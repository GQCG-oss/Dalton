!> @file 
!> MPI routines for main DEC driver.
!> \author Kasper Kristensen

#ifdef VAR_MPI
module dec_driver_slave_module


  use precision
  use lstiming
  use infpar_module
  use lsmpi_type
  use typedeftype
  use dec_typedef_module
  use DALTONINFO, only: ls_free
  use BUILDAOBATCH


  ! DEC DEPENDENCIES (within deccc directory) 
  ! *****************************************
  use dec_fragment_utils
  use dec_settings_mod!, only: dec_set_default_config
  use full_molecule!, only: molecule_finalize
  use decmpi_module
  use atomic_fragment_operations
  use mp2_module
  use ccsd_module
  use mp2_gradient_module
  use fragment_energy_module

public :: main_fragment_driver_slave
private

contains

!> \brief Handling of local master ranks in the main fragment driver for DEC.
!> Local masters (infpar%lg_comm=0): Receive fragment jobs from main master and carry these out.
!> Local slaves (infpar%lg_com/=0) : Help local master using MPI parallelization of
!> integral/amplitude routine (MP2) or coupled-cluster residual (CC2 or CCSD).
!> \author Kasper Kristensen
!> \date March 2012
  subroutine main_fragment_driver_slave()


    implicit none
    integer :: job,i
    INTEGER(kind=ls_mpik) :: MPISTATUS(MPI_STATUS_SIZE), DUMMYSTAT(MPI_STATUS_SIZE)
    integer :: nunocc,nocc,natoms,siz
    type(lsitem) :: MyLsitem
    type(decfrag),pointer :: AtomicFragments(:),EstAtomicFragments(:)
    type(fullmolecule) :: MyMolecule
    type(decorbital),pointer :: OccOrbitals(:), UnoccOrbitals(:)
    logical,pointer :: dofrag(:)
    type(joblist) :: jobs
    integer :: groups,signal,step
    integer(kind=ls_mpik) :: groupsize,ierr
    logical :: dens_save,FO_save,grad_save, localslave, in_master_group,esti
    integer(kind=ls_mpik) :: master
    master=0


    ! DEC settings
    ! ************
    ! Just in case, initialize using default settings
    call dec_set_default_config(0)
    ! Get DEC settings for current calculation from master
    call mpibcast_dec_settings(DECinfo,MPI_COMM_LSDALTON)

    ! Set output unit number to 0 for slaves
    DECinfo%output=0


    ! GET FULL MOLECULAR INFO FROM MASTER
    ! ***********************************

    ! Get very basic dimension information from main master
    call mpi_dec_fullinfo_master_to_slaves_precursor(esti,nocc,nunocc,master)

    ! Allocate arrays needed for fragment calculations
    call mem_alloc(OccOrbitals,nocc)
    call mem_alloc(UnoccOrbitals,nunocc)

    ! Get remaining full molecular information needed for all fragment calculations
    call mpi_dec_fullinfo_master_to_slaves(nocc,nunocc,&
         & OccOrbitals, UnoccOrbitals, MyMolecule, MyLsitem)
    natoms=MyMolecule%natoms

    ! Get list of which atoms have orbitals assigned
    call mem_alloc(dofrag,natoms)
    call which_atoms_have_orbitals_assigned(MyMolecule%ncore,nocc,nunocc,&
         & natoms,OccOrbitals,UnoccOrbitals,dofrag,MyMolecule%PhantomAtom)

    ! Internal control of first order property keywords
    ! (Necessary because these must be false during fragment optimization.)
    ! Better solution should be implemented at some point...
    dens_save = DECinfo%MP2density
    FO_save = DECinfo%first_order
    grad_save = DECinfo%gradient
    DECinfo%MP2density=.false.
    DECinfo%first_order=.false.
    DECinfo%gradient=.false.

    ! Atomic fragments
    call mem_alloc(AtomicFragments,natoms)
    do i=1,natoms
       call atomic_fragment_nullify(AtomicFragments(i))
    end do


    ! *************************************************************
    ! *                   FRAGMENT CALCULATIONS                   *
    ! *************************************************************

    ! Two steps
    ! ---------
    ! step=1: Atomic fragment optimizations (and possibly estimated pair fragment calculations)
    ! step=2: Pair fragment calculations for given FOT (and possibly repeating atomic fragments)

    StepLoop: do step=1,2


       ! Special for step 1: Get estimated fragments
       ! *******************************************
       if(step==1 .and. esti) then
          ! Get optimized atomic fragments from master node
          call mem_alloc(EstAtomicFragments,natoms)
          do i=1,natoms
             call atomic_fragment_nullify(EstAtomicFragments(i))
          end do
          call mpi_bcast_many_fragments(natoms,dofrag,EstAtomicFragments,MPI_COMM_LSDALTON)          

          ! Receive CC models to use for each pair based on estimates
          call ls_mpibcast(MyMolecule%ccmodel,natoms,natoms,master,MPI_COMM_LSDALTON)
       end if
       

       ! Special for step 2: Reset first-order keywords and get optimized fragments
       ! **************************************************************************
       Step2: if(step==2) then

          ! Restore first order keywords for second step
          DECinfo%MP2density=dens_save
          DECinfo%first_order = FO_save
          DECinfo%gradient = grad_save

          ! Get FOT-optimized atomic fragments from master node
          call mpi_bcast_many_fragments(natoms,dofrag,AtomicFragments,MPI_COMM_LSDALTON)

       end if Step2


       ! Initialize local MPI groups
       ! ***************************

       ! Receive signal from master to init group (sanity precaution)
       call ls_mpibcast(signal,infpar%master,MPI_COMM_LSDALTON)
       if(signal /= GROUPINIT) then
          call lsquit('main_fragment_driver_slave: Expected signal to init group',-1)
       end if

       ! Get groupsize from master
       call ls_mpibcast(groupsize,master,MPI_COMM_LSDALTON)
       if(DECinfo%PL>0) print *, 'node/groupsize', infpar%mynum, groupsize

       ! Initialize new MPI groups
       call init_mpi_groups(groupsize,DECinfo%output)

       ! Local master if rank=0 WITHIN the local group (might be different in steps 1 and 2)
       if(infpar%lg_mynum/=master) then ! local slave
          localslave=.true.
       else ! local master
          localslave = .false.
       end if

       !  Get fragment job list (includes initialization of pointers in job list)
       call bcast_dec_fragment_joblist(jobs,MPI_COMM_LSDALTON)


       ! REDEFINING LOCAL GROUPS
       ! =======================
       infpar%lg_morejobs = .true.   ! more jobs for local group
       do while(localslave)

          ! Send local slave to local waiting position until local master contacts it
          call DEC_lsmpi_slave(infpar%lg_comm)

          ! Slaves were kicked out of waiting position - two possibilities:
          !
          ! 1. lg_morejobs=.true.
          ! There are more jobs to be done - divide local group, and then go back to local slave position
          ! unless this rank has become a new local master (infpar%lg_mynum=master)
          !
          ! 2. lg_morejobs=.false.
          ! There are no more jobs to be done - quit do-loop regardless of current rank number

          if(.not. infpar%lg_morejobs) exit

          ! Get new local groups
          ! (redefine globals infpar%lg_mynum, infpar%lg_nodtot, and infpar%lg_comm)
          call dec_half_local_group

          ! Check if the current rank has become a local master (rank=master within local group)
          if(infpar%lg_mynum==master) localslave=.false.

       end do


       ! Receive  fragment jobs from master and carry those out
       ! ======================================================
       if(step==1 .and. esti) then
          ! Calculate estimated pair fragment energies
          call fragments_slave(natoms,nocc,nunocc,OccOrbitals,&
               & UnoccOrbitals,MyMolecule,MyLsitem,AtomicFragments,jobs,&
               & EstAtomicFragments=EstAtomicFragments)
       else
          call fragments_slave(natoms,nocc,nunocc,OccOrbitals,&
               & UnoccOrbitals,MyMolecule,MyLsitem,AtomicFragments,jobs)
       end if

       ! Remaining local slaves should exit local slave routine for good (infpar%lg_morejobs=.false.)
       job=QUITNOMOREJOBS
       if(infpar%lg_mynum==master .and. infpar%lg_nodtot>1) then
          call ls_mpibcast(job,master,infpar%lg_comm)
       end if

       ! Done with existing job list
       call free_joblist(jobs)

       ! Free existing group communicators
       call MPI_COMM_FREE(infpar%lg_comm, IERR)

       ! Clean up estimated fragments used in step 1 and receive CC models to use for all pairs
       CleanupAndUpdateCCmodel: if(step==1 .and. esti) then
          do i=1,natoms
             if(dofrag(i)) then
                call atomic_fragment_free_simple(EstAtomicFragments(i))
             end if
          end do
          call mem_dealloc(EstAtomicFragments)

          ! Receive CC models to use for each pair based on estimates
          call ls_mpibcast(MyMolecule%ccmodel,natoms,natoms,master,MPI_COMM_LSDALTON)

       end if CleanupAndUpdateCCmodel

    end do StepLoop


    ! Free stuff
    ! **********
    do i=1,natoms
       if(dofrag(i)) then
          call atomic_fragment_free_simple(AtomicFragments(i))
       end if
    end do
    call mem_dealloc(AtomicFragments)
    do i=1,nOcc
       call orbital_free(OccOrbitals(i))
    end do
    do i=1,nUnocc
       call orbital_free(UnoccOrbitals(i))
    end do
    call mem_dealloc(OccOrbitals)
    call mem_dealloc(UnoccOrbitals)
    call mem_dealloc(dofrag)
    call ls_free(MyLsitem)
    call molecule_finalize(MyMolecule)

  end subroutine main_fragment_driver_slave



!> \brief For each local master: Carry out  fragment calculations (both single and pair)
!> for those  fragments requested by main master.
!> \author Kasper Kristensen
!> \date May 2012
  subroutine fragments_slave(natoms,nocc,nunocc,OccOrbitals,&
       & UnoccOrbitals,MyMolecule,MyLsitem,AtomicFragments,jobs,EstAtomicFragments)

    implicit none

    !> Number of atoms in the molecule
    integer,intent(in) :: natoms
    !> Number of occupied orbitals in the molecule
    integer,intent(in) :: nocc
    !> Number of unoccupied orbitals in the molecule
    integer,intent(in) :: nunocc
    !> Occupied orbitals, DEC format
    type(decorbital),intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals, DEC format
    type(decorbital),intent(in) :: UnoccOrbitals(nunocc)
    !> Full molecular info (not changed at output)
    type(fullmolecule),intent(inout) :: MyMolecule
    !> LS item structure
    type(lsitem),intent(inout) :: MyLsitem
    !> Atomic fragments
    type(decfrag),dimension(natoms),intent(inout) :: AtomicFragments
    !> Job list for  fragments
    type(joblist),intent(in) :: jobs
    !> Estimated atomic fragments
    type(decfrag),dimension(natoms),intent(inout),optional :: EstAtomicFragments
    type(joblist) :: singlejob
    type(decfrag) :: PairFragment
    type(decfrag) :: MyFragment
    integer :: job,atomA,atomB,i,slavejob,ntasks
    real(realk) :: flops_slaves
    type(mp2grad) :: grad
    logical :: morejobs, divide, split,only_update
    real(realk) :: fragenergy(ndecenergies),tottime
    real(realk) :: t1cpu, t2cpu, t1wall, t2wall
    real(realk) :: t1cpuacc, t2cpuacc, t1wallacc, t2wallacc
    real(realk) :: flops
    integer(kind=ls_mpik) :: master
    master = 0
    fragenergy=0.0_realk
    only_update=.true.

    call LSTIMER('START',t1cpuacc,t1wallacc,DECinfo%output)

    if(any(jobs%esti)) then
       if(.not. present(EstAtomicFragments)) then
          call lsquit('fragment_slaves: Estimated pair fragments requested, but estimated & 
               & atomic fragments are not present!',-1)
       end if
    end if



    ! ********************************************************************
    ! *                        FRAGMENT CALCULATIONS                     *
    ! ********************************************************************

    ! Local master asks main master for fragment jobs
    slavejob=QUITMOREJOBS
    divide=.true.

    ! Send empty job in first round
    job=0
    morejobs=.true.
    ! Job list containing only a single job which will be overwritten in each new fragment calculation
    call init_joblist(1,singlejob)


    AskForJob: do while(morejobs)
      


       ! Send finished job to master
       ! ***************************
       ! (If job=0 it is an empty job not containing any real information)
       call ls_mpisendrecv(job,MPI_COMM_LSDALTON,infpar%mynum,master)

       ! If real job: Send info to master
       SendToMaster: if(job>0) then

          ! Update stuff (same for single and pair)
          ! ***************************************

          FragoptCheck1: if(jobs%dofragopt(job)) then  
             ! Job is fragment optimization --> send atomic fragment info to master
             call mpi_send_recv_single_fragment(MyFragment,MPI_COMM_LSDALTON,&
                  & infpar%mynum,master,singlejob)
             call atomic_fragment_free_simple(MyFragment)

          else

             ! No fragment optimization done --> send fragment energy to master

             if(DECinfo%first_order) then ! also send density,gradient
                call mpi_send_recv_mp2grad(MPI_COMM_LSDALTON,infpar%mynum,master,&
                     & grad,fragenergy,singlejob,only_update)
                call free_mp2grad(grad)
             else ! just energy
                call mpi_send_recv_fragmentenergy(MPI_COMM_LSDALTON,infpar%mynum,master,&
                     & fragenergy,singlejob)
             end if

          end if FragoptCheck1


       end if SendToMaster


       ! Receive new job task
       ! ********************
       call ls_mpisendrecv(job,MPI_COMM_LSDALTON,master,infpar%mynum)

       ! Carry out fragment optimization (job>0), or finish if all jobs are done (job=-1)
       ! ********************************************************************************
       DoJob: IF(job==-1) then
          call LSTIMER('START',t2cpuacc,t2wallacc,DECinfo%output)
          print '(a,i6,a,g14.6)', 'Slave ', infpar%mynum, ' exits  fragments. Time: ',&
               & t2wallacc-t1wallacc
          morejobs=.false.
       else

          ! Timing and flops
          call LSTIMER('START',t1cpu,t1wall,DECinfo%output)
          call start_flop_counter()

          atomA=jobs%atom1(job)
          atomB=jobs%atom2(job)

          ! Find number of tasks required by integral program
          ! *************************************************
          if(atomA==atomB) then ! single fragment
             FragoptCheck2: if(.not. jobs%dofragopt(job)) then
                ! Init basis for fragment
                call atomic_fragment_init_basis_part(nunocc, nocc, OccOrbitals,&
                     & UnoccOrbitals,MyMolecule,mylsitem,AtomicFragments(atomA))

                call get_number_of_integral_tasks_for_mpi(AtomicFragments(atomA),ntasks)
  
             else
                ! Set ntasks to be zero to initialize it to something, although it is not used for
                ! fragment optimizations.
                ! Note that the groups never divide for fragment optimizations below.
                ntasks=0
             end if FragoptCheck2
          else ! pair fragment

             ! init pair
             if(jobs%esti(job)) then
                ! Estimated pair fragment
                call merged_fragment_init(EstAtomicFragments(atomA), EstAtomicFragments(atomB),&
                     & nunocc, nocc, natoms,OccOrbitals,UnoccOrbitals, &
                     & MyMolecule,mylsitem,.true.,PairFragment,esti=.true.)
             else
                call merged_fragment_init(AtomicFragments(atomA), AtomicFragments(atomB),&
                     & nunocc, nocc, natoms,OccOrbitals,UnoccOrbitals, &
                     & MyMolecule,mylsitem,.true.,PairFragment)
             end if
             call get_number_of_integral_tasks_for_mpi(PairFragment,ntasks)

          end if

          ! Divide local group into two smaller groups if:
          ! number of tasks < (#nodes in local group)*(magic factor)
          !
          ! The magic factor (DECinfo%MPIsplit) should in general be larger than 1, 
          ! and it can be changed by .MPIsplit keyword.
          !
          ! Thus, we obtain optimal MPI performance for the individual fragment
          ! calculations if the number of tasks is significantly larger
          ! than the number of nodes in the local group.
          ! Also check that the local group is larger than just 1 node.
          !

          split=.false.
          divide=.true.
          ! Special case: Never divide for fragment optimization
          if(jobs%dofragopt(job)) divide=.false.

          DoDivide: do while(divide)
             if( (ntasks < infpar%lg_nodtot*DECinfo%MPIsplit) .and. (infpar%lg_nodtot>1)  ) then

                print '(a,4i8)', 'DIVIDE! Job, tasks, mynode, #nodes: ', &
                     & job,ntasks,infpar%mynum,infpar%lg_nodtot

                ! Kick slaves out of local group
                call ls_mpibcast(slavejob,master,infpar%lg_comm)

                ! Divide local group into two smaller groups
                call dec_half_local_group
                split=.true.
             else
                divide=.false.
             end if
          end do DoDivide

          print '(a,i8,a,i8)', 'Slave ', infpar%mynum, ' will do  job ', job

          ! Communicator in setting may have changed due to division of local group
          ! Ensure that the correct local communicator is used.
          if(split) then
             if(atomA==atomB) then
                AtomicFragments(atomA)%mylsitem%setting%comm = infpar%lg_comm
             else
                PairFragment%mylsitem%setting%comm = infpar%lg_comm
             end if
          end if


          ! RUN FRAGMENT CALCULATION
          ! ************************

          FragoptCheck3: if(jobs%dofragopt(job)) then
             ! Sanity check
             if(atomA/=atomB) then
                call lsquit('fragments_slave: Fragment optimization requested for pair fragment!',-1)
             end if

             ! Fragment optimization
             call optimize_atomic_fragment(atomA,MyFragment,MyMolecule%nAtoms, &
                  & OccOrbitals,nOcc,UnoccOrbitals,nUnocc, &
                  & MyMolecule,mylsitem,.true.)
             flops_slaves = MyFragment%flops_slaves
             tottime = MyFragment%slavetime ! time used by all local slaves
             fragenergy = MyFragment%energies
             call copy_fragment_info_job(MyFragment,singlejob)

          else


             SingleOrPair: if( atomA == atomB ) then ! single  fragment

                call atomic_driver(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals,&
                     & AtomicFragments(atomA),grad=grad)
                flops_slaves = AtomicFragments(atomA)%flops_slaves
                tottime = AtomicFragments(atomA)%slavetime ! time used by all local slaves
                fragenergy = AtomicFragments(atomA)%energies
                call copy_fragment_info_job(AtomicFragments(atomA),singlejob)
                call atomic_fragment_free_basis_info(AtomicFragments(atomA))

             else ! pair  fragment

                if(jobs%esti(job)) then
                   ! Estimated pair fragment
                   call pair_driver(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals,&
                        & EstAtomicFragments(atomA), EstAtomicFragments(atomB), &
                        & natoms,PairFragment,grad)
                else
                   ! Pair fragment according to FOT precision
                   call pair_driver(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals,&
                        & AtomicFragments(atomA), AtomicFragments(atomB), &
                        & natoms,PairFragment,grad)
                end if
                flops_slaves = PairFragment%flops_slaves
                tottime = PairFragment%slavetime ! time used by all local slaves
                fragenergy=PairFragment%energies

                call copy_fragment_info_job(PairFragment,singlejob)
                ! Free pair
                call atomic_fragment_free(PairFragment)

             end if SingleOrPair

          end if FragoptCheck3


          ! Set fragment job info
          ! *********************
          call LSTIMER('START',t2cpu,t2wall,DECinfo%output)
          call end_flop_counter(flops) ! flops for local master
          singlejob%LMtime(1) = t2wall - t1wall  ! wall time used by local master
          tottime = tottime + singlejob%LMtime(1) ! total time for all local slaves and local master
          singlejob%flops(1) = flops + flops_slaves  ! FLOPS for local master + local slaves
          singlejob%nslaves(1) = infpar%lg_nodtot ! Sizes of local slot (local master + local slaves)
          ! load distribution: { tottime / time(local master) } / number of nodes (ideally 1.0)
          singlejob%load(1) = (tottime/singlejob%LMtime(1))/real(singlejob%nslaves(1))
          singlejob%jobsdone(1) = .true.
          singlejob%esti(1) = jobs%esti(job)
          singlejob%dofragopt(1) = jobs%dofragopt(job)

          print '(a,i8,a,i8,g14.6)', 'Slave ', infpar%mynum, ' is done with  job/time ', &
               & job, singlejob%LMtime(1)
       end if DoJob
   


    end do AskForJob


    call free_joblist(singlejob)

  end subroutine fragments_slave


!> \brief Get number of tasks in integral loop (nalpha*ngamma)
!> \author Kasper Kristensen
!> \date May 2012
subroutine get_number_of_integral_tasks_for_mpi(MyFragment,ntasks)

  implicit none
  !> Atomic fragment
  type(decfrag),intent(inout) :: MyFragment
  !> Number of tasks
  integer,intent(inout) :: ntasks
  type(mp2_batch_construction) :: bat
  integer :: MaxActualDimAlpha,nbatchesAlpha
  integer :: MaxActualDimGamma,nbatchesGamma
  integer, pointer :: orb2batchAlpha(:), batchdimAlpha(:), batchsizeAlpha(:), batchindexAlpha(:)
  integer, pointer :: orb2batchGamma(:), batchdimGamma(:), batchsizeGamma(:), batchindexGamma(:)
  logical :: mpi_split

  ! Initialize stuff (just dummy arguments here)
  nullify(orb2batchAlpha)
  nullify(batchdimAlpha)
  nullify(batchsizeAlpha)
  nullify(batchindexAlpha)
  nullify(orb2batchGamma)
  nullify(batchdimGamma)
  nullify(batchsizeGamma)
  nullify(batchindexGamma)

  ntasks = 0

  ! Determine optimal batchsizes with available memory
  if(MyFragment%ccmodel==MODEL_MP2) then ! MP2
     call get_optimal_batch_sizes_for_mp2_integrals(MyFragment,DECinfo%first_order,bat,.false.)
  else  ! CC2 or CCSD
     mpi_split = .true.
     call wrapper_get_ccsd_batch_sizes(MyFragment,bat,mpi_split,ntasks)
     ! In case of MO-based algorithm the number of tasks depends on the
     ! number of MO batches and is directly returned by the routine.
     if (ntasks>0) return
  end if


  ! Get number of gamma batches
  call mem_alloc(orb2batchGamma,MyFragment%nbasis)
  call build_batchesofAOS(DECinfo%output,MyFragment%mylsitem%setting,bat%MaxAllowedDimGamma,&
       & MyFragment%nbasis,MaxActualDimGamma,batchsizeGamma,batchdimGamma,&
       & batchindexGamma,nbatchesGamma,orb2BatchGamma,'R')
  call mem_dealloc(orb2batchGamma)
  call mem_dealloc(batchdimGamma)
  call mem_dealloc(batchsizeGamma)
  call mem_dealloc(batchindexGamma)

  ! Get number of alpha batches
  call mem_alloc(orb2batchAlpha,MyFragment%nbasis)
  call build_batchesofAOS(DECinfo%output,MyFragment%mylsitem%setting,bat%MaxAllowedDimAlpha,&
       & MyFragment%nbasis,MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,&
       & batchindexAlpha,nbatchesAlpha,orb2BatchAlpha,'R')
  call mem_dealloc(orb2batchAlpha)
  call mem_dealloc(batchdimAlpha)
  call mem_dealloc(batchsizeAlpha)
  call mem_dealloc(batchindexAlpha)


  ! Number of tasks = nalpha*ngamma
  ntasks = nbatchesGamma*nbatchesAlpha

end subroutine get_number_of_integral_tasks_for_mpi


end module dec_driver_slave_module


#else
module dec_driver_slave_module

contains

!Added to avoid "has no symbols" linking warning
subroutine dec_driver_slave_void()
end subroutine dec_driver_slave_void

end module dec_driver_slave_module

#endif

