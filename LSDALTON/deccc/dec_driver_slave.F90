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
  type(ccatom),pointer :: AtomicFragments(:)
  type(fullmolecule) :: MyMolecule
  real(realk),pointer :: DistanceTable(:,:)
  type(ccorbital),pointer :: OccOrbitals(:), UnoccOrbitals(:)
  logical,pointer :: dofrag(:)
  type(joblist) :: jobs
  integer :: groups,signal
  integer(kind=ls_mpik) :: groupsize,ierr
  logical :: dens_save,FO_save,grad_save, localslave, in_master_group, morejobs
  integer(kind=ls_mpik) :: master
  master=0


  ! DEC settings
  ! ************
  ! Just in case, initialize using default settings
  call dec_set_default_config(0)
  siz=sizeof(DECinfo)
  ! Get DEC settings for current calculation from master (quick and dirty bcast...)
  ierr=0

    call MPI_BCAST(DECinfo,siz,MPI_CHARACTER,master,MPI_COMM_LSDALTON,ierr)

  if(ierr/=0) then
     call  lsquit('main_fragment_driver_slave: Something went wrong when &
          & bcasting DEC structure!',-1)
  end if
  ! Set output unit number to 0 for slaves
  DECinfo%output=0

  ! Internal control of first order property keywords
  ! (Necessary because these must be false during fragment optimization.)
  ! Better solution should be implemented at some point...
  dens_save = DECinfo%MP2density
  FO_save = DECinfo%first_order
  grad_save = DECinfo%gradient
  DECinfo%MP2density=.false.
  DECinfo%first_order=.false.
  DECinfo%gradient=.false.


  ! GET FULL MOLECULAR INFO FROM MASTER
  ! ***********************************

  ! Get very basic dimension information from main master
  call ls_mpibcast(natoms,master,MPI_COMM_LSDALTON)
  call ls_mpibcast(nocc,master,MPI_COMM_LSDALTON)
  call ls_mpibcast(nunocc,master,MPI_COMM_LSDALTON)

  ! Allocate arrays needed for fragment calculations
  call mem_alloc(OccOrbitals,nocc)
  call mem_alloc(UnoccOrbitals,nunocc)
  call mem_alloc(DistanceTable,natoms,natoms)

  ! Get remaining full molecular information needed for all fragment calculations
  call mpi_dec_fullinfo_master_to_slaves(natoms,nocc,nunocc,DistanceTable,&
       & OccOrbitals, UnoccOrbitals, MyMolecule, MyLsitem)


  ! Local master if rank=0 WITHIN the local group
  if(infpar%lg_mynum/=master) then ! local slave
     localslave=.true.
  else ! local master
     localslave = .false.
  end if

  if(localslave) then
     ! Send local slave to local waiting position until local master contacts it
     call DEC_lsmpi_slave(infpar%lg_comm)
  end if



  ! **********************************************
  ! *        ATOMIC FRAGMENT OPTIMIZATION        *
  ! **********************************************

  call atomic_fragments_slave(natoms,nocc,nunocc,DistanceTable,OccOrbitals,&
       & UnoccOrbitals,MyMolecule,MyLsitem)

  ! Remaining local slaves should exit local slave routine (but there will be more jobs)
  job=QUITMOREJOBS
  if(infpar%lg_mynum==0 .and. infpar%lg_nodtot>1) then
     call ls_mpibcast(job,master,infpar%lg_comm)
  end if



  ! ********************************************************************
  ! *       Get new optimized atomic fragments from master node        *
  ! ********************************************************************
  ! 1. Get list of which atoms are EOS atoms
  call mem_alloc(dofrag,natoms)
  dofrag=.false.
  call ls_mpibcast(dofrag,natoms,master,MPI_COMM_LSDALTON)

  ! 2. Get  fragments
  call mem_alloc(AtomicFragments,natoms)
  call mpi_bcast_many_fragments(natoms,dofrag,AtomicFragments,MPI_COMM_LSDALTON)




     ! *************************************************************
     ! *         POST-FRAGOPT FRAGMENT CALCULATIONS                *
     ! *************************************************************

     ! Init stuff
     ! **********

     ! Restore first order keywords
     DECinfo%MP2density=dens_save
     DECinfo%first_order = FO_save
     DECinfo%gradient = grad_save


     ! Continue as long as there are more jobs to be done
     ! (we might need to add more pairs than in original  job list to adapt to precision)
     morejobs=.true.
     MorePairs: do while(morejobs)

        !  Get  fragment job list (includes initalization of pointers in job list)
        call bcast_post_fragopt_joblist(jobs,MPI_COMM_LSDALTON)
        ! Receive pair distance threshold to make sure we are consistent
        call ls_mpibcast(DECinfo%pair_distance_threshold,master,MPI_COMM_LSDALTON)


        ! Redefine MPI groups for  fragments
        ! ***************************************
        call lsmpi_barrier(MPI_COMM_LSDALTON)

        ! Free existing group communicators
        call MPI_COMM_FREE(infpar%lg_comm, IERR)

        ! Receive signal from master to init group/sanity check
        call ls_mpibcast(signal,infpar%master,MPI_COMM_LSDALTON)
        if(signal /= GROUPINIT) then
           call lsquit('main_fragment_driver_slave: Expected signal to init group',-1)
        end if
     
        ! Get groupsize from master
        call ls_mpibcast(groupsize,master,MPI_COMM_LSDALTON)
        if(DECinfo%PL>0) print *, 'node/groupsize', infpar%mynum, groupsize

        ! Initialize new group
        call init_mpi_groups(groupsize,DECinfo%output)

        ! Local master/slave assigment might have changed - reset
        if(infpar%lg_mynum/=0) then ! local slave
           localslave=.true.
        else ! local master
           localslave = .false.
        end if

        
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
           ! unless this rank has become a new local master (infpar%lg_mynum=0)
           !
           ! 2. lg_morejobs=.false.
           ! There are no more jobs to be done - quit do-loop regardless of current rank number
        
           if(.not. infpar%lg_morejobs) exit

           ! Get new local groups
           ! (redefine globals infpar%lg_mynum, infpar%lg_nodtot, and infpar%lg_comm)
           call dec_half_local_group
           
           ! Check if the current rank has become a local master (rank=0 within local group)
           if(infpar%lg_mynum==0) localslave=.false.
        
        end do

  
        ! Receive  fragment jobs from master and carry those out
        ! ======================================================
        call fragments_slave(natoms,nocc,nunocc,DistanceTable,OccOrbitals,&
             & UnoccOrbitals,MyMolecule,MyLsitem,AtomicFragments,jobs)

        ! Remaining local slaves should exit local slave routine for good (infpar%lg_morejobs=.false.)
        job=QUITNOMOREJOBS
        if(infpar%lg_mynum==0 .and. infpar%lg_nodtot>1) then
           call ls_mpibcast(job,master,infpar%lg_comm)
        end if
     
        ! Receive information telling whether there are more jobs or not
        call ls_mpibcast(morejobs,master,MPI_COMM_LSDALTON)
        
        ! Done with existing job list
        call free_joblist(jobs)
        
     end do MorePairs


  ! Free stuff
  ! **********
  do i=1,natoms
     if(dofrag(i)) then
        call atomic_fragment_free_simple(AtomicFragments(i))
     end if
  end do
  do i=1,nOcc
     call orbital_free(OccOrbitals(i))
  end do
  do i=1,nUnocc
     call orbital_free(UnoccOrbitals(i))
  end do
  call mem_dealloc(dofrag)
  call mem_dealloc(AtomicFragments)
  call mem_dealloc(OccOrbitals)
  call mem_dealloc(UnoccOrbitals)
  call mem_dealloc(DistanceTable)
  call ls_free(MyLsitem)
  call molecule_finalize(MyMolecule)

  ! Free existing group communicators
  call MPI_COMM_FREE(infpar%lg_comm, IERR)


end subroutine main_fragment_driver_slave


!> \brief For each local master: Carry out fragment optimization for
!> those atomic fragments requested by main master.
!> \author Kasper Kristensen
!> \date May 2012
subroutine atomic_fragments_slave(natoms,nocc,nunocc,DistanceTable,OccOrbitals,&
     & UnoccOrbitals,MyMolecule,MyLsitem)


  implicit none

  !> Number of atoms in the molecule
  integer,intent(in) :: natoms
  !> Number of occupied orbitals in the molecule
  integer,intent(in) :: nocc
  !> Number of unoccupied orbitals in the molecule
  integer,intent(in) :: nunocc
  !> Distance table with atomic distances
  real(realk) :: DistanceTable(natoms,natoms)
  !> Occupied orbitals, DEC format
  type(ccorbital),intent(in) :: OccOrbitals(nocc)
  !> Unoccupied orbitals, DEC format
  type(ccorbital),intent(in) :: UnoccOrbitals(nunocc)
  !> Full molecular info
  type(fullmolecule),intent(in) :: MyMolecule
  !> LS item structure
  type(lsitem),intent(inout) :: MyLsitem
  !> Atomic fragments to be determined
  type(ccatom) :: AtomicFragment
  type(joblist) :: job  
  real(realk) :: flops
  integer :: jobidx
  real(realk) :: t1cpu, t2cpu, t1wall, t2wall, tottime
  real(realk) :: t1cpuacc, t2cpuacc, t1wallacc, t2wallacc
  integer(kind=ls_mpik) :: master
  master = 0

  call LSTIMER('START',t1cpuacc,t1wallacc,DECinfo%output)

  ! Job list containing only a single job which will be overwritten in each new fragment calculation
  call init_joblist(1,job)
  

  ! ********************************************************************
  ! *                    ATOMIC FRAGMENT OPTIMIZATION                  *
  ! ********************************************************************
  ! Local master asks main master for fragment jobs


  ! Send empty job in first round
  jobidx=0

  FragmentOptimization: do while(.true.)


     ! Send finished job to master
     ! ***************************
     ! (If jobidx=0 it is an empty job not containing any real information)
     call ls_mpisendrecv(jobidx,MPI_COMM_LSDALTON,infpar%mynum,master)

     ! Real job: Send fragment info to master
     if(jobidx/=0) then
        call mpi_send_recv_single_fragment(AtomicFragment,MPI_COMM_LSDALTON,&
             & infpar%mynum,master,job)
        call atomic_fragment_free_simple(AtomicFragment)
     end if

     ! Receive new job task
     ! ********************
     call ls_mpisendrecv(jobidx,MPI_COMM_LSDALTON,master,infpar%mynum)
     print '(a,i6,a,i6)', 'Slave ', infpar%mynum, ' will do atomic job ', jobidx

     ! Carry out fragment optimization (jobidx>0), or finish if all jobs are done (jobidx=-1)
     ! ********************************************************************************
     IF(jobidx==-1) then
        call LSTIMER('START',t2cpuacc,t2wallacc,DECinfo%output)
        print '(a,i6,a,g14.6)', 'Slave ', infpar%mynum, ' exits fragment optimization. Time: ',&
             & t2wallacc-t1wallacc
        exit
     else
        ! Timing and flops
        call LSTIMER('START',t1cpu,t1wall,DECinfo%output)
        call start_flop_counter()

        ! Carry out fragment optimization job task
        call optimize_atomic_fragment(jobidx,AtomicFragment,nAtoms, &
             & OccOrbitals,nOcc,UnoccOrbitals,nUnocc,DistanceTable, &
             & MyMolecule,mylsitem,.true.)

        ! Set job info (statistics)
        ! *************************
        call LSTIMER('START',t2cpu,t2wall,DECinfo%output)
        call end_flop_counter(flops) ! flops for local master

        call copy_fragment_info_job(AtomicFragment,job)

        ! Job is done...
        job%jobsdone(1) = .true.
        ! wall time used by local master
        job%LMtime(1) = t2wall - t1wall
        ! total time for all local slaves and local master
        tottime = AtomicFragment%slavetime + job%LMtime(1)
        ! FLOPS for local master + local slaves
        job%flops(1) = flops + AtomicFragment%flops_slaves
        job%nslaves(1) = infpar%lg_nodtot ! Sizes of local slot (local master + local slaves)
        ! load distribution: { tottime / time(local master) } / number of nodes (ideally 1.0)
        job%load(1)= (tottime/job%LMtime(1))/job%nslaves(1)

     end if
     print '(a,i6,a,i6,g14.6)', 'Local master ', infpar%mynum, ' is done with atomic job/time ', &
          & jobidx, job%LMtime(1)

  end do FragmentOptimization


  call free_joblist(job)

end subroutine atomic_fragments_slave




!> \brief For each local master: Carry out  fragment calculations (both single and pair)
!> for those  fragments requested by main master.
!> \author Kasper Kristensen
!> \date May 2012
subroutine fragments_slave(natoms,nocc,nunocc,DistanceTable,OccOrbitals,&
     & UnoccOrbitals,MyMolecule,MyLsitem,AtomicFragments,jobs)

  implicit none

  !> Number of atoms in the molecule
  integer,intent(in) :: natoms
  !> Number of occupied orbitals in the molecule
  integer,intent(in) :: nocc
  !> Number of unoccupied orbitals in the molecule
  integer,intent(in) :: nunocc
  !> Distances between atoms
  real(realk) :: DistanceTable(natoms,natoms)
  !> Occupied orbitals, DEC format
  type(ccorbital),intent(in) :: OccOrbitals(nocc)
  !> Unoccupied orbitals, DEC format
  type(ccorbital),intent(in) :: UnoccOrbitals(nunocc)
  !> Full molecular info (not changed at output)
  type(fullmolecule),intent(inout) :: MyMolecule
  !> LS item structure
  type(lsitem),intent(inout) :: MyLsitem
  !> Atomic fragments
  type(ccatom),dimension(natoms),intent(inout) :: AtomicFragments
  !> Job list for  fragments
  type(joblist),intent(inout) :: jobs
  type(joblist) :: singlejob
  type(ccatom) :: PairFragment
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



  ! ********************************************************************
  ! *               POST-FRAGOPT FRAGMENT CALCULATIONS                 *
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
        if(DECinfo%first_order) then ! density,gradient
           call mpi_send_recv_mp2grad(MPI_COMM_LSDALTON,infpar%mynum,master,&
                & grad,fragenergy,singlejob,only_update)
           call free_mp2grad(grad)
        else ! just energy
           call mpi_send_recv_fragmentenergy(MPI_COMM_LSDALTON,infpar%mynum,master,&
                & fragenergy,singlejob)
        end if

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

           ! Init basis for fragment
           call atomic_fragment_init_basis_part(nunocc, nocc, OccOrbitals,&
                & UnoccOrbitals,MyMolecule,mylsitem,AtomicFragments(atomA))

           call get_number_of_integral_tasks_for_mpi(AtomicFragments(atomA),ntasks)
        else ! pair fragment

           ! init pair
           call merged_fragment_init(AtomicFragments(atomA), AtomicFragments(atomB),&
                & nunocc, nocc, natoms,OccOrbitals,UnoccOrbitals, DistanceTable, &
                & MyMolecule,mylsitem,.true.,PairFragment)
           call get_number_of_integral_tasks_for_mpi(PairFragment,ntasks)

        end if

        ! Divide local group into two smaller groups if:
        ! number of tasks < (#nodes in local group)*(magic factor)
        !
        ! The magic factor is 10 by default and can be changed by .MPIsplit keyword.
        !
        ! Thus, we obtain optimal MPI performance for the individual fragment
        ! calculations if the number of tasks is significantly larger
        ! than the number of nodes in the local group.
        ! (The factor 10 is a "rule-of-thumb").
        ! Also check that the local group is larger than just 1 node.

        split=.false.
        divide=.true.
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
        SingleOrPair: if( atomA == atomB ) then ! single  fragment

           call atomic_driver(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals,&
                & AtomicFragments(atomA),grad=grad)
           flops_slaves = AtomicFragments(atomA)%flops_slaves
           tottime = AtomicFragments(atomA)%slavetime ! time used by all local slaves
           fragenergy = AtomicFragments(atomA)%energies
           call copy_fragment_info_job(AtomicFragments(atomA),singlejob)
           call atomic_fragment_free_basis_info(AtomicFragments(atomA))

        else ! pair  fragment

           call pair_driver(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals,&
                & AtomicFragments(atomA), AtomicFragments(atomB), &
                & natoms,DistanceTable,PairFragment,grad)
           flops_slaves = PairFragment%flops_slaves
           tottime = PairFragment%slavetime ! time used by all local slaves
           fragenergy=PairFragment%energies

           call copy_fragment_info_job(PairFragment,singlejob)
           ! Free pair
           call atomic_fragment_free(PairFragment)

        end if SingleOrPair


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
  type(ccatom),intent(inout) :: MyFragment
  !> Number of tasks
  integer,intent(inout) :: ntasks
  type(mp2_batch_construction) :: bat
  integer :: MaxActualDimAlpha,nbatchesAlpha
  integer :: MaxActualDimGamma,nbatchesGamma
  integer, pointer :: orb2batchAlpha(:), batchdimAlpha(:), batchsizeAlpha(:), batchindexAlpha(:)
  integer, pointer :: orb2batchGamma(:), batchdimGamma(:), batchsizeGamma(:), batchindexGamma(:)
  integer :: scheme,nocc,nunocc,MinAObatch,iter
  real(realk) :: MemFree

  ! Initialize stuff (just dummy arguments here)
  nullify(orb2batchAlpha)
  nullify(batchdimAlpha)
  nullify(batchsizeAlpha)
  nullify(batchindexAlpha)
  nullify(orb2batchGamma)
  nullify(batchdimGamma)
  nullify(batchsizeGamma)
  nullify(batchindexGamma)

  ! For fragment with local orbitals where we really want to use the fragment-adapted orbitals
  ! we need to set nocc and nvirt equal to the fragment-adapted dimensions
  if(DECinfo%fragadapt .and. (.not. MyFragment%fragmentadapted) ) then
     nocc=MyFragment%noccFA
     nunocc=MyFragment%nunoccFA
  else
     nocc=MyFragment%noccAOS
     nunocc=MyFragment%nunoccAOS
  end if

  ! Determine optimal batchsizes with available memory
  if(DECinfo%ccmodel==1) then ! MP2
     call get_optimal_batch_sizes_for_mp2_integrals(MyFragment,DECinfo%first_order,bat,.false.)
  else  ! CC2 or CCSD
     iter=1
     call determine_maxBatchOrbitalsize(DECinfo%output,MyFragment%MyLsItem%setting,MinAObatch)
     call get_currently_available_memory(MemFree)
     call get_max_batch_sizes(scheme,MyFragment%number_basis,nunocc,nocc,bat%MaxAllowedDimAlpha,&
          & bat%MaxAllowedDimGamma,MinAObatch,DECinfo%manual_batchsizes,iter,MemFree,.true.)
  end if


  ! Get number of gamma batches
  call mem_alloc(orb2batchGamma,MyFragment%number_basis)
  call build_batchesofAOS(DECinfo%output,MyFragment%mylsitem%setting,bat%MaxAllowedDimGamma,&
       & MyFragment%number_basis,MaxActualDimGamma,batchsizeGamma,batchdimGamma,&
       & batchindexGamma,nbatchesGamma,orb2BatchGamma)
  call mem_dealloc(orb2batchGamma)
  call mem_dealloc(batchdimGamma)
  call mem_dealloc(batchsizeGamma)
  call mem_dealloc(batchindexGamma)

  ! Get number of alpha batches
  call mem_alloc(orb2batchAlpha,MyFragment%number_basis)
  call build_batchesofAOS(DECinfo%output,MyFragment%mylsitem%setting,bat%MaxAllowedDimAlpha,&
       & MyFragment%number_basis,MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,&
       & batchindexAlpha,nbatchesAlpha,orb2BatchAlpha)
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

