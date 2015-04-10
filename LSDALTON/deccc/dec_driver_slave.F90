!> @file 
!> MPI routines for main DEC driver.
!> \author Kasper Kristensen

#ifdef VAR_MPI
module dec_driver_slave_module


  use precision
  use lstiming
  use infpar_module
  use lsmpi_type
  use lsmpi_op, only: init_slave_timers, get_slave_timers
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
    integer :: nvirt,nocc,siz,nfrags
    type(lsitem) :: MyLsitem
    type(decfrag),pointer :: AtomicFragments(:),EstAtomicFragments(:)
    type(fullmolecule) :: MyMolecule
    type(decorbital),pointer :: OccOrbitals(:), virtOrbitals(:)
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

    ! Allocate BackGround Buffer
    call BackGroundBufferAllocate()

    ! GET FULL MOLECULAR INFO FROM MASTER
    ! ***********************************

    ! Get very basic dimension information from main master
    call mpi_dec_fullinfo_master_to_slaves_precursor(esti,nocc,nvirt,master)

    ! Allocate arrays needed for fragment calculations
    call mem_alloc(OccOrbitals,nocc)
    call mem_alloc(virtOrbitals,nvirt)

    ! Get remaining full molecular information needed for all fragment calculations
    call mpi_dec_fullinfo_master_to_slaves(nocc,nvirt,&
         & OccOrbitals, virtOrbitals, MyMolecule, MyLsitem)
    nfrags=MyMolecule%nfrags

    ! Get list of which atoms have orbitals assigned
    call mem_alloc(dofrag,nfrags)
    call which_fragments_to_consider(MyMolecule%ncore,nocc,nvirt,&
         & nfrags,OccOrbitals,virtOrbitals,dofrag,MyMolecule%PhantomAtom)

    IF(DECinfo%StressTest)THEN
     call StressTest_mod_dofrag(MyMolecule%nfrags,nocc,nvirt,MyMolecule%ncore,&
          & MyMolecule%DistanceTable,OccOrbitals,virtOrbitals,dofrag,mylsitem)
    ENDIF

    ! Internal control of first order property keywords
    ! (Necessary because these must be false during fragment optimization.)
    ! Better solution should be implemented at some point...
    dens_save = DECinfo%density
    FO_save = DECinfo%first_order
    grad_save = DECinfo%gradient
    DECinfo%density=.false.
    DECinfo%first_order=.false.
    DECinfo%gradient=.false.

    ! Atomic fragments
    call mem_alloc(AtomicFragments,nfrags)
    do i=1,nfrags
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
          call mem_alloc(EstAtomicFragments,nfrags)
          do i=1,nfrags
             call atomic_fragment_nullify(EstAtomicFragments(i))
          end do
          call mpi_bcast_many_fragments(nfrags,dofrag,EstAtomicFragments,MPI_COMM_LSDALTON)          

          ! Receive CC models to use for each pair based on estimates
          call ls_mpibcast(MyMolecule%ccmodel,nfrags,nfrags,master,MPI_COMM_LSDALTON)
       end if
       

       ! Special for step 2: Reset first-order keywords and get optimized fragments
       ! **************************************************************************
       Step2: if(step==2) then

          ! Restore first order keywords for second step
          DECinfo%density=dens_save
          DECinfo%first_order = FO_save
          DECinfo%gradient = grad_save

          ! Get FOT-optimized atomic fragments from master node
          call mpi_bcast_many_fragments(nfrags,dofrag,AtomicFragments,MPI_COMM_LSDALTON)

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
          call dec_half_local_group( print_ = DECinfo%print_small_calc )

          ! Check if the current rank has become a local master (rank=master within local group)
          if(infpar%lg_mynum==master) localslave=.false.

       end do


       ! Receive  fragment jobs from master and carry those out
       ! ======================================================
       if(step==1 .and. esti) then
          ! Calculate estimated pair fragment energies
          call fragments_slave(nfrags,nocc,nvirt,OccOrbitals,&
               & virtOrbitals,MyMolecule,MyLsitem,AtomicFragments,jobs,&
               & EstAtomicFragments=EstAtomicFragments)
       else
          call fragments_slave(nfrags,nocc,nvirt,OccOrbitals,&
               & virtOrbitals,MyMolecule,MyLsitem,AtomicFragments,jobs)
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
       infpar%lg_comm  = MPI_COMM_LSDALTON
       infpar%lg_mynum = infpar%mynum

       ! Clean up estimated fragments used in step 1 and receive CC models to use for all pairs
       CleanupAndUpdateCCmodel: if(step==1 .and. esti) then
          do i=1,nfrags
             if(dofrag(i)) then
                call atomic_fragment_free_simple(EstAtomicFragments(i))
             end if
          end do
          call mem_dealloc(EstAtomicFragments)

          ! Receive CC models and pair FOT levels to use for each pair based on estimates
          call ls_mpibcast(MyMolecule%ccmodel,nfrags,nfrags,master,MPI_COMM_LSDALTON)
          call ls_mpibcast(MyMolecule%pairfotlevel,nfrags,nfrags,master,MPI_COMM_LSDALTON)

       end if CleanupAndUpdateCCmodel


    end do StepLoop


    ! Free stuff
    ! **********
    do i=1,nfrags
       if(dofrag(i)) then
          call atomic_fragment_free_simple(AtomicFragments(i))
       end if
    end do
    call mem_dealloc(AtomicFragments)
    do i=1,nOcc
       call orbital_free(OccOrbitals(i))
    end do
    do i=1,nvirt
       call orbital_free(virtOrbitals(i))
    end do
    call mem_dealloc(OccOrbitals)
    call mem_dealloc(virtOrbitals)
    call mem_dealloc(dofrag)
    call ls_free(MyLsitem)

    call molecule_finalize(MyMolecule,.false.)
    ! Deallocate BackGround Buffer
    call BackGroundBufferDeallocate()

  end subroutine main_fragment_driver_slave

!> \brief allocate the background buffer framework
!> \author Thomas Kjaergaard
!> \date April 2015
  subroutine BackGroundBufferAllocate()
    implicit none
    logical :: use_bg_buf
    real(realk) :: bytes_to_alloc
    use_bg_buf = .FALSE.
#ifdef VAR_MPI
    IF(DECinfo%use_bg_buffer)use_bg_buf = .TRUE.
#endif
    IF(use_bg_buf)THEN
       !DECinfo%memory is in GB we need it in bytes 
       IF(DECinfo%PL.GT.0)THEN
          print*,'Background buffer initilizes with ',DECinfo%bg_memory,' GB'
          !WRITE(DECinfo%output,'(A,F10.2,A)')'Background buffer initilizes with ',DECinfo%bg_memory,' GB'
          IF(DECinfo%bg_memory.LT.0.OR.DECinfo%bg_memory.GT.DECinfo%memory)THEN
             print*,'DECinfo%bg_memory=',DECinfo%bg_memory
             print*,'DECinfo%memory   =',DECinfo%memory
             call lsquit('.BGMEMORY not set correctly',-1)
          ENDIF
          bytes_to_alloc = DECinfo%bg_memory*1.0E+9_realk 
          call mem_init_background_alloc(bytes_to_alloc)     
       ENDIF
    ENDIF
  end subroutine BackGroundBufferAllocate

!> \brief deallocate the background buffer framework
!> \author Thomas Kjaergaard
!> \date April 2015
  subroutine BackGroundBufferDeallocate()
    implicit none    
    IF(mem_is_background_buf_init())THEN
       call mem_free_background_alloc()
    ENDIF
  end subroutine BackGroundBufferDeallocate

!> \brief For each local master: Carry out  fragment calculations (both single and pair)
!> for those  fragments requested by main master.
!> \author Kasper Kristensen
!> \date May 2012
  subroutine fragments_slave(nfrags,nocc,nvirt,OccOrbitals,&
       & virtOrbitals,MyMolecule,MyLsitem,AtomicFragments,jobs,EstAtomicFragments)

    implicit none

    !> Number of fragments (# atoms for atom-based DEC, #occ orbitals for DECCO)
    integer,intent(in) :: nfrags
    !> Number of occupied orbitals in the molecule
    integer,intent(in) :: nocc
    !> Number of virtupied orbitals in the molecule
    integer,intent(in) :: nvirt
    !> Occupied orbitals, DEC format
    type(decorbital),intent(in) :: OccOrbitals(nocc)
    !> virtupied orbitals, DEC format
    type(decorbital),intent(in) :: virtOrbitals(nvirt)
    !> Full molecular info (not changed at output)
    type(fullmolecule),intent(inout) :: MyMolecule
    !> LS item structure
    type(lsitem),intent(inout) :: MyLsitem
    !> Atomic fragments
    type(decfrag),dimension(nfrags),intent(inout) :: AtomicFragments
    !> Job list for  fragments
    type(joblist),intent(in) :: jobs
    !> Estimated atomic fragments
    type(decfrag),dimension(nfrags),intent(inout),optional :: EstAtomicFragments
    type(joblist) :: singlejob
    type(decfrag) :: PairFragment
    type(decfrag) :: MyFragment
    integer :: job,atomA,atomB,i,slavejob,ntasks,no,nv,mymodel
    real(realk) :: flops_slaves,gpu_flops_slaves
    type(mp2grad) :: grad
    logical :: morejobs, continuedivide, split,only_update
    real(realk) :: fragenergy(ndecenergies),tottime
    real(realk) :: t1cpu, t2cpu, t1wall, t2wall
    real(realk) :: tot_work_time, tot_comm_time, tot_idle_time
    real(realk) :: test_work_time, test_comm_time, test_idle_time, test_master, testtime
    real(realk) :: t1cpuacc, t2cpuacc, t1wallacc, t2wallacc, mwork, midle, mcomm
    real(realk) :: flops,gpu_flops
    real(realk), pointer :: slave_times(:)
    integer(kind=ls_mpik) :: master
    logical :: dividegroup
    master = 0
    fragenergy=0.0_realk
    only_update=.true.
    slave_times => null() 

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
    continuedivide=.true.

    ! Send empty job in first round
    job=0
    morejobs=.true.
    ! Job list containing only a single job which will be overwritten in each new fragment calculation
    call init_joblist(1,singlejob)


    AskForJob: do while(morejobs)
      


       ! Send finished job to master
       ! ***************************
       call time_start_phase(PHASE_COMM)
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
       call time_start_phase(PHASE_WORK)

       ! Carry out fragment optimization (job>0), or finish if all jobs are done (job=-1)
       ! ********************************************************************************
       DoJob: IF(job==-1) then
          call LSTIMER('START',t2cpuacc,t2wallacc,DECinfo%output)
          if(t2wallacc-t1wallacc > 5.0E0_realk.or.DECinfo%PL>2)then
             print '(X,a,i5,a,g14.6)', 'Slave ', infpar%mynum, ' exits  fragments. Time: ',&
               & t2wallacc-t1wallacc
          endif
          morejobs=.false.
       else

          ! Timing and flops
          call time_start_phase(PHASE_WORK, twall = t1wall, swwork = mwork, swcomm = mcomm, swidle = midle )
          call start_flop_counter()

          atomA=jobs%atom1(job)
          atomB=jobs%atom2(job)

          ! Find number of tasks required by integral program
          ! *************************************************
          if(atomA==atomB) then ! single fragment
             FragoptCheck2: if(.not. jobs%dofragopt(job)) then
                ! Init basis for fragment
                call atomic_fragment_init_basis_part(nvirt, nocc, OccOrbitals,&
                     & virtOrbitals,MyMolecule,mylsitem,AtomicFragments(atomA))

                call get_number_of_integral_tasks_for_mpi(AtomicFragments(atomA),ntasks)
                no = AtomicFragments(atomA)%noccAOS
                nv = AtomicFragments(atomA)%nvirtAOS
                mymodel = AtomicFragments(atomA)%ccmodel
  
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
                     & nvirt, nocc, nfrags,OccOrbitals,virtOrbitals, &
                     & MyMolecule,mylsitem,.true.,PairFragment,esti=.true.)
             else
                call merged_fragment_init(AtomicFragments(atomA), AtomicFragments(atomB),&
                     & nvirt, nocc, nfrags,OccOrbitals,virtOrbitals, &
                     & MyMolecule,mylsitem,.true.,PairFragment)
             end if
             call get_number_of_integral_tasks_for_mpi(PairFragment,ntasks)
             no = PairFragment%noccAOS
             nv = PairFragment%nvirtAOS
             mymodel = PairFragment%ccmodel

          end if

          ! Divide local group?
          ! (see details in divideMPIgroup function)
          split=.false.
          continuedivide=.true.
          ! Special case: Never divide for fragment optimization
          call time_start_phase(PHASE_COMM)
          if(jobs%dofragopt(job)) continuedivide=.false.

         

          DoDivide: do while(continuedivide)
             ! Divide group?
             dividegroup = divideMPIgroup(ntasks,no,nv,mymodel)

             if( dividegroup ) then

                if(DECinfo%print_small_calc) print '(a,4i8)', 'DIVIDE! Job, tasks, mynode, #nodes: ', &
                     & job,ntasks,infpar%mynum,infpar%lg_nodtot

                ! Kick slaves out of local group
                call ls_mpibcast(slavejob,master,infpar%lg_comm)

                ! Divide local group into two smaller groups
                call dec_half_local_group
                split=.true.
             else
                continuedivide=.false.
             end if
          end do DoDivide

          call time_start_phase(PHASE_WORK)
          !print '(a,i8,a,i8)', 'Slave ', infpar%mynum, ' will do  job ', job

          ! Communicator in setting may have changed due to division of local group
          ! Ensure that the correct local communicator is used.
          if(split) then
             if(atomA==atomB) then
                AtomicFragments(atomA)%mylsitem%setting%comm = infpar%lg_comm
             else
                PairFragment%mylsitem%setting%comm = infpar%lg_comm
             end if
          end if


          call init_slave_timers(slave_times,infpar%lg_comm)

          ! RUN FRAGMENT CALCULATION
          ! ************************

          FragoptCheck3: if(jobs%dofragopt(job)) then
             ! Sanity check
             if(atomA/=atomB) then
                call lsquit('fragments_slave: Fragment optimization requested for pair fragment!',-1)
             end if

             write(*, '(1X,a,i5,a,i6,a,i15,a,i8)') 'Slave ',infpar%mynum,' will do job: ', &
                & job, ' of size ', jobs%jobsize(job),&
                &  ' single fragment optimization: ', atomA

             ! Fragment optimization
             call optimize_atomic_fragment(atomA,MyFragment,MyMolecule%nfrags, &
                & OccOrbitals,nOcc,virtOrbitals,nvirt, MyMolecule,mylsitem,.true.)


             flops_slaves = MyFragment%flops_slaves
             gpu_flops_slaves = MyFragment%gpu_flops_slaves
             !   tottime = MyFragment%slavetime_work(MyFragment%ccmodel) + &
             !      & MyFragment%slavetime_comm(MyFragment%ccmodel) + &
             !      & MyFragment%slavetime_idle(MyFragment%ccmodel) 
             fragenergy = MyFragment%energies

             call copy_fragment_info_job(MyFragment,singlejob)

          else


             SingleOrPair: if( atomA == atomB ) then ! single  fragment
                write(*, '(1X,a,i5,a,i6,a,i15,a,i4,a,i4,a,i4,a,i6)') 'Slave ',infpar%mynum,'will do job: ', &
                   &job, ' of size ', jobs%jobsize(job),&
                   &  ' with #o',AtomicFragments(atomA)%noccAOS,' #v', AtomicFragments(atomA)%nvirtAOS,&
                   &' #b',AtomicFragments(atomA)%nbasis,' single fragment: ', atomA

                call atomic_driver(MyMolecule,mylsitem,OccOrbitals,virtOrbitals,&
                   & AtomicFragments(atomA),grad=grad)

                flops_slaves = AtomicFragments(atomA)%flops_slaves
                gpu_flops_slaves = AtomicFragments(atomA)%gpu_flops_slaves
                !tottime = AtomicFragments(atomA)%slavetime_work(AtomicFragments(atomA)%ccmodel) + &
                !   & AtomicFragments(atomA)%slavetime_comm(AtomicFragments(atomA)%ccmodel) + &
                !   & AtomicFragments(atomA)%slavetime_idle(AtomicFragments(atomA)%ccmodel) 
                fragenergy = AtomicFragments(atomA)%energies

                call copy_fragment_info_job(AtomicFragments(atomA),singlejob)
                call atomic_fragment_free_basis_info(AtomicFragments(atomA))

             else ! pair  fragment

                if(jobs%esti(job)) then
                   ! Estimated pair fragment
                   call pair_driver(EstAtomicFragments(atomA), EstAtomicFragments(atomB), &
                      & PairFragment,grad)

                   write(*, '(1X,a,i5,a,i6,a,i15,a,i4,a,i4,a,i4,a,i6,i6)')'Slave ',infpar%mynum, &
                      &' will do job: ', job, ' of size ', jobs%jobsize(job),&
                      &  ' with #o',PairFragment%noccAOS,' #v', PairFragment%nvirtAOS,&
                      &' #b',PairFragment%nbasis,' pair   estimate: ', atomA,atomB
                else
                   ! Pair fragment according to FOT precision
                   call pair_driver(AtomicFragments(atomA), AtomicFragments(atomB), &
                      & PairFragment,grad)

                   write(*, '(1X,a,i5,a,i6,a,i15,a,i4,a,i4,a,i4,a,i6,i6)')'Slave ',infpar%mynum, &
                      &' will do job: ', job, ' of size ', jobs%jobsize(job),&
                      &  ' with #o',PairFragment%noccAOS,' #v', PairFragment%nvirtAOS,&
                      &' #b',PairFragment%nbasis,' pair   fragment: ', atomA,atomB

                end if

                flops_slaves = PairFragment%flops_slaves
                gpu_flops_slaves = PairFragment%gpu_flops_slaves
                !tottime = PairFragment%slavetime_work(PairFragment%ccmodel) + &
                !   & PairFragment%slavetime_comm(PairFragment%ccmodel) + &
                !   & PairFragment%slavetime_idle(PairFragment%ccmodel) 
                fragenergy=PairFragment%energies

                call copy_fragment_info_job(PairFragment,singlejob)
                ! Free pair
                call atomic_fragment_free(PairFragment)

             end if SingleOrPair

          end if FragoptCheck3




          ! Set fragment job info
          ! *********************

          !TIME
          call get_slave_timers(slave_times,infpar%lg_comm)
          call time_start_phase(PHASE_WORK, twall = t2wall , dwwork = mwork, dwcomm = mcomm, dwidle = midle )
          !update the times the slaves have spent in the calclulation with the
          !master time-> this gives rise to the time spent in the initialization
          !and finalization of the fragment calculation
          tot_work_time  = mwork
          tot_comm_time  = mcomm
          tot_idle_time  = midle
          test_work_time = slave_times(PHASE_WORK_IDX)
          test_comm_time = slave_times(PHASE_COMM_IDX)
          test_idle_time = slave_times(PHASE_IDLE_IDX)
          test_master    = test_work_time + test_comm_time + test_idle_time 
          do i = 2, infpar%lg_nodtot
             tot_work_time  = tot_work_time  + slave_times((i-1)*nphases+PHASE_WORK_IDX)
             tot_comm_time  = tot_comm_time  + slave_times((i-1)*nphases+PHASE_COMM_IDX)
             tot_idle_time  = tot_idle_time  + slave_times((i-1)*nphases+PHASE_IDLE_IDX)
             test_work_time = test_work_time + slave_times((i-1)*nphases+PHASE_WORK_IDX)
             test_comm_time = test_comm_time + slave_times((i-1)*nphases+PHASE_COMM_IDX)
             test_idle_time = test_idle_time + slave_times((i-1)*nphases+PHASE_IDLE_IDX)
          enddo
          tottime  = tot_work_time  + tot_comm_time  + tot_idle_time
          testtime = test_work_time + test_comm_time + test_idle_time

          singlejob%LMtime(1)  = t2wall - t1wall  ! wall time used by local master
          singlejob%nslaves(1) = infpar%lg_nodtot ! Sizes of local slot (local master + local slaves)
          singlejob%workt(1)   = tot_work_time    ! collective work time
          singlejob%commt(1)   = tot_comm_time    ! collective comm time 
          singlejob%idlet(1)   = tot_idle_time    ! collective idle time

          !check if timings are as expected, within 1% at least
          if( (abs(mwork + mcomm + midle - singlejob%LMtime(1)) > 1.0E-2_realk*singlejob%LMtime(1)).or.&
             &( abs((testtime/singlejob%nslaves(1)) - test_master)>1.0E-2_realk*test_master) )then
             print *,"WARNING(fragments_slave):timing for fragment weird",mwork,mcomm,&
                &midle,singlejob%LMtime(1),testtime,test_master,singlejob%nslaves(1)
          endif

          !FLOPS
          call end_flop_counter(flops,gpu_flops) ! flops for local master
          singlejob%flops(1)     = flops + flops_slaves  ! FLOPS for local master + local slaves
          singlejob%gpu_flops(1)     = gpu_flops + gpu_flops_slaves  ! GPU FLOPS for local master + local slaves
          singlejob%jobsdone(1)  = .true.
          singlejob%esti(1)      = jobs%esti(job)
          singlejob%dofragopt(1) = jobs%dofragopt(job)

          print '(X,a,i5,a,i8,g14.6)', 'Slave ', infpar%mynum, ' is done with  job/time ', &
               & job, singlejob%LMtime(1)

          call mem_dealloc(slave_times)
          slave_times => null() 


       end if DoJob
   

    end do AskForJob


    call free_joblist(singlejob)

  end subroutine fragments_slave

  !> \brief Set logical defining whether MPI group should be divided or not
  !> \author Kasper Kristensen
  !> \date March 2015
  function divideMPIgroup(ntasks,no,nv,ccmodel) result(dividegroup)
    implicit none
    !> Number of tasks (not used for RI-MP2)
    integer,intent(in) :: ntasks
    !> Number of occupied/virtual AOS orbitals (only used for RI-MP2)
    integer,intent(in) :: no,nv
    !> CC model 
    integer,intent(in) :: ccmodel
    logical :: dividegroup
    real(realk) :: rifactor

    ! Divide local group into two smaller groups if:
    ! number of tasks < (#nodes in local group)*(magic factor)
    !
    ! The magic factor (DECinfo%MPIsplit) should in general be larger than 1, 
    ! and it can be changed by .MPIsplit keyword.
    ! 
    ! Special case for RI MP2, divide if:
    ! #nodes > nocc*nvirt/DECinfo%RIMPIsplit
    !
    ! Thus, we obtain optimal MPI performance for the individual fragment
    ! calculations if the number of tasks is significantly larger
    ! than the number of nodes in the local group.
    ! Also check that the local group is larger than just 1 node.
    !

    dividegroup = .false.
    WHichModel: if(ccmodel==MODEL_RIMP2) then 

       ! RIMP2
       rifactor = real(no)*real(nv)/real(DECinfo%RIMPIsplit)
       if(real(infpar%lg_nodtot) > rifactor ) then
          dividegroup=.true.
       end if

    else

       ! Not RIMP2
       if( (ntasks .LE. infpar%lg_nodtot*DECinfo%MPIsplit) ) then
          dividegroup=.true.
       end if

    end if WHichModel


    ! Never split if there is already just one node in the group.
    if(infpar%lg_nodtot==1) dividegroup=.false.

    ! Never split if MPIsplit is set to zero (even for RIMP2)
    if(DECinfo%MPIsplit==0) dividegroup=.false.

  end function divideMPIgroup

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
  real(realk) :: memoryneeded

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

  !MODIFY FOR NEW MODEL

  ! Determine optimal batchsizes with available memory
  if(MyFragment%ccmodel==MODEL_MP2) then ! MP2
     call get_optimal_batch_sizes_for_mp2_integrals(MyFragment,DECinfo%first_order,bat,.false.,.false.,memoryneeded)
  elseif(MyFragment%ccmodel==MODEL_RIMP2) then ! RIMP2
     !do nothing ntasks is not used anyway
  elseif(MyFragment%ccmodel==MODEL_LSTHCRIMP2) then ! LS-THC-RIMP2
     !do nothing
  else  ! CC2 or CCSD
     mpi_split = .true.
     call wrapper_get_ccsd_batch_sizes(MyFragment,bat,mpi_split,ntasks)
     ! In case of MO-based algorithm the number of tasks depends on the
     ! number of MO batches and is directly returned by the routine.
     if (ntasks>0) return
  end if

  if (.not. MyFragment%ccmodel==MODEL_RIMP2) then
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
  end if

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

