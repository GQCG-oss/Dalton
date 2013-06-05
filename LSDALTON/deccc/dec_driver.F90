!> @file
!> Main DEC fragment driver.
!> \author Marcin Ziolkowski and Kasper Kristensen

module dec_driver_module

  use fundamental
  use memory_handling
  use precision
  use lstiming!, only: lstimer
  use typedeftype!, only: lsitem
  use matrix_module!, only:matrix
  use matrix_operations!, only: mat_init, mat_zero,mat_free
  use dec_typedef_module


  ! DEC DEPENDENCIES (within deccc directory) 
  ! *****************************************
  use dec_fragment_utils
  use array2_simple_operations !,only: array2_init, array2_free, array2_copy
  use orbital_operations!,only: get_number_of_orbitals_per_atom
  use full_molecule!,only:molecule_get_fock
  use ccsd_module
  use atomic_fragment_operations!,only: write_fragment_info_to_output,&
 !      & get_fragmentt1_aosaos_from_full,merged_fragment_init,&
 !      & redo_fragment_calculations,extract_specific_fragmentt1,&
 !      & update_full_t1_from_pair_frag,update_full_t1_from_atomic_frag,&
 !      & atomic_fragment_init_basis_part,&
 !      & atomic_fragment_nullify,add_fragment_to_file,restart_atomic_fragments_from_file,&
 !      & estimate_atomic_fragment_sizes,free_joblist,init_joblist,put_job_into_joblist
  use mp2_gradient_module!,only:free_mp2grad,&
!       & init_fullmp2grad,update_full_mp2gradient,get_mp2gradient_main,free_fullmp2grad,&
!       & write_gradient_and_energies_for_restart,read_gradient_and_energies_for_restart
  use fragment_energy_module

#ifdef VAR_MPI
  use infpar_module
  use dec_driver_slave_module
#endif

public:: DEC_wrapper,main_fragment_driver
private

contains

  !> \brief Main DEC driver
  !> \author Marcin Ziokowski and Kasper Kristensen
  !> \date November 2010
  subroutine DEC_wrapper(MyMolecule,mylsitem,D)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(inout) :: MyMolecule
    !> LS item
    type(lsitem), intent(inout) :: mylsitem
    !> Density matrix (only stored on master node, not needed for fragment calculations)
    type(matrix),intent(in) :: D
    type(ccorbital), pointer :: OccOrbitals(:)
    type(ccorbital), pointer :: UnoccOrbitals(:)
    real(realk), pointer :: DistanceTable(:,:)
    real(realk) :: Ecorr, Ehf, Eerr
    real(realk), pointer :: mp2gradient(:,:)
    integer :: nBasis,nOcc,nUnocc,nAtoms,i


    ! Print DEC info
    Eerr = 0.0_realk
    call print_dec_info()

    nOcc = MyMolecule%numocc
    nUnocc = MyMolecule%numvirt
    nBasis = MyMolecule%nbasis
    nAtoms = MyMolecule%natoms

    ! -- Calculate distance matrix
    call mem_alloc(DistanceTable,nAtoms,nAtoms)
    DistanceTable=0.0E0_realk
    call GetDistances(DistanceTable,nAtoms,mylsitem,DECinfo%output) ! distances in atomic units


    ! -- Analyze basis and create orbitals
    call mem_alloc(OccOrbitals,nOcc)
    call mem_alloc(UnoccOrbitals,nUnocc)
    call GenerateOrbitals_driver(MyMolecule,mylsitem,nocc,nunocc,natoms, &
         & OccOrbitals, UnoccOrbitals, DistanceTable)


    ! Calculation type
    ! ****************
    ! Possibilities:
    ! * Single atomic fragment calculation
    ! * Single atomic pair fragment calculation
    ! * Simulate full calculation and use full molecule for all fragments
    ! * Calculate all atomic fragment and pairs to get total correlation energy

    CalculationType: if(DECinfo%mp2energydebug) then
       ! DEC calculation for full molecule where exact single and pair energies are calculated.
       ! Only for testing and only for MP2
       call Full_DEC_calculation(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals, &
            & natoms,nocc,nunocc,DistanceTable,Ecorr)
       Ehf = get_HF_energy_fullmolecule(MyMolecule,Mylsitem,D) 
       call print_total_energy_summary(EHF,Ecorr,Eerr)

    else ! Optimize all fragments and calculate pairs to get full correlation energy

       call mem_alloc(mp2gradient,3,natoms)
       call main_fragment_driver(MyMolecule,mylsitem,D,&
            &OccOrbitals,UnoccOrbitals, &
            & natoms,nocc,nunocc,DistanceTable,EHF,Ecorr,mp2gradient,Eerr)
       call mem_dealloc(mp2gradient)

    end if CalculationType

    call mem_dealloc(DistanceTable)

    ! Delete orbitals
    do i=1,nOcc
       call orbital_free(OccOrbitals(i))
    end do

    do i=1,nUnocc
       call orbital_free(UnoccOrbitals(i))
    end do
    if(DECinfo%PL>0) write(DECinfo%output,'(a)') 'Orbitals removed'

    call mem_dealloc(OccOrbitals)
    call mem_dealloc(UnoccOrbitals)

    ! Check that file handling went OK
    if(files_opened /= 0) then
       write(DECinfo%output,*)
       write(DECinfo%output,'(1X,a)') 'WARNING: Something went wrong with C file handling!'
       write(DECinfo%output,'(1X,a,i5,a)') 'There are ',&
            & files_opened, ' unclosed files!'
       write(DECinfo%output,*)
    end if

  end subroutine DEC_wrapper



  !> \brief Calculate all atomic fragment energies
  !> and all pair interaction energies, and add them
  !> to get the correlation energy for the full molecule
  !> Also calculate density and/or gradient if requested (only for MP2)
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine main_fragment_driver(MyMolecule,mylsitem,D,&
       & OccOrbitals,UnoccOrbitals, &
       & natoms,nocc,nunocc,DistanceTable,EHF,Ecorr,molgrad,Eerr)

    implicit none
    !> Number of occupied orbitals in full molecule
    integer,intent(in) :: nOcc
    !> Number of unoccupied orbitals in full molecule
    integer,intent(in) :: nUnocc
    !> Number of atoms in full molecule
    integer,intent(in) :: nAtoms
    !> Full molecule info (will be unchanged at output)
    type(fullmolecule), intent(inout) :: MyMolecule
    !> LSDalton info
    type(lsitem), intent(inout) :: mylsitem
    !> HF density matrix
    type(matrix),intent(in) :: D
    !> Occupied orbitals in DEC format (not changed at output, is intent(inout) for MPI purposes)
    type(ccorbital), intent(inout) :: OccOrbitals(nocc)
    !> Unoccupied orbitals in DEC format (not changed at output, is intent(inout) for MPI purposes)
    type(ccorbital), intent(inout) :: UnoccOrbitals(nunocc)
    !> Distances between all atoms (not changed at output, is intent(inout) for MPI purposes)
    real(realk),intent(inout) :: DistanceTable(natoms,natoms)
    !> Hartree-Fock energy
    real(realk),intent(inout) :: EHF
    !> MP2 correlation energy
    real(realk),intent(inout) :: Ecorr
    !> Molecular gradient (only calculated if DECinfo%gradient is set!)
    real(realk),intent(inout) :: molgrad(3,natoms)
    !> Estimated intrinsic energy error of DEC calculation
    real(realk),intent(inout) :: Eerr
    ! Fragment energies for Lagrangian (:,:,1), occupied (:,:,2), and virtual (:,:,3) schemes
    real(realk) :: FragEnergies(natoms,natoms,ndecenergies)
    type(ccatom),pointer :: AtomicFragments(:)
    integer :: i,j,k,dims(2),nbasis,counter
    real(realk) :: energies(ndecenergies)
    logical :: dens_save ! Internal control of MP2 density keyword
    logical :: FO_save  ! Internal control of first order property keyword
    logical :: grad_save  ! Internal control of MP2 gradient keyword
    logical :: dofrag(natoms), fragdone(natoms)
    integer, dimension(natoms) :: nocc_per_atom, nunocc_per_atom
    type(array2) :: t1old,fockt1,t1new
    logical :: redo,post_fragopt_restart
    logical,dimension(natoms) :: whichfrags
    type(mp2grad) :: grad
    type(fullmp2grad) :: fullgrad
    integer :: jobdone,newjob, nworkers, njobs, siz, jobidx
    integer(kind=ls_mpik) :: groupsize
    integer :: af_list(natoms),MPIdatatype
    type(joblist) :: jobs,singlejob
    real(realk) :: tcpu,twall,oldpaircut,newpaircut,tcpu1,tcpu2,twall1,twall2,mastertime
    ! (T) contribution to fragment energies for occupied (:,:,1), and virtual (:,:,2) schemes 
    !> (:,:,3): Occupied E[4] contribution;  (:,:,4): Virtual E[4] contribution
    !> (:,:,5): Occupied E[5] contribution;  (:,:,6): Virtual E[5] contribution
    logical :: morejobs,atomic_fragment_restart
    integer(kind=ls_mpik) :: master,IERR,comm,sender
#ifdef VAR_MPI
    INTEGER(kind=ls_mpik) :: MPISTATUS(MPI_STATUS_SIZE), DUMMYSTAT(MPI_STATUS_SIZE)
#endif
    master=0

    call LSTIMER('START',tcpu,twall,DECinfo%output)



    ! ************************************************************************
    ! *                           Initialize stuff                           *
    ! ************************************************************************

    Eerr=0.0_realk
    whichfrags=.false.
    redo=.false.
    nbasis = MyMolecule%nbasis
    call mem_alloc(AtomicFragments,natoms)
    do i=1,natoms
       call atomic_fragment_nullify(AtomicFragments(i))
    end do

    ! Number of orbitals per atom
    nocc_per_atom =  get_number_of_orbitals_per_atom(OccOrbitals,nocc,natoms)
    nunocc_per_atom =  get_number_of_orbitals_per_atom(UnOccOrbitals,nunocc,natoms)

    ! Which fragments to consider
    dofrag=.true.
    do i=1,natoms
       if( (nocc_per_atom(i)==0) .and. (nunocc_per_atom(i)==0) ) dofrag(i)=.false.
    end do
    njobs = count(dofrag)


    ! Special treatment of singles amplitudes
    ! ***************************************
    ! For CC2 and CCSD, collect singles amplitudes from fragment calculations
    ! to construct approximate singles amplitudes for full molecule
    ! to include singles polarization effects in a simple manner.
    ! Initialize appriximate full singles amplitudes here (order: virt,occ)
    if(DECinfo%SinglesPolari) then
       dims(1) = nunocc
       dims(2) = nocc
       t1old = array2_init_plain(dims)
    end if


    ! Restart option: In case some fragments are already done and stored in atomicfragments.info.
    fragdone=.false.
    post_fragopt_restart=.false.
    if(DECinfo%restart) then
       call restart_atomic_fragments_from_file(natoms,MyMolecule,MyLsitem,OccOrbitals,&
            & UnoccOrbitals,.false.,AtomicFragments,jobs,atomic_fragment_restart)


       AFrestartfileexist: if(atomic_fragment_restart) then
          ! Atomic fragment restart files exist - we restart...
          
          write(DECinfo%output,'(a,2i8)') 'RESTARTING STANDARD FRAGMENTS - jobs to do ', &
               & count(dofrag)-count(jobs%jobsdone)

          ! Restart pair fragments only if all atomic fragments are done AND fragment file exists
          if(count(dofrag) == count(jobs%jobsdone)) then
             post_fragopt_restart = fragment_restart_file_exist(DECinfo%first_order)
             if(post_fragopt_restart) then
                write(DECinfo%output,'(a)') 'Fragment optimization is done, restart remaining fragments'
             else
                write(DECinfo%output,'(a)') 'Fragment optimization is done, but no restart file &
                     & for remaining fragments'
                write(DECinfo%output,'(a)') '--> We will calculate remaining fragments from scratch!'
             end if
          end if

          ! Make list of dimension natoms telling which atomic fragments are already done
          ! (note that number of jobs in job list is equal to count(dofrag) which
          !  in general is smaller than natoms because not all atoms have orbitals assigned.
          !  When all calculations are done dofrag=fragdone).
          do i=1,jobs%njobs
             if(jobs%jobsdone(i)) then  ! job number "i" is done
                ! Job number "i" in job list corresponds to atom number jobs%atom1(i)
                fragdone(jobs%atom1(i)) = .true.
             end if
          end do
       
          ! Sanity checks
          if(count(fragdone) /= count(jobs%jobsdone) ) then
             write(DECinfo%output,*) 'jobsdone / fragdone: ', count(jobs%jobsdone),count(fragdone)
             call lsquit('Main driver: Inconsistency in atomic fragment restart job list 1',-1)
          end if
          if(count(dofrag) /= jobs%njobs) then
             write(DECinfo%output,*) 'dofrag / njobs: ', count(dofrag),jobs%njobs
             call lsquit('Main driver: Inconsistency in atomic fragment restart job list 2',-1)
          end if

          ! Subtract number of fragments already done from job count
          njobs = njobs - count(jobs%jobsdone)
          jobidx = count(jobs%jobsdone) ! counter used for putting new jobs into job list

       else
          ! Atomic fragment restart files do not exist --> calculate from scratch
          call init_joblist(njobs,jobs)
          jobidx=0
       end if AFrestartfileexist

    else ! create new empty job list to be updated
       call init_joblist(njobs,jobs)
       jobidx=0
    end if


    ! ************************************************************************
    ! *                         MPI intialization                            *
    ! ************************************************************************


#ifdef VAR_MPI

    ! Number of workers (slaves) = Number of nodes minus master itself
    nworkers = infpar%nodtot -1
    if(nworkers<1) then
       call lsquit('DEC MPI requires at least two nodes!',-1)
    end if

    ! MPI local group size
    if(DECinfo%MPIgroupsize>0) then ! group size was defined explicitly in input
       groupsize=DECinfo%MPIgroupsize
    else 
       ! If nworkers > njobs, then half the atomic fragments will start immediately
       if(njobs>0) then
          groupsize = ceiling(real(nworkers)/real(njobs))
          groupsize = min(2*groupsize,nworkers)
       else  ! sanity precaution to not divide by zero
          groupsize=nworkers
       end if
    end if

    write(DECinfo%output,*) '*** MPI GROUPSIZE SET TO ', groupsize

    ! Initialize local groups 
    call init_mpi_groups(groupsize,DECinfo%output)

    ! Wake up local masters for fragment jobs
    call ls_mpibcast(DECDRIVER,master,MPI_COMM_LSDALTON)

    ! Send DEC input information to slaves
    siz=sizeof(DECinfo)
    call MPI_BCAST(DECinfo,siz,MPI_CHARACTER,master,MPI_COMM_LSDALTON,ierr)

    ! Pass very basic information on dimensions to local masters (necessary to allocate arrays)
    call ls_mpibcast(natoms,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(nocc,master,MPI_COMM_LSDALTON)
    call ls_mpibcast(nunocc,master,MPI_COMM_LSDALTON)

    ! Pass remaining full molecular info to local masters
    call mpi_dec_fullinfo_master_to_slaves(natoms,nocc,nunocc,DistanceTable,&
         & OccOrbitals, UnoccOrbitals, MyMolecule, MyLsitem)

#else
    nworkers=0   ! master node does all jobs
    siz=1
    call init_joblist(siz,singlejob)
#endif


    ! Internal control of first order property keywords
    ! (Necessary because these must be false during fragment optimization.)
    dens_save = DECinfo%MP2density
    FO_save = DECinfo%first_order
    grad_save = DECinfo%gradient
    DECinfo%MP2density=.false.
    DECinfo%first_order=.false.
    DECinfo%gradient=.false.


    ! Sort atomic fragments according to estimated size
    call estimate_atomic_fragment_sizes(natoms,nocc,nunocc,DistanceTable,&
         & OccOrbitals, UnoccOrbitals, mylsitem,af_list)

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '*** Initiating optimization of atomic fragments ***'
    write(DECinfo%output,*)
    write(DECinfo%output,*)

    call LSTIMER('DEC INIT',tcpu,twall,DECinfo%output)
    call LSTIMER('START',tcpu1,twall1,DECinfo%output)

    ! loop over all atoms
    ! MPI: Extra loop for each local worker to finalize things correctly
    counter=0
    jobdone=0
    DoAtomicFragments: do k=1,natoms+nworkers

       if(k<=natoms) then ! start up new calculation
          i = af_list(k)  ! atom to consider
          ! cycle if no orbitals assigned to atom or fragment is already done
          if( (.not. dofrag(i)) .or. fragdone(i)) cycle
       else
          i=k  ! not a real calculation, used for MPI quit signal
       end if
       ! if i>natoms then this is a dummy calculation used to finalize local MPI masters


       ! ************************************************************************
       ! *               DEC-MPI ATOMIC FRAGMENT OPTIMIZATION SCHEME            *
       ! ***********************************************************************
#ifdef VAR_MPI
       ! Counter used to distinquish real (counter<=njobs) and quit (counter>njobs) calculations
       counter=counter+1

       ! Check for local masters to do a job
       CALL MPI_PROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_LSDALTON,MPISTATUS,IERR)
       ! Receive information about which job was done by local master
       ! (if jobdone=0 it is an empty job used just for starting up the scheme)
       call ls_mpisendrecv(jobdone,MPI_COMM_LSDALTON,MPISTATUS(MPI_SOURCE),master)

       ! If job done is NOT an empty job -> book keeping and receive fragment info
       if(jobdone/=0) then
          write(DECinfo%output,*) 'Received atomic fragment: ', jobdone, &
               & ' -- missing #jobs: ', jobs%njobs-count(jobs%jobsdone)-1
          comm = MPI_COMM_LSDALTON
          sender = MPISTATUS(MPI_SOURCE)
          ! Receive fragment info
          call mpi_send_recv_single_fragment(AtomicFragments(jobdone),comm,&
               & sender,master,singlejob)

          ! Increase job counter and put received job info into big job list
          jobidx = jobidx+1
          call put_job_into_joblist(singlejob,jobidx,jobs)
          ! Done with singlejob for now
          call free_joblist(singlejob)
          ! Save fragment info to file atomicfragments.info
          call add_fragment_to_file(AtomicFragments(jobdone),jobs)
       end if

       ! Send new job task to local master
       if(counter <= njobs) then
          newjob=i
          write(DECinfo%output,*) 'Optimizing single fragment:', newjob
       else  ! send quit signal to local master
          newjob=-1
       end if
       call ls_mpisendrecv(newjob,MPI_COMM_LSDALTON,master,MPISTATUS(MPI_SOURCE))
#else

       ! No MPI: Master nodes does all jobs.
       if(DECinfo%SinglesPolari) then ! singles effects for full molecule
          call optimize_atomic_fragment(i,AtomicFragments(i),nAtoms, &
               & OccOrbitals,nOcc,UnoccOrbitals,nUnocc,DistanceTable, &
               & MyMolecule,mylsitem,.true.,t1full=t1old)
       else
          call optimize_atomic_fragment(i,AtomicFragments(i),nAtoms, &
               & OccOrbitals,nOcc,UnoccOrbitals,nUnocc,DistanceTable, &
               & MyMolecule,mylsitem,.true.)
       end if

       ! Copy fragment information into joblist
       call copy_fragment_info_job(AtomicFragments(i),singlejob)
       ! Increase job counter and put received job info into big job list
       jobidx = jobidx+1
       singlejob%jobsdone(1) = .true.
       call put_job_into_joblist(singlejob,jobidx,jobs)
       ! Save fragment info to file atomicfragments.info
       call add_fragment_to_file(AtomicFragments(i),jobs)

#endif


    end do DoAtomicFragments


    ! Save fragment energies
    FragEnergies=0E0_realk
    do i=1,natoms
       if( dofrag(i) ) then
          do j=1,ndecenergies
             FragEnergies(i,i,j) = AtomicFragments(i)%energies(j)
          end do
       end if
    end do

    ! Now all atomic fragment energies have been calculated and the
    ! fragment information has been stored in AtomicFragments.
    call LSTIMER('START',tcpu2,twall2,DECinfo%output)
    mastertime = twall2-twall1
    call LSTIMER('DEC ATOMFRAG',tcpu,twall,DECinfo%output)

#ifdef VAR_MPI
    ! Print MPI statistics
    call print_MPI_fragment_statistics(jobs,mastertime,'ATOMIC FRAGMENTS')
#else
    call free_joblist(singlejob)
#endif
    ! Done with atomic fragment job list
    call free_joblist(jobs)




    ! ************************************************************************
    ! *             Construct job list for remaining fragments               *
    ! ************************************************************************

    ! Restore first order
    DECinfo%MP2density=dens_save
    DECinfo%first_order = FO_save
    DECinfo%gradient = grad_save

    ! Get job list 
    call create_dec_joblist_driver(MyMolecule,mylsitem,natoms,nocc,nunocc,DistanceTable,&
         &OccOrbitals,UnoccOrbitals,AtomicFragments,dofrag,jobs)
    njobs = jobs%njobs


    ! Zero fragment energies because they are recalculated
    ! (except for simple DEC-MP2 with only energies)
    if(DECinfo%ccmodel/=1 .or. DECinfo%first_order .or. DECinfo%InclFullMolecule ) then
       FragEnergies = 0.0E0_realk
    end if

    ! Initialize matrices used for MP2 gradient or density if requested
    ! Note: For density we use subset of gradient structure.
    if(DECinfo%first_order) then
       call init_fullmp2grad(MyMolecule,fullgrad)
    end if


#ifdef VAR_MPI

    ! Fragment MPI communication
    ! **************************

    ! 1. Communicate list of which atoms are EOS atoms
    call ls_mpibcast(dofrag,natoms,master,MPI_COMM_LSDALTON)

    ! 2. Communicate fragments (single fragments) to slaves
    call mpi_bcast_many_fragments(natoms,dofrag,AtomicFragments,MPI_COMM_LSDALTON)

#endif
    call LSTIMER('DEC JOBLIST',tcpu,twall,DECinfo%output)
    call LSTIMER('START',tcpu1,twall1,DECinfo%output)




    ! ***************************************************************************
    ! *                    SINGLE+PAIR FRAGMENTS                                *
    ! ***************************************************************************
    

    morejobs=.true.

    ! Continue as long as there are more jobs to be done
    ! (we might need to add more pairs than in original job list to adapt to precision)
    MorePairs: do while(morejobs)


       redo=.true.
       counter=0
       MainFragLoop: do while(redo)

          ! ***************************************************************************
          ! *               Get full molecular T1-transformed Fock matrix             *
          ! ***************************************************************************
          ! For CC2 and CCSD, use full singles amplitude array (constructed from
          ! fragment calculations) to calculate appriximate T1-transformed Fock matrix
          ! which effectively describes singles polarization effects in a simple manner.
          if(DECinfo%SinglesPolari) then

             fockt1 = array2_init([nbasis,nbasis])
             ! Get T1 transformed Fock matrix in AO basis
             call fullmolecular_get_AOt1Fock(mylsitem,MyMolecule,t1old,fockt1)
             
             ! Set fock matrix associated with fullmolecule structure
             ! equal to T1 transformed Fock matrix
             MyMolecule%fock(1:nbasis,1:nbasis) = fockt1%val(1:nbasis,1:nbasis)
             call array2_free(fockt1)

             ! Init improved full molecular singles constructed from
             ! single and pair fragment calculations.
             dims(1) = nunocc
             dims(2) = nocc
             t1new = array2_init_plain(dims)
          end if

          write(DECinfo%output,*)
          write(DECinfo%output,*) '*** Calculating single and pair interaction energies ***'
          write(DECinfo%output,*)

          ! Calculate all fragments (both single and pairs)
          call fragment_jobs(post_fragopt_restart,nocc,nunocc,natoms,MyMolecule,mylsitem,&
               & OccOrbitals,UnoccOrbitals,DistanceTable,jobs,AtomicFragments,&
               & FragEnergies,fullgrad,t1old,t1new)
          
          ! Compare amplitudes from single fragment calculations
          ! to amplitudes from  fragment (single+pair) calculations
          redo=.false.  ! never redo unless singles correction is requested
          if(DECinfo%SinglesPolari) then
             call redo_fragment_calculations(t1old,t1new,redo)

             ! If we need to redo, then the current t1new takes the role of
             ! t1old in the next round
             if(redo) then
                print *, 'Rerunning ALL  fragment calculations!'
                counter=counter+1
                call array2_copy(t1old,t1new)
                call array2_free(t1new)
             end if
          end if

       end do MainFragLoop


       ! DEC energy control center
       ! *************************

       ! Determine whether more pairs are required, if so determine new pair cutoff
       oldpaircut = DECinfo%pair_distance_threshold
       call dec_energy_control_center(natoms,oldpaircut,dofrag,DistanceTable,&
            & FragEnergies,AtomicFragments,morejobs,newpaircut,Eerr)

       if(morejobs) then ! update pair cutoff
          DECinfo%pair_distance_threshold = newpaircut
          
          ! Expand existing job list (where are jobs are already done) to include additional pair jobs
          ! ------------------------------------------------------------------------------------------
          ! Append new jobs to existing job list
          call expand_joblist_to_include_more_pairs(nocc,nunocc,natoms,DistanceTable,dofrag,&
               & AtomicFragments,oldpaircut,newpaircut,MyMolecule,OccOrbitals,&
               & UnoccOrbitals,jobs)

          ! total number of jobs (already done + new jobs to be done)
          njobs = jobs%njobs 
          
       end if

       if(morejobs .and. DECinfo%SinglesPolari) then
          call lsquit('DEC singles polarization not implemented when job list in increased!',-1)
       end if

#ifdef VAR_MPI
       ! Bcast information telling whether there are more jobs or not
       call ls_mpibcast(morejobs,master,MPI_COMM_LSDALTON)
#endif

       ! If we need to take one more round we should never restart because the 
       ! restart file would not contain the new pairs.
       post_fragopt_restart=.false.
       
    end do MorePairs


    call LSTIMER('START',tcpu2,twall2,DECinfo%output)
    mastertime = twall2-twall1
    call LSTIMER('DEC ALLFRAG',tcpu,twall,DECinfo%output)

#ifdef VAR_MPI
    call MPI_COMM_FREE(infpar%lg_comm,IERR)
    ! Set all MPI groups equal to the world group
    call lsmpi_default_mpi_group
    ! Print MPI statistics
    call print_MPI_fragment_statistics(jobs,mastertime,'ALL FRAGMENTS')
#endif


    ! Total correlation energy 
    do j=1,ndecenergies
       call add_dec_energies(natoms,FragEnergies(:,:,j),dofrag,energies(j))
    end do


    ! Print all fragment energies (should be turned of before release...)
    call print_all_fragment_energies(natoms,FragEnergies,dofrag,AtomicFragments,&
         & DistanceTable,energies)

    ! Set output energy
    select case(DECinfo%ccmodel)
    case(1)
       ! MP2, use occ energy
       Ecorr = energies(1)
    case(2)
       ! CC2, use occ energy
       Ecorr = energies(4)
    case(3)
       ! CCSD, use occ energy
       Ecorr = energies(6)
    case(4)
       ! CCSD(T), use occ energy - of course include both CCSD and (T) contributions
       Ecorr = energies(6) + energies(8)
    end select

    ! If singles polarization was considered, we need to
    ! ensure that the fullmolecule structure contains the standard
    ! (NOT T1 transformed) Fock matrix at output
    if(DECinfo%SinglesPolari) then
       call mem_dealloc(MyMolecule%fock)
       call molecule_get_fock(MyMolecule,mylsitem)  ! Get standard Fock matrix again
       call array2_free(t1old)
       call array2_free(t1new)
    end if
    do i=1,nAtoms
       if(.not. dofrag(i)) cycle
       call atomic_fragment_free_simple(AtomicFragments(i))
    end do
    call mem_dealloc(AtomicFragments)

    ! HF energy
    Ehf = get_HF_energy_fullmolecule(MyMolecule,Mylsitem,D) 

    ! If requested, calculate MP2 gradient and MP2 density (or just density) for
    ! full molecule from fragment contributions
    if(DECinfo%SkipFull) then
       write(DECinfo%output,*) 'WARNING: SKIPPING FULL MOLECULAR PART &
            & OF DENSITY/GRADIENT CALCULATION!'
    else
       if(DECinfo%first_order) then
          fullgrad%EHF = EHF
          ! Calculate MP2 density/gradient and save MP2 density matrices to file
          call get_mp2gradient_main(MyMolecule,mylsitem,D,molgrad,fullgrad)
       end if
    end if

    if(DECinfo%first_order) then
       call free_fullmp2grad(fullgrad)
    end if

    call free_joblist(jobs)

    ! Print short summary
    call print_total_energy_summary(EHF,Ecorr,Eerr)

    call LSTIMER('DEC FINAL',tcpu,twall,DECinfo%output)

  end subroutine main_fragment_driver


  !> \brief Print info about DEC calculation.
  !> \author Kasper Kristensen
  !> \date November 2010
  subroutine print_dec_info()
    implicit none

    write(DECinfo%output,'(/,a)') ' ================================================ '
    write(DECinfo%output,'(a)')   '                  DEC-CC driver                   '
    write(DECinfo%output,'(a,/)') ' ================================================ '

    ! print dec input parameters
    write(DECinfo%output,'(/,a)') '--------------------------'
    write(DECinfo%output,'(a)')   '   DEC input parameters   '
    write(DECinfo%output,'(a,/)') '--------------------------'
    write(DECinfo%output,'(a,g15.2)')  'FOT                    = ',DECinfo%FOT
    write(DECinfo%output,'(a,g15.3)')  'Pair distance thresh.  = ',DECinfo%pair_distance_threshold
    write(DECinfo%output,'(a,g15.3)')  'Pair reduction thresh. = ',DECinfo%PairReductionDistance
    write(DECinfo%output,'(a,g15.3)')  'Simple orbital thr.    = ',DECinfo%simple_orbital_threshold
    write(DECinfo%output,'(a,i4)')     'Expansion step size    = ',DECinfo%FragmentExpansionSize
    write(DECinfo%output,'(a,i4)')     'Print level            = ',DECinfo%PL

    ! print cc parameters
    write(DECinfo%output,'(/,a)') '--------------------------'
    write(DECinfo%output,'(a)')   '  Coupled-cluster input  '
    write(DECinfo%output,'(a,/)') '--------------------------'
    write(DECinfo%output,'(a,a)')      'Wave function          = ',DECinfo%cc_models(DECinfo%ccModel)
    write(DECinfo%output,'(a,i4)')     'MaxIter                = ',DECinfo%ccMaxIter
    write(DECinfo%output,'(a,e8.1e2)') 'Convergence thr        = ',DECinfo%ccConvergenceThreshold
    write(DECinfo%output,'(a,l1)')     'Use CROP               = ',DECinfo%use_crop
    write(DECinfo%output,'(a,i4)')     'CROP subspace          = ',DECinfo%ccMaxDIIS
    write(DECinfo%output,'(a,l1)')     'Preconditioner         = ',DECinfo%use_preconditioner
    write(DECinfo%output,'(a,l1)')     'Precond. B             = ',DECinfo%use_preconditioner_in_b
    write(DECinfo%output,'(a,l1)')     'Debug mode             = ',DECinfo%cc_driver_debug

  end subroutine print_dec_info


  !> \brief Print all fragment energies for given CC model.
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine print_all_fragment_energies(natoms,FragEnergies,dofrag,AtomicFragments,&
       & DistanceTable,energies)
    implicit none
    !> Number of atoms in full molecule
    integer,intent(in) :: nAtoms
    ! Fragment energies as listed in ccatom type def "energies"
    real(realk),intent(in) :: FragEnergies(natoms,natoms,ndecenergies)
    !> Which atoms are associated with a fragment?
    logical,intent(in) :: dofrag(natoms)
    !> Atomic fragments
    type(ccatom),intent(inout) :: AtomicFragments(natoms)
    !> Distances between all atoms (not changed at output, is intent(inout) for MPI purposes)
    real(realk),intent(inout) :: DistanceTable(natoms,natoms)
    !> Total DEC energies (sum of frag energies)
    real(realk),intent(in) :: energies(ndecenergies)



    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '============================================================================='



    select case(DECinfo%ccmodel)
    case(1)
       call print_atomic_fragment_energies(natoms,FragEnergies(:,:,1),dofrag,&
            & 'MP2 Lagrangian single energies')
       call print_atomic_fragment_energies(natoms,FragEnergies(:,:,2),dofrag,&
            & 'MP2 occupied single energies')
       call print_atomic_fragment_energies(natoms,FragEnergies(:,:,3),dofrag,&
            & 'MP2 virtual single energies')

       call print_pair_fragment_energies(natoms,FragEnergies(:,:,1),dofrag,AtomicFragments,&
            & DistanceTable, 'MP2 Lagrangian pair energies')
       call print_pair_fragment_energies(natoms,FragEnergies(:,:,2),dofrag,AtomicFragments,&
            & DistanceTable, 'MP2 occupied pair energies')
       call print_pair_fragment_energies(natoms,FragEnergies(:,:,3),dofrag,AtomicFragments,&
            & DistanceTable, 'MP2 virtual pair energies')

       write(DECinfo%output,*)
       write(DECinfo%output,'(1X,a,g20.10)') 'MP2 Lagrangian correlation energy : ', energies(1)
       write(DECinfo%output,'(1X,a,g20.10)') 'MP2 occupied   correlation energy : ', energies(2)
       write(DECinfo%output,'(1X,a,g20.10)') 'MP2 virtual    correlation energy : ', energies(3)
       write(DECinfo%output,*)

    case(2)
       call print_atomic_fragment_energies(natoms,FragEnergies(:,:,4),dofrag,&
            & 'CC2 occupied single energies')
       call print_atomic_fragment_energies(natoms,FragEnergies(:,:,5),dofrag,&
            & 'CC2 virtual single energies')

       call print_pair_fragment_energies(natoms,FragEnergies(:,:,4),dofrag,AtomicFragments,&
            & DistanceTable, 'CC2 occupied pair energies')
       call print_pair_fragment_energies(natoms,FragEnergies(:,:,5),dofrag,AtomicFragments,&
            & DistanceTable, 'CC2 virtual pair energies')

       write(DECinfo%output,*)
       write(DECinfo%output,'(1X,a,g20.10)') 'CC2 occupied   correlation energy : ', energies(4)
       write(DECinfo%output,'(1X,a,g20.10)') 'CC2 virtual    correlation energy : ', energies(5)
       write(DECinfo%output,*)

    case(3)
       if(.not.DECinfo%CCDhack)then
         call print_atomic_fragment_energies(natoms,FragEnergies(:,:,6),dofrag,&
              & 'CCSD occupied single energies')
         call print_atomic_fragment_energies(natoms,FragEnergies(:,:,7),dofrag,&
              & 'CCSD virtual single energies')
        
         call print_pair_fragment_energies(natoms,FragEnergies(:,:,6),dofrag,AtomicFragments,&
              & DistanceTable, 'CCSD occupied pair energies')
         call print_pair_fragment_energies(natoms,FragEnergies(:,:,7),dofrag,AtomicFragments,&
              & DistanceTable, 'CCSD virtual pair energies')
        
         write(DECinfo%output,*)
         write(DECinfo%output,'(1X,a,g20.10)') 'CCSD occupied   correlation energy : ', energies(6)
         write(DECinfo%output,'(1X,a,g20.10)') 'CCSD virtual    correlation energy : ', energies(7)
         write(DECinfo%output,*)
       else
         call print_atomic_fragment_energies(natoms,FragEnergies(:,:,6),dofrag,&
              & 'CCD occupied single energies')
         call print_atomic_fragment_energies(natoms,FragEnergies(:,:,7),dofrag,&
              & 'CCD virtual single energies')
        
         call print_pair_fragment_energies(natoms,FragEnergies(:,:,6),dofrag,AtomicFragments,&
              & DistanceTable, 'CCD occupied pair energies')
         call print_pair_fragment_energies(natoms,FragEnergies(:,:,7),dofrag,AtomicFragments,&
              & DistanceTable, 'CCD virtual pair energies')
        
         write(DECinfo%output,*)
         write(DECinfo%output,'(1X,a,g20.10)') 'CCD occupied   correlation energy : ', energies(6)
         write(DECinfo%output,'(1X,a,g20.10)') 'CCD virtual    correlation energy : ', energies(7)
         write(DECinfo%output,*)
       endif

    case(4)

       call print_atomic_fragment_energies(natoms,FragEnergies(:,:,6),dofrag,&
            & 'CCSD occupied single energies')
       call print_atomic_fragment_energies(natoms,FragEnergies(:,:,7),dofrag,&
            & 'CCSD virtual single energies')

       call print_pair_fragment_energies(natoms,FragEnergies(:,:,6),dofrag,AtomicFragments,&
            & DistanceTable, 'CCSD occupied pair energies')
       call print_pair_fragment_energies(natoms,FragEnergies(:,:,7),dofrag,AtomicFragments,&
            & DistanceTable, 'CCSD virtual pair energies')

       write(DECinfo%output,*)
       write(DECinfo%output,'(1X,a,g20.10)') 'CCSD occupied correlation energy : ', energies(6)
       write(DECinfo%output,'(1X,a,g20.10)') 'CCSD virtual  correlation energy : ', energies(7)
       write(DECinfo%output,*)

       call print_atomic_fragment_energies(natoms,FragEnergies(:,:,8),dofrag,&
            & '(T) occupied single energies')
       call print_atomic_fragment_energies(natoms,FragEnergies(:,:,9),dofrag,&
            & '(T) virtual single energies')
       call print_atomic_fragment_energies(natoms,FragEnergies(:,:,10),dofrag,&
            & '(T) occupied single energies (fourth order)')
       call print_atomic_fragment_energies(natoms,FragEnergies(:,:,11),dofrag,&
            & '(T) virtual single energies (fourth order)')
       call print_atomic_fragment_energies(natoms,FragEnergies(:,:,12),dofrag,&
            & '(T) occupied single energies (fifth order)')
       call print_atomic_fragment_energies(natoms,FragEnergies(:,:,13),dofrag,&
            & '(T) virtual single energies (fifth order)')

       call print_pair_fragment_energies(natoms,FragEnergies(:,:,8),dofrag,AtomicFragments,&
            & DistanceTable, '(T) occupied pair energies')
       call print_pair_fragment_energies(natoms,FragEnergies(:,:,9),dofrag,AtomicFragments,&
            & DistanceTable, '(T) virtual pair energies')
       call print_pair_fragment_energies(natoms,FragEnergies(:,:,10),dofrag,AtomicFragments,&
            & DistanceTable, '(T) occupied pair energies (fourth order)')
       call print_pair_fragment_energies(natoms,FragEnergies(:,:,11),dofrag,AtomicFragments,&
            & DistanceTable, '(T) virtual pair energies (fourth order)')
       call print_pair_fragment_energies(natoms,FragEnergies(:,:,12),dofrag,AtomicFragments,&
            & DistanceTable, '(T) occupied pair energies (fifth order)')
       call print_pair_fragment_energies(natoms,FragEnergies(:,:,13),dofrag,AtomicFragments,&
            & DistanceTable, '(T) virtual pair energies (fifth order)')

       write(DECinfo%output,*)
       write(DECinfo%output,'(1X,a,g20.10)') '(T) occupied correlation energy : ', energies(8)
       write(DECinfo%output,'(1X,a,g20.10)') '(T) virtual  correlation energy : ', energies(9)
       write(DECinfo%output,'(1X,a,g20.10)') '(T) occupied 4th order energy   : ', energies(10)
       write(DECinfo%output,'(1X,a,g20.10)') '(T) virtual  4th order energy   : ', energies(11)
       write(DECinfo%output,'(1X,a,g20.10)') '(T) occupied 5th order energy   : ', energies(12)
       write(DECinfo%output,'(1X,a,g20.10)') '(T) virtual  5th order energy   : ', energies(13)
       write(DECinfo%output,*)
       write(DECinfo%output,'(1X,a,g20.10)') 'Total CCSD(T) occupied correlation energy : ', &
            & energies(6)+energies(8)
       write(DECinfo%output,'(1X,a,g20.10)') 'Total CCSD(T) virtual  correlation energy : ', &
            & energies(7)+energies(9)
       write(DECinfo%output,*)

    case default
       ! MODIFY FOR NEW MODEL
       ! If you implement new model, please print the fragment energies here,
       ! see ccatom type def. to determine the number for your model (e.g. 1,2, and 3 for MP2).
       write(DECinfo%output,*) 'WARNING: print_all_fragment_energies needs implementation &
            & for model: ', DECinfo%ccmodel
    end select

    ! MODIFY FOR NEW CORRECTION
    ! E.g. for F12:
    if(DECInfo%F12) then
       
       print *, "(DEC_driver) Total energy for MP2-F12: ", energies(14)
       write(DECinfo%output,*)
       write(DECinfo%output,'(1X,a,g20.10)') 'MP2F12-V_gr_term occupied correlation energy : ', energies(14)
       write(DECinfo%output,*)       
         
    endif

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '============================================================================='

  end subroutine print_all_fragment_energies


  !> \brief Print atomic fragment energies
  !> \author Kasper Kristensen
  !> \date September 2012
  subroutine print_atomic_fragment_energies(natoms,FragEnergies,dofrag,string)

    implicit none
    !> Number of atoms in full molecule
    integer,intent(in) :: nAtoms
    ! Fragment energies 
    real(realk),intent(in) :: FragEnergies(natoms,natoms)
    !> Which atoms are associated with a fragment?
    logical,intent(in) :: dofrag(natoms)
    !> Character string to print
    character(*),intent(in) :: string
    integer :: i


    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '================================================================='
    write(DECinfo%output,*) trim(string)
    write(DECinfo%output,*) '================================================================='
    write(DECinfo%output,*)
    write(DECinfo%output,*) 'Fragment       Energy'
    do i=1,natoms
       if(.not. dofrag(i)) cycle
       write(DECinfo%output,'(I6,3X,g20.10)') i,FragEnergies(i,i)
    end do


  end subroutine print_atomic_fragment_energies


  !> \brief Print pair fragment energies
  !> \author Kasper Kristensen
  !> \date September 2012
  subroutine print_pair_fragment_energies(natoms,FragEnergies,dofrag,AtomicFragments,&
       & DistanceTable, string)

    implicit none

    !> Number of atoms in full molecule
    integer,intent(in) :: nAtoms
    ! Fragment energies 
    real(realk),intent(in) :: FragEnergies(natoms,natoms)
    !> Which atoms are associated with a fragment?
    logical,intent(in) :: dofrag(natoms)
    !> Atomic fragments
    type(ccatom),intent(inout) :: AtomicFragments(natoms)
    !> Distances between all atoms (a.u.)
    real(realk),intent(in) :: DistanceTable(natoms,natoms)
    !> Character string to print
    character(*),intent(in) :: string
    integer :: i,j
    real(realk) :: pairdist


    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '================================================================='
    write(DECinfo%output,*) trim(string)
    write(DECinfo%output,*) '================================================================='
    write(DECinfo%output,*)
    write(DECinfo%output,'(2X,a)') 'Frag1  Frag2     Dist(Ang)        Energy'

    do i=1,natoms
       do j=i+1,natoms

          ! Skip if no fragment
          if(.not. dofrag(i)) cycle
          if(.not. dofrag(j)) cycle

          ! Pair distance between fragment centers
          pairdist = get_distance_between_fragments(atomicfragments(i),&
               & atomicfragments(j),natoms,DistanceTable)

          DistanceCheck: if(pairdist < DECinfo%pair_distance_threshold ) then
             write(DECinfo%output,'(I6,2X,I6,2X,g14.5,2X,g18.10)') &
                  & i,j,pairdist*bohr_to_angstrom, FragEnergies(i,j)
          end if DistanceCheck

       end do
    end do

  end subroutine print_pair_fragment_energies


  !> \brief Full DEC calculation where exact single and pair energies are calculated. Mainly for testing.
  !> Only implemented for MP2.
  !> \author Kasper Kristensen
  !> \date April 2011
  subroutine Full_DEC_calculation(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals, &
       & natoms,nocc,nunocc, DistanceTable,Ecorr)


    implicit none
    !> Number of atoms in full molecule
    integer,intent(in) :: natoms
    !> Number of occupied orbitals in molecule
    integer,intent(in) :: nocc
    !> Number of unoccupied orbitals in molecule
    integer,intent(in) :: nunocc
    !> Molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> LS Dalton info
    type(lsitem), intent(inout) :: mylsitem
    !> Occupied MOs
    type(ccorbital), intent(in) :: OccOrbitals(nocc)
    !> UnOccupied MOs
    type(ccorbital), intent(in) :: UnoccOrbitals(nunocc)
    !> Distance table with interatomic distances for atoms in molecule
    real(realk),intent(in) :: DistanceTable(natoms,natoms)
    !> Correlation energy
    real(realk),intent(inout) :: Ecorr

    ! Only MP2, Lagrangian scheme
    if(DECinfo%ccmodel/=1) then
       call lsquit('Full DEC calculation only implemented for MP2!',DECinfo%output)
    end if

    call Full_DEC_calculation_Lagrangian(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals, &
         & natoms,nocc,nunocc, DistanceTable,Ecorr)

  end subroutine Full_DEC_calculation





  !> \brief Driver for doing all fragment jobs and updating
  !> the relevant quantities (energy, density, gradient, t1).
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine fragment_jobs(post_fragopt_restart,nocc,nunocc,natoms,MyMolecule,mylsitem,OccOrbitals,&
       & UnoccOrbitals,DistanceTable,jobs,AtomicFragments,FragEnergies,&
       & fullgrad,t1old,t1new)

    implicit none
    !> Restart fragments
    logical,intent(in) :: post_fragopt_restart
    !> Number of occupied orbitals in full molecule
    integer,intent(in) :: nOcc
    !> Number of unoccupied orbitals in full molecule
    integer,intent(in) :: nUnocc
    !> Number of atoms in full molecule
    integer,intent(in) :: nAtoms
    !> Full molecule info (will be unchanged at output)
    type(fullmolecule), intent(inout) :: MyMolecule
    !> LSDalton info
    type(lsitem), intent(inout) :: mylsitem
    !> Occupied orbitals in DEC format 
    type(ccorbital), intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals in DEC format 
    type(ccorbital), intent(in) :: UnoccOrbitals(nunocc)
    !> Distances between atoms
    real(realk),intent(in) :: DistanceTable(natoms,natoms)
    !>  fragment job list
    type(joblist),intent(inout) :: jobs
    !> Atomic fragments
    type(ccatom),dimension(natoms),intent(inout) :: AtomicFragments
    !> Fragment energies, see "energies" in ccatom type def
    real(realk),intent(inout) :: FragEnergies(natoms,natoms,ndecenergies)
    !> MP2 gradient structure (only used if DECinfo%first_order is set)
    type(fullmp2grad),intent(inout) :: fullgrad
    !> Old t1 amplitudes (see main_fragment_driver for details)
    type(array2),intent(in) :: t1old
    !> New t1 amplitudes (see main_fragment_driver for details)
    type(array2),intent(inout) :: t1new
    type(ccatom) :: PairFragment
    integer :: k,atomA,atomB,i,j,counter,jobdone,nworkers,newjob,jobstodo,siz
    type(mp2grad) :: grad
    type(joblist) :: singlejob,jobsold
    real(realk) :: fragenergy(ndecenergies)
    real(realk) :: t1cpu, t2cpu, t1wall, t2wall, dt
    integer(kind=ls_mpik) ::master, sender, groupsize,IERR
    logical :: only_update
#ifdef VAR_MPI
    INTEGER(kind=ls_mpik) :: MPISTATUS(MPI_STATUS_SIZE), DUMMYSTAT(MPI_STATUS_SIZE)
#endif
    master = 0
    fragenergy=0.0_realk
    only_update=.true.


#ifdef VAR_MPI
    ! Number of workers = Number of nodes minus master itself
    nworkers = infpar%nodtot -1
#else
    nworkers=0   ! master node does all jobs
    siz=1
    call init_joblist(siz,singlejob)
#endif

    ! Restart in case some fragments are already done
    if(post_fragopt_restart) then

       ! Copy existing job list (only used for sanity check)
       call copy_joblist(jobs,jobsold)

       ! Free existing job list 
       call free_joblist(jobs)

       ! Overwrite with job list read from file
       if(DECinfo%first_order) then ! density or gradient
          call read_gradient_and_energies_for_restart(natoms,FragEnergies,jobs,fullgrad)
       else
          call read_fragment_energies_for_restart(natoms,FragEnergies,jobs)
       end if

       ! Sanity check for job list read from file
       call fragment_sanity_check(jobs,jobsold)
       ! Done with old job list
       call free_joblist(jobsold)

       ! Only do jobs which are not already finished
       jobstodo = jobs%njobs - count(jobs%jobsdone)
       write(DECinfo%output,'(a,2i8)') 'RESTARTING FRAGMENTS - jobs to do ',  jobstodo
    else
       ! Do all jobs
       jobstodo = jobs%njobs
    end if


    ! Send fragment job list to slaves and redefine MPI groups
#ifdef VAR_MPI

       ! Send fragment job list
       call bcast_post_fragopt_joblist(jobs,MPI_COMM_LSDALTON)
       ! Send pair distance threshold to make sure we are consistent
       call ls_mpibcast(DECinfo%pair_distance_threshold,master,MPI_COMM_LSDALTON)


       ! Reset MPI group structure for post-fragopt fragments
       ! ****************************************************
       ! New groupsize: 64 or number of slaves
       if(DECinfo%MPIgroupsize>0) then ! group size was defined explicitly in input
          groupsize=DECinfo%MPIgroupsize
       else
          groupsize=min(64,nworkers)
       end if
       write(DECinfo%output,'(1X,a,i8)') 'REDEFINE GROUPSIZE TO ', groupsize
       call lsmpi_barrier(MPI_COMM_LSDALTON)
       call MPI_COMM_FREE(infpar%lg_comm,IERR)

       ! Initialize new group
       call init_mpi_groups(groupsize,DECinfo%output)

#endif


    ! Timing for backup
    call LSTIMER('START',t1cpu,t1wall,DECinfo%output)

    counter=0
    JobLoop: do k=1,jobs%njobs+nworkers

       ! Counter used to distinquish real (counter<=jobstodo) and quit (counter>jobstodo) calculations
       counter=counter+1
       
       if(k<=jobs%njobs) then
          if(jobs%jobsdone(k)) cycle  ! job is already done
       end if

       ! *********************************************
       ! *    MPI PARALLELIZATION OF FRAGMENT JOBS   *
       ! *********************************************

#ifdef VAR_MPI

       ! Check for local masters to do a job
       CALL MPI_PROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_LSDALTON,MPISTATUS,IERR)
       ! Receive information about which job was done by local master
       ! (if jobdone=0 it is an empty job used just for starting up the scheme)


       call ls_mpisendrecv(jobdone,MPI_COMM_LSDALTON,MPISTATUS(MPI_SOURCE),master)

       ! Local master "MPISTATUS(MPI_SOURCE)" did job "jobdone"
       ! ******************************************************
       ! Receive stuff from local master and update
       ReceiveJob: if(jobdone>0) then
       write(DECinfo%output,'(a,i8,a,i8,a,i8)') 'Master received job', jobdone, &
            & ' from ', MPISTATUS(MPI_SOURCE), ' -- missing #jobs: ', &
            & jobs%njobs-count(jobs%jobsdone)-1

          sender = MPISTATUS(MPI_SOURCE)

          ! Receive fragment info (and update density or gradient)
          ! ------------------------------------------------------
          if(DECinfo%first_order) then ! gradient
             call nullify_mp2grad(grad)
             call mpi_send_recv_mp2grad(MPI_COMM_LSDALTON,sender,master,&
                  & grad,fragenergy,singlejob,only_update)
             call update_full_mp2gradient(grad,fullgrad)
             call free_mp2grad(grad)
          else ! just energy
             call mpi_send_recv_fragmentenergy(MPI_COMM_LSDALTON,sender,master,&
                  & fragenergy, singlejob)
          end if

          ! Put received job into total job list at position "jobdone"
          call put_job_into_joblist(singlejob,jobdone,jobs)
          ! Done with single job for now
          call free_joblist(singlejob)

          ! Save fragment energy (symmetric for pairs)
          atomA = jobs%atom1(jobdone)
          atomB = jobs%atom2(jobdone)
          do j=1,ndecenergies
             if(atomA==atomB) then
                FragEnergies(atomA,atomA,j) = fragenergy(j)
             else
                FragEnergies(atomA,atomB,j) = fragenergy(j)
                FragEnergies(atomB,atomA,j) = fragenergy(j)
             end if
          end do


          ! Backup files to be able to do restart
          ! *************************************

          ! Time dt since last backup
          call LSTIMER('START',t2cpu,t2wall,DECinfo%output)
          dt = t2wall - t1wall

          ! Backup if time passed is more than DECinfo%TimeBackup or if all jobs are done
          Backup: if( (dt > DECinfo%TimeBackup) .or. all(jobs%jobsdone) ) then

             ! Save info for restart
             if(DECinfo%first_order) then  ! density and/or gradient 
                call write_gradient_and_energies_for_restart(natoms,FragEnergies,jobs,fullgrad)
             else ! just energy
                call write_fragment_energies_for_restart(natoms,FragEnergies,jobs)
             end if

             ! Reset timer
             call LSTIMER('START',t1cpu,t1wall,DECinfo%output)

          end if Backup

       end if ReceiveJob


       ! Send new job task to local master
       if(counter <= jobs%njobs) then
          newjob=k
       else  ! send quit signal to local master
          newjob=-1
       end if
       call ls_mpisendrecv(newjob,MPI_COMM_LSDALTON,master,MPISTATUS(MPI_SOURCE))

#else
! NO MPI --> "Master" does all fragment jobs

       ! Identify job "k"
       atomA = jobs%atom1(k)
       atomB = jobs%atom2(k)
       ! atomA=atomB  : Single fragment job
       ! atomA/=atomB : Pair fragment job

       if(atomA==atomB) then ! single
          print '(1X,a,i8,a,i15,a,i8)', 'Job: ', k, ' of size ', jobs%jobsize(k),&
               &  ' is single fragment: ', atomA

          ! Fragment "atomA" is stored in AtomicFragments(atomA).
          ! However, the basis information has not yet been initialized

          ! Get single fragment energy (and possibly density or gradient)
          ! *******************************************************************
          if(DECinfo%SinglesPolari) then
             call atomic_driver_advanced(nocc,nunocc,&
                  & OccOrbitals,UnoccOrbitals,MyLsitem,MyMolecule,&
                  & AtomicFragments(atomA),grad,t1old=t1old,t1new=t1new)
          else

             ! Init fragment basis information
             call atomic_fragment_init_basis_part(nunocc, nocc, OccOrbitals,&
                  & UnoccOrbitals,MyMolecule,mylsitem,AtomicFragments(atomA))

             ! Call main driver to get energy (and possibly density or gradient)
             call atomic_driver(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals,&
                  & AtomicFragments(atomA),grad=grad)

             ! Free basis info again
             call atomic_fragment_free_basis_info(AtomicFragments(atomA))

          end if

          ! Update energies
          do j=1,ndecenergies
             FragEnergies(atomA,atomA,j) = AtomicFragments(atomA)%energies(j)
          end do

          ! Copy fragment information into joblist
          call copy_fragment_info_job(AtomicFragments(atomA),singlejob)
          singlejob%jobsdone(1) = .true.

       else ! pair calculation

          print '(1X,a,i8,a,i15,a,2i8)', 'Job: ', k, ' of size ', jobs%jobsize(k),&
               &  ' is pair fragment: ', atomA,atomB

          ! Get energy (and possibly density or gradient)
          ! *********************************************

          ! Init pair
          call merged_fragment_init(AtomicFragments(atomA), AtomicFragments(atomB),&
               & nunocc, nocc, natoms,OccOrbitals,UnoccOrbitals, DistanceTable, &
               & MyMolecule,mylsitem,.true.,PairFragment)


          if(DECinfo%SinglesPolari) then
             call pair_driver_singles(natoms,nocc,nunocc,DistanceTable,&
                  & OccOrbitals,UnoccOrbitals,MyLsitem,MyMolecule,&
                  & AtomicFragments(atomA), AtomicFragments(atomB),PairFragment,t1old,t1new)
          else
             call pair_driver(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals,&
                  & AtomicFragments(atomA), AtomicFragments(atomB),&
                  & natoms,DistanceTable,PairFragment,grad)
          end if

          ! Update pair energies
          do j=1,ndecenergies
             fragenergies(atomA,atomB,j) = PairFragment%energies(j)
             fragenergies(atomB,atomA,j) = PairFragment%energies(j)
          end do

          ! Copy fragment information into joblist
          call copy_fragment_info_job(PairFragment,singlejob)
          singlejob%jobsdone(1) = .true.

          call atomic_fragment_free(PairFragment)

       end if


       ! Update stuff (same for single and pair)
       ! ***************************************
       if(DECinfo%first_order) then  ! density and/or gradient
          call update_full_mp2gradient(grad,fullgrad)
          call free_mp2grad(grad)
       end if


       ! Save info for restart
       ! *********************
       ! Put received job into total job list
       call put_job_into_joblist(singlejob,k,jobs)
       if(DECinfo%first_order) then  ! density and/or gradient 
          call write_gradient_and_energies_for_restart(natoms,FragEnergies,jobs,fullgrad)
       else ! just energy
          call write_fragment_energies_for_restart(natoms,FragEnergies,jobs)
       end if

#endif

    end do JobLoop

#ifndef VAR_MPI
    call free_joblist(singlejob)
#endif


  end subroutine fragment_jobs


end module dec_driver_module


