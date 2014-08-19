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

public:: DEC_wrapper
private

contains

  !> \brief Main DEC driver
  !> \author Marcin Ziokowski and Kasper Kristensen
  !> \date November 2010
  subroutine DEC_wrapper(MyMolecule,mylsitem,D,EHF,Ecorr,molgrad,Eerr)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(inout) :: MyMolecule
    !> LS item
    type(lsitem), intent(inout) :: mylsitem
    !> Density matrix (only stored on master node, not needed for fragment calculations)
    type(matrix),intent(in) :: D
    !> Hartree-Fock energy
    real(realk),intent(inout) :: EHF
    !> MP2 correlation energy
    real(realk),intent(inout) :: Ecorr
    !> Molecular gradient (only calculated if DECinfo%gradient is set!)
    real(realk),intent(inout) :: molgrad(3,MyMolecule%natoms)
    !> Estimated energy error
    real(realk),intent(inout) :: Eerr
    type(decorbital), pointer :: OccOrbitals(:)
    type(decorbital), pointer :: UnoccOrbitals(:)
    real(realk),pointer :: FragEnergiesOcc(:,:)
    integer :: nBasis,nOcc,nUnocc,nAtoms,i


    ! Print DEC info
    call print_dec_info()

    nOcc   = MyMolecule%nocc
    nUnocc = MyMolecule%nunocc
    nBasis = MyMolecule%nbasis
    nAtoms = MyMolecule%natoms

    ! -- Analyze basis and create orbitals
    call mem_alloc(OccOrbitals,nOcc)
    call mem_alloc(UnoccOrbitals,nUnocc)
    call GenerateOrbitals_driver(MyMolecule,mylsitem,nocc,nunocc,natoms, &
         & OccOrbitals, UnoccOrbitals)


    ! *************************************************
    ! Optimize all atomic fragments and calculate pairs
    ! *************************************************
    call mem_alloc(FragEnergiesOcc,MyMolecule%natoms,MyMolecule%natoms)
    call main_fragment_driver(MyMolecule,mylsitem,D,&
         &OccOrbitals,UnoccOrbitals, &
         & natoms,nocc,nunocc,EHF,Ecorr,molgrad,Eerr,FragEnergiesOcc)


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
    call mem_dealloc(FragEnergiesOcc)

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
       & natoms,nocc,nunocc,EHF,Ecorr,molgrad,Eerr,FragEnergiesOcc)

    implicit none
    !> Number of occupied orbitals in full molecule (not changed, inout for MPI reasons)
    integer,intent(inout) :: nOcc
    !> Number of unoccupied orbitals in full molecule (not changed, inout for MPI reasons)
    integer,intent(inout) :: nUnocc
    !> Number of atoms in full molecule
    integer,intent(in) :: nAtoms
    !> Full molecule info (will be unchanged at output)
    type(fullmolecule), intent(inout) :: MyMolecule
    !> LSDalton info
    type(lsitem), intent(inout) :: mylsitem
    !> HF density matrix
    type(matrix),intent(in) :: D
    !> Occupied orbitals in DEC format (not changed at output, is intent(inout) for MPI purposes)
    type(decorbital), intent(inout) :: OccOrbitals(nocc)
    !> Unoccupied orbitals in DEC format (not changed at output, is intent(inout) for MPI purposes)
    type(decorbital), intent(inout) :: UnoccOrbitals(nunocc)
    !> Hartree-Fock energy
    real(realk),intent(inout) :: EHF
    !> MP2 correlation energy
    real(realk),intent(inout) :: Ecorr
    !> Molecular gradient (only calculated if DECinfo%gradient is set!)
    real(realk),intent(inout) :: molgrad(3,natoms)
    !> Estimated energy error
    real(realk),intent(inout) :: Eerr
    !> Fragment energies for occupied partitioning scheme
    real(realk),dimension(natoms,natoms),intent(inout) :: FragEnergiesOcc
    logical :: esti
    ! Fragment energies
    !real(realk) :: FragEnergies(natoms,natoms,ndecenergies)
    real(realk),pointer :: FragEnergies(:,:,:) !(natoms,natoms,ndecenergies)
    type(decfrag),pointer :: AtomicFragments(:)
    integer :: i,j,k,dims(2),nbasis,counter
    !real(realk) :: energies(ndecenergies)
    !real(realk) :: Interactionenergies(ndecenergies)
    real(realk),pointer :: energies(:)!(ndecenergies)
    real(realk),pointer :: Interactionenergies(:)!(ndecenergies)
    real(realk) :: InteractionEcorr
    real(realk) :: InteractionEerr
    logical :: dens_save ! Internal control of MP2 density keyword
    logical :: FO_save  ! Internal control of first order property keyword
    logical :: grad_save  ! Internal control of MP2 gradient keyword
    logical :: dofrag(natoms)
    type(array2) :: t1old,fockt1,t1new
    logical :: redo
    type(mp2grad) :: grad
    type(fullmp2grad) :: fullgrad
    integer :: jobdone,newjob, nworkers, siz, jobidx
    integer(kind=ls_mpik) :: groupsize
    integer :: af_list(natoms),MPIdatatype,MyAtom
    type(joblist) :: jobs
    real(realk) :: tcpu,twall,oldpaircut,newpaircut,tcpu1,tcpu2,twall1,twall2,mastertime
    ! (T) contribution to fragment energies for occupied (:,:,1), and virtual (:,:,2) schemes 
    !> (:,:,3): Occupied E[4] contribution;  (:,:,4): Virtual E[4] contribution
    !> (:,:,5): Occupied E[5] contribution;  (:,:,6): Virtual E[5] contribution
    logical :: calcAF,ForcePrintTime
    integer(kind=ls_mpik) :: master,IERR,comm,sender
#ifdef VAR_MPI
    INTEGER(kind=ls_mpik) :: MPISTATUS(MPI_STATUS_SIZE), DUMMYSTAT(MPI_STATUS_SIZE)
#endif
    master=0
    ForcePrintTime = .TRUE.

    call LSTIMER('START',tcpu,twall,DECinfo%output)



    ! ************************************************************************
    ! *                           Initialize stuff                           *
    ! ************************************************************************

    redo=.false.
    nbasis = MyMolecule%nbasis
    call mem_alloc(AtomicFragments,natoms)
    do i=1,natoms
       call atomic_fragment_nullify(AtomicFragments(i))
    end do
    call mem_alloc(FragEnergies,natoms,natoms,ndecenergies)
    FragEnergies=0E0_realk

    ! Find out which atoms have one or more orbitals assigned
    call which_atoms_have_orbitals_assigned(MyMolecule%ncore,nocc,nunocc,natoms,&
         & OccOrbitals,UnoccOrbitals,dofrag,MyMolecule%PhantomAtom)

    if(DECinfo%StressTest)then
       call StressTest_mod_dofrag(MyMolecule%natoms,nocc,nunocc,&
            & MyMolecule%DistanceTable,OccOrbitals, UnoccOrbitals, dofrag, mylsitem)
    endif

    if(DECinfo%PairEstimate .and. count(dofrag)>1) then
       ! Use estimated pair fragments to determine which pair fragments to calculate at the FOT level
       ! (only if there actually are any pair fragments)
       esti = .true.
    else
       esti = .false.
    end if

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




    ! ************************************************************************
    ! *                         MPI intialization                            *
    ! ************************************************************************


#ifdef VAR_MPI

    ! Wake up all slaves for fragment jobs
    call ls_mpibcast(DECDRIVER,master,MPI_COMM_LSDALTON)

    ! Send DEC input information to slaves
    call mpibcast_dec_settings(DECinfo,MPI_COMM_LSDALTON)

    ! Pass very basic information on dimensions to local masters (necessary to allocate arrays)
    call mpi_dec_fullinfo_master_to_slaves_precursor(esti,nocc,nunocc,master)

    ! Pass remaining full molecular info to local masters
    call mpi_dec_fullinfo_master_to_slaves(nocc,nunocc,&
         & OccOrbitals, UnoccOrbitals, MyMolecule, MyLsitem)

#else
    nworkers=0   ! master node does all jobs
    siz=1
#endif


    ! Internal control of first order property keywords
    ! (Necessary because these must be false during fragment optimization.)
    dens_save           = DECinfo%density
    FO_save             = DECinfo%first_order
    grad_save           = DECinfo%gradient
    DECinfo%density  = .false.
    DECinfo%first_order = .false.
    DECinfo%gradient    = .false.
    call LSTIMER('DEC INIT',tcpu,twall,DECinfo%output,ForcePrintTime)


    ! FRAGMENT OPTIMIZATION AND (POSSIBLY) ESTIMATED FRAGMENTS
    ! ********************************************************
    call fragopt_and_estimated_frags(nOcc,nUnocc,OccOrbitals,UnoccOrbitals, &
         & MyMolecule,mylsitem,dofrag,esti,AtomicFragments,FragEnergies)


    !Crash calculation on purpose to test restart option
    IF(DECinfo%CRASHCALC)THEN
       print*,'Calculation was intentionally crashed due to keyword .CRASHCALC'
       print*,'This keyword is only used for debug and testing purposes'
       print*,'We want to be able to test the .RESTART keyword'
       WRITE(DECinfo%output,*)'Calculation was intentionally crashed due to keyword .CRASHCALC'
       WRITE(DECinfo%output,*)'This keyword is only used for debug and testing purposes'
       WRITE(DECinfo%output,*)'We want to be able to test the .RESTART keyword'
       call lsquit('Crashed Calculation due to .CRASHCALC keyword',DECinfo%output)
    ENDIF


    ! Send CC models to use for all pairs based on estimates
    if(esti) then
#ifdef VAR_MPI
       call ls_mpibcast(MyMolecule%ccmodel,natoms,natoms,master,MPI_COMM_LSDALTON)
#endif
    end if

    IF(esti)THEN
       call LSTIMER('DEC Atomic Frags and Estimates',tcpu,twall,DECinfo%output,ForcePrintTime)
    ELSE
       call LSTIMER('DEC Atomic Fragment Calculation',tcpu,twall,DECinfo%output,ForcePrintTime)
    ENDIF

    ! Done with estimates
    esti=.false.
    ! Save fragment energies and set model for atomic fragments appropriately
    do i=1,natoms
       if( dofrag(i) ) then
          do j=1,ndecenergies
             FragEnergies(i,i,j) = AtomicFragments(i)%energies(j)
          end do

          ! If atomic fragments are to be repeated we need to reset original model
          ! for fragments here! For example, if fragment optimization was done
          ! at the MP2 level, then AtomicFragments(i)%ccmodel is MODEL_MP2 now.
          ! However, if the target model is CCSD, then we of course need to use the original
          ! model for the subsequent atomic fragment calculations.
          if(DECinfo%RepeatAF) then
             AtomicFragments(i)%ccmodel = DECinfo%ccmodel
          end if

       end if
    end do

    ! Now all atomic fragment energies have been calculated and the
    ! fragment information has been stored in AtomicFragments.

    IF(DECinfo%InteractionEnergy)THEN
       if(DECinfo%ccmodel.NE.DECinfo%fragopt_red_model) then
          !the energies was calculated with a wrong wave function and
          !should be recalculated but they are not needed for 
          !interaction energies so we set the energy to zero to avoid 
          !confusion. 
          FragEnergies = 0.0E0_realk
       endif
    ENDIF
    ! ************************************************************************
    ! *             Construct job list for remaining fragments               *
    ! ************************************************************************

    ! Restore first order
    DECinfo%density=dens_save
    DECinfo%first_order = FO_save
    DECinfo%gradient = grad_save

    ! Get job list 
    calcAF = DECinfo%RepeatAF
    call create_dec_joblist_driver(calcAF,MyMolecule,mylsitem,natoms,nocc,nunocc,&
         &OccOrbitals,UnoccOrbitals,AtomicFragments,dofrag,.false.,jobs)

    ! Zero fragment energies if they are recalculated
    if(calcAF) then
       FragEnergies = 0.0E0_realk
    end if

    ! Initialize matrices used for MP2 gradient or density if requested
    ! Note: For density we use subset of gradient structure.
    if(DECinfo%first_order) then
       call init_fullmp2grad(MyMolecule,fullgrad)
    end if

#ifdef VAR_MPI

    ! Communicate atomic fragments to slaves
    call mpi_bcast_many_fragments(natoms,dofrag,AtomicFragments,MPI_COMM_LSDALTON)

#endif
    call LSTIMER('DEC JOBLIST',tcpu,twall,DECinfo%output)
    call LSTIMER('START',tcpu1,twall1,DECinfo%output)




    ! ***************************************************************************
    ! *                    SINGLE+PAIR FRAGMENTS                                *
    ! ***************************************************************************


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


       ! ************************************************************************
       ! Calculate pair fragments (and also atomic fragments if DECinfo%RepeatAF)
       ! ************************************************************************

       ! Restart?
       if(DECinfo%DECrestart) then
          if(DECinfo%first_order) then ! density or gradient
             write(DECinfo%output,*) 'Restarting pair fragments - energy and first order prop...'
             call read_gradient_and_energies_for_restart(natoms,FragEnergies,jobs,fullgrad)
          else
             write(DECinfo%output,*) 'Restarting pair fragments - energy...'
             call read_fragment_energies_for_restart(natoms,FragEnergies,jobs,esti)
          end if
       end if
       if(DECinfo%only_n_frag_jobs > 0)then
          jobs%jobsdone  = .false.
          jobs%dofragopt = .false.
       endif

       call fragment_jobs(nocc,nunocc,natoms,MyMolecule,mylsitem,&
            & OccOrbitals,UnoccOrbitals,jobs,AtomicFragments,&
            & FragEnergies,esti,fullgrad=fullgrad,t1old=t1old,t1new=t1new)

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


    ! Plot pair interaction energies using occ. partitioning scheme
    ! *************************************************************
    IF(DECinfo%onlyVirtPart)THEN
       call get_virtfragenergies(natoms,DECinfo%ccmodel,FragEnergies,FragEnergiesOcc)
    ELSE
       call get_occfragenergies(natoms,DECinfo%ccmodel,FragEnergies,FragEnergiesOcc)
    ENDIF
    call plot_pair_energies(natoms,DECinfo%pair_distance_threshold,FragEnergiesOcc,&
         & MyMolecule,dofrag)

    call LSTIMER('START',tcpu2,twall2,DECinfo%output)
    mastertime = twall2-twall1
    IF(DECinfo%RepeatAF)THEN
       call LSTIMER('DEC Atomic and Pair Fragcalc',tcpu,twall,DECinfo%output,ForcePrintTime)
    ELSE
       call LSTIMER('DEC PAIR Fragment calc',tcpu,twall,DECinfo%output,ForcePrintTime)
    ENDIF

#ifdef VAR_MPI
    ! Set all MPI groups equal to the world group
    call lsmpi_default_mpi_group
    ! Print MPI statistics
    if(DECinfo%RepeatAF) then
       call print_MPI_fragment_statistics(jobs,mastertime,'ALL FRAGMENTS')
    else
       call print_MPI_fragment_statistics(jobs,mastertime,'PAIR FRAGMENTS')
    end if
#endif

    call mem_alloc(energies,ndecenergies)
    IF(DECinfo%InteractionEnergy)THEN
       ! Interaction correlation energy 
       do j=1,ndecenergies
          call add_dec_Interactionenergies(natoms,FragEnergies(:,:,j),dofrag,&
               & energies(j),mymolecule%SubSystemIndex)
       end do
    ELSE
       ! Total correlation energy 
       do j=1,ndecenergies
          call add_dec_energies(natoms,FragEnergies(:,:,j),dofrag,energies(j))
       end do
       if(DECinfo%PrintInteractionEnergy)then
         call mem_alloc(Interactionenergies,ndecenergies)
          do j=1,ndecenergies
             call add_dec_Interactionenergies(natoms,FragEnergies(:,:,j),dofrag,&
                  & Interactionenergies(j),mymolecule%SubSystemIndex)
          end do
       endif
    endif

    ! Print all fragment energies
    call print_all_fragment_energies(natoms,FragEnergies,dofrag,&
         & mymolecule%DistanceTable,energies)
     call mem_dealloc(FragEnergies)
    !Obtain The Correlation Energy from the list energies
    InteractionEcorr = 0.0E0_realk
    call ObtainModelEnergyFromEnergies(DECinfo%ccmodel,energies,Ecorr)
    IF(DECinfo%InteractionEnergy)THEN
       InteractionEcorr = Ecorr
    ENDIF
    IF(DECinfo%PrintInteractionEnergy)THEN
       call ObtainModelEnergyFromEnergies(DECinfo%ccmodel,interactionenergies,InteractionEcorr)
    ENDIF
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

    ! Estimate energy error
    call get_estimated_energy_error(natoms,energies,Eerr)
    call mem_dealloc(energies)

    ! Print short summary
    call print_total_energy_summary(EHF,Ecorr,Eerr)
    IF(DECinfo%PrintInteractionEnergy)THEN
       call get_estimated_energy_error(natoms,Interactionenergies,InteractionEerr)
       call print_Interaction_energy(InteractionEcorr,InteractionEerr)   
       call mem_dealloc(Interactionenergies)
    ENDIF
    call LSTIMER('DEC FINAL',tcpu,twall,DECinfo%output,ForcePrintTime)

  end subroutine main_fragment_driver

  subroutine ObtainModelEnergyFromEnergies(ccmodel,energies,Ecorr)
    implicit none
    integer,intent(in)        ::  ccmodel
    real(realk),intent(in)    ::  energies(ndecenergies)
    real(realk),intent(inout) ::  Ecorr
    ! MODIFY FOR NEW MODEL
    ! MODIFY FOR NEW CORRECTION: Add correction to output energy
    ! Set output energy: We choose occupied partitioning scheme energy as general output
    select case(DECinfo%ccmodel)
    case(MODEL_MP2)
       ! MP2, use occ energy
       IF(DECinfo%onlyVirtPart)THEN
          Ecorr = energies(FRAGMODEL_VIRTMP2)
       ELSE
          Ecorr = energies(FRAGMODEL_OCCMP2)
       ENDIF
#ifdef MOD_UNRELEASED
       if(DECinfo%F12) then
          IF(DECinfo%onlyVirtPart)THEN
             Ecorr = energies(FRAGMODEL_MP2f12) + energies(FRAGMODEL_VIRTMP2)
          ELSE
             Ecorr = energies(FRAGMODEL_MP2f12) + energies(FRAGMODEL_OCCMP2)
          ENDIF
       endif
#endif 
    case(MODEL_RPA)
       ! RPA, use occ energy
       IF(DECinfo%onlyVirtPart)THEN
          Ecorr = energies(FRAGMODEL_VIRTRPA)
       ELSE
          Ecorr = energies(FRAGMODEL_OCCRPA)
       ENDIF
    case(MODEL_CC2)
       ! CC2, use occ energy
       IF(DECinfo%onlyVirtPart)THEN
          Ecorr = energies(FRAGMODEL_VIRTCC2)
       ELSE
          Ecorr = energies(FRAGMODEL_OCCCC2)
       ENDIF
    case(MODEL_CCSD)
       ! CCSD, use occ energy
       IF(DECinfo%onlyVirtPart)THEN
          Ecorr = energies(FRAGMODEL_VIRTCCSD)
       ELSE
          Ecorr = energies(FRAGMODEL_OCCCCSD)
       ENDIF
    case(MODEL_CCSDpT)
       ! CCSD(T), use occ energy - of course include both CCSD and (T) contributions
       IF(DECinfo%onlyVirtPart)THEN
          Ecorr = energies(FRAGMODEL_VIRTCCSD) + energies(FRAGMODEL_VIRTpT)
       ELSE
          Ecorr = energies(FRAGMODEL_OCCCCSD) + energies(FRAGMODEL_OCCpT)
       ENDIF
    case default
       write(DECinfo%output,*) 'main_fragment_driver: Needs implementation for model ', DECinfo%ccmodel
       call lsquit('main_fragment_driver: Needs implementation for model!',-1)
    end select
  end subroutine ObtainModelEnergyFromEnergies
!> \brief Print info about DEC calculation.
!> \author Kasper Kristensen
!> \date November 2010
subroutine print_dec_info()
  implicit none
  integer :: LU
  character(len=5) :: LogicString(2)
  LogicString(1) = ' TRUE'
  LogicString(2) = 'FALSE'

  LU = DECinfo%output
  write(LU,'(/,a)') ' ================================================ '
  write(LU,'(a)')   '                  DEC-CC driver                   '
  write(LU,'(a,/)') ' ================================================ '
  
  IF(.NOT.DECinfo%full_molecular_cc)THEN
   ! print dec input parameters
   write(LU,'(/,a)') '--------------------------'
   write(LU,'(a)')   '   DEC input parameters   '
   write(LU,'(a,/)') '--------------------------'
   write(LU,'(a,g15.2)') 'FOT (Fragment Optimization Threshold)               = ',DECinfo%FOT
   write(LU,'(a,A5)')    'Use Pair Estimates to screen pairs                  =           ',&
        & LogicString(Log2It(DECinfo%PairEstimate))
   
   if(DECinfo%PairEstimate) then
        write(LU,'(a,g15.3)') 'Pair distance cutoff threshold (Angstrom)           = ',&
          & DECinfo%pair_distance_threshold*bohr_to_angstrom
        if (DECinfo%orb_based_fragopt) then 
           write(LU,'(a,i5)') 'Use Pair Estimate initialisation number of atom     = ',&
             & DECinfo%estimateInitAtom
        else
           write(LU,'(a,e15.3)') 'Use Pair Estimate initialisation radius (Angstrom)  = ',&
             & DECinfo%estimateINITradius*bohr_to_angstrom
        end if
   endif

   write(LU,'(a,g15.3)') 'Simple orbital thr.                                 = ',&
        & DECinfo%simple_orbital_threshold
   write(LU,'(a,i5)')    'Expansion step size                                 =           ',&
        & DECinfo%Frag_Exp_Size
   write(LU,'(a,A5)')    'Use RI for Expansion step                           =           ',&
        & LogicString(Log2It(DECinfo%FragmentExpansionRI))
   write(LU,'(a,i5)')    'Print level                                         =           ',DECinfo%PL
   write(LU,'(a,A5)')    'Fragment-adapted orbitals                           =           ',&
        & LogicString(Log2It(DECinfo%FragAdapt))
   write(LU,'(a,A5)')    'Calculate Interaction Energies                      =           ',&
        & LogicString(Log2It(DECinfo%InteractionEnergy))
   write(LU,'(a,A5)')    'Fit Molecular Orbitals                              =           ',&
        & LogicString(Log2It(DECinfo%FitOrbitals))
   !Oribtal Assignment 
   write(LU,'(A)') ' '
   if(DECinfo%BoughtonPulay) then
      write(LU,'(A)') 'DEC orbitals will be generated using Boughton-Pulay criteria'
      if(DECinfo%Distance) then
       write(LU,'(A)') 'Assignment of Molecular Orbitals to Atoms will be based Distance criteria'
      else 
         IF(DECinfo%Mulliken)then
          write(LU,'(A)') 'Assignment of Molecular Orbitals to Atoms will be based on Mulliken population analysis'
         ELSE
          write(LU,'(A)') 'Assignment of Molecular Orbitals to Atoms will be based on Lowdin charge analysis'
         ENDIF
      endif
   else 
      write(LU,'(A)') 'DEC orbitals will be generated using  simple Lowdin charge analysis'
      if(DECinfo%Distance) then
       write(LU,'(A)') 'Assignment of Molecular Orbitals to Atoms will be based Distance criteria'
      else 
       write(LU,'(A)') 'Assignment of Molecular Orbitals to Atoms will be based on Lowdin charge analysis'
      endif
   endif

   !Fragment Optimization 
   IF(DECinfo%ccModel.NE.MODEL_MP2)THEN
      write(LU,'(A)') ' '
      ! for CC models beyond MP2 (e.g. CCSD), option to use MP2 optimized fragments
      write(LU,'(A,A)')'The wave function Model used for Atomic Fragment expansion scheme = ',&
           & DECinfo%cc_models(DECinfo%fragopt_exp_model)
      IF(DECinfo%fragopt_exp_model.NE.DECinfo%ccModel)THEN
         write(LU,'(A)')'This wave function model can be changed using the .FRAGEXPMODEL keyword'
         write(LU,'(A)')'.FRAGEXPMODEL'
         write(LU,'(A1,A)')'.',DECinfo%cc_models(DECinfo%ccModel)
      ENDIF
      write(LU,'(A,A)')'Wave function Model used for Atomic Fragment reduction scheme = ',&
           & DECinfo%cc_models(DECinfo%fragopt_red_model)
      IF(DECinfo%fragopt_red_model.NE.DECinfo%ccModel)THEN
         write(LU,'(A)')'This wave function model can be changed using the .FRAGREDMODEL keyword'
         write(LU,'(A)')'.FRAGREDMODEL'
         write(LU,'(A1,A)')'.',DECinfo%cc_models(DECinfo%ccModel)
      ENDIF
   ENDIF
  ENDIF
  ! print cc parameters
  write(LU,'(/,a)') '--------------------------'
  write(LU,'(a)')   '  Coupled-cluster input   '
  write(LU,'(a,/)') '--------------------------'
  write(LU,'(a,A)')      'Wave function                =           ',DECinfo%cc_models(DECinfo%ccModel)
  write(LU,'(a,i5)')     'Maximum number of iterations =           ',DECinfo%ccMaxIter
  write(LU,'(a,g15.3)')  'Convergence threshold        = ',DECinfo%ccConvergenceThreshold
  write(LU,'(a,A5)')     'Use CROP                     =           ',LogicString(Log2It(DECinfo%use_crop))
  write(LU,'(a,i5)')     'CROP subspace                =           ',DECinfo%ccMaxDIIS
  write(LU,'(a,A5)')     'Use Preconditioner           =           ',LogicString(Log2It(DECinfo%use_preconditioner))
  write(LU,'(a,A5)')     'Precond. B                   =           ',LogicString(Log2It(DECinfo%use_preconditioner_in_b))
  write(LU,'(a,A5)')     'Debug mode                   =           ',LogicString(Log2It(DECinfo%cc_driver_debug))
  write(LU,'(a,A5)')     'CC Solver distributes Memory =           ',LogicString(Log2It(DECinfo%solver_par))
  write(LU,'(a,F16.2)')  'Memory Usage Allowed(GB)     =',DECinfo%memory
  write(LU,'(a,F16.2)')  'Backup Time Interval(s)      =',DECinfo%TimeBackup
  write(LU,'(a,A5)')     'Use F12 correction           =           ',LogicString(Log2It(DECinfo%F12))
  write(LU,'(A)')' '
  end subroutine print_dec_info

  function Log2It(LogicInput)
    implicit none
    logical,intent(in) :: LogicInput
    integer :: Log2It
    IF(LogicInput)THEN
       Log2It=1
    ELSE
       Log2It=2
    ENDIF
  end function Log2It


  !> \brief Driver for doing all fragment jobs and updating
  !> the relevant quantities (energy, density, gradient, t1).
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine fragment_jobs(nocc,nunocc,natoms,MyMolecule,mylsitem,OccOrbitals,&
       & UnoccOrbitals,jobs,AtomicFragments,FragEnergies,esti,&
       & fullgrad,t1old,t1new,fragoptjobs,estijobs,EstAtomicFragments)

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
    !> Occupied orbitals in DEC format 
    type(decorbital), intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals in DEC format 
    type(decorbital), intent(in) :: UnoccOrbitals(nunocc)
    !>  fragment job list
    type(joblist),intent(inout) :: jobs
    !> Atomic fragments with orbital spaces according to the FOT
    type(decfrag),dimension(natoms),intent(inout) :: AtomicFragments
    !> Fragment energies, see FRAGMODEL_* in dec_typedef.F90
    real(realk),intent(inout) :: FragEnergies(natoms,natoms,ndecenergies)
    !> Is this a calculation with predefined small orbital spaces used to estimate frag energies?
    logical,intent(in) :: esti
    !> MP2 gradient structure (only if DECinfo%first_order is set)
    type(fullmp2grad),intent(inout),optional :: fullgrad
    !> Old t1 amplitudes (see main_fragment_driver for details)
    type(array2),intent(in),optional :: t1old
    !> New t1 amplitudes (see main_fragment_driver for details)
    type(array2),intent(inout),optional :: t1new
    !> Job list for atomic fragment optimization  (only when esti=true)
    !> Identical to the first part of "jobs" job list but necessary for restart reasons.
    type(joblist),intent(inout),optional :: fragoptjobs
    !> Job list for estimated pair fragment calculations  (only when esti=true)
    !> Identical to the second part of "jobs" job list but necessary for restart reasons.
    type(joblist),intent(inout),optional :: estijobs
    !> Atomic fragments with estimated orbital spaces (intent(in) for practical purposes)
    type(decfrag),dimension(natoms),intent(inout),optional :: EstAtomicFragments
    type(decfrag) :: PairFragment
    integer :: k,atomA,atomB,i,j,counter,jobdone,nworkers,newjob,siz,nfragopt,estipos
    type(mp2grad) :: grad
    type(joblist) :: singlejob
    real(realk) :: fragenergy(ndecenergies)
    real(realk) :: t1cpu, t2cpu, t1wall, t2wall, dt
    integer(kind=ls_mpik) ::master, sender, groupsize,IERR
    logical :: only_update,dofragopt,backup_files
#ifdef VAR_MPI
    INTEGER(kind=ls_mpik) :: MPISTATUS(MPI_STATUS_SIZE)
#endif
    master = 0
    fragenergy=0.0_realk
    only_update=.true.

    ! Do any fragment optimizations?
    dofragopt=any(jobs%dofragopt)
    ! Number of fragment optimizations
    nfragopt = count(jobs%dofragopt)

    ! Sanity check for atomic fragment optimization 
    if(dofragopt .and. esti) then
       if(.not. present(fragoptjobs) ) then
          call lsquit('fragment_jobs: Atomic fragment optimization in job list, but &
               & fragment opt. job list is not present!',-1)
       end if
    end if

    ! Sanity check for estimated fragments
    if(any(jobs%esti)) then
       if(.not. present(EstAtomicFragments)) then
          call lsquit('fragment_jobs: Estimated pair fragments requested, but estimated &
               & atomic fragments are not present!',-1)
       end if
       if(.not. present(estijobs) ) then
          call lsquit('fragment_jobs: Estimated pair fragments requested, but estimated &
               & pair fragment job list is not present!',-1)
       end if
    end if

    ! Sanity check for singles polarization effects
    if(DECinfo%SinglesPolari) then
       if( (.not. present(t1old)) .or. (.not. present(t1new)) ) then
          call lsquit('fragment_jobs: Singles polarization requested but no &
               & t1 amplitudes present!',-1)
       end if
       if(any(jobs%esti)) then
          call lsquit('fragment_jobs: Singles polarization not implemented for estimated fragments',-1)
       end if
    end if

    ! Sanity check for first order properties
    if(DECinfo%first_order) then
       if(.not. present(fullgrad)) then
          call lsquit('fragment_jobs: First order properties requested but gradient structure &
               & is not present!',-1)
       end if
    end if


#ifdef VAR_MPI
    ! Number of workers = Number of nodes minus master itself
    nworkers = infpar%nodtot -1
    if(nworkers<1) then
       call lsquit('DEC calculations using MPI require at least two MPI processes!',-1)
    end if
#else
    nworkers=0   ! master node does all jobs
    siz=1
    call init_joblist(siz,singlejob)
#endif





    ! Define MPI groups and send fragment job list to slaves
#ifdef VAR_MPI

    ! MPI local group size: 64 or number of slaves
    ! --> this should done in a more sophisticated manner...
    if(DECinfo%MPIgroupsize>0) then
       ! group size was defined explicitly in input
       groupsize=DECinfo%MPIgroupsize
    else
       groupsize=min(64,nworkers)
    end if

    write(DECinfo%output,*) '*** MPI GROUPSIZE SET TO ', groupsize


    ! Initialize local groups 
    call init_mpi_groups(groupsize,DECinfo%output)

    ! Send fragment job list
    call bcast_dec_fragment_joblist(jobs,MPI_COMM_LSDALTON)

#endif

    ! Timing for backup
    call LSTIMER('START',t1cpu,t1wall,DECinfo%output)

    counter=0
    JobLoop: do k=1,jobs%njobs+nworkers

       ! Counter used to distinquish real (counter<=jobs%njobs) and quit (counter>jobs%njobs) 
       ! calculations
       counter=counter+1

       if(k<=jobs%njobs) then
          if(jobs%jobsdone(k)) cycle JobLoop ! job is already done
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

          ! Which atoms are involved in "jobdone"
          atomA = jobs%atom1(jobdone)
          atomB = jobs%atom2(jobdone)

          ! Receive fragment info (and update density or gradient)
          ! ------------------------------------------------------
          FragoptCheck1: if(jobs%dofragopt(jobdone)) then

             ! Sanity precaution - only fragopt for atomic fragments.
             if(atomA/=atomB) call lsquit('fragment_jobs: Fragment opt. only for atomic frags!',-1)

             ! Job is fragment optimization --> receive atomic fragment info from slave
             call mpi_send_recv_single_fragment(AtomicFragments(atomA),MPI_COMM_LSDALTON,&
                  & sender,master,singlejob)

          else

             ! No fragment optimization - receive energy (and possibly gradient) information
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

          end if FragoptCheck1


          ! Save fragment energy (symmetric for pairs)
          do j=1,ndecenergies
             if(atomA==atomB) then
                FragEnergies(atomA,atomA,j) = fragenergy(j)
             else
                FragEnergies(atomA,atomB,j) = fragenergy(j)
                FragEnergies(atomB,atomA,j) = fragenergy(j)
             end if
          end do

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

       ! Job to do: "k"
       jobdone=k


       ! Identify job
       atomA = jobs%atom1(jobdone)
       atomB = jobs%atom2(jobdone)
       ! atomA=atomB  : Single fragment job
       ! atomA/=atomB : Pair fragment job

       if(atomA==atomB) then ! single

          if(jobs%dofragopt(jobdone)) then
             write(*, '(1X,a,i6,a,i15,a,i8)') 'Job: ', jobdone, ' of size ', jobs%jobsize(jobdone),&
                &  ' is single fragment optimization: ', atomA
          else
             write(*, '(1X,a,i6,a,i15,a,i4,a,i4,a,i4,a,i6)') 'Job: ', jobdone, ' of size ', jobs%jobsize(jobdone),&
                &  ' with #o',AtomicFragments(atomA)%noccAOS,' #v', AtomicFragments(atomA)%nunoccAOS,&
                &' #b',AtomicFragments(atomA)%nbasis,' is single fragment: ', atomA
          endif

          ! Fragment "atomA" is stored in AtomicFragments(atomA).
          ! However, the basis information has not yet been initialized

          ! Get single fragment energy (and possibly density or gradient)
          ! *******************************************************************
          if(DECinfo%SinglesPolari) then
             call atomic_driver_advanced(nocc,nunocc,&
                  & OccOrbitals,UnoccOrbitals,MyLsitem,MyMolecule,&
                  & AtomicFragments(atomA),grad,t1old=t1old,t1new=t1new)
          else

             FragoptCheck2: if(jobs%dofragopt(jobdone)) then

                ! Fragment optimization
                call optimize_atomic_fragment(atomA,AtomicFragments(atomA),MyMolecule%nAtoms, &
                     & OccOrbitals,nOcc,UnoccOrbitals,nUnocc, &
                     & MyMolecule,mylsitem,.true.)

             else

                ! Init fragment basis information
                call atomic_fragment_init_basis_part(nunocc, nocc, OccOrbitals,&
                     & UnoccOrbitals,MyMolecule,mylsitem,AtomicFragments(atomA))

                ! Call main driver to get energy (and possibly density or gradient)
                call atomic_driver(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals,&
                     & AtomicFragments(atomA),grad=grad)

                ! Free basis info again
                call atomic_fragment_free_basis_info(AtomicFragments(atomA))

             end if FragoptCheck2

          end if

          ! Update energies
          do j=1,ndecenergies
             FragEnergies(atomA,atomA,j) = AtomicFragments(atomA)%energies(j)
          end do

          ! Copy fragment information into joblist
          call copy_fragment_info_job(AtomicFragments(atomA),singlejob)
          singlejob%jobsdone(1) = .true.
          singlejob%esti(1) = jobs%esti(jobdone)
          singlejob%dofragopt(1) = jobs%dofragopt(jobdone)

       else ! pair calculation


          ! Get energy (and possibly density or gradient)
          ! *********************************************

          ! Init pair
          if(jobs%esti(jobdone)) then
             ! Estimated pair fragment
             call merged_fragment_init(EstAtomicFragments(atomA), EstAtomicFragments(atomB),&
                & nunocc, nocc, natoms,OccOrbitals,UnoccOrbitals, &
                & MyMolecule,mylsitem,.true.,PairFragment,esti=.true.)

             write(*, '(1X,a,i6,a,i15,a,i4,a,i4,a,i4,a,i6,i6)') 'Job: ', jobdone, ' of size ', jobs%jobsize(jobdone),&
                &  ' with #o',PairFragment%noccAOS,' #v', PairFragment%nunoccAOS,&
                &' #b',PairFragment%nbasis,' is pair   estimate: ', atomA,atomB

          else
             ! Pair fragment according to FOT precision
             call merged_fragment_init(AtomicFragments(atomA), AtomicFragments(atomB),&
                & nunocc, nocc, natoms,OccOrbitals,UnoccOrbitals, &
                & MyMolecule,mylsitem,.true.,PairFragment)

             write(*, '(1X,a,i6,a,i15,a,i4,a,i4,a,i4,a,i6,i6)') 'Job: ', jobdone, ' of size ', jobs%jobsize(jobdone),&
                &  ' with #o',PairFragment%noccAOS,' #v', PairFragment%nunoccAOS,&
                &' #b',PairFragment%nbasis,' is pair   fragment: ', atomA,atomB

          end if


          if(DECinfo%SinglesPolari) then

             call pair_driver_singles(natoms,nocc,nunocc,&
                  & OccOrbitals,UnoccOrbitals,MyLsitem,MyMolecule,&
                  & AtomicFragments(atomA), AtomicFragments(atomB),PairFragment,t1old,t1new)

          else

             if(jobs%esti(jobdone)) then
                call pair_driver(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals,&
                     & EstAtomicFragments(atomA), EstAtomicFragments(atomB),&
                     & natoms,PairFragment,grad)
             else
                call pair_driver(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals,&
                     & AtomicFragments(atomA), AtomicFragments(atomB),&
                     & natoms,PairFragment,grad)
             end if

          end if

          ! Update pair energies
          do j=1,ndecenergies
             fragenergies(atomA,atomB,j) = PairFragment%energies(j)
             fragenergies(atomB,atomA,j) = PairFragment%energies(j)
          end do

          ! Copy fragment information into joblist
          call copy_fragment_info_job(PairFragment,singlejob)
          singlejob%jobsdone(1) = .true.
          singlejob%esti(1) = jobs%esti(jobdone)
          singlejob%dofragopt(1) = jobs%dofragopt(jobdone)

          call atomic_fragment_free(PairFragment)

       end if

       ! Update stuff (same for single and pair)
       ! ***************************************
       if(DECinfo%first_order) then  ! density and/or gradient
          call update_full_mp2gradient(grad,fullgrad)
          call free_mp2grad(grad)
       end if

#endif



       RestartStuff: if(jobdone>0) then

          ! Backup files to be able to do restart
          ! *************************************

          ! Put received job into total job list at position "jobdone"
          call put_job_into_joblist(singlejob,jobdone,jobs)

          ! For estimated fragments we also need to store in estijobs to be able to make restart
          ! files below
          if(jobs%esti(jobdone)) then
             ! Compensate for the fact that estijobs contains only estimated pair fragments
             ! and NOT atomic fragments when putting job info into joblist.
             estipos = jobdone-nfragopt
             call put_job_into_joblist(singlejob,estipos,estijobs)
          end if


          FragoptCheck3: if(jobs%dofragopt(jobdone)) then
             if(esti) then
                ! Also put job into fragopt job list to enable restart when using estimated fragments
                call put_job_into_joblist(singlejob,jobdone,fragoptjobs)
                ! Save fragment info to file atomicfragments.info
                call add_fragment_to_file(AtomicFragments(atomA),fragoptjobs)
             else
                call add_fragment_to_file(AtomicFragments(atomA),jobs)
             end if
          end if FragoptCheck3

          ! Time dt since last backup
          call LSTIMER('START',t2cpu,t2wall,DECinfo%output)
          dt = t2wall - t1wall

                    !BACKUP IF:
          !1) We are doing a fragopt, safe every converged fragment

          !2) OR 
          !   a) the finished job is one of the first quarter according to the
          !      job-list, i.e. if it is one of the largest calcluations.
          !   b) the time since the last backup is larger than DECinfo%TimeBackup

          !3) DECinfo%only_one_frag_job is requested -- this is only for debugging

          backup_files =  (((float(jobdone) < 0.25*float(jobs%njobs)) .or. &
              &(dt > DECinfo%TimeBackup) .or. all(jobs%jobsdone) ) .and. & 
              & (.not. all(jobs%dofragopt))) .or. (DECinfo%only_n_frag_jobs>0) 

          ! Backup if time passed is more than DECinfo%TimeBackup or if all jobs are done
          Backup: if( backup_files )then
             ! Note: If only fragment optimization jobs are requested this is not necessary
             ! because the fragment energies are anyway stored in add_fragment_to_file above.

             if(esti) then
                ! Save info for estimated pair fragments restart
                call write_fragment_energies_for_restart(natoms,FragEnergies,estijobs,esti)
             else
                ! Standard fragments, save info for restart
                if(DECinfo%first_order) then  ! density and/or gradient 
                   call write_gradient_and_energies_for_restart(natoms,FragEnergies,jobs,fullgrad)
                else ! just energy
                   call write_fragment_energies_for_restart(natoms,FragEnergies,jobs,esti)
                end if

             end if

             ! Reset timer
             call LSTIMER('START',t1cpu,t1wall,DECinfo%output)

          end if Backup

#ifdef VAR_MPI
          call free_joblist(singlejob)
#endif

       end if RestartStuff


    end do JobLoop

#ifdef VAR_MPI
    call MPI_COMM_FREE(infpar%lg_comm,IERR)
#else
    call free_joblist(singlejob)
#endif


  end subroutine fragment_jobs

  !> \brief Carry out fragment optimizations and (possibly) calculate 
  !> estimated pair fragment energies.
  !> \author Kasper Kristensen
  !> \date November 2013
  subroutine fragopt_and_estimated_frags(nOcc,nUnocc,OccOrbitals,UnoccOrbitals, &
       & MyMolecule,mylsitem,dofrag,esti,AtomicFragments,FragEnergies)

    implicit none
    !> Full molecule info (CC model for each pair fragment resulting from estimate analysis is stored 
    !> in MyMolecule%ccmodel)
    type(fullmolecule), intent(inout) :: MyMolecule
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> Number of occupied orbitals in full molecule
    integer, intent(in) :: nocc
    !> Number of unoccupied orbitals in full molecule
    integer, intent(in) :: nunocc
    !> Information about DEC occupied orbitals
    type(decorbital), dimension(nOcc), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(decorbital), dimension(nUnocc), intent(in) :: UnoccOrbitals
    !> List of which atoms have orbitals assigned
    logical,intent(in),dimension(MyMolecule%natoms) :: dofrag
    !> Is this a calculation with predefined small orbital spaces used to estimate frag energies?
    !> (This is effectively intent(in), only intent(inout) in practice for MPI purposes).
    logical,intent(inout) :: esti
    !> Atomic Fragments with orbital spaces determined to the FOT precision
    type(decfrag), intent(inout),dimension(MyMolecule%natoms) :: AtomicFragments
    !> Fragment energies 
    real(realk),intent(inout) :: FragEnergies(MyMolecule%natoms,MyMolecule%natoms,ndecenergies)
    real(realk),pointer :: FragEnergiesPart(:,:)
    type(decfrag),pointer :: EstAtomicFragments(:)
    logical :: DoBasis,calcAF
    real(realk) :: init_radius,tcpu1,twall1,tcpu2,twall2,mastertime,Epair_est,Eskip_est
    integer :: natoms,i,j,k
    type(joblist) :: jobs,fragoptjobs,estijobs
    integer(kind=ls_mpik) :: master
    master=0

    call LSTIMER('START',tcpu1,twall1,DECinfo%output)


    natoms= MyMolecule%natoms

    ! Initialize job list for atomic fragment optimizations
    call create_dec_joblist_fragopt(natoms,nocc,nunocc,MyMolecule%DistanceTable,&
       & OccOrbitals, UnoccOrbitals, dofrag, mylsitem,fragoptjobs)

    if(DECinfo%DECrestart) then
       write(DECinfo%output,*) 'Restarting atomic fragment optimizations....'
       call restart_atomic_fragments_from_file(natoms,MyMolecule,MyLsitem,OccOrbitals,&
          & UnoccOrbitals,.false.,AtomicFragments,fragoptjobs)
    end if

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '*** Initiating optimization of atomic fragments ***'
    write(DECinfo%output,*)
    write(DECinfo%output,*)


    if(esti) then

       ! fragment opt. AND estimate pair fragments
       ! *****************************************

       ! All estimates are done using MP2, set model for all fragments to MP2
       ! (MyMolecule%ccmodel is then redefined below based on these estimates)
       do i=1,natoms
          do j=i+1,natoms
             MyMolecule%ccmodel(i,j) = MODEL_MP2
             MyMolecule%ccmodel(j,i) = MODEL_MP2
          end do
       end do
       IF(DECinfo%InteractionEnergy)THEN
          !We are only intrested in certain fragments. 
          !So we only do estimates on those we are 
          !intrested in. 
          do i=1,natoms
             do j=1,natoms
                if(MyMolecule%SubSystemIndex(i).EQ.MyMolecule%SubSystemIndex(j))then
                   MyMolecule%ccmodel(i,j) = MODEL_NONE
                endif
             end do
          end do
       endif

!
       ! Init estimated atomic fragments by including orbitals assigned to neighbour atoms
       ! within 2 Angstrom.
       init_radius = DECinfo%EstimateINITradius
       DoBasis = .false.
       call mem_alloc(EstAtomicFragments,natoms)
       call init_estimated_atomic_fragments(nOcc,nUnocc,OccOrbitals,UnoccOrbitals, &
            & MyMolecule,mylsitem,DoBasis,init_radius,dofrag,EstAtomicFragments)

#ifdef VAR_MPI
       ! Send estimated fragment information to slaves
       call mpi_bcast_many_fragments(natoms,dofrag,EstAtomicFragments,MPI_COMM_LSDALTON)

       ! Send CC models to use for each pair based on estimates (always MP2, but keep it general)
       call ls_mpibcast(MyMolecule%ccmodel,natoms,natoms,master,MPI_COMM_LSDALTON)
#endif

       ! Get job list for estimated pair fragments
       calcAF = .false.  ! No atomic fragments, just pairs
       call create_dec_joblist_driver(calcAF,MyMolecule,mylsitem,natoms,nocc,nunocc,&
            &OccOrbitals,UnoccOrbitals,EstAtomicFragments,dofrag,esti,estijobs)
       if(DECinfo%DECrestart) then
          write(DECinfo%output,*) 'Restarting pair fragment estimate calculations...'
          call read_fragment_energies_for_restart(natoms,FragEnergies,estijobs,esti)
       end if

       ! Merge job list of atomic fragment optimization and estimated fragments (in this order)
       call concatenate_joblists(fragoptjobs,estijobs,jobs)

       call fragment_jobs(nocc,nunocc,natoms,MyMolecule,mylsitem,OccOrbitals,&
            & UnoccOrbitals,jobs,AtomicFragments,FragEnergies,esti,&
            & fragoptjobs=fragoptjobs,estijobs=estijobs,EstAtomicFragments=EstAtomicFragments)

       do i=1,natoms
          if(dofrag(i)) then
             call atomic_fragment_free_simple(EstAtomicFragments(i))
          end if
       end do
       call mem_dealloc(EstAtomicFragments)
    else

       ! Just fragment opt.
       ! ******************
       call fragment_jobs(nocc,nunocc,natoms,MyMolecule,mylsitem,OccOrbitals,&
            & UnoccOrbitals,fragoptjobs,AtomicFragments,FragEnergies,esti)

    end if



    if(esti) then
       ! Get estimated pair fragment energies for occupied partitioning scheme
       call mem_alloc(FragEnergiesPart,natoms,natoms)
       IF(DECinfo%onlyVirtPart)THEN
          call get_virtfragenergies(natoms,MODEL_MP2,FragEnergies,FragEnergiesPart)
       ELSE
          call get_occfragenergies(natoms,MODEL_MP2,FragEnergies,FragEnergiesPart)
       ENDIF
       ! We do not want to consider atomic fragment energies now so zero them
       ! (they might be zero already but in this way we avoid wrong print out below).
       do i=1,natoms
          FragEnergiesPart(i,i)=0.0_realk
       end do
       ! Define which model to use in each pair calculation (info stored in MyMolecule%ccmodel)
       ! (Also calculate estimated correlation energy and estimated error by skipping pairs).
       call define_pair_calculations(natoms,dofrag,FragEnergiesPart,MyMolecule,Epair_est,Eskip_est)
       call mem_dealloc(FragEnergiesPart)
    end if


    ! Zero all pair fragment energies to make sure that there are no leftovers from estimated pairs.
    do k=1,ndecenergies
       do j=1,natoms
          do i=j+1,natoms
             FragEnergies(i,j,k) = 0.0_realk
             FragEnergies(j,i,k) = 0.0_realk
          end do
       end do
    end do


    ! Timings and statistics print
    call LSTIMER('START',tcpu2,twall2,DECinfo%output)
    mastertime = twall2-twall1
    if(esti) then
#ifdef VAR_MPI
       call print_MPI_fragment_statistics(jobs,mastertime,'FRAGMENT OPT. + ESTIM.')
#endif
    else
#ifdef VAR_MPI
       call print_MPI_fragment_statistics(fragoptjobs,mastertime,'FRAGMENT OPT.')
#endif
    end if

    ! Free job list(s)
    call free_joblist(fragoptjobs)
    if(esti) then
       call free_joblist(jobs)
       call free_joblist(estijobs)
    end if



  end subroutine fragopt_and_estimated_frags


end module dec_driver_module


