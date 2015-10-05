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
  use crayio_tools_module
  use tensor_type_def_module
  use tensor_interface_module
  use background_buffer_module
  use lsparameters

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
  use ccsd_gradient_module

#ifdef VAR_MPI
  use decmpi_module
  use infpar_module
  use dec_driver_slave_module
  use lsmpi_module
  use lsmpi_param
  use lsmpi_type
#endif


public:: DEC_wrapper,main_fragment_driver
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
    type(decorbital), pointer :: virtOrbitals(:)
    real(realk),pointer :: FragEnergiesOcc(:,:)
    integer :: nBasis,nOcc,nvirt,nAtoms,i
    type(decfrag),pointer :: AtomicFragments(:)


    ! Print DEC info
    call print_dec_info()

    nOcc   = MyMolecule%nocc
    nvirt = MyMolecule%nvirt
    nBasis = MyMolecule%nbasis
    nAtoms = MyMolecule%natoms

    ! -- Analyze basis and create orbitals
    call mem_alloc(OccOrbitals,nOcc)
    call mem_alloc(virtOrbitals,nvirt)
    call GenerateOrbitals_driver(MyMolecule,mylsitem,nocc,nvirt,natoms, &
         & OccOrbitals, virtOrbitals)

    if (DECinfo%only_generate_DECorbs) then
       print*, 'Quits after DEC orbitals have been written to file'
       write(DECinfo%output,*)
       write(DECinfo%output,'(1X,a)') 'Quits after DEC orbitals have been written to file'
       call lsquit('Quitting after DEC orbitals have been written to file',DECinfo%output)
    end if



    ! *************************************************
    ! Optimize all atomic fragments and calculate pairs
    ! *************************************************
    call mem_alloc(FragEnergiesOcc,MyMolecule%nfrags,MyMolecule%nfrags)
    call mem_alloc(AtomicFragments,MyMolecule%nfrags)
    call main_fragment_driver(MyMolecule,mylsitem,D,&
         &OccOrbitals,virtOrbitals,natoms,nocc,nvirt,EHF,Ecorr,molgrad,&
         & Eerr,FragEnergiesOcc,AtomicFragments,.false.)


    ! Delete orbitals
    do i=1,nOcc
       call orbital_free(OccOrbitals(i))
    end do

    do i=1,nvirt
       call orbital_free(virtOrbitals(i))
    end do
    if(DECinfo%PL>0) write(DECinfo%output,'(a)') 'Orbitals removed'

    call mem_dealloc(OccOrbitals)
    call mem_dealloc(virtOrbitals)
    call mem_dealloc(FragEnergiesOcc)

    do i=1,MyMolecule%nfrags
       if(.not. associated(AtomicFragments(i)%EOSatoms)) cycle
       call atomic_fragment_free_simple(AtomicFragments(i))
    end do
    call mem_dealloc(AtomicFragments)

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
       & OccOrbitals,virtOrbitals, &
       & natoms,nocc,nvirt,EHF,Ecorr,molgrad,dE_est1,FragEnergiesOcc,AtomicFragments,AFset)

    implicit none
    !> Number of occupied orbitals in full molecule (not changed, inout for MPI reasons)
    integer,intent(inout) :: nOcc
    !> Number of virtupied orbitals in full molecule (not changed, inout for MPI reasons)
    integer,intent(inout) :: nvirt
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
    !> virtupied orbitals in DEC format (not changed at output, is intent(inout) for MPI purposes)
    type(decorbital), intent(inout) :: virtOrbitals(nvirt)
    !> Hartree-Fock energy
    real(realk),intent(inout) :: EHF
    !> MP2 correlation energy
    real(realk),intent(inout) :: Ecorr
    !> Molecular gradient (only calculated if DECinfo%gradient is set!)
    real(realk),intent(inout) :: molgrad(3,natoms)
    !> Estimated energy error
    real(realk),intent(inout) :: dE_est1
    !> Fragment energies for occupied partitioning scheme
    real(realk),dimension(MyMolecule%nfrags,MyMolecule%nfrags),intent(inout) :: FragEnergiesOcc
    logical :: esti
    ! Fragment energies
    !real(realk) :: FragEnergies(natoms,natoms,ndecenergies)
    real(realk),pointer :: FragEnergies(:,:,:) !(natoms,natoms,ndecenergies)
    !> Atomic fragments to be optimized (AFset=false) or already set by input (AFset=true)
    type(decfrag),intent(inout) :: AtomicFragments(MyMolecule%nfrags)
    !> Have atomic fragments already been set by input (true) 
    !> or do they need to be optimized (false)
    !> Note: If AFset=true then atomic fragment energies will always be calculated here!
    logical,intent(in) :: AFset
    integer :: i,j,k,dims(2),nbasis,counter
    real(realk) :: Esos,Eerrs
    !real(realk) :: energies(ndecenergies)
    real(realk),pointer :: energies(:)!(ndecenergies)
    logical :: dens_save ! Internal control of MP2 density keyword
    logical :: FO_save  ! Internal control of first order property keyword
    logical :: grad_save  ! Internal control of MP2 gradient keyword
    logical :: dofrag(MyMolecule%nfrags)
    type(array2) :: t1old,fockt1,t1new
    logical :: redo
    type(mp2grad) :: grad
    type(fullmp2grad) :: fullgrad
    integer :: jobdone,newjob, nworkers, jobidx
    integer(kind=ls_mpik) :: groupsize
    integer :: MPIdatatype,MyAtom,nfrags
    type(joblist) :: jobs
    real(realk) :: tcpu,twall,oldpaircut,newpaircut,tcpu1,tcpu2,twall1,twall2,mastertime
    ! (T) contribution to fragment energies for occupied (:,:,1), and virtual (:,:,2) schemes 
    !> (:,:,3): Occupied E[4] contribution;  (:,:,4): Virtual E[4] contribution
    !> (:,:,5): Occupied E[5] contribution;  (:,:,6): Virtual E[5] contribution
    logical :: calcAF,ForcePrintTime
    integer(kind=ls_mpik) :: master,IERR,comm,sender
    real(realk) :: Edft,dE_est2,dE_est3,Eskip_est
#ifdef VAR_MPI
    INTEGER(kind=ls_mpik) :: MPISTATUS(MPI_STATUS_SIZE), DUMMYSTAT(MPI_STATUS_SIZE)
#endif
    Eskip_est = 0.0E0_realk
    dE_est3 = 0.0E0_realk
    master=0
    ForcePrintTime = .TRUE.

    ! Number of potential fragments
    nfrags = MyMolecule%nfrags


    call LSTIMER('START',tcpu,twall,DECinfo%output)



    ! ************************************************************************
    ! *                           Initialize stuff                           *
    ! ************************************************************************

    redo=.false.
    nbasis = MyMolecule%nbasis
    if(.not. AFset) then
       do i=1,nfrags
          call atomic_fragment_nullify(AtomicFragments(i))
       end do
    end if

    call mem_alloc(FragEnergies,nfrags,nfrags,ndecenergies)
    FragEnergies=0E0_realk

    ! Which fragments should we consider
    call which_fragments_to_consider(MyMolecule%ncore,nocc,nvirt,nfrags,&
         & OccOrbitals,virtOrbitals,dofrag,MyMolecule%PhantomAtom)

    if(DECinfo%StressTest)then
       call StressTest_mod_dofrag(MyMolecule%natoms,nocc,nvirt,MyMolecule%ncore,&
            & MyMolecule%DistanceTable,OccOrbitals, virtOrbitals, dofrag, mylsitem)
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
       dims(1) = nvirt
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
    call mpi_dec_fullinfo_master_to_slaves_precursor(esti,nocc,nvirt,master)

    ! Pass remaining full molecular info to local masters
    call mpi_dec_fullinfo_master_to_slaves(nocc,nvirt,&
         & OccOrbitals, virtOrbitals, MyMolecule, MyLsitem)

#else
    nworkers=0   ! master node does all jobs
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
    call fragopt_and_estimated_frags(nOcc,nvirt,OccOrbitals,virtOrbitals, &
         & MyMolecule,mylsitem,dofrag,esti,AtomicFragments,FragEnergies, &
         & AFset,dE_est3,Eskip_est)

    ! Send CC models and pair FOTs to use for all pairs based on estimates
    if(esti) then
#ifdef VAR_MPI
       call ls_mpibcast(MyMolecule%ccmodel,nfrags,nfrags,master,MPI_COMM_LSDALTON)
       call ls_mpibcast(MyMolecule%pairfotlevel,nfrags,nfrags,master,MPI_COMM_LSDALTON)
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
    IsAfSet: if(.not. AFset) then
       do i=1,nfrags
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
    end if IsAfSet

    ! Now all atomic fragment energies have been calculated and the
    ! fragment information has been stored in AtomicFragments.

    ! ************************************************************************
    ! *             Construct job list for remaining fragments               *
    ! ************************************************************************

    ! Restore first order
    DECinfo%density=dens_save
    DECinfo%first_order = FO_save
    DECinfo%gradient = grad_save

    ! Get job list 
    calcAF = DECinfo%RepeatAF
    !This is a hack to specify that only pair fragment jobs should be done
    if(DECinfo%only_pair_frag_jobs) calcAF = .false.
    ! Always calculate atomic fragment energies if AFset is true
    if(AFset) calcAF=.true.

    call create_dec_joblist_driver(calcAF,MyMolecule,mylsitem,nfrags,nocc,nvirt,&
         &OccOrbitals,virtOrbitals,AtomicFragments,dofrag,.false.,jobs)

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
    call mpi_bcast_many_fragments(nfrags,dofrag,AtomicFragments,MPI_COMM_LSDALTON)

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
          if(MyMolecule%mem_distributed)then
             call lsquit("ERROR(main_fragment_driver): singles polarization does &
             &not work with distributed molecular structure",-1)
          endif

          fockt1 = array2_init([nbasis,nbasis])
          ! Get T1 transformed Fock matrix in AO basis
          call fullmolecular_get_AOt1Fock(mylsitem,MyMolecule,t1old,fockt1)

          ! Set fock matrix associated with fullmolecule structure
          ! equal to T1 transformed Fock matrix
          MyMolecule%fock%elm2(1:nbasis,1:nbasis) = fockt1%val(1:nbasis,1:nbasis)
          call array2_free(fockt1)

          ! Init improved full molecular singles constructed from
          ! single and pair fragment calculations.
          dims(1) = nvirt
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
             call read_gradient_and_energies_for_restart(nfrags,FragEnergies,jobs,fullgrad)
          else
             write(DECinfo%output,*) 'Restarting pair fragments - energy...'
             call read_fragment_energies_for_restart_driver(nfrags,FragEnergies,jobs,esti)
          end if
       end if
       if(DECinfo%only_n_frag_jobs > 0)then
          jobs%jobsdone  = .false.
          jobs%dofragopt = .false.
       endif

       call fragment_jobs(nocc,nvirt,nfrags,MyMolecule,mylsitem,&
            & OccOrbitals,virtOrbitals,jobs,AtomicFragments,&
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
    if (.not. DECinfo%no_pairs) then 
       IF(DECinfo%onlyVirtPart)THEN
          call get_virtfragenergies(nfrags,DECinfo%ccmodel,FragEnergies,FragEnergiesOcc)
       ELSE
          call get_occfragenergies(nfrags,DECinfo%ccmodel,FragEnergies,FragEnergiesOcc)
       ENDIF
       call plot_pair_energies(nfrags,DECinfo%pair_distance_threshold,FragEnergiesOcc,&
            & MyMolecule,dofrag)
    end if

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
    else if (.not.DECinfo%no_pairs) then
       call print_MPI_fragment_statistics(jobs,mastertime,'PAIR FRAGMENTS')
    end if
#endif

    call mem_alloc(energies,ndecenergies)
    ! Total correlation energy 
    do j=1,ndecenergies
       call add_dec_energies(nfrags,FragEnergies(:,:,j),dofrag,energies(j))
    end do

    ! Print all fragment energies
    call print_all_fragment_energies(nfrags,FragEnergies,dofrag,mymolecule%DistanceTable,energies)

    call mem_dealloc(FragEnergies)
    !Obtain The Correlation Energy from the list energies
    call ObtainModelEnergyFromEnergies(DECinfo%ccmodel,energies,Ecorr)
    if(DECinfo%ccmodel == MODEL_RPA) then
       call ObtainModelEnergyFromEnergies(DECinfo%ccmodel,energies,Esos,.true.)
    endif
    ! If singles polarization was considered, we need to
    ! ensure that the fullmolecule structure contains the standard
    ! (NOT T1 transformed) Fock matrix at output
    if(DECinfo%SinglesPolari) then
       call tensor_free(MyMolecule%fock)
       call molecule_get_fock(MyMolecule,mylsitem)  ! Get standard Fock matrix again
       call array2_free(t1old)
       call array2_free(t1new)
    end if

    ! HF energy
    Ehf = get_HF_energy_fullmolecule(MyMolecule,Mylsitem,D) 
    Edft =0.0_realk
    !DFT energy
    if(DECinfo%DFTreference) then
       Edft = get_dft_energy_fullmolecule(MyMolecule,Mylsitem,D) 
    endif

    ! If requested, calculate MP2 gradient and MP2 density (or just density) for
    ! full molecule from fragment contributions
    if(DECinfo%SkipFull) then
       write(DECinfo%output,*) 'WARNING: SKIPPING FULL MOLECULAR PART &
            & OF DENSITY/GRADIENT CALCULATION!'
    else
       if(DECinfo%first_order) then
         if(DECinfo%ccmodel == MODEL_CCSD)then
           fullgrad%EHF = EHF
           ! Calculate CCSD density/gradient and save density matrices to file
           call get_CCSDgradient_main(MyMolecule,mylsitem,D,molgrad,fullgrad)
         end if
         if(DECinfo%ccmodel == MODEL_MP2 .or. DECinfo%ccmodel == MODEL_RIMP2) then
           fullgrad%EHF = EHF
           ! Calculate MP2 density/gradient and save MP2 density matrices to file
           call get_mp2gradient_main(MyMolecule,mylsitem,D,molgrad,fullgrad)
         end if

       end if
    end if

    if(DECinfo%first_order) then
       call free_fullmp2grad(fullgrad)
    end if

    call free_joblist(jobs)

    ! Estimate energy error
    call get_estimated_energy_error(nfrags,energies,dE_est1,dE_est2)
    if(DECinfo%ccmodel == MODEL_RPA) then
       call get_estimated_energy_error(nfrags,energies,Eerrs,dE_est2,.true.)
    endif
    call mem_dealloc(energies)

    if(DECinfo%NaturalLinearScalingF12Terms) then
       print*,'NaturalLinearScalingTerms:'
       print*,'add NaturalLinearScalingTerms: Ecorr               =',Ecorr
       print*,'add NaturalLinearScalingTerms: MyMolecule%EF12NLSV1=',MyMolecule%EF12NLSV1
       print*,'add NaturalLinearScalingTerms: MyMolecule%EF12NLSB1=',MyMolecule%EF12NLSB1
       print*,'add NaturalLinearScalingTerms: MyMolecule%EF12NLSX1=',MyMolecule%EF12NLSX1
       Ecorr = Ecorr + MyMolecule%EF12NLSB1 + MyMolecule%EF12NLSV1 + MyMolecule%EF12NLSX1
       print*,'========================================================================='
       print*,'add NaturalLinearScalingTerms: New Ecorr         =',Ecorr
    endif

    ! Print short summary
    call print_total_energy_summary(EHF,Edft,Ecorr,MyMolecule%EF12singles,&
         & dE_est1,dE_est2,dE_est3,Eskip_est)
    if(DECinfo%ccmodel==MODEL_RPA)then
       call print_total_energy_summary(EHF,Edft,Esos,MyMolecule%EF12singles,&
            &  Eerrs,dE_est2,dE_est3,Eskip_est,doSOS=.true.)
    endif
    call LSTIMER('DEC FINAL',tcpu,twall,DECinfo%output,ForcePrintTime)

  end subroutine main_fragment_driver

  subroutine ObtainModelEnergyFromEnergies(ccmodel,energies,Ecorr,doSOS)
    implicit none
    integer,intent(in)        ::  ccmodel
    real(realk),intent(in)    ::  energies(ndecenergies)
    real(realk),intent(inout) ::  Ecorr
    logical,intent(in),optional :: doSOS
    logical ::                     SOS
    ! MODIFY FOR NEW MODEL
    ! MODIFY FOR NEW CORRECTION: Add correction to output energy
    ! Set output energy: We choose occupied partitioning scheme energy as general output
    SOS = .false.
    if(present(doSOS)) SOS = doSOS
    select case(DECinfo%ccmodel)
    case(MODEL_MP2)
       ! MP2, use occ energy
       IF(DECinfo%onlyVirtPart)THEN
          Ecorr = energies(FRAGMODEL_VIRTMP2)
       ELSE
          Ecorr = energies(FRAGMODEL_OCCMP2)
       ENDIF
       if(DECinfo%F12) then
          Ecorr = Ecorr + energies(FRAGMODEL_MP2f12)
       endif
    case(MODEL_RIMP2)
       ! RI-MP2, use occ energy
       IF(DECinfo%onlyVirtPart)THEN
          Ecorr = energies(FRAGMODEL_VIRTRIMP2)
       ELSE
          Ecorr = energies(FRAGMODEL_OCCRIMP2)
       ENDIF
       if(DECinfo%F12) then
          Ecorr = Ecorr + energies(FRAGMODEL_RIMP2f12)
       endif
    case(MODEL_LSTHCRIMP2)
       ! LS-THC-RI-MP2, use occ energy
       IF(DECinfo%onlyVirtPart)THEN
          Ecorr = energies(FRAGMODEL_VIRTLSTHCRIMP2)
       ELSE
          Ecorr = energies(FRAGMODEL_OCCLSTHCRIMP2)
       ENDIF
    case(MODEL_RPA)
       ! RPA, use occ energy
       if(SOS) then
         IF(DECinfo%onlyVirtPart)THEN
           Ecorr = energies(FRAGMODEL_VIRTSOS)
         ELSE
           Ecorr = energies(FRAGMODEL_OCCSOS)
         ENDIF
       ELSE
         IF(DECinfo%onlyVirtPart)THEN
           Ecorr = energies(FRAGMODEL_VIRTRPA)
         ELSE
           Ecorr = energies(FRAGMODEL_OCCRPA)
         ENDIF
       ENDIF
    case(MODEL_SOSEX)
      IF(DECinfo%onlyVirtPart)THEN
        Ecorr = energies(FRAGMODEL_VIRTSOS)
      ELSE
        Ecorr = energies(FRAGMODEL_OCCSOS)
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
       if(DECinfo%F12) then
          Ecorr = Ecorr + energies(FRAGMODEL_CCSDf12)
       endif
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
        if (DECinfo%decco) then 
           write(LU,'(a,i5)') 'Use Pair Estimate initialisation number of atom     = ',&
             & DECinfo%estimateInitAtom
        else
           write(LU,'(a,e15.3)') 'Use Pair Estimate initialisation radius (Angstrom)  = ',&
             & DECinfo%estimateINITradius*bohr_to_angstrom
        end if
        ! Pair estimates
        write(LU,'(a,a)')     'Pair estimate model                                 =     ', &
             & DECinfo%cc_models(DECinfo%PairEstimateModel)
   endif

   write(LU,'(a,g15.3)') 'Simple orbital thr.                                 = ',&
        & DECinfo%simple_orbital_threshold
   write(LU,'(a,i5)')    'Expansion step size                                 =           ',&
        & DECinfo%Frag_Exp_Size
   write(LU,'(a,i5)')    'Print level                                         =           ',DECinfo%PL
   write(LU,'(a,A5)')    'Fragment-adapted orbitals                           =           ',&
        & LogicString(Log2It(DECinfo%FragAdapt))
   write(LU,'(a,ES12.4)')'The Integral Screening threshold                    =           ',&
        & DECinfo%IntegralThreshold

   !Oribtal Assignment 
   write(LU,'(A)') ' '
   write(LU,'(A)') 'DEC orbitals will be generated using  simple Lowdin charge analysis'
   if(DECinfo%Distance) then
      write(LU,'(A)') 'Assignment of Molecular Orbitals to Atoms will be based Distance criteria'
   else 
      write(LU,'(A)') 'Assignment of Molecular Orbitals to Atoms will be based on Lowdin charge analysis'
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
  write(LU,*)

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
  subroutine fragment_jobs(nocc,nvirt,nfrags,MyMolecule,mylsitem,OccOrbitals,&
       & virtOrbitals,jobs,AtomicFragments,FragEnergies,esti,&
       & fullgrad,t1old,t1new,fragoptjobs,estijobs,EstAtomicFragments)

    implicit none
    !> Number of occupied orbitals in full molecule
    integer,intent(in) :: nOcc
    !> Number of virtupied orbitals in full molecule
    integer,intent(in) :: nvirt
    !> Number of fragments (# atoms for atom-based DEC, #occ orbitals for DECCO)
    integer,intent(in) :: nfrags
    !> Full molecule info (will be unchanged at output)
    type(fullmolecule), intent(inout) :: MyMolecule
    !> LSDalton info
    type(lsitem), intent(inout) :: mylsitem
    !> Occupied orbitals in DEC format 
    type(decorbital), intent(in) :: OccOrbitals(nocc)
    !> virtupied orbitals in DEC format 
    type(decorbital), intent(in) :: virtOrbitals(nvirt)
    !>  fragment job list
    type(joblist),intent(inout) :: jobs
    !> Atomic fragments with orbital spaces according to the FOT
    type(decfrag),dimension(nfrags),intent(inout) :: AtomicFragments
    !> Fragment energies, see FRAGMODEL_* in dec_typedef.F90
    real(realk),intent(inout) :: FragEnergies(nfrags,nfrags,ndecenergies)
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
    type(decfrag),dimension(nfrags),intent(inout),optional :: EstAtomicFragments
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
    dofragopt=.false.
    nfragopt = 0

    if (jobs%njobs>0) then
       ! Do any fragment optimizations?
       dofragopt=any(jobs%dofragopt)
       ! Number of fragment optimizations
       nfragopt = count(jobs%dofragopt)

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
    end if 

    ! Sanity check for atomic fragment optimization 
    if(dofragopt .and. esti) then
       if(.not. present(fragoptjobs) ) then
          call lsquit('fragment_jobs: Atomic fragment optimization in job list, but &
               & fragment opt. job list is not present!',-1)
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
                &  ' with #o',AtomicFragments(atomA)%noccAOS,' #v', AtomicFragments(atomA)%nvirtAOS,&
                &' #b',AtomicFragments(atomA)%nbasis,' is single fragment: ', atomA
          endif

          ! Fragment "atomA" is stored in AtomicFragments(atomA).
          ! However, the basis information has not yet been initialized

          ! Get single fragment energy (and possibly density or gradient)
          ! *******************************************************************
          if(DECinfo%SinglesPolari) then
             call atomic_driver_advanced(nocc,nvirt,&
                  & OccOrbitals,virtOrbitals,MyLsitem,MyMolecule,&
                  & AtomicFragments(atomA),grad,t1old=t1old,t1new=t1new)
          else

             FragoptCheck2: if(jobs%dofragopt(jobdone)) then

                ! Fragment optimization
                call optimize_atomic_fragment(atomA,AtomicFragments(atomA),MyMolecule%nfrags, &
                     & OccOrbitals,nOcc,virtOrbitals,nvirt, &
                     & MyMolecule,mylsitem,.true.)

             else

                ! Init fragment basis information
                call atomic_fragment_init_basis_part(nvirt, nocc, OccOrbitals,&
                     & virtOrbitals,MyMolecule,mylsitem,AtomicFragments(atomA))

                ! Call main driver to get energy (and possibly density or gradient)
                call atomic_driver(MyMolecule,mylsitem,OccOrbitals,virtOrbitals,&
                     & AtomicFragments(atomA),grad=grad)

                ! Free basis info again
                call atomic_fragment_free_basis_info(AtomicFragments(atomA))

             end if FragoptCheck2

          end if

          ! Update energies
          do j=1,ndecenergies
             FragEnergies(atomA,atomA,j) = AtomicFragments(atomA)%energies(j)
          end do

       if(DECinfo%first_order) then  ! density and/or gradient
        !Update the CCSD density
        if(DECinfo%ccmodel == MODEL_CCSD)then
          call update_CCSDdensity(AtomicFragments(atomA),grad,fullgrad)
        endif
       endif

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
                & nvirt, nocc, nfrags,OccOrbitals,virtOrbitals, &
                & MyMolecule,mylsitem,.true.,PairFragment,esti=.true.)

             write(*, '(1X,a,i6,a,i15,a,i4,a,i4,a,i4,a,i6,i6)') 'Job: ', jobdone, ' of size ', jobs%jobsize(jobdone),&
                &  ' with #o',PairFragment%noccAOS,' #v', PairFragment%nvirtAOS,&
                &' #b',PairFragment%nbasis,' is pair   estimate: ', atomA,atomB

          else
             ! Pair fragment according to FOT precision
             call merged_fragment_init(AtomicFragments(atomA), AtomicFragments(atomB),&
                & nvirt, nocc, nfrags,OccOrbitals,virtOrbitals, &
                & MyMolecule,mylsitem,.true.,PairFragment)

             write(*, '(1X,a,i6,a,i15,a,i4,a,i4,a,i4,a,i6,i6)') 'Job: ', jobdone, ' of size ', jobs%jobsize(jobdone),&
                &  ' with #o',PairFragment%noccAOS,' #v', PairFragment%nvirtAOS,&
                &' #b',PairFragment%nbasis,' is pair   fragment: ', atomA,atomB

          end if


          if(DECinfo%SinglesPolari) then

             call pair_driver_singles(nfrags,nocc,nvirt,&
                  & OccOrbitals,virtOrbitals,MyLsitem,MyMolecule,&
                  & AtomicFragments(atomA), AtomicFragments(atomB),PairFragment,t1old,t1new)

          else

             if(jobs%esti(jobdone)) then
                call pair_driver(EstAtomicFragments(atomA), EstAtomicFragments(atomB),&
                     & PairFragment,grad)
             else
                call pair_driver(AtomicFragments(atomA), AtomicFragments(atomB),&
                     & PairFragment,grad)
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

       if(DECinfo%first_order) then  ! density and/or gradient
        !Update the CCSD density
        if(DECinfo%ccmodel == MODEL_CCSD)then
          call update_CCSDdensity(PairFragment,grad,fullgrad)
        endif
       endif

          call atomic_fragment_free(PairFragment)

       end if

       ! Update stuff (same for single and pair)
       ! ***************************************
       if(DECinfo%first_order) then  ! density and/or gradient
!        if(DECinfo%ccmodel == MODEL_CCSD)then
!          call update_full_CCSDgradient(grad,fullgrad)
!          call free_mp2grad(grad)
!        endif
        if(DECinfo%ccmodel == MODEL_MP2 .or. DECinfo%ccmodel == MODEL_RIMP2) then
          call update_full_mp2gradient(grad,fullgrad)
          call free_mp2grad(grad)
        end if
       end if

#endif


       ! Crash calculation on purpose to test restart option
       if (( jobdone >= 2*nfragopt-1).and.(DECinfo%CrashEsti)) then
          print*,'Calculation was intentionally crashed due to keyword .CRASHESTI'
          print*,'This keyword is only used for debug and testing purposes'
          print*,'We want to be able to test the .RESTART keyword'
          WRITE(DECinfo%output,*)'Calculation was intentionally crashed due to keyword .CRASHESTI'
          WRITE(DECinfo%output,*)'This keyword is only used for debug and testing purposes'
          WRITE(DECinfo%output,*)'We want to be able to test the .RESTART keyword'
          call lsquit('Crashed Calculation due to .CRASHESTI keyword',DECinfo%output)
       end if


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

          ! BACKUP IF:
          ! ----------
          !
          ! 1) It is time to back up files base on DECinfo%TimeBacku
          ! 2) All jobs are done
          ! 3) DECinfo%CRASHESTI is true (Debug and testing of estimate restart files)
          ! 4) (DECinfo%only_n_frag_jobs>0) Mainly for debugging 
          !
          backup_files =  (dt > DECinfo%TimeBackup) .or. all(jobs%jobsdone) & 
             & .or. DECinfo%CRASHESTI .or. (DECinfo%only_n_frag_jobs>0)

          Backup: if( backup_files )then
             ! Note: If only fragment optimization jobs are requested this is not necessary
             ! because the fragment energies are anyway stored in add_fragment_to_file above.

             if(esti) then
                ! Save info for estimated pair fragments restart
                call write_fragment_energies_for_restart_driver(nfrags,FragEnergies,estijobs,esti)
             else
                ! Standard fragments, save info for restart
                if(DECinfo%first_order) then  ! density and/or gradient 
                   call write_gradient_and_energies_for_restart(nfrags,FragEnergies,jobs,fullgrad)
                else ! just energy
                   call write_fragment_energies_for_restart_driver(nfrags,FragEnergies,jobs,esti)
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
    infpar%lg_comm  = MPI_COMM_LSDALTON
    infpar%lg_mynum = infpar%mynum
#else
    call free_joblist(singlejob)
#endif


  end subroutine fragment_jobs

  !> \brief Carry out fragment optimizations and (possibly) calculate 
  !> estimated pair fragment energies.
  !> \author Kasper Kristensen
  !> \date November 2013
  subroutine fragopt_and_estimated_frags(nOcc,nvirt,OccOrbitals,virtOrbitals, &
       & MyMolecule,mylsitem,dofrag,esti,AtomicFragments,FragEnergies, &
       & AFset,dE_est3,Eskip_est)

    implicit none
    !> Full molecule info (CC model for each pair fragment resulting from estimate analysis is stored 
    !> in MyMolecule%ccmodel)
    type(fullmolecule), intent(inout) :: MyMolecule
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> Number of occupied orbitals in full molecule
    integer, intent(in) :: nocc
    !> Number of virtupied orbitals in full molecule
    integer, intent(in) :: nvirt
    !> Information about DEC occupied orbitals
    type(decorbital), dimension(nOcc), intent(in) :: OccOrbitals
    !> Information about DEC virtupied orbitals
    type(decorbital), dimension(nvirt), intent(in) :: virtOrbitals
    !> List of which atoms have orbitals assigned
    logical,intent(in),dimension(MyMolecule%nfrags) :: dofrag
    !> Is this a calculation with predefined small orbital spaces used to estimate frag energies?
    !> (This is effectively intent(in), only intent(inout) in practice for MPI purposes).
    logical,intent(inout) :: esti
    !> Atomic Fragments with orbital spaces determined to the FOT precision
    type(decfrag), intent(inout),dimension(MyMolecule%nfrags) :: AtomicFragments
    !> Fragment energies 
    real(realk),intent(inout) :: FragEnergies(MyMolecule%nfrags,MyMolecule%nfrags,ndecenergies)
    !> Have atomic fragments already been set by input (true) 
    !> or do they need to be optimized (false)
    !> NOTE: If AFset=false, nothing (except some MPI communication) is effectively done here.
    logical,intent(in) :: AFset
    !> error energy estimate based on information from the fragment optimization
    real(realk), intent(out) :: dE_est3, Eskip_est
    real(realk),pointer :: FragEnergiesPart(:,:)
    type(decfrag),pointer :: EstAtomicFragments(:)
    logical :: DoBasis,calcAF
    real(realk) :: init_radius,tcpu1,twall1,tcpu2,twall2,mastertime,Epair_est,deltaE_p
    integer :: nfrags,i,j,k
    type(joblist) :: jobs,fragoptjobs,estijobs
    integer(kind=ls_mpik) :: master


    master=0

    call LSTIMER('START',tcpu1,twall1,DECinfo%output)

    nfrags= MyMolecule%nfrags


    ! Initialize job list for atomic fragment optimizations
    call create_dec_joblist_fragopt(nfrags,nocc,nvirt,MyMolecule%ncore,MyMolecule%DistanceTable,&
         & OccOrbitals, virtOrbitals, dofrag, mylsitem,fragoptjobs)

    ! If atomic fragments have already been set, we do not carry out fragment optimization
    ! We therefore say that the job has already been done
    ! KK fixme - this could be done more elegantly...
    if(AFset) then
       do i=1,fragoptjobs%njobs
          fragoptjobs%jobsdone(i)=.true.
       end do
    end if

    if(DECinfo%DECrestart) then
       write(DECinfo%output,*) 'Restarting atomic fragment optimizations....'
       call restart_atomic_fragments_from_file(MyMolecule,MyLsitem,OccOrbitals,&
            & virtOrbitals,.false.,AtomicFragments,fragoptjobs)
    end if

    ! Sanity Check:
    if (DECinfo%DECNP .and. esti) call lsquit("ERROR(fragopt_and_estimated_frags): &
       &Pair estimates should be turned off in combination with DECNP!!",DECinfo%output)


    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '*** Initiating optimization of atomic fragments ***'
    write(DECinfo%output,*)
    write(DECinfo%output,*)


    UseEstimatedFragments: if(esti) then

       ! fragment opt. AND estimate pair fragments
       ! *****************************************

       ! Which model to use for pair estimate calculations?
       ! (MyMolecule%ccmodel is later redefined based on these estimates)
       do i=1,nfrags
          do j=i+1,nfrags
             MyMolecule%ccmodel(i,j) = DECinfo%PairEstimateModel
             MyMolecule%ccmodel(j,i) = DECinfo%PairEstimateModel
          end do
       end do

       !
       ! Init estimated atomic fragments by including orbitals assigned to neighbour atoms
       ! within 2 Angstrom.
       init_radius = DECinfo%EstimateINITradius
       DoBasis = .false.
       call mem_alloc(EstAtomicFragments,nfrags)
       call init_estimated_atomic_fragments(nOcc,nvirt,OccOrbitals,virtOrbitals, &
            & MyMolecule,mylsitem,DoBasis,init_radius,dofrag,EstAtomicFragments)

#ifdef VAR_MPI
       ! Send estimated fragment information to slaves
       call mpi_bcast_many_fragments(nfrags,dofrag,EstAtomicFragments,MPI_COMM_LSDALTON)

       ! Send CC models to use for each pair based on estimates (always MP2, but keep it general)
       call ls_mpibcast(MyMolecule%ccmodel,nfrags,nfrags,master,MPI_COMM_LSDALTON)
#endif

       ! Get job list for estimated pair fragments
       calcAF = .false.  ! No atomic fragments, just pairs
       call create_dec_joblist_driver(calcAF,MyMolecule,mylsitem,nfrags,nocc,nvirt,&
            &OccOrbitals,virtOrbitals,EstAtomicFragments,dofrag,esti,estijobs)

       if(DECinfo%DECrestart) then
          write(DECinfo%output,*) 'Restarting pair fragment estimate calculations...'
          call read_fragment_energies_for_restart_driver(nfrags,FragEnergies,estijobs,esti)
       end if

       ! Merge job list of atomic fragment optimization and estimated fragments (in this order)
       call concatenate_joblists(fragoptjobs,estijobs,jobs)

       call fragment_jobs(nocc,nvirt,nfrags,MyMolecule,mylsitem,OccOrbitals,&
            & virtOrbitals,jobs,AtomicFragments,FragEnergies,esti,&
            & fragoptjobs=fragoptjobs,estijobs=estijobs,EstAtomicFragments=EstAtomicFragments)
       do i=1,nfrags
          if(dofrag(i)) then
             call atomic_fragment_free_simple(EstAtomicFragments(i))
          end if
       end do
       call mem_dealloc(EstAtomicFragments)
    else

       ! Just fragment opt.
       ! ******************
       call fragment_jobs(nocc,nvirt,nfrags,MyMolecule,mylsitem,OccOrbitals,&
            & virtOrbitals,fragoptjobs,AtomicFragments,FragEnergies,esti)

    end if UseEstimatedFragments


    call evaluate_AtomicFragment_energy_error(AtomicFragments,nfrags,deltaE_p)


    if(esti) then
       ! Get estimated pair fragment energies for occupied partitioning scheme
       call mem_alloc(FragEnergiesPart,nfrags,nfrags)
       IF(DECinfo%onlyVirtPart)THEN
          call get_virtfragenergies(nfrags,DECinfo%PairEstimateModel,FragEnergies,FragEnergiesPart)
       ELSE
          call get_occfragenergies(nfrags,DECinfo%PairEstimateModel,FragEnergies,FragEnergiesPart)
       ENDIF
       ! We do not want to consider atomic fragment energies now so zero them
       ! (they might be zero already but in this way we avoid wrong print out below).
       do i=1,nfrags
          FragEnergiesPart(i,i)=0.0_realk
       end do
       ! Define which model to use in each pair calculation (info stored in MyMolecule%ccmodel)
       ! (Also calculate estimated correlation energy and estimated error by skipping pairs).
       call define_pair_calculations(nfrags,dofrag,FragEnergiesPart,MyMolecule,Epair_est,Eskip_est)
       call mem_dealloc(FragEnergiesPart)
    else 
       Eskip_est = 0.0E0_realk
    end if

    !Get the estimated total error of the calculation
    dE_est3 = deltaE_p + Eskip_est

    ! Zero all pair fragment energies to make sure that there are no leftovers from estimated pairs.
    do k=1,ndecenergies
       do j=1,nfrags
          do i=j+1,nfrags
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

  !> \author Patrick Ettenhuber
  !> extracting the approximate error made in a fragment calculation. we may use
  !this as a complementary measure to estimate the error made in a DEC
  !calculation
  subroutine evaluate_AtomicFragment_energy_error(AtomicFragments,nfrags,deltaE_p)
     implicit none
     integer, intent(in) :: nfrags
     type(decfrag), intent(in) :: AtomicFragments(nfrags)
     real(realk), intent(inout) :: deltaE_p
     integer :: frag, maxis
     real(realk) :: dE_occ, dE_vir, dE_Lag

     if(.not.DECinfo%InclFullMolecule .and. .not. DECinfo%simulate_full)then
        dE_occ = 0.0E0_realk
        dE_vir = 0.0E0_realk
        dE_Lag = 0.0E0_realk


        !Note, that %*_exp always contains the fragment energies of the expanded
        !fragment on the level of theory performed in the fragment reduction.
        !Therefore we may always get the difference here
        do frag=1,nfrags
           dE_occ = dE_occ + abs( AtomicFragments(frag)%Eocc_err )
           dE_vir = dE_vir + abs( AtomicFragments(frag)%Evir_err )
           dE_Lag = dE_Lag + abs( AtomicFragments(frag)%Elag_err )
        enddo

        maxis = 0
        if(DECinfo%OnlyOccPart)then
           maxis = 1
        else if(DECinfo%OnlyVirtPart)then
           maxis = 2
        else
           if(abs(dE_occ) > abs(dE_vir))then
              if( abs(dE_occ) > abs(dE_Lag))then
                 !occ > (virt and lag)
                 maxis = 1
              else
                 !lag > occ > virt
                 maxis = 3
              endif
           else
              if( abs(dE_vir) > abs(dE_Lag))then
                 !virt > (occ and lag)
                 maxis = 2
              else
                 !lag > virt > occ
                 maxis = 3
              endif
           endif
        endif


        deltaE_p =  0.0E0_realk
        select case(maxis)
        case(1)
           deltaE_p =  dE_occ
        case(2)
           deltaE_p =  dE_vir
        case(3)
           deltaE_p =  dE_Lag
        end select

        write(DECinfo%output,*)
        write(DECinfo%output,*)
        write(DECinfo%output,'(1X,a)') '******************************************************************'
        write(DECinfo%output,'(1X,a)') '*                SUMMARY ATOMIC FRAGMENT ERROR                   *'
        write(DECinfo%output,'(1X,a)') '******************************************************************'
        write(DECinfo%output,*)
        write(DECinfo%output,'(1X,a,i10)')    'Total number of fragments:                                      ', nfrags
        write(DECinfo%output,*)
        write(DECinfo%output,'(1X,a,g20.10)') 'Estimated contribution beyond the fragment energy contribution: ',deltaE_p
        if(DECinfo%PL>0)then
           write(DECinfo%output,*)
           write(DECinfo%output,'(1X,a,g20.10)') '-- individual estimated errors:'
           write(DECinfo%output,'(1X,a,g20.10)') '   occupied fragment energy:   ',dE_occ
           write(DECinfo%output,'(1X,a,g20.10)') '   virtual fragment energy:    ',dE_vir
           write(DECinfo%output,'(1X,a,g20.10)') '   Lagrangian fragment energy: ',dE_Lag
        endif
        write(DECinfo%output,*)
     else
        deltaE_p=0.0E0_realk
     endif
  end subroutine evaluate_AtomicFragment_energy_error

end module dec_driver_module


