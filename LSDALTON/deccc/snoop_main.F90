!> @file

!> Module to calculate interaction energies using SNOOP scheme.
!> \author Kasper Kristensen

module snoop_main_module

  use fundamental
  use precision
  use lstiming!,only: lstimer
  use typedeftype!,only: lsitem
  use matrix_module!,only:matrix,matrixp
  use files!,only:lsopen,lsclose
  use matrix_util!,only:util_MO_to_AO_different_trans
  use matrix_operations!,only: mat_dotproduct, mat_init, mat_zero, &
!       & mat_free, mat_set_from_full, mat_trans, mat_daxpy, mat_mul,&
!       & mat_add, mat_assign
  use matrix_operations_aux
  use memory_handling!,only: mem_alloc,mem_dealloc, mem_turnonthread_memory,&
!       & init_threadmemvar, collect_thread_memory, mem_turnoffthread_memory
  use dec_typedef_module
  use lsparameters
  use IntegralInterfaceMOD!,only: ii_get_twoelectron_gradient, ii_get_reorthonormalization, &
!       & ii_get_oneelectron_gradient, ii_get_nn_gradient
  use optimlocMOD, only: optimloc
  use configurationType
  use DALTONINFO


  ! DEC DEPENDENCIES (within deccc directory)   
  ! *****************************************
  use dec_fragment_utils
  use snoop_tools_module
  use crop_tools_module
  use ccintegrals
  use full_molecule
  use full, only: full_driver
  use dec_driver_module
  use array2_simple_operations
  use orbital_operations
  use atomic_fragment_operations
  use f12_routines_module

  public :: snoop_driver
  private

contains

  !> Driver for calculation interaction enegy using local orbitals. 
  !> \author Kasper Kristensen
  subroutine snoop_driver(Lsfull,config,MyMoleculeFULL,D)
    implicit none
    !> LSitem for full system
    type(lsitem), intent(inout) :: lsfull
    !> Config info
    type(configItem),intent(in)  :: config
    !> Molecule info for full system
    type(fullmolecule),intent(inout) :: MyMoleculeFULL
    !> Density matrix for full system
    type(matrix),intent(in) :: D
    real(realk),pointer :: S(:,:)

    ! Overlap matrix
    call mem_alloc(S,MyMoleculeFULL%nbasis,MyMoleculeFULL%nbasis)
    call II_get_mixed_overlap_full(DECinfo%output,DECinfo%output,lsfull%SETTING,&
         & S,MyMoleculeFULL%nbasis,MyMoleculeFULL%nbasis,AORdefault,AORdefault)

    call snoop_workhorse(Lsfull,config,MyMoleculeFULL,D,S)

    call mem_dealloc(S)

  end subroutine snoop_driver



  !> Driver for calculation interaction enegy using local orbitals. 
  !> \author Kasper Kristensen
  subroutine snoop_workhorse(Lsfull,config,MyMoleculeFULL,D,S)
    implicit none
    !> LSitem for full system
    type(lsitem), intent(inout) :: lsfull
    !> Config info
    type(configItem),intent(in)  :: config
    !> Molecule info
    type(fullmolecule),intent(inout) :: MyMoleculeFULL
    !> Density matrix for full system
    type(matrix),intent(in) :: D
    !> Overlap matrix
    real(realk),intent(in) :: S(MyMoleculeFULL%nbasis,MyMoleculeFULL%nbasis)
    real(realk),pointer :: EHFsnoop(:),Ecorrsnoop(:),EHFiso(:) 
    real(realk) :: EHFfull,Ecorrfull,Eerr
    integer :: nsub,nbasis,nvirtfull,noccfull,i,this,noccsnoop,nvirtsnoop,nbasissnoop,nelsnoop
    type(matrix),pointer :: Cocciso(:),Cvirtiso(:)
    type(matrix) :: FAOsnoop, FAOiso,Coccsnoop,Cvirtsnoop
    type(matrix) :: C
    type(lsitem) :: lssnoop,lsiso
    integer :: nocciso,nvirtiso,nbasisiso,neliso,nMO
    type(decorbital),pointer :: OccOrbitals(:), VirtOrbitals(:)
    real(realk) :: dummy(3,MyMoleculeFULL%natoms)  ! gradient input, just a dummy for now
    real(realk),pointer :: FragEnergiesOcc(:,:)
    type(decfrag),pointer :: AFfull(:)
    logical :: file_exist,onesub

    ! Sanity zeroing
    EHFfull=0.0_realk
    Ecorrfull=0.0_realk

    !> Calculate just one subsystem? 
    if(DECinfo%SNOOPonesub==-1) then
       ! Calculate all subsystems
       onesub=.false.
    else
       ! Calculate just one subsystem (possibly the full molecule)
       onesub=.true.
    end if

    ! Determine DEC orbital structures
    call mem_alloc(OccOrbitals,MyMoleculeFULL%nocc)
    call mem_alloc(VirtOrbitals,MyMoleculeFULL%nvirt)
    call GenerateOrbitals_driver(MyMoleculeFULL,lsfull,MyMoleculeFULL%nocc,MyMoleculeFULL%nvirt,&
         & MyMoleculeFULL%natoms, OccOrbitals, VirtOrbitals)
    call mem_alloc(FragEnergiesOcc,MyMoleculeFULL%nfrags,MyMoleculeFULL%nfrags) ! init frag energy array
    call mem_alloc(AFfull,MyMoleculeFULL%nfrags)

    ! Number of subsystems
    nsub = lsfull%input%molecule%nSubSystems

    ! HF energy and correlation energy for full molecular system
    ! **********************************************************

    ! Restart file exists?
    this=0
    if(DECinfo%SNOOPrestart) then
       call snoop_read_restart_file(this,EHFfull,Ecorrfull,file_exist)
    else
       file_exist=.false.
    end if

    RestartFull: if(file_exist) then
       write(DECinfo%output,*) 'SNOOP: Read full molecular energies from file!'
    else
       ! Restart file did not exist, we calculate full HF and corr energies from scratch

       OneSub1: if(onesub .and. DECinfo%SNOOPonesub/=0) then
          ! We are only interested in one subsystem, which is not the full molecule
          write(DECinfo%output,*) 'SNOOP: Skipping full molecular calculation!'

       else

          ! Perform full calculations
          if(DECinfo%full_molecular_cc) then
             ! Full calculation
             write(DECinfo%output,*) 'SNOOP: Starting full system calculation using full driver'
             call full_driver(MyMoleculeFULL,lsfull,D,EHFfull,Ecorrfull)
          else
             ! DEC calculation
             write(DECinfo%output,*) 'SNOOP: Starting full system calculation using DEC driver'
             call main_fragment_driver(MyMoleculeFULL,lsfull,D,&
                  & OccOrbitals,VirtOrbitals,MyMoleculeFULL%natoms,MyMoleculeFULL%nocc,MyMoleculeFULL%nvirt,&
                  & EHFfull,Ecorrfull,dummy,Eerr,FragEnergiesOcc,AFfull,.false.)
          end if

       end if OneSub1

    end if RestartFull


    ! F12 singles contribution
    if(DECinfo%F12 .and. DECinfo%F12singles) then
       ! Add F12 singles correction to HF energy
       EHFfull = EHFfull + MyMoleculeFULL%EF12singles
    endif


    ! Make restart file for full molecular energies
    this=0
    call snoop_write_restart_file(this,EHFfull,Ecorrfull)


    ! Notation
    ! --------
    ! "iso": Isolated subsystem (monomer) using only its own basis function
    ! "snoop": Subsystem (monomer) using also some effective "virtual ghost functions" (the SNOOP way)

    call mem_alloc(EHFsnoop,nsub)
    call mem_alloc(EHFiso,nsub)
    call mem_alloc(Ecorrsnoop,nsub)
    EHFsnoop = 0.0_realk  
    EHFiso = 0.0_realk  
    Ecorrsnoop=0.0_realk

    ! Number of basis functions and virt orbitals for full system
    nbasis = MyMoleculeFULL%nbasis
    nvirtfull = MyMoleculeFULL%nvirt
    noccfull = MyMoleculeFULL%nocc

    write(DECinfo%output,'(1X,a,i6,a)') 'Starting SNOOP subsystem calculations for ', nsub, &
         & ' subsystems.'


    ! **************************************************************************************
    ! Starting orbitals : HF calculations on isolated monomers OR localized full orbitals *
    ! **************************************************************************************
    call mem_alloc(Cocciso,nsub)
    call mem_alloc(Cvirtiso,nsub)

    ! Starting/reference orbitals come from isolated HF calculations

    IsoLOOP: do i=1,nsub
       this = i

       ! LSitem for SNOOP subsystem
       call build_subsystem_lsitem_no_ghost(this,MyMoleculeFULL,lsfull,lsiso)

       ! Number of electrons/orbitals for subsystem 
       ! ==========================================
       nbasisiso = lsiso%input%molecule%nbastREG
       neliso = get_num_electrons(lsiso)
       ! Currently SNOOP assumes closed-shell subsystem
       if(mod(neliso,2)/=0) then
          call lsquit('SNOOP only implemented for closed-shell systems',-1)
       else
          nocciso = neliso/2
       end if
       nvirtiso = nbasisiso-nocciso


       write(DECinfo%output,'(1X,a,i6,a,3i8)') 'ISOLATED SUBSYSTEM: ', this, &
            & ' -- #occ, #virt, #basis ', nocciso,nvirtiso,nbasisiso

       ! Initial occ and virt MO coefficients for subsystem - simple h1 diagonalization for now
       ! --> This should be made more efficient!
       call mat_init(Cocciso(this),nbasisiso,nocciso)
       call mat_init(Cvirtiso(this),nbasisiso,nvirtiso)
       call starting_orbitals_from_h1_split(lsiso,Cocciso(this),Cvirtiso(this))

       ! SCF optimization for isolated subsystem "this"
       call mat_init(FAOiso,nbasisiso,nbasisiso)
       write(DECinfo%output,*) 'Starting isolated subsystem SCF solver for subsystem', this
       call solve_subsystem_scf_rh(lsiso,Cocciso(this),Cvirtiso(this),FAOiso,EHFiso(this))

       print '(1X,a,i5,a,g20.12)', 'Isolated subsystem: ', this, &
            & '  *** HF energy: ', EHFiso(this)
       write(DECinfo%output,'(1X,a,i5,a,g20.12)') 'Isolated subsystem: ', this, &
            & '  *** HF energy: ', EHFiso(this)

       ! Free stuff for subsystem
       call mat_free(FAOiso)
       call ls_free(lsiso)

    end do IsoLOOP


    ! ****************************************************************************************
    ! HF calculations using AO basis from isolated monomer + virtual space from other monomers
    ! ****************************************************************************************

    ! SNOOP HF energy and correlation energy (if requested) for all subsystems
    call mat_init(FAOsnoop,nbasis,nbasis)
    SNOOPLOOP: do i=1,nsub

       this = i

       ! Restart file exists?
       if(DECinfo%SNOOPrestart) then
          call snoop_read_restart_file(this,EHFsnoop(this),Ecorrsnoop(this),file_exist)
       else
          file_exist=.false.
       end if

       ! If file exist, the HF and corr energies for this subsystem were read from file
       ! and we proceed to the next subsystem
       if(file_exist) then
          write(DECinfo%output,'(1X,a,i7,a)') 'SNOOP: SUBSYSTEM ', this, &
               & ' restart info read from file '
          cycle SNOOPLOOP
       end if

       ! Calculate just one subsystem?
       OneSub2: if(onesub .and. DECinfo%SNOOPonesub/=this) then
          ! We are only interested in one subsystem, and it is not this subsystem
          write(DECinfo%output,*) 'SNOOP: Skipping subsystem', this
          cycle SNOOPLOOP
       end if OneSub2


       IdenticalSubsystems: if(config%samesubsystems .and. this>1) then
          ! Same subsystems, we can skip calculation for all but subsystem 1 and 
          ! just copy HF and correlation energies.
          write(DECinfo%output,*) 'IDENTICAL SUBSYSTEMS - NO CALCULATION FOR SUBSYSTEM', this
          EHFsnoop(this) = EHFsnoop(1)
          Ecorrsnoop(this) = Ecorrsnoop(1)

       else

          write(DECinfo%output,*) 'STARTING SNOOP FOR SUBSYSTEM', this

          ! LSitem for subsystem using ghost functions on the other subsystems
          call build_subsystem_lsitem_ghost(this,lsfull,lssnoop)

          ! Number of electrons/orbitals for subsystem 
          ! ==========================================
          nelsnoop = get_num_electrons(lssnoop)
          ! Currently SNOOP assumes closed-shell subsystem
          if(mod(nelsnoop,2)/=0) then
             call lsquit('SNOOP only implemented for closed-shell systems with zero ',-1)
          else
             noccsnoop = nelsnoop/2
          end if


          ! Initial occ and virt MO coefficients for subsystem
          ! --------------------------------------------------
          ! 1. Take occ and virt MOs from isolated monomer calculations augmented 
          !    with zeros for basis functions on the other subsystems
          ! 2. Add virtual orbitals (Cvirtother) from other subsystems augmented with zeros
          !    for basis functions on this subsystem
          ! 3. Orthogonalize Cvirtother againts MOs on this subsystem while keeping occ and virt
          !    MOs for this subsystem fixed (they are already orthogonal).

          ! NOTE: Coccsnoop is initialized here, but Cvirtsnoop is
          !       initialized inside subroutine because we
          !       do not yet know the dimensions!

          call mat_init(Coccsnoop,nbasis,noccsnoop)

          call get_orthogonal_basis_for_subsystem(this,nsub,&
               & MyMoleculeFULL,Cocciso,Cvirtiso,Coccsnoop,Cvirtsnoop,S)

          ! Sanity check for initial orbitals
          call subsystem_orbitals_sanity_check(Coccsnoop,&
               & Cvirtsnoop,MyMoleculeFULL,S)

          ! SCF optimization for subsystem "this"
          call solve_subsystem_scf_rh(lssnoop,Coccsnoop,&
               & Cvirtsnoop,FAOsnoop,EHFsnoop(this))

          ! Sanity check for optimized orbitals
          call subsystem_orbitals_sanity_check(Coccsnoop,&
               & Cvirtsnoop,MyMoleculeFULL,S)

          ! Determine orbitals for correlated SNOOP monomer calculations

          ! SNOOP-DEC stuff
          DECcalc: if(.not. DECinfo%full_molecular_cc) then
             if(DECinfo%SNOOPlocalize) then
                ! Localize orbitals
                ! -----------------
                nMO = Coccsnoop%ncol + Cvirtsnoop%ncol
                call mat_init(C,nbasis,nMO)
                call collect_MO_coeff_in_one_matrix(Coccsnoop,Cvirtsnoop,C) 
                call optimloc(C,noccsnoop,config%decomp%cfg_mlo_m,lssnoop,&
                     & config%davidOrbLoc)
                call partition_MO_coeff_into_two_matrices(C,Coccsnoop,Cvirtsnoop)
                call mat_free(C)
             else
                ! Rotate subsystem orbitals using natural connection such that they are as
                ! close as possible to the full orbitals in a least-squares sense.
                call rotate_subsystem_orbitals_to_mimic_FULL_orbitals(MyMoleculeFULL,this,&
                     & OccOrbitals,VirtOrbitals,lssnoop,Coccsnoop, Cvirtsnoop,S)
             end if
          end if DECcalc

          ! Correlation energy for subsystem
          if(.not. DECinfo%SNOOPjustHF) then
             call subsystem_correlation_energy(this,MyMoleculeFULL,OccOrbitals,VirtOrbitals,AFfull,&
                  & Coccsnoop,Cvirtsnoop,FAOsnoop,lssnoop,EHFsnoop(this),Ecorrsnoop(this))
          end if

          ! Free stuff for subsystem
          call ls_free(lssnoop)
          call mat_free(Coccsnoop)
          call mat_free(Cvirtsnoop)


       end if IdenticalSubsystems

       ! Make restart file
       call snoop_write_restart_file(this,EHFsnoop(this),Ecorrsnoop(this))

       print '(1X,a,i5,a,3g20.12)', 'SNOOP subsystem: ', this, &
            & '  *** HF/corr/HFdiff energy: ', EHFsnoop(this), Ecorrsnoop(this),EHFsnoop(this)-EHFiso(this)

       write(DECinfo%output,'(1X,a,i5,a,3g20.12)') 'SNOOP subsystem: ', this, &
            & '  *** HF/corr/HFdiff energy: ', EHFsnoop(this), Ecorrsnoop(this),EHFsnoop(this)-EHFiso(this)


    End do SNOOPLOOP


    ! Print interaction energy summary
    call SNOOP_interaction_energy_print(nsub,EHFsnoop,Ecorrsnoop,EHFfull,Ecorrfull)

    call mat_free(FAOsnoop)
    call mem_dealloc(EHFsnoop)
    call mem_dealloc(EHFiso)
    call mem_dealloc(Ecorrsnoop)
    do i=1,nsub
       call mat_free(Cocciso(i))
       call mat_free(Cvirtiso(i))
    end do
    call mem_dealloc(Cocciso)
    call mem_dealloc(Cvirtiso)

    do i=1,MyMoleculeFULL%nocc
       call orbital_free(OccOrbitals(i))
    end do
    do i=1,MyMoleculeFULL%nvirt
       call orbital_free(VirtOrbitals(i))
    end do
    call mem_dealloc(OccOrbitals)
    call mem_dealloc(VirtOrbitals)
    call mem_dealloc(FragEnergiesOcc)

    do i=1,MyMoleculeFULL%nfrags
       if(.not. associated(AFfull(i)%EOSatoms)) cycle
       call atomic_fragment_free_simple(AFfull(i))
    end do
    call mem_dealloc(AFfull)


  end subroutine snoop_workhorse




  !> \brief Solve SCF equations for subsystem to get optimized occ and virt
  !> orbitals Cocc and Cvirt for subsystem using simple RH/DIIS scheme.
  !> The equations to solve is that the virt-occ Fock matrix elements should be zero:
  !>
  !> Fsub(Dsub)_ai = 0          i : occupied in this subsystem, a: virtual
  !>
  !> Dsub_{mu,nu} = sum_{i \in subsystem} C_{mu i} C_{nu i}
  !> Fsub is Fock matrix where only interactions between electrons and nuclei in
  !> subsystem are included.
  !>
  !> Subsystem is defined solely by lssub input, so in principle this subroutine
  !> could also be used for the full molecular system (i.e., one subsystem which is
  !> the full molecule), although this is not the intention.
  !>
  !> \author Kasper Kristensen
  subroutine solve_subsystem_scf_rh(lssub,CoccAO,CvirtAO,FAOsub,EHF)

    implicit none
    !> Integral info
    type(lsitem), intent(inout) :: lssub
    !> Initial occupied MO coefficients for subsystem
    !> At output: Optimized occupied MOs for subsystem
    type(matrix), intent(inout) :: CoccAO
    !> Inital virtual MO coefficients
    !> At output: Optimized virtual MOs
    type(matrix), intent(inout) :: CvirtAO
    !> AO Fock matrix for density corresponding to optimized orbitals
    type(matrix),intent(inout) :: FAOsub
    !> Hf energy for converged MOs for subsystem
    real(realk),intent(inout) :: EHF
    integer :: nocc, nvirt,nbasis
    type(matrix), pointer :: residual(:),FAO(:)
    type(matrix) :: FAOopt,DAO,CMOi,FMOi,SMOi,CAO,CAOold
    reaL(realk) :: tcpu,twall
    real(realk), pointer :: B(:,:),c(:),eival(:)
    logical :: converged,dodiis
    real(realk) :: prev_norm, resnorm, convthr,diisnorm
    integer :: iter, last_iter, i,j,n


    ! ************************************************************
    !               Note on orbital basis notation
    ! ************************************************************
    ! We use two different bases:
    !
    ! AO  : Atomic orbital basis
    ! MOi : MO basis defined by input orbitals


    call LSTIMER('START',tcpu,twall,DECinfo%output)



    ! Initialize stuff
    ! ****************
    nocc   = CoccAO%ncol         
    nvirt = CvirtAO%ncol         
    nbasis = CoccAO%nrow         
    n = nocc + nvirt             
    ! Note: n < nbasis if routine is for SNOOP subsystem, 
    !       n=nbasis if routine is used for isolated subsystem (with no ghost functions).


    ! CROP/DIIS matrices
    call mem_alloc(B,DECinfo%SNOOPMaxIter,DECinfo%SNOOPMaxIter)
    call mem_alloc(c,DECinfo%SNOOPMaxIter)

    ! Fock matrices in AO basis and residuals for each iteration
    call mem_alloc(FAO,DECinfo%SNOOPMaxIter)
    call mem_alloc(residual,DECinfo%SNOOPMaxIter)

    ! Other matrices
    call mat_init(DAO,nbasis,nbasis)
    call mat_init(FMOi,n,n)
    call mat_init(CMOi,n,n)
    call mat_init(CAO,nbasis,n)
    call mat_init(SMOi,n,n)    ! overlap matrix in MOi basis = identity matrix
    call mat_identity(SMOi)
    call mem_alloc(eival,n)
    call mat_init(FAOopt,nbasis,nbasis)

    ! MO coeff for initial MOs expressed in AO basis,
    ! i.e., CAOold define MOi basis
    call mat_init(CAOold,nbasis,n)
    call collect_MO_coeff_in_one_matrix(CoccAO,CvirtAO,CAOold) 


    ! Solver information
    dodiis=.true.
    if(DECinfo%SNOOPMaxDIIS==0) dodiis=.false.
    convthr = DECinfo%SNOOPTHR      ! The convergence threshold
    converged=.false.
    prev_norm = huge(1.0_realk)

    write(DECinfo%output,*)
    write(DECinfo%output,*) 'Starting subsystem SCF iterations'
    write(DECinfo%output,*) '---------------------------------'
    write(DECinfo%output,'(1X,a)')  '###  Iteration     Residual norm            HF energy'
    write(6,'(1X,a)')  '###  Iteration     Residual norm            HF energy'


    ! Start Roothaan-Hall iterations
    ! ******************************

    RHiteration : do iter=1,DECinfo%SNOOPMaxIter

       ! Current iteration
       last_iter = iter

       ! remove old vectors
       RemoveOldVectors : if(iter > DECinfo%SNOOPMaxDIIS) then
          call mat_free(residual(iter-DECinfo%SNOOPMaxDIIS))
          call mat_free(FAO(iter-DECinfo%SNOOPMaxDIIS))
       end if RemoveOldVectors


       ! Check for convergence of current MOs
       ! ************************************

       ! AO Density matrix for current occupied orbitals
       call get_density_from_occ_orbitals_mat(CoccAO,DAO)

       ! Fock matrix in AO basis for current density
       call mat_init(FAO(iter),nbasis,nbasis)
       call dec_fock_transformation(FAO(iter),DAO,Lssub,.true.,incl_h=.true.)

       ! HF energy
       EHF = get_HF_energy(DAO,FAO(iter),Lssub)

       ! Residual = Virt-occ block of Fock matrix with current MO orbitals (before diagonalization),
       ! i.e., the orbitals corresponding to the density from which FAO(iter) was constructed
       call mat_init(residual(iter),nvirt,nocc)
       call util_AO_to_MO_different_trans(CvirtAO,FAO(iter),CoccAO,residual(iter))

       ! Residual norm
       resnorm = mat_sqnorm2(residual(iter))  ! norm squared
       resnorm = sqrt(resnorm)                ! norm itself

       write(DECinfo%output,'(1X,a,2X,i4,5X,g20.10,3X,g22.12)')  '### ',iter, resnorm,EHF
       write(*,'(1X,a,2X,i4,5X,g20.10,3X,g22.12)')  '### ',iter, resnorm,EHF

       ! Sanity warning
       if(resnorm > prev_norm) then
          write(DECinfo%output,'(a)') ' warning :: norm was smaller in previous iteration !!! '
       end if
       prev_norm=resnorm

       ! check if this is the last iteration
       if(resnorm < convthr) converged=.true.
       if( (iter == DECinfo%SNOOPMaxIter) .or. converged) then
          call mat_assign(FAOsub,FAO(iter))
          exit RHiteration
       end if



       ! DIIS ACCELERATION
       ! *****************

       DIISacc: if(dodiis) then

          ! calculate diis matrix
          B=0.0E0_realk; c=0.0E0_realk
          do i=iter,max(iter-DECinfo%SNOOPMaxDIIS+1,1),-1
             do j=iter,i,-1
                B(i,j) = mat_dotproduct(residual(i),residual(j))
                B(j,i) = B(i,j)
             end do
          end do

          ! solve crop/diis equation
          call CalculateDIIScoefficients(DECinfo%SNOOPMaxDIIS,DECinfo%SNOOPMaxIter,iter,B,c, &
               DECinfo%SNOOPdebug)

          ! mixing to get optimal Fock matrix
          call mat_zero(FAOopt)
          do i=iter,max(iter-DECinfo%SNOOPMaxDIIS+1,1),-1
             call mat_daxpy(c(i),FAO(i),FAOopt)
          end do

       else
          ! No DIIS: Optimal Fock matrix is just the one from current iteration
          call mat_assign(FAOopt, FAO(iter))
       end if DIISacc


       ! Diagonalize Fock matrix and get corresponding MOs
       ! *************************************************

       ! Construct Fock matrix in MOi basis 
       ! (for dimensionality reasons we use the "different_trans" subroutine
       !  even though we actually transform with CAOold both from left and right)
       call util_AO_to_MO_different_trans(CAOold,FAOopt,CAOold,FMOi)

       ! Diagonalize Fock matrix in MOi basis
       ! NOTE: It is crucial that this is done in MOi basis and NOT in AO basis
       !       to ensure that we do not mix with occ orbitals on other subsystems!
       call mat_diag_f(FMOi,SMOi,eival,CMOi)

       ! CMOi are MO coeffecients expressed in the MOi basis, we want them in AO basis.
       ! We get this by multiplying with initial MO coefficients expressed in AO basis:
       ! "New MOs in AO basis" = "MOis in AO basis" times "Current MOs in MOi basis"
       ! CAO = CAOold CMOi
       call mat_mul(CAOold, CMOi, 'N', 'N', 1.0_realk, 0.0_realk, CAO)

       ! Put new MO coefficients into CoccAO and CvirtAO
       call partition_MO_coeff_into_two_matrices(CAO,CoccAO,CvirtAO)

    end do RHiteration


    if(converged) then
       write(DECinfo%output,'(A41,I5,A12)') &
            & 'SCF subsystem equation solved in', last_iter, ' iterations!'
    else
       call lsquit('solve_subsystem_scf: Equations not solved!',-1)
    endif

    ! deallocate stuff
    call mem_dealloc(B)
    call mem_dealloc(c)
    call mem_dealloc(eival)
    do i=last_iter,max(last_iter-DECinfo%SNOOPMaxDIIS+1,1),-1
       call mat_free(FAO(i))
       call mat_free(residual(i))
    end do
    call mem_dealloc(residual)
    call mem_dealloc(FAO)
    call mat_free(DAO)
    call mat_free(FMOi)
    call mat_free(CMOi)
    call mat_free(CAO)
    call mat_free(SMOi)
    call mat_free(FAOopt)
    call mat_free(CAOold)

    call LSTIMER('SUBSYSTEM SCF',tcpu,twall,DECinfo%output)

  end subroutine solve_subsystem_scf_rh


  !> Check that subsystem orbitals are properly orthogonal and normalized.
  subroutine subsystem_orbitals_sanity_check(Coccsub_mat,&
       & Cvirtsub_mat,MyMoleculeFULL,S)
    implicit none

    !> Occ and virt MO coefficients 
    type(matrix),intent(in) :: Coccsub_mat, Cvirtsub_mat
    !> Full molecule info
    type(fullmolecule),intent(in) :: MyMoleculeFULL
    !> Overlap matrix
    real(realk),intent(in) :: S(MyMoleculeFULL%nbasis,MyMoleculeFULL%nbasis)
    real(realk),pointer :: tmp(:,:),Coccsub(:,:),Cvirtsub(:,:)
    integer :: nocc,nvirt,nbasis,i,j
    real(realk) :: thr,one

    thr = 1.0E-8_realk
    one = 1.0_realk

    ! Dimensions
    nbasis = Coccsub_mat%nrow
    nocc = Coccsub_mat%ncol
    nvirt = Cvirtsub_mat%ncol

    ! Work with Fortran arrays - occ and virt MOs for subsystem
    call mem_alloc(Coccsub,nbasis,nocc)
    call mem_alloc(Cvirtsub,nbasis,nvirt)
    call mat_to_full(Coccsub_mat, 1.0_realk, Coccsub)
    call mat_to_full(Cvirtsub_mat, 1.0_realk, Cvirtsub)


    ! 1. Check overlap between Coccsub and Cvirtsub: tmp = Coccsub^T SAO Cvirtsub
    call mem_alloc(tmp,nocc,nvirt)
    call dec_diff_basis_transform1(nbasis,nocc,nvirt,&
         & Coccsub,Cvirtsub,S,tmp)
    do j=1,nvirt
       do i=1,nocc
          if(abs(tmp(i,j))>thr) then
             print *, 'Check1: i,j,value,thr',i,j,tmp(i,j),thr
             call lsquit('subsystem_orbitals_sanity_check1: Orbitals are not orthogonal!',-1)
          end if
       end do
    end do
    call mem_dealloc(tmp)


    ! 2. Check that Coccsub orbitals are orthonormal
    call mem_alloc(tmp,nocc,nocc)
    call dec_simple_basis_transform1(nbasis,nocc,Coccsub,S,tmp)
    do j=1,nocc

       do i=1,nocc

          if(i==j) then

             ! Is orbital "i" normalized?
             if(abs(one-tmp(i,j)) > thr) then
                print *, 'Check2a: i,j,value,thr',i,j,one-tmp(i,j),thr
                call lsquit('subsystem_orbitals_sanity_check2: Orbitals are not normalized!',-1)
             end if

          else

             ! Are orbitals "i" and "j" orthogonal
             if(abs(tmp(i,j)) > thr) then
                print *, 'Check2b: i,j,value,thr',i,j,tmp(i,j),thr
                call lsquit('subsystem_orbitals_sanity_check2: Orbitals are not orthogonal!',-1)
             end if

          end if

       end do

    end do
    call mem_dealloc(tmp)



    ! 3. Check that Cvirtsub orbitals are orthonormal
    call mem_alloc(tmp,nvirt,nvirt)
    call dec_simple_basis_transform1(nbasis,nvirt,Cvirtsub,S,tmp)
    do j=1,nvirt

       do i=1,nvirt

          if(i==j) then

             ! Is orbital "i" normalized?
             if(abs(one-tmp(i,j)) > thr) then
                print *, 'Check3a: i,j,value,thr',i,j,one-tmp(i,j),thr
                call lsquit('subsystem_orbitals_sanity_check3: Orbitals are not normalized!',-1)
             end if

          else

             ! Are orbitals "i" and "j" orthogonal
             if(abs(tmp(i,j)) > thr) then
                print *, 'Check3b: i,j,value,thr',i,j,tmp(i,j),thr
                call lsquit('subsystem_orbitals_sanity_check3: Orbitals are not orthogonal!',-1)
             end if

          end if

       end do

    end do
    call mem_dealloc(tmp)


    call mem_dealloc(Coccsub)
    call mem_dealloc(Cvirtsub)

  end subroutine subsystem_orbitals_sanity_check


  !> Print energy summary for SNOOP interaction energy calculation
  subroutine SNOOP_interaction_energy_print(nsub,EHFsub,Ecorrsub,&
       & EHFfull,Ecorrfull)
    implicit none
    !> Number of subsystems
    integer,intent(in) :: nsub
    !> HF energies for subsystems
    real(realk),intent(in) :: EHFsub(nsub)
    !> Correlation energies for subsystems
    real(realk),intent(in) :: Ecorrsub(nsub)
    !> HF and correlation energy for total system
    real(realk),intent(in) :: EHFfull,Ecorrfull
    real(realk) :: EHFint,Ecorrint
    integer :: i

    write(DECinfo%output,'(1X,a)') ''
    write(DECinfo%output,'(1X,a)') ''
    write(DECinfo%output,'(1X,a)') '***************************************************************'
    write(DECinfo%output,'(1X,a)') '            SNOOP interaction energies - summary               '
    write(DECinfo%output,'(1X,a)') '***************************************************************'
    write(DECinfo%output,'(1X,a,g22.12)') 'Full system     --- HF   energy : ', EHFfull
    if(.not. DECinfo%SNOOPjustHF) then
       write(DECinfo%output,'(1X,a,g22.12)') 'Full system     --- corr energy : ', Ecorrfull
       write(DECinfo%output,'(1X,a,g22.12)') 'Full system     --- tot  energy : ', EHFfull+Ecorrfull
    end if
    EHFint = EHFfull
    Ecorrint = Ecorrfull
    do i=1,nsub
       write(DECinfo%output,'(1X,a)') ''
       write(DECinfo%output,'(1X,a,i5,a,g21.12)') 'Subsystem: ', i, &
            & ' --- HF   energy: ', EHFsub(i)
       if(.not. DECinfo%SNOOPjustHF) then
          write(DECinfo%output,'(1X,a,i5,a,g21.12)') 'Subsystem: ', i, &
               & ' --- corr energy: ', Ecorrsub(i)
          write(DECinfo%output,'(1X,a,i5,a,g21.12)') 'Subsystem: ', i, &
               & ' --- tot  energy: ', EHFsub(i)+Ecorrsub(i)
       end if
       EHFint = EHFint - EHFsub(i)
       Ecorrint = Ecorrint - Ecorrsub(i)
    end do
    write(DECinfo%output,'(1X,a)') ''
    write(DECinfo%output,'(1X,a)') '---------------------------------------------------------------'
    write(DECinfo%output,'(1X,a,g22.12)') 'HF Interaction energy      = ', EHFint
    if(.not. DECinfo%SNOOPjustHF) then
       write(DECinfo%output,'(1X,a,g22.12)') 'Corr Interaction energy    = ', Ecorrint
       write(DECinfo%output,'(1X,a,g22.12)') 'Total Interaction energy   = ', EHFint+Ecorrint
    end if
    write(DECinfo%output,'(1X,a)') '---------------------------------------------------------------'
    write(DECinfo%output,'(1X,a)') ''
    write(DECinfo%output,'(1X,a)') ''


  end subroutine SNOOP_interaction_energy_print





  !> Calculate correlation energy for subsystem. F12 singles correction is also added to HF energy
  !> if requested.
  !> \author Kasper Kristensen
  !> \date October 2014
  subroutine subsystem_correlation_energy(this,MyMoleculeFULL,OccOrbitalsFULL,&
       & VirtOrbitalsFULL,AFfull,Cocc,Cvirt,F,lssub,EHF,Ecorr)
    implicit none

    !> Which subsystem
    integer,intent(in) :: this
    !> Molecule info for full system
    type(fullmolecule),intent(in) :: MyMoleculeFULL
    !> Occ orbitals for full molecule
    type(decorbital),intent(in) :: OccOrbitalsFULL(MyMoleculeFULL%nocc)
    !> Virt orbitals for full molecule
    type(decorbital),intent(in) :: VirtOrbitalsFULL(MyMoleculeFULL%nvirt)
    !> Atomic fragments for full molecule
    type(decfrag),intent(in) :: AFfull(MyMoleculeFULL%nfrags)
    !> Occupied and virtual MO coefficients for subsystem
    type(matrix),intent(in) :: Cocc, Cvirt
    !> Fock matrix in AO basis for subsystem
    type(matrix),intent(in) :: F
    !> LSitem for subsystem
    type(lsitem), intent(inout) :: lssub
    !> Subsystem HF energy (F12 singles correction is added if requested)
    real(realk),intent(inout) :: EHF
    !> Subsystem correlation energy
    real(realk),intent(inout) :: Ecorr
    type(matrix) ::D,C
    type(fullmolecule) :: MySubsystem
    real(realk) :: dummyE,Eerr,simple_orbital_threshold_save
    logical :: simulate_full_save,InclFullMolecule_save,full_molecular_cc_save
    integer :: nMO,nbasis,i
    real(realk),pointer :: dummyG(:,:)

    ! Dimensions
    nbasis = F%nrow
    nMO= Cocc%ncol + Cvirt%ncol     ! #MOs = #occ + #virt   (different from nbasis in general)


    ! Make fullmolecule structure for subsystem
    ! *****************************************

    ! Density matrix for subsystem
    call mat_init(D,nbasis,nbasis)
    call get_density_from_occ_orbitals_mat(Cocc,D)

    ! Collect MO coefficients in one matrix
    call mat_init(C,nbasis,nMO)
    call collect_MO_coeff_in_one_matrix(Cocc,Cvirt,C)

    ! Molecule structure for subsystem
    call molecule_init_from_inputs(MySubsystem,lssub,F,C,D)
    call mat_free(C)

    ! F12 singles contribution
    if(DECinfo%F12 .and. DECinfo%F12singles ) then
       call F12singles_driver(MySubsystem,lssub,D)
       ! Add F12 singles correction to HF energy
       EHF = EHF + MySubsystem%EF12singles
    endif


    ! Correlation energy for subsystem
    ! ********************************
    WhichScheme: if(DECinfo%full_molecular_cc) then

       ! Full calculation
       write(DECinfo%output,'(1X,a,i7,a)') 'SNOOP: Starting subsystem ', this, &
            & ' calculation using full driver'

       if(DECinfo%SNOOPdecfrag) then

          ! Consider SNOOP monomer to be one DEC fragment
          ! *********************************************

          ! Quick and dirty solution to ensure that we only 
          ! have one fragment representing the SNOOP monomer.
          simple_orbital_threshold_save = DECinfo%simple_orbital_threshold
          simulate_full_save = DECinfo%simulate_full
          InclFullMolecule_save = DECinfo%InclFullMolecule
          full_molecular_cc_save = DECinfo%full_molecular_cc
          DECinfo%simulate_full = .true.
          DECinfo%simple_orbital_threshold = -0.1_realk
          DECinfo%InclFullMolecule = .true.
          DECinfo%full_molecular_cc = .false.
          call mem_alloc(dummyG,3,MySubsystem%natoms) 

          ! DEC calculation with simulatefull where the one fragment
          ! is the SNOOP monomer
          call DEC_wrapper(MySubsystem,lssub,D,dummyE,Ecorr,dummyG,Eerr)

          ! Free and reset
          call mem_dealloc(dummyG)
          DECinfo%simulate_full = simulate_full_save
          DECinfo%simple_orbital_threshold = simple_orbital_threshold_save
          DECinfo%InclFullMolecule = InclFullMolecule_save 
          DECinfo%full_molecular_cc = full_molecular_cc_save 

       else
          call full_driver(MySubsystem,lssub,D,dummyE,Ecorr)
       end if

    else

       ! DEC calculation
       call DECsubsystem_correlation_energy(this,MyMoleculeFULL,OccOrbitalsFULL,&
            & VirtOrbitalsFULL,AFfull,Cocc,Cvirt,F,D,MySubsystem,lssub,Ecorr)

    end if WhichScheme


    call mat_free(D)
    call molecule_finalize(MySubsystem,.true.)

  end subroutine subsystem_correlation_energy


  !> Calculate DEC correlation energy for subsystem u
  !> \author Kasper Kristensen
  !> \date October 2014
  subroutine DECsubsystem_correlation_energy(sub,MyMoleculeFULL,OccOrbitalsFULL,&
       & VirtOrbitalsFULL,AFfull,Cocc,Cvirt,F,D,MySubsystem,lssub,Ecorr)
    implicit none

    !> Which subsystem
    integer,intent(in) :: sub
    !> Molecule info for full system
    type(fullmolecule),intent(in) :: MyMoleculeFULL
    !> Occ orbitals for full molecule
    type(decorbital),intent(in) :: OccOrbitalsFULL(MyMoleculeFULL%nocc)
    !> Virt orbitals for full molecule
    type(decorbital),intent(in) :: VirtOrbitalsFULL(MyMoleculeFULL%nvirt)
    !> Atomic fragments for full molecule
    type(decfrag),intent(in) :: AFfull(MyMoleculeFULL%nfrags)
    !> Occupied and virtual MO coefficients for subsystem
    type(matrix),intent(in) :: Cocc, Cvirt
    !> Fock matrix in AO basis for subsystem
    type(matrix),intent(in) :: F
    !> Density matrix in AO basis for subsystem
    type(matrix),intent(in) :: D
    !> Full molecule structure for subsystem
    type(fullmolecule),intent(inout) :: MySubsystem
    !> LSitem for subsystem
    type(lsitem), intent(inout) :: lssub
    !> Subsystem correlation energy
    real(realk),intent(inout) :: Ecorr
    type(decfrag),pointer :: AFsub(:)
    logical :: dofragSUB(MySubsystem%nfrags), dofragFULL(MyMoleculeFULL%nfrags)
    real(realk) :: EHF,Eerr
    integer :: i
    real(realk),pointer :: dummy(:,:),FragEnergiesOcc(:,:)
    type(decorbital),pointer :: OccOrbitalsSUB(:), VirtOrbitalsSUB(:)
    real(realk) :: energies(ndecenergies)


    ! DEC calculation
    ! ---------------

    call mem_alloc(dummy,3,MySubsystem%natoms) ! gradient input, just a dummy for now

    if(DECinfo%SNOOPsamespace) then
       ! Use "same" orbital space for subsystem as for full system, as defined 
       ! by natural connection

       write(DECinfo%output,'(1X,a,i7,a)') 'SNOOP: Starting subsystem ', sub, &
            & ' calculation using full DEC space'

       ! Fragment energies
       call mem_alloc(FragEnergiesOcc,MySubsystem%nfrags,MySubsystem%nfrags)
       FragEnergiesOcc=0.0_realk

       ! Determine DEC orbital structures
       call mem_alloc(OccOrbitalsSUB,MySubsystem%nocc)
       call mem_alloc(VirtOrbitalsSUB,MySubsystem%nvirt)
       call GenerateOrbitals_driver(MySubsystem,lssub,MySubsystem%nocc,MySubsystem%nvirt,&
            & MySubsystem%natoms, OccOrbitalsSUB, VirtOrbitalsSUB)

       ! Check that orbital assignment for subsystem is consistent with that for full system
       call Orbitals_subsystem_vs_FULL_sanity_check(sub,MyMoleculeFULL,MySubsystem,&
            & OccOrbitalsFULL,VirtOrbitalsFULL,OccOrbitalsSUB,VirtOrbitalsSUB)

       !  List of which fragments to consider for subsystem
       call which_fragments_to_consider(MySubsystem%ncore,MySubsystem%nocc,MySubsystem%nvirt,&
            & MySubsystem%nfrags,OccOrbitalsSUB,VirtOrbitalsSUB,dofragSUB,MySubsystem%PhantomAtom)

       !  List of which fragments to consider for full system
       call which_fragments_to_consider(MyMoleculeFULL%ncore,MyMoleculeFULL%nocc,&
            & MyMoleculeFULL%nvirt,MyMoleculeFULL%nfrags,OccOrbitalsFULL,VirtOrbitalsFULL,&
            & dofragFULL,MyMoleculeFULL%PhantomAtom)

       ! Get atomic fragments for subsystem with one-to-one correspondence to
       ! full atomic fragments
       call mem_alloc(AFsub,MySubsystem%nfrags)
       do i=1,MySubsystem%nfrags
          call atomic_fragment_nullify(AFsub(i))
       end do
       call subsystemAOS_equals_FULLAOS(sub,MyMoleculeFULL,&
            & OccOrbitalsFULL,MySubsystem,lssub,OccOrbitalsSUB,&
            & VirtOrbitalsSUB,dofragFULL,dofragsub,AFfull,AFsub)

       ! Run DEC fragment calculations with atomic fragment defined above
        call main_fragment_driver(MySubsystem,lssub,D,OccOrbitalsSUB,&
            & VirtOrbitalsSUB,MySubsystem%natoms,MySubsystem%nocc,&
            & MySubsystem%nvirt,EHF,Ecorr,dummy,Eerr,FragEnergiesOcc,AFsub,.true.)


       ! Free stuff
       do i=1,MySubsystem%nocc
          call orbital_free(OccOrbitalsSUB(i))
       end do
       do i=1,MySubsystem%nvirt
          call orbital_free(VirtOrbitalsSUB(i))
       end do
       call mem_dealloc(OccOrbitalsSUB)
       call mem_dealloc(VirtOrbitalsSUB)
       call mem_dealloc(FragEnergiesOcc)
       do i=1,MySubsystem%nfrags
          if(.not. associated(AFsub(i)%EOSatoms)) cycle
          call atomic_fragment_free_simple(AFsub(i))
       end do
       call mem_dealloc(AFsub)

    else

       ! Do independent DEC fragment optimization for monomer
       write(DECinfo%output,'(1X,a,i7,a)') 'SNOOP: Starting subsystem ', sub, &
            & ' calculation using DEC driver'
       call DEC_wrapper(MySubsystem,lssub,D,EHF,Ecorr,dummy,Eerr)

    end if

    call mem_dealloc(dummy)
    call molecule_finalize(MySubsystem,.true.)

  end subroutine DECsubsystem_correlation_energy


end module snoop_main_module
