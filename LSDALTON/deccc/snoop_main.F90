!> @file

!> Module to handle interaction energy calculations using local orbitals.
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
  use IntegralInterfaceMOD!,only: ii_get_twoelectron_gradient, ii_get_reorthonormalization, &
!       & ii_get_oneelectron_gradient, ii_get_nn_gradient


  ! DEC DEPENDENCIES (within deccc directory)   
  ! *****************************************
  use dec_fragment_utils
  use snoop_tools_module
  use crop_tools_module
  use full_molecule
  use full
!!$  use array2_simple_operations!,only: array2_init, array2_transpose, array2_free,&
!!$ !      & extract_occupied_EOS_MO_indices
!!$  use array4_simple_operations!,only: array4_init, array4_contract1, array4_free,&
!!$ !      & array4_reorder, array4_contract3
!!$  use ccintegrals!,only:dec_fock_transformation
!!$  use orbital_operations
!!$  use ccdriver
!!$  use full
!!$  use mp2_module
!!$  use fragment_energy_module
!!$  use atomic_fragment_operations
!!$  use tensor_type_def_module
!!$  use full_molecule
!!$  use dec_driver_module
!!$  use mp2_gradient_module
  

  public :: snoop_driver
  private

contains

  !> Driver for calculation interaction enegy using local orbitals. 
  !> \author Kasper Kristensen
  subroutine snoop_driver(Lsfull,MyMolecule,D)
    implicit none
    !> LSitem for full system
    type(lsitem), intent(inout) :: lsfull
    !> Molecule info
    type(fullmolecule),intent(inout) :: MyMolecule
    !> Density matrix for full system
    type(matrix),intent(in) :: D
!!$    type(decorbital), pointer :: OccOrbitals(:),VirtOrbitals(:)
!!$    integer :: this,i,nsub,noccfull,nvirtfull,nbasis,noccsub,nvirtsub,natomssub,nelsub,nbasissub,j
!!$    real(realk),pointer :: EHFsub(:),Ecorrsub(:),EHFiso(:)

!!$    type(matrix) :: FAOsub,tmp
!!$
!!$
    real(realk),pointer :: EHFsnoop(:),Ecorrsnoop(:),EHFiso(:) 
    real(realk) :: EHFfull,Ecorrfull
    integer :: nsub,nbasis,nvirtfull,noccfull,i,this,noccsnoop,nvirtsnoop,nbasissnoop,nelsnoop
    type(matrix),pointer :: Cocciso(:),Cvirtiso(:),Coccsnoop(:),Cvirtsnoop(:)
    type(lsitem) :: lssnoop


    ! Number of subsystems
    nsub = lsfull%input%molecule%nSubSystems
    if(nsub /= 2) then
       print *, 'Number of subsystems = ',nsub
       call lsquit('snoop_driver: Only implemented for two subsystems!',-1)
    end if

    ! HF energy and correlation energy for full molecular system
    call full_driver(MyMolecule,lsfull,D,EHFfull,Ecorrfull)

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
    nbasis = MyMolecule%nbasis
    nvirtfull = MyMolecule%nunocc
    noccfull = MyMolecule%nocc

    write(DECinfo%output,'(1X,a,i6,a)') 'Starting SNOOP subsystem calculations for ', nsub, &
         & ' subsystems.'


    ! ************************************
    ! HF calculations on isolated monomers
    ! ************************************
    call mem_alloc(Cocciso,nsub)
    call mem_alloc(Cvirtiso,nsub)
    call mem_alloc(Coccsnoop,nsub)
    call mem_alloc(Cvirtsnoop,nsub)

    ! HF energy for all subsystems
    do i=1,nsub

       this = i

       ! LSitem for SNOOP subsystem
       call build_subsystem_lsitem_no_ghost(this,MyMolecule,lsfull,lssnoop)
       
       ! Number of electrons/orbitals for subsystem 
       ! ==========================================
       nbasissnoop = lssnoop%input%molecule%nbastREG
       nelsnoop = get_num_electrons(lssnoop)
       ! Currently SNOOP assumes closed-shell subsystem
       if(mod(nelsnoop,2)/=0) then
          call lsquit('SNOOP only implemented for closed-shell systems with zero ',-1)
       else
          noccsnoop = nelsnoop/2
       end if
       nvirtsnoop = nbasissnoop-noccsnoop


       write(DECinfo%output,'(1X,a,i6,a,3i8)') 'SUBSYSTEM: ', this, &
            & ' -- #occ, #virt, #basis ', noccsnoop,nvirtsnoop,nbasissnoop

       ! Initial occ and virt MO coefficients for subsystem - simple h1 diagonalization for now
       ! --> This should be made more efficient!
       call mat_init(Cocciso(this),nbasissnoop,noccsnoop)
       call mat_init(Cvirtiso(this),nbasissnoop,nvirtsnoop)
       call starting_orbitals_from_h1_split(lssnoop,Cocciso(this),Cvirtiso(this))
!!$
!!$       ! SCF optimization for isolated subsystem "this"
!!$       call mat_init(FAOsnoop,nbasissnoop,nbasissnoop)
!!$       call solve_subsystem_scf_rh(this,lssnoop,Cocciso(this),Cvirtiso(this),FAOsnoop,EHFiso(this))
!!$
!!$       print '(1X,a,i5,a,g20.12)', 'Isolated subsystem: ', this, &
!!$            & '  *** HF energy: ', EHFiso(this)
!!$       write(DECinfo%output,'(1X,a,i5,a,g20.12)') 'Isolated subsystem: ', this, &
!!$            & '  *** HF energy: ', EHFiso(this)
!!$
!!$
!!$       if(DECinfo%hack3) then
!!$          call get_mp2_density_for_subsystem(this,lssnoop,Cocciso(this),Cvirtiso(this),FAOsnoop,&
!!$               & CoccMP2(this),CvirtMP2(this))
!!$       end if
!!$
!!$
!!$       ! Free stuff for subsystem
!!$       call mat_free(FAOsnoop)
!!$       call ls_free(lssnoop)

    end do
!!$
!!$
!!$
!!$    ! ****************************************************************************************
!!$    ! HF calculations using AO basis from isolated monomer + virtual space from other monomers
!!$    ! ****************************************************************************************
!!$
!!$    SNOOPSCF: do j=1,2
!!$
!!$       print *, '******** SNOOPSCF round ',j
!!$
!!$       ! HF energy and correlation energy (if requested) for all subsystems
!!$       do i=1,nsub
!!$
!!$          this = i
!!$
!!$          ! LSitem for subsystem using ghost functions on the other subsystems
!!$          call build_subsystem_lsitem_ghost(this,lsfull,lssnoop)
!!$
!!$          ! Number of electrons/orbitals for subsystem 
!!$          ! (NOTE: assumes closed-shell and uncharged subsystem!)
!!$          nbasissnoop = lssnoop%input%molecule%nbastREG
!!$          nelsnoop = get_num_electrons(lssnoop,zerocharge=.true.)
!!$          noccsnoop = nelsnoop/2
!!$
!!$          ! Initial occ and virt MO coefficients for subsystem
!!$          ! --------------------------------------------------
!!$          ! 1. Take occ and virt MOs from isolated monomer calculations augmented 
!!$          !    with zeros for basis functions on the other subsystems
!!$          ! 2. Add virtual orbitals (Cvirtother) from other subsystems augmented with zeros
!!$          !    for basis functions on this subsystem
!!$          ! 3. Orthogonalize Cvirtother againts MOs on this subsystem while keeping occ and virt
!!$          !    MOs for this subsystem fixed (they are already orthogonal).
!!$
!!$          ! NOTE: Coccsnoop is initialized here, but Cvirtsnoop is
!!$          !       initialized inside subroutine because we
!!$          !       do not yet know the dimensions!
!!$
!!$          if(j==1) then
!!$             call mat_init(Coccsnoop(i),nbasis,noccsnoop)
!!$          else
!!$             call mat_free(Cvirtsnoop(i))
!!$          end if
!!$
!!$          if(j==1) then
!!$             call get_orthogonal_basis_for_subsystem(this,nsub,&
!!$                  & MyMolecule,Cocciso,Cvirtiso,Coccsnoop(i),Cvirtsnoop(i),.false.)
!!$          else
!!$             call get_orthogonal_basis_for_subsystem(this,nsub,&
!!$                  & MyMolecule,Coccsnoop,Cvirtsnoop,Coccsnoop(i),Cvirtsnoop(i),.false.)
!!$          end if
!!$             
!!$          ! call mat_print_all('COCCSUB1',Coccsub,DECinfo%output)
!!$          ! call mat_print_all('CVIRTSUB1',Cvirtsub,DECinfo%output)
!!$
!!$          ! Sanity check for orbitals
!!$          call subsystem_orbitals_sanity_check(Coccsnoop(i),&
!!$               & Cvirtsnoop(i),MyMolecule)
!!$
!!$          ! SCF optimization for subsystem "this"
!!$          call mat_init(FAOsnoop,nbasis,nbasis)
!!$          call solve_subsystem_scf_rh(this,lssnoop,Coccsnoop(i),&
!!$               & Cvirtsnoop(i),FAOsnoop,EHFsnoop(this))
!!$
!!$          ! call mat_print_all('COCCSUB2',Coccsub,DECinfo%output)
!!$          ! call mat_print_all('CVIRTSUB2',Cvirtsub,DECinfo%output)
!!$
!!$          call subsystem_orbitals_sanity_check(Coccsnoop(i),&
!!$               & Cvirtsnoop(i),MyMolecule)
!!$
!!$          ! Correlation energy for subsystem
!!$          if(.not. DECinfo%hack2) then
!!$             call subsystem_correlation_energy(this,MyMolecule,&
!!$                  & OccOrbitals,Coccsnoop(i),Cvirtsnoop(i),&
!!$                  & FAOsnoop,lssnoop,Ecorrsnoop(this))
!!$          end if
!!$
!!$          print '(1X,a,i5,a,3g20.12)', 'Subsystem: ', this, &
!!$               & '  *** HF/corr/HFdiff energy: ', EHFsnoop(this), Ecorrsnoop(this),EHFsnoop(this)-EHFiso(this)
!!$
!!$          write(DECinfo%output,'(1X,a,i5,a,3g20.12)') 'Subsystem: ', this, &
!!$               & '  *** HF/corr/HFdiff energy: ', EHFsnoop(this), Ecorrsnoop(this),EHFsnoop(this)-EHFiso(this)
!!$
!!$
!!$          ! Free stuff for subsystem
!!$          call mat_free(FAOsnoop)
!!$          !       call mat_free(Coccsnoop)
!!$          !       call mat_free(Cvirtsnoop)
!!$          call ls_free(lssnoop)
!!$
!!$       End do
!!$
!!$    end do SNOOPSCF
!!$
!!$
!!$    ! Print interaction energy summary
!!$    call local_interaction_energy_print(nsub,EHFsnoop,Ecorrsnoop,EHFfull,Ecorrfull)
!!$
!!$    do i=1,MyMolecule%nocc
!!$       call orbital_free(OccOrbitals(i))
!!$    end do
!!$    do i=1,MyMolecule%nunocc
!!$       call orbital_free(VirtOrbitals(i))
!!$    end do
!!$    call mem_dealloc(OccOrbitals)
!!$    call mem_dealloc(VirtOrbitals)
!!$    call mem_dealloc(EHFsnoop)
!!$    call mem_dealloc(EHFiso)
!!$    call mem_dealloc(Ecorrsnoop)
!!$
!!$    do i=1,nsub
!!$       call mat_free(Cocciso(i))
!!$       call mat_free(Cvirtiso(i))
!!$       call mat_free(CoccMP2(i))
!!$       call mat_free(CvirtMP2(i))
!!$       call mat_free(Coccsnoop(i))
!!$       call mat_free(Cvirtsnoop(i))
!!$    end do
!!$    call mem_dealloc(CoccMP2)
!!$    call mem_dealloc(CvirtMP2)
!!$    call mem_dealloc(Cocciso)
!!$    call mem_dealloc(Cvirtiso)
!!$    call mem_dealloc(Coccsnoop)
!!$    call mem_dealloc(Cvirtsnoop)

  end subroutine snoop_driver




!!$  !> \brief Solve SCF equations for subsystem to get optimized occ and virt
!!$  !> orbitals Cocc and Cvirt for subsystem.
!!$  !> The equations to solve is that the virt-occ Fock matrix elements should be zero:
!!$  !>
!!$  !> Fsub(Dsub)_ai = 0          i : occupied in this subsystem, a: virtual
!!$  !>
!!$  !> Dsub_{mu,nu} = sum_{i \in subsystem} C_{mu i} C_{nu i}
!!$  !> Fsub is Fock matrix where only interactions between electrons and nuclei in
!!$  !> subsystem are included.
!!$  !>
!!$  !> We do not include occupied MO coefficients from other subsystems,
!!$  !> such that our final MO coefficients have no components from other subsystems.
!!$  !> 
!!$  !> \author Kasper Kristensen
!!$  subroutine solve_subsystem_scf_rh(this,Lssub,CoccAO,CvirtAO,FAOsub,EHF)
!!$
!!$    implicit none
!!$    !> Which subsystem to consider
!!$    integer,intent(in) :: this
!!$    !> Integral info
!!$    type(lsitem), intent(inout) :: lssub
!!$    !> Initial occupied MO coefficients ONLY for subsystem "this"
!!$    !> At output: Optimized occupied MOs for subsystem "this"
!!$    type(matrix), intent(inout) :: CoccAO
!!$    !> Inital virtual MO coefficients for all subsystems 
!!$    !> At output: Optimized virtual MOs for subsystem "this"
!!$    type(matrix), intent(inout) :: CvirtAO
!!$    !> AO Fock matrix for density corresponding to optimized orbitals
!!$    type(matrix),intent(inout) :: FAOsub
!!$    !> Hf energy for converged MOs for subsystem
!!$    real(realk),intent(inout) :: EHF
!!$    integer :: nocc, nvirt,nbasis
!!$    type(matrix), pointer :: residual(:),FAO(:)
!!$    type(matrix) :: FAOopt,DAO,CLO,FLO,SLO,CAO,CAOold
!!$    reaL(realk) :: tcpu,twall
!!$    real(realk), pointer :: B(:,:),c(:),eival(:)
!!$    logical :: converged,dodiis
!!$    real(realk) :: prev_norm, resnorm, convthr,diisnorm
!!$    integer :: iter, last_iter, i,j,n
!!$
!!$
!!$    ! ************************************************************
!!$    !               Note on orbital basis notation
!!$    ! ************************************************************
!!$    ! We use two different bases:
!!$    !
!!$    ! AO  : Atomic orbital basis
!!$    ! LO : Localized molecular orbital for Large system basis
!!$
!!$
!!$    call LSTIMER('START',tcpu,twall,DECinfo%output)
!!$    write(DECinfo%output,*) 'Starting subsystem SCF solver for subsystem', this
!!$
!!$
!!$    ! Initialize stuff
!!$    ! ****************
!!$    nocc   = CoccAO%ncol          ! number of occupied orbitals for subsystem
!!$    nvirt = CvirtAO%ncol         ! number of virtupied orbitals for full molecule
!!$    nbasis = CoccAO%nrow          ! number of basis functions (same for subsystem and full molecule)
!!$    n = nocc + nvirt             ! Number of occ for subsystem + total number of virt
!!$    ! Note: n < nbasis if routine is used for subsystem, 
!!$    !       n=nbasis if routine is used for full molecule
!!$
!!$
!!$    ! CROP/DIIS matrices
!!$    call mem_alloc(B,DECinfo%kappaMaxIter,DECinfo%kappaMaxIter)
!!$    call mem_alloc(c,DECinfo%kappaMaxIter)
!!$
!!$    ! Fock matrices in AO basis and residuals for each iteration
!!$    call mem_alloc(FAO,DECinfo%kappaMaxIter)
!!$    call mem_alloc(residual,DECinfo%kappaMaxIter)
!!$
!!$    ! Other matrices
!!$    call mat_init(DAO,nbasis,nbasis)
!!$    call mat_init(FLO,n,n)
!!$    call mat_init(CLO,n,n)
!!$    call mat_init(CAO,nbasis,n)
!!$    call mat_init(SLO,n,n)    ! overlap matrix in LO basis = identity matrix
!!$    call mat_identity(SLO)
!!$    call mem_alloc(eival,n)
!!$    call mat_init(FAOopt,nbasis,nbasis)
!!$
!!$    ! MO coeff for localized MOs expressed in AO basis,
!!$    ! i.e., CAOold define LO basis
!!$    call mat_init(CAOold,nbasis,n)
!!$    call collect_MO_coeff_in_one_matrix(CoccAO,CvirtAO,CAOold) 
!!$
!!$
!!$    ! Solver information
!!$    dodiis=.true.
!!$    convthr = DECinfo%kappaTHR      ! The convergence threshold
!!$    converged=.false.
!!$    prev_norm = huge(1.0_realk)
!!$
!!$    write(DECinfo%output,*)
!!$    write(DECinfo%output,*) 'Starting subsystem SCF iterations'
!!$    write(DECinfo%output,*) '---------------------------------'
!!$    write(DECinfo%output,'(1X,a)')  '###  Iteration     Residual norm            HF energy'
!!$    write(6,'(1X,a)')  '###  Iteration     Residual norm            HF energy'
!!$
!!$
!!$    ! Start Roothaan-Hall iterations
!!$    ! ******************************
!!$
!!$    RHiteration : do iter=1,DECinfo%kappaMaxIter
!!$
!!$       ! Current iteration
!!$       last_iter = iter
!!$
!!$       ! remove old vectors
!!$       RemoveOldVectors : if(iter > DECinfo%kappaMaxDIIS) then
!!$          call mat_free(residual(iter-DECinfo%kappaMaxDIIS))
!!$          call mat_free(FAO(iter-DECinfo%kappaMaxDIIS))
!!$       end if RemoveOldVectors
!!$
!!$
!!$       ! Check for convergence of current MOs
!!$       ! ************************************
!!$
!!$       ! AO Density matrix for current occupied orbitals
!!$       call get_density_from_occ_orbitals_mat(CoccAO,DAO)
!!$
!!$       ! Fock matrix in AO basis for current density
!!$       call mat_init(FAO(iter),nbasis,nbasis)
!!$       call dec_fock_transformation(FAO(iter),DAO,Lssub,.true.,incl_h=.true.)
!!$
!!$       ! HF energy
!!$       EHF = get_HF_energy(DAO,FAO(iter),Lssub)
!!$
!!$       ! Residual = Virt-occ block of Fock matrix with current MO orbitals (before diagonalization),
!!$       ! i.e., the orbitals corresponding to the density from which FAO(iter) was constructed
!!$       call mat_init(residual(iter),nvirt,nocc)
!!$       call util_AO_to_MO_different_trans(CvirtAO,FAO(iter),CoccAO,residual(iter))
!!$
!!$       ! Residual norm
!!$       resnorm = mat_sqnorm2(residual(iter))  ! norm squared
!!$       resnorm = sqrt(resnorm)                ! norm itself
!!$
!!$       write(DECinfo%output,'(1X,a,2X,i4,5X,g20.10,3X,g22.12)')  '### ',iter, resnorm,EHF
!!$       write(6,'(1X,a,2X,i4,5X,g20.10,3X,g22.12)')  '### ',iter, resnorm,EHF
!!$
!!$       ! Sanity warning
!!$       if(resnorm > prev_norm) then
!!$          write(DECinfo%output,'(a)') ' warning :: norm was smaller in previous iteration !!! '
!!$       end if
!!$       prev_norm=resnorm
!!$
!!$       ! check if this is the last iteration
!!$       if(resnorm < convthr) converged=.true.
!!$       if( (iter == DECinfo%kappaMaxIter) .or. converged) then
!!$          call mat_assign(FAOsub,FAO(iter))
!!$          exit RHiteration
!!$       end if
!!$
!!$
!!$
!!$       ! DIIS ACCELERATION
!!$       ! *****************
!!$
!!$       DIISacc: if(dodiis) then
!!$
!!$          ! calculate diis matrix
!!$          B=0.0E0_realk; c=0.0E0_realk
!!$          do i=iter,max(iter-DECinfo%kappaMaxDIIS+1,1),-1
!!$             do j=iter,i,-1
!!$                B(i,j) = mat_dotproduct(residual(i),residual(j))
!!$                B(j,i) = B(i,j)
!!$             end do
!!$          end do
!!$
!!$          ! solve crop/diis equation
!!$          call CalculateDIIScoefficients(DECinfo%kappaMaxDIIS,DECinfo%kappaMaxIter,iter,B,c, &
!!$               DECinfo%kappa_driver_debug)
!!$
!!$          ! mixing to get optimal Fock matrix
!!$          call mat_zero(FAOopt)
!!$          do i=iter,max(iter-DECinfo%kappaMaxDIIS+1,1),-1
!!$             call mat_daxpy(c(i),FAO(i),FAOopt)
!!$          end do
!!$
!!$       else
!!$          ! No DIIS: Optimal Fock matrix is just the one from current iteration
!!$          call mat_assign(FAOopt, FAO(iter))
!!$       end if DIISacc
!!$
!!$
!!$       ! Diagonalize Fock matrix and get corresponding MOs
!!$       ! *************************************************
!!$
!!$       ! Construct Fock matrix in LO basis 
!!$       ! (for dimensionality reasons we use the "different_trans" subroutine
!!$       !  even though we actually transform with CAOold both from left and right)
!!$       call util_AO_to_MO_different_trans(CAOold,FAOopt,CAOold,FLO)
!!$
!!$       ! Diagonalize Fock matrix in LO basis
!!$       ! NOTE: It is crucial that this is done in LO basis and NOT in AO basis
!!$       !       to ensure that we do not mix with occ orbitals on other subsystems!
!!$       call mat_diag_f(FLO,SLO,eival,CLO)
!!$
!!$       ! CLO are MO coeffecients expressed in the LO basis, we want them in AO basis.
!!$       ! We get this by multiplying with initial MO coefficients expressed in AO basis:
!!$       ! "New MOs in AO basis" = "LOs in AO basis" times "Current MOs in LO basis"
!!$       ! CAO = CAOold CLO
!!$       call mat_mul(CAOold, CLO, 'N', 'N', 1.0_realk, 0.0_realk, CAO)
!!$
!!$       ! Put new MO coefficients into CoccAO and CvirtAO
!!$       call partition_MO_coeff_into_two_matrices(CAO,CoccAO,CvirtAO)
!!$
!!$    end do RHiteration
!!$
!!$
!!$    if(converged) then
!!$       write(DECinfo%output,'(A41,I5,A12)') &
!!$            & 'SCF subsystem equation solved in', last_iter, ' iterations!'
!!$    else
!!$       call lsquit('solve_subsystem_scf: Equations not solved!',-1)
!!$    endif
!!$
!!$    ! deallocate stuff
!!$    call mem_dealloc(B)
!!$    call mem_dealloc(c)
!!$    call mem_dealloc(eival)
!!$    do i=last_iter,max(last_iter-DECinfo%kappaMaxDIIS+1,1),-1
!!$       call mat_free(FAO(i))
!!$       call mat_free(residual(i))
!!$    end do
!!$    call mem_dealloc(residual)
!!$    call mem_dealloc(FAO)
!!$    call mat_free(DAO)
!!$    call mat_free(FLO)
!!$    call mat_free(CLO)
!!$    call mat_free(CAO)
!!$    call mat_free(SLO)
!!$    call mat_free(FAOopt)
!!$    call mat_free(CAOold)
!!$
!!$    call LSTIMER('SUBSYSTEM SCF',tcpu,twall,DECinfo%output)
!!$
!!$  end subroutine solve_subsystem_scf_rh
!!$
!!$
!!$  !> Check that subsystem orbitals are orthogonal.
!!$  subroutine subsystem_orbitals_sanity_check(Coccsub_mat,&
!!$       & Cvirtsub_mat,MyMolecule)
!!$    implicit none
!!$
!!$    !> Occ and virt MO coefficients 
!!$    type(matrix),intent(in) :: Coccsub_mat, Cvirtsub_mat
!!$    !> Full molecule info
!!$    type(fullmolecule),intent(in) :: MyMolecule
!!$    real(realk),pointer :: tmp(:,:),Coccsub(:,:),Cvirtsub(:,:)
!!$    integer :: nocc,nvirt,nbasis,i,j
!!$    real(realk) :: thr,one
!!$    type(matrix) :: A
!!$
!!$    thr = 1.0E-8_realk
!!$    one = 1.0_realk
!!$
!!$    ! Dimensions
!!$    nbasis = Coccsub_mat%nrow
!!$    nocc = Coccsub_mat%ncol
!!$    nvirt = Cvirtsub_mat%ncol
!!$
!!$    ! Work with Fortran arrays - occ and virt MOs for subsystem
!!$    call mem_alloc(Coccsub,nbasis,nocc)
!!$    call mem_alloc(Cvirtsub,nbasis,nvirt)
!!$    call mat_to_full(Coccsub_mat, 1.0_realk, Coccsub)
!!$    call mat_to_full(Cvirtsub_mat, 1.0_realk, Cvirtsub)
!!$
!!$
!!$    ! 1. Check overlap between Coccsub and Cvirtsub: tmp = Coccsub^T SAO Cvirtsub
!!$    call mem_alloc(tmp,nocc,nvirt)
!!$    call dec_diff_basis_transform1(nbasis,nocc,nvirt,&
!!$         & Coccsub,Cvirtsub,MyMolecule%overlap,tmp)
!!$
!!$    ! HACK
!!$    call mat_init(A,nocc,nvirt)
!!$    call mat_set_from_full(tmp,1.0_realk,A)
!!$
!!$    call mat_free(A)
!!$    ! END HACK
!!$
!!$    do j=1,nvirt
!!$       do i=1,nocc
!!$          if(abs(tmp(i,j))>thr) then
!!$             print *, 'Check1: i,j,value,thr',i,j,tmp(i,j),thr
!!$             call lsquit('subsystem_orbitals_sanity_check1: Orbitals are not orthogonal!',-1)
!!$          end if
!!$       end do
!!$    end do
!!$    call mem_dealloc(tmp)
!!$
!!$
!!$    ! 2. Check that Coccsub orbitals are orthonormal
!!$    call mem_alloc(tmp,nocc,nocc)
!!$    call dec_simple_basis_transform1(nbasis,nocc,Coccsub,MyMolecule%overlap,tmp)
!!$
!!$    ! HACK
!!$    call mat_init(A,nocc,nocc)
!!$    call mat_set_from_full(tmp,1.0_realk,A)
!!$
!!$    call mat_free(A)
!!$    ! END HACK
!!$
!!$    do j=1,nocc
!!$
!!$       do i=1,nocc
!!$
!!$          if(i==j) then
!!$
!!$             ! Is orbital "i" normalized?
!!$             if(abs(one-tmp(i,j)) > thr) then
!!$                print *, 'Check2a: i,j,value,thr',i,j,one-tmp(i,j),thr
!!$                call lsquit('subsystem_orbitals_sanity_check2: Orbitals are not normalized!',-1)
!!$             end if
!!$
!!$          else
!!$
!!$             ! Are orbitals "i" and "j" orthogonal
!!$             if(abs(tmp(i,j)) > thr) then
!!$                print *, 'Check2b: i,j,value,thr',i,j,tmp(i,j),thr
!!$                call lsquit('subsystem_orbitals_sanity_check2: Orbitals are not orthogonal!',-1)
!!$             end if
!!$
!!$          end if
!!$
!!$       end do
!!$
!!$    end do
!!$    call mem_dealloc(tmp)
!!$
!!$
!!$
!!$    ! 3. Check that Cvirtsub orbitals are orthonormal
!!$    call mem_alloc(tmp,nvirt,nvirt)
!!$    call dec_simple_basis_transform1(nbasis,nvirt,Cvirtsub,MyMolecule%overlap,tmp)
!!$
!!$    ! HACK
!!$    call mat_init(A,nvirt,nvirt)
!!$    call mat_set_from_full(tmp,1.0_realk,A)
!!$
!!$    call mat_free(A)
!!$    ! END HACK
!!$
!!$    do j=1,nvirt
!!$
!!$       do i=1,nvirt
!!$
!!$          if(i==j) then
!!$
!!$             ! Is orbital "i" normalized?
!!$             if(abs(one-tmp(i,j)) > thr) then
!!$                print *, 'Check3a: i,j,value,thr',i,j,one-tmp(i,j),thr
!!$                write(DECinfo%output,*) 'SNOOPWARNING3: Orbitals are not normalized!'
!!$!                call lsquit('subsystem_orbitals_sanity_check3: Orbitals are not normalized!',-1)
!!$             end if
!!$
!!$          else
!!$
!!$             ! Are orbitals "i" and "j" orthogonal
!!$             if(abs(tmp(i,j)) > thr) then
!!$                print *, 'Check3b: i,j,value,thr',i,j,tmp(i,j),thr
!!$                write(DECinfo%output,*) 'SNOOPWARNING3: Orbitals are not orthogonal!'
!!$!                call lsquit('subsystem_orbitals_sanity_check3: Orbitals are not orthogonal!',-1)
!!$             end if
!!$
!!$          end if
!!$
!!$       end do
!!$
!!$    end do
!!$    call mem_dealloc(tmp)
!!$
!!$
!!$    call mem_dealloc(Coccsub)
!!$    call mem_dealloc(Cvirtsub)
!!$
!!$  end subroutine subsystem_orbitals_sanity_check
!!$
!!$
!!$
!!$
!!$  !> Check that subsystem orbitals are normalized and orthogonal to orbitals from other subsystems.
!!$  subroutine subsystem_orbitals_sanity_check_old(this,Coccsub_mat,Cvirtsub_mat,MyMolecule,OccOrbitals)
!!$    implicit none
!!$
!!$    !> Subsystem label
!!$    integer,intent(in) :: this
!!$    !> Occ and virt MO coefficients 
!!$    type(matrix),intent(in) :: Coccsub_mat, Cvirtsub_mat
!!$    !> Full molecule info
!!$    type(fullmolecule),intent(in) :: MyMolecule
!!$    !> Occ orbitals in DEC format
!!$    type(decorbital),intent(in) :: OccOrbitals(MyMolecule%nocc)
!!$    real(realk),pointer :: tmp(:,:),Coccsub(:,:),Cvirtsub(:,:)
!!$    integer :: noccfull,nvirt,noccsub,nbasis,i,j,atomj,subj
!!$    real(realk) :: thr,one
!!$
!!$    thr = 1.0E-10_realk
!!$    one = 1.0_realk
!!$
!!$    ! Dimensions
!!$    nbasis = MyMolecule%nbasis
!!$    noccfull = MyMolecule%nocc   ! Number of occ orbitals in full molecule
!!$    noccsub = Coccsub_mat%ncol       ! Number of occ orbitals in subsystem
!!$    nvirt = MyMolecule%nunocc    ! Number of virt orbitals (same for subsystem and full system)
!!$
!!$    ! Work with Fortran arrays - occ and virt MOs for subsystem
!!$    call mem_alloc(Coccsub,nbasis,noccsub)
!!$    call mem_alloc(Cvirtsub,nbasis,nvirt)
!!$    call mat_to_full(Coccsub_mat, 1.0_realk, Coccsub)
!!$    call mat_to_full(Cvirtsub_mat, 1.0_realk, Cvirtsub)
!!$
!!$
!!$    ! 1. Check overlap between Coccsub and Cvirtsub: tmp = Coccsub^T SAO Cvirtsub
!!$    call mem_alloc(tmp,noccsub,nvirt)
!!$    call dec_diff_basis_transform1(nbasis,noccsub,nvirt,&
!!$         & Coccsub,Cvirtsub,MyMolecule%overlap,tmp)
!!$    do j=1,nvirt
!!$       do i=1,noccsub
!!$          if(abs(tmp(i,j))>thr) then
!!$             print *, 'Check1: i,j,value,thr',i,j,tmp(i,j),thr
!!$             call lsquit('subsystem_orbitals_sanity_check1: Orbitals are not orthogonal!',-1)
!!$          end if
!!$       end do
!!$    end do
!!$    call mem_dealloc(tmp)
!!$
!!$    ! 2. Check overlap between Coccsub and occ orbitals for other subsystems
!!$    call mem_alloc(tmp,noccsub,noccfull)
!!$    call dec_diff_basis_transform1(nbasis,noccsub,noccfull,&
!!$         & Coccsub,MyMolecule%Co,MyMolecule%overlap,tmp)
!!$    do j=1,noccfull
!!$
!!$       ! Subsystem to which orbital "j" is assigned                  
!!$       atomj = OccOrbitals(j)%centralatom
!!$       subj = MyMolecule%SubSystemIndex(atomj)
!!$
!!$       ! Only consider orbital "j" if assigned to another subsystem
!!$       if(subj/=this) then
!!$          do i=1,noccsub
!!$             if(abs(tmp(i,j))>thr) then
!!$                print *, 'Check2: i,j,value,thr',i,j,tmp(i,j),thr
!!$                call lsquit('subsystem_orbitals_sanity_check2: Orbitals are not orthogonal!',-1)
!!$             end if
!!$          end do
!!$       end if
!!$
!!$    end do
!!$    call mem_dealloc(tmp)
!!$
!!$
!!$
!!$    ! 3. Check that Coccsub orbitals are orthonormal
!!$    call mem_alloc(tmp,noccsub,noccsub)
!!$    call dec_simple_basis_transform1(nbasis,noccsub,Coccsub,MyMolecule%overlap,tmp)
!!$    do j=1,noccsub
!!$
!!$       do i=1,noccsub
!!$
!!$          if(i==j) then
!!$
!!$             ! Is orbital "i" normalized?
!!$             if(abs(one-tmp(i,j)) > thr) then
!!$                print *, 'Check3a: i,j,value,thr',i,j,one-tmp(i,j),thr
!!$                call lsquit('subsystem_orbitals_sanity_check3: Orbitals are not normalized!',-1)
!!$             end if
!!$
!!$          else
!!$
!!$             ! Are orbitals "i" and "j" orthogonal
!!$             if(abs(tmp(i,j)) > thr) then
!!$                print *, 'Check3b: i,j,value,thr',i,j,tmp(i,j),thr
!!$                call lsquit('subsystem_orbitals_sanity_check3: Orbitals are not orthogonal!',-1)
!!$             end if
!!$
!!$          end if
!!$
!!$       end do
!!$
!!$    end do
!!$    call mem_dealloc(tmp)
!!$
!!$
!!$
!!$    ! 4. Check that Cvirtsub orbitals are orthonormal
!!$    call mem_alloc(tmp,nvirt,nvirt)
!!$    call dec_simple_basis_transform1(nbasis,nvirt,Cvirtsub,MyMolecule%overlap,tmp)
!!$    do j=1,nvirt
!!$
!!$       do i=1,nvirt
!!$
!!$          if(i==j) then
!!$
!!$             ! Is orbital "i" normalized?
!!$             if(abs(one-tmp(i,j)) > thr) then
!!$                print *, 'Check4a: i,j,value,thr',i,j,one-tmp(i,j),thr
!!$                call lsquit('subsystem_orbitals_sanity_check4: Orbitals are not normalized!',-1)
!!$             end if
!!$
!!$          else
!!$
!!$             ! Are orbitals "i" and "j" orthogonal
!!$             if(abs(tmp(i,j)) > thr) then
!!$                print *, 'Check4b: i,j,value,thr',i,j,tmp(i,j),thr
!!$                call lsquit('subsystem_orbitals_sanity_check4: Orbitals are not orthogonal!',-1)
!!$             end if
!!$
!!$          end if
!!$
!!$       end do
!!$
!!$    end do
!!$    call mem_dealloc(tmp)
!!$
!!$
!!$
!!$
!!$
!!$    ! 5. Check overlap between Cvirtsub and occ orbitals for other subsystems
!!$    call mem_alloc(tmp,nvirt,noccfull)
!!$    call dec_diff_basis_transform1(nbasis,nvirt,noccfull,&
!!$         & Cvirtsub,MyMolecule%Co,MyMolecule%overlap,tmp)
!!$    do j=1,noccfull
!!$
!!$       ! Subsystem to which orbital "j" is assigned
!!$       atomj = OccOrbitals(j)%centralatom
!!$       subj = MyMolecule%SubSystemIndex(atomj)
!!$
!!$       ! Only consider orbital "j" if assigned to another subsystem
!!$       if(subj/=this) then
!!$          do i=1,nvirt
!!$             if(abs(tmp(i,j))>thr) then
!!$                print *, 'Check5: i,j,value,thr',i,j,tmp(i,j),thr
!!$                call lsquit('subsystem_orbitals_sanity_check5: Orbitals are not orthogonal!',-1)
!!$             end if
!!$          end do
!!$       end if
!!$
!!$    end do
!!$    call mem_dealloc(tmp)
!!$
!!$
!!$
!!$    call mem_dealloc(Coccsub)
!!$    call mem_dealloc(Cvirtsub)
!!$
!!$  end subroutine subsystem_orbitals_sanity_check_old
!!$
!!$
!!$
!!$
!!$
!!$  !> Print energy summary for interaction energy calculation
!!$  subroutine local_interaction_energy_print(nsub,EHFsub,Ecorrsub,EHFfull,Ecorrfull)
!!$    implicit none
!!$    !> Number of subsystems
!!$    integer,intent(in) :: nsub
!!$    !> HF energies for subsystems
!!$    real(realk),intent(in) :: EHFsub(nsub)
!!$    !> Correlation energies for subsystems
!!$    real(realk),intent(in) :: Ecorrsub(nsub)
!!$    !> HF and correlation energy for total system
!!$    real(realk),intent(in) :: EHFfull,Ecorrfull
!!$    real(realk) :: EHFint,Ecorrint
!!$    integer :: i
!!$
!!$    write(DECinfo%output,'(1X,a)') ''
!!$    write(DECinfo%output,'(1X,a)') ''
!!$    write(DECinfo%output,'(1X,a)') '***************************************************************'
!!$    write(DECinfo%output,'(1X,a)') '            Local interaction energies - summary               '
!!$    write(DECinfo%output,'(1X,a)') '***************************************************************'
!!$    write(DECinfo%output,'(1X,a,g22.12)') 'Full system     --- HF   energy : ', EHFfull
!!$    write(DECinfo%output,'(1X,a,g22.12)') 'Full system     --- corr energy : ', Ecorrfull
!!$    write(DECinfo%output,'(1X,a,g22.12)') 'Full system     --- tot  energy : ', EHFfull+Ecorrfull
!!$    EHFint = EHFfull
!!$    Ecorrint = Ecorrfull
!!$    do i=1,nsub
!!$       write(DECinfo%output,'(1X,a)') ''
!!$       write(DECinfo%output,'(1X,a,i5,a,g21.12)') 'Subsystem: ', i, &
!!$            & ' --- HF   energy: ', EHFsub(i)
!!$       write(DECinfo%output,'(1X,a,i5,a,g21.12)') 'Subsystem: ', i, &
!!$            & ' --- corr energy: ', Ecorrsub(i)
!!$       write(DECinfo%output,'(1X,a,i5,a,g21.12)') 'Subsystem: ', i, &
!!$            & ' --- tot  energy: ', EHFsub(i)+Ecorrsub(i)
!!$       EHFint = EHFint - EHFsub(i)
!!$       Ecorrint = Ecorrint - Ecorrsub(i)
!!$    end do
!!$    write(DECinfo%output,'(1X,a)') ''
!!$    write(DECinfo%output,'(1X,a)') '---------------------------------------------------------------'
!!$    write(DECinfo%output,'(1X,a,g22.12)') 'HF Interaction energy      = ', EHFint
!!$    write(DECinfo%output,'(1X,a,g22.12)') 'Corr Interaction energy    = ', Ecorrint
!!$    write(DECinfo%output,'(1X,a,g22.12)') 'Total Interaction energy   = ', EHFint+Ecorrint
!!$    write(DECinfo%output,'(1X,a)') '---------------------------------------------------------------'
!!$    write(DECinfo%output,'(1X,a)') ''
!!$    write(DECinfo%output,'(1X,a)') ''
!!$
!!$  end subroutine local_interaction_energy_print
!!$
!!$
!!$
!!$
!!$
!!$  !> Calculate MP2 correlation energy for subsystem
!!$  subroutine subsystem_correlation_energy(this,MyMolecule,OccOrbitals,Cocc_mat,&
!!$       & Cvirt_mat,F_mat,lssub,Ecorr)
!!$    implicit none
!!$
!!$    !> Subsystem label
!!$    integer,intent(in) :: this
!!$    !> Molecule info
!!$    type(fullmolecule),intent(in) :: MyMolecule
!!$    !> Occ orbitals in DEC format
!!$    type(decorbital),intent(in) :: OccOrbitals(MyMolecule%nocc)
!!$    !> Occupied and virtual MO coefficients for subsystem
!!$    type(matrix),intent(in) :: Cocc_mat, Cvirt_mat
!!$    !> Fock matrix in AO basis for subsystem
!!$    type(matrix),intent(in) :: F_mat
!!$    !> LSitem for subsystem
!!$    type(lsitem), intent(inout) :: lssub
!!$    !> Subsystem correlation energy
!!$    real(realk),intent(inout) :: Ecorr
!!$    real(realk),pointer :: Cocc(:,:), Cvirt(:,:), tmp(:,:),ppfock(:,:),qqfock(:,:),F(:,:)
!!$    integer :: nbasis,nocc,nvirt,ncore,offset,noccall,i,a,j,b
!!$    logical :: local
!!$    type(array2) :: t1
!!$    type(array4) :: t2,VOVO,t2_tmp
!!$
!!$    local=.true.
!!$#ifdef VAR_MPI
!!$    if(infpar%lg_nodtot>1)local=.false.
!!$#endif
!!$
!!$    nbasis = Cocc_mat%nrow  ! basis functions 
!!$    nvirt = Cvirt_mat%ncol  ! virtual orbitals
!!$    ncore = get_number_of_core_orbitals_in_subsystem(this,MyMolecule,OccOrbitals)  ! core orbitals
!!$    noccall = Cocc_mat%ncol   ! Core + valence orbitals for subsystem
!!$    if(DECinfo%frozencore) then
!!$       ! nocc is the number of valence orbitals for subsystem
!!$       nocc = noccall - ncore
!!$       offset=ncore
!!$    else
!!$       ! nocc is the number of core+valence orbitals for subsystem
!!$       nocc = noccall
!!$       offset=0
!!$    end if
!!$
!!$
!!$    ! Virtual MO coefficients in Fortran format (simple copy)
!!$    call mem_alloc(Cvirt,nbasis,nvirt)
!!$    call mat_to_full(Cvirt_mat, 1.0_realk, Cvirt)    
!!$
!!$    ! Occupied MO coefficients in Fortran format
!!$    ! ******************************************
!!$    ! (caution to take only valence orbitals for frozen core approx)
!!$    call mem_alloc(tmp,nbasis,noccall)
!!$    call mat_to_full(Cocc_mat, 1.0_realk, tmp)    ! core+valence orbitals
!!$    call mem_alloc(Cocc,nbasis,nocc)
!!$    do i=1,nocc
!!$       Cocc(:,i)  = tmp(:,i+offset)
!!$    end do
!!$    call mem_dealloc(tmp)
!!$
!!$    ! AO Fock matrix
!!$    call mem_alloc(F,nbasis,nbasis)
!!$    call mat_to_full(F_mat, 1.0_realk, F)
!!$
!!$    ! Occ-occ and virt-virt blocks of Fock matrix
!!$    call mem_alloc(ppfock,nocc,nocc)
!!$    call mem_alloc(qqfock,nvirt,nvirt)
!!$    call dec_simple_basis_transform1(nbasis,nocc,Cocc,F,ppfock)
!!$    call dec_simple_basis_transform1(nbasis,nvirt,Cvirt,F,qqfock)
!!$
!!$    if(DECinfo%ccmodel==MODEL_MP2) then
!!$       call get_VOVO_integrals(lssub,nbasis,nocc,nvirt,&
!!$            & Cvirt,Cocc,VOVO)
!!$       call mp2_solver(nocc,nvirt,ppfock,qqfock,VOVO,t2)
!!$    elseif(DECinfo%ccmodel==MODEL_CCSD) then
!!$
!!$       call ccsolver_par(DECinfo%ccmodel,Cocc,Cvirt,F,nbasis,nocc,nvirt,&
!!$            & lssub,DECinfo%PL,ppfock,qqfock,Ecorr,&
!!$            & t1,t2_tmp,VOVO,.false.,local,DECinfo%use_pnos)
!!$       call get_combined_SingleDouble_amplitudes(t1,t2_tmp,t2)
!!$       call array4_free(t2_tmp)
!!$       call array2_free(t1)
!!$    else
!!$       print *, 'Model ', DECinfo%ccmodel
!!$       call lsquit('subsystem_correlation_energy: Not implemented for model',-1)
!!$    end if
!!$
!!$
!!$    Ecorr=0.0_realk
!!$
!!$    do j=1,nocc
!!$       do b=1,nvirt
!!$          do i=1,nocc
!!$             do a=1,nvirt
!!$                Ecorr = Ecorr + &
!!$                     & t2%val(a,i,b,j)*( 2.0_realk*VOVO%val(a,i,b,j) - VOVO%val(b,i,a,j)  )
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$
!!$
!!$
!!$    call mem_dealloc(F)
!!$    call mem_dealloc(Cocc)
!!$    call mem_dealloc(Cvirt)
!!$    call mem_dealloc(ppfock)
!!$    call mem_dealloc(qqfock)
!!$    call array4_free(t2)
!!$    call array4_free(VOVO)
!!$
!!$  end subroutine subsystem_correlation_energy
!!$
!!$
!!$
!!$  !> Calculate correlation energy for full system
!!$  subroutine fullsystem_correlation_energy(MyMolecule,OccOrbitals,VirtOrbitals,&
!!$       & lsfull,Ecorr)
!!$    implicit none
!!$
!!$    !> Molecule info
!!$    type(fullmolecule),intent(in) :: MyMolecule
!!$    !> Orbitals in DEC format
!!$    type(decorbital) :: OccOrbitals(MyMolecule%nocc),VirtOrbitals(MyMolecule%nunocc)
!!$    !> LSitem for full system
!!$    type(lsitem), intent(inout) :: lsfull
!!$    !> correlation energy
!!$    real(realk),intent(inout) :: Ecorr
!!$    real(realk),pointer :: Cocc(:,:), Cvirt(:,:), tmp(:,:),ppfock(:,:),qqfock(:,:)
!!$    integer :: nbasis,nocc,nvirt,ncore,offset,noccall,i,a,j,b
!!$    type(array2) :: t1
!!$    type(array4) :: t2,VOVO
!!$    type(decfrag) :: FullFrag
!!$    logical :: occ_list(MyMolecule%nocc), unocc_list(MyMolecule%nunocc)
!!$    type(array) :: t1_arr,t2_arr,VOVO_arr,u_arr
!!$
!!$    nbasis = MyMolecule%nbasis
!!$    nvirt = MyMolecule%nunocc
!!$    ncore = MyMolecule%ncore
!!$    noccall = MyMolecule%nocc
!!$
!!$    if(DECinfo%frozencore) then
!!$       ! nocc is the number of valence orbitals for subsystem
!!$       nocc = noccall - ncore
!!$       offset=ncore
!!$    else
!!$       ! nocc is the number of core+valence orbitals for subsystem
!!$       nocc = noccall
!!$       offset=0
!!$    end if
!!$
!!$    ! Virtual MO coefficients in Fortran format (simple copy)
!!$    call mem_alloc(Cvirt,nbasis,nvirt)
!!$    Cvirt = MyMolecule%Cv
!!$
!!$    ! Occupied MO coefficients in Fortran format
!!$    ! ******************************************
!!$    ! (caution to take only valence orbitals for frozen core approx)
!!$    call mem_alloc(Cocc,nbasis,nocc)
!!$    do i=1,nocc
!!$       Cocc(:,i)  = MyMolecule%Co(:,i+offset)
!!$    end do
!!$
!!$    ! Occ-occ and virt-virt blocks of Fock matrix
!!$    call mem_alloc(ppfock,nocc,nocc)
!!$    call mem_alloc(qqfock,nvirt,nvirt)
!!$    call dec_simple_basis_transform1(nbasis,nocc,Cocc,mymolecule%fock,ppfock)
!!$    call dec_simple_basis_transform1(nbasis,nvirt,Cvirt,mymolecule%fock,qqfock)
!!$
!!$
!!$    WhichModel: if(DECinfo%ccmodel==MODEL_MP2) then
!!$
!!$
!!$       call get_VOVO_integrals(lsfull,nbasis,nocc,nvirt,&
!!$            & Cvirt,Cocc,VOVO)
!!$       call mp2_solver(nocc,nvirt,ppfock,qqfock,VOVO,t2)
!!$
!!$    elseif(DECinfo%ccmodel==MODEL_CCSD) then
!!$
!!$       ! Init fragment which is the full molecule
!!$       occ_list=.true.
!!$       unocc_list=.true.
!!$       call atomic_fragment_init_orbital_specific(1,MyMolecule%nunocc,&
!!$            & MyMolecule%nocc, unocc_list, &
!!$            & occ_list,OccOrbitals,VirtOrbitals,MyMolecule,lsfull,&
!!$            & FullFrag,.true.,.false.)
!!$
!!$       ! Solve CCSD eq
!!$       call fragment_ccsolver(FullFrag,t1_arr,t2_arr,VOVO_arr)
!!$
!!$       ! Integrals in array4 format
!!$       VOVO = array4_init(VOVO_arr%dims)       
!!$       call array_convert(VOVO_arr,VOVO%val)
!!$       call array_free(VOVO_arr)
!!$
!!$       ! Combined amplitudes in array4 format
!!$       call get_combined_SingleDouble_amplitudes(t1_arr,t2_arr,u_arr)
!!$       call array_free(t1_arr)
!!$       call array_free(t2_arr)
!!$       t2 = array4_init(u_arr%dims)       
!!$       call array_convert(u_arr,t2%val)
!!$       call array_free(u_arr)
!!$
!!$       call atomic_fragment_free(FullFrag)
!!$
!!$
!!$    else
!!$       print *, 'Model ', DECinfo%ccmodel
!!$       call lsquit('fullsystem_correlation_energy: Not implemented for model',-1)
!!$    end if WhichModel
!!$
!!$    Ecorr=0.0_realk
!!$
!!$    do j=1,nocc
!!$       do b=1,nvirt
!!$          do i=1,nocc
!!$             do a=1,nvirt
!!$                Ecorr = Ecorr + &
!!$                     & t2%val(a,i,b,j)*( 2.0_realk*VOVO%val(a,i,b,j) - VOVO%val(b,i,a,j)  )
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    call interactionenergy_test(VOVO,t2,MyMolecule,OccOrbitals,VirtOrbitals)
!!$
!!$    call mem_dealloc(Cocc)
!!$    call mem_dealloc(Cvirt)
!!$    call mem_dealloc(ppfock)
!!$    call mem_dealloc(qqfock)
!!$
!!$    call array4_free(t2)
!!$    call array4_free(VOVO)
!!$
!!$  end subroutine fullsystem_correlation_energy
!!$
!!$
!!$  !> \brief Calculate MP2 density matrix for subsystem
!!$  subroutine get_mp2_density_for_subsystem(this,mylsitem,Cocc,Cvirt,F,CoccMP2,CvirtMP2)
!!$    implicit none
!!$
!!$    !> Which subsystem to consider
!!$    integer,intent(in) :: this
!!$    !> Integral info
!!$    type(lsitem), intent(inout) :: mylsitem
!!$    !> Occupied and virtual MOs
!!$    type(matrix), intent(in) :: Cocc
!!$    type(matrix), intent(in) :: Cvirt
!!$    !> AO Fock matrix
!!$    type(matrix),intent(in) :: F
!!$    !> First occ/virt MP2 orbitals
!!$    type(matrix),intent(inout) :: CoccMP2,CvirtMP2
!!$    type(matrix) :: DHF,S,C,rho,DMP2ao,CMP2,U,Sunit,DMP2mo
!!$    integer :: nocc,nvirt,nbasis,i,j
!!$    type(fullmolecule) :: SubMolecule
!!$    type(fullmp2grad) :: fullgrad
!!$    real(realk) :: EHF,Ecorr,Eerr
!!$    real(realk),pointer :: molgrad(:,:),eival(:)
!!$    
!!$
!!$    DECinfo%InclFullMolecule=.true.
!!$    DECinfo%simulate_full=.true.
!!$    DECinfo%first_order=.true.
!!$    DECinfo%density=.true.
!!$
!!$    nbasis = Cocc%nrow
!!$    nocc = Cocc%ncol
!!$    nvirt = Cvirt%ncol
!!$
!!$    ! HF density matrix for subsystem
!!$    call mat_init(DHF,nbasis,nbasis)
!!$    call get_density_from_occ_orbitals_mat(Cocc,DHF)
!!$
!!$    ! Overlap matrix
!!$    call mat_init(S,nbasis,nbasis)
!!$    call II_get_overlap(DECinfo%output,DECinfo%output,mylsitem%setting,S)
!!$
!!$    ! MOs collected
!!$    call mat_init(C,nbasis,nbasis)
!!$    call collect_MO_coeff_in_one_matrix(Cocc,Cvirt,C)
!!$
!!$    ! Molecule for subsystem
!!$    call molecule_init_from_inputs(SubMolecule,mylsitem,F,S,C,DHF)
!!$
!!$    ! Calculate MP2 density structure
!!$    call mem_alloc(molgrad,3,SubMolecule%natoms)
!!$    call DEC_wrapper(SubMolecule,mylsitem,DHF,EHF,Ecorr,molgrad,Eerr,fullgrad=fullgrad)
!!$
!!$    ! Correlation density matrix
!!$    call mat_init(rho,nbasis,nbasis)
!!$    call mat_set_from_full(fullgrad%rho,1.0_realk,rho)
!!$
!!$    ! DMP2 = 2*DHF + rho  (in AO basis)
!!$    call mat_init(DMP2ao,nbasis,nbasis)
!!$    call mat_add(2.0_realk, DHF, 1.0_realk, rho, DMP2ao)
!!$
!!$    ! Express DMP2 in MO basis
!!$    call mat_init(DMP2mo,nbasis,nbasis)
!!$    call util_AO_to_MO_2(S,C,DMP2ao,DMP2mo,.false.)
!!$
!!$    ! Diagonalize DMP2
!!$    call mat_init(U,nbasis,nbasis)
!!$    call mat_init(Sunit,nbasis,nbasis)
!!$    call mat_identity(Sunit)
!!$    call mem_alloc(eival,nbasis)
!!$    call mat_diag_f(DMP2mo,Sunit,eival,U)
!!$    do i=1,nbasis
!!$       write(DECinfo%output,*) 'Subsystem, eival', this, i,eival(i)
!!$    end do
!!$    ! call mat_print_all('U',U,DECinfo%output)
!!$
!!$    ! CMP2 = C U
!!$    call mat_init(CMP2,nbasis,nbasis)
!!$    call mat_mul(C, U, 'N', 'N', 1.0_realk, 0.0_realk, CMP2)
!!$
!!$    ! Save MP2 coefficients in two matrices
!!$    ! ("virtual ones" come BEFORE "occupied ones" in CMP2 matrix)
!!$    call partition_MO_coeff_into_two_matrices2(CMP2,CoccMP2,CvirtMP2)
!!$
!!$
!!$    call free_fullmp2grad(fullgrad)
!!$    call mat_free(DHF)
!!$    call mat_free(S)
!!$    call mat_free(C)
!!$    call mat_free(rho)
!!$    call mat_free(DMP2mo)
!!$    call mat_free(DMP2ao)
!!$    call mat_free(CMP2)
!!$    call mat_free(U)
!!$    call mat_free(Sunit)
!!$    call molecule_finalize(SubMolecule)
!!$    call mem_dealloc(molgrad)
!!$    call mem_dealloc(eival)
!!$
!!$    DECinfo%simulate_full=.false.
!!$    DECinfo%InclFullMolecule=.false.
!!$
!!$  end subroutine get_mp2_density_for_subsystem
!!$
!!$
!!$

end module snoop_main_module
