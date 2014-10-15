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
  use ccintegrals
  use full_molecule
  use full

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
    real(realk),pointer :: EHFsnoop(:),Ecorrsnoop(:),EHFiso(:) 
    real(realk) :: EHFfull,Ecorrfull
    integer :: nsub,nbasis,nvirtfull,noccfull,i,this,noccsnoop,nvirtsnoop,nbasissnoop,nelsnoop
    type(matrix),pointer :: Cocciso(:),Cvirtiso(:),Coccsnoop(:),Cvirtsnoop(:)
    type(matrix) :: FAOsnoop, FAOiso
    type(lsitem) :: lssnoop,lsiso
    integer :: nocciso,nvirtiso,nbasisiso,neliso


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
    IsoLOOP: do i=1,nsub

       this = i

       ! LSitem for SNOOP subsystem
       call build_subsystem_lsitem_no_ghost(this,MyMolecule,lsfull,lsiso)
       
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
       SNOOPLOOP: do i=1,nsub

          this = i

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

          call mat_init(Coccsnoop(i),nbasis,noccsnoop)

          call get_orthogonal_basis_for_subsystem(this,nsub,&
               & MyMolecule,Cocciso,Cvirtiso,Coccsnoop(i),Cvirtsnoop(i))

          ! Sanity check for initial orbitals
          call subsystem_orbitals_sanity_check(Coccsnoop(i),&
               & Cvirtsnoop(i),MyMolecule)

          ! SCF optimization for subsystem "this"
          call mat_init(FAOsnoop,nbasis,nbasis)
          call solve_subsystem_scf_rh(lssnoop,Coccsnoop(i),&
               & Cvirtsnoop(i),FAOsnoop,EHFsnoop(this))

          ! Sanity check for optimized orbitals
          call subsystem_orbitals_sanity_check(Coccsnoop(i),&
               & Cvirtsnoop(i),MyMolecule)

          ! Correlation energy for subsystem
          if(.not. DECinfo%SNOOPjustHF) then
!             call subsystem_correlation_energy(this,MyMolecule,&
!                  & OccOrbitals,Coccsnoop(i),Cvirtsnoop(i),&
!                  & FAOsnoop,lssnoop,Ecorrsnoop(this))
          end if

          print '(1X,a,i5,a,3g20.12)', 'Subsystem: ', this, &
               & '  *** HF/corr/HFdiff energy: ', EHFsnoop(this), Ecorrsnoop(this),EHFsnoop(this)-EHFiso(this)

          write(DECinfo%output,'(1X,a,i5,a,3g20.12)') 'Subsystem: ', this, &
               & '  *** HF/corr/HFdiff energy: ', EHFsnoop(this), Ecorrsnoop(this),EHFsnoop(this)-EHFiso(this)


          ! Free stuff for subsystem
          call mat_free(FAOsnoop)
          call ls_free(lssnoop)

       End do SNOOPLOOP

    ! Print interaction energy summary
    call SNOOP_interaction_energy_print(nsub,EHFsnoop,Ecorrsnoop,EHFfull,Ecorrfull)


    call mem_dealloc(EHFsnoop)
    call mem_dealloc(EHFiso)
    call mem_dealloc(Ecorrsnoop)
    do i=1,nsub
       call mat_free(Cocciso(i))
       call mat_free(Cvirtiso(i))
       call mat_free(Coccsnoop(i))
       call mat_free(Cvirtsnoop(i))
    end do
    call mem_dealloc(Cocciso)
    call mem_dealloc(Cvirtiso)
    call mem_dealloc(Coccsnoop)
    call mem_dealloc(Cvirtsnoop)

  end subroutine snoop_driver




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
    call mem_alloc(B,DECinfo%kappaMaxIter,DECinfo%kappaMaxIter)
    call mem_alloc(c,DECinfo%kappaMaxIter)

    ! Fock matrices in AO basis and residuals for each iteration
    call mem_alloc(FAO,DECinfo%kappaMaxIter)
    call mem_alloc(residual,DECinfo%kappaMaxIter)

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
    convthr = DECinfo%kappaTHR      ! The convergence threshold
    converged=.false.
    prev_norm = huge(1.0_realk)

    write(DECinfo%output,*)
    write(DECinfo%output,*) 'Starting subsystem SCF iterations'
    write(DECinfo%output,*) '---------------------------------'
    write(DECinfo%output,'(1X,a)')  '###  Iteration     Residual norm            HF energy'
    write(6,'(1X,a)')  '###  Iteration     Residual norm            HF energy'


    ! Start Roothaan-Hall iterations
    ! ******************************

    RHiteration : do iter=1,DECinfo%kappaMaxIter

       ! Current iteration
       last_iter = iter

       ! remove old vectors
       RemoveOldVectors : if(iter > DECinfo%kappaMaxDIIS) then
          call mat_free(residual(iter-DECinfo%kappaMaxDIIS))
          call mat_free(FAO(iter-DECinfo%kappaMaxDIIS))
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
       if( (iter == DECinfo%kappaMaxIter) .or. converged) then
          call mat_assign(FAOsub,FAO(iter))
          exit RHiteration
       end if



       ! DIIS ACCELERATION
       ! *****************

       DIISacc: if(dodiis) then

          ! calculate diis matrix
          B=0.0E0_realk; c=0.0E0_realk
          do i=iter,max(iter-DECinfo%kappaMaxDIIS+1,1),-1
             do j=iter,i,-1
                B(i,j) = mat_dotproduct(residual(i),residual(j))
                B(j,i) = B(i,j)
             end do
          end do

          ! solve crop/diis equation
          call CalculateDIIScoefficients(DECinfo%kappaMaxDIIS,DECinfo%kappaMaxIter,iter,B,c, &
               DECinfo%kappa_driver_debug)

          ! mixing to get optimal Fock matrix
          call mat_zero(FAOopt)
          do i=iter,max(iter-DECinfo%kappaMaxDIIS+1,1),-1
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
    do i=last_iter,max(last_iter-DECinfo%kappaMaxDIIS+1,1),-1
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
       & Cvirtsub_mat,MyMolecule)
    implicit none

    !> Occ and virt MO coefficients 
    type(matrix),intent(in) :: Coccsub_mat, Cvirtsub_mat
    !> Full molecule info
    type(fullmolecule),intent(in) :: MyMolecule
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
         & Coccsub,Cvirtsub,MyMolecule%overlap,tmp)
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
    call dec_simple_basis_transform1(nbasis,nocc,Coccsub,MyMolecule%overlap,tmp)
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
    call dec_simple_basis_transform1(nbasis,nvirt,Cvirtsub,MyMolecule%overlap,tmp)
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
  subroutine SNOOP_interaction_energy_print(nsub,EHFsub,Ecorrsub,EHFfull,Ecorrfull)
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
    write(DECinfo%output,'(1X,a)') '            Local interaction energies - summary               '
    write(DECinfo%output,'(1X,a)') '***************************************************************'
    write(DECinfo%output,'(1X,a,g22.12)') 'Full system     --- HF   energy : ', EHFfull
    write(DECinfo%output,'(1X,a,g22.12)') 'Full system     --- corr energy : ', Ecorrfull
    write(DECinfo%output,'(1X,a,g22.12)') 'Full system     --- tot  energy : ', EHFfull+Ecorrfull
    EHFint = EHFfull
    Ecorrint = Ecorrfull
    do i=1,nsub
       write(DECinfo%output,'(1X,a)') ''
       write(DECinfo%output,'(1X,a,i5,a,g21.12)') 'Subsystem: ', i, &
            & ' --- HF   energy: ', EHFsub(i)
       write(DECinfo%output,'(1X,a,i5,a,g21.12)') 'Subsystem: ', i, &
            & ' --- corr energy: ', Ecorrsub(i)
       write(DECinfo%output,'(1X,a,i5,a,g21.12)') 'Subsystem: ', i, &
            & ' --- tot  energy: ', EHFsub(i)+Ecorrsub(i)
       EHFint = EHFint - EHFsub(i)
       Ecorrint = Ecorrint - Ecorrsub(i)
    end do
    write(DECinfo%output,'(1X,a)') ''
    write(DECinfo%output,'(1X,a)') '---------------------------------------------------------------'
    write(DECinfo%output,'(1X,a,g22.12)') 'HF Interaction energy      = ', EHFint
    write(DECinfo%output,'(1X,a,g22.12)') 'Corr Interaction energy    = ', Ecorrint
    write(DECinfo%output,'(1X,a,g22.12)') 'Total Interaction energy   = ', EHFint+Ecorrint
    write(DECinfo%output,'(1X,a)') '---------------------------------------------------------------'
    write(DECinfo%output,'(1X,a)') ''
    write(DECinfo%output,'(1X,a)') ''

  end subroutine SNOOP_interaction_energy_print





  !> Calculate MP2 correlation energy for subsystem
  subroutine subsystem_correlation_energy(Cocc_mat,&
       & Cvirt_mat,F_mat,lssub,Ecorr)
    implicit none

    !> Occupied and virtual MO coefficients for subsystem
    type(matrix),intent(in) :: Cocc_mat, Cvirt_mat
    !> Fock matrix in AO basis for subsystem
    type(matrix),intent(in) :: F_mat
    !> LSitem for subsystem
    type(lsitem), intent(inout) :: lssub
    !> Subsystem correlation energy
    real(realk),intent(inout) :: Ecorr


  end subroutine subsystem_correlation_energy


end module snoop_main_module
