!> @file
!> Contains all informations about molecule and HF calculation
!> \author Marcin Ziolkowski

!> All information about a molecule and HF calculations.
module full_molecule


  ! Outside DEC directory
  use TYPEDEF!,only: count_ncore
  use memory_handling
  use precision
  use lstiming!, only: lstimer
  use typedeftype!,only: lsitem
  use typedef!,only: lsitem
  use files!,only:lsopen,lsclose
  use tensor_interface_module
  use matrix_operations_pdmm, only: BLOCK_SIZE_PDM
  use matrix_module!, only:matrix
  use matrix_util!,only:mat_init,mat_zero,mat_daxpy,mat_free,mat_to_full,&
!       & mat_diag_f, mat_write_to_disk
  use matrix_operations!, only: mat_set_from_full,mat_free,mat_init,mat_to_full
  use dec_typedef_module
  use IntegralInterfaceMod
  use II_XC_interfaceModule

  ! CABS
  use CABS_operations
#ifdef MOD_UNRELEASED
  ! F12 MO-matrices
  use f12_routines_module!,only: get_F12_mixed_MO_Matrices, MO_transform_AOMatrix
#endif
  ! DEC DEPENDENCIES (within deccc directory) 
  ! *****************************************
  use dec_fragment_utils
  use array2_simple_operations

  integer, save :: mol_block_size = -14938343

contains

  !> \brief Initialize informations about full molecule by reading HF info from file.
  !> NOTE: If this routine is modified, then molecule_init_from_inputs must be modified
  !> in the same manner!
  !> \author Marcin Ziolkowski
  !> \param molecule Full molecule info
  !> \param mylsitem Integral program input
  subroutine molecule_init_from_files(molecule,mylsitem,D)

    implicit none
    type(fullmolecule), intent(inout) :: molecule
    type(lsitem), intent(inout) :: mylsitem
    !> Density Matrix 
    type(matrix), optional, intent(in) :: D  ! Needed for creating the hJir MO-matrix
    real(realk) :: tcpu, twall
    
    call LSTIMER('START',tcpu,twall,DECinfo%output)

    ! Init basic info (molecular dimensions etc.)
    call molecule_init_basics(molecule,mylsitem)       
    
    ! Skip read-in of info for molecule if requested (only for testing)
    if(DECinfo%SkipReadIn) then
       write(DECinfo%output,*) 'WARNING: I do NOT read in the molecular info files &
            & as requested in the input!'
       return
    end if

    ! Get Fock, overlap, and MO coefficient matrices.
    call molecule_get_reference_state(molecule,mylsitem)

    call molecule_mo_fock(molecule)

    if(DECinfo%use_canonical) then ! overwrite local orbitals and use canonical orbitals
       call dec_get_canonical_orbitals(molecule,mylsitem)
    end if

    call molecule_get_carmom(molecule,mylsitem)

    !> Interatomic distances in atomic units
    call mem_alloc(molecule%DistanceTable,molecule%nfrags,molecule%nfrags)
    call GetDistances(molecule,mylsitem,DECinfo%output) 

    call mem_alloc(molecule%PhantomAtom,molecule%nAtoms)
    call getPhantomAtoms(mylsitem,molecule%PhantomAtom,molecule%nAtoms)

    if(DECinfo%F12) then 
#ifdef MOD_UNRELEASED
       !> Sanity check 
       if(.NOT. present(D)) then
          call lsquit("ERROR: (molecule_init_from_files) : Density needs to be present for F12 calc",-1)
       end if
       IF(DECinfo%full_molecular_cc)THEN
          call dec_get_CABS_orbitals(molecule,mylsitem)
          call dec_get_RI_orbitals(molecule,mylsitem)
       ELSE
          !> F12 Fock matrices in MO basis
          call molecule_mo_f12(molecule,mylsitem,D)
       ENDIF
#endif
    end if

    
    ! Do not store AO Fock matrix if requested
    if(DECinfo%noaofock) then
       write(DECinfo%output,*) 'Warning: Not storing AO Fock matrix to minimize memory usage - use at own risk!'
       call tensor_free(molecule%fock)
    end if

    call LSTIMER('DEC: MOL INIT',tcpu,twall,DECinfo%output)

  end subroutine molecule_init_from_files


  !> \brief Initialize informations about full molecule
  !> using input Fock and MO matrices.
  !> NOTE: If this routine is modified, then molecule_init_from_files must be modified
  !> in the same manner!
  !> \author Kasper Kristensen
  !> \date November 2011
  subroutine molecule_init_from_inputs(molecule,mylsitem,F,C,D)

    implicit none
    !> Full molecule structure to be initialized
    type(fullmolecule), intent(inout) :: molecule
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> Fock matrix
    type(matrix),intent(in) :: F
    !> MO coefficients
    type(matrix),intent(in) :: C
    !> Density Matrix 
    type(matrix),intent(in) :: D  ! Needed for creating the hJir MO-matrix
    real(realk) :: tcpu, twall
    integer :: nMO
    
    call LSTIMER('START',tcpu,twall,DECinfo%output)

    ! Number of MOs (keep it general such that nMO can be different from nbasis for subsystems)
    nMO = C%ncol

     ! Init basic info (molecular dimensions etc.)
    call molecule_init_basics(molecule,mylsitem,nMO=nMO)

    ! Copy Fock and MO overlap matrices to molecule structure
    call molecule_copy_FC_matrices(molecule,F,C)

    ! Get absolute overlap matrix for space selection
    call molecule_init_abs_overlap(molecule,mylsitem)

    ! Fock matrix in MO basis
    call molecule_mo_fock(molecule)

 
    if(DECinfo%use_canonical) then
       ! overwrite local orbitals and use canonical orbitals
       call dec_get_canonical_orbitals(molecule,mylsitem)
    end if
     
    call molecule_get_carmom(molecule,mylsitem)

    !> Interatomic distances in atomic units
    call mem_alloc(molecule%DistanceTable,molecule%nfrags,molecule%nfrags)
    call GetDistances(molecule,mylsitem,DECinfo%output) 

    call mem_alloc(molecule%PhantomAtom,molecule%nAtoms)
    call getPhantomAtoms(mylsitem,molecule%PhantomAtom,molecule%nAtoms)

    if(DECinfo%F12) then ! overwrite local orbitals and use CABS orbitals
#ifdef MOD_UNRELEASED
       IF(DECinfo%full_molecular_cc)THEN
          call dec_get_CABS_orbitals(molecule,mylsitem)
          call dec_get_RI_orbitals(molecule,mylsitem)
       ELSE
          !> F12 Fock matrices in MO basis
          call molecule_mo_f12(molecule,mylsitem,D)
       ENDIF
#endif
    end if

    ! Do not store AO Fock matrix if requested
    if(DECinfo%noaofock) then
       write(DECinfo%output,*) 'Warning: Not storing AO Fock matrix to minimize memory usage - use at own risk!'
       call tensor_free(molecule%fock)
    end if
    
    call LSTIMER('DEC: MOL INIT',tcpu,twall,DECinfo%output)

  end subroutine molecule_init_from_inputs



  !> \brief Initialize basic informations about full molecule (number of basis functions,
  !> number of occipied orbitals etc.).
  !> \author Marcin Ziolkowski/Kasper Kristensen
  !> \param molecule Full molecule info
  !> \param mylsitem Integral program input
  subroutine molecule_init_basics(molecule,mylsitem,nMO)

    implicit none
    type(fullmolecule), intent(inout) :: molecule
    type(lsitem), intent(inout) :: mylsitem
    !> Number of MOs (need to be input only if it is different from nbasis which 
    !> can be the case for subsystems)
    integer,intent(in),optional :: nMO
    real(realk) :: memory_use, tcpu, twall

    call LSTIMER('START',tcpu,twall,DECinfo%output)

    molecule%EF12singles = 0.0_realk
    molecule%Edisp = 0.0_realk
    molecule%Ect = 0.0_realk
    molecule%Esub = 0.0_realk
    molecule%natoms = get_num_atoms(mylsitem)
    molecule%nelectrons = get_num_electrons(mylsitem)
    molecule%nbasis = get_num_basis_functions(mylsitem)
    molecule%nauxbasis = get_num_aux_basis_functions(mylsitem)

    ! Number of MOs can be different from nbasis if specified by input
    if(present(nMO)) then
       molecule%nMO=nMO
    else
       molecule%nMO = molecule%nbasis
    end if

    molecule%nocc = molecule%nelectrons/2
    molecule%nvirt = molecule%nMO - molecule%nocc
    molecule%ncore = count_ncore(mylsitem)
    molecule%nval = molecule%nocc - molecule%ncore
    molecule%nCabsAO = 0
    molecule%nCabsMO = 0

    ! Number of possible fragments:
    ! natoms for atom-based approach
    ! nocc for orbital-based approach
    if(DECinfo%DECCO) then
       molecule%nfrags=molecule%nocc
    else
       molecule%nfrags=molecule%natoms
    end if

    ! Which basis functions are on which atoms?
    call molecule_get_atomic_sizes(molecule,mylsitem)

    ! Memory use for full molecule
    call calculate_fullmolecule_memory(molecule,memory_use)
    DECinfo%fullmolecule_memory = memory_use

    !> SubSystem index
    call mem_alloc(molecule%SubSystemIndex,molecule%natoms)
    call GetSubSystemIndex(molecule%SubSystemIndex,molecule%natoms,mylsitem,DECinfo%output) 

    !> Which model to use for different pair calculations?
    !> At this initialization step - use the input CC model for all pairs
    call mem_alloc(molecule%ccmodel,molecule%nfrags,molecule%nfrags)
    molecule%ccmodel = DECinfo%ccmodel
    
     !> FOT level to use for each pair calculation
     !>  0: Use input FOT
     !>  n>0: Use AOS information from fragment%REDfrags(n)
    call mem_alloc(molecule%PairFOTlevel,molecule%nfrags,molecule%nfrags)
    molecule%PairFOTlevel=0


    ! Print some info about the molecule
    write(DECinfo%output,*)
    if(molecule%nMO /= molecule%nbasis) then ! subsystem

       write(DECinfo%output,'(/,a)') '-- Subsystem info --'
       write(DECinfo%output,'(/,a,i6)') 'SUB: Overall charge of molecule : ',nint(mylsitem%input%molecule%charge)
       write(DECinfo%output,'(/,a,i6)') 'SUB: Number of electrons        : ',molecule%nelectrons
       write(DECinfo%output,'(a,i6)')   'SUB: Number of atoms            : ',molecule%natoms
       write(DECinfo%output,'(a,i6)')   'SUB: Number of basis func.      : ',molecule%nbasis
       write(DECinfo%output,'(a,i6)')   'SUB: Number of aux. basis func. : ',molecule%nauxbasis
       write(DECinfo%output,'(a,i6)')   'SUB: Number of core orbitals    : ',molecule%ncore
       write(DECinfo%output,'(a,i6)')   'SUB: Number of valence orbitals : ',molecule%nval
       write(DECinfo%output,'(a,i6)')   'SUB: Number of occ. orbitals    : ',molecule%nocc
       write(DECinfo%output,'(a,i6)')   'SUB: Number of virt. orbitals   : ',molecule%nvirt

    else      ! full molecule

       write(DECinfo%output,'(/,a)') '-- Full molecular info --'
       write(DECinfo%output,'(/,a,i6)') 'FULL: Overall charge of molecule : ',nint(mylsitem%input%molecule%charge)
       write(DECinfo%output,'(/,a,i6)') 'FULL: Number of electrons        : ',molecule%nelectrons
       write(DECinfo%output,'(a,i6)')   'FULL: Number of atoms            : ',molecule%natoms
       write(DECinfo%output,'(a,i6)')   'FULL: Number of basis func.      : ',molecule%nbasis
       write(DECinfo%output,'(a,i6)')   'FULL: Number of aux. basis func. : ',molecule%nauxbasis
       write(DECinfo%output,'(a,i6)')   'FULL: Number of core orbitals    : ',molecule%ncore
       write(DECinfo%output,'(a,i6)')   'FULL: Number of valence orbitals : ',molecule%nval
       write(DECinfo%output,'(a,i6)')   'FULL: Number of occ. orbitals    : ',molecule%nocc
       write(DECinfo%output,'(a,i6)')   'FULL: Number of virt. orbitals   : ',molecule%nvirt
       write(DECinfo%output,'(a,g9.2)') 'FULL: Local memory use type full : ',memory_use
       write(DECinfo%output,'(a,L)')    'FULL: Distribute matrices        : ',molecule%mem_distributed

    end if
    write(DECinfo%output,*)

  end subroutine molecule_init_basics


  !> \brief Copy Fock and MO matrices in type(matrix) format to molecule structure.
  !> \author Kasper Kristensen
  !> \date November 2011
  subroutine molecule_copy_FC_matrices(molecule,F,C)

    implicit none
    !> Full molecule structure to be initialized
    type(fullmolecule), intent(inout) :: molecule
    !> Fock matrix
    type(matrix),intent(in) :: F
    !> MO coefficients
    type(matrix),intent(in) :: C
    real(realk),pointer :: basis(:,:)
    logical :: loc
    integer :: tdim(2)

    loc = .not.molecule%mem_distributed
    tdim = [mol_block_size,mol_block_size]

    !FIXME: Avoid the mat_to_full and distribute it directly as tensor
    ! Fock matrix
    call tensor_minit(molecule%fock,[F%nrow,F%ncol],2, local=loc, atype="TDPD",tdims=tdim)
    call mat_to_full(F, 1.0_realk, molecule%fock%elm2)
    if(.not.loc)then
       call tensor_mv_dense2tiled(molecule%fock,.true.,dealloc_local=.true.)
    endif


    ! MO coefficient matrix
    call mem_alloc(basis,C%nrow,C%ncol)
    call mat_to_full(C, 1.0_realk, basis)
    call molecule_generate_basis(molecule,basis)
    call mem_dealloc(basis)

  end subroutine molecule_copy_FC_matrices


  !> \brief Copy Fock and MO matrices in molecule strucutre to type(matrix) format.
  !> (included intialization of matrices)
  !> \author Kasper Kristensen
  !> \date December 2012
  subroutine molecule_copyback_FC_matrices(mylsitem,molecule,F,C)

    implicit none
    !> LS integral info
    type(lsitem), intent(inout) :: mylsitem
    !> Full molecule structure
    type(fullmolecule), intent(inout) :: molecule
    !> Fock matrix
    type(matrix),intent(inout) :: F
    !> MO coefficients
    type(matrix),intent(inout) :: C
    integer :: nbasis,i
    real(realk),pointer :: tmp(:,:)

    nbasis = molecule%nbasis
    call mat_init(C,nbasis,nbasis)
    if(DECinfo%noaofock) then
       ! We need to read Fock file first because it 
       ! is not stored in full molecule structure.
       call molecule_get_fock(molecule,mylsitem)
    end if
    if( Molecule%mem_distributed )then
       call tensor_cp_tiled2dense(Molecule%fock,.false.)
    endif

    call mat_init(F,nbasis,nbasis)
    call mat_set_from_full(Molecule%fock%elm2,1.0_realk,F)



    call mem_alloc(tmp,nbasis,nbasis)
    if( Molecule%mem_distributed )then
       call tensor_deallocate_dense(Molecule%fock)
       call tensor_cp_tiled2dense(Molecule%Co,.false.)
    endif

    ! Put occ orbitals into tmp
    do i=1,Molecule%nocc
       tmp(1:nbasis,i) = Molecule%Co%elm2(1:nbasis,i)
    end do

    if( Molecule%mem_distributed )then
       call tensor_deallocate_dense(Molecule%Co)
       call tensor_cp_tiled2dense(Molecule%Cv,.false.)
    endif

    ! Put virt orbitals into tmp
    do i=1,Molecule%nvirt
       tmp(1:nbasis,i+Molecule%nocc) = Molecule%Cv%elm2(1:nbasis,i)
    end do

    if( Molecule%mem_distributed )then
       call tensor_deallocate_dense(Molecule%Cv)
    endif
    
    ! All orbitals into C
    call mat_set_from_full(tmp,1.0_realk,C)
    call mem_dealloc(tmp)


  end subroutine molecule_copyback_FC_matrices



  !> \brief Get number of atoms in the molecule
  !> \param mylsitem Integral program input
  !> \return Number of atoms
  function get_num_atoms(mylsitem) result(natoms)

    implicit none
    type(lsitem), intent(inout) :: mylsitem
    integer :: natoms

    natoms = mylsitem%input%molecule%nAtoms

  end function get_num_atoms

  !> \brief Number of electrons (as a sum of nuclear charges minus overall charge)
  !> \param mylsitem Integral program input
  !> \return Number of electrons
  function get_num_electrons(mylsitem) result(electrons)

    implicit none
    type(lsitem), intent(inout) :: mylsitem
    integer :: i,electrons,natoms,charge

    electrons = 0
    natoms = mylsitem%input%molecule%natoms
    do i=1,natoms
       IF(mylsitem%input%molecule%Atom(i)%Phantom)CYCLE
       electrons = electrons + mylsitem%input%molecule%Atom(i)%Charge
    end do
    charge = nint(mylsitem%input%molecule%charge)
    electrons = electrons - charge

    return
  end function get_num_electrons

  subroutine getPhantomAtoms(mylsitem,PhantomAtom,nAtoms)
    implicit none
    integer,intent(in) :: nAtoms
    logical,intent(inout) :: PhantomAtom(nAtoms)
    type(lsitem), intent(inout) :: mylsitem
    !
    integer :: i
    do i=1,natoms
       PhantomAtom(i) = mylsitem%input%molecule%Atom(i)%Phantom
    end do
  end subroutine getPhantomAtoms

  !> \brief Get number of regular basis functions
  !> \param mylsitem Integral program input
  !> \return Number of basis functions
  function get_num_basis_functions(mylsitem) result(nbas)
    implicit none
    type(lsitem), intent(inout) :: mylsitem
    integer :: nbas
    nbas = 0
    nbas = mylsitem%input%molecule%nbastREG
    return
  end function get_num_basis_functions

  !> \brief Get number of auxiliary basis functions if used
  !> \param mylsitem Integral program input
  !> \return Number of auxiliary basis functions
  function get_num_aux_basis_functions(mylsitem) result(naux)

    implicit none
    type(lsitem), intent(inout) :: mylsitem
    integer :: naux
    naux = 0
    naux = mylsitem%input%molecule%nbastAUX
    return
  end function get_num_aux_basis_functions


  !> \brief Read or construct (if it does not already exist) the fock matrix
  subroutine molecule_get_fock(molecule,mylsitem)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(inout) :: molecule
    !> LSitem
    type(lsitem), intent(inout) :: mylsitem
    integer :: nbasis
    logical :: fock_exist, loc
    integer :: tdim(2)

    loc = .not.molecule%mem_distributed
    tdim = [mol_block_size,mol_block_size]

    ! Init stuff
    nbasis = molecule%nbasis
    call tensor_minit(molecule%fock,[nbasis,nbasis],2, local=loc, atype="TDPD",tdims=tdim)
    inquire(file='fock.restart',exist=fock_exist)

    ! Read or construct fock matrix
    ! *****************************
    if(fock_exist) then

       ! Read Fock matrix from file
       if( .not. loc) print *,"WARNING(molecule_get_fock): reading to dense, this should be MPI I/O"
       write(DECinfo%output,*) 'Reading Fock matrix from file fock.restart...'
       write(DECinfo%output,*)
       call dec_read_mat_from_file('fock.restart',nbasis,nbasis,molecule%fock%elm2)

    else

       call lsquit('molecule_get_fock: Fock matrix does not exist!',-1)

    end if

    if(.not.loc)then
       call tensor_mv_dense2tiled(molecule%fock,.true.,dealloc_local=.true.)
    endif

  end subroutine molecule_get_fock

  !> \brief Set orbitals used in DEC to canonical orbitals (only for testing).
  !> Assumes that Fock matrix is stored in molecule%fock when this subroutine is called.
  !> After this call molecule%Co and molecule%Cv will contain the occupied and virtual MO
  !> coefficients, respectively, while molecule%ppfock and molecule%qqfock will be the occ-occ
  !> and virt-virt blocks of the diagonal canonical MO Fock matrix.
  !> \author Kasper Kristensen
  !> \date September 2011
  subroutine dec_get_canonical_orbitals(molecule,mylsitem)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(inout) :: molecule
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    integer :: nbasis,i,nocc,nvirt
    real(realk), pointer :: eival(:), C(:,:), S(:,:)

    ! KK quick fix: Do NOT call this subroutine when nMO/=nbasis,
    !               this is the case for SNOOP subsystems. 
    !               A proper fix is needed here, on my todo-list... 
    if(molecule%nMO/=molecule%nbasis) then
       write(DECinfo%output,*) 'WARNING: Quitting dec_get_canonical_orbitals because nMO/=nbasis!'
       write(DECinfo%output,*) 'WARNING: Proper solution is required!'
       return
    end if

    if(DECinfo%noaofock) then
       call lsquit('ERROR(dec_get_canonical_orbitals): You cannot use canonical orbitals &
          &in combination with .NOAOFOCK keyword!',-1)
    end if
    if( molecule%mem_distributed )then
       call lsquit('ERROR(dec_get_canonical_orbitals): You cannot use canonical orbitals &
          &when the matrices in the molecule structure are distributed',-1)
    endif

    nbasis = molecule%nbasis
    nocc = molecule%nocc
    nvirt = molecule%nvirt

    ! AO overlap
    call mem_alloc(S,nbasis,nbasis)
    call II_get_mixed_overlap_full(DECinfo%output,DECinfo%output,MyLsitem%SETTING,&
         & S,nbasis,nbasis,AORdefault,AORdefault)

    ! Canonical MO coefficients
    call mem_alloc(C,nbasis,nbasis)

    ! Orbital energies
    call mem_alloc(eival,nbasis)

    ! Diagonalize Fock matrix
    call solve_eigenvalue_problem(nbasis,molecule%fock%elm2,S,eival,C)

    ! Set MO coefficients
    Molecule%Co%elm2 = C(:,1:nocc)   ! occupied 
    Molecule%Cv%elm2 = C(:,nocc+1:nbasis)   ! virtupied

    ! Set Fock matrix in canonical MO basis 
    Molecule%oofock%elm2=0.0_realk
    Molecule%vvfock%elm2=0.0_realk
    do i=1,nocc
       molecule%oofock%elm2(i,i) = eival(i)
    end do
    do i=1,nvirt
       molecule%vvfock%elm2(i,i) = eival(i+nocc)
    end do

    write(DECinfo%output,*) 'Orbital energies:'
    do i=1,nbasis
       write(DECinfo%output,*) i, eival(i)
    end do

    call mem_dealloc(C)
    call mem_dealloc(eival)
    call mem_dealloc(S)

  end subroutine dec_get_canonical_orbitals


  !> \brief Get information about HF solution (fock and orbitals)
  !> \param molecule Full molecule info
  !> \param mylsitem Integral program input
  !> \author Marcin Ziolkowski
  subroutine molecule_get_reference_state(molecule,mylsitem)

    implicit none
    type(fullmolecule), intent(inout) :: molecule
    type(lsitem), intent(inout) :: mylsitem
    integer :: nbasis
    real(realk),pointer :: C(:,:)

    nbasis = molecule%nbasis

    ! Get Fock matrix
    call molecule_get_fock(molecule,mylsitem)

    ! Get orbitals
    call mem_alloc(C,nbasis,nbasis)
    IF(DECinfo%use_canonical)THEN
       call dec_read_mat_from_file('cmo_orbitals.u',nbasis,nbasis,C)       
    ELSE
       call dec_read_mat_from_file('lcm_orbitals.u',nbasis,nbasis,C)
    ENDIF
    call molecule_generate_basis(molecule,C)
    call mem_dealloc(C)

  end subroutine molecule_get_reference_state




  !> \author Patrick Ettenhuber
  !> \date Feb 2015
  !> \calclate the numerical overlap between all occ and virt orbitals for prioritorizing
  subroutine molecule_init_abs_overlap(molecule,mylsitem)
    implicit none
    !> Full molecule info
    type(fullmolecule), intent(inout) :: molecule
    !> LSitem
    type(lsitem), intent(inout) :: mylsitem


    if(DECinfo%use_abs_overlap)then

       if(molecule%mem_distributed)then
          call lsquit("ERROR(molecule_init_abs_overlap): not working with PDM yet",-1)
       endif

       call mem_alloc(molecule%ov_abs_overlap,molecule%nocc,molecule%nvirt)

       molecule%ov_abs_overlap=0.0E0_realk

       call II_get_AbsoluteValue_overlap(DECinfo%output,6,Mylsitem%SETTING,molecule%nbasis,&
          &molecule%nocc,molecule%nvirt,molecule%Co%elm2,molecule%Cv%elm2,molecule%ov_abs_overlap)

    endif

 end subroutine molecule_init_abs_overlap

  !> \brief Read or construct (if it does not already exist) the carmom matrices
  !> \author Thomas Kjaergaard
  !> \date January 2013
  subroutine molecule_get_carmom(molecule,mylsitem)
    implicit none
    !> Full molecule info
    type(fullmolecule), intent(inout) :: molecule
    !> LSitem
    type(lsitem), intent(inout) :: mylsitem
    type(matrix) :: XYZmat(4),Cocc,Cvirt,Xocc,Xvirt
    integer :: nbasis,nocc,nvirt,natoms,nmat,nderiv,XYZ,I
    real(realk) :: CenterX,CenterY,CenterZ


    ! Init stuff
    nbasis = molecule%nbasis
    nocc = molecule%nocc
    nvirt = molecule%nvirt
    natoms = molecule%natoms

!    inquire(file='carmommatrix',exist=carmom_exist)
    ! Read or construct carmom matrices
    ! ********************************
!    if(carmom_exist) then
!       call lsquit('not implemented',-1)
       ! Read overlap matrix from file
!       write(DECinfo%output,*) 'Reading carmom matrices from file carmommatrix...'
!       write(DECinfo%output,*)
!       call dec_read_mat_tensor_from_file('carmommatrix',nbasis,nbasis,molecule%carmom,3)
!    else
       ! Calculate carmom matrix from scratch
       write(DECinfo%output,*) 'Calculating carmom matrix for DEC calculation...'
       call mat_init(XYZmat(1),nbasis,nbasis) !overlap matrix
       call mat_init(XYZmat(2),nbasis,nbasis) !X matrix
       call mat_init(XYZmat(3),nbasis,nbasis) !Y matrix
       call mat_init(XYZmat(4),nbasis,nbasis) !Z matrix
!       call mat_zero(XYZmat)
       nMat = 4
       nDeriv = 1
       CenterX = 0.0E0_realk
       CenterY = 0.0E0_realk
       CenterZ = 0.0E0_realk
       call II_get_carmom(DECinfo%output,DECinfo%output,mylsitem%setting,XYZmat,&
            & nMat,nDeriv,CenterX,CenterY,CenterZ)

       call mat_free(XYZmat(1)) !no need for overlap matrix
       
       ! Set MO orbitals
       ! ****************
       if( Molecule%mem_distributed )then
          call tensor_cp_tiled2dense(Molecule%Co,.false.)
          call tensor_cp_tiled2dense(Molecule%Cv,.false.)
       endif

       call mat_init(Cocc,nbasis,nocc)
       call mat_init(Cvirt,nbasis,nvirt)

       call mat_init(Xocc,nocc,nocc)
       call mat_init(Xvirt,nvirt,nvirt)

       call mat_set_from_full(Molecule%Co%elm2(1:nbasis,1:nocc), 1E0_realk,Cocc)
       call mat_set_from_full(Molecule%Cv%elm2(1:nbasis,1:nvirt), 1E0_realk,Cvirt)

       if( Molecule%mem_distributed )then
          call tensor_deallocate_dense(Molecule%Co)
          call tensor_deallocate_dense(Molecule%Cv)
       endif

       call mem_alloc(molecule%carmomocc,3,nocc)
       call mem_alloc(molecule%carmomvirt,3,nvirt)

       do XYZ=1,3
          call util_AO_to_MO_different_trans(Cocc, XYZmat(XYZ+1), Cocc, Xocc)
          call util_AO_to_MO_different_trans(Cvirt, XYZmat(XYZ+1), Cvirt, Xvirt)
          call mat_free(XYZmat(XYZ+1)) 
          !extract diagonal for X
          do I = 1,nocc
             molecule%carmomocc(XYZ,I) = Xocc%elms(I+(I-1)*nocc)
          enddo
          do I = 1,nvirt
             molecule%carmomvirt(XYZ,I) = Xvirt%elms(I+(I-1)*nvirt)
          enddo
       enddo

       call mat_free(Cocc)
       call mat_free(Cvirt)
       call mat_free(Xocc)
       call mat_free(Xvirt)
       
       call mem_alloc(molecule%AtomCenters,3,nAtoms)
       call getAtomicCenters(mylsitem%setting,molecule%AtomCenters,nAtoms)

  end subroutine molecule_get_carmom


  !> \brief Destroy fullmolecule structure
  !> \param molecule Full molecular info
  subroutine molecule_finalize(molecule,free_pdm)

    implicit none
    type(fullmolecule), intent(inout) :: molecule
    logical, intent(in) :: free_pdm
    logical :: free_tensors

    ! Do not free anything if it has not been initiated (mainly for testing)
    if(DECinfo%SkipReadIn) return

    if(molecule%mem_distributed)then
       free_tensors = free_pdm
    else
       free_tensors = .true.
    endif
    
    ! Delete transformation matrices for general basis
    if(molecule%Co%initialized.and.free_tensors) then
       call tensor_free(molecule%Co)
    end if

    if(molecule%Cv%initialized.and.free_tensors) then
       call tensor_free(molecule%Cv)
    end if

    ! Delete AO fock matrix
    if(molecule%fock%initialized.and.free_tensors) then
       call tensor_free(molecule%fock)
    end if

    ! OOFock
    if(molecule%oofock%initialized.and.free_tensors) then
       call tensor_free(molecule%oofock)
    end if

    ! VVFock
    if(molecule%vvfock%initialized.and.free_tensors) then
       call tensor_free(molecule%vvfock)
    end if

    !Deallocate CABS MO!
!    if(associated(molecule%Ccabs)) then
!       call mem_dealloc(molecule%Ccabs)
!    end if

    !Deallocate CABS RI MO!
!    if(associated(molecule%Cri)) then
!       call mem_dealloc(molecule%Cri)
!    end if

    ! Delete F12-Fock and K and hJir info
    if(associated(molecule%Fij)) then
       call mem_dealloc(molecule%Fij)
    end if

    if(associated(molecule%hJir)) then
       call mem_dealloc(molecule%hJir)
    end if

    if(associated(molecule%Krs)) then
       call mem_dealloc(molecule%Krs)
    end if

    if(associated(molecule%Frs)) then
       call mem_dealloc(molecule%Frs)
    end if

    if(associated(molecule%Fac)) then
       call mem_dealloc(molecule%Fac)
    end if

    if(associated(molecule%Frm)) then
       call mem_dealloc(molecule%Frm)
    end if

    if(associated(molecule%Fcp)) then
       call mem_dealloc(molecule%Fcp)
    end if

!!$    if(associated(molecule%Fcd)) then
!!$       call mem_dealloc(molecule%Fcd)
!!$    end if

    ! Delete atomic info
    if(associated(molecule%atom_size)) then
       call mem_dealloc(molecule%atom_size)
    end if

    if(associated(molecule%atom_start)) then
       call mem_dealloc(molecule%atom_start)
    end if

    if(associated(molecule%atom_end)) then
       call mem_dealloc(molecule%atom_end)
    end if

    if(associated(molecule%bas_start)) then
       call mem_dealloc(molecule%bas_start)
    end if

    if(associated(molecule%bas_end)) then
       call mem_dealloc(molecule%bas_end)
    end if

    if(associated(molecule%atom_cabssize)) then
       call mem_dealloc(molecule%atom_cabssize)
    end if

    if(associated(molecule%atom_cabsstart)) then
       call mem_dealloc(molecule%atom_cabsstart)
    end if

    if(associated(molecule%carmomocc)) then
       call mem_dealloc(molecule%carmomocc)
    end if

    if(associated(molecule%carmomvirt)) then
       call mem_dealloc(molecule%carmomvirt)
    end if

    if(associated(molecule%AtomCenters)) then
       call mem_dealloc(molecule%AtomCenters)
    end if

    if(associated(molecule%PhantomAtom)) then
       call mem_dealloc(molecule%PhantomAtom)
    end if

    if(associated(molecule%SubSystemIndex)) then
       call mem_dealloc(molecule%SubSystemIndex)
    end if

    if(associated(molecule%DistanceTable)) then
       call mem_dealloc(molecule%DistanceTable)
    end if

    if(associated(molecule%ccmodel)) then
       call mem_dealloc(molecule%ccmodel)
    end if

    if(associated(molecule%PairFOTlevel)) then
       call mem_dealloc(molecule%PairFOTlevel)
    end if

    if(associated(molecule%ov_abs_overlap)) then
       call mem_dealloc(molecule%ov_abs_overlap)
    end if

    call free_cabs()

  end subroutine molecule_finalize

  !> \brief Get number of atomic orbitals on atoms, first and last index in AO basis for full molecular matrices
  !> \param molecule Full molecule info
  !> \param mylsitem Integral program input
  subroutine molecule_get_atomic_sizes(molecule,mylsitem)

    implicit none
    type(fullmolecule), intent(inout) :: molecule
    type(lsitem), intent(inout) :: mylsitem
    integer :: natoms,r,i,iset,itype,basis,icharge,nbasis,iOrbitalIndex
    integer :: kmult,nAngmom,ang,norb,j,k,nCont,iOrbitalIndexSave
    logical :: status_info

    natoms = molecule%natoms
    nbasis = molecule%nbasis
    call mem_alloc(molecule%atom_size,natoms)
    molecule%atom_size=0

    r = mylsitem%input%basis%binfo(RegBasParam)%labelindex

    ! loop over atoms
    do i=1,natoms

       if(r == 0) then
          icharge = int(mylsitem%input%molecule%atom(i)%charge)
          itype = mylsitem%input%basis%binfo(RegBasParam)%chargeindex(icharge)
       else
          itype = mylsitem%input%molecule%atom(i)%idtype(1)
       end if

       molecule%atom_size(i) = mylsitem%input%basis%binfo(RegBasParam)%&
            atomtype(itype)%TotNOrb

    end do

    ! get first and last index of an atom in ao matrix
    call mem_alloc(molecule%atom_start,natoms)
    molecule%atom_start = 0
    call mem_alloc(molecule%atom_end,natoms)
    molecule%atom_end = 0

    molecule%atom_start(1) = 1
    molecule%atom_end(1) = molecule%atom_start(1) + molecule%atom_size(1)-1
    basis=1
    do i=1,natoms-1
       basis = basis + molecule%atom_size(i)
       molecule%atom_start(i+1) = basis
       molecule%atom_end(i+1) = molecule%atom_start(i+1) &
            + molecule%atom_size(i+1)-1
    end do

    ! get first and last index of an basis shell in ao matrix
    call mem_alloc(molecule%bas_start,nbasis)
    molecule%bas_start = 0
    call mem_alloc(molecule%bas_end,nbasis)
    molecule%bas_end = 0

    ! loop over atoms
    iOrbitalIndex = 0
    do i=1,natoms
       if(r == 0) then
          icharge = int(mylsitem%input%molecule%atom(i)%charge)
          itype = mylsitem%input%basis%binfo(RegBasParam)%chargeindex(icharge)
       else
          itype = mylsitem%input%molecule%atom(i)%idtype(1)
       end if

       nAngmom = mylsitem%input%basis%binfo(RegBasParam)%ATOMTYPE(itype)%nAngmom
       kmult = 1
       iOrbitalIndexSave = iOrbitalIndex
       nCont = mylsitem%input%basis%binfo(RegBasParam)%ATOMTYPE(itype)%ToTnorb
       do ang = 0,nAngmom-1
          norb = mylsitem%input%basis%binfo(RegBasParam)%ATOMTYPE(itype)%SHELL(ang+1)%norb
          IF(kmult.EQ.1)THEN
             !Include all S orbitals if one S orbital is included  
             do j = 1, norb
                if(DECinfo%AtomicExtent)then
                   molecule%bas_start(iOrbitalIndex+j) = iOrbitalIndexSave+1
                   molecule%bas_end(iOrbitalIndex+j) = iOrbitalIndexSave+nCont
                else
                   molecule%bas_start(iOrbitalIndex+j) = iOrbitalIndex+1
                   molecule%bas_end(iOrbitalIndex+j) = iOrbitalIndex+norb
                endif
             enddo
             iOrbitalIndex = iOrbitalIndex + norb
             kmult = kmult + 2
          ELSE
             do j = 1, norb
                !Include all Orbital Components of a given function (For P include Px,Py,Pz)
                do k=1,kmult
                   if(DECinfo%AtomicExtent)then
                      molecule%bas_start(iOrbitalIndex+k) = iOrbitalIndexSave+1
                      molecule%bas_end(iOrbitalIndex+k) = iOrbitalIndexSave+nCont
                   else
                      molecule%bas_start(iOrbitalIndex+k) = iOrbitalIndex+1
                      molecule%bas_end(iOrbitalIndex+k) = iOrbitalIndex+kmult
                   endif
                enddo
                iOrbitalIndex = iOrbitalIndex + kmult
             enddo
             kmult = kmult + 2
          ENDIF
       enddo
    enddo

    IF(decinfo%F12)THEN
     call mem_alloc(molecule%atom_cabssize,natoms)
     molecule%atom_cabssize=0

     r = mylsitem%input%basis%binfo(CABBasParam)%labelindex
       
     ! loop over atoms
     do i=1,natoms
      if(r == 0) then
         icharge = int(mylsitem%input%molecule%atom(i)%charge)
         itype = mylsitem%input%basis%binfo(CABBasParam)%chargeindex(icharge)
      else
         itype = mylsitem%input%molecule%atom(i)%idtype(r)
      end if
      molecule%atom_cabssize(i) = &
           & mylsitem%input%basis%binfo(CABBasParam)%atomtype(itype)%TotNOrb
     end do

     ! get first and last index of an atom in ao matrix
     call mem_alloc(molecule%atom_cabsstart,natoms)
     molecule%atom_cabsstart = 0
     molecule%atom_cabsstart(1) = 1
     basis=1
     do i=1,natoms-1
        basis = basis + molecule%atom_cabssize(i)
        molecule%atom_cabsstart(i+1) = basis
     end do
    ENDIF
  end subroutine molecule_get_atomic_sizes

  !> \brief Set occupied and virtual MO orbitals in molecule type
  subroutine molecule_generate_basis(molecule,C)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(inout) :: molecule
    !> All MO coefficients (occupied and virtual)
    real(realk),dimension(molecule%nbasis,molecule%nbasis),intent(in) :: C
    integer :: nbasis,nocc,nvirt,i,j,k
    integer :: tdim(2)
    logical :: loc

    if( molecule%mem_distributed )then
       print *,"WARNING(molecule_generate_basis): going to full, this is not scalable"
    endif
    loc = .not.molecule%mem_distributed
    tdim = [mol_block_size,mol_block_size]

    nbasis = molecule%nbasis
    nocc = molecule%nocc
    nvirt = molecule%nvirt
    call tensor_minit(molecule%Co,[nbasis,nocc ],2, local=loc, atype="TDPD",tdims=tdim)

    ! assign
    !molecule%Co = C(1:nbasis,1:nocc)
    !molecule%Cv = C(1:nbasis,nocc+1:nbasis)
    do j = 1,nocc
       do i = 1,nbasis
          molecule%Co%elm2(i,j) = C(i,j)
       enddo
    enddo

    if(.not.loc)then
       call tensor_mv_dense2tiled(molecule%Co,.true.,dealloc_local=.true.)
    endif


    call tensor_minit(molecule%Cv,[nbasis,nvirt],2, local=loc, atype="TDPD",tdims=tdim)
    k=nocc+1
    do j=1,nvirt
       do i=1,nbasis
          molecule%Cv%elm2(i,j) = C(i,k)
       enddo
       k=k+1
    enddo

    if(.not.loc)then
       call tensor_mv_dense2tiled(molecule%Cv,.true.,dealloc_local=.true.)
    endif

  end subroutine molecule_generate_basis

  !> \brief Get full molecular Fock matrix in MO basis (occupied and virtupied)
  !> \param molecule Full molecule info
  subroutine molecule_mo_fock(molecule)

     implicit none
     type(fullmolecule), intent(inout) :: molecule
     integer :: nocc, nvirt,nbasis
     logical :: loc
     integer :: tdim(2),ord(2)
     type(tensor) :: tmp

     nocc = molecule%nocc
     nvirt = molecule%nvirt
     nbasis = molecule%nbasis

     ord = [1,2]

     loc = .not.molecule%mem_distributed
     tdim = [mol_block_size,mol_block_size]

     ! Occ-occ Fock matrix
     call tensor_minit(molecule%oofock,[nocc,nocc],2, local=loc, atype="TDAR",tdims=tdim)
     !call dec_simple_basis_transform1(nbasis,nocc,molecule%Co,&
     !     & molecule%fock,molecule%oofock)

     call tensor_minit(tmp, [nbasis,nocc], 2, local=loc, atype='TDAR',tdims=tdim  )
     call tensor_contract(1.0E0_realk,molecule%fock,molecule%Co,[2],[1],1,0.0E0_realk,tmp,ord,force_sync=.true.)
     call tensor_contract(1.0E0_realk,molecule%Co,tmp,[1],[1],1,0.0E0_realk,molecule%oofock,ord,force_sync=.true.)
     call tensor_free(tmp)


     ! Virt-virt Fock matrix
     call tensor_minit(molecule%vvfock,[nvirt,nvirt],2, local=loc, atype="TDAR",tdims=tdim)
     !call dec_simple_basis_transform1(nbasis,nvirt,molecule%Cv,&
     !     & molecule%fock,molecule%qqfock)

     call tensor_minit(tmp, [nbasis,nvirt], 2, local=loc, atype='TDAR',tdims=tdim  )
     call tensor_contract(1.0E0_realk,molecule%fock,molecule%Cv,[2],[1],1,0.0E0_realk,tmp,ord,force_sync=.true.)
     call tensor_contract(1.0E0_realk,molecule%Cv,tmp,[1],[1],1,0.0E0_realk,molecule%vvfock,ord,force_sync=.true.)
     call tensor_free(tmp)

  end subroutine molecule_mo_fock

  
  subroutine molecule_mo_f12(MyMolecule,MyLsitem,D)
     type(fullmolecule), intent(inout) :: MyMolecule
     type(lsitem), intent(inout) :: MyLsitem
     type(matrix), intent(in) :: D
#ifdef MOD_UNRELEASED

     integer :: nbasis,nocc,nvirt,noccfull,ncabsAO,nocvfull,ncabsMO

     nbasis   = MyMolecule%nbasis
     nocc     = MyMolecule%nocc
     nvirt    = MyMolecule%nvirt
     noccfull = nocc

     !HACK we do call Fcp for Fcp - indicating 
     !     that this is a Fock(nCabsMO,nbasis)
     !     However all AO -> MO transformations
     !     realted to CABS and RI is postponed
     !     so Fock(nCabsMO,nbasis) will actually be 
     !     Fock(nCabsAO,nbasis) and be a 
     !     half transfomed matrix

     call determine_CABS_nbast(ncabsAO,ncabsMO,MyLsitem%setting,DECinfo%output)
     MyMolecule%nCabsAO = ncabsAO
     MyMolecule%nCabsMO = ncabsMO

     nocvfull = nocc + nvirt

     if(DECinfo%F12debug) then
       ! print *, "--------------------------"
       ! print *, "Molecule_mo_f12"
       ! print *, "--------------------------"
       ! print *, "nbasis:   ", nbasis
       ! print *, "nocc:     ", nocc
       ! print *, "nvirt:    ", nvirt
       ! print *, "--------------------------"
       ! print *, "ncabsAO:  ", ncabsAO
       ! print *, "ncabsMO:  ", ncabsMO
       ! print *, "nocvfull: ", nocc+nvirt
       ! print *, "--------------------------"
     end if

     ! Mixed regular/CABS one-electron  and Coulomb matrix (h+J) combination in AO basis
     call mem_alloc(MyMolecule%hJir,nocc,ncabsAO)    !HACK not RI MO orbitals (AO basis)
     call mem_alloc(MyMolecule%Krs,ncabsAO,ncabsAO)  !HACK not RI MO orbitals (AO basis)
     call mem_alloc(MyMolecule%Fac,nvirt,ncabsAO)    !HACK not nvirt,ncabsMO - not CABS MOs
     call mem_alloc(MyMolecule%Frs,ncabsAO,ncabsAO)  !HACK not RI MO orbitals (AO basis)
     call mem_alloc(MyMolecule%Frm,ncabsAO,noccfull) !HACK not RI MO orbitals (AO basis)
     call mem_alloc(MyMolecule%Fcp,ncabsAO,nbasis)   !HACK not ncabsMO,nbasis - not CABS MOs
     call mem_alloc(MyMolecule%Fij,nocc,nocc)
     !call mem_alloc(MyMolecule%Fcd,ncabsAO,ncabsAO)

     ! Constructing the F12 MO matrices from F12_routines.F90
     call get_F12_mixed_MO_Matrices_real(MyLsitem,MyMolecule,D,nbasis,ncabsAO,&
        & nocc,noccfull,nvirt,MyMolecule%hJir,MyMolecule%Krs,MyMolecule%Frs,&
        & MyMolecule%Fac,MyMolecule%Fij,MyMolecule%Frm,MyMolecule%Fcp)

#endif
  end subroutine molecule_mo_f12


  !> \brief Calculate how much memory is used for the fullmolecule type (in GB).
  !> Only two-dimensional arrays are considered and memory for
  !> vectors and single elements are ignored.
  !> \author Kasper Kristensen
  !> \date February 2011
  subroutine calculate_fullmolecule_memory(MyMolecule,molmem)

    implicit none
    !> Full molecule information
    type(fullmolecule), intent(inout) :: MyMolecule
    real(realk), intent(inout) :: molmem
    real(realk) :: O,V,A,tmp
    integer :: nnod, nblocks, nlocal_blocks
    ! GB conversion
    real(realk), parameter :: GB=1.024E3_realk**3! 1GB


    ! Number of occupied (O), Virtual (V), atomic basis functions (A)
    ! ***************************************************************
    O = MyMolecule%nocc
    V = MyMolecule%nvirt
    A = MyMolecule%nbasis


    ! Use type fullmolecule to calculate memory use
    ! *********************************************

    ! MO coefficients, Fock
    molmem = 2E0_realk*A*A

    ! ppfock and qqfock
    tmp = O*O + V*V
    molmem = molmem + tmp


#ifdef VAR_MPI
    !Do we need to distribute the arrays? -> if more than 10% of the available memory, the memory will be distributed
    if(molmem>((0.1*DECinfo%memory*GB)/realk))then
       MyMolecule%mem_distributed = .true.
    else
       MyMolecule%mem_distributed = .false.
    endif

    if(Decinfo%force_distribution)then
       MyMolecule%mem_distributed = DECinfo%distribute_fullmolecule
    endif

    nnod = infpar%nodtot
#else
    MyMolecule%mem_distributed = .false.
    nnod = 1
#endif

    if(MyMolecule%mem_distributed)then
       if( matrix_type == mtype_pdmm )then
          mol_block_size = BLOCK_SIZE_PDM
       else
          !Block size to about 100MB, this is a compile--time constant
          mol_block_size = int(sqrt((0.1*GB)/dble(realk)))
       endif

       !re-evaluate the local memory requirements
       nblocks       = molmem/mol_block_size
       !this depends on the distribution, with the one chosen for now, the
       !following measure holds
       nlocal_blocks = ceiling(float(nblocks)/float(nnod))
       
       molmem        = nlocal_blocks
    endif

    ! Convert to GB
    molmem = realk*molmem/GB

  end subroutine calculate_fullmolecule_memory


  !> \brief Get density matrix in type matrix form by reading from file.
  !> \author Kasper Kristensen
  !> \date April 2013
  subroutine dec_get_density_matrix_from_file(nbasis,D)
    implicit none
    !> Number of basis functions in molecule
    integer,intent(in) :: nbasis
    !> Density matrix (will be intialized here)
    type(matrix),intent(inout) :: D
    integer :: funit,dim1,dim2
    integer(kind=long) :: longdim1,longdim2
    real(realk),pointer :: Dfull(:,:)
    logical :: gcbasis
    logical(kind=8) :: gcbasis8
    ! Open density file
    funit=-1
    call lsopen(funit,'dens.restart','OLD','UNFORMATTED')

    ! Allocate real vector to temporarily hold density values
    call mem_alloc(Dfull,nbasis,nbasis)

    READ(funit) longdim1,longdim2 !dens.restart always written using kind=8
    dim1 = int(longdim1)
    dim2 = int(longdim2)

    ! Sanity check
    if( (dim1/=nbasis) .or. (dim2/=nbasis) ) then
       print *, 'nbasis: ', nbasis
       print *, 'dim1,dim2', dim1,dim2
       call lsquit('dec_get_density_matrix_from_file: Error in density dimensions!',-1)
    end if

    ! Read density elements
    read(funit) Dfull

    read(funit) gcbasis8
    gcbasis = gcbasis8

    ! Basis set Sanity check
    if (gcbasis .and. .not. DECinfo%GCBASIS) then
        WRITE(DECinfo%output,*) 'Your dens.restart was constructed using the grand-canonical (GC) basis,'
        WRITE(DECinfo%output,*) 'while your LSDALTON.INP uses the standard input basis. '
        WRITE(DECinfo%output,*) 'The GC basis is default unless you use a dunnings basis set,'
        WRITE(DECinfo%output,*) 'or you specify .NOGCBASIS under *GENERAL'
        call lsquit('Calculation in standard basis, dens.restart in GC basis!',DECinfo%output)
     else if (DECinfo%GCBASIS .and. .not. gcbasis) then
        WRITE(DECinfo%output,*) 'Your dens.restart was constructed using the standard input basis, while your'
        WRITE(DECinfo%output,*) 'LSDALTON.INP uses the grand-canonical (GC) basis.'
        WRITE(DECinfo%output,*) 'The GC basis is default unless you use a dunnings basis set,'
        WRITE(DECinfo%output,*) 'or you specify .NOGCBASIS under *GENERAL'
        call lsquit('Calculation in GC basis, dens.restart in standard input basis!',DECinfo%output)
     else if (DECinfo%GCBASIS .and. gcbasis) then
        WRITE(DECinfo%output,*) 'Basis check ok: Using GC basis consistently'
     else if (.not. DECinfo%GCBASIS .and. .not. gcbasis) then
        WRITE(DECinfo%output,*) 'Basis check ok: Using standard basis consistently'
     else
        call lsquit('Basis check is messed up!!',DECinfo%output)
     endif

     call lsclose(funit,'KEEP')

    ! Init density matrix
    call mat_init(D,nbasis,nbasis)

    call mat_set_from_full(Dfull,1.0_realk,D)
    call mem_dealloc(Dfull)

  end subroutine dec_get_density_matrix_from_file
  

  ! THIS ROUTINE SHOULD BE RECONSIDERED IF WE FIND A GOOD ORBITAL INTERACTION MATRIX TO USE
  ! FOR FRAGMENT EXPANSION:   
  !> Calculate occ and virt interaction matrices which are used for atomic fragment
  !> optimization as a measure of orbital interactions:
  !> G_ia = sum_{mu nu} C_{mu i} C_{nu a} G_{mu nu}
  !> G_{mu nu} is approximation to (mu nu | mu nu)
  !> More precisely G_{mu nu} > (mu nu | mu nu), see II_get_2int_ScreenMat.
  !> Thus G_{ia} is a rough measure of (i a | i a) integrals which enter the CC energy, and therefore
  !> G_{ia} is a measure of the interaction between orbital "i" and orbital "a" 
  !> in a correlation energy context.
!!$  subroutine molecule_get_interaction_matrices(molecule,mylsitem)
!!$    implicit none
!!$
!!$    !> Molecule info
!!$    type(fullmolecule), intent(inout) :: molecule
!!$    !> Integral info
!!$    type(lsitem), intent(inout) :: mylsitem
!!$    real(realk),pointer :: AOint(:,:)
!!$    type(matrix),target :: AOint_mat
!!$    integer :: i,j,idx
!!$
!!$
!!$    ! Get interaction matrix in AO basis
!!$    call mat_init(AOint_mat,molecule%nbasis,molecule%nbasis)
!!$    call mat_zero(AOint_mat)
!!$    call II_get_2int_ScreenMat(DECinfo%output,DECinfo%output,mylsitem%setting,AOint_mat)
!!$
!!$    ! Temporary fix until Thomas makes routine which gives Fortran arrays as output
!!$    call mem_alloc(AOint,molecule%nbasis,molecule%nbasis)
!!$    idx=0
!!$    do j=1,molecule%nbasis
!!$       do i=1,molecule%nbasis
!!$          idx = idx+1
!!$          AOint(i,j) = AOint_mat%elms(idx)
!!$       end do
!!$    end do
!!$    call mat_free(AOint_mat)
!!$
!!$    ! Init interaction matrix: Occupied,virtual dimension
!!$    call mem_alloc(molecule%orbint,molecule%nocc,molecule%nvirt)
!!$
!!$    ! Transform to MO basis
!!$    call dec_diff_basis_transform1(molecule%nbasis,molecule%nocc,molecule%nvirt,&
!!$         & molecule%Co, molecule%Cv, AOint, molecule%orbint)
!!$
!!$    ! Take absolute value (should not be necessary but do it to be on the safe side)
!!$    do j=1,molecule%nvirt
!!$       do i=1,molecule%nocc
!!$          molecule%orbint(i,j) = abs(molecule%orbint(i,j))
!!$       end do
!!$    end do
!!$
!!$    call mem_dealloc(AOint)
!!$
!!$  end subroutine molecule_get_interaction_matrices


end module full_molecule
