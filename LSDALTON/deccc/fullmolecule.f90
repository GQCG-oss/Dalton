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
  use files!,only:lsopen,lsclose
  use matrix_module!, only:matrix
  use matrix_util!,only:mat_init,mat_zero,mat_daxpy,mat_free,mat_to_full,&
!       & mat_diag_f, mat_write_to_disk
  use matrix_operations!, only: mat_set_from_full,mat_free,mat_init,mat_to_full
  use dec_typedef_module
  use IntegralInterfaceMod


  ! DEC DEPENDENCIES (within deccc directory) 
  ! *****************************************
  use dec_fragment_utils
  use array2_simple_operations

contains

  !> \brief Initialize informations about full molecule by reading HF info from file.
  !> NOTE: If this routine is modified, then molecule_init_from_inputs must be modified
  !> in the same manner!
  !> \author Marcin Ziolkowski
  !> \param molecule Full molecule info
  !> \param mylsitem Integral program input
  subroutine molecule_init_from_files(molecule,mylsitem)

    implicit none
    type(fullmolecule), intent(inout) :: molecule
    type(lsitem), intent(inout) :: mylsitem
    integer :: natoms,basis,i
    integer :: r,iset,itype
    logical :: status_info
    real(realk) :: memory_use, tcpu, twall

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
    call molecule_get_overlap(molecule,mylsitem)
    call molecule_mo_fock(molecule)

    if(DECinfo%use_canonical) then ! overwrite local orbitals and use canonical orbitals
       call dec_get_canonical_orbitals(molecule)
    end if

    call LSTIMER('DEC: MOL INIT',tcpu,twall,DECinfo%output)

  end subroutine molecule_init_from_files



  !> \brief Initialize informations about full molecule
  !> using input Fock,MO, and overlap matrix 
  !> NOTE: If this routine is modified, then molecule_init_from_files must be modified
  !> in the same manner!
  !> \author Kasper Kristensen
  !> \date November 2011
  subroutine molecule_init_from_inputs(molecule,mylsitem,F,S,C)

    implicit none
    !> Full molecule structure to be initialized
    type(fullmolecule), intent(inout) :: molecule
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> Fock matrix
    type(matrix),intent(in) :: F
    !> Overlap matrix
    type(matrix),intent(in) :: S
    !> MO coefficients
    type(matrix),intent(in) :: C
    real(realk) :: memory_use, tcpu, twall

    call LSTIMER('START',tcpu,twall,DECinfo%output)

    ! Init basic info (molecular dimensions etc.)
    call molecule_init_basics(molecule,mylsitem)

    ! Copy Fock, density, MO, and overlap matrices to molecule structure
    call molecule_copy_FSC_matrices(molecule,F,S,C)

    ! Fock matrix in MO basis
    call molecule_mo_fock(molecule)

    if(DECinfo%use_canonical) then ! overwrite local orbitals and use canonical orbitals
       call dec_get_canonical_orbitals(molecule)
    end if

    call LSTIMER('DEC: MOL INIT',tcpu,twall,DECinfo%output)

  end subroutine molecule_init_from_inputs



  !> \brief Initialize basic informations about full molecule (number of basis functions,
  !> number of occipied orbitals etc.).
  !> \author Marcin Ziolkowski/Kasper Kristensen
  !> \param molecule Full molecule info
  !> \param mylsitem Integral program input
  subroutine molecule_init_basics(molecule,mylsitem)

    implicit none
    type(fullmolecule), intent(inout) :: molecule
    type(lsitem), intent(inout) :: mylsitem
    integer :: natoms,basis,i
    integer :: r,iset,itype
    logical :: status_info
    real(realk) :: memory_use, tcpu, twall

    call LSTIMER('START',tcpu,twall,DECinfo%output)

    molecule%natoms = get_num_atoms(mylsitem)
    molecule%nelectrons = get_num_electrons(mylsitem)
    molecule%nbasis = get_num_basis_functions(mylsitem)
    molecule%nauxbasis = get_num_aux_basis_functions(mylsitem)
    molecule%numocc = molecule%nelectrons/2
    molecule%numvirt = molecule%nbasis - molecule%numocc
    molecule%ncore = count_ncore(mylsitem)
    molecule%nval = molecule%numocc - molecule%ncore

    ! Which basis functions are on which atoms?
    call molecule_get_atomic_sizes(molecule,mylsitem)

    ! Memory use for full molecule
    call calculate_fullmolecule_memory(molecule,memory_use)
    DECinfo%fullmolecule_memory = memory_use

    ! Print some info about the molecule
    write(DECinfo%output,*)
    write(DECinfo%output,'(/,a)') '-- Full moleculecular info --'
    write(DECinfo%output,'(/,a,i6)') 'FULL: Overall charge of molecule : ',nint(mylsitem%input%molecule%charge)
    write(DECinfo%output,'(/,a,i6)') 'FULL: Number of electrons        : ',molecule%nelectrons
    write(DECinfo%output,'(a,i6)')   'FULL: Number of atoms            : ',molecule%natoms
    write(DECinfo%output,'(a,i6)')   'FULL: Number of basis func.      : ',molecule%nbasis
    write(DECinfo%output,'(a,i6)')   'FULL: Number of aux. basis func. : ',molecule%nauxbasis
    write(DECinfo%output,'(a,i6)')   'FULL: Number of core orbitals    : ',molecule%ncore
    write(DECinfo%output,'(a,i6)')   'FULL: Number of valence orbitals : ',molecule%nval
    write(DECinfo%output,'(a,i6)')   'FULL: Number of occ. orbitals    : ',molecule%numocc
    write(DECinfo%output,'(a,i6)')   'FULL: Number of virt. orbitals   : ',molecule%numvirt
    write(DECinfo%output,*)

  end subroutine molecule_init_basics


  !> \brief Copy Fock,density,MO, and overlap matrices in type(matrix) format to molecule structure.
  !> \author Kasper Kristensen
  !> \date November 2011
  subroutine molecule_copy_FSC_matrices(molecule,F,S,C)

    implicit none
    !> Full molecule structure to be initialized
    type(fullmolecule), intent(inout) :: molecule
    !> Fock matrix
    type(matrix),intent(in) :: F
    !> Overlap matrix
    type(matrix),intent(in) :: S
    !> MO coefficients
    type(matrix),intent(in) :: C
    integer :: nbasis,i
    real(realk),pointer :: basis(:,:)

    nbasis = molecule%nbasis

    ! Fock matrix
    call mem_alloc(molecule%fock,nbasis,nbasis)
    call mat_to_full(F, 1.0_realk, molecule%fock)

    ! Overlap matrix
    call mem_alloc(molecule%overlap,nbasis,nbasis)
    call mat_to_full(S, 1.0_realk, molecule%overlap)

    ! MO coefficient matrix
    call mem_alloc(basis,nbasis,nbasis)
    call mat_to_full(C, 1.0_realk, basis)
    call molecule_generate_basis(molecule,basis)
    call mem_dealloc(basis)

  end subroutine molecule_copy_FSC_matrices


  !> \brief Copy Fock,density,MO, and overlap matrices in molecule strucutre to type(matrix) format.
  !> (included intialization of matrices)
  !> \author Kasper Kristensen
  !> \date December 2012
  subroutine molecule_copyback_FSC_matrices(molecule,F,S,C)

    implicit none
    !> Full molecule structure
    type(fullmolecule), intent(in) :: molecule
    !> Fock matrix
    type(matrix),intent(inout) :: F
    !> Overlap matrix
    type(matrix),intent(inout) :: S
    !> MO coefficients
    type(matrix),intent(inout) :: C
    integer :: nbasis,i
    real(realk),pointer :: tmp(:,:)

    nbasis = molecule%nbasis
    call mat_init(F,nbasis,nbasis)
    call mat_init(S,nbasis,nbasis)
    call mat_init(C,nbasis,nbasis)
    call mat_set_from_full(Molecule%fock,1.0_realk,F)
    call mat_set_from_full(Molecule%overlap,1.0_realk,S)

    call mem_alloc(tmp,nbasis,nbasis)
    ! Put occ orbitals into tmp
    do i=1,Molecule%numocc
       tmp(1:nbasis,i) = Molecule%ypo(1:nbasis,i)
    end do
    ! Put virt orbitals into tmp
    do i=1,Molecule%numvirt
       tmp(1:nbasis,i+Molecule%numocc) = Molecule%ypv(1:nbasis,i)
    end do
    
    ! All orbitals into C
    call mat_set_from_full(tmp,1.0_realk,C)
    call mem_dealloc(tmp)


  end subroutine molecule_copyback_FSC_matrices



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
       electrons = electrons + mylsitem%input%molecule%Atom(i)%Charge
    end do
    charge = nint(mylsitem%input%molecule%charge)
    electrons = electrons - charge

    return
  end function get_num_electrons

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


  !> \brief Read or construct (if it does not already exit) the fock matrix
  subroutine molecule_get_fock(molecule,mylsitem)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(inout) :: molecule
    !> LSitem
    type(lsitem), intent(inout) :: mylsitem
    integer :: nbasis
    logical :: fock_exist

    ! Init stuff
    nbasis = molecule%nbasis
    call mem_alloc(molecule%fock,nbasis,nbasis)
    inquire(file='fock.restart',exist=fock_exist)

    ! Read or construct fock matrix
    ! *****************************
    if(fock_exist) then

       ! Read Fock matrix from file
       write(DECinfo%output,*) 'Reading Fock matrix from file fock.restart...'
       write(DECinfo%output,*)
       call dec_read_mat_from_file('fock.restart',nbasis,nbasis,molecule%fock)

    else

       call lsquit('molecule_get_fock: Fock matrix does not exist!',-1)

    end if


  end subroutine molecule_get_fock

  !> \brief Set orbitals used in DEC to canonical orbitals (only for testing).
  !> Assumes that Fock matrix is stored in molecule%fock when this subroutine is called.
  !> After this call molecule%ypo and molecule%ypv will contain the occupied and virtual MO
  !> coefficients, respectively, while molecule%ppfock and molecule%qqfock will be the occ-occ
  !> and virt-virt blocks of the diagonal canonical MO Fock matrix.
  !> \author Kasper Kristensen
  !> \date September 2011
  subroutine dec_get_canonical_orbitals(molecule)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(inout) :: molecule
    integer :: nbasis,i,nocc,nunocc
    real(realk), pointer :: eival(:), C(:,:)

    nbasis = molecule%nbasis
    nocc = molecule%numocc
    nunocc = molecule%numvirt

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) 'Using canonical orbitals as requested in input!'
    write(DECinfo%output,*) 'Warning: This may cause meaningless results when the DEC scheme is used!'
    write(DECinfo%output,*)

    ! Canonical MO coefficients
    call mem_alloc(C,nbasis,nbasis)

    ! Orbital energies
    call mem_alloc(eival,nbasis)

    ! Diagonalize Fock matrix
    call solve_eigenvalue_problem(nbasis,molecule%fock,molecule%overlap,eival,C)

    ! Set MO coefficients
    Molecule%ypo = C(:,1:nocc)   ! occupied 
    Molecule%ypv = C(:,nocc+1:nbasis)   ! unoccupied

    ! Set Fock matrix in canonical MO basis 
    Molecule%ppfock=0.0_realk
    Molecule%qqfock=0.0_realk
    do i=1,nocc
       molecule%ppfock(i,i) = eival(i)
    end do
    do i=1,nunocc
       molecule%qqfock(i,i) = eival(i+nocc)
    end do

    write(DECinfo%output,*) 'Orbital energies:'
    do i=1,nbasis
       write(DECinfo%output,*) i, eival(i)
    end do

    call mem_dealloc(C)
    call mem_dealloc(eival)

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
    call dec_read_mat_from_file('lcm_orbitals.u',nbasis,nbasis,C)
    call molecule_generate_basis(molecule,C)
    call mem_dealloc(C)

  end subroutine molecule_get_reference_state





  !> \brief Read or construct (if it does not already exist) the overlap matrix
  !> \author Kasper Kristensen
  !> \date January 2012
  subroutine molecule_get_overlap(molecule,mylsitem)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(inout) :: molecule
    !> LSitem
    type(lsitem), intent(inout) :: mylsitem
    type(matrix) :: S
    integer :: nbasis
    logical :: overlap_exist

    ! Init stuff
    nbasis = molecule%nbasis
    call mem_alloc(molecule%overlap,nbasis,nbasis)
    inquire(file='overlapmatrix',exist=overlap_exist)

    ! Read or construct overlap matrix
    ! ********************************
    if(overlap_exist) then

       ! Read overlap matrix from file
       write(DECinfo%output,*) 'Reading overlap matrix from file overlapmatrix...'
       write(DECinfo%output,*)
       call dec_read_mat_from_file('overlapmatrix',nbasis,nbasis,molecule%overlap)

    else

       ! Calculate overlap matrix from scratch
       write(DECinfo%output,*) 'Calculating overlap matrix for DEC calculation...'
       call mat_init(s,nbasis,nbasis)
       call mat_zero(s)
       call II_get_overlap(DECinfo%output,DECinfo%output,mylsitem%setting,s)
       call mem_alloc(molecule%overlap,nbasis,nbasis)
       molecule%overlap=0.0E0_realk
       call mat_to_full(s,1.0E0_realk,molecule%overlap)
       call mat_free(s)

    end if


  end subroutine molecule_get_overlap


  !> \brief Destroy fullmolecule structure
  !> \param molecule Full molecular info
  subroutine molecule_finalize(molecule)

    implicit none
    type(fullmolecule), intent(inout) :: molecule

    ! Do not free anything if it has not been initiated (mainly for testing)
    if(DECinfo%SkipReadIn) return

    ! Delete transformation matrices for general basis
    if(associated(molecule%ypo)) then
       call mem_dealloc(molecule%ypo)
    end if

    if(associated(molecule%ypv)) then
       call mem_dealloc(molecule%ypv)
    end if

    ! Delete AO fock matrix
    if(associated(molecule%fock)) then
       call mem_dealloc(molecule%fock)
    end if

    ! P^Fock
    if(associated(molecule%ppfock)) then
       call mem_dealloc(molecule%ppfock)
    end if

    ! Q^Fock
    if(associated(molecule%qqfock)) then
       call mem_dealloc(molecule%qqfock)
    end if

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

    if(associated(molecule%overlap)) then
       call mem_dealloc(molecule%overlap)
    end if

  end subroutine molecule_finalize

  !> \brief Get number of atomic orbitals on atoms, first and last index in AO basis for full molecular matrices
  !> \param molecule Full molecule info
  !> \param mylsitem Integral program input
  subroutine molecule_get_atomic_sizes(molecule,mylsitem)

    implicit none
    type(fullmolecule), intent(inout) :: molecule
    type(lsitem), intent(inout) :: mylsitem
    integer :: natoms,r,i,iset,itype,basis,icharge
    logical :: status_info

    natoms = molecule%natoms
    call mem_alloc(molecule%atom_size,natoms)
    molecule%atom_size=0

    r = mylsitem%input%basis%regular%labelindex

    ! loop over atoms
    do i=1,natoms

       if(r == 0) then
          icharge = int(mylsitem%input%molecule%atom(i)%charge)
          itype = mylsitem%input%basis%regular%chargeindex(icharge)
       else
          itype = mylsitem%input%molecule%atom(i)%idtype(1)
       end if

       molecule%atom_size(i) = mylsitem%input%basis%regular%&
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

    return
  end subroutine molecule_get_atomic_sizes

  !> \brief Set occupied and virtual MO orbitals in molecule type
  subroutine molecule_generate_basis(molecule,C)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(inout) :: molecule
    !> All MO coefficients (occupied and virtual)
    real(realk),dimension(molecule%nbasis,molecule%nbasis),intent(in) :: C
    integer :: nbasis,nocc,nvirt

    nbasis = molecule%nbasis
    nocc = molecule%numocc
    nvirt = molecule%numvirt
    call mem_alloc(molecule%ypo,nbasis,nocc)
    call mem_alloc(molecule%ypv,nbasis,nvirt)

    ! assign
    molecule%ypo = C(1:nbasis,1:nocc)
    molecule%ypv = C(1:nbasis,nocc+1:nbasis)

  end subroutine molecule_generate_basis

  !> \brief Get full molecular Fock matrix in MO basis (occupied and unoccupied)
  !> \param molecule Full molecule info
  subroutine molecule_mo_fock(molecule)

    implicit none
    type(fullmolecule), intent(inout) :: molecule
    type(array2) :: ppfock, qqfock, ypo,ypv,yho,yhv,fock
    integer :: nocc, nvirt, oo(2), bo(2), bv(2), vv(2), bb(2),nbasis

    nocc = molecule%numocc
    nvirt = molecule%numvirt
    nbasis = molecule%nbasis
    oo(1)=nocc
    oo(2)=nocc
    vv(1)=nvirt
    vv(2)=nvirt
    bo(1)=nbasis
    bo(2)=nocc
    bv(1)=nbasis
    bv(2)=nvirt
    bb(1)=nbasis
    bb(2)=nbasis

    ! Fock matrix in AO basis
    fock = array2_init(bb,molecule%fock)

    ! Occ-occ block
    ypo = array2_init(bo,molecule%ypo)
    yho = array2_init(bo,molecule%ypo)
    ppfock = array2_similarity_transformation(ypo,fock,yho,oo)
    call array2_free(ypo)
    call array2_free(yho)
    call mem_alloc(molecule%ppfock,nocc,nocc)
    molecule%ppfock(1:nocc,1:nocc) = ppfock%val(1:nocc,1:nocc)
    call array2_free(ppfock)

    ! Virt-virt block
    ypv = array2_init(bv,molecule%ypv)
    yhv = array2_init(bv,molecule%ypv)
    qqfock = array2_similarity_transformation(ypv,fock,yhv,vv)
    call array2_free(ypv)
    call array2_free(yhv)
    call array2_free(fock)
    call mem_alloc(molecule%qqfock,nvirt,nvirt)
    molecule%qqfock(1:nvirt,1:nvirt) = qqfock%val(1:nvirt,1:nvirt)
    call array2_free(qqfock)

  end subroutine molecule_mo_fock



  !> \brief Calculate how much memory is used for the fullmolecule type (in GB).
  !> Only two-dimensional arrays are considered and memory for
  !> vectors and single elements are ignored.
  !> \author Kasper Kristensen
  !> \date February 2011
  subroutine calculate_fullmolecule_memory(MyMolecule,molmem)

    implicit none
    !> Full molecule information
    type(fullmolecule), intent(in) :: MyMolecule
    real(realk), intent(inout) :: molmem
    real(realk) :: O,V,A,tmp,GB

    ! GB conversion
    GB = 1.000E9_realk ! 1 GB


    ! Number of occupied (O), Virtual (V), atomic basis functions (A)
    ! ***************************************************************
    O = MyMolecule%numocc
    V = MyMolecule%numvirt
    A = MyMolecule%nbasis


    ! Use type fullmolecule to calculate memory use
    ! *********************************************

    ! MO coefficients, Fock, overlap
    molmem = 3E0_realk*A*A

    ! ppfock and qqfock
    tmp = O*O + V*V
    molmem = molmem + tmp

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
    integer(kind=8) :: longdim1,longdim2
    integer(kind=4) :: dim1_32,dim2_32
    real(realk),pointer :: Dfull(:,:)

    ! Open density file
    funit=-1
    call lsopen(funit,'dens.restart','OLD','UNFORMATTED')

    ! Allocate real vector to temporarily hold density values
    call mem_alloc(Dfull,nbasis,nbasis)


    ! Safe handling of integers 
    if(DECinfo%convert64to32) then   ! convert from 64 to 32 bit integers?
       READ(funit) longdim1,longdim2
       dim1 = int(longdim1)
       dim2 = int(longdim2)
    elseif(DECinfo%convert32to64) then   ! convert from 32 to 64 bit integers?
       READ(funit) dim1_32,dim2_32
       dim1 = dim1_32
       dim2 = dim2_32
    else
       READ(funit) dim1,dim2
    end if

    ! Sanity check
    if( (dim1/=nbasis) .or. (dim2/=nbasis) ) then
       print *, 'nbasis: ', nbasis
       print *, 'dim1,dim2', dim1,dim2
       call lsquit('dec_get_density_matrix_from_file: Error in density dimensions!',-1)
    end if

    ! Read density elements
    read(funit) Dfull
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
!!$    call mem_alloc(molecule%orbint,molecule%numocc,molecule%numvirt)
!!$
!!$    ! Transform to MO basis
!!$    call dec_diff_basis_transform1(molecule%nbasis,molecule%numocc,molecule%numvirt,&
!!$         & molecule%ypo, molecule%ypv, AOint, molecule%orbint)
!!$
!!$    ! Take absolute value (should not be necessary but do it to be on the safe side)
!!$    do j=1,molecule%numvirt
!!$       do i=1,molecule%numocc
!!$          molecule%orbint(i,j) = abs(molecule%orbint(i,j))
!!$       end do
!!$    end do
!!$
!!$    call mem_dealloc(AOint)
!!$
!!$  end subroutine molecule_get_interaction_matrices


end module full_molecule