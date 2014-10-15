!> @file
!> Utils for DEC subroutines
!> \author Marcin Ziolkowski (modified by Kasper Kristensen)
module snoop_tools_module

  use fundamental
  use precision
  use lstiming
  use ls_util!,only: dgemm_ts
  use typedeftype!,only: lsitem
  use molecule_module!, only: get_geometry
  use files!,only:lsopen,lsclose
  use DALTONINFO!, only: ls_free
  use dec_typedef_module
  use memory_handling!, only: mem_alloc, mem_dealloc, mem_allocated_global,&
  !       & stats_mem, get_avaiLable_memory
  use,intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
  use matrix_module!, only:matrix
  use matrix_operations
  use matrix_util
  use matrix_operations_aux
  use tensor_interface_module
  use array4_simple_operations
  use IntegralInterfaceMOD!, only: ii_get_h1, ii_get_nucpot
  use BUILDAOBATCH
  use DALTONINFO

#ifdef VAR_MPI
  use infpar_module
  use lsmpi_op
#endif

#ifdef VAR_PAPI
  use papi_module
#endif


contains


  !> Build lsitem for subsystem with no ghost functions from the other subsystems
  !> \author Kasper Kristensen
  !> \date August 2014
  subroutine build_subsystem_lsitem_no_ghost(this,MyMolecule,lsfull,lssub)
    implicit none
    !> Subsystem label
    integer,intent(in) :: this
    !> Full molecule info
    type(fullmolecule),intent(in) :: MyMolecule
    !> LSitem for full system
    type(lsitem),intent(in) :: lsfull
    !> lsitem for subsystem
    type(lsitem),intent(inout) :: lssub
    integer,pointer :: subatoms(:)
    integer :: natomssub


    ! Number of atoms for subsystem and list of these atoms
    call get_natoms_for_subsystem(this,MyMolecule,natomssub)
    call mem_alloc(subatoms,natomssub)
    call list_of_atoms_for_subsystem(this,MyMolecule,natomssub,subatoms)

    ! Build lssub with atoms for subsystem
    call build_ccfragmentlsitem(lsfull,lssub,subatoms,natomssub,DECinfo%output,DECinfo%output)
    call mem_dealloc(subatoms)

  end subroutine build_subsystem_lsitem_no_ghost


  !> Get number of number of atoms for subsystem
  !> \author Kasper Kristensen
  !> \date August 2014
  subroutine get_natoms_for_subsystem(this,MyMolecule,natomssub)
    implicit none
    !> Subsystem label
    integer,intent(in) :: this
    !> Molecule info
    type(fullmolecule),intent(in) :: MyMolecule
    !> Number of atoms for subsystem "this"
    integer,intent(inout) :: natomssub
    integer :: i

    ! Number of atoms in subsystem "this"
    ! ***********************************
    natomssub=0
    do i=1,MyMolecule%natoms
       if(MyMolecule%SubSystemIndex(i)==this) then
          natomssub = natomssub + 1
       end if
    end do

  end subroutine get_natoms_for_subsystem



  !> Get list of atoms for a given subsystem
  !> \author Kasper Kristensen
  !> \date August 2014
  subroutine list_of_atoms_for_subsystem(this,MyMolecule,natomssub,subatoms)
    implicit none

    !> Subsystem label
    integer,intent(in) :: this
    !> Molecule info
    type(fullmolecule),intent(in) :: MyMolecule
    !> Number of atoms in subsystem
    integer,intent(in) :: natomssub
    !> List of atoms in subsystem
    integer,intent(inout) :: subatoms(natomssub)
    integer :: i,idx

    idx=0
    do i=1,MyMolecule%natoms
       if(MyMolecule%SubSystemIndex(i)==this) then
          idx=idx+1
          subatoms(idx) = i
       end if
    end do

    ! Sanity check
    if(idx/=natomssub) then
       print *, 'idx,natomssub',idx,natomssub
       call lsquit('list_of_atoms_for_subsystem: Counter mismatch',-1)
    end if

  end subroutine list_of_atoms_for_subsystem


  !> \brief Obtain initial orbital guess by diagonalizing one-el. Hamiltonian
  !> Can be used also for subsystems as long as mylsitem corresponds to subsystem
  !> \author Kasper Kristensen
  !> \date 2014
  subroutine starting_orbitals_from_h1(mylsitem,CMO)
    implicit none
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> MO coefficients corresponding to diagonalization of H1
    type(matrix),intent(inout)  :: CMO
    type(matrix) :: H1,S
    real(realk), pointer    :: eival(:)
    integer :: nbasis

    ! Get one-electron matrix and overlap
    nbasis = CMO%nrow
    call mat_init(H1,nbasis,nbasis)
    call mat_init(S,nbasis,nbasis)
    call ii_get_h1(DECinfo%output,DECinfo%output,mylsitem%setting,H1)
    call II_get_overlap(DECinfo%output,DECinfo%output,mylsitem%setting,S)

    ! Diagonalize one-electron matrix
    call mem_alloc(eival,nbasis)
    call mat_diag_f(H1,S,eival,Cmo)
    call mem_dealloc(eival)
    call mat_free(H1)
    call mat_free(S)

  end subroutine starting_orbitals_from_h1


  !> \brief Obtain initial orbital guess by diagonalizing one-el. Hamiltonian, where
  !> the orbitals are split into occ and virt matrices.
  !> Can be used also for subsystems as long as mylsitem corresponds to subsystem
  !> \author Kasper Kristensen
  !> \date 2014
  subroutine starting_orbitals_from_h1_split(mylsitem,Cocc,Cvirt)
    implicit none
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> MO coefficients corresponding to diagonalization of H1
    type(matrix),intent(inout)  :: Cocc,Cvirt
    type(matrix) :: CMO
    real(realk), pointer    :: eival(:)
    integer :: nbasis

    ! MO coefficients collected
    nbasis = Cocc%nrow
    call mat_init(CMO,nbasis,nbasis)
    call starting_orbitals_from_h1(mylsitem,CMO)

    ! Split into two matrices
    call partition_MO_coeff_into_two_matrices(CMO,Cocc,Cvirt)
    call mat_free(CMO)

  end subroutine starting_orbitals_from_h1_split


  !> Collect occ and virt MO coefficients in one type(matrix), where the number of
  !> basis functions is not necessarily equal to nocc+nvirt.
  !> \author Kasper Kristensen
  !> \date August 2014
  subroutine collect_MO_coeff_in_one_matrix(Cocc,Cvirt,C)
    implicit none
    !> Occ and virt MO coeff
    type(matrix),intent(in) :: Cocc,Cvirt
    !> Occ and virt MO coeff combined
    type(matrix),intent(inout) :: C
    integer :: mu,p,idx,offset

    if(matrix_type/=mtype_dense) then
       call lsquit('collect_MO_coeff_in_one_matrix: Only implemented for dense matrices!',-1)
    end if

    ! Copy occ orbitals
    idx=0
    do mu=1,Cocc%nrow  ! AO loop
       do p=1,Cocc%ncol   ! occ MO loop
          idx=idx+1
          C%elms(idx) = Cocc%elms(idx)
       end do
    end do
    offset = idx   ! Number of elements for occ orbitals

    ! Copy virt orbitals
    idx=0
    do mu=1,Cvirt%nrow  ! AO loop
       do p=1,Cvirt%ncol   ! virt MO loop
          idx=idx+1
          C%elms(idx+offset) = Cvirt%elms(idx)
       end do
    end do

  end subroutine collect_MO_coeff_in_one_matrix



  !> Partion MO coefficients collected in one matrix into an occupied and a virtual MO matrix.
  !> (the reverse of collect_MO_coeff_in_one_matrix).
  !> \author Kasper Kristensen
  !> \date August 2014
  subroutine partition_MO_coeff_into_two_matrices(C,Cocc,Cvirt)
    implicit none
    !> Occ and virt MO coeff combined
    type(matrix),intent(in) :: C
    !> Occ and virt MO coeff
    type(matrix),intent(inout) :: Cocc,Cvirt
    integer :: mu,p,idx,offset

    if(matrix_type/=mtype_dense) then
       call lsquit('partition_MO_coeff_into_two_matrices: Only implemented for dense matrices!',-1)
    end if

    ! Copy occ orbitals
    idx=0
    do mu=1,Cocc%nrow  ! AO loop
       do p=1,Cocc%ncol   ! occ MO loop
          idx=idx+1
          Cocc%elms(idx) = C%elms(idx) 
       end do
    end do
    offset = idx   ! Number of elements for occ orbitals

    ! Copy virt orbitals
    idx=0
    do mu=1,Cvirt%nrow  ! AO loop
       do p=1,Cvirt%ncol   ! virt MO loop
          idx=idx+1
          Cvirt%elms(idx) = C%elms(idx+offset) 
       end do
    end do

  end subroutine partition_MO_coeff_into_two_matrices


end module snoop_tools_module
