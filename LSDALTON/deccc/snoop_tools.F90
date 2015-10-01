!> @file
!> SNOOP utilities
!> \author Kasper Kristensen
module snoop_tools_module

  use fundamental
  use precision
  use lstiming
  use ls_util!,only: dgemm_ts
  use TYPEDEF
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

  ! DEC dependencies
  use dec_fragment_utils
  use orbital_operations
  use atomic_fragment_operations

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
  subroutine build_subsystem_lsitem_no_ghost(this,MyMoleculeFULL,lsfull,lssub)
    implicit none
    !> Subsystem label
    integer,intent(in) :: this
    !> Full molecule info
    type(fullmolecule),intent(in) :: MyMoleculeFULL
    !> LSitem for full system
    type(lsitem),intent(in) :: lsfull
    !> lsitem for subsystem
    type(lsitem),intent(inout) :: lssub
    integer,pointer :: subatoms(:)
    integer :: natomssub


    ! Number of atoms for subsystem and list of these atoms
    call get_natoms_for_subsystem(this,MyMoleculeFULL,natomssub)

    call mem_alloc(subatoms,natomssub)
    call list_of_atoms_for_subsystem(this,MyMoleculeFULL,natomssub,subatoms)


    ! Build lssub with atoms for subsystem
    call build_AtomSpecfragmentlsitem(lsfull,lssub,subatoms,natomssub,DECinfo%output,DECinfo%output)
    call mem_dealloc(subatoms)

  end subroutine build_subsystem_lsitem_no_ghost


  !> Get number of number of atoms for subsystem
  !> \author Kasper Kristensen
  !> \date August 2014
  subroutine get_natoms_for_subsystem(this,MyMoleculeFULL,natomssub)
    implicit none
    !> Subsystem label
    integer,intent(in) :: this
    !> Molecule info
    type(fullmolecule),intent(in) :: MyMoleculeFULL
    !> Number of atoms for subsystem "this"
    integer,intent(inout) :: natomssub
    integer :: i

    ! Number of atoms in subsystem "this"
    ! ***********************************
    natomssub=0
    do i=1,MyMoleculeFULL%natoms
       if(MyMoleculeFULL%SubSystemIndex(i)==this .and. (.not. MyMoleculeFUll%PhantomAtom(i))) then
          natomssub = natomssub + 1
       end if
    end do

  end subroutine get_natoms_for_subsystem



  !> Get list of atoms for a given subsystem
  !> \author Kasper Kristensen
  !> \date August 2014
  subroutine list_of_atoms_for_subsystem(this,MyMoleculeFULL,natomssub,subatoms)
    implicit none

    !> Subsystem label
    integer,intent(in) :: this
    !> Molecule info
    type(fullmolecule),intent(in) :: MyMoleculeFULL
    !> Number of atoms in subsystem
    integer,intent(in) :: natomssub
    !> List of atoms in subsystem
    integer,intent(inout) :: subatoms(natomssub)
    integer :: i,idx

    idx=0
    do i=1,MyMoleculeFULL%natoms
       if(MyMoleculeFULL%SubSystemIndex(i)==this .and. (.not. MyMoleculeFUll%PhantomAtom(i)) ) then
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



  !> Get density D = Cocc Cocc^T from occupied orbitals - type(matrix)
  !> \author Kasper Kristensen
  !> \date July 2014
  subroutine get_density_from_occ_orbitals_mat(Cocc,dens)
    implicit none
    !> Occupied MO coefficients 
    type(matrix),intent(in) :: Cocc
    !> Density
    type(matrix),intent(inout) :: dens
    type(matrix) :: Cocc_copy

    ! Cocc copy (avoid passing the same element into dgemm twice)
    call mat_init(Cocc_copy,Cocc%nrow,Cocc%ncol)
    call mat_assign(Cocc_copy,Cocc)

    ! density = Cocc Cocc^T 
    call mat_mul(Cocc, Cocc_copy, 'N', 'T', 1.0_realk, 0.0_realk, dens)
    call mat_free(Cocc_copy)

  end subroutine get_density_from_occ_orbitals_mat


  !> Build lsitem for subsystem with ghost functions on the other subsystems
  !> \author Kasper Kristensen
  !> \date August 2014
  subroutine build_subsystem_lsitem_ghost(this,lsfull,lssub)
    implicit none
    !> Subsystem label
    integer,intent(in) :: this
    !> LSitem for full system
    type(lsitem),intent(in) :: lsfull
    !> lsitem for subsystem
    type(lsitem),intent(inout) :: lssub
    integer :: natoms,i,nelectronsSub
    integer,pointer :: atoms(:)


    ! Strategy
    ! --------
    ! 1. Build lssub to be identical to lsfull
    ! 2. Make modifications such that atoms on other centers are considered as "ghost atoms"

    ! In this way basis function for full system are included, but electron-nuclear 
    ! and nuclear-nuclear interactions involving nuclei outside subsystem "this" are not
    ! considered when Fock matrix and energy are calculated.


    ! 1. Build lssub to be identical to lsfull
    ! ****************************************
    ! List of all atoms
    natoms = lsfull%input%molecule%nAtoms
    call mem_alloc(atoms,natoms)
    do i=1,natoms
       atoms(i)=i
    end do

    ! Build lssub with all atoms included
    call build_AtomSpecfragmentlsitem(lsfull,lssub,atoms,natoms,DECinfo%output,DECinfo%output)


    ! 2. Make modifications such that atoms on other centers are considered as "ghost atoms"
    ! **************************************************************************************


    ! Set number of electrons for subsystem and set atoms on other subsystems to be ghost atoms
    nelectronsSub = 0 
    Do i = 1,natoms
       IF(lssub%input%Molecule%Atom(i)%SubSystemIndex.NE.this)THEN
          lssub%input%Molecule%Atom(i)%Phantom = .TRUE.
       ELSE
          nelectronsSub = nelectronsSub + INT(lssub%input%Molecule%Atom(i)%Charge)
       ENDIF
    Enddo
    lssub%input%Molecule%nelectrons = nelectronsSub

    ! Impose that we only have one subsystem!
    lssub%input%molecule%nSubSystems=1

    call mem_dealloc(atoms)

  end subroutine build_subsystem_lsitem_ghost



  !> Get initial orthogonal basis for subsystem
  !> 1. Take occ and virt MOs from isolated monomer calculations augmented                          
  !>    with zeros for basis functions on the other subsystems                                      
  !> 2. Add virtual orbitals (Cvirtother) from other subsystems augmented with zeros                
  !>   for basis functions on this subsystem                                                       
  !> 3. Orthogonalize Cvirtother againts MOs on this subsystem while keeping occ and virt           
  !>    MOs for this subsystem fixed (they are already orthogonal).                  
  !> \author Kasper Kristensen
  subroutine get_orthogonal_basis_for_subsystem(this,nsub,&
       & MyMoleculeFULL,Cocciso,Cvirtiso,Coccsub,Cvirtsub,S)
    implicit none
    !> Subsystem under consideration
    integer,intent(in) :: this
    !> Number of subsystems
    integer,intent(in) :: nsub
    !> Full molecule info
    type(fullmolecule),intent(in) :: MyMoleculeFULL
    !> MO coefficients for isolated monomer
    type(matrix),intent(in) :: Cocciso(nsub), Cvirtiso(nsub)
    !> Starting guess for occ MO coefficients for subsystem 
    !> expressed in full basis
    !> IMPORTANT: Should be initialized before this routine 
    !             is called!
    type(matrix),intent(inout) :: Coccsub
    !> Starting guess for virt MO coefficients for subsystem 
    !> expressed in full basis
    !> IMPORTANT: Should NOT be initialized before this routine
    !             is called because we do not yet the dimensions.
    !             (We do not know how many redundancies are present
    !              in the virtual MOs for other subsystems,
    !              this is investigated inside this subroutine).
    type(matrix),intent(inout) :: Cvirtsub
    !> AO Overlap matrix
    real(realk),intent(in) :: S(MyMoleculeFULL%nbasis,MyMoleculeFULL%nbasis)
    integer :: nocciso,nvirtiso,nbasisiso,nvirtmax,i,nbasisfull,j,atom
    integer :: idx,nbasissub,sub,nbasisother,nvirtother,norb,nvirtsub
    integer,pointer :: basisidx(:),basisidx_other(:)
    real(realk),pointer :: C(:,:), Cother(:,:)
    real(realk) :: Cother_norm,thr,testnorm(1,1)


    ! Dimensions
    ! ----------
    nbasisfull = MyMoleculeFULL%nbasis  ! nbasis for full system
    nocciso = Cocciso(this)%ncol ! nocc for isolated monomer
    nvirtiso = Cvirtiso(this)%ncol     ! nvirt for isolated monomer
    nbasisiso = Cocciso(this)%nrow     ! nbasis for isolated monomer
    ! Maximum number of virtual orbitals (if there are no redundancies)
    ! - sum of virtual orbitals for all subsystems + possibly ghost atoms
    nvirtmax=0
    do i=1,nsub
       nvirtmax = nvirtmax + Cvirtiso(i)%ncol
    end do
    do atom=1,MyMoleculeFULL%natoms  ! loop over atoms
       if(MyMoleculeFULL%PhantomAtom(atom)) then
          nvirtmax = nvirtmax +MyMoleculeFULL%atom_size(atom)
       end if
    end do


    ! Indices for basis functions in subsystem "this"
    call mem_alloc(basisidx,nbasisiso)
    call basis_idx_subsystem(this,MyMoleculeFULL,nbasisiso,basisidx)    

    ! We discard a virtual orbital from other subsystem if its
    ! norm after projection agaist existing orbitals is below 
    ! this threshold
    thr = 1.0e-5_realk




    ! *************************************************************
    ! Put occ and virt MOs for isolated subsystem into array
    ! with dimensions for the full system, and where we
    ! put zeros on entries corresponding to basis functions
    ! for other subsystems, schematically:
    !        
    ! C =    Cocc(this)   Cvirt(this)   0
    !        0            0             0
    ! *************************************************************
    call mem_alloc(C,nbasisfull,nbasisfull)
    C = 0.0_realk
    ! (C is on purpose allocated to be a little too large here)

    ! Occ orbitals for isolated system "this"
    idx=0
    do j=1,nocciso
       do i=1,nbasisiso
          idx = idx+1
          C(basisidx(i),j) = Cocciso(this)%elms(idx)
       end do
    end do

    ! Virt orbitals for isolated system "this"
    idx=0
    do j=nocciso+1,nvirtiso+nocciso
       do i=1,nbasisiso
          idx = idx+1
          C(basisidx(i),j) = Cvirtiso(this)%elms(idx)
       end do
    end do

    ! Number of orbitals in C array is currently the ones
    ! from subsystem this
    norb = nocciso +nvirtiso
    ! Below we add orbitals from other subsystems


    ! *************************************************************
    ! Put virtual orbitals for other subsystems into C array,
    ! but only after they have been orthogonalized against the
    ! existing orbitals, and only if they are not redundant.
    !        
    ! C =    Cocc(this)   Cvirt(this)   nonzero tails   0
    !        0            0             Cvirt(other)    0
    ! *************************************************************
    call mem_alloc(Cother,nbasisfull,1)
    SubsystemLoop: do sub=1,nsub  ! loop over subsystem

       DifferentSub: if(sub/=this) then 
          ! consider other subsystems than "this"

          ! Number of basis functions for other subsystem
          nbasisother = Cvirtiso(sub)%nrow
          nvirtother = Cvirtiso(sub)%ncol

          ! Basis indices for other subsystem
          call mem_alloc(basisidx_other,nbasisother)
          call basis_idx_subsystem(sub,MyMoleculeFULL,&
               & nbasisother,basisidx_other)

          ! Loop over virtual orbitals on other subsystem
          ! *********************************************
          idx=0
          OtherVirt: do i=1,nvirtother

             ! Virtual orbital "i" on other subsystem
             ! --------------------------------------
             Cother=0.0_realk
             OtherBasis: do j=1,nbasisother
                idx=idx+1
                ! Copy virtual orbital "j" into Cother
                ! where we put zeros on entries corresponding
                ! to other subsystems than "sub"
                Cother(basisidx_other(j),1) = Cvirtiso(sub)%elms(idx)

             end do OtherBasis

             ! Project out components from the norb orbitals already
             ! included in C array
             call project_out_MOs(nbasisfull,norb,C(:,1:norb),&
                  & S,Cother,Cother_norm)

             ! Include orbital "i" only if its norm after
             ! projection is above threshold
             if(Cother_norm > thr) then
                norb = norb + 1

                ! Normalize before putting Cother into array
                do j=1,nbasisfull
                   C(j,norb) = (1.0_realk/Cother_norm)*Cother(j,1)
                end do
                if(DECinfo%SNOOPdebug) then
                   call dec_simple_basis_transform1(nbasisfull,1,&
                        & C(:,norb),S,testnorm)
                   write(DECinfo%output,*) 'sub,norb,norm1,norm2',sub,&
                        & norb,Cother_norm,testnorm(1,1)
                end if

             end if

          end do OtherVirt
          call mem_dealloc(basisidx_other)

       end if DifferentSub

    end do SubsystemLoop



    ! Include possible ghost atoms (phantoms)
    AtomLoop: do atom=1,MyMoleculeFULL%natoms  ! loop over atoms

       Phantom: if(MyMoleculeFULL%PhantomAtom(atom)) then

          ! Number of basis functions for phantom
          nbasisother = MyMoleculeFULL%atom_size(atom)

          ! Basis indices for phantom
          call mem_alloc(basisidx_other,nbasisother)
          call basis_idx_phantom(atom,MyMoleculeFULL,nbasisother,basisidx_other)

          ! Loop over AO basis functions on phantom and include them
          ! ********************************************************
          idx=0
          PhantomVirt: do i=1,nbasisother

             ! Make MO coefficient vector which is just phantom AO number "i"
             ! --------------------------------------------------------------
             Cother=0.0_realk
             Cother(basisidx_other(i),1) = 1.0_realk

             ! Project out components from phantom AO already included in C array
             call project_out_MOs(nbasisfull,norb,C(:,1:norb),&
                  & S,Cother,Cother_norm)

             ! Include orbital "i" only if its norm after
             ! projection is above threshold
             if(Cother_norm > thr) then
                norb = norb + 1

                ! Normalize before putting Cother into array
                do j=1,nbasisfull
                   C(j,norb) = (1.0_realk/Cother_norm)*Cother(j,1)
                end do
                if(DECinfo%SNOOPdebug) then
                   call dec_simple_basis_transform1(nbasisfull,1,&
                        & C(:,norb),S,testnorm)
                   write(DECinfo%output,*) 'PHANTOM,norb,norm1,norm2',atom,&
                        & norb,Cother_norm,testnorm(1,1)
                end if

             end if

          end do PhantomVirt
          call mem_dealloc(basisidx_other)

       end if Phantom

    end do AtomLoop


    call mem_dealloc(Cother)

    ! ***********************************************************
    ! Now C(:,1:norb) contains the orthogonal MOs which can 
    ! be used as starting guess for the calculation on subsystem
    ! "this" using the extended basis.
    ! Finally, put this information into the outputs of type(matrix)
    ! ***********************************************************

    nvirtsub = norb-nocciso
    write(DECinfo%output,'(a,i7,a,2i7)') 'Subsystem: ', this, &
         & '  *** nvirt actual/max ', nvirtsub,nvirtmax

    ! Occupied MOs (see comments in declaration above)
    call mat_set_from_full(C(:,1:nocciso),1E0_realk, Coccsub)

    ! Virtual MOs (see comments in declaration above)
    call mat_init(Cvirtsub,nbasisfull,nvirtsub)
    call mat_set_from_full(C(:,nocciso+1:norb),1E0_realk, Cvirtsub)


    call mem_dealloc(C)
    call mem_dealloc(basisidx)

  end subroutine get_orthogonal_basis_for_subsystem



  !> Get initial orthogonal basis for subsystem be orthogonalizing all virtual orbitals
  !> against the occupied orbitals on subsystem assuming that all MOs
  !> use all are expanded in terms of all basis functions.
  !> \author Kasper Kristensen
  subroutine get_orthogonal_basis_for_subsystem_allvirt(MyMoleculeFULL,Cocc,Cvirt,S)
    implicit none
    !> Full molecule info
    type(fullmolecule),intent(in) :: MyMoleculeFULL
    !> Occupied MO coefficients (not changed)
    type(matrix),intent(in) :: Cocc
    !> Virtual MO coefficients (orthogonalized against Cocc -and possibly remove redundancies)
    type(matrix),intent(inout) :: Cvirt
    !> AO Overlap matrix
    real(realk),intent(in) :: S(MyMoleculeFULL%nbasis,MyMoleculeFULL%nbasis)
    integer :: nocc,nbasis,i,j,nvirtmax,nvirtsub,idx,norb
    real(realk),pointer :: C(:,:), Ctmp(:,:)
    real(realk) :: Ctmp_norm,thr,testnorm(1,1)

    if(Cocc%nrow/=Cvirt%nrow) then
       print *, '#basis function for occ,virt', Cocc%nrow,Cvirt%nrow
       call lsquit('get_orthogonal_basis_for_subsystem_allvirt: Dimension mismatch!',-1)
    end if

    ! Dimensions
    ! ----------
    nocc = Cocc%ncol
    nvirtmax = Cvirt%ncol
    nbasis = Cocc%nrow

    ! We discard a virtual orbital from other subsystem if its
    ! norm after projection agaist existing orbitals is below 
    ! this threshold
    thr = 1.0e-5_realk


    ! Matrix for keep
    call mem_alloc(C,nbasis,nbasis)
    C = 0.0_realk

    ! Copy occ orbitals
    idx=0
    do j=1,nocc
       do i=1,nbasis
          idx = idx+1
          C(i,j) = Cocc%elms(idx)
       end do
    end do
    norb = nocc

    ! *************************************************************
    ! Put virtual orbitals into C array,
    ! but only after they have been orthogonalized against the
    ! existing orbitals, and only if they are not redundant.
    ! *************************************************************

    call mem_alloc(Ctmp,nbasis,1)

    ! Loop over virtual orbitals
    idx=0
    VirtLoop: do i=1,nvirtmax

       ! Virtual orbital "i" on other subsystem
       ! --------------------------------------
       Ctmp=0.0_realk
       do j=1,nbasis
          idx=idx+1
          ! Copy virtual orbital "j" into Ctmp
          ! where we put zeros on entries corresponding
          ! to other subsystems than "sub"
          Ctmp(j,1) = Cvirt%elms(idx)
       end do

       ! Project out components already included in C array
       call project_out_MOs(nbasis,norb,C(:,1:norb),&
            & S,Ctmp,Ctmp_norm)

       ! Include orbital "i" only if its norm after
       ! projection is above threshold
       if(Ctmp_norm > thr) then
          norb = norb + 1

          ! Normalize before putting Ctmp into array
          do j=1,nbasis
             C(j,norb) = (1.0_realk/Ctmp_norm)*Ctmp(j,1)
          end do
          if(DECinfo%SNOOPdebug) then
             call dec_simple_basis_transform1(nbasis,1,&
                  & C(:,norb),S,testnorm)
             write(DECinfo%output,*) 'i,norb,norm1,norm2',i,&
                  & norb,Ctmp_norm,testnorm(1,1)
          end if

       end if

    end do VirtLoop
    call mem_dealloc(Ctmp)


    ! ***********************************************************
    ! Now C(:,1:norb) contains the orthogonal MOs which can 
    ! be used as starting guess for the calculation on subsystem
    ! Finally, put this information into the outputs of type(matrix)
    ! ***********************************************************

    nvirtsub = norb-nocc
    write(DECinfo%output,'(a,2i7)') 'Subsystem: nvirt actual/max ', nvirtsub,nvirtmax

    ! Virtual MOs (see comments in declaration above)
    call mat_free(Cvirt)
    call mat_init(Cvirt,nbasis,nvirtsub)
    call mat_set_from_full(C(:,nocc+1:norb),1E0_realk, Cvirt)


    call mem_dealloc(C)

  end subroutine get_orthogonal_basis_for_subsystem_allvirt




  !> Get basis indices for the basis function
  !> associated with a given subsystem
  subroutine basis_idx_subsystem(this,MyMoleculeFULL,nbasissub,basisidx)
    implicit none
    !> Subsystem index
    integer,intent(in) :: this
    !> Full molecule info
    type(fullmolecule),intent(in) :: MyMoleculeFULL
    !> Number of basis functions for subsystem
    integer,intent(in) :: nbasissub
    !> Basis function indices for subsystem basis functions
    integer,intent(inout) :: basisidx(nbasissub)
    integer :: atomsub,idx,i,j

    basisidx=0

    FullDimension: if(nbasissub==MyMoleculeFULL%nbasis) then

       ! Basis functions for full molecule are all used for subsystem (ghost functions are used),
       ! so basisidx is simply 1,2,3,...,nbasis
       do i=1,nbasissub
          basisidx(i) = i
       end do

    else
       ! Only basis functions on atoms in subsystem are used - extract these basis function indices

       idx=0
       do i=1,MyMoleculeFULL%natoms
          ! Subsystem to which atom "i" belongs
          atomsub = MyMoleculeFULL%SubSystemIndex(i)

          if(atomsub == this .and. (.not. MyMoleculeFULL%PhantomAtom(i)) ) then 
             ! Atom "i" is in subsystem 
             ! --> include basis function indices for atom "i"
             do j=1,MyMoleculeFULL%atom_size(i)
                idx=idx+1
                basisidx(idx) = MyMoleculeFULL%atom_start(i) + j-1
             end do
          end if

       end do

       ! Sanity check
       if(idx/=nbasissub) then
          print *, 'idx, nbasissub', idx, nbasissub
          call lsquit('basis_idx_subsystem: Counter mismatch',-1)
       end if

    end if FullDimension


  end subroutine basis_idx_subsystem



  !> Get basis indices for the basis function associated with a phantom atom
  subroutine basis_idx_phantom(this,MyMoleculeFULL,nbasisphantom,basisidx)
    implicit none
    !> Phantom atom index
    integer,intent(in) :: this
    !> Full molecule info
    type(fullmolecule),intent(in) :: MyMoleculeFULL
    !> Number of basis functions for subsystem
    integer,intent(in) :: nbasisphantom
    !> Basis function indices for subsystem basis functions
    integer,intent(inout) :: basisidx(nbasisphantom)
    integer :: j,idx

    basisidx=0
    idx=0
    do j=1,MyMoleculeFULL%atom_size(this)
       idx=idx+1
       basisidx(idx) = MyMoleculeFULL%atom_start(this) + j-1
    end do

    ! Sanity check
    if(idx/=nbasisphantom) then
       print *, 'idx, nbasisphantom',idx, nbasisphantom
       call lsquit('basis_idx_phantom: Index mismatch',-1)
    end if

  end subroutine basis_idx_phantom



  ! Project out components of a molecular orbital Cb
  ! defined by an input matrix of other MOs Ca:
  !
  ! Cb --> (1 - Ca Ca^T S ) Cb
  ! 
  ! Cb is one column vector, while Ca is in general a matrix
  ! representing a set of MOs, and S is the overlap matrix
  ! in the AO basis.
  ! Also calculate norm of resulting Cb.
  subroutine project_out_MOs(nbasis,na,Ca,S,Cb,Cb_norm)
    implicit none
    !> Number of basis functions
    integer,intent(in) :: nbasis
    !> Number of orbitals in Ca matrix
    integer,intent(in) :: na
    !> Ca and S as defined above
    real(realk),intent(in) :: Ca(nbasis,na), S(nbasis,nbasis)
    !> Output Cb vector
    real(realk),intent(inout) :: Cb(nbasis,1)
    !> Norm of projected Cb (before renormalization)
    real(realk),intent(inout) :: Cb_norm
    real(realk),pointer :: tmp(:,:),tmp2(:,:),tmp3(:)
    integer :: one,i
    real(realk) :: norm_squared(1,1)

    ! tmp = Ca^T S
    call mem_alloc(tmp,na,nbasis)
    call dec_simple_dgemm(na,nbasis,nbasis,Ca,S,tmp,'T','N')

    ! tmp2 = Ca Ca^T S = Ca tmp
    call mem_alloc(tmp2,nbasis,nbasis)
    call dec_simple_dgemm(nbasis,na,nbasis,Ca,tmp,tmp2,'N','N')
    call mem_dealloc(tmp)
    
    ! tmp3 = Ca Ca^T S Cb = tmp2 Cb
    one = 1
    call mem_alloc(tmp3,nbasis)
    call dec_simple_dgemm(nbasis,nbasis,one,tmp2,Cb,tmp3,'N','N')
    call mem_dealloc(tmp2)
    
    ! Output Cb = Cb - Ca Ca^T S Cb = Cb - tmp3
    do i=1,nbasis
       Cb(i,1) = Cb(i,1) - tmp3(i)
    end do
    call mem_dealloc(tmp3)
    

    ! Norm of final Cb
    call dec_simple_basis_transform1(nbasis,one,Cb,S,norm_squared)
    Cb_norm = sqrt(norm_squared(1,1))

  end subroutine project_out_MOs


  !> Find initial MOs for subsystem calculation from (localized) full MOs.
  subroutine initial_subsystem_MOs_from_full(MyMoleculeFULL,mylsitem,nsub,Coccsub,Cvirtall)
    implicit none
    !> Full molecule info
    type(fullmolecule),intent(in) :: MyMoleculeFULL
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> Number of subsystems
    integer,intent(in) :: nsub
    !> Occupied MOs for each subsystem by simply copying localized full orbitals
    !> (matrices are also initialized in here because we do not know their dimensions beforehand)
    type(matrix),dimension(nsub),intent(inout) :: Coccsub 
    !> All virtual orbitals put into matrix structure (initialized BEFORE calling this)
    type(matrix),intent(inout) :: Cvirtall
    real(realk),pointer :: Cocc(:,:),Cvirt(:,:)
    integer :: i,this,j,noccsub,nvirtsub,atom,idx
    type(decorbital),pointer :: OccOrbitals(:), VirtOrbitals(:)

    if(MyMoleculeFULL%mem_distributed)then
       call lsquit("ERROR(initial_subsystem_MOs_from_full) not implemented for pdm",-1)
    endif

    ! Determine DEC orbital structures for local orbital analysis below
    call mem_alloc(OccOrbitals,MyMoleculeFULL%nocc)
    call mem_alloc(VirtOrbitals,MyMoleculeFULL%nvirt)
    call GenerateOrbitals_driver(MyMoleculeFULL,mylsitem,MyMoleculeFULL%nocc,MyMoleculeFULL%nvirt,&
         & MyMoleculeFULL%natoms, OccOrbitals, VirtOrbitals)

    ! Allocate temporary arrays to be to big to make sure we have enough space
    call mem_alloc(Cocc,MyMoleculeFULL%nbasis,MyMoleculeFULL%nocc)


    ! Occupied orbitals
    ! =================
    do i=1,nsub

       ! Subsystem under consideration
       this = i

       ! Number of occupied orbitals for subsystem "this"
       idx=0
       do j=1,MyMoleculeFULL%nocc
          atom = OccOrbitals(j)%centralatom
          ! Does atom to which orbital "j" is assigned belong to subsystem "this"
          if(MyMoleculeFULL%SubSystemIndex(atom)==this) then  
             idx = idx + 1
             Cocc(:,idx) = MyMoleculeFULL%Co%elm2(:,j)
          end if
       end do
       noccsub=idx

       ! Set output matrix
       call mat_init(Coccsub(this),MyMoleculeFULL%nbasis,noccsub)
       call mat_set_from_full(Cocc(:,1:noccsub),1E0_realk, Coccsub(this))

    end do


    ! Virtual orbitals
    ! ================
    ! Simply copy from molecule structure
    call mat_set_from_full(MyMoleculeFULL%Cv%elm2,1E0_realk, Cvirtall)

  
    call mem_dealloc(Cocc)

    do i=1,MyMoleculeFULL%nocc
       call orbital_free(OccOrbitals(i))
    end do
    do i=1,MyMoleculeFULL%nvirt
       call orbital_free(VirtOrbitals(i))
    end do
    call mem_dealloc(OccOrbitals)
    call mem_dealloc(VirtOrbitals)

  end subroutine initial_subsystem_MOs_from_full



  !> \brief Carry out unitary rotation of occupied as well as virtual subsystem orbitals
  !> such that they mimic the FULL occupied and virtual orbitals as much as possible
  !> in a least squares sense (defined by the natural connection).
  subroutine rotate_subsystem_orbitals_to_mimic_FULL_orbitals(MyMoleculeFULL,sub,&
       & OccOrbitals,VirtOrbitals,lssub,CoccSUBmat,CvirtSUBmat,S)
    implicit none
    !> Full molecule info for FULL
    type(fullmolecule),intent(in) :: MyMoleculeFULL
    !> Which subsystem
    integer,intent(in) :: sub
    !> Occupied and virtual orbitals in DEC format
    type(decorbital),intent(in) :: OccOrbitals(MyMoleculeFULL%nvirt), &
         & VirtOrbitals(MyMoleculeFULL%nvirt)
    !> LSitem for subsystem 
    type(lsitem), intent(inout) :: lssub
    !> Occupied and virtual orbitals for subsystem to be rotated
    type(matrix),intent(inout) :: CoccSUBmat, CvirtSUBmat
    !> AO Overlap matrix
    real(realk),intent(in) :: S(MyMoleculeFULL%nbasis,MyMoleculeFULL%nbasis)
    real(realk),pointer :: CcoreSUB(:,:),CvalSUB(:,:),CvirtSUB(:,:),CoccSUB(:,:)
    real(realk),pointer :: CcoreFULL(:,:),CvalFULL(:,:),CvirtFULL(:,:)
    integer :: noccSUB,nbasis,nvirt,ncoreFULL,ncoreSUB,nvalSUB,nvalFULL,noccFULL
    integer :: i

    if(MyMoleculeFULL%mem_distributed)then
       call lsquit("ERROR(rotate_subsystem_orbitals_to_mimic_FULL_orbitals) not implemented for pdm",-1)
    endif

    nvirt = MyMoleculeFULL%nvirt  ! #virt orbitals same for SNOOP subsystem and FULL
    nbasis = MyMoleculeFULL%nbasis
    noccFULL = MyMoleculeFULL%nocc    ! Occupied full molecule
    ncoreFULL = MyMoleculeFULL%ncore
    nvalFULL = noccFULL - ncoreFULL
    noccSUB = CoccSUBmat%ncol       ! # occupied subsystem orbitals
    ncoreSUB = count_ncore(lssub)
    nvalSUB = noccSUB-ncoreSUB

    ! For simplicity, we consider everything as fortran arrays now

    ! Initial occupied and virtual subsystem orbitals
    ! -----------------------------------------------
    call mem_alloc(CoccSUB,nbasis,noccSUB)
    call mat_to_full(CoccSUBmat, 1.0_realk, CoccSUB)  ! Occupied subsystem orbitals
    call mem_alloc(CcoreSUB,nbasis,ncoreSUB)
    call mem_alloc(CvalSUB,nbasis,nvalSUB)
    ! Partition into core and valence orbitals
    do i=1,ncoreSUB
       CcoreSUB(:,i) = CoccSUB(:,i)
    end do
    do i=1,nvalSUB
       CvalSUB(:,i) = CoccSUB(:,i+ncoreSUB)
    end do

    ! Virtual orbitals
    call mem_alloc(CvirtSUB,nbasis,nvirt)
    call mat_to_full(CvirtSUBmat, 1.0_realk, CvirtSUB)



    ! FULL occupied and virtual orbitals relevant for subsystem
    ! ----------------------------------------------------------
    ! Extract occupied (core+valence) orbitals assigned to subsystem
    call mem_alloc(CcoreFULL,nbasis,ncoreSUB)
    call mem_alloc(CvalFULL,nbasis,nvalSUB)
    call extract_subsystem_occorbitals_from_FULL(MyMoleculeFULL,OccOrbitals,&
         & sub,ncoreSUB,nvalSUB,CcoreFULL,CvalFULL)

    ! Simply copy ALL virtual FULL orbitals because virtual dimensions are
    ! the same for subsystem and FULL
    call mem_alloc(CvirtFULL,nbasis,nvirt)
    CvirtFULL = MyMoleculeFULL%Cv%elm2

    
    ! Rotate core subsystem orbitals to mimic FULL orbitals
    call get_natural_connection_subsystem_matrix(nbasis,ncoreSUB,S,&
         & CcoreFULL,CcoreSUB)

    ! Rotate valence subsystem orbitals to mimic FULL orbitals
    call get_natural_connection_subsystem_matrix(nbasis,nvalSUB,S,&
         & CvalFULL,CvalSUB)

    ! Rotate virtual subsystem orbitals to mimic FULL orbitals
    call get_natural_connection_subsystem_matrix(nbasis,nvirt,S,&
         & CvirtFULL,CvirtSUB)

    ! Collect occupied orbitals in one matrix with core before valence
    do i=1,ncoreSUB
       CoccSUB(:,i) = CcoreSUB(:,i)
    end do
    do i=1,nvalSUB
       CoccSUB(:,i+ncoreSUB) = CvalSUB(:,i)
    end do

    ! Set output in type(matrix) form
    call mat_set_from_full(CoccSUB,1E0_realk,CoccSUBmat)
    call mat_set_from_full(CvirtSUB,1E0_realk,CvirtSUBmat)


    call mem_dealloc(CoccSUB)
    call mem_dealloc(CcoreSUB)
    call mem_dealloc(CvalSUB)
    call mem_dealloc(CvirtSUB)
    call mem_dealloc(CcoreFULL)
    call mem_dealloc(CvalFULL)
    call mem_dealloc(CvirtFULL)

  end subroutine rotate_subsystem_orbitals_to_mimic_FULL_orbitals



  !> \brief Extract occupied (partitioned into core and valence) orbitals 
  !> assigned to a given subsystem from orbitals for full system stored in MyMoleculeFULL.
  subroutine extract_subsystem_occorbitals_from_FULL(MyMoleculeFULL,OccOrbitals,&
       & sub,ncoreSUB,nvalSUB,Ccore,Cval)
    implicit none
    !> Full molecule info for FULL
    type(fullmolecule),intent(in) :: MyMoleculeFULL
    !> Occupied  orbitals in DEC format
    type(decorbital),intent(in) :: OccOrbitals(MyMoleculeFULL%nvirt)
    !> Which subsystem
    integer,intent(in) :: sub
    !> Number of core and valence orbitals for subsystem
    integer,intent(in) :: ncoreSUB,nvalSUB
    !> Extracted core and valence MO coefficients corresponding to orbitals
    !> assigned to subsystem
    real(realk) :: Ccore(MyMoleculeFULL%nbasis,ncoreSUB), &
         & Cval(MyMoleculeFULL%nbasis,nvalSUB)
    integer :: j,idx,atom


    if(MyMoleculeFULL%mem_distributed)then
       call lsquit("ERROR(extract_subsystem_occorbitals_from_FULL) not implemented for pdm",-1)
    endif

    ! Extract core orbitals
    ! *********************
    idx=0
    do j=1,MyMoleculeFULL%ncore
       atom = OccOrbitals(j)%centralatom
       ! Does atom to which orbital "j" belong to subsystem "sub"
       if(MyMoleculeFULL%SubSystemIndex(atom)==sub) then  
          idx = idx + 1
          if(idx>ncoreSUB) then
             print *, 'Core error: idx,nMO',idx,ncoreSUB
             call lsquit('extract_subsystem_occorbitals_from_FULL: &
                  & Core idx is too large',-1)
          end if
          Ccore(:,idx) = MyMoleculeFULL%Co%elm2(:,j)
       end if
    end do

    if(idx/=ncoreSUB) then
       print *, 'Core: idx,ncoreSUB',idx,ncoreSUB
       call lsquit('extract_subsystem_occorbitals_from_FULL: core idx error!',-1)
    end if



    ! Extract valence orbitals
    ! ************************
    idx=0
    do j=MyMoleculeFULL%ncore+1,MyMoleculeFULL%nocc
       atom = OccOrbitals(j)%centralatom
       ! Does atom to which orbital "j" belong to subsystem "sub"
       if(MyMoleculeFULL%SubSystemIndex(atom)==sub) then  
          idx = idx + 1
          if(idx>nvalSUB) then
             print *, 'Valence error: idx,nvalSUB',idx,nvalSUB
             call lsquit('extract_subsystem_occorbitals_from_FULL: &
                  & Valence idx is too large',-1)
          end if
          Cval(:,idx) = MyMoleculeFULL%Co%elm2(:,j)
       end if
    end do

    if(idx/=nvalSUB) then
       print *, 'Valence: idx,nvalSUB',idx,nvalSUB
       call lsquit('extract_subsystem_occorbitals_from_FULL: valence idx error!',-1)
    end if


  end subroutine extract_subsystem_occorbitals_from_FULL


  !> Get natural connection matrix which rotates subsystem orbitals to
  !> mimic FULL orbitals assigned to that subsystem as much as possible.
  !> (Theor Chim Acta 90,421 (1995) applied to this problem).
  subroutine get_natural_connection_subsystem_matrix(nbasis,nMO,S,CFULL,Csub)
    implicit none
    !> Number of atomic basis functions
    integer,intent(in) :: nbasis
    !> Number of molecular orbitals (either nocc or nvirt)
    integer,intent(in) :: nMO
    !> AO overlap matrix
    real(realk),intent(in) :: S(nbasis,nbasis)
    !> FULL orbitals
    real(realk),intent(in) :: CFULL(nbasis,nMO)
    !> Subsystem orbitals
    real(realk),intent(inout) :: Csub(nbasis,nMO)
    real(realk),pointer :: W(:,:),T(:,:),tmp(:,:)
    real(realk) :: funcval1,funcval2   ! KKFIXME delete this when properly tested


    ! Strategy
    ! ********
    !
    ! We want to rotate the subsystem orbitals such that they resemble the
    ! FULL orbitals as much as possible in a least-squares-sense.
    ! This is done using the natural connection (Theor Chim Acta 90,421 (1995)).
    ! 
    ! Effectively we replace the subsystem orbitals:
    !
    ! Csub --> Csub T
    ! 
    ! where
    ! 
    ! T = W^-1 (W W^T)^{1/2}
    ! W = CFULL^T S Csub

    ! Calculate W
    call mem_alloc(W,nMO,nMO)
    call dec_diff_basis_transform1(nbasis,nMO,nMO,CFULL,Csub,S,W)

    ! Get T matrix: T = W^-1 (W W^T)^{1/2} 
    call mem_alloc(T,nMO,nMO)
    call get_natural_connection_T_matrix(nMO,W,T)
    ! Value of difference measure function (m
    funcval1 = snoop_natural_orbitals_diff_measure(nbasis,nMO,S,CFULL,W) 

    ! Csub --> Csub T 
    call mem_alloc(tmp,nbasis,nMO)
    tmp = Csub
    call dec_simple_dgemm(nbasis,nMO,nMO,tmp,T,Csub,'n','n')
    call mem_dealloc(tmp)

    ! Value of difference measure function after changing Csub
    call dec_diff_basis_transform1(nbasis,nMO,nMO,CFULL,Csub,S,W)
    funcval2 = snoop_natural_orbitals_diff_measure(nbasis,nMO,S,CFULL,W)
    write(DECinfo%output,'(1X,a,i7,3g15.5)') 'SNOOP funcval ',&
         & nMO,funcval1,funcval2,funcval2-funcval1

    call mem_dealloc(W)
    call mem_dealloc(T)

  end subroutine get_natural_connection_subsystem_matrix


  
  !> Get value of function which is minimized by choosing natural connection
  !> for determining unitary transformation in get_natural_connection_subsystem_matrix.
  !> (Mainly for debugging).
  function snoop_natural_orbitals_diff_measure(nbasis,nMO,S,CFULL,W) result(funcval)
    implicit none
    !> Number of atomic basis functions
    integer,intent(in) :: nbasis
    !> Number of molecular orbitals (either nocc or nvirt)
    integer,intent(in) :: nMO
    !> AO overlap matrix
    real(realk),intent(in) :: S(nbasis,nbasis)
    !> FULL orbitals
    real(realk),intent(in) :: CFULL(nbasis,nMO)
    !> W matrix defined in get_natural_connection_subsystem_matrix
    real(realk),dimension(nMO,nMO),intent(in) :: W
    real(realk) :: funcval
    real(realk),pointer :: SMO(:,:)
    integer :: i

    !> FULL MO overlap (has to be unit matrix but calculate it here to be general)
    call mem_alloc(SMO,nMO,nMO)
    call dec_simple_basis_transform1(nbasis,nMO,CFULL,S,SMO)

    ! Calculate function value:
    ! Eq. 9 in Theor Chim Acta 90,421 (1995) for T=1 since this is point where we are standing
    funcval = 0.0_realk
    do i=1,nMO
       funcval = funcval + SMO(i,i) +1.0_realk - 2.0_realk*W(i,i)
    end do

    call mem_dealloc(SMO)

  end function snoop_natural_orbitals_diff_measure


  !> Calculate natural connection T matrix, Eq. 21 in Theor Chim Acta 90,421 (1995):
  !> T = W^-1 (W W^T)^{1/2} 
  subroutine get_natural_connection_T_matrix(nMO,W,T)
    implicit none
    !> Number of MOs
    integer,intent(in) :: nMO
    !> W matrix as defined in get_natural_connection_subsystem_matrix
    real(realk),dimension(nMO,nMO),intent(in) :: W
    !> T = W^-1 (W W^T)^{1/2}
    real(realk),dimension(nMO,nMO),intent(inout) :: T
    real(realk),pointer :: Winv(:,:),tmp(:,:),tmp2(:,:)

    ! Invert W
    call mem_alloc(Winv,nMO,nMO)
    call invert_matrix(W,Winv,nMO)

    ! (W W^T)^{1/2} 
    ! -------------------

    ! tmp2 = W W^T
    call mem_alloc(tmp,nMO,nMO)
    call mem_alloc(tmp2,nMO,nMO)
    tmp = W
    call dec_simple_dgemm(nMO,nMO,nMO,W,tmp,tmp2,'n','t')

    ! tmp = (W W^T)^{1/2}
    call get_power_of_symmetric_matrix(nMO,0.5_realk,tmp2,tmp)

    ! T = W^-1 (W W^T)^{1/2} 
    call dec_simple_dgemm(nMO,nMO,nMO,Winv,tmp,T,'n','n')

    call mem_dealloc(tmp)
    call mem_dealloc(tmp2)
    call mem_dealloc(Winv)


  end subroutine get_natural_connection_T_matrix


  !> Initialize atomic fragments for subsystem to use FULL orbital spaces
  !> using the one-one mapping defined by the natural connection 
  !> (see rotate_subsystem_orbitals_to_mimic_FULL_orbitals).
  subroutine subsystemAOS_equals_FULLAOS(sub,MyMoleculeFULL,OccOrbitalsFULL,&
       & MySubsystem,lssub,OccOrbitalsSUB,&
       & VirtOrbitalsSUB,dofragFULL,dofragsub,AFfull,AFsub)
    implicit none
    !> Which subsystem
    integer,intent(in) :: sub
    !> Molecule structure for full system
    type(fullmolecule),intent(in) :: MyMoleculeFULL
    !> Occupied orbitals for full system in DEC format
    type(decorbital),intent(in) :: OccOrbitalsFULL(MyMoleculeFULL%nocc)
    !> Molecule structure for subsystem
    type(fullmolecule),intent(in) :: MySubsystem
    !> LSitem for subsystem 
    type(lsitem), intent(inout) :: lssub
    !> Occ and virt orbitals for subsystem in DEC format
    type(decorbital),intent(in) :: OccOrbitalsSUB(MySubsystem%nocc), VirtOrbitalsSUB(MySubsystem%nvirt)
    !> Which atoms to consider for FULL
    logical,intent(in) :: dofragFULL(MyMoleculeFULL%natoms)
    !> Which atoms to consider for subsystem
    logical,intent(in) :: dofragsub(MyMoleculeFULL%natoms)
    !>  Atomic fragments for FULL
    type(decfrag),intent(in) :: AFfull(MyMoleculeFULL%natoms)
    !> Atomic fragments for subsystem
    type(decfrag),intent(inout) :: AFsub(MyMoleculeFULL%natoms)
    integer :: i,j,noccAOSsub,fullidx,counter,subidx,natoms
    integer :: fulltosub(MyMoleculeFULL%nocc)
    integer,pointer :: occAOS(:)
    logical :: DoBasis,pairfrag

    ! Number of atoms (KK fixme: Needs to be modified for DECCO)
    natoms = MyMoleculeFULL%natoms

    ! For occupied orbitals, get list which define the one-to-one correspondance from
    ! full indices to subsystem indices.
    call get_fulloccAOS_to_suboccAOS(sub,MyMoleculeFULL,OccOrbitalsFULL,fulltosub)
    ! Note: Since the number of virtual orbitals for full molecule and subsystem
    ! is the same, the corresponding list for virtual orbitals would just be a 
    ! trivial 1,2,3...,nvirt list so we do not construct it.

    do i=1,natoms
       SubsystemFragment: if(dofragsub(i)) then ! atomic fragment considered for subsystem
          ! Sanity check: If atom is not considered for FULL we are in trouble...
          if(.not. dofragFULL(i)) then
             print *, 'Atom = ',i
             print *, 'dofragsub  ', dofragsub
             print *, 'dofragFULL', dofragFULL
             call lsquit('subsystemAOS_equals_FULLAOS: Atom considered for subsystem but not FULL!',-1)
          end if

          ! Set indices for occupied AOS
          ! ----------------------------
          ! Here care must be exercised because the number of occupied orbitals
          ! is different for subsystem and full system. We therefore use the fulltosub
          ! conversion to translate the full occupied AOS index into the
          ! corresponding subsystem occupied AOS index.
          ! However, if the full occupied AOS includes orbitals which are not assigned 
          ! to the subsystem in question, we do not include them. 
          ! (Full occupied orbitals assigned to subsystem B simply do not have a counterpart
          !  in the subsystem A calculation etc.).  (*)
          ! 
          counter=0
          call mem_alloc(occAOS,AFfull(i)%noccAOS)  ! Subsystem occ AOS
          do j=1,AFfull(i)%noccAOS  ! loop over FULL occ AOS
             fullidx = AFfull(i)%occAOSidx(j)   ! Index in full list
             subidx = fulltosub(fullidx)         ! corresponding subsystem index

             ! if subidx/=0 --> include occ AOS orbital because it is assigned to subsystem
             if(subidx/=0) then
                counter=counter+1
                occAOS(counter) = subidx
             end if
          end do
          noccAOSsub=counter   ! occAOS dimension for subsystem

          ! Initialize atomic fragment for subsystem with "same" orbitals as in FULL
          ! calculation in the sense defined in the subroutine 
          ! rotate_subsystem_orbitals_to_mimic_FULL_orbitals 
          ! (although comment (*) above means that there might be a slight difference
          !  in occupied AOS).
          DoBasis=.false.
          pairfrag=.false.
          call atomic_fragment_init_integer_list(i,MySubsystem%nvirt, MySubsystem%nocc, &
               & AFfull(i)%nvirtAOS,noccAOSsub,&
               & AFfull(i)%virtAOSidx,occAOS(1:noccAOSsub),OccOrbitalsSUB,VirtOrbitalsSUB,&
               & MySubsystem,lssub,AFsub(i),DoBasis,pairfrag)
          call mem_dealloc(occAOS)

       end if SubsystemFragment

    end do

    ! Sanity check that orbital spaces have the correct sizes
    call SUBatomic_fragments_sanity_check(sub,natoms,AFfull,AFsub,dofragsub)

  end subroutine subsystemAOS_equals_FULLAOS


  !> For each occupied orbital in full molecule, get corresponding index in subsystem.
  !> If the full occupied index is not present in the subsystem, the corresponding index
  !> is set to zero.
  subroutine get_fulloccAOS_to_suboccAOS(sub,MyMoleculeFULL,OccOrbitals,fulltosub)
    implicit none
    !> Which subsystem
    integer,intent(in) :: sub
    !> Molecule structure for full molecule
    type(fullmolecule),intent(in) :: MyMoleculeFULL
    !> Occ orbitals for full molecule in DEC format
    type(decorbital),intent(in) :: OccOrbitals(MyMoleculeFULL%nocc)
    !> Index conversion from full system to subsystem
    integer,intent(inout) :: fulltosub(MyMoleculeFULL%nocc)
    integer :: j,idx,atom

    ! Zero for all indices corresponding to other subsystems
    fulltosub=0

    idx=0
    do j=1,MyMoleculeFULL%nocc
       atom = OccOrbitals(j)%centralatom
       ! Does atom to which orbital "j" is assigned belong to subsystem "sub"
       if(MyMoleculeFULL%SubSystemIndex(atom)==sub) then
          idx = idx + 1
          ! Orbital "j" in full system corresponds to orbital "idx" in subsystem
          fulltosub(j)=idx
       end if
    end do

  end subroutine get_fulloccAOS_to_suboccAOS


  !> Check that central atom for all subsystem orbitals are the same as for full system.
  !> If necessary, reassign orbitals to ensure that this is the case.
  subroutine Orbitals_subsystem_vs_FULL_sanity_check(sub,MyMoleculeFULL,MySubsystem,&
       & OccOrbitalsFULL,VirtOrbitalsFULL,OccOrbitalsSUB,VirtOrbitalsSUB)
    implicit none
    !> Which subsystem
    integer,intent(in) :: sub
    !> Molecule structure for full system
    type(fullmolecule),intent(in) :: MyMoleculeFULL
    !> Molecule structure for subsystem
    type(fullmolecule),intent(in) :: MySubsystem
    !> Occupied orbitals for full system in DEC format
    type(decorbital),intent(in) :: OccOrbitalsFULL(MyMoleculeFULL%nocc)
    !> Virtual orbitals for full system in DEC format
    type(decorbital),intent(in) :: VirtOrbitalsFULL(MyMoleculeFULL%nvirt)
    !> Occ and virt orbitals for subsystem in DEC format
    type(decorbital),intent(inout) :: OccOrbitalsSUB(MySubsystem%nocc), VirtOrbitalsSUB(MySubsystem%nvirt)
    integer :: fulltosub(MyMoleculeFULL%nocc)
    integer :: i,subidx,subatom,fullatom
    logical :: dofragSUB1(MySubsystem%nfrags), dofragSUB2(MySubsystem%nfrags)

    !  List of which fragments to consider for subsystem (for extra check)
    call which_fragments_to_consider(MySubsystem%ncore,MySubsystem%nocc,MySubsystem%nvirt,&
         & MySubsystem%nfrags,OccOrbitalsSUB,VirtOrbitalsSUB,&
         & dofragSUB1,MySubsystem%PhantomAtom)

    ! Get list taking us FROM full occupied index TO subsystem occupied index
    call get_fulloccAOS_to_suboccAOS(sub,MyMoleculeFULL,OccOrbitalsFULL,fulltosub)

    do i=1,MyMoleculeFULL%nocc
       subidx = fulltosub(i)

       OrbitalInSubsystem: if(subidx/=0) then  ! "orbital i is included for subsystem"

          ! Atom to which orbital "i" is assigned for full system and subsystem
          fullatom = OccOrbitalsFULL(i)%centralatom
          subatom = OccOrbitalsSUB(subidx)%centralatom

          if(fullatom /= subatom) then
             ! This occupied orbital is assigned differently for subsystem,
             ! reassign such that it is assigned to same atom as for full system.
             OccOrbitalsSUB(subidx)%centralatom = fullatom
          end if

       end if OrbitalInSubsystem

    end do


    ! Same for virtual - except now we do not need a "fulltosub" list because
    ! there is a one-to-one correspondence between the orbitals.
    ! (The number of virtual orbitals is the same for full system and subsystem).
    do i=1,MyMoleculeFULL%nvirt

       ! Atom to which orbital "i" is assigned for full system and subsystem
       fullatom = VirtOrbitalsFULL(i)%centralatom
       subatom = VirtOrbitalsSUB(i)%centralatom

       if(fullatom /= subatom) then
          VirtOrbitalsSUB(i)%centralatom = fullatom
       end if

    end do


    !  List of which fragments to consider for subsystem after reassigning
    call which_fragments_to_consider(MySubsystem%ncore,MySubsystem%nocc,MySubsystem%nvirt,&
         & MySubsystem%nfrags,OccOrbitalsSUB,VirtOrbitalsSUB,&
         & dofragSUB2,MySubsystem%PhantomAtom)
    
    ! Check that reassign did not modify list of which fragments to consider
    do i=1,MySubsystem%nfrags
       if(dofragSUB1(i) .neqv. dofragSUB2(i)) then
          print *, 'Fragment ',i
          print *, 'dofragSUB1 ', dofragSUB1(i)
          print *, 'dofragSUB2 ', dofragSUB2(i)
          call lsquit('Orbitals_subsystem_vs_FULL_sanity_check: Final check failed!')
       end if
    end do

  end subroutine Orbitals_subsystem_vs_FULL_sanity_check



  !> Check that occupied and virtual EOS as well as virtual AOS are the same
  !> for subsystem and full system.
  !> (Occupied AOS is not checked because we cannot guarantee that they are identical,
  !   see subsystemAOS_equals_FULLAOS).
  subroutine SUBatomic_fragments_sanity_check(sub,natoms,AFfull,AFsub,dofragsub)
    implicit none
    !> Which subsystem?
    integer,intent(in) :: sub
    !> Number of atoms in FULL molecule
    integer,intent(in) :: natoms
    !>  Atomic fragments for FULL
    type(decfrag),intent(in) :: AFfull(natoms)
    !> Atomic fragments for subsystem
    type(decfrag),intent(in) :: AFsub(natoms)
    !> Which atoms to consider for subsystem
    logical,intent(in) :: dofragsub(natoms)
    integer :: atom
    logical :: something_wrong

    do atom=1,natoms
       if(dofragsub(atom)) then

          something_wrong=.false.

          ! Occupied EOS
          if(AFfull(atom)%noccEOS /= AFsub(atom)%noccEOS) then
             something_wrong=.true.
          end if

          ! Virtual EOS
          if(AFfull(atom)%nvirtEOS /= AFsub(atom)%nvirtEOS) then
             something_wrong=.true.
          end if

          ! Virtual AOS
          if(AFfull(atom)%nvirtAOS /= AFsub(atom)%nvirtAOS) then
             something_wrong=.true.
          end if

          if(something_wrong) then
             print *, 'Subsystem: ',sub
             print *, 'Atom: ', atom
             print *, 'Occ  EOS: Full/sub', AFfull(atom)%noccEOS,AFsub(atom)%noccEOS
             print *, 'Virt EOS: Full/sub', AFfull(atom)%nvirtEOS,AFsub(atom)%nvirtEOS
             print *, 'Virt AOS: Full/sub', AFfull(atom)%nvirtAOS,AFsub(atom)%nvirtAOS
             call lsquit('SUBatomic_fragments_sanity_check: Orbital mismatch!',-1)
          end if

       end if
    end do

  end subroutine SUBatomic_fragments_sanity_check


  !> \brief Write SNOOP restart file for full system (sub=0) or
  !> one of the subsystems
  !> \author Kasper Kristensen
  !> \date August 2015
  subroutine snoop_write_restart_file(sub,EHF,Ecorr)
    implicit none
    !> Which subsystem? (If sub=0, then we consider the full system)
    integer,intent(in) :: sub 
    !> HF and correlation energy 
    real(realk),intent(in) :: EHF,Ecorr
    integer :: funit
    character(len=17) :: FileName    

    ! Get file name
    call SNOOP_get_filename(sub,FileName)

    ! Open file (always delete existing file if it exists)
    funit = -1
    call lsopen(funit,FileName,'REPLACE','UNFORMATTED')

    ! Write energies to file and close file
    write(funit) EHF
    write(funit) Ecorr
    call lsclose(funit,'KEEP')

  end subroutine snoop_write_restart_file


  !> \brief Read SNOOP restart file for full system (sub=0) or
  !> one of the subsystems
  !> \author Kasper Kristensen
  !> \date August 2015
  subroutine snoop_read_restart_file(sub,EHF,Ecorr,file_exist)
    implicit none
    !> Which subsystem? (If sub=0, then we consider the full system)
    integer,intent(in) :: sub 
    !> HF and correlation energy 
    real(realk),intent(inout) :: EHF,Ecorr
    !> Did restart file exist? (If not, EHF and Ecorr are not modified)
    logical,intent(inout) :: file_exist
    integer :: funit
    character(len=17) :: FileName    

    ! Set file name
    call SNOOP_get_filename(sub,FileName)
    
    ! Does restart file exist?
    inquire(file=FileName,exist=file_exist)

    if(.not. file_exist) return

    ! Open file (always delete existing file if it exists)
    funit = -1
    call lsopen(funit,FileName,'OLD','UNFORMATTED')

    ! Read energies from file and close file
    read(funit) EHF
    read(funit) Ecorr
    call lsclose(funit,'KEEP')

  end subroutine snoop_read_restart_file


  !> Get filename for snoop restart file
  !> \author Kasper Kristensen
  !> \date August 2015  
  subroutine SNOOP_get_filename(sub,FileName)
    implicit none
    !> Which subsystem? (If sub=0, then we consider the full system)
    integer,intent(in) :: sub 
    !> File name for restart file
    character(len=17),intent(inout) :: FileName    

    ! Set file name
    FileName(1:5)='SNOOP'
    write(FileName(6:9),'(i4.4)') sub
    FileName(10:17)='.restart'

  end subroutine SNOOP_get_filename


end module snoop_tools_module
