!> @file
!> Module contains operations to create an atomic fragment and a pair atomic fragment
!> \author Marcin Ziolkowski and Kasper Kristensen

module atomic_fragment_operations

  use fundamental
  use precision
  use lstiming!, only: lstimer
  use typedeftype!, only: lsitem
  use DALTONINFO!, only: build_ccfragmentlsitem
  use files!,only:lsopen,lsclose
  use matrix_module!, only:matrix
  use matrix_operations!, only: mat_init, mat_zero, mat_to_full, mat_free
  use dec_typedef_module
  use memory_handling!, only: mem_alloc,mem_dealloc
  use IntegralInterfaceMod!, only: ii_get_mixed_overlap_full
  use Integralparameters!,only: AORdefault
#ifdef VAR_LSMPI
  use infpar_module
#endif


  ! DEC DEPENDENCIES (within deccc directory)   
  ! *****************************************
  use dec_fragment_utils
  use array2_simple_operations!, only: array2_init,array2_matmul,array2_free,array2_add, operator(*)
  use array4_simple_operations!, only: array4_init,array4_extract_eos_indices_both_schemes,&
  use orbital_operations!,only: get_number_of_orbitals_per_atom,orbital_init,&
!       & copy_orbital
  use ao_contractions!,only: max_batch_dimension,get_vovo_integrals,&
!       & mp2_integrals_memory_for_updating_arrays,&
!       & mp2_integrals_and_amplitudes
!       & array4_free, operator(*)


contains

  !> \brief Nullify everything and set variables to zero
  !> \author Marcin Ziolkowski
  !> \param fragment Molecular fragment
  subroutine atomic_fragment_nullify(fragment)

    implicit none
    type(ccatom), intent(inout) :: fragment


    fragment%occEOSidx => null()
    fragment%unoccEOSidx => null()
    fragment%occAOSidx => null()
    fragment%unoccAOSidx => null()
    fragment%coreidx => null()

    fragment%REDoccAOSidx => null()
    fragment%REDunoccAOSidx => null()

    fragment%occAOSorb => null()
    fragment%unoccAOSorb => null()

    fragment%idxo => null()
    fragment%idxu => null()

    fragment%atoms_idx => null()

    fragment%ypo => null()
    fragment%ypv => null()

    fragment%fock => null()
    fragment%ppfock => null()
    fragment%ccfock => null()
    fragment%qqfock => null()

    fragment%OccContribs => null()
    fragment%VirtContribs => null()

    fragment%OccMat => null()
    fragment%VirtMat => null()

    fragment%basisinfoisset=.false.
    fragment%atomic_number = 0
    fragment%atomic_number2 = 0
    fragment%noccEOS = 0
    fragment%nunoccEOS = 0
    fragment%noccAOS = 0
    fragment%nunoccAOS = 0
    fragment%number_atoms = 0
    fragment%number_basis = 0

    fragment%t1dims=0
    fragment%t1 => null()
    fragment%t1_occidx => null()
    fragment%t1_virtidx => null()


  end subroutine atomic_fragment_nullify



  !> \brief Initialize atomic fragment based on a list of atoms for the occupied and
  !> unoccupied orbital spaces. 
  !> For a super fragment calculation it is assumed that SF_atomlist and SF_nfrags
  !> are set BEFORE calling this routine!
  !> \author Kasper Kristensen
  subroutine atomic_fragment_init_atom_specific(MyAtom,natoms,Unocc_atoms, &
       & Occ_atoms,nOcc,nUnocc,OccOrbitals,UnoccOrbitals, &
       & MyMolecule,mylsitem,fragment,DoBasis,FA)

    implicit none
    !> Number of atom to build a fragment on
    integer, intent(in) :: MyAtom
    !> Number of atoms in full molecule
    integer, intent(in) :: nAtoms
    !> Logical vector telling which atoms are in unocc AOS
    logical, dimension(natoms), intent(in) :: Unocc_atoms
    !> Logical vector telling which atoms are in occ AOS
    logical, dimension(natoms), intent(in) :: Occ_atoms
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> Number of occupied orbitals in full molecule
    integer, intent(in) :: nocc
    !> Number of unoccupied orbitals in full molecule
    integer, intent(in) :: nunocc
    !> Information about DEC occupied orbitals
    type(ccorbital), dimension(nOcc), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(ccorbital), dimension(nUnocc), intent(in) :: UnoccOrbitals
    !> Fragment to construct
    type(ccatom), intent(inout) :: fragment
    !> Make fragment basis (MO coeff and Fock matrix for fragment)?
    logical, intent(in) :: DoBasis
    !> Use fragment-adapted orbitals?
    !> NOTE: If FA=.true. then the MO coefficients and MO Fock matrix
    !> elements in the fragment are NOT set here. Rather they should be
    !> set AFTER calling this routine (see init_fragment_adapted).
    logical,intent(in),optional :: FA
    integer :: j,idx, CentralAtom
    logical,pointer :: occ_list(:),unocc_list(:)

    ! Determine list of occupied orbitals assigned to one of the atoms in occ_atoms
    call mem_alloc(occ_list,nocc)
    occ_list=.false.
    do j=1,nocc
       CentralAtom=OccOrbitals(j)%centralatom
       if( Occ_atoms(CentralAtom) ) then  ! occupied orbital j is included in fragment
          occ_list(j)=.true.
       end if
    end do


    ! Determine list of unoccupied orbitals assigned to one of the atoms in unocc_atoms
    call mem_alloc(unocc_list,nunocc)
    unocc_list=.false.
    do j=1,nunocc
       CentralAtom=UnoccOrbitals(j)%centralatom
       if( Unocc_atoms(CentralAtom) ) then
          ! unoccupied orbital j is included in fragment
          unocc_list(j)=.true.
       end if
    end do


    ! Create fragment based on logical vectors for occupied and virtual AOS
    call atomic_fragment_init_orbital_specific(MyAtom,nunocc, nocc, unocc_list, &
       & occ_list,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,fragment,DoBasis)
    call mem_dealloc(occ_list)
    call mem_dealloc(unocc_list)


  end subroutine atomic_fragment_init_atom_specific



  !> \brief Initialize atomic fragment based on a list of specific AOS orbitals.
  !> For a super fragment calculation it is assumed that SF_atomlist and SF_nfrags
  !> are set BEFORE calling this routine!
  !> \author Kasper Kristensen
  subroutine atomic_fragment_init_orbital_specific(MyAtom,nunocc, nocc, unocc_list, &
       & occ_list,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,fragment,DoBasis,FA)

    implicit none
    !> Number of atom to build a fragment on
    integer, intent(in) :: MyAtom
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> Number of occupied orbitals in full molecule
    integer, intent(in) :: nocc
    !> Number of unoccupied orbitals in full molecule
    integer, intent(in) :: nunocc
    !> Logical vector telling which unoccupied AOS orbitals are included in fragment
    logical, dimension(nunocc), intent(in) :: Unocc_list
    !> Logical vector telling which occupied AOS orbitals are included in fragment
    logical, dimension(nocc), intent(in) :: Occ_list
    !> Information about DEC occupied orbitals
    type(ccorbital), dimension(nOcc), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(ccorbital), dimension(nUnocc), intent(in) :: UnoccOrbitals
    !> Fragment to construct
    type(ccatom), intent(inout) :: fragment
    !> Make fragment basis (MO coeff and Fock matrix for fragment)?
    logical, intent(in) :: DoBasis
    !> Use fragment-adapted orbitals?
    !> NOTE: If FA=.true. then the MO coefficients and MO Fock matrix
    !> elements in the fragment are NOT set here. Rather they should be
    !> set AFTER calling this routine (see init_fragment_adapted).
    logical,intent(in),optional :: FA
    integer :: j,idx,i,listidx,natoms,startidx
    integer :: CentralAtom
    real(realk) :: tcpu, twall
    logical,pointer :: occ_listEFF(:)


    ! Use fragment-adapted orbitals?
    fragment%fragmentadapted=.false.
    if(present(FA)) then
       if(FA) then
          fragment%fragmentadapted=.true.
       end if
    end if

    natoms = MyMolecule%natoms

    call LSTIMER('START',tcpu,twall,DECinfo%output)

    ! To be on the safe side, initialize everything to zero or null
    call atomic_fragment_nullify(fragment)
    nullify(fragment%parenthesis_t)

    ! -- Assign things --
    fragment%atomic_number = MyAtom


    ! SUPER FRAGMENT SETTINGS
    ! ***********************
    IsThisASuperFragment: if(DECinfo%SF) then ! Super fragment
       if(DECinfo%PL>0) write(DECinfo%output,*) '-- Initializing super fragment for atoms : ', fragment%SF_atomlist

       ! Sanity checks
       if( (.not. associated(fragment%SF_atomlist)) ) then
          call lsquit('atomic_fragment_init_orbital_specific: Super fragment calculation, &
               & but SF_atomlist has not been initiated!',-1)
       end if
       if( fragment%SF_nfrags /= size(fragment%SF_atomlist) ) then
          call lsquit('atomic_fragment_init_orbital_specific: Super fragment calculation, &
               & Numbers of atoms in atomlist is incorrect!',-1)
       end if

       ! Set remaining super fragment info
       fragment%SF=.true. ! super fragment

    else ! NOT a super fragment

       if(DECinfo%PL>0) write(DECinfo%output,*) '-- Initializing standard fragment on atom : ',MyAtom
       fragment%SF=.false. ! not super fragment
       fragment%SF_nfrags=1 ! contains just one fragment (itself)

       if(associated(fragment%SF_atomlist)) then  ! dealloc and realloc with the correct size of 1
          call mem_dealloc(fragment%SF_atomlist)
          fragment%SF_atomlist => null()
       end if
       call mem_alloc(fragment%SF_atomlist,fragment%SF_nfrags)
       fragment%SF_atomlist(1) = MyAtom
    end if IsThisASuperFragment

    ! ---- Done with super fragment settings



    ! Size of occupied AOS - number of "true" elements in logical occupied vector
    ! ***************************************************************************
    call mem_alloc(occ_listEFF,nocc)
    occ_listEFF = occ_list
    if(DECinfo%frozencore) then  ! for frozen core, core orbitals are removed from AOS orbital list
       occ_listEFF(1:MyMolecule%ncore) = .false.
    end if
    fragment%noccAOS = count(occ_listEFF)

    ! Size of unoccupied AOS - number of "true" elements in logical unoccupied vector
    ! *******************************************************************************
    fragment%nunoccAOS = count(unocc_list)

    ! Fragment-adapted information, for now set equal to local dimensions
    fragment%noccFA = fragment%noccAOS
    fragment%nunoccFA = fragment%nunoccAOS
    

    ! Size of occupied EOS
    ! ********************

    ! Occupied EOS orbitals are those assigned to central atom in fragment
    ! (or central atomS in case of a super fragment)

    ! Avoid including core orbitals in EOS if frozen core approx. is used
    if(DECinfo%frozencore) then
       startidx=MyMolecule%ncore+1 ! loop from first valence orbital
    else
       startidx=1 ! loop over all occ orbitals
    end if

    fragment%noccEOS=0
    ! loop over occupied orbitals (only valence for frozen core)
    OccEOSSize: do j=startidx,nOcc  
       CentralAtom=OccOrbitals(j)%centralatom

       ! Loop over atoms in super fragment list (just one atom, MyAtom, if it is not a superfragment)
       do i=1,fragment%SF_nfrags
          listidx=fragment%SF_atomlist(i)
          if( CentralAtom==listidx ) then ! Orbital is included in the EOS
             fragment%noccEOS = fragment%noccEOS + 1
          end if
       end do

    end do OccEOSSize


    ! Size of unoccupied EOS
    ! **********************

    ! Virtual EOS orbitals are those assigned to the central atom
    fragment%nunoccEOS=0
    do j=1,nunocc
       CentralAtom=UnoccOrbitals(j)%centralatom

       ! Loop over atoms in super fragment list (just one atom, Myatom, if it is not a superfragment)
       do i=1,fragment%SF_nfrags
          listidx=fragment%SF_atomlist(i)
          if( CentralAtom==listidx ) then ! Orbital is included in the EOS
             fragment%nunoccEOS = fragment%nunoccEOS + 1
          end if
       end do

    end do


    ! -- Assign occupied AOS indices
    call mem_alloc(fragment%occAOSidx,fragment%noccAOS)
    idx=0
    do i=1,nocc
       if(occ_listEFF(i)) then
          idx=idx+1
          fragment%occAOSidx(idx) = i
       end if
    end do
    if(idx /= fragment%noccAOS) &
         & call lsquit('atomic_fragment_init_orbital_specific: idx /= fragment%noccAOS',-1)

    ! -- Assign unoccupied AOS indices
    call mem_alloc(fragment%unoccAOSidx,fragment%nunoccAOS)
    idx=0
    do i=1,nunocc
       if(unocc_list(i)) then
          idx=idx+1
          fragment%unoccAOSidx(idx) = i
       end if
    end do
    if(idx /= fragment%nunoccAOS) &
         & call lsquit('atomic_fragment_init_orbital_specific: idx /= fragment%nunoccAOS',-1)


    ! Occupied EOS orbital indices
    ! ****************************
    call mem_alloc(fragment%occEOSidx,fragment%noccEOS)

    ! Occupied EOS orbitals are those assigned to central atom in fragment
    ! (or central atomS in case of a super fragment)

    fragment%occEOSidx = 0
    idx=0
    do j=startidx,nOcc   ! no core orbitals in EOS for frozen core approx
       CentralAtom=OccOrbitals(j)%centralatom

       ! Loop over atoms in super fragment list (just one atom - MyAtom - if it is not a superfragment)
       do i=1,fragment%SF_nfrags
          listidx=fragment%SF_atomlist(i)
          if( CentralAtom==listidx ) then ! Orbital is included in the EOS
             idx=idx+1
             fragment%occEOSidx(idx)=j
          end if
       end do

    end do
    if(idx /= fragment%noccEOS) then
       call lsquit('atomic_fragment_init_orbital_specific: idx /= fragment%noccEOS',-1)
    end if



    ! Virtual EOS orbital indices
    ! ***************************
    call mem_alloc(fragment%unoccEOSidx,fragment%nunoccEOS)

    ! Virtual EOS orbitals are those assigned to the central atom
    ! (or central atomS in case of a super fragment)

    fragment%unoccEOSidx = 0
    idx=0
    do j=1,nunocc
       CentralAtom=UnOccOrbitals(j)%centralatom

       ! Loop over atoms in super fragment list (just one atom - MyAtom - if it is not a superfragment)
       do i=1,fragment%SF_nfrags
          listidx=fragment%SF_atomlist(i)
          if( CentralAtom==listidx ) then ! Orbital is included in the EOS
             idx=idx+1
             fragment%unoccEOSidx(idx)=j
          end if
       end do

    end do
    if(idx /= fragment%nunoccEOS) then
       call lsquit('atomic_fragment_init_orbital_specific: idx /= fragment%nunoccEOS',DECinfo%output)
    end if



    ! Initially, set reduced fragment of lower accuracy to have the same orbitals as the original fragment
    ! ****************************************************************************************************
    fragment%REDnoccAOS = fragment%noccAOS
    fragment%REDnunoccAOS = fragment%nunoccAOS
    call mem_alloc(fragment%REDoccAOSidx,fragment%REDnoccAOS)
    call mem_alloc(fragment%REDunoccAOSidx,fragment%REDnunoccAOS)
    fragment%REDoccAOSidx = fragment%occAOSidx
    fragment%REDunoccAOSidx = fragment%unoccAOSidx

    ! Set core orbital info (redundant if we do not use frozen core approx, but do it anyway)
    call set_Core_orbitals_for_fragment(MyMolecule,nocc,OccOrbitals,Fragment) 

    ! Set EOS indices in the total list of fragment indices (idxo and idxu)
    call target_indices(fragment)


    ! Fragment-adapted orbital info
    !******************************
    fragment%RejectThr=0.0_realk
    ! correlation density matrices not set
    fragment%CDset=.false.  
    ! Transformation matrices between local/fragment-adapted bases not set
    fragment%FAtransSet=.false.

    ! Set atomic fragment extent info
    call init_atomic_fragment_extent(OccOrbitals,UnoccOrbitals,MyMolecule,fragment)

    ! Only create fragment basis (the information in the "expensive box" in ccatom type definition)
    ! if DoBasis is true.
    if(DoBasis) then
       call atomic_fragment_init_basis_part(nunocc, nocc, OccOrbitals,UnoccOrbitals,&
            & MyMolecule,mylsitem,fragment)
    else
       fragment%BasisInfoIsSet=.false.
    end if

    ! Information only used for pair, set to zero
    fragment%atomic_number2=0
    fragment%pairdist=0.0E0_realk


    ! Occ and virt energy contributions
    call mem_alloc(fragment%OccContribs,fragment%noccAOS)
    call mem_alloc(fragment%VirtContribs,fragment%nunoccAOS)
    fragment%OccContribs=0.0E0_realk
    fragment%VirtContribs=0.0E0_realk
    fragment%energies = 0.0E0_realk
    fragment%EoccFOP = 0.0_realk
    fragment%EvirtFOP = 0.0_realk
    fragment%LagFOP = 0.0_realk


    ! Information related to singles amplitudes - only relevant for CC2 and CCSD
    ! **************************************************************************
    ! By default we do not initialize any space for singles amplitudes
    fragment%t1_stored=.false.
    fragment%t1dims=0
    nullify(fragment%t1)
    nullify(fragment%t1_occidx)
    nullify(fragment%t1_virtidx)

    ! FLOP ACCOUNTING
    fragment%flops_slaves=0.0E0_realk
    fragment%ntasks=0

    ! TIME FOR LOCAL MPI SLAVES
    fragment%slavetime = 0.0E0_realk

    ! Free stuff
    call mem_dealloc(occ_listEFF)

    if(DECinfo%PL>0) then
       write(DECinfo%output,*)
       write(DECinfo%output,*) 'FRAGINIT: Initialized fragment :', Fragment%atomic_number
       write(DECinfo%output,'(a)')      ' -- Target --'
       write(DECinfo%output,'(a,i6)')   ' FRAGINIT: Occ EOS     : ',Fragment%noccEOS
       write(DECinfo%output,'(a,i6)')   ' FRAGINIT: Unocc EOS   : ',Fragment%nunoccEOS
       write(DECinfo%output,'(a)')      ' -- Target + Buffer --'
       write(DECinfo%output,'(a,i6)')   ' FRAGINIT: Occ AOS    : ',Fragment%noccAOS
       write(DECinfo%output,'(a,i6)')   ' FRAGINIT: Unocc AOS  : ',Fragment%nunoccAOS
       write(DECinfo%output,'(a,i6)')   ' FRAGINIT: Basis      : ',Fragment%number_basis
       write(DECinfo%output,*)
    end if

    call LSTIMER('FRAGMENT INIT',tcpu,twall,DECinfo%output)

  end subroutine atomic_fragment_init_orbital_specific





  !> Generate fragment where all quantites are expressed in the "fragment-adapted basis".
  !> Specifically, (1) the correlation density matrices in the local matrix are diagonalized, 
  !> (2) the unitary transformation matrices define the fragment-adapted basis,
  !> (3) fragment-adapted orbitals with eigenvalues smaller than 
  !> abs(LocalFragment%Rejectthr(1)) (occupied) or abs(LocalFragment%Rejectthr(2)) (virtual).
  !> are NOT included in the FAfragment (thereby the dimensions are reduced significantly
  !> compared to the fragment in the local basis).
  !> NOTE: The EOS orbitals (assigned to central atom in fragment) are left untouched, 
  !> so only the subset of AOS orbitals OUTSIDE the EOS are rotated.
  !> which are NOT in the EOS
  !> \author Kasper Kristensen
  !> \date February 2013
  subroutine fragment_adapted_driver(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals,&
       & LocalFragment,FAfragment)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> Information about DEC occupied orbitals
    type(ccorbital), dimension(MyMolecule%numocc), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(ccorbital), dimension(MyMolecule%numvirt), intent(in) :: UnoccOrbitals
    !> Atomic fragment where all quantities are expressed in local basis
    !> The transformation coefficients between AO and fragment-adapted basis
    !> will be stored in LocalFragment%CoccFA and LocalFragment%CunoccFA.
    type(ccatom), intent(inout) :: LocalFragment
    !> Atomic fragment where all quantities are expressed in fragment-adapted basis
    type(ccatom), intent(inout) :: FAfragment
    real(realk) :: tcpu,twall

    call LSTIMER('START',tcpu,twall,DECinfo%output)

    ! Get unitary transformation matrices from local to fragment-adapted basis 
    call fragment_adapted_transformation_matrices(LocalFragment)

    ! Initialize fragment-adapted fragment
    call init_fragment_adapted(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals,&
         & LocalFragment,FAfragment)

    call LSTIMER('SITE FRAGMENT',tcpu,twall,DECinfo%output)

  end subroutine fragment_adapted_driver



  !> Get transformation matrices from local basis to fragment-adapted basis
  !> (both for occ and virt spaces, see fragment_adapted_driver for details).
  !> \author Kasper Kristensen
  !> \date February 2013
  subroutine fragment_adapted_transformation_matrices(LocalFragment)

    implicit none
    !> Atomic fragment where all quantities are expressed in local basis
    type(ccatom), intent(inout) :: LocalFragment
    real(realk),pointer :: VirtMat(:,:),OccMat(:,:),tmpeival(:),tmpU(:,:)
    integer :: i,j,ix,jx
    integer :: noccTRANS,nvirtTRANS,noccEOS,nvirtEOS,offset,nocc,nvirt
    real(realk),pointer :: OU(:,:),VU(:,:),Oeival(:),Veival(:),OUred(:,:),VUred(:,:)
    logical,pointer :: virtEOS(:),occEOS(:), VirtOrbs(:), OccOrbs(:)

    ! Init stuff
    nocc = LocalFragment%noccAOS
    nvirt = LocalFragment%nunoccAOS
    noccEOS = LocalFragment%noccEOS
    nvirtEOS = LocalFragment%nunoccEOS
    call mem_alloc(OU,nocc,nocc)
    call mem_alloc(VU,nvirt,nvirt)
    call mem_alloc(Oeival,nocc)
    call mem_alloc(Veival,nvirt)

    ! Number of orbitals to transform (leave EOS untouched in standard DEC calculation)
    noccTRANS = nocc-noccEOS
    nvirtTRANS = nvirt-nvirtEOS

    ! Special debug case: AOS=EOS, transform all orbitals
    if(LocalFragment%noccEOS == LocalFragment%noccAOS) then
       noccTRANS = nocc
    end if
    if(LocalFragment%nunoccEOS == LocalFragment%nunoccAOS) then
       nvirtTRANS = nvirt
    end if

    ! Sanity check
    if( (noccTRANS < 1) .or. (nvirtTRANS < 1) ) then
       print *, 'EOS: nocc, nvirt', noccEOS, nvirtEOS
       print *, 'AOS: nocc, nvirt', nocc, nvirt
       call lsquit('fragment_adapted_driver: No orbitals to rotate!',-1)
    end if




    ! ============================================================
    !                          VIRTUAL SPACE                     !
    ! ============================================================

    ! Logical vector for virtual EOS
    call mem_alloc(virtEOS,nvirt)
    virtEOS=.false.
    do i=1,nvirtEOS
       virtEOS(LocalFragment%idxu(i)) = .true.
    end do

    ! Virtual correlation density matrix
    call mem_alloc(VirtMat,nvirtTRANS,nvirtTRANS)
    call mem_alloc(tmpeival,nvirtTRANS)
    call mem_alloc(tmpU,nvirtTRANS,nvirtTRANS)

    if(LocalFragment%nunoccEOS ==LocalFragment%nunoccAOS) then
       ! Transform all orbitals
       VirtMat = LocalFragment%VirtMat
    else

       ! Consider only AOS orbitals that are NOT in EOS
       ix=0
       do i=1,nvirt
          if(virtEOS(i)) cycle
          ix=ix+1

          jx=0
          do j=1,nvirt
             if(virtEOS(j)) cycle
             jx = jx+1
             VirtMat(jx,ix) = LocalFragment%VirtMat(j,i)
          end do
       end do

    end if

    ! Diagonalize virtual correlation density matrix to define fragment-adapted virtual AOS orbitals
    call solve_eigenvalue_problem_unitoverlap(nvirtTRANS,VirtMat,tmpeival,tmpU)


    ! Now the fragment-adapted orbitals (psi) are given from local orbitals (phi) as:
    ! psi(c) = sum_a tmpU(a,c) phi(a)

    ! However, only orbitals OUTSIDE the EOS are described by tmpU while the local
    ! EOS orbitals have not been touched. To understand the ordering in the
    ! output VU matrix consider the example where there are 5 AOS orbitals,
    ! and where 2 of these are EOS orbitals.
    ! In the AOS list of orbitals the EOS orbitals are number 2 and 4, i.e. 
    ! LocalFragment%idxu = [2,4].
    !
    ! The VU matrix has dimensions (local basis, fragment-adapted basis).
    !
    ! For the local basis the EOS orbitals thus correspond to rows 2 and 4.
    ! Later on we need to remove some of the fragment-adapted orbitals, and therefore
    ! it is most convenient to place all EOS orbitals BEFORE the remaining AOS
    ! orbitals. Thus, the EOS orbitals for the fragment-adapted basis (identical
    ! to the EOS orbitals in the local basis) are described by columns 1 and 2.
    ! VU therefore has the structure:
    ! 
    ! 
    !         0     0     X     X     X    
    !         1     0     0     0     0    
    ! VU =    0     0     X     X     X     
    !         0     1     0     0     0    
    !         0     0     X     X     X
    !
    ! where X is, in general, NOT 0 or 1.
    VU = 0.0_realk

    ! Set EOS indices of VU
    ! ---------------------
    ! Also set eigenvalue output for the EOS orbitals to 1, 
    ! although in principle these are not defined 
    ! - but at least they are then initialized to something.
    do i=1,nvirtEOS
       ! VU index 1: Local index   (correspond to 2 and 4 in example above)
       ! VU index 2: Fragment-adapted index   (correspond to 1 and 2 in example above)
       VU(LocalFragment%idxu(i),i) = 1.0_realk    
       Veival(i) = 1.0_realk
    end do


    ! Set remaining AOS  indices (not EOS)
    ! ------------------------------------
    offset = nvirt-nvirtTRANS    ! offset to put all non-EOS orbitals after the EOS orbitals
    do j=1,nvirtTRANS  ! loop over fragment-adapted orbital indices
       Veival(j+offset) = tmpeival(j)
       ix=0
       localloop1: do i=1,nvirt  ! loop over local orbital indices
          if(virtEOS(i)) cycle localloop1   ! not consider local EOS, which have already been set
          ix=ix+1
          VU(i,j+offset) = tmpU(ix,j)
       end do localloop1
    end do


    call mem_dealloc(tmpeival)
    call mem_dealloc(tmpU)




    ! ============================================================
    !                         OCCUPIED SPACE                     !
    ! ============================================================
    ! Same strategy as for virtual space above.

    ! Logical vector for occupied EOS
    call mem_alloc(occEOS,nocc)
    occEOS=.false.
    do i=1,noccEOS
       occEOS(LocalFragment%idxo(i)) = .true.
    end do

    ! Occupied correlation density matrix
    call mem_alloc(OccMat,noccTRANS,noccTRANS)
    call mem_alloc(tmpeival,noccTRANS)
    call mem_alloc(tmpU,noccTRANS,noccTRANS)

    if(LocalFragment%noccEOS ==LocalFragment%noccAOS) then
       ! Transform all orbitals
       OccMat = LocalFragment%OccMat
    else

       ! Consider only AOS orbitals that are NOT in EOS
       ix=0
       do i=1,nocc
          if(occEOS(i)) cycle
          ix=ix+1

          jx=0
          do j=1,nocc
             if(occEOS(j)) cycle
             jx = jx+1
             OccMat(jx,ix) = LocalFragment%OccMat(j,i)
          end do
       end do

    end if

    ! Diagonalize occupied correlation density matrix to define fragment-adapted occupied AOS orbitals
    call solve_eigenvalue_problem_unitoverlap(noccTRANS,OccMat,tmpeival,tmpU)

    ! Set blocks of OU in the same way as for VU above.
    OU = 0.0_realk

    do i=1,noccEOS
       OU(LocalFragment%idxo(i),i) = 1.0_realk
       Oeival(i) = 1.0_realk
    end do

    offset = nocc-noccTRANS
    do j=1,noccTRANS  ! loop over fragment-adapted orbital indices
       Oeival(j+offset) = tmpeival(j)
       ix=0
       localloop2: do i=1,nocc  ! loop over local orbital indices
          if(occEOS(i)) cycle localloop2
          ix=ix+1
          OU(i,j+offset) = tmpU(ix,j)
       end do localloop2
    end do

    call mem_dealloc(tmpeival)
    call mem_dealloc(tmpU)
    call mem_dealloc(VirtMat)
    call mem_dealloc(OccMat)
    call mem_dealloc(occEOS)
    call mem_dealloc(virtEOS)



    !=====================================================================================
    ! At this point OU and VU contains the transformation matrices between local and
    ! fragment-adapted bases (first index: local, second index: fragment-adapted).
    ! We now want to remove the fragment-adapted orbitals with small eigenvalues
    ! as defined by LocalFragment%RejectThr.
    !=====================================================================================


    ! Find which fragment-adapted orbitals to include
    call mem_alloc(VirtOrbs,nvirt)
    call mem_alloc(OccOrbs,nocc)
    call which_fragment_adapted_orbitals(LocalFragment,nocc,nvirt,Oeival,Veival,OccOrbs,VirtOrbs)
    write(DECinfo%output,'(1X,a,i6,a,i6,a)') 'FOP Removed ', nvirt-count(VirtOrbs), ' of ', &
         & nvirt-LocalFragment%nunoccEOS, ' fragment-adapted virtual orbitals'
    write(DECinfo%output,'(1X,a,i6,a,i6,a)') 'FOP Removed ', nocc-count(OccOrbs), ' of ', &
         & nocc-LocalFragment%noccEOS, ' fragment-adapted occupied orbitals'


    ! Set fragment-adapted MO coefficients
    ! ------------------------------------

    ! The fragment-adapted occ orbitals (psi) are given from local orbitals (phi) as:
    !
    ! psi(i) = sum_k OU(k,i) phi(k)         (only for certain "i" defined by OccOrbs)
    ! 
    ! We want to include only the orbitals "i" defined by OccOrbs to generate a
    ! matrix OUred where:
    ! Dimension 1: Number of occupied AOS orbitals in LOCAL basis
    ! Dimension 2: Number of occupied AOS orbitals in SITE SPECIFIC basis
    ! In general dimension 1 is larger than dimension 2.
    LocalFragment%noccFA = count(OccOrbs)
    call mem_alloc(OUred,LocalFragment%noccAOS,LocalFragment%noccFA)
    ix=0
    do i=1,nocc
       if(OccOrbs(i)) then
          ix=ix+1
          OUred(:,ix) = OU(:,i)
       end if
    end do


    ! Same for virtual space
    LocalFragment%nunoccFA = count(VirtOrbs)
    call mem_alloc(VUred,LocalFragment%nunoccAOS,LocalFragment%nunoccFA)
    ix=0
    do i=1,nvirt
       if(VirtOrbs(i)) then
          ix=ix+1
          VUred(:,ix) = VU(:,i)
       end if
    end do


    ! Set fragment-adapted (FA) orbital coefficients: AO-->FA basis
    ! -------------------------------------------------------------

    ! The FA MO coefficient matrix for the occupied space
    ! can be found as (dimensions given below):
    !
    ! (MO fragment-adapted basis)    =     (MO local basis)      *      OUred
    !     (nbasis,noccFA)            (nbasis,noccLOCAL)      (noccLOCAL,noccFA)
    !
    if(LocalFragment%FAtransSet) then
       call mem_dealloc(LocalFragment%CoccFA)
       call mem_dealloc(LocalFragment%CunoccFA)
    end if
    nullify(LocalFragment%CoccFA,LocalFragment%CunoccFA)
    call mem_alloc(LocalFragment%CoccFA,LocalFragment%number_basis,LocalFragment%noccFA)
    call dec_simple_dgemm(LocalFragment%number_basis,LocalFragment%noccAOS,&
         & LocalFragment%noccFA,LocalFragment%ypo,OUred,LocalFragment%CoccFA,'n','n')

    ! Virtual space (same strategy as for occ space)
    call mem_alloc(LocalFragment%CunoccFA,LocalFragment%number_basis,LocalFragment%nunoccFA)
    call dec_simple_dgemm(LocalFragment%number_basis,LocalFragment%nunoccAOS,&
         & LocalFragment%nunoccFA,LocalFragment%ypv,VUred,LocalFragment%CunoccFA,'n','n')
    LocalFragment%FAtransSet=.true.

    call mem_dealloc(OccOrbs)
    call mem_dealloc(VirtOrbs)
    call mem_dealloc(OU)
    call mem_dealloc(VU)
    call mem_dealloc(OUred)
    call mem_dealloc(VUred)
    call mem_dealloc(Oeival)
    call mem_dealloc(Veival)

  end subroutine fragment_adapted_transformation_matrices



  !> Initialize fragment where all quantites are expressed in the "fragment-adapted basis".
  !> See fragment_adapted_driver for details.
  !> \author Kasper Kristensen
  !> \date February 2013
  subroutine init_fragment_adapted(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals,&
       & LocalFragment,FAfragment)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> Information about DEC occupied orbitals
    type(ccorbital), dimension(MyMolecule%numocc), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(ccorbital), dimension(MyMolecule%numvirt), intent(in) :: UnoccOrbitals
    !> Atomic fragment where all quantities are expressed in local basis
    !> and where site specific MO coefficients have been stored in 
    !> CoccFA and CunoccFA (see fragment_adapted_transformation_matrices).
    type(ccatom),intent(inout) :: LocalFragment
    !> Atomic fragment where all quantities are expressed in fragment-adapted basis
    !> (should be an "empty" fragment at input)
    type(ccatom), intent(inout) :: FAfragment
    integer :: i,ix,MyAtom,nocc,nvirt
    logical,pointer :: VirtLocal(:), OccLocal(:)
    logical :: pairfrag,SFsave

    ! Dims
    nocc = LocalFragment%noccAOS
    nvirt = LocalFragment%nunoccAOS

    ! Is this a pair fragment?
    if(LocalFragment%atomic_number2/=0) then
       pairfrag = .true.  ! yes it is, two central atoms have been defined
    else
       pairfrag = .false.  ! no it is not, only 1 central atom
    end if


    ! =================================================
    !            Init fragment-adapted fragment          !
    ! =================================================
    ! First, we initialize fragment-adapted fragment identical to local fragment.
    ! Second, we modify only the properties of the fragment-adapted fragment 
    ! which differ from local fragment.

    ! Init fragment identical to local fragment
    ! *****************************************
    call mem_alloc(VirtLocal,MyMolecule%numvirt)
    VirtLocal=.false.
    do i=1,LocalFragment%nunoccAOS
       VirtLocal(LocalFragment%unoccAOSidx(i)) = .true.
    end do
    call mem_alloc(OccLocal,MyMolecule%numocc)
    OccLocal=.false.
    do i=1,LocalFragment%noccAOS
       OccLocal(LocalFragment%occAOSidx(i)) = .true.
    end do

    ! Just to be sure, nullify site fragment to zero
    call atomic_fragment_nullify(FAfragment)

    ! Special info required for pair fragment (similar to merged_fragment_init subroutine)
    ! (Also for fragment which is formally a single super fragment)
    if(pairfrag .or. LocalFragment%SF) then
       ! Ensure super fragment keyword is set (should be replaced by more elegant solution)
       SFsave = DECinfo%SF
       DECinfo%SF = .true.
       
       ! Copy super fragment info (pair fragment is effectively considered to be a super fragment)
       FAfragment%SF_nfrags = localfragment%SF_nfrags
       call mem_alloc(FAfragment%SF_atomlist,FAfragment%SF_nfrags)
       do i=1,FAfragment%SF_nfrags
          FAfragment%SF_atomlist(i) = localfragment%SF_atomlist(i)
       end do

       MyAtom=0  ! no central atom because it is a pair fragment

    else
       MyAtom = LocalFragment%atomic_number
    end if


    call atomic_fragment_init_orbital_specific(MyAtom,MyMolecule%numvirt,&
         & MyMolecule%numocc, VirtLocal, OccLocal,OccOrbitals,UnoccOrbitals,MyMolecule,&
         & mylsitem,FAfragment,.true.,FA=.true.)

    ! Copy special pair fragment info
    if(pairfrag) then
       FAfragment%atomic_number = LocalFragment%atomic_number
       FAfragment%atomic_number2 = LocalFragment%atomic_number2
       FAfragment%pairdist = LocalFragment%pairdist
    end if
    

    ! Reset information for site fragment that differs from local fragment
    ! ********************************************************************
    ! AOS dimensions
    FAfragment%noccAOS = LocalFragment%noccFA
    FAfragment%nunoccAOS = LocalFragment%nunoccFA

    ! Total number of occupied orbitals
    if(DECinfo%frozencore) then
       ! core + valence (AOS)
       FAfragment%nocctot = FAfragment%ncore + FAfragment%noccAOS
    else 
       ! AOS already contains core orbitals
       FAfragment%nocctot = FAfragment%noccAOS
    end if


    ! AOS indices are now ill-defined because the fragment-adapted orbitals
    ! are not defined in the full molecular local basis!
    ! We now set the entries to -1 to hopefully detect if someone by mistake tries to use the indices!
    FAfragment%occAOSidx = -1
    FAfragment%unoccAOSidx = -1

    ! Do the same for reduced space indices. However, these should be removed at some point!
    FAfragment%REDoccAOSidx = -1
    FAfragment%REDunoccAOSidx = -1
    


    ! EOS orbitals in the AOS list of orbitals
    ! ----------------------------------------
    ! Now the orbitals are listed with EOS orbitals before  the remaining AOS orbitals
    ! (see fragment_adapted_transformation_matrices), so the EOS indices in the total list
    ! of AOS orbitals are simply the first noccEOS (or nunoccEOS) orbitals.
    do i=1,FAfragment%noccEOS
       FAfragment%idxo(i) = i
    end do
    do i=1,FAfragment%nunoccEOS
      FAfragment%idxu(i) = i
    end do


    ! Set fragment-adapted MO coefficients
    ! ------------------------------------

    ! These should be stored in LocalFragment, quit if this is not the case
    if(.not. LocalFragment%FAtransSet) then
       call lsquit('init_fragment_adapted: Fragment-adapted MO coefficients &
            & have not been set!',-1)
    end if

    ! Copy MOs
    call mem_alloc(FAfragment%ypo,FAfragment%number_basis,FAfragment%noccAOS)
    FAfragment%ypo = LocalFragment%CoccFA
    call mem_alloc(FAfragment%ypv,FAfragment%number_basis,FAfragment%nunoccAOS)
    FAfragment%ypv = LocalFragment%CunoccFA



    ! Set fragment-adapted MO Fock matrix
    ! -----------------------------------
    ! Transform local FOck matrix:  F_fragmentadapted = C^T F_AO C
    ! C: MO coefficient matrix from AO to FA basis.

    ! Occ space
    call mem_alloc(FAfragment%ppfock,FAfragment%noccAOS,FAfragment%noccAOS)
    call dec_simple_basis_transform1(FAfragment%number_basis,FAfragment%noccAOS,&
         & FAfragment%ypo,FAfragment%fock,FAfragment%ppfock)

    ! Virt space
    call mem_alloc(FAfragment%qqfock,FAfragment%nunoccAOS,FAfragment%nunoccAOS)
    call dec_simple_basis_transform1(FAfragment%number_basis,FAfragment%nunoccAOS,&
         & FAfragment%ypv,FAfragment%fock,FAfragment%qqfock)


    ! Reset keyword (more elegant solution with soon be introduced)
    if(pairfrag .or. LocalFragment%SF) then
       ! Reset super fragment keyword (more elegant solution is coming up...)
       DECinfo%SF = SFsave
    end if

    call mem_dealloc(VirtLocal)
    call mem_dealloc(OccLocal)


  end subroutine init_fragment_adapted




  !> Set logical vectors describing which fragment-adapted orbitals to include in
  !> occupied and virtual AOS for fragment.
  !> This is done by including only those orbitals where the eigenvalue is
  !> above the threshold store in LocalFragment%RejectThr(1) (occupied space)
  !> or LocalFragment%RejectThr(2) (virtual space).
  !> \author Kasper Kristensen
  !> \date February 2013
  subroutine which_fragment_adapted_orbitals(LocalFragment,nocc,nvirt,Oeival,Veival,OccOrbs,VirtOrbs)

    implicit none
    !> Atomic fragment where all quantities are expressed in local basis
    type(ccatom),intent(inout) :: LocalFragment
    !> Number of occupied AOS orbitals for local fragment
    integer, intent(in) :: nocc
    !> Number of virtual AOS orbitals for local fragment
    integer, intent(in) :: nvirt
    !> Eigenvalues for occupied fragment-adapted orbitals (see fragment_adapted_driver)
    real(realk),intent(in) :: Oeival(nocc)
    !> Eigenvalues for virtual fragment-adapted orbitals (see fragment_adapted_driver)
    real(realk),intent(in) :: Veival(nvirt)
    !> Logical vector describing which fragment-adapted occ orbitals to include
    logical,intent(inout) :: OccOrbs(nocc)
    !> Logical vector describing which fragment-adapted virt orbitals to include
    logical,intent(inout) :: VirtOrbs(nvirt)
    integer :: i,ix,maxvirtidx,maxoccidx
    real(realk) :: maxvirt,maxocc
   

    ! ================================================
    !    Find which virtual AOS orbitals to include  !
    ! ================================================
    VirtOrbs=.false.
    ! Always include all EOS orbitals 
    ! (listed before remaining orbitals, see example in fragment_adapted_transformation_matrices)
    do i=1,LocalFragment%nunoccEOS
       VirtOrbs(i) = .true.
    end do

    ! Include fragment-adapted orbitals in the remaining part of AOS (NOT EOS) for which
    ! eigenvalue is above threshold value.
    maxvirt=-1.0_realk
    maxvirtidx=0
    do i=1,nvirt  

       ! not consider EOS orbitals already included
       if(.not. VirtOrbs(i)) then 

          if( abs(Veival(i)) > LocalFragment%RejectThr(2) ) then
             VirtOrbs(i)=.true.
          end if
          if(abs(Veival(i)) > maxvirt) then
             maxvirt = abs(Veival(i))
             maxvirtidx=i
          end if
       end if
    end do
    ! Sanity precaution: Always include at least one orbital outside EOS,
    !                    even if it falls below the threshold.
    if(LocalFragment%nunoccEOS/=LocalFragment%nunoccAOS) VirtOrbs(maxvirtidx) = .true.


    ! =================================================
    !    Find which occupied AOS orbitals to include  !
    ! =================================================
    ! Same strategy as for virtual orbitals above
    OccOrbs=.false.
    do i=1,LocalFragment%noccEOS
       OccOrbs(i) = .true.
    end do

    maxocc=-1.0_realk
    maxoccidx=0
    do i=1,nocc  

       if(.not. OccOrbs(i) .and. LocalFragment%noccEOS/=LocalFragment%noccAOS) then 
          if( abs(Oeival(i)) > LocalFragment%RejectThr(1) ) then
             OccOrbs(i)=.true.
          end if
          if(abs(Oeival(i)) > maxocc) then
             maxocc = abs(Oeival(i))
             maxoccidx=i
          end if
       end if
    end do
    ! Sanity precaution: Always include at least one orbital outside EOS,
    !                    even if it falls below the threshold.
    if(LocalFragment%noccEOS/=LocalFragment%noccAOS) OccOrbs(maxoccidx) = .true.


  end subroutine which_fragment_adapted_orbitals





  !> \brief Initialize the basis part of fragment, i.e. the information listed
  !> in the "expensive box" in the ccatom type definition.
  !> \author Kasper Kristensen
  !> \date March 2012
  subroutine atomic_fragment_init_basis_part(nunocc, nocc, OccOrbitals,UnoccOrbitals,&
       & MyMolecule,mylsitem,fragment)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> Number of occupied orbitals in full molecule
    integer, intent(in) :: nocc
    !> Number of unoccupied orbitals in full molecule
    integer, intent(in) :: nunocc
    !> Information about DEC occupied orbitals
    type(ccorbital), dimension(nOcc), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(ccorbital), dimension(nUnocc), intent(in) :: UnoccOrbitals
    !> Fragment to construct
    type(ccatom), intent(inout) :: fragment
    integer :: j,idx

    ! Basis info
    call atomic_fragment_basis(fragment,MyMolecule)


    ! lsitem
#ifdef VAR_LSMPI
    ! Quick fix such that lsitem is never constructed for global master
    ! as this will destroy the overall MPI framework.
    if(infpar%mynum/=infpar%master) then
       call build_ccfragmentlsitem(mylsitem,fragment%mylsitem,fragment%atoms_idx,&
         fragment%number_atoms,DECinfo%output,0)
    end if
#else
       call build_ccfragmentlsitem(mylsitem,fragment%mylsitem,fragment%atoms_idx,&
         fragment%number_atoms,DECinfo%output,0)
#endif

    ! normalize if requested
    if(DECinfo%NormalizeFragment) then
       call atomic_fragment_normalize(fragment)
    end if

    ! Basis info has now been set
    fragment%BasisInfoIsSet=.true.

  end subroutine atomic_fragment_init_basis_part



  !> \brief Normalize atomic fragment basis
  !> \author Marcin Ziolkowski
  !> \param fragment Atomic fragment/atomic pair fragment

  ! Normalize basis for an atomic fragment. truncated basis in not normalized if
  ! molecular orbitals expansion coefficients are taken only for selected atoms.
  ! this subroutine normalized occupied and unoccupied orbitals for the
  ! fragment.
  subroutine atomic_fragment_normalize(fragment)


    implicit none
    type(ccatom), intent(inout) :: fragment
    integer :: nocc,nvirt,nbasis
    real(realk), pointer :: aoS(:,:), moSp(:,:)
    real(realk) :: norm_factor_particle
    integer :: mu,i

    write(DECinfo%output,*) 'WARNING: Normalization of basis for fragment is not recommended!'
    write(DECinfo%output,*) 'WARNING: It is both slow and may give unstable results!'

    nocc = fragment%noccAOS
    nvirt = fragment%nunoccAOS
    nbasis = fragment%number_basis

    ! calculate overlap matrix
    call mem_alloc(aoS,nbasis,nbasis)
    call ii_get_mixed_overlap_full(DECinfo%output,DECinfo%output,&
         & fragment%mylsitem%setting,aoS,nbasis,nbasis,AORdefault,AORdefault)

    ! -- occupied

    ! transform overlap to occupied
    call mem_alloc(moSp,nocc,nocc)  ! particle
    moSp = 0.0E0_realk
    moSp = matmul(matmul(transpose(fragment%ypo),aoS),fragment%ypo)

    ! normalize MOs coefficients
    do i=1,nocc
       norm_factor_particle = 1.0E0_realk/sqrt(moSp(i,i))
       do mu=1,nbasis
          fragment%ypo(mu,i) = fragment%ypo(mu,i) * norm_factor_particle
       end do
    end do

    call mem_dealloc(moSp)

    ! -- virtuals

    call mem_alloc(moSp,nvirt,nvirt)
    moSp = 0.0E0_realk
    moSp = matmul(matmul(transpose(fragment%ypv),aoS),fragment%ypv)

    do i=1,nvirt
       norm_factor_particle = 1.0E0_realk/sqrt(moSp(i,i))
       do mu=1,nbasis
          fragment%ypv(mu,i) = fragment%ypv(mu,i) * norm_factor_particle
       end do
    end do

    call mem_dealloc(moSp)


    call mem_dealloc(aoS)

    return
  end subroutine atomic_fragment_normalize




  !> \brief Initialize all atomic fragments such that they all contain the entire molecule.
  !> This is used only to simulate a full calculation within the DEC framework.
  !> \author Kasper Kristensen
  !> \date August 2012
  subroutine fragment_init_simulate_full(myAtom,nunocc, nocc, OccOrbitals,UnoccOrbitals,&
       & MyMolecule,mylsitem,fragment,DoBasis)

    implicit none
    !> Central atom in fragment
    integer,intent(in) :: MyAtom
    !> Number of occupied orbitals in full molecule
    integer, intent(in) :: nocc
    !> Number of unoccupied orbitals in full molecule
    integer, intent(in) :: nunocc
    !> Information about DEC occupied orbitals
    type(ccorbital), dimension(nOcc), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(ccorbital), dimension(nUnocc), intent(in) :: UnoccOrbitals
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> Fragment to construct
    type(ccatom), intent(inout) :: fragment
    !> Make fragment basis (MO coeff and Fock matrix for fragment)?
    logical, intent(in) :: DoBasis
    logical, dimension(nunocc) :: Unocc_list
    logical, dimension(nocc) :: Occ_list

    ! All orbitals included in fragment
    unocc_list=.true.
    occ_list=.true.

    call atomic_fragment_init_orbital_specific(MyAtom,nunocc, nocc, unocc_list, &
         & occ_list,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,fragment,DoBasis)

  end subroutine fragment_init_simulate_full


  !> \brief Merge two single fragment into one pair fragment
  !> NOTE: The pair fragment is considered as a super fragment.
  !> \author Kasper Kristensen
  !> \date December 2011
  subroutine merged_fragment_init(Fragment1,Fragment2,nunocc, nocc, natoms, &
       & OccOrbitals,UnoccOrbitals,DistanceTable,MyMolecule,mylsitem,DoBasis,pairfragment)


    implicit none
    !> Fragment 1 in pair
    type(ccatom),intent(inout) :: fragment1
    !> Fragment 2 in pair
    type(ccatom),intent(inout) :: fragment2
    !> Number of occupied orbitals in full molecule
    integer, intent(in) :: nocc
    !> Number of unoccupied orbitals in full molecule
    integer, intent(in) :: nunocc
    !> Number of atoms for full molecule
    integer, intent(in) :: natoms
    !> Information about DEC occupied orbitals
    type(ccorbital), dimension(nOcc), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(ccorbital), dimension(nUnocc), intent(in) :: UnoccOrbitals
    !> Distance table for all atoms in the molecule
    real(realk), intent(in) :: DistanceTable(natoms,natoms)
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> Make pair fragment basis (MO coeff and Fock matrix for fragment)?
    logical, intent(in) :: DoBasis
    !> Pair fragment to be constructed
    type(ccatom),intent(inout) :: pairfragment
    logical, dimension(nunocc) :: Unocc_list
    logical, dimension(nocc) :: Occ_list
    logical :: SFsave,pairreduction
    logical,pointer :: EOSatoms(:)
    integer :: i,j,idx
    real(realk) :: pairdist



    ! Sanity check : No overlap between EOS atoms in fragment1 and fragment2
    ! **********************************************************************

    ! We also set a logical vector for EOS atoms to be used below
    call mem_alloc(EOSatoms,natoms)
    EOSatoms=.false.

    do i=1,fragment1%SF_nfrags
       EOSatoms(fragment1%SF_atomlist(i))=.true.
       do j=1,fragment2%SF_nfrags
          EOSatoms(fragment2%SF_atomlist(j))=.true.
          if(fragment1%SF_atomlist(i) == fragment2%SF_atomlist(j)) then
             call lsquit('merged_fragment_init: Overlap between EOS atoms for fragment 1 and 2 ',DECinfo%output)
          end if
       end do
    end do

    ! Save super fragment info - the super fragment keyword must be true
    ! when atomic_fragment_init_orbital_specific is called.
    ! At the end, the keyword is restored
    SFsave = DECinfo%SF
    DECinfo%SF = .true.

    ! just to be sure, initialize everything to zero or null
    call atomic_fragment_nullify(pairfragment)

    ! Find pair distance
    pairdist = get_distance_between_superfragments(fragment1,&
         & fragment2,natoms,DistanceTable)
    ! Use reduced fragments for large distances
    if(DECinfo%PL>0) write(DECinfo%output,'(1X,a,F12.3)') 'Pair distance (Angstrom)      = '&
         &, PairDist*bohr_to_angstrom
    if(DECinfo%PL>0) write(DECinfo%output,'(1X,a,F12.3)') 'Pair reduction thr (Angstrom) = '&
         &, DECinfo%PairReductionDistance*bohr_to_angstrom
    if(pairdist > DECinfo%PairReductionDistance) then
       pairreduction=.true.
       if(DECinfo%PL>0) write(DECinfo%output,*) '--> Initiating pair fragment using REDUCED fragments'
    else
       pairreduction=.false.
       if(DECinfo%PL>0) write(DECinfo%output,*) '--> Initiating pair fragment using STANDARD fragments'
    end if



    ! Find EOS atoms for pair: Union of EOS atoms for fragment1 and fragment2
    ! ***********************************************************************

    ! Since there is no overlap between the EOS atoms (sanity check 2),
    ! the pair EOS size is just the sum of the individual sizes
    pairfragment%SF_nfrags = fragment1%SF_nfrags + fragment2%SF_nfrags
    call mem_alloc(pairfragment%SF_atomlist,pairfragment%SF_nfrags)
    if(pairfragment%SF_nfrags /= count(EOSatoms)) then
       write(DECinfo%output,*) 'Atoms in fragment1 = ', fragment1%SF_atomlist
       write(DECinfo%output,*) 'Atoms in fragment2 = ', fragment2%SF_atomlist
       call lsquit('merged_fragment_init: pairfragment%SF_nfrags /= count(EOSatoms)',DECinfo%output)
    end if

    ! Set indices for EOS atoms for pair
    idx=0
    do i=1,natoms
       if(EOSatoms(i)) then
          idx=idx+1
          pairfragment%SF_atomlist(idx) = i
       end if
    end do


    ! Occupied AOS and unoccupied AOS for pair: Union of input fragments
    ! ******************************************************************

    occ_list=.false.
    unocc_list=.false.

    if(pairreduction) then  ! use reduced orbital space

       ! Occupied
       do i=1,fragment1%REDnoccAOS
          occ_list(fragment1%REDoccAOSidx(i))=.true.
       end do
       do i=1,fragment2%REDnoccAOS
          occ_list(fragment2%REDoccAOSidx(i))=.true.
       end do

       ! Unoccupied
       do i=1,fragment1%REDnunoccAOS
          unocc_list(fragment1%REDunoccAOSidx(i))=.true.
       end do
       do i=1,fragment2%REDnunoccAOS
          unocc_list(fragment2%REDunoccAOSidx(i))=.true.
       end do

    else

       ! Occupied
       do i=1,fragment1%noccAOS
          occ_list(fragment1%occAOSidx(i))=.true.
       end do
       do i=1,fragment2%noccAOS
          occ_list(fragment2%occAOSidx(i))=.true.
       end do

       ! Unoccupied
       do i=1,fragment1%nunoccAOS
          unocc_list(fragment1%unoccAOSidx(i))=.true.
       end do
       do i=1,fragment2%nunoccAOS
          unocc_list(fragment2%unoccAOSidx(i))=.true.
       end do

    end if


    ! Construct pair fragment as union of input fragments
    ! ***************************************************
    ! Simply set MyAtom=0 for pairs
    call atomic_fragment_init_orbital_specific(0,nunocc, nocc, unocc_list, &
         & occ_list,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,pairfragment,DoBasis)


    ! Set remaining information special for pair
    ! ******************************************
    PairFragment%atomic_number = fragment1%atomic_number ! atom 1 in pair
    PairFragment%atomic_number2 = fragment2%atomic_number ! atom 2 in pair
    ! Distance between original (super) fragments
    PairFragment%pairdist =  pairdist
    call mem_dealloc(EOSatoms)

    ! Print out summary
    if(DECinfo%PL>0) then
       write(DECinfo%output,*)
       write(DECinfo%output,*) 'PAIR CONSTRUCTION SUMMARY'
       write(DECinfo%output,*) '*************************'
       write(DECinfo%output,'(1X,a,20i4)')  'PCS: EOS atoms in fragment 1          :', Fragment1%SF_atomlist
       write(DECinfo%output,'(1X,a,20i4)')  'PCS: EOS atoms in fragment 2          :', Fragment2%SF_atomlist
       write(DECinfo%output,'(1X,a,i4)')    'PCS: Number of orbitals in virt total :', PairFragment%nunoccAOS
       write(DECinfo%output,'(1X,a,i4)')    'PCS: Number of orbitals in occ total  :', PairFragment%noccAOS
       write(DECinfo%output,'(1X,a,i4)')    'PCS: Number of basis functions        :', PairFragment%number_basis
       write(DECinfo%output,'(1X,a,F12.3)') 'PCS: Pair distance (Angstrom)         :', PairDist*bohr_to_angstrom
       write(DECinfo%output,*)
    end if

    if(DECinfo%fragadapt) then 
       call pair_fragment_adapted_transformation_matrices(MyMolecule,nocc,nunocc,&
            & OccOrbitals,UnoccOrbitals, Fragment1,Fragment2,PairFragment)
   end if


    ! Set DECinfo%SF back to the input value
    DECinfo%SF = SFsave


  end subroutine merged_fragment_init




  !> Set transformation matrices between fragment-adapted and local orbital bases for
  !> pair fragment (FragmentPQ%CoccFA andFragmentPQ%CunoccFA).
  !> \author Kasper Kristensen
  !> \date February 2013
  subroutine pair_fragment_adapted_transformation_matrices(MyMolecule,nocctot,nunocctot,&
       & OccOrbitals,UnoccOrbitals,FragmentP,FragmentQ,FragmentPQ)
    implicit none

    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule    
    !> Number of occupied orbitals in full molecule
    integer, intent(in) :: nocctot
    !> Number of unoccupied orbitals in full molecule
    integer, intent(in) :: nunocctot
    !> Information about DEC occupied orbitals for full molecule
    type(ccorbital), dimension(nOcctot), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals for full molecule
    type(ccorbital), dimension(nUnocctot), intent(in) :: UnoccOrbitals
    !> Fragment P
    type(ccatom),intent(inout) :: fragmentP
    !> Fragment Q
    type(ccatom),intent(inout) :: fragmentQ
    !> Pair fragment PQ (it is assumed that everything except fragment-adapted
    !> stuff has already been initialized).
    type(ccatom),intent(inout) :: FragmentPQ
    real(realk) :: lambdathr
    real(realk),pointer :: CoccPQ(:,:), CunoccPQ(:,:),aoS(:,:),moS(:,:),CoccPQ_tmp(:,:),SinvEOS(:,:)
    real(realk),pointer :: T(:,:),lambda(:),MOtmp(:,:),tmp(:,:),tmp2(:,:),CEOS(:,:),CunoccPQ_tmp(:,:)
    logical,pointer :: whichOrbitals(:)
    integer :: noccPQ,nunoccPQ,nbasisPQ,noccPQ_FA, nunoccPQ_FA,i,mu,idx, EOSidx,j
    logical :: debugprint,keepon  ! temporary debug prints
    real(realk) :: diagdev, nondiagdev

    debugprint=.true.

    ! Dimensions for FA space (remove EOS)
    noccPQ = fragmentP%noccFA + fragmentQ%noccFA &
         & - fragmentP%noccEOS - fragmentQ%noccEOS
    nunoccPQ = fragmentP%nunoccFA + fragmentQ%nunoccFA &
         & - fragmentP%nunoccEOS - fragmentQ%nunoccEOS
    nbasisPQ = fragmentPQ%number_basis

    ! Sanity checks
    if( (.not. fragmentP%FAtransSet) .or. (.not. fragmentQ%FAtransSet) ) then
       call lsquit('pair_fragment_adapted_transformation_matrices: Transformation &
            & matrices for atomic fragments are not set!',-1)
    end if
    if( (fragmentPQ%noccEOS == fragmentPQ%noccAOS) .or. &
         & (fragmentPQ%nunoccEOS == fragmentPQ%nunoccAOS) ) then
       ! Special case: No FA orbitals to consider.
       call pair_fragment_adapted_transformation_matrices_justEOS(fragmentP,fragmentQ,fragmentPQ)
       return
    end if


    ! Set pair PQ MO-coefficients by simply copying P and Q coefficients
    ! ******************************************************************
    call mem_alloc(CoccPQ_tmp,NbasisPQ,noccPQ)
    call mem_alloc(CunoccPQ_tmp,NbasisPQ,nunoccPQ)
    call set_pair_fragment_adapted_redundant_orbitals(MyMolecule,FragmentP,FragmentQ,&
         & FragmentPQ,noccPQ,nunoccPQ,CoccPQ_tmp,CunoccPQ_tmp)
    ! Note: At this stage the PQ orbitals are redundant and not orthogonal.
    !       Furthermore, the EOS orbitals have been taken out of the orbital pool.
    !       However, FA orbitals originally used for fragment P may still contain
    !       components of EOS orbitals on fragment Q (and vice versa).
    !       These EOS components must be projected out.


    ! Calculate AO overlap matrix for pair fragment
    call mem_alloc(aoS,nbasisPQ,nbasisPQ)
    call adjust_square_matrix(MyMolecule%overlap,aoS,fragmentPQ%atoms_idx, &
         & MyMolecule%atom_size,MyMolecule%atom_start,MyMolecule%atom_end, &
         & MyMolecule%nbasis,MyMolecule%natoms,fragmentPQ%number_basis,FragmentPQ%number_atoms)
!    call ii_get_mixed_overlap_full(DECinfo%output,DECinfo%output,&
 !        & fragmentPQ%mylsitem%setting,aoS,nbasisPQ,nbasisPQ,AORdefault,AORdefault)




    ! ******************************************************************
    !                          OCCUPIED SPACE                          *
    ! ******************************************************************

    ! Project out EOS components
    ! ==========================
    ! MO coefficients where orbital components are projected out:
    !
    ! CoccPQ = CoccPQ_tmp - CEOS SinvEOS CEOS^T aoS CoccPQ_tmp
    !
    ! where CoccPQ_tmp are the MO coefficients determined above 
    ! (which currently contain EOS components),
    ! CEOS are the EOS MO coefficients for the pair fragment, 
    ! and SinvEOS is the inverse overlap matrix for EOS orbitals.
    ! (SinvEOS is simply the identity matrix in the special case where the atomic extent
    !  contains all basis functions in the molecule.)


    ! Extract EOS indices
    ! -------------------
    call mem_alloc(CEOS,nbasisPQ,fragmentPQ%noccEOS)
    do i=1,fragmentPQ%noccEOS
       CEOS(:,i) = fragmentPQ%ypo(:,fragmentPQ%idxo(i))
    end do


    ! Get inverse overlap for EOS orbitals
    ! ------------------------------------
    ! EOS overlap: tmp = CEOS^T aoS CEOS
    call mem_alloc(tmp,fragmentPQ%noccEOS,fragmentPQ%noccEOS)
    call dec_simple_basis_transform1(nbasisPQ,fragmentPQ%noccEOS,CEOS,aoS,tmp)

    ! Inverse EOS overlap: SinvEOS = tmp^-1
    call mem_alloc(SinvEOS,fragmentPQ%noccEOS,fragmentPQ%noccEOS)
    call invert_matrix(tmp,SinvEOS,fragmentPQ%noccEOS)
    call mem_dealloc(tmp)


    ! tmp = aoS CoccPQ_tmp
    ! --------------------
    call mem_alloc(tmp,nbasisPQ,noccPQ)
    call dec_simple_dgemm(nbasisPQ,nbasisPQ,noccPQ,aoS,CoccPQ_tmp,tmp,'n','n')    


    ! tmp2 = CEOS^T aoS CoccPQ_tmp
    ! ----------------------------
    call mem_alloc(tmp2,fragmentPQ%noccEOS,noccPQ)
    call dec_simple_dgemm(fragmentPQ%noccEOS,nbasisPQ,noccPQ,CEOS,tmp,tmp2,'t','n')
    call mem_dealloc(tmp)


    ! tmp = SinvEOS CEOS^T aoS CoccPQ_tmp
    ! -----------------------------------
    call mem_alloc(tmp,fragmentPQ%noccEOS,noccPQ)
    call dec_simple_dgemm(fragmentPQ%noccEOS,fragmentPQ%noccEOS,noccPQ,SinvEOS,tmp2,tmp,'n','n')
    call mem_dealloc(tmp2)


    ! tmp2 = CEOS SinvEOS CEOS^T aoS CoccPQ_tmp
    ! -----------------------------------------
    call mem_alloc(tmp2,nbasisPQ,noccPQ)
    call dec_simple_dgemm(nbasisPQ,fragmentPQ%noccEOS,noccPQ,CEOS,tmp,tmp2,'n','n')
    call mem_dealloc(tmp)


    ! CoccPQ = CoccPQ_tmp - CEOS SinvEOS CEOS^T aoS CoccPQ_tmp  
    ! --------------------------------------------------------
    call mem_alloc(CoccPQ,nbasisPQ,noccPQ)
    CoccPQ = CoccPQ_tmp - tmp2
    call mem_dealloc(tmp2)
    call mem_dealloc(CEOS)
    call mem_dealloc(CoccPQ_tmp)
    call mem_dealloc(SinvEOS)


    ! Pair MO overlap matrix
    ! ======================
    call mem_alloc(moS,noccPQ,noccPQ)
    call dec_simple_basis_transform1(nbasisPQ,noccPQ,CoccPQ,aoS,moS)


    ! Diagonalize overlap matrix for redundant pair molecular orbitals
    ! ================================================================

    call mem_alloc(lambda,noccPQ)
    call mem_alloc(T,noccPQ,noccPQ)
    call solve_eigenvalue_problem_unitoverlap(noccPQ,moS,lambda,T)
    lambda=abs(lambda)
    call mem_alloc(whichorbitals,noccPQ)

    ! Determine which orbitals to include
    keepon=.true.
    lambdathr = 1.0e-3_realk
    OccCheck: do while(keepon)
       whichorbitals=.false.
       do i=1,noccPQ
          if(lambda(i) > lambdathr) whichorbitals(i)=.true.
       end do

       ! Make sanity check that number of orbitals is not larger than in local basis
       if(count(whichorbitals)>fragmentPQ%noccAOS - fragmentPQ%noccEOS) then
          lambdathr = lambdathr*10.0_realk ! increase lambda threshold
       else ! done          
          keepon=.false.
       end if

    end do OccCheck
    ! Number of FA orbitals for pair (not including EOS orbitals)
    noccPQ_FA = count(whichorbitals)
       
    if(debugprint) then
       do i=1,noccPQ
          print *, 'lambda occ', i, lambda(i)
       end do
    end if



    ! Get orthogonal non-redundant MOs for pair
    ! =========================================

    ! Orthogonal fragment-adapted molecular orbitals (psi) are given as:
    ! psi(i) = lambda(i)^{-1/2} * sum_{mu} (CoccPQ T)_{mu i} chi(mu)
    !
    ! where chi(mu) is atomic orbital "mu" in pair atomic fragment extent.
    ! 

    ! MOtmp = CoccPQ T
    call mem_alloc(MOtmp,nbasisPQ,noccPQ)
    call dec_simple_dgemm(nbasisPQ,noccPQ,noccPQ,CoccPQ,T,MOtmp,'n','n')

    ! Final MO coefficents (EOS orbitals + FA orbitals) are stored in fragmentPQ%CoccFA and found by:
    ! 1. Copying existing EOS MO coefficients into the first "noccEOS" columns
    ! 2. Copying the non-redundant orthogonal orbitals (defined by whichorbitals) 
    !    in MOtmp into the remaining columns.

    ! Total number of pair orbitals (EOS +FA)
    fragmentPQ%noccFA = fragmentPQ%noccEOS + noccPQ_FA

    call mem_alloc(fragmentPQ%CoccFA,nbasisPQ,fragmentPQ%noccFA)
    fragmentPQ%CoccFA = 0.0_realk

    ! 1. Copy EOS
    ! Note: For FA orbitals we always put EOS orbitals before remaining orbitals,
    !      this is not the case for local orbitals (see fragment_adapted_transformation_matrices)
    do i=1,fragmentPQ%noccEOS
       ! EOS index for local orbitals in AOS list for local orbitals
       EOSidx = fragmentPQ%idxo(i)
       do mu=1,nbasisPQ
          fragmentPQ%CoccFA(mu,i) = fragmentPQ%ypo(mu,EOSidx)
       end do
    end do
    idx= fragmentPQ%noccEOS  ! counter

    ! 2. Copy FA orbitals (also divide my lambda^{-1/2} as shown above)
    do i=1,noccPQ
       if(whichorbitals(i)) then
          idx = idx+1
          do mu=1,nbasisPQ
             fragmentPQ%CoccFA(mu,idx) = MOtmp(mu,i)/sqrt(lambda(i))
          end do
       end if
    end do
    if(idx/=fragmentPQ%noccFA) then
       call lsquit('pair_fragment_adapted_transformation_matrices: &
            & Counter is different from number of occupied FA orbitals!',-1)
    end if


    ! JUST TESTING
    if(debugprint) then
       call mem_dealloc(moS)
       call mem_alloc(moS,fragmentPQ%noccFA,fragmentPQ%noccFA)
       call dec_simple_basis_transform1(nbasisPQ,fragmentPQ%noccFA,&
            &  fragmentPQ%CoccFA,aoS,moS)
       diagdev=0.0_realk
       nondiagdev=0.0_realk
       do j=1,fragmentPQ%noccFA
          diagdev = max(diagdev,abs(moS(j,j)-1.0_realk))
          do i=1,fragmentPQ%noccFA
             if(i/=j) nondiagdev = max(nondiagdev,abs(moS(i,j)))
          end do
       end do
       print *, 'DIAG / NON-DIAG: ', diagdev, nondiagdev
    end if

    call mem_dealloc(MOtmp)
    call mem_dealloc(lambda)
    call mem_dealloc(T)
    call mem_dealloc(whichorbitals)
    call mem_dealloc(moS)
    call mem_dealloc(CoccPQ)





    ! ******************************************************************
    !                         UNOCCUPIED SPACE                         *
    ! ******************************************************************
    ! Do exactly the same as for occupied space (see comments above)

    call mem_alloc(CEOS,nbasisPQ,fragmentPQ%nunoccEOS)
    do i=1,fragmentPQ%nunoccEOS
       CEOS(:,i) = fragmentPQ%ypv(:,fragmentPQ%idxu(i))
    end do
    call mem_alloc(tmp,fragmentPQ%nunoccEOS,fragmentPQ%nunoccEOS)
    call dec_simple_basis_transform1(nbasisPQ,fragmentPQ%nunoccEOS,CEOS,aoS,tmp)
    call mem_alloc(SinvEOS,fragmentPQ%nunoccEOS,fragmentPQ%nunoccEOS)
    call invert_matrix(tmp,SinvEOS,fragmentPQ%nunoccEOS)
    call mem_dealloc(tmp)
    call mem_alloc(tmp,nbasisPQ,nunoccPQ)
    call dec_simple_dgemm(nbasisPQ,nbasisPQ,nunoccPQ,aoS,CunoccPQ_tmp,tmp,'n','n')    
    call mem_alloc(tmp2,fragmentPQ%nunoccEOS,nunoccPQ)
    call dec_simple_dgemm(fragmentPQ%nunoccEOS,nbasisPQ,nunoccPQ,CEOS,tmp,tmp2,'t','n')
    call mem_dealloc(tmp)
    call mem_alloc(tmp,fragmentPQ%nunoccEOS,nunoccPQ)
    call dec_simple_dgemm(fragmentPQ%nunoccEOS,fragmentPQ%nunoccEOS,nunoccPQ,SinvEOS,tmp2,tmp,'n','n')
    call mem_dealloc(tmp2)
    call mem_alloc(tmp2,nbasisPQ,nunoccPQ)
    call dec_simple_dgemm(nbasisPQ,fragmentPQ%nunoccEOS,nunoccPQ,CEOS,tmp,tmp2,'n','n')
    call mem_dealloc(tmp)
    call mem_alloc(CunoccPQ,nbasisPQ,nunoccPQ)
    CunoccPQ = CunoccPQ_tmp - tmp2
    call mem_dealloc(tmp2)
    call mem_dealloc(CEOS)
    call mem_dealloc(CunoccPQ_tmp)
    call mem_dealloc(SinvEOS)


    ! Pair MO overlap matrix
    call mem_alloc(moS,nunoccPQ,nunoccPQ)
    call dec_simple_basis_transform1(nbasisPQ,nunoccPQ,CunoccPQ,aoS,moS)


    ! Diagonalize overlap matrix for redundant pair molecular orbitals
    call mem_alloc(lambda,nunoccPQ)
    call mem_alloc(T,nunoccPQ,nunoccPQ)
    call solve_eigenvalue_problem_unitoverlap(nunoccPQ,moS,lambda,T)
    lambda=abs(lambda)
    call mem_alloc(whichorbitals,nunoccPQ)
    whichorbitals=.false.

    ! Determine which orbitals to include
    keepon=.true.
    lambdathr = 1.0e-3_realk
    UnoccCheck: do while(keepon)
       whichorbitals=.false.
       do i=1,nunoccPQ
          if(lambda(i) > lambdathr) whichorbitals(i)=.true.
       end do

       ! Make sanity check that number of orbitals is not larger than in local basis
       if(count(whichorbitals)>fragmentPQ%nunoccAOS - fragmentPQ%nunoccEOS) then
          lambdathr = lambdathr*10.0_realk ! increase lambda threshold
       else ! done          
          keepon=.false.
       end if

    end do UnoccCheck
    nunoccPQ_FA = count(whichorbitals)

    if(debugprint) then
       print *
       do i=1,nunoccPQ
          if(lambda(i)>lambdathr) whichorbitals(i)=.true.
          print *, 'lambda unocc', i, lambda(i)
       end do
    end if


    ! Get orthogonal non-redundant MOs for pair
    call mem_alloc(MOtmp,nbasisPQ,nunoccPQ)
    call dec_simple_dgemm(nbasisPQ,nunoccPQ,nunoccPQ,CunoccPQ,T,MOtmp,'n','n')
    fragmentPQ%nunoccFA = fragmentPQ%nunoccEOS + nunoccPQ_FA
    call mem_alloc(fragmentPQ%CunoccFA,nbasisPQ,fragmentPQ%nunoccFA)
    fragmentPQ%CunoccFA = 0.0_realk
    do i=1,fragmentPQ%nunoccEOS
       EOSidx = fragmentPQ%idxu(i)
       do mu=1,nbasisPQ
          fragmentPQ%CunoccFA(mu,i) = fragmentPQ%ypv(mu,EOSidx)
       end do
    end do
    idx= fragmentPQ%nunoccEOS  
    do i=1,nunoccPQ
       if(whichorbitals(i)) then
          idx = idx+1
          do mu=1,nbasisPQ
             fragmentPQ%CunoccFA(mu,idx) = MOtmp(mu,i)/sqrt(lambda(i))
          end do
       end if
    end do
    if(idx/=fragmentPQ%nunoccFA) then
       call lsquit('pair_fragment_adapted_transformation_matrices: &
            & Counter is different from number of unoccupied FA orbitals!',-1)
    end if



    ! JUST TESTING
    if(debugprint) then
       call mem_dealloc(moS)
       call mem_alloc(moS,fragmentPQ%nunoccFA,fragmentPQ%nunoccFA)
       call dec_simple_basis_transform1(nbasisPQ,fragmentPQ%nunoccFA,&
            &  fragmentPQ%CunoccFA,aoS,moS)
       diagdev=0.0_realk
       nondiagdev=0.0_realk
       do j=1,fragmentPQ%nunoccFA
          diagdev = max(diagdev,abs(moS(j,j)-1.0_realk))
          do i=1,fragmentPQ%nunoccFA
             if(i/=j) nondiagdev = max(nondiagdev,abs(moS(i,j)))
          end do
       end do
       print *, 'DIAG / NON-DIAG: ', diagdev, nondiagdev
    end if

    call mem_dealloc(MOtmp)
    call mem_dealloc(lambda)
    call mem_dealloc(T)
    call mem_dealloc(whichorbitals)
    call mem_dealloc(moS)
    call mem_dealloc(CunoccPQ)
    call mem_dealloc(aoS)


    ! Transformation matrices set
    fragmentPQ%FAtransSet=.true.

  end subroutine pair_fragment_adapted_transformation_matrices



  !> Special case for pair_fragment_adapted_transformation_matrices where
  !> there are only EOS orbitals.
  !> ONLY relevant for small test systems.
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine pair_fragment_adapted_transformation_matrices_justEOS(fragmentP,&
       & fragmentQ, fragmentPQ)
    implicit none

    !> Fragment P
    type(ccatom),intent(inout) :: fragmentP
    !> Fragment Q
    type(ccatom),intent(inout) :: fragmentQ
    !> Pair fragment PQ (it is assumed that everything except fragment-adapted
    !> stuff has already been initialized).
    type(ccatom),intent(inout) :: FragmentPQ

    ! Sanity check
    if( (fragmentPQ%noccEOS /= fragmentPQ%noccAOS) .or. &
         & (fragmentPQ%nunoccEOS /= fragmentPQ%nunoccAOS) ) then
       call lsquit('pair_fragment_adapted_transformation_matrices_justEOS:  &
            & EOS must be identical to AOS! ',-1)
    end if


    ! Only occ EOS orbitals
    fragmentPQ%noccFA = fragmentPQ%noccEOS 
    call mem_alloc(fragmentPQ%CoccFA,fragmentPQ%number_basis,fragmentPQ%noccFA)
    fragmentPQ%CoccFA = fragmentPQ%ypo

    ! Only unocc EOS orbitals
    fragmentPQ%nunoccFA = fragmentPQ%nunoccEOS 
    call mem_alloc(fragmentPQ%CunoccFA,fragmentPQ%number_basis,fragmentPQ%nunoccFA)
    fragmentPQ%CunoccFA = fragmentPQ%ypv

    fragmentPQ%FAtransSet=.true.

  end subroutine pair_fragment_adapted_transformation_matrices_justEOS



  !> \brief Get dimensions for fragment-adapted AOS for pair fragment.
  !> Only to be used for job scheduling, pair fragment is deleted again.
  !> \author Kasper Kristensen
  !> \date February 2013
  subroutine get_fragment_adapted_dimensions_for_pair(Fragment1,Fragment2,nunoccFULL, noccFULL, natoms, &
       & OccOrbitals,UnoccOrbitals,DistanceTable,MyMolecule,mylsitem,&
       & noccFA,nunoccFA,nbasisFA)

    implicit none
    !> Fragment 1 in pair
    type(ccatom),intent(inout) :: fragment1
    !> Fragment 2 in pair
    type(ccatom),intent(inout) :: fragment2
    !> Number of occupied orbitals in full molecule
    integer, intent(in) :: noccFULL
    !> Number of unoccupied orbitals in full molecule
    integer, intent(in) :: nunoccFULL
    !> Number of atoms for full molecule
    integer, intent(in) :: natoms
    !> Information about DEC occupied orbitals
    type(ccorbital), dimension(noccFULL), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(ccorbital), dimension(nunoccFULL), intent(in) :: UnoccOrbitals
    !> Distance table for all atoms in the molecule
    real(realk), intent(in) :: DistanceTable(natoms,natoms)
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> Number of occupied fragment-adapted orbitals for pair fragment
    integer,intent(inout) :: noccFA
    !> Number of unoccupied fragment-adapted orbitals for pair fragment
    integer,intent(inout) :: nunoccFA
    !> Number of basis functions for pair fragment
    integer,intent(inout) :: nbasisFA
    type(ccatom) :: pairfragment


    ! Init pair fragment
    ! ------------------
    call merged_fragment_init(Fragment1,Fragment2,nunoccFULL, noccFULL, natoms, &
       & OccOrbitals,UnoccOrbitals,DistanceTable,MyMolecule,mylsitem,.true.,pairfragment)

    noccFA = pairfragment%noccFA
    nunoccFA = pairfragment%nunoccFA
    nbasisFA = pairfragment%number_basis

    call atomic_fragment_free(pairfragment)

  end subroutine get_fragment_adapted_dimensions_for_pair



  !> \brief Initialize atomic fragment extent info for fragment
  !> (fragment%number_atoms, fragment%number_basis, fragment%atoms_idx)
  !> and orbital info in DEC format (fragment%occAOSorb and fragment%unoccAOSorb)
  !> \author Kasper Kristensen
  !> \date Matrch 2013
  subroutine init_atomic_fragment_extent(OccOrbitals,UnoccOrbitals,MyMolecule,fragment)

    implicit none

    !> FUll molecular info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Information about DEC occupied orbitals
    type(ccorbital), dimension(MyMolecule%numocc), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(ccorbital), dimension(MyMolecule%numvirt), intent(in) :: UnoccOrbitals
    !> Fragment info
    type(ccatom), intent(inout) :: fragment
    integer :: i,j,idx
    logical,pointer :: which_atoms(:)


    ! -- Copy occupied orbitals for total occ space
    call mem_alloc(fragment%occAOSorb,Fragment%noccAOS)
    do j=1,Fragment%noccAOS
       idx = fragment%occAOSidx(j)
       fragment%occAOSorb(j) = orbital_init(idx,OccOrbitals(idx)%centralatom, &
            OccOrbitals(idx)%numberofatoms,OccOrbitals(idx)%atoms)
    end do
    ! --

    ! -- Copy unoccupied orbitals for total unocc space
    call mem_alloc(fragment%unoccAOSorb,Fragment%nunoccAOS)
    do j=1,Fragment%nunoccAOS
       idx = fragment%unoccAOSidx(j)
       fragment%unoccAOSorb(j) = orbital_init(idx,UnoccOrbitals(idx)%centralatom, &
            UnoccOrbitals(idx)%numberofatoms,UnoccOrbitals(idx)%atoms)
    end do
    ! --


    ! Which atoms are included in atomic fragment extent?
    call mem_alloc(which_atoms,MyMolecule%natoms)
    which_atoms=.false.

    do i=1,fragment%nunoccAOS  ! loop over unocc orbitals
       ! Loop over atoms used to span MO "i"
       do j=1,fragment%unoccAOSorb(i)%numberofatoms  
          idx = fragment%unoccAOSorb(i)%atoms(j) ! full space index for atom j in orbital extent
          which_atoms(idx)=.true. ! atom "idx" included in orbital extent for orbital "i"
       end do
    end do

    ! Same for occ orbitals
    do i=1,fragment%noccAOS  
       do j=1,fragment%occAOSorb(i)%numberofatoms  
          idx = fragment%occAOSorb(i)%atoms(j)
          which_atoms(idx)=.true.
       end do
    end do

    ! assign atomic extent
    fragment%number_atoms = count(which_atoms)
    call mem_alloc(fragment%atoms_idx,Fragment%number_atoms)
    idx=0
    do i=1,MyMolecule%natoms
       if(which_atoms(i)) then
          idx=idx+1
          fragment%atoms_idx(idx)=i
       end if
    end do
    call mem_dealloc(which_atoms)

    ! count number of basis functions on selected atoms
    fragment%number_basis = 0
    do i=1,fragment%number_atoms
       fragment%number_basis = fragment%number_basis &
            & + MyMolecule%atom_size(fragment%atoms_idx(i))
    end do

    ! Set basis function indices
    call mem_alloc(fragment%basis_idx,fragment%number_basis)
    idx=0
    do i=1,fragment%number_atoms
       do j=1,MyMolecule%atom_size(fragment%atoms_idx(i))
          idx=idx+1
          ! Basis function index 
          fragment%basis_idx(idx) = MyMolecule%atom_start(fragment%atoms_idx(i)) + j-1
       end do
    end do

    ! Sanity check
    if(idx/=fragment%number_basis) then
       print *, 'Counter,nbasis: ', idx,fragment%number_basis
       call lsquit('init_atomic_fragment_extent: Counter does not equal number &
            & of basis functions!',-1)
    end if

  end subroutine init_atomic_fragment_extent



  !> \brief Adjust molecular basis to orbitals included in the molecular fragment
  !> \author Marcin Ziolkowski and Kasper Kristensen
  !> \param fragment Atomic fragment/atomic pair fragment
  !> \param MyMolecule Full molecule info
  subroutine atomic_fragment_basis(fragment,MyMolecule)

    implicit none
    type(ccatom), intent(inout) :: fragment
    type(fullmolecule), intent(in) :: MyMolecule
    type(array2) :: S,tmp1,tmp2
    real(realk), pointer :: smallS(:,:)
    real(realk), pointer :: correct_vector_moS(:), approximated_orbital(:)
    integer :: i,j,idx,atom_k,nocc,nunocc,nbasis,natoms,k,bas_offset,offset, &
         full_orb_idx
    integer, dimension(2) :: dims, dimsAO, dimsMO
    logical,pointer :: which_atoms(:)


    if(.not. fragment%fragmentadapted) then
       ! allocate C^o(nbasis,occ) C^v(nbasis,unocc)
       call mem_alloc(fragment%ypo, fragment%number_basis,  fragment%noccAOS   )
       call mem_alloc(fragment%ypv, fragment%number_basis,  fragment%nunoccAOS )
       ! Core
       fragment%ypo=0.0E0_realk
       fragment%ypv=0.0E0_realk
    end if

    call mem_alloc(fragment%CoreMO, fragment%number_basis, fragment%ncore)
    fragment%CoreMO=0.0E0_realk

    ! truncate basis to this set of atoms
    nocc = MyMolecule%numocc
    nunocc = MyMolecule%numvirt
    nbasis = MyMolecule%nbasis
    natoms = MyMolecule%natoms
    dimsAO(1)=nbasis
    dimsAO(2)=nbasis

    FitOrbitalsForFragment: if(DECinfo%FitOrbitals) then ! fit orbitals for fragment to exact orbitals

       ! fit orbitals
       dimsMO(1) = nbasis
       dimsMO(2) = nocc
       dims(1) = nocc
       dims(2) = nbasis
       S = array2_init(dims)
       call mem_alloc(correct_vector_moS,fragment%number_basis)
       call mem_alloc(approximated_orbital,fragment%number_basis)
       call mem_alloc(smallS,fragment%number_basis,fragment%number_basis)

       call adjust_square_matrix(MyMolecule%overlap,smallS,fragment%atoms_idx, &
            MyMolecule%atom_size,MyMolecule%atom_start,MyMolecule%atom_end, &
            nbasis,natoms,fragment%number_basis,Fragment%number_atoms)

       ! Fragment YPO
       ! half transform overlap
       tmp1 = array2_init(dimsMO,MyMolecule%ypo)
       tmp2 = array2_init(dimsAO,MyMolecule%overlap)
       call array2_matmul(tmp1,tmp2,S,'T','N',1E0_realk,0E0_realk)

       ! Occ orbitals (only valence if frozen core approx is used)
       SkipOcc: if(.not. fragment%fragmentadapted) then ! these are set inside init_fragment_adapted
          do i=1,fragment%noccAOS

             full_orb_idx = fragment%occAOSidx(i)

             ! for each orbital truncate half transformed overlap orbital
             correct_vector_moS = 0.0E0_realk
             bas_offset = 1
             do k=1,Fragment%number_atoms
                atom_k = fragment%atoms_idx(k)
                correct_vector_moS(bas_offset:bas_offset + MyMolecule%atom_size(atom_k) - 1) = &
                     S%val(full_orb_idx,MyMolecule%atom_start(atom_k):MyMolecule%atom_end(atom_k))
                bas_offset = bas_offset + MyMolecule%atom_size(atom_k)
             end do

             ! fit orbital
             approximated_orbital = 0.0E0_realk
             call solve_linear_equations(smallS,approximated_orbital, &
                  correct_vector_moS,fragment%number_basis)
             fragment%ypo(:,i) = approximated_orbital

          end do

       end if SkipOcc


       ! Core orbitals
       do i=1,fragment%ncore

          full_orb_idx = fragment%coreidx(i)

          ! for each orbital truncate half transformed overlap orbital
          correct_vector_moS = 0.0E0_realk
          bas_offset = 1
          do k=1,Fragment%number_atoms
             atom_k = fragment%atoms_idx(k)
             correct_vector_moS(bas_offset:bas_offset + MyMolecule%atom_size(atom_k) - 1) = &
                  S%val(full_orb_idx,MyMolecule%atom_start(atom_k):MyMolecule%atom_end(atom_k))
             bas_offset = bas_offset + MyMolecule%atom_size(atom_k)
          end do

          ! fit orbital
          approximated_orbital = 0.0E0_realk
          call solve_linear_equations(smallS,approximated_orbital, &
               correct_vector_moS,fragment%number_basis)
          fragment%coreMO(:,i) = approximated_orbital

       end do


       SkipUnocc: if(.not. fragment%fragmentadapted) then ! these are set inside init_fragment_adapted
          call array2_free(tmp1)
          call array2_free(S)
          dimsMO(1) = nbasis
          dimsMO(2) = nunocc
          dims(1) = nunocc
          dims(2) = nbasis
          S = array2_init(dims)
          tmp1 = array2_init(dimsMO,MyMolecule%ypv)


          ! Fragmant YPV
          call array2_matmul(tmp1,tmp2,S,'T','N',1E0_realk,0E0_realk)

          do i=1,fragment%nunoccAOS

             full_orb_idx = fragment%unoccAOSidx(i)

             correct_vector_moS = 0.0E0_realk
             bas_offset = 1
             do k=1,Fragment%number_atoms
                atom_k = fragment%atoms_idx(k)
                correct_vector_moS(bas_offset:bas_offset + MyMolecule%atom_size(atom_k) - 1) = &
                     S%val(full_orb_idx,MyMolecule%atom_start(atom_k):MyMolecule%atom_end(atom_k))
                bas_offset = bas_offset + MyMolecule%atom_size(atom_k)
             end do

             ! fit orbital
             approximated_orbital = 0.0E0_realk
             call solve_linear_equations(smallS,approximated_orbital, &
                  correct_vector_moS,fragment%number_basis)
             fragment%ypv(:,i) = approximated_orbital
          end do

       end if SkipUnocc

       ! remove stuff
       call mem_dealloc(correct_vector_moS)
       call mem_dealloc(approximated_orbital)
       call mem_dealloc(smallS)
       call array2_free(tmp1)
       call array2_free(tmp2)
       call array2_free(S)



    else ! Simply extract the exact MO coefficients for atoms in fragment without fitting

       if(DECinfo%frozencore .or. fragment%fragmentadapted) then
          call lsquit('atomic_fragment_basis: Needs implementation!',-1)
       end if

       ! Fragment YPO
       call adjust_basis_matrix(MyMolecule%ypo,fragment%ypo,fragment%occAOSidx, &
            fragment%atoms_idx,MyMolecule%atom_size,MyMolecule%atom_start,MyMolecule%atom_end,nbasis,nocc,natoms, &
            fragment%number_basis,fragment%noccAOS,Fragment%number_atoms)

       ! Fragment YPV
       call adjust_basis_matrix(MyMolecule%ypv,fragment%ypv,fragment%unoccAOSidx, &
            fragment%atoms_idx,MyMolecule%atom_size,MyMolecule%atom_start,MyMolecule%atom_end,nbasis,nunocc,natoms, &
            fragment%number_basis,fragment%nunoccAOS,Fragment%number_atoms)

    end if FitOrbitalsForFragment


 ! adjust fock matrix in ao basis
 call mem_alloc(fragment%fock,fragment%number_basis,fragment%number_basis)
 fragment%fock=0.0E0_realk
 call adjust_square_matrix(MyMolecule%fock,fragment%fock,fragment%atoms_idx,MyMolecule%atom_size, &
      MyMolecule%atom_start,MyMolecule%atom_end,MyMolecule%nbasis,MyMolecule%natoms, &
      fragment%number_basis,Fragment%number_atoms)


 ! adjust fock matrices in mo basis
 ! --------------------------------

 ! For fragment-adapted orbitals, MO Fock matrix is set elsewhere (init_fragment_adapted)
 SkipFock: if(.not. fragment%fragmentadapted) then

    ! Occ-occ block  (valence-valence for frozen core)
    call mem_alloc(fragment%ppfock,fragment%noccAOS,fragment%noccAOS)
    call dec_simple_basis_transform1(fragment%number_basis,fragment%noccAOS,&
         & fragment%ypo,fragment%fock,fragment%ppfock) 

    ! Virtual-virtual block
    call mem_alloc(fragment%qqfock,fragment%nunoccAOS,fragment%nunoccAOS)
    call dec_simple_basis_transform1(fragment%number_basis,fragment%nunoccAOS,&
         & fragment%ypv,fragment%fock,fragment%qqfock) 

 end if SkipFock

 ! Core-core block
 if(fragment%ncore>0) then
    call mem_alloc(fragment%ccfock,fragment%ncore,fragment%ncore)
    call dec_simple_basis_transform1(fragment%number_basis,fragment%ncore,&
         & fragment%coreMO,fragment%fock,fragment%ccfock) 
 end if
 
end subroutine atomic_fragment_basis




  !> \brief Indices of target orbitals
  !> \author Marcin Ziolkowski
  !> \param fragment Atomic fragment/atomic pair fragment
  subroutine target_indices(fragment)

    implicit none
    type(ccatom), intent(inout) :: fragment
    integer :: i,j

    ! don't do that if no occ assigned
    if(fragment%noccEOS == 0) return

    call mem_alloc(fragment%idxo,fragment%noccEOS)
    do i=1,fragment%noccEOS
       do j=1,fragment%noccAOS
          if(fragment%occAOSidx(j) == fragment%occEOSidx(i)) then
             fragment%idxo(i) = j
             exit
          end if
       end do
    end do

    call mem_alloc(fragment%idxu,fragment%nunoccEOS)
    do i=1,fragment%nunoccEOS
       do j=1,fragment%nunoccAOS
          if(fragment%unoccAOSidx(j) == fragment%unoccEOSidx(i)) then
             fragment%idxu(i) = j
             exit
          end if
       end do
    end do

    return
  end subroutine target_indices



  !> \brief Write info for a given fragment to the common file
  !> atomicfragments.info, which, eventually will contain the information
  !> for all atomic fragments (not super fragments).
  !> Furthermore, the file atomicfragmentsdone.info containing the job list info 
  !> (e.g. which fragments have been calculated at this stage) is updated.
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine add_fragment_to_file(fragment,jobs)
    implicit none
    !> Fragment to add to common fragment file
    type(ccatom),intent(inout) :: fragment
    !> Job list
    type(joblist),intent(in) :: jobs
    character(len=40) :: FileName
    integer :: funit,idx,i
    logical :: file_exist

    ! Sanity check: Job has been done for fragment in question
    do i=1,jobs%njobs
       if(jobs%atom1(i)==fragment%atomic_number) then
          if(.not. jobs%jobsdone(i)) then
             call lsquit('add_fragment_to_file: Fragment calculation is not done1',-1)
          end if
       end if
    end do


    ! Write job list info
    ! *******************

    ! Init stuff
    funit = -1
    FileName='atomicfragmentsdone.info'

    ! Replace existing file
    call lsopen(funit,FileName,'REPLACE','UNFORMATTED')

    ! Write fragment job info to file
    call write_fragment_joblist_to_file(jobs,funit)
    call lsclose(funit,'KEEP')



    ! Write minimum (but necessary) fragment info by appending to existing file
    ! *************************************************************************

    ! Init stuff
    funit = 99
    FileName='atomicfragments.info'

    inquire(file=FileName,exist=file_exist)
    if(file_exist) then ! read file and update it

       ! Need to open without lsopen to be able to append to end of file
       open(funit,FILE=FileName,STATUS='OLD',FORM='UNFORMATTED',POSITION='APPEND')

    else ! create new file
       open(funit,FILE=FileName,STATUS='NEW',FORM='UNFORMATTED')
    end if

    ! Write fragment data to file
    call fragment_write_data(funit,fragment)
    close(funit,STATUS='KEEP')


  end subroutine add_fragment_to_file



  !> \brief Write optimized fragment to file.
  !> Assumes that file unit "wunit" has been opened and does NOT close this file.
  !> \author Kasper Kristensen
  !> \param fragment Atomic fragment/atomic pair fragment
  subroutine fragment_write_data(wunit,fragment)
    implicit none
    !> File unit number to write to
    integer, intent(in) :: wunit
    type(ccatom),intent(inout) :: fragment
    integer :: i


    write(wunit) fragment%atomic_number

    ! Occupied AOS orbitals
    write(wunit) fragment%noccAOS
    write(wunit) fragment%occAOSidx

    ! Unoccupied AOS orbitals
    write(wunit) fragment%nunoccAOS
    write(wunit) fragment%unoccAOSidx


    ! Info used mainly for super fragments
    write(wunit) fragment%SF
    write(wunit) fragment%SF_nfrags
    do i=1,fragment%SF_nfrags
       write(wunit) fragment%SF_atomlist(i)
    end do

    ! Energies
    write(wunit) fragment%energies

    ! Occupied AOS orbitals for reduced fragment of lower accuracy
    write(wunit) fragment%REDnoccAOS
    write(wunit) fragment%REDoccAOSidx

    ! Unoccupied AOS orbitals for reduced fragment of lower accuracy
    write(wunit) fragment%REDnunoccAOS
    write(wunit) fragment%REDunoccAOSidx

    ! Correlation density matrices
    write(wunit) fragment%CDset
    if(fragment%CDset) then
       write(wunit) fragment%occmat
       write(wunit) fragment%virtmat
       write(wunit) fragment%RejectThr
       write(wunit) fragment%noccFA
       write(wunit) fragment%nunoccFA
    end if

  end subroutine fragment_write_data



  !> \brief Read fragment info for the fragments which are done
  !> from file atomicfragments.info. Used for restart.
  !> The simple logical list of which fragments are done is read from
  !> file atomicfragmentsdone.info.
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine restart_fragments_from_file(natoms,MyMolecule,MyLsitem,OccOrbitals,UnoccOrbitals,&
       & DoBasis,fragments,jobs)
    implicit none
    !> Number of atoms in the full molecule
    integer,intent(in) :: natoms
    !> Full molecule info
    type(fullmolecule),intent(in)  :: MyMolecule
    !> LSitem infi
    type(lsitem),intent(inout)        :: MyLsitem
    !> Occupied orbitals in DEC format
    type(ccorbital),intent(in)     :: OccOrbitals(MyMolecule%numocc)
    !> Unoccupied orbitals in DEC format
    type(ccorbital),intent(in)     :: UnoccOrbitals(MyMolecule%numvirt)
    !> Construct Fock matrix and MO coeff for fragments?
    logical,intent(in) :: DoBasis
    !> Atomic fragments. NOTE: Only those fragments specified by bookkeeping will be initialized here.
    type(ccatom), intent(inout),dimension(natoms) :: fragments
    !> Job list
    type(joblist),intent(inout) :: jobs
    character(len=40) :: FileName
    integer :: funit, i, MyAtom, ndone,idx,j
    logical :: file_exist


    ! Read list of finished fragments
    ! *******************************

    ! Init stuff
    funit = -1
    FileName='atomicfragmentsdone.info'

    ! Sanity check
    inquire(file=FileName,exist=file_exist)
    if(.not. file_exist) then ! something wrong
       call lsquit('restart_fragments_from_file: File atomicfragmentsdone.info &
            & does not exist!',-1)
    end if

    call lsopen(funit,FileName,'OLD','UNFORMATTED')

    ! Read job list (includes initalization of pointers in job list)
    call read_fragment_joblist_from_file(jobs,funit)
    call lsclose(funit,'KEEP')

    ! Number of fragments done
    ndone = count(jobs%jobsdone)
    write(DECinfo%output,'(a,i8,a)') 'DEC restart: Reading ', ndone, ' fragments'


    ! Read and initialize fragments specified by bookkeeping
    ! ******************************************************

    ! Init stuff
    funit = -1
    FileName='atomicfragments.info'

    ! Sanity check
    inquire(file=FileName,exist=file_exist)
    if(.not. file_exist) then ! something wrong
       call lsquit('restart_fragments_from_file: File atomicfragments.info &
            & does not exist!',-1)
    end if

    ! Open file
    call lsopen(funit,FileName,'OLD','UNFORMATTED')


    ! Read the fragments which were done
    do i=1,ndone

       ! Atom index for the i'th fragment in atomicfragments.info
       if(DECinfo%convert64to32) then
          call read_64bit_to_32bit(funit,MyAtom)
       elseif(DECinfo%convert32to64) then
          call read_32bit_to_64bit(funit,MyAtom)
       else
          read(funit) MyAtom
       end if
       backspace(funit)  ! Go one step in file again to read properly in fragment_read_data


       ! Sanity check
       ! ------------
       ! Find index in job list corresponding to atomic index just read from file
       idx=0
       FindAtom: do j=1,jobs%njobs
          if(jobs%atom1(j)==MyAtom) then
             idx = j
             exit FindAtom
          end if
       end do FindAtom
       ! Check that we found atom in job list
       if(idx==0) then
          call lsquit('restart_fragments_from_file: Atom not found in job list!',-1)
       end if

       ! Consistency check that atom MyAtom (represented by job list index "idx") is 
       ! indeed in the list of finished jobs.
       if(.not. jobs%jobsdone(idx)) then
          print *, 'MyAtom = ', MyAtom
          print *, 'idx in job list = ', idx
          print *, 'jobsdone=', jobs%jobsdone
          call lsquit('restart_fragments_from_file: Central atom read from file &
               & is not in bookkeeping list',-1)
       end if

       call fragment_read_data(funit,fragments(MyAtom),&
            & OccOrbitals,UnoccOrbitals,MyMolecule,Mylsitem,DoBasis)

    end do

    call lsclose(funit,'KEEP')


  end subroutine restart_fragments_from_file



  !> \brief Read optimized fragment from file.
  !> Assumes that file unit "runit" has been opened and does NOT close file.
  !> \author Kasper Kristensen
  !> \param mylsitem Full molecular lsitem
  !> \param runit File unit number to read from
  !> \param fragment Atomic fragment
  subroutine fragment_read_data(runit,fragment,&
       &OccOrbitals,UnoccOrbitals,MyMolecule,MyLsitem,DoBasis)
    implicit none
    type(fullmolecule),intent(in)  :: MyMolecule
    type(lsitem),intent(inout)        :: MyLsitem
    type(ccorbital),intent(in)     :: OccOrbitals(MyMolecule%numocc)
    type(ccorbital),intent(in)     :: UnoccOrbitals(MyMolecule%numvirt)
    type(ccatom),intent(inout)        :: fragment
    integer, intent(in) :: runit
    !> Construct Fock matrix and MO coeff for fragment?
    logical,intent(in) :: DoBasis
    integer             :: i,MyAtom
    logical,pointer :: Occ_list(:),virt_list(:)
    integer :: noccAOS,nvirtAOS
    integer,pointer :: occAOSidx(:), virtAOSidx(:)


    ! Central atom
    if(DECinfo%convert64to32) then
       call read_64bit_to_32bit(runit,MyAtom)
    elseif(DECinfo%convert32to64) then
       call read_32bit_to_64bit(runit,MyAtom)
    else
       read(runit) MyAtom
    end if

    ! Number of occupied AOS orbitals
    if(DECinfo%convert64to32) then
       call read_64bit_to_32bit(runit,noccAOS)
    elseif(DECinfo%convert32to64) then
       call read_32bit_to_64bit(runit,noccAOS)
    else
       read(runit) noccAOS
    end if

    ! Occupied AOS indices
    call mem_alloc(occAOSidx,noccAOS)
    occAOSidx=0
    if(DECinfo%convert64to32) then
       call read_64bit_to_32bit(runit,noccAOS,occAOSidx)
    elseif(DECinfo%convert32to64) then
       call read_32bit_to_64bit(runit,noccAOS,occAOSidx)
    else
       read(runit) occAOSidx
    end if

    ! Logical vector keeping track of which occupied AOS orbitals are included in fragment
    call mem_alloc(occ_list,MyMolecule%numocc)
    occ_list=.false.
    do i=1,noccAOS
       occ_list(occAOSidx(i)) = .true.
    end do


    ! Number of virtual AOS orbitals
    if(DECinfo%convert64to32) then
       call read_64bit_to_32bit(runit,nvirtAOS)
    elseif(DECinfo%convert32to64) then
       call read_64bit_to_32bit(runit,nvirtAOS)
    else
       read(runit) nvirtAOS
    end if

    ! Virtual AOS indices
    call mem_alloc(virtAOSidx,nvirtAOS)
    virtAOSidx=0
    if(DECinfo%convert64to32) then
       call read_64bit_to_32bit(runit,nvirtAOS,virtAOSidx)
    elseif(DECinfo%convert32to64) then
       call read_32bit_to_64bit(runit,nvirtAOS,virtAOSidx)
    else
       read(runit) virtAOSidx
    end if

    ! Logical vector keeping track of which virtual AOS orbitals are included in fragment
    call mem_alloc(virt_list,MyMolecule%numvirt)
    virt_list=.false.
    do i=1,nvirtAOS
       virt_list(virtAOSidx(i)) = .true.
    end do


    ! Info used mainly for super fragment
    if(DECinfo%convert64to32) then
       call read_64bit_to_32bit(runit,fragment%SF)
    elseif(DECinfo%convert32to64) then
       call read_32bit_to_64bit(runit,fragment%SF)
    else
       read(runit) fragment%SF
    end if

    if(DECinfo%convert64to32) then
       call read_64bit_to_32bit(runit,fragment%SF_nfrags)
    elseif(DECinfo%convert32to64) then
       call read_32bit_to_64bit(runit,fragment%SF_nfrags)
    else
       read(runit) fragment%SF_nfrags
    end if

    call mem_alloc(fragment%SF_atomlist,fragment%SF_nfrags)
    do i=1,fragment%SF_nfrags
       if(DECinfo%convert64to32) then
          call read_64bit_to_32bit(runit,fragment%SF_atomlist(i))
       elseif(DECinfo%convert32to64) then
          call read_32bit_to_64bit(runit,fragment%SF_atomlist(i))
       else
          read(runit) fragment%SF_atomlist(i)
       end if
    end do

    ! Initialize fragment
    call atomic_fragment_init_orbital_specific(MyAtom,MyMolecule%numvirt, MyMolecule%numocc,&
         & virt_list,occ_list,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,fragment,DoBasis)

    ! if DECinfo%use_ccsd_frag is .true., then we can safely initialize the parenthesis_t ccatom type
    ! from the CCSD atomic fragment restart files
    if (DECinfo%ccModel .eq. 4 .and. DECinfo%use_ccsd_frag) then

       ! allocate the parenthesis_t ccatom type
       allocate(fragment%parenthesis_t)

       ! Initialize fragment
       call atomic_fragment_init_orbital_specific(MyAtom,MyMolecule%numvirt,MyMolecule%numocc,&
            & virt_list,occ_list,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,fragment%parenthesis_t,DoBasis)

    end if

    ! Fragment energies (only for the Lagrangian scheme are all three used)
    read(runit) fragment%energies

    if (DECinfo%ccModel .eq. 4 .and. DECinfo%use_ccsd_frag) then
       fragment%parenthesis_t%energies = fragment%energies
    end if


    ! Occupied AOS orbitals for reduced fragment of lower accuracy
    if(DECinfo%convert64to32) then
       call read_64bit_to_32bit(runit,fragment%REDnoccAOS)
    elseif(DECinfo%convert32to64) then
       call read_32bit_to_64bit(runit,fragment%REDnoccAOS)
    else
       read(runit) fragment%REDnoccAOS
    end if
    call mem_dealloc(fragment%REDoccAOSidx)
    call mem_alloc(fragment%REDoccAOSidx,fragment%REDnoccAOS)
    if(DECinfo%convert64to32) then
       call read_64bit_to_32bit(runit,fragment%REDnoccAOS,fragment%REDoccAOSidx)
    elseif(DECinfo%convert32to64) then
       call read_32bit_to_64bit(runit,fragment%REDnoccAOS,fragment%REDoccAOSidx)
    else
       read(runit) fragment%REDoccAOSidx
    end if

    ! Unoccupied AOS orbitals for reduced fragment of lower accuracy
    if(DECinfo%convert64to32) then
       call read_64bit_to_32bit(runit,fragment%REDnunoccAOS)
    elseif(DECinfo%convert32to64) then
       call read_32bit_to_64bit(runit,fragment%REDnunoccAOS)
    else
       read(runit) fragment%REDnunoccAOS
    end if
    call mem_dealloc(fragment%REDunoccAOSidx)
    call mem_alloc(fragment%REDunoccAOSidx,fragment%REDnunoccAOS)
    if(DECinfo%convert64to32) then
       call read_64bit_to_32bit(runit,fragment%REDnunoccAOS,fragment%REDunoccAOSidx)
    elseif(DECinfo%convert32to64) then
       call read_32bit_to_64bit(runit,fragment%REDnunoccAOS,fragment%REDunoccAOSidx)
    else
       read(runit) fragment%REDunoccAOSidx
    end if


    ! Correlation density matrices
    if(DECinfo%convert64to32) then
       call read_64bit_to_32bit(runit,fragment%CDset)
    elseif(DECinfo%convert32to64) then
       call read_32bit_to_64bit(runit,fragment%CDset)
    else
       read(runit) fragment%CDset
    end if

    if(fragment%CDset) then
       call mem_alloc(Fragment%OccMat,Fragment%noccAOS,Fragment%noccAOS)
       call mem_alloc(Fragment%VirtMat,Fragment%nunoccAOS,Fragment%nunoccAOS)
       read(runit) fragment%occmat
       read(runit) fragment%virtmat
       read(runit) fragment%RejectThr

       if(DECinfo%convert64to32) then
          call read_64bit_to_32bit(runit,fragment%noccFA)
          call read_64bit_to_32bit(runit,fragment%nunoccFA)
       elseif(DECinfo%convert32to64) then
          call read_32bit_to_64bit(runit,fragment%noccFA)
          call read_32bit_to_64bit(runit,fragment%nunoccFA)
       else
          read(runit) fragment%noccFA
          read(runit) fragment%nunoccFA
       end if

    end if


    call mem_dealloc(Occ_list)
    call mem_dealloc(virt_list)
    call mem_dealloc(occAOSidx)
    call mem_dealloc(virtAOSidx)

  end subroutine fragment_read_data




  !> \brief Get file name for atomic (super) fragment, depends on partitioning scheme and atomic number:
  !> FileName = prefixMyAtom  (e.g. LagFragment17 for standard fragment 17, see prefixes inside)
  !> \author Kasper Kristensen
  !> \date December 2011
  subroutine get_fragment_name(MyAtom,FileName)
    implicit none
    !> Atomic number for fragment
    integer,intent(in) :: MyAtom
    character(len=40),intent(inout) :: FileName
    character(len=40) :: prefix
    integer :: length


    ! Different prefix depending on scheme and standard vs. super fragments
    if(DECinfo%SF) then ! Super fragment
       prefix = "SF_LagFragment"
    else
       prefix = "LagFragment"
    end if

    ! Filename = prefixMyAtom   (e.g. LagFragment17 for standard fragment number 17)
    FileName=prefix
    length=LEN_TRIM(prefix)

    ! Not very elegant - but quick - solution to annoying Fortran issue
    if( (MyAtom > 0) .and. (MyAtom < 10) ) then
       write(FileName(length+1:length+1),'(i1.1)') MyAtom
    elseif( (MyAtom > 9) .and. (MyAtom < 100) ) then
       write(FileName(length+1:length+2),'(i2.2)') MyAtom
    elseif( (MyAtom > 99) .and. (MyAtom < 1000) ) then
       write(FileName(length+1:length+3),'(i3.3)') MyAtom
    elseif( (MyAtom > 999) .and. (MyAtom < 10000) ) then
       write(FileName(length+1:length+4),'(i4.4)') MyAtom
    elseif( (MyAtom > 9999) .and. (MyAtom < 100000) ) then
       write(FileName(length+1:length+5),'(i5.5)') MyAtom
    else
       if(DECinfo%PL>0) write(DECinfo%output,*) 'MyAtom = ', MyAtom
       call lsquit('get_fragment_name: Illegal atom number!',DECinfo%output)
    end if

  end subroutine get_fragment_name



  !> \brief Transform orbital matrix (dimension nrow X ncol) to atom matrix
  !> (natoms X natoms) which contains the largest (absolute)
  !> element for each set of atoms to which the orbitals were assigned.
  !> E.g. if the orbital matrix is the occupied-virtual overlap matrix S
  !> then the (P,Q)th element of the atom matrix contains
  !> the largest (absolute) value of S, where the occupied orbital was assigned
  !> to P and the unoccupied orbital was assigned to Q.
  !> \author Kasper Kristensen
  !> \date 2010-10
  subroutine get_atom_matrix_from_orbital_matrix(Orbitals1,Orbitals2, &
       & orbital_mat, atom_mat, natoms, dim1, dim2)


    implicit none
    !> Dimension 1 for orbital matrix (number of rows)
    integer, intent(in) :: dim1
    !> Dimension 2 for orbital matrix (number of columns)
    integer, intent(in) :: dim2
    !> Number of atoms in the fragment (which may be the full molecule)
    integer, intent(in) :: natoms
    !> Orbital information for orbitals describing the first index of orbital_mat
    type(ccorbital), dimension(dim1), intent(in) :: Orbitals1
    !> Orbital information for orbitals describing the second index of orbital_mat
    type(ccorbital), dimension(dim2), intent(in) :: Orbitals2
    !> Orbital matrix
    real(realk), dimension(dim1,dim2),intent(in) :: orbital_mat
    !> Atom matrix
    real(realk), dimension(natoms,natoms), intent(inout) :: atom_mat
    real(realk) :: tstart,tend,tmp
    integer :: m,n,atomM,atomN


    ! Initialize stuff
    ! ****************
    call cpu_time(tstart)
    atom_mat = 0E0_realk


    Orbloop1: do n=1,dim2

       ! Find atomN to which orbital n has been assigned
       AtomN = Orbitals2(n)%centralatom

       Orbloop2: do m=1,dim1

          ! Find atomM to which orbital m has been assigned
          AtomM = Orbitals1(m)%centralatom

          tmp = abs(orbital_mat(m,n))

          ! Set element of atom matrix
          ReplaceElement: if( tmp > atom_mat(AtomM,AtomN) ) then
             atom_mat(AtomM,AtomN) = tmp
          end if ReplaceElement

       end do Orbloop2
    end do Orbloop1


    call cpu_time(tend)
    if(DECinfo%PL>0) then
       if(DECinfo%show_time) write(DECinfo%output,'(1X,a,g18.8)') &
            & 'Time used in get_atom_matrix_from_orbital_matrix', tend-tstart
    end if


  end subroutine get_atom_matrix_from_orbital_matrix



  !> \brief Estimate memory consumption before fragment calculation is carried out.
  !> If memory consumption is larger than DECinfo%memory, then
  !> the largest arrays are stored on file (Set: DECinfo%array4OnFile=.true.),
  !> else we keep everything in memory (Set: DECinfo%array4OnFile=.false.).
  !> \author Kasper Kristensen
  !> \date February 2011
  subroutine estimate_memory_consumption(MyFragment,offset_memory)


    implicit none
    !> Original fragment information
    type(ccatom), intent(inout) :: MyFragment
    !> Memory already in use before calling estimate_memory_consumption (in GB)
    real(realk), intent(in) :: offset_memory
    real(realk) :: O,V,A,B,Oeos,Veos
    integer :: Bint,intStep
    real(realk) :: intMEM, solMEM, mem_required,LthetaMEM
    real(realk) :: fragmem,mem_in_use
#ifdef VAR_OMP
    integer, external :: OMP_GET_MAX_THREADS
#endif
    integer :: nthreads

    ! If array4OnFile was specified in the input we skip this routine
    if(DECinfo%array4OnFile_specified) then
       if(DECinfo%PL>0) then
          write(DECinfo%output,*)
          write(DECinfo%output,*) 'Arrays are stored on file as requested in input.'
          write(DECinfo%output,*) 'Therefore memory check in estimate_memory_consumption is skipped.'
          write(DECinfo%output,*)
       end if
       return
    end if

    ! Sanity check: Only implemented for MP2
    if(DECinfo%ccModel/=1) then
       call lsquit('estimate_memory_consumption: &
            & Only implemented for the MP2 model!',-1)
    end if


    ! Init stuff
    ! **********

    ! Number of threads
#ifdef VAR_OMP
    ! LSDALTON compiled with OMP
    nthreads=OMP_GET_MAX_THREADS()
#else
    ! No OMP, set number of threads to one
    nthreads=1
#endif

    ! Fragment orbital space
    O = MyFragment%noccAOS ! Number of occupied orbitals
    V = MyFragment%nunoccAOS ! Number of virtual orbitals
    A = MyFragment%number_basis ! Number of atomic orbitals
    ! Maximum batch dimension in integral subroutines
    Bint = max_batch_dimension(MyFragment%mylsitem,MyFragment%number_basis)
    B=Bint
    Oeos = MyFragment%noccEOS ! Number of occupied EOS orbitals
    Veos = MyFragment%nunoccEOS ! Number of unoccupied EOS orbitals

    ! Memory
    intMEM = 0E0_realk
    solMEM = 0E0_realk
    LthetaMEM = 0E0_realk


    ! Fragment memory use
    ! *******************
    call calculate_fragment_memory(MyFragment,fragmem)


    ! Memory in use now
    ! *****************
    ! Memory in use now is the offset_memory (input) plus the fragment memory plus the
    ! memory associated with the full molecule type.
    mem_in_use = offset_memory + fragmem + DECinfo%fullmolecule_memory


    ! NOTE: Currently the memory estimates are hard-coded for the MP2 energy and
    ! first order property calculations.
    if(DECinfo%first_order) then ! First order properties requested
       call estimate_memory_for_mp2_firstorder(nthreads,O,V,A,B,Oeos,Veos,intMEM,intStep,solMEM)
    else
       call estimate_memory_for_mp2_energy(nthreads,O,V,A,B,intMEM,intStep,solMEM)
    end if


    ! Maximum memory required (also add memory in use now)
    mem_required = max(intMEM,solMEM)
    mem_required = mem_required + mem_in_use

    ! Set DECinfo%array4OnFile=.true. if more memory is required than available
    if(mem_required > DECinfo%memory) then ! arrays on file
       DECinfo%array4OnFile = .true.
    else ! arrays in memory
       DECinfo%array4OnFile = .false.
    end if


    ! Print memory summary
    ! ********************
if(DECinfo%PL>0) then
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '========================================================'
    if(DECinfo%first_order) then
       write(DECinfo%output,*) '  FIRST ORDER PROP: MEMORY ESTIMATE SUMMARY'
    else
       write(DECinfo%output,*) '            ENERGY: MEMORY ESTIMATE SUMMARY'
    end if
    write(DECinfo%output,*) '--------------------------------------------------------'
    write(DECinfo%output,'(1X,a,i9)') 'ME: Number of OMP threads         = ', &
         & nthreads
    write(DECinfo%output,'(1X,a,i9)') 'ME: Number of occupied orbitals   = ', &
         & MyFragment%noccAOS
    write(DECinfo%output,'(1X,a,i9)') 'ME: Number of virtual orbitals    = ', &
         & MyFragment%nunoccAOS
    write(DECinfo%output,'(1X,a,i9)') 'ME: Number of atomic orbitals     = ', &
         & MyFragment%number_basis
    write(DECinfo%output,'(1X,a,i9)') 'ME: Maximum number of batches     = ', Bint
    write(DECinfo%output,'(1X,a,g16.6)') 'ME: Offset memory  (GB)           = ', offset_memory
    write(DECinfo%output,'(1X,a,g16.6)') 'ME: Memory in use, fragment  (GB) = ', fragmem
    write(DECinfo%output,'(1X,a,g16.6)') 'ME: Memory in use, full mol  (GB) = ', &
         & DECinfo%fullmolecule_memory
    write(DECinfo%output,'(1X,a,g16.6,a,i3)') &
         & 'ME: Integral memory required (GB) = ', intMEM, ' in step ', intStep
    write(DECinfo%output,'(1X,a,g16.6)') 'ME: Solver memory required   (GB) = ', solMEM
    write(DECinfo%output,'(1X,a,g16.6)') 'ME: Maximum memory required  (GB) = ', &
         & mem_required
    write(DECinfo%output,'(1X,a,g16.6)') 'ME: Maximum memory available (GB) = ', DECinfo%memory
    write(DECinfo%output,*)
    if(intMEM > solMEM) then
       if(DECinfo%array4OnFile) then
          write(DECinfo%output,*) 'ME: INTEGRAL ROUTINE DEFINES MEMORY CONSUMPTION &
               & --> STORE ON FILE'
       else
          write(DECinfo%output,*) 'ME: INTEGRAL ROUTINE DEFINES MEMORY CONSUMPTION &
               & --> STORE IN MEMORY'
       end if

    else
       if(DECinfo%array4OnFile) then
          write(DECinfo%output,*) 'ME: SOLVER ROUTINE DEFINES MEMORY CONSUMPTION &
               & --> STORE ON FILE'
       else
          write(DECinfo%output,*) 'ME: SOLVER ROUTINE DEFINES MEMORY CONSUMPTION &
               & --> STORE IN MEMORY'
       end if
    end if

    write(DECinfo%output,*) '--------------------------------------------------------'
    write(DECinfo%output,*)
    write(DECinfo%output,*)
 end if

  end subroutine estimate_memory_consumption




  !> \brief Estimate memory used in MP2 calculation for first order properties
  !> \author Kasper Kristensen
  !> \date February 2011
  subroutine estimate_memory_for_mp2_firstorder(nthreads,O,V,A,B,Oeos,Veos,intMEM,intStep,solMEM)

    implicit none
    !> Number of OMP threads
    integer, intent(in) :: nthreads
    !> Number of occupied orbitals (as a real)
    real(realk), intent(in) :: O
    !> Number of virtual orbitals (as a real)
    real(realk), intent(in) :: V
    !> Number of atomic orbitals (as a real)
    real(realk), intent(in) :: A
    !> Maximum batch dimension (as a real)
    real(realk), intent(in) :: B
    !> Number of occupied EOS orbitals (as real)
    real(realk), intent(in) :: Oeos
    !> Number of virtual EOS orbitals (as real)
    real(realk), intent(in) :: Veos
    !> Maximum memory for integrals (in GB)
    real(realk),intent(inout) :: intMEM
    !> Step in integral routine where memory use is greatest
    integer, intent(inout) :: intStep
    !> Maximum memory for solver (in GB)
    real(realk),intent(inout) :: solMEM
    real(realk) :: tmp,GB

    ! Initialize stuff
    ! ****************
    GB = 1.000E9_realk ! 1 GB


    ! Memory for integrals (in GB)
    ! ****************************

    ! In different places of the get_VOVO_integrals_mem routine, different amounts of memory
    ! are allocated. Here we find the maximum.
    ! Roughly speaking, we can talk about 5 different memory comsumptions.

    intMEM = O*V*O*A + Veos*Veos*O*A + V*O*A*B + (A*A*B*B + V*A*B*B + V*O*B*B)*nthreads  ! 1
    intStep=1

    tmp = O*V*O*A + Veos*Veos*O*A + 2E0_realk*V*O*A*B  ! 2
    if(tmp > intMEM) then
       intMEM=tmp
       intStep=2
    end if

    tmp = O*V*O*A + Veos*Veos*O*A + V*O*A*B + O*V*O*B  ! 3
    if(tmp > intMEM) then
       intMEM=tmp
       intStep=3
    end if

    tmp = O*V*O*A + Veos*Veos*O*A + V*O*A*B + Veos*V*O*B ! 4
    if(tmp > intMEM) then
       intMEM=tmp
       intStep=4
    end if

    tmp = 2E0_realk*O*V*O*A + V*O*Oeos*Oeos   ! 5
    if(tmp > intMEM) then
       intMEM=tmp
       intStep=5
    end if

    tmp = O*V*O*A + V*O*Oeos*Oeos + Oeos*O*V*O   ! 6
    if(tmp > intMEM) then
       intMEM=tmp
       intStep=6
    end if

    tmp = O*V*O*A + 2E0_realk*V*O*Oeos*Oeos   ! 7
    if(tmp > intMEM) then
       intMEM=tmp
       intStep=7
    end if

    tmp = O*V*O*A + V*O*Oeos*Oeos + V*O*V*O   ! 8
    if(tmp > intMEM) then
       intMEM=tmp
       intStep=8
    end if

    tmp = 2E0_realk*Veos*Veos*O*A    ! 9
    if(tmp > intMEM) then
       intMEM=tmp
       intStep=9
    end if

    tmp = Veos*Veos*O*A + 2E0_realk*V*O*Veos*Veos    ! 10
    if(tmp > intMEM) then
       intMEM=tmp
       intStep=10
    end if

    tmp =V*O*Oeos*Oeos + V*O*V*O + V*O*Veos*Veos  ! 11
    if(tmp > intMEM) then
       intMEM=tmp
       intStep=11
    end if

    ! Multiply intMEM by realk (8) and divide by GB to get memory in GB
    intMEM = realk*intMEM/GB


    ! Memory for solver (in GB)
    ! *************************
    ! Maximum memory allocated in solver is three arrays of dimensions (V,O,V,O)
    solMEM = 3E0_realk*(V*O*V*O)

    ! Special case, the dominating contribution may be when calculating energy - only relevant for testing
    tmp = V*O*V*O + 2*Oeos*Oeos*V*V + 2*O*O*Veos*Veos
    if(tmp > solMEM) then
       write(DECinfo%output,*) 'WARNING: Max memory consumption for energy contractions!'
       solMEM=tmp
    end if
    solMEM = realk*solMEM/GB


  end subroutine estimate_memory_for_mp2_firstorder



  !> \brief Get list of central atoms for occupied EOS orbitals.
  !> Used for pair fragments to tell which orbitals are assigned to which atoms.
  !> (Useful to avoid double counting).
  !> \author Kasper Kristensen
  !> \date November 2010
  subroutine get_central_atoms_for_occ_EOS(occ_atoms,nocc_eos,Fragment)

    implicit none
    !> Number of occupied orbitals in EOS space
    integer, intent(in) :: nocc_eos
    !> Vector listing central atoms for all occupied EOS orbitals
    integer, dimension(nocc_eos), intent(inout) :: occ_atoms
    !> Fragment info
    type(ccatom),intent(inout) :: Fragment
    integer :: i,ix

    ! Initialize stuff
    ! ****************
    occ_atoms=0

    do i=1,nocc_eos
       ! Index of occupied EOS orbital in total AOS list
       ix = Fragment%idxo(i)
       ! Central atom for occupied EOS orbital
       occ_atoms(i) = Fragment%occAOSorb(ix)%centralatom

       ! Sanity check
       if(occ_atoms(i) < 1) then
          write(DECinfo%output,*) 'get_central_atoms_for_occ_EOS:&
               & Atom index for occupied EOS orbital is invalid!'
          write(DECinfo%output,*) 'Occupied EOS orbital index', i
          write(DECinfo%output,*) 'Atom index', occ_atoms(i)
          call lsquit('get_central_atoms_for_occ_EOS:&
               & Atom index for occupied EOS orbital is invalid!',-1)
       end if

    end do

  end subroutine get_central_atoms_for_occ_EOS



  !> \brief Get list of central atoms for virtual EOS orbitals.
  !> Used for pair fragments to tell which orbitals are assigned to which atoms.
  !> (Useful to avoid double counting).
  !> \author Kasper Kristensen
  !> \date November 2010
  subroutine get_central_atoms_for_virt_EOS(unocc_atoms,nunocc_eos,Fragment)

    implicit none
    !> Number of unoccupied orbitals in EOS space
    integer, intent(in) :: nunocc_eos
    !> Vector listing central atoms for all unoccupied EOS orbitals
    integer, dimension(nunocc_eos), intent(inout) :: unocc_atoms
    !> Fragment info
    type(ccatom),intent(inout) :: Fragment
    integer :: a,ax

    ! Initialize stuff
    ! ****************
    unocc_atoms=0


    ! Set indices
    ! ***********

    do a=1,nunocc_eos
       ! Index of unoccupied EOS orbital in total AOS list
       ax = Fragment%idxu(a)
       ! Central atom for unoccupied EOS orbital
       unocc_atoms(a) = Fragment%unoccAOSorb(ax)%centralatom

       ! Sanity check
       if(unocc_atoms(a) < 1) then
          write(DECinfo%output,*) 'get_central_atoms_for_virt_EOS:&
               & Atom index for virtual EOS orbital is invalid!'
          write(DECinfo%output,*) 'Virtual EOS orbital index', a
          write(DECinfo%output,*) 'Atom index', unocc_atoms(a)
          call lsquit('get_central_atoms_for_virt_EOS:&
               & Atom index for virtual EOS orbital is invalid!',-1)
       end if

    end do

  end subroutine get_central_atoms_for_virt_EOS


  !> Function that checks if fragment info file exists
  !> Returns true if file exists, and false if it does not exist.
  !> Author: Ida-Marie Hoeyvik
  function DoesFragFileExist(MyAtom) result(FileExists)
    implicit none

    integer :: MyAtom
    logical :: FileExists
    character(len=80) :: FileName

    call get_fragment_name(myatom,FileName)
    inquire(file=FileName,exist=FileExists)

  end function DoesFragFileExist



  !> Given an atomic fragment which is one of two atomic fragments which
  !> make up a pair fragment:
  !> Set an integer vector frag_in_pair such that frag_in_pair(i) gives the position of occupied
  !> AOS orbital "i" in the atomic fragment list in the pair fragment occupied AOS list.
  !> \author Kasper Kristensen
  !> \date February 2013
  subroutine atomicfragAOS_in_pairfragAOS_occ(atomicfragment,pairfragment,frag_in_pair)
    implicit none

    !> Atomic fragment (one of the two fragments in the pair)
    type(ccatom),intent(inout) :: atomicfragment
    !> Pair fragment
    type(ccatom),intent(inout) :: pairfragment
    !> For each occupied atomic fragment AOS index, what is the position of the same
    !> index in the pair fragment occupied AOS index list
    integer,intent(inout) :: frag_in_pair(atomicfragment%noccAOS)
    integer :: i,j

    frag_in_pair = 0
    fragloop: do i=1,atomicfragment%noccAOS  

       ! Find index "i" in atomic fragment AOS in the pair AOS list
       do j=1,pairfragment%noccAOS
          if(pairfragment%occAOSidx(j) == atomicfragment%occAOSidx(i)) then
             frag_in_pair(i) = j
             cycle fragloop
          end if
       end do

    end do fragloop

  end subroutine atomicfragAOS_in_pairfragAOS_occ



  !> Same as atomicfragAOS_in_pairfragAOS_occ but for the unocc AOS.
  !> \author Kasper Kristensen
  !> \date February 2013
  subroutine atomicfragAOS_in_pairfragAOS_unocc(atomicfragment,pairfragment,frag_in_pair)
    implicit none

    !> Atomic fragment (one of the two fragments in the pair)
    type(ccatom),intent(inout) :: atomicfragment
    !> Pair fragment
    type(ccatom),intent(inout) :: pairfragment
    !> For each unoccupied atomic fragment AOS index, what is the position of the same
    !> index in the pair fragment unoccupied AOS index list
    integer,intent(inout) :: frag_in_pair(atomicfragment%nunoccAOS)
    integer :: i,j

    frag_in_pair = 0
    fragloop: do i=1,atomicfragment%nunoccAOS  

       ! Find index "i" in atomic fragment AOS in the pair AOS list
       do j=1,pairfragment%nunoccAOS
          if(pairfragment%unoccAOSidx(j) == atomicfragment%unoccAOSidx(i)) then
             frag_in_pair(i) = j
             cycle fragloop
          end if
       end do

    end do fragloop

  end subroutine atomicfragAOS_in_pairfragAOS_unocc



  !> Given an atomic fragment which is one of two atomic fragments which
  !> make up a pair fragment:
  !> Set an integer vector frag_in_pair such that frag_in_pair(i) gives the position of 
  !> atomic orbital "i" in the atomic fragment orbital extent 
  !> in the pair fragment orbital extent.
  !> \author Kasper Kristensen
  !> \date February 2013
  subroutine atomicfragAE_in_pairfragAE(nAF,nPF,AF,PF,frag_in_pair)
    implicit none

    !> Number of basis functions in atomic fragment?
    integer,intent(in) :: nAF
    !> Number of basis functions in pair fragment?
    integer,intent(in) :: nPF
    !> Which orbitals in atomic fragment (see which_orbitals_in_atomic_extents)
    integer,intent(in) :: AF(nAF)
    !> Which orbitals in pair fragment (see which_orbitals_in_atomic_extents)
    integer,intent(in) :: PF(nPF)
    integer,intent(inout) :: frag_in_pair(nAF)
    integer :: i,j

    frag_in_pair = 0
    fragloop: do i=1,nAF

       ! Find index "i" in atomic fragment orbital extent in the pair orbital extent list
       do j=1,nPF
          if(PF(j) == AF(i)) then
             frag_in_pair(i) = j
             cycle fragloop
          end if
       end do
       ! Should never go here -something wrong
       call lsquit('atomicfragAE_in_pairfragAE: Index not found!',-1)

    end do fragloop

  end subroutine atomicfragAE_in_pairfragAE


  !> \brief Update full t1 singles amplitude array with
  !> contribution from singles amplitudes for a given atomic fragment
  !> (NOT pair fragment).
  !> \author Kasper Kristensen
  !> \date January 2012
  subroutine update_full_t1_from_atomic_frag(MyFragment,t1full)
    implicit none
    !> Atomic fragment info
    type(ccatom),intent(inout) :: MyFragment
    !> t1 amplitudes for the full molecule, to be updated (order virtual,occupied)
    type(array2), intent(inout) :: t1full
    integer :: i,a,ix,ax

    ! Sanity check
    if(.not. MyFragment%t1_stored) then
       call lsquit('update_full_t1_from_atomic_frag: &
            & Singles amplitudes for fragment are not stored!',DECinfo%output)
    end if


    do i=1,MyFragment%t1dims(2)    ! loop over virtual orbitals in fragment
       ! Occupied index in full list of orbitals
       ix = MyFragment%t1_occidx(i)

       do a=1,MyFragment%t1dims(1)       ! loop over virtual orbitals in fragment

          ! Virtual index in full list of orbitals
          ax = MyFragment%t1_virtidx(a)

          ! Update full t1 array
          t1full%val(ax,ix) = t1full%val(ax,ix) + MyFragment%t1(a,i)

       end do
    end do



  end subroutine update_full_t1_from_atomic_frag


  !> \brief Update full t1 singles amplitude array with
  !> contribution from singles amplitudes for a given pair fragment
  !> (NOT atomic fragment).
  !> \author Kasper Kristensen
  !> \date January 2012
  subroutine update_full_t1_from_pair_frag(PairFragment,nocc,nvirt,&
       & natoms,dopair,OccOrbitals,VirtOrbitals,t1full)
    implicit none
    !> Pair fragment info
    type(ccatom),intent(inout) :: PairFragment
    !> Number of occupied orbitals in molecule
    integer,intent(in) :: nocc
    !> Number of virtual orbitals in molecule
    integer,intent(in) :: nvirt
    !> Number of atoms in molecule
    integer,intent(in) :: natoms
    !> Logical vector telling which orbital pairs to include
    !> (Needed to avoid double counting, see subroutine which_pairs)
    logical,dimension(natoms,natoms),intent(in) :: dopair
    !> Occupied orbitals in DEC format
    type(ccorbital),dimension(nocc) :: OccOrbitals
    !> Virtual orbitals in DEC format
    type(ccorbital),dimension(nvirt) :: VirtOrbitals
    !> t1 amplitudes for the full molecule to be updated (order virtual,occupied)
    type(array2), intent(inout) :: t1full
    integer :: i,a,ix,ax, atomA,atomI

    ! Sanity check
    if(.not. PairFragment%t1_stored) then
       call lsquit('update_full_t1_from_pair_frag: &
            & Singles amplitudes for fragment are not stored!',DECinfo%output)
    end if


    ! Loop over orbital indices for fragment
    do i=1,PairFragment%t1dims(2)    ! loop over occupied orbitals in fragment
       ! Occupied index in full list of orbitals
       ix = PairFragment%t1_occidx(i)
       ! Atom I to which orbital "ix" is assigned
       atomI = OccOrbitals(ix)%centralatom

       do a=1,PairFragment%t1dims(1)       ! loop over virtual orbitals in fragment
          ! Virtual index in full list of orbitals
          ax = PairFragment%t1_virtidx(a)
          ! Atom A to which orbital "ax" is assigned
          atomA = VirtOrbitals(ax)%centralatom

          ! Update full t1 array
          ! Only if orbitals "ix" and "ax" are assigned to different
          ! super fragments (to avoid double counting).
          if(dopair(atomA,atomI)) then
             t1full%val(ax,ix) = t1full%val(ax,ix) + PairFragment%t1(a,i)
          end if

       end do
    end do



  end subroutine update_full_t1_from_pair_frag




  !> \brief From full set of t1 amplitudes, extract those amplitudes
  !> where the virtual/occupied indices belong to the (virtual AOS, occupied AOS)
  !> for the input fragment.
  !> The information is stored in MyFragment%t1.
  !> \author Kasper Kristensen
  !> \date January 2012
  subroutine get_fragmentt1_AOSAOS_from_full(MyFragment,t1full)
    implicit none
    !> Fragment info (only t1-related information will be modified here)
    type(ccatom), intent(inout) :: MyFragment
    !> Singles amplitudes for full molecular system (stored as virtual,occupied)
    type(array2),intent(in) :: t1full
    integer :: noccfrag,nvirtfrag,noccfull,nvirtfull,i,a,ix,ax


    ! Init stuff
    ! **********
    noccfrag = MyFragment%noccAOS  ! occ dimension, fragment
    nvirtfrag = MyFragment%nunoccAOS ! virt dimension, fragment
    noccfull = t1full%dims(2) ! occ dimension, full molecule
    nvirtfull = t1full%dims(1) ! virt dimension, full molecule

    ! Free t1 stuff (in case old amplitudes are already stored)
    ! *********************************************************
    call free_fragment_t1(MyFragment)

    ! Fragment t1 info
    ! ****************
    Myfragment%t1_stored=.true.   ! singles amplitudes will now be stored
    Myfragment%t1dims(1)=nvirtfrag
    Myfragment%t1dims(2)=noccfrag

    ! Initiate fragment t1
    ! *********************
    call mem_alloc(MyFragment%t1_occidx,noccfrag)
    call mem_alloc(MyFragment%t1_virtidx,nvirtfrag)
    call mem_alloc(MyFragment%t1,nvirtfrag,noccfrag)
    MyFragment%t1_occidx = MyFragment%occAOSidx ! occupied AOS indices
    MyFragment%t1_virtidx = MyFragment%unoccAOSidx ! virtual AOS indices


    ! Set fragment t1 amplitude equal to the corresponding full amplitudes
    ! ********************************************************************
    do i=1,noccfrag
       ix=MyFragment%t1_occidx(i)
       do a=1,nvirtfrag
          ax=MyFragment%t1_virtidx(a)

          ! Amplitude elements
          MyFragment%t1(a,i) = t1full%val(ax,ix)

       end do
    end do

  end subroutine get_fragmentt1_AOSAOS_from_full




  !> \brief Extract specific set of indices from t1 indices in fragment.
  !> It is assumed that the fragment t1 is stored in the form
  !> (virt AOS,occ AOS) when this routine is called.
  !> It is then possible for example to extract (virt EOS, occ AOS) etc.
  !> such that at output the fragment t1 contains only these specific indices.
  !> \author Kasper Kristensen
  !> \date January 2012
  subroutine extract_specific_fragmentt1(MyFragment,occEOS,virtEOS)
    implicit none
    !> Fragment info (only t1-related information will be modified here)
    type(ccatom), intent(inout) :: MyFragment
    !> Reduce occupied indices from AOS to EOS size?
    logical,intent(in) :: occEOS
    !> Reduce virtual indices from AOS to EOS size?
    logical,intent(in) :: virtEOS
    integer :: i,a,ix,ax,noccAOS,nvirtAOS,newocc,newvirt
    real(realk),pointer :: oldt1(:,:)

    ! Init stuff
    ! **********
    noccAOS = MyFragment%noccAOS
    nvirtAOS = MyFragment%nunoccAOS


    ! Sanity check: MyFragment%t1 contains amplitudes stored as (virt AOS, occ AOS)
    ! *****************************************************************************
    if(.not. MyFragment%t1_stored) then
       call lsquit('extract_specific_fragmentt1: t1 amplitudes are not stored!', DECinfo%output)
    end if
    if( (nvirtAOS /= MyFragment%t1dims(1)) .or. &
         & (noccAOS /= MyFragment%t1dims(2)) ) then
       write(DECinfo%output,*) 't1 dimensions  = ', MyFragment%t1dims
       write(DECinfo%output,*) 'virtAOS/occAOS = ', nvirtAOS,noccAOS
       call lsquit('extract_specific_fragmentt1: t1 dimensions do not match fragment!', DECinfo%output)
    end if
    if(size(MyFragment%t1) /= noccAOS*nvirtAOS) then
       write(DECinfo%output,*) 'Expected size = ', noccAOS*nvirtAOS
       write(DECinfo%output,*) 'Actual size  = ', size(MyFragment%t1)
       call lsquit('extract_specific_fragmentt1: Wrong dimension of t1 amplitudes!', DECinfo%output)
    end if


    ! Save existing t1 amplitudes
    ! ***************************
    call mem_alloc(oldt1,nvirtAOS,noccAOS)
    do i=1,noccAOS
       do a=1,nvirtAOS
          oldt1(a,i) = MyFragment%t1(a,i)
       end do
    end do


    ! Free t1 information in fragment structure
    ! *****************************************
    call free_fragment_t1(MyFragment)


    ! Re-initialize index bookkeeping according to input (EOS or AOS)
    ! ***************************************************************

    ! Virtual dimension and indices
    if(virtEOS) then
       newvirt = MyFragment%nunoccEOS ! virtual EOS dimension
       call mem_alloc(MyFragment%t1_virtidx,newvirt)
       MyFragment%t1_virtidx = MyFragment%unoccEOSidx ! EOS indices
    else
       newvirt = MyFragment%nunoccAOS ! virtual AOS dimension
       call mem_alloc(MyFragment%t1_virtidx,newvirt)
       MyFragment%t1_virtidx = MyFragment%unoccAOSidx ! AOS indices
    end if

    ! Occupied dimension and indices
    if(occEOS) then
       newocc = MyFragment%noccEOS ! occupied EOS dimension
       call mem_alloc(MyFragment%t1_occidx,newocc)
       MyFragment%t1_occidx = MyFragment%occEOSidx ! EOS indices
    else
       newocc = MyFragment%noccAOS ! occupied AOS dimension
       call mem_alloc(MyFragment%t1_occidx,newocc)
       MyFragment%t1_occidx = MyFragment%occAOSidx ! AOS indices
    end if

    ! Store new dimensions in fragment structure
    MyFragment%t1dims(1) = newvirt
    MyFragment%t1dims(2) = newocc


    ! Init amplitudes and store amplitude elements requested in input
    ! ***************************************************************
    call mem_alloc(MyFragment%t1,newvirt,newocc)

    do i=1,newocc ! occupied

       ! Use EOS or AOS space for occupied indices?
       if(occEOS) then
          ! ix: position of EOS orbital "i" in AOS list of orbitals
          ix = MyFragment%idxo(i)
       else
          ! Simply use all occupied AOS orbitals
          ix = i
       end if

       do a=1,newvirt ! virtual

          ! Use EOS or AOS space for virtual indices?
          if(virtEOS) then
             ! ax: position of EOS orbital "a" in AOS list of orbitals
             ax = MyFragment%idxu(a)
          else
             ! Simply use all virtual AOS orbitals
             ax = a
          end if

          ! Set amplitude values
          MyFragment%t1(a,i) = oldt1(ax,ix)

       end do
    end do


    ! Free stuff
    call mem_dealloc(oldt1)



  end subroutine extract_specific_fragmentt1


  !> \brief Determine whether a new set of fragment calculations
  !> should be carried out with improved full molecular singles amplitudes.
  !> If the relative difference between the new and old singles amplitudes
  !> is larger than the DECinfo%SinglesThr, then "redo=TRUE"; otherwise "redo=FALSE".
  !> \author Kasper Kristensen
  !> \date January 2012
  subroutine redo_fragment_calculations(t1old,t1new,redo)
    implicit none
    !> Old full molecular singles amplitudes (virt,occ) currently to describe long-range effects
    type(array2),intent(in) :: t1old
    !> Improved full molecular singles amplitudes
    type(array2),intent(in) :: t1new
    !> Redo calculation or not?
    logical,intent(inout) :: redo
    type(array2) ::deltat1
    real(realk) :: deltat1_norm, t1new_norm, t1old_norm, reldiff

    ! Difference between new and old singles amplitudes
    deltat1 = array2_add(1.0E0_realk,t1new,-1.0E0_realk,t1old)
    deltat1_norm = sqrt(deltat1*deltat1)
    t1new_norm = sqrt(t1new*t1new)
    t1old_norm = sqrt(t1old*t1old)
    reldiff = deltat1_norm/t1new_norm   ! relative difference

    ! Print out summary
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) 'Singles amplitude summary'
    write(DECinfo%output,*) '*************************'
    write(DECinfo%output,'(1X,a,g12.5)') 'Norm of new ampl = ', t1new_norm
    write(DECinfo%output,'(1X,a,g12.5)') 'Norm of old ampl = ', t1old_norm
    write(DECinfo%output,'(1X,a,g12.5)') 'Norm of delta-t1 = ', deltat1_norm
    write(DECinfo%output,'(1X,a,g12.5)') 'Relative diff    = ', reldiff
    write(DECinfo%output,'(1X,a,g12.5)') 'Threshold value  = ', DECinfo%SinglesThr

    ! Set redo parameter
    if(reldiff > DECinfo%SinglesThr) then
       write(DECinfo%output,*) 'Relative singles difference is larger than threhold!'
       write(DECinfo%output,*) '--> We redo the super fragment calculations'
       print *, 'Relative singles difference is larger than threhold!'
       print *, '--> We redo the super fragment calculations'
       redo=.true.
    else
       write(DECinfo%output,*) 'Relative singles difference is smaller than threhold!'
       write(DECinfo%output,*) '--> Singles treatment is converged.'
       print *, 'Relative singles difference is smaller than threhold!'
       print *, '--> Singles treatment is converged.'
       redo=.false.
    end if
    write(DECinfo%output,*)

    call array2_free(deltat1)


  end subroutine redo_fragment_calculations



  !> Very rough estimate of sizes of atomic fragments.
  !> Measure of size for a given atom "i":
  !> (Number of orbitals within 4 Angstrom from "i") * (AtomNumber) * (EOS orbitals for "i")
  !> In this way the enviroment (first factor), the nature of the atom
  !> (second factor), and the direct orbital assignments (third factor) become determining factors.
  !> The atoms are sorted according to estimated size, with the largest first.
  subroutine estimate_atomic_fragment_sizes(natoms,nocc,nunocc,DistanceTable,&
       & OccOrbitals, UnoccOrbitals, mylsitem,af_list)

    implicit none

    !> Number of atoms in full molecule
    integer, intent(in) :: nAtoms
    !> Number of occupied orbitals in full molecule
    integer, intent(in) :: nocc
    !> Number of unoccupied orbitals in full molecule
    integer, intent(in) :: nunocc
    !> Distance table for all atoms in the molecule
    real(realk), dimension(natoms,natoms), intent(in) :: DistanceTable
    !> Information about DEC occupied orbitals
    type(ccorbital), dimension(nOcc), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(ccorbital), dimension(nUnocc), intent(in) :: UnoccOrbitals
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> List of atoms, sorted according to estimated size (largest first)
    integer,dimension(natoms),intent(inout) :: af_list
    integer :: i,j, atomtype
    integer,dimension(natoms) :: norb,interactions,SizeMeasure
    real(realk) :: dist


    ! Get number of orbitals assigned to each atom "i"
    ! ************************************************
    ! norb is the number of EOS orbitals
    norb=0
    do i=1,natoms

       do j=1,nocc
          ! If occ orbital "j" is assigned to atom "i", then increase counter for atom "i"
          if(OccOrbitals(j)%CentralAtom==i) norb(i) = norb(i)+1
       end do
       do j=1,nunocc
          ! If unocc orbital "j" is assigned to atom "i", then increase counter for atom "i"
          if(UnoccOrbitals(j)%CentralAtom==i) norb(i) = norb(i)+1
       end do

    end do


    ! For each atom "i" count number of "interaction orbitals" assigned to atoms within 4 Angstrom
    ! ********************************************************************************************
    ! Interactions is a rough estimate of AOS size
    dist = 4.0/bohr_to_angstrom
    interactions=0
    do i=1,natoms
       do j=1,natoms

          ! Atom "i": Count atom j's orbitals if j is within distance
          if(DistanceTable(j,i)<dist) then
             interactions(i) = interactions(i) + norb(j)
          end if

       end do
    end do


    ! Size measure = interaction orbitals * atomtype * EOSorbitals
    ! ************************************************************
    SizeMeasure=0
    do i=1,natoms
       atomtype = MyLsitem%input%molecule%atom(i)%Atomic_number ! atom type for atom "i"
       SizeMeasure(i) = interactions(i)*atomtype*norb(i)
    end do

    ! Sort according to size measure (largest first)
    ! **********************************************
    call integer_inv_sort_with_tracking(SizeMeasure,af_list,natoms)

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '******************************************************************'
    write(DECinfo%output,*) '       ATOMIC FRAGMENTS SORTED ACCORDING TO ESTIMATED SIZES       '
    write(DECinfo%output,*) '******************************************************************'
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)') 'JobNumber  AtomType  AtomIndex  SizeMeasure'
    do i=1,natoms
       j=af_list(i)  ! atom index for job number "i"
       write(DECinfo%output,'(1X,i6,7X,a4,2X,i6,2X,i10)') i, MyLsitem%input%molecule%atom(j)%name, &
            & j, SizeMeasure(i)
    end do
    write(DECinfo%output,*)
    write(DECinfo%output,*)

  end subroutine estimate_atomic_fragment_sizes



  !> \brief Set core orbital info in fragment structure.
  !> \author Kasper Kristensen
  !> \date November 2012
  subroutine set_Core_orbitals_for_fragment(MyMolecule,nocc,OccOrbitals,MyFragment) 
    implicit none
    !> Information about the full molecule
    type(fullmolecule),intent(in) :: MyMolecule
    !> Number of occupied orbitals for full molecule
    integer,intent(in) :: nocc
    !> Occupied orbitals
    type(ccorbital), intent(in) :: OccOrbitals(nocc)
    !> Fragment where core information is to be set
    type(ccatom),intent(inout) :: MyFragment 
    integer :: i,idx,atom
    logical,pointer :: which_atoms(:), which_core_orbitals(:)


    ! Sanity check
    if( (MyFragment%noccAOS==0) .or. (MyFragment%REDnoccAOS==0) ) then
       write(DECinfo%output,'(1X,a,2i7)') 'nocc AOS/ reduced AOS', &
            & MyFragment%noccAOS, MyFragment%REDnoccAOS
       call lsquit('set_CoreVal_orbitals_for_fragment: Occupied AOS is not set!',-1)
    end if


    ! Find list of atoms where one or more occupied AOS orbitals are assigned
    call mem_alloc(which_atoms,MyMolecule%natoms)
    which_atoms=.false.
    do i=1,MyFragment%noccAOS
       idx = MyFragment%occAOSidx(i)  ! orbital index in full list
       ! central atom for orbital "i" in AOS list (which is orbital "idx" in full list)
       atom = OccOrbitals(idx)%centralatom  
       which_atoms(atom)=.true.
    end do


    ! Include core orbitals for atoms which already contribute with one or more valence orbitals
    ! ******************************************************************************************
    call mem_alloc(which_core_orbitals,MyMolecule%ncore)
    which_core_orbitals=.false.
    do i=1,MyMolecule%ncore   ! loop over core orbitals
       idx = OccOrbitals(i)%centralatom
       if(which_atoms(idx)) then   ! atom "idx" is already included in occupied valence AOS
          which_core_orbitals(i)=.true.  
       end if
    end do

    ! Set number of core orbitals for fragment
    MyFragment%ncore = count(which_core_orbitals)

    ! Total number of occupied orbitals
    if(DECinfo%frozencore) then
       ! core + valence (AOS)
       MyFragment%nocctot = MyFragment%ncore + MyFragment%noccAOS
    else 
       ! AOS already contains core orbitals
       MyFragment%nocctot = MyFragment%noccAOS
    end if


    ! Set core orbital indices
    call mem_alloc(MyFragment%coreidx,MyFragment%ncore)
    idx=0
    do i=1,MyMolecule%ncore
       if(which_core_orbitals(i)) then
          idx=idx+1
          MyFragment%coreidx(idx) = i
       end if
    end do


    call mem_dealloc(which_atoms)
    call mem_dealloc(which_core_orbitals)

  end subroutine set_Core_orbitals_for_fragment


  !> Get density matrix from fragment as Cocc Cocc^T 
  !> where Cocc are the occupied MO coefficients for the fragment
  !> (Cocc are only valence orbitals if the frozen core approx is used)
  !> \author Kasper Kristensen
  !> \date November 2012
  subroutine get_density_for_fragment(MyFragment,dens)
    implicit none
    !> Fragment info
    type(ccatom),intent(inout) :: MyFragment
    !> Density = Cocc Cocc^T 
    real(realk),intent(inout),dimension(MyFragment%number_basis,MyFragment%number_basis) :: dens

    call get_density_from_occ_orbitals(MyFragment%number_basis,MyFragment%noccAOS,MyFragment%ypo,dens)

  end subroutine get_density_for_fragment

  !> Get core density matrix from fragment as Ccore Ccore^T 
  !> where Ccore are the core MO coefficients for the fragment
  !> \author Kasper Kristensen
  !> \date November 2012
  subroutine get_coredensity_for_fragment(MyFragment,dens)
    implicit none
    !> Fragment info
    type(ccatom),intent(inout) :: MyFragment
    !> Density = Cocc Cocc^T 
    real(realk),intent(inout),dimension(MyFragment%number_basis,MyFragment%number_basis) :: dens

    call get_density_from_occ_orbitals(MyFragment%number_basis,MyFragment%ncore,MyFragment%coreMO,dens)

  end subroutine get_coredensity_for_fragment



  ! SUPER FRAGMENT STUFF - TO BE MODIFIED ASAP
  ! ******************************************



  !> \brief Create job list for DEC calculations remaining after fragment optimization.
  !> \author Kasper Kristensen
  !> \date January 2013
  subroutine create_dec_joblist(MyMolecule,mylsitem,natoms,nocc,nunocc,&
       &DistanceTable,OccOrbitals,UnoccOrbitals,AtomicFragments,which_superfragments,jobs)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> LS item info for full molecule
    type(lsitem), intent(inout) :: mylsitem
    !> Number of atoms for full molecule
    integer, intent(in) :: natoms
    !> Number of occupied orbitals for full molecule
    integer, intent(in) :: nocc
    !> Number of unoccupied orbitals for full molecule
    integer, intent(in) :: nunocc
    !> Table with atomic distances
    real(realk), intent(in) :: DistanceTable(natoms,natoms)
    !> Occupied orbitals
    type(ccorbital), intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals
    type(ccorbital), intent(in) :: UnoccOrbitals(nunocc)
    !> Atomic fragments
    type(ccatom),dimension(natoms),intent(inout) :: AtomicFragments
    !> which_superfragments(i) is true if atom "i" is central in one of the super fragments
    logical, dimension(natoms),intent(inout) :: which_superfragments
    !> Job list of super fragments listed according to size
    type(joblist),intent(inout) :: jobs

    call create_superfragments(MyMolecule,mylsitem,natoms,nocc,nunocc,&
         &DistanceTable,OccOrbitals,UnoccOrbitals, &
         & AtomicFragments,which_superfragments,jobs)

  end subroutine create_dec_joblist


  !> \brief Create super fragments based on information for optimized standard fragments.
  !> This is done if:
  !> (i)  The size of the super fragment and all its associated pairs are smaller than
  !>      the maximum dimensions for existing fragments.
  !> (ii) The distance between any atoms within a super fragment is below
  !>      the maximum accepted super fragment distance (default: 2.5 Angstrom).
  !> \author Kasper Kristensen
  !> \date August 2011
  subroutine create_superfragments(MyMolecule,mylsitem,natoms,nocc,nunocc,&
       &DistanceTable,OccOrbitals,UnoccOrbitals,AtomicFragments,which_superfragments,jobs)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> LS item info for full molecule
    type(lsitem), intent(inout) :: mylsitem
    !> Number of atoms for full molecule
    integer, intent(in) :: natoms
    !> Number of occupied orbitals for full molecule
    integer, intent(in) :: nocc
    !> Number of unoccupied orbitals for full molecule
    integer, intent(in) :: nunocc
    !> Table with atomic distances
    real(realk), intent(in) :: DistanceTable(natoms,natoms)
    !> Occupied orbitals
    type(ccorbital), intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals
    type(ccorbital), intent(in) :: UnoccOrbitals(nunocc)
    !> Standard fragments
    !> Note 1: Only standard atoms with orbitals assigned are associated with
    !>         a standard fragment.
    !> Note 2: Standard fragments will be deleted and replaced by super fragments
    !> Note 3: Both at input (standard fragments) and at output (super fragments)
    !>         only the information in the ccatom type outside the "expensive box"
    !>         is initiated (see type definition of ccatom).
    type(ccatom),dimension(natoms),intent(inout) :: AtomicFragments
    !> which_superfragments(i) is true if atom "i" is central in one of the super fragments
    logical, dimension(natoms),intent(inout) :: which_superfragments
    !> Job list of super fragments listed according to size
    type(joblist),intent(inout) :: jobs
    real(realk), dimension(nAtoms,nAtoms) :: SortedDistTable
    integer, dimension(nAtoms,nAtoms) :: DistanceTrackMatrix
    integer, dimension(natoms) :: nocc_per_atom, nunocc_per_atom
    logical, dimension(natoms,natoms) :: Superfragment
    logical,pointer :: occAOS(:,:), unoccAOS(:,:), REDoccAOS(:,:), REDunoccAOS(:,:)
    logical,pointer :: occAOSpt(:,:), unoccAOSpt(:,:)
    logical,pointer :: occEOS(:,:), unoccEOS(:,:)
    logical, dimension(nocc) :: occpairAOS, occpairEOS, REDoccpairAOS
    logical, dimension(nunocc) :: UNOCCpairAOS, unoccpairEOS, REDunoccpairAOS
    logical,dimension(natoms) :: tmpsuperEOS
    integer :: i,j,k,l,idx,nsuperfrag,myatom
    integer :: atomA, atomB
    logical :: distance_accepted, size_accepted
    integer, pointer :: atomlist(:)
    integer :: npairs, nsuperpairs, REDmaxocc, REDmaxunocc
    integer :: occtmp, unocctmp, EOStmp,atom, nfrags,tmp,maxoccFA,maxunoccFA
    real(realk) :: dist,avoccFA,avunoccFA
    real(realk) :: pairdist(natoms,natoms)
    integer,pointer :: fragsize(:),fragTrack(:)
    real(realk) :: avOCC, avUNOCC, REDavOCC, REDavUNOCC
    integer :: maxoccAOS, maxunoccAOS, noccAOS, nunoccAOS, noccEOS, nunoccEOS
    integer :: maxoccEOS, maxunoccEOS,noccFA,nunoccFA,nbasisFA
    real(realk) :: memavailable, maxmem, memtmp, tcpu,twall
    integer :: njobs

    call LSTIMER('START',tcpu,twall,DECinfo%output)

    ! Get available memory
    call get_currently_available_memory(MemAvailable)

    ! Sanity check
    if(.not. DECinfo%SF) then
       call lsquit('create_superfragments_orbital_specific: DECinfo%SF not set!', DECinfo%output)
    end if

    write(DECinfo%output,*) 'Creating super fragments...'
    write(DECinfo%output,*)

    ! Init stuff
    call mem_alloc(occAOS,nocc,natoms)
    call mem_alloc(unoccAOS,nunocc,natoms)
    call mem_alloc(occEOS,nocc,natoms)
    call mem_alloc(unoccEOS,nunocc,natoms)
    call mem_alloc(REDoccAOS,nocc,natoms)
    call mem_alloc(REDunoccAOS,nunocc,natoms)
    occAOS=.false.
    unoccAOS=.false.
    occEOS=.false.
    unoccEOS=.false.
    REDoccAOS=.false.
    REDunoccAOS=.false.
    if(DECinfo%ccmodel==4) then
       call mem_alloc(occAOSpt,nocc,natoms)
       call mem_alloc(unoccAOSpt,nunocc,natoms)
       occAOSpt=.false.
       unoccAOSpt=.false.
    end if

    ! Get number of occupied/unoccupied orbitals assigned to each atom
    nocc_per_atom =  get_number_of_orbitals_per_atom(OccOrbitals,nocc,natoms)
    nunocc_per_atom =  get_number_of_orbitals_per_atom(UnoccOrbitals,nunocc,natoms)


    ! Notation to keep track of fragments and superfragments
    ! ******************************************************
    ! Initially the diagonal elements of the logical array superfragment tells which atoms are
    ! associated with a standard atomic fragment, i.e. all atoms with at least one orbital assigned.
    ! Before constructing any SFs, this is usually all non-hydrogenic atoms.
    ! During the subroutine small atomic fragments are absorbed into larger atomic fragments to
    ! create super fragments. This in controlled by the off-diagonal elements of the superfragment array.
    !
    ! Example:
    ! --------
    ! If orbitals are originaly assigned to both atom 4 and atom 5 we set:
    ! superfragment(4,4)=.true.   superfragment(j,4)=.false.   for j/=4
    ! superfragment(5,5)=.true.   superfragment(j,5)=.false.   for j/=5
    !
    ! If fragment 5 is absorbed in fragment 4 we set:
    ! superfragment(4,4)=.true.   superfragment(5,4)=.true.
    ! superfragment(j,5)=.false.  for all j (including j=5)
    !
    !
    ! Similar bookkeeping is used for the orbital spaces:
    ! ---------------------------------------------------
    !
    ! occAOS: Occupied AOS orbital space - dimension (nocc,natoms)
    ! Thus, if occAOS(j,i)=.true. then occupied orbital "j" is included in the (super) fragment for atom "i"
    !
    ! unoccAOS: Unoccupied AOS orbital space - dimension (nunocc,natoms)
    ! Thus, if unoccAOS(j,i)=.true. then unoccupied orbital "j" is included in the (super) fragment for atom "i"
    !
    ! etc.



    ! Set initial fragment information (no SF so far)
    ! ***********************************************
    superfragment=.false.
    nfrags=0
    do i=1,natoms

       if( (nunocc_per_atom(i) /= 0) .or. (nocc_per_atom(i) /= 0) ) then
          nfrags=nfrags+1
          superfragment(i,i)=.true.
       end if

    end do


    ! Read information from standard fragments
    ! ****************************************

    ! We need to turn off the superfragment keyword when handling the standard fragments
    ! It is turned back on after the standard fragment information has been obtained.
    maxoccFA=0
    maxunoccFA=0
    avoccFA=0.0_realk
    avunoccFA=0.0_realk
    maxmem=0E0_realk
    maxoccAOS=0
    maxunoccAOS=0
    REDmaxocc=0
    REDmaxunocc=0
    avOCC=0.0E0_realk
    avUNOCC=0.0E0_realk
    REDavOCC=0.0E0_realk
    REDavUNOCC=0.0E0_realk
    call mem_alloc(fragsize,natoms)
    fragsize=0

    GetStandardFrag: do atom=1,natoms

       if(.not. superfragment(atom,atom)) cycle

       ! Set occupied EOS logical vector
       ! ===============================
       do j=1,AtomicFragments(atom)%noccEOS
          idx=AtomicFragments(atom)%occEOSidx(j)  ! index for occupied EOS orbital
          occEOS(idx,atom) = .true.  ! idx is included in "atom" fragment
       end do

       ! Set unoccupied EOS logical vector
       ! =================================
       do j=1,AtomicFragments(atom)%nunoccEOS
          idx=AtomicFragments(atom)%unoccEOSidx(j)  ! index for unoccupied EOS orbital
          unoccEOS(idx,atom) = .true.  ! idx is included in "atom" fragment
       end do

       ! Set occupied AOS logical vector
       ! ===============================
       do j=1,AtomicFragments(atom)%noccAOS
          idx=AtomicFragments(atom)%occAOSidx(j)  ! index for occupied AOS orbital
          occAOS(idx,atom) = .true.  ! idx is included in "atom" fragment
       end do
       ! Same for (T) fragment
       if(DECinfo%ccmodel==4) then
          do j=1,AtomicFragments(atom)%parenthesis_t%noccAOS
             idx=AtomicFragments(atom)%parenthesis_t%occAOSidx(j)  
             occAOSpt(idx,atom) = .true. 
          end do
       end if

       ! Set unoccupied AOS logical vector
       ! =================================
       do j=1,AtomicFragments(atom)%nunoccAOS
          idx=AtomicFragments(atom)%unoccAOSidx(j)  ! index for unoccupied AOS orbital
          unoccAOS(idx,atom) = .true.  ! idx is included in "atom" fragment
       end do
       ! Same for (T) fragment
       if(DECinfo%ccmodel==4) then
          do j=1,AtomicFragments(atom)%parenthesis_t%nunoccAOS
             idx=AtomicFragments(atom)%parenthesis_t%unoccAOSidx(j)  
             unoccAOSpt(idx,atom) = .true.  
          end do
       end if

       ! Also set occ and unocc AOS vectors for reduced fragment of low accuracy
       ! =======================================================================
       do j=1,AtomicFragments(atom)%REDnoccAOS
          idx=AtomicFragments(atom)%REDoccAOSidx(j)  ! index for occupied AOS orbital
          REDoccAOS(idx,atom) = .true.  ! idx is included in "atom" fragment
       end do
       do j=1,AtomicFragments(atom)%REDnunoccAOS
          idx=AtomicFragments(atom)%REDunoccAOSidx(j)  ! index for unoccupied AOS orbital
          REDunoccAOS(idx,atom) = .true.  ! idx is included in "atom" fragment
       end do


       ! Statistics: Fragment sizes
       ! ==========================

       ! Max occupied AOS (measured in no. orbitals)
       occtmp = count(occAOS(1:nocc,atom))
       if(occtmp > maxoccAOS) maxoccAOS=occtmp
       avOCC = avOCC +real(occtmp)
       if(AtomicFragments(atom)%noccFA > maxoccFA) maxoccFA=AtomicFragments(atom)%noccFA
       avoccFA = avoccFA + real(AtomicFragments(atom)%noccFA)
       
       ! Max unoccupied AOS (measured in no. orbitals)
       unocctmp = count(unoccAOS(1:nunocc,atom))
       if(unocctmp > maxunoccAOS) maxunoccAOS=unocctmp
       avUNOCC = avUNOCC +real(unocctmp)
       if(AtomicFragments(atom)%nunoccFA > maxunoccFA) maxunoccFA=AtomicFragments(atom)%nunoccFA
       avunoccFA = avunoccFA + real(AtomicFragments(atom)%nunoccFA)

       ! Fragment size for "atom" measured as (occ AOS size)*(unocc AOS size)
       fragsize(atom) = occtmp*unocctmp

       ! Info for reduced (lower-accuracy) fragment
       tmp = count(REDoccAOS(1:nocc,atom))
       if(tmp > REDmaxocc) REDmaxocc=tmp
       REDavOCC = REDavOCC +real(tmp)

       tmp = count(REDunoccAOS(1:nunocc,atom))
       if(tmp > REDmaxunocc) REDmaxunocc=tmp
       REDavUNOCC = REDavUNOCC +real(tmp)

       ! Free standard fragment (not for simulated super fragments)
       if(.not. DECinfo%simulateSF) then
          if(DECinfo%ccmodel==4) then  ! also free (T) fragments
             call atomic_fragment_free_simple(AtomicFragments(atom)%parenthesis_t)
             deallocate(AtomicFragments(atom)%parenthesis_t)
             nullify(AtomicFragments(atom)%parenthesis_t)
          end if
          call atomic_fragment_free_simple(AtomicFragments(atom))
       end if

    end do GetStandardFrag


    ! Average AOS sizes
    avOCC = avOCC/real(nfrags)
    avUNOCC = avUNOCC/real(nfrags)
    REDavOCC = REDavOCC/real(nfrags)
    REDavUNOCC = REDavUNOCC/real(nfrags)
    avoccFA = avoccFA/real(nfrags)
    avunoccFA = avunoccFA/real(nfrags)


    ! Sort standard fragments according to size and print out
    ! *******************************************************
    call mem_alloc(fragtrack,natoms)
    call integer_inv_sort_with_tracking(fragsize, fragtrack, natoms)

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)') '***************************************************************&
         &****************'
    write(DECinfo%output,*) '             Standard fragments listed according to total size'
    write(DECinfo%output,'(1X,a)') '***************************************************************&
         &****************'

    if(DECinfo%ccmodel==4) then
       write(DECinfo%output,*) '   Index     Occupied (no. orb)      Virtual (no. orb)   Occ(T)  Vir(T)'
    else
       write(DECinfo%output,*) '   Index     Occupied (no. orb)      Virtual (no. orb)'
    end if

    do i=1,natoms
       myatom = fragtrack(i)
       if(.not. SuperFragment(myatom,myatom)) cycle

       occtmp = count(occAOS(1:nocc,myatom))
       unocctmp = count(unoccAOS(1:nunocc,myatom))
       if(DECinfo%ccmodel==4) then
          write(DECinfo%output,'(1X,i6,10X,i6,17X,i6,10X,2i6)') myatom, occtmp, unocctmp,&
               & count(occAOSpt(1:nocc,myatom)), count(unoccAOSpt(1:nunocc,myatom))
       else
          write(DECinfo%output,'(1X,i6,10X,i6,17X,i6)') myatom, occtmp, unocctmp
       end if

    end do
    write(DECinfo%output,*)

    if(DECinfo%fragadapt) then
       write(DECinfo%output,'(1X,a,i8)')    'STANDARD: MAXIMUM OCCUPIED FA SPACE = ', maxoccFA
       write(DECinfo%output,'(1X,a,F12.3)') 'STANDARD: AVERAGE OCCUPIED FA SPACE = ', avoccFA
       write(DECinfo%output,'(1X,a,i8)')    'STANDARD: MAXIMUM VIRTUAL  FA SPACE = ', maxunoccFA
       write(DECinfo%output,'(1X,a,F12.3)') 'STANDARD: AVERAGE VIRTUAL  FA SPACE = ', avunoccFA
    else
       write(DECinfo%output,'(1X,a,i8)')    'STANDARD: MAXIMUM OCCUPIED SPACE = ', maxoccAOS
       write(DECinfo%output,'(1X,a,F12.3)') 'STANDARD: AVERAGE OCCUPIED SPACE = ', avOCC
       write(DECinfo%output,'(1X,a,i8)')    'STANDARD: MAXIMUM VIRTUAL  SPACE = ', maxunoccAOS
       write(DECinfo%output,'(1X,a,F12.3)') 'STANDARD: AVERAGE VIRTUAL  SPACE = ', avUNOCC
       write(DECinfo%output,'(1X,a,i8)')    'STANDARD: REDUCED MAX OCC  SPACE = ', REDmaxOCC
       write(DECinfo%output,'(1X,a,F12.3)') 'STANDARD: REDUCED AVE OCC  SPACE = ', REDavOCC
       write(DECinfo%output,'(1X,a,i8)')    'STANDARD: REDUCED MAX VIRT SPACE = ', REDmaxUNOCC
       write(DECinfo%output,'(1X,a,F12.3)') 'STANDARD: REDUCED AVE VIRT SPACE = ', REDavUNOCC
    end if
    write(DECinfo%output,*)



    ! Find maximum sizes for standard pair fragments
    ! **********************************************
    maxoccAOS=0
    maxunoccAOS=0
    maxoccEOS=0
    maxunoccEOS=0
    maxoccFA=0
    maxunoccFA=0
    avoccFA=0.0_realk
    avunoccFA=0.0_realk
    avocc=0.0_realk
    avunocc=0.0_realk
    npairs=0
    do i=1,natoms
       if(.not. superfragment(i,i)) cycle

       do j=i+1,natoms

          if(.not. superfragment(j,j)) cycle

          PairDistCheck: if(DistanceTable(i,j)<DECinfo%pair_distance_threshold) then
             ! Pair (i,j) is needed - check its size now
             npairs=npairs+1


             if(DECinfo%fragadapt) then  ! different dimensions

                ! Get fragment-adapted dimensions
                call get_fragment_adapted_dimensions_for_pair(AtomicFragments(i),AtomicFragments(j),&
                     & nunocc, nocc, natoms, OccOrbitals,UnoccOrbitals,DistanceTable,MyMolecule,&
                     & mylsitem,noccFA,nunoccFA,nbasisFA)


                maxoccFA = max(maxoccFA,noccFA)
                maxunoccFA = max(maxunoccFA,nunoccFA)
                avoccFA = avoccFA + real(noccFA)
                avunoccFA = avunoccFA + real(nunoccFA)
                
             else
                
                if(DistanceTable(i,j) < DECinfo%PairReductionDistance) then ! standard pair
                   
                   call get_main_pair_info(nocc,nunocc,occAOS(1:nocc,i),unoccAOS(1:nunocc,i),&
                        & occEOS(1:nocc,i),unoccEOS(1:nunocc,i),&
                        & occAOS(1:nocc,j),unoccAOS(1:nunocc,j),&
                        & occEOS(1:nocc,j),unoccEOS(1:nunocc,j),&
                        & noccAOS,nunoccAOS,noccEOS,nunoccEOS,memtmp)
                   
                else ! reduced pair
                   
                   call get_main_pair_info(nocc,nunocc,REDoccAOS(1:nocc,i),REDunoccAOS(1:nunocc,i),&
                        & occEOS(1:nocc,i),unoccEOS(1:nunocc,i),&
                        & REDoccAOS(1:nocc,j),REDunoccAOS(1:nunocc,j),&
                        & occEOS(1:nocc,j),unoccEOS(1:nunocc,j),&
                        & noccAOS,nunoccAOS,noccEOS,nunoccEOS,memtmp)
                   
                end if

                ! Max memory
                if(DECinfo%ccModel==1) maxmem = max(maxmem,memtmp)

                ! Max sizes
                maxoccAOS = max(maxoccAOS,noccAOS)
                maxunoccAOS = max(maxunoccAOS,nunoccAOS)
                maxoccEOS = max(maxoccEOS,noccEOS)
                maxunoccEOS = max(maxunoccEOS,nunoccEOS)

                ! Av sizes
                avocc = avocc + real(noccAOS)
                avunocc = avunocc + real(nunoccAOS)

             end if

          end if PairDistCheck

       end do
    end do

    if(DECinfo%fragadapt) then

       if(npairs==0) then  ! avoid dividing by zero
          avoccFA = 0.0_realk
          avunoccFA = 0.0_realk
       else
          avoccFA = avoccFA/real(npairs)
          avunoccFA = avunoccFA/real(npairs)
       end if

       write(DECinfo%output,'(1X,a,i8)')    'STANDARD PAIR: MAXIMUM OCCUPIED FA  = ', maxoccFA
       write(DECinfo%output,'(1X,a,F12.3)') 'STANDARD PAIR: AVERAGE OCCUPIED FA  = ', avoccFA
       write(DECinfo%output,'(1X,a,i8)')    'STANDARD PAIR: MAXIMUM VIRTUAL  FA  = ', maxunoccFA
       write(DECinfo%output,'(1X,a,F12.3)') 'STANDARD PAIR: AVERAGE VIRTUAL  FA  = ', avunoccFA

    else

       if(npairs==0) then  ! avoid dividing by zero
          avocc = 0.0_realk
          avunocc = 0.0_realk
       else
          avocc = avocc/real(npairs)
          avunocc = avunocc/real(npairs)
       end if

       write(DECinfo%output,'(a,2i7)') 'STANDARD PAIR: Max occ (AOS,EOS)     ', maxoccAOS, maxoccEOS
       write(DECinfo%output,'(a,2i7)') 'STANDARD PAIR: Max unocc (AOS,EOS)   ', maxunoccAOS, maxunoccEOS
       write(DECinfo%output,'(a,F12.3)') 'STANDARD PAIR: Average occ (AOS)   ', avocc
       write(DECinfo%output,'(a,F12.3)') 'STANDARD PAIR: Average unocc (AOS) ', avunocc
       if(DECinfo%ccModel==1) then
          write(DECinfo%output,'(a,2g12.3)') 'STANDARD PAIR: Max memory, mem available (GB)  ', &
               & maxmem, memavailable
          if(maxmem > memavailable) then
             call lsquit('Maximum memory required is larger than memory available!',DECinfo%output)
          end if
       end if

    end if


    ! Create sorted distance table
    ! ****************************
    SortedDistTable=DistanceTable
    call sort_track(SortedDistTable,DistanceTrackMatrix,nAtoms)

    ! At this stage no superfragments have yet been created but we are now in a position to do that...


    ! If super fragments are just simulated, we skip super fragment construction part
    if(DECinfo%SimulateSF) then
       write(DECinfo%output,*) 'SIMULATING SUPER FRAGMENTS, EFFECTIVELY USING STANDARD FRAGMENTS!'
       goto 100
    end if


    ! *************************************************************************************
    ! *************************************************************************************
    ! *                             CREATE SUPER FRAGMENTS                                *
    ! *************************************************************************************
    ! *************************************************************************************

    ! Merge two atomic fragments into a super fragment if:
    ! (i)  The size of the super fragment and all its associated pairs are smaller than
    !      the maximum dimensions maxoccAOS, maxunoccAOS, and mem (for MP2) determined above.
    ! (ii) The distance between the two atomic fragments is below
    !      the maximum accepted super fragment distance (default: 2.5 Angstrom).
    !      (More generally, this distance criterion must be satisfied for any two atoms
    !       in a super fragment).


    ! Loop over current atomic fragments
    CurrentFrags: do i=1,natoms

       myatom=fragtrack(i)  ! largest atomic fragment first
       if(.not. Superfragment(myatom,myatom)) cycle CurrentFrags
       !       write(DECinfo%output,*)
       !       write(DECinfo%output,*) 'Super fragment analysis for atom: ', myatom

       ! Loop over all atomic fragments to see if any can be absorbed in super fragment "myatom",
       ! starting with the atomic fragment where the central atom is closest to "myatom"
       PossibleToAbsorb: do j=1,natoms

          ! Index "idx" for atom number "j" according to distance from myatom
          idx = DistanceTrackMatrix(j,myatom)

          ! Fragment "idx" has no orbitals assigned or is the same as fragment "myatom" --> No merging
          if( (.not. Superfragment(idx,idx)) .or. (idx==myatom) ) cycle PossibleToAbsorb


          ! Distance between fragments is larger than distance threshold --> No merging
          ! ---------------------------------------------------------------------------
          distance_accepted=.true.
          DistanceCheck: do k=1,natoms
             if(Superfragment(k,myatom)) then ! Is atom k included in super fragment "myatom"?
                do l=1,natoms
                   if(Superfragment(l,idx)) then ! Is atom l included in super fragment "idx"?

                      ! (k,l) distance must be below threshold for any atoms k and l
                      ! in (super) atomic fragments myatom and idx.
                      if(DistanceTable(k,l) > DECinfo%SF_maxdist) then
                         distance_accepted=.false.
                         exit
                      end if

                   end if
                end do
             end if
             if(.not. distance_accepted) exit
          end do DistanceCheck
          if(.not. distance_accepted) cycle PossibleToAbsorb


          ! Distance condition accepted: Next, check size of possible new super fragment
          ! ****************************************************************************

          ! Occ space
          ! ---------
          ! Note: occpairAOS etc. describe the possible new super fragment
          ! which is equivalent to a pair fragment formed between current fragments
          ! "myatom" and "idx" (hence, the "pair" name).
          Call get_logical_pair_vector(nocc,occAOS(1:nocc,myatom),occAOS(1:nocc,idx),occpairAOS)
          Call get_logical_pair_vector(nocc,REDoccAOS(1:nocc,myatom),REDoccAOS(1:nocc,idx),&
               & REDoccpairAOS)
          call get_logical_pair_vector(nocc,occEOS(1:nocc,myatom),occEOS(1:nocc,idx),occpairEOS)
          ! Skip super fragment if occupied AOS of super fragment (myatom,idx) is larger than
          ! the maximum size of standard (single and pair) fragments
          if(count(occpairAOS) > maxoccAOS) cycle PossibleToAbsorb


          ! Unocc space
          ! -----------
          call get_logical_pair_vector(nunocc,unoccAOS(1:nunocc,myatom),&
               & unoccAOS(1:nunocc,idx),unoccpairAOS)
          call get_logical_pair_vector(nunocc,REDunoccAOS(1:nunocc,myatom),&
               & REDunoccAOS(1:nunocc,idx),REDunoccpairAOS)
          call get_logical_pair_vector(nunocc,unoccEOS(1:nunocc,myatom),&
               & unoccEOS(1:nunocc,idx),unoccpairEOS)
          ! Skip super fragment if unoccupied AOS of super fragment (myatom,idx) is larger than
          ! the maximum size for standard (single and pair) fragments
          if(count(unoccpairAOS) > maxunoccAOS) cycle PossibleToAbsorb


          ! EOS atoms
          ! ---------
          call get_logical_pair_vector(natoms,superfragment(1:natoms,myatom),&
               & superfragment(1:natoms,idx), tmpsuperEOS)


          ! Check sizes of all relevant pairs associated with possible new super fragment
          size_accepted=.true.
          kloop: do k=1,natoms

             if(.not. superfragment(k,k) .or. (k==myatom) .or. (k==idx) ) cycle kloop

             ! Distance for possible super pair fragment([myatom,idx], k)
             dist = get_superpair_distance_simple(natoms,&
                  & tmpsuperEOS,superfragment(1:natoms,k), DistanceTable)

             CheckPair: if(dist < DECinfo%pair_distance_threshold) then

                if(dist < DECinfo%PairReductionDistance) then ! standard pair

                   call get_main_pair_info(nocc,nunocc,occpairAOS,unoccpairAOS,&
                        & occpairEOS,unoccpairEOS,&
                        & occAOS(1:nocc,k),unoccAOS(1:nunocc,k),&
                        & occEOS(1:nocc,k),unoccEOS(1:nunocc,k),&
                        & noccAOS,nunoccAOS,noccEOS,nunoccEOS,memtmp)

                else ! reduced pair

                   call get_main_pair_info(nocc,nunocc,REDoccpairAOS,REDunoccpairAOS,&
                        & occpairEOS,unoccpairEOS,&
                        & REDoccAOS(1:nocc,k),REDunoccAOS(1:nunocc,k),&
                        & occEOS(1:nocc,k),unoccEOS(1:nunocc,k),&
                        & noccAOS,nunoccAOS,noccEOS,nunoccEOS,memtmp)

                end if

                ! Is pair larger than largest existing pair? If so, exit k-loop.
                if(noccAOS > maxoccAOS) size_accepted=.false.
                if(nunoccAOS > maxunoccAOS) size_accepted=.false.
                if(DECinfo%ccModel==1) then
                   if(memtmp > maxmem) size_accepted=.false.
                end if

                if(.not. size_accepted) exit kloop

             end if CheckPair

          end do kloop


          ! No super fragments if size is not accepted (cycle)
          if(.not. size_accepted) cycle PossibleToAbsorb


          ! ABSORB FRAGMENT "idx" in FRAGMENT "myatom"
          ! ------------------------------------------
          !          write(DECinfo%output,'(1X,a,i5,a,i5)') 'Absorbing fragment', idx, &
          !               & ' in fragment ', myatom

          ! 1. Absorb central atom
          call absorb_logical_vector(natoms, Superfragment(1:natoms,myatom), Superfragment(1:natoms,idx) )
          ! 2. Absorb occupied EOS vector
          call absorb_logical_vector(nocc, occEOS(1:nocc,myatom), occEOS(1:nocc,idx) )
          ! 3. Absorb virtual EOS vector
          call absorb_logical_vector(nunocc, unoccEOS(1:nunocc,myatom), unoccEOS(1:nunocc,idx) )
          ! 4. Absorb occupied AOS vector
          call absorb_logical_vector(nocc, occAOS(1:nocc,myatom), occAOS(1:nocc,idx) )
          ! 5. Absorb virtual AOS vector
          call absorb_logical_vector(nunocc, unoccAOS(1:nunocc,myatom), unoccAOS(1:nunocc,idx) )
          ! 6. Absorb occupied reduced AOS vector
          call absorb_logical_vector(nocc, REDoccAOS(1:nocc,myatom), REDoccAOS(1:nocc,idx) )
          ! 7. Absorb virtual reduced AOS vector
          call absorb_logical_vector(nunocc, REDunoccAOS(1:nunocc,myatom), REDunoccAOS(1:nunocc,idx) )


       end do PossibleToAbsorb
    end do CurrentFrags



    ! Set output vector for superfragments
    ! ************************************
100 which_superfragments=.false.
    nsuperfrag=0
    do i=1,natoms
       if(SuperFragment(i,i)) then
          which_superfragments(i) = .true.
          nsuperfrag = nsuperfrag + 1
       end if
    end do


    ! Set super fragment output (ccatom type)
    ! ***************************************
    maxoccAOS=0
    maxunoccAOS=0
    REDmaxocc=0
    REDmaxunocc=0
    avOCC=0.0E0_realk
    avUNOCC=0.0E0_realk
    REDavOCC=0.0E0_realk
    REDavUNOCC=0.0E0_realk
    nsuperfrag=0
    fragsize=0
    SetSuperFrag: do atom=1,natoms
       if(.not. which_superfragments(atom)) cycle SetSuperFrag

       InitNewFrag: if(.not. DECinfo%simulateSF) then

          ! Set information required to initialize super fragment
          ! *****************************************************

          ! Count number of EOS atoms that have been merged to create SF
          AtomicFragments(atom)%SF_nfrags = count(superfragment(:,atom))

          ! Set list of original atoms now merged into SF
          call mem_alloc(AtomicFragments(atom)%SF_atomlist,AtomicFragments(atom)%SF_nfrags)
          idx=0
          do j=1,natoms
             if(superfragment(j,atom)) then
                idx = idx+1
                AtomicFragments(atom)%SF_atomlist(idx)=j
             end if
          end do

          ! Initialize super fragment (but not basis info)
! HACK
          call atomic_fragment_init_orbital_specific(Atom,nunocc, nocc, unoccAOS(1:nunocc,atom), &
               & occAOS(1:nocc,atom),OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,&
               & AtomicFragments(atom),.false.)

          if (DECinfo%ccModel==4) then

             allocate(AtomicFragments(atom)%parenthesis_t)

             ! Use same super fragments for (T) fragments as for CCSD in CCSD(T) calculations
             AtomicFragments(atom)%parenthesis_t%SF_nfrags = AtomicFragments(atom)%SF_nfrags
             call mem_alloc(AtomicFragments(atom)%parenthesis_t%SF_atomlist,&
                  & AtomicFragments(atom)%SF_nfrags)
             do j=1,AtomicFragments(atom)%SF_nfrags
                AtomicFragments(atom)%parenthesis_t%SF_atomlist(j) = AtomicFragments(atom)%SF_atomlist(j)
             end do
! HACK
             call atomic_fragment_init_orbital_specific(Atom,nunocc, nocc, unoccAOSpt(1:nunocc,atom), &
                  & occAOSpt(1:nocc,atom),OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,&
                  & AtomicFragments(atom)%parenthesis_t,.false.)

          end if


          ! Set reduced fragment information correctly
          ! ******************************************
          AtomicFragments(atom)%REDnoccAOS = count(REDoccAOS(1:nocc,atom))
          AtomicFragments(atom)%REDnunoccAOS = count(REDunoccAOS(1:nunocc,atom))
          call mem_dealloc(AtomicFragments(atom)%REDoccAOSidx)
          call mem_dealloc(AtomicFragments(atom)%REDunoccAOSidx)
          call mem_alloc(AtomicFragments(atom)%REDoccAOSidx,&
               & AtomicFragments(atom)%REDnoccAOS)
          call mem_alloc(AtomicFragments(atom)%REDunoccAOSidx,&
               & AtomicFragments(atom)%REDnunoccAOS)
          idx=0
          do j=1,nocc  ! occupied AOS indices for reduced superfragment
             if(REDoccAOS(j,atom)) then
                idx=idx+1
                AtomicFragments(atom)%REDoccAOSidx(idx) = j
             end if
          end do

          idx=0
          do j=1,nunocc  ! unoccupied AOS indices for reduced superfragment
             if(REDunoccAOS(j,atom)) then
                idx=idx+1
                AtomicFragments(atom)%REDunoccAOSidx(idx) = j
             end if
          end do

       else ! keep fragment as it is but formally consider it to be a super fragment
          AtomicFragments(atom)%SF = .true.
       end if InitNewFrag


       ! Statistics: Fragment sizes
       ! ==========================
       nsuperfrag=nsuperfrag+1

       ! Max occupied AOS (measured in no. orbitals)
       occtmp = count(occAOS(1:nocc,atom))
       if(occtmp > maxoccAOS) maxoccAOS=occtmp
       avOCC = avOCC +real(occtmp)

       ! Max unoccupied AOS (measured in no. orbitals)
       unocctmp = count(unoccAOS(1:nunocc,atom))
       if(unocctmp > maxunoccAOS) maxunoccAOS=unocctmp
       avUNOCC = avUNOCC +real(unocctmp)

       ! Fragment size for "atom" measured as (occ AOS size)*(unocc AOS size)
       fragsize(atom) = occtmp*unocctmp

       ! Info for reduced (lower-accuracy) fragment
       tmp = count(REDoccAOS(1:nocc,atom))
       if(tmp > REDmaxocc) REDmaxocc=tmp
       REDavOCC = REDavOCC +real(tmp)

       tmp = count(REDunoccAOS(1:nunocc,atom))
       if(tmp > REDmaxunocc) REDmaxunocc=tmp
       REDavUNOCC = REDavUNOCC +real(tmp)


    end do SetSuperFrag


    ! Average sizes
    avOCC = avOCC/real(nsuperfrag)
    avUNOCC = avUNOCC/real(nsuperfrag)
    REDavOCC = REDavOCC/real(nsuperfrag)
    REDavUNOCC = REDavUNOCC/real(nsuperfrag)


    ! Sort list according to size
    call integer_inv_sort_with_tracking(fragsize, fragtrack, natoms)


    ! Print superfragment info to output
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)') '**************************************************************&
         &*******************'
    write(DECinfo%output,*) '                 Super fragments listed according to total size'
    write(DECinfo%output,'(1X,a)') '**************************************************************&
         &*******************'
    write(DECinfo%output,'(5X,a)') 'Index   Occ(orb)   Virt(orb)    EOS(atoms)  --  Atoms in EOS'

    do i=1,natoms
       myatom = fragtrack(i)
       if(.not. SuperFragment(myatom,myatom)) cycle

       occtmp = count(occAOS(1:nocc,myatom))
       unocctmp = count(unoccAOS(1:nunocc,myatom))
       EOStmp = count(superfragment(1:natoms,myatom))


       ! Set list of original atoms now merged into SF
       call mem_alloc(atomlist,EOStmp)
       idx=0
       do j=1,natoms
          if(superfragment(j,myatom)) then
             idx = idx+1
             atomlist(idx)=j
          end if
       end do

       write(DECinfo%output,'(2X,i6,5X,i6,5X,i6,6X,i6,3X,3X,a,20i6)') myatom, occtmp, unocctmp, &
            & EOStmp, ' -- ', atomlist

       call mem_dealloc(atomlist)
    end do
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a,i8)')    'SUPER: MAXIMUM OCCUPIED SPACE = ', maxoccAOS
    write(DECinfo%output,'(1X,a,F12.3)') 'SUPER: AVERAGE OCCUPIED SPACE = ', avOCC
    write(DECinfo%output,'(1X,a,i8)')    'SUPER: MAXIMUM VIRTUAL  SPACE = ', maxunoccAOS
    write(DECinfo%output,'(1X,a,F12.3)') 'SUPER: AVERAGE VIRTUAL  SPACE = ', avUNOCC
    write(DECinfo%output,'(1X,a,i8)')    'SUPER: REDUCED MAX OCC  SPACE = ', REDmaxOCC
    write(DECinfo%output,'(1X,a,F12.3)') 'SUPER: REDUCED AVE OCC  SPACE = ', REDavOCC
    write(DECinfo%output,'(1X,a,i8)')    'SUPER: REDUCED VIRTUAL  SPACE = ', REDmaxUNOCC
    write(DECinfo%output,'(1X,a,F12.3)') 'SUPER: REDUCED VIRTUAL  SPACE = ', REDavUNOCC


    write(DECinfo%output,*)




    ! Find maximum sizes for super pair fragments
    ! **********************************************
    maxoccAOS=0
    maxunoccAOS=0
    maxoccEOS=0
    maxunoccEOS=0
    do i=1,natoms
       if(.not. superfragment(i,i)) cycle

       do j=i+1,natoms

          if(.not. superfragment(j,j)) cycle

          dist = get_superpair_distance_simple(natoms,&
               & superfragment(1:natoms,i),superfragment(1:natoms,j), DistanceTable)

          SuperPairDistCheck: if(dist<DECinfo%pair_distance_threshold) then

             ! Pair (i,j) is needed - check its size now
             if(dist < DECinfo%PairReductionDistance) then ! super pair

                call get_main_pair_info(nocc,nunocc,occAOS(1:nocc,i),unoccAOS(1:nunocc,i),&
                     & occEOS(1:nocc,i),unoccEOS(1:nunocc,i),&
                     & occAOS(1:nocc,j),unoccAOS(1:nunocc,j),&
                     & occEOS(1:nocc,j),unoccEOS(1:nunocc,j),&
                     & noccAOS,nunoccAOS,noccEOS,nunoccEOS,memtmp)

             else ! reduced pair

                call get_main_pair_info(nocc,nunocc,REDoccAOS(1:nocc,i),REDunoccAOS(1:nunocc,i),&
                     & occEOS(1:nocc,i),unoccEOS(1:nunocc,i),&
                     & REDoccAOS(1:nocc,j),REDunoccAOS(1:nunocc,j),&
                     & occEOS(1:nocc,j),unoccEOS(1:nunocc,j),&
                     & noccAOS,nunoccAOS,noccEOS,nunoccEOS,memtmp)

             end if

             ! Max memory
             if(DECinfo%ccModel==1) maxmem = max(maxmem,memtmp)

             ! Max sizes
             maxoccAOS = max(maxoccAOS,noccAOS)
             maxunoccAOS = max(maxunoccAOS,nunoccAOS)
             maxoccEOS = max(maxoccEOS,noccEOS)
             maxunoccEOS = max(maxunoccEOS,nunoccEOS)

          end if SuperPairDistCheck

       end do
    end do

    write(DECinfo%output,'(a,2i7)') 'SUPER PAIR: Max occ (AOS,EOS)  ', maxoccAOS, maxoccEOS
    write(DECinfo%output,'(a,2i7)') 'SUPER PAIR: Max unocc (AOS,EOS)', maxunoccAOS, maxunoccEOS
    if(DECinfo%ccModel==1) then
       write(DECinfo%output,'(a,2g12.3)') 'SUPER PAIR: Max memory, mem available (GB)  ', &
            & maxmem, memavailable
       if(maxmem > memavailable) then
          call lsquit('Super pair: Maximum memory required is larger than memory available!',DECinfo%output)
       end if
    end if




    ! Single fragment list - for easy grepping
    ! ****************************************
    if(DECinfo%PL>0) then
       write(DECinfo%output,*)
       write(DECinfo%output,*)
       write(DECinfo%output,*) 'SINGLE SUPER FRAGMENT LIST'
       write(DECinfo%output,*) '**************************'
       write(DECinfo%output,*)
       do atomA=1,natoms
          if(which_superfragments(atomA)) then
             write(DECinfo%output,*) 'SF SINGLE', atomA
          end if
       end do
    end if
    call mem_dealloc(fragsize)
    call mem_dealloc(fragtrack)

    ! Pair fragment list - for easy grepping
    ! **************************************



    ! Find distances between super fragments:
    ! =======================================
    ! Super fragments distance between SF1 and SF2 =
    ! Minimum distance between an EOS atom in SF1 and an EOS atom in SF2
    pairdist=0E0_realk

    do k=1,natoms

       if(.not. Superfragment(k,k)) cycle
       do l=k+1,natoms
          if(.not. Superfragment(l,l)) cycle

          pairdist(k,l)= get_superpair_distance_simple(natoms,&
               & superfragment(1:natoms,k), superfragment(1:natoms,l), DistanceTable)
          pairdist(l,k) = pairdist(k,l)

       end do
    end do


    if(DECinfo%PL>0) then
       write(DECinfo%output,*)
       write(DECinfo%output,*)
       write(DECinfo%output,*) 'PAIR SUPER FRAGMENT LIST'
       write(DECinfo%output,*) '************************'
       write(DECinfo%output,*)
    end if
    nsuperpairs=0
    do atomA=1,natoms
       do atomB=atomA+1,natoms

          if(which_superfragments(atomA) .and. which_superfragments(atomB)) then
             ! AtomA and atomB are both superfragments

             if(pairdist(AtomA,AtomB) < DECinfo%pair_distance_threshold) then
                ! Distance between SF A and B is smaller than pair cut-off threshold

                if(DECinfo%PL>0) then
                   write(DECinfo%output,'(1X,a,F12.3,2i5)') 'SF PAIR (distance[Ang],atom1,atom2)', &
                        & bohr_to_angstrom*pairdist(AtomA,AtomB), atomA,atomB
                end if
                nsuperpairs=nsuperpairs+1

             end if

          end if

       end do
    end do

    ! Set job list for super fragments
    ! ---------------------------------
    if(DECinfo%first_order .or. (.not. DECinfo%simulateSF) .or. DECinfo%InclFullMolecule &
         & .or. count(which_superfragments)==1 ) then
       ! First order calculation or actual super fragment calculation requires
       ! atomic fragment calculations to be repeated
       njobs = nsuperfrag + nsuperpairs
    else
       ! For simple energy calculation with no super fragments it is not necessary to repeat atomic
       ! fragment calculations so we only do the pairs
       njobs = nsuperpairs
    end if
    call init_joblist(njobs,jobs)

    call set_superfragment_joblist(natoms,nocc,nunocc,occAOS,unoccAOS,&
         & REDoccAOS,REDunoccAOS,superfragment, DistanceTable, MyMolecule, &
         & OccOrbitals,UnoccOrbitals,jobs)

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '*****************************************************'
    write(DECinfo%output,*) '*              SUPER FRAGMENT JOB LIST              *'
    write(DECinfo%output,*) '*****************************************************'
    write(DECinfo%output,*) 'Number of jobs = ', njobs
    write(DECinfo%output,*)
    write(DECinfo%output,*) 'JobIndex            Jobsize         Atom(s) involved '
    do i=1,njobs
       if(jobs%atom1(i)==jobs%atom2(i)) then ! single
          write(DECinfo%output,'(1X,i8,4X,i15,7X,i8)') i,jobs%jobsize(i),jobs%atom1(i)
       else ! pair
          write(DECinfo%output,'(1X,i8,4X,i15,7X,2i8)') i,jobs%jobsize(i),jobs%atom1(i),jobs%atom2(i)
       end if
    end do
    write(DECinfo%output,*)
    write(DECinfo%output,*)

    ! Summary print out
    write(DECinfo%output,'(1X,a,2i8,F10.2,i8)') 'SF SUMMARY (single superfrags, single frags, ratio)', &
         & nsuperfrag, nfrags, real(nsuperfrag)/real(nfrags)
    if(npairs>0) then
       write(DECinfo%output,'(1X,a,2i8,F10.2,i8)') 'SF SUMMARY (pair superfrags,   pair frags,   ratio)', &
            & nsuperpairs, npairs, real(nsuperpairs)/real(npairs)
    end if
    write(DECinfo%output,*)
    write(DECinfo%output,*)

    call mem_dealloc(occAOS)
    call mem_dealloc(unoccAOS)
    call mem_dealloc(occEOS)
    call mem_dealloc(unoccEOS)
    call mem_dealloc(REDoccAOS)
    call mem_dealloc(REDunoccAOS)
    if(DECinfo%ccmodel==4) then
       call mem_dealloc(occAOSpt)
       call mem_dealloc(unoccAOSpt)
    end if
    call LSTIMER('GET SUPERFRAGS',tcpu,twall,DECinfo%output)


  end subroutine create_superfragments


  ! \brief Get main info about size for pair fragment from sets of logical arrays
  !> without constructing pair fragment explicitly.
  !> \author Kasper Kristensen
  !> \date March 2012
  subroutine get_main_pair_info(nocc,nunocc,occAOS1,unoccAOS1,occEOS1,unoccEOS1,&
       & occAOS2,unoccAOS2,occEOS2,unoccEOS2,noccAOS,nunoccAOS,noccEOS,nunoccEOS,mem)

    implicit none
    !> Number of occupied orbitals for full molecule
    integer, intent(in) :: nocc
    !> Number of unoccupied orbitals for full molecule
    integer, intent(in) :: nunocc
    !> Logical vector giving occupied AOS orbitals for fragment 1
    logical,intent(in) :: occAOS1(nocc)
    !> Logical vector giving unoccupied AOS orbitals for fragment 1
    logical,intent(in) :: unoccAOS1(nunocc)
    !> Logical vector giving unoccupied EOS orbitals for fragment 1
    logical,intent(in) :: occEOS1(nocc)
    !> Logical vector giving unoccupied EOS orbitals for fragment 1
    logical,intent(in) :: unoccEOS1(nunocc)
    !> Logical vector giving occupied AOS orbitals for fragment 2
    logical,intent(in) :: occAOS2(nocc)
    !> Logical vector giving unoccupied AOS orbitals for fragment 2
    logical,intent(in) :: unoccAOS2(nunocc)
    !> Logical vector giving unoccupied EOS orbitals for fragment 2
    logical,intent(in) :: occEOS2(nocc)
    !> Logical vector giving unoccupied EOS orbitals for fragment 2
    logical,intent(in) :: unoccEOS2(nunocc)
    !> Number of occupied AOS orbitals in pair fragment
    integer,intent(inout) :: noccAOS
    !> Number of unoccupied AOS orbitals in pair fragment
    integer,intent(inout) :: nunoccAOS
    !> Number of occupied EOS orbitals in pair fragment
    integer,intent(inout) :: noccEOS
    !> Number of unoccupied EOS orbitals in pair fragment
    integer,intent(inout) :: nunoccEOS
    !> KK fix me: Currently returns zero, will soon be redundant!
    real(realk),intent(inout) :: mem
    logical, dimension(nocc) :: occpairAOS, occpairEOS
    logical, dimension(nunocc) :: unoccpairAOS, unoccpairEOS


    ! Merge occupied AOS for fragment 1 and 2
    call get_logical_pair_vector(nocc,occAOS1,occAOS2,occpairAOS)
    noccAOS = count(occpairAOS)

    ! Merge unoccupied AOS for fragment 1 and 2
    call get_logical_pair_vector(nunocc,unoccAOS1,unoccAOS2,unoccpairAOS)
    nunoccAOS = count(unoccpairAOS)

    ! Merge occupied EOS for fragment 1 and 2
    call get_logical_pair_vector(nocc,occEOS1,occEOS2,occpairEOS)
    noccEOS = count(occpairEOS)

    ! Merge unoccupied EOS for fragment 1 and 2
    call get_logical_pair_vector(nunocc,unoccEOS1,unoccEOS2,unoccpairEOS)
    nunoccEOS = count(unoccpairEOS)

    ! Memory used (fix me - generalize to CC2 and CCSD)
    mem=0E0_realk

end subroutine get_main_pair_info



  !> \brief Set superfragment job list. The jobs are listed according to size
  !> with the largest jobs first.
  !> Note: MPI fragment statistics is not modified here.
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine set_superfragment_joblist(natoms,nocc,nunocc,occAOS,unoccAOS,&
       & REDoccAOS,REDunoccAOS,superfragment, DistanceTable, MyMolecule, &
       & OccOrbitals,UnoccOrbitals, jobs)

    implicit none
    !> Number of atoms
    integer,intent(in) :: natoms
    !> Number of occupied orbitals
    integer,intent(in) :: nocc
    !> Number of unoccupied orbitals
    integer,intent(in) :: nunocc
    !> Logical vector describing occupied AOS (see create_superfragments)
    logical,dimension(nocc,natoms),intent(in) :: occAOS
    !> Logical vector describing unoccupied AOS (see create_superfragments)
    logical,dimension(nunocc,natoms),intent(in) :: unoccAOS
    !> Logical vector describing reduced occupied AOS (see create_superfragments)
    logical,dimension(nocc,natoms),intent(in) :: REDoccAOS
    !> Logical vector describing reduced unoccupied AOS (see create_superfragments)
    logical,dimension(nunocc,natoms),intent(in) :: REDunoccAOS
    !> Logical vector describing super fragment composition (see create_superfragments)
    logical,dimension(natoms,natoms),intent(in) :: superfragment
    !> Distance table with interatomic distances
    real(realk),dimension(natoms,natoms),intent(in) :: DistanceTable
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Occupied orbitals
    type(ccorbital), intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals
    type(ccorbital), intent(in) :: UnoccOrbitals(nunocc)
    !> Job list for super fragments
    type(joblist),intent(inout) :: jobs
    logical :: occpairAOS(nocc), unoccpairAOS(nunocc)
    integer :: i,j,k,njobs,nbasisFragment,nsingle
    real(realk) :: dist,tcpu,twall
    integer,pointer :: atom1(:),atom2(:),order(:)

    call LSTIMER('START',tcpu,twall,DECinfo%output)

    ! Init stuff
    k=0
    njobs = jobs%njobs
    if(njobs < 1) then ! Sanity check
       call lsquit('set_superfragment_joblist : Number of jobs must be positive!',-1)
    end if
    call mem_alloc(atom1,njobs)
    call mem_alloc(atom2,njobs)
    atom1=0
    atom2=0

    nsingle=0
    do i=1,natoms
       if(superfragment(i,i)) nsingle=nsingle+1
    end do

    ! **************************
    ! MAIN LOOP TO GET JOB SIZES
    ! **************************

    do i=1,natoms  ! Loop over atoms

       if(.not. superfragment(i,i)) cycle  ! No superfragment for atom i

       if(DECinfo%first_order .or. (.not. DECinfo%simulateSF) .or. DECinfo%InclFullMolecule &
            & .or. nsingle==1 ) then

          ! First order calculation or actual super fragment calculation requires
          ! atomic fragment calculations to be repeated, while this is not necessary for
          ! simple energy calculation.

          ! Set atom indices for single super fragment job
          k=k+1
          atom1(k) = i
          atom2(k) = i   ! Same index for both atoms to distinguish single from pair jobs

          ! Job size is defined as occupied AOS * unoccupied AOS dimensions * nbasis
          call get_nbasis_for_fragment(nocc,nunocc,occAOS(1:nocc,i),unoccAOS(1:nunocc,i),&
               & OccOrbitals,UnoccOrbitals,MyMolecule,nbasisFragment)
          jobs%jobsize(k) = count(occAOS(1:nocc,i))*count(unoccAOS(1:nunocc,i))*nbasisFragment

       end if

       ! Pair loop
       do j=i+1,natoms
          if(.not. superfragment(j,j)) cycle ! No super fragment for atom j

          ! Distance between atoms i and j
          dist = get_superpair_distance_simple(natoms,&
               & superfragment(1:natoms,i),superfragment(1:natoms,j), DistanceTable)

          CheckPair: if(dist < DECinfo%pair_distance_threshold) then  ! Pair needs to be computed
             k=k+1

             if(dist < DECinfo%PairReductionDistance) then

                ! Merge AOS for fragment 1 and 2 for standard pair
                call get_logical_pair_vector(nocc,occAOS(1:nocc,i),occAOS(1:nocc,j),occpairAOS)
                call get_logical_pair_vector(nunocc,unoccAOS(1:nunocc,i),&
                     &unoccAOS(1:nunocc,j),unoccpairAOS)

             else

                ! Merge AOS for fragment 1 and 2 for reduced pair
                call get_logical_pair_vector(nocc,REDoccAOS(1:nocc,i),REDoccAOS(1:nocc,j),occpairAOS)
                call get_logical_pair_vector(nunocc,REDunoccAOS(1:nunocc,i),&
                     & REDunoccAOS(1:nunocc,j),unoccpairAOS)

             end if

             ! Atomic indices and jobsize for pair
             atom1(k) = i
             atom2(k) = j
             call get_nbasis_for_fragment(nocc,nunocc,occpairAOS,unoccpairAOS,&
                  & OccOrbitals,UnoccOrbitals,MyMolecule,nbasisFragment)
             jobs%jobsize(k) = count(occpairAOS(1:nocc))*count(unoccpairAOS(1:nunocc))*nbasisFragment

             if(jobs%jobsize(k)<1) then
                print *, 'dist', dist
                print *, 'k=',k
                print *, 'pair=', i,j
                print *, 'job size: ', jobs%jobsize(k)
                call lsquit('set_superfragment_joblist: Non-positive job size, something wrong',DECinfo%output)
             end if

          end if CheckPair

       end do
    end do

    ! Sanity check: k must equal number of jobs
    if(k/=njobs) then
       print *, 'k     = ',k
       print *, 'njobs = ', njobs
       call lsquit('set_superfragment_joblist: Something wrong in superfragment bookkeeping',-1)
    end if

    ! Sort job list according to size (largest elements first)
    call mem_alloc(order,njobs)
    call integer_inv_sort_with_tracking(jobs%jobsize,order,njobs)

    ! Arrange atoms in correct order and put into job structure
    do i=1,njobs
       jobs%atom1(i) = atom1(order(i))
       jobs%atom2(i) = atom2(order(i))
    end do

    ! No jobs have been done
    jobs%jobsdone=.false.

    ! Clean up
    call mem_dealloc(atom1)
    call mem_dealloc(atom2)
    call mem_dealloc(order)

    call LSTIMER('SF JOBLIST',tcpu,twall,DECinfo%output)


  end subroutine set_superfragment_joblist



  !> \brief Write fragment energies for super fragment for easy restart to
  !> file superenergies.info.
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine write_superfragment_energies_for_restart(natoms,FragEnergies,jobs)

    implicit none
    !> Number of atoms in the molecule
    integer,intent(in) :: natoms
    !> Fragment energies (see ccatom type def)
    real(realk),dimension(natoms,natoms,ndecenergies),intent(in) :: FragEnergies
    !> Job list of super fragments 
    type(joblist),intent(in) :: jobs
    character(len=40) :: FileName
    integer :: funit
    logical :: file_exist

    ! Init stuff
    funit = -1
    FileName='superenergies.info'

    ! Delete existing file
    inquire(file=FileName,exist=file_exist)
    if(file_exist) then  ! backup exisiting file
#ifdef SYS_AIX
       call rename('superenergies.info\0','superenergies.backup\0')
#else
       call rename('superenergies.info','superenergies.backup')
#endif
    end if

    ! Create a new file superenergies.info
    call lsopen(funit,FileName,'NEW','UNFORMATTED')

    ! Write job list info and fragment energies
    call write_fragment_joblist_to_file(jobs,funit)
    write(funit) FragEnergies

    call lsclose(funit,'KEEP')

  end subroutine write_superfragment_energies_for_restart



  !> \brief Read fragment energies for super fragment from file superenergies.info.
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine read_superfragment_energies_for_restart(natoms,FragEnergies,jobs)

    implicit none
    !> Number of atoms in the molecule
    integer,intent(in) :: natoms
    !> Fragment energies (see ccatom type def)
    real(realk),dimension(natoms,natoms,ndecenergies),intent(inout) :: FragEnergies
    !> Job list of super fragments 
    type(joblist),intent(inout) :: jobs
    character(len=40) :: FileName
    integer :: funit
    logical :: file_exist

    ! Init stuff
    funit = -1
    FileName='superenergies.info'

    ! Sanity check
    inquire(file=FileName,exist=file_exist)
    if(.not. file_exist) then
       call lsquit('read_superfragment_energies_for_restart: &
            & File superenergies.info does not exist!',-1)
    end if

    ! Create a new file superenergies.info
    call lsopen(funit,FileName,'OLD','UNFORMATTED')

    ! Read job list and fragment energies from file
    call read_fragment_joblist_from_file(jobs,funit)
    read(funit) FragEnergies
    call lsclose(funit,'KEEP')

  end subroutine read_superfragment_energies_for_restart



  !> \brief Expand superfragment job list to include additional pair fragments
  !> (the additional pairs to include are determined in dec_energy_control_center).
  !> Note: The order of the existing jobs will be unchanged, the new jobs will
  !> simply be appended at the end of the job list.
  !> \author Kasper Kristensen
  !> \data October 2012
  subroutine expand_joblist_to_include_more_pairs(nocc,nunocc,natoms,SuperDistanceTable,dofrag,&
       & Fragments,oldpaircut,newpaircut,MyMolecule,OccOrbitals,UnoccOrbitals,jobs)

    implicit none
    !> Number of occupied orbitals in full molecule
    integer,intent(in) :: nocc
    !> Number of unoccupied orbitals in full molecule
    integer,intent(in) :: nunocc
    !> Number of atoms in full molecule 
    integer,intent(in) :: nAtoms
    !> Distances between super fragments, NOT standard atoms (see set_super_distance_table) 
    real(realk),intent(in) :: SuperDistanceTable(natoms,natoms)
    !> dofrag(P) is true if P is central atom for a super fragment (see main_fragment_driver)
    logical,intent(in) :: dofrag(natoms)
    !> Single super fragments
    type(ccatom),intent(inout) :: Fragments(natoms)
    !> Pair cutoff distance used for current job list
    real(realk),intent(in) :: oldpaircut
    !> Pair cutoff distance to be used for new job list
    real(realk),intent(in) :: newpaircut
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Occupied orbitals
    type(ccorbital), intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals
    type(ccorbital), intent(in) :: UnoccOrbitals(nunocc)
    !> Current job list which will be appended with new pairs
    type(joblist),intent(inout) :: jobs
    type(joblist) :: oldjobs
    integer :: i,j,nsingle,npairold,npairnew,npairdelta, nold,nnew,k,nbasisFragment
    integer,pointer :: atom1(:), atom2(:), jobsize(:),order(:)
    logical,pointer :: occAOS(:,:), unoccAOS(:,:), REDoccAOS(:,:), REDunoccAOS(:,:)
    logical,pointer :: occpairAOS(:), unoccpairAOS(:)
    real(realk) :: dist

    ! Counting: Number of single fragments, old pairs, and new pairs
    ! **************************************************************
    nsingle=0
    npairold=0
    npairnew=0
    do i=1,natoms
       if(.not. dofrag(i)) cycle ! atom "i" is not central atom in a super fragment
       nsingle = nsingle +1 ! number of single super fragments

       do j=i+1,natoms
          if(.not. dofrag(j)) cycle ! atom "j" is not central atom in a super fragment

          ! Number of pairs for old cutoff
          if(SuperDistanceTable(i,j) < oldpaircut) npairold = npairold+1

          ! Number of pairs for new cutoff
          if(SuperDistanceTable(i,j) < newpaircut) npairnew = npairnew+1

       end do
    end do


    ! Sanity check 1: Number current jobs should be sum of single jobs + old pairs
    if( jobs%njobs /= (nsingle+npairold) ) then
       write(DECinfo%output,*) 'Number of single super frags: ', nsingle
       write(DECinfo%output,*) 'Number of pair   super frags: ', npairold
       write(DECinfo%output,*) 'Number of jobs in job list  : ', jobs%njobs
       call lsquit('expand_joblist_to_include_more_pairs: Error for job bookkeeping!',-1)
    end if

    ! Sanity check 2: New pairs to include
    if(npairold >= npairnew) then
       write(DECinfo%output,*) 'Old: Paircut/npairs: ', oldpaircut,npairold
       write(DECinfo%output,*) 'New: Paircut/npairs: ', newpaircut,npairnew
       call lsquit('expand_joblist_to_include_more_pairs: No new pairs to include!',-1)
    end if


    ! Save info for current job list
    ! ******************************

    ! Copy old job list
    call copy_joblist(jobs,oldjobs)
    nold = jobs%njobs


    ! Delete current job list and reallocate with new dimensions
    ! **********************************************************
    call free_joblist(jobs)
    nnew = nsingle + npairnew  ! New jobs: Same single fragments as before + a larger number of pairs
    call init_joblist(nnew,jobs)    
    npairdelta = nnew - nold  ! Number of pairs not included in old list


    ! Set the first "nold" job in new job list to be identical to old job list
    ! ************************************************************************
    jobs%atom1(1:nold) = oldjobs%atom1(1:nold)
    jobs%atom2(1:nold) = oldjobs%atom2(1:nold)
    jobs%jobsize(1:nold) = oldjobs%jobsize(1:nold)
    jobs%jobsdone(1:nold) = oldjobs%jobsdone(1:nold)
    jobs%nslaves(1:nold) = oldjobs%nslaves(1:nold)
    jobs%nocc(1:nold) = oldjobs%nocc(1:nold)
    jobs%nvirt(1:nold) = oldjobs%nvirt(1:nold)
    jobs%nbasis(1:nold) = oldjobs%nbasis(1:nold)
    jobs%ntasks(1:nold) = oldjobs%ntasks(1:nold)
    jobs%flops(1:nold) = oldjobs%flops(1:nold)
    jobs%LMtime(1:nold) = oldjobs%LMtime(1:nold)
    jobs%load(1:nold) = oldjobs%load(1:nold)


    ! Set logical vectors with fragment orbital information
    ! *****************************************************
    ! occAOS(i,j) is true (false) if occupied orbital "i" is (not) included in super fragment "j"
    call mem_alloc(occAOS,nocc,natoms)
    call mem_alloc(unoccAOS,nunocc,natoms)
    ! same for reduced fragment spaces (see ccatom type)
    call mem_alloc(REDoccAOS,nocc,natoms)
    call mem_alloc(REDunoccAOS,nunocc,natoms)
    call get_logical_vectors_for_AOS(nocc,nunocc,natoms,dofrag,Fragments,&
         & occAOS, unoccAOS, REDoccAOS, REDunoccAOS)




    ! Set new job list
    ! ****************

    k = 0
    call mem_alloc(occpairAOS,nocc)
    call mem_alloc(unoccpairAOS,nunocc)
    call mem_alloc(atom1,npairdelta)
    call mem_alloc(atom2,npairdelta)
    call mem_alloc(jobsize,npairdelta)

    do i=1,natoms  ! Loop over atoms
       if(.not. dofrag(i)) cycle  ! No superfragment for atom i

       do j=i+1,natoms
          if(.not. dofrag(j) ) cycle ! No super fragment for atom j

          ! Distance between atoms i and j
          dist = SuperDistanceTable(i,j)

          ! Add pairs with distances between old and new paircut
          AddPairToJobList: if( (dist < newpaircut) .and. (dist >= oldpaircut) ) then  

             if(dist < DECinfo%PairReductionDistance) then  
                ! Merge AOS for fragment 1 and 2 for standard pair
                call get_logical_pair_vector(nocc,occAOS(1:nocc,i),occAOS(1:nocc,j),occpairAOS)
                call get_logical_pair_vector(nunocc,unoccAOS(1:nunocc,i),&
                     &unoccAOS(1:nunocc,j),unoccpairAOS)

             else
                ! Merge AOS for fragment 1 and 2 for reduced pair
                call get_logical_pair_vector(nocc,REDoccAOS(1:nocc,i),REDoccAOS(1:nocc,j),occpairAOS)
                call get_logical_pair_vector(nunocc,REDunoccAOS(1:nunocc,i),&
                     & REDunoccAOS(1:nunocc,j),unoccpairAOS)

             end if


             ! Set pair job info
             ! *****************
             k=k+1 
             atom1(k) = i
             atom2(k) = j
             ! Job size is defined as occupied AOS * unoccupied AOS dimensions * nbasis
             call get_nbasis_for_fragment(nocc,nunocc,occpairAOS,unoccpairAOS,&
                  & OccOrbitals,UnoccOrbitals,MyMolecule,nbasisFragment)
             jobsize(k) = count(occpairAOS(1:nocc))*count(unoccpairAOS(1:nunocc))*nbasisFragment

             if(jobsize(k)<1) then
                print *, 'dist', dist
                print *, 'k=',k
                print *, 'pair=', i,j
                print *, 'jobsize: ', jobsize(k)
                call lsquit('expand_joblist_to_include_more_pairs: Non-positive job size!',DECinfo%output)
             end if

          end if AddPairToJobList

       end do
    end do


    ! Sanity check
    if(k/=npairdelta) then
       call lsquit('expand_joblist_to_include_more_pairs: &
            & Counter does not match number of new jobs',-1)
    end if


    ! Sort newly added pair jobs according to size (largest elements first)
    call mem_alloc(order,npairdelta)
    call integer_inv_sort_with_tracking(jobsize,order,npairdelta)

    ! Put job info into job structure
    ! (note that jobsize has already been sorted, while atom1 and atom2 need to be sorted)
    do i=1,npairdelta
       jobs%jobsize(nold+i) = jobsize(i)
       jobs%atom1(nold+i) = atom1(order(i))
       jobs%atom2(nold+i) = atom2(order(i))
    end do


    call free_joblist(oldjobs)
    call mem_dealloc(jobsize)
    call mem_dealloc(order)
    call mem_dealloc(occpairAOS)
    call mem_dealloc(unoccpairAOS)
    call mem_dealloc(atom1)
    call mem_dealloc(atom2)
    call mem_dealloc(occAOS)
    call mem_dealloc(unoccAOS)
    call mem_dealloc(REDoccAOS)
    call mem_dealloc(REDunoccAOS)

  end subroutine expand_joblist_to_include_more_pairs


  !> Determine logical vector for occ and unocc AOS where e.g.
  !> occAOS(i,P) is true (false) if occupied orbital "i" is (not) included in super fragment "P".
  !> The same is done for fragment spaces of reduced size.
  !> \author Kasper Kristensen
  !> \date October 2012
  subroutine get_logical_vectors_for_AOS(nocc,nunocc,natoms,dofrag,Fragments,&
       & occAOS, unoccAOS, REDoccAOS, REDunoccAOS)

    implicit none
    !> Number of occupied orbitals in full molecule
    integer,intent(in) :: nocc
    !> Number of unoccupied orbitals in full molecule
    integer,intent(in) :: nunocc
    !> Number of atoms in full molecule 
    integer,intent(in) :: nAtoms
    !> dofrag(P) is true if P is central atom for a super fragment (see main_fragment_driver)
    logical,intent(in) :: dofrag(natoms)
    !> Single super fragments
    type(ccatom),intent(inout) :: Fragments(natoms)
    !> Logical vector for occupied AOS
    logical,intent(inout) :: occAOS(nocc,natoms)
    !> Logical vector for unoccupied AOS
    logical,intent(inout) :: unoccAOS(nunocc,natoms)
    !> Logical vector for reduced occupied AOS
    logical,intent(inout) :: REDoccAOS(nocc,natoms)
    !> Logical vector for reduced unoccupied AOS
    logical,intent(inout) :: REDunoccAOS(nunocc,natoms)
    integer :: i,idx,P

    ! Init
    occAOS = .false.
    unoccAOS = .false.
    REDoccAOS = .false.
    REDunoccAOS = .false.

    
    ! Set logical vectors
    ! *******************

    do P=1,natoms  ! loop over atoms
       if(.not. dofrag(P)) cycle  ! skip if P is not a super fragment

       ! Set occupied AOS for P
       do i=1,Fragments(P)%noccAOS
          ! index "idx" in total occupied orbital list for orbital "i" in fragment orbital list
          idx = Fragments(P)%occAOSidx(i)    
          occAOS(idx,P) = .true.
       end do

       ! Set reduced occupied AOS for P
       do i=1,Fragments(P)%REDnoccAOS
          idx = Fragments(P)%REDoccAOSidx(i)  
          REDoccAOS(idx,P) = .true.
       end do

       ! Set unoccupied AOS for P
       do i=1,Fragments(P)%nunoccAOS
          idx = Fragments(P)%unoccAOSidx(i)  
          unoccAOS(idx,P) = .true.
       end do

       ! Set reduced unoccupied AOS for P
       do i=1,Fragments(P)%REDnunoccAOS
          idx = Fragments(P)%REDunoccAOSidx(i)  
          REDunoccAOS(idx,P) = .true.
       end do

    end do
    

  end subroutine get_logical_vectors_for_AOS


  !> Copy super fragment job list (the new copied joblist is also initialized here)
  !> \author Kasper Kristensen
  !> \date November 2012
  subroutine copy_joblist(jobs,jobscopy)
    implicit none
    !> Original job list
    type(joblist),intent(in) :: jobs
    !> Copy of job list
    type(joblist),intent(inout) :: jobscopy
    integer :: i

    ! Init new job list
    call init_joblist(jobs%njobs,jobscopy)


    ! Copy old job list information
    jobscopy%atom1 = jobs%atom1
    jobscopy%atom2 = jobs%atom2
    jobscopy%jobsize = jobs%jobsize
    jobscopy%jobsdone = jobs%jobsdone
    jobscopy%nslaves = jobs%nslaves
    jobscopy%nocc = jobs%nocc
    jobscopy%nvirt = jobs%nvirt
    jobscopy%nbasis = jobs%nbasis
    jobscopy%ntasks = jobs%ntasks
    jobscopy%flops = jobs%flops
    jobscopy%LMtime = jobs%LMtime
    jobscopy%load = jobs%load

  end subroutine copy_joblist

  !> \brief Compare list of single super fragments 
  !> for two different super fragment job lists.
  !> If they do not match, the program prints error message and quits.
  !> (deliberately we do not compare pairs since they may be different if
  !>  the pair cutoff distance is increased).
  subroutine superfragment_sanity_check(jobs1,jobs2)
    implicit none
    !> The first job list
    type(joblist),intent(in) :: jobs1
    !> The second job list
    type(joblist),intent(in) :: jobs2
    integer :: i,n1,n2,max1,max2
    logical,pointer :: singlejobs1(:),singlejobs2(:)


    ! Check 1: The same number of single jobs in both job lists
    ! *********************************************************
    
    ! Job list 1
    n1=0
    max1=0
    do i=1,jobs1%njobs
       ! Single fragment if atom1=atom2
       if(jobs1%atom1(i)==jobs1%atom2(i)) then
          n1 = n1+1
          if( jobs1%atom1(i) > max1) max1=jobs1%atom1(i)   ! Max atomic index
       end if
    end do

    ! Job list 2
    n2=0
    max2=0
    do i=1,jobs2%njobs
       if(jobs2%atom1(i)==jobs2%atom2(i)) then
          n2 = n2+1
          if( jobs2%atom1(i) > max2) max2=jobs2%atom1(i)   ! Max atomic index
       end if
    end do

    ! Same number of single fragments
    if(n1/=n2) then
       write(DECinfo%output,*) 'n1,n2: ',n1,n2
       call lsquit('superfragment_sanity_check: The number of single fragments differ!',-1)
    end if

    ! Same max index
    if(max1/=max2) then
       write(DECinfo%output,*) 'max1,max2: ',max1,max2
       call lsquit('superfragment_sanity_check: Maximum atomic indices differ!',-1)
    end if


    ! Sanity check 2: The single fragment indices are identical
    ! *********************************************************
    call mem_alloc(singlejobs1,max1)
    call mem_alloc(singlejobs2,max1)  ! we know that max1=max2
    singlejobs1=.false.
    singlejobs2=.false.
    
    ! Set logical atomic index vector for jobs1
    do i=1,jobs1%njobs
       if(jobs1%atom1(i)==jobs1%atom2(i)) then
          ! Single fragment "jobs1%atom1(i)" is present in job list 1
          singlejobs1(jobs1%atom1(i)) = .true.
       end if
    end do

    ! Set logical atomic index vector for jobs2
    do i=1,jobs2%njobs
       if(jobs2%atom1(i)==jobs2%atom2(i)) then
          ! Single fragment "jobs2%atom1(i)" is present in job list 2
          singlejobs2(jobs2%atom1(i)) = .true.
       end if
    end do

    ! Check that logical atomic index vectors are identical
    do i=1,max1
       if(singlejobs1(i) .neqv. singlejobs2(i)) then
          write(DECinfo%output,*)
          write(DECinfo%output,*) 'Single jobs1: ', singlejobs1
          write(DECinfo%output,*)
          write(DECinfo%output,*) 'Single jobs2: ', singlejobs2
          call lsquit('superfragment_sanity_check: Single super fragment indices differ!',-1)
       end if
    end do


    call mem_dealloc(singlejobs1)
    call mem_dealloc(singlejobs2)

  end subroutine superfragment_sanity_check

  
  !> \brief Get number of basis functions for fragment based on logical
  !> list describing occ and virt AOS orbitals.
  !> \author Kasper Kristensen
  !> \date January 2013
  subroutine get_nbasis_for_fragment(nocc,nunocc,occAOS,unoccAOS,&
       & OccOrbitals,UnoccOrbitals,MyMolecule,nbasisFragment)
    implicit none

    !> Number of occupied orbitals for full molecule
    integer, intent(in) :: nocc
    !> Number of unoccupied orbitals for full molecule
    integer, intent(in) :: nunocc
    !> Logical vector describing occ AOS for fragment
    logical,intent(in) :: occAOS(nocc)
    !> Logical vector describing unocc AOS for fragment
    logical,intent(in) :: unoccAOS(nunocc)
    !> Occupied orbitals
    type(ccorbital), intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals
    type(ccorbital), intent(in) :: UnoccOrbitals(nunocc)
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Number of basis functions in fragment
    integer,intent(inout) :: nbasisFragment
    logical,pointer :: FragmentAtoms(:)
    integer :: i,j


    ! FragmentAtoms: Defines which atoms are included in atomic fragment extent describing
    !                molecular orbitals for fragment.
    call mem_alloc(FragmentAtoms,MyMolecule%natoms)
    FragmentAtoms=.false.


    ! Occupied orbitals
    do i=1,nocc  

       if(occAOS(i)) then ! occupied orbital "i" is in fragment AOS

          do j=1,OccOrbitals(i)%numberofatoms 
             ! Atom "j" in list is included in atomic fragment extent
             FragmentAtoms(OccOrbitals(i)%atoms(j)) = .true.
          end do

       end if
    end do


    ! Unoccupied orbitals
    do i=1,nunocc  

       if(unoccAOS(i)) then ! unoccupied orbital "i" is in fragment AOS

          do j=1,UnoccOrbitals(i)%numberofatoms 
             ! Atom "j" in list is included in atomic fragment extent
             FragmentAtoms(UnoccOrbitals(i)%atoms(j)) = .true.
          end do

       end if
    end do


    ! Count number of basis function based on which atoms are included
    nbasisFragment=0
    do i=1,MyMolecule%natoms
       if(FragmentAtoms(i)) then
          nbasisFragment = nbasisFragment + MyMolecule%atom_size(i)
       end if
    end do

    call mem_dealloc(FragmentAtoms)


  end subroutine get_nbasis_for_fragment


  !> \brief Set pair fragment-adapted MO coefficients, simply by
  !> taking the union of fragment-adapted orbitals
  !> for the incoming atomic fragments.
  !> Note 1: These orbitals are NOT orthogonal, and they
  !> are also redundant! Redundancies need to be removed
  !> and orbitals must be orthogonalized before they can be used
  !> in DEC fragment calculation (see pair_fragment_adapted_transformation_matrices)
  !> Note 2: The EOS orbitals (assigned to atoms P and Q) are taken out of the orbital
  !> pool here since they should not be modified.
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine set_pair_fragment_adapted_redundant_orbitals(MyMolecule,FragmentP,FragmentQ,&
       & FragmentPQ,noccPQ,nunoccPQ,CoccPQ,CunoccPQ)
    implicit none

    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule    
    !> Fragment P
    type(ccatom),intent(inout) :: fragmentP
    !> Fragment Q
    type(ccatom),intent(inout) :: fragmentQ
    !> Pair Fragment PQ
    type(ccatom),intent(inout) :: FragmentPQ
    !> Number of occupied redundant pair orbitals
    !> (should be sum of fragmentP%noccFA and fragmentQ%noccFA)
    integer,intent(in) :: noccPQ
    !> Number of unoccupied redundant pair orbitals
    !> (should be sum of fragmentP%nunoccFA and fragmentQ%nunoccFA)
    integer,intent(in) :: nunoccPQ
    !> Occupied MO coefficients
    real(realk),intent(inout) :: CoccPQ(FragmentPQ%number_basis,noccPQ)
    !> Unoccupied MO coefficients
    real(realk),intent(inout) :: CunoccPQ(FragmentPQ%number_basis,nunoccPQ)
    integer :: i,j,ix,jx,offset,orbidx,nP,nQ
    integer,pointer :: Pbasis(:), Qbasis(:), PQbasis(:)
    integer,pointer :: fragP_in_pair(:),fragQ_in_pair(:)

    ! Sanity check
    if( (noccPQ /= fragmentP%noccFA + fragmentQ%noccFA - fragmentP%noccEOS - fragmentQ%noccEOS) .or. &
         & (nunoccPQ /= fragmentP%nunoccFA + fragmentQ%nunoccFA &
         & - fragmentP%nunoccEOS - fragmentQ%nunoccEOS) ) then
       print *, 'noccP: ', fragmentP%noccFA - fragmentP%noccEOS
       print *, 'noccQ: ', fragmentQ%noccFA - fragmentQ%noccEOS
       print *, 'noccPQ: ', noccPQ
       print *, 'nunoccP: ', fragmentP%nunoccFA - fragmentP%nunoccEOS
       print *, 'nunoccQ: ', fragmentQ%nunoccFA - fragmentQ%nunoccEOS
       print *, 'nunoccPQ: ', nunoccPQ
       call lsquit('set_pair_fragment_adapted_redundant_orbitals: Dimension mismatch',-1)
    end if

    ! Set list of orbitals in full molecular list
    call mem_alloc(Pbasis,FragmentP%number_basis)
    call mem_alloc(Qbasis,FragmentQ%number_basis)
    call mem_alloc(PQbasis,FragmentPQ%number_basis)
    call which_orbitals_in_atomic_extents(MyMolecule,FragmentP,Pbasis)
    call which_orbitals_in_atomic_extents(MyMolecule,FragmentQ,Qbasis)
    call which_orbitals_in_atomic_extents(MyMolecule,FragmentPQ,PQbasis)

    ! Set list of atomic orbital extent in pair orbital extent list
    call mem_alloc(fragP_in_pair,FragmentP%number_basis)
    call mem_alloc(fragQ_in_pair,FragmentQ%number_basis)
    call atomicfragAE_in_pairfragAE(FragmentP%number_basis,FragmentPQ%number_basis,&
         & Pbasis,PQbasis,fragP_in_pair)
    call atomicfragAE_in_pairfragAE(FragmentQ%number_basis,FragmentPQ%number_basis,&
         & Qbasis,PQbasis,fragQ_in_pair)


    ! Set occupied FA-orbitals for pair
    ! *********************************

    ! MO coefficients for pair PQ
    !
    ! C_PQ  =  ( C_P C_Q )     (*)
    ! 
    ! i.e. with fragment P orbitals before fragment Q orbitals.

    ! Note: Atomic extent for pair is in union of atomic extent for P and Q.
    !       Therefore atomic extent for pair is (in general) larger than atomic extent 
    !       for either P or Q.
    !       In C_PQ matrix, the orbitals for fragment P are simply copied into PQ fragment coefficient
    !       matrix with the modification that MO coefficients are set to zero for
    !       atomic orbitals which are included in the atomic fragment extent
    !       for Q but not for P (and vice versa).

    ! Put FA orbitals for P into pair fragment orbital matrix
    ! Note: We start counting at the first orbital outside EOS, and use that for FA orbitals
    !       that EOS orbitals are listed before the remaining FA orbitals.
    !       (See fragment_adapted_transformation_matrices).
    CoccPQ=0.0_realk
    offset=fragmentP%noccEOS
    nP = fragmentP%noccFA - fragmentP%noccEOS
    do j=1,nP
       orbidx = j
       do i=1,FragmentP%number_basis
          ix = fragP_in_pair(i)
          CoccPQ(ix,orbidx) = FragmentP%CoccFA(i,j+offset)
       end do
    end do

    ! Put FA orbitals for Q into pair fragment orbital matrix
    offset=fragmentQ%noccEOS
    nQ = fragmentQ%noccFA - fragmentQ%noccEOS
    do j=1,nQ
       orbidx = j + nP
       do i=1,FragmentQ%number_basis
          ix = fragQ_in_pair(i)
          CoccPQ(ix,orbidx) = FragmentQ%CoccFA(i,j+offset)
       end do
    end do


    ! Set unoccupied FA-orbitals for pair
    ! ***********************************
    ! Same strategy as for occupied space.


    CunoccPQ=0.0_realk
    offset=fragmentP%nunoccEOS
    nP = fragmentP%nunoccFA - fragmentP%nunoccEOS
    do j=1,nP
       orbidx = j
       do i=1,FragmentP%number_basis
          ix = fragP_in_pair(i)
          CunoccPQ(ix,orbidx) = FragmentP%CunoccFA(i,j+offset)
       end do
    end do

    offset=fragmentQ%nunoccEOS
    nQ = fragmentQ%nunoccFA - fragmentQ%nunoccEOS
    do j=1,nQ
       orbidx = j + nP
       do i=1,FragmentQ%number_basis
          ix = fragQ_in_pair(i)
          CunoccPQ(ix,orbidx) = FragmentQ%CunoccFA(i,j+offset)
       end do
    end do


    call mem_dealloc(fragP_in_pair)
    call mem_dealloc(fragQ_in_pair)
    call mem_dealloc(Pbasis)
    call mem_dealloc(Qbasis)
    call mem_dealloc(PQbasis)


  end subroutine set_pair_fragment_adapted_redundant_orbitals


  !> \brief Find out which atomic basis functions are included in
  !> the atomic fragment extent.
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine which_orbitals_in_atomic_extents(MyMolecule,Fragment,which_orbitals)

    implicit none

    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule    
    !> Fragment
    type(ccatom),intent(inout) :: fragment
    !> Which atomic orbitals in full AO list are included in atomic fragment extent?
    integer,intent(inout) :: which_orbitals(fragment%number_basis)
    integer :: i,j,idx,atom

    idx=0
    do i=1,fragment%number_atoms  ! loop over atoms in fragment
       atom = fragment%atoms_idx(i)  ! atom index in full atom list

       do j=1,MyMolecule%atom_size(atom)  ! loop over basis functions for atom
          idx =idx+1
          which_orbitals(idx) = MyMolecule%atom_start(atom) + j -1
       end do

    end do

    ! Sanity check
    if(idx/=fragment%number_basis) then
       call lsquit('which_orbitals_in_atomic_extents: &
            & counter does not match number of basis functions in atomic &
            & fragment extent!',-1)
    end if
    

  end subroutine which_orbitals_in_atomic_extents

end module atomic_fragment_operations
