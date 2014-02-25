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
#ifdef VAR_MPI
  use infpar_module
#endif

  ! DEC DEPENDENCIES (within deccc directory)
  ! *****************************************
  use dec_fragment_utils
  use array2_simple_operations!, only: array2_init,array2_matmul,array2_free,array2_add, operator(*)
  use array4_simple_operations!, only: array4_init,array4_extract_eos_indices_both_schemes,&
  use orbital_operations!,only: get_number_of_orbitals_per_atom,orbital_init,&
  !       & copy_orbital
  use mp2_module!,only: max_batch_dimension,get_vovo_integrals,&
  !       & mp2_integrals_memory_for_updating_arrays,&
  !       & mp2_integrals_and_amplitudes
  !       & array4_free, operator(*)

  ! F12 DEPENDENCIES 
  ! *****************************************
  use CABS_operations
  use f12_routines_module

contains



  !> \brief Initialize atomic fragment based on a list of specific AOS orbitals.
  !> For a pair fragment calculation it is assumed that EOSatoms and nEOSatoms
  !> are set and that fragment%pairfrag=.true. before calling this routine
  !> \author Kasper Kristensen
  subroutine atomic_fragment_init_orbital_specific(MyAtom,nunocc, nocc, unocc_list, &
       & occ_list,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,fragment,DoBasis,pairfrag)

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
    type(decorbital), dimension(nOcc), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(decorbital), dimension(nUnocc), intent(in) :: UnoccOrbitals
    !> Fragment to construct
    type(decfrag), intent(inout) :: fragment
    !> Make fragment basis (MO coeff and Fock matrix for fragment)?
    logical, intent(in) :: DoBasis
    !> Is it a pair fragment?
    logical,intent(in) :: pairfrag
    integer :: j,idx,i,natoms,startidx
    integer :: CentralAtom, cntr
    real(realk) :: tcpu, twall
    logical,pointer :: occ_listEFF(:),occEOS(:),unoccEOS(:)
    !> list of atoms with AOS orbitals assigned
    logical,pointer :: all_atoms(:)

    all_atoms  => null()

    ! Integer pointers for some dimensions
    call fragment_init_dimension_pointers(fragment)

    ! Use of fragment-adapted orbitals set to false here, they can be set later using
    ! fragment_adapted_transformation_matrices.
    fragment%fragmentadapted=.false.

    natoms = MyMolecule%natoms

    call LSTIMER('START',tcpu,twall,DECinfo%output)

    call mem_alloc( all_atoms,  natoms )
    all_atoms  = .false.



    ! To be on the safe side, initialize everything to zero or null
    call atomic_fragment_nullify(fragment)


    ! ATOMIC VS. PAIR FRAGMENT
    ! ************************
    IsThisAPairFragment: if(pairfrag) then
       if(DECinfo%PL>0) write(DECinfo%output,*) '-- Initializing pair fragment for atoms : ', fragment%EOSatoms

       ! Sanity checks
       if( (.not. associated(fragment%EOSatoms)) ) then
          call lsquit('atomic_fragment_init_orbital_specific: Pair fragment calculation, &
               & but EOSatoms has not been initiated!',-1)
       end if
       if( fragment%nEOSatoms /= size(fragment%EOSatoms) ) then
          call lsquit('atomic_fragment_init_orbital_specific: Pair fragment calculation, &
               & Numbers of atoms in atomlist is incorrect!',-1)
       end if
       fragment%pairfrag=.true.

    else

       if(DECinfo%PL>0) write(DECinfo%output,*) '-- Initializing standard fragment on atom : ',MyAtom
       fragment%nEOSatoms=1 ! atomic fragment, 1 EOS atom

       if(associated(fragment%EOSatoms)) then  ! dealloc and realloc with the correct size of 1
          call mem_dealloc(fragment%EOSatoms)
          fragment%EOSatoms => null()
       end if
       call mem_alloc(fragment%EOSatoms,fragment%nEOSatoms)
       fragment%EOSatoms(1) = MyAtom
       fragment%pairfrag=.false.
    end if IsThisAPairFragment

    ! Set model to use for fragment calculation (see define_pair_calculations for details)
    if(fragment%pairfrag) then
       fragment%ccmodel = MyMolecule%ccmodel(fragment%EOSatoms(1),fragment%EOSatoms(2))
    else
       fragment%ccmodel = MyMolecule%ccmodel(MyAtom,MyAtom)
    end if


    ! Size of occupied EOS
    ! ********************

    ! Occupied EOS orbitals are those assigned to central atom in fragment
    ! (or central atomS in case of a pair fragment)

    ! Avoid including core orbitals in EOS if frozen core approx. is used
    if(DECinfo%frozencore) then
       startidx=MyMolecule%ncore+1 ! loop from first valence orbital
    else
       startidx=1 ! loop over all occ orbitals
    end if

    ! Logical vector keeping track of EOS orbitals
    call mem_alloc(occEOS,nocc)
    occEOS=.false.
    fragment%noccEOS=0
    ! loop over occupied orbitals (only valence for frozen core)
    OccEOSSize: do j=startidx,nOcc
       CentralAtom=OccOrbitals(j)%centralatom

       ! Loop over atoms in pair fragment list (just one atom, MyAtom, if it is not a pairfragment)
       do i=1,fragment%nEOSatoms
          if( CentralAtom == fragment%EOSatoms(i) ) then ! Orbital is included in the EOS
             fragment%noccEOS         = fragment%noccEOS + 1
             occEOS(j)                = .true.
             all_atoms( CentralAtom ) = .true.
          end if
       end do


    end do OccEOSSize

    ! Size of unoccupied EOS
    ! **********************

    ! Logical vector keeping track of EOS orbitals
    call mem_alloc(unoccEOS,nunocc)
    unoccEOS=.false.
    ! Virtual EOS orbitals are those assigned to the central atom
    fragment%nunoccEOS=0
    UnoccEOSLoop: do j=1,nunocc
       CentralAtom=UnoccOrbitals(j)%centralatom

       ! Loop over atoms in pair fragment list (just one atom, Myatom, if it is not a pairfragment)
       do i=1,fragment%nEOSatoms
          if( CentralAtom==fragment%EOSatoms(i) ) then ! Orbital is included in the EOS
             fragment%nunoccEOS       = fragment%nunoccEOS + 1
             unoccEOS(j)              = .true.
             all_atoms( CentralAtom ) = .true.

             ! Special case: Only occupied partitioning
             ! ----------------------------------------
             ! When we are only interested in the occupied partitioning scheme,
             ! there are effectively zero virtual EOS orbitals.
             ! However, setting fragment%nunoccEOS to 0 would cause numerous
             ! problems many places in the code, since some arrays would have zero size.
             ! For now, we solve this problem in a pragmatic and dirty manner:
             ! We simply initialize an unocc EOS containing a single dummy orbital.
             ! At some point we want to separate out the occ and virt
             ! partitioning scheme, and when this is done, this temporary solution
             ! will be superfluous.
             if(DECinfo%onlyoccpart) then
                exit UnoccEOSLoop
             end if

          end if
       end do


    end do UnoccEOSLoop


    ! Size of occupied AOS - number of "true" elements in logical occupied vector
    ! ***************************************************************************
    call mem_alloc(occ_listEFF,nocc)
    occ_listEFF = occ_list
    if(DECinfo%frozencore) then  ! for frozen core, core orbitals are removed from AOS orbital list
       occ_listEFF(1:MyMolecule%ncore) = .false.
    end if

    fragment%noccLOC = count(occ_listEFF)
    fragment%noccAOS => fragment%noccLOC ! Local orbitals = AOS orbitals

    ! Size of unoccupied AOS - number of "true" elements in logical unoccupied vector
    ! *******************************************************************************
    fragment%nunoccLOC = count(unocc_list)
    fragment%nunoccAOS => fragment%nunoccLOC  ! AOS orbitals = local orbitals

    ! Fragment-adapted information, for now set equal to local dimensions
    fragment%noccFA = fragment%noccLOC
    fragment%nunoccFA = fragment%nunoccLOC

    ! Occupied orbital indices
    ! ************************
    call mem_alloc(fragment%occEOSidx,fragment%noccEOS)
    call mem_alloc(fragment%occAOSidx,fragment%noccAOS)
    fragment%occEOSidx=0
    fragment%occAOSidx=0

    ! Note that the procedure below ensures that the EOS orbitals are always
    ! listed before the remaining AOS orbitals

    ! Loop over EOS orbitals
    idx=0
    do j=startidx,nOcc   ! no core orbitals in EOS for frozen core approx

       if(occEOS(j)) then

          idx=idx+1
          fragment%occEOSidx(idx)=j
          fragment%occAOSidx(idx)=j   ! also set
          all_atoms( OccOrbitals(j)%centralatom ) = .true.

       end if


    end do

    if(idx /= fragment%noccEOS) then
       call lsquit('atomic_fragment_init_orbital_specific: idx /= fragment%noccEOS',-1)
    end if

    ! Loop over remaining AOS orbitals (not EOS because they have already been set)
    do i=1,nocc

       if(occ_listEFF(i) .and. (.not. occEOS(i)) ) then
          idx=idx+1
          fragment%occAOSidx(idx) = i
          all_atoms( OccOrbitals(i)%centralatom ) = .true.
       end if

    end do
    if(idx /= fragment%noccAOS) &
         & call lsquit('atomic_fragment_init_orbital_specific: idx /= fragment%noccAOS',-1)



    ! Virtual EOS orbital indices
    ! ***************************
    call mem_alloc(fragment%unoccEOSidx,fragment%nunoccEOS)
    call mem_alloc(fragment%unoccAOSidx,fragment%nunoccAOS)
    fragment%unoccEOSidx = 0
    fragment%unoccAOSidx = 0

    ! Loop over EOS orbitals
    idx=0
    do j=1,nunocc

       if(unoccEOS(j)) then
          idx=idx+1
          fragment%unoccEOSidx(idx)=j
          fragment%unoccAOSidx(idx)=j
          all_atoms( UnoccOrbitals(j)%centralatom ) = .true.
       end if

    end do
    if(idx /= fragment%nunoccEOS) then
       call lsquit('atomic_fragment_init_orbital_specific: idx /= fragment%nunoccEOS',DECinfo%output)
    end if

    ! Loop over remaining AOS orbitals (not EOS)
    do i=1,nunocc

       if(unocc_list(i) .and. (.not. unoccEOS(i)) ) then
          idx=idx+1
          fragment%unoccAOSidx(idx) = i
          all_atoms( UnoccOrbitals(i)%centralatom ) = .true.
       end if

    end do
    if(idx /= fragment%nunoccAOS) &
         & call lsquit('atomic_fragment_init_orbital_specific: idx /= fragment%nunoccAOS',-1)


    !set max distance in AOS space
    fragment%RmaxAOS = 0.0E0_realk
    fragment%RaveAOS = 0.0E0_realk
    fragment%RsdvAOS = 0.0E0_realk
    cntr = 0
    do i=1,natoms
      if( all_atoms(i) )then
        do j=i+1,natoms
          if( all_atoms(j) )then
            fragment%RmaxAOS = max(fragment%RmaxAOS,MyMolecule%DistanceTable(j,i))
            fragment%RaveAOS = fragment%RaveAOS + MyMolecule%DistanceTable(j,i)
            cntr = cntr + 1
          endif
        enddo
      endif
    enddo
    fragment%RaveAOS = fragment%RaveAOS / float( cntr )
    do i=1,natoms
      if( all_atoms(i) )then
        do j=i+1,natoms
          if( all_atoms(j) )then
            fragment%RsdvAOS = fragment%RsdvAOS + (MyMolecule%DistanceTable(j,i) - fragment%RaveAOS )**2
          endif
        enddo
      endif
    enddo
    fragment%RsdvAOS = ( fragment%RsdvAOS / float( cntr - 1 ) )**(0.5E0_realk)

    call mem_dealloc( all_atoms )

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
    fragment%FAset=.false.

    ! Set atomic fragment extent info
    call init_atomic_fragment_extent(OccOrbitals,UnoccOrbitals,MyMolecule,fragment)

    ! Only create fragment basis (the information in the "expensive box" in decfrag type definition)
    ! if DoBasis is true.
    if(DoBasis) then
       call atomic_fragment_init_basis_part(nunocc, nocc, OccOrbitals,UnoccOrbitals,&
            & MyMolecule,mylsitem,fragment)
    else
       fragment%BasisInfoIsSet=.false.
    end if

    ! Information only used for pair, set to zero (set correctly in merged_fragment_init)
    fragment%pairdist=0.0E0_realk


    ! Occ and virt energy contributions
    call mem_alloc(fragment%OccContribs,fragment%noccAOS)
    call mem_alloc(fragment%VirtContribs,fragment%nunoccAOS)
    fragment%OccContribs  = 0.0E0_realk
    fragment%VirtContribs = 0.0E0_realk
    fragment%energies     = 0.0E0_realk
    fragment%EoccFOP      = 0.0_realk
    fragment%EvirtFOP     = 0.0_realk
    fragment%LagFOP       = 0.0_realk


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
    call mem_dealloc(occEOS)
    call mem_dealloc(unoccEOS)

    if(DECinfo%PL>0) then
       write(DECinfo%output,*)
       write(DECinfo%output,*) 'FRAGINIT: Initialized fragment :', Fragment%EOSatoms(1)
       write(DECinfo%output,'(a)')      ' -- Target --'
       write(DECinfo%output,'(a,i6)')   ' FRAGINIT: Occ EOS     : ',Fragment%noccEOS
       write(DECinfo%output,'(a,i6)')   ' FRAGINIT: Unocc EOS   : ',Fragment%nunoccEOS
       write(DECinfo%output,'(a)')      ' -- Target + Buffer --'
       write(DECinfo%output,'(a,i6)')   ' FRAGINIT: Occ AOS     : ',Fragment%noccAOS
       write(DECinfo%output,'(a,i6)')   ' FRAGINIT: Unocc AOS   : ',Fragment%nunoccAOS
       write(DECinfo%output,'(a,i6)')   ' FRAGINIT: Basis       : ',Fragment%nbasis
       write(DECinfo%output,'(a,2f10.3)')   ' FRAGINIT: Dist AOS/AE : ',&
       &Fragment%RmaxAOS*bohr_to_angstrom,Fragment%RmaxAE*bohr_to_angstrom
       write(DECinfo%output,*)
    end if

    call LSTIMER('FRAGMENT INIT',tcpu,twall,DECinfo%output)
     
  end subroutine atomic_fragment_init_orbital_specific

  
  !> \brief Initialize atomic fragment based on a list of specific AOS orbitals for F12-related matrices
  !> \author Yang Min Wang
  subroutine atomic_fragment_init_f12(fragment, MyMolecule)
    type(fullmolecule), intent(in) :: MyMolecule
    type(decfrag), intent(inout) :: fragment

    !> F12 Specific Variables
    integer :: nbasis, noccEOS, nunoccEOS, noccfull, nocvAOS, nvirtAOS, ncabsAO, ncabsMO
    integer :: noccAOS, nunoccAOS
    integer :: ix, iy

    ! ============================================================
    !                        F12-Specific                        !
    ! ============================================================
    !> F12 Specific Variables
    nbasis   = fragment%nbasis
    noccEOS  = fragment%noccEOS
    nunoccEOS = fragment%nunoccEOS

    noccAOS  = fragment%noccAOS
    nunoccAOS = fragment%nunoccAOS  
    nocvAOS  = fragment%noccAOS + fragment%nunoccAOS
    nvirtAOS = fragment%nunoccAOS
    
    if(DECinfo%f12debug) then
       !CABS MO and RI AO (Fragment ncabsMO = Molecule ncabsMO) Needs to be changed
       !ncabsAO = size(MyMolecule%Ccabs,1)
       !ncabsMO = size(MyMolecule%Ccabs,2)
       !print *, "associated(MyMolecule%Ccabs)", associated(MyMolecule%Ccabs)
       !print *, "size(MyMolecule%Ccabs,2)", size(MyMolecule%Ccabs,2)
    endif

    if(DECinfo%f12debug) then
       print *, "associated(fragment%Ccabs)", associated(fragment%Ccabs)
       print *, "size(fragment%Ccabs,2)", size(fragment%Ccabs,2)
    endif

    ncabsAO  = size(fragment%Ccabs,1)
    ncabsMO = size(fragment%Ccabs,2)

    if(DECinfo%F12debug) then
       print *, "---------------------------------------"
       print *, " atomic_fragment_init_f12 dec_atom.F90 "
       print *, "---------------------------------------"
       print *, "nbasis: ", nbasis
       print *, "noccEOS: ", noccEOS
       print *, "nunoccEOS: ", nunoccEOS
       print *, "---------------------------------------"
       print *, "nocvAOS", nocvAOS
       print *, "noccAOS", noccAOS
       print *, "nvirtAOS", nvirtAOS
       print *, "ncabsAO", ncabsAO
       print *, "ncabsMO", ncabsMO
       print *, "---------------------------------------"
    end if 
    
    ! hJir
    call mem_alloc(fragment%hJir, noccEOS, ncabsAO)
    do j=1,ncabsAO
       do i=1, fragment%noccEOS
          ix = fragment%occEOSidx(i)
          fragment%hJir(i,:) = MyMolecule%hJir(ix,:)
       enddo
    enddo
  
    ! Krs
    call mem_alloc(fragment%Krs, ncabsAO, ncabsAO)
    call dcopy(ncabsAO*ncabsAO, MyMolecule%Krs, 1, fragment%Krs, 1)

    ! Frs
    call mem_alloc(fragment%Frs, ncabsAO, ncabsAO)
    call dcopy(ncabsAO*ncabsAO, MyMolecule%Frs, 1, fragment%Frs, 1)

    ! Frm
    call mem_alloc(fragment%Frm, ncabsAO, noccAOS)
    do i=1, fragment%noccAOS
       iy = fragment%occAOSidx(i)
       fragment%Frm(:,i) = MyMolecule%Frm(:,iy)
    enddo

    ! Fcp in the order of the index (occ to virt)
    call mem_alloc(fragment%Fcp, ncabsMO, nocvAOS)
    do i=1, fragment%noccAOS
       iy = fragment%occAOSidx(i)
       fragment%Fcp(:,i) = MyMolecule%Fcp(:,iy)
    enddo
  
    do i=fragment%noccAOS+1, fragment%nunoccAOS+fragment%noccAOS
       iy = fragment%unoccAOSidx(i-fragment%noccAOS)
       fragment%Fcp(:,i) = MyMolecule%Fcp(:,iy+MyMolecule%nocc)
    enddo

    if(DECinfo%F12debug) then
       print *, "---------------------------------------"
       print *, " atomic_fragment_init_f12 dec_atom.F90 "
       print *, "---------------------------------------"
       print *,"norm2D(fragment%hJir)", norm2D(fragment%hJir)
       print *,"norm2D(fragment%Krs)" , norm2D(fragment%Krs) 
       print *,"norm2D(fragment%Frs)" , norm2D(fragment%Frs) 
       print *,"norm2D(fragment%Frm)" , norm2D(fragment%Frm) 
       print *,"norm2D(fragment%Fcp)" , norm2D(fragment%Fcp) 
       print *, "---------------------------------------"
    endif

  end subroutine atomic_fragment_init_f12

  !> \brief Initialize atomic fragments by simply including neighbouring atoms within
  !> a certain distance.
  !> \author Kasper Kristensen
  !> \date November 2013
  subroutine init_estimated_atomic_fragments(nOcc,nUnocc,OccOrbitals,UnoccOrbitals, &
       & MyMolecule,mylsitem,DoBasis,init_radius,dofrag,AtomicFragments)

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
    type(decorbital), dimension(nOcc), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(decorbital), dimension(nUnocc), intent(in) :: UnoccOrbitals
    !> Make fragment basis (MO coeff and Fock matrix for fragment)?
    logical, intent(in) :: DoBasis
    !> Distance beyond which to include neighbour atoms
    real(realk),intent(in) :: init_radius
    !> List of which atoms have orbitals assigned
    logical,intent(in),dimension(MyMolecule%natoms) :: dofrag
    !> Fragments to construct
    type(decfrag), intent(inout), dimension(MyMolecule%natoms) :: AtomicFragments
    integer :: i,MyAtom

    do i=1,MyMolecule%natoms
       if(.not. dofrag(i)) cycle
       MyAtom=i
       call atomic_fragment_init_within_distance(MyAtom,&
            & nOcc,nUnocc,OccOrbitals,UnoccOrbitals, &
            & MyMolecule,mylsitem,DoBasis,init_radius,AtomicFragments(MyAtom))
    end do

  end subroutine init_estimated_atomic_fragments

  !> \brief Nullify pointers
  !> \author Marcin Ziolkowski
  !> \param fragment Molecular fragment
  subroutine atomic_fragment_nullify(fragment)

    implicit none
    type(decfrag), intent(inout) :: fragment

    fragment%occEOSidx => null()
    fragment%unoccEOSidx => null()
    fragment%occAOSidx => null()
    fragment%unoccAOSidx => null()
    fragment%coreidx => null()


    fragment%occAOSorb => null()
    fragment%unoccAOSorb => null()

    fragment%idxo => null()
    fragment%idxu => null()

    fragment%atoms_idx => null()

    fragment%Co => null()
    fragment%Cv => null()

    !Free CABS MO F12
    fragment%Ccabs => null()
    
    !Free RI MO F12
    fragment%Cri => null()

    !Free F12-matrices
    fragment%hJir => null()
    fragment%Krs => null()
    fragment%Frs => null()
    fragment%Fac => null()
     
    fragment%fock => null()
    fragment%ppfock => null()
    fragment%ccfock => null()
    fragment%qqfock => null()
    fragment%ppfockLOC => null()
    fragment%qqfockLOC => null()
    fragment%ppfockFA => null()
    fragment%qqfockFA => null()

    fragment%OccContribs => null()
    fragment%VirtContribs => null()

    fragment%OccMat => null()
    fragment%VirtMat => null()
    fragment%CoFA => null()
    fragment%CvFA => null()
    fragment%CDocceival => null()
    fragment%CDunocceival => null()
    fragment%t1 => null()
    fragment%t1_occidx => null()
    fragment%t1_virtidx => null()


  end subroutine atomic_fragment_nullify



  !> \brief Initialize atomic fragment based on a list of atoms for the occupied and
  !> unoccupied orbital spaces.
  !> For a pair fragment calculation it is assumed that EOSatoms and nEOSatoms
  !> are set and that fragment%pairfrag=.true. before calling this routine
  !> \author Kasper Kristensen
  subroutine atomic_fragment_init_atom_specific(MyAtom,natoms,Unocc_atoms, &
       & Occ_atoms,nOcc,nUnocc,OccOrbitals,UnoccOrbitals, &
       & MyMolecule,mylsitem,fragment,DoBasis,Pairfrag)

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
    type(decorbital), dimension(nOcc), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(decorbital), dimension(nUnocc), intent(in) :: UnoccOrbitals
    !> Fragment to construct
    type(decfrag), intent(inout) :: fragment
    !> Make fragment basis (MO coeff and Fock matrix for fragment)?
    logical, intent(in) :: DoBasis
    !> Is it a pair fragment?
    logical,intent(in) :: pairfrag
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
         & occ_list,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,fragment,DoBasis,pairfrag)
    call mem_dealloc(occ_list)
    call mem_dealloc(unocc_list)


  end subroutine atomic_fragment_init_atom_specific



  !> \brief Initialize atomic fragment by simply including all local orbitals assigned
  !> to atoms which are closer to central atom than the input radius.
  !> \author Kasper Kristensen
  !> \date October 2013
  subroutine atomic_fragment_init_within_distance(MyAtom,&
       & nOcc,nUnocc,OccOrbitals,UnoccOrbitals, &
       & MyMolecule,mylsitem,DoBasis,init_radius,Fragment)

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
    !> Information about DEC occupied orbitals
    type(decorbital), dimension(nOcc), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(decorbital), dimension(nUnocc), intent(in) :: UnoccOrbitals
    !> Make fragment basis (MO coeff and Fock matrix for fragment)?
    logical, intent(in) :: DoBasis
    !> Distance beyond which to include neighbour atoms
    real(realk),intent(in) :: init_radius
    !> Fragment to construct
    type(decfrag), intent(inout) :: fragment
    integer :: natoms
    logical, pointer :: Unocc_atoms(:), Occ_atoms(:)
    integer, pointer :: nocc_per_atom(:),nunocc_per_atom(:)
    logical :: pairfrag

    ! How many occ/unocc orbitals assigned to each atom
    natoms = MyMolecule%natoms
    call mem_alloc( nocc_per_atom,   natoms )
    call mem_alloc( nunocc_per_atom, natoms )

    nocc_per_atom   = get_number_of_orbitals_per_atom(OccOrbitals,nocc,natoms)
    nunocc_per_atom = get_number_of_orbitals_per_atom(UnoccOrbitals,nunocc,natoms)

    ! Determine logical vectors describing which atoms to include in fragment,
    ! i.e., atoms where the distance to MyAtom is smaller than init_radius
    call mem_alloc( occ_atoms,   natoms )
    call mem_alloc( unocc_atoms, natoms )

    call InitialFragment(natoms,nocc_per_atom,nunocc_per_atom,MyMolecule%distancetable(:,MyAtom),&
         & init_radius,Occ_atoms, Unocc_atoms)

    ! Init atomic fragment
    pairfrag=.false.
    call atomic_fragment_init_atom_specific(MyAtom,natoms,Unocc_atoms, &
         & Occ_atoms,nOcc,nUnocc,OccOrbitals,UnoccOrbitals, &
         & MyMolecule,mylsitem,fragment,DoBasis,pairfrag)

    call mem_dealloc(nocc_per_atom)
    call mem_dealloc(nunocc_per_atom)
    call mem_dealloc(occ_atoms)
    call mem_dealloc(unocc_atoms)

  end subroutine atomic_fragment_init_within_distance




  !> Get transformation matrices from local basis to fragment-adapted basis.
  !> Note: This routine also makes the pointers for MO coefficients and MO Fock matrices
  !> point to the newly generated MO coefficients and MO Fock matrices in the FO basis.
  !> \author Kasper Kristensen
  !> \date February 2013
  subroutine fragment_adapted_transformation_matrices(Myfragment)

    implicit none
    !> Atomic fragment where all quantities are expressed in local basis
    type(decfrag), intent(inout) :: Myfragment
    real(realk),pointer :: VirtMat(:,:),OccMat(:,:),tmpU(:,:),occeival(:),virteival(:)
    integer :: i,j,ix,jx
    integer :: noccTRANS,nvirtTRANS,noccEOS,nvirtEOS,offset,nocc,nvirt
    real(realk),pointer :: OU(:,:),VU(:,:),Oeival(:),Veival(:),OUred(:,:),VUred(:,:)
    logical,pointer :: virtEOS(:),occEOS(:), VirtOrbs(:), OccOrbs(:)

    ! Note: When we do this, the fragment has to use MO coeff and dimensions for LOCAL orbitals
    call fragment_basis_point_to_LOs(Myfragment)
    ! Then, after the fragment-adapted orbitals have been generated we direct the
    ! fragment data pointers to the fragment-adapted data rather than local.
    ! (See the end of this subroutine)

    ! Init stuff
    nocc = Myfragment%noccAOS
    nvirt = Myfragment%nunoccAOS
    noccEOS = Myfragment%noccEOS
    nvirtEOS = Myfragment%nunoccEOS
    call mem_alloc(OU,nocc,nocc)
    call mem_alloc(VU,nvirt,nvirt)
    call mem_alloc(Oeival,nocc)
    call mem_alloc(Veival,nvirt)

    ! Number of orbitals to transform (leave EOS untouched in standard DEC calculation)
    noccTRANS = nocc-noccEOS
    nvirtTRANS = nvirt-nvirtEOS

    ! Special debug case: AOS=EOS, transform all orbitals
    if(Myfragment%noccEOS == Myfragment%noccAOS) then
       noccTRANS = nocc
    end if
    if(Myfragment%nunoccEOS == Myfragment%nunoccAOS) then
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
       virtEOS(Myfragment%idxu(i)) = .true.
    end do

    ! Virtual correlation density matrix
    call mem_alloc(VirtMat,nvirtTRANS,nvirtTRANS)
    call mem_alloc(tmpU,nvirtTRANS,nvirtTRANS)
    call mem_alloc(virteival,nvirtTRANS)

    if(Myfragment%nunoccEOS ==Myfragment%nunoccAOS) then
       ! Transform all orbitals
       VirtMat = Myfragment%VirtMat
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
             VirtMat(jx,ix) = Myfragment%VirtMat(j,i)
          end do
       end do

    end if

    ! Diagonalize virtual correlation density matrix to define fragment-adapted virtual AOS orbitals
    call solve_eigenvalue_problem_unitoverlap(nvirtTRANS,VirtMat,virteival,tmpU)

    ! Now the fragment-adapted orbitals (psi) are given from local orbitals (phi) as:
    ! psi(c) = sum_a tmpU(a,c) phi(a)

    ! However, only orbitals OUTSIDE the EOS are described by tmpU while the local
    ! EOS orbitals have not been touched. To understand the ordering in the
    ! output VU matrix consider the example where there are 5 AOS orbitals,
    ! and where 2 of these are EOS orbitals.
    ! In the AOS list of orbitals the EOS orbitals are number 2 and 4, i.e.
    ! Myfragment%idxu = [2,4].
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
       VU(Myfragment%idxu(i),i) = 1.0_realk
       Veival(i) = 1.0_realk
    end do


    ! Set remaining AOS  indices (not EOS)
    ! ------------------------------------
    offset = nvirt-nvirtTRANS    ! offset to put all non-EOS orbitals after the EOS orbitals
    do j=1,nvirtTRANS  ! loop over fragment-adapted orbital indices
       Veival(j+offset) = virteival(j)
       ix=0
       localloop1: do i=1,nvirt  ! loop over local orbital indices
          if(virtEOS(i)) cycle localloop1   ! not consider local EOS, which have already been set
          ix=ix+1
          VU(i,j+offset) = tmpU(ix,j)
       end do localloop1
    end do

    call mem_dealloc(tmpU)




    ! ============================================================
    !                         OCCUPIED SPACE                     !
    ! ============================================================
    ! Same strategy as for virtual space above.

    ! Logical vector for occupied EOS
    call mem_alloc(occEOS,nocc)
    occEOS=.false.
    do i=1,noccEOS
       occEOS(Myfragment%idxo(i)) = .true.
    end do

    ! Occupied correlation density matrix
    call mem_alloc(OccMat,noccTRANS,noccTRANS)
    call mem_alloc(occeival,noccTRANS)
    call mem_alloc(tmpU,noccTRANS,noccTRANS)


    if(Myfragment%noccEOS ==Myfragment%noccAOS) then
       ! Transform all orbitals
       OccMat = Myfragment%OccMat
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
             OccMat(jx,ix) = Myfragment%OccMat(j,i)
          end do
       end do

    end if

    ! Diagonalize occupied correlation density matrix to define fragment-adapted occupied AOS orbitals
    call solve_eigenvalue_problem_unitoverlap(noccTRANS,OccMat,occeival,tmpU)

    ! Set blocks of OU in the same way as for VU above.
    OU = 0.0_realk

    do i=1,noccEOS
       OU(Myfragment%idxo(i),i) = 1.0_realk
       Oeival(i) = 1.0_realk
    end do

    offset = nocc-noccTRANS
    do j=1,noccTRANS  ! loop over fragment-adapted orbital indices
       Oeival(j+offset) = occeival(j)
       ix=0
       localloop2: do i=1,nocc  ! loop over local orbital indices
          if(occEOS(i)) cycle localloop2
          ix=ix+1
          OU(i,j+offset) = tmpU(ix,j)
       end do localloop2
    end do

    call mem_dealloc(tmpU)
    call mem_dealloc(VirtMat)
    call mem_dealloc(OccMat)
    call mem_dealloc(occEOS)
    call mem_dealloc(virtEOS)



    !=====================================================================================
    ! At this point OU and VU contains the transformation matrices between local and
    ! fragment-adapted bases (first index: local, second index: fragment-adapted).
    ! We now want to remove the fragment-adapted orbitals with small eigenvalues
    ! as defined by Myfragment%RejectThr.
    !=====================================================================================


    ! Find which fragment-adapted orbitals to include
    call mem_alloc(VirtOrbs,nvirt)
    call mem_alloc(OccOrbs,nocc)
    call which_fragment_adapted_orbitals(Myfragment,nocc,nvirt,Oeival,Veival,OccOrbs,VirtOrbs)
    write(DECinfo%output,'(1X,a,i6,a,i6,a)') 'FOP Removed ', nvirt-count(VirtOrbs), ' of ', &
         & nvirtTRANS, ' fragment-adapted virtual orbitals'
    write(DECinfo%output,'(1X,a,i6,a,i6,a)') 'FOP Removed ', nocc-count(OccOrbs), ' of ', &
         & noccTRANS, ' fragment-adapted occupied orbitals'


    ! Set fragment-adapted MO coefficients
    ! ------------------------------------

    ! The fragment-adapted occ orbitals (psi) are given from local orbitals (phi) as:
    !
    ! psi(i) = sum_k OU(k,i) phi(k)         (only for certain "i" defined by OccOrbs)
    !
    ! We want to include only the orbitals "i" defined by OccOrbs to generate a
    ! matrix OUred where:
    ! Dimension 1: Number of occupied AOS orbitals in LOCAL basis
    ! Dimension 2: Number of occupied AOS orbitals in FRAGMENT-ADAPTED basis
    ! In general dimension 1 is larger than dimension 2.
    Myfragment%noccFA = count(OccOrbs)
    call mem_alloc(OUred,Myfragment%noccAOS,Myfragment%noccFA)
    if(associated(Myfragment%CDocceival)) then
       call mem_dealloc(Myfragment%CDocceival)
    end if
    call mem_alloc(Myfragment%CDocceival,Myfragment%noccFA)
    ix=0
    do i=1,nocc
       if(OccOrbs(i)) then
          ix=ix+1
          OUred(:,ix) = OU(:,i)
          Myfragment%CDocceival(ix) = Oeival(i)
       end if
    end do


    ! Same for virtual space
    Myfragment%nunoccFA = count(VirtOrbs)
    call mem_alloc(VUred,Myfragment%nunoccAOS,Myfragment%nunoccFA)
    if(associated(Myfragment%CDunocceival)) then
       call mem_dealloc(Myfragment%CDunocceival)
    end if
    call mem_alloc(Myfragment%CDunocceival,Myfragment%nunoccFA)
    ix=0
    do i=1,nvirt
       if(VirtOrbs(i)) then
          ix=ix+1
          VUred(:,ix) = VU(:,i)
          Myfragment%CDunocceival(ix) = Veival(i)
       end if
    end do

    call mem_dealloc(occeival)
    call mem_dealloc(virteival)


    ! Set fragment-adapted (FA) orbital coefficients: AO-->FA basis
    ! -------------------------------------------------------------

    ! The FA MO coefficient matrix for the occupied space
    ! can be found as (dimensions given below):
    !
    ! (MO fragment-adapted basis)    =     (MO local basis)      *      OUred
    !     (nbasis,noccFA)            (nbasis,noccLOCAL)      (noccLOCAL,noccFA)
    !
    if(Myfragment%FAset) then
       call mem_dealloc(Myfragment%CoFA)
       call mem_dealloc(Myfragment%CvFA)
       call mem_dealloc(Myfragment%ppfockFA)
       call mem_dealloc(Myfragment%qqfockFA)
    end if
    nullify(Myfragment%CoFA,Myfragment%CvFA,Myfragment%ppfockFA,Myfragment%qqfockFA)
    call mem_alloc(Myfragment%CoFA,Myfragment%nbasis,Myfragment%noccFA)
    call dec_simple_dgemm(Myfragment%nbasis,Myfragment%noccAOS,&
         & Myfragment%noccFA,Myfragment%Co,OUred,Myfragment%CoFA,'n','n')

    ! Virtual space (same strategy as for occ space)
    call mem_alloc(Myfragment%CvFA,Myfragment%nbasis,Myfragment%nunoccFA)
    call dec_simple_dgemm(Myfragment%nbasis,Myfragment%nunoccAOS,&
         & Myfragment%nunoccFA,Myfragment%Cv,VUred,Myfragment%CvFA,'n','n')
    MyFragment%FAset=.true.

    ! Make fragment MO coefficients point to FOs (but not Fock matrices)
    call fragment_basis_point_to_FOs(Myfragment,skipfock=.true.)

    if(DECinfo%PurifyMOs) then
       call fragment_purify(Myfragment)
    end if

    ! FO Fock matrices
    call get_fragment_FA_fock(MyFragment)

    ! Make fragment pointers point to FOs for both MO coefficients and Fock matrices
    call fragment_basis_point_to_FOs(Myfragment)

    call mem_dealloc(OccOrbs)
    call mem_dealloc(VirtOrbs)
    call mem_dealloc(OU)
    call mem_dealloc(VU)
    call mem_dealloc(OUred)
    call mem_dealloc(VUred)
    call mem_dealloc(Oeival)
    call mem_dealloc(Veival)

  end subroutine fragment_adapted_transformation_matrices




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
    type(decfrag),intent(inout) :: LocalFragment
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
    integer :: i,ix,maxvirtidx,maxoccidx,noccEOS,nvirtEOS
    real(realk) :: maxvirt,maxocc

    ! EOS dimensions
    noccEOS = LocalFragment%noccEOS
    nvirtEOS = LocalFragment%nunoccEOS


    ! ================================================
    !    Find which virtual AOS orbitals to include  !
    ! ================================================
    VirtOrbs=.false.
    ! Always include all EOS orbitals
    ! (listed before remaining orbitals, see example in fragment_adapted_transformation_matrices)
    do i=1,nvirtEOS
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
    if(nvirtEOS/=LocalFragment%nunoccAOS) VirtOrbs(maxvirtidx) = .true.


    ! =================================================
    !    Find which occupied AOS orbitals to include  !
    ! =================================================
    ! Same strategy as for virtual orbitals above
    OccOrbs=.false.
    do i=1,noccEOS
       OccOrbs(i) = .true.
    end do

    maxocc=-1.0_realk
    maxoccidx=0
    do i=1,nocc

       if(.not. OccOrbs(i) .and. noccEOS/=LocalFragment%noccAOS) then
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
    if(noccEOS/=LocalFragment%noccAOS) OccOrbs(maxoccidx) = .true.


  end subroutine which_fragment_adapted_orbitals




  !> \brief Initialize the basis part of fragment, i.e. the information listed
  !> in the "expensive box" in the decfrag type definition.
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
    type(decorbital), dimension(nOcc), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(decorbital), dimension(nUnocc), intent(in) :: UnoccOrbitals
    !> Fragment to construct
    type(decfrag), intent(inout) :: fragment
    integer :: j,idx

    ! Basis info
    call atomic_fragment_basis(fragment,MyMolecule)

    ! lsitem
#ifdef VAR_MPI
    ! Quick fix such that lsitem is never constructed for global master
    ! as this will destroy the overall MPI framework.
    if(infpar%mynum/=infpar%master) then
       call build_ccfragmentlsitem(mylsitem,fragment%mylsitem,fragment%atoms_idx,&
            fragment%natoms,DECinfo%output,0)

       !Build F12 Ccabs and Ri for a fragment
       if(DECinfo%F12) then
          print *, "create_f12_cabs_and_ri_fragment_info(fragment)"
          call create_f12_cabs_and_ri_fragment_info(fragment)
       end if
    end if

    !F12-calculation F12-Fock terms
    if(DECinfo%F12) then     
       call atomic_fragment_init_f12(fragment,MyMolecule)
    endif !F12

#else
    call build_ccfragmentlsitem(mylsitem,fragment%mylsitem,fragment%atoms_idx,&
         fragment%natoms,DECinfo%output,0)

    !Build F12 Cabs and Ri for a fragment
    if(DECinfo%F12) then
       print *, "create_f12_cabs_and_ri_fragment_info(fragment)"
       call create_f12_cabs_and_ri_fragment_info(fragment)
    endif

    !F12-calculation F12-Fock terms
    if(DECinfo%F12) then     
       call atomic_fragment_init_f12(fragment,MyMolecule)
    endif !F12

#endif

    ! Basis info has now been set
    fragment%BasisInfoIsSet=.true.
    
  end subroutine atomic_fragment_init_basis_part


  subroutine create_f12_cabs_and_ri_fragment_info(fragment)
    implicit none 
    type(decfrag), intent(inout) :: fragment

    ! ============================================================
    !                F12-Specific for Cabs and Cri               !
    ! ============================================================    
    if(DECinfo%F12) then 
       call dec_get_CABS_orbitals_fragment(fragment, fragment%mylsitem)
       call dec_get_RI_orbitals_fragment(fragment, fragment%mylsitem)
    endif
  end subroutine create_f12_cabs_and_ri_fragment_info
  
  
  subroutine  dec_get_CABS_orbitals_fragment(fragment,mylsitem)
    implicit none

    !> Fragment molecule structure to be initialized
    type(decfrag), intent(inout) :: fragment
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem

    type(matrix) :: CMO_cabs
    integer :: ncabsAO,ncabs

    call determine_CABS_nbast(ncabsAO,ncabs,mylsitem%setting,DECinfo%output)
    call mat_init(CMO_cabs,nCabsAO,nCabs)

    call init_cabs()
    call build_CABS_MO(CMO_cabs,ncabsAO,mylsitem%SETTING,DECinfo%output)
    call free_cabs()

    ! NB! Memory leak need to be freed somewhere
    call mem_alloc(fragment%Ccabs,ncabsAO,nCabs)
    call mat_to_full(CMO_cabs,1.0E0_realk,fragment%Ccabs)
    call mat_free(CMO_cabs)

  end subroutine dec_get_CABS_orbitals_fragment

  subroutine  dec_get_RI_orbitals_fragment(fragment,mylsitem)
    implicit none

    !> Fragment molecule structure to be initialized
    type(decfrag), intent(inout) :: fragment
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem

    type(matrix) :: CMO_RI
    integer :: ncabsAO,ncabs,lupri

    call determine_CABS_nbast(ncabsAO,ncabs,mylsitem%setting,DECinfo%output)

    call mat_init(CMO_RI,ncabsAO,ncabsAO)

    call init_cabs()
    call build_RI_MO(CMO_RI,ncabsAO,mylsitem%SETTING,lupri)
    call free_cabs()

    ! NB! Memory leak need to be freed somewhere
    call mem_alloc(fragment%Cri,ncabsAO,ncabsAO) 
    call mat_to_full(CMO_RI,1.0E0_realk,fragment%Cri)

    call mat_free(CMO_RI)

  end subroutine dec_get_RI_orbitals_fragment


  !\ brief Purify fragment MO coefficients by (i) projecting out possible occupied
  !> components from the unoccupied MOs (and vice versa), (ii) orthogonalize orbitalsæ.
  !> \author Kasper Kristensen
  !> \date April 2013
  subroutine fragment_purify(fragment)
    implicit none
    !> Fragment where the MOs (fragment%Co, fragment%Cv, fragment%coremo) will be purified
    type(decfrag),intent(inout) :: fragment
    integer :: nXOS,nbasis,noccAOS,nunoccAOS,ncore,i
    real(realk),pointer :: XOS(:,:)

    ! Dimensions
    nbasis = fragment%nbasis
    noccAOS = fragment%noccAOS
    nunoccAOS = fragment%nunoccAOS
    ncore = fragment%ncore

    ! We purify twice to be on the safe side
    Purify: do i=1,2

       ! OCCUPIED ORBITALS
       ! *****************
       ! Note: We only consider the XOS (AOS-EOS) because the EOS orbitals must remain untouched!

       ! Extract occupied XOS orbitals
       nXOS = fragment%noccAOS - fragment%noccEOS
       if(nXOS>0) then
          call mem_alloc(XOS,nbasis,nXOS)
          call extract_XOS_orbitals_occ(Fragment,nbasis,nXOS,XOS)

          ! (i) Project out possible unoccupied components from the occupied XOS
          call project_out_MO_space(nunoccAOS,nXOS,nbasis,Fragment%Cv,fragment%S,XOS)
          ! For frozen core also project out possible core orbital components
          if(DECinfo%frozencore .and. ncore>0) then
             call project_out_MO_space(ncore,nXOS,nbasis,Fragment%coreMO,fragment%S,XOS)
          end if

          ! (ii) Orthogonalize XOS orbitals
          call orthogonalize_MOs(nXOS,nbasis,fragment%S,XOS)

          ! Put purified XOS orbitals back into fragment structure
          call put_XOS_orbitals_occ(nbasis,nXOS,XOS,Fragment)
          call mem_dealloc(XOS)
       end if


       ! UNOCCUPIED ORBITALS
       ! *******************
       ! Same procedure as for occ orbitals

       ! Extract unoccupied XOS orbitals
       nXOS = fragment%nunoccAOS - fragment%nunoccEOS
       if(nXOS>0) then
          call mem_alloc(XOS,nbasis,nXOS)
          call extract_XOS_orbitals_unocc(Fragment,nbasis,nXOS,XOS)

          ! (i) Project out possible occupied components from the unoccupied XOS
          call project_out_MO_space(noccAOS,nXOS,nbasis,Fragment%Co,fragment%S,XOS)
          ! For frozen core also project out possible core orbital components
          if(DECinfo%frozencore .and. ncore>0) then
             call project_out_MO_space(ncore,nXOS,nbasis,Fragment%coreMO,fragment%S,XOS)
          end if

          ! (ii) Orthogonalize XOS orbitals
          call orthogonalize_MOs(nXOS,nbasis,fragment%S,XOS)

          ! Put purified XOS orbitals back into fragment structure
          call put_XOS_orbitals_unocc(nbasis,nXOS,XOS,Fragment)
          call mem_dealloc(XOS)
       end if


       ! CORE ORBITALS
       ! *************
       if(ncore>0) then

          ! (i) Project out possible unoccupied components from core space
          call project_out_MO_space(nunoccAOS,ncore,nbasis,Fragment%Cv,fragment%S,fragment%coreMO)
          ! For frozen core also project out possible valence orbital components
          if(DECinfo%frozencore) then
             call project_out_MO_space(noccAOS,ncore,nbasis,fragment%Co,fragment%S,Fragment%coreMO)
          end if

          ! (ii) Orthogonalize core orbitals
          call orthogonalize_MOs(ncore,nbasis,fragment%S,fragment%coreMO)

       end if

    end do Purify


  end subroutine fragment_purify


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
    type(decorbital), dimension(nOcc), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(decorbital), dimension(nUnocc), intent(in) :: UnoccOrbitals
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> Fragment to construct
    type(decfrag), intent(inout) :: fragment
    !> Make fragment basis (MO coeff and Fock matrix for fragment)?
    logical, intent(in) :: DoBasis
    logical, dimension(nunocc) :: Unocc_list
    logical, dimension(nocc) :: Occ_list

    ! All orbitals included in fragment
    unocc_list=.true.
    occ_list=.true.

    call atomic_fragment_init_orbital_specific(MyAtom,nunocc, nocc, unocc_list, &
         & occ_list,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,fragment,DoBasis,.false.)

  end subroutine fragment_init_simulate_full


  !> \brief Merge two single fragment into one pair fragment
  !> \author Kasper Kristensen
  !> \date December 2011
  subroutine merged_fragment_init(Fragment1,Fragment2,nunocc, nocc, natoms, &
       & OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,DoBasis,pairfragment,esti)


    implicit none
    !> Fragment 1 in pair
    type(decfrag),intent(inout) :: fragment1
    !> Fragment 2 in pair
    type(decfrag),intent(inout) :: fragment2
    !> Number of occupied orbitals in full molecule
    integer, intent(in) :: nocc
    !> Number of unoccupied orbitals in full molecule
    integer, intent(in) :: nunocc
    !> Number of atoms for full molecule
    integer, intent(in) :: natoms
    !> Information about DEC occupied orbitals
    type(decorbital), dimension(nOcc), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(decorbital), dimension(nUnocc), intent(in) :: UnoccOrbitals
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> Make pair fragment basis (MO coeff and Fock matrix for fragment)?
    logical, intent(in) :: DoBasis
    !> Pair fragment to be constructed
    type(decfrag),intent(inout) :: pairfragment
    !> Use estimated fragments? (If this is the case FOs are effectively turned off)
    logical,intent(in),optional :: esti
    logical, dimension(nunocc) :: Unocc_list
    logical, dimension(nocc) :: Occ_list
    logical :: pairfrag,estimated_frags
    logical,pointer :: EOSatoms(:)
    integer :: i,j,idx
    real(realk) :: pairdist

    ! Note: When we do this, the atomic fragment has to use MO coeff and dimensions for LOCAL orbitals
    call fragment_basis_point_to_LOs(fragment1)
    call fragment_basis_point_to_LOs(fragment2)
    ! Then, if fragment-adapted orbitals are requested we direct the
    ! fragment data pointers to the fragment-adapted data rather than local.
    ! This is done at the end of pair_fragment_adapted_transformation_matrices for the pair fragment
    ! and at the end of the routine for the incoming atomic fragments.


    ! Use estimated fragments?
    estimated_frags=.false.
    if(present(esti)) then
       if(esti) estimated_frags=.true.
    end if

    pairfrag=.true.

    ! Sanity check : No overlap between EOS atoms in fragment1 and fragment2
    ! **********************************************************************

    ! We also set a logical vector for EOS atoms to be used below
    call mem_alloc(EOSatoms,natoms)
    EOSatoms=.false.

    do i=1,fragment1%nEOSatoms
       EOSatoms(fragment1%EOSatoms(i))=.true.
       do j=1,fragment2%nEOSatoms
          EOSatoms(fragment2%EOSatoms(j))=.true.
          if(fragment1%EOSatoms(i) == fragment2%EOSatoms(j)) then
             call lsquit('merged_fragment_init: Overlap between EOS atoms for fragment 1 and 2 ',DECinfo%output)
          end if
       end do
    end do

    ! just to be sure, initialize everything to zero or null
    call atomic_fragment_nullify(pairfragment)

    ! Find pair distance
    pairdist = get_distance_between_fragments(fragment1,&
         & fragment2,natoms,MyMolecule%DistanceTable)
    if(DECinfo%PL>0) write(DECinfo%output,*) '--> Initiating pair fragment...'
    if(DECinfo%PL>0) write(DECinfo%output,'(1X,a,F12.3)') 'Pair distance (Angstrom)      = '&
         &, PairDist*bohr_to_angstrom




    ! Find EOS atoms for pair: Union of EOS atoms for fragment1 and fragment2
    ! ***********************************************************************

    ! Since there is no overlap between the EOS atoms (sanity check 2),
    ! the pair EOS size is just the sum of the individual sizes
    pairfragment%nEOSatoms = fragment1%nEOSatoms + fragment2%nEOSatoms
    call mem_alloc(pairfragment%EOSatoms,pairfragment%nEOSatoms)
    if(pairfragment%nEOSatoms /= count(EOSatoms)) then
       write(DECinfo%output,*) 'Atoms in fragment1 = ', fragment1%EOSatoms
       write(DECinfo%output,*) 'Atoms in fragment2 = ', fragment2%EOSatoms
       call lsquit('merged_fragment_init: pairfragment%nEOSatoms /= count(EOSatoms)',DECinfo%output)
    end if

    ! Set indices for EOS atoms for pair
    idx=0
    do i=1,natoms
       if(EOSatoms(i)) then
          idx=idx+1
          pairfragment%EOSatoms(idx) = i
       end if
    end do


    ! Occupied AOS and unoccupied AOS for pair: Union of input fragments
    ! ******************************************************************

    occ_list=.false.
    unocc_list=.false.

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


    ! Construct pair fragment as union of input fragments
    ! ***************************************************
    ! Simply set MyAtom=0 for pairs
    pairfragment%pairfrag=.true.
    call atomic_fragment_init_orbital_specific(0,nunocc, nocc, unocc_list, &
         & occ_list,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,pairfragment,DoBasis,pairfrag)


    ! Set remaining information special for pair
    ! ******************************************
    PairFragment%EOSatoms(1) = fragment1%EOSatoms(1) ! atom 1 in pair
    ! Distance between original fragments
    PairFragment%pairdist =  pairdist
    call mem_dealloc(EOSatoms)

    ! Print out summary
    if(DECinfo%PL>0) then
       write(DECinfo%output,*)
       write(DECinfo%output,*) 'PAIR CONSTRUCTION SUMMARY'
       write(DECinfo%output,*) '*************************'
       write(DECinfo%output,'(1X,a,20i4)')  'PCS: EOS atoms in fragment 1          :', Fragment1%EOSatoms
       write(DECinfo%output,'(1X,a,20i4)')  'PCS: EOS atoms in fragment 2          :', Fragment2%EOSatoms
       write(DECinfo%output,'(1X,a,i4)')    'PCS: Number of orbitals in virt total :', PairFragment%nunoccAOS
       write(DECinfo%output,'(1X,a,i4)')    'PCS: Number of orbitals in occ total  :', PairFragment%noccAOS
       write(DECinfo%output,'(1X,a,i4)')    'PCS: Number of basis functions        :', PairFragment%nbasis
       write(DECinfo%output,'(1X,a,F12.3)') 'PCS: Pair distance (Angstrom)         :', PairDist*bohr_to_angstrom
       write(DECinfo%output,*)
    end if

    if(DECinfo%fragadapt .and. (.not. estimated_frags)  ) then
       call pair_fragment_adapted_transformation_matrices(MyMolecule,nocc,nunocc,&
            & OccOrbitals,UnoccOrbitals, Fragment1,Fragment2,PairFragment)
       ! Set atomic fragment pointers back to FOs data.
       call fragment_basis_point_to_FOs(fragment1)
       call fragment_basis_point_to_FOs(fragment2)
    end if

  end subroutine merged_fragment_init




  !> Set fragment-adapted MO coefficients for
  !> pair fragment (FragmentPQ%CoFA andFragmentPQ%CvFA).
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
    type(decorbital), dimension(nOcctot), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals for full molecule
    type(decorbital), dimension(nUnocctot), intent(in) :: UnoccOrbitals
    !> Fragment P
    type(decfrag),intent(inout) :: fragmentP
    !> Fragment Q
    type(decfrag),intent(inout) :: fragmentQ
    !> Pair fragment PQ (it is assumed that everything except fragment-adapted
    !> stuff has already been initialized).
    type(decfrag),intent(inout) :: FragmentPQ
    real(realk) :: lambdathr,lambdathr_default
    real(realk),pointer :: CoccPQ(:,:), CunoccPQ(:,:),moS(:,:),Sinv(:,:)
    real(realk),pointer :: T(:,:),lambda(:),MOtmp(:,:),tmp(:,:),tmp2(:,:),CEOS(:,:)
    logical,pointer :: whichOrbitals(:)
    integer :: noccPQ,nunoccPQ,nbasisPQ,noccPQ_FA, nunoccPQ_FA,i,mu,idx, EOSidx,j
    logical :: debugprint,keepon  ! temporary debug prints
    real(realk) :: diagdev, nondiagdev
    logical,pointer :: WhichOccP(:), WhichOccQ(:), WhichUnoccP(:), WhichUnoccQ(:)


    if(DECinfo%PL>0) then
       debugprint=.true.
    else
       debugprint=.false.
    end if
    lambdathr_default = 1.0e-3_realk

    ! Threshold used to remove certain FOs
    FragmentPQ%RejectThr = set_pairFOthr(FragmentPQ%pairdist)

    ! Dimensions of pair fragment PQ as union of spaces for atomic fragments P and Q,
    ! where only FOs with eigenvalues larger than pairFOthr in the original atomic
    ! fragments are kept.
    call mem_alloc(WhichOccP,fragmentP%noccFA)
    call mem_alloc(WhichOccQ,fragmentQ%noccFA)
    call mem_alloc(WhichUnoccP,fragmentP%nunoccFA)
    call mem_alloc(WhichUnoccQ,fragmentQ%nunoccFA)
    call get_pairFO_union(fragmentP,fragmentQ,fragmentPQ,noccPQ,nunoccPQ,nbasisPQ,&
         & WhichOccP, WhichOccQ, WhichUnoccP, WhichUnoccQ)


    ! Sanity checks
    if( (.not. fragmentP%FAset) .or. (.not. fragmentQ%FAset) ) then
       call lsquit('pair_fragment_adapted_transformation_matrices: Transformation &
            & matrices for atomic fragments are not set!',-1)
    end if
    if( (fragmentPQ%noccEOS == fragmentPQ%noccAOS) .or. &
         & (fragmentPQ%nunoccEOS == fragmentPQ%nunoccAOS) ) then
       ! Special case: No FA orbitals to consider.
       call pair_fragment_adapted_transformation_matrices_justEOS(fragmentP,fragmentQ,fragmentPQ)
       call mem_dealloc(WhichOccP)
       call mem_dealloc(WhichOccQ)
       call mem_dealloc(WhichUnoccP)
       call mem_dealloc(WhichUnoccQ)
       return
    end if


    ! *****************************************************************************************
    ! *    Fragment-adapted orbitals (FOs) for pair PQ are set by the following procedure:   *
    ! *****************************************************************************************
    !
    ! 1. Set up MO coefficient matrix using union of existing FOs for atomic fragments P and Q
    !    (but where only atomic fragment FOs with eigenvalues larger than fragmentPQ%rejecthr
    !     are kept as described above).
    ! 2. Project onto local AOS for pair (cleanup).
    ! 3. Project out orbitals assignsed to P or Q.
    ! 4. Setup MO overlap matrix for resulting orbitals.
    ! 5. Diagonalize MO overlap matrix: SMO = T^T lambda T
    ! 6. Each lambda (and its associated eigenvector) represents a FO
    !    --> Throw away FOs where lambda<lambdathr_default.
    ! 7. Calculate FOs in terms of lambda eigenvalues and T eigenvectors.
    !
    ! (the purpose and more details for each step is given below)
    !


    ! 1. Set pair PQ MO-coefficients by simply copying P and Q coefficients
    ! *********************************************************************
    call mem_alloc(CoccPQ,NbasisPQ,noccPQ)
    call mem_alloc(CunoccPQ,NbasisPQ,nunoccPQ)
    call set_pair_fragment_adapted_redundant_orbitals(MyMolecule,FragmentP,FragmentQ,&
         & FragmentPQ,noccPQ,nunoccPQ,WhichOccP,WhichOccQ,WhichUnoccP,WhichUnoccQ,CoccPQ,CunoccPQ)
    ! Note: At this stage the PQ orbitals are redundant and not orthogonal.
    !       Furthermore, the EOS orbitals have been taken out of the orbital pool.
    !       However, FOs originally used for fragment P may still contain
    !       components of EOS orbitals on fragment Q (and vice versa).
    !       These EOS components must be projected out.
    !
    ! We consider first the occupied space and then do the same for the unoccupied space.



    ! 2. Project onto local AOS
    ! *************************
    ! We need to project against the local AOS since the current CoccPQ coefficients
    ! are simply the FO coefficients for fragments P and Q. These FO coefficients have been
    ! determined by diagonalizing the correlation density matrix expressed in terms of local orbitals.
    ! However, in general P and Q have different atomic extents, which again are different
    ! from the atomic extent for pair fragment PQ. This means that the FOs from P and Q may
    ! contain (artificial) components outside the local AOS for pair fragment PQ.
    ! Therefore we need to project against the local AOS for pair fragment PQ to ensure that
    ! CoccPQ really does not contain components outside the local AOS for pair PQ.
    !
    ! The MOs projected against the local AOS for PQ are:
    !
    ! CoccPQ --> PROJ CoccPQ
    !
    ! PROJ = Clocal Sinv Clocal^T aoS    (projector)
    !
    ! where Clocal are the local MOs for pair PQ (fitted in the atomic extent for PQ),
    ! Sinv is the inverse overlap matrix for the local MOs for pair PQ (the identity matrix in the
    ! limit where the atomic extent for PQ contains all atomic basis functions in the molecule),
    ! aoS is the AO overlap matrix.
    call project_onto_MO_space(fragmentPQ%noccAOS,noccPQ,nbasisPQ,fragmentPQ%Co,&
         & fragmentPQ%S,CoccPQ)


    ! 3. Project out EOS components
    ! *****************************
    !
    ! MO coefficients where EOS orbital components are projected out:
    !
    ! CoccPQ = CoccPQ - CEOS Sinv CEOS^T aoS CoccPQ
    !
    ! where CoccPQ are the MO coefficients determined above
    ! (which currently contain EOS components),
    ! CEOS are the EOS MO coefficients for the pair fragment (orbitals assigned to P or Q),
    ! and Sinv is the inverse overlap matrix for EOS orbitals.


    ! Extract EOS indices
    call mem_alloc(CEOS,nbasisPQ,fragmentPQ%noccEOS)
    do i=1,fragmentPQ%noccEOS
       CEOS(:,i) = fragmentPQ%Co(:,fragmentPQ%idxo(i))
    end do

    ! Project out EOS componets
    call project_out_MO_space(fragmentPQ%noccEOS,noccPQ,nbasisPQ,CEOS,fragmentPQ%S,CoccPQ)
    call mem_dealloc(CEOS)


    ! 4. Setup MO overlap matrix
    ! **************************
    call mem_alloc(moS,noccPQ,noccPQ)
    call dec_simple_basis_transform1(nbasisPQ,noccPQ,CoccPQ,fragmentPQ%S,moS)



    ! 5. Diagonalize MO overlap matrix: moS = T^T lambda T
    ! ****************************************************
    call mem_alloc(lambda,noccPQ)
    call mem_alloc(T,noccPQ,noccPQ)
    call solve_eigenvalue_problem_unitoverlap(noccPQ,moS,lambda,T)
    lambda=abs(lambda)
    call mem_alloc(whichorbitals,noccPQ)



    ! 6. Throw away eigenvalues smaller than threshold
    ! ************************************************

    ! Determine which orbitals to include
    keepon=.true.
    lambdathr = lambdathr_default
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
    ! Number of FOs for pair (not including EOS orbitals)
    noccPQ_FA = count(whichorbitals)

    if(debugprint) then
       do i=1,noccPQ
          print *, 'lambda occ', i, lambda(i)
       end do
    end if



    ! 7. Calculate FOs from lambda and T
    ! ***********************************

    ! Orthogonal fragment-adapted molecular orbitals (psi) are given as:
    !
    ! psi(i) = lambda(i)^{-1/2} * sum_{mu} (CoccPQ T)_{mu i} chi(mu)
    !
    ! where chi(mu) is atomic orbital "mu" in pair atomic fragment extent.
    !

    ! MOtmp = CoccPQ T
    call mem_alloc(MOtmp,nbasisPQ,noccPQ)
    call dec_simple_dgemm(nbasisPQ,noccPQ,noccPQ,CoccPQ,T,MOtmp,'n','n')

    ! Final MO coefficents (EOS orbitals + FA orbitals) are stored in fragmentPQ%CoFA and found by:
    ! (i)  Copying existing EOS MO coefficients into the first "noccEOS" columns
    ! (ii) Copying the non-redundant orthogonal orbitals (defined by whichorbitals)
    !      in MOtmp into the remaining columns.

    ! Total number of pair orbitals (EOS +FA)
    fragmentPQ%noccFA = fragmentPQ%noccEOS + noccPQ_FA

    call mem_alloc(fragmentPQ%CoFA,nbasisPQ,fragmentPQ%noccFA)
    fragmentPQ%CoFA = 0.0_realk

    ! (i) Copy EOS
    ! Note: For FA orbitals we always put EOS orbitals before remaining orbitals,
    !      this is not the case for local orbitals (see fragment_adapted_transformation_matrices)
    do i=1,fragmentPQ%noccEOS
       ! EOS index for local orbitals in AOS list for local orbitals
       EOSidx = fragmentPQ%idxo(i)
       do mu=1,nbasisPQ
          fragmentPQ%CoFA(mu,i) = fragmentPQ%Co(mu,EOSidx)
       end do
    end do
    idx= fragmentPQ%noccEOS  ! counter

    ! (ii) Copy FA orbitals (also divide my lambda^{-1/2} as shown above)
    do i=1,noccPQ
       if(whichorbitals(i)) then
          idx = idx+1
          do mu=1,nbasisPQ
             fragmentPQ%CoFA(mu,idx) = MOtmp(mu,i)/sqrt(lambda(i))
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
            &  fragmentPQ%CoFA,fragmentPQ%S,moS)
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


    ! 2
    call project_onto_MO_space(fragmentPQ%nunoccAOS,nunoccPQ,nbasisPQ,fragmentPQ%Cv,fragmentPQ%S,CunoccPQ)


    ! 3
    call mem_alloc(CEOS,nbasisPQ,fragmentPQ%nunoccEOS)
    do i=1,fragmentPQ%nunoccEOS
       CEOS(:,i) = fragmentPQ%Cv(:,fragmentPQ%idxu(i))
    end do
    call project_out_MO_space(fragmentPQ%nunoccEOS,nunoccPQ,nbasisPQ,CEOS,fragmentPQ%S,CunoccPQ)
    call mem_dealloc(CEOS)

    ! 4
    call mem_alloc(moS,nunoccPQ,nunoccPQ)
    call dec_simple_basis_transform1(nbasisPQ,nunoccPQ,CunoccPQ,fragmentPQ%S,moS)


    ! 5
    call mem_alloc(lambda,nunoccPQ)
    call mem_alloc(T,nunoccPQ,nunoccPQ)
    call solve_eigenvalue_problem_unitoverlap(nunoccPQ,moS,lambda,T)
    lambda=abs(lambda)
    call mem_alloc(whichorbitals,nunoccPQ)


    ! 6
    keepon=.true.
    lambdathr = lambdathr_default
    UnoccCheck: do while(keepon)
       whichorbitals=.false.
       do i=1,nunoccPQ
          if(lambda(i) > lambdathr) whichorbitals(i)=.true.
       end do

       if(count(whichorbitals)>fragmentPQ%nunoccAOS - fragmentPQ%nunoccEOS) then
          lambdathr = lambdathr*10.0_realk ! increase lambda threshold
       else ! done
          keepon=.false.
       end if

    end do UnoccCheck
    nunoccPQ_FA = count(whichorbitals)
    if(debugprint) then
       do i=1,nunoccPQ
          print *, 'lambda unocc', i, lambda(i)
       end do
    end if



    ! 7
    call mem_alloc(MOtmp,nbasisPQ,nunoccPQ)
    call dec_simple_dgemm(nbasisPQ,nunoccPQ,nunoccPQ,CunoccPQ,T,MOtmp,'n','n')
    fragmentPQ%nunoccFA = fragmentPQ%nunoccEOS + nunoccPQ_FA
    call mem_alloc(fragmentPQ%CvFA,nbasisPQ,fragmentPQ%nunoccFA)
    fragmentPQ%CvFA = 0.0_realk
    do i=1,fragmentPQ%nunoccEOS
       EOSidx = fragmentPQ%idxu(i)
       do mu=1,nbasisPQ
          fragmentPQ%CvFA(mu,i) = fragmentPQ%Cv(mu,EOSidx)
       end do
    end do
    idx= fragmentPQ%nunoccEOS  ! counter
    do i=1,nunoccPQ
       if(whichorbitals(i)) then
          idx = idx+1
          do mu=1,nbasisPQ
             fragmentPQ%CvFA(mu,idx) = MOtmp(mu,i)/sqrt(lambda(i))
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
            &  fragmentPQ%CvFA,fragmentPQ%S,moS)
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
    call mem_dealloc(WhichOccP)
    call mem_dealloc(WhichOccQ)
    call mem_dealloc(WhichUnoccP)
    call mem_dealloc(WhichUnoccQ)


    ! Transformation matrices have been set!
    fragmentPQ%FAset=.true.

    ! Make fragment MO coefficients point to FOs (but not Fock matrices)
    call fragment_basis_point_to_FOs(FragmentPQ,skipfock=.true.)

    if(DECinfo%PurifyMOs) then
       call fragment_purify(FragmentPQ)
    end if

    ! FO Fock matrices
    call get_fragment_FA_fock(fragmentPQ)

    ! Make pointers for MO coeff and Fock matrices point to FO quantities
    call fragment_basis_point_to_FOs(fragmentPQ)

  end subroutine pair_fragment_adapted_transformation_matrices


  !> \brief Allocate and calculate occ-occ and virt-virt blocks of Fock matrix in FO basis
  !> as CoccFA^T FAO CoccFA
  !> \author Kasper Kristensen
  !> \date November 2013
  subroutine get_fragment_FA_fock(myfragment)
    implicit none
    type(decfrag),intent(inout) :: myfragment
   
    ! Occ space, check the allocation status of the corresponding arrays
    !********************************************************************
    if(.not.associated(Myfragment%ppfockFA))then

      call mem_alloc(Myfragment%ppfockFA,Myfragment%noccFA,Myfragment%noccFA)

    else

      if(size(Myfragment%ppfockFA)/=Myfragment%noccFA**2)then

        call mem_dealloc(Myfragment%ppfockFA)
        call mem_alloc(Myfragment%ppfockFA,Myfragment%noccFA,Myfragment%noccFA)

      endif

    endif

    call dec_simple_basis_transform1(Myfragment%nbasis,Myfragment%noccFA,&
         & Myfragment%CoFA,Myfragment%fock,Myfragment%ppfockFA)




    ! Virt space, same check
    !************************
    if(.not.associated(Myfragment%qqfockFA))then

      call mem_alloc(Myfragment%qqfockFA,Myfragment%nunoccFA,Myfragment%nunoccFA)

    else

      if(size(Myfragment%qqfockFA)/=Myfragment%nunoccFA**2)then

        call mem_dealloc(Myfragment%qqfockFA)
        call mem_alloc(Myfragment%qqfockFA,Myfragment%nunoccFA,Myfragment%nunoccFA)

      endif

    endif

    call dec_simple_basis_transform1(Myfragment%nbasis,Myfragment%nunoccFA,&
         & Myfragment%CvFA,Myfragment%fock,Myfragment%qqfockFA)
 
  end subroutine get_fragment_FA_fock


  !> Special case for pair_fragment_adapted_transformation_matrices where
  !> there are only EOS orbitals.
  !> ONLY relevant for small test systems.
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine pair_fragment_adapted_transformation_matrices_justEOS(fragmentP,&
       & fragmentQ, fragmentPQ)
    implicit none

    !> Fragment P
    type(decfrag),intent(inout) :: fragmentP
    !> Fragment Q
    type(decfrag),intent(inout) :: fragmentQ
    !> Pair fragment PQ (it is assumed that everything except fragment-adapted
    !> stuff has already been initialized).
    type(decfrag),intent(inout) :: FragmentPQ

    ! Occupied space
    fragmentPQ%noccFA = fragmentPQ%noccAOS
    call mem_alloc(fragmentPQ%CoFA,fragmentPQ%nbasis,fragmentPQ%noccFA)
    fragmentPQ%CoFA = fragmentPQ%Co

    ! Unoccupied space
    fragmentPQ%nunoccFA = fragmentPQ%nunoccAOS
    call mem_alloc(fragmentPQ%CvFA,fragmentPQ%nbasis,fragmentPQ%nunoccFA)
    fragmentPQ%CvFA = fragmentPQ%Cv

    ! FO Fock matrices
    call get_fragment_FA_fock(fragmentPQ)

    fragmentPQ%FAset=.true.

    ! Make pointers for MO coeff and Fock matrices point to FO quantities
    call fragment_basis_point_to_FOs(fragmentPQ)

  end subroutine pair_fragment_adapted_transformation_matrices_justEOS



  !> \brief Get dimensions for fragment-adapted AOS for pair fragment.
  !> Only to be used for job scheduling, pair fragment is deleted again.
  !> \author Kasper Kristensen
  !> \date February 2013
  subroutine get_fragment_adapted_dimensions_for_pair(Fragment1,Fragment2,nunoccFULL, noccFULL, natoms, &
       & OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,&
       & noccFA,nunoccFA,nbasisFA)

    implicit none
    !> Fragment 1 in pair
    type(decfrag),intent(inout) :: fragment1
    !> Fragment 2 in pair
    type(decfrag),intent(inout) :: fragment2
    !> Number of occupied orbitals in full molecule
    integer, intent(in) :: noccFULL
    !> Number of unoccupied orbitals in full molecule
    integer, intent(in) :: nunoccFULL
    !> Number of atoms for full molecule
    integer, intent(in) :: natoms
    !> Information about DEC occupied orbitals
    type(decorbital), dimension(noccFULL), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(decorbital), dimension(nunoccFULL), intent(in) :: UnoccOrbitals
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
    type(decfrag) :: pairfragment


    ! Init pair fragment
    ! ------------------
    call merged_fragment_init(Fragment1,Fragment2,nunoccFULL, noccFULL, natoms, &
         & OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,.true.,pairfragment)

    noccFA = pairfragment%noccFA
    nunoccFA = pairfragment%nunoccFA
    nbasisFA = pairfragment%nbasis

    call atomic_fragment_free(pairfragment)

  end subroutine get_fragment_adapted_dimensions_for_pair



  !> \brief Initialize atomic fragment extent info for fragment
  !> (fragment%natoms, fragment%nbasis, fragment%atoms_idx)
  !> and orbital info in DEC format (fragment%occAOSorb and fragment%unoccAOSorb)
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine init_atomic_fragment_extent(OccOrbitals,UnoccOrbitals,MyMolecule,fragment)

    implicit none

    !> FUll molecular info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Information about DEC occupied orbitals
    type(decorbital), dimension(MyMolecule%nocc), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(decorbital), dimension(MyMolecule%nunocc), intent(in) :: UnoccOrbitals
    !> Fragment info
    type(decfrag), intent(inout) :: fragment
    integer :: i,j,idx,cntr
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
    fragment%RmaxAE = 0.0E0_realk
    fragment%RaveAE = 0.0E0_realk
    fragment%RsdvAE = 0.0E0_realk
    fragment%natoms = count(which_atoms)
    call mem_alloc(fragment%atoms_idx,Fragment%natoms)
    idx  = 0
    cntr = 0
    do i=1,MyMolecule%natoms
       if(which_atoms(i)) then
          idx=idx+1
          fragment%atoms_idx(idx)=i

          !set RmaxAE
          do j=i+1,MyMolecule%natoms
            if(which_atoms(j)) then
              fragment%RmaxAE = max(fragment%RmaxAE,MyMolecule%DistanceTable(j,i))
              fragment%RaveAE = fragment%RaveAE + MyMolecule%DistanceTable(j,i)
              cntr = cntr + 1
            endif
          enddo

       end if
    end do
    fragment%RaveAE = fragment%RaveAE / float( cntr )
    do i=1,MyMolecule%natoms
       if(which_atoms(i)) then
          do j=i+1,MyMolecule%natoms
            if(which_atoms(j)) then
              fragment%RsdvAE = fragment%RsdvAE + ( MyMolecule%DistanceTable(j,i)- fragment%RaveAE)**2
            endif
          enddo

       end if
    end do
    fragment%RsdvAE = ( fragment%RsdvAE / float( cntr - 1 ) )**(0.5E0_realk)
    call mem_dealloc(which_atoms)

    ! count number of basis functions on selected atoms
    fragment%nbasis = 0
    do i=1,fragment%natoms
       fragment%nbasis = fragment%nbasis &
            & + MyMolecule%atom_size(fragment%atoms_idx(i))
    end do

    ! Set basis function indices
    call mem_alloc(fragment%basis_idx,fragment%nbasis)
    idx=0
    do i=1,fragment%natoms
       do j=1,MyMolecule%atom_size(fragment%atoms_idx(i))
          idx=idx+1
          ! Basis function index
          fragment%basis_idx(idx) = MyMolecule%atom_start(fragment%atoms_idx(i)) + j-1
       end do
    end do

    ! Sanity check
    if(idx/=fragment%nbasis) then
       print *, 'Counter,nbasis: ', idx,fragment%nbasis
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
    type(decfrag), intent(inout) :: fragment
    type(fullmolecule), intent(in) :: MyMolecule
    type(array2) :: S,tmp1,tmp2
    real(realk), pointer :: correct_vector_moS(:), approximated_orbital(:)
    integer :: i,j,idx,atom_k,nocc,nunocc,nbasis,natoms,k,bas_offset,offset, &
         full_orb_idx
    integer, dimension(2) :: dims, dimsAO, dimsMO
    logical,pointer :: which_atoms(:)

    ! allocate C^o(nbasis,occ) C^v(nbasis,unocc)
    call mem_alloc(fragment%CoLOC, fragment%nbasis,  fragment%noccLOC   )
    call mem_alloc(fragment%CvLOC, fragment%nbasis,  fragment%nunoccLOC )
    ! Core
    fragment%CoLOC=0.0E0_realk
    fragment%CvLOC=0.0E0_realk

    call mem_alloc(fragment%CoreMO, fragment%nbasis, fragment%ncore)
    fragment%CoreMO=0.0E0_realk

    ! truncate basis to this set of atoms
    nocc = MyMolecule%nocc
    nunocc = MyMolecule%nunocc
    nbasis = MyMolecule%nbasis
    natoms = MyMolecule%natoms
    dimsAO(1)=nbasis
    dimsAO(2)=nbasis

    ! Overlap matrix for fragment
    call mem_alloc(fragment%S,fragment%nbasis,fragment%nbasis)
    call adjust_square_matrix(MyMolecule%overlap,fragment%S,fragment%atoms_idx, &
         & MyMolecule%atom_size,MyMolecule%atom_start,MyMolecule%atom_end, &
         & nbasis,natoms,fragment%nbasis,Fragment%natoms)

    FitOrbitalsForFragment: if(DECinfo%FitOrbitals) then ! fit orbitals for fragment to exact orbitals

       ! fit orbitals
       dimsMO(1) = nbasis
       dimsMO(2) = nocc
       dims(1) = nocc
       dims(2) = nbasis
       S = array2_init(dims)
       call mem_alloc(correct_vector_moS,fragment%nbasis)
       call mem_alloc(approximated_orbital,fragment%nbasis)

       ! Fragment Co
       ! half transform overlap
       tmp1 = array2_init(dimsMO,MyMolecule%Co)
       tmp2 = array2_init(dimsAO,MyMolecule%overlap)
       call array2_matmul(tmp1,tmp2,S,'T','N',1E0_realk,0E0_realk)

       ! Occ orbitals (only valence if frozen core approx is used)
       do i=1,fragment%noccLOC

          full_orb_idx = fragment%occAOSidx(i)

          ! for each orbital truncate half transformed overlap orbital
          correct_vector_moS = 0.0E0_realk
          bas_offset = 1
          do k=1,Fragment%natoms
             atom_k = fragment%atoms_idx(k)
             correct_vector_moS(bas_offset:bas_offset + MyMolecule%atom_size(atom_k) - 1) = &
                  S%val(full_orb_idx,MyMolecule%atom_start(atom_k):MyMolecule%atom_end(atom_k))
             bas_offset = bas_offset + MyMolecule%atom_size(atom_k)
          end do

          ! fit orbital
          approximated_orbital = 0.0E0_realk
          call solve_linear_equations(fragment%S,approximated_orbital, &
               correct_vector_moS,fragment%nbasis)
          fragment%CoLOC(:,i) = approximated_orbital

       end do


       ! Core orbitals
       do i=1,fragment%ncore

          full_orb_idx = fragment%coreidx(i)

          ! for each orbital truncate half transformed overlap orbital
          correct_vector_moS = 0.0E0_realk
          bas_offset = 1
          do k=1,Fragment%natoms
             atom_k = fragment%atoms_idx(k)
             correct_vector_moS(bas_offset:bas_offset + MyMolecule%atom_size(atom_k) - 1) = &
                  S%val(full_orb_idx,MyMolecule%atom_start(atom_k):MyMolecule%atom_end(atom_k))
             bas_offset = bas_offset + MyMolecule%atom_size(atom_k)
          end do

          ! fit orbital
          approximated_orbital = 0.0E0_realk
          call solve_linear_equations(fragment%S,approximated_orbital, &
               correct_vector_moS,fragment%nbasis)
          fragment%coreMO(:,i) = approximated_orbital

       end do


          call array2_free(tmp1)
          call array2_free(S)
          dimsMO(1) = nbasis
          dimsMO(2) = nunocc
          dims(1) = nunocc
          dims(2) = nbasis
          S = array2_init(dims)
          tmp1 = array2_init(dimsMO,MyMolecule%Cv)


          ! Fragmant Cv
          call array2_matmul(tmp1,tmp2,S,'T','N',1E0_realk,0E0_realk)

          do i=1,fragment%nunoccLOC

             full_orb_idx = fragment%unoccAOSidx(i)

             correct_vector_moS = 0.0E0_realk
             bas_offset = 1
             do k=1,Fragment%natoms
                atom_k = fragment%atoms_idx(k)
                correct_vector_moS(bas_offset:bas_offset + MyMolecule%atom_size(atom_k) - 1) = &
                     S%val(full_orb_idx,MyMolecule%atom_start(atom_k):MyMolecule%atom_end(atom_k))
                bas_offset = bas_offset + MyMolecule%atom_size(atom_k)
             end do

             ! fit orbital
             approximated_orbital = 0.0E0_realk
             call solve_linear_equations(fragment%S,approximated_orbital, &
                  correct_vector_moS,fragment%nbasis)
             fragment%CvLOC(:,i) = approximated_orbital
          end do

       ! remove stuff
       call mem_dealloc(correct_vector_moS)
       call mem_dealloc(approximated_orbital)
       call array2_free(tmp1)
       call array2_free(tmp2)
       call array2_free(S)



    else ! Simply extract the exact MO coefficients for atoms in fragment without fitting

       if(DECinfo%frozencore) then
          call lsquit('atomic_fragment_basis: Needs implementation!',-1)
       end if

       ! Fragment Co
       call adjust_basis_matrix(MyMolecule%Co,fragment%CoLOC,fragment%occAOSidx, &
            fragment%atoms_idx,MyMolecule%atom_size,MyMolecule%atom_start,MyMolecule%atom_end,nbasis,nocc,natoms, &
            fragment%nbasis,fragment%noccLOC,Fragment%natoms)

       ! Fragment Cv
       call adjust_basis_matrix(MyMolecule%Cv,fragment%CvLOC,fragment%unoccAOSidx, &
            fragment%atoms_idx,MyMolecule%atom_size,MyMolecule%atom_start,MyMolecule%atom_end,nbasis,nunocc,natoms, &
            fragment%nbasis,fragment%nunoccLOC,Fragment%natoms)

    end if FitOrbitalsForFragment


    ! KK: Purification can be problematic for local orbitals in the context of fragment optimization,
    ! so currently it is only used in combination with fragment-adapted orbitals.
!!$    ! Purify fragment MO coefficients
!!$    ! (not for fragment-adapted orbitals since these are automatically orthogonal
!!$    !  by virtue of the purification of the local orbitals).
!!$    if(DECinfo%PurifyMOs) then
!!$       call fragment_purify(fragment)
!!$    end if


    ! adjust fock matrix in ao basis
    call mem_alloc(fragment%fock,fragment%nbasis,fragment%nbasis)
    fragment%fock=0.0E0_realk
    call adjust_square_matrix(MyMolecule%fock,fragment%fock,fragment%atoms_idx,MyMolecule%atom_size, &
         MyMolecule%atom_start,MyMolecule%atom_end,MyMolecule%nbasis,MyMolecule%natoms, &
         fragment%nbasis,Fragment%natoms)


    ! adjust fock matrices in mo basis
    ! --------------------------------

    ! For fragment-adapted orbitals, MO Fock matrix is set elsewhere (fragment_adapted_transformation_matrices)

       ! Occ-occ block  (valence-valence for frozen core)
       call mem_alloc(fragment%ppfockLOC,fragment%noccLOC,fragment%noccLOC)
       call dec_simple_basis_transform1(fragment%nbasis,fragment%noccLOC,&
            & fragment%CoLOC,fragment%fock,fragment%ppfockLOC)

       ! Virtual-virtual block
       call mem_alloc(fragment%qqfockLOC,fragment%nunoccLOC,fragment%nunoccLOC)
       call dec_simple_basis_transform1(fragment%nbasis,fragment%nunoccLOC,&
            & fragment%CvLOC,fragment%fock,fragment%qqfockLOC)

    ! Core-core block
    if(fragment%ncore>0) then
       call mem_alloc(fragment%ccfock,fragment%ncore,fragment%ncore)
       call dec_simple_basis_transform1(fragment%nbasis,fragment%ncore,&
            & fragment%coreMO,fragment%fock,fragment%ccfock)
    end if

    ! Make MO coeff and Fock matrices point to local orbital quantities unless we use fragment-adapted
    ! orbitals
    if(fragment%fragmentadapted .and. fragment%FAset) then

       if(.not.associated(fragment%ppfock).or..not.associated(fragment%qqfock))then
         call get_fragment_FA_fock(fragment)
       endif
       ! Fragment-adapted
       call fragment_basis_point_to_FOs(Fragment)

    else
       ! Local orbitals
       call fragment_basis_point_to_LOs(Fragment)
    end if

  end subroutine atomic_fragment_basis




  !> \brief Indices of target orbitals
  !> \author Marcin Ziolkowski
  !> \param fragment Atomic fragment/atomic pair fragment
  subroutine target_indices(fragment)

    implicit none
    type(decfrag), intent(inout) :: fragment
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
  !> for all atomic fragments (not pair fragments).
  !> Furthermore, the file atomicfragmentsdone.info containing the job list info
  !> (e.g. which fragments have been calculated at this stage) is updated.
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine add_fragment_to_file(fragment,jobs)
    implicit none
    !> Fragment to add to common fragment file
    type(decfrag),intent(inout) :: fragment
    !> Job list
    type(joblist),intent(in) :: jobs
    character(len=40) :: FileName
    integer :: funit,idx,i
    logical :: file_exist

    ! Sanity check: Job has been done for fragment in question
    do i=1,jobs%njobs
       if(jobs%dofragopt(i) .and. jobs%atom1(i)==fragment%EOSatoms(1)) then
          if(.not. jobs%jobsdone(i)) then
             call lsquit('add_fragment_to_file: Fragment calculation is not done!',-1)
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

    ! Init stuff. Yes, I know this hardcored funit is not pretty, but I need to circumvent lsopen...
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
    type(decfrag),intent(inout) :: fragment
    integer :: i


    write(wunit) fragment%EOSatoms(1)

    ! Occupied AOS orbitals
    write(wunit) fragment%noccAOS
    write(wunit) fragment%occAOSidx

    ! Unoccupied AOS orbitals
    write(wunit) fragment%nunoccAOS
    write(wunit) fragment%unoccAOSidx


    ! EOS atom(s)
    write(wunit) fragment%nEOSatoms
    do i=1,fragment%nEOSatoms
       write(wunit) fragment%EOSatoms(i)
    end do

    ! Energies
    write(wunit) fragment%energies

    ! Correlation density matrices
    write(wunit) fragment%CDset
    write(wunit) fragment%FAset
    write(wunit) fragment%noccFA
    write(wunit) fragment%nunoccFA
    if(fragment%CDset) then
       write(wunit) fragment%occmat
       write(wunit) fragment%virtmat
       write(wunit) fragment%RejectThr
    end if
    if(fragment%FAset) then
       write(wunit) fragment%CoFA
       write(wunit) fragment%CvFA
       write(wunit) fragment%CDocceival
       write(wunit) fragment%CDunocceival
    end if

  end subroutine fragment_write_data



  !> \brief Read fragment info for the fragments which are done
  !> from file atomicfragments.info. Used for restart.
  !> The simple logical list of which fragments are done is read from
  !> file atomicfragmentsdone.info.
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine restart_atomic_fragments_from_file(natoms,MyMolecule,MyLsitem,OccOrbitals,UnoccOrbitals,&
       & DoBasis,fragments,jobs)
    implicit none
    !> Number of atoms in the full molecule
    integer,intent(in) :: natoms
    !> Full molecule info
    type(fullmolecule),intent(in)  :: MyMolecule
    !> LSitem info
    type(lsitem),intent(inout)        :: MyLsitem
    !> Occupied orbitals in DEC format
    type(decorbital),intent(in)     :: OccOrbitals(MyMolecule%nocc)
    !> Unoccupied orbitals in DEC format
    type(decorbital),intent(in)     :: UnoccOrbitals(MyMolecule%nunocc)
    !> Construct Fock matrix and MO coeff for fragments?
    logical,intent(in) :: DoBasis
    !> Atomic fragments. NOTE: Only those fragments specified by bookkeeping will be initialized here.
    type(decfrag), intent(inout),dimension(natoms) :: fragments
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
    if(.not. file_exist) then
       write(DECinfo%output,*) 'You requested DEC atomic fragment restart but no restart files exist!'
       write(DECinfo%output,*) '--> I will calculate all atomic fragments from scratch...'
       return
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
       call lsquit('restart_atomic_fragments_from_file: File atomicfragments.info &
            & does not exist!',-1)
    end if

    ! Open file
    call lsopen(funit,FileName,'OLD','UNFORMATTED')

    write(*,             '(" Restarting with ",I4," converged atomic fragments")')ndone
    write(DECinfo%output,'(" Restarting with ",I4," converged atomic fragments")')ndone

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
          call lsquit('restart_atomic_fragments_from_file: Atom not found in job list!',-1)
       end if

       ! Consistency check that atom MyAtom (represented by job list index "idx") is
       ! indeed in the list of finished jobs.
       if(.not. jobs%jobsdone(idx)) then
          print *, 'MyAtom = ', MyAtom
          print *, 'idx in job list = ', idx
          print *, 'jobsdone=', jobs%jobsdone
          call lsquit('restart_atomic_fragments_from_file: Central atom read from file &
               & is not in bookkeeping list',-1)
       end if

       !F12 restart from file
       if(DECinfo%F12) then     
          print *, "Restart from F12 file" 
          call atomic_fragment_init_f12(fragments(MyAtom),MyMolecule)
       endif
    
       call fragment_read_data(funit,fragments(MyAtom),&
            & OccOrbitals,UnoccOrbitals,MyMolecule,Mylsitem,DoBasis)

    end do

    call lsclose(funit,'KEEP')
 
  end subroutine restart_atomic_fragments_from_file



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
    type(decorbital),intent(in)     :: OccOrbitals(MyMolecule%nocc)
    type(decorbital),intent(in)     :: UnoccOrbitals(MyMolecule%nunocc)
    type(decfrag),intent(inout)        :: fragment
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
    call mem_alloc(occ_list,MyMolecule%nocc)
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
    call mem_alloc(virt_list,MyMolecule%nunocc)
    virt_list=.false.
    do i=1,nvirtAOS
       virt_list(virtAOSidx(i)) = .true.
    end do


    if(DECinfo%convert64to32) then
       call read_64bit_to_32bit(runit,fragment%nEOSatoms)
    elseif(DECinfo%convert32to64) then
       call read_32bit_to_64bit(runit,fragment%nEOSatoms)
    else
       read(runit) fragment%nEOSatoms
    end if

    call mem_alloc(fragment%EOSatoms,fragment%nEOSatoms)
    do i=1,fragment%nEOSatoms
       if(DECinfo%convert64to32) then
          call read_64bit_to_32bit(runit,fragment%EOSatoms(i))
       elseif(DECinfo%convert32to64) then
          call read_32bit_to_64bit(runit,fragment%EOSatoms(i))
       else
          read(runit) fragment%EOSatoms(i)
       end if
    end do

    ! Initialize fragment
    call atomic_fragment_init_orbital_specific(MyAtom,MyMolecule%nunocc, MyMolecule%nocc,&
         & virt_list,occ_list,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,fragment,DoBasis,.false.)


    ! Fragment energies
    read(runit) fragment%energies


    ! Correlation density matrices and fragment-adapted orbitals
    if(DECinfo%convert64to32) then
       call read_64bit_to_32bit(runit,fragment%CDset)
       call read_64bit_to_32bit(runit,fragment%FAset)
       call read_64bit_to_32bit(runit,fragment%noccFA)
       call read_64bit_to_32bit(runit,fragment%nunoccFA)
    elseif(DECinfo%convert32to64) then
       call read_32bit_to_64bit(runit,fragment%CDset)
       call read_32bit_to_64bit(runit,fragment%FAset)
       call read_32bit_to_64bit(runit,fragment%noccFA)
       call read_32bit_to_64bit(runit,fragment%nunoccFA)
    else
       read(runit) fragment%CDset
       read(runit) fragment%FAset
       read(runit) fragment%noccFA
       read(runit) fragment%nunoccFA
    end if


    ! Correlation density
    if(fragment%CDset) then
       call mem_alloc(Fragment%OccMat,Fragment%noccAOS,Fragment%noccAOS)
       call mem_alloc(Fragment%VirtMat,Fragment%nunoccAOS,Fragment%nunoccAOS)
       read(runit) fragment%occmat
       read(runit) fragment%virtmat
       read(runit) fragment%RejectThr
    end if

    ! Fragment-adapted orbitals
    if(fragment%FAset) then
       call mem_alloc(Fragment%CoFA,Fragment%nbasis,Fragment%noccFA)
       call mem_alloc(Fragment%CvFA,Fragment%nbasis,Fragment%nunoccFA)
       call mem_alloc(Fragment%CDocceival,Fragment%noccFA)
       call mem_alloc(Fragment%CDunocceival,Fragment%nunoccFA)
       read(runit) fragment%CoFA
       read(runit) fragment%CvFA
       read(runit) fragment%CDocceival
       read(runit) fragment%CDunocceival
    end if


    call mem_dealloc(Occ_list)
    call mem_dealloc(virt_list)
    call mem_dealloc(occAOSidx)
    call mem_dealloc(virtAOSidx)

  end subroutine fragment_read_data




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
    type(decorbital), dimension(dim1), intent(in) :: Orbitals1
    !> Orbital information for orbitals describing the second index of orbital_mat
    type(decorbital), dimension(dim2), intent(in) :: Orbitals2
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
       write(DECinfo%output,'(1X,a,g18.8)') &
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
    type(decfrag), intent(inout) :: MyFragment
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
    if(MyFragment%ccModel/=MODEL_MP2) then
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
    A = MyFragment%nbasis ! Number of atomic orbitals
    ! Maximum batch dimension in integral subroutines
    Bint = max_batch_dimension(MyFragment%mylsitem,MyFragment%nbasis)
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
            & MyFragment%nbasis
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
    type(decfrag),intent(inout) :: MyFragment
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
    type(decfrag),intent(inout) :: PairFragment
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
    type(decorbital),dimension(nocc) :: OccOrbitals
    !> Virtual orbitals in DEC format
    type(decorbital),dimension(nvirt) :: VirtOrbitals
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
          ! fragments (to avoid double counting).
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
    type(decfrag), intent(inout) :: MyFragment
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
    type(decfrag), intent(inout) :: MyFragment
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
       write(DECinfo%output,*) '--> We redo the fragment calculations'
       print *, 'Relative singles difference is larger than threhold!'
       print *, '--> We redo the fragment calculations'
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
    type(decorbital), dimension(nOcc), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(decorbital), dimension(nUnocc), intent(in) :: UnoccOrbitals
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
       atomtype = MyLsitem%input%molecule%atom(i)%atomic_number ! atom type for atom "i"
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
    type(decorbital), intent(in) :: OccOrbitals(nocc)
    !> Fragment where core information is to be set
    type(decfrag),intent(inout) :: MyFragment
    integer :: i,idx,atom
    logical,pointer :: which_atoms(:), which_core_orbitals(:)


    ! Sanity check
    if( MyFragment%noccAOS==0 ) then
       write(DECinfo%output,'(1X,a,i7)') 'nocc AOS', MyFragment%noccAOS
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
    type(decfrag),intent(inout) :: MyFragment
    !> Density = Cocc Cocc^T
    real(realk),intent(inout),dimension(MyFragment%nbasis,MyFragment%nbasis) :: dens

    call get_density_from_occ_orbitals(MyFragment%nbasis,MyFragment%noccAOS,MyFragment%Co,dens)

  end subroutine get_density_for_fragment

  !> Get core density matrix from fragment as Ccore Ccore^T
  !> where Ccore are the core MO coefficients for the fragment
  !> \author Kasper Kristensen
  !> \date November 2012
  subroutine get_coredensity_for_fragment(MyFragment,dens)
    implicit none
    !> Fragment info
    type(decfrag),intent(inout) :: MyFragment
    !> Density = Cocc Cocc^T
    real(realk),intent(inout),dimension(MyFragment%nbasis,MyFragment%nbasis) :: dens

    call get_density_from_occ_orbitals(MyFragment%nbasis,MyFragment%ncore,MyFragment%coreMO,dens)

  end subroutine get_coredensity_for_fragment



  !> \brief Create job list for DEC fragment optimization calculations.
  !> \author Kasper Kristensen
  !> \date November 2013
  subroutine create_dec_joblist_fragopt(natoms,nocc,nunocc,DistanceTable,&
       & OccOrbitals, UnoccOrbitals, dofrag, mylsitem,jobs)
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
    type(decorbital), dimension(nOcc), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(decorbital), dimension(nUnocc), intent(in) :: UnoccOrbitals
    !> Logical vector telling which atoms have orbitals assigned
    logical,dimension(natoms),intent(in) :: dofrag
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> Job list of fragments listed according to size (all pointers in the type are initialized here)
    type(joblist),intent(inout) :: jobs
    integer,dimension(natoms) :: af_list
    integer :: njobs,i

    ! Get list of atomic fragment ordered according to their (very roughly) estimated sizes
    call estimate_atomic_fragment_sizes(natoms,nocc,nunocc,DistanceTable,&
         & OccOrbitals, UnoccOrbitals, mylsitem,af_list)

    ! Init job list with number of jobs corresponding to number of atoms with orbitals assigned
    njobs = count(dofrag)
    call init_joblist(njobs,jobs)


    ! Set members of job list
    ! ***********************

    do i=1,njobs

       ! Atoms according to list based on estimated sizes
       jobs%atom1(i) = af_list(i)
       jobs%atom2(i) = af_list(i)  ! atomic fragment, so atom1=atom2

       ! Job size is not known because fragment size is not known - simply set it to 1 for all atoms.
       jobs%jobsize(i) = 1

       ! Do fragment optimization!
       jobs%dofragopt(i) = .true.
    end do

    ! All other members of job list were set appropriately by the call to init_joblist.

  end subroutine create_dec_joblist_fragopt


  !> \brief Create job list for DEC calculations remaining after fragment optimization.
  !> \author Kasper Kristensen
  !> \date January 2013
  subroutine create_dec_joblist_driver(calcAF,MyMolecule,mylsitem,natoms,nocc,nunocc,&
       &OccOrbitals,UnoccOrbitals,AtomicFragments,which_fragments,esti,jobs)

    implicit none
    !> Calculate atomic fragments (true) or just pair fragments (false)
    logical,intent(in) :: calcAF
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
    !> Occupied orbitals
    type(decorbital), intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals
    type(decorbital), intent(in) :: UnoccOrbitals(nunocc)
    !> Optimized atomic fragments
    type(decfrag),dimension(natoms),intent(in) :: AtomicFragments
    !> which_fragments(i) is true if atom "i" is central in one of the fragments
    logical, dimension(natoms),intent(in) :: which_fragments
    !> Use estimated fragments?
    logical,intent(in) :: esti
    !> Job list of fragments listed according to size
    type(joblist),intent(inout) :: jobs
    integer :: maxocc,maxunocc,occdim,unoccdim,basisdim,nfrags, minocc,minunocc,minbasis
    integer:: maxbasis, nbasis,atom,idx,i,j,myatom,nsingle,npair,njobs
    real(realk) :: avocc,avunocc,tcpu,twall,avbasis,avRmaxAOS, avRmaxAE, maxRmaxAOS, maxRmaxAE, minRmaxAOS, minRmaxAE
    logical,pointer :: occAOS(:,:),unoccAOS(:,:),fragbasis(:,:)
    integer,pointer :: fragsize(:),fragtrack(:),occsize(:),unoccsize(:),basissize(:)
    real(realk), pointer :: RmaxAOS(:), RmaxAE(:)

    call LSTIMER('START',tcpu,twall,DECinfo%output)

    write(DECinfo%output,*) 'Preparing job list...'
    write(DECinfo%output,*)
    nbasis = MyMolecule%nbasis


    ! Fragment dimension statistics
    ! *****************************
    nsingle    = count(which_fragments)
    maxocc     = 0
    maxunocc   = 0
    maxbasis   = 0
    avocc      = 0.0_realk
    avunocc    = 0.0_realk
    avbasis    = 0.0_realk
    minocc     = huge(minocc)
    minunocc   = huge(minunocc)
    minbasis   = huge(minbasis)
    maxRmaxAOS = 0.0E0_realk
    maxRmaxAE  = 0.0E0_realk
    avRmaxAOS  = 0.0E0_realk
    avRmaxAE   = 0.0E0_realk
    minRmaxAOS = huge(minRmaxAOS)
    minRmaxAE  = huge(minRmaxAE)

    call mem_alloc(occAOS,nocc,natoms)
    call mem_alloc(unoccAOS,nunocc,natoms)
    call mem_alloc(Fragbasis,nbasis,natoms)
    call mem_alloc(fragsize,natoms)
    call mem_alloc(occsize,natoms)
    call mem_alloc(unoccsize,natoms)
    call mem_alloc(basissize,natoms)
    call mem_alloc(RmaxAOS,natoms)
    call mem_alloc(RmaxAE,natoms)
    occAOS    = .false.
    unoccAOS  = .false.
    fragbasis = .false.
    fragsize  = 0
    occsize   = 0
    unoccsize = 0
    basissize = 0


    GetStandardFrag: do atom=1,natoms

       if(.not. which_fragments(atom)) cycle


       ! Set occupied AOS logical vector
       ! ===============================
       do j=1,AtomicFragments(atom)%noccAOS
          idx=AtomicFragments(atom)%occAOSidx(j)  ! index for local occupied AOS orbital
          occAOS(idx,atom) = .true.  ! idx is included in "atom" fragment
       end do


       ! Set unoccupied AOS logical vector
       ! =================================
       do j=1,AtomicFragments(atom)%nunoccAOS
          idx=AtomicFragments(atom)%unoccAOSidx(j)  ! index for unoccupied AOS orbital
          unoccAOS(idx,atom) = .true.  ! idx is included in "atom" fragment
       end do


       ! Set AO basis logical vector
       ! ===========================
       do j=1,AtomicFragments(atom)%nbasis
          idx=AtomicFragments(atom)%basis_idx(j)  ! index for AO basis function
          fragbasis(idx,atom) = .true.  ! idx is included in "atom" fragment
       end do


       ! Statistics: Fragment sizes
       ! ==========================
       if(DECinfo%fragadapt) then ! use dimensions for fragment-adapted orbitals
          occdim   = AtomicFragments(atom)%noccFA
          unoccdim = AtomicFragments(atom)%nunoccFA
       else ! use dimensions for local orbitals
          occdim   = AtomicFragments(atom)%noccAOS
          unoccdim = AtomicFragments(atom)%nunoccAOS
       end if
       basisdim = AtomicFragments(atom)%nbasis

       ! Max and average dimensions
       maxocc   = max(maxocc,occdim)
       maxunocc = max(maxunocc,unoccdim)
       maxbasis = max(maxbasis,basisdim)
       avOCC    = avOCC   + real(occdim)
       avUNOCC  = avUNOCC + real(unoccdim)
       avbasis  = avbasis + real(basisdim)
       minocc   = min(minocc,occdim)
       minunocc = min(minunocc,unoccdim)
       minbasis = min(minbasis,basisdim)
       !Get max distances in Fragment
       RmaxAE(atom)  = AtomicFragments(atom)%RmaxAE
       RmaxAOS(atom) = AtomicFragments(atom)%RmaxAOS
       maxRmaxAE     = max(maxRmaxAE,RmaxAE(atom))
       maxRmaxAOS    = max(maxRmaxAOS,RmaxAOS(atom))
       avRmaxAE      = avRmaxAE  + RmaxAE(atom)
       avRmaxAOS     = avRmaxAOS + RmaxAOS(atom)
       minRmaxAE     = min(minRmaxAE,RmaxAE(atom))
       minRmaxAOS    = min(minRmaxAOS,RmaxAOS(atom))

       ! Store dimensions
       occsize(atom)   = occdim
       unoccsize(atom) = unoccdim
       basissize(atom) = basisdim

       ! Fragment size measure: occ*unocc*basis
       fragsize(atom) = occdim*unoccdim*basisdim


    end do GetStandardFrag


    ! Average dimensions
    avOCC     = avOCC     / real(nsingle)
    avUNOCC   = avUNOCC   / real(nsingle)
    avbasis   = avbasis   / real(nsingle)
    avRmaxAOS = avRmaxAOS / real(nsingle)
    avRmaxAE  = avRmaxAE  / real(nsingle)



    ! Sort standard fragments according to size and print out
    ! *******************************************************
    call mem_alloc(fragtrack,natoms)
    call integer_inv_sort_with_tracking(fragsize, fragtrack, natoms)

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)') '***************************************************************&
         &****************'
    write(DECinfo%output,*) '             Atomic fragments listed according to total size'
    write(DECinfo%output,'(1X,a)') '***************************************************************&
         &****************'

    write(DECinfo%output,*) '   Index     Occupied    Virtual    #Basis funcs.  Rmax(AOS/AE) Rave(AOS/AE) Rsdv(AOS/AE)'
    write(DECinfo%output,*) '             (no. orb)  (no. orb)  '

    do i=1,natoms
       myatom = fragtrack(i)

       PrintFragInfo: if(which_fragments(myatom)) then

          write(DECinfo%output,'(1X,i6,6X,i6,5X,i6,8X,i5,7X,g9.4,"/",g9.4,5X,g9.4,"/",g9.4,5X,g9.4,"/",g9.4)') &
               & myatom, &
               & occsize(myatom), &
               & unoccsize(myatom), &
               & basissize(myatom), &
               & RmaxAOS(myatom)*bohr_to_angstrom, &
               & RmaxAE(myatom)*bohr_to_angstrom

       end if PrintFragInfo

    end do
    write(DECinfo%output,*)

    write(DECinfo%output,'(1X,a,i8,7X,"/",g15.5,"/",i8)')&
    &'FRAGANALYSIS: Max/Ave/Min occ     : ', maxocc,avocc,minocc
    write(DECinfo%output,'(1X,a,i8,7X,"/",g15.5,"/",i8)')&
    &'FRAGANALYSIS: Max/Ave/Min unocc   : ', maxunocc,avunocc,minunocc
    write(DECinfo%output,'(1X,a,i8,7X,"/",g15.5,"/",i8)')&
    &'FRAGANALYSIS: Max/Ave/Min basis   : ', maxbasis,avbasis,minbasis
    write(DECinfo%output,'(1X,a,g15.5,"/",g15.5,"/",g15.5)')&
    &'FRAGANALYSIS: Max/Ave/Min dist AE : ', &
    &maxRmaxAE*bohr_to_angstrom,avRmaxAE*bohr_to_angstrom,minRmaxAE*bohr_to_angstrom
    write(DECinfo%output,'(1X,a,g15.5,"/",g15.5,"/",g15.5)')&
    &'FRAGANALYSIS: Max/Ave/Min dist AOS: ',&
    &maxRmaxAOS*bohr_to_angstrom,avRmaxAOS*bohr_to_angstrom,minRmaxAOS*bohr_to_angstrom
    write(DECinfo%output,*)


    ! Number of pair fragments
    npair=0
    do i=1,natoms
       if(.not. which_fragments(i)) cycle
       do j=i+1,natoms
          if(.not. which_fragments(j)) cycle
          CheckPair: if(MyMolecule%ccmodel(i,j) /= MODEL_NONE .and. &
               & MyMolecule%distancetable(i,j)<DECinfo%pair_distance_threshold ) then
             ! Pair needs to be computed
             npair = npair+1
          end if CheckPair
       end do
    end do


    ! Set job list for fragments
    ! --------------------------
    if(calcAF) then
       ! Calculate atomic fragments: Both atomic and pair frags
       njobs = nsingle+npair
    else
       ! Atomic fragments are already done - only do pairs
       njobs = npair
    end if

    call init_joblist(njobs,jobs)

    call set_dec_joblist(MyMolecule,calcAF,natoms,nocc,nunocc,nbasis,occAOS,unoccAOS,&
         & FragBasis,which_fragments, mymolecule%DistanceTable, esti,jobs)

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '*****************************************************'
    if(esti) then
       write(DECinfo%output,*) '*      DEC ESTIMATED FRAGMENT JOB LIST              *'
    else
       write(DECinfo%output,*) '*               DEC FRAGMENT JOB LIST               *'
    end if
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
    write(DECinfo%output,'(1X,a,i10)') 'DEC JOB SUMMARY: Number of single jobs = ', nsingle
    write(DECinfo%output,'(1X,a,i10)') 'DEC JOB SUMMARY: Number of pair jobs   = ', npair
    write(DECinfo%output,'(1X,a,i10)') 'DEC JOB SUMMARY: Total number of jobs  = ', njobs
    write(DECinfo%output,*)
    write(DECinfo%output,*)

    call mem_dealloc( occAOS    )
    call mem_dealloc( unoccAOS  )
    call mem_dealloc( Fragbasis )
    call mem_dealloc( fragsize  )
    call mem_dealloc( fragtrack )
    call mem_dealloc( occsize   )
    call mem_dealloc( unoccsize )
    call mem_dealloc( basissize )
    call mem_dealloc( RmaxAOS   )
    call mem_dealloc( RmaxAE    )


  end subroutine create_dec_joblist_driver



  !> \brief Construct a job list, joblist12, which contains the jobs of two existing job lists,
  !> joblist1 and joblist2 in the order [joblist1 joblist2], with the internal job orderings within
  !> joblist1 or joblist2 unchanged.
  !> \author Kasper Kristensen
  !> \date November 2013
  subroutine concatenate_joblists(joblist1,joblist2,joblist12)
    implicit none
    !> Joblists 1 and 2
    type(joblist),intent(in) :: joblist1, joblist2
    !> Output joblist containing [joblist1 joblist2]
    type(joblist),intent(inout) :: joblist12
    integer :: njobs,startidx

    ! Number of jobs in combined job list: Sum of original job lists
    njobs = joblist1%njobs + joblist2%njobs
    call init_joblist(njobs,joblist12)


    ! Copy information
    ! ****************

    ! Joblist 1 --> Joblist12
    joblist12%atom1(1:joblist1%njobs) = joblist1%atom1
    joblist12%atom2(1:joblist1%njobs) = joblist1%atom2
    joblist12%jobsize(1:joblist1%njobs) = joblist1%jobsize
    joblist12%jobsdone(1:joblist1%njobs) = joblist1%jobsdone
    joblist12%dofragopt(1:joblist1%njobs) = joblist1%dofragopt
    joblist12%esti(1:joblist1%njobs) = joblist1%esti
    joblist12%nslaves(1:joblist1%njobs) = joblist1%nslaves
    joblist12%nocc(1:joblist1%njobs) = joblist1%nocc
    joblist12%nunocc(1:joblist1%njobs) = joblist1%nunocc
    joblist12%nbasis(1:joblist1%njobs) = joblist1%nbasis
    joblist12%ntasks(1:joblist1%njobs) = joblist1%ntasks
    joblist12%flops(1:joblist1%njobs) = joblist1%flops
    joblist12%LMtime(1:joblist1%njobs) = joblist1%LMtime
    joblist12%load(1:joblist1%njobs) = joblist1%load

    ! Joblist 2 --> Joblist12
    startidx=joblist1%njobs+1
    joblist12%atom1(startidx:njobs) = joblist2%atom1
    joblist12%atom2(startidx:njobs) = joblist2%atom2
    joblist12%jobsize(startidx:njobs) = joblist2%jobsize
    joblist12%jobsdone(startidx:njobs) = joblist2%jobsdone
    joblist12%dofragopt(startidx:njobs) = joblist2%dofragopt
    joblist12%esti(startidx:njobs) = joblist2%esti
    joblist12%nslaves(startidx:njobs) = joblist2%nslaves
    joblist12%nocc(startidx:njobs) = joblist2%nocc
    joblist12%nunocc(startidx:njobs) = joblist2%nunocc
    joblist12%nbasis(startidx:njobs) = joblist2%nbasis
    joblist12%ntasks(startidx:njobs) = joblist2%ntasks
    joblist12%flops(startidx:njobs) = joblist2%flops
    joblist12%LMtime(startidx:njobs) = joblist2%LMtime
    joblist12%load(startidx:njobs) = joblist2%load

  end subroutine concatenate_joblists


  ! \brief Get main info about size for pair fragment from sets of logical arrays
  !> without constructing pair fragment explicitly.
  !> \author Kasper Kristensen
  !> \date March 2012
  subroutine get_main_pair_info(nocc,nunocc,occAOS1,unoccAOS1,occEOS1,unoccEOS1,&
       & occAOS2,unoccAOS2,occEOS2,unoccEOS2,noccAOS,nunoccAOS,noccEOS,nunoccEOS)

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


  end subroutine get_main_pair_info




  !> \brief Set fragment job list. The jobs are listed according to size
  !> with the largest jobs first.
  !> Note: MPI fragment statistics is not modified here.
  !> \author Kasper Kristensen
  !> \date April 2013
  subroutine set_dec_joblist(MyMolecule,calcAF,natoms,nocc,nunocc,nbasis,occAOS,unoccAOS,&
       & FragBasis,which_fragments, DistanceTable, esti,jobs)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Calculate atomic fragments (true) or just pair fragments (false)
    logical,intent(in) :: calcAF
    !> Number of atoms in full molecule
    integer,intent(in) :: natoms
    !> Number of occupied orbitals in full molecule
    integer,intent(in) :: nocc
    !> Number of unoccupied orbitals in full molecule
    integer,intent(in) :: nunocc
    !> Number of basis functions in full molecule
    integer,intent(in) :: nbasis
    !> Logical vector describing occupied AOS (see create_dec_joblist_driver)
    logical,dimension(nocc,natoms),intent(in) :: occAOS
    !> Logical vector describing unoccupied AOS (see create_dec_joblist_driver)
    logical,dimension(nunocc,natoms),intent(in) :: unoccAOS
    !> Logical vector describing which atomic basis functions to include for each fragment
    logical,dimension(nbasis,natoms) :: FragBasis
    !> Logical vector describing which atoms have orbitals assigned
    logical,dimension(natoms),intent(in) :: which_fragments
    !> Distance table with interatomic distances
    real(realk),dimension(natoms,natoms),intent(in) :: DistanceTable
    !> Use estimated fragments?
    logical,intent(in) :: esti
    !> Job list for fragments
    type(joblist),intent(inout) :: jobs
    logical,pointer :: occpairAOS(:), unoccpairAOS(:),basispair(:)
    integer :: i,j,k,njobs,nsingle
    real(realk) :: dist
    integer,pointer :: atom1(:),atom2(:),order(:)

    ! Init stuff
    k=0
    njobs = jobs%njobs
    if(njobs < 1) then ! Sanity check
       call lsquit('set_dec_joblist : Number of jobs must be positive!',-1)
    end if
    call mem_alloc(atom1,njobs)
    call mem_alloc(atom2,njobs)
    atom1=0
    atom2=0
    nsingle=count(which_fragments)
    call mem_alloc(occpairAOS,nocc)
    call mem_alloc(unoccpairAOS,nunocc)
    call mem_alloc(basispair,nbasis)



    ! **************************
    ! MAIN LOOP TO GET JOB SIZES
    ! **************************

    do i=1,natoms  ! Loop over atoms

       if(.not. which_fragments(i)) cycle  ! No fragment for atom i

       ! Repeat atomic fragments if requested
       if(calcAF) then

          ! Set atom indices for single fragment job
          k=k+1
          atom1(k) = i
          atom2(k) = i   ! Same index for both atoms to distinguish single from pair jobs

          ! Job size is defined as occupied AOS * unoccupied AOS dimensions * nbasis
          jobs%jobsize(k) = count(occAOS(1:nocc,i))*count(unoccAOS(1:nunocc,i))&
               &*count(fragbasis(1:nbasis,i))

       end if

       ! Pair loop
       do j=i+1,natoms
          if(.not. which_fragments(j)) cycle ! No fragment for atom j

          ! Distance between atoms i and j
          dist = DistanceTable(i,j)

          CheckPair: if(MyMolecule%ccmodel(i,j) /= MODEL_NONE .and. &
               & MyMolecule%distancetable(i,j)<DECinfo%pair_distance_threshold ) then
             k=k+1

             ! Merge AOS for fragment 1 and 2 for standard pair
             call get_logical_pair_vector(nocc,occAOS(1:nocc,i),occAOS(1:nocc,j),occpairAOS)
             call get_logical_pair_vector(nunocc,unoccAOS(1:nunocc,i),&
                  &unoccAOS(1:nunocc,j),unoccpairAOS)

             ! Logical vector for basis functions
             call get_logical_pair_vector(nbasis,FragBasis(1:nbasis,i),FragBasis(1:nbasis,j),&
                  & basispair)

             ! Atomic indices and jobsize for pair
             atom1(k) = i
             atom2(k) = j

             ! Job size is defined as occupied AOS * unoccupied AOS dimensions * nbasis
             jobs%jobsize(k) = count(occpairAOS)*count(unoccpairAOS)*count(basispair)

             if(jobs%jobsize(k)<1) then
                print *, 'dist', dist
                print *, 'k=',k
                print *, 'pair=', i,j
                print *, 'job size: ', jobs%jobsize(k)
                call lsquit('set_dec_joblist: Non-positive job size, something wrong',DECinfo%output)
             end if

          end if CheckPair

       end do
    end do

    ! Sanity check: k must equal number of jobs
    if(k/=njobs) then
       print *, 'k     = ',k
       print *, 'njobs = ', njobs
       call lsquit('set_dec_joblist: Something wrong in fragment bookkeeping',-1)
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

    ! Estimated fragments?
    jobs%esti=esti

    ! Clean up
    call mem_dealloc(atom1)
    call mem_dealloc(atom2)
    call mem_dealloc(order)
    call mem_dealloc(occpairAOS)
    call mem_dealloc(unoccpairAOS)
    call mem_dealloc(basispair)

  end subroutine set_dec_joblist



  !> \brief Write fragment energies for fragment for easy restart to
  !> file fragenergies.info (or estimate_fragenergies.info for estimated energies).
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine write_fragment_energies_for_restart(natoms,FragEnergies,jobs,esti)

    implicit none
    !> Number of atoms in the molecule
    integer,intent(in) :: natoms
    !> Fragment energies (see decfrag type def)
    real(realk),dimension(natoms,natoms,ndecenergies),intent(in) :: FragEnergies
    !> Job list of fragment jobs
    type(joblist),intent(in) :: jobs
    !> Read estimated fragment energies for restart?
    logical,intent(in) :: esti
    character(len=40) :: FileName, filebackup
    integer :: funit
    logical :: file_exist

    ! Init stuff
    funit = -1
    filename = get_fragenergy_restart_filename(esti)

    ! backup exisiting file
    inquire(file=FileName,exist=file_exist)
    if(file_exist) then
       filebackup = get_fragenergy_restart_filename_backup(esti)
       call rename(filename,filebackup)
    end if

    ! Create a new file fragenergies.info
    call lsopen(funit,FileName,'NEW','UNFORMATTED')

    ! Write job list info and fragment energies
    call write_fragment_joblist_to_file(jobs,funit)
    write(funit) FragEnergies

    call lsclose(funit,'KEEP')

  end subroutine write_fragment_energies_for_restart



  !> \brief Read fragment energies for fragment from file fragenergies.info.
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine read_fragment_energies_for_restart(natoms,FragEnergies,jobs,esti)

    implicit none
    !> Number of atoms in the molecule
    integer,intent(in) :: natoms
    !> Fragment energies (see decfrag type def)
    real(realk),dimension(natoms,natoms,ndecenergies),intent(inout) :: FragEnergies
    !> Job list of fragments
    type(joblist),intent(inout) :: jobs
    !> Read estimated fragment energies for restart?
    logical,intent(in) :: esti
    character(len=40) :: FileName
    integer :: funit
    logical :: file_exist

    ! Init stuff
    funit = -1
    filename = get_fragenergy_restart_filename(esti)

    ! Sanity check
    inquire(file=FileName,exist=file_exist)
    if(.not. file_exist) then
       write(DECinfo%output,*) 'WARNING: Restart file: ', FileName
       write(DECinfo%output,*) 'does not exist! I therefore calculate it from scratch...'
       return
    end if

    ! Create a new file fragenergies.info
    call lsopen(funit,FileName,'OLD','UNFORMATTED')

    ! Read job list and fragment energies from file
    call read_fragment_joblist_from_file(jobs,funit)
    read(funit) FragEnergies
    call lsclose(funit,'KEEP')

  end subroutine read_fragment_energies_for_restart



  !> \brief Expand fragment job list to include additional pair fragments
  !> (the additional pairs to include are determined in dec_energy_control_center).
  !> Note: The order of the existing jobs will be unchanged, the new jobs will
  !> simply be appended at the end of the job list.
  !> \author Kasper Kristensen
  !> \data October 2012
  subroutine expand_joblist_to_include_more_pairs(nocc,nunocc,natoms,dofrag,&
       & Fragments,oldpaircut,newpaircut,MyMolecule,OccOrbitals,UnoccOrbitals,jobs)

    implicit none
    !> Number of occupied orbitals in full molecule
    integer,intent(in) :: nocc
    !> Number of unoccupied orbitals in full molecule
    integer,intent(in) :: nunocc
    !> Number of atoms in full molecule
    integer,intent(in) :: nAtoms
    !> dofrag(P) is true if P is central atom for a fragment (see main_fragment_driver)
    logical,intent(in) :: dofrag(natoms)
    !> Single fragments
    type(decfrag),intent(inout) :: Fragments(natoms)
    !> Pair cutoff distance used for current job list
    real(realk),intent(in) :: oldpaircut
    !> Pair cutoff distance to be used for new job list
    real(realk),intent(in) :: newpaircut
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Occupied orbitals
    type(decorbital), intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals
    type(decorbital), intent(in) :: UnoccOrbitals(nunocc)
    !> Current job list which will be appended with new pairs
    type(joblist),intent(inout) :: jobs
    type(joblist) :: oldjobs
    integer :: i,j,nsingle,npairold,npairnew,npairdelta, nold,nnew,k,nbasisFragment,n
    integer,pointer :: atom1(:), atom2(:), jobsize(:),order(:)
    logical,pointer :: occAOS(:,:), unoccAOS(:,:)
    logical,pointer :: occpairAOS(:), unoccpairAOS(:)
    real(realk) :: dist

    ! Counting: Number of single fragments, old pairs, and new pairs
    ! **************************************************************
    nsingle=0
    npairold=0
    npairnew=0
    do i=1,natoms
       if(.not. dofrag(i)) cycle ! atom "i" is not central atom in a  fragment
       nsingle = nsingle +1 ! number of single  fragments

       do j=i+1,natoms
          if(.not. dofrag(j)) cycle ! atom "j" is not central atom in a  fragment

          ! Number of pairs for old cutoff
          if(mymolecule%DistanceTable(i,j) < oldpaircut) npairold = npairold+1

          ! Number of pairs for new cutoff
          if(mymolecule%DistanceTable(i,j) < newpaircut) npairnew = npairnew+1

       end do
    end do

    ! Number of jobs in job list (for MP2 energy calculations atomic frags are not
    ! included in job list because they have already been determined during fragment optimization)
    if(DECinfo%ccmodel==MODEL_MP2 .and. .not. DECinfo%first_order) then
       n = npairold
    else
       n = nsingle+npairold
    end if


    ! Sanity check 1: Number current jobs should be sum of single jobs + old pairs
    if( jobs%njobs /= n ) then
       write(DECinfo%output,*) 'Number of single  frags: ', nsingle
       write(DECinfo%output,*) 'Number of pair    frags: ', npairold
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
    jobs%dofragopt(1:nold) = oldjobs%dofragopt(1:nold)
    jobs%esti(1:nold) = oldjobs%esti(1:nold)
    jobs%nslaves(1:nold) = oldjobs%nslaves(1:nold)
    jobs%nocc(1:nold) = oldjobs%nocc(1:nold)
    jobs%nunocc(1:nold) = oldjobs%nunocc(1:nold)
    jobs%nbasis(1:nold) = oldjobs%nbasis(1:nold)
    jobs%ntasks(1:nold) = oldjobs%ntasks(1:nold)
    jobs%flops(1:nold) = oldjobs%flops(1:nold)
    jobs%LMtime(1:nold) = oldjobs%LMtime(1:nold)
    jobs%load(1:nold) = oldjobs%load(1:nold)


    ! Set logical vectors with fragment orbital information
    ! *****************************************************
    ! occAOS(i,j) is true (false) if occupied orbital "i" is (not) included in  fragment "j"
    call mem_alloc(occAOS,nocc,natoms)
    call mem_alloc(unoccAOS,nunocc,natoms)
    ! same for reduced fragment spaces (see decfrag type)
    call get_logical_vectors_for_AOS(nocc,nunocc,natoms,dofrag,Fragments,&
         & occAOS, unoccAOS)



    ! Set new job list
    ! ****************

    k = 0
    call mem_alloc(occpairAOS,nocc)
    call mem_alloc(unoccpairAOS,nunocc)
    call mem_alloc(atom1,npairdelta)
    call mem_alloc(atom2,npairdelta)
    call mem_alloc(jobsize,npairdelta)

    do i=1,natoms  ! Loop over atoms
       if(.not. dofrag(i)) cycle  ! No fragment for atom i

       do j=i+1,natoms
          if(.not. dofrag(j) ) cycle ! No  fragment for atom j

          ! Distance between atoms i and j
          dist = mymolecule%DistanceTable(i,j)

          ! Add pairs with distances between old and new paircut
          AddPairToJobList: if( (dist < newpaircut) .and. (dist >= oldpaircut) ) then

             ! Merge AOS for fragment 1 and 2 for standard pair
             call get_logical_pair_vector(nocc,occAOS(1:nocc,i),occAOS(1:nocc,j),occpairAOS)
             call get_logical_pair_vector(nunocc,unoccAOS(1:nunocc,i),&
                  &unoccAOS(1:nunocc,j),unoccpairAOS)

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

  end subroutine expand_joblist_to_include_more_pairs


  !> Determine logical vector for occ and unocc AOS where e.g.
  !> occAOS(i,P) is true (false) if occupied orbital "i" is (not) included in  fragment "P".
  !> The same is done for fragment spaces of reduced size.
  !> \author Kasper Kristensen
  !> \date October 2012
  subroutine get_logical_vectors_for_AOS(nocc,nunocc,natoms,dofrag,Fragments,&
       & occAOS, unoccAOS)

    implicit none
    !> Number of occupied orbitals in full molecule
    integer,intent(in) :: nocc
    !> Number of unoccupied orbitals in full molecule
    integer,intent(in) :: nunocc
    !> Number of atoms in full molecule
    integer,intent(in) :: nAtoms
    !> dofrag(P) is true if P is central atom for a  fragment (see main_fragment_driver)
    logical,intent(in) :: dofrag(natoms)
    !> Single  fragments
    type(decfrag),intent(inout) :: Fragments(natoms)
    !> Logical vector for occupied AOS
    logical,intent(inout) :: occAOS(nocc,natoms)
    !> Logical vector for unoccupied AOS
    logical,intent(inout) :: unoccAOS(nunocc,natoms)
    integer :: i,idx,P

    ! Init
    occAOS = .false.
    unoccAOS = .false.


    ! Set logical vectors
    ! *******************

    do P=1,natoms  ! loop over atoms
       if(.not. dofrag(P)) cycle  ! skip if P is not a  fragment

       ! Set occupied AOS for P
       do i=1,Fragments(P)%noccAOS
          ! index "idx" in total occupied orbital list for orbital "i" in fragment orbital list
          idx = Fragments(P)%occAOSidx(i)
          occAOS(idx,P) = .true.
       end do

       ! Set unoccupied AOS for P
       do i=1,Fragments(P)%nunoccAOS
          idx = Fragments(P)%unoccAOSidx(i)
          unoccAOS(idx,P) = .true.
       end do

    end do


  end subroutine get_logical_vectors_for_AOS


  !> Copy  fragment job list (the new copied joblist is also initialized here)
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
    jobscopy%dofragopt = jobs%dofragopt
    jobscopy%esti = jobs%esti
    jobscopy%nslaves = jobs%nslaves
    jobscopy%nocc = jobs%nocc
    jobscopy%nunocc = jobs%nunocc
    jobscopy%nbasis = jobs%nbasis
    jobscopy%ntasks = jobs%ntasks
    jobscopy%flops = jobs%flops
    jobscopy%LMtime = jobs%LMtime
    jobscopy%load = jobs%load

  end subroutine copy_joblist

  !> \brief Compare list of single  fragments
  !> for two different  fragment job lists.
  !> If they do not match, the program prints error message and quits.
  !> (deliberately we do not compare pairs since they may be different if
  !>  the pair cutoff distance is increased).
  subroutine fragment_sanity_check(jobs1,jobs2)
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
       call lsquit('fragment_sanity_check: The number of single fragments differ!',-1)
    end if

    ! Same max index
    if(max1/=max2) then
       write(DECinfo%output,*) 'max1,max2: ',max1,max2
       call lsquit('fragment_sanity_check: Maximum atomic indices differ!',-1)
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
          call lsquit('fragment_sanity_check: Single  fragment indices differ!',-1)
       end if
    end do


    call mem_dealloc(singlejobs1)
    call mem_dealloc(singlejobs2)

  end subroutine fragment_sanity_check


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
    type(decorbital), intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals
    type(decorbital), intent(in) :: UnoccOrbitals(nunocc)
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


  !> \brief Set pair fragment-adapted MO coefficients by
  !> taking the union of fragment-adapted orbitals for the incoming atomic fragments,
  !> but where only FOs in the original atomic fragments with eigenvalues larger than some
  !> threshold are kept (see get_pairFO_union).
  !> Note 1: These orbitals are NOT orthogonal, and they
  !> are also redundant! Redundancies need to be removed
  !> and orbitals must be orthogonalized before they can be used
  !> in DEC fragment calculation (see pair_fragment_adapted_transformation_matrices)
  !> Note 2: The EOS orbitals (assigned to atoms P and Q) are taken out of the orbital
  !> pool here since they should not be modified.
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine set_pair_fragment_adapted_redundant_orbitals(MyMolecule,FragmentP,FragmentQ,&
       & FragmentPQ,noccPQ,nunoccPQ,WhichOccP,WhichOccQ,WhichUnoccP,WhichUnoccQ,CoccPQ,CunoccPQ)
    implicit none

    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Fragment P
    type(decfrag),intent(inout) :: fragmentP
    !> Fragment Q
    type(decfrag),intent(inout) :: fragmentQ
    !> Pair Fragment PQ
    type(decfrag),intent(inout) :: FragmentPQ
    !> Number of occupied redundant pair orbitals (see get_pairFO_union)
    integer,intent(in) :: noccPQ
    !> Number of unoccupied redundant pair orbitals (see get_pairFO_union)
    integer,intent(in) :: nunoccPQ
    !> Which occupied FOs to include from atomic fragments P and Q (see get_pairFO_union)
    logical,intent(in) :: WhichOccP(fragmentP%noccFA), WhichOccQ(fragmentQ%noccFA)
    !> Which unoccupied FOs to include from atomic fragments P and Q (see get_pairFO_union)
    logical,intent(in) :: WhichUnoccP(fragmentP%nunoccFA), WhichUnoccQ(fragmentQ%nunoccFA)
    !> Occupied MO coefficients
    real(realk),intent(inout) :: CoccPQ(FragmentPQ%nbasis,noccPQ)
    !> Unoccupied MO coefficients
    real(realk),intent(inout) :: CunoccPQ(FragmentPQ%nbasis,nunoccPQ)
    integer :: i,j,ix,jx,orbidx
    integer,pointer :: Pbasis(:), Qbasis(:), PQbasis(:)
    integer,pointer :: fragP_in_pair(:),fragQ_in_pair(:)


    ! Set list of orbitals in full molecular list
    call mem_alloc(Pbasis,FragmentP%nbasis)
    call mem_alloc(Qbasis,FragmentQ%nbasis)
    call mem_alloc(PQbasis,FragmentPQ%nbasis)
    call which_orbitals_in_atomic_extents(MyMolecule,FragmentP,Pbasis)
    call which_orbitals_in_atomic_extents(MyMolecule,FragmentQ,Qbasis)
    call which_orbitals_in_atomic_extents(MyMolecule,FragmentPQ,PQbasis)

    ! Set list of atomic fragment extent in pair fragment extent list
    call mem_alloc(fragP_in_pair,FragmentP%nbasis)
    call mem_alloc(fragQ_in_pair,FragmentQ%nbasis)
    call atomicfragAE_in_pairfragAE(FragmentP%nbasis,FragmentPQ%nbasis,&
         & Pbasis,PQbasis,fragP_in_pair)
    call atomicfragAE_in_pairfragAE(FragmentQ%nbasis,FragmentPQ%nbasis,&
         & Qbasis,PQbasis,fragQ_in_pair)


    ! Set occupied FA-orbitals for pair
    ! *********************************

    ! MO coefficients for pair PQ
    !
    ! C_PQ  =  ( C_P C_Q )     (*)
    !
    ! i.e. with fragment P orbitals before fragment Q orbitals and where orbitals
    ! in the original P and Q atomic fragments are only included if they have eigenvalues
    ! below threshold (see get_pairFO_union)

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
    orbidx=0
    do j=1,fragmentP%noccFA
       if(WhichOccP(j)) then ! "j" is XOS orbital that we want to include
          orbidx = orbidx+1  ! Orbital index for pair MOs
          do i=1,FragmentP%nbasis
             ix = fragP_in_pair(i)
             CoccPQ(ix,orbidx) = FragmentP%CoFA(i,j)
          end do
       end if
    end do

    ! Put FA orbitals for Q into pair fragment orbital matrix
    do j=1,fragmentQ%noccFA
       if(WhichOccQ(j)) then ! "j" is XOS orbital that we want to include
          orbidx = orbidx +1   ! Orbital index for pair MOs
          do i=1,FragmentQ%nbasis
             ix = fragQ_in_pair(i)
             CoccPQ(ix,orbidx) = FragmentQ%CoFA(i,j)
          end do
       end if
    end do



    ! Set unoccupied FA-orbitals for pair
    ! ***********************************
    ! Same strategy as for occupied space.


    ! Put FA occupied orbitals for P into pair fragment orbital matrix
    CunoccPQ=0.0_realk
    orbidx=0
    do j=1,fragmentP%nunoccFA
       if(WhichUnoccP(j)) then ! "j" is XOS orbital that we want to include
          orbidx = orbidx+1  ! Orbital index for pair MOs
          do i=1,FragmentP%nbasis
             ix = fragP_in_pair(i)
             CunoccPQ(ix,orbidx) = FragmentP%CvFA(i,j)
          end do
       end if
    end do

    ! Put FA unoccupied orbitals for Q into pair fragment orbital matrix
    do j=1,fragmentQ%nunoccFA
       if(WhichUnoccQ(j)) then ! "j" is XOS orbital that we want to include
          orbidx = orbidx +1   ! Orbital index for pair MOs
          do i=1,FragmentQ%nbasis
             ix = fragQ_in_pair(i)
             CunoccPQ(ix,orbidx) = FragmentQ%CvFA(i,j)
          end do
       end if
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
    type(decfrag),intent(inout) :: fragment
    !> Which atomic orbitals in full AO list are included in atomic fragment extent?
    integer,intent(inout) :: which_orbitals(fragment%nbasis)
    integer :: i,j,idx,atom

    idx=0
    do i=1,fragment%natoms  ! loop over atoms in fragment
       atom = fragment%atoms_idx(i)  ! atom index in full atom list

       do j=1,MyMolecule%atom_size(atom)  ! loop over basis functions for atom
          idx =idx+1
          which_orbitals(idx) = MyMolecule%atom_start(atom) + j -1
       end do

    end do

    ! Sanity check
    if(idx/=fragment%nbasis) then
       call lsquit('which_orbitals_in_atomic_extents: &
            & counter does not match number of basis functions in atomic &
            & fragment extent!',-1)
    end if


  end subroutine which_orbitals_in_atomic_extents



  !> \brief Extract XOS orbitals from occupied AOS orbitals. The set of XOS orbitals is the
  !> Xtra Orbital Space needed in addition to the EOS orbitals. For example,
  !> for fragment P the XOS orbitals are the subset of AOS orbitals NOT assigned to atom P.
  !> \author Kasper Kristensen
  !> \date April 2013
  subroutine extract_XOS_orbitals_occ(MyFragment,nbasis,nXOS,XOS)
    implicit none
    !> Fragment info (atomic fragment or pair fragment)
    type(decfrag),intent(inout) :: MyFragment
    !> Number of basis functions in fragment
    integer,intent(in) :: nbasis
    !> Number of XOS orbitals
    integer,intent(in) :: nXOS
    !> XOS MO coefficients
    real(realk),intent(inout),dimension(nbasis,nXOS) :: XOS
    integer :: i,idx,nAOS,nEOS
    logical,pointer :: which_XOS(:)

    ! Dimensions
    nAOS = MyFragment%noccAOS
    nEOS = MyFragment%noccEOS


    ! Sanity checks for dimensions
    ! ****************************
    if(nXOS /= nAOS - nEOS) then
       print '(a,3i8)','EOS,AOS,XOS', nEOS,nAOS,nXOS
       call lsquit('extract_XOS_orbitals_occ: XOS dimension mismatch!',-1)
    end if

    ! AO basis dimension consistent
    if(MyFragment%nbasis /= nbasis) then
       print '(a,i8)','basis input', nbasis
       print '(a,i8)','basis frag ', MyFragment%nbasis
       call lsquit('extract_XOS_orbitals_occ: AO basis dimension mismatch!',-1)
    end if

    ! Special case: No XOS orbitals, just return
    if(nXOS==0) return


    ! Which AOS orbitals are XOS orbitals?
    call mem_alloc(which_XOS,nAOS)
    which_XOS=.true.
    do i=1,nEOS
       ! Index MyFragment%idxo(i) is an EOS orbital and therefore not an XOS orbital
       which_XOS(MyFragment%idxo(i)) = .false.
    end do

    ! Extract XOS orbitals and put them into XOS output array
    idx=0
    do i=1,nAOS
       if(which_XOS(i)) then
          idx=idx+1
          XOS(:,idx) = MyFragment%Co(:,i)
       end if
    end do
    if(idx/=nXOS) then
       call lsquit('extract_XOS_orbitals_occ: XOS book keeping error',-1)
    end if
    call mem_dealloc(which_XOS)

  end subroutine extract_XOS_orbitals_occ


  !> \brief Extract XOS orbitals from unoccupied AOS orbitals, see extract_XOS_orbitals_occ.
  !> \author Kasper Kristensen
  !> \date April 2013
  subroutine extract_XOS_orbitals_unocc(MyFragment,nbasis,nXOS,XOS)
    implicit none
    !> Fragment info (atomic fragment or pair fragment)
    type(decfrag),intent(inout) :: MyFragment
    !> Number of basis functions in fragment
    integer,intent(in) :: nbasis
    !> Number of XOS orbitals
    integer,intent(in) :: nXOS
    !> XOS MO coefficients
    real(realk),intent(inout),dimension(nbasis,nXOS) :: XOS
    integer :: i,idx,nAOS,nEOS
    logical,pointer :: which_XOS(:)

    ! Dimensions
    nAOS = MyFragment%nunoccAOS
    nEOS = MyFragment%nunoccEOS


    ! Sanity checks for dimensions
    ! ****************************
    if(nXOS /= nAOS - nEOS) then
       print '(a,3i8)','EOS,AOS,XOS', nEOS,nAOS,nXOS
       call lsquit('extract_XOS_orbitals_unocc: XOS dimension mismatch!',-1)
    end if

    ! AO basis dimension consistent
    if(MyFragment%nbasis /= nbasis) then
       print '(a,i8)','basis input', nbasis
       print '(a,i8)','basis frag ', MyFragment%nbasis
       call lsquit('extract_XOS_orbitals_unocc: AO basis dimension mismatch!',-1)
    end if

    ! Special case: No XOS orbitals, just return
    if(nXOS==0) return


    ! Which AOS orbitals are XOS orbitals?
    call mem_alloc(which_XOS,nAOS)
    which_XOS=.true.
    do i=1,nEOS
       ! Index MyFragment%idxu(i) is an EOS orbital and therefore not an XOS orbital
       which_XOS(MyFragment%idxu(i)) = .false.
    end do

    ! Extract XOS orbitals and put them into XOS output array
    idx=0
    do i=1,nAOS
       if(which_XOS(i)) then
          idx=idx+1
          XOS(:,idx) = MyFragment%Cv(:,i)
       end if
    end do
    if(idx/=nXOS) then
       call lsquit('extract_XOS_orbitals_occ: XOS book keeping error',-1)
    end if
    call mem_dealloc(which_XOS)

  end subroutine extract_XOS_orbitals_unocc


  !> \brief Put occupied XOS orbitals into fragment structure, effectively copying elements
  !> "the inverse way" of what extract_XOS_orbitals_occ is doing.
  !> \author Kasper Kristensen
  !> \date April 2013
  subroutine put_XOS_orbitals_occ(nbasis,nXOS,XOS,MyFragment)
    implicit none
    !> Number of basis functions in fragment
    integer,intent(in) :: nbasis
    !> Number of XOS orbitals
    integer,intent(in) :: nXOS
    !> XOS MO coefficients
    real(realk),intent(in),dimension(nbasis,nXOS) :: XOS
    !> Fragment info (MyFragment%Co will be modified here)
    type(decfrag),intent(inout) :: MyFragment
    integer :: i,idx,nAOS,nEOS
    logical,pointer :: which_XOS(:)


    ! Dimensions
    nAOS = MyFragment%noccAOS
    nEOS = MyFragment%noccEOS

    ! Sanity checks for dimensions
    ! ****************************
    if(nXOS /= nAOS - nEOS) then
       print '(a,3i8)','EOS,AOS,XOS', nEOS,nAOS,nXOS
       call lsquit('put_XOS_orbitals_occ: XOS dimension mismatch!',-1)
    end if

    ! AO basis dimension consistent
    if(MyFragment%nbasis /= nbasis) then
       print '(a,i8)','basis input', nbasis
       print '(a,i8)','basis frag ', MyFragment%nbasis
       call lsquit('put_XOS_orbitals_occ: AO basis dimension mismatch!',-1)
    end if

    ! Special case: No XOS orbitals, just return
    if(nXOS==0) return


    ! Which AOS orbitals are XOS orbitals?
    call mem_alloc(which_XOS,nAOS)
    which_XOS=.true.
    do i=1,nEOS
       ! Index MyFragment%idxo(i) is an EOS orbital and therefore not an XOS orbital
       which_XOS(MyFragment%idxo(i)) = .false.
    end do

    ! Put XOS orbitals into fragment structure
    idx=0
    do i=1,nAOS
       if(which_XOS(i)) then
          idx=idx+1
          MyFragment%Co(:,i) = XOS(:,idx)
       end if
    end do
    if(idx/=nXOS) then
       call lsquit('extract_XOS_orbitals_occ: XOS book keeping error',-1)
    end if
    call mem_dealloc(which_XOS)

  end subroutine put_XOS_orbitals_occ



  !> \brief Put unoccupied XOS orbitals into fragment structure, effectively copying elements
  !> "the inverse way" of what extract_XOS_orbitals_unocc is doing.
  !> \author Kasper Kristensen
  !> \date April 2013
  subroutine put_XOS_orbitals_unocc(nbasis,nXOS,XOS,MyFragment)
    implicit none
    !> Number of basis functions in fragment
    integer,intent(in) :: nbasis
    !> Number of XOS orbitals
    integer,intent(in) :: nXOS
    !> XOS MO coefficients
    real(realk),intent(in),dimension(nbasis,nXOS) :: XOS
    !> Fragment info (MyFragment%Cv will be modified here)
    type(decfrag),intent(inout) :: MyFragment
    integer :: i,idx,nAOS,nEOS
    logical,pointer :: which_XOS(:)

    ! Dimensions
    nAOS = MyFragment%nunoccAOS
    nEOS = MyFragment%nunoccEOS

    ! Sanity checks for dimensions
    ! ****************************
    if(nXOS /= nAOS - nEOS) then
       print '(a,3i8)','EOS,AOS,XOS', nEOS,nAOS,nXOS
       call lsquit('put_XOS_orbitals_unocc: XOS dimension mismatch!',-1)
    end if

    ! AO basis dimension consistent
    if(MyFragment%nbasis /= nbasis) then
       print '(a,i8)','basis input', nbasis
       print '(a,i8)','basis frag ', MyFragment%nbasis
       call lsquit('put_XOS_orbitals_unocc: AO basis dimension mismatch!',-1)
    end if

    ! Special case: No XOS orbitals, just return
    if(nXOS==0) return


    ! Which AOS orbitals are XOS orbitals?
    call mem_alloc(which_XOS,nAOS)
    which_XOS=.true.
    do i=1,nEOS
       ! Index MyFragment%idxu(i) is an EOS orbital and therefore not an XOS orbital
       which_XOS(MyFragment%idxu(i)) = .false.
    end do

    ! Put XOS orbitals into fragment structure
    idx=0
    do i=1,nAOS
       if(which_XOS(i)) then
          idx=idx+1
          MyFragment%Cv(:,i) = XOS(:,idx)
       end if
    end do
    if(idx/=nXOS) then
       call lsquit('extract_XOS_orbitals_unocc: XOS book keeping error',-1)
    end if
    call mem_dealloc(which_XOS)

  end subroutine put_XOS_orbitals_unocc


  !> \brief Copy basic fragment information into job structure
  !> \author Kasper Kristensen
  !> \date May 2013
  subroutine copy_fragment_info_job(myfragment,myjob)
    implicit none
    !> Fragent info
    type(decfrag),intent(in) :: myfragment
    !> Job list of length 1
    type(joblist) :: myjob

    ! Sanity check
    if(myjob%njobs/=1) then
       call lsquit('copy_fragment_info_job: Length of job list is not 1!',-1)
    end if


    ! Copy info
    if(myfragment%nEOSatoms==1) then ! atomic fragment
       myjob%atom1(1) = myfragment%EOSatoms(1)
       myjob%atom2(1) = myfragment%EOSatoms(1)
    else
       myjob%atom1(1) = myfragment%EOSatoms(1)
       myjob%atom2(1) = myfragment%EOSatoms(2)
    end if

    myjob%nocc(1) = myfragment%noccAOS
    myjob%nunocc(1) = myfragment%nunoccAOS
    myjob%nbasis(1) = myfragment%nbasis
    myjob%ntasks(1) = myfragment%ntasks
    myjob%jobsize(1) = myjob%nocc(1)*myjob%nunocc(1)*myjob%nbasis(1)

  end subroutine copy_fragment_info_job


  !> \brief Calculate occ-occ and virt-virt blocks of correlation density matrix
  !> for atomic fragment from AOS amplitude input.
  !> Different definitions for the density matrices are given depending on the model
  !> and on whether we use both occupied and virtual partitioning schemes (details inside subroutine).
  !> \author Kasper Kristensen
  !> \date August 2013
  subroutine calculate_corrdens_frag(t2,MyFragment)
    implicit none
    !> Doubles amplitudes in AOS stored as (a,i,b,j)
    type(array4),intent(in) :: t2
    !> MyFragment - output occ and virt density matrices are stored in
    !> MyFragment%occmat and MyFragment%virtmat, respectively.
    type(decfrag),intent(inout) :: MyFragment

    ! Delete existing correlation density matrix (if present)
    if(MyFragment%CDset) then
       call mem_dealloc(MyFragment%occmat)
       call mem_dealloc(MyFragment%virtmat)
    end if
    call mem_alloc(MyFragment%occmat,MyFragment%noccAOS,MyFragment%noccAOS)
    call mem_alloc(MyFragment%virtmat,MyFragment%nunoccAOS,MyFragment%nunoccAOS)

    ! Different density matrix definitions depending on scheme
    ! - this is work in progress and will probably be modified.
    ! See DECsettings type definition for details.
    print *, 'CORRDENS: Using scheme ', DECinfo%CorrDensScheme

    CorrDensDefinition: select case(DECinfo%CorrDensScheme)
    case(1)
       ! Construct density matrix based only on EOS amplitudes
       call calculate_corrdens_EOS(t2,MyFragment)
    case(2)
       ! Use AOS amplitudes but with special emphasis on EOS for virtual FOs,
       ! while, for occupied FOs, we put equal weight on all AOS amplitudes.
       call calculate_corrdens_semiEOS(t2,MyFragment)
    case(3)
       ! Use AOS amplitudes with equal weight on all amplitudes.
       call calculate_corrdens_AOS(t2,MyFragment)
    case default
       call lsquit('calculate_corrdens: Invalid corrdens scheme',-1)
    end select CorrDensDefinition

    MyFragment%CDset=.true.

  end subroutine calculate_corrdens_frag


  !> \brief Calculate occ-occ and virt-virt blocks of correlation density matrix
  !> for atomic fragment using EOS subset of AOS amplitudes (see inside subroutine).
  !> \author Kasper Kristensen
  !> \date August 2013
  subroutine calculate_corrdens_EOS(t2,MyFragment)
    implicit none
    !> Doubles amplitudes in AOS stored as (a,i,b,j)
    type(array4),intent(in) :: t2
    !> MyFragment - output occ and virt density matrices are stored in
    !> MyFragment%occmat and MyFragment%virtmat, respectively.
    type(decfrag),intent(inout) :: MyFragment
    integer :: i,j,k,a,b,c,ax,bx,ix,jx
    real(realk) :: i2,i4

    i2 = 2.0_realk
    i4 = 4.0_realk

    ! Occ-occ block of density matrix
    ! *******************************
    ! OccMat(i,j) = sum_{abk} t_{ik}^{ab}  tbar_{jk}^{ab}
    ! where ab are assigned to central atom in fragment, and ijk are in occupied AOS
    ! tbar_{ij}^{ab} = 4t_{ij}^{ab} - 2t_{ij}^{ba}
    MyFragment%OccMat = 0.0_realk
    do b=1,MyFragment%nunoccEOS
       bx = MyFragment%idxu(b)  ! index for EOS orbital b in AOS list of orbitals
       do a=1,MyFragment%nunoccEOS
          ax = MyFragment%idxu(a)
          do k=1,MyFragment%noccAOS
             do i=1,MyFragment%noccAOS
                do j=1,MyFragment%noccAOS
                   MyFragment%OccMat(i,j) = MyFragment%OccMat(i,j) + &
                        & t2%val(ax,i,bx,k)*(i4*t2%val(ax,j,bx,k) - i2*t2%val(bx,j,ax,k))
                end do
             end do
          end do
       end do
    end do


    ! Virt-virt block of density matrix
    ! *********************************

    ! VirtMat(a,b) = sum_{ijc} t_{ij}^{ac}  tbar_{ij}^{bc}
    ! where ij are assigned to central atom in fragment, and abc are in virtual AOS
    MyFragment%VirtMat = 0.0_realk
    do i=1,MyFragment%noccEOS
       ix = MyFragment%idxo(i)
       do j=1,MyFragment%noccEOS
          jx = MyFragment%idxo(j)
          do a=1,MyFragment%nunoccAOS
             do b=1,MyFragment%nunoccAOS
                do c=1,MyFragment%nunoccAOS
                   MyFragment%VirtMat(a,b) = MyFragment%VirtMat(a,b) + &
                        & t2%val(a,ix,c,jx)*(i4*t2%val(b,ix,c,jx) - i2*t2%val(c,ix,b,jx))
                end do
             end do
          end do
       end do
    end do

  end subroutine calculate_corrdens_EOS


  !> \brief Calculate occ-occ and virt-virt blocks of correlation density matrix
  !> for atomic fragment using AOS amplitudes - but where the EOS amplitudes are given more weight
  !> (see details inside subroutine).
  !> \author Kasper Kristensen
  !> \date August 2013
  subroutine calculate_corrdens_semiEOS(t2,MyFragment)
    implicit none
    !> Doubles amplitudes in AOS stored as (a,i,b,j)
    type(array4),intent(in) :: t2
    !> MyFragment - output occ and virt density matrices are stored in
    !> MyFragment%occmat and MyFragment%virtmat, respectively.
    type(decfrag),intent(inout) :: MyFragment
    integer :: i,j,a,b,c
    real(realk) :: i2,i4
    logical,pointer :: OccEOS(:)

    i2 = 2.0_realk
    i4 = 4.0_realk

    ! Occ-occ block of density matrix
    ! *******************************
    ! OccMat(i,j) = sum_{abk} t_{ik}^{ab}  tbar_{jk}^{ab}
    ! tbar_{ij}^{ab} = 4t_{ij}^{ab} - 2t_{ij}^{ba}
    ! - for now we let all indices be AOS indices here
    ! (meaning that we do not give the virtual orbitals assigned to P
    ! special weight, and thus that we do not consider the virtual part. scheme).
    call calculate_corrdens_AOS_occocc(t2,MyFragment)



    ! Virt-virt block of density matrix
    ! *********************************

    ! VirtMat(a,b) = sum_{ijc} t_{ij}^{ac}  tbar_{ij}^{bc}
    ! where EITHER i OR j is assigned to central atom in fragment,
    ! and the other occupied index is in the AOS (including the EOS).
    ! The virtual abc indices are in the AOS.
    ! In this way we put emphasis on the EOS amplitudes and the remaining AOS
    ! amplitudes which couple most strongly to the EOS.

    ! Logical vector for occ EOS
    call mem_alloc(OccEOS,MyFragment%noccAOS)
    OccEOS=.false.
    do i=1,MyFragment%noccEOS
       OccEOS(MyFragment%idxo(i)) =.true.
    end do


    MyFragment%VirtMat = 0.0_realk
    do i=1,MyFragment%noccAOS
       do j=1,MyFragment%noccAOS
          ! Only consider contribution if either orbital "i" or "j"
          ! is an EOS orbital (assigned to central atom)
          AddContribution: if(OccEOS(i) .or. OccEOS(j)) then
             do a=1,MyFragment%nunoccAOS
                do b=1,MyFragment%nunoccAOS
                   do c=1,MyFragment%nunoccAOS
                      MyFragment%VirtMat(a,b) = MyFragment%VirtMat(a,b) + &
                           & t2%val(a,i,c,j)*(i4*t2%val(b,i,c,j) - i2*t2%val(c,i,b,j))
                   end do
                end do
             end do
          end if AddContribution
       end do
    end do

    call mem_dealloc(OccEOS)

  end subroutine calculate_corrdens_semiEOS





  !> \brief Calculate occ-occ and virt-virt blocks of correlation density matrix
  !> for atomic fragment using all AOS amplitudes (see details inside subroutine).
  !> \author Kasper Kristensen
  !> \date August 2013
  subroutine calculate_corrdens_AOS(t2,MyFragment)
    implicit none
    !> Doubles amplitudes in AOS stored as (a,i,b,j)
    type(array4),intent(in) :: t2
    !> MyFragment - output occ and virt density matrices are stored in
    !> MyFragment%occmat and MyFragment%virtmat, respectively.
    type(decfrag),intent(inout) :: MyFragment

    ! Occ-occ block of density matrix
    call calculate_corrdens_AOS_occocc(t2,MyFragment)

    ! Virt-virt block of density matrix
    call calculate_corrdens_AOS_virtvirt(t2,MyFragment)

  end subroutine calculate_corrdens_AOS



  !> \brief Calculate occ-occ block of correlation density matrix
  !> for atomic fragment using all AOS amplitudes (see details inside subroutine).
  !> \author Kasper Kristensen
  !> \date August 2013
  subroutine calculate_corrdens_AOS_occocc(t2,MyFragment)
    implicit none
    !> Doubles amplitudes in AOS stored as (a,i,b,j)
    type(array4),intent(in) :: t2
    !> MyFragment - output occ and virt density matrices are stored in
    !> MyFragment%occmat and MyFragment%virtmat, respectively.
    type(decfrag),intent(inout) :: MyFragment
    integer :: i,j,k,a,b
    real(realk) :: i2,i4

    i2 = 2.0_realk
    i4 = 4.0_realk

    ! Occ-occ block of density matrix
    ! *******************************

    ! OccMat(i,j) = sum_{abk} t_{ik}^{ab}  tbar_{jk}^{ab}
    ! where all indices are in AOS
    ! tbar_{ij}^{ab} = 4t_{ij}^{ab} - 2t_{ij}^{ba}
    MyFragment%OccMat = 0.0_realk
    do b=1,t2%dims(1)
       do a=1,t2%dims(1)
          do k=1,t2%dims(2)
             do i=1,t2%dims(2)
                do j=1,t2%dims(2)
                   MyFragment%OccMat(i,j) = MyFragment%OccMat(i,j) + &
                        & t2%val(a,i,b,k)*(i4*t2%val(a,j,b,k) - i2*t2%val(b,j,a,k))
                end do
             end do
          end do
       end do
    end do


  end subroutine calculate_corrdens_AOS_occocc




  !> \brief Calculate virt-virt block of correlation density matrix
  !> for atomic fragment using all AOS amplitudes (see details inside subroutine).
  !> \author Kasper Kristensen
  !> \date August 2013
  subroutine calculate_corrdens_AOS_virtvirt(t2,MyFragment)
    implicit none
    !> Doubles amplitudes in AOS stored as (a,i,b,j)
    type(array4),intent(in) :: t2
    !> MyFragment - output occ and virt density matrices are stored in
    !> MyFragment%occmat and MyFragment%virtmat, respectively.
    type(decfrag),intent(inout) :: MyFragment
    integer :: i,j,a,b,c
    real(realk) :: i2,i4

    i2 = 2.0_realk
    i4 = 4.0_realk


    ! Virt-virt block of density matrix
    ! *********************************

    ! VirtMat(a,b) = sum_{ijc} t_{ij}^{ac}  tbar_{ij}^{bc}
    ! where all indices are in AOS.
    MyFragment%VirtMat = 0.0_realk
    do i=1,t2%dims(2)
       do j=1,t2%dims(2)
          do a=1,t2%dims(1)
             do b=1,t2%dims(1)
                do c=1,t2%dims(1)
                   MyFragment%VirtMat(a,b) = MyFragment%VirtMat(a,b) + &
                        & t2%val(a,i,c,j)*(i4*t2%val(b,i,c,j) - i2*t2%val(c,i,b,j))
                end do
             end do
          end do
       end do
    end do


  end subroutine calculate_corrdens_AOS_virtvirt

  !> \brief Set threshold for removing individual FOs from atomic fragments P and Q
  !> before merging the atomic fragments to generate pair fragment PQ.
  !> This is an active research area, and the optimal procedure is yet unknown!
  !> \author Kasper Kristensen
  !> \date September 2013
  function set_pairFOthr(pairdist) result(pairFOthr)
    implicit none
    !> Distance between atoms P and Q (a.u.)
    real(realk), intent(in) :: pairdist
    !> Threshold for removing FOs on atomic fragments P and Q for occupied (first entry)
    !> and virtual (second entry) FOs.
    real(realk) :: pairFOthr(2)

    ! In the long run, pairdist and the FOT should be used to dictate
    ! the value of pairFOthr. For now, we simply set using input keyword for virtual space
    ! and 0 for occupied space
    pairFOthr(1) = 0.0_realk
    pairFOthr(2) = DECinfo%pairFOthr
    write(DECinfo%output,'(a,2g20.5)') 'Setting pairFOthr to ', pairFOthr

  end function set_pairFOthr


  !> \brief When using fragment-adapted orbitals: Get dimensions of pair fragment PQ as unions
  !> of spaces for atomic fragments P and Q,
  !> where only FOs in the original atomic fragments with eigenvalues larger than pairFOthr are kept.
  !> NOTE: We take the EOS orbitals completely out of the treatment at this stage since they are not
  !> allowed to mix with the XOS orbitals at later stages. Thus, the dimensions do NOT include EOS
  !> orbitals and the logical vectors are false for the EOS orbitals.
  !> \author Kasper Kristensen
  !> \date September 2013
  subroutine get_pairFO_union(fragmentP,fragmentQ,fragmentPQ,noccPQ,nunoccPQ,nbasisPQ,&
       & WhichOccP, WhichOccQ, WhichUnoccP, WhichUnoccQ)
    implicit none
    !> Fragment P
    type(decfrag),intent(in) :: fragmentP
    !> Fragment Q
    type(decfrag),intent(in) :: fragmentQ
    !> Pair fragment PQ using local orbitals
    type(decfrag),intent(in) :: FragmentPQ
    !> Number of occupied, unoccupied MOs and number of atomic basis functions for pair fragment
    integer,intent(inout) :: noccPQ,nunoccPQ,nbasisPQ
    !> Which occupied FOs to include from atomic fragments P and Q
    logical,intent(inout) :: WhichOccP(fragmentP%noccFA), WhichOccQ(fragmentQ%noccFA)
    !> Which unoccupied FOs to include from atomic fragments P and Q
    logical,intent(inout) :: WhichUnoccP(fragmentP%nunoccFA), WhichUnoccQ(fragmentQ%nunoccFA)
    integer :: i, noccP,noccQ,nunoccP,nunoccQ,maxidx
    real(realk) :: themax


    ! Number of FOs for atomic fragment P with eigenvalues below threshold
    ! ********************************************************************

    ! Threshold is stored in FragmentPQ%RejectThr (first/second entry refer to occ/unocc threshold)


    ! OCCUPIED
    ! ========

    ! Fragment P, occupied
    noccP=0
    WhichOccP=.false.
    themax = 0.0_realk
    do i=fragmentP%noccEOS+1,fragmentP%noccFA  ! Skip EOS orbitals in loop (first noccEOS orbitals)
       if( abs(fragmentP%CDocceival(i)) > fragmentPQ%rejectthr(1) ) then
          noccP = noccP + 1
          WhichOccP(i) = .true.
       end if
       ! Find max eigenvalue in case sanity check below needs to be invoked
       if( abs(fragmentP%CDocceival(i)) > themax ) then
          themax = abs(fragmentP%CDocceival(i))
          maxidx = i
       end if
    end do
    ! Sanity check, ensure that we have at least one pair FO, choose the one with largest eigenvalue
    if(count(WhichOccP)==0) then
       WhichOccP(maxidx) = .true.
    end if

    ! Fragment Q, occupied
    noccQ=0
    WhichOccQ=.false.
    themax = 0.0_realk
    do i=fragmentQ%noccEOS+1,fragmentQ%noccFA
       if( abs(fragmentQ%CDocceival(i)) > fragmentPQ%rejectthr(1) ) then
          noccQ = noccQ + 1
          WhichOccQ(i) = .true.
       end if
       ! Find max eigenvalue in case sanity check below needs to be invoked
       if( abs(fragmentQ%CDocceival(i)) > themax ) then
          themax = abs(fragmentQ%CDocceival(i))
          maxidx = i
       end if
    end do
    ! Sanity check, ensure that we have at least one pair FO, choose the one with largest eigenvalue
    if(count(WhichOccQ)==0) then
       WhichOccQ(maxidx) = .true.
    end if

    ! Number of occupied PQ XOS orbitals is the union of XOS orbitals above threshold for P and Q.
    noccPQ = count(WhichOccP) + count(WhichOccQ)



    ! UNOCCUPIED
    ! ==========

    ! Fragment P, unoccupied
    nunoccP=0
    WhichUnoccP=.false.
    themax = 0.0_realk
    do i=fragmentP%nunoccEOS+1,fragmentP%nunoccFA  ! Skip EOS orbitals in loop (first nunoccEOS orbitals)
       if( abs(fragmentP%CDunocceival(i)) > fragmentPQ%rejectthr(2) ) then
          nunoccP = nunoccP + 1
          WhichUnoccP(i) = .true.
       end if
       ! Find max eigenvalue in case sanity check below needs to be invoked
       if( abs(fragmentP%CDunocceival(i)) > themax ) then
          themax = abs(fragmentP%CDunocceival(i))
          maxidx = i
       end if
    end do
    ! Sanity check, ensure that we have at least one pair FO, choose the one with largest eigenvalue
    if(count(WhichUnoccP)==0) then
       WhichUnoccP(maxidx) = .true.
    end if

    ! Fragment Q, unoccupied
    nunoccQ=0
    WhichUnoccQ=.false.
    themax = 0.0_realk
    do i=fragmentQ%nunoccEOS+1,fragmentQ%nunoccFA
       if( abs(fragmentQ%CDunocceival(i)) > fragmentPQ%rejectthr(2) ) then
          nunoccQ = nunoccQ + 1
          WhichUnoccQ(i) = .true.
       end if
       ! Find max eigenvalue in case sanity check below needs to be invoked
       if( abs(fragmentQ%CDunocceival(i)) > themax ) then
          themax = abs(fragmentQ%CDunocceival(i))
          maxidx = i
       end if
    end do
    ! Sanity check, ensure that we have at least one pair FO, choose the one with largest eigenvalue
    if(count(WhichUnoccQ)==0) then
       WhichUnoccQ(maxidx) = .true.
    end if


    ! Number of unoccupied PQ XOS orbitals is the union of XOS orbitals above threshold for P and Q.
    nunoccPQ = count(WhichUnoccP) + count(WhichUnoccQ)



    ! Number of basis functions
    ! --> same as for local fragment, simply copy dimension
    nbasisPQ = fragmentPQ%nbasis


  end subroutine get_pairFO_union


  !> \brief For each pair, determine whether pair calculation should be calculated
  !> should be calculated (i) using input CC model, (ii) MP2 model, or (iii) simply be skipped.
  !> This information is stored in MyMolecule%ccmodel(P,Q) in terms of an integer
  !> defining the model to be used for pair (P,Q), according to the model.
  !> See details defining the model inside subroutine.
  !> definitions in dec_typedef.F90.
  !> \author Kasper Kristensen
  !> \date October 2013
  subroutine define_pair_calculations(natoms,dofrag,FragEnergies,MyMolecule,Epair_est,Eskip_est)
    implicit none
    !> Number of atoms in molecule
    integer,intent(in) :: natoms
    !> Logical vector describing which atoms have orbitals assigned
    !> (i.e., which atoms to consider in atomic fragment calculations)
    logical,intent(in) :: dofrag(natoms)
    !> Estimated fragment energies
    real(realk),intent(in) :: FragEnergies(natoms,natoms)
    !> Full molecule structure, where only MyMolecule%ccmodel is modified
    type(fullmolecule),intent(inout) :: MyMolecule
    !> Estimated correlation energy from estimated pair fragments
    real(realk),intent(inout) :: Epair_est
    !> Estimated correlation energy from the pairs that will be skipped in the calculation
    real(realk),intent(inout) :: Eskip_est
    integer :: P,Q,Qidx,Qstart,Nskip,NMP2,NCC,Qquit,Pmax,Qmax
    real(realk),dimension(natoms) :: EPQ
    real(realk) :: Eacc,Epairmax
    integer,dimension(natoms) :: atomindices


    ! ****************************************************************************************
    !            DIFFERENT TREATMENT OF PAIR FRAGMENTS BASED ON ENERGY ESTIMATES
    ! ****************************************************************************************
    !
    ! The basic philosophy is that, since we have an error of ~FOT in the atomic fragment energy E_P,
    ! we also allow error(s) of size ~FOT in the pair fragment energies dE_PQ associated with
    ! atomic site P.
    ! We use usually write the DEC correlation energy like this:
    !
    ! Ecorr = sum_P E_P  +  sum_P sum_{Q<P} dE_PQ
    !
    ! For this analysis it is more convinient to rewrite it like:
    !
    ! Ecorr = sum_P [ E_P + 1/2 sum_P sum_{Q .neq. P} dE_PQ ]
    !
    ! In this way everything inside [...] is associated with atom P, and all atoms have
    ! the same number of associated pair fragments. The procedure use is then:
    !
    ! 1. For each atomic site P, sort estimated pair energies (1/2)dE_PQ according to size.
    ! 2. Add up the smallest contributions until they add up to the FOT. Skip those pairs.
    ! 3. Add up the "second-smallest" contributions until they add up to the FOT.
    !    Calculate these pairs at the MP2 level.
    !
    ! 4. Thus, at the end we get this where the (absolute)
    !    pair energies are arranged in decreasing order:
    !
    !
    !    (1/2) dE_PQ: (*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*,*)
    !                 |          CC MODEL             |    MP2 LEVEL  |    SKIP PAIRS   |
    !
    !
    ! In this way we (2) skip pairs for which the contribution to the energy is so small
    ! that it anyways drowns in numerical noise from the E_P FOT error; and then (3) we allow an
    ! additional error of order FOT by treating the next pairs only at the MP2 level. Point (3) is
    ! a very conservative strategy, where we basically assume that the error introduced by doing MP2
    ! rather than the actual CC model is of the same order as the actual pair interaction energy.
    ! In practice this error is typically 1%-10% of the pair interaction energy itself,
    ! but for now we choose to be very conservative here...


    ! Init stuff
    MyMolecule%ccmodel = DECinfo%ccmodel  ! use original CC model as initialization model


    ! Loop over all atoms P
    Ploop: do P=1,natoms

       ! If no orbitals assigned to P, no pairs for this atom
       if(.not. dofrag(P)) then
          Qloop1: do Q=1,natoms
             MyMolecule%ccmodel(P,Q)=MODEL_NONE
             MyMolecule%ccmodel(Q,P)=MODEL_NONE
          end do Qloop1
          cycle Ploop
       end if

       ! Sort absolute pair interaction energies dE_PQ for given P and all Q's.
       ! Multiply by 1/2 for the reasons discussed above.
       EPQ = 0.5_realk*abs(FragEnergies(:,P))
       call real_inv_sort_with_tracking(EPQ,atomindices,natoms)


       ! Step 2 above: Find which pairs to skip
       ! **************************************
       Eacc = 0.0_realk
       Qloop2: do Q=natoms,1,-1

          ! Skip if no orbitals assigned to atom Q
          Qidx = atomindices(Q)
          if(.not. dofrag(Qidx) .or. Qidx==P) cycle Qloop2

          ! Accumulated estimated energy from smallest pair interaction energies
          Eacc = Eacc + EPQ(Q)
          if(Eacc < DECinfo%FOT) then
             ! Accumulated value smaller than FOT --> (P,Qidx) can be skipped!
             MyMolecule%ccmodel(P,Qidx)=MODEL_NONE
          else
             ! Not possible to skip (P,Q). Save current counter value.
             Qquit = Q
             exit Qloop2
          end if

       end do Qloop2


       ! Step 3 above: Find which pairs to treat at MP2 level
       ! ****************************************************
       ! (Only do this if DECinfo%PairMP2 is turned on)

       NotMP2model: if(DECinfo%PairMP2 .and. DECinfo%ccmodel/=MODEL_MP2) then

          ! Same procedure as step 2 above, now the loop starts where we ended above...
          Eacc = 0.0_realk
          Qstart=Qquit
          Qloop3: do Q=Qstart,1,-1
             Qidx = atomindices(Q)
             if(.not. dofrag(Qidx) .or. Qidx==P) cycle Qloop3

             Eacc = Eacc + EPQ(Q)
             if(Eacc < DECinfo%FOT) then
                MyMolecule%ccmodel(P,Qidx)=MODEL_MP2
             else
                Qquit = Q
                exit Qloop3
             end if

          end do Qloop3

       end if NotMP2model

       ! Sanity precaution: Atomic fragment energies should always use input CC model
       MyMolecule%ccmodel(P,P) = DECinfo%ccmodel

    end do Ploop

    ! Fix inconsistencies which may arise if MyMolecule%ccmodel(P,Q) is not
    ! identical to MyMolecule%ccmodel(Q,P) from the above procedure
    call fix_inconsistencies_for_pair_model(MyMolecule)

    ! Count number of pairs to be treated at the different levels
    NCC=0
    NMP2=0
    Nskip=0
    Epairmax=0.0_realk
    Pmax=0
    Qmax=0
    do P=1,natoms
       if(.not. dofrag(P)) cycle
       do Q=1,P-1
          if(.not. dofrag(Q)) cycle
          if(MyMolecule%ccmodel(P,Q)==DECinfo%ccmodel) then
             NCC=NCC+1
          elseif(MyMolecule%ccmodel(P,Q)==MODEL_MP2) then
             NMP2=NMP2+1
          elseif(MyMolecule%ccmodel(P,Q)==MODEL_NONE) then
             Nskip=Nskip+1
          end if

          ! Find max pair energy and associated indices (for sanity check below)
          if(Epairmax < abs(FragEnergies(P,Q))) then
             Epairmax = abs(FragEnergies(P,Q))
             Pmax=P
             Qmax=Q
          end if

       end do
    end do

    ! Sanity check to avoid empty job list --> Always have at least one pair
    ! (only relevant for debug molecules)
    if(NCC+NMP2==0) then
       ! Simply calculate the pair with the largest (absolute) estimated energy
       MyMolecule%ccmodel(Pmax,Qmax) = DECinfo%ccmodel
       MyMolecule%ccmodel(Qmax,Pmax) = DECinfo%ccmodel
       NCC=1
       Nskip = Nskip-1
    end if

    ! Estimate correlation energy by adding all estimated atomic and pair fragment energy contributions
    call add_dec_energies(natoms,FragEnergies,dofrag,Epair_est)

    ! Estimate energy contribution from pairs which will be skipped from the calculation
    call estimate_energy_of_skipped_pairs(natoms,FragEnergies,dofrag,MyMolecule,Eskip_est)

    ! Print summary for pair analysis and estimated energies
    call print_pair_estimate_summary(natoms,Nskip,NMP2,NCC,dofrag,Epair_est,Eskip_est,&
         & FragEnergies,MyMolecule%DistanceTable)

  end subroutine define_pair_calculations


  !> \brief Print summary of analysis used to define pair calculations (see define_pair_calculations).
  !> \author Kasper Kristensen
  !> \date October 2013
  subroutine print_pair_estimate_summary(natoms,Nskip,NMP2,NCC,dofrag,Epair_est,Eskip_est,&
       & FragEnergies,DistanceTable)
    implicit none
    !> Number of atoms in molecule
    integer,intent(in) :: natoms
    !> Number of pairs skipped
    integer,intent(in) :: Nskip
    !> Number of pairs treated at MP2 level
    integer,intent(in) :: NMP2
    !> Number of pairs treated at input CC level
    integer,intent(in) :: NCC
    !> Which atomic sites have orbitals assigned?
    logical,intent(in),dimension(natoms) :: dofrag
    !> Estimated correlation energy from estimated pair fragments
    real(realk),intent(in) :: Epair_est
    !> Estimated correlation energy from the pairs that will be skipped in the calculation
    real(realk),intent(in) :: Eskip_est
    !> Estimated fragment energies
    real(realk),intent(in) :: FragEnergies(natoms,natoms)
    !> Interatomic distance table
    real(realk),intent(in) :: DistanceTable(natoms,natoms)
    integer :: npairs

    ! Total number of pairs
    npairs = count(dofrag)*(count(dofrag)-1)/2

    write(DECinfo%output,'(1X,a)') '******************************************************************'
    write(DECinfo%output,'(1X,a)') '*            SUMMARY FOR PAIR ESTIMATE ANALYSIS                  *'
    write(DECinfo%output,'(1X,a)') '******************************************************************'
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a,i10)') 'Total number of pairs:                ', npairs
    write(DECinfo%output,'(1X,a,i10)') 'Pairs to be treated at input CC level ', NCC
    if(DECinfo%PairMP2 .and. DECinfo%ccmodel/=MODEL_MP2) then
       write(DECinfo%output,'(1X,a,i10)') 'Pairs to be treated at MP2 level      ', NMP2
    end if
    write(DECinfo%output,'(1X,a,i10)') 'Pairs to be skipped from calculation  ', Nskip
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a,g20.10)') 'Estimated pair correlation energy:         ', Epair_est
    write(DECinfo%output,'(1X,a,g20.10)') 'Estimated contribution from skipped pairs: ', Eskip_est
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    ! Also print estimated pair fragment energies. Maybe this should be removed at some point
    call print_pair_fragment_energies(natoms,FragEnergies,dofrag,&
         & DistanceTable, 'Estimated occupied pair energies','PF_ESTIMATE')
    write(DECinfo%output,*)
    write(DECinfo%output,*)

  end subroutine print_pair_estimate_summary

  !> Once the CC model to be used for each pair has been defined in define_pair_calculations,
  !> there will in general be inconsistencies. For example, it might be that from P's point
  !> of view the pair (P,Q) should be treated at the CCSD level [MyMolecule%ccmodel(P,Q)=MODEL_CCSD],
  !> while from the point of view of Q, this pair should be skipped
  !> [MyMolecule%ccmodel(Q,P)=MODEL_NONE]. In such cases we always choose the more accurate model.
  !> (Thus, for this example we would do a CCSD calculation for pair (P,Q)).
  !> \author Kasper Kristensen
  !> \date October 2013
  subroutine fix_inconsistencies_for_pair_model(MyMolecule)
    implicit none
    !> Full molecule structure, where only MyMolecule%ccmodel is modified
    type(fullmolecule),intent(inout) :: MyMolecule
    integer :: P,Q

    do P=1,MyMolecule%natoms
       do Q=1,MyMolecule%natoms
          ! Check whether model defined by (P,Q) is different from that defined by (Q,P)
          Inconsistency: if(MyMolecule%ccmodel(P,Q) /= MyMolecule%ccmodel(Q,P)) then

             ! Possibility 1: Either (P,Q) or (Q,P) requires original CC model, the
             ! other one requires a less accurate model
             if(MyMolecule%ccmodel(P,Q)==DECinfo%ccmodel .or. &
                  & MyMolecule%ccmodel(Q,P)==DECinfo%ccmodel) then
                MyMolecule%ccmodel(P,Q)=DECinfo%ccmodel
                MyMolecule%ccmodel(Q,P)=DECinfo%ccmodel

                ! Possibility 2: Either (P,Q) or (Q,P) requires MP2 model, the
                ! other one dictates that the pair should be skipped
             elseif(MyMolecule%ccmodel(P,Q)==MODEL_MP2 .or. &
                  & MyMolecule%ccmodel(Q,P)==MODEL_MP2) then
                MyMolecule%ccmodel(P,Q)=MODEL_MP2
                MyMolecule%ccmodel(Q,P)=MODEL_MP2

                ! This should never happen! Something very wrong...
             else
                print *, 'P,Q: ', P,Q
                print *, 'Model(P,Q): ', MyMolecule%ccmodel(P,Q)
                print *, 'Model(Q,P): ', MyMolecule%ccmodel(Q,P)
                call lsquit('fix_inconsistencies_for_pair_model: Something wrong &
                     & with definitions of CC model for pair!',-1)
             end if

          end if Inconsistency

       end do
    end do

  end subroutine fix_inconsistencies_for_pair_model

end module atomic_fragment_operations
