!> @file
!> Contains orbital type and subroutines need to approimate it
!> \author Marcin Ziolkowski
module orbital_operations


  ! Outside DEC directory
  use fundamental
  use memory_handling
  use precision
  use matrix_module!, only:matrix
  use matrix_operations
  use matrix_util!,only: util_AO_to_MO_different_trans
  use lstiming!, only: lstimer
  use files!,only: lsopen,lsclose
  use typedeftype!,only:lsitem
  use dec_typedef_module
  use lsparameters
  use IntegralInterfaceMOD


  ! DEC DEPENDENCIES (within deccc directory)    
  ! *****************************************
  use dec_fragment_utils


contains


  !> \brief Generate ocupied and virtual orbitals.
  !> \author Kasper Kristensen
  !> \date November 2011
  subroutine GenerateOrbitals_driver(MyMolecule,mylsitem,nocc,nvirt,natoms, &
       & OccOrbitals, UnoccOrbitals)

    implicit none
    !> Full molecule structure ( Note: MyMolecule is only changed if modbasis=.true.!)
    type(fullmolecule), intent(in) :: MyMolecule
    !> LSitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> Number of occupied orbitals
    integer,intent(in) :: nocc
    !> Number of virtupied orbitals
    integer,intent(in) :: nvirt
    !> Number of atoms in the molecule
    integer,intent(in) :: natoms
    !> Occupied orbitals to be generated
    type(decorbital),intent(inout) :: OccOrbitals(nocc)
    !> Unoccupied orbitals to be generated
    type(decorbital),intent(inout) :: UnoccOrbitals(nvirt)


    DECCO: if(DECinfo%DECCO) then
       ! Orbital-based DEC
       call GenerateOrbitals_DECCO(nocc,nvirt,natoms, &
          & MyMolecule,MyLsitem,DECinfo%simple_orbital_threshold,OccOrbitals,UnoccOrbitals)
       return

    else

       ! Atom-based
       OrbitalGeneration: if(DECinfo%read_dec_orbitals) then ! Read DEC orbitals form file
          call read_DECorbitals_from_file(nocc,nvirt,&
             &OccOrbitals,UnoccOrbitals)
       else ! Generate DEC orbitals

          write(DECinfo%output,*) 'Generating DEC orbitals using simple Lowdin charge analysis'
          
          call GenerateOrbitals_simple(nocc,nvirt,natoms, &
               & MyMolecule,MyLsitem,DECinfo%simple_orbital_threshold,OccOrbitals,UnoccOrbitals)

       end if OrbitalGeneration

    end if DECCO

    ! If we want to simulate full calculation, we simply assign ALL orbitals to atom 1
    if(DECinfo%simulate_full) then
       write(DECinfo%output,*) 'THIS IS A SIMULATED FULL MOLECULAR CALCULATION!'
       call adjust_orbitals_for_full_simulation(nocc,nvirt,&
            &OccOrbitals,UnoccOrbitals,natoms,DECinfo%simulate_natoms)
    end if

    ! Set secondary atom for orbitals (see set_secondary_atom_for_orbitals)
    call set_secondary_atom_for_orbitals(natoms,nocc,nvirt,MyMolecule,mylsitem,&
         & OccOrbitals,UnoccOrbitals)

    ! Print orbital info
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) 'Simple Lowdin-based orbital assignment: Summary'
    write(DECinfo%output,*) '-----------------------------------------------'
    call print_orbital_info(mylsitem,nocc,natoms,nvirt,OccOrbitals,UnoccOrbitals,ncore=MyMolecule%ncore)


    ! Check that assigment is meaningful
    call dec_orbital_sanity_check(natoms,nocc,nvirt,OccOrbitals,&
         & UnoccOrbitals,MyMolecule)

    ! Write orbitals to file
    call write_DECorbitals_to_file(nocc,nvirt,&
         &OccOrbitals,UnoccOrbitals)


  end subroutine GenerateOrbitals_driver


  !> \brief Set secondary central atom for orbitals. 
  !> For example, if virtual orbital "a" is assigned to a hydrogen atom "H_A" for which 
  !> there are no occupied orbitals assigned, the secondary          
  !> central atom for orbital "a" will be the atom closest to "H_A" which has a nonzero number
  !> of occupied AND virtual orbitals in the original assignment.
  !> \author Kasper Kristensen
  !> \date June 2014
  subroutine set_secondary_atom_for_orbitals(natoms,nocc,nvirt,MyMolecule,mylsitem,&
       & OccOrbitals,UnoccOrbitals)
    implicit none

    !> Number of atoms in the molecule
    integer,intent(in) :: natoms
    !> Number of occupied orbitals
    integer,intent(in) :: nocc
    !> Number of virtupied orbitals
    integer,intent(in) :: nvirt
    !> Full molecule structure ( Note: MyMolecule is only changed if modbasis=.true.!)
    type(fullmolecule), intent(in) :: MyMolecule
    !> LSitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> Occupied orbitals to be generated
    type(decorbital),intent(inout) :: OccOrbitals(nocc)
    !> Unoccupied orbitals to be generated
    type(decorbital),intent(inout) :: UnoccOrbitals(nvirt)
    integer :: nocc_per_atom(natoms), nvirt_per_atom(natoms),offset,P,i,central,secondary
    logical :: reset_secondary_atom(natoms),found
    real(realk) :: dummy(natoms)
    integer,pointer :: DistSortedIdx(:,:)


    !> Set offset for frozen core to not consider core orbitals in bookkeeping.
    if(DECinfo%frozencore) then
       offset=MyMolecule%ncore
    else
       offset=0
    end if

    ! Number of orbitals assigned to each atom
    nocc_per_atom =  get_number_of_orbitals_per_atom(OccOrbitals,nocc,natoms,.true.,offset=offset)
    nvirt_per_atom =  get_number_of_orbitals_per_atom(UnoccOrbitals,nvirt,natoms,.true.)

    ! First, set secondary atom equal to central atom
    do i=1,nocc
       OccOrbitals(i)%secondaryatom=OccOrbitals(i)%centralatom
    end do
    do i=1,nvirt
       UnoccOrbitals(i)%secondaryatom=UnoccOrbitals(i)%centralatom
    end do

    ! Now we determine if any atoms have zero occ and nonzero virtual orbitals assigned (or vice versa)
    ! to determine which orbitals need a new secondary atom.
    ! Also, for each atom we sort the other atoms according to distance.
    reset_secondary_atom=.false.
    call mem_alloc(DistSortedIdx,natoms,natoms)
    do P=1,natoms
       ! Reset secondary atom?
       if(nocc_per_atom(P)==0 .and. nvirt_per_atom(P)/=0) reset_secondary_atom(P)=.true.
       if(nocc_per_atom(P)/=0 .and. nvirt_per_atom(P)==0) reset_secondary_atom(P)=.true.
       
       ! Set sorted distance table
       call GetSortedList(dummy,DistSortedIdx(:,P),mymolecule%DistanceTable,natoms,natoms,P)

    end do


    ! Reassign occupied orbitals
    do i=1,nocc
       central = OccOrbitals(i)%centralatom
       if(reset_secondary_atom(central)) then

          ! Find atom closest to current central atom which has both occ and virt orbitals assigned
          found=.false.
          Ploop: do P=1,natoms
             ! Secondary atom from distance list
             secondary = DistSortedIdx(P,central)
             if(nocc_per_atom(secondary)/=0 .and. nvirt_per_atom(secondary)/=0) then
                found=.true.
                OccOrbitals(i)%secondaryatom = secondary
                exit Ploop
             end if
          end do Ploop
          ! Sanity check
          if(.not. found) then
             call print_orbital_info(mylsitem,nocc,natoms,nvirt,OccOrbitals,UnoccOrbitals,&
                  & ncore=MyMolecule%ncore)
             call lsquit('set_secondary_atom_for_orbitals: Occupied orbital reassigning failed!')
          end if

       end if
    end do

    ! Reassign virtupied orbitals
    do i=1,nvirt
       central = UnoccOrbitals(i)%centralatom
       if(reset_secondary_atom(central)) then

          ! Find atom closest to current central atom which has both occ and virt orbitals assigned
          found=.false.
          Ploop2: do P=1,natoms
             ! Secondary atom from distance list
             secondary = DistSortedIdx(P,central)
             if(nocc_per_atom(secondary)/=0 .and. nvirt_per_atom(secondary)/=0) then
                found=.true.
                UnoccOrbitals(i)%secondaryatom = secondary
                exit Ploop2
             end if
          end do Ploop2
          ! Sanity check
          if(.not. found) then
             call print_orbital_info(mylsitem,nocc,natoms,nvirt,OccOrbitals,UnoccOrbitals,&
                  & ncore=MyMolecule%ncore)
             call lsquit('set_secondary_atom_for_orbitals: Unoccupied orbital reassigning failed!')
          end if

       end if
    end do

    call mem_dealloc(DistSortedIdx)

  end subroutine set_secondary_atom_for_orbitals



  !> \brief Initialize orbital, assign center and list of significant atoms.
  function orbital_init(orb_number,central_atom,num_aos,list_aos) result(myorbital)

    implicit none
    type(decorbital) :: myorbital
    integer, intent(in) :: orb_number
    integer, intent(in) :: central_atom
    integer, intent(in) :: num_aos
    integer, dimension(num_aos), intent(in) :: list_aos
    integer :: i

    myorbital%orbitalnumber = orb_number
    myorbital%centralatom = central_atom
    myorbital%numberofaos = num_aos

    ! Set secondary atom equal to central atom in initialization (may be changed later)
    myorbital%secondaryatom = central_atom

    call mem_alloc(myorbital%aos,num_aos)
    myorbital%aos = 0

    do i=1,num_aos
       myorbital%aos(i) = list_aos(i)
    end do

  end function orbital_init


!!$  !> \brief Initialize orbital, assign center and list of significant atoms from logical array.
!!$  !> (Small wrapper for orbital_init)
!!$  function orbital_init_from_logical_array(orb_number,central_atom,natoms,which_atoms) result(myorbital)
!!$
!!$    implicit none
!!$    !> Orbital to initialize
!!$    type(decorbital) :: myorbital
!!$    !> Orbital index
!!$    integer, intent(in) :: orb_number
!!$    !> Central atom for orbital
!!$    integer, intent(in) :: central_atom
!!$    !> Number of atoms in FULL molecule
!!$    integer,intent(in) :: natoms
!!$    !> Which atoms are included in orbital extent?
!!$    logical,dimension(natoms) :: which_atoms
!!$    integer :: num_atoms
!!$    integer, pointer :: list_atoms(:)
!!$    integer :: i,idx
!!$
!!$    ! Number of atoms in orbital extent
!!$    num_atoms = count(which_atoms)
!!$
!!$    ! Sanity check
!!$    if(num_atoms ==0) then
!!$       write(DECinfo%output,*) 'No atoms for orbital: ', orb_number
!!$       call lsquit('orbital_init_from_logical_array: No atoms in orbital extent!',-1)
!!$    end if
!!$
!!$    ! Convert logical atom list to integer list of atomic indices
!!$    call mem_alloc(list_atoms,num_atoms)
!!$    idx=0
!!$    do i=1,natoms
!!$       if(which_atoms(i)) then
!!$          idx =idx+1
!!$          list_atoms(idx) = i
!!$       end if
!!$    end do
!!$
!!$    ! Init orbital
!!$    myorbital = orbital_init(orb_number,central_atom,num_atoms,list_atoms)
!!$
!!$    call mem_dealloc(list_atoms)
!!$
!!$  end function orbital_init_from_logical_array
!!$
!!$
!!$
  !> \brief Copy orbital (orbital also intialized here)
  subroutine copy_orbital(OriginalOrbital,OrbitalCopy)
    implicit none

    !> Original orbital
    type(decorbital),intent(in) :: OriginalOrbital
    !> Copy of original orbitals
    type(decorbital),intent(inout) :: OrbitalCopy

    OrbitalCopy%orbitalnumber = OriginalOrbital%orbitalnumber
    OrbitalCopy%centralatom = OriginalOrbital%centralatom
    OrbitalCopy%numberofaos = OriginalOrbital%numberofaos
    OrbitalCopy%secondaryatom = OriginalOrbital%secondaryatom

    call mem_alloc(OrbitalCopy%aos,OrbitalCopy%numberofaos)
    OrbitalCopy%aos = OriginalOrbital%aos

  end subroutine copy_orbital

  !> \brief Generate DEC orbitals for both occ and virt orbitals. For each orbital:
  !> 1. List atoms according to Lowdin charge for that orbital.
  !> 2. Include atoms from this list until "1 minus the sum of these Lowdin charges"
  !>    is smaller than the input approximated_norm_threshold.
  !> 3. Regarding orbital assigning we ensure that we have one atomic fragment for each heavy atom
  !>    + one atomic fragment for all hydrogens (if any) assigned to that heavy atom.
  !> \author Kasper Kristensen
  !> \date September 2011
  subroutine GenerateOrbitals_simple(nocc,nvirt,natoms, &
       & MyMolecule,MyLsitem,approximated_norm_threshold,OccOrbitals,UnoccOrbitals)

    implicit none
    !> Number of occupied orbitals
    integer, intent(in) :: nocc
    !> Number of virtupied orbitals
    integer, intent(in) :: nvirt
    !> Number of atoms in the molecule
    integer,intent(in) :: natoms
    !> Molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> General LS info
    type(lsitem),intent(inout) :: MyLsitem
    !> Threshold for orbital norm (see above)
    real(realk),intent(in) :: approximated_norm_threshold
    !> Occupied orbitals to create
    type(decorbital), intent(inout), dimension(nocc) :: OccOrbitals
    !> Unoccupied orbitals to create
    type(decorbital), intent(inout), dimension(nvirt) :: UnoccOrbitals
    integer :: i,j,central_atom,n,norbital_extent,nbasis,heavyatom,ni,bas,atom
    integer, pointer :: list_of_aos_to_consider(:)
    real(realk) :: error,charge,mindist,twall,tcpu,maxlowdin
    logical :: keepon,ReAssignVirtHydrogenOrbs
    real(realk), pointer :: ShalfC(:,:)
    real(realk), pointer :: lowdin_charge(:,:)
    integer, pointer :: basis_idx(:,:), atomic_idx(:,:), countocc(:), countvirt(:),central_atom2(:)
    integer :: maxidx,offset,changedatom,nreass
    logical,pointer ::which_hydrogens(:), dofrag(:)
    real(realk),pointer :: tmplowdin_charge(:),atomic_gross_charge(:)
    integer,pointer :: tmpatomic_idx(:),basToatom(:)
    integer :: nvirtperatom,II,IDX,nu,k,kk,endidx,nMO
    logical,pointer :: WhichAOs(:)

    call LSTIMER('START',tcpu,twall,DECinfo%output)

    ! Init stuff
    nbasis = MyMolecule%nbasis
    nMO = MyMolecule%nMO

    ! For orbital assignment bookkeeping, only consider valence orbitals
    offset = MyMolecule%ncore  

    call mem_alloc(lowdin_charge,nbasis,nMO)
    call mem_alloc(ShalfC,nbasis,nMO)
    call mem_alloc(basis_idx,nbasis,nMO)
    call mem_alloc(atomic_idx,natoms,nMO)
    if(DECinfo%Distance) then
       call mem_alloc(central_atom2,nMO)
    endif

    ! Get Lowdin matrix S^{1/2} C
    call Get_matrix_for_lowdin_analysis(MyMolecule, MyLsitem, ShalfC)
    
    ! *********************************************
    ! Get Lowdin charges for all molecular orbitals
    ! *********************************************
    GetLowdinCharges: do i=1,nMO

       if(DECinfo%Distance) then
          ! Use Distance criteria to determine central atom for orbital "i"  
          call GetDistanceCentralAtom(i,natoms,MyMolecule,central_atom2(i))
       endif

       ! Get vector with Lowdin charges for all atoms for orbital "i"
       call GetLowdinVector(i,nbasis,nMO,ShalfC,lowdin_charge(:,i) )

       call mem_alloc(atomic_gross_charge,natoms)
       atomic_gross_charge = 0.0E0_realk
       do atom=1,natoms
          do nu=MyMolecule%atom_start(atom),MyMolecule%atom_end(atom) ! nu \in atom
             atomic_gross_charge(atom)=atomic_gross_charge(atom)+abs(lowdin_charge(nu,i))
          end do
       enddo

       ! Sort Lowdin charges
       call real_inv_sort_with_tracking(lowdin_charge(:,i),basis_idx(:,i),nbasis)
       call real_inv_sort_with_tracking(atomic_gross_charge(:),atomic_idx(:,i),natoms)

       ! Reassign central atom for phantom atoms - however not for SNOOP
       IF(MyMolecule%PhantomAtom(atomic_idx(1,i)) .and. (.not. DECinfo%SNOOP)  ) THEN
          !the first atom in the atomic_gross_charge(:,i) list is a Phantom atom
          !we reorder the atomic_gross_charge(:,i) list to put a non Phantom atom on top
          call mem_alloc(tmplowdin_charge,nAtoms)
          call mem_alloc(tmpatomic_idx,nAtoms)
          ni=0
          do n=1,natoms
             tmplowdin_charge(n) = atomic_gross_charge(n)
             tmpatomic_idx(n) = atomic_idx(n,i)
             IF(.NOT.MyMolecule%PhantomAtom(atomic_idx(n,i)))THEN
                IF(ni.EQ.0)THEN
                   !found the first non Phantom atom in lowdin list
                   ni = n
                ENDIF
             ENDIF
          enddo
          !place the first non Phantom atom first
          atomic_gross_charge(1) = tmplowdin_charge(ni) 
          atomic_idx(1,i) = tmpatomic_idx(ni)
          !place the other atoms in list
          do n=1,ni-1                               
             atomic_gross_charge(1+n) = tmplowdin_charge(n)
             atomic_idx(1+n,i) = tmpatomic_idx(n)                
          enddo
          !number ni should not be placed because it is now number 1
          !place the rest of the atoms in the list
          do n=ni+1,nAtoms
             atomic_gross_charge(n) = tmplowdin_charge(n)
             atomic_idx(n,i) = tmpatomic_idx(n)                
          enddo
          call mem_dealloc(tmplowdin_charge)
          call mem_dealloc(tmpatomic_idx)
       ENDIF
       call mem_dealloc(atomic_gross_charge)
    end do GetLowdinCharges
    call mem_dealloc(ShalfC)

    ! *********************************************
    ! Initial orbital assignment and orbital extent
    ! *********************************************
    call mem_alloc(basToatom,nbasis)
    do atom=1,natoms
       do nu=MyMolecule%atom_start(atom),MyMolecule%atom_end(atom) ! nu \in atom
          basToatom(nu)=atom
       enddo
    enddo
    call mem_alloc(WhichAOs,nbasis)
    OrbitalLoop: do i=1,nMO
       ! Central atom is the one with largest Lowdin charge (although this may be modified below)
       if(DECinfo%Distance) then
          !central atom determined from shortest distance from MO
          central_atom = central_atom2(i)
          IF(atomic_idx(1,i).EQ.central_atom)THEN
             !do nothing the closest atom is also the one with the largest Lowdin charge
          ELSE
             !we could reorder the lowdin_charge(:,i) list to put the central atom on top
             !but we do not care
          ENDIF
       ELSE
          central_atom = atomic_idx(1,i)
       ENDIF

       !include all basis functions attached to new central atom
       charge = 0E0_realk
       WhichAOs = .FALSE.
       bas = 0
       do nu=MyMolecule%atom_start(central_atom),MyMolecule%atom_end(central_atom)
          WhichAOs(nu) = .TRUE.          
          bas = bas + 1
          kloop22: do k=1,nbasis
             IF(basis_idx(k,i).EQ.nu)THEN 
                charge = charge + lowdin_charge(k,i)                
                exit kloop22
             ENDIF
          enddo kloop22
       enddo
       IF(bas.NE.COUNT(WhichAOs))THEN
          call lsquit('AO dim mismatch in GenerateOrbitals_simple',-1)
       ENDIF
       error = 1E0_realk - charge
       norbital_extent = 0
          
       if(error < approximated_norm_threshold .or. COUNT(WhichAOs)==nbasis ) then 
          ! ao list converged for orbital i
          ! Set list of aos to consider for orbital and exit loop
          norbital_extent = bas
          call mem_alloc(list_of_aos_to_consider,norbital_extent)
          !include all basis functions attached to central atom
          bas = 0
          do nu=MyMolecule%atom_start(central_atom),MyMolecule%atom_end(central_atom)
             bas = bas + 1
             list_of_aos_to_consider(bas) = nu
          enddo
       else
          ! Add basis functions until sum of Lowdin charges is close enough to 1
          k=0
          LowdinAddLoop3: do n=1,nbasis
             !ensure the AO is not already included 
             IF(.NOT.WhichAOs(basis_idx(n,i)))THEN
                WhichAOs(basis_idx(n,i)) = .TRUE.
                charge = charge + lowdin_charge(n,i)
                !also include all components of this shell
                IDX = basis_idx(n,i)
                DO J=MyMolecule%bas_Start(IDX),MyMolecule%bas_End(IDX)
                   IF(.NOT.WhichAOs(J))THEN !if not already included 
                      WhichAOs(J) = .TRUE.
                      kloop23: do k=1,nbasis
                         IF(basis_idx(k,i).EQ.J)THEN 
                            charge = charge + lowdin_charge(k,i)
                            exit kloop23
                         ENDIF
                      enddo kloop23
                   ENDIF
                ENDDO
                error = 1E0_realk - charge
                if(error < approximated_norm_threshold .or. COUNT(WhichAOs)==nbasis) then 
                   ! atom list converged for orbital i
                   exit LowdinAddLoop3
                elseif(error .GT. 6.0E0_realk*approximated_norm_threshold) then 
                   !include all basis functions attached to this atom
                   atom = basToatom(basis_idx(n,i))
                   do nu=MyMolecule%atom_start(atom),MyMolecule%atom_end(atom)
                      IF(.NOT.WhichAOs(nu))THEN !if not already included 
                         WhichAOs(nu) = .TRUE.
                         kloop24: do k=1,nbasis
                            IF(basis_idx(k,i).EQ.nu)THEN 
                               charge = charge + lowdin_charge(k,i)
                               exit kloop24
                            ENDIF
                         enddo kloop24
                      ENDIF
                   enddo
                end if
             endif
          end do LowdinAddLoop3                   
          ! Set list of aos to consider for orbital and exit loop
          norbital_extent = COUNT(WhichAOs)
          call mem_alloc(list_of_aos_to_consider,norbital_extent)
          bas = 0
          do k=1,nbasis
             IF(WhichAOs(k))THEN
                bas = bas + 1
                list_of_aos_to_consider(bas) = k
             ENDIF
          enddo
       endif
       !SANITY CHECK Include full Shell (all 3 P x,y,z components) if one part of shell have been included
       DO II=1,norbital_extent
          IDX = list_of_aos_to_consider(II)
          WhichAOs(IDX) = .TRUE.
          DO J=MyMolecule%bas_Start(IDX),MyMolecule%bas_End(IDX)
             WhichAOs(J) = .TRUE.
          ENDDO
       ENDDO
       IF(norbital_extent.NE.COUNT(WhichAOs))THEN
          if(DECinfo%PL>1) then
             charge=0.0E0_realk
          endif
          !add functions
          call mem_dealloc(list_of_aos_to_consider)
          norbital_extent = COUNT(WhichAOs)
          call mem_alloc(list_of_aos_to_consider,norbital_extent)
          II = 1
          do IDX=1,nbasis
             IF(WhichAOs(IDX))THEN
                list_of_aos_to_consider(II) = IDX
                II = II + 1
                if(DECinfo%PL>1) then
                   kloop21: do k=1,nbasis
                      IF(basis_idx(k,i).EQ.IDX)THEN 
                         charge = charge + lowdin_charge(k,i)
                         exit kloop21
                      ENDIF
                   enddo kloop21
                endif
             ENDIF
          enddo
          IF(II-1.NE.norbital_extent)call lsquit('DEC error in AO shell sanitycheck',-1)
       ENDIF

       ! Print orbital info for high print levels
       if(DECinfo%PL>1) then
          write(DECinfo%output,'(1X,a,i10)') 'ORBITAL: ', i
          write(DECinfo%output,*) '-------------------------------------'
          write(DECinfo%output,'(1X,a,1i5)')        '#AOS   : ', norbital_extent
          do j=1,norbital_extent,10
             endidx=min(j+9,norbital_extent)
             write(DECinfo%output,'(1X,10i7)') list_of_aos_to_consider(j:endidx)
          end do

!          write(DECinfo%output,'(1X,a,100f10.3)')   'LOWDIN : ', lowdin_charge(1:norbital_extent,i) !not true anymore
          write(DECinfo%output,'(1X,a,f12.5)')      'TOTAL LOWDIN : ', charge 
          write(DECinfo%output,'(1X,a,i6)')         'Central Atom : ', central_atom
          if(DECinfo%Distance) then
             write(DECinfo%output,'(1X,a,i6)') 'Central Atom : ', central_atom2(i)
          endif
          write(DECinfo%output,*)
       end if

       ! -- Create orbital
       if(i<=nocc) then   ! Occupied orbital
          OccOrbitals(i) = orbital_init(i,central_atom, &
               norbital_extent,list_of_aos_to_consider)
       else  ! virtupied orbital
          UnoccOrbitals(i-nocc) = orbital_init(i,central_atom, &
               norbital_extent,list_of_aos_to_consider)
       end if

       call mem_dealloc(list_of_aos_to_consider)

    end do OrbitalLoop
    call mem_dealloc(WhichAOs)
    call mem_dealloc(basToatom)


    ! Which atoms are hydrogens?
    call mem_alloc(which_hydrogens,natoms)
    which_hydrogens=.false.
    do atom=1,natoms
       if(myLsitem%input%molecule%Atom(atom)%atomic_number==1) then
          which_hydrogens(atom)=.true.
       end if
    end do

    if(count(which_hydrogens)==natoms .and. (.NOT.decinfo%PureHydrogendebug) &
         & .and. (.not. DECinfo%OnlyOccPart.AND..not. DECinfo%OnlyVirtPart) ) then
       print*,'Orbital assignment failed because there are only hydrogen atoms!'
       print*,'For development & debug purposes the keyword PUREHYDROGENDEBUG can be used.'
       call lsquit('Orbital assignment failed because there are only hydrogen atoms!',-1)
    end if

    ! Assign each hydrogen atom to a heavy atom
    AbsorbHydrogenAtoms: if(DECinfo%AbsorbHatoms .and. (.NOT.decinfo%PureHydrogendebug)) then 
       ! Reassign orbitals originally assigned to hydrogen
       call reassign_orbitals(nocc,OccOrbitals,natoms,MyMolecule%DistanceTable,mylsitem)
       call reassign_orbitals(nvirt,UnOccOrbitals,natoms,MyMolecule%DistanceTable,mylsitem)
    end if AbsorbHydrogenAtoms


    HydrogenDebug: if( decinfo%PureHydrogendebug ) THEN

       do i=1,nocc
          OccOrbitals(i)%centralatom = i
       enddo
       nvirtperatom = nvirt/nocc
       do i=1,nvirt
          ReAssignVirtHydrogenOrbs = .TRUE.
          do j=1,nocc
             IF(UnOccOrbitals(i)%centralatom.EQ.j)THEN
                ReAssignVirtHydrogenOrbs = .FALSE.                
             ENDIF
          enddo
          IF(ReAssignVirtHydrogenOrbs)THEN
             do j=1,nocc   
                UnOccOrbitals(i)%centralatom = (nvirt-1)/nvirtperatom + 1
             enddo
          ENDIF
       enddo

    end if HydrogenDebug


    ! ****************************************************************************************
    ! * Reassign: Ensure that all atoms have both occupied and virtupied orbitals assigned  *
    ! ****************************************************************************************

    REASSIGNING: if( .not. decinfo%PureHydrogendebug ) then

       ! Count # orbitals assigned to each atom
       call mem_alloc(countocc,natoms)
       call mem_alloc(countvirt,natoms)
       countocc=0
       countvirt=0
       do i=offset+1,nocc
          countocc(OccOrbitals(i)%centralatom) = countocc(OccOrbitals(i)%centralatom)+1
       end do
       do i=1,nvirt
          countvirt(UnoccOrbitals(i)%centralatom) = countvirt(UnoccOrbitals(i)%centralatom)+1
       end do

       ! The atoms which will be central in atomic fragments are those which
       ! at this point have some orbitals assigned
       call mem_alloc(dofrag,natoms)
       dofrag=.false.
       do i=1,natoms
          if( countocc(i)/=0 .or. countvirt(i)/=0 ) dofrag(i)=.true.
       end do

       ! Now reassign to ensure that all atomic fragment have both occupied and virtupied orbitals
       keepon = .true.
       nreass=0  ! number of reassigment steps
       ContinueReassign: do while(keepon)
          nreass = nreass+1

          ReassignAtomLoop: do atom=1,natoms

             ! Reassign occupied orbitals          
             ! Never reassign occupied orbitals for only occupied partitioning
             OccReassign: if(dofrag(atom) .and. countocc(atom)==0 &
                  & .and. (.not. DECinfo%onlyoccpart) ) then

                ! Atom is supposed to be central in an atomic fragment but
                ! it has no occupied orbitals assigned:
                ! 1. Check orbitals for which atom is number 2 in the Lowdin priority list
                !    and find orbital with largest Lowdin charge in this set of orbitals
                ! 2. Reassign that orbital to atom under consideration
                maxlowdin = 0.0_realk
                maxidx = 0
                changedatom=0
                do j=offset+1,nocc

                   ! Only consider reassigning if:
                   ! (i)   Atom under consideration is number 2 in Lowdin list
                   ! (ii)  Lowdin charge for atom is larger than current max value
                   ! (iii) Orbital is not a core orbital
                   if(atomic_idx(2,j)==atom .and. lowdin_charge(2,j) > maxlowdin ) then
                      maxlowdin = lowdin_charge(2,j)
                      maxidx = j
                      changedatom=OccOrbitals(j)%centralatom  ! atom which will "loose" an orbital
                   end if
                end do

                ! Reassign orbital
                if(maxidx/=0) then
                   if(DECinfo%PL>1) write(DECinfo%output,'(1X,a,i8,a,i8)') &
                        & 'Reassign occ orbital ', maxidx, ' to atom ', atom
                   OccOrbitals(maxidx)%centralatom = atom
                   countocc(changedatom) = countocc(changedatom) - 1
                   countocc(atom) = countocc(atom) + 1
                end if

             end if OccReassign


             ! Reassign virtupied orbitals (same procedure as for occ space)
             ! Never reassign virtual orbitals for only virtual partitioning
             UnoccReassign: If(dofrag(atom) .and. countvirt(atom)==0 &
                  & .and. (.not. DECinfo%onlyvirtpart) ) then

                maxlowdin = 0.0_realk
                maxidx = 0
                changedatom=0
                do j=nocc+1,nMO
                   if(atomic_idx(2,j)==atom .and. lowdin_charge(2,j) > maxlowdin ) then
                      maxlowdin = lowdin_charge(2,j)
                      maxidx = j-nocc   ! virtupied index in list of virtupied orbitals
                      changedatom=UnoccOrbitals(j-nocc)%centralatom
                   end if
                end do

                ! Reassign orbital
                if(maxidx/=0) then
                   if(DECinfo%PL>1) write(DECinfo%output,'(1X,a,i8,a,i8)') &
                        & 'Reassign virt orbital ', maxidx, ' to atom ', atom
                   UnoccOrbitals(maxidx)%centralatom = atom
                   countvirt(changedatom) = countvirt(changedatom) - 1
                   countvirt(atom) = countvirt(atom) + 1
                end if

             end if UnoccReassign


             ! Check that all atoms have either have
             ! (i) both occupied AND virtupied orbitals assigned   OR
             ! (ii) zero orbitals assigned
             keepon=.false.
             CheckAssignment: do i=1,natoms

                WhichScheme: if(DECinfo%onlyoccpart) then

                   if( (countocc(i)/=0 .and. countvirt(i)==0) ) then
                      ! Still not acceptable orbital distribution - keep on
                      keepon=.true.
                      exit CheckAssignment
                   end if

                elseif(DECinfo%onlyvirtpart) then

                   if( (countocc(i)==0 .and. countvirt(i)/=0) ) then
                      ! Still not acceptable orbital distribution - keep on
                      keepon=.true.
                      exit CheckAssignment
                   end if

                else

                   if( (countocc(i)/=0 .and. countvirt(i)==0) .or. &
                        & (countocc(i)==0 .and. countvirt(i)/=0) ) then
                      ! Still not acceptable orbital distribution - keep on
                      keepon=.true.
                      exit CheckAssignment
                   end if

                end if WhichScheme

             end do CheckAssignment

          end do ReassignAtomLoop

          ! Avoid infinite loop
          if(keepon .and. nreass>5) then
             if(count(which_hydrogens)==natoms) then
                print*,'Orbital assignment failed because there are only hydrogen atoms!'
                print*,'For development & debug purposes the keyword PUREHYDROGENDEBUG can be used.'
                call lsquit('Orbital assignment failed because there are only hydrogen atoms!',-1)
             else 
                print *, 'Reassignment procedure failed!'
                print *, 'Suggestion: Remove .NOTABSORBH keyword'
                print *, 'If you are not using .NOTABSORBH - then DEC cannot be used for this system!'
                call lsquit('Reassignment procedure failed!',DECinfo%output)
             end if
          end if

       end do ContinueReassign

       call mem_dealloc(dofrag)
       call mem_dealloc(countocc)
       call mem_dealloc(countvirt)

    end if REASSIGNING


    call mem_dealloc(which_hydrogens)
    call mem_dealloc(basis_idx)
    call mem_dealloc(atomic_idx)
    call mem_dealloc(lowdin_charge)
    if(DECinfo%Distance) then
       call mem_dealloc(central_atom2)
    endif
    call LSTIMER('GenerateOrb',tcpu,twall,DECinfo%output)

  end subroutine GenerateOrbitals_simple





  !> \brief Generate DEC orbitals for both occ and virt orbitals using DECCO scheme. For each orbital:
  !> 1. List atoms according to Lowdin charge for that orbital.
  !> 2. Include orbitals from this list until "1 minus the sum of these Lowdin charges"
  !>    is smaller than the input approximated_norm_threshold.
  !> \author Kasper Kristensen
  !> \date September 2014
  subroutine GenerateOrbitals_DECCO(nocc,nvirt,natoms, &
       & MyMolecule,MyLsitem,approximated_norm_threshold,OccOrbitals,UnoccOrbitals)

    implicit none
    !> Number of occupied orbitals
    integer, intent(in) :: nocc
    !> Number of virtupied orbitals
    integer, intent(in) :: nvirt
    !> Number of atoms in the molecule
    integer,intent(in) :: natoms
    !> Molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> General LS info
    type(lsitem),intent(inout) :: MyLsitem
    !> Threshold for orbital norm (see above)
    real(realk),intent(in) :: approximated_norm_threshold
    !> Occupied orbitals to create
    type(decorbital), intent(inout), dimension(nocc) :: OccOrbitals
    !> Unoccupied orbitals to create
    type(decorbital), intent(inout), dimension(nvirt) :: UnoccOrbitals
    integer :: i,central_orbital,n,norbital_extent,nbasis,atom,j,bas
    integer, pointer :: list_of_aos_to_consider(:)
    real(realk) :: error,charge,twall,tcpu
    real(realk), pointer :: ShalfC(:,:)
    real(realk), pointer :: lowdin_charge(:,:)
    integer, pointer :: basis_idx(:,:),basToatom(:)
    integer :: offset,II,IDX,k,nu,nMO
    real(realk),pointer :: DistoccUnocc(:,:),DistOccOcc(:,:),sorted_dists(:)
    integer :: sorted_orbitals(nocc)
    logical,pointer :: WhichAOs(:)

    call LSTIMER('START',tcpu,twall,DECinfo%output)

    ! Init stuff
    nbasis = MyMolecule%nbasis
    nMO = MyMolecule%nMO

    offset = MyMolecule%ncore  

    call mem_alloc(lowdin_charge,nbasis,nMO)
    call mem_alloc(ShalfC,nbasis,nMO)
    call mem_alloc(basis_idx,nbasis,nMO)

    ! Get Lowdin matrix S^{1/2} C
    call Get_matrix_for_lowdin_analysis(MyMolecule, MyLsitem, ShalfC)


    ! ***********************************
    ! Get Lowdin charges for all orbitals
    ! ***********************************
    GetLowdinCharges: do i=1,nbasis

       ! Get vector with Lowdin charges for all atoms for orbital "i"
       call GetLowdinVector(i,nbasis,nMO,ShalfC,lowdin_charge(:,i) )

       ! Sort Lowdin charges
       call real_inv_sort_with_tracking(lowdin_charge(:,i),basis_idx(:,i),nbasis)

    end do GetLowdinCharges
    call mem_dealloc(ShalfC)


    ! Distances between occ and virt orbitals
    call mem_alloc(DistOccUnocc,nocc,nvirt)
    call general_distance_table(nocc,nvirt,MyMolecule%carmomocc,MyMolecule%carmomvirt,DistOccUnocc)
    ! .. and between occ and occ orbitals
    call mem_alloc(DistOccOcc,nocc,nocc)
    call general_distance_table(nocc,nocc,MyMolecule%carmomocc,MyMolecule%carmomocc,DistOccOcc)
    call mem_alloc(sorted_dists,nocc)

    ! *************************************
    ! Orbital assignment and orbital extent
    ! *************************************
    OrbitalLoop: do i=1,nMO

       charge = 0E0_realk

       ! Add basis functions until sum of Lowdin charges is close enough to 1
       call mem_alloc(basToatom,nbasis)
       do atom=1,natoms
          do nu=MyMolecule%atom_start(atom),MyMolecule%atom_end(atom) ! nu \in atom
             basToatom(nu)=atom
          enddo
       enddo
       call mem_alloc(WhichAOs,nbasis)
       WhichAOs = .FALSE.
       LowdinAddLoop: do n=1,nbasis
          IF(.NOT.WhichAOs(basis_idx(n,i)))THEN !if not already included
             charge = charge + lowdin_charge(n,i)
             error = 1E0_realk - charge
             WhichAOs(basis_idx(n,i)) = .TRUE.
             if(error < approximated_norm_threshold .or. COUNT(WhichAOs)==nbasis ) then 
                ! atom list converged for orbital i
                exit LowdinAddLoop
             elseif(error .GT. 6.0E0_realk*approximated_norm_threshold) then 
                !include all basis functions attached to this atom
                atom = basToatom(basis_idx(n,i))
                do nu=MyMolecule%atom_start(atom),MyMolecule%atom_end(atom)
                   IF(.NOT.WhichAOs(nu))THEN !if not already included 
                      WhichAOs(nu) = .TRUE.
                      !modify the charge accordingly
                      kloop25: do k=1,nbasis
                         IF(basis_idx(k,i).EQ.nu)THEN 
                            charge = charge + lowdin_charge(k,i)
                            exit kloop25
                         ENDIF
                      enddo kloop25
                   ENDIF
                enddo
             end if
          ENDIF
       end do LowdinAddLoop
       !SANITY CHECK Include full Shell (all 3 P x,y,z components) if one part of shell have been included
       do bas=1,nbasis
          IF(WhichAOs(basis_idx(bas,i)))THEN
             DO J=MyMolecule%bas_Start(basis_idx(bas,i)),MyMolecule%bas_End(basis_idx(bas,i))
                WhichAOs(J) = .TRUE.
             ENDDO
          ENDIF
       enddo
       ! Set list of aos to consider
       norbital_extent = COUNT(WhichAOs)
       call mem_alloc(list_of_aos_to_consider,norbital_extent)
       II = 1
       do IDX=1,nbasis
          IF(WhichAOs(IDX))THEN
             list_of_aos_to_consider(II) = IDX
             II = II + 1
          ENDIF
       enddo
       IF(II-1.NE.norbital_extent)call lsquit('DEC error in AO shell sanitycheck',-1)
       call mem_dealloc(WhichAOs)
       call mem_dealloc(basToatom)

       ! Print orbital info for high print levels
       ! ****************************************
       if(DECinfo%PL>1) then
          write(DECinfo%output,'(1X,a,i10)') 'ORBITAL: ', i
          write(DECinfo%output,*) '-------------------------------------'
          write(DECinfo%output,'(1X,a,100i5)')    'BASISFUNCTIONS  : ', list_of_aos_to_consider
          write(DECinfo%output,'(1X,a,100f10.3)') 'LOWDIN          : ', lowdin_charge(1:norbital_extent,i)
          write(DECinfo%output,*)
       end if


       ! -- Create orbital
       ! *****************

       if(  (i>offset)  .and.  (i<=nocc) ) then   ! Valence orbital

          ! Central orbital is the occupied orbital itself
          ! ----------------------------------------------
          central_orbital = i
          OccOrbitals(i) = orbital_init(i,central_orbital, &
               norbital_extent,list_of_aos_to_consider)

       else  ! virtupied orbital or core orbital

          ! Sort occ orbitals according to distance to orbital "i"
          if(i.le.offset) then  ! core orbital
             sorted_dists = DistOccOcc(:,i)
          else ! virt orbital
             sorted_dists = DistOccUnocc(:,i-nocc)
          end if
          call real_inv_sort_with_tracking(sorted_dists,sorted_orbitals,nocc)

          ! Assign to nearest valence orbital (ensure that we do not assign to core orbitals)
          Assigning: do n=nocc,1,-1
             if(sorted_orbitals(n)>offset) then
                central_orbital = sorted_orbitals(n)
                exit Assigning
             end if
          end do Assigning

          if(DECinfo%PL>1) then
             write(DECinfo%output,'(1X,a,i10)') 'Sorted (occ orbs,dist) for orbital: ', i
             write(DECinfo%output,*) '--------------------------------------------------------'
             do j=1,nocc
                write(DECinfo%output,*) sorted_orbitals(j), sorted_dists(j)
             end do
             write(DECinfo%output,*) 'Central orbital: ', central_orbital
          end if

          if(i.le.offset) then  ! core orbital
             OccOrbitals(i) = orbital_init(i,central_orbital, &
                  norbital_extent,list_of_aos_to_consider)
          else ! Unocc orbital
             UnoccOrbitals(i-nocc) = orbital_init(i,central_orbital, &
                  norbital_extent,list_of_aos_to_consider)
          end if

       end if

       call mem_dealloc(list_of_aos_to_consider)

    end do OrbitalLoop

    ! Sanity check
    call DECCO_assignment_sanity_check(MyMolecule,DistOccUnocc,OccOrbitals,UnoccOrbitals)

    ! Print orbital info
    call print_orbital_info_DECCO(nocc,nvirt,OccOrbitals,&
         & UnoccOrbitals,MyMolecule%ncore)

    call mem_dealloc(sorted_dists)
    call mem_dealloc(DistOccUnocc)
    call mem_dealloc(DistOccOcc)
    call mem_dealloc(lowdin_charge)
    call mem_dealloc(basis_idx)

    call LSTIMER('GenerateOrb',tcpu,twall,DECinfo%output)

  end subroutine GenerateOrbitals_DECCO




  !> Sanity check for orbital assignment in DECCO
  !> Includes reassigning of virt orbitals, if necessary
  subroutine DECCO_assignment_sanity_check(MyMolecule,DistOccUnocc,OccOrbitals,UnoccOrbitals)

    implicit none
    !> Molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Distances between occ and virt orbitals
    real(realk),intent(in) :: DistOccUnocc(MyMolecule%nocc,MyMolecule%nvirt)
    !> Occupied orbitals
    type(decorbital), intent(in), dimension(MyMolecule%nocc) :: OccOrbitals
    !> Unoccupied orbitals
    type(decorbital), intent(inout), dimension(MyMolecule%nvirt) :: UnoccOrbitals
    integer :: i, nvirt_per_occ(MyMolecule%nocc),a,sorted_orbitals(MyMolecule%nvirt),idx,ax
    real(realk) :: sorted_dists(MyMolecule%nvirt)
    real(realk),parameter :: maxdist = 2.0/bohr_to_angstrom
    logical :: did_reassign


    call DECCO_get_nvirt_per_occ_orbital(MyMolecule%nocc,MyMolecule%nvirt,&
         & UnoccOrbitals,nvirt_per_occ)

    ! Else we check that each valence orbital has at least on virt orbital assigned
    ! and that core orbitals have nothing assigned

    ! Valence check
    ValenceCheck: do i=MyMolecule%ncore+1,MyMolecule%nocc

       ValenceProblem: if(nvirt_per_occ(i)==0) then

          ! No virt orbitals assigned to orbital "i"
          ! --> we try to fix the problem by reassigning...

          ! Sort virt orbitals according to distance from orbital "i"
          do a=1,MyMolecule%nvirt
             sorted_dists(a) = DistOccUnocc(i,a)
          end do
          call real_inv_sort_with_tracking(sorted_dists,sorted_orbitals,MyMolecule%nvirt)

          ! Reassigning procedure
          ! *********************
          did_reassign=.false.

          ! Loop over virt orbitals, start with those closest to "i"
          ReAssign: do a=MyMolecule%nvirt,1,-1

             ! Occupied orbital to which virt orbital "a" is currently assigned
             idx = UnoccOrbitals(sorted_orbitals(a))%centralatom

             ! Unsorted virt orbital index
             ax = sorted_orbitals(a)

             ! Reassign virt orbital "a" to occ orbital "i" IF:
             ! 
             ! 1. The occ orbital "idx" to which "a" is currently assigned 
             !    has more than one virt orbital assigned,
             !    such that orbital "idx" does not get into the same
             !    problem that orbital "i" is currently facing.
             !
             !    AND
             ! 
             ! 2. The distance from "a" to "i" is smaller than 2 Angstrom
             ! 
             ! 
             if(nvirt_per_occ(idx)>1 .and. sorted_dists(a)<maxdist) then
                UnoccOrbitals(ax)%centralatom = i
                nvirt_per_occ(idx) = nvirt_per_occ(idx) - 1
                did_reassign=.true.
                exit ReAssign
             end if
          end do ReAssign

          if(did_reassign) then
             write(DECinfo%output,*) 'WARNING: Reassigning virtupied orbital ', ax, &
                  & ' from occ orbital ', idx, ' to occ orbital ',i
          else
             print *, 'Error for orbital: ',i
             print *, 'You cannot use DECCO for this system - unless you are willing to&
                  & change the source code!'
             print *, 'nvirt_per_occ: ', nvirt_per_occ

             write(DECinfo%output,*) 'Print DECCO orbital info before quitting...'
             call print_orbital_info_DECCO(MyMolecule%nocc,MyMolecule%nvirt,OccOrbitals,&
                  & UnoccOrbitals,MyMolecule%ncore)
             call lsquit('DECCO_assignment_sanity_check: &
                  &Valence orbital has no virtupied orbital(s) assigned!',-1)
          end if

       end if ValenceProblem

    end do ValenceCheck

    ! Core check
    do i=1,MyMolecule%ncore
       if(nvirt_per_occ(i)/=0) then
          ! This should never happen!
          ! Something went wrong somewhere...
          print *, 'Error for orbital: ',i
          print *, 'You cannot use DECCO for this system - unless you are willing to&
               & change the source code!'
          print *, 'nvirt_per_occ: ', nvirt_per_occ

          write(DECinfo%output,*) 'Print DECCO orbital info before quitting...'
          call print_orbital_info_DECCO(MyMolecule%nocc,MyMolecule%nvirt,OccOrbitals,&
               & UnoccOrbitals,MyMolecule%ncore)
          call lsquit('DECCO_assignment_sanity_check: &
               &Core orbital has virtupied orbital(s) assigned!',-1)
       end if
    end do


  end subroutine DECCO_assignment_sanity_check


  !> \brief Get matrix [S^{1/2} C] used for Lowdin population analysis
  !> \author Kasper Kristensen
  subroutine Get_matrix_for_lowdin_analysis(MyMolecule, MyLsitem, ShalfC)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> General LS info
    type(lsitem),intent(inout) :: MyLsitem
    !> S^{1/2} C matrix
    real(realk), dimension(MyMolecule%nbasis,MyMolecule%nMO),intent(inout) :: ShalfC
    real(realk), pointer :: Shalf(:,:)
    integer :: nb,no,nv,i,j,k,nMO
    real(realk),pointer :: basis(:,:),S(:,:), Co(:,:), Cv(:,:)

    nb = MyMolecule%nbasis
    no = MyMolecule%nocc
    nv = MyMolecule%nvirt
    nMO = MyMolecule%nMO
    call mem_alloc(basis,nb,nb)

    !basis(1:nb,1:MyMolecule%nocc) = MyMolecule%Co(1:nb,1:MyMolecule%nocc)
    !basis(1:nb,MyMolecule%nocc+1:nb) = MyMolecule%Cv(1:nb,1:MyMolecule%nvirt)
    if( MyMolecule%mem_distributed )then
       call mem_alloc(Co,nb,MyMolecule%nocc)
       call mem_alloc(Cv,nb,MyMolecule%nvirt)

       call tensor_gather(1.0E0_realk,MyMolecule%Co,0.0E0_realk,Co,i8*nb*no)
       call tensor_gather(1.0E0_realk,MyMolecule%Cv,0.0E0_realk,Cv,i8*nb*nv)
    else
       Co => MyMolecule%Co%elm2
       Cv => MyMolecule%Cv%elm2
    endif

    do j=1,MyMolecule%nocc
       do i =1,nb
          basis(i,j) = Co(i,j)
       enddo
    enddo
    k = MyMolecule%nocc+1
    do j = 1,MyMolecule%nvirt
       do i =1,nb
          basis(i,k) = Cv(i,j)
       enddo
       k=k+1
    enddo

    if( MyMolecule%mem_distributed )then
       call mem_dealloc(Co)
       call mem_dealloc(Cv)
    endif

    Co => null()
    Cv => null()

    ! AO overlap
    call mem_alloc(S,nb,nb)
    call II_get_mixed_overlap_full(DECinfo%output,DECinfo%output,MyLsitem%SETTING,&
         & S,nb,nb,AORdefault,AORdefault)

    ! Get S^{1/2} matrix
    ! ******************
    call mem_alloc(Shalf,nb,nb)
    call get_power_of_symmetric_matrix(nb,0.5E0_realk,S,Shalf)
    call mem_dealloc(S)

    ! S^{1/2} C
    ! *********
    call dec_simple_dgemm(nb,nb,nMO,Shalf,basis,ShalfC,'n','n')

    call mem_dealloc(Shalf)
    call mem_dealloc(basis)


  end subroutine Get_matrix_for_lowdin_analysis



  !> \brief Get Mulliken gross charges for a given orbital on all basisfunctions
  subroutine GetMullikenVector(orbI,nbasis,MyMolecule,charges,CtS)

    implicit none
    integer, intent(in) :: orbI ! orbital number
    integer, intent(in) :: nbasis ! Number of basis functions
    type(fullmolecule), intent(in) :: MyMolecule
    real(realk), dimension(nbasis), intent(inout) :: charges
    real(realk), dimension(nbasis,nbasis), intent(in) :: CtS
#ifdef VAR_PTR_RESHAPE
    real(realk), pointer, contiguous :: C(:,:)
#else
    real(realk), pointer :: C(:,:)
#endif
    integer :: nu
    integer :: nocc,nvirt

    ! Init stuff
    nocc=MyMolecule%nocc
    nvirt=MyMolecule%nvirt
    call mem_alloc(C,nbasis,nbasis)

    if( MyMolecule%mem_distributed )then
       print *,"WARNING(GetMullikenVector): getting matrix in full, this should&
       & be replaced by a PDM operation, plus this is not tested at all"
       call tensor_gather(1.0E0_realk,MyMolecule%Co,0.0E0_realk,C(1:nbasis,1:nocc) ,i8*nbasis*nocc)
       call tensor_gather(1.0E0_realk,MyMolecule%Cv,0.0E0_realk,C(1:nbasis,nocc+1:nbasis),i8*nbasis*nvirt)
    else
       C(1:nbasis,1:nocc) = MyMolecule%Co%elm2(1:nbasis,1:nocc)
       C(1:nbasis,nocc+1:nbasis) = MyMolecule%Cv%elm2(1:nbasis,1:nvirt)
    endif

!    charges=0.0E0_realk

    do nu=1,nbasis
       ! Mulliken charge for orbital I, on the AO "nu":
       !
       ! chargeI(nu) =  [C^T S]_{I,nu} C_{nu,I}
       !
       ! where C are MO coefficients and S is AO overlap matrix,
       ! nu is an AO index and I is a MO index for the orbital in question.
       !
       ! It is seen that:
       ! sum_{all basisfunctions} chargeI(nu) = sum_{all nu} [C^T S]_{I,nu} C_{nu,I}
       !                               = [C^T S C]_{I,I}
       !                               = 1
       ! as it must be.
       charges(nu)=CtS(orbI,nu)*C(nu,orbI)
    end do
    call mem_dealloc(C)

  end subroutine GetMullikenVector

  !> \brief Get Mulliken gross charges for a given orbital on all atoms
  subroutine GetDistanceCentralAtom(orbI,natoms,MyMolecule,central_atom)
    implicit none
    integer, intent(in) :: orbI ! orbital number
    integer, intent(in) :: natoms ! number of atoms
    type(fullmolecule), intent(in) :: MyMolecule
    integer, intent(inout) :: central_atom
    !
    integer :: atom,catom,nocc,nvirt
    real(realk) :: XMO,YMO,ZMO,XATOM,YATOM,ZATOM,XDIST,YDIST,ZDIST
    real(realk) :: SQDIST(nAtoms),SQDISTVAL    
    nocc=MyMolecule%nocc
    ! Init stuff
    IF(orbI.GT.nocc)THEN
       !virtual orbital
       XMO = - MyMolecule%carmomvirt(1,orbI-nocc)
       YMO = - MyMolecule%carmomvirt(2,orbI-nocc)
       ZMO = - MyMolecule%carmomvirt(3,orbI-nocc)
    ELSE
       !occupied orbital
       XMO = - MyMolecule%carmomocc(1,orbI)
       YMO = - MyMolecule%carmomocc(2,orbI)
       ZMO = - MyMolecule%carmomocc(3,orbI)
    endif
    ! Loop over atoms
    do atom=1,natoms
       XATOM = MyMolecule%AtomCenters(1,atom)
       YATOM = MyMolecule%AtomCenters(2,atom)
       ZATOM = MyMolecule%AtomCenters(3,atom)
       XDIST = XATOM + XMO
       YDIST = YATOM + YMO
       ZDIST = ZATOM + ZMO
       SQDIST(atom) = XDIST*XDIST + YDIST*YDIST + ZDIST*ZDIST
    enddo
    SQDISTVAL=HUGE(SQDIST(1))
    catom = -101
    do atom=1,natoms
       IF(SQDIST(atom).LT.SQDISTVAL)THEN
          IF(.NOT.MyMolecule%PhantomAtom(atom))THEN
             SQDISTVAL = SQDIST(atom)
             catom = atom
          ENDIF
       ENDIF
    enddo
    central_atom = catom
  end subroutine GetDistanceCentralAtom


  !> \brief Lowdin population analysis: Get Lowdin charges on all atoms for a given orbital
  !> \author Kasper Kristensen
  !> \date September 2011
  subroutine GetLowdinVector(orb_idx,nbasis,nMO,ShalfC,charges)

    implicit none
    !> Orbital number
    integer, intent(in) :: orb_idx
!    !> Number of atoms
!    integer, intent(in) :: natoms
    !> Number of basis functions
    integer, intent(in) :: nbasis
    !> Number of MOs (usually equal to nbasis but can be different)
    integer,intent(in) :: nMO
!    !> Full molecule info
!    type(fullmolecule), intent(in) :: MyMolecule
    !> Overlap matrix to power 1/2 multiplied by MO coefficients: S^{1/2} C
    real(realk), dimension(nbasis,nMO), intent(in) :: ShalfC
    !> Lowdin charges
    real(realk), dimension(nbasis), intent(inout) :: charges
!    integer, pointer :: atom_start(:) => null()
!    integer, pointer :: atom_end(:) => null()
    integer :: atom,mu


    ! Init stuff
!    charges=0.0E0_realk
!    atom_start => MyMolecule%atom_start
!    atom_end => MyMolecule%atom_end

    do mu=1,nbasis

       ! Lowdin charge for orbital I, on the atom "atom":
       !
       ! chargeI(atom) = sum_{mu \in atom} { [S^{1/2} C]_{mu,I} }^2
       !
       ! where C are MO coefficients, S^{1/2} is AO overlap matrix to the power 1/2,
       ! mu is an AO index and I is a MO index for the orbital in question.
       !
       ! It is seen that:
       ! sum_{all atoms} chargeI(atom) = 1

       charges(mu)= ShalfC(mu,orb_idx)**2
    end do

!    atom_start => null()
!    atom_end => null()

  end subroutine GetLowdinVector



  !> \brief Lowdin population analysis: Get Lowdin charges on all orbitals for a given orbital
  !> \author Kasper Kristensen
  !> \date September 2011
!!$  subroutine GetLowdinVector_orbitals(orb_idx,nbasis,ShalfC,charges)
!!$
!!$    implicit none
!!$    !> Orbital number
!!$    integer, intent(in) :: orb_idx
!!$    !> Number of basis functions
!!$    integer, intent(in) :: nbasis
!!$    !> Overlap matrix to power 1/2 multiplied by MO coefficients: S^{1/2} C
!!$    real(realk), dimension(nbasis,nbasis), intent(in) :: ShalfC
!!$    !> Lowdin charges
!!$    real(realk), dimension(nbasis), intent(inout) :: charges
!!$    integer :: mu
!!$
!!$
!!$    ! Loop over orbitals
!!$    do mu=1,nbasis
!!$
!!$       ! Lowdin charge for AO orbital mu: { [S^{1/2} C]_{mu,I} }^2
!!$       !
!!$       ! where C are MO coefficients, S^{1/2} is AO overlap matrix to the power 1/2,
!!$       ! mu is an AO index and I is a MO index for the orbital in question.
!!$       charges(mu) = ShalfC(mu,orb_idx)**2
!!$
!!$    end do
!!$
!!$
!!$  end subroutine GetLowdinVector_orbitals



  !> \brief Atomic norms for expansion coefficients
!!$  subroutine GetOrbitalAtomicNorm(num,MyMolecule,natoms,atomic_norms)
!!$
!!$    implicit none
!!$    integer, intent(in) :: num
!!$    type(fullmolecule), intent(in) :: MyMolecule
!!$    integer, intent(in) :: natoms
!!$    real(realk), dimension(natoms), intent(inout) :: atomic_norms
!!$
!!$    integer, pointer :: atom_start(:) => null()
!!$    integer, pointer :: atom_end(:) => null()
!!$    real(realk), pointer :: C(:,:)
!!$    integer :: i,nbasis,nocc,nvirt
!!$
!!$    ! Init stuff
!!$    nbasis=MyMolecule%nbasis
!!$    nocc=MyMolecule%nocc
!!$    nvirt=MyMolecule%nvirt
!!$    call mem_alloc(C,nbasis,nbasis)
!!$    C(1:nbasis,1:nocc) = MyMolecule%Co(1:nbasis,1:nocc)
!!$    C(1:nbasis,nocc+1:nbasis) = MyMolecule%Cv(1:nbasis,1:nvirt)
!!$    atom_start => MyMolecule%atom_start
!!$    atom_end => MyMolecule%atom_end
!!$
!!$    atomic_norms=0.0E0_realk
!!$    do i=1,natoms
!!$       atomic_norms(i)=dot_product(C(atom_start(i):atom_end(i),num), &
!!$            C(atom_start(i):atom_end(i),num))
!!$    end do
!!$
!!$    atom_start => null()
!!$    atom_end => null()
!!$    call mem_dealloc(C)
!!$
!!$  end subroutine GetOrbitalAtomicNorm


  !> \brief Write orbital to file
  subroutine orbital_write(orb,iunit)
    implicit none
    type(decorbital), intent(in) :: orb
    integer, intent(in) :: iunit
    integer(kind=long) :: orbitalnumber64,centralatom64,numberofaos64,aos64(orb%numberofaos)

    orbitalnumber64 = orb%orbitalnumber
    centralatom64   = orb%centralatom
    numberofaos64   = orb%numberofaos
    aos64           = orb%aos

    write(iunit) orbitalnumber64
    write(iunit) centralatom64
    write(iunit) numberofaos64
    write(iunit) aos64

  end subroutine orbital_write

  !> \brief Read orbital from file
  function orbital_read(iunit) result(orb)
    implicit none
    type(decorbital) :: orb
    integer, intent(in) :: iunit
    integer :: i
    integer(kind=long) :: orbitalnumber64,centralatom64,numberofaos64
    integer(kind=4) :: orbitalnumber32,centralatom32,numberofaos32
    integer(kind=long),pointer :: aos64(:)
    integer(kind=4),pointer :: aos32(:)

    !ConvertFrom64Bit: if(DECinfo%convert64to32) then
    !files always written in 64 bit integers
       read(iunit) orbitalnumber64
       read(iunit) centralatom64
       read(iunit) numberofaos64
       call mem_alloc(aos64,numberofaos64)
       read(iunit) aos64

       orb%orbitalnumber = int(orbitalnumber64)
       orb%centralatom   = int(centralatom64)
       orb%numberofaos = int(numberofaos64)
       call mem_alloc(orb%aos,orb%numberofaos)
       do i=1,orb%numberofaos
          orb%aos(i) = int(aos64(i))
       end do
       call mem_dealloc(aos64)

    !elseif(DECinfo%convert32to64)then

    !   read(iunit) orbitalnumber32
    !   read(iunit) centralatom32
    !   read(iunit) numberofatoms32
    !   call mem_alloc(atoms32,numberofatoms32)
    !   read(iunit) atoms32
    !   ! Convert 32 bit integers to 32 bit and store in decorbital type
    !   orb%orbitalnumber = orbitalnumber32
    !   orb%centralatom = centralatom32
    !   orb%numberofatoms = numberofatoms32
    !   call mem_alloc(orb%atoms,orb%numberofatoms)
    !   do i=1,orb%numberofatoms
    !      orb%atoms(i) = atoms32(i)
    !   end do
    !   call mem_dealloc(atoms32)

    !else

    !   read(iunit) orb%orbitalnumber
    !   read(iunit) orb%centralatom
    !   read(iunit) orb%numberofatoms
    !   call mem_alloc(orb%atoms,orb%numberofatoms)
    !   read(iunit) orb%atoms

    !end if ConvertFrom64Bit

  end function orbital_read


  !> \brief Check that LCM orbitals are correct by projecting
  !> against canonical orbitals.
  !> \author Kasper Kristensen
  !> \date November 2010
  subroutine check_lcm_against_canonical(MyMolecule,MyLsitem)


    implicit none
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> LSDALTON info
    type(lsitem), intent(inout) :: MyLsitem
    type(matrix) :: Cocc_can, Cocc_lcm, Cvirt_can, Cvirt_lcm
    type(matrix) :: Ccan,F,S, SCLocc, SCLvirt, UnitMatrix, Clcm
    type(matrix) :: tmp1
    real(realk), pointer :: eival(:)
    integer :: nb,no,nv, nstart, nend,nvirt_end,i,j
    integer :: idx,a,b,funit
    real(realk) :: nocc_calc, nvirt_calc, tmp
    real(realk) :: can_norm, lcm_norm
    real(realk),pointer :: occ_vector(:), virt_vector(:)

    ! Sanity check: This is just a debugging routine and it only works for dense matrix type
    if(matrix_type/=mtype_dense) then
       call lsquit('check_lcm_against_canonical: Only implemented for dense matrix type!',-1)
    end if


    ! Initialize stuff
    ! ****************
    nb = MyMolecule%nbasis
    no = MyMolecule%nocc
    nv = MyMolecule%nvirt
    call mat_init(Ccan, nb, nb )
    call mat_init(F,    nb, nb )
    call mat_init(S,    nb, nb )
    call mem_alloc(eival,   nb )


    ! Set Fock and overlap matrices
    ! *****************************
    
    if( MyMolecule%mem_distributed )then
       !if(matrix_type == mtype_dense) then
       !   call tensor_convert(MyMolecule%fock,F%elms)
       !if(matrix_type == mtype_pdmm) then !pdmm_Marray is not known here!!!
       !   call tensor_copy_data(MyMolecule%fock,pdmm_Marray(F%PDMID)%p)
       !else
          call lsquit("ERROR(check_lcm_against_canonical)this does not work with the selected matrix type (yet)",-1)
       !endif
    else
       call mat_set_from_full(MyMolecule%fock%elm2(1:nb,1:nb), 1E0_realk,F)
    endif
    call II_get_overlap(DECinfo%output,DECinfo%output,mylsitem%SETTING,S)


    ! Get canonical orbitals, sorted in order of increasing orbital energies
    ! **********************************************************************
    call mat_diag_f(F,S,eival,Ccan)



    ! Occupied canonical orbitals
    ! ***************************

    ! These are the first nb*no elements in Ccan
    call mat_init(Cocc_can,nb,no)
    nstart=1
    nend = nb*no
    Cocc_can%elms(nstart:nend) = Ccan%elms(nstart:nend)

    ! Unoccupied canonical orbitals
    ! *****************************
    ! These are the last nb*nv elements,
    ! i.e. from element nb*no+1 to nb*nb
    call mat_init(Cvirt_can,nb,nv)
    nstart = nb*no + 1
    nvirt_end = nb*nv
    nend = nb*nb
    Cvirt_can%elms(1:nvirt_end) = Ccan%elms(nstart:nend)



    ! Set LCM orbitals
    ! ****************
    call mat_init(Cocc_lcm,nb,no)
    call mat_init(Cvirt_lcm,nb,nv)
    if( MyMolecule%mem_distributed )then
       if(matrix_type == mtype_dense) then
          call tensor_gather(1.0E0_realk,MyMolecule%Co,0.0E0_realk,Cocc_lcm%elms ,i8*nb*no)
          call tensor_gather(1.0E0_realk,MyMolecule%Cv,0.0E0_realk,Cvirt_lcm%elms,i8*nb*nv)
       else
          call lsquit("ERROR(check_lcm_against_canonical)this does not work with the selected matrix type (yet)",-1)
       endif
    else
       call mat_set_from_full(MyMolecule%Co%elm2(1:nb,1:no), 1E0_realk,Cocc_lcm)
       call mat_set_from_full(MyMolecule%Cv%elm2(1:nb,1:nv), 1E0_realk,Cvirt_lcm)
    endif


    ! Construct canonical/LCM overlap matrix for occupied space
    ! *********************************************************
    ! SCLocc = Cocc_can^T S Cocc_lcm
    call mat_init(SCLocc,no,no)
    call util_AO_to_MO_different_trans(Cocc_can, S, Cocc_lcm, SCLocc)


    ! Construct canonical/LCM overlap matrix for virtual space
    ! ********************************************************
    ! SCLvirt = Cvirt_can^T S Cvirt_lcm
    call mat_init(SCLvirt,nv,nv)
    call util_AO_to_MO_different_trans(Cvirt_can, S, Cvirt_lcm, SCLvirt)


    ! Check orthonormality
    ! ********************

    ! Unit matrices
    call mat_init(UnitMatrix,nb,nb)
    call mat_zero(UnitMatrix)
    do i=1,nb
       idx = get_matrix_position(i,i,nb,nb)
       UnitMatrix%elms(idx) = 1E0_realk
    end do

    ! Overlap for canonical orbitals
    ! ''''''''''''''''''''''''''''''
    call mat_init(tmp1,nb,nb)
    call util_AO_to_MO_different_trans(Ccan, S, Ccan, tmp1)
    ! Subtract unit matrix from occupied-occupied overlap
    call mat_daxpy(-1E0_realk,UnitMatrix,tmp1)
    ! Find norm of overlap minus unit matrix (this should be zero)
    can_norm = mat_sqnorm2(tmp1)
    can_norm = sqrt(can_norm)

    ! Overlap for LCM orbitals
    ! ''''''''''''''''''''''''
    call mat_init(Clcm,nb,nb)
    ! Set occupied orbitals
    nstart=1
    nend = nb*no
    Clcm%elms(nstart:nend) = Cocc_lcm%elms(nstart:nend)
    ! Set virtual orbitals
    nstart = nb*no + 1
    nend = nb*nb
    nvirt_end = nb*nv
    Clcm%elms(nstart:nend) = Cvirt_lcm%elms(1:nvirt_end)
    ! Get LCM overlap
    call util_AO_to_MO_different_trans(Clcm, S, Clcm, tmp1)
    ! Subtract unit matrix from occupied-occupied overlap
    call mat_daxpy(-1E0_realk,UnitMatrix,tmp1)
    ! Find norm of overlap minus unit matrix (this should be zero)
    lcm_norm = mat_sqnorm2(tmp1)
    lcm_norm = sqrt(lcm_norm)



    call mat_free(F)
    call mat_free(Ccan)
    call mat_free(Clcm)
    call mat_free(tmp1)
    call mat_free(UnitMatrix)


    ! Create occupied and virtupied vectors
    ! **************************************

    ! The occupied vector contains each occupied canonical orbital projected
    ! against the full set of occupied LCM orbitals:
    ! occupied_vector(i) = sum_j SCLocc(i,j)**2
    ! and similarly for the virtual vector.
    ! Thus, each element in the occupied_vector and the virtupied_vector
    ! has to equal one. And the sum of elements equals no/nv for
    ! occupied/virtupied vectors.

    call mem_alloc(occ_vector,no)
    call mem_alloc(virt_vector,nv)
    occ_vector=0E0_realk
    virt_vector=0E0_realk

    ! Occupied vector
    ! '''''''''''''''
    do i=1,no
       do j=1,no
          idx = get_matrix_position(i,j,no,no)
          tmp = SCLocc%elms(idx)
          occ_vector(i) = occ_vector(i) + tmp**2
       end do
    end do

    ! Total number of occupied orbitals
    nocc_calc=0E0_realk
    do i=1,no
       nocc_calc = nocc_calc + occ_vector(i)
    end do


    ! Virtual vector
    ! ''''''''''''''
    do a=1,nv
       do b=1,nv
          idx = get_matrix_position(a,b,nv,nv)
          tmp = SCLvirt%elms(idx)
          virt_vector(a) = virt_vector(a) + tmp**2
       end do
    end do

    ! Total number of virtual orbitals
    nvirt_calc=0E0_realk
    do a=1,nv
       nvirt_calc = nvirt_calc + virt_vector(a)
    end do


    ! Print out
    ! *********

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)') '================================================'
    write(DECinfo%output,'(1X,a)') '               TESTING LCM BASIS'
    write(DECinfo%output,'(1X,a)') '================================================'

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(3X,a)') '  Occupied projections  '
    write(DECinfo%output,'(3X,a)') '------------------------'
    do i=1,no
       write(DECinfo%output,'(3X,i8,g18.8)') i, occ_vector(i)
    end do
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(3X,a)') '  Virtual projections   '
    write(DECinfo%output,'(3X,a)') '------------------------'
    do a=1,nv
       write(DECinfo%output,'(3X,i8,g18.8)') a, virt_vector(a)
    end do
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(3X,a,g18.8)') 'Calculated number of occ. orbitals  :', nocc_calc
    write(DECinfo%output,'(3X,a,g18.8)') 'Calculated number of virt. orbitals :', nvirt_calc
    write(DECinfo%output,'(3X,a,g18.8)') 'Exact number of occ. orbitals       :', real(no)
    write(DECinfo%output,'(3X,a,g18.8)') 'Exact number of virt. orbitals      :', real(nv)
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(3X,a)') '     Orthonormality     '
    write(DECinfo%output,'(3X,a)') '------------------------'
    write(DECinfo%output,'(3X,a,g18.8)') &
         & 'Norm of canonical overlap minus unit matrix :', can_norm
    write(DECinfo%output,'(3X,a,g18.8)') &
         & 'Norm of LCM overlap minus unit matrix       :', lcm_norm
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)


    ! Free remaining matrices and vectors
    ! ***********************************
    call mat_free(Cocc_can)
    call mat_free(Cocc_lcm)
    call mat_free(Cvirt_can)
    call mat_free(Cvirt_lcm)
    call mat_free(S)
    call mat_free(SCLocc)
    call mat_free(SCLvirt)
    call mem_dealloc(eival)
    call mem_dealloc(occ_vector)
    call mem_dealloc(virt_vector)

  end subroutine check_lcm_against_canonical


  !> \brief Check the norm of the occupied orbitals on the other Subsystem
  !> \author Thomas Kjaergaard
  !> \date November 2014
  subroutine check_Occupied_SubSystemLocality(MyMolecule,MyLsitem)
    implicit none
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> LSDALTON info
    type(lsitem), intent(inout) :: MyLsitem
    !
    type(decorbital), pointer :: OccOrbitals(:)
    type(decorbital), pointer :: UnoccOrbitals(:)
    real(realk),pointer :: Cocc1(:,:),Cocc2(:,:)
    integer :: nbasis,nocc,nvirt,natoms,noccsub(2),nsize
    integer,pointer :: nOrb(:),centralatom(:)
    integer :: i,isys,iatom,ibasissub(2),nbasisSub(2)
    real(realk), external :: ddot

    nbasis = MyMolecule%nbasis
    nocc = MyMolecule%nocc
    nvirt = MyMolecule%nvirt
    nBasis = MyMolecule%nbasis
    nAtoms = MyMolecule%natoms

    call mem_alloc(OccOrbitals,nOcc)
    call mem_alloc(UnoccOrbitals,nvirt)
    call GenerateOrbitals_driver(MyMolecule,mylsitem,nocc,nvirt,natoms, &
         & OccOrbitals, UnoccOrbitals)
    !I do not need the UnoccOrbitals
    do i=1,nvirt
       call orbital_free(UnoccOrbitals(i))
    end do
    call mem_dealloc(UnoccOrbitals)

    call mem_alloc(CentralAtom,nocc)
    do i=1,nocc
       CentralAtom(i) = OccOrbitals(i)%centralatom
    enddo
    do i=1,nOcc
       call orbital_free(OccOrbitals(i))
    end do
    call mem_dealloc(OccOrbitals)

    !determine number of occupied orbitals assigned to which subsystems
    noccsub = 0 
    do i=1,nocc
       isys = Mymolecule%SubSystemIndex(centralatom(i))
       noccsub(isys) = noccsub(isys) + 1       
    enddo
    !determine the number of AO orbitals for each atom
    call mem_alloc(nOrb,nAtoms)
    do iAtom=1,nAtoms
       nOrb(iAtom) = mylsitem%input%molecule%ATOM(IAtom)%nContOrbREG
    enddo
    !determine the number of AO orbitals in each subsystem
    ibasisSub = 0     
    do iAtom=1,nAtoms
       isys = Mymolecule%SubSystemIndex(iAtom)
       ibasisSub(isys) = ibasisSub(isys) + nOrb(iAtom)
    enddo
    nbasisSub = ibasisSub
    !make 2 seperate MO coefficient matrices
    !MOs assigned to SubSystem 1 and AOs belonging to Subsystem 2
    call mem_alloc(Cocc1,nbasisSub(2),noccsub(1))
    call BuildSubsystemCMO(1,nbasisSub(2),noccsub(1),nocc,nAtoms,nbasis,&
         & nOrb,MyMolecule%Co%elm2,Cocc1,CentralAtom,MyMolecule%SubSystemIndex)      

    !MOs assigned to SubSystem 2 and AOs belonging to Subsystem 1
    call mem_alloc(Cocc2,nbasisSub(1),noccsub(2))
    call BuildSubsystemCMO(2,nbasisSub(1),noccsub(2),nocc,nAtoms,nbasis,&
         & nOrb,MyMolecule%Co%elm2,Cocc2,CentralAtom,MyMolecule%SubSystemIndex)      

    ! Delete orbitals
    call mem_dealloc(nOrb)
    call mem_dealloc(CentralAtom)

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)') '================================================'
    write(DECinfo%output,'(1X,a)') '         Occupied SubSystem Locality            '
    write(DECinfo%output,'(1X,a)') '================================================'

    write(DECinfo%output,'(1X,a,A)') ' The Occupied Orbitals Assigned to SubSystem: ',&
         & mylsitem%input%molecule%SubsystemLabel(1)
    nsize = nbasisSub(2)*noccsub(1)
    write(DECinfo%output,'(1X,A,ES16.8)') ' Norm = ',sqrt(ddot(nsize,Cocc1,1,Cocc1,1))
    write(DECinfo%output,'(1X,A,I5,I5)')  ' Dim  = ',nbasisSub(2),noccsub(1)
    call ls_output(Cocc1,1,nbasisSub(2),1,noccsub(1),nbasisSub(2),noccsub(1),1,DECinfo%output)

    write(DECinfo%output,'(1X,a,A)') ' The Occupied Orbitals Assigned to SubSystem: ',&
         & mylsitem%input%molecule%SubsystemLabel(2)
    nsize = nbasisSub(1)*noccsub(2)
    write(DECinfo%output,'(1X,A,ES16.8)') ' Norm = ',sqrt(ddot(nsize,Cocc1,1,Cocc1,1))
    write(DECinfo%output,'(1X,A,I5,I5)')  ' Dim  = ',nbasisSub(1),noccsub(2)
    call ls_output(Cocc2,1,nbasisSub(1),1,noccsub(2),nbasisSub(1),noccsub(2),1,DECinfo%output)

    WRITE(DECinfo%output,*)'Full CMO '
    call ls_output(MyMolecule%Co%elm2,1,nbasis,1,nocc,nbasis,nocc,1,DECinfo%output)

    call mem_dealloc(Cocc1)
    call mem_dealloc(Cocc2)
  end subroutine check_Occupied_SubSystemLocality

  !> \brief Check the norm of the occupied orbitals on the other Subsystem
  !> \author Thomas Kjaergaard
  !> \date November 2014
  subroutine force_Occupied_SubSystemLocality(MyMolecule,MyLsitem)
    implicit none
    !> Full molecule info
    type(fullmolecule), intent(inout) :: MyMolecule
    !> LSDALTON info
    type(lsitem), intent(inout) :: MyLsitem
    !
    type(decorbital), pointer :: OccOrbitals(:)
    type(decorbital), pointer :: UnoccOrbitals(:)
    real(realk),pointer :: Cocc1(:,:),Cocc2(:,:)
    integer :: nbasis,nocc,nvirt,natoms,noccsub(2),nsize
    integer,pointer :: nOrb(:),centralatom(:)
    integer :: i,isys,iatom,ibasissub(2),nbasisSub(2)
    real(realk), external :: ddot

    nbasis = MyMolecule%nbasis
    nocc = MyMolecule%nocc
    nvirt = MyMolecule%nvirt
    nBasis = MyMolecule%nbasis
    nAtoms = MyMolecule%natoms

    call mem_alloc(OccOrbitals,nOcc)
    call mem_alloc(UnoccOrbitals,nvirt)
    call GenerateOrbitals_driver(MyMolecule,mylsitem,nocc,nvirt,natoms, &
         & OccOrbitals, UnoccOrbitals)
    !I do not need the UnoccOrbitals
    do i=1,nvirt
       call orbital_free(UnoccOrbitals(i))
    end do
    call mem_dealloc(UnoccOrbitals)

    call mem_alloc(CentralAtom,nocc)
    do i=1,nocc
       CentralAtom(i) = OccOrbitals(i)%centralatom
    enddo
    do i=1,nOcc
       call orbital_free(OccOrbitals(i))
    end do
    call mem_dealloc(OccOrbitals)

    !determine number of occupied orbitals assigned to which subsystems
    noccsub = 0 
    do i=1,nocc
       isys = Mymolecule%SubSystemIndex(centralatom(i))
       noccsub(isys) = noccsub(isys) + 1       
    enddo
    !determine the number of AO orbitals for each atom
    call mem_alloc(nOrb,nAtoms)
    do iAtom=1,nAtoms
       nOrb(iAtom) = mylsitem%input%molecule%ATOM(IAtom)%nContOrbREG
    enddo
    !determine the number of AO orbitals in each subsystem
    ibasisSub = 0     
    do iAtom=1,nAtoms
       isys = Mymolecule%SubSystemIndex(iAtom)
       ibasisSub(isys) = ibasisSub(isys) + nOrb(iAtom)
    enddo
    nbasisSub = ibasisSub
    !make 2 seperate MO coefficient matrices
    !MOs assigned to SubSystem 1 and AOs belonging to Subsystem 2

    call ZeroSubsystemCMO(1,nbasisSub(2),noccsub(1),nocc,nAtoms,nbasis,&
         & nOrb,MyMolecule%Co%elm2,CentralAtom,MyMolecule%SubSystemIndex)      

    call ZeroSubsystemCMO(2,nbasisSub(1),noccsub(2),nocc,nAtoms,nbasis,&
         & nOrb,MyMolecule%Co%elm2,CentralAtom,MyMolecule%SubSystemIndex)      

    WRITE(DECinfo%output,*)'CMO after force_Occupied_SubSystemLocality'
    call ls_output(MyMolecule%Co%elm2,1,nbasis,1,nocc,nbasis,nocc,1,DECinfo%output)

    ! Delete orbitals
    call mem_dealloc(nOrb)
    call mem_dealloc(CentralAtom)

  end subroutine force_Occupied_SubSystemLocality

  !suboutine can be used for both subsystems combinations but made (name conventions)
  !for the case where:
  !this occupied orbital have been assigned to subsystem 1
  !we therefore loop over all AO basis functions that belog to 
  !subsystem 2     
  subroutine BuildSubsystemCMO(isysInput,nbas2,nocc1,nocc,nAtoms,nbasis,&
       & nOrb,Co,CoccSub,CentralAtom,SubSystemIndex)      
    implicit none
    integer :: nbas2,nocc1,nocc,isysInput,nAtoms,nbasis
    integer,intent(in) :: nOrb(nAtoms)
    real(realk),intent(in) :: Co(nbasis,nocc)
    real(realk),intent(inout) :: CoccSub(nbas2,nocc1)
    integer,intent(in) :: CentralAtom(nOcc),SubSystemIndex(nAtoms)
    !
    integer :: noccsub,i,iatom,isys
    integer :: ibasis,ibasisSub,iAO,isys2
    noccsub = 0 
    do i=1,nocc
       isys = SubSystemIndex(CentralAtom(i))
       IF(isys.EQ.isysInput)THEN
          noccsub = noccsub + 1       
          ibasis = 0 
          ibasisSub = 0 
          do iAtom=1,nAtoms
             isys2 = SubSystemIndex(iAtom)
             IF(isys2.NE.isys)THEN
                !include this in the Cocc (assuming isysinput=1)
                !this occupied orbital have been assigned to subsystem 1
                !we therefore loop over all AO basis functions that belog to 
                !subsystem 2 
                do iAO = 1,nOrb(iAtom)
                   CoccSub(ibasisSub+iAO,noccsub)=Co(ibasis+iAO,i)
                enddo
                ibasisSub=ibasisSub+nOrb(iAtom)
             ENDIF
             ibasis=ibasis+nOrb(iAtom)
          end do
       ENDIF
    enddo
  end subroutine BuildSubsystemCMO

  !suboutine can be used for both subsystems combinations but made (name conventions)
  !for the case where:
  !this occupied orbital have been assigned to subsystem 1
  !we therefore loop over all AO basis functions that belog to 
  !subsystem 2     
  subroutine ZeroSubsystemCMO(isysInput,nbas2,nocc1,nocc,nAtoms,nbasis,&
       & nOrb,Co,CentralAtom,SubSystemIndex)      
    implicit none
    integer :: nbas2,nocc1,nocc,isysInput,nAtoms,nbasis
    integer,intent(in) :: nOrb(nAtoms)
    real(realk),intent(inout) :: Co(nbasis,nocc)
    integer,intent(in) :: CentralAtom(nOcc),SubSystemIndex(nAtoms)
    !
    integer :: noccsub,i,iatom,isys
    integer :: ibasis,ibasisSub,iAO,isys2
    noccsub = 0 
    do i=1,nocc
       isys = SubSystemIndex(CentralAtom(i))
       IF(isys.EQ.isysInput)THEN
          noccsub = noccsub + 1       
          ibasis = 0 
          ibasisSub = 0 
          do iAtom=1,nAtoms
             isys2 = SubSystemIndex(iAtom)
             IF(isys2.NE.isys)THEN
                !include this in the Cocc (assuming isysinput=1)
                !this occupied orbital have been assigned to subsystem 1
                !we therefore loop over all AO basis functions that belog to 
                !subsystem 2 
                do iAO = 1,nOrb(iAtom)
                   Co(ibasis+iAO,i) = 0.0E0_realk
                enddo
                ibasisSub=ibasisSub+nOrb(iAtom)
             ENDIF
             ibasis=ibasis+nOrb(iAtom)
          end do
       ENDIF
    enddo
  end subroutine ZeroSubsystemCMO

  !> \brief Write all orbitals to file "DECorbitals.info"
  !> with occupied orbitals before virtupied orbitals.
  !> \author Kasper Kristensen
  !> \date January 2011
  subroutine write_DECorbitals_to_file(nocc,nvirt,&
       &OccOrbitals,UnoccOrbitals)

    implicit none
    !> Number of occupied orbitals in full molecule
    integer, intent(in) :: nocc
    !> Number of virtupied orbitals in full molecule
    integer, intent(in) :: nvirt
    !> Occupied orbital info for full molecule
    type(decorbital), dimension(nocc), intent(in) :: OccOrbitals
    !> Unoccupied orbital info for full molecule
    type(decorbital), dimension(nvirt), intent(in) :: UnoccOrbitals
    character(len=16) :: FileName
    integer :: funit,i


    ! Open file "DECorbitals.info"
    ! ****************************
    if(DECinfo%PL>0) write(DECinfo%output,*) 'Writing DEC orbitals to file DECorbitals.info...'
    FileName = "DECorbitals.info"
    funit=-1
    call lsopen(funit,FileName,'REPLACE','UNFORMATTED')


    ! Write occupied orbitals to file
    ! *******************************
    do i=1,nocc
       call orbital_write(OccOrbitals(i),funit)
    end do


    ! Write virtupied orbitals to file
    ! *********************************
    do i=1,nvirt
       call orbital_write(UnoccOrbitals(i),funit)
    end do


    ! Close file
    ! **********
    call lsclose(funit,'KEEP')


  end subroutine write_DECorbitals_to_file



  !> \brief Read all orbitals from file "DECorbitals.info"
  !> with occupied orbitals before virtupied orbitals.
  !> \author Kasper Kristensen
  !> \date January 2011
  subroutine read_DECorbitals_from_file(nocc,nvirt,&
       &OccOrbitals,UnoccOrbitals)

    implicit none
    !> Number of occupied orbitals in full molecule
    integer, intent(in) :: nocc
    !> Number of virtupied orbitals in full molecule
    integer, intent(in) :: nvirt
    !> Occupied orbital info for full molecule
    type(decorbital), dimension(nocc), intent(inout) :: OccOrbitals
    !> Unoccupied orbital info for full molecule
    type(decorbital), dimension(nvirt), intent(inout) :: UnoccOrbitals
    character(len=16) :: FileName
    integer :: funit,i


    ! Open file "DECorbitals.info"
    ! ****************************
    if(DECinfo%PL>0) write(DECinfo%output,*) 'Reading DEC orbitals from file DECorbitals.info...'
    FileName = "DECorbitals.info"

    funit=-1
    call lsopen(funit,FileName,'OLD','UNFORMATTED')


    ! Read occupied orbitals from file
    ! ********************************
    do i=1,nocc
       OccOrbitals(i) = orbital_read(funit)
    end do


    ! Read virtupied orbitals from file
    ! **********************************
    do i=1,nvirt
       UnoccOrbitals(i) = orbital_read(funit)
    end do


    ! Close file
    ! **********
    call lsclose(funit,'KEEP')



  end subroutine read_DECorbitals_from_file



  !> \brief Reassign orbitals such that no hydrogen atoms have orbitals assigned.
  !> This is done by reassigning orbitals originally assigned to an H atom to the nearest neighbour atom.
  !> (Special case: If the nearest neighbour is also H, then no reassignment is made).
  !> \author Kasper Kristensen
  !> \date May 2011
  subroutine reassign_orbitals(norb,Orbitals,natoms,DistanceTable,mylsitem)

    implicit none
    !> Number of orbitals
    integer, intent(in) :: norb
    !> Orbital info for full molecule
    type(decorbital), dimension(norb), intent(inout) :: Orbitals
    !> Number of atoms in molecule
    integer, intent(in) :: natoms
    !> Distance table for atoms in molecule
    real(realk), dimension(natoms,natoms), intent(in) :: DistanceTable
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !real(realk), dimension(natoms,natoms) :: SortedDistTable
    real(realk), pointer :: SortedDistTable(:,:)
    integer, dimension(nAtoms,nAtoms) :: TrackMatrix
    integer :: i,j,centralatom,neighbor,atomnumber,neighbor_atomnumber
    logical :: included, reassign
    real(realk) :: maxdist, bohr

    ! Hardcode max acceptable distance for reassignment to 1.5 Angstrom
    ! (Should perhaps be a keyword when this subroutine is made more general)
    bohr=bohr_to_angstrom
    maxdist = 1.5E0_realk/bohr

    ! Sort atoms according to distance, and keep track of original indices in TrackMatrix
    call mem_alloc(SortedDistTable,natoms,natoms)
    do i = 1,natoms
     do j=1,natoms
     SortedDistTable(i,j)=DistanceTable(i,j)
     enddo
    enddo
    !SortedDistTable(:,:)=DistanceTable(:,:)
    call sort_track(SortedDistTable,TrackMatrix,nAtoms)


    ! Loop over all orbitals
    OrbitalLoop: do i=1,norb

       ! Init
       reassign=.false.

       ! Central atom
       centralatom = Orbitals(i)%centralatom

       ! Atomic number for central atom
       atomnumber = myLsitem%input%molecule%Atom(centralatom)%atomic_number

       ! Reassign if centralatom is hydrogen and IF:
       ! Nearest neighbor is not H (CHECK 2)
       ! Distance to nearest neighbor is smaller than maxdist (CHECK 3)
       Hatom: if(atomnumber ==1) then

          ! Reassign orbital
          reassign=.true.

          ! Neighbour: Second element in Trackmatrix
          neighbor = TrackMatrix(2,centralatom)

          ! Neighbor is not hydrogen
          neighbor_atomnumber = myLsitem%input%molecule%Atom(neighbor)%atomic_number
          CheckNeighbor: if(neighbor_atomnumber == 1) then ! No reassignment
             if(DECinfo%PL>0) then
                write(DECinfo%output,*)
                write(DECinfo%output,*) 'WARNING, No reassignment: Nearest neighbor is H (atom,neighbor)', &
                     & centralatom,neighbor
                write(DECinfo%output,'(1X,a,F12.3)') 'Distance to neighbor (Angstrom)    =', &
                     & SortedDistTable(2,centralatom)*bohr
             end if
             reassign=.false.
          end if CheckNeighbor

          ! CHECK 3: Distance to neighbor is smaller than maximum distance (currently hardcoded to 1.5 Angstrom)
          CheckDist: if(SortedDistTable(2,centralatom) > maxdist) then ! No reassignment
             if(DECinfo%PL>0) then
                write(DECinfo%output,*)
                write(DECinfo%output,*) 'WARNING, No reassignment: Neighbor distance too large (atom,neighbor)',&
                     & centralatom,neighbor
                write(DECinfo%output,*) 'Max acceptable distance (Angstrom) =', maxdist*bohr
                write(DECinfo%output,'(1X,a,F12.3)') 'Distance to neighbor (Angstrom)    =', &
                     & SortedDistTable(2,centralatom)*bohr
             end if
             reassign=.false.
          end if CheckDist


          if(reassign) then ! Reassign atom to nearest neighbor
             Orbitals(i)%centralatom = neighbor
          end if

       end if Hatom

    end do OrbitalLoop
    call mem_dealloc(SortedDistTable)

  end subroutine reassign_orbitals


  !> \brief Reassign orbitals such that no hydrogen atoms have orbitals assigned using a simple 
  !> integer input rather than the decorbital input used in reassign_orbitals.
  !> This is done by reassigning orbitals originally assigned to an H atom to the nearest neighbour atom.
  !> (Special case: If the nearest neighbour is also H, then no reassignment is made).
  !> \author Kasper Kristensen
  !> \date November 2012
  subroutine reassign_orbitals_simple(nbasis,natoms,DistanceTable,mylsitem,centralatoms)

    implicit none
    !> Number of orbitals (nocc+nvirt)
    integer, intent(in) :: nbasis
    !> Number of atoms in molecule
    integer, intent(in) :: natoms
    !> Distance table for atoms in molecule
    real(realk), dimension(natoms,natoms), intent(in) :: DistanceTable
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> Central atom for each orbital
    integer,intent(inout) :: centralatoms(nbasis)
    real(realk), dimension(natoms,natoms) :: SortedDistTable
    integer, dimension(nAtoms,nAtoms) :: TrackMatrix
    integer :: i,centralatom,neighbor,atomnumber,neighbor_atomnumber
    logical :: reassign
    real(realk) :: maxdist, bohr

    ! Hardcode max acceptable distance for reassignment to 1.5 Angstrom
    bohr=bohr_to_angstrom
    maxdist = 1.5E0_realk/bohr

    ! Sort atoms according to distance, and keep track of original indices in TrackMatrix
    SortedDistTable(:,:)=DistanceTable(:,:)
    call sort_track(SortedDistTable,TrackMatrix,nAtoms)


    ! Loop over all orbitals
    OrbitalLoop: do i=1,nbasis

       ! Init
       reassign=.false.

       ! Central atom for orbital "i"
       centralatom = centralatoms(i)

       ! Atomic number for central atom
       atomnumber = myLsitem%input%molecule%Atom(centralatom)%atomic_number

       ! Reassign if centralatom is hydrogen and IF:
       ! Nearest neighbor is not H (CHECK 1)
       ! Distance to nearest neighbor is smaller than maxdist (CHECK 2)
       Hatom: if(atomnumber ==1) then

          ! Reassign orbital
          reassign=.true.

          ! Neighbour: Second element in Trackmatrix
          neighbor = TrackMatrix(2,centralatom)

          ! CHECK 1: Neighbor is not hydrogen
          neighbor_atomnumber = myLsitem%input%molecule%Atom(neighbor)%atomic_number
          if(neighbor_atomnumber == 1) reassign=.false.

          ! CHECK 2: Distance to neighbor is smaller than maximum accepted distance
          if(SortedDistTable(2,centralatom) > maxdist) reassign=.false.

          if(reassign) then ! Reassign atom to nearest neighbor
             centralatoms(i) = neighbor
          end if

       end if Hatom

    end do OrbitalLoop

  end subroutine reassign_orbitals_simple



  !> \brief Reassign orbitals from donor atom to acceptor atom.
  !> \author Kasper Kristensen
  !> \date August 2011
  subroutine reassign_orbitals_from_atom_to_atom(norb,Orbitals,donor,acceptor)

    implicit none
    !> Number of orbitals
    integer, intent(in) :: norb
    !> Orbital info for full molecule (either occupied or virtual)
    type(decorbital), dimension(norb), intent(inout) :: Orbitals
    !> Donor atom
    integer, intent(in) :: donor
    !> Acceptor atom
    integer, intent(in) :: acceptor
    integer :: i,centralatom

    do i=1,norb

       ! Central atom
       centralatom = Orbitals(i)%centralatom

       ! Reassign if the central atom is the donor
       if(centralatom == donor) then
          Orbitals(i)%centralatom = acceptor
       end if

    end do


  end subroutine reassign_orbitals_from_atom_to_atom




  !> \brief Assign orbitals to a few (possibly just one) atom(s) and extend orbital extent
  !> to include the full molecule.
  !> In this way we automatically simulate a full calculation by
  !> going through the usual DEC routines.
  !> \author Kasper Kristensen
  !> \date April 2011
  subroutine adjust_orbitals_for_full_simulation(nocc,nvirt,&
       &OccOrbitals,UnoccOrbitals,natoms,n)

    implicit none
    !> Number of occupied orbitals in full molecule
    integer, intent(in) :: nocc
    !> Number of virtupied orbitals in full molecule
    integer, intent(in) :: nvirt
    !> Occupied orbital info for full molecule
    type(decorbital), dimension(nocc), intent(inout) :: OccOrbitals
    !> Unoccupied orbital info for full molecule
    type(decorbital), dimension(nvirt), intent(inout) :: UnoccOrbitals
    !> Number of atoms in the molecule
    integer, intent(in) :: natoms
    !> Number of atomic sites with orbitals assigned in the simulation (default: 1)
    !> Thus, the first n atoms will get all orbitals evenly assigned.
    integer,intent(in) :: n
    integer :: i,j,orb,nocc_per_atom,nvirt_per_atom,atom,nbasis

    ! Number of occupied/virtual orbitals per atom evenly distributed over n atoms
    nocc_per_atom = floor(real(nocc)/real(n))
    nvirt_per_atom = floor(real(nvirt)/real(n))
    nbasis = nocc + nvirt

    ! Sanity check: This routine only makes sense to call
    ! if the full molecule if included for the fragments.
    if(.not. DECinfo%InclFullMolecule) then
       call lsquit('adjust_orbitals_for_full_simulation: &
            & DECinfo%InclFullMolecule must be true when this routine is called!',DECinfo%output)
    end if
    ! Check that  n input is meaningful
    if(n > natoms) then
       write(DECinfo%output,*) 'n, natoms', n,natoms
       call lsquit('Input n is larger than number of atoms!', DECinfo%output)
    end if


    ! Assign occupied orbitals
    ! ************************

    orb=0
    do atom=1,n  ! loop over the n first atoms
       do i=1,nocc_per_atom  ! loop over nocc_per_atom orbitals per atom

          ! Increase orbital counter
          orb=orb+1

          ! Free existing orbital info
          call orbital_free(OccOrbitals(orb))
          call mem_alloc(OccOrbitals(orb)%aos,nbasis)

          ! Orbital extent is the full molecule
          OccOrbitals(orb)%numberofaos=nbasis
          do j=1,nbasis
             OccOrbitals(orb)%aos(j) = j
          end do

          ! Central atom for this occupied orbital
          OccOrbitals(orb)%centralatom = atom

       end do
    end do
    i=orb+1

    ! Assign remaining occupied orbitals (if any) to the last atom n
    do orb=i,nocc
       call orbital_free(OccOrbitals(orb))
       call mem_alloc(OccOrbitals(orb)%aos,nbasis)
       OccOrbitals(orb)%numberofaos=nbasis
       do j=1,nbasis
          OccOrbitals(orb)%aos(j) = j
       end do
       OccOrbitals(orb)%centralatom = n
    end do


    ! Assign virtupied orbitals
    ! **************************

    orb=0
    do atom=1,n  ! loop over the n first atoms
       do i=1,nvirt_per_atom  ! loop over nvirt_per_atom orbitals per atom

          ! Increase orbital counter
          orb=orb+1

          ! Free existing orbital info
          call orbital_free(UnoccOrbitals(orb))
          call mem_alloc(UnoccOrbitals(orb)%aos,nbasis)

          ! Orbital extent is the full molecule
          UnoccOrbitals(orb)%numberofaos=nbasis
          do j=1,nbasis
             UnoccOrbitals(orb)%aos(j) = j
          end do

          ! Central atom for this virtupied orbital
          UnoccOrbitals(orb)%centralatom = atom

       end do
    end do
    i=orb+1

    ! Assign remaining virtupied orbitals (if any) to the last atom n
    do orb=i,nvirt
       call orbital_free(UnoccOrbitals(orb))
       call mem_alloc(UnoccOrbitals(orb)%aos,nbasis)
       UnoccOrbitals(orb)%numberofaos=nbasis
       do j=1,nbasis
          UnoccOrbitals(orb)%aos(j) = j
       end do
       UnoccOrbitals(orb)%centralatom = n
    end do


  end subroutine adjust_orbitals_for_full_simulation


  !> \brief Count number of orbitals assigned to each atom.
  !> \author Kasper Kristensen
  !> \date October 2010
  function get_number_of_orbitals_per_atom(Orbitals,norb,natoms,mainass,offset) &
       & result(norb_per_atom)

    implicit none
    !> Total number of orbitals in Orbitals vector
    integer, intent(in) :: norb
    !> Number of atoms in the molecule
    integer, intent(in) :: natoms
    !> Orbital vector (may either be occupied or virtual orbitals)
    type(decorbital), dimension(norb), intent(in) :: Orbitals
    !> Count orbital based on main assignment (true) or secondary assignment (false)
    logical,intent(in) :: mainass
    !> Not consider orbitals with indices 1:offset (default: consider all orbitals)
    integer,intent(in),optional :: offset
    !> Number of orbitals assigned to each atom
    integer, dimension(natoms) :: norb_per_atom
    integer :: j,MyAtom,startidx

    norb_per_atom=0

    if(present(offset)) then
       startidx = offset+1  ! not consider orbitals 1:offset
    else
       startidx=1 ! consider all orbitals
    end if

    do j=startidx,norb
          if(mainass) then
             ! main assignment
             MyAtom=Orbitals(j)%centralatom
          else ! secondary assigment
             MyAtom=Orbitals(j)%secondaryatom
          end if
          norb_per_atom(MyAtom) = norb_per_atom(MyAtom) +1
    end do

  end function get_number_of_orbitals_per_atom



  !> \brief For the set of ntot integers in tot_idx, this function contructs
  !> a logical vector of length ntot where each entry is
  !> .true. if one of the integers in target_idx is contained in tot_idx (target index)
  !> .false. if none of the integers in target_idx are contained in tot_idx (buffer index).
  !> \author Kasper Kristensen
  !> \param ntarget Number of target indices
  !> \param ntot Total number of indices to compare with
  !> \param target_idx List of target indices
  !> \param tot_idx List of target indices
  !> \return is_in_target True if index in total list i contained in target list, false otherwise.
  function which_orbitals_are_target_orbitals(ntarget,ntot,target_idx,tot_idx) &
       result(is_in_target)

    implicit none
    logical, dimension(ntot) :: is_in_target
    integer, intent(in) :: ntarget, ntot
    integer, dimension(ntarget), intent(in) :: target_idx
    integer, dimension(ntot), intent(in) :: tot_idx
    integer :: i,j

    is_in_target(:) = .false.

    do i=1,ntot

       ! Check if the index in the total list (tot_idx) equals one
       ! of the indices in the target list.
       do j=1,ntarget
          if(tot_idx(i) == target_idx(j)) then
             is_in_target(i) = .true.
          end if
       end do

    end do

  end function which_orbitals_are_target_orbitals


  !> \brief Print number of occupied and virtupied orbitals assigned to each atom for
  !> all atoms in the molecule.
  !> \author Kasper Kristensen
  subroutine print_orbital_info(mylsitem,nocc,natoms,nvirt,OccOrbitals,&
       & UnoccOrbitals,ncore)

    implicit none
    !> Dalton LSITEM (just for printing atom type)
    type(lsitem), intent(inout) :: mylsitem
    !> Number of atoms
    integer, intent(in) :: natoms
    !> Number of occupied orbitals
    integer, intent(in) :: nocc
    !> Number of virtupied orbitals
    integer, intent(in) :: nvirt
    !> Occupied orbitals
    type(decorbital), intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals
    type(decorbital), intent(in) :: UnoccOrbitals(nvirt)
    !> Number of of core orbitals. If present and frozen core approx is used,
    !> the first ncore orbitals will not be printed.
    integer,intent(in),optional :: ncore
    integer :: nocc_per_atom(natoms), nvirt_per_atom(natoms)
    integer :: SECnocc_per_atom(natoms), SECnvirt_per_atom(natoms)
    integer :: i, occ_max_orbital_extent, virt_max_orbital_extent, occ_idx, virt_idx,j,offset
    real(realk) :: occ_av_orbital_extent, virt_av_orbital_extent

    if(present(ncore) .and. DECinfo%frozencore) then
       offset=ncore
    else
       offset=0
    end if

    ! Find average and maximum orbital extent
    ! ***************************************

    ! Occupied
    call get_orbital_extent_info(nocc,OccOrbitals,occ_max_orbital_extent,&
         & occ_av_orbital_extent, occ_idx)

    ! Unoccupied
    call get_orbital_extent_info(nvirt,UnoccOrbitals,virt_max_orbital_extent,&
         & virt_av_orbital_extent, virt_idx)


    ! Number of orbitals per atom
    ! ***************************

    ! Occupied
    ! --------
    ! Main assigning
    nocc_per_atom =  get_number_of_orbitals_per_atom(OccOrbitals,nocc,natoms,.true.,offset=offset)
    ! Secondary assigning
    SECnocc_per_atom =  get_number_of_orbitals_per_atom(OccOrbitals,nocc,natoms,.false.,offset=offset)

    ! Unoccupied
    ! ----------
    ! Main assigninig
    nvirt_per_atom =  get_number_of_orbitals_per_atom(UnoccOrbitals,nvirt,natoms,.true.)
    ! Secondary assigning
    SECnvirt_per_atom =  get_number_of_orbitals_per_atom(UnoccOrbitals,nvirt,natoms,.false.)

    ! Print out number of orbitals
    ! ****************************
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) 'ORBITAL DISTRIBUTION INFORMATION'
    write(DECinfo%output,*) '********************************'
    write(DECinfo%output,*)
    if(DECinfo%frozencore) then
       write(DECinfo%output,*) '--- using frozen core approximation, only valence orbitals are printed'
       write(DECinfo%output,*)
    end if
    write(DECinfo%output,*) '   Atom type     Occ Orbitals     Unocc Orbitals     Occ Secondary    Unocc Secondary'
    do i=1,natoms
       write(DECinfo%output,'(1X,I5,4X,A4,4X,I6,11X,I6,11X,I6,11X,i6)') i, MyLsitem%input%molecule%atom(i)%name,&
            & nocc_per_atom(i), nvirt_per_atom(i), SECnocc_per_atom(i), SECnvirt_per_atom(i)
    end do
    write(DECinfo%output,'(1X,A,11X,I6,11X,I6)') 'Total:', nocc-offset, nvirt
    write(DECinfo%output,*)

    ! Print out orbital extent summary
    ! ********************************
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) 'ORBITAL EXTENT INFORMATION (NUMBER OF ATOMS USED TO SPAN EACH ORBITAL)'
    write(DECinfo%output,*) '**********************************************************************'
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a,i6,a,i6)') 'Maximum occ orbital extent (#AOs) = ', occ_max_orbital_extent, &
         & '   -- Orbital index', occ_idx
    write(DECinfo%output,'(1X,a,i6,a,i6)') 'Maximum virt orbital extent (#AOs) = ', virt_max_orbital_extent, &
         & '   -- Orbital index', virt_idx
    write(DECinfo%output,'(1X,a,f12.4)') 'Average occ orbital extent (#AOs)  = ', occ_av_orbital_extent
    write(DECinfo%output,'(1X,a,f12.4)') 'Average virt orbital extent (#AOs) = ', virt_av_orbital_extent
    write(DECinfo%output,*)
    write(DECinfo%output,*)

    if(DECinfo%PL>0) then ! print specific info for each orbital
       do i=1,nocc
          write(DECinfo%output,*) 'Occupied orbital: ', i
          write(DECinfo%output,*) '************************************************'
          write(DECinfo%output,*) 'Central atom: ', OccOrbitals(i)%centralatom
          write(DECinfo%output,*) '# AOs in orbital extent: ', OccOrbitals(i)%numberofaos
          write(DECinfo%output,*) 'AOs in orbital extent  : '
          do j=1,OccOrbitals(i)%numberofaos
             write(DECinfo%output,*) OccOrbitals(i)%aos(j)
          end do
          write(DECinfo%output,*)
       end do
       write(DECinfo%output,*)
       write(DECinfo%output,*)

       do i=1,nvirt
          write(DECinfo%output,*) 'Virtual orbital: ', i
          write(DECinfo%output,*) '************************************************'
          write(DECinfo%output,*) 'Central atom: ', UnOccOrbitals(i)%centralatom
          write(DECinfo%output,*) '# AOs in orbital extent: ', UnOccOrbitals(i)%numberofaos
          write(DECinfo%output,*) 'AOs in orbital extent  : ' 
          do j=1,UnOccOrbitals(i)%numberofaos
             write(DECinfo%output,*) UnOccOrbitals(i)%aos(j)
          end do
          write(DECinfo%output,*)
       end do
       write(DECinfo%output,*)
       write(DECinfo%output,*)

    end if


  end subroutine print_orbital_info






  !> \brief Print number of occupied and virtupied orbitals assigned to each atom for
  !> all atoms in the molecule.
  !> \author Kasper Kristensen
  subroutine print_orbital_info_DECCO(nocc,nvirt,OccOrbitals,&
       & UnoccOrbitals,ncore)

    implicit none
    !> Number of occupied orbitals
    integer, intent(in) :: nocc
    !> Number of virtupied orbitals
    integer, intent(in) :: nvirt
    !> Occupied orbitals
    type(decorbital), intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals
    type(decorbital), intent(in) :: UnoccOrbitals(nvirt)
    !> Number of of core orbitals. If present and frozen core approx is used,
    !> the first ncore orbitals will not be printed.
    integer,intent(in) :: ncore
    integer :: nvirt_per_occ(nocc)
    integer :: i, occ_max_orbital_extent, virt_max_orbital_extent, occ_idx, virt_idx,j
    real(realk) :: occ_av_orbital_extent, virt_av_orbital_extent


    ! Find average and maximum orbital extent
    ! ***************************************
    ! Occupied
    call get_orbital_extent_info(nocc,OccOrbitals,occ_max_orbital_extent,&
         & occ_av_orbital_extent, occ_idx)

    ! Unoccupied
    call get_orbital_extent_info(nvirt,UnoccOrbitals,virt_max_orbital_extent,&
         & virt_av_orbital_extent, virt_idx)

    ! Number of virtupied orbitals assigned to each occupied orbital
    call DECCO_get_nvirt_per_occ_orbital(nocc,nvirt,UnoccOrbitals,nvirt_per_occ)


    ! Print out number of orbitals
    ! ****************************
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '**************************************'
    write(DECinfo%output,*) 'DECCO ORBITAL DISTRIBUTION INFORMATION'
    write(DECinfo%output,*) '**************************************'
    write(DECinfo%output,*)
    if(DECinfo%frozencore) then
       write(DECinfo%output,*) '--- one fragment per valence orbital!'
       write(DECinfo%output,*) '--- core orbitals are absorbed into valence fragments.'
       write(DECinfo%output,*)
    end if
    write(DECinfo%output,*) '   Occ index     #Unocc Orbitals'
    do i=ncore+1,nocc
       write(DECinfo%output,'(1X,I5,12X,I10)') i,nvirt_per_occ(i)
    end do
    write(DECinfo%output,'(1X,A,11X,I6,11X,I6)') 'Total:', nocc-ncore, nvirt
    write(DECinfo%output,*)


    ! Print out orbital extent summary
    ! ********************************
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) 'ORBITAL EXTENT INFORMATION (NUMBER OF ATOMS USED TO SPAN EACH ORBITAL)'
    write(DECinfo%output,*) '**********************************************************************'
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a,i6,a,i6)') 'Maximum occ orbital extent (#AOs) = ', occ_max_orbital_extent, &
         & '   -- Orbital index', occ_idx
    write(DECinfo%output,'(1X,a,i6,a,i6)') 'Maximum virt orbital extent (#AOs) = ', virt_max_orbital_extent, &
         & '   -- Orbital index', virt_idx
    write(DECinfo%output,'(1X,a,f12.4)') 'Average occ orbital extent (#AOs) = ', occ_av_orbital_extent
    write(DECinfo%output,'(1X,a,f12.4)') 'Average virt orbital extent (#AOs) = ', virt_av_orbital_extent
    write(DECinfo%output,*)
    write(DECinfo%output,*)


    ! print specific info for each orbital
    if(DECinfo%PL>0) then 
       do i=1,nocc
          write(DECinfo%output,*) 'Occupied orbital: ', i
          write(DECinfo%output,*) '************************************************'
          write(DECinfo%output,*) '# AOs in orbital extent: ', OccOrbitals(i)%numberofaos
          write(DECinfo%output,*) 'AOs in orbital extent  : '
          do j=1,OccOrbitals(i)%numberofaos
             write(DECinfo%output,*) OccOrbitals(i)%aos(j)
          end do
          write(DECinfo%output,*)
       end do
       write(DECinfo%output,*)
       write(DECinfo%output,*)

       do i=1,nvirt
          write(DECinfo%output,*) 'Virtual orbital: ', i
          write(DECinfo%output,*) '************************************************'
          write(DECinfo%output,*) '# AOs in orbital extent: ', UnOccOrbitals(i)%numberofaos
          write(DECinfo%output,*) 'AOs in orbital extent  : ' 
          do j=1,UnOccOrbitals(i)%numberofaos
             write(DECinfo%output,*) UnOccOrbitals(i)%aos(j)
          end do
          write(DECinfo%output,*)
       end do
       write(DECinfo%output,*)
       write(DECinfo%output,*)

    end if


  end subroutine print_orbital_info_DECCO


  !> Get number of virtual orbital assigned to each occupied orbital
  subroutine DECCO_get_nvirt_per_occ_orbital(nocc,nvirt,UnoccOrbitals,nvirt_per_occ)
    implicit none
    !> Number of occupied orbitals
    integer, intent(in) :: nocc
    !> Number of virtupied orbitals
    integer, intent(in) :: nvirt
    !> Unoccupied orbitals
    type(decorbital), intent(in) :: UnoccOrbitals(nvirt)
    !> Number of virtupied orbitals assigned to each occupied orbital
    integer, intent(inout) :: nvirt_per_occ(nocc)
    integer :: i,j

    nvirt_per_occ=0
    do j=1,nvirt
       nvirt_per_occ(UnoccOrbitals(j)%centralatom) = nvirt_per_occ(UnoccOrbitals(j)%centralatom) + 1
    end do


  end subroutine DECCO_get_nvirt_per_occ_orbital

  !> \brief The DEC scheme only works if for each atom either zero occupied AND zero virtual
  !> orbitals are assigned - or if nonzero occupied AND nonzero virtual orbitals are assigned.
  !> If this is not the case, the system under consideration is presumably a debug molecule
  !> and we quit here, rather than encountering uninitialized pointers later on...
  !> \author Kasper Kristensen 
  !> \date December 2011
  subroutine dec_orbital_sanity_check(natoms,nocc,nvirt,OccOrbitals,&
       & UnoccOrbitals,MyMolecule)

    implicit none
    !> Number of atoms
    integer, intent(in) :: natoms
    !> Number of occupied orbitals
    integer, intent(in) :: nocc
    !> Number of virtupied orbitals
    integer, intent(in) :: nvirt
    !> Occupied orbitals
    type(decorbital), intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals
    type(decorbital), intent(in) :: UnoccOrbitals(nvirt)
    !> Full molecule info
    type(fullmolecule),intent(in) :: MyMolecule
    integer :: nocc_per_atom(natoms), nvirt_per_atom(natoms)
    integer :: i,nfrags,offset
    logical :: something_wrong


    if(DECinfo%frozencore) then
       offset=MyMolecule%ncore
    else
       offset=0
    end if

    ! Number of orbitals per atom
    ! ***************************

    ! Occupied
    nocc_per_atom =  get_number_of_orbitals_per_atom(OccOrbitals,nocc,natoms,.true.,&
         & offset=offset)

    ! Unoccupied
    nvirt_per_atom =  get_number_of_orbitals_per_atom(UnoccOrbitals,nvirt,natoms,.true.)


    something_wrong=.false.
    nfrags=0
    do i=1,natoms

       if(.not.(DECinfo%onlyoccpart.or.DECinfo%onlyvirtpart)&
            & .and.(.NOT.decinfo%PureHydrogendebug))then
          if( (nocc_per_atom(i) == 0) .and. (nvirt_per_atom(i)/=0) ) something_wrong=.true.
          if( (nocc_per_atom(i) /= 0) .and. (nvirt_per_atom(i)==0) ) something_wrong=.true.
          if(something_wrong) then
             write(DECinfo%output,*) 'Atom = ',i
             write(DECinfo%output,*) 'Number of occupied orbitals   assigned = ', nocc_per_atom(i)
             write(DECinfo%output,*) 'Number of virtupied orbitals assigned = ', nvirt_per_atom(i)
             write(DECinfo%output,*) 'If you use Phantom Atoms try .ONLYOCCPART keyword'
             print*,'If you use Phantom Atoms try .ONLYOCCPART keyword'
             call lsquit('Orbital assigment is inconsistent &
                  & with DEC scheme',DECinfo%output)
          end if
       end if

       ! Count number of atomic fragments
       if( (nocc_per_atom(i) /= 0) .and. (nvirt_per_atom(i)/=0) ) then
          nfrags=nfrags+1
       end if
    end do


    ! Secondary assignment check
    ! **************************
    nocc_per_atom =  get_number_of_orbitals_per_atom(OccOrbitals,nocc,natoms,.false.,&
         & offset=offset)
    nvirt_per_atom =  get_number_of_orbitals_per_atom(UnoccOrbitals,nvirt,natoms,.false.)
    something_wrong=.false.
    do i=1,natoms
       if( (nocc_per_atom(i) == 0) .and. (nvirt_per_atom(i)/=0) ) something_wrong=.true.
       if( (nocc_per_atom(i) /= 0) .and. (nvirt_per_atom(i)==0) ) something_wrong=.true.
       if(something_wrong) then
          write(DECinfo%output,*) 'Atom = ',i
          write(DECinfo%output,*) 'Number of occupied orbitals   assigned = ', nocc_per_atom(i)
          write(DECinfo%output,*) 'Number of virtupied orbitals assigned = ', nvirt_per_atom(i)
          call lsquit('Secondary orbital assigment is inconsistent &
               & with DEC scheme',DECinfo%output)
       end if
    end do


  end subroutine dec_orbital_sanity_check



  !> \brief Get maximum and average orbital extent (number of atoms used to span each orbital)
  !> for input orbitals.
  !> \author Kasper Kristensen
  subroutine get_orbital_extent_info(norb,Orbitals,max_orbital_extent,av_orbital_extent,orb_index)

    implicit none
    !> Number of orbitals
    integer, intent(in) :: norb
    !> List of orbitals
    type(decorbital), intent(in) :: Orbitals(norb)
    !> Maximum orbital extent
    integer, intent(inout) :: max_orbital_extent
    !> Average orbital extent
    real(realk), intent(inout) :: av_orbital_extent
    !> Atom where orbital with the maximum orbital extent is assigned
    integer, intent(inout) :: orb_index
    integer :: i

    ! Init stuff
    max_orbital_extent = 0
    av_orbital_extent = 0e0_realk
    orb_index=0

    ! Get average and maximum number orbital extent
    do i=1,norb

       ! Max
       if(Orbitals(i)%numberofaos > max_orbital_extent) then
          max_orbital_extent = Orbitals(i)%numberofaos
          orb_index = i
       end if

       ! Average
       av_orbital_extent = av_orbital_extent + real(Orbitals(i)%numberofaos)

    end do

    ! Average orbital extent
    av_orbital_extent = av_orbital_extent / real(norb)


  end subroutine get_orbital_extent_info



  !> \brief Determine which fragments to consider
  !> For atom-based DEC - the atoms with orbitals assigned
  !> For DECCO - all occupied orbitals (but not the core orbitals if frozencore is used)
  !> \author Kasper Kristensen
  !> \date October 2013
  subroutine which_fragments_to_consider(ncore,nocc,nvirt,natoms,&
       & OccOrbitals,UnoccOrbitals,dofrag,PhantomAtom)
    !> Number of core orbitals in full molecule
    integer,intent(in) :: ncore
    !> Number of occupied orbitals in full molecule
    integer,intent(in) :: nocc
    !> Number of virtupied orbitals in full molecule
    integer,intent(in) :: nvirt
    !> Number of atoms in full molecule
    integer,intent(in) :: nAtoms
    !> Occupied orbitals in DEC format
    type(decorbital), intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals in DEC format
    type(decorbital), intent(in) :: UnoccOrbitals(nvirt)
    !> dofrag(P) is true/false if atom P has one or more/zero orbitals assigned.
    logical,intent(inout) :: dofrag(natoms)
    !> Which atoms are Phantom Atoms
    logical, intent(in) :: PhantomAtom(natoms)

    if(DECinfo%DECCO) then
       ! Never make a fragment for a core orbital!
       ! Core orbitals are absorbed into valence fragments, see GenerateOrbitals_DECCO
       dofrag=.true.
       dofrag(1:ncore)=.false.
    else
       call which_atoms_have_orbitals_assigned(ncore,nocc,nvirt,natoms,&
            & OccOrbitals,UnoccOrbitals,dofrag,PhantomAtom)
    end if

  end subroutine which_fragments_to_consider


  !> \brief Determine which atoms have one or more orbitals assigned.
  !> \author Kasper Kristensen
  !> \date October 2013
  subroutine which_atoms_have_orbitals_assigned(ncore,nocc,nvirt,natoms,&
       & OccOrbitals,UnoccOrbitals,dofrag,PhantomAtom)
    implicit none
    !> Number of core orbitals in full molecule
    integer,intent(in) :: ncore
    !> Number of occupied orbitals in full molecule
    integer,intent(in) :: nocc
    !> Number of virtupied orbitals in full molecule
    integer,intent(in) :: nvirt
    !> Number of atoms in full molecule
    integer,intent(in) :: nAtoms
    !> Occupied orbitals in DEC format
    type(decorbital), intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals in DEC format
    type(decorbital), intent(in) :: UnoccOrbitals(nvirt)
    !> dofrag(P) is true/false if atom P has one or more/zero orbitals assigned.
    logical,intent(inout) :: dofrag(natoms)
    !> Which atoms are Phantom Atoms
    logical, intent(in) :: PhantomAtom(natoms)
    !local 
    integer, dimension(natoms) :: nocc_per_atom, nvirt_per_atom
    integer :: i

    ! Number of orbitals per atom
    if(DECinfo%frozencore) then ! only count valence orbitals
       nocc_per_atom =  get_number_of_orbitals_per_atom(OccOrbitals,nocc,natoms,.true.,offset=ncore)
    else
       nocc_per_atom =  get_number_of_orbitals_per_atom(OccOrbitals,nocc,natoms,.true.)
    end if

    nvirt_per_atom =  get_number_of_orbitals_per_atom(UnOccOrbitals,nvirt,natoms,.true.)

    ! Which fragments to consider
    dofrag=.true.
    do i=1,natoms
       if(DECinfo%onlyoccpart) then
          ! Only consider occupied orbitals

          if( (nocc_per_atom(i)==0) ) then
             dofrag(i)=.false.
          else
             if(PhantomAtom(i))then
                print*,'ERROR   nocc_per_atom',nocc_per_atom(i),'nvirt_per_atom(i)',nvirt_per_atom(i)
                print*,'ERROR   i',i,'PhantomAtom',PhantomAtom(i)
                dofrag(i)=.false.
                print*,'Setting dofrag to false'
             endif
          endif
       elseif(DECinfo%onlyvirtpart) then
          ! Only consider virtual orbitals
          if( (nvirt_per_atom(i)==0) ) then
             dofrag(i)=.false.
          else
             if(PhantomAtom(i))then
                print*,'ERROR   nocc_per_atom',nocc_per_atom(i),'nvirt_per_atom(i)',nvirt_per_atom(i)
                print*,'ERROR   i',i,'PhantomAtom',PhantomAtom(i)
                dofrag(i)=.false.
                print*,'Setting dofrag to false'
             endif
          endif
       else
          ! Consider occupied as well as virtupied orbitals
          if( (nocc_per_atom(i)==0) .and. (nvirt_per_atom(i)==0) )then
             dofrag(i)=.false.
          else
             if(PhantomAtom(i))then
                print*,'ERROR   nocc_per_atom',nocc_per_atom(i),'nvirt_per_atom(i)',nvirt_per_atom(i)
                print*,'ERROR   i',i,'PhantomAtom',PhantomAtom(i)
                dofrag(i)=.false.
                print*,'Setting dofrag to false'
             endif
          endif
       end if
    end do

    if( DECinfo%only_n_frag_jobs>0)then

       if(DECinfo%only_n_frag_jobs == 1)then
          print *,"HACK TO ONLY DO ONE FRAGMENT JOB"
       else
          print *,"HACK TO ONLY DO",DECinfo%only_n_frag_jobs," FRAGMENT JOBS"
       endif

       dofrag = .false.
       do i = 1, DECinfo%only_n_frag_jobs
          dofrag(DECinfo%frag_job_nr(i)) = .true.
       enddo

    endif

  end subroutine which_atoms_have_orbitals_assigned

end module orbital_operations
