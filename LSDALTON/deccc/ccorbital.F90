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


  ! DEC DEPENDENCIES (within deccc directory)    
  ! *****************************************
  use dec_fragment_utils


contains


  !> \brief Generate ocupied and virtual orbitals.
  !> \author Kasper Kristensen
  !> \date November 2011
  subroutine GenerateOrbitals_driver(MyMolecule,mylsitem,nocc,nunocc,natoms, &
       & OccOrbitals, UnoccOrbitals)

    implicit none
    !> Full molecule structure ( Note: MyMolecule is only changed if modbasis=.true.!)
    type(fullmolecule), intent(in) :: MyMolecule
    !> LSitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> Number of occupied orbitals
    integer,intent(in) :: nocc
    !> Number of unoccupied orbitals
    integer,intent(in) :: nunocc
    !> Number of atoms in the molecule
    integer,intent(in) :: natoms
    !> Occupied orbitals to be generated
    type(decorbital),intent(inout) :: OccOrbitals(nocc)
    !> Unoccupied orbitals to be generated
    type(decorbital),intent(inout) :: UnoccOrbitals(nunocc)


    DECCO: if(DECinfo%DECCO) then
       ! Orbital-based DEC
       call GenerateOrbitals_DECCO(nocc,nunocc,natoms, &
            & MyMolecule,DECinfo%simple_orbital_threshold,OccOrbitals,UnoccOrbitals)
       call print_orbital_info_DECCO(mylsitem,nocc,nunocc,OccOrbitals,&
            & UnoccOrbitals,MyMolecule%ncore)
       return

    else

       ! Atom-based
       OrbitalGeneration: if(DECinfo%read_dec_orbitals) then ! Read DEC orbitals form file
          call read_DECorbitals_from_file(nocc,nunocc,&
               &OccOrbitals,UnoccOrbitals)
       else ! Generate DEC orbitals

          if(DECinfo%BoughtonPulay) then
             write(DECinfo%output,*) 'Generating occupied DEC orbitals using Boughton-Pulay criteria'
             call GenerateOrbitals_BP(OccOrbitals,nOcc,0,MyMolecule,&
                  & DECinfo%mulliken_threshold, DECinfo%simple_mulliken_threshold,&
                  & DECinfo%approximated_norm_threshold,.FALSE.,DECinfo%output)
             if(DECinfo%PL>0) call PrintOrbitalsInfo(OccOrbitals,nOcc,DECinfo%output)

             write(DECinfo%output,*) 'Generating unoccupied DEC orbitals using Boughton-Pulay criteria'
             call GenerateOrbitals_BP(UnoccOrbitals,nUnocc,nOcc,MyMolecule,&
                  & DECinfo%mulliken_threshold, DECinfo%simple_mulliken_threshold,&
                  & DECinfo%approximated_norm_threshold,.FALSE.,DECinfo%output)
             if(DECinfo%PL>0) call PrintOrbitalsInfo(UnoccOrbitals,nUnocc,DECinfo%output)

             ! For Boughton-Pulay, reassign orbitals originally assigned to hydrogen
             if(DECinfo%AbsorbHatoms) then
                call reassign_orbitals(nocc,OccOrbitals,natoms,MyMolecule%DistanceTable,mylsitem)
                call reassign_orbitals(nunocc,UnOccOrbitals,natoms,MyMolecule%DistanceTable,mylsitem)
             end if

          else ! Simple Lowdin charge procedure to determine atomic extent

             write(DECinfo%output,*) 'Generating DEC orbitals using simple Lowdin charge analysis'

             call GenerateOrbitals_simple(nocc,nunocc,natoms, &
                  & MyMolecule,MyLsitem,DECinfo%simple_orbital_threshold,OccOrbitals,UnoccOrbitals)
             if(DECinfo%PL>0) call PrintOrbitalsInfo(OccOrbitals,nocc,DECinfo%output)
             if(DECinfo%PL>0) call PrintOrbitalsInfo(UnoccOrbitals,nUnocc,DECinfo%output)

          end if

       end if OrbitalGeneration

    end if DECCO

    ! If we want to simulate full calculation, we simply assign ALL orbitals to atom 1
    if(DECinfo%simulate_full) then
       write(DECinfo%output,*) 'THIS IS A SIMULATED FULL MOLECULAR CALCULATION!'
       call adjust_orbitals_for_full_simulation(nocc,nunocc,&
            &OccOrbitals,UnoccOrbitals,natoms,DECinfo%simulate_natoms)
    end if

    ! Set secondary atom for orbitals (see set_secondary_atom_for_orbitals)
    call set_secondary_atom_for_orbitals(natoms,nocc,nunocc,MyMolecule,mylsitem,&
         & OccOrbitals,UnoccOrbitals)

    ! Print orbital info
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    if(DECinfo%BoughtonPulay) then
       write(DECinfo%output,*) 'Boughton-Pulay based orbital assignment: Summary'
       write(DECinfo%output,*) '------------------------------------------------'
    else
       write(DECinfo%output,*) 'Simple Lowdin-based orbital assignment: Summary'
       write(DECinfo%output,*) '-----------------------------------------------'
    end if
    call print_orbital_info(mylsitem,nocc,natoms,nunocc,OccOrbitals,UnoccOrbitals,ncore=MyMolecule%ncore)


    ! Check that assigment is meaningful
    call dec_orbital_sanity_check(natoms,nocc,nunocc,OccOrbitals,&
         & UnoccOrbitals,MyMolecule)

    ! Write orbitals to file
    call write_DECorbitals_to_file(nocc,nunocc,&
         &OccOrbitals,UnoccOrbitals)


  end subroutine GenerateOrbitals_driver


  !> \brief Set secondary central atom for orbitals. 
  !> For example, if virtual orbital "a" is assigned to a hydrogen atom "H_A" for which 
  !> there are no occupied orbitals assigned, the secondary          
  !> central atom for orbital "a" will be the atom closest to "H_A" which has a nonzero number
  !> of occupied AND virtual orbitals in the original assignment.
  !> \author Kasper Kristensen
  !> \date June 2014
  subroutine set_secondary_atom_for_orbitals(natoms,nocc,nunocc,MyMolecule,mylsitem,&
       & OccOrbitals,UnoccOrbitals)
    implicit none

    !> Number of atoms in the molecule
    integer,intent(in) :: natoms
    !> Number of occupied orbitals
    integer,intent(in) :: nocc
    !> Number of unoccupied orbitals
    integer,intent(in) :: nunocc
    !> Full molecule structure ( Note: MyMolecule is only changed if modbasis=.true.!)
    type(fullmolecule), intent(in) :: MyMolecule
    !> LSitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> Occupied orbitals to be generated
    type(decorbital),intent(inout) :: OccOrbitals(nocc)
    !> Unoccupied orbitals to be generated
    type(decorbital),intent(inout) :: UnoccOrbitals(nunocc)
    integer :: nocc_per_atom(natoms), nunocc_per_atom(natoms),offset,P,i,central,secondary
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
    nunocc_per_atom =  get_number_of_orbitals_per_atom(UnoccOrbitals,nunocc,natoms,.true.)

    ! First, set secondary atom equal to central atom
    do i=1,nocc
       OccOrbitals(i)%secondaryatom=OccOrbitals(i)%centralatom
    end do
    do i=1,nunocc
       UnoccOrbitals(i)%secondaryatom=UnoccOrbitals(i)%centralatom
    end do

    ! Now we determine if any atoms have zero occ and nonzero virtual orbitals assigned (or vice versa)
    ! to determine which orbitals need a new secondary atom.
    ! Also, for each atom we sort the other atoms according to distance.
    reset_secondary_atom=.false.
    call mem_alloc(DistSortedIdx,natoms,natoms)
    do P=1,natoms
       ! Reset secondary atom?
       if(nocc_per_atom(P)==0 .and. nunocc_per_atom(P)/=0) reset_secondary_atom(P)=.true.
       if(nocc_per_atom(P)/=0 .and. nunocc_per_atom(P)==0) reset_secondary_atom(P)=.true.
       
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
             if(nocc_per_atom(secondary)/=0 .and. nunocc_per_atom(secondary)/=0) then
                found=.true.
                OccOrbitals(i)%secondaryatom = secondary
                exit Ploop
             end if
          end do Ploop
          ! Sanity check
          if(.not. found) then
             call print_orbital_info(mylsitem,nocc,natoms,nunocc,OccOrbitals,UnoccOrbitals,&
                  & ncore=MyMolecule%ncore)
             call lsquit('set_secondary_atom_for_orbitals: Occupied orbital reassigning failed!')
          end if

       end if
    end do

    ! Reassign unoccupied orbitals
    do i=1,nunocc
       central = UnoccOrbitals(i)%centralatom
       if(reset_secondary_atom(central)) then

          ! Find atom closest to current central atom which has both occ and virt orbitals assigned
          found=.false.
          Ploop2: do P=1,natoms
             ! Secondary atom from distance list
             secondary = DistSortedIdx(P,central)
             if(nocc_per_atom(secondary)/=0 .and. nunocc_per_atom(secondary)/=0) then
                found=.true.
                UnoccOrbitals(i)%secondaryatom = secondary
                exit Ploop2
             end if
          end do Ploop2
          ! Sanity check
          if(.not. found) then
             call print_orbital_info(mylsitem,nocc,natoms,nunocc,OccOrbitals,UnoccOrbitals,&
                  & ncore=MyMolecule%ncore)
             call lsquit('set_secondary_atom_for_orbitals: Unoccupied orbital reassigning failed!')
          end if

       end if
    end do

    call mem_dealloc(DistSortedIdx)

  end subroutine set_secondary_atom_for_orbitals



  !> \brief Initialize orbital, assign center and list of significant atoms.
  function orbital_init(orb_number,central_atom,num_atoms,list_atoms) result(myorbital)

    implicit none
    type(decorbital) :: myorbital
    integer, intent(in) :: orb_number
    integer, intent(in) :: central_atom
    integer, intent(in) :: num_atoms
    integer, dimension(num_atoms), intent(in) :: list_atoms
    integer :: i

    myorbital%orbitalnumber = orb_number
    myorbital%centralatom = central_atom
    myorbital%numberofatoms = num_atoms

    ! Set secondary atom equal to central atom in initialization (may be changed later)
    myorbital%secondaryatom = central_atom

    call mem_alloc(myorbital%atoms,num_atoms)
    myorbital%atoms = 0

    do i=1,num_atoms
       myorbital%atoms(i) = list_atoms(i)
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
    OrbitalCopy%numberofatoms = OriginalOrbital%numberofatoms
    OrbitalCopy%secondaryatom = OriginalOrbital%secondaryatom

    call mem_alloc(OrbitalCopy%atoms,OrbitalCopy%numberofatoms)
    OrbitalCopy%atoms = OriginalOrbital%atoms

  end subroutine copy_orbital



  !> \brief Generate orbitals for a given set (occ or virt) using Boughton-Pulay criteria
  subroutine GenerateOrbitals_BP(orbital_set,size_of_set,offset,MyMolecule,&
       & mulliken_threshold2,simple_mulliken_threshold2,&
       & approximated_norm_threshold,modbasis,lu_output)

    implicit none
    integer, intent(in) :: size_of_set,lu_output
    logical,intent(in) :: modbasis !modify basis
    type(decorbital), intent(inout), dimension(size_of_set) :: orbital_set
    ! Note: MyMolecule is only changed if modbasis=.true.!
    type(fullmolecule), intent(in) :: MyMolecule
    integer, intent(in) :: offset  ! 0 - for occupied, nocc - for virtual set
    real(realk), intent(in) :: mulliken_threshold2
    logical, intent(in) :: simple_mulliken_threshold2
    real(realk),intent(in) :: approximated_norm_threshold

    real(realk), pointer :: CtS(:,:), smallS(:,:)
    real(realk), pointer :: gross_charge(:)
    integer, pointer :: atomic_idx(:)
    integer :: num_atoms_to_consider,i,j,nbasis,natoms,atom,central_atom,k,n, &
         selected_nbasis,bas_offset,atom_k
    integer, pointer :: list_of_atoms_to_consider(:)
    real(realk) :: approx_vector_norm, full_vector_norm, approximated_norm, error, half_norm
    real(realk), pointer :: correct_vector_moS(:), approximated_orbital(:), tmp1(:)
    real(realk), pointer :: ShalfC(:,:)
    integer :: tki,atom_start_i,atom_end_i,gamma,catom(1),central_atom2
    integer, pointer :: atom_start(:) => null()
    integer, pointer :: atom_end(:) => null()
    real(realk) :: twall,tcpu
    real(realk),pointer :: basis(:,:)


    call LSTIMER('START',tcpu,twall,DECinfo%output)
    nbasis = MyMolecule%nbasis
    natoms = MyMolecule%natoms
    call mem_alloc(basis,nbasis,nbasis)
    basis(1:nbasis,1:MyMolecule%nocc) = MyMolecule%Co(1:nbasis,1:MyMolecule%nocc)
    basis(1:nbasis,MyMolecule%nocc+1:nbasis) = MyMolecule%Cv(1:nbasis,1:MyMolecule%nunocc)
    call mem_alloc(gross_charge,natoms)

    ! half transformed overlap
    call mem_alloc(CtS,nbasis,nbasis)
    CtS = 0.0E0_realk

    if(.not. DECinfo%Mulliken) then ! Construct matrix used for Lowdin population analysis
       call mem_alloc(ShalfC,nbasis,nbasis)
       call Get_matrix_for_lowdin_analysis(MyMolecule, ShalfC)
    end if

    ! Get matrix C^T S
    call dgemm('t','n',nbasis,nbasis,nbasis,1.0E0_realk,basis,nbasis, &
         & MyMolecule%overlap,nbasis,0.0E0_realk,CtS,nbasis)

    do i=1,size_of_set


       if(DECinfo%Distance) then
          ! Use Distance for orbital "i" to determine central atom 
          call GetDistanceCentralAtom(offset+i,nbasis,natoms,MyMolecule,central_atom2)
       endif
       
       if(DECinfo%Mulliken) then
          ! Use Mulliken population analysis for orbital "i": gross_charge -> Mulliken
          call GetMullikenVector(offset+i,nbasis,natoms,MyMolecule,gross_charge,CtS)
       else
          ! Use Lowdin population analysis for orbital "i": gross_charge -> Lowdin
          call GetLowdinVector(offset+i,nbasis,natoms,MyMolecule,ShalfC,gross_charge)
       end if

       if(simple_mulliken_threshold2) then

          ! -- Count number of atoms with Mulliken charge above the threshold
          num_atoms_to_consider = 0
          do j=1,natoms
             if(abs(gross_charge(j)) > mulliken_threshold2) &
                  num_atoms_to_consider = num_atoms_to_consider + 1
          end do

          if(DECinfo%Distance) then
             if( abs(gross_charge(central_atom2)) > mulliken_threshold2)THEN
                !the new central atom according to distance criteria
                !is already included in the list
             else
                !add the central atom 
                num_atoms_to_consider = num_atoms_to_consider + 1                
             endif
          endif

          call mem_alloc(list_of_atoms_to_consider,num_atoms_to_consider)
          list_of_atoms_to_consider = 0
          atom=1
          do j=1,natoms
             if(abs(gross_charge(j)) > mulliken_threshold2) then
                list_of_atoms_to_consider(atom) = j
                atom = atom + 1
             end if
          end do
          ! --
          if(DECinfo%Distance) then
             if( abs(gross_charge(central_atom2)) > mulliken_threshold2)THEN
                !the new central atom according to distance criteria
                !is already included in the list
             else
                !add the central atom
                list_of_atoms_to_consider(atom) = central_atom2 
             endif
             central_atom = central_atom2 
          else
             ! -- number of atom with max value of the Mulliken charge             
             catom = maxloc(abs(gross_charge),dim=1,MASK=.NOT.MyMolecule%PhantomAtom)
             central_atom = catom(1)
          endif
       else


          ! sort mulliken charge/atomic density keeping track on the original index
          call mem_alloc(atomic_idx,natoms)
          do k=1,natoms
             atomic_idx(k) = k
          end do

          gross_charge = abs(gross_charge)
          !obtain the central atom - as the atom (which is not Phantom) with biggest 
          !mulliken charge
          catom = maxloc(gross_charge,dim=1,MASK=.NOT.MyMolecule%PhantomAtom)
          call real_inv_sort_with_tracking(gross_charge,atomic_idx,natoms)

          if(DECinfo%PL>0) then
             do n=1,natoms
                write(lu_output,'(2i4,f16.10)') n,atomic_idx(n),gross_charge(n)
             end do
          end if

          ! get an optimal set of atoms by taking them one by one until the
          ! convergence i.e. until orbital is correctly approximated
          do n=1,natoms

             ! take 'n' first atoms
             num_atoms_to_consider = n
             call mem_alloc(list_of_atoms_to_consider,n)
             do k=1,n
                list_of_atoms_to_consider(k) = atomic_idx(k)
             end do

             ! sort list of atoms
             call int_sort(list_of_atoms_to_consider,n,n)

             ! number of basis function on selected atoms
             selected_nbasis = 0
             do k=1,n
                selected_nbasis = selected_nbasis + &
                     MyMolecule%atom_size(list_of_atoms_to_consider(k))
             end do

             ! ao overlap on selected atoms
             call mem_alloc(smallS,selected_nbasis,selected_nbasis)
             call adjust_square_matrix(MyMolecule%overlap,smallS,&
                  list_of_atoms_to_consider,MyMolecule%atom_size, &
                  MyMolecule%atom_start,MyMolecule%atom_end, &
                  nbasis,natoms,selected_nbasis,n)

             ! mo overlap vector for given orbital
             call mem_alloc(correct_vector_moS,selected_nbasis)
             correct_vector_moS = 0.0E0_realk
             bas_offset = 1
             do k=1,n
                atom_k = list_of_atoms_to_consider(k)
                correct_vector_moS(bas_offset:bas_offset + MyMolecule%atom_size(atom_k) - 1) = &
                     CtS(offset+i,MyMolecule%atom_start(atom_k):MyMolecule%atom_end(atom_k))
                bas_offset = bas_offset + MyMolecule%atom_size(atom_k)
             end do

             ! approximated orbital
             call mem_alloc(approximated_orbital,selected_nbasis)
             approximated_orbital = 0.0E0_realk
             call solve_linear_equations(smallS,approximated_orbital, &
                  correct_vector_moS,selected_nbasis)

             ! transform ao overlap using approximated orbital to check its norm
             call mem_alloc(tmp1,selected_nbasis)
             call dgemm('n','n',selected_nbasis,1,selected_nbasis,1.0E0_realk,smallS,selected_nbasis,&
                  & approximated_orbital,selected_nbasis,0.0E0_realk,tmp1,nbasis)
             approximated_norm = dot_product(approximated_orbital,tmp1)
             call mem_dealloc(tmp1)

             ! Boughton-Pulay: Error = 1 - norm of approximated orbital
             error = abs(1.0E0_realk - approximated_norm)


             ! debug info
             DebugInfo : if(DECinfo%PL>0) then
                call mem_alloc(tmp1,selected_nbasis)
                tmp1 = matmul(smallS,approximated_orbital)
                write(lu_output,'(a,f16.6)') 'New orbital norm :', approximated_norm
                write(lu_output,'(a,f16.6)') 'Error in norm    :', error
                write(lu_output,'(a,f16.6)') 'Old orbital norm :', &
                     dot_product(basis(:,offset+i), &
                     matmul(MyMolecule%overlap,basis(:,offset+i)))
                write(lu_output,'(a)') '-------'
                write(lu_output,'(a)') 'mu  FullMoS(mu) ApproxMoS(mu)'
                do k=1,selected_nbasis
                   write(lu_output,'(i4,3f16.6)') &
                        k,correct_vector_moS(k),tmp1(k),correct_vector_moS(k)-tmp1(k)
                end do
                write(lu_output,'(a)') '-------'
                call mem_dealloc(tmp1)
             end if DebugInfo

             ! delete stuff
             call mem_dealloc(correct_vector_moS)
             call mem_dealloc(smallS)

             ! break procedure keeping list of atoms
             if( error < approximated_norm_threshold .or. n == natoms) then
                if(modbasis)then
                   tki = 1
                   do k=1,n
                      atom_k = list_of_atoms_to_consider(k)
                      atom_start_i = MyMolecule%atom_start(atom_k)
                      atom_end_i = MyMolecule%atom_end(atom_k)
                      do gamma=atom_start_i,atom_end_i
                         basis(gamma,i) = approximated_orbital(tki)
                         tki = tki+1
                      enddo
                   enddo
                endif
                call mem_dealloc(approximated_orbital)
                exit
             else
                call mem_dealloc(approximated_orbital)
                call mem_dealloc(list_of_atoms_to_consider)
             endif
          end do

          ! -- number of atom with max value of the Mulliken charge
          central_atom = atomic_idx(1)
          iF(MyMolecule%PhantomAtom(atomic_idx(1)))THEN
             central_atom = catom(1)
          ENDIF
          call mem_dealloc(atomic_idx)
       endif

       ! -- print extended orbital info
       if(DECinfo%PL>0) then
          write(DECinfo%output,'(/,a,i4)') 'Orbital                 : ',i
          write(DECinfo%output,'(a,i4)')   'Orbital # in full basis : ',offset+i
          write(DECinfo%output,'(a,i4,/)') 'Central atom            : ',central_atom
          write(DECinfo%output,'(a,i4,/)') 'Central atom2           : ',central_atom2
          write(DECinfo%output,'(a,i4,/)') 'Number of atoms         : ',num_atoms_to_consider
          do j=1,num_atoms_to_consider
             write(DECinfo%output,*) list_of_atoms_to_consider(j)
          end do

          write(DECinfo%output,'(2/,a)') 'i     Mul.Ch.'
          do j=1,natoms
             write(DECinfo%output,'(i4,3f16.10)') &
                  j,gross_charge(j)
          end do

       end if
       ! --

       ! -- Create orbital
       orbital_set(i) = orbital_init(i,central_atom, &
            num_atoms_to_consider,list_of_atoms_to_consider)

       call mem_dealloc(list_of_atoms_to_consider)

    end do

    call mem_dealloc(basis)
    call mem_dealloc(gross_charge)
    call mem_dealloc(CtS)
    if(.not. DECinfo%Mulliken) then
       call mem_dealloc(ShalfC)
    end if

    call LSTIMER('GenerateOrb',tcpu,twall,DECinfo%output)

  end subroutine GenerateOrbitals_BP



  !> \brief Generate DEC orbitals for both occ and virt orbitals. For each orbital:
  !> 1. List atoms according to Lowdin charge for that orbital.
  !> 2. Include atoms from this list until "1 minus the sum of these Lowdin charges"
  !>    is smaller than the input approximated_norm_threshold.
  !> 3. Regarding orbital assigning we ensure that we have one atomic fragment for each heavy atom
  !>    + one atomic fragment for all hydrogens (if any) assigned to that heavy atom.
  !> \author Kasper Kristensen
  !> \date September 2011
  subroutine GenerateOrbitals_simple(nocc,nunocc,natoms, &
       & MyMolecule,MyLsitem,approximated_norm_threshold,OccOrbitals,UnoccOrbitals)

    implicit none
    !> Number of occupied orbitals
    integer, intent(in) :: nocc
    !> Number of unoccupied orbitals
    integer, intent(in) :: nunocc
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
    type(decorbital), intent(inout), dimension(nunocc) :: UnoccOrbitals
    integer :: i,j,atom,central_atom,n,norbital_extent,nbasis,heavyatom,ni
    integer, pointer :: list_of_atoms_to_consider(:)
    real(realk) :: error,charge,mindist,twall,tcpu,maxlowdin
    logical :: keepon,ReAssignVirtHydrogenOrbs
    real(realk), pointer :: ShalfC(:,:)
    real(realk), pointer :: lowdin_charge(:,:)
    integer, pointer :: atomic_idx(:,:), countOcc(:), countUnocc(:),central_atom2(:)
    integer :: maxidx,offset,changedatom,nreass
    logical,pointer ::which_hydrogens(:), dofrag(:)
    real(realk),pointer :: tmplowdin_charge(:)
    integer,pointer :: tmpatomic_idx(:)
    integer :: nunoccperatom


    call LSTIMER('START',tcpu,twall,DECinfo%output)

    ! Init stuff
    nbasis = MyMolecule%nbasis

    ! For orbital assignment bookkeeping, only consider valence orbitals
    offset = MyMolecule%ncore  

    call mem_alloc(lowdin_charge,natoms,nbasis)
    call mem_alloc(ShalfC,nbasis,nbasis)
    call mem_alloc(atomic_idx,natoms,nbasis)
    if(DECinfo%Distance) then
       call mem_alloc(central_atom2,nbasis)
    endif

    ! Get Lowdin matrix S^{1/2} C
    call Get_matrix_for_lowdin_analysis(MyMolecule, ShalfC)


    ! ***********************************
    ! Get Lowdin charges for all orbitals
    ! ***********************************
    GetLowdinCharges: do i=1,nbasis

       if(DECinfo%Distance) then
          ! Use Distance criteria to determine central atom for orbital "i"  
          call GetDistanceCentralAtom(i,nbasis,natoms,MyMolecule,central_atom2(i))
       endif

       ! Get vector with Lowdin charges for all atoms for orbital "i"
       call GetLowdinVector(i,nbasis,natoms,MyMolecule,ShalfC,lowdin_charge(:,i) )

       ! Sort Lowdin charges
       call real_inv_sort_with_tracking(lowdin_charge(:,i),atomic_idx(:,i),natoms)

       IF(MyMolecule%PhantomAtom(atomic_idx(1,i)))THEN
          !the first atom in the lowdin_charge(:,i) list is a Phantom atom
          !we reorder the lowdin_charge(:,i) list to put a non Phantom atom on top
          call mem_alloc(tmplowdin_charge,nAtoms)
          call mem_alloc(tmpatomic_idx,nAtoms)
          ni=0
          do n=1,natoms
             tmplowdin_charge(n) = lowdin_charge(n,i)
             tmpatomic_idx(n) = atomic_idx(n,i)
             IF(.NOT.MyMolecule%PhantomAtom(atomic_idx(n,i)))THEN
                IF(ni.EQ.0)THEN
                   !found the first non Phantom atom in lowdin list
                   ni = n
                ENDIF
             ENDIF
          enddo
          !place the first non Phantom atom first
          lowdin_charge(1,i) = tmplowdin_charge(ni) 
          atomic_idx(1,i) = tmpatomic_idx(ni)
          !place the other atoms in list
          do n=1,ni-1                               
             lowdin_charge(1+n,i) = tmplowdin_charge(n)
             atomic_idx(1+n,i) = tmpatomic_idx(n)                
          enddo
          !number ni should not be placed because it is now number 1
          !place the rest of the atoms in the list
          do n=ni+1,nAtoms
             lowdin_charge(n,i) = tmplowdin_charge(n)
             atomic_idx(n,i) = tmpatomic_idx(n)                
          enddo
          call mem_dealloc(tmplowdin_charge)
          call mem_dealloc(tmpatomic_idx)
       ENDIF
    end do GetLowdinCharges
    call mem_dealloc(ShalfC)



    ! *********************************************
    ! Initial orbital assignment and orbital extent
    ! *********************************************
    OrbitalLoop: do i=1,nbasis

       charge = 0E0_realk

       ! Central atom is the one with largest Lowdin charge (although this may be modified below)
       if(DECinfo%Distance) then
          !central atom determined from shortest distance from MO
          central_atom = central_atom2(i)
          IF(atomic_idx(1,i).EQ.central_atom)THEN
             !do nothing the closest atom is also the one with the largest Lowdin charge
          ELSE
             !we reorder the lowdin_charge(:,i) list to put the central atom on top
             call mem_alloc(tmplowdin_charge,nAtoms)
             call mem_alloc(tmpatomic_idx,nAtoms)
             do n=1,natoms
                tmplowdin_charge(n) = lowdin_charge(n,i)
                tmpatomic_idx(n) = atomic_idx(n,i)
                IF(atomic_idx(n,i).EQ.central_atom)THEN
                   !found the central_atom in lowdin list
                   ni = n
                ENDIF
             enddo
             !place the central atom first
             lowdin_charge(1,i) = tmplowdin_charge(ni) 
             atomic_idx(1,i) = tmpatomic_idx(ni)
             !place the non central atoms in list
             do n=1,ni-1                               
                lowdin_charge(1+n,i) = tmplowdin_charge(n)
                atomic_idx(1+n,i) = tmpatomic_idx(n)                
             enddo
             !number ni should not be placed because it is now number 1
             !place the rest of the atoms in the list
             do n=ni+1,nAtoms
                lowdin_charge(n,i) = tmplowdin_charge(n)
                atomic_idx(n,i) = tmpatomic_idx(n)                
             enddo
             call mem_dealloc(tmplowdin_charge)
             call mem_dealloc(tmpatomic_idx)
          ENDIF
       ELSE
          central_atom = atomic_idx(1,i)
       ENDIF
       ! Add atoms until sum of Lowdin charges is close enough to 1
       LowdinAddLoop: do n=1,natoms
          charge = charge + lowdin_charge(n,i)
          error = 1E0_realk - charge

          if(error < approximated_norm_threshold .or. n==natoms ) then 
             ! atom list converged for orbital i

             ! Set list of atoms to consider for orbital and exit loop
             norbital_extent = n
             call mem_alloc(list_of_atoms_to_consider,norbital_extent)
             do atom=1,norbital_extent
                list_of_atoms_to_consider(atom) = atomic_idx(atom,i)
             end do
             exit LowdinAddLoop

          end if

       end do LowdinAddLoop


       ! Print orbital info for high print levels
       if(DECinfo%PL>1) then
          write(DECinfo%output,'(1X,a,i10)') 'ORBITAL: ', i
          write(DECinfo%output,*) '-------------------------------------'
          write(DECinfo%output,'(1X,a,100i5)')    'ATOMS  : ', list_of_atoms_to_consider
          write(DECinfo%output,'(1X,a,100f10.3)') 'LOWDIN : ', lowdin_charge(1:norbital_extent,i)
          write(DECinfo%output,'(1X,a,i6)') 'Central Atom : ', central_atom
          if(DECinfo%Distance) then
             write(DECinfo%output,'(1X,a,i6)') 'Central Atom : ', central_atom2(i)
          endif
          write(DECinfo%output,*)
       end if

       ! -- Create orbital
       if(i<=nocc) then   ! Occupied orbital
          OccOrbitals(i) = orbital_init(i,central_atom, &
               norbital_extent,list_of_atoms_to_consider)
       else  ! unoccupied orbital
          UnoccOrbitals(i-nocc) = orbital_init(i,central_atom, &
               norbital_extent,list_of_atoms_to_consider)
       end if

       call mem_dealloc(list_of_atoms_to_consider)

    end do OrbitalLoop


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
       call reassign_orbitals(nunocc,UnOccOrbitals,natoms,MyMolecule%DistanceTable,mylsitem)
    end if AbsorbHydrogenAtoms


    HydrogenDebug: if( decinfo%PureHydrogendebug ) THEN

       do i=1,nocc
          OccOrbitals(i)%centralatom = i
       enddo
       nunoccperatom = nunocc/nocc
       do i=1,nunocc
          ReAssignVirtHydrogenOrbs = .TRUE.
          do j=1,nocc
             IF(UnOccOrbitals(i)%centralatom.EQ.j)THEN
                ReAssignVirtHydrogenOrbs = .FALSE.                
             ENDIF
          enddo
          IF(ReAssignVirtHydrogenOrbs)THEN
             do j=1,nocc   
                UnOccOrbitals(i)%centralatom = (nunocc-1)/nunoccperatom + 1
             enddo
          ENDIF
       enddo

    end if HydrogenDebug


    ! ****************************************************************************************
    ! * Reassign: Ensure that all atoms have both occupied and unoccupied orbitals assigned  *
    ! ****************************************************************************************

    REASSIGNING: if( .not. decinfo%PureHydrogendebug ) then

       ! Count # orbitals assigned to each atom
       call mem_alloc(countOcc,natoms)
       call mem_alloc(countUnocc,natoms)
       countOcc=0
       countUnocc=0
       do i=offset+1,nocc
          countocc(OccOrbitals(i)%centralatom) = countocc(OccOrbitals(i)%centralatom)+1
       end do
       do i=1,nunocc
          countunocc(UnoccOrbitals(i)%centralatom) = countunocc(UnoccOrbitals(i)%centralatom)+1
       end do

       ! The atoms which will be central in atomic fragments are those which
       ! at this point have some orbitals assigned
       call mem_alloc(dofrag,natoms)
       dofrag=.false.
       do i=1,natoms
          if( countocc(i)/=0 .or. countunocc(i)/=0 ) dofrag(i)=.true.
       end do

       ! Now reassign to ensure that all atomic fragment have both occupied and unoccupied orbitals
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


             ! Reassign unoccupied orbitals (same procedure as for occ space)
             ! Never reassign virtual orbitals for only virtual partitioning
             UnoccReassign: If(dofrag(atom) .and. countunocc(atom)==0 &
                  & .and. (.not. DECinfo%onlyvirtpart) ) then

                maxlowdin = 0.0_realk
                maxidx = 0
                changedatom=0
                do j=nocc+1,nbasis
                   if(atomic_idx(2,j)==atom .and. lowdin_charge(2,j) > maxlowdin ) then
                      maxlowdin = lowdin_charge(2,j)
                      maxidx = j-nocc   ! unoccupied index in list of unoccupied orbitals
                      changedatom=UnoccOrbitals(j-nocc)%centralatom
                   end if
                end do

                ! Reassign orbital
                if(maxidx/=0) then
                   if(DECinfo%PL>1) write(DECinfo%output,'(1X,a,i8,a,i8)') &
                        & 'Reassign unocc orbital ', maxidx, ' to atom ', atom
                   UnoccOrbitals(maxidx)%centralatom = atom
                   countunocc(changedatom) = countunocc(changedatom) - 1
                   countunocc(atom) = countunocc(atom) + 1
                end if

             end if UnoccReassign


             ! Check that all atoms have either have
             ! (i) both occupied AND unoccupied orbitals assigned   OR
             ! (ii) zero orbitals assigned
             keepon=.false.
             CheckAssignment: do i=1,natoms

                WhichScheme: if(DECinfo%onlyoccpart) then

                   if( (countocc(i)/=0 .and. countunocc(i)==0) ) then
                      ! Still not acceptable orbital distribution - keep on
                      keepon=.true.
                      exit CheckAssignment
                   end if

                elseif(DECinfo%onlyvirtpart) then

                   if( (countocc(i)==0 .and. countunocc(i)/=0) ) then
                      ! Still not acceptable orbital distribution - keep on
                      keepon=.true.
                      exit CheckAssignment
                   end if

                else

                   if( (countocc(i)/=0 .and. countunocc(i)==0) .or. &
                        & (countocc(i)==0 .and. countunocc(i)/=0) ) then
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
       call mem_dealloc(countOcc)
       call mem_dealloc(countUnocc)

    end if REASSIGNING

    call mem_dealloc(which_hydrogens)
    call mem_dealloc(lowdin_charge)
    call mem_dealloc(atomic_idx)
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
  subroutine GenerateOrbitals_DECCO(nocc,nunocc,natoms, &
       & MyMolecule,approximated_norm_threshold,OccOrbitals,UnoccOrbitals)

    implicit none
    !> Number of occupied orbitals
    integer, intent(in) :: nocc
    !> Number of unoccupied orbitals
    integer, intent(in) :: nunocc
    !> Number of atoms in the molecule
    integer,intent(in) :: natoms
    !> Molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Threshold for orbital norm (see above)
    real(realk),intent(in) :: approximated_norm_threshold
    !> Occupied orbitals to create
    type(decorbital), intent(inout), dimension(nocc) :: OccOrbitals
    !> Unoccupied orbitals to create
    type(decorbital), intent(inout), dimension(nunocc) :: UnoccOrbitals
    integer :: i,central_orbital,n,norbital_extent,nbasis,atom,j
    integer, pointer :: list_of_atoms_to_consider(:)
    real(realk) :: error,charge,twall,tcpu
    real(realk), pointer :: ShalfC(:,:)
    real(realk), pointer :: lowdin_charge(:,:)
    integer, pointer :: atomic_idx(:,:)
    integer :: offset
    real(realk),pointer :: DistoccUnocc(:,:),DistOccOcc(:,:),sorted_dists(:)
    integer :: sorted_orbitals(nocc)


    call LSTIMER('START',tcpu,twall,DECinfo%output)

    ! Init stuff
    nbasis = MyMolecule%nbasis

    offset = MyMolecule%ncore  

    call mem_alloc(lowdin_charge,natoms,nbasis)
    call mem_alloc(ShalfC,nbasis,nbasis)
    call mem_alloc(atomic_idx,natoms,nbasis)

    ! Get Lowdin matrix S^{1/2} C
    call Get_matrix_for_lowdin_analysis(MyMolecule, ShalfC)


    ! ***********************************
    ! Get Lowdin charges for all orbitals
    ! ***********************************
    GetLowdinCharges: do i=1,nbasis

       ! Get vector with Lowdin charges for all atoms for orbital "i"
       call GetLowdinVector(i,nbasis,natoms,MyMolecule,ShalfC,lowdin_charge(:,i) )

       ! Sort Lowdin charges
       call real_inv_sort_with_tracking(lowdin_charge(:,i),atomic_idx(:,i),natoms)

    end do GetLowdinCharges
    call mem_dealloc(ShalfC)


    ! Distances between occ and unocc orbitals
    call mem_alloc(DistOccUnocc,nocc,nunocc)
    call general_distance_table(nocc,nunocc,MyMolecule%carmomocc,MyMolecule%carmomvirt,DistOccUnocc)
    ! .. and between occ and occ orbitals
    call mem_alloc(DistOccOcc,nocc,nocc)
    call general_distance_table(nocc,nocc,MyMolecule%carmomocc,MyMolecule%carmomocc,DistOccOcc)
    call mem_alloc(sorted_dists,nocc)

    ! *************************************
    ! Orbital assignment and orbital extent
    ! *************************************
    OrbitalLoop: do i=1,nbasis

       charge = 0E0_realk

       ! Add atoms until sum of Lowdin charges is close enough to 1
       ! **********************************************************
       LowdinAddLoop: do n=1,natoms
          charge = charge + lowdin_charge(n,i)
          error = 1E0_realk - charge

          if(error < approximated_norm_threshold .or. n==natoms ) then 
             ! atom list converged for orbital i

             ! Set list of atoms to consider for orbital and exit loop
             norbital_extent = n
             call mem_alloc(list_of_atoms_to_consider,norbital_extent)
             do atom=1,norbital_extent
                list_of_atoms_to_consider(atom) = atomic_idx(atom,i)
             end do
             exit LowdinAddLoop

          end if

       end do LowdinAddLoop


       ! Print orbital info for high print levels
       ! ****************************************
       if(DECinfo%PL>1) then
          write(DECinfo%output,'(1X,a,i10)') 'ORBITAL: ', i
          write(DECinfo%output,*) '-------------------------------------'
          write(DECinfo%output,'(1X,a,100i5)')    'ATOMS  : ', list_of_atoms_to_consider
          write(DECinfo%output,'(1X,a,100f10.3)') 'LOWDIN : ', lowdin_charge(1:norbital_extent,i)
          write(DECinfo%output,*)
       end if


       ! -- Create orbital
       ! *****************

       if(  (i>offset)  .and.  (i<=nocc) ) then   ! Valence orbital

          ! Central orbital is the occupied orbital itself
          ! ----------------------------------------------
          central_orbital = i
          OccOrbitals(i) = orbital_init(i,central_orbital, &
               norbital_extent,list_of_atoms_to_consider)

       else  ! unoccupied orbital or core orbital

          ! Sort occ orbitals according to distance to orbital "i"
          if(i.le.offset) then  ! core orbital
             sorted_dists = DistOccOcc(:,i)
          else ! unocc orbital
             sorted_dists = DistOccUnocc(:,i-nocc)
          end if
          call real_inv_sort_with_tracking(sorted_dists,sorted_orbitals,nocc)

          ! Assign to nearest valence orbital (ensure that we do not assign to core orbitals)
          Assigning: do n=1,nocc
             if(sorted_orbitals(n)>offset) then
                central_orbital = sorted_orbitals(n)
                exit Assigning
             end if
          end do Assigning

          if(DECinfo%PL>1) then
             write(DECinfo%output,'(1X,a,i10)') 'Sorted occ orbs for orbital: ', i
             write(DECinfo%output,*) '--------------------------------------------------------'
             do j=1,nocc
                write(DECinfo%output,*) sorted_orbitals(i)
             end do
             write(DECinfo%output,*) 'Central orbital: ', central_orbital
          end if

          if(i.le.offset) then  ! core orbital
             OccOrbitals(i) = orbital_init(i,central_orbital, &
                  norbital_extent,list_of_atoms_to_consider)
          else ! Unocc orbital
             UnoccOrbitals(i-nocc) = orbital_init(i,central_orbital, &
                  norbital_extent,list_of_atoms_to_consider)
          end if

       end if

       call mem_dealloc(list_of_atoms_to_consider)

    end do OrbitalLoop



    call mem_dealloc(sorted_dists)
    call mem_dealloc(DistOccUnocc)
    call mem_dealloc(DistOccOcc)
    call mem_dealloc(lowdin_charge)
    call mem_dealloc(atomic_idx)


    ! Sanity check
    call DECCO_assignment_sanity_check(MyMolecule,OccOrbitals,UnoccOrbitals)

    call LSTIMER('GenerateOrb',tcpu,twall,DECinfo%output)

  end subroutine GenerateOrbitals_DECCO




  !> Sanity check for orbital assignment in DECCO
  subroutine DECCO_assignment_sanity_check(MyMolecule,OccOrbitals,UnoccOrbitals)

    implicit none
    !> Molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Occupied orbitals
    type(decorbital), intent(in), dimension(MyMolecule%nocc) :: OccOrbitals
    !> Unoccupied orbitals
    type(decorbital), intent(in), dimension(MyMolecule%nunocc) :: UnoccOrbitals
    integer :: i, nunocc_per_occ(MyMolecule%nocc)

    ! Only occupied partitioning scheme - nothing to check
    if(DECinfo%OnlyOccPart) then
       return
    end if

    call DECCO_get_nvirt_per_occ_orbital(MyMolecule%nocc,MyMolecule%nunocc,&
         & UnoccOrbitals,nunocc_per_occ)

    ! Else we check that each valence orbital has at least on unocc orbital assigned
    ! and that core orbitals have nothing assigned

    ! Valence check
    do i=MyMolecule%ncore+1,MyMolecule%nocc
       if(nunocc_per_occ(i)==0) then
          call lsquit('DECCO_assignment_sanity_check: &
               &Valence orbital has no unoccupied orbital(s) assigned!',-1)
          print *, 'Error for orbital: ',i
          print *, 'You cannot use DECCO for this system - unless you are willing to&
               & change the source code!'
          print *, 'nunocc_per_occ: ', nunocc_per_occ
       end if
    end do

    ! Core check
    do i=1,MyMolecule%ncore
       if(nunocc_per_occ(i)/=0) then
          call lsquit('DECCO_assignment_sanity_check: &
               &Core orbital has unoccupied orbital(s) assigned!',-1)
          print *, 'Error for orbital: ',i
          print *, 'You cannot use DECCO for this system - unless you are willing to&
               & change the source code!'
          print *, 'nunocc_per_occ: ', nunocc_per_occ
       end if
    end do


  end subroutine DECCO_assignment_sanity_check



  !> \brief Get norm of approximate orbital.
  !> \author Kasper Kristensen
  !> \date September 2011
  subroutine GetApproximateOrbitalNorm(MyMolecule,which_atoms,orb_idx,OrbitalNorm)

    implicit none
    !> Information for molecule
    type(fullmolecule), intent(in) :: MyMolecule
    !> which_atoms(A) is true if atom A is used to approximate orbital, false otherwise
    logical, dimension(MyMolecule%natoms), intent(in) :: which_atoms
    !> Orbital index
    integer, intent(in) :: Orb_idx
    !> Norm of approximat orbital
    real(realk), intent(inout) :: OrbitalNorm
    real(realk), pointer :: C(:,:)
    real(realk), pointer :: S(:,:) => null()
    integer, pointer :: atom_start(:) => null()
    integer, pointer :: atom_end(:) => null()
    integer :: atom1,atom2,mu,nu,natoms,nbasis,nocc,nvirt

    ! Init stuff
    nbasis=MyMolecule%nbasis
    nocc=MyMolecule%nocc
    nvirt=MyMolecule%nunocc
    call mem_alloc(C,nbasis,nbasis)
    C(1:nbasis,1:nocc) = MyMolecule%Co(1:nbasis,1:nocc)
    C(1:nbasis,nocc+1:nbasis) = MyMolecule%Cv(1:nbasis,1:nvirt)
    S => MyMolecule%overlap ! AO overlap
    atom_start => MyMolecule%atom_start
    atom_end => MyMolecule%atom_end
    OrbitalNorm = 0E0_realk
    natoms = MyMolecule%natoms


    ! Approximate orbital norm for molecular orbital "i":
    !
    ! < phi | phi > = \sum_{mu in SET} sum_{nu in SET}  C_{mu i} S_{mu nu} C_{nu i}
    !
    ! where SET is the set of atoms defined by which_atoms.
    do atom1=1,natoms
       do atom2=1,natoms

          ! Only add contribution if both mu AND nu are in set
          if( which_atoms(atom1) .and. which_atoms(atom2) ) then
             do mu=atom_start(atom1),atom_end(atom1) ! mu \in atom
                do nu=atom_start(atom2),atom_end(atom2) ! nu \in atom
                   OrbitalNorm = OrbitalNorm + C(mu,orb_idx) * S(mu,nu) * C(nu,orb_idx)
                end do
             end do
          end if

       end do
    end do


    call mem_dealloc(C)
    S => null()
    atom_start => null()
    atom_end => null()

  end subroutine GetApproximateOrbitalNorm




  !> \brief Get matrix [S^{1/2} C] used for Lowdin population analysis
  !> \author Kasper Kristensen
  subroutine Get_matrix_for_lowdin_analysis(MyMolecule, ShalfC)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> S^{1/2} C matrix
    real(realk), dimension(MyMolecule%nbasis,MyMolecule%nbasis) :: ShalfC
    real(realk), pointer :: Shalf(:,:)
    integer :: nbasis,i,j,k
    real(realk),pointer :: basis(:,:)

    nbasis = MyMolecule%nbasis
    call mem_alloc(basis,nbasis,nbasis)
    !basis(1:nbasis,1:MyMolecule%nocc) = MyMolecule%Co(1:nbasis,1:MyMolecule%nocc)
    !basis(1:nbasis,MyMolecule%nocc+1:nbasis) = MyMolecule%Cv(1:nbasis,1:MyMolecule%nunocc)
    do j=1,MyMolecule%nocc
    do i =1,nbasis
    basis(i,j) =MyMolecule%Co(i,j)
    enddo
    enddo
    k = MyMolecule%nocc+1
    do j = 1,MyMolecule%nunocc
    do i =1,nbasis
    basis(i,k) = MyMolecule%Cv(i,j)
    enddo
    k=k+1
    enddo

    ! Get S^{1/2} matrix
    ! ******************
    call mem_alloc(Shalf,nbasis,nbasis)
    call get_power_of_symmetric_matrix(nbasis,0.5E0_realk,MyMolecule%overlap(1:nbasis,1:nbasis),Shalf)

    ! S^{1/2} C
    ! *********
    call dgemm('n','n',nbasis,nbasis,nbasis,1.0E0_realk,Shalf,nbasis, &
         & basis(1:nbasis,1:nbasis),nbasis,0.0E0_realk,ShalfC,nbasis)

    call mem_dealloc(Shalf)
    call mem_dealloc(basis)

  end subroutine Get_matrix_for_lowdin_analysis



  !> \brief Get Mulliken gross charges for a given orbital on all atoms
  subroutine GetMullikenVector(orbI,nbasis,natoms,MyMolecule,charges,CtS)

    implicit none
    integer, intent(in) :: orbI ! orbital number
    integer, intent(in) :: natoms ! number of atoms
    integer, intent(in) :: nbasis ! Number of basis functions
    type(fullmolecule), intent(in) :: MyMolecule
    real(realk), dimension(natoms), intent(inout) :: charges
    real(realk), dimension(nbasis,nbasis), intent(in) :: CtS
    real(realk), pointer :: C(:,:)
    integer, pointer :: atom_start(:) => null()
    integer, pointer :: atom_end(:) => null()
    integer :: atom,nu
    integer :: nocc,nvirt

    ! Init stuff
    nocc=MyMolecule%nocc
    nvirt=MyMolecule%nunocc
    call mem_alloc(C,nbasis,nbasis)
    C(1:nbasis,1:nocc) = MyMolecule%Co(1:nbasis,1:nocc)
    C(1:nbasis,nocc+1:nbasis) = MyMolecule%Cv(1:nbasis,1:nvirt)
    charges=0.0E0_realk
    atom_start => MyMolecule%atom_start
    atom_end => MyMolecule%atom_end

    ! Loop over atoms
    do atom=1,natoms

       ! Mulliken charge for orbital I, on the atom "atom":
       !
       ! chargeI(atom) = sum_{nu \in atom} [C^T S]_{I,nu} C_{nu,I}
       !
       ! where C are MO coefficients and S is AO overlap matrix,
       ! nu is an AO index and I is a MO index for the orbital in question.
       !
       ! It is seen that:
       ! sum_{all atoms} chargeI(atom) = sum_{all nu} [C^T S]_{I,nu} C_{nu,I}
       !                               = [C^T S C]_{I,I}
       !                               = 1
       ! as it must be.

       do nu=atom_start(atom),atom_end(atom) ! nu \in atom
          charges(atom)=charges(atom)+CtS(orbI,nu)*C(nu,orbI)
       end do

    end do


    call mem_dealloc(C)
    atom_start => null()
    atom_end => null()

  end subroutine GetMullikenVector

  !> \brief Get Mulliken gross charges for a given orbital on all atoms
  subroutine GetDistanceCentralAtom(orbI,nbasis,natoms,MyMolecule,central_atom)
    implicit none
    integer, intent(in) :: orbI ! orbital number
    integer, intent(in) :: natoms ! number of atoms
    integer, intent(in) :: nbasis ! Number of basis functions
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
  subroutine GetLowdinVector(orb_idx,nbasis,natoms,MyMolecule,ShalfC,charges)

    implicit none
    !> Orbital number
    integer, intent(in) :: orb_idx
    !> Number of atoms
    integer, intent(in) :: natoms
    !> Number of basis functions
    integer, intent(in) :: nbasis
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Overlap matrix to power 1/2 multiplied by MO coefficients: S^{1/2} C
    real(realk), dimension(nbasis,nbasis), intent(in) :: ShalfC
    !> Lowdin charges
    real(realk), dimension(natoms), intent(inout) :: charges
    integer, pointer :: atom_start(:) => null()
    integer, pointer :: atom_end(:) => null()
    integer :: atom,mu


    ! Init stuff
    charges=0.0E0_realk
    atom_start => MyMolecule%atom_start
    atom_end => MyMolecule%atom_end

    ! Loop over atoms
    do atom=1,natoms

       ! Lowdin charge for orbital I, on the atom "atom":
       !
       ! chargeI(atom) = sum_{mu \in atom} { [S^{1/2} C]_{mu,I} }^2
       !
       ! where C are MO coefficients, S^{1/2} is AO overlap matrix to the power 1/2,
       ! mu is an AO index and I is a MO index for the orbital in question.
       !
       ! It is seen that:
       ! sum_{all atoms} chargeI(atom) = 1

       do mu=atom_start(atom),atom_end(atom) ! nu \in atom
          charges(atom)=charges(atom) + ShalfC(mu,orb_idx)**2
       end do

    end do

    atom_start => null()
    atom_end => null()

  end subroutine GetLowdinVector



  !> \brief Lowdin population analysis: Get Lowdin charges on all orbitals for a given orbital
  !> \author Kasper Kristensen
  !> \date September 2011
  subroutine GetLowdinVector_orbitals(orb_idx,nbasis,ShalfC,charges)

    implicit none
    !> Orbital number
    integer, intent(in) :: orb_idx
    !> Number of basis functions
    integer, intent(in) :: nbasis
    !> Overlap matrix to power 1/2 multiplied by MO coefficients: S^{1/2} C
    real(realk), dimension(nbasis,nbasis), intent(in) :: ShalfC
    !> Lowdin charges
    real(realk), dimension(nbasis), intent(inout) :: charges
    integer :: mu


    ! Loop over orbitals
    do mu=1,nbasis

       ! Lowdin charge for AO orbital mu: { [S^{1/2} C]_{mu,I} }^2
       !
       ! where C are MO coefficients, S^{1/2} is AO overlap matrix to the power 1/2,
       ! mu is an AO index and I is a MO index for the orbital in question.
       charges(mu) = ShalfC(mu,orb_idx)**2

    end do


  end subroutine GetLowdinVector_orbitals



  !> \brief Atomic norms for expansion coefficients
  subroutine GetOrbitalAtomicNorm(num,MyMolecule,natoms,atomic_norms)

    implicit none
    integer, intent(in) :: num
    type(fullmolecule), intent(in) :: MyMolecule
    integer, intent(in) :: natoms
    real(realk), dimension(natoms), intent(inout) :: atomic_norms

    integer, pointer :: atom_start(:) => null()
    integer, pointer :: atom_end(:) => null()
    real(realk), pointer :: C(:,:)
    integer :: i,nbasis,nocc,nvirt

    ! Init stuff
    nbasis=MyMolecule%nbasis
    nocc=MyMolecule%nocc
    nvirt=MyMolecule%nunocc
    call mem_alloc(C,nbasis,nbasis)
    C(1:nbasis,1:nocc) = MyMolecule%Co(1:nbasis,1:nocc)
    C(1:nbasis,nocc+1:nbasis) = MyMolecule%Cv(1:nbasis,1:nvirt)
    atom_start => MyMolecule%atom_start
    atom_end => MyMolecule%atom_end

    atomic_norms=0.0E0_realk
    do i=1,natoms
       atomic_norms(i)=dot_product(C(atom_start(i):atom_end(i),num), &
            C(atom_start(i):atom_end(i),num))
    end do

    atom_start => null()
    atom_end => null()
    call mem_dealloc(C)

  end subroutine GetOrbitalAtomicNorm

  !> \brief Print some info about orbitals
  subroutine PrintOrbitalsInfo(orbitals,num_orbitals,lu_output)

    implicit none
    integer, intent(in) :: num_orbitals,lu_output
    type(decorbital), dimension(num_orbitals), intent(in) :: orbitals
    integer :: i,j

    write(lu_output,'(/,a)') 'Orbital Atom #Atoms List of atoms'
    do i=1,num_orbitals
       write(lu_output,*) orbitals(i)%orbitalnumber, &
            orbitals(i)%centralatom, orbitals(i)%numberofatoms, ' -> '
       do j=1,orbitals(i)%numberofatoms
          write(lu_output,*) orbitals(i)%atoms(j)
       end do
       write(lu_output,'(a)') ''
    end do

    return
  end subroutine PrintOrbitalsInfo


  !> \brief Write orbital to file
  subroutine orbital_write(orb,iunit)
    implicit none
    type(decorbital), intent(in) :: orb
    integer, intent(in) :: iunit
    integer(kind=8) :: orbitalnumber64,centralatom64,numberofatoms64,atoms64(orb%numberofatoms)

    orbitalnumber64 = orb%orbitalnumber
    centralatom64   = orb%centralatom
    numberofatoms64 = orb%numberofatoms
    atoms64         = orb%atoms

    write(iunit) orbitalnumber64
    write(iunit) centralatom64
    write(iunit) numberofatoms64
    write(iunit) atoms64

  end subroutine orbital_write

  !> \brief Read orbital from file
  function orbital_read(iunit) result(orb)
    implicit none
    type(decorbital) :: orb
    integer, intent(in) :: iunit
    integer :: i
    integer(kind=8) :: orbitalnumber64,centralatom64,numberofatoms64
    integer(kind=4) :: orbitalnumber32,centralatom32,numberofatoms32
    integer(kind=8),pointer :: atoms64(:)
    integer(kind=4),pointer :: atoms32(:)

    !ConvertFrom64Bit: if(DECinfo%convert64to32) then
       read(iunit) orbitalnumber64
       read(iunit) centralatom64
       read(iunit) numberofatoms64
       call mem_alloc(atoms64,numberofatoms64)
       read(iunit) atoms64

       orb%orbitalnumber = int(orbitalnumber64)
       orb%centralatom   = int(centralatom64)
       orb%numberofatoms = int(numberofatoms64)
       call mem_alloc(orb%atoms,orb%numberofatoms)
       do i=1,orb%numberofatoms
          orb%atoms(i) = int(atoms64(i))
       end do
       call mem_dealloc(atoms64)

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
    integer :: nbasis,nocc,nunocc, nstart, nend,nunocc_end,i,j
    integer :: idx,a,b,funit
    real(realk) :: nocc_calc, nunocc_calc, tmp
    real(realk) :: can_norm, lcm_norm
    real(realk),pointer :: occ_vector(:), unocc_vector(:)

    ! Sanity check: This is just a debugging routine and it only works for dense matrix type
    if(matrix_type/=mtype_dense) then
       call lsquit('check_lcm_against_canonical: Only implemented for dense matrix type!',-1)
    end if


    ! Initialize stuff
    ! ****************
    nbasis = MyMolecule%nbasis
    nocc = MyMolecule%nocc
    nunocc = MyMolecule%nunocc
    call mat_init(Ccan,nbasis,nbasis)
    call mat_init(F,nbasis,nbasis)
    call mat_init(S,nbasis,nbasis)
    call mem_alloc(eival,nbasis)


    ! Set Fock and overlap matrices
    ! *****************************
    call mat_set_from_full(MyMolecule%fock(1:nbasis,1:nbasis), 1E0_realk,F)
    call mat_set_from_full(MyMolecule%overlap(1:nbasis,1:nbasis), 1E0_realk,S)


    ! Get canonical orbitals, sorted in order of increasing orbital energies
    ! **********************************************************************
    call mat_diag_f(F,S,eival,Ccan)



    ! Occupied canonical orbitals
    ! ***************************

    ! These are the first nbasis*nocc elements in Ccan
    call mat_init(Cocc_can,nbasis,nocc)
    nstart=1
    nend = nbasis*nocc
    Cocc_can%elms(nstart:nend) = Ccan%elms(nstart:nend)

    ! Unoccupied canonical orbitals
    ! *****************************
    ! These are the last nbasis*nunocc elements,
    ! i.e. from element nbasis*nocc+1 to nbasis*nbasis
    call mat_init(Cvirt_can,nbasis,nunocc)
    nstart = nbasis*nocc + 1
    nunocc_end = nbasis*nunocc
    nend = nbasis*nbasis
    Cvirt_can%elms(1:nunocc_end) = Ccan%elms(nstart:nend)



    ! Set LCM orbitals
    ! ****************
    call mat_init(Cocc_lcm,nbasis,nocc)
    call mat_init(Cvirt_lcm,nbasis,nunocc)
    call mat_set_from_full(MyMolecule%Co(1:nbasis,1:nocc), 1E0_realk,Cocc_lcm)
    call mat_set_from_full(MyMolecule%Cv(1:nbasis,1:nunocc), 1E0_realk,Cvirt_lcm)


    ! Construct canonical/LCM overlap matrix for occupied space
    ! *********************************************************
    ! SCLocc = Cocc_can^T S Cocc_lcm
    call mat_init(SCLocc,nocc,nocc)
    call util_AO_to_MO_different_trans(Cocc_can, S, Cocc_lcm, SCLocc)


    ! Construct canonical/LCM overlap matrix for virtual space
    ! ********************************************************
    ! SCLvirt = Cvirt_can^T S Cvirt_lcm
    call mat_init(SCLvirt,nunocc,nunocc)
    call util_AO_to_MO_different_trans(Cvirt_can, S, Cvirt_lcm, SCLvirt)


    ! Check orthonormality
    ! ********************

    ! Unit matrices
    call mat_init(UnitMatrix,nbasis,nbasis)
    call mat_zero(UnitMatrix)
    do i=1,nbasis
       idx = get_matrix_position(i,i,nbasis,nbasis)
       UnitMatrix%elms(idx) = 1E0_realk
    end do

    ! Overlap for canonical orbitals
    ! ''''''''''''''''''''''''''''''
    call mat_init(tmp1,nbasis,nbasis)
    call util_AO_to_MO_different_trans(Ccan, S, Ccan, tmp1)
    ! Subtract unit matrix from occupied-occupied overlap
    call mat_daxpy(-1E0_realk,UnitMatrix,tmp1)
    ! Find norm of overlap minus unit matrix (this should be zero)
    can_norm = mat_sqnorm2(tmp1)
    can_norm = sqrt(can_norm)

    ! Overlap for LCM orbitals
    ! ''''''''''''''''''''''''
    call mat_init(Clcm,nbasis,nbasis)
    ! Set occupied orbitals
    nstart=1
    nend = nbasis*nocc
    Clcm%elms(nstart:nend) = Cocc_lcm%elms(nstart:nend)
    ! Set virtual orbitals
    nstart = nbasis*nocc + 1
    nend = nbasis*nbasis
    nunocc_end = nbasis*nunocc
    Clcm%elms(nstart:nend) = Cvirt_lcm%elms(1:nunocc_end)
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


    ! Create occupied and unoccupied vectors
    ! **************************************

    ! The occupied vector contains each occupied canonical orbital projected
    ! against the full set of occupied LCM orbitals:
    ! occupied_vector(i) = sum_j SCLocc(i,j)**2
    ! and similarly for the virtual vector.
    ! Thus, each element in the occupied_vector and the unoccupied_vector
    ! has to equal one. And the sum of elements equals nocc/nunocc for
    ! occupied/unoccupied vectors.

    call mem_alloc(occ_vector,nocc)
    call mem_alloc(unocc_vector,nunocc)
    occ_vector=0E0_realk
    unocc_vector=0E0_realk

    ! Occupied vector
    ! '''''''''''''''
    do i=1,nocc
       do j=1,nocc
          idx = get_matrix_position(i,j,nocc,nocc)
          tmp = SCLocc%elms(idx)
          occ_vector(i) = occ_vector(i) + tmp**2
       end do
    end do

    ! Total number of occupied orbitals
    nocc_calc=0E0_realk
    do i=1,nocc
       nocc_calc = nocc_calc + occ_vector(i)
    end do


    ! Virtual vector
    ! ''''''''''''''
    do a=1,nunocc
       do b=1,nunocc
          idx = get_matrix_position(a,b,nunocc,nunocc)
          tmp = SCLvirt%elms(idx)
          unocc_vector(a) = unocc_vector(a) + tmp**2
       end do
    end do

    ! Total number of virtual orbitals
    nunocc_calc=0E0_realk
    do a=1,nunocc
       nunocc_calc = nunocc_calc + unocc_vector(a)
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
    do i=1,nocc
       write(DECinfo%output,'(3X,i8,g18.8)') i, occ_vector(i)
    end do
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(3X,a)') '  Virtual projections   '
    write(DECinfo%output,'(3X,a)') '------------------------'
    do a=1,nunocc
       write(DECinfo%output,'(3X,i8,g18.8)') a, unocc_vector(a)
    end do
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(3X,a,g18.8)') 'Calculated number of occ. orbitals  :', nocc_calc
    write(DECinfo%output,'(3X,a,g18.8)') 'Calculated number of virt. orbitals :', nunocc_calc
    write(DECinfo%output,'(3X,a,g18.8)') 'Exact number of occ. orbitals       :', real(nocc)
    write(DECinfo%output,'(3X,a,g18.8)') 'Exact number of virt. orbitals      :', real(nunocc)
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
    call mem_dealloc(unocc_vector)

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
    nvirt = MyMolecule%nunocc
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
         & nOrb,MyMolecule%Co,Cocc1,CentralAtom,MyMolecule%SubSystemIndex)      

    !MOs assigned to SubSystem 2 and AOs belonging to Subsystem 1
    call mem_alloc(Cocc2,nbasisSub(1),noccsub(2))
    call BuildSubsystemCMO(2,nbasisSub(1),noccsub(2),nocc,nAtoms,nbasis,&
         & nOrb,MyMolecule%Co,Cocc2,CentralAtom,MyMolecule%SubSystemIndex)      

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
    call ls_output(MyMolecule%Co,1,nbasis,1,nocc,nbasis,nocc,1,DECinfo%output)

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
    nvirt = MyMolecule%nunocc
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
         & nOrb,MyMolecule%Co,CentralAtom,MyMolecule%SubSystemIndex)      

    call ZeroSubsystemCMO(2,nbasisSub(1),noccsub(2),nocc,nAtoms,nbasis,&
         & nOrb,MyMolecule%Co,CentralAtom,MyMolecule%SubSystemIndex)      

    WRITE(DECinfo%output,*)'CMO after force_Occupied_SubSystemLocality'
    call ls_output(MyMolecule%Co,1,nbasis,1,nocc,nbasis,nocc,1,DECinfo%output)

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
  !> with occupied orbitals before unoccupied orbitals.
  !> \author Kasper Kristensen
  !> \date January 2011
  subroutine write_DECorbitals_to_file(nocc,nunocc,&
       &OccOrbitals,UnoccOrbitals)

    implicit none
    !> Number of occupied orbitals in full molecule
    integer, intent(in) :: nocc
    !> Number of unoccupied orbitals in full molecule
    integer, intent(in) :: nunocc
    !> Occupied orbital info for full molecule
    type(decorbital), dimension(nocc), intent(in) :: OccOrbitals
    !> Unoccupied orbital info for full molecule
    type(decorbital), dimension(nunocc), intent(in) :: UnoccOrbitals
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


    ! Write unoccupied orbitals to file
    ! *********************************
    do i=1,nunocc
       call orbital_write(UnoccOrbitals(i),funit)
    end do


    ! Close file
    ! **********
    call lsclose(funit,'KEEP')


  end subroutine write_DECorbitals_to_file



  !> \brief Read all orbitals from file "DECorbitals.info"
  !> with occupied orbitals before unoccupied orbitals.
  !> \author Kasper Kristensen
  !> \date January 2011
  subroutine read_DECorbitals_from_file(nocc,nunocc,&
       &OccOrbitals,UnoccOrbitals)

    implicit none
    !> Number of occupied orbitals in full molecule
    integer, intent(in) :: nocc
    !> Number of unoccupied orbitals in full molecule
    integer, intent(in) :: nunocc
    !> Occupied orbital info for full molecule
    type(decorbital), dimension(nocc), intent(inout) :: OccOrbitals
    !> Unoccupied orbital info for full molecule
    type(decorbital), dimension(nunocc), intent(inout) :: UnoccOrbitals
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


    ! Read unoccupied orbitals from file
    ! **********************************
    do i=1,nunocc
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

          ! CHECK 1: Is neighbor is included in orbital extent?
          ! (This is not crucial, but nice to know for statistics).
          ! Therefore, even if neighbor is not in orbital extent, we reassign anyway.
          included=.false.
          do j=1,Orbitals(i)%numberofatoms
             if(Orbitals(i)%atoms(j) == neighbor) then
                included=.true.
                exit
             end if
          end do
          if(.not. included) then
             if(DECinfo%PL>0) then
                write(DECinfo%output,*)
                write(DECinfo%output,*) 'WARNING in reassign_orbitals: &
                     &Atom neighbouring hydrogen is not included in the orbital extent'
                write(DECinfo%output,*) 'Orbital is reassigned anyway...'
                write(DECinfo%output,*) 'Orbital number =', i
                write(DECinfo%output,*) 'Central atom =', centralatom
                write(DECinfo%output,*) 'Neighbor =', neighbor
                write(DECinfo%output,*) 'Orbital extent:', Orbitals(i)%atoms(:)
                write(DECinfo%output,'(1X,a,F12.3)') 'Distance to neighbor (Angstrom)    =', &
                     & SortedDistTable(2,centralatom)*bohr
             end if
          end if

          ! CHECK 2: Neighbor is not hydrogen
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
  subroutine adjust_orbitals_for_full_simulation(nocc,nunocc,&
       &OccOrbitals,UnoccOrbitals,natoms,n)

    implicit none
    !> Number of occupied orbitals in full molecule
    integer, intent(in) :: nocc
    !> Number of unoccupied orbitals in full molecule
    integer, intent(in) :: nunocc
    !> Occupied orbital info for full molecule
    type(decorbital), dimension(nocc), intent(inout) :: OccOrbitals
    !> Unoccupied orbital info for full molecule
    type(decorbital), dimension(nunocc), intent(inout) :: UnoccOrbitals
    !> Number of atoms in the molecule
    integer, intent(in) :: natoms
    !> Number of atomic sites with orbitals assigned in the simulation (default: 1)
    !> Thus, the first n atoms will get all orbitals evenly assigned.
    integer,intent(in) :: n
    integer :: i,j,orb,nocc_per_atom,nunocc_per_atom,atom

    ! Number of occupied/virtual orbitals per atom evenly distributed over n atoms
    nocc_per_atom = floor(real(nocc)/real(n))
    nunocc_per_atom = floor(real(nunocc)/real(n))


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
          call mem_alloc(OccOrbitals(orb)%atoms,natoms)

          ! Orbital extent is the full molecule
          OccOrbitals(orb)%numberofatoms=natoms
          do j=1,natoms
             OccOrbitals(orb)%atoms(j) = j
          end do

          ! Central atom for this occupied orbital
          OccOrbitals(orb)%centralatom = atom

       end do
    end do
    i=orb+1

    ! Assign remaining occupied orbitals (if any) to the last atom n
    do orb=i,nocc
       call orbital_free(OccOrbitals(orb))
       call mem_alloc(OccOrbitals(orb)%atoms,natoms)
       OccOrbitals(orb)%numberofatoms=natoms
       do j=1,natoms
          OccOrbitals(orb)%atoms(j) = j
       end do
       OccOrbitals(orb)%centralatom = n
    end do


    ! Assign unoccupied orbitals
    ! **************************

    orb=0
    do atom=1,n  ! loop over the n first atoms
       do i=1,nunocc_per_atom  ! loop over nunocc_per_atom orbitals per atom

          ! Increase orbital counter
          orb=orb+1

          ! Free existing orbital info
          call orbital_free(UnoccOrbitals(orb))
          call mem_alloc(UnoccOrbitals(orb)%atoms,natoms)

          ! Orbital extent is the full molecule
          UnoccOrbitals(orb)%numberofatoms=natoms
          do j=1,natoms
             UnoccOrbitals(orb)%atoms(j) = j
          end do

          ! Central atom for this unoccupied orbital
          UnoccOrbitals(orb)%centralatom = atom

       end do
    end do
    i=orb+1

    ! Assign remaining unoccupied orbitals (if any) to the last atom n
    do orb=i,nunocc
       call orbital_free(UnoccOrbitals(orb))
       call mem_alloc(UnoccOrbitals(orb)%atoms,natoms)
       UnoccOrbitals(orb)%numberofatoms=natoms
       do j=1,natoms
          UnoccOrbitals(orb)%atoms(j) = j
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


  !> \brief Print number of occupied and unoccupied orbitals assigned to each atom for
  !> all atoms in the molecule.
  !> \author Kasper Kristensen
  subroutine print_orbital_info(mylsitem,nocc,natoms,nunocc,OccOrbitals,&
       & UnoccOrbitals,ncore)

    implicit none
    !> Dalton LSITEM (just for printing atom type)
    type(lsitem), intent(inout) :: mylsitem
    !> Number of atoms
    integer, intent(in) :: natoms
    !> Number of occupied orbitals
    integer, intent(in) :: nocc
    !> Number of unoccupied orbitals
    integer, intent(in) :: nunocc
    !> Occupied orbitals
    type(decorbital), intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals
    type(decorbital), intent(in) :: UnoccOrbitals(nunocc)
    !> Number of of core orbitals. If present and frozen core approx is used,
    !> the first ncore orbitals will not be printed.
    integer,intent(in),optional :: ncore
    integer :: nocc_per_atom(natoms), nunocc_per_atom(natoms)
    integer :: SECnocc_per_atom(natoms), SECnunocc_per_atom(natoms)
    integer :: i, occ_max_orbital_extent, unocc_max_orbital_extent, occ_idx, unocc_idx,j,offset
    real(realk) :: occ_av_orbital_extent, unocc_av_orbital_extent

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
    call get_orbital_extent_info(nunocc,UnoccOrbitals,unocc_max_orbital_extent,&
         & unocc_av_orbital_extent, unocc_idx)


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
    nunocc_per_atom =  get_number_of_orbitals_per_atom(UnoccOrbitals,nunocc,natoms,.true.)
    ! Secondary assigning
    SECnunocc_per_atom =  get_number_of_orbitals_per_atom(UnoccOrbitals,nunocc,natoms,.false.)

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
            & nocc_per_atom(i), nunocc_per_atom(i), SECnocc_per_atom(i), SECnunocc_per_atom(i)
    end do
    write(DECinfo%output,'(1X,A,11X,I6,11X,I6)') 'Total:', nocc-offset, nunocc
    write(DECinfo%output,*)

    ! Print out orbital extent summary
    ! ********************************
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) 'ORBITAL EXTENT INFORMATION (NUMBER OF ATOMS USED TO SPAN EACH ORBITAL)'
    write(DECinfo%output,*) '**********************************************************************'
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a,i6,a,i6)') 'Maximum occ orbital extent   = ', occ_max_orbital_extent, &
         & '   -- Orbital index', occ_idx
    write(DECinfo%output,'(1X,a,i6,a,i6)') 'Maximum unocc orbital extent = ', unocc_max_orbital_extent, &
         & '   -- Orbital index', unocc_idx
    write(DECinfo%output,'(1X,a,f12.4)') 'Average occ orbital extent   = ', occ_av_orbital_extent
    write(DECinfo%output,'(1X,a,f12.4)') 'Average unocc orbital extent = ', unocc_av_orbital_extent
    write(DECinfo%output,*)
    write(DECinfo%output,*)

    if(DECinfo%PL>0) then ! print specific info for each orbital
       do i=1,nocc
          write(DECinfo%output,*) 'Occupied orbital: ', i
          write(DECinfo%output,*) '************************************************'
          write(DECinfo%output,*) 'Central atom: ', OccOrbitals(i)%centralatom
          write(DECinfo%output,*) '# atoms in orbital extent: ', OccOrbitals(i)%numberofatoms
          write(DECinfo%output,*) 'Atoms in orbital extent: '
          do j=1,OccOrbitals(i)%numberofatoms
             write(DECinfo%output,*) OccOrbitals(i)%atoms(j)
          end do
          write(DECinfo%output,*)
       end do
       write(DECinfo%output,*)
       write(DECinfo%output,*)

       do i=1,nunocc
          write(DECinfo%output,*) 'Virtual orbital: ', i
          write(DECinfo%output,*) '************************************************'
          write(DECinfo%output,*) 'Central atom: ', UnOccOrbitals(i)%centralatom
          write(DECinfo%output,*) '# atoms in orbital extent: ', UnOccOrbitals(i)%numberofatoms
          write(DECinfo%output,*) 'Atoms in orbital extent: ' 
          do j=1,UnOccOrbitals(i)%numberofatoms
             write(DECinfo%output,*) UnOccOrbitals(i)%atoms(j)
          end do
          write(DECinfo%output,*)
       end do
       write(DECinfo%output,*)
       write(DECinfo%output,*)

    end if


  end subroutine print_orbital_info






  !> \brief Print number of occupied and unoccupied orbitals assigned to each atom for
  !> all atoms in the molecule.
  !> \author Kasper Kristensen
  subroutine print_orbital_info_DECCO(mylsitem,nocc,nunocc,OccOrbitals,&
       & UnoccOrbitals,ncore)

    implicit none
    !> Dalton LSITEM (just for printing atom type)
    type(lsitem), intent(inout) :: mylsitem
    !> Number of occupied orbitals
    integer, intent(in) :: nocc
    !> Number of unoccupied orbitals
    integer, intent(in) :: nunocc
    !> Occupied orbitals
    type(decorbital), intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals
    type(decorbital), intent(in) :: UnoccOrbitals(nunocc)
    !> Number of of core orbitals. If present and frozen core approx is used,
    !> the first ncore orbitals will not be printed.
    integer,intent(in) :: ncore
    integer :: nunocc_per_occ(nocc)
    integer :: i, occ_max_orbital_extent, unocc_max_orbital_extent, occ_idx, unocc_idx,j
    real(realk) :: occ_av_orbital_extent, unocc_av_orbital_extent


    ! Find average and maximum orbital extent
    ! ***************************************
    ! Occupied
    call get_orbital_extent_info(nocc,OccOrbitals,occ_max_orbital_extent,&
         & occ_av_orbital_extent, occ_idx)

    ! Unoccupied
    call get_orbital_extent_info(nunocc,UnoccOrbitals,unocc_max_orbital_extent,&
         & unocc_av_orbital_extent, unocc_idx)

    ! Number of unoccupied orbitals assigned to each occupied orbital
    call DECCO_get_nvirt_per_occ_orbital(nocc,nunocc,UnoccOrbitals,nunocc_per_occ)


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
       write(DECinfo%output,'(1X,I5,12X,I10)') i,nunocc_per_occ(i)
    end do
    write(DECinfo%output,'(1X,A,11X,I6,11X,I6)') 'Total:', nocc-ncore, nunocc
    write(DECinfo%output,*)


    ! Print out orbital extent summary
    ! ********************************
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) 'ORBITAL EXTENT INFORMATION (NUMBER OF ATOMS USED TO SPAN EACH ORBITAL)'
    write(DECinfo%output,*) '**********************************************************************'
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a,i6,a,i6)') 'Maximum occ orbital extent   = ', occ_max_orbital_extent, &
         & '   -- Orbital index', occ_idx
    write(DECinfo%output,'(1X,a,i6,a,i6)') 'Maximum unocc orbital extent = ', unocc_max_orbital_extent, &
         & '   -- Orbital index', unocc_idx
    write(DECinfo%output,'(1X,a,f12.4)') 'Average occ orbital extent   = ', occ_av_orbital_extent
    write(DECinfo%output,'(1X,a,f12.4)') 'Average unocc orbital extent = ', unocc_av_orbital_extent
    write(DECinfo%output,*)
    write(DECinfo%output,*)


    ! print specific info for each orbital
    if(DECinfo%PL>0) then 
       do i=1,nocc
          write(DECinfo%output,*) 'Occupied orbital: ', i
          write(DECinfo%output,*) '************************************************'
          write(DECinfo%output,*) '# atoms in orbital extent: ', OccOrbitals(i)%numberofatoms
          write(DECinfo%output,*) 'Atoms in orbital extent: '
          do j=1,OccOrbitals(i)%numberofatoms
             write(DECinfo%output,*) OccOrbitals(i)%atoms(j)
          end do
          write(DECinfo%output,*)
       end do
       write(DECinfo%output,*)
       write(DECinfo%output,*)

       do i=1,nunocc
          write(DECinfo%output,*) 'Virtual orbital: ', i
          write(DECinfo%output,*) '************************************************'
          write(DECinfo%output,*) '# atoms in orbital extent: ', UnOccOrbitals(i)%numberofatoms
          write(DECinfo%output,*) 'Atoms in orbital extent: ' 
          do j=1,UnOccOrbitals(i)%numberofatoms
             write(DECinfo%output,*) UnOccOrbitals(i)%atoms(j)
          end do
          write(DECinfo%output,*)
       end do
       write(DECinfo%output,*)
       write(DECinfo%output,*)

    end if


  end subroutine print_orbital_info_DECCO


  !> Get number of virtual orbital assigned to each occupied orbital
  subroutine DECCO_get_nvirt_per_occ_orbital(nocc,nunocc,UnoccOrbitals,nunocc_per_occ)
    implicit none
    !> Number of occupied orbitals
    integer, intent(in) :: nocc
    !> Number of unoccupied orbitals
    integer, intent(in) :: nunocc
    !> Unoccupied orbitals
    type(decorbital), intent(in) :: UnoccOrbitals(nunocc)
    !> Number of unoccupied orbitals assigned to each occupied orbital
    integer, intent(inout) :: nunocc_per_occ(nocc)
    integer :: i,j

    nunocc_per_occ=0
    do j=1,nunocc
       nunocc_per_occ(UnoccOrbitals(j)%centralatom) = nunocc_per_occ(UnoccOrbitals(j)%centralatom) + 1
    end do


  end subroutine DECCO_get_nvirt_per_occ_orbital

  !> \brief The DEC scheme only works if for each atom either zero occupied AND zero virtual
  !> orbitals are assigned - or if nonzero occupied AND nonzero virtual orbitals are assigned.
  !> If this is not the case, the system under consideration is presumably a debug molecule
  !> and we quit here, rather than encountering uninitialized pointers later on...
  !> We also check whether it is necessary to repeat the atomic fragment calcs
  !> after the fragment optimization and set DECinfo%RepeatAF accordingly.
  !> \author Kasper Kristensen
  !> \date December 2011
  subroutine dec_orbital_sanity_check(natoms,nocc,nunocc,OccOrbitals,&
       & UnoccOrbitals,MyMolecule)

    ! NOTE!!! 
    ! DECinfo%RepeatAF is also intent(inout) here, although it is part of a global structure.
    ! I know this is not pretty but this solution will have to do for now...

    implicit none
    !> Number of atoms
    integer, intent(in) :: natoms
    !> Number of occupied orbitals
    integer, intent(in) :: nocc
    !> Number of unoccupied orbitals
    integer, intent(in) :: nunocc
    !> Occupied orbitals
    type(decorbital), intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals
    type(decorbital), intent(in) :: UnoccOrbitals(nunocc)
    !> Full molecule info
    type(fullmolecule),intent(in) :: MyMolecule
    integer :: nocc_per_atom(natoms), nunocc_per_atom(natoms)
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
    nunocc_per_atom =  get_number_of_orbitals_per_atom(UnoccOrbitals,nunocc,natoms,.true.)


    something_wrong=.false.
    nfrags=0
    do i=1,natoms

       if(.not.(DECinfo%onlyoccpart.or.DECinfo%onlyvirtpart)&
            & .and.(.NOT.decinfo%PureHydrogendebug))then
          if( (nocc_per_atom(i) == 0) .and. (nunocc_per_atom(i)/=0) ) something_wrong=.true.
          if( (nocc_per_atom(i) /= 0) .and. (nunocc_per_atom(i)==0) ) something_wrong=.true.
          if(something_wrong) then
             write(DECinfo%output,*) 'Atom = ',i
             write(DECinfo%output,*) 'Number of occupied orbitals   assigned = ', nocc_per_atom(i)
             write(DECinfo%output,*) 'Number of unoccupied orbitals assigned = ', nunocc_per_atom(i)
             write(DECinfo%output,*) 'If you use Phantom Atoms try .ONLYOCCPART keyword'
             print*,'If you use Phantom Atoms try .ONLYOCCPART keyword'
             call lsquit('Orbital assigment is inconsistent &
                  & with DEC scheme',DECinfo%output)
          end if
       end if

       ! Count number of atomic fragments
       if( (nocc_per_atom(i) /= 0) .and. (nunocc_per_atom(i)/=0) ) then
          nfrags=nfrags+1
       end if
    end do


    ! Secondary assignment check
    ! **************************
    nocc_per_atom =  get_number_of_orbitals_per_atom(OccOrbitals,nocc,natoms,.false.,&
         & offset=offset)
    nunocc_per_atom =  get_number_of_orbitals_per_atom(UnoccOrbitals,nunocc,natoms,.false.)
    something_wrong=.false.
    do i=1,natoms
       if( (nocc_per_atom(i) == 0) .and. (nunocc_per_atom(i)/=0) ) something_wrong=.true.
       if( (nocc_per_atom(i) /= 0) .and. (nunocc_per_atom(i)==0) ) something_wrong=.true.
       if(something_wrong) then
          write(DECinfo%output,*) 'Atom = ',i
          write(DECinfo%output,*) 'Number of occupied orbitals   assigned = ', nocc_per_atom(i)
          write(DECinfo%output,*) 'Number of unoccupied orbitals assigned = ', nunocc_per_atom(i)
          call lsquit('Secondary orbital assigment is inconsistent &
               & with DEC scheme',DECinfo%output)
       end if
    end do




    ! Repeat atomic fragment calcs after fragment optimization???
    ! -----------------------------------------------------------
    ! Cases where it is necessary to repeat atomic fragment calcs:
    ! - first order properties are requested
    ! - only one fragment 
    ! - fragment opt where reduction step is done for a different model than the target CC model.
    if(DECinfo%first_order .or. nfrags==1 .or. DECinfo%InclFullMolecule &
         & .or. (DECinfo%ccmodel/=DECinfo%fragopt_red_model ) ) then
       DECinfo%RepeatAF=.true.
    else
       DECinfo%RepeatAF=.false.
    end if
    IF(DECinfo%FragmentExpansionRI)THEN
       DECinfo%RepeatAF=.true.
    ENDIF
    IF(DECinfo%InteractionEnergy)THEN
       IF(DECinfo%ccmodel.NE.DECinfo%fragopt_red_model)THEN
          !turn of recalculation of Atomic fragment calculations
          DECinfo%RepeatAF = .FALSE.
       ENDIF
    ENDIF

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
       if(Orbitals(i)%numberofatoms > max_orbital_extent) then
          max_orbital_extent = Orbitals(i)%numberofatoms
          orb_index = i
       end if

       ! Average
       av_orbital_extent = av_orbital_extent + real(Orbitals(i)%numberofatoms)

    end do

    ! Average orbital extent
    av_orbital_extent = av_orbital_extent / real(norb)


  end subroutine get_orbital_extent_info



  !> \brief Determine which fragments to consider
  !> For atom-based DEC - the atoms with orbitals assigned
  !> For DECCO - all occupied orbitals (but not the core orbitals if frozencore is used)
  !> \author Kasper Kristensen
  !> \date October 2013
  subroutine which_fragments_to_consider(ncore,nocc,nunocc,natoms,&
       & OccOrbitals,UnoccOrbitals,dofrag,PhantomAtom)
    !> Number of core orbitals in full molecule
    integer,intent(in) :: ncore
    !> Number of occupied orbitals in full molecule
    integer,intent(in) :: nOcc
    !> Number of unoccupied orbitals in full molecule
    integer,intent(in) :: nUnocc
    !> Number of atoms in full molecule
    integer,intent(in) :: nAtoms
    !> Occupied orbitals in DEC format
    type(decorbital), intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals in DEC format
    type(decorbital), intent(in) :: UnoccOrbitals(nunocc)
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
       call which_atoms_have_orbitals_assigned(ncore,nocc,nunocc,natoms,&
            & OccOrbitals,UnoccOrbitals,dofrag,PhantomAtom)
    end if

  end subroutine which_fragments_to_consider


  !> \brief Determine which atoms have one or more orbitals assigned.
  !> \author Kasper Kristensen
  !> \date October 2013
  subroutine which_atoms_have_orbitals_assigned(ncore,nocc,nunocc,natoms,&
       & OccOrbitals,UnoccOrbitals,dofrag,PhantomAtom)
    implicit none
    !> Number of core orbitals in full molecule
    integer,intent(in) :: ncore
    !> Number of occupied orbitals in full molecule
    integer,intent(in) :: nOcc
    !> Number of unoccupied orbitals in full molecule
    integer,intent(in) :: nUnocc
    !> Number of atoms in full molecule
    integer,intent(in) :: nAtoms
    !> Occupied orbitals in DEC format
    type(decorbital), intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals in DEC format
    type(decorbital), intent(in) :: UnoccOrbitals(nunocc)
    !> dofrag(P) is true/false if atom P has one or more/zero orbitals assigned.
    logical,intent(inout) :: dofrag(natoms)
    !> Which atoms are Phantom Atoms
    logical, intent(in) :: PhantomAtom(natoms)
    !local 
    integer, dimension(natoms) :: nocc_per_atom, nunocc_per_atom
    integer :: i

    ! Number of orbitals per atom
    if(DECinfo%frozencore) then ! only count valence orbitals
       nocc_per_atom =  get_number_of_orbitals_per_atom(OccOrbitals,nocc,natoms,.true.,offset=ncore)
    else
       nocc_per_atom =  get_number_of_orbitals_per_atom(OccOrbitals,nocc,natoms,.true.)
    end if

    nunocc_per_atom =  get_number_of_orbitals_per_atom(UnOccOrbitals,nunocc,natoms,.true.)

    ! Which fragments to consider
    dofrag=.true.
    do i=1,natoms
       if(DECinfo%onlyoccpart) then
          ! Only consider occupied orbitals

          if( (nocc_per_atom(i)==0) ) then
             dofrag(i)=.false.
          else
             if(PhantomAtom(i))then
                print*,'ERROR   nocc_per_atom',nocc_per_atom(i),'nunocc_per_atom(i)',nunocc_per_atom(i)
                print*,'ERROR   i',i,'PhantomAtom',PhantomAtom(i)
                dofrag(i)=.false.
                print*,'Setting dofrag to false'
             endif
          endif
       elseif(DECinfo%onlyvirtpart) then
          ! Only consider virtual orbitals
          if( (nunocc_per_atom(i)==0) ) then
             dofrag(i)=.false.
          else
             if(PhantomAtom(i))then
                print*,'ERROR   nocc_per_atom',nocc_per_atom(i),'nunocc_per_atom(i)',nunocc_per_atom(i)
                print*,'ERROR   i',i,'PhantomAtom',PhantomAtom(i)
                dofrag(i)=.false.
                print*,'Setting dofrag to false'
             endif
          endif
       else
          ! Consider occupied as well as unoccupied orbitals
          if( (nocc_per_atom(i)==0) .and. (nunocc_per_atom(i)==0) )then
             dofrag(i)=.false.
          else
             if(PhantomAtom(i))then
                print*,'ERROR   nocc_per_atom',nocc_per_atom(i),'nunocc_per_atom(i)',nunocc_per_atom(i)
                print*,'ERROR   i',i,'PhantomAtom',PhantomAtom(i)
                dofrag(i)=.false.
                print*,'Setting dofrag to false'
             endif
          endif
       end if
    end do
    if( DECinfo%only_n_frag_jobs>0)then
       print *,"HACK TO ONLY DO ONE FRAGMENT JOB"
       dofrag = .false.
       do i = 1, DECinfo%only_n_frag_jobs
          dofrag(DECinfo%frag_job_nr(i)) = .true.
       enddo
    endif

  end subroutine which_atoms_have_orbitals_assigned

end module orbital_operations
