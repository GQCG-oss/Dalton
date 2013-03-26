!> @file
!> Main driver for DEC program. 

!> \author Marcin Ziolkowski (modified for Dalton by Kasper Kristensen)
!> \date 2010-09

module dec_main_mod

#define version 0
#define subversion 2

  use precision
  use lstiming!,only: lstimer
  use typedeftype!,only: lsitem
  use matrix_module!,only:matrix
  use matrix_operations !,only: mat_read_from_disk,mat_set_from_full
  use memory_handling!,only: mem_alloc, mem_dealloc
  use dec_typedef_module
  use files !,only:lsopen,lsclose


  ! DEC DEPENDENCIES (within deccc directory) 
  ! *****************************************
  use dec_fragment_utils
  use array3_memory_manager!,only: print_memory_currents_3d
  use array4_memory_manager!,only: print_memory_currents4
  use full_molecule!,only: molecule_init, molecule_finalize,molecule_init_from_inputs
  use orbital_operations!,only: check_lcm_against_canonical
  use array_operations!,only:test_array_struct,test_array_reorderings
  use full_molecule!,only: molecule_copyback_FSC_matrices
  use mp2_gradient_module!,only: dec_get_error_difference
  use dec_driver_module,only: dec_wrapper, main_fragment_driver
  use full,only: full_driver

public :: dec_main_prog, get_mp2gradient_and_energy_from_inputs, get_total_mp2energy_from_inputs
private

contains

  !> \brief Main DEC program.
  !> \author Marcin Ziolkowski (modified for Dalton by Kasper Kristensen)
  !> \date September 2010
  subroutine dec_main_prog(MyLsitem)

    implicit none
    type(lsitem), intent(inout) :: mylsitem
    type(fullmolecule) :: molecule
    character(len=10) :: program_version
    character(len=50) :: MyHostname
    integer, dimension(8) :: values
    real(realk) :: tcpu1, twall1, tcpu2, twall2
    type(matrix) :: D
    integer :: funit,dim1,dim2
    integer(8) :: longdim1,longdim2
    integer(kind=4) :: dim132,dim232
    logical :: onmaster
    real(realk),pointer :: Dfull(:,:)

    onmaster=.true.

    !Array test
    if (DECinfo%array_test)then
      print *,"TEST ARRAY MODULE"
      call test_array_struct()
      return
    endif
    if (DECinfo%reorder_test)then
      print *,"TEST REORDERINGS"
      call test_array_reorderings()
      return
    endif

    call LSTIMER('START',tcpu1,twall1,DECinfo%output)

    write(program_version,'("v:",i2,".",i2.2)') version,subversion
    ! =============================================================
    ! Main program 
    ! =============================================================
    write(DECinfo%output,*) 
    write(DECinfo%output,*)
    write(DECinfo%output,'(a)') '============================================================================='
    write(DECinfo%output,'(a,a)') &
         '             -- Divide, Expand & Consolidate Coupled-Cluster -- ',program_version
    write(DECinfo%output,'(a)') '============================================================================='
    write(DECinfo%output,'(a)') ' Authors: Marcin Ziolkowski @ AU 2009,2010' 
    write(DECinfo%output,'(a)') '          Kasper Kristensen (kasperk@chem.au.dk)' 
    write(DECinfo%output,'(a)') '          Ida-Marie Hoeyvik (idamh@chem.au.dk)' 
    write(DECinfo%output,'(a)') '          Patrick Ettenhuber (pett@chem.au.dk)'
    write(DECinfo%output,'(a)') '          Janus Juul Eriksen (janusje@chem.au.dk)'
    write(DECinfo%output,*)
    write(DECinfo%output,*)

    ! Set DEC memory
    call get_memory_for_dec_calculation()

    ! Get informations about full molecule
    call molecule_init(molecule,mylsitem,DECinfo%output,DECinfo%output,DECinfo%output)

    ! Get density matrix in type(matrix) form (more elegant solution is seeked)
    funit=-1
    call lsopen(funit,'dens.restart','OLD','UNFORMATTED')
    call mat_init(D,molecule%nbasis,molecule%nbasis)
    call mem_alloc(Dfull,D%nrow,D%ncol)
    if(DECinfo%convert64to32) then
       READ(funit) longdim1,longdim2
       dim1 = int(longdim1)
       dim2 = int(longdim2)
    elseif(DECinfo%convert32to64) then
       READ(funit) dim132,dim232
       dim1 = dim132
       dim2 = dim232
    else
       READ(funit) dim1,dim2
    end if
    if(dim1 /= D%nrow) stop 'DEC read dens: dim1 /= A%nrow'
    if(dim1 /= D%ncol) stop 'DEC read dens: dim2 /= A%ncol'
    read(funit) Dfull
    call mat_set_from_full(Dfull,1E0_realk,D)
    call mem_dealloc(Dfull)
    call lsclose(funit,'KEEP')

    ! Sanity check: LCM orbitals span the same space as canonical orbitals 
    if(DECinfo%check_lcm_orbitals) then
       call check_lcm_against_canonical(molecule,MyLsitem)
    end if

    if(DECinfo%full_molecular_cc) then
       ! -- Call full molecular CC
       write(DECinfo%output,'(/,a,/)') 'Full molecular calculation is carried out...'
       call full_driver(molecule,mylsitem,D)
       ! --
    else
       ! -- Initialize DEC driver for energy calculation
       write(DECinfo%output,'(/,a,/)') 'DEC fragment calculation is carried out...'
       call DEC_wrapper(molecule,mylsitem,D)
       ! --
    end if

    ! Finalize everything
    call molecule_finalize(molecule)
    call mat_free(D)
    write(DECinfo%output,'(a)') 'Full molecular data deleted'

    call LSTIMER('START',tcpu2,twall2,DECinfo%output)


    ! Print memory summary
    write(DECinfo%output,*)
    write(DECinfo%output,*) 'DEC memory summary'
    write(DECinfo%output,*) '------------------'
    call print_memory_currents4(DECinfo%output)
    write(DECinfo%output,*) '------------------'
    call print_memory_currents_3d(DECinfo%output)
    write(DECinfo%output,*) '------------------'
    write(DECinfo%output,*)


    write(DECinfo%output,*)
    write(DECinfo%output,'(/,a)') '------------------------------------------------------'
    write(DECinfo%output,'(a,g20.6,a)') 'Total CPU  time used in DEC           :',tcpu2-tcpu1,' s'
    write(DECinfo%output,'(a,g20.6,a)') 'Total Wall time used in DEC           :',twall2-twall1,' s'
    write(DECinfo%output,'(a,/)') '------------------------------------------------------'
    write(DECinfo%output,*)

#ifdef __GNUC__  
    call hostnm(MyHostname)
    write(DECinfo%output,'(a,a)')   'Hostname       : ',MyHostname
#endif
    call date_and_time(VALUES=values)
    write(DECinfo%output,'(2(a,i2.2),a,i4.4,3(a,i2.2))') &
         'Job finished   : Date: ',values(3),'/',values(2),'/',values(1), &
         '   Time: ',values(5),':',values(6),':',values(7)

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(/,a)') '============================================================================='
    write(DECinfo%output,'(a)')   '                          -- end of DEC program --                           '
    write(DECinfo%output,'(a,/)') '============================================================================='
    write(DECinfo%output,*)
    write(DECinfo%output,*)

    ! Update number of DEC calculations for given FOT level
    DECinfo%ncalc(DECinfo%FOTlevel) = DECinfo%ncalc(DECinfo%FOTlevel) +1
    call print_calculation_bookkeeping()

  end subroutine dec_main_prog


  !> \brief Calculate MP2 energy and molecular gradient using DEC scheme
  !> using the Fock, density, and MO matrices and input - rather than reading them from file
  !> as is done in a "conventional" DEC calculation which uses the dec_main_prog subroutine.
  !> Intended to be used for MP2 geometry optimizations.
  !> \author Kasper Kristensen
  !> \date November 2011
  subroutine get_mp2gradient_and_energy_from_inputs(MyLsitem,F,D,S,C,natoms,MP2gradient,EMP2,Eerr)

    implicit none
    !> LSitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> Fock matrix 
    type(matrix),intent(inout) :: F
    !> HF density matrix 
    type(matrix),intent(in) :: D
    !> Overlap matrix
    type(matrix),intent(inout) :: S
    !> MO coefficients 
    type(matrix),intent(inout) :: C
    !> Number of atoms in molecule
    integer,intent(in) :: natoms
    !> MP2 molecular gradient
    real(realk), intent(inout) :: mp2gradient(3,natoms)
    !> Total MP2 energy (Hartree-Fock + correlation contribution)
    real(realk),intent(inout) :: EMP2
    !> Difference between intrinsic energy error at this and the previous geometry
    !> (zero for single point calculation or first geometry step)
    real(realk),intent(inout) :: Eerr
    real(realk) :: Ecorr,EHF
    type(fullmolecule) :: MyMolecule
    integer :: nBasis,nOcc,nUnocc
    type(ccorbital), pointer :: OccOrbitals(:)
    type(ccorbital), pointer :: UnoccOrbitals(:)
    integer :: i
    real(realk), pointer :: DistanceTable(:,:)

    write(DECinfo%output,*) 'Calculating DEC-MP2 gradient, FOT = ', DECinfo%FOT

    ! Sanity checks
    ! *************
    if(.not. DECinfo%gradient) then
       call lsquit('get_mp2gradient_from_inputs called - but gradient keyword is not set!  &
            & Suggstion: Insert .gradient keyword in **DEC section',DECinfo%output)
    end if
    if(DECinfo%single_calculation) then
       call lsquit('get_mp2gradient_from_inputs do not work for single fragment calculations! &
       & Suggstion: Remove keywords .singleFragment and .singlePair in **DEC section',DECinfo%output)
    end if

    write(DECinfo%output,*) 'Calculating MP2 energy and gradient from Fock, density, overlap, and MO inputs...'


    ! Get informations about full molecule
    ! ************************************
    call molecule_init_from_inputs(MyMolecule,mylsitem,F,S,C)
    nOcc = MyMolecule%numocc
    nUnocc = MyMolecule%numvirt
    nBasis = MyMolecule%nbasis

    ! No reason to save F,S and C twice. Delete the ones in matrix format and reset at the end
    call mat_free(F)
    call mat_free(S)
    call mat_free(C)


    ! Calculate distance matrix 
    ! *************************
    call mem_alloc(DistanceTable,nAtoms,nAtoms)
    DistanceTable=0.0E0_realk
    call GetDistances(DistanceTable,nAtoms,mylsitem,DECinfo%output) ! distances in atomic units


    ! Analyze basis and create orbitals
    ! *********************************
    call mem_alloc(OccOrbitals,nOcc)
    call mem_alloc(UnoccOrbitals,nUnocc)
    call GenerateOrbitals_driver(MyMolecule,mylsitem,nocc,nunocc,natoms, &
         & OccOrbitals, UnoccOrbitals, DistanceTable)


    ! -- Calculate molecular MP2 gradient
    call main_fragment_driver(MyMolecule,mylsitem,D,&
         & OccOrbitals,UnoccOrbitals, &
         & natoms,nocc,nunocc,DistanceTable,EHF,Ecorr,mp2gradient,Eerr)

    ! Total MP2 energy: EHF + Ecorr
    EMP2 = EHF + Ecorr

    ! Restore input matrices
    call molecule_copyback_FSC_matrices(MyMolecule,F,S,C)

    ! Free molecule structure and other stuff
    call molecule_finalize(MyMolecule)
    call mem_dealloc(DistanceTable)
    
    ! Delete orbitals 
    do i=1,nOcc
       call orbital_free(OccOrbitals(i))
    end do

    do i=1,nUnocc
       call orbital_free(UnoccOrbitals(i))
    end do

    call mem_dealloc(OccOrbitals)
    call mem_dealloc(UnOccOrbitals)

    ! Set Eerr equal to the difference between the intrinsic error at this geometry
    ! (the current value of Eerr) and the intrinsic error at the previous geometry.
    call dec_get_error_difference(Eerr)

    ! Update number of DEC calculations for given FOT level
    DECinfo%ncalc(DECinfo%FOTlevel) = DECinfo%ncalc(DECinfo%FOTlevel) +1
    call print_calculation_bookkeeping()

  end subroutine get_mp2gradient_and_energy_from_inputs




  !> \brief Calculate MP2 energy using DEC scheme
  !> using the Fock, density, and MO matrices and input - rather than reading them from file
  !> as is done in a "conventional" DEC calculation which uses the dec_main_prog subroutine.
  !> \author Kasper Kristensen
  !> \date November 2011
  subroutine get_total_mp2energy_from_inputs(MyLsitem,F,D,S,C,EMP2,Eerr)

    implicit none
    !> LSitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> Fock matrix 
    type(matrix),intent(inout) :: F
    !> HF density matrix 
    type(matrix),intent(in) :: D
    !> Overlap matrix
    type(matrix),intent(inout) :: S
    !> MO coefficients 
    type(matrix),intent(inout) :: C
    !> Total MP2 energy (Hartree-Fock + correlation contribution)
    real(realk),intent(inout) :: EMP2
    !> Difference between intrinsic energy error at this and the previous geometry
    !> (zero for single point calculation or first geometry step)
    real(realk),intent(inout) :: Eerr
    real(realk) :: Ecorr,EHF
    type(fullmolecule) :: MyMolecule
    integer :: nBasis,nOcc,nUnocc,natoms
    type(ccorbital), pointer :: OccOrbitals(:)
    type(ccorbital), pointer :: UnoccOrbitals(:)
    real(realk), pointer :: DistanceTable(:,:), dummy(:,:)
    integer :: i
    logical :: save_first_order, save_grad, save_dens
    ! Quick solution to ensure that the MP2 gradient contributions are not set
    save_first_order = DECinfo%first_order
    save_grad = DECinfo%gradient
    save_dens = DECinfo%MP2density
    DECinfo%first_order = .false.
    DECinfo%gradient = .false.
    DECinfo%MP2density=.false.
    
    write(DECinfo%output,*) 'Calculating DEC-MP2 energy, FOT = ', DECinfo%FOT

    ! Get informations about full molecule
    ! ************************************
    call molecule_init_from_inputs(MyMolecule,mylsitem,F,S,C)

    ! No reason to save F,S and C twice. Delete the ones in matrix format and reset at the end
    call mat_free(F)
    call mat_free(S)
    call mat_free(C)

    nOcc = MyMolecule%numocc
    nUnocc = MyMolecule%numvirt
    nBasis = MyMolecule%nbasis
    natoms = MyMolecule%natoms


    ! Calculate distance matrix 
    ! *************************
    call mem_alloc(DistanceTable,nAtoms,nAtoms)
    DistanceTable=0.0E0_realk
    call GetDistances(DistanceTable,nAtoms,mylsitem,DECinfo%output) ! distances in atomic units


    ! Analyze basis and create orbitals
    ! *********************************
    call mem_alloc(OccOrbitals,nOcc)
    call mem_alloc(UnoccOrbitals,nUnocc)
    call GenerateOrbitals_driver(MyMolecule,mylsitem,nocc,nunocc,natoms, &
         & OccOrbitals, UnoccOrbitals, DistanceTable)

    ! -- Calculate correlation energy
    call mem_alloc(dummy,3,natoms)
    call main_fragment_driver(MyMolecule,mylsitem,D,&
         & OccOrbitals,UnoccOrbitals, &
         & natoms,nocc,nunocc,DistanceTable,EHF,Ecorr,dummy,Eerr)
    call mem_dealloc(dummy)

    ! Total MP2 energy: EHF + Ecorr
    EMP2 = EHF + Ecorr


    ! Restore input matrices
    call molecule_copyback_FSC_matrices(MyMolecule,F,S,C)

    ! Free molecule structure and other stuff
    call molecule_finalize(MyMolecule)
    call mem_dealloc(DistanceTable)

    ! Delete orbitals 
    do i=1,nOcc
       call orbital_free(OccOrbitals(i))
    end do

    do i=1,nUnocc
       call orbital_free(UnoccOrbitals(i))
    end do

    call mem_dealloc(OccOrbitals)
    call mem_dealloc(UnOccOrbitals)


    ! Reset DEC parameters to the same as they were at input
    DECinfo%first_order = save_first_order
    DECinfo%gradient = save_grad
    DECinfo%MP2density = save_dens

    ! Set Eerr equal to the difference between the intrinsic error at this geometry
    ! (the current value of Eerr) and the intrinsic error at the previous geometry.
    call dec_get_error_difference(Eerr)

    ! Update number of DEC calculations for given FOT level
    DECinfo%ncalc(DECinfo%FOTlevel) = DECinfo%ncalc(DECinfo%FOTlevel) +1
    call print_calculation_bookkeeping()

  end subroutine get_total_mp2energy_from_inputs

  !> \brief Print number of DEC calculations for each FOT level (only interesting for geometry opt)
  !> \author Kasper Kristensen
  !> \date December 2012
  subroutine print_calculation_bookkeeping()
    implicit none
    integer :: i

    write(DECinfo%output,*) 
    write(DECinfo%output,*) 'DEC CALCULATION BOOK KEEPING'
    write(DECinfo%output,*) '============================'
    write(DECinfo%output,*) 
    do i=1,7
       write(DECinfo%output,'(1X,a,i6,a,i6)') '# calculations done for FOTlevel ',i, ' is ', DECinfo%ncalc(i)
    end do
    write(DECinfo%output,*) 

  end subroutine print_calculation_bookkeeping

end module dec_main_mod
