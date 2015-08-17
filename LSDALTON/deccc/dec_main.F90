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
  use reorder_frontend_module 
  use tensor_tester_module
  use Matrix_util!, only: get_AO_gradient
  use configurationType

  !> F12routines
  !use f12_routines_module
  
  ! DEC DEPENDENCIES (within deccc directory) 
  ! *****************************************
  use snoop_main_module
  use dec_fragment_utils
  use array3_memory_manager!,only: print_memory_currents_3d
  use array4_memory_manager!,only: print_memory_currents4
  use full_molecule!,only: molecule_init, molecule_finalize,molecule_init_from_inputs
  use orbital_operations!,only: check_lcm_against_canonical
  use full_molecule
  use mp2_gradient_module
  use f12_routines_module
  use dec_driver_module,only: dec_wrapper
  use full,only: full_driver

public :: dec_main_prog_input, dec_main_prog_file, &
     & get_mp2gradient_and_energy_from_inputs, get_total_CCenergy_from_inputs
private

contains

  !> Wrapper for main DEC program to use when Fock,density, and MO coefficient
  !> matrices are available from HF calculation.
  !> \author Kasper Kristensen
  !> \date April 2013
  subroutine dec_main_prog_input(mylsitem,config,F,D,C,E)
    implicit none

    !> Integral info
    type(lsitem), intent(inout) :: mylsitem
    !> Config info
    type(configItem),intent(in)  :: config
    !> Fock matrix 
    type(matrix),intent(inout) :: F
    !> HF density matrix 
    type(matrix),intent(in) :: D
    !> MO coefficients 
    type(matrix),intent(inout) :: C
    !> CC Energy (intent out) 
    real(realk),intent(inout) :: E
    !local variables
    type(fullmolecule) :: Molecule
    !SCF gradient norm
    real(realk) :: gradnorm

    print *, 'Hartree-Fock info comes directly from HF calculation...'


    ! Get informations about full molecule
    ! ************************************
    call molecule_init_from_inputs(Molecule,mylsitem,F,C,D)
    
    !> F12
    !call molecule_init_f12(molecule,mylsitem,D)

    ! Fock, overlap, and MO coefficient matrices are now stored
    ! in Molecule, and there is no reason to store them twice.
    ! So we delete them now and reset them at the end.
    call mat_free(F)
    call mat_free(C)
    
    call dec_main_prog(MyLsitem,config,molecule,D,E)

    ! Restore input matrices
    call molecule_copyback_FC_matrices(mylsitem,Molecule,F,C)

    ! Delete molecule structure
    call molecule_finalize(molecule,.true.)


  end subroutine dec_main_prog_input

  subroutine SCFgradError(gradnorm)
    implicit none
    real(realk) :: gradnorm
    WRITE(DECinfo%output,'(A,ES16.8)')'WARNING: The SCF gradient norm = ',gradnorm
    WRITE(DECinfo%output,'(A,ES16.8)')'WARNING: Is greater then the FOT=',DECinfo%FOT
    WRITE(DECinfo%output,'(A,ES16.8)')'WARNING: The Error of the DEC calculation would be determined by the SCF error.'
    WRITE(DECinfo%output,'(A,ES16.8)')'WARNING: Tighten the SCF threshold!'
!    call lsquit('DEC ERROR: SCF gradient too large',-1)
  end subroutine SCFgradError

  !> Wrapper for main DEC program to use when Fock,density,overlap, and MO coefficient
  !> matrices are not directly available from current HF calculation, but rather needs to
  !> be read in from previous HF calculation.
  !> \author Kasper Kristensen
  !> \date April 2013
  subroutine dec_main_prog_file(mylsitem,config)

    implicit none
    !> Integral info
    type(lsitem), intent(inout) :: mylsitem
    !> Config info
    type(configItem),intent(in)  :: config
    type(matrix) :: D
    type(fullmolecule) :: Molecule
    integer :: nbasis
    real(realk) :: E
    E = 0.0E0_realk
    
    ! Minor tests
    ! ***********
    !Array test
    if (DECinfo%tensor_test)then
      print *,"TEST ARRAY MODULE"
      call test_tensor_struct(DECinfo%output)
      return
    endif
    ! Reorder test
    if (DECinfo%reorder_test)then
      print *,"TEST REORDERINGS"
      call test_array_reorderings(DECinfo%output)
      return
    endif

    print *, 'Hartree-Fock info is read from file...'

    ! Get density matrix
    Molecule%nbasis = get_num_basis_functions(mylsitem)
    call dec_get_density_matrix_from_file(Molecule%nbasis,D)
    
    ! Get informations about full molecule by reading from file
    call molecule_init_from_files(molecule,mylsitem,D)
 
    !> F12
    !call molecule_init_f12(molecule,mylsitem,D)
       
    ! Main DEC program
    call dec_main_prog(MyLsitem,config,molecule,D,E)
     
    ! Delete molecule structure and density
    call molecule_finalize(molecule,.true.)
    call mat_free(D)

  end subroutine dec_main_prog_file



  !> \brief Main DEC program.
  !> \author Marcin Ziolkowski (modified for Dalton by Kasper Kristensen)
  !> \date September 2010
  subroutine dec_main_prog(MyLsitem,config,molecule,D,E)

    implicit none
    !> Integral info
    type(lsitem), intent(inout) :: mylsitem
    !> Config info
    type(configItem),intent(in)  :: config
    !> Molecule info
    type(fullmolecule),intent(inout) :: molecule
    !> HF density matrix
    type(matrix),intent(in) :: D
    !> Energy (maybe HF energy as input, CC energy as output) 
    real(realk),intent(inout) :: E
    character(len=10) :: program_version
    character(len=50) :: MyHostname
    integer, dimension(8) :: values
    real(realk) :: tcpu1, twall1, tcpu2, twall2, EHF,Ecorr,Eerr,ES2
    real(realk) :: molgrad(3,Molecule%natoms)  
    
    ! Set DEC memory
    call get_memory_for_dec_calculation()

    ! Perform SNOOP calculation and skip DEC calculation
    ! (at some point SNOOP and DEC might be merged)
    if(DECinfo%SNOOP) then
       write(DECinfo%output,*) '***********************************************************'
       write(DECinfo%output,*) '      Performing SNOOP interaction energy calculation...'
       write(DECinfo%output,*) '***********************************************************'
       call snoop_driver(mylsitem,config,Molecule,D)
       return
    end if

    ! Sanity check: LCM orbitals span the same space as canonical orbitals 
    if(DECinfo%check_lcm_orbitals) then
       call check_lcm_against_canonical(molecule,MyLsitem)
       return
    end if

    if(DECinfo%force_Occ_SubSystemLocality)then
       call force_Occupied_SubSystemLocality(molecule,MyLsitem)
    endif

    if(DECinfo%check_Occ_SubSystemLocality)then
       call check_Occupied_SubSystemLocality(molecule,MyLsitem)
    endif


    ! Actual DEC calculation
    ! **********************

    call LSTIMER('START',tcpu1,twall1,DECinfo%output)

    if(.not. DECinfo%full_molecular_cc) then
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
       write(DECinfo%output,*)
    end if

    if(DECinfo%use_canonical) then
       write(DECinfo%output,*) 'Using canonical orbitals as requested in input!'
       if(.not. DECinfo%full_molecular_cc) then
          write(DECinfo%output,*) 'Warning: This may cause meaningless results for the DEC scheme'
       end if
       write(DECinfo%output,*)
       write(DECinfo%output,*)
    end if
    
    
    if(DECinfo%F12 .and. DECinfo%F12singles ) then
       call F12singles_driver(Molecule,MyLsitem,D)
    endif

    
    if(DECinfo%full_molecular_cc) then
       ! -- Call full molecular CC
       write(DECinfo%output,'(/,a,/)') 'Full molecular calculation is carried out...'
       call full_driver(molecule,mylsitem,D,EHF,Ecorr)
       ! --
    else
       ! -- Initialize DEC driver for energy calculation
       write(DECinfo%output,'(/,a,/)') 'DEC calculation is carried out...'
       call DEC_wrapper(molecule,mylsitem,D,EHF,Ecorr,molgrad,Eerr)
       ! --
    end if
    E = EHF + Ecorr 

    ! Update number of DEC calculations for given FOT level
    DECinfo%ncalc(DECinfo%FOTlevel) = DECinfo%ncalc(DECinfo%FOTlevel) +1
    call print_calculation_bookkeeping()

    call LSTIMER('START',tcpu2,twall2,DECinfo%output)

    ! Print memory summary
    write(DECinfo%output,*)
    write(DECinfo%output,*) 'CC Memory summary'
    write(DECinfo%output,*) '-----------------'
    call print_memory_currents4(DECinfo%output)
    write(DECinfo%output,*) '------------------'
    call print_memory_currents_3d(DECinfo%output)
    write(DECinfo%output,*) '------------------'
    write(DECinfo%output,*)
    write(DECinfo%output,'(/,a)') '------------------------------------------------------'
    if(DECinfo%full_molecular_cc) then
       write(DECinfo%output,'(a,g20.6,a)') 'Total CPU  time used in CC           :',tcpu2-tcpu1,' s'
       write(DECinfo%output,'(a,g20.6,a)') 'Total Wall time used in CC           :',twall2-twall1,' s'
    else
       write(DECinfo%output,'(a,g20.6,a)') 'Total CPU  time used in DEC          :',tcpu2-tcpu1,' s'
       write(DECinfo%output,'(a,g20.6,a)') 'Total Wall time used in DEC          :',twall2-twall1,' s'
    end if
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
    if(DECinfo%full_molecular_cc) then
       write(DECinfo%output,'(a)')   '                          -- end of CC program --'
    else
       write(DECinfo%output,'(a)')   '                          -- end of DEC program --'
    end if
    write(DECinfo%output,'(a,/)') '============================================================================='
    write(DECinfo%output,*)
    write(DECinfo%output,*)

  end subroutine dec_main_prog


  !> \brief Calculate MP2 energy and molecular gradient using DEC scheme
  !> using the Fock, density, and MO matrices and input - rather than reading them from file
  !> as is done in a "conventional" DEC calculation which uses the dec_main_prog subroutine.
  !> Intended to be used for MP2 geometry optimizations.
  !> \author Kasper Kristensen
  !> \date November 2011
  subroutine get_mp2gradient_and_energy_from_inputs(MyLsitem,F,D,C,natoms,MP2gradient,EMP2,Eerr)

    implicit none
    !> LSitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> Fock matrix 
    type(matrix),intent(inout) :: F
    !> HF density matrix 
    type(matrix),intent(in) :: D
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
    type(fullmolecule) :: Molecule

    write(DECinfo%output,*) 'Calculating DEC-MP2 gradient, FOT = ', DECinfo%FOT

    ! Sanity check
    ! ************
    if( (.not. DECinfo%gradient) .or. (.not. DECinfo%first_order) &
         &  .or. (DECinfo%ccmodel/=MODEL_MP2 .and. DECinfo%ccmodel/=MODEL_RIMP2 ) ) then
       call lsquit("ERROR(get_mp2gradient_and_energy_from_inputs): Inconsitent input!!",DECinfo%output)
    end if
    write(DECinfo%output,*) 'Calculating MP2 energy and gradient from Fock, density, overlap, and MO inputs...'

    ! Set DEC memory
    call get_memory_for_dec_calculation()

    ! Get informations about full molecule
    ! ************************************
    call molecule_init_from_inputs(Molecule,mylsitem,F,C,D)

    ! No reason to save F and C twice. Delete the ones in matrix format and reset at the end
    call mat_free(F)
    call mat_free(C)

    ! -- Calculate molecular MP2 gradient and correlation energy
    call DEC_wrapper(Molecule,mylsitem,D,EHF,Ecorr,mp2gradient,Eerr)

    ! Total MP2 energy: EHF + Ecorr
    EMP2 = EHF + Ecorr

    ! Restore input matrices
    call molecule_copyback_FC_matrices(mylsitem,Molecule,F,C)

    ! Free molecule structure and other stuff
    call molecule_finalize(Molecule,.true.)
    

    ! Set Eerr equal to the difference between the intrinsic error at this geometry
    ! (the current value of Eerr) and the intrinsic error at the previous geometry.
    call dec_get_error_for_geoopt(Eerr)

    ! Update number of DEC calculations for given FOT level
    DECinfo%ncalc(DECinfo%FOTlevel) = DECinfo%ncalc(DECinfo%FOTlevel) +1
    call print_calculation_bookkeeping()

  end subroutine get_mp2gradient_and_energy_from_inputs




  !> \brief Calculate MP2 energy using DEC scheme
  !> using the Fock, density, and MO matrices and input - rather than reading them from file
  !> as is done in a "conventional" DEC calculation which uses the dec_main_prog subroutine.
  !> \author Kasper Kristensen
  !> \date November 2011
  subroutine get_total_CCenergy_from_inputs(MyLsitem,F,D,C,ECC,Eerr)

    implicit none
    !> LSitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> Fock matrix 
    type(matrix),intent(inout) :: F
    !> HF density matrix 
    type(matrix),intent(in) :: D
    !> MO coefficients 
    type(matrix),intent(inout) :: C
    !> Total CC energy (Hartree-Fock + correlation contribution)
    real(realk),intent(inout) :: ECC
    !> Difference between intrinsic energy error at this and the previous geometry
    !> (zero for single point calculation or first geometry step)
    real(realk),intent(inout) :: Eerr
    real(realk) :: Ecorr,EHF
    type(fullmolecule) :: Molecule
    real(realk), pointer :: dummy(:,:)
    integer :: i
    logical :: save_first_order, save_grad, save_dens
    ! Quick solution to ensure that the MP2 gradient contributions are not set
    save_first_order = DECinfo%first_order
    save_grad = DECinfo%gradient
    save_dens = DECinfo%density
    DECinfo%first_order = .false.
    DECinfo%gradient = .false.
    DECinfo%density=.false.
    
    write(DECinfo%output,*) 'Calculating DEC correlation energy, FOT = ', DECinfo%FOT

    ! Set DEC memory
    call get_memory_for_dec_calculation()

    ! Get informations about full molecule
    ! ************************************
    call molecule_init_from_inputs(Molecule,mylsitem,F,C,D)

    ! No reason to save F and C twice. Delete the ones in matrix format and reset at the end
    call mat_free(F)
    call mat_free(C)

    ! -- Calculate correlation energy
    if(DECinfo%full_molecular_cc) then
       ! -- Call full molecular CC
       write(DECinfo%output,'(/,a,/)') 'Full molecular calculation is carried out ...'
       call full_driver(molecule,mylsitem,D,EHF,Ecorr)
       ! --
    else
       ! -- Initialize DEC driver for energy calculation
       write(DECinfo%output,'(/,a,/)') 'DEC calculation is carried out...'
       call mem_alloc(dummy,3,Molecule%natoms)
       call DEC_wrapper(Molecule,mylsitem,D,EHF,Ecorr,dummy,Eerr)
       call mem_dealloc(dummy)
       ! --
    end if

    ! Total CC energy: EHF + Ecorr
    ECC = EHF + Ecorr
    ! Restore input matrices
    call molecule_copyback_FC_matrices(mylsitem,Molecule,F,C)

    ! Free molecule structure and other stuff
    call molecule_finalize(Molecule,.true.)

    ! Reset DEC parameters to the same as they were at input
    DECinfo%first_order = save_first_order
    DECinfo%gradient = save_grad
    DECinfo%density = save_dens

    ! Set Eerr equal to the difference between the intrinsic error at this geometry
    ! (the current value of Eerr) and the intrinsic error at the previous geometry.
    call dec_get_error_for_geoopt(Eerr)

    ! Update number of DEC calculations for given FOT level
    DECinfo%ncalc(DECinfo%FOTlevel) = DECinfo%ncalc(DECinfo%FOTlevel) +1
    call print_calculation_bookkeeping()

  end subroutine get_total_CCenergy_from_inputs

  !> \brief Print number of DEC calculations for each FOT level (only interesting for geometry opt)
  !> \author Kasper Kristensen
  !> \date December 2012
  subroutine print_calculation_bookkeeping()
    implicit none
    integer :: i

    ! This only be done for geometry optimization and DEC scheme
    if( (.not. DECinfo%full_molecular_cc) .and. DECinfo%gradient) then
       write(DECinfo%output,*) 
       write(DECinfo%output,*) 'DEC CALCULATION BOOK KEEPING'
       write(DECinfo%output,*) '============================'
       write(DECinfo%output,*) 
       do i=1,nFOTs
          if(DECinfo%ncalc(i)/=0) then
             write(DECinfo%output,'(1X,a,i6,a,i6)') '# calculations done for FOTlevel ',i, ' is ', DECinfo%ncalc(i)
          end if
       end do
       write(DECinfo%output,*) 
    end if

  end subroutine print_calculation_bookkeeping

end module dec_main_mod
