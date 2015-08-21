!> @file
!> Full calculation
!> This file is mainly a playground for new developments, not intended to be included in a release.

module full_mp3_module

  use fundamental
  use precision
  use typedeftype!,only:lsitem
  use typedef
  use dec_typedef_module
  use matrix_module
  use matrix_operations
  use memory_handling
  use MemoryLeakToolMod

  !  DEC DEPENDENCIES (within deccc directory)   
  !  *****************************************
  use dec_tools_module
  use dec_fragment_utils
  use array4_simple_operations
  use array3_simple_operations
  use array2_simple_operations
  use fullmp2

  public :: full_canonical_mp3
  private

contains

  !> Calculate MP3 energy using canonical orbital.
  !> MP3 energy is the sum of the MP2 energy (Eq. 14.4.56 in purple book)
  !> and the third-order correction (Eq. 14.4.61 in purple book).
  !> \author Kasper Kristensen
  !> \date August 2015
  subroutine full_canonical_mp3(MyMolecule,MyLsitem,E3)
    implicit none
    !> Full molecule info
    type(fullmolecule), intent(inout) :: MyMolecule
    !> Lsitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> Canonical MP3 correlation energy (i.e. sum of second and third order energies)
    real(realk),intent(inout) :: E3
    real(realk) :: E2, dE3

    ! Second-order energy (MP2 energy)
    call full_canonical_mp2(MyMolecule,MyLsitem,E2)

    ! Third order-correction to energy
    dE3 = 0.0_realk

    ! Total third-order energy
    E3 = E2 + dE3

  end subroutine full_canonical_mp3



  !> Calculate third-order energy correction ((Eq. 14.4.61 in purple book)
  !> \author Kasper Kristensen
  !> \date August 2015  
  subroutine third_order_energy_correction(MyMolecule,MyLsitem,dE3)
    implicit none
    !> Full molecule info
    type(fullmolecule), intent(inout) :: MyMolecule
    !> Lsitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> Third-order energy correction
    real(realk),intent(inout) :: dE3

    dE3=0.0_realk

  end subroutine third_order_energy_correction

end module full_mp3_module

