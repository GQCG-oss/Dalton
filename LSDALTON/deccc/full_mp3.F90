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
  use ccintegrals

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
    call third_order_energy_correction(MyMolecule,MyLsitem,dE3)

    ! Total third-order energy
    E3 = E2 + dE3

  end subroutine full_canonical_mp3



  !> Calculate third-order energy correction (Eq. 14.4.61 in purple book)
  !> Noddy implementation!
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
    type(array4) :: govov,gvvoo,gvvvv,goooo
    integer :: i,j,k,l,a,b,c,d
    real(realk) :: tbar,ttilde,eps,X
    real(realk),pointer :: t(:,:,:,:)
    integer :: nocc,nvirt,nbasis

    ! Dimensions
    nocc = MyMolecule%nocc
    nvirt = MyMolecule%nvirt
    nbasis = MyMolecule%nbasis


    ! Get necessary integrals
    call third_order_energy_correction_integrals(MyMolecule,mylsitem,govov,gvvoo,gvvvv,goooo)


    ! *********************************
    ! Canonical first order amplitudes
    ! ********************************
    call mem_alloc(t,nvirt,nocc,nvirt,nocc)
    do j=1,nocc
       do b=1,nvirt
          do i=1,nocc
             do a=1,nvirt

                ! Orbital energy difference
                eps = MyMolecule%vvfock%elm2(a,a) + MyMolecule%vvfock%elm2(b,b) &
                     & - MyMolecule%oofock%elm2(i,i) -  MyMolecule%oofock%elm2(j,j)

                ! Amplitude t (Eq. 14.4.41)
                t(a,i,b,j) = - govov%val(i,a,j,b) / eps

             end do
          end do
       end do
    end do



    ! ***************************************
    ! Energy using Eq. 14.4.61 in purple book
    ! ***************************************
    dE3=0.0_realk
    do i=1,nocc
       do j=1,nocc
          do a=1,nvirt
             do b=1,nvirt


                ! Construct tTILDE_{ij}^{ab} 
                ! **************************

                ! Multiplier (Eq. 14.4.42 combined with 14.4.41 + symmetry of integrals)
                tbar = 4.0_realk*t(a,i,b,j) - 2.0_realk*t(b,i,a,j)

                ! Multiplier using unrestricted summation (Eq. 14.4.59)
                if(a==b .and. i==j) then
                   ttilde = 2.0_realk*tbar
                else
                   ttilde = tbar
                end if

                ! Construct X_{ij}^{ab}
                ! *********************
                X=0.0_realk

                ! 1: 1/2 sum_{cd} t_ij^cd g_acbd
                do c=1,nvirt
                   do d=1,nvirt
                   end do
                end do


             end do
          end do
       end do
    end do

    call mem_dealloc(t)
    call array4_free(govov)
    call array4_free(gvvoo)
    call array4_free(gvvvv)
    call array4_free(goooo)

  end subroutine third_order_energy_correction


  !> Calculate integrals needed for third-order energy correction 
  !> Noddy implementation!
  !> \author Kasper Kristensen
  !> \date August 2015
  subroutine third_order_energy_correction_integrals(MyMolecule,mylsitem,govov,gvvoo,gvvvv,goooo)
    implicit none
    !> Full molecule info
    type(fullmolecule), intent(inout) :: MyMolecule
    !> Lsitem structure
    type(lsitem), intent(inout) :: mylsitem
    !> 2-electron MO integrals with occupied and virtual indices as shown
    type(array4),intent(inout) :: govov,gvvoo,gvvvv,goooo
    type(array4) :: gao
    type(array2) :: Co,Cv
    integer :: nocc,nvirt,nbasis

    ! Dimensions
    nocc = MyMolecule%nocc
    nvirt = MyMolecule%nvirt
    nbasis = MyMolecule%nbasis
    call get_full_eri(mylsitem,nbasis,gao)

    Co = array2_init([nbasis,nocc],MyMolecule%Co%elm2)
    Cv = array2_init([nbasis,nvirt],MyMolecule%Cv%elm2)

    govov = get_gmo_simple(gao,Co,Cv,Co,Cv)
    gvvoo = get_gmo_simple(gao,Cv,Cv,Co,Co)
    gvvvv = get_gmo_simple(gao,Cv,Cv,Cv,Cv)
    goooo = get_gmo_simple(gao,Co,Co,Co,Co)


    call array2_free(Co)
    call array2_free(Cv)
    call array4_free(gao)

  end subroutine third_order_energy_correction_integrals


end module full_mp3_module

