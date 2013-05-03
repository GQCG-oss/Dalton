
!> @file
!> association routines for pointers

module ptr_assoc_module
 use precision
 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  REALK ASSOCIATIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !with integers in dims on comilation size
  subroutine ass_D1to1(elms,ptr,dims)
    implicit none
    integer,intent(in) :: dims(1)
    real(realk), target,intent(in) :: elms(dims(1))
    real(realk), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_D1to1
  subroutine ass_D1to2(elms,ptr,dims)
    implicit none
    integer,intent(in) :: dims(2)
    real(realk), target,intent(in) :: elms(dims(1),dims(2))
    real(realk), pointer:: ptr(:,:)
    ptr => elms
  end subroutine ass_D1to2
  subroutine ass_D1to3(elms,ptr,dims)
    implicit none
    integer,intent(in) :: dims(3)
    real(realk), target,intent(in) :: elms(dims(1),dims(2),dims(3))
    real(realk), pointer:: ptr(:,:,:)
    ptr => elms
  end subroutine ass_D1to3
  subroutine ass_D1to4(elms,ptr,dims)
    implicit none
    integer,intent(in) :: dims(4)
    real(realk), target,intent(in) :: elms(dims(1),dims(2),dims(3),dims(4))
    real(realk), pointer:: ptr(:,:,:,:)
    ptr => elms
  end subroutine ass_D1to4
  subroutine ass_D1to5(elms,ptr,dims)
    implicit none
    integer,intent(in) :: dims(5)
    real(realk), target,intent(in) :: elms(dims(1),dims(2),dims(3),dims(4),dims(5))
    real(realk), pointer:: ptr(:,:,:,:,:)
    ptr => elms
  end subroutine ass_D1to5
  subroutine ass_D1to6(elms,ptr,dims)
    implicit none
    integer,intent(in) :: dims(6)
    real(realk), target,intent(in) :: elms(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6))
    real(realk), pointer:: ptr(:,:,:,:,:,:)
    ptr => elms
  end subroutine ass_D1to6
  subroutine ass_D1to7(elms,ptr,dims)
    implicit none
    integer,intent(in) :: dims(7)
    real(realk), target,intent(in) :: elms(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7))
    real(realk), pointer:: ptr(:,:,:,:,:,:,:)
    ptr => elms
  end subroutine ass_D1to7
  subroutine ass_D2to1(elms,ptr,dims)
    implicit none
    integer,intent(in) :: dims(2)
    real(realk), target,intent(in) :: elms(dims(1)*dims(2))
    real(realk), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_D2to1
  subroutine ass_D3to1(elms,ptr,dims)
    implicit none
    integer,intent(in) :: dims(3)
    real(realk), target,intent(in) :: elms(dims(1)*dims(2)*dims(3))
    real(realk), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_D3to1
  subroutine ass_D4to1(elms,ptr,dims)
    implicit none
    integer,intent(in) :: dims(4)
    real(realk), target,intent(in) :: elms(dims(1)*dims(2)*dims(3)*dims(4))
    real(realk), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_D4to1
  subroutine ass_D5to1(elms,ptr,dims)
    implicit none
    integer,intent(in) :: dims(5)
    real(realk), target,intent(in) :: elms(dims(1)*dims(2)*dims(3)*dims(4)*dims(5))
    real(realk), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_D5to1
  subroutine ass_D6to1(elms,ptr,dims)
    implicit none
    integer,intent(in) :: dims(6)
    real(realk), target,intent(in) :: elms(dims(1)*dims(2)*dims(3)*dims(4)*dims(5)*dims(6))
    real(realk), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_D6to1
  subroutine ass_D7to1(elms,ptr,dims)
    implicit none
    integer,intent(in) :: dims(7)
    real(realk), target,intent(in) :: elms(dims(1)*dims(2)*dims(3)*dims(4)*dims(5)*dims(6)*dims(7))
    real(realk), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_D7to1

  !WITH 32bit dims
  subroutine ass_4D1to2(elms,ptr,dims)
    implicit none
    integer(kind=4),intent(in) :: dims(2)
    real(realk), target,intent(in) :: elms(dims(1),dims(2))
    real(realk), pointer:: ptr(:,:)
    ptr => elms
  end subroutine ass_4D1to2
  subroutine ass_4D1to3(elms,ptr,dims)
    implicit none
    integer(kind=4),intent(in) :: dims(3)
    real(realk), target,intent(in) :: elms(dims(1),dims(2),dims(3))
    real(realk), pointer:: ptr(:,:,:)
    ptr => elms
  end subroutine ass_4D1to3
  subroutine ass_4D1to4(elms,ptr,dims)
    implicit none
    integer(kind=4),intent(in) :: dims(4)
    real(realk), target,intent(in) :: elms(dims(1),dims(2),dims(3),dims(4))
    real(realk), pointer:: ptr(:,:,:,:)
    ptr => elms
  end subroutine ass_4D1to4
  subroutine ass_4D1to5(elms,ptr,dims)
    implicit none
    integer(kind=4),intent(in) :: dims(5)
    real(realk), target,intent(in) :: elms(dims(1),dims(2),dims(3),dims(4),dims(5))
    real(realk), pointer:: ptr(:,:,:,:,:)
    ptr => elms
  end subroutine ass_4D1to5
  subroutine ass_4D1to6(elms,ptr,dims)
    implicit none
    integer(kind=4),intent(in) :: dims(6)
    real(realk), target,intent(in) :: elms(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6))
    real(realk), pointer:: ptr(:,:,:,:,:,:)
    ptr => elms
  end subroutine ass_4D1to6
  subroutine ass_4D1to7(elms,ptr,dims)
    implicit none
    integer(kind=4),intent(in) :: dims(7)
    real(realk), target,intent(in) :: elms(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7))
    real(realk), pointer:: ptr(:,:,:,:,:,:,:)
    ptr => elms
  end subroutine ass_4D1to7
  subroutine ass_4D2to1(elms,ptr,dims)
    implicit none
    integer(kind=4),intent(in) :: dims(2)
    real(realk), target,intent(in) :: elms(dims(1)*dims(2))
    real(realk), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_4D2to1
  subroutine ass_4D3to1(elms,ptr,dims)
    implicit none
    integer(kind=4),intent(in) :: dims(3)
    real(realk), target,intent(in) :: elms(dims(1)*dims(2)*dims(3))
    real(realk), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_4D3to1
  subroutine ass_4D4to1(elms,ptr,dims)
    implicit none
    integer(kind=4),intent(in) :: dims(4)
    real(realk), target,intent(in) :: elms(dims(1)*dims(2)*dims(3)*dims(4))
    real(realk), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_4D4to1
  subroutine ass_4D5to1(elms,ptr,dims)
    implicit none
    integer(kind=4),intent(in) :: dims(5)
    real(realk), target,intent(in) :: elms(dims(1)*dims(2)*dims(3)*dims(4)*dims(5))
    real(realk), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_4D5to1
  subroutine ass_4D6to1(elms,ptr,dims)
    implicit none
    integer(kind=4),intent(in) :: dims(6)
    real(realk), target,intent(in) :: elms(dims(1)*dims(2)*dims(3)*dims(4)*dims(5)*dims(6))
    real(realk), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_4D6to1
  subroutine ass_4D7to1(elms,ptr,dims)
    implicit none
    integer(kind=4),intent(in) :: dims(7)
    real(realk), target,intent(in) :: elms(dims(1)*dims(2)*dims(3)*dims(4)*dims(5)*dims(6)*dims(7))
    real(realk), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_4D7to1

  !WITH 64bit dims
  subroutine ass_8D1to2(elms,ptr,dims)
    implicit none
    integer(kind=8),intent(in) :: dims(2)
    real(realk), target,intent(in) :: elms(dims(1),dims(2))
    real(realk), pointer:: ptr(:,:)
    ptr => elms
  end subroutine ass_8D1to2
  subroutine ass_8D1to3(elms,ptr,dims)
    implicit none
    integer(kind=8),intent(in) :: dims(3)
    real(realk), target,intent(in) :: elms(dims(1),dims(2),dims(3))
    real(realk), pointer:: ptr(:,:,:)
    ptr => elms
  end subroutine ass_8D1to3
  subroutine ass_8D1to4(elms,ptr,dims)
    implicit none
    integer(kind=8),intent(in) :: dims(4)
    real(realk), target,intent(in) :: elms(dims(1),dims(2),dims(3),dims(4))
    real(realk), pointer:: ptr(:,:,:,:)
    ptr => elms
  end subroutine ass_8D1to4
  subroutine ass_8D1to5(elms,ptr,dims)
    implicit none
    integer(kind=8),intent(in) :: dims(5)
    real(realk), target,intent(in) :: elms(dims(1),dims(2),dims(3),dims(4),dims(5))
    real(realk), pointer:: ptr(:,:,:,:,:)
    ptr => elms
  end subroutine ass_8D1to5
  subroutine ass_8D1to6(elms,ptr,dims)
    implicit none
    integer(kind=8),intent(in) :: dims(6)
    real(realk), target,intent(in) :: elms(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6))
    real(realk), pointer:: ptr(:,:,:,:,:,:)
    ptr => elms
  end subroutine ass_8D1to6
  subroutine ass_8D1to7(elms,ptr,dims)
    implicit none
    integer(kind=8),intent(in) :: dims(7)
    real(realk), target,intent(in) :: elms(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7))
    real(realk), pointer:: ptr(:,:,:,:,:,:,:)
    ptr => elms
  end subroutine ass_8D1to7
  subroutine ass_8D2to1(elms,ptr,dims)
    implicit none
    integer(kind=8),intent(in) :: dims(2)
    real(realk), target,intent(in) :: elms(dims(1)*dims(2))
    real(realk), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_8D2to1
  subroutine ass_8D3to1(elms,ptr,dims)
    implicit none
    integer(kind=8),intent(in) :: dims(3)
    real(realk), target,intent(in) :: elms(dims(1)*dims(2)*dims(3))
    real(realk), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_8D3to1
  subroutine ass_8D4to1(elms,ptr,dims)
    implicit none
    integer(kind=8),intent(in) :: dims(4)
    real(realk), target,intent(in) :: elms(dims(1)*dims(2)*dims(3)*dims(4))
    real(realk), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_8D4to1
  subroutine ass_8D5to1(elms,ptr,dims)
    implicit none
    integer(kind=8),intent(in) :: dims(5)
    real(realk), target,intent(in) :: elms(dims(1)*dims(2)*dims(3)*dims(4)*dims(5))
    real(realk), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_8D5to1
  subroutine ass_8D6to1(elms,ptr,dims)
    implicit none
    integer(kind=8),intent(in) :: dims(6)
    real(realk), target,intent(in) :: elms(dims(1)*dims(2)*dims(3)*dims(4)*dims(5)*dims(6))
    real(realk), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_8D6to1
  subroutine ass_8D7to1(elms,ptr,dims)
    implicit none
    integer(kind=8),intent(in) :: dims(7)
    real(realk), target,intent(in) :: elms(dims(1)*dims(2)*dims(3)*dims(4)*dims(5)*dims(6)*dims(7))
    real(realk), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_8D7to1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  INTEGER ASSOCIATIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !BOTH 32bit INTEGERS
  subroutine ass_44I2to1(elms,ptr,dims)
    implicit none
    integer(kind=4),intent(in) :: dims(2)
    integer(kind=4), target,intent(in) :: elms(dims(1)*dims(2))
    integer(kind=4), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_44I2to1
  subroutine ass_44I3to1(elms,ptr,dims)
    implicit none
    integer(kind=4),intent(in) :: dims(3)
    integer(kind=4), target,intent(in) :: elms(dims(1)*dims(2)*dims(3))
    integer(kind=4), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_44I3to1
  subroutine ass_44I4to1(elms,ptr,dims)
    implicit none
    integer(kind=4),intent(in) :: dims(4)
    integer(kind=4), target,intent(in) :: elms(dims(1)*dims(2)*dims(3)*dims(4))
    integer(kind=4), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_44I4to1


  !DIMENSION ON 64bit, ELEMENTS ON 32bit
  subroutine ass_48I2to1(elms,ptr,dims)
    implicit none
    integer(kind=8),intent(in) :: dims(2)
    integer(kind=4), target,intent(in) :: elms(dims(1)*dims(2))
    integer(kind=4), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_48I2to1
  subroutine ass_48I3to1(elms,ptr,dims)
    implicit none
    integer(kind=8),intent(in) :: dims(3)
    integer(kind=4), target,intent(in) :: elms(dims(1)*dims(2)*dims(3))
    integer(kind=4), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_48I3to1
  subroutine ass_48I4to1(elms,ptr,dims)
    implicit none
    integer(kind=8),intent(in) :: dims(4)
    integer(kind=4), target,intent(in) :: elms(dims(1)*dims(2)*dims(3)*dims(4))
    integer(kind=4), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_48I4to1

  !DIMENSION ON 32bit, ELEMENTS ON 64bit
  subroutine ass_84I2to1(elms,ptr,dims)
    implicit none
    integer(kind=4),intent(in) :: dims(2)
    integer(kind=8), target,intent(in) :: elms(dims(1)*dims(2))
    integer(kind=8), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_84I2to1
  subroutine ass_84I3to1(elms,ptr,dims)
    implicit none
    integer(kind=4),intent(in) :: dims(3)
    integer(kind=8), target,intent(in) :: elms(dims(1)*dims(2)*dims(3))
    integer(kind=8), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_84I3to1
  subroutine ass_84I4to1(elms,ptr,dims)
    implicit none
    integer(kind=4),intent(in) :: dims(4)
    integer(kind=8), target,intent(in) :: elms(dims(1)*dims(2)*dims(3)*dims(4))
    integer(kind=8), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_84I4to1

  !DIMENSION ON 64bit, ELEMENTS ON 64bit
  subroutine ass_88I2to1(elms,ptr,dims)
    implicit none
    integer(kind=8),intent(in) :: dims(2)
    integer(kind=8), target,intent(in) :: elms(dims(1)*dims(2))
    integer(kind=8), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_88I2to1
  subroutine ass_88I3to1(elms,ptr,dims)
    implicit none
    integer(kind=8),intent(in) :: dims(3)
    integer(kind=8), target,intent(in) :: elms(dims(1)*dims(2)*dims(3))
    integer(kind=8), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_88I3to1
  subroutine ass_88I4to1(elms,ptr,dims)
    implicit none
    integer(kind=8),intent(in) :: dims(4)
    integer(kind=8), target,intent(in) :: elms(dims(1)*dims(2)*dims(3)*dims(4))
    integer(kind=8), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_88I4to1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  LOGICAL ASSOCIATIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ass_44L2to1(elms,ptr,dims)
    implicit none
    integer(kind=4),intent(in) :: dims(2)
    logical(kind=4), target,intent(in) :: elms(dims(1)*dims(2))
    logical(kind=4), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_44L2to1

  subroutine ass_48L2to1(elms,ptr,dims)
    implicit none
    integer(kind=8),intent(in) :: dims(2)
    logical(kind=4), target,intent(in) :: elms(dims(1)*dims(2))
    logical(kind=4), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_48L2to1

  subroutine ass_84L2to1(elms,ptr,dims)
    implicit none
    integer(kind=4),intent(in) :: dims(2)
    logical(kind=8), target,intent(in) :: elms(dims(1)*dims(2))
    logical(kind=8), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_84L2to1

  subroutine ass_88L2to1(elms,ptr,dims)
    implicit none
    integer(kind=8),intent(in) :: dims(2)
    logical(kind=8), target,intent(in) :: elms(dims(1)*dims(2))
    logical(kind=8), pointer:: ptr(:)
    ptr => elms
  end subroutine ass_88L2to1
end module ptr_assoc_module
