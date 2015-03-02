module mat3d_mod
use precision
TYPE MAT3D
real(realk), pointer :: elements(:,:,:)
INTEGER              :: dim1,dim2,dim3
END TYPE MAT3D

CONTAINS
!> \brief Writes a mat3d-item to file
!> \author S. Reine
!> \date 2011-01-14
!> \param mat The MAT3D to be written to disk  
!> \param iunit The unit number of the file
SUBROUTINE write_mat3d_to_disk(mat,iunit)
implicit none
type(mat3d),intent(in) :: mat
integer,intent(in)     :: iunit

write(iunit) mat%dim1,mat%dim2,mat%dim3
write(iunit) mat%elements

END SUBROUTINE write_mat3d_to_disk

!> \brief Reads a mat3d-item from file
!> \author S. Reine
!> \date 2011-01-14
!> \param mat The MAT3D to be written to disk  
!> \param iunit The unit number of the file
SUBROUTINE read_mat3d_from_disk(mat,iunit,lupri)
implicit none
type(mat3d),intent(inout) :: mat
integer,intent(in)        :: iunit,lupri
!
integer :: dim1,dim2,dim3

Read(iunit) dim1,dim2,dim3
IF (dim1.NE.mat%dim1 .OR. dim2.NE.mat%dim2 .OR. dim3.NE.mat%dim3) THEN
  write(lupri,'(1X,A)')     'Dimension mismatch in read_mat3d_from_disk'
  write(lupri,'(5X,A,3I6)') '- on file',dim1,dim2,dim3
  write(lupri,'(5X,A,3I6)') '- in mat ',mat%dim1,mat%dim2,mat%dim3
  call lsquit('Dimension mismatch in read_mat3d_from_disk',lupri)
ENDIF
Read(iunit) mat%elements

END SUBROUTINE read_mat3d_from_disk

!> \brief This subroutine allocates memory and zeros a 3 dimensional matrix type
!> \author P. Merlot, S. Reine
!> \date 2010-05-04
!> \param threeIndex the 3D matrix pointer
!> \param m first dimension of the matrix
!> \param n second dimension of the matrix
!> \param k third dimension of the matrix
SUBROUTINE init_MAT3D(threeIndex,m,n,k)
use memory_handling
  implicit none
  type(MAT3D),intent(inout)	:: threeIndex
  integer,intent(in)		:: m,n,k
  threeIndex%dim1 = m
  threeIndex%dim2 = n
  threeIndex%dim3 = k
  call mem_alloc(threeIndex%elements,m,n,k)
  threeIndex%elements = 0E0_realk
END SUBROUTINE init_MAT3D

!> \brief This subroutine deallocate memory for a 3 dimensional matrix type
!> \author P. Merlot, S. Reine
!> \date 2010-05-04
!> \param threeIndex the 3D matrix pointer
!> \param m first dimension of the matrix
!> \param n second dimension of the matrix
!> \param k third dimension of the matrix
SUBROUTINE free_MAT3D(threeIndex)
use memory_handling
  implicit none
  type(MAT3D),intent(inout)	:: threeIndex
  threeIndex%dim1 = 0
  threeIndex%dim2 = 0
  threeIndex%dim3 = 0
  call mem_dealloc(threeIndex%elements)
END SUBROUTINE free_MAT3D

end module mat3d_mod

module mat2d_mod
use precision
TYPE MAT2D
real(realk), pointer :: elements(:,:)
INTEGER              :: dim1,dim2
END TYPE MAT2D

contains
subroutine init_mat2d(mat,m,n)
  use memory_handling
  implicit none
  type(mat2d),intent(inout) :: mat
  integer,intent(in)        :: m,n
  mat%dim1=m
  mat%dim2=n
  call mem_alloc(mat%elements,m,n)
  mat%elements=0E0_realk
end subroutine

subroutine free_mat2d(mat)
  use memory_handling
  implicit none
  type(mat2d),intent(inout) :: mat
  mat%dim1=0
  mat%dim2=0
  call mem_dealloc(mat%elements)
end subroutine free_mat2d

end module mat2d_mod
