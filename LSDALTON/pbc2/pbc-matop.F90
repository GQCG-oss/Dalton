#ifdef MOD_UNRELEASED
MODULE pbc_matrix_operations
USE precision
!use matrix_operations
use typedef
use lattice_type
use lattice_vectors
use files
use fundamental

CONTAINS

!This routine writes a matrix to screen in a consistent form.
SUBROUTINE write_matrix(A,nrows,mcol,lupri)!,mattxt)
IMPLICIT NONE
INTEGER, INTENT(IN) :: nrows, mcol
INTEGER,INTENT(IN),OPTIONAL :: lupri
!CHARACTER(len=40),INTENT(IN),OPTIONAL :: mattxt
REAL(realk), INTENT(IN) :: A(nrows,mcol)
!LOCAL VARIABLES
INTEGER :: i,j,k
CHARACTER(LEN=3) :: nline

 if(PRESENT(lupri)) then
   k=0
   DO i=1,nrows
    DO j=1,mcol
     k=k+1
     nline='no'
     if(k .eq. mcol) THEN
       k=0
       nline='yes'
     ENDIF
     write(lupri,'(E25.12)',advance=nline) A(i,j)
    ENDDO
   ENDDO
 else
   k=0
   DO i=1,nrows
    DO j=1,mcol
     k=k+1
     nline='no'
     if(k .eq. mcol) THEN
       k=0
       nline='yes'
     ENDIF
     write(*,'(E25.12)',advance=nline)  A(i,j)
    ENDDO
   ENDDO
 endif

END SUBROUTINE write_matrix

SUBROUTINE write_tmatrix(A,nrows,mcol,lupri)
IMPLICIT NONE
INTEGER, INTENT(IN) :: nrows, mcol
INTEGER,INTENT(IN),OPTIONAL :: lupri
!CHARACTER(len=40),INTENT(IN),OPTIONAL :: mattxt
REAL(realk), INTENT(IN) :: A(nrows,mcol)
!LOCAL VARIABLES
INTEGER :: i,j,k
CHARACTER(LEN=3) :: nline

 if(PRESENT(lupri)) then
   k=0
   DO i=1,nrows
    DO j=1,mcol
     k=k+1
     nline='no'
     if(k .eq. mcol) THEN
       k=0
       nline='yes'
     ENDIF
     write(lupri,'(E25.12)',advance=nline) A(j,i)
    ENDDO
   ENDDO
 else
   k=0
   DO i=1,nrows
    DO j=1,mcol
     k=k+1
     nline='no'
     if(k .eq. mcol) THEN
       k=0
       nline='yes'
     ENDIF
     write(*,'(E25.12)',advance=nline)  A(j,i)
    ENDDO
   ENDDO
 endif

END SUBROUTINE write_tmatrix


!This routine writes a complex matrix to screen in a consistent form.
SUBROUTINE write_zmatrix(A,nrows,mcol,lupri)
IMPLICIT NONE
INTEGER, INTENT(IN) :: nrows, mcol
COMPLEX(complexk), INTENT(IN) :: A(nrows,mcol)
INTEGER, INTENT(IN),OPTIONAL :: lupri
!LOCAL VARIABLES
INTEGER :: i,j,k
CHARACTER(LEN=3) :: nline

 if(PRESENT(lupri)) then
   k=0
   DO i=1,nrows
    DO j=1,mcol
     k=k+1
     nline='no'
     if(k .eq. mcol) THEN
       k=0
       nline='yes'
     ENDIF
     write(lupri,101,advance=nline) A(i,j)
     101 FORMAT('(',2E25.12,')')
    ENDDO
   ENDDO
 else
   k=0
   DO i=1,nrows
    DO j=1,mcol
     k=k+1
     nline='no'
     if(k .eq. mcol) THEN
       k=0
       nline='yes'
     ENDIF
     write(*,100,advance=nline) A(i,j)
     100 FORMAT('(',2E25.12,')')
    ENDDO
   ENDDO
 endif


END SUBROUTINE write_zmatrix

SUBROUTINE write_tzmatrix(A,nrows,mcol,lupri)
IMPLICIT NONE
INTEGER, INTENT(IN) :: nrows, mcol
COMPLEX(complexk), INTENT(IN) :: A(nrows,mcol)
INTEGER, INTENT(IN),OPTIONAL :: lupri
!LOCAL VARIABLES
INTEGER :: i,j,k
CHARACTER(LEN=3) :: nline

 if(PRESENT(lupri)) then
   k=0
   DO i=1,nrows
    DO j=1,mcol
     k=k+1
     nline='no'
     if(k .eq. mcol) THEN
       k=0
       nline='yes'
     ENDIF
     write(lupri,101,advance=nline) A(j,i)
     101 FORMAT('(',2E18.9,')')
    ENDDO
   ENDDO
 else
   k=0
   DO i=1,nrows
    DO j=1,mcol
     k=k+1
     nline='no'
     if(k .eq. mcol) THEN
       k=0
       nline='yes'
     ENDIF
     write(*,100,advance=nline) A(j,i)
     100 FORMAT('(',2E18.9,')')
    ENDDO
   ENDDO
 endif


END SUBROUTINE write_tzmatrix


!SUBROUTINE to write the fock matrix to disk
SUBROUTINE pbc_fockmat_write(Aop,nrows,ncols,oper1,oper2,diis,lu)
IMPLICIT NONE 
TYPE(lvec_list_t),INTENT(IN) :: Aop
INTEGER,INTENT(IN) :: nrows,ncols,oper1,oper2,lu
CHARACTER(len=12)          :: diis
!LOCAL VARIABLES
INTEGER :: layer
INTEGER :: l1,l2,l3

do layer =1,size(Aop%lvec)
   if(Aop%lvec(layer)%f1_computed .or. Aop%lvec(layer)%g2_computed) then
     l1=int(Aop%lvec(layer)%lat_coord(1))
     l2=int(Aop%lvec(layer)%lat_coord(2))
     l3=int(Aop%lvec(layer)%lat_coord(3))

#ifdef DEBUGPBC
     write(lu,*) 'Fock matrix',l1
     call mat_print(Aop%lvec(layer)%oper(oper2),1,nrows,1,ncols,lu)
#endif
     !write(*,*) 'layer', layer,l1,l2,l3,Aop%lvec(layer)%f1_computed
     call pbc_get_file_and_write(Aop,nrows,ncols,layer,oper1,oper2,diis)
   endif
enddo

END SUBROUTINE pbc_fockmat_write

SUBROUTINE pbc_densitymat_write(dmat,Aop,nrows,ncols,oper1,diis)
IMPLICIT NONE 
TYPE(lvec_list_t),INTENT(IN) :: Aop
TYPE(matrix),INTENT(IN) :: dmat(size(Aop%lvec))
INTEGER,INTENT(IN) :: nrows,ncols,oper1
CHARACTER(len=12)          :: diis
!LOCAL VARIABLES
INTEGER :: layer,iunit
INTEGER :: l1,l2,l3
CHARACTER(len=10) :: stl1,stl2,stl3
CHARACTER(LEN=100) :: filename

iunit=-1

do layer =1,size(Aop%lvec)
   l1=int(Aop%lvec(layer)%lat_coord(1))
   l2=int(Aop%lvec(layer)%lat_coord(2))
   l3=int(Aop%lvec(layer)%lat_coord(3))
   write(stl1,'(I5)')  l1
   write(stl2,'(I5)')  l2
   write(stl3,'(I5)')  l3
   stl1=adjustl(stl1)
   stl2=adjustl(stl2)
   stl3=adjustl(stl3)

   if((abs(l1) .le. Aop%ndmat .and. abs(l2) .le. Aop%ndmat) .and. abs(l3) .le. Aop%ndmat)then
   !    Now we are doing DIIS scf iterations and we need to store 
   !    the matrix for some of the iterations
   filename = trim(diis)//trim(adjustl(Aop%opdat(oper1)%basename))//'_'//trim(stl1)//'_'//trim(stl2)//'_'//trim(stl3)

    filename=adjustl(filename)
    filename=trim(filename)

#ifdef DEBUGPBC
    write(6,*) 'dmat, for l1 = ' , l1
    call mat_print(dmat(layer),1,nrows,1,ncols,6)
#endif

    call LSOPEN(iunit,filename,'unknown','UNFORMATTED')

    !write(iunit,*) ncols
    !We write the matrix to a binary file
    call write_binary_matrix(dmat(layer),nrows,ncols,iunit)

    call LSCLOSE(iunit,'KEEP')


   endif
enddo

END SUBROUTINE pbc_densitymat_write

! This function prints a T or W matrix in a format that allows it
! to be cut-n-pasted into MATLAB. It is for debugging purposes only.
!
!
subroutine pbc_matlab_print(matrix,maxrows,maxcols,text,lupri)
  implicit none

  character*(*) :: text
  integer,intent(in) :: maxrows, maxcols,lupri
  real(realk) :: matrix(maxrows,maxcols)

  integer :: row, col, col_max

  write(LUPRI,*) text,' = ['
  do row = 1,maxrows
     do col = 1,maxcols,5
        col_max = col + 4
        if (col_max .lt. maxcols) then
           write(LUPRI,'(5G18.8,A)') matrix(row,col:col_max),' ...'
        else if (col .eq. maxcols) then
           write(LUPRI,'(G18.8,A)') matrix(row,col),' ;'
        else if (col+1 .eq. maxcols) then
           write(LUPRI,'(2G18.8,A)') matrix(row,col:maxcols),' ;'
        else if (col+2 .eq. maxcols) then
           write(LUPRI,'(3G18.8,A)') matrix(row,col:maxcols),' ;'
        else if (col+3 .eq. maxcols) then
           write(LUPRI,'(4G18.8,A)') matrix(row,col:maxcols),' ;'
        else if (col+4 .eq. maxcols) then
           write(LUPRI,'(5G18.8,A)') matrix(row,col:maxcols),' ;'
        else
           call lsquit('pbc_matlab_print: programmer error',lupri)
        end if
     end do
  end do
  write(LUPRI,*) '   ]; '
end subroutine pbc_matlab_print



!SUBROUTINE to get filename for the operator, open the file and write 
!to disk
SUBROUTINE pbc_get_file_and_write(Aop,nrows,ncols,cell,oper1,oper2,diis)
IMPLICIT NONE 
TYPE(lvec_list_t),INTENT(IN) :: Aop
INTEGER,INTENT(IN) :: nrows,ncols,cell,oper1,oper2
CHARACTER(len=12)  :: diis
!LOCAL VARIABLES
INTEGER :: l1,l2,l3
CHARACTER(len=10) :: stl1,stl2,stl3
CHARACTER(LEN=100) :: filename
INTEGER :: iunit

  l1=int(Aop%lvec(cell)%lat_coord(1))
  l2=int(Aop%lvec(cell)%lat_coord(2))
  l3=int(Aop%lvec(cell)%lat_coord(3))
  write(stl1,'(I5)')  l1
  write(stl2,'(I5)')  l2
  write(stl3,'(I5)')  l3
  stl1=adjustl(stl1)
  stl2=adjustl(stl2)
  stl3=adjustl(stl3)

  iunit=-1

  !    Now we are doing DIIS scf iterations and we need to store 
  !    the matrix for some of the iterations
  filename = trim(diis)//trim(adjustl(Aop%opdat(oper1)%basename))//'_'//trim(stl1)//'_'//trim(stl2)//'_'//trim(stl3)
  
  filename=adjustl(filename)
  filename=trim(filename)

  call LSOPEN(iunit,filename,'unknown','UNFORMATTED')

  !We write the matrix to a binary file
  call write_binary_matrix(Aop%lvec(cell)%oper(oper2),nrows,ncols,iunit)

  call LSCLOSE(iunit,'KEEP')


END SUBROUTINE pbc_get_file_and_write


!!!Writes a matrix, or a vector, if mcol=1, to 
!  an unformatted file
SUBROUTINE write_binary_matrix(A,nrows,mcol,lupri)
IMPLICIT NONE
INTEGER, INTENT(IN) :: nrows, mcol
INTEGER,INTENT(IN) :: lupri
TYPE(matrix), INTENT(IN) :: A
!LOCAL VARIABLES
INTEGER :: i,j,nbast

    nbast=nrows*mcol

    !write(lupri,*) nbast!just debug
    call mat_write_to_disk(lupri,A,.true.)
    !DO j=1,mcol
    !  write(lupri,*) (A%elms(i+(j-1)*mcol),i=1,nrows)
    !ENDDO


END SUBROUTINE write_binary_matrix


SUBROUTINE pbc_read_matrix(Aop,nrows,ncols,oper1,oper2,diis)
IMPLICIT NONE 
TYPE(lvec_list_t),INTENT(INOUT) :: Aop
INTEGER,INTENT(IN) :: nrows,ncols,oper1,oper2
CHARACTER(len=12)          :: diis
!LOCAL Variables
INTEGER :: layer,l1,l2,l3

Do layer=1,size(Aop%lvec)
   l1=int(Aop%lvec(layer)%lat_coord(1))
   l2=int(Aop%lvec(layer)%lat_coord(2))
   l3=int(Aop%lvec(layer)%lat_coord(3))

   select case(oper1)
   
   case(7) !fock
   if(Aop%lvec(layer)%g2_computed .or.Aop%lvec(layer)%f1_computed )then
     !  call mat_init(Aop%lvec(layer)%oper(oper2),nrows,ncols)
       !write(*,*) 'DEBUG get fock',layer, l1
       !write(*,*) 'DEBUG false true',Aop%lvec(layer)%f1_computed,Aop%lvec(layer)%g2_computed
       call mat_zero(Aop%lvec(layer)%oper(oper2))
       call pbc_get_file_and_read(Aop,nrows,ncols,layer,7,oper2,diis)
     endif

   case(1) !overlap
     if(Aop%lvec(layer)%f1_computed)then
        if(Aop%lvec(layer)%oper(oper2)%init_magic_tag .ne. mat_init_magic_value)then
          call mat_init(Aop%lvec(layer)%oper(oper2),nrows,ncols)
        endif
       call mat_zero(Aop%lvec(layer)%oper(oper2))
       call pbc_get_file_and_read(Aop,nrows,ncols,layer,1,oper2,diis)
     endif
   case(2)!kinetic
     if(Aop%lvec(layer)%f1_computed)then
     !  call mat_init(Aop%lvec(layer)%oper(oper2),nrows,ncols)
       call mat_zero(Aop%lvec(layer)%oper(oper2))
       call pbc_get_file_and_read(Aop,nrows,ncols,layer,2,oper2,diis)
     endif
   case(3) !el-nuc
     if(Aop%lvec(layer)%f1_computed)then
     !  call mat_init(Aop%lvec(layer)%oper(oper2),nrows,ncols)
       call mat_zero(Aop%lvec(layer)%oper(oper2))
       call pbc_get_file_and_read(Aop,nrows,ncols,layer,3,1,diis)
     endif
   case(4) !Coulomb
     if(Aop%lvec(layer)%J_computed)then
     !  call mat_init(Aop%lvec(layer)%oper(oper2),nrows,ncols)
       call mat_zero(Aop%lvec(layer)%oper(oper2))
       call pbc_get_file_and_read(Aop,nrows,ncols,layer,4,oper2,diis)
     endif
   case(5) !exact exchange
     if(Aop%lvec(layer)%Kx_computed)then
     !  call mat_init(Aop%lvec(layer)%oper(oper2),nrows,ncols)
       call mat_zero(Aop%lvec(layer)%oper(oper2))
       call pbc_get_file_and_read(Aop,nrows,ncols,layer,5,oper2,diis)
     endif
   case(6) !exchange correlation
     if(Aop%lvec(layer)%Kx_computed)then
     !  call mat_init(Aop%lvec(layer)%oper(oper2),nrows,ncols)
       call mat_zero(Aop%lvec(layer)%oper(oper2))
       call pbc_get_file_and_read(Aop,nrows,ncols,layer,6,oper2,diis)
     endif
   !case(9) !nuc-farfield
   !  if((abs(l1) .le. Aop%oneop1 .and. abs(l2) .le. Aop%oneop2)&
   !  & .and. abs(l3) .le. Aop%oneop3)then
   !  !  call mat_init(Aop%lvec(layer)%oper(oper2),nrows,ncols)
   !    call mat_zero(Aop%lvec(layer)%oper(oper2))
   !    call pbc_get_file_and_read(Aop,nrows,ncols,layer,9,oper2,diis)
   !  endif
   end select

ENDDO

END SUBROUTINE pbc_read_matrix

!ROUTINE for retrieving the fock matrix
!reads each part of the fock matrix from disk
!and adds them together.
SUBROUTINE pbc_read_fock_matrix(Aop,nrows,ncols,diis)
IMPLICIT NONE 
TYPE(lvec_list_t),INTENT(INOUT) :: Aop
INTEGER,INTENT(IN) :: nrows,ncols
CHARACTER(len=12)          :: diis
!LOCAL Variables
INTEGER :: layer,l1,l2,l3,indred
INTEGER :: fdim(3)

Do layer=1,size(Aop%lvec)

   l1=int(Aop%lvec(layer)%lat_coord(1))
   l2=int(Aop%lvec(layer)%lat_coord(2))
   l3=int(Aop%lvec(layer)%lat_coord(3))
   if(.not.Aop%lvec(layer)%is_redundant .or.&
    &             l1**2+l2**2+l3**2 .eq. 0) then

     if(Aop%lvec(layer)%f1_computed .or. Aop%lvec(layer)%g2_computed)then
       call mat_init(Aop%lvec(layer)%oper(1),nrows,ncols)
       call mat_zero(Aop%lvec(layer)%oper(1))
       call mat_init(Aop%lvec(layer)%oper(2),nrows,ncols)
       write(*,*) 'init fock matrix',layer,l1
       call mat_zero(Aop%lvec(layer)%oper(2))
     endif

     if(Aop%lvec(layer)%f1_computed)then
         !overlap
         call mat_zero(Aop%lvec(layer)%oper(1))
         call pbc_get_file_and_read(Aop,nrows,ncols,layer,1,1,"            ")
         !kinetic
         call mat_zero(Aop%lvec(layer)%oper(1))
         call pbc_get_file_and_read(Aop,nrows,ncols,layer,2,1,"            ")
         call mat_daxpy(1.D0,Aop%lvec(layer)%oper(1),Aop%lvec(layer)%oper(2))
         !nuclear attraction near-field
         call mat_zero(Aop%lvec(layer)%oper(1))
         call pbc_get_file_and_read(Aop,nrows,ncols,layer,3,1,"            ")
         call mat_daxpy(1.D0,Aop%lvec(layer)%oper(1),Aop%lvec(layer)%oper(2))
         !nuclear attraction far-field
         !call mat_zero(Aop%lvec(layer)%oper(1))
         !call pbc_get_file_and_read(Aop,nrows,ncols,layer,9,1,"            ")
         !call mat_daxpy(1.D0,Aop%lvec(layer)%oper(1),Aop%lvec(layer)%oper(2))
     endif
     if(Aop%lvec(layer)%J_computed)then
       !Coulomb matrix
       call mat_zero(Aop%lvec(layer)%oper(1))
       call pbc_get_file_and_read(Aop,nrows,ncols,layer,4,1,diis)
       call mat_daxpy(1.D0,Aop%lvec(layer)%oper(1),Aop%lvec(layer)%oper(2))
     endif
     if(Aop%lvec(layer)%kx_computed)then
       !Exact exchange matrix
       call mat_zero(Aop%lvec(layer)%oper(1))
       call pbc_get_file_and_read(Aop,nrows,ncols,layer,5,1,diis)
       call mat_daxpy(1.D0,Aop%lvec(layer)%oper(1),Aop%lvec(layer)%oper(2))
     endif
   endif
   if(Aop%lvec(layer)%is_redundant .and. &
    &             l1**2+l2**2+l3**2 .gt. 0) then
     l1=int(Aop%lvec(layer)%lat_coord(1))
     l2=int(Aop%lvec(layer)%lat_coord(2))
     l3=int(Aop%lvec(layer)%lat_coord(3))
     call find_latt_index(indred,-l1,-l2,-l3,fdim,Aop,Aop%max_layer)
     if(Aop%lvec(indred)%oper(2)%init_magic_tag.EQ.mat_init_magic_value)then
       call mat_init(Aop%lvec(layer)%oper(2),nrows,ncols)
       call mat_zero(Aop%lvec(layer)%oper(2))
       call mat_trans(Aop%lvec(indred)%oper(2),Aop%lvec(layer)%oper(2))
     endif
   endif
ENDDO

END SUBROUTINE pbc_read_fock_matrix

SUBROUTINE pbc_add_fock_matrix(f_1,g_2,ll,nrows,ncols,numvecs)
IMPLICIT NONE 
INTEGER,INTENT(IN) :: nrows,ncols,numvecs
TYPE(lvec_list_t),INTENT(INOUT) :: ll
TYPE(matrix),intent(inout) ::f_1(numvecs),g_2(numvecs) !g_2 is freed so have to input 
!LOCAL Variables
INTEGER :: layer,l1,l2,l3
INTEGER :: g1,g2,g3
INTEGER :: fdim(3),indred

g1=max(ll%col1,ll%Kx1)
g2=max(ll%col2,ll%Kx2)
g3=max(ll%col3,ll%Kx3)

Do layer=1,size(ll%lvec)
   l1=int(ll%lvec(layer)%lat_coord(1))
   l2=int(ll%lvec(layer)%lat_coord(2))
   l3=int(ll%lvec(layer)%lat_coord(3))
   if(.not.ll%lvec(layer)%is_redundant .or. &
                &  l1**2+l2**2+l3**2 .eq. 0) then

     if(ll%lvec(layer)%g2_computed .or.ll%lvec(layer)%f1_computed )then
       call mat_init(ll%lvec(layer)%oper(2),nrows,ncols)
       call mat_zero(ll%lvec(layer)%oper(2))
     endif

     if(ll%lvec(layer)%f1_computed )then !ONE PARTICLE PART
#ifdef DEBUGPBC
         write(*,*) 'layer',layer,l1,'f_1'
         call mat_print(f_1(layer),1,nrows,1,nrows,6)
#endif
         call mat_daxpy(1.D0,f_1(layer),ll%lvec(layer)%oper(2))
     endif
     if(ll%lvec(layer)%g2_computed) THEN !TWO PARTICLE PART
       !write(*,*) 'layer',layer, l1,l2,l3,g_2(layer)%init_magic_tag
#ifdef DEBUGPBC
       write(*,*) 'layer',layer,l1,'g_2'
       call mat_print(g_2(layer),1,nrows,1,nrows,6)
#endif
       call mat_daxpy(1.D0,g_2(layer),ll%lvec(layer)%oper(2))
       call mat_free(g_2(layer))
     endif
   endif
   if(ll%lvec(layer)%is_redundant .and. &
    &             l1**2+l2**2+l3**2 .gt. 0) then
     l1=int(ll%lvec(layer)%lat_coord(1))
     l2=int(ll%lvec(layer)%lat_coord(2))
     l3=int(ll%lvec(layer)%lat_coord(3))
     call find_latt_index(indred,-l1,-l2,-l3,fdim,ll,ll%max_layer)
     if(ll%lvec(indred)%g2_computed .or.ll%lvec(indred)%f1_computed )then
       call mat_init(ll%lvec(layer)%oper(2),nrows,ncols)
       call mat_zero(ll%lvec(layer)%oper(2))
       call mat_trans(ll%lvec(indred)%oper(2),ll%lvec(layer)%oper(2))
     endif
     if(g_2(layer)%init_magic_tag.EQ.mat_init_magic_value)then
       call mat_free(g_2(layer))
     endif
   endif

enddo

END SUBROUTINE pbc_add_fock_matrix

SUBROUTINE pbc_free_read_matrices(Aop)
IMPLICIT NONE 
TYPE(lvec_list_t),INTENT(INOUT) :: Aop
!LOCAL Variables
INTEGER :: layer,l1,l2,l3


Do layer=1,size(Aop%lvec)
   l1=int(Aop%lvec(layer)%lat_coord(1))
   l2=int(Aop%lvec(layer)%lat_coord(2))
   l3=int(Aop%lvec(layer)%lat_coord(3))

   if((abs(l1) .le. Aop%fc1 .and. abs(l2) .le. Aop%fc2) .and. abs(l3) .le. Aop%fc3)then
     if(Aop%lvec(layer)%oper(1)%init_magic_tag.EQ.mat_init_magic_value)then
       call mat_free(Aop%lvec(layer)%oper(1))
     endif
     if(Aop%lvec(layer)%oper(2)%init_magic_tag.EQ.mat_init_magic_value)then
       call mat_free(Aop%lvec(layer)%oper(2))
     endif
   endif

enddo

END SUBROUTINE pbc_free_read_matrices

!SUBROUTINE to get filename for the operator, open the file and read 
!from disk
SUBROUTINE pbc_get_file_and_read(Aop,nrows,ncols,cell,oper1,oper2,diis)
IMPLICIT NONE 
TYPE(lvec_list_t),INTENT(INOUT) :: Aop
INTEGER,INTENT(IN) :: nrows,ncols,cell,oper1,oper2
CHARACTER(len=12)          :: diis
!LOCAL VARIABLES
INTEGER :: l1,l2,l3
CHARACTER(len=7) :: stl1,stl2,stl3
CHARACTER(LEN=100) :: filename
INTEGER :: iunit,FILELEN,IOS
LOGICAL :: fileexists

  l1=int(Aop%lvec(cell)%lat_coord(1))
  l2=int(Aop%lvec(cell)%lat_coord(2))
  l3=int(Aop%lvec(cell)%lat_coord(3))
  write(stl1,'(I5)')  l1
  write(stl2,'(I5)')  l2
  write(stl3,'(I5)')  l3
  stl1=adjustl(stl1)
  stl2=adjustl(stl2)
  stl3=adjustl(stl3)

  iunit=-1

! write(*,*) 'debug 3',oper1
  !    Now we are doing DIIS scf iterations and we read the stored
  !    matrices 
  filename = trim(diis)//trim(adjustl(Aop%opdat(oper1)%basename))//'_'//trim(stl1)//'_'//trim(stl2)//'_'//trim(stl3)
  
  filename=adjustl(filename)
  filename=trim(filename)
  FILELEN=len(filename)
  INQUIRE(FILE=filename(1:FILELEN),EXIST=fileexists,IOSTAT=IOS)
  IF(fileexists) then
    call LSOPEN(iunit,filename,'OLD','UNFORMATTED')

    !We write the matrix to a binary file
    call read_binary_matrix(Aop%lvec(cell)%oper(oper2),nrows,ncols,iunit)

    call LSCLOSE(iunit,'KEEP')
  ENDIF


END SUBROUTINE pbc_get_file_and_read




!!! READS a matrix or vector from unformatted file
SUBROUTINE read_binary_matrix(A,nrows,mcol,lupri)
IMPLICIT NONE
INTEGER, INTENT(IN) :: nrows, mcol
INTEGER,INTENT(IN) :: lupri
TYPE(matrix), INTENT(INOUT) :: A
!LOCAL VARIABLES
INTEGER :: i,j,nbast

    nbast=nrows*mcol

    rewind(lupri)
    call mat_read_from_disk(lupri,A,.true.)
    !DO j=1,mcol
    !  read(lupri,*) (A(i+(j-1)*mcol),i=1,nrows)
    !ENDDO


END SUBROUTINE read_binary_matrix


!Creates an upper triagonal matrix A consisting of ones.
!I do not know why I made it.
SUBROUTINE make_uptriag1_mat(A,nrow,ncol)
IMPLICIT NONE
INTEGER, INTENT(IN) :: nrow, ncol
REAL(realk), INTENT(INOUT)    :: A(nrow,ncol)
!local variables
INTEGER             :: i
A=0.0
DO i=1, nrow-1
   A(i,i+1)=1.0
ENDDO

END SUBROUTINE make_uptriag1_mat
END MODULE pbc_matrix_operations
#endif
