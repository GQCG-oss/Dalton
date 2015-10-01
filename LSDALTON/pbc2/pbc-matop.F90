MODULE pbc_matrix_operations
	USE precision
	USE typedef
	USE lattice_type
	USE lattice_vectors
	USE files
	USE fundamental

CONTAINS

!> \author johannes rekkedal
!> \date 2013
!> \brief this routine writes a matrix to screen in a consistent form.
!> \brief Writes to print unit lupri if lupri is present.
!> \param a 		matrix to be printed.
!> \param nrows 	mtrx dim.
!> \param ncols  	mtrx dim.
!> \param lupri 	logical print unit for output file. (optional)
SUBROUTINE write_matrix(A,nrows,mcol,lupri)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: nrows, mcol
	INTEGER,INTENT(IN),OPTIONAL :: lupri
	REAL(realk), INTENT(IN) :: A(nrows,mcol)
	! local
	INTEGER :: i,j,k
	CHARACTER(LEN=3) :: nline

	if (PRESENT(lupri)) then
		k=0
		do i=1,nrows
			do j=1,mcol
				k=k+1
				nline='no'
				if(k .eq. mcol) then
					k=0
					nline='yes'
				endif
				write(lupri,'(E25.12)',advance=nline) A(i,j)
			enddo
		enddo
	else
		k=0
		do i=1,nrows
			do j=1,mcol
				k=k+1
				nline='no'
				if(k .eq. mcol) then
					k=0
					nline='yes'
				endif
				write(*,'(E25.12)',advance=nline)  A(i,j)
			enddo
		enddo
	endif

END SUBROUTINE write_matrix

!> \author johannes rekkedal
!> \date 2013
!> \brief This routine writes a complex matrix to screen in a consistent form.
!> \brief Writes to print unit lupri if lupri is present.
!> \param a 		matrix to be printed.
!> \param nrows 	mtrx dim.
!> \param ncols  	mtrx dim.
!> \param lupri 	logical print unit for output file. (optional)
SUBROUTINE write_zmatrix(A,nrows,mcol,lupri)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: nrows, mcol
	COMPLEX(COMPLEXK), INTENT(IN) :: A(nrows,mcol)
	INTEGER, INTENT(IN),OPTIONAL :: lupri
	! local
	INTEGER :: i,j,k
	CHARACTER(LEN=3) :: nline

	if (PRESENT(lupri)) then
		k=0
		do i=1,nrows
			do j=1,mcol
				k=k+1
				nline='no'
				if (k .eq. mcol) then
					k=0
					nline='yes'
				endif
				write(lupri,101,advance=nline) A(i,j)
				101 FORMAT('(',2E25.12,')')
			enddo
		enddo
	else
		k=0
		do i=1,nrows
			do j=1,mcol
				k=k+1
				nline='no'
				if (k .eq. mcol) then
					k=0
					nline='yes'
				endif
				write(*,100,advance=nline) A(i,j)
				100 format('(',2E25.12,')')
			enddo
		enddo
	endif

END SUBROUTINE write_zmatrix

!> \author johannes rekkedal
!> \date 2013
!> \brief This routine writes the Fock matrix to disk.
!> \param Aop 		Lattice info including Fock matrix.
!> \param nrows 	mtrx dim.
!> \param ncols  	mtrx dim.
!> \param oper1  	mtrx ???.
!> \param oper2  	mtrx ???.
!> \param diis 	???
!> \param lu 	 	logical print unit for output file.
SUBROUTINE pbc_fockmat_write(Aop,nrows,ncols,oper1,oper2,diis,lu)
	IMPLICIT NONE 
	TYPE(lvec_list_t),INTENT(IN) :: Aop
	INTEGER,INTENT(IN) :: nrows,ncols,oper1,oper2,lu
	CHARACTER(LEN=12)  :: diis
	! local
	INTEGER :: layer
	INTEGER :: l1,l2,l3

	do layer =1,size(Aop%lvec)
		if (Aop%lvec(layer)%f1_computed .or. Aop%lvec(layer)%g2_computed) then
			l1=int(Aop%lvec(layer)%lat_coord(1))
			l2=int(Aop%lvec(layer)%lat_coord(2))
			l3=int(Aop%lvec(layer)%lat_coord(3))
			call pbc_get_file_and_write(Aop,nrows,ncols,layer,oper1,oper2,diis)
		endif
	enddo

END SUBROUTINE pbc_fockmat_write

!> \author johannes rekkedal
!> \date 2013
!> \brief This routine writes the density matrix to disk.
!> \param dmat 	Density matrix.
!> \param Aop 		Lattice info.
!> \param nrows 	mtrx dim.
!> \param ncols  	mtrx dim.
!> \param oper1  	???.
!> \param diis 	???
SUBROUTINE pbc_densitymat_write(dmat,Aop,nrows,ncols,oper1,diis)
	IMPLICIT NONE 
	TYPE(lvec_list_t),INTENT(IN) :: Aop
	TYPE(matrix),INTENT(IN) :: dmat(size(Aop%lvec))
	INTEGER,INTENT(IN) :: nrows,ncols,oper1
	CHARACTER(len=12)          :: diis
	! local
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

		if(Aop%lvec(layer)%dm_computed)then
			! Now we are doing DIIS scf iterations and we need to store 
			! the matrix for some of the iterations
			filename = trim(diis)//trim(adjustl(Aop%opdat(oper1)%basename))//'_'//trim(stl1)//'_'//trim(stl2)//'_'//trim(stl3)

			filename=adjustl(filename)
			filename=trim(filename)
			!We write the matrix to a binary file
			call LSOPEN(iunit,filename,'unknown','UNFORMATTED')
			call write_binary_matrix(dmat(layer),nrows,ncols,iunit)
			call LSCLOSE(iunit,'KEEP')

		endif
	enddo

END SUBROUTINE pbc_densitymat_write

!author Johannes Rekkedal
!writes max density element from each layer
SUBROUTINE print_maxdens(dmat,ll,lupri)
IMPLICIT NONE
	TYPE(lvec_list_t),INTENT(IN) :: ll
	TYPE(matrix),INTENT(IN) :: dmat(size(ll%lvec))
        INTEGER,INTENT(IN)      :: lupri
        !LOCAL
        INTEGER                 :: layer,l1,l2,l3
        real(realk)             :: maxdens

  DO layer = 1,size(ll%lvec)
     l1=int(ll%lvec(layer)%lat_coord(1))
     l2=int(ll%lvec(layer)%lat_coord(2))
     l3=int(ll%lvec(layer)%lat_coord(3))
     if(ll%lvec(layer)%dm_computed) then
       call mat_abs_max_elm(dmat(layer),maxdens)
       write(lupri,*) maxdens,l1,l2,l3
     endif
  END DO


END SUBROUTINE print_maxdens

! todo for debug purpuses (delete??)
!> \author johannes rekkedal
!> \date 2013
!> \brief This function prints a T or W matrix in a format that allows it
!> \brief to be cut-n-pasted into MATLAB. It is for debugging purposes only.
!> \param matrix  	Mtrx to be printed.
!> \param maxrows 	Mtrx dim.
!> \param maxcols  	Mtrx dim.
!> \param text  		Name of matrix.
!> \param lupri 		Logical print unit.
SUBROUTINE pbc_matlab_print(matrix,maxrows,maxcols,text,lupri)
	IMPLICIT NONE
	CHARACTER*(*) :: text
	INTEGER,INTENT(IN) :: maxrows, maxcols,lupri
	REAL(realk) :: matrix(maxrows,maxcols)
	! local
	INTEGER :: row, col, col_max

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

!> \author johannes rekkedal
!> \date 2013
!> \brief Get filename for the operator, open the file and write.
!> \brief to disk
!> \param Aop  		Lattice info incl. operators (matrices).
!> \param nrows 		Mtrx dim.
!> \param ncols  		Mtrx dim.
!> \param cell  		Index of cell.
!> \param oper1  		???
!> \param oper2  		???
!> \param diis 		???
SUBROUTINE pbc_get_file_and_write(Aop,nrows,ncols,cell,oper1,oper2,diis)
	IMPLICIT NONE 
	TYPE(lvec_list_t),INTENT(IN) :: Aop
	INTEGER,INTENT(IN) :: nrows,ncols,cell,oper1,oper2
	CHARACTER(len=12)  :: diis
	! local
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
	call lsclose(iunit,'KEEP')

END SUBROUTINE pbc_get_file_and_write

! todo Really necc?? only calls an other routine. delete.
!> \author johannes rekkedal
!> \date 2013
! \brief Writes a matrix, or a vector, if mcol=1, to an unformatted file
!> \param A 	 		Mtrx to be printesd.
!> \param nrows 		Mtrx dim.
!> \param ncols  		Mtrx dim.
!> \param lupri 		Logical print unit.
SUBROUTINE write_binary_matrix(A,nrows,mcol,lupri)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: nrows, mcol
	INTEGER,INTENT(IN) :: lupri
	TYPE(matrix), INTENT(IN) :: A

	call mat_write_to_disk(lupri,A,.true.)

END SUBROUTINE write_binary_matrix

!> \author johannes rekkedal
!> \date 2013
! \brief Read a matrix, or a vector, if mcol=1, to an unformatted file
!> \param A 	 		Mtrx to be printesd.
!> \param nrows 		Mtrx dim.
!> \param ncols  		Mtrx dim.
!> \param oper1  		???
!> \param oper2  		???
!> \param diis 		???
SUBROUTINE pbc_read_matrix(Aop,nrows,ncols,oper1,oper2,diis)
	IMPLICIT NONE 
	TYPE(lvec_list_t),INTENT(INOUT) :: Aop
	INTEGER,INTENT(IN) :: nrows,ncols,oper1,oper2
	CHARACTER(LEN=12) :: diis
	! local
	INTEGER :: layer,l1,l2,l3

	Do layer=1,size(Aop%lvec)
		l1=int(Aop%lvec(layer)%lat_coord(1))
		l2=int(Aop%lvec(layer)%lat_coord(2))
		l3=int(Aop%lvec(layer)%lat_coord(3))

		select case(oper1)

		case(7) !fock
			if(Aop%lvec(layer)%g2_computed .or.Aop%lvec(layer)%f1_computed )then
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
		case(2) !kinetic
			if(Aop%lvec(layer)%f1_computed)then
				call mat_zero(Aop%lvec(layer)%oper(oper2))
				call pbc_get_file_and_read(Aop,nrows,ncols,layer,2,oper2,diis)
			endif
		case(3) !el-nuc
			if(Aop%lvec(layer)%f1_computed)then
				call mat_zero(Aop%lvec(layer)%oper(oper2))
				call pbc_get_file_and_read(Aop,nrows,ncols,layer,3,1,diis)
			endif
		case(4) !coulomb
			if(Aop%lvec(layer)%J_computed)then
				call mat_zero(Aop%lvec(layer)%oper(oper2))
				call pbc_get_file_and_read(Aop,nrows,ncols,layer,4,oper2,diis)
			endif
		case(5) !exact exchange
			if(Aop%lvec(layer)%Kx_computed)then
				call mat_zero(Aop%lvec(layer)%oper(oper2))
				call pbc_get_file_and_read(Aop,nrows,ncols,layer,5,oper2,diis)
			endif
		case(6) !exchange correlation
			if(Aop%lvec(layer)%Kx_computed)then
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

!> \author johannes rekkedal
!> \date 2013
!> \brief Retrieving the fock matrix reads each part of the fock matrix from 
!> \brief disk and adds them together.
!> \param Aop 			Lattice info incl. operators (matrices).
!> \param nrows 		Mtrx dim.
!> \param ncols 		Mtrx dim.
!> \param diis 		???
SUBROUTINE pbc_read_fock_matrix(Aop,nrows,ncols,diis)
	IMPLICIT NONE 
	TYPE(lvec_list_t),INTENT(INOUT) :: Aop
	INTEGER,INTENT(IN) :: nrows,ncols
	CHARACTER(LEN=12) :: diis
	! local
	INTEGER :: layer,l1,l2,l3,indred

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
				call mat_daxpy(1.0_realk,Aop%lvec(layer)%oper(1),Aop%lvec(layer)%oper(2))
				!nuclear attraction near-field
				call mat_zero(Aop%lvec(layer)%oper(1))
				call pbc_get_file_and_read(Aop,nrows,ncols,layer,3,1,"            ")
				call mat_daxpy(1.0_realk,Aop%lvec(layer)%oper(1),Aop%lvec(layer)%oper(2))
				!nuclear attraction far-field
				!call mat_zero(Aop%lvec(layer)%oper(1))
				!call pbc_get_file_and_read(Aop,nrows,ncols,layer,9,1,"            ")
				!call mat_daxpy(1.0_realk,Aop%lvec(layer)%oper(1),Aop%lvec(layer)%oper(2))
			endif
			if(Aop%lvec(layer)%J_computed)then
				!Coulomb matrix
				call mat_zero(Aop%lvec(layer)%oper(1))
				call pbc_get_file_and_read(Aop,nrows,ncols,layer,4,1,diis)
				call mat_daxpy(1.0_realk,Aop%lvec(layer)%oper(1),Aop%lvec(layer)%oper(2))
			endif
			if(Aop%lvec(layer)%kx_computed)then
				!Exact exchange matrix
				call mat_zero(Aop%lvec(layer)%oper(1))
				call pbc_get_file_and_read(Aop,nrows,ncols,layer,5,1,diis)
				call mat_daxpy(1.0_realk,Aop%lvec(layer)%oper(1),Aop%lvec(layer)%oper(2))
			endif
		endif
		if(Aop%lvec(layer)%is_redundant .and. l1**2+l2**2+l3**2 .gt. 0) then
			l1=int(Aop%lvec(layer)%lat_coord(1))
			l2=int(Aop%lvec(layer)%lat_coord(2))
			l3=int(Aop%lvec(layer)%lat_coord(3))
			call find_latt_index(indred,-l1,-l2,-l3,Aop,Aop%max_layer)
			if(Aop%lvec(indred)%oper(2)%init_magic_tag.EQ.mat_init_magic_value)then
				call mat_init(Aop%lvec(layer)%oper(2),nrows,ncols)
				call mat_zero(Aop%lvec(layer)%oper(2))
				call mat_trans(Aop%lvec(indred)%oper(2),Aop%lvec(layer)%oper(2))
			endif
		endif
	enddo

END SUBROUTINE pbc_read_fock_matrix

!> \author johannes rekkedal
!> \date 2013
!> \brief ??? 
!> \param ??? 		???
SUBROUTINE pbc_add_fock_matrix(f_1,g_2,ll,nrows,ncols,numvecs)
	IMPLICIT NONE 
	INTEGER,INTENT(IN) :: nrows,ncols,numvecs
	TYPE(lvec_list_t),INTENT(INOUT) :: ll
	!g_2 is freed so have to input 
	TYPE(matrix),INTENT(INOUT) ::f_1(numvecs),g_2(numvecs) 
	! local
	INTEGER :: layer,l1,l2,l3
	INTEGER :: g1,g2,g3
	INTEGER :: indred

	g1=max(ll%col1,ll%Kx1)
	g2=max(ll%col2,ll%Kx2)
	g3=max(ll%col3,ll%Kx3)

	Do layer=1,size(ll%lvec)
		l1=int(ll%lvec(layer)%lat_coord(1))
		l2=int(ll%lvec(layer)%lat_coord(2))
		l3=int(ll%lvec(layer)%lat_coord(3))
		if(.not.ll%lvec(layer)%is_redundant .or. l1**2+l2**2+l3**2 .eq. 0) then
			if(ll%lvec(layer)%g2_computed .or.ll%lvec(layer)%f1_computed )then
				call mat_init(ll%lvec(layer)%oper(2),nrows,ncols)
				call mat_zero(ll%lvec(layer)%oper(2))
			endif
			if(ll%lvec(layer)%f1_computed ) then !ONE PARTICLE PART
				call mat_daxpy(1.0_realk,f_1(layer),ll%lvec(layer)%oper(2))
			endif
			if(ll%lvec(layer)%g2_computed) then !two particle part
				call mat_daxpy(1.0_realk,g_2(layer),ll%lvec(layer)%oper(2))
				call mat_free(g_2(layer))
			endif
		endif
		if(ll%lvec(layer)%is_redundant .and. l1**2+l2**2+l3**2 .gt. 0) then
			l1=int(ll%lvec(layer)%lat_coord(1))
			l2=int(ll%lvec(layer)%lat_coord(2))
			l3=int(ll%lvec(layer)%lat_coord(3))
			call find_latt_index(indred,-l1,-l2,-l3,ll,ll%max_layer)
			if(ll%lvec(indred)%g2_computed .or.ll%lvec(indred)%f1_computed )then
				call mat_init(ll%lvec(layer)%oper(2),nrows,ncols)
				call mat_zero(ll%lvec(layer)%oper(2))
				call mat_trans(ll%lvec(indred)%oper(2),ll%lvec(layer)%oper(2))
			endif
			if(g_2(layer)%init_magic_tag.eq.mat_init_magic_value)then
				call mat_free(g_2(layer))
			endif
		endif
	enddo

END SUBROUTINE pbc_add_fock_matrix

!> \author johannes rekkedal
!> \date 2013
!> \brief ??? 
!> \param ??? 		???
SUBROUTINE pbc_free_read_matrices(Aop)
	IMPLICIT NONE 
	TYPE(lvec_list_t),INTENT(INOUT) :: Aop
	! local
	INTEGER :: layer,l1,l2,l3


	do layer=1,size(Aop%lvec)
		l1=int(Aop%lvec(layer)%lat_coord(1))
		l2=int(Aop%lvec(layer)%lat_coord(2))
		l3=int(Aop%lvec(layer)%lat_coord(3))

		if ( abs(l1).le.Aop%fc1 & 
			& .and.abs(l2).le.Aop%fc2 & 
			& .and.abs(l3).le.Aop%fc3) then
			if(Aop%lvec(layer)%oper(1)%init_magic_tag.eq.mat_init_magic_value)then
				call mat_free(Aop%lvec(layer)%oper(1))
			endif
			if(Aop%lvec(layer)%oper(2)%init_magic_tag.eq.mat_init_magic_value)then
				call mat_free(Aop%lvec(layer)%oper(2))
			endif
		endif
	enddo

END SUBROUTINE pbc_free_read_matrices

!> \author johannes rekkedal
!> \date 2013
!> \brief Get filename for the operator, open the file and read 
!> \brief from disk
!> \param Aop  		Lattice info incl. operators (matrices).
!> \param nrows 		Mtrx dim.
!> \param ncols  		Mtrx dim.
!> \param cell  		Index of cell.
!> \param oper1  		???
!> \param oper2  		???
!> \param diis 		???
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

	! now we are doing diis scf iterations and we read the stored matrices 
	filename = trim(diis)//trim(adjustl(Aop%opdat(oper1)%basename))//'_'//trim(stl1)//'_'//trim(stl2)//'_'//trim(stl3)

	filename=adjustl(filename)
	filename=trim(filename)
	filelen=len(filename)
	INQUIRE(FILE=filename(1:FILELEN),EXIST=fileexists,IOSTAT=IOS)
	if (fileexists) then
		call LSOPEN(iunit,filename,'OLD','UNFORMATTED')

		!We write the matrix to a binary file
		call read_binary_matrix(Aop%lvec(cell)%oper(oper2),nrows,ncols,iunit)
		call lsclose(iunit,'KEEP')
	endif

END SUBROUTINE pbc_get_file_and_read

! todo really necc?
!> \author johannes rekkedal
!> \date 2013
!> \brief reads a matrix or vector from unformatted file
!> \param A 			Matrix outfut. 
!> \param nrows 		Mtrx dim.
!> \param ncols 		Mtrx dim.
!> \param lupri 		Logical print unit.
SUBROUTINE read_binary_matrix(A,nrows,mcol,lupri)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: nrows, mcol
	INTEGER,INTENT(IN) :: lupri
	TYPE(matrix), INTENT(INOUT) :: A

	rewind(lupri)
	call mat_read_from_disk(lupri,A,.true.)

END SUBROUTINE read_binary_matrix


!> \author johannes rekkedal
!> \date 2013
!> \brief Creates an upper triagonal matrix A consisting of ones.
!> \param A 		Matrix.
!> \param nrow 	Mtrx dim.
!> \param nrow 	Mtrx dim.
SUBROUTINE make_uptriag1_mat(A,nrow,ncol)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: nrow, ncol
	REAL(realk), INTENT(INOUT)    :: A(nrow,ncol)
	! local
	INTEGER :: i
	a=0.0
	do i=1, nrow-1
		a(i,i+1)=1.0
	enddo

END SUBROUTINE make_uptriag1_mat

!! todo not in use
!!> \author 	Karl R. Leikanger
!!> \date 	2013
!!> \brief 	Orthogonalize complex matrix c
!!> \param 	c 		mxn matrix to be orthogonalized
!!> \param 	m 		mtrx dim
!!> \param 	n 		mtrx dim
!!> \param 	lupri logical print unit
!SUBROUTINE pbcmo_orthogonalize(c,m,n,lupri)
!	IMPLICIT NONE
!	INTEGER, INTENT(IN) :: m,n,lupri
!	COMPLEX(complexk), INTENT(INOUT) :: c(m,n)
!	!local
!	INTEGER :: lda, info, lwork
!	COMPLEX(complexk), POINTER :: work(:), tau(:)
!
!	lda = m
!	lwork = 64*n
!
!	call mem_alloc(tau, min(m,n))
!	call mem_alloc(work, lwork)
!
!	call zgeqrf(m,n,c,lda,tau,work,lwork,info)
!	if (info .lt. 0) then
!		write (*,*) 'error zgeqrf in pbcmo_orthogonalize, info=', info
!		write (lupri,*) 'error zgeqrf in pbcmo_orthogonalize, info=', info
!	end if
!	call zungqr(m,n,n,c,lda,tau,work,lwork,info)
!	if (info .lt. 0) then
!		write (*,*) 'error zungqr in pbcmo_orthogonalize, info=', info
!		write (lupri,*) 'error zungqr in pbcmo_orthogonalize, info=', info
!	end if
!END SUBROUTINE pbcmo_orthogonalize


!This routine makes sure that C^dagger S C = I
!SUBROUTINE zggram_schmidt(C,S,nrowsc,mcolsc,nrowss,mcolss,lupri)
!IMPLICIT NONE
!INTEGER, INTENT(IN) :: nrowsc, mcolsc,nrowss,mcolss
!COMPLEX(complexk), INTENT(INOUT) :: C(nrowsc,mcolsc)
!COMPLEX(complexk), INTENT(IN) :: S(nrowss,mcolss)
!INTEGER, INTENT(IN), OPTIONAL :: lupri
!!local variables
!INTEGER :: i,j,k
!COMPLEX(complexk) :: alpha,beta
!COMPLEX(complexk),pointer :: ctmp1(:,:),ctmp2(:,:),vtmp1(:,:)
!COMPLEX(complexk),pointer :: vtmp2(:,:)
!COMPLEX(complexk),external :: zdotc
!REAL(realk) :: radii,angle
!COMPLEX(complexk) :: factor1,factor2,fac
!
! alpha=CMPLX(1._realk,0._realk,complexk)
! beta=CMPLX(0._realk,0._realk,complexk)
! factor1=CMPLX(0._realk,0._realk,complexk)
! factor2=CMPLX(0._realk,0._realk,complexk)
!
! !vtmp1(:)=C(:,1)
!
! call mem_alloc(ctmp1,nrowsc,1)
! call mem_alloc(ctmp2,nrowsc,1)
! call mem_alloc(vtmp1,nrowsc,1)
! call mem_alloc(vtmp2,nrowsc,1)
!
!! call zgemm('C','N',nrowsc,mcolss,1,alpha,vtmp1,nrowsc,S,nrowss,nrowss,beta, &
!!           & vtmp2,nrowsc)
! !call zgemm('N','N',nrowsc,mcolss,1,alpha,vtmp2,nrowsc,S,nrowss,nrowss,beta &
! !          & vtmp1,nrowsc)
!
! !c(:,1)=vtmp1(:)!/zdotc(nrowsc,vtmp2,1,vtmp1,1)
! !c(:,1)=vtmp1(:)
!
! !Orthogonalize the eigenvector C_i S C_j= 0
! do i=2,mcolsc
!    vtmp1(:)=c(:,i)
!    vtmp2(:)=c(:,i)
!    call zgemm('C','N',1,mcolss,nrowsc,alpha,vtmp2,nrowss,S,nrowss,beta, &
!            & ctmp1,1)
!    do j=1,i-1
!
!      beta=CMPLX(0._realk,0._realk,complexk)
!      !ctmp1(:)=c(:,j)
!      ctmp2(:)=c(:,j)
!
!      factor1=zdotc(nrowsc,ctmp1,1,ctmp2,1)
!
!      call zgemm('C','N',1,mcolss,nrowsc,alpha,ctmp2,nrowsc,S,nrowss,beta, &
!            & ctmp1,1)
!
!      factor2=zdotc(nrowsc,ctmp1,1,ctmp2,1)
!
!      vtmp1=vtmp1-ctmp2*factor2/factor1
!
!    enddo
!    c(:,i)=vtmp1(:)
! enddo
!
! !Normalize C_i S C_i = 1
!
! Do i=1,mcolsc
!    vtmp1(:)=c(:,i)
!    vtmp2(:)=c(:,i)
!    call zgemm('C','N',1,mcolss,nrowsc,alpha,vtmp2,nrowsc,S,nrowss,beta, &
!            & ctmp1,1)
!    factor1=sqrt(zdotc(nrowsc,ctmp1,1,vtmp2,1))
!    c(:,i)=vtmp1(:)/factor1
! enddo
!
! call mem_dealloc(ctmp1)
! call mem_dealloc(ctmp2)
! call mem_dealloc(vtmp1)
! call mem_dealloc(vtmp2)
!
!END SUBROUTINE zggram_schmidt
END MODULE pbc_matrix_operations
