!> @file
!> Contains sparse1 (compressed sparse row, CSR) matrix module.

!----------------------------------------------------------------------
!---------- sparse matrix class (linked-list implementation) ----------
!---------- Andreas Hesselmann (2005)                        ----------

!----------------------------------------------------------------------
module row_entry
!----------------------------------------------------------------------
  type row_item
    integer :: j
    real(8) :: value 
  end type row_item

  interface operator (<) ! for sorting
    module procedure less_than_row  
  end interface
  interface operator (==) ! for sorting
    module procedure equal_to_row
  end interface

contains                                  ! overload definitions only
  function less_than_row(row1, row2) result(Boolean)
    type(row_item), intent(in) :: row1, row2
    logical                    :: Boolean
    Boolean = row1%j < row2%j   ! standard (<) here
  end function less_than_row
  function equal_to_row(row1, row2) result (Boolean)
    type(row_item), intent(in) :: row1, row2
    logical                    :: Boolean
    Boolean = row1%j == row2%j  ! standard (==) here
  end function equal_to_row
end module row_entry


!----------------------------------------------------------------------
module row
!----------------------------------------------------------------------
  use row_entry
  implicit none

  type :: row_node                        ! Singly Linked Node
     type(row_item)          :: row       ! Object attribute
     type(row_node), pointer :: next      ! Pointer to next node
  end type row_node

  type row_list                       ! Singly Linked List of Nodes
     type(row_node), pointer :: first ! Dummy first object in list
  end type row_list
 

contains
  
  subroutine row_delete(links, Obj, found)
    type(row_list), intent (inout) :: links
    type(row_item), intent (in)    :: Obj
    logical,        intent (out)   :: found
    type(row_node), pointer        :: previous, current

    ! find location of Obj
    previous => links%first        ! begin at top of list
    if ( .not. associated (previous) ) print *, 'NA first'   
    current => previous%next
    if ( .not. associated ( current) ) print *, 'NA current ' 
    found = .false.
    do
       if( found .or. (.not. associated (current)) ) return ! end of list
       if( Obj == current%row ) then
          found = .true. ; exit ! this location search
       else ! move the next node in list
          previous => previous%next 
          current  => current%next
       end if
    end do ! to find location of node with Obj
    ! delete if found
    if ( found )  then
       previous%next => current%next  ! redirect pointer
       deallocate ( current )         ! free space for node
    end if
  end subroutine row_delete
   
  
  subroutine row_insert(links, Obj)
    type(row_list), intent(inout)  :: links
    type(row_item), intent(in)     :: Obj
    type(row_node), pointer        :: previous, current
    integer                        :: istat

    if ( .not. associated(links%first) ) then   
       print *, 'allocating first'              
       allocate (links%first, stat=istat)       
       if (istat > 0) print *, 'allocate failed'
    end if    
    !  Find location to insert a new object
    previous => links%first
    current => previous%next
    if ( .not. associated (current) ) then ! insert at end 
       allocate ( previous%next )      ! get new node space
       previous%next%row = Obj         ! new object inserted
       previous%next%next => current   ! new next pointer
       return
    endif

    do
       if ( .not. associated (current) ) exit ! insert at end 
       if ( Obj < current%row ) exit      ! insert before current
       previous => current                ! move to next node
       current  => current%next           ! move to next node
    end do ! to locate insert node
    !  Insert before current (duplicates not allowed)
    if( Obj == previous%row ) then
       previous%row%value=Obj%value
    else 
       allocate ( previous%next )      ! get new node space
       previous%next%row = Obj         ! new object inserted
       previous%next%next => current   ! new next pointer
    endif
  end subroutine row_insert
    

  function is_row_empty(links) result(t_or_f)
    type(row_list), intent(in) :: links
    logical                    :: t_or_f
    t_or_f = .not. associated ( links%first%next )
  end function is_row_empty  


  function row_new() result(new_list)
    type(row_list) :: new_list
    allocate( new_list%first )        ! get memory for the object
    nullify( new_list%first%next )    ! begin with empty list
  end function row_new


  function number_of_nodes(links) result(number)
    type(row_list), intent(in) :: links
    integer :: number
    type(row_node), pointer :: current
    current => links%first%next
    number=0
    do 
       if ( .not. associated (current) ) exit ! list end
       number=number+1
       current => current%next
    enddo
  end function number_of_nodes
  

  subroutine print_row_list(links)
    type(row_list), intent(in) :: links
    type(row_node), pointer    :: current
    integer                    :: counter
    current => links%first%next
    counter = 0 
    write(*,*) 'Link    Object Value'
    do 
       if ( .not. associated (current) ) exit ! list end
       counter = counter + 1 
       !print *, counter, '     ', current%row
       write(*,'(i4,2x,i6,f16.8)') counter,current%row%j,current%row%value
       current => current%next
    end do
  end subroutine print_row_list
  
end module row


!----------------------------------------------------------------------
module sparse_matrix
!----------------------------------------------------------------------
  use row
  implicit none

  type spmatrix
     integer :: ncol
     TYPE(row_list), DIMENSION(:),pointer :: col
  end type spmatrix
  
contains

  subroutine spmatrix_construct(a,n)
    integer, intent(in)           :: n
    type(spmatrix), intent(inout) :: a
    integer :: i

    a%ncol=n
    allocate(a%col(n))
    !initialise the row linked lists
    do i=1,n
       a%col(i)=row_new()
    enddo
  end subroutine spmatrix_construct
  

  subroutine spmatrix_destruct(a)
    type(spmatrix) :: a

    deallocate(a%col)
  end subroutine spmatrix_destruct


  function spmatrix_number_of_elements(a) result(number)
    type(spmatrix), intent(in) :: a
    integer :: number,i

    number=0
    do i=1,a%ncol
       number=number+number_of_nodes(a%col(i))
    enddo
  end function spmatrix_number_of_elements


  subroutine spmatrix_print(a)
    type(spmatrix), intent(in) :: a
    integer :: i

    write(*,*)
    do i=1,a%ncol
       write(*,'(1x,a,i6)') 'i=',i
       call print_row_list(a%col(i))
    enddo
  end subroutine spmatrix_print
  
end module sparse_matrix

!---------- end of sparse matrix class --------------------------------
!----------------------------------------------------------------------


!> \brief Contains matrix operation routines for type(matrix) = sparse1
module matrix_operations_sparse1
  use matrix_module
  use memory_handling
  private :: nidata, nelm, symm, pack, asymm, usymm, unzipd, zipd, zero, zero2,waste_frac
  private :: iopt, jopt, ijopt, jiopt, sort, nosort
  !> Size of the integer array idata in type(matrix)
  integer, parameter :: nidata = 3
  !> First position in idata is nelm = number of non-zero elements in the sparse matrix
  integer, parameter :: nelm = 1 
  !> Second position in idata is symm; indicates whether the sparse matrix has any symmetry
  integer, parameter :: symm = 2 
  !> Third position in idata is pack; indicates whether the sparse matrix is packed
  integer, parameter :: pack = 3
  !> idata(2) = idata(symm) is set to 1 if sparse matrix is asymmetric
  integer, parameter :: asymm = 1 
  !> idata(2) = idata(symm) is set to 3 if sparse matrix is nonsymmetric
  integer, parameter :: usymm = 3 
  !> idata(3) = idata(pack) is set to 1 if sparse matrix is unzipped
  integer, parameter :: unzipd = 1
  !> idata(3) = idata(pack) is set to 2 if sparse matrix is zipped
  integer, parameter :: zipd = 2
  !> If an array has nelm non-zero elements, we allow it to occupy at most nelm*waste_frac space
  real(realk), parameter :: waste_frac = 1.3
  !> See documentation for sp_sort in matop_sparse1.f90
  integer, parameter :: iopt = 1 
  !> See documentation for sp_sort in matop_sparse1.f90
  integer, parameter :: jopt = 2 
  !> See documentation for sp_sort in matop_sparse1.f90
  integer, parameter :: ijopt = 3 
  !> See documentation for sp_sort in matop_sparse1.f90
  integer, parameter :: jiopt = 4
  !> Parameter used for sorting internally in sparse1 matrix module 
  integer, parameter :: sort = 1
  !> Parameter used for sorting internally in sparse1 matrix module
  integer, parameter :: nosort = 2
  !> Threshold for elements to neglect
  real(realk) :: zero = 1E-9_realk   
  !> Threshold for elements to neglect, squared
  real(realk) :: zero2 = 1E-18_realk
  !> Overload: The '*' sign may be used to scale a type(smat1)
  interface operator(*)
     module procedure sp_scale_smat1
  end interface

contains

!> See mat_inquire_cutoff in mat-operations.f90
  subroutine mat_sparse1_inquire_cutoff(cutoff)
    implicit none
    real(realk), intent(out) :: cutoff

    cutoff=zero

  end subroutine mat_sparse1_inquire_cutoff

!> See mat_zero_cutoff in mat-operations.f90
  subroutine mat_sparse1_zero_cutoff(cutoff)
    implicit none
    real(realk), intent(in) :: cutoff

    zero=cutoff
    zero2 = zero**2

  end subroutine mat_sparse1_zero_cutoff
  
!> \brief See mat_init in mat-operations.f90
  subroutine mat_sparse1_init(a,nrow,ncol)
     implicit none
     TYPE(Matrix) :: a 
     integer, intent(in) :: nrow, ncol
     NULLIFY(a%selm1)
     NULLIFY(a%idata)
     ALLOCATE(a%idata(nidata))
     a%nrow = nrow
     a%ncol = ncol
     a%idata(nelm) = 0
     a%idata(symm) = usymm
     a%idata(pack) = unzipd

  end subroutine mat_sparse1_init

!> \brief See mat_free in mat-operations.f90
  subroutine mat_sparse1_free(a)
     implicit none
     type(matrix) :: a
     integer :: nsize

     if (.not.ASSOCIATED(a%selm1)) then
       if (a%idata(Nelm) /= 0) then
         STOP 'Error in mat_sparse1_free - memory previously released'
       endif
     else
       nsize = SIZE(a%selm1)*mem_realsize
       call sp1_stat_deallocated_memory(nsize)
       DEALLOCATE(a%selm1)
     endif
     DEALLOCATE(A%idata)
     NULLIFY(a%selm1,a%idata)
  end subroutine mat_sparse1_free

!> \brief Puts a full matrix on sparse matrix form and multiplies it with a scalar alpha
!> \param afull The full matrix
!> \param alpha The scalar
!> \param a afull on sparse matrix form
  subroutine mat_sparse1_set_from_full(afull,alpha,a)
     implicit none
     real(realk), INTENT(IN) :: afull(*)
     real(realk), intent(in) :: alpha
     TYPE(Matrix)            :: a  !output
     integer                 :: i, j,num_elm

     a%idata(nelm) = 0
     do i = 1,a%nrow*a%ncol
       if (ABS(afull(i)) > zero) a%idata(nelm) = a%idata(nelm) + 1
     enddo
     call sp_allocate(a,a%idata(nelm))
     num_elm = 0
     do j = 1,a%ncol
       do i = 1,a%nrow
         if (ABS(afull(a%nrow*(j-1)+i)) > zero) then
           num_elm = num_elm + 1
           a%selm1(num_elm)%i = i
           a%selm1(num_elm)%j = j
           a%selm1(num_elm)%val = alpha*afull(a%nrow*(j-1)+i)
         endif
       enddo
     enddo
     if (num_elm /= a%idata(nelm)) then
       print*,num_elm,a%idata(nelm)
       STOP 'programming error in mat_sparse1_set_from_full'
     endif
      
  end subroutine mat_sparse1_set_from_full

!> \brief Puts a sparse matrix on full matrix form and multiplies it with a scalar.
!> \param a The sparse matrix
!> \param alpha The scalar
!> \param afull a on full matrix form
  subroutine mat_sparse1_to_full(a, alpha, afull)
     implicit none
     TYPE(Matrix), intent(in):: a
     real(realk), intent(in) :: alpha
     real(realk), intent(out):: afull(*)  !output
     integer                 :: i,ndim

     if (.not.ASSOCIATED(a%selm1)) STOP 'a in mat_sparse_to_full non-existant'
     if (a%idata(pack) /= unzipd) then
       STOP 'implement me in mat_sparse1_to_full'
     endif
     ndim = a%nrow*a%ncol
     afull(1:ndim) = 0.0E0_realk
     do i = 1,a%idata(nelm)
       !index in afull: nrow(j-1)+i
       afull(a%nrow*(a%selm1(i)%j - 1) + a%selm1(i)%i) = a%selm1(i)%val*alpha
     enddo

  end subroutine mat_sparse1_to_full

!> \brief Prints the full matrix section a(irow1:i_rown,j_col1:j_coln), 4 columns at the time
!> \param A The matrix to be written
!> \param i_row1 Start printing at this row
!> \param i_rown End printing at this row
!> \param j_col1 Start printing at this column
!> \param j_coln End printing at this column
!> \param lu The logical unit number of the file
  subroutine mat_sparse1_print(A, i_row1, i_rown, j_col1, j_coln, lu)
     implicit none
     TYPE(Matrix),intent(in) :: A
     integer, intent(in)     :: i_row1, i_rown, j_col1, j_coln, lu
     type(matrix)            :: ax
     real(realk)             :: val4(4)
     integer                 :: ip(a%nrow+1), ncols, icol, i, j, k, &
                              & rest

     if(.not.ASSOCIATED(A%selm1)) STOP 'A in mat_sparse1_print non-existant'
     call mat_sparse1_init(ax,a%nrow,a%ncol)
     call mat_sparse1_assign(ax,a)
     call sp_sort(ijopt,ax)
     call sp_pointers(iopt,nosort,ax,ip)
     !number of columns to print = j_coln - j_col1 + 1
     !number of rows to print = i_rown - i_row1 + 1
     ncols = j_coln - j_col1 + 1
     icol = j_col1 - 1
     do i = 1,ncols/4
       WRITE(lu,*)
       WRITE(lu,'(9x,4("   Column",i6))') icol+1,icol+2,icol+3,icol+4
       do j = i_row1,i_rown
         !find the nonzero elements beloning to this row
         do k = 1,4
           if (IP(j)+1 <= a%idata(nelm)) then
             if(Ax%selm1(IP(j)+1)%j==icol+k.and.Ax%selm1(IP(j)+1)%i==j)then
               IP(j)=IP(j)+1
               val4(k)=Ax%selm1(IP(j))%val
             else
               val4(k)=0.0000000000
             endif
           else
             val4(k)=0.0000000000
           endif
         enddo
         WRITE(lu,'(i7,2x,4F15.8)') j,(val4(k),k=1,4)
       enddo
       icol = icol + 4
     enddo
     rest=MOD(ncols,4)
     if (rest>0) then
       WRITE(lu,*)
       WRITE(lu,'(9x,4("   Column",i6))') icol+1,icol+2,icol+3,icol+4
       do j=i_row1,i_rown
         !find the nonzero elements belonging to this row
         do k=1,rest
           if (IP(j)+1 <= a%idata(nelm)) then
             if (Ax%selm1(IP(j)+1)%j==icol+k.and.Ax%selm1(IP(j)+1)%i==j)then
                IP(j)=IP(j)+1
                val4(k)=Ax%selm1(IP(j))%val  
             else
                val4(k)=0.0000000000
             endif
           else
             val4(k)=0.00000000000
           endif
         enddo
         WRITE(lu,'(i7,2x,4F15.8)') j,(val4(k),k=1,rest)
       enddo
     endif
     call mat_sparse1_free(ax)
  end subroutine mat_sparse1_print

!> \brief See mat_trans in mat-operations.f90
!> \param A The matrix leaves the routine unchanged. 
!> \param B The transpose of A
!>
!> Nothing is assumed about the order of the matrix elements.
!> All variables in A should be defined upon entry. \n
!>
  subroutine mat_sparse1_trans(A,B)
    implicit none
    type(matrix), intent(in) :: A
    type(matrix)             :: B !output
    integer :: i

    if(.not. ASSOCIATED(A%selm1)) STOP 'input in mat_sparse1_trans non-existant'
    b%idata(nelm) = a%idata(nelm)
    b%idata(symm) = a%idata(symm)
    b%idata(pack) = a%idata(pack)    
    call sp_allocate(b,a%idata(nelm))

    do i=1,A%idata(Nelm)
      b%selm1(i)%i = A%selm1(i)%j
      b%selm1(i)%j = A%selm1(i)%i
      b%selm1(i)%val = A%selm1(i)%val
    enddo
  end subroutine mat_sparse1_trans

!> \brief Make an exact copy of a matrix.
!> \param b The matrix you want a copy of
!> \param a The copy of b
  subroutine mat_sparse1_assign(a,b)
    implicit none
    TYPE(Matrix), INTENT(INOUT) :: a
    TYPE(Matrix), INTENT(IN)    :: b
    integer                     :: i

    if(.not.ASSOCIATED(b%selm1)) STOP 'input in mat_sparse1_assign non-existant'
    call sp_allocate(a,b%idata(nelm))
    A%idata(Nelm)=b%idata(Nelm)
    A%idata(pack)=b%idata(pack)
    A%idata(symm)=b%idata(symm)
    do i=1,A%idata(Nelm)
      A%selm1(i)=b%selm1(i)
    enddo
  end subroutine mat_sparse1_assign

!> \brief Make a MPI broadcast of matrix.
!> \param a The matrix
!> \param slave logical if slave process
!> \param master integer of master process
  subroutine mat_sparse1_mpicopy(a,slave,master)
! NOT TESTET
    use lsmpi_type
    implicit none
    TYPE(Matrix), TARGET, INTENT(INOUT) :: a
    logical                     :: slave
    integer                     :: i,nrow,ncol
    integer(kind=ls_mpik)       :: master

    nrow=a%nrow
    ncol=a%ncol
    call LS_MPI_BUFFER(nrow,Master)
    call LS_MPI_BUFFER(ncol,Master)
    IF(SLAVE)THEN
       NULLIFY(a%iaux, a%raux)
       nullify(A%elms)
       nullify(A%elmsb)
       a%init_self_ptr => a
       a%init_magic_tag = mat_init_magic_value
       call mat_sparse1_init(a,nrow,ncol)
    ENDIF
    IF(.NOT.SLAVE)i = a%idata(nelm)
    call LS_MPI_BUFFER(i,Master)
    IF(SLAVE)call sp_allocate(a,i)
    call LS_MPI_BUFFER(a%idata(Nelm),Master)
    call LS_MPI_BUFFER(a%idata(pack),Master)
    call LS_MPI_BUFFER(a%idata(symm),Master)
    do i=1,A%idata(Nelm)
      call LS_MPI_BUFFER(A%selm1(i)%i,Master)
      call LS_MPI_BUFFER(A%selm1(i)%j,Master)
      call LS_MPI_BUFFER(A%selm1(i)%val,Master)
    enddo

  end subroutine mat_sparse1_mpicopy

!> \brief Make a copy of a matrix multiplied by a scalar
!> \param alpha The scalar
!> \param a The matrix you want to multiply by alpha and copy to b
!> \param b The copy of a*alpha
  subroutine mat_sparse1_copy(alpha,a,b)
     implicit none
     REAL(REALK),  INTENT(IN)    :: alpha
     TYPE(Matrix), INTENT(IN)    :: a
     TYPE(Matrix), INTENT(INOUT) :: b
     integer                     :: i
     
     if (.not.ASSOCIATED(a%selm1)) STOP 'input in mat_sparse1_copy non-existant'
     call sp_allocate(b,a%idata(nelm))
     b%idata(Nelm)=a%idata(Nelm)
     b%idata(pack)=a%idata(pack)
     b%idata(symm)=a%idata(symm)
     do i=1,A%idata(Nelm)
       b%selm1(i)=alpha*a%selm1(i)
     enddo
  end subroutine mat_sparse1_copy

!> \brief Makes the trace of a sparse matrix.
!> \param A The sparse matrix
!> \return The trace of A
  function mat_sparse1_Tr(A)
     implicit none
     TYPE(Matrix), intent(IN) :: A
     REAL(realk) :: mat_sparse1_tr
     integer :: i

     if (.not.ASSOCIATED(a%selm1)) STOP 'input in mat_sparse1_Tr non-existant'
     mat_SParse1_Tr = 0.0E0_realk
     do i=1,A%idata(Nelm)
       if (A%selm1(i)%i==A%selm1(i)%j) then
         mat_SParse1_Tr = mat_SParse1_Tr + A%selm1(i)%val
       endif
     enddo
  end function mat_sparse1_Tr

!> \brief See mat_TrAB in mat-operations.f90
!>
!> sum_ij A_ij B_ji = sum_ij A_ij B^T_ij
!>
  function mat_sparse1_TrAB(a,b)
     implicit none
     TYPE(Matrix), intent(IN) :: a,b
     REAL(realk) :: mat_sparse1_trAB
     type(matrix) :: abx,bT
     integer :: i,index

     if(.not.ASSOCIATED(a%selm1).or..not.ASSOCIATED(b%selm1)) then
        STOP 'input in mat_sparse1_TrAB non-existant'
     endif
     if (a%idata(pack) /= unzipd .or. b%idata(pack) /= unzipd) then
       STOP 'mat_sparse1_trAB is not implemented for packed stuff'
     endif
     call mat_sparse1_init(bT,b%ncol,b%nrow)
     call mat_sparse1_trans(b,bT)
     call mat_sparse1_init(abx,a%nrow,a%ncol)
     abx%idata(nelm) = a%idata(nelm) + b%idata(nelm)
     call sp_allocate(abx,abx%idata(nelm))
     do i = 1,a%idata(nelm)
       abx%selm1(i) = a%selm1(i)
     enddo
     do i = a%idata(nelm)+1,abx%idata(nelm)
       abx%selm1(i) = bT%selm1(i-a%idata(nelm))
     enddo 
     call mat_sparse1_free(bT)
     call sp_sort(ijopt,abx)
     mat_sparse1_TrAB = 0.0E0_realk
     index = 1
     do
       if (index >= abx%idata(nelm)) exit
       if (abx%selm1(index)%i == abx%selm1(index+1)%i .and.&
          &abx%selm1(index)%j == abx%selm1(index+1)%j) then
         mat_sparse1_TrAB = mat_sparse1_TrAB + &
                     & abx%selm1(index)%val*abx%selm1(index+1)%val
         index = index + 2
       else
         index = index + 1
       endif
     enddo  
     call mat_sparse1_free(abx)
  end function mat_sparse1_TrAB

!> \brief Computes C = alpha*AB + beta*C
!> 
!> Maximum allocation during routine:~2*NCX~=2*(2*C%Nelm)>A%Nelm+B%Nelm \n
!> Allocated on exit: (allocation on entry)+C%Nelm \n
!> Nothing is assumed about the order of the matrix elements.
!> 
!> \param A First sparse matrix to be multiplied.
!> \param B Second sparse matrix to be multiplied,
!> \param transa 'n'/'N': The matrix A is multiplied as is. 't'/'T': The transposed of A is multiplied.
!> \param transb 'n'/'N': The matrix B is multiplied as is. 't'/'T': The transposed of B is multiplied.
!> \param alpha The multiplied matrices AB is scaled with alpha
!> \param beta If beta /= 0, the input matrix C is scaled by beta and added to the matrix alpha*AB
!> \param C The product of A and B scaled by alpha and added with beta*C(in)
  subroutine mat_sparse1_mul(A,B,transa, transb,alpha,beta,C)
     !c = alpha*ab + beta*c
     !transa = 'T'/'t' - transposed, 'N'/'n' - normal
     implicit none
     TYPE(Matrix), intent(IN) :: A, B
     character, intent(in)    :: transa, transb
     REAL(realk), INTENT(IN)  :: alpha, beta
     TYPE(Matrix), intent(inout) :: C   
     type(matrix)             :: xa, xb, xxc,temp
     type(smat1), allocatable:: xc(:),xcx(:)
     real(realk)              :: zerosqrt,zeroN
     integer                  :: i,j,k,nindex,klo,khi,&
                               & nc,nxc,ncon,jpc(c%ncol+1),flops,ngrand,ipc(c%nrow+1),NiC(c%nrow)
     integer, allocatable     :: ipa(:),ipb(:),NiB(:)
     logical                  :: aT,bT,Asmall

     if (.not.ASSOCIATED(A%selm1)) STOP 'input mat1 in MAT_sparse1_MUL non-existant'
     if (.not.ASSOCIATED(B%selm1)) STOP 'input mat2 in MAT_sparse1_MUL non-existant'
     if (ABS(beta) > zero2 .and. .not. ASSOCIATED(c%selm1))Then
        STOP 'input mat3 in MAT_sparse1_mul non-existant'
     ENDIF
     if (a%idata(pack)/=unzipd .or. b%idata(pack)/=unzipd) then
       STOP 'mat_sparse1_mul not implemented for packed matrices'
     endif
!** Initializing
     zeroN = zero/a%nrow
     zerosqrt = sqrt(zeroN)
     ngrand = 0
     flops = 0
     NiC = 0
     if (transa == 'n'.or.transa == 'N') then
       aT = .false.
       ALLOCATE(ipa(a%nrow+1))
     elseif (transa == 't'.or.transa == 'T') then
       aT = .true.
       ALLOCATE(ipa(a%ncol+1))
     else
       STOP ' wrong input for transa in mat_sparse1_mul'
     endif
     if (transb == 'n' .or. transb == 'N') then
       bT = .false.
       ALLOCATE(ipb(b%nrow+1),NiB(b%nrow))
     elseif (transb == 't' .or. transb == 'T') then
       bT = .true.
       ALLOCATE(ipb(b%ncol+1),NiB(b%ncol))
     else
       STOP 'wrong input for transb in mat_sparse1_mul'
     endif
     call mat_sparse1_init(xxc,a%nrow,b%ncol)
     call mat_sparse1_init(temp,a%nrow,b%ncol)
!*********************************************
!If beta /= 0, sort C in the row-index
!*********************************************
     if (ABS(beta) > zero2) then
       call mat_sparse1_assign(xxc,c)
       call sp_pointers(iopt,sort,xxc,ipc)
       do i = 1,c%nrow
         NiC(i) = ipc(i+1) - ipc(i)
       enddo
     endif
!***********************************************
!Sort A in the row-index if .not. aT
!Sort A in the column-index if aT
!***********************************************
     call mat_sparse1_init(xa,a%nrow,a%ncol)
     call mat_sparse1_assign(xa,a)
     if (.not.aT) then
       call SP_POINTERS(iopt,sort,XA,IPA)
     else
       call SP_POINTERS(jopt,sort,XA,IPA)
     endif
!*****************************
!Sort B in the row-index if .not. bT
!Sort B in the column-index if bT
!*****************************
     call mat_sparse1_init(xb,b%nrow,b%ncol)
     call mat_sparse1_assign(xb,b)
     if (.not.bT) then
       call SP_POINTERS(iopt,sort,XB,IPB)
       nindex = b%nrow
     else
       call sp_pointers(jopt,sort,xb,ipb)
       nindex = b%ncol
     endif
!***********************************************
!If alfa /= 1, multiply A or B by alfa
!***********************************************
     if (alpha > 1.0E0_realk + zero .or. alpha < 1E0_realk - zero) then
       if (a%idata(nelm) < b%idata(nelm)) then
         do i = 1,a%idata(nelm)
           xa%selm1(i)%val = xa%selm1(i)%val * alpha
         enddo
       else
         do i = 1,b%idata(nelm)
           xb%selm1(i)%val = xb%selm1(i)%val * alpha
         enddo
       endif
     endif
!******************************************************************
!Make array with information about the number of B-elements with
!certain i-index if .not. bT or j-index if bT
!******************************************************************
     NiB = 0
     do i = 1,nindex
       NiB(i) = IPB(i+1) - IPB(i)
     enddo
!*******************************************************************
!Find the max number of contributions to one C-row 
!*******************************************************************
     NXC = 0
     if (.not.aT) then
       do i = 1,a%nrow
         NC = NiC(i)  !NiC(i) = 0 if beta = 0
         do j = IPA(i)+1,IPA(i+1)
           NC = NC + NiB(XA%selm1(j)%j)
         enddo
         NXC = MAX(NXC,NC)
       enddo
     else
       do i = 1,a%ncol
         NC = NiC(i)  !NiC(i) = 0 if beta = 0
         do j = IPA(i)+1,IPA(i+1)
           NC = NC + NiB(XA%selm1(j)%i)
         enddo
         NXC = MAX(NXC,NC)
       enddo
     endif
!******************************************************************
!Allocate space for the contributions. If there are more than
!60EE6 contributions, it is not possible to store all on a 32-bit
!machine.
!******************************************************************
     C%idata(nelm) = 0
     if (NXC /= 0) then
       ALLOCATE(xc(nxc),xcx(nxc))
  !******************************************************************
  !Compute the contributions to C. Row by row. For efficient use of
  !memory, the computation is split in pieces
  !******************************************************************
       Ncon = 0
       JPC = 0
  !** Allocate room for C
       call sp_allocate(c,nxc)
       if (.not.aT) then
         nindex = a%nrow
       else
         nindex = a%ncol
       endif
       do i = 1,nindex
  !Make certain row in C
         if (IPA(i+1)-IPA(i) > 0) then
           do j = IPA(i)+1,IPA(i+1)
             if (.not.aT) then
               klo = ipb(xa%selm1(j)%j)+1
               khi = ipb(xa%selm1(j)%j+1)
             else
               klo = ipb(xa%selm1(j)%i)+1
               khi = ipb(xa%selm1(j)%i+1)
             endif
             if (ABS(XA%selm1(j)%val) < zerosqrt) then
               Asmall = .true.
             else
               Asmall = .false.
             endif
             do k = klo,khi
               if (Asmall .and. ABS(XB%selm1(k)%val) < zerosqrt) then
                 cycle
               endif
               XC(Ncon+1)%val = XA%selm1(j)%val*XB%selm1(k)%val  
               flops = flops + 1
               if (ABS(XC(Ncon+1)%val) > zeroN) then
                 ncon = ncon + 1
                 if (.not.aT) then
                   XC(Ncon)%i   = XA%selm1(j)%i
                 else
                   XC(Ncon)%i   = XA%selm1(j)%j
                 endif
                 if (.not.bT) then
                   XC(Ncon)%j   = XB%selm1(k)%j
                 else
                   XC(Ncon)%j   = XB%selm1(k)%i
                 endif
  !** On the go we also make array of pointers to beg. of elements with certain
  !** j-index in Cx
                 JPC(XC(Ncon)%j+1) = JPC(XC(Ncon)%j+1) + 1  
               endif             
             enddo
           enddo
           if (ABS(beta) > zero2) then
             do j = ipc(i)+1,ipc(i+1)
               ncon = ncon + 1
               xc(ncon) = xxc%selm1(j)
               JPC(xc(ncon)%j+1) = JPC(xc(ncon)%j+1) + 1
             enddo
           endif
           if (Ncon == 0) cycle
  !*************************************************************
  !Evaluate row by row
  !**************************************************************
  ! print*,"ncon",Ncon
  !** Sort XC in the j-index. All i-indexes are the same
           do j = 2,C%Ncol+1
             JPC(j) = JPC(j-1) + JPC(j)
           enddo
           do j = 1,Ncon
             JPC(XC(j)%j) = JPC(XC(j)%j) + 1
             XCX(JPC(XC(j)%j)) = XC(j)
           enddo
  !** Allocate room for C?
           NC = 1
           do j = 2,Ncon
             !nc = number of new contributions to c
             if (XCX(j)%j /= XCX(j-1)%j) then
               nc = nc + 1
             endif 
           enddo
           if (c%idata(nelm) + nc > SIZE(c%selm1)) then
             call mat_sparse1_assign(temp,c)
             call sp_allocate(c,c%idata(nelm)+(nindex-i+1)*nc)
             !if mat_sparse1_assign(c,temp) was called, it could resize c again
             do j = 1,temp%idata(nelm)
               c%selm1(j) = temp%selm1(j)
             enddo
             c%idata(nelm) = temp%idata(nelm)
             ngrand = ngrand + 1
           endif
  !** Sum and merge elements with similar i- and j-index
  !** Make C-array larger if necessary
           C%idata(Nelm) = C%idata(Nelm) + 1
           C%selm1(C%idata(Nelm)) = XCX(1)
           do j = 2,Ncon
             if(XCX(j)%j==XCX(j-1)%j)then
               C%selm1(C%idata(Nelm))%val = C%selm1(C%idata(Nelm))%val + XCX(j)%val
             else
               if (ABS(C%selm1(C%idata(Nelm))%val) > zero) then
                 C%idata(Nelm) = C%idata(Nelm) + 1
               endif
               C%selm1(C%idata(Nelm)) = XCX(j)
             endif
           enddo
           if (ABS(C%selm1(C%idata(Nelm))%val) <= zero) then
             !overwrite
             C%idata(Nelm) = C%idata(Nelm) - 1
           endif
           Ncon = 0
           JPC = 0
         endif
       enddo
       DEALLOCATE(xc,xcx)
     endif
     if (C%idata(Nelm) == 0) then
       C%idata(Nelm) = 1
       call sp_allocate(c,1)
       c%selm1(1)%i = 2
       c%selm1(1)%j = 1
       c%selm1(1)%val = 0.0E0_realk
     endif
     c%idata(symm) = usymm
     c%idata(pack) = unzipd
     call mat_sparse1_free(xa)
     call mat_sparse1_free(xb)
     if (ASSOCIATED(xxc%selm1)) call mat_sparse1_free(xxc)
     if (ASSOCIATED(temp%selm1)) call mat_sparse1_free(temp)
     if (mat_info) then
       WRITE(mat_lu,*) 'mat_sparse1_mul:'
       WRITE(mat_LU,*) 'flops : ',flops,' c-resizing: ',ngrand, ' c%nelm/size(c): ', &
                       &REAL(c%idata(nelm))/REAL(SIZE(c%selm1))
       WRITE(mat_lu,*) 'a%nelm/size: ',rEAL(a%idata(nelm))/(a%nrow*a%ncol),' b%nelm/size: ',&
                       &real(b%idata(nelm))/(b%nrow*b%ncol),' c%nelm/size: ',&
                       &real(c%idata(nelm))/(c%nrow*c%ncol)
     endif
!** TODO - Make C-array smaller?
  end subroutine mat_sparse1_mul                                                                  

!> \brief Computes C = alpha*A + beta*B
!>
!> Nothing is assumed about the order of the matrix elements. \n
!> Maximum allocation during routine: 2*(A%Nelm+B%Nelm)+MAX(Nrow,Ncol) \n
!> Allocated on exit: (allocation on entry)+C%Nelm 
!>
!> \param A First matrix to be added
!> \param B Second matrix to be added
!> \param alpha A is scaled by alpha before addition
!> \param beta B is scaled by beta before addition
!> \param C The sum of A and B
  subroutine mat_sparse1_add(alpha,A,beta,B,C)
     implicit none
     TYPE(Matrix), intent(IN) :: A, B
     REAL(realk), INTENT(IN)  :: alpha, beta
     TYPE(Matrix)             :: C   !output
     type(matrix)             :: xc
     integer                  :: i,maxelm, nc
     if (.not.ASSOCIATED(A%selm1)) STOP 'A in MAT_sparse1_add non-existant'
     if (.not.ASSOCIATED(B%selm1)) STOP 'B in MAT_sparse1_add non-existant'
     if (a%idata(pack) /= b%idata(pack)) THEN
        STOP 'mat_sparse1_add not implemented for different packed input matrices'
     endif
!   
!**************************************************************************
!When adding the matrices A and B, the maximum number of elements in the
!resultmatrix will be maxelm:
!**************************************************************************
     maxelm=A%idata(Nelm)+B%idata(Nelm)
     if (maxelm == 0) then
       !C is a zero matrix
       call mat_sparse1_zero(C)
       return
     endif
     call mat_sparse1_init(xc,a%nrow,a%ncol)
     xc%idata(nelm) = maxelm
     call sp_allocate(xc,maxelm)
!** Make XC-sparsematrix with all the elements from A and B
     do i=1,A%idata(Nelm)
       XC%selm1(i) = alpha*A%selm1(i)
     enddo
     do i=A%idata(Nelm)+1,maxelm
       XC%selm1(i) = beta*B%selm1(i-A%idata(Nelm))
     enddo
!** Sort the elements in both indexes
     call SP_SORT(jiopt,XC)
!** How large should C be?
     nc = maxelm
     i = 2
     do 
       if (i > maxelm) exit
       if (XC%selm1(i)%i == XC%selm1(i-1)%i .and. XC%selm1(i)%j == XC%selm1(i-1)%j) then
!FIXME: How much would it cost to check the size of the elements here too?
         nc = nc - 1
         i = i+2
       else
         i = i+1 
       endif
     enddo
!************************
!** If C exists but possibly is too short, it is de- and re-allocated
!** If C does not exist, it is allocated
     call sp_allocate(c,nc)
!************************
!** Sum and merge elements with similar i- and j-index
     nc = 1
     C%selm1(1)=XC%selm1(1)
     do i=2,maxelm
       if(XC%selm1(i)%i==XC%selm1(i-1)%i.and.XC%selm1(i)%j==XC%selm1(i-1)%j)then
         C%selm1(nc)%val=C%selm1(nc)%val+XC%selm1(i)%val
       else
         if (ABS(c%selm1(nc)%val) > zero) then
           nc=nc+1
         endif
         C%selm1(nc)=XC%selm1(i)
       endif
     enddo
     if (ABS(c%selm1(nc)%val) < zero) nc = nc - 1
!** Make output array as small as possible?
     call mat_sparse1_free(xc)
!** Assign values to the variables of C
     C%idata(Nelm) = nc
     if (A%idata(Symm) == B%idata(Symm)) then
       C%idata(Symm) = A%idata(Symm)
     else
       C%idata(Symm) = usymm
     endif
     C%idata(pack) = A%idata(pack)
     if (mat_info) then
       WRITE(mat_LU,*) 'mat_sparse1_add: c%nelm/size(c) = ',REAL(c%idata(nelm))/SIZE(c%selm1)
       WRITE(mat_lu,*) 'a%nelm/size: ',REAL(a%idata(nelm))/(a%nrow*a%ncol),' b%nelm/size: ',&
                       &REAL(b%idata(nelm))/(b%nrow*b%ncol),' c%nelm/size: ',&
                       &REAL(c%idata(nelm))/(c%nrow*c%ncol)
     endif

  end subroutine mat_sparse1_add

!> \brief Y = alpha*X + Y
!> \param alpha scalar 
!> \param X Matrix to be multiplied with alpha and added to Y
!> \param Y Matrix to be added to alpha*X. On output, the sum of alpha*X and Y
  subroutine mat_sparse1_daxpy(alpha,x,y)
     implicit none
     real(realk),intent(in)       :: alpha
     TYPE(Matrix), intent(IN)     :: X
     TYPE(Matrix), intent(INOUT)  :: Y
     type(matrix)                 :: y2
     integer                      :: i,maxelm,ny

     if (.not.ASSOCIATED(X%selm1)) STOP 'X in mat_SParse1_DAXPY non-existant'
     if (.not.ASSOCIATED(Y%selm1)) STOP 'Y in mat_SParse1_DAXPY non-existant'
     if (x%idata(pack) /= y%idata(pack)) then
        STOP 'different packed input mats in mat_sparse1_daxpy not implemented'
     endif
!**************************************************************************
!When adding the matrices X and Y, the maximum number of elements in the
!resultmatrix will be maxelm:
!**************************************************************************
     maxelm=X%idata(Nelm) + Y%idata(Nelm)
     call mat_sparse1_init(y2,y%nrow,y%ncol)
     call sp_allocate(y2,maxelm)
     Y2%idata(Nelm) = maxelm
!** Make Y2-sparsematrix with all the elements from aX and Y
     do i=1,X%idata(Nelm)
       Y2%selm1(i)=alpha*X%selm1(i)
     enddo
     do i=X%idata(Nelm)+1,maxelm
       Y2%selm1(i)=Y%selm1(i-X%idata(Nelm))
     enddo
!** Sort the elements in both indexes
     call SP_SORT(jiopt,Y2)
!** How large should Y be?
     ny = maxelm
     i = 2
     do 
       if (i > maxelm) exit
       if (Y2%selm1(i)%i == Y2%selm1(i-1)%i .and. Y2%selm1(i)%j == Y2%selm1(i-1)%j) then
!FIXME: How much would it cost to check the size of the elements here too?
         ny = ny - 1
         i = i+2
       else
         i = i+1 
       endif
     enddo
!************************
!** If Y is too short, it is de- and re-allocated
     call sp_deallocate(y)
     call sp_allocate(y,ny)
!************************
!** Sum and merge elements with similar i- and j-index
     ny=1
     Y%selm1(1)=Y2%selm1(1)
     do i=2,maxelm
       if((Y2%selm1(i)%i==Y2%selm1(i-1)%i).and.(Y2%selm1(i)%j==Y2%selm1(i-1)%j))then
         Y%selm1(ny)%val=Y%selm1(ny)%val+Y2%selm1(i)%val
       else
         if (ABS(y%selm1(ny)%val) > zero) then
           ny=ny+1
         endif
         Y%selm1(ny)=Y2%selm1(i)
       endif
     enddo
     if (ABS(y%selm1(ny)%val) < zero) ny = ny - 1
     y%idata(nelm) = ny
!** Make output array as small as possible ?
    call mat_sparse1_free(y2)
!** Assign values to the variables of C
    if (x%idata(symm) /= y%idata(symm)) then
      y%idata(symm) = usymm
    endif

  end subroutine mat_sparse1_daxpy

!> \brief Makes the dot-product of two sparse matrices
!> \param A The first sparse matrix
!> \param B The second sparse matrix
!> \return The dotproduct of the matrices
  function mat_sparse1_dotproduct(A,B)
     implicit none
     TYPE(Matrix), intent(IN) :: A,B
     REAL(realk) :: mat_sparse1_dotproduct
     integer  :: i,totelm
     type(matrix) :: ab

     if (.not.ASSOCIATED(A%selm1)) STOP 'A in mat_SParse1_DOTPRODUCT non-existant'
     if (.not.ASSOCIATED(B%selm1)) STOP 'B in mat_SParse1_DOTPRODUCT non-existant'
     totelm=A%idata(Nelm)+B%idata(Nelm)
     call mat_sparse1_init(ab,a%nrow,a%ncol)
     call sp_allocate(ab,totelm)
     AB%idata(Nelm) = totelm
 
 !** Make sparsematrix with all the elements from A and B
     do i=1,A%idata(Nelm)
       AB%selm1(i) = A%selm1(i)
     enddo
     do i=A%idata(Nelm)+1,totelm
       AB%selm1(i) = B%selm1(i-A%idata(Nelm))
     enddo
 !** Sort the elements in both indexes
     call SP_SORT(ijopt,AB)
 !** Elements with similar indexes should be multiplied and summed up 
     mat_sparse1_dotproduct = 0E0_realk
     i = 2
     do 
       if (i > totelm) exit
       if (AB%selm1(i)%i == AB%selm1(i-1)%i.and.&
          &AB%selm1(i)%j == AB%selm1(i-1)%j) then
         mat_sparse1_dotproduct = mat_sparse1_dotproduct + &
                                & AB%selm1(i)%val*AB%selm1(i-1)%val
         i = i + 2
       else
         i = i + 1
       endif
     enddo
     call mat_sparse1_free(ab)

  end function mat_sparse1_dotproduct

!> \brief Makes the square norm of a sparse matrix
!> \param A The sparse matrix
!> \return The square norm of A
  function mat_sparse1_sqnorm2(A)
     implicit none
     TYPE(Matrix), intent(IN) :: A
     REAL(realk) :: mat_sparse1_sqnorm2
     integer :: i

     if(.not.ASSOCIATED(A%selm1)) STOP 'A in mat_SParse1_sqNORM2 non-existant'
     mat_sparse1_sqnorm2 = 0.0E0_realk
     do i=1,A%idata(Nelm)
       mat_sparse1_sqnorm2 = mat_sparse1_sqnorm2 + A%selm1(i)%val**2
     enddo

  end function mat_sparse1_sqnorm2

!> \brief Find the largest element of a sparse matrix.
!> \param A The sparse matrix
!> \return The largest element of A
   REAL(realk) function mat_sparse1_max_elm(A)
     implicit none
     TYPE(Matrix), intent(IN) :: A
     integer :: i

     if(.not.ASSOCIATED(A%selm1)) STOP 'A in mat_SParse1_max_elm non-existant'
     mat_sparse1_max_elm = MAXVAL(A%selm1(:)%val)
  end function mat_sparse1_max_elm

!> \brief  Makes the square norm of the off-diagonal elements of a sparse matrix.
!> \param A The sparse matrix
!> \return The square norm of the off-diagonal elements of A
  function mat_sparse1_outdia_sqnorm2(A)
     implicit none
     TYPE(Matrix), intent(IN) :: A
     REAL(realk) :: mat_sparse1_outdia_sqnorm2
     integer :: i

     if(.not.ASSOCIATED(A%selm1)) STOP 'A in mat_SParse1_outdia_sqNORM2 non-existant'
     mat_sparse1_outdia_sqnorm2 = 0.0E0_realk
     do i=1,A%idata(Nelm)
       if (A%selm1(i)%i /= A%selm1(i)%j) then
         mat_sparse1_outdia_sqnorm2 = mat_sparse1_outdia_sqnorm2 + A%selm1(i)%val**2
       endif
     enddo

  end function mat_sparse1_outdia_sqnorm2

!> \brief See mat_diag_f in mat-operations.f90
  subroutine mat_sparse1_diag_f(F,S,eival,Cmo)
     !solves FC = SCe 
     implicit none
     TYPE(Matrix), intent(IN) :: F,S
     type(matrix)             :: Cmo  !output
     real(realk), intent(OUT) :: eival(:)
     integer                  :: ndim
     real(realk), allocatable :: wrk(:),tmp(:),Cmo_full(:)
     real(realk) :: dummy
    
     ndim = s%nrow
     ALLOCATE(tmp(Ndim*Ndim),Cmo_full(Ndim*Ndim))
  ! diagonalizaton carried out in full basis
     call mat_sparse1_to_full(S, 1E0_realk, tmp)
     call mat_sparse1_to_full(F,1E0_realk,Cmo_full)
  !
     call my_DSYGV(ndim,Cmo_full,tmp,eival,"SPARSE1_GET_DENS_SOL")
  !
     call mat_sparse1_set_from_full(Cmo_full,1E0_realk,Cmo)
  !
     DEALLOCATE(tmp,Cmo_full)
  end subroutine mat_sparse1_diag_f

  !subroutine mat_sparse1_d_from_f(X,S,Dnew,eival,Cmo)
  !   use matrix_operations_dense
     !solves FdampC = SCe and construct new density (Dnew) from the 
     !eigenvectors C
  !   implicit none
  !   TYPE(Matrix), intent(IN) :: X,S
  !   type(matrix)             :: Dnew,Cmo  !output
  !   type(matrix)             :: eival     !intent(out)
  !   integer                  :: ndim
  !   real(realk), allocatable :: wrk(:),tmp(:),Cmo_full(:),Dnew_full(:)
  !   real(realk), allocatable :: eig(:)
  !   real(realk) :: dummy
  !   logical :: getd, getdv
     
  !   ndim = s%nrow  
  !   ALLOCATE(wrk((8+2)*Ndim),tmp(Ndim*Ndim),Cmo_full(Ndim*Ndim), &
  !          & Dnew_full(Ndim*Ndim),eig(Ndim))
  !   call mat_sparse1_to_full(S, 1E0_realk, tmp)
  !   call mat_sparse1_to_full(X,1E0_realk,Cmo_full)
  !   call my_DSYGV(ndim,Cmo_full,tmp,eig,wrk,size(wrk),"SPARSE1_GET_DENS_SOL")
     !** create new density from eivec
  !   getd = .true.; getdv = .false.  !calculate DCAO, not DVAO
  !   call fckden(getd,getdv,Dnew_full,dummy,Cmo_full,dummy,scratch_area,scratch_size)
     !because it is a closed shell case, the diagonal = 2 and we would like
     !it to be 1 for idempotency reasons
  !   call mat_sparse1_set_from_full(Cmo_full,1E0_realk,Cmo)
  !   call mat_sparse1_set_from_full(Dnew_full,0.5E0_realk,Dnew)
  !   call mat_sparse1_set_from_full(eig,1.0E0_realk,eival)  
  !   DEALLOCATE(wrk,tmp,Cmo_full,Dnew_full,eig)

  ! end subroutine mat_sparse1_d_from_f

!> \brief See mat_dE_dmu in mat-operations.f90
  FUNCTION mat_sparse1_dE_dmu(Fnew,Fdamp,SDS,Cmo,nocc)
     !Computes dE_SCF/dmu where mu is the damping in damped roothan
     !Find dE/dmu = sum_nu,I dE/dC_nu,I dC_nu,I/dmu
     implicit none
     type(matrix),intent(in) :: Fnew,Fdamp,SDS,Cmo
     integer, intent(in)     :: nocc
     real(realk)             :: mat_sparse1_dE_dmu
     type(matrix)            :: CTFC, tSDS, Cmo_occ, X, X2
     real(realk)             :: tFdamp_dia(Fnew%nrow)
     integer                 :: ndim,i,index

     ndim = Fnew%nrow
   ![dE/dC]ij = 4[FC]ij  !FIXME: is this correct for DFT?????
     call mat_sparse1_init(CTFC,ndim,nocc)
     call mat_sparse1_init(X,ndim,nocc)
     call mat_sparse1_init(cmo_occ,ndim,nocc)
     call mat_sparse1_section(Cmo,1,ndim,1,nocc,Cmo_occ)
 write(mat_lu,*) 'nocc',nocc
 write(mat_lu,*) 'Cmo:'
call mat_sparse1_print(Cmo, 1, ndim, 1, ndim, mat_lu)
write(mat_lu,*) 'Cmo_occ:'
call mat_sparse1_print(Cmo_occ, 1, ndim, 1, nocc, mat_lu)
     call mat_sparse1_mul(Fnew,Cmo_occ,'n','n',1E0_realk,0E0_realk,X)
     call mat_sparse1_mul(Cmo,X,'t','n',4.0E0_realk,0E0_realk,CTFC)

   ![dC/dmu]ij = -[dC/dmu]ji => [dC/dmu]ij = 0 for i=j
   ![dC/dmu]ij = [C^TSDSC]ij/(C^TFdampC_ii - C^TFdampC_jj))
     call mat_sparse1_init(tSDS,ndim,nocc)
     call mat_sparse1_mul(SDS,Cmo_occ,'n','n',1E0_realk,0E0_realk,X)
     call mat_sparse1_mul(Cmo,X,'t','n',1E0_realk,0E0_realk,tSDS)
     call mat_sparse1_free(Cmo_occ)
     !construct the C^TFdampC diagonal - maybe some day a mat_get_dia routine 
     !will be made and some of this would be replaced with a call for that
     !tFdamp_ii = sum_km (Fdamp_mk*C_mi*C_ki) = sum_k C_ki*[Fdamp C]_ki
     call mat_sparse1_free(X)
     call mat_sparse1_init(X,ndim,ndim)
     call mat_sparse1_mul(Fdamp,Cmo,'n','n',1E0_realk,0E0_realk,X)
     !collect [Fdamp C] and C in same matrix and sort (standard procedure)
     call mat_sparse1_init(X2,ndim,ndim)
     x2%idata(nelm) = x%idata(nelm) + Cmo%idata(nelm)
     call sp_allocate(x2,x2%idata(nelm))
     do i = 1,x%idata(nelm)
       x2%selm1(i) = x%selm1(i)
     enddo
     do i = x%idata(nelm)+1,x2%idata(nelm)
       x2%selm1(i) = Cmo%selm1(i-x%idata(nelm))
     enddo 
     call mat_sparse1_free(X)
     call sp_sort(ijopt,x2)
     tFdamp_dia = 0.0E0_realk
     index = 1
     do
       if (index >= x2%idata(nelm)) exit
       if (x2%selm1(index)%i == x2%selm1(index+1)%i .and.&
          &x2%selm1(index)%j == x2%selm1(index+1)%j) then
         tFdamp_dia(x2%selm1(index)%j) = tFdamp_dia(x2%selm1(index)%j) + & 
                     & x2%selm1(index)%val*x2%selm1(index+1)%val
         index = index + 2
       else
         index = index + 1
       endif
     enddo
     call mat_sparse1_free(X2)    
     do i = 1,tSDS%idata(nelm)
       if (tSDS%selm1(i)%i == tSDS%selm1(i)%j) then
         tSDS%selm1(i)%val = 0E0_realk
       else
         tSDS%selm1(i)%val =                tSDS%selm1(i)%val/&
                             &(tFdamp_dia(tSDS%selm1(i)%i) - tFdamp_dia(tSDS%selm1(i)%j))
       endif
     enddo
   !dE/dmu = sum_nu,I dE/dC_nu,I dC_nu,I/dmu
     mat_sparse1_dE_dmu = mat_sparse1_dotproduct(CTFC,tSDS)
     call mat_sparse1_free(CTFC)
     call mat_sparse1_free(tSDS)
  END FUNCTION mat_sparse1_dE_dmu

!> \brief See mat_column_norm in mat-operations.f90
  function mat_sparse1_column_norm(Mat,ncol,from_row,to_row)
     implicit none
     type(Matrix), intent(in) :: Mat
     integer, intent(in) :: ncol, from_row, to_row
     real(realk) :: mat_sparse1_column_norm
     type(Matrix) :: xmat
     integer :: i,jp(Mat%ncol +1)

     if (.not.ASSOCIATED(Mat%selm1)) STOP 'input non-existant in mat_sparse1_column_norm'
     call mat_sparse1_init(xmat,mat%nrow,mat%ncol)
     call mat_sparse1_assign(xmat,mat)
     call sp_pointers(jopt,sort,xMat,jp)
     mat_sparse1_column_norm = 0E0_realk
     !run over all elements with ncol as column index
     do i = jp(ncol)+1,jp(ncol+1)
       if (xMat%selm1(i)%i >= from_row .and. xMat%selm1(i)%i <= to_row) then
         mat_sparse1_column_norm = mat_sparse1_column_norm + (xMat%selm1(i)%val)**2
       endif
     enddo
     call mat_sparse1_free(xmat)

  end function mat_sparse1_column_norm

!> \brief See mat_section in mat-operations.f90
  subroutine mat_sparse1_section(A,from_row,to_row,from_col,to_col,Asec)
     implicit none
     type(Matrix), intent(in) :: A
     integer, intent(in) :: from_row, to_row, from_col, to_col
     type(Matrix), intent(inout) :: Asec  !output
     integer :: i, index, irow, icol

     if (.not.ASSOCIATED(A%selm1)) STOP 'input non-existant in mat_sparse1_section'
     Asec%nrow = to_row - from_row + 1
     Asec%ncol = to_col - from_col + 1
     index = 0 
     do i = 1,A%idata(nelm)
       if (A%selm1(i)%i >= from_row .and. A%selm1(i)%i <= to_row &
          & .and. A%selm1(i)%j >= from_col .and. A%selm1(i)%j <= to_col) then 
         index = index + 1
       endif
     enddo
     call sp_allocate(Asec,index)
     Asec%idata(nelm) = index
     index = 0
     do i = 1,A%idata(nelm)
       if (A%selm1(i)%i >= from_row .and. A%selm1(i)%i <= to_row &
          & .and. A%selm1(i)%j >= from_col .and. A%selm1(i)%j <= to_col) then 
         index = index + 1
         Asec%selm1(index)%i = A%selm1(i)%i - from_row + 1
         Asec%selm1(index)%j = A%selm1(i)%j - from_col + 1
         Asec%selm1(index)%val = A%selm1(i)%val
       endif
     enddo
  end subroutine mat_sparse1_section
 
!> \brief See mat_mo_precond in mat-operations.f90
  subroutine mat_sparse1_mo_precond(nocc,omega,Eorb_final,X_MO)
    implicit none
    integer, intent(in) :: nocc
    real(realk), intent(in) :: omega, Eorb_final(:)
    type(Matrix), intent(inout) :: X_MO
    real(realk) :: dia
    integer :: i

    do i = 1,X_MO%idata(nelm)
      if (X_MO%selm1(i)%i > nocc .and. X_MO%selm1(i)%j <= nocc) then
        !lower left block (ai block)
        !dia = E[2]_dia + omega*S[2]_dia
        dia = 2.0E0_realk*(Eorb_final(X_MO%selm1(i)%i) - Eorb_final(X_MO%selm1(i)%j)) &
            & + omega * 2E0_realk
        X_MO%selm1(i)%val = X_MO%selm1(i)%val/dia
      elseif (X_MO%selm1(i)%j > nocc .and. X_MO%selm1(i)%i <= nocc) then
        !upper right block (ia block)
        !dia = E[2]_dia + omega*S[2]_dia
        dia = 2.0E0_realk*(Eorb_final(X_MO%selm1(i)%j) - Eorb_final(X_MO%selm1(i)%i)) &
            & - omega * 2E0_realk
        X_MO%selm1(i)%val = X_MO%selm1(i)%val/dia
      endif
    enddo

  end subroutine mat_sparse1_mo_precond
   
!> \brief See mat_identity in mat-operations.f90
  subroutine mat_sparse1_identity(I)
    implicit none
    type(Matrix), intent(inout) :: I
    integer :: n
    call sp_allocate(I,I%nrow)
    do n = 1,I%nrow
      I%selm1(n)%i = n
      I%selm1(n)%j = n
      I%selm1(n)%val = 1.0E0_realk
    enddo
    I%idata(nelm) = I%nrow

  end subroutine mat_sparse1_identity

!> \brief See mat_create_block in mat-operations.f90
  subroutine mat_sparse1_create_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
    implicit none
    integer, intent(in) :: fullrow,fullcol,insertrow,insertcol
    real(Realk), intent(in) :: fullmat(fullrow,fullcol)
    type(Matrix), intent(inout) :: A
    type(Matrix) :: ax
    integer :: i, j, k, newsize, extra_elms ! = number of extra elms to be allocated
    logical :: found

    !A. Very simple and unefficient implementation:
    !do i = insertrow, insertrow+fullrow-1
    !   do j = insertcol, insertcol+fullcol-1
    !      !write(mat_lu,*) 'i,j,val =', i,j,fullmat(i-insertrow+1,j-insertcol+1)
    !      call mat_sparse1_create_elm(i,j,fullmat(i-insertrow+1,j-insertcol+1),A)
    !   enddo
    !enddo

    !B. A little more efficient implementation - could still be improved by
    !sorting:
    !Check how much extra space must be allocated:
    extra_elms = 0
    do i = insertrow, insertrow+fullrow-1
       do j = insertcol, insertcol+fullcol-1
          found = .FALSE.
          do k = 1,A%idata(nelm)
             if ((A%selm1(k)%i == i) .and. (A%selm1(k)%j == j)) then
                !nothing - elm already exists
                found = .TRUE.
             endif
          enddo
          if(.NOT.found)THEN
             if (abs(fullmat(i-insertrow+1,j-insertcol+1)) > zero) then
                extra_elms = extra_elms + 1
             endif
          endif
       enddo
    enddo

    !If extra space is needed, allocate it:
    if (extra_elms /= 0) then
       newsize = A%idata(nelm) + extra_elms
       call mat_sparse1_init(ax,a%nrow,a%ncol)
       IF (ASSOCIATED(A%selm1)) THEN 
         call mat_sparse1_assign(ax,A)
         call sp_deallocate(A)
         call sp_allocate(A,newsize)
         call mat_sparse1_assign(A,ax)
         call mat_sparse1_free(ax)
       ELSE
         call sp_allocate(A,newsize)
       ENDIF
    endif

    !We have now allocated needed space, insert elements:
    do i = insertrow, insertrow+fullrow-1
       do j = insertcol, insertcol+fullcol-1
          if (abs(fullmat(i-insertrow+1,j-insertcol+1)) > zero) then
             call mat_sparse1_create_elm(i,j,fullmat(i-insertrow+1,j-insertcol+1),A)
          endif
       enddo
    enddo

  end subroutine mat_sparse1_create_block

!> \brief See mat_retrieve_block in mat-operations.f90
  subroutine mat_sparse1_retrieve_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
    implicit none
    integer, intent(in) :: fullrow,fullcol,insertrow,insertcol
    real(Realk), intent(out) :: fullmat(fullrow,fullcol)
    type(Matrix), intent(in) :: A
    type(Matrix) :: ax
    integer :: i, j, k, newsize, extra_elms ! = number of extra elms to be allocated
    logical :: found

    extra_elms = 0
    do i = insertrow, insertrow+fullrow-1
       do j = insertcol, insertcol+fullcol-1
          found = .FALSE.
          do k = 1,A%idata(nelm)
             if ((A%selm1(k)%i == i) .and. (A%selm1(k)%j == j)) then
                !elm exists
                fullmat(i-insertrow+1,j-insertcol+1) = A%selm1(k)%val
                found = .TRUE.
             endif
          enddo
          if(.NOT.found)THEN
                fullmat(i-insertrow+1,j-insertcol+1) = 0E0_realk
          endif
       enddo
    enddo

  end subroutine mat_sparse1_retrieve_block

!> \brief See mat_add_block in mat-operations.f90
  subroutine mat_sparse1_add_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
    implicit none
    integer, intent(in) :: fullrow,fullcol,insertrow,insertcol
    real(Realk), intent(inout) :: fullmat(fullrow,fullcol)
    type(Matrix), intent(inout) :: A
    type(Matrix) :: ax
    integer :: i, j, k, newsize, extra_elms ! = number of extra elms to be allocated
    logical :: found

    !A. Very simple and unefficient implementation:
    !do i = insertrow, insertrow+fullrow-1
    !   do j = insertcol, insertcol+fullcol-1
    !      !write(mat_lu,*) 'i,j,val =', i,j,fullmat(i-insertrow+1,j-insertcol+1)
    !      call mat_sparse1_create_elm(i,j,fullmat(i-insertrow+1,j-insertcol+1),A)
    !   enddo
    !enddo

    !B. A little more efficient implementation - could still be improved by
    !sorting:
    !Check how much extra space must be allocated:
    extra_elms = 0
    do i = insertrow, insertrow+fullrow-1
       do j = insertcol, insertcol+fullcol-1
          found = .FALSE.
          do k = 1,A%idata(nelm)
             if ((A%selm1(k)%i == i) .and. (A%selm1(k)%j == j)) then
                !elm already exists so we add to block in order to overwrite later
                fullmat(i-insertrow+1,j-insertcol+1) = fullmat(i-insertrow+1,j-insertcol+1)+A%selm1(k)%val
                found = .TRUE.
             endif
          enddo
          if(.NOT.found)THEN
             if (abs(fullmat(i-insertrow+1,j-insertcol+1)) > zero) then
                extra_elms = extra_elms + 1
             endif
          endif
       enddo
    enddo

    !If extra space is needed, allocate it:
    if (extra_elms /= 0) then
       newsize = A%idata(nelm) + extra_elms
       call mat_sparse1_init(ax,a%nrow,a%ncol)
       IF (ASSOCIATED(A%selm1)) THEN
         call mat_sparse1_assign(ax,A)
         call sp_deallocate(A)
         call sp_allocate(A,newsize)
         call mat_sparse1_assign(A,ax)
         call mat_sparse1_free(ax)
       ELSE
         call sp_allocate(A,newsize)
       ENDIF
    endif

    !We have now allocated needed space, insert elements:
    do i = insertrow, insertrow+fullrow-1
       do j = insertcol, insertcol+fullcol-1
          if (abs(fullmat(i-insertrow+1,j-insertcol+1)) > zero) then
             call mat_sparse1_create_elm(i,j,fullmat(i-insertrow+1,j-insertcol+1),A)
          endif
       enddo
    enddo

  end subroutine mat_sparse1_add_block

!> \brief See mat_create_elm in mat-operations.f90
  subroutine mat_sparse1_create_elm(i,j,val,A)
    implicit none
    integer, intent(in) :: i,j
    real(realk), intent(in) :: val
    type(Matrix), intent(inout) :: A
    type(Matrix) :: ax
    integer :: n
    logical :: found
    found = .false.
    do n = 1,A%idata(nelm)
      if (A%selm1(n)%i == i .and. A%selm1(n)%j == j) then
        A%selm1(n)%val = val
        found = .true.
        exit
      endif
    enddo

    if (.not.found) then
      if (A%idata(nelm) == 0 ) then
        call sp_allocate(A,1)
      endif

      if (SIZE(A%selm1) < A%idata(nelm)+1) then
        call mat_sparse1_init(ax,a%nrow,a%ncol)
        call mat_sparse1_assign(ax,A)
        call sp_deallocate(A)
        call sp_allocate(A,A%idata(nelm)+1)
        call mat_sparse1_assign(A,ax)
        call mat_sparse1_free(ax)
      endif
      A%idata(nelm) = A%idata(nelm) + 1
      A%selm1(A%idata(nelm))%i = i
      A%selm1(A%idata(nelm))%j = j
      A%selm1(A%idata(nelm))%val = val
    endif
  end subroutine mat_sparse1_create_elm

!> \brief Scale a sparse matrix.
!> \param alpha The scaling parameter which is multiplied with every element in A
!> \param A The sparse matrix
  subroutine mat_sparse1_scal(alpha,A)
    implicit none
    real(realk), intent(in) :: alpha
    type(Matrix), intent(inout) :: A
    integer :: i

    if (.not.ASSOCIATED(A%selm1)) STOP 'input in mat_sparse1_SCAL non-existant'
    do i = 1,A%idata(Nelm)
      A%selm1(i)%val = alpha*A%selm1(i)%val
    enddo

  end subroutine mat_sparse1_scal

!> \brief See mat_scal_dia in mat-operations.f90
  subroutine mat_sparse1_scal_dia(alpha,A)
    implicit none
    real(realk), intent(in) :: alpha
    type(matrix), intent(inout) :: A
    integer :: i
    if (.not.ASSOCIATED(A%selm1)) STOP 'input in mat_sparse1_SCAL_DIA non-existant'
    do i=1,A%idata(Nelm)
      if (A%selm1(i)%i == A%selm1(i)%j) then
        A%selm1(i)%val = alpha*A%selm1(i)%val
      endif
    enddo

  end subroutine mat_sparse1_scal_dia

!> \brief See mat_zero in mat-operations.f90
!>
!>    The elements of A should all be zero, in sparse terms that
!>    gives a matrix with no elements.
!>    Since a matrix doesn't really exist before it is associated
!>    with some memory, we make one zero element.
!>
  subroutine mat_sparse1_zero(A)
    implicit none
    type(Matrix), intent(inout) :: A
    A%idata(nelm) = 1
    if (.not.ASSOCIATED(A%selm1)) then
      call sp_allocate(A,1)
    endif
    if (A%nrow > 1) then
      A%selm1(1)%i = 2
      A%selm1(1)%j = 1
    elseif (A%ncol > 1) then
      A%selm1(1)%i = 1
      A%selm1(1)%j = 2
    else
      A%selm1(1)%i = 1
      A%selm1(1)%j = 1
    endif
   A%selm1(1)%val = 0.0E0_realk

  end subroutine mat_sparse1_zero

!> \brief Set the lower or upper triangle to zero.
!> \param part If 'LT'/'lt': The lower triangle+diagonal is set to 0. If 'UT'/'ut': The upper triangle is set to 0, diagonal is kept
!> \param A The sparse matrix. On output: without the upper or lower triangle
  subroutine mat_sparse1_zerohalf(part,A)
    implicit none
    character(len=2),intent(in) :: part
    type(matrix), intent(inout) :: A
    type(Matrix) :: A_x
    integer :: i,nlower, nupper

    if (.not.ASSOCIATED(A%selm1)) STOP 'input in mat_sparse1_zerohalf non-existant'
    nlower = 0; nupper = 0
    do i = 1,A%idata(Nelm)
      if (A%selm1(i)%i < A%selm1(i)%j) then
        nupper = nupper + 1
      else
        nlower = nlower + 1
      endif
    enddo
    call mat_sparse1_init(A_x,A%nrow,A%ncol)
    call mat_sparse1_assign(A_x,A)
    if (part == 'LT' .or. part == 'lt') then
      call sp_ALLOCATE(A,nupper)
      A%idata(Nelm) = 0
      do i = 1,A_x%idata(Nelm)
        if (A_x%selm1(i)%i < A_x%selm1(i)%j) then
          A%idata(Nelm) = A%idata(Nelm) + 1
          A%selm1(A%idata(Nelm)) = A_x%selm1(i)
        endif
      enddo 
    elseif (part == 'UT' .or. part == 'ut') then
      call sp_ALLOCATE(A,nlower)
      A%idata(Nelm) = 0
      do i = 1,A_x%idata(Nelm)
        if (A_x%selm1(i)%i >= A_x%selm1(i)%j) then
          A%idata(Nelm) = A%idata(Nelm) + 1
          A%selm1(A%idata(Nelm)) = A_x%selm1(i)
        endif
      enddo
    else
      STOP ' unknown option given to mat_sparse1_ZEROHALF'
    endif
    call mat_sparse1_free(A_x)
  end subroutine mat_sparse1_zerohalf

!> \brief Writes matrix A to file.
!> \param iunit The logical unit number for the file
!> \param A The sparse matrix that should be written to disk
  subroutine mat_sparse1_write_to_disk(iunit,A)
    implicit none
    integer, intent(in) :: iunit
    type(Matrix), intent(in) :: A
    character(len=20) :: filename

    if (.not.ASSOCIATED(A%selm1)) STOP 'A in mat_sparse1_WRITE_TO_DISK non-existant'

    WRITE(iunit) A%Nrow, A%Ncol, A%idata(nelm)
    WRITE(iunit) A%selm1
  end subroutine mat_sparse1_write_to_disk


!> \brief Reads a matrix from disk.
!> \param iunit The logical unit number for the file
!> \param A The matrix read from disk
  subroutine mat_sparse1_read_from_disk(iunit,A)
    implicit none
    integer, intent(in) :: iunit
    type(Matrix), intent(inout) :: A
    integer :: i

    READ(iunit) A%Nrow, A%Ncol, A%idata(Nelm)
    call sp_allocate(A,A%idata(nelm))
    READ(iunit) A%selm1(1:A%idata(nelm)) 

  end subroutine mat_sparse1_read_from_disk

!> \brief See mat_vec_to_mat in mat-operations.f90
  subroutine mat_sparse1_vec_to_mat(symmetry, VEC, MAT)
     implicit none
     integer                  :: n, k, nelmmat, xint
     real(realk)              :: xreal
     character, intent(in)    :: symmetry
     type(matrix), intent(in) :: VEC
     type(Matrix)             :: MAT !Intent(out)
     LOGICAL                  :: DIAG

   if (.not.ASSOCIATED(VEC%selm1)) STOP 'input in mat_sparse1_vec_to_mat non-existant'
   if (symmetry == 's' .OR. symmetry == 'S') then
      call sp_allocate(MAT,VEC%idata(nelm)*2)
   else if (symmetry == 'a' .OR. symmetry == 'A') then
      call sp_allocate(MAT,VEC%idata(nelm)*2 + VEC%nrow) !Diagonal of matrix not included in vector size
   endif

!      do n = 1, dimension_of_matrix     
!         do m = 1, n
!            i = n*(n-1)/2 + m
!            MAT(n,m) = VEC(i)
!            MAT(m,n) = -VEC(i)		 ! +VEC(i) for symmetric matrices
!         enddo
!      enddo

   !If no elements in VEC, MAT is a zero matrix:
   if (VEC%idata(nelm) == 0)  then
      call mat_sparse1_zero(MAT)
   endif

   nelmmat = 0
   if (symmetry == 's' .OR. symmetry == 'S') then
      do n = 1, VEC%idata(nelm)
      !Before putting the n'th element in MAT, we must determine it's position:
         xreal = 2*VEC%selm1(n)%i                     !Find largest x so that x**2 =< 2k, k is vector index
         xreal = SQRT(xreal)
         xint = INT(xreal)
         if (xint*(xint+1)/2 + 1 > VEC%selm1(n)%i) then                    !If x*(x+1)/2 + 1 > n
            nelmmat = nelmmat + 1
            MAT%selm1(nelmmat)%i = xint                                    !Mat(i,j) = (x, n - x*(x-1)/2)
            MAT%selm1(nelmmat)%j = VEC%selm1(n)%i - xint*(xint-1)/2
            !Now, put element in MAT:
            MAT%selm1(nelmmat)%val = VEC%selm1(n)%val
            !If element is not diagonal, also do transposed element:
            if (MAT%selm1(nelmmat)%i /= MAT%selm1(nelmmat)%j) then
               nelmmat = nelmmat + 1
               MAT%selm1(nelmmat)%j = xint                               !Mat(j,i) = (n - x*(x-1)/2, x)
               MAT%selm1(nelmmat)%i = VEC%selm1(n)%i - xint*(xint-1)/2
               !Now, put transposed element in MAT:
               MAT%selm1(nelmmat)%val = VEC%selm1(n)%val
            endif
         else
            nelmmat = nelmmat + 1
            MAT%selm1(nelmmat)%i = xint + 1                                !Mat(i,j) = (x + 1, n - x*(x-1)/2)
            MAT%selm1(nelmmat)%j = VEC%selm1(n)%i - (xint+1)*((xint+1)-1)/2
            !Now, put element in MAT:
            MAT%selm1(nelmmat)%val = VEC%selm1(n)%val
            !If element is not diagonal, also do transposed element:
            if (MAT%selm1(nelmmat)%i /= MAT%selm1(nelmmat)%j) then
               nelmmat = nelmmat + 1
               MAT%selm1(nelmmat)%j = xint + 1                           !Mat(j,i) = (n - x*(x-1)/2,x + 1)
               MAT%selm1(nelmmat)%i = VEC%selm1(n)%i - (xint+1)*((xint+1)-1)/2
               !Now, put transposed element in MAT:
               MAT%selm1(nelmmat)%val = VEC%selm1(n)%val
            endif
         endif
      enddo
   endif
   if (symmetry == 'a' .OR. symmetry == 'A') then
      do n = 1, VEC%idata(nelm)
      !Before putting the n'th element in MAT, we must determine it's position:
         xreal = 2*VEC%selm1(n)%i                     !Find largest x so that x**2 =< 2k, k is vector index
         xreal = SQRT(xreal)
         xint = INT(xreal)
         if (xint*(xint+1)/2 + 1 > VEC%selm1(n)%i) then                    !If x*(x+1)/2 + 1 > n
            nelmmat = nelmmat + 1
            MAT%selm1(nelmmat)%i = xint + 1 !RETTET!!!!                    !Mat(i,j) = (x+1, n - x*(x-1)/2)
            MAT%selm1(nelmmat)%j = VEC%selm1(n)%i - xint*(xint-1)/2
            !Now, put element in MAT:
            MAT%selm1(nelmmat)%val = VEC%selm1(n)%val
            !Element cannot be diagonal, also do transposed element:
            nelmmat = nelmmat + 1
            MAT%selm1(nelmmat)%j = xint + 1                                !Mat(j,i) = (n - x*(x-1)/2, x+1)
            MAT%selm1(nelmmat)%i = VEC%selm1(n)%i - xint*(xint-1)/2
            !Now, put transposed element in MAT:
            MAT%selm1(nelmmat)%val = -VEC%selm1(n)%val
         else
            nelmmat = nelmmat + 1
            MAT%selm1(nelmmat)%i = xint + 2  !RETTET!!!!                   !Mat(i,j) = (x+2, n - x*(x-1)/2)
            MAT%selm1(nelmmat)%j = VEC%selm1(n)%i - (xint+1)*((xint+1)-1)/2
            !Now, put element in MAT:
            MAT%selm1(nelmmat)%val = VEC%selm1(n)%val
            !Element cannot be diagonal, also do transposed element:
            nelmmat = nelmmat + 1
            MAT%selm1(nelmmat)%j = xint + 2                                 !Mat(j,i) = (n - x*(x-1)/2, x+2)
            MAT%selm1(nelmmat)%i = VEC%selm1(n)%i - (xint+1)*((xint+1)-1)/2
            !Now, put transposed element in MAT:
            MAT%selm1(nelmmat)%val = -VEC%selm1(n)%val
         endif
      enddo
   endif

  MAT%idata(nelm) = nelmmat
  end subroutine mat_sparse1_vec_to_mat

!> \brief See mat_to_vec in mat-operations.f90
  subroutine mat_sparse1_mat_to_vec(symmetry, MAT, VEC)
     implicit none
     character, intent(in)    :: symmetry
     integer                  :: n, k, nelmvec
     type(matrix), intent(in) :: MAT
     type(Matrix)             :: VEC !Intent(out)

   if (.not.ASSOCIATED(MAT%selm1)) STOP 'input in mat_sparse1_mat_to_vec non-existant'
   !#elements in vec cannot exceede nelm(mat)/nelm(mat) + ndim
   call sp_allocate(VEC,MAT%idata(nelm)/2+MAT%nrow) 

   !do n = 1, dimension_of_matrix     !Algorithm concept. 
   !   do m = 1, n
   !      i = n*(n-1)/2 + m
   !      !VEC(i) = MAT(n,m)
   !   enddo
   !enddo

   nelmvec = 0

   if (symmetry == 's' .OR. symmetry == 'S') then
      do n = 1, MAT%idata(nelm)
         if (MAT%selm1(n)%i >= MAT%selm1(n)%j) then
            nelmvec = nelmvec + 1
            k = MAT%selm1(n)%i*(MAT%selm1(n)%i - 1)/2 + MAT%selm1(n)%j !k = n*(n-1)/2 + m
            VEC%selm1(nelmvec)%i = k
            VEC%selm1(nelmvec)%j = 1
            VEC%selm1(nelmvec)%val = MAT%selm1(n)%val                   !VEC(k) = MAT(n,m)
         endif
      enddo
   endif

   if (symmetry == 'a' .OR. symmetry == 'A') then
      do n = 1, MAT%idata(nelm)
         if (MAT%selm1(n)%i > MAT%selm1(n)%j) then
            nelmvec = nelmvec + 1
            k = MAT%selm1(n)%i*(MAT%selm1(n)%i - 1)/2 + MAT%selm1(n)%j  & !k = n*(n-1)/2 + m
              & - MAT%selm1(n)%i + 1
            VEC%selm1(nelmvec)%i = k
            VEC%selm1(nelmvec)%j = 1
            VEC%selm1(nelmvec)%val = MAT%selm1(n)%val                   !VEC(k) = MAT(n,m)
         endif
      enddo
   endif

   VEC%idata(nelm) = nelmvec
  end subroutine mat_sparse1_mat_to_vec

!> \brief See mat_report_sparsity in mat-operations.f90
  subroutine mat_sparse1_report_sparsity(A,sparsity)
  implicit none
  type(Matrix) :: A
  real(realk)  :: nnz,sparsity
  integer      :: i

       nnz=0E0_realk
       do i=1, A%idata(nelm)
         nnz = nnz+1E0_realk
       enddo

       sparsity = nnz/(A%ncol*A%nrow)
  end subroutine mat_sparse1_report_sparsity

!*****************************************************************************
!Routines private to this module

!>  \brief Local allocation routine for this module.
!>
!> Check that matrix A is allocated and has room for nsize elements. If not,
!> deallocate it and reallocate it with proper size. It is more efficient to
!> deallocate and allocate the whole matrix instead of just single elements, 
!> because this ensures consecutive memory, which is crucial for efficiency.
!> Allocated memory is counted using sp1_stat_allocated_memory.
!>
!> \param A The sparse matrix
!> \param nsize The requested size of the matrix
  subroutine sp_allocate(A,nsize)
     implicit none
     type(matrix) :: A
     integer, intent(in) :: nsize
     integer :: sz

     if (ASSOCIATED(a%selm1)) then
       if (SIZE(a%selm1) < nsize .or. waste_frac*nsize > SIZE(a%selm1) ) then
         call sp_deallocate(a)
         ALLOCATE(a%selm1(nsize))
         sz = nsize*mem_realsize
         call sp1_stat_allocated_memory(sz)
       endif
     else
       ALLOCATE(a%selm1(nsize))
       sz = nsize*mem_realsize
       call sp1_stat_allocated_memory(sz)
     endif
  end subroutine sp_allocate

!> \brief Local deallocation routine for this module.
!>
!> Deallocate a matrix A and count the deallocated memory using sp1_stat_deallocated_memory.
!>
!> \param A The sparse matrix
  subroutine sp_deallocate(A)
     implicit none
     type(matrix) :: A
     integer :: nsize
 
     nsize = SIZE(a%selm1)*mem_realsize
     call sp1_stat_deallocated_memory(nsize)
     DEALLOCATE(a%selm1)
     NULLIFY(a%selm1)

  end subroutine sp_deallocate

!> \brief See stat_allocated_memory in mat-operations.f90
  subroutine sp1_stat_allocated_memory(size)
     implicit none
     integer, intent(in) :: size
     integer(kind=long) :: nsize
     nsize = size
     call mem_allocated_mem_type_matrix(nsize)
  
  end subroutine sp1_stat_allocated_memory

!> \brief See stat_deallocated_memory in mat-operations.f90
  subroutine sp1_stat_deallocated_memory(size)
     implicit none
     integer, intent(in) :: size
     integer(kind=long) :: nsize
     nsize = size  
     call mem_deallocated_mem_type_matrix(nsize)

  end subroutine sp1_stat_deallocated_memory

!> \brief Local scaling routine used for type(smat1).
!> \param alpha The scaling parameter
!> \param A The type(smat1) matrix to be scaled
  function sp_scale_smat1(alpha,A)
     type(smat1), intent(in) :: A
     real(realk), intent(in) :: alpha
     type(smat1)             :: sp_scale_smat1

     sp_scale_smat1%i = a%i
     sp_scale_smat1%j = a%j
     sp_scale_smat1%val = alpha*a%val
  end function sp_scale_smat1

!> \brief Local sorting routine for this module.
!>
!> Maximum allocation during routine:  A%idata(Nelm)+MAX(Nrow,Ncol) \n
!> Allocated on exit: same as allocation on entry. \n
!> Nothing is assumed about the order of the matrix elements. \n
!> All variables in A should be defined upon entry.
!>
!> \param opt If opt = iopt, A is only sorted in the i-index. If opt = jopt, A is only sorted in the j-index. If opt = ijopt, A is sorted in both indices, but with i as super-index. If opt = jiopt, A is sorted in both indices, but with j as super-index.
!> \param A The matrix is sorted according to the opt-value
  subroutine sp_sort(opt,A)
    implicit none
    integer, intent(in) :: opt
    type(matrix), intent(inout) :: A
    type(smat1), allocatable :: xa(:)
    if(.not.ASSOCIATED(A%selm1)) STOP 'A in SP_SORT non-existant'
    ALLOCATE(xa(a%idata(nelm)))
    select case(opt)
    case(iopt)
        call sp_isort(a,xa)
    case(jopt)
        call sp_jsort(a,xa)
    case(ijopt)
        call sp_jsort(a,xa)
        call sp_isort(a,xa)
    case(jiopt)
        call sp_isort(a,xa)
        call sp_jsort(a,xa)
    case default
      STOP 'sp_sort called with undefined sorting-rule'
    end select
    DEALLOCATE(xa)
  end subroutine sp_sort

!> \brief Local sorting routine for this module - sort a sparse1 type(matrix) in i-index.
!> \param A The matrix to be sorted
!> \param xa Scratch array
  subroutine sp_isort(A,xa)
    implicit none
    type(matrix), intent(inout) :: A
    type(smat1), intent(inout) :: xa(:) !scratch array
!    integer :: IP(A%Nrow+1),i
    integer              :: i
    integer, allocatable :: IP(:)

    allocate(IP(A%Nrow+1))

!*******************************
!Set up pointers
!*******************************
    IP=0
    !write(mat_lu,*) 'isort, nelm =',a%idata(nelm)
    do i=1,A%idata(Nelm)
      !write(mat_lu,'(2i10,f10.3," elm ",i4)') a%selm1(i)%i, a%selm1(i)%j, a%selm1(i)%val, i
      if (A%selm1(i)%i>0) then
        IP(A%selm1(i)%i+1)=IP(A%selm1(i)%i+1)+1
      else
        print*,'i',a%selm1(i)%i
        print*,'j',a%selm1(i)%j
        print*,'val',a%selm1(i)%val
        STOP 'wrong i-indexes in sp_isort'
      endif
    enddo
    do i=2,A%Nrow+1
      IP(i)=IP(i-1)+IP(i)   !points to [entry-1] of first element with
    enddo                   !certain i-index in A
!***************
    do i=1,A%idata(Nelm)
      XA(i)=A%selm1(i)
    enddo
    do i=1,A%idata(Nelm)
      IP(XA(i)%i)=IP(XA(i)%i)+1
      A%selm1(IP(XA(i)%i))=XA(i)
    enddo

    deallocate(IP)
  end subroutine sp_isort

!> \brief Local sorting routine for this module - sort a sparse1 type(matrix) in j-index.
!> \param A The matrix to be sorted
!> \param xa Scratch array
  subroutine sp_jsort(A,xa)
    implicit none
    type(matrix), intent(inout) :: A
    type(smat1), intent(inout) :: xa(:) !scratch array
    integer :: JP(A%Ncol+1),i

!*******************************
!Set up pointers
!*******************************
    JP=0
    do i=1,A%idata(Nelm)
      if (A%selm1(i)%j>0) then
        JP(A%selm1(i)%j+1)=JP(A%selm1(i)%j+1)+1
      else
        STOP 'wrong j-indexes in sp_jsort'
      endif
    enddo
    do i=2,A%Ncol+1
      JP(i)=JP(i-1)+JP(i)   !points to [entry-1] of first element with
    enddo                   !certain i-index in A
!***************
    do i=1,A%idata(Nelm)
      XA(i)=A%selm1(i)
    enddo
    do i=1,A%idata(Nelm)
      JP(XA(i)%j)=JP(XA(i)%j)+1
      A%selm1(JP(XA(i)%j))=XA(i)
    enddo
  end subroutine sp_jsort

!> \brief Local routine for this module. Set up pointers for sparse matrices.
  subroutine sp_pointers(opt,optsort,A,ip)
    implicit none
    integer, intent(in) :: opt,optsort
    type(matrix), intent(inout) :: A
    integer, intent(out) :: ip(:)
    if(.not.ASSOCIATED(A%selm1)) STOP 'A in SP_pointers non-existant'
    select case(opt)
    case(iopt)
        call sp_ipointers(optsort,a,ip)
    case(jopt)
        call sp_jpointers(optsort,a,ip)
    case default
      STOP 'sp_pointers called with undefined option'
    end select
  end subroutine sp_pointers

!> \brief Local routine for this module. Set up pointers for sparse matrices.
  subroutine sp_ipointers(optsort,a,ip)
    implicit none
    integer, intent(in) :: optsort
    type(matrix), intent(inout) :: a
    integer, intent(out) :: ip(a%nrow+1)
    type(smat1),allocatable :: xa(:) !scratch array
    integer :: i
!*******************************
!Set up pointers
!*******************************
    IP=0
    do i=1,A%idata(Nelm)
      if (A%selm1(i)%i>0) then
        IP(A%selm1(i)%i+1)=IP(A%selm1(i)%i+1)+1
      else
        STOP 'wrong i-indexes in sp_ipointers'
      endif
    enddo
    do i=2,A%Nrow+1
      IP(i)=IP(i-1)+IP(i)   !points to [entry-1] of first element with
    enddo                   !certain i-index in A
!***************
    if (optsort == sort) then
      ALLOCATE(xa(a%idata(nelm)))
      do i=1,A%idata(Nelm)
        XA(i)=A%selm1(i)
      enddo
      do i=1,A%idata(Nelm)
        IP(XA(i)%i)=IP(XA(i)%i)+1
        A%selm1(IP(XA(i)%i))=XA(i)
      enddo
      DEALLOCATE(xa)
      IP(2:A%Nrow+1)=IP(1:A%Nrow) !REESTABLISH IP
      IP(1)=0
    endif
  end subroutine sp_ipointers

!> \brief Local routine for this module. Set up pointers for sparse matrices.
  subroutine sp_jpointers(optsort,a,jp)
    implicit none
    integer, intent(in) :: optsort
    type(matrix), intent(inout) :: a
    integer, intent(out) :: jp(a%ncol+1)
    type(smat1), allocatable :: xa(:) !scratch array
    integer :: i

!*******************************
!Set up pointers
!*******************************
    JP=0
    do i=1,A%idata(Nelm)
      if (A%selm1(i)%j>0) then
        JP(A%selm1(i)%j+1)=JP(A%selm1(i)%j+1)+1
      else
        STOP 'wrong j-indexes in sp_jsort'
      endif
    enddo
    do i=2,A%Ncol+1
      JP(i)=JP(i-1)+JP(i)   !points to [entry-1] of first element with
    enddo                   !certain j-index in A
!***************
    if (optsort == sort) then
      ALLOCATE(xa(a%idata(nelm)))
      do i=1,A%idata(Nelm)
        XA(i)=A%selm1(i)
      enddo
      do i=1,A%idata(Nelm)
        JP(XA(i)%j)=JP(XA(i)%j)+1
        A%selm1(JP(XA(i)%j))=XA(i)
      enddo
      DEALLOCATE(xa)
      JP(2:A%Ncol+1)=JP(1:A%Ncol) !REESTABLISH JP
      JP(1)=0
    endif
  end subroutine sp_jpointers


!*****************************************************************************
!Routines needed for purification
! - commented out because purification is not documented and no one really knows
! if purification works!

!  subroutine mat_sparse1_cholesky(A,B)
!  !returns cholesky decomposition of matrix A in B
!    use sparse_matrix
!    implicit none
!    type(Matrix), intent(in) :: A
!    type(Matrix), intent(inout) :: B
!    integer :: i,j,k,n,ndim,jj
!    real(realk) :: sum
!    real(realk), dimension(:), allocatable :: p
!    type(spmatrix) :: C
!    type(row_item) :: j_and_val
!    type(row_node), pointer :: current,current_i,current_j
!    logical :: found
!
!    !initialise and fill a new (linked-list) sparse matrix with A
!    ndim=A%nrow
!    call spmatrix_construct(C,ndim)
!    do n=1,A%idata(nelm)
!       i=A%selm1(n)%i
!       j_and_val%j=A%selm1(n)%j
!       j_and_val%value=A%selm1(n)%val
!       call row_insert(C%col(i),j_and_val)
!    enddo
!
!    !allocate p for the diagonal elements
!    allocate(p(ndim))
!
!    !cholesky decomposition of C
!    do i=1,ndim
!       do j=i,ndim
!          sum=0E0_realk
!          current => C%col(i)%first%next
!          do
!             if(.not. associated(current)) exit
!             jj=current%row%j
!             if(jj.eq.j) then
!                sum=current%row%value
!                exit
!             endif
!             current=>current%next
!          enddo
!
!          current_i => C%col(i)%first%next
!          current_j => C%col(j)%first%next
!
!          do
!             if(.not. associated(current_i)) exit
!             if(.not. associated(current_j)) exit
!             if(current_i%row%j.lt.current_j%row%j) then
!                current_i => current_i%next
!             elseif(current_i%row%j.gt.current_j%row%j) then
!                current_j => current_j%next
!             else
!                k=current_i%row%j
!                if(k.lt.i) then
!                   sum=sum-current_i%row%value*current_j%row%value
!                endif
!                current_i => current_i%next
!                current_j => current_j%next
!             endif
!          enddo
!
!          if(i.eq.j) then
!             if(sum.le. 0E0_realk) stop 'choldc failed'
!             j_and_val%j=i
!             j_and_val%value=sqrt(sum)
!             p(i)=j_and_val%value
!             call row_insert(C%col(i),j_and_val)
!          else
!             j_and_val%j=i
!             j_and_val%value=sum/p(i)
!             call row_insert(C%col(j),j_and_val)
!          endif
!       enddo
!    enddo    
!    deallocate(p)
!
!    !set upper triangle to zero
!    do i=1,ndim
!       current => C%col(i)%first%next
!       do
!          if(.not. associated(current)) exit
!          j=current%row%j
!          if(j.gt.i) call row_delete(C%col(i),current%row,found)
!          current=>current%next
!       enddo
!    enddo
!
!    !call spmatrix_print(C)
!
!    !transfer result to B
!    B%idata(nelm)=spmatrix_number_of_elements(C)
!    call sp_allocate(B,B%idata(nelm))
!    n=0
!    do i=1,ndim
!       current => C%col(i)%first%next
!       do
!          if(.not. associated(current)) exit
!          n=n+1
!          B%selm1(n)%i = i
!          B%selm1(n)%j = current%row%j
!          B%selm1(n)%val = current%row%value
!
!          current => current%next
!       enddo
!    enddo
!    
!    call spmatrix_destruct(C)
!
!    return
!  end subroutine mat_sparse1_cholesky
!
!
!  subroutine mat_sparse1_inverse_triang(A,B)
!  !returns inverse of lower triangular matrix A in B
!    use sparse_matrix
!    implicit none
!    type(Matrix), intent(in) :: A
!    type(Matrix), intent(inout) :: B
!    integer :: i,j,k,ii,kk,n,ndim
!    real(realk) :: sum
!    type(spmatrix) :: C
!    type(row_item) :: j_and_val
!    type(row_node), pointer :: current,current_i,current_j,current_k
!
!    !initialise and fill a new (linked-list) sparse matrix with A
!    ndim=A%nrow
!    call spmatrix_construct(C,ndim)
!    do n=1,A%idata(nelm)
!       i=A%selm1(n)%i
!       j_and_val%j=A%selm1(n)%j
!       j_and_val%value=A%selm1(n)%val
!       call row_insert(C%col(i),j_and_val)
!    enddo
!
!    !inverse of cholesky matrix
!    do i=1,ndim
!       current => C%col(i)%first%next
!       do
!          if(.not. associated(current)) exit
!          j=current%row%j
!          if(j.eq.i) then
!             j_and_val%j=i
!             j_and_val%value=1E0_realk/current%row%value
!             call row_insert(C%col(i),j_and_val)
!             exit
!          endif
!          current=>current%next
!       enddo
!
!       do j=i+1,ndim
!          sum=0E0_realk
!          do k=i,j-1
!             current_j => C%col(j)%first%next
!             do
!                current_k => C%col(k)%first%next
!                if(.not. associated(current_j)) exit
!                do
!                   if(.not. associated(current_k)) exit
!
!                   kk=current_j%row%j
!                   ii=current_k%row%j
!                   if(kk.eq.k.and.ii.eq.i) then
!                      sum=sum-current_j%row%value*current_k%row%value
!                   endif
!                   current_k => current_k%next 
!                enddo
!                current_j => current_j%next
!             enddo
!
!          enddo
!          if(sum.eq. 0E0_realk) cycle    
!          current_j => C%col(j)%first%next
!          do
!             if(.not. associated(current_j)) exit
!             if(current_j%row%j.eq.j) then
!                j_and_val%j=i
!                j_and_val%value=sum/current_j%row%value     
!                exit
!             endif
!             current_j => current_j%next
!          enddo
!          call row_insert(C%col(j),j_and_val)
!       enddo
!
!    enddo
!
!    !transfer result to B
!    B%idata(nelm)=spmatrix_number_of_elements(C)
!    call sp_allocate(B,B%idata(nelm))
!    n=0
!    do i=1,ndim
!       current => C%col(i)%first%next
!       do
!          if(.not. associated(current)) exit
!          n=n+1
!          B%selm1(n)%i = i
!          B%selm1(n)%j = current%row%j
!          B%selm1(n)%val = current%row%value
!
!          current => current%next
!       enddo
!    enddo    
!
!    call spmatrix_destruct(C)
!
!    return
!  end subroutine mat_sparse1_inverse_triang
!  
!
!  subroutine mat_sparse1_simtran(A,B,transb,C)
!  !returns similarity transformation B.A.B^T in C
!    implicit none
!    character, intent(in) :: transb
!    type(Matrix), intent(in) :: A,B
!    type(Matrix), intent(inout) :: C
!    type(Matrix) :: D
!
!    call mat_sparse1_init(D,A%nrow,A%ncol)
!
!    if(transb.eq.'t'.or.transb.eq.'T') then
!       call mat_sparse1_mul(B,A,'T','N',1E0_realk,0E0_realk,D)
!       call mat_sparse1_mul(D,B,'N','N',1E0_realk,0E0_realk,C)
!    else
!       call mat_sparse1_mul(B,A,'N','N',1E0_realk,0E0_realk,D)
!       call mat_sparse1_mul(D,B,'N','T',1E0_realk,0E0_realk,C)
!    endif
!
!    call mat_sparse1_free(D)
!
!    return
!  end subroutine mat_sparse1_simtran
!
!
!  subroutine mat_sparse1_gershgorin_minmax(A,min,max)
!    implicit none
!    type(Matrix), intent(in) :: A
!    real(realk), intent(out) :: min,max
!    integer :: i,j,n
!    real(realk) :: tmin,tmax
!    real(realk), dimension(:), allocatable :: offsum,p
!
!    allocate(offsum(A%nrow),p(A%nrow))
!    offsum=0E0_realk
!
!    do n=1,A%idata(nelm)
!       i=A%selm1(n)%i
!       j=A%selm1(n)%j
!       if(i.eq.j) then
!          p(i)=A%selm1(n)%val
!       else
!          offsum(i)=offsum(i)+abs(A%selm1(n)%val)
!       endif       
!    enddo
!
!    min=p(1)-offsum(1)
!    max=p(1)+offsum(1)
!    do i=2,A%nrow
!       tmin=p(i)-offsum(i)
!       tmax=p(i)+offsum(i)
!       if(min.gt.tmin) min=tmin
!       if(max.lt.tmax) max=tmax
!    enddo
!    
!    deallocate(offsum,p)
!
!    return
!  end subroutine mat_sparse1_gershgorin_minmax
!  
!
!  subroutine mat_sparse1_clean(A,lowcut)
!  !symmetrise and filter small values of matrix
!    use sparse_matrix
!    implicit none
!    type(Matrix), intent(inout) :: A
!    real(realk), intent(in) :: lowcut
!    integer :: i,j,n,ndim
!    type(spmatrix) :: B
!    type(row_item) :: j_and_val
!    type(row_node), pointer :: current
!    logical :: found
!
!    !initialise and fill a new (linked-list) sparse matrix with A
!    ndim=A%nrow
!    call spmatrix_construct(B,ndim)
!    do n=1,A%idata(nelm)
!       i=A%selm1(n)%i
!       j_and_val%j=A%selm1(n)%j
!       j_and_val%value=A%selm1(n)%val
!       call row_insert(B%col(i),j_and_val)
!    enddo    
!
!    !delete zeros and copy lower to upper triangle
!    do i=1,ndim
!       current => B%col(i)%first%next
!       do
!          if(.not. associated(current)) exit
!          j=current%row%j
!          if(j.gt.i) exit
!          j_and_val%j=i
!          j_and_val%value=current%row%value
!          if(abs(j_and_val%value).gt.lowcut) then
!             call row_insert(B%col(j),j_and_val)
!          else
!             call row_delete(B%col(i),current%row,found)
!          endif
!          current=>current%next
!       enddo
!    enddo
!
!    !transfer result to A
!    call sp_deallocate(A)
!    A%idata(nelm)=spmatrix_number_of_elements(B)
!    call sp_allocate(A,A%idata(nelm))
!    n=0
!    do i=1,ndim
!       current => B%col(i)%first%next
!       do
!          if(.not. associated(current)) exit
!          n=n+1
!          A%selm1(n)%i = i
!          A%selm1(n)%j = current%row%j
!          A%selm1(n)%val = current%row%value
!
!          current => current%next
!       enddo
!    enddo
!
!    call spmatrix_destruct(B)
!
!    return
!  end subroutine mat_sparse1_clean
!
!> \brief See mat_sum in mat-operations.f90
  function mat_sparse1_sum(A)
  !returns sum of all elements of matrix
    implicit none
    type(Matrix), intent(in) :: A
    real(realk) :: mat_sparse1_sum
    integer :: n

    mat_sparse1_sum=0E0_realk
    do n=1,A%idata(nelm)
       mat_sparse1_sum=mat_sparse1_sum+A%selm1(n)%val
    enddo

    return
  end function mat_sparse1_sum

end module matrix_operations_sparse1
