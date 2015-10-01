!> @file
!> Simple structure to store and handle three-dimensional arrays
!> \author Janus Juul Eriksen (after 4d template by Marcin Ziolkowski)

!> Three-dimensional array operations
module array3_simple_operations


  ! Outside DEC directory
  use precision
  use dec_typedef_module
  use files!,only: lsopen,lsclose
  use LSTIMING!,only:lstimer
  use memory_handling!, only: mem_alloc,mem_dealloc
  use reorder_frontend_module


  ! DEC DEPENDENCIES (within deccc directory)                                                         
  ! *****************************************
  use array3_memory_manager!, only: memory_allocate_3d, memory_deallocate_3d
!  use dec_fragment_utils

  !> Number of array
  integer(kind=long) :: ArrayNumber=0
  !> Number of created arrays
  integer(kind=long) :: CreatedArrays=0
  !> Number of destroyed arrays
  integer(kind=long) :: DestroyedArrays=0

  ! for timing
  real(realk) :: time_array3_init = 0.0E0_realk
  real(realk) :: time_array3_free = 0.0E0_realk

  private

contains

  !> \brief Initialize three-dimensional array
  !> \author Janus Eriksen (after template by Marcin Ziolkowski)
  !> \param dims Dimensions
  !> \return Three-dimensional array
  function array3_init_standard(dims) result(array)

    implicit none
    type(array3) :: array
    integer, dimension(3), intent(in) :: dims
    real(realk) :: t0,t1

    call cpu_time(t0)

    CreatedArrays = CreatedArrays+1
    array%dims=dims
    array%order=[1,2,3]

    call memory_allocate_3d(array%val,dims)
    call ls_dzero(array%val,size(array%val))

    call cpu_time(t1)
    time_array3_init = time_array3_init + (t1-t0)

    return
  end function array3_init_standard

  subroutine array3_free(array)

    implicit none
    type(array3), intent(inout) :: array
    real(realk) :: t0,t1

    call cpu_time(t0)

    DestroyedArrays = DestroyedArrays + 1

    call memory_deallocate_3d(array%val)

    call cpu_time(t1)
    time_array3_free  = time_array3_free + (t1-t0)

    return
  end subroutine array3_free

  !> \brief Create a copy of an array
  !> \author Marcin Ziolkowski (modified by Kasper Kristensen (+ Janus Juul Eriksen))
  !> \param this three-dimensional array to be copied
  !> \return New array that is a duplicate of the requested array
  function array3_duplicate(this) result(res)
    implicit none
    type(array3), intent(in) :: this
    type(array3) :: res
    integer :: vec_size
    integer(kind=long) :: vec_size64

    res = array3_init_standard(this%dims)

    ! Same ordering for the two arrays
    res%order = this%order

    ! copy data
    DataInMemory: if(associated(this%val)) then

       vec_size64 = int(this%dims(1)*this%dims(2)*this%dims(3),kind=8)
       if(vec_size64>MAXINT)then
          call lsquit('ERROR(array4_duplicate): size of array cannot be &
                       &described by current integer type, please try another &
                       &compilation or fix this routine', DECinfo%output)
       endif 
       vec_size = this%dims(1)*this%dims(2)*this%dims(3)
       if(vec_size < 1) then
          call lsquit('array3_duplicate: Something wrong with vec_size. &
               & Perhaps vec_size is larger than maxium allowed integer!', DECinfo%output)
       end if
       call dcopy(vec_size,this%val,1,res%val,1)

    end if DataInMemory

  end function array3_duplicate

  !> \brief Add scaled arrays when array elements are kept in memory
  subroutine array3_add_memory(alpha,A,beta,B,res_array,zero,zeta,Z)
    implicit none
    type(array3), intent(inout) :: res_array
    type(array3), intent(in) :: A,B
    type(array3), optional, intent(in) :: Z
    real(realk), intent(in) :: alpha,beta
    real(realk), optional, intent(in) :: zeta
    logical, intent(in) :: zero
    integer :: nelements
    integer(kind=long) :: nelements64

    nelements64=int(res_array%dims(1)*res_array%dims(2)*res_array%dims(3),kind=8)
    if(nelements64>MAXINT)then
       call lsquit('ERROR(array3_add_memory): size of array cannot be &
                    &described by current integer type, please try another &
                    &compilation or fix this routine', DECinfo%output)
    endif
    nelements=res_array%dims(1)*res_array%dims(2)*res_array%dims(3)
    if(zero) call ls_dzero(res_array%val,nelements)
    if (.not. present(zeta) .and. (.not. present(Z))) then
       !#ifdef USE_BLAS
       call daxpy(nelements,alpha,A%val,1,res_array%val,1)
       call daxpy(nelements,beta,B%val,1,res_array%val,1)
       !#else
       !    res_array%val=alpha*A%val+beta*B%val
       !#endif
    else if (present(zeta) .and. present(Z)) then
       !#ifdef USE_BLAS
       call daxpy(nelements,alpha,A%val,1,res_array%val,1)
       call daxpy(nelements,beta,B%val,1,res_array%val,1)
       call daxpy(nelements,zeta,Z%val,1,res_array%val,1)
       !#else
       !    res_array%val=alpha*A%val+beta*B%val+zeta*Z%val
       !#endif
    end if

  end subroutine array3_add_memory

  !> \brief Scale B and add to A when array elements are stored in memory
  subroutine array3_add_to_memory(A,beta,B)
    implicit none
    type(array3), intent(inout) :: A
    type(array3), intent(in) :: B
    real(realk), intent(in) :: beta
    integer :: nelements
    integer(kind=long) :: nelements64

    nelements64=int(A%dims(1)*A%dims(2)*A%dims(3),kind=8)
    if(nelements64>MAXINT)then
       call lsquit('ERROR(array3_add_to_memory): size of array cannot be &
                    &described by current integer type, please try another &
                    &compilation or fix this routine', DECinfo%output)
    endif
    nelements=A%dims(1)*A%dims(2)*A%dims(3)
    call daxpy(nelements,beta,B%val,1,A%val,1)

  end subroutine array3_add_to_memory

  !> \brief Contraction of array3 and array2 summing the first (trans) or second index (.not. trans)
  subroutine array3_contract1(A,B,C,zero,trans)

    implicit none
    type(array3), intent(in) :: A
    type(array3), intent(inout) :: C
    type(array2), intent(in) :: B
    logical, intent(in) :: zero, trans
    real(realk) :: MaxElement1,MaxElement2
    real(realk) :: beta
    real(realk) :: tcpu1,twall1,tcpu2,twall2
    integer :: k,l
    integer :: vec_size
    integer(kind=long) :: vec_size64

    call LSTIMER('START',tcpu1,twall1,DECinfo%output)

    if (trans) then
       if(A%dims(1) /= B%dims(1)) then
          print *,'index 1 -> ',A%dims(1),B%dims(1)
          stop 'Error :: Contraction array3_contract1 :: Dimensions do not match'
       end if
    else if(.not. trans) then
       if(A%dims(1) /= B%dims(2)) then
          print *,'index 1 -> ',A%dims(1),B%dims(2)
          stop 'Error :: Contraction array3_contract1 :: Dimensions do not match'
       end if
    end if

    vec_size64 = int(C%dims(1)*C%dims(2)*C%dims(3),kind=8)
    if(vec_size64>MAXINT)then
       call lsquit('ERROR(array3_contract1): size of array cannot be &
                    &described by current integer type, please try another &
                    &compilation or fix this routine', DECinfo%output)
    endif
    vec_size = C%dims(1)*C%dims(2)*C%dims(3)
    if(zero) call ls_dzero(C%val,vec_size)

    if( (A%dims(2) /= C%dims(2)) .or. (A%dims(3) /= C%dims(3))) then
       write(DECinfo%output,*) 'A%dims', A%dims
       write(DECinfo%output,*) 'C%dims', C%dims
       call lsquit('array3_contract1: A and C dimensions do not match',-1)
    end if

    if (trans) then
       do l=1,A%dims(3)
          call dgemm('t','n',B%dims(2),A%dims(2),B%dims(1), &
                   & 1.0E0_realk,B%val,B%dims(1),A%val(:,:,l),A%dims(1), &
                   & 1.0E0_realk,C%val(:,:,l),C%dims(1))
       end do
    else if (.not. trans) then
       do l=1,A%dims(3)
          call dgemm('n','n',B%dims(1),A%dims(2),B%dims(2), &
             &  1.0E0_realk,B%val,B%dims(1),A%val(:,:,l),A%dims(1), &
             &  1.0E0_realk,C%val(:,:,l),C%dims(1))
       end do
    end if

    call LSTIMER('START',tcpu2,twall2,DECinfo%output)

  end subroutine array3_contract1


  !> \brief Contract with two indices (two first)
  subroutine array3_contract2(A,B,C)

    implicit none
    type(array3), intent(in) :: A,B
    type(array2), intent(inout) :: C
    integer :: i,j,k,l,m,n
    integer :: dim1,dim2,dim3
    real(realk) :: tcpu1,twall1,tcpu2,twall2

    call LSTIMER('START',tcpu1,twall1,DECinfo%output)

    ! contractions
    if(A%dims(1) /= B%dims(1) .or. A%dims(2) /= B%dims(2)) then
       do i=1,2
          print *,'index ',i,' -> ',A%dims(i),B%dims(i)
       end do
       stop 'error :: contraction 2 :: dimensions do not match'
    end if

    ! check the dimensions of the final array
    if(C%dims(1) /= A%dims(3) .or. C%dims(2) /= B%dims(3)) then
       print *,'-- Final array --'
       print *,' C(1) ',C%dims(1), ' /= A(3) ',A%dims(3)
       print *,' C(2) ',C%dims(2), ' /= B(3) ',B%dims(3)
       stop 'error :: contraction 2 :: final arrays dimensions do not match'
    end if

    !#ifdef USE_BLAS
    dim1=A%dims(3)
    dim2=A%dims(1)*A%dims(2)
    dim3=B%dims(3)
    call dgemm('t','n',dim1,dim3,dim2, &
         1.0E0_realk,A%val,dim2,B%val,dim2,0.0E0_realk,C%val,dim1)

!!$#else
!!$       do m=1,B%dims(3)
!!$             do k=1,A%dims(3)
!!$
!!$                ! Contraction
!!$                do j=1,A%dims(2)
!!$                   do i=1,A%dims(1)
!!$                      C%val(k,m) = C%val(k,m) + &
!!$                           A%val(i,j,k)*B%val(i,j,m)
!!$                   end do
!!$                end do
!!$
!!$             end do
!!$       end do
!!$#endif

    call LSTIMER('START',tcpu2,twall2,DECinfo%output)

  end subroutine array3_contract2

  !> \brief Reorder indices with additional memory allocation
  !> \author Marcin Ziolkowski
  subroutine array3_reorder(array,order)

    implicit none
    type(array3), intent(inout) :: array
    integer, dimension(3), intent(in) :: order
    integer, dimension(3) :: new_dims,new_order,order1,order2
    real(realk), pointer :: new_data(:,:,:)
    integer :: a,b,c,d
    integer :: dim1,dim2,dim3
    integer :: i,j
    integer :: aa,bb,cc,dd
    integer :: order_type,m,n
    real(realk) :: tcpu1,twall1,tcpu2,twall2

    call LSTIMER('START',tcpu1,twall1,DECinfo%output)

    ! get new order
    do i=1,3
!!$#ifdef VAR_LSDEBUG
!!$       print *, array%order(i),'(',array%dims(i), &
!!$            ') => ',array%order(order(i)),' (',array%dims(order(i)),')'
!!$#endif
       new_order(i) = array%order(order(i))
       new_dims(i) = array%dims(order(i))
    end do

    ! Allocate space for reordered data

    ! Note: To make things work when array3_reorder is called in parallel
    !       we need to make a manual allocation of new_data (not using memory_allocate_3d).
    call mem_alloc(new_data,new_dims(1),new_dims(2),new_dims(3) )
    call ls_dzero(new_data,size(new_data))

    call array_reorder_3d(1.0E0_realk,array%val,array%dims(1),array%dims(2),&
    &array%dims(3),order,0.0E0_realk,new_data)

    ! assign new data
    call memory_deallocate_3d(array%val)
    array%dims=new_dims
    array%order(1)=1
    array%order(2)=2
    array%order(3)=3
    call memory_allocate_3d(array%val,array%dims)
    if(size(array%val) /= size(new_data)) then
       write(DECinfo%output,*) 'array%dims', array%dims
       write(DECinfo%output,*) 'new_dims  ', new_dims
       write(DECinfo%output,*) 'size(array%val)', size(array%val)
       write(DECinfo%output,*) 'size(new_data)', size(new_data)
       call lsquit('array3_reorder: size of array and new_data do no match!',-1)
    end if
    array%val = new_data
    call mem_dealloc(new_data)

    call LSTIMER('START',tcpu2,twall2,DECinfo%output)

  end subroutine array3_reorder


end module array3_simple_operations
