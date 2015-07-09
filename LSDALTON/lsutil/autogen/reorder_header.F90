
  interface get_midx
    module procedure get_mode_idx8,get_mode_idx4
  end interface get_midx

  interface get_cidx
    module procedure get_comp_idx
  end interface get_cidx

  !interface array_reorder
  !  module procedure array_reorder_4d,&
  !                  &array_reorder_3d,&
  !                  &array_reorder_2d
  !end interface array_reorder

#ifdef VAR_OPENACC
  interface array_reorder_acc
    module procedure array_reorder_4d_acc,array_reorder_3d_acc
  end interface array_reorder_acc
#endif

  contains

  !> \author Patrick Ettenhuber
  !> \date June 2012
  !> \brief get mode index from composite index
  subroutine get_mode_idx8(a,inds,dims,modes)
    implicit none
    integer(kind=8),intent(in):: a
    integer,intent(in)        :: dims(*),modes
    integer,intent(inout)     :: inds(*)
    integer(kind=8)           :: i,cind,ndim
    select case(modes)
    case default
      cind=a
      do i=1,modes
        ndim = dims(i)
        inds(i)=mod(cind-1,ndim)+1
        cind=(cind-inds(i))/dims(i) + 1
      enddo
    end select
  end subroutine get_mode_idx8
  !> \author Patrick Ettenhuber
  !> \date June 2012
  !> \brief get mode index from composite index
  subroutine get_mode_idx4(a,inds,dims,modes)
    implicit none
    integer(kind=4),intent(in):: a
    integer,intent(in)        :: dims(*),modes
    integer,intent(inout)     :: inds(*)
    integer(kind=4)           :: i,cind
    select case(modes)
    case default
      cind=a
      do i=1,modes
        inds(i)=mod(cind-1,dims(i))+1
        cind=(cind-inds(i))/dims(i) + 1
      enddo
    end select
  end subroutine get_mode_idx4

  !> \author Patrick Ettenhuber
  !> \date June 2012
  !> \brief get composite index from mode index
  function get_comp_idx(inds,dims,modes) result(a)
    implicit none
    integer,intent(in) :: inds(modes),dims(modes),modes
    integer :: i,j,cdim,cd(modes-2)
    integer :: a,b
    select case(modes)
    case(2)
      a=inds(1)+(inds(2)-1)*dims(1)
    case(3)
      cd(1) = dims(1)*dims(2)
      a=inds(1)+(inds(2)-1)*dims(1)+(inds(3)-1)*cd(1)
    case(4)
      cd(1) = dims(1)*dims(2)
      cd(2) = dims(1)*dims(2)*dims(3)
      a=inds(1)+(inds(2)-1)*dims(1)+(inds(3)-1)*cd(1)+(inds(4)-1)*cd(2)
    case default
      a=1
      do i=1,modes
        cdim=1
        do j=i-1,1,-1
          cdim=cdim*dims(j)
        enddo
        a=a+(inds(i)-1)*cdim
      enddo
    end select
  end function get_comp_idx
  !> \brief general array reordering routine, can add to destination matrix
  !> \author Patrick Ettenhuber, adapted scheme from Marcin Ziolkowski
  subroutine array_reorder_4d(pre1,array_in,d1,d2,d3,d4,order,pre2,array_out)
    implicit none
    integer,intent(in) ::        d1,d2,d3,d4
    real(realk), intent(in)::    array_in(i8*d1*d2*d3*d4),pre1,pre2
    real(realk), intent(inout):: array_out(i8*d1*d2*d3*d4)
    integer, dimension(4), intent(in) :: order

    integer, dimension(4) :: new_order,order1,order2,dims
    integer :: a,b,c,d,maxdim
    integer :: i,j,l
    integer :: aa,bb,cc,dd,fina,finb,finc,find
    integer :: order_type,m,n
    integer :: di3(3), di2(2)
    real(realk) :: tcpu1,twall1,tcpu2,twall2
    integer :: vec_size
    integer(kind=long) :: vec_size64

    vec_size64 = int(d1*d2*d3*d4,kind=8)
    if(vec_size64>MAXINT)then
       call lsquit('ERROR(array_reorder_4d): size of array cannot be &
                    &described by current integer type, please try another &
                    &compilation or fix this routine', -1)
    endif
    vec_size = d1*d2*d3*d4

    call LSTIMER('START',tcpu1,twall1,-1)
    dims(1)=d1
    dims(2)=d2
    dims(3)=d3
    dims(4)=d4

    ! select  invalid reordering type
    order_type = -1

    !MAPPING TO LOWER ORDER REORDERINGS
    if(order(1)==1 .and. order(2)==2 .and. &
         order(3)==3 .and. order(4)==4) order_type = 0
    !CASE 2 D REORDERINGS
    ! rephrase to 2 1
    if(order(1)==3 .and. order(2)==4 .and. &
         order(3)==1 .and. order(4)==2) order_type = 1
    if(order(1)==4 .and. order(2)==1 .and. &
         order(3)==2 .and. order(4)==3) order_type = 2
    if(order(1)==2 .and. order(2)==3 .and. &
         order(3)==4 .and. order(4)==1) order_type = 3
    !CASE 3 D REORDERINGS
    ! rephrase to 1 3 2
    if(order(1)==1 .and. order(2)==2 .and. &
         order(3)==4 .and. order(4)==3) order_type = 4
    if(order(1)==1 .and. order(2)==4 .and. &
         order(3)==2 .and. order(4)==3) order_type = 5
    if(order(1)==1 .and. order(2)==3 .and. &
         order(3)==4 .and. order(4)==2) order_type = 6
    ! rephrase to 2 1 3
    if(order(1)==3 .and. order(2)==1 .and. &
         order(3)==2 .and. order(4)==4) order_type = 7
    if(order(1)==2 .and. order(2)==3 .and. &
         order(3)==1 .and. order(4)==4) order_type = 8
    if(order(1)==2 .and. order(2)==1 .and. &
         order(3)==3 .and. order(4)==4) order_type = 9
    ! rephrase to 3 2 1
    if(order(1)==4 .and. order(2)==3 .and. &
         order(3)==1 .and. order(4)==2) order_type = 10
    if(order(1)==4 .and. order(2)==2 .and. &
         order(3)==3 .and. order(4)==1) order_type = 11
    if(order(1)==3 .and. order(2)==4 .and. &
         order(3)==2 .and. order(4)==1) order_type = 12
    !CASE REAL 4 D REORDERINGS
    if(order(1)==2 .and. order(2)==4 .and. &
         order(3)==1 .and. order(4)==3) order_type = 13
    if(order(1)==3 .and. order(2)==2 .and. &
         order(3)==1 .and. order(4)==4) order_type = 14
    if(order(1)==1 .and. order(2)==3 .and. &
         order(3)==2 .and. order(4)==4) order_type = 15
    if(order(1)==4 .and. order(2)==1 .and. &
         order(3)==3 .and. order(4)==2) order_type = 16
    if(order(1)==2 .and. order(2)==1 .and. &
         order(3)==4 .and. order(4)==3) order_type = 17
    if(order(1)==4 .and. order(2)==3 .and. &
         order(3)==2 .and. order(4)==1) order_type = 18
    if(order(1)==2 .and. order(2)==4 .and. &
         order(3)==3 .and. order(4)==1) order_type = 19
    if(order(1)==1 .and. order(2)==4 .and. &
         order(3)==3 .and. order(4)==2) order_type = 20
    if(order(1)==3 .and. order(2)==1 .and. &
         order(3)==4 .and. order(4)==2) order_type = 21
    if(order(1)==3 .and. order(2)==2 .and. &
         order(3)==4 .and. order(4)==1) order_type = 22
    if(order(1)==4 .and. order(2)==2 .and. &
         order(3)==1 .and. order(4)==3) order_type = 23

    ! do the reordering
    TypeOfReordering: select case(order_type)
    case(0)
      ! CASE 1 2 3 4
       if (pre2 /= 0.0E0_realk) then
         call dscal(vec_size,pre2,array_out,1)
         call daxpy(vec_size,pre1,array_in,1,array_out,1)
       else
         call dcopy(vec_size,array_in,1,array_out,1)
         if (pre1 /= 1.0E0_realk) then
           call dscal(vec_size,pre1,array_out,1)
         endif
       endif

    case(1)
       ! CASE 3 4 1 2
       di2(1) = dims(1)*dims(2)
       di2(2) = dims(3)*dims(4)
       call manual_21_reordering(di2,pre1,array_in,pre2,array_out)
    case(2)
       ! CASE 4 1 2 3
       di2(1) = dims(1)*dims(2)*dims(3)
       di2(2) = dims(4)
       call manual_21_reordering(di2,pre1,array_in,pre2,array_out)
    case(3)
       ! CASE 2 3 4 1
       di2(1) = dims(1)
       di2(2) = dims(2)*dims(3)*dims(4)
       call manual_21_reordering(di2,pre1,array_in,pre2,array_out)

    case(4)
       ! CASE  1 2 4 3
       di3(1) = dims(1)*dims(2)
       di3(2) = dims(3)
       di3(3) = dims(4)
       call manual_132_reordering(di3,pre1,array_in,pre2,array_out)
    case(5)
       ! CASE  1 4 2 3
       di3(1) = dims(1)
       di3(2) = dims(2)*dims(3)
       di3(3) = dims(4)
       call manual_132_reordering(di3,pre1,array_in,pre2,array_out)
    case(6)
       ! CASE  1 3 4 2
       di3(1) = dims(1)
       di3(2) = dims(2)
       di3(3) = dims(3)*dims(4)
       call manual_132_reordering(di3,pre1,array_in,pre2,array_out)


    case(7)
       ! CASE 3 1 2 4
       di3(1) = dims(1)*dims(2)
       di3(2) = dims(3)
       di3(3) = dims(4)
       call manual_213_reordering(di3,pre1,array_in,pre2,array_out)
    case(8)
       ! CASE  2 3 1 4
       di3(1) = dims(1)
       di3(2) = dims(2)*dims(3)
       di3(3) = dims(4)
       call manual_213_reordering(di3,pre1,array_in,pre2,array_out)
    case(9)
       ! CASE 2 1 3 4
       di3(1) = dims(1)
       di3(2) = dims(2)
       di3(3) = dims(3)*dims(4)
       call manual_213_reordering(di3,pre1,array_in,pre2,array_out)

    case(10)
       ! CASE  4 3 1 2
       di3(1) = dims(1)*dims(2)
       di3(2) = dims(3)
       di3(3) = dims(4)
       call manual_321_reordering(di3,pre1,array_in,pre2,array_out)
    case(11)
       ! CASE  4 2 3 1
       di3(1) = dims(1)
       di3(2) = dims(2)*dims(3)
       di3(3) = dims(4)
       call manual_321_reordering(di3,pre1,array_in,pre2,array_out)
    case(12)
       ! CASE  3 4 2 1
       di3(1) = dims(1)
       di3(2) = dims(2)
       di3(3) = dims(3)*dims(4)
       call manual_321_reordering(di3,pre1,array_in,pre2,array_out)


    case(13)
       call manual_2413_reordering(dims,pre1,array_in,pre2,array_out)
    case(14)
       call manual_3214_reordering(dims,pre1,array_in,pre2,array_out)
    case(15)
       call manual_1324_reordering(dims,pre1,array_in,pre2,array_out)
    case(16)
       call manual_4132_reordering(dims,pre1,array_in,pre2,array_out)
    case(17)
       call manual_2143_reordering(dims,pre1,array_in,pre2,array_out)
    case(18)
       call manual_4321_reordering(dims,pre1,array_in,pre2,array_out)
    case(19)
       call manual_2431_reordering(dims,pre1,array_in,pre2,array_out)
    case(20)
       call manual_1432_reordering(dims,pre1,array_in,pre2,array_out)
    case(21)
       call manual_3142_reordering(dims,pre1,array_in,pre2,array_out)
    case(22)
       call manual_3241_reordering(dims,pre1,array_in,pre2,array_out)
    case(23)
       call manual_4213_reordering(dims,pre1,array_in,pre2,array_out)
    case default
       print *,'4d_reordering case does not exist, THIS IS IMPOSSIBLE UNLESS&
       & SOMEBODY DID SOMETHING STUPID'
       call lsquit("ERROR(array_reorder_4d):invalid case",-1)
    end select TypeOfReordering


    call LSTIMER('START',tcpu2,twall2,-1)
  end subroutine array_reorder_4d

#ifdef VAR_OPENACC
  !> \brief general gpu array reordering routine, can add to destination matrix
  !> \author Janus Juul Eriksen, adapted scheme from Patrick Ettenhuber and Marcin Ziolkowski
  subroutine array_reorder_4d_acc(pre1,array_in,d1,d2,d3,d4,order,pre2,array_out,async_idx,async_wait)

    use openacc
    implicit none

    integer,intent(in) ::        d1,d2,d3,d4
    real(realk), intent(in)::    array_in(i8*d1*d2*d3*d4),pre1,pre2
    real(realk), intent(inout):: array_out(i8*d1*d2*d3*d4)
    integer, dimension(4), intent(in) :: order
    integer(kind=acc_handle_kind), intent(in) :: async_idx
    integer(kind=acc_handle_kind), intent(in), optional :: async_wait

    integer, dimension(4) :: new_order,order1,order2,dims
    integer :: a,b,c,d,maxdim
    integer :: i,j,l
    integer :: aa,bb,cc,dd,fina,finb,finc,find
    integer :: order_type,m,n
    integer :: di3(3), di2(2)
    real(realk) :: tcpu1,twall1,tcpu2,twall2
    integer(kind=long) :: vec_size64
    logical :: wait_arg
    integer(kind=acc_handle_kind) :: async_idx2

    vec_size64 = int(d1*d2*d3*d4,kind=8)
    if(vec_size64>MAXINT)then
       call lsquit('ERROR(array_reorder_4d_acc): size of array cannot be &
                    &described by current integer type, please try another &
                    &compilation or fix this routine', -1)
    endif

    wait_arg = .false.
    if (present(async_wait)) wait_arg = .true.

    if (wait_arg) then
       async_idx2 = async_wait
    else
       async_idx2 = async_idx
    endif

    call LSTIMER('START',tcpu1,twall1,-1)
    dims(1)=d1
    dims(2)=d2
    dims(3)=d3
    dims(4)=d4

    ! select  invalid reordering type
    order_type = -1

    !MAPPING TO LOWER ORDER REORDERINGS
    if(order(1)==1 .and. order(2)==2 .and. &
         order(3)==3 .and. order(4)==4) order_type = 0
    !CASE 2 D REORDERINGS
    ! rephrase to 2 1
    if(order(1)==3 .and. order(2)==4 .and. &
         order(3)==1 .and. order(4)==2) order_type = 1
    if(order(1)==4 .and. order(2)==1 .and. &
         order(3)==2 .and. order(4)==3) order_type = 2
    if(order(1)==2 .and. order(2)==3 .and. &
         order(3)==4 .and. order(4)==1) order_type = 3
    !CASE 3 D REORDERINGS
    ! rephrase to 1 3 2
    if(order(1)==1 .and. order(2)==2 .and. &
         order(3)==4 .and. order(4)==3) order_type = 4
    if(order(1)==1 .and. order(2)==4 .and. &
         order(3)==2 .and. order(4)==3) order_type = 5
    if(order(1)==1 .and. order(2)==3 .and. &
         order(3)==4 .and. order(4)==2) order_type = 6
    ! rephrase to 2 1 3
    if(order(1)==3 .and. order(2)==1 .and. &
         order(3)==2 .and. order(4)==4) order_type = 7
    if(order(1)==2 .and. order(2)==3 .and. &
         order(3)==1 .and. order(4)==4) order_type = 8
    if(order(1)==2 .and. order(2)==1 .and. &
         order(3)==3 .and. order(4)==4) order_type = 9
    ! rephrase to 3 2 1
    if(order(1)==4 .and. order(2)==3 .and. &
         order(3)==1 .and. order(4)==2) order_type = 10
    if(order(1)==4 .and. order(2)==2 .and. &
         order(3)==3 .and. order(4)==1) order_type = 11
    if(order(1)==3 .and. order(2)==4 .and. &
         order(3)==2 .and. order(4)==1) order_type = 12
    !CASE REAL 4 D REORDERINGS
    if(order(1)==2 .and. order(2)==4 .and. &
         order(3)==1 .and. order(4)==3) order_type = 13
    if(order(1)==3 .and. order(2)==2 .and. &
         order(3)==1 .and. order(4)==4) order_type = 14
    if(order(1)==1 .and. order(2)==3 .and. &
         order(3)==2 .and. order(4)==4) order_type = 15
    if(order(1)==4 .and. order(2)==1 .and. &
         order(3)==3 .and. order(4)==2) order_type = 16
    if(order(1)==2 .and. order(2)==1 .and. &
         order(3)==4 .and. order(4)==3) order_type = 17
    if(order(1)==4 .and. order(2)==3 .and. &
         order(3)==2 .and. order(4)==1) order_type = 18
    if(order(1)==2 .and. order(2)==4 .and. &
         order(3)==3 .and. order(4)==1) order_type = 19
    if(order(1)==1 .and. order(2)==4 .and. &
         order(3)==3 .and. order(4)==2) order_type = 20
    if(order(1)==3 .and. order(2)==1 .and. &
         order(3)==4 .and. order(4)==2) order_type = 21
    if(order(1)==3 .and. order(2)==2 .and. &
         order(3)==4 .and. order(4)==1) order_type = 22
    if(order(1)==4 .and. order(2)==2 .and. &
         order(3)==1 .and. order(4)==3) order_type = 23

    ! do the reordering
    TypeOfReordering4d_acc: select case(order_type)
    case(0)
       ! CASE 1 2 3 4
       print *,'4d_acc_reordering case 1234 - no reordering - do not call this routine'
       call lsquit("ERROR(array_reorder_4d_acc):case 1234 - no reordering - do not call this routine",-1)

    case(1)
       ! CASE 3 4 1 2
       di2(1) = dims(1)*dims(2)
       di2(2) = dims(3)*dims(4)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_21_reordering_0(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_21_reordering_1(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_21_reordering_2(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_21_reordering_3(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_21_reordering_4(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_21_reordering_5(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case(2)
       ! CASE 4 1 2 3
       di2(1) = dims(1)*dims(2)*dims(3)
       di2(2) = dims(4)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_21_reordering_0(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_21_reordering_1(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_21_reordering_2(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_21_reordering_3(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_21_reordering_4(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_21_reordering_5(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case(3)
       ! CASE 2 3 4 1
       di2(1) = dims(1)
       di2(2) = dims(2)*dims(3)*dims(4)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_21_reordering_0(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_21_reordering_1(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_21_reordering_2(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_21_reordering_3(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_21_reordering_4(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_21_reordering_5(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if

    case(4)
       ! CASE  1 2 4 3
       di3(1) = dims(1)*dims(2)
       di3(2) = dims(3)
       di3(3) = dims(4)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_132_reordering_0(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_132_reordering_1(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_132_reordering_2(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_132_reordering_3(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_132_reordering_4(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_132_reordering_5(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case(5)
       ! CASE  1 4 2 3
       di3(1) = dims(1)
       di3(2) = dims(2)*dims(3)
       di3(3) = dims(4)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_132_reordering_0(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_132_reordering_1(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_132_reordering_2(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_132_reordering_3(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_132_reordering_4(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_132_reordering_5(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case(6)
       ! CASE  1 3 4 2
       di3(1) = dims(1)
       di3(2) = dims(2)
       di3(3) = dims(3)*dims(4)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_132_reordering_0(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_132_reordering_1(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_132_reordering_2(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_132_reordering_3(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_132_reordering_4(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_132_reordering_5(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case(7)
       ! CASE 3 1 2 4
       di3(1) = dims(1)*dims(2)
       di3(2) = dims(3)
       di3(3) = dims(4)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_213_reordering_0(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_213_reordering_1(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_213_reordering_2(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_213_reordering_3(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_213_reordering_4(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_213_reordering_5(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case(8)
       ! CASE  2 3 1 4
       di3(1) = dims(1)
       di3(2) = dims(2)*dims(3)
       di3(3) = dims(4)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_213_reordering_0(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_213_reordering_1(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_213_reordering_2(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_213_reordering_3(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_213_reordering_4(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_213_reordering_5(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case(9)
       ! CASE 2 1 3 4
       di3(1) = dims(1)
       di3(2) = dims(2)
       di3(3) = dims(3)*dims(4)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_213_reordering_0(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_213_reordering_1(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_213_reordering_2(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_213_reordering_3(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_213_reordering_4(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_213_reordering_5(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case(10)
       ! CASE  4 3 1 2
       di3(1) = dims(1)*dims(2)
       di3(2) = dims(3)
       di3(3) = dims(4)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_321_reordering_0(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_321_reordering_1(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_321_reordering_2(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_321_reordering_3(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_321_reordering_4(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_321_reordering_5(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case(11)
       ! CASE  4 2 3 1
       di3(1) = dims(1)
       di3(2) = dims(2)*dims(3)
       di3(3) = dims(4)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_321_reordering_0(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_321_reordering_1(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_321_reordering_2(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_321_reordering_3(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_321_reordering_4(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_321_reordering_5(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case(12)
       ! CASE  3 4 2 1
       di3(1) = dims(1)
       di3(2) = dims(2)
       di3(3) = dims(3)*dims(4)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_321_reordering_0(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_321_reordering_1(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_321_reordering_2(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_321_reordering_3(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_321_reordering_4(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_321_reordering_5(di3,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if

    case(13)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_2413_reordering_0(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_2413_reordering_1(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_2413_reordering_2(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_2413_reordering_3(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_2413_reordering_4(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_2413_reordering_5(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case(14)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_3214_reordering_0(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_3214_reordering_1(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_3214_reordering_2(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_3214_reordering_3(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_3214_reordering_4(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_3214_reordering_5(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case(15)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_1324_reordering_0(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_1324_reordering_1(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_1324_reordering_2(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_1324_reordering_3(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_1324_reordering_4(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_1324_reordering_5(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case(16)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_4132_reordering_0(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_4132_reordering_1(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_4132_reordering_2(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_4132_reordering_3(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_4132_reordering_4(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_4132_reordering_5(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case(17)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_2143_reordering_0(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_2143_reordering_1(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_2143_reordering_2(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_2143_reordering_3(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_2143_reordering_4(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_2143_reordering_5(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case(18)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_4321_reordering_0(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_4321_reordering_1(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_4321_reordering_2(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_4321_reordering_3(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_4321_reordering_4(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_4321_reordering_5(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case(19)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_2431_reordering_0(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_2431_reordering_1(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_2431_reordering_2(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_2431_reordering_3(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_2431_reordering_4(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_2431_reordering_5(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case(20)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_1432_reordering_0(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_1432_reordering_1(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_1432_reordering_2(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_1432_reordering_3(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_1432_reordering_4(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_1432_reordering_5(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case(21)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_3142_reordering_0(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_3142_reordering_1(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_3142_reordering_2(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_3142_reordering_3(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_3142_reordering_4(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_3142_reordering_5(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case(22)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_3241_reordering_0(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_3241_reordering_1(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_3241_reordering_2(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_3241_reordering_3(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_3241_reordering_4(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_3241_reordering_5(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case(23)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_4213_reordering_0(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_4213_reordering_1(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_4213_reordering_2(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_4213_reordering_3(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_4213_reordering_4(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_4213_reordering_5(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case default
       print *,'4d_acc_reordering case does not exist, THIS IS IMPOSSIBLE UNLESS&
       & SOMEBODY DID SOMETHING STUPID'
       call lsquit("ERROR(array_reorder_4d_acc):invalid case",-1)
    end select TypeOfReordering4d_acc


    call LSTIMER('START',tcpu2,twall2,-1)
  end subroutine array_reorder_4d_acc
#endif


  !> \brief general 3d array reordering routine, can add to destination matrix
  !> \author Janus Juul Eriksen, adapted scheme from Marcin Ziolkowski (Patrick Ettenhuber)
  subroutine array_reorder_3d(pre1,array_in,d1,d2,d3,order,pre2,array_out)
    implicit none
    integer,intent(in) ::        d1,d2,d3
    real(realk), intent(in)::    array_in(i8*d1*d2*d3),pre1,pre2
    real(realk), intent(inout):: array_out(i8*d1*d2*d3)
    integer, dimension(3), intent(in) :: order

    integer, dimension(3) :: new_order,order1,order2,dims
    integer :: a,b,c,fina,finb,finc
    integer :: aa,bb,cc
    integer :: order_type
    integer :: vec_size
    integer :: di2(2)
    integer(kind=long) :: vec_size64
    real(realk) :: tcpu1,twall1,tcpu2,twall2

    call LSTIMER('START',tcpu1,twall1,-1)

    vec_size64 = int(d1*d2*d3,kind=8)
    if(vec_size64>MAXINT)then
       call lsquit('ERROR(array_reorder_3d): size of array cannot be &
                    &described by current integer type, please try another &
                    &compilation or fix this routine', -1)
    endif
    vec_size = d1*d2*d3

    dims(1)=d1
    dims(2)=d2
    dims(3)=d3

    ! select  general type of the reordering
    order_type = -1
    if(order(1)==1 .and. order(2)==2 .and. &
         order(3)==3) order_type = 0
    !CASE 2 D REORDERINGS
    ! rephrase to 2 1
    if(order(1)==3 .and. order(2)==1 .and. &
         order(3)==2) order_type = 1
    if(order(1)==2 .and. order(2)==3 .and. &
         order(3)==1) order_type = 2
    !REAL 3D REORDERINGS
    if(order(1)==1 .and. order(2)==3 .and. &
         order(3)==2) order_type = 3
    if(order(1)==2 .and. order(2)==1 .and. &
         order(3)==3) order_type = 4
    if(order(1)==3 .and. order(2)==2 .and. &
         order(3)==1) order_type = 5

    ! do the reordering
    TypeOfReordering: select case(order_type)
    case(0)
       if (pre2 .ne. 0.0E0_realk) then
         call dscal(vec_size,pre2,array_out,1)
         call daxpy(vec_size,pre1,array_in,1,array_out,1)
       else
         call dcopy(vec_size,array_in,1,array_out,1)
         if (pre1 .ne. 1.0E0_realk) then
           call dscal(vec_size,pre1,array_out,1)
         endif
       endif
    case(1)
       ! CASE 3 1 2
       di2(1) = dims(1)*dims(2)
       di2(2) = dims(3)
       call manual_21_reordering(di2,pre1,array_in,pre2,array_out)
    case(2)
       ! CASE 2 3 1
       di2(1) = dims(1)
       di2(2) = dims(2) * dims(3)
       call manual_21_reordering(di2,pre1,array_in,pre2,array_out)

    case(3)
       call manual_132_reordering(dims,pre1,array_in,pre2,array_out)
    case(4)
       call manual_213_reordering(dims,pre1,array_in,pre2,array_out)
    case(5)
       call manual_321_reordering(dims,pre1,array_in,pre2,array_out)
    case default
       print *,'3d_reordering case does not exist, THIS IS IMPOSSIBLE UNLESS&
       & SOMEBODY DID SOMETHING STUPID'
       call lsquit("ERROR(array_reorder_3d):invalid case",-1)

    end select TypeOfReordering

    call LSTIMER('START',tcpu2,twall2,-1)

  end subroutine array_reorder_3d


#ifdef VAR_OPENACC
  !> \brief general gpu 3d array reordering routine, can add to destination matrix
  !> \author Janus Juul Eriksen, adapted scheme from Patrick Ettenhuber and Marcin Ziolkowski
  subroutine array_reorder_3d_acc(pre1,array_in,d1,d2,d3,order,pre2,array_out,async_idx,async_wait)

    use openacc
    implicit none

    integer,intent(in) ::        d1,d2,d3
    real(realk), intent(in)::    array_in((i8*d1)*d2*d3),pre1,pre2
    real(realk), intent(inout):: array_out((i8*d1)*d2*d3)
    integer, dimension(3), intent(in) :: order
    integer(kind=acc_handle_kind), intent(in) :: async_idx
    integer(kind=acc_handle_kind), intent(in), optional :: async_wait

    integer, dimension(3) :: new_order,order1,order2,dims
    integer :: a,b,c,fina,finb,finc
    integer :: aa,bb,cc
    integer :: order_type
    integer :: di2(2)
    integer(kind=long) :: vec_size64
    real(realk) :: tcpu1,twall1,tcpu2,twall2
    logical :: wait_arg
    integer(kind=acc_handle_kind) :: async_idx2

    wait_arg = .false.
    if (present(async_wait)) wait_arg = .true.

    if (wait_arg) then 
       async_idx2 = async_wait
    else
       async_idx2 = async_idx
    endif

    vec_size64 = int(d1*d2*d3,kind=8)
    if(vec_size64>MAXINT)then
       call lsquit('ERROR(array_reorder_3d_acc): size of array cannot be &
                    &described by current integer type, please try another &
                    &compilation or fix this routine', -1)
    endif

    call LSTIMER('START',tcpu1,twall1,-1)

    dims(1)=d1
    dims(2)=d2
    dims(3)=d3

    ! select  general type of the reordering
    order_type = -1
    if(order(1)==1 .and. order(2)==2 .and. &
         order(3)==3) order_type = 0
    !CASE 2 D REORDERINGS
    ! rephrase to 2 1
    if(order(1)==3 .and. order(2)==1 .and. &
         order(3)==2) order_type = 1
    if(order(1)==2 .and. order(2)==3 .and. &
         order(3)==1) order_type = 2
    !REAL 3D REORDERINGS
    if(order(1)==1 .and. order(2)==3 .and. &
         order(3)==2) order_type = 3
    if(order(1)==2 .and. order(2)==1 .and. &
         order(3)==3) order_type = 4
    if(order(1)==3 .and. order(2)==2 .and. &
         order(3)==1) order_type = 5

    ! do the reordering
    TypeOfReordering3d_acc: select case(order_type)
    case(0)
       print *,'3d_acc_reordering case 123 - no reordering - do not call this routine'
       call lsquit("ERROR(array_reorder_3d_acc):case 123 - no reordering - do not call this routine",-1)
    case(1)
       ! CASE 3 1 2
       di2(1) = dims(1)*dims(2)
       di2(2) = dims(3)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_21_reordering_0(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_21_reordering_1(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_21_reordering_2(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_21_reordering_3(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_21_reordering_4(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_21_reordering_5(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case(2)
       ! CASE 2 3 1
       di2(1) = dims(1)
       di2(2) = dims(2) * dims(3)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_21_reordering_0(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_21_reordering_1(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_21_reordering_2(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_21_reordering_3(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_21_reordering_4(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_21_reordering_5(di2,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if

    case(3)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_132_reordering_0(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_132_reordering_1(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_132_reordering_2(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_132_reordering_3(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_132_reordering_4(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_132_reordering_5(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case(4)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_213_reordering_0(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_213_reordering_1(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_213_reordering_2(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_213_reordering_3(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_213_reordering_4(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_213_reordering_5(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case(5)
       if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_321_reordering_0(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 0.0E0_realk) then
          call manual_acc_321_reordering_1(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_321_reordering_2(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .eq. 1.0E0_realk) then
          call manual_acc_321_reordering_3(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .eq. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_321_reordering_4(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       else if (pre1 .ne. 1.0E0_realk .and. pre2 .ne. 1.0E0_realk) then
          call manual_acc_321_reordering_5(dims,pre1,array_in,pre2,array_out,async_idx,async_idx2,wait_arg)
       end if
    case default
       print *,'3d_reordering_acc case does not exist, THIS IS IMPOSSIBLE UNLESS&
       & SOMEBODY DID SOMETHING STUPID'
       call lsquit("ERROR(array_reorder_3d_acc):invalid case",-1)

    end select TypeOfReordering3d_acc

    call LSTIMER('START',tcpu2,twall2,-1)

  end subroutine array_reorder_3d_acc
#endif


  !> \brief 2d array reordering routine, for debugging and testing
  !> \author Patrick Ettenhuber
  subroutine array_reorder_2d(pre1,array_in,d1,d2,order,pre2,array_out)
     implicit none
     integer,intent(in) ::        d1,d2
     real(realk), intent(in)::    array_in((i8*d1)*d2),pre1,pre2
     real(realk), intent(inout):: array_out((i8*d1)*d2)
     integer, dimension(2), intent(in) :: order
     if (order(1) == 1 .and. order(2) == 2 )then
        if (pre2 .ne. 0.0E0_realk) then
           !call dscal(d1*d2,pre2,array_out,1)
           !call daxpy(d1*d2,pre1,array_in,1,array_out,1)
           array_out = pre2 * array_out + pre1*array_in
        else
           array_out = pre1*array_in
        endif
        elseif (order(1) == 2 .and. order(2) == 1) then
        call mat_transpose(d1,d2,pre1,array_in,pre2,array_out)
     else
        call lsquit("ERROR(array_reorder_2d): reordering not defined",-1)
     endif
  end subroutine array_reorder_2d
  !\>  \brief another transposition routine, intended to replace the others
  !\>  \> author Patrick Ettenhuber
  subroutine mat_transpose(r,c,p1,x,p2,y)
    implicit none
    integer,intent(in) ::        r,c
    real(realk), intent(in)::    x((i8*r)*c),p1,p2
    real(realk), intent(inout):: y((i8*c)*r)

    integer, dimension(2) :: dims
    real(realk) :: tcpu1,twall1,tcpu2,twall2

    call LSTIMER('START',tcpu1,twall1,-1)

    dims(1)=r
    dims(2)=c

    call manual_21_reordering(dims,p1,x,p2,y)

    call LSTIMER('START',tcpu2,twall2,-1)

  end subroutine mat_transpose

  !> \brief Transpose any (esp. rectangular)  matrix in place  -> algorithm 513
  subroutine alg513(a, m, n, mn, move, iwrk, iok)                    !tra   10
  !*****
  ! algorithm 380 - revised
  !*****
  ! a is a one-dimensional array of length mn=m*n, which
  ! contains the mxn matrix to be transposed (stored
  ! columwise). move is a one-dimensional array of length iwrk
  ! used to store information to speed up the process.  the
  ! value iwrk=(m+n)/2 is recommended. iok indicates the
  ! success or failure of the routine.
  ! normal return  iok=0
  ! errors         iok=-1 ,mn not equal to m*n
  !                iok=-2 ,iwrk negative or zero
  !                iok.gt.0, (should never occur),in this case
  ! we set iok equal to the final value of i when the search
  ! is completed but some loops have not been moved
  ! note * move(i) will stay zero for fixed points
  !      dimension a(mn), move(iwrk)
        implicit none
        integer,intent(in) :: m, n, mn,iwrk
        integer,intent(inout) :: iok 
        real(realk), dimension(mn), intent(inout) :: a
        integer, dimension(iwrk), intent(inout) :: move
        integer :: ncount, k, i, ir1, ir2
        integer :: ir0, im, i2, kmi, i1c, i2c, i1, j, max, n1, j1
        real(realk) :: b, c, d
  !check arguments and initialize.
        if (m.lt.2 .or. n.lt.2) go to 120
        if (mn.ne.m*n) go to 180
        if (iwrk.lt.1) go to 190
        if (m.eq.n) go to 130
        ncount = 2
        k = mn - 1
        do 10 i=1,iwrk
          move(i) = 0
     10 continue
        if (m.lt.3 .or. n.lt.3) go to 30
  !calculate the number of fixed points, euclids algorithm
  !for gcd(m-1,n-1).
        ir2 = m - 1
        ir1 = n - 1
     20 ir0 = mod(ir2,ir1)
        ir2 = ir1
        ir1 = ir0
        if (ir0.ne.0) go to 20
        ncount = ncount + ir2 - 1
  !set initial values for search
     30 i = 1
        im = m
  !at least one loop must be re-arranged
        go to 80
  !search for loops to rearrange
     40 max = k - i
        i = i + 1
        if (i.gt.max) go to 160
        im = im + m
        if (im.gt.k) im = im - k
        i2 = im
        if (i.eq.i2) go to 40
        if (i.gt.iwrk) go to 60
        if (move(i).eq.0) go to 80
        go to 40
     50 i2 = m*i1 - k*(i1/n)
     60 if (i2.le.i .or. i2.ge.max) go to 70
        i1 = i2
        go to 50
     70 if (i2.ne.i) go to 40
  !rearrange the elements of a loop and its companion loop
     80 i1 = i
        kmi = k - i
        b = a(i1+1)
        i1c= kmi
        c = a(i1c+1)
     90 i2 = m*i1 - k*(i1/n)
        i2c = k - i2
        if (i1.le.iwrk) move(i1) = 2
        if (i1c.le.iwrk) move(i1c) = 2
        ncount = ncount + 2
        if (i2.eq.i) go to 110
        if (i2.eq.kmi) go to 100
        a(i1+1) = a(i2+1)
        a(i1c+1) = a(i2c+1)
        i1 = i2
        i1c = i2c
        go to 90
  !final store and test for finished
    100 d = b
        b = c
        c = d
    110 a(i1+1) = b
        a(i1c+1) = c
        if (ncount.lt.mn) go to 40
  !normal return
    120 iok = 0
       return
  !if matrix is square,exchange elements a(i,j) and a(j,i).
    130 n1 = n - 1
        do 150 i=1,n1
          j1 = i + 1
          do 140 j=j1,n
            i1 = i + (j-1)*n
            i2 = j + (i-1)*m
            b = a(i1)
            a(i1) = a(i2)
            a(i2) = b
    140   continue
    150 continue
        go to 120
  !error returns.
    160 iok = i
    170 return
    180 iok = -1
        go to 170
    190 iok = -2
        go to 170
  end subroutine alg513

  !> \autor Patrick Ettenhuber
  !> \brief Subroutine to put a specific tile in a general matrix
  !> \date January 2013
  subroutine tile_in_fort(pre1,tilein,tnr,tdims,pre2,fort,full_arr_dim,mode,o)
     implicit none
     !> scaling factors
     real(realk) :: pre1,pre2
     !> input array with data in arbtirary order
     real(realk), intent(inout) :: fort(*)
     !> mode infortmation about how to interpret data
     integer, intent(in) :: mode
     !> the tile number in column major ordering of the tiles
     integer, intent(in) :: tnr
     !> dimension information for the mew array
     integer, intent(in) :: full_arr_dim(mode)
     !> batch information for the tiles
     integer, intent(in) :: tdims(mode)
     !> specify how to reorder the tile to the new full array
     integer,intent(in) :: o(mode)
     !> tile output
     real(realk), intent(in) :: tilein(*)
     integer :: i,nelms,k
     integer :: tmodeidx(mode)
     integer :: idxintile(mode),ro(mode),full_dim_tiled_arr(mode)
     integer :: ccels,ntimes,acttdim(mode),fels(mode)
     integer :: pos1,ntpm(mode),glbmodeidx(mode)
     integer :: order_type

     order_type=0
     do i=1,mode
        if(o(i)/=i)order_type=-1
        ro(o(i))=i
     enddo


     do i=1,mode
        full_dim_tiled_arr(i) = full_arr_dim(ro(i))
        ntpm(i) = full_dim_tiled_arr(i)/tdims(i)
        if(mod(full_dim_tiled_arr(i),tdims(i))>0)ntpm(i) = ntpm(i) + 1
     enddo


     select case(mode)
     case(4)
        if(o(1)==1.and.o(2)==2.and.o(3)==3.and.o(4)==4)order_type = 0
        if(o(1)==3.and.o(2)==4.and.o(3)==1.and.o(4)==2)order_type = 1
        if(o(1)==4.and.o(2)==1.and.o(3)==2.and.o(4)==3)order_type = 2
        if(o(1)==2.and.o(2)==3.and.o(3)==4.and.o(4)==1)order_type = 3
        if(o(1)==1.and.o(2)==2.and.o(3)==4.and.o(4)==3)order_type = 4
        if(o(1)==1.and.o(2)==4.and.o(3)==2.and.o(4)==3)order_type = 5
        if(o(1)==1.and.o(2)==3.and.o(3)==4.and.o(4)==2)order_type = 6
        if(o(1)==3.and.o(2)==1.and.o(3)==2.and.o(4)==4)order_type = 7
        if(o(1)==2.and.o(2)==3.and.o(3)==1.and.o(4)==4)order_type = 8
        if(o(1)==2.and.o(2)==1.and.o(3)==3.and.o(4)==4)order_type = 9
        if(o(1)==4.and.o(2)==3.and.o(3)==1.and.o(4)==2)order_type = 10
        if(o(1)==4.and.o(2)==2.and.o(3)==3.and.o(4)==1)order_type = 11
        if(o(1)==3.and.o(2)==4.and.o(3)==2.and.o(4)==1)order_type = 12
        if(o(1)==2.and.o(2)==4.and.o(3)==1.and.o(4)==3)order_type = 13
        if(o(1)==3.and.o(2)==2.and.o(3)==1.and.o(4)==4)order_type = 14
        if(o(1)==1.and.o(2)==3.and.o(3)==2.and.o(4)==4)order_type = 15
        if(o(1)==4.and.o(2)==1.and.o(3)==3.and.o(4)==2)order_type = 16
        if(o(1)==2.and.o(2)==1.and.o(3)==4.and.o(4)==3)order_type = 17
        if(o(1)==4.and.o(2)==3.and.o(3)==2.and.o(4)==1)order_type = 18
        if(o(1)==2.and.o(2)==4.and.o(3)==3.and.o(4)==1)order_type = 19
        if(o(1)==1.and.o(2)==4.and.o(3)==3.and.o(4)==2)order_type = 20
        if(o(1)==3.and.o(2)==1.and.o(3)==4.and.o(4)==2)order_type = 21
        if(o(1)==3.and.o(2)==2.and.o(3)==4.and.o(4)==1)order_type = 22
        if(o(1)==4.and.o(2)==2.and.o(3)==1.and.o(4)==3)order_type = 23
     case(3)
        if(o(1)==1.and.o(2)==2.and.o(3)==3)            order_type = 0
        if(o(1)==3.and.o(2)==1.and.o(3)==2)            order_type = 24
        if(o(1)==2.and.o(2)==3.and.o(3)==1)            order_type = 25
        if(o(1)==1.and.o(2)==3.and.o(3)==2)            order_type = 26
        if(o(1)==2.and.o(2)==1.and.o(3)==3)            order_type = 27
        if(o(1)==3.and.o(2)==2.and.o(3)==1)            order_type = 28
     case(2)
        if(o(1)==1.and.o(2)==2)                        order_type = 0
        if(o(1)==2.and.o(2)==1)                        order_type = 29
     end select

     call get_midx(tnr,tmodeidx,ntpm,mode)

     ntimes=1
     do i=1,mode
        fels(i) = (tmodeidx(i)-1) * tdims(i) + 1
        if(tmodeidx(i)*tdims(i)>full_dim_tiled_arr(i))then
           acttdim(i)=mod(full_dim_tiled_arr(i),tdims(i))
        else
           acttdim(i)=tdims(i)
        endif
        if(i>1)ntimes=ntimes*acttdim(i)
     enddo


     select case(order_type)
     case(0)
        ccels=acttdim(1)
        do i=1,ntimes
           call get_midx(i,idxintile(2:mode),acttdim(2:mode),mode-1)
           idxintile(1)=1
           do k=1,mode
              glbmodeidx(k)=idxintile(k) +(tmodeidx(k)-1) *tdims(k)
           enddo
           pos1=get_cidx(glbmodeidx,full_arr_dim,mode)
           if(pre1==1.0E0_realk.and.pre2==0.0E0_realk)then
              call dcopy(ccels,tilein(1+(i-1)*ccels),1,fort(pos1),1)
           else
              call dscal(ccels,pre2,fort(pos1),1)
              call daxpy(ccels,pre1,tilein(1+(i-1)*ccels),1,fort(pos1),1)
           endif
        enddo
     case(1)
        call manual_3412_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
     case(2)
        call manual_4123_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
     case(3)
        call manual_2341_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
     case(4)
        call manual_1243_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
     case(5)
        call manual_1423_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
     case(6)
        call manual_1342_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
     case(7)
        call manual_3124_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
     case(8)
        call manual_2314_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
     case(9)
        call manual_2134_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
     case(10)
        call manual_4312_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
     case(11)
        call manual_4231_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
     case(12)
        call manual_3421_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
     case(13)
        call manual_2413_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
     case(14)
        call manual_3214_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
     case(15)
        call manual_1324_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
     case(16)
        call manual_4132_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
     case(17)
        call manual_2143_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
     case(18)
        call manual_4321_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
     case(19)
        call manual_2431_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
     case(20)
        call manual_1432_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
     case(21)
        call manual_3142_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
     case(22)
        call manual_3241_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
     case(23)
        call manual_4213_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
      case(24)
        call manual_312_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
      case(25)
        call manual_231_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
      case(26)
        call manual_132_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
      case(27)
        call manual_213_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
      case(28)
        call manual_321_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
      case(29)
        call manual_21_reordering_t2f(acttdim,full_dim_tiled_arr,fels,pre1,tilein,pre2,fort)
     case default
        print *,"expensive default tile_in_fort",o
        !print *,"order  :",o
        !print *,"rorder :",ro
        !print *,"fad    :",full_arr_dim
        !print *,"tad    :",full_dim_tiled_arr
        !print *,"atd    :",acttdim
        !count elements in the current tile for loop over elements
        !identify their original position and put them in tile
        nelms=1
        do i=1,mode
           nelms = nelms * acttdim(i)
        enddo

        do i = 1,nelms
           !get mode index of element in tile
           call get_midx(i,idxintile,acttdim,mode)
           !get full index in new array from indices referencing the tile
           do k=1,mode
              glbmodeidx(ro(k))=idxintile(k) + (tmodeidx(k)-1)*tdims(k)
           enddo
           !get index in new array
           pos1=get_cidx(glbmodeidx,full_arr_dim,mode)
           fort(pos1)=pre2*fort(pos1) + pre1 * tilein(i)
        enddo
     end select
  end subroutine tile_in_fort



  !> \autor Patrick Ettenhuber
  !> \brief Subroutine to extract a specific tile from a general matrix
  !> \date November 2012
  subroutine tile_from_fort(pre1,fort,full_arr_dim,mode,pre2,tileout,tnr,tdims,o)
    implicit none
    !> scaling factors
    real(realk) :: pre1,pre2
    !> input array with data in arbtirary order
    real(realk), intent(in) :: fort(*)
    !> mode infortmation about how to interpret data
    integer, intent(in) :: mode
    !> the tile number in column major ordering of the tiles
    integer, intent(in) :: tnr
    !> dimension infortmation about how to interpret data
    integer, intent(in) :: full_arr_dim(mode)
    !> batch information for the tiles
    integer, intent(in) :: tdims(mode)
    !> reorder information for the array with respect to the original array,
    !> if optorder is given, then the dimensions, tdims and tnr are with
    !reference to the tile to calculate
    integer, intent(in) :: o(mode)
    !> tile output
    real(realk), intent(inout) :: tileout(*)
    integer :: i,nccblcks,nels,k
    integer :: tmodeidx(mode)
    integer :: idxintile(mode),glbidx
    integer :: ccels,ntimes,el,acttdim(mode),full_dim_tiled_arr(mode),nelms
    integer :: pos1,pos2,ntpm(mode),glbmodeidx(mode),ro(mode),rtd(mode),fels(mode)
    integer :: order_type,a,b,c,d


    order_type=0
    do i=1,mode
      if(o(i)/=i)order_type=-1
    enddo


    do i=1,mode
      !get the reverse order information
      ro(o(i))=i
    enddo

    !calculate number of tiles per mode
    nels=1
    do i=1,mode
      full_dim_tiled_arr(i)=full_arr_dim(o(i))
      ntpm(i) = full_dim_tiled_arr(i)/tdims(i)
      if(mod(full_dim_tiled_arr(i),tdims(i))>0)ntpm(i) = ntpm(i) + 1
      nels=nels*full_dim_tiled_arr(i)
    enddo

    !print *,"fad    :",full_arr_dim
    !print *,"tad    :",full_dim_tiled_arr
    !print *,"ntpm   :",ntpm
    !print *,"td     :",tdims

    select case(mode)
    case(4)
       if(o(1)==1.and.o(2)==2.and.o(3)==3.and.o(4)==4)order_type = 0
       if(o(1)==3.and.o(2)==4.and.o(3)==1.and.o(4)==2)order_type = 1
       if(o(1)==4.and.o(2)==1.and.o(3)==2.and.o(4)==3)order_type = 2
       if(o(1)==2.and.o(2)==3.and.o(3)==4.and.o(4)==1)order_type = 3
       if(o(1)==1.and.o(2)==2.and.o(3)==4.and.o(4)==3)order_type = 4
       if(o(1)==1.and.o(2)==4.and.o(3)==2.and.o(4)==3)order_type = 5
       if(o(1)==1.and.o(2)==3.and.o(3)==4.and.o(4)==2)order_type = 6
       if(o(1)==3.and.o(2)==1.and.o(3)==2.and.o(4)==4)order_type = 7
       if(o(1)==2.and.o(2)==3.and.o(3)==1.and.o(4)==4)order_type = 8
       if(o(1)==2.and.o(2)==1.and.o(3)==3.and.o(4)==4)order_type = 9
       if(o(1)==4.and.o(2)==3.and.o(3)==1.and.o(4)==2)order_type = 10
       if(o(1)==4.and.o(2)==2.and.o(3)==3.and.o(4)==1)order_type = 11
       if(o(1)==3.and.o(2)==4.and.o(3)==2.and.o(4)==1)order_type = 12
       if(o(1)==2.and.o(2)==4.and.o(3)==1.and.o(4)==3)order_type = 13
       if(o(1)==3.and.o(2)==2.and.o(3)==1.and.o(4)==4)order_type = 14
       if(o(1)==1.and.o(2)==3.and.o(3)==2.and.o(4)==4)order_type = 15
       if(o(1)==4.and.o(2)==1.and.o(3)==3.and.o(4)==2)order_type = 16
       if(o(1)==2.and.o(2)==1.and.o(3)==4.and.o(4)==3)order_type = 17
       if(o(1)==4.and.o(2)==3.and.o(3)==2.and.o(4)==1)order_type = 18
       if(o(1)==2.and.o(2)==4.and.o(3)==3.and.o(4)==1)order_type = 19
       if(o(1)==1.and.o(2)==4.and.o(3)==3.and.o(4)==2)order_type = 20
       if(o(1)==3.and.o(2)==1.and.o(3)==4.and.o(4)==2)order_type = 21
       if(o(1)==3.and.o(2)==2.and.o(3)==4.and.o(4)==1)order_type = 22
       if(o(1)==4.and.o(2)==2.and.o(3)==1.and.o(4)==3)order_type = 23
    case(3)
       if(o(1)==1.and.o(2)==2.and.o(3)==3)            order_type = 0
       if(o(1)==3.and.o(2)==1.and.o(3)==2)            order_type = 24
       if(o(1)==2.and.o(2)==3.and.o(3)==1)            order_type = 25
       if(o(1)==1.and.o(2)==3.and.o(3)==2)            order_type = 26
       if(o(1)==2.and.o(2)==1.and.o(3)==3)            order_type = 27
       if(o(1)==3.and.o(2)==2.and.o(3)==1)            order_type = 28
    case(2)
       if(o(1)==1.and.o(2)==2)                        order_type = 0
       if(o(1)==2.and.o(2)==1)                        order_type = 29
    end select

    call get_midx(tnr,tmodeidx,ntpm,mode)

    ntimes=1
    do i=1,mode
      fels(o(i)) = (tmodeidx(i)-1) * tdims(i) + 1
      if(tmodeidx(i)*tdims(i)>full_dim_tiled_arr(i))then
        acttdim(i)=mod(full_dim_tiled_arr(i),tdims(i))
      else
        acttdim(i)=tdims(i)
      endif
      if(i>1)ntimes=ntimes*acttdim(i)
      rtd(o(i))     = acttdim(i)
    enddo

    select case(order_type)
      case(0)
        ccels=acttdim(1)
        !loop over the remaining not-consecutive dimensions
        do i=1,ntimes
          !get the mode-index in the remaining dimensions
          call get_midx(i,idxintile(2:mode),acttdim(2:mode),mode-1)
          !get the position of the first element in the consecutive stretch
          idxintile(1)=1
          do k=1,mode
            glbmodeidx(k)=idxintile(k) +(tmodeidx(k)-1) *tdims(k)
          enddo
          pos1=get_cidx(glbmodeidx,full_dim_tiled_arr,mode)
          if(pre1==1.0E0_realk.and.pre2==0.0E0_realk)then
            call dcopy(ccels,fort(pos1),1,tileout(1+(i-1)*ccels),1)
          else
            call dscal(ccels,pre2,tileout(1+(i-1)*ccels),1)
            call daxpy(ccels,pre1,fort(pos1),1,tileout(1+(i-1)*ccels),1)
          endif
        enddo
      case(1)
        call manual_3412_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(2)
        call manual_4123_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(3)
        call manual_2341_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(4)
        call manual_1243_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(5)
        call manual_1423_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(6)
        call manual_1342_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(7)
        call manual_3124_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(8)
        call manual_2314_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(9)
        call manual_2134_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(10)
        call manual_4312_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(11)
        call manual_4231_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(12)
        call manual_3421_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(13)
        call manual_2413_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(14)
        call manual_3214_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(15)
        call manual_1324_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(16)
        call manual_4132_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(17)
        call manual_2143_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(18)
        call manual_4321_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(19)
        call manual_2431_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(20)
        call manual_1432_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(21)
        call manual_3142_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(22)
        call manual_3241_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(23)
        call manual_4213_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(24)
        call manual_312_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(25)
        call manual_231_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(26)
        call manual_132_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(27)
        call manual_213_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(28)
        call manual_321_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case(29)
        call manual_21_reordering_f2t(rtd,full_arr_dim,fels,pre1,fort,pre2,tileout)
      case default
        print *,"expensive default tile_from_fort",o
        print *,"order  :",o
        print *,"rorder :",ro
        print *,"atd    :",acttdim
        !count elements in the current tile for loop over elements
        !identify their original position and put them in tile
        nelms=1
        do i=1,mode
          nelms = nelms * acttdim(i)
        enddo
        do i = 1,nelms
          !get mode index of element in tile
          call get_midx(i,idxintile,acttdim,mode)
          !get global index of element, example order = 2 3 1 4 of new array with
          !respect to old --> element 54 3 27 8 of old goes to 3 27 54 8 of new -
          ! old with respect to new 3 1 2 4 
          do k=1,mode
            glbmodeidx(o(k))=idxintile(k) + (tmodeidx(k)-1)*tdims(k)
          enddo

          pos1=get_cidx(glbmodeidx,full_arr_dim,mode)

          if(pre2==0)then
             tileout(i)=pre1*fort(pos1)
          else
             tileout(i)=pre2*tileout(i)+pre1*fort(pos1)
          endif

        enddo
    end select

  end subroutine tile_from_fort
