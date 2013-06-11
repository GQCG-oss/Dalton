!\> \brief this module is inteded to contain high performance reorderings for
!mutlidimensional arrays.
! first part: 4d arrays
! next part: 3d arrays
!\> \author Patrick Ettenhuber & Janus Juul Eriksen
!\> \date November 2012
module manual_reorderings_module
  use precision
  use memory_handling
  use LSTIMING!,only:lstimer
  contains
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
    integer :: dim1,dim2,dim3,dim4,dim1b,dim2b,dim3b,vdim
    integer :: i,j,l
    integer :: aa,bb,cc,dd,block_size,fina,finb,finc,find
    integer :: order_type,m,n
    real(realk) :: tcpu1,twall1,tcpu2,twall2


    call LSTIMER('START',tcpu1,twall1,-1)

    ! assuming available cache memory is 8 MB and we need to store two blocks at a time 
    block_size = int(((8000.0*1000.0)/(8.0*2.0))**(1.0/4.0))
    print *,"block size",block_size, d1,d2,d3,d4
    dims(1)=d1
    dims(2)=d2
    dims(3)=d3
    dims(4)=d4
    dim1 = dims(1)
    dim2 = dims(2)*dims(1)
    dim3 = dims(3)*dims(2)*dims(1)
    dim4 = dims(4)

    ! select  general type of the reordering
    order_type = -1
    if(order(1)==1 .and. order(2)==2 .and. &
         order(3)==3 .and. order(4)==4) order_type = 0
    if(order(1)==3 .and. order(2)==4 .and. &
         order(3)==1 .and. order(4)==2) order_type = 1
    if(order(1)==4 .and. order(2)==1 .and. &
         order(3)==2 .and. order(4)==3) order_type = 1
    !equal to case 9
    if(order(1)==2 .and. order(2)==3 .and. &
         order(3)==4 .and. order(4)==1) order_type = 1
    !equal to case 16
    if(order(1)==2 .and. order(2)==1 .and. &
         order(3)==3 .and. order(4)==4) order_type = 2
    if(order(1)==3 .and. order(2)==1 .and. &
         order(3)==2 .and. order(4)==4) order_type = 2
    if(order(1)==2 .and. order(2)==3 .and. &
         order(3)==1 .and. order(4)==4) order_type = 2
    !if necessary introduce special treatments for slow
    ! reorderings here and put them into manual_reorderings.f90
    if(order(1)==4 .and. order(2)==2 .and. &
         order(3)==3 .and. order(4)==1) order_type = 3
    if(order(1)==2 .and. order(2)==4 .and. &
         order(3)==1 .and. order(4)==3) order_type = 4
    if(order(1)==3 .and. order(2)==2 .and. &
         order(3)==1 .and. order(4)==4) order_type = 5
    if(order(1)==1 .and. order(2)==3 .and. &
         order(3)==4 .and. order(4)==2) order_type = 6
    if(order(1)==1 .and. order(2)==3 .and. &
         order(3)==2 .and. order(4)==4) order_type = 7
    if(order(1)==2 .and. order(2)==3 .and. &
         order(3)==1 .and. order(4)==4) order_type = 8
    ! this already appears in transp case 1!!!
    if(order(1)==2 .and. order(2)==3 .and. &
         order(3)==4 .and. order(4)==1) order_type = 9
    if(order(1)==4 .and. order(2)==1 .and. &
         order(3)==2 .and. order(4)==3) order_type = 10
    if(order(1)==4 .and. order(2)==1 .and. &
         order(3)==3 .and. order(4)==2) order_type = 11
    if(order(1)==1 .and. order(2)==2 .and. &
         order(3)==4 .and. order(4)==3) order_type = 12
    if(order(1)==2 .and. order(2)==1 .and. &
         order(3)==4 .and. order(4)==3) order_type = 13
    if(order(1)==1 .and. order(2)==4 .and. &
         order(3)==2 .and. order(4)==3) order_type = 14
    ! new janus routines
    if(order(1)==4 .and. order(2)==3 .and. &
         order(3)==2 .and. order(4)==1) order_type = 15
    ! this already appears in transp case 2!!!
    if(order(1)==2 .and. order(2)==1 .and. &
         order(3)==3 .and. order(4)==4) order_type = 16
    if(order(1)==2 .and. order(2)==4 .and. &
         order(3)==3 .and. order(4)==1) order_type = 17
    if(order(1)==3 .and. order(2)==4 .and. &
         order(3)==2 .and. order(4)==1) order_type = 18
    if(order(1)==1 .and. order(2)==4 .and. &
         order(3)==3 .and. order(4)==2) order_type = 19
    if(order(1)==3 .and. order(2)==1 .and. &
         order(3)==4 .and. order(4)==2) order_type = 20

    ! do the reordering
    TypeOfReordering: select case(order_type)
    case(0)
       if (pre2 /= 0.0E0_realk) then
         call dscal(d1*d2*d3*d4,pre2,array_out,1)
         call daxpy(d1*d2*d3*d4,pre1,array_in,1,array_out,1)
       else
         call dcopy(d1*d2*d3*d4,array_in,1,array_out,1)
         if (pre1 /= 1.0E0_realk) then
           call dscal(d1*d2*d3*d4,pre1,array_out,1)
         endif
       endif
    case(1)
       if(order(1) == 2) then
         m = dims(1)
         n = dims(3)*dims(4)*dims(2)
       endif
       if(order(1) == 3) then
         m = dims(1)*dims(2)
         n = dims(3)*dims(4)
       endif
       if(order(1) == 4) then
         m = dims(1)*dims(2)*dims(3)
         n = dims(4)
       endif
       !print *,"WARNING(array_reorder4d):case 1 deprecated"
       call mat_transpose_pl(m,n,pre1,array_in,pre2,array_out)
    case(2)
       if(order(1) == 2 .and. order(2) == 1) then
         l = dims(3)
         m = dims(1)
         n = dims(2)
       endif
       if(order(1) == 2 .and. order(2) == 3) then
         l = 1
         m = dims(1)
         n = dims(2)*dims(3)
       endif
       if(order(1) == 3) then
         l = 1
         m = dims(1)*dims(2)
         n = dims(3)
       endif
       !print *,"WARNING(array_reorder4d):case 2 deprecated"
       do j=1,dims(4)*l,1
         call mat_transpose_pl(m,n,pre1,array_in((j-1)*m*n+1), pre2,array_out((j-1)*m*n+1))
       end do
    case(3)
       call manual_4231_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case(4)
       call manual_2413_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case(5)
       call manual_3214_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case(6)
       call manual_1342_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case(7)
       call manual_1324_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case(8)
       call manual_2314_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case(9)
       call manual_2341_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case(10)
       call manual_4123_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case(11)
       call manual_4132_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case(12)
       call manual_1243_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case(13)
       call manual_2143_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case(14)
       call manual_1423_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case(15)
       call manual_4321_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case(16)
       call manual_2134_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case(17)
       call manual_2431_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case(18)
       call manual_3421_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case(19)
       call manual_1432_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case(20)
       call manual_3142_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case default
       print *,'4d_reordering case DEFAULT'
       print *,order
       !loops over blocks
       do dd = 1, dim4, block_size
         find= min(dd + block_size - 1, dims(4))
         do cc = 1, dim3, block_size
           finc= min(cc + block_size - 1, dims(3))
           do bb = 1, dim2, block_size
             finb=min(bb + block_size - 1, dims(2))
             do aa = 1, dim1, block_size
               fina=min(aa + block_size - 1, dims(1))
               !loops over elements
               do d = dd,find
                 do c = cc,finc
                   do b = bb,finb
                     do a = aa,fina 

                       ! old order
                       order1 = [a,b,c,d]
                       ! new order
                       order2 = [order1(order(1)),order1(order(2)), &
                            order1(order(3)),order1(order(4))]
                       ! reorder
                       dim1b=dims(order(1))
                       dim2b=dims(order(1))*dims(order(2))
                       dim3b=dims(order(1))*dims(order(2))*dims(order(3))
                       array_out(order2(1)+(order2(2)-1)*dim1b+(order2(3)-1)*&
                            &dim2b+(order2(4)-1)*dim3b) = &
                            &pre2*array_out(order2(1)+(order2(2)-1)*dim1b+&
                            &(order2(3)-1)*dim2b+(order2(4)-1)*dim3b) + &
                         pre1*array_in(a+(b-1)*dim1+(c-1)*dim2+(d-1)*dim3)
                     end do
                   end do
                 end do
               end do

             end do
           end do
         end do
       end do
    end select TypeOfReordering


    call LSTIMER('START',tcpu2,twall2,-1)
  end subroutine array_reorder_4d

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
    integer :: dim1,dim2,dim3,dim1b,dim2b
    integer :: aa,bb,cc,block_size
    integer :: order_type
    integer :: vec_size
    integer(kind=long) :: vec_size64
    real(realk) :: tcpu1,twall1,tcpu2,twall2

    call LSTIMER('START',tcpu1,twall1,-1)

    ! assuming available cache memory is 8 MB and we need to store two blocks at a time 
    block_size = int(((8000.0*1000.0)/(8.0*2.0))**(1.0/3.0))

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
    dim1 = dims(1)
    dim2 = dims(2)*dims(1)
    dim3 = dims(3)

    ! select  general type of the reordering
    order_type = -1
    if(order(1)==1 .and. order(2)==2 .and. &
         order(3)==3) order_type = 0
    if(order(1)==1 .and. order(2)==3 .and. &
         order(3)==2) order_type = 1
    if(order(1)==2 .and. order(2)==1 .and. &
         order(3)==3) order_type = 2
    if(order(1)==2 .and. order(2)==3 .and. &
         order(3)==1) order_type = 3
    if(order(1)==3 .and. order(2)==1 .and. &
         order(3)==2) order_type = 4
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
       call manual_132_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case(2)
       call manual_213_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case(3)
       call manual_231_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case(4)
       call manual_312_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case(5)
       call manual_321_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case default
    ! ******************************
    ! we should never enter in here!
    ! ******************************
       print *,'3d_reordering case DEFAULT'
       !loops over blocks
       do cc = 1, dim3, block_size
         finc= min(cc + block_size - 1, dims(3))
         do bb = 1, dim2, block_size
           finb=min(bb + block_size - 1, dims(2))
           do aa = 1, dim1, block_size
             fina=min(aa + block_size - 1, dims(1))
             !loops over elements
             do c = cc,finc
               do b = bb,finb
                 do a = aa,fina

                   ! old order
                   order1 = [a,b,c]
                   ! new order
                   order2 = [order1(order(1)),order1(order(2)), &
                        order1(order(3))]
                   ! reorder
                   dim1b=dims(order(1))
                   dim2b=dims(order(1))*dims(order(2))
                   array_out(order2(1)+(order2(2)-1)*dim1b+(order2(3)-1)*dim2b) = &
                        &pre2*array_out(order2(1)+(order2(2)-1)*dim1b+&
                        &(order2(3)-1)*dim2b) + &
                     pre1*array_in(a+(b-1)*dim1+(c-1)*dim2)
                 end do
               end do
             end do

           end do
         end do
       end do

    end select TypeOfReordering

    call LSTIMER('START',tcpu2,twall2,-1)

  end subroutine array_reorder_3d

  !\>  \brief another transposition routine, intended to replace the others
  !\>  \> author Patrick Ettenhuber
  subroutine mat_transpose_pl(r,c,p1,x,p2,y)
    integer,intent(in) :: r, c
    real(realk), intent(in) :: p1,p2
    real(realk), dimension(r,c), intent(in) :: x
    real(realk), dimension(c,r), intent(inout) :: y
    integer :: i,j, s, p, k, a, b,e,f
    integer :: nr, nc, rr, rc, IOK, kb,d, iwrk
    real(realk) :: val
    real(realk), dimension(:)  , pointer :: sm
    integer,     dimension(:)  , pointer :: h
    kb = 512 !cache size of cpu --> /proc/cpuinfo
    s = kb * 1024 / 8
    k = int(0.9*sqrt(float(s)))
    nr=r/k
    nc=c/k
    rr=mod(r,k)
    rc=mod(c,k)
    d = max(rr,rc)
    iwrk = (d+k) / 4
    call mem_alloc(sm,k*k)
    call mem_alloc(h,iwrk)
    if(p2/=1.0E0_realk)call dscal(r*c,p2,y,1)
    !if(p1/=1.0E0_realk)call dscal(r*c,p1,x,1)
    !call dscal(r*c,p2,y,1)
    !call dscal(r*c,p1,x,1)
    do a = 0, nr-1, 1
      do b = 0, nc-1, 1
        do j=1,k,1
          do i=1,k,1
            sm((j-1)*k+i) = x(i+a*k,j+b*k)
          enddo
        enddo
!        call dcopy(k*k,x(a*k,j+b*k),1,sm((j-1)*k+i),1)
        call alg513(sm,k,k,k*k,h,iwrk,IOK)
        call dscal(k*k,p1,sm,1)
!        call daxpy(k*k,p1,sm,1,y(b*k:k+k*b,a*k:k+a*k),1)
        do j=1,k,1
          do i=1,k,1
            y(i+b*k,j+a*k) = y(i+b*k,j+a*k) + sm((j-1)*k+i)
          enddo
        enddo
      enddo
    enddo

    if (rc>0) then
      iwrk = (rc+k)/4
      do a = 0, nr-1, 1
        do j=1,rc,1
          do i=1,k,1
            sm((j-1)*k+i) = x(i+a*k,j+nc*k)
          enddo
        enddo
        call alg513(sm,k,rc,k*rc,h,iwrk,IOK)
        call dscal(k*rc,p1,sm,1)
        do i=1,rc,1
          do j=1,k,1
            y(i+nc*k,j+a*k) = y(i+b*k,j+a*k) + sm((j-1)*rc+i)
          enddo
        enddo
      enddo
    endif

    if (rr>0) then
      iwrk = (rr+k)/4
      do b = 0, nc-1, 1
        do j=1,k,1
          do i=1,rr,1
            sm((j-1)*rr+i) = x(i+nr*k,j+b*k)
          enddo
        enddo
        call alg513(sm,rr,k,k*rr,h,iwrk,IOK)
        call dscal(k*rr,p1,sm,1)
        do j=1,rr,1
          do i=1,k,1
            y(i+b*k,j+nr*k) = y(i+b*k,j+nr*k) + sm((j-1)*k+i)
          enddo
        enddo
      enddo
    endif

    if (rr>0 .and. rc>0) then
      iwrk = (rc+rr)/4
      do j=1,rc,1
        do i=1,rr,1
          sm((j-1)*rr+i) = x(i+nr*k,j+nc*k)
        enddo
      enddo
      call alg513(sm,rr,rc,rc*rr,h,iwrk,IOK)
      call dscal(rc*rr,p1,sm,1)
      do j=1,rr,1
        do i=1,rc,1
          y(i+nc*k,j+nr*k) = y(i+nc*k,j+nr*k) + sm((j-1)*rc+i)
        enddo
      enddo
    endif
    call mem_dealloc(h)
    call mem_dealloc(sm)
  end subroutine mat_transpose_pl

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

  !
  ! *******************
  ! **** 4d arrays ****
  ! *******************
  !
  !\> \brief reorder a 4 diensional array  to get the indices
  !   in the order 3 1 4 2 , this is a quite expensive reordering
  !   and thus requires additional attention
  !\> \autor Janus Juul Eriksen
  !\> \date February 2013
  subroutine manual_3142_reordering(bs,dims,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the original array
    integer, intent(in) :: dims(4)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1,pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3),dims(4))
    !> reordered array
    real(realk),intent(inout) :: array_out(dims(3),dims(1),dims(4),dims(2))
    integer :: a,b,c,d,ba,bb,bc,bd,da,db,dc,dd,da2,db2,dc2,dd2,bcntr
    logical :: moda,modb,modc,modd

    da=dims(1)
    db=dims(2)
    dc=dims(3)
    dd=dims(4)

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs
    dd2=(dd/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)
    modd=(mod(dd,bs)>0)

    bcntr=bs-1

    if (pre2 .eq. 0.0E0_realk .and. pre1 .eq. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
              do bc=1,dc2,bs

                do b=0,bcntr
                  do d=0,bcntr
                    do a=0,bcntr
                      do c=0,bcntr
                        array_out(bc+c,ba+a,bd+d,bb+b)=array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo

              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do bc=1,dc2,bs

                do b=0,bcntr
                  do d=0,bcntr
                    do a=da2+1,da
                      do c=0,bcntr
                        array_out(bc+c,a,bd+d,bb+b)=array_in(a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo

            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bc=1,dc2,bs

              do b=db2+1,db
                do d=0,bcntr
                  do a=0,bcntr
                    do c=0,bcntr
                      array_out(bc+c,ba+a,bd+d,b)= array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo

            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs

              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(c,ba+a,bd+d,bb+b)= array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo

            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs

              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bc+c,ba+a,d,bb+b)= array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo

            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs

            do b=db2+1,db
              do d=0,bcntr
                do a=da2+1,da
                  do c=0,bcntr
                    array_out(bc+c,a,bd+d,b)= array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs

            do d=0,bcntr
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(c,a,bd+d,bb+b)= array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs

            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(bc+c,a,d,bb+b)= array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs

            do b=db2+1,db
              do d=0,bcntr
                do a=0,bcntr
                  do c=dc2+1,dc
                    array_out(c,ba+a,bd+d,b)= array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs

            do b=db2+1,db
              do d=dd2+1,dd
                do a=0,bcntr
                  do c=0,bcntr
                    array_out(bc+c,ba+a,d,b)= array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs

            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(c,ba+a,d,bb+b)= array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(da2>0.and.modb.and.modc.and.modd)then
        do ba=1,da2,bs

          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(c,ba+a,d,b)= array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo

        enddo
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        do bb=1,db2,bs

          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(c,a,d,bb+b)= array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo

        enddo
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        do bc=1,dc2,bs

          do b=db2+1,db
            do d=dd2+1,dd
              do a=da2+1,da
                do c=0,bcntr
                  array_out(bc+c,a,d,b)= array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo

        enddo
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        do bd=1,dd2,bs

          do b=db2+1,db
            do d=0,bcntr
              do a=da2+1,da
                do c=dc2+1,dc
                  array_out(c,a,bd+d,b)= array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo

        enddo
      endif
      if(moda.and.modb.and.modc.and.modd)then
        do b=db2+1,db
          do d=dd2+1,dd
            do a=da2+1,da
              do c=dc2+1,dc
                array_out(c,a,d,b)= array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    else if (pre2 .eq. 0.0E0_realk .and. pre1 .ne. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
              do bc=1,dc2,bs

                do b=0,bcntr
                  do d=0,bcntr
                    do a=0,bcntr
                      do c=0,bcntr
                        array_out(bc+c,ba+a,bd+d,bb+b)=pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo

              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do bc=1,dc2,bs

                do b=0,bcntr
                  do d=0,bcntr
                    do a=da2+1,da
                      do c=0,bcntr
                        array_out(bc+c,a,bd+d,bb+b)=pre1*array_in(a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo

            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bc=1,dc2,bs

              do b=db2+1,db
                do d=0,bcntr
                  do a=0,bcntr
                    do c=0,bcntr
                      array_out(bc+c,ba+a,bd+d,b)= pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo

            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs

              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(c,ba+a,bd+d,bb+b)= pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo

            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs

              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bc+c,ba+a,d,bb+b)= pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo

            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs

            do b=db2+1,db
              do d=0,bcntr
                do a=da2+1,da
                  do c=0,bcntr
                    array_out(bc+c,a,bd+d,b)= pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs

            do d=0,bcntr
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(c,a,bd+d,bb+b)= pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs

            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(bc+c,a,d,bb+b)= pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs

            do b=db2+1,db
              do d=0,bcntr
                do a=0,bcntr
                  do c=dc2+1,dc
                    array_out(c,ba+a,bd+d,b)= pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs

            do b=db2+1,db
              do d=dd2+1,dd
                do a=0,bcntr
                  do c=0,bcntr
                    array_out(bc+c,ba+a,d,b)= pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs

            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(c,ba+a,d,bb+b)= pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(da2>0.and.modb.and.modc.and.modd)then
        do ba=1,da2,bs

          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(c,ba+a,d,b)= pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo

        enddo
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        do bb=1,db2,bs

          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(c,a,d,bb+b)= pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo

        enddo
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        do bc=1,dc2,bs

          do b=db2+1,db
            do d=dd2+1,dd
              do a=da2+1,da
                do c=0,bcntr
                  array_out(bc+c,a,d,b)= pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo

        enddo
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        do bd=1,dd2,bs

          do b=db2+1,db
            do d=0,bcntr
              do a=da2+1,da
                do c=dc2+1,dc
                  array_out(c,a,bd+d,b)= pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo

        enddo
      endif
      if(moda.and.modb.and.modc.and.modd)then
        do b=db2+1,db
          do d=dd2+1,dd
            do a=da2+1,da
              do c=dc2+1,dc
                array_out(c,a,d,b)= pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    else if (pre2 .ne. 0.0E0_realk .and. pre1 .eq. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
              do bc=1,dc2,bs

                do b=0,bcntr
                  do d=0,bcntr
                    do a=0,bcntr
                      do c=0,bcntr
                        array_out(bc+c,ba+a,bd+d,bb+b)=pre2*array_out(bc+c,ba+a,bd+d,bb+b) &
                                                      & + array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo

              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do bc=1,dc2,bs

                do b=0,bcntr
                  do d=0,bcntr
                    do a=da2+1,da
                      do c=0,bcntr
                        array_out(bc+c,a,bd+d,bb+b)=pre2*array_out(bc+c,a,bd+d,bb+b) &
                                                   & + array_in(a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo

            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bc=1,dc2,bs

              do b=db2+1,db
                do d=0,bcntr
                  do a=0,bcntr
                    do c=0,bcntr
                      array_out(bc+c,ba+a,bd+d,b)= pre2*array_out(bc+c,ba+a,bd+d,b) &
                                                 & + array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo

            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs

              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(c,ba+a,bd+d,bb+b)= pre2*array_out(c,ba+a,bd+d,bb+b) &
                                                 & + array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo

            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs

              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bc+c,ba+a,d,bb+b)= pre2*array_out(bc+c,ba+a,d,bb+b) &
                                                 & + array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo

            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs

            do b=db2+1,db
              do d=0,bcntr
                do a=da2+1,da
                  do c=0,bcntr
                    array_out(bc+c,a,bd+d,b)= pre2*array_out(bc+c,a,bd+d,b) &
                                            & + array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs

            do d=0,bcntr
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(c,a,bd+d,bb+b)= pre2*array_out(c,a,bd+d,bb+b) &
                                            & + array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs

            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(bc+c,a,d,bb+b)= pre2*array_out(bc+c,a,d,bb+b) &
                                            & + array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs

            do b=db2+1,db
              do d=0,bcntr
                do a=0,bcntr
                  do c=dc2+1,dc
                    array_out(c,ba+a,bd+d,b)= pre2*array_out(c,ba+a,bd+d,b) &
                                            & + array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs

            do b=db2+1,db
              do d=dd2+1,dd
                do a=0,bcntr
                  do c=0,bcntr
                    array_out(bc+c,ba+a,d,b)= pre2*array_out(bc+c,ba+a,d,b) &
                                            & + array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs

            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(c,ba+a,d,bb+b)= pre2*array_out(c,ba+a,d,bb+b) &
                                            & + array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(da2>0.and.modb.and.modc.and.modd)then
        do ba=1,da2,bs

          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(c,ba+a,d,b)= pre2*array_out(c,ba+a,d,b) &
                                       & + array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo

        enddo
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        do bb=1,db2,bs

          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(c,a,d,bb+b)= pre2*array_out(c,a,d,bb+b) &
                                       & + array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo

        enddo
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        do bc=1,dc2,bs

          do b=db2+1,db
            do d=dd2+1,dd
              do a=da2+1,da
                do c=0,bcntr
                  array_out(bc+c,a,d,b)= pre2*array_out(bc+c,a,d,b) &
                                       & + array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo

        enddo
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        do bd=1,dd2,bs

          do b=db2+1,db
            do d=0,bcntr
              do a=da2+1,da
                do c=dc2+1,dc
                  array_out(c,a,bd+d,b)= pre2*array_out(c,a,bd+d,b) &
                                       & + array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo

        enddo
      endif
      if(moda.and.modb.and.modc.and.modd)then
        do b=db2+1,db
          do d=dd2+1,dd
            do a=da2+1,da
              do c=dc2+1,dc
                array_out(c,a,d,b)= pre2*array_out(c,a,d,b) &
                                  & + array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    else if (pre2 .ne. 0.0E0_realk .and. pre1 .ne. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
              do bc=1,dc2,bs

                do b=0,bcntr
                  do d=0,bcntr
                    do a=0,bcntr
                      do c=0,bcntr
                        array_out(bc+c,ba+a,bd+d,bb+b)=pre2*array_out(bc+c,ba+a,bd+d,bb+b) &
                                                      & + pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo

              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do bc=1,dc2,bs

                do b=0,bcntr
                  do d=0,bcntr
                    do a=da2+1,da
                      do c=0,bcntr
                        array_out(bc+c,a,bd+d,bb+b)=pre2*array_out(bc+c,a,bd+d,bb+b) &
                                                   & + pre1*array_in(a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo

            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bc=1,dc2,bs

              do b=db2+1,db
                do d=0,bcntr
                  do a=0,bcntr
                    do c=0,bcntr
                      array_out(bc+c,ba+a,bd+d,b)= pre2*array_out(bc+c,ba+a,bd+d,b) &
                                                 & + pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo

            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs

              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(c,ba+a,bd+d,bb+b)= pre2*array_out(c,ba+a,bd+d,bb+b) &
                                                 & + pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo

            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs

              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bc+c,ba+a,d,bb+b)= pre2*array_out(bc+c,ba+a,d,bb+b) &
                                                 & + pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo

            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs

            do b=db2+1,db
              do d=0,bcntr
                do a=da2+1,da
                  do c=0,bcntr
                    array_out(bc+c,a,bd+d,b)= pre2*array_out(bc+c,a,bd+d,b) &
                                            & + pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs

            do d=0,bcntr
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(c,a,bd+d,bb+b)= pre2*array_out(c,a,bd+d,bb+b) &
                                            & + pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs

            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(bc+c,a,d,bb+b)= pre2*array_out(bc+c,a,d,bb+b) &
                                            & + pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs

            do b=db2+1,db
              do d=0,bcntr
                do a=0,bcntr
                  do c=dc2+1,dc
                    array_out(c,ba+a,bd+d,b)= pre2*array_out(c,ba+a,bd+d,b) &
                                            & + pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs

            do b=db2+1,db
              do d=dd2+1,dd
                do a=0,bcntr
                  do c=0,bcntr
                    array_out(bc+c,ba+a,d,b)= pre2*array_out(bc+c,ba+a,d,b) &
                                            & + pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs

            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(c,ba+a,d,bb+b)= pre2*array_out(c,ba+a,d,bb+b) &
                                            & + pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo

          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(da2>0.and.modb.and.modc.and.modd)then
        do ba=1,da2,bs

          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(c,ba+a,d,b)= pre2*array_out(c,ba+a,d,b) &
                                       & + pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo

        enddo
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        do bb=1,db2,bs

          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(c,a,d,bb+b)= pre2*array_out(c,a,d,bb+b) &
                                       & + pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo

        enddo
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        do bc=1,dc2,bs

          do b=db2+1,db
            do d=dd2+1,dd
              do a=da2+1,da
                do c=0,bcntr
                  array_out(bc+c,a,d,b)= pre2*array_out(bc+c,a,d,b) &
                                       & + pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo

        enddo
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        do bd=1,dd2,bs

          do b=db2+1,db
            do d=0,bcntr
              do a=da2+1,da
                do c=dc2+1,dc
                  array_out(c,a,bd+d,b)= pre2*array_out(c,a,bd+d,b) &
                                       & + pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo

        enddo
      endif
      if(moda.and.modb.and.modc.and.modd)then
        do b=db2+1,db
          do d=dd2+1,dd
            do a=da2+1,da
              do c=dc2+1,dc
                array_out(c,a,d,b)= pre2*array_out(c,a,d,b) &
                                  & + pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    else
      call lsquit("MANUAL REORDERING:Case not found",-1)
    endif

!    print *,"check 3142"
!    if(pre2==0.0E0_realk)then
!      do b=1,db
!        do d=1,dd
!          do a=1,da
!            do c=1,dc
!              if(array_out(c,a,d,b)/=pre1*array_in(a,b,c,d))then
!                print *,"1432 reordering not correct",a,b,c,d,da,db,dc,dd
!                print *,array_out(c,a,d,b),array_in(a,b,c,d)
!                stop 0
!              endif
!            enddo
!          enddo
!        enddo
!      enddo
!    endif
  end subroutine manual_3142_reordering

  !\> \brief reorder a 4 diensional array  to get the indices
  !   in the order 1 4 3 2 , this is a quite expensive reordering
  !   and thus requires additional attnetion
  !\> \autor Patrick Ettenhuber
  !\> \date Movember 2012
  subroutine manual_1432_reordering(bs,dims,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the original array
    integer, intent(in) :: dims(4)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1,pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3),dims(4))
    !> reordered array
    real(realk),intent(inout) :: array_out(dims(1),dims(4),dims(3),dims(2))
    integer :: a,b,c,d,ba,bb,bc,bd,da,db,dc,dd,da2,db2,dc2,dd2,bcntr
    logical :: moda,modb,modc,modd
    
    da=dims(1)
    db=dims(2)
    dc=dims(3)
    dd=dims(4)

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs
    dd2=(dd/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)
    modd=(mod(dd,bs)>0)

    bcntr=bs-1

    if(pre2==0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
              do ba=1,da2,bs
     
                do b=0,bcntr
                  do c=0,bcntr
                    do d=0,bcntr
                      do a=0,bcntr
                        array_out(ba+a,bd+d,bc+c,bb+b)=array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
     
                do b=0,bcntr
                  do c=0,bcntr
                    do d=0,bcntr
                      do a=da2+1,da
                        array_out(a,bd+d,bc+c,bb+b)=array_in(a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
     
              do b=db2+1,db
                do c=0,bcntr
                  do d=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bd+d,bc+c,b)= array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
     
              do b=0,bcntr
                do c=dc2+1,dc
                  do d=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bd+d,c,bb+b)= array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,d,bc+c,bb+b)= array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do b=db2+1,db
              do c=0,bcntr
                do d=0,bcntr
                  do a=da2+1,da
                    array_out(a,bd+d,bc+c,b)= array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
     
            do b=0,bcntr
              do c=dc2+1,dc
                do d=0,bcntr
                  do a=da2+1,da
                    array_out(a,bd+d,c,bb+b)= array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,d,bc+c,bb+b)= array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do b=db2+1,db
              do c=dc2+1,dc
                do d=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,bd+d,c,b)= array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(ba+a,d,bc+c,b)= array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,d,c,bb+b)= array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(da2>0.and.modb.and.modc.and.modd)then
        do ba=1,da2,bs
     
          do b=db2+1,db
            do c=dc2+1,dc
              do d=dd2+1,dd
                do a=0,bcntr
                  array_out(ba+a,d,c,b)= array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(a,d,c,bb+b)= array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(a,d,bc+c,b)= array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        do bd=1,dd2,bs
     
          do b=db2+1,db
            do c=dc2+1,dc
              do d=0,bcntr
                do a=da2+1,da
                  array_out(a,bd+d,c,b)= array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
      endif
      if(moda.and.modb.and.modc.and.modd)then
        do b=db2+1,db
          do c=dc2+1,dc
            do d=dd2+1,dd
              do a=da2+1,da
                array_out(a,d,c,b)= array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2==0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
              do ba=1,da2,bs
     
                do b=0,bcntr
                  do c=0,bcntr
                    do d=0,bcntr
                      do a=0,bcntr
                        array_out(ba+a,bd+d,bc+c,bb+b)=pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
     
                do b=0,bcntr
                  do c=0,bcntr
                    do d=0,bcntr
                      do a=da2+1,da
                        array_out(a,bd+d,bc+c,bb+b)=pre1*array_in(a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
     
              do b=db2+1,db
                do c=0,bcntr
                  do d=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bd+d,bc+c,b)= pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
     
              do b=0,bcntr
                do c=dc2+1,dc
                  do d=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bd+d,c,bb+b)= pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,d,bc+c,bb+b)= pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do b=db2+1,db
              do c=0,bcntr
                do d=0,bcntr
                  do a=da2+1,da
                    array_out(a,bd+d,bc+c,b)= pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
     
            do b=0,bcntr
              do c=dc2+1,dc
                do d=0,bcntr
                  do a=da2+1,da
                    array_out(a,bd+d,c,bb+b)= pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,d,bc+c,bb+b)= pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do b=db2+1,db
              do c=dc2+1,dc
                do d=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,bd+d,c,b)= pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(ba+a,d,bc+c,b)= pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,d,c,bb+b)= pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(da2>0.and.modb.and.modc.and.modd)then
        do ba=1,da2,bs
     
          do b=db2+1,db
            do c=dc2+1,dc
              do d=dd2+1,dd
                do a=0,bcntr
                  array_out(ba+a,d,c,b)= pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(a,d,c,bb+b)= pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(a,d,bc+c,b)= pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        do bd=1,dd2,bs
     
          do b=db2+1,db
            do c=dc2+1,dc
              do d=0,bcntr
                do a=da2+1,da
                  array_out(a,bd+d,c,b)= pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
      endif
      if(moda.and.modb.and.modc.and.modd)then
        do b=db2+1,db
          do c=dc2+1,dc
            do d=dd2+1,dd
              do a=da2+1,da
                array_out(a,d,c,b)= pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2/=0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
              do ba=1,da2,bs
     
                do b=0,bcntr
                  do c=0,bcntr
                    do d=0,bcntr
                      do a=0,bcntr
                        array_out(ba+a,bd+d,bc+c,bb+b)=pre2*array_out(ba+a,bd+d,bc+c,bb+b)&
                                                     &+array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
     
                do b=0,bcntr
                  do c=0,bcntr
                    do d=0,bcntr
                      do a=da2+1,da
                        array_out(a,bd+d,bc+c,bb+b)=pre2*array_out(a,bd+d,bc+c,bb+b)&
                                                     &+array_in(a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
     
              do b=db2+1,db
                do c=0,bcntr
                  do d=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bd+d,bc+c,b)=pre2*array_out(ba+a,bd+d,bc+c,b)&
                                                   &+array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
     
              do b=0,bcntr
                do c=dc2+1,dc
                  do d=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bd+d,c,bb+b)=pre2*array_out(ba+a,bd+d,c,bb+b)&
                                                   &+array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,d,bc+c,bb+b)=pre2*array_out(ba+a,d,bc+c,bb+b)&
                                                   &+array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do b=db2+1,db
              do c=0,bcntr
                do d=0,bcntr
                  do a=da2+1,da
                    array_out(a,bd+d,bc+c,b)=pre2*array_out(a,bd+d,bc+c,b)&
                                                 &+array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
     
            do b=0,bcntr
              do c=dc2+1,dc
                do d=0,bcntr
                  do a=da2+1,da
                    array_out(a,bd+d,c,bb+b)=pre2*array_out(a,bd+d,c,bb+b)&
                                                 &+array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,d,bc+c,bb+b)=pre2*array_out(a,d,bc+c,bb+b)&
                                                 &+array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do b=db2+1,db
              do c=dc2+1,dc
                do d=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,bd+d,c,b)=pre2*array_out(ba+a,bd+d,c,b)&
                                                 &+array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(ba+a,d,bc+c,b)=pre2*array_out(ba+a,d,bc+c,b)&
                                                 &+array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,d,c,bb+b)=pre2*array_out(ba+a,d,c,bb+b)&
                                                 &+array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(da2>0.and.modb.and.modc.and.modd)then
        do ba=1,da2,bs
     
          do b=db2+1,db
            do c=dc2+1,dc
              do d=dd2+1,dd
                do a=0,bcntr
                  array_out(ba+a,d,c,b)=pre2*array_out(ba+a,d,c,b)&
                                               &+array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(a,d,c,bb+b)=pre2*array_out(a,d,c,bb+b)&
                                               &+array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(a,d,bc+c,b)=pre2*array_out(a,d,bc+c,b)&
                                               &+array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        do bd=1,dd2,bs
     
          do b=db2+1,db
            do c=dc2+1,dc
              do d=0,bcntr
                do a=da2+1,da
                  array_out(a,bd+d,c,b)=pre2*array_out(a,bd+d,c,b)&
                                               &+array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
      endif
      if(moda.and.modb.and.modc.and.modd)then
        do b=db2+1,db
          do c=dc2+1,dc
            do d=dd2+1,dd
              do a=da2+1,da
                array_out(a,d,c,b)=pre2*array_out(a,d,c,b)&
                                             &+array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2/=0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
              do ba=1,da2,bs
     
                do b=0,bcntr
                  do c=0,bcntr
                    do d=0,bcntr
                      do a=0,bcntr
                        array_out(ba+a,bd+d,bc+c,bb+b)=pre2*array_out(ba+a,bd+d,bc+c,bb+b)&
                                                     &+pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
     
                do b=0,bcntr
                  do c=0,bcntr
                    do d=0,bcntr
                      do a=da2+1,da
                        array_out(a,bd+d,bc+c,bb+b)=pre2*array_out(a,bd+d,bc+c,bb+b)&
                                                     &+pre1*array_in(a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
     
              do b=db2+1,db
                do c=0,bcntr
                  do d=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bd+d,bc+c,b)=pre2*array_out(ba+a,bd+d,bc+c,b)&
                                                   &+pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
     
              do b=0,bcntr
                do c=dc2+1,dc
                  do d=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bd+d,c,bb+b)=pre2*array_out(ba+a,bd+d,c,bb+b)&
                                                   &+pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,d,bc+c,bb+b)=pre2*array_out(ba+a,d,bc+c,bb+b)&
                                                   &+pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do b=db2+1,db
              do c=0,bcntr
                do d=0,bcntr
                  do a=da2+1,da
                    array_out(a,bd+d,bc+c,b)=pre2*array_out(a,bd+d,bc+c,b)&
                                                 &+pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
     
            do b=0,bcntr
              do c=dc2+1,dc
                do d=0,bcntr
                  do a=da2+1,da
                    array_out(a,bd+d,c,bb+b)=pre2*array_out(a,bd+d,c,bb+b)&
                                                 &+pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,d,bc+c,bb+b)=pre2*array_out(a,d,bc+c,bb+b)&
                                                 &+pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do b=db2+1,db
              do c=dc2+1,dc
                do d=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,bd+d,c,b)=pre2*array_out(ba+a,bd+d,c,b)&
                                                 &+pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(ba+a,d,bc+c,b)=pre2*array_out(ba+a,d,bc+c,b)&
                                                 &+pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,d,c,bb+b)=pre2*array_out(ba+a,d,c,bb+b)&
                                                 &+pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do b=db2+1,db
            do c=dc2+1,dc
              do d=dd2+1,dd
                do a=0,bcntr
                  array_out(ba+a,d,c,b)=pre2*array_out(ba+a,d,c,b)&
                                               &+pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(a,d,c,bb+b)=pre2*array_out(a,d,c,bb+b)&
                                               &+pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(a,d,bc+c,b)=pre2*array_out(a,d,bc+c,b)&
                                               &+pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do b=db2+1,db
            do c=dc2+1,dc
              do d=0,bcntr
                do a=da2+1,da
                  array_out(a,bd+d,c,b)=pre2*array_out(a,bd+d,c,b)&
                                               &+pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do b=db2+1,db
          do c=dc2+1,dc
            do d=dd2+1,dd
              do a=da2+1,da
                array_out(a,d,c,b)=pre2*array_out(a,d,c,b)&
                                             &+pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    else
      call lsquit("MANUAL REORDERING:Case not found",-1)
    endif

    !print *,"check 1432"
    !if(pre2==0.0E0_realk)then
    !  do b=1,db
    !    do c=1,dc
    !      do d=1,dd
    !        do a=1,da
    !          if(array_out(a,d,c,b)/=pre1*array_in(a,b,c,d))then
    !            print *,"1432 reordering not correct",a,b,c,d,da,db,dc,dd
    !            print *,array_out(a,d,c,b),array_in(a,b,c,d)
    !            stop 0
    !          endif
    !        enddo
    !      enddo
    !    enddo
    !  enddo
    !endif
  end subroutine manual_1432_reordering
  !\> \brief reorder a 4 diensional array  to get the indices
  !   in the order 1 4 2 3 , this is a quite expensive reordering
  !   and thus requires additional attnetion
  !\> \autor Patrick Ettenhuber
  !\> \date Movember 2012
  subroutine manual_1423_reordering(bs,dims,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the original array
    integer, intent(in) :: dims(4)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1,pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3),dims(4))
    !> reordered array
    real(realk),intent(inout) :: array_out(dims(1),dims(4),dims(2),dims(3))
    integer :: a,b,c,d,ba,bb,bc,bd,da,db,dc,dd,da2,db2,dc2,dd2,bcntr
    logical :: moda,modb,modc,modd
    
    da=dims(1)
    db=dims(2)
    dc=dims(3)
    dd=dims(4)

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs
    dd2=(dd/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)
    modd=(mod(dd,bs)>0)

    bcntr=bs-1

    if(pre2==0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do bd=1,dd2,bs
              do ba=1,da2,bs
     
                do c=0,bcntr
                  do b=0,bcntr
                    do d=0,bcntr
                      do a=0,bcntr
                        array_out(ba+a,bd+d,bb+b,bc+c)=array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do bd=1,dd2,bs
     
                do c=0,bcntr
                  do b=0,bcntr
                    do d=0,bcntr
                      do a=da2+1,da
                        array_out(a,bd+d,bb+b,bc+c)= array_in(a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
     
              do c=0,bcntr
                do b=db2+1,db
                  do d=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bd+d,b,bc+c)=array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
     
              do c=dc2+1,dc
                do b=0,bcntr
                  do d=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bd+d,bb+b,c)=array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,d,bb+b,bc+c)=array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do c=0,bcntr
              do b=db2+1,db
                do d=0,bcntr
                  do a=da2+1,da
                    array_out(a,bd+d,b,bc+c)=array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
     
            do c=dc2+1,dc
              do b=0,bcntr
                do d=0,bcntr
                  do a=da2+1,da
                    array_out(a,bd+d,bb+b,c)=array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,d,bb+b,bc+c)=array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do c=dc2+1,dc
              do b=db2+1,db
                do d=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,bd+d,b,c)=array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(ba+a,d,b,bc+c)=array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,d,bb+b,c)=array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do d=dd2+1,dd
                do a=0,bcntr
                  array_out(ba+a,d,b,c)=array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(a,d,bb+b,c)=array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(a,d,b,bc+c)=array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do d=0,bcntr
                do a=da2+1,da
                  array_out(a,bd+d,b,c)=array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do c=dc2+1,dc
          do b=db2+1,db
            do d=dd2+1,dd
              do a=da2+1,da
                array_out(a,d,b,c)=array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2==0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do bd=1,dd2,bs
              do ba=1,da2,bs
     
                do c=0,bcntr
                  do b=0,bcntr
                    do d=0,bcntr
                      do a=0,bcntr
                        array_out(ba+a,bd+d,bb+b,bc+c)=pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do bd=1,dd2,bs
     
                do c=0,bcntr
                  do b=0,bcntr
                    do d=0,bcntr
                      do a=da2+1,da
                        array_out(a,bd+d,bb+b,bc+c)= pre1*array_in(a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
     
              do c=0,bcntr
                do b=db2+1,db
                  do d=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bd+d,b,bc+c)=pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
     
              do c=dc2+1,dc
                do b=0,bcntr
                  do d=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bd+d,bb+b,c)=pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,d,bb+b,bc+c)=pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do c=0,bcntr
              do b=db2+1,db
                do d=0,bcntr
                  do a=da2+1,da
                    array_out(a,bd+d,b,bc+c)=pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
     
            do c=dc2+1,dc
              do b=0,bcntr
                do d=0,bcntr
                  do a=da2+1,da
                    array_out(a,bd+d,bb+b,c)=pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,d,bb+b,bc+c)=pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do c=dc2+1,dc
              do b=db2+1,db
                do d=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,bd+d,b,c)=pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(ba+a,d,b,bc+c)=pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,d,bb+b,c)=pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do d=dd2+1,dd
                do a=0,bcntr
                  array_out(ba+a,d,b,c)=pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(a,d,bb+b,c)=pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(a,d,b,bc+c)=pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do d=0,bcntr
                do a=da2+1,da
                  array_out(a,bd+d,b,c)=pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do c=dc2+1,dc
          do b=db2+1,db
            do d=dd2+1,dd
              do a=da2+1,da
                array_out(a,d,b,c)=pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2/=0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do bd=1,dd2,bs
              do ba=1,da2,bs
     
                do c=0,bcntr
                  do b=0,bcntr
                    do d=0,bcntr
                      do a=0,bcntr
                        array_out(ba+a,bd+d,bb+b,bc+c)=pre2*array_out(ba+a,bd+d,bb+b,bc+c)&
                                                     &+array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do bd=1,dd2,bs
     
                do c=0,bcntr
                  do b=0,bcntr
                    do d=0,bcntr
                      do a=da2+1,da
                        array_out(a,bd+d,bb+b,bc+c)=pre2*array_out(a,bd+d,bb+b,bc+c)&
                                                     &+array_in(a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
     
              do c=0,bcntr
                do b=db2+1,db
                  do d=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bd+d,b,bc+c)=pre2*array_out(ba+a,bd+d,b,bc+c)&
                                                   &+array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
     
              do c=dc2+1,dc
                do b=0,bcntr
                  do d=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bd+d,bb+b,c)=pre2*array_out(ba+a,bd+d,bb+b,c)&
                                                   &+array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,d,bb+b,bc+c)=pre2*array_out(ba+a,d,bb+b,bc+c)&
                                                   &+array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do c=0,bcntr
              do b=db2+1,db
                do d=0,bcntr
                  do a=da2+1,da
                    array_out(a,bd+d,b,bc+c)=pre2*array_out(a,bd+d,b,bc+c)&
                                                 &+array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
     
            do c=dc2+1,dc
              do b=0,bcntr
                do d=0,bcntr
                  do a=da2+1,da
                    array_out(a,bd+d,bb+b,c)=pre2*array_out(a,bd+d,bb+b,c)&
                                                 &+array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,d,bb+b,bc+c)=pre2*array_out(a,d,bb+b,bc+c)&
                                                 &+array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do c=dc2+1,dc
              do b=db2+1,db
                do d=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,bd+d,b,c)=pre2*array_out(ba+a,bd+d,b,c)&
                                                 &+array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(ba+a,d,b,bc+c)=pre2*array_out(ba+a,d,b,bc+c)&
                                                 &+array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,d,bb+b,c)=pre2*array_out(ba+a,d,bb+b,c)&
                                                 &+array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do d=dd2+1,dd
                do a=0,bcntr
                  array_out(ba+a,d,b,c)=pre2*array_out(ba+a,d,b,c)&
                                               &+array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(a,d,bb+b,c)=pre2*array_out(a,d,bb+b,c)&
                                               &+array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(a,d,b,bc+c)=pre2*array_out(a,d,b,bc+c)&
                                               &+array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do d=0,bcntr
                do a=da2+1,da
                  array_out(a,bd+d,b,c)=pre2*array_out(a,bd+d,b,c)&
                                               &+array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do c=dc2+1,dc
          do b=db2+1,db
            do d=dd2+1,dd
              do a=da2+1,da
                array_out(a,d,b,c)=pre2*array_out(a,d,b,c)&
                                             &+array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2/=0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do bd=1,dd2,bs
              do ba=1,da2,bs
     
                do c=0,bcntr
                  do b=0,bcntr
                    do d=0,bcntr
                      do a=0,bcntr
                        array_out(ba+a,bd+d,bb+b,bc+c)=pre2*array_out(ba+a,bd+d,bb+b,bc+c)&
                                                     &+pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do bd=1,dd2,bs
     
                do c=0,bcntr
                  do b=0,bcntr
                    do d=0,bcntr
                      do a=da2+1,da
                        array_out(a,bd+d,bb+b,bc+c)=pre2*array_out(a,bd+d,bb+b,bc+c)&
                                                     &+pre1*array_in(a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
     
              do c=0,bcntr
                do b=db2+1,db
                  do d=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bd+d,b,bc+c)=pre2*array_out(ba+a,bd+d,b,bc+c)&
                                                   &+pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
     
              do c=dc2+1,dc
                do b=0,bcntr
                  do d=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bd+d,bb+b,c)=pre2*array_out(ba+a,bd+d,bb+b,c)&
                                                   &+pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,d,bb+b,bc+c)=pre2*array_out(ba+a,d,bb+b,bc+c)&
                                                   &+pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do c=0,bcntr
              do b=db2+1,db
                do d=0,bcntr
                  do a=da2+1,da
                    array_out(a,bd+d,b,bc+c)=pre2*array_out(a,bd+d,b,bc+c)&
                                                 &+pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
     
            do c=dc2+1,dc
              do b=0,bcntr
                do d=0,bcntr
                  do a=da2+1,da
                    array_out(a,bd+d,bb+b,c)=pre2*array_out(a,bd+d,bb+b,c)&
                                                 &+pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,d,bb+b,bc+c)=pre2*array_out(a,d,bb+b,bc+c)&
                                                 &+pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do c=dc2+1,dc
              do b=db2+1,db
                do d=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,bd+d,b,c)=pre2*array_out(ba+a,bd+d,b,c)&
                                                 &+pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(ba+a,d,b,bc+c)=pre2*array_out(ba+a,d,b,bc+c)&
                                                 &+pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,d,bb+b,c)=pre2*array_out(ba+a,d,bb+b,c)&
                                                 &+pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do d=dd2+1,dd
                do a=0,bcntr
                  array_out(ba+a,d,b,c)=pre2*array_out(ba+a,d,b,c)&
                                               &+pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(a,d,bb+b,c)=pre2*array_out(a,d,bb+b,c)&
                                               &+pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(a,d,b,bc+c)=pre2*array_out(a,d,b,bc+c)&
                                               &+pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do d=0,bcntr
                do a=da2+1,da
                  array_out(a,bd+d,b,c)=pre2*array_out(a,bd+d,b,c)&
                                               &+pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do c=dc2+1,dc
          do b=db2+1,db
            do d=dd2+1,dd
              do a=da2+1,da
                array_out(a,d,b,c)=pre2*array_out(a,d,b,c)&
                                             &+pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    else
      call lsquit("MANUAL REORDERING:Case not found",-1)
    endif

    !print *,"check 1423"
    !if(pre2==0.0E0_realk)then
    !  do c=1,dc
    !    do a=1,da
    !      do d=1,dd
    !        do b=1,db
    !          if(array_out(a,d,b,c)/=pre1*array_in(a,b,c,d))then
    !            print *,"2143 reordering not correct",a,b,c,d,da,db,dc,dd
    !            stop 0
    !          endif
    !        enddo
    !      enddo
    !    enddo
    !  enddo
    !endif
  end subroutine manual_1423_reordering
  !\> \brief reorder a 4 diensional array  to get the indices
  !   in the order 2 1 4 3 , this is a quite expensive reordering
  !   and thus requires additional attnetion
  !\> \autor Patrick Ettenhuber
  !\> \date Movember 2012
  subroutine manual_2143_reordering(bs,dims,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the original array
    integer, intent(in) :: dims(4)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1,pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3),dims(4))
    !> reordered array
    real(realk),intent(inout) :: array_out(dims(2),dims(1),dims(4),dims(3))
    integer :: a,b,c,d,ba,bb,bc,bd,da,db,dc,dd,da2,db2,dc2,dd2,bcntr
    logical :: moda,modb,modc,modd
    
    da=dims(1)
    db=dims(2)
    dc=dims(3)
    dd=dims(4)

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs
    dd2=(dd/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)
    modd=(mod(dd,bs)>0)

    bcntr=bs-1
    if(pre2==0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
              do bb=1,db2,bs
     
                do c=0,bcntr
                  do d=0,bcntr
                    do a=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,ba+a,bd+d,bc+c)=array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
     
              do c=0,bcntr
                do d=0,bcntr
                  do a=da2+1,da
                    do b=0,bcntr
                      array_out(bb+b,a,bd+d,bc+c)=array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
                do d=0,bcntr
                  do c=0,bcntr
                    do b=db2+1,db
                      do a=0,bcntr
                        array_out(b,ba+a,bd+d,bc+c)=array_in(ba+a,b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bb=1,db2,bs
     
              do c=dc2+1,dc
                do d=0,bcntr
                  do a=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,ba+a,bd+d,c)=array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bb+b,ba+a,d,bc+c)=array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do c=0,bcntr
              do d=0,bcntr
                do a=da2+1,da
                  do b=db2+1,db
                    array_out(b,a,bd+d,bc+c)=array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do c=dc2+1,dc
              do d=0,bcntr
                do a=da2+1,da
                  do b=0,bcntr
                    array_out(bb+b,a,bd+d,c)= array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do c=0,bcntr
              do d=dd2+1,dd
                do a=da2+1,da
                  do b=0,bcntr
                    array_out(bb+b,a,d,bc+c)=array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,ba+a,bd+d,c)=array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,ba+a,d,bc+c)=array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
     
            do c=dc2+1,dc
              do d=dd2+1,dd
                do a=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,ba+a,d,c)=array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,ba+a,d,c)=array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do c=dc2+1,dc
            do d=dd2+1,dd
              do a=da2+1,da
                do b=0,bcntr
                  array_out(bb+b,a,d,c)=array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(b,a,d,bc+c)=array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do c=dc2+1,dc
            do d=0,bcntr
              do a=da2+1,da
                do b=db2+1,db
                  array_out(b,a,bd+d,c)=array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do c=dc2+1,dc
          do d=dd2+1,dd
            do a=da2+1,da
              do b=db2+1,db 
                array_out(b,a,d,c)=array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2==0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
              do bb=1,db2,bs
     
                do c=0,bcntr
                  do d=0,bcntr
                    do a=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,ba+a,bd+d,bc+c)=pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
     
              do c=0,bcntr
                do d=0,bcntr
                  do a=da2+1,da
                    do b=0,bcntr
                      array_out(bb+b,a,bd+d,bc+c)=pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
                do d=0,bcntr
                  do c=0,bcntr
                    do b=db2+1,db
                      do a=0,bcntr
                        array_out(b,ba+a,bd+d,bc+c)=pre1*array_in(ba+a,b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bb=1,db2,bs
     
              do c=dc2+1,dc
                do d=0,bcntr
                  do a=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,ba+a,bd+d,c)=pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bb+b,ba+a,d,bc+c)=pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do c=0,bcntr
              do d=0,bcntr
                do a=da2+1,da
                  do b=db2+1,db
                    array_out(b,a,bd+d,bc+c)=pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do c=dc2+1,dc
              do d=0,bcntr
                do a=da2+1,da
                  do b=0,bcntr
                    array_out(bb+b,a,bd+d,c)= pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do c=0,bcntr
              do d=dd2+1,dd
                do a=da2+1,da
                  do b=0,bcntr
                    array_out(bb+b,a,d,bc+c)=pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,ba+a,bd+d,c)=pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,ba+a,d,bc+c)=pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
     
            do c=dc2+1,dc
              do d=dd2+1,dd
                do a=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,ba+a,d,c)=pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,ba+a,d,c)=pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do c=dc2+1,dc
            do d=dd2+1,dd
              do a=da2+1,da
                do b=0,bcntr
                  array_out(bb+b,a,d,c)=pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(b,a,d,bc+c)=pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do c=dc2+1,dc
            do d=0,bcntr
              do a=da2+1,da
                do b=db2+1,db
                  array_out(b,a,bd+d,c)=pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do c=dc2+1,dc
          do d=dd2+1,dd
            do a=da2+1,da
              do b=db2+1,db 
                array_out(b,a,d,c)=pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2/=0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
              do bb=1,db2,bs
     
                do c=0,bcntr
                  do d=0,bcntr
                    do a=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,ba+a,bd+d,bc+c)=pre2*array_out(bb+b,ba+a,bd+d,bc+c)&
                                                     &+array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
     
              do c=0,bcntr
                do d=0,bcntr
                  do a=da2+1,da
                    do b=0,bcntr
                      array_out(bb+b,a,bd+d,bc+c)=pre2*array_out(bb+b,a,bd+d,bc+c)&
                                                   &+array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
                do d=0,bcntr
                  do c=0,bcntr
                    do b=db2+1,db
                      do a=0,bcntr
                        array_out(b,ba+a,bd+d,bc+c)=pre2*array_out(b,ba+a,bd+d,bc+c)&
                                                     &+array_in(ba+a,b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bb=1,db2,bs
     
              do c=dc2+1,dc
                do d=0,bcntr
                  do a=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,ba+a,bd+d,c)=pre2*array_out(bb+b,ba+a,bd+d,c)&
                                                   &+array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bb+b,ba+a,d,bc+c)=pre2*array_out(bb+b,ba+a,d,bc+c)&
                                                   &+array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do c=0,bcntr
              do d=0,bcntr
                do a=da2+1,da
                  do b=db2+1,db
                    array_out(b,a,bd+d,bc+c)=pre2*array_out(b,a,bd+d,bc+c)&
                                                 &+array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do c=dc2+1,dc
              do d=0,bcntr
                do a=da2+1,da
                  do b=0,bcntr
                    array_out(bb+b,a,bd+d,c)=pre2*array_out(bb+b,a,bd+d,c)&
                                                 &+array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do c=0,bcntr
              do d=dd2+1,dd
                do a=da2+1,da
                  do b=0,bcntr
                    array_out(bb+b,a,d,bc+c)=pre2*array_out(bb+b,a,d,bc+c)&
                                                 &+array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,ba+a,bd+d,c)=pre2*array_out(b,ba+a,bd+d,c)&
                                                 &+array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,ba+a,d,bc+c)=pre2*array_out(b,ba+a,d,bc+c)&
                                                 &+array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
     
            do c=dc2+1,dc
              do d=dd2+1,dd
                do a=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,ba+a,d,c)=pre2*array_out(bb+b,ba+a,d,c)&
                                                 &+array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,ba+a,d,c)=pre2*array_out(b,ba+a,d,c)&
                                               &+array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do c=dc2+1,dc
            do d=dd2+1,dd
              do a=da2+1,da
                do b=0,bcntr
                  array_out(bb+b,a,d,c)=pre2*array_out(bb+b,a,d,c)&
                                               &+array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(b,a,d,bc+c)=pre2*array_out(b,a,d,bc+c)&
                                               &+array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do c=dc2+1,dc
            do d=0,bcntr
              do a=da2+1,da
                do b=db2+1,db
                  array_out(b,a,bd+d,c)=pre2*array_out(b,a,bd+d,c)&
                                               &+array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do c=dc2+1,dc
          do d=dd2+1,dd
            do a=da2+1,da
              do b=db2+1,db 
                array_out(b,a,d,c)=pre2*array_out(b,a,d,c)&
                                             &+array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2/=0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
              do bb=1,db2,bs
     
                do c=0,bcntr
                  do d=0,bcntr
                    do a=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,ba+a,bd+d,bc+c)=pre2*array_out(bb+b,ba+a,bd+d,bc+c)&
                                                     &+pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
     
              do c=0,bcntr
                do d=0,bcntr
                  do a=da2+1,da
                    do b=0,bcntr
                      array_out(bb+b,a,bd+d,bc+c)=pre2*array_out(bb+b,a,bd+d,bc+c)&
                                                   &+pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
                do d=0,bcntr
                  do c=0,bcntr
                    do b=db2+1,db
                      do a=0,bcntr
                        array_out(b,ba+a,bd+d,bc+c)=pre2*array_out(b,ba+a,bd+d,bc+c)&
                                                     &+pre1*array_in(ba+a,b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bb=1,db2,bs
     
              do c=dc2+1,dc
                do d=0,bcntr
                  do a=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,ba+a,bd+d,c)=pre2*array_out(bb+b,ba+a,bd+d,c)&
                                                   &+pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bb+b,ba+a,d,bc+c)=pre2*array_out(bb+b,ba+a,d,bc+c)&
                                                   &+pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do c=0,bcntr
              do d=0,bcntr
                do a=da2+1,da
                  do b=db2+1,db
                    array_out(b,a,bd+d,bc+c)=pre2*array_out(b,a,bd+d,bc+c)&
                                                 &+pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do c=dc2+1,dc
              do d=0,bcntr
                do a=da2+1,da
                  do b=0,bcntr
                    array_out(bb+b,a,bd+d,c)=pre2*array_out(bb+b,a,bd+d,c)&
                                                 &+pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do c=0,bcntr
              do d=dd2+1,dd
                do a=da2+1,da
                  do b=0,bcntr
                    array_out(bb+b,a,d,bc+c)=pre2*array_out(bb+b,a,d,bc+c)&
                                                 &+pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,ba+a,bd+d,c)=pre2*array_out(b,ba+a,bd+d,c)&
                                                 &+pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,ba+a,d,bc+c)=pre2*array_out(b,ba+a,d,bc+c)&
                                                 &+pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
     
            do c=dc2+1,dc
              do d=dd2+1,dd
                do a=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,ba+a,d,c)=pre2*array_out(bb+b,ba+a,d,c)&
                                                 &+pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,ba+a,d,c)=pre2*array_out(b,ba+a,d,c)&
                                               &+pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do c=dc2+1,dc
            do d=dd2+1,dd
              do a=da2+1,da
                do b=0,bcntr
                  array_out(bb+b,a,d,c)=pre2*array_out(bb+b,a,d,c)&
                                               &+pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(b,a,d,bc+c)=pre2*array_out(b,a,d,bc+c)&
                                               &+pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do c=dc2+1,dc
            do d=0,bcntr
              do a=da2+1,da
                do b=db2+1,db
                  array_out(b,a,bd+d,c)=pre2*array_out(b,a,bd+d,c)&
                                               &+pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do c=dc2+1,dc
          do d=dd2+1,dd
            do a=da2+1,da
              do b=db2+1,db 
                array_out(b,a,d,c)=pre2*array_out(b,a,d,c)&
                                             &+pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    else
      call lsquit("MANUAL REORDERING:Case not found",-1)
    endif

    !print *,"check 2143"
    !if(pre2==0.0E0_realk)then
    !  do c=1,dc
    !    do a=1,da
    !      do d=1,dd
    !        do b=1,db
    !          if(array_out(b,a,d,c)/=pre1*array_in(a,b,c,d))then
    !            print *,"2143 reordering not correct",a,b,c,d,da,db,dc,dd
    !            stop 0
    !          endif
    !        enddo
    !      enddo
    !    enddo
    !  enddo
    !endif
  end subroutine manual_2143_reordering
  !\> \brief reorder a 4 diensional array  to get the indices
  !   in the order 1 2 4 3 , this is a quite expensive reordering
  !   and thus requires additional attnetion
  !\> \autor Patrick Ettenhuber
  !\> \date Movember 2012
  subroutine manual_1243_reordering(bs,dims,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the original array
    integer, intent(in) :: dims(4)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1,pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3),dims(4))
    !> reordered array
    real(realk),intent(inout) :: array_out(dims(1),dims(2),dims(4),dims(3))
    integer :: a,b,c,d,ba,bb,bc,bd,da,db,dc,dd,da2,db2,dc2,dd2,bcntr
    logical :: moda,modb,modc,modd
    
    da=dims(1)
    db=dims(2)
    dc=dims(3)
    dd=dims(4)

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs
    dd2=(dd/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)
    modd=(mod(dd,bs)>0)

    bcntr=bs-1
    if(pre2==0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
              do ba=1,da2,bs
     
                do c=0,bcntr
                  do d=0,bcntr
                    do b=0,bcntr
                      do a=0,bcntr
                        array_out(ba+a,bb+b,bd+d,bc+c)=array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
     
              do c=0,bcntr
                do d=0,bcntr
                  do b=0,bcntr
                    do a=da2+1,da
                      array_out(a,bb+b,bd+d,bc+c)=array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
     
                do c=0,bcntr
                  do d=0,bcntr
                    do b=db2+1,db
                      do a=0,bcntr
                        array_out(ba+a,b,bd+d,bc+c)=array_in(ba+a,b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do c=dc2+1,dc
                do d=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bb+b,bd+d,c)=array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bb+b,d,bc+c)=array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do c=0,bcntr
              do d=0,bcntr
                do b=db2+1,db
                  do a=da2+1,da
                    array_out(a,b,bd+d,bc+c)=array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do c=dc2+1,dc
              do d=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,bb+b,bd+d,c)=array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,bb+b,d,bc+c)=array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do c=dc2+1,dc
              do d=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(ba+a,b,bd+d,c)=array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(ba+a,b,d,bc+c)=array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do c=dc2+1,dc
              do d=dd2+1,dd
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,bb+b,d,c)=array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do c=dc2+1,dc
            do d=dd2+1,dd
              do b=db2+1,db
                do a=0,bcntr
                  array_out(ba+a,b,d,c)=array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do c=dc2+1,dc
            do d=dd2+1,dd
              do b=0,bcntr
                do a=da2+1,da
                  array_out(a,bb+b,d,c)=array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(a,b,d,bc+c)=array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do c=dc2+1,dc
            do d=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(a,b,bd+d,c)=array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do c=dc2+1,dc
          do d=dd2+1,dd
            do b=db2+1,db 
              do a=da2+1,da
                array_out(a,b,d,c)=array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2==0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
              do ba=1,da2,bs
     
                do c=0,bcntr
                  do d=0,bcntr
                    do b=0,bcntr
                      do a=0,bcntr
                        array_out(ba+a,bb+b,bd+d,bc+c)=pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
     
              do c=0,bcntr
                do d=0,bcntr
                  do b=0,bcntr
                    do a=da2+1,da
                      array_out(a,bb+b,bd+d,bc+c)=pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
     
                do c=0,bcntr
                  do d=0,bcntr
                    do b=db2+1,db
                      do a=0,bcntr
                        array_out(ba+a,b,bd+d,bc+c)=pre1*array_in(ba+a,b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do c=dc2+1,dc
                do d=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bb+b,bd+d,c)=pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bb+b,d,bc+c)=pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do c=0,bcntr
              do d=0,bcntr
                do b=db2+1,db
                  do a=da2+1,da
                    array_out(a,b,bd+d,bc+c)=pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do c=dc2+1,dc
              do d=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,bb+b,bd+d,c)=pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,bb+b,d,bc+c)=pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do c=dc2+1,dc
              do d=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(ba+a,b,bd+d,c)=pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(ba+a,b,d,bc+c)=pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do c=dc2+1,dc
              do d=dd2+1,dd
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,bb+b,d,c)=pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do c=dc2+1,dc
            do d=dd2+1,dd
              do b=db2+1,db
                do a=0,bcntr
                  array_out(ba+a,b,d,c)=pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do c=dc2+1,dc
            do d=dd2+1,dd
              do b=0,bcntr
                do a=da2+1,da
                  array_out(a,bb+b,d,c)=pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(a,b,d,bc+c)=pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do c=dc2+1,dc
            do d=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(a,b,bd+d,c)=pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do c=dc2+1,dc
          do d=dd2+1,dd
            do b=db2+1,db 
              do a=da2+1,da
                array_out(a,b,d,c)=pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2/=0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
              do ba=1,da2,bs
     
                do c=0,bcntr
                  do d=0,bcntr
                    do b=0,bcntr
                      do a=0,bcntr
                        array_out(ba+a,bb+b,bd+d,bc+c)=pre2*array_out(ba+a,bb+b,bd+d,bc+c)&
                                                     &+array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
     
              do c=0,bcntr
                do d=0,bcntr
                  do b=0,bcntr
                    do a=da2+1,da
                      array_out(a,bb+b,bd+d,bc+c)=pre2*array_out(a,bb+b,bd+d,bc+c)&
                                                   &+array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
     
                do c=0,bcntr
                  do d=0,bcntr
                    do b=db2+1,db
                      do a=0,bcntr
                        array_out(ba+a,b,bd+d,bc+c)=pre2*array_out(ba+a,b,bd+d,bc+c)&
                                                     &+array_in(ba+a,b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do c=dc2+1,dc
                do d=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bb+b,bd+d,c)=pre2*array_out(ba+a,bb+b,bd+d,c)&
                                                   &+array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bb+b,d,bc+c)=pre2*array_out(ba+a,bb+b,d,bc+c)&
                                                   &+array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do c=0,bcntr
              do d=0,bcntr
                do b=db2+1,db
                  do a=da2+1,da
                    array_out(a,b,bd+d,bc+c)=pre2*array_out(a,b,bd+d,bc+c)&
                                                 &+array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do c=dc2+1,dc
              do d=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,bb+b,bd+d,c)=pre2*array_out(a,bb+b,bd+d,c)&
                                                 &+array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,bb+b,d,bc+c)=pre2*array_out(a,bb+b,d,bc+c)&
                                                 &+array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do c=dc2+1,dc
              do d=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(ba+a,b,bd+d,c)=pre2*array_out(ba+a,b,bd+d,c)&
                                                 &+array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(ba+a,b,d,bc+c)=pre2*array_out(ba+a,b,d,bc+c)&
                                                 &+array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do c=dc2+1,dc
              do d=dd2+1,dd
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,bb+b,d,c)=pre2*array_out(ba+a,bb+b,d,c)&
                                                 &+array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do c=dc2+1,dc
            do d=dd2+1,dd
              do b=db2+1,db
                do a=0,bcntr
                  array_out(ba+a,b,d,c)=pre2*array_out(ba+a,b,d,c)&
                                               &+array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do c=dc2+1,dc
            do d=dd2+1,dd
              do b=0,bcntr
                do a=da2+1,da
                  array_out(a,bb+b,d,c)=pre2*array_out(a,bb+b,d,c)&
                                               &+array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(a,b,d,bc+c)=pre2*array_out(a,b,d,bc+c)&
                                               &+array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do c=dc2+1,dc
            do d=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(a,b,bd+d,c)=pre2*array_out(a,b,bd+d,c)&
                                               &+array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do c=dc2+1,dc
          do d=dd2+1,dd
            do b=db2+1,db 
              do a=da2+1,da
                array_out(a,b,d,c)=pre2*array_out(a,b,d,c)&
                                             &+array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2/=0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
              do ba=1,da2,bs
     
                do c=0,bcntr
                  do d=0,bcntr
                    do b=0,bcntr
                      do a=0,bcntr
                        array_out(ba+a,bb+b,bd+d,bc+c)=pre2*array_out(ba+a,bb+b,bd+d,bc+c)&
                                                     &+pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
     
              do c=0,bcntr
                do d=0,bcntr
                  do b=0,bcntr
                    do a=da2+1,da
                      array_out(a,bb+b,bd+d,bc+c)=pre2*array_out(a,bb+b,bd+d,bc+c)&
                                                   &+pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do ba=1,da2,bs
     
                do c=0,bcntr
                  do d=0,bcntr
                    do b=db2+1,db
                      do a=0,bcntr
                        array_out(ba+a,b,bd+d,bc+c)=pre2*array_out(ba+a,b,bd+d,bc+c)&
                                                     &+pre1*array_in(ba+a,b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do c=dc2+1,dc
                do d=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bb+b,bd+d,c)=pre2*array_out(ba+a,bb+b,bd+d,c)&
                                                   &+pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bb+b,d,bc+c)=pre2*array_out(ba+a,bb+b,d,bc+c)&
                                                   &+pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do c=0,bcntr
              do d=0,bcntr
                do b=db2+1,db
                  do a=da2+1,da
                    array_out(a,b,bd+d,bc+c)=pre2*array_out(a,b,bd+d,bc+c)&
                                                 &+pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do c=dc2+1,dc
              do d=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,bb+b,bd+d,c)=pre2*array_out(a,bb+b,bd+d,c)&
                                                 &+pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,bb+b,d,bc+c)=pre2*array_out(a,bb+b,d,bc+c)&
                                                 &+pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do c=dc2+1,dc
              do d=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(ba+a,b,bd+d,c)=pre2*array_out(ba+a,b,bd+d,c)&
                                                 &+pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(ba+a,b,d,bc+c)=pre2*array_out(ba+a,b,d,bc+c)&
                                                 &+pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do c=dc2+1,dc
              do d=dd2+1,dd
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,bb+b,d,c)=pre2*array_out(ba+a,bb+b,d,c)&
                                                 &+pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do c=dc2+1,dc
            do d=dd2+1,dd
              do b=db2+1,db
                do a=0,bcntr
                  array_out(ba+a,b,d,c)=pre2*array_out(ba+a,b,d,c)&
                                               &+pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do c=dc2+1,dc
            do d=dd2+1,dd
              do b=0,bcntr
                do a=da2+1,da
                  array_out(a,bb+b,d,c)=pre2*array_out(a,bb+b,d,c)&
                                               &+pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(a,b,d,bc+c)=pre2*array_out(a,b,d,bc+c)&
                                               &+pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do c=dc2+1,dc
            do d=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(a,b,bd+d,c)=pre2*array_out(a,b,bd+d,c)&
                                               &+pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do c=dc2+1,dc
          do d=dd2+1,dd
            do b=db2+1,db 
              do a=da2+1,da
                array_out(a,b,d,c)=pre2*array_out(a,b,d,c)&
                                             &+pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    else
      call lsquit("MANUAL REORDERING:Case not found",-1)
    endif

    !print *,"check 1243"
    !if(pre2==0.0E0_realk)then
    !  do c=1,dc
    !    do a=1,da
    !      do d=1,dd
    !        do b=1,db
    !          if(array_out(a,b,d,c)/=pre1*array_in(a,b,c,d))then
    !            print *,"1243 reordering not correct",a,b,c,d,da,db,dc,dd
    !            stop 0
    !          endif
    !        enddo
    !      enddo
    !    enddo
    !  enddo
    !endif
  end subroutine manual_1243_reordering
  !\> \brief reorder a 4 diensional array  to get the indices
  !   in the order 4 1 3 2 , this is a quite expensive reordering
  !   and thus requires additional attnetion
  !\> \autor Patrick Ettenhuber
  !\> \date Movember 2012
  subroutine manual_4132_reordering(bs,dims,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the original array
    integer, intent(in) :: dims(4)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1,pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3),dims(4))
    !> reordered array
    real(realk),intent(inout) :: array_out(dims(4),dims(1),dims(3),dims(2))
    integer :: a,b,c,d,ba,bb,bc,bd,da,db,dc,dd,da2,db2,dc2,dd2,bcntr
    logical :: moda,modb,modc,modd
    
    da=dims(1)
    db=dims(2)
    dc=dims(3)
    dd=dims(4)

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs
    dd2=(dd/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)
    modd=(mod(dd,bs)>0)

    bcntr=bs-1

    if(pre2==0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
              do bd=1,dd2,bs
     
                do b=0,bcntr
                  do c=0,bcntr
                    do a=0,bcntr
                      do d=0,bcntr
                        array_out(bd+d,ba+a,bc+c,bb+b)=array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
     
              do b=0,bcntr
                do c=0,bcntr
                  do a=da2+1,da
                    do d=0,bcntr
                      array_out(bd+d,a,bc+c,bb+b)=array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
     
              do b=db2+1,db
                do c=0,bcntr
                  do a=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,ba+a,bc+c,b)=array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
            do bd=1,dd2,bs
     
              do b=0,bcntr
                do c=dc2+1,dc
                  do a=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,ba+a,c,bb+b)=array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(d,ba+a,bc+c,bb+b)=array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do b=db2+1,db
              do c=0,bcntr
                do a=da2+1,da
                  do d=0,bcntr
                    array_out(bd+d,a,bc+c,b)=array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
     
            do b=0,bcntr
              do c=dc2+1,dc
                do a=da2+1,da
                  do d=0,bcntr
                    array_out(bd+d,a,c,bb+b)=array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(d,a,bc+c,bb+b)=array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
     
            do b=db2+1,db
              do c=dc2+1,dc
                do a=0,bcntr
                  do d=0,bcntr
                    array_out(bd+d,ba+a,c,b)=array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(d,ba+a,bc+c,b)=array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(d,ba+a,c,bb+b)=array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(d,ba+a,c,b)=array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(d,a,c,bb+b)=array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(d,a,bc+c,b)=array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do b=db2+1,db
            do c=dc2+1,dc
              do a=da2+1,da
                do d=0,bcntr
                  array_out(bd+d,a,c,b)=array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do b=db2+1,db
          do c=dc2+1,dc
            do a=da2+1,da
              do d=dd2+1,dd
                array_out(d,a,c,b)=array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2==0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
              do bd=1,dd2,bs
     
                do b=0,bcntr
                  do c=0,bcntr
                    do a=0,bcntr
                      do d=0,bcntr
                        array_out(bd+d,ba+a,bc+c,bb+b)=pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
     
              do b=0,bcntr
                do c=0,bcntr
                  do a=da2+1,da
                    do d=0,bcntr
                      array_out(bd+d,a,bc+c,bb+b)=pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
     
              do b=db2+1,db
                do c=0,bcntr
                  do a=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,ba+a,bc+c,b)=pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
            do bd=1,dd2,bs
     
              do b=0,bcntr
                do c=dc2+1,dc
                  do a=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,ba+a,c,bb+b)=pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(d,ba+a,bc+c,bb+b)=pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do b=db2+1,db
              do c=0,bcntr
                do a=da2+1,da
                  do d=0,bcntr
                    array_out(bd+d,a,bc+c,b)=pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
     
            do b=0,bcntr
              do c=dc2+1,dc
                do a=da2+1,da
                  do d=0,bcntr
                    array_out(bd+d,a,c,bb+b)=pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(d,a,bc+c,bb+b)=pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
     
            do b=db2+1,db
              do c=dc2+1,dc
                do a=0,bcntr
                  do d=0,bcntr
                    array_out(bd+d,ba+a,c,b)=pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(d,ba+a,bc+c,b)=pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(d,ba+a,c,bb+b)=pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(d,ba+a,c,b)=pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(d,a,c,bb+b)=pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(d,a,bc+c,b)=pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do b=db2+1,db
            do c=dc2+1,dc
              do a=da2+1,da
                do d=0,bcntr
                  array_out(bd+d,a,c,b)=pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do b=db2+1,db
          do c=dc2+1,dc
            do a=da2+1,da
              do d=dd2+1,dd
                array_out(d,a,c,b)=pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2/=0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
              do bd=1,dd2,bs
     
                do b=0,bcntr
                  do c=0,bcntr
                    do a=0,bcntr
                      do d=0,bcntr
                        array_out(bd+d,ba+a,bc+c,bb+b)=pre2*array_out(bd+d,ba+a,bc+c,bb+b)&
                                                     &+array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
     
              do b=0,bcntr
                do c=0,bcntr
                  do a=da2+1,da
                    do d=0,bcntr
                      array_out(bd+d,a,bc+c,bb+b)=pre2*array_out(bd+d,a,bc+c,bb+b)&
                                                   &+array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
     
              do b=db2+1,db
                do c=0,bcntr
                  do a=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,ba+a,bc+c,b)=pre2*array_out(bd+d,ba+a,bc+c,b)&
                                                   &+array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
            do bd=1,dd2,bs
     
              do b=0,bcntr
                do c=dc2+1,dc
                  do a=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,ba+a,c,bb+b)=pre2*array_out(bd+d,ba+a,c,bb+b)&
                                                   &+array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(d,ba+a,bc+c,bb+b)=pre2*array_out(d,ba+a,bc+c,bb+b)&
                                                   &+array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do b=db2+1,db
              do c=0,bcntr
                do a=da2+1,da
                  do d=0,bcntr
                    array_out(bd+d,a,bc+c,b)=pre2*array_out(bd+d,a,bc+c,b)&
                                                 &+array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
     
            do b=0,bcntr
              do c=dc2+1,dc
                do a=da2+1,da
                  do d=0,bcntr
                    array_out(bd+d,a,c,bb+b)=pre2*array_out(bd+d,a,c,bb+b)&
                                                 &+array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(d,a,bc+c,bb+b)=pre2*array_out(d,a,bc+c,bb+b)&
                                                 &+array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
     
            do b=db2+1,db
              do c=dc2+1,dc
                do a=0,bcntr
                  do d=0,bcntr
                    array_out(bd+d,ba+a,c,b)=pre2*array_out(bd+d,ba+a,c,b)&
                                                 &+array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(d,ba+a,bc+c,b)=pre2*array_out(d,ba+a,bc+c,b)&
                                                 &+array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(d,ba+a,c,bb+b)=pre2*array_out(d,ba+a,c,bb+b)&
                                                 &+array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(d,ba+a,c,b)=pre2*array_out(d,ba+a,c,b)&
                                               &+array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(d,a,c,bb+b)=pre2*array_out(d,a,c,bb+b)&
                                               &+array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(d,a,bc+c,b)=pre2*array_out(d,a,bc+c,b)&
                                               &+array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do b=db2+1,db
            do c=dc2+1,dc
              do a=da2+1,da
                do d=0,bcntr
                  array_out(bd+d,a,c,b)=pre2*array_out(bd+d,a,c,b)&
                                               &+array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do b=db2+1,db
          do c=dc2+1,dc
            do a=da2+1,da
              do d=dd2+1,dd
                array_out(d,a,c,b)=pre2*array_out(d,a,c,b)&
                                             &+array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2/=0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
              do bd=1,dd2,bs
     
                do b=0,bcntr
                  do c=0,bcntr
                    do a=0,bcntr
                      do d=0,bcntr
                        array_out(bd+d,ba+a,bc+c,bb+b)=pre2*array_out(bd+d,ba+a,bc+c,bb+b)&
                                                     &+pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
     
              do b=0,bcntr
                do c=0,bcntr
                  do a=da2+1,da
                    do d=0,bcntr
                      array_out(bd+d,a,bc+c,bb+b)=pre2*array_out(bd+d,a,bc+c,bb+b)&
                                                   &+pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
     
              do b=db2+1,db
                do c=0,bcntr
                  do a=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,ba+a,bc+c,b)=pre2*array_out(bd+d,ba+a,bc+c,b)&
                                                   &+pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
            do bd=1,dd2,bs
     
              do b=0,bcntr
                do c=dc2+1,dc
                  do a=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,ba+a,c,bb+b)=pre2*array_out(bd+d,ba+a,c,bb+b)&
                                                   &+pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(d,ba+a,bc+c,bb+b)=pre2*array_out(d,ba+a,bc+c,bb+b)&
                                                   &+pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do b=db2+1,db
              do c=0,bcntr
                do a=da2+1,da
                  do d=0,bcntr
                    array_out(bd+d,a,bc+c,b)=pre2*array_out(bd+d,a,bc+c,b)&
                                                 &+pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
     
            do b=0,bcntr
              do c=dc2+1,dc
                do a=da2+1,da
                  do d=0,bcntr
                    array_out(bd+d,a,c,bb+b)=pre2*array_out(bd+d,a,c,bb+b)&
                                                 &+pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(d,a,bc+c,bb+b)=pre2*array_out(d,a,bc+c,bb+b)&
                                                 &+pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
     
            do b=db2+1,db
              do c=dc2+1,dc
                do a=0,bcntr
                  do d=0,bcntr
                    array_out(bd+d,ba+a,c,b)=pre2*array_out(bd+d,ba+a,c,b)&
                                                 &+pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(d,ba+a,bc+c,b)=pre2*array_out(d,ba+a,bc+c,b)&
                                                 &+pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(d,ba+a,c,bb+b)=pre2*array_out(d,ba+a,c,bb+b)&
                                                 &+pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(d,ba+a,c,b)=pre2*array_out(d,ba+a,c,b)&
                                               &+pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(d,a,c,bb+b)=pre2*array_out(d,a,c,bb+b)&
                                               &+pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(d,a,bc+c,b)=pre2*array_out(d,a,bc+c,b)&
                                               &+pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do b=db2+1,db
            do c=dc2+1,dc
              do a=da2+1,da
                do d=0,bcntr
                  array_out(bd+d,a,c,b)=pre2*array_out(bd+d,a,c,b)&
                                               &+pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do b=db2+1,db
          do c=dc2+1,dc
            do a=da2+1,da
              do d=dd2+1,dd
                array_out(d,a,c,b)=pre2*array_out(d,a,c,b)&
                                             &+pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    else
      call lsquit("MANUAL REORDERING:Case not found",-1)
    endif

    !print *,"check 4132"
    !if(pre2==0.0E0_realk)then
    !  do c=1,dc
    !    do a=1,da
    !      do d=1,dd
    !        do b=1,db
    !          if(array_out(d,a,c,b)/=pre1*array_in(a,b,c,d))then
    !            print *,"4123 reordering not correct",a,b,c,d,da,db,dc,dd
    !            stop 0
    !          endif
    !        enddo
    !      enddo
    !    enddo
    !  enddo
    !endif
  end subroutine manual_4132_reordering
  !\> \brief reorder a 4 diensional array  to get the indices
  !   in the order 4 1 2 3 , this is a quite expensive reordering
  !   and thus requires additional attnetion
  !\> \autor Patrick Ettenhuber
  !\> \date Movember 2012
  subroutine manual_4123_reordering(bs,dims,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the original array
    integer, intent(in) :: dims(4)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1,pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3),dims(4))
    !> reordered array
    real(realk),intent(inout) :: array_out(dims(4),dims(1),dims(2),dims(3))
    integer :: a,b,c,d,ba,bb,bc,bd,da,db,dc,dd,da2,db2,dc2,dd2,bcntr
    logical :: moda,modb,modc,modd
    
    da=dims(1)
    db=dims(2)
    dc=dims(3)
    dd=dims(4)

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs
    dd2=(dd/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)
    modd=(mod(dd,bs)>0)

    bcntr=bs-1
    if(pre2==0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
              do bd=1,dd2,bs
     
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      do d=0,bcntr
                        array_out(bd+d,ba+a,bb+b,bc+c)=array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do bd=1,dd2,bs
     
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    do d=0,bcntr
                      array_out(bd+d,a,bb+b,bc+c)=array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
            do bd=1,dd2,bs
     
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,ba+a,b,bc+c)=array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
            do bd=1,dd2,bs
     
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,ba+a,bb+b,c)=array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(d,ba+a,bb+b,bc+c)=array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  do d=0,bcntr
                    array_out(bd+d,a,b,bc+c)=array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
     
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  do d=0,bcntr
                    array_out(bd+d,a,bb+b,c)=array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(d,a,bb+b,bc+c)=array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
     
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  do d=0,bcntr
                    array_out(bd+d,ba+a,b,c)=array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(d,ba+a,b,bc+c)=array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(d,ba+a,bb+b,c)=array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(d,ba+a,b,c)=array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(d,a,bb+b,c)=array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(d,a,b,bc+c)=array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do a=da2+1,da
                do d=0,bcntr
                  array_out(bd+d,a,b,c)=array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do c=dc2+1,dc
          do b=db2+1,db
            do a=da2+1,da
              do d=dd2+1,dd
                array_out(d,a,b,c)=array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2==0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
              do bd=1,dd2,bs
     
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      do d=0,bcntr
                        array_out(bd+d,ba+a,bb+b,bc+c)=pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do bd=1,dd2,bs
     
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    do d=0,bcntr
                      array_out(bd+d,a,bb+b,bc+c)=pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
            do bd=1,dd2,bs
     
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,ba+a,b,bc+c)=pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
            do bd=1,dd2,bs
     
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,ba+a,bb+b,c)=pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(d,ba+a,bb+b,bc+c)=pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  do d=0,bcntr
                    array_out(bd+d,a,b,bc+c)=pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
     
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  do d=0,bcntr
                    array_out(bd+d,a,bb+b,c)=pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(d,a,bb+b,bc+c)=pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
     
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  do d=0,bcntr
                    array_out(bd+d,ba+a,b,c)=pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(d,ba+a,b,bc+c)=pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(d,ba+a,bb+b,c)=pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(d,ba+a,b,c)=pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(d,a,bb+b,c)=pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(d,a,b,bc+c)=pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do a=da2+1,da
                do d=0,bcntr
                  array_out(bd+d,a,b,c)=pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do c=dc2+1,dc
          do b=db2+1,db
            do a=da2+1,da
              do d=dd2+1,dd
                array_out(d,a,b,c)=pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2/=0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
              do bd=1,dd2,bs
     
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      do d=0,bcntr
                        array_out(bd+d,ba+a,bb+b,bc+c)=pre2*array_out(bd+d,ba+a,bb+b,bc+c)&
                                                     &+array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do bd=1,dd2,bs
     
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    do d=0,bcntr
                      array_out(bd+d,a,bb+b,bc+c)=pre2*array_out(bd+d,a,bb+b,bc+c)&
                                                   &+array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
            do bd=1,dd2,bs
     
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,ba+a,b,bc+c)=pre2*array_out(bd+d,ba+a,b,bc+c)&
                                                   &+array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
            do bd=1,dd2,bs
     
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,ba+a,bb+b,c)=pre2*array_out(bd+d,ba+a,bb+b,c)&
                                                   &+array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(d,ba+a,bb+b,bc+c)=pre2*array_out(d,ba+a,bb+b,bc+c)&
                                                   &+array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  do d=0,bcntr
                    array_out(bd+d,a,b,bc+c)=pre2*array_out(bd+d,a,b,bc+c)&
                                                 &+array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
     
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  do d=0,bcntr
                    array_out(bd+d,a,bb+b,c)=pre2*array_out(bd+d,a,bb+b,c)&
                                                 &+array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(d,a,bb+b,bc+c)=pre2*array_out(d,a,bb+b,bc+c)&
                                                 &+array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
     
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  do d=0,bcntr
                    array_out(bd+d,ba+a,b,c)=pre2*array_out(bd+d,ba+a,b,c)&
                                                 &+array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(d,ba+a,b,bc+c)=pre2*array_out(d,ba+a,b,bc+c)&
                                                 &+array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(d,ba+a,bb+b,c)=pre2*array_out(d,ba+a,bb+b,c)&
                                                 &+array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(d,ba+a,b,c)=pre2*array_out(d,ba+a,b,c)&
                                               &+array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(d,a,bb+b,c)=pre2*array_out(d,a,bb+b,c)&
                                               &+array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(d,a,b,bc+c)=pre2*array_out(d,a,b,bc+c)&
                                               &+array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do a=da2+1,da
                do d=0,bcntr
                  array_out(bd+d,a,b,c)=pre2*array_out(bd+d,a,b,c)&
                                               &+array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do c=dc2+1,dc
          do b=db2+1,db
            do a=da2+1,da
              do d=dd2+1,dd
                array_out(d,a,b,c)=pre2*array_out(d,a,b,c)&
                                             &+array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2/=0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
              do bd=1,dd2,bs
     
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      do d=0,bcntr
                        array_out(bd+d,ba+a,bb+b,bc+c)=pre2*array_out(bd+d,ba+a,bb+b,bc+c)&
                                                     &+pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do bd=1,dd2,bs
     
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    do d=0,bcntr
                      array_out(bd+d,a,bb+b,bc+c)=pre2*array_out(bd+d,a,bb+b,bc+c)&
                                                   &+pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
            do bd=1,dd2,bs
     
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,ba+a,b,bc+c)=pre2*array_out(bd+d,ba+a,b,bc+c)&
                                                   &+pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
            do bd=1,dd2,bs
     
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,ba+a,bb+b,c)=pre2*array_out(bd+d,ba+a,bb+b,c)&
                                                   &+pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(d,ba+a,bb+b,bc+c)=pre2*array_out(d,ba+a,bb+b,bc+c)&
                                                   &+pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  do d=0,bcntr
                    array_out(bd+d,a,b,bc+c)=pre2*array_out(bd+d,a,b,bc+c)&
                                                 &+pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
     
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  do d=0,bcntr
                    array_out(bd+d,a,bb+b,c)=pre2*array_out(bd+d,a,bb+b,c)&
                                                 &+pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(d,a,bb+b,bc+c)=pre2*array_out(d,a,bb+b,bc+c)&
                                                 &+pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
     
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  do d=0,bcntr
                    array_out(bd+d,ba+a,b,c)=pre2*array_out(bd+d,ba+a,b,c)&
                                                 &+pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(d,ba+a,b,bc+c)=pre2*array_out(d,ba+a,b,bc+c)&
                                                 &+pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(d,ba+a,bb+b,c)=pre2*array_out(d,ba+a,bb+b,c)&
                                                 &+pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(d,ba+a,b,c)=pre2*array_out(d,ba+a,b,c)&
                                               &+pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(d,a,bb+b,c)=pre2*array_out(d,a,bb+b,c)&
                                               &+pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(d,a,b,bc+c)=pre2*array_out(d,a,b,bc+c)&
                                               &+pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do a=da2+1,da
                do d=0,bcntr
                  array_out(bd+d,a,b,c)=pre2*array_out(bd+d,a,b,c)&
                                               &+pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do c=dc2+1,dc
          do b=db2+1,db
            do a=da2+1,da
              do d=dd2+1,dd
                array_out(d,a,b,c)=pre2*array_out(d,a,b,c)&
                                             &+pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    else
      call lsquit("MANUAL REORDERING:Case not found",-1)
    endif

    !print *,"check 4123"
    !if(pre2==0.0E0_realk)then
    !  do c=1,dc
    !    do a=1,da
    !      do d=1,dd
    !        do b=1,db
    !          if(array_out(d,a,b,c)/=pre1*array_in(a,b,c,d))then
    !            print *,"4123 reordering not correct",a,b,c,d,da,db,dc,dd
    !            stop 0
    !          endif
    !        enddo
    !      enddo
    !    enddo
    !  enddo
    !endif
  end subroutine manual_4123_reordering

  !\> \brief reorder a 4 diensional array  to get the indices
  !   in the order 2 3 4 1 , this is a quite expensive reordering
  !   and thus requires additional attnetion
  !\> \autor Patrick Ettenhuber
  !\> \date Movember 2012
  subroutine manual_2341_reordering(bs,dims,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the original array
    integer, intent(in) :: dims(4)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1,pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3),dims(4))
    !> reordered array
    real(realk),intent(inout) :: array_out(dims(2),dims(3),dims(4),dims(1))
    integer :: a,b,c,d,ba,bb,bc,bd,da,db,dc,dd,da2,db2,dc2,dd2,bcntr
    logical :: moda,modb,modc,modd
    
    da=dims(1)
    db=dims(2)
    dc=dims(3)
    dd=dims(4)

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs
    dd2=(dd/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)
    modd=(mod(dd,bs)>0)

    bcntr=bs-1


    if(pre2==0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
            do bc=1,dc2,bs
              do bb=1,db2,bs
     
                do a=0,bcntr
                  do d=0,bcntr
                    do c=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,bc+c,bd+d,ba+a)=array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do bb=1,db2,bs
     
              do a=da2+1,da
                do d=0,bcntr
                  do c=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bc+c,bd+d,a)=array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do a=0,bcntr
                      array_out(b,bc+c,bd+d,ba+a)=array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bb+b,c,bd+d,ba+a)=array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bb+b,bc+c,d,ba+a)=array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do a=da2+1,da
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    array_out(b,bc+c,bd+d,a)=array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do a=da2+1,da
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    array_out(bb+b,c,bd+d,a)=array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do a=da2+1,da
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,bc+c,d,a)=array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,c,bd+d,ba+a)=array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,bc+c,d,ba+a)=array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(bb+b,c,d,ba+a)=array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,c,d,ba+a)=array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do a=da2+1,da
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  array_out(bb+b,c,d,a)=array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do a=da2+1,da
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  array_out(b,bc+c,d,a)=array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do a=da2+1,da
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  array_out(b,c,bd+d,a)=array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do a=da2+1,da
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                array_out(b,c,d,a)=array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2==0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
            do bc=1,dc2,bs
              do bb=1,db2,bs
     
                do a=0,bcntr
                  do d=0,bcntr
                    do c=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,bc+c,bd+d,ba+a)=pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do bb=1,db2,bs
     
              do a=da2+1,da
                do d=0,bcntr
                  do c=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bc+c,bd+d,a)=pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do a=0,bcntr
                      array_out(b,bc+c,bd+d,ba+a)=pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bb+b,c,bd+d,ba+a)=pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bb+b,bc+c,d,ba+a)=pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do a=da2+1,da
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    array_out(b,bc+c,bd+d,a)=pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do a=da2+1,da
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    array_out(bb+b,c,bd+d,a)=pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do a=da2+1,da
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,bc+c,d,a)=pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,c,bd+d,ba+a)=pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,bc+c,d,ba+a)=pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(bb+b,c,d,ba+a)=pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,c,d,ba+a)=pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do a=da2+1,da
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  array_out(bb+b,c,d,a)=pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do a=da2+1,da
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  array_out(b,bc+c,d,a)=pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do a=da2+1,da
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  array_out(b,c,bd+d,a)=pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do a=da2+1,da
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                array_out(b,c,d,a)=pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2/=0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
            do bc=1,dc2,bs
              do bb=1,db2,bs
     
                do a=0,bcntr
                  do d=0,bcntr
                    do c=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,bc+c,bd+d,ba+a)=pre2*array_out(bb+b,bc+c,bd+d,ba+a)&
                                                     &+array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do bb=1,db2,bs
     
              do a=da2+1,da
                do d=0,bcntr
                  do c=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bc+c,bd+d,a)=pre2*array_out(bb+b,bc+c,bd+d,a)&
                                                   &+array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do a=0,bcntr
                      array_out(b,bc+c,bd+d,ba+a)=pre2*array_out(b,bc+c,bd+d,ba+a)&
                                                   &+array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bb+b,c,bd+d,ba+a)=pre2*array_out(bb+b,c,bd+d,ba+a)&
                                                   &+array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bb+b,bc+c,d,ba+a)=pre2*array_out(bb+b,bc+c,d,ba+a)&
                                                   &+array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do a=da2+1,da
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    array_out(b,bc+c,bd+d,a)=pre2*array_out(b,bc+c,bd+d,a)&
                                                 &+array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do a=da2+1,da
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    array_out(bb+b,c,bd+d,a)=pre2*array_out(bb+b,c,bd+d,a)&
                                                 &+array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do a=da2+1,da
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,bc+c,d,a)=pre2*array_out(bb+b,bc+c,d,a)&
                                                 &+array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,c,bd+d,ba+a)=pre2*array_out(b,c,bd+d,ba+a)&
                                                 &+array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,bc+c,d,ba+a)=pre2*array_out(b,bc+c,d,ba+a)&
                                                 &+array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(bb+b,c,d,ba+a)=pre2*array_out(bb+b,c,d,ba+a)&
                                                 &+array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,c,d,ba+a)=pre2*array_out(b,c,d,ba+a)&
                                               &+array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do a=da2+1,da
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  array_out(bb+b,c,d,a)=pre2*array_out(bb+b,c,d,a)&
                                               &+array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do a=da2+1,da
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  array_out(b,bc+c,d,a)=pre2*array_out(b,bc+c,d,a)&
                                               &+array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do a=da2+1,da
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  array_out(b,c,bd+d,a)=pre2*array_out(b,c,bd+d,a)&
                                               &+array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do a=da2+1,da
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                array_out(b,c,d,a)=pre2*array_out(b,c,d,a)&
                                             &+array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2/=0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
            do bc=1,dc2,bs
              do bb=1,db2,bs
     
                do a=0,bcntr
                  do d=0,bcntr
                    do c=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,bc+c,bd+d,ba+a)=pre2*array_out(bb+b,bc+c,bd+d,ba+a)&
                                                     &+pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do bb=1,db2,bs
     
              do a=da2+1,da
                do d=0,bcntr
                  do c=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bc+c,bd+d,a)=pre2*array_out(bb+b,bc+c,bd+d,a)&
                                                   &+pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do a=0,bcntr
                      array_out(b,bc+c,bd+d,ba+a)=pre2*array_out(b,bc+c,bd+d,ba+a)&
                                                   &+pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bb+b,c,bd+d,ba+a)=pre2*array_out(bb+b,c,bd+d,ba+a)&
                                                   &+pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bb+b,bc+c,d,ba+a)=pre2*array_out(bb+b,bc+c,d,ba+a)&
                                                   &+pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do a=da2+1,da
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    array_out(b,bc+c,bd+d,a)=pre2*array_out(b,bc+c,bd+d,a)&
                                                 &+pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do a=da2+1,da
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    array_out(bb+b,c,bd+d,a)=pre2*array_out(bb+b,c,bd+d,a)&
                                                 &+pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do a=da2+1,da
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,bc+c,d,a)=pre2*array_out(bb+b,bc+c,d,a)&
                                                 &+pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,c,bd+d,ba+a)=pre2*array_out(b,c,bd+d,ba+a)&
                                                 &+pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,bc+c,d,ba+a)=pre2*array_out(b,bc+c,d,ba+a)&
                                                 &+pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(bb+b,c,d,ba+a)=pre2*array_out(bb+b,c,d,ba+a)&
                                                 &+pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,c,d,ba+a)=pre2*array_out(b,c,d,ba+a)&
                                               &+pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do a=da2+1,da
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  array_out(bb+b,c,d,a)=pre2*array_out(bb+b,c,d,a)&
                                               &+pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do a=da2+1,da
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  array_out(b,bc+c,d,a)=pre2*array_out(b,bc+c,d,a)&
                                               &+pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do a=da2+1,da
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  array_out(b,c,bd+d,a)=pre2*array_out(b,c,bd+d,a)&
                                               &+pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do a=da2+1,da
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                array_out(b,c,d,a)=pre2*array_out(b,c,d,a)&
                                             &+pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    else
      call lsquit("MANUAL REORDERING:Case not found",-1)
    endif

    !print *,"check 2341"
    !if(pre2==0.0E0_realk)then
    !  do c=1,dc
    !    do a=1,da
    !      do d=1,dd
    !        do b=1,db
    !          if(array_out(b,c,d,a)/=pre1*array_in(a,b,c,d))then
    !            print *,"2341 reordering not correct",a,b,c,d,da,db,dc,dd
    !            stop 0
    !          endif
    !        enddo
    !      enddo
    !    enddo
    !  enddo
    !endif
  end subroutine manual_2341_reordering
  !\> \brief reorder a 4 diensional array  to get the indices
  !   in the order 2 3 1 4 , this is a quite expensive reordering
  !   and thus requires additional attnetion
  !\> \autor Patrick Ettenhuber
  !\> \date Movember 2012
  subroutine manual_2314_reordering(bs,dims,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the original array
    integer, intent(in) :: dims(4)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1,pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3),dims(4))
    !> reordered array
    real(realk),intent(inout) :: array_out(dims(2),dims(3),dims(1),dims(4))
    integer :: a,b,c,d,ba,bb,bc,bd,da,db,dc,dd,da2,db2,dc2,dd2,bcntr
    logical :: moda,modb,modc,modd
    
    da=dims(1)
    db=dims(2)
    dc=dims(3)
    dd=dims(4)

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs
    dd2=(dd/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)
    modd=(mod(dd,bs)>0)

    bcntr=bs-1


    if(pre2==0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bc=1,dc2,bs
              do bb=1,db2,bs
     
                do d=0,bcntr
                  do a=0,bcntr
                    do c=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,bc+c,ba+a,bd+d)=array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do bb=1,db2,bs
     
              do d=0,bcntr
                do a=da2+1,da
                  do c=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bc+c,a,bd+d)=array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do a=0,bcntr
                      array_out(b,bc+c,ba+a,bd+d)=array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bb+b,c,ba+a,bd+d)=array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bb=1,db2,bs
     
              do d=dd2+1,dd
                do a=0,bcntr
                  do c=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bc+c,ba+a,d)=array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do d=0,bcntr
              do a=da2+1,da
                do c=0,bcntr
                  do b=db2+1,db
                    array_out(b,bc+c,a,bd+d)=array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do d=0,bcntr
              do a=da2+1,da
                do c=dc2+1,dc
                  do b=0,bcntr
                    array_out(bb+b,c,a,bd+d)=array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do a=da2+1,da
                do c=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,bc+c,a,d)=array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,c,ba+a,bd+d)=array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,bc+c,ba+a,d)=array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(bb+b,c,ba+a,d)=array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,c,ba+a,d)=array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do a=da2+1,da
              do c=dc2+1,dc
                do b=0,bcntr
                  array_out(bb+b,c,a,d)=array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do a=da2+1,da
              do c=0,bcntr
                do b=db2+1,db
                  array_out(b,bc+c,a,d)=array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do d=0,bcntr
            do a=da2+1,da
              do c=dc2+1,dc
                do b=db2+1,db
                  array_out(b,c,a,bd+d)=array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do d=dd2+1,dd
          do a=da2+1,da
            do c=dc2+1,dc
              do b=db2+1,db
                array_out(b,c,a,d)=array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2==0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bc=1,dc2,bs
              do bb=1,db2,bs
     
                do d=0,bcntr
                  do a=0,bcntr
                    do c=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,bc+c,ba+a,bd+d)=pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do bb=1,db2,bs
     
              do d=0,bcntr
                do a=da2+1,da
                  do c=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bc+c,a,bd+d)=pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do a=0,bcntr
                      array_out(b,bc+c,ba+a,bd+d)=pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bb+b,c,ba+a,bd+d)=pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bb=1,db2,bs
     
              do d=dd2+1,dd
                do a=0,bcntr
                  do c=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bc+c,ba+a,d)=pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do d=0,bcntr
              do a=da2+1,da
                do c=0,bcntr
                  do b=db2+1,db
                    array_out(b,bc+c,a,bd+d)=pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do d=0,bcntr
              do a=da2+1,da
                do c=dc2+1,dc
                  do b=0,bcntr
                    array_out(bb+b,c,a,bd+d)=pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do a=da2+1,da
                do c=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,bc+c,a,d)=pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,c,ba+a,bd+d)=pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,bc+c,ba+a,d)=pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(bb+b,c,ba+a,d)=pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,c,ba+a,d)=pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do a=da2+1,da
              do c=dc2+1,dc
                do b=0,bcntr
                  array_out(bb+b,c,a,d)=pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do a=da2+1,da
              do c=0,bcntr
                do b=db2+1,db
                  array_out(b,bc+c,a,d)=pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do d=0,bcntr
            do a=da2+1,da
              do c=dc2+1,dc
                do b=db2+1,db
                  array_out(b,c,a,bd+d)=pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do d=dd2+1,dd
          do a=da2+1,da
            do c=dc2+1,dc
              do b=db2+1,db
                array_out(b,c,a,d)=pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2/=0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bc=1,dc2,bs
              do bb=1,db2,bs
     
                do d=0,bcntr
                  do a=0,bcntr
                    do c=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,bc+c,ba+a,bd+d)=pre2*array_out(bb+b,bc+c,ba+a,bd+d)&
                                                     &+array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do bb=1,db2,bs
     
              do d=0,bcntr
                do a=da2+1,da
                  do c=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bc+c,a,bd+d)=pre2*array_out(bb+b,bc+c,a,bd+d)&
                                                   &+array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do a=0,bcntr
                      array_out(b,bc+c,ba+a,bd+d)=pre2*array_out(b,bc+c,ba+a,bd+d)&
                                                   &+array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bb+b,c,ba+a,bd+d)=pre2*array_out(bb+b,c,ba+a,bd+d)&
                                                   &+array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bb=1,db2,bs
     
              do d=dd2+1,dd
                do a=0,bcntr
                  do c=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bc+c,ba+a,d)=pre2*array_out(bb+b,bc+c,ba+a,d)&
                                                   &+array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do d=0,bcntr
              do a=da2+1,da
                do c=0,bcntr
                  do b=db2+1,db
                    array_out(b,bc+c,a,bd+d)=pre2*array_out(b,bc+c,a,bd+d)&
                                                 &+array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do d=0,bcntr
              do a=da2+1,da
                do c=dc2+1,dc
                  do b=0,bcntr
                    array_out(bb+b,c,a,bd+d)=pre2*array_out(bb+b,c,a,bd+d)&
                                                 &+array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do a=da2+1,da
                do c=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,bc+c,a,d)=pre2*array_out(bb+b,bc+c,a,d)&
                                                 &+array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,c,ba+a,bd+d)=pre2*array_out(b,c,ba+a,bd+d)&
                                                 &+array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,bc+c,ba+a,d)=pre2*array_out(b,bc+c,ba+a,d)&
                                                 &+array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(bb+b,c,ba+a,d)=pre2*array_out(bb+b,c,ba+a,d)&
                                                 &+array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,c,ba+a,d)=pre2*array_out(b,c,ba+a,d)&
                                               &+array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do a=da2+1,da
              do c=dc2+1,dc
                do b=0,bcntr
                  array_out(bb+b,c,a,d)=pre2*array_out(bb+b,c,a,d)&
                                               &+array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do a=da2+1,da
              do c=0,bcntr
                do b=db2+1,db
                  array_out(b,bc+c,a,d)=pre2*array_out(b,bc+c,a,d)&
                                               &+array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do d=0,bcntr
            do a=da2+1,da
              do c=dc2+1,dc
                do b=db2+1,db
                  array_out(b,c,a,bd+d)=pre2*array_out(b,c,a,bd+d)&
                                               &+array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do d=dd2+1,dd
          do a=da2+1,da
            do c=dc2+1,dc
              do b=db2+1,db
                array_out(b,c,a,d)=pre2*array_out(b,c,a,d)&
                                             &+array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2/=0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bc=1,dc2,bs
              do bb=1,db2,bs
     
                do d=0,bcntr
                  do a=0,bcntr
                    do c=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,bc+c,ba+a,bd+d)=pre2*array_out(bb+b,bc+c,ba+a,bd+d)&
                                                     &+pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do bb=1,db2,bs
     
              do d=0,bcntr
                do a=da2+1,da
                  do c=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bc+c,a,bd+d)=pre2*array_out(bb+b,bc+c,a,bd+d)&
                                                   &+pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do a=0,bcntr
                      array_out(b,bc+c,ba+a,bd+d)=pre2*array_out(b,bc+c,ba+a,bd+d)&
                                                   &+pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bb+b,c,ba+a,bd+d)=pre2*array_out(bb+b,c,ba+a,bd+d)&
                                                   &+pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bb=1,db2,bs
     
              do d=dd2+1,dd
                do a=0,bcntr
                  do c=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bc+c,ba+a,d)=pre2*array_out(bb+b,bc+c,ba+a,d)&
                                                   &+pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do d=0,bcntr
              do a=da2+1,da
                do c=0,bcntr
                  do b=db2+1,db
                    array_out(b,bc+c,a,bd+d)=pre2*array_out(b,bc+c,a,bd+d)&
                                                 &+pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do d=0,bcntr
              do a=da2+1,da
                do c=dc2+1,dc
                  do b=0,bcntr
                    array_out(bb+b,c,a,bd+d)=pre2*array_out(bb+b,c,a,bd+d)&
                                                 &+pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do a=da2+1,da
                do c=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,bc+c,a,d)=pre2*array_out(bb+b,bc+c,a,d)&
                                                 &+pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,c,ba+a,bd+d)=pre2*array_out(b,c,ba+a,bd+d)&
                                                 &+pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,bc+c,ba+a,d)=pre2*array_out(b,bc+c,ba+a,d)&
                                                 &+pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(bb+b,c,ba+a,d)=pre2*array_out(bb+b,c,ba+a,d)&
                                                 &+pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,c,ba+a,d)=pre2*array_out(b,c,ba+a,d)&
                                               &+pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do a=da2+1,da
              do c=dc2+1,dc
                do b=0,bcntr
                  array_out(bb+b,c,a,d)=pre2*array_out(bb+b,c,a,d)&
                                               &+pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do a=da2+1,da
              do c=0,bcntr
                do b=db2+1,db
                  array_out(b,bc+c,a,d)=pre2*array_out(b,bc+c,a,d)&
                                               &+pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do d=0,bcntr
            do a=da2+1,da
              do c=dc2+1,dc
                do b=db2+1,db
                  array_out(b,c,a,bd+d)=pre2*array_out(b,c,a,bd+d)&
                                               &+pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do d=dd2+1,dd
          do a=da2+1,da
            do c=dc2+1,dc
              do b=db2+1,db
                array_out(b,c,a,d)=pre2*array_out(b,c,a,d)&
                                             &+pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    else
      call lsquit("MANUAL REORDERING:Case not found",-1)
    endif

    !print *,"check 2314"
    !if(pre2==0.0E0_realk)then
    !  do c=1,dc
    !    do a=1,da
    !      do d=1,dd
    !        do b=1,db
    !          if(array_out(b,c,a,d)/=pre1*array_in(a,b,c,d))then
    !            print *,"2314 reordering not correct",a,b,c,d,da,db,dc,dd
    !            stop 0
    !          endif
    !        enddo
    !      enddo
    !    enddo
    !  enddo
    !endif
  end subroutine manual_2314_reordering
  !\> \brief reorder a 4 diensional array  to get the indices
  !   in the order 1 3 2 4 , this is a quite expensive reordering
  !   and thus requires additional attnetion
  !\> \autor Patrick Ettenhuber
  !\> \date Movember 2012
  subroutine manual_1324_reordering(bs,dims,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the original array
    integer, intent(in) :: dims(4)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1,pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3),dims(4))
    !> reordered array
    real(realk),intent(inout) :: array_out(dims(1),dims(3),dims(2),dims(4))
    integer :: a,b,c,d,ba,bb,bc,bd,da,db,dc,dd,da2,db2,dc2,dd2,bcntr
    logical :: moda,modb,modc,modd
    
    da=dims(1)
    db=dims(2)
    dc=dims(3)
    dd=dims(4)

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs
    dd2=(dd/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)
    modd=(mod(dd,bs)>0)

    bcntr=bs-1
    if(pre2==0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do bc=1,dc2,bs
              do ba=1,da2,bs
     
                do d=0,bcntr
                  do b=0,bcntr
                    do c=0,bcntr
                      do a=0,bcntr
                        array_out(ba+a,bc+c,bb+b,bd+d)=array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do bc=1,dc2,bs
     
                do d=0,bcntr
                  do b=0,bcntr
                    do c=0,bcntr
                      do a=da2+1,da
                        array_out(a,bc+c,bb+b,bd+d)=array_in(a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do b=db2+1,db
                  do c=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bc+c,b,bd+d)=array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,c,bb+b,bd+d)=array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bc+c,bb+b,d)=array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do d=0,bcntr
              do b=db2+1,db
                do c=0,bcntr
                  do a=da2+1,da
                    array_out(a,bc+c,b,bd+d)=array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,c,bb+b,bd+d)=array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,bc+c,bb+b,d)=array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do b=db2+1,db
                do c=dc2+1,dc
                  do a=0,bcntr
                    array_out(ba+a,c,b,bd+d)=array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do b=db2+1,db
                do c=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,bc+c,b,d)=array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,c,bb+b,d)=array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do b=db2+1,db
              do c=dc2+1,dc
                do a=0,bcntr
                  array_out(ba+a,c,b,d)=array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(a,c,bb+b,d)=array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do b=db2+1,db
              do c=0,bcntr
                do a=da2+1,da
                  array_out(a,bc+c,b,d)=array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do d=0,bcntr
            do b=db2+1,db
              do c=dc2+1,dc
                do a=da2+1,da
                  array_out(a,c,b,bd+d)=array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do b=db2+1,db
          do d=dd2+1,dd
            do c=dc2+1,dc
              do a=da2+1,da
                array_out(a,c,b,d)=array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2==0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do bc=1,dc2,bs
              do ba=1,da2,bs
     
                do d=0,bcntr
                  do b=0,bcntr
                    do c=0,bcntr
                      do a=0,bcntr
                        array_out(ba+a,bc+c,bb+b,bd+d)=pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do bc=1,dc2,bs
     
                do d=0,bcntr
                  do b=0,bcntr
                    do c=0,bcntr
                      do a=da2+1,da
                        array_out(a,bc+c,bb+b,bd+d)=pre1*array_in(a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do b=db2+1,db
                  do c=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bc+c,b,bd+d)=pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,c,bb+b,bd+d)=pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bc+c,bb+b,d)=pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do d=0,bcntr
              do b=db2+1,db
                do c=0,bcntr
                  do a=da2+1,da
                    array_out(a,bc+c,b,bd+d)=pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,c,bb+b,bd+d)=pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,bc+c,bb+b,d)=pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do b=db2+1,db
                do c=dc2+1,dc
                  do a=0,bcntr
                    array_out(ba+a,c,b,bd+d)=pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do b=db2+1,db
                do c=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,bc+c,b,d)=pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,c,bb+b,d)=pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do b=db2+1,db
              do c=dc2+1,dc
                do a=0,bcntr
                  array_out(ba+a,c,b,d)=pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(a,c,bb+b,d)=pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do b=db2+1,db
              do c=0,bcntr
                do a=da2+1,da
                  array_out(a,bc+c,b,d)=pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do d=0,bcntr
            do b=db2+1,db
              do c=dc2+1,dc
                do a=da2+1,da
                  array_out(a,c,b,bd+d)=pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do b=db2+1,db
          do d=dd2+1,dd
            do c=dc2+1,dc
              do a=da2+1,da
                array_out(a,c,b,d)=pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2/=0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do bc=1,dc2,bs
              do ba=1,da2,bs
     
                do d=0,bcntr
                  do b=0,bcntr
                    do c=0,bcntr
                      do a=0,bcntr
                        array_out(ba+a,bc+c,bb+b,bd+d)=pre2*array_out(ba+a,bc+c,bb+b,bd+d)&
                                                     &+array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do bc=1,dc2,bs
     
                do d=0,bcntr
                  do b=0,bcntr
                    do c=0,bcntr
                      do a=da2+1,da
                        array_out(a,bc+c,bb+b,bd+d)=pre2*array_out(a,bc+c,bb+b,bd+d)&
                                                  &+array_in(a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do b=db2+1,db
                  do c=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bc+c,b,bd+d)=pre2*array_out(ba+a,bc+c,b,bd+d)&
                                                   &+array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,c,bb+b,bd+d)=pre2*array_out(ba+a,c,bb+b,bd+d)&
                                                   &+array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bc+c,bb+b,d)=pre2*array_out(ba+a,bc+c,bb+b,d)&
                                                   &+array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do d=0,bcntr
              do b=db2+1,db
                do c=0,bcntr
                  do a=da2+1,da
                    array_out(a,bc+c,b,bd+d)=pre2*array_out(a,bc+c,b,bd+d)&
                                                 &+array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,c,bb+b,bd+d)=pre2*array_out(a,c,bb+b,bd+d)&
                                                 &+array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,bc+c,bb+b,d)=pre2*array_out(a,bc+c,bb+b,d)&
                                                 &+array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do b=db2+1,db
                do c=dc2+1,dc
                  do a=0,bcntr
                    array_out(ba+a,c,b,bd+d)=pre2*array_out(ba+a,c,b,bd+d)&
                                                 &+array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do b=db2+1,db
                do c=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,bc+c,b,d)=pre2*array_out(ba+a,bc+c,b,d)&
                                                 &+array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,c,bb+b,d)=pre2*array_out(ba+a,c,bb+b,d)&
                                                 &+array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do b=db2+1,db
              do c=dc2+1,dc
                do a=0,bcntr
                  array_out(ba+a,c,b,d)=pre2*array_out(ba+a,c,b,d)&
                                               &+array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(a,c,bb+b,d)=pre2*array_out(a,c,bb+b,d)&
                                               &+array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do b=db2+1,db
              do c=0,bcntr
                do a=da2+1,da
                  array_out(a,bc+c,b,d)=pre2*array_out(a,bc+c,b,d)&
                                               &+array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do d=0,bcntr
            do b=db2+1,db
              do c=dc2+1,dc
                do a=da2+1,da
                  array_out(a,c,b,bd+d)=pre2*array_out(a,c,b,bd+d)&
                                               &+array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do b=db2+1,db
          do d=dd2+1,dd
            do c=dc2+1,dc
              do a=da2+1,da
                array_out(a,c,b,d)=pre2*array_out(a,c,b,d)&
                                             &+array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2/=0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do bc=1,dc2,bs
              do ba=1,da2,bs
     
                do d=0,bcntr
                  do b=0,bcntr
                    do c=0,bcntr
                      do a=0,bcntr
                        array_out(ba+a,bc+c,bb+b,bd+d)=pre2*array_out(ba+a,bc+c,bb+b,bd+d)&
                                                     &+pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do bc=1,dc2,bs
     
                do d=0,bcntr
                  do b=0,bcntr
                    do c=0,bcntr
                      do a=da2+1,da
                        array_out(a,bc+c,bb+b,bd+d)=pre2*array_out(a,bc+c,bb+b,bd+d)&
                                                  &+pre1*array_in(a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do b=db2+1,db
                  do c=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bc+c,b,bd+d)=pre2*array_out(ba+a,bc+c,b,bd+d)&
                                                   &+pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,c,bb+b,bd+d)=pre2*array_out(ba+a,c,bb+b,bd+d)&
                                                   &+pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bc+c,bb+b,d)=pre2*array_out(ba+a,bc+c,bb+b,d)&
                                                   &+pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do d=0,bcntr
              do b=db2+1,db
                do c=0,bcntr
                  do a=da2+1,da
                    array_out(a,bc+c,b,bd+d)=pre2*array_out(a,bc+c,b,bd+d)&
                                                 &+pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,c,bb+b,bd+d)=pre2*array_out(a,c,bb+b,bd+d)&
                                                 &+pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,bc+c,bb+b,d)=pre2*array_out(a,bc+c,bb+b,d)&
                                                 &+pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do b=db2+1,db
                do c=dc2+1,dc
                  do a=0,bcntr
                    array_out(ba+a,c,b,bd+d)=pre2*array_out(ba+a,c,b,bd+d)&
                                                 &+pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do b=db2+1,db
                do c=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,bc+c,b,d)=pre2*array_out(ba+a,bc+c,b,d)&
                                                 &+pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,c,bb+b,d)=pre2*array_out(ba+a,c,bb+b,d)&
                                                 &+pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do b=db2+1,db
              do c=dc2+1,dc
                do a=0,bcntr
                  array_out(ba+a,c,b,d)=pre2*array_out(ba+a,c,b,d)&
                                               &+pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(a,c,bb+b,d)=pre2*array_out(a,c,bb+b,d)&
                                               &+pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do b=db2+1,db
              do c=0,bcntr
                do a=da2+1,da
                  array_out(a,bc+c,b,d)=pre2*array_out(a,bc+c,b,d)&
                                               &+pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do d=0,bcntr
            do b=db2+1,db
              do c=dc2+1,dc
                do a=da2+1,da
                  array_out(a,c,b,bd+d)=pre2*array_out(a,c,b,bd+d)&
                                               &+pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do b=db2+1,db
          do d=dd2+1,dd
            do c=dc2+1,dc
              do a=da2+1,da
                array_out(a,c,b,d)=pre2*array_out(a,c,b,d)&
                                             &+pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    else
      call lsquit("MANUAL REORDERING:Case not found",-1)
    endif

    !print *,"check 1324"
    !if(pre2==0.0E0_realk)then
    !  do c=1,dc
    !    do a=1,da
    !      do d=1,dd
    !        do b=1,db
    !          if(array_out(a,c,b,d)/=pre1*array_in(a,b,c,d))then
    !            print *,"1324 reordering not correct",a,b,c,d,da,db,dc,dd
    !            stop 0
    !          endif
    !        enddo
    !      enddo
    !    enddo
    !  enddo
    !endif
  end subroutine manual_1324_reordering

  !\> \brief reorder a 4 diensional array  to get the indices
  !   in the order 1 3 4 2 , this is a quite expensive reordering
  !   and thus requires additional attnetion
  !\> \autor Patrick Ettenhuber
  !\> \date Movember 2012
  subroutine manual_1342_reordering(bs,dims,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the original array
    integer, intent(in) :: dims(4)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1,pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3),dims(4))
    !> reordered array
    real(realk),intent(inout) :: array_out(dims(1),dims(3),dims(4),dims(2))
    integer :: a,b,c,d,ba,bb,bc,bd,da,db,dc,dd,da2,db2,dc2,dd2,bcntr
    logical :: moda,modb,modc,modd
    
    !routine go get [a b c d] --> [d b c a]
    da=dims(1)
    db=dims(2)
    dc=dims(3)
    dd=dims(4)

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs
    dd2=(dd/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)
    modd=(mod(dd,bs)>0)

    bcntr=bs-1
    if(pre2==0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do bc=1,dc2,bs
              do ba=1,da2,bs
     
                do b=0,bcntr
                  do d=0,bcntr
                    do c=0,bcntr
                      do a=0,bcntr
                        array_out(ba+a,bc+c,bd+d,bb+b)=array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do bc=1,dc2,bs
     
                do b=0,bcntr
                  do d=0,bcntr
                    do c=0,bcntr
                      do a=da2+1,da
                        array_out(a,bc+c,bd+d,bb+b)=array_in(a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do b=db2+1,db
                do d=0,bcntr
                  do c=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bc+c,bd+d,b)=array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,c,bd+d,bb+b)=array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bc+c,d,bb+b)=array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do b=db2+1,db
              do d=0,bcntr
                do c=0,bcntr
                  do a=da2+1,da
                    array_out(a,bc+c,bd+d,b)=array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,c,bd+d,bb+b)=array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,bc+c,d,bb+b)=array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do b=db2+1,db
              do d=0,bcntr
                do c=dc2+1,dc
                  do a=0,bcntr
                    array_out(ba+a,c,bd+d,b)=array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do b=db2+1,db
              do d=dd2+1,dd
                do c=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,bc+c,d,b)=array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,c,d,bb+b)=array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do b=db2+1,db
            do d=dd2+1,dd
              do c=dc2+1,dc
                do a=0,bcntr
                  array_out(ba+a,c,d,b)=array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(a,c,d,bb+b)=array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do b=db2+1,db
            do d=dd2+1,dd
              do c=0,bcntr
                do a=da2+1,da
                  array_out(a,bc+c,d,b)=array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do b=db2+1,db
            do d=0,bcntr
              do c=dc2+1,dc
                do a=da2+1,da
                  array_out(a,c,bd+d,b)=array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do b=db2+1,db
          do d=dd2+1,dd
            do c=dc2+1,dc
              do a=da2+1,da
                array_out(a,c,d,b)=array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2==0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do bc=1,dc2,bs
              do ba=1,da2,bs
     
                do b=0,bcntr
                  do d=0,bcntr
                    do c=0,bcntr
                      do a=0,bcntr
                        array_out(ba+a,bc+c,bd+d,bb+b)=pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do bc=1,dc2,bs
     
                do b=0,bcntr
                  do d=0,bcntr
                    do c=0,bcntr
                      do a=da2+1,da
                        array_out(a,bc+c,bd+d,bb+b)=pre1*array_in(a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do b=db2+1,db
                do d=0,bcntr
                  do c=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bc+c,bd+d,b)=pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,c,bd+d,bb+b)=pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bc+c,d,bb+b)=pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do b=db2+1,db
              do d=0,bcntr
                do c=0,bcntr
                  do a=da2+1,da
                    array_out(a,bc+c,bd+d,b)=pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,c,bd+d,bb+b)=pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,bc+c,d,bb+b)=pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do b=db2+1,db
              do d=0,bcntr
                do c=dc2+1,dc
                  do a=0,bcntr
                    array_out(ba+a,c,bd+d,b)=pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do b=db2+1,db
              do d=dd2+1,dd
                do c=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,bc+c,d,b)=pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,c,d,bb+b)=pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do b=db2+1,db
            do d=dd2+1,dd
              do c=dc2+1,dc
                do a=0,bcntr
                  array_out(ba+a,c,d,b)=pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(a,c,d,bb+b)=pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do b=db2+1,db
            do d=dd2+1,dd
              do c=0,bcntr
                do a=da2+1,da
                  array_out(a,bc+c,d,b)=pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do b=db2+1,db
            do d=0,bcntr
              do c=dc2+1,dc
                do a=da2+1,da
                  array_out(a,c,bd+d,b)=pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do b=db2+1,db
          do d=dd2+1,dd
            do c=dc2+1,dc
              do a=da2+1,da
                array_out(a,c,d,b)=pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2/=0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do bc=1,dc2,bs
              do ba=1,da2,bs
     
                do b=0,bcntr
                  do d=0,bcntr
                    do c=0,bcntr
                      do a=0,bcntr
                        array_out(ba+a,bc+c,bd+d,bb+b)=pre2*array_out(ba+a,bc+c,bd+d,bb+b)&
                                                     &+array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do bc=1,dc2,bs
     
                do b=0,bcntr
                  do d=0,bcntr
                    do c=0,bcntr
                      do a=da2+1,da
                        array_out(a,bc+c,bd+d,bb+b)=pre2*array_out(a,bc+c,bd+d,bb+b)&
                                                  &+array_in(a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do b=db2+1,db
                do d=0,bcntr
                  do c=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bc+c,bd+d,b)=pre2*array_out(ba+a,bc+c,bd+d,b)&
                                                   &+array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,c,bd+d,bb+b)=pre2*array_out(ba+a,c,bd+d,bb+b)&
                                                   &+array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bc+c,d,bb+b)=pre2*array_out(ba+a,bc+c,d,bb+b)&
                                                   &+array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do b=db2+1,db
              do d=0,bcntr
                do c=0,bcntr
                  do a=da2+1,da
                    array_out(a,bc+c,bd+d,b)=pre2*array_out(a,bc+c,bd+d,b)&
                                                 &+array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,c,bd+d,bb+b)=pre2*array_out(a,c,bd+d,bb+b)&
                                                 &+array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,bc+c,d,bb+b)=pre2*array_out(a,bc+c,d,bb+b)&
                                                 &+array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do b=db2+1,db
              do d=0,bcntr
                do c=dc2+1,dc
                  do a=0,bcntr
                    array_out(ba+a,c,bd+d,b)=pre2*array_out(ba+a,c,bd+d,b)&
                                                 &+array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do b=db2+1,db
              do d=dd2+1,dd
                do c=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,bc+c,d,b)=pre2*array_out(ba+a,bc+c,d,b)&
                                                 &+array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,c,d,bb+b)=pre2*array_out(ba+a,c,d,bb+b)&
                                                 &+array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do b=db2+1,db
            do d=dd2+1,dd
              do c=dc2+1,dc
                do a=0,bcntr
                  array_out(ba+a,c,d,b)=pre2*array_out(ba+a,c,d,b)&
                                               &+array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(a,c,d,bb+b)=pre2*array_out(a,c,d,bb+b)&
                                               &+array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do b=db2+1,db
            do d=dd2+1,dd
              do c=0,bcntr
                do a=da2+1,da
                  array_out(a,bc+c,d,b)=pre2*array_out(a,bc+c,d,b)&
                                               &+array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do b=db2+1,db
            do d=0,bcntr
              do c=dc2+1,dc
                do a=da2+1,da
                  array_out(a,c,bd+d,b)=pre2*array_out(a,c,bd+d,b)&
                                               &+array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do b=db2+1,db
          do d=dd2+1,dd
            do c=dc2+1,dc
              do a=da2+1,da
                array_out(a,c,d,b)=pre2*array_out(a,c,d,b)&
                                             &+array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2/=0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do bc=1,dc2,bs
              do ba=1,da2,bs
     
                do b=0,bcntr
                  do d=0,bcntr
                    do c=0,bcntr
                      do a=0,bcntr
                        array_out(ba+a,bc+c,bd+d,bb+b)=pre2*array_out(ba+a,bc+c,bd+d,bb+b)&
                                                     &+pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do bc=1,dc2,bs
     
                do b=0,bcntr
                  do d=0,bcntr
                    do c=0,bcntr
                      do a=da2+1,da
                        array_out(a,bc+c,bd+d,bb+b)=pre2*array_out(a,bc+c,bd+d,bb+b)&
                                                  &+pre1*array_in(a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do b=db2+1,db
                do d=0,bcntr
                  do c=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bc+c,bd+d,b)=pre2*array_out(ba+a,bc+c,bd+d,b)&
                                                   &+pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,c,bd+d,bb+b)=pre2*array_out(ba+a,c,bd+d,bb+b)&
                                                   &+pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(ba+a,bc+c,d,bb+b)=pre2*array_out(ba+a,bc+c,d,bb+b)&
                                                   &+pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do b=db2+1,db
              do d=0,bcntr
                do c=0,bcntr
                  do a=da2+1,da
                    array_out(a,bc+c,bd+d,b)=pre2*array_out(a,bc+c,bd+d,b)&
                                                 &+pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,c,bd+d,bb+b)=pre2*array_out(a,c,bd+d,bb+b)&
                                                 &+pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(a,bc+c,d,bb+b)=pre2*array_out(a,bc+c,d,bb+b)&
                                                 &+pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do b=db2+1,db
              do d=0,bcntr
                do c=dc2+1,dc
                  do a=0,bcntr
                    array_out(ba+a,c,bd+d,b)=pre2*array_out(ba+a,c,bd+d,b)&
                                                 &+pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do b=db2+1,db
              do d=dd2+1,dd
                do c=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,bc+c,d,b)=pre2*array_out(ba+a,bc+c,d,b)&
                                                 &+pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,c,d,bb+b)=pre2*array_out(ba+a,c,d,bb+b)&
                                                 &+pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do b=db2+1,db
            do d=dd2+1,dd
              do c=dc2+1,dc
                do a=0,bcntr
                  array_out(ba+a,c,d,b)=pre2*array_out(ba+a,c,d,b)&
                                               &+pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(a,c,d,bb+b)=pre2*array_out(a,c,d,bb+b)&
                                               &+pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do b=db2+1,db
            do d=dd2+1,dd
              do c=0,bcntr
                do a=da2+1,da
                  array_out(a,bc+c,d,b)=pre2*array_out(a,bc+c,d,b)&
                                               &+pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do b=db2+1,db
            do d=0,bcntr
              do c=dc2+1,dc
                do a=da2+1,da
                  array_out(a,c,bd+d,b)=pre2*array_out(a,c,bd+d,b)&
                                               &+pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do b=db2+1,db
          do d=dd2+1,dd
            do c=dc2+1,dc
              do a=da2+1,da
                array_out(a,c,d,b)=pre2*array_out(a,c,d,b)&
                                             &+pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    else
      call lsquit("MANUAL REORDERING:Case not found",-1)
    endif

    !print *,"check 1342"
    !if(pre2==0.0E0_realk)then
    !  do c=1,dc
    !    do a=1,da
    !      do d=1,dd
    !        do b=1,db
    !          if(array_out(a,c,d,b)/=pre1*array_in(a,b,c,d))then
    !            print *,"1342 reordering not correct",a,b,c,d,da,db,dc,dd
    !            stop 0
    !          endif
    !        enddo
    !      enddo
    !    enddo
    !  enddo
    !endif
  end subroutine manual_1342_reordering
  !\> \brief reorder a 4 diensional array  to get the indices
  !   in the order 2 4 1 3 , this is a quite expensive reordering
  !   and thus requires additional attnetion
  !\> \autor Patrick Ettenhuber
  !\> \date Movember 2012
  subroutine manual_3214_reordering(bs,dims,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the original array
    integer, intent(in) :: dims(4)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1,pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3),dims(4))
    !> reordered array
    real(realk),intent(inout) :: array_out(dims(3),dims(2),dims(1),dims(4))
    integer :: a,b,c,d,ba,bb,bc,bd,da,db,dc,dd,da2,db2,dc2,dd2,bcntr
    logical :: moda,modb,modc,modd
    
    !routine go get [a b c d] --> [d b c a]
    da=dims(1)
    db=dims(2)
    dc=dims(3)
    dd=dims(4)

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs
    dd2=(dd/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)
    modd=(mod(dd,bs)>0)

    bcntr=bs-1
    if(pre2==0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bb=1,db2,bs
              do bc=1,dc2,bs
     
                do d=0,bcntr
                  do a=0,bcntr
                    do b=0,bcntr
                      do c=0,bcntr
                        array_out(bc+c,bb+b,ba+a,bd+d)=array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do bc=1,dc2,bs
     
              do d=0,bcntr
                do a=da2+1,da
                  do b=0,bcntr
                    do c=0,bcntr
                      array_out(bc+c,bb+b,a,bd+d)=array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bc=1,dc2,bs
     
              do d=0,bcntr
                do a=0,bcntr
                  do b=db2+1,db
                    do c=0,bcntr
                      array_out(bc+c,b,ba+a,bd+d)=array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bb=1,db2,bs
     
                do d=0,bcntr
                  do c=dc2+1,dc
                    do b=0,bcntr
                      do a=0,bcntr
                        array_out(c,bb+b,ba+a,bd+d)=array_in(ba+a,bb+b,c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
            do bc=1,dc2,bs
     
              do d=dd2+1,dd
                do a=0,bcntr
                  do b=0,bcntr
                    do c=0,bcntr
                      array_out(bc+c,bb+b,ba+a,d)=array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do d=0,bcntr
              do a=da2+1,da
                do b=db2+1,db
                  do c=0,bcntr
                    array_out(bc+c,b,a,bd+d)=array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do d=0,bcntr
              do a=da2+1,da
                do b=0,bcntr
                  do c=dc2+1,dc
                    array_out(c,bb+b,a,bd+d)=array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
     
            do d=dd2+1,dd
              do a=da2+1,da
                do b=0,bcntr
                  do c=0,bcntr
                    array_out(bc+c,bb+b,a,d)=array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(c,b,ba+a,bd+d)=array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
     
            do d=dd2+1,dd
              do a=0,bcntr
                do b=db2+1,db
                  do c=0,bcntr
                    array_out(bc+c,b,ba+a,d)=array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(c,bb+b,ba+a,d)=array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(c,b,ba+a,d)=array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do a=da2+1,da
              do b=0,bcntr
                do c=dc2+1,dc
                  array_out(c,bb+b,a,d)=array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     

        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
       do bc=1,dc2,bs
     
         do d=dd2+1,dd
           do a=da2+1,da
             do b=db2+1,db
               do c=0,bcntr
                 array_out(bc+c,b,a,d)=array_in(a,b,bc+c,d)
               enddo
             enddo
           enddo
         enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do d=0,bcntr
            do a=da2+1,da
              do b=db2+1,db
                do c=dc2+1,dc
                  array_out(c,b,a,bd+d)=array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do d=dd2+1,dd
          do a=da2+1,da
            do b=db2+1,db
              do c=dc2+1,dc
                array_out(c,b,a,d)=array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2==0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bb=1,db2,bs
              do bc=1,dc2,bs
     
                do d=0,bcntr
                  do a=0,bcntr
                    do b=0,bcntr
                      do c=0,bcntr
                        array_out(bc+c,bb+b,ba+a,bd+d)=pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do bc=1,dc2,bs
     
              do d=0,bcntr
                do a=da2+1,da
                  do b=0,bcntr
                    do c=0,bcntr
                      array_out(bc+c,bb+b,a,bd+d)=pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bc=1,dc2,bs
     
              do d=0,bcntr
                do a=0,bcntr
                  do b=db2+1,db
                    do c=0,bcntr
                      array_out(bc+c,b,ba+a,bd+d)=pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bb=1,db2,bs
     
                do d=0,bcntr
                  do c=dc2+1,dc
                    do b=0,bcntr
                      do a=0,bcntr
                        array_out(c,bb+b,ba+a,bd+d)=pre1*array_in(ba+a,bb+b,c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
            do bc=1,dc2,bs
     
              do d=dd2+1,dd
                do a=0,bcntr
                  do b=0,bcntr
                    do c=0,bcntr
                      array_out(bc+c,bb+b,ba+a,d)=pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do d=0,bcntr
              do a=da2+1,da
                do b=db2+1,db
                  do c=0,bcntr
                    array_out(bc+c,b,a,bd+d)=pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do d=0,bcntr
              do a=da2+1,da
                do b=0,bcntr
                  do c=dc2+1,dc
                    array_out(c,bb+b,a,bd+d)=pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
     
            do d=dd2+1,dd
              do a=da2+1,da
                do b=0,bcntr
                  do c=0,bcntr
                    array_out(bc+c,bb+b,a,d)=pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(c,b,ba+a,bd+d)=pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
     
            do d=dd2+1,dd
              do a=0,bcntr
                do b=db2+1,db
                  do c=0,bcntr
                    array_out(bc+c,b,ba+a,d)=pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(c,bb+b,ba+a,d)=pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(c,b,ba+a,d)=pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do a=da2+1,da
              do b=0,bcntr
                do c=dc2+1,dc
                  array_out(c,bb+b,a,d)=pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     

        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
       do bc=1,dc2,bs
     
         do d=dd2+1,dd
           do a=da2+1,da
             do b=db2+1,db
               do c=0,bcntr
                 array_out(bc+c,b,a,d)=pre1*array_in(a,b,bc+c,d)
               enddo
             enddo
           enddo
         enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do d=0,bcntr
            do a=da2+1,da
              do b=db2+1,db
                do c=dc2+1,dc
                  array_out(c,b,a,bd+d)=pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do d=dd2+1,dd
          do a=da2+1,da
            do b=db2+1,db
              do c=dc2+1,dc
                array_out(c,b,a,d)=pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2/=0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bb=1,db2,bs
              do bc=1,dc2,bs
     
                do d=0,bcntr
                  do a=0,bcntr
                    do b=0,bcntr
                      do c=0,bcntr
                        array_out(bc+c,bb+b,ba+a,bd+d)=pre2*array_out(bc+c,bb+b,ba+a,bd+d)&
                                                     &+array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do bc=1,dc2,bs
     
              do d=0,bcntr
                do a=da2+1,da
                  do b=0,bcntr
                    do c=0,bcntr
                      array_out(bc+c,bb+b,a,bd+d)=pre2*array_out(bc+c,bb+b,a,bd+d)&
                                                   &+array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bc=1,dc2,bs
     
              do d=0,bcntr
                do a=0,bcntr
                  do b=db2+1,db
                    do c=0,bcntr
                      array_out(bc+c,b,ba+a,bd+d)=pre2*array_out(bc+c,b,ba+a,bd+d)&
                                                   &+array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bb=1,db2,bs
     
                do d=0,bcntr
                  do c=dc2+1,dc
                    do b=0,bcntr
                      do a=0,bcntr
                        array_out(c,bb+b,ba+a,bd+d)=pre2*array_out(c,bb+b,ba+a,bd+d)&
                                                     &+array_in(ba+a,bb+b,c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
            do bc=1,dc2,bs
     
              do d=dd2+1,dd
                do a=0,bcntr
                  do b=0,bcntr
                    do c=0,bcntr
                      array_out(bc+c,bb+b,ba+a,d)=pre2*array_out(bc+c,bb+b,ba+a,d)&
                                                   &+array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do d=0,bcntr
              do a=da2+1,da
                do b=db2+1,db
                  do c=0,bcntr
                    array_out(bc+c,b,a,bd+d)=pre2*array_out(bc+c,b,a,bd+d)&
                                                 &+array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do d=0,bcntr
              do a=da2+1,da
                do b=0,bcntr
                  do c=dc2+1,dc
                    array_out(c,bb+b,a,bd+d)=pre2*array_out(c,bb+b,a,bd+d)&
                                                 &+array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
     
            do d=dd2+1,dd
              do a=da2+1,da
                do b=0,bcntr
                  do c=0,bcntr
                    array_out(bc+c,bb+b,a,d)=pre2*array_out(bc+c,bb+b,a,d)&
                                                 &+array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(c,b,ba+a,bd+d)=pre2*array_out(c,b,ba+a,bd+d)&
                                                 &+array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
     
            do d=dd2+1,dd
              do a=0,bcntr
                do b=db2+1,db
                  do c=0,bcntr
                    array_out(bc+c,b,ba+a,d)=pre2*array_out(bc+c,b,ba+a,d)&
                                                 &+array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(c,bb+b,ba+a,d)=pre2*array_out(c,bb+b,ba+a,d)&
                                                 &+array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(c,b,ba+a,d)=pre2*array_out(c,b,ba+a,d)&
                                               &+array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do a=da2+1,da
              do b=0,bcntr
                do c=dc2+1,dc
                  array_out(c,bb+b,a,d)=pre2*array_out(c,bb+b,a,d)&
                                               &+array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     

        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
       do bc=1,dc2,bs
     
         do d=dd2+1,dd
           do a=da2+1,da
             do b=db2+1,db
               do c=0,bcntr
                 array_out(bc+c,b,a,d)=pre2*array_out(bc+c,b,a,d)&
                                              &+array_in(a,b,bc+c,d)
               enddo
             enddo
           enddo
         enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do d=0,bcntr
            do a=da2+1,da
              do b=db2+1,db
                do c=dc2+1,dc
                  array_out(c,b,a,bd+d)=pre2*array_out(c,b,a,bd+d)&
                                               &+array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do d=dd2+1,dd
          do a=da2+1,da
            do b=db2+1,db
              do c=dc2+1,dc
                array_out(c,b,a,d)=pre2*array_out(c,b,a,d)&
                                             &+array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
     
      endif
    elseif(pre2/=0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bb=1,db2,bs
              do bc=1,dc2,bs
     
                do d=0,bcntr
                  do a=0,bcntr
                    do b=0,bcntr
                      do c=0,bcntr
                        array_out(bc+c,bb+b,ba+a,bd+d)=pre2*array_out(bc+c,bb+b,ba+a,bd+d)&
                                                     &+pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do bc=1,dc2,bs
     
              do d=0,bcntr
                do a=da2+1,da
                  do b=0,bcntr
                    do c=0,bcntr
                      array_out(bc+c,bb+b,a,bd+d)=pre2*array_out(bc+c,bb+b,a,bd+d)&
                                                   &+pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bc=1,dc2,bs
     
              do d=0,bcntr
                do a=0,bcntr
                  do b=db2+1,db
                    do c=0,bcntr
                      array_out(bc+c,b,ba+a,bd+d)=pre2*array_out(bc+c,b,ba+a,bd+d)&
                                                   &+pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bb=1,db2,bs
     
                do d=0,bcntr
                  do c=dc2+1,dc
                    do b=0,bcntr
                      do a=0,bcntr
                        array_out(c,bb+b,ba+a,bd+d)=pre2*array_out(c,bb+b,ba+a,bd+d)&
                                                     &+pre1*array_in(ba+a,bb+b,c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
            do bc=1,dc2,bs
     
              do d=dd2+1,dd
                do a=0,bcntr
                  do b=0,bcntr
                    do c=0,bcntr
                      array_out(bc+c,bb+b,ba+a,d)=pre2*array_out(bc+c,bb+b,ba+a,d)&
                                                   &+pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do d=0,bcntr
              do a=da2+1,da
                do b=db2+1,db
                  do c=0,bcntr
                    array_out(bc+c,b,a,bd+d)=pre2*array_out(bc+c,b,a,bd+d)&
                                                 &+pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do d=0,bcntr
              do a=da2+1,da
                do b=0,bcntr
                  do c=dc2+1,dc
                    array_out(c,bb+b,a,bd+d)=pre2*array_out(c,bb+b,a,bd+d)&
                                                 &+pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0)then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
     
            do d=dd2+1,dd
              do a=da2+1,da
                do b=0,bcntr
                  do c=0,bcntr
                    array_out(bc+c,bb+b,a,d)=pre2*array_out(bc+c,bb+b,a,d)&
                                                 &+pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(c,b,ba+a,bd+d)=pre2*array_out(c,b,ba+a,bd+d)&
                                                 &+pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
     
            do d=dd2+1,dd
              do a=0,bcntr
                do b=db2+1,db
                  do c=0,bcntr
                    array_out(bc+c,b,ba+a,d)=pre2*array_out(bc+c,b,ba+a,d)&
                                                 &+pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(c,bb+b,ba+a,d)=pre2*array_out(c,bb+b,ba+a,d)&
                                                 &+pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(c,b,ba+a,d)=pre2*array_out(c,b,ba+a,d)&
                                               &+pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do a=da2+1,da
              do b=0,bcntr
                do c=dc2+1,dc
                  array_out(c,bb+b,a,d)=pre2*array_out(c,bb+b,a,d)&
                                               &+pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     

        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
       do bc=1,dc2,bs
     
         do d=dd2+1,dd
           do a=da2+1,da
             do b=db2+1,db
               do c=0,bcntr
                 array_out(bc+c,b,a,d)=pre2*array_out(bc+c,b,a,d)&
                                              &+pre1*array_in(a,b,bc+c,d)
               enddo
             enddo
           enddo
         enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do d=0,bcntr
            do a=da2+1,da
              do b=db2+1,db
                do c=dc2+1,dc
                  array_out(c,b,a,bd+d)=pre2*array_out(c,b,a,bd+d)&
                                               &+pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do d=dd2+1,dd
          do a=da2+1,da
            do b=db2+1,db
              do c=dc2+1,dc
                array_out(c,b,a,d)=pre2*array_out(c,b,a,d)&
                                             &+pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    else
      call lsquit("MANUAL REORDERING:Case not found",-1)
    endif


    !print *,"checking 3214"
    !if(pre2==0.0E0_realk)then
    !  do c=1,dc
    !    do a=1,da
    !      do d=1,dd
    !        do b=1,db
    !          if(array_out(c,b,a,d)/=pre1*array_in(a,b,c,d))then
    !            print *,"3214 reordering not correct",a,b,c,d,da,db,dc,dd
    !            stop 0
    !          endif
    !        enddo
    !      enddo
    !    enddo
    !  enddo
    !endif
  end subroutine manual_3214_reordering

  !\> \brief reorder a 4 diensional array  to get the indices
  !   in the order 2 4 1 3 , this is a quite expensive reordering
  !   and thus requires additional attnetion
  !\> \autor Patrick Ettenhuber
  !\> \date Movember 2012
  subroutine manual_2413_reordering(bs,dims,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the original array
    integer, intent(in) :: dims(4)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1,pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3),dims(4))
    !> reordered array
    real(realk),intent(inout) :: array_out(dims(2),dims(4),dims(1),dims(3))
    integer :: a,b,c,d,ba,bb,bc,bd,da,db,dc,dd,da2,db2,dc2,dd2,bcntr
    logical :: moda,modb,modc,modd
    
    !routine go get [a b c d] --> [d b c a]
    da=dims(1)
    db=dims(2)
    dc=dims(3)
    dd=dims(4)

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs
    dd2=(dd/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)
    modd=(mod(dd,bs)>0)

    bcntr=bs-1

    if(pre2==0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
            do bd=1,dd2,bs
              do bb=1,db2,bs
     
                do c=0,bcntr
                  do a=0,bcntr
                    do d=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,bd+d,ba+a,bc+c)=array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
     
              do c=0,bcntr
                do a=da2+1,da
                  do d=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bd+d,a,bc+c)=array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do a=0,bcntr
                        array_out(b,bd+d,ba+a,bc+c)=array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
     
              do c=dc2+1,dc
                do a=0,bcntr
                  do d=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bd+d,ba+a,c)=array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bb+b,d,ba+a,bc+c)=array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do c=0,bcntr
              do a=da2+1,da
                do d=0,bcntr
                  do b=db2+1,db
                    array_out(b,bd+d,a,bc+c)=array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do c=dc2+1,dc
              do a=da2+1,da
                do d=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,bd+d,a,c)=array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do c=0,bcntr
              do a=da2+1,da
                do d=dd2+1,dd
                  do b=0,bcntr
                    array_out(bb+b,d,a,bc+c)=array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,bd+d,ba+a,c)=array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,d,ba+a,bc+c)=array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(bb+b,d,ba+a,c)=array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,d,ba+a,c)=array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do c=dc2+1,dc
            do a=da2+1,da
              do d=dd2+1,dd
                do b=0,bcntr
                  array_out(bb+b,d,a,c)=array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do c=0,bcntr
            do a=da2+1,da
              do d=dd2+1,dd
                do b=db2+1,db
                  array_out(b,d,a,bc+c)=array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do c=dc2+1,dc
            do a=da2+1,da
              do d=0,bcntr
                do b=db2+1,db
                  array_out(b,bd+d,a,c)=array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do c=dc2+1,dc
          do a=da2+1,da
            do d=dd2+1,dd
              do b=db2+1,db
                array_out(b,d,a,c)=array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2==0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
            do bd=1,dd2,bs
              do bb=1,db2,bs
     
                do c=0,bcntr
                  do a=0,bcntr
                    do d=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,bd+d,ba+a,bc+c)=pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
     
              do c=0,bcntr
                do a=da2+1,da
                  do d=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bd+d,a,bc+c)=pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do a=0,bcntr
                        array_out(b,bd+d,ba+a,bc+c)=pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
     
              do c=dc2+1,dc
                do a=0,bcntr
                  do d=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bd+d,ba+a,c)=pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bb+b,d,ba+a,bc+c)=pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do c=0,bcntr
              do a=da2+1,da
                do d=0,bcntr
                  do b=db2+1,db
                    array_out(b,bd+d,a,bc+c)=pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do c=dc2+1,dc
              do a=da2+1,da
                do d=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,bd+d,a,c)=pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do c=0,bcntr
              do a=da2+1,da
                do d=dd2+1,dd
                  do b=0,bcntr
                    array_out(bb+b,d,a,bc+c)=pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,bd+d,ba+a,c)=pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,d,ba+a,bc+c)=pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(bb+b,d,ba+a,c)=pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,d,ba+a,c)=pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do c=dc2+1,dc
            do a=da2+1,da
              do d=dd2+1,dd
                do b=0,bcntr
                  array_out(bb+b,d,a,c)=pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do c=0,bcntr
            do a=da2+1,da
              do d=dd2+1,dd
                do b=db2+1,db
                  array_out(b,d,a,bc+c)=pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do c=dc2+1,dc
            do a=da2+1,da
              do d=0,bcntr
                do b=db2+1,db
                  array_out(b,bd+d,a,c)=pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do c=dc2+1,dc
          do a=da2+1,da
            do d=dd2+1,dd
              do b=db2+1,db
                array_out(b,d,a,c)=pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2/=0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
            do bd=1,dd2,bs
              do bb=1,db2,bs
     
                do c=0,bcntr
                  do a=0,bcntr
                    do d=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,bd+d,ba+a,bc+c)=pre2*array_out(bb+b,bd+d,ba+a,bc+c)&
                                                     &+array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
     
              do c=0,bcntr
                do a=da2+1,da
                  do d=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bd+d,a,bc+c)=pre2*array_out(bb+b,bd+d,a,bc+c)&
                                                   &+array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do a=0,bcntr
                        array_out(b,bd+d,ba+a,bc+c)=pre2*array_out(b,bd+d,ba+a,bc+c)&
                                                     &+array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
     
              do c=dc2+1,dc
                do a=0,bcntr
                  do d=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bd+d,ba+a,c)=pre2*array_out(bb+b,bd+d,ba+a,c)&
                                                &+array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bb+b,d,ba+a,bc+c)=pre2*array_out(bb+b,d,ba+a,bc+c)&
                                                   &+array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do c=0,bcntr
              do a=da2+1,da
                do d=0,bcntr
                  do b=db2+1,db
                    array_out(b,bd+d,a,bc+c)=pre2*array_out(b,bd+d,a,bc+c)&
                                                 &+array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do c=dc2+1,dc
              do a=da2+1,da
                do d=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,bd+d,a,c)=pre2*array_out(bb+b,bd+d,a,c)&
                                                 &+array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do c=0,bcntr
              do a=da2+1,da
                do d=dd2+1,dd
                  do b=0,bcntr
                    array_out(bb+b,d,a,bc+c)=pre2*array_out(bb+b,d,a,bc+c)&
                                                 &+array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,bd+d,ba+a,c)=pre2*array_out(b,bd+d,ba+a,c)&
                                                 &+array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,d,ba+a,bc+c)=pre2*array_out(b,d,ba+a,bc+c)&
                                                 &+array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(bb+b,d,ba+a,c)=pre2*array_out(bb+b,d,ba+a,c)&
                                                 &+array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,d,ba+a,c)=pre2*array_out(b,d,ba+a,c)&
                                               &+array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do c=dc2+1,dc
            do a=da2+1,da
              do d=dd2+1,dd
                do b=0,bcntr
                  array_out(bb+b,d,a,c)=pre2*array_out(bb+b,d,a,c)&
                                               &+array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do c=0,bcntr
            do a=da2+1,da
              do d=dd2+1,dd
                do b=db2+1,db
                  array_out(b,d,a,bc+c)=pre2*array_out(b,d,a,bc+c)&
                                               &+array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do c=dc2+1,dc
            do a=da2+1,da
              do d=0,bcntr
                do b=db2+1,db
                  array_out(b,bd+d,a,c)=pre2*array_out(b,bd+d,a,c)&
                                               &+array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do c=dc2+1,dc
          do a=da2+1,da
            do d=dd2+1,dd
              do b=db2+1,db
                array_out(b,d,a,c)=pre2*array_out(b,d,a,c)&
                                             &+array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2/=0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
            do bd=1,dd2,bs
              do bb=1,db2,bs
     
                do c=0,bcntr
                  do a=0,bcntr
                    do d=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,bd+d,ba+a,bc+c)=pre2*array_out(bb+b,bd+d,ba+a,bc+c)&
                                                     &+pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
     
              do c=0,bcntr
                do a=da2+1,da
                  do d=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bd+d,a,bc+c)=pre2*array_out(bb+b,bd+d,a,bc+c)&
                                                   &+pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do a=0,bcntr
                        array_out(b,bd+d,ba+a,bc+c)=pre2*array_out(b,bd+d,ba+a,bc+c)&
                                                     &+pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
     
              do c=dc2+1,dc
                do a=0,bcntr
                  do d=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bd+d,ba+a,c)=pre2*array_out(bb+b,bd+d,ba+a,c)&
                                                &+pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bb+b,d,ba+a,bc+c)=pre2*array_out(bb+b,d,ba+a,bc+c)&
                                                   &+pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do c=0,bcntr
              do a=da2+1,da
                do d=0,bcntr
                  do b=db2+1,db
                    array_out(b,bd+d,a,bc+c)=pre2*array_out(b,bd+d,a,bc+c)&
                                                 &+pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do c=dc2+1,dc
              do a=da2+1,da
                do d=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,bd+d,a,c)=pre2*array_out(bb+b,bd+d,a,c)&
                                                 &+pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do c=0,bcntr
              do a=da2+1,da
                do d=dd2+1,dd
                  do b=0,bcntr
                    array_out(bb+b,d,a,bc+c)=pre2*array_out(bb+b,d,a,bc+c)&
                                                 &+pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,bd+d,ba+a,c)=pre2*array_out(b,bd+d,ba+a,c)&
                                                 &+pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,d,ba+a,bc+c)=pre2*array_out(b,d,ba+a,bc+c)&
                                                 &+pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(bb+b,d,ba+a,c)=pre2*array_out(bb+b,d,ba+a,c)&
                                                 &+pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,d,ba+a,c)=pre2*array_out(b,d,ba+a,c)&
                                               &+pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.db2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
     
          do c=dc2+1,dc
            do a=da2+1,da
              do d=dd2+1,dd
                do b=0,bcntr
                  array_out(bb+b,d,a,c)=pre2*array_out(bb+b,d,a,c)&
                                               &+pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.dc2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
     
          do c=0,bcntr
            do a=da2+1,da
              do d=dd2+1,dd
                do b=db2+1,db
                  array_out(b,d,a,bc+c)=pre2*array_out(b,d,a,bc+c)&
                                               &+pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.modb.and.modc.and.dd2>0)then
        !$OMP DO
        do bd=1,dd2,bs
     
          do c=dc2+1,dc
            do a=da2+1,da
              do d=0,bcntr
                do b=db2+1,db
                  array_out(b,bd+d,a,c)=pre2*array_out(b,bd+d,a,c)&
                                               &+pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
        do c=dc2+1,dc
          do a=da2+1,da
            do d=dd2+1,dd
              do b=db2+1,db
                array_out(b,d,a,c)=pre2*array_out(b,d,a,c)&
                                             &+pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    else
      call lsquit("MANUAL REORDERING:Case not found",-1)
    endif

    !print *,"checking the result 2413"
    !if(pre2==0.0E0_realk)then
    !  do c=1,dc
    !    do a=1,da
    !      do d=1,dd
    !        do b=1,db
    !          if(array_out(b,d,a,c)/=pre1*array_in(a,b,c,d))then
    !            print *,"2413 reordering not correct",a,b,c,d,da,db,dc,dd
    !            stop 0
    !          endif
    !        enddo
    !      enddo
    !    enddo
    !  enddo
    !endif
  end subroutine manual_2413_reordering
  
  !\> \brief reorder a 4 diensional array  to get the indices
  !   in the order 4 2 3 1, this is a quite expensive reordering
  !   and thus requires additional attnetion
  !\> \autor Patrick Ettenhuber
  !\> \date Movember 2012
  subroutine manual_4231_reordering(bs,dims,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the original array
    integer, intent(in) :: dims(4)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1,pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3),dims(4))
    !> reordered array
    real(realk),intent(inout) :: array_out(dims(4),dims(2),dims(3),dims(1))
    integer :: a,b,c,d,ba,bb,bc,bd,da,db,dc,dd,da2,db2,dc2,dd2,bcntr
    logical :: moda,modb,modc,modd
    
    !routine go get [a b c d] --> [d b c a]
    da=dims(1)
    db=dims(2)
    dc=dims(3)
    dd=dims(4)

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs
    dd2=(dd/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)
    modd=(mod(dd,bs)>0)

    bcntr=bs-1


    if(pre2==0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bb=1,db2,bs
              do bd=1,dd2,bs
     
                do a=0,bcntr
                  do c=0,bcntr
                    do b=0,bcntr
                      do d=0,bcntr
                        array_out(bd+d,bb+b,bc+c,ba+a)=array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.dc2>0.and.db2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do bd=1,dd2,bs
      
              do a=da2+1,da
                do c=0,bcntr
                  do b=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,bb+b,bc+c,a)=array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modc.and.db2>0.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
            do bd=1,dd2,bs
      
              do a=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,bb+b,c,ba+a)=array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.dc2>0.and.modb.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
      
              do a=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do d=0,bcntr
                      array_out(bd+d,b,bc+c,ba+a)=array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      !may be better to use the loop over a as the fastest, since more elements than in d
      if(dc2>0.and.db2>0.and.da2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
      
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(d,bb+b,bc+c,ba+a)=array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      !two at the ends
      if(db2>0.and.dd2>0.and.modc.and.moda)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
      
            do a=da2+1,da
              do c=dc2+1,dc
                do b=0,bcntr
                  do d=0,bcntr
                    array_out(bd+d,bb+b,c,a)=array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(dc2>0.and.dd2>0.and.modb.and.moda)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
      
            do a=da2+1,da
              do c=0,bcntr
                do b=db2+1,db
                  do d=0,bcntr
                    array_out(bd+d,b,bc+c,a)=array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(dc2>0.and.db2>0.and.moda.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
      
            do a=da2+1,da
              do c=0,bcntr
                do b=0,bcntr
                  do d=dd2+1,dd
                    array_out(d,bb+b,bc+c,a)=array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.dd2>0.and.modb.and.modc)then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
      
            do a=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do d=0,bcntr
                    array_out(bd+d,b,c,ba+a)=array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(db2>0.and.da2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
      
              do d=dd2+1,dd
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(d,bb+b,c,ba+a)=array_in(ba+a,bb+b,c,d)
                    enddo
                  enddo
                enddo
              enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(dc2>0.and.da2>0.and.modd.and.modb)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
      
             do d=dd2+1,dd
               do c=0,bcntr
                 do b=db2+1,db
                   do a=0,bcntr
                     array_out(d,b,bc+c,ba+a)=array_in(ba+a,b,bc+c,d)
                   enddo
                 enddo
               enddo
             enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(dd2>0.and.moda.and.modb.and.modc)then
        !$OMP DO
        do bd=1,dd2,bs
      
          do a=da2+1,da
            do c=dc2+1,dc
              do b=db2+1,db
                do d=0,bcntr
                  array_out(bd+d,b,c,a)=array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      if(db2>0.and.moda.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
      
          do a=da2+1,da
            do c=dc2+1,dc
              do b=0,bcntr
                do d=dd2+1,dd
                  array_out(d,bb+b,c,a)=array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      if(dc2>0.and.moda.and.modb.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
      
          do a=da2+1,da
            do c=0,bcntr
              do b=db2+1,db
                do d=dd2+1,dd
                  array_out(d,b,bc+c,a)=array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
      
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(d,b,c,ba+a)=array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
      !last block
        do a=da2+1,da
          do c=dc2+1,dc
            do b=db2+1,db
              do d=dd2+1,dd
                array_out(d,b,c,a)=array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif

    elseif(pre2==0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bb=1,db2,bs
              do bd=1,dd2,bs
     
                do a=0,bcntr
                  do c=0,bcntr
                    do b=0,bcntr
                      do d=0,bcntr
                        array_out(bd+d,bb+b,bc+c,ba+a)=pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.dc2>0.and.db2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do bd=1,dd2,bs
      
              do a=da2+1,da
                do c=0,bcntr
                  do b=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,bb+b,bc+c,a)=pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modc.and.db2>0.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
            do bd=1,dd2,bs
      
              do a=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,bb+b,c,ba+a)=pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.dc2>0.and.modb.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
      
              do a=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do d=0,bcntr
                      array_out(bd+d,b,bc+c,ba+a)=pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      !may be better to use the loop over a as the fastest, since more elements than in d
      if(dc2>0.and.db2>0.and.da2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
      
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(d,bb+b,bc+c,ba+a)=pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      !two at the ends
      if(db2>0.and.dd2>0.and.modc.and.moda)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
      
            do a=da2+1,da
              do c=dc2+1,dc
                do b=0,bcntr
                  do d=0,bcntr
                    array_out(bd+d,bb+b,c,a)=pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(dc2>0.and.dd2>0.and.modb.and.moda)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
      
            do a=da2+1,da
              do c=0,bcntr
                do b=db2+1,db
                  do d=0,bcntr
                    array_out(bd+d,b,bc+c,a)=pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(dc2>0.and.db2>0.and.moda.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
      
            do a=da2+1,da
              do c=0,bcntr
                do b=0,bcntr
                  do d=dd2+1,dd
                    array_out(d,bb+b,bc+c,a)=pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.dd2>0.and.modb.and.modc)then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
      
            do a=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do d=0,bcntr
                    array_out(bd+d,b,c,ba+a)=pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(db2>0.and.da2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
      
              do d=dd2+1,dd
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(d,bb+b,c,ba+a)=pre1*array_in(ba+a,bb+b,c,d)
                    enddo
                  enddo
                enddo
              enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(dc2>0.and.da2>0.and.modd.and.modb)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
      
             do d=dd2+1,dd
               do c=0,bcntr
                 do b=db2+1,db
                   do a=0,bcntr
                     array_out(d,b,bc+c,ba+a)=pre1*array_in(ba+a,b,bc+c,d)
                   enddo
                 enddo
               enddo
             enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(dd2>0.and.moda.and.modb.and.modc)then
        !$OMP DO
        do bd=1,dd2,bs
      
          do a=da2+1,da
            do c=dc2+1,dc
              do b=db2+1,db
                do d=0,bcntr
                  array_out(bd+d,b,c,a)=pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      if(db2>0.and.moda.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
      
          do a=da2+1,da
            do c=dc2+1,dc
              do b=0,bcntr
                do d=dd2+1,dd
                  array_out(d,bb+b,c,a)=pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      if(dc2>0.and.moda.and.modb.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
      
          do a=da2+1,da
            do c=0,bcntr
              do b=db2+1,db
                do d=dd2+1,dd
                  array_out(d,b,bc+c,a)=pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
      
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(d,b,c,ba+a)=pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
      !last block
        do a=da2+1,da
          do c=dc2+1,dc
            do b=db2+1,db
              do d=dd2+1,dd
                array_out(d,b,c,a)=pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif

    elseif(pre2/=0.0E0_realk.and.pre1==1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bb=1,db2,bs
              do bd=1,dd2,bs
     
                do a=0,bcntr
                  do c=0,bcntr
                    do b=0,bcntr
                      do d=0,bcntr
                        array_out(bd+d,bb+b,bc+c,ba+a)=pre2*array_out(bd+d,bb+b,bc+c,ba+a)&
                                                     &+array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.dc2>0.and.db2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do bd=1,dd2,bs
      
              do a=da2+1,da
                do c=0,bcntr
                  do b=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,bb+b,bc+c,a)=pre2*array_out(bd+d,bb+b,bc+c,a)&
                                                &+array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modc.and.db2>0.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
            do bd=1,dd2,bs
      
              do a=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,bb+b,c,ba+a)=pre2*array_out(bd+d,bb+b,c,ba+a)&
                                                &+array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.dc2>0.and.modb.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
      
              do a=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do d=0,bcntr
                      array_out(bd+d,b,bc+c,ba+a)=pre2*array_out(bd+d,b,bc+c,ba+a)&
                                                &+array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      !may be better to use the loop over a as the fastest, since more elements than in d
      if(dc2>0.and.db2>0.and.da2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
      
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(d,bb+b,bc+c,ba+a)=pre2*array_out(d,bb+b,bc+c,ba+a)&
                                                &+array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      !two at the ends
      if(db2>0.and.dd2>0.and.modc.and.moda)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
      
            do a=da2+1,da
              do c=dc2+1,dc
                do b=0,bcntr
                  do d=0,bcntr
                    array_out(bd+d,bb+b,c,a)=pre2*array_out(bd+d,bb+b,c,a)&
                                           &+array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(dc2>0.and.dd2>0.and.modb.and.moda)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
      
            do a=da2+1,da
              do c=0,bcntr
                do b=db2+1,db
                  do d=0,bcntr
                    array_out(bd+d,b,bc+c,a)=pre2*array_out(bd+d,b,bc+c,a)&
                                           &+array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(dc2>0.and.db2>0.and.moda.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
      
            do a=da2+1,da
              do c=0,bcntr
                do b=0,bcntr
                  do d=dd2+1,dd
                    array_out(d,bb+b,bc+c,a)=pre2*array_out(d,bb+b,bc+c,a)&
                                           &+array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.dd2>0.and.modb.and.modc)then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
      
            do a=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do d=0,bcntr
                    array_out(bd+d,b,c,ba+a)=pre2*array_out(bd+d,b,c,ba+a)&
                                           &+array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(db2>0.and.da2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
      
              do d=dd2+1,dd
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(d,bb+b,c,ba+a)=pre2*array_out(d,bb+b,c,ba+a)&
                                             &+array_in(ba+a,bb+b,c,d)
                    enddo
                  enddo
                enddo
              enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(dc2>0.and.da2>0.and.modd.and.modb)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
      
             do d=dd2+1,dd
               do c=0,bcntr
                 do b=db2+1,db
                   do a=0,bcntr
                     array_out(d,b,bc+c,ba+a)=pre2*array_out(d,b,bc+c,ba+a)&
                                            &+array_in(ba+a,b,bc+c,d)
                   enddo
                 enddo
               enddo
             enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(dd2>0.and.moda.and.modb.and.modc)then
        !$OMP DO
        do bd=1,dd2,bs
      
          do a=da2+1,da
            do c=dc2+1,dc
              do b=db2+1,db
                do d=0,bcntr
                  array_out(bd+d,b,c,a)=pre2*array_out(bd+d,b,c,a)&
                                      &+array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      if(db2>0.and.moda.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
      
          do a=da2+1,da
            do c=dc2+1,dc
              do b=0,bcntr
                do d=dd2+1,dd
                  array_out(d,bb+b,c,a)=pre2*array_out(d,bb+b,c,a)&
                                      &+array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      if(dc2>0.and.moda.and.modb.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
      
          do a=da2+1,da
            do c=0,bcntr
              do b=db2+1,db
                do d=dd2+1,dd
                  array_out(d,b,bc+c,a)=pre2*array_out(d,b,bc+c,a)&
                                      &+array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
      
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(d,b,c,ba+a)=pre2*array_out(d,b,c,ba+a)&
                                      &+array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
      !last block
        do a=da2+1,da
          do c=dc2+1,dc
            do b=db2+1,db
              do d=dd2+1,dd
                array_out(d,b,c,a)=pre2*array_out(d,b,c,a)&
                                 &+array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif(pre2/=0.0E0_realk.and.pre1/=1.0E0_realk)then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bb=1,db2,bs
              do bd=1,dd2,bs
     
                do a=0,bcntr
                  do c=0,bcntr
                    do b=0,bcntr
                      do d=0,bcntr
                        array_out(bd+d,bb+b,bc+c,ba+a)=pre2*array_out(bd+d,bb+b,bc+c,ba+a)&
                                                     &+pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(moda.and.dc2>0.and.db2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do bd=1,dd2,bs
      
              do a=da2+1,da
                do c=0,bcntr
                  do b=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,bb+b,bc+c,a)=pre2*array_out(bd+d,bb+b,bc+c,a)&
                                                &+pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modc.and.db2>0.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
            do bd=1,dd2,bs
      
              do a=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,bb+b,c,ba+a)=pre2*array_out(bd+d,bb+b,c,ba+a)&
                                                &+pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.dc2>0.and.modb.and.dd2>0)then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
      
              do a=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do d=0,bcntr
                      array_out(bd+d,b,bc+c,ba+a)=pre2*array_out(bd+d,b,bc+c,ba+a)&
                                                &+pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      !may be better to use the loop over a as the fastest, since more elements than in d
      if(dc2>0.and.db2>0.and.da2>0.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
      
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(d,bb+b,bc+c,ba+a)=pre2*array_out(d,bb+b,bc+c,ba+a)&
                                                &+pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      !two at the ends
      if(db2>0.and.dd2>0.and.modc.and.moda)then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
      
            do a=da2+1,da
              do c=dc2+1,dc
                do b=0,bcntr
                  do d=0,bcntr
                    array_out(bd+d,bb+b,c,a)=pre2*array_out(bd+d,bb+b,c,a)&
                                           &+pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(dc2>0.and.dd2>0.and.modb.and.moda)then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
      
            do a=da2+1,da
              do c=0,bcntr
                do b=db2+1,db
                  do d=0,bcntr
                    array_out(bd+d,b,bc+c,a)=pre2*array_out(bd+d,b,bc+c,a)&
                                           &+pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(dc2>0.and.db2>0.and.moda.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
      
            do a=da2+1,da
              do c=0,bcntr
                do b=0,bcntr
                  do d=dd2+1,dd
                    array_out(d,bb+b,bc+c,a)=pre2*array_out(d,bb+b,bc+c,a)&
                                           &+pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.dd2>0.and.modb.and.modc)then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
      
            do a=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do d=0,bcntr
                    array_out(bd+d,b,c,ba+a)=pre2*array_out(bd+d,b,c,ba+a)&
                                           &+pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(db2>0.and.da2>0.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
      
              do d=dd2+1,dd
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(d,bb+b,c,ba+a)=pre2*array_out(d,bb+b,c,ba+a)&
                                             &+pre1*array_in(ba+a,bb+b,c,d)
                    enddo
                  enddo
                enddo
              enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(dc2>0.and.da2>0.and.modd.and.modb)then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
      
             do d=dd2+1,dd
               do c=0,bcntr
                 do b=db2+1,db
                   do a=0,bcntr
                     array_out(d,b,bc+c,ba+a)=pre2*array_out(d,b,bc+c,ba+a)&
                                            &+pre1*array_in(ba+a,b,bc+c,d)
                   enddo
                 enddo
               enddo
             enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(dd2>0.and.moda.and.modb.and.modc)then
        !$OMP DO
        do bd=1,dd2,bs
      
          do a=da2+1,da
            do c=dc2+1,dc
              do b=db2+1,db
                do d=0,bcntr
                  array_out(bd+d,b,c,a)=pre2*array_out(bd+d,b,c,a)&
                                      &+pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      if(db2>0.and.moda.and.modc.and.modd)then
        !$OMP DO
        do bb=1,db2,bs
      
          do a=da2+1,da
            do c=dc2+1,dc
              do b=0,bcntr
                do d=dd2+1,dd
                  array_out(d,bb+b,c,a)=pre2*array_out(d,bb+b,c,a)&
                                      &+pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      if(dc2>0.and.moda.and.modb.and.modd)then
        !$OMP DO
        do bc=1,dc2,bs
      
          do a=da2+1,da
            do c=0,bcntr
              do b=db2+1,db
                do d=dd2+1,dd
                  array_out(d,b,bc+c,a)=pre2*array_out(d,b,bc+c,a)&
                                      &+pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2>0.and.modb.and.modc.and.modd)then
        !$OMP DO
        do ba=1,da2,bs
      
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(d,b,c,ba+a)=pre2*array_out(d,b,c,ba+a)&
                                      &+pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if(moda.and.modb.and.modc.and.modd)then
      !last block
        do a=da2+1,da
          do c=dc2+1,dc
            do b=db2+1,db
              do d=dd2+1,dd
                array_out(d,b,c,a)=pre2*array_out(d,b,c,a)&
                                 &+pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif

    else
      call lsquit("MANUAL REORDERING:Case not found",-1)
    endif
    !if(pre2==0.0E0_realk)then
    !  do a=1,da
    !    do c=1,dc
    !      do b=1,db
    !        do d=1,dd
    !          if(array_out(d,b,c,a)/=pre1*array_in(a,b,c,d))then
    !            print *,"reordering not correct",a,b,c,d,da,db,dc,dd
    !            stop 0
    !          endif
    !        enddo
    !      enddo
    !    enddo
    !  enddo
    !endif
    
  end subroutine manual_4231_reordering

  !\> \brief reorder a 4 diensional array  to get the indices
  !   in the order 4 3 2 1, this is a quite expensive reordering
  !   and thus requires additional attention
  !\> \author Janus Juul Eriksen & Patrick Ettenhuber
  !\> \date November 2012
  subroutine manual_4321_reordering(bs,dims,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the original array
    integer, intent(in) :: dims(4)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1,pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3),dims(4))
    !> reordered array
    real(realk),intent(inout) :: array_out(dims(4),dims(3),dims(2),dims(1))
    integer :: a,b,c,d,ba,bb,bc,bd,da,db,dc,dd,da2,db2,dc2,dd2,bcntr
    logical :: moda,modb,modc,modd
    
    da=dims(1)
    db=dims(2)
    dc=dims(3)
    dd=dims(4)

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs
    dd2=(dd/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)
    modd=(mod(dd,bs)>0)

    bcntr=bs-1

    if (pre2 .eq. 0.0E0_realk .and. pre1 .eq. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
            do bc=1,dc2,bs
              do bd=1,dd2,bs
     
                do a=0,bcntr
                  do b=0,bcntr
                    do c=0,bcntr
                      do d=0,bcntr
                        array_out(bd+d,bc+c,bb+b,ba+a)=array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
      
              do a=da2+1,da
                do b=0,bcntr
                  do c=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,bc+c,bb+b,a)=array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
      
              do a=0,bcntr
                do b=db2+1,db
                  do c=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,bc+c,b,ba+a)=array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
      
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bd+d,c,bb+b,ba+a)=array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (dc2 .gt. 0 .and. db2 .gt. 0 .and. da2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
      
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(d,bc+c,bb+b,ba+a)=array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      !two at the ends
      if (dc2 .gt. 0 .and. dd2 .gt. 0 .and. modb .and. moda) then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
      
            do a=da2+1,da
              do b=db2+1,db
                do c=0,bcntr
                  do d=0,bcntr
                    array_out(bd+d,bc+c,b,a)=array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (db2 .gt. 0 .and. dd2 .gt. 0 .and. modc .and. moda) then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
      
            do a=da2+1,da
              do b=0,bcntr
                do c=dc2+1,dc
                  do d=0,bcntr
                    array_out(bd+d,c,bb+b,a)=array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (dc2 .gt. 0 .and. db2 .gt. 0 .and. moda .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
      
            do a=da2+1,da
              do b=0,bcntr
                do c=0,bcntr
                  do d=dd2+1,dd
                    array_out(d,bc+c,bb+b,a)=array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. dd2 .gt. 0 .and. modb .and. modc) then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
      
            do a=0,bcntr
              do b=db2+1,db
                do c=dc2+1,dc
                  do d=0,bcntr
                    array_out(bd+d,c,b,ba+a)=array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (db2 .gt. 0 .and. da2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
      
              do d=dd2+1,dd
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(d,c,bb+b,ba+a)=array_in(ba+a,bb+b,c,d)
                    enddo
                  enddo
                enddo
              enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (dc2 .gt. 0 .and. da2 .gt. 0 .and. modd .and. modb) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
      
             do d=dd2+1,dd
               do c=0,bcntr
                 do b=db2+1,db
                   do a=0,bcntr
                     array_out(d,bc+c,b,ba+a)=array_in(ba+a,b,bc+c,d)
                   enddo
                 enddo
               enddo
             enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (dd2 .gt. 0 .and. moda .and. modb .and. modc) then
        !$OMP DO
        do bd=1,dd2,bs
      
          do a=da2+1,da
            do b=db2+1,db
              do c=dc2+1,dc
                do d=0,bcntr
                  array_out(bd+d,c,b,a)=array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      if (db2 .gt. 0 .and. moda .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
      
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(d,c,bb+b,a)=array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      if (dc2 .gt. 0 .and. moda .and. modb .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
      
          do a=da2+1,da
            do b=db2+1,db
              do c=0,bcntr
                do d=dd2+1,dd
                  array_out(d,bc+c,b,a)=array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. modc .and. modd) then
        !$OMP DO
        do ba=1,da2,bs
      
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(d,c,b,ba+a)=array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. modc .and. modd) then
      !last block
        do a=da2+1,da
          do b=db2+1,db
            do c=dc2+1,dc
              do d=dd2+1,dd
                array_out(d,c,b,a)=array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif

    elseif (pre2 .eq. 0.0E0_realk .and. pre1 .ne. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
            do bc=1,dc2,bs
              do bd=1,dd2,bs
     
                do a=0,bcntr
                  do b=0,bcntr
                    do c=0,bcntr
                      do d=0,bcntr
                        array_out(bd+d,bc+c,bb+b,ba+a)=pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
      
              do a=da2+1,da
                do b=0,bcntr
                  do c=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,bc+c,bb+b,a)=pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
      
              do a=0,bcntr
                do b=db2+1,db
                  do c=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,bc+c,b,ba+a)=pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
      
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bd+d,c,bb+b,ba+a)=pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (dc2 .gt. 0 .and. db2 .gt. 0 .and. da2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
      
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(d,bc+c,bb+b,ba+a)=pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      !two at the ends
      if (dc2 .gt. 0 .and. dd2 .gt. 0 .and. modb .and. moda) then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
      
            do a=da2+1,da
              do b=db2+1,db
                do c=0,bcntr
                  do d=0,bcntr
                    array_out(bd+d,bc+c,b,a)=pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (db2 .gt. 0 .and. dd2 .gt. 0 .and. modc .and. moda) then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
      
            do a=da2+1,da
              do b=0,bcntr
                do c=dc2+1,dc
                  do d=0,bcntr
                    array_out(bd+d,c,bb+b,a)=pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (dc2 .gt. 0 .and. db2 .gt. 0 .and. moda .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
      
            do a=da2+1,da
              do b=0,bcntr
                do c=0,bcntr
                  do d=dd2+1,dd
                    array_out(d,bc+c,bb+b,a)=pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. dd2 .gt. 0 .and. modb .and. modc) then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
      
            do a=0,bcntr
              do b=db2+1,db
                do c=dc2+1,dc
                  do d=0,bcntr
                    array_out(bd+d,c,b,ba+a)=pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (db2 .gt. 0 .and. da2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
      
              do d=dd2+1,dd
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(d,c,bb+b,ba+a)=pre1*array_in(ba+a,bb+b,c,d)
                    enddo
                  enddo
                enddo
              enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (dc2 .gt. 0 .and. da2 .gt. 0 .and. modd .and. modb) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
      
             do d=dd2+1,dd
               do c=0,bcntr
                 do b=db2+1,db
                   do a=0,bcntr
                     array_out(d,bc+c,b,ba+a)=pre1*array_in(ba+a,b,bc+c,d)
                   enddo
                 enddo
               enddo
             enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (dd2 .gt. 0 .and. moda .and. modb .and. modc) then
        !$OMP DO
        do bd=1,dd2,bs
      
          do a=da2+1,da
            do b=db2+1,db
              do c=dc2+1,dc
                do d=0,bcntr
                  array_out(bd+d,c,b,a)=pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      if (db2 .gt. 0 .and. moda .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
      
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(d,c,bb+b,a)=pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      if (dc2 .gt. 0 .and. moda .and. modb .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
      
          do a=da2+1,da
            do b=db2+1,db
              do c=0,bcntr
                do d=dd2+1,dd
                  array_out(d,bc+c,b,a)=pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. modc .and. modd) then
        !$OMP DO
        do ba=1,da2,bs
      
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(d,c,b,ba+a)=pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. modc .and. modd) then
      !last block
        do a=da2+1,da
          do b=db2+1,db
            do c=dc2+1,dc
              do d=dd2+1,dd
                array_out(d,c,b,a)=pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif

    elseif (pre2 .ne. 0.0E0_realk .and. pre1 .eq. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
            do bc=1,dc2,bs
              do bd=1,dd2,bs
     
                do a=0,bcntr
                  do b=0,bcntr
                    do c=0,bcntr
                      do d=0,bcntr
                        array_out(bd+d,bc+c,bb+b,ba+a)=pre2*array_out(bd+d,bc+c,bb+b,ba+a)&
                                                      & + array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
      
              do a=da2+1,da
                do b=0,bcntr
                  do c=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,bc+c,bb+b,a)=pre2*array_out(bd+d,bc+c,bb+b,a)&
                                                 & + array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
      
              do a=0,bcntr
                do b=db2+1,db
                  do c=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,bc+c,b,ba+a)=pre2*array_out(bd+d,bc+c,b,ba+a)&
                                                 & + array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
      
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bd+d,c,bb+b,ba+a)=pre2*array_out(bd+d,c,bb+b,ba+a)&
                                                 & + array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (dc2 .gt. 0 .and. db2 .gt. 0 .and. da2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
      
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(d,bc+c,bb+b,ba+a)=pre2*array_out(d,bc+c,bb+b,ba+a)&
                                                 & + array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      !two at the ends
      if (dc2 .gt. 0 .and. dd2 .gt. 0 .and. modb .and. moda) then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
      
            do a=da2+1,da
              do b=db2+1,db
                do c=0,bcntr
                  do d=0,bcntr
                    array_out(bd+d,bc+c,b,a)=pre2*array_out(bd+d,bc+c,b,a)&
                                            & + array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (db2 .gt. 0 .and. dd2 .gt. 0 .and. modc .and. moda) then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
      
            do a=da2+1,da
              do b=0,bcntr
                do c=dc2+1,dc
                  do d=0,bcntr
                    array_out(bd+d,c,bb+b,a)=pre2*array_out(bd+d,c,bb+b,a)&
                                            & + array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (dc2 .gt. 0 .and. db2 .gt. 0 .and. moda .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
      
            do a=da2+1,da
              do b=0,bcntr
                do c=0,bcntr
                  do d=dd2+1,dd
                    array_out(d,bc+c,bb+b,a)=pre2*array_out(d,bc+c,bb+b,a)&
                                            & + array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. dd2 .gt. 0 .and. modb .and. modc) then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
      
            do a=0,bcntr
              do b=db2+1,db
                do c=dc2+1,dc
                  do d=0,bcntr
                    array_out(bd+d,c,b,ba+a)=pre2*array_out(bd+d,c,b,ba+a)&
                                            & + array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (db2 .gt. 0 .and. da2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
      
              do d=dd2+1,dd
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(d,c,bb+b,ba+a)=pre2*array_out(d,c,bb+b,ba+a)&
                                              & + array_in(ba+a,bb+b,c,d)
                    enddo
                  enddo
                enddo
              enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (dc2 .gt. 0 .and. da2 .gt. 0 .and. modd .and. modb) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
      
             do d=dd2+1,dd
               do c=0,bcntr
                 do b=db2+1,db
                   do a=0,bcntr
                     array_out(d,bc+c,b,ba+a)=pre2*array_out(d,bc+c,b,ba+a)&
                                             & + array_in(ba+a,b,bc+c,d)
                   enddo
                 enddo
               enddo
             enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (dd2 .gt. 0 .and. moda .and. modb .and. modc) then
        !$OMP DO
        do bd=1,dd2,bs
      
          do a=da2+1,da
            do b=db2+1,db
              do c=dc2+1,dc
                do d=0,bcntr
                  array_out(bd+d,c,b,a)=pre2*array_out(bd+d,c,b,a)&
                                       & + array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      if (db2 .gt. 0 .and. moda .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
      
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(d,c,bb+b,a)=pre2*array_out(d,c,bb+b,a)&
                                       & + array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      if (dc2 .gt. 0 .and. moda .and. modb .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
      
          do a=da2+1,da
            do b=db2+1,db
              do c=0,bcntr
                do d=dd2+1,dd
                  array_out(d,bc+c,b,a)=pre2*array_out(d,bc+c,b,a)&
                                       & + array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. modc .and. modd) then
        !$OMP DO
        do ba=1,da2,bs
      
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(d,c,b,ba+a)=pre2*array_out(d,c,b,ba+a)&
                                       & + array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. modc .and. modd) then
      !last block
        do a=da2+1,da
          do b=db2+1,db
            do c=dc2+1,dc
              do d=dd2+1,dd
                array_out(d,c,b,a)=pre2*array_out(d,c,b,a)&
                                  & + array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif

    elseif (pre2 .ne. 0.0E0_realk .and. pre1 .ne. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
            do bc=1,dc2,bs
              do bd=1,dd2,bs
     
                do a=0,bcntr
                  do b=0,bcntr
                    do c=0,bcntr
                      do d=0,bcntr
                        array_out(bd+d,bc+c,bb+b,ba+a)=pre2*array_out(bd+d,bc+c,bb+b,ba+a)&
                                                      & + pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
      
              do a=da2+1,da
                do b=0,bcntr
                  do c=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,bc+c,bb+b,a)=pre2*array_out(bd+d,bc+c,bb+b,a)&
                                                 & + pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
      
              do a=0,bcntr
                do b=db2+1,db
                  do c=0,bcntr
                    do d=0,bcntr
                      array_out(bd+d,bc+c,b,ba+a)=pre2*array_out(bd+d,bc+c,b,ba+a)&
                                                 & + pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
      
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bd+d,c,bb+b,ba+a)=pre2*array_out(bd+d,c,bb+b,ba+a)&
                                                 & + pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (dc2 .gt. 0 .and. db2 .gt. 0 .and. da2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
      
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(d,bc+c,bb+b,ba+a)=pre2*array_out(d,bc+c,bb+b,ba+a)&
                                                 & + pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
      
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      !two at the ends
      if (dc2 .gt. 0 .and. dd2 .gt. 0 .and. modb .and. moda) then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
      
            do a=da2+1,da
              do b=db2+1,db
                do c=0,bcntr
                  do d=0,bcntr
                    array_out(bd+d,bc+c,b,a)=pre2*array_out(bd+d,bc+c,b,a)&
                                            & + pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (db2 .gt. 0 .and. dd2 .gt. 0 .and. modc .and. moda) then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
      
            do a=da2+1,da
              do b=0,bcntr
                do c=dc2+1,dc
                  do d=0,bcntr
                    array_out(bd+d,c,bb+b,a)=pre2*array_out(bd+d,c,bb+b,a)&
                                            & + pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (dc2 .gt. 0 .and. db2 .gt. 0 .and. moda .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
      
            do a=da2+1,da
              do b=0,bcntr
                do c=0,bcntr
                  do d=dd2+1,dd
                    array_out(d,bc+c,bb+b,a)=pre2*array_out(d,bc+c,bb+b,a)&
                                            & + pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. dd2 .gt. 0 .and. modb .and. modc) then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
      
            do a=0,bcntr
              do b=db2+1,db
                do c=dc2+1,dc
                  do d=0,bcntr
                    array_out(bd+d,c,b,ba+a)=pre2*array_out(bd+d,c,b,ba+a)&
                                            & + pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (db2 .gt. 0 .and. da2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
      
              do d=dd2+1,dd
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(d,c,bb+b,ba+a)=pre2*array_out(d,c,bb+b,ba+a)&
                                              & + pre1*array_in(ba+a,bb+b,c,d)
                    enddo
                  enddo
                enddo
              enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (dc2 .gt. 0 .and. da2 .gt. 0 .and. modd .and. modb) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
      
             do d=dd2+1,dd
               do c=0,bcntr
                 do b=db2+1,db
                   do a=0,bcntr
                     array_out(d,bc+c,b,ba+a)=pre2*array_out(d,bc+c,b,ba+a)&
                                             & + pre1*array_in(ba+a,b,bc+c,d)
                   enddo
                 enddo
               enddo
             enddo
      
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (dd2 .gt. 0 .and. moda .and. modb .and. modc) then
        !$OMP DO
        do bd=1,dd2,bs
      
          do a=da2+1,da
            do b=db2+1,db
              do c=dc2+1,dc
                do d=0,bcntr
                  array_out(bd+d,c,b,a)=pre2*array_out(bd+d,c,b,a)&
                                       & + pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      if (db2 .gt. 0 .and. moda .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
      
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(d,c,bb+b,a)=pre2*array_out(d,c,bb+b,a)&
                                       & + pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      if (dc2 .gt. 0 .and. moda .and. modb .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
      
          do a=da2+1,da
            do b=db2+1,db
              do c=0,bcntr
                do d=dd2+1,dd
                  array_out(d,bc+c,b,a)=pre2*array_out(d,bc+c,b,a)&
                                       & + pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. modc .and. modd) then
        !$OMP DO
        do ba=1,da2,bs
      
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(d,c,b,ba+a)=pre2*array_out(d,c,b,ba+a)&
                                       & + pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
      
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. modc .and. modd) then
      !last block
        do a=da2+1,da
          do b=db2+1,db
            do c=dc2+1,dc
              do d=dd2+1,dd
                array_out(d,c,b,a)=pre2*array_out(d,c,b,a)&
                                  & + pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif

    else
      call lsquit("MANUAL REORDERING:Case not found",-1)
    endif
    ! if(pre2==0.0E0_realk)then
    !   do a=1,da
    !     do b=1,db
    !       do c=1,dc
    !         do d=1,dd
    !           if(array_out(d,c,b,a)/=pre1*array_in(a,b,c,d))then
    !             print *,"reordering (4321) not correct",a,b,c,d,da,db,dc,dd
    !             stop 0
    !           endif
    !         enddo
    !       enddo
    !     enddo
    !   enddo
    ! endif

  end subroutine manual_4321_reordering

  !\> \brief reorder a 4 diensional array  to get the indices
  !   in the order 2 1 3 4 , this is a quite expensive reordering
  !   and thus requires additional attention
  !\> \author Janus Juul Eriksen & Patrick Ettenhuber
  !\> \date November 2012
  subroutine manual_2134_reordering(bs,dims,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the original array
    integer, intent(in) :: dims(4)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1,pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3),dims(4))
    !> reordered array
    real(realk),intent(inout) :: array_out(dims(2),dims(1),dims(3),dims(4))
    integer :: a,b,c,d,ba,bb,bc,bd,da,db,dc,dd,da2,db2,dc2,dd2,bcntr
    logical :: moda,modb,modc,modd
    
    da=dims(1)
    db=dims(2)
    dc=dims(3)
    dd=dims(4)

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs
    dd2=(dd/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)
    modd=(mod(dd,bs)>0)

    bcntr=bs-1
    if (pre2 .eq. 0.0E0_realk .and. pre1 .eq. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
              do bb=1,db2,bs
     
                do d=0,bcntr
                  do c=0,bcntr
                    do a=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,ba+a,bc+c,bd+d)=array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do bb=1,db2,bs
     
              do d=0,bcntr
                do c=0,bcntr
                  do a=da2+1,da
                    do b=0,bcntr
                      array_out(bb+b,a,bc+c,bd+d)=array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
                do d=0,bcntr
                  do c=0,bcntr
                    do b=db2+1,db
                      do a=0,bcntr
                        array_out(b,ba+a,bc+c,bd+d)=array_in(ba+a,b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bb=1,db2,bs
     
              do d=0,bcntr
                do c=dc2+1,dc
                  do a=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,ba+a,c,bd+d)=array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
            do bb=1,db2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do a=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,ba+a,bc+c,d)=array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do d=0,bcntr
              do c=0,bcntr
                do a=da2+1,da
                  do b=db2+1,db
                    array_out(b,a,bc+c,bd+d)=array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do a=da2+1,da
                  do b=0,bcntr
                    array_out(bb+b,a,c,bd+d)= array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do a=da2+1,da
                  do b=0,bcntr
                    array_out(bb+b,a,bc+c,d)=array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,ba+a,c,bd+d)=array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,ba+a,bc+c,d)=array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do a=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,ba+a,c,d)=array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2 .gt. 0 .and. modb .and. modc .and. modd) then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,ba+a,c,d)=array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do a=da2+1,da
                do b=0,bcntr
                  array_out(bb+b,a,c,d)=array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do a=da2+1,da
                do b=db2+1,db
                  array_out(b,a,bc+c,d)=array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
     
          do d=0,bcntr
            do c=dc2+1,dc
              do a=da2+1,da
                do b=db2+1,db
                  array_out(b,a,c,bd+d)=array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. modc .and. modd) then
        do d=dd2+1,dd
          do c=dc2+1,dc
            do a=da2+1,da
              do b=db2+1,db 
                array_out(b,a,c,d)=array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif (pre2 .eq. 0.0E0_realk .and. pre1 .ne. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
              do bb=1,db2,bs
     
                do d=0,bcntr
                  do c=0,bcntr
                    do a=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,ba+a,bc+c,bd+d)=pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do bb=1,db2,bs
     
              do d=0,bcntr
                do c=0,bcntr
                  do a=da2+1,da
                    do b=0,bcntr
                      array_out(bb+b,a,bc+c,bd+d)=pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
                do d=0,bcntr
                  do c=0,bcntr
                    do b=db2+1,db
                      do a=0,bcntr
                        array_out(b,ba+a,bc+c,bd+d)=pre1*array_in(ba+a,b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bb=1,db2,bs
     
              do d=0,bcntr
                do c=dc2+1,dc
                  do a=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,ba+a,c,bd+d)=pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
            do bb=1,db2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do a=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,ba+a,bc+c,d)=pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do d=0,bcntr
              do c=0,bcntr
                do a=da2+1,da
                  do b=db2+1,db
                    array_out(b,a,bc+c,bd+d)=pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do a=da2+1,da
                  do b=0,bcntr
                    array_out(bb+b,a,c,bd+d)=pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do a=da2+1,da
                  do b=0,bcntr
                    array_out(bb+b,a,bc+c,d)=pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,ba+a,c,bd+d)=pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,ba+a,bc+c,d)=pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do a=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,ba+a,c,d)=pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2 .gt. 0 .and. modb .and. modc .and. modd) then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,ba+a,c,d)=pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do a=da2+1,da
                do b=0,bcntr
                  array_out(bb+b,a,c,d)=pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do a=da2+1,da
                do b=db2+1,db
                  array_out(b,a,bc+c,d)=pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
     
          do d=0,bcntr
            do c=dc2+1,dc
              do a=da2+1,da
                do b=db2+1,db
                  array_out(b,a,c,bd+d)=pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. modc .and. modd) then
        do d=dd2+1,dd
          do c=dc2+1,dc
            do a=da2+1,da
              do b=db2+1,db 
                array_out(b,a,c,d)=pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif (pre2 .ne. 0.0E0_realk .and. pre1 .eq. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
              do bb=1,db2,bs
     
                do d=0,bcntr
                  do c=0,bcntr
                    do a=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,ba+a,bc+c,bd+d)=pre2*array_out(bb+b,ba+a,bc+c,bd+d)&
                                                      & + array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do bb=1,db2,bs
     
              do d=0,bcntr
                do c=0,bcntr
                  do a=da2+1,da
                    do b=0,bcntr
                      array_out(bb+b,a,bc+c,bd+d)=pre2*array_out(bb+b,a,bc+c,bd+d)&
                                                 & + array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
                do d=0,bcntr
                  do c=0,bcntr
                    do b=db2+1,db
                      do a=0,bcntr
                        array_out(b,ba+a,bc+c,bd+d)=pre2*array_out(b,ba+a,bc+c,bd+d)&
                                                   & + array_in(ba+a,b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bb=1,db2,bs
     
              do d=0,bcntr
                do c=dc2+1,dc
                  do a=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,ba+a,c,bd+d)=pre2*array_out(bb+b,ba+a,c,bd+d)&
                                                 & + array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
            do bb=1,db2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do a=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,ba+a,bc+c,d)=pre2*array_out(bb+b,ba+a,bc+c,d)&
                                                 & + array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do d=0,bcntr
              do c=0,bcntr
                do a=da2+1,da
                  do b=db2+1,db
                    array_out(b,a,bc+c,bd+d)=pre2*array_out(b,a,bc+c,bd+d)&
                                            & + array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do a=da2+1,da
                  do b=0,bcntr
                    array_out(bb+b,a,c,bd+d)=pre2*array_out(bb+b,a,c,bd+d)&
                                            & + array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do a=da2+1,da
                  do b=0,bcntr
                    array_out(bb+b,a,bc+c,d)=pre2*array_out(bb+b,a,bc+c,d)&
                                            & + array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,ba+a,c,bd+d)=pre2*array_out(b,ba+a,c,bd+d)&
                                            & + array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,ba+a,bc+c,d)=pre2*array_out(b,ba+a,bc+c,d)&
                                            & + array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do a=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,ba+a,c,d)=pre2*array_out(bb+b,ba+a,c,d)&
                                            & + array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2 .gt. 0 .and. modb .and. modc .and. modd) then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,ba+a,c,d)=pre2*array_out(b,ba+a,c,d)&
                                       & + array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do a=da2+1,da
                do b=0,bcntr
                  array_out(bb+b,a,c,d)=pre2*array_out(bb+b,a,c,d)&
                                       & + array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do a=da2+1,da
                do b=db2+1,db
                  array_out(b,a,bc+c,d)=pre2*array_out(b,a,bc+c,d)&
                                       & + array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
     
          do d=0,bcntr
            do c=dc2+1,dc
              do a=da2+1,da
                do b=db2+1,db
                  array_out(b,a,c,bd+d)=pre2*array_out(b,a,c,bd+d)&
                                       & + array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. modc .and. modd) then
        do d=dd2+1,dd
          do c=dc2+1,dc
            do a=da2+1,da
              do b=db2+1,db 
                array_out(b,a,c,d)=pre2*array_out(b,a,c,d)&
                                  & + array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif (pre2 .ne. 0.0E0_realk .and. pre1 .ne. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
              do bb=1,db2,bs
     
                do d=0,bcntr
                  do c=0,bcntr
                    do a=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,ba+a,bc+c,bd+d)=pre2*array_out(bb+b,ba+a,bc+c,bd+d)&
                                                      & + pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do bb=1,db2,bs
     
              do d=0,bcntr
                do c=0,bcntr
                  do a=da2+1,da
                    do b=0,bcntr
                      array_out(bb+b,a,bc+c,bd+d)=pre2*array_out(bb+b,a,bc+c,bd+d)&
                                                 & + pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
                do d=0,bcntr
                  do c=0,bcntr
                    do b=db2+1,db
                      do a=0,bcntr
                        array_out(b,ba+a,bc+c,bd+d)=pre2*array_out(b,ba+a,bc+c,bd+d)&
                                                   & + pre1*array_in(ba+a,b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
            do bb=1,db2,bs
     
              do d=0,bcntr
                do c=dc2+1,dc
                  do a=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,ba+a,c,bd+d)=pre2*array_out(bb+b,ba+a,c,bd+d)&
                                                 & + pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
            do bb=1,db2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do a=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,ba+a,bc+c,d)=pre2*array_out(bb+b,ba+a,bc+c,d)&
                                                 & + pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do d=0,bcntr
              do c=0,bcntr
                do a=da2+1,da
                  do b=db2+1,db
                    array_out(b,a,bc+c,bd+d)=pre2*array_out(b,a,bc+c,bd+d)&
                                            & + pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do a=da2+1,da
                  do b=0,bcntr
                    array_out(bb+b,a,c,bd+d)=pre2*array_out(bb+b,a,c,bd+d)&
                                            & + pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do a=da2+1,da
                  do b=0,bcntr
                    array_out(bb+b,a,bc+c,d)=pre2*array_out(bb+b,a,bc+c,d)&
                                            & + pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,ba+a,c,bd+d)=pre2*array_out(b,ba+a,c,bd+d)&
                                            & + pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,ba+a,bc+c,d)=pre2*array_out(b,ba+a,bc+c,d)&
                                            & + pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do a=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,ba+a,c,d)=pre2*array_out(bb+b,ba+a,c,d)&
                                            & + pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if(da2 .gt. 0 .and. modb .and. modc .and. modd) then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,ba+a,c,d)=pre2*array_out(b,ba+a,c,d)&
                                       & + pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do a=da2+1,da
                do b=0,bcntr
                  array_out(bb+b,a,c,d)=pre2*array_out(bb+b,a,c,d)&
                                       & + pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
     
          do d=dd2+1,dd
            do c=0,bcntr
              do a=da2+1,da
                do b=db2+1,db
                  array_out(b,a,bc+c,d)=pre2*array_out(b,a,bc+c,d)&
                                       & + pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
     
          do d=0,bcntr
            do c=dc2+1,dc
              do a=da2+1,da
                do b=db2+1,db
                  array_out(b,a,c,bd+d)=pre2*array_out(b,a,c,bd+d)&
                                       & + pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. modc .and. modd) then
        do d=dd2+1,dd
          do c=dc2+1,dc
            do a=da2+1,da
              do b=db2+1,db 
                array_out(b,a,c,d)=pre2*array_out(b,a,c,d)&
                                  & + pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif

    else
      call lsquit("MANUAL REORDERING:Case not found",-1)
    endif

    ! print *,"check 2134"
    ! if(pre2==0.0E0_realk)then
    !   do d=1,dd
    !     do c=1,dc
    !       do a=1,da
    !         do b=1,db
    !           if(array_out(b,a,c,d)/=pre1*array_in(a,b,c,d))then
    !             print *,"2134 reordering not correct",a,b,c,d,da,db,dc,dd
    !             stop 0
    !           endif
    !         enddo
    !       enddo
    !     enddo
    !   enddo
    ! endif
  end subroutine manual_2134_reordering

  !\> \brief reorder a 4 diensional array  to get the indices
  !   in the order 2 4 3 1 , this is a quite expensive reordering
  !   and thus requires additional attention
  !\> \author Janus Juul Eriksen & Patrick Ettenhuber
  !\> \date November 2012
  subroutine manual_2431_reordering(bs,dims,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the original array
    integer, intent(in) :: dims(4)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1,pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3),dims(4))
    !> reordered array
    real(realk),intent(inout) :: array_out(dims(2),dims(4),dims(3),dims(1))
    integer :: a,b,c,d,ba,bb,bc,bd,da,db,dc,dd,da2,db2,dc2,dd2,bcntr
    logical :: moda,modb,modc,modd
    
    da=dims(1)
    db=dims(2)
    dc=dims(3)
    dd=dims(4)

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs
    dd2=(dd/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)
    modd=(mod(dd,bs)>0)

    bcntr=bs-1

    if (pre2 .eq. 0.0E0_realk .and. pre1 .eq. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
              do bb=1,db2,bs
     
                do a=0,bcntr
                  do c=0,bcntr
                    do d=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,bd+d,bc+c,ba+a)=array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
     
              do a=da2+1,da
                do c=0,bcntr
                  do d=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bd+d,bc+c,a)=array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do a=0,bcntr
                        array_out(b,bd+d,bc+c,ba+a)=array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
     
              do a=0,bcntr
                do c=dc2+1,dc
                  do d=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bd+d,c,ba+a)=array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bb+b,d,bc+c,ba+a)=array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do a=da2+1,da
              do c=0,bcntr
                do d=0,bcntr
                  do b=db2+1,db
                    array_out(b,bd+d,bc+c,a)=array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do a=da2+1,da
              do c=dc2+1,dc
                do d=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,bd+d,c,a)=array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do a=da2+1,da
              do c=0,bcntr
                do d=dd2+1,dd
                  do b=0,bcntr
                    array_out(bb+b,d,bc+c,a)=array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,bd+d,c,ba+a)=array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,d,bc+c,ba+a)=array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(bb+b,d,c,ba+a)=array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. modc .and. modd) then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,d,c,ba+a)=array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
     
          do a=da2+1,da
            do c=dc2+1,dc
              do d=dd2+1,dd
                do b=0,bcntr
                  array_out(bb+b,d,c,a)=array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
     
          do a=da2+1,da
            do c=0,bcntr
              do d=dd2+1,dd
                do b=db2+1,db
                  array_out(b,d,bc+c,a)=array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
     
          do a=da2+1,da
            do c=dc2+1,dc
              do d=0,bcntr
                do b=db2+1,db
                  array_out(b,bd+d,c,a)=array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. modc .and. modd) then
        do a=da2+1,da
          do c=dc2+1,dc
            do d=dd2+1,dd
              do b=db2+1,db
                array_out(b,d,c,a)=array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif (pre2 .eq. 0.0E0_realk .and. pre1 .ne. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
              do bb=1,db2,bs
     
                do a=0,bcntr
                  do c=0,bcntr
                    do d=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,bd+d,bc+c,ba+a)=pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
     
              do a=da2+1,da
                do c=0,bcntr
                  do d=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bd+d,bc+c,a)=pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do a=0,bcntr
                        array_out(b,bd+d,bc+c,ba+a)=pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
     
              do a=0,bcntr
                do c=dc2+1,dc
                  do d=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bd+d,c,ba+a)=pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bb+b,d,bc+c,ba+a)=pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do a=da2+1,da
              do c=0,bcntr
                do d=0,bcntr
                  do b=db2+1,db
                    array_out(b,bd+d,bc+c,a)=pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do a=da2+1,da
              do c=dc2+1,dc
                do d=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,bd+d,c,a)=pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do a=da2+1,da
              do c=0,bcntr
                do d=dd2+1,dd
                  do b=0,bcntr
                    array_out(bb+b,d,bc+c,a)=pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,bd+d,c,ba+a)=pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,d,bc+c,ba+a)=pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(bb+b,d,c,ba+a)=pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. modc .and. modd) then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,d,c,ba+a)=pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
     
          do a=da2+1,da
            do c=dc2+1,dc
              do d=dd2+1,dd
                do b=0,bcntr
                  array_out(bb+b,d,c,a)=pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
     
          do a=da2+1,da
            do c=0,bcntr
              do d=dd2+1,dd
                do b=db2+1,db
                  array_out(b,d,bc+c,a)=pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
     
          do a=da2+1,da
            do c=dc2+1,dc
              do d=0,bcntr
                do b=db2+1,db
                  array_out(b,bd+d,c,a)=pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. modc .and. modd) then
        do a=da2+1,da
          do c=dc2+1,dc
            do d=dd2+1,dd
              do b=db2+1,db
                array_out(b,d,c,a)=pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif (pre2 .ne. 0.0E0_realk .and. pre1 .eq. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
              do bb=1,db2,bs
     
                do a=0,bcntr
                  do c=0,bcntr
                    do d=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,bd+d,bc+c,ba+a)=pre2*array_out(bb+b,bd+d,bc+c,ba+a)&
                                                      & + array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
     
              do a=da2+1,da
                do c=0,bcntr
                  do d=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bd+d,bc+c,a)=pre2*array_out(bb+b,bd+d,bc+c,a)&
                                                 & + array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do a=0,bcntr
                        array_out(b,bd+d,bc+c,ba+a)=pre2*array_out(b,bd+d,bc+c,ba+a)&
                                                   & + array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
     
              do a=0,bcntr
                do c=dc2+1,dc
                  do d=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bd+d,c,ba+a)=pre2*array_out(bb+b,bd+d,c,ba+a)&
                                                 & + array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bb+b,d,bc+c,ba+a)=pre2*array_out(bb+b,d,bc+c,ba+a)&
                                                 & + array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do a=da2+1,da
              do c=0,bcntr
                do d=0,bcntr
                  do b=db2+1,db
                    array_out(b,bd+d,bc+c,a)=pre2*array_out(b,bd+d,bc+c,a)&
                                            & + array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do a=da2+1,da
              do c=dc2+1,dc
                do d=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,bd+d,c,a)=pre2*array_out(bb+b,bd+d,c,a)&
                                            & + array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do a=da2+1,da
              do c=0,bcntr
                do d=dd2+1,dd
                  do b=0,bcntr
                    array_out(bb+b,d,bc+c,a)=pre2*array_out(bb+b,d,bc+c,a)&
                                            & + array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,bd+d,c,ba+a)=pre2*array_out(b,bd+d,c,ba+a)&
                                            & + array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,d,bc+c,ba+a)=pre2*array_out(b,d,bc+c,ba+a)&
                                            & + array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(bb+b,d,c,ba+a)=pre2*array_out(bb+b,d,c,ba+a)&
                                            & + array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. modc .and. modd) then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,d,c,ba+a)=pre2*array_out(b,d,c,ba+a)&
                                       & + array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
     
          do a=da2+1,da
            do c=dc2+1,dc
              do d=dd2+1,dd
                do b=0,bcntr
                  array_out(bb+b,d,c,a)=pre2*array_out(bb+b,d,c,a)&
                                       & + array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
     
          do a=da2+1,da
            do c=0,bcntr
              do d=dd2+1,dd
                do b=db2+1,db
                  array_out(b,d,bc+c,a)=pre2*array_out(b,d,bc+c,a)&
                                       & + array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
     
          do a=da2+1,da
            do c=dc2+1,dc
              do d=0,bcntr
                do b=db2+1,db
                  array_out(b,bd+d,c,a)=pre2*array_out(b,bd+d,c,a)&
                                       & + array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. modc .and. modd) then
        do a=da2+1,da
          do c=dc2+1,dc
            do d=dd2+1,dd
              do b=db2+1,db
                array_out(b,d,c,a)=pre2*array_out(b,d,c,a)&
                                  & + array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif (pre2 .ne. 0.0E0_realk .and. pre1 .ne. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bd=1,dd2,bs
              do bb=1,db2,bs
     
                do a=0,bcntr
                  do c=0,bcntr
                    do d=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,bd+d,bc+c,ba+a)=pre2*array_out(bb+b,bd+d,bc+c,ba+a)&
                                                      & + pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
     
              do a=da2+1,da
                do c=0,bcntr
                  do d=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bd+d,bc+c,a)=pre2*array_out(bb+b,bd+d,bc+c,a)&
                                                 & + pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do a=0,bcntr
                        array_out(b,bd+d,bc+c,ba+a)=pre2*array_out(b,bd+d,bc+c,ba+a)&
                                                   & + pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
            do bb=1,db2,bs
     
              do a=0,bcntr
                do c=dc2+1,dc
                  do d=0,bcntr
                    do b=0,bcntr
                      array_out(bb+b,bd+d,c,ba+a)=pre2*array_out(bb+b,bd+d,c,ba+a)&
                                                 & + pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bb+b,d,bc+c,ba+a)=pre2*array_out(bb+b,d,bc+c,ba+a)&
                                                 & + pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do bd=1,dd2,bs
     
            do a=da2+1,da
              do c=0,bcntr
                do d=0,bcntr
                  do b=db2+1,db
                    array_out(b,bd+d,bc+c,a)=pre2*array_out(b,bd+d,bc+c,a)&
                                            & + pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
     
            do a=da2+1,da
              do c=dc2+1,dc
                do d=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,bd+d,c,a)=pre2*array_out(bb+b,bd+d,c,a)&
                                            & + pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do a=da2+1,da
              do c=0,bcntr
                do d=dd2+1,dd
                  do b=0,bcntr
                    array_out(bb+b,d,bc+c,a)=pre2*array_out(bb+b,d,bc+c,a)&
                                            & + pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,bd+d,c,ba+a)=pre2*array_out(b,bd+d,c,ba+a)&
                                            & + pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(b,d,bc+c,ba+a)=pre2*array_out(b,d,bc+c,ba+a)&
                                            & + pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(bb+b,d,c,ba+a)=pre2*array_out(bb+b,d,c,ba+a)&
                                            & + pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. modc .and. modd) then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,d,c,ba+a)=pre2*array_out(b,d,c,ba+a)&
                                       & + pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
     
          do a=da2+1,da
            do c=dc2+1,dc
              do d=dd2+1,dd
                do b=0,bcntr
                  array_out(bb+b,d,c,a)=pre2*array_out(bb+b,d,c,a)&
                                       & + pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
     
          do a=da2+1,da
            do c=0,bcntr
              do d=dd2+1,dd
                do b=db2+1,db
                  array_out(b,d,bc+c,a)=pre2*array_out(b,d,bc+c,a)&
                                       & + pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
     
          do a=da2+1,da
            do c=dc2+1,dc
              do d=0,bcntr
                do b=db2+1,db
                  array_out(b,bd+d,c,a)=pre2*array_out(b,bd+d,c,a)&
                                       & + pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. modc .and. modd) then
        do a=da2+1,da
          do c=dc2+1,dc
            do d=dd2+1,dd
              do b=db2+1,db
                array_out(b,d,c,a)=pre2*array_out(b,d,c,a)&
                                  & + pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif

    else
      call lsquit("MANUAL REORDERING:Case not found",-1)
    endif

   ! print *,"checking the result 2431"
   ! if(pre2==0.0E0_realk)then
   !   do a=1,da
   !     do c=1,dc
   !       do d=1,dd
   !         do b=1,db
   !           if(array_out(b,d,c,a)/=pre1*array_in(a,b,c,d))then
   !             print *,"2431 reordering not correct",a,b,c,d,da,db,dc,dd
   !             stop 0
   !           endif
   !         enddo
   !       enddo
   !     enddo
   !   enddo
   ! endif
  end subroutine manual_2431_reordering

  !\> \brief reorder a 4 diensional array  to get the indices
  !   in the order 3 4 2 1 , this is a quite expensive reordering
  !   and thus requires additional attention
  !\> \author Janus Juul Eriksen & Patrick Ettenhuber
  !\> \date November 2012
  subroutine manual_3421_reordering(bs,dims,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the original array
    integer, intent(in) :: dims(4)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1,pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3),dims(4))
    !> reordered array
    real(realk),intent(inout) :: array_out(dims(3),dims(4),dims(2),dims(1))
    integer :: a,b,c,d,ba,bb,bc,bd,da,db,dc,dd,da2,db2,dc2,dd2,bcntr
    logical :: moda,modb,modc,modd
    
    da=dims(1)
    db=dims(2)
    dc=dims(3)
    dd=dims(4)

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs
    dd2=(dd/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)
    modd=(mod(dd,bs)>0)

    bcntr=bs-1

    if (pre2 .eq. 0.0E0_realk .and. pre1 .eq. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
            do bd=1,dd2,bs
              do bc=1,dc2,bs
     
                do a=0,bcntr
                  do b=0,bcntr
                    do d=0,bcntr
                      do c=0,bcntr
                        array_out(bc+c,bd+d,bb+b,ba+a)=array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do bc=1,dc2,bs
     
              do a=da2+1,da
                do b=0,bcntr
                  do d=0,bcntr
                    do c=0,bcntr
                      array_out(bc+c,bd+d,bb+b,a)=array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
            do bc=1,dc2,bs
     
              do a=0,bcntr
                do b=db2+1,db
                  do d=0,bcntr
                    do c=0,bcntr
                        array_out(bc+c,bd+d,b,ba+a)=array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(c,bd+d,bb+b,ba+a)=array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bc+c,d,bb+b,ba+a)=array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do a=da2+1,da
              do b=db2+1,db
                do d=0,bcntr
                  do c=0,bcntr
                    array_out(bc+c,bd+d,b,a)=array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
     
            do a=da2+1,da
              do b=0,bcntr
                do d=0,bcntr
                  do c=dc2+1,dc
                    array_out(c,bd+d,bb+b,a)=array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
     
            do a=da2+1,da
              do b=0,bcntr
                do d=dd2+1,dd
                  do c=0,bcntr
                    array_out(bc+c,d,bb+b,a)=array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(c,bd+d,b,ba+a)=array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(bc+c,d,b,ba+a)=array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(c,d,bb+b,ba+a)=array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. modc .and. modd) then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(c,d,b,ba+a)=array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(c,d,bb+b,a)=array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
     
          do a=da2+1,da
            do b=db2+1,db
              do d=dd2+1,dd
                do c=0,bcntr
                  array_out(bc+c,d,b,a)=array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
     
          do a=da2+1,da
            do b=db2+1,db
              do d=0,bcntr
                do c=dc2+1,dc
                  array_out(c,bd+d,b,a)=array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. modc .and. modd) then
        do a=da2+1,da
          do b=db2+1,db
            do d=dd2+1,dd
              do c=dc2+1,dc
                array_out(c,d,b,a)=array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif (pre2 .eq. 0.0E0_realk .and. pre1 .ne. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
            do bd=1,dd2,bs
              do bc=1,dc2,bs
     
                do a=0,bcntr
                  do b=0,bcntr
                    do d=0,bcntr
                      do c=0,bcntr
                        array_out(bc+c,bd+d,bb+b,ba+a)=pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do bc=1,dc2,bs
     
              do a=da2+1,da
                do b=0,bcntr
                  do d=0,bcntr
                    do c=0,bcntr
                      array_out(bc+c,bd+d,bb+b,a)=pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
            do bc=1,dc2,bs
     
              do a=0,bcntr
                do b=db2+1,db
                  do d=0,bcntr
                    do c=0,bcntr
                        array_out(bc+c,bd+d,b,ba+a)=pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(c,bd+d,bb+b,ba+a)=pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bc+c,d,bb+b,ba+a)=pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do a=da2+1,da
              do b=db2+1,db
                do d=0,bcntr
                  do c=0,bcntr
                    array_out(bc+c,bd+d,b,a)=pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
     
            do a=da2+1,da
              do b=0,bcntr
                do d=0,bcntr
                  do c=dc2+1,dc
                    array_out(c,bd+d,bb+b,a)=pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
     
            do a=da2+1,da
              do b=0,bcntr
                do d=dd2+1,dd
                  do c=0,bcntr
                    array_out(bc+c,d,bb+b,a)=pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(c,bd+d,b,ba+a)=pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(bc+c,d,b,ba+a)=pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(c,d,bb+b,ba+a)=pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. modc .and. modd) then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(c,d,b,ba+a)=pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(c,d,bb+b,a)=pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
     
          do a=da2+1,da
            do b=db2+1,db
              do d=dd2+1,dd
                do c=0,bcntr
                  array_out(bc+c,d,b,a)=pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
     
          do a=da2+1,da
            do b=db2+1,db
              do d=0,bcntr
                do c=dc2+1,dc
                  array_out(c,bd+d,b,a)=pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. modc .and. modd) then
        do a=da2+1,da
          do b=db2+1,db
            do d=dd2+1,dd
              do c=dc2+1,dc
                array_out(c,d,b,a)=pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif (pre2 .ne. 0.0E0_realk .and. pre1 .eq. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
            do bd=1,dd2,bs
              do bc=1,dc2,bs
     
                do a=0,bcntr
                  do b=0,bcntr
                    do d=0,bcntr
                      do c=0,bcntr
                        array_out(bc+c,bd+d,bb+b,ba+a)=pre2*array_out(bc+c,bd+d,bb+b,ba+a)&
                                                      & + array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do bc=1,dc2,bs
     
              do a=da2+1,da
                do b=0,bcntr
                  do d=0,bcntr
                    do c=0,bcntr
                      array_out(bc+c,bd+d,bb+b,a)=pre2*array_out(bc+c,bd+d,bb+b,a)&
                                                 & + array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
            do bc=1,dc2,bs
     
              do a=0,bcntr
                do b=db2+1,db
                  do d=0,bcntr
                    do c=0,bcntr
                        array_out(bc+c,bd+d,b,ba+a)=pre2*array_out(bc+c,bd+d,b,ba+a)&
                                                   & + array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(c,bd+d,bb+b,ba+a)=pre2*array_out(c,bd+d,bb+b,ba+a)&
                                                 & + array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bc+c,d,bb+b,ba+a)=pre2*array_out(bc+c,d,bb+b,ba+a)&
                                                 & + array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do a=da2+1,da
              do b=db2+1,db
                do d=0,bcntr
                  do c=0,bcntr
                    array_out(bc+c,bd+d,b,a)=pre2*array_out(bc+c,bd+d,b,a)&
                                            & + array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
     
            do a=da2+1,da
              do b=0,bcntr
                do d=0,bcntr
                  do c=dc2+1,dc
                    array_out(c,bd+d,bb+b,a)=pre2*array_out(c,bd+d,bb+b,a)&
                                            & + array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
     
            do a=da2+1,da
              do b=0,bcntr
                do d=dd2+1,dd
                  do c=0,bcntr
                    array_out(bc+c,d,bb+b,a)=pre2*array_out(bc+c,d,bb+b,a)&
                                            & + array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(c,bd+d,b,ba+a)=pre2*array_out(c,bd+d,b,ba+a)&
                                            & + array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(bc+c,d,b,ba+a)=pre2*array_out(bc+c,d,b,ba+a)&
                                            & + array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(c,d,bb+b,ba+a)=pre2*array_out(c,d,bb+b,ba+a)&
                                            & + array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. modc .and. modd) then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(c,d,b,ba+a)=pre2*array_out(c,d,b,ba+a)&
                                       & + array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(c,d,bb+b,a)=pre2*array_out(c,d,bb+b,a)&
                                       & + array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
     
          do a=da2+1,da
            do b=db2+1,db
              do d=dd2+1,dd
                do c=0,bcntr
                  array_out(bc+c,d,b,a)=pre2*array_out(bc+c,d,b,a)&
                                       & + array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
     
          do a=da2+1,da
            do b=db2+1,db
              do d=0,bcntr
                do c=dc2+1,dc
                  array_out(c,bd+d,b,a)=pre2*array_out(c,bd+d,b,a)&
                                       & + array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. modc .and. modd) then
        do a=da2+1,da
          do b=db2+1,db
            do d=dd2+1,dd
              do c=dc2+1,dc
                array_out(c,d,b,a)=pre2*array_out(c,d,b,a)&
                                  & + array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif (pre2 .ne. 0.0E0_realk .and. pre1 .ne. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
            do bd=1,dd2,bs
              do bc=1,dc2,bs
     
                do a=0,bcntr
                  do b=0,bcntr
                    do d=0,bcntr
                      do c=0,bcntr
                        array_out(bc+c,bd+d,bb+b,ba+a)=pre2*array_out(bc+c,bd+d,bb+b,ba+a)&
                                                      & + pre1*array_in(ba+a,bb+b,bc+c,bd+d)
                      enddo
                    enddo
                  enddo
                enddo
     
              enddo
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
            do bc=1,dc2,bs
     
              do a=da2+1,da
                do b=0,bcntr
                  do d=0,bcntr
                    do c=0,bcntr
                      array_out(bc+c,bd+d,bb+b,a)=pre2*array_out(bc+c,bd+d,bb+b,a)&
                                                 & + pre1*array_in(a,bb+b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bd=1,dd2,bs
            do bc=1,dc2,bs
     
              do a=0,bcntr
                do b=db2+1,db
                  do d=0,bcntr
                    do c=0,bcntr
                        array_out(bc+c,bd+d,b,ba+a)=pre2*array_out(bc+c,bd+d,b,ba+a)&
                                                   & + pre1*array_in(ba+a,b,bc+c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=0,bcntr
                do c=dc2+1,dc
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(c,bd+d,bb+b,ba+a)=pre2*array_out(c,bd+d,bb+b,ba+a)&
                                                 & + pre1*array_in(ba+a,bb+b,c,bd+d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bc+c,d,bb+b,ba+a)=pre2*array_out(bc+c,d,bb+b,ba+a)&
                                                 & + pre1*array_in(ba+a,bb+b,bc+c,d)
                    enddo
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. dc2 .gt. 0 .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do bc=1,dc2,bs
     
            do a=da2+1,da
              do b=db2+1,db
                do d=0,bcntr
                  do c=0,bcntr
                    array_out(bc+c,bd+d,b,a)=pre2*array_out(bc+c,bd+d,b,a)&
                                            & + pre1*array_in(a,b,bc+c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bd=1,dd2,bs
     
            do a=da2+1,da
              do b=0,bcntr
                do d=0,bcntr
                  do c=dc2+1,dc
                    array_out(c,bd+d,bb+b,a)=pre2*array_out(c,bd+d,bb+b,a)&
                                            & + pre1*array_in(a,bb+b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
     
            do a=da2+1,da
              do b=0,bcntr
                do d=dd2+1,dd
                  do c=0,bcntr
                    array_out(bc+c,d,bb+b,a)=pre2*array_out(bc+c,d,bb+b,a)&
                                            & + pre1*array_in(a,bb+b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
          do ba=1,da2,bs
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(c,bd+d,b,ba+a)=pre2*array_out(c,bd+d,b,ba+a)&
                                            & + pre1*array_in(ba+a,b,c,bd+d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(bc+c,d,b,ba+a)=pre2*array_out(bc+c,d,b,ba+a)&
                                            & + pre1*array_in(ba+a,b,bc+c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
     
            do d=dd2+1,dd
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(c,d,bb+b,ba+a)=pre2*array_out(c,d,bb+b,ba+a)&
                                            & + pre1*array_in(ba+a,bb+b,c,d)
                  enddo
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. modc .and. modd) then
        !$OMP DO
        do ba=1,da2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(c,d,b,ba+a)=pre2*array_out(c,d,b,ba+a)&
                                       & + pre1*array_in(ba+a,b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. modc .and. modd) then
        !$OMP DO
        do bb=1,db2,bs
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=0,bcntr
                do a=da2+1,da
                  array_out(c,d,bb+b,a)=pre2*array_out(c,d,bb+b,a)&
                                       & + pre1*array_in(a,bb+b,c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. dc2 .gt. 0 .and. modd) then
        !$OMP DO
        do bc=1,dc2,bs
     
          do a=da2+1,da
            do b=db2+1,db
              do d=dd2+1,dd
                do c=0,bcntr
                  array_out(bc+c,d,b,a)=pre2*array_out(bc+c,d,b,a)&
                                       & + pre1*array_in(a,b,bc+c,d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. modb .and. modc .and. dd2 .gt. 0) then
        !$OMP DO
        do bd=1,dd2,bs
     
          do a=da2+1,da
            do b=db2+1,db
              do d=0,bcntr
                do c=dc2+1,dc
                  array_out(c,bd+d,b,a)=pre2*array_out(c,bd+d,b,a)&
                                       & + pre1*array_in(a,b,c,bd+d)
                enddo
              enddo
            enddo
          enddo
     
        enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. modc .and. modd) then
        do a=da2+1,da
          do b=db2+1,db
            do d=dd2+1,dd
              do c=dc2+1,dc
                array_out(c,d,b,a)=pre2*array_out(c,d,b,a)&
                                  & + pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif

    else
      call lsquit("MANUAL REORDERING:Case not found",-1)
    endif
   
   ! print *,"checking the result 3421"
   ! if(pre2==0.0E0_realk)then
   !   do a=1,da
   !     do b=1,db
   !       do d=1,dd
   !         do c=1,dc
   !           if(array_out(c,d,b,a)/=pre1*array_in(a,b,c,d))then
   !             print *,"3421 reordering not correct",a,b,c,d,da,db,dc,dd
   !             stop 0
   !           endif
   !         enddo
   !       enddo
   !     enddo
   !   enddo
   ! endif
  end subroutine manual_3421_reordering
  !
  ! *******************
  ! **** 3d arrays ****
  ! *******************
  !

  !\> \brief reorder a 3 diensional array  to get the indices
  !   in the order 3 2 1 , this is a quite expensive reordering
  !   and thus requires additional attention
  !\> \author Janus Juul Eriksen & Patrick Ettenhuber
  !\> \date November 2012
  subroutine manual_321_reordering(bs,dims,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the original array
    integer, intent(in) :: dims(3)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1,pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3))
    !> reordered array
    real(realk),intent(inout) :: array_out(dims(3),dims(2),dims(1))
    integer :: a,b,c,ba,bb,bc,da,db,dc,da2,db2,dc2,bcntr
    logical :: moda,modb,modc
    
    da=dims(1)
    db=dims(2)
    dc=dims(3)

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)

    bcntr=bs-1

    if (pre2 .eq. 0.0E0_realk .and. pre1 .eq. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,ba,bb,bc),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,da2,db2,dc2,bcntr,moda,modb,modc,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
            do bc=1,dc2,bs
     
              do a=0,bcntr
                do b=0,bcntr
                  do c=0,bcntr
                    array_out(bc+c,bb+b,ba+a)=array_in(ba+a,bb+b,bc+c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
     
            do a=da2+1,da
              do b=0,bcntr
                do c=0,bcntr
                  array_out(bc+c,bb+b,a)=array_in(a,bb+b,bc+c)
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
   
            do a=0,bcntr
              do b=db2+1,db
                do c=0,bcntr
                  array_out(bc+c,b,ba+a)=array_in(ba+a,b,bc+c)
                enddo
              enddo
            enddo
   
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc) then
        !$OMP DO
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(c,bb+b,ba+a)=array_in(ba+a,bb+b,c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. dc2 .gt. 0) then
        do bc=1,dc2,bs
     
          do a=da2+1,da
            do b=db2+1,db
              do c=0,bcntr
                array_out(bc+c,b,a)=array_in(a,b,bc+c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. db2 .gt. 0 .and. modc) then
        do bb=1,db2,bs
     
          do a=da2+1,da
            do b=0,bcntr
              do c=dc2+1,dc
                array_out(c,bb+b,a)=array_in(a,bb+b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (da2 .gt. 0 .and. modb .and. modc) then
        do ba=1,da2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do a=0,bcntr
                array_out(c,b,ba+a)=array_in(ba+a,b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. modb .and. modc) then
        do a=da2+1,da
          do b=db2+1,db
            do c=dc2+1,dc
              array_out(c,b,a)=array_in(a,b,c)
            enddo
          enddo
        enddo
      endif
    elseif (pre2 .eq. 0.0E0_realk .and. pre1 .ne. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,ba,bb,bc),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,da2,db2,dc2,bcntr,moda,modb,modc,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
            do bc=1,dc2,bs
     
              do a=0,bcntr
                do b=0,bcntr
                  do c=0,bcntr
                    array_out(bc+c,bb+b,ba+a)=pre1*array_in(ba+a,bb+b,bc+c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
     
            do a=da2+1,da
              do b=0,bcntr
                do c=0,bcntr
                  array_out(bc+c,bb+b,a)=pre1*array_in(a,bb+b,bc+c)
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
   
            do a=0,bcntr
              do b=db2+1,db
                do c=0,bcntr
                  array_out(bc+c,b,ba+a)=pre1*array_in(ba+a,b,bc+c)
                enddo
              enddo
            enddo
   
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc) then
        !$OMP DO
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(c,bb+b,ba+a)=pre1*array_in(ba+a,bb+b,c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. dc2 .gt. 0) then
        do bc=1,dc2,bs
     
          do a=da2+1,da
            do b=db2+1,db
              do c=0,bcntr
                array_out(bc+c,b,a)=pre1*array_in(a,b,bc+c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. db2 .gt. 0 .and. modc) then
        do bb=1,db2,bs
     
          do a=da2+1,da
            do b=0,bcntr
              do c=dc2+1,dc
                array_out(c,bb+b,a)=pre1*array_in(a,bb+b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (da2 .gt. 0 .and. modb .and. modc) then
        do ba=1,da2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do a=0,bcntr
                array_out(c,b,ba+a)=pre1*array_in(ba+a,b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. modb .and. modc) then
        do a=da2+1,da
          do b=db2+1,db
            do c=dc2+1,dc
              array_out(c,b,a)=pre1*array_in(a,b,c)
            enddo
          enddo
        enddo
      endif
    elseif (pre2 .ne. 0.0E0_realk .and. pre1 .eq. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,ba,bb,bc),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,da2,db2,dc2,bcntr,moda,modb,modc,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
            do bc=1,dc2,bs
     
              do a=0,bcntr
                do b=0,bcntr
                  do c=0,bcntr
                    array_out(bc+c,bb+b,ba+a)=pre2*array_out(bc+c,bb+b,ba+a)&
                                             & + array_in(ba+a,bb+b,bc+c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
     
            do a=da2+1,da
              do b=0,bcntr
                do c=0,bcntr
                  array_out(bc+c,bb+b,a)=pre2*array_out(bc+c,bb+b,a)&
                                        & + array_in(a,bb+b,bc+c)
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
   
            do a=0,bcntr
              do b=db2+1,db
                do c=0,bcntr
                  array_out(bc+c,b,ba+a)=pre2*array_out(bc+c,b,ba+a)&
                                        & + array_in(ba+a,b,bc+c)
                enddo
              enddo
            enddo
   
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc) then
        !$OMP DO
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(c,bb+b,ba+a)=pre2*array_out(c,bb+b,ba+a)&
                                          & + array_in(ba+a,bb+b,c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. dc2 .gt. 0) then
        do bc=1,dc2,bs
     
          do a=da2+1,da
            do b=db2+1,db
              do c=0,bcntr
                array_out(bc+c,b,a)=pre2*array_out(bc+c,b,a)&
                                   & + array_in(a,b,bc+c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. db2 .gt. 0 .and. modc) then
        do bb=1,db2,bs
     
          do a=da2+1,da
            do b=0,bcntr
              do c=dc2+1,dc
                array_out(c,bb+b,a)=pre2*array_out(c,bb+b,a)&
                                   & + array_in(a,bb+b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (da2 .gt. 0 .and. modb .and. modc) then
        do ba=1,da2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do a=0,bcntr
                array_out(c,b,ba+a)=pre2*array_out(c,b,ba+a)&
                                   & + array_in(ba+a,b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. modb .and. modc) then
        do a=da2+1,da
          do b=db2+1,db
            do c=dc2+1,dc
              array_out(c,b,a)=pre2*array_out(c,b,a)&
                              & + array_in(a,b,c)
            enddo
          enddo
        enddo
      endif
    elseif (pre2 .ne. 0.0E0_realk .and. pre1 .ne. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,ba,bb,bc),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,da2,db2,dc2,bcntr,moda,modb,modc,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bb=1,db2,bs
            do bc=1,dc2,bs
     
              do a=0,bcntr
                do b=0,bcntr
                  do c=0,bcntr
                    array_out(bc+c,bb+b,ba+a)=pre2*array_out(bc+c,bb+b,ba+a)&
                                             & + pre1*array_in(ba+a,bb+b,bc+c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
     
            do a=da2+1,da
              do b=0,bcntr
                do c=0,bcntr
                  array_out(bc+c,bb+b,a)=pre2*array_out(bc+c,bb+b,a)&
                                        & + pre1*array_in(a,bb+b,bc+c)
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
   
            do a=0,bcntr
              do b=db2+1,db
                do c=0,bcntr
                  array_out(bc+c,b,ba+a)=pre2*array_out(bc+c,b,ba+a)&
                                        & + pre1*array_in(ba+a,b,bc+c)
                enddo
              enddo
            enddo
   
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc) then
        !$OMP DO
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(c,bb+b,ba+a)=pre2*array_out(c,bb+b,ba+a)&
                                          & + pre1*array_in(ba+a,bb+b,c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. dc2 .gt. 0) then
        do bc=1,dc2,bs
     
          do a=da2+1,da
            do b=db2+1,db
              do c=0,bcntr
                array_out(bc+c,b,a)=pre2*array_out(bc+c,b,a)&
                                   & + pre1*array_in(a,b,bc+c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. db2 .gt. 0 .and. modc) then
        do bb=1,db2,bs
     
          do a=da2+1,da
            do b=0,bcntr
              do c=dc2+1,dc
                array_out(c,bb+b,a)=pre2*array_out(c,bb+b,a)&
                                   & + pre1*array_in(a,bb+b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (da2 .gt. 0 .and. modb .and. modc) then
        do ba=1,da2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do a=0,bcntr
                array_out(c,b,ba+a)=pre2*array_out(c,b,ba+a)&
                                   & + pre1*array_in(ba+a,b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. modb .and. modc) then
        do a=da2+1,da
          do b=db2+1,db
            do c=dc2+1,dc
              array_out(c,b,a)=pre2*array_out(c,b,a)&
                              & + pre1*array_in(a,b,c)
            enddo
          enddo
        enddo
      endif

    else
      call lsquit("MANUAL REORDERING:Case not found",-1)
    endif
   
   ! print *,"checking the result 321"
   ! if(pre2==0.0E0_realk)then
   !   do a=1,da
   !     do b=1,db
   !       do c=1,dc
   !         if(array_out(c,b,a)/=pre1*array_in(a,b,c))then
   !           print *,"321 reordering not correct",a,b,c,da,db,dc
   !           stop 0
   !         endif
   !       enddo
   !     enddo
   !   enddo
   ! endif
  end subroutine manual_321_reordering

  !\> \brief reorder a 3 diensional array  to get the indices
  !   in the order 3 1 2 , this is a quite expensive reordering
  !   and thus requires additional attention
  !\> \author Janus Juul Eriksen & Patrick Ettenhuber
  !\> \date November 2012
  subroutine manual_312_reordering(bs,dims,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the original array
    integer, intent(in) :: dims(3)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1,pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3))
    !> reordered array
    real(realk),intent(inout) :: array_out(dims(3),dims(1),dims(2))
    integer :: a,b,c,ba,bb,bc,da,db,dc,da2,db2,dc2,bcntr
    logical :: moda,modb,modc
    
    da=dims(1)
    db=dims(2)
    dc=dims(3)

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)

    bcntr=bs-1

    if (pre2 .eq. 0.0E0_realk .and. pre1 .eq. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,ba,bb,bc),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,da2,db2,dc2,bcntr,moda,modb,modc,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
            do bc=1,dc2,bs
     
              do b=0,bcntr
                do a=0,bcntr
                  do c=0,bcntr
                    array_out(bc+c,ba+a,bb+b)=array_in(ba+a,bb+b,bc+c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
     
            do b=0,bcntr
              do a=da2+1,da
                do c=0,bcntr
                  array_out(bc+c,a,bb+b)=array_in(a,bb+b,bc+c)
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
   
            do b=db2+1,db
              do a=0,bcntr
                do c=0,bcntr
                  array_out(bc+c,ba+a,b)=array_in(ba+a,b,bc+c)
                enddo
              enddo
            enddo
   
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc) then
        !$OMP DO
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(c,ba+a,bb+b)=array_in(ba+a,bb+b,c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. dc2 .gt. 0) then
        do bc=1,dc2,bs
     
          do b=db2+1,db
            do a=da2+1,da
              do c=0,bcntr
                array_out(bc+c,a,b)=array_in(a,b,bc+c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. db2 .gt. 0 .and. modc) then
        do bb=1,db2,bs
     
          do c=dc2+1,dc
            do b=0,bcntr
              do a=da2+1,da
                array_out(c,a,bb+b)=array_in(a,bb+b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (da2 .gt. 0 .and. modb .and. modc) then
        do ba=1,da2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do a=0,bcntr
                array_out(c,ba+a,b)=array_in(ba+a,b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. modb .and. modc) then
        do b=db2+1,db
          do a=da2+1,da
            do c=dc2+1,dc
              array_out(c,a,b)=array_in(a,b,c)
            enddo
          enddo
        enddo
      endif
    elseif (pre2 .eq. 0.0E0_realk .and. pre1 .ne. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,ba,bb,bc),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,da2,db2,dc2,bcntr,moda,modb,modc,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
            do bc=1,dc2,bs
     
              do b=0,bcntr
                do a=0,bcntr
                  do c=0,bcntr
                    array_out(bc+c,ba+a,bb+b)=pre1*array_in(ba+a,bb+b,bc+c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
     
            do b=0,bcntr
              do a=da2+1,da
                do c=0,bcntr
                  array_out(bc+c,a,bb+b)=pre1*array_in(a,bb+b,bc+c)
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
   
            do b=db2+1,db
              do a=0,bcntr
                do c=0,bcntr
                  array_out(bc+c,ba+a,b)=pre1*array_in(ba+a,b,bc+c)
                enddo
              enddo
            enddo
   
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc) then
        !$OMP DO
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(c,ba+a,bb+b)=pre1*array_in(ba+a,bb+b,c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. dc2 .gt. 0) then
        do bc=1,dc2,bs
     
          do b=db2+1,db
            do a=da2+1,da
              do c=0,bcntr
                array_out(bc+c,a,b)=pre1*array_in(a,b,bc+c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. db2 .gt. 0 .and. modc) then
        do bb=1,db2,bs
     
          do c=dc2+1,dc
            do b=0,bcntr
              do a=da2+1,da
                array_out(c,a,bb+b)=pre1*array_in(a,bb+b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (da2 .gt. 0 .and. modb .and. modc) then
        do ba=1,da2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do a=0,bcntr
                array_out(c,ba+a,b)=pre1*array_in(ba+a,b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. modb .and. modc) then
        do b=db2+1,db
          do a=da2+1,da
            do c=dc2+1,dc
              array_out(c,a,b)=pre1*array_in(a,b,c)
            enddo
          enddo
        enddo
      endif
    elseif (pre2 .ne. 0.0E0_realk .and. pre1 .eq. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,ba,bb,bc),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,da2,db2,dc2,bcntr,moda,modb,modc,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
            do bc=1,dc2,bs
     
              do b=0,bcntr
                do a=0,bcntr
                  do c=0,bcntr
                    array_out(bc+c,ba+a,bb+b)=pre2*array_out(bc+c,ba+a,bb+b)&
                                             & + array_in(ba+a,bb+b,bc+c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
     
            do b=0,bcntr
              do a=da2+1,da
                do c=0,bcntr
                  array_out(bc+c,a,bb+b)=pre2*array_out(bc+c,a,bb+b)&
                                        & + array_in(a,bb+b,bc+c)
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
   
            do b=db2+1,db
              do a=0,bcntr
                do c=0,bcntr
                  array_out(bc+c,ba+a,b)=pre2*array_out(bc+c,ba+a,b)&
                                        & + array_in(ba+a,b,bc+c)
                enddo
              enddo
            enddo
   
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc) then
        !$OMP DO
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(c,ba+a,bb+b)=pre2*array_out(c,ba+a,bb+b)&
                                          & + array_in(ba+a,bb+b,c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. dc2 .gt. 0) then
        do bc=1,dc2,bs
     
          do b=db2+1,db
            do a=da2+1,da
              do c=0,bcntr
                array_out(bc+c,a,b)=pre2*array_out(bc+c,a,b)&
                                   & + array_in(a,b,bc+c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. db2 .gt. 0 .and. modc) then
        do bb=1,db2,bs
     
          do c=dc2+1,dc
            do b=0,bcntr
              do a=da2+1,da
                array_out(c,a,bb+b)=pre2*array_out(c,a,bb+b)&
                                   & + array_in(a,bb+b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (da2 .gt. 0 .and. modb .and. modc) then
        do ba=1,da2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do a=0,bcntr
                array_out(c,ba+a,b)=pre2*array_out(c,ba+a,b)&
                                   & + array_in(ba+a,b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. modb .and. modc) then
        do b=db2+1,db
          do a=da2+1,da
            do c=dc2+1,dc
              array_out(c,a,b)=pre2*array_out(c,a,b)&
                              & + array_in(a,b,c)
            enddo
          enddo
        enddo
      endif
    elseif (pre2 .ne. 0.0E0_realk .and. pre1 .ne. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,ba,bb,bc),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,da2,db2,dc2,bcntr,moda,modb,modc,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do ba=1,da2,bs
            do bc=1,dc2,bs
     
              do b=0,bcntr
                do a=0,bcntr
                  do c=0,bcntr
                    array_out(bc+c,ba+a,bb+b)=pre2*array_out(bc+c,ba+a,bb+b)&
                                             & + pre1*array_in(ba+a,bb+b,bc+c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
     
            do b=0,bcntr
              do a=da2+1,da
                do c=0,bcntr
                  array_out(bc+c,a,bb+b)=pre2*array_out(bc+c,a,bb+b)&
                                        & + pre1*array_in(a,bb+b,bc+c)
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
   
            do b=db2+1,db
              do a=0,bcntr
                do c=0,bcntr
                  array_out(bc+c,ba+a,b)=pre2*array_out(bc+c,ba+a,b)&
                                        & + pre1*array_in(ba+a,b,bc+c)
                enddo
              enddo
            enddo
   
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc) then
        !$OMP DO
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(c,ba+a,bb+b)=pre2*array_out(c,ba+a,bb+b)&
                                          & + pre1*array_in(ba+a,bb+b,c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. dc2 .gt. 0) then
        do bc=1,dc2,bs
     
          do b=db2+1,db
            do a=da2+1,da
              do c=0,bcntr
                array_out(bc+c,a,b)=pre2*array_out(bc+c,a,b)&
                                   & + pre1*array_in(a,b,bc+c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. db2 .gt. 0 .and. modc) then
        do bb=1,db2,bs
     
          do c=dc2+1,dc
            do b=0,bcntr
              do a=da2+1,da
                array_out(c,a,bb+b)=pre2*array_out(c,a,bb+b)&
                                   & + pre1*array_in(a,bb+b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (da2 .gt. 0 .and. modb .and. modc) then
        do ba=1,da2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do a=0,bcntr
                array_out(c,ba+a,b)=pre2*array_out(c,ba+a,b)&
                                   & + pre1*array_in(ba+a,b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. modb .and. modc) then
        do b=db2+1,db
          do a=da2+1,da
            do c=dc2+1,dc
              array_out(c,a,b)=pre2*array_out(c,a,b)&
                              & + pre1*array_in(a,b,c)
            enddo
          enddo
        enddo
      endif

    else
      call lsquit("MANUAL REORDERING:Case not found",-1)
    endif

   ! print *,"checking the result 312"
   ! if(pre2==0.0E0_realk)then
   !   do b=1,db
   !     do a=1,da
   !       do c=1,dc
   !         if(array_out(c,a,b)/=pre1*array_in(a,b,c))then
   !           print *,"312 reordering not correct",a,b,c,da,db,dc
   !           stop 0
   !         endif
   !       enddo
   !     enddo
   !   enddo
   ! endif
  end subroutine manual_312_reordering

  !\> \brief reorder a 3 diensional array  to get the indices
  !   in the order 1 3 2 , this is a quite expensive reordering
  !   and thus requires additional attention
  !\> \author Janus Juul Eriksen & Patrick Ettenhuber
  !\> \date November 2012
  subroutine manual_132_reordering(bs,dims,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the original array
    integer, intent(in) :: dims(3)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1,pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3))
    !> reordered array
    real(realk),intent(inout) :: array_out(dims(1),dims(3),dims(2))
    integer :: a,b,c,ba,bb,bc,da,db,dc,da2,db2,dc2,bcntr
    logical :: moda,modb,modc
    
    da=dims(1)
    db=dims(2)
    dc=dims(3)

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)

    bcntr=bs-1

    if (pre2 .eq. 0.0E0_realk .and. pre1 .eq. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,ba,bb,bc),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,da2,db2,dc2,bcntr,moda,modb,modc,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do b=0,bcntr
                do c=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,bc+c,bb+b)=array_in(ba+a,bb+b,bc+c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
     
            do b=0,bcntr
              do c=0,bcntr
                do a=da2+1,da
                  array_out(a,bc+c,bb+b)=array_in(a,bb+b,bc+c)
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
   
            do b=db2+1,db
              do c=0,bcntr
                do a=0,bcntr
                  array_out(ba+a,bc+c,b)=array_in(ba+a,b,bc+c)
                enddo
              enddo
            enddo
   
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc) then
        !$OMP DO
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,c,bb+b)=array_in(ba+a,bb+b,c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. dc2 .gt. 0) then
        do bc=1,dc2,bs
     
          do b=db2+1,db
            do c=0,bcntr
              do a=da2+1,da
                array_out(a,bc+c,b)=array_in(a,b,bc+c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. db2 .gt. 0 .and. modc) then
        do bb=1,db2,bs
     
          do c=dc2+1,dc
            do b=0,bcntr
              do a=da2+1,da
                array_out(a,c,bb+b)=array_in(a,bb+b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (da2 .gt. 0 .and. modb .and. modc) then
        do ba=1,da2,bs
     
          do b=db2+1,db
            do c=dc2+1,dc
              do a=0,bcntr
                array_out(ba+a,c,b)=array_in(ba+a,b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. modb .and. modc) then
        do b=db2+1,db
          do c=dc2+1,dc
            do a=da2+1,da
              array_out(a,c,b)=array_in(a,b,c)
            enddo
          enddo
        enddo
      endif
    elseif (pre2 .eq. 0.0E0_realk .and. pre1 .ne. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,ba,bb,bc),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,da2,db2,dc2,bcntr,moda,modb,modc,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do b=0,bcntr
                do c=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,bc+c,bb+b)=pre1*array_in(ba+a,bb+b,bc+c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
     
            do b=0,bcntr
              do c=0,bcntr
                do a=da2+1,da
                  array_out(a,bc+c,bb+b)=pre1*array_in(a,bb+b,bc+c)
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
   
            do b=db2+1,db
              do c=0,bcntr
                do a=0,bcntr
                  array_out(ba+a,bc+c,b)=pre1*array_in(ba+a,b,bc+c)
                enddo
              enddo
            enddo
   
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc) then
        !$OMP DO
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,c,bb+b)=pre1*array_in(ba+a,bb+b,c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. dc2 .gt. 0) then
        do bc=1,dc2,bs
     
          do b=db2+1,db
            do c=0,bcntr
              do a=da2+1,da
                array_out(a,bc+c,b)=pre1*array_in(a,b,bc+c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. db2 .gt. 0 .and. modc) then
        do bb=1,db2,bs
     
          do c=dc2+1,dc
            do b=0,bcntr
              do a=da2+1,da
                array_out(a,c,bb+b)=pre1*array_in(a,bb+b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (da2 .gt. 0 .and. modb .and. modc) then
        do ba=1,da2,bs
     
          do b=db2+1,db
            do c=dc2+1,dc
              do a=0,bcntr
                array_out(ba+a,c,b)=pre1*array_in(ba+a,b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. modb .and. modc) then
        do b=db2+1,db
          do c=dc2+1,dc
            do a=da2+1,da
              array_out(a,c,b)=pre1*array_in(a,b,c)
            enddo
          enddo
        enddo
      endif
    elseif (pre2 .ne. 0.0E0_realk .and. pre1 .eq. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,ba,bb,bc),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,da2,db2,dc2,bcntr,moda,modb,modc,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do b=0,bcntr
                do c=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,bc+c,bb+b)=pre2*array_out(ba+a,bc+c,bb+b)&
                                             & + array_in(ba+a,bb+b,bc+c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
     
            do b=0,bcntr
              do c=0,bcntr
                do a=da2+1,da
                  array_out(a,bc+c,bb+b)=pre2*array_out(a,bc+c,bb+b)&
                                        & + array_in(a,bb+b,bc+c)
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
   
            do b=db2+1,db
              do c=0,bcntr
                do a=0,bcntr
                  array_out(ba+a,bc+c,b)=pre2*array_out(ba+a,bc+c,b)&
                                        & + array_in(ba+a,b,bc+c)
                enddo
              enddo
            enddo
   
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc) then
        !$OMP DO
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,c,bb+b)=pre2*array_out(ba+a,c,bb+b)&
                                          & + array_in(ba+a,bb+b,c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. dc2 .gt. 0) then
        do bc=1,dc2,bs
     
          do b=db2+1,db
            do c=0,bcntr
              do a=da2+1,da
                array_out(a,bc+c,b)=pre2*array_out(a,bc+c,b)&
                                   & + array_in(a,b,bc+c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. db2 .gt. 0 .and. modc) then
        do bb=1,db2,bs
     
          do c=dc2+1,dc
            do b=0,bcntr
              do a=da2+1,da
                array_out(a,c,bb+b)=pre2*array_out(a,c,bb+b)&
                                   & + array_in(a,bb+b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (da2 .gt. 0 .and. modb .and. modc) then
        do ba=1,da2,bs
     
          do b=db2+1,db
            do c=dc2+1,dc
              do a=0,bcntr
                array_out(ba+a,c,b)=pre2*array_out(ba+a,c,b)&
                                   & + array_in(ba+a,b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. modb .and. modc) then
        do b=db2+1,db
          do c=dc2+1,dc
            do a=da2+1,da
              array_out(a,c,b)=pre2*array_out(a,c,b)&
                              & + array_in(a,b,c)
            enddo
          enddo
        enddo
      endif
    elseif (pre2 .ne. 0.0E0_realk .and. pre1 .ne. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,ba,bb,bc),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,da2,db2,dc2,bcntr,moda,modb,modc,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
            do ba=1,da2,bs
     
              do b=0,bcntr
                do c=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,bc+c,bb+b)=pre2*array_out(ba+a,bc+c,bb+b)&
                                             & + pre1*array_in(ba+a,bb+b,bc+c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bb=1,db2,bs
          do bc=1,dc2,bs
     
            do b=0,bcntr
              do c=0,bcntr
                do a=da2+1,da
                  array_out(a,bc+c,bb+b)=pre2*array_out(a,bc+c,bb+b)&
                                        & + pre1*array_in(a,bb+b,bc+c)
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
   
            do b=db2+1,db
              do c=0,bcntr
                do a=0,bcntr
                  array_out(ba+a,bc+c,b)=pre2*array_out(ba+a,bc+c,b)&
                                        & + pre1*array_in(ba+a,b,bc+c)
                enddo
              enddo
            enddo
   
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc) then
        !$OMP DO
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(ba+a,c,bb+b)=pre2*array_out(ba+a,c,bb+b)&
                                          & + pre1*array_in(ba+a,bb+b,c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. dc2 .gt. 0) then
        do bc=1,dc2,bs
     
          do b=db2+1,db
            do c=0,bcntr
              do a=da2+1,da
                array_out(a,bc+c,b)=pre2*array_out(a,bc+c,b)&
                                   & + pre1*array_in(a,b,bc+c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. db2 .gt. 0 .and. modc) then
        do bb=1,db2,bs
     
          do c=dc2+1,dc
            do b=0,bcntr
              do a=da2+1,da
                array_out(a,c,bb+b)=pre2*array_out(a,c,bb+b)&
                                   & + pre1*array_in(a,bb+b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (da2 .gt. 0 .and. modb .and. modc) then
        do ba=1,da2,bs
     
          do b=db2+1,db
            do c=dc2+1,dc
              do a=0,bcntr
                array_out(ba+a,c,b)=pre2*array_out(ba+a,c,b)&
                                   & + pre1*array_in(ba+a,b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. modb .and. modc) then
        do b=db2+1,db
          do c=dc2+1,dc
            do a=da2+1,da
              array_out(a,c,b)=pre2*array_out(a,c,b)&
                              & + pre1*array_in(a,b,c)
            enddo
          enddo
        enddo
      endif

    else
      call lsquit("MANUAL REORDERING:Case not found",-1)
    endif

   ! print *,"checking the result 132"
   ! if(pre2==0.0E0_realk)then
   !   do b=1,db
   !     do c=1,dc
   !       do a=1,da
   !         if(array_out(a,c,b)/=pre1*array_in(a,b,c))then
   !           print *,"132 reordering not correct",a,b,c,da,db,dc
   !           stop 0
   !         endif
   !       enddo
   !     enddo
   !   enddo
   ! endif
  end subroutine manual_132_reordering

  !\> \brief reorder a 3 diensional array  to get the indices
  !   in the order 2 3 1 , this is a quite expensive reordering
  !   and thus requires additional attention
  !\> \author Janus Juul Eriksen & Patrick Ettenhuber
  !\> \date November 2012
  subroutine manual_231_reordering(bs,dims,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the original array
    integer, intent(in) :: dims(3)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1,pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3))
    !> reordered array
    real(realk),intent(inout) :: array_out(dims(2),dims(3),dims(1))
    integer :: a,b,c,ba,bb,bc,da,db,dc,da2,db2,dc2,bcntr
    logical :: moda,modb,modc
    
    da=dims(1)
    db=dims(2)
    dc=dims(3)

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)

    bcntr=bs-1

    if (pre2 .eq. 0.0E0_realk .and. pre1 .eq. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,ba,bb,bc),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,da2,db2,dc2,bcntr,moda,modb,modc,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bb=1,db2,bs
     
              do a=0,bcntr
                do c=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,bc+c,ba+a)=array_in(ba+a,bb+b,bc+c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do a=da2+1,da
              do c=0,bcntr
                do b=0,bcntr
                  array_out(bb+b,bc+c,a)=array_in(a,bb+b,bc+c)
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
   
            do c=0,bcntr
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,bc+c,ba+a)=array_in(ba+a,b,bc+c)
                enddo
              enddo
            enddo
   
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc) then
        !$OMP DO
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(bb+b,c,ba+a)=array_in(ba+a,bb+b,c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. dc2 .gt. 0) then
        do bc=1,dc2,bs
     
          do a=da2+1,da
            do c=0,bcntr
              do b=db2+1,db
                array_out(b,bc+c,a)=array_in(a,b,bc+c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. db2 .gt. 0 .and. modc) then
        do bb=1,db2,bs
     
          do a=da2+1,da
            do c=dc2+1,dc
              do b=0,bcntr
                array_out(bb+b,c,a)=array_in(a,bb+b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (da2 .gt. 0 .and. modb .and. modc) then
        do ba=1,da2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do a=0,bcntr
                array_out(b,c,ba+a)=array_in(ba+a,b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. modb .and. modc) then
        do a=da2+1,da
          do c=dc2+1,dc
            do b=db2+1,db
              array_out(b,c,a)=array_in(a,b,c)
            enddo
          enddo
        enddo
      endif
    elseif (pre2 .eq. 0.0E0_realk .and. pre1 .ne. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,ba,bb,bc),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,da2,db2,dc2,bcntr,moda,modb,modc,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bb=1,db2,bs
     
              do a=0,bcntr
                do c=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,bc+c,ba+a)=pre1*array_in(ba+a,bb+b,bc+c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do a=da2+1,da
              do c=0,bcntr
                do b=0,bcntr
                  array_out(bb+b,bc+c,a)=pre1*array_in(a,bb+b,bc+c)
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
   
            do c=0,bcntr
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,bc+c,ba+a)=pre1*array_in(ba+a,b,bc+c)
                enddo
              enddo
            enddo
   
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc) then
        !$OMP DO
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(bb+b,c,ba+a)=pre1*array_in(ba+a,bb+b,c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. dc2 .gt. 0) then
        do bc=1,dc2,bs
     
          do a=da2+1,da
            do c=0,bcntr
              do b=db2+1,db
                array_out(b,bc+c,a)=pre1*array_in(a,b,bc+c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. db2 .gt. 0 .and. modc) then
        do bb=1,db2,bs
     
          do a=da2+1,da
            do c=dc2+1,dc
              do b=0,bcntr
                array_out(bb+b,c,a)=pre1*array_in(a,bb+b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (da2 .gt. 0 .and. modb .and. modc) then
        do ba=1,da2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do a=0,bcntr
                array_out(b,c,ba+a)=pre1*array_in(ba+a,b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. modb .and. modc) then
        do a=da2+1,da
          do c=dc2+1,dc
            do b=db2+1,db
              array_out(b,c,a)=pre1*array_in(a,b,c)
            enddo
          enddo
        enddo
      endif
    elseif (pre2 .ne. 0.0E0_realk .and. pre1 .eq. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,ba,bb,bc),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,da2,db2,dc2,bcntr,moda,modb,modc,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bb=1,db2,bs
     
              do a=0,bcntr
                do c=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,bc+c,ba+a)=pre2*array_out(bb+b,bc+c,ba+a)&
                                             & + array_in(ba+a,bb+b,bc+c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do a=da2+1,da
              do c=0,bcntr
                do b=0,bcntr
                  array_out(bb+b,bc+c,a)=pre2*array_out(bb+b,bc+c,a)&
                                        & + array_in(a,bb+b,bc+c)
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
   
            do c=0,bcntr
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,bc+c,ba+a)=pre2*array_out(b,bc+c,ba+a)&
                                        & + array_in(ba+a,b,bc+c)
                enddo
              enddo
            enddo
   
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc) then
        !$OMP DO
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(bb+b,c,ba+a)=pre2*array_out(bb+b,c,ba+a)&
                                          & + array_in(ba+a,bb+b,c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. dc2 .gt. 0) then
        do bc=1,dc2,bs
     
          do a=da2+1,da
            do c=0,bcntr
              do b=db2+1,db
                array_out(b,bc+c,a)=pre2*array_out(b,bc+c,a)&
                                   & + array_in(a,b,bc+c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. db2 .gt. 0 .and. modc) then
        do bb=1,db2,bs
     
          do a=da2+1,da
            do c=dc2+1,dc
              do b=0,bcntr
                array_out(bb+b,c,a)=pre2*array_out(bb+b,c,a)&
                                   & + array_in(a,bb+b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (da2 .gt. 0 .and. modb .and. modc) then
        do ba=1,da2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do a=0,bcntr
                array_out(b,c,ba+a)=pre2*array_out(b,c,ba+a)&
                                   & + array_in(ba+a,b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. modb .and. modc) then
        do a=da2+1,da
          do c=dc2+1,dc
            do b=db2+1,db
              array_out(b,c,a)=pre2*array_out(b,c,a)&
                              & + array_in(a,b,c)
            enddo
          enddo
        enddo
      endif
    elseif (pre2 .ne. 0.0E0_realk .and. pre1 .ne. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,ba,bb,bc),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,da2,db2,dc2,bcntr,moda,modb,modc,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do ba=1,da2,bs
          do bc=1,dc2,bs
            do bb=1,db2,bs
     
              do a=0,bcntr
                do c=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,bc+c,ba+a)=pre2*array_out(bb+b,bc+c,ba+a)&
                                             & + pre1*array_in(ba+a,bb+b,bc+c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do a=da2+1,da
              do c=0,bcntr
                do b=0,bcntr
                  array_out(bb+b,bc+c,a)=pre2*array_out(bb+b,bc+c,a)&
                                        & + pre1*array_in(a,bb+b,bc+c)
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
   
            do c=0,bcntr
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,bc+c,ba+a)=pre2*array_out(b,bc+c,ba+a)&
                                        & + pre1*array_in(ba+a,b,bc+c)
                enddo
              enddo
            enddo
   
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc) then
        !$OMP DO
          do bb=1,db2,bs
            do ba=1,da2,bs
     
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=0,bcntr
                    array_out(bb+b,c,ba+a)=pre2*array_out(bb+b,c,ba+a)&
                                          & + pre1*array_in(ba+a,bb+b,c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. dc2 .gt. 0) then
        do bc=1,dc2,bs
     
          do a=da2+1,da
            do c=0,bcntr
              do b=db2+1,db
                array_out(b,bc+c,a)=pre2*array_out(b,bc+c,a)&
                                   & + pre1*array_in(a,b,bc+c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. db2 .gt. 0 .and. modc) then
        do bb=1,db2,bs
     
          do a=da2+1,da
            do c=dc2+1,dc
              do b=0,bcntr
                array_out(bb+b,c,a)=pre2*array_out(bb+b,c,a)&
                                   & + pre1*array_in(a,bb+b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (da2 .gt. 0 .and. modb .and. modc) then
        do ba=1,da2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do a=0,bcntr
                array_out(b,c,ba+a)=pre2*array_out(b,c,ba+a)&
                                   & + pre1*array_in(ba+a,b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. modb .and. modc) then
        do a=da2+1,da
          do c=dc2+1,dc
            do b=db2+1,db
              array_out(b,c,a)=pre2*array_out(b,c,a)&
                              & + pre1*array_in(a,b,c)
            enddo
          enddo
        enddo
      endif

    else
      call lsquit("MANUAL REORDERING:Case not found",-1)
    endif

   ! print *,"checking the result 231"
   ! if(pre2==0.0E0_realk)then
   !   do a=1,da
   !     do c=1,dc
   !       do b=1,db
   !         if(array_out(b,c,a)/=pre1*array_in(a,b,c))then
   !           print *,"231 reordering not correct",a,b,c,da,db,dc
   !           stop 0
   !         endif
   !       enddo
   !     enddo
   !   enddo
   ! endif
  end subroutine manual_231_reordering

  !\> \brief reorder a 3 diensional array  to get the indices
  !   in the order 2 1 3 , this is a quite expensive reordering
  !   and thus requires additional attention
  !\> \author Janus Juul Eriksen & Patrick Ettenhuber
  !\> \date November 2012
  subroutine manual_213_reordering(bs,dims,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the original array
    integer, intent(in) :: dims(3)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1,pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3))
    !> reordered array
    real(realk),intent(inout) :: array_out(dims(2),dims(1),dims(3))
    integer :: a,b,c,ba,bb,bc,da,db,dc,da2,db2,dc2,bcntr
    logical :: moda,modb,modc
    
    da=dims(1)
    db=dims(2)
    dc=dims(3)

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)

    bcntr=bs-1

    if (pre2 .eq. 0.0E0_realk .and. pre1 .eq. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,ba,bb,bc),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,da2,db2,dc2,bcntr,moda,modb,modc,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
            do bb=1,db2,bs
     
              do c=0,bcntr
                do a=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,ba+a,bc+c)=array_in(ba+a,bb+b,bc+c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do c=0,bcntr
              do a=da2+1,da
                do b=0,bcntr
                  array_out(bb+b,a,bc+c)=array_in(a,bb+b,bc+c)
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
   
            do c=0,bcntr
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,ba+a,bc+c)=array_in(ba+a,b,bc+c)
                enddo
              enddo
            enddo
   
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc) then
        !$OMP DO
          do ba=1,da2,bs
            do bb=1,db2,bs
     
              do c=dc2+1,dc
                do a=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,ba+a,c)=array_in(ba+a,bb+b,c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. dc2 .gt. 0) then
        do bc=1,dc2,bs
     
          do c=0,bcntr
            do a=da2+1,da
              do b=db2+1,db
                array_out(b,a,bc+c)=array_in(a,b,bc+c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. db2 .gt. 0 .and. modc) then
        do bb=1,db2,bs
     
          do c=dc2+1,dc
            do a=da2+1,da
              do b=0,bcntr
                array_out(bb+b,a,c)=array_in(a,bb+b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (da2 .gt. 0 .and. modb .and. modc) then
        do ba=1,da2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do a=0,bcntr
                array_out(b,ba+a,c)=array_in(ba+a,b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. modb .and. modc) then
        do c=dc2+1,dc
          do a=da2+1,da
            do b=db2+1,db
              array_out(b,a,c)=array_in(a,b,c)
            enddo
          enddo
        enddo
      endif
    elseif (pre2 .eq. 0.0E0_realk .and. pre1 .ne. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,ba,bb,bc),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,da2,db2,dc2,bcntr,moda,modb,modc,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
            do bb=1,db2,bs
     
              do c=0,bcntr
                do a=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,ba+a,bc+c)=pre1*array_in(ba+a,bb+b,bc+c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do c=0,bcntr
              do a=da2+1,da
                do b=0,bcntr
                  array_out(bb+b,a,bc+c)=pre1*array_in(a,bb+b,bc+c)
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
   
            do c=0,bcntr
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,ba+a,bc+c)=pre1*array_in(ba+a,b,bc+c)
                enddo
              enddo
            enddo
   
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc) then
        !$OMP DO
          do ba=1,da2,bs
            do bb=1,db2,bs
     
              do c=dc2+1,dc
                do a=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,ba+a,c)=pre1*array_in(ba+a,bb+b,c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. dc2 .gt. 0) then
        do bc=1,dc2,bs
     
          do c=0,bcntr
            do a=da2+1,da
              do b=db2+1,db
                array_out(b,a,bc+c)=pre1*array_in(a,b,bc+c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. db2 .gt. 0 .and. modc) then
        do bb=1,db2,bs
     
          do c=dc2+1,dc
            do a=da2+1,da
              do b=0,bcntr
                array_out(bb+b,a,c)=pre1*array_in(a,bb+b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (da2 .gt. 0 .and. modb .and. modc) then
        do ba=1,da2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do a=0,bcntr
                array_out(b,ba+a,c)=pre1*array_in(ba+a,b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. modb .and. modc) then
        do c=dc2+1,dc
          do a=da2+1,da
            do b=db2+1,db
              array_out(b,a,c)=pre1*array_in(a,b,c)
            enddo
          enddo
        enddo
      endif
    elseif (pre2 .ne. 0.0E0_realk .and. pre1 .eq. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,ba,bb,bc),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,da2,db2,dc2,bcntr,moda,modb,modc,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
            do bb=1,db2,bs
     
              do c=0,bcntr
                do a=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,ba+a,bc+c)=pre2*array_out(bb+b,ba+a,bc+c)&
                                             & + array_in(ba+a,bb+b,bc+c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do c=0,bcntr
              do a=da2+1,da
                do b=0,bcntr
                  array_out(bb+b,a,bc+c)=pre2*array_out(bb+b,a,bc+c)&
                                        & + array_in(a,bb+b,bc+c)
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
   
            do c=0,bcntr
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,ba+a,bc+c)=pre2*array_out(b,ba+a,bc+c)&
                                        & + array_in(ba+a,b,bc+c)
                enddo
              enddo
            enddo
   
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc) then
        !$OMP DO
          do ba=1,da2,bs
            do bb=1,db2,bs
     
              do c=dc2+1,dc
                do a=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,ba+a,c)=pre2*array_out(bb+b,ba+a,c)&
                                          & + array_in(ba+a,bb+b,c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. dc2 .gt. 0) then
        do bc=1,dc2,bs
     
          do c=0,bcntr
            do a=da2+1,da
              do b=db2+1,db
                array_out(b,a,bc+c)=pre2*array_out(b,a,bc+c)&
                                   & + array_in(a,b,bc+c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. db2 .gt. 0 .and. modc) then
        do bb=1,db2,bs
     
          do c=dc2+1,dc
            do a=da2+1,da
              do b=0,bcntr
                array_out(bb+b,a,c)=pre2*array_out(bb+b,a,c)&
                                   & + array_in(a,bb+b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (da2 .gt. 0 .and. modb .and. modc) then
        do ba=1,da2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do a=0,bcntr
                array_out(b,ba+a,c)=pre2*array_out(b,ba+a,c)&
                                   & + array_in(ba+a,b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. modb .and. modc) then
        do c=dc2+1,dc
          do a=da2+1,da
            do b=db2+1,db
              array_out(b,a,c)=pre2*array_out(b,a,c)&
                              & + array_in(a,b,c)
            enddo
          enddo
        enddo
      endif
    elseif (pre2 .ne. 0.0E0_realk .and. pre1 .ne. 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,ba,bb,bc),SHARED(array_in,array_out,&
      !$OMP& da,db,dc,da2,db2,dc2,bcntr,moda,modb,modc,pre1,pre2,bs)
      if ( da2 .gt. 0 .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
            do bb=1,db2,bs
     
              do c=0,bcntr
                do a=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,ba+a,bc+c)=pre2*array_out(bb+b,ba+a,bc+c)&
                                             & + pre1*array_in(ba+a,bb+b,bc+c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (moda .and. db2 .gt. 0 .and. dc2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do bb=1,db2,bs
     
            do c=0,bcntr
              do a=da2+1,da
                do b=0,bcntr
                  array_out(bb+b,a,bc+c)=pre2*array_out(bb+b,a,bc+c)&
                                        & + pre1*array_in(a,bb+b,bc+c)
                enddo
              enddo
            enddo
     
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. modb .and. dc2 .gt. 0) then
        !$OMP DO
        do bc=1,dc2,bs
          do ba=1,da2,bs
   
            do c=0,bcntr
              do b=db2+1,db
                do a=0,bcntr
                  array_out(b,ba+a,bc+c)=pre2*array_out(b,ba+a,bc+c)&
                                        & + pre1*array_in(ba+a,b,bc+c)
                enddo
              enddo
            enddo
   
          enddo
        enddo
        !$OMP END DO NOWAIT
      endif
      if (da2 .gt. 0 .and. db2 .gt. 0 .and. modc) then
        !$OMP DO
          do ba=1,da2,bs
            do bb=1,db2,bs
     
              do c=dc2+1,dc
                do a=0,bcntr
                  do b=0,bcntr
                    array_out(bb+b,ba+a,c)=pre2*array_out(bb+b,ba+a,c)&
                                          & + pre1*array_in(ba+a,bb+b,c)
                  enddo
                enddo
              enddo
     
            enddo
          enddo
        !$OMP END DO NOWAIT
      endif
      !$OMP END PARALLEL
      if (moda .and. modb .and. dc2 .gt. 0) then
        do bc=1,dc2,bs
     
          do c=0,bcntr
            do a=da2+1,da
              do b=db2+1,db
                array_out(b,a,bc+c)=pre2*array_out(b,a,bc+c)&
                                   & + pre1*array_in(a,b,bc+c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. db2 .gt. 0 .and. modc) then
        do bb=1,db2,bs
     
          do c=dc2+1,dc
            do a=da2+1,da
              do b=0,bcntr
                array_out(bb+b,a,c)=pre2*array_out(bb+b,a,c)&
                                   & + pre1*array_in(a,bb+b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (da2 .gt. 0 .and. modb .and. modc) then
        do ba=1,da2,bs
     
          do c=dc2+1,dc
            do b=db2+1,db
              do a=0,bcntr
                array_out(b,ba+a,c)=pre2*array_out(b,ba+a,c)&
                                   & + pre1*array_in(ba+a,b,c)
              enddo
            enddo
          enddo
     
        enddo
      endif
      if (moda .and. modb .and. modc) then
        do c=dc2+1,dc
          do a=da2+1,da
            do b=db2+1,db
              array_out(b,a,c)=pre2*array_out(b,a,c)&
                              & + pre1*array_in(a,b,c)
            enddo
          enddo
        enddo
      endif

    else
      call lsquit("MANUAL REORDERING:Case not found",-1)
    endif

   ! print *,"checking the result 213"
   ! if(pre2==0.0E0_realk)then
   !   do c=1,dc
   !     do a=1,da
   !       do b=1,db
   !         if(array_out(b,a,c)/=pre1*array_in(a,b,c))then
   !           print *,"213 reordering not correct",a,b,c,da,db,dc
   !           stop 0
   !         endif
   !       enddo
   !     enddo
   !   enddo
   ! endif
  end subroutine manual_213_reordering
end module manual_reorderings_module
