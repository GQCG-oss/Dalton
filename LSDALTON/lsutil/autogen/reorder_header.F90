
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
    if(order(1)==3 .and. order(2)==2 .and. &
         order(3)==4 .and. order(4)==1) order_type = 21
    if(order(1)==4 .and. order(2)==3 .and. &
         order(3)==1 .and. order(4)==2) order_type = 22
    if(order(1)==4 .and. order(2)==2 .and. &
         order(3)==1 .and. order(4)==3) order_type = 23

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
    case(21)
       call manual_3241_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case(22)
       call manual_4312_reordering(block_size,dims,pre1,array_in,pre2,array_out)
    case(23)
       call manual_4213_reordering(block_size,dims,pre1,array_in,pre2,array_out)
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
