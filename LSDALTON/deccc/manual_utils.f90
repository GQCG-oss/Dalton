!\> \brief additional stuff to manual reorderings, generalizations of
!some routines 
!> \author Patrick Ettenhuber
!\> \date March 2013
module manual_utils_module
  use precision
  use dec_typedef_module
  contains
  
  !\> \brief reorder a 4 diensional tile into the full matrix  to get the indices
  !   in the order 2 1 4 3 
  !\> \autor Patrick Ettenhuber
  !\> \date March 2013
  subroutine manual_2143_reordering_full2tile(bs,dims,fdims,fels,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the tile
    integer, intent(in) :: dims(4)
    !>  the dimensions of the different modes in the full array
    integer, intent(in) :: fdims(4)
    !> first elements in the tile corresponding to the full array 
    integer, intent(in) :: fels(4)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1
    real(realk),intent(in) :: pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3),dims(4))
    !> reordered array
    real(realk),intent(inout) :: array_out(fdims(2),fdims(1),fdims(4),fdims(3))
    integer :: a,b,c,d,ba,bb,bc,bd,da,db,dc,dd,da2,db2,dc2,dd2,bcntr
    integer :: baf,bbf,bcf,bdf,fa,fb,fc,fd
    logical :: moda,modb,modc,modd

    da=dims(1)
    db=dims(2)
    dc=dims(3)
    dd=dims(4)

    fa=fels(1)-1
    fb=fels(2)-1
    fc=fels(3)-1
    fd=fels(4)-1

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs
    dd2=(dd/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)
    modd=(mod(dd,bs)>0)

    bcntr=bs-1

    !elseif (pre2 /= 0.0E0_realk .and. pre1 /= 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd,baf,bbf,bcf,bdf),SHARED(array_in,array_out,&
      !$OMP& !da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs,fa,fb,fc,fd)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          bcf = fc + bc
          do bd=1,dd2,bs
            bdf = fd + bd
            do ba=1,da2,bs
              baf = fa + ba
              do bb=1,db2,bs
                bbf = fb + bb
     
                do c=0,bcntr
                  do d=0,bcntr
                    do a=0,bcntr
                      do b=0,bcntr
                        array_out(bb+b,ba+a,bd+d,bc+c)=pre2*array_out(bb+b,ba+a,bd+d,bc+c)&
                                                     &+pre1*array_in(baf+a,bbf+b,bcf+c,bdf+d)
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
          bcf = fc + bc
          do bd=1,dd2,bs
            bdf = fd + bd
            do bb=1,db2,bs
              bbf = fb + bb
     
              do c=0,bcntr
                do d=0,bcntr
                  do a=da2+1,da
                    do b=0,bcntr
                      array_out(bb+b,a,bd+d,bc+c)=pre2*array_out(bb+b,a,bd+d,bc+c)&
                                                   &+pre1*array_in(fa+a,bbf+b,bcf+c,bdf+d)
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
          bdf = fd + bd
          do bc=1,dc2,bs
            bcf = fc + bc
            do ba=1,da2,bs
              baf = fa + ba
     
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do a=0,bcntr
                      array_out(fb+b,baf+a,bdf+d,bcf+c)=pre2*array_out(b,ba+a,bd+d,bc+c)&
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
          bdf = fd + bd
          do ba=1,da2,bs
            baf = fa + ba
            do bb=1,db2,bs
              bbf = fb + bb
     
              do c=dc2+1,dc
                do d=0,bcntr
                  do a=0,bcntr
                    do b=0,bcntr
                      array_out(bbf+b,baf+a,bdf+d,fc+c)=pre2*array_out(bbf+b,baf+a,bdf+d,fc+c)&
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
          bcf = fc + bc
          do bb=1,db2,bs
            bbf = fb + bb
            do ba=1,da2,bs
              baf = fa + ba
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bbf+b,baf+a,fd+d,bcf+c)=pre2*array_out(bbf+b,baf+a,fd+d,bcf+c)&
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
          bdf = fd + bd
          do bc=1,dc2,bs
            bcf = fc + bc
     
            do d=0,bcntr
              do c=0,bcntr
                do b=db2+1,db
                  do a=da2+1,da
                    array_out(fb+b,fa+a,bdf+d,bcf+c)=pre2*array_out(fb+b,fa+a,bdf+d,bcf+c)&
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
          bdf = fd + bd
          do bb=1,db2,bs
            bbf = fb + bb
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(bbf+b,fa+a,bdf+d,fc+c)=pre2*array_out(bbf+b,fa+a,bdf+d,fc+c)&
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
          bcf = fc + bc
          do bb=1,db2,bs
            bbf = fb + bb
     
            do c=0,bcntr
              do d=dd2+1,dd
                do a=da2+1,da
                  do b=0,bcntr
                    array_out(bbf+b,fa+a,fd+d,bcf+c)=pre2*array_out(bbf+b,fa+a,fd+d,bcf+c)&
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
          bdf = fd + bd
          do ba=1,da2,bs
            baf = fa + ba
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(fb+b,baf+a,bdf+d,fc+c)=pre2*array_out(fb+b,baf+a,bdf+d,fc+c)&
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
          bcf = fc + bc
          do ba=1,da2,bs
            baf = fa + ba
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(fb+b,baf+a,fd+d,bcf+c)=pre2*array_out(fb+b,baf+a,fd+d,bcf+c)&
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
          baf = fa + ba
          do bb=1,db2,bs
            bbf = fb + bb
     
            do c=dc2+1,dc
              do d=dd2+1,dd
                do a=0,bcntr
                  do b=0,bcntr
                    array_out(bbf+b,baf+a,fd+d,fc+c)=pre2*array_out(bbf+b,baf+a,fd+d,fc+c)&
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
          baf = fa + ba
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(fb+b,baf+a,fd+d,fc+c)=pre2*array_out(fb+b,baf+a,fd+d,fc+c)&
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
          bbf = fb + bb
     
          do c=dc2+1,dc
            do d=dd2+1,dd
              do a=da2+1,da
                do b=0,bcntr
                  array_out(bbf+b,fa+a,fd+d,fc+c)=pre2*array_out(bbf+b,fa+a,fd+d,fc+c)&
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
          bcf = fc + bc
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(fb+b,fa+a,fd+d,bcf+c)=pre2*array_out(fb+b,fa+a,fd+d,bcf+c)&
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
          bdf = fd + bd
     
          do c=dc2+1,dc
            do d=0,bcntr
              do a=da2+1,da
                do b=db2+1,db
                  array_out(fb+b,fa+a,bdf+d,fc+c)=pre2*array_out(fb+b,fa+a,bdf+d,fc+c)&
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
                array_out(fb+b,fa+a,fd+d,fc+c)=pre2*array_out(fb+b,fa+a,fd+d,fc+c)&
                                             &+pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    !else
    !  call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
    !endif

    print *,"check 2143"
    if(pre2==0.0E0_realk)then
      do b=1,db
        do d=1,dd
          do a=1,da
            do c=1,dc
              if(array_out(fb+b,fa+a,fd+d,fc+c)/=pre1*array_in(a,b,c,d))then
                print *,"1432 reordering not correct",a,b,c,d,da,db,dc,dd
                print *,array_out(fb+b,fa+a,fd+d,fc+c),array_in(a,b,c,d)
                stop 0
              endif
            enddo
          enddo
        enddo
      enddo
    endif
  end subroutine manual_2143_reordering_full2tile

  !\> \brief reorder a 4 diensional tile into the full matrix  to get the indices
  !   in the order 2 1 4 3 
  !\> \autor Patrick Ettenhuber
  !\> \date March 2013
  subroutine manual_2143_reordering_tile2full(bs,dims,fdims,fels,pre1,array_in,pre2,array_out)
    implicit none
    !> input for the block size in tiled reordering
    integer, intent(in) :: bs
    !>  the dimensions of the different modes in the tile
    integer, intent(in) :: dims(4)
    !>  the dimensions of the different modes in the full array
    integer, intent(in) :: fdims(4)
    !> first elements in the tile corresponding to the full array 
    integer, intent(in) :: fels(4)
    !> as this routine can be used for adding and scaling these are the prefactors
    real(realk),intent(in) :: pre1
    real(realk),intent(in) :: pre2
    !> array to be reordered
    real(realk),intent(in) :: array_in(dims(1),dims(2),dims(3),dims(4))
    !> reordered array
    real(realk),intent(inout) :: array_out(fdims(2),fdims(1),fdims(4),fdims(3))
    integer :: a,b,c,d,ba,bb,bc,bd,da,db,dc,dd,da2,db2,dc2,dd2,bcntr
    integer :: baf,bbf,bcf,bdf,fa,fb,fc,fd
    logical :: moda,modb,modc,modd

    da=dims(1)
    db=dims(2)
    dc=dims(3)
    dd=dims(4)

    fa=fels(1)-1
    fb=fels(2)-1
    fc=fels(3)-1
    fd=fels(4)-1

    da2=(da/bs)*bs
    db2=(db/bs)*bs
    dc2=(dc/bs)*bs
    dd2=(dd/bs)*bs

    moda=(mod(da,bs)>0)
    modb=(mod(db,bs)>0)
    modc=(mod(dc,bs)>0)
    modd=(mod(dd,bs)>0)

    bcntr=bs-1

    if (pre2 == 0.0E0_realk .and. pre1 == 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd,baf,bbf,bcf,bdf),SHARED(array_in,array_out,&
      !$OMP& !da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs,fa,fb,fc,fd)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          bcf = fc + bc
          do bd=1,dd2,bs
            bdf = fd + bd
            do ba=1,da2,bs
              baf = fa + ba
              do bb=1,db2,bs
                bbf = fb + bb
     
                do c=0,bcntr
                  do d=0,bcntr
                    do a=0,bcntr
                      do b=0,bcntr
                        array_out(bbf+b,baf+a,bdf+d,bcf+c)=array_in(ba+a,bb+b,bc+c,bd+d)
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
          bcf = fc + bc
          do bd=1,dd2,bs
            bdf = fd + bd
            do bb=1,db2,bs
              bbf = fb + bb
     
              do c=0,bcntr
                do d=0,bcntr
                  do a=da2+1,da
                    do b=0,bcntr
                      array_out(bbf+b,fa+a,bdf+d,bcf+c)=array_in(a,bb+b,bc+c,bd+d)
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
          bdf = fd + bd
          do bc=1,dc2,bs
            bcf = fc + bc
            do ba=1,da2,bs
              baf = fa + ba
     
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do a=0,bcntr
                      array_out(fb+b,baf+a,bdf+d,bcf+c)=array_in(ba+a,b,bc+c,bd+d)
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
          bdf = fd + bd
          do ba=1,da2,bs
            baf = fa + ba
            do bb=1,db2,bs
              bbf = fb + bb
     
              do c=dc2+1,dc
                do d=0,bcntr
                  do a=0,bcntr
                    do b=0,bcntr
                      array_out(bbf+b,baf+a,bdf+d,fc+c)=array_in(ba+a,bb+b,c,bd+d)
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
          bcf = fc + bc
          do bb=1,db2,bs
            bbf = fb + bb
            do ba=1,da2,bs
              baf = fa + ba
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bbf+b,baf+a,fd+d,bcf+c)=array_in(ba+a,bb+b,bc+c,d)
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
          bdf = fd + bd
          do bc=1,dc2,bs
            bcf = fc + bc
     
            do d=0,bcntr
              do c=0,bcntr
                do b=db2+1,db
                  do a=da2+1,da
                    array_out(fb+b,fa+a,bdf+d,bcf+c)=array_in(a,b,bc+c,bd+d)
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
          bdf = fd + bd
          do bb=1,db2,bs
            bbf = fb + bb
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(bbf+b,fa+a,bdf+d,fc+c)=array_in(a,bb+b,c,bd+d)
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
          bcf = fc + bc
          do bb=1,db2,bs
            bbf = fb + bb
     
            do c=0,bcntr
              do d=dd2+1,dd
                do a=da2+1,da
                  do b=0,bcntr
                    array_out(bbf+b,fa+a,fd+d,bcf+c)=array_in(a,bb+b,bc+c,d)
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
          bdf = fd + bd
          do ba=1,da2,bs
            baf = fa + ba
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(fb+b,baf+a,bdf+d,fc+c)=array_in(ba+a,b,c,bd+d)
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
          bcf = fc + bc
          do ba=1,da2,bs
            baf = fa + ba
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(fb+b,baf+a,fd+d,bcf+c)=array_in(ba+a,b,bc+c,d)
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
          baf = fa + ba
          do bb=1,db2,bs
            bbf = fb + bb
     
            do c=dc2+1,dc
              do d=dd2+1,dd
                do a=0,bcntr
                  do b=0,bcntr
                    array_out(bbf+b,baf+a,fd+d,fc+c)=array_in(ba+a,bb+b,c,d)
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
          baf = fa + ba
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(fb+b,baf+a,fd+d,fc+c)=array_in(ba+a,b,c,d)
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
          bbf = fb + bb
     
          do c=dc2+1,dc
            do d=dd2+1,dd
              do a=da2+1,da
                do b=0,bcntr
                  array_out(bbf+b,fa+a,fd+d,fc+c)=array_in(a,bb+b,c,d)
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
          bcf = fc + bc
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(fb+b,fa+a,fd+d,bcf+c)=array_in(a,b,bc+c,d)
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
          bdf = fd + bd
     
          do c=dc2+1,dc
            do d=0,bcntr
              do a=da2+1,da
                do b=db2+1,db
                  array_out(fb+b,fa+a,bdf+d,fc+c)=array_in(a,b,c,bd+d)
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
                array_out(fb+b,fa+a,fd+d,fc+c)=array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif (pre2 == 0.0E0_realk .and. pre1 /= 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd,baf,bbf,bcf,bdf),SHARED(array_in,array_out,&
      !$OMP& !da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs,fa,fb,fc,fd)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          bcf = fc + bc
          do bd=1,dd2,bs
            bdf = fd + bd
            do ba=1,da2,bs
              baf = fa + ba
              do bb=1,db2,bs
                bbf = fb + bb
     
                do c=0,bcntr
                  do d=0,bcntr
                    do a=0,bcntr
                      do b=0,bcntr
                        array_out(bbf+b,baf+a,bdf+d,bcf+c)=pre1*array_in(ba+a,bb+b,bc+c,bd+d)
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
          bcf = fc + bc
          do bd=1,dd2,bs
            bdf = fd + bd
            do bb=1,db2,bs
              bbf = fb + bb
     
              do c=0,bcntr
                do d=0,bcntr
                  do a=da2+1,da
                    do b=0,bcntr
                      array_out(bbf+b,fa+a,bdf+d,bcf+c)=pre1*array_in(a,bb+b,bc+c,bd+d)
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
          bdf = fd + bd
          do bc=1,dc2,bs
            bcf = fc + bc
            do ba=1,da2,bs
              baf = fa + ba
     
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do a=0,bcntr
                      array_out(fb+b,baf+a,bdf+d,bcf+c)=pre1*array_in(ba+a,b,bc+c,bd+d)
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
          bdf = fd + bd
          do ba=1,da2,bs
            baf = fa + ba
            do bb=1,db2,bs
              bbf = fb + bb
     
              do c=dc2+1,dc
                do d=0,bcntr
                  do a=0,bcntr
                    do b=0,bcntr
                      array_out(bbf+b,baf+a,bdf+d,fc+c)=pre1*array_in(ba+a,bb+b,c,bd+d)
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
          bcf = fc + bc
          do bb=1,db2,bs
            bbf = fb + bb
            do ba=1,da2,bs
              baf = fa + ba
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bbf+b,baf+a,fd+d,bcf+c)=pre1*array_in(ba+a,bb+b,bc+c,d)
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
          bdf = fd + bd
          do bc=1,dc2,bs
            bcf = fc + bc
     
            do d=0,bcntr
              do c=0,bcntr
                do b=db2+1,db
                  do a=da2+1,da
                    array_out(fb+b,fa+a,bdf+d,bcf+c)=pre1*array_in(a,b,bc+c,bd+d)
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
          bdf = fd + bd
          do bb=1,db2,bs
            bbf = fb + bb
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(bbf+b,fa+a,bdf+d,fc+c)=pre1*array_in(a,bb+b,c,bd+d)
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
          bcf = fc + bc
          do bb=1,db2,bs
            bbf = fb + bb
     
            do c=0,bcntr
              do d=dd2+1,dd
                do a=da2+1,da
                  do b=0,bcntr
                    array_out(bbf+b,fa+a,fd+d,bcf+c)=pre1*array_in(a,bb+b,bc+c,d)
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
          bdf = fd + bd
          do ba=1,da2,bs
            baf = fa + ba
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(fb+b,baf+a,bdf+d,fc+c)=pre1*array_in(ba+a,b,c,bd+d)
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
          bcf = fc + bc
          do ba=1,da2,bs
            baf = fa + ba
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(fb+b,baf+a,fd+d,bcf+c)=pre1*array_in(ba+a,b,bc+c,d)
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
          baf = fa + ba
          do bb=1,db2,bs
            bbf = fb + bb
     
            do c=dc2+1,dc
              do d=dd2+1,dd
                do a=0,bcntr
                  do b=0,bcntr
                    array_out(bbf+b,baf+a,fd+d,fc+c)=pre1*array_in(ba+a,bb+b,c,d)
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
          baf = fa + ba
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(fb+b,baf+a,fd+d,fc+c)=pre1*array_in(ba+a,b,c,d)
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
          bbf = fb + bb
     
          do c=dc2+1,dc
            do d=dd2+1,dd
              do a=da2+1,da
                do b=0,bcntr
                  array_out(bbf+b,fa+a,fd+d,fc+c)=pre1*array_in(a,bb+b,c,d)
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
          bcf = fc + bc
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(fb+b,fa+a,fd+d,bcf+c)=pre1*array_in(a,b,bc+c,d)
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
          bdf = fd + bd
     
          do c=dc2+1,dc
            do d=0,bcntr
              do a=da2+1,da
                do b=db2+1,db
                  array_out(fb+b,fa+a,bdf+d,fc+c)=pre1*array_in(a,b,c,bd+d)
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
                array_out(fb+b,fa+a,fd+d,fc+c)=pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif (pre2 /= 0.0E0_realk .and. pre1 == 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd,baf,bbf,bcf,bdf),SHARED(array_in,array_out,&
      !$OMP& !da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs,fa,fb,fc,fd)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          bcf = fc + bc
          do bd=1,dd2,bs
            bdf = fd + bd
            do ba=1,da2,bs
              baf = fa + ba
              do bb=1,db2,bs
                bbf = fb + bb
     
                do c=0,bcntr
                  do d=0,bcntr
                    do a=0,bcntr
                      do b=0,bcntr
                        array_out(bbf+b,baf+a,bdf+d,bcf+c)=pre2*array_out(bbf+b,baf+a,bdf+d,bcf+c)&
                                                     &+ array_in(ba+a,bb+b,bc+c,bd+d)
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
          bcf = fc + bc
          do bd=1,dd2,bs
            bdf = fd + bd
            do bb=1,db2,bs
              bbf = fb + bb
     
              do c=0,bcntr
                do d=0,bcntr
                  do a=da2+1,da
                    do b=0,bcntr
                      array_out(bbf+b,fa+a,bdf+d,bcf+c)=pre2*array_out(bbf+b,fa+a,bdf+d,bcf+c)&
                                                   &+ array_in(a,bb+b,bc+c,bd+d)
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
          bdf = fd + bd
          do bc=1,dc2,bs
            bcf = fc + bc
            do ba=1,da2,bs
              baf = fa + ba
     
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do a=0,bcntr
                      array_out(fb+b,baf+a,bdf+d,bcf+c)=pre2*array_out(fb+b,baf+a,bdf+d,bcf+c)&
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
          bdf = fd + bd
          do ba=1,da2,bs
            baf = fa + ba
            do bb=1,db2,bs
              bbf = fb + bb
     
              do c=dc2+1,dc
                do d=0,bcntr
                  do a=0,bcntr
                    do b=0,bcntr
                      array_out(bbf+b,baf+a,bdf+d,fc+c)=pre2*array_out(bbf+b,baf+a,bdf+d,fc+c)&
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
          bcf = fc + bc
          do bb=1,db2,bs
            bbf = fb + bb
            do ba=1,da2,bs
              baf = fa + ba
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bbf+b,baf+a,fd+d,bcf+c)=pre2*array_out(bbf+b,baf+a,fd+d,bcf+c)&
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
          bdf = fd + bd
          do bc=1,dc2,bs
            bcf = fc + bc
     
            do d=0,bcntr
              do c=0,bcntr
                do b=db2+1,db
                  do a=da2+1,da
                    array_out(fb+b,fa+a,bdf+d,bcf+c)=pre2*array_out(fb+b,fa+a,bdf+d,bcf+c)&
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
          bdf = fd + bd
          do bb=1,db2,bs
            bbf = fb + bb
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(bbf+b,fa+a,bdf+d,fc+c)=pre2*array_out(bbf+b,fa+a,bdf+d,fc+c)&
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
          bcf = fc + bc
          do bb=1,db2,bs
            bbf = fb + bb
     
            do c=0,bcntr
              do d=dd2+1,dd
                do a=da2+1,da
                  do b=0,bcntr
                    array_out(bbf+b,fa+a,fd+d,bcf+c)=pre2*array_out(bbf+b,fa+a,fd+d,bcf+c)&
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
          bdf = fd + bd
          do ba=1,da2,bs
            baf = fa + ba
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(fb+b,baf+a,bdf+d,fc+c)=pre2*array_out(fb+b,baf+a,bdf+d,fc+c)&
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
          bcf = fc + bc
          do ba=1,da2,bs
            baf = fa + ba
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(fb+b,baf+a,fd+d,bcf+c)=pre2*array_out(fb+b,baf+a,fd+d,bcf+c)&
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
          baf = fa + ba
          do bb=1,db2,bs
            bbf = fb + bb
     
            do c=dc2+1,dc
              do d=dd2+1,dd
                do a=0,bcntr
                  do b=0,bcntr
                    array_out(bbf+b,baf+a,fd+d,fc+c)=pre2*array_out(bbf+b,baf+a,fd+d,fc+c)&
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
          baf = fa + ba
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(fb+b,baf+a,fd+d,fc+c)=pre2*array_out(fb+b,baf+a,fd+d,fc+c)&
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
          bbf = fb + bb
     
          do c=dc2+1,dc
            do d=dd2+1,dd
              do a=da2+1,da
                do b=0,bcntr
                  array_out(bbf+b,fa+a,fd+d,fc+c)=pre2*array_out(bbf+b,fa+a,fd+d,fc+c)&
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
          bcf = fc + bc
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(fb+b,fa+a,fd+d,bcf+c)=pre2*array_out(fb+b,fa+a,fd+d,bcf+c)&
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
          bdf = fd + bd
     
          do c=dc2+1,dc
            do d=0,bcntr
              do a=da2+1,da
                do b=db2+1,db
                  array_out(fb+b,fa+a,bdf+d,fc+c)=pre2*array_out(fb+b,fa+a,bdf+d,fc+c)&
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
                array_out(fb+b,fa+a,fd+d,fc+c)=pre2*array_out(fb+b,fa+a,fd+d,fc+c)&
                                             &+array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    elseif (pre2 /= 0.0E0_realk .and. pre1 /= 1.0E0_realk) then
      !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,c,d,ba,bb,bc,bd,baf,bbf,bcf,bdf),SHARED(array_in,array_out,&
      !$OMP& !da,db,dc,dd,da2,db2,dc2,dd2,bcntr,moda,modb,modc,modd,pre1,pre2,bs,fa,fb,fc,fd)
      if(da2>0.and.db2>0.and.dc2>0.and.dd2>0)then
        !$OMP DO
        do bc=1,dc2,bs
          bcf = fc + bc
          do bd=1,dd2,bs
            bdf = fd + bd
            do ba=1,da2,bs
              baf = fa + ba
              do bb=1,db2,bs
                bbf = fb + bb
     
                do c=0,bcntr
                  do d=0,bcntr
                    do a=0,bcntr
                      do b=0,bcntr
                        array_out(bbf+b,baf+a,bdf+d,bcf+c)=pre2*array_out(bbf+b,baf+a,bdf+d,bcf+c)&
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
          bcf = fc + bc
          do bd=1,dd2,bs
            bdf = fd + bd
            do bb=1,db2,bs
              bbf = fb + bb
     
              do c=0,bcntr
                do d=0,bcntr
                  do a=da2+1,da
                    do b=0,bcntr
                      array_out(bbf+b,fa+a,bdf+d,bcf+c)=pre2*array_out(bbf+b,fa+a,bdf+d,bcf+c)&
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
          bdf = fd + bd
          do bc=1,dc2,bs
            bcf = fc + bc
            do ba=1,da2,bs
              baf = fa + ba
     
              do d=0,bcntr
                do c=0,bcntr
                  do b=db2+1,db
                    do a=0,bcntr
                      array_out(fb+b,baf+a,bdf+d,bcf+c)=pre2*array_out(fb+b,baf+a,bdf+d,bcf+c)&
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
          bdf = fd + bd
          do ba=1,da2,bs
            baf = fa + ba
            do bb=1,db2,bs
              bbf = fb + bb
     
              do c=dc2+1,dc
                do d=0,bcntr
                  do a=0,bcntr
                    do b=0,bcntr
                      array_out(bbf+b,baf+a,bdf+d,fc+c)=pre2*array_out(bbf+b,baf+a,bdf+d,fc+c)&
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
          bcf = fc + bc
          do bb=1,db2,bs
            bbf = fb + bb
            do ba=1,da2,bs
              baf = fa + ba
     
              do d=dd2+1,dd
                do c=0,bcntr
                  do b=0,bcntr
                    do a=0,bcntr
                      array_out(bbf+b,baf+a,fd+d,bcf+c)=pre2*array_out(bbf+b,baf+a,fd+d,bcf+c)&
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
          bdf = fd + bd
          do bc=1,dc2,bs
            bcf = fc + bc
     
            do d=0,bcntr
              do c=0,bcntr
                do b=db2+1,db
                  do a=da2+1,da
                    array_out(fb+b,fa+a,bdf+d,bcf+c)=pre2*array_out(fb+b,fa+a,bdf+d,bcf+c)&
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
          bdf = fd + bd
          do bb=1,db2,bs
            bbf = fb + bb
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=0,bcntr
                  do a=da2+1,da
                    array_out(bbf+b,fa+a,bdf+d,fc+c)=pre2*array_out(bbf+b,fa+a,bdf+d,fc+c)&
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
          bcf = fc + bc
          do bb=1,db2,bs
            bbf = fb + bb
     
            do c=0,bcntr
              do d=dd2+1,dd
                do a=da2+1,da
                  do b=0,bcntr
                    array_out(bbf+b,fa+a,fd+d,bcf+c)=pre2*array_out(bbf+b,fa+a,fd+d,bcf+c)&
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
          bdf = fd + bd
          do ba=1,da2,bs
            baf = fa + ba
     
            do d=0,bcntr
              do c=dc2+1,dc
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(fb+b,baf+a,bdf+d,fc+c)=pre2*array_out(fb+b,baf+a,bdf+d,fc+c)&
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
          bcf = fc + bc
          do ba=1,da2,bs
            baf = fa + ba
     
            do d=dd2+1,dd
              do c=0,bcntr
                do b=db2+1,db
                  do a=0,bcntr
                    array_out(fb+b,baf+a,fd+d,bcf+c)=pre2*array_out(fb+b,baf+a,fd+d,bcf+c)&
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
          baf = fa + ba
          do bb=1,db2,bs
            bbf = fb + bb
     
            do c=dc2+1,dc
              do d=dd2+1,dd
                do a=0,bcntr
                  do b=0,bcntr
                    array_out(bbf+b,baf+a,fd+d,fc+c)=pre2*array_out(bbf+b,baf+a,fd+d,fc+c)&
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
          baf = fa + ba
     
          do d=dd2+1,dd
            do c=dc2+1,dc
              do b=db2+1,db
                do a=0,bcntr
                  array_out(fb+b,baf+a,fd+d,fc+c)=pre2*array_out(fb+b,baf+a,fd+d,fc+c)&
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
          bbf = fb + bb
     
          do c=dc2+1,dc
            do d=dd2+1,dd
              do a=da2+1,da
                do b=0,bcntr
                  array_out(bbf+b,fa+a,fd+d,fc+c)=pre2*array_out(bbf+b,fa+a,fd+d,fc+c)&
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
          bcf = fc + bc
     
          do d=dd2+1,dd
            do c=0,bcntr
              do b=db2+1,db
                do a=da2+1,da
                  array_out(fb+b,fa+a,fd+d,bcf+c)=pre2*array_out(fb+b,fa+a,fd+d,bcf+c)&
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
          bdf = fd + bd
     
          do c=dc2+1,dc
            do d=0,bcntr
              do a=da2+1,da
                do b=db2+1,db
                  array_out(fb+b,fa+a,bdf+d,fc+c)=pre2*array_out(fb+b,fa+a,bdf+d,fc+c)&
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
                array_out(fb+b,fa+a,fd+d,fc+c)=pre2*array_out(fb+b,fa+a,fd+d,fc+c)&
                                             &+pre1*array_in(a,b,c,d)
              enddo
            enddo
          enddo
        enddo
      endif
    else
      call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
    endif

    print *,"check 2143"
    if(pre2==0.0E0_realk)then
      do b=1,db
        do d=1,dd
          do a=1,da
            do c=1,dc
              if(array_out(fb+b,fa+a,fd+d,fc+c)/=pre1*array_in(a,b,c,d))then
                print *,"1432 reordering not correct",a,b,c,d,da,db,dc,dd
                print *,array_out(fb+b,fa+a,fd+d,fc+c),array_in(a,b,c,d)
                stop 0
              endif
            enddo
          enddo
        enddo
      enddo
    endif
  end subroutine manual_2143_reordering_tile2full
end module manual_utils_module
