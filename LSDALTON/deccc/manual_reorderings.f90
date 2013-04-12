!\> \brief this module is inteded to contain high performance reorderings for
!mutlidimensional arrays.
! first part: 4d arrays
! next part: 3d arrays
!\> \author Patrick Ettenhuber & Janus Juul Eriksen
!\> \date November 2012
module manual_reorderings_module
  use precision
  use dec_typedef_module
  contains
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
      call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
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
      call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
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
      call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
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
      call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
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
      call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
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
      call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
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
      call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
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
      call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
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
      call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
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
      call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
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
      call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
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
      call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
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
      call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
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
      call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
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
      call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
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
      call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
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
      call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
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
      call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
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
      call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
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
      call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
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
      call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
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
      call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
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
      call lsquit("MANUAL REORDERING:Case not found",DECinfo%output)
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
