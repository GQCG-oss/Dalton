! -------------------------------------------------------------------------
!
! Program:      Interest
!
! File:         module_interest_selftest.f90 
!
! Description:  Performing selftest calculations of all integrals implemented in src/
!
! Licensing:    This code is distributed under the GNU LGPL license
!
! Author:       Michal Repisky (michal.repisky@uit.no)
!
! Revisions:    
!
! -------------------------------------------------------------------------
module module_interest_selftest

  implicit none

  public interest_selftest_overlap
  public interest_selftest_dipole
  public interest_selftest_nuclear_attraction_point_nucleus 

  private

  integer, parameter :: lmx1 = 6
  integer, parameter :: lmx2 = 2*lmx1+1 

  real(8), parameter :: pi   = 4.0d0*datan(1.0d0)
!  real(8), parameter :: pi34 = (0.5d0/pi)**0.75d0 
  real(8), parameter :: pi34 = 0.25197943553838073021D0
  
  integer, parameter :: ncc(lmx1) = (/1,3,6,10,15,21/)

  integer, parameter, dimension(lmx1)   :: lshell = (/1,2,3,4,5,6   /)
  real(8), parameter, dimension(lmx1)   :: zeta_a = (/8.630d0,8.630d0,8.630d0,8.630d0,8.630d0,8.630d0/)
  real(8), parameter, dimension(lmx1)   :: zeta_b = (/2.910d0,2.910d0,2.910d0,2.910d0,2.910d0,2.910d0/)

  real(8), parameter, dimension(lmx1,3) :: cent_a = RESHAPE((/3.089244D0,3.089244D0,3.089244D0,3.089244D0,3.089244D0,3.089244D0,&
                                                      -0.17483D0,-0.17483D0,-0.17483D0,-0.17483D0,-0.17483D0,-0.17483D0,&
                                                      1.4532D0,1.4532D0,1.4532D0,1.4532D0,1.4532D0,1.4532D0/),(/lmx1,3/))

  real(8), parameter, dimension(lmx1,3) :: cent_b = RESHAPE((/5.036335D0,5.036335D0,5.036335D0,5.036335D0,5.036335D0,5.036335D0,&
                                                      -1.27970D0,-1.27970D0,-1.27970D0,-1.27970D0,-1.27970D0,-1.27970D0,&
                                                      3.4782D0,3.4782D0,3.4782D0,3.4782D0,3.4782D0,3.4782D0/),(/lmx1,3/))

  integer, parameter :: ncentr = 2
  integer, parameter :: natomc(ncentr) = (/17,39/)
  real(8), parameter, dimension(3,ncentr) :: cent = &
     &RESHAPE((/5.036335D0,-1.27970D0,3.4782D0,5.036335D0,-1.27970D0,3.4782D0/),(/3,ncentr/))

contains

! -------------------------------------------------------------------------
  subroutine  interest_selftest_overlap

    integer :: k
    integer :: i,li
    integer :: j,lj
    integer :: ij_batch
    real(8) :: value_max 
    real(8) :: alpha,cntra,ax,ay,az
    real(8) :: beta ,cntrb,bx,by,bz
    real(8) :: gout(ncc(lmx1)*ncc(lmx1)*10)
    real(8) :: gref(ncc(lmx1)*ncc(lmx1)*10)

    write(6,*)
    write(6,*)
    write(6,'(1x,a)') '*** SELFTEST OF OVERLAP_RKB INTEGRALS ***'
    write(6,*)

    open(unit=5,file='reference/overlap_rkb.ref',status='old',form='formatted')
    rewind(5)

    value_max = 0.0d0

    do j=1,lmx1
      lj    = lshell(j)
      beta  = zeta_b(j)
      cntrb = (4.0d0*beta)**(0.5d0*lj+0.25d0)*pi34
      bx    = cent_b(j,1)
      by    = cent_b(j,2)
      bz    = cent_b(j,3)

      do i=j,lmx1
        li    = lshell(i)
        alpha = zeta_a(i)
        cntra = (4.0d0*alpha)**(0.5d0*li+0.25d0)*pi34
        ax    = cent_a(i,1)
        ay    = cent_a(i,2)
        az    = cent_a(i,3)

        call interest_overlap(gout,ij_batch,          &
                              li,alpha,ax,ay,az,cntra,&
                              lj,beta, bx,by,bz,cntrb,&
                              .true.)         

!       do k=1,ij_batch
!           write(6,'(2(2x,i3),5x,3(2x,i4),8x,1(5x,d20.10))') i,j,li,lj,k,gout(k)
!       enddo

!       do k=1,ij_batch
!         write(6,'(d30.16)') gout(k)
!       enddo

        do k=1,ij_batch
          read(5,'(d30.16)') gref(k) 
        enddo
        do k=1,ij_batch
          if((abs(gref(k)-gout(k))/abs(gout(k))) > 1.d-10)then
            write(6,'(2(2x,i3),5x,3(2x,i4),8x,3(5x,d20.10))') i,j,li,lj,k,gref(k),gout(k),(abs(gref(k)-gout(k)))/abs(gout(k))
            value_max = max( value_max,(abs(gref(k)-gout(k)))/abs(gout(k)) )
          endif
        enddo

      enddo
    enddo 

    if( value_max < 1.d-10 )then
      write(6,*)
      write(6,'(1x,a,d20.14)') 'with the result: --passed-- with the relative error = ',value_max 
      write(6,*)
    else
      write(6,*)
      write(6,'(1x,a,d20.14)') 'with the result: --error-- with the relative error = ',value_max 
      write(6,*)
    endif

    close(5,status='keep')
  end subroutine

! -------------------------------------------------------------------------
  subroutine  interest_selftest_dipole

    integer :: k,l
    integer :: i,li
    integer :: j,lj
    integer :: ij_batch
    real(8) :: value_max(3) 
    real(8) :: alpha,cntra,ax,ay,az
    real(8) :: beta ,cntrb,bx,by,bz
    real(8) :: gout(ncc(lmx1)*ncc(lmx1)*10,3)
    real(8) :: gref(ncc(lmx1)*ncc(lmx1)*10)

    write(6,*)
    write(6,*)
    write(6,'(1x,a)') '*** SELFTEST OF DIPOLE_RKB INTEGRALS ***'
    write(6,*)

    open(unit=5,file='reference/dipole_rkb.ref',status='old',form='formatted')
    rewind(5)

    value_max = 0.0d0

    do j=1,lmx1
      lj    = lshell(j)
      beta  = zeta_b(j)
      cntrb = (4.0d0*beta)**(0.5d0*lj+0.25d0)*pi34
      bx    = cent_b(j,1)
      by    = cent_b(j,2)
      bz    = cent_b(j,3)

      do i=j,lmx1
        li    = lshell(i)
        alpha = zeta_a(i)
        cntra = (4.0d0*alpha)**(0.5d0*li+0.25d0)*pi34
        ax    = cent_a(i,1)
        ay    = cent_a(i,2)
        az    = cent_a(i,3)

        call interest_dipole(gout(1,1),              &
                             gout(1,2),              &
                             gout(1,3),ij_batch,     &
                             li,alpha,ax,ay,az,cntra,&
                             lj,beta, bx,by,bz,cntrb,&
                             0.5d0,1.0d0,1.5d0,      &
                             .true.)

!       do l=1,3
!         do k=1,ij_batch
!           write(6,'(2(2x,i3),5x,3(2x,i4),8x,1(5x,d20.10))') i,j,li,lj,k,gout(k,l)
!         enddo
!       enddo

!       do l=1,3
!         do k=1,ij_batch
!           write(6,'(d30.16)') gout(k,l)
!         enddo
!       enddo

        do l=1,3
          do k=1,ij_batch
            read(5,'(d30.16)') gref(k) 
          enddo
          do k=1,ij_batch
            if((abs(gref(k)-gout(k,l))/abs(gout(k,l))) > 1.d-10)then
              write(6,'(2(2x,i3),5x,3(2x,i4),8x,3(5x,d20.10))') &
     &i,j,li,lj,k,gref(k),gout(k,l),(abs(gref(k)-gout(k,l)))/abs(gout(k,l))
              value_max(l) = max( value_max(l),(abs(gref(k)-gout(k,l)))/abs(gout(k,l)) )
            endif
          enddo
        enddo

      enddo
    enddo 

    if( maxval(value_max) < 1.d-10 )then
      write(6,*)
      write(6,'(1x,a,3(d20.14,2x))') 'with the result: --passed-- with the relative errors = ',(value_max(l),l=1,3) 
      write(6,*)
    else
      write(6,*)
      write(6,'(1x,a,3(d20.14,2x))') 'with the result: --error-- with the relative errors = ',(value_max(l),l=1,3)
      write(6,*)
    endif

    close(5,status='keep')
  end subroutine

! -------------------------------------------------------------------------
  subroutine interest_selftest_nuclear_attraction_point_nucleus 

    integer :: k
    integer :: i,li
    integer :: j,lj
    integer :: ij_batch
    real(8) :: value_max 
    real(8) :: alpha,cntra,ax,ay,az
    real(8) :: beta ,cntrb,bx,by,bz
    real(8) :: gout(ncc(lmx1)*ncc(lmx1)*10)
    real(8) :: gref(ncc(lmx1)*ncc(lmx1)*10)

    write(6,*)
    write(6,*)
    write(6,'(1x,a)') '*** SELFTEST OF NUCLEAR_ATTRACTION_RKB_POINT_NUCLEUS INTEGRALS ***'
    write(6,*)

    open(unit=5,file='reference/nuclear_attraction_point_nucleus.ref',status='old',form='formatted')
    rewind(5)

    value_max = 0.0d0

    do j=1,lmx1
      lj    = lshell(j)
      beta  = zeta_b(j)
      cntrb = (4.0d0*beta)**(0.5d0*lj+0.25d0)*pi34
      bx    = cent_b(j,1)
      by    = cent_b(j,2)
      bz    = cent_b(j,3)

      do i=j,lmx1
        li    = lshell(i)
        alpha = zeta_a(i)
        cntra = (4.0d0*alpha)**(0.5d0*li+0.25d0)*pi34
        ax    = cent_a(i,1)
        ay    = cent_a(i,2)
        az    = cent_a(i,3)

        call interest_nuclear_attraction_point_nucleus(gout,ij_batch,          &
                                                       li,alpha,ax,ay,az,cntra,&
                                                       lj,beta, bx,by,bz,cntrb,&
                                                       ncentr,natomc,cent,     &
                                                       .true.)

!       do k=1,ij_batch
!           write(6,'(2(2x,i3),5x,3(2x,i4),8x,1(5x,d20.10))') i,j,li,lj,k,gout(k)
!       enddo

!       do k=1,ij_batch
!         write(6,'(d30.16)') gout(k)
!       enddo

        do k=1,ij_batch
          read(5,'(d30.16)') gref(k) 
        enddo
        do k=1,ij_batch
          if((abs(gref(k)-gout(k))/abs(gout(k))) > 1.d-10)then
            write(6,'(2(2x,i3),5x,3(2x,i4),8x,3(5x,d20.10))') i,j,li,lj,k,gref(k),gout(k),(abs(gref(k)-gout(k)))/abs(gout(k))
            value_max = max( value_max,(abs(gref(k)-gout(k)))/abs(gout(k)) )
          endif
        enddo

      enddo
    enddo 

    if( value_max < 1.d-10 )then
      write(6,*)
      write(6,'(1x,a,d20.14)') 'with the result: --passed-- with the relative error = ',value_max 
      write(6,*)
    else
      write(6,*)
      write(6,'(1x,a,d20.14)') 'with the result: --error-- with the relative error = ',value_max 
      write(6,*)
    endif

    close(5,status='keep')
  end subroutine

! -------------------------------------------------------------------------
end module
