
!------------------------------------------------------------------------------

      subroutine get_Tk_integrals(Tk_ints, nints, k, coord, work, nwrk)

#include <mxcent.h>
#include <qm3.h>
#include <orgcom.h>

        character(8), dimension(9*mxcent) :: labint
        integer*8, dimension(9*mxcent) :: intrep, intadr

        integer, intent(in) :: k, nints, nwrk
        real*8, dimension(nints), intent(inout) :: Tk_ints
        real*8, dimension(3), intent(in) :: coord
        real*8, dimension(nwrk), intent(inout) :: work

        logical :: trimat
        integer :: i, j, m
        character(7) :: inttype
        real*8, dimension(3) :: backup

        backup = diporg
        diporg = coord
        runqm3 = .true.

        trimat = .false.

        if (k == 0) then
            ncomp = 1
            inttype = 'NPETES '
        else if (k == 1) then
            ncomp = 3
            inttype = 'NEFIELD'
        else if (k == 2) then
            ncomp = 6
            inttype = 'ELFGRDC'
        else if (k == 3) then
            if (nwrk < (24 * nints / 10)) then
                print *, 'Not enough work space for T^(3) integrals!'
            end if
            ncomp = 6
            inttype = 'ELFGRDC'
        end if

        Tk_ints = 0.0d0

        if (k <= 2) then
            call get1in(Tk_ints(1), inttype, ncomp, work(1), nwrk,
     &                  labint, intrep, intadr, -1, .false., 0, trimat,
     &                  work, .false., work, 1)
        else if (k == 3) then
            m = nints / 10
            i = 6 * nints / 10
            l = 0
            do j = 1, 3
                diporg(j) = diporg(j) + 0.01d0
                call get1in(work(l*i+1), inttype, ncomp,
     &                      work(j*i+1), nwrk, labint, intrep,
     &                      intadr, -1, .false., 0, trimat, work,
     &                      .false., work, 1)
                diporg(j) = diporg(j) - 2.0d0 * 0.01d0
                call get1in(work(j*i+1), inttype, ncomp,
     &                      work((j+1)*i+1), nwrk, labint, intrep,
     &                      intadr, -1, .false., 0, trimat, work,
     &                      .false., work, 1)
                diporg(j) = coord(j)
                work(l*i+1:j*i) = (work(l*i+1:j*i)
     &                            - work(j*i+1:(j+1)*i))
     &                            / (2.0d0 * 0.01d0)
                l = l + 1
            end do
            Tk_ints(1:m) = work(1:m)
            Tk_ints(m+1:2*m) = work(i+1:i+1+m)
            Tk_ints(2*m+1:3*m) = work(2*i+1:2*i+1+m)
            Tk_ints(3*m+1:4*m) = work(i+1+m:i+1+2*m)
            Tk_ints(4*m+1:5*m) = work(2*i+1+m:2*i+1+2*m)
            Tk_ints(5*m+1:6*m) = work(2*i+1+2*m:2*i+1+3*m)
            Tk_ints(6*m+1:7*m) = work(i+1+3*m:i+1+4*m)
            Tk_ints(7*m+1:8*m) = work(2*i+1+3*m:2*i+1+4*m)
            Tk_ints(8*m+1:9*m) = work(2*i+1+4*m:2*i+1+5*m)
            Tk_ints(9*m+1:10*m) = work(2*i+1+5*m:2*i+1+6*m)
            Tk_ints = -1.0d0 * Tk_ints
        end if

        runqm3 = .false.
        diporg = backup

      end subroutine get_Tk_integrals


