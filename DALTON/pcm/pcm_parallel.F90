module pcm_parallel

implicit none

public j1intp_pcm
public j1ints_pcm

private

logical, public :: pcm_parallel_initialized = .false.

contains

   subroutine j1intp_pcm(nbast, nosim, ksymp, dmat, exp1vl, expval, &
                         nr_points, nr_points_irr, centers, tofile, &
                         iprtyp, matdim, work, lwork)
   !
   !     Master routine for distributing tesseraes to slaves
   !     K.Ruud, June 10 2005, Pisa
   !
#include "maxorb.h"
#include "mxcent.h"
#include "priunit.h"
#include "mtags.h"
#include "infpar.h"
#include "mpif.h"

   ! Passed variables
   integer(8) :: nr_points, nr_points_irr
   integer    :: nbast, nosim, ksymp, lwork, matdim
   real(8)    :: centers(3, nr_points), expval(nr_points, nosim)
   real(8)    :: dmat(matdim, nosim)
   real(8)    :: work(lwork)
   logical    :: exp1vl, tofile
   integer    :: iprtyp

   ! Local variables
   integer    :: iprpcm = 0
   integer    :: its, iwho, nwho, islave, ierr, nnbasx

   if (tofile) call quit('Parallel calculations do not allow for storing PCM-integrals on disk')

   !
   !     Wake up the slaves
   !
   nnbasx = nbast*(nbast + 1)/2
   call mpixbcast(iprtyp, 1, 'INTEGER', master)
   call mpixbcast(iprpcm, 1, 'INTEGER', master)

   call mpixbcast(ksymp, 1, 'INTEGER', master)
   call mpixbcast(exp1vl, 1, 'LOGICAL', master)
   call mpixbcast(nbast, 1, 'INTEGER', master)
   call mpixbcast(nr_points_irr, 1, 'INTEGER', master)
   call mpixbcast(nosim, 1, 'INTEGER', master)
   if (exp1vl) then
      call mpixbcast(dmat, nnbasx, 'DOUBLE', master)
   else
      call mpixbcast(expval, nr_points*nosim, 'DOUBLE', master)
      call dzero(dmat, matdim*nosim)
   end if
   call mpixbcast(centers(1, :), nr_points_irr, 'DOUBLE', master)
   call mpixbcast(centers(2, :), nr_points_irr, 'DOUBLE', master)
   call mpixbcast(centers(3, :), nr_points_irr, 'DOUBLE', master)
   call mpixbcast(tofile, 1, 'LOGICAL', master)
   !
   !  Loop over all tesserae
   !
   do its = 1, nr_points_irr
      iwho = -1
      call mpixrecv(nwho, 1, 'INTEGER', iwho, mptag1)
      call mpixsend(its, 1, 'INTEGER', nwho, mptag2)
   end do
   !
   !  Send end message to all slaves
   !
   its = -1
   do islave = 1, nodtot
      iwho = -1
      call mpixrecv(nwho, 1, 'INTEGER', iwho, mptag1)
      call mpixsend(its , 1, 'INTEGER', nwho, mptag2)
   end do
   !
   !  Collect data from all slaves
   !
   if (exp1vl) then
      call dzero(work, nr_points)
      call mpi_reduce(work, expval, nr_points, mpi_double_precision,  &
                      mpi_sum, 0, mpi_comm_world, ierr)
   else if (.not. tofile) then
   !
   !     work(1) = dmat(1)
   !
      call dzero(work(nosim*matdim+1), matdim*nosim)
      call mpi_reduce(work(nosim*matdim+1), dmat, matdim*nosim,  &
                      mpi_double_precision, mpi_sum, 0, mpi_comm_world,  &
                      ierr)
   end if

   end subroutine j1intp_pcm

   subroutine j1ints_pcm(work, lwork, iprtyp, iprtmp)

#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "iratdef.h"
#include "priunit.h"
#include "orgcom.h"
#include "symmet.h"
#include "infpar.h"
#include "mtags.h"
#include "mpif.h"

      ! Passed variables
      integer :: lwork, iprtyp, iprtmp
      real(8) :: work(lwork)
      ! Local variables
      integer(8) :: nr_points, nr_points_irr
      logical :: tofile, exp1vl, trimat
      character(7) :: intlab
      character(8) :: labint(3*mxcoor)
      integer      :: i, iadr, ierr, iosim, iprpcm
      integer      :: its, jadr, kadr, kden, kexpvl, kintad
      integer      :: kintrp, klast, klst1, kmat, kpatom
      integer      :: ksymp, ktmp, l, lwrk, madr, matdim
      integer      :: nnbasx, n2basx, ncomp, ntesp
      integer      :: nosim, nbast
      real(8), allocatable :: x_centers(:), y_centers(:), z_centers(:)


      iprpcm = iprtmp
      call mpixbcast(ksymp, 1, 'INTEGER', master)
      call mpixbcast(exp1vl, 1, 'LOGICAL', master)
      call mpixbcast(nbast, 1, 'INTEGER', master)
      call mpixbcast(nr_points_irr, 1, 'INTEGER', master)
      call mpixbcast(nosim, 1, 'INTEGER', master)
      nnbasx = nbast*(nbast + 1)/2
      n2basx = nbast*nbast
      if (iprtyp .eq. 133) then
         matdim = nnbasx
         trimat = .true.
         ncomp  = 1
         intlab = 'NPETES '
      else
         matdim = n2basx
         trimat = .false.
         ncomp = 3
         intlab = 'PCMBSOL'
      end if

      nr_points = nr_points_irr * (maxrep + 1)
      allocate(x_centers(nr_points_irr))
      x_centers = 0.0d0
      allocate(y_centers(nr_points_irr))
      y_centers = 0.0d0
      allocate(z_centers(nr_points_irr))
      z_centers = 0.0d0

      kden   = 1
      kintrp = kden + matdim*nosim
      kintad = kintrp + (3*mxcoor + 1)/irat
      kexpvl = kintad + (3*mxcoor + 1)/irat
      klst1  = kexpvl + nr_points_irr*(maxrep + 1)*nosim

      if (exp1vl) then
         call mpixbcast(work(kden), nnbasx, 'DOUBLE', master)
         call dzero(work(kexpvl), nr_points_irr*(maxrep + 1))
      else
         call mpixbcast(work(kexpvl), nr_points_irr*(maxrep + 1)*nosim, 'DOUBLE', master)
         call dzero(work(kden), matdim*nosim)
      end if
      call mpixbcast(x_centers(1:), nr_points_irr, 'DOUBLE', master)
      call mpixbcast(y_centers(1:), nr_points_irr, 'DOUBLE', master)
      call mpixbcast(z_centers(1:), nr_points_irr, 'DOUBLE', master)
      call mpixbcast(tofile, 1, 'LOGICAL', master)
      !
      !     Loop over tesseraes to be calculated on this node
      !
 10   continue
      call mpixsend(mynum, 1, 'INTEGER', master, mptag1)
      call mpixrecv(its, 1, 'INTEGER', master, mptag2)

      if (its .gt. 0) then
         diporg(1) = x_centers(its)
         diporg(2) = y_centers(its)
         diporg(3) = z_centers(its)
         ntesp = 1
         kpatom = 0
         !
         !  Calculates nuclear potential energy integrals (in AO basis) for
         !  the given tessera
         !
         l=1
         ktmp = klst1
         if (.not. tofile .and. .not. exp1vl) then
            kmat = ktmp + 8
            if (iprtyp .eq. 133) then
               klast = kmat + (maxrep + 1)*matdim
            else
               klast = kmat + (maxrep + 1)*matdim*ncomp
            end if
            ncomp = (maxrep + 1)
         else
            kmat = ktmp + 8
            klast = kmat
            ncomp = 0
         end if
         lwrk = lwork - klast + 1
         call get1in(work(kmat), intlab, ncomp, work(klast), lwrk, labint,  &
                     work(kintrp), work(kintad), l, tofile, kpatom, trimat, &
                     work(ktmp), exp1vl, work(kden), iprpcm)
         if (iprtyp .eq. 134) then
            do iosim = 1,  nosim
               iadr = kexpvl + its - 1 + nr_points*(iosim - 1)
               jadr = kmat + (iosim - 1)*matdim
               madr = kden + (iosim - 1)*matdim
               call daxpy(matdim, work(iadr), work(jadr), 1, &
                          work(madr), 1)
            end do
         else if (exp1vl) then
            do i = 1,  ncomp
               iadr = kexpvl + its - 1 + (i-1)*nr_points_irr
               work(iadr) = -work(ktmp+i-1)
            end do
         else
            do iosim = 1,  nosim
               iadr = kexpvl + its - 1 + nr_points_irr*(maxrep + 1)*(iosim - 1)
               kadr = kmat + (ksymp - 1)*matdim
               madr = kden + (iosim - 1)*matdim
               call daxpy(matdim, -work(iadr), work(kadr), 1, &
                          work(madr), 1)
            end do
         end if
         go to 10
      end if
      !
      !   No more tesseraes to calculate
      !
      if (exp1vl) then
         call mpi_reduce(work(kexpvl), mpi_in_place, nr_points_irr*(maxrep + 1),  &
                         mpi_double_precision, mpi_sum, 0, mpi_comm_world, &
                         ierr)
      else if (.not. tofile) then
         call mpi_reduce(work(kden), mpi_in_place, matdim*nosim,           &
                         mpi_double_precision, mpi_sum, 0, mpi_comm_world, &
                         ierr)
      end if

      deallocate(x_centers)
      deallocate(y_centers)
      deallocate(z_centers)

   end subroutine j1ints_pcm

end module pcm_parallel
