module pcm_integrals

use pcm_config, only: pcm_configuration, pcm_cfg

implicit none

public get_mep
public get_nuclear_mep
public get_electronic_mep
public pot_int_tess
public j1x_pcm

contains

   subroutine get_mep(nr_points, nr_points_irr, centers, mep, density, work, lwork)

   integer(8), intent(in)  :: nr_points
   integer(8), intent(in)  :: nr_points_irr
   real(8), intent(in)  :: centers(3, nr_points)
   real(8), intent(out) :: mep(nr_points)
   real(8)              :: density(*)
   real(8)              :: work(*)
   integer              :: lwork

   ! For some reason that I don't understand, the potential vector is zeroed out
   ! somewhere in j1int_pcm. So the mep must be calculated in this order.
   call get_electronic_mep(nr_points, nr_points_irr, centers, mep, density, work, lwork, .false.)
   call get_nuclear_mep(nr_points, nr_points_irr, centers, mep)

   end subroutine get_mep

   subroutine get_nuclear_mep(nr_points, nr_points_irr, centers, mep)

#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "symmet.h"
#include "nuclei.h"

   integer(8), intent(in)  :: nr_points
   integer(8), intent(in)  :: nr_points_irr
   real(8), intent(in)  :: centers(3, nr_points)
   real(8), intent(out) :: mep(nr_points)

   real(8)              :: coora(3, nucdep), charges(nucdep)
   real(8), allocatable :: scratch(:)
   real(8)              :: dist, temp
   integer              :: i, j, k, ipoint
   real(8)              :: renorm

   renorm = dble(maxrep + 1)
   ! Get coordinates and charges of all centers, not only the ones
   ! that are symmetry independent
   call getacord(coora)
   i = 0
   do j = 1, nucind
      do k = 1, nucdeg(j)
         i = i + 1
         charges(i) = charge(j)
      enddo
   enddo

   ! In order to renormalize for the number of symmetry operations we have
   ! to branch based on the type of algorithms chosen: .SEPARATE or
   ! .TOTAL MEP/ASC
   ! In the second case, in fact, the nuclear potential is accumulated on
   ! top of the electronic and the renormalization will then affect also
   ! the electronic part (which we do not want)
   if (pcm_cfg%separate) then
     do i = 1, nucdep
        do ipoint = 1, nr_points_irr
           dist = (coora(1, i) - centers(1, ipoint))**2 +      &
                  (coora(2, i) - centers(2, ipoint))**2 +      &
                  (coora(3, i) - centers(3, ipoint))**2
           dist = sqrt(dist)
           mep(ipoint) = mep(ipoint) + (charges(i) / dist)
        end do
     end do

     do ipoint = 1, nr_points_irr
        temp = mep(ipoint) * renorm
        mep(ipoint) = temp
        !write(lupri, *) "v_nuc(",ipoint,") = ", mep(ipoint)
     end do
   else
     allocate(scratch(nr_points_irr))
     scratch = 0.0d0
     do i = 1, nucdep
        do ipoint = 1, nr_points_irr
           dist = (coora(1, i) - centers(1, ipoint))**2 +      &
                  (coora(2, i) - centers(2, ipoint))**2 +      &
                  (coora(3, i) - centers(3, ipoint))**2
           dist = sqrt(dist)
           scratch(ipoint) = scratch(ipoint) + (charges(i) / dist)
        end do
     end do

     do ipoint = 1, nr_points_irr
        temp = scratch(ipoint) * renorm
        scratch(ipoint) = temp
        !write(lupri, *) "v_nuc(",ipoint,") = ", scratch(ipoint)
     end do
     do ipoint = 1, nr_points_irr
       mep(ipoint) = mep(ipoint) + scratch(ipoint)
     end do
     deallocate(scratch)
   end if

   !print *, "Debug print of v_nuc at cavity points"
   !do ipoint = 1, nr_points
   !   print *, "v_nuc(",ipoint,") = ", mep(ipoint)
   !end do

   end subroutine get_nuclear_mep

   subroutine get_electronic_mep(nr_points, nr_points_irr, centers, vector, matrix, work, lwork, get_matrix)
   !
   ! Driver routine for the calculation of the electronic part of the molecular
   ! electrostatic potential on a certain grid of points {r_i}:
   !     v_el(r_i) = tr(DV_i)
   ! tr is the trace operator, D is the density matrix, V^i is the matrix of the
   ! "nuclear attraction" integrals calculated at point r_i of the grid:
   !     V_mu,nu,i =  - <mu|1/|r-r_i||nu>
   !
   ! Written, tested, debugged: R. Bast, R. Di Remigio, L. Frediani, K. Ruud
   !
   !  array ADER(mu, nu, r_i) contains these integrals.
   !  array vc(r_i) contains the electronic MEP at point r_i.
   !
   ! RDR 0512.
   !
   ! RDR 050312 CANNOT yet handle symmetry.
   ! RDR 220312 This routine will be used to form both potentials and Fock
   !            matrix contribution for PCM.
   !            matrix is Fock or density matrix, vector is potentials
   !            or charges vector.
   !            get_matrix logical is present and TRUE:
   !            charges vector as input, Fock matrix contribution as output.
   !            get_matrix logical is absent or is present and FALSE:
   !            density matrix as input, potentials vector as output.
   !
   integer(8), intent(in)  :: nr_points
   integer(8), intent(in)  :: nr_points_irr
   real(8), intent(in)  :: centers(3, nr_points)
   real(8), intent(out) :: vector(nr_points)
   real(8)              :: matrix(*)
   real(8)              :: work(*)
   integer              :: lwork
   logical, optional    :: get_matrix
   ! Local variables
   logical              :: do_matrix
   integer :: ipoint
   integer :: nosim, ksymp
   logical :: tofile
   character(7) :: integral_type

   ! Fock matrix contribution or potential calculation?
   do_matrix = .false.
   if (present(get_matrix)) then
     if (get_matrix) then
       do_matrix = .true.
     end if
   end if

   ! Decide which integration routines to use
   if (pcm_cfg%fast_integration) then
      ksymp = 1
      call vectorized_integration_pcm(nr_points, centers, vector, matrix, ksymp, work, lwork, do_matrix)
   else
      nosim  = 1
      tofile = .false.
      integral_type = 'NPETES '
      ksymp = 1
      call j1int_pcm(vector, nr_points, nr_points_irr, centers, &
    (.not.do_matrix), matrix, nosim, tofile, integral_type, ksymp, work, lwork)
   end if

   !if (do_matrix) then
   !   print *, "Called with get_matrix"
   !else
   !   print *, "Debug print of v_ele at cavity points"
   !   do ipoint = 1, nr_points
   !     print *, "v_ele(", ipoint,") = ", vector(ipoint)
   !   end do
   !end if

   end subroutine get_electronic_mep

   subroutine j1int_pcm(expval, nr_points, nr_points_irr, centers, exp1vl, denmat, &
                        nosim, tofile, intlab, ksymp, work, lwork)
   !
   ! This subroutines is used both for the calculation of the electronic
   ! MEP and the PCM Fock matrix contribution.
   ! if (exp1vl) then
   !             density matrix as input and vector of electronic MEP as
   !             output
   ! else
   !             vector of ASC and Fock matrix contribution as output
   ! end if
   !
#if defined (VAR_MPI)
   use pcm_parallel, only: j1intp_pcm
#endif

#include "dummy.h"
#include "mxcent.h"
#include "nuclei.h"
#include "inforb.h"
#include "orgcom.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "symmet.h"
#include "infpar.h"

   ! Passed variables
   integer(8)   :: nr_points, nr_points_irr
   integer      :: nosim, ksymp, lwork
   real(8)      :: centers(3, nr_points), expval(nr_points, nosim)
   real(8)      :: denmat(*)
   real(8)      :: work(*)
   logical      :: exp1vl, tofile
   character(7) :: intlab

   ! Local variables
   character(8) :: labint(9*mxcent)
   integer      :: isum
   integer      :: i, j, iaddr, iosim, iprpcm
   integer      :: iadr, iprtyp, isymd, its, jmat, kden, klast
   integer      :: kmat, kpatom, ktmp, l, lwrk, matdim, nbastold
   integer      :: lwrk1, klast1
   integer      :: ncomp, nnbasxold, ntesp
   real(8)      :: intrep(9*mxcent), intadr(9*mxcent)
   real(8)      :: xdiporg, ydiporg, zdiporg
   logical      :: trimat, vcheck

#include "iprtyp.h"

   nbastold  = nbast
   nnbasxold = nnbasx
   nbast  = isum(maxrep + 1, naos, 1)
   nnbasx = nbast * (nbast + 1)/2
   n2basx = nbast*nbast
   if (intlab .eq. 'PCMBSOL') then
      matdim = n2basx
      iprtyp = 134
      ncomp  = 3
      trimat = .false.
   else
      matdim = nnbasx
      iprtyp = 133
      ncomp  = 1
      trimat = .true.
   end if

   !
   ! We use as a quick way of transfering tessera coordinates to hermit:
   ! the dipole origin. Need to be restored.
   !
   xdiporg = diporg(1)
   ydiporg = diporg(2)
   zdiporg = diporg(3)
   !
   ! 2) calculation of apparent charges generated by the solute's nuclei.
   !
   iprpcm=0
   if (exp1vl) then
      kden = 1
      klast = kden + nnbasx
      lwrk  = lwork - klast
      isymd = ksymp - 1
      if (ksymp .eq. 1) then
         call pksym1(denmat,work(klast),nbas,nsym,1)
         call dsym1(work(kden),dummy,work(klast),dummy,.false.,nbast,iprpcm)
      else
         call dcopy(nnbasx,denmat,1,work(kden),1)
      end if
      if (nosim .gt. 1) call quit('nosim .gt. 1 and exp1vl not permitted in j1int')
   else
      kden = 1
      if (nodtot .ge. 1) then
         klast = kden + matdim*nosim
         call dcopy(matdim*nosim,denmat,1,work(kden),1)
      else
         klast = kden
      end if
      lwrk = lwork - klast
   end if

   ! RDR 210515 Copied verbatim from sirpcm.F
#if defined (VAR_MPI)
      if (nodtot .ge. 1) then
         call j1intp_pcm(nbast, nosim, ksymp, work(kden), exp1vl, expval, &
                     nr_points, nr_points_irr, centers, tofile,       &
                     iprtyp, matdim, work(klast), lwrk)
         if (.not. exp1vl) call daxpy(matdim*nosim, 1.0d0, work(kden), 1, &
                                     denmat, 1)
      else
#endif

   do  its = 1, nr_points_irr
      diporg(1) = centers(1,its)
      diporg(2) = centers(2,its)
      diporg(3) = centers(3,its)

      ntesp = 1
      kpatom = 0
   !
   ! calculates nuclear potential energy integrals (in ao basis) for the given tessera
   !
      l=1
      ktmp = klast
      if (.not. tofile .and. .not. exp1vl) then
         kmat = ktmp + 8
         if (iprtyp .eq. 11) then
            klast1 = kmat + (maxrep + 1)*matdim
         else
            klast1 = kmat + (maxrep + 1)*matdim*ncomp
         end if
         ncomp = nsym
      else
         kmat  = ktmp + 8
         klast1 = kmat
         ncomp = 0
      end if
      lwrk1 = lwork - klast1
      call get1in(work(kmat),intlab,ncomp,work(klast1),lwrk1,labint, &
                  intrep,intadr,l,tofile,kpatom,trimat,work(ktmp), &
                  exp1vl,work(kden),iprpcm)
      if (iprtyp .eq. 13) then
         jmat = kmat
         do iosim = 1, nosim
            call daxpy(matdim,expval(its,iosim),work(jmat),1,denmat(matdim*(iosim - 1) + 1),1)
            jmat = jmat + matdim
         end do
      else if (exp1vl) then
         do i = 1, ncomp
           expval(its+(i-1)*nr_points_irr,1) = -work(ktmp+i-1)
         end do
      else if (.not. tofile) then
         do iosim = 1, nosim
            iadr = kmat + (ksymp - 1)*matdim
            call daxpy(matdim,-expval(its,iosim),work(iadr),1,denmat(matdim*(iosim - 1) + 1),1)
         end do
      end if
   enddo

   ! RDR 210515 Copied verbatim from sirpcm.F
#if defined (VAR_MPI)
      END IF
#endif

   diporg(1) = xdiporg
   diporg(2) = ydiporg
   diporg(3) = zdiporg
   nbast  = nbastold
   nnbasx = nnbasxold

   end subroutine j1int_pcm

   subroutine j1x_pcm(nr_points, nr_points_irr, centers, nosim, vtex, &
                      ucmo, ubo, udv, udvtr, jwopsy,                  &
                      work, lwork)
   ! Calculates one-index transformed potentials at tesserae centers
#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxash.h"
#include "iratdef.h"
#include "priunit.h"
#include "inforb.h"
#include "infind.h"
#include "orgcom.h"
#include "symmet.h"
#include "wrkrsp.h"
#include "infrsp.h"

   ! Parameters
   real(8), parameter :: d0 = 0.0d0
   ! Passed variables
   integer(8) :: nr_points
   integer(8) :: nr_points_irr
   real(8) :: centers(3, nr_points)
   integer :: nosim
   real(8) :: vtex(nr_points, nosim)
   real(8) :: ucmo(norbt * nbast)
   real(8) :: ubo(nosim * n2orbx)
   real(8) :: udv(n2ashx)
   real(8) :: udvtr(n2ashx)
   integer :: jwopsy
   integer :: lwork
   real(8) :: work(lwork)
   ! Local variables
   character(8) :: labint(3 * mxcoor)
   integer      :: its, iadr, iosim, iprpcm
   real(8) :: dipx, dipy, dipz
   real(8) :: factor
   logical :: exp1vl, trimat
   integer :: jj1ao, jj1x, jj1xac, jpcmx, jts, jubo
   integer :: kexpvl, kintad, kintrp, kj1, kj1ao, kj1sq, kj1xac, kj1xsq
   integer :: klast, kpatom, kpcmx, kubo, kucmo, kudv, l, lwrk, ncomp
   integer :: ndimz
   real(8) :: tj1xac, slvtlm, slvqlm
   logical :: tofile

   ndimz = max(n2orbx, nr_points)
   call dzero(work, nosim * ndimz)
   call dzero(vtex, nr_points * nosim)
   iprpcm = 0
   kintrp = 1
   kintad = kintrp + (3 * mxcoor + 1) / irat
   kexpvl = kintad + (3 * mxcoor + 1) / irat
   kj1ao  = kexpvl + nr_points_irr * (maxrep + 1) * nosim
   kj1    = kj1ao  + nnbasx * (maxrep + 1)
   kj1sq  = kj1    + nnorbx
   kj1xsq = kj1sq  + n2orbx
   kj1xac = kj1xsq + n2orbx * nosim
   kucmo  = kj1xac + n2ashx * nosim
   kudv   = kucmo  + norbt  * nbast
   kubo   = kudv   + n2ashx
   kpcmx  = kubo   + n2orbx * nosim
   klast  = kpcmx  + n2orbx * nosim
   lwrk   = lwork - klast + 1
   !
   !  Loop over tesseraes to be calculated on this node
   !
   call dzero(work(kpcmx),  nosim * n2orbx)
   call dzero(work(kexpvl), nosim * nr_points_irr * (maxrep + 1))

   jubo = 1

   dipx = diporg(1)
   dipy = diporg(2)
   dipz = diporg(3)
   IrreducibleCavityPoints: do its = 1, nr_points_irr
      l = 1
      ncomp     = maxrep + 1
      diporg(1) = centers(1, its)
      diporg(2) = centers(2, its)
      diporg(3) = centers(3, its)
      exp1vl    = .false.
      tofile    = .false.
      kpatom    = 0
      trimat    = .true.

      call get1in(work(kj1ao), 'NPETES ', ncomp, work(klast), lwrk,      &
                  labint, work(kintrp), work(kintad), l, tofile, kpatom, &
                  trimat, dummy, exp1vl, dummy, iprpcm)
      jj1ao = kj1ao
      !  Transform AO pot. int. into MO basis  V(AO) --> V(MO)
      call uthu(work(jj1ao), work(kj1), ucmo, work(klast), nbast, norbt)
      ! Transform V(MO) from triangular to square format
      call dsptsi(norbt, work(kj1), work(kj1sq))
      call dzero(work(kj1xsq), n2orbx * nosim)
      TildeVSymmetries: do iosim = 1, nosim
         jubo = 1 + (iosim - 1) * n2orbx
         jj1x = kj1xsq + (iosim - 1) * n2orbx
         jj1xac = kj1xac + (iosim - 1) * n2ashx
         jpcmx  = kpcmx + (iosim - 1)*n2orbx
         call onexh1(ubo(jubo), work(kj1sq), work(jj1x))
         if (nasht .gt. 0) call getacq(work(jj1x), work(jj1xac))
         if (iprrsp .ge. 15) then
            write (lupri,'(/a,i5)') 'J1X_mo matrix tess:', its
            call output(work(jj1x), 1, norbt, 1, norbt, norbt, norbt, 1, lupri)
            if (nasht .gt. 0) then
               write (lupri,'(/a)') ' J1X_ac matrix:'
               call output(work(jj1xac), 1, nasht, 1, nasht, nasht, nasht, 1, lupri)
            end if
         end if
         !
         ! Expectation value of transformed potential on tesserae:
         !               <0|\tilde{V}|0>
         !
         if (ksymop .eq. 1) then
            if (trplet) then
               factor = slvtlm(work(kudv), work(jj1xac), work(jj1x), tj1xac)
            else
               factor = slvqlm(work(kudv), work(jj1xac), work(jj1x), tj1xac)
            end if
            iadr = kexpvl + nr_points_irr * (maxrep + 1) * (iosim - 1) + its - 1
            vtex(its, iosim) = factor
            if (iprrsp .ge. 6) then
               write (lupri,'(a,f17.8)') ' --- J1X expectation value :', work(iadr)
               write (lupri,'(a,f17.8)') ' --- active part of J1X    :', tj1xac
            end if
         end if
      end do TildeVSymmetries
      ! Non-totally symmetric perturbation operators
      if (ksymop .gt. 1) then
         jj1ao = kj1ao + (ksymop - 1) * nnbasx
         jts = (ksymop - 1) * nr_points_irr + its
         ! Transform AO pot. int. into MO basis  V(AO) --> V(MO)
         call uthu(work(jj1ao), work(kj1), ucmo, work(klast), nbast, norbt)
         ! Transform V(MO) from triangular to square format
         call dsptsi(norbt, work(kj1), work(kj1sq))
         call dzero(work(kj1xsq), n2orbx * nosim)
         Symmetries: do iosim = 1, nosim
            jubo = 1 + (iosim - 1) * n2orbx
            jj1x = kj1xsq + (iosim - 1) * n2orbx
            jj1xac = kj1xac + (iosim - 1) * n2ashx
            call onexh1(ubo(jubo), work(kj1sq), work(jj1x))
            if (nasht .gt. 0) call getacq(work(jj1x), work(jj1xac))
            if (iprrsp .ge. 15) then
               write (lupri,'(/a)') ' J1X_mo matrix (KSYMOP):'
               call output(work(jj1x), 1, norbt, 1, norbt, norbt, norbt, 1, lupri)
               if (nasht .gt. 0) then
                  write (lupri,'(/a)') ' J1X_ac matrix (KSYMOP):'
                  call output(work(jj1xac), 1, nasht, 1, nasht, nasht, nasht, 1, lupri)
               end if
            end if
            ! Expectation value of transformed potential on tesserae:
            !         <0|\tilde{V}|0>
            !
            if (trplet) then
               factor = slvtlm(work(kudv), work(jj1xac), work(jj1x), tj1xac)
            else
               factor = slvqlm(work(kudv), work(jj1xac), work(jj1x), tj1xac)
            endif
            iadr = kexpvl + nr_points_irr * (maxrep + 1) * (iosim - 1) + jts - 1
            vtex(jts, iosim) = factor
            if (iprrsp .ge. 6) then
               write (lupri,'(a,f17.8)') ' --- J1X expectation value :', work(iadr)
               write (lupri,'(a,f17.8)') ' --- active part of J1X    :', tj1xac
            end if
         enddo Symmetries
      end if
   end do IrreducibleCavityPoints
   diporg(1) = dipx
   diporg(2) = dipy
   diporg(3) = dipz

   end subroutine j1x_pcm

   subroutine vectorized_integration_pcm(nr_points, centers, vector, matrix, ksymp, work, lwork, do_matrix)

#include "dummy.h"
#include "mxcent.h"
#include "nuclei.h"
#include "inforb.h"
#include "orgcom.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "symmet.h"
#include "infpar.h"
#include "inftap.h"
#include "aovec.h"
#include "ccom.h"
#include "shells.h"
#include "onecom.h"
#include "primit.h"
#include "lmns.h"
#include "priunit.h"

   ! Parameters
   real(8), parameter :: pi = acos(-1.0d0)

   integer(8), intent(in)  :: nr_points
   real(8), intent(in)  :: centers(3, nr_points)
   real(8), intent(out) :: vector(nr_points)
   real(8)              :: matrix(nnbasx)
   real(8)              :: work(*)
   integer              :: lwork
   logical              :: do_matrix
   integer, intent(in)  :: ksymp

   ! Local variables
   real(8), allocatable :: ader(:, :, :)
   ! DALTON packs the density/Fock matrix.
   ! matrix_internal is an "alias" to that, stored in full format
   real(8), allocatable :: matrix_internal(:, :)
   real(8)              :: tolog, tols, factor = 1.0d0
   integer              :: ishela, ishelb, ica, icb, ia, ib
   integer              :: multa, multb, nhktab, kab
   integer              :: ipoint, idxmax, iprint, i, j
   integer              :: ierr

   iprint = 0

   allocate(matrix_internal(nbast, nbast))
   matrix_internal = 0.0d0

   ! If calculating the potential, transfer matrix to matrix_internal
   if (.not. do_matrix) then
      call dsptge(nbast, matrix, matrix_internal)
   end if

   tols = thrs**2
   tolog = 2 * log(thrs)
   idxmax = kmax
   ! Loop over bras <mu|
   idena = 0
   do ishela = 1, idxmax
     ica   = lclass(ishela)
     nhkta = nhkt(ishela)
     khkta = khkt(ishela)
     kckta = kckt(ishela)
     sphra = sphr(ishela)
     call lmnval(nhkta, kckta, lvalua, mvalua, nvalua)
     ncenta = ncent(ishela)
     icenta = nucnum(ncenta, 1)
     mula   = istbao(ishela)
     multa  = mult(mula)
     nuca   = nuco(ishela)
     numcfa = numcf(ishela)
     jsta   = jstrt(ishela)
     corax  = cent(ishela, 1, 1)
     coray  = cent(ishela, 2, 1)
     coraz  = cent(ishela, 3, 1)
     ! Loop over kets |nu>
     idenb = 0
     do ishelb = 1, ishela
       icb   = lclass(ishelb)
       ldiag = ishela .eq. ishelb
       nhktb = nhkt(ishelb)
       khktb = khkt(ishelb)
       kcktb = kckt(ishelb)
       sphrb = sphr(ishelb)
       call lmnval(nhktb, kcktb, lvalub, mvalub, nvalub)
       ncentb = ncent(ishelb)
       nhktab = nhkta + nhktb
       mulb   = istbao(ishelb)
       multb  = mult(mulb)
       nucb   = nuco(ishelb)
       numcfb = numcf(ishelb)
       jstb   = jstrt(ishelb)
       corbx  = cent(ishelb, 1, 1)
       corby  = cent(ishelb, 2, 1)
       corbz  = cent(ishelb, 3, 1)
       khktab = khkta * khktb
       kcktab = kckta * kcktb
       mab    = ior(mula, mulb)
       kab    = iand(mula ,mulb)
       hkab   = fmult(kab)

       sphrab = sphr(ishela) .or. sphr(ishelb)
       ! Calculate -<mu|1/|r-r_i||nu> integrals for each shell pair mu, nu and point r_i,
       ! then contract with the density matrix. The minus sign accounts for the
       ! charge of the electron.
       if (ica == icb) then
         allocate(ader(kckta, kcktb, nr_points), stat=ierr)
         if (ierr .ne. 0) then
             call quit("Allocating ader failed in vectorized_integration_pcm")
         end if
         ader = 0.0d0
         call vc_shell(nr_points, ader, tolog, tols, centers, iprint, work, lwork)
         if (sphrab) then
            call sphrm1(ader, ader, nr_points, work, lwork, .false., iprint)
         endif
       ! do ipoint = 1, nr_points
       !    write(lupri, *) "ELECTROSTATIC_POTENTIAL_MATRIX @point", ipoint
       !    do ib = 1, kcktb
       !      do ia = 1, kckta
       !        ! write(lupri, *) "ELEMENT", ia, ib
       !        ! write(lupri, *) ader(ia, ib, ipoint)
       !    call output(ader(:,:, ipoint),1,kckta,1,kcktb,kckta,kcktb,2,lupri)
       !      enddo
       !    enddo
       ! enddo
       do ipoint = 1, nr_points
         if (do_matrix) then
           do ib = 1, khktb
             do ia = 1, khkta
             matrix_internal(idena+ia, idenb+ib) = matrix_internal(idena+ia, idenb+ib) &
                  - ader(ia, ib, ipoint) * vector(ipoint)
             end do
           end do
         else
           vector(ipoint) = vector(ipoint) +                                  &
              sum(matrix_internal(idena+1:idena+khkta, idenb+1:idenb+khktb) * &
                        ader(1:khkta, 1:khktb, ipoint))
         end if
       end do
       deallocate(ader)
       end if
       idenb = idenb + khktb * multb
    end do
    idena = idena + khkta * multa
   end do

   if (do_matrix) then
      ! Multiply by 2.0 off-diagonal elements
      do i = 1, nbast
         do j = 1, nbast
           if (i /= j) then
              matrix_internal(i, j) = 2.0d0 * matrix_internal(i, j)
           end if
         end do
      end do
      ! Transfer from matrix_internal to matrix and clean-up
      call dgetsp(nbast, matrix_internal, matrix)
   end if

   deallocate(matrix_internal)

   !if (.not.do_matrix) then
   !   do ipoint = 1, nr_points
   !      write(lupri, *) "PRINT vector(",ipoint,")", vector(ipoint)
   !   end do
   !end if

   end subroutine vectorized_integration_pcm

   subroutine vc_shell(nr_points, ader, tolog, tols, points, iprint, work, lwork)
   !
   ! Calculates the contribution for one primitive orbital set.
   !
   ! Written, tested, debugged: R. Bast and R. Di Remigio
   !
   ! RDR 090312 Clean-up.
   !

#include "mxcent.h"
#include "maxaqn.h"
#include "aovec.h"
#include "maxorb.h"
#include "onecom.h"
#include "ader.h"
#include "primit.h"

   ! Parameters
   real(8), parameter :: pi = acos(-1.0d0)
   ! Passed variables
   integer(8), intent(in) :: nr_points
   integer, intent(in) :: iprint, lwork
   real(8) :: points(3, nr_points), ader(kckta, kcktb, nr_points)
   real(8) :: work(*)
   real(8) :: tolog, tols, factor=1.0
   ! Local variables
   real(8), allocatable :: ahgtf(:, :, :, :), odc(:, :, :, :, :, :)
   real(8), allocatable :: r(:, :, :, :, :)
   real(8) :: difab(3), corp(3), origin(3), cora(3), corb(3)
   real(8) :: difcp(3, nr_points)
   real(8) :: distab, conta, expa, contb, expb, expp, exppi
   real(8) :: expabq, saab, asaab, saab13, expapi, expbpi
   integer :: jmaxa, jmaxb, jmaxd, jmaxt, jmaxm
   integer :: ipoint, iprima, iprimb, jprima, jprimb
   integer :: idummy
   real(8) :: pval

   ! Allocation
   jmaxd = 2
   jmaxa = nhkta - 1
   jmaxb = nhktb - 1
   jmaxt = jmaxa + jmaxb + jmaxd
   jmaxm = 0
   jmax = jmaxa + jmaxb
   ! Initialization
   allocate(odc(0:jmaxa,0:jmaxb,0:jmaxt,0:jmaxd,0:jmaxm,3))
   odc = 0.0d0
   allocate(ahgtf(nr_points, 0:jmax, 0:jmax, 0:jmax))
   ahgtf = 0.0d0
   allocate(r(nr_points, 0:jmax, 0:jmax, 0:jmax, 0:jmax))
   r = 0.0d0

   cora(1) = corax
   cora(2) = coray
   cora(3) = coraz
   corb(1) = corbx
   corb(2) = corby
   corb(3) = corbz
   difab(:) = cora(:) - corb(:)
   distab = difab(1) * difab(1) + difab(2) * difab(2) + difab(3) * difab(3)

   ! Loop over primitive orbitals
   ! Shell a
   do iprima = 1, nuca
     jprima = jsta + iprima
     conta = priccf(jprima, numcfa)
     expa = priexp(jprima)
     ! Shell b
     do iprimb = 1, nucb
       jprimb = jstb + iprimb
       contb = priccf(jprimb, numcfb)
       expb = priexp(jprimb)
       expp = expa + expb
       exppi = 1.0d0 / expp
       expabq = expa * expb * exppi
       saab = conta * contb * exp(-expabq * distab)
       asaab = abs(saab)
       if (expabq * distab < tolog) then
             cycle
       end if
       if (asaab < tols) then
             cycle
       end if
       saab13 = sign(asaab**(1.0d0/3.0d0), saab)
       expapi  = expa * exppi
       expbpi  = expb * exppi
       corp(:) = expapi * cora(:) + expbpi * corb(:)
       do ipoint = 1, nr_points
         difcp(:, ipoint) =  points(:, ipoint) - corp(:)
       end do
       ! Calculate the Overlap Distribution Coefficients
       idummy = 0
       call getodc(odc, jmaxa, jmaxb, jmaxt, jmaxd, jmaxm, .false.,  &
       &      .false., onecen, expa, expb, iprint, saab13, exppi,    &
       &      work, lwork, corp(1), corp(2), corp(3),                &
       &     .true.,.false.,origin,idummy)
       call vnuc_vec(ahgtf, r, nr_points, jmax, expp, difcp)
       call cart_vc_vec(odc, jmaxa, jmaxb, jmaxt, jmaxd, jmaxm,      &
       &            ader, ahgtf, nr_points)
     end do ! Close loop over second shell
   end do ! Close loop over first shell

   deallocate(odc)
   deallocate(ahgtf)
   deallocate(r)

   end subroutine


   subroutine vnuc_vec(ahgtf, r, nr_points, jmax, pval, cp)
   !
   ! This subroutine calculates the R integrals as defined by
   ! McMurchie and Davidson in J. Comp. Phys. 26 (1978) 218.
   ! The recursion formulas (4.6) - (4.8) are used.
   !
   ! Written, tested, debugged: R. Bast, R. Di Remigio, J. Sikkema
   !
   ! JHS 260308 Only slightly slower than the implementation hernai
   !            in abacus/her1car.f of TUH
   ! RDR 090312 Clean-up.
   !
#include "maxaqn.h"
#include "gamcom.h"

   ! Parameters
   real(8), parameter :: pi = acos(-1.0d0)
   ! Passed variables
   integer :: jmax
   integer(8) :: nr_points
   real(8), intent(out) :: ahgtf(nr_points, 0:jmax, 0:jmax, 0:jmax)
   real(8), intent(in)  :: cp(3, nr_points)
   real(8) :: r(nr_points, 0:jmax, 0:jmax, 0:jmax, 0:jmax)
   real(8) :: pval
   ! Local variables
   real(8) :: factor, prod
   integer :: jval, t, u, v, ipoint

   do ipoint = 1, nr_points
     ! Incomplete gamma function
     wval = pval * (cp(1,ipoint)**2 + cp(2,ipoint)**2 + cp(3,ipoint)**2)
     jmax0 = jmax
     call gamma_function
     ! Calculate r(ipoint, 0, 0, 0, jval)
     factor = (2.0d0 * pi) / pval
     do jval = 0, jmax
       fjw(jval)         =   factor * fjw(jval)
       factor            = - 2.0d0 * pval * factor
       r(ipoint, 0, 0, 0, jval)  =   fjw(jval)
     end do
     ! Calculate r(t, u, v, jval)
     do jval = jmax, 1, -1
       do v = 0, jmax - jval
         prod = -cp(3, ipoint) * r(ipoint, 0, 0, v, jval)
         if (v > 0) then
           prod = prod + v * r(ipoint, 0, 0, v - 1, jval)
         end if
         r(ipoint, 0, 0, v + 1, jval - 1) = prod
         do u = 0, jmax - jval - v
           prod = -cp(2, ipoint) * r(ipoint, 0, u, v, jval)
           if (u > 0) then
             prod = prod + u * r(ipoint, 0, u - 1, v, jval)
           end if
           r(ipoint, 0, u + 1, v, jval - 1) = prod
           do t = 0, jmax - jval - u - v
             prod = -cp(1, ipoint) * r(ipoint, t, u, v, jval)
             if (t > 0) then
               prod = prod + t * r(ipoint, t - 1, u, v, jval)
             end if
             r(ipoint, t + 1, u, v, jval - 1) = prod
           end do
         end do
       end do
     end do
   end do
   ! The nuclear attraction integrals are given as r(ipoint, t, u, v, 0)
   do v = 0, jmax
     do u = 0, jmax
       do t = 0, jmax
         do ipoint = 1, nr_points
           ahgtf(ipoint, t, u, v) = - r(ipoint, t, u, v, 0)
         end do
       end do
     end do
   end do

   end subroutine

   subroutine cart_vc_vec(odc, jmaxa, jmaxb, jmaxt, jmaxd, jmaxm, ader, ahgtf, nr_points)
   !
   ! Written, tested, debugged: R. Bast, R. Di Remigio, J. Sikkema
   !
   ! RDR 060312 Clean-up.
   !

#include "maxaqn.h"
#include "onecom.h"
#include "lmns.h"

   ! Passed variables
   integer :: jmaxa, jmaxb, jmaxt, jmaxd, jmaxm
   integer(8) :: nr_points
   real(8) :: ahgtf(nr_points, 0:jmax, 0:jmax, 0:jmax)
   real(8) :: ader(kckta, kcktb, nr_points)
   real(8) :: odc(0:jmaxa, 0:jmaxb, 0:jmaxt, 0:jmaxd, 0:jmaxm, 3)
   ! Local variables
   real(8) :: ev, ee, eee
   integer :: icompa, lvala, mvala, nvala
   integer :: icompb, lvalb, mvalb, nvalb
   integer :: t, u, v, ipoint

   do icompa = 1,kckta
       lvala = lvalua(icompa)
       mvala = mvalua(icompa)
       nvala = nvalua(icompa)
     do icompb = 1,kcktb
         lvalb = lvalub(icompb)
         mvalb = mvalub(icompb)
         nvalb = nvalub(icompb)
       do v = 0, nvala + nvalb
           ev = odc(nvala,nvalb,v,0,0,3)
         do u = 0, mvala + mvalb
             ee = odc(mvala,mvalb,u,0,0,2)*ev
           do t = 0, lvala + lvalb
               eee = odc(lvala,lvalb,t,0,0,1)*ee
             do ipoint = 1, nr_points
               ader(icompa, icompb, ipoint) = ader(icompa, icompb, ipoint) + eee * ahgtf(ipoint, t, u, v)
             end do
           end do
         end do
       end do
     end do
   end do

   end subroutine

   subroutine gamma_function
   !
   ! Trygve Ulf Helgaker fall 1984
   !
   ! This subroutine calculates the incomplete gamma function as
   ! described by McMurchie & Davidson, J. Comp. Phys. 26 (1978) 218.
   !
   ! Roberto Di Remigio May 2012
   ! Purified from the evil implicit.h and all the other common blocks.
   !
#include "maxaqn.h"

   real(8), parameter :: d1 = 1.0d0,  d2 = 2.0d0, d10 = 10.0d0
   real(8), parameter :: half = 0.5d0, tenth = 0.1d0, ten6 = 1.0d6
   real(8), parameter :: pi = acos(-1.0d0)
   real(8), parameter :: sqrtpi = sqrt(pi)
   real(8), parameter :: pi2 = pi * pi
   real(8), parameter :: sqrpih = sqrtpi/d2
   real(8), parameter :: coef2 = half,  coef3 = - d1/6.0d0, coef4 = d1/24.0d0
   real(8), parameter :: coef5 = - d1/120.0d0, coef6 = d1/720.0d0
   real(8), parameter :: gfac30 = 0.4999489092d0, gfac31 = -0.2473631686d0,        &
   &   gfac32 = 0.321180909d0, gfac33 = -0.3811559346d0, gfac20 = 0.4998436875d0,  &
   &   gfac21 = -0.24249438d0, gfac22 = 0.24642845d0, gfac10 = 0.499093162d0,      &
   &   gfac11 = -0.2152832d0, gfac00 = 0.490d0


   ! Local variables
   real(8)            :: tabjfw, wdif, d2wal, rexpw, denom, rwval, summ, term
   real(8)            :: r2max1, d2max1, gval, factor
   integer            :: jmax, jmx, j, iadr, jadr, maxj0, istart, ipoint, iorder

#include "gamcom.h"
   !
   save maxj0
   data maxj0 /-1/
   !
   ipoint = d10 * min(wval, ten6) + half
   !     have seen problems with NINT intrinsic function here (rarely)
   !     therefore the "+ HALF" before integer truncation
   if (ipoint < 0) then
      call quit('Fatal error in gammafun')
   else if (ipoint < 120) then
      istart = 1 + 121 * jmax0 + ipoint
      wdif = wval - tenth * ipoint
      fjw(jmax0) = (((((coef6 * tabfjw(istart + 726) * wdif    &  ! 726 = 6*121
   &                   + coef5 * tabfjw(istart + 605)) * wdif   &
   &                    + coef4 * tabfjw(istart + 484)) * wdif  &
   &                     + coef3 * tabfjw(istart + 363)) * wdif   &
   &                      + coef2 * tabfjw(istart + 242)) * wdif  &
   &                       - tabfjw(istart + 121)) * wdif       &
   &                        + tabfjw(istart)
      d2wal = d2 * wval
      rexpw = exp(-wval)
      denom = 2.0d0 * jmax0 + 1.0d0
      do j = jmax0, 1, -1
        denom = denom - d2
        fjw(j - 1) = (d2wal * fjw(j) + rexpw) / denom
      end do
   else if (ipoint <= (20 * jmax0 + 360)) then
      rwval = d1 / wval
      rexpw = exp(-wval)
      gval = gfac30 + rwval * (gfac31 + rwval * (gfac32 + rwval * gfac33))
      fjw(0) = sqrpih * sqrt(rwval) - rexpw * gval * rwval
      factor = half * rwval
      term = factor * rexpw
      do j = 1, jmax0
        fjw(j) = factor * fjw(j - 1) - term
        factor = rwval + factor
      end do
   else
      rwval  = d1 / wval
      fjw(0) = sqrpih * sqrt(rwval)
      factor = half * rwval
      do j = 1, jmax0
        fjw(j) = factor * fjw(j - 1)
        factor = rwval + factor
      end do
   end if
   return
   !
   !     ***** Tabulation of incomplete gamma function *****
   !
   entry gamtab(jmx)
   !
   !     For j = jmx a power series expansion is used, see for
   !     example Eq.(39) given by V. Saunders in "Computational
   !     Techniques in Quantum Chemistry and Molecular Physics",
   !     Reidel 1975.  For j < jmx the values are calculated
   !     using downward recursion in j.
   !
   !
   if (jmx > maxj) then
      call quit('Gamtab error: jmx greater than limit.')
   end if
   jmax = jmx + 6
   maxj0 = jmax
   !
   !     WVAL = 0.0
   !
   iadr = 1
   denom = d1
   do j = 0, jmax
     tabfjw(iadr) = d1 / denom
     iadr = iadr + 121
     denom = denom + d2
   end do
   !
   !     WVAL = 0.1, 0.2, 0.3,... 12.0
   !
   iadr = iadr - 121
   d2max1 = 2.0d0 * jmax + 1.0d0
   r2max1 = d1 / d2max1
   do ipoint = 1, 120
      wval = tenth * ipoint
      d2wal = wval + wval
      iadr = iadr + 1
      term = r2max1
      summ = term
      denom = d2max1
      do iorder = 2, 200
        denom = denom + d2
        term = term * d2wal / denom
        summ = summ + term
        if (term .le. 1.0d-15) exit
      end do
      rexpw = exp(-wval)
      tabfjw(iadr) = rexpw * summ
      denom = d2max1
      jadr = iadr
      do j = 1, jmax
         denom = denom - d2
         tabfjw(jadr - 121) = (tabfjw(jadr) * d2wal + rexpw) / denom
         jadr = jadr - 121
      end do
   end do

   end subroutine

   subroutine pot_int_tess(potint, tessera, trimat, work, lwork)

#include "dummy.h"
#include "maxorb.h"
#include "mxcent.h"
#include "priunit.h"
#include "orgcom.h"
#include "inforb.h"

   real(8), intent(out) :: potint(*)
   real(8), intent(in)  :: tessera(3)
   real(8)              :: work(*)
   logical              :: trimat, tofile, exp1vl
   character(7)         :: intlab
   character(8)         :: labint(9*mxcent)
   integer              :: intrep(9*mxcent), intadr(9*mxcent), lwork
   integer              :: il, iprint, j, kfree, kpatom, lfree, ncomp
   real(8)              :: xdiporg, ydiporg, zdiporg

   xdiporg = diporg(1)
   ydiporg = diporg(2)
   zdiporg = diporg(3)
   diporg(1) = tessera(1)
   diporg(2) = tessera(2)
   diporg(3) = tessera(3)

   intlab = 'NPETES '
   ncomp = nsym
   il = 1
   tofile = .false.
   kpatom = 0
   exp1vl = .false.
   iprint = 0

   kfree = 1
   lfree = lwork - kfree + 1
   if (lfree .lt. 0) call errwrk('pot_int_tess', kfree, lwork)

   call get1in(potint,intlab,ncomp,work(kfree),lfree,labint, &
               intrep,intadr,il,tofile,kpatom,trimat,dummy,  &
               exp1vl,dummy,iprint)

   write(lupri, *) "Potentials for tessera rsp",1,(diporg(j), j=1,3),mxcent
   call outpak(potint, norbt, 1, lupri)

   diporg(1) = xdiporg
   diporg(2) = ydiporg
   diporg(3) = zdiporg

   end subroutine pot_int_tess

end module pcm_integrals
