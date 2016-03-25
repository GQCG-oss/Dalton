module pcm_linrsp

  use, intrinsic :: iso_c_binding
  use pcmsolver
  use pcm_config
  use pcm_write
  use pcm_integrals

  implicit none

  public pcm_lr_driver

  private

contains

  subroutine pcm_lr_driver(nosim, bovecs, cmo, indxci, udv, dv, udvtr, dvtr, dtv, dtvtr, sovecs, wrk, lwrk)
    ! Shamelessly copied from the PCMLTR subroutine in rsp/rspief.F
    ! Originally written by Luca Frediani, Benedetta Mennucci and Roberto Cammi (24.10.01)
    ! RDR (10.11.14) remove everything configurational from this subroutine

#include "dummy.h"
#include "inforb.h"
#include "infrsp.h"
#include "inftap.h"
#include "wrkrsp.h"
#include "infvar.h"

    integer(4) :: nosim, indxci(*), lwrk
    real(8) :: bovecs(*), cmo(*)
    real(8) :: udv(*), dv(*)
    real(8) :: udvtr(*), dvtr(*)
    real(8) :: dtv(*), dtvtr(*)
    real(8) :: sovecs(*)
    real(8) :: wrk(lwrk)

    if ( nosim .gt. 0 ) then
      if (iprrsp.gt.101) then
        write(pcm_global%print_unit,*)' LINEAR TRANSFORMED ORBITAL VECTOR'
        write(pcm_global%print_unit,*)' **** BEFORE orbital_hessian **** ',kzyvar,nosim
        call output(sovecs,1,kzyvar,1,nosim,kzyvar,nosim,1,pcm_global%print_unit)
      end if

      call orbital_hessian(pcm_global%nr_points, pcm_global%nr_points_irr, pcm_global%tess_cent,  &
                           nosim, bovecs, cmo, indxci, udv, dv, udvtr, dvtr, sovecs, wrk, lwrk)

      if (iprrsp.gt.101) then
        write(pcm_global%print_unit,*)' LINEAR TRANSFORMED ORBITAL VECTOR'
        write(pcm_global%print_unit,*)' **** AFTER orbital_hessian **** '
        call output(sovecs,1,kzyvar,1,nosim,kzyvar,nosim,1,pcm_global%print_unit)
      end if
    end if

  end subroutine pcm_lr_driver

  subroutine orbital_hessian(nr_points, nr_points_irr, centers, nosim, bovecs, cmo, &
                             indxci, udv, dv, udvtr, dvtr, sovecs, wrk, lwrk)
    !  Purpose:  Calculate MCSCF E2 contribution from a
    !            surrounding ief-pcm medium to an orbital trial vector.

    real(8), external :: slvtlm, slvqlm

    integer(8), intent(in) :: nr_points, nr_points_irr
    real(8),    intent(in) :: centers(3, nr_points)
    integer(4) :: nosim, indxci(*), lwrk
    real(8) :: bovecs(*), cmo(*)
    real(8) :: udv(*), dv(*)
    real(8) :: udvtr(*), dvtr(*)
    real(8) :: wrk(lwrk)
#include "iratdef.h"
    !
    real(8), parameter :: d0 = 0.0d0, d1 = 1.0d0, d2 = 2.0d0, dp5 = 0.5d0
    logical :: fndlab, tofile, exp1vl, trimat
#include "dummy.h"
    !
    !  Used from common blocks:
    !    INFINP : NLMSOL, LSOLMX, INERSI
    !    INFORB : NNASHX, NNORBX, NNBASX, etc.
    !    INFVAR : JWOP
    !    INFRSP :
    !    WRKRSP :
    !    INFTAP : LUSOL,  LBSYMB
    !
#include "maxash.h"
#include "maxorb.h"
#include "mxcent.h"
    integer(4) :: intrep(9*mxcent), intadr(9*mxcent)
    character*8 :: labint(9*mxcent)
    logical :: file_exist
#include "orgcom.h"
#include "infinp.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "inftap.h"
#include "infpri.h"
#include "infpar.h"
    real(8) :: sovecs(kzyvar, *)

    ! Local variables
    integer(4) :: i, iosim, its, jaop, jj1ao, jj1x, jj1xac, jpcmx
    integer(4) :: jubo, kintad, kintrp, kj1, kj1ao, kj1sq, kj1xac, kj1xact
    integer(4) :: kj1xsq, kjpxac, kjpxact, kneed, kovlp, kpatom
    integer(4) :: kpcmx, kpcmxa, kpcmxat, kpcmxt, kqtex, ktssim
    integer(4) :: ktssymop, kubo, kucmo, kvtex, kw10, kw20, kw30, kw50
    integer(4) :: l, lw50, ncomp
    logical    :: laddmm = .true.
    real(c_double), allocatable :: oit_mep(:, :), oit_asc(:, :)
    real(c_double), allocatable :: oit_mep_slice(:), oit_asc_slice(:)
    real(c_double), allocatable :: shuffled_oit_asc(:, :)
    real(c_double), allocatable :: asc(:)
    real(8) :: factor, qfactor, qtexs, tj1xac, xi, yi, zi
    integer(c_int) :: irrep
    integer :: jts

    character(7) :: totasc
    character(7) :: oitmep, oitasc

    !
    !     Core allocation
    !
    kintrp = 1
    kintad = kintrp + (3*mxcoor + 1)/irat
    kvtex  = kintad + (3*mxcoor + 1)/irat
    kqtex  = kvtex + nr_points*nosim
    kubo   = kqtex + nr_points*nosim
    kw10   = kubo  + nosim*n2orbx
    !
    kucmo  = kw10
    kpcmx  = kucmo  + norbt*nbast

    if (trplet) then
      kpcmxt  = kpcmx  + nosim*n2orbx
      kpcmxa  = kpcmxt + nosim*n2orbx
      kpcmxat = kpcmxa + nosim*n2ashx
      kjpxac  = kpcmxat + nosim*n2ashx
    else
      kpcmxa = kpcmx  + nosim*n2orbx
      kjpxac = kpcmxa + nosim*n2ashx
    endif
    if (trplet) then
      kjpxact = kjpxac + nosim
      kw20    = kjpxact + nosim
    else
      kw20   = kjpxac + nosim
    endif
    !
    kj1ao  = kw20
    kj1sq  = kj1ao  + nnbasx*nsym
    kj1    = kj1sq  + n2orbx
    kj1xsq = kj1    + nnorbx
    kj1xac = kj1xsq + n2orbx*nosim
    if (trplet) then
      kj1xact = kj1xac   + nosim * n2ashx
      kw30    = kj1xact  + nosim * n2ashx
    else
      kw30    = kj1xac   + nosim * n2ashx
    endif
    !
    !  3.0 SOLSC
    kovlp = kw30
    kw50   = kovlp  + nosim
    lw50   = lwrk  + 1 - kw50
    !
    kneed = max(kw30, kw50)
    if (kneed .gt. lwrk) call errwrk('orbital_hessian', kneed, lwrk)

    ! ksymop is the index of the irrep the perturbation operator belongs to
    ! the indexing is 1-based while the PCM module uses a 0-based indexing
    ! Transfer the information in ksymop to irrep
    irrep = ksymop - 1

    !  Unpack symmetry blocked CMO
    call upkcmo(cmo, wrk(kucmo))

    !  Calculate unpacked orbital trial vectors in UBO
    if (nosim .gt. 0) then
      call rspzym(nosim, bovecs, wrk(kubo))
      call dscal(nosim * n2orbx, -1.0d0, wrk(kubo), 1)
      if (iprrsp .ge. 55) then
        do iosim = 1, nosim
          jubo = kubo + (iosim-1)*n2orbx
          write(pcm_global%print_unit,*) 'Orbital trial vector unpacked to matrix form (no.', IOSIM, ' of', NOSIM
          call output(wrk(jubo), 1, norbt, 1, norbt, norbt, norbt, 1, pcm_global%print_unit)
        end do
      end if
    end if

    ! Contributions from all tesserae are included.
    call dzero(wrk(kpcmx), nosim * n2orbx)
    if (trplet) call dzero(wrk(kpcmxt), nosim * n2orbx)
    xi = diporg(1)
    yi = diporg(2)
    zi = diporg(3)

    ! Set names for surface functions
    totasc = 'TotASC'//char(0)
    oitmep = 'oitMEP'//char(0)
    oitasc = 'oitASC'//char(0)
    ! Retrieve total unperturbed surface charges from module
    allocate(asc(nr_points))
    asc = 0.0d0
    call pcmsolver_get_surface_function(context_, nr_points, asc, totasc)

    allocate(oit_mep(nr_points, nosim))
    oit_mep = 0.0d0
    allocate(oit_asc(nr_points, nosim))
    oit_asc = 0.0d0

    NrPointsIrr: do i = 1, nr_points_irr
      ! Read AO potential integrals on tesserae
      l = 1
      ncomp     = nsym
      diporg(1) = centers(1, i)
      diporg(2) = centers(2, i)
      diporg(3) = centers(3, i)
      exp1vl    = .false.
      tofile    = .false.
      kpatom    = 0
      trimat    = .true.
      ! Calculate V_{mu, nu}^i (AO basis)
      call get1in(wrk(kj1ao), 'NPETES ', ncomp, wrk(kw50), lw50,        &
        labint, wrk(kintrp), wrk(kintad), l, tofile ,kpatom,  &
        trimat, dummy, exp1vl, dummy, iprrsp)
      jj1ao = kj1ao
      ! Transform to MO basis: V_{pq}^i
      call uthu(wrk(jj1ao), wrk(kj1), wrk(kucmo), wrk(kw30), nbast, norbt)
      ! Transform V_{pq}^i from triangular to square format
      call dsptsi(norbt, wrk(kj1), wrk(kj1sq))
      call dzero(wrk(kj1xsq), n2orbx * nosim)
      do iosim = 1, nosim
        jubo = kubo + (iosim - 1) * n2orbx
        jj1x = kj1xsq + (iosim - 1) * n2orbx
        jj1xac = kj1xac + (iosim - 1) * n2ashx
        call onexh1(wrk(jubo), wrk(kj1sq), wrk(jj1x))
        if (nasht .gt. 0) then
          call getacq(wrk(jj1x), wrk(jj1xac))
        end if
        if (iprrsp .ge. 15) then
          write(pcm_global%print_unit,'(/A,I5)') 'J1X_mo matrix: tess', i
          call output(wrk(jj1x), 1, norbt, 1, norbt, norbt, norbt, 1, pcm_global%print_unit)
          if (nasht .gt. 0) then
            write(pcm_global%print_unit,'(/A)') ' J1X_ac matrix:'
            call output(wrk(jj1xac), 1, nasht, 1, nasht, nasht, nasht, 1, pcm_global%print_unit)
          end if
        end if
        !     Expectation value of transformed potential on tesserae:
        !                     <0|\tilde{V}|0>
        if (irrep .eq. 0) then
          if (trplet) then
            factor = slvtlm(udvtr, wrk(jj1xac), wrk(jj1x), tj1xac)
          else
            factor = slvqlm(udv, wrk(jj1xac), wrk(jj1x), tj1xac)
          end if
          oit_mep(i, iosim) = factor
          if (iprrsp .ge. 6) then
            write(pcm_global%print_unit, *) ' oit_mep(',i,', ', iosim,') = ', oit_mep(i, iosim)
            write(pcm_global%print_unit,'(A,F17.8)') ' --- active part of J1X    :',TJ1XAC
          end if
        end if
      end do
      qfactor = asc(i)
      if (trplet) then
        call daxpy(nosim * n2orbx, qfactor, wrk(kj1xsq), 1, wrk(kpcmxt), 1)
      else
        call daxpy(nosim * n2orbx, qfactor, wrk(kj1xsq), 1, wrk(kpcmx), 1)
      endif
      ! KPCMX: \tilde{J} + \tilde{X}
      !
      !     For non-totally symmetric perturbation operators
      !
      if (irrep .ne. 0) then
        its = irrep * nr_points_irr + i
        ! Transform AO pot. int. into MO basis  V(AO) --> V(MO)
        jj1ao = kj1ao + (ksymop - 1) * nnbasx
        call uthu(wrk(jj1ao), wrk(kj1), wrk(kucmo), wrk(kw50), nbast, norbt)
        ! Transform V(MO) from triangular to square format
        call dsptsi(norbt, wrk(kj1), wrk(kj1sq))
        call dzero(wrk(kj1xsq), n2orbx * nosim)
        do iosim = 1, nosim
          jubo = kubo + (iosim - 1) * n2orbx
          jj1x = kj1xsq + (iosim - 1) * n2orbx
          jj1xac = kj1xac + (iosim - 1) * n2ashx
          call onexh1(wrk(jubo), wrk(kj1sq), wrk(jj1x))
          if (nasht .gt. 0) then
            call getacq(wrk(jj1x), wrk(jj1xac))
          end if
          if (iprrsp .ge. 15) then
            write(pcm_global%print_unit,'(/A)') ' J1X_mo matrix :'
            call output(wrk(jj1x), 1, norbt, 1, norbt, norbt, norbt, 1, pcm_global%print_unit)
            if (nasht .gt. 0) then
              write(pcm_global%print_unit,'(/A)') ' J1X_ac matrix :'
              call output(wrk(jj1xac), 1, nasht, 1, nasht, nasht, nasht, 1, pcm_global%print_unit)
            end if
          end if
          !     Expectation value of transformed potential on tesserae:
          !                     <0|\tilde{V}|0>
          if (trplet) then
            factor = slvtlm(udvtr, wrk(jj1xac), wrk(jj1x), tj1xac)
          else
            factor = slvqlm(udv, wrk(jj1xac), wrk(jj1x), tj1xac)
          end if
          oit_mep(its, iosim) = factor
          if (iprrsp .ge. 6) then
            write(pcm_global%print_unit,'(A,F17.8)') ' oit_mep(',i,',',iosim,') = ', oit_mep(i, iosim)
            write(pcm_global%print_unit,'(A,F17.8)') ' --- active part of J1X    :',TJ1XAC
          end if
        end do
      end if
    end do NrPointsIrr

    if (iprrsp .ge. 50) then
      do iosim = 1, nosim
        jpcmx = kpcmx + (iosim - 1) * n2orbx
        write(pcm_global%print_unit,'(/A,I3,A,I3)') ' >>> orbital_hessian - (JPCMX) half.',IOSIM,' of',NOSIM
        call output(wrk(jpcmx), 1, norbt, 1, norbt, norbt, norbt, 1, pcm_global%print_unit)
      end do
    end if
    allocate(oit_mep_slice(nr_points))
    allocate(oit_asc_slice(nr_points))
    do iosim = 1, nosim
      ! Set a cavity surface function with the MEP
      oit_mep_slice = 0.0d0
      oit_mep_slice = oit_mep(:, iosim)
      oit_asc_slice = 0.0d0
      call pcmsolver_set_surface_function(context_, nr_points, oit_mep_slice, oitmep)
      call pcmsolver_compute_response_asc(context_, oitmep, oitasc, irrep)
      ! Get polarization charges @tesserae centers
      call pcmsolver_get_surface_function(context_, nr_points, oit_asc_slice, oitasc)
      oit_asc(:, iosim) = oit_asc_slice
      if (iprrsp .ge. 6) then
        do i = 1, nr_points
          write(pcm_global%print_unit, *) ' oit_asc(',i,', ', iosim,') = ', oit_asc(i, iosim)
        end do
      end if
    end do
    deallocate(oit_mep_slice)
    deallocate(oit_asc_slice)
    ! TRANSFORMED CHARGES MULTIPIED TO POTENTIALS.
    ! There are only one-index-transformed charges for the totally symmetric irrep.
    call dzero(wrk(kj1xsq), nosim * nnbasx)
    ! The one-index transformed ASC is stored in oit_asc according to its irrep.
    ! Move it to the top segment (associated with irrep 0 <--> totally symmetric)
    ! if we are handling the non totally-symmetric irrep of the point group
    if (irrep .gt. 0) then
      allocate(shuffled_oit_asc(nr_points, nosim))
      shuffled_oit_asc = 0.0d0
      do i = 1, nr_points_irr
        shuffled_oit_asc(i, :) = oit_asc(irrep * nr_points_irr + i, :)
      end do
      call j1int_pcm(shuffled_oit_asc, nr_points, nr_points_irr, pcm_global%tess_cent, .false., wrk(kj1xsq), &
        nosim, .false., 'NPETES ', ksymop, wrk(kw50), lw50)
    else
      call j1int_pcm(oit_asc, nr_points, nr_points_irr, pcm_global%tess_cent, .false., wrk(kj1xsq),          &
        nosim, .false., 'NPETES ', ksymop, wrk(kw50), lw50)
    end if
    do iosim = 1, nosim
      jaop  = kj1xsq + (iosim - 1) * nnbasx
      jpcmx = kpcmx  + (iosim - 1) * n2orbx
      call uthu(wrk(jaop), wrk(kj1), wrk(kucmo), wrk(kw50), nbast, norbt)
      call dsptsi(norbt, wrk(kj1), wrk(kj1sq))
      call daxpy(n2orbx, -d1, wrk(kj1sq), 1, wrk(jpcmx), 1)
    end do
    ! Restore dipole origin
    diporg(1) = xi
    diporg(2) = yi
    diporg(3) = zi

    if (iprrsp .ge. 50) then
      WRITE(pcm_global%print_unit,*) ' >>> orbital_hessian - (KJ1XSQ)'
      call output(wrk(kj1xsq),1,nbast,1,nbast,nbast,nbast,1,pcm_global%print_unit)
      do iosim = 1,nosim
        jpcmx = kpcmx + (iosim-1)*n2orbx
        WRITE(pcm_global%print_unit,'(/A,I3,A,I3)') ' >>> IEFLNO - (JPCMX) matrix no.',IOSIM,' of',NOSIM
        call output(wrk(jpcmx),1,norbt,1,norbt,norbt,norbt,1,pcm_global%print_unit)
      end do
    end if
    ! \sigma_{oo} = <0|[q_j,\tilde{Z}]|0>
    laddmm = .true. ! Related to MM/PCM
    if (.not. laddmm) then
      call dzero(wrk(kpcmx),nosim*n2orbx)
      if (trplet) call dzero(wrk(kpcmxt),nosim*n2orbx)
    endif

    if(trplet) then
      call slvsor(.true.,.false.,nosim,udvtr,sovecs(1,1),wrk(kpcmx))
      call slvsor(.true.,.true., nosim,udv,  sovecs(1,1),wrk(kpcmxt))
    else
      call slvsor(.true.,.true., nosim,udv,  sovecs(1,1),wrk(kpcmx))
    endif

  end subroutine orbital_hessian

end module pcm_linrsp
