      subroutine pe_grad(cref, cmo, cindx, dv, grd, energy, wrk, nwrk)
!
!     Written by Erik Donovan Hedegård (edh) and Jógvan Magnus H. Olsen
!                based on PCMGRAD 
!
!     Purpose:  calculate (MCSCF) energy and gradient contribution
!               from a PE potential using the PE module
!
!     Output:
!     grd       MCSCF gradient with PE contribution added
!     energy    total PE energy
!
! Used from common blocks:
!   INFVAR: NCONF,  NWOPT,  NVAR,   NVARH
!   INFORB: NNASHX, NNBASX, NNORBX, etc.
!   INFTAP: LUIT2
!

      use polarizable_embedding, only: pe_master

      implicit none

#include "priunit.h"
#include "infvar.h"
#include "inforb.h"
#include "inftap.h"

      integer :: nwrk
      real*8 :: energy
      real*8, dimension(*) :: cref, cmo, cindx, dv, grd
      real*8, dimension(nwrk) :: wrk
      character*8 :: star8 = '********'
      character*8 :: solvdi = 'SOLVDIAG'
      character*8 :: eodata = 'EODATA  '

      logical :: fndlab
      integer :: nc4, nw4, i
      real*8 :: solelm, ddot
      real*8 :: tmo, tac, test
      real*8, dimension(:), allocatable :: fckmo, fckac
      real*8, dimension(:), allocatable :: pegrd, diape

      call qenter('pe_grad')

      allocate(fckmo(nnorbx))
      call pe_fckmo(cmo, fckmo, energy, dv, wrk, nwrk)

      allocate(fckac(nnashx))
      if (nasht > 0) call getac2(fckmo, fckac)

      tmo = solelm(dv, fckac, fckmo, tac)

      allocate(pegrd(nvarh))
      pegrd = 0.0d0
      if (nconf > 1) then
         ! edh: SOLGC calc. < u | Fg | 0 > + < 0 | Fg | 0 > c_u
        call solgc(cref, fckac, tac, pegrd, cindx, wrk, nwrk)
      end if
      if (nwopt > 0) then
        ! edh: SOLGO calc. 2 < 0 | [Ers, Fg] | 0 >
        call solgo(2.0d0, dv, fckmo, pegrd(1+nconf:nvarh))
      end if

      allocate(diape(nvar))
      diape = 0.0d0
      call soldia(tac, fckac, cindx, fckmo, dv, diape, wrk, nwrk)
      diape = - diape
      deallocate(fckmo, fckac)

! ******************* Orthogonality test **********************
       test = ddot(nconf, cref, 1, pegrd, 1)
       if (abs(test) > 1.0d-8) then
          nwarn = nwarn + 1
          write(lupri,*) ' >>> PE GRADIENT WARNING <<< '
          write(lupri,*) ' < CREF | GRAD > =', test
       end if
! *************************************************************
     
!     Add PE gradient contribution to MCSCF gradient
      call daxpy(nvarh, 1.0d0, pegrd, 1, grd, 1)
      deallocate(pegrd)

      if (luit2 > 0) then
         nc4 = max(nconf, 4)
         nw4 = max(nwopt, 4)
         rewind luit2
         if (fndlab(eodata,luit2)) backspace luit2
         write(luit2) star8, star8, star8, solvdi
         if (nconf > 1) call writt(luit2, nc4, diape)
         write(luit2) star8, star8, star8, eodata
      end if

      call qexit('pe_grad')

      end subroutine pe_grad

!------------------------------------------------------------------------------

      subroutine pe_fckmo(cmo, fckmo, energy, dv, wrk, nwrk)
!
!     Written by Erik Donovan Hedegård (edh) and Jógvan Magnus H. Olsen
!
!     Purpose:  Get PE operator in MO basis using PE module
!
!     Output:
!     fckmo     PE operator in MO basis
!     energy    PE energy
!
! Used from common blocks:
!   INFVAR: NCONF,  NWOPT,  NVAR,   NVARH
!   INFORB: NNASHX, NNBASX, NNORBX, etc.
!   INFTAP: LUIT2
!

      use polarizable_embedding, only: pe_master

      implicit none

#include "priunit.h"
#include "infvar.h"
#include "inforb.h"
#include "inftap.h"

      integer :: nwrk
      real*8, dimension(1) :: energy
      real*8, dimension(*) :: fckmo
      real*8, dimension(*) :: cmo, dv
      real*8, dimension(nwrk) :: wrk

      integer :: i
      real*8, dimension(:), allocatable :: dcao, dvao, fdtao, fckao

      call qenter('pe_fckmo')

      allocate(dcao(n2basx), dvao(n2basx))
      call fckden((nisht>0), (nasht>0), dcao, dvao, cmo, dv, wrk, nwrk)
      if (nasht > 0) dcao = dcao + dvao
      deallocate(dvao)

      allocate(fdtao(nnbasx))
      call dgefsp(nbast, dcao, fdtao)
      deallocate(dcao)

      allocate(fckao(nnbasx))
      call pe_master(runtype='fock', denmats=fdtao, fckmats=fckao,
     &               nmats=1, energies=energy, dalwrk=wrk)
      deallocate(fdtao)

      call uthu(fckao, fckmo, cmo, wrk, nbast, norbt)
      deallocate(fckao)

      call qexit('pe_fckmo')

      end subroutine pe_fckmo

!------------------------------------------------------------------------------

      subroutine pe_lin(ncsim, nosim, bcvecs, bovecs, cref, cmo, cindx,
     &                  dv, dtv, scvecs, sovecs, orblin, wrk, nwrk)
!
! Written by Erik Donovan Hedegård and Jógvan Magnus H. Olsen
!            after original code by  Hans Joergen Aa. Jensen
!
! Common driver for pe_lnc and pe_lno
!
!   Used from common blocks:
!   INFLIN : NWOPPT,NVARPT

      implicit none

#include "priunit.h"
#include "inflin.h"
#include "infvar.h"
#include "inforb.h"

      logical :: orblin
      integer :: ncsim, nosim, nwrk
      real*8, dimension(*) :: bcvecs, bovecs, scvecs, sovecs
      real*8, dimension(*) :: cmo, cref, cindx, dv, dtv
      real*8, dimension(nwrk) :: wrk

      integer :: i, nso

      call qenter('pe_lin')

      if (ncsim > 0) then
        call pe_lnc(ncsim, bcvecs, cref, cmo, cindx, dv, dtv, scvecs,
     &              wrk, nwrk)
      end if

      if (nosim > 0) then
        if (.not. orblin) then
            nso = nvarpt
        else
            nso = nwoppt
        end if
        call pe_lno(nosim, bovecs, cref, cmo, cindx, dv, sovecs, nso,
     &              wrk, nwrk)
      end if

      call qexit('pe_lin')

      end subroutine pe_lin

!------------------------------------------------------------------------------

      subroutine pe_lnc(ncsim, bcvecs, cref, cmo, cindx, dv, dtv,
     &                  scvecs, wrk, nwrk)
!
!  Written by Erik Donovan Hedegaard and Jógvan Magnus H. Olsen
!             after original routine by Hans Jørgen Aa. Jensen
!
!  Purpose:  Calculate Hessian contribution from a polarizable
!            embedding potantial to a csf trial vector.
!
!
!  Used from common blocks:
!    INFORB : NNASHX, NNORBX, NNBASX, etc.
!    INFVAR : NWOPH
!    INFLIN : NCONST, NVARPT, NWOPPT
!

      use polarizable_embedding, only: pe_master, pe_polar

      implicit none

#include "priunit.h"
#include "dummy.h"
#include "inforb.h"
#include "infvar.h"
#include "inflin.h"
#include "infdim.h"

      integer :: ncsim, nwrk
      real*8, dimension(*) :: bcvecs, cref, cmo, cindx, dv
      real*8, dimension(nnashx,*) :: dtv
      real*8, dimension(nvarpt,*) :: scvecs
      real*8, dimension(nwrk) :: wrk

      logical :: fndlab
      integer :: i, j, jscvec, mwoph
      real*8 :: tfxc, tfyc, tfycac, energy, solelm
      real*8, dimension(:), allocatable :: udtvao, fdtvaos, fxcaos
      real*8, dimension(:), allocatable :: tfxcacs, fyc, fycac
      real*8, dimension(:,:), allocatable :: fxcs, fxcacs

      call qenter('pe_lnc')

      allocate(fxcs(nnorbx,ncsim))
      allocate(fxcacs(nnashx,ncsim))
      allocate(tfxcacs(ncsim))
      if (pe_polar) then
        allocate(udtvao(n2basx))
        allocate(fdtvaos(ncsim*nnbasx))
        do i = 1, ncsim
            j = (i - 1) * nnbasx + 1
            call fckden(.false., .true., dummy, udtvao, cmo,
     &                  dtv(:,i), wrk, nwrk)
            call dgefsp(nbast, udtvao, fdtvaos(j))
        end do
        deallocate(udtvao)

        allocate(fxcaos(ncsim*nnbasx))
        call pe_master(runtype='response', denmats=fdtvaos,
     &                 fckmats=fxcaos, nmats=ncsim, dalwrk=wrk(1:nwrk))
        deallocate(fdtvaos)

        do i = 1, ncsim
            j = (i - 1) * nnbasx + 1
            call uthu(fxcaos(j), fxcs(:,i), cmo, wrk, nbast, norbt)
            if (nasht > 0) call getac2(fxcs(:,i), fxcacs(:,i))
            tfxc = solelm(dv, fxcacs(:,i), fxcs(:,i), tfxcacs(i))
        end do
        deallocate(fxcaos)
      end if

      allocate(fyc(nnorbx), fycac(nnashx))
      call pe_fckmo(cmo, fyc, energy, dv, wrk, nwrk)
      if (nasht > 0) call getac2(fyc, fycac)
      tfyc = solelm(dv, fycac, fyc, tfycac)

!    ...CSF part of sigma vectors    
      call solsc(ncsim, 0, bcvecs, cref, scvecs, fxcacs, fycac, tfxcacs,
     &          tfycac, cindx, wrk, nwrk)
      deallocate(fxcacs, fycac, tfxcacs)

      if (nwoppt > 0) then
        mwoph = nwoph
        nwoph = nwoppt
        jscvec = 1 + nconst
        do i = 1, ncsim
            if (pe_polar) then
                call solgo(2.0d0, dv, fxcs(:,i), scvecs(jscvec,i))
            end if
            call solgo(0.0d0, dtv(:,i), fyc, scvecs(jscvec,i))
        end do
        nwoph = mwoph
      end if
      deallocate(fyc, fxcs)

      call qexit('pe_lnc')

      end subroutine pe_lnc

!------------------------------------------------------------------------------

      subroutine pe_lno(nosim, bovecs, cref, cmo, cindx, dv,
     &                  sovecs, nso, wrk, nwrk)
!
!  Written by Erik Donovan Hedegaard and Jógvan Magnus H. Olsen
!             after original code by Hans Jorgen Aa. Jensen
!
!  Purpose:  Calculate Hessian contribution from a
!            PE potential to an orbital trial vector.
!
!  NSVEC     may be NVAR or NWOPT, dependent on LINTRN
!
!
!  Used from common blocks:
!    INFORB : NNASHX, NNORBX, NNBASX, etc.
!    INFVAR : JWOP
!    INFLIN : NWOPPT, NVARPT, NCONST, NCONRF
!

      use polarizable_embedding, only: pe_master, pe_polar

      implicit none

#include "priunit.h"
#include "dummy.h"
#include "inforb.h"
#include "infvar.h"
#include "inflin.h"

      integer :: nosim, nso, nwrk
      real*8, dimension(*) :: cref, cmo, dv, cindx
      real*8, dimension(nwrk) :: wrk
      real*8, dimension(nwoppt,*) :: bovecs
      real*8, dimension(nso,*) :: sovecs

      integer :: i, j, jsovec, mwoph, ncolim
      logical :: fulhes, fndlab
      real*8 :: solelm
      real*8 :: txyo, energy
      real*8, dimension(:), allocatable :: txyoacs
      real*8, dimension(:), allocatable :: ubodcao, ubodvao
      real*8, dimension(:), allocatable :: bodtaos, fxoaos
      real*8, dimension(:), allocatable :: fckmo, fyo, ufyo
      real*8, dimension(:,:), allocatable :: ubovecs, fxos
      real*8, dimension(:,:), allocatable :: fxyos, fxyoacs


      call qenter('pe_lno')

      allocate(ubovecs(n2orbx,nosim))
      if (nosim > 0) then
          do i = 1, nosim
            call upkwop(nwoppt, jwop, bovecs(:,i), ubovecs(:,i))
          end do
      end if

      if (pe_polar) then
        allocate(ubodcao(n2basx), ubodvao(n2basx))
        allocate(bodtaos(nosim*nnbasx))
        do i = 1, nosim
            j = (i - 1) * nnbasx + 1
            call tr1den(cmo, ubovecs(:,i), dv, ubodcao, ubodvao,
     &                  wrk, nwrk)
            if (nasht > 0) ubodcao = ubodcao + ubodvao
            call dgefsp(nbast, ubodcao, bodtaos(j))
        end do
        deallocate(ubodcao, ubodvao)

        allocate(fxoaos(nosim*nnbasx))
        call pe_master(runtype='response', denmats=bodtaos,
     &                 fckmats=fxoaos, nmats=nosim, dalwrk=wrk(1:nwrk))
        deallocate(bodtaos)

        allocate(fxos(nnorbx,nosim)) 
        do i = 1, nosim
            j = (i - 1) * nnbasx + 1
            call uthu(fxoaos(j), fxos(:,i), cmo, wrk, nbast, norbt)
        end do
        deallocate(fxoaos)
      end if

      allocate(fckmo(nnorbx))
      call pe_fckmo(cmo, fckmo, energy, dv, wrk, nwrk)

      allocate(fyo(n2orbx))
      call dsptsi(norbt, fckmo, fyo)
      deallocate(fckmo)

      allocate(ufyo(n2orbx), txyoacs(nosim))
      allocate(fxyos(nnorbx,nosim), fxyoacs(nnashx,nosim))
      do i = 1, nosim
        ufyo = 0.0d0
        call tr1uh1(ubovecs(:,i), fyo, ufyo, 1)
        call dgetsp(norbt, ufyo, fxyos(:,i))
        if (pe_polar) then
            call daxpy(nnorbx, 1.0d0, fxos(:,i), 1, fxyos(:,i), 1)
        end if
        if (nasht > 0) then
            call getac2(fxyos(:,i), fxyoacs(:,i))
        end if
        txyo = solelm(dv, fxyoacs(:,i), fxyos(:,i), txyoacs(i))
      end do

!     ... CSF part of sigma vectors
      if (lsymrf == lsymst) then
        ncolim = 1
      else
        ncolim = 0
      end if

!     Determine if full Hessian or only orbital Hessian
      fulhes = (nso == nvarpt)
      if (fulhes) then
        jsovec = 1 + nconst
      else
        jsovec = 1
      end if

      if (fulhes .and. (nconst > ncolim)) then
        call solsc(0, nosim, dummy, cref, sovecs, fxyoacs, dummy,
     &             txyoacs, dummy, cindx, wrk, nwrk)
      end if

!     ... orbital part of sigma vectors
      mwoph = nwoph
      nwoph = nwoppt
!    ... tell SOLGO only to use the NWOPPT first JWOP entries
      do i = 1, nosim
        call solgo(2.0d0, dv, fxyos(:,i), sovecs(jsovec,i))
      end do
      nwoph = mwoph

      call qexit('pe_lno')

      end subroutine pe_lno

!------------------------------------------------------------------------------

      subroutine pe_rsplin(ncsim, nosim, bcvecs, bovecs, cref, cmo,
     &                     cindx, udv, dv, udvtr, dvtr, dtv, dtvtr,
     &                     scvecs, sovecs, wrk, nwrk)

        implicit none

#include "priunit.h"
#include "dummy.h"
#include "wrkrsp.h"
#include "infrsp.h"
#include "inftap.h"

        integer :: ncsim, nosim, nwrk
        real*8, dimension(*) :: bcvecs, bovecs
        real*8, dimension(*) :: cref, cmo, cindx, udv, dv
        real*8, dimension(*) :: udvtr, dvtr, dtv, dtvtr
        real*8, dimension(*) :: scvecs, sovecs
        real*8, dimension(nwrk) :: wrk

        call qenter('pe_rsplin')

        if (kzconf >  0) then
            stop 'ERROR: PE-MCSCF response not implemented.'
        end if

!        if (ncsim > 0) then
!            call pe_rsplnc(ncsim, bcvecs, cref, cmo, cindx, udv, dv,
!     &                     udvtr, dvtr, dtv, dtvtr, scvecs, wrk, nwrk)
!        end if

        if (nosim > 0) then
            call pe_rsplno(nosim, bovecs, cref, cmo, cindx, udv, dv,
     &                     udvtr, dvtr, sovecs, wrk, nwrk)
        end if

        call qexit('pe_rsplin')



      end subroutine pe_rsplin

!------------------------------------------------------------------------------
!
!      subroutine pe_rsplnc(ncsim, bcvecs, cref, cmo, cindx, udv, dv,
!     &                     udvtr, dvtr, dtv, dtvtr, scvecs, wrk, nwrk)
!
!        implicit none
!
!#include "priunit.h"
!#include "dummy.h"
!#include "infrsp.h"
!#include "inftap.h"
!
!
!      end subroutine pe_rsplnc
!
!------------------------------------------------------------------------------

      subroutine pe_rsplno(nosim, bovecs, cref, cmo, cindx, udv, dv,
     &                     udvtr, dvtr, sovecs, wrk, nwrk)

        use polarizable_embedding, only: pe_master, pe_polar, pe_gspol

        implicit none

#include "priunit.h"
#include "dummy.h"
#include "wrkrsp.h"
#include "inforb.h"
#include "infrsp.h"
!#include "inftap.h"

        integer :: nosim, nwrk
        real*8, dimension(*) :: bovecs
        real*8, dimension(kzyvar,*) :: sovecs
        real*8, dimension(*) :: cref, cmo, cindx, udv, dv, udvtr, dvtr
        real*8, dimension(nwrk) :: wrk

        integer :: i, j
        real*8, dimension(:), allocatable :: dcao, dvao
!        real*8, dimension(:), allocatable :: dcaotr, dvaotr
        real*8, dimension(:), allocatable :: daos, fckaos
!        real*8, dimension(:), allocatable :: daotrs, fckaotrs
        real*8, dimension(:), allocatable :: evec
        real*8, dimension(:,:), allocatable :: ubovecs, evecs, eacs
!        real*8, dimension(:,:), allocatable :: evectrs, eactrs

        ! return if no polarization and not MCSCF
        if ((kzconf == 0) .and. .not. pe_polar) then
            return
        ! no polarization for triplet excitations in closed shell
        else if ((nasht == 0) .and. trplet) then
            return
        ! ground state polarization
        else if (pe_gspol) then
            return
        ! triplet response for open shell systems not ready yet
        else if ((nasht > 0) .and. trplet) then
            stop 'ERROR: triplet response for open shell systems not
     & ready yet'
        end if

        call qenter('pe_rsplno')

        allocate(ubovecs(n2orbx,nosim))
        call rspzym(nosim, bovecs, ubovecs)

        ubovecs = - 2.0d0 * ubovecs

        allocate(dcao(n2basx), dvao(n2basx), daos(nosim*nnbasx))
!        if (trplet) then
!            allocate(dcaotr(n2basx), dvaotr(n2basx),
!     &               daotrs(nosim*nnbasx))
!        end if
        do i = 1, nosim
            j = 1 + (i - 1) * nnbasx
            call deq27(cmo, ubovecs(:,i), udv, dcao, dvao, wrk, nwrk)
!            if (trplet) then
!                call deq27(cmo, ubovecs(:,i), udvtr, dcaotr, dvaotr,
!     &                     wrk, nwrk)
!            end if
            ! factor 1/2 probably related to factor 2 above?
            if (nasht > 0) then
                dcao = dcao + 0.5d0 * dvao
!                if (trplet) then
!                    dcaotr = dcaotr + 0.50d0 * dvaotr
!                end if
            end if
            call dgefsp(nbast, dcao, daos(j))
!            if (trplet) then
!                call dgefsp(nbast, dcaotr, daotrs(j))
!            end if
        end do
        deallocate(dcao, dvao)
!        if (trplet) then
!            deallocate(dcaotr, dvaotr)
!        end if

        allocate(fckaos(nosim*nnbasx))
        call pe_master(runtype='response', denmats=daos,
     &                 fckmats=fckaos, nmats=nosim,
     &                 dalwrk=wrk)
        deallocate(daos)
!        if (trplet) then
!            allocate(fckaotrs(nosim*nnbasx))
!            call pe_master(runtype='response', denmats=daotrs,
!     &                     fckmats=fckaotrs, nmats=nosim,
!     &                     dalwrk=wrk)
!            deallocate(daotrs)
!        end if

        allocate(evec(nnorbx))
        allocate(evecs(n2orbx,nosim), eacs(n2ashx,nosim))
!        if (trplet) then
!            allocate(evectrs(n2orbx,nosim), eactrs(n2ashx,nosim))
!        end if
        do i = 1, nosim
            j = 1 + (i - 1) * nnbasx
            call uthu(fckaos(j), evec, cmo, wrk, nbast, norbt)
            call dsptsi(norbt, evec, evecs(:,i))
!            if (trplet) then
!                call uthu(fckaotrs(j), evec, cmo, wrk, nbast, norbt)
!                call dsptsi(norbt, evec, evectrs(:,i))
!            end if
!            if (nasht > 0) then
!                call getacq(evecs(:,i), eacs(:,i))
!                if (trplet) then
!                    call getacq(evectrs(:,i), eactrs(:,i))
!                end if
!            end if
!            tr = solelm(dv, fxyoacs(:,i), fxyos(:,i), txyoacs(i))
        end do
        deallocate(evec)

!        if (trplet) then
!            call slvsor(.true., .false., nosim, udvtr, sovecs, evectrs)
!            call slvsor(.true., .true., nosim, udv, sovecs, evecs)
!        else
            call slvsor(.true., .true., nosim, udv, sovecs, evecs)
!        end if

        deallocate(evecs, eacs)
!        if (trplet) then
!            deallocate(evectrs, eactrs)
!        end if

        call qexit('pe_rsplno')

      end subroutine pe_rsplno

!------------------------------------------------------------------------------
!
!      subroutine pe_rspqro()
!
!      end subroutine pe_rspqro
!
