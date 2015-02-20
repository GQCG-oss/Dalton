!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2015 (2015), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
#if defined (BUILD_PELIB)
      subroutine pe_ifc_input_reader(word)
        use polarizable_embedding, only: pe_input_reader
        implicit none
#include "priunit.h"
        character(len=7), intent(inout) :: word
        call qenter('pe_ifc_input_reader')
        call pe_input_reader(word, lucmd)
        call qexit('pe_ifc_input_reader')
      end subroutine pe_ifc_input_reader

!-----------------------------------------------------------------------------

      subroutine pe_ifc_init()
        use polarizable_embedding, only: pe_init
        implicit none
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
        call qenter('pe_ifc_init')
        call pe_init(lupri, cord(1:3,1:natoms), charge(1:natoms))
        call qexit('pe_ifc_init')
      end subroutine pe_ifc_init

!-----------------------------------------------------------------------------

      subroutine pe_ifc_finalize()
        use polarizable_embedding, only: pe_finalize
        implicit none
        call qenter('pe_ifc_finalize')
        call pe_finalize()
        call qexit('pe_ifc_finalize')
      end subroutine pe_ifc_finalize

!-----------------------------------------------------------------------------

      subroutine pe_ifc_fock(denmats, fckmats, nmats, energies)
        use polarizable_embedding, only: pe_master
        implicit none
#include "inforb.h"
        integer, intent(in) :: nmats
        real*8, dimension(nmats*nnbasx), intent(in) :: denmats
        real*8, dimension(nmats*nnbasx), intent(out) :: fckmats
        real*8, dimension(nmats), intent(out) :: energies
        call qenter('pe_ifc_fock')
#if defined(VAR_MPI)
        call pe_ifc_start_slaves(1)
#endif
        call pe_master(runtype='full_fock',
     &                 triang=.true.,
     &                 ndim=nbast,
     &                 nmats=nmats,
     &                 denmats=denmats,
     &                 fckmats=fckmats,
     &                 expvals=energies)
        call qexit('pe_ifc_fock')
      end subroutine pe_ifc_fock

!-----------------------------------------------------------------------------

      subroutine pe_ifc_energy(denmats, nmats)
        use polarizable_embedding, only: pe_master
        implicit none
#include "inforb.h"
        integer, intent(in) :: nmats
        real*8, dimension(nmats*nnbasx), intent(in) :: denmats
        call qenter('pe_ifc_energy')
#if defined(VAR_MPI)
        call pe_ifc_start_slaves(2)
#endif
        call pe_master(runtype='print_energy',
     &                 triang=.true.,
     &                 ndim=nbast,
     &                 nmats=nmats,
     &                 denmats=denmats)
        call qexit('pe_ifc_energy')
      end subroutine pe_ifc_energy

!-----------------------------------------------------------------------------

      subroutine pe_ifc_london(fckmats)
        use polarizable_embedding, only: pe_master
        implicit none
#include "inforb.h"
        real*8, dimension(3*n2basx), intent(out) :: fckmats
        integer :: i, j, k, l, m
        real*8, dimension(:), allocatable :: fckmats_packed
        call qenter('pe_ifc_london')
        allocate(fckmats_packed(3*nnbasx))
#if defined(VAR_MPI)
        call pe_ifc_start_slaves(4)
#endif
        call pe_master('magnetic_gradient',
     &                 triang=.true.,
     &                 ndim=nbast,
     &                 fckmats=fckmats_packed)
        do i = 1, 3
            j = (i - 1) * nnbasx + 1
            k = i * nnbasx
            l = (i - 1) * n2basx + 1
            m = i * n2basx
            call daptge(nbas, fckmats_packed(j:k), fckmats(l:m))
        end do
        deallocate(fckmats_packed)
        call qexit('pe_ifc_london')
      end subroutine pe_ifc_london

!-----------------------------------------------------------------------------

#if defined(VAR_MPI)
      subroutine pe_ifc_slave(runtype)
        use polarizable_embedding, only: pe_slave
        implicit none
        integer, intent(in) :: runtype
        call qenter('pe_ifc_slave')
        if (runtype == 1) then
            call pe_slave('full_fock')
        else if (runtype == 2) then
            call pe_slave('print_energy')
        else if (runtype == 3) then
            call pe_slave('dynamic_response')
        else if (runtype == 4) then
            call pe_slave('magnetic_gradient')
        end if
        call qexit('pe_ifc_slave')
      end subroutine pe_ifc_slave

!-----------------------------------------------------------------------------

      subroutine pe_ifc_start_slaves(runtyp)
        implicit none
        integer :: runtyp
#include "iprtyp.h"
#include "maxorb.h"
#include "infpar.h"
        integer, parameter :: iprtyp = POLARIZABLE_EMBEDDING
        call qenter('pe_ifc_start_slaves')
        if (master /= 0) call quit('ERROR: PElib assumes master id 0')
        if (nodtot >= 1) then
            call mpixbcast(iprtyp, 1, 'INTEGER', master)
            call mpixbcast(runtyp, 1, 'INTEGER', master)
        end if
        call qexit('pe_ifc_start_slaves')
      end subroutine pe_ifc_start_slaves
#endif

!-----------------------------------------------------------------------------

      subroutine pe_grad(cref, cmo, cindx, dv, grd, energy, wrk, nwrk)
!
!     Written by Erik Donovan Hedegård (edh) and Jogvan Magnus H. Olsen
!                based on PCMGRAD
!
!     Purpose:  calculate (MCSCF) energy and gradient contribution
!               from an embedding potential using the PE library
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
        real*8, dimension(1) :: etmp
        real*8, dimension(:), allocatable :: fckmo, fckac
        real*8, dimension(:), allocatable :: pegrd, diape

        call qenter('pe_grad')

        allocate(fckmo(nnorbx))
        call pe_fckmo(cmo, dv, fckmo, etmp, wrk, nwrk)
        energy = etmp(1)

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

        !--------------- Orthogonality test ----------------
        test = ddot(nconf, cref, 1, pegrd, 1)
        if (abs(test) > 1.0d-8) then
            nwarn = nwarn + 1
            write(lupri,*) ' >>> PE GRADIENT WARNING <<< '
            write(lupri,*) ' < CREF | GRAD > =', test
        end if

        ! Add PE gradient contribution to MCSCF gradient
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

      subroutine pe_fckmo(cmo, dv, fckmo, energy, wrk, nwrk)
!
!     Written by Erik Donovan Hedegård (edh) and Jogvan Magnus H. Olsen
!
!     Purpose:  Compute PE operator in MO basis using PE library
!
!     Output:
!     fckmo     PE operator in MO basis
!     energy    PE energy
!
!     Used from common blocks:
!     INFVAR: NCONF,  NWOPT,  NVAR,   NVARH
!     INFORB: NNASHX, NNBASX, NNORBX, etc.
!     INFTAP: LUIT2
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

        real*8, dimension(:), allocatable :: dcao, dvao, fdtao, fckao

        call qenter('pe_fckmo')

        allocate(dcao(n2basx), dvao(n2basx))
        call fckden((nisht > 0), (nasht > 0), dcao, dvao, cmo, dv,
     &              wrk, nwrk)
        if (nisht == 0) dcao = 0.0d0
        if (nasht > 0) dcao = dcao + dvao
        deallocate(dvao)

        allocate(fdtao(nnbasx))
        call dgefsp(nbast, dcao, fdtao)
        deallocate(dcao)

        allocate(fckao(nnbasx))
#if defined(VAR_MPI)
        call pe_ifc_start_slaves(1)
#endif
        call pe_master(runtype='full_fock',
     &                 triang=.true.,
     &                 ndim=nbast,
     &                 nmats=1,
     &                 denmats=fdtao,
     &                 fckmats=fckao,
     &                 expvals=energy)
        deallocate(fdtao)

        call uthu(fckao, fckmo, cmo, wrk, nbast, norbt)
        deallocate(fckao)

        call qexit('pe_fckmo')

      end subroutine pe_fckmo

!------------------------------------------------------------------------------

      subroutine pe_lin(ncsim, nosim, bcvecs, bovecs, cref, cmo, cindx,
     &                  dv, dtv, scvecs, sovecs, orblin, wrk, nwrk)
!
! Written by Erik Donovan Hedegård and Jogvan Magnus H. Olsen
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
            call pe_lnc(ncsim, bcvecs, cref, cmo, cindx, dv, dtv,
     &                  scvecs, wrk, nwrk)
        end if

        if (nosim > 0) then
            if (.not. orblin) then
                nso = nvarpt
            else
                nso = nwoppt
            end if
            call pe_lno(nosim, bovecs, cref, cmo, cindx, dv, sovecs,
     &                  nso, wrk, nwrk)
        end if

        call qexit('pe_lin')

      end subroutine pe_lin

!------------------------------------------------------------------------------

      subroutine pe_lnc(ncsim, bcvecs, cref, cmo, cindx, dv, dtv,
     &                  scvecs, wrk, nwrk)
!
!  Written by Erik Donovan Hedegaard and Jogvan Magnus H. Olsen
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

        use pe_variables, only: pe_polar
        use polarizable_embedding, only: pe_master

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
        real*8 :: tfxc, tfyc, tfycac, solelm
        real*8, dimension(1) :: energy
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
     &                      dtv(:,i), wrk, nwrk)
                call dgefsp(nbast, udtvao, fdtvaos(j))
            end do
            deallocate(udtvao)

            ! Fxc = -R<0|Fe|B>Fe in fxcaos
            allocate(fxcaos(ncsim*nnbasx))
#if defined(VAR_MPI)
            call pe_ifc_start_slaves(3)
#endif
            call pe_master(runtype='dynamic_response',
     &                     triang=.true.,
     &                     ndim=nbast,
     &                     nmats=ncsim,
     &                     denmats=fdtvaos,
     &                     fckmats=fxcaos)
            deallocate(fdtvaos)

            do i = 1, ncsim
                j = (i - 1) * nnbasx + 1
                call uthu(fxcaos(j), fxcs(:,i), cmo, wrk, nbast, norbt)
                if (nasht > 0) call getac2(fxcs(:,i), fxcacs(:,i))
                tfxc = solelm(dv, fxcacs(:,i), fxcs(:,i), tfxcacs(i))
            end do
            deallocate(fxcaos)
        end if

        ! Fg = Vmul -R<0|Fe|0>Fe in fyc
        allocate(fyc(nnorbx), fycac(nnashx))
        call pe_fckmo(cmo, dv, fyc, energy, wrk, nwrk)
        if (nasht > 0) call getac2(fyc, fycac)
        tfyc = solelm(dv, fycac, fyc, tfycac)

!       ...CSF part of sigma vectors
        call solsc(ncsim, 0, bcvecs, cref, scvecs, fxcacs, fycac,
     &             tfxcacs, tfycac, cindx, wrk, nwrk)
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
!  Written by Erik Donovan Hedegaard and Jogvan Magnus H. Olsen
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

        use pe_variables, only: pe_polar
        use polarizable_embedding, only: pe_master

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
        real*8 :: txyo
        real*8, dimension(1) :: energy
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

        ! 1. Calculation of Fxo = R*<0|Fe(k)|O>Fe
        !    Store in fxos
        if (pe_polar) then
            allocate(ubodcao(n2basx), ubodvao(n2basx))
            allocate(bodtaos(nosim*nnbasx))
            do i = 1, nosim
                j = (i - 1) * nnbasx + 1
                call tr1den(cmo, ubovecs(:,i), dv, ubodcao, ubodvao,
     &                      wrk, nwrk)
                if (nasht > 0) ubodcao = ubodcao + ubodvao
                call dgefsp(nbast, ubodcao, bodtaos(j))
            end do
            deallocate(ubodcao, ubodvao)

            allocate(fxoaos(nosim*nnbasx))
#if defined(VAR_MPI)
            call pe_ifc_start_slaves(3)
#endif
            call pe_master(runtype='dynamic_response',
     &                     triang=.true.,
     &                     ndim=nbast,
     &                     nmats=nosim,
     &                     denmats=bodtaos,
     &                     fckmats=fxoaos)
            deallocate(bodtaos)


            allocate(fxos(nnorbx,nosim))
            do i = 1, nosim
                j = (i - 1) * nnbasx + 1
                call uthu(fxoaos(j), fxos(:,i), cmo, wrk, nbast, norbt)
            end do
            deallocate(fxoaos)
        end if

        ! 2. Calculation of Fyo = V(k) + R<0|F|0>Fe(k)
        !    Store in fyos
        allocate(fckmo(nnorbx))
        call pe_fckmo(cmo, dv, fckmo, energy, wrk, nwrk)

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
        ! 3.   /        <0[Epq,Fxo + Fyo]|0>      \  orbital part
        !      \ 2<0|Fyo + Fxo|mu> - <0|Fyo|0>*c0 /  CSF part
        !     ... CSF part of sigma vectors
        if (lsymrf == lsymst) then
            ncolim = 1
        else
            ncolim = 0
        end if

        ! Determine if full Hessian or only orbital Hessian
        fulhes = (nso == nvarpt)
        if (fulhes) then
            jsovec = 1 + nconst
        else
            jsovec = 1
        end if

        if (fulhes .and. (nconst > ncolim)) then
            call solsc(0, nosim, dummy, cref, sovecs, fxyoacs, dummy,
     &                 txyoacs, dummy, cindx, wrk, nwrk)
        end if

        ! ... orbital part of sigma vectors
        mwoph = nwoph
        nwoph = nwoppt
        ! ... tell SOLGO only to use the NWOPPT first JWOP entries
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

        if (ncsim > 0 .and. .not. soppa) then
            call pe_rsplnc(ncsim, bcvecs, cref, cmo, cindx, udv, dv,
     &                     udvtr, dvtr, dtv, dtvtr, scvecs, wrk, nwrk)
        end if

        if (nosim > 0) then
            call pe_rsplno(ncsim, nosim, bovecs, cref, cmo, cindx,
     &                     udv, dv, udvtr, dvtr, sovecs, wrk, nwrk)
        end if

        call qexit('pe_rsplin')

      end subroutine pe_rsplin

!------------------------------------------------------------------------------
!
       subroutine pe_rsplnc(ncsim, bcvecs, cref, cmo, cindx, udv, dv,
     &                      udvtr, dvtr, dtv, dtvtr, scvecs, wrk, nwrk)

        use pe_variables, only: pe_polar
        use polarizable_embedding, only: pe_master

        implicit none
!
#include "priunit.h"
#include "dummy.h"
#include "infrsp.h"
#include "inftap.h"
#include "wrkrsp.h"
#include "inforb.h"
#include "qrinf.h"
#include "infvar.h"
! Used from common: n2ashx, kzyvar, kzconf, ncref, maxwop

        integer :: i, j
        integer :: ncsim, nwrk

        real*8, dimension(*) :: bcvecs, cref, cmo, cindx, udv, dv
        real*8, dimension(*) :: udvtr, dvtr
        real*8, dimension(n2ashx,*) :: dtv
        real*8, dimension(n2ashx,*) :: dtvtr
        real*8, dimension(kzyvar,*) :: scvecs
        real*8, dimension(nwrk) :: wrk

        real*8 :: ovlap, solelm, tfpeac, tfpe
        real*8, dimension(:,:), allocatable :: udtv!, udtvao
        real*8, dimension(:,:), allocatable :: fxcs, fxcacs
        real*8, dimension(:), allocatable :: udtvao, fuxcs
        real*8, dimension(:), allocatable :: fdtvaos, fxcaos
        real*8, dimension(:), allocatable :: tfxc, tfxcacs
        real*8, dimension(:), allocatable :: fpe, fupe, fpeac

        logical :: lexist, lopen, locdeb
        logical :: fndlab
        logical :: tdm, norho2

        call qenter('pe_rsplnc')

        locdeb = .false.

        lopen = .false.
        tdm = .true.
        norho2 = .true.

        allocate(fxcs(nnorbx,ncsim))
        allocate(fuxcs(n2orbx))
        allocate(fxcacs(nnashx,ncsim))
        allocate(tfxc(ncsim))
        allocate(tfxcacs(ncsim))
        fxcs    = 0.0d0
        fuxcs   = 0.0d0
        fxcacs  = 0.0d0
        tfxc    = 0.0d0
        tfxcacs = 0.0d0

        ! Fxc = R*(<0(L)|Fe|0> + <0|Fe|0(R)>)Fe
        if (pe_polar) then
            call getref(cref, ncref)
            ! ...Construct <0(L)|...|0> + <0|...|0(R)>
            allocate(udtv(n2ashx,ncsim))
            udtv = 0.0d0
            call rsptdm(ncsim,irefsy,ksymst,ncref,kzconf,cref,
     &                  bcvecs,udtv,dummy,0,0,tdm,norho2,
     &                  cindx,wrk,1,nwrk)
            udtv = -1.0d0*udtv

            if ( ncsim > 0 ) then
                allocate(fdtvaos(nnbasx*ncsim))
                fdtvaos = 0.0d0
                allocate(udtvao(n2basx))
                udtvao = 0.0d0
                do i = 1, ncsim
                   j = (i - 1) * nnbasx + 1
                    call fckden2(.false.,.true., dummy, udtvao, cmo,
     &                           udtv(:,i), wrk, nwrk)
                    call dgefsp(nbast, udtvao, fdtvaos(j))
                end do
                deallocate(udtv,udtvao)
            end if

            allocate(fxcaos(ncsim*nnbasx))
            fxcaos = 0.0d0
#if defined(VAR_MPI)
            call pe_ifc_start_slaves(3)
#endif
            call pe_master(runtype='dynamic_response',
     &                     triang=.true.,
     &                     ndim=nbast,
     &                     nmats=ncsim,
     &                     denmats=fdtvaos,
     &                     fckmats=fxcaos)
            deallocate(fdtvaos)

            do i = 1, ncsim
                j = (i - 1) * nnbasx + 1
                call uthu(fxcaos(j), fxcs(:,i), cmo, wrk, nbast, norbt)
                if (nasht > 0) call getac2(fxcs(:,i), fxcacs(:,i))
                tfxc = solelm(dv, fxcacs(:,i), fxcs(:,i), tfxcacs(i))
            end do
            deallocate(fxcaos)
        end if

        ! Fg = V - <0|F|0>Fe -unpack into fupe
        if (.not. tdhf) then
            allocate(fpe(nnorbx))
            fpe = 0.0d0
            if (lusifc <= 0) then
                call gpopen(lusifc, 'SIRIFC', 'OLD', ' ', 'UNFORMATTED',
     &                      idummy, .false.)
                lopen = .true.
            end if
            rewind(lusifc)
            call mollab('PEFMAT  ', lusifc, lupri)
            call readt(lusifc, nnorbx, fpe)
            if (lopen) call gpclose(lusifc, 'KEEP')
            allocate(fupe(n2orbx), fpeac(nnashx))
            fupe = 0.0d0
            fpeac = 0.0d0
            call dsptsi(norbt, fpe, fupe)
            if (nasht > 0) call getac2(fpe, fpeac)
                tfpe = solelm(dv, fpeac, fpe, tfpeac)
            deallocate(fpe)
        end if

!       Calculate Fxc(Rxc) and Fg(Ryc) contributions to SCVECS(NVAR,NCSIM)
!       ... CSF part of sigma vectors
        if (locdeb) then
           write(lupri,*)' Linear transformed configuration vector'
           write(lupri,*)' **** Before slvsc in pe_rsplnc **** '
           call output(scvecs,1,kzyvar,1,ncsim,kzyvar,ncsim,1,lupri)
        endif

        call slvsc(ncsim, 0, nnashx, bcvecs, cref, scvecs, fxcacs,
     &             fpeac, tfxcacs, tfpeac, cindx, wrk, nwrk)
        deallocate(fxcacs, tfxcacs, fpeac)

        if (locdeb) then
            write(lupri,*)' Linear transformed configuration vector'
            write(lupri,*)' **** After slvsc in pe_rsplnc **** '
            call output(scvecs,1,kzyvar,1,ncsim,kzyvar,ncsim,1,lupri)
        end if

!       ... orbital part of sigma vector(s)
        if (kzwopt .gt. 0) then
            do i = 1,ncsim
                fuxcs = 0.0d0
                call dsptsi(norbt,fxcs(:,i), fuxcs)
                call slvsor(.true.,.true., 1, udv, scvecs(1,i), fuxcs)
                if (locdeb) then
                    write(lupri,*)' *** After slvsor in pe_rsplnc *** '
                    write(lupri,*)
     &              'Orbital part of linear transformed conf vec no ', i
                    write(lupri,*)' Txc contribution'
                    call output(scvecs(1,i), 1, kzyvar, 1, 1, kzyvar,
     &                          1, 1,lupri)
                end if

                call slvsor(.false., .false., 1, dtv(1,i), scvecs(1,i),
     &                      fupe)
                if (locdeb) then
                    write(lupri,*)
     &              'Orbital part of linear transformed conf vec no ', i
                    write(lupri,*)' Tg contribution'
                    call output(scvecs(1,i), 1, kzyvar, 1, 1, kzyvar,
     &                          1, 1, lupri)
                end if
            end do
            deallocate(fupe, fuxcs)

            if (locdeb) then
                write(lupri,*)' linear transformed conf. vector'
                write(lupri,*)' *** after slvsor in pe_rsplnc *** '
                call output(scvecs, 1, kzyvar, 1, ncsim, kzyvar, ncsim,
     &                     1, lupri)
            end if
        end if

        if (ncref /= kzconf) then
            call quit('ERROR in pe_rsplnc: ncref /= kzconf')
        end if

        call qexit('pe_rsplnc')

      end subroutine pe_rsplnc

!-------------------------------------------------------------------------------

      subroutine pe_rsplno(ncsim, nosim, bovecs, cref, cmo, cindx,
     &                     udv, dv, udvtr, dvtr, sovecs, wrk, nwrk)

        use pe_variables, only: pe_polar, pe_gspol
        use polarizable_embedding, only: pe_master

        implicit none

#include "priunit.h"
#include "dummy.h"
#include "wrkrsp.h"
#include "inforb.h"
#include "infrsp.h"
#include "inftap.h"

        integer :: nosim, ncsim, nwrk
        real*8, dimension(*) :: bovecs
        real*8, dimension(kzyvar,*) :: sovecs
        real*8, dimension(*) :: cref, cmo, cindx, udv, dv, udvtr, dvtr
        real*8, dimension(nwrk) :: wrk

        integer :: i, j
        real*8 :: txyo
        real*8 :: ddot, slvqlm
        real*8, dimension(:), allocatable :: dcao, dvao
        real*8, dimension(:), allocatable :: daos, fckaos
        real*8, dimension(:), allocatable :: daotrs
        real*8, dimension(:), allocatable :: evec
        real*8, dimension(:,:), allocatable :: ubovecs, evecs, eacs
        real*8, dimension(:), allocatable :: fpemo,fupemo
        real*8, dimension(:), allocatable :: txyoacs
        real*8, dimension(:), allocatable :: ovlp
        logical :: lexist, lopen

        call qenter('pe_rsplno')

        ! return if no polarization and not MCSCF
        if (tdhf .and. .not. pe_polar) then
            call qexit('pe_rsplno')
            return
        ! no polarization for triplet excitations in closed shell SCF
        else if ((nasht == 0) .and. trplet) then
            call qexit('pe_rsplno')
            return
        ! ground state polarization approximation
        else if (pe_gspol) then
            call qexit('pe_rsplno')
            return
        ! triplet response for open shell systems not ready yet
        else if ((nasht > 0) .and. trplet) then
            call quit('ERROR: triplet operators for open shell'//
     &                ' systems not implemented')
        end if

        lopen = .false.

        if (.not. tdhf) then
            ! Read Fg = V - <0|F|0>Fe from file
            allocate(fpemo(nnorbx))
            if (lusifc <= 0) then
                call gpopen(lusifc, 'SIRIFC', 'OLD', ' ', 'UNFORMATTED',
     &                      idummy, .false.)
                lopen = .true.
            end if
            rewind(lusifc)
            call mollab('PEFMAT  ', lusifc, lupri)
            call readt(lusifc, nnorbx, fpemo)
            if (lopen) call gpclose(lusifc, 'KEEP')
            allocate(fupemo(n2orbx))
            call dsptsi(norbt, fpemo, fupemo)
            deallocate(fpemo)
        end if

        allocate(ubovecs(n2orbx,nosim))
        call rspzym(nosim, bovecs, ubovecs)

        ubovecs = - ubovecs

        allocate(dcao(n2basx), dvao(n2basx), daos(nosim*nnbasx))
        ! Calculate Fxo = <0|Fe(k)|0>Fe
        do i = 1, nosim
            j = (i - 1) * nnbasx + 1
            call deq27(cmo, ubovecs(:,i), udv, dcao, dvao,
     &                 wrk, nwrk)
            if (nasht > 0) then
                dcao = dcao + dvao
            end if
            call dgefsp(nbast, dcao, daos(j))
        end do
        deallocate(dcao, dvao)

        allocate(fckaos(nosim*nnbasx))
#if defined(VAR_MPI)
        call pe_ifc_start_slaves(3)
#endif
        call pe_master(runtype='dynamic_response',
     &                 triang=.true.,
     &                 ndim=nbast,
     &                 nmats=nosim,
     &                 denmats=daos,
     &                 fckmats=fckaos)
        deallocate(daos)

        allocate(evec(nnorbx))
        allocate(evecs(n2orbx,nosim))
        evecs = 0.0d0
        if (.not. tdhf) then
            allocate(eacs(n2ashx,nosim))
            allocate(txyoacs(nosim))
            eacs = 0.0d0
            txyoacs = 0.0d0
        end if

        do i = 1, nosim
            j = (i - 1) * nnbasx + 1
            call uthu(fckaos(j), evec, cmo, wrk, nbast, norbt)
            call dsptsi(norbt, evec, evecs(:,i))

            ! Fyo = V(k) - <0|F|0>Fe(k)
            if (.not. tdhf) then
                call onexh1(ubovecs(:,i), fupemo, evecs(:,i))
                call getacq(evecs(:,i), eacs(:,i))
                txyo = slvqlm(udv, eacs(:,i), evecs(:,i), txyoacs(i))
            end if
        end do

        deallocate(evec)
        if (.not. tdhf) then
            deallocate(fupemo)
            call slvsc(0, nosim, n2ashx, dummy, cref, sovecs, eacs,
     &                 dummy, txyoacs, dummy, cindx, wrk, nwrk)
            deallocate(eacs)
            deallocate(txyoacs)
        end if

        call slvsor(.true., .true., nosim, udv, sovecs, evecs)

        deallocate(evecs)

        call qexit('pe_rsplno')

      end subroutine pe_rsplno

!------------------------------------------------------------------------------

      subroutine pe_rspqro(vecb, vecc, etrs, xindx, zymb, zymc,
     &                 udv, wrk, nwrk, kzyva, kzyvb, kzyvc,
     &                 isyma, isymb, isymc, cmo, mjwop)

        use polarizable_embedding, only: pe_master

        implicit none

#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "infdim.h"
#include "qrinf.h"

        integer :: kzyva, kzyvb, kzyvc
        integer :: isyma, isymb, isymc
        integer :: nwrk
        real*8, dimension(nwrk) :: wrk
        real*8, dimension(kzyva) :: etrs
        real*8, dimension(kzyvb) :: vecb
        real*8, dimension(kzyvc) :: vecc
        real*8, dimension(ncmot) :: cmo
        real*8, dimension(norbt,norbt) :: zymb, zymc
        real*8, dimension(nashdi,nashdi) :: udv
        real*8, dimension(lcindx) :: xindx
        integer, dimension(2,maxwop,8) :: mjwop

        integer :: i, j, k
        integer :: idum = 1
        real*8, dimension(:), allocatable :: udcao, ufcmo
        real*8, dimension(:), allocatable :: dcaos, fcaos
        real*8, dimension(:), allocatable :: fcmo

        call qenter('pe_rspqro')

        call gtzymt(1, vecb, kzyvb, isymb, zymb, mjwop)
        call gtzymt(1, vecc, kzyvc, isymc, zymc, mjwop)

        allocate(udcao(n2basx))
        allocate(ufcmo(n2orbx))
        allocate(dcaos(4*nnbasx))
        dcaos = 0.0d0

        !  D(1k)
        udcao = 0.0d0
        call cdens1(isymb, cmo, zymb, udcao, wrk, nwrk)
        call dgefsp(nbast, udcao, dcaos(1:nnbasx))

        ! D(1k,2k)
        udcao = 0.0d0
        call cdens2(isymb, isymc, cmo, zymb, zymc, udcao,
     &              wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
        call dgefsp(nbast, udcao, dcaos(nnbasx+1:2*nnbasx))

        !  D(2k)
        udcao = 0.0d0
        call cdens1(isymc, cmo, zymc, udcao, wrk, nwrk)
        call dgefsp(nbast, udcao, dcaos(2*nnbasx+1:3*nnbasx))

        !  D(2k,1k)
        udcao = 0.0d0
        call cdens2(isymc, isymb, cmo, zymc, zymb, udcao,
     &              wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
        call dgefsp(nbast, udcao, dcaos(3*nnbasx+1:4*nnbasx))

        deallocate(udcao)

        allocate(fcaos(4*nnbasx))
#if defined(VAR_MPI)
         call pe_ifc_start_slaves(3)
#endif
         call pe_master(runtype='dynamic_response',
     &                  triang=.true.,
     &                  ndim=nbast,
     &                  nmats=4,
     &                  denmats=dcaos,
     &                  fckmats=fcaos)
        deallocate(dcaos)

        allocate(fcmo(nnorbx))
        ufcmo = 0.0d0

        i = 1
        j = nnbasx
        call uthu(2.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
        wrk(1:n2orbx) = 0.0d0
        call dsptsi(norbt, fcmo, wrk(1:n2orbx))
        call oith1(isymc, zymc, wrk(1:n2orbx), ufcmo, isyma)

        i = i + nnbasx
        j = j + nnbasx
        call uthu(fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
        wrk(1:n2orbx) = 0.0d0
        call dsptsi(norbt, fcmo, wrk(1:n2orbx))
        ufcmo = ufcmo + wrk(1:n2orbx)

        i = i + nnbasx
        j = j + nnbasx
        call uthu(2.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
        wrk(1:n2orbx) = 0.0d0
        call dsptsi(norbt, fcmo, wrk(1:n2orbx))
        call oith1(isymb, zymb, wrk(1:n2orbx), ufcmo, isyma)

        i = i + nnbasx
        j = j + nnbasx
        call uthu(fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
        wrk(1:n2orbx) = 0.0d0
        call dsptsi(norbt, fcmo, wrk(1:n2orbx))
        ufcmo = ufcmo + wrk(1:n2orbx)

        call rsp1gr(1, kzyva, idum, 0, isyma, 0, 1, etrs,
     &              wrk, idum, idum, 1.0d0, 1, udv, ufcmo, xindx,
     &              mjwop, wrk, nwrk, .true., .false., .false.)

        deallocate(fcaos, fcmo, ufcmo)

        call qexit('pe_rspqro')

      end subroutine pe_rspqro

!------------------------------------------------------------------------------
#ifdef MOD_UNRELEASED

      subroutine pe_rspcro(vecb, vecc, vecd, etrs, xindx, zymb, zymc,
     &                     zymd, udv, wrk, nwrk, kzyva, kzyvb, kzyvc,
     &                     kzyvd, isyma, isymb, isymc, isymd, cmo,mjwop)

        use polarizable_embedding, only: pe_master

        implicit none

#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "infdim.h"
#include "qrinf.h"

        integer :: kzyva, kzyvb, kzyvc, kzyvd
        integer :: isyma, isymb, isymc, isymd
        integer :: nwrk
        real*8, dimension(nwrk) :: wrk
        real*8, dimension(kzyva) :: etrs
        real*8, dimension(kzyvb) :: vecb
        real*8, dimension(kzyvc) :: vecc
        real*8, dimension(kzyvd) :: vecd
        real*8, dimension(ncmot) :: cmo
        real*8, dimension(norbt,norbt) :: zymb, zymc, zymd
        real*8, dimension(nashdi,nashdi) :: udv
        real*8, dimension(lcindx) :: xindx
        integer, dimension(2,maxwop,8) :: mjwop

        integer :: i, j, k
        integer :: idum = 1
        real*8, dimension(:), allocatable :: udcao, ufcmo
        real*8, dimension(:), allocatable :: dcaos, fcaos
        real*8, dimension(:), allocatable :: fcmo

        call qenter('pe_rspcro')

        call gtzymt(1, vecb, kzyvb, isymb, zymb, mjwop)
        call gtzymt(1, vecc, kzyvc, isymc, zymc, mjwop)
        call gtzymt(1, vecd, kzyvd, isymd, zymd, mjwop)

        allocate(udcao(n2basx))
        allocate(ufcmo(n2orbx))
        allocate(dcaos(15*nnbasx))
        dcaos = 0.0d0

        udcao = 0.0d0
        call cdens1(isymb, cmo, zymb, udcao, wrk, nwrk)
        call dgefsp(nbast, udcao, dcaos(1:nnbasx))
        udcao = 0.0d0
        call cdens1(isymc, cmo, zymc, udcao, wrk, nwrk)
        call dgefsp(nbast, udcao, dcaos(nnbasx+1:2*nnbasx))
        udcao = 0.0d0
        call cdens1(isymd, cmo, zymd, udcao, wrk, nwrk)
        call dgefsp(nbast, udcao, dcaos(2*nnbasx+1:3*nnbasx))

        udcao = 0.0d0
        call cdens2(isymb, isymc, cmo, zymb, zymc, udcao,
     &              wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
        call dgefsp(nbast, udcao, dcaos(3*nnbasx+1:4*nnbasx))
        udcao = 0.0d0
        call cdens2(isymc, isymb, cmo, zymc, zymb, udcao,
     &              wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
        call dgefsp(nbast, udcao, dcaos(4*nnbasx+1:5*nnbasx))
        udcao = 0.0d0
        call cdens2(isymb, isymd, cmo, zymb, zymd, udcao,
     &              wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
        call dgefsp(nbast, udcao, dcaos(5*nnbasx+1:6*nnbasx))
        udcao = 0.0d0
        call cdens2(isymd, isymb, cmo, zymd, zymb, udcao,
     &              wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
        call dgefsp(nbast, udcao, dcaos(6*nnbasx+1:7*nnbasx))
        udcao = 0.0d0
        call cdens2(isymc, isymd, cmo, zymc, zymd, udcao,
     &              wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
        call dgefsp(nbast, udcao, dcaos(7*nnbasx+1:8*nnbasx))
        udcao = 0.0d0
        call cdens2(isymd, isymc, cmo, zymd, zymc, udcao,
     &              wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
        call dgefsp(nbast, udcao, dcaos(8*nnbasx+1:9*nnbasx))

        udcao = 0.0d0
        call cdens3(isymb, isymc, isymd, cmo, zymb, zymc, zymd, udcao,
     &              wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
        call dgefsp(nbast, udcao, dcaos(9*nnbasx+1:10*nnbasx))
        udcao = 0.0d0
        call cdens3(isymd, isymb, isymc, cmo, zymd, zymb, zymc, udcao,
     &              wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
        call dgefsp(nbast, udcao, dcaos(10*nnbasx+1:11*nnbasx))
        udcao = 0.0d0
        call cdens3(isymc, isymd, isymb, cmo, zymc, zymd, zymb, udcao,
     &              wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
        call dgefsp(nbast, udcao, dcaos(11*nnbasx+1:12*nnbasx))
        udcao = 0.0d0
        call cdens3(isymb, isymd, isymc, cmo, zymb, zymd, zymc, udcao,
     &              wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
        call dgefsp(nbast, udcao, dcaos(12*nnbasx+1:13*nnbasx))
        udcao = 0.0d0
        call cdens3(isymc, isymb, isymd, cmo, zymc, zymb, zymd, udcao,
     &              wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
        call dgefsp(nbast, udcao, dcaos(13*nnbasx+1:14*nnbasx))
        udcao = 0.0d0
        call cdens3(isymd, isymc, isymb, cmo, zymd, zymc, zymb, udcao,
     &              wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
        call dgefsp(nbast, udcao, dcaos(14*nnbasx+1:15*nnbasx))

        deallocate(udcao)

        allocate(fcaos(15*nnbasx))
#if defined(VAR_MPI)
         call pe_ifc_start_slaves(3)
#endif
         call pe_master(runtype='dynamic_response',
     &                  triang=.true.,
     &                  ndim=nbast,
     &                  nmats=15,
     &                  denmats=dcaos,
     &                  fckmats=fcaos)
        deallocate(dcaos)

        allocate(fcmo(nnorbx))
        ufcmo = 0.0d0

        i = 1
        j = nnbasx
        call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
        wrk(1:2*n2orbx) = 0.0d0
        call dsptsi(norbt, fcmo, wrk(1:n2orbx))
        call oith1(isymc, zymc, wrk(1:n2orbx),
     &             wrk(n2orbx+1:2*n2orbx), isyma)
        call oith1(isymd, zymd, wrk(n2orbx+1:2*n2orbx), ufcmo, isyma)
        i = i + nnbasx
        j = j + nnbasx
        call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
        wrk(1:2*n2orbx) = 0.0d0
        call dsptsi(norbt, fcmo, wrk(1:n2orbx))
        call oith1(isymb, zymb, wrk(1:n2orbx),
     &             wrk(n2orbx+1:2*n2orbx), isyma)
        call oith1(isymd, zymd, wrk(n2orbx+1:2*n2orbx), ufcmo, isyma)
        i = i + nnbasx
        j = j + nnbasx
        call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
        wrk(1:2*n2orbx) = 0.0d0
        call dsptsi(norbt, fcmo, wrk(1:n2orbx))
        call oith1(isymb, zymb, wrk(1:n2orbx),
     &             wrk(n2orbx+1:2*n2orbx), isyma)
        call oith1(isymc, zymc, wrk(n2orbx+1:2*n2orbx), ufcmo, isyma)
        i = 1
        j = nnbasx
        call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
        wrk(1:2*n2orbx) = 0.0d0
        call dsptsi(norbt, fcmo, wrk(1:n2orbx))
        call oith1(isymd, zymd, wrk(1:n2orbx),
     &             wrk(n2orbx+1:2*n2orbx), isyma)
        call oith1(isymc, zymc, wrk(n2orbx+1:2*n2orbx), ufcmo, isyma)
        i = i + nnbasx
        j = j + nnbasx
        call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
        wrk(1:2*n2orbx) = 0.0d0
        call dsptsi(norbt, fcmo, wrk(1:n2orbx))
        call oith1(isymd, zymd, wrk(1:n2orbx),
     &             wrk(n2orbx+1:2*n2orbx), isyma)
        call oith1(isymb, zymb, wrk(n2orbx+1:2*n2orbx), ufcmo, isyma)
        i = i + nnbasx
        j = j + nnbasx
        call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
        wrk(1:2*n2orbx) = 0.0d0
        call dsptsi(norbt, fcmo, wrk(1:n2orbx))
        call oith1(isymc, zymc, wrk(1:n2orbx),
     &             wrk(n2orbx+1:2*n2orbx), isyma)
        call oith1(isymb, zymb, wrk(n2orbx+1:2*n2orbx), ufcmo, isyma)

        i = i + nnbasx
        j = j + nnbasx
        call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
        wrk(1:n2orbx) = 0.0d0
        call dsptsi(norbt, fcmo, wrk(1:n2orbx))
        call oith1(isymd, zymd, wrk(1:n2orbx), ufcmo, isyma)
        i = i + nnbasx
        j = j + nnbasx
        call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
        wrk(1:n2orbx) = 0.0d0
        call dsptsi(norbt, fcmo, wrk(1:n2orbx))
        call oith1(isymd, zymd, wrk(1:n2orbx), ufcmo, isyma)
        i = i + nnbasx
        j = j + nnbasx
        call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
        wrk(1:n2orbx) = 0.0d0
        call dsptsi(norbt, fcmo, wrk(1:n2orbx))
        call oith1(isymc, zymc, wrk(1:n2orbx), ufcmo, isyma)
        i = i + nnbasx
        j = j + nnbasx
        call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
        wrk(1:n2orbx) = 0.0d0
        call dsptsi(norbt, fcmo, wrk(1:n2orbx))
        call oith1(isymc, zymc, wrk(1:n2orbx), ufcmo, isyma)
        i = i + nnbasx
        j = j + nnbasx
        call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
        wrk(1:n2orbx) = 0.0d0
        call dsptsi(norbt, fcmo, wrk(1:n2orbx))
        call oith1(isymb, zymb, wrk(1:n2orbx), ufcmo, isyma)
        i = i + nnbasx
        j = j + nnbasx
        call uthu(1.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
        wrk(1:n2orbx) = 0.0d0
        call dsptsi(norbt, fcmo, wrk(1:n2orbx))
        call oith1(isymb, zymb, wrk(1:n2orbx), ufcmo, isyma)

        i = i + nnbasx
        j = j + nnbasx
        call uthu(1.0d0/3.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
        wrk(1:n2orbx) = 0.0d0
        call dsptsi(norbt, fcmo, wrk(1:n2orbx))
        ufcmo = ufcmo + wrk(1:n2orbx)
        i = i + nnbasx
        j = j + nnbasx
        call uthu(1.0d0/3.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
        wrk(1:n2orbx) = 0.0d0
        call dsptsi(norbt, fcmo, wrk(1:n2orbx))
        ufcmo = ufcmo + wrk(1:n2orbx)
        i = i + nnbasx
        j = j + nnbasx
        call uthu(1.0d0/3.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
        wrk(1:n2orbx) = 0.0d0
        call dsptsi(norbt, fcmo, wrk(1:n2orbx))
        ufcmo = ufcmo + wrk(1:n2orbx)
        i = i + nnbasx
        j = j + nnbasx
        call uthu(1.0d0/3.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
        wrk(1:n2orbx) = 0.0d0
        call dsptsi(norbt, fcmo, wrk(1:n2orbx))
        ufcmo = ufcmo + wrk(1:n2orbx)
        i = i + nnbasx
        j = j + nnbasx
        call uthu(1.0d0/3.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
        wrk(1:n2orbx) = 0.0d0
        call dsptsi(norbt, fcmo, wrk(1:n2orbx))
        ufcmo = ufcmo + wrk(1:n2orbx)
        i = i + nnbasx
        j = j + nnbasx
        call uthu(1.0d0/3.0d0*fcaos(i:j), fcmo, cmo, wrk, nbast, norbt)
        wrk(1:n2orbx) = 0.0d0
        call dsptsi(norbt, fcmo, wrk(1:n2orbx))
        ufcmo = ufcmo + wrk(1:n2orbx)

        call rsp1gr(1, kzyva, idum, 0, isyma, 0, 1, etrs,
     &              wrk, idum, idum, 1.0d0, 1, udv, ufcmo, xindx,
     &              mjwop, wrk, nwrk, .true., .false., .false.)

        deallocate(fcaos, fcmo, ufcmo)

        call qexit('pe_rspcro')

      end subroutine pe_rspcro
#endif

#endif
