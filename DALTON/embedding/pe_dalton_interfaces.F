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

      subroutine pe_ifc_infld
        use polarizable_embedding, only: pe_master
        implicit none
        call qenter('pe_ifc_infld')
        call pe_master(runtype='infld')
        call qexit('pe_ifc_infld')
      end subroutine

!-----------------------------------------------------------------------------

      subroutine pe_ifc_grad(denmats)
        use polarizable_embedding, only: pe_master
        implicit none
#include "inforb.h"
#include "mxcent.h"
#include "nuclei.h"
#include "energy.h"
        real*8, dimension(nnbasx), intent(in) :: denmats
        integer :: i, j
        call qenter('pe_ifc_grad')
#if defined(VAR_MPI)
        call pe_ifc_start_slaves(5)
#endif
        call pe_master(runtype='molecular_gradient',
     &                 triang=.true.,
     &                 ndim=nbast,
     &                 nmats=1,
     &                 denmats=denmats,
     &                 expvals=pegrad)
        call qexit('pe_ifc_grad')
      end subroutine pe_ifc_grad

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

      subroutine pe_ifc_response(denmats, fckmats, nmats)
        use polarizable_embedding, only: pe_master
        implicit none
#include "inforb.h"
        integer, intent(in) :: nmats
        real*8, dimension(nmats*nnbasx), intent(in) :: denmats
        real*8, dimension(nmats*nnbasx), intent(out) :: fckmats
        call qenter('pe_ifc_response')
#if defined(VAR_MPI)
        call pe_ifc_start_slaves(3)
#endif
        call pe_master(runtype='dynamic_response',
     &                 triang=.true.,
     &                 ndim=nbast,
     &                 nmats=nmats,
     &                 denmats=denmats,
     &                 fckmats=fckmats)
        call qexit('pe_ifc_response')
      end subroutine pe_ifc_response

!-----------------------------------------------------------------------------

      function pe_ifc_mep()
        use pe_variables, only: pe_mep
        implicit none
        integer :: pe_ifc_mep
        call qenter('pe_ifc_mep')
        if (pe_mep) then
            pe_ifc_mep = 1
        else
            pe_ifc_mep = 0
        end if
        call qexit('pe_ifc_mep')
      end function pe_ifc_mep

!-----------------------------------------------------------------------------

      function pe_ifc_mep_noqm()
        use pe_variables, only: pe_mep, mep_qmcube
        implicit none
        integer :: pe_ifc_mep_noqm
        call qenter('pe_ifc_mep_noqm')
        if (pe_mep .and. .not. mep_qmcube) then
            pe_ifc_mep_noqm = 1
        else
            pe_ifc_mep_noqm = 0
        end if
        call qexit('pe_ifc_mep_noqm')
      end function pe_ifc_mep_noqm

!-----------------------------------------------------------------------------

      subroutine pe_ifc_mepcalc(denmats)
        use polarizable_embedding, only: pe_master
        use pe_variables, only: mep_qmcube
        implicit none
#include "inforb.h"
        real*8, dimension(nnbasx), intent(in) :: denmats
        call qenter('pe_ifc_mepcalc')
#if defined(VAR_MPI)
        call pe_ifc_start_slaves(6)
#endif
        if (mep_qmcube) then
            call pe_master(runtype='mep',
     &                     triang=.true.,
     &                     ndim=nbast,
     &                     nmats=1,
     &                     denmats=denmats)
        end if
        call qexit('pe_ifc_mepcalc')
      end subroutine pe_ifc_mepcalc

!-----------------------------------------------------------------------------

      function pe_ifc_docube()
        use pe_variables, only: pe_cube
        implicit none
        integer :: pe_ifc_docube
        call qenter('pe_ifc_docube')
        if (pe_cube) then
            pe_ifc_docube = 1
        else
            pe_ifc_docube = 0
        end if
        call qexit('pe_ifc_docube')
      end function pe_ifc_docube

!-----------------------------------------------------------------------------

      subroutine pe_ifc_cube(denmats, idx_state)
        use polarizable_embedding, only: pe_master
        use pe_variables, only: pe_cube
        implicit none
#include "inforb.h"
        real*8, dimension(nnbasx), intent(in) :: denmats
        integer, intent(in) :: idx_state
        call qenter('pe_ifc_cube')
#if defined(VAR_MPI)
        call pe_ifc_start_slaves(7)
#endif
        if (pe_cube) then
          call pe_master(runtype='cube',
     &                   triang=.true.,
     &                   ndim=nbast,
     &                   nmats=1,
     &                   denmats=denmats,
     &                   idx = idx_state )
        end if
        call qexit('pe_ifc_cube')
      end subroutine pe_ifc_cube

!-----------------------------------------------------------------------------

      subroutine pe_ifc_mepcalc_noqm()
        use polarizable_embedding, only: pe_master
        implicit none
        call qenter('pe_ifc_mepcalc_noqm')
#if defined(VAR_MPI)
        call pe_ifc_start_slaves(6)
#endif
        call pe_master(runtype='mep',
     &                 triang=.true.,
     &                 ndim=0,
     &                 nmats=0)
        call qexit('pe_ifc_mepcalc_noqm')
      end subroutine pe_ifc_mepcalc_noqm

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
        call pe_master(runtype='magnetic_gradient',
     &                 triang=.true.,
     &                 ndim=nbast,
     &                 nmats=0,
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
        else if (runtype == 5) then
            call pe_slave('molecular_gradient')
        else if (runtype == 6) then
            call pe_slave('mep')
        else if (runtype == 7) then
            call pe_slave('cube')
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
            call pe_ifc_response(fdtvaos, fxcaos, ncsim)
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
            call pe_ifc_response(bodtaos, fxoaos, nosim)
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
            call pe_ifc_response(fdtvaos, fxcaos, ncsim)
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
!        real*8, dimension(:), allocatable :: dcaotr, dvaotr
        real*8, dimension(:), allocatable :: daos, fckaos
        real*8, dimension(:), allocatable :: daotrs
!        real*8, dimension(:), allocatable :: daotrs, fckaotrs
        real*8, dimension(:), allocatable :: evec
        real*8, dimension(:,:), allocatable :: ubovecs, evecs, eacs
!        real*8, dimension(:,:), allocatable :: evectrs, eactrs
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
!        if (trplet) then
!            allocate(dcaotr(n2basx), dvaotr(n2basx),
!     &               daotrs(nosim*nnbasx))
!        end if
        ! Calculate Fxo = <0|Fe(k)|0>Fe
        do i = 1, nosim
            j = (i - 1) * nnbasx + 1
            call deq27(cmo, ubovecs(:,i), udv, dcao, dvao,
     &                 wrk, nwrk)
!            edh: Does deq27 assume symmetric matrix? is this
!            reasonable?
!            if (trplet) then
!                call deq27(cmo, ubovecs(:,i), udvtr, dcaotr, dvaotr,
!    &                     wrk, nwrk)
!            end if
            if (nasht > 0) then
                dcao = dcao + dvao
!                if (trplet) then
!                    dcaotr = dcaotr + dvaotr
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
        call pe_ifc_response(daos, fckaos, nosim)
        deallocate(daos)
!        if (trplet) then
!            allocate(fckaotrs(nosim*nnbasx))
!            call pe_master(runtype='response', denmats=daotrs,
!     &                     fckmats=fckaotrs, nmats=nosim)
!            deallocate(daotrs)
!        end if

        allocate(evec(nnorbx))
        allocate(evecs(n2orbx,nosim))
        evecs = 0.0d0
        if (.not. tdhf) then
            allocate(eacs(n2ashx,nosim))
            allocate(txyoacs(nosim))
            eacs = 0.0d0
            txyoacs = 0.0d0
        end if
!        if (trplet) then
!            allocate(evectrs(n2orbx,nosim), eactrs(n2ashx,nosim))
!        end if

        do i = 1, nosim
            j = (i - 1) * nnbasx + 1
            call uthu(fckaos(j), evec, cmo, wrk, nbast, norbt)
            call dsptsi(norbt, evec, evecs(:,i))
!            if (trplet) then
!                uthu(fckaotrs(j), evectrs, cmo, wrk, nbast, norbt)
!                call dsptsi(norbt, evectrs, evecstrs(:,i))
!            end if

            ! Fyo = V(k) - <0|F|0>Fe(k)
            if (.not. tdhf) then
                call onexh1(ubovecs(:,i), fupemo, evecs(:,i))
                call getacq(evecs(:,i), eacs(:,i))
                txyo = slvqlm(udv, eacs(:,i), evecs(:,i), txyoacs(i))
!                if (trplet) then
!                    call getacq(evecstrs, eacstrs)
!                    txyot  = slvqlm(udvtr, eacstrs, evecstrs, fyoat(i))
!                end if
            end if
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
        if (.not. tdhf) then
            deallocate(fupemo)
            call slvsc(0, nosim, n2ashx, dummy, cref, sovecs, eacs,
     &                 dummy, txyoacs, dummy, cindx, wrk, nwrk)
            deallocate(eacs)
            deallocate(txyoacs)
        end if

!        if (trplet) then
!            call slvsor(.true., .false., nosim, udvtr, sovecs, evectrs)
!            call slvsor(.true., .true., nosim, udv, sovecs, evecs)
!        else
            call slvsor(.true., .true., nosim, udv, sovecs, evecs)
!        end if

        deallocate(evecs)
!        if (trplet) then
!            deallocate(evectrs, eactrs)
!        end if

!        write(lupri,*)' Linear transformed orbital vector'
!        call output(sovecs,1,kzyvar,1,nosim,kzyvar,nosim,1,lupri)

        call qexit('pe_rsplno')

      end subroutine pe_rsplno

!------------------------------------------------------------------------------


      subroutine pe_rspmcqr(vecb, vecc, etrs, xindx, zymb, zymc,
     &                      den1, udv, wrk, lfree, kzyva, kzyvb, kzyvc,
     &                      isyma, isymb, isymc, cmo, mjwop)

         use pe_variables, only: pe_polar
         use polarizable_embedding, only: pe_master

         implicit none

#include "inforb.h"
#include "infvar.h"
#include "infdim.h"
#include "qrinf.h"
#include "priunit.h"
#include "dummy.h"
#include "inftap.h"
#include "infrsp.h"
#include "wrkrsp.h"

         integer :: kzyva, kzyvb, kzyvc
         integer :: isyma, isymb, isymc, isymbc
         integer :: lfree
         integer :: ilsym, irsym, ncl, ncr, kzvarl, kzvarr
         integer :: isymdn, isymst
         integer :: kcref, nzyvec, nzcvec
         integer :: iprone, nzconf, nzvar
         integer :: n2ash

         real*8 :: ovlap
         real*8 :: fact
         real*8 :: ddot

         real*8, dimension(*) :: wrk
         real*8, dimension(1) :: tmpwrk
         real*8, dimension(*) :: cmo, xindx

         real*8, dimension(kzyva) :: etrs
         real*8, dimension(kzyvb) :: vecb
         real*8, dimension(kzyvc) :: vecc

         real*8, dimension(norbt,norbt) :: zymb, zymc
         real*8, dimension(nashdi,nashdi) :: udv, den1
         real*8, dimension(nnashx) :: dv

         integer, dimension(2,maxwop,8) :: mjwop

         real*8, dimension(:), allocatable :: fpe
         real*8, dimension(:), allocatable :: cref
         real*8, dimension(:), allocatable :: dcaos, fcaos
         real*8, dimension(:), allocatable :: udtv, udtvao
         real*8, dimension(:), allocatable :: dvaao, dvbao, dvatr
         real*8, dimension(:), allocatable :: udcao, udvao
         real*8, dimension(:), allocatable :: udcmo, udvmo
         real*8, dimension(:), allocatable :: fcmo

         real*8, dimension(:,:), allocatable :: dva, dvb
         real*8, dimension(:,:), allocatable :: fupe
         real*8, dimension(:,:), allocatable :: fxpeb,fxpec, fx2pe
         real*8, dimension(:,:), allocatable :: fxo1k, fxc1s
         real*8, dimension(:,:), allocatable :: fxo2k, fxc2s
         real*8, dimension(:,:), allocatable :: fcas2_1, fcas2_2
         real*8, dimension(:,:), allocatable :: fcas3_1, fcas3_2
         real*8, dimension(:,:), allocatable :: fxo, fxo1k2k, fxo2k1k
         real*8, dimension(:,:), allocatable :: fxc1s2s
         real*8, dimension(:,:), allocatable :: fxc1s2k, fxc1k2s

         logical :: lexist, lopen, lcon, lorb
         logical :: fndlab

         lopen = .false.

         call qenter('pe_rspmcqr')

         call gtzymt(1, vecb, kzyvb, isymb, zymb, mjwop)
         call gtzymt(1, vecc, kzyvc, isymc, zymc, mjwop)

!         CALL GTZYMT(NSIM,VEC1,KZYV1,ISYMV1,ZYM1,MJWOP)
!         CALL GTZYMT(NSIM,VEC2,KZYV2,ISYMV2,ZYM2,MJWOP)

         !-----------------------------------------------------------
         ! Get Fg = Vmul - R*<0|F|>Fe from file
         !-----------------------------------------------------------
         if (.not. tdhf) then
            allocate(fpe(nnorbx))
            if (lusifc <= 0) then
                call gpopen(lusifc, 'SIRIFC', 'OLD', ' ', 'UNFORMATTED',
     &       idummy, .false.)
                    lopen = .true.
            end if
            rewind(lusifc)
            call mollab('PEFMAT  ', lusifc, lupri)
            call readt(lusifc, nnorbx, fpe)
            if (lopen) call gpclose(lusifc, 'KEEP')
            allocate(fupe(norbt,norbt))
            call dsptsi(norbt, fpe, fupe)
            deallocate(fpe)
         end if
         !-----------------------------------------------------------
         ! Density Factory ...
         !-----------------------------------------------------------

         allocate(cref(mzconf(1)))
         call getref(cref, mzconf(1))

         if (pe_polar) then
            allocate(udcao(n2basx))
            allocate(udvao(n2basx))
            if (.not. tdhf) then
               allocate(dcaos(10*nnbasx))
               else
               allocate(dcaos(4*nnbasx))
            end if
            dcaos = 0.0d0

            !  DTX = D_pq(k1) = <0|[k1,Epq]|0>
            allocate(udcmo(n2orbx),udvmo(n2orbx))
            udcmo = 0.0d0
            udvmo = 0.0d0
            call deq27mo(isymb, zymb, udv, udcmo, udvmo,
     &                   wrk, lfree)
            if (nasht > 0) then
               udcmo = udcmo + udvmo
            end if
            udcao = 0.0d0
            call motoao(udcmo,udcao,cmo,isymb,wrk,lfree)
            call dgefsp(nbast, udcao, dcaos(1:nnbasx))
            ! needed to fit with HF code
            dcaos(1:nnbasx) = 0.5d0*dcaos(1:nnbasx)

            !  DT2X = D_pq(k2,k1) = <0|[k2,[k1,Epq]|0>
            udvmo = 0.0d0
            call oitd1(isymc,zymc,udcmo,udvmo,isymb)
            ! DT2X in udvmo (re-used to save memory)
            udcao = 0.0d0
            isymbc = muld2h(isymb,isymc)
            call motoao(udvmo,udcao,cmo,isymbc,wrk,lfree)
            call dgefsp(nbast, udcao, dcaos(nnbasx+1:2*nnbasx))
            ! needed to fit with HF code
            dcaos(nnbasx+1:2*nnbasx) = 0.5d0*dcaos(nnbasx+1:2*nnbasx)

            !  DTX = D_pq(k2) = <0|[k2,Epq]|0>
            udcmo = 0.0d0
            udvmo = 0.0d0
            call deq27mo(isymc, zymc, udv, udcmo, udvmo,
     &                   wrk, lfree)
            if (nasht > 0) then
               udcmo = udcmo + udvmo
            end if
            udcao = 0.0d0
            call motoao(udcmo,udcao,cmo,isymb,wrk,lfree)
            call dgefsp(nbast, udcao, dcaos(2*nnbasx+1:3*nnbasx))
            ! needed to fit with HF code
            dcaos(2*nnbasx+1:3*nnbasx) =
     &      0.5d0*dcaos(2*nnbasx+1:3*nnbasx)

            !  DT2X = D_pq(k1,k2) = <0|[k1,[k2,Epq]|0>
            udvmo = 0.0d0
            call oitd1(isymb,zymb,udcmo,udvmo,isymc)
            ! DT2X in udvmo (re-used to save memory)
            udcao = 0.0d0
            isymbc = muld2h(isymc,isymb)
            call motoao(udvmo,udcao,cmo,isymbc,wrk,lfree)
            call dgefsp(nbast, udcao, dcaos(3*nnbasx+1:4*nnbasx))
            ! needed to fit with HF code
            dcaos(3*nnbasx+1:4*nnbasx) =
     &      0.5d0*dcaos(3*nnbasx+1:4*nnbasx)
            deallocate(udcmo,udvmo)

            if (tdhf) then
            write(lupri,*) 'PE-DFT or HF QR detected: Skipping CI dens.'
            end if
            if (.not. tdhf ) then
            write(lupri,*) 'PE-MCSCF QR detected: Constructing CI dens.'

               ! Construct the density matrix <02L|..|0> + <0|..|02R>
               ilsym  = irefsy
               irsym  = muld2h(irefsy,isymc)
               ncl    = mzconf(1)
               ncr    = mzconf(isymc)
               kzvarl = mzconf(1)
               kzvarr = mzyvar(isymc)

               den1 = 0.0d0 ! edh: This is equal to udtv later...
               allocate(udtv(n2ashx), udtvao(n2basx))

               udtv = 0.0d0
               udtvao = 0.0d0
               call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,
     &                     cref, vecc, ovlap, udtv, dummy, 0 ,0, .true.,
     &                    .true., xindx, wrk, 1, lfree, .true.)
               call fckden2(.false.,.true., dummy, udtvao, cmo,
     &                   udtv, wrk, lfree)
               call dgefsp(nbast, udtvao, dcaos(4*nnbasx+1:5*nnbasx))
               dcaos(4*nnbasx+1:5*nnbasx) =
     &         1.0d0*dcaos(4*nnbasx+1:5*nnbasx)

               ! Construct the density matrix <01L|..|0> + <0|..|01R>
               ilsym  = irefsy
               irsym  = muld2h(irefsy,isymb)
               ncl    = mzconf(1)
               ncr    = mzconf(isymb)
               kzvarl = mzconf(1)
               kzvarr = mzyvar(isymb)

               udtv = 0.0d0
               udtvao = 0.0d0
               call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,
     &                     cref, vecb, ovlap, udtv, dummy, 0 ,0, .true.,
     &                     .true., xindx, wrk, 1, lfree, .true.)
               call fckden2(.false.,.true., dummy, udtvao, cmo,
     &                      udtv, wrk, lfree)
               call dgefsp(nbast, udtvao, dcaos(5*nnbasx+1:6*nnbasx))
                dcaos(5*nnbasx+1:6*nnbasx) =
     &         1.0d0*dcaos(5*nnbasx+1:6*nnbasx)

               if (mzconf(isymb) .gt. 0 .and. mzconf(isymc) .gt. 0) then

                  ! Construct <01L|..|02R> + <02L|..|01R> density
                  ilsym  = muld2h(irefsy,isymb)
                  irsym  = muld2h(irefsy,isymc)
                  ncl    = mzconf(isymb)
                  ncr    = mzconf(isymc)
                  kzvarl = mzyvar(isymb)
                  kzvarr = mzyvar(isymc)
                  isymdn = muld2h(ilsym,irsym)

                  udtv = 0.0d0
                  udtvao = 0.0d0
                  call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,
     &                        vecb, vecc, ovlap, udtv, dummy, 0 ,0,
     &                       .true., .true., xindx, wrk, 1, lfree,
     &                       .false.)
                  call fckden2(.false.,.true., dummy, udtvao, cmo,
     &                       udtv, wrk, lfree)
                  call dgefsp(nbast, udtvao, dcaos(6*nnbasx+1:7*nnbasx))
               end if
               dcaos(6*nnbasx+1:7*nnbasx) =
     &         1.0d0*dcaos(6*nnbasx+1:7*nnbasx)

               ! D_pq = <0|Epq|0>
               udcao = 0.0d0
               udvao = 0.0d0
               call dgefsp(nasht, udv, dv)
               call fckden((nisht>0), (nasht>0), udcao, udvao,
     &                      cmo, dv, wrk, lfree)
               if (nisht==0) udcao = 0.0d0
               udcao = udcao + udvao
               call dgefsp(nbast, udcao, dcaos(7*nnbasx+1:8*nnbasx))
               dcaos(7*nnbasx+1:8*nnbasx) =
     &         0.5d0*dcaos(7*nnbasx+1:8*nnbasx)

               ! D_pq(S1,k2) = <01L|[k2,Epq]|0> + <0|[k2,Epq]|01R>
               allocate(dva(norbt,nasht), dvb(norbt,nasht))
               allocate(dvaao(n2basx), dvbao(n2basx), dvatr(n2basx))
               dva    = 0.0d0
               dvb    = 0.0d0
               dvaao  = 0.0d0
               dvbao  = 0.0d0
               dvatr  = 0.0d0
               udtvao = 0.0d0
               call rsptr1(1, udv, zymb, dva, dvb)
               call fckden2(.false.,.true., dummy, dvaao, cmo,
     &                      dva, wrk, lfree)
               call fckden2(.false.,.true., dummy, dvbao, cmo,
     &                      dvb, wrk, lfree)
               call mtrsp(nbast, nbast, dvaao, nbast, dvatr, nbast)
               udtvao = dvbao - dvatr
               call dgefsp(nbast, udtvao, dcaos(8*nnbasx+1:9*nnbasx))
               dcaos(8*nnbasx+1:9*nnbasx) =
     &         2.0d0*dcaos(8*nnbasx+1:9*nnbasx)

               ! D_pq(k1,S2) = <02L|[k1,Epq]|0> + <0|[k1,Epq]|02R>
               dva    = 0.0d0
               dvb    = 0.0d0
               dvaao  = 0.0d0
               dvbao  = 0.0d0
               dvatr  = 0.0d0
               udtvao = 0.0d0
               call rsptr1(1, udv, zymc, dva, dvb)
               call fckden2(.false.,.true., dummy, dvaao, cmo,
     &                      dva, wrk, lfree)
               call fckden2(.false.,.true., dummy, dvbao, cmo,
     &                      dvb, wrk, lfree)
               call mtrsp(nbast, nbast, dvaao, nbast, dvatr, nbast)
               udtvao = dvbao - dvatr
               call dgefsp(nbast, udtvao, dcaos(9*nnbasx+1:10*nnbasx))
               dcaos(9*nnbasx+1:10*nnbasx) =
     &         2.0d0*dcaos(9*nnbasx+1:10*nnbasx)

               deallocate(dva, dvb, dvaao, dvbao, dvatr)
               deallocate(udtv,udtvao)

            end if
            deallocate(udcao, udvao)
            !-----------------------------------------------------------
            ! Calculate PE response operators in AO basis
            !-----------------------------------------------------------

            if (.not. tdhf) then
               allocate(fcaos(10*nnbasx))
               fcaos = 0.0d0
#if defined(VAR_MPI)
               call pe_ifc_start_slaves(3)
#endif
               call pe_master(runtype='dynamic_response',
     &                        triang=.true.,
     &                        ndim=nbast,
     &                        nmats=10,
     &                        denmats=dcaos,
     &                        fckmats=fcaos)
            else
               allocate(fcaos(4*nnbasx))
               fcaos = 0.0d0
#if defined(VAR_MPI)
               call pe_ifc_start_slaves(3)
#endif
               call pe_master(runtype='dynamic_response',
     &                        triang=.true.,
     &                        ndim=nbast,
     &                        nmats=4,
     &                        denmats=dcaos,
     &                        fckmats=fcaos)
            end if
            deallocate(dcaos)

         end if

         if ( .not. tdhf ) then
            !-----------------------------------------------------------
            !case 1
            !-----------------------------------------------------------
            if ( mzconf(isymb) .eq. 0 .or. mzconf(isymc) .eq. 0 ) return

            !/   <01L| [qj,TB] |02R>  + <02L| [qj,TB] |01R>  \
            !|                       0                       |
            !|   <01L| [qj+,TB] |02R> + <02L| [qj+,TB] |01R> |
            !\                       0                       /

            ! ionstruct <01L|..|02R> + <02L|..|01R> density
            ilsym  = muld2h(irefsy,isymb)
            irsym  = muld2h(irefsy,isymc)
            ncl    = mzconf(isymb)
            ncr    = mzconf(isymc)
            kzvarl = mzyvar(isymb)
            kzvarr = mzyvar(isymc)

            den1 = 0.0d0
            call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,
     &                  vecb, vecc, ovlap, den1, dummy, 0, 0, .true.,
     &                  .true., xindx, wrk, 1, lfree, .false.)

            ! Make the gradient
            isymdn = muld2h(ilsym,irsym)

            if ( mzwopt(isyma) .gt. 0 ) then
               call orbsx(1, isyma, kzyva, etrs, fupe, ovlap,
     &                    isymdn, den1, mjwop, 1, lfree)
           end if
        end if
        !-----------------------------------------------------------
        !case 2
        !-----------------------------------------------------------

        if (pe_polar) then

           allocate(fcmo(nnorbx))
           allocate(fxo1k(norbt,norbt))
           fxo1k = 0.0d0
           if (.not. tdhf) then
              allocate(fxc1s(norbt,norbt))
              fxc1s = 0.0d0
           end if

           ! Fxo = R*<0|Fe(1k)|0>Fe
           fcmo = 0.0d0
           call uthu(2.0d0*fcaos(1:nnbasx), fcmo, cmo,
     &               wrk, nbast, norbt)
           call dsptsi(norbt, fcmo, fxo1k)

           if (.not. tdhf) then
              ! Fxc(1S) = ( R*<01lE|0>+<0|E01R> )Fe
              ! edh: Should it be 1.0d0 or 2.0d0 ???
              fcmo = 0.0d0
              call uthu(1.0d0*fcaos(4*nnbasx+1:5*nnbasx), fcmo, cmo,
     &                  wrk, nbast, norbt)
              call dsptsi(norbt, fcmo, fxc1s)
           end if

           ! fcas2_1 = Fa[1](1k)
           allocate(fcas2_1(norbt,norbt))
           fcas2_1 = 0.0d0
           if (.not. tdhf) then
              fcas2_1 = fxo1k + fxc1s
              deallocate(fxc1s)
           else
              fcas2_1 = fxo1k
           end if
           deallocate(fxo1k)

           if (.not. tdhf) then
              if (mzconf(isymc) .le. 0) return

             !/   0    \
             !| Sj(2)  | * <0| Fa[1](1k) |0>
             !|   0    |
             !\ Sj(2)* /

              if (isyma .eq. isymc) then
                  ovlap = 1.0d0
                  call melone(fcas2_1, 1, udv, ovlap, fact,
     &                        200,'fact for Fxo(1k) + Fxc(1S) ')
                  nzconf = mzconf(isyma)
                  nzvar  = mzvar(isyma)
                  call daxpy(nzconf, fact, vecc, 1, etrs, 1)
                  call daxpy(nzconf,fact,
     &                         vecc(nzvar+1), 1, etrs(nzvar+1), 1)
              end if
           end if

           allocate(fxo2k(norbt,norbt))
           fxo2k = 0.0d0
           if (.not. tdhf) then
              allocate(fxc2s(norbt,norbt))
              fxc2s = 0.0d0
           end if

           ! Fxo(2k) = R*<0|[2k,Epq]|0>Fe
           fcmo = 0.0d0
           call uthu(2.0d0*fcaos(2*nnbasx+1:3*nnbasx),
     &                fcmo, cmo, wrk, nbast, norbt)
           call dsptsi(norbt, fcmo, fxo2k)

           if (.not. tdhf) then
              ! Fxc(1S) = ( R*<01lE|0>+<0|E01R> )Fe
              ! edh: Should it be 1.0d0 or 2.0d0 ???
              fcmo = 0.0d0
              call uthu(1.0d0*fcaos(5*nnbasx+1:6*nnbasx), fcmo, cmo,
     &                   wrk, nbast, norbt)
               call dsptsi(norbt, fcmo, fxc2s)
           end if

           ! fcas2_2 = Fa[1](2k)
           allocate(fcas2_2(norbt,norbt))
           fcas2_2 = 0.0d0
           if (.not. tdhf) then
              fcas2_2 = fxo2k + fxc2s
              deallocate(fxc2s)
           else
              fcas2_2 = fxo2k
           end if
           deallocate(fxo2k)

           if (.not. tdhf) then
              if (mzconf(isymb) .le. 0) return

             !/   0    \
             !| Sj(1)  | * <0| Fa[1](2k) |0>
             !|   0    |
             !\ Sj(1)* /

              if (isyma .eq. isymb) then
                  ovlap = 1.0d0
                  call melone(fcas2_2, 1, udv, ovlap, fact,
     &                        200,'fact for Fxo(1k) + Fxc(1S) ')
                    nzconf = mzconf(isyma)
                    nzvar  = mzvar(isyma)
                    call daxpy(nzconf, fact, vecb, 1, etrs, 1)
                    call daxpy(nzconf,fact,
     &                         vecb(nzvar+1), 1, etrs(nzvar+1), 1)
              end if
           end if
        end if
!       !-----------------------------------------------------------
!       !case 3
!       !-----------------------------------------------------------
        allocate(fxpeb(norbt,norbt))
        fxpeb = 0.0d0
        call oith1(isymb,zymb,fupe,fxpeb,1)

        ! fcas3_1 = Fg(k1) + F[a1]
        allocate(fcas3_1(norbt,norbt))
        if (pe_polar) then
            fcas3_1 = fxpeb + fcas2_1
        else
            fcas3_1 = fxpeb
        end if

        if (.not. tdhf) then

           if (mzconf(isymc) .le. 0) return

          !/   <0| [qj,TD1] |02R>  + <02L| [qj,TD1] |0>  \
          !|   <j| TD1 |02R>                             |
          !|   <0| [qj+,TD1] |02R> + <02L| [qj+,TD1] |0> |
          !\  -<02L| TD1 |j>                             /

          ! 1a. Construct the density matrix <02L|..|0> + <0|..|02R>
           ilsym  = irefsy
           irsym  = muld2h(irefsy,isymc)
           ncl    = mzconf(1)
           ncr    = mzconf(isymc)
           kzvarl = mzconf(1)
           kzvarr = mzyvar(isymc)

           den1 = 0.0d0
           call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,
     &                 cref, vecc, ovlap, den1, dummy, 0 ,0, .true.,
     &                 .true., xindx, wrk, 1, lfree, .true.)

!          1b. Make the gradient
           isymdn = muld2h(ilsym,irsym)
           isymst = muld2h(isyma,irefsy)
           if ( isymst .eq. irefsy ) then
              lcon = ( mzconf(isyma) .gt. 1 )
           else
              lcon = ( mzconf(isyma) .gt. 0 )
           end if
           lorb    = ( mzwopt(isyma) .gt. 0 )
           nzyvec = mzyvar(isymc)
           nzcvec = mzconf(isymc)

           call rsp1gr(1, kzyva, idummy, 0 , isyma, 0, isymc, etrs,
     &               vecc, nzyvec, nzcvec, ovlap, isymdn, den1, fcas3_1,
     &               xindx, mjwop, wrk(1), lfree, lorb, lcon, .false.)
        end if
        deallocate(fcas3_1)

        allocate(fxpec(norbt,norbt))
        fxpec = 0.0d0
        call oith1(isymc,zymc,fupe,fxpec,1)

        ! fcas3_2 = Fg(2k) + Fa[1]
        allocate(fcas3_2(norbt,norbt))
        if (pe_polar) then
            fcas3_2 = fxpec + fcas2_2
        else
            fcas3_2 = fxpec
        end if

        if (.not. tdhf) then

           if (mzconf(isymb) .le. 0) return

          !/   <0| [qj,TD2] |01R>  + <01L| [qj,TD2] |0>  \
          !|   <j| TD2 |01R>                             |
          !|   <0| [qj+,TD2] |01R> + <01L| [qj+,TD2] |0> |
          !\  -<01L| TD2 |j>                             /

          ! 2a. Construct the density matrix <01L|..|0> + <0|..|01R>
           ilsym  = irefsy
           irsym  = muld2h(irefsy,isymb)
           ncl    = mzconf(1)
           ncr    = mzconf(isymb)
           kzvarl = mzconf(1)
           kzvarr = mzyvar(isymb)

           den1 = 0.0d0
           call rspgdm(1, ilsym, irsym, ncl, ncr, kzvarl, kzvarr,
     &               cref, vecb, ovlap, den1, dummy, 0 ,0, .true.,
     &               .true., xindx, wrk, 1, lfree, .true.)

          ! 2b. Make the gradient
           isymdn = muld2h(ilsym,irsym)
           isymst = muld2h(isyma,irefsy)
           if ( isymst .eq. irefsy ) then
              lcon = ( mzconf(isyma) .gt. 1 )
           else
              lcon = ( mzconf(isyma) .gt. 0 )
           end if
           lorb    = ( mzwopt(isyma) .gt. 0 )
           nzyvec = mzyvar(isymb)
           nzcvec = mzconf(isymb)

           call rsp1gr(1, kzyva, idummy, 0 , isyma, 0, isymb, etrs,
     &               vecb, nzyvec, nzcvec, ovlap, isymdn, den1, fcas3_2,
     &               xindx, mjwop, wrk(1), lfree, lorb, lcon, .false.)
        end if
        deallocate(fcas3_2)
!       !-----------------------------------------------------------
!       !case 4
!       !-----------------------------------------------------------

!       !fx2pe = 0.5*Fg(1k,2k) + 0.5*Fg(2k,1k) + ...
        allocate(fx2pe(norbt,norbt))
        fx2pe = 0.0d0
        if (.not. tdhf) then
           call oith1(isymc, zymc, fxpeb, fx2pe, isymb)
           call oith1(isymb, zymb, fxpec, fx2pe, isymc)
           deallocate(fxpeb,fxpec)
        end if
        fx2pe = 0.5d0 * fx2pe

        ! ... + fcas2_1(2k) + fcas2_2(k1)
        if (pe_polar) then
            call oith1(isymc,zymc,fcas2_1, fx2pe,  isymb)
            call oith1(isymb,zymb,fcas2_2, fx2pe, isymc)
            deallocate(fcas2_1,fcas2_2)

            if (.not. tdhf ) then
            ! + ( S(1)S*(2) + S(2)S*(1) ) * fxo  + ...
               if ((isymb. eq. isymc) .and. (mzconf(isymb) .gt. 0)) then
                  allocate(fxo(norbt,norbt))
                  fcmo = 0.0d0
                  fxo = 0.0d0
                 call uthu(0.25d0*fcaos(7*nnbasx+1:8*nnbasx), fcmo, cmo,
     &                      wrk, nbast, norbt)
                  call dsptsi(norbt, fcmo, fxo)
                  nzconf = mzconf(isymb)
                  nzvar  = mzvar(isymb)
                  fact   = ddot(nzconf, vecb, 1, vecc(nzvar+1), 1) +
     &                     ddot(nzconf, vecc, 1, vecb(nzvar+1), 1)
                           call daxpy(n2orbx, fact, fxo, 1, fx2pe, 1)
                  deallocate(fxo)
               end if
            end if

            allocate(fxo1k2k(norbt,norbt), fxo2k1k(norbt,norbt))
            fxo1k2k = 0.0d0
            fxo2k1k = 0.0d0
            if (.not. tdhf) then
               allocate(fxc1s2s(norbt,norbt))
               allocate(fxc1k2s(norbt,norbt), fxc1s2k(norbt,norbt))
               fxc1s2s = 0.0d0
               fxc1s2k = 0.0d0
               fxc1k2s = 0.0d0
            end if

            ! ... + fxo(1k,2k) + fxo(2k,1k)
            if (mzwopt(isymb).gt.0 .and. mzwopt(isymc).gt.0) then
               !fxo(1k,2k) = <0|Fe(1k,2k)|0>Fe
               fcmo = 0.0d0
               call uthu(1.0d0*fcaos(nnbasx+1:2*nnbasx),
     &                   fcmo, cmo, wrk, nbast, norbt)
               call dsptsi(norbt, fcmo, fxo1k2k)

               !fxo(2k,1k) = R*<0|Fe(2k,1k)|0>Fe
               fcmo    = 0.0d0
               call uthu(1.0d0*fcaos(3*nnbasx+1:4*nnbasx),
     &                   fcmo, cmo, wrk, nbast, norbt)
               call dsptsi(norbt, fcmo, fxo2k1k)
            end if

            ! + ... fxc(1s2s) + fxc(2s1s) + 2fxc(1s2k) + 2fxc(1k2S)
            if (.not. tdhf) then
               ! fxc(1s2s) = R* ( <01L|..|02R> + <02L|..|01R> )Fe
               fcmo    = 0.0d0
               call uthu(0.25d0*fcaos(6*nnbasx+1:7*nnbasx), fcmo, cmo,
     &                   wrk, nbast, norbt)
               call dsptsi(norbt, fcmo, fxc1s2s)

               ! + ... fxc(1s2k) + fxc(1k2S)
               ! fxc(1s2k) = ( <01L|[k2,Epq]|0> + <0|[k2,Epq]|01R> )Fe
               call uthu(2.0d0*fcaos(8*nnbasx+1:9*nnbasx), fcmo, cmo,
     &                   wrk, nbast, norbt)
               call dsptsi(norbt, fcmo, fxc1s2k)

               ! fxc(1k2s) = ( <02L|[k1,Epq]|0> + <0|[k1,Epq]|02R> )Fe
               call uthu(2.0d0*fcaos(9*nnbasx+1:10*nnbasx), fcmo, cmo,
     &                   wrk, nbast, norbt)
               call dsptsi(norbt, fcmo, fxc1k2s)
            end if
            deallocate(fcmo)

            if (.not. tdhf) then
               fx2pe = fx2pe + fxo1k2k + fxo2k1k
     &               + fxc1s2s + fxc1s2k + fxc1k2s
               deallocate(fxc1s2s, fxc1s2k, fxc1k2s)
            else
               fx2pe = fx2pe + fxo1k2k + fxo2k1k
            end if
            deallocate(fxo1k2k, fxo2k1k)

       end if

       !/ <0| [qj ,TE] |0> \
       !| <j| TE |0>       |
       !| <0| [qj+,TE] |0> |
       !\ -<0| TE |j>      /

        isymdn = 1
        ovlap  = 1.0d0
        isymst = muld2h(isyma, irefsy)
        if ( isymst .eq. irefsy ) then
           lcon = ( mzconf(isyma) .gt. 1 )
        else
           lcon = ( mzconf(isyma) .gt. 0 )
        end if
        lorb   = ( mzwopt(isyma) .gt. 0 )
        nzyvec = mzconf(1)
        nzcvec = mzconf(1)

        call rsp1gr(1 ,kzyva, idummy,0, isyma, 0, irefsy, etrs,
     &              cref, nzyvec, nzcvec, ovlap, isymdn, udv, fx2pe,
     &              xindx, mjwop, wrk(1), lfree, lorb, lcon, .true.)
        deallocate(cref)
        deallocate(fupe)

        call qexit('pe_rspmcqr')

        end subroutine pe_rspmcqr

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
        call pe_ifc_response(dcaos, fcaos, 4)
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
        call pe_ifc_response(dcaos, fcaos, 15)
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

!------------------------------------------------------------------------------

      subroutine dalton_cavity(all_centers, all_charges, ncents,
     &                         wrk, nwrk)

!        implicit none

#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "mxcent.h"
#include "pcmdef.h"
#include "pcm.h"
#include "pcmlog.h"
#include "nuclei.h"

        integer :: ncents, nwrk
        real*8, dimension(ncents) :: all_charges
        real*8, dimension(3,ncents) :: all_centers
        real*8, dimension(nwrk) :: wrk

        integer :: i
        integer :: ntess, lusurf, lucav
        integer :: bck_nucind, bck_nucdep, bck_natoms
        integer, dimension(:), allocatable :: isphere
        real*8, dimension(:), allocatable :: areas, rvdw
        real*8, dimension(:,:), allocatable :: tess_coords

        call qenter('dalton_cavity')

C       van der Waals radii taken from "The Elements", 2nd edition,
C       John Emsley, Clarendon Press, Oxford, 1991.  Unknown values
C       are simply set to D0, rather than trying to guess them.
C
C       A.Bondi, J.Phys.Chem. 68: 441-451(1964) gives alternate
C       values, and a few transition metals.
        allocate(rvdw(99))
        rvdw = (/ 1.20d0, 1.22d0, 0.00d0, 0.00d0, 2.08d0, 1.85d0,
     &            1.54d0, 1.40d0, 1.35d0, 1.60d0, 2.31d0, 0.00d0,
     &            2.05d0, 2.00d0, 1.90d0, 1.85d0, 1.81d0, 1.91d0,
     &            2.31d0, 13*0.0d0, 2.00d0, 2.00d0, 1.95d0, 1.98d0,
     &            2.44d0, 13*0.0d0, 2.20d0, 2.20d0, 2.15d0, 0.00d0,
     &            2.62d0, 27*0.0d0, 2.40d0, 16*0.0d0 /)
C     override the above table with U. Pisa's experience
C     as to what works best for singly bonded C,N,O
        rvdw(6) = 1.70d0
        rvdw(7) = 1.60d0
        rvdw(8) = 1.50d0

        if (ncents > mxsp) then
            stop 'ERROR: ncents > mxsp'
        end if

        nesfp = ncents
        nesf = nesfp

        do i = 1, nesfp
            xe(i) = all_centers(1,i)
            ye(i) = all_centers(2,i)
            ze(i) = all_centers(3,i)
            if (nint(all_charges(i)) > 99) then
                stop 'ERROR: charge is too big'
            end if
            rin(i) = rvdw(nint(all_charges(i)))
            if (rin(i) == 0.0d0) then
                stop 'ERROR: no vdw radius defined'
            end if
        end do

        bck_nucind = nucind
        bck_nucdep = nucdep
        bck_natoms = natoms

        nucind = ncents
        nucdep = ncents
        natoms = ncents

        omega = 40.0d0
        fro = 0.7d0
        ret = 0.2d0
        nrwcav = 2
        areats = 0.3d0
        oldcen = .false.
        iprpcm = 0
        dr = 1.0d-4
        alpha = 1.2d0
        rsolv = 1.385d0

        call pedram(wrk, nwrk)

        nucind = bck_nucind
        nucdep = bck_nucdep
        natoms = bck_natoms

        lucav = -999
        call gpopen(lucav, 'CAVDATA', 'OLD', 'SEQUENTIAL',
     &              'UNFORMATTED', idummy, .FALSE.)
        read(lucav) ntess
        allocate(isphere(ntess))
        allocate(tess_coords(3,ntess))
        allocate(areas(ntess))
        read(lucav) isphere
        read(lucav) tess_coords(1,:)
        read(lucav) tess_coords(2,:)
        read(lucav) tess_coords(3,:)
        read(lucav) areas
        call gpclose(lucav, 'KEEP')

        lusurf = -999
        call gpopen(lusurf, 'SURFACE.INP', 'NEW', 'SEQUENTIAL',
     &              'FORMATTED', idummy, .FALSE.)

        write(lusurf,'(i6)') ntess
        write(lusurf,'(a)') 'AA'
        do i = 1, ntess
            write(lusurf,'(4f12.6)') tess_coords(:,i), areas(i)
        end do

        call gpclose(lusurf,'KEEP')

        call qexit('dalton_cavity')

      end subroutine dalton_cavity

!------------------------------------------------------------------------------

#endif
