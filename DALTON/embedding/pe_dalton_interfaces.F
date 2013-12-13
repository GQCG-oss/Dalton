!
!...   Copyright (c) 2013 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2013 (2013), see http://daltonprogram.org"
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
      subroutine pe_ifc_dalinp_reader(word)
        use pe_variables, only: peqm
        use polarizable_embedding, only: pe_dalton_input
#include "priunit.h"
        character(len=7), intent(inout) :: word
        peqm = .true.
        call pe_dalton_input(word, lucmd, lupri)
      end subroutine pe_ifc_dalinp_reader

      subroutine pe_ifc_init()
        use polarizable_embedding, only: pe_init
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
        call pe_init(lupri, cord(1:3,1:natoms), charge(1:natoms))
      end subroutine pe_ifc_init

      subroutine pe_ifc_fock(denmats, fckmats, nmats, energies,
     &                       work, nwrk)
        use polarizable_embedding, only: pe_master
#include "inforb.h"
        integer, intent(in) :: nmats, nwrk
        real*8, dimension(nmats*nnbasx), intent(in) :: denmats
        real*8, dimension(nmats*nnbasx), intent(out) :: fckmats
        real*8, dimension(nmats), intent(out) :: energies
        real*8, dimension(nwrk), intent(inout) :: work
        call pe_master(runtype='fock', denmats=denmats, fckmats=fckmats,
     &                 nmats=nmats, energies=energies, dalwrk=work)
      end subroutine pe_ifc_fock

      subroutine pe_ifc_energy(denmats, nmats, work, nwrk)
        use polarizable_embedding, only: pe_master
#include "inforb.h"
        integer, intent(in) :: nmats, nwrk
        real*8, dimension(nmats*nnbasx), intent(in) :: denmats
        real*8, dimension(nwrk), intent(inout) :: work
        call pe_master(runtype='energy', denmats=denmats, nmats=nmats,
     &                 dalwrk=work)
      end subroutine pe_ifc_energy

#if defined(VAR_MPI)
      subroutine pe_ifc_mpi(work, nwrk, runtype)
        use polarizable_embedding, only: pe_mpi
        integer, intent(in) :: nwrk, runtype
        real*8, dimension(nwrk), intent(inout) :: work
        call pe_mpi(work, runtype)
      end subroutine pe_ifc_mpi
#endif

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

        if (.not. tdhf) then
            write(lupri, *) 'ERROR: PE-MCSCF response not implemented.'
            call quit('ERROR: PE-MCSCF response not implemented.')
        end if

        if (nosim > 0) then
            call pe_rsplno(ncsim, nosim, bovecs, cref, cmo, cindx, 
     &                     udv, dv, udvtr, dvtr, sovecs, wrk, nwrk)
        end if

        call qexit('pe_rsplin')

      end subroutine pe_rsplin

!------------------------------------------------------------------------------

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
        real*8, dimension(:), allocatable :: dcao, dvao
!        real*8, dimension(:), allocatable :: dcaotr, dvaotr
        real*8, dimension(:), allocatable :: daos, fckaos
!        real*8, dimension(:), allocatable :: daotrs, fckaotrs
        real*8, dimension(:), allocatable :: evec
        real*8, dimension(:,:), allocatable :: ubovecs, evecs
!        real*8, dimension(:,:), allocatable :: evectrs, eactrs
        logical :: lexist

        ! return if no polarization
        if (.not. pe_polar) then
            return
        ! no polarization for triplet excitations in closed shell
        else if ((nasht == 0) .and. trplet) then
            return
        ! ground state polarization approximation
        else if (pe_gspol) then
            return
        ! triplet response for open shell systems not ready yet
        else if ((nasht > 0) .and. trplet) then
            call quit('ERROR: triplet operators for open shell'//
     &                ' systems not implemented')
        end if

        call qenter('pe_rsplno')

        allocate(ubovecs(n2orbx,nosim))
        call rspzym(nosim, bovecs, ubovecs)

        ubovecs = - ubovecs

        allocate(dcao(n2basx), dvao(n2basx), daos(nosim*nnbasx))
!        if (trplet) then
!            allocate(dcaotr(n2basx), dvaotr(n2basx),
!     &               daotrs(nosim*nnbasx))
!        end if

        do i = 1, nosim
            j = (i - 1) * nnbasx + 1
            call deq27(cmo, ubovecs(:,i), udv, dcao, dvao,
     &                 wrk, nwrk)
!            if (trplet) then
!                call deq27(cmo, ubovecs(:,i), udvtr, dcaotr, dvaotr,
!     &                     wrk, nwrk)
!            end if
            if (nasht > 0) then
                dcao = dcao + dvao
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
        allocate(evecs(n2orbx,nosim))
        evecs = 0.0d0
!        if (trplet) then
!            allocate(evectrs(n2orbx,nosim), eactrs(n2ashx,nosim))
!        end if

        do i = 1, nosim
            j = (i - 1) * nnbasx + 1
            call uthu(fckaos(j), evec, cmo, wrk, nbast, norbt)
            call dsptsi(norbt, evec, evecs(:,i))
!            if (trplet) then
!                call uthu(fckaotrs(j), evec, cmo, wrk, nbast, norbt)
!                call dsptsi(norbt, evec, evectrs(:,i))
!            end if
        end do

        deallocate(evec)
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

         call gtzymt(1, vecb, kzyvb, isymb, zymb, mjwop)
         call gtzymt(1, vecc, kzyvc, isymc, zymc, mjwop)

         allocate(udcao(n2basx))
         allocate(ufcmo(n2orbx))
         allocate(dcaos(4*nnbasx))
         dcaos = 0.0d0

         udcao = 0.0d0
         call cdens1(isymb, cmo, zymb, udcao, wrk, nwrk)
         call dgefsp(nbast, udcao, dcaos(1:nnbasx))

         udcao = 0.0d0
         call cdens2(isymb, isymc, cmo, zymb, zymc, udcao,
     &               wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
         call dgefsp(nbast, udcao, dcaos(nnbasx+1:2*nnbasx))

         udcao = 0.0d0
         call cdens1(isymc, cmo, zymc, udcao, wrk, nwrk)
         call dgefsp(nbast, udcao, dcaos(2*nnbasx+1:3*nnbasx))

         udcao = 0.0d0
         call cdens2(isymc, isymb, cmo, zymc, zymb, udcao,
     &               wrk(1:n2basx), wrk(n2basx+1:2*n2basx), ufcmo)
         call dgefsp(nbast, udcao, dcaos(3*nnbasx+1:4*nnbasx))

         deallocate(udcao)

         allocate(fcaos(4*nnbasx))
         fcaos = 0.0d0
         call pe_master(runtype='response',
     &                  denmats=dcaos, fckmats=fcaos,
     &                  nmats=4, dalwrk=wrk(1:nwrk))
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
     &               wrk, idum, idum, 1.0d0, 1, udv, ufcmo, xindx,
     &               mjwop, wrk, nwrk, .true., .false., .false.)

         deallocate(fcaos, fcmo, ufcmo)

      end subroutine pe_rspqro
#endif
