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
!
      SUBROUTINE QMNPMM_LNO(NOSIM,BOVECS,CREF,CMO,XINDX,UDV,DV,UDVTR,   &
     &                      DVTR,EVECS,WORK,LWORK)
!
! Purpose:
!     Computes QM/NP/MM contribution to linear response vector.
!
! Input:
!   NOSIM  - Number of first Fock/Kohn-Sham matrices
!   BOVECs - Linear response vectors
!   CREF   - Reference state CI coeficients (Currently not used)
!   CMO    - Molecular orbitals
!   XINDX  - Indexing array of response vectors
!   UDV    - Density matrix
!   DV     - Active density matrix
!   UDVTR  - Triplet density matrix
!   DVTR   - Active triplet density matrix
!   WORK   - Temporary memory array
!   LWORK  - Size of temporary memory array
! Output:
!   EVECS  - Response vectors (XY).
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!

      use qmcmm, only: getdim_relmat, read_relmat
      use qmcmm_lr, only: get_fvvec_lr
      use qmcmm_qr, only: get_xyvec

      implicit none

#include "inforb.h"
#include "infdim.h"
#include "wrkrsp.h"
#include "infrsp.h"
#include "priunit.h"
#include "qmnpmm.h"

      integer, intent(in) :: nosim
      integer, intent(in) :: lwork

      real(8) :: BOVECS(*), CMO(*), XINDX(*), UDV(NASHDI,NASHDI)
      real(8) :: UDVTR(N2ASHX), DVTR(*), EVECS(KZYVAR,*)
      real(8) :: WORK(LWORK), DV(*), CREF(*)
      integer :: i, j, idim, idimx, kfree, lfree, ioff
      integer :: kbov, kcmo, kfvvec, kmqvec, krelmat, krxy, krxyt
!
      KFREE = 1
      LFREE = LWORK
!
!     NP/MM contribution to response matrix is zero for triplet
!     perturbations applied to singlet reference state
      IF ((NASHT.EQ.0).AND.(TRPLET)) RETURN
!     Non-iterative method
      IF (.NOT.MQITER) THEN
        CALL GETDIM_RELMAT(IDIM,.FALSE.)
!       Allocate FV and MQ vectors for all perturbed densities
        CALL MEMGET('REAL',KMQVEC,NOSIM*IDIM,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KFVVEC,NOSIM*IDIM,WORK,KFREE,LFREE)
!       Allocate unpacked temporary CMO, etc vectors
        CALL MEMGET('REAL',KCMO,NBAST*NORBT,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KBOV,NOSIM*N2ORBX,WORK,KFREE,LFREE)
!       Zero & unpack CMO and ZY vectors
        CALL DZERO(WORK(KCMO),NORBT*NBAST)
        CALL UPKCMO(CMO,WORK(KCMO))
        CALL DZERO(WORK(KBOV),NOSIM*N2ORBX)
        IF (NOSIM.GT.0) THEN
           CALL RSPZYM(NOSIM,BOVECS,WORK(KBOV))
           CALL DSCAL(NOSIM*N2ORBX,-1.0d0,WORK(KBOV),1)
        END IF
!       Determine electric field/potential vector for perturbed
!       density matrices
        CALL GET_FVVEC_LR(WORK(KFVVEC), idim, nosim, UDV,UDVTR,WORK(KCMO),WORK(KBOV), &
     &                    WORK(KFREE),LFREE)
!       Allocate and compute Relay matrix
        CALL GETDIM_RELMAT(IDIMX,.TRUE.)
        CALL MEMGET('REAL',KRELMAT,IDIMX,WORK,KFREE,LFREE)
        CALL READ_RELMAT(WORK(KRELMAT))
!       Determine induced induced dipoles and charges
        CALL DZERO(WORK(KMQVEC),NOSIM*IDIM)
        DO I=1,NOSIM
           IOFF = (I-1)*IDIM
           CALL DGEMV('N',IDIM,IDIM,1.0d0,WORK(KRELMAT),IDIM,              &
     &                WORK(KFVVEC+IOFF),1,0.0d0,WORK(KMQVEC+IOFF),1)
          if (iprtlvl > 14) then
             write(lupri, '(/,2x,a,i0)') &
                 '*** Computed MQ vector start 1st-order density ', i
             do j = 1, idim
                write(lupri, '(i8, f18.8)') j, WORK(KMQVEC+IOFF - 1 + j)
             end do
             write(lupri, '(/,2x,a)') '*** Computed MQ vector end ***'
          end if
        END DO
!       Allocate temporary XY vector contributions
        CALL MEMGET('REAL',KRXY,NOSIM*N2ORBX,WORK,KFREE,LFREE)
        CALL DZERO(WORK(KRXY),NOSIM*N2ORBX)
        IF (TRPLET) THEN
           CALL MEMGET('REAL',KRXYT,NOSIM*N2ORBX,WORK,KFREE,LFREE)
           CALL DZERO(WORK(KRXYT),NOSIM*N2ORBX)
        END IF

         ! compute xy contributions from induced dipoles and charges
         if (trplet) then
            call get_xyvec(work(kcmo),    &
                           idim,          &
                           nosim,         &
                           work(krxyt),   &
                           work(kfree),   &
                           lfree,         &
                           work(kmqvec))
         else
            call get_xyvec(work(kcmo),    &
                           idim,          &
                           nosim,         &
                           work(krxy),    &
                           work(kfree),   &
                           lfree,         &
                           work(kmqvec))
         end if
      ELSE
!        Fix me
      END IF
!     Add QM/NP/MM contributions to transformed resp. vectors
      IF (TRPLET) THEN
         CALL SLVSOR(.TRUE.,.FALSE.,NOSIM,UDVTR,EVECS(1,1),WORK(KRXY))
         CALL SLVSOR(.TRUE.,.TRUE.,NOSIM,UDV,EVECS(1,1),WORK(KRXYT))
      ELSE
         CALL SLVSOR(.TRUE.,.TRUE.,NOSIM,UDV,EVECS(1,1),WORK(KRXY))
      ENDIF
!
      CALL MEMREL('QMNPMM_LNO',WORK,1,1,KFREE,LFREE)
!
      end subroutine
      SUBROUTINE QMNPMMQRO(VEC1,VEC2,ETRS,XINDX,ZYM1,ZYM2,              &
     &                     UDV,WORK,LWORK,KZYVR,KZYV1,KZYV2,            &
     &                     IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP,              &
     &                     ISPIN0,ISPIN1,ISPIN2)

      use qmcmm, only: getdim_relmat, read_relmat
      use qmcmm_qr, only: get_fvvec_qr, get_xyvec

      implicit none

#include "dummy.h"
#include "maxorb.h"
#include "inforb.h"
#include "infdim.h"
#include "infinp.h"
#include "infvar.h"
#include "infrsp.h"
#include "infpri.h"
#include "rspprp.h"
#include "infcr.h"
#include "inftap.h"
#include "qrinf.h"
#include "mxcent.h"
#include "priunit.h"
#include "wrkrsp.h"
#include "orgcom.h"
#include "ccinftap.h"
#include "nuclei.h"
#include "infpar.h"
#include "qmnpmm.h"

      integer :: kzyvr
      integer :: kzyv1
      integer :: kzyv2
      integer :: lwork,igrsym,isymv1,isymv2,ispin0,ispin1,ispin2
      real(8) :: etrs(kzyvr),xindx(*)
      real(8) :: udv(nashdi,nashdi)
      real(8) :: zym1(*),zym2(*),work(lwork),cmo(*)
      real(8) :: vec1(kzyv1),vec2(kzyv2)
      integer :: mjwop(2,maxwop,8)
      integer :: i, j, ioff
      logical   lcon, lorb, lref
      integer :: kmqvec1
      integer :: kmqvec2
      integer :: kfvvec1
      integer :: kfvvec2
      integer :: kcref
      integer :: ktres
      integer :: kucmo
      integer :: ktlma
      integer :: ktlmb
      integer :: ktrmo
      integer :: krelmat
      integer :: kutr
      integer :: idim, idimx
      integer :: isymt, isymv, isymst, jspin, nsim
      integer :: lfree, kfree
      integer :: nzyvec, nzcvec
      integer :: ISYMDN, idum
      real(8) :: ovlap

      CALL QENTER('QMNPMMQRO')
!
      KFREE = 1
      LFREE = LWORK
!     Allocate arrays for response
      CALL MEMGET('REAL',KCREF,NCREF,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTRES,N2ORBX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KUCMO,NORBT*NBAST,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTLMA,N2ORBX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTLMB,N2ORBX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTRMO,NNORBX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KUTR,N2ORBX,WORK,KFREE,LFREE)
!     Initialize allocated arrays
      CALL DZERO(WORK(KCREF),NCREF)
      CALL DZERO(WORK(KTRES),N2ORBX)
      CALL DZERO(WORK(KUCMO),NORBT*NBAST)
      CALL DZERO(WORK(KTLMA),N2ORBX)
      CALL DZERO(WORK(KTLMB),N2ORBX)
      CALL DZERO(WORK(KTRMO),NNORBX)
      CALL DZERO(WORK(KUTR),N2ORBX)
!     Reset symmetry variables
      NSIM  = 1
      ISYMT = 1
!     Get the reference state
      CALL GETREF(WORK(KCREF),MZCONF(1))
!     Unpack the response vectors
      CALL GTZYMT(NSIM,VEC1,KZYV1,ISYMV1,ZYM1,MJWOP)
      CALL GTZYMT(NSIM,VEC2,KZYV2,ISYMV2,ZYM2,MJWOP)
!     Unpack symmetry blocked CMO
      CALL UPKCMO(CMO,WORK(KUCMO))
!     Non-iterative method
      IF (.NOT.MQITER) THEN
        CALL GETDIM_RELMAT(IDIM,.FALSE.)
!       Allocate FV and MQ vectors for all perturbed densities
        CALL MEMGET('REAL',KMQVEC1,NSIM*IDIM,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KMQVEC2,NSIM*IDIM,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KFVVEC1,NSIM*IDIM,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KFVVEC2,NSIM*IDIM,WORK,KFREE,LFREE)
!       Compute FV vectors for second order pertubed densities
        CALL GET_FVVEC_QR(WORK(KFVVEC1),WORK(KFVVEC2),IDIM,NSIM,        &
     &                    UDV,WORK(KUCMO),ISYMT,ISYMV1,ISYMV2,          &
     &                    ZYM1,ZYM2,WORK(KFREE),LFREE)
!       Allocate and compute Relay matrix
        CALL GETDIM_RELMAT(IDIMX,.TRUE.)
        CALL MEMGET('REAL',KRELMAT,IDIMX,WORK,KFREE,LFREE)
        CALL READ_RELMAT(WORK(KRELMAT))
!       Determine induced induced dipoles and charges
        CALL DZERO(WORK(KMQVEC1),NSIM*IDIM)
        CALL DZERO(WORK(KMQVEC2),NSIM*IDIM)
        DO I=1,NSIM
           IOFF = (I-1)*IDIM
           CALL DGEMV('N',IDIM,IDIM,1.0d0,WORK(KRELMAT),IDIM,              &
     &                WORK(KFVVEC1+IOFF),1,0.0d0,WORK(KMQVEC1+IOFF),1)
           CALL DGEMV('N',IDIM,IDIM,1.0d0,WORK(KRELMAT),IDIM,              &
     &                WORK(KFVVEC2+IOFF),1,0.0d0,WORK(KMQVEC2+IOFF),1)
          if (iprtlvl > 14) then
             write(lupri, '(/,2x,a,i0)') &
                 '*** Computed MQ vector start v1 1st-order density ', i
             do j = 1, idim
                write(lupri, '(i8, f18.8)') j, WORK(KMQVEC1+IOFF-1+j)
             end do
             write(lupri, '(/,2x,a)') '*** Computed MQ vector end ***'
             write(lupri, '(/,2x,a,i0)') &
                 '*** Computed MQ vector start v2 1st-order density ', i
             do j = 1, idim
                write(lupri, '(i8, f18.8)') j, WORK(KMQVEC2+IOFF-1+j)
             end do
             write(lupri, '(/,2x,a)') '*** Computed MQ vector end ***'
          end if
        END DO

        ! determine mm region contribution to qm region potential from
        ! second order density
        call get_xyvec(work(kucmo),   &
                       idim,          &
                       nsim,          &
                       work(ktres),   &
                       work(kfree),   &
                       lfree,         &
                       work(kfvvec1), &
                       work(kfvvec2), &
                       isymt,         &
                       isymv2,        &
                       zym2)

      ELSE
!       FIX ME: ITERATIVE METHOD
      END IF
!     Set up paramterers for quadratic response gradient formation
      ISYMDN = 1
      OVLAP  = 1.0d0
      JSPIN  = 0
      ISYMV  = IREFSY
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF ( ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LORB   = ( MZWOPT(IGRSYM) .GT. 0 )
      LREF = .TRUE.
      NZYVEC = NCREF
      NZCVEC = NCREF
!     Compute gradient
      CALL RSP1GR(NSIM,KZYVR,IDUM,JSPIN,IGRSYM,JSPIN,ISYMV,ETRS,        &
     &            WORK(KCREF),NZYVEC,NZCVEC,OVLAP,ISYMDN,UDV,           &
     &            WORK(KTRES),XINDX,MJWOP,WORK(KFREE),LFREE,            &
     &            LORB,LCON,LREF)
!
      CALL MEMREL('QMNPMM_QRO',WORK,1,1,KFREE,LFREE)
!
      CALL QEXIT('QMNPMMQRO')
!
      end subroutine
