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
!  /* Deck qmnpmm_nlo */
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

#include "implicit.h"
#include "inforb.h"
#include "infdim.h"
#include "wrkrsp.h"
#include "infrsp.h"
#include "priunit.h"
#include "qmnpmm.h"
!
      DIMENSION BOVECS(*), CMO(*), XINDX(*), UDV(NASHDI,NASHDI)
      DIMENSION UDVTR(N2ASHX), DVTR(*), EVECS(KZYVAR,*)
      DIMENSION WORK(LWORK), DV(*), CREF(*)
!
      PARAMETER (DM1 = -1.0D0, D1 = 1.0D0, D0 = 0.0D0)
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
           CALL DSCAL(NOSIM*N2ORBX,DM1,WORK(KBOV),1)
        END IF
!       Determine electric field/potential vector for perturbed
!       density matrices
        CALL GET_FVVEC_LR(UDV,UDVTR,WORK(KCMO),WORK(KBOV),WORK(KFVVEC), &
     &                    IDIM,NOSIM,WORK(KFREE),LFREE)
!       Allocate and compute Relay matrix
        CALL GETDIM_RELMAT(IDIMX,.TRUE.)
        CALL MEMGET('REAL',KRELMAT,IDIMX,WORK,KFREE,LFREE)
        CALL READ_RELMAT(WORK(KRELMAT))
!       Determine induced induced dipoles and charges
        CALL DZERO(WORK(KMQVEC),NOSIM*IDIM)
        DO I=1,NOSIM
           IOFF = (I-1)*IDIM
           CALL DGEMV('N',IDIM,IDIM,D1,WORK(KRELMAT),IDIM,              &
     &                WORK(KFVVEC+IOFF),1,D0,WORK(KMQVEC+IOFF),1)
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
!       Compute XY contributions from induced dipoles and charges
        IF (TRPLET) THEN
           CALL GET_XYVEC_LR(WORK(KMQVEC),WORK(KCMO),IDIM,NOSIM,        &
     &                       WORK(KRXYT),WORK(KFREE),LFREE)
        ELSE
           CALL GET_XYVEC_LR(WORK(KMQVEC),WORK(KCMO),IDIM,NOSIM,        &
     &                       WORK(KRXY),WORK(KFREE),LFREE)
        END IF
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
      RETURN
      end subroutine
!  /* Deck qmnpmmqro */
      SUBROUTINE QMNPMMQRO(VEC1,VEC2,ETRS,XINDX,ZYM1,ZYM2,              &
     &                     UDV,WORK,LWORK,KZYVR,KZYV1,KZYV2,            &
     &                     IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP,              &
     &                     ISPIN0,ISPIN1,ISPIN2)

      use qmcmm, only: getdim_relmat, read_relmat

#include "implicit.h"
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
!
      DIMENSION ETRS(KZYVR),XINDX(*)
      DIMENSION UDV(NASHDI,NASHDI)
      DIMENSION ZYM1(*),ZYM2(*),WORK(LWORK),CMO(*)
      DIMENSION VEC1(KZYV1),VEC2(KZYV2)
      DIMENSION MJWOP(2,MAXWOP,8)
      LOGICAL   LCON, LORB, LREF
!
      PARAMETER (D1 = 1.0D0, D0 = 0.0D0)
!
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
           CALL DGEMV('N',IDIM,IDIM,D1,WORK(KRELMAT),IDIM,              &
     &                WORK(KFVVEC1+IOFF),1,D0,WORK(KMQVEC1+IOFF),1)
           IF (IPRTLVL.GE.15) THEN
              WRITE(LUPRI,'(/,2X,A,I2,A)') '*** Computed MQ1 vector ',  &
     &              I,' perturbed first order density ***'
              CALL OUTPUT(WORK(KMQVEC1+IOFF),1,IDIM,1,1,IDIM,1,1,LUPRI)
           END IF
           CALL DGEMV('N',IDIM,IDIM,D1,WORK(KRELMAT),IDIM,              &
     &                WORK(KFVVEC2+IOFF),1,D0,WORK(KMQVEC2+IOFF),1)
           IF (IPRTLVL.GE.15) THEN
              WRITE(LUPRI,'(/,2X,A,I2,A)') '*** Computed MQ2 vector ',  &
     &              I,' perturbed first order density ***'
              CALL OUTPUT(WORK(KMQVEC2+IOFF),1,IDIM,1,1,IDIM,1,1,LUPRI)
           END IF
        END DO
!       Determine MM region contribution to QM region potential from
!       second order density
        CALL GET_XYVEC_QR(WORK(KFVVEC1),WORK(KFVVEC2),WORK(KUCMO),      &
     &                    IDIM,NSIM,WORK(KTRES),ISYMT,ISYMV2,ZYM2,      &
     &                    WORK(KFREE),LFREE)

      ELSE
!       FIX ME: ITERATIVE METHOD
      END IF
!     Set up paramterers for quadratic response gradient formation
      ISYMDN = 1
      OVLAP  = D1
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
      RETURN
      end subroutine

!  /* Deck get_fvvec_lr */
      SUBROUTINE GET_FVVEC_LR(UDV,UDVTR,CMO,BOVECS,FVVEC,IDIM,NOSIM,    &
     &                        WORK,LWORK)
!
! Purpose:
!     Computes electric field/potential vector generated by first
!     order perturbed density matrix.
!
! Input:
!   WORK   - Temporary memory array
!   LWORK  - Size of temporary memory array.
! Output:
!   FVVEC  - Electric field/potential vector at NP/MM centers.
!   IDIM   - Size of electric field/potential vector
!   NOSIM  - Number of perturbed density matrices-
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "qmnpmm.h"
#include "inforb.h"
#include "infdim.h"
#include "mxcent.h"
#include "iratdef.h"
#include "nuclei.h"
#include "orgcom.h"
#include "qm3.h"
#include "infrsp.h"
!
      DIMENSION UDV(NASHDI,NASHDI), UDVTR(N2ASHX), CMO(NORBT*NBAST)
      DIMENSION BOVECS(NOSIM*N2ORBX), FVVEC(*), WORK(LWORK)
!
      LOGICAL TOFILE,TRIMAT,EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      DIMENSION RSAVORG(3)
!
      KFREE = 1
      LFREE = LWORK
!
      CALL DZERO(FVVEC,IDIM*NOSIM)
!     Save origin coordinates
      RSAVORG(1) = DIPORG(1)
      RSAVORG(2) = DIPORG(2)
      RSAVORG(3) = DIPORG(3)
!     Compute electric field if needed
      IF (DONPPOL) THEN
!       Allocate integrals buffer
        CALL MEMGET('REAL',KINTAO,3*NNBASX,WORK,KFREE,LFREE)
!       Allocate temp. matrices used for integrals tranformations
        CALL MEMGET('REAL',KTRMO,NNORBX,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KUTR,N2ORBX,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KUTRX,N2ORBX,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KURXAC,N2ASHX,WORK,KFREE,LFREE)
!       Set up integrals calculation parameters
        KPATOM = 0
        NOCOMP = 3
        TOFILE = .FALSE.
        TRIMAT = .TRUE.
        EXP1VL = .FALSE.
        RUNQM3 = .TRUE.
!       Loop over perturbed first order densities
        DO I=1,NOSIM
           IOFF = (I-1)*IDIM
           ISIMOFF = (I-1)*N2ORBX+1
!          Loop over NP centers
           DO J=1,TNPATM
             JOFF = IOFF+(J-1)*3
             DIPORG(1) = NPCORD(1,J)
             DIPORG(2) = NPCORD(2,J)
             DIPORG(3) = NPCORD(3,J)
             CALL DZERO(WORK(KINTAO),3*NNBASX)
             CALL GET1IN(WORK(KINTAO),'NEFIELD',NOCOMP,WORK(KFREE),     &
     &                   LFREE,LABINT,INTREP,INTADR,J,TOFILE,KPATOM,    &
     &                   TRIMAT,DUMMY,EXP1VL,DUMMY,0)
!            Compute X component of electric field
!            zero temp. arrays
             CALL DZERO(WORK(KTRMO),NNORBX)
             CALL DZERO(WORK(KUTR),N2ORBX)
             CALL DZERO(WORK(KUTRX),N2ORBX)
             CALL DZERO(WORK(KURXAC),N2ASHX)
!            transfor integrals
             CALL UTHU(WORK(KINTAO),WORK(KTRMO),CMO,WORK(KFREE),NBAST,  &
     &                 NORBT)
             CALL DSPTSI(NORBT,WORK(KTRMO),WORK(KUTR))
             CALL ONEXH1(BOVECS(ISIMOFF),WORK(KUTR),WORK(KUTRX))
!
             IF (NASHT.GT.0) CALL GETACQ(WORK(KUTRX),WORK(KURXAC))
             IF (TRPLET) THEN
                 FVVEC(JOFF+1) = SLVTLM(UDVTR,WORK(KURXAC),WORK(KUTRX), &
     &                                  TAC)
             ELSE
                 FVVEC(JOFF+1) = SLVQLM(UDV,WORK(KURXAC),WORK(KUTRX),   &
     &                                  TAC)
             ENDIF
!            Compute Y component of electric field
!            zero temp. arrays
             CALL DZERO(WORK(KTRMO),NNORBX)
             CALL DZERO(WORK(KUTR),N2ORBX)
             CALL DZERO(WORK(KUTRX),N2ORBX)
             CALL DZERO(WORK(KURXAC),N2ASHX)
!            transfor integrals
             CALL UTHU(WORK(KINTAO+NNBASX),WORK(KTRMO),CMO,WORK(KFREE), &
     &                 NBAST,NORBT)
             CALL DSPTSI(NORBT,WORK(KTRMO),WORK(KUTR))
             CALL ONEXH1(BOVECS(ISIMOFF),WORK(KUTR),WORK(KUTRX))
!
             IF (NASHT.GT.0) CALL GETACQ(WORK(KUTRX),WORK(KURXAC))
             IF (TRPLET) THEN
                 FVVEC(JOFF+2) = SLVTLM(UDVTR,WORK(KURXAC),WORK(KUTRX), &
     &                                  TAC)
             ELSE
                 FVVEC(JOFF+2) = SLVQLM(UDV,WORK(KURXAC),WORK(KUTRX),   &
     &                                  TAC)
             ENDIF
!            Compute Z component of electric field
!            zero temp. arrays
             CALL DZERO(WORK(KTRMO),NNORBX)
             CALL DZERO(WORK(KUTR),N2ORBX)
             CALL DZERO(WORK(KUTRX),N2ORBX)
             CALL DZERO(WORK(KURXAC),N2ASHX)
!            transfor integrals
             CALL UTHU(WORK(KINTAO+2*NNBASX),WORK(KTRMO),CMO,           &
     &                 WORK(KFREE),NBAST,NORBT)
             CALL DSPTSI(NORBT,WORK(KTRMO),WORK(KUTR))
             CALL ONEXH1(BOVECS(ISIMOFF),WORK(KUTR),WORK(KUTRX))
!
             IF (NASHT.GT.0) CALL GETACQ(WORK(KUTRX),WORK(KURXAC))
             IF (TRPLET) THEN
                 FVVEC(JOFF+3) = SLVTLM(UDVTR,WORK(KURXAC),WORK(KUTRX), &
     &                                  TAC)
             ELSE
                 FVVEC(JOFF+3) = SLVQLM(UDV,WORK(KURXAC),WORK(KUTRX),   &
     &                                  TAC)
             ENDIF
           END DO
        END DO
        RUNQM3 = .FALSE.
        CALL MEMREL('GET_FVVEC_LR',WORK,1,1,KFREE,LFREE)
      END IF
      IF (DONPCAP) THEN
        ISTART = 0
        IF (DONPPOL) ISTART = 3*TNPATM
!       Fix me MM region shift
!       Allocate integrals buffer
        CALL MEMGET('REAL',KINTAO,NNBASX,WORK,KFREE,LFREE)
!       Allocate temp. matrices used for integrals tranformations
        CALL MEMGET('REAL',KTRMO,NNORBX,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KUTR,N2ORBX,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KUTRX,N2ORBX,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KURXAC,N2ASHX,WORK,KFREE,LFREE)
!       Set up integrals calculation parameters
        KPATOM = 0
        NOCOMP = 1
        TOFILE = .FALSE.
        TRIMAT = .TRUE.
        EXP1VL = .FALSE.
        RUNQM3 = .TRUE.
!       Loop over perturbed first order densities
        DO I=1,NOSIM
           IOFF = (I-1)*IDIM
           ISIMOFF = (I-1)*N2ORBX+1
!          Loop over NP centers
           DO J=1,TNPATM
             JOFF = IOFF+ISTART+J
             DIPORG(1) = NPCORD(1,J)
             DIPORG(2) = NPCORD(2,J)
             DIPORG(3) = NPCORD(3,J)
             CALL DZERO(WORK(KINTAO),NNBASX)
             CALL GET1IN(WORK(KINTAO),'NPETES ',NOCOMP,WORK(KFREE),     &
     &                   LFREE,LABINT,INTREP,INTADR,J,TOFILE,KPATOM,    &
     &                   TRIMAT,DUMMY,EXP1VL,DUMMY,0)
!            Zero temporary arrays
             CALL DZERO(WORK(KTRMO),NNORBX)
             CALL DZERO(WORK(KUTR),N2ORBX)
             CALL DZERO(WORK(KUTRX),N2ORBX)
             CALL DZERO(WORK(KURXAC),N2ASHX)
!            Transform integrals
             CALL UTHU(WORK(KINTAO),WORK(KTRMO),CMO,WORK(KFREE),NBAST,  &
     &                 NORBT)
             CALL DSPTSI(NORBT,WORK(KTRMO),WORK(KUTR))
             CALL ONEXH1(BOVECS(ISIMOFF),WORK(KUTR),WORK(KUTRX))
!
             IF (NASHT.GT.0) CALL GETACQ(WORK(KUTRX),WORK(KURXAC))
             IF (TRPLET) THEN
                 FVVEC(JOFF) = SLVTLM(UDVTR,WORK(KURXAC),WORK(KUTRX),   &
     &                                TAC)
             ELSE
                 FVVEC(JOFF) = SLVQLM(UDV,WORK(KURXAC),WORK(KUTRX),     &
     &                                TAC)
             ENDIF
            END DO
!           Set Lagrangian for charge equilibration
            DO IBLK=1,TNPBLK
              FVVEC(IOFF+IDIM) = FVVEC(IOFF+IDIM)+NPCHRG(IBLK)
            END DO
        END DO
        RUNQM3 = .FALSE.
        CALL MEMREL('GET_FVVEC_LR',WORK,1,1,KFREE,LFREE)
      END IF
!     Restore origin coordinates
      DIPORG(1) = RSAVORG(1)
      DIPORG(2) = RSAVORG(2)
      DIPORG(3) = RSAVORG(3)
!     Print final FV vector
      DO I=1,NOSIM
         if (iprtlvl > 14 .and. .not. mqiter) then
            IOFF = (I-1)*IDIM+1
            write(lupri, '(/,2x,a,i0)') &
                '*** Computed FV vector start 1st-order density ', i
            do j = 1, idim
               write(lupri, '(i8, f18.8)') j, fvvec(IOFF + j - 1)
            end do
            write(lupri, '(/,2x,a)') '*** Computed FV vector end ***'
         end if
      END DO
!
      RETURN
      end subroutine

!  /* Deck get_fvvec_qr */
      SUBROUTINE GET_FVVEC_QR(FVVEC1,FVVEC2,IDIM,NSIM,UDV,UCMO,         &
     &                        ISYMT,ISYMV1,ISYMV2,ZYM1,ZYM2,WORK,LWORK)
!
! Purpose:
!     Computes electric field/potential vector generated by second
!     order perturbed density matrix.
!
! Input:
!   WORK   - Temporary memory array
!   LWORK  - Size of temporary memory array.
! Output:
!   FVVEC1 - Electric field/potential vector at NP/MM centers.
!   FVVEC2 - Electric field/potential vector at NP/MM centers.
!   IDIM   - Size of electric field/potential vector
!   NSIM   - Number of perturbed density matrices.
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "qmnpmm.h"
#include "inforb.h"
#include "infdim.h"
#include "mxcent.h"
#include "iratdef.h"
#include "nuclei.h"
#include "orgcom.h"
#include "qm3.h"
#include "infrsp.h"
!
      DIMENSION FVVEC1(*), FVVEC2(*), WORK(LWORK)
      DIMENSION UDV(NASHDI,NASHDI),UCMO(*),ZYM1(*),ZYM2(*)
!
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
!
      LOGICAL TOFILE,TRIMAT,EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      DIMENSION RSAVORG(3)
!
      KFREE = 1
      LFREE = LWORK
!
      CALL DZERO(FVVEC1,IDIM*NSIM)
      CALL DZERO(FVVEC2,IDIM*NSIM)
!     Save origin coordinates
      RSAVORG(1) = DIPORG(1)
      RSAVORG(2) = DIPORG(2)
      RSAVORG(3) = DIPORG(3)
!     Compute electric field if needed
      IF (DONPPOL) THEN
!       Allocate integrals buffer
        CALL MEMGET('REAL',KINTAO,3*NNBASX,WORK,KFREE,LFREE)
!       Allocate temp. matrices used for integrals tranformations
        CALL MEMGET('REAL',KTRMO,NNORBX,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KUTR,N2ORBX,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KTLMA,N2ORBX,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KTLMB,N2ORBX,WORK,KFREE,LFREE)
!       Set up integrals calculation parameters
        KPATOM = 0
        NOCOMP = 3
        TOFILE = .FALSE.
        TRIMAT = .TRUE.
        EXP1VL = .FALSE.
        RUNQM3 = .TRUE.
!       Loop over perturbed first order densities
        DO I=1,NSIM
           IOFF = (I-1)*IDIM
!          Loop over NP centers
           DO J=1,TNPATM
             JOFF = IOFF+(J-1)*3
             DIPORG(1) = NPCORD(1,J)
             DIPORG(2) = NPCORD(2,J)
             DIPORG(3) = NPCORD(3,J)
             CALL DZERO(WORK(KINTAO),3*NNBASX)
             CALL GET1IN(WORK(KINTAO),'NEFIELD',NOCOMP,WORK(KFREE),     &
     &                   LFREE,LABINT,INTREP,INTADR,J,TOFILE,KPATOM,    &
     &                   TRIMAT,DUMMY,EXP1VL,DUMMY,0)
!            Compute X component of electric field
             CALL DZERO(WORK(KTRMO),NNORBX)
             CALL DZERO(WORK(KUTR),N2ORBX)
!            Transform integrals
             CALL UTHU(WORK(KINTAO),WORK(KTRMO),UCMO,WORK(KFREE),NBAST, &
     &                 NORBT)
             CALL DSPTSI(NORBT,WORK(KTRMO),WORK(KUTR))
!            Determine electric field component size
             F1VAL = D0
             F2VAL = D0
             IF (ISYMT.EQ.ISYMV1) THEN
                CALL DZERO(WORK(KTLMA),N2ORBX)
                CALL OITH1(ISYMV1,ZYM1,WORK(KUTR),WORK(KTLMA),ISYMT)
                CALL MELONE(WORK(KTLMA),1,UDV,D1,F1VAL,200,             &
     &                      'QMNPQRO')
                FVVEC1(JOFF+1) = F1VAL
             END IF
             IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
                CALL DZERO(WORK(KTLMA),N2ORBX)
                CALL DZERO(WORK(KTLMB),N2ORBX)
                CALL OITH1(ISYMV1,ZYM1,WORK(KUTR),WORK(KTLMA),ISYMT)
                CALL OITH1(ISYMV2,ZYM2,WORK(KTLMA),WORK(KTLMB),ISYMV2)
                CALL MELONE(WORK(KTLMB),1,UDV,D1,F2VAL,200,             &
     &                      'QMNPQRO')
                FVVEC2(JOFF+1) = F2VAL
             ENDIF
!            Compute Y component of electric field
             CALL DZERO(WORK(KTRMO),NNORBX)
             CALL DZERO(WORK(KUTR),N2ORBX)
!            Transform integrals
             CALL UTHU(WORK(KINTAO+NNBASX),WORK(KTRMO),UCMO,WORK(KFREE),&
     &                 NBAST,NORBT)
             CALL DSPTSI(NORBT,WORK(KTRMO),WORK(KUTR))
!            Determine electric field component size
             F1VAL = D0
             F2VAL = D0
             IF (ISYMT.EQ.ISYMV1) THEN
                CALL DZERO(WORK(KTLMA),N2ORBX)
                CALL OITH1(ISYMV1,ZYM1,WORK(KUTR),WORK(KTLMA),ISYMT)
                CALL MELONE(WORK(KTLMA),1,UDV,D1,F1VAL,200,             &
     &                      'QMNPQRO')
                FVVEC1(JOFF+2) = F1VAL
             END IF
             IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
                CALL DZERO(WORK(KTLMA),N2ORBX)
                CALL DZERO(WORK(KTLMB),N2ORBX)
                CALL OITH1(ISYMV1,ZYM1,WORK(KUTR),WORK(KTLMA),ISYMT)
                CALL OITH1(ISYMV2,ZYM2,WORK(KTLMA),WORK(KTLMB),ISYMV2)
                CALL MELONE(WORK(KTLMB),1,UDV,D1,F2VAL,200,             &
     &                      'QMNPQRO')
                FVVEC2(JOFF+2) = F2VAL
             ENDIF
!            Compute Z component of electric field
             CALL DZERO(WORK(KTRMO),NNORBX)
             CALL DZERO(WORK(KUTR),N2ORBX)
!            Transform integrals
             CALL UTHU(WORK(KINTAO+2*NNBASX),WORK(KTRMO),UCMO,          &
     &                 WORK(KFREE),NBAST,NORBT)
             CALL DSPTSI(NORBT,WORK(KTRMO),WORK(KUTR))
!            Determine electric field component size
             F1VAL = D0
             F2VAL = D0
             IF (ISYMT.EQ.ISYMV1) THEN
                CALL DZERO(WORK(KTLMA),N2ORBX)
                CALL OITH1(ISYMV1,ZYM1,WORK(KUTR),WORK(KTLMA),ISYMT)
                CALL MELONE(WORK(KTLMA),1,UDV,D1,F1VAL,200,             &
     &                      'QMNPQRO')
                FVVEC1(JOFF+3) = F1VAL
             END IF
             IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
                CALL DZERO(WORK(KTLMA),N2ORBX)
                CALL DZERO(WORK(KTLMB),N2ORBX)
                CALL OITH1(ISYMV1,ZYM1,WORK(KUTR),WORK(KTLMA),ISYMT)
                CALL OITH1(ISYMV2,ZYM2,WORK(KTLMA),WORK(KTLMB),ISYMV2)
                CALL MELONE(WORK(KTLMB),1,UDV,D1,F2VAL,200,             &
     &                      'QMNPQRO')
                FVVEC2(JOFF+3) = F2VAL
             ENDIF
           END DO
        END DO
        RUNQM3 = .FALSE.
        CALL MEMREL('GET_FVVEC_QR',WORK,1,1,KFREE,LFREE)
      END IF
      IF (DONPCAP) THEN
        ISTART = 0
        IF (DONPPOL) ISTART = 3*TNPATM
!       Fix me MM region shift
!       Allocate integrals buffer
        CALL MEMGET('REAL',KINTAO,NNBASX,WORK,KFREE,LFREE)
!       Allocate temp. matrices used for integrals tranformations
        CALL MEMGET('REAL',KTRMO,NNORBX,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KUTR,N2ORBX,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KTLMA,N2ORBX,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KTLMB,N2ORBX,WORK,KFREE,LFREE)
!       Set up integrals calculation parameters
        KPATOM = 0
        NOCOMP = 1
        TOFILE = .FALSE.
        TRIMAT = .TRUE.
        EXP1VL = .FALSE.
        RUNQM3 = .TRUE.
!       Loop over perturbed first order densities
        DO I=1,NSIM
           IOFF = (I-1)*IDIM
!          Loop over NP centers
           DO J=1,TNPATM
             JOFF = IOFF+ISTART+J
             DIPORG(1) = NPCORD(1,J)
             DIPORG(2) = NPCORD(2,J)
             DIPORG(3) = NPCORD(3,J)
             CALL DZERO(WORK(KINTAO),NNBASX)
             CALL GET1IN(WORK(KINTAO),'NPETES ',NOCOMP,WORK(KFREE),     &
     &                   LFREE,LABINT,INTREP,INTADR,J,TOFILE,KPATOM,    &
     &                   TRIMAT,DUMMY,EXP1VL,DUMMY,0)
!            Zero temporary arrays
             CALL DZERO(WORK(KTRMO),NNORBX)
             CALL DZERO(WORK(KUTR),N2ORBX)
!            Transform integrals
             CALL UTHU(WORK(KINTAO),WORK(KTRMO),UCMO,WORK(KFREE),NBAST, &
     &                 NORBT)
             CALL DSPTSI(NORBT,WORK(KTRMO),WORK(KUTR))
!            Determine electric field component size
             F1VAL = D0
             F2VAL = D0
             IF (ISYMT.EQ.ISYMV1) THEN
                CALL DZERO(WORK(KTLMA),N2ORBX)
                CALL OITH1(ISYMV1,ZYM1,WORK(KUTR),WORK(KTLMA),ISYMT)
                CALL MELONE(WORK(KTLMA),1,UDV,D1,F1VAL,200,             &
     &                      'QMNPQRO')
                FVVEC1(JOFF) = F1VAL
             END IF
             IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
                CALL DZERO(WORK(KTLMA),N2ORBX)
                CALL DZERO(WORK(KTLMB),N2ORBX)
                CALL OITH1(ISYMV1,ZYM1,WORK(KUTR),WORK(KTLMA),ISYMT)
                CALL OITH1(ISYMV2,ZYM2,WORK(KTLMA),WORK(KTLMB),ISYMV2)
                CALL MELONE(WORK(KTLMB),1,UDV,D1,F2VAL,200,             &
     &                      'QMNPQRO')
                FVVEC2(JOFF) = F2VAL
             ENDIF
           END DO
!          Set Lagrangian for charge equilibration
           DO IBLK=1,TNPBLK
              FVVEC1(IOFF+IDIM) = FVVEC1(IOFF+IDIM)+NPCHRG(IBLK)
              FVVEC2(IOFF+IDIM) = FVVEC2(IOFF+IDIM)+NPCHRG(IBLK)
           END DO
        END DO
        RUNQM3 = .FALSE.
        CALL MEMREL('GET_FVVEC_QR',WORK,1,1,KFREE,LFREE)
      END IF
!     Restore origin coordinates
      DIPORG(1) = RSAVORG(1)
      DIPORG(2) = RSAVORG(2)
      DIPORG(3) = RSAVORG(3)
!     Print final FV vector
      DO I=1,NSIM
         IF ((IPRTLVL.GE.15).AND.(.NOT.MQITER)) THEN
            WRITE(LUPRI,'(/,2X,A,I2,A)') '*** Computed FV1 vector for ',&
     &            I, ' pertubed second order density matrix ***'
            IOFF = (I-1)*IDIM+1
            CALL OUTPUT(FVVEC1(IOFF),1,IDIM,1,1,IDIM,1,1,LUPRI)
            WRITE(LUPRI,'(/,2X,A,I2,A)') '*** Computed FV2 vector for ',&
     &            I, ' pertubed second order density matrix ***'
            CALL OUTPUT(FVVEC2(IOFF),1,IDIM,1,1,IDIM,1,1,LUPRI)
         END IF
      END DO
!
      RETURN
      end subroutine

!  /* Deck get_xyvec_lr */
      SUBROUTINE GET_XYVEC_LR(FMQVEC,CMO,IDIM,NOSIM,FXYVEC,WORK,LWORK)
!
! Purpose:
!     Computes contribution to XY vector from induced dipoles moments
!     and charges.
!
! Input:
!   FMQVEC - Induced dipole moments and charges
!   CMO    - Molecular orbitals
!   WORK   - Temporary memory array
!   LWORK  - Size of temporary memory array
! Output:
!   FXYVEC - Contribution to XY vector from induced dipole moments and
!            charges
!   IDIM   - Size of electric field/potential vector
!   NOSIM  - Number of perturbed density matrices
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "qmnpmm.h"
#include "inforb.h"
#include "infdim.h"
#include "mxcent.h"
#include "iratdef.h"
#include "nuclei.h"
#include "orgcom.h"
#include "qm3.h"
#include "infrsp.h"
!
      DIMENSION FMQVEC(NOSIM*IDIM), CMO(NORBT*NBAST)
      DIMENSION FXYVEC(NOSIM*N2ORBX), WORK(LWORK)
!
      LOGICAL TOFILE,TRIMAT,EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      DIMENSION RSAVORG(3)
!
      KFREE = 1
      LFREE = LWORK
!     Save origin coordinates
      RSAVORG(1) = DIPORG(1)
      RSAVORG(2) = DIPORG(2)
      RSAVORG(3) = DIPORG(3)
!     Induced dipole moment in NP region interaction with QM region
      IF (DONPPOL.AND.NOVDAMP) THEN
!        Electronic interaction part
         CALL MEMGET('REAL',KINTAO,3*NNBASX,WORK,KFREE,LFREE)
!        Allocate temp. matrices used for integrals tranformations
         CALL MEMGET('REAL',KTRMO,NNORBX,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KUTR,N2ORBX,WORK,KFREE,LFREE)
!        Set integrals evaluation flags
         KPATOM = 0
         NOCOMP = 3
         TOFILE = .FALSE.
         TRIMAT = .TRUE.
         EXP1VL = .FALSE.
         RUNQM3 = .TRUE.
!        Compute contributions
         DO I=1,NOSIM
            IOFF = (I-1)*IDIM
            ISIMOFF = (I-1)*N2ORBX+1
            DO J=1,TNPATM
               JOFF = IOFF+(J-1)*3
               DIPORG(1) = NPCORD(1,J)
               DIPORG(2) = NPCORD(2,J)
               DIPORG(3) = NPCORD(3,J)
               CALL DZERO(WORK(KINTAO),3*NNBASX)
               CALL GET1IN(WORK(KINTAO),'NEFIELD',NOCOMP,WORK(KFREE),   &
     &                     LFREE,LABINT,INTREP,INTADR,J,TOFILE,KPATOM,  &
     &                     TRIMAT,DUMMY,EXP1VL,DUMMY,0)
!              X-component
               CALL DZERO(WORK(KTRMO),NNORBX)
               CALL DZERO(WORK(KUTR),N2ORBX)
               CALL UTHU(WORK(KINTAO),WORK(KTRMO),CMO,WORK(KFREE),      &
     &                   NBAST,NORBT)
               CALL DSPTSI(NORBT,WORK(KTRMO),WORK(KUTR))
               FACT = -FMQVEC(JOFF+1)
               CALL DAXPY(N2ORBX,FACT,WORK(KUTR),1,FXYVEC(ISIMOFF),1)
!              Y-component
               CALL DZERO(WORK(KTRMO),NNORBX)
               CALL DZERO(WORK(KUTR),N2ORBX)
               CALL UTHU(WORK(KINTAO+NNBASX),WORK(KTRMO),CMO,           &
     &                   WORK(KFREE),NBAST,NORBT)
               CALL DSPTSI(NORBT,WORK(KTRMO),WORK(KUTR))
               FACT = -FMQVEC(JOFF+2)
               CALL DAXPY(N2ORBX,FACT,WORK(KUTR),1,FXYVEC(ISIMOFF),1)
!              Z-component
               CALL DZERO(WORK(KTRMO),NNORBX)
               CALL DZERO(WORK(KUTR),N2ORBX)
               CALL UTHU(WORK(KINTAO+2*NNBASX),WORK(KTRMO),CMO,         &
     &                   WORK(KFREE),NBAST,NORBT)
               CALL DSPTSI(NORBT,WORK(KTRMO),WORK(KUTR))
               FACT = -FMQVEC(JOFF+3)
               CALL DAXPY(N2ORBX,FACT,WORK(KUTR),1,FXYVEC(ISIMOFF),1)
            END DO
         END DO
         RUNQM3 = .FALSE.
         CALL MEMREL('GET_XYVEC_LR',WORK,1,1,KFREE,LFREE)
      END IF
!     Induced dipole moment in NP region interaction with QM region
      IF (DONPCAP.AND.NOVDAMP) THEN
         ISTART = 0
         IF (DONPPOL) ISTART = 3*TNPATM
!        Fix me MM region shift
!        Electronic interaction part
         CALL MEMGET('REAL',KINTAO,NNBASX,WORK,KFREE,LFREE)
!        Allocate temp. matrices used for integrals tranformations
         CALL MEMGET('REAL',KTRMO,NNORBX,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KUTR,N2ORBX,WORK,KFREE,LFREE)
!        Set integrals evaluation flags
         KPATOM = 0
         NOCOMP = 1
         TOFILE = .FALSE.
         TRIMAT = .TRUE.
         EXP1VL = .FALSE.
         RUNQM3 = .TRUE.
!        Compute contributions
         DO I=1,NOSIM
            IOFF = (I-1)*IDIM
            ISIMOFF = (I-1)*N2ORBX+1
            DO J=1,TNPATM
               JOFF = IOFF+ISTART+J
               DIPORG(1) = NPCORD(1,J)
               DIPORG(2) = NPCORD(2,J)
               DIPORG(3) = NPCORD(3,J)
               CALL DZERO(WORK(KINTAO),NNBASX)
               CALL GET1IN(WORK(KINTAO),'NPETES ',NOCOMP,WORK(KFREE),   &
     &                     LFREE,LABINT,INTREP,INTADR,J,TOFILE,KPATOM,  &
     &                     TRIMAT,DUMMY,EXP1VL,DUMMY,0)
!              Zero temporary arrays
               CALL DZERO(WORK(KTRMO),NNORBX)
               CALL DZERO(WORK(KUTR),N2ORBX)
!              Transfor integrals
               CALL UTHU(WORK(KINTAO),WORK(KTRMO),CMO,WORK(KFREE),      &
     &                   NBAST,NORBT)
               CALL DSPTSI(NORBT,WORK(KTRMO),WORK(KUTR))
               FACT = FMQVEC(JOFF)
               CALL DAXPY(N2ORBX,FACT,WORK(KUTR),1,FXYVEC(ISIMOFF),1)
           END DO
         END DO
         RUNQM3 = .FALSE.
         CALL MEMREL('GET_XYVEC_LR',WORK,1,1,KFREE,LFREE)
      END IF
!     Restore origin coordinates
      DIPORG(1) = RSAVORG(1)
      DIPORG(2) = RSAVORG(2)
      DIPORG(3) = RSAVORG(3)
!
      RETURN
      end subroutine

!  /* Deck get_xyvec_qr */
      SUBROUTINE GET_XYVEC_QR(FMQVEC1,FMQVEC2,UCMO,IDIM,NSIM,FVEC,      &
     &                        ISYMT,ISYMV2,ZYM2,WORK,LWORK)
!
! Purpose:
!     Computes contribution to XY vector from induced dipoles moments
!     and charges.
!
! Input:
!   FMQVEC - Induced dipole moments and charges
!   CMO    - Molecular orbitals
!   WORK   - Temporary memory array
!   LWORK  - Size of temporary memory array
! Output:
!   FXYVEC - Contribution to XY vector from induced dipole moments and
!            charges
!   IDIM   - Size of electric field/potential vector
!   NOSIM  - Number of perturbed density matrices
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "qmnpmm.h"
#include "inforb.h"
#include "infdim.h"
#include "mxcent.h"
#include "iratdef.h"
#include "nuclei.h"
#include "orgcom.h"
#include "qm3.h"
#include "infrsp.h"
!
      DIMENSION FMQVEC1(NSIM*IDIM), FMQVEC2(NSIM*IDIM)
      DIMENSION UCMO(NORBT*NBAST), ZYM2(*)
      DIMENSION FVEC(NSIM*N2ORBX), WORK(LWORK)
!
      LOGICAL TOFILE,TRIMAT,EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      DIMENSION RSAVORG(3)
!
      KFREE = 1
      LFREE = LWORK
!     Save origin coordinates
      RSAVORG(1) = DIPORG(1)
      RSAVORG(2) = DIPORG(2)
      RSAVORG(3) = DIPORG(3)
!     Induced dipole moment in NP region interaction with QM region
      IF (DONPPOL.AND.NOVDAMP) THEN
!        Electronic interaction part
         CALL MEMGET('REAL',KINTAO,3*NNBASX,WORK,KFREE,LFREE)
!        Allocate temp. matrices used for integrals tranformations
         CALL MEMGET('REAL',KTRMO,NNORBX,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KUTR,N2ORBX,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KTLMA,N2ORBX,WORK,KFREE,LFREE)
!        Set integrals evaluation flags
         KPATOM = 0
         NOCOMP = 3
         TOFILE = .FALSE.
         TRIMAT = .TRUE.
         EXP1VL = .FALSE.
         RUNQM3 = .TRUE.
!        Compute contributions
         DO I=1,NSIM
            IOFF = (I-1)*IDIM
            ISIMOFF = (I-1)*N2ORBX+1
            DO J=1,TNPATM
              JOFF = IOFF+(J-1)*3
              DIPORG(1) = NPCORD(1,J)
              DIPORG(2) = NPCORD(2,J)
              DIPORG(3) = NPCORD(3,J)
              CALL DZERO(WORK(KINTAO),3*NNBASX)
              CALL GET1IN(WORK(KINTAO),'NEFIELD',NOCOMP,WORK(KFREE),    &
     &                    LFREE,LABINT,INTREP,INTADR,J,TOFILE,KPATOM,   &
     &                    TRIMAT,DUMMY,EXP1VL,DUMMY,0)
!             X-component
              CALL DZERO(WORK(KTRMO),NNORBX)
              CALL DZERO(WORK(KUTR),N2ORBX)
!             Transform integrals
              CALL UTHU(WORK(KINTAO),WORK(KTRMO),UCMO,WORK(KFREE),NBAST,&
     &                  NORBT)
              CALL DSPTSI(NORBT,WORK(KTRMO),WORK(KUTR))
!             Determine MM region contribution
              F1VAL = -FMQVEC1(JOFF+1)
              F2VAL = -0.50D0*FMQVEC2(JOFF+1)
              CALL DZERO(WORK(KTLMA),N2ORBX)
              CALL OITH1(ISYMV2,ZYM2,WORK(KUTR),WORK(KTLMA),ISYMT)
              CALL DAXPY(N2ORBX,F1VAL,WORK(KTLMA),1,FVEC,1)
              CALL DAXPY(N2ORBX,F2VAL,WORK(KTLMA),1,FVEC,1)
!             Y-component
              CALL DZERO(WORK(KTRMO),NNORBX)
              CALL DZERO(WORK(KUTR),N2ORBX)
!             Transform integrals
              CALL UTHU(WORK(KINTAO+NNBASX),WORK(KTRMO),UCMO,           &
     &                  WORK(KFREE),NBAST,NORBT)
              CALL DSPTSI(NORBT,WORK(KTRMO),WORK(KUTR))
!             Determine MM region contribution
              F1VAL = -FMQVEC1(JOFF+2)
              F2VAL = -0.50D0*FMQVEC2(JOFF+2)
              CALL DZERO(WORK(KTLMA),N2ORBX)
              CALL OITH1(ISYMV2,ZYM2,WORK(KUTR),WORK(KTLMA),ISYMT)
              CALL DAXPY(N2ORBX,F1VAL,WORK(KTLMA),1,FVEC,1)
              CALL DAXPY(N2ORBX,F2VAL,WORK(KTLMA),1,FVEC,1)
!             Z-component
              CALL DZERO(WORK(KTRMO),NNORBX)
              CALL DZERO(WORK(KUTR),N2ORBX)
!             Transform integrals
              CALL UTHU(WORK(KINTAO+2*NNBASX),WORK(KTRMO),UCMO,         &
     &                  WORK(KFREE),NBAST,NORBT)
              CALL DSPTSI(NORBT,WORK(KTRMO),WORK(KUTR))
!             Determine MM region contribution
              F1VAL = -FMQVEC1(JOFF+3)
              F2VAL = -0.50D0*FMQVEC2(JOFF+3)
              CALL DZERO(WORK(KTLMA),N2ORBX)
              CALL OITH1(ISYMV2,ZYM2,WORK(KUTR),WORK(KTLMA),ISYMT)
              CALL DAXPY(N2ORBX,F1VAL,WORK(KTLMA),1,FVEC,1)
              CALL DAXPY(N2ORBX,F2VAL,WORK(KTLMA),1,FVEC,1)
            END DO
         END DO
         RUNQM3 = .FALSE.
         CALL MEMREL('GET_XYVEC_QR',WORK,1,1,KFREE,LFREE)
      END IF
!     Induced dipole moment in NP region interaction with QM region
      IF (DONPCAP.AND.NOVDAMP) THEN
         ISTART = 0
         IF (DONPPOL) ISTART = 3*TNPATM
!        Fix me MM region shift
!        Electronic interaction part
         CALL MEMGET('REAL',KINTAO,NNBASX,WORK,KFREE,LFREE)
!        Allocate temp. matrices used for integrals tranformations
         CALL MEMGET('REAL',KTRMO,NNORBX,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KUTR,N2ORBX,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KTLMA,N2ORBX,WORK,KFREE,LFREE)
!        Set integrals evaluation flags
         KPATOM = 0
         NOCOMP = 1
         TOFILE = .FALSE.
         TRIMAT = .TRUE.
         EXP1VL = .FALSE.
         RUNQM3 = .TRUE.
!        Compute contributions
         DO I=1,NSIM
            IOFF = (I-1)*IDIM
            DO J=1,TNPATM
               JOFF = IOFF+ISTART+J
               DIPORG(1) = NPCORD(1,J)
               DIPORG(2) = NPCORD(2,J)
               DIPORG(3) = NPCORD(3,J)
               CALL DZERO(WORK(KINTAO),NNBASX)
               CALL GET1IN(WORK(KINTAO),'NPETES ',NOCOMP,WORK(KFREE),   &
     &                     LFREE,LABINT,INTREP,INTADR,J,TOFILE,KPATOM,  &
     &                     TRIMAT,DUMMY,EXP1VL,DUMMY,0)
!             Zero integral buffers
              CALL DZERO(WORK(KTRMO),NNORBX)
              CALL DZERO(WORK(KUTR),N2ORBX)
!             Transform integrals
              CALL UTHU(WORK(KINTAO),WORK(KTRMO),UCMO,WORK(KFREE),NBAST,&
     &                  NORBT)
              CALL DSPTSI(NORBT,WORK(KTRMO),WORK(KUTR))
!             Determine MM region contribution
              F1VAL = FMQVEC1(JOFF)
              F2VAL = 0.50D0*FMQVEC2(JOFF)
              CALL DZERO(WORK(KTLMA),N2ORBX)
              CALL OITH1(ISYMV2,ZYM2,WORK(KUTR),WORK(KTLMA),ISYMT)
              CALL DAXPY(N2ORBX,F1VAL,WORK(KTLMA),1,FVEC,1)
              CALL DAXPY(N2ORBX,F2VAL,WORK(KTLMA),1,FVEC,1)
           END DO
         END DO
         RUNQM3 = .FALSE.
         CALL MEMREL('GET_XYVEC_QR',WORK,1,1,KFREE,LFREE)
      END IF
!     Restore origin coordinates
      DIPORG(1) = RSAVORG(1)
      DIPORG(2) = RSAVORG(2)
      DIPORG(3) = RSAVORG(3)
!
      RETURN
      end subroutine

