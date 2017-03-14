!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2016 (2015), see http://daltonprogram.org"
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
C
      SUBROUTINE T4DRV(IBCDEQ,ISYMA,ISYMB,ISYMC,ISYMD,VECA,VECB,VECC,
     *                 VECD,FREQB,FREQC,FREQD,XINDX,UDV,PV,MJWOP,
     &                 WRK,LWRK,CMO,FC)
C
C     Purpose:
C     Compute E[4] times three vectors
C     Compute S[4] times three vectors, multiply with omega, and add to result
C
C     OUTPUT is put in WRK:
C     First  KZYVA elements are result vector
C
#include "implicit.h"
C
      PARAMETER ( D0 = 0.0D0 , DM1 = -1.0D0 , DEM10 = 1.0D-10 )
C
#include "maxorb.h"
#include "priunit.h"
#include "thrzer.h"
#include "infvar.h"
#include "qrinf.h"
#include "rspprp.h"
#include "infhyp.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "infpri.h"
#include "infspi.h"
#include "infcr.h"
C
      DIMENSION WRK(*)
      DIMENSION XINDX(*)
      DIMENSION UDV(*), MJWOP(2,MAXWOP,8)
      DIMENSION PV(*)
      DIMENSION VECA(*),VECB(*),VECC(*),VECD(*)
      DIMENSION CMO(*),FC(*)
C
C     Compute the length of storage needed
C     KZYVA is the length of the resulting vector
C
      KZYVA  = MZYVAR(ISYMA)
      KZYVB  = MZYVAR(ISYMB)
      KZYVC  = MZYVAR(ISYMC)
      KZYVD  = MZYVAR(ISYMD)
C
      KE4    = 1
      KS4    = KE4  + KZYVA
      KWRK   = KS4  + KZYVA
      LWRKE  = LWRK - KS4
      LWRKS  = LWRK - KWRK
      IF (LWRKS.LT.0) CALL ERRWRK('T4DRV',KWRK-1,LWRK)
C
C     Initialise the result vector
C
      CALL DZERO(WRK(KE4),KZYVA)
C
C
      IF( IPRRSP .GT. 100) THEN
         WRITE(LUPRI, '(//A)')
     *   ' Characteristics in T4DRV routine'
         WRITE(LUPRI, '( //A,3I10 )')
     *   ' Symmetry of vectors b, c and d: ', ISYMB,ISYMC,ISYMD
         WRITE(LUPRI, '( A,3I10 )')
     *   ' Length of orbital part        : ',
     *     MZWOPT(ISYMB),MZWOPT(ISYMC),MZWOPT(ISYMD)
         WRITE(LUPRI, '( A,3I10 )')
     *   ' Length of vectors             : ', KZYVB,KZYVC,KZYVD
         WRITE(LUPRI, '( A,3F10.7 )')
     *   ' Frequencies of vectors        : ', FREQB,FREQC,FREQD
         WRITE(LUPRI, '( A,I10 )')
     *   ' vector b equal to vector c or d ?  : ', IBCDEQ
      END IF
C
      IF ( IBCDEQ .EQ. 0 ) RETURN
C
C
C     Compute ( E4(jklm) + permutations ) N(k) N(l) N(m)
C
      CALL E4INIT(VECA,VECB,VECC,VECD,IBCDEQ,
     *            WRK(KE4),XINDX,UDV,PV,WRK(KS4),LWRKE,
     *            KZYVA,KZYVB,KZYVC,KZYVD,
     *            ISYMA,ISYMB,ISYMC,ISYMD,CMO,FC,MJWOP)
C
C
C     Compute ( S(jklm) + S(jkml) ) N(k) N(l) N(m)
C
      IF (FREQB .NE. D0 ) THEN
         CALL S4INIT(KZYVA,KZYVB,KZYVC,KZYVD,
     *               ISYMA,ISYMB,ISYMC,ISYMD,
     *               VECB,VECC,VECD,
     *               WRK(KS4),XINDX,UDV,MJWOP,WRK(KWRK),LWRKS)
C
C        Multiply with omega-b
C
         CALL DAXPY(KZYVA,FREQB,WRK(KS4),1,WRK(KE4),1)
      END IF
C
C     Compute ( S(jlkm) + S(jlmk) ) N(k) N(l) N(m)
C
      IF (FREQC .NE. D0 ) THEN
         CALL S4INIT(KZYVA,KZYVC,KZYVB,KZYVD,
     *               ISYMA,ISYMC,ISYMB,ISYMD,
     *               VECC,VECB,VECD,
     *               WRK(KS4),XINDX,UDV,MJWOP,WRK(KWRK),LWRKS)
C
C        Multiply with omega-c
C
         CALL DAXPY(KZYVA,FREQC,WRK(KS4),1,WRK(KE4),1)
      END IF
C
C     Compute ( S(jmkl) + S(jmlk) ) N(k) N(l) N(m)
C
      IF (FREQD .NE. D0 ) THEN
         CALL S4INIT(KZYVA,KZYVD,KZYVB,KZYVC,
     *               ISYMA,ISYMD,ISYMB,ISYMC,
     *               VECD,VECB,VECC,
     *               WRK(KS4),XINDX,UDV,MJWOP,WRK(KWRK),LWRKS)
C
C        Multiply with omega-d
C
         CALL DAXPY(KZYVA,FREQD,WRK(KS4),1,WRK(KE4),1)
      END IF
C
C
C     Result is now in WRK(KE4)
C
C
      IF (IPRRSP .GT. 100 ) THEN
         WRITE(LUPRI,'(//A)') ' Final result in T4DRV'
         WRITE(LUPRI,'(A)') ' ======================'
         CALL OUTPUT(WRK(KE4),1,KZYVA,1,1,KZYVA,1,1,LUPRI)
      END IF
C
      RETURN
      END
      SUBROUTINE E4INIT(VECA,VECB,VECC,VECD,IBCDEQ,
     *                  ETRS,XINDX,UDV,PV,WRK,LWRK,
     *                  KZYVA,KZYVB,KZYVC,KZYVD,
     *                  ISYMA,ISYMB,ISYMC,ISYMD,CMO,FC,MJWOP)
C
C
      use pelib_interface, only: use_pelib
#include "implicit.h"
C
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0 )
C
#include "priunit.h"
#include "infdim.h"
#include "maxorb.h"
#include "maxash.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "rspprp.h"
#include "infcr.h"
#include "infinp.h"
#include "infpri.h"
#include "dftcom.h"
#include "gnrinf.h"
C
      DIMENSION VECA(*),VECB(*),VECC(*),VECD(*)
      DIMENSION ETRS(KZYVA), MJWOP(2,MAXWOP,8)
      DIMENSION XINDX(*)
      DIMENSION UDV(*)
      DIMENSION PV(*)
      DIMENSION WRK(*)
      DIMENSION CMO(*),FC(*)

      LOGICAL   ADDFOCK, DFTSAV

C
C     Initialise variables and layout some workspace
C
      KFREE = 1
      LFREE = LWRK
      CALL MEMGET('REAL',KFI,NORBT*NORBT,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KZYMB,NORBT*NORBT,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KZYMC,NORBT*NORBT,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KZYMD,NORBT*NORBT,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDEN1,NASHT*NASHT,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDEN2,N2ASHX*N2ASHX,WRK,KFREE,LFREE)
C
C     Direct one-index transformation check
C
      IF (DIROIT.AND.TDHF) THEN
         DFTSAV = DFTADD
         DFTADD = .FALSE.
         CALL C4FOCK(IBCDEQ,WRK(KFI),VECB,VECC,VECD,
     *              WRK(KZYMB),WRK(KZYMC),WRK(KZYMD),
     *              KZYVB,KZYVC,KZYVD,
     *              ISYMA,ISYMB,ISYMC,ISYMD,
     *              CMO,FC,MJWOP,WRK(KFREE),LFREE)
         KDUM = KFREE
         CALL DZERO(WRK(KDUM),NORBT*NORBT)
         CALL ORBEX(ISYMA,1,KZYVA,WRK(KFI),WRK(KDUM),WRK(KDUM),
     *              WRK(KDUM),ETRS,D1,WRK(KDUM),MJWOP)
         VAL = DDOT(KZYVA,VECA,1,ETRS,1)
           WRITE (LUPRI,'(A,F20.12)') 'Contribution from C4FOCK:',
     *          VAL
         DFTADD = DFTSAV
      ELSE
         IF(DFTADD) CALL QUIT('DFT not implemented for this path.')
         CALL E4DRV(VECA,VECB,VECC,VECD,IBCDEQ,ETRS,XINDX,WRK(KFI),
     *       WRK(KZYMB),WRK(KZYMC),WRK(KZYMD),WRK(KDEN1),WRK(KDEN2),
     *              UDV,PV,WRK(KFREE),LFREE,KZYVA,KZYVB,KZYVC,KZYVD,
     *              ISYMA,ISYMB,ISYMC,ISYMD,CMO,FC,MJWOP)
      END IF

      IF (DFTADD) THEN
         IF (DOHFSRDFT .OR. DOMCSRDFT) THEN
            CALL QUIT('srDFT not implemented yet for cubic response')
         END IF
         CALL DZERO(WRK(KFI),N2ORBX)
         CALL DFT4DRV(VECB, VECC, VECD,
     &                   WRK(KFI),WRK(KZYMB),WRK(KZYMC),WRK(KZYMD),
     &                   KZYVB,KZYVC, KZYVD,
     &                   ISYMB,ISYMC, ISYMD,
     &                   ISPINB,ISPINC,ISPIND,
     &                   CMO,MJWOP,WRK(KFREE),LFREE)
           KDUM = KFREE
           CALL DZERO(WRK(KDUM),NORBT*NORBT)
           CALL ORBEX(ISYMA,1,KZYVA,WRK(KFI),WRK(KDUM),WRK(KDUM),
     *          WRK(KDUM),ETRS,D1,WRK(KDUM),MJWOP)
           WRITE (LUPRI,'(A,F20.1)') 'Contribution from DFTCRCS:',
     *          DDOT(KZYVA,VECA,1,ETRS,1) - VAL
      END IF
C
      IF (USE_PELIB()) THEN
        CALL QUIT('ERROR: cubic response not implemented'//
     &             ' for PEQM.')
      END IF
C
      IF (QMMM) THEN
         CALL QUIT('ERROR: cubic response not implemented'//
     &             ' for QMMM.')
      END IF

      IF (FLAG(16)) THEN
         CALL MEMGET('INTE',KSYRLM,(2*LSOLMX+1),WRK,KFREE,LFREE)
         VAL = DDOT(KZYVA,VECA,1,ETRS,1)
         IF (TDHF) THEN
            CALL C4SOL(VECA, VECB, VECC, VECD,
     *                 ETRS,XINDX,WRK(KZYMB),WRK(KZYMC),WRK(KZYMD),
     *                 UDV,WRK(KFREE),LFREE,
     *                 KZYVA,KZYVB,KZYVC,KZYVD,
     *                 ISYMA,ISYMB,ISYMC,ISYMD,CMO,MJWOP,WRK(KSYRLM))
            CALL C4SOL(VECA, VECB, VECD, VECC,
     *                 ETRS,XINDX,WRK(KZYMB),WRK(KZYMD),WRK(KZYMC),
     *                 UDV,WRK(KFREE),LFREE,
     *                 KZYVA,KZYVB,KZYVD,KZYVC,
     *                 ISYMA,ISYMB,ISYMD,ISYMC,CMO,MJWOP,WRK(KSYRLM))
            CALL C4SOL(VECA, VECC, VECB, VECD,
     *                 ETRS,XINDX,WRK(KZYMC),WRK(KZYMB),WRK(KZYMD),
     *                 UDV,WRK(KFREE),LFREE,
     *                 KZYVA,KZYVC,KZYVB,KZYVD,
     *                 ISYMA,ISYMC,ISYMB,ISYMD,CMO,MJWOP,WRK(KSYRLM))
            CALL C4SOL(VECA, VECC, VECD, VECB,
     *                 ETRS,XINDX,WRK(KZYMC),WRK(KZYMD),WRK(KZYMB),
     *                 UDV,WRK(KFREE),LFREE,
     *                 KZYVA,KZYVC,KZYVD,KZYVB,
     *                 ISYMA,ISYMC,ISYMD,ISYMB,CMO,MJWOP,WRK(KSYRLM))
            CALL C4SOL(VECA, VECD, VECC, VECB,
     *                 ETRS,XINDX,WRK(KZYMD),WRK(KZYMC),WRK(KZYMB),
     *                 UDV,WRK(KFREE),LFREE,
     *                 KZYVA,KZYVD,KZYVC,KZYVB,
     *                 ISYMA,ISYMD,ISYMC,ISYMB,CMO,MJWOP,WRK(KSYRLM))
            CALL C4SOL(VECA, VECD, VECB, VECC,
     *                 ETRS,XINDX,WRK(KZYMD),WRK(KZYMB),WRK(KZYMC),
     *                 UDV,WRK(KFREE),LFREE,
     *                 KZYVA,KZYVD,KZYVB,KZYVC,
     *                 ISYMA,ISYMD,ISYMB,ISYMC,CMO,MJWOP,WRK(KSYRLM))
         ELSE
            CALL E4SOL(VECA, VECB, VECC, VECD,
     *                 ETRS,XINDX,WRK(KZYMB),WRK(KZYMC),WRK(KZYMD),
     *                 WRK(KDEN1),UDV,WRK(KFREE),LFREE,
     *                 KZYVA,KZYVB,KZYVC,KZYVD,
     *                 ISYMA,ISYMB,ISYMC,ISYMD,CMO,MJWOP,WRK(KSYRLM))
         END IF
         WRITE (LUPRI,'(A,F20.12)') 'Total contribution from E4SOL:',
     *        DDOT(KZYVA,VECA,1,ETRS,1) - VAL
      END IF
C
C     overall minus sign for E[4] (missing in Eq. 34 jcp 105 (1996) 6401)
C
      CALL DSCAL(KZYVA,DM1,ETRS,1)
C
      RETURN
      END
      SUBROUTINE BCDCHK(DOCAL,IBCDEQ,LURSPRES,DIPLEN,GAMMA,
     *                  ISYMA,ISYMB,ISYMC,ISYMD,ISYMBC,ISYMBD,ISYMCD,
     *                  ALAB,BLAB,CLAB,DLAB,IAOP,IBOP,ICOP,IDOP,
     *                  IBFR,ICFR,IDFR,FREQA,FREQB,FREQC,FREQD,
     *                  KZYVA,KZYVB,KZYVC,KZYVD,KZYVBC,KZYVBD,KZYVCD,
     *                  VECA,VECB,VECC,VECD,VECBC,VECBD,VECCD,IEQTO)
C
C     called from CRHYP or CRTMO or CRTPA.
C
C     initialize variables.
C
C     check if an equivalent calculation already has been done,
C     DOCAL indicates the result. Save uniqe calculations in
C     LAB(N,4) and FRQ(N,4).
C
C     check if some of the response vectors are equal or zero,
C     IBCDEQ indicates the result:
C       0  at least one of VECB VECC VECD equal to zero
C       1  all of VECB, VECC, VECD unequal
C       2  VECB = VECC
C       3  VECB = VECD
C       4  VECC = VECD
C      24  VECB = VECC = VECD
C
C
#include "implicit.h"
C
      PARAMETER ( D0 = 0.0D0, THRSML = 1.0D-9 )
C
#include "thrzer.h"
#include "maxorb.h"
#include "priunit.h"
#include "infdim.h"
#include "inforb.h"
#include "wrkrsp.h"
#include "infrsp.h"
#include "infpri.h"
#include "infvar.h"
#include "qrinf.h"
#include "infspi.h"
#include "infden.h"
#include "rspprp.h"
#include "infcr.h"
#include "indcr.h"
#include "inftmo.h"
#include "inftpa.h"
C
      PARAMETER (THD = 1.0D-8, MXCALC = 200)
C
      LOGICAL DOCAL,DIPLEN,ONFIL
C
      CHARACTER*8 LAB, ALABEL, BLABEL, CLABEL, DLABEL, LABEL
      CHARACTER*8 ALAB,BLAB,CLAB,DLAB
C
      DIMENSION VECA(*), VECB(*),VECC(*),VECD(*),VECBC(*),
     &          VECBD(*),VECCD(*)
      DIMENSION LAB(MXCALC,4), FRQ(MXCALC,4), ISYM(MXCALC,4)
      DIMENSION GAMMA(3,3,3,3)
      DIMENSION IEQTO(5)
C
      SAVE LAB, FRQ, N, ISYM
C
      DATA N/0/
C
C     set variables
C
      ONFIL = .FALSE.
      ISYMBC = MULD2H(ISYMB,ISYMC)
      ISYMBD = MULD2H(ISYMB,ISYMD)
      ISYMCD = MULD2H(ISYMC,ISYMD)
C
      KZYVA  = MZYVAR(ISYMA)
      KZYVB  = MZYVAR(ISYMB)
      KZYVC  = MZYVAR(ISYMC)
      KZYVD  = MZYVAR(ISYMD)
      KZYVBC = MZYVAR(ISYMBC)
      KZYVBD = MZYVAR(ISYMBD)
      KZYVCD = MZYVAR(ISYMCD)
C
      IF (CRCAL) THEN
         FREQB = BCRFR(IBFR)
         FREQC = CCRFR(ICFR)
         FREQD = DCRFR(IDFR)
         FREQA = FREQB + FREQC + FREQD
C
         ALAB = ACRLB(ISYMA,IAOP)
         BLAB = BCRLB(ISYMB,IBOP)
         CLAB = CCRLB(ISYMC,ICOP)
         DLAB = DCRLB(ISYMD,IDOP)
      END IF
C
      IF (TOMOM) THEN
         FREQB = -BTMFR(IBFR)
         FREQC = -CTMFR(ICFR)
         FREQD = EXCIT2(ISYMD,IDFR)
         FREQA = FREQB + FREQC + FREQD
C
         ALAB = ATMLB(ISYMA,IAOP)
         BLAB = BTMLB(ISYMB,IBOP)
         CLAB = CTMLB(ISYMC,ICOP)
         DLAB = 'EXCITLAB'
      END IF
C
      IF (TPAMP) THEN
         FREQB = -BTPFR(IBFR)
         FREQC = -EXCIT2(ISYMC,ICFR)
         FREQD = EXCIT2(ISYMD,IDFR)
         FREQA = FREQB + FREQC + FREQD
C
         ALAB = ATPLB(ISYMA,IAOP)
         BLAB = BTPLB(ISYMB,IBOP)
         CLAB = 'EXCITLAB'
         DLAB = 'EXCITLAB'
      END IF
C
C If dipole operators only, store result in variable GAMMA
C
      DIPLEN = .FALSE.
      IF (ALAB(2:7) .EQ. 'DIPLEN' .AND.
     &    BLAB(2:7) .EQ. 'DIPLEN' .AND.
     &    CLAB(2:7) .EQ. 'DIPLEN' .AND.
     &    DLAB(2:7) .EQ. 'DIPLEN') THEN
          DIPLEN = .TRUE.
      END IF
C
C In case of gamma, only calculate contributions to average gamma
C
      IF (.NOT. GAMALL .AND. DIPLEN) THEN
         IF (.NOT.((ALAB.EQ.BLAB .OR. ALAB.EQ.CLAB .OR. ALAB.EQ.DLAB)
     &       .AND.(BLAB.EQ.ALAB .OR. BLAB.EQ.CLAB .OR. BLAB.EQ.DLAB)
     &       .AND.(CLAB.EQ.ALAB .OR. CLAB.EQ.BLAB .OR. CLAB.EQ.DLAB)
     &       .AND.(DLAB.EQ.ALAB .OR. DLAB.EQ.BLAB .OR. DLAB.EQ.CLAB)))
     &        THEN
            DOCAL = .FALSE.
            WRITE (LUPRI,'(/8A,/A)') ' Component ',ALAB,', ',BLAB,', ',
     &           CLAB,', ',DLAB,'   not calculated because it does'//
     &           ' not contribute to average gamma'
            RETURN
         END IF
      END IF
C
C     check if an equivalent calculation already has been done.
C
      DOCAL = .TRUE.
      DO 110 I = 1,N
      DO 120 J = 1,4
      DO 130 K = 1,4
      IF (K.NE.J) THEN
         DO 140 L = 1,4
         IF (L.NE.K .AND. L.NE.J) THEN
            DO 150 M = 1,4
            IF (M.NE.L .AND. M.NE.K .AND. M.NE.J) THEN
C
               IF ( ALAB.EQ.LAB(I,J) .AND.
     *              BLAB.EQ.LAB(I,K) .AND.
     *              CLAB.EQ.LAB(I,L) .AND.
     *              DLAB.EQ.LAB(I,M) .AND.
     *              ISYMA.EQ.ISYM(I,J) .AND.
     *              ISYMB.EQ.ISYM(I,K) .AND.
     *              ISYMC.EQ.ISYM(I,L) .AND.
     *              ISYMD.EQ.ISYM(I,M) .AND.
     *              ABS(-FREQA-FRQ(I,J)).LT.THD .AND.
     *              ABS( FREQB-FRQ(I,K)).LT.THD .AND.
     *              ABS( FREQC-FRQ(I,L)).LT.THD .AND.
     *              ABS( FREQD-FRQ(I,M)).LT.THD ) THEN
                   IF (DOCAL.AND.CRCAL) THEN
C
C Store the previously calculated result in gamma(A,B,C,D)
C
                      IF (DIPLEN) THEN
                         CALL DIPLAB(ALAB,I1)
                         CALL DIPLAB(BLAB,I2)
                         CALL DIPLAB(CLAB,I3)
                         CALL DIPLAB(DLAB,I4)
                         CALL DIPLAB(LAB(I,1),I5)
                         CALL DIPLAB(LAB(I,2),I6)
                         CALL DIPLAB(LAB(I,3),I7)
                         CALL DIPLAB(LAB(I,4),I8)
                         GAMMA(I1,I2,I3,I4)=GAMMA(I5,I6,I7,I8)
                      ELSE
      WRITE(LUPRI,'(//A,4(/A,A10,I4,F10.6))')
     *     ' Cubic response function value in a.u. for',
     *     ' A operator, symmetry, frequency: ',ALAB,ISYMA,-FREQA,
     *     ' B operator, symmetry, frequency: ',BLAB,ISYMB,FREQB,
     *     ' C operator, symmetry, frequency: ',CLAB,ISYMC,FREQC,
     *     ' D operator, symmetry, frequency: ',DLAB,ISYMD,FREQD
      WRITE(LUPRI,'(/A)') ' is equal to the already calculated'
      WRITE(LUPRI,'(/A,4(/A,A10,I4,F10.6))')
     *' Cubic response function value in a.u. for',
     *' A operator, symmetry, frequency: ',LAB(I,1),ISYM(I,1),FRQ(I,1),
     *' B operator, symmetry, frequency: ',LAB(I,2),ISYM(I,2),FRQ(I,2),
     *' C operator, symmetry, frequency: ',LAB(I,3),ISYM(I,3),FRQ(I,3),
     *' D operator, symmetry, frequency: ',LAB(I,4),ISYM(I,4),FRQ(I,4)
                      END IF
                   END IF
                   IF (DOCAL.AND.TOMOM) THEN
      WRITE(LUPRI,'(//A,3(/A,A10,I4,F10.6),/A,I7,F10.6)')
     *     ' Third order transition moment in a.u. for',
     *     ' A operator, symmetry, frequency: ',ALAB,ISYMA,FREQA,
     *     ' B operator, symmetry, frequency: ',BLAB,ISYMB,-FREQB,
     *     ' C operator, symmetry, frequency: ',CLAB,ISYMC,-FREQC,
     *     ' State no.,  symmetry, excitation energy:',ISYMD,FREQD
      WRITE(LUPRI,'(/A)') ' is equal to the already calculated'
      WRITE(LUPRI,'(/A,4(/A,A10,I4,F10.6))')
     *' Third order transition moment in a.u. for',
     *' A operator, symmetry, frequency: ',LAB(I,1),ISYM(I,1),-FRQ(I,1),
     *' B operator, symmetry, frequency: ',LAB(I,2),ISYM(I,2),-FRQ(I,2),
     *' C operator, symmetry, frequency: ',LAB(I,3),ISYM(I,3),-FRQ(I,3),
     *' Symmetry, excitation energy:     ','  ',ISYM(I,4),FRQ(I,4)
                      CALL DIPLAB(LAB(I,1),I1)
                      CALL DIPLAB(LAB(I,2),I2)
                      CALL DIPLAB(LAB(I,3),I3)
                      IEQTO(1)=I1
                      IEQTO(2)=I2
                      IEQTO(3)=I3
                      IEQTO(4)=IDFR
                      IEQTO(5)=ISYM(I,4)
                   END IF
                   IF (DOCAL.AND.TPAMP) THEN
      WRITE(LUPRI,'(//A,2(/A,A10,I4,F10.6),2(/A,I7,F10.6))')
     *     ' Third order transition moment in a.u. for',
     *     ' A operator, symmetry, frequency: ',ALAB,ISYMA,FREQA,
     *     ' B operator, symmetry, frequency: ',BLAB,ISYMB,-FREQB,
     *     ' State no.,  symmetry, excitation energy:',ISYMC,-FREQC,
     *     ' State no.,  symmetry, excitation energy:',ISYMD,FREQD
      WRITE(LUPRI,'(/A)') ' is equal to the already calculated'
      WRITE(LUPRI,'(/A,4(/A,A10,I4,F10.6))')
     *' Third order transition moment in a.u. for',
     *' A operator, symmetry, frequency: ',LAB(I,1),ISYM(I,1),-FRQ(I,1),
     *' B operator, symmetry, frequency: ',LAB(I,2),ISYM(I,2),-FRQ(I,2),
     *' Symmetry, excitation energy:     ','  ',ISYM(I,3),-FRQ(I,3),
     *' Symmetry, excitation energy:     ','  ',ISYM(I,4),FRQ(I,4)
C    This will set up IEQTO used for output analysis
                      CALL DIPLAB(LAB(I,1),I1)
                      CALL DIPLAB(LAB(I,2),I2)
                      CALL DIPLAB(LAB(I,3),I3)
                      IEQTO(1)=I1
                      IEQTO(2)=I2
                      IEQTO(3)=I3
                      IEQTO(4)=IDFR
                      IEQTO(5)=ISYM(I,4)
C    End of setting up IEQTO
                   END IF
                   DOCAL = .FALSE.
               END IF
C
            END IF
 150        CONTINUE
         END IF
 140     CONTINUE
      END IF
 130  CONTINUE
 120  CONTINUE
 110  CONTINUE
C
C     finally we check if the result is available from a previous
C     calculation
C     
      IF (DOCAL) THEN
         REWIND (LURSPRES)
 346     READ (LURSPRES,'(A8)',END=347) LABEL
         IF (LABEL(2:6) .EQ. 'Cubic') THEN
            READ (LURSPRES,'(36X,A8,I5,F10.3)') ALABEL,ISYMAF,FRQA
            READ (LURSPRES,'(36X,A8,I5,F10.3)') BLABEL,ISYMBF,FRQB
            READ (LURSPRES,'(36X,A8,I5,F10.3)') CLABEL,ISYMCF,FRQC
            READ (LURSPRES,'(36X,A8,I5,F10.3)') DLABEL,ISYMDF,FRQD
            READ (LURSPRES,'(A8)') LABEL
            READ (LURSPRES,'(21X,F20.8)') CUBVAL
         ELSE IF (LABEL(2:6) .EQ. 'Third') THEN
            READ (LURSPRES,'(36X,A8,I5,F10.3)') ALABEL,ISYMAF,FRQA
            READ (LURSPRES,'(36X,A8,I5,F10.3)') BLABEL,ISYMBF,FRQB
            READ (LURSPRES,'(36X,A8,I5,F10.3)') CLABEL,ISYMCF,FRQC
            IF (CLABEL .EQ. '        ') THEN
               CLABEL = 'EXCITLAB'
               READ (LURSPRES,'(41X,2I4,F10.3)') IC,ISYMCF,FRQC
            ELSE
               READ (LURSPRES,'(A8)') LABEL
            END IF
            DLABEL = 'EXCITLAB'
            READ (LURSPRES,'(41X,2I4,F10.3)') ID,ISYMDF,FRQD
            READ (LURSPRES,'(A8)') LABEL
            READ (LURSPRES,'(21X,F20.8)') CUBVAL
         ELSE
            GOTO 346
         END IF
         IF (CRCAL) THEN
            IF ((ALABEL .EQ. ALAB) .AND. (BLABEL .EQ. BLAB) .AND. 
     &           (CLABEL .EQ. CLAB) .AND. (DLABEL .EQ. DLAB) .AND.  
     &           (ISYMAF .EQ. ISYMA) .AND. (ISYMBF .EQ. ISYMB) .AND.
     &           (ISYMCF .EQ. ISYMC) .AND. (ISYMDF .EQ. ISYMD) .AND.
     &         ABS(FREQA+FRQA).LT.THD .AND. ABS(FREQB-FRQB).LT.THD .AND.
     &         ABS(FREQC-FRQC).LT.THD .AND. ABS(FREQD-FRQD).LT.THD) THEN
               WRITE(LUPRI,'(/A)') ' The following cubic response '//
     &              'function has already been calculated'  
               WRITE(LUPRI,'(/A,4(/A,A10,I4,F10.6))')
     *             ' Cubic response function value in a.u. for',
     *             ' A operator, symmetry, frequency: ',ALAB,ISYMA,FRQA,
     *             ' B operator, symmetry, frequency: ',BLAB,ISYMB,FRQB,
     *             ' C operator, symmetry, frequency: ',CLAB,ISYMC,FRQC,
     *             ' D operator, symmetry, frequency: ',DLAB,ISYMD,FRQD
               WRITE(LUPRI,'(/A,F20.8)') ' << A; B, C, D >>  = ', CUBVAL
               CALL FLSHFO(LUPRI)
               DOCAL = .FALSE.
               ONFIL = .TRUE.
               IF (DIPLEN) THEN
                  CALL DIPLAB(ALAB,I1)
                  CALL DIPLAB(BLAB,I2)
                  CALL DIPLAB(CLAB,I3)
                  CALL DIPLAB(DLAB,I4)
                  CALL DIPLAB(LAB(I,1),I5)
                  CALL DIPLAB(LAB(I,2),I6)
                  CALL DIPLAB(LAB(I,3),I7)
                  CALL DIPLAB(LAB(I,4),I8)
                  GAMMA(I1,I2,I3,I4)=-CUBVAL
               END IF
            END IF
         ELSE IF (TOMOM) THEN
            IF ((ALABEL .EQ. ALAB) .AND. (BLABEL .EQ. BLAB) .AND. 
     &           (CLABEL .EQ. CLAB) .AND. (DLABEL .EQ. DLAB) .AND.  
     &           (ISYMAF .EQ. ISYMA) .AND. (ISYMBF .EQ. ISYMB) .AND.
     &           (ISYMCF .EQ. ISYMC) .AND. (ISYMDF .EQ. ISYMD) .AND.
     &         ABS(FREQA-FRQA).LT.THD .AND. ABS(FREQB+FRQB).LT.THD .AND.
     &         ABS(FREQC+FRQC).LT.THD .AND. ABS(FREQD-FRQD).LT.THD) THEN
               WRITE(LUPRI,'(/A)') ' The following third order '//
     &              'transition moment has already been calculated'  
               WRITE(LUPRI,'(/A,3(/A,A10,I4,F10.6))')
     *              ' Third order transition moment in a.u. for',
     *           ' A operator, symmetry, frequency: ',ALAB,ISYMA,FREQA,
     *           ' B operator, symmetry, frequency: ',BLAB,ISYMB,-FREQB,
     *           ' C operator, symmetry, frequency: ',CLAB,ISYMC,-FREQC
               WRITE(LUPRI,'(/A,2I4,F10.6)')
     *        ' State no., symmetry, excitation energy:',ID,ISYMD,FREQD
               WRITE(LUPRI,'(/A,F20.8)') ' < 0 | ABC | f >  = ', CUBVAL
               CALL FLSHFO(LUPRI)
               DOCAL = .FALSE.
               ONFIL = .TRUE.
            END IF
         ELSE IF (TPAMP) THEN
            IF ((ALABEL .EQ. ALAB) .AND. (BLABEL .EQ. BLAB) .AND. 
     &           (CLABEL .EQ. CLAB) .AND. (DLABEL .EQ. DLAB) .AND.  
     &           (ISYMAF .EQ. ISYMA) .AND. (ISYMBF .EQ. ISYMB) .AND.
     &           (ISYMCF .EQ. ISYMC) .AND. (ISYMDF .EQ. ISYMD) .AND.
     &         ABS(FREQA-FRQA).LT.THD .AND. ABS(FREQB+FRQB).LT.THD .AND.
     &         ABS(FREQC+FRQC).LT.THD .AND. ABS(FREQD-FRQD).LT.THD) THEN
               WRITE(LUPRI,'(/A)') ' The following third order '//
     &              'transition moment has already been calculated'  
               WRITE(LUPRI,'(/A,2(/A,A10,I4,F10.6))')
     *            ' Third order transition moment in a.u. for',
     *            ' A operator, symmetry, frequency: ',ALAB,ISYMA,FREQA,
     *            ' B operator, symmetry, frequency: ',BLAB,ISYMB,-FREQB
               WRITE(LUPRI,'(2(/A,2I4,F10.6))')
     *       ' State no., symmetry, excitation energy:',IC,ISYMC,-FREQC,
     *       ' State no., symmetry, excitation energy:',ID,ISYMD,FREQD
               WRITE(LUPRI,'(/A,F20.8)') ' < e | AB | f >  = ', CUBVAL
               CALL FLSHFO(LUPRI)
               DOCAL = .FALSE.
               ONFIL = .TRUE.
            END IF
         END IF
         GOTO 346
 347     CONTINUE
#if defined (VAR_MFDS)
         BACKSPACE (LURSPRES)
#endif
      END IF

C
C     if not done before put in list with unique calculations
C
      IF ((DOCAL .OR. ONFIL).AND.(CRCAL.OR.TOMOM.OR.TPAMP)) THEN
         N = N + 1
         IF (N .GT. MXCALC) THEN
            WRITE (LUPRI,'(/A,I3)')
     &         'BCDCHK: # unique calculations .gt. MXCALC =',MXCALC
            CALL QUIT('BCDCHK: # unique calculations .gt. MXCALC')
         END IF
         LAB(N,1) = ALAB
         LAB(N,2) = BLAB
         LAB(N,3) = CLAB
         LAB(N,4) = DLAB
         ISYM(N,1) = ISYMA
         ISYM(N,2) = ISYMB
         ISYM(N,3) = ISYMC
         ISYM(N,4) = ISYMD
         FRQ(N,1) = -FREQA
         FRQ(N,2) = FREQB
         FRQ(N,3) = FREQC
         FRQ(N,4) = FREQD
      ELSE
         RETURN
      END IF
C
C
C     read in response vectors
C
      CALL READVE(ISYMA,ISYMB,ISYMC,ISYMD,
     *            ALAB,BLAB,CLAB,DLAB,
     *            FREQA,FREQB,FREQC,FREQD,
     *            KZYVA,KZYVB,KZYVC,KZYVD,KZYVBC,KZYVBD,KZYVCD,
     *            VECA,VECB,VECC,VECD,VECBC,VECBD,VECCD)
C
C     check if response vectors are zero or equal.
C
      IBCDEQ = 1
      IF ( ISYMB .EQ. ISYMC .AND. FREQB .EQ. FREQC ) THEN
         VMIN = D0
         DO 210 I = 1,KZYVB
            VMIN = VMIN + ABS(VECB(I) - VECC(I))
  210    CONTINUE
         IF (VMIN .LE. THRZER) IBCDEQ = 2
      END IF
C
      IF ( ISYMB .EQ. ISYMD .AND. FREQB .EQ. FREQD ) THEN
         VMIN = D0
         DO 230 I = 1,KZYVB
            VMIN = VMIN + ABS(VECB(I) - VECD(I))
  230    CONTINUE
         IF (VMIN .LE. THRZER) IBCDEQ = IBCDEQ * 3
      END IF
C
      IF ( ISYMC .EQ. ISYMD .AND. FREQD .EQ. FREQC ) THEN
         VMIN = D0
         DO 250 I = 1,KZYVC
            VMIN = VMIN + ABS(VECC(I) - VECD(I))
  250    CONTINUE
         IF (VMIN .LE. THRZER) IBCDEQ = IBCDEQ * 4
      END IF
C
      VNORM = DNRM2(KZYVA,VECA,1)
      IF ( VNORM .LT. THRSML ) IBCDEQ = 0
C
      VNORM = DNRM2(KZYVB,VECB,1)
      IF ( VNORM .LT. THRSML ) IBCDEQ = 0
C
      VNORM = DNRM2(KZYVC,VECC,1)
      IF ( VNORM .LT. THRSML ) IBCDEQ = 0
C
      VNORM = DNRM2(KZYVD,VECD,1)
      IF ( VNORM .LT. THRSML ) IBCDEQ = 0
C
      RETURN
      END

      SUBROUTINE DIPLAB(LAB,I)
C
#include "implicit.h"
C
      CHARACTER*8 LAB
C
C Map dipole operators to integer number
C X => 1, Y => 2, and Z => 3
C
      IF (LAB(1:1).EQ.'X') I=1
      IF (LAB(1:1).EQ.'Y') I=2
      IF (LAB(1:1).EQ.'Z') I=3
C
      RETURN
      END

      SUBROUTINE DFT4DRV(VECB, VECC, VECD,
     &                   FI,ZYMB,ZYMC,ZYMD,
     &                   KZYVB,KZYVC, KZYVD,
     &                   ISYMB,ISYMC, ISYMD,
     &                   ISPINB,ISPINC,ISPIND,
     &                   CMO,MJWOP,WRK,LFREE)
#include "implicit.h"
#include "infvar.h"
#include "inforb.h"
#include "maxorb.h"
#include "infinp.h"
#include "dftcom.h"
 
      DIMENSION VECB(*),VECC(*),VECD(*)
      DIMENSION ZYMB(*),ZYMC(*),ZYMD(*)
      DIMENSION WRK(*),CMO(*)
      DIMENSION FI(NORBT,NORBT), MJWOP(2,MAXWOP,8)
      LOGICAL ADDFOCK
c     unpack response vectors.
      NSIM = 1
      CALL GTZYMT(NSIM,VECB,KZYVB,ISYMB,ZYMB,MJWOP)
      CALL GTZYMT(NSIM,VECC,KZYVC,ISYMC,ZYMC,MJWOP)
      CALL GTZYMT(NSIM,VECD,KZYVD,ISYMD,ZYMD,MJWOP)
C     ADDFOCK = .NOT.DIRFCK
      ADDFOCK = .FALSE.
c     compute dft contribution
      CALL DFTCRCF(FI,CMO,
     &             ZYMB,ISYMB,
     &             ZYMC,ISYMC,
     &             ZYMD,ISYMD,
     &             WRK,LFREE,IPRDFT)
c     pack stuff back.
      END
