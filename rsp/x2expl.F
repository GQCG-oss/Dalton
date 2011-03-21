!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
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
#include <single.h>
C
#ifdef REV_LOG
===========================================================================
Revision 1.2  2000/05/24 19:13:57  hjj
new getref calls with appropriate NDREF instead of NCREF
(fixing error for triplet with CSF)
===========================================================================
#endif
C  /* Deck x2expl */
      SUBROUTINE X2EXPL(NSIM,MZYVA,MZYVC,ISYMA,ISPINA,
     *               ISYMC,ISPINC,DUMMY,CVEC,X2N,XINDX,
     *               UDV,PV,OPLBL,ISYMB,ISPINB,CMO,MJWOP,
     *               WRK,LWRK)
C
#include <implicit.h>
C
C  THIS IS A TEST ROUTINE THAT
C  1)  CALCULATE THE X(2) MATRIX EXPLICITLY BY CARRYING
C      OUT LINEAR TRANSFORMATIONS ON UNIT VECTORS
C
#include <priunit.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <infpri.h>
#include <maxorb.h>
#include <infvar.h>
#include <qrinf.h>
C
      PARAMETER ( D0 = 0.0D0 , D1 = 1.0D0 , DTOL = 1.0D-8 )
      DIMENSION CVEC(MZYVC), X2N(MZYVA), XINDX(*), UDV(*), PV(*)
      DIMENSION MJWOP(2,MAXWOP,8)
      DIMENSION CMO(*), WRK(*)
C
C  ALLOCATE WORK SPACE
C
      KX2 = 1
      KCREF = KX2 + MZYVA*MZYVC
      KWRK1 = KCREF  + MZCONF(1)
      LWRK1 = LWRK   - KWRK1
      IF (LWRK1.LT.0) CALL ERRWRK('X2EXPL',KWRK1-1,LWRK)
C
      IF (ISYMA.EQ.1 .AND. MZCONF(ISYMA).GT.0
     *                .OR.
     *    ISYMC.EQ.1 .AND. MZCONF(ISYMC).GT.0) THEN
         CALL GETREF(WRK(KCREF),MZCONF(1))
      END IF
      CALL DZERO(CVEC,MZYVAR(ISYMC))
      DO 100 I = 1,MZYVAR(ISYMC)
         IOFF = I
         CVEC(IOFF) = D1
         IF ( ISYMC.EQ.1 ) THEN
            IF (IOFF.LE.MZCONF(ISYMC)) THEN
               ICOFF = 0
               ICREF  = IOFF
            END IF
            IF (IOFF.GT.MZVAR(ISYMC)
     *           .AND. IOFF.LE.(MZVAR(ISYMC)+MZCONF(ISYMC))) THEN
               ICOFF = MZVAR(ISYMC)
               ICREF = IOFF - MZVAR(ISYMC)
            END IF
            CALL DAXPY(MZCONF(ISYMC),-WRK(KCREF+ICREF-1),WRK(KCREF),1,
     *                 CVEC(1+ICOFF),1)
         END IF
         IF (IPRRSP.GT.110) THEN
            KZYDIM = MZYVAR(ISYMC)
            WRITE(LUPRI,'(/A)') 'X2EXPL: TRIAL VECTOR'
            CALL OUTPUT(CVEC,1,KZYDIM/2,1,2,KZYDIM,2,1,LUPRI)
         END IF
         CALL X2INIT(1,MZYVAR(ISYMA),MZYVAR(ISYMC),ISYMA,ISPINA,
     *               ISYMC,ISPINC,DUMMY,CVEC,WRK(KWRK1),XINDX,
     *               UDV,PV,OPLBL,ISYMB,ISPINB,CMO,MJWOP,
     *               WRK(KWRK1+MZYVAR(ISYMA)),LWRK1-MZYVAR(ISYMA))
C     CALL X2INIT(NSIM,KZYVR,KZYV,IGRSYM,IGRSPI,ISYMV,ISPINV,
C    *                 GRVEC,VEC,X2TRS,XINDX,UDV,PV,OPLBL,IOPSYM,IOPSPI,
C    *                 CMO,MJWOP,WRK,LWRK)
C
C PROJECT OUT RERERENCE STATE COMPONENTS FROM LINEAR TRANSFORMED
C X2
C
         IF ((.NOT.TDHF ).AND.(ISYMA.EQ.1)) THEN
            X2 = DDOT(MZCONF(ISYMA),WRK(KCREF),1,WRK(KWRK1),1)
            CALL DAXPY(MZCONF(ISYMA),-X2,WRK(KCREF),1,WRK(KWRK1),1)
            X2 = DDOT(MZCONF(ISYMA),WRK(KCREF),1,
     *           WRK(KWRK1+MZVAR(ISYMA)),1)
            CALL DAXPY(MZCONF(ISYMA),-X2,WRK(KCREF),1,
     *           WRK(KWRK1+MZVAR(ISYMA)),1)
         END IF
         CALL DCOPY(MZYVAR(ISYMA),WRK(KWRK1),1,
     *              WRK(KX2+(I-1)*MZYVAR(ISYMA)),1)
         IF (ISYMA.EQ.1) THEN
            CALL DZERO(CVEC,MZYVAR(ISYMC))
         ELSE
            CVEC(IOFF) = D0
         END IF
 100  CONTINUE
      WRITE(LUPRI,'(A,2I8)')' X(2) MATRIX : DIMENSION ',
     *    MZYVAR(ISYMA),MZYVAR(ISYMC)
      CALL OUTPUT(WRK(KX2),1,MZYVAR(ISYMA),1,MZYVAR(ISYMC),
     *            MZYVAR(ISYMA),MZYVAR(ISYMC),1,LUPRI)
      RETURN
#ifdef DEBUGOUTPUT
      ZAMAX = D0
      IA    = 0
      JA    = 0
      ZBMAX = D0
      IB    = 0
      JB    = 0
      ZSMAX = D0
      IS    = 0
      JS    = 0
      ZDMAX = D0
      ID    = 0
      JD    = 0
      DO 150 IZ = 1,KZVAR
         IY = IZ +KZVAR
         DO 160 JZ = 1,IZ
            JY = JZ + KZVAR
            ZAIJ  = WRK(KX2-1+(IZ-1)*KZYVAR+JZ)
            ZAJI  = WRK(KX2-1+(JZ-1)*KZYVAR+IZ)
            YAIJ  = WRK(KX2-1+(IY-1)*KZYVAR+JY)
            YAJI  = WRK(KX2-1+(JY-1)*KZYVAR+IY)
            ZDEV  = MAX(ABS(ZAIJ-ZAJI),ABS(ZAIJ-YAIJ))
            ZDEV  = MAX(ZDEV,ABS(ZAIJ-YAJI))
            IF (ZDEV.GT.ZAMAX) THEN
               ZAMAX = ZDEV
               IA    = IZ
               JA    = JZ
            END IF
            ZBIJ  = WRK(KX2-1+(IZ-1)*KZYVAR+JY)
            ZBJI  = WRK(KX2-1+(JZ-1)*KZYVAR+IY)
            YBIJ  = WRK(KX2-1+(IY-1)*KZYVAR+JZ)
            YBJI  = WRK(KX2-1+(JY-1)*KZYVAR+IZ)
            ZDEV  = MAX(ABS(ZBIJ-ZBJI),ABS(ZBIJ-YBIJ))
            ZDEV  = MAX(ZDEV,ABS(ZBIJ-YBJI))
            IF (ZDEV.GT.ZBMAX) THEN
               ZBMAX = ZDEV
               IB    = IZ
               JB    = JZ
            END IF
 160     CONTINUE
 150  CONTINUE
      IZ=IA
      JZ=JA
      IY=IZ+KZVAR
      JY=JZ+KZVAR
      WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *' A(I,J)-MATRIX : MAXIMUM DEVIATION, I=',IA,' J=',JA,' DEV=',ZAMAX
      IF ( ZAMAX.GT.DTOL)
     *WRITE(LUPRI,'(/4(/A,I5,A,I5,A,1P,G16.8))')
     *' IZ=',IZ,' JZ=',JZ,' X2(IZ,JZ)',WRK(KX2-1+(IZ-1)*KZYVAR+JZ),
     *' IZ=',JZ,' JZ=',IZ,' X2(IZ,JZ)',WRK(KX2-1+(JZ-1)*KZYVAR+IZ),
     *' IZ=',IY,' JZ=',JY,' X2(IZ,JZ)',WRK(KX2-1+(IY-1)*KZYVAR+JY),
     *' IZ=',JY,' JZ=',IY,' X2(IZ,JZ)',WRK(KX2-1+(JY-1)*KZYVAR+IY)
      IZ=IB
      JZ=JB
      IY=IZ+KZVAR
      JY=JZ+KZVAR
      WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *' B(I,J)-MATRIX : MAXIMUM DEVIATION, I=',IB,' J=',JB,' DEV=',ZBMAX
      IF (ZBMAX.GT.DTOL)
     *WRITE(LUPRI,'(/4(/A,I5,A,I5,A,1P,G16.8))')
     *' IZ=',IZ,' JZ=',JY,' X2(IZ,JZ)',WRK(KX2-1+(IZ-1)*KZYVAR+JY),
     *' IZ=',JZ,' JZ=',IY,' X2(IZ,JZ)',WRK(KX2-1+(JZ-1)*KZYVAR+IY),
     *' IZ=',IY,' JZ=',JZ,' X2(IZ,JZ)',WRK(KX2-1+(IY-1)*KZYVAR+JZ),
     *' IZ=',JY,' JZ=',IZ,' X2(IZ,JZ)',WRK(KX2-1+(JY-1)*KZYVAR+IZ)
      RETURN
#endif
      END
