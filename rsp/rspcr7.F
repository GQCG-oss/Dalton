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
      SUBROUTINE S4INIT(KZYVA,KZYVB,KZYVC,KZYVD,
     *                  ISYMA,ISYMB,ISYMC,ISYMD,
     *                  VECB,VECC,VECD,
     *                  S4TRS,XINDX,UDV,MJWOP,WRK,LWRK)
C
C     Layout the core for the calculation of S4 times three vectors
C
#include <implicit.h>
#include <infdim.h>
#include <inforb.h>
#include <maxorb.h>
#include <maxash.h>
#include <infvar.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <rspprp.h>
#include <infhyp.h>
#include <qrinf.h>
#include <infpri.h>
#include <indcr.h>
#include <infcr.h>
C
      DIMENSION WRK(*)
      DIMENSION S4TRS(KZYVA), MJWOP(2,MAXWOP,8)
      DIMENSION VECB(KZYVB), VECC(KZYVC), VECD(KZYVD)
      DIMENSION XINDX(*)
      DIMENSION UDV(NASHDI,NASHDI)
C
C     Initialise the gradient to zero.
C
      CALL DZERO(S4TRS,KZYVA)
C
      KZYMAT  = 1
      KDEN1   = KZYMAT + NORBT * NORBT
      KFREE   = KDEN1  + NASHT * NASHT
      LWRKF   = LWRK - KFREE
      IF (LWRKF.LT.0) CALL ERRWRK('S4INIT',KFREE-1,LWRK)
C
      CALL S4DRV(KZYVA,KZYVB,KZYVC,KZYVD,
     *           ISYMA,ISYMB,ISYMC,ISYMD,
     *           S4TRS,VECB,VECC,VECD,UDV,
     *           WRK(KZYMAT),WRK(KDEN1),XINDX,MJWOP,WRK(KFREE),LWRKF)
C
      RETURN
      END
      SUBROUTINE TRZYM2(VEC1,VEC2,VEC3,KZYV1,KZYV2,KZYV3,
     *                  ISYM1,ISYM2,ISYM3,ZYMAT,MJWOP,
     *                  WRK,LWRK)
C
C     This subroutine unpacks the ZY matrices from the three vectors and
C     does the transformation
C
#include <implicit.h>
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <infvar.h>
#include <inforb.h>
#include <infind.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <rspprp.h>
#include <qrinf.h>
#include <infpri.h>
#include <infspi.h>
#include <indcr.h>
#include <infcr.h>
C
      DIMENSION VEC1(KZYV1), VEC2(KZYV2), VEC3(KZYV3), MJWOP(2,MAXWOP,8)
      DIMENSION ZYMAT(NORBT,NORBT)
      DIMENSION WRK(*)
C
C     Layout workspace
C
      KZY1  = 1
      KZY2  = KZY1  + NORBT * NORBT
      KZY3  = KZY2  + NORBT * NORBT
      KTEMP = KZY3  + NORBT * NORBT
      KFREE = KTEMP + NORBT * NORBT
      LFREE = LWRK  - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('TRZYMT',KFREE-1,LWRK)
C
C     Unpack the kappa(1), kappa(2), and kappa(3)
C
      CALL GTZYMT(1,VEC1,KZYV1,ISYM1,WRK(KZY1),MJWOP )
      CALL GTZYMT(1,VEC2,KZYV2,ISYM2,WRK(KZY2),MJWOP )
      CALL GTZYMT(1,VEC3,KZYV3,ISYM3,WRK(KZY3),MJWOP )
C
C     Calculate the commutator [k3,[k2,k1]]
C
         CALL DGEMM('N','N',NORBT,NORBT,NORBT,1.D0,
     &              WRK(KZY2),NORBT,
     &              WRK(KZY1),NORBT,0.D0,
     &              WRK(KTEMP),NORBT)
C
         CALL DGEMM('N','N',NORBT,NORBT,NORBT,-1.D0,
     &              WRK(KZY1),NORBT,
     &              WRK(KZY2),NORBT,1.D0,
     &              WRK(KTEMP),NORBT)
C
         CALL DGEMM('N','N',NORBT,NORBT,NORBT,1.D0,
     &              WRK(KZY3),NORBT,
     &              WRK(KTEMP),NORBT,0.D0,
     &              ZYMAT,NORBT)
C
         CALL DGEMM('N','N',NORBT,NORBT,NORBT,-1.D0,
     &              WRK(KTEMP),NORBT,
     &              WRK(KZY3),NORBT,1.D0,
     &              ZYMAT,NORBT)
C
         IF( IPRRSP .GT. 100 ) THEN
            WRITE(LUPRI,'(/A)') ' Final result in TRZYM2'
            WRITE(LUPRI,'(A)')  ' ======================'
            CALL OUTPUT(ZYMAT,1,NORBT,1,NORBT,NORBT,NORBT,
     *                  1,LUPRI)
         END IF
C
      RETURN
      END
      SUBROUTINE GETNXY(OMEGA,RESVEC,CMO,UDV,PV,FC,FV,FCAC,H2AC,
     *                  XINDX,WRK,LWRK)
C
#include <implicit.h>
C
C  THIS IS A TEST ROUTINE THAT
C  CALCULATES THE E(2) AND S(2) MATRICES EXPLICITLY BY CARRYING
C  OUT LINEAR TRANSFORMATIONS ON UNIT VECTORS
C  CONSTRUCTS E(2) - w*S(2) AND COMPUTES THE INVERSE
C
#include <priunit.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <infpri.h>
C
      PARAMETER ( D0 = 0.0D0 , D1 = 1.0D0 )
      DIMENSION RESVEC(*)
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
C
C  ALLOCATE WORK SPACE
C
      KE2   = 1
      KS2   = KE2    + KZYVAR*KZYVAR
      KBVEC = KS2    + KZYVAR*KZYVAR
      KWRK1 = KBVEC  + KZYVAR
      IF (KSYMOP .EQ. 1) THEN
         KCREF = KWRK1
         KWRK1 = KCREF  + KZCONF
      ELSE
         KCREF = -999 999 999
      END IF
      LWRK1 = LWRK   - KWRK1
      IF (LWRK1.LT.0) CALL ERRWRK('RSPES2',KWRK1-1,LWRK)
C
      IF (KSYMOP .EQ. 1 .AND. KZCONF.GT.0) THEN
         CALL GETREF(WRK(KCREF),KZCONF)
      END IF
      CALL DZERO(WRK(KBVEC),KZYVAR)
      DO 100 I = 1,KZYVAR
         IF (I.LE.KZCONF) THEN
            NCSIM = 1
            NOSIM = 0
            IOFF  = I
         ELSE IF (I.LE.KZVAR) THEN
            NCSIM = 0
            NOSIM = 1
            IOFF  = I - KZCONF
         ELSE IF (I.LE.KZVAR+KZCONF) THEN
            GO TO 100
         ELSE
            NCSIM = 0
            NOSIM = 1
            IOFF  = I - KZVAR - KZCONF + KZWOPT
         ENDIF
         WRK(KBVEC-1+IOFF) = D1
         IF (( NCSIM.GT.0 ).AND.( KSYMOP.EQ.1)) THEN
            IBOFF = 0
            ICREF  = IOFF
            IF (IOFF.GT.KZCONF) THEN
               IBOFF = KZCONF
               ICREF = ICREF - KZCONF
            END IF
            CALL DAXPY(KZCONF,-WRK(KCREF+ICREF-1),WRK(KCREF),1,
     *                 WRK(KBVEC+IBOFF),1)
         END IF
         IF (IPRRSP.GT.110) THEN
            IF(NOSIM.GT.0)  THEN
               KZYDIM = KZYWOP
               WRITE(LUPRI,'(/A)')' ORBITAL TRIAL VECTOR'
            END IF
            IF(NCSIM.GT.0) THEN
               KZYDIM = KZCONF
               WRITE(LUPRI,'(/A)')' CONFIGURATION TRIAL VECTOR'
            END IF
            CALL OUTPUT(WRK(KBVEC),1,KZYDIM,1,1,KZYDIM,1,1,LUPRI)
         END IF
         CALL RSPLIN(NCSIM,NOSIM,WRK(KBVEC),WRK(KBVEC),
     *               CMO,UDV,PV,FC,FV,FCAC,H2AC,
     *               XINDX,WRK(KWRK1),LWRK1)
C
C        CALL RSPLIN(NCSIM,NOSIM,ZYCVEC,ZYOVEC,
C    *               CMO,UDV,PV,FC,FV,FCAC,H2AC,
C    *               XINDX,WRK,LWRK)
C
C PROJECT OUT RERERENCE STATE COMPONENTS FROM LINEAR TRANSFORMED
C E2 AND S2 VECTORS
C
         IF ((.NOT.TDHF ).AND.( KSYMOP.EQ.1)) THEN
            E2OVL = DDOT(KZCONF,WRK(KCREF),1,WRK(KWRK1),1)
            CALL DAXPY(KZCONF,-E2OVL,WRK(KCREF),1,WRK(KWRK1),1)
            S2OVL = DDOT(KZCONF,WRK(KCREF),1,WRK(KWRK1+KZYVAR),1)
            CALL DAXPY(KZCONF,-S2OVL,WRK(KCREF),1,WRK(KWRK1+KZYVAR),1)
            E2OVL = DDOT(KZCONF,WRK(KCREF),1,WRK(KWRK1+KZVAR),1)
            CALL DAXPY(KZCONF,-E2OVL,WRK(KCREF),1,WRK(KWRK1+KZVAR),1)
            S2OVL = DDOT(KZCONF,WRK(KCREF),1,
     *                   WRK(KWRK1+KZYVAR+KZVAR),1)
            CALL DAXPY(KZCONF,-S2OVL,WRK(KCREF),1,
     *                 WRK(KWRK1+KZYVAR+KZVAR),1)
         END IF
         CALL DCOPY(KZYVAR,WRK(KWRK1),1,WRK(KE2+(I-1)*KZYVAR),1)
         CALL DCOPY(KZYVAR,WRK(KWRK1+KZYVAR),1,
     *              WRK(KS2+(I-1)*KZYVAR),1)
         IF (( NCSIM.GT.0 ).AND.( KSYMOP.EQ.1)) THEN
            CALL DZERO(WRK(KBVEC),KZCONF)
         ELSE
            WRK(KBVEC-1+IOFF) = D0
         END IF
 100  CONTINUE
      KTOT = KZVAR * KZYVAR
      DO 105 J = 1,KZCONF
         JTOT = (J-1)*KZYVAR
         CALL DCOPY(KZVAR,WRK(KE2+JTOT),1,WRK(KE2+JTOT+KTOT+KZVAR),1)
         CALL DCOPY(KZVAR,WRK(KE2+JTOT+KZVAR),1,WRK(KE2+JTOT+KTOT),1)
         CALL DCOPY(KZVAR,WRK(KS2+JTOT),1,WRK(KS2+JTOT+KTOT+KZVAR),1)
         CALL DCOPY(KZVAR,WRK(KS2+JTOT+KZVAR),1,WRK(KS2+JTOT+KTOT),1)
         CALL DSCAL(KZYVAR,-1.0D0,WRK(KS2+JTOT+KTOT),1)
 105  CONTINUE
      IF (IPRRSP.GT.100) THEN
         WRITE(LUPRI,'(A,I8)')' E(2) MATRIX : DIMENSION ',KZYVAR
         CALL OUTPUT(WRK(KE2),1,KZYVAR,1,KZYVAR,KZYVAR,KZYVAR,1,LUPRI)
         WRITE(LUPRI,'(A,I8)')' S(2) MATRIX : DIMENSION ',KZYVAR
         CALL OUTPUT(WRK(KS2),1,KZYVAR,1,KZYVAR,KZYVAR,KZYVAR,1,LUPRI)
      END IF
C
C     Compute inv[ E(2) - w*S(2) ] * XYvec
C
      CALL DAXPY(KZYVAR*KZYVAR,-OMEGA,WRK(KS2),1,WRK(KE2),1)
      CALL DGESOL(1,KZYVAR,KZYVAR,KZYVAR,WRK(KE2),RESVEC,WRK(KS2),INFO)
      IF (INFO.NE.0) WRITE (LUPRI,'(/A)')'**** ERROR IN GETNXY '
      RETURN
      END

