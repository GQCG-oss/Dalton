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
C  /* Deck eriexp */
      SUBROUTINE ERIEXP(EXPAB,EXPCD,ALPHA,FACPQ,HEXPP,HEXPQ,
     &                  RODPE1,RODPE2,IPNTPP,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxorb.h>
#include <r12int.h>
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, DP5 = 0.5D0)
      INTEGER AB, CD, PQ, RST
      DIMENSION EXPAB(NPQBCS,MLTPZ,NPRFAB,3),
     &          EXPCD(NPQBCS,MLTPZ,NPRFCD,3),
     &          ALPHA(NPQBCS,MLTPZ,NPRFPQ),
     &          HEXPP(NPQBCS,MLTPZ,NPRFPQ),
     &          HEXPQ(NPQBCS,MLTPZ,NPRFPQ),
     &          FACPQ(NPQBCS,NPRFPQ),
     &          RODPE1(NODPP1,  2), RODPE2(NODPP2,  2),
     &          IPNTPP(MAXBCH,2)
#include <cbieri.h>
#include <odclss.h>
#include <ericom.h>
#include <hertop.h>
#include <aobtch.h>
#include <symmet.h>
#include <ibtfun.h>
C
      IF (IPRINT .GT. 10) THEN
         CALL TITLER('Output from ERIEXP','*',103)
      END IF
C
      SFAC = SQRT(DP5)*FMULT(
     &       IBTAND(IBTAND(ISTBLA,ISTBLB),IBTAND(ISTBLC,ISTBLD)))
C
      DO 100 I = 1, NPRFAB
      DO 100 PQ = 1, NPQBCS
         AB = IPNTPP(PQ,1) + I - 1
         EXPA = RODPE1(AB,1)
         EXPB = RODPE1(AB,2)
         PINV = D1/(EXPA + EXPB)
         EXPAP = EXPA*PINV
         EXPBP = EXPB*PINV
         DO 150 RST = 1, MLTPZ
            EXPAB(PQ,RST,I,1) = EXPAP
            EXPAB(PQ,RST,I,2) = EXPBP
            EXPAB(PQ,RST,I,3) = PINV
  150    CONTINUE
  100 CONTINUE
C
      DO I = 1, NPRFCD
      DO PQ = 1, NPQBCS
         CD = IPNTPP(PQ,2) + I - 1
         EXPC = RODPE2(CD,1)
         EXPD = RODPE2(CD,2)
         QINV = D1/(EXPC + EXPD)
         EXPCQ = EXPC*QINV
         EXPDQ = EXPD*QINV
         DO RST = 1, MLTPZ
            EXPCD(PQ,RST,I,1) = EXPCQ
            EXPCD(PQ,RST,I,2) = EXPDQ
            EXPCD(PQ,RST,I,3) = QINV
         END DO
      END DO
      END DO
C
      IJ = 0
      DO J = 1, NPRFCD
      DO I = 1, NPRFAB
         IJ = IJ + 1
         DO PQ = 1, NPQBCS
            PINV = EXPAB(PQ,1,I,3) 
            QINV = EXPCD(PQ,1,J,3) 
            ALFA = D2/(PINV+QINV)
            IF (R12EIN) THEN
               FACPQ(PQ,IJ) = SFAC*SQRT(D2)
            ELSE
               FACPQ(PQ,IJ) = SFAC*SQRT(PINV*QINV*ALFA)
            END IF
            ALPHA(PQ,1,IJ) = -ALFA
            HEXPP(PQ,1,IJ) = D1/PINV
            HEXPQ(PQ,1,IJ) = D1/QINV
         END DO
         DO RST = 2, MLTPZ
         DO PQ = 1, NPQBCS
            ALPHA(PQ,RST,IJ) = ALPHA(PQ,1,IJ)
            HEXPP(PQ,RST,IJ) = HEXPP(PQ,1,IJ)
            HEXPQ(PQ,RST,IJ) = HEXPQ(PQ,1,IJ)
         END DO
         END DO
      END DO
      END DO
C
      IF (IPRINT .GT. 25) THEN
         CALL HEADER('EXPAB in ERIEXP',-1)
         CALL OUTPUT(EXPAB,1,NPPAB,1,3,NPPAB,3,1,LUPRI)
         CALL HEADER('EXPCD in ERIEXP',-1)
         CALL OUTPUT(EXPCD,1,NPPCD,1,3,NPPCD,3,1,LUPRI)
         CALL HEADER('FACPQ in ERIEXP',-1)
         CALL OUTPUT(FACPQ,1,NPQBCS, 1,MLTPZ,NPQBCX,MLTPZ,1,LUPRI)
         CALL HEADER('ALPHA in ERIEXP',-1)
         CALL OUTPUT(ALPHA,1,NPPX, 1,1,NPPX, 1,1,LUPRI)
         CALL HEADER('HEXPP in ERIEXP',-1)
         CALL OUTPUT(HEXPP,1,NPPX, 1,1,NPPX, 1,1,LUPRI)
         CALL HEADER('HEXPQ in ERIEXP',-1)
         CALL OUTPUT(HEXPQ,1,NPPX, 1,1,NPPX, 1,1,LUPRI)
      END IF
C
      RETURN
      END
C  /* Deck ericor */
      SUBROUTINE ERICOR(COORTR,IPNTCR,CORBCH,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxorb.h>
      INTEGER PQ
      DIMENSION COORTR(NPQBCS,3,4),
     &          IPNTCR(MAXBCH,4), CORBCH(NPQBCS,3)
#include <cbieri.h>
#include <odclss.h>
#include <ericom.h>
#include <hertop.h>
#include <aobtch.h>
#include <symmet.h>
#include <ibtfun.h>
C
      DO IATOM = 1,4
#if defined (VAR_CRY)
         CALL GATHER(NPQBCS,CORBCH(1,1),CORXBT,IPNTCR(1,IATOM))
         CALL GATHER(NPQBCS,CORBCH(1,2),CORYBT,IPNTCR(1,IATOM))
         CALL GATHER(NPQBCS,CORBCH(1,3),CORZBT,IPNTCR(1,IATOM))
#else
         DO PQ = 1, NPQBCS
            CORBCH(PQ,1) = CORXBT(IPNTCR(PQ,IATOM))
            CORBCH(PQ,2) = CORYBT(IPNTCR(PQ,IATOM))
            CORBCH(PQ,3) = CORZBT(IPNTCR(PQ,IATOM))
         END DO
#endif
         DO PQ = 1, NPQBCS
            COORTR(PQ,1,IATOM) = CORBCH(PQ,1)
            COORTR(PQ,2,IATOM) = CORBCH(PQ,2)
            COORTR(PQ,3,IATOM) = CORBCH(PQ,3)
         END DO
      END DO
C
      IF (IPRINT .GT. 25) THEN
         CALL TITLER('Output from ERICOR','*',103)
         CALL HEADER('Coordinates A in ERICOR',-1)
         CALL OUTPUT(COORTR(1,1,1),1,NPQBCS,1,3,NPQBCS,3,1,LUPRI)
         CALL HEADER('Coordinates B in ERICOR',-1)
         CALL OUTPUT(COORTR(1,1,2),1,NPQBCS,1,3,NPQBCS,3,1,LUPRI)
         CALL HEADER('Coordinates C in ERICOR',-1)
         CALL OUTPUT(COORTR(1,1,3),1,NPQBCS,1,3,NPQBCS,3,1,LUPRI)
         CALL HEADER('Coordinates D in ERICOR',-1)
         CALL OUTPUT(COORTR(1,1,4),1,NPQBCS,1,3,NPQBCS,3,1,LUPRI)
      END IF
C
      RETURN
      END
C  /* Deck erisso */
      SUBROUTINE ERISSO(IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxorb.h>
      PARAMETER (D1 = 1.0D0)
      INTEGER PQ, R, S, T, RST, X
#include <cbieri.h>
#include <odclss.h>
#include <ericom.h>
#include <hertop.h>
#include <aobtch.h>
#include <symmet.h>
#include <erisop.h>
#include <ibtfun.h>
      XPAR(I) = PT(IBTAND(ISYMAX(1,1),I))
      YPAR(I) = PT(IBTAND(ISYMAX(2,1),I))
      ZPAR(I) = PT(IBTAND(ISYMAX(3,1),I))
C
      RST = 0
      IT = 0
      DO T = 0, MAXOPR
      IF (IBTAND(T,ISTBLT) .EQ. 0) THEN
         IT = IT + 1
         XC0 = XPAR(T)
         YC0 = YPAR(T)
         ZC0 = ZPAR(T)
         IS = 0
         DO S = 0, MAXOPR
         IF (IBTAND(S,ISTBLS) .EQ. 0) THEN
            IS = IS + 1
            XD0 = XPAR(IBTXOR(S,T))
            YD0 = YPAR(IBTXOR(S,T))
            ZD0 = ZPAR(IBTXOR(S,T))
            IR = 0
            DO R = 0, MAXOPR
            IF (IBTAND(R,ISTBLR) .EQ. 0) THEN
               IR = IR + 1
               RST = RST + 1
               NSOP(RST,1)  = IR
               NSOP(RST,2)  = IS
               NSOP(RST,3)  = IT
               NSOQ(RST,1)  = R
               NSOQ(RST,2)  = T
               NSOQ(RST,3)  = IBTXOR(S,T)
               SGNX(RST,1)  = D1
               SGNX(RST,2)  = D1
               SGNX(RST,3)  = D1
               SGNX(RST,4)  = XPAR(R)
               SGNX(RST,5)  = YPAR(R)
               SGNX(RST,6)  = ZPAR(R)
               SGNX(RST,7)  = XC0
               SGNX(RST,8)  = YC0
               SGNX(RST,9)  = ZC0
               SGNX(RST,10) = XD0
               SGNX(RST,11) = YD0
               SGNX(RST,12) = ZD0
            END IF
            END DO
         END IF
         END DO
      END IF
      END DO
C
      IF (IPRINT .GT. 25) THEN
         CALL TITLER('Output from ERISSO','*',103)
      END IF
C
      RETURN
      END
C  /* Deck eripfs */
      SUBROUTINE ERIPFS(FACINT,FACPQ,RODAB,RODCD,RODPF1,RODPF2,IPNTPP,
     &                  ROD1,ROD2,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxorb.h>
#include <symmet.h>
      INTEGER AB, CD, PQ, RST
      DIMENSION FACINT(NPQBCS,MLTPZ,NPRFPQ),
     &          FACPQ(NPQBCS,NPRFPQ), 
     &          RODAB(NPQBCS,MLTPZ,NPRFAB),
     &          RODCD(NPQBCS,MLTPZ,NPRFCD),
     &          RODPF1(MAXOPR+1,*), RODPF2(MAXOPR+1,*),
     &          ROD1(NPQBCS,NPRFAB,3), ROD2(NPQBCS,NPRFCD,3),
     &          IPNTPP(MAXBCH,2)
#include <cbieri.h>
#include <odclss.h>
#include <ericom.h>
#include <hertop.h>
#include <aobtch.h>
#include <erisop.h>
C
      IF (IPRINT .GT. 10) THEN
         CALL TITLER('Output from ERIPFS','*',103)
      END IF
C
C     Collect data 
C     ============
C
      DO RST = 1, MLTPZ
         IR = NSOP(IOFRST+RST,1)
         DO I  = 1, NPRFAB
         DO PQ = 1, NPQBCS
            RODAB(PQ,RST,I) = RODPF1(IR,IPNTPP(PQ,1) + I - 1)
         END DO
         END DO
      END DO
      DO RST = 1, MLTPZ
         IR = NSOP(IOFRST+RST,2)
         DO I  = 1, NPRFCD
         DO PQ = 1, NPQBCS
            RODCD(PQ,RST,I) = RODPF2(IR,IPNTPP(PQ,2) + I - 1)
         END DO
         END DO
      END DO
C
C     Preexponential and other overall factors
C     ========================================
C
      IJ = 0
      DO J = 1, NPRFCD
      DO I = 1, NPRFAB
         IJ = IJ + 1
         DO RST = 1, MLTPZ
         DO PQ  = 1, NPQBCS
            FACINT(PQ,RST,IJ) = 
     &         FACPQ(PQ,IJ)*RODAB(PQ,RST,I)*RODCD(PQ,RST,J)
         END DO
         END DO
      END DO
      END DO
C
C     Print
C     =====
C
      IF (IPRINT .GT. 25) THEN
         CALL HEADER('FACINT in ERIPFS',-1)
         CALL OUTPUT(FACINT,1,NPQBCS*MLTPZ,1,NPRFPQ,
     &               NPQBCS*MLTPZ,NPRFPQ,1,LUPRI)
      END IF
C
      RETURN
      END
C  /* Deck ericrs */
      SUBROUTINE ERICRS(COORAO,COORTR,IPRINT)
#include <implicit.h>
#include <priunit.h>
      INTEGER PQ, RST, X
      DIMENSION COORTR(NPQBCS,12),
     &          COORAO(NPQBCS,MLTPZ,12)
#include <cbieri.h>
#include <ericom.h>
#include <erisop.h>
C
      DO X = 1, 12
         DO RST = 1, MLTPZ
            FAC = SGNX(IOFRST+RST,X)
            DO PQ = 1, NPQBCS
               COORAO(PQ,RST,X) = FAC*COORTR(PQ,X)
            END DO
         END DO
      END DO
C
      IF (IPRINT .GT. 25) THEN
         CALL TITLER('Output from ERICRS','*',103)
         CALL HEADER('Coordinates A in ERICRS',-1)
         CALL OUTPUT(COORAO(1,1,1),1,NPQBCX,1,3,NPQBCX,3,1,LUPRI)
         CALL HEADER('Coordinates B in ERICRS',-1)
         CALL OUTPUT(COORAO(1,1,4),1,NPQBCX,1,3,NPQBCX,3,1,LUPRI)
         CALL HEADER('Coordinates C in ERICRS',-1)
         CALL OUTPUT(COORAO(1,1,7),1,NPQBCX,1,3,NPQBCX,3,1,LUPRI)
         CALL HEADER('Coordinates D in ERICRS',-1)
         CALL OUTPUT(COORAO(1,1,10),1,NPQBCX,1,3,NPQBCX,3,1,LUPRI)
      END IF
C
      RETURN
      END
C  /* Deck erivcs */
      SUBROUTINE ERIVCS(EXPAB,EXPCD,COORAO,COORAB,COORCD,COORPQ,IPRINT)
#include <implicit.h>
#include <priunit.h>
      INTEGER X
      DIMENSION COORAO(NPQBCX,3,4),
     &          EXPAB (NPQBCX,NPRFAB,3), 
     &          EXPCD (NPQBCX,NPRFCD,3), 
     &          COORAB(NPQBCX,NPRFAB,3,2), 
     &          COORCD(NPQBCX,NPRFCD,3,2), 
     &          COORPQ(NPQBCX,NPRFAB,NPRFCD,3)
#include <ericom.h>
#include <eriao.h>
#include <erithr.h>
#include <hertop.h>
#include <ibtfun.h>
C
C     AB coordinates
C     ==============
C
      DO X = 1, 3
         DO K = 1, NPQBCX 
            AX = COORAO(K,X,1)
            BX = COORAO(K,X,2)
            ABX = AX - BX
            DO I = 1, NPRFAB
               COORAB(K,I,X,1) = ABX
               COORAB(K,I,X,2) = EXPAB(K,I,1)*AX + EXPAB(K,I,2)*BX
            END DO
         END DO
      END DO
C
C     CD coordinates
C     ==============
C
      DO X = 1, 3
         DO K = 1, NPQBCX 
            CX = COORAO(K,X,3)
            DX = COORAO(K,X,4)
            CDX = CX - DX
            DO I = 1, NPRFCD
               COORCD(K,I,X,1) = CDX 
               COORCD(K,I,X,2) = EXPCD(K,I,1)*CX + EXPCD(K,I,2)*DX
            END DO
         END DO
      END DO
C
C     PQ coordinates
C     ==============
C
      DO X = 1, 3
         DO J = 1, NPRFCD
         DO I = 1, NPRFAB
            DO K = 1, NPQBCX 
               COORPQ(K,I,J,X) = COORAB(K,I,X,2) - COORCD(K,J,X,2)
            END DO
         END DO
         END DO
      END DO
C
C     Local symmetries
C     ================
C
      CALL LOCSYM(COORAB,NPQBCX,NPRFAB,IAB0,IABXYZ)
      CALL LOCSYM(COORCD,NPQBCX,NPRFCD,ICD0,ICDXYZ)
      CALL LOCSYM(COORPQ,NPPX,1,IPQ0,IPQXYZ)
      IF (BDER) THEN
         CALL LOCSYB(COORAO(1,1,1),IAB0,IABXYZ)
         CALL LOCSYB(COORAO(1,1,3),ICD0,ICDXYZ)
      END IF
      IHHXYZ = IPQXYZ
      IHCXYZ = IBTAND(IPQXYZ,IABXYZ)
      ICHXYZ = IBTAND(IPQXYZ,ICDXYZ)
      ICCXYZ = IBTAND(IPQXYZ,IBTAND(IABXYZ,ICDXYZ))
C
      IF (IPATH .EQ. 1) THEN
         I12XYZ  = IABXYZ
         I120(1) = IAB0(1)
         I120(2) = IAB0(2)
         I120(3) = IAB0(3)
C
         I34XYZ  = ICDXYZ
         I340(1) = ICD0(1)
         I340(2) = ICD0(2)
         I340(3) = ICD0(3)
C
         IHHBIT = IHHXYZ
         IHCBIT = IHCXYZ
         ICHBIT = ICHXYZ
         ICCBIT = ICCXYZ
      ELSE
         I12XYZ  = ICDXYZ
         I120(1) = ICD0(1)
         I120(2) = ICD0(2)
         I120(3) = ICD0(3)
C
         I34XYZ  = IABXYZ
         I340(1) = IAB0(1)
         I340(2) = IAB0(2)
         I340(3) = IAB0(3)
C
         IHHBIT = IHHXYZ
         IHCBIT = ICHXYZ
         ICHBIT = IHCXYZ
         ICCBIT = ICCXYZ
      END IF
C
C     Print section
C     =============
C
      IF (IPRINT .GT. 25) THEN
         CALL HEADER('COORPQ in ERIVCS',1)
         CALL OUTPUT(COORPQ,1,NPPX,1,3,NPPX,3,1,LUPRI)
         IF (JMAXAB .GT. 0) THEN
            CALL HEADER('COORAB (A-B) in ERIVCS',1)
            CALL OUTPUT(COORAB(1,1,1,1),1,NPPAB,1,3,NPPAB,3,1,LUPRI)
         END IF
         IF (JMAXCD .GT. 0) THEN
            CALL HEADER('COORCD (C-D) in ERIVCS',1)
            CALL OUTPUT(COORCD(1,1,1,1),1,NPPCD,1,3,NPPCD,3,1,LUPRI)
         END IF
         WRITE (LUPRI,'(/2X,A,E20.10,/)') 'THRSH', THRSH
         WRITE (LUPRI,'(2X,A,11X,I5,5X,3I5)')
     &      'IPQXYZ, IPQ0(3):', IPQXYZ, IPQ0(1), IPQ0(2), IPQ0(3)
         WRITE (LUPRI,'(2X,A,11X,I5,5X,3I5)')
     &      'IABXYZ, IAB0(3):', IABXYZ, IAB0(1), IAB0(2), IAB0(3)
         WRITE (LUPRI,'(2X,A,11X,I5,5X,3I5)')
     &      'ICDXYZ, ICD0(3):', ICDXYZ, ICD0(1), ICD0(2), ICD0(3)
         WRITE (LUPRI,'(2X,A,4I5)')
     &      'IHHXYZ,IHCXYZ,ICHXYZ,ICCXYZ',
     &       IHHXYZ,IHCXYZ,ICHXYZ,ICCXYZ
      END IF
      RETURN
      END
C  /* Deck locsym */
      SUBROUTINE LOCSYM(VECTOR,NINNER,NDIM,I0,IXYZ)
#include <implicit.h>
      PARAMETER (D0 = 0.0D0)
      INTEGER X
      DIMENSION VECTOR(NINNER,NDIM,3), I0(3)
#include <erithr.h>
#include <cbieri.h>
C
      IF (NOLOCS) THEN
         I0(1) = 0
         I0(2) = 0
         I0(3) = 0
      ELSE
#if defined (SYS_CRAY)
         I0(1) = 1
         I0(2) = 1
         I0(3) = 1
         DX = D0
         DY = D0
         DZ = D0
         DO I = 1, NINNER
            DX = MAX(DX,ABS(VECTOR(I,1,1)))
            DY = MAX(DY,ABS(VECTOR(I,1,2)))
            DZ = MAX(DZ,ABS(VECTOR(I,1,3)))
         END DO
         IF (DX .GT. THRSH) I0(1) = 0
         IF (DY .GT. THRSH) I0(2) = 0
         IF (DZ .GT. THRSH) I0(3) = 0
#else
         DO X = 1, 3
            I0(X) = 1
            DO I = 1, NINNER
               IF (ABS(VECTOR(I,1,X)) .GT. THRSH) THEN
                  I0(X) = 0
                  GO TO 100
               END IF
            END DO
  100       CONTINUE
         END DO
#endif
      END IF
      IXYZ = I0(1) + 2*I0(2) + 4*I0(3)
      RETURN
      END
C  /* Deck locsyb */
      SUBROUTINE LOCSYB(COORAO,I0,IXYZ)
#include <implicit.h>
      INTEGER X
      DIMENSION COORAO(NPQBCX,3), I0(3)
#include <orgcom.h>
#include <erithr.h>
#include <cbieri.h>
#include <ericom.h>
C
      IF (.NOT.NOLOCS) THEN
         DO X = 1, 3
            DO I = 1, NPQBCX
               IF (ABS(COORAO(I,X)-ORIGIN(X)) .GT. THRSH) THEN
                  I0(X) = 0
                  GO TO 100
               END IF
            END DO
  100       CONTINUE
         END DO
      END IF
      IXYZ = I0(1) + 2*I0(2) + 4*I0(3)
      RETURN
      END
C  /* Deck ericnt */
      SUBROUTINE ERICNT(NCENTR,IPNTCR,IPRINT)
C
C     tuh
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
      INTEGER RST, PQ, A
      DIMENSION NCENTR(MLTPZ,NPQBCS,4), IPNTCR(MAXBCH,4)
#include <nuclei.h>
#include <aobtch.h>
#include <ericom.h>
#include <erisop.h>
#include <ibtfun.h>
      IBTAXO(I,J) = IBTAND(I,IBTXOR(I,J))
C
      IF (IPRINT .GT. 10) THEN
         CALL TITLER('Output from ERICNT','*',103)
      END IF
C
      DO PQ = 1, NPQBCS
         IA = NCNTBT(IPNTCR(PQ,1))
         IB = NCNTBT(IPNTCR(PQ,2))
         IC = NCNTBT(IPNTCR(PQ,3))
         ID = NCNTBT(IPNTCR(PQ,4))
         ISTB = ISTBNU(IB)
         ISTC = ISTBNU(IC)
         ISTD = ISTBNU(ID)
         DO RST = 1, MLTPZ
            NCENTR(RST,PQ,1) = NUCNUM(IA,1)
            NCENTR(RST,PQ,2) = NUCNUM(IB,IBTAXO(NSOQ(RST,1),ISTB)+1)
            NCENTR(RST,PQ,3) = NUCNUM(IC,IBTAXO(NSOQ(RST,2),ISTC)+1)
            NCENTR(RST,PQ,4) = NUCNUM(ID,IBTAXO(NSOQ(RST,3),ISTD)+1)
         END DO
      END DO
      IF (IPRINT .GT. 10) THEN
         DO A = 1, 4
            WRITE (LUPRI,'(1X,A,I5)') ' NCENTR for orbital ',A
            DO PQ = 1, NPQBCS
               WRITE (LUPRI,'(2X,39I2)') (NCENTR(RST,PQ,A),RST=1,MLTPZ)
            END DO
         END DO
      END IF
      RETURN
      END
C  /* Deck ericls */
      SUBROUTINE ERICLS(FACINT,COORPQ,EXPAB,EXPCD,IPRINT)
#include <implicit.h>
#include <priunit.h>
C
C     tuh
C
#include <maxaqn.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxorb.h>
      PARAMETER (D0 = 0.0D0)
      LOGICAL CLASIC 
      DIMENSION FACINT(NPQBCX,NPRFAB,NPRFCD),
     &          COORPQ(NPQBCX,NPRFAB,NPRFCD,3), 
     &          EXPAB (NPQBCX,NPRFAB,3),
     &          EXPCD (NPQBCX,NPRFCD,3)
#include <cbieri.h>
#include <odclss.h>
#include <ericom.h>
#include <hertop.h>
#include <aobtch.h>
#include <erisop.h>
#include <erithr.h>
#include <symmet.h>
#include <clsfmm.h>
C
      IF (IPRINT .GT. 10) THEN
         CALL TITLER('Output from ERICLS','*',103)
      END IF
C
      FAC = ERFCIV(THRCLS)
      NNCLS = 0
      DO K  = 1, NPQBCX
         CLASIC = .TRUE.
         DO J = 1, NPRFCD
         DO I = 1, NPRFAB
            RPQ = SQRT(COORPQ(K,I,J,1)**2 + COORPQ(K,I,J,2)**2
     &                                    + COORPQ(K,I,J,3)**2)
            EPQ = FAC*(SQRT(EXPAB(K,I,3)) + SQRT(EXPCD(K,J,3)))
            IF (RPQ .LE. EPQ) CLASIC = .FALSE. 
         END DO
         END DO
         IF (CLASIC) THEN
            NNCLS = NNCLS + 1
            DO J = 1, NPRFCD
            DO I = 1, NPRFAB
               FACINT(K,I,J) = D0
            END DO
            END DO
         END IF
      END DO
C
      IF (IPRINT .GT. 25) THEN
         CALL HEADER('FACINT in ERICLS',-1)
         WRITE (LUPRI,'(2X,A,2I5)') ' NPQBCX, NNCLS ', NPQBCX, NNCLS
         CALL OUTPUT(FACINT,1,NPQBCX,1,NPRFPQ,NPQBCX,NPRFPQ,1,LUPRI)
      END IF
C
      RETURN
      END
