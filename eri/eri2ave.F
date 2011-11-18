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
C  /* Deck eriave */
      SUBROUTINE ERIAVE(SO,PMAT,DMAT,D2MAT,ID2MAT,IPNTCR,IODDCC,IPNTUV,
     &                  INDXBT,WORK,LWORK,IPRINT)
C
C     tuh march 99 - april 01
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
C
      DIMENSION SO(*), PMAT(*), DMAT(*), D2MAT(*), ID2MAT(*),
     &          IPNTCR(*), IPNTUV(*), IODDCC(*),
     &          INDXBT(*), WORK(*)
#include <ericom.h>
C
C     Allocations
C
      KPOINT = 1 
      KLAST  = KPOINT + (5*NCCS - 1)/IRAT + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('ERIAVE',' ',KLAST,LWORK)
      CALL ERIAV1(SO,DMAT,D2MAT,ID2MAT,IPNTCR,IODDCC,IPNTUV,
     &            INDXBT,WORK(KPOINT),IPRINT)
      RETURN
      END
C  /* Deck eriav1 */
      SUBROUTINE ERIAV1(SO,DMAT,D2MAT,ID2MAT,IPNTCR,IODDCC,IPNTUV,
     &                  INDXBT,IPOINT,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
      PARAMETER (D1=1.0D0, DP5=0.5D0)
      INTEGER A, B, C, D, R, S, T
      LOGICAL DOREP(0:7,4), CTRIAB, CTRICD, DAB, DCD
      DIMENSION IPOINT(NCCS,5), 
     &          IPNTCR(MAXBCH,4), INDXBT(MXSHEL*MXCONT,0:7),
     &          IODDCC(NRTOP), IPNTUV(KC2MAX,0:NRDER,2),
     &          IPNRST(0:7,3), IODXYZ(3),
     &          IADCMP(MXAQN,MXAQN,2),
     &          SO(NCCS,MLTPR,MLTPS,MLTPT,KHKTAB,KHKTCD,12),
     &          DMAT(NBASE,NBASE),
     &          D2MAT(NBASE,NBASE,NPDIMB,NPDIMA),
     &          ID2MAT(MXCORB,2)
#include <cbieri.h>
#include <ericom.h>
#include <erithr.h>
#include <eribuf.h>
#include <aobtch.h>
#include <hertop.h>
#include <symmet.h>
#include <nuclei.h>
#include <energy.h>
C
#include <ibtfun.h>
      IBTEST(I,J,K,L) = IBTAND(I,IBTXOR(J,ISYMAO(K,L)))
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine ERIAV1',-1)
C
      IF (GRDZER) CALL DZERO(GRADEE,MXCOOR)
C
      CALL PRPREP(DOREP(0,1),NHKTA,KHKTA,ISTBLA)
      CALL PRPREP(DOREP(0,2),NHKTB,KHKTB,ISTBLB)
      CALL PRPREP(DOREP(0,3),NHKTC,KHKTC,ISTBLC)
      CALL PRPREP(DOREP(0,4),NHKTD,KHKTD,ISTBLD)
C
      IF (IPRINT .GT. 10) THEN
         WRITE (LUPRI,'(/,2X,A,8L2)')'DOREP A  ',(DOREP(I,1),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,8L2)')  'DOREP B  ',(DOREP(I,2),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,8L2)')  'DOREP C  ',(DOREP(I,3),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,8L2)')  'DOREP D  ',(DOREP(I,4),I=0,MAXREP)
      END IF
C
      CALL CMPADR(IADCMP(1,1,1),KHKTA,KHKTB,TKMPAB)
      CALL CMPADR(IADCMP(1,1,2),KHKTC,KHKTD,TKMPCD)
C
      CALL GETRST(IPNRST(0,1),ISTBLR)
      CALL GETRST(IPNRST(0,2),ISTBLS)
      CALL GETRST(IPNRST(0,3),ISTBLT)
C
      DO A = 0, MAXREP
      IF (DOREP(A,1)) THEN
      DO B = 0, MAXREP
      IF (DOREP(B,2)) THEN
      DO C = 0, MAXREP
      IF (DOREP(C,3) .AND. DOREP(IBTXOR(IBTXOR(A,B),C),4)) THEN
         D = IBTXOR(IBTXOR(A,B),C)
         IF (DIAGAB .AND. B.GT.A) GO TO 100
         IF (DIAGCD .AND. D.GT.C) GO TO 100
C
         CTRIAB = DIAGAB .AND. A.EQ.B 
         CTRICD = DIAGCD .AND. C.EQ.D 
C
         R = IPNRST(B,1)
         S = IPNRST(D,2)
         T = IPNRST(IBTXOR(C,D),3)
C
         CALL ERIPNT(IPOINT,A,B,C,D,IPNTCR,INDXBT,0)
C
         IA   = -1
         MAXB = KHKTB
         MAXD = KHKTD
         DO ICMPA = 1, KHKTA
         IVARA = IBTEST(ISTBLA,A,NHKTA,ICMPA)
         IF (IVARA.EQ.0) THEN
            IA = IA + 1
            IB = -1
            IF (CTRIAB) MAXB = ICMPA
            DO ICMPB = 1, MAXB
            IVARB = IBTEST(ISTBLB,B,NHKTB,ICMPB)
            IF (IVARB.EQ.0) THEN
               IB = IB + 1
               IC = -1
               ICMPAB = IADCMP(ICMPA,ICMPB,1)
               IODXYZ(1) = IODDCC(IPNTUV(ICMPAB,1,1))
               IODXYZ(2) = IODDCC(IPNTUV(ICMPAB,2,1))
               IODXYZ(3) = IODDCC(IPNTUV(ICMPAB,3,1))
               DO ICMPC = 1, KHKTC
               IVARC = IBTEST(ISTBLC,C,NHKTC,ICMPC)
               IF (IVARC.EQ.0) THEN
                  IC = IC + 1
                  ID = -1
                  IF (CTRICD) MAXD = ICMPC
                  DO ICMPD = 1, MAXD
                  IVARD = IBTEST(ISTBLD,D,NHKTD,ICMPD)
                  IF (IVARD.EQ.0) THEN
                     ID = ID + 1
                     ICMPCD = IADCMP(ICMPC,ICMPD,2)
                     IODDCD = IODDCC(IPNTUV(ICMPCD,0,2))
                     IF(IODDCD.EQ.IODXYZ(1).OR.IODDCD.EQ.IODXYZ(2) 
     &                                     .OR.IODDCD.EQ.IODXYZ(3))THEN
                        DAB = GCONAB .AND. CTRIAB.AND.ICMPA.EQ.ICMPB
                        DCD = GCONCD .AND. CTRICD.AND.ICMPC.EQ.ICMPD
                        FAC = D1
                        IF (DIAGPQ) FAC = DP5*FAC
                        CALL ERIFRC(IA,IB,IC,ID,
     &                              SO(1,R,S,T,ICMPAB,ICMPCD,1),
     &                              IPOINT,DMAT,D2MAT,ID2MAT,FAC,
     &                              IODDCD,IODXYZ,DAB,DCD)
                     END IF
                  END IF
                  END DO
               END IF
               END DO
            END IF
            END DO
         END IF
         END DO
  100    CONTINUE
      END IF
      END DO
      END IF
      END DO
      END IF
      END DO
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('GRADEE in ERIAV1',-1) 
         CALL OUTPUT(GRADEE,1,3,1,NUCDEP,3,NUCDEP,1,LUPRI)
      END IF
C 
      RETURN
      END
C  /* Deck erifrc */
      SUBROUTINE ERIFRC(IA,IB,IC,ID,SO,IPOINT,DMAT,D2MAT,ID2MAT,
     &                  DFAC,IODDCD,IODXYZ,DAB,DCD)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
      PARAMETER (D4=4.0D0, D2=2.0D0, D1=1.0D0, DP5=0.5D0, DP25=0.25D0)
      LOGICAL DAB, DCD, NOTRAS
      DIMENSION SO(NAOINT,3,4), IPOINT(NCCS,5), DMAT(NBASE,NBASE), 
     &          D2MAT(NBASE,NBASE,NPDIMB,NPDIMA), ID2MAT(MXCORB,2),
     &          IODXYZ(3)
#include <cbieri.h>
#include <ericom.h>
#include <eribuf.h>
#include <aobtch.h>
#include <hertop.h>
#include <nuclei.h>
#include <symmet.h>
#include <energy.h>
C    
      NOTRAS = .TRUE.
      DO I = 1, NCCS
         KA = IPOINT(I,1) + IA
         KB = IPOINT(I,2) + IB
         KC = IPOINT(I,3) + IC
         KD = IPOINT(I,4) + ID
         IF (CCRUN) THEN
            PMAT = DFAC*D2MAT(KD,KC,ID2MAT(KB,2),ID2MAT(KA,1))
            IF (DCD .OR. KC.EQ.KD) PMAT = DP5*PMAT
         ELSE
            PMAT = DFAC*D4*(DMAT(KA,KB)*DMAT(KC,KD) 
     &               - DP25*DMAT(KA,KC)*DMAT(KB,KD)
     &               - DP25*DMAT(KA,KD)*DMAT(KB,KC))
            IF (DAB .OR. KA .EQ. KB) PMAT = DP5*PMAT
            IF (DCD .OR. KC .EQ. KD) PMAT = DP5*PMAT
         END IF
C
         JA = 3*ICNTAO(KA) - 3
         JB = 3*ICNTAO(KB) - 3
         JC = 3*ICNTAO(KC) - 3
         JD = 3*ICNTAO(KD) - 3
C
         IF (NOTRAS) THEN
            J = IPOINT(I,5) 
            DO ICOOR = 1, 3
            IF (IODDCD .EQ. IODXYZ(ICOOR)) THEN
               ISCORA = IPTCNT(JA+ICOOR,0,1)
               ISCORB = IPTCNT(JB+ICOOR,0,1)
               ISCORC = IPTCNT(JC+ICOOR,0,1)
               ISCORD = IPTCNT(JD+ICOOR,0,1)
               IF (ISCORA.GT.0) GRADEE(ISCORA) 
     &                        = GRADEE(ISCORA) + PMAT*SO(J,ICOOR,1)
               IF (ISCORB.GT.0) GRADEE(ISCORB) 
     &                        = GRADEE(ISCORB) + PMAT*SO(J,ICOOR,2)
               IF (ISCORC.GT.0) GRADEE(ISCORC) 
     &                        = GRADEE(ISCORC) + PMAT*SO(I,ICOOR,3)
               IF (ISCORD.GT.0) GRADEE(ISCORD) 
     &                        = GRADEE(ISCORD) + PMAT*SO(I,ICOOR,4)
            END IF
            END DO
         ELSE
            DO ICOOR = 1, 3
            IF (IODDCD .EQ. IODXYZ(ICOOR)) THEN
               FA = PMAT*SO(J,ICOOR,1)
               FB = PMAT*SO(I,ICOOR,2)
               FC = PMAT*SO(I,ICOOR,3)
               FD = - FA - FB - FC
               ISCORA = IPTCNT(JA+ICOOR,0,1)
               ISCORB = IPTCNT(JB+ICOOR,0,1)
               ISCORC = IPTCNT(JC+ICOOR,0,1)
               ISCORD = IPTCNT(JD+ICOOR,0,1)
               IF (ISCORA.GT.0) GRADEE(ISCORA) = GRADEE(ISCORA) + FA 
               IF (ISCORB.GT.0) GRADEE(ISCORB) = GRADEE(ISCORB) + FB 
               IF (ISCORC.GT.0) GRADEE(ISCORC) = GRADEE(ISCORC) + FC 
               IF (ISCORD.GT.0) GRADEE(ISCORD) = GRADEE(ISCORD) + FD 
            END IF
            END DO
         END IF
      END DO
      RETURN
      END
