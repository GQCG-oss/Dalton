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
C  /* Deck ccsd_triple */
      SUBROUTINE CCSD_TRIPLE(OMEGA1,OMEGA2,T1AM,T2AM,FOCK,
     *                       XLAMDP,XLAMDH,WORK,LWORK)
C
C     Written by Henrik Koch 19-Sep-1994
C
C     Version 1.0
C
C     Purpose:
C
C     Calculate the iterative triples corrections to the CCSD
C     equations.
C
C     NB! The T2 amplitudes are assumed to be a full square.
C
C
C     NB! It is assumed that the vectors are allocated in the following
C     order:
C           T1AM(*), OMEGA1(*), OMEGA2(*), T2AM(*), SCR(*), WRK(*).
C
C
      IMPLICIT NONE  
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <ccorb.h>
#include <inftap.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccfield.h>
C
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

#if defined (SYS_CRAY)
      REAL ZERO, HALF, ONE, TWO
#else
      DOUBLE PRECISION ZERO, HALF, ONE, TWO
#endif
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      INTEGER LWORK
#if defined (SYS_CRAY)
      REAL OMEGA1(*),OMEGA2(*),T1AM(*),T2AM(*),FOCK(*)
      REAL XLAMDP(*),XLAMDH(*),WORK(LWORK)
      REAL DDOT, FF
#else
      DOUBLE PRECISION OMEGA1(*),OMEGA2(*),T1AM(*),T2AM(*),FOCK(*)
      DOUBLE PRECISION XLAMDP(*),XLAMDH(*),WORK(LWORK)
      DOUBLE PRECISION DDOT, FF
#endif
C
      INTEGER INDEX, KSCR1, KFOCKD, KINT1T, KINT2T,
     &        KINT1S, KINT2S, KT3AM, KOME1, KOME2, KIAJB, KYIAJB,
     &        KEND1, LWRK1, KT3AM2, KONEP, KEND2, 
     &        ISYMD1, ILLL, ISYMD, KXINT, LWRK2,
     &        IJ, NIJ, ISYDIS, IDEL
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      IF (DIRECT) CALL QUIT('Direct not implemented in CCSD_TRIPLE')

      IF (NSYM.NE.1) CALL QUIT('No symmetry in CCSD_TRIPLE')
C
C------------------------
C     Dynamic Allocation.
C------------------------
C
      KSCR1  = 1
      KFOCKD = KSCR1  + NT1AMX
      KINT1T = KFOCKD + NORBT
      KINT2T = KINT1T + NT1AMX*NVIRT*NVIRT
      KINT1S = KINT2T + NRHFT*NRHFT*NT1AMX
      KINT2S = KINT1S + NT1AMX*NVIRT*NVIRT
      KT3AM  = KINT2S + NRHFT*NRHFT*NT1AMX
      KOME1  = KT3AM  + NT1AMX*NT1AMX*NT1AMX
      KOME2  = KOME1  + NT1AMX
      KIAJB  = KOME2  + NT1AMX*NT1AMX
      KYIAJB = KIAJB  + NT1AMX*NT1AMX
      KEND1  = KYIAJB + NT1AMX*NT1AMX
C
      IF ((NONHF) .AND. (NFIELD .GT. 0)) THEN
         KONEP  = KEND1
         KT3AM2 = KONEP + N2BST(ISYMOP)
         KEND1  = KT3AM2 + NT1AMX*NT1AMX*NT1AMX
      ENDIF
C
      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSD_TRIPLE')
      ENDIF
C
C--------------------------------
C     Initialize integral arrays.
C--------------------------------
C
      CALL DZERO(WORK(KINT1T),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2T),NT1AMX*NRHFT*NRHFT)
      CALL DZERO(WORK(KINT1S),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2S),NT1AMX*NRHFT*NRHFT)
      CALL DZERO(WORK(KOME1),NT1AMX)
      CALL DZERO(WORK(KOME2),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KIAJB),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KYIAJB),NT1AMX*NT1AMX)
C
C--------------------------------------------------------
C     Read canonical orbital energies from interface file
C--------------------------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBT)
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C=======================================================
C     Get the one electron integrals for the fields
C=======================================================
C
      IF ((NONHF) .AND. NFIELD .GT. 0) THEN
         CALL DZERO(WORK(KONEP),N2BST(ISYMOP))
         DO I = 1, NFIELD
            FF = EFIELD(I)
            CALL CC_ONEP(WORK(KONEP),WORK(KEND1),LWRK1,FF,1,LFIELD(I))
         ENDDO
         CALL CC_FCKMO(WORK(KONEP),XLAMDP,XLAMDH,
     *                 WORK(KEND1),LWRK1,1,1,1)
      ENDIF
C
C====================================================
C     Start the loop over distributions of integrals.
C====================================================
C
      DO ISYMD1 = 1,NSYM
C
         DO ILLL = 1, NBAS(ISYMD1)
C
            IDEL  = IBAS(ISYMD1) + ILLL
            ISYMD = ISYMD1
C
            ISYDIS = MULD2H(ISYMD,ISYMOP)
C
C           ----------------------------
C           Work space allocation no. 2.
C           ----------------------------
            KXINT  = KEND1
            KEND2  = KXINT + NDISAO(ISYDIS)
            LWRK2  = LWORK - KEND2
C
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
               CALL QUIT('Insufficient space in CCSD_TRIPLE')
            ENDIF
C
C           ---------------------------
C           Read in batch of integrals.
C           ---------------------------
            CALL CCRDAO(WORK(KXINT),IDEL,1,WORK(KEND2),LWRK2,
     *                  IDUMMY,DIRECT)
C
C           --------------------------
C           Calculate integrals needed
C           --------------------------
            CALL CCSDT_TRAN1(WORK(KINT1T),WORK(KINT2T),XLAMDP,
     *                       XLAMDH,WORK(KXINT),IDEL)
C
            CALL CC3_TRAN2(WORK(KIAJB),WORK(KYIAJB),XLAMDP,
     *                       XLAMDH,WORK(KXINT),IDEL)
C
            CALL CCSDT_TRAN3(WORK(KINT1S),WORK(KINT2S),XLAMDP,
     *                       XLAMDH,WORK(KXINT),IDEL)
C
         END DO
      END DO
C
      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'CCSD_TRIPLE> norm2(xiajb):',
     &     DDOT(NT1AMX*NT1AMX,WORK(KIAJB),1,WORK(KIAJB),1)
        WRITE(LUPRI,*) 'CCSD_TRIPLE> norm2(int1t):',
     &     DDOT(NT1AMX*NVIRT*NVIRT,WORK(KINT1T),1,WORK(KINT1T),1)
        WRITE(LUPRI,*) 'CCSD_TRIPLE> norm2(int2t):',
     &     DDOT(NT1AMX*NRHFT*NRHFT,WORK(KINT2T),1,WORK(KINT2T),1)
      END IF
C
C-------------------------------------
C     Calculate the triple amplitudes.
C-------------------------------------
C
      CALL CCSDT_T03AM(WORK(KT3AM),WORK(KINT1S),WORK(KINT2S),T2AM,
     *               WORK(KSCR1),WORK(KFOCKD),WORK(KONEP),WORK(KT3AM2))

COMMENT COMMENT
COMMENT COMMENT
C        LUSIFC = -1
C        CALL GPOPEN(LUSIFC,'T3AMFD','UNKNOWN',' ','FORMATTED',IDUMMY,
C    &               .FALSE.)
C        REWIND LUSIFC
C        DO I = 1, NT1AMX*NT1AMX*NT1AMX
C          WRITE (LUSIFC,*) -WORK(KT3AM-1+I)
C        END DO
C        CALL GPCLOSE(LUSIFC,'KEEP')
COMMENT COMMENT
COMMENT COMMENT
C-----------------------------------------
C     Calculate the triples corrections.
C-----------------------------------------
C
      CALL CCSDT_OMEGA1(WORK(KOME1),WORK(KIAJB),WORK(KT3AM))
C
      CALL CCSDT_OMEGA2(WORK(KOME2),WORK(KINT1T),WORK(KINT2T),
     *                  WORK(KT3AM),FOCK)
C
      DO 300 I = 1,NT1AMX
         OMEGA1(I) = OMEGA1(I) + WORK(KOME1+I-1)
  300 CONTINUE

C
      DO 310 I = 1,NT1AMX
         DO 320 J = 1,I
            IJ = NT1AMX*(I-1) + J
            NIJ = INDEX(I,J)
            OMEGA2(NIJ) = OMEGA2(NIJ) + WORK(KOME2+IJ-1)
  320    CONTINUE
  310 CONTINUE
C
      RETURN
      END
C  /* Deck ccsdt_tran1 */
      SUBROUTINE CCSDT_TRAN1(XINT1,XINT2,XLAMDP,XLAMDH,AOINT,IDEL)
C
C
C
#include <implicit.h>
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
      DIMENSION XINT1(NT1AMX,NVIRT,NVIRT), XINT2(NT1AMX,NRHFT,NRHFT)
      DIMENSION XLAMDP(NBAST,NORBT), XLAMDH(NBAST,NORBT)
      DIMENSION AOINT(NNBAST,NBAST)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO 100 G = 1,NBAST
         DO 110 IB = 1,NBAST
            DO 120 A = 1,NBAST
               NAB = INDEX(A,IB)
C
               if (aoint(nab,g) .eq. 0.0d0) goto 120
               DO 200 D = 1,NVIRT
                  DO 210 B = 1,NVIRT
                     DO 220 K = 1,NRHFT
                        DO 230 C = 1,NVIRT
C
                           NCK = NVIRT*(K-1) + C
C
                           XINT1(NCK,B,D) = XINT1(NCK,B,D)
     *               + AOINT(NAB,G)*XLAMDH(A,NRHFT+C)*XLAMDP(IB,K)
     *                      *XLAMDP(G,NRHFT+B)*XLAMDH(IDEL,NRHFT+D)
C
  230                   CONTINUE
  220                CONTINUE
  210             CONTINUE
  200          CONTINUE
C
               DO 300 J = 1,NRHFT
                  DO 310 L = 1,NRHFT
                     DO 320 K = 1,NRHFT
                        DO 330 C = 1,NVIRT
C
                           NCK = NVIRT*(K-1) + C
C
                           XINT2(NCK,L,J) = XINT2(NCK,L,J)
     *                  + AOINT(NAB,G)*XLAMDH(A,NRHFT+C)*XLAMDP(IB,K)
     *                                *XLAMDP(G,L)*XLAMDH(IDEL,J)
C
  330                   CONTINUE
  320                CONTINUE
  310             CONTINUE
  300          CONTINUE
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccsdt_tran2 */
      SUBROUTINE CCSDT_TRAN2(XIAJB,YIAJB,CMO,AOINT,IDEL)
C
C
C
#include <implicit.h>
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>

      PARAMETER (TWO = 2.0D0)
      DIMENSION XIAJB(NT1AMX,NT1AMX), AOINT(NNBAST,NBAST)
      DIMENSION YIAJB(NT1AMX,NT1AMX)
      DIMENSION CMO(NORBT,NORBT)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO 100 G = 1,NBAST
         DO 110 B = 1,NBAST
            DO 120 A = 1,NBAST
               NAB = INDEX(A,B)
C
               if (aoint(nab,g) .eq. 0.0d0) goto 120
               DO 200 L = 1,NRHFT
                  DO 210 D = 1,NVIRT
                     NDL = NVIRT*(L-1) + D
                     DO 220 K = 1,NRHFT
                        DO 230 C = 1,NVIRT
C
                           NCK = NVIRT*(K-1) + C
C
                           XIAJB(NCK,NDL) = XIAJB(NCK,NDL)
     *                             + AOINT(NAB,G)*
     *        (TWO*CMO(A,NRHFT+C)*CMO(B,K)*CMO(G,NRHFT+D)*CMO(IDEL,L)
     *         - CMO(A,NRHFT+C)*CMO(B,L)*CMO(G,NRHFT+D)*CMO(IDEL,K))
C
                           YIAJB(NCK,NDL) = YIAJB(NCK,NDL)
     *                             + AOINT(NAB,G)*
     *           CMO(A,NRHFT+C)*CMO(B,K)*CMO(G,NRHFT+D)*CMO(IDEL,L)
C
  230                   CONTINUE
  220                CONTINUE
  210             CONTINUE
  200          CONTINUE
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccsdt_tran2 */
      SUBROUTINE CC3_TRAN2(XIAJB,YIAJB,XLAMDP,XLAMDH,AOINT,IDEL)
C
C
C
#include <implicit.h>
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>

      PARAMETER (TWO = 2.0D0)
      DIMENSION XIAJB(NT1AMX,NT1AMX), AOINT(NNBAST,NBAST)
      DIMENSION YIAJB(NT1AMX,NT1AMX)
      DIMENSION XLAMDP(NBAST,NORBT), XLAMDH(NBAST,NORBT)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO 100 G = 1,NBAST
         DO 110 B = 1,NBAST
            DO 120 A = 1,NBAST
               NAB = INDEX(A,B)
C
               if (aoint(nab,g) .eq. 0.0d0) goto 120
               DO 200 L = 1,NRHFT
                  DO 210 D = 1,NVIRT
                     NDL = NVIRT*(L-1) + D
                     DO 220 K = 1,NRHFT
                        DO 230 C = 1,NVIRT
C
                           NCK = NVIRT*(K-1) + C
C
                           XIAJB(NCK,NDL) = XIAJB(NCK,NDL)
     *                             + AOINT(NAB,G)*
     *       (TWO*XLAMDP(A,K)*XLAMDH(B,NRHFT+C)
     *           *XLAMDP(G,L)*XLAMDH(IDEL,NRHFT+D)
     *         -  XLAMDP(A,K)*XLAMDH(B,NRHFT+D)
     *           *XLAMDP(G,L)*XLAMDH(IDEL,NRHFT+C))
C
                           YIAJB(NCK,NDL) = YIAJB(NCK,NDL)
     *                             + AOINT(NAB,G)*
     *           XLAMDP(A,K)*XLAMDH(B,NRHFT+C)
     *          *XLAMDP(G,L)*XLAMDH(IDEL,NRHFT+D)
C
  230                   CONTINUE
  220                CONTINUE
  210             CONTINUE
  200          CONTINUE
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccsdt_omega1 */
      SUBROUTINE CCSDT_OMEGA1(OMEGA1,XIAJB,T3AM)
C
C
C
#include <implicit.h>
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>

      PARAMETER (TWO = 2.0D0)
      DIMENSION OMEGA1(NT1AMX),XIAJB(NT1AMX,NT1AMX)
      DIMENSION T3AM(NT1AMX,NT1AMX,NT1AMX)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO 100 I = 1,NRHFT
         DO 110 A = 1,NVIRT
            NAI = NVIRT*(I-1) + A
C
            DO 120 NBJ = 1,NT1AMX
C
               DO 130 K = 1,NRHFT
                  NAK = NVIRT*(K-1) + A
                  DO 140 C = 1,NVIRT
                     NCK = NVIRT*(K-1) + C
                     NCI = NVIRT*(I-1) + C
                     OMEGA1(NAI) = OMEGA1(NAI) - XIAJB(NCK,NBJ)*
     *               ( T3AM(NAI,NBJ,NCK) - T3AM(NAK,NBJ,NCI) )
  140             CONTINUE
  130          CONTINUE
C
  120       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccsdt_omega2 */
      SUBROUTINE CCSDT_OMEGA2(OMEGA2,XINT1T,XINT2T,T3AM,FOCK)
C
C
C
#include <implicit.h>
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
      PARAMETER (TWO = 2.0D0)
      DIMENSION XINT1T(NT1AMX,NVIRT,NVIRT)
      DIMENSION XINT2T(NT1AMX,NRHFT,NRHFT)
      DIMENSION OMEGA2(NT1AMX,NT1AMX)
      DIMENSION T3AM(NT1AMX,NT1AMX,NT1AMX),FOCK(NORBT,NORBT)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO 100 I = 1,NRHFT
         DO 110 A = 1,NVIRT
            NAI = NVIRT*(I-1) + A
C
            DO 120 J = 1,NRHFT
               DO 130 B = 1,NVIRT
                  NBJ = NVIRT*(J-1) + B
C
                  DO 140 K = 1,NRHFT
                     NBK = NVIRT*(K-1) + B
                     NAK = NVIRT*(K-1) + A
                     DO 150 C = 1,NVIRT
C
                        NCK = NVIRT*(K-1) + C
                        NCJ = NVIRT*(J-1) + C
                        NCI = NVIRT*(I-1) + C
C
                        OMEGA2(NBJ,NAI) = OMEGA2(NBJ,NAI) -
     *      (T3AM(NAI,NBJ,NCK) - T3AM(NAI,NBK,NCJ))*FOCK(K,NRHFT+C)
C
                        DO 160 D = 1,NVIRT
C
                           NDJ = NVIRT*(J-1) + D
                           NDK = NVIRT*(K-1) + D
                           NDI = NVIRT*(I-1) + D
C
                           OMEGA2(NBJ,NAI) = OMEGA2(NBJ,NAI) +
     *  (T3AM(NBK,NCI,NDJ) - TWO*T3AM(NBJ,NCI,NDK) + T3AM(NBJ,NCK,NDI))
     *   *XINT1T(NDK,A,C)
C
  160                   CONTINUE
C
                        DO 170 L = 1,NRHFT
C
                           NBL = NVIRT*(L-1) + B
                           NCL = NVIRT*(L-1) + C
                           NAL = NVIRT*(L-1) + A
C
                           OMEGA2(NBJ,NAI) = OMEGA2(NBJ,NAI) -
     *  (T3AM(NBL,NAK,NCJ) - TWO*T3AM(NBJ,NAK,NCL) + T3AM(NBJ,NAL,NCK))
     *   *XINT2T(NCL,K,I)
C
  170                   CONTINUE
C
  150                CONTINUE
  140             CONTINUE
C
  130          CONTINUE
  120       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      DO 200 NAI = 1,NT1AMX
         DO 210 NBJ = 1,NAI
C
            XAIBJ = OMEGA2(NAI,NBJ) + OMEGA2(NBJ,NAI)
            OMEGA2(NAI,NBJ) = XAIBJ
            OMEGA2(NBJ,NAI) = XAIBJ
C
  210    CONTINUE
  200 CONTINUE
C
      RETURN
      END
C  /* Deck ccsdt_tran3 */
      SUBROUTINE CCSDT_TRAN3(XINT1,XINT2,XLAMDP,XLAMDH,AOINT,IDEL)
C
C
C
#include <implicit.h>
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
      DIMENSION XINT1(NT1AMX,NVIRT,NVIRT), XINT2(NT1AMX,NRHFT,NRHFT)
      DIMENSION XLAMDP(NBAST,NORBT), XLAMDH(NBAST,NORBT)
      DIMENSION AOINT(NNBAST,NBAST)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO 100 G = 1,NBAST
         DO 110 IB = 1,NBAST
            DO 120 A = 1,NBAST
               NAB = INDEX(A,IB)
C
               if (aoint(nab,g) .eq. 0.0d0) goto 120
               DO 200 D = 1,NVIRT
                  DO 210 B = 1,NVIRT
                     DO 220 K = 1,NRHFT
                        DO 230 C = 1,NVIRT
C
                           NCK = NVIRT*(K-1) + C
C
                           XINT1(NCK,B,D) = XINT1(NCK,B,D)
     *               + AOINT(NAB,G)*XLAMDP(A,NRHFT+C)*XLAMDH(IB,K)
     *                      *XLAMDP(G,NRHFT+B)*XLAMDH(IDEL,NRHFT+D)
C
  230                   CONTINUE
  220                CONTINUE
  210             CONTINUE
  200          CONTINUE
C
               DO 300 J = 1,NRHFT
                  DO 310 L = 1,NRHFT
                     DO 320 K = 1,NRHFT
                        DO 330 C = 1,NVIRT
C
                           NCK = NVIRT*(K-1) + C
C
                           XINT2(NCK,L,J) = XINT2(NCK,L,J)
     *                  + AOINT(NAB,G)*XLAMDP(A,NRHFT+C)*XLAMDH(IB,K)
     *                                *XLAMDP(G,L)*XLAMDH(IDEL,J)
C
  330                   CONTINUE
  320                CONTINUE
  310             CONTINUE
  300          CONTINUE
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
*---------------------------------------------------------------------*
C  /* Deck ccsdt_t03am */
*=====================================================================*
      SUBROUTINE CCSDT_T03AM(T3AM,XINT1S,XINT2S,T2AM,SCR1,FOCKD,
     *                       FCKFLD,T3AM2)
*---------------------------------------------------------------------*
*
*     Purpose: compute zero-order triples cluster amplitudes
*
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccfield.h>
#include <maxorb.h>
#include <ccorb.h>

#if defined (SYS_CRAY)
      REAL XINT1S(NT1AMX,NVIRT,NVIRT)
      REAL XINT2S(NT1AMX,NRHFT,NRHFT)
      REAL T3AM(NT1AMX,NT1AMX,NT1AMX),SCR1(NT1AMX),FOCKD(NORBT)
      REAL T2AM(NT1AMX,NT1AMX),T3AM2(NT1AMX,NT1AMX,NT1AMX)
      REAL FCKFLD(NORBT,NORBT)
#else
      DOUBLE PRECISION XINT1S(NT1AMX,NVIRT,NVIRT) 
      DOUBLE PRECISION XINT2S(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION T3AM(NT1AMX,NT1AMX,NT1AMX),SCR1(NT1AMX)
      DOUBLE PRECISION T2AM(NT1AMX,NT1AMX),T3AM2(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION FCKFLD(NORBT,NORBT),FOCKD(NORBT)
#endif
 
      CALL DZERO(T3AM,NT1AMX*NT1AMX*NT1AMX)

      CALL CCSDT_T3AM_R(T3AM,0.0D0,XINT1S,XINT2S,T2AM,SCR1,FOCKD,
     *                  NONHF,FCKFLD,.FALSE.)

      CALL CCSDT_3AM(T3AM,0.0D0,SCR1,FOCKD,NONHF,FCKFLD,.FALSE.,T3AM2)

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_T03AM                          *
*---------------------------------------------------------------------*
*---------------------------------------------------------------------*
C  /* Deck ccsdt_t3am_r */
*=====================================================================*
      SUBROUTINE CCSDT_T3AM_R(T3AM,FREQ,XINT1,XINT2,T2AM,SCR1,FOCKD,
     *                        FCKINCL,FCKFLD,DIVIDE)
*---------------------------------------------------------------------*
*
*     Purpose: compute contributions to triples part of the coupled
*              cluster vector function (but without <mu_3|[F,T3]|HF>)
*
*         FCKINCL=.true.: include external field contribution
*                           <mu_3|[[V,T2],T2]|HF>
*
*---------------------------------------------------------------------*
#include <implicit.h>
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
      DIMENSION XINT1(NT1AMX,NVIRT,NVIRT), XINT2(NT1AMX,NRHFT,NRHFT)
      DIMENSION T3AM(NT1AMX,NT1AMX,NT1AMX),SCR1(NT1AMX),FOCKD(NORBT)
      DIMENSION T2AM(NT1AMX,NT1AMX)
      DIMENSION FCKFLD(NORBT,NORBT)
C 
      LOGICAL DIVIDE, FCKINCL
C
      PARAMETER (HALF = 0.5D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO 50 I = 1,NRHFT
         DO 60 A = 1,NVIRT
            NAI = NVIRT*(I-1) + A
            SCR1(NAI) = FOCKD(NRHFT+A) - FOCKD(I)
   60    CONTINUE
   50 CONTINUE
C
      DO 100 NCK = 1,NT1AMX
C
         DO 110 J = 1,NRHFT
            DO 120 B = 1,NVIRT
C
               NBJ = NVIRT*(J-1) + B
C
               DO 130 NAI = 1,NT1AMX
C
                  AIBJCK = 0.0D0
                  DO 140 D = 1,NVIRT
C
                     NDJ = NVIRT*(J-1) + D
C
                     AIBJCK = AIBJCK + XINT1(NCK,B,D)*T2AM(NDJ,NAI)
C
  140             CONTINUE
C
                  DO 150 L = 1,NRHFT
C
                     NBL = NVIRT*(L-1) + B
C
                     AIBJCK = AIBJCK - XINT2(NCK,L,J)*T2AM(NBL,NAI)
C
  150             CONTINUE
C
                  IF (FCKINCL) THEN
                     DO D = 1, NVIRT
                        NDJ = NVIRT*(J-1) + D
                        DO L = 1, NRHFT
                           NBL = NVIRT*(L-1) + B
                           AIBJCK = AIBJCK
     *                            - T2AM(NAI,NBL)
     *                             *T2AM(NCK,NDJ)
     *                             *FCKFLD(L,NRHFT+D)
                        ENDDO
                     ENDDO
                  ENDIF

                  IF (DIVIDE) THEN
                   AIBJCK = AIBJCK/(SCR1(NAI)+SCR1(NBJ)+SCR1(NCK)-FREQ)
                  END IF
C
                  T3AM(NAI,NBJ,NCK) = T3AM(NAI,NBJ,NCK) + AIBJCK
                  T3AM(NAI,NCK,NBJ) = T3AM(NAI,NCK,NBJ) + AIBJCK
                  T3AM(NBJ,NAI,NCK) = T3AM(NBJ,NAI,NCK) + AIBJCK
                  T3AM(NCK,NAI,NBJ) = T3AM(NCK,NAI,NBJ) + AIBJCK
                  T3AM(NBJ,NCK,NAI) = T3AM(NBJ,NCK,NAI) + AIBJCK
                  T3AM(NCK,NBJ,NAI) = T3AM(NCK,NBJ,NAI) + AIBJCK
C
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE

C----------------------------------------
C     Remove the forbidden elements.
C----------------------------------------
  
      CALL CCSDT_CLEAN_T3(T3AM,NT1AMX,NVIRT,NRHFT)
 
      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_T3AM_R                         *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_3AM(T3AM,FREQ,SCR1,FOCKD,
     &                     NONHF,FIELD,TRANS,T3AM2)
*---------------------------------------------------------------------*
*
*    Purpose: solve triples equations
*
*        NONHF=.FALSE.  canonical case:
*                 divide triples vector by orbital energy 
*                 differences and a frequency
*
*                 T_nu_3 = T_nu_3 / (eps_nu_3 - freq)
*
*                 FOCKD diagonal of Fock matrix (eps_i)
*                 FIELD, T3AM2 are not used in this case 
*
*        NONHF=.TRUE.  non-canonical case (external fields):
*                      solve triples equations iteratively
*               
*                  0 = R_nu_3 - (eps_nu_3-freq) T_nu_3 - [V,T]_nu_3
*
*                  --> T3AM on input  R_nu_3 
*                           on output T_nu_3 
*                      FOCKD diagonal of Fock matrix
*                      FIELD the potential V
*                      T3AM2 scratch array for triples
*                      TRANS flag, if set use V transposed
*                            (needed for "left" equations)
*        
*     Written by Christof Haettig, April 2002 
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <dummy.h>

      LOGICAL NONHF, TRANS, LOCDBG
      PARAMETER ( LOCDBG = .FALSE. )

#if defined (SYS_CRAY)
      REAL FOCKD(NORBT), SCR1(NT1AMX)
      REAL T3AM(NT1AMX,NT1AMX,NT1AMX), FREQ
      REAL AIBJCK, DDOT, XNORM, HALF
      REAL T3AM2(NT1AMX,NT1AMX,NT1AMX), FIELD(NORBT,NORBT)
#else
      DOUBLE PRECISION FOCKD(NORBT), SCR1(NT1AMX)
      DOUBLE PRECISION T3AM(NT1AMX,NT1AMX,NT1AMX), FREQ
      DOUBLE PRECISION T3AM2(NT1AMX,NT1AMX,NT1AMX), FIELD(NORBT,NORBT)
      DOUBLE PRECISION AIBJCK, DDOT, XNORM, HALF
#endif
      PARAMETER ( HALF = 0.5D0 )

      LOGICAL NONCONV
      INTEGER NI, NA, NAI, NBJ, NCK, NDJ, NBL,
     &        LUTEMP, MAXITE

      DO NI = 1,NRHFT
         DO NA = 1,NVIRT
            NAI = NVIRT*(NI-1) + NA
            SCR1(NAI) = FOCKD(NRHFT+NA) - FOCKD(NI)
         END DO
      END DO
 
      DO NCK = 1, NT1AMX
        DO NBJ = 1, NT1AMX
          DO NAI = 1, NT1AMX
            T3AM(NAI,NBJ,NCK) = T3AM(NAI,NBJ,NCK) /
     &         (SCR1(NAI) + SCR1(NBJ) + SCR1(NCK) - FREQ)
          END DO
        END DO
      END DO

      CALL CCSDT_CLEAN_T3(T3AM,NT1AMX,NVIRT,NRHFT)

      IF (NONHF) THEN
         LUTEMP = -1
         CALL GPOPEN(LUTEMP,'CC3AMTMP','UNKNOWN',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         REWIND LUTEMP
         WRITE (LUTEMP) (((T3AM(NAI,NBJ,NCK),NAI=1,NT1AMX),
     &                    NBJ=1,NT1AMX),NCK=1,NT1AMX)

         CALL DCOPY(NT1AMX*NT1AMX*NT1AMX,T3AM,1,T3AM2,1)

         MAXITE  = 1
         NONCONV = .TRUE.

         DO WHILE (NONCONV)

           DO NCK = 1, NT1AMX
            DO J = 1, NRHFT
             DO B = 1, NVIRT
              NBJ = NVIRT*(J-1) + B
              DO NAI = 1, NT1AMX
                AIBJCK = 0.0D0
                IF (TRANS) THEN
                  DO D = 1, NVIRT
                    NDJ = NVIRT*(J-1) + D
                    AIBJCK = AIBJCK - HALF*T3AM2(NAI,NDJ,NCK)
     *                           *FIELD(NRHFT+D,NRHFT+B)
                  ENDDO
                  DO L = 1, NRHFT
                    NBL = NVIRT*(L-1) + B
                    AIBJCK = AIBJCK + HALF*T3AM2(NAI,NBL,NCK)
     *                        *FIELD(J,L)
                  ENDDO
                ELSE 
                  DO D = 1, NVIRT
                    NDJ = NVIRT*(J-1) + D
                    AIBJCK = AIBJCK - HALF*T3AM2(NAI,NDJ,NCK)
     *                           *FIELD(NRHFT+B,NRHFT+D)
                  ENDDO
                  DO L = 1, NRHFT
                    NBL = NVIRT*(L-1) + B
                    AIBJCK = AIBJCK + HALF*T3AM2(NAI,NBL,NCK)
     *                        *FIELD(L,J)
                  ENDDO
                END IF
                AIBJCK = AIBJCK/(SCR1(NAI)+SCR1(NBJ)+SCR1(NCK)-FREQ)
                T3AM(NAI,NBJ,NCK) = T3AM(NAI,NBJ,NCK) + AIBJCK
                T3AM(NAI,NCK,NBJ) = T3AM(NAI,NCK,NBJ) + AIBJCK
                T3AM(NBJ,NAI,NCK) = T3AM(NBJ,NAI,NCK) + AIBJCK
                T3AM(NCK,NAI,NBJ) = T3AM(NCK,NAI,NBJ) + AIBJCK
                T3AM(NBJ,NCK,NAI) = T3AM(NBJ,NCK,NAI) + AIBJCK
                T3AM(NCK,NBJ,NAI) = T3AM(NCK,NBJ,NAI) + AIBJCK
              ENDDO
             ENDDO
            ENDDO
           ENDDO

           CALL CCSDT_CLEAN_T3(T3AM,NT1AMX,NVIRT,NRHFT)

           CALL DAXPY(NT1AMX*NT1AMX*NT1AMX,-1.0D0,T3AM,1,T3AM2,1)
           XNORM = DSQRT(DDOT(NT1AMX*NT1AMX*NT1AMX,T3AM2,1,T3AM2,1))

           IF (LOCDBG.OR.DEBUG) 
     &       WRITE(LUPRI,*)'CCSDT_3AM> Norm for iteration ',maxite,xnorm
           IF (XNORM .LT. 1.0D-14) THEN
             NONCONV = .FALSE.
             IF (LOCDBG.OR.DEBUG) WRITE(LUPRI,*)
     &         'CCSDT_3AM> converged in iteration',maxite
           ELSE
             CALL DCOPY(NT1AMX*NT1AMX*NT1AMX,T3AM,1,T3AM2,1)
             IF(LOCDBG.OR.DEBUG)WRITE(LUPRI,*)'present norm of t3am:',
     &        DSQRT(DDOT(NT1AMX*NT1AMX*NT1AMX,T3AM,1,T3AM,1)) 
             MAXITE = MAXITE + 1
             IF (MAXITE .GT. 30)
     &            CALL QUIT('MAX ITERATIONS EXCEEDED in CCSDT_3AM')
             REWIND LUTEMP
             READ (LUTEMP) (((T3AM(NAI,NBJ,NCK),NAI=1,NT1AMX),
     &                        NBJ=1,NT1AMX),NCK=1,NT1AMX)
           ENDIF

         ENDDO

         CALL GPCLOSE(LUTEMP,'KEEP')
      END IF

      RETURN
      END 
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_3AM                            *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_CLEAN_T3(T3AM,NT1AMX,NVIRT,NRHFT)
*---------------------------------------------------------------------*
*     Purpose: remove forbidden elements in a triples vector
*---------------------------------------------------------------------*
      IMPLICIT NONE  

      INTEGER NT1AMX, NVIRT, NRHFT

#if defined (SYS_CRAY)
      REAL T3AM(NT1AMX,NT1AMX,NT1AMX)
#else
      DOUBLE PRECISION T3AM(NT1AMX,NT1AMX,NT1AMX)
#endif

      INTEGER NA, NI, NAI, NB, NBI, NC, NCI, NJ, NAJ, NK, NAK
    
      DO NA = 1, NVIRT
       DO NI = 1, NRHFT
         NAI = NVIRT*(NI-1) + NA
         DO NB = 1, NVIRT
            NBI = NVIRT*(NI-1) + NB
            DO NC = 1, NVIRT
               NCI = NVIRT*(NI-1) + NC
               T3AM(NAI,NBI,NCI) = 0.0d0
               T3AM(NAI,NCI,NBI) = 0.0d0
               T3AM(NBI,NAI,NCI) = 0.0d0
               T3AM(NBI,NCI,NAI) = 0.0d0
               T3AM(NCI,NAI,NBI) = 0.0d0
               T3AM(NCI,NBI,NAI) = 0.0d0
            ENDDO
         ENDDO
         DO NJ = 1, NRHFT
            NAJ = NVIRT*(NJ-1) + NA
            DO NK = 1, NRHFT
               NAK = NVIRT*(NK-1) + NA
               T3AM(NAI,NAJ,NAK) = 0.0d0
               T3AM(NAI,NAK,NAJ) = 0.0d0
               T3AM(NAJ,NAI,NAK) = 0.0d0
               T3AM(NAJ,NAK,NAI) = 0.0d0
               T3AM(NAK,NAI,NAJ) = 0.0d0
               T3AM(NAK,NAJ,NAI) = 0.0d0
            ENDDO
         ENDDO
       ENDDO
      ENDDO
 
      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_CLEAN_T3                       *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE PRINT_PT3_NODDY(T3AM)
*---------------------------------------------------------------------*
*
*    Purpose: Print a T3AM triples vector from a noddy code routine
*             in such a way that it can directly be compared 
*             to the printout from PRINT_PT3 routine.
*
*     Written by Filip Pawlowski, 26-09-2002 
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER NAI, NBJ, NCK
C
#if defined (SYS_CRAY)
      REAL T3AM(NT1AMX,NT1AMX,NT1AMX)
#else
      DOUBLE PRECISION T3AM(NT1AMX,NT1AMX,NT1AMX)
#endif
C
      DO A = 1, NVIRT
        DO B = 1, NVIRT
          DO C = 1, NVIRT
            DO I = 1, NRHFT
              DO J = 1, NRHFT
                DO K = 1, NRHFT
C
                 NAI = NVIRT*(I-1) + A
                 NBJ = NVIRT*(J-1) + B
                 NCK = NVIRT*(K-1) + C
C
                 IF (ABS(T3AM(NAI,NBJ,NCK)) .GT. 1.0D-12) THEN
C
                     WRITE(LUPRI,1) 'noddy T3AM(',A,',',B,',',
     *                                            C,',',I,',',
     *                                            J,',',K,') = ',
     *               T3AM(NAI,NBJ,NCK)
C
                 END IF

                END DO 
              END DO 
            END DO 
          END DO 
        END DO 
      END DO 
C
    1 FORMAT(1X,A11,I3,A1,I3,A1,I3,A1,I3,A1,I3,A1,I3,A4,E20.10)
C
      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE PRINT_PT3_NODDY                      *
*---------------------------------------------------------------------*
