C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
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
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <maxash.h>
#include <mxorb.h>
#include <aovec.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION INDEXA(MXCORB)
      DIMENSION OMEGA1(*),OMEGA2(*),T1AM(*),T2AM(*),FOCK(*)
      DIMENSION XLAMDP(*),XLAMDH(*),WORK(LWORK)
C
#include <inforb.h>
#include <infind.h>
#include <blocks.h>
#include <inftap.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccsdio.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      IF (NSYM .NE. 1) THEN
         CALL QUIT('No symmetry in this part yet')
      ENDIF
C
C------------------------
C     Dynamic Allocation.
C------------------------
C
      KCMO   = 1
      KSCR1  = KCMO   + NLAMDT
      KFOCKD = KSCR1  + NT1AMX
      KINT1  = KFOCKD + NORBT
      KINT2  = KINT1  + NT1AMX*NVIRT*NVIRT
      KINT1T = KINT2  + NRHFT*NRHFT*NT1AMX
      KINT2T = KINT1T + NT1AMX*NVIRT*NVIRT
      KINT1S = KINT2T + NRHFT*NRHFT*NT1AMX
      KINT2S = KINT1S + NT1AMX*NVIRT*NVIRT
      KT3AM  = KINT2S + NRHFT*NRHFT*NT1AMX
      KOME1  = KT3AM  + NT1AMX*NT1AMX*NT1AMX
      KOME2  = KOME1  + NT1AMX
      KIAJB  = KOME2  + NT1AMX*NT1AMX
      KSIG2  = KOME2  + NT1AMX*NT1AMX
      KYIAJB = KSIG2  + NT1AMX*NT1AMX
      KEND1  = KYIAJB + NT1AMX*NT1AMX
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSD_TRIPLE')
      ENDIF
C
C--------------------------------
C     Initialize integral arrays.
C--------------------------------
C
      CALL DZERO(WORK(KINT1),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2),NT1AMX*NRHFT*NRHFT)
      CALL DZERO(WORK(KINT1T),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2T),NT1AMX*NRHFT*NRHFT)
      CALL DZERO(WORK(KINT1S),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2S),NT1AMX*NRHFT*NRHFT)
      CALL DZERO(WORK(KT3AM),NT1AMX*NT1AMX*NT1AMX)
      CALL DZERO(WORK(KOME1),NT1AMX)
      CALL DZERO(WORK(KOME2),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KIAJB),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KYIAJB),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KSIG2),NT1AMX*NT1AMX)
C
C----------------------------------------------
C     Read MO-coefficients from interface file.
C----------------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
C
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KCMO+I-1), I=1,NLAMDT)
C
C-------------------------------------
C     Read canonical orbital energies.
C-------------------------------------
C
      REWIND LUSIFC
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBT)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C====================================================
C     Start the loop over distributions of integrals.
C====================================================
C
      IF (DIRECT) THEN
         NTOSYM = 1
         DTIME  = SECOND()
c        CALL HERDI1(WORK(KEND1),LWRK1,IPRINT)
         CALL QUIT('Direct not implemented')
         DTIME  = SECOND() - DTIME
         TIMHER1 = TIMHER1 + DTIME
      ELSE
         NTOSYM = NSYM
      ENDIF
C
      DO 100 ISYMD1 = 1,NTOSYM
C
         IF (DIRECT) THEN
            NTOT = MAXSHL
         ELSE
            NTOT = NBAS(ISYMD1)
         ENDIF
C
         DO 110 ILLL = 1,NTOT
C
C-----------------------------------------------------------------
C           If direct calculate the integrals and transposed t2am.
C-----------------------------------------------------------------
C
            IF (DIRECT) THEN
               DTIME  = SECOND()
c              CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,IPRINT)
               CALL QUIT('Direct not implemented')
               DTIME  = SECOND() - DTIME
               TIMHER2 = TIMHER2 + DTIME
            ELSE
               NUMDIS = 1
            ENDIF
C
C-----------------------------------------------------
C           Loop over number of distributions in disk.
C-----------------------------------------------------
C
            DO 120 IDEL2 = 1,NUMDIS
C
               IF (DIRECT) THEN
                  IDEL  = INDEXA(IDEL2)
                  ISYMD = ISAO(IDEL)
               ELSE
                  IDEL  = IBAS(ISYMD1) + ILLL
                  ISYMD = ISYMD1
               ENDIF
C
               ISYDIS = MULD2H(ISYMD,ISYMOP)
C
C------------------------------------------
C              Work space allocation no. 2.
C------------------------------------------
C
               KXINT  = KEND1
               KEND2  = KXINT + NDISAO(ISYDIS)
               LWRK2  = LWORK - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
                  CALL QUIT('Insufficient space in CCSD_TRIPLE')
               ENDIF
C
C
C-----------------------------------------
C              Read in batch of integrals.
C-----------------------------------------
C
               DTIME   = SECOND()
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                     DIRECT)
               DTIME   = SECOND() - DTIME
               TIMRDAO = TIMRDAO  + DTIME
C
C----------------------------------------------------
C              Calculate integrals needed in CCSD[T].
C----------------------------------------------------
C
c              CALL CCSDT_TRAN1(WORK(KINT1),WORK(KINT2),WORK(KCMO),
c    *                          WORK(KCMO),WORK(KXINT),IDEL)
C
               CALL CCSDT_TRAN1(WORK(KINT1T),WORK(KINT2T),XLAMDP,
     *                          XLAMDH,WORK(KXINT),IDEL)
C
               CALL CCSDT_TRAN2(WORK(KIAJB),WORK(KYIAJB),WORK(KCMO),
     *                          WORK(KXINT),IDEL)
C
               CALL CCSDT_TRAN3(WORK(KINT1S),WORK(KINT2S),XLAMDP,
     *                          XLAMDH,WORK(KXINT),IDEL)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C-------------------------------------
C     Calculate the triple amplitudes.
C-------------------------------------
C
      CALL CCSDT_T3AM(WORK(KT3AM),WORK(KINT1S),WORK(KINT2S),T2AM,
     *                WORK(KSCR1),WORK(KFOCKD))
C
C------------------------------------
C     Calculate the CCSD corrections.
C------------------------------------
C
      CALL CCSDT_OMEGA1(WORK(KOME1),WORK(KIAJB),WORK(KT3AM))
C
      CALL CCSDT_OMEGA2(WORK(KOME2),WORK(KINT1T),WORK(KINT2T),
     *                  WORK(KT3AM),FOCK)
C
      DO 300 I = 1,NT1AMX
         OMEGA1(I) = OMEGA1(I) + WORK(KOME1+I-1)
  300 CONTINUE
      ENERT1 = ddot(nt1amx,t1am,1,work(kome1),1)
      write(LUPRI,*) 'The E5(ST) contribution :',ENERT1
      CALL ECCSDT(T2AM,WORK(KOME2))

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
      SUBROUTINE ECCSDT(T2AM,OMEGA2)
C
C
C
#include <implicit.h>
      PARAMETER (TWO=2.0d0)
      DIMENSION T2AM(NT1AMX,NT1AMX),OMEGA2(NT1AMX,NT1AMX)
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      ENERGY = 0.0d0
      DO 100 I = 1,NRHFT
         DO 110 A = 1,NVIRT
            NAI = NVIRT*(I-1) + A
            DO 120 J = 1,NRHFT
               NAJ = NVIRT*(J-1) + A
               DO 130 B = 1,NVIRT
                  NBJ = NVIRT*(J-1) + B
                  NBI = NVIRT*(I-1) + B
                  DELTA = OMEGA2(NAI,NBJ)
     *                   *(TWO*T2AM(NAI,NBJ) - T2AM(NAJ,NBI))
                  ENERGY = ENERGY + DELTA
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      WRITE(LUPRI,*) 'Energy contribution E4(TT) :',ENERGY
C
      RETURN
      END
C  /* Deck ccsdt_tran1 */
      SUBROUTINE CCSDT_TRAN1(XINT1,XINT2,XLAMDP,XLAMDH,AOINT,IDEL)
C
C
C
#include <implicit.h>
      DIMENSION XINT1(NT1AMX,NVIRT,NVIRT), XINT2(NT1AMX,NRHFT,NRHFT)
      DIMENSION XLAMDP(NBAST,NORBT), XLAMDH(NBAST,NORBT)
      DIMENSION AOINT(NNBAST,NBAST)
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
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
      PARAMETER (TWO = 2.0D0)
      DIMENSION XIAJB(NT1AMX,NT1AMX), AOINT(NNBAST,NBAST)
      DIMENSION YIAJB(NT1AMX,NT1AMX)
      DIMENSION CMO(NORBT,NORBT)
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
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
C  /* Deck ccsdt_t3am */
      SUBROUTINE CCSDT_T3AM(T3AM,XINT1,XINT2,T2AM,SCR1,FOCKD)
C
C
C
#include <implicit.h>
      DIMENSION XINT1(NT1AMX,NVIRT,NVIRT), XINT2(NT1AMX,NRHFT,NRHFT)
      DIMENSION T3AM(NT1AMX,NT1AMX,NT1AMX),SCR1(NT1AMX),FOCKD(*)
      DIMENSION T2AM(NT1AMX,NT1AMX)
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
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
                  AIBJCK = AIBJCK/(SCR1(NAI) + SCR1(NBJ) + SCR1(NCK))
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
C
      RETURN
      END
C  /* Deck ccsdt_omega1 */
      SUBROUTINE CCSDT_OMEGA1(OMEGA1,XIAJB,T3AM)
C
C
C
#include <implicit.h>
      PARAMETER (TWO = 2.0D0)
      DIMENSION OMEGA1(NT1AMX),XIAJB(NT1AMX,NT1AMX)
      DIMENSION T3AM(NT1AMX,NT1AMX,NT1AMX)
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
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
      PARAMETER (TWO = 2.0D0)
      DIMENSION XINT1T(NT1AMX,NVIRT,NVIRT)
      DIMENSION XINT2T(NT1AMX,NRHFT,NRHFT)
      DIMENSION OMEGA2(NT1AMX,NT1AMX)
      DIMENSION T3AM(NT1AMX,NT1AMX,NT1AMX),FOCK(NORBT,NORBT)
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
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
C  /* Deck ccsdt_sigma2 */
      SUBROUTINE CCSDT_SIGMA2(SIGMA2,XIAJB,YIAJB,T1AM,T3AM)
C
C
C
#include <implicit.h>
      PARAMETER (TWO = 2.0D0)
      DIMENSION T1AM(NT1AMX)
      DIMENSION XIAJB(NT1AMX,NT1AMX), YIAJB(NT1AMX,NT1AMX)
      DIMENSION SIGMA2(NT1AMX,NT1AMX)
      DIMENSION T3AM(NT1AMX,NT1AMX,NT1AMX)
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
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
C
                        DO 160 D = 1,NVIRT
C
                           NDI = NVIRT*(I-1) + D
                           NDK = NVIRT*(K-1) + D
                           NDJ = NVIRT*(J-1) + D
C
                           DO 170 L = 1,NRHFT
C
                              NCL = NVIRT*(L-1) + C
                              NDL = NVIRT*(L-1) + D
                              NAL = NVIRT*(L-1) + A
                              NBL = NVIRT*(L-1) + B
C
                   SIGMA2(NAI,NBJ) = SIGMA2(NAI,NBJ)
     *         + T3AM(NAI,NBJ,NCK)*T1AM(NDL)*XIAJB(NDL,NCK)
     *         - T3AM(NAI,NDJ,NCK)*T1AM(NBL)*XIAJB(NDL,NCK)
     *         - T3AM(NAI,NBL,NCK)*T1AM(NDJ)*XIAJB(NDL,NCK)
     *         - T3AM(NAI,NCJ,NBK)*T1AM(NDL)*XIAJB(NDL,NCK)
     *         + T3AM(NBL,NCJ,NAK)*T1AM(NDI)*YIAJB(NDK,NCL)
     *         + T3AM(NBK,NDJ,NCI)*T1AM(NAL)*YIAJB(NDK,NCL)
C
  170                      CONTINUE
C
  160                   CONTINUE
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
            ZAIBJ = SIGMA2(NAI,NBJ) + SIGMA2(NBJ,NAI)
            SIGMA2(NAI,NBJ) = ZAIBJ
            SIGMA2(NBJ,NAI) = ZAIBJ
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
      DIMENSION XINT1(NT1AMX,NVIRT,NVIRT), XINT2(NT1AMX,NRHFT,NRHFT)
      DIMENSION XLAMDP(NBAST,NORBT), XLAMDH(NBAST,NORBT)
      DIMENSION AOINT(NNBAST,NBAST)
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
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
