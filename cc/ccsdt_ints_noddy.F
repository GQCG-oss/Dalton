C  /* Deck ccsdt_tran1_r */
      SUBROUTINE CCSDT_TRAN1_R(XINT1,XINT2,XLAMDP0,XLAMDH0,
     &                         XLAMDP1,XLAMDH1,AOINT,IDEL)
C
C     XINT1 = XINT1 + (K^p0 C^h0|B^p1 D^h0)
C     XINT2 = XINT2 + (K^p0 C^h0|L^p0 J^h1)
C
C
#include <implicit.h>
      DIMENSION XINT1(NT1AMX,NVIRT,NVIRT), XINT2(NT1AMX,NRHFT,NRHFT)
      DIMENSION XLAMDP0(NBAST,NORBT), XLAMDH0(NBAST,NORBT)
      DIMENSION XLAMDP1(NBAST,NORBT), XLAMDH1(NBAST,NORBT)
      DIMENSION AOINT(NNBAST,NBAST)
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO 100 IGAM = 1,NBAST
         DO 110 IBET = 1,NBAST
            DO 120 IALP = 1,NBAST
               NAB = INDEX(IALP,IBET)
C
               if (aoint(nab,IGAM) .eq. 0.0d0) goto 120
               DO 200 D = 1,NVIRT
                  DO 210 B = 1,NVIRT
                     DO 220 K = 1,NRHFT
                        DO 230 C = 1,NVIRT
C
                           NCK = NVIRT*(K-1) + C
C
                           XINT1(NCK,B,D)=XINT1(NCK,B,D)+AOINT(NAB,IGAM)
     &                     *XLAMDP0(IBET,K)      *XLAMDH0(IALP,NRHFT+C)
     &                     *XLAMDP1(IGAM,NRHFT+B)*XLAMDH0(IDEL,NRHFT+D)
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
                           XINT2(NCK,L,J)=XINT2(NCK,L,J)+AOINT(NAB,IGAM)
     &                      * XLAMDP0(IBET,K) * XLAMDH0(IALP,NRHFT+C)
     &                      * XLAMDP0(IGAM,L) * XLAMDH1(IDEL,J)
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
C  /* Deck ccsdt_tran3_r */
      SUBROUTINE CCSDT_TRAN3_R(XINT1,XINT2,XLAMDP0,XLAMDH0,
     &                         XLAMDP1,XLAMDH1,XLAMDP2,XLAMDH2,
     &                         AOINT,IDEL)
C
C     XINT1 = XINT1 + (C^p1 K^h1|B^p2 D^h0)
C     XINT2 = XINT2 + (C^p1 K^h1|L^p0 J^h2)
C
#include <implicit.h>
      DIMENSION XINT1(NT1AMX,NVIRT,NVIRT), XINT2(NT1AMX,NRHFT,NRHFT)
      DIMENSION XLAMDP0(NBAST,NORBT), XLAMDH0(NBAST,NORBT)
      DIMENSION XLAMDP1(NBAST,NORBT), XLAMDH1(NBAST,NORBT)
      DIMENSION XLAMDP2(NBAST,NORBT), XLAMDH2(NBAST,NORBT)
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
                           XINT1(NCK,B,D) = XINT1(NCK,B,D)+AOINT(NAB,G)
     &                      * XLAMDP1(A,NRHFT+C) * XLAMDH1(IB,K)
     &                      * XLAMDP2(G,NRHFT+B) * XLAMDH0(IDEL,NRHFT+D)
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
                           XINT2(NCK,L,J) = XINT2(NCK,L,J)+AOINT(NAB,G)
     &                      * XLAMDP1(A,NRHFT+C) * XLAMDH1(IB,K)
     &                      * XLAMDP0(G,L)       * XLAMDH2(IDEL,J)
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
C  /* Deck ccfop_tran1 */
      SUBROUTINE CCFOP_TRAN1_R(XINT1,XINT2,XINT3,XINT4,
     &                         XLAMDP0,XLAMDH0,
     &                         XLAMDP1,XLAMDH1,
     &                         XLAMDP2,XLAMDH2,
     &                         AOINT,IDEL)
C
C     XINT1 = (K^p0 C^h0|D^p1 L^h2)    "(O-0 V-0|V-1 O-2)"
C     XINT2 = (K^p0 L^h1|C^p2 D^h0)    "(O-0 O-1|V-2 V-0)"
C     XINT3 = (K^p0 L^h1|M^p0 N^h2)    "(O-0 O-1|O-0 O-2)"
C     XINT4 = (C^p1 D^h0|E^p2 F^h0)    "(V-1 V-0|V-2 V-0)"
C
#include <implicit.h>
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
      DIMENSION XINT1(NRHFT,NVIRT,NVIRT,NRHFT)
      DIMENSION XINT2(NRHFT,NRHFT,NVIRT,NVIRT)
      DIMENSION XINT3(NRHFT,NRHFT,NRHFT,NRHFT)
      DIMENSION XINT4(NVIRT,NVIRT,NVIRT,NVIRT)
      DIMENSION XLAMDP0(NBAST,NORBT), XLAMDH0(NBAST,NORBT)
      DIMENSION XLAMDP1(NBAST,NORBT), XLAMDH1(NBAST,NORBT)
      DIMENSION XLAMDP2(NBAST,NORBT), XLAMDH2(NBAST,NORBT)
      DIMENSION AOINT(NNBAST,NBAST)
C
      LOGICAL LDEBUG
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      LDEBUG = .TRUE.
C
C----------------------------------------
C     Calculate integrals :
C----------------------------------------
C
      DO 100 G = 1,NBAST
         DO 110 IB = 1,NBAST
            DO 120 A = 1,NBAST
               NAB = INDEX(A,IB)
C
               if (aoint(nab,g) .eq. 0.0d0) goto 120
C
               DO NC = 1,NVIRT
                  DO ND = 1,NVIRT
                     DO NK = 1,NRHFT
                        DO NL = 1,NRHFT
C
                           XINT1(NK,NC,ND,NL) = XINT1(NK,NC,ND,NL)
     &                       + AOINT(NAB,G) *
     &                          XLAMDP0(A,NK) * 
     &                          XLAMDH0(IB,NRHFT+NC) *
     &                          XLAMDP1(G,NRHFT+ND) * 
     &                          XLAMDH2(IDEL,NL)
C
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
C
               DO NC = 1,NVIRT
                  DO ND = 1,NVIRT
                     DO NK = 1,NRHFT
                        DO NL = 1,NRHFT
C
                           XINT2(NK,NL,NC,ND) = XINT2(NK,NL,NC,ND)
     &                       + AOINT(NAB,G) *
     &                          XLAMDP0(A,NK) * 
     &                          XLAMDH1(IB,NL) *
     &                          XLAMDP2(G,NRHFT+NC) *
     &                          XLAMDH0(IDEL,NRHFT+ND)
C
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
C
               DO NK = 1,NRHFT
                  DO NL = 1,NRHFT
                     DO NM = 1,NRHFT
                        DO NN = 1,NRHFT
C
                           XINT3(NK,NL,NM,NN) = XINT3(NK,NL,NM,NN)
     &                       + AOINT(NAB,G) *
     &                          XLAMDP0(A,NK) *
     &                          XLAMDH1(IB,NL) *
     &                          XLAMDP0(G,NM) *
     &                          XLAMDH2(IDEL,NN)
C
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
C
               DO NC = 1,NVIRT
                  DO ND = 1,NVIRT
                     DO NE = 1,NVIRT
                        DO NF = 1,NVIRT
C
                           XINT4(NC,ND,NE,NF) = XINT4(NC,ND,NE,NF)
     &                       + AOINT(NAB,G) * 
     &                          XLAMDP1(A,NRHFT+NC) *
     &                          XLAMDH0(IB,NRHFT+ND) *
     &                          XLAMDP2(G,NRHFT+NE) *
     &                          XLAMDH0(IDEL,NRHFT+NF)
C
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
