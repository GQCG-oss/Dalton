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
#include <single.h>
C
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
*=====================================================================*
      SUBROUTINE CCSDT_INTS0_NODDY(DO_INTXY,XIAJB,YIAJB,
     &                             DO_INTS0,XINT1S,XINT2S,
     &                             DO_INTT0,XINT1T,XINT2T,
     &                             XLAMDP,XLAMDH,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: compute some standard integrals for CC3
*
*             if (do_intxy) compute 
*                 XIAJB  = 2(kc|ld) - (kd|lc)
*                 YIAJB  =  (kc|ld)
*
*             if (do_ints0) compute 
*                 XINT1S = (ck|bd) and 
*                 XINT2S = (ck|lj)
*
*             if (do_intt0) compute 
*                 XINT1T = (kc|bd) and 
*                 XINT2T = (kc|lj)
*
* Written by Christof Haettig, November 2002, based on CCSDT_XI3_NODDY.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccfield.h>
#include <ccorb.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG=.FALSE.)
      INTEGER ISYM0
      PARAMETER (ISYM0=1)

      LOGICAL DO_INTXY, DO_INTS0, DO_INTT0
      INTEGER LWORK
#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL XINT1T(NT1AMX*NVIRT*NVIRT)
      REAL XINT2T(NRHFT*NRHFT*NT1AMX)
      REAL XINT1S(NT1AMX*NVIRT*NVIRT)
      REAL XINT2S(NRHFT*NRHFT*NT1AMX)
      REAL XIAJB(NT1AMX*NT1AMX)
      REAL YIAJB(NT1AMX*NT1AMX)
      REAL XLAMDP(*), XLAMDH(*)
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION XINT1T(NT1AMX*NVIRT*NVIRT)
      DOUBLE PRECISION XINT2T(NRHFT*NRHFT*NT1AMX)
      DOUBLE PRECISION XINT1S(NT1AMX*NVIRT*NVIRT)
      DOUBLE PRECISION XINT2S(NRHFT*NRHFT*NT1AMX)
      DOUBLE PRECISION XIAJB(NT1AMX*NT1AMX)
      DOUBLE PRECISION YIAJB(NT1AMX*NT1AMX)
      DOUBLE PRECISION XLAMDP(*), XLAMDH(*)
#endif

      INTEGER ISYMD, ILLL, IDEL, ISYDIS, KXINT, KEND1, LWRK1, KRECNR

      CALL QENTER('CCSDT_INTS0_NODDY')

      IF (DIRECT) 
     &  CALL QUIT('DIRECT NOT IMPLEMENTED IN CCSDT_INTS0_NODDY')

*---------------------------------------------------------------------*
*     Loop over distributions of integrals:
*---------------------------------------------------------------------*
      IF (DO_INTS0) THEN
        CALL DZERO(XINT1S,NT1AMX*NVIRT*NVIRT)
        CALL DZERO(XINT2S,NT1AMX*NRHFT*NRHFT)
      END IF
      IF (DO_INTT0) THEN
         CALL DZERO(XINT1T,NT1AMX*NVIRT*NVIRT)
         CALL DZERO(XINT2T,NT1AMX*NRHFT*NRHFT)
      END IF
      IF (DO_INTXY) THEN
         CALL DZERO(XIAJB,NT1AMX*NT1AMX)
         CALL DZERO(YIAJB,NT1AMX*NT1AMX)
      END IF

      DO ISYMD = 1, NSYM
         DO ILLL = 1,NBAS(ISYMD)
            IDEL   = IBAS(ISYMD) + ILLL
            ISYDIS = MULD2H(ISYMD,ISYM0)
 
C           ----------------------------
C           Work space allocation no. 2.
C           ----------------------------
            KXINT  = 1
            KEND1  = KXINT + NDISAO(ISYDIS)
            LWRK1  = LWORK - KEND1
            IF (LWRK1 .LT. 0) THEN
               WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
               CALL QUIT('Insufficient space in CCSDT_INTS0_NODDY')
            ENDIF
 
C           ---------------------------
C           Read in batch of integrals.
C           ---------------------------
            CALL CCRDAO(WORK(KXINT),IDEL,1,WORK(KEND1),LWRK1,
     *                  WORK(KRECNR),DIRECT)
 
C           ----------------------------------
C           Calculate integrals needed in CC3:
C           ----------------------------------
            IF (DO_INTT0) THEN
               CALL CCSDT_TRAN1(XINT1T,XINT2T,XLAMDP,XLAMDH,
     *                          WORK(KXINT),IDEL)
            END IF

            IF (DO_INTXY) THEN
               CALL CC3_TRAN2(XIAJB,YIAJB,XLAMDP,XLAMDH,
     *                        WORK(KXINT),IDEL)
            END IF
 
            IF (DO_INTS0) THEN
              CALL CCSDT_TRAN3(XINT1S,XINT2S,XLAMDP,XLAMDH,
     *                         WORK(KXINT),IDEL)
            END IF
 
         END DO   
      END DO  
*---------------------------------------------------------------------*
*     End Loop over distributions of integrals.
*---------------------------------------------------------------------*
 
      CALL QEXIT('CCSDT_INTS0_NODDY')
      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_INTS0_NODDY                    *
*---------------------------------------------------------------------*
