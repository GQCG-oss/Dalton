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
*---------------------------------------------------------------------*
C  /* Deck cc_lhtr_noddy */
*---------------------------------------------------------------------*
      SUBROUTINE CC_LHTR_NODDY(FREQ,OMEGA1,OMEGA2,T1AM,T2AM,C1AM,C2AM,
     *                         XLAMDP,XLAMDH,WORK,LWORK,IOPT)
*=====================================================================*
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
C     IOPT = 0 No storage of T3, L3
C     IOPT = 1 Store T3 and L3 in 'T3AMPL' and 'L3AMPL'
C
C
*---------------------------------------------------------------------*
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <maxash.h>
#include <aovec.h>
#include <inforb.h>
#include <infind.h>
#include <blocks.h>
#include <inftap.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <mxcent.h>
#include <eribuf.h>
#include <ccfield.h>
#include <ccnoddy.h>
C
      LOGICAL LOCDBG
      PARAMETER(LOCDBG = .FALSE.)
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION OMEGA1(*),OMEGA2(*),T1AM(*),T2AM(*), C1AM(*), C2AM(*)
      DIMENSION XLAMDP(*),XLAMDH(*),WORK(LWORK)
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
      KSCR1  = 1
      KFOCK  = KSCR1  + NT1AMX
      KFOCKD = KFOCK  + N2BST(ISYMOP)
      KEND1  = KFOCKD + NORBT

      IF ((NONHF) .AND. (NFIELD .GT. 0)) THEN
         KONEP  = KEND1
         KONEPAO= KONEP + N2BST(ISYMOP)
         KEND1  = KONEPAO+N2BST(ISYMOP)
         LWRK1  = LWORK  - KEND1
      ENDIF

      KINT1S = KEND1
      KINT2S = KINT1S + NT1AMX*NVIRT*NVIRT
      KEND1  = KINT2S + NRHFT*NRHFT*NT1AMX

      KIOVVO = KEND1
      KIOOVV = KIOVVO + NRHFT*NVIRT*NVIRT*NRHFT
      KIOOOO = KIOOVV + NRHFT*NRHFT*NVIRT*NVIRT
      KIVVVV = KIOOOO + NRHFT*NRHFT*NRHFT*NRHFT
      KEND1  = KIVVVV + NVIRT*NVIRT*NVIRT*NVIRT

      KL3AM  = KEND1
      KEND1  = KL3AM  + NT1AMX*NT1AMX*NT1AMX

      ! what is above has to be kept until the end...
      ! everything below might be overwritten in CC_LHPART_NODDY
      KEND1A  = KEND1
      LWRK1A  = LWORK  - KEND1A  

      KINT1T = KEND1 
      KINT2T = KINT1T + NT1AMX*NVIRT*NVIRT
      KEND1  = KINT2T + NRHFT*NRHFT*NT1AMX

      KT3AM  = KEND1 
      KEND1  = KT3AM  + NT1AMX*NT1AMX*NT1AMX

      KOME1  = KEND1 
      KOME2  = KOME1  + NT1AMX
      KEND1  = KOME2  + NT1AMX*NT1AMX

      KIAJB  = KEND1 
      KYIAJB = KIAJB  + NT1AMX*NT1AMX
      KEND1  = KYIAJB + NT1AMX*NT1AMX

      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CC_LHTR_NODDY')
      ENDIF
C
C--------------------------------
C     Initialize integral arrays.
C--------------------------------
C
      CALL DZERO(WORK(KL3AM),NT1AMX*NT1AMX*NT1AMX)

      CALL CCSDT_READ_NODDY(.TRUE.,WORK(KFOCKD),WORK(KFOCK),
     &                             WORK(KONEP), WORK(KONEPAO),
     &                      .TRUE.,WORK(KIAJB), WORK(KYIAJB),
     &                      .TRUE.,WORK(KINT1S),WORK(KINT2S),
     &                      .TRUE.,WORK(KINT1T),WORK(KINT2T),
     &                      .TRUE.,WORK(KIOVVO),WORK(KIOOVV),
     &                             WORK(KIOOOO),WORK(KIVVVV),
     &                      NORBT,NLAMDT,NRHFT,NVIRT,NT1AMX)
C
C-------------------------------------
C     Get the T3 amplitudes.
C-------------------------------------
      LUTEMP = -1
      CALL GPOPEN(LUTEMP,FILNODT30,'UNKNOWN',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      READ(LUTEMP) (WORK(KT3AM+I-1), I=1,NT1AMX*NT1AMX*NT1AMX)
      CALL GPCLOSE(LUTEMP,'KEEP')
C
      IF (IOPT .EQ. 1 .OR. (NONHF .AND. NFIELD.GT.0)) THEN
         LUSIFC = -1
         CALL GPOPEN(LUSIFC,'T3AMPL','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND LUSIFC
         WRITE (LUSIFC) (WORK(KT3AM-1+I), I=1,NT1AMX*NT1AMX*NT1AMX)
         CALL GPCLOSE(LUSIFC,'KEEP')
      ENDIF
C
C------------------------------------
C     Calculate the CC3 corrections.
C------------------------------------
C
      CALL DZERO(WORK(KOME1),NT1AMX)
C
      CALL T3_LEFT1(WORK(KOME1),WORK(KYIAJB),WORK(KIAJB),
     *              C2AM,WORK(KT3AM))
C
      CALL T3_LEFT2(WORK(KOME1),WORK(KYIAJB),WORK(KIAJB),
     *              C2AM,WORK(KT3AM))
C
      CALL T3_LEFT3(WORK(KOME1),WORK(KIAJB),C2AM,WORK(KT3AM))
C
      DO I = 1,NT1AMX
         OMEGA1(I) = OMEGA1(I) + WORK(KOME1+I-1)
      ENDDO
C
C---------------------------------------------------------------
C     Calculate rhs vector for triples left amplitude equations:
C---------------------------------------------------------------
C
      CALL DZERO(WORK(KL3AM),NT1AMX*NT1AMX*NT1AMX)

      CALL CCSDT_L3AM_R(WORK(KL3AM),FREQ,WORK(KINT1T),WORK(KINT2T),
     *                  WORK(KIAJB),WORK(KFOCK),C1AM,C2AM,
     *                  WORK(KSCR1),WORK(KFOCKD),.FALSE.)
C
C----------------------------------------------------------
C     Solve triples equations and partition solution vector
C     into the singles and doubles results vectors:
C----------------------------------------------------------
C
      CALL CC_LHPART_NODDY(OMEGA1,OMEGA2,WORK(KL3AM),FREQ,
     &                     WORK(KFOCKD),WORK(KONEP),
     &                     WORK(KIOOOO),WORK(KIOVVO),
     &                     WORK(KIOOVV),WORK(KIVVVV),
     &                     T2AM,WORK(KINT1S),WORK(KINT2S),
     &                     WORK(KEND1A),LWRK1A)


      IF (IOPT .EQ. 1) THEN
         LUSIFC = -1
         CALL GPOPEN(LUSIFC,'L3AMPL','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND LUSIFC
         WRITE (LUSIFC) (WORK(KL3AM-1+I), I=1,NT1AMX*NT1AMX*NT1AMX)
         CALL GPCLOSE(LUSIFC,'KEEP')
      ENDIF

      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'Result vector after CC_LHTR_NODDY:'
        CALL CC_PRP(OMEGA1,OMEGA2,1,1,1)
      END IF
 
      RETURN
      END
*---------------------------------------------------------------------*
C  /* Deck t3_left1 */
*---------------------------------------------------------------------*
      SUBROUTINE T3_LEFT1(OMEGA1,YIAJB,XIAJB,C2AM,T3AM)
*=====================================================================*
C
C
C
*---------------------------------------------------------------------*
#include <implicit.h>
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
      DIMENSION OMEGA1(NT1AMX), YIAJB(NT1AMX,NT1AMX)
      DIMENSION XIAJB(NT1AMX,NT1AMX)
      DIMENSION C2AM(NT1AMX,NT1AMX), T3AM(NT1AMX,NT1AMX,NT1AMX)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO A = 1, NVIRT
      DO I = 1, NRHFT
        NAI = NVIRT*(I-1) + A
        DO D = 1, NVIRT
        DO L = 1, NRHFT
          NDI = NVIRT*(I-1) + D
          NDL = NVIRT*(L-1) + D
          NAL = NVIRT*(L-1) + A
          DO E = 1, NVIRT
          DO M = 1, NRHFT
            NEM = NVIRT*(M-1) + E
            DO F = 1, NVIRT
            DO N = 1, NRHFT
              NEN = NVIRT*(N-1) + E
              NFN = NVIRT*(N-1) + F
              NFM = NVIRT*(M-1) + F
C
              OMEGA1(NAI) = OMEGA1(NAI)
     *                    + T3AM(NDL,NEM,NFN)
     *                    * C2AM(NFM,NDI)*YIAJB(NEN,NAL)
     *                    - T3AM(NDL,NEN,NFM)
     *                    * C2AM(NFM,NDI)*XIAJB(NEN,NAL)
            ENDDO
            ENDDO
          ENDDO
          ENDDO
        ENDDO
        ENDDO
      ENDDO
      ENDDO
C
      RETURN
      END
*---------------------------------------------------------------------*
C  /* Deck t3_left2 */
*---------------------------------------------------------------------*
      SUBROUTINE T3_LEFT2(OMEGA1,YIAJB,XIAJB,C2AM,T3AM)
*=====================================================================*
C
C
C
*---------------------------------------------------------------------*
#include <implicit.h>
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
      DIMENSION OMEGA1(NT1AMX), YIAJB(NT1AMX,NT1AMX)
      DIMENSION XIAJB(NT1AMX,NT1AMX)
      DIMENSION C2AM(NT1AMX,NT1AMX), T3AM(NT1AMX,NT1AMX,NT1AMX)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO A = 1, NVIRT
      DO I = 1, NRHFT
        NAI = NVIRT*(I-1) + A
        DO D = 1, NVIRT
        DO L = 1, NRHFT
          NDL = NVIRT*(L-1) + D
          DO E = 1, NVIRT
          DO M = 1, NRHFT
            NAM = NVIRT*(M-1) + A
            NEI = NVIRT*(I-1) + E
            NEM = NVIRT*(M-1) + E
            DO F = 1, NVIRT
            DO N = 1, NRHFT
              NDN = NVIRT*(N-1) + D
              NFL = NVIRT*(L-1) + F
              NFN = NVIRT*(N-1) + F
C
              OMEGA1(NAI) = OMEGA1(NAI)
     *                    + T3AM(NDL,NEM,NFN)
     *                    * C2AM(NAM,NDN)*YIAJB(NEI,NFL)
     *                    - T3AM(NDN,NEM,NFL)
     *                    * C2AM(NAM,NDN)*XIAJB(NEI,NFL)
            ENDDO
            ENDDO
          ENDDO
          ENDDO
        ENDDO
        ENDDO
      ENDDO
      ENDDO
C
      RETURN
      END
*---------------------------------------------------------------------*
C  /* Deck t3_left3 */
*---------------------------------------------------------------------*
      SUBROUTINE T3_LEFT3(OMEGA1,XIAJB,C2AM,T3AM)
*=====================================================================*
C
C     Note : XIAJB is L and not g.
C
C
*---------------------------------------------------------------------*
#include <implicit.h>
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
      DIMENSION OMEGA1(NT1AMX), XIAJB(NT1AMX,NT1AMX)
      DIMENSION C2AM(NT1AMX,NT1AMX), T3AM(NT1AMX,NT1AMX,NT1AMX)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO A = 1, NVIRT
      DO I = 1, NRHFT
        NAI = NVIRT*(I-1) + A
        DO D = 1, NVIRT
        DO L = 1, NRHFT
          NDL = NVIRT*(L-1) + D
          DO E = 1, NVIRT
          DO M = 1, NRHFT
            NEM = NVIRT*(M-1) + E
            DO F = 1, NVIRT
            DO N = 1, NRHFT
              NEN = NVIRT*(N-1) + E
              NFM = NVIRT*(M-1) + F
              NFN = NVIRT*(N-1) + F
C
              OMEGA1(NAI) = OMEGA1(NAI)
     *                    + (T3AM(NDL,NEM,NFN)-T3AM(NDL,NEN,NFM))
     *                    * C2AM(NEM,NDL)*XIAJB(NFN,NAI)
            ENDDO
            ENDDO
          ENDDO
          ENDDO
        ENDDO
        ENDDO
      ENDDO
      ENDDO
C
      RETURN
      END
*---------------------------------------------------------------------*
C  /* Deck ccsdt_l03am */
*=====================================================================*
      SUBROUTINE CCSDT_L03AM(L3AM,XINT1T,XINT2T,XIAJB,FOCK,
     *                       L1AM,L2AM,SCR1,FOCKD,FCKFLD,L3AM2)
*---------------------------------------------------------------------*
*
*     Purpose: compute zero-order triples lagrangian multipliers
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
      REAL XINT1T(*), XINT2T(*), XIAJB(*), FOCK(*)
      REAL L3AM(*), L2AM(*), L1AM(*), L3AM2(*)
      REAL SCR1(*), FOCKD(*), FCKFLD(*)
#else
      DOUBLE PRECISION XINT1T(*), XINT2T(*), XIAJB(*), FOCK(*)
      DOUBLE PRECISION L3AM(*), L2AM(*), L1AM(*), L3AM2(*)
      DOUBLE PRECISION SCR1(*), FOCKD(*), FCKFLD(*)
#endif
 
      CALL DZERO(L3AM,NT1AMX*NT1AMX*NT1AMX)

      CALL CCSDT_L3AM_R(L3AM,0.0D0,XINT1T,XINT2T,XIAJB,FOCK,
     *                  L1AM,L2AM,SCR1,FOCKD,.FALSE.)

      CALL CCSDT_3AM(L3AM,0.0D0,SCR1,FOCKD,NONHF,FCKFLD,.TRUE.,L3AM2)

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_L03AM                          *
*---------------------------------------------------------------------*
*---------------------------------------------------------------------*
C  /* Deck ccsdt_l3am_r */
*=====================================================================*
      SUBROUTINE CCSDT_L3AM_R(T3BAR,FREQ,XINT1,XINT2,XIAJB,FOCK,T1AM,
     *                        T2AM,SCR1,FOCKD,DIVIDE)
*---------------------------------------------------------------------*
*
*     Purpose: compute contributions to left triples amplitudes
*
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit.h>
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccfield.h>
C
#if defined (SYS_CRAY) 
      REAL XINT1(NT1AMX,NVIRT,NVIRT)
      REAL XINT2(NT1AMX,NRHFT,NRHFT)
      REAL XIAJB(NT1AMX,NT1AMX), FOCK(NORBT,NORBT)
      REAL T3BAR(NT1AMX,NT1AMX,NT1AMX), SCR1(NT1AMX)
      REAL T1AM(NT1AMX), T2AM(NT1AMX,NT1AMX), FOCKD(*)
      REAL AIBJCK, FREQ, TWO
#else
      DOUBLE PRECISION XINT1(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT2(NT1AMX,NRHFT,NRHFT) 
      DOUBLE PRECISION XIAJB(NT1AMX,NT1AMX), FOCK(NORBT,NORBT)
      DOUBLE PRECISION T3BAR(NT1AMX,NT1AMX,NT1AMX), SCR1(NT1AMX)
      DOUBLE PRECISION T1AM(NT1AMX), T2AM(NT1AMX,NT1AMX), FOCKD(*)
      DOUBLE PRECISION AIBJCK, FREQ, TWO
#endif
 
      PARAMETER (TWO = 2.0D0)

      LOGICAL DIVIDE
 
      INTEGER NI, NA, NAI, NK, NC, NCK, NJ, NCJ, NBJ, NBK, NBI, NCI,
     &        NAJ, NAK, ND, NDI, NDJ, NDK, NL, NAL, NBL, NB
C
      ! we can just divide by orbital energy denominators if we
      ! have non-HF external fields
      IF (DIVIDE .AND. NONHF) 
     &    CALL QUIT('NONHF Problem in CCSDT_L3AM_R!')
C
      DO NI = 1,NRHFT
         DO NA = 1,NVIRT
            NAI = NVIRT*(NI-1) + NA
            SCR1(NAI) = FOCKD(NRHFT+NA) - FOCKD(NI)
         END DO
      END DO
C
      DO 100 NK = 1,NRHFT
      DO 105 NC = 1,NVIRT
C
         NCK = NVIRT*(NK-1) + NC
C
         DO 110 NJ = 1,NRHFT
            NCJ = NVIRT*(NJ-1) + NC
            DO 120 NB = 1,NVIRT
C
               NBJ = NVIRT*(NJ-1) + NB
               NBK = NVIRT*(NK-1) + NB
C
               DO 130 NI = 1,NRHFT
                  NBI = NVIRT*(NI-1) + NB
                  NCI = NVIRT*(NI-1) + NC
               DO 135 NA = 1,NVIRT
C
                  NAI = NVIRT*(NI-1) + NA
                  NAJ = NVIRT*(NJ-1) + NA
                  NAK = NVIRT*(NK-1) + NA
C
                  AIBJCK = 0.0D0
C
C                 T1* = TWO T1 => Factor of two
                  AIBJCK = AIBJCK + XIAJB(NBJ,NCK)*T1AM(NAI)
C
C                 T1* = TWO T1 => Factor of two
                  AIBJCK = AIBJCK - XIAJB(NBJ,NCI)*T1AM(NAK)
C
                  AIBJCK = AIBJCK + FOCK(NK,NRHFT+NC)*
     *                           T2AM(NAI,NBJ)
C
                  AIBJCK = AIBJCK - FOCK(NJ,NRHFT+NC)*
     *                           T2AM(NAI,NBK)
C
                  DO 140 ND = 1,NVIRT
C
                     NDI = NVIRT*(NI-1) + ND
                     NDJ = NVIRT*(NJ-1) + ND
                     NDK = NVIRT*(NK-1) + ND
C
                     AIBJCK = AIBJCK + 
     *                        (TWO*XINT1(NCK,ND,NB)-XINT1(NBK,ND,NC))*
     *                        T2AM(NDJ,NAI)
C
                     AIBJCK = AIBJCK - 
     *                        XINT1(NBI,ND,NC)*
     *                        T2AM(NDK,NAJ)
C
  140             CONTINUE
C
                  DO 150 NL = 1,NRHFT
C
                     NAL = NVIRT*(NL-1) + NA
                     NBL = NVIRT*(NL-1) + NB
C
                     AIBJCK = AIBJCK - 
     *                        (TWO*XINT2(NCK,NJ,NL)-XINT2(NCJ,NK,NL))*
     *                        T2AM(NBL,NAI)
C
                     AIBJCK = AIBJCK + 
     *                        XINT2(NCJ,NI,NL)*
     *                        T2AM(NBK,NAL)
C
  150             CONTINUE
C
                  IF (DIVIDE) THEN
                    AIBJCK = AIBJCK/(SCR1(NAI)+SCR1(NBJ)+SCR1(NCK)-FREQ)
                  END IF
C
                  T3BAR(NAI,NBJ,NCK) = T3BAR(NAI,NBJ,NCK) + AIBJCK
                  T3BAR(NAI,NCK,NBJ) = T3BAR(NAI,NCK,NBJ) + AIBJCK
                  T3BAR(NBJ,NAI,NCK) = T3BAR(NBJ,NAI,NCK) + AIBJCK
                  T3BAR(NCK,NAI,NBJ) = T3BAR(NCK,NAI,NBJ) + AIBJCK
                  T3BAR(NBJ,NCK,NAI) = T3BAR(NBJ,NCK,NAI) + AIBJCK
                  T3BAR(NCK,NBJ,NAI) = T3BAR(NCK,NBJ,NAI) + AIBJCK
C
  135          CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  105 CONTINUE
  100 CONTINUE
C
C------------------------------------------------------
C     Get rid of amplitudes that are not allowed.
C------------------------------------------------------
C
      CALL CCSDT_CLEAN_T3(T3BAR,NT1AMX,NVIRT,NRHFT) 
C
      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_L3AM_R                         *
*---------------------------------------------------------------------*
C  /* Deck cc_lhpart_noddy */
*=====================================================================*
      SUBROUTINE CC_LHPART_NODDY(OMEGA1,OMEGA2,L3AM,FREQ,
     &                           FOCKD,FIELD,
     &                           XOOOO,XOVVO,XOOVV,XVVVV,
     &                           T2AM,XINT1S,XINT2S,
     &                           WORK,LWORK)
*---------------------------------------------------------------------*
*
*     Purpose: solve 'left' triples equations and partition the
*              triples solution vector into an effective lhs vector
*
*              Note: in case of non-HF external fields zero-order
*              cluster amplitudes must be available on file 'T3AMPL'
*
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <dummy.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccfield.h>
#include <maxorb.h>
#include <ccorb.h>

#if defined (SYS_CRAY)
      REAL OMEGA1(*), OMEGA2(*), L3AM(*), FREQ, FOCKD(*)
      REAL FIELD(*), T2AM(*), XINT1S(*), XINT2S(*), WORK(*)
      REAL XOOOO(*), XOVVO(*), XOOVV(*), XVVVV(*)
#else
      DOUBLE PRECISION OMEGA1(*), OMEGA2(*), L3AM(*), FREQ, FOCKD(*)
      DOUBLE PRECISION FIELD(*), T2AM(*), XINT1S(*), XINT2S(*), WORK(*)
      DOUBLE PRECISION XOOOO(*), XOVVO(*), XOOVV(*), XVVVV(*)
#endif

      LOGICAL TRANSPOSE
      INTEGER LWORK, KSCR1, KEND1, KL3SCR, LWRK1, KOME1, KOME2, KT03AM,
     &        LUTEMP, IJ, NIJ, INDEX
 
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J

      CALL QENTER('CCLPARTN')

*---------------------------------------------------------------------*
*     Solve triples equations:
*---------------------------------------------------------------------*
      KSCR1 = 1
      KEND1 = KSCR1 + NT1AMX
      IF (NONHF) THEN
         KL3SCR = KEND1
         KEND1  = KL3SCR + NT1AMX*NT1AMX*NT1AMX
      END IF
 
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CC_LHPART_NODDY (1)')
      ENDIF

      TRANSPOSE = .TRUE.
      CALL CCSDT_3AM(L3AM,FREQ,WORK(KSCR1),FOCKD,
     *               NONHF,FIELD,TRANSPOSE,WORK(KL3SCR))

      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0d0,L3AM,1)
 
*---------------------------------------------------------------------*
*     Add contribution to singles doubles result vectors:
*---------------------------------------------------------------------*
      KOME1 = 1
      KOME2 = KOME1 + NT1AMX
      KEND1 = KOME2 + NT1AMX*NT1AMX

      IF (NONHF) THEN
         KT03AM = KEND1
         KEND1  = KT03AM + NT1AMX*NT1AMX*NT1AMX
      END IF
 
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CC_LHPART_NODDY (2)')
      ENDIF

      CALL DZERO(WORK(KOME1),NT1AMX)
      CALL DZERO(WORK(KOME2),NT1AMX*NT1AMX)
 
      CALL CC3_L3_OMEGA1_NODDY(WORK(KOME1),L3AM,
     *                         XOOOO,XOVVO,XOOVV,XVVVV,T2AM)
 
      CALL CC3_L3_OMEGA2_NODDY(WORK(KOME2),L3AM,XINT1S,XINT2S)

c     ------------------------------------------
c     contributions from non-HF external fields:
c     ------------------------------------------
      IF (NONHF) THEN
         LUTEMP = -1
         CALL GPOPEN(LUTEMP,'T3AMPL','UNKNOWN',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         REWIND LUTEMP
         READ (LUTEMP) (WORK(KT03AM-1+I),I=1,NT1AMX*NT1AMX*NT1AMX)
         CALL GPCLOSE(LUTEMP,'KEEP')

         CALL CCSDT_E1AM(WORK(KOME1),L3AM,WORK(KT03AM),FIELD)

         CALL CCSDT_E2AM(WORK(KOME2),L3AM,T2AM,FIELD)
      ENDIF
 

      DO I = 1,NT1AMX
         OMEGA1(I) = OMEGA1(I) + WORK(KOME1+I-1)
      END DO  
 
      DO I = 1,NT1AMX
         DO J = 1,I
            IJ = NT1AMX*(I-1) + J
            NIJ = INDEX(I,J)
            OMEGA2(NIJ) = OMEGA2(NIJ) + WORK(KOME2+IJ-1)
         ENDDO
      ENDDO

      CALL QEXIT('CCLPARTN')
      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_LHPART_NODDY                      *
*---------------------------------------------------------------------*
C  /* Deck cc3__l3_omega */
      SUBROUTINE CC3_L3_OMEGA2_NODDY(OMEGA2,L3AM,XINT1,XINT2)
C
C     Kasper Hald, Spring 2002.
C
C     Calculate the doubles of the L3 part of CC3 left hand side doubles
C
#include <implicit.h>
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INTEGER INDEX, NAI, NBJ, NBK, NAK, NCK, NCJ, NCI, NDJ, NDK, NDI
      INTEGER NBL, NCL, NAL, NAJ, NBI, AI
C
#if defined (SYS_CRAY)
      REAL OMEGA2(NT1AMX,NT1AMX)
      REAL L3AM(NT1AMX,NT1AMX,NT1AMX)
      REAL XINT1(NT1AMX,NVIRT,NVIRT)
      REAL XINT2(NT1AMX,NRHFT,NRHFT)
      REAL ONE, TWO, XAIBJ
#else
      DOUBLE PRECISION OMEGA2(NT1AMX,NT1AMX)
      DOUBLE PRECISION L3AM(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION XINT1(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT2(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION ONE, TWO, XAIBJ
#endif
      PARAMETER (ONE = 1.0D0, TWO = 2.0D0)
C
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
C===================================================
C     Calculate the L3 contributions to omega2
C===================================================
C
      CALL DZERO(OMEGA2,NT1AMX*NT1AMX)
C
      DO NI = 1,NRHFT
         DO NA = 1,NVIRT
            NAI = NVIRT*(NI-1) + NA
C
            DO NJ = 1,NRHFT
               DO NB = 1,NVIRT
                  NBJ = NVIRT*(NJ-1) + NB
C
                  DO NK = 1,NRHFT
                     NBK = NVIRT*(NK-1) + NB
                     NAK = NVIRT*(NK-1) + NA
                     DO NC = 1,NVIRT
C
                        NCK = NVIRT*(NK-1) + NC
                        NCJ = NVIRT*(NJ-1) + NC
                        NCI = NVIRT*(NI-1) + NC
C
                        DO ND = 1,NVIRT
C
                           NDJ = NVIRT*(NJ-1) + ND
                           NDK = NVIRT*(NK-1) + ND
                           NDI = NVIRT*(NI-1) + ND
C
                           OMEGA2(NBJ,NAI) = OMEGA2(NBJ,NAI) +
     *              L3AM(NBJ,NCI,NDK)*XINT1(NDK,NC,NA)
C
C
                        ENDDO
C
                        DO NL = 1,NRHFT
C
                           NBL = NVIRT*(NL-1) + NB
                           NCL = NVIRT*(NL-1) + NC
                           NAL = NVIRT*(NL-1) + NA
C
                           OMEGA2(NBJ,NAI) = OMEGA2(NBJ,NAI)
     *                     - L3AM(NBJ,NAK,NCL)*XINT2(NCL,NI,NK)
C
                        ENDDO
C
                     ENDDO
                  ENDDO
C
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C
C----------------------------------------------------
C     Make P^{ab}_{ij} for the T3BAR contributions.
C----------------------------------------------------
C
      DO NAI = 1,NT1AMX
         DO NBJ = 1,NAI
C
            XAIBJ = OMEGA2(NAI,NBJ) + OMEGA2(NBJ,NAI)
            OMEGA2(NAI,NBJ) = XAIBJ
            OMEGA2(NBJ,NAI) = XAIBJ
C
         ENDDO
      ENDDO
C
      RETURN
      END
C  /* Deck cc3__l3_omega1_noddy */
      SUBROUTINE CC3_L3_OMEGA1_NODDY(OMEGA1,L3AM,
     *                                XOOOO,XOVVO,XOOVV,XVVVV,T2AM)
C
C     Kasper Hald, Spring 2002.
C
C     Calculate the doubles of the L3 part of CC3 left hand side singles
C
#include <implicit.h>
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
C
      INTEGER INDEX, NAI, NBJ, NBK, NAK, NCK, NCJ, NCI, NDJ, NDK, NDI
      INTEGER NBL, NCL, NAL, NAJ, NBI, AI
C
#if defined (SYS_CRAY)
      REAL OMEGA1(NT1AMX)
      REAL L3AM(NT1AMX,NT1AMX,NT1AMX)
      REAL XOOOO(NRHFT,NRHFT,NRHFT,NRHFT)
      REAL XOVVO(NRHFT,NVIRT,NVIRT,NRHFT)
      REAL XOOVV(NRHFT,NRHFT,NVIRT,NVIRT)
      REAL XVVVV(NVIRT,NVIRT,NVIRT,NVIRT)
      REAL T2AM(NT1AMX,NT1AMX)
      REAL ONE, TWO, XAIBJ
#else
      DOUBLE PRECISION OMEGA1(NT1AMX)
      DOUBLE PRECISION L3AM(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION XOOOO(NRHFT,NRHFT,NRHFT,NRHFT)
      DOUBLE PRECISION XOVVO(NRHFT,NVIRT,NVIRT,NRHFT)
      DOUBLE PRECISION XOOVV(NRHFT,NRHFT,NVIRT,NVIRT)
      DOUBLE PRECISION XVVVV(NVIRT,NVIRT,NVIRT,NVIRT)
      DOUBLE PRECISION T2AM(NT1AMX,NT1AMX)
      DOUBLE PRECISION ONE, TWO, XAIBJ
#endif
      PARAMETER (ONE = 1.0D0, TWO = 2.0D0)
C
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
C===================================================
C     Calculate the L3BAR contributions to omega1
C===================================================
C
C-------------------------------------------
C     4 occupied integrals in integrals
C-------------------------------------------
C
      if (.true.) then
      DO NA = 1, NVIRT
       DO NB = 1, NVIRT
        DO NC = 1, NVIRT
         DO NI = 1, NRHFT
          NAI = NVIRT*(NI-1) + NA
          DO NJ = 1, NRHFT
           NBJ = NVIRT*(NJ-1) + NB
           DO NK = 1, NRHFT
            NCK = NVIRT*(NK-1) + NC
            DO NL = 1, NRHFT
             NAL = NVIRT*(NL-1) + NA
             DO NM = 1, NRHFT
                NCM = NVIRT*(NM-1) + NC
                OMEGA1(NAI) = OMEGA1(NAI)
     *           + L3AM(NBJ,NAL,NCM)*T2AM(NBJ,NCK)*XOOOO(NI,NL,NK,NM)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      endif
C
C-------------------------------------------
C     4 virtual integrals in integrals
C-------------------------------------------
C
      if (.true.) then
      DO NA = 1, NVIRT
       DO NB = 1, NVIRT
        DO NC = 1, NVIRT
         DO ND = 1, NVIRT
          DO NE = 1, NVIRT
           DO NI = 1, NRHFT
            NAI = NVIRT*(NI-1) + NA
            NDI = NVIRT*(NI-1) + ND
            DO NJ = 1, NRHFT
             NBJ = NVIRT*(NJ-1) + NB
             DO NK = 1, NRHFT
              NEK = NVIRT*(NK-1) + NE
              NCK = NVIRT*(NK-1) + NC
C
              OMEGA1(NAI) = OMEGA1(NAI)
     *          + L3AM(NBJ,NDI,NEK)*T2AM(NBJ,NCK)*XVVVV(ND,NA,NE,NC)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      endif
C
C--------------------------------------------------------
C     2 terms with g_{ovvo} and 2 terms with g_{oovv}
C--------------------------------------------------------
C
      if (.true.) then
      DO NA = 1, NVIRT
       DO NB = 1, NVIRT
        DO NC = 1, NVIRT
         DO ND = 1, NVIRT
          DO NI = 1, NRHFT
           NAI = NVIRT*(NI-1) + NA
           NCI = NVIRT*(NI-1) + NC
           NDI = NVIRT*(NI-1) + ND
           DO NJ = 1, NRHFT
            NBJ = NVIRT*(NJ-1) + NB
            DO NK = 1, NRHFT
             NAK = NVIRT*(NK-1) + NA
             NCK = NVIRT*(NK-1) + NC
             NDK = NVIRT*(NK-1) + ND
             DO NL = 1, NRHFT
              NAL = NVIRT*(NL-1) + NA
              NCL = NVIRT*(NL-1) + NC
              NDL = NVIRT*(NL-1) + ND
C
C - L3^{daf}_{lmn} * t^{de}_{lm} * g_{iefn}
              OMEGA1(NAI) = OMEGA1(NAI)
     *       - L3AM(NBJ,NAK,NDL)*T2AM(NBJ,NCK)*XOVVO(NI,NC,ND,NL)
C
C - L3^{daf}_{lnm} * t^{de}_{lm} * g_{infe}
              OMEGA1(NAI) = OMEGA1(NAI)
     *       - L3AM(NBJ,NAL,NDK)*T2AM(NBJ,NCK)*XOOVV(NI,NL,ND,NC)
C
C - L3^{def}_{lin} * t^{de}_{lm} * g_{mafn}
              OMEGA1(NAI) = OMEGA1(NAI)
     *       - L3AM(NBJ,NCI,NDL)*T2AM(NBJ,NCK)*XOVVO(NK,NA,ND,NL)
C
C - L3^{def}_{lni} * t^{de}_{lm} * g_{mnfa}
              OMEGA1(NAI) = OMEGA1(NAI)
     *       - L3AM(NBJ,NCL,NDI)*T2AM(NBJ,NCK)*XOOVV(NK,NL,ND,NA)
C
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      endif
C
C=====================================
C     Write out all results
C=====================================
C
      if (locdbg) then
       write(lupri,*) '        '
       do nai = 1, nt1amx
          if (abs(omega1(nai)) .gt. 1.0D-12) then
             write(lupri,*) 'Omega1_noddy(',nai,')             = ',
     *                       omega1(nai)
          endif
       enddo
       write(lupri,*) '        '
      end if
C
      RETURN
      END
C  /* Deck ccfop_tran1 */
      SUBROUTINE CCFOP_TRAN1(XINT1,XINT2,XINT3,XINT4,
     *                       XLAMDP,XLAMDH,AOINT,IDEL)
C
C
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
      DIMENSION XLAMDP(NBAST,NORBT), XLAMDH(NBAST,NORBT)
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
     *                  + AOINT(NAB,G)*XLAMDP(A,NK)*XLAMDH(IB,NRHFT+NC)
     *                         *XLAMDP(G,NRHFT+ND)*XLAMDH(IDEL,NL)
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
     *                  + AOINT(NAB,G)*XLAMDP(A,NK)*XLAMDH(IB,NL)
     *                         *XLAMDP(G,NRHFT+NC)*XLAMDH(IDEL,NRHFT+ND)
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
     *                  + AOINT(NAB,G)*XLAMDP(A,NK)*XLAMDH(IB,NL)
     *                         *XLAMDP(G,NM)*XLAMDH(IDEL,NN)
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
     *           + AOINT(NAB,G)*XLAMDP(A,NRHFT+NC)*XLAMDH(IB,NRHFT+ND)
     *                         *XLAMDP(G,NRHFT+NE)*XLAMDH(IDEL,NRHFT+NF)
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
C  /* Deck cc_lhtr_noddy2 */
      SUBROUTINE CC_LHTR_NODDY2(OMEGA1,OMEGA2,OMEGA3,T1AM,T2AM,T3AM,
     *                          C1AM,C2AM,C3AM,XLAMDP,XLAMDH,FOCK,
     *                          WORK,LWORK)
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
C     IOPT = 0 No storage of T3, L3
C     IOPT = 1 Store T3 and L3 in 'T3AMPL' and 'L3AMPL'
C
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <maxorb.h>
#include <maxash.h>
#include <aovec.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION INDEXA(MXCORB)
      DIMENSION OMEGA1(*),OMEGA2(*), OMEGA3(*)
      DIMENSION T1AM(*),T2AM(*), T3AM(*)
      DIMENSION C1AM(*), C2AM(*), C3AM(*)
      DIMENSION XLAMDP(*),XLAMDH(*), FOCK(*), WORK(LWORK)
      LOGICAL   L2INCL
C
#include <ccfield.h>
#include <inforb.h>
CCN #include <infind.h> not compatible with R12!
#include <ccisao.h>
#include <r12int.h>
#include <blocks.h>
#include <inftap.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccsdio.h>
COMMENT COMMENT
COMMENT COMMENT
#include <mxcent.h>
#include <eribuf.h>
COMMENT COMMENT
COMMENT COMMENT
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      IF (NSYM .NE. 1) THEN
         CALL QUIT('CC_LHTR_NODDY2: No symmetry in this part yet')
      ENDIF
C
C------------------------
C     Dynamic Allocation.
C------------------------
C
      KCMO   = 1
      KFOCKD = KCMO   + NLAMDT
      KINT1  = KFOCKD + NORBT
      KINT2  = KINT1  + NT1AMX*NVIRT*NVIRT
      KINT1T = KINT2  + NRHFT*NRHFT*NT1AMX
      KINT2T = KINT1T + NT1AMX*NVIRT*NVIRT
      KINT1S = KINT2T + NRHFT*NRHFT*NT1AMX
      KINT2S = KINT1S + NT1AMX*NVIRT*NVIRT
      KINT3S = KINT2S + NRHFT*NRHFT*NT1AMX
      KINT4S = KINT3S + NT1AMX*NVIRT*NVIRT
      KOME1  = KINT4S + NRHFT*NRHFT*NT1AMX
      KOME2  = KOME1  + NT1AMX
      KOME3  = KOME2  + NT1AMX*NT1AMX
      KIAJB  = KOME3  + NT1AMX*NT1AMX*NT1AMX
      KYIAJB = KIAJB  + NT1AMX*NT1AMX
      KEND1  = KYIAJB + NT1AMX*NT1AMX
      LWRK1  = LWORK  - KEND1
C
C     New Integrals.
      KIOVVO = KEND1
      KIOOVV = KIOVVO + NRHFT*NVIRT*NVIRT*NRHFT
      KIOOOO = KIOOVV + NRHFT*NRHFT*NVIRT*NVIRT
      KIVVVV = KIOOOO + NRHFT*NRHFT*NRHFT*NRHFT
      KEND1  = KIVVVV + NVIRT*NVIRT*NVIRT*NVIRT
      LWRK1  = LWORK  - KEND1
C
      IF ((NONHF) .AND. (NFIELD .GT. 0)) THEN
         KONEP  = KEND1
         KT3SCR = KONEP  + NORBT*NORBT
         KEND1  = KT3SCR + NT1AMX*NT1AMX*NT1AMX
         LWRK1  = LWORK  - KEND1
      ENDIF
C
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
      CALL DZERO(WORK(KINT3S),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT4S),NT1AMX*NRHFT*NRHFT)
      CALL DZERO(WORK(KOME1),NT1AMX)
      CALL DZERO(WORK(KOME2),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KOME3),NT1AMX*NT1AMX*NT1AMX)
      CALL DZERO(WORK(KIAJB),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KYIAJB),NT1AMX*NT1AMX)
C
      CALL DZERO(WORK(KIOVVO),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KIOOVV),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KIOOOO),NRHFT*NRHFT*NRHFT*NRHFT)
      CALL DZERO(WORK(KIVVVV),NVIRT*NVIRT*NVIRT*NVIRT)
C
C=======================================================
C     Get the one electron integrals from the fields
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
      IF (DIRECT) THEN
         CALL QUIT('Direct not implemented')
#ifdef NOT_IMPLEMENTED
         NTOSYM = 1
         DTIME  = SECOND()
         CALL HERDI1(WORK(KEND1),LWRK1,IPRINT)
         DTIME  = SECOND() - DTIME
         TIMHER1 = TIMHER1 + DTIME
#endif
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
               CALL QUIT('Direct not implemented')
#ifdef NOT_IMPLEMENTED
               DTIME  = SECOND()
               CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,IPRINT)
               DTIME  = SECOND() - DTIME
               TIMHER2 = TIMHER2 + DTIME
COMMENT COMMENT
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient core in CC_FOP3')
               END IF
COMMENT COMMENT
#endif
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
#ifdef NOT_IMPLEMENTED
                  IDEL  = INDEXA(IDEL2)
                  IF (NOAUXB) THEN
                     IDUM = 1
                     CALL IJKAUX(IDEL,IDUM,IDUM,IDUM)
                  END IF
                  ISYMD = ISAO(IDEL)
#endif
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
c              DTIME   = SECOND()
      call quit('CC_LHTR_NODDY2: KRECNR not defined, FIXME')
      krecnr = 1 ! hjaaj Aug 07: to silence ftnchek
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                     WORK(KRECNR),DIRECT)
c              DTIME   = SECOND() - DTIME
c              TIMRDAO = TIMRDAO  + DTIME
C
C----------------------------------------------------
C              Calculate integrals needed in CCSD[T].
C----------------------------------------------------
C
               CALL CCSDT_TRAN1(WORK(KINT1T),WORK(KINT2T),XLAMDP,
     *                          XLAMDH,WORK(KXINT),IDEL)
C
               CALL CC3_TRAN2(WORK(KIAJB),WORK(KYIAJB),XLAMDP,
     *                          XLAMDH,WORK(KXINT),IDEL)
C
               CALL CCSDT_TRAN3(WORK(KINT1S),WORK(KINT2S),XLAMDP,
     *                          XLAMDH,WORK(KXINT),IDEL)
C
               CALL CCFOP_TRAN1(WORK(KIOVVO),WORK(KIOOVV),
     *                          WORK(KIOOOO),WORK(KIVVVV),
     *                          XLAMDP,XLAMDH,WORK(KXINT),IDEL)
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C------------------------------------
C     Calculate the corrections from T3.
C------------------------------------
C
      CALL T3_LEFT1(WORK(KOME1),WORK(KYIAJB),WORK(KIAJB),
     *              C2AM,T3AM)
C
      CALL T3_LEFT2(WORK(KOME1),WORK(KYIAJB),WORK(KIAJB),
     *              C2AM,T3AM)
C
      CALL T3_LEFT3(WORK(KOME1),WORK(KIAJB),C2AM,T3AM)
C
      DO I = 1,NT1AMX
         OMEGA1(I) = OMEGA1(I) + WORK(KOME1+I-1)
      ENDDO
C
C----------------------------------------------------------
C         Calculate the contributions to Omega2 from C3
C         The result is written out inside the routine.
C----------------------------------------------------------
C
      CALL DZERO(WORK(KOME1),NT1AMX)
      CALL DZERO(WORK(KOME2),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KOME3),NT1AMX*NT1AMX*NT1AMX)
 
      CALL CC3_L3_OMEGA1_NODDY(WORK(KOME1),C3AM,
     *                         WORK(KIOOOO),WORK(KIOVVO),
     *                         WORK(KIOOVV),WORK(KIVVVV),T2AM)
 
      CALL CC3_L3_OMEGA2_NODDY(WORK(KOME2),C3AM,
     *                         WORK(KINT1S),WORK(KINT2S))

      IF (NONHF .AND. NFIELD.GT.0) THEN
         CALL CCSDT_E1AM(WORK(KOME1),C3AM,T3AM,WORK(KONEP))

         CALL CCSDT_E2AM(WORK(KOME2),C3AM,T2AM,WORK(KONEP))
      ENDIF
 
      DO 300 I = 1,NT1AMX
         OMEGA1(I) = OMEGA1(I) + WORK(KOME1+I-1)
  300 CONTINUE

      CALL DAXPY(NT1AMX*NT1AMX,1.0D0,WORK(KOME2),1,OMEGA2,1)
C
C---------------------------------------------------------
C     Calculate the contributions to omega3 from C1, C2
C---------------------------------------------------------
      DO I = 1, NRHFT
         WORK(KFOCKD-1+I) = 0.0d0
      ENDDO
      DO I = 1, NVIRT
         WORK(KFOCKD-1+NRHFT+I) = 1.0d0/3.0d0
      ENDDO
C
C     for non-HF external fields the above trick for the 
C     orbital energies does solve all problems...
C     --> change to CCSDT_L3AM_R
      IF (NONHF) CALL QUIT('fix me!')
C
C hjaaj-aug-2007: changed undefined KSCR1 to KEND1, hope it is OK.
      CALL CCSDT_L03AM(WORK(KOME3),WORK(KINT1T),WORK(KINT2T),
     *                WORK(KIAJB),FOCK,C1AM,C2AM,
     *                WORK(KEND1),WORK(KFOCKD),
     *                WORK(KONEP),WORK(KT3SCR))

      IF (NONHF .AND. NFIELD.GT.0) THEN
          L2INCL = .TRUE.
          CALL CCSDT_E3AM(WORK(KOME3),C2AM,C3AM,WORK(KONEP),L2INCL)
      END IF

      CALL DAXPY(NT1AMX*NT1AMX*NT1AMX,1.0D0,WORK(KOME3),1,OMEGA3,1)
C
C----------------------------------------
C     Print norm of result vectors.
C----------------------------------------
C
      XNORM = DDOT(NT1AMX,OMEGA1,1,OMEGA1,1)
      write(lupri,*) 'Norm Omega1 at the end of noddy :',xnorm
      XNORM = DDOT(NT1AMX**2,OMEGA2,1,OMEGA2,1)
      write(lupri,*) 'Norm Omega2 at the end of noddy :',xnorm
      XNORM = DDOT(NT1AMX**3,OMEGA3,1,OMEGA3,1)
      write(lupri,*) 'Norm Omega3 at the end of noddy :',xnorm
C
      RETURN
      END
