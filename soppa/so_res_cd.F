C  /* Deck so_res_cd */
      SUBROUTINE SO_RES_CD(RES2E,LRES2E,RES2D,LRES2D,TR2E,LTR2E,TR2D,
     &                     LTR2D,FOCKD,LFOCKD,ISYRES,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, April 1996
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Calculate contribution from D matrix to the 2p2h
C              result vectors.
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION RES2E(LRES2E), RES2D(LRES2D)
      DIMENSION TR2E(LTR2E),   TR2D(LTR2D)
      DIMENSION FOCKD(LFOCKD), WORK(LWORK)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <soppinf.h>
C
C------------------------------
C     Statement function INDEX.
C------------------------------
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_RES_CD')
C
C------------------------------------------------
C     Loop over the combined symmetry of B and J.
C------------------------------------------------
C
      DO 100 ISYMBJ = 1,NSYM
C
         ISYMAI = MULD2H(ISYMBJ,ISYRES)
C
C---------------------------------
C        Allocation of work space.
C---------------------------------
C
         LFBJ   = NT1AM(ISYMBJ)
         LFAI   = NT1AM(ISYMAI)
C
         KFBJ    = 1
         KFAI    = KFBJ  + LFBJ
         KEND    = KFAI  + LFAI
         LWORK1  = LWORK - KEND
C
         CALL SO_MEMMAX ('SO_RES_CD.1',LWORK1)
         IF (LWORK1 .LT. 0) CALL STOPIT('SO_RES_CD.1',' ',KEND,LWORK)
C
C----------------------------------------------------------------
C        Make difference of fock-diagonals B and J in WORK(KFBJ).
C----------------------------------------------------------------
C
         DO 201 ISYMJ = 1,NSYM
C
            ISYMB = MULD2H(ISYMJ,ISYMBJ)
C
            DO 202 J = 1,NRHF(ISYMJ)
C
               KOFFJ = IRHF(ISYMJ) + J
C
               DO 203 B = 1,NVIR(ISYMB)
C
                  NBJ   = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
                  KOFFB = IVIR(ISYMB) + B
C
                  WORK(KFBJ+NBJ) =  FOCKD(KOFFB) - FOCKD(KOFFJ)
C
  203          CONTINUE
C
  202       CONTINUE
C
  201    CONTINUE
C
C----------------------------------------------------------------
C        Make difference of fock-diagonals A and I in WORK(KFAI).
C----------------------------------------------------------------
C
         DO  301 ISYMI = 1,NSYM
C
            ISYMA = MULD2H(ISYMI,ISYMAI)
C
            DO 302 I = 1,NRHF(ISYMI)
C
               KOFFI = IRHF(ISYMI) + I
C
               DO 303 A = 1,NVIR(ISYMA)
C
                  NAI   = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
C
                  KOFFA = IVIR(ISYMA) + A
C
                  WORK(KFAI+NAI) =  FOCKD(KOFFA) - FOCKD(KOFFI)
C
  303          CONTINUE
C
  302       CONTINUE
C
  301    CONTINUE
C
C---------------------------------------------------------------
C        Multiply energy-differences EAIBJ and 2p2h trialvectors
C        to obtain the D-matrix contribution to the 2p2h result-
C        vectors.
C---------------------------------------------------------------
C
         IF ( ISYMAI .EQ. ISYMBJ) THEN
C
            DO 401 NBJ = 1,NT1AM(ISYMBJ)
C
               DO 402 NAI = 1,NBJ
C
                  NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     &                  + INDEX(NAI,NBJ)
C
                  EAIBJ = WORK(KFAI+NAI) + WORK(KFBJ+NBJ)
C
                  RES2E(NAIBJ) = RES2E(NAIBJ)
     &                         + EAIBJ * TR2E(NAIBJ)
                  RES2D(NAIBJ) = RES2D(NAIBJ)
     &                         + EAIBJ * TR2D(NAIBJ)
C
  402          CONTINUE
C
  401       CONTINUE
C
         ELSE IF ( ISYMAI .LT. ISYMBJ) THEN
C
            DO 501 NBJ = 1,NT1AM(ISYMBJ)
C
               DO 502 NAI = 1,NT1AM(ISYMAI)
C
                  NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     &                  + NT1AM(ISYMAI)*(NBJ - 1) + NAI
C
                  EAIBJ = WORK(KFAI+NAI) + WORK(KFBJ+NBJ)
C
                  RES2E(NAIBJ) = RES2E(NAIBJ)
     &                         + EAIBJ * TR2E(NAIBJ)
                  RES2D(NAIBJ) = RES2D(NAIBJ)
     &                         + EAIBJ * TR2D(NAIBJ)
C
  502          CONTINUE
C
  501       CONTINUE
C
         END IF
C
  100 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_RES_CD')
C
      RETURN
      END
