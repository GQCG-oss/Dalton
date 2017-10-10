C  /* Deck so_res_cdt */
      SUBROUTINE SO_RES_CDT(RES2E,LRES2E,RES2D,LRES2D,
     &                      FOCKD,LFOCKD,ISYRES,
     &                      DO_DEX,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, April 1996
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C     Pi A. B. Haase 29.04.2016: Triplet version
C     Rasmus Faber, 2017: Now done in-place
C
C     PURPOSE: Calculate contribution from D matrix to the 2p2h
C              result vectors.
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION RES2E(LRES2E), RES2D(LRES2D)
      DIMENSION FOCKD(LFOCKD), WORK(LWORK)
      LOGICAL, INTENT(IN) :: DO_DEX
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "soppinf.h"
CPi moved triplet variables to soppinf.h
C#include "infsop.h"
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_RES_CDT')
C
C------------------------------------------------
C     Loop over the combined symmetry of I AND J
C------------------------------------------------
C
CPi   Outer loop with occupied as in so_init
      DO 100 ISYMIJ = 1,NSYM
C
         ISYMAB = MULD2H(ISYMIJ,ISYRES)
C
C         PRINT *,'so_res_cdt: sym',ISYRES,ISYMIJ,ISYMAB
C---------------------------------
C        Allocation of work space.
C---------------------------------
C
         LFIJ1   = NTOO(ISYMIJ)
         LFIJ2   = NTOO(ISYMIJ)
         LFIJ3   = NSOO(ISYMIJ)
         LFAB1   = NTVV(ISYMAB)
         LFAB2   = NSVV(ISYMAB)
         LFAB3   = NTVV(ISYMAB)
C
         KFIJ1  = 1
         KFIJ2  = KFIJ1 + LFIJ1
         KFIJ3  = KFIJ2 + LFIJ2
         KFAB1  = KFIJ3 + LFIJ3
         KFAB2  = KFAB1 + LFAB1
         KFAB3  = KFAB2 + LFAB2
C
         KEND    = KFAB3  + LFAB3
         LWORK1  = LWORK - KEND
C
         CALL SO_MEMMAX ('SO_RES_CDT.1',LWORK1)
         IF (LWORK1 .LT. 0) CALL STOPIT('SO_RES_CDT.1',' ',KEND,LWORK)
C
C----------------------------------------------------------------
C        Make sum of fock-diagonals I and J in WORK(KFIJ).
C----------------------------------------------------------------
C
         DO 201 ISYMI = 1,NSYM
C
            ISYMJ = MULD2H(ISYMI,ISYMIJ)
C
            IOFFSYM1 = ITOO(MAX(ISYMI,ISYMJ),MIN(ISYMI,ISYMJ))
            IOFFSYM2 = ITOO(MAX(ISYMI,ISYMJ),MIN(ISYMI,ISYMJ))
            IOFFSYM3 = ISOO(MAX(ISYMI,ISYMJ),MIN(ISYMI,ISYMJ))
C
            IOFFIJ1 = 1
            IOFFIJ2 = 1
            IOFFIJ3 = 1
C
            IF (ISYMI .EQ. ISYMJ) THEN
C
               DO 202 I = 1,NRHF(ISYMI)
C
                  KOFFI = IRHF(ISYMI) + I
C
                  DO 203 J = 1,I ! i>=j T3
C
                     KOFFJ = IRHF(ISYMJ) + J
C
                     NIJ3 = IOFFSYM3 + IOFFIJ3 - 1
                     WORK(KFIJ3+NIJ3) =  FOCKD(KOFFI) + FOCKD(KOFFJ)
                     IOFFIJ3 = IOFFIJ3 + 1
C
  203             CONTINUE
C
                  DO 204 J = 1,(I-1) ! i>j T1,T2
C
                     KOFFJ = IRHF(ISYMJ) + J
C
                     NIJ1 = IOFFSYM1 + IOFFIJ1 - 1
                     WORK(KFIJ1+NIJ1) =  FOCKD(KOFFI) + FOCKD(KOFFJ)
                     IOFFIJ1 = IOFFIJ1 + 1
C
                     NIJ2 = IOFFSYM2 + IOFFIJ2 - 1
                     WORK(KFIJ2+NIJ2) =  FOCKD(KOFFI) + FOCKD(KOFFJ)
                     IOFFIJ2 = IOFFIJ2 + 1
C
  204             CONTINUE
C
  202          CONTINUE
C
            ELSE IF (ISYMI .GT. ISYMJ) THEN
C
               DO 205 I = 1,NRHF(ISYMI)
C
                  KOFFI = IRHF(ISYMI) + I
C
                  DO 206 J = 1,NRHF(ISYMJ)
C
                     KOFFJ = IRHF(ISYMJ) + J
C
                     NIJ1 = IOFFSYM1 + IOFFIJ1 - 1
                     WORK(KFIJ1+NIJ1) =  FOCKD(KOFFI) + FOCKD(KOFFJ)
                     IOFFIJ1 = IOFFIJ1 + 1
C
                     NIJ2 = IOFFSYM2 + IOFFIJ2 - 1
                     WORK(KFIJ2+NIJ2) =  FOCKD(KOFFI) + FOCKD(KOFFJ)
                     IOFFIJ2 = IOFFIJ2 + 1
C
                     NIJ3 = IOFFSYM3 + IOFFIJ3 - 1
                     WORK(KFIJ3+NIJ3) =  FOCKD(KOFFI) + FOCKD(KOFFJ)
                     IOFFIJ3 = IOFFIJ3 + 1
C
  206             CONTINUE
C
  205          CONTINUE
C
            END IF
C
  201    CONTINUE
C
C----------------------------------------------------------------
C        Make sum of fock-diagonals A and B in WORK(KFAB).
C----------------------------------------------------------------
C
         DO 301 ISYMA = 1,NSYM
C
            ISYMB = MULD2H(ISYMA,ISYMAB)
C
            IOFFSYM1 = ITVV(MAX(ISYMA,ISYMB),MIN(ISYMA,ISYMB))
            IOFFSYM2 = ISVV(MAX(ISYMA,ISYMB),MIN(ISYMA,ISYMB))
            IOFFSYM3 = ITVV(MAX(ISYMA,ISYMB),MIN(ISYMA,ISYMB))
C
            IOFFAB1 = 1
            IOFFAB2 = 1
            IOFFAB3 = 1
C
            IF (ISYMA .EQ. ISYMB) THEN
C
               DO 302 A = 1,NVIR(ISYMA)
C
                  KOFFA = IVIR(ISYMA) + A
C
                  DO 303 B = 1,A
C
                     KOFFB = IVIR(ISYMB) + B
C
                     NAB2 = IOFFSYM2 + IOFFAB2 - 1
                     WORK(KFAB2+NAB2) = FOCKD(KOFFA) + FOCKD(KOFFB)
                     IOFFAB2 = IOFFAB2 + 1
C
  303             CONTINUE
C
                  DO 304 B = 1,(A-1)
C
                     KOFFB = IVIR(ISYMB) + B
C
                     NAB1 = IOFFSYM1 + IOFFAB1 - 1
                     WORK(KFAB1+NAB1) = FOCKD(KOFFA) + FOCKD(KOFFB)
                     IOFFAB1 = IOFFAB1 + 1
C
                     NAB3 = IOFFSYM3 + IOFFAB3 - 1
                     WORK(KFAB3+NAB3) = FOCKD(KOFFA) + FOCKD(KOFFB)
                     IOFFAB3 = IOFFAB3 + 1
C
  304             CONTINUE
C
  302          CONTINUE
C
            ELSE IF (ISYMA .GT. ISYMB) THEN
C
               DO 305 A = 1,NVIR(ISYMA)
C
                  KOFFA = IVIR(ISYMA) + A
C
                  DO 306 B = 1,NVIR(ISYMB)
C
                     KOFFB = IVIR(ISYMB) + B
C
                     NAB1 = IOFFSYM1 + IOFFAB1 - 1
                     WORK(KFAB1+NAB1) = FOCKD(KOFFA) + FOCKD(KOFFB)
                     IOFFAB1 = IOFFAB1 + 1
C
                     NAB2 = IOFFSYM2 + IOFFAB2 - 1
                     WORK(KFAB2+NAB2) = FOCKD(KOFFA) + FOCKD(KOFFB)
                     IOFFAB2 = IOFFAB2 + 1
C
                     NAB3 = IOFFSYM3 + IOFFAB3 - 1
                     WORK(KFAB3+NAB3) = FOCKD(KOFFA) + FOCKD(KOFFB)
                     IOFFAB3 = IOFFAB3 + 1
C
  306             CONTINUE
C
  305          CONTINUE
C
            END IF
C
  301    CONTINUE
C
C---------------------------------------------------------------
C        Multiply energy-differences EABIJ and 2p2h trialvectors
C        to obtain the D-matrix contribution to the 2p2h result-
C        vectors.
C---------------------------------------------------------------
C
         IOFF1 = IT2AMT1(ISYMIJ,ISYMAB)
         IOFF2 = IT2AMT2(ISYMIJ,ISYMAB)
         IOFF3 = IT2AMT3(ISYMIJ,ISYMAB)
C
         IF (DO_DEX) THEN
C
            DO 501 NIJ = 1,NTOO(ISYMIJ) ! Nr of occ/occ pairs with i>j
C
               IOFFIJ1 = IOFF1 + NTVV(ISYMAB) * (NIJ - 1)
               IOFFIJ2 = IOFF2 + NSVV(ISYMAB) * (NIJ - 1)
C
C----------------------------------------------------------------
C              T1 contribution with i>j and a>b
C----------------------------------------------------------------
C
               DO 502 NAB = 1,NTVV(ISYMAB) !Nr of vir/vir pairs with a>b
C
                  NABIJ1 = IOFFIJ1 + NAB
C
                  EABIJ1 = WORK(KFAB1 + NAB - 1) - WORK(KFIJ1 + NIJ - 1)
C
                  RES2E(NABIJ1) = EABIJ1 * RES2E(NABIJ1)
                  RES2D(NABIJ1) = EABIJ1 * RES2D(NABIJ1)
C
  502          CONTINUE
C
C----------------------------------------------------------------
C              T2 contribution with i>j and a>=b
C----------------------------------------------------------------
C
               DO 503 NAB = 1,NSVV(ISYMAB) !Nr of vir/vir pairs with a>=b
C
                  NABIJ2 = IOFFIJ2 + NAB
C
                  EABIJ2 = WORK(KFAB2 + NAB - 1) - WORK(KFIJ2 + NIJ - 1)
C
                  RES2E(NABIJ2) = EABIJ2 * RES2E(NABIJ2)
                  RES2D(NABIJ2) = EABIJ2 * RES2D(NABIJ2)
C
  503          CONTINUE
C
  501       CONTINUE ! Loop over occ/occ pairs with i>j
C
C----------------------------------------------------------------
C              T3 contribution with i>=j and a>b
C----------------------------------------------------------------
C
            DO 504 NIJ = 1,NSOO(ISYMIJ) ! Nr of occ/occ pairs with i>=j
C
               IOFFIJ3 = IOFF3 + NTVV(ISYMAB) * (NIJ - 1)
C
               DO 505 NAB = 1,NTVV(ISYMAB) !Nr of vir/vir pairs with a>b
C
                  NABIJ3 = IOFFIJ3 + NAB
C
                  EABIJ3 = WORK(KFAB3 + NAB - 1) - WORK(KFIJ3 + NIJ - 1)
C
                  RES2E(NABIJ3) = EABIJ3 * RES2E(NABIJ3)
                  RES2D(NABIJ3) = EABIJ3 * RES2D(NABIJ3)
C
  505          CONTINUE
C
  504       CONTINUE
C
         ELSE ! Static case, one vector only
C            
            DO 601 NIJ = 1,NTOO(ISYMIJ) ! Nr of occ/occ pairs with i>j
C
               IOFFIJ1 = IOFF1 + NTVV(ISYMAB) * (NIJ - 1)
               IOFFIJ2 = IOFF2 + NSVV(ISYMAB) * (NIJ - 1)
C
C----------------------------------------------------------------
C              T1 contribution with i>j and a>b
C----------------------------------------------------------------
C
               DO 602 NAB = 1,NTVV(ISYMAB) !Nr of vir/vir pairs with a>b
C
                  NABIJ1 = IOFFIJ1 + NAB
C
                  EABIJ1 = WORK(KFAB1 + NAB - 1) - WORK(KFIJ1 + NIJ - 1)
C
                  RES2E(NABIJ1) = EABIJ1 * RES2E(NABIJ1)
C
  602          CONTINUE
C
C----------------------------------------------------------------
C              T2 contribution with i>j and a>=b
C----------------------------------------------------------------
C
               DO 603 NAB = 1,NSVV(ISYMAB) !Nr of vir/vir pairs with a>=b
C
                  NABIJ2 = IOFFIJ2 + NAB
C
                  EABIJ2 = WORK(KFAB2 + NAB - 1) - WORK(KFIJ2 + NIJ - 1)
C
                  RES2E(NABIJ2) = EABIJ2 * RES2E(NABIJ2)
C
  603          CONTINUE
C
  601       CONTINUE ! Loop over occ/occ pairs with i>j
C
C----------------------------------------------------------------
C              T3 contribution with i>=j and a>b
C----------------------------------------------------------------
C
            DO 604 NIJ = 1,NSOO(ISYMIJ) ! Nr of occ/occ pairs with i>=j
C
               IOFFIJ3 = IOFF3 + NTVV(ISYMAB) * (NIJ - 1)
C
               DO 605 NAB = 1,NTVV(ISYMAB) !Nr of vir/vir pairs with a>b
C
                  NABIJ3 = IOFFIJ3 + NAB
C
                  EABIJ3 = WORK(KFAB3 + NAB - 1) - WORK(KFIJ3 + NIJ - 1)
C
                  RES2E(NABIJ3) = EABIJ3 * RES2E(NABIJ3)
C
  605          CONTINUE
C
  604       CONTINUE


         END IF

  100 CONTINUE ! Loop over ISYMIJ (and corresponding ISYMAB)

C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_RES_CDT')
C
      RETURN
      END
