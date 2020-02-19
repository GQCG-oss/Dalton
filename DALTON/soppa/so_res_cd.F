C  /* Deck so_res_cd */
      SUBROUTINE SO_RES_CD(RES2E,LRES2E,RES2D,LRES2D,TR2E,LTR2E,TR2D,
     &                     LTR2D,FOCKD,LFOCKD,ISYRES,DO_DEX,
     &                     WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, April 1996
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C     Rasmus Faber: April 2016 - Now sets rather than increments
C                                Result-Vector
C
C     PURPOSE: Calculate contribution from D matrix to the 2p2h
C              result vectors.
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION RES2E(LRES2E), RES2D(LRES2D)
      DIMENSION TR2E(LTR2E),   TR2D(LTR2D)
      DIMENSION FOCKD(LFOCKD), WORK(LWORK)
      LOGICAL, INTENT(IN) :: DO_DEX
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "soppinf.h"
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
C        Since nothing is done in the end, skip this it.
         IF (ISYMAI.GT.ISYMBJ) GOTO 100
C
C---------------------------------
C        Allocation of work space.
C---------------------------------
C
         LFBJ   = NT1AM(ISYMBJ)
         LFAI   = NT1AM(ISYMAI)
C
         KFBJ    = 1
         IF(ISYMAI.LT.ISYMBJ) THEN
            KFAI    = KFBJ  + LFBJ
         ELSE ! They're the same!!!
            KFAI = KFBJ
         ENDIF
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
            IOFFSYM = IT1AM(ISYMB,ISYMJ)
C
            DO 202 J = 1,NRHF(ISYMJ)
C
               KOFFJ = IRHF(ISYMJ) + J
C
               DO 203 B = 1,NVIR(ISYMB)
C
                  NBJ   = IOFFSYM + NVIR(ISYMB)*(J - 1) + B - 1
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
CRF      Only if it is different than the calculated above!!!
         IF (ISYMAI.LT.ISYMBJ) THEN
            DO  301 ISYMI = 1,NSYM
C
               ISYMA = MULD2H(ISYMI,ISYMAI)
C
               IOFFSYM = IT1AM(ISYMA,ISYMI)
               DO 302 I = 1,NRHF(ISYMI)
C
                  KOFFI = IRHF(ISYMI) + I
C
                  DO 303 A = 1,NVIR(ISYMA)
C
                     NAI = IOFFSYM + NVIR(ISYMA)*(I - 1) + A - 1
C
                     KOFFA = IVIR(ISYMA) + A
C
                     WORK(KFAI+NAI) =  FOCKD(KOFFA) - FOCKD(KOFFI)
C
  303             CONTINUE
C
  302          CONTINUE
C
  301       CONTINUE
         ENDIF
C
C---------------------------------------------------------------
C        Multiply energy-differences EAIBJ and 2p2h trialvectors
C        to obtain the D-matrix contribution to the 2p2h result-
C        vectors.
C---------------------------------------------------------------
C
         IF ( ISYMAI .EQ. ISYMBJ) THEN
C
            IOFF = IT2AM(ISYMAI, ISYMBJ) + 1
            CALL EXTR2_TOTSYM( TR2E(IOFF), RES2E(IOFF),
     &                         WORK(KFBJ), ISYMBJ)
            IF (DO_DEX) THEN
               CALL EXTR2_TOTSYM( TR2D(IOFF), RES2D(IOFF),
     &                            WORK(KFBJ), ISYMBJ)
            END IF
         ELSE IF ( ISYMAI .LT. ISYMBJ) THEN
C
            IOFF = IT2AM(ISYMAI,ISYMBJ) + 1
            CALL EXTR2_ASYM(TR2E(IOFF),RES2E(IOFF),WORK(KFAI),
     &                      WORK(KFBJ),ISYMAI,ISYMBJ)
            IF (DO_DEX) THEN
               CALL EXTR2_ASYM(TR2D(IOFF),RES2D(IOFF),WORK(KFAI),
     &                         WORK(KFBJ),ISYMAI,ISYMBJ)
            END IF
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
      CONTAINS
         pure SUBROUTINE EXTR2_TOTSYM ( TR2, RES2, EPAIR, ISYM)
            DOUBLE PRECISION, INTENT(IN) :: TR2(*), EPAIR(*)
            DOUBLE PRECISION, INTENT(INOUT) :: RES2(*)
            integer, intent(in) :: isym
C
            integer :: nbj, nai, naibj, ioffbj
            double precision :: eaibj
C
C            ioffsym = it2am(isym,isym)

            do nbj = 1, nt1am(isym)
               ioffbj = nbj*(nbj-1)/2
               do nai = 1, nbj
                  naibj = ioffbj + nai
                  eaibj = epair(nbj) + epair(nai)
C                  res2(naibj) = res2(naibj) + eaibj*tr2(naibj)
                  res2(naibj) =  eaibj*tr2(naibj)
               end do
            end do
            return
         end subroutine

         pure subroutine extr2_asym ( tr2, res2, epair1, epair2,
     &                                isym1, isym2 )
            double precision, intent(in) :: tr2(*), epair1(*), epair2(*)
            double precision, intent(inout) :: res2(*)
            integer, intent(in) :: isym1, isym2
C
            integer :: nbj, nai, naibj, ioffbj
            double precision :: eaibj

            do nbj = 1, nt1am(isym2)
               ioffbj = nt1am(isym1)*(nbj-1)
               do nai = 1, nt1am(isym1)
                  naibj = ioffbj + nai
                  eaibj = epair1(nai) + epair2(nbj)
C                  res2(naibj) = res2(naibj) + eaibj*tr2(naibj)
                  res2(naibj) = eaibj*tr2(naibj)
               end do
            end do

            return
         end subroutine
      END
