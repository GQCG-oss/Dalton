C  /* Deck so_t2m1 */
      SUBROUTINE SO_T2M1(T2M1,LT2M1,T2MP,LT2MP,CMO,LCMO,IDEL,
     &                   ISYMD,ISYDIS,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak and Henrik Koch, September 1995
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Calculate MP2 T2-amplitudes with one back-transformed
C              index.
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION T2M1(LT2M1), T2MP(LT2MP), CMO(LCMO)
      DIMENSION WORK(LWORK)
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
C
C------------------------------
C     Statement function INDEX.
C------------------------------
C
C      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_T2M1')
C
      ISYMB = ISYMD
C
      LCDB  = NVIR(ISYMB)
CRF   If there's nothing to do here, skip straight to exit
      IF ( LCDB .LE. 0 ) GOTO 4321
C
      KCDB    = 1
      KEND1   = KCDB  + LCDB
      LWORK1  = LWORK - KEND1
C
      CALL SO_MEMMAX ('SO_T2M1.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_T2M1.1',' ',KEND1,LWORK)
C
      KOFF1 = ILMVIR(ISYMD) + IDEL - IBAS(ISYMD)
C
C--------------------------------------------------
C     Copy delta MO-coefficients to the vector CDB.
C--------------------------------------------------
C
      CALL DCOPY(NVIR(ISYMB),CMO(KOFF1),NBAS(ISYMD),WORK(KCDB),1)
C
C------------------------------------------------------------------
C     Loop over symmetry-blocks, do reductions
C------------------------------------------------------------------
C
C$OMP PARALLEL PRIVATE(ISYMBJ,ISYMJ,ISYMAI,NAIBJ1,NAIBJ,FACT,NBJ,NBJ1,
C$OMP&                 J,B,NAI,NAIBJ2,ISYMI,ISYMA,I,A,KOFFBJ,TMP,ISTART)
      DO 100 ISYMJ = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMJ,ISYMB)
C
         ISYMAI  = ISYMBJ
C
         NAIBJ1 = IT2AM(ISYMAI,ISYMBJ)
C$OMP DO
         DO 120 J = 1,NRHF(ISYMJ)
C
            KOFF = IT2BCD(ISYMAI,ISYMJ) + NT1AM(ISYMAI)*(J-1)
            NBJ1 = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J-1)
C
C           Zero output array 
C           (only the part accessed in ai<bj loop)
            DO NAI = 1, NBJ1+NVIR(ISYMB)
               T2M1(KOFF+NAI) = ZERO
            END DO
C
C           Do reduction for this symmetry-block
C           T2 (ai<=bj) * C(b) -> T2M1 (ai,j)
C
            DO 130 B = 1,NVIR(ISYMB)
C
               NBJ  = NBJ1 + B
               FACT = WORK(KCDB-1+B)
               NAIBJ2 = NAIBJ1 + NBJ*(NBJ-1)/2
C
C              Stride 1 loop when ai <= bj
               DO NAI = 1, NBJ
                   NAIBJ = NAIBJ2 + NAI
                   T2M1(KOFF+NAI) = FACT*T2MP(NAIBJ)+T2M1(KOFF+NAI)
               END DO
C
  130       CONTINUE ! LOOP B
C
C           Do the same for ai>bi :
C           T2(bj<ai)*C(b) -> T2M1 (ai,j)
            DO ISYMI = ISYMJ, NSYM
               ISYMA = MULD2H(ISYMAI,ISYMI)
               IF (ISYMI .EQ. ISYMJ) THEN
                  ! Same symmetries: two cases I>J and I=J
                  ! DO I = J Seperately
                  I = J
                  DO A = 2, NVIR(ISYMA)
                     NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1)+A
                     TMP = ZERO
                     KOFFBJ = NAIBJ1 + NAI*(NAI-1)/2 + NBJ1
                     ! In this case B<A, The B>=A has allready been
                     ! handled
                     DO B = 1, A-1
                        NAIBJ = KOFFBJ + B
                        TMP = TMP + T2MP(NAIBJ)*WORK(KCDB-1+B)
                     END DO
                     T2M1(KOFF+NAI) = T2M1(KOFF+NAI) + TMP
                  END DO
                  ! I > J taken care of in the following
                  ISTART = J + 1
               ELSE
                  ! ISYMI > ISYMJ, so I>J is implicit, start loop from
                  ! one
                  ISTART = 1
               END IF

               DO I = ISTART, NRHF(ISYMI)
                  ! Since I>J, we have no restrictions on A,B
                  DO A = 1, NVIR(ISYMA)
                     NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1)+A
                     TMP = ZERO
                     KOFFBJ = NAIBJ1 + NAI*(NAI-1)/2 + NBJ1
                     ! This could be a DDOT
                     DO B = 1, NVIR(ISYMB)
                        NAIBJ = KOFFBJ + B
                        TMP = TMP + T2MP(NAIBJ)*WORK(KCDB-1+B)
                     END DO
                     T2M1(KOFF+NAI) = TMP
                  END DO
               END DO
            END DO
C
  120    CONTINUE ! LOOP J
C$OMP END DO NOWAIT
  100 CONTINUE ! LOOP ISYMJ
C$OMP END PARALLEL
C
C-----------------------
C     Remove from trace.
C-----------------------
C
4321  CALL QEXIT('SO_T2M1')
C
      RETURN
      END
