C  /* Deck so_t2m1 */
      SUBROUTINE SO_T2M1(T2M1,LT2M1,T2MP,LT2MP,CMO,LCMO,IDEL,
     &                   ISYMD,ISYDIS,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak and Henrik Koch, September 1995
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C     Rasmus Faber 2016: Rewrite for better memory access
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
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_T2M1')
C
      ISYMB = ISYMD
C
      LCDB  = NVIR(ISYMB)
CPi If no virtual orbitals in this symmetry, skip this transformation
      IF (LCDB .EQ. 0) THEN
        CALL QEXIT('SO_T2M1')
        RETURN
      END IF
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
C  /* Deck so_x2m1 */
      SUBROUTINE SO_X2M1(X2M1,LX2M1,X2MP,LX2MP,CMO,LCMO,IDEL,
     &                   ISYMD,ISYMTR,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Rasmus Faber 2016: Based one SO_T2M1, handles non-totally
C                        symmetric vectors
C                        ~
C     PURPOSE: Calculate X2-trial vectors with one back-transformed
C              index.
C
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DOUBLE PRECISION, PARAMETER :: INVSQ2 = SQRT(HALF)
C
      DIMENSION X2M1(LX2M1), X2MP(LX2MP), CMO(LCMO)
      DIMENSION WORK(LWORK)
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_X2M1')
C
      ISYMB = ISYMD
C
      LCDB  = NVIR(ISYMB)
C
      KCDB    = 1
      KEND1   = KCDB  + LCDB
      LWORK1  = LWORK - KEND1
C
      CALL SO_MEMMAX ('SO_X2M1.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_X2M1.1',' ',KEND1,LWORK)
C
      KOFF1 = ILMVIR(ISYMD) + IDEL - IBAS(ISYMD)
C
C--------------------------------------------------
C     Copy delta MO-coefficients to the vector CDB.
C--------------------------------------------------
C
      CALL DCOPY(NVIR(ISYMB),CMO(KOFF1),NBAS(ISYMD),WORK(KCDB),1)
C
C     Incorporate factor 1/sqrt(2)
      CALL DSCAL(NVIR(ISYMB),INVSQ2,WORK(KCDB),1)
C
C------------------------------------------------------------------
C     Loop over symmetry-blocks, do reductions
C------------------------------------------------------------------
C
      IF (ISYMTR .EQ. 1) THEN
C$OMP PARALLEL PRIVATE(ISYMBJ,ISYMJ,ISYMAI,ISYMI,ISYMA,
C$OMP&                 NAIBJ1,NAIBJ,NAIBJ2,
C$OMP&                 FACT,KOFF,TMP,ISTART,
C$OMP&                 NBJ,NBJ1,NAI,
C$OMP&                 J,B,I,A)
         DO 100 ISYMJ = 1,NSYM
C
            ISYMBJ = MULD2H(ISYMJ,ISYMB)
C
            ISYMAI  = ISYMBJ
C
            NAIBJ1 = IT2AM(ISYMAI,ISYMBJ)
C$OMP    DO
            DO 120 J = 1,NRHF(ISYMJ)
C
               KOFF = IT2BCD(ISYMAI,ISYMJ) + NT1AM(ISYMAI)*(J-1)
               NBJ1 = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J-1)
C
C              Zero output array
C           (only the part accessed in ai<bj loop)
               DO NAI = 1, NBJ1+NVIR(ISYMB)
                  X2M1(KOFF+NAI) = ZERO
               END DO
C
C              Do reduction for this symmetry-block
C              X2 (ai<=bj) * C(b) -> X2M1 (ai,j)
C
               DO 130 B = 1,NVIR(ISYMB)
C
                  NBJ  = NBJ1 + B
                  FACT = WORK(KCDB-1+B)
                  NAIBJ2 = NAIBJ1 + NBJ*(NBJ-1)/2
C
C                 Stride 1 loop when ai <= bj
                  DO NAI = 1, NBJ
                      NAIBJ = NAIBJ2 + NAI
                      X2M1(KOFF+NAI) = FACT*X2MP(NAIBJ)+X2M1(KOFF+NAI)
                  END DO
C
  130          CONTINUE ! LOOP B
C
C              Do the same for ai>bi :
C              X2(bj<ai)*C(b) -> X2M1 (ai,j)
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
                           TMP = TMP + X2MP(NAIBJ)*WORK(KCDB-1+B)
                        END DO
                        X2M1(KOFF+NAI) = X2M1(KOFF+NAI) + TMP
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
                           TMP = TMP + X2MP(NAIBJ)*WORK(KCDB-1+B)
                        END DO
                        X2M1(KOFF+NAI) = TMP
                     END DO
                  END DO
               END DO
C
  120       CONTINUE ! LOOP J
C$OMP       END DO NOWAIT
  100    CONTINUE ! LOOP ISYMJ
C$OMP    END PARALLEL
      ELSE
C$OMP    PARALLEL PRIVATE(ISYMBJ,ISYMJ,ISYMAI,
C$OMP&                    FACT,DSUM,
C$OMP&                    KOFF,IOFFAIBJ,KOFFBJ,KOUT,KINP,KOFFINP,
C$OMP&                    NBJ,NBJ1,NAI,
C$OMP&                    J,B)
         KOFF = 0
         DO ISYMJ = 1, NSYM
            ISYMBJ = MULD2H(ISYMJ,ISYMB)
            ISYMAI = MULD2H(ISYMBJ,ISYMTR)
            IOFFAIBJ = IT2AM(ISYMAI,ISYMBJ)
            KOFFBJ = IT1AM(ISYMB,ISYMJ)

            IF (ISYMAI.LT.ISYMBJ) THEN

               KOFFINP = IOFFAIBJ + NT1AM(ISYMAI)*KOFFBJ
C$OMP          DO
               DO J = 1, NRHF(ISYMJ)
                  KOUT = KOFF + NT1AM(ISYMAI)*(J-1)+1
                  KINP = KOFFINP + NT1AM(ISYMAI)*NVIR(ISYMB)*(J-1)+1
                  CALL DGEMV('N',NT1AM(ISYMAI),NVIR(ISYMB),ONE,
     &                       X2MP(KINP),MAX(1,NT1AM(ISYMAI)),
     &                       WORK(KCDB),1,ZERO,X2M1(KOUT),1)

               END DO
C$OMP          END DO NOWAIT

            ELSE ! ISYMAI > ISYMBJ

C$OMP          DO
               DO J = 1,NRHF(ISYMJ)
                  KOUT = KOFF + NT1AM(ISYMAI)*(J-1)
                  KOFFINP = IOFFAIBJ + KOFFBJ + NVIR(ISYMB)*(J-1)
                  DO NAI = 1,NT1AM(ISYMAI)
                     DSUM = 0.0D0
                     DO B = 1, NVIR(ISYMB)
                        DSUM = DSUM + X2MP(KOFFINP+B)*WORK(KCDB-1+B)
                     END DO
                     X2M1(KOUT+NAI) = DSUM
                     KOFFINP = KOFFINP + NT1AM(ISYMBJ)
                  END DO
               END DO
C$OMP          END DO NOWAIT
            END IF

            KOFF = KOFF + NT1AM(ISYMAI)*NRHF(ISYMJ)

         END DO
C$OMP    END PARALLEL
      ENDIF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_X2M1')
C
      RETURN
      END
C  /* Deck so_x2sm1 */
      SUBROUTINE SO_X2SM1(X2M1,LX2M1,X2SQ,LX2SQ,CMO,LCMO,IDEL,
     &                    ISYMD,ISYMTR,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Rasmus Faber 2017: Version of SO_X2M1, that handles squared, non-totally
C                        symmetric vectors
C                        ~
C     PURPOSE: Calculate X2-trial vectors with one back-transformed
C              index.
C
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DOUBLE PRECISION, PARAMETER :: INVSQ2 = SQRT(HALF)
      character(len=*), parameter :: myname = 'SO_X2SM1'
C
      DIMENSION X2M1(LX2M1), X2SQ(LX2SQ), CMO(LCMO)
      DIMENSION WORK(LWORK)
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER(myname)
C
      ISYMB = ISYMD
C
      LCDB  = NVIR(ISYMB)
C
      KCDB    = 1
      KEND1   = KCDB  + LCDB
      LWORK1  = LWORK - KEND1
C
      CALL SO_MEMMAX (myname//'.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT(myname//'.1',' ',KEND1,LWORK)
C
      KOFF1 = ILMVIR(ISYMD) + IDEL - IBAS(ISYMD)
C
C--------------------------------------------------
C     Copy delta MO-coefficients to the vector CDB.
C--------------------------------------------------
C
      CALL DCOPY(NVIR(ISYMB),CMO(KOFF1),NBAS(ISYMD),WORK(KCDB),1)
C
C     Incorporate factor 1/sqrt(2)
      CALL DSCAL(NVIR(ISYMB),INVSQ2,WORK(KCDB),1)
C
C------------------------------------------------------------------
C     Loop over symmetry-blocks, do reductions
C------------------------------------------------------------------
C
      KOFF = 0
      DO ISYMJ = 1, NSYM
         ISYMBJ = MULD2H(ISYMJ,ISYMB)
         ISYMAI = MULD2H(ISYMBJ,ISYMTR)
         IOFFAIBJ = IT2SQ(ISYMAI,ISYMBJ)
         KOFFBJ = IT1AM(ISYMB,ISYMJ)


         KOFFINP = IOFFAIBJ + NT1AM(ISYMAI)*KOFFBJ
         DO J = 1, NRHF(ISYMJ)
            KOUT = KOFF + NT1AM(ISYMAI)*(J-1)+1
            KINP = KOFFINP + NT1AM(ISYMAI)*NVIR(ISYMB)*(J-1)+1
            CALL DGEMV('N',NT1AM(ISYMAI),NVIR(ISYMB),ONE,
     &                 X2SQ(KINP),MAX(1,NT1AM(ISYMAI)),
     &                 WORK(KCDB),1,ONE,X2M1(KOUT),1)
         END DO


         KOFF = KOFF + NT1AM(ISYMAI)*NRHF(ISYMJ)

      END DO
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT(myname)
C
      RETURN
      END

C  /* Deck so_t2x1 */
      SUBROUTINE SO_T2X1(X2M1,LX2M1,T2MP,LT2MP,TJ1,LTJ1,IDEL,
     &                   ISYMD,ISYMTR,WORK,LWORK)
C
C     Rasmus Faber, 2016: Bases on SO_T2M1
C
C     PURPOSE: Calculates a T2(ai, bj) * x1( b, delta) intermediate and
C              and adds it to the back-transformed x2 trial-vectors.
C              This allows term (1) and (5) of the B-matrix to be
C              calculated using SO_RES_TCB
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION X2M1(LX2M1), T2MP(LT2MP), TJ1(LTJ1)
      DIMENSION WORK(LWORK)
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_T2X1')
C
      ISYMB = MULD2H(ISYMD,ISYMTR)
C
      LCDB  = NVIR(ISYMB)
CRF   If there's nothing to do here, skip straight to exit
      IF ( LCDB .LE. 0 ) GOTO 4321
C
      KCDB   = 1
      KEND1  = KCDB  + LCDB
      LWORK1 = LWORK - KEND1
C
      CALL SO_MEMMAX ('SO_T2X1.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_T2X1.1',' ',KEND1,LWORK)
C
      KOFF1 = IMATAV(ISYMD,ISYMB) + IDEL - IBAS(ISYMD)
C
C--------------------------------------------------
C     Copy delta X1-coefficients to the vector CDB.
C--------------------------------------------------
C
      CALL DCOPY(NVIR(ISYMB),TJ1(KOFF1),NBAS(ISYMD),WORK(KCDB),1)
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
         ISYMAI = ISYMBJ
C
         NAIBJ1 = IT2AM(ISYMAI,ISYMBJ)
C$OMP DO
         DO 120 J = 1,NRHF(ISYMJ)
C
            KOFF = IT2BCD(ISYMAI,ISYMJ) + NT1AM(ISYMAI)*(J-1)
            NBJ1 = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J-1)
C
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
                   X2M1(KOFF+NAI) = FACT*T2MP(NAIBJ)+X2M1(KOFF+NAI)
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
                     X2M1(KOFF+NAI) = X2M1(KOFF+NAI) + TMP
                  END DO
                  ! I > J taken care of in the following
                  ISTART = J + 1
               ELSE
                  ! ISYMI > ISYMJ, so I>J is implicit, start loop from
                  ! one
                  ISTART = 1
               END IF
C
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
                     X2M1(KOFF+NAI) = X2M1(KOFF+NAI) + TMP
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
4321  CALL QEXIT('SO_T2X1')
C
      RETURN
      END

      SUBROUTINE SO_M1SHUF(T2M1,LT2M1,ISYMD,ISYMTR)
C
C     Calculate
C        v                ~           ~
C        t(ai,j) = -1/3 ( t(ai,j) + 2 t(aj,i) )
C     in place.
      use so_info, only: sop_dp
      implicit none

#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
      real(sop_dp), intent(inout) :: T2M1(LT2M1)
      integer, intent(in) :: LT2M1, ISYMD, ISYMTR

      character(len=*), parameter :: myname = 'SO_M1SHUF'
      real(sop_dp), parameter :: inv3m = -1.0D0/3.0D0, two = 2.0D0

      integer :: isymi, isymj, isymai, isyma, isymjd, isymaj
      integer :: ioffj, ioffi, ioffaij, ioffaji
      integer :: nrhfi
      real(sop_dp) :: xaij, xaji

      CALL QENTER(MYNAME)

      DO ISYMJ = 1, NSYM
         isymjd = muld2h(isymj,isymd)
         do ISYMI = 1, isymj ! Ensure i < j
            isymai = muld2h(isymjd,isymtr)
            isyma = muld2h(isymai,isymi)
            isymaj = muld2h(isymj,isyma)

            do j = 1, nrhf(isymj)
               ioffj = it2bcd(isymai,isymj) + nt1am(isymai)*(j-1)
               if (isymi.eq.isymj) then
                  nrhfi = j
               else
                  nrhfi = nrhf(isymi)
               end if
               do i = 1, nrhfi
                  ioffi = it2bcd(isymaj,isymi)+nt1am(isymaj)*(i-1)
                  ioffaij = ioffj + it1am(isyma,isymi)
     &                    + nvir(isyma)*(i-1)
                  ioffaji = ioffi + it1am(isyma,isymj)
     &                    + nvir(isyma)*(j-1)
                  do a = 1, nvir(isyma)
                     xaij = t2m1(ioffaij+a)
                     xaji = t2m1(ioffaji+a)
                     t2m1(ioffaij+a) = inv3m * (xaij+two*xaji)
                     t2m1(ioffaji+a) = inv3m * (xaji+two*xaij)
                  end do ! a
               end do ! i
            end do ! j
         end do ! isymi
      END DO ! isymj

      CALL QEXIT(MYNAME)
      END SUBROUTINE
