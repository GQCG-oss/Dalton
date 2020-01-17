C  /* Deck cc_t2sq3a */
      SUBROUTINE CC_T2SQ3A(T2AM,T2SQ,ISYM,FACT)
!
!--------------------------------------------------------
!     Rasmus Faber 2017:
!
!     Adds "FACT" times the content of the antisymmetric,
!     packed array T2AM to the square array T2SQ
!
!     Based on CC_T2SQ3 by Kasper Hald
!--------------------------------------------------------
!
      IMPLICIT NONE
!if defined (SYS_CRAY)
!     REAL T2AM(*), T2SQ(*)
!else
      DOUBLE PRECISION, PARAMETER :: ONEM = -1.0D0
      DOUBLE PRECISION, INTENT(IN) :: T2AM(*), FACT
      DOUBLE PRECISION, INTENT(INOUT) :: T2SQ(*)
!endif
      INTEGER, INTENT(IN) :: ISYM
      INTEGER KOFF1, KOFF2, ISYMBJ, KOFF, ISYMAI, NAMP, NAI
      INTEGER NBJ
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
!
      CALL QENTER('CC_T2SQ3')
!
      IF (ISYM.EQ.1) THEN
         KOFF1 = 1
         KOFF2 = 1
         DO 100 ISYMBJ = 1,NSYM
            IF (NT1AM(ISYMBJ) .GT. 0) THEN
               CALL SQMATR3A(NT1AM(ISYMBJ),T2AM(KOFF1),T2SQ(KOFF2),FACT)
               KOFF1 = KOFF1 + NT1AM(ISYMBJ)*(NT1AM(ISYMBJ)+1)/2
               KOFF2 = KOFF2 + NT1AM(ISYMBJ)*NT1AM(ISYMBJ)
            ENDIF
  100    CONTINUE
!
      ELSE
!
         KOFF = 1
         DO 200 ISYMBJ = 1,NSYM
            ISYMAI = MULD2H(ISYM,ISYMBJ)
!
            IF (ISYMBJ.LT.ISYMAI) CYCLE
!
            NAMP = NT1AM(ISYMAI)*NT1AM(ISYMBJ)
            KOFF = IT2AM(ISYMAI,ISYMBJ) + 1
!
            IF (NAMP .GT. 0) THEN
               KOFF1 = IT2SQ(ISYMAI,ISYMBJ) + 1
               CALL DAXPY(NAMP,FACT,T2AM(KOFF),1,T2SQ(KOFF1),1)
               NAI = MAX(NT1AM(ISYMAI),1)
               NBJ = MAX(NT1AM(ISYMBJ),1)
               KOFF2 = IT2SQ(ISYMBJ,ISYMAI) + 1
               CALL TRM3A(T2AM(KOFF),NT1AM(ISYMAI),NT1AM(ISYMBJ),
     *                    T2SQ(KOFF2),ONEM*FACT)
!
            ENDIF
!
  200    CONTINUE
!
      ENDIF
!
      CALL QEXIT('CC_T2SQ3')
!
      RETURN
      CONTAINS
!
      SUBROUTINE TRM3A(A,M,N,B,FACTM)
!
!---------------------------------------------------------------
!
!     Adds FACT times the transpose of A to the array B
!
!     Based on trm3 by Kasper Hald
!     Based on TRM by Ove Christiansen.
!---------------------------------------------------------------
!
      INTEGER, INTENT(IN) :: M, N
      DOUBLE PRECISION, INTENT(IN) :: A(M,N), FACTM
      DOUBLE PRECISION, INTENT(INOUT) :: B(N,M)
!
      DOUBLE PRECISION :: XMONE
      INTEGER :: I
!
      DO 100 I = 1, M
!
C         CALL DAXPY(N,FACT,A(I,1),M,B(1,I),1)
         DO J = 1, N
            B(J,I) = B(J,I) + FACTM*A(I,J)
         END DO
!
 100  CONTINUE
!
      RETURN
      END SUBROUTINE
C
      PURE SUBROUTINE SQMATR3A(NDIM,PKMAT,SQMAT,FACT)
!
!-----------------------------------------------------
!
!     This subroutine adds the packed antisymmetric matrix
!     PKMAT*FACT to the square matrix SQMAT
!
!     Based on SQMATR3 by Kasper Hald.
!     Based on SQMATR by Henrik Koch.
!-----------------------------------------------------
!
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: FACT, PKMAT(*)
      DOUBLE PRECISION, INTENT(INOUT) :: SQMAT(NDIM,NDIM)
!
      INTEGER I, J, IJ
      DOUBLE PRECISION ZERO,XMONE
!
      PARAMETER(XMONE = -1.0D00)
!
      IJ = 0
      DO 100 I = 1,NDIM
         DO 110 J = 1,I-1
               IJ = IJ + 1
               SQMAT(I,J) = -FACT*PKMAT(IJ) + SQMAT(I,J)
               SQMAT(J,I) =  FACT*PKMAT(IJ) + SQMAT(J,I)
  110    CONTINUE
         IJ = IJ + 1
         SQMAT(I,I) = 0.0D0
  100 CONTINUE
!
      RETURN
      END SUBROUTINE

      END
C  /* Deck cc_12c3 */
      SUBROUTINE CC_12C3(RHO2,CTR1,ISYMV,XLAMDH,ISYMH,
     *                   ISINT,WORK,LWORK,
     *                   LUO3,O3FIL)
C
C     Rasmus Faber - 2017
C
C     Based on CC_21B12C by Asger Halkier & Henrik Koch 13/9 - 1995.
C
C     Note: the file-read performed here is also done in
C           CC_21G: we might want to merge those routines.
C
      IMPLICIT NONE
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, ONE = 1.0D0,
     &                               ONEM = -1.0D0
      CHARACTER(len=*), PARAMETER :: myname = 'CC_12C3'
      DOUBLE PRECISION, INTENT(INOUT) :: RHO2(*)
      DOUBLE PRECISION, INTENT(IN) :: CTR1(*),XLAMDH(*)
      DOUBLE PRECISION, INTENT(OUT) :: WORK(LWORK)
      CHARACTER, INTENT(IN) :: O3FIL*(*)
      INTEGER, INTENT(IN) :: ISYMV, ISYMH, ISINT, LWORK, LUO3
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
#include "cclr.h"
C
      INTEGER :: ISYMA, ISYMB, ISYMD, ISYMI, ISYMJ, ISYMK,
     &           ISYMAI, ISYMBJ, ISYMIK, ISYIKJ,
     &           ISYMLI, ISYRES
      INTEGER :: NTOT, NTOTA, NTOTD, NTOTI, NTOIKJ
      INTEGER :: KSCR1, KSCR3, KOFF1, KOFF2, KOFF3, KEND1
      INTEGER :: LWRK1, IOFF, NBJ
C
      CALL QENTER(myname)
C
      ISYMLI = MULD2H(ISINT,ISYMH)
      ISYRES = MULD2H(ISYMV,ISYMLI)
C
      DO 100 ISYMD = 1,NSYM
C
         IF (NBAS(ISYMD) .EQ. 0) GOTO 100
C
         ISYIKJ = MULD2H(ISINT,ISYMD)
         ISYMB  = MULD2H(ISYMH,ISYMD)
C
C---------------------------------
C        Allocation of work space.
C---------------------------------
C
         KSCR1 = 1
         KSCR3 = KSCR1 + NMAIJK(ISYIKJ)*NVIR(ISYMB)
         KEND1 = KSCR3 + NMAIJK(ISYIKJ)*NBAS(ISYMD)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient work space area in '//myname)
         ENDIF
C
C------------------------------------
C        Read integrals from disc.
C------------------------------------
C
         NTOT = NMAIJK(ISYIKJ)*NBAS(ISYMD)
         IOFF = I3ODEL(ISYIKJ,ISYMD) + 1
C
         CALL GETWA2(LUO3,O3FIL,WORK(KSCR3),IOFF,NTOT)
C
C---------------------------
C        Transform AO index.
C---------------------------
C
         KOFF1  = ILMVIR(ISYMB) + 1
C
         NTOIKJ = MAX(NMAIJK(ISYIKJ),1)
         NTOTD  = MAX(NBAS(ISYMD),1)
C
         CALL DGEMM('N','N',NMAIJK(ISYIKJ),NVIR(ISYMB),NBAS(ISYMD),
     *              ONE,WORK(KSCR3),NTOIKJ,XLAMDH(KOFF1),NTOTD,
     *              ZERO,WORK(KSCR1),NTOIKJ)
C
         DO 110 B = 1,NVIR(ISYMB)
C
            DO 170 ISYMJ = 1,NSYM
C
               ISYMBJ = MULD2H(ISYMJ,ISYMB)
               ISYMAI = MULD2H(ISYMBJ,ISYRES)
               ISYMIK = MULD2H(ISYMJ,ISYIKJ)
C
               DO J = 1, NRHF(ISYMJ)
C
C-----------------------------------------------------------
C                 Contract integrals with trial vector CTR1.
C-----------------------------------------------------------
C
                  DO 190 ISYMK = 1,NSYM
C
                     ISYMI  = MULD2H(ISYMK,ISYMIK)
                     ISYMA  = MULD2H(ISYMK,ISYMV)
                     NBJ = IT1AM(ISYMB,ISYMJ)
     *                   + NVIR(ISYMB)*(J-1) + B
C
                     KOFF1 = IT1AM(ISYMA,ISYMK) + 1
                     KOFF2 = KSCR1 + NMAIJK(ISYIKJ)*(B-1)
     *                             + IMAIJK(ISYMIK,ISYMJ)
     *                             + NMATIJ(ISYMIK)*(J-1)
     *                             + IMATIJ(ISYMI,ISYMK)
                     KOFF3 = IT2SQ(ISYMAI,ISYMBJ)
     *                     + NT1AM(ISYMAI)*(NBJ-1)
     *                     + IT1AM(ISYMA,ISYMI) + 1
C
                     NTOTA = MAX(NVIR(ISYMA),1)
                     NTOTI = MAX(NRHF(ISYMI),1)
C
                     CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),
     *                          NRHF(ISYMK),ONEM,CTR1(KOFF1),NTOTA,
     *                          WORK(KOFF2),NTOTI,ONE,RHO2(KOFF3),
     *                          NTOTA)
C
  190             CONTINUE
               END DO
  170       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT(myname)
C
      RETURN
      END
C  /* Deck cc_xd3 */
      SUBROUTINE CC_XD3(XOUT,XMAT,ISYMX,XLAMDH,XLAMDP,ISYMLH,
     *                 WORK,LWORK)

C
C     Rasmus Faber 2017
C
C     Purpose: Transform X-matrix to AO-basis!
C
C     Based on CC_YD by Asger Halkier 8/12 - 1995.
C
      IMPLICIT NONE
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, ONE = 1.0D0,
     &                               ONEM = -1.0D0
      CHARACTER(LEN=*), PARAMETER :: myname = 'CC_XD3'
      DOUBLE PRECISION, INTENT(IN) :: XMAT(*), XLAMDH(*), XLAMDP(*)
      DOUBLE PRECISION, INTENT(INOUT) :: XOUT(*), WORK(LWORK)
      INTEGER, INTENT(IN) :: ISYMX, ISYMLH, LWORK

      INTEGER :: ISYMM, ISYMN, ISYMAL, ISYMBE, ISYMXD
      INTEGER :: NTOTAL, NTOTBE, NTOTM
      INTEGER :: KOFF1, KOFF2, LWRKCH

C#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "cclr.h"
C
      ISYMXD = MULD2H(ISYMLH,ISYMX)
C
C------------------------------------
C     Transform X-matrix to AO-basis.
C------------------------------------
C
      DO 100 ISYMAL = 1,NSYM
C
         ISYMM  = ISYMAL
         ISYMBE = MULD2H(ISYMAL,ISYMXD)
         ISYMN  = MULD2H(ISYMBE,ISYMLH)
C
         LWRKCH = LWORK - NBAS(ISYMAL)*NVIR(ISYMN)
C
         IF (LWRKCH .LT. 0) THEN
            CALL QUIT('Insufficient work space in '//myname)
         ENDIF
C
         KOFF1 = ILMRHF(ISYMM) + 1
         KOFF2 = IMATIJ(ISYMM,ISYMN) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTM  = MAX(NRHF(ISYMM),1)
C
C        lambdp (alpha,m) * x(m,n) => work(alpha,n)
         CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMN),NRHF(ISYMM),
     *              ONE,XLAMDP(KOFF1),NTOTAL,XMAT(KOFF2),NTOTM,
     *              ZERO,WORK,NTOTAL)
C
         KOFF1 = IGLMRH(ISYMBE,ISYMN) + 1
         KOFF2 = IAODIS(ISYMAL,ISYMBE) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
C        work(alpha,n) * lambdh(beta,n) => xout(alpha,beta)
         CALL DGEMM('N','T',NBAS(ISYMAL),NBAS(ISYMBE),NRHF(ISYMN),
     *              ONEM,WORK,NTOTAL,XLAMDH(KOFF1),NTOTBE,
     *              ONE,XOUT(KOFF2),NTOTAL)
C
  100 CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE CC_T2SQTRANSP(X2SQ,ISYMTR)
C
C     Transpose the squared doubles amplitude-like array,
C     x2 (ai,bj) -> (bj,ai)
C     Needed in the triplet case, when the squared array
C     typically contain a linear combination of (+) and (-)
C     components
      IMPLICIT NONE

#include "ccorb.h"
#include "ccsdsym.h"

      INTEGER, INTENT(IN) :: ISYMTR
      DOUBLE PRECISION, INTENT(INOUT) :: X2SQ(*)
      INTEGER :: ISYMAI, ISYMBJ, IXPOS, IXPOS2
C
      IF (ISYMTR.EQ.1) THEN

         DO ISYMBJ = 1, NSYM
            ISYMAI = ISYMBJ
            IXPOS = IT2SQ(ISYMAI,ISYMBJ) + 1

            CALL TRANSP(X2SQ(IXPOS),NT1AM(ISYMAI))

         END DO

      ELSE

          DO ISYMBJ = 1, NSYM
             ISYMAI = MULD2H(ISYMBJ,ISYMTR)
             IF (ISYMAI .GT. ISYMBJ) CYCLE
             IXPOS = IT2SQ(ISYMAI,ISYMBJ) + 1
             IXPOS2 = IT2SQ(ISYMBJ,ISYMAI) + 1
             CALL TRANSPA(X2SQ(IXPOS),X2SQ(IXPOS2),
     &                    NT1AM(ISYMAI),NT1AM(ISYMBJ))
          END DO

      END IF
C
      CONTAINS
C
         PURE SUBROUTINE TRANSP(A,N)
C
         DOUBLE PRECISION, INTENT(INOUT) :: A(N,N)
         INTEGER, INTENT(IN) :: N
         DOUBLE PRECISION :: TMP
         INTEGER I, J
C
         DO I = 1, N
            DO J = 1, I-1
               TMP = A(I,J)
               A(I,J) = A(J,I)
               A(J,I) = TMP
            END DO
         END DO
C
         END SUBROUTINE
C
         PURE SUBROUTINE TRANSPA(A,B,M,N)
C
         DOUBLE PRECISION, INTENT(INOUT) :: A(M,N), B(N,M)
         INTEGER, INTENT(IN) :: N, M
         INTEGER :: I, J
         DOUBLE PRECISION TMP
C
         DO I = 1, M
            DO J = 1, N
               TMP = A(I,J)
               A(I,J) = B(J,I)
               B(J,I) = TMP
            END DO
         END DO
C
         END SUBROUTINE

      END SUBROUTINE
C
      SUBROUTINE CC_T2SQSYMSCAL(X2SQ,ISYMTR,FACT)
C
C     Scale the symmetric part of the squared doubles array by FACT:
C     x2 (ai,bj) <- (1+FACT)/2 (ai,bj) - (1-FACT)/2 (bj,ai)
C     Needed in the triplet case, when the squared array
C     typically contain a linear combination of (+) and (-)
C     components
      IMPLICIT NONE

#include "ccorb.h"
#include "ccsdsym.h"

      INTEGER, INTENT(IN) :: ISYMTR
      DOUBLE PRECISION, INTENT(INOUT) :: X2SQ(*)
      DOUBLE PRECISION, INTENT(IN) :: FACT
      DOUBLE PRECISION, PARAMETER :: HALF = 0.5D0, ONE = 1.0D0
      INTEGER :: ISYMAI, ISYMBJ, IXPOS, IXPOS2
      DOUBLE PRECISION :: FN, FT
C
      FN = HALF*(ONE+FACT)
      FT = HALF*(ONE-FACT)
C
      IF (ISYMTR.EQ.1) THEN

         DO ISYMBJ = 1, NSYM
            ISYMAI = ISYMBJ
            IXPOS = IT2SQ(ISYMAI,ISYMBJ) + 1

            CALL SCAL1(X2SQ(IXPOS),NT1AM(ISYMAI))

         END DO

      ELSE

          DO ISYMBJ = 1, NSYM
             ISYMAI = MULD2H(ISYMBJ,ISYMTR)
             IF (ISYMAI .GT. ISYMBJ) CYCLE
             IXPOS = IT2SQ(ISYMAI,ISYMBJ) + 1
             IXPOS2 = IT2SQ(ISYMBJ,ISYMAI) + 1
             CALL SCAL2(X2SQ(IXPOS),X2SQ(IXPOS2),
     &                  NT1AM(ISYMAI),NT1AM(ISYMBJ))
          END DO

      END IF
C
      CONTAINS
C
         PURE SUBROUTINE SCAL1(A,N)
C
         DOUBLE PRECISION, INTENT(INOUT) :: A(N,N)
         INTEGER, INTENT(IN) :: N
         DOUBLE PRECISION :: TMP1, TMP2
         INTEGER I, J
C
         DO I = 1, N
            DO J = 1, I-1
               TMP1 = A(I,J)
               TMP2 = A(J,I)
               A(I,J) = FN*TMP1 - FT*TMP2
               A(J,I) = FN*TMP2 - FT*TMP1
            END DO
         END DO
C
         END SUBROUTINE
C
         PURE SUBROUTINE SCAL2(A,B,M,N)
C
         DOUBLE PRECISION, INTENT(INOUT) :: A(M,N), B(N,M)
         INTEGER, INTENT(IN) :: N, M
         INTEGER :: I, J
         DOUBLE PRECISION TMP1, TMP2
C
         DO I = 1, M
            DO J = 1, N
               TMP1 = A(I,J)
               TMP2 = B(J,I)
               A(I,J) = FN*TMP1 - FT*TMP2
               B(J,I) = FN*TMP2 - FT*TMP1
            END DO
         END DO
C
         END SUBROUTINE

      END SUBROUTINE
C
      SUBROUTINE CCSD_TCMEPK3(T2AM,ISYOPE)
C
C     Swap I and J in the antisymmetric doubles array
C     T2AM(ai,bj) -> T2AM(aj,bi)
C     Since only the lower half is stored,
C     in practice A and B is swapped instead:
C     T2AM(ai,bj) -> - T2AM(bi,aj)
C     only for i != j
C
C     Based on CCSD_TCMEPK by
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C
#include "implicit.h"
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, ONEM = -1.0D0)
C
      DIMENSION T2AM(*)
C
#include "priunit.h"
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CCSD_TCMEPK3')
C
      DO 100 ISYMIJ = 1,NSYM
C
         ISYMAB = MULD2H(ISYMIJ,ISYOPE)
C
         DO 110 ISYMJ = 1,NSYM
C
            ISYMI = MULD2H(ISYMJ,ISYMIJ)
C
            IF (ISYMI .GT. ISYMJ) GOTO 110
C
            DO 120 ISYMB = 1,NSYM
C
               ISYMA = MULD2H(ISYMB,ISYMAB)
C
               IF (ISYMA .GT. ISYMB) GOTO 120
C
               ISYMAI = MULD2H(ISYMA,ISYMI)
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYMBI = MULD2H(ISYMB,ISYMI)
               ISYMAJ = MULD2H(ISYMA,ISYMJ)
C
C
               DO 130 J = 1,NRHF(ISYMJ)
C
                  IF (ISYMI .EQ. ISYMJ) THEN
                     NRHFI = J - 1
                  ELSE
                     NRHFI = NRHF(ISYMI)
                  ENDIF
C
                  DO 140 I = 1,NRHFI
C
                     DO 150 B = 1,NVIR(ISYMB)
C
                        IF (ISYMB .EQ. ISYMA) THEN
                           NVIRA = B
                        ELSE
                           NVIRA = NVIR(ISYMA)
                        ENDIF
C
                        NBI = IT1AM(ISYMB,ISYMI)
     *                      + NVIR(ISYMB)*(I - 1) + B
                        NBJ = IT1AM(ISYMB,ISYMJ)
     *                      + NVIR(ISYMB)*(J - 1) + B
C
                        DO 160 A = 1,NVIRA
C
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I - 1) + A
                           NAJ = IT1AM(ISYMA,ISYMJ)
     *                         + NVIR(ISYMA)*(J - 1) + A
C
                         FACT = ONE
                         IF (ISYMAI.EQ.ISYMBJ) THEN
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                           + INDEX(NAI,NBJ)
                         ELSE IF (ISYMAI.LT.ISYMBJ) THEN
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                           + NT1AM(ISYMAI)*(NBJ-1)+NAI
                         ELSE IF (ISYMBJ.LT.ISYMAI) THEN
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                           + NT1AM(ISYMBJ)*(NAI-1)+NBJ
                           FACT = FACT*ONEM
                         END IF
C
                         IF (ISYMAJ.EQ.ISYMBI) THEN
                           NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                           + INDEX(NAJ,NBI)
                           FACT = FACT*ONEM
                         ELSE IF (ISYMAJ.LT.ISYMBI) THEN
                           NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                           + NT1AM(ISYMAJ)*(NBI-1)+NAJ
                         ELSE IF (ISYMBI.LT.ISYMAJ) THEN
                           NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                           + NT1AM(ISYMBI)*(NAJ-1)+NBI
                           FACT = FACT*ONEM
                         END IF
C
                         XAIBJ = T2AM(NAJBI)
                         XAJBI = T2AM(NAIBJ)
C
                         T2AM(NAJBI) = FACT*XAJBI
                         T2AM(NAIBJ) = FACT*XAIBJ
C
  160                   CONTINUE
  150                CONTINUE
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C---------------------------------------
C     Scale diagonal elements of result.
C---------------------------------------
C
      IF (ISYOPE .NE. 1) THEN
         CALL QEXIT('CCSD_TCMEPK3')
         RETURN
      ENDIF
C
      DO 200 ISYMAI = 1,NSYM
         DO 210 NAI = 1,NT1AM(ISYMAI)
            NAIAI = IT2AM(ISYMAI,ISYMAI) + INDEX(NAI,NAI)
            T2AM(NAIAI) = 0.0D0
  210    CONTINUE
  200 CONTINUE
C
      CALL QEXIT('CCSD_TCMEPK3')
C
      RETURN
      END
C  /* Deck cc_zwvi3 */
      SUBROUTINE CC_ZWVI3(ZINT,CTR2,ISYMC2 ,TINT,ISYTIN,
     *                   WORK,LWORK)
C
C     Written by Asger Halkier 26/10 - 1995
C
C     Version: 1.0
C
C     Purpose: To calculate the intermediates entering some of the
C              terms in the 2.1-block.
C     Ove Christiansen 1-10-1996:
C              general symmetry of ctr2 (isymc2)
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ZINT(*), CTR2(*), TINT(*), WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "cclr.h"
C
C-----------------------------
C     Initialize result array.
C-----------------------------
C
      ISYZIN = MULD2H(ISYMC2,ISYTIN)
C
      CALL DZERO(ZINT,NT2BCD(ISYZIN))
C
C------------------------
C     Do the calculation.
C------------------------
C
      DO 100 ISYMDK = 1,NSYM
C
         ISYMEI = MULD2H(ISYMDK,ISYMC2)
         ISYMJ  = MULD2H(ISYMDK,ISYTIN)
C
         KOFF1  = IT2SQ(ISYMDK,ISYMEI) + 1
         KOFF2  = IT2BCD(ISYMDK,ISYMJ) + 1
         KOFF3  = IT2BCD(ISYMEI,ISYMJ) + 1
C
         NTOTEI = MAX(NT1AM(ISYMEI),1)
         NTOTDK = MAX(NT1AM(ISYMDK),1)
C
         CALL DGEMM('T','N',NT1AM(ISYMEI),NRHF(ISYMJ),NT1AM(ISYMDK),
     *              ONE,CTR2(KOFF1),NTOTDK,TINT(KOFF2),NTOTDK,ZERO,
     *              ZINT(KOFF3),NTOTEI)
C
  100 CONTINUE

      RETURN
      END

C  /* Deck cc_t2ao3 */
      SUBROUTINE CC_T2AO3(T2AM,XLAMDH,ISYMLH,SCRM,WORK,LWORK,
     *                    IDEL,ISYMD,ISYMTR)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Written by Rasmus Faber, 2017
C     Based on CC_T2AO by H. Kock, A. Sanchez, O. Christiansen and
C     A. Halkier.
C     PURPOSE:
C             Tdjci -> Tci,j (delta)
C     i.e. transform the other index than CC_T2AO, needed in case
C     of non-symmetric T2AM
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      implicit none
#include "priunit.h"
#include "iratdef.h"
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, ONE = 1.0D0
      CHARACTER(LEN=*), PARAMETER :: myname = 'CC_T2AO3'
C
      INTEGER, INTENT(IN) :: ISYMLH, LWORK, IDEL, ISYMD, ISYMTR
      DOUBLE PRECISION, INTENT(IN) ::  T2AM(*),XLAMDH(*)
      DOUBLE PRECISION, INTENT(OUT) :: SCRM(*)
      DOUBLE PRECISION, INTENT(INOUT) :: WORK(LWORK)
C
#include "ccorb.h"
#include "ccsdsym.h"
C
      INTEGER :: ISYDVI, ISYMM, ISYMCI, ISYMJ, ISYMDJ
      INTEGER :: KOFF1, KOFF2, KOFF3
      INTEGER :: NTOTCI, NCI, NTOTDVI

C
C-----------------------------------------------------
C     Calculate the transformed t2-amplitude and save.
C-----------------------------------------------------
C
      ISYDVI = MULD2H(ISYMLH,ISYMD)
      ISYMM = MULD2H(ISYMTR,ISYDVI)
      CALL DZERO(SCRM,NT2BCD(ISYMM))
C
      IF ( LWORK .LT. NVIR(ISYDVI)) THEN
         CALL QUIT('Insufficient core in '//myname)
      ENDIF
C
      CALL DZERO(WORK,NVIR(ISYDVI))
C
      KOFF1 = IGLMVI(ISYMD,ISYDVI) + IDEL - IBAS(ISYMD)
C
      CALL DCOPY(NVIR(ISYDVI),XLAMDH(KOFF1),NBAS(ISYMD),WORK,1)
C
      DO 100 ISYMCI = 1,NSYM
C
         ISYMDJ = MULD2H(ISYMTR,ISYMCI)
         ISYMJ  = MULD2H(ISYMDJ,ISYDVI)
C
         NTOTCI = MAX(1,NT1AM(ISYMCI))
         NTOTDVI = MAX(1,NVIR(ISYDVI))
C
         DO 110 NCI = 1, NT1AM(ISYMCI)
C
            KOFF2 = IT2SQ(ISYMDJ,ISYMCI)
     *            + NT1AM(ISYMDJ)*(NCI - 1)
     &            + IT1AM(ISYDVI,ISYMJ) + 1
            KOFF3 = IT2BCD(ISYMCI,ISYMJ)
     *            + NCI
C
            CALL DGEMV('T',NVIR(ISYDVI),NRHF(ISYMJ),ONE,
     *                 T2AM(KOFF2),NTOTDVI,WORK,1,ZERO,
     *                 SCRM(KOFF3),NTOTCI)
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_t2mot */
      SUBROUTINE CC_T2MOT(RHO1,CTR2,ISYMC2,OMEGA2,RHO2,GAMMA,XLAMDP,
     *                    XLAMPC,ISYMPC,WORK,LWORK,ISYMBF,ICON)
C
C     Rasmus Faber
C
C     Based on cc_t2mo by Henrik Koch and Alfredo Sanchez.
C
C     Transform the Omega2 vector from the AO basis to the MO
C     basis. Triplet version.
C
C
      implicit none
C
#include "priunit.h"
#include "maxorb.h"
      CHARACTER(len=*), PARAMETER :: myname = "CC_T2MOT"
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0
      DOUBLE PRECISION, PARAMETER :: FOURTH = 0.25D0
      DOUBLE PRECISION, INTENT(IN):: CTR2(*), OMEGA2(*),
     &                               GAMMA(*), XLAMDP(*), XLAMPC(*)
      DOUBLE PRECISION, INTENT(OUT):: RHO1(*), RHO2(*), WORK(*)
#include "ccorb.h"
#include "ccsdsym.h"
#include "symsq.h"
C
      INTEGER, INTENT(IN) :: ISYMC2, ISYMPC, ISYMBF, ICON, LWORK

      integer :: index
      INTEGER :: ISYMI, ISYMJ, ISYAL, ISYBE, ISALBE, ISYALI, ISYBEJ,
     &           ISYMA, ISYMB, ISYMAI, ISYMBJ, ISYMO1, ISYMO2,
     &           ISYMIJ, ISYMAB
      INTEGER :: NVA, NRA, NVB, NRB, NIJP, NIJM, NABP, NABIJP, NABIJM,
     &           NAB, NAI, NBJ, NBASA, NBASB, NVIRA, NAIBJ
      INTEGER :: KSCR1, KSCR2, KSCR3, KSCR4, KSCR5
      INTEGER :: KEND1, LWRK1
      INTEGER :: KOFF1, KOFF2
      DOUBLE PRECISION :: FACP, FACM
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      ISYMO2 = MULD2H(ISYMBF,ISYMPC)
      ISYMO1 = MULD2H(ISYMO2,ISYMC2)
C
      IF ((ICON.EQ.1).OR.(ICON.EQ.2)) THEN
         CALL DZERO(RHO2,NT2SQ(ISYMO2))
      ENDIF
C
      DO 100 ISYMJ = 1,NSYM
         DO 110 ISYMI = 1,NSYM
C
            ISYMIJ = MULD2H(ISYMI,ISYMJ)
            ISALBE = MULD2H(ISYMIJ,ISYMBF)
            ISYMAB = MULD2H(ISYMIJ,ISYMO2)
C
            DO 120 ISYBE = 1,NSYM
C
               ISYAL  = MULD2H(ISYBE,ISALBE)
               ISYALI = MULD2H(ISYAL,ISYMI)
               ISYBEJ = MULD2H(ISYBE,ISYMJ)
C
C-----------------------------------------------
C              Dynamic allocation of work space.
C-----------------------------------------------
C
               ISYMA = MULD2H(ISYAL,ISYMPC)
               NVA = MAX(NVIR(ISYMA),NVIR(ISYAL))
               NRA = MAX(NRHF(ISYMA),NRHF(ISYAL))
               ISYMB = MULD2H(ISYBE,ISYMPC)
               NVB = MAX(NVIR(ISYMB),NVIR(ISYBE),NRHF(ISYBE))
               NRB = MAX(NRHF(ISYMB),NRHF(ISYBE))
C
               KSCR1 = 1
               IF (ICON.NE.4) THEN
                  KSCR2 = KSCR1 + NBAS(ISYAL)*NBAS(ISYBE)
                  KSCR3 = KSCR2 + NBAS(ISYAL)*NVB
                  KEND1 = KSCR3 + NVA*NVB
               ELSE
                  KSCR4 = KSCR1 + NBAS(ISYAL)*NBAS(ISYBE)
                  KEND1 = KSCR4 + NBAS(ISYAL)*NRB
               END IF
               LWRK1 = LWORK - KEND1
C
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Not enough space in '//myname)
               END IF
C
               DO 130 J = 1,NRHF(ISYMJ)
                  DO 140 I = 1,NRHF(ISYMI)
C
C------------------------------------------
C                    Squareup the AB block.
C------------------------------------------
C                    Plus version is stored for i<j with an implicit
C                    symmetry i<j = - i>j
C
                     IF (ISYMI .EQ. ISYMJ) THEN
                        NIJP = IMIJP(ISYMI,ISYMJ) + INDEX(I,J)
                        FACP = ONE
                        IF (I .EQ. J) FACP = ZERO
                        IF (I .GT. J) FACP = -ONE
                     ELSE IF (ISYMI .LT. ISYMJ) THEN
                        NIJP = IMIJP(ISYMI,ISYMJ)
     *                       + NRHF(ISYMI)*(J - 1) + I
                        FACP = ONE
                     ELSE
                        NIJP = IMIJP(ISYMJ,ISYMI)
     *                       + NRHF(ISYMJ)*(I - 1) + J
                        FACP = -ONE
                     ENDIF
C                    Minus version: i,j stored as square array
                     NIJM = IMATIJ(ISYMI,ISYMJ)
     &                    + NRHF(ISYMI)*(J-1) + I
C
                     DO 166 B = 1, NBAS(ISYBE)
                        DO 167 A = 1, NBAS(ISYAL)
C
                           IF (ISYAL .EQ. ISYBE) THEN
                              NABP = IAODPK(ISYAL,ISYBE)
     *                             + INDEX(A,B)
                              FACM = ONE
                              IF (A .EQ. B) FACM = ZERO
                              IF (A .GT. B) FACM = -ONE
                           ELSE IF (ISYAL .LT. ISYBE) THEN
                              NABP = IAODPK(ISYAL,ISYBE)
     *                             + NBAS(ISYAL)*(B - 1) + A
                              FACM = ONE
                           ELSE
                              NABP = IAODPK(ISYBE,ISYAL)
     *                             + NBAS(ISYBE)*(A - 1) + B
                              FACM = -ONE
                           ENDIF
C
                           NABIJP = IT2ORT(ISALBE,ISYMIJ)
     *                            + NNBST(ISALBE)*(NIJP - 1) + NABP
C
                           NABIJM = NT2ORT(ISYMBF)
     *                            + IT2ORT3(ISALBE,ISYMIJ)
     *                            + NNBST(ISALBE)*(NIJM - 1) + NABP
C
                           NAB   = KSCR1 + NBAS(ISYAL)*(B - 1) + A - 1
C
                           WORK(NAB) = FOURTH*
     &                       (FACP*OMEGA2(NABIJP) +
     &                        FACM*OMEGA2(NABIJM))
C
  167                   CONTINUE
  166                CONTINUE
C
C------------------------------------------------------------
C                    Transform the AB block to virtual space.
C------------------------------------------------------------
C
                     IF ((ICON.EQ.1).OR.(ICON.EQ.2)) THEN
C
                     ISYMA = MULD2H(ISYAL,ISYMPC)
                     ISYMB = ISYBE
                     ISYMAI = MULD2H(ISYMA,ISYMI)
                     ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
                     NBASA = MAX(NBAS(ISYAL),1)
                     NBASB = MAX(NBAS(ISYBE),1)
                     NVIRA = MAX(NVIR(ISYMA),1)
C
                     KOFF1 = ILMVIR(ISYBE) + 1
C
                     CALL DGEMM('N','N',NBAS(ISYAL),NVIR(ISYMB),
     *                          NBAS(ISYBE),ONE,WORK(KSCR1),NBASA,
     *                          XLAMDP(KOFF1),NBASB,ZERO,WORK(KSCR2),
     *                          NBASA)
C
                     KOFF2 = IGLMVI(ISYAL,ISYMA) + 1
C
                     CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),
     *                          NBAS(ISYAL),ONE,XLAMPC(KOFF2),NBASA,
     *                          WORK(KSCR2),NBASA,ZERO,WORK(KSCR3),
     *                          NVIRA)
C
C--------------------------------------------
C                    Store the omega2 vector.
C--------------------------------------------
C
                     DO 170 B = 1,NVIR(ISYMB)
                        NBJ   = IT1AM(ISYMB,ISYMJ)
     *                        + NVIR(ISYMB)*(J-1) + B
                        DO 180 A = 1,NVIR(ISYMA)
C
                           NAI   = IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1) + A
                           NAB   = KSCR3 + NVIR(ISYMA)*(B - 1) + A - 1
C
                           NAIBJ = IT2SQ(ISYMAI,ISYMBJ) +
     &                             NT1AM(ISYMAI)*(NBJ-1) + NAI
C
                           RHO2(NAIBJ) = RHO2(NAIBJ) + WORK(NAB)
C
  180                   CONTINUE
  170                CONTINUE
C
                     ENDIF
C
C--------------------------------------
C                    CCLR contribution.
C--------------------------------------
C
                     IF (ICON .EQ. 2 ) THEN
C
                        CALL DZERO(WORK(KSCR2),NVA*NVB)
                        ISYMA = ISYAL
                        ISYMB = MULD2H(ISYBE,ISYMPC)
                        ISYMAI = MULD2H(ISYMA,ISYMI)
                        ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
                        NBASA = MAX(NBAS(ISYAL),1)
                        NBASB = MAX(NBAS(ISYBE),1)
                        NVIRA = MAX(NVIR(ISYMA),1)
C
                        KOFF1 = IGLMVI(ISYBE,ISYMB) + 1
C
                        CALL DGEMM('N','N',NBAS(ISYAL),NVIR(ISYMB),
     *                             NBAS(ISYBE),ONE,WORK(KSCR1),NBASA,
     *                             XLAMPC(KOFF1),NBASB,ZERO,WORK(KSCR2),
     *                             NBASA)
C
                        KOFF2 = ILMVIR(ISYAL) + 1
C
                        CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),
     *                             NBAS(ISYAL),ONE,XLAMDP(KOFF2),NBASA,
     *                             WORK(KSCR2),NBASA,ZERO,WORK(KSCR3),
     *                             NVIRA)
C
C--------------------------------------------
C                    Store the omega2 vector.
C--------------------------------------------
                        DO 171 B = 1,NVIR(ISYMB)
                           NBJ   = IT1AM(ISYMB,ISYMJ)
     *                           + NVIR(ISYMB)*(J-1) + B
                           DO 181 A = 1,NVIR(ISYMA)
C
                              NAI   = IT1AM(ISYMA,ISYMI)
     *                              + NVIR(ISYMA)*(I-1) + A
                              NAB  = KSCR3 + NVIR(ISYMA)*(B - 1) + A - 1
C
                              NAIBJ = IT2SQ(ISYMAI,ISYMBJ) +
     &                                NT1AM(ISYMAI)*(NBJ-1) + NAI
C
                              RHO2(NAIBJ) = RHO2(NAIBJ) + WORK(NAB)
C
  181                      CONTINUE
  171                   CONTINUE
C
C
                     ENDIF
C
C
  999                CONTINUE
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END

      SUBROUTINE CC_TRIP_EXTRACT(X2SQ,X2AMP,ISYMTR,ANTI)
C
C     Rasmus Faber, Aug 2017
C
C     Extracts the components of the "+" triplet doubles vector from the
C     square array "X2SQ" and puts it in X2AMP ordered AI<BJ
C     If anti is true, extract (-) component instead.
C
      IMPLICIT NONE
#include "ccorb.h"
#include "ccsdsym.h"
C
      DOUBLE PRECISION, INTENT(IN) :: X2SQ(*)
      DOUBLE PRECISION, INTENT(OUT) :: X2AMP(*)
      INTEGER, INTENT(IN) :: ISYMTR
      LOGICAL, INTENT(IN) :: ANTI
C
      INTEGER :: ISYMAI, ISYMBJ
      INTEGER :: NAI, NBJ, NAIBJS, NAIBJT, NBJAIS
      DOUBLE PRECISION :: FACTOR, XAIBJ, XBJAI
C
      FACTOR = 1.0D0
      IF (ANTI) FACTOR = -1.0D0
C
      DO ISYMBJ = 1, NSYM
         ISYMAI = MULD2H(ISYMBJ,ISYMTR)

         IF (ISYMTR.EQ.1) THEN
            DO NBJ = 1, NT1AM(ISYMBJ)
               DO NAI = 1, NBJ
C
                  NAIBJS = IT2SQ(ISYMAI,ISYMBJ)
     &                   + NT1AM(ISYMAI)*(NBJ-1) + NAI
                  NBJAIS = IT2SQ(ISYMBJ,ISYMAI)
     %                   + NT1AM(ISYMBJ)*(NAI-1) + NBJ
                  NAIBJT = IT2AM(ISYMAI,ISYMBJ)
     &                   + NBJ*(NBJ-1)/2 + NAI
C
                  XAIBJ  = X2SQ(NAIBJS)
                  XBJAI  = X2SQ(NBJAIS)
                  X2AMP(NAIBJT) = XAIBJ + FACTOR*XBJAI
               END DO
            END DO
         ELSE IF ( ISYMAI .LT. ISYMBJ) THEN
            DO NBJ = 1, NT1AM(ISYMBJ)
               DO NAI = 1, NT1AM(ISYMAI)
                  NAIBJS = IT2SQ(ISYMAI,ISYMBJ)
     &                   + NT1AM(ISYMAI)*(NBJ-1) + NAI
                  NBJAIS = IT2SQ(ISYMBJ,ISYMAI)
     %                   + NT1AM(ISYMBJ)*(NAI-1) + NBJ
                  NAIBJT = IT2AM(ISYMAI,ISYMBJ)
     &                   + NT1AM(ISYMAI)*(NBJ-1) + NAI
C
                  XAIBJ  = X2SQ(NAIBJS)
                  XBJAI  = X2SQ(NBJAIS)
                  X2AMP(NAIBJT) = XAIBJ + FACTOR*XBJAI
               END DO
            END DO
         END IF

      END DO
C
      END SUBROUTINE
C
C  /* Deck cc_22am3 */
      SUBROUTINE CC_22AM3(RHO2,XIAJB,XMINT,ISYMIM,WORK,LWORK)
C
C     Written by Rasmus Faber - 2017
C     Based on cc_22am by Asger Halkier
C
C     Purpose: Contract the M intermediate with the (ma|nb) integrals
C     to obtain the E contribution to the left transformed doubles
C     vector.
C
C
      IMPLICIT NONE
C
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, ONE = 1.0D0,
     &                               HALF = 0.5D0
      DOUBLE PRECISION :: RHO2(*), XIAJB(*), XMINT(*), WORK(LWORK)
      INTENT(IN) :: XIAJB
      INTENT(INOUT) :: XMINT
      INTENT(OUT) :: RHO2, WORK
      INTEGER, INTENT(IN) :: ISYMIM, LWORK
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "cclr.h"
C
      INTEGER :: INDEX
      INTEGER :: ISYMJ, ISYMB, ISYMA, ISYMI, ISYMM, ISYMN,
     &           ISYRES, ISYMIN, ISYMBJ, ISYMAI,
     &           ISYMAM, ISYMAN, ISYMMI, ISYMBN
      INTEGER :: KSCR1, KSCR2
      INTEGER :: LWRK1, LWRK2
      INTEGER :: KEND1, KEND2
      INTEGER :: KOFF1, KOFF2, KOFF3
      INTEGER :: NBJ, NBN, NAM, NAI, NAMBN, NAIBJ
      INTEGER :: NTOTA, NTOTM
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      ISYRES = ISYMIM
C
      CALL TRANSP_MI(XMINT,ISYRES,WORK(1))
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMIN = MULD2H(ISYMJ,ISYRES)
C
         IF (NRHF(ISYMJ) .EQ. 0) GOTO 100
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            DO 120 ISYMB = 1,NSYM
C
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYMAI = MULD2H(ISYMBJ,ISYRES)
C
               IF (NVIR(ISYMB) .EQ. 0) GOTO 120
C
C----------------------------------------
C              Work space allocation one.
C----------------------------------------
C
               KSCR1 = 1
               KEND1 = KSCR1 + NT1AM(ISYMAI)
               LWRK1 = LWORK - KEND1
C
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient work space in CC_22AM3')
               ENDIF
C
               DO 130 B = 1,NVIR(ISYMB)
C
                  CALL DZERO(WORK(KSCR1),NT1AM(ISYMAI))
C
                  NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
                  DO 140 ISYMN = 1,NSYM
C
                     ISYMBN = MULD2H(ISYMN,ISYMB)
                     ISYMAM = MULD2H(ISYMBN,ISYMOP)
                     ISYMMI = MULD2H(ISYMN,ISYMIN)
C
C----------------------------------------------
C                    Work space allocation two.
C----------------------------------------------
C
                     KSCR2 = KEND1
                     KEND2 = KSCR2 + NT1AM(ISYMAM)
                     LWRK2 = LWORK - KEND2
C
                     IF (LWRK2 .LT. 0) THEN
                        CALL QUIT('Insufficient work space in CC_22AM3')
                     ENDIF
C
                     DO 150 N = 1,NRHF(ISYMN)
C
                        NBN = IT1AM(ISYMB,ISYMN)
     &                      + NVIR(ISYMB)*(N - 1) + B
C
C-------------------------------------------------------
C                       Copy submatrix out of integrals.
C-------------------------------------------------------
C
                        DO 160 NAM = 1,NT1AM(ISYMAM)
C
                           NAMBN = IT2AM(ISYMAM,ISYMBN)
     &                           + INDEX(NAM,NBN)
C
                           WORK(KSCR2 + NAM - 1) = XIAJB(NAMBN)
C
  160                   CONTINUE
C
C--------------------------------------------------------------------
C                       Contraction of integrals with M-intermediate.
C--------------------------------------------------------------------
C
                        DO 170 ISYMA = 1,NSYM
C
                           ISYMM = MULD2H(ISYMA,ISYMAM)
                           ISYMI = MULD2H(ISYMA,ISYMAI)
C
                           KOFF1 = KSCR2 + IT1AM(ISYMA,ISYMM)
                           KOFF2 = I3ORHF(ISYMIN,ISYMJ)
     &                           + NMAIJK(ISYMIN)*(J - 1)
     &                           + IMAIJK(ISYMMI,ISYMN)
     &                           + NMATIJ(ISYMMI)*(N - 1)
     &                           + IMATIJ(ISYMM,ISYMI) + 1
                           KOFF3 = KSCR1 + IT1AM(ISYMA,ISYMI)
C
                           NTOTA = MAX(NVIR(ISYMA),1)
                           NTOTM = MAX(NRHF(ISYMM),1)
C
                           CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),
     &                                NRHF(ISYMM),ONE,WORK(KOFF1),
     &                                NTOTA,XMINT(KOFF2),NTOTM,ONE,
     &                                WORK(KOFF3),NTOTA)
C
  170                   CONTINUE
C
  150                CONTINUE
  140             CONTINUE
C
C-----------------------------------------------
C                 Storage in RHO2 result vector.
C-----------------------------------------------
C
                  DO NAI = 1, NT1AM(ISYMAI)
                     NAIBJ = IT2SQ(ISYMAI,ISYMBJ)
     &                     + NT1AM(ISYMAI)*(NBJ-1) + NAI
                     RHO2(NAIBJ) = RHO2(NAIBJ)
     &                           + HALF*WORK(KSCR1+NAI-1)
                  END DO
C
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
C
      CONTAINS
         SUBROUTINE TRANSP_MI(MINT,ISYMINT,WORK)
C
C        Transposes the M intermediate:
C        M(M,I,N,J) -> M(N,J,M,I)
C
         DOUBLE PRECISION, INTENT(INOUT) :: MINT(*), WORK(*)
         INTEGER, INTENT(IN) :: ISYMINT

         INTEGER :: I, J, N
         INTEGER :: ISYMI, ISYMJ, ISYMM, ISYMN, ISYMIN, ISYNJM,
     &              ISYMMI, ISYMNJ, ISYMMN
         INTEGER :: KOFF1, KOFF2
         INTEGER :: NTOTI, NTOTM

         DO ISYMJ = 1, NSYM
            ISYMIN = MULD2H(ISYMJ,ISYMINT)
            DO ISYMI = 1, ISYMJ
               ISYMMN = MULD2H(ISYMI,ISYMIN)
               ISYNJM = MULD2H(ISYMI,ISYMINT)
               NTOTI = NRHF(ISYMI)
               DO J = 1, NRHF(ISYMJ)
                  IF (ISYMI.EQ.ISYMJ) NTOTI = J - 1
                  DO I = 1, NTOTI
                     DO ISYMN = 1, NSYM
                        ISYMM = MULD2H(ISYMN,ISYMMN)
                        ISYMMI = MULD2H(ISYMIN,ISYMN)
                        ISYMNJ = MULD2H(ISYNJM,ISYMM)
                        DO N = 1, NRHF(ISYMN)
                           KOFF1 = I3ORHF(ISYMIN,ISYMJ)
     &                           + NMAIJK(ISYMIN)*(J-1)
     &                           + IMAIJK(ISYMMI,ISYMN)
     &                           + NMATIJ(ISYMMI)*(N-1)
     &                           + IMATIJ(ISYMM,ISYMI)
     &                           + NRHF(ISYMM)*(I-1) + 1

                           KOFF2 = I3ORHF(ISYNJM,ISYMI)
     &                           + NMAIJK(ISYNJM)*(I-1)
     &                           + IMAIJK(ISYMNJ,ISYMM)
     &                           + IMATIJ(ISYMN,ISYMJ)
     &                           + NRHF(ISYMN)*(J-1) + N

                           CALL DCOPY(NRHF(ISYMM),
     &                                MINT(KOFF2),NMATIJ(ISYMNJ),
     &                                WORK(1),1)

                           CALL DCOPY(NRHF(ISYMM),
     &                                MINT(KOFF1),1,
     &                                MINT(KOFF2),NMATIJ(ISYMNJ))

                           CALL DCOPY(NRHF(ISYMM),
     &                                WORK(1),1,
     &                                MINT(KOFF1),1)
                        END DO ! N
                     END DO ! ISYMN
                  END DO ! I
!
!       Handle I = J case
                  IF (ISYMI.EQ.ISYMJ) THEN
                     I = J
                     DO ISYMN = 1, NSYM
                        ISYMM = MULD2H(ISYMMN,ISYMN)
                        IF (ISYMM.GT.ISYMN) CYCLE
                        NTOTM = NRHF(ISYMM)
                        ISYMMI = MULD2H(ISYMIN,ISYMN)
                        ISYMNJ = MULD2H(ISYNJM,ISYMM)
                        DO N = 1, NRHF(ISYMN)
                           IF (ISYMM.EQ.ISYMN) NTOTM = N - 1
                           KOFF1 = I3ORHF(ISYMIN,ISYMJ)
     &                           + NMAIJK(ISYMIN)*(J-1)
     &                           + IMAIJK(ISYMMI,ISYMN)
     &                           + NMATIJ(ISYMMI)*(N-1)
     &                           + IMATIJ(ISYMM,ISYMI)
     &                           + NRHF(ISYMM)*(I-1) + 1

                           KOFF2 = I3ORHF(ISYNJM,ISYMI)
     &                           + NMAIJK(ISYNJM)*(I-1)
     &                           + IMAIJK(ISYMNJ,ISYMM)
     &                           + IMATIJ(ISYMN,ISYMJ)
     &                           + NRHF(ISYMN)*(J-1) + N

                           CALL DCOPY(NTOTM,
     &                                MINT(KOFF2),NMATIJ(ISYMNJ),
     &                                WORK(1),1)

                           CALL DCOPY(NTOTM,
     &                                MINT(KOFF1),1,
     &                                MINT(KOFF2),NMATIJ(ISYMNJ))

                           CALL DCOPY(NTOTM,
     &                                WORK(1),1,
     &                                MINT(KOFF1),1)
                        END DO
                     END DO
                  END IF
               END DO ! J
            END DO ! ISYMI
         END DO ! ISYMJ

         END SUBROUTINE
      END
C
      SUBROUTINE CCRHS_ASQ(OMEGA2,T2AM,GAMMA,WORK,LWORK,ISYGAM,ISYVEC,
     *                     IOPT)
C
C     Rasmus Faber 2017
C     Based on CCRHS_A by Henrik Koch & Ove Christiansen
C
C     Purpose: Calculates the doubles A term
C     Contract the "Gamma" intermediate with the quadratic
C     trial-vector "T2AM".
C     Use
C        IOPT = 1 for right transformation
C        IOPT = 2 for left transformation
C
C     This version uses a quadratic storage on the output vector, OMEGA2
C
      IMPLICIT NONE
C
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D0, ONE=1.0D0, HALF=0.5D0
      CHARACTER(LEN=*), PARAMETER :: myname = "CCRHS_ASQ"
C
      DOUBLE PRECISION, INTENT(INOUT) :: OMEGA2(*)
      DOUBLE PRECISION, INTENT(IN) :: T2AM(*), GAMMA(*)
      DOUBLE PRECISION, INTENT(OUT) :: WORK(LWORK)
      INTEGER, INTENT(IN) :: LWORK, ISYGAM, ISYVEC, IOPT
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
      INTEGER :: INDEX
      INTEGER :: ISYMA, ISYMB, ISYMI, ISYMJ, ISYMK, ISYML,
     &           ISYMAI, ISYMAK,  ISYMBJ, ISYMBL, ISYMKI, ISYMLJ, ISAIBJ
      INTEGER :: NAI, NBJ, NBL, NLJ, NKI, NAIBJ, NKILJ, NSTO
      INTEGER :: NRHFK, NVIRA
      INTEGER :: KSCR1, KSCR2
      INTEGER :: KEND1, KEND2
      INTEGER :: LWRK1, LWRK2
      INTEGER :: KOFF1, KOFF2, KOFF3
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
C----------------------------
C     Calculate contribution.
C----------------------------
C
      ISAIBJ = MULD2H(ISYGAM,ISYVEC)
C
      DO 100 ISYMLJ = 1,NSYM
C
         ISYMKI = MULD2H(ISYMLJ,ISYGAM)
C
         KSCR1 = 1
         KEND1 = KSCR1 + NMATIJ(ISYMKI)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient space for allocation in '//myname)
         END IF
C
         DO 110 ISYMJ = 1,NSYM
C
            ISYML = MULD2H(ISYMJ,ISYMLJ)
C
            DO 120 J = 1,NRHF(ISYMJ)
C
               DO 130 L = 1,NRHF(ISYML)
C
                  IF (IOPT .EQ. 1) THEN
C
                     NLJ = IMATIJ(ISYML,ISYMJ)
     *                   + NRHF(ISYML)*(J - 1) + L
C
                  ELSE IF (IOPT .EQ. 2) THEN
C
                     NLJ = IMATIJ(ISYMJ,ISYML)
     *                   + NRHF(ISYMJ)*(L - 1) + J
C
                  ENDIF
C
                  DO 140 ISYMK = 1,NSYM
C
                     ISYMI = MULD2H(ISYMK,ISYMKI)
C
                     DO 150 I = 1,NRHF(ISYMI)
C
                        DO 160 K = 1,NRHF(ISYMK)
C
                           IF (IOPT .EQ. 1) THEN
C
                              NKI = IMATIJ(ISYMK,ISYMI)
     *                            + NRHF(ISYMK)*(I - 1) + K
C
                           ELSE IF (IOPT .EQ. 2) THEN
C
                              NKI = IMATIJ(ISYMI,ISYMK)
     *                            + NRHF(ISYMI)*(K - 1) + I
C
                           ENDIF
C
                           IF (ISYMKI .EQ. ISYMLJ) THEN
                              NKILJ = IGAMMA(ISYMKI,ISYMLJ)
     *                              + INDEX(NKI,NLJ)
                           ELSE
                              IF (ISYMKI .LT. ISYMLJ) THEN
                                 NKILJ = IGAMMA(ISYMKI,ISYMLJ)
     *                                 + NMATIJ(ISYMKI)*(NLJ - 1) + NKI
                              ELSE
                                 NKILJ = IGAMMA(ISYMLJ,ISYMKI)
     *                                 + NMATIJ(ISYMLJ)*(NKI - 1) + NLJ
                              ENDIF
                           ENDIF
C
                           NSTO = IMATIJ(ISYMK,ISYMI)
     *                          + NRHF(ISYMK)*(I - 1) + K
C
                           WORK(KSCR1 + NSTO - 1) = GAMMA(NKILJ)
C
  160                   CONTINUE
  150                CONTINUE
  140             CONTINUE
C
                  DO 170 ISYMB = 1,NSYM
C
                     ISYMBJ = MULD2H(ISYMB,ISYMJ)
                     ISYMAI = MULD2H(ISYMBJ,ISAIBJ)
                     ISYMBL = MULD2H(ISYMB,ISYML)
                     ISYMAK = MULD2H(ISYVEC,ISYMBL)
C
                     KSCR2 = KEND1
                     KEND2 = KSCR2 + NT1AM(ISYMAI)
                     LWRK2 = LWORK - KEND2
C
                     IF (LWRK2 .LT. 0) THEN
                        CALL QUIT('Insufficient space in '//myname)
                     END IF
C
C                     IF (ISYMAI .GT. ISYMBJ) GOTO 170
C
                     DO 180 B = 1,NVIR(ISYMB)
C
                        NBJ = IT1AM(ISYMB,ISYMJ)
     *                      + NVIR(ISYMB)*(J - 1) + B
                        NBL = IT1AM(ISYMB,ISYML)
     *                      + NVIR(ISYMB)*(L - 1) + B
C
                        CALL DZERO(WORK(KSCR2),NT1AM(ISYMAI))
C
                        DO 190 ISYMI = 1,NSYM
C
                           ISYMK = MULD2H(ISYMI,ISYMKI)
                           ISYMA = MULD2H(ISYMK,ISYMAK)
C
                           NVIRA = MAX(NVIR(ISYMA),1)
                           NRHFK = MAX(NRHF(ISYMK),1)
C
                           KOFF1 = IT2SQ(ISYMAK,ISYMBL)
     *                           + NT1AM(ISYMAK)*(NBL - 1)
     *                           + IT1AM(ISYMA,ISYMK) + 1
                           KOFF2 = KSCR1 + IMATIJ(ISYMK,ISYMI)
                           KOFF3 = KSCR2 + IT1AM(ISYMA,ISYMI)
C
                           CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),
     *                                NRHF(ISYMK),ONE,T2AM(KOFF1),
     *                                NVIRA,WORK(KOFF2),NRHFK,ZERO,
     *                                WORK(KOFF3),NVIRA)
C
  190                   CONTINUE
C
                        DO 200 NAI = 1, NT1AM(ISYMAI)
C
                           NAIBJ = IT2SQ(ISYMAI,ISYMBJ)
     &                           + NT1AM(ISYMAI)*(NBJ-1) + NAI
C
                           OMEGA2(NAIBJ) = OMEGA2(NAIBJ)
     *                                   + HALF*WORK(KSCR2 + NAI - 1)
C
  200                   CONTINUE
C
  180                CONTINUE
  170             CONTINUE
C
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_22e3 */
      SUBROUTINE CC_22E3(RHO2,T2AM,XMAT,YMAT,ISYMXY,WORK,LWORK)
C
C     Rasmus Faber - 2017
C     Based on CC_22EC by Asger Halkier 21/11 - 1995
C
C     Purpose: To calculate the triplet doubles F-contribution
C              to RHO2.
C
C
C
      IMPLICIT NONE
C
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, ONE = 1.0D0,
     &                               ONEM = -ONE
      DOUBLE PRECISION, INTENT(INOUT) :: RHO2(*)
      DOUBLE PRECISION, INTENT(IN) :: T2AM(*), XMAT(*), YMAT(*)
      DOUBLE PRECISION, INTENT(OUT) :: WORK(LWORK)
      INTEGER, INTENT(IN) :: LWORK, ISYMXY
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "cclr.h"
C
      INTEGER :: INDEX
      INTEGER :: ISYMB, ISYMJ, ISYMAI, ISYMAN, ISYMAT, ISYMAX, ISYMAY,
     &           ISYMBJ, ISYMFI, ISYMFY, ISYMIX, ISYMIY, ISYMNX, ISYRES
      INTEGER :: KOFF1, KOFF2, KOFF3, KOFF4, KOFF5, KOFF6
      INTEGER :: NBJ, NAN, NANBJ
      INTEGER :: NTOTA, NTOTF, NTOTN
      INTEGER :: KSCR1, KSCR2, KEND1, LWRK1
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      ISYRES = MULD2H(ISYMOP,ISYMXY)
      ISYMAT = ISYRES
C
      DO 100 ISYMJ = 1,NSYM
C
         IF (NRHF(ISYMJ) .EQ. 0) GOTO 100
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            DO 120 ISYMB = 1,NSYM
C
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYMAI = MULD2H(ISYMBJ,ISYRES)
               ISYMAN = MULD2H(ISYMBJ,ISYMOP)
               ISYMFI = ISYMAN
C
               IF (NVIR(ISYMB) .EQ. 0) GOTO 120
C
C------------------------------------
C              Work space allocation.
C------------------------------------
C
               KSCR1 = 1
               KSCR2 = KSCR1 + NT1AM(ISYMAI)
               KEND1 = KSCR2 + NT1AM(ISYMAN)
               LWRK1 = LWORK - KEND1
C
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient work space in CC_22E3')
               ENDIF
C
               DO 130 B = 1,NVIR(ISYMB)
C
                  CALL DZERO(WORK(KSCR1),NT1AM(ISYMAI))
C
                  NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
C-------------------------------------------------
C                 Copy submatrix out of integrals.
C-------------------------------------------------
C
                  DO 140 NAN = 1,NT1AM(ISYMAN)
C
                     NANBJ = IT2AM(ISYMAN,ISYMBJ) + INDEX(NAN,NBJ)
C
                     WORK(KSCR2 + NAN - 1) = T2AM(NANBJ)
C
  140             CONTINUE
C
C-----------------------------------------------------------------
C                 Contraction of integrals with X- and Y-matrices.
C-----------------------------------------------------------------
C
                  DO 150 ISYMAX = 1,NSYM
C
                     ISYMNX = MULD2H(ISYMAX,ISYMAN)
                     ISYMIX = MULD2H(ISYMNX,ISYMAT)
                     ISYMFY = ISYMAX
                     ISYMIY = MULD2H(ISYMFY,ISYMFI)
                     ISYMAY = MULD2H(ISYMFY,ISYMAT)
C
                     KOFF1 = KSCR2 + IT1AM(ISYMAX,ISYMNX)
                     KOFF2 = IMATIJ(ISYMNX,ISYMIX) + 1
                     KOFF3 = KSCR1 + IT1AM(ISYMAX,ISYMIX)
C
                     NTOTA = MAX(NVIR(ISYMAX),1)
                     NTOTN = MAX(NRHF(ISYMNX),1)
C
                     CALL DGEMM('N','N',NVIR(ISYMAX),NRHF(ISYMIX),
     &                          NRHF(ISYMNX),ONE,WORK(KOFF1),NTOTA,
     &                          XMAT(KOFF2),NTOTN,ONE,WORK(KOFF3),
     &                          NTOTA)
C
                     KOFF4 = IMATAB(ISYMFY,ISYMAY) + 1
                     KOFF5 = KSCR2 + IT1AM(ISYMFY,ISYMIY)
                     KOFF6 = KSCR1 + IT1AM(ISYMAY,ISYMIY)
C
                     NTOTF = MAX(NVIR(ISYMFY),1)
                     NTOTA = MAX(NVIR(ISYMAY),1)
C
                     CALL DGEMM('T','N',NVIR(ISYMAY),NRHF(ISYMIY),
     &                          NVIR(ISYMFY),ONE,YMAT(KOFF4),NTOTF,
     &                          WORK(KOFF5),NTOTF,ONE,WORK(KOFF6),
     &                          NTOTA)
C
  150             CONTINUE
C
C-----------------------------------------------
C                 Storage in RHO2 result vector.
C-----------------------------------------------
C                 Try to rearrange the above DGEMMs to
C                 do it there directly
                  KOFF1 = IT2SQ(ISYMAI,ISYMBJ)
     &                  + NT1AM(ISYMAI)*(NBJ-1) + 1
                  CALL DAXPY(NT1AM(ISYMAI),ONEM,WORK(KSCR1),1,
     &                       RHO2(KOFF1),1)

  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_22cd3 */
      SUBROUTINE CC_22CD3(RHO2,CTR2,ISYVEC,XLAMDH,WORK,LWORK,
     *                    ISYDIM,LUF,FIL,IV,IOPT,FACT)
C
C     Written by Rasmus Faber, 2017
C     Based on CC_22D by Asger Halkier
C
C     Purpose: To calculate the triplet 22C or 22D contribution to RHO2.
C              If 22D: CTR2 needs to be transposed before this routine.
C              If 22C: RHO2 must be transposed.
C
C                 if iopt = 1 the D intermediate is assumed
C                 to be as in energy calc.
C
C                 if iopt ne. 1 we use the intermediate
C                 on lud with address given according to
C                 transformed vector nr iv (iv is not 1
C                 if several vectors are transformed
C                 simultaneously.)
C
C                 in Lambda vector calc: iv=1,iopt=1
C
      IMPLICIT NONE
      CHARACTER(LEN=*), PARAMETER :: myname = 'CC_22CD3'
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, ONE = 1.0D0
      DOUBLE PRECISION, INTENT(INOUT) ::  RHO2(*),  WORK(LWORK)
      DOUBLE PRECISION, INTENT(IN) :: CTR2(*), XLAMDH(*), FACT
      CHARACTER, INTENT(IN) :: FIL*(*)
      INTEGER, INTENT(IN) :: ISYVEC, LWORK, ISYDIM, LUF, IV, IOPT
C
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "cclr.h"
#include "maxorb.h"
#include "ccsdio.h"
C
      INTEGER :: ISYMB, ISYMD, ISYMJ, ISYMAI, ISYMBJ, ISYMEM,
     &           ISYAIJ, ISYJEM, ISYRES
      INTEGER :: IDEL, ID, IOFF, IOFFD, LEN
      INTEGER :: NTOTD, NTOTJ, NTOTAI, NTOTEM
      INTEGER :: KOFF1, KOFF2, KOFF3, KOFF4, KOFF5
      INTEGER :: KPINT, KSCRN, KLAMD, KEND1
      INTEGER :: LWRK1


C
      ISYRES = MULD2H(ISYVEC,ISYDIM)
C
      DO 100 ISYMD = 1,NSYM
C
         ISYMB  = ISYMD
         ISYJEM = MULD2H(ISYMB,ISYDIM)
         ISYAIJ = MULD2H(ISYMB,ISYRES)
C
C------------------------------
C        Work space allocation.
C------------------------------
C
         KPINT = 1
         KSCRN = KPINT + NT2BCD(ISYJEM)
         KLAMD = KSCRN + NT2BCD(ISYAIJ)
         KEND1 = KLAMD + NVIR(ISYMB)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient work space in '//myname)
         ENDIF
C
         LEN = NT2BCD(ISYJEM)
         IF (LEN .LE. 0) CYCLE ! Nothing to do for this sym

         DO 110 IDEL = 1,NBAS(ISYMD)
C
C-----------------------------------------
C           Copy row "IDEL" out of XLAMDH.
C-----------------------------------------
C
            IOFFD = ILMVIR(ISYMB) + IDEL
            NTOTD = MAX(NBAS(ISYMD),1)
            CALL DCOPY(NVIR(ISYMB),XLAMDH(IOFFD),NTOTD,
     &                 WORK(KLAMD),1)

C
C---------------------------------------------------
C           Read P-intermediate submatrix from disc.
C---------------------------------------------------
C
            ID = IDEL + IBAS(ISYMD)
C
            IF (IOPT .EQ. 1) THEN
               IOFF = IT2DEL(ID) + 1
            ELSE
               IOFF = IT2DLR(ID,IV) + 1
            ENDIF
C
            CALL GETWA2(LUF,FIL,WORK(KPINT),IOFF,LEN)
C
            DO 120 ISYMJ = 1,NSYM
C
C--------------------------------------------------------
C              Contraction of intermediate and trial vector.
C--------------------------------------------------------
C
               ISYMAI = MULD2H(ISYMJ,ISYAIJ)
               ISYMEM = MULD2H(ISYMAI,ISYVEC)
               ISYMBJ = MULD2H(ISYMJ,ISYMB)
C
               KOFF1  = IT2SQ(ISYMEM,ISYMAI) + 1
               KOFF2  = KPINT + IT2BCT(ISYMJ,ISYMEM)
               KOFF3  = KSCRN + IT2BCD(ISYMAI,ISYMJ)
C
               NTOTAI = MAX(NT1AM(ISYMAI),1)
               NTOTEM = MAX(NT1AM(ISYMEM),1)
               NTOTJ  = MAX(NRHF(ISYMJ),1)
C
               CALL DGEMM('T','T',NT1AM(ISYMAI),NRHF(ISYMJ),
     &                    NT1AM(ISYMEM),ONE,CTR2(KOFF1),NTOTEM,
     &                    WORK(KOFF2),NTOTJ,ZERO,WORK(KOFF3),NTOTAI)
C
C-----------------------------------------------------------------
C              Final scaling with XLAMDH-element and storage in RHO2.
C-----------------------------------------------------------------
C
               DO J = 1, NRHF(ISYMJ)

                  KOFF4 = KOFF3 + NT1AM(ISYMAI)*(J-1)

                  KOFF5 = IT2SQ(ISYMAI,ISYMBJ) +
     &           NT1AM(ISYMAI)* (IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J-1))
     &                  + 1

                  CALL DGER(NT1AM(ISYMAI),NVIR(ISYMB),FACT,
     &                      WORK(KOFF4),1,WORK(KLAMD),1,
     &                      RHO2(KOFF5),NTOTAI)

               END DO
C
  120       CONTINUE

  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_L1FCK3 */
      SUBROUTINE CC_L1FCK3(RHO2,C1AM,FCKMO,ISYMV,ISYMF,WORK,LWORK)
C
C     Rasmus Faber 2017,
C     based on CC_L1FCK by Asger Halkier
C
C     Purpose: To calculate the disconnected (first) term
C              of the H contribution to RHO2. Triplet version
C
C     rho2(ai,bj) = L1am(ai)*FCKMO(jb)
C
      IMPLICIT NONE
      DOUBLE PRECISION, PARAMETER ::ONE = 1.0D0, TWO = 2.0D0
      CHARACTER(LEN=*), PARAMETER :: myname = 'CC_L1FCK3'
C
      DOUBLE PRECISION, INTENT(INOUT) :: RHO2(*), WORK(LWORK)
      DOUBLE PRECISION, INTENT(IN) :: C1AM(*), FCKMO(*)
      INTEGER, INTENT(IN) :: ISYMV, ISYMF, LWORK
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "cclr.h"
      INTEGER :: ISYMC, ISYMK, ISYMAI, ISYMBJ, ISYRES
      INTEGER :: KOFF1, KOFF2, KOFF3, NTOTAI
C
      CALL QENTER(myname)
C
      ISYRES = MULD2H(ISYMV,ISYMF)
C
      IF (LWORK .LT. NT1AM(ISYMF)) THEN
         CALL QUIT('Insufficient work-space-area in '//myname)
      ENDIF
C
C-----------------------------------
C     Extract the fock matrix F(jb).
C-----------------------------------
C
      DO 50 ISYMC = 1,NSYM
C
         ISYMK = MULD2H(ISYMC,ISYMF)
C
         DO 60 K = 1,NRHF(ISYMK)
C
            DO 70 C = 1,NVIR(ISYMC)
C
               KOFF1 = IFCVIR(ISYMK,ISYMC) + NORB(ISYMK)*(C - 1) + K
               KOFF2 = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1) + C
C
               WORK(KOFF2) = FCKMO(KOFF1)
C
   70       CONTINUE
   60    CONTINUE
   50 CONTINUE
C
C--------------------------------------------------------------
C     Calculate the contraction of C1AM with FCKMO.
C--------------------------------------------------------------
C
      ISYMAI = ISYMV
      ISYMBJ = ISYMF
      NTOTAI = MAX(1,NT1AM(ISYMAI))

      KOFF3 = IT2SQ(ISYMAI,ISYMBJ) + 1

      CALL DGER(NT1AM(ISYMAI),NT1AM(ISYMBJ),ONE,C1AM,1,WORK,1,
     &          RHO2(KOFF3),NTOTAI)
C
C     All done
      CALL QEXIT(myname)
C
      RETURN
      END

C  /* Deck cc_L1FCK3P */
      SUBROUTINE CC_L1FCK3P(RHO2,C1AM,FCKMO,ISYMV,ISYMF)
C
C     Rasmus Faber 2017,
C     based on CC_L1FCK by Asger Halkier
C
C     Purpose: To calculate the disconnected (first) term
C              of the H contribution to RHO2.
C              Triplet version, with packed RHO2
C
C     rho2+(ai,bj) = L1am(ai)*FCKMO(jb) + FCKMO(ai)*L1am(bj)
C     rho2-(ai,bj) = L1am(ai)*FCKMO(jb) - FCKMO(ai)*L1am(bj)
C
      IMPLICIT NONE
      CHARACTER(LEN=*), PARAMETER :: myname = 'CC_L1FCK3P'
C
      DOUBLE PRECISION, INTENT(INOUT) :: RHO2(*)
      DOUBLE PRECISION, INTENT(IN) :: C1AM(*), FCKMO(*)
      INTEGER, INTENT(IN) :: ISYMV, ISYMF
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "cclr.h"
      INTEGER :: ISYMAI, ISYMBJ, ISYRES
      INTEGER :: NAI, NBJ
      INTEGER :: KP, KM, KRHOP, KRHOM
      INTEGER :: KOFFP, KOFFM, KOFF3, NTOTAI
      DOUBLE PRECISION :: TEMP, TEMP1, TEMP2
C
      CALL QENTER(myname)
C
      ISYRES = MULD2H(ISYMV,ISYMF)
      KRHOP = 1
      KRHOM = NT2AM(ISYRES) + 1
C
C--------------------------------------------------------------
C     Calculate the contraction of C1AM with FCKMO.
C--------------------------------------------------------------
C
      ISYMAI = ISYMV
      ISYMBJ = ISYMF

      IF (ISYMV.EQ.ISYMF) THEN
         ISYMAI = ISYMV
         ISYMBJ = ISYMAI
         DO NBJ = 1, NT1AM(ISYMBJ)
            KOFFP = KRHOP + IT2AM(ISYMAI,ISYMBJ) + NBJ*(NBJ-1)/2
            KOFFM = KRHOM + IT2AM(ISYMAI,ISYMBJ) + NBJ*(NBJ-1)/2
            DO NAI = 1, NBJ-1
               KP = KOFFP + (NAI-1)
               KM = KOFFM + (NAI-1)
               TEMP1 = C1AM(NAI)*FCKMO(NBJ)
               TEMP2 = FCKMO(NAI)*C1AM(NBJ)
               RHO2(KP) = RHO2(KP) + TEMP1 + TEMP2
               RHO2(KM) = RHO2(KM) + TEMP1 - TEMP2
            END DO
            RHO2(KOFFP+NBJ-1) = 0.0D0
            RHO2(KOFFM+NBJ-1) = 0.0D0
        END DO
      ELSE IF (ISYMV.LT.ISYMF) THEN ! only first term is non-zero
         ISYMAI = ISYMV
         ISYMBJ = ISYMF
         DO NBJ = 1, NT1AM(ISYMBJ)
            KOFFP = KRHOP + IT2AM(ISYMAI,ISYMBJ) + NT1AM(ISYMAI)*(NBJ-1)
            KOFFM = KRHOM + IT2AM(ISYMAI,ISYMBJ) + NT1AM(ISYMAI)*(NBJ-1)
            DO NAI = 1, NT1AM(ISYMAI)
               KP = KOFFP + (NAI-1)
               KM = KOFFM + (NAI-1)
               TEMP = C1AM(NAI)*FCKMO(NBJ)
               RHO2(KP) = RHO2(KP) + TEMP
               RHO2(KM) = RHO2(KM) + TEMP
            END DO
         END DO
      ELSE ! Only the second term is stored
         ISYMAI = ISYMF
         ISYMBJ = ISYMV
         DO NBJ = 1, NT1AM(ISYMBJ)
            KOFFP = KRHOP + IT2AM(ISYMAI,ISYMBJ) + NT1AM(ISYMAI)*(NBJ-1)
            KOFFM = KRHOM + IT2AM(ISYMAI,ISYMBJ) + NT1AM(ISYMAI)*(NBJ-1)
            DO NAI = 1, NT1AM(ISYMAI)
               KP = KOFFP + (NAI-1)
               KM = KOFFM + (NAI-1)
               TEMP = C1AM(NBJ)*FCKMO(NAI)
               RHO2(KP) = RHO2(KP) + TEMP
               RHO2(KM) = RHO2(KM) - TEMP
            END DO
         END DO
      END IF
C
C     All done
      CALL QEXIT(myname)
C
      RETURN
      END
