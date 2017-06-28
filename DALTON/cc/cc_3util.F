C  /* Deck cc_t2sq3a */
      SUBROUTINE CC_T2SQ3A(T2AM,T2SQ,ISYM,FACT)
!
!--------------------------------------------------------
!     Rasmus Faber 2017:
!
!     Adds "FACT" times the content of the antisymmetric,
!     packed array T2AM to the square array T2SQ
!
!     Based on CC_T2SQ by Kasper Hald
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
!
            IF (NAMP .GT. 0) THEN
               KOFF1 = IT2SQ(ISYMAI,ISYMBJ) + 1
               CALL DAXPY(NAMP,FACT,T2AM(KOFF),1,T2SQ(KOFF1),1)
               NAI = MAX(NT1AM(ISYMAI),1)
               NBJ = MAX(NT1AM(ISYMBJ),1)
               KOFF2 = IT2SQ(ISYMBJ,ISYMAI) + 1
               CALL TRM3A(T2AM(KOFF),NT1AM(ISYMAI),NT1AM(ISYMBJ),
     *                    T2SQ(KOFF2),ONEM*FACT)
               KOFF = KOFF + NAMP
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
      SUBROUTINE TRM3A(A,M,N,B,FACT)
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
      DOUBLE PRECISION, INTENT(IN) :: A(M,N), FACT
      DOUBLE PRECISION, INTENT(INOUT) :: B(N,M)
!
      DOUBLE PRECISION :: XMONE
      INTEGER :: I
!
      DO 100 I = 1, M
!
         CALL DAXPY(N,FACT,A(I,1),M,B(1,I),1)
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
               SQMAT(I,J) = -FACT * PKMAT(IJ) + SQMAT(I,J)
               SQMAT(J,I) =  FACT*PKMAT(IJ)
  110    CONTINUE
         IJ = IJ + 1
         SQMAT(I,I) = 0.0D0
  100 CONTINUE
!
      RETURN
      END SUBROUTINE

      END
C  /* Deck cc_21b12ca */
      SUBROUTINE CC_21B12C3(RHO1,RHO2,CTR1,ISYMV,XLAMDH,ISYMH,
     *                      XMAT,ISYMX,ISINT,WORK,LWORK,
     *                      LUO3,O3FIL,IOPT)
C
C     Rasmus Faber - 2017
C
C
C
C     Based on CC_21B12C by Asger Halkier & Henrik Koch 13/9 - 1995.
C     Ove Christiansen 13-6-1996:
C
C     Generalised to calculate contributions to F matrix * vector.
C     ISYMV is symmetry of CTR1.
C     ISINT is symmetry of intermediate on disk.
C     ISYMX is symmetry of intermediate on disk and XMAT.
C     ISYMH is symmetry of lambda matrix.
C
C     rho2(ai,bj) = P(ai,bj)(-sum(k)CTR1(a,k)*L(jbik)
C                   where k may be a barred index.
C
C     IOPT controls if both 21B and 12 C are calculated:
C     if IOPT = 2 both are calculated if IOPT = 1 only
C     21B is calculated.
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION RHO1(*),RHO2(*),CTR1(*),XLAMDH(*),WORK(LWORK),XMAT(*)
      CHARACTER O3FIL*(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
#include "cclr.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC_21B12C')
C
      ISYMLI = MULD2H(ISINT,ISYMH)
      ISYRES = MULD2H(ISYMV,ISYMLI)
      IF (IOPT .EQ. 1) ISYRES = MULD2H(ISYMX,ISYMLI)
      IF ((IOPT .EQ. 2) .AND. (ISYMV .NE. ISYMX)) THEN
         WRITE(LUPRI,*) ' Symmetry mismatch in CC_21B12C '
         CALL QUIT(' Symmetry mismatch in CC_21B12C ')
      ENDIF
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
         KSCR2 = KSCR1 + NMAIJK(ISYIKJ)*NVIR(ISYMB)
         KSCR3 = KSCR2 + NMAIJK(ISYIKJ)
         KEND1 = KSCR3 + NMAIJK(ISYIKJ)*NBAS(ISYMD)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient work space area in CC_21B12C')
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
C-------------------------------------------------------------
C        Calculate 2 coulomb - exchange for (ik|jb) integrals.
C-------------------------------------------------------------
C
         DO 110 B = 1,NVIR(ISYMB)
C
            DO 120 ISYMJ = 1,NSYM
C
               ISYMIK = MULD2H(ISYIKJ,ISYMJ)
C
               DO 130 J = 1,NRHF(ISYMJ)
C
                  DO 140 ISYMK = 1,NSYM
C
                     ISYMI  = MULD2H(ISYMIK,ISYMK)
                     ISYMJK = MULD2H(ISYMJ,ISYMK)
C
                     DO 150 K = 1,NRHF(ISYMK)
C
                        DO 160 I = 1,NRHF(ISYMI)
C
                           NIKJ  = IMAIJK(ISYMIK,ISYMJ)
     *                           + NMATIJ(ISYMIK)*(J - 1)
     *                           + IMATIJ(ISYMI,ISYMK)
     *                           + NRHF(ISYMI)*(K - 1) + I
                           NJKI  = IMAIJK(ISYMJK,ISYMI)
     *                           + NMATIJ(ISYMJK)*(I - 1)
     *                           + IMATIJ(ISYMJ,ISYMK)
     *                           + NRHF(ISYMJ)*(K - 1) + J
                           NBIKJ = NMAIJK(ISYIKJ)*(B - 1) + NIKJ
                           NBJKI = NMAIJK(ISYIKJ)*(B - 1) + NJKI
C
                           WORK(KSCR2 + NIKJ - 1) =
C     *                              TWO*WORK(KSCR1 + NBIKJ - 1) -
     *                              WORK(KSCR1 + NBJKI - 1)
C
  160                   CONTINUE
  150                CONTINUE
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
C
            DO 170 ISYMJ = 1,NSYM
C
               ISYMBJ = MULD2H(ISYMJ,ISYMB)
               ISYMAI = MULD2H(ISYMBJ,ISYRES)
               ISYMIK = MULD2H(ISYMJ,ISYIKJ)
C
               LAVL   = LWORK - KSCR3 - NT1AM(ISYMAI)
C
               IF (LAVL .LT. 0) THEN
                  CALL QUIT('Insufficient space in CC_21B12C')
               ENDIF
C
               DO 180 J = 1,NRHF(ISYMJ)
C
C------------------------------------------------
C                 Calculate the 21B contribution.
C------------------------------------------------
C
                  NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
                  IF (.NOT.(CC2)) THEN
C
                     IF (ISYMBJ .EQ. ISYRES) THEN
C
                        KOFF1 = KSCR2 + IMAIJK(ISYMIK,ISYMJ)
     *                                + NMATIJ(ISYMIK)*(J - 1)
C
                        RHO1(NBJ) = RHO1(NBJ) - DDOT(NMATIJ(ISYMIK),
     *                              XMAT,1,WORK(KOFF1),1)
C
                     ENDIF
C
                  ENDIF
C
C-----------------------------------------------------------------
C                 IF only 21B contributions then goto end of loop.
C-----------------------------------------------------------------
C
                  IF ( IOPT .EQ. 1 ) GOTO 180
C
C-----------------------------------------------------------
C                 Contract integrals with trial vector CTR1.
C-----------------------------------------------------------
C
                  DO 190 ISYMK = 1,NSYM
C
                     ISYMI  = MULD2H(ISYMK,ISYMIK)
                     ISYMA  = MULD2H(ISYMK,ISYMV)
                     ISYMAJ = MULD2H(ISYMA,ISYMJ)
                     ISYMBI = MULD2H(ISYMB,ISYMI)
C
                     KOFF1 = IT1AM(ISYMA,ISYMK) + 1
                     KOFF2 = KSCR2 + IMAIJK(ISYMIK,ISYMJ)
     *                             + NMATIJ(ISYMIK)*(J - 1)
     *                             + IMATIJ(ISYMI,ISYMK)
                     KOFF3 = KSCR3 + IT1AM(ISYMA,ISYMI)
C
                     NTOTA = MAX(NVIR(ISYMA),1)
                     NTOTI = MAX(NRHF(ISYMI),1)
C
                     CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),
     *                          NRHF(ISYMK),-ONE,CTR1(KOFF1),NTOTA,
     *                          WORK(KOFF2),NTOTI,ZERO,WORK(KOFF3),
     *                          NTOTA)
C
  190             CONTINUE
C
C------------------------------------------------
C                 Calculate the 12C contribution.
C------------------------------------------------
C
                  IF (ISYMAI .EQ. ISYMBJ)
     *                      WORK(KSCR3+NBJ-1) = TWO*WORK(KSCR3+NBJ-1)
C
                  DO 200 NAI = 1,NT1AM(ISYMAI)
C
                     IF (ISYMAI .EQ. ISYMBJ) THEN
                        NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
                     ELSE IF (ISYMAI .LT. ISYMBJ) THEN
                        NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                        + NT1AM(ISYMAI)*(NBJ - 1) + NAI
                     ELSE IF (ISYMAI .GT. ISYMBJ) THEN
                        NAIBJ = IT2AM(ISYMBJ,ISYMAI)
     *                        + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
                     ENDIF
C
                     RHO2(NAIBJ) = RHO2(NAIBJ) + WORK(KSCR3+NAI-1)
C
  200             CONTINUE
  180          CONTINUE
  170       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CC_21B12C')
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

#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "cclr.h"
C
      ISYMXD = MULD2H(ISYMLH,ISYMX)
C
C------------------------------------
C     Transform Y-matrix to AO-basis.
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
         NTOTM  = MAX(NVIR(ISYMM),1)
C
         CALL DGEMM('N','N',NBAS(ISYMAL),NVIR(ISYMN),NVIR(ISYMM),
     *              ONE,XLAMDP(KOFF1),NTOTAL,XMAT(KOFF2),NTOTM,
     *              ZERO,WORK,NTOTAL)
C
         KOFF1 = IGLMRH(ISYMAL,ISYMN) + 1
         KOFF2 = IAODIS(ISYMAL,ISYMBE) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
C
         CALL DGEMM('N','T',NBAS(ISYMAL),NBAS(ISYMBE),NVIR(ISYMN),
     *              ONEM,WORK,NTOTAL,XLAMDH(KOFF1),NTOTBE,
     *              ONE,XOUT(KOFF2),NTOTAL)
C
  100 CONTINUE
C
      RETURN
      END
