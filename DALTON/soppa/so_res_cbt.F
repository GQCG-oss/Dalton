C
C  /* Deck so_res_cbt */
      SUBROUTINE SO_RES_CBT(RES2E,LRES2E,RES2D,LRES2D,
     &                      DSRHF,LDSRHF,BTR1E,LBTR1E,
     &                      BTR1D,LBTR1D,BTJ1E,LBTJ1E,
     &                      BTJ1D,LBTJ1D,CMO,LCMO,
     &                      IDEL,ISDEL,ISYDIS,ISYMTR,DO_DEX,
     &                      WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, February 1996
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C     Stephan P. A. Sauer: April 2006: Triplet version
C
C     Pi 29.03.16: Copied from dalton20aosoppa
C
C     PURPOSE: Calculate C times b contribution to 2p2h resultvectors
C              as described in eq. (62) and (63).
C
      use so_info, only: sop_dp
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION RES2E(LRES2E), RES2D(LRES2D)
      DIMENSION DSRHF(LDSRHF), BTR1E(LBTR1E), BTR1D(LBTR1D)
      DIMENSION BTJ1E(LBTJ1E), BTJ1D(LBTJ1D), CMO(LCMO)
      DIMENSION WORK(LWORK)
      LOGICAL, INTENT(IN) :: DO_DEX
      REAL(SOP_DP), PARAMETER :: SQ2 = DSQRT(TWO), ONESQ2 = ONE/SQ2
C
      character(len=*), parameter :: myname = 'SO_RES_CBT'
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "soppinf.h"
C
C------------------------------
C     Statement function INDEX.
C------------------------------
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
      INDLT(I,J) = (MAX(I,J) - 1)*(MAX(I,J) - 2)/2 + MIN(I,J)
C
C------------------
C     Add to trace.
C------------------
C
C
      CALL QENTER('SO_RES_CBT')
C
      ISYMB = ISDEL
C
      LCDB  = NVIR(ISYMB)
CPi If no virtual orbitals in this symmetry, skip this transformation
C      IF (LCDB .EQ. 0) THEN
C        CALL QEXIT(myname)
C        RETURN
C      END IF
C
      KCDB    = 1
      KEND0   = KCDB  + LCDB
      LWORK0  = LWORK - KEND0
C
      CALL SO_MEMMAX (myname//'.0',LWORK0)
      IF (LWORK0 .LT. 0) CALL STOPIT(myname//'.0',' ',KEND0,LWORK)
      KOFF1 = ILMVIR(ISDEL) + IDEL - IBAS(ISDEL)
C
C--------------------------------------------------
C     Copy delta MO-coefficients to the vector CDB.
C--------------------------------------------------
C
      CALL DCOPY(NVIR(ISYMB),CMO(KOFF1),NBAS(ISDEL),WORK(KCDB),1)
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMJ,ISYMB)
         ISALBE = MULD2H(ISYMJ,ISYDIS)
         ISYMAI = MULD2H(ISALBE,ISYMTR)
C
         LSCR1  = N2BST(ISALBE)
         LSCR2  = NT1AM(ISYMAI)
         LSCR3  = NT1AM(ISYMAI)
C
         KSCR1   = KEND0
         KSCR2   = KSCR1 + LSCR1
         KSCR3   = KSCR2 + LSCR2
         KEND1   = KSCR3 + LSCR3
         LWORK1  = LWORK - KEND1
C
         CALL SO_MEMMAX ('SO_RES_CBT.1',LWORK1)
         IF (LWORK1 .LT. 0) CALL STOPIT('SO_RES_CBT.1',' ',KEND1,LWORK)
C
         DO 200 J = 1,NRHF(ISYMJ)
C
            KOFF1 = IDSRHF(ISALBE,ISYMJ) + NNBST(ISALBE)*(J - 1) + 1
C
C-----------------------------------------------------------------------
C           Get a squared set of ( alfa beta | j delta ) for given j and
C           delta.
C-----------------------------------------------------------------------
C
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KSCR1))
C
            DO 300 ISYMI = 1,NSYM
C
               ISYMA  = MULD2H(ISYMI,ISYMAI)
C
C----------------------------------------------------------------------
C                                       ~ ~
C              Generate first part of ( a i | j delta ) for given j and
C              delta and given symmetry of ai in KSCR2 and KSCR3 for
C              excitations and de-excitations, respectively.
C----------------------------------------------------------------------
C
               ISBETA = ISYMA
               ISALFA = MULD2H(ISBETA,ISALBE)
C
               LSCR4  = NBAS(ISBETA)*NRHF(ISYMI)
               LSCR5  = NBAS(ISBETA)*NRHF(ISYMI)
C
               KSCR4   = KEND1
               KSCR5   = KSCR4 + LSCR4
               KEND2   = KSCR5 + LSCR5
C
               LWORK2  = LWORK - KEND2
C
               CALL SO_MEMMAX ('SO_RES_CBT.2',LWORK2)
               IF (LWORK2 .LT. 0)
     &             CALL STOPIT('SO_RES_CBT.2',' ',KEND2,LWORK)
C
               NTOTAL = MAX(NBAS(ISALFA),1)
               NTOTBE = MAX(NBAS(ISBETA),1)
               NTOTA  = MAX(NVIR(ISYMA),1)
C
               KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
               KOFF3  = IT1AO(ISALFA,ISYMI) + 1
               KOFF4  = ILMVIR(ISYMA) + 1
               KOFF5  = KSCR2 + IT1AM(ISYMA,ISYMI)
               KOFF6  = KSCR3 + IT1AM(ISYMA,ISYMI)
C
C------------------------------
C              For excitations.
C------------------------------
C                                                    ~
C              (alpha, beta | * x (alpha,i) => (beta,i|
C
               CALL DGEMM('T','N',NBAS(ISBETA),NRHF(ISYMI),
     &                    NBAS(ISALFA),ONE,
     &                    WORK(KOFF2),NTOTAL,
     &                    BTR1E(KOFF3),NTOTAL,
     &                    ZERO,WORK(KSCR4),NTOTBE)
C                                  ~         ~
C              c (beta,a) * (beta, i | => (a,i|
C
               CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     &                    NBAS(ISBETA),ONE,CMO(KOFF4),NTOTBE,
     &                    WORK(KSCR4),NTOTBE,ZERO,WORK(KOFF5),NTOTA)
C
C---------------------------------
C              For de-excitations.
C---------------------------------
C
               IF (DO_DEX) THEN
                  CALL DGEMM('T','N',NBAS(ISBETA),NRHF(ISYMI),
     &                       NBAS(ISALFA),ONE,
     &                       WORK(KOFF2),NTOTAL,
     &                       BTR1D(KOFF3),NTOTAL,
     &                       ZERO,WORK(KSCR5),NTOTBE)
C
                  CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     &                       NBAS(ISBETA),ONE,CMO(KOFF4),NTOTBE,
     &                       WORK(KSCR5),NTOTBE,ZERO,WORK(KOFF6),NTOTA)
               END IF
C
C-----------------------------------------------------------------------
C                                        ~ ~
C              Generate second part of ( a i | j delta ) for given j and
C              delta and given symmetry of ai and add to KSCR2 and KSCR3
C              for excitations and de-excitations, respectively.
C-----------------------------------------------------------------------
C
               ISALFA = ISYMI
               ISBETA = MULD2H(ISALFA,ISALBE)
C
               LSCR4  = NBAS(ISBETA)*NRHF(ISYMI)
C
               KSCR4   = KEND1
               KEND3   = KSCR4 + LSCR4
               LWORK3  = LWORK - KEND3
C
               CALL SO_MEMMAX ('SO_RES_CBT.3',LWORK3)
               IF (LWORK3 .LT. 0)
     &             CALL STOPIT('SO_RES_CBT.3',' ',KEND3,LWORK)
C
               NTOTAL = MAX(NBAS(ISALFA),1)
               NTOTBE = MAX(NBAS(ISBETA),1)
               NTOTA  = MAX(NVIR(ISYMA),1)
C
               KOFF2  = KSCR1 + IAODIS(ISALFA,ISBETA)
               KOFF3  = ILMRHF(ISYMI) + 1
               KOFF4  = IMATAV(ISBETA,ISYMA) + 1
               KOFF5  = KSCR2 + IT1AM(ISYMA,ISYMI)
               KOFF6  = KSCR3 + IT1AM(ISYMA,ISYMI)
C
C              ( alpha, beta | -- ) * c(alpha, i) => ( beta, i | --)
               CALL DGEMM('T','N',NBAS(ISBETA),NRHF(ISYMI),
     &                    NBAS(ISALFA),ONE,WORK(KOFF2),NTOTAL,
     &                    CMO(KOFF3),NTOTAL,ZERO,WORK(KSCR4),NTOTBE)
C
C------------------------------
C              For excitations.
C------------------------------
C                                                  ~
C              ( beta, i | -- ) * x(beta, a ) => ( a, i | --)
               CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     &                    NBAS(ISBETA),-ONE,BTJ1E(KOFF4),NTOTBE,
     &                    WORK(KSCR4),NTOTBE,ONE,WORK(KOFF5),NTOTA)
C
C---------------------------------
C              For de-excitations.
C---------------------------------
C
               IF (DO_DEX) THEN
                  CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     &                       NBAS(ISBETA),-ONE,BTJ1D(KOFF4),NTOTBE,
     &                       WORK(KSCR4),NTOTBE,ONE,WORK(KOFF6),NTOTA)
               END IF


               ISYMAB = MULD2H(ISYMA,ISYMB)
               ISYMIJ = MULD2H(ISYMI,ISYMJ)

              
               KPOSIJ1 = IT2AMT1(ISYMIJ,ISYMAB)
     &                 + NTVV(ISYMAB)*ITOO(ISYMI,ISYMJ)
     &                 + ITVV(ISYMA,ISYMB) + 1
               KPOSIJ2 = IT2AMT2(ISYMIJ,ISYMAB)
     &                 + NSVV(ISYMAB)*ITOO(ISYMI,ISYMJ)
     &                 + ISVV(ISYMA,ISYMB) + 1
               KPOSIJ3 = IT2AMT3(ISYMIJ,ISYMAB)
     &                 + NTVV(ISYMAB)*ISOO(ISYMI,ISYMJ)
     &                 + ITVV(ISYMA,ISYMB) + 1

               CALL OUTPACK(RES2E(KPOSIJ1),RES2E(KPOSIJ2),
     &                      RES2E(KPOSIJ3),WORK(KCDB),WORK(KOFF5),J,
     &                      ISYMA,ISYMB,ISYMI,ISYMJ)

               IF (DO_DEX) THEN
                  CALL OUTPACK(RES2D(KPOSIJ1),RES2D(KPOSIJ2),
     &                         RES2D(KPOSIJ3),WORK(KCDB),WORK(KOFF6),J,
     &                         ISYMA,ISYMB,ISYMI,ISYMJ)
               END IF

C
  300       CONTINUE
C
  200    CONTINUE
C
  100 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_RES_CBT')
C
      RETURN

      CONTAINS

         SUBROUTINE OUTPACK(RES1,RES2,RES3,CDB,SCR,J,
     &                      ISYMA,ISYMB,ISYMI,ISYMJ)
            IMPLICIT NONE
            REAL(SOP_DP), INTENT(INOUT) :: RES1(*), RES2(*), RES3(*)
            REAL(SOP_DP), INTENT(IN) :: CDB(*), SCR(*)
            INTEGER, INTENT(IN) :: J, ISYMA, ISYMB, ISYMI, ISYMJ
             
            INTEGER :: NIJ1, NIJ2, NIJ3, KSCR
            REAL(SOP_DP) :: FACTIJ, FACTAB
            INTEGER :: KV1, KV10, KV2, KV20, LOFFV2, LOFFV1
            INTEGER :: NVIR1, NVIR2 
            INTEGER :: I, ISYMAB
            INTEGER :: NRHFI, NRHFJ

            ISYMAB = MULD2H(ISYMA,ISYMB)

            IF (ISYMI.EQ.ISYMJ) THEN
               IF (ISYMA.EQ.ISYMB) THEN
C  Case I<J
                  FACTIJ = -1.0D0
                  DO I = 1, J-1
                     NIJ1 = NTVV(ISYMAB)*(((J-2)*(J-1))/2+I-1)+1
                     NIJ2 = NSVV(ISYMAB)*(((J-2)*(J-1))/2+I-1)+1
                     NIJ3 = NTVV(ISYMAB)*((J*(J-1))/2+I-1)+1
                     KSCR = NVIR(ISYMA)*(I-1) + 1
                     CALL EXPAND_SYM(RES1(NIJ1),RES2(NIJ2),
     &                               RES3(NIJ3),CDB,
     &                               SCR(KSCR),FACTIJ,NVIR(ISYMB))
                  END DO
C  Handle I==J
                  NIJ3 = NTVV(ISYMAB)*((J*(J+1))/2-1)+1
                  KSCR = NVIR(ISYMA)*(J-1) + 1
                  CALL EXPAND_SYM1(RES3(NIJ3),CDB,
     &                             SCR(KSCR),NVIR(ISYMB))
C  Case I>I
                  FACTIJ = 1.D0
                  DO I = J+1, NRHF(ISYMI)
                     NIJ1 = NTVV(ISYMAB)*(((I-2)*(I-1))/2+J-1)+1
                     NIJ2 = NSVV(ISYMAB)*(((I-2)*(I-1))/2+J-1)+1
                     NIJ3 = NTVV(ISYMAB)*((I*(I-1))/2+J-1)+1
                     KSCR = NVIR(ISYMA)*(I-1) + 1
                     CALL EXPAND_SYM(RES1(NIJ1),RES2(NIJ2),
     &                               RES3(NIJ3),CDB,
     &                               SCR(KSCR),FACTIJ,NVIR(ISYMB))
                  END DO

               ELSE
C
                  IF (ISYMA.GT.ISYMB) THEN
                     FACTAB = 1.0D0
C  Case I<J
                     FACTIJ = -1.0D0
                     DO I = 1, J-1
                        NIJ1 = NTVV(ISYMAB)*(((J-2)*(J-1))/2+I-1)+1
                        NIJ2 = NSVV(ISYMAB)*(((J-2)*(J-1))/2+I-1)+1
                        NIJ3 = NTVV(ISYMAB)*((J*(J-1))/2+I-1)+1
                        KSCR = NVIR(ISYMA)*(I-1) + 1
                        CALL EXPAND_GEN(RES1(NIJ1),RES2(NIJ2),
     &                                  RES3(NIJ3),
     &                                  SCR(KSCR),CDB,FACTIJ,FACTAB,
     &                                  NVIR(ISYMA),NVIR(ISYMB))
                     END DO
C  Handle I==J
                     NIJ3 = NTVV(ISYMAB)*((J*(J+1))/2-1)+1
                     KSCR = NVIR(ISYMA)*(J-1) + 1
                     CALL EXPAND_GEN1(RES3(NIJ3),
     &                             SCR(KSCR),CDB,FACTAB,
     &                             NVIR(ISYMA),NVIR(ISYMB))
C  Case I>I
                     FACTIJ = 1.D0
                     DO I = J+1, NRHF(ISYMI)
                        NIJ1 = NTVV(ISYMAB)*(((I-2)*(I-1))/2+J-1)+1
                        NIJ2 = NSVV(ISYMAB)*(((I-2)*(I-1))/2+J-1)+1
                        NIJ3 = NTVV(ISYMAB)*((I*(I-1))/2+J-1)+1
                        KSCR = NVIR(ISYMA)*(I-1) + 1
                        CALL EXPAND_GEN(RES1(NIJ1),RES2(NIJ2),
     &                                  RES3(NIJ3),
     &                                  SCR(KSCR),CDB,FACTIJ,FACTAB,
     &                                  NVIR(ISYMA),NVIR(ISYMB))
                     END DO
                  ELSE
                     FACTAB = -1.0D0
                     FACTIJ = -1.0D0
                     DO I = 1, J-1
                        NIJ1 = NTVV(ISYMAB)*(((J-2)*(J-1))/2+I-1)+1
                        NIJ2 = NSVV(ISYMAB)*(((J-2)*(J-1))/2+I-1)+1
                        NIJ3 = NTVV(ISYMAB)*((J*(J-1))/2+I-1)+1
                        KSCR = NVIR(ISYMA)*(I-1) + 1
                        CALL EXPAND_GEN(RES1(NIJ1),RES2(NIJ2),
     &                                  RES3(NIJ3),CDB,
     &                                  SCR(KSCR),FACTIJ,FACTAB,
     &                                  NVIR(ISYMB),NVIR(ISYMA))
                     END DO
C  Handle I==J
                     NIJ3 = NTVV(ISYMAB)*((J*(J+1))/2-1)+1
                     KSCR = NVIR(ISYMA)*(J-1) + 1
                     CALL EXPAND_GEN1(RES3(NIJ3),CDB,
     &                             SCR(KSCR),FACTAB,
     &                             NVIR(ISYMB),NVIR(ISYMA))
C  Case I>I
                     FACTIJ = 1.D0
                     DO I = J+1, NRHF(ISYMI)
                        NIJ1 = NTVV(ISYMAB)*(((I-2)*(I-1))/2+J-1)+1
                        NIJ2 = NSVV(ISYMAB)*(((I-2)*(I-1))/2+J-1)+1
                        NIJ3 = NTVV(ISYMAB)*((I*(I-1))/2+J-1)+1
                        KSCR = NVIR(ISYMA)*(I-1) + 1
                        CALL EXPAND_GEN(RES1(NIJ1),RES2(NIJ2),
     &                                  RES3(NIJ3),CDB,
     &                                  SCR(KSCR),FACTIJ,FACTAB,
     &                                  NVIR(ISYMB),NVIR(ISYMA))
                     END DO
                  END IF
               END IF
            ELSE
               IF (ISYMI.GT.ISYMJ) THEN
                  FACTIJ = 1.0D0
                  NRHFI = 1
                  NRHFJ = NRHF(ISYMJ)
               ELSE
                  FACTIJ =-1.0D0
                  NRHFI = NRHF(ISYMI)
                  NRHFJ = 1
               END IF
C
               IF (ISYMA.EQ.ISYMB) THEN
                  DO I = 1, NRHF(ISYMI)
                     NIJ1 = NTVV(ISYMAB)*(NRHFI*(J-1)+NRHFJ*(I-1))+1
                     NIJ2 = NSVV(ISYMAB)*(NRHFI*(J-1)+NRHFJ*(I-1))+1
                     NIJ3 = NTVV(ISYMAB)*(NRHFI*(J-1)+NRHFJ*(I-1))+1
                     KSCR = NVIR(ISYMA)*(I-1) + 1
                     CALL EXPAND_SYM(RES1(NIJ1),RES2(NIJ2),
     &                               RES3(NIJ3),CDB,
     &                               SCR(KSCR),FACTIJ,NVIR(ISYMB))
                  END DO
               ELSEIF (ISYMA.GT.ISYMB) THEN
                  FACTAB = 1.0D0
                  DO I = 1, NRHF(ISYMI)
                     NIJ1 = NTVV(ISYMAB)*(NRHFI*(J-1)+NRHFJ*(I-1))+1
                     NIJ2 = NSVV(ISYMAB)*(NRHFI*(J-1)+NRHFJ*(I-1))+1
                     NIJ3 = NTVV(ISYMAB)*(NRHFI*(J-1)+NRHFJ*(I-1))+1
                     KSCR = NVIR(ISYMA)*(I-1) + 1
                     CALL EXPAND_GEN(RES1(NIJ1),RES2(NIJ2),
     &                               RES3(NIJ3),SCR(KSCR),
     &                               CDB,FACTIJ,FACTAB,
     &                               NVIR(ISYMA),NVIR(ISYMB))
                  END DO
               ELSE
                  FACTAB = -1.0D0
                  DO I = 1, NRHF(ISYMI)
                     NIJ1 = NTVV(ISYMAB)*(NRHFI*(J-1)+NRHFJ*(I-1))+1
                     NIJ2 = NSVV(ISYMAB)*(NRHFI*(J-1)+NRHFJ*(I-1))+1
                     NIJ3 = NTVV(ISYMAB)*(NRHFI*(J-1)+NRHFJ*(I-1))+1
                     KSCR = NVIR(ISYMA)*(I-1) + 1
                     CALL EXPAND_GEN(RES1(NIJ1),RES2(NIJ2),
     &                               RES3(NIJ3),CDB,
     &                               SCR(KSCR),FACTIJ,FACTAB,
     &                               NVIR(ISYMB),NVIR(ISYMA))
                  END DO
               END IF

            END IF
         END SUBROUTINE
C 
         SUBROUTINE EXPAND_SYM(RES1,RES2,RES3,CMO,INTAI,FACTIJ,NB)
            IMPLICIT NONE
            REAL(SOP_DP),INTENT(INOUT) :: RES1(*), RES2(*), RES3(*)
            REAL(SOP_DP),INTENT(IN) :: CMO(*), INTAI(*), FACTIJ
            INTEGER, INTENT(IN) :: NB

            REAL(SOP_DP) :: TMP1, TMP2
            INTEGER :: POS1, POS2, A, B
            POS1 = 0
            POS2 = 0
            DO B = 1, NB
               DO A = 1, B-1
                  POS1 = POS1 + 1
                  POS2 = POS2 + 1
                  TMP1 = CMO(B)*INTAI(A)
                  TMP2 = CMO(A)*INTAI(B)
                  RES1(POS1) = RES1(POS1) - FACTIJ*(TMP2-TMP1)
                  RES2(POS2) = RES2(POS2) + FACTIJ*(TMP2+TMP1)*ONESQ2
                  RES3(POS1) = RES3(POS1) + ONESQ2*(TMP2-TMP1)
               END DO
C
               POS2 = POS2 + 1
               RES2(POS2) = RES2(POS2) + FACTIJ*CMO(B)*INTAI(B)
            END DO
            RETURN
         END SUBROUTINE
C
         SUBROUTINE EXPAND_GEN(RES1,RES2,RES3,V1,V2,FACTIJ,FACTAB,NB,NA)
            IMPLICIT NONE
            REAL(SOP_DP),INTENT(INOUT) :: RES1(*), RES2(*), RES3(*)
            REAL(SOP_DP),INTENT(IN) :: V1(*), V2(*), FACTIJ, FACTAB
            INTEGER, INTENT(IN) :: NA, NB

            REAL(SOP_DP) :: TMP
            INTEGER :: POS1, A, B
            POS1 = 0
            DO B = 1, NB
               DO A = 1, NA
                  POS1 = POS1 + 1
                  TMP = V1(B)*V2(A)
                  RES1(POS1) = RES1(POS1) - FACTAB*FACTIJ*TMP
                  RES2(POS1) = RES2(POS1) + FACTIJ*ONESQ2*TMP
                  RES3(POS1) = RES3(POS1) + FACTAB*ONESQ2*TMP
               END DO
            END DO
            RETURN
         END SUBROUTINE
C
         SUBROUTINE EXPAND_SYM1(RES3,CMO,INTAI,NB)
            IMPLICIT NONE
            REAL(SOP_DP),INTENT(INOUT) :: RES3(*)
            REAL(SOP_DP),INTENT(IN) :: CMO(*), INTAI(*)
            INTEGER, INTENT(IN) :: NB

            REAL(SOP_DP) :: TMP1, TMP2
            INTEGER :: POS1, A, B
            POS1 = 0
            DO B = 1, NB
               DO A = 1, B-1
                  POS1 = POS1 + 1
                  TMP1 = CMO(B)*INTAI(A)
                  TMP2 = CMO(A)*INTAI(B)
                  RES3(POS1) = RES3(POS1) + (TMP2-TMP1)
               END DO
            END DO
            RETURN
         END SUBROUTINE
C         
         SUBROUTINE EXPAND_GEN1(RES3,V1,V2,FACTAB,NB,NA)
            IMPLICIT NONE
            REAL(SOP_DP),INTENT(INOUT) :: RES3(*)
            REAL(SOP_DP),INTENT(IN) :: V1(*), V2(*), FACTAB
            INTEGER, INTENT(IN) :: NA, NB

            REAL(SOP_DP) :: TMP
            INTEGER :: POS1, A, B
            POS1 = 0
            DO B = 1, NB
               DO A = 1, NA
                  POS1 = POS1 + 1
                  TMP = V1(B)*V2(A)
                  RES3(POS1) = RES3(POS1) + FACTAB*TMP
               END DO
            END DO
            RETURN
         END SUBROUTINE
      END
