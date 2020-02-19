C  /* Deck so_res_cb */
      SUBROUTINE SO_RES_CB(RES2E,LRES2E,RES2D,LRES2D,
     &                     DSRHF,LDSRHF,BTR1E,LBTR1E,BTR1D,LBTR1D,
     &                     BTJ1E,LBTJ1E,BTJ1D,LBTJ1D,CMO,LCMO,
     &                     IDEL,ISDEL,ISYDIS,ISYMTR,DO_DEX,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, February 1996
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Calculate C times b contribution to 2p2h resultvectors
C              as described in eq. (62) and (63).
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
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_RES_CB')
C
      SQ2   = DSQRT(TWO)
C
      ISYMB = ISDEL
      KCDB  = 1
      KEND1A = KCDB + NVIR(ISYMB)
      KOFF1 = ILMVIR(ISYMB) + IDEL - IBAS(ISDEL)
C
C--------------------------------------------------
C     Copy delta MO-coefficients to the vector CDB.
C--------------------------------------------------
C
      CALL DCOPY(NVIR(ISYMB),CMO(KOFF1),NBAS(ISDEL),WORK(KCDB),1)
C     The terms in this routine carries a factor of SQRT(2),
C     this can be taken care of by scaling the MO-coefficients
      CALL DSCAL(NVIR(ISYMB),SQ2,WORK(KCDB),1)
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
         KSCR1  = KEND1A
         KSCR2  = KSCR1 + LSCR1
         KSCR3  = KSCR2 + LSCR2
         KEND1  = KSCR3 + LSCR3
         LWORK1 = LWORK - KEND1
C
         CALL SO_MEMMAX ('SO_RES_CB.1',LWORK1)
         IF (LWORK1 .LT. 0) CALL STOPIT('SO_RES_CB.1',' ',KEND1,LWORK)
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
               KSCR4  = KEND1
               KSCR5  = KSCR4 + LSCR4
               KEND2  = KSCR5 + LSCR5
               LWORK2 = LWORK - KEND2
C
               CALL SO_MEMMAX ('SO_RES_CB.2',LWORK2)
               IF (LWORK2 .LT. 0)
     &             CALL STOPIT('SO_RES_CB.2',' ',KEND2,LWORK)
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
     &                    NBAS(ISALFA),ONE,
     &                    WORK(KOFF2),NTOTAL,
     &                    BTR1D(KOFF3),NTOTAL,
     &                    ZERO,WORK(KSCR5),NTOTBE)
C
               CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     &                    NBAS(ISBETA),ONE,CMO(KOFF4),NTOTBE,
     &                    WORK(KSCR5),NTOTBE,ZERO,WORK(KOFF6),NTOTA)
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
               LSCR4  = NBAS(ISALFA)*NRHF(ISYMI)
C
               KSCR4  = KEND1
               KEND3  = KSCR4 + LSCR4
               LWORK3 = LWORK - KEND3
C
               CALL SO_MEMMAX ('SO_RES_CB.3',LWORK3)
               IF (LWORK3 .LT. 0)
     &             CALL STOPIT('SO_RES_CB.3',' ',KEND3,LWORK)
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
     &                    NBAS(ISBETA),-ONE,BTJ1D(KOFF4),NTOTBE,
     &                    WORK(KSCR4),NTOTBE,ONE,WORK(KOFF6),NTOTA)
               END IF
C
  300       CONTINUE
C                       ~~
C           Calculate ( ai, j, del)*C(del,b) => X(ai,bj)
C
            NBJ1 = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1)
            KOFFRES = IT2AM(ISYMAI,ISYMBJ) + 1

            CALL OUTPACK(RES2E(KOFFRES),WORK(KCDB),WORK(KSCR2))
            IF(DO_DEX) 
     &         CALL OUTPACK(RES2D(KOFFRES),WORK(KCDB),WORK(KSCR3))


  200    CONTINUE
C
  100 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_RES_CB')
C
      RETURN
      CONTAINS
      
         SUBROUTINE OUTPACK(RES,CDB,SCR)!,ISYMAI,ISYMBJ,NBJ1)
            USE SO_INFO, ONLY: SOP_DP
            REAL(SOP_DP), INTENT(INOUT) :: RES(*)
            REAL(SOP_DP), INTENT(IN) :: CDB(*), SCR(*)
C            INTEGER, INTENT(IN) :: ISYMAI, ISYMBJ, NBJ1
            IF (ISYMAI.EQ.ISYMBJ) THEN
C              Triangular storage => complicated case
C              Calculate contributions to from ai=<bj
C              -> loop b first, then ai

               DO B = 1, NVIR(ISYMB)

                  NBJ   = NBJ1 + B
                  KOFFBJ = NBJ*(NBJ-1)/2
                  FAC   = CDB(B)

                  DO NAI = 1, NBJ
                     NAIBJ = KOFFBJ + NAI
                     RES(NAIBJ) = RES(NAIBJ)+FAC*SCR(NAI)
                  END DO
               END DO
C
C              Second term, ai>bj
C
               DO ISYMI = ISYMJ, NSYM
                  ISYMA = MULD2H(ISYMAI,ISYMI)
                  IF (ISYMI .EQ. ISYMJ) THEN
                  ! Same symmetries: two cases I>J and I=J
                  ! DO I = J Seperately
                     I = J
                     DO A = 2, NVIR(ISYMA)
                        NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1)+A

                        KOFFAI = NAI*(NAI-1)/2 + NBJ1
                     ! In this case B<A, The B>=A has allready been
                     ! handled
                        DO B = 1, A-1
                           NAIBJ = KOFFAI + B
                           RES(NAIBJ) = RES(NAIBJ) + CDB(B)*SCR(NAI)
                        END DO

                     END DO
                     ! Then I>J
                     ISTART = J + 1
                  ELSE
                     ! ISYMI > ISYMJ, so I>J already ensured
                     ISTART = 1
                  END IF

                  DO I = ISTART, NRHF(ISYMI)
                  ! Since I > J, we have no restrictions on A,B
                     DO A = 1, NVIR(ISYMA)
                        NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1)+A

                        KOFFAI = NAI*(NAI-1)/2 + NBJ1
                        ! In this case B<A, The B>=A has allready been
                        ! handled
                        DO B = 1, NVIR(ISYMB)
                           NAIBJ = KOFFAI + B
                           RES(NAIBJ) = RES(NAIBJ) + CDB(B)*SCR(NAI)
                        END DO
                     END DO
                  END DO

               END DO ! ISYMI

            ELSE IF (ISYMAI .LT. ISYMBJ) THEN
C              ai < bj ensured, only one term
C
C              Loop B first, then AI
               DO B = 1, NVIR(ISYMB)
                  NBJ   = NBJ1 + B
                  KOFFBJ = NT1AM(ISYMAI)*(NBJ-1)
                  FAC   =  CDB(B)
                  DO NAI = 1, NT1AM(ISYMAI)
                     NAIBJ = KOFFBJ + NAI
                     RES(NAIBJ) = RES(NAIBJ)+FAC*SCR(NAI)
                  END DO
               END DO

            ELSE ! ISYMAI > ISYMBJ
C              ai > bj ensured, one term!
               DO NAI = 1, NT1AM(ISYMAI)
                  KOFFAI = NT1AM(ISYMBJ)*(NAI-1) + NBJ1
                  DO B = 1, NVIR(ISYMB)
                     NAIBJ = KOFFAI + B
                     RES(NAIBJ) = RES(NAIBJ)+CDB(B)*SCR(NAI)
                  END DO
               END DO
            END IF
            
            RETURN
         END SUBROUTINE

      END
