C  /* Deck so_res_cb */
      SUBROUTINE SO_RES_CB(RES2E,LRES2E,RES2D,LRES2D,
     &                     DSRHF,LDSRHF,BTR1E,LBTR1E,BTR1D,LBTR1D,
     &                     BTJ1E,LBTJ1E,BTJ1D,LBTJ1D,CMO,LCMO,
     &                     IDEL,ISDEL,ISYDIS,ISYMTR,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, February 1996
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Calculate C times b contribution to 2p2h resultvectors
C              as described in eq. (62) and (63).
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION RES2E(LRES2E), RES2D(LRES2D)
      DIMENSION DSRHF(LDSRHF), BTR1E(LBTR1E), BTR1D(LBTR1D)
      DIMENSION BTJ1E(LBTJ1E), BTJ1D(LBTJ1D), CMO(LCMO)
      DIMENSION WORK(LWORK)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
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
      CALL QENTER('SO_RES_CB')
C
      SQ2   = DSQRT(TWO)
C
      ISYMB = ISDEL
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
         KSCR1   = 1
         KSCR2   = KSCR1 + LSCR1
         KSCR3   = KSCR2 + LSCR2
         KEND1   = KSCR3 + LSCR3
         LWORK1  = LWORK - KEND1
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
               ISALFA = ISYMA
               ISBETA = MULD2H(ISALFA,ISALBE)
C
               LSCR4  = NBAS(ISALFA)*NRHF(ISYMI)
               LSCR5  = NBAS(ISALFA)*NRHF(ISYMI)
C
               KSCR4   = KEND1
               KSCR5   = KSCR4 + LSCR4
               KEND2   = KSCR5 + LSCR5
               LWORK2  = LWORK - KEND2
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
               KOFF3  = IT1AO(ISBETA,ISYMI) + 1
               KOFF4  = ILMVIR(ISYMA) + 1
               KOFF5  = KSCR2 + IT1AM(ISYMA,ISYMI)
               KOFF6  = KSCR3 + IT1AM(ISYMA,ISYMI)
C
C------------------------------
C              For excitations.
C------------------------------
C
               CALL DGEMM('N','N',NBAS(ISALFA),NRHF(ISYMI),
     &                    NBAS(ISBETA),ONE,WORK(KOFF2),NTOTAL,
     &                    BTR1E(KOFF3),NTOTBE,ZERO,WORK(KSCR4),NTOTAL)
C
               CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     &                    NBAS(ISALFA),ONE,CMO(KOFF4),NTOTAL,
     &                    WORK(KSCR4),NTOTAL,ZERO,WORK(KOFF5),NTOTA)
C
C---------------------------------
C              For de-excitations.
C---------------------------------
C
               CALL DGEMM('N','N',NBAS(ISALFA),NRHF(ISYMI),
     &                    NBAS(ISBETA),ONE,WORK(KOFF2),NTOTAL,
     &                    BTR1D(KOFF3),NTOTBE,ZERO,WORK(KSCR5),NTOTAL)
C
               CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     &                    NBAS(ISALFA),ONE,CMO(KOFF4),NTOTAL,
     &                    WORK(KSCR5),NTOTAL,ZERO,WORK(KOFF6),NTOTA)
C
C-----------------------------------------------------------------------
C                                        ~ ~
C              Generate second part of ( a i | j delta ) for given j and
C              delta and given symmetry of ai and add to KSCR2 and KSCR3
C              for excitations and de-excitations, respectively.
C-----------------------------------------------------------------------
C
               ISBETA = ISYMI
               ISALFA = MULD2H(ISBETA,ISALBE)
C
               LSCR4  = NBAS(ISALFA)*NRHF(ISYMI)
C
               KSCR4   = KEND1
               KEND3   = KSCR4 + LSCR4
               LWORK3  = LWORK - KEND3
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
               KOFF4  = IMATAV(ISALFA,ISYMA) + 1
               KOFF5  = KSCR2 + IT1AM(ISYMA,ISYMI)
               KOFF6  = KSCR3 + IT1AM(ISYMA,ISYMI)
C
               CALL DGEMM('N','N',NBAS(ISALFA),NRHF(ISYMI),
     &                    NBAS(ISBETA),ONE,WORK(KOFF2),NTOTAL,
     &                    CMO(KOFF3),NTOTBE,ZERO,WORK(KSCR4),NTOTAL)
C
C------------------------------
C              For excitations.
C------------------------------
C
               CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     &                    NBAS(ISALFA),-ONE,BTJ1E(KOFF4),NTOTAL,
     &                    WORK(KSCR4),NTOTAL,ONE,WORK(KOFF5),NTOTA)
C
C---------------------------------
C              For de-excitations.
C---------------------------------
C
               CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     &                    NBAS(ISALFA),-ONE,BTJ1D(KOFF4),NTOTAL,
     &                    WORK(KSCR4),NTOTAL,ONE,WORK(KOFF6),NTOTA)
C
  300       CONTINUE
C
C--------------------------------------
C           4. allocation of workspace.
C--------------------------------------
C
            LSCR4  = NT1AM(ISYMAI)
            LSCR5  = NT1AM(ISYMAI)
C
            KSCR4   = KEND1
            KSCR5   = KSCR4 + LSCR4
            KEND4   = KSCR5 + LSCR5
            LWORK4  = LWORK - KEND4
C
            CALL SO_MEMMAX ('SO_RES_CB.4',LWORK4)
            IF (LWORK4 .LT. 0) 
     &          CALL STOPIT('SO_RES_CB.4',' ',KEND4,LWORK)
C
            DO 400 B = 1,NVIR(ISYMB)
C
               KOFF6 = ILMVIR(ISYMB) + NBAS(ISDEL)*(B - 1)
     &               + ( IDEL - IBAS(ISDEL) )
C
               NBJ   = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
C--------------------------------------------------------------
C                         ~ ~
C              Generate ( a i | j b ) for given j and b for
C              excitations and de-excitations in SCR4 and SCR5.
C--------------------------------------------------------------
C
               CALL DCOPY(NT1AM(ISYMAI),WORK(KSCR2),1,WORK(KSCR4),1)
               CALL DCOPY(NT1AM(ISYMAI),WORK(KSCR3),1,WORK(KSCR5),1)
C
               FAC   = SQ2 * CMO(KOFF6)
C
               CALL DSCAL(NT1AM(ISYMAI),FAC,WORK(KSCR4),1)
               CALL DSCAL(NT1AM(ISYMAI),FAC,WORK(KSCR5),1)
C
C------------------------------------------------------------
C                          ~ ~
C              Add SQ2 * ( a i | j b ) to 2p2h excitation and
C              de-excitation result-vectors.
C------------------------------------------------------------
C
               IF ( ISYMAI .EQ. ISYMBJ ) THEN
C
                  DO 501 NAI = 1,NT1AM(ISYMAI)
C
                     NAIBJ        = IT2AM(ISYMAI,ISYMBJ)
     &                            + INDEX(NAI,NBJ)
C
                     RES2E(NAIBJ) = RES2E(NAIBJ) + WORK(KSCR4 + NAI - 1)
                     RES2D(NAIBJ) = RES2D(NAIBJ) + WORK(KSCR5 + NAI - 1)
C
  501             CONTINUE
C
               ELSE IF ( ISYMAI .LT. ISYMBJ ) THEN
C
                  DO 502 NAI = 1,NT1AM(ISYMAI)
C
                     NAIBJ        = IT2AM(ISYMAI,ISYMBJ)
     &                            + NT1AM(ISYMAI)*(NBJ - 1) + NAI
C
                     RES2E(NAIBJ) = RES2E(NAIBJ) + WORK(KSCR4 + NAI - 1)
                     RES2D(NAIBJ) = RES2D(NAIBJ) + WORK(KSCR5 + NAI - 1)
C
  502             CONTINUE
C
               ELSE IF ( ISYMAI .GT. ISYMBJ ) THEN
C
                  DO 503 NAI = 1,NT1AM(ISYMAI)
C
                     NAIBJ        = IT2AM(ISYMBJ,ISYMAI)
     &                            + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
C
                     RES2E(NAIBJ) = RES2E(NAIBJ) + WORK(KSCR4 + NAI - 1)
                     RES2D(NAIBJ) = RES2D(NAIBJ) + WORK(KSCR5 + NAI - 1)
C
  503             CONTINUE
C
               END IF
C
  400       CONTINUE
C
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
      END
