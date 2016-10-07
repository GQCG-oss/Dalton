C
C  /* Deck so_res_cbt */
      SUBROUTINE SO_RES_CBT(RES2E,LRES2E,RES2D,LRES2D,
     &                      DSRHF,LDSRHF,BTR1E,LBTR1E,
     &                      BTR1D,LBTR1D,BTJ1E,LBTJ1E,
     &                      BTJ1D,LBTJ1D,CMO,LCMO,
     &                      IDEL,ISDEL,ISYDIS,ISYMTR,WORK,LWORK)
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
#include "implicit.h"
#include "priunit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION RES2E(LRES2E), RES2D(LRES2D)
      DIMENSION DSRHF(LDSRHF), BTR1E(LBTR1E), BTR1D(LBTR1D)
      DIMENSION BTJ1E(LBTJ1E), BTJ1D(LBTJ1D), CMO(LCMO)
      DIMENSION WORK(LWORK)
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "soppinf.h"
CPi moved triplet variables to soppinf.h
C#include "infsop.h"
CPFP
c#include <inforb.h>
Cend-PFP
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
      SQ2   = DSQRT(TWO)
      ONESQ2= ONE/SQ2
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
               ISALFA = ISYMA
               ISBETA = MULD2H(ISALFA,ISALBE)
C
               LSCR4  = NBAS(ISALFA)*NRHF(ISYMI)
               LSCR5  = NBAS(ISALFA)*NRHF(ISYMI)
               LSCR6  = NBAS(ISALFA)*NRHF(ISYMI)
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
            CALL SO_MEMMAX ('SO_RES_CBT.4',LWORK4)
            IF (LWORK4 .LT. 0)
     &          CALL STOPIT('SO_RES_CBT.4',' ',KEND4,LWORK)
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
               FAC   = CMO(KOFF6)
C
               CALL DSCAL(NT1AM(ISYMAI),FAC,WORK(KSCR4),1)
               CALL DSCAL(NT1AM(ISYMAI),FAC,WORK(KSCR5),1)
C
CSPAS:26.04.2006: from here on I made modifications already
C
C------------------------------------------------------------
C                          ~ ~
C              Add SQ2 * ( a i | j b ) to 2p2h excitation and
C              de-excitation result-vectors.
C------------------------------------------------------------
C
               DO 500 ISYMI = 1, NSYM
C
                  ISYMA  = MULD2H(ISYMI,ISYMAI)
CSPASwas here
c                  IF (ISYMA .GE. ISYMB) THEN
CSPAShas left
                  ISYMAB = MULD2H(ISYMA,ISYMB)
                  ISYMIJ = MULD2H(ISYMI,ISYMJ)
C
                  IABT   = ITVV(MAX(ISYMA,ISYMB),MIN(ISYMA,ISYMB))
                  IABS   = ISVV(MAX(ISYMA,ISYMB),MIN(ISYMA,ISYMB))
                  IAB1   = IABT + B
                  IAB2   = IABS + B
C
                  IIJT   = ITOO(MAX(ISYMI,ISYMJ),MIN(ISYMI,ISYMJ))
                  IIJS   = ISOO(MAX(ISYMI,ISYMJ),MIN(ISYMI,ISYMJ))
                  IIJ1   = IIJT + J
                  IIJ3   = IIJS + J
C
                  DO 600 I = 1, NRHF(ISYMI)
C
                     NAI1 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA) * (I - 1)
C
                     IF (ISYMI .EQ. ISYMJ) THEN
                        NIJ1  = IIJT + INDLT(I,J)
                        NIJ2  = NIJ1
                        NIJ3  = IIJS + INDEX(I,J)
                        IF(I.EQ.J) THEN
                          SGNIJ1 = ZERO
                          SGNIJ2 = ZERO
                          SGNIJ3 = SQ2 ! Pi: In T3 part elimnate ONESQ2
                        ELSE IF(I.GT.J) THEN
                          SGNIJ1 = ONE
                          SGNIJ2 = ONE
                          SGNIJ3 = ONE
                        ELSE
                          SGNIJ1 = -ONE
                          SGNIJ2 = -ONE
                          SGNIJ3 = ONE
                        END IF
C
                     ELSE IF (ISYMI .GT. ISYMJ) THEN
                        KJOFF = NRHF(ISYMJ) * (I - 1)
                        NIJ1  = IIJ1 + KJOFF
                        NIJ2  = NIJ1
                        NIJ3  = IIJ3 + KJOFF
                        SGNIJ1= ONE
                        SGNIJ2= ONE
                        SGNIJ3= ONE
                     ELSE IF (ISYMI .LT. ISYMJ) THEN
                        KJOFF = NRHF(ISYMI) * (J - 1)
                        NIJ1  = IIJT + I + KJOFF
                        NIJ2  = NIJ1
                        NIJ3  = IIJS + I + KJOFF
                        SGNIJ1= -ONE
                        SGNIJ2= -ONE
                        SGNIJ3= ONE
                     END IF
C
                     IABIJ1 = IT2AMT1(ISYMIJ,ISYMAB)
     &                        + NTVV(ISYMAB) * (NIJ1 - 1)
                     IABIJ2 = IT2AMT2(ISYMIJ,ISYMAB)
     &                        + NSVV(ISYMAB) * (NIJ2 - 1)
                     IABIJ3 = IT2AMT3(ISYMIJ,ISYMAB)
     &                        + NTVV(ISYMAB) * (NIJ3 - 1)
C
                     DO 700 A = 1, NVIR(ISYMA)
C
                        NAI       = NAI1 + A
C
                        IF (ISYMA .EQ. ISYMB) THEN
                           NAB1  = IABT + INDLT(A,B)
                           NAB2  = IABS + INDEX(A,B)
                           NAB3  = NAB1
C
                           IF(A.EQ.B) THEN
                             SGN1 = ZERO
                             SGN2 = SGNIJ2*SQ2 ! CPi eliminate ONESQ2
                             SGN3 = ZERO
                           ELSE IF(A.GT.B) THEN
                             SGN1 = -SGNIJ1
                             SGN2 = SGNIJ2
                             SGN3 = SGNIJ3
                           ELSE
                             SGN1 = SGNIJ1
                             SGN2 = SGNIJ2
                             SGN3 = -SGNIJ3
                           END IF
C
                        ELSE IF (ISYMA .GT. ISYMB) THEN
C
                           KBOFF = NVIR(ISYMB) * (A - 1)
                           NAB1  = IAB1 + KBOFF
                           NAB2  = IAB2 + KBOFF
                           NAB3  = NAB1
                           SGN1  = -SGNIJ1
                           SGN2  = SGNIJ2
                           SGN3  = SGNIJ3
                        ELSE IF (ISYMA .LT. ISYMB) THEN
                           KBOFF = NVIR(ISYMA) * (B - 1)
                           NAB1  = IABT + A + KBOFF
                           NAB2  = IABS + A + KBOFF
                           NAB3  = NAB1
                           SGN1  = SGNIJ1
                           SGN2  = SGNIJ2
                           SGN3  = -SGNIJ3
                        END IF
C
                        NABIJ1    = IABIJ1 + NAB1
                        NABIJ2    = IABIJ2 + NAB2
                        NABIJ3    = IABIJ3 + NAB3
C
                 RES2E(NABIJ1) = RES2E(NABIJ1) + SGN1*WORK(KSCR4+NAI-1)
C     &                           *ONESQ2
C
                 RES2D(NABIJ1) = RES2D(NABIJ1) + SGN1*WORK(KSCR5+NAI-1)
C     &                           *ONESQ2
C
                 RES2E(NABIJ2) = RES2E(NABIJ2) + SGN2*WORK(KSCR4+NAI-1)
     &                           *ONESQ2
C
                 RES2D(NABIJ2) = RES2D(NABIJ2) + SGN2*WORK(KSCR5+NAI-1)
     &                           *ONESQ2
C
                 RES2E(NABIJ3) = RES2E(NABIJ3)
     &                         + SGN3*WORK(KSCR4+NAI-1)*ONESQ2
C
                 RES2D(NABIJ3) = RES2D(NABIJ3)
     &                         + SGN3*WORK(KSCR5+NAI-1)*ONESQ2
C
  700         CONTINUE
C
  600      CONTINUE
CSPASwas here
c                  END IF
CSPAShas left
C
  500   CONTINUE
C
CKeinSPASmehr
C
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
      CALL QEXIT('SO_RES_CBT')
C
      RETURN
      END
