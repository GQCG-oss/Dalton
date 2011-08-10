!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
#include <single.h>
C
C  /* Deck cclr_dia */
      SUBROUTINE CCLR_DIA(DIAA,ISYMTR,TRIPLET,WORK,LWORK)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Written by Ove Christiansen 4-9-1994
C     Triplet extensions by Christof Haettig & Kasper Hald, May 1999
C     R12 extension by Christof Haettig, Jun 2003
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      IMPLICIT NONE
#include <dummy.h>
#include <priunit.h>
#include <inftap.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
Cholesky
#include <maxorb.h>
#include <ccdeco.h>
Cholesky

      LOGICAL TRIPLET, LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER LWORK

#if defined (SYS_CRAY)
      REAL WORK(*), DIAA(*)
      REAL ONE, HALF, TMPKH
#else
      DOUBLE PRECISION WORK(*), DIAA(*)
      DOUBLE PRECISION ONE, HALF, TMPKH
#endif
      PARAMETER (ONE = 1.0D0, HALF = 0.5D00, TMPKH= 1.0D08)
C
      INTEGER KOFF1, KOFF2, KOFF2P, KOFF2M, KOFF12, KFOCKD, KEND1, 
     &        LWRK1, ISYMA, ISYMI, MI, MA, KAI, ISYMBJ, ISYMAI, 
     &        ISYMB, ISYMJ, MJ, MB, NAI, NBJ, NAIBJ, NAIAI, ISYMTR, 
     &        INDEX
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) -3 )/2 + I + J
C
      CALL QENTER('CCLR_DIA')

C     ----------------------------------------------------
C     set start addresses for the individual blocks parts:
C     (singles, doubles, R12 doubles)
C     ----------------------------------------------------
      KOFF1  = 1
      KOFF2  = KOFF1  + NT1AM(ISYMTR)
      IF (TRIPLET) THEN
        KOFF2P = KOFF2
        KOFF2M = KOFF2P + NT2AM(ISYMTR)
        KOFF12 = KOFF2M + NT2AMA(ISYMTR)
      ELSE
        KOFF12 = KOFF2  + NT2AM(ISYMTR)
      END IF
C
C-----------------------------
C     Allocation of workspace.
C-----------------------------
C
      KFOCKD  = 1
      KEND1   = KFOCKD  + NORBTS
      LWRK1   = LWORK   - KEND1
C
      IF ( LWRK1 .LT. 0 ) THEN
         CALL QUIT('Insufficient space in CCLR_DIA ')
      ENDIF
C
C-------------------------------------
C     Read canonical orbital energies.
C-------------------------------------
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBTS)
      CALL GPCLOSE(LUSIFC,'KEEP')
C
      IF (IPRINT.GT.100) THEN
         CALL AROUND('FOCK DIAGONAL ')
         CALL OUTPUT(WORK(KFOCKD),1,NORBTS,1,1,NORBTS,1,1,LUPRI)
      ENDIF
C
      IF (FROIMP .OR. FROEXP)
     *   CALL CCSD_DELFRO(WORK(KFOCKD),WORK(KEND1),LWRK1)

C------------------------------------------------
C     FOCK based  approximative  A(1,1) diagonal.
C------------------------------------------------
      DO ISYMA = 1,NSYM
         ISYMI = MULD2H(ISYMA,ISYMTR)
         DO I = 1,NRHF(ISYMI)
            MI = IORB(ISYMI) + I
            DO A = 1,NVIR(ISYMA)
                MA = IORB(ISYMA) + NRHF(ISYMA) +  A
                KAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
                DIAA(KAI) = WORK(MA) - WORK(MI)
            END DO
         END DO
      END DO
C
      IF (IPRINT.GT.20 .OR. DEBUG .OR. LOCDBG) THEN
         CALL AROUND('FOCK BASED DIAGONAL single exc. ')
         CALL OUTPUT(DIAA,1,NT1AM(ISYMTR),1,1,NT1AM(ISYMTR),1,1,
     &        LUPRI)
      ENDIF
C
      IF (CCS .OR. CIS) GOTO 9999
Cholesky
      IF (CHOINT .AND. CC2) THEN
         CALL QEXIT('CCLR_DIA')
         RETURN
      END IF
Cholesky
C
C------------------------------------------------
C     FOCK based approximative A(2,2) diagonal.
C------------------------------------------------
      DO 200 ISYMBJ = 1,NSYM
C
         ISYMAI = MULD2H(ISYMBJ,ISYMTR)
C
         DO 210 ISYMJ = 1,NSYM
C
            ISYMB = MULD2H(ISYMJ,ISYMBJ)
C
            DO 220 ISYMI = 1,NSYM
C
               ISYMA = MULD2H(ISYMI,ISYMAI)
C
               DO 230 J = 1,NRHF(ISYMJ)
C
                  MJ = IORB(ISYMJ) + J
C
                  DO 240 B = 1,NVIR(ISYMB)
C
                     NBJ = IT1AM(ISYMB,ISYMJ)
     *                + NVIR(ISYMB)*(J - 1) + B
C
                     MB = IORB(ISYMB) + NRHF(ISYMB) + B
C
                     DO 250 I = 1,NRHF(ISYMI)
C
                        MI = IORB(ISYMI) + I
C
                        DO 260 A = 1,NVIR(ISYMA)
C
                           NAI = IT1AM(ISYMA,ISYMI)
     *                      + NVIR(ISYMA)*(I - 1) + A
C
                           MA = IORB(ISYMA) + NRHF(ISYMA) +  A
C
                           IF (ISYMAI.EQ.ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                            + INDEX(NAI,NBJ)
                           ELSE IF (ISYMAI.LT.ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                            + NT1AM(ISYMAI)*(NBJ-1) + NAI
                           ELSE IF (ISYMAI.GT.ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMBJ,ISYMAI)
     *                            + NT1AM(ISYMBJ)*(NAI-1) + NBJ
                           ENDIF
C
                           DIAA(NAIBJ+NT1AM(ISYMTR)) =
     *                       WORK(MA)+WORK(MB)-WORK(MI)-WORK(MJ)
C
  260                   CONTINUE
  250                CONTINUE
  240             CONTINUE
  230          CONTINUE
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
C
      IF (TRIPLET) THEN
C        --------------------------------------------------------
C        2- block as for singlet vector, but with the ai,ai 
C        diagonal set to a huge value to ensure that the
C        corresponding elements in the trial vectors are zero.
C        --------------------------------------------------------
         CALL DCOPY(NT2AM(ISYMTR),DIAA(KOFF2P),1,DIAA(KOFF2M),1)
         IF (ISYMTR.EQ.1) THEN
          DO ISYMAI = 1,NSYM
            DO ISYMI = 1,NSYM
               ISYMA = MULD2H(ISYMI,ISYMAI)
               DO I = 1,NRHF(ISYMI)
                  MI = IORB(ISYMI) + I
                  DO A = 1,NVIR(ISYMA)
                     MA  = IORB(ISYMA) + NRHF(ISYMA) +  A
                     NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1)+A
                     NAIAI = IT2AM(ISYMAI,ISYMAI) + INDEX(NAI,NAI)
                     DIAA(KOFF2M-1+NAIAI) = 1.0D08
               END DO
               END DO
            END DO
          END DO
         END IF

C        -----------------------------------------------------------
C        2+ block as for singlet vector, but times a factor of 0.5
C        and with the elements a<=b or i<=j set to huge value
C        to ensure that the corresponding elements in the trial
C        vectors are zero 
C        -----------------------------------------------------------
         CALL DCOPY(NT2AM(ISYMTR),TMPKH,0,DIAA(KOFF2P),1)
         DO ISYMBJ = 1,NSYM
            ISYMAI = MULD2H(ISYMBJ,ISYMTR)
C
            IF (ISYMAI .LE. ISYMBJ) THEN
!
            DO ISYMI = 1,NSYM
            ISYMA = MULD2H(ISYMI,ISYMAI)
            DO ISYMJ = 1,NSYM
               ISYMB = MULD2H(ISYMJ,ISYMBJ)
C
               DO I = 1,NRHF(ISYMI)
               DO J = 1,NRHF(ISYMJ)
C
                  MI = IORB(ISYMI) + I
                  MJ = IORB(ISYMJ) + J
                  IF (MI .NE. MJ) THEN
C
                  DO A = 1,NVIR(ISYMA)
                  DO B = 1,NVIR(ISYMB)
C
                     MA = IORB(ISYMA) + NRHF(ISYMA) +  A
                     MB = IORB(ISYMB) + NRHF(ISYMB) + B
                     IF (MA .NE. MB) THEN
C
                     NBJ = IT1AM(ISYMB,ISYMJ)+NVIR(ISYMB)*(J-1)+B
                     NAI = IT1AM(ISYMA,ISYMI)+NVIR(ISYMA)*(I-1)+A
C
                        IF (ISYMAI.EQ.ISYMBJ) THEN
C
                          IF ((NAI .LT. NBJ) .AND.
     *                        (MA .LT. MB)) THEN
C
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                            + INDEX(NAI,NBJ)
                           DIAA(NT1AM(ISYMTR)+NAIBJ) =
     *                          DIAA(KOFF2M-1+NAIBJ)
!
                          ENDIF
                        ELSE IF (ISYMAI.LT.ISYMBJ) THEN
!
                          IF (((MA .LT. MB) .AND.
     *                         (MI .LT. MJ)) .OR.
     *                        ((MA .GT. MB) .AND.
     *                         (MI .GT. MJ))) THEN
!
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                            + NT1AM(ISYMAI)*(NBJ-1) + NAI
                              DIAA(NT1AM(ISYMTR)+NAIBJ) =
     *                             DIAA(KOFF2M-1+NAIBJ)
!
                          ENDIF
                        ENDIF

                     ENDIF
!
                  END DO  
                  END DO  
!
                  ENDIF
C
               END DO  
               END DO  
C
            END DO  
            END DO  
C
!
         ENDIF
!
         END DO  
C
      ENDIF
C
      IF (CCR12) THEN
C       CALL CCLR_DIAR12(DIAA(KOFF12),WORK,LWORK)
C       CALL QUIT('Missing subroutine in CCLR_DIA...')
c       Write(LUPRI,*)'Warning Missing subroutine in CCLR_DIA...'
c       Write(LUPRI,*)'Set R12 Jacobian diagonal to 1.0D01...'
        CALL DCOPY(NTR12AM(ISYMTR),1.0D01,0,DIAA(KOFF12),1)
      END IF
C
      IF (IPRINT.GT. 20 .OR. DEBUG .OR. LOCDBG) THEN
         CALL AROUND('APPROXIMATIVE  DIAGONAL 22 BLOCK ')
         IF (TRIPLET) THEN
           CALL OUTPUT(DIAA(1+NT1AM(ISYMTR)),1,2*NT2AM(ISYMTR),
     *                 1,1,2*NT2AM(ISYMTR),1,1,LUPRI)
         ELSE
           CALL OUTPUT(DIAA(1+NT1AM(ISYMTR)),1,1*NT2AM(ISYMTR),
     *                 1,1,1*NT2AM(ISYMTR),1,1,LUPRI)
         END IF
         IF (CCR12) THEN
           WRITE(LUPRI,*) 'R12 doubles part:'
           CALL OUTPUT(DIAA(KOFF12),1,NTR12AM(ISYMTR),
     *                 1,1,NTR12AM(ISYMTR),1,1,LUPRI)
         END IF
      ENDIF
C
 9999 CONTINUE
C
      CALL QEXIT('CCLR_DIA')
      RETURN
      END
C  /* Deck cclr_diascl */
      SUBROUTINE CCLR_DIASCL(OMEGA2,X,ISYMTR)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Written by Ove Christiansen 6-5-1994
C     symmetry 16-2-1995.
C
C     Purpose: Scale diagonal of omega2 with x
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
      DIMENSION OMEGA2(*)
#include <ccorb.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      IF (ISYMTR.EQ.1) THEN
C
         DO 50 ISYMBJ = 1, NSYM
C
            ISYMAI = MULD2H(ISYMBJ,ISYMTR)
C
            DO 100 NBJ = 1,NT1AM(ISYMBJ)
C
               NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NBJ,NBJ)
               OMEGA2(NAIBJ) = OMEGA2(NAIBJ)*X
C
  100       CONTINUE
C
   50    CONTINUE

      ENDIF
C
      RETURN
      END
C  /* Deck cclr_1c1f */
      SUBROUTINE CCLR_1C1F(RHO1,FOCKC,ISYMTR)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Ove Christiansen 19-Jan-1994
C
C     Purpose: Add Sum-k-Laikk-bar to rho1(ai)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
      PARAMETER (ONE = 1.0D0, ZERO = 0.0D0)
      DIMENSION RHO1(*),FOCKC(*)
#include <ccorb.h>
#include <ccsdsym.h>
C
C----------------------
C     Add contribution.
C----------------------
C
      ISYMAI = ISYMTR
C
      KOFF1 = 1
C
      DO 110 ISYMI = 1,NSYM
C
         ISYMA = MULD2H(ISYMI,ISYMAI)
C
         DO 120 I = 1, NRHF(ISYMI)
C
            K1 = KOFF1 + NORB(ISYMA)*(I-1) + NRHF(ISYMA)
            K2 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + 1
            CALL DAXPY(NVIR(ISYMA),ONE,FOCKC(K1),1,RHO1(K2),1)
C
  120    CONTINUE
C
         KOFF1 = KOFF1 + NORB(ISYMA)*NRHF(ISYMI)
C
  110 CONTINUE
C
      END
C  /* Deck cclr_e1c1 */
      SUBROUTINE CCLR_E1C1(RHO1,C1AM,EMAT1,EMAT2,WORK,LWORK,
     *                     ISYMV,ISYMIM,TRANS)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Written by Ove Christiansen 31-Jan-1995
C
C     input: Ematrices transformed in ccrhs_E
C     Purpose: Calculate E-terms in 1C1 part of linear transformation.
C
C     Calculates rho(ai) = Sum(c)C1AM(c,i)*E1(a,c) (If trans then E1(c,a)
C                        - Sum(k)C1AM(a,k)*E2(k,i) (If trans then E2(i,k)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
      PARAMETER (ONE = 1.0D0, XMONE=-1.0D0,TWO = 2.0D0)
      DIMENSION EMAT1(*), EMAT2(*),C1AM(*)
      DIMENSION WORK(LWORK),RHO1(*)
      CHARACTER*1 TRANS
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
C
C----------------------------------------
C     Calculate 1C1 contribution from E1.
C----------------------------------------
C
      ISYMAI = MULD2H(ISYMIM,ISYMV)
C
      DO 100 ISYMA = 1, NSYM
C
         ISYMI = MULD2H(ISYMAI,ISYMA)
         ISYMC = MULD2H(ISYMV ,ISYMI)
C
         NVIRA = MAX(NVIR(ISYMA),1)
         NVIRC = MAX(NVIR(ISYMC),1)
         KOFF2 = IT1AM(ISYMC,ISYMI) + 1
C
         KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
         IF (TRANS.EQ.'N') THEN
C
            KOFF1 = IMATAB(ISYMA,ISYMC) + 1
C
            CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMC),
     *                 ONE,EMAT1(KOFF1),NVIRA,C1AM(KOFF2),NVIRC,ONE,
     *                 RHO1(KOFF3),NVIRA)
C
         ELSE IF (TRANS.EQ.'T') THEN
C
            KOFF1 = IMATAB(ISYMC,ISYMA) + 1
C
            CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMC),
     *                 ONE,EMAT1(KOFF1),NVIRC,C1AM(KOFF2),NVIRC,ONE,
     *                 RHO1(KOFF3),NVIRA)
         ENDIF        
  100 CONTINUE
C
C----------------------------------------
C     Calculate 1C1 contribution from E2.
C----------------------------------------
C
      DO 200 ISYMA = 1, NSYM
C
         ISYMI = MULD2H(ISYMAI,ISYMA)
         ISYMK = MULD2H(ISYMV ,ISYMA)
C
         KOFF1 = IT1AM(ISYMA,ISYMK) + 1
         KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
         NVIRA = MAX(NVIR(ISYMA),1)
         NRHFK = MAX(NRHF(ISYMK),1)
         NRHFI = MAX(NRHF(ISYMI),1)
C
         IF (TRANS.EQ.'N') THEN
C
            KOFF2 = IMATIJ(ISYMK,ISYMI) + 1
C
            CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMK),
     *                 XMONE,C1AM(KOFF1),NVIRA,EMAT2(KOFF2),NRHFK,ONE,
     *                 RHO1(KOFF3),NVIRA)
C
         ELSE IF (TRANS.EQ.'T') THEN
C
            KOFF2 = IMATIJ(ISYMI,ISYMK) + 1
C
            CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMK),
     *                 XMONE,C1AM(KOFF1),NVIRA,EMAT2(KOFF2),NRHFI,ONE,
     *                 RHO1(KOFF3),NVIRA)
C
         ENDIF 
C
  200 CONTINUE
C
      RETURN
C
      END
c*DECK CCLR_JAC
      SUBROUTINE CCLR_JAC(ISIDE,WORK,LWORK,APROXR12)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Written by Ove Christiansen May/November 1994,
C     to calculate the Jacobian by finite difference and by
C     construction with transformation from the right on the CC jacobian.
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
#include <maxorb.h>
#include <iratdef.h>
#include <ccorb.h>
#include <aovec.h>
C
      LOGICAL   RSPIM2
      DIMENSION WORK(LWORK)
      CHARACTER*24 CBLOCK
      CHARACTER*(*) APROXR12
C
#include <ccsdinp.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <leinf.h>
#include <priunit.h>
C
      IF (IPRINT .GT. 10) THEN
         CALL AROUND(' START OF CCLR_JAC ')
      ENDIF
C
      IF (NSYM.GT.1) THEN
         WRITE(LUPRI,*) 
     *   'Finite difference Joacobian does only work without symmetry'
         CALL QUIT(
     *   'Finite difference Joacobian does only work without symmetry')
      ENDIF
C
C--------------------
C     Initialization.
C--------------------
C
      ISYMTR = 1
      IPRLE  = IPRINT
      IF (JACTST) THEN
         NC1VEC = NT1AM(ISYMTR)
         NC2VEC = NT2AM(ISYMTR)
c        IF (NT1AMX.GT.3) THEN
c           NC1VEC = 3
c        ENDIF
c        IF (NT2AMX.GT.5) THEN
c           NC2VEC = 5
c        ENDIF
      ENDIF
      IF (FDJAC ) THEN
         NC1VEC = NT1AM(ISYMTR)
         NC2VEC = NT2AM(ISYMTR)
      ENDIF
      IF (CCR12) THEN
        NTEMP  = NT1AM(ISYMTR) + NT2AM(ISYMTR) + NTR12AM(ISYMTR)
        NTEMP2 = NTEMP*(NC1VEC + NC2VEC + NTR12AM(ISYMTR))
      ELSE
        NTEMP  = NT1AM(ISYMTR) + NT2AM(ISYMTR)
        NTEMP2 = NTEMP*(NC1VEC + NC2VEC)
      END IF
C
      KEND1 = 1
      LWRK1 = LWORK
C
C------------------------------------------------------------
C     Calculate Jacobian by Transformation with unit vectors.
C     First elements in workspace contains the jacobian.
C------------------------------------------------------------
C
      IF (IPRINT.GT.10) WRITE(LUPRI,*) 'CCLR_JAC Workspace:',LWORK
C     
      IF (JACTST) THEN
         KJACAN  = KEND1
         KEND1   = KJACAN + NTEMP2
         LWRK1   = LWORK  - KEND1
         IF ( LWRK1 .LE. 0 ) THEN
           WRITE(LUPRI,*) 'out of memory:'
           WRITE(LUPRI,*) 'need:',KEND1
           WRITE(LUPRI,*) 'have:',LWORK
           CALL QUIT('Too little workspace in CCLR_JAC')
         END IF
         CALL AROUND( 'Calculation of analytical CC Jacobian')
         CALL CCLR_JACAN(NC1VEC,NC2VEC,ISIDE,WORK(KJACAN),LWORK,
     &                   APROXR12)
      ENDIF
C
C-------------------------------------------------------
C     Calculate Jacobian by finite difference.
C     both for jactst and jacexp.
C     First elements in workspace contains the jacobian.
C-------------------------------------------------------
C
      KFDJAC  = KEND1
      KEND2   = KFDJAC + NTEMP2
      LWRK2   = LWORK - KEND2
      IF ( LWRK2 .LE. 0 ) THEN
           WRITE(LUPRI,*) 'out of memory:'
           WRITE(LUPRI,*) 'need:',KEND1
           WRITE(LUPRI,*) 'have:',LWORK
           CALL QUIT('Too little workspace in CCLR_JAC-2')
      END IF
      CALL AROUND( 'Calculation of finite difference  CC Jacobian')
      RSPIM2 =  RSPIM
      RSPIM  = .FALSE.
      CALL CCLR_FDJAC(NC1VEC,NC2VEC,WORK(KFDJAC),LWRK1,APROXR12)
      RSPIM  = RSPIM2
C
C--------------------------------------------------------
C        calculate difference between the two jacobiants.
C--------------------------------------------------------
C
      IF (JACTST) THEN
         CALL DAXPY(NTEMP2,-1.0D00,WORK(KFDJAC),1,WORK(KJACAN),1)
         IDXMX = IDAMAX(NTEMP2,WORK(KJACAN),1) 
         IDXCOL = 1 + (IDXMX-1)/NTEMP
         IDXROW = IDXMX - ((IDXCOL-1)*NTEMP)
         IF (IDXCOL.LE.NC1VEC) THEN
            CBLOCK = 'singles block:'
         ELSE IF (IDXCOL.LE.(NC1VEC+NC2VEC)) THEN
            CBLOCK = 'doubles block:'
            IDXCOL = IDXCOL - NC1VEC
         ELSE IF (IDXCOL.LE.NCCVAR) THEN
            CBLOCK = 'R12 doubles block:'
            IDXCOL = IDXCOL - NC1VEC - NC2VEC
         ELSE
            CALL QUIT('Error in CCLR_JAC')
         END IF
         IF ( IPRINT .GE. 5) THEN
            CALL AROUND( 'DIFFERENCE. CC JACOBIAN - singles PART  ' )
            CALL OUTPUT(WORK(KJACAN),1,NTEMP,1,NC1VEC,NTEMP,NC1VEC,1,
     &           LUPRI)
            CALL AROUND( 'DIFFERENCE. CC JACOBIAN - doubles PART  ' )
            CALL OUTPUT(WORK(KJACAN+NTEMP*NC1VEC),1,NTEMP,1,NC2VEC,
     *                  NTEMP,NC2VEC,1,LUPRI)
            IF (CCR12) THEN
              CALL AROUND('DIFFERENCE. CC JACOBIAN - R12 doubles PART')
              CALL OUTPUT(WORK(KJACAN+NTEMP*(NC1VEC+NC2VEC)),
     *                    1,NTEMP,1,NTR12AM(1),NTEMP,NTR12AM(1),1,LUPRI)
            END IF
         ENDIF
         DIFNOR  = DDOT(NTEMP2,WORK(KJACAN),1,WORK(KJACAN),1)
         WRITE(LUPRI,'(/A,E20.10/)') ' The norm of the difference'//
     *        ' between fd and analytical jacobian is ',SQRT(DIFNOR)
         WRITE(LUPRI,'(/2A/,3I5,E20.10/)') 
     *        ' Largest element in difference is in ',
     *         CBLOCK, IDXMX, IDXCOL, IDXROW, WORK(KJACAN-1+IDXMX)
         WRITE(LUPRI,'(A/)') ' END OF JACTST TEST'
      ENDIF
C
      IF (IPRINT .GT. 10) THEN
         CALL AROUND(' END OF CCLR_JAC ')
      ENDIF
C
      RETURN
      END
c*DECK CCLR_JACAN
      SUBROUTINE CCLR_JACAN(NC1VEC,NC2VEC,ISIDE,WORK,LWORK,APROXR12)
C
C-------------------------------------------------------------------------
C     Test routine for calculating the CC Jacobian by Transformation of
C     unit vectors.
C     calculates the first nc1vec rows of the C1 blocks
C     and nc2vec of the rows of the c2 blocks.
C     Written by Ove Christiansen 26-5-1994
C-------------------------------------------------------------------------
C
#include <implicit.h>
#include <maxorb.h>
#include <iratdef.h>
#include <ccorb.h>
#include <aovec.h>
#include <ccsdinp.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <leinf.h>
#include <priunit.h>
C
      DIMENSION WORK(LWORK)
      PARAMETER (XHALF = 0.5D00,XMTWO = -2.0D00, DELTA = 1.0D-08)
      CHARACTER*(*) APROXR12
C
      IF (IPRINT .GT. 10) THEN
         CALL AROUND(' START OF CCLR_JACAN')
      ENDIF
C
C----------------------------
C     Work space allocations.
C----------------------------
C
      ISYMTR     = 1
      IF (CCR12) THEN
        NTAMP      = NT1AM(ISYMTR) + NT2AM(ISYMTR) + NTR12AM(ISYMTR)
        NTAMP2     = NTAMP*(NC1VEC + NC2VEC + NTR12AM(ISYMTR))
      ELSE
        NTAMP      = NT1AMX + NT2AMX
        NTAMP2     = NTAMP*(NC1VEC + NC2VEC )
      END IF

      KJACOBI    = 1
      KEND1      = 1 + NTAMP2
      IF (CCR12) THEN
        KMETRIC = KEND1
        KEND1   = KMETRIC + NTAMP2
      END IF
      LWRK1      = LWORK - KEND1
      IF (LWRK1 .LT. NT1AMX + NT2AMX) 
     *        CALL QUIT('INSUFFICIENT WORK SPACE IN CCLR_JACAN')
C
C--------------------------------------------------------------
C     Put up nc1vec C1 unit vectors and nc2vec C2 unit vectors.
C--------------------------------------------------------------
C
      KC1  = KJACOBI
      KC2  = KC1 + NTAMP*NC1VEC
      KC12 = KC2 + NTAMP*NC2VEC
      CALL CCLR_UNIT(NTAMP,NC1VEC,NC2VEC,
     &               WORK(KC1),WORK(KC2),WORK(KC12))

      IF (CCR12) THEN
        ! for R12 initialize Metric with the same unit vectors:
        CALL DCOPY(NTAMP2,WORK(KJACOBI),1,WORK(KMETRIC),1)
      END IF
C
C---------------------------------------------------------
C     Make the transformation of the unit vectors.
C     Transformed vectors are returned as first
C     elements in work.
C     Work follows immediately after vectors as it should.
C---------------------------------------------------------
C
      DO I = 1, NC1VEC + NC2VEC + NTR12AM(1)
         CALL DCOPY(NTAMP,WORK(KC1+NTAMP*(I-1)),1,WORK(KEND1),1)
         CALL CC_ATRR(0.0D0,ISYMOP,ISIDE,WORK(KEND1),LWRK1,
     &                .FALSE.,DUMMY,APROXR12)
         CALL DCOPY(NTAMP,WORK(KEND1),1,WORK(KC1+NTAMP*(I-1)),1)
  
         IF (CCR12) THEN  
           KS12AM = KEND1 + NTAMP
           KM12AM = KMETRIC + NT1AM(ISYMTR)+NT2AM(ISYMTR)+NTAMP*(I-1)
           CALL DCOPY(NTR12AM(1),WORK(KS12AM),1,WORK(KM12AM),1)
         END IF
      END DO
C
C----------------------------------------------------------------
C     Calculate norms to compare with finite difference jacobian.
C----------------------------------------------------------------
C
      XNJ = DDOT(NTAMP2,WORK(KJACOBI),1,WORK(KJACOBI),1)
C
      WRITE(LUPRI,*) '  '
      WRITE(LUPRI,*) ' NORM OF UNIT VEC. TRA. JAC.', SQRT(XNJ)
      WRITE(LUPRI,*) '  '
C
      CALL CCLR_NORMS(X11,X21,XR1,X12,X22,XR2,X1R,X2R,XRR,
     &                WORK(KJACOBI),NTAMP,NC1VEC,NC2VEC)
C
      WRITE(LUPRI,*) 'NORM OF 11 PART OF UNIT VECT. JAC. ', SQRT(X11)
      WRITE(LUPRI,*) 'NORM OF 21 PART OF UNIT VECT. JAC. ', SQRT(X21)
      WRITE(LUPRI,*) 'NORM OF R1 PART OF UNIT VECT. JAC. ', SQRT(XR1)
      WRITE(LUPRI,*) 'NORM OF 12 PART OF UNIT VECT. JAC. ', SQRT(X12)
      WRITE(LUPRI,*) 'NORM OF 22 PART OF UNIT VECT. JAC. ', SQRT(X22)
      WRITE(LUPRI,*) 'NORM OF R2 PART OF UNIT VECT. JAC. ', SQRT(XR2)
      WRITE(LUPRI,*) 'NORM OF 1R PART OF UNIT VECT. JAC. ', SQRT(X1R)
      WRITE(LUPRI,*) 'NORM OF 2R PART OF UNIT VECT. JAC. ', SQRT(X2R)
      WRITE(LUPRI,*) 'NORM OF RR PART OF UNIT VECT. JAC. ', SQRT(XRR)
C
C------------------------------------------------
C     Print the columns of the jacobian obtained.
C------------------------------------------------
C
      IF (IPRINT .GE. 5) THEN
         CALL AROUND( 'CC JACOBIANT I PRESUME - A*C1 PART ' )
         CALL OUTPUT(WORK(KJACOBI),1,NTAMP,1,NC1VEC,NTAMP,NC1VEC,1,
     &        LUPRI)
         CALL AROUND( 'CC JACOBIANT I PRESUME - A*C2 PART ' )
         CALL OUTPUT(WORK(KJACOBI+NTAMP*NC1VEC),1,NTAMP,1,NC2VEC,
     *               NTAMP,NC2VEC,1,LUPRI)
         IF ( CCR12) THEN
           CALL AROUND( 'CC JACOBIANT I PRESUME - A*C12 PART ' )
           CALL OUTPUT(WORK(KJACOBI+NTAMP*(NC1VEC+NC2VEC)),1,NTAMP,
     *                 1,NTR12AM(1),NTAMP,NTR12AM(1),1,LUPRI)
           CALL AROUND( 'CC-R12 METRIC') 
           CALL OUTPUT(WORK(KMETRIC),1,NTAMP,1,NC1VEC+NC2VEC+NTR12AM(1),
     *                 NTAMP,NC1VEC+NC2VEC+NTR12AM(1),1,LUPRI)
         END IF

         CALL AROUND( 'FULL ANALYTICAL CC JACOBIANT' )
         IF (ISIDE.EQ.-1) THEN
           CALL CCLR_TRANSSQ(WORK(KJACOBI),NTAMP)
         ENDIF
         if (.not. CCR12) then
         call output(WORK(KJACOBI),1,NTAMP,1,NC1VEC+NC2VEC,NTAMP,
     &               NC1VEC+NC2VEC,1,LUPRI)
         else
         call output(WORK(KJACOBI),1,NTAMP,1,NC1VEC+NC2VEC+NTR12AM(1),
     &               NTAMP,NC1VEC+NC2VEC+NTR12AM(1),1,LUPRI)
         end if

         CALL AROUND(' END OF CCLR_JACAN ')
      ENDIF
C
      RETURN
      END
C
c*DECK CCLR_UNIT
      SUBROUTINE CCLR_UNIT(NTAMP,NC1VEC,NC2VEC,C1AM,C2AM,C12AM)
C
C     Put up the first NC1VEC and NC2VEC unit vectors in the single and double
C     excitation space respectively.
C
#include <implicit.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <cclr.h>
C
      DIMENSION C1AM(NTAMP,NC1VEC),C2AM(NTAMP,NC2VEC)
      DIMENSION C12AM(NTAMP,*)
C
C--------------------------------------
C     Put up the required unit vectors.
C--------------------------------------
C
      DO J = 1, NC1VEC
        CALL DZERO(C1AM(1,J),NTAMP)
        C1AM(J,J) = 1.0D00
      END DO

      DO J = 1, NC2VEC
        CALL DZERO(C2AM(1,J),NTAMP)
        C2AM(J + NT1AMX,J) = 1.0D00
      END DO

      IF (CCR12) THEN
        DO J = 1, NTR12AM(1)
          CALL DZERO(C12AM(1,J),NTAMP)
          C12AM(J + NT1AMX+NT2AMX,J) = 1.0D00
        END DO
      END IF
C
      RETURN
      END
c*DECK CCLR_NORMS
      SUBROUTINE CCLR_NORMS(X11,X21,XR1,X12,X22,XR2,X1R,X2R,XRR,
     &                      XJAC,NTAMP,NC1VEC,NC2VEC)
#include <implicit.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <cclr.h>
C
      DIMENSION XJAC(NTAMP,*)
C
      X11 = 0.0D00
      X12 = 0.0D00
      X1R = 0.0D00
      X21 = 0.0D00
      X22 = 0.0D00
      X2R = 0.0D00
      XR1 = 0.0D00
      XR2 = 0.0D00
      XRR = 0.0D00
C
      DO 100 I = 1, NC1VEC
C
         X11 = X11 + DDOT(NT1AMX,XJAC(1,I),1,XJAC(1,I),1)
         X21 = X21 + DDOT(NT2AMX,XJAC(1 + NT1AMX,I),1,
     *                    XJAC(1 + NT1AMX,I),1)
         IF (CCR12) THEN
           XR1 = XR1 + DDOT(NTR12AM(1),XJAC(1+NT1AMX+NT2AMX,I),1,
     *                                XJAC(1+NT1AMX+NT2AMX,I),1)
         END IF
C
 100  CONTINUE
C
      DO 200 I = NC1VEC+1, NC1VEC + NC2VEC
C
         X12 = X12 + DDOT(NT1AMX,XJAC(1,I),1,XJAC(1,I),1)
         X22 = X22 + DDOT(NT2AMX,XJAC(1 + NT1AMX,I),1,
     *                    XJAC(1 + NT1AMX,I),1)
         IF (CCR12) THEN
           XR2 = XR2 + DDOT(NTR12AM(1),XJAC(1+NT1AMX+NT2AMX,I),1,
     *                                XJAC(1+NT1AMX+NT2AMX,I),1)
         END IF
C
 200  CONTINUE

      IF (CCR12) THEN
       DO I = NC1VEC+NC2VEC+1,NC1VEC+NC2VEC+NTR12AM(1)
         X1R = X1R + DDOT(NT1AMX,XJAC(1,I),1,XJAC(1,I),1)
         X2R = X2R + DDOT(NT2AMX,XJAC(1 + NT1AMX,I),1,
     *                           XJAC(1 + NT1AMX,I),1)
         XRR = XRR + DDOT(NTR12AM(1),XJAC(1+NT1AMX+NT2AMX,I),1,
     *                              XJAC(1+NT1AMX+NT2AMX,I),1)
       END DO
      END IF
C
      RETURN
      END
c*DECK CCLR_DX
      SUBROUTINE CCLR_DX(VEC,X,N)
C
C     Written by Ove Christiansen 26-5-1994
C
#include <implicit.h>
C
      DIMENSION VEC(*)
      DO 100 I = 1, N
         VEC(I) = X
100   CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
c*DECK CCLR_FDJAC
      SUBROUTINE CCLR_FDJAC(NC1VEC,NC2VEC,WORK,LWORK,APROXR12)
C-----------------------------------------------------------------------
C     Test routine for calculating the CC Jacobian by finite 
C     difference on the CC vector function.
C     IF JACTST then first elements in work contains the jacobian.
C     calculates the first nc1vec rows of the C1 blocks
C     and nc2vec of the rows of the c2 blocks.
C     Written by Ove Christiansen 26-5-1994
C     Changes for R12 by Christof Haettig 21-5-2003
C-----------------------------------------------------------------------
      IMPLICIT NONE
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
#include <dummy.h>
#include <second.h>
#include <ccsdinp.h>
#include <cclr.h>
#include <r12int.h>
#include <ccr12int.h>
 
      INTEGER NC1VEC, NC2VEC, LWORK
   
#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL XHALF, XMTWO, DELTA
      REAL TI,X11,X12,X21,X22,XNJ,X1R,XR1,X2R,XR2,XRR,DDOT
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION XHALF, XMTWO, DELTA, DELTAI
      DOUBLE PRECISION TI,X11,X12,X21,X22,XNJ,X1R,XR1,X2R,XR2,XRR,DDOT
#endif 
      PARAMETER (XHALF = 0.5D00,XMTWO = -2.0D00, 
C    &           DELTA = 1.0D-05, DELTAI = 1.0D00/DELTA)
     &           DELTA = 1.0D-07, DELTAI = 1.0D00/DELTA)
      CHARACTER*10 MODEL,MODELR12
      CHARACTER*8 FC12AM,FS12AM,FC2AM,FS2AM
      CHARACTER*(*) APROXR12
      PARAMETER (FC12AM='CCR_C12M',FS12AM='CCR_S12M')
      PARAMETER (FC2AM='CCR_C2AM',FS2AM='CCR_S2DM')

      INTEGER NTAMP,NTAMP2,KJACOBI,KJACCO,KRHO1,KRHO2,KRHO1D,KRHO2D,
     &        KC1AM,KC2AM,KEND1,LWRK1,KJACOBI2,IOPT,LUJAC, KJACOBIR,
     &        NAI, NBJ, NKI, NLJ, KCR12AM, KRHOR12, KRHOR12D, 
     &        LUNIT, LUFC12, LUFS12, LUMET,
     &        LUFS2,LUFC2,LUMET2,LUMET3,KC2AMPCK,KSCR1
C
      INTEGER INDEX
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      CALL QENTER('FDJAC')
 
      LUFC2  = -1
      LUFS2  = -1
      LUFS12 = -1
      LUFC12 = -1

      IF (OMEGSQ) CALL QUIT('CCLR_FDJAC incompatible with OMEGSQ.')
      IF (LCOR.OR.LSEC) 
     &   CALL QUIT('CCLR_FDJAC incompatible with FCORE and FSECON.')
      IF (NSYM.NE.1) CALL QUIT('CCLR_FDJAC must run in C1 symmetry.')
C
      IF (IPRINT.GT.5) THEN
         CALL AROUND( 'IN CCLR_FDJAC: MAKING FINITE DIFF. CC JACOBIANT')
         WRITE(LUPRI,*) 'Number of single excitations:',NC1VEC
         WRITE(LUPRI,*) 'Number of double excitations:',NC2VEC
         IF (CCR12) THEN
           WRITE(LUPRI,*) 'This is a CC-R12 calculation...'
           WRITE(LUPRI,*) 'Number of R12 double excitations:',NTR12AM(1)
         END IF
      ENDIF
C 
C--------------------------------------------------
C     Write the jacobian to disk if FDJAC is done.
C--------------------------------------------------
C
      IF (FDJAC ) THEN
         LUJAC = -1
         CALL GPOPEN(LUJAC,'CCLR_JAC','UNKNOWN',' ','UNFORMATTED',
     *               IDUMMY,.FALSE.)
         IF (IPRINT.GT.5) WRITE(LUPRI,'(3X,A,//)') 
     &        'FINITE DIFF. JACOBIAN WILL BE WRITTEN TO DISK'
         IF (CCR12) THEN
           LUMET = -1
           CALL GPOPEN(LUMET,'CCLR_MET','UNKNOWN',' ','UNFORMATTED',
     *                 IDUMMY,.FALSE.)
           IF (IPRINT.GT.5) WRITE(LUPRI,'(3X,A,//)') 
     &        'METRIC MATRIX WILL BE WRITTEN TO DISK'
           CALL WOPEN2(LUFC12,FC12AM,64,0)
           CALL WOPEN2(LUFS12,FS12AM,64,0)
           IF (IANR12.EQ.2) THEN
             LUMET2 = -1
             CALL GPOPEN(LUMET2,'CCLR_MET2','UNKNOWN',' ','UNFORMATTED',
     *                 IDUMMY,.FALSE.)
             LUMET3 = -1
             CALL GPOPEN(LUMET3,'CCLR_MET3','UNKNOWN',' ','UNFORMATTED',
     *                 IDUMMY,.FALSE.)
             CALL WOPEN2(LUFS2,FS2AM,64,0)
             CALL WOPEN2(LUFC2,FC2AM,64,0)
           END IF
         END IF
      ENDIF
C
C----------------------------
C     Work space allocations.
C----------------------------
C
      ISYMTR     = 1
      IF (CCR12) THEN
        NTAMP      = NT1AM(ISYMTR) + NT2AM(ISYMTR) + NTR12AM(1)
        NTAMP2     = NTAMP*(NC1VEC + NC2VEC + NTR12AM(1))
      ELSE
        NTAMP      = NT1AM(ISYMTR) + NT2AM(ISYMTR)
        NTAMP2     = NTAMP*(NC1VEC + NC2VEC)
      END IF

      KEND1 = 1

      IF (JACTST) THEN
         KJACOBI = KEND1
         KEND1   = KJACOBI  + NTAMP2
      ELSE
         KJACOBI = -999999
      ENDIF

      IF (FDJAC ) THEN
         KJACCO  = KEND1
         KEND1   = KJACCO   + NTAMP
      ELSE
         KJACCO  = -999999
      ENDIF

      KRHO1      = KEND1
      KRHO2      = KRHO1    + NT1AMX
      KC1AM      = KRHO2    + MAX(NT2AMX,NT2AO(ISYMOP),2*NT2ORT(ISYMOP))
      KRHO1D     = KC1AM    + NT1AMX
      KRHO2D     = KRHO1D   + NT1AM(ISYMOP)
      KC2AM      = KRHO2D   + MAX(NT2AMX,NT2AO(ISYMOP),2*NT2ORT(ISYMOP))
      KEND1      = KC2AM    + MAX(NT2SQ(ISYMOP),NT2R12(1))

      IF (CCR12) THEN
        KCR12AM  = KEND1 
        KRHOR12  = KCR12AM  + NTR12AM(1)
        KRHOR12D = KRHOR12  + NTR12AM(1)
        KEND1    = KRHOR12D + NTR12AM(1)
        IF (IANR12.EQ.2) THEN
          KC2AMPCK = KEND1
          KSCR1  = KC2AMPCK + NT2AM(ISYMTR)
          KEND1  = KSCR1 + MAX(NT2AM(ISYMTR),NTR12AM(ISYMTR))
        END IF
      ELSE
        KCR12AM  = -999999
        KRHOR12  = -999999
        KRHOR12D = -999999
      ENDIF 

      LWRK1      = LWORK    - KEND1

      IF (LWRK1.LT.0) THEN
         WRITE (LUPRI,*) 'Too little work space in cc_jacfd '
         WRITE (LUPRI,*) 'AVAILABLE: LWORK   =  ',LWORK
         WRITE (LUPRI,*) 'NEEDED (AT LEAST)  =  ',KEND1
         CALL QUIT('TOO LITTLE WORKSPACE IN CCLR_FDJAC ')
      ENDIF

      IF (IPRINT .GT. 15) THEN
         WRITE (LUPRI,*) ' IN CCLR_FDJAC: KJACOBI =  ',KJACOBI
         WRITE (LUPRI,*) ' IN CCLR_FDJAC: KJACCO  =  ',KJACCO
         WRITE (LUPRI,*) ' IN CCLR_FDJAC: KRHO1   =  ',KRHO1
         WRITE (LUPRI,*) ' IN CCLR_FDJAC: KRHO2   =  ',KRHO2
         WRITE (LUPRI,*) ' IN CCLR_FDJAC: KC1AM   =  ',KC1AM
         WRITE (LUPRI,*) ' IN CCLR_FDJAC: KRHO1D  =  ',KRHO1D
         WRITE (LUPRI,*) ' IN CCLR_FDJAC: KRHO2D  =  ',KRHO2D
         WRITE (LUPRI,*) ' IN CCLR_FDJAC: KC2AM   =  ',KC2AM
         WRITE (LUPRI,*) ' IN CCLR_FDJAC: KCR12AM =  ',KCR12AM
         WRITE (LUPRI,*) ' IN CCLR_FDJAC: KRHOR12 =  ',KRHOR12
         WRITE (LUPRI,*) ' IN CCLR_FDJAC: KRHOR12D=  ',KRHOR12D
         WRITE (LUPRI,*) ' IN CCLR_FDJAC: KEND1   =  ',KEND1
      ENDIF

      KJACOBI2   = KJACOBI  + NC1VEC*NTAMP
      KJACOBIR   = KJACOBI2 + NC2VEC*NTAMP
C
C---------------------
C     Initializations.
C---------------------
C
      IF (JACTST) CALL DZERO(WORK(KJACOBI),NTAMP2)
      IF (FDJAC ) CALL DZERO(WORK(KJACCO),NTAMP)

      X11 = 0.0D00
      X12 = 0.0D00
      X1R = 0.0D00
      X21 = 0.0D00
      X22 = 0.0D00
      X2R = 0.0D00
      XR1 = 0.0D00
      XR2 = 0.0D00
      XRR = 0.0D00
C
C------------------------------------------------
C     Read the CC reference amplitudes From disk.
C------------------------------------------------
C
      IOPT = 3
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KC1AM),WORK(KC2AM))
C
C--------------------------------------
C     Calculate reference omega vector.
C--------------------------------------
C
      RSPIM = .FALSE.
C
      IF (.FALSE.) THEN
      CALL CCRHSN(WORK(KRHO1D),WORK(KRHO2D),WORK(KC1AM),
     &            WORK(KC2AM),
     *            WORK(KEND1),LWRK1,APROXR12)
C
      CALL DCOPY(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1),1)
      CALL DCOPY(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2),1)
      IF (CCR12) THEN
        LUNIT = -1
        CALL GPOPEN(LUNIT,'CC_OMEGAR12','UNKNOWN',' ','UNFORMATTED',
     *                     IDUMMY,.FALSE.)
        READ(LUNIT) (WORK(KRHOR12+N), N = 0, NTR12AM(1)-1 )
        CALL GPCLOSE(LUNIT,'KEEP')
      END IF
C
      IF (IPRINT.GT.125) THEN
         CALL AROUND( 'RHO1' )
         CALL OUTPUT(WORK(KRHO1),1,NVIRT,1,NRHFT,NVIRT,NRHFT,1,LUPRI)
         CALL AROUND( 'RHO2' )
         CALL OUTPAK(WORK(KRHO2),NT1AMX,1,LUPRI)
         IF (CCR12) THEN
           CALL AROUND( 'RHO-R12' )
           CALL OUTPAK(WORK(KRHOR12),NMATKI(1),1,LUPRI)
         END IF
      ENDIF
      END IF
C
C---------------------------------------------
C     Calculate Jacobiant by finite difference.
C---------------------------------------------
C
      DO I = 1, NC1VEC
         TI   = SECOND()

C        Compute Omega(t - 0.5 delta)
         WORK(KC1AM +I -1) = WORK(KC1AM +I -1 ) - 0.5D0*DELTA

         CALL CCRHSN(WORK(KRHO1),WORK(KRHO2),WORK(KC1AM),
     &               WORK(KC2AM),
     *               WORK(KEND1),LWRK1,APROXR12)
         IF (CCR12) THEN
           LUNIT = -1
           CALL GPOPEN(LUNIT,'CC_OMEGAR12','UNKNOWN',' ','UNFORMATTED',
     &                     IDUMMY,.FALSE.)
           READ(LUNIT) (WORK(KRHOR12+N), N = 0, NTR12AM(1)-1 )
           CALL GPCLOSE(LUNIT,'KEEP')
         END IF

C        Compute Omega(t + 0.5 delta)
         WORK(KC1AM +I -1) = WORK(KC1AM +I -1 ) + DELTA

         CALL CCRHSN(WORK(KRHO1D),WORK(KRHO2D),
     &               WORK(KC1AM),WORK(KC2AM),
     *               WORK(KEND1),LWRK1,APROXR12)

C        Restore t amplitudes
         WORK(KC1AM +I -1) = WORK(KC1AM +I -1) - 0.5D0*DELTA

         IF (IPRINT.GT.125) THEN
           CALL AROUND( 'RHO1D ' )
           CALL OUTPUT(WORK(KRHO1D),1,NVIRT,1,NRHFT,NVIRT,NRHFT,1,LUPRI)
           CALL AROUND( 'RHO2D ' )
           CALL OUTPAK(WORK(KRHO2D),NT1AMX,1,LUPRI)
         ENDIF

C        Calculate [Omega(t+0.5delta)-(omega(t-0.5delta)]/delta
         CALL DAXPY(NT1AMX,-1.0D00,WORK(KRHO1),1,WORK(KRHO1D),1)
         CALL DAXPY(NT2AMX,-1.0D00,WORK(KRHO2),1,WORK(KRHO2D),1)
         CALL CCLR_DIASCL(WORK(KRHO2D),XHALF,ISYMTR)
         CALL DSCAL(NT1AMX,DELTAI,WORK(KRHO1D),1)
         CALL DSCAL(NT2AMX,DELTAI,WORK(KRHO2D),1)

         IF (CCR12) THEN
           LUNIT = -1
           CALL GPOPEN(LUNIT,'CC_OMEGAR12','UNKNOWN',' ','UNFORMATTED',
     &                     IDUMMY,.FALSE.)
           READ(LUNIT) (WORK(KRHOR12D+N), N = 0, NTR12AM(1)-1 )
           CALL GPCLOSE(LUNIT,'KEEP') 
           IF (IPRINT.GT.125) THEN
             CALL AROUND( 'RHO-R12' )
             CALL OUTPAK(WORK(KRHO2),NMATKI(1),1,LUPRI)
           END IF
           CALL DAXPY(NTR12AM(1),-1.0D0,WORK(KRHOR12),1,
     &                WORK(KRHOR12D),1)
           CALL CCLR_DIASCLR12(WORK(KRHOR12D),BRASCL,1)
           CALL DSCAL(NTR12AM(1),DELTAI,WORK(KRHOR12D),1)
         END IF

         IF (JACTST) THEN
            CALL DCOPY(NT1AMX,WORK(KRHO1D),1,
     *                 WORK(KJACOBI+NTAMP*(I-1)),1)
            CALL DCOPY(NT2AMX,WORK(KRHO2D),1,
     *                WORK(KJACOBI+NTAMP*(I-1)+NT1AMX),1)
            IF (CCR12) CALL DCOPY(NTR12AM(1),WORK(KRHOR12D),1,
     *              WORK(KJACOBI+NTAMP*(I-1)+NT1AMX+NT2AMX),1)
         ENDIF
         IF (FDJAC ) THEN
            CALL DCOPY(NT1AMX,WORK(KRHO1D),1,WORK(KJACCO),1)
            CALL DCOPY(NT2AMX,WORK(KRHO2D),1,WORK(KJACCO+NT1AMX),1)
            IF (CCR12) CALL DCOPY(NTR12AM(1),WORK(KRHOR12D),1,
     *                                WORK(KJACCO+NT1AMX+NT2AMX),1)
            WRITE(LUJAC) (WORK(KJACCO +J-1), J = 1, NTAMP)
         ENDIF
         X11 = X11 + DDOT(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1D),1)
         X21 = X21 + DDOT(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2D),1)
         IF (CCR12) 
     *    XR1 = XR1 + DDOT(NTR12AM(1),WORK(KRHOR12D),1,WORK(KRHOR12D),1)
C
         TI   = SECOND() - TI
         IF (IPRINT.GT.5 ) THEN
            WRITE(LUPRI,*) '  '
            WRITE(LUPRI,*) 'FDJAC ROW NR. ',I,' DONE IN ',TI,' SEC.'
         ENDIF
C
      END DO
C
C----------------------------------------------------------------
C     Loop over T2 amplitudes. Take care of diagonal t2 elements
C     is in a different convention in the energy code.
C     Factor 1/2 from right , and factor 2 from left.
C----------------------------------------------------------------
C
      DO NAI = 1, NT1AMX
        DO NBJ = 1, NAI
         I = INDEX(NAI,NBJ)
C
         IF (I.LE.NC2VEC) THEN
           TI   = SECOND()

C          Compute Omega(t - 0.5 delta)
           WORK(KC2AM + I -1) = WORK(KC2AM+I -1) - 0.5D0*DELTA

           CALL CCRHSN(WORK(KRHO1),WORK(KRHO2),WORK(KC1AM),
     *               WORK(KC2AM),WORK(KEND1),LWRK1,APROXR12)
           IF (CCR12) THEN
            LUNIT = -1
            CALL GPOPEN(LUNIT,'CC_OMEGAR12','UNKNOWN',' ','UNFORMATTED',
     &                      IDUMMY,.FALSE.)
            READ(LUNIT) (WORK(KRHOR12+N), N = 0, NTR12AM(1)-1 )
            CALL GPCLOSE(LUNIT,'KEEP') 
           END IF
 
C          Compute Omega(t + 0.5 delta)
           WORK(KC2AM + I -1) = WORK(KC2AM+I -1) + DELTA

           CALL CCRHSN(WORK(KRHO1D),WORK(KRHO2D),
     &               WORK(KC1AM),
     *               WORK(KC2AM),WORK(KEND1),LWRK1,APROXR12)
 
C          Restore t amplitudes
           WORK(KC2AM + I -1) = WORK(KC2AM+I -1) - 0.5D0*DELTA
C
           IF (IPRINT.GT.125) THEN
              CALL AROUND( 'RHO1D ' )
              CALL OUTPUT(WORK(KRHO1D),1,NVIRT,1,NRHFT,NVIRT,NRHFT,1,
     &             LUPRI)
              CALL AROUND( 'RHO2D ' )
              CALL OUTPAK(WORK(KRHO2D),NT1AMX,1,LUPRI)
           ENDIF
C
C          Calculate [Omega(t+0.5delta)-(omegat-0.5delta)]/delta
           CALL DAXPY(NT1AMX,-1.0D00,WORK(KRHO1),1,WORK(KRHO1D),1)
           CALL DAXPY(NT2AMX,-1.0D00,WORK(KRHO2),1,WORK(KRHO2D),1)
           CALL CCLR_DIASCL(WORK(KRHO2D),XHALF,ISYMTR)
           CALL DSCAL(NT1AMX,DELTAI,WORK(KRHO1D),1)
           CALL DSCAL(NT2AMX,DELTAI,WORK(KRHO2D),1)

           IF (CCR12) THEN
             LUNIT = -1
             CALL GPOPEN(LUNIT,'CC_OMEGAR12','UNKNOWN',' ',
     &                     'UNFORMATTED',IDUMMY,.FALSE.)
             READ(LUNIT) (WORK(KRHOR12D+N), N = 0, NTR12AM(1)-1 )
             CALL GPCLOSE(LUNIT,'KEEP') 
             IF (IPRINT.GT.125) THEN
               CALL AROUND( 'RHO-R12' )
               CALL OUTPAK(WORK(KRHO2),NMATKI(1),1,LUPRI)
             END IF
             CALL DAXPY(NTR12AM(1),-1.0D0,WORK(KRHOR12),1,
     *                                   WORK(KRHOR12D),1)
             CALL CCLR_DIASCLR12(WORK(KRHOR12D),BRASCL,1)
             CALL DSCAL(NTR12AM(1),DELTAI,WORK(KRHOR12D),1)
           END IF
C
           IF (NAI.EQ.NBJ) THEN
              CALL DSCAL(NT1AMX,2.0D00,WORK(KRHO1D),1)
              CALL DSCAL(NT2AMX,2.0D00,WORK(KRHO2D),1)
              IF (CCR12) CALL DSCAL(NTR12AM(1),2.0D00,WORK(KRHOR12D),1)
           ENDIF
C
           IF (JACTST) THEN
              CALL DCOPY(NT1AMX,WORK(KRHO1D),1,
     *                 WORK(KJACOBI2+NTAMP*(I-1)),1)
              CALL DCOPY(NT2AMX,WORK(KRHO2D),1,
     *                 WORK(KJACOBI2+NTAMP*(I-1)+NT1AMX),1)
              IF (CCR12) CALL DCOPY(NTR12AM(1),WORK(KRHOR12D),1,
     *                WORK(KJACOBI2+NTAMP*(I-1)+NT1AMX+NT2AMX),1)
           ENDIF
C
           IF (FDJAC ) THEN
              CALL DCOPY(NT1AMX,WORK(KRHO1D),1,WORK(KJACCO ),1)
              CALL DCOPY(NT2AMX,WORK(KRHO2D),1,
     *                  WORK(KJACCO+NT1AMX),1)
              IF (CCR12) CALL DCOPY(NTR12AM(1),WORK(KRHOR12D),1,
     *                  WORK(KJACCO+NT1AMX+NT2AMX),1)
              WRITE(LUJAC) (WORK(KJACCO +J-1), J = 1, NTAMP)
              IF (CCR12.AND.IANR12.EQ.2) THEN
c               make unit vector for this element (0 1 0) and
c               compute lower nondiagonal part for metric and
c               store it on file CCLR_MET3               
                CALL DZERO(WORK(KC2AMPCK),NT2AM(1))
                WORK(KC2AMPCK+I-1) = 1.0d0
                CALL CC_WVEC(LUFC2,FC2AM,NT2AM(1),NT2AM(1),1,
     &                       WORK(KC2AMPCK))
                CALL DZERO(WORK(KSCR1),NTR12AM(1))
                CALL CC_WVEC(LUFC12,FC12AM,NTR12AM(1),NTR12AM(1),1,
     &                       WORK(KSCR1))
                CALL CC_R12METRIC(1,BRASCL,KETSCL,WORK(KEND1),LWRK1,
     &                            FC2AM,LUFC2,FC12AM,LUFC12,FS12AM,
     &                            LUFS12,FS2AM,LUFS2,1,.FALSE.,DUMMY)
                CALL CC_RVEC(LUFS12,FS12AM,NTR12AM(1),NTR12AM(1),1,
     &                       WORK(KSCR1))
                WRITE(LUMET3) (WORK(KSCR1 +J-1), J = 1, NTR12AM(1))
              END IF
           ENDIF
C
           X12 = X12 + DDOT(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1D),1)
           X22 = X22 + DDOT(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2D),1)
           IF (CCR12)
     *      XR2=XR2 + DDOT(NTR12AM(1),WORK(KRHOR12D),1,WORK(KRHOR12D),1)

           TI   = SECOND() - TI
           IF (IPRINT.GT.5 ) THEN
              WRITE(LUPRI,*) '  '
              WRITE(LUPRI,*) 'FDJAC ROW NR. ',I+NT1AMX,
     *                  ' DONE IN ',TI,' SEC.'
           ENDIF
C
         ENDIF
C
        END DO
      END DO
C
C----------------------------------------------------------------
C     Loop over R12 doubles amplitudes. 
C       Take care of the diagonal elements, which are in a 
C       different convention than in the energy code:
C          factor 1/2 from right , and factor 2 from left.
C     (Analogous to the conventional doubles amplitudes t^ab_ij)
C----------------------------------------------------------------
C
      IF (CCR12) THEN
       ! read R12 doubles amplitudes from file
       IOPT = 32
       CALL CC_RDRSP('R0 ',0,1,IOPT,MODELR12,DUMMY,WORK(KCR12AM))

       ! make copy from which correct amplitudes can be 
       ! restored after finite difference calculation
       LUNIT = -1
       CALL GPOPEN(LUNIT,'CCR12_TEMP','UNKNOWN',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
       WRITE(LUNIT) (WORK(KCR12AM+N), N = 0, NTR12AM(1)-1 )
       CALL GPCLOSE(LUNIT,'KEEP') 

       DO NKI = 1, NMATKI(1)
        DO NLJ = 1, NKI
         I = INDEX(NKI,NLJ)
C
           TI   = SECOND()

C          Compute Omega(t - 0.5 delta)
           WORK(KCR12AM + I -1) = WORK(KCR12AM+I -1) - 0.5D0*DELTA
           IOPT = 32
           CALL CC_WRRSP('R0 ',0,1,IOPT,MODELR12,DUMMY,DUMMY,
     *                   WORK(KCR12AM),WORK(KEND1),LWRK1)

           CALL CCRHSN(WORK(KRHO1),WORK(KRHO2),WORK(KC1AM),
     *               WORK(KC2AM),WORK(KEND1),LWRK1,APROXR12)
           LUNIT = -1
           CALL GPOPEN(LUNIT,'CC_OMEGAR12','UNKNOWN',' ','UNFORMATTED',
     &             IDUMMY,.FALSE.)
           READ(LUNIT) (WORK(KRHOR12+N), N = 0, NTR12AM(1)-1 )
           CALL GPCLOSE(LUNIT,'KEEP') 

C          Compute Omega(t + 0.5 delta)
           WORK(KCR12AM + I -1) = WORK(KCR12AM+I -1) + DELTA
           IOPT = 32
           CALL CC_WRRSP('R0 ',0,1,IOPT,MODELR12,DUMMY,DUMMY,
     *                   WORK(KCR12AM),WORK(KEND1),LWRK1)

           CALL CCRHSN(WORK(KRHO1D),WORK(KRHO2D),
     &                 WORK(KC1AM),
     *                 WORK(KC2AM),WORK(KEND1),LWRK1,APROXR12)
 
C          Restore t amplitudes
           WORK(KCR12AM + I -1) = WORK(KCR12AM+I -1) - 0.5D0*DELTA
C
           IF (IPRINT.GT.125) THEN
              CALL AROUND( 'RHO1D ' )
              CALL OUTPUT(WORK(KRHO1D),1,NVIRT,1,NRHFT,NVIRT,NRHFT,1,
     &             LUPRI)
              CALL AROUND( 'RHO2D ' )
              CALL OUTPAK(WORK(KRHO2D),NT1AMX,1,LUPRI)
           ENDIF
 
C          Calculate [Omega(t+0.5delta)-(omegat-0.5delta)]/delta
           CALL DAXPY(NT1AMX,-1.0D00,WORK(KRHO1),1,WORK(KRHO1D),1)
           CALL DAXPY(NT2AMX,-1.0D00,WORK(KRHO2),1,WORK(KRHO2D),1)
           CALL CCLR_DIASCL(WORK(KRHO2D),XHALF,ISYMTR)
           CALL DSCAL(NT1AMX,DELTAI,WORK(KRHO1D),1)
           CALL DSCAL(NT2AMX,DELTAI,WORK(KRHO2D),1)

           LUNIT = -1
           CALL GPOPEN(LUNIT,'CC_OMEGAR12','UNKNOWN',' ','UNFORMATTED',
     &             IDUMMY,.FALSE.)
           READ(LUNIT) (WORK(KRHOR12D+N), N = 0, NTR12AM(1)-1 )
           CALL GPCLOSE(LUNIT,'KEEP') 
           IF (IPRINT.GT.125) THEN
             CALL AROUND( 'RHO-R12' )
             CALL OUTPAK(WORK(KRHO2),NMATKI(1),1,LUPRI)
           END IF
           CALL DAXPY(NTR12AM(1),-1.0D0,WORK(KRHOR12),1,
     *                                 WORK(KRHOR12D),1)
           CALL CCLR_DIASCLR12(WORK(KRHOR12D),BRASCL,1)
           CALL DSCAL(NTR12AM(1),DELTAI,WORK(KRHOR12D),1)
C
           IF (NKI.EQ.NLJ) THEN
              CALL DSCAL(NT1AMX,KETSCL,WORK(KRHO1D),1)
              CALL DSCAL(NT2AMX,KETSCL,WORK(KRHO2D),1)
              IF (CCR12) CALL DSCAL(NTR12AM(1),KETSCL,WORK(KRHOR12D),1)
           ENDIF
C
           IF (JACTST) THEN
              CALL DCOPY(NT1AMX,WORK(KRHO1D),1,
     *                 WORK(KJACOBIR+NTAMP*(I-1)),1)
              CALL DCOPY(NT2AMX,WORK(KRHO2D),1,
     *                 WORK(KJACOBIR+NTAMP*(I-1)+NT1AMX),1)
              IF (CCR12) CALL DCOPY(NTR12AM(1),WORK(KRHOR12D),1,
     *                WORK(KJACOBIR+NTAMP*(I-1)+NT1AMX+NT2AMX),1)
           ENDIF
C
           IF (FDJAC ) THEN
              CALL DCOPY(NT1AMX,WORK(KRHO1D),1,WORK(KJACCO ),1)
              CALL DCOPY(NT2AMX,WORK(KRHO2D),1,
     *                  WORK(KJACCO+NT1AMX),1)
              CALL DCOPY(NTR12AM(1),WORK(KRHOR12D),1,
     *                  WORK(KJACCO+NT1AMX+NT2AMX),1)
              WRITE(LUJAC) (WORK(KJACCO +J-1), J = 1, NTAMP)

              ! make unit vector for this element and compute
              ! S * R_12 and store on file CCLR_MET:
              CALL DZERO(WORK(KRHOR12D),NTR12AM(1))
              WORK(KRHOR12D + I -1) = 1.0D0
              CALL CC_WVEC(LUFC12,FC12AM,NTR12AM(1),NTR12AM(1),1,
     &                     WORK(KRHOR12D))
              IF (IANR12.EQ.2) THEN
                CALL DZERO(WORK(KC2AMPCK),NT2AM(1))
                CALL CC_WVEC(LUFC2,FC2AM,NT2AM(1),NT2AM(1),1,
     &                       WORK(KC2AMPCK))
              END IF
              CALL CC_R12METRIC(1,BRASCL,KETSCL,WORK(KEND1),LWRK1,
     &                          FC2AM,LUFC2,FC12AM,LUFC12,FS12AM,
     &                          LUFS12,FS2AM,LUFS2,1,.FALSE.,DUMMY)
              CALL CC_RVEC(LUFS12,FS12AM,NTR12AM(1),NTR12AM(1),1,
     &                     WORK(KRHOR12D))
              IF (IANR12.EQ.2) THEN
                CALL CC_RVEC(LUFS2,FS2AM,NT2AM(1),NT2AM(1),1,
     &                       WORK(KC2AMPCK))
c               write mixed S matrix on file
                WRITE(LUMET2) (WORK(KC2AMPCK +J-1), J = 1, NT2AM(1))
              END IF
        
              WRITE(LUMET) (WORK(KRHOR12D +J-1), J = 1, NTR12AM(1))
           ENDIF
C
           X1R = X1R +DDOT(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1D),1)
           X2R = X2R +DDOT(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2D),1)
           XRR = XRR +DDOT(NTR12AM(1),WORK(KRHOR12D),1,WORK(KRHOR12D),1)

           TI   = SECOND() - TI
           IF (IPRINT.GT.5 ) THEN
              WRITE(LUPRI,*) '  '
              WRITE(LUPRI,*) 'FDJAC ROW NR. ',I+NT1AMX+NT2AMX,
     *                  ' DONE IN ',TI,' SEC.'
           ENDIF
C
        END DO
       END DO

       ! restore R12 amplitudes
       LUNIT = -1
       CALL GPOPEN(LUNIT,'CCR12_TEMP','UNKNOWN',' ','UNFORMATTED',
     &         IDUMMY,.FALSE.)
       READ(LUNIT) (WORK(KCR12AM+N), N = 0, NTR12AM(1)-1 )
       CALL GPCLOSE(LUNIT,'DELETE') 
       IOPT = 32
       CALL CC_WRRSP('R0 ',0,1,IOPT,MODELR12,DUMMY,DUMMY,
     *               WORK(KCR12AM),WORK(KEND1),LWRK1) 
      END IF
C
      WRITE(LUPRI,*) '    '
      WRITE(LUPRI,*) '**  FINITE DIFF WITH DELTA ',DELTA, '**'
      WRITE(LUPRI,*) '    '
      IF (FDJAC) CALL GPCLOSE(LUJAC,'KEEP')
      IF (FDJAC .AND. CCR12) THEN
          CALL GPCLOSE(LUMET,'KEEP')
          CALL WCLOSE2(LUFC12,FC12AM,'DELETE')
          CALL WCLOSE2(LUFS12,FS12AM,'DELETE')
          IF (IANR12.EQ.2) THEN
            CALL WCLOSE2(LUFS2,FS2AM,'DELETE')
            CALL WCLOSE2(LUFC2,FC2AM,'DELETE')
            CALL GPCLOSE(LUMET2,'KEEP')
            CALL GPCLOSE(LUMET3,'KEEP')
          END IF
      END IF
      IF (IPRINT.GE.5 .AND. JACTST) THEN
         CALL AROUND( 'FINITE DIFF. CC JACOBIANT - singles PART  ' )
         CALL OUTPUT(WORK(KJACOBI),1,NTAMP,1,NC1VEC,NTAMP,NC1VEC,1,
     &        LUPRI)
         CALL AROUND( 'FINITE DIFF. CC JACOBIANT - doubles PART  ' )
         CALL OUTPUT(WORK(KJACOBI2),1,NTAMP,1,NC2VEC,
     *               NTAMP,NC2VEC,1,LUPRI)
         IF (CCR12) THEN
           CALL AROUND( 'FINITE DIFF. CC JACOBIANT - R12 doubles PART')
           CALL OUTPUT(WORK(KJACOBIR),1,NTAMP,1,NTR12AM(1),
     *                 NTAMP,NTR12AM(1),1,LUPRI)
         END IF
      ENDIF
      IF (JACTST) THEN
         XNJ = X11 + X12 + X1R + X21 + X22 + X2R + XR1 + XR2 + XRR
         WRITE(LUPRI,*) '  '
         WRITE(LUPRI,*) ' NORM OF FIN. DIFF. JAC.', SQRT(XNJ)
         WRITE(LUPRI,*) '  '
         WRITE(LUPRI,*) 'NORM OF 11 PART OF FIN. DIFF. JAC. ', SQRT(X11)
         WRITE(LUPRI,*) 'NORM OF 21 PART OF FIN. DIFF. JAC. ', SQRT(X21)
         WRITE(LUPRI,*) 'NORM OF R1 PART OF FIN. DIFF. JAC. ', SQRT(XR1)
         WRITE(LUPRI,*) 'NORM OF 12 PART OF FIN. DIFF. JAC. ', SQRT(X12)
         WRITE(LUPRI,*) 'NORM OF 22 PART OF FIN. DIFF. JAC. ', SQRT(X22)
         WRITE(LUPRI,*) 'NORM OF R2 PART OF FIN. DIFF. JAC. ', SQRT(XR2)
         WRITE(LUPRI,*) 'NORM OF 1R PART OF FIN. DIFF. JAC. ', SQRT(X1R)
         WRITE(LUPRI,*) 'NORM OF 2R PART OF FIN. DIFF. JAC. ', SQRT(X2R)
         WRITE(LUPRI,*) 'NORM OF RR PART OF FIN. DIFF. JAC. ', SQRT(XRR)
      ENDIF
C
C-----------------------------------
C     restore intermediates on disk:
C-----------------------------------
C
      IOPT = 3
      RSPIM = .TRUE.
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KC1AM),WORK(KC2AM))
      CALL CCRHSN(WORK(KRHO1D),WORK(KRHO2D),WORK(KC1AM),
     &            WORK(KC2AM),
     *            WORK(KEND1),LWRK1,APROXR12)
C
      IF (IPRINT .GT. 10) THEN
         CALL AROUND(' END OF CCLR_FDJAC ')
      ENDIF
C
      CALL QEXIT('FDJAC')
C
      RETURN
      END
C-----------------------------------------------------------------------
c*DECK CCLR_DUMTRR
      SUBROUTINE CCLR_DUMTRR(CVEC1,RHO1,RHO2,RHO12,S12AM,S2AM,
     &                       WORK,LWORK)
C-------------------------------------------------------------------
C     Written by Ove Christiansen 21-11-1994
C
C     Makes the transformation from vectors onto the jacobian
C     by using a finite difference jacobian read in from disk.
C
C-------------------------------------------------------------------
C
#include <implicit.h>
#include <dummy.h>
C
#include <priunit.h>
#include <ccorb.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <cclr.h>
#include <ccsdio.h>
#include <aovec.h>
#include <leinf.h>
#include <r12int.h>
      PARAMETER ( TWO = 2.0D00,XHALF=0.5D00 )
      CHARACTER*1 TRANS
      DIMENSION CVEC1(*),WORK(LWORK)
      DIMENSION RHO1(*),RHO2(*), RHO12(*), S12AM(*), S2AM(*)
C
      IF (IPRINT .GT. 10) THEN
         CALL AROUND(' START OF CCLR_DUMTRR ')
      ENDIF
C
C----------------------------------
C        Print transformed vectors.
C----------------------------------
C
      IF (IPRINT.GT.10) THEN
         CALL AROUND(' CCLR_DUMTRR: vector to be transformed ')
         CALL CC_PRP(CVEC1,CVEC1(1+NT1AM(ISYMTR)),ISYMTR,1,1)
         CALL AROUND('R12 double excitation part of vector:')
         CALL OUTPUT(CVEC1(1+NT1AM(ISYMTR)+NT2AM(ISYMTR)),1,
     *               NTR12AM(ISYMTR),1,1,NTR12AM(ISYMTR),1,1,LUPRI)
      ENDIF
C
      NTAMP = NT1AM(ISYMTR) + NT2AM(ISYMTR)
      IF (CCR12) NTAMP = NTAMP + NTR12AM(ISYMTR)
C
C------------------------
C     Open Jacobian file.
C------------------------
C
      LUJAC = -1
      CALL GPOPEN(LUJAC,'CCLR_JAC','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUJAC)
C
C---------------------------
C     Work space allocation.
C---------------------------
C
      KJACCR= 1
      KEND1 = KJACCR+ NTAMP
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0)
     &     CALL QUIT(' INSUFFICIENT WORK SPACE IN CCLR_DUMTRR ')
C
C---------------------------------------------------
C     Transform vectors by one column/row at a time.
C---------------------------------------------------
C
      IF (IPRINT.GT.110) THEN
         CALL AROUND( 'PRINTOUT FROM CCLR_DUMTRR')
      ENDIF
      DO 100 I = 1, NT1AM(ISYMTR)
C
         IF (NSIDE.EQ.-1) THEN
C
            READ(LUJAC) (WORK(KJACCR+K-1),K=1,NTAMP)
            IF (IPRINT.GT.110) THEN
               WRITE (LUPRI,*) 
     *             'THE ',I,'-th coulumn of the jacobian is '
               CALL OUTPUT(WORK(KJACCR),1,NTAMP,1,1,NTAMP,1,1,LUPRI)
            ENDIF
C
         ELSE IF (NSIDE.EQ.1) THEN
C
            CALL CCLR_JACCR(WORK(KJACCR),I,LUJAC,WORK(KEND1),LWRK1)
            IF (IPRINT.GT.110) THEN
               WRITE (LUPRI,*) 'THE ',I,'-th row of the jacobian is '
               CALL OUTPUT(WORK(KJACCR),1,NTAMP,1,1,NTAMP,1,1,LUPRI)
            ENDIF
C
         ENDIF
C
         RHO1(I) = DDOT(NTAMP,WORK(KJACCR),1,CVEC1,1)
C
 100  CONTINUE
C
      DO 200 I = 1, NT2AM(ISYMTR)
C
         IF (NSIDE.EQ.-1) THEN
C
            READ(LUJAC) (WORK(KJACCR+K-1),K=1,NTAMP)
            IF (IPRINT.GT.110) THEN
               WRITE (LUPRI,*) 'THE',I+NT1AM(ISYMTR),
     *                    '-th coulumn of the jacobian '
               CALL OUTPUT(WORK(KJACCR),1,NTAMP,1,1,NTAMP,1,1,LUPRI)
            ENDIF
C
         ELSE IF (NSIDE.EQ.1) THEN
C
            CALL CCLR_JACCR(WORK(KJACCR),I+NT1AMX,LUJAC,WORK(KEND1),
     *                      LWRK1)
            IF (IPRINT.GT.110) THEN
               WRITE (LUPRI,*) 'THE ',I+NT1AM(ISYMTR),
     *                    '-th row of the jacobian '
               CALL OUTPUT(WORK(KJACCR),1,NTAMP,1,1,NTAMP,1,1,LUPRI)
            ENDIF
C
         ENDIF
C
         RHO2(I) = DDOT(NTAMP,WORK(KJACCR),1,CVEC1,1)
C
 200  CONTINUE

      IF (CCR12) THEN
       DO I = 1, NTR12AM(ISYMTR)
         IF (NSIDE.EQ.-1) THEN
            READ(LUJAC) (WORK(KJACCR+K-1),K=1,NTAMP)
            IF (IPRINT.GT.110) THEN
               WRITE (LUPRI,*) 'THE',I+NT1AM(ISYMTR)+NT2AM(ISYMTR),
     *                    '-th coulumn of the jacobian '
               CALL OUTPUT(WORK(KJACCR),1,NTAMP,1,1,NTAMP,1,1,LUPRI)
            ENDIF
         ELSE IF (NSIDE.EQ.1) THEN
            CALL CCLR_JACCR(WORK(KJACCR),I+NT1AMX+NT2AMX,LUJAC,
     *                      WORK(KEND1),LWRK1)
            IF (IPRINT.GT.110) THEN
               WRITE (LUPRI,*) 'THE ',I+NT1AM(ISYMTR)+NT2AM(ISYMTR),
     *                    '-th row of the jacobian '
               CALL OUTPUT(WORK(KJACCR),1,NTAMP,1,1,NTAMP,1,1,LUPRI)
            ENDIF
         ENDIF
         RHO12(I) = DDOT(NTAMP,WORK(KJACCR),1,CVEC1,1)
       END DO
      
       KSMAT = 1
       KEND1 = KSMAT + NTR12AM(ISYMTR)*NTR12AM(ISYMTR)
       IF (IANR12.EQ.2) THEN
         KS2MAT = KEND1
         KEND1  = KS2MAT + NT2AM(ISYMTR)*NT2AM(ISYMTR)
       END IF
       LWRK1 = LWORK - KEND1
       IF (LWRK1 .LT. 0)
     &     CALL QUIT(' INSUFFICIENT WORK SPACE IN CCLR_DUMTRR ')

       LUMET = -1
       CALL GPOPEN(LUMET,'CCLR_MET','UNKNOWN',' ','UNFORMATTED',
     *             IDUMMY,.FALSE.)
       CALL DZERO(WORK(KSMAT),NTR12AM(ISYMTR)*NTR12AM(ISYMTR))
       DO I = 1, NTR12AM(1)
         IOFFS = KSMAT-1+NTR12AM(ISYMTR)*(I-1)
         READ(LUMET) (WORK(IOFFS+J),J=1,NTR12AM(ISYMTR))
       END DO
       CALL GPCLOSE(LUMET,'KEEP')

       IF (IANR12.EQ.2) THEN
         LUMET2 = -1
         CALL GPOPEN(LUMET2,'CCLR_MET2','UNKNOWN',' ','UNFORMATTED',
     *               IDUMMY,.FALSE.)
         CALL DZERO(WORK(KS2MAT),NT2AM(ISYMTR)*NT2AM(ISYMTR))
         DO I = 1, NT2AM(ISYMTR)
           IOFFS2 = KS2MAT -1+NT2AM(ISYMTR)*(I-1)
           READ(LUMET2)(WORK(IOFFS2+J),J=1,NT2AM(ISYMTR))
         END DO
         CALL GPCLOSE(LUMET2,'KEEP')
       END IF

       IF (NSIDE.EQ.1) THEN
         TRANS = 'N'
       ELSE IF (NSIDE.EQ.-1) THEN
         TRANS = 'T'
       ELSE
         CALL QUIT('Illegal NSIDE in CCLR_DUMTRR.')
       END IF
       KR12 = NT1AM(ISYMTR)+NT2AM(ISYMTR)+1
       CALL DGEMV(TRANS,NTR12AM(ISYMTR),NTR12AM(ISYMTR),1.0D0,
     *            WORK(KSMAT),NTR12AM(ISYMTR),CVEC1(KR12),1,0.0D0,
     *            S12AM,1)
       IF (IANR12.EQ.2) THEN
         KR12 = NT1AM(ISYMTR)+1
         CALL DGEMV(TRANS,NT2AM(ISYMTR),NT2AM(ISYMTR),1.0D0,
     &              WORK(KS2MAT),NT2AM(ISYMTR),CVEC1(KR12),1,0.0D0,
     &              S2AM,1)

       END IF
      END IF
C
C----------------------------------
C        Print transformed vectors.
C----------------------------------
C
      IF (IPRINT .GT. 10) THEN
         CALL AROUND(' CCLR_DUMTRR: transformed vectors ')
         CALL CC_PRP(RHO1,RHO2,ISYMTR,1,1)
         CALL AROUND('R12 double excitation part of vector:')
         CALL OUTPUT(RHO12,1,NTR12AM(ISYMTR),1,1,
     *               NTR12AM(ISYMTR),1,1,LUPRI)
         CALL AROUND('R12 double excitation part of metric:')
         CALL OUTPUT(S12AM,1,NTR12AM(ISYMTR),1,1,
     *               NTR12AM(ISYMTR),1,1,LUPRI)
      ENDIF
C
      IF (IPRINT .GT. 10) THEN
         CALL AROUND(' END OF CCLR_DUMTRR ')
      ENDIF
C
      CALL GPCLOSE(LUJAC,'KEEP')
C
      RETURN
      END
c*DECK CCLR_JACCR
      SUBROUTINE CCLR_JACCR(XJACCR,I,LU,WORK,LWORK)
C
C---------------------------------------------------------------
C     Reads the I-th row of the jacobian stored in the file
C     specified by unit number LU.
C----------------------------------------------------------------
C
#include <implicit.h>
C
#include <priunit.h>
#include <ccorb.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <cclr.h>
#include <ccsdio.h>
#include <aovec.h>
C
      DIMENSION XJACCR(*),WORK(LWORK)
C
      NTAMP = NT1AM(ISYMTR) + NT2AM(ISYMTR)
      IF (CCR12) NTAMP = NTAMP + NTR12AM(ISYMTR)
C
C---------------------------
C     Work space allocation.
C---------------------------
C
      KROW  = 1
      KEND1 = KROW  + NTAMP
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0)
     &     CALL QUIT('INSUFFICIENT WORK SPACE IN CCLR_JACCR')
C
C---------------------------------------------------
C     Read in vectors by one column at a time.
C---------------------------------------------------
C
      REWIND(LU)
C
      DO 100 J = 1, NTAMP
C
         READ(LU) (WORK(KROW +K -1),K= 1, NTAMP)
C
         XJACCR(J) = WORK(KROW + I -1)
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cclr_lamtra */
      SUBROUTINE CCLR_LAMTRA(XLAMDP,XLAMPC,XLAMDH,XLAMHC,C1AM,ISYMTR)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Written by Ove Christiansen 10-2-1995
C
C     PURPOSE:
C             transform lamda matrices wtih C1AM
C             debugged 10-2-1995
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
#include <iratdef.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, XMONE= -1.0D00)
      DIMENSION C1AM(*),XLAMDP(*),XLAMDH(*),XLAMPC(*),XLAMHC(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
C
      IF (IPRINT .GT.25) THEN
         CALL AROUND('IN CCLR_LAMTRA  ')
      ENDIF
C
C-----------------------------------------
C     Transform lamda particle matrix.
C     LaP~(al,a) = -sum(k)[ La(al,k)*C(a,k)]
C NB!! notice the minus sign.
C-----------------------------------------
C
      CALL DZERO(XLAMPC,NGLMDT(ISYMTR))
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMK  = MULD2H(ISYMTR,ISYMA)
         ISYMAL = ISYMK
C
         NBASAL = MAX(NBAS(ISYMAL),1)
         NBASA  = MAX(NVIR(ISYMA),1)
C
         KOFF1  = IGLMRH(ISYMAL,ISYMK) + 1
         KOFF2  = IT1AM(ISYMA,ISYMK) + 1
         KOFF3  = IGLMVI(ISYMAL,ISYMA) + 1
C
         CALL DGEMM('N','T',NBAS(ISYMAL),NVIR(ISYMA),NRHF(ISYMK),
     *              XMONE,XLAMDP(KOFF1),NBASAL,C1AM(KOFF2),NBASA,
     *              ZERO,XLAMPC(KOFF3),NBASAL)
C
  100 CONTINUE
C
C-----------------------------------------
C     Transform lamda hole matrix.
C     LaH~(al,i) = sum(c)[ La(al,c)*C(c,i)]
C-----------------------------------------
C
      CALL DZERO(XLAMHC,NGLMDT(ISYMTR))
C
      DO 200 ISYMI = 1,NSYM
C
         ISYMC  = MULD2H(ISYMTR,ISYMI)
         ISYMAL = ISYMC
C
         NBASAL = MAX(NBAS(ISYMAL),1)
         NBASC  = MAX(NVIR(ISYMC),1)
C
         KOFF1  = IGLMVI(ISYMAL,ISYMC) + 1
         KOFF2  = IT1AM(ISYMC,ISYMI) + 1
         KOFF3  = IGLMRH(ISYMAL,ISYMI) + 1
C
         CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMI),NVIR(ISYMC),
     *              ONE,XLAMDH(KOFF1),NBASAL,C1AM(KOFF2),NBASC,
     *              ZERO,XLAMHC(KOFF3),NBASAL)
C
  200 CONTINUE
C
      RETURN
      END
C  /* Deck cc_prtspv */
      SUBROUTINE CC_PRTSPV(T1AM,T2AM)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Ove Christiansen 24-1-1994
C     PRint Total Symmetric Packed Vector (PRTSPV) (T1,T2)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
      DIMENSION T1AM(NT1AMX),T2AM(NT2AMX)
C
C------------------------------------
C     Write single excitation vector.
C------------------------------------
C
      CALL AROUND('single excitation part of vector ')
      DO 300 ISYM = 1,NSYM
         WRITE(LUPRI,*) 'Symmetry block number : ',ISYM
           KOFF = IT1AM(ISYM,ISYM) + 1
           CALL OUTPUT(T1AM(KOFF),1,NVIR(ISYM),1,NRHF(ISYM),
     *                 NVIR(ISYM),NRHF(ISYM),1,LUPRI)
         WRITE(LUPRI,*) ' '
  300 CONTINUE
C
C------------------------------------
C     Write double excitation vector.
C------------------------------------
C
      CALL AROUND('double excitation part of vector ')
      DO 310 ISYM = 1,NSYM
         WRITE(LUPRI,*) 'Symmetry block number : ',ISYM
         KOFF = IT2AM(ISYM,ISYM) + 1
         CALL OUTPAK(T2AM(KOFF),NT1AM(ISYM),1,LUPRI)
         WRITE(LUPRI,*) ' '
  310 CONTINUE
C
      RETURN
      END
C  /* Deck cc_prp */
      SUBROUTINE CC_PRP(T1AM,T2AM,ISYM,NT1,NT2)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Ove Christiansen 24-1-1994
C     PRint Packed Vector - PRP (in general non. tot. sym.) (t1,t2).
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
C
      DIMENSION T1AM(*),T2AM(*)
C
C------------------------------------
C     set dimensions
C------------------------------------
C
      LT1AM = NT1AM(ISYM)
      LT2AM = NT2AM(ISYM)

C
C------------------------------------
C     Write single excitation vector.
C------------------------------------
C
      DO 10 I = 1, NT1
C
         CALL AROUND('single excitation part of vector ')
         IF (NT1.GT.1) WRITE(LUPRI,*) 'Vector Number:',I
C
         DO 100 ISYMI = 1,NSYM
C
            ISYMA = MULD2H(ISYMI,ISYM)
C
            WRITE(LUPRI,*) 'Symmetry block number(a,i): ',ISYMA,ISYMI
            KOFF = LT1AM*(I-1) + IT1AM(ISYMA,ISYMI) + 1
            CALL OUTPUT(T1AM(KOFF),1,NVIR(ISYMA),1,NRHF(ISYMI),
     *                  NVIR(ISYMA),NRHF(ISYMI),1,LUPRI)
            WRITE(LUPRI,*) ' '
C
  100    CONTINUE
C
   10 CONTINUE
C
C------------------------------------
C     Write double excitation vector.
C------------------------------------
C
      DO 20 I = 1,NT2
C
         CALL AROUND('double excitation part of vector ')
         IF (NT2.GT.1) WRITE(LUPRI,*) 'Vector Number:',I
C
         DO 200 ISYMBJ = 1,NSYM
C
            ISYMAI = MULD2H(ISYMBJ,ISYM)
C
            KOFF = LT2AM*(I-1) + IT2AM(ISYMAI,ISYMBJ) + 1

            IF (ISYMAI.EQ.ISYMBJ) THEN
C
               WRITE(LUPRI,*) 'Symmetry block number(ai,bj): ',
     &              ISYMAI,ISYMBJ
               CALL OUTPAK(T2AM(KOFF),NT1AM(ISYMAI),1,LUPRI)
               WRITE(LUPRI,*) ' '
C
            ELSE IF (ISYMBJ.GT.ISYMAI) THEN
C
               WRITE(LUPRI,*) 'Symmetry block number(ai,bj): ',
     &              ISYMAI,ISYMBJ
               CALL OUTPUT(T2AM(KOFF),1,NT1AM(ISYMAI),1,NT1AM(ISYMBJ),
     *                     NT1AM(ISYMAI),NT1AM(ISYMBJ),1,LUPRI)
               WRITE(LUPRI,*) ' '
C
            ENDIF
C
  200    CONTINUE
C
   20 CONTINUE
C
      RETURN
      END
C  /* Deck cc_prpao */
      SUBROUTINE CC_PRPAO(T1AO,T2AO,ISYM,NT1,NT2)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Ove Christiansen 24-1-1994
C     PRint Packed Vector - PRP (in general non. tot. sym.) (t1,t2).
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
C
      DIMENSION T1AO(*),T2AO(*)
C
C------------------------------------
C     set dimensions
C------------------------------------
C
      LT1AO = NT1AO(ISYM)
      LT2AO = NT2AO(ISYM)

C
C------------------------------------
C     Write single excitation vector.
C------------------------------------
C
      DO 10 I = 1, NT1
C
         CALL AROUND('single excitation part of vector ')
         IF (NT1.GT.1) WRITE(LUPRI,*) 'Vector Number:',I
C
         DO 100 ISYMI = 1,NSYM
C
            ISYMA = MULD2H(ISYMI,ISYM)
C
            WRITE(LUPRI,*) 'Symmetry block number(a,i): ',ISYMA,ISYMI
            KOFF = LT1AO*(I-1) + IT1AO(ISYMA,ISYMI) + 1
            CALL OUTPUT(T1AO(KOFF),1,NBAS(ISYMA),1,NRHF(ISYMI),
     *                  NBAS(ISYMA),NRHF(ISYMI),1,LUPRI)
            WRITE(LUPRI,*) ' '
C
  100    CONTINUE
C
   10 CONTINUE
C
C------------------------------------
C     Write double excitation vector.
C------------------------------------
C
      DO 20 I = 1,NT2
C
         CALL AROUND('double excitation part of vector ')
         IF (NT2.GT.1) WRITE(LUPRI,*) 'Vector Number:',I
C
         DO 200 ISYMBJ = 1,NSYM
C
            ISYMAI = MULD2H(ISYMBJ,ISYM)
C
            KOFF = LT2AO*(I-1) + IT2AO(ISYMAI,ISYMBJ) + 1

            IF (ISYMAI.EQ.ISYMBJ) THEN
C
               WRITE(LUPRI,*) 'Symmetry block number(ai,bj): ',
     &              ISYMAI,ISYMBJ
               CALL OUTPAK(T2AO(KOFF),NT1AO(ISYMAI),1,LUPRI)
               WRITE(LUPRI,*) ' '
C
            ELSE IF (ISYMBJ.GT.ISYMAI) THEN
C
               WRITE(LUPRI,*) 'Symmetry block number(ai,bj): ',
     &              ISYMAI,ISYMBJ
               CALL OUTPUT(T2AO(KOFF),1,NT1AO(ISYMAI),1,NT1AO(ISYMBJ),
     *                     NT1AO(ISYMAI),NT1AO(ISYMBJ),1,LUPRI)
               WRITE(LUPRI,*) ' '
C
            ENDIF
C
  200    CONTINUE
C
   20 CONTINUE
C
      RETURN
      END
C  /* Deck cc_prsq */
      SUBROUTINE CC_PRSQ(T1AM,T2AM,ISYM,NT1,NT2)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Ove Christiansen 24-1-1994
C     PRint SQuared vector. general non. tot. sym. (t1,t2).
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
C
      DIMENSION T1AM(*),T2AM(*)

C
C------------------------------------
C     Set dimensions                  
C------------------------------------
C
      LT1AM = NT1AM(ISYM)
      LT2AM = NT2SQ(ISYM)
C
C------------------------------------
C     Write single excitation vector.
C------------------------------------
C
      DO 10 I = 1, NT1
C
         CALL AROUND('single excitation part of vector ')
         IF (NT1.GT.1) WRITE(LUPRI,*) 'Vector Number:',I
C
         DO 100 ISYMI = 1,NSYM
C
            ISYMA = MULD2H(ISYMI,ISYM)
C
            WRITE(LUPRI,*) 'Symmetry block number(a,i): ',ISYMA,ISYMI

            KOFF = LT1AM*(I-1) + IT1AM(ISYMA,ISYMI) + 1
            CALL OUTPUT(T1AM(KOFF),1,NVIR(ISYMA),1,NRHF(ISYMI),
     *               NVIR(ISYMA),NRHF(ISYMI),1,LUPRI)
            WRITE(LUPRI,*) ' '
C
  100    CONTINUE
C
   10 CONTINUE
C
C------------------------------------
C     Write double excitation vector.
C------------------------------------
C
      DO 20 I = 1, NT2
C
         CALL AROUND('double excitation part of vector ')
         IF (NT2.GT.1) WRITE(LUPRI,*) 'Vector Number:',I
C
         DO 200 ISYMBJ = 1,NSYM
C
            ISYMAI = MULD2H(ISYMBJ,ISYM)
C
            WRITE(LUPRI,*) 'Symmetry block number(ai,bj): ',
     &           ISYMAI,ISYMBJ
C
            KOFF = LT2AM*(I-1) + IT2SQ(ISYMAI,ISYMBJ) + 1
            CALL OUTPUT(T2AM(KOFF),1,NT1AM(ISYMAI),1,NT1AM(ISYMBJ),
     *                  NT1AM(ISYMAI),NT1AM(ISYMBJ),1,LUPRI)
C
            WRITE(LUPRI,*) ' '
C
  200    CONTINUE
C
   20 CONTINUE
C
      RETURN
      END
c*DECK CCLR_FDEXCI
      SUBROUTINE CCLR_FDEXCI(WORK,LWORK)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Purpose: Calculation of Excitation energies from
C              explicit jacobian constructed from finite
C              difference.
C
C     Written by Ove Christiansen November 1994
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
#include <maxorb.h>
#include <ccorb.h>
#include <iratdef.h>
#include <codata.h>
#include <priunit.h>
#include <ccsdinp.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccsdio.h>
c
#include <r12int.h>
C
      LOGICAL OPNSTA, LOCDBG
      DIMENSION WORK(LWORK)
      PARAMETER (THRZER = 1.0D-99)
C
#include <leinf.h>
C
C---------------------------------------------
C     Header of Excitation Energy calculation.
C---------------------------------------------
C
      CALL QENTER('FDEXCI')
c
      CALL TIMER('START ',TIMEIN,TIMOUT)

      IF (IPRINT .GT. 0) THEN
         WRITE (LUPRI,'(A,/)') '  '
         WRITE (LUPRI,'(A,/)')
     *    '1 <<<<<<<<<< Output from CCLR_FDEXCI >>>>>>>>> '
         IF (CCR12) WRITE (LUPRI,'(A,/)') 'This is a R12 run.'
         WRITE (LUPRI,'(A,/)') '  '
      ENDIF
C
C---------------------------------
C     Initialization of variables.
C---------------------------------
C
      MATZ    = 1
      NREDH   = NT1AMX + NT2AMX
      NTAMP   = NT1AMX + NT2AMX
      IF (CCR12) THEN
        NREDH = NREDH + NTR12AM(1)
        NTAMP = NTAMP + NTR12AM(1)
      END IF
C
C----------------------------
C     Alloction of Workspace.
C----------------------------
C
      IF (NTAMP.GT.MAXRED .OR. NREDH.GT.MAXRED) THEN
        WRITE (LUPRI,*)'MAXRED, NTAMP, NREDH', MAXRED, NTAMP, NREDH
        CALL QUIT('MAXRED too small in CCLR_FDEXCI.')
      END IF
C
      KWI     = 1
      KAMAT   = KWI  + MAXRED
      KIV1    = KAMAT+ MAXRED*MAXRED
      KFV1    = KIV1 + MAXRED
      KEND    = KFV1 + MAXRED
      LEND    = LWORK - KEND
      KEIVAL  = KEND
      KSOLEQ  = KEIVAL + MAXRED
      KAMAT2  = KSOLEQ + MAXRED*MAXRED
      KWRK1   = KAMAT2 + MAXRED*MAXRED
      IF (CCR12) THEN
        KDENOM = KWRK1
        KSMAT  = KDENOM + MAXRED 
        KWRK1  = KSMAT  + MAXRED*MAXRED
        IF (IANR12.EQ.2) THEN
          KSEIVAL  = KWRK1
          KSWI     = KSEIVAL + MAXRED
          KSOL     = KSWI + MAXRED
          KSOLEQS  = KSOL + MAXRED*MAXRED
          KSCOPY   = KSOLEQS + MAXRED*MAXRED
          KWRK1    = KSCOPY + MAXRED*MAXRED
        END IF
      END IF
      LWRK1   = LWORK - KWRK1
      IF (LWRK1.LT. 0)
     &     CALL QUIT('TOO LITTLE WORK SPACE IN CCLR_FDEXCI ')
C
C---------------------
C     Readin Jacobian.
C---------------------
C
      CALL CCLR_JACIN(WORK(KSOLEQ),WORK(KSMAT))
C
      IF (IPRINT .GE. 10) THEN
         CALL AROUND( 'FIN. DIF. CC JACOBIANT - A*C1 PART ' )
         CALL OUTPUT(WORK(KSOLEQ),1,NTAMP,1,NT1AMX,MAXRED,NTAMP,1,LUPRI)
         CALL AROUND( 'FIN. DIF. CC JACOBIANT - A*C2 PART ' )
         CALL OUTPUT(WORK(KSOLEQ +MAXRED*NT1AMX),1,NTAMP,1,NT2AMX,
     *              MAXRED,NT2AMX,1,LUPRI)
         IF (CCR12) THEN
          CALL AROUND( 'FIN. DIF. CC JACOBIANT - A*CR12 PART ' )
          CALL OUTPUT(WORK(KSOLEQ +MAXRED*(NT1AMX+NT2AMX)),1,NTAMP,
     *               1,NTR12AM(1),MAXRED,NTR12AM(1),1,LUPRI)
          CALL AROUND( 'CC METRIC' )
          CALL OUTPUT(WORK(KSMAT),1,NTAMP,1,NTAMP,MAXRED,NTAMP,1,LUPRI)
         END IF
      ENDIF
C
C---------------------------
C     Solve for eigenvalues.
C---------------------------
C
      IF (CCR12) THEN
        CALL DZERO(WORK(KAMAT),MAXRED*MAXRED)
        CALL DCOPY(MAXRED*MAXRED,WORK(KSOLEQ),1,WORK(KAMAT),1)
        WRITE(LUPRI,'(/A)')'Jacobi matrix:' 
        CALL OUTPUT(WORK(KAMAT),1,MAXRED,1,NREDH,MAXRED,MAXRED,1,LUPRI)

C       lujacobi = 0
C       call gpopen(lujacobi,'JACOBIMAT','UNKNOWN',' ','FORMATTED',
C    &              IDUMMY,.FALSE.)

C       do i = 1, nrhf(1)
C       do a = 1, nvir(1)
C          nai = nvir(1)*(i-1) + a
C          do j = 1, nrhf(1)
C          do b = 1, nvir(1)
C            nbj = nvir(1)*(j-1) + b 
C            if (nbj.le.nai) then
C              naibj = nai*(nai-1)/2 + nbj
C              idxaibj = naibj + nvir(1)*nrhf(1)

C              do k = 1, nrhf(1)
CCN            do m = 1, nrhf(1)
C              do m = 1, nrhfb(1)
CCN              nmk = nrhf(1)*(k-1) + m
C                nmk = nrhfb(1)*(k-1) + m
C                  do l = 1, nrhf(1)
CCN                do n = 1, nrhf(1)
C                  do n = 1, nrhfb(1)
CCN                  nnl = nrhf(1)*(l-1) + n
C                    nnl = nrhfb(1)*(l-1) + n
C                    if (nnl.le.nmk) then
C                      nmknl = nmk*(nmk-1)/2 + nnl
C                      idxmknl = nmknl + nvir(1)*nrhf(1) +
C    &                         (nrhf(1)*nvir(1))*(nrhf(1)*nvir(1)+1)/2

C                      idx23 = maxred*(idxaibj-1) + idxmknl 
C                      idx32 = maxred*(idxmknl-1) + idxaibj 

C                      if ((dabs(work(kamat-1+idx23))+
C    &                      dabs(work(ksmat-1+idx23))+
C    &                      dabs(work(kamat-1+idx32))+
C    &                      dabs(work(ksmat-1+idx32))).gt.1.0d-8) then

C                        write(lujacobi,'(8i2,4e20.10)') 
C    &                     i,a,j,b,k,m,l,n,
C    &                      work(kamat-1+idx23),
C    &                      work(ksmat-1+idx23),
C    &                      work(kamat-1+idx32),
C    &                      work(ksmat-1+idx32)
C                       end if

C                    end if 
C                  end do
C                  end do
C               end do
C               end do
C            end if
C          end do
C          end do
C       end do
C       end do
          
C       call gpclose(lujacobi,'KEEP')
        
c       get idea of stability matrix
        call cc_r12stabmat(work(kamat),maxred,nt1am(1),nt2am(1),
     &                     ntr12am(1),work(kwrk1),lwrk1)
       
        write(lupri,*)'A(2,3) Block:'
        call output(work(kamat),nt1am(1)+1,
     &              nt1am(1)+nt2am(1),nt1am(1)+nt2am(1)+1,
     &              nt1am(1)+nt2am(1)+ntr12am(1),maxred,maxred,1,lupri)
        write(lupri,*)'A(3,2) Block:'
        call output(work(kamat),nt1am(1)+nt2am(1)+1,
     &              nt1am(1)+nt2am(1)+ntr12am(1),nt1am(1)+1,
     &              nt1am(1)+nt2am(1),maxred,maxred,1,lupri)
c
        IF (IANR12.EQ.2) THEN
          WRITE(LUPRI,*)'S matrix:' 
          CALL OUTPUT(WORK(KSMAT),1,NTAMP,1,NTAMP,MAXRED,NTAMP,1,LUPRI)
          CALL DCOPY(MAXRED*MAXRED,WORK(KSMAT),1,WORK(KSCOPY),1)
chf
          write(lupri,*)'S(2,3) Block:'
          call output(work(ksmat),nt1am(1)+1,
     &              nt1am(1)+nt2am(1),nt1am(1)+nt2am(1)+1,
     &              nt1am(1)+nt2am(1)+ntr12am(1),maxred,maxred,1,lupri)
          write(lupri,*)'S(3,2) Block:'
          call output(work(ksmat),nt1am(1)+nt2am(1)+1,
     &                nt1am(1)+nt2am(1)+ntr12am(1),nt1am(1)+1,
     &                nt1am(1)+nt2am(1),maxred,maxred,1,lupri)
c         get eigenvalues of S matrix to check if it is positive definite
c         caution rg routine destroys old S matrix
          CALL RG(MAXRED,NREDH,WORK(KSCOPY),WORK(KSEIVAL),WORK(KSWI),
     &            MATZ,WORK(KSOLEQS),WORK(KIV1),WORK(KFV1),IERR)
          WRITE(LUPRI,'(/A)')'Eigenvalues of S matrix:' 
          CALL OUTPUT(WORK(KSEIVAL),1,NREDH,1,1,NREDH,MAXRED,1,LUPRI)
        END IF
c
        CALL RGG(MAXRED,NREDH,WORK(KAMAT),WORK(KSMAT),WORK(KEIVAL),
     *           WORK(KWI),WORK(KDENOM),MATZ,WORK(KSOLEQ),IERR)
        DO I = 1, NREDH
          IF (ABS(WORK(KDENOM-1+I)).GT.THRZER) THEN
            WORK(KEIVAL-1+I)  = WORK(KEIVAL-1+I)/WORK(KDENOM-1+I)
            WORK(KWI-1+I)     = WORK(KWI-1+I)   /WORK(KDENOM-1+I)
          ELSE
            WORK(KEIVAL-1+I)  = 1.0D0/THRZER
            WORK(KWI-1+I)     = WORK(KWI-1+I)/THRZER
          END IF
        END DO
      ELSE
        CALL DCOPY(MAXRED*MAXRED,WORK(KSOLEQ),1,WORK(KAMAT),1)
        CALL DCOPY(MAXRED*MAXRED,WORK(KSOLEQ),1,WORK(KAMAT2),1)
        CALL RG(MAXRED,NREDH,WORK(KAMAT),WORK(KEIVAL),WORK(KWI),MATZ,
     *          WORK(KSOLEQ),WORK(KIV1),WORK(KFV1),IERR)
      END IF
C
      IF (IPRINT .GE. 70 .OR. IERR .NE. 0) THEN
         WRITE (LUPRI,'(/A)') ' EIGENVALUES real part:'
         WRITE (LUPRI,'(A)') ' before sort of eigenvalues'
         CALL OUTPUT(WORK(KEIVAL),1,NREDH,1,1,NREDH,MAXRED,1,LUPRI)
         WRITE (LUPRI,'(/A)')
     *   ' EIGENVALUES imaginary part:'
         WRITE (LUPRI,'(A)') ' before sort of eigenvalues'
         CALL OUTPUT(WORK(KWI),1,NREDH,1,1,NREDH,MAXRED,1,LUPRI)
         WRITE (LUPRI,'(/A)') ' EIGENVECTORS :'
         WRITE (LUPRI,'(A)') ' before sort of eigenvalues'
         CALL OUTPUT(WORK(KSOLEQ),1,NREDH,1,NREDH,MAXRED,MAXRED,1,LUPRI)
      END IF
C
      CALL RGORD(MAXRED,NREDH,WORK(KEIVAL),WORK(KWI),WORK(KSOLEQ),
     *           .FALSE.)
C
      IF (IPRINT .GE. 70 ) THEN
         WRITE (LUPRI,'(/A)') ' EIGENVECTORS :'
         WRITE (LUPRI,'(A)') ' After  sort of eigenvalues'
         CALL OUTPUT(WORK(KSOLEQ),1,NREDH,1,NREDH,MAXRED,MAXRED,1,LUPRI)
      END IF
      IF ( IERR.NE.0 ) THEN
         WRITE(LUPRI,'(/A,I5)')
     *   '  EIGENVALUE PROBLEM NOT CONVERGED IN RG, IERR =',IERR
         CALL QUIT(' CCLR_EXPJAC:  EIGENVALUE EQUATION NOT CONVERGED ')
      END IF
C
c     CALL RGTST(LUPRI,IPRINT,MAXRED,NREDH,WORK(KAMAT2),WORK(KSOLEQ),
c    *          WORK(KWRK1),LWRK1,INFO)
C
C----------------------------------
C    Write out Excitation energies.
C----------------------------------
C
      WRITE (LUPRI,'(//A/A/A//A/A/)')
     *'  CCSD excitation energies :',
     *' ============================',
     *' (conversion factor used: 1 au = 27.2113957 eV)',
     *' Excitation no.       Hartree               eV',
     *' --------------       -------               --'
      DO 400 I = 1,NREDH
         WRITE (LUPRI,'(I10,2F20.6)')
     *      I,WORK(KEIVAL-1+I),WORK(KEIVAL-1+I)*XTEV
  400 CONTINUE
C
C---------------------------------
C    Analysis of solution vectors.
C---------------------------------
C
      THRESH = 0.05
      MAXLIN = 100
C
      DO 500 I = 1,NREDH
         WRITE(LUPRI,'(//5X,A,I2)')
     *'Analysis of the Coupled Cluster Excitation Vector Number : ',I
         WRITE(LUPRI,'(5X,A)')
     *'-------------------------------------------------------------'
         WRITE(LUPRI,'(/15X,A,F10.5,A)')
     *   'Excitation Energy   :  ',WORK(KEIVAL+I-1)*XTEV,
     *   ' eV'
         CALL CC_PRAM(WORK(KSOLEQ + (I-1)*MAXRED),PT1,1,.FALSE.)
  500 CONTINUE
C
      CALL QEXIT('FDEXCI')
c
      RETURN
      END
c*DECK CCLR_JACIN
      SUBROUTINE CCLR_JACIN(FDJACO,SMATRIX)
C
C-------------------------------------------------------------------
C     Written by Ove Christiansen 21-11-1994
C
C     Reads the jacobian from disk and put it into JAC.
C
C-------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <cclr.h>
#include <ccsdio.h>
#include <aovec.h>
#include <cbirea.h>
#include <r12int.h>
C
      DIMENSION FDJACO(MAXRED,MAXRED)
      DIMENSION SMATRIX(MAXRED,MAXRED)
C
      NTAMP = NT1AMX + NT2AMX
      IF (LMULBS) NTAMP = NTAMP + NTR12AM(1)
C
C------------------------
C     Open Jacobian file.
C------------------------
C
      LUJAC = -1
      CALL GPOPEN(LUJAC,'CCLR_JAC','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUJAC)
C
      DO I = 1, NTAMP
         READ(LUJAC) (FDJACO(J,I),J=1,NTAMP)
      END DO
C
      CALL GPCLOSE(LUJAC,'KEEP')


      IF (CCR12) THEN
        LUMET = -1
        CALL GPOPEN(LUMET,'CCLR_MET','UNKNOWN',' ','UNFORMATTED',
     *              IDUMMY,.FALSE.)
        CALL DZERO(SMATRIX,MAXRED*MAXRED)
        DO I = 1, NT1AMX+NT2AMX
          SMATRIX(I,I) = 1.0D0
        END DO
chf
        IF (IANR12.EQ.2) THEN
          LUMET2 = -1
          CALL GPOPEN(LUMET2,'CCLR_MET2','UNKNOWN',' ','UNFORMATTED',
     *              IDUMMY,.FALSE.)
          LUMET3 = -1
          CALL GPOPEN(LUMET3,'CCLR_MET3','UNKNOWN',' ','UNFORMATTED',
     *              IDUMMY,.FALSE.)

          DO I = NT1AM(1)+NT2AM(1)+1, NT1AM(1)+NT2AM(1)+NTR12AM(1)
            READ(LUMET2)(SMATRIX(NT1AM(1)+J,I),J=1,NT2AM(1))
          END DO
c
          DO I = NT1AM(1)+1, NT1AM(1)+NT2AM(1)
            READ(LUMET3)(SMATRIX(NT1AM(1)+NT2AM(1)+J,I),J=1,NTR12AM(1))
          END DO
          CALL GPCLOSE(LUMET2,'KEEP')
          CALL GPCLOSE(LUMET3,'KEEP')
        END IF
c
        DO I = NT1AMX+NT2AMX+1, NT1AMX+NT2AMX+NTR12AM(1)
          READ(LUMET) (SMATRIX(NT1AMX+NT2AMX+J,I),J=1,NTR12AM(1))
        END DO
        CALL GPCLOSE(LUMET,'KEEP')
c       WRITE(LUPRI,*)'S-Matrix out of cclr_jacin'
c       CALL OUTPUT(SMATRIX,1,NT1AMX+NT2AMX+NTR12AM(1),1,
c    &              NT1AMX+NT2AMX+NTR12AM(1),MAXRED,MAXRED,1,LUPRI)
      END IF
C
      RETURN
C
      END
C  /* Deck cc_core */
      SUBROUTINE CC_CORE(T1AM,T2AM,ISYMTR)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     30-5-1995 Ove Christiansen
C
C     Purpose: To zero core and secondary part of
C              CC vector.
C              ISYMTR is symmetry of vector.
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
C
      PARAMETER ( ZERO= 0.0D0)
C
      DIMENSION T1AM(*),T2AM(*)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      IF (IPRINT .GT. 45) THEN
         CALL AROUND('START OF CC_CORE: (T1,T2) vector packed ')
         CALL CC_PRP(T1AM,T2AM,ISYMTR,1,1)
      ENDIF
C
      ISYMAI = MULD2H(ISYMTR,ISYMOP)
C
C-------------------------------------------
C     Loop through single excitation vector.
C-------------------------------------------
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMI = MULD2H(ISYMAI,ISYMA)
C
         IF (LCOR) THEN
C
            DO 110 I = 1,NRHF(ISYMI)
C
               IF (I .LE. ICOR(ISYMI)) THEN
                  KOFF = IT1AM(ISYMA,ISYMI)
     *                 + NVIR(ISYMA)*(I-1) + 1
                  CALL DZERO(T1AM(KOFF),NVIR(ISYMA))
               ENDIF
C
  110       CONTINUE
C
         ENDIF
         IF (LSEC) THEN
C
            DO 120 A=1,NVIR(ISYMA)
C
               IA = NVIR(ISYMA) - A + 1
C
               IF (IA .LE. ISEC(ISYMA)) THEN
                  DO 130 I = 1, NRHF(ISYMI)
                     NAI = IT1AM(ISYMA,ISYMI)
     *                    + NVIR(ISYMA)*(I-1) + A
                     T1AM(NAI) = ZERO
  130             CONTINUE
               ENDIF
C
  120       CONTINUE
C
         ENDIF
C
  100 CONTINUE
C
C--------------------------------------------
C     Loop through Doublee excitation vector.
C     If not ccs or ccp2
C--------------------------------------------
C
      IF ( CCS ) RETURN
C
      ISAIBJ = MULD2H(ISYMTR,ISYMOP)
C
      DO 200 ISYMAI = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMAI,ISAIBJ)
C
         DO 210 ISYMJ = 1,NSYM
C
            ISYMB = MULD2H(ISYMJ,ISYMBJ)
C
            DO 220 ISYMI = 1,NSYM
C
               ISYMA = MULD2H(ISYMI,ISYMAI)
C
               DO 230 J = 1,NRHF(ISYMJ)
C
                  DO 240 B = 1,NVIR(ISYMB)
C
                     NBJ = IT1AM(ISYMB,ISYMJ)
     *                   + NVIR(ISYMB)*(J - 1) + B
C
                     DO 250 I = 1,NRHF(ISYMI)
C
                        DO 260 A = 1,NVIR(ISYMA)
C
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I - 1) + A
C
                           IF (((ISYMAI.EQ.ISYMBJ).AND.
     *                         (NAI .GT. NBJ)).OR.(ISYMAI.GT.ISYMBJ))
     *                          GOTO 260
C
                           IF (ISYMAI.EQ.ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                             + INDEX(NAI,NBJ)
                           ELSE
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                            + NT1AM(ISYMAI)*(NBJ-1) + NAI
                           ENDIF
C
                           IF (LCOR) THEN
C
                              IF (I .LE. ICOR(ISYMI)) THEN
C
                                 T2AM(NAIBJ) = ZERO
C
                              ELSE IF (J .LE. ICOR(ISYMJ)) THEN
C
                                 T2AM(NAIBJ) = ZERO
C
                              ENDIF
C
                           ENDIF
C
                           IF (LSEC) THEN
C
                              IA = NVIR(ISYMA) - A + 1
                              IB = NVIR(ISYMB) - B + 1
                              IF (IA .LE. ISEC(ISYMA)) THEN
C
                                 T2AM(NAIBJ) = ZERO
C
                              ELSE IF (IB .LE. ISEC(ISYMB)) THEN
C
                                 T2AM(NAIBJ) = ZERO
C
                              ENDIF
C
                           ENDIF
C
  260                   CONTINUE
  250                CONTINUE
  240             CONTINUE
  230          CONTINUE
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
C
      IF (IPRINT .GT. 45) THEN
         CALL AROUND('END OF CC_CORE: (T1,T2) vector packed ')
         CALL CC_PRP(T1AM,T2AM,ISYMTR,1,1)
      ENDIF
C
      RETURN
      END
C  /* Deck cc_pram */
      SUBROUTINE CC_PRAM(CAM,PT1,ISYMTR,LGRS)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     30-5-1995 Ove Christiansen
C
C     Purpose: Writes out vector:
C              %T1 and %T2 and ||T1||/||T2||
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
C
      PARAMETER (TWO = 2.0D0, THPRT = 1.0D-9, HUNDRED = 100.0D0)
C
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
Cholesky
#include <maxorb.h>
#include <ccdeco.h>
C
      LOGICAL CCSEFF
Cholesky
C
C
      LOGICAL LGRS
      DIMENSION CAM(*)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
Cholesky
      CCSEFF = CCS .OR. (CHOINT.AND.CC2)
Cholesky
C
C------------------------
C     Add up the vectors.
C------------------------
C
      C1NOSQ = 0.0D0
      C2NOSQ = 0.0D0
      KC1 = 1
      KC2 = KC1 + NT1AM(ISYMTR)
      C1NOSQ = DDOT(NT1AM(ISYMTR),CAM(KC1),1,CAM(KC1),1)
Chol  IF (.NOT. CCS) C2NOSQ = DDOT(NT2AM(ISYMTR),CAM(KC2),1,CAM(KC2),1)
      IF (.NOT. CCSEFF)
     &   C2NOSQ = DDOT(NT2AM(ISYMTR),CAM(KC2),1,CAM(KC2),1)
      CNOSQ  = C1NOSQ + C2NOSQ
C
      IF (.NOT. (CCSEFF.OR.CCP2) .AND. CNOSQ.NE.0.0D0) THEN
C
         WRITE(LUPRI,'(//10X,A)')
     *     'CC_PRAM:Overall Contribution of the Different Components'
         WRITE(LUPRI,'(10X,A//)')
     *     '--------------------------------------------------------'
         WRITE(LUPRI,'(/10X,A,10X,F10.4,A)')
     *     'Single Excitation Contribution : ',
     *     (C1NOSQ/CNOSQ)*HUNDRED,' %'
         WRITE(LUPRI,'(/10X,A,10X,F10.4,A)')
     *     'Double Excitation Contribution : ',
     *     (C2NOSQ/CNOSQ)*HUNDRED,' %'
         WRITE(LUPRI,'(/10X,A,10X,F10.4)')
     *     '||T1||/||T2||                  : ',
     *      SQRT(C1NOSQ/C2NOSQ)
         IF (LGRS) WRITE(LUPRI,'(/10X,A,10X,F10.4)')
     *     'Tau1 diagnostic                : ',
     *      SQRT(C1NOSQ/(TWO*DFLOAT(NRHFT)))
         PT1 = (C1NOSQ/CNOSQ)*HUNDRED
      ELSE
         PT1 = HUNDRED
      ENDIF
      WRITE(LUPRI,'(/10X,A,10X,F10.4)')
     *  'Norm of Total Amplitude Vector : ',SQRT(CNOSQ)
C
      CALL FLSHFO(LUPRI)
C
C----------------------------------------------
C     Initialize threshold etc from Printlevel.
C----------------------------------------------
C
      NL = MAX(1,2*IPRINT)
C
      CNOSQ = MAX(CNOSQ,THPRT)
C
      THR1 = SQRT(C1NOSQ/CNOSQ)/NL
      THR2 = SQRT(C2NOSQ/CNOSQ)/NL
      THR1 = MAX(THR1,1.0D-02)
      THR2 = MAX(THR2,1.0D-02)
      SUMOFP = 0.0D00
C
      IF (DEBUG) THR1 = 0.0D0
C
C-------------------------------------
C     Loop until a few is Printed out.
C-------------------------------------
C
C
C---------------------------------------
C     Loop through One excitation part.
C---------------------------------------
C
      WRITE(LUPRI,'(//A)')
     *     ' +=============================================='
     *    //'===============================+'
      WRITE(LUPRI,'(1X,A)')
     *     '| symmetry|  orbital index  |   Excitation Numbers'
     *     //'             |   Amplitude  |'
      WRITE(LUPRI,'(1X,A)')
     *     '|  Index  |   a   b   i   j |      NAI      NBJ |'
     *     //'     NAIBJ    |              |'
      WRITE(LUPRI,'(A)')
     *     ' +=============================================='
     *    //'===============================+'
C
      ISYMAI = MULD2H(ISYMTR,ISYMOP)
C
  1   CONTINUE
      N1 = 0
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMI = MULD2H(ISYMAI,ISYMA)
C
         DO 110 I = 1,NRHF(ISYMI)
C
            MI = IORB(ISYMI) + I
C
            DO 120 A=1,NVIR(ISYMA)
C
               NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
C
               MA = IORB(ISYMA) + NRHF(ISYMA) +  A
C
               IF (ABS(CAM(NAI)) .GE. THR1 ) THEN
C
                  WRITE(LUPRI,9990) ISYMA,ISYMI,A,I,NAI,CAM(NAI)
C
                  N1 = N1 + 1
                  SUMOFP = SUMOFP + CAM(NAI)*CAM(NAI)
C
               ENDIF
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      IF ((N1.LT.1).AND.(SQRT(C1NOSQ/CNOSQ).GT.1.0D-3)) THEN
         THR1 = THR1/5.0D0
         GOTO 1
      ENDIF
C
      CALL FLSHFO(LUPRI)
C
C--------------------------------------------
C     Loop through Double excitation vector.
C     If not ccs or ccp2
C--------------------------------------------
C
      IF (.NOT. ( CCSEFF .OR. CCP2 )) THEN
C
      WRITE(LUPRI,'(A)')
     *     ' +----------------------------------------------'
     *    //'-------------------------------+'
C

 2    CONTINUE
      N2 = 0
C
      DO 200 ISYMAI = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMAI,ISYMTR)
C
         DO 210 ISYMJ = 1,NSYM
C
            ISYMB = MULD2H(ISYMJ,ISYMBJ)
C
            DO 220 ISYMI = 1,NSYM
C
               ISYMA = MULD2H(ISYMI,ISYMAI)
C
               DO 230 J = 1,NRHF(ISYMJ)
C
                  MJ = IORB(ISYMJ) + J
C
                  DO 240 B = 1,NVIR(ISYMB)
C
                     NBJ = IT1AM(ISYMB,ISYMJ)
     *                   + NVIR(ISYMB)*(J - 1) + B
C
                     MB = IORB(ISYMB) + NRHF(ISYMB) + B
C
                     DO 250 I = 1,NRHF(ISYMI)
C
                        MI = IORB(ISYMI) + I
C
                        DO 260 A = 1,NVIR(ISYMA)
C
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I - 1) + A
C
                           MA = IORB(ISYMA) + NRHF(ISYMA) +  A
C
                           IF (((ISYMAI.EQ.ISYMBJ).AND.
     *                         (NAI .LT. NBJ)).OR.(ISYMAI.LT.ISYMBJ))
     *                          GOTO 260
C
                           IF (ISYMAI.EQ.ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                             + INDEX(NAI,NBJ)
                           ELSE
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                            + NT1AM(ISYMAI)*(NBJ-1) + NAI
                           ENDIF
C
                           KAIBJ = NAIBJ + NT1AM(ISYMTR)
                           IF (ABS(CAM(KAIBJ)) .GT. THR2 ) THEN
C
                              WRITE(LUPRI,9991) ISYMA,ISYMB,ISYMI,ISYMJ,
     *                                      A,B,I,J,NAI,NBJ,NAIBJ,
     *                                      CAM(KAIBJ)
                              N2 = N2 + 1
C
                              SUMOFP = SUMOFP + CAM(KAIBJ)*CAM(KAIBJ)
C
                           ENDIF
C
  260                   CONTINUE
  250                CONTINUE
  240             CONTINUE
  230          CONTINUE
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
C
      IF ((N2.LT.1).AND.(SQRT(C2NOSQ/CNOSQ).GT.1.0D-3)) THEN
         THR2 = THR2/5D00
         GOTO 2
      ENDIF
C
      ENDIF
C
      WRITE(LUPRI,'(A)')
     *     ' +=============================================='
     *    //'===============================+'
C
      WRITE(LUPRI,'(//10X,A,8X,F10.4)')
     *     'Norm of Printed Amplitude Vector : ',SQRT(SUMOFP)
      WRITE(LUPRI,'(/10X,A43,F9.6)')
     *     'Printed all single excitations greater than',THR1
      IF (.NOT. (CCSEFF.OR.CCP2)) THEN
         WRITE(LUPRI,'(/10X,A43,F9.6)')
     *     'Printed all double excitations greater than',THR2
      ENDIF
C
 9990 FORMAT(1X,'| ',I1,3X,I1,2X,' | ',I3,5X,I3,4X,' | ',I8,9x,
     *       ' | ',12x,' | ',1x, F10.6,'  |')
 9991 FORMAT(1X,'| ',I1,1X,I1,1X,I1,1X,I1,' | ',
     *       I3,1X,I3,1X,I3,1X,I3,' | ',
     *       I8,1x,I8,' | ',I12,' | ',1x,F10.6,'  |')
C
      RETURN
      END
C  /* Deck cc_prtm */
      SUBROUTINE CC_PRTM(TM,ISYMD,ISYMTR)
C
#include <implicit.h>
#include <iratdef.h>
      DIMENSION TM(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMDJ = MULD2H(ISYMJ,ISYMD)
         ISYMCI = MULD2H(ISYMTR,ISYMDJ)
C
         NTOTCI = MAX(NT1AM(ISYMCI),1)
C
         WRITE(LUPRI,*) 'Transformede double exc. ampitudes Tm '
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'ISYMCI,ISYMJ = ',ISYMCI,ISYMJ
         DO 110 J = 1,NRHF(ISYMJ)
C
            KOFF3 = IT2BCD(ISYMCI,ISYMJ)
     *            + NT1AM(ISYMCI)*(J-1) + 1
C
            WRITE(LUPRI,*) 'FOR J= ',J
            CALL OUTPUT(TM(KOFF3),1,NT1AM(ISYMCI),1,1,
     *                  NT1AM(ISYMCI),1,1,LUPRI)
C
  110    CONTINUE
  100 CONTINUE
C
      END
C  /* Deck cc_prlam */
      SUBROUTINE CC_PRLAM(XLAMDP,XLAMDH,ISYML)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Ove Christiansen 13-7-1995 based on lammat by Henrik Koch.
C
C     PURPOSE:
C             Printout lambda matrices
C             for usual lambda in CC opt ISYML = 1
C             for c1 transformed ISYML = ISYMTR
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
      DIMENSION XLAMDH(*),XLAMDP(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
C
      WRITE(LUPRI,*)
      WRITE(LUPRI,*) 'In cc_prlam, symmetry of lambda matrices is :'
     *            ,ISYML
      CALL AROUND('Lambda Particle matrix ')
      DO 200 ISYMAO = 1,NSYM
         ISYMMO = MULD2H(ISYML,ISYMAO)
         WRITE(LUPRI,1) ISYMAO,ISYMMO
         WRITE(LUPRI,2)
         WRITE(LUPRI,3)
         IF (NRHF(ISYMMO) .EQ. 0) THEN
            WRITE(LUPRI,4)
            GOTO 210
         ENDIF
         KOFF1 = IGLMRH(ISYMAO,ISYMMO) + 1
         CALL OUTPUT(XLAMDP(KOFF1),1,NBAS(ISYMAO),1,NRHF(ISYMMO),
     *               NBAS(ISYMAO),NRHF(ISYMMO),1,LUPRI)
  210    WRITE(LUPRI,5)
         WRITE(LUPRI,6)
         IF (NVIR(ISYMMO) .EQ. 0) THEN
            WRITE(LUPRI,4)
            GOTO 220
         ENDIF
         KOFF2 = IGLMVI(ISYMAO,ISYMMO) + 1
         CALL OUTPUT(XLAMDP(KOFF2),1,NBAS(ISYMAO),1,NVIR(ISYMMO),
     *               NBAS(ISYMAO),NVIR(ISYMMO),1,LUPRI)
C
  220    CONTINUE
  200 CONTINUE
C
      CALL AROUND('Lambda Hole matrix ')
      DO 300 ISYMAO = 1,NSYM
         ISYMMO = MULD2H(ISYML,ISYMAO)
         WRITE(LUPRI,1) ISYMAO,ISYMMO
         WRITE(LUPRI,7)
         WRITE(LUPRI,8)
         IF (NRHF(ISYMMO) .EQ. 0) THEN
            WRITE(LUPRI,4)
            GOTO 310
         ENDIF
         KOFF1 = IGLMRH(ISYMAO,ISYMMO) + 1
         CALL OUTPUT(XLAMDH(KOFF1),1,NBAS(ISYMAO),1,NRHF(ISYMMO),
     *               NBAS(ISYMAO),NRHF(ISYMMO),1,LUPRI)
  310    WRITE(LUPRI,9)
         WRITE(LUPRI,10)
         IF (NVIR(ISYMMO) .EQ. 0) THEN
            WRITE(LUPRI,4)
            GOTO 320
         ENDIF
         KOFF2 = IGLMVI(ISYMAO,ISYMMO) + 1
         CALL OUTPUT(XLAMDH(KOFF2),1,NBAS(ISYMAO),1,NVIR(ISYMMO),
     *               NBAS(ISYMAO),NVIR(ISYMMO),1,LUPRI)
C
  320    CONTINUE
  300 CONTINUE
C
      RETURN
C
    1 FORMAT(/,/,7X,'Symmetry numbers ao,mo:',I5,I5)
    2 FORMAT(/,/,7X,'Lambda particle occupied part')
    3 FORMAT(7X,'-----------------------------')
    4 FORMAT(/,/,7X,'This symmetry is empty')
    5 FORMAT(/,/,7X,'Lambda particle virtual part')
    6 FORMAT(7X,'----------------------------')
    7 FORMAT(/,/,7X,'Lambda hole occupied part')
    8 FORMAT(7X,'-------------------------')
    9 FORMAT(/,/,7X,'Lambda hole virtual part')
   10 FORMAT(7X,'------------------------')
C
      END
C  /* Deck cc_print */
      SUBROUTINE CC_PRINT(XINT,DSRHF,ISYDIS)
C
C     Written by Ove Christiansen 24-7-1995
C
C     Purpose: Write out integrals.
C
#include <implicit.h>
C
      DIMENSION XINT(*),DSRHF(*)
C
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
C
      KOFF1 = 1
      KOFF2 = 1
      DO 100 ISYMJ = 1,NSYM
C
         ISYMG  = ISYMJ
         ISYMAB = MULD2H(ISYMG,ISYDIS)
C
         NNBSAB = MAX(NNBST(ISYMAB),1)
         NBASG  = MAX(NBAS(ISYMG),1)
C
         CALL AROUND( ' AO - INTEGRALS ')
         WRITE(LUPRI,*) 'ISYMG: ',ISYMG
         CALL OUTPUT(XINT(KOFF1),1,NNBST(ISYMAB),1,NBAS(ISYMG),
     *               NNBSAB,NBASG,1,LUPRI)
C
         CALL AROUND( 'DSRHF ')
         WRITE(LUPRI,*) 'ISYMJ: ',ISYMJ
         CALL OUTPUT(DSRHF(KOFF2),1,NNBST(ISYMAB),1,NRHF(ISYMJ),
     *               NNBSAB,NBASG,1,LUPRI)
C
         KOFF1 = KOFF1 + NNBST(ISYMAB)*NBAS(ISYMG)
         KOFF2 = KOFF2 + NNBST(ISYMAB)*NRHF(ISYMJ)
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_prei */
      SUBROUTINE CC_PREI(EI1,EI2,ISYMEI,LEI1MO)
C
#include <implicit.h>
#include <iratdef.h>
      DIMENSION EI1(*),EI2(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
C
      CALL AROUND( 'EI1 -intermediate ')
C
      DO 100 ISYMB = 1,NSYM
C
         ISYMA = MULD2H(ISYMB,ISYMEI)
C
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'ISYMA,ISYMB = ',ISYMA,ISYMB
         WRITE(LUPRI,*) ' '
C
         IF (LEI1MO.EQ.1) THEN
            KOFF1 = IMATAB(ISYMA,ISYMB) + 1
            LB    = NVIR(ISYMB)
         ELSE
            KOFF1 = IEMAT1(ISYMA,ISYMB) + 1
            LB    = NBAS(ISYMB)
         ENDIF
         CALL OUTPUT(EI1(KOFF1),1,NVIR(ISYMA),1,LB,
     *               NVIR(ISYMA),LB,1,LUPRI)
C
  100 CONTINUE
C
      CALL AROUND( 'EI2 -intermediate ')
C
      DO 200 ISYMJ = 1,NSYM
C
         ISYMI = MULD2H(ISYMJ,ISYMEI)
C
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'ISYMI,ISYMJ = ',ISYMI,ISYMJ
         WRITE(LUPRI,*) ' '
C
         KOFF1 = 1 + IMATIJ(ISYMI,ISYMJ)
         CALL OUTPUT(EI2(KOFF1),1,NRHF(ISYMI),1,NRHF(ISYMJ),
     *               NRHF(ISYMI),NRHF(ISYMJ),1,LUPRI)
C
  200 CONTINUE
C
      END
C  /* Deck cc_praoden */
      SUBROUTINE CC_PRAODEN(DENS,ISYDEN)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Ove Christiansen 13-7-1995
C
C     Purpose: Print density matrices
C	       calculated from ordinary lambda matrices
C              and from C1 transformed lambda hole matrix.
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
      DIMENSION DENS(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
C
      CALL AROUND('Lamda density matrix')
      KOFF1 = 1
      DO 110 ISYMB = 1,NSYM
C
         ISYMA = MULD2H(ISYMB,ISYDEN)
         WRITE(LUPRI,*) 'Symmetry number alfa,beta: ',ISYMA,ISYMB
         NBASA = NBAS(ISYMA)
         NBASB = NBAS(ISYMB)
         CALL OUTPUT(DENS(KOFF1),1,NBASA,1,NBASB,NBASA,NBASB,1,LUPRI)
         KOFF1 = KOFF1 + NBASA*NBASB
C
  110 CONTINUE
C
      END
C  /* Deck cc_prfckao */
      SUBROUTINE CC_PRFCKAO(FOCK,ISYFCK)
C
C     Ove Christiansen 14-7-1994
C
C     Purpose: Print AO Fock matrix.
C
#include <implicit.h>
      DIMENSION FOCK(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
C
      CALL AROUND('The AO Fock matrix')
      KOFF1 = 1
      DO 50 ISYMB = 1,NSYM
         ISYMA = MULD2H(ISYMB,ISYFCK)
         WRITE(LUPRI,*) 'Symmetry of A,B : ',ISYMA,ISYMB
         NBASA = NBAS(ISYMA)
         NBASB = NBAS(ISYMB)
         CALL OUTPUT(FOCK(KOFF1),1,NBASA,1,NBASB,NBASA,NBASB,1,LUPRI)
         KOFF1 = KOFF1 + NBAS(ISYMA)*NBAS(ISYMB)
   50 CONTINUE
C
      END
C  /* Deck cc_prfckmo */
      SUBROUTINE CC_PRFCKMO(FOCK,ISYFCK)
C
C     Ove Christiansen 14-7-1994
C
C     Purpose: Print MO Fock matrix.
C
#include <implicit.h>
      DIMENSION FOCK(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
C
      CALL AROUND('The MO Fock matrix ')
C
      KOFF1 = 1
      KOFF2 = NLRHFR(ISYFCK) + 1
      DO 300 ISYMQ = 1,NSYM
C
         ISYMP = MULD2H(ISYMQ,ISYFCK)
C
         WRITE(LUPRI,1) ISYMP,ISYMQ
         WRITE(LUPRI,2)
         WRITE(LUPRI,3)
         IF (NRHF(ISYMQ) .EQ. 0) THEN
            WRITE(LUPRI,4)
            GOTO 310
         ENDIF
         CALL OUTPUT(FOCK(KOFF1),1,NORB(ISYMP),1,NRHF(ISYMQ),
     *               NORB(ISYMP),NRHF(ISYMQ),1,LUPRI)
  310    CONTINUE
C
         WRITE(LUPRI,5)
         WRITE(LUPRI,6)
         IF (NVIR(ISYMQ) .EQ. 0) THEN
            WRITE(LUPRI,4)
            GOTO 320
         ENDIF
         CALL OUTPUT(FOCK(KOFF2),1,NORB(ISYMP),1,NVIR(ISYMQ),
     *               NORB(ISYMP),NVIR(ISYMQ),1,LUPRI)
C
  320    CONTINUE
C
         KOFF1 = KOFF1 + NORB(ISYMP)*NRHF(ISYMQ)
         KOFF2 = KOFF2 + NORB(ISYMP)*NVIR(ISYMQ)
C
  300 CONTINUE
C
    1 FORMAT(/,/,7X,'Symmetry of P,Q :',I5,I5)
    2 FORMAT(/,/,7X,'Occupied part')
    3 FORMAT(7X,'-------------')
    4 FORMAT(/,/,7X,'This symmetry is empty')
    5 FORMAT(/,/,7X,'Virtual part')
    6 FORMAT(7X,'------------')
C
      END
C  /* Deck cc_rvec */
      SUBROUTINE CC_RVEC(LU,FIL,LLEN,LEN,NR,VEC)
C
C     Ove Christiansen 22-1-1996:
C
C             Read vector NR on file FIL with unit number LU and
C             put into VEC. The position is calculated relative to
C             LLEN which is the length of the vectors according to
C             which the vectors was put there in the first place.
C
C     tbp: Use fortran I/O for Cholesky....
C
C
#include <implicit.h>
Cholesky
#include <maxorb.h>
#include <ccdeco.h>
#include <priunit.h>
Cholesky
      DIMENSION VEC(*)
      CHARACTER FIL*(*)
C
      IF (CHOINT) THEN

C        Fortran I/O section.
C        --------------------

         IF (LLEN .NE. LEN) THEN
            WRITE(LUPRI,'(//,5X,A,/,5X,A,I10,/,5X,A,I10,/)')
     &     'Fatal error in CC_RVEC: LLEN .ne. LEN:',
     &     'LLEN = ',LLEN,
     &     'LEN  = ',LEN
           CALL QUIT('Fatal [FORTRAN] I/O error in CC_RVEC')
         ENDIF

         CALL CHO_MOREAD(VEC,LEN,1,NR,LU)

      ELSE

C        crayio section.
C        ---------------
C
         IOFF = 1 + LLEN*(NR-1)
C
         IF (LEN .GT. 0) THEN
            CALL GETWA2(LU,FIL,VEC,IOFF,LEN)
C
         ENDIF
C
      END IF
C
      RETURN
      END
C  /* Deck cc_wvec */
      SUBROUTINE CC_WVEC(LU,FIL,LLEN,LEN,NR,VEC)
C
C     Ove Christiansen 22-1-1996:
C
C             Write vector NR given in VEC on file FIL with unit number LU.
C             The position is calculated relative to
C             LLEN which is the length of the vectors according to
C             which the vectors was put there in the first place.
C
#include <implicit.h>
Cholesky
#include <maxorb.h>
#include <ccdeco.h>
#include <priunit.h>
Cholesky
      DIMENSION VEC(*)
      CHARACTER FIL*(*)
C
      IF (CHOINT) THEN

C        Fortran I/O section.
C        --------------------

         IF (LLEN .NE. LEN) THEN
            WRITE(LUPRI,'(//,5X,A,/,5X,A,I10,/,5X,A,I10,/)')
     &     'Fatal error in CC_WVEC: LLEN .ne. LEN:',
     &     'LLEN = ',LLEN,
     &     'LEN  = ',LEN
           CALL QUIT('Fatal [FORTRAN] I/O error in CC_WVEC')
         ENDIF

         CALL CHO_MOWRITE(VEC,LEN,1,NR,LU)

      ELSE

C        crayio section.
C        ---------------

         IOFF = 1 + LLEN*(NR-1)
C
         IF (LEN .GT. 0) THEN
            CALL PUTWA2(LU,FIL,VEC,IOFF,LEN)
         ENDIF
C
      ENDIF
C
      RETURN
      END
c*DECK CC_JACEXP
      SUBROUTINE CC_JACEXP(WORK,LWORK)
C
C----------------------------------------------------------------------
C     Calculates  the CC Jacobian by Transformation of
C     unit vectors and writes it to disk in CCLR_JAC.
C
C     Written by Ove Christiansen 13-10-1995
C
C     Christof Haettig, October 1998:
C       changed to a dummy routine, because call of CCLR_RHTR and
C       CCLHST via CCLR_TRR was switched off since quite some time...
C       after change of F matrix implementation CCLR_TRR was 
C       completly dummy and removed...
C
C----------------------------------------------------------------------
C
#include <implicit.h>
#include <dummy.h>
#include <maxorb.h>
#include <iratdef.h>
#include <ccorb.h>
#include <aovec.h>
#include <ccsdinp.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <leinf.h>
#include <priunit.h>
C
      DIMENSION WORK(LWORK)
      PARAMETER (XHALF = 0.5D00,XMTWO = -2.0D00, DELTA = 1.0D-08)
C
CCH   IF (IPRINT .GT. 10) THEN
         CALL AROUND(' START OF CC_JACEXP')
         CALL AROUND(' The routine is outdated... Leaving CC_JACEXP.')
         RETURN
CCH   ENDIF
#ifdef OUTDATED_ROUTINE
C
C-------------------
C     Open jac file.
C-------------------
C
      LUJAC = -1
      CALL GPOPEN(LUJAC,'CCLR_JAC','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     *            .FALSE.)
      REWIND(LUJAC)
C
C----------------------------
C     Work space allocations.
C----------------------------
C
      NTAMP      = NT1AM(ISYMTR) + NT2AM(ISYMTR)
C
      KV         = 1
      KEND1      = 1 + NTAMP
      LWRK1      = LWORK - KEND1
      IF (LWRK1 .LT. 0 )
     &     CALL QUIT('INSUFFICIENT WORK SPACE IN CC_JACEXP')
C
C---------------------------------------------------------------------
C     Loop over vectors and save on disk the coloumns of the jacobian.
C---------------------------------------------------------------------
C
      CALL DZERO(WORK(KV),NTAMP)
C
      JACEXP = .FALSE.
C
      DO 100 I = 1,NTAMP
C
         WORK(KV+I-1) = 1.0D00
CCH      CALL CCLR_TRR(1,0,WORK(KV),DUMMY,
CCH  *                 XLINLE,WORK(KEND1),LWRK1)
         WRITE(LUJAC) (WORK(KEND1+J-1),J=1,NTAMP)
         WORK(KV+I-1) = 0.0D00
C
 100  CONTINUE
C
      JACEXP = .TRUE.
C
      IF (IPRINT .GT. 10) THEN
         CALL AROUND(' END OF CC_JACEXP ')
      ENDIF
C
      CALL GPCLOSE(LUJAC,'KEEP')
C
      RETURN
#endif
      END
C
C  /* Deck cc_pronelao */
      SUBROUTINE CC_PRONELAO(FOCK,ISYFCK)
C
C     Purpose: Print AO one-electron matrix.
C
#include <implicit.h>
      DIMENSION FOCK(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
C
      KOFF1 = 1
      DO ISYMB = 1,NSYM
         ISYMA = MULD2H(ISYMB,ISYFCK)
         WRITE(LUPRI,*) 'Symmetry of A,B : ',ISYMA,ISYMB
         NBASA = NBAS(ISYMA)
         NBASB = NBAS(ISYMB)
         CALL OUTPUT(FOCK(KOFF1),1,NBASA,1,NBASB,NBASA,NBASB,1,LUPRI)
         KOFF1 = KOFF1 + NBAS(ISYMA)*NBAS(ISYMB)
      END DO
C
      END
C
C  /* Deck cc_prpr12 */
      SUBROUTINE CC_PRPR12(TR12AM,ISYM,NTR12,LHEAD)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Christian Neiss, Feb. 2005
C     PRint Packed R12-vector - PRPR12 (in general non. tot. sym.)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
C
      LOGICAL   LHEAD
      DIMENSION TR12AM(*)

C
C------------------------------------
C     Write R12 vector.
C------------------------------------
C
      DO 20 I = 1,NTR12
C
         IF (LHEAD) CALL AROUND('R12 part of vector ')
         IF (NTR12.GT.1) WRITE(LUPRI,*) 'Vector Number:',I
C
         DO 200 ISYMLJ = 1,NSYM
C
            ISYMKI = MULD2H(ISYMLJ,ISYM)
C
            KOFF = ITR12AM(ISYMKI,ISYMLJ) + 1

            IF (ISYMKI.EQ.ISYMLJ) THEN
C
               WRITE(LUPRI,*) 'Symmetry block number(ki,lj): ',
     &              ISYMKI,ISYMLJ
               IF (NMATKI(ISYMKI).EQ.0) THEN
                  WRITE(LUPRI,*) 'This symmetry is empty'
               ELSE
                  CALL OUTPAK(TR12AM(KOFF),NMATKI(ISYMKI),1,LUPRI)
               END IF
               WRITE(LUPRI,*) ' '
C
            ELSE IF (ISYMLJ.GT.ISYMKI) THEN
C
               WRITE(LUPRI,*) 'Symmetry block number(ki,lj): ',
     &              ISYMKI,ISYMLJ
               IF (NMATKI(ISYMKI).EQ.0 .OR. NMATKI(ISYMLJ).EQ.0) THEN
                  WRITE(LUPRI,*) 'This symmetry is empty'
               ELSE
                  CALL OUTPUT(TR12AM(KOFF),1,NMATKI(ISYMKI),
     &                        1,NMATKI(ISYMLJ),
     *                        NMATKI(ISYMKI),NMATKI(ISYMLJ),1,LUPRI)
               END IF
               WRITE(LUPRI,*) ' '
C
            ENDIF
C
  200    CONTINUE
C
   20 CONTINUE
C
      RETURN
      END
C
C  /* Deck cc_prsqr12 */
      SUBROUTINE CC_PRSQR12(TR12AM,ISYM,TRANS,NTR12,LHEAD)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Christian Neiss  Feb. 2005
C     PRint SQuared R12 vector. general non. tot. sym.
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
C
      LOGICAL     LHEAD
      CHARACTER*1 TRANS
      DIMENSION   TR12AM(*)
C
      IF ((TRANS.NE.'T').AND.(TRANS.NE.'N')) 
     &   CALL QUIT('Illegal value for "TRANS" in CC_PRSQR12')
C
C
C------------------------------------
C     Write R12 vector.
C------------------------------------
C
      DO 20 I = 1, NTR12
C
         IF (LHEAD) CALL AROUND('R12 part of vector ')
        
         IF (NTR12.GT.1) WRITE(LUPRI,*) 'Vector Number:',I
C
         DO 200 ISYMIJ = 1,NSYM
C
            ISYMKL = MULD2H(ISYMIJ,ISYM)
C
            IF (TRANS.EQ.'T') THEN
             WRITE(LUPRI,*) 'Symmetry block (ij,kl):',
     &             ISYMIJ,ISYMKL
            ELSE IF (TRANS.EQ.'N') THEN
             WRITE(LUPRI,*) 'Symmetry block (kl,ij):',
     &             ISYMKL,ISYMIJ
            END IF
            IF (NMATKL(ISYMKL).EQ.0 .OR. NMATIJ(ISYMIJ).EQ.0) THEN
            WRITE(LUPRI,*) 'This symmetry is empty'
            ELSE
              IF (TRANS.EQ.'T') THEN
               KOFF = ITR12SQT(ISYMIJ,ISYMKL) + 1
               CALL OUTPUT(TR12AM(KOFF),1,NMATIJ(ISYMIJ),1,
     *                   NMATKL(ISYMKL),NMATIJ(ISYMIJ),NMATKL(ISYMKL),
     *                   1,LUPRI)
              ELSE IF (TRANS.EQ.'N') THEN
               KOFF = ITR12SQ(ISYMKL,ISYMIJ) + 1
               CALL OUTPUT(TR12AM(KOFF),1,NMATKL(ISYMKL),1,
     *                   NMATIJ(ISYMIJ),NMATKL(ISYMKL),NMATIJ(ISYMIJ),
     *                   1,LUPRI)
              END IF
            END IF
C
            WRITE(LUPRI,*) ' '
C
  200    CONTINUE
C
   20 CONTINUE
C
      RETURN
      END
C
C  /* Deck cclr_transsq */
      SUBROUTINE CCLR_TRANSSQ(MATRIX,NROWS)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Christian Neiss  Mar. 2005
C     Transpose a square matrix 'in place'
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      implicit none
#include <priunit.h>
C
      INTEGER IDXIJ,IDXJI,NROWS,I,J
#if defined (SYS_CRAY)
      REAL X1,X2,matrix(*)
#else
      DOUBLE PRECISION X1,X2,matrix(*) 
#endif
C
      do i = 1, nrows
        do j = 1, i
          idxij = nrows*(j-1) + i
          idxji = nrows*(i-1) + j
          x1 = matrix(idxij)
          x2 = matrix(idxji)
          matrix(idxij) = x2
          matrix(idxji) = x1
        end do
      end do

      RETURN
      END
C
C  /* Deck cclr_trsqr12 */
      SUBROUTINE CCLR_TRSQR12(MATRIX,ISYM)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Christian Neiss  April 2005
C     Transpose a symmetry packed square matrix X_kl^ij:
C     X(kl,ij) --> X(ij,kl) 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      implicit none
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER IDXKLIJ,IDXIJKL,IDXKL,IDXIJ,ISYM,ISYMKL,ISYMIJ
#if defined (SYS_CRAY)
      REAL X1,X2,matrix(*)
#else
      DOUBLE PRECISION X1,X2,matrix(*)
#endif
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
C
      if (locdbg) then
        do I = 1, NSYM
        write(lupri,*) 'ITR12SQ = ', (ITR12SQ(I,J), J=1,NSYM)
        write(lupri,*) 'ITR12SQT= ', (ITR12SQT(I,J),J=1,NSYM)
        end do
        call around('Input Matrix in CCLR_TRSQR12')
        call cc_prsqr12(matrix,isym,'N',1,.false.)
      end if
C
      do isymkl = 1, nsym
        isymij = muld2h(isym,isymkl)
        do idxkl = 1, nmatkl(isymkl)  
          do idxij = 1, nmatij(isymij)
            idxklij = itr12sq(isymkl,isymij) + 
     &                nmatkl(isymkl)*(idxij-1) + idxkl
            idxijkl = itr12sqt(isymij,isymkl) + 
     &                nmatij(isymij)*(idxkl-1) + idxij
            if (idxklij.lt.idxijkl) then
              x1 = matrix(idxklij)
              x2 = matrix(idxijkl)
              matrix(idxklij) = x2
              matrix(idxijkl) = x1
            end if
          end do
        end do
      end do
C
      if (locdbg) then
        call around('Transposed Matrix in CCLR_TRSQR12')
        call cc_prsqr12(matrix,isym,'T',1,.false.)
      end if
C
      RETURN
      END
C
