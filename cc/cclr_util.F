C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
C  /* Deck cclr_dia */
      SUBROUTINE CCLR_DIA(DIAA,DIAB,TRIPLET,WORK,LWORK)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Written by Ove Christiansen 4-9-1994
C     Modified 22-11 to use Ematrices from
C     calculation of the CCSD vector.
C     NB no symmetry in E part 26-1-1995
C     For CCS and CC2 it uses Fock diagonal.
C     Triplet extensions by Christof Haettig & Kasper Hald, May 1999
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
#include <dummy.h>
      PARAMETER (ONE = 1.0D0, IZERO = 0, HALF = 0.5D00, TMPKH= 1.0D08)
#include <priunit.h>
#include <inftap.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <cclr.h>
#include <iratdef.h>
#include <leinf.h>
      DIMENSION WORK(LWORK),DIAA(*),DIAB(*)
      LOGICAL FOCKDI,EMATDI,TRIPLET
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) -3 )/2 + I + J
C
C-------------------------------------------
C     Testoptions.
C     Use fock diagonal for the moment.
C-------------------------------------------
C
      CALL QENTER('CCLR_DIA')
      EMATDI = .FALSE.
      FOCKDI = .TRUE.
      IF (CCS.OR.CCP2.OR.CC2 ) THEN
         FOCKDI = .TRUE.
         EMATDI = .FALSE.
      ENDIF
C
      IF (FOCKDI) THEN
C
C--------------------------------
C        Allocation of workspace.
C--------------------------------
C
         KFOCKD  = 1
         KEND1   = KFOCKD  + NORBTS
         LWRK1   = LWORK   - KEND1
C
         IF ( LWRK1 .LT. 0 ) THEN
            CALL QUIT('Insufficient space in CCLR_DIA ')
         ENDIF
C
C----------------------------------------
C        Read canonical orbital energies.
C----------------------------------------
c
         CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND LUSIFC
C
         CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
         READ (LUSIFC) POTNUC,EMY,EACTIV,EMCSCF,ISTATE,ISPIN,NACTEL,
     *                 LSYM
C
         CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
         READ (LUSIFC)
         READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBTS)
C
         CALL GPCLOSE(LUSIFC,'KEEP')
C
         IF (IPRINT.GT.100) THEN
            CALL AROUND('FOCK DIAGONAL ')
            CALL OUTPUT(WORK(KFOCKD),1,NORBTS,1,1,NORBTS,1,1,LUPRI)
         ENDIF
C
         IF (FROIMP .OR. FROEXP)
     *      CALL CCSD_DELFRO(WORK(KFOCKD),WORK(KEND1),LWRK1)
C
C---------------------------------------------------
C        FOCK based  approximative  A(1,1) diagonal.
C---------------------------------------------------
C
         DO 100 ISYMA = 1,NSYM
C
            ISYMI  = MULD2H(ISYMA,ISYMTR)
C
            DO 110 I = 1,NRHF(ISYMI)
C
               MI = IORB(ISYMI) + I
C
               DO 120 A = 1,NVIR(ISYMA)
C
                   MA = IORB(ISYMA) + NRHF(ISYMA) +  A
C
                   KAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
C
                   DIAA(KAI) = WORK(MA) - WORK(MI)
C
  120          CONTINUE
C
  110       CONTINUE
C
  100    CONTINUE
C
         IF (IPRINT.GT. 20 .OR. DEBUG) THEN
            CALL AROUND('FOCK BASED DIAGONAL single exc. ')
            CALL OUTPUT(DIAA,1,NT1AM(ISYMTR),1,1,NT1AM(ISYMTR),1,1,
     &           LUPRI)
         ENDIF
C
C-----------------------------------
C     Approximative A(2,2) diagonal.
C     Skip This for CCS!
C-----------------------------------
C
         IF (CCS.OR.CCP2) THEN
            DO 190 I =1,LETYPA-NT1AM(ISYMTR)
               DIAA(NT1AM(ISYMTR)+I) = 1.0D08
  190       CONTINUE
            CALL QEXIT('CCLR_DIA')
            RETURN
         ENDIF
C
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

CCH
CCH                           old stuff... seems to be bugy...
CCH
C                             IF ((ISYMAI.EQ.ISYMBJ).AND.
C    *                            (NAI .LT. NBJ)) GOTO 260
C
C                             IF (ISYMAI.EQ.ISYMBJ) THEN
C                                NAIBJ = IT2AM(ISYMAI,ISYMBJ)
C    *                                + INDEX(NAI,NBJ)
C                             ELSE
C                                NAIBJ = IT2AM(ISYMAI,ISYMBJ)
C    *                               + NT1AM(ISYMAI)*(NBJ-1) + NAI
C                             ENDIF

CCH
CCH                           new version... hopefully correct...
CCH
                              IF (ISYMAI.EQ.ISYMBJ) THEN
                                 NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                               + INDEX(NAI,NBJ)
                              ELSE IF (ISYMAI.LT.ISYMBJ) THEN
                                 NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                               + NT1AM(ISYMAI)*(NBJ-1) + NAI
                              ELSE IF (ISYMAI.GT.ISYMBJ) THEN
                                 NAIBJ = IT2AM(ISYMBJ,ISYMAI)
     *                               + NT1AM(ISYMBJ)*(NAI-1) + NBJ
                              ENDIF
C
                              DIAA(NAIBJ+NT1AM(ISYMTR)) =
     *                          WORK(MA)+WORK(MB)-WORK(MI)-WORK(MJ)
C
  260                      CONTINUE
  250                   CONTINUE
  240                CONTINUE
  230             CONTINUE
  220          CONTINUE
  210       CONTINUE
  200    CONTINUE
C
         IF (TRIPLET) THEN
            KOFF1  = 1
            KOFF2P = KOFF1  + NT1AM(ISYMTR)
            KOFF2M = KOFF2P + NT2AM(ISYMTR)

C           --------------------------------------------------------
C           2- block as for singlet vector, but with the ai,ai 
C           diagonal set to a huge value to ensure that the
C           corresponding elements in the trial vectors are zero.
C           --------------------------------------------------------
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

C           -----------------------------------------------------------
C           2+ block as for singlet vector, but times a factor of 0.5
C           and with the elements a<=b or i<=j set to huge value
C           to ensure that the corresponding elements in the trial
C           vectors are zero 
C           -----------------------------------------------------------
!
            CALL DCOPY(NT2AM(ISYMTR),TMPKH,0,DIAA(KOFF2P),1)
!
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
     *                           (MA .LT. MB)) THEN
C
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                               + INDEX(NAI,NBJ)
                              DIAA(NT1AM(ISYMTR)+NAIBJ) =
     *                             DIAA(KOFF2M-1+NAIBJ)
!
                             ENDIF
                           ELSE IF (ISYMAI.LT.ISYMBJ) THEN
!
                             IF (((MA .LT. MB) .AND.
     *                            (MI .LT. MJ)) .OR.
     *                           ((MA .GT. MB) .AND.
     *                            (MI .GT. MJ))) THEN
!
                                 NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                               + NT1AM(ISYMAI)*(NBJ-1) + NAI
                                 DIAA(NT1AM(ISYMTR)+NAIBJ) =
     *                                DIAA(KOFF2M-1+NAIBJ)
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
         END IF
C
      ENDIF
C
      IF (EMATDI) THEN
C
C--------------------------------
C        Allocation of workspace.
C--------------------------------
C
         KEI1    = 1
         KEI2   = KEI1 + NMATAB(ISYMOP)
         KEND1  = KEI2 + NMATIJ(ISYMOP)
         LWRK1  = LWORK - KEND1
C
         IF ( LWRK1 .LT. 0 ) THEN
            CALL QUIT('Insufficient spaces in CCLR_DIA ')
         ENDIF
C
C-------------------------
C        Readin Ematrices.
C-------------------------
C
         LUE1 = 0
         CALL GPOPEN(LUE1,'CC_E1IM','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LUE1)
         READ (LUE1)(WORK(KEI1+ J-1),J = 1,NMATAB(ISYMOP))
         CALL GPCLOSE(LUE1,'KEEP')
C
         LUE2 = 0
         CALL GPOPEN(LUE2,'CC_E2IM','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LUE2)
         READ (LUE2)(WORK(KEI2+ J-1),J = 1,NMATIJ(ISYMOP))
         CALL GPCLOSE(LUE2,'KEEP' )
C
C---------------------------------------------------
C        EMAT based  approximative  A(1,1) diagonal.
C---------------------------------------------------
C
C        old crap deleted - write new if you want E diagonal.
C
C-----------------------------------
C     Approximative A(2,2) diagonal.
C-----------------------------------
C
      ENDIF
C
      IF (IPRINT.GT. 20 .OR. DEBUG) THEN
         CALL AROUND('APPROXIMATIVE  DIAGONAL 22 BLOCK ')
         IF (TRIPLET) THEN
           CALL OUTPUT(DIAA(1+NT1AM(ISYMTR)),1,2*NT2AM(ISYMTR),
     *                 1,1,2*NT2AM(ISYMTR),1,1,LUPRI)
         ELSE
           CALL OUTPUT(DIAA(1+NT1AM(ISYMTR)),1,1*NT2AM(ISYMTR),
     *                 1,1,1*NT2AM(ISYMTR),1,1,LUPRI)
         END IF
      ENDIF
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
      SUBROUTINE CCLR_JAC(WORK,LWORK)
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
      NTEMP  = NT1AM(ISYMTR) + NT2AM(ISYMTR)
      NTEMP2 = NTEMP*(NC1VEC + NC2VEC)
C
C------------------------------------------------------------
C     Calculate Jacobian by Transformation with unit vectors.
C     First elements in workspace contains the jacobian.
C------------------------------------------------------------
C
      IF (IPRINT.GT.10) WRITE(LUPRI,*) 'CCLR_JAC Workspace:',LWORK
C     
      KEND1 = 1
      LWRK1 = LWORK
      IF (JACTST) THEN
         KJACAN  = 1
         KEND1   = 1 + NTEMP2
         LWRK1   = LWORK - KEND1
         IF ( LWRK1 .LE. 0 )
     &        CALL QUIT('Too little workspace in CCLR_JAC')
         CALL AROUND( 'Calculation of analytical CC Jacobian')
         CALL CCLR_JACAN(NC1VEC,NC2VEC,WORK(KJACAN),LWORK)
      ENDIF
C
C-------------------------------------------------------
C     Calculate Jacobian by finite difference.
C     both for jactst and jacexp.
C     First elements in workspace contains the jacobian.
C-------------------------------------------------------
C
      write(LUPRI,*) 'lwork,ntemp2,ken1',lwork,ntemp2,kend1
      KFDJAC  = KEND1
      KEND2   = KFDJAC + NTEMP2
      LWRK2   = LWORK - KEND2
      write(LUPRI,*) 'lwrk2,ken2',lwrk2,kend2
      IF ( LWRK2 .LE. 0 )
     &     CALL QUIT('Too little workspace in CCLR_JAC-2')
      CALL AROUND( 'Calculation of finite difference  CC Jacobian')
      RSPIM2 =  RSPIM
      RSPIM  = .FALSE.
      CALL CCLR_FDJAC(NC1VEC,NC2VEC,WORK(KFDJAC),LWRK1)
      RSPIM  = RSPIM2
C
C--------------------------------------------------------
C        calculate difference between the two jacobiants.
C--------------------------------------------------------
C
      IF (JACTST) THEN
         CALL DAXPY(NTEMP2,-1.0D00,WORK(KFDJAC),1,WORK(KJACAN),1)
         IF ( IPRINT .GT. 40) THEN
            CALL AROUND( 'DIFFERENCE. CC JACOBIANT - 11 & 21 PART  ' )
            CALL OUTPUT(WORK(KJACAN),1,NTEMP,1,NC1VEC,NTEMP,NC1VEC,1,
     &           LUPRI)
            CALL AROUND( 'DIFFERENCE. CC JACOBIANT - 12 & 22 PART  ' )
            CALL OUTPUT(WORK(KJACAN+NTEMP*NC1VEC),1,NTEMP,1,NC2VEC,
     *                  NTEMP,NC2VEC,1,LUPRI)
         ENDIF
         DIFNOR  = DDOT(NTEMP2,WORK(KJACAN),1,WORK(KJACAN),1)
         WRITE(LUPRI,*) '  '
         WRITE(LUPRI,*) ' The norm of the difference between fd and '//
     *               ' ana jacobiant is ',SQRT(DIFNOR)
         WRITE(LUPRI,*) '  '
         WRITE(LUPRI,*) ' END OF JACTST TEST'
         WRITE(LUPRI,*) '  '
      ENDIF
C
      IF (IPRINT .GT. 10) THEN
         CALL AROUND(' END OF CCLR_JAC ')
      ENDIF
C
      RETURN
      END
c*DECK CCLR_JACAN
      SUBROUTINE CCLR_JACAN(NC1VEC,NC2VEC,WORK,LWORK)
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
C
      IF (IPRINT .GT. 10) THEN
         CALL AROUND(' START OF CCLR_JACAN')
      ENDIF
C
C----------------------------
C     Work space allocations.
C----------------------------
C
      NTAMP      = NT1AMX + NT2AMX
      NTAMP2     = NTAMP*(NC1VEC + NC2VEC )
      KJACOBI    = 1
      KEND1      = 1 + NTAMP2
      LWRK1      = LWORK - KEND1
      IF (LWRK1 .LT. NT1AMX + NT2AMX) 
     *        CALL QUIT('INSUFFICIENT WORK SPACE IN CCLR_JACAN')
C
C--------------------------------------------------------------
C     Put up nc1vec C1 unit vectors and nc2vec C2 unit vectors.
C--------------------------------------------------------------
C
      KC1 = KJACOBI
      KC2 = KC1 + (NT1AMX + NT2AMX)*NC1VEC
      CALL CCLR_UNIT(NC1VEC,NC2VEC,WORK(KC1),WORK(KC2))
C
C---------------------------------------------------------
C     Make the transformation of the unit vectors.
C     Transformed vectors are returned as first
C     elements in work.
C     Work follows immediately after vectors as it should.
C---------------------------------------------------------
C
      DO 100 I = 1, NC1VEC + NC2VEC
         CALL DCOPY(NTAMP,WORK(KC1 + (NT1AMX + NT2AMX)*(I-1)),1,
     *              WORK(KEND1),1)
         ISIDE = 1
         CALL CC_ATRR(ISYMOP,ISIDE,WORK(KEND1),LWRK1)
         CALL DCOPY(NTAMP,WORK(KEND1),1,
     *              WORK(KC1 + (NT1AMX + NT2AMX)*(I-1)),1) 
  100 CONTINUE
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
      CALL CCLR_NORMS(X11,X21,X12,X22,WORK(KJACOBI),NC1VEC,NC2VEC)
C
      WRITE(LUPRI,*) 'NORM OF 11 PART OF UNIT VECT. JAC. ', SQRT(X11)
      WRITE(LUPRI,*) 'NORM OF 21 PART OF UNIT VECT. JAC. ', SQRT(X21)
      WRITE(LUPRI,*) 'NORM OF 12 PART OF UNIT VECT. JAC. ', SQRT(X12)
      WRITE(LUPRI,*) 'NORM OF 22 PART OF UNIT VECT. JAC. ', SQRT(X22)
C
C------------------------------------------------
C     Print the columbs of the jacobian obtained.
C------------------------------------------------
C
      IF (IPRINT .GT. 40) THEN
C
         CALL AROUND( 'CC JACOBIANT I PRESUME - A*C1 PART ' )
         CALL OUTPUT(WORK(KJACOBI),1,NTAMP,1,NC1VEC,NTAMP,NC1VEC,1,
     &        LUPRI)
         CALL AROUND( 'CC JACOBIANT I PRESUME - A*C2 PART ' )
         CALL OUTPUT(WORK(KJACOBI+NTAMP*NC1VEC),1,NTAMP,1,NC2VEC,
     *               NTAMP,NC2VEC,1,LUPRI)
C
      ENDIF
C
      IF (IPRINT .GT. 10) THEN
         CALL AROUND(' END OF CCLR_JACAN ')
      ENDIF
C
      RETURN
      END
C
c*DECK CCLR_UNIT
      SUBROUTINE CCLR_UNIT(NC1VEC,NC2VEC,C1AM,C2AM)
C
C     Put up the first NC1VEC and NC2VEC unit vectors in the single and double
C     excitation space respectively.
C
#include <implicit.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <cclr.h>
C
      DIMENSION C1AM(NT1AMX+ NT2AMX,NC1VEC),C2AM(NT1AMX + NT2AMX,NC2VEC)
C
C--------------------------------------
C     Put up the required unit vectors.
C--------------------------------------
C
      DO 100 J = 1, NC1VEC
C
            CALL DZERO(C1AM(1,J),NT1AMX + NT2AMX)
            C1AM(J,J) = 1.0D00
C
 100  CONTINUE
C
      DO 200 J = 1, NC2VEC
C
            CALL DZERO(C2AM(1,J),NT1AMX + NT2AMX)
            C2AM(J + NT1AMX,J) = 1.0D00
C
 200  CONTINUE
C
      RETURN
      END
c*DECK CCLR_NORMS
      SUBROUTINE CCLR_NORMS(X11,X21,X12,X22,XJAC,NC1VEC,NC2VEC)
#include <implicit.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <cclr.h>
C
      DIMENSION XJAC(NT1AMX + NT2AMX,NC1VEC + NC2VEC)
C
      NTAMP = NT1AMX + NT2AMX
      X11 = 0.0D00
      X12 = 0.0D00
      X21 = 0.0D00
      X22 = 0.0D00
C
      DO 100 I = 1, NC1VEC
C
         X11 = X11 + DDOT(NT1AMX,XJAC(1,I),1,XJAC(1,I),1)
         X21 = X21 + DDOT(NT2AMX,XJAC(1 + NT1AMX,I),1,
     *                    XJAC(1 + NT1AMX,I),1)
C
 100  CONTINUE
C
      DO 200 I = NC1VEC+1, NC1VEC + NC2VEC
C
         X12 = X12 + DDOT(NT1AMX,XJAC(1,I),1,XJAC(1,I),1)
         X22 = X22 + DDOT(NT2AMX,XJAC(1 + NT1AMX,I),1,
     *                    XJAC(1 + NT1AMX,I),1)
C
 200  CONTINUE
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
c*DECK CCLR_FDJAC
      SUBROUTINE CCLR_FDJAC(NC1VEC,NC2VEC,WORK,LWORK)
C
C-------------------------------------------------------------------------
C     Test routine for calculating the CC Jacobian by finite difference on
C     the cc-vector.
C     IF JACTST then first elements in work contains the jacobian.
C     calculates the first nc1vec rows of the C1 blocks
C     and nc2vec of the rows of the c2 blocks.
C     Written by Ove Christiansen 26-5-1994
C-------------------------------------------------------------------------
C
#include <implicit.h>
#include <maxorb.h>
#include <dummy.h>
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
      CHARACTER*10 MODEL
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
C
      IF (IPRINT .GT. 10) THEN
         CALL AROUND(' START OF CCLR_FDJAC ')
      ENDIF
C
      IF (IPRINT.GT.5) THEN
         CALL AROUND( 'IN CCLR_FDJAC: MAKING FINITE DIFF. CC JACOBIANT')
      ENDIF
C
      IF (FDJAC .AND.(IPRINT.GT.5)) THEN
         WRITE(LUPRI,'(3x,a40,//)') 
     &        'FINITE DIFF. JACOBIAN WRITTEN TO DISK'
      ENDIF
C
C----------------------------
C     Work space allocations.
C----------------------------
C
      NTAMP      = NT1AM(ISYMTR) + NT2AM(ISYMTR)
      NTAMP2     = NTAMP*(NC1VEC + NC2VEC )
      KJACOBI    = 1
      IF (JACTST) THEN
         KJACCO  = KJACOBI  + NTAMP2
      ELSE
         KJACCO  = KJACOBI
      ENDIF
      IF (FDJAC ) THEN
         KRHO1   = KJACCO   + NT1AMX + NT2AMX
      ELSE
         KRHO1   = KJACCO
      ENDIF
      KRHO2      = KRHO1    + NT1AMX
      IF (OMEGSQ) THEN
         KC1AM   = KRHO2
     *         + MAX(NT2AMX,NT2AM(ISYMOP),NT2AO(ISYMOP),NT2AOS(ISYMOP))
      ELSE
         KC1AM   = KRHO2
     *        + MAX(NT2AMX,NT2AM(ISYMOP),NT2AO(ISYMOP),2*NT2ORT(ISYMOP))
      ENDIF
      KRHO1D     = KC1AM    + NT1AMX
      KRHO2D     = KRHO1D   + NT1AM(ISYMOP)
      IF (OMEGSQ) THEN
         KC2AM  = KRHO2D
     *         + MAX(NT2AMX,NT2AM(ISYMOP),NT2AO(ISYMOP),NT2AOS(ISYMOP))
      ELSE
         KC2AM   = KRHO2D
     *        + MAX(NT2AMX,NT2AM(ISYMOP),NT2AO(ISYMOP),2*NT2ORT(ISYMOP))
      ENDIF
      KSCR1      = KC2AM      + NT2AMX
      KEND1      = KSCR1      + NT2AMX + NT1AMX
      LWRK1      = LWORK      - KEND1
      IF (IPRINT .GT. 100 ) THEN
         WRITE (LUPRI,*) ' IN CCLR_FDJAC: KJACOBI =  ',KJACOBI
         WRITE (LUPRI,*) ' IN CCLR_FDJAC: KJACCO  =  ',KJACCO
         WRITE (LUPRI,*) ' IN CCLR_FDJAC: KRHO1   =  ',KRHO1
         WRITE (LUPRI,*) ' IN CCLR_FDJAC: KRHO2   =  ',KRHO2
         WRITE (LUPRI,*) ' IN CCLR_FDJAC: KC1AM   =  ',KC1AM
         WRITE (LUPRI,*) ' IN CCLR_FDJAC: KRHO1D  =  ',KRHO1D
         WRITE (LUPRI,*) ' IN CCLR_FDJAC: KRHO2D  =  ',KRHO2D
         WRITE (LUPRI,*) ' IN CCLR_FDJAC: KC2AM   =  ',KC2AM
         WRITE (LUPRI,*) ' IN CCLR_FDJAC: KSCR1   =  ',KSCR1
         WRITE (LUPRI,*) ' IN CCLR_FDJAC: KEND1   =  ',KEND1
      ENDIF
      IF (LWRK1.LT.0 ) THEN
         WRITE (LUPRI,*) 'Too little work space in cc_jacfd '
         WRITE (LUPRI,*) 'AVAILABLE: LWORK   =  ',LWORK
         WRITE (LUPRI,*) 'NEEDED (AT LEAST)  =  ',KEND1
         CALL QUIT('TOO LITTLE WORKSPACE IN CCLR_FDJAC ')
      ENDIF
      KJACOBI2   = KJACOBI + NC1VEC*NTAMP
C
C---------------------
C     Initializations.
C---------------------
C
      CALL DZERO(WORK(KC1AM),NT1AMX)
      CALL DZERO(WORK(KC2AM),NT2AMX)
      CALL DZERO(WORK(KJACOBI),NTAMP2)
      DELTAI = 1.0D00/DELTA
      X11 = 0.0D00
      X12 = 0.0D00
      X21 = 0.0D00
      X22 = 0.0D00
      XNJ = 0.0D00
C
C------------------------------------------------
C     Read the CC reference amplitudes From disk.
C------------------------------------------------
C
      IOPT = 3
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KC1AM),WORK(KC2AM))
C
C--------------------------------------------------
C     Write the jacobian to disk if FDJAC  is done.
C--------------------------------------------------
C
      IF (FDJAC ) THEN
         LUJAC = 0
         CALL GPOPEN(LUJAC,'CCLR_JAC','UNKNOWN',' ','UNFORMATTED',
     *               IDUMMY,.FALSE.)
      ENDIF
C
C----------------------------
C     Different test options.
C----------------------------
C
c     CALL DZERO(WORK(KC1AM),NT1AMX)
c     CALL DZERO(WORK(KC2AM),NT2AMX)
c     CALL CCLR_DX(WORK(KC1AM),1.0D-2,NT1AMX)
c     CALL CCLR_DX(WORK(KC2AM),1.0D-2,NT2AMX)
C
C--------------------------------------
C     Calculate reference omega vector.
C--------------------------------------
C
      IF (IPRINT.GT.125) THEN
         CALL AROUND( 'C1AM ' )
         CALL OUTPUT(WORK(KC1AM),1,NVIRT,1,NRHFT,NVIRT,NRHFT,1,LUPRI)
         CALL AROUND( 'C2AM ' )
         CALL OUTPAK(WORK(KC2AM),NT1AMX,1,LUPRI)
      ENDIF
C
      CALL CCRHSN(WORK(KRHO1D),WORK(KRHO2D),WORK(KC1AM),WORK(KC2AM),
     *            WORK(KEND1),LWRK1)
C
C----------------------------------------------
C     Zero out components in the trial vectors.
C----------------------------------------------
C
      IF (LCOR .OR. LSEC) THEN
C
         CALL CC_CORE(WORK(KRHO1D),WORK(KRHO2D),ISYMTR)
C
      ENDIF
C
      CALL DCOPY(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1),1)
      CALL DCOPY(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2),1)
C
      IF (IPRINT.GT.125) THEN
         CALL AROUND( 'RHO1 ' )
         CALL OUTPUT(WORK(KRHO1),1,NVIRT,1,NRHFT,NVIRT,NRHFT,1,LUPRI)
         CALL AROUND( 'RHO2 ' )
         CALL OUTPAK(WORK(KRHO2),NT1AMX,1,LUPRI)
      ENDIF
C
C---------------------------------------------
C     Calculate Jacobiant by finite difference.
C---------------------------------------------
C
      DO 100 I = 1, NC1VEC
C
         TI   = SECOND()
         WORK(KC1AM +I -1) = WORK(KC1AM +I -1 ) + DELTA
         IF (LCOR .OR. LSEC) THEN
            CALL CC_CORE(WORK(KC1AM),WORK(KC2AM),ISYMTR)
         ENDIF
         CALL CCRHSN(WORK(KRHO1D),WORK(KRHO2D),WORK(KC1AM),WORK(KC2AM),
     *               WORK(KEND1),LWRK1)
C
         IF (LCOR .OR. LSEC) THEN
            CALL CC_CORE(WORK(KRHO1D),WORK(KRHO2D),ISYMTR)
         ENDIF
C
         WORK(KC1AM +I -1) = WORK(KC1AM +I -1) - DELTA
         IF (IPRINT.GT.125) THEN
            CALL AROUND( 'RHO1D ' )
            CALL OUTPUT(WORK(KRHO1D),1,NVIRT,1,NRHFT,NVIRT,NRHFT,1,
     &           LUPRI)
            CALL AROUND( 'RHO2D ' )
            CALL OUTPAK(WORK(KRHO2D),NT1AMX,1,LUPRI)
         ENDIF
         CALL DAXPY(NT1AMX,-1.0D00,WORK(KRHO1),1,WORK(KRHO1D),1)
         CALL DAXPY(NT2AMX,-1.0D00,WORK(KRHO2),1,WORK(KRHO2D),1)
         CALL CCLR_DIASCL(WORK(KRHO2D),XHALF,ISYMTR)
         CALL DSCAL(NT1AMX,DELTAI,WORK(KRHO1D),1)
         CALL DSCAL(NT2AMX,DELTAI,WORK(KRHO2D),1)
         IF (JACTST) THEN
            CALL DCOPY(NT1AMX,WORK(KRHO1D),1,
     *                 WORK(KJACOBI+NTAMP*(I-1)),1)
            CALL DCOPY(NT2AMX,WORK(KRHO2D),1,
     *                WORK(KJACOBI+NTAMP*(I-1)+NT1AMX),1)
         ENDIF
         IF (FDJAC ) THEN
            CALL DCOPY(NT1AMX,WORK(KRHO1D),1,
     *                 WORK(KJACCO),1)
            CALL DCOPY(NT2AMX,WORK(KRHO2D),1,
     *                WORK(KJACCO+NT1AMX),1)
            WRITE(LUJAC) (WORK(KJACCO +J-1), J = 1, NTAMP)
         ENDIF
         X11 = X11 + DDOT(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1D),1)
         X21 = X21 + DDOT(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2D),1)
C
         TI   = SECOND() - TI
         IF (IPRINT.GT.5 ) THEN
            WRITE(LUPRI,*) '  '
            WRITE(LUPRI,*) 'FDJAC ROW NR. ',I,' DONE IN ',TI,' SEC.'
         ENDIF
C
 100  CONTINUE
C
C----------------------------------------------------------------
C     Loop over T2 amplitudes. Take care of diagonal t2 elements
C     is in a different convention in the energy code.
C     Factor 1/2 from right , and factor 2 from left.
C----------------------------------------------------------------
C
      DO 200 NAI = 1, NT1AMX
        DO 300 NBJ = 1, NAI
         I = INDEX(NAI,NBJ)
C
         IF (I.LE.NC2VEC) THEN
C
           TI   = SECOND()
           WORK(KC2AM + I -1) = WORK(KC2AM+I -1) + DELTA
           IF (LCOR .OR. LSEC) THEN
             CALL CC_CORE(WORK(KC1AM),WORK(KC2AM),ISYMTR)
           ENDIF
           CALL CCRHSN(WORK(KRHO1D),WORK(KRHO2D),WORK(KC1AM),
     *               WORK(KC2AM),WORK(KEND1),LWRK1)
C
           IF (LCOR .OR. LSEC) THEN
              CALL CC_CORE(WORK(KRHO1D),WORK(KRHO2D),ISYMTR)
           ENDIF
           WORK(KC2AM + I -1) = WORK(KC2AM+I -1) - DELTA
C
           IF (IPRINT.GT.125) THEN
              CALL AROUND( 'RHO1D ' )
              CALL OUTPUT(WORK(KRHO1D),1,NVIRT,1,NRHFT,NVIRT,NRHFT,1,
     &             LUPRI)
              CALL AROUND( 'RHO2D ' )
              CALL OUTPAK(WORK(KRHO2D),NT1AMX,1,LUPRI)
           ENDIF
C
           CALL DAXPY(NT1AMX,-1.0D00,WORK(KRHO1),1,WORK(KRHO1D),1)
           CALL DAXPY(NT2AMX,-1.0D00,WORK(KRHO2),1,WORK(KRHO2D),1)
           CALL CCLR_DIASCL(WORK(KRHO2D),XHALF,ISYMTR)
           CALL DSCAL(NT1AMX,DELTAI,WORK(KRHO1D),1)
           CALL DSCAL(NT2AMX,DELTAI,WORK(KRHO2D),1)
C
           IF (NAI.EQ.NBJ) THEN
C
              CALL DSCAL(NT1AMX,2.0D00,WORK(KRHO1D),1)
              CALL DSCAL(NT2AMX,2.0D00,WORK(KRHO2D),1)
C
           ENDIF
C
           IF (JACTST) THEN
              CALL DCOPY(NT1AMX,WORK(KRHO1D),1,
     *                 WORK(KJACOBI2+NTAMP*(I-1)),1)
              CALL DCOPY(NT2AMX,WORK(KRHO2D),1,
     *                 WORK(KJACOBI2+NTAMP*(I-1)+NT1AMX),1)
           ENDIF
C
           IF (FDJAC ) THEN
              CALL DCOPY(NT1AMX,WORK(KRHO1D),1,
     *                   WORK(KJACCO ),1)
              CALL DCOPY(NT2AMX,WORK(KRHO2D),1,
     *                  WORK(KJACCO +NT1AMX),1)
              WRITE(LUJAC) (WORK(KJACCO +J-1), J = 1, NTAMP)
           ENDIF
C
           X12 = X12 + DDOT(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1D),1)
           X22 = X22 + DDOT(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2D),1)
           TI   = SECOND() - TI
           IF (IPRINT.GT.5 ) THEN
              WRITE(LUPRI,*) '  '
              WRITE(LUPRI,*) 'FDJAC ROW NR. ',I+NT1AMX,
     *                  ' DONE IN ',TI,' SEC.'
           ENDIF
C
         ENDIF
C
 300    CONTINUE
 200  CONTINUE
C
      WRITE(LUPRI,*) '    '
      WRITE(LUPRI,*) '**  FINITE DIFF WITH DELTA ',DELTA, '**'
      WRITE(LUPRI,*) '    '
      IF (FDJAC) CALL GPCLOSE(LUJAC,'KEEP')
      IF ((IPRINT .GT. 40).AND.(JACTST)) THEN
         CALL AROUND( 'FINITE DIFF. CC JACOBIANT - 11 & 21 PART  ' )
         CALL OUTPUT(WORK(KJACOBI),1,NTAMP,1,NC1VEC,NTAMP,NC1VEC,1,
     &        LUPRI)
         CALL AROUND( 'FINITE DIFF. CC JACOBIANT - 12 & 22 PART  ' )
         CALL OUTPUT(WORK(KJACOBI+NTAMP*NC1VEC),1,NTAMP,1,NC2VEC,
     *               NTAMP,NC2VEC,1,LUPRI)
      ENDIF
      IF (JACTST) THEN
         XNJ = X11 + X12 + X21 + X22
         WRITE(LUPRI,*) '  '
         WRITE(LUPRI,*) ' NORM OF FIN. DIFF. JAC.', SQRT(XNJ)
         WRITE(LUPRI,*) '  '
         WRITE(LUPRI,*) 'NORM OF 11 PART OF FIN. DIFF. JAC. ', SQRT(X11)
         WRITE(LUPRI,*) 'NORM OF 21 PART OF FIN. DIFF. JAC. ', SQRT(X21)
         WRITE(LUPRI,*) 'NORM OF 12 PART OF FIN. DIFF. JAC. ', SQRT(X12)
         WRITE(LUPRI,*) 'NORM OF 22 PART OF FIN. DIFF. JAC. ', SQRT(X22)
      ENDIF
C
      IF (IPRINT .GT. 10) THEN
         CALL AROUND(' END OF CCLR_FDJAC ')
      ENDIF
C
      RETURN
      END
c*DECK CCLR_DUMTRR
      SUBROUTINE CCLR_DUMTRR(CVEC1,RHO1,RHO2,WORK,LWORK)
C
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
      PARAMETER ( TWO = 2.0D00,XHALF=0.5D00 )
      DIMENSION CVEC1(*),WORK(LWORK)
      DIMENSION RHO1(*),RHO2(*)
C
      IF (IPRINT .GT. 10) THEN
         CALL AROUND(' START OF CCLR_DUMTRR ')
      ENDIF
C
C----------------------------------
C        Print transformed vectors.
C----------------------------------
C
      IF (IPRINT .GT.110) THEN
         CALL AROUND(' CCLR_DUMTRR: vectors to be transformed ')
         CALL CC_PRP(CVEC1,CVEC1(1+NT1AM(ISYMTR)),ISYMTR,1,1)
      ENDIF
C
      NTAMP = NT1AM(ISYMTR) + NT2AM(ISYMTR)
C
C------------------------
C     Open Jacobian file.
C------------------------
C
      LUJAC = 0
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
C
C----------------------------------
C        Print transformed vectors.
C----------------------------------
C
      IF (IPRINT .GT.110) THEN
         CALL AROUND(' CCLR_DUMTRR: transformed vectors ')
         CALL CC_PRP(RHO1,RHO2,ISYMTR,1,1)
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
CCH   DIMENSION T1AM(NT1AM(ISYM),NT1),T2AM(NT2AM(ISYM),NT2)
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
CCH         KOFF = IT2AM(ISYMAI,ISYMBJ) + 1
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
CCH   DIMENSION T1AM(NT1AM(ISYM),NT1),T2AM(NT2SQ(ISYM),NT2)
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
C
      LOGICAL OPNSTA
      DIMENSION WORK(LWORK)
C
#include <leinf.h>
C
C---------------------------------------------
C     Header of Excitation Energy calculation.
C---------------------------------------------
C
      CALL TIMER('START ',TIMEIN,TIMOUT)
      IF (IPRINT .GT. 0) THEN
         WRITE (LUPRI,'(A,/)') '  '
         WRITE (LUPRI,'(A,/)')
     *    '1 <<<<<<<<<< Output from CCLR_FDEXCI >>>>>>>>> '
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
C
C----------------------------
C     Alloction of Workspace.
C----------------------------
C
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
      LWRK1   = LWORK - KWRK1
      IF (LWRK1.LT. 0)
     &     CALL QUIT('TOO LITTLE WORK SPACE IN CCLR_FDEXCI ')
C
C---------------------
C     Readin Jacobian.
C---------------------
C
      CALL CCLR_JACIN(WORK(KSOLEQ))
C
      IF (IPRINT .GT. 100) THEN
C
         CALL AROUND( 'FIN. DIF. CC JACOBIANT - A*C1 PART ' )
         CALL OUTPUT(WORK(KSOLEQ),1,NTAMP,1,NT1AMX,MAXRED,NTAMP,1,LUPRI)
         CALL AROUND( 'FIN. DIF. CC JACOBIANT - A*C2 PART ' )
         CALL OUTPUT(WORK(KSOLEQ +MAXRED*NT1AMX),1,NTAMP,1,NT2AMX,
     *              MAXRED,NT2AMX,1,LUPRI)
C
      ENDIF
C
C---------------------------
C     Solve for eigenvalues.
C---------------------------
C
      CALL DCOPY(MAXRED*MAXRED,WORK(KSOLEQ),1,WORK(KAMAT),1)
      CALL DCOPY(MAXRED*MAXRED,WORK(KSOLEQ),1,WORK(KAMAT2),1)
      CALL RG(MAXRED,NREDH,WORK(KAMAT),WORK(KEIVAL),WORK(KWI),MATZ,
     *        WORK(KSOLEQ),WORK(KIV1),WORK(KFV1),IERR)
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
      CALL RGORD(MAXRED,NREDH,WORK(KEIVAL),WORK(KWI),WORK(KSOLEQ))
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
         CALL CC_PRAM(WORK(KSOLEQ + (I-1)*MAXRED),PT1,1)
  500 CONTINUE
C
      RETURN
      END
c*DECK CCLR_JACIN
      SUBROUTINE CCLR_JACIN(FDJACO)
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
C
      DIMENSION FDJACO(MAXRED,MAXRED)
C
      NTAMP = NT1AMX + NT2AMX
C
C------------------------
C     Open Jacobian file.
C------------------------
C
      LUJAC = 0
      CALL GPOPEN(LUJAC,'CCLR_JAC','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUJAC)
C
      DO 100 I = 1, NTAMP
C
         READ(LUJAC) (FDJACO(J,I),J=1,NTAMP)
C
  100 CONTINUE
C
      CALL GPCLOSE(LUJAC,'KEEP')
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
C     Purpose: To zero core and secondary part og
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
      SUBROUTINE CC_PRAM(CAM,PT1,ISYMTR)
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
C
C
      DIMENSION CAM(*)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C------------------------
C     Add up the vectors.
C------------------------
C
      C1NOSQ = 0
      C2NOSQ = 0
      KC1 = 1
      KC2 = KC1 + NT1AM(ISYMTR)
      C1NOSQ = DDOT(NT1AM(ISYMTR),CAM(KC1),1,CAM(KC1),1)
      IF (.NOT. CCS) C2NOSQ = DDOT(NT2AM(ISYMTR),CAM(KC2),1,CAM(KC2),1)
      CNOSQ  = C1NOSQ + C2NOSQ
C
      IF (.NOT. (CCS.OR.CCP2)) THEN
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
         WRITE(LUPRI,'(/10X,A,10X,F10.4)')
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
C     Loop through Doublee excitation vector.
C     If not ccs or ccp2
      IF (.NOT. ( CCS .OR. CCP2 )) THEN
C--------------------------------------------
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
      WRITE(LUPRI,'(//10X,A43,1X,F8.6)')
     *     'Printed all single excitations greater than',THR1
      IF (.NOT. (CCS.OR.CCP2)) THEN
         WRITE(LUPRI,'(//10X,A43,1X,F8.6)')
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
#include <implicit.h>
      DIMENSION VEC(*)
      CHARACTER FIL*(*)
C
      IOFF = 1 + LLEN*(NR-1)
C
      IF (LEN .GT. 0) THEN
         CALL GETWA2(LU,FIL,VEC,IOFF,LEN)
C
      ENDIF
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
      DIMENSION VEC(*)
      CHARACTER FIL*(*)
C
      IOFF = 1 + LLEN*(NR-1)
C
      IF (LEN .GT. 0) THEN
         CALL PUTWA2(LU,FIL,VEC,IOFF,LEN)
      END IF
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
C       CCLHST via CCLR_TRR was switched of since quite some time...
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
C
C-------------------
C     Open jac file.
C-------------------
C
      LUJAC = 0
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
