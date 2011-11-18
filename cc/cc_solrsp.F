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
*=====================================================================*
      SUBROUTINE CC_SOLRSP(WORK,LWORK,APROXR12)
C---------------------------------------------------------------------*
C
C     Purpose: call equation driver for the various types of
C              response equations we have to solve
C
C      Written by Ove Christiansen October 1996
C      Restructured by Christof Haettig, 1997/1998
C      PL1. Sonia 2000
C      SLV98,OC: Self consistent solution of 1-order response equations.
C
C---------------------------------------------------------------------*
C
      IMPLICIT NONE
#include <priunit.h>
#include <maxorb.h>
#include <dummy.h>
#include <iratdef.h>
#include <cclr.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <ccexci.h>
#include <cclrmrsp.h>
#include <ccer1rsp.h>
#include <ccer2rsp.h>
#include <ccel1rsp.h>
#include <ccel2rsp.h>
#include <ccr1rsp.h>
#include <cco1rsp.h>
#include <ccx1rsp.h>
#include <ccl1rsp.h>
#include <ccr2rsp.h>
#include <ccl2rsp.h>
#include <cco2rsp.h>
#include <ccx2rsp.h>
#include <ccr3rsp.h>
#include <ccl3rsp.h>
#include <cco3rsp.h>
#include <ccx3rsp.h>
#include <ccr4rsp.h>
#include <ccl4rsp.h>
#include <cco4rsp.h>
#include <ccx4rsp.h>
#include <ccn2rsp.h>
#include <ccrc1rsp.h>
#include <cclc1rsp.h>
#include <cccr2rsp.h>
#include <cccl2rsp.h>
#include <ccco2rsp.h>
#include <cccx2rsp.h>
#include <ccexgr.h>
#include <leinf.h>
#include <ccpl1rsp.h>
#include <ccslvinf.h>
Cholesky
#include <ccdeco.h>
Cholesky
C
C
      LOGICAL LOCDBG, lold
      PARAMETER (LOCDBG = .FALSE., lold = .false.)
C
      CHARACTER*(*) APROXR12
      CHARACTER*8 FC1AM,FC2AM,FRHO1,FRHO2
      PARAMETER (FC1AM='CCR_C1AM',FC2AM='CCR_C2AM')
      PARAMETER (FRHO1='CCR_RHO1',FRHO2='CCR_RHO2')

      INTEGER LWORK

#if defined (SYS_CRAY)
      REAL TOL, ONE, TWO, WORK(LWORK), DDOT, X1, X2, FREK
#else
      DOUBLE PRECISION TOL, ONE, TWO, WORK(LWORK), DDOT, X1, X2, FREK
#endif
      LOGICAL LDUMMY
      INTEGER ISYM0, IZERO

      PARAMETER(TOL=1.0D-12, ISYM0 = 1, LDUMMY = .FALSE.)
      PARAMETER(ONE = 1.0D0, IZERO = 0, TWO = 2.0D0)

      CHARACTER MODEL*10,MODEL1*10
      CHARACTER*3 LIST 
      INTEGER MAXFDIM
      PARAMETER (MAXFDIM = 500)
      INTEGER IFTRAN(3,MAXFDIM)
      INTEGER ORDER, I1DXORD, IOPREL, NSTAT, ISIDE, IDXR1, IOPT,
     &        ILSTNR, ISYM, KT1, KT2, KEND1, LEND1, IDXM, IDXRE, IDXF,
     &        NTAMP, KETA, LWRK1, KBEVC, KEND2, LWRK2, IDXRC, KVEC1,
     &        KVEC2, NCAU, IVEC, INUM, IDXLC, IDXL2, IDXR2, IDXL3, 
     &        IDXR3, NLRCLBLSAVE, KTR12
C
      LOGICAL R1HFSKIP, O1SKIP, X1SKIP
      SAVE    R1HFSKIP, O1SKIP, X1SKIP
      DATA R1HFSKIP / .FALSE. /
      DATA O1SKIP   / .FALSE. /
      DATA X1SKIP   / .FALSE. /

* external functions:
      INTEGER IR1TAMP, ILRCAMP, ILC1AMP, ICR2AMP, IRHSCR2, IR2TAMP,
     &        ICL2AMP, IETACL2, IR3TAMP, ILSTSYM

C-------------------------------------
C     Some initializations:
C-------------------------------------
C
      CALL QENTER('CC_SOLRSP')
C
      MODEL = 'CCSD      '
      IF (CCS) MODEL = 'CCS       '
      IF (CC2) MODEL = 'CC2       '
C
      IF (.NOT.(CCS.OR.CC2.OR.CCSD.OR.CC3)) THEN
         WRITE(LUPRI,'(A)') ' CC_SOLRSP: Do not want to calculate '
     *           //'anything else than CCS, CC2 and CCSD properties '
         CALL QUIT('Model not CCS or CC2 or CCSD in CC_SOLRSP')
      ENDIF
 
      IF (CHOINT .AND. (.NOT. CC2)) THEN
         WRITE(LUPRI,'(A,A,A)')
     &   ' CC_SOLRSP: Cholesky decomposed integrals not implemented',
     &   ' for ',MODEL
         CALL QUIT('Cholesky '//MODEL//' not implemented in CC_SOLRSP')
      ENDIF

      WRITE (LUPRI,'(4X,A,/)') '  '
      WRITE (LUPRI,'(4X,A)')
     *'*********************************************************'//
     *'**********'
      WRITE (LUPRI,'(4X,A)')
     *'*          SOLVING COUPLED CLUSTER RESPONSE EQUATIONS    '//
     *'         *'
      WRITE (LUPRI,'(4X,A,///)')
     *'*********************************************************'//
     *'**********'

      IF (CHOINT) THEN
         WRITE(LUPRI,'(4X,A,///)')
     &   '         (Calculation using Cholesky decomposed integrals)'
      ELSE
         WRITE(LUPRI,'(A)') ' '
      ENDIF

C======================================================================

      WRITE(LUPRI,'(1X,A1,70("="),A1)') '+','+'
      WRITE(LUPRI,'(1X,A1,21X,A,20X,A1)')
     &  '!','RHS & ETA VECTORS TO COMPUTE:','!'
      WRITE(LUPRI,'(1X,A1,70("="),A1)') '+','+'

      WRITE(LUPRI, '(1X,A1,1X,A4,1X,A1,1X,A6,1X,A1,1X,A52,1X,A1)')
     &   '|','TYPE','|','# VEC.','|',
     &   'NEEDED FOR:                                        ','|'

      WRITE(LUPRI,'(1X,A1,70("-"),A1)') '+','+'

      IF (NO1LBL.GT.0) WRITE(LUPRI,2) '|',' O1 ','|',   NO1LBL,'|',
     &   'first-order amplitude equations                    ','|'
      IF (NX1LBL.GT.0) WRITE(LUPRI,2) '|',' X1 ','|',   NX1LBL,'|',
     &   'first-order multipliers equations                  ','|'
      IF (NO2LBL.GT.0) WRITE(LUPRI,2) '|',' O2 ','|',   NO2LBL,'|',
     &   'second-order amplitude equations                   ','|'
      IF (NX2LBL.GT.0) WRITE(LUPRI,2) '|',' X2 ','|',   NX2LBL,'|',
     &   'second-order multipliers equations                 ','|'
      IF (NCO2LBL.GT.0) WRITE(LUPRI,2) '|','CO2','|',  NCO2LBL,'|',
     &   'second-order right Cauchy equations                ','|'
      IF (NCX2LBL.GT.0) WRITE(LUPRI,2) '|','CX2','|',  NCX2LBL,'|',
     &   'second-order left Cauchy equations                 ','|'
      IF (NO3LBL.GT.0) WRITE(LUPRI,2) '|',' O3 ','|',   NO3LBL,'|',
     &   'third-order amplitude equations                    ','|'
      IF (NX3LBL.GT.0) WRITE(LUPRI,2) '|',' X3 ','|',   NX3LBL,'|',
     &   'third-order multiplier equations                   ','|'
      IF (NXGRST.GT.0) WRITE(LUPRI,2) '|',' BE ','|',   NXGRST,'|',
     &   'zeroth-order excited-state Lagrange multiplier     ','|'
      IF (NQRN2.GT.0)  WRITE(LUPRI,2) '|',' BR ','|',    NQRN2,'|',
     &   'zeroth-order excited-state transition Lagr. multip.','|'
      IF (NER1LBL.GT.0) WRITE(LUPRI,2)'|',' EO1','|',  NER1LBL,'|',
     &   'first-order excited state right response equations ','|'
      IF (NEL1LBL.GT.0) WRITE(LUPRI,2)'|',' EX1','|',  NEL1LBL,'|',
     &   'first-order excited state left response equations  ','|'
      IF (NEL2LBL.GT.0) WRITE(LUPRI,2)'|',' EO2','|',  NER2LBL,'|',
     &   'second-order excited state right response equations','|'
      IF (NER2LBL.GT.0) WRITE(LUPRI,2)'|',' EX2','|',  NEL2LBL,'|',
     &   'second-order excited state left response equations ','|'

      WRITE(LUPRI,'(1X,A1,70("="),A1,///)') '+','+'

C======================================================================

      WRITE(LUPRI,'(1X,A1,70("="),A1)') '+','+'
      WRITE(LUPRI,'(1X,A1,21X,A,20X,A1)')
     &  '!','  LINEAR EQUATIONS TO SOLVE: ','!'
      WRITE(LUPRI,'(1X,A1,70("="),A1)') '+','+'

      WRITE(LUPRI, '(1X,A1,1X,A4,1X,A1,1X,A6,1X,A1,1X,A52,1X,A1)')
     &   '|','TYPE','|','# VEC.','|',
     &   'EQUATION:                                          ','|'

      WRITE(LUPRI,'(1X,A1,70("-"),A1)') '+','+'

      IF (NLRTHFLBL.GT.0)WRITE(LUPRI,2)'|',' R1 ','|',NLRTHFLBL,'|',
     &   'first-order orbtial response                       ','|'
      IF (NLRTLBL.GT.0) WRITE(LUPRI,2) '|',' R1 ','|', NLRTLBL,'|',
     &   'first-order amplitude response                     ','|'
      IF (NLRZLBL.GT.0) WRITE(LUPRI,2) '|',' L1 ','|', NLRZLBL,'|',
     &   'first-order multiplier response                    ','|'
      IF (NLRCLBL.GT.0) WRITE(LUPRI,2) '|',' RC ','|', NLRCLBL,'|',
     &   'first-order right Cauchy vectors                   ','|'
      IF (NLC1LBL.GT.0) WRITE(LUPRI,2) '|',' LC ','|', NLC1LBL,'|',
     &   'first-order left Cauchy vectors                    ','|'
      IF (NR2TLBL.GT.0) WRITE(LUPRI,2) '|',' R2 ','|', NR2TLBL,'|',
     &   'second-order amplitude response                    ','|'
      IF (NCR2LBL.GT.0) WRITE(LUPRI,2) '|',' CR2','|', NCR2LBL,'|',
     &   'second-order right Cauchy vectors                  ','|'
      IF (NL2LBL.GT.0)  WRITE(LUPRI,2) '|',' L2 ','|',  NL2LBL,'|',
     &   'second-order multiplier response                   ','|'
      IF (NCL2LBL.GT.0) WRITE(LUPRI,2) '|',' CL2','|', NCL2LBL,'|',
     &   'second-order left Cauchy vectors                   ','|'
      IF (NR3TLBL.GT.0) WRITE(LUPRI,2) '|',' R3 ','|', NR3TLBL,'|',
     &   'third-order amplitude response                     ','|'
      IF (NL3LBL.GT.0)  WRITE(LUPRI,2) '|',' L3 ','|',  NL3LBL,'|',
     &   'third-order multiplier response                    ','|'
      IF (NLRM.GT.0)    WRITE(LUPRI,2) '|',' M1 ','|',    NLRM,'|',
     &   'first-order ground - excited state multipliers     ','|'
      IF (NXGRST.GT.0)  WRITE(LUPRI,2) '|',' E0 ','|',  NXGRST,'|',
     &   'zeroth-order excited-state Lagrange multiplier     ','|'
      IF (NQRN2.GT.0)   WRITE(LUPRI,2) '|',' N2 ','|',   NQRN2,'|',
     &   'zeroth-order excited-state transition Lagr. multip.','|'
      IF (NER1LBL.GT.0) WRITE(LUPRI,2) '|',' ER1','|', NER1LBL,'|',
     &   'first-order right excited state response           ','|'
      IF (NEL1LBL.GT.0) WRITE(LUPRI,2) '|',' EL1','|', NEL1LBL,'|',
     &   'first-order left excited state response            ','|'
      IF (NER2LBL.GT.0) WRITE(LUPRI,2) '|',' ER2','|', NER2LBL,'|',
     &   'second-order right excited state response          ','|'
      IF (NEL2LBL.GT.0) WRITE(LUPRI,2) '|',' EL2','|', NEL2LBL,'|',
     &   'second-order left excited state response           ','|'
      IF (NPL1LBL.GT.0) WRITE(LUPRI,2) '|',' PL1','|', NPL1LBL,'|',
     &   'first-order projected multiplier response          ','|'

      WRITE(LUPRI,'(1X,A1,70("="),A1,///)') '+','+'

C======================================================================

      WRITE(LUPRI,'(1X,A1,70("="),A1)') '+','+'
      WRITE(LUPRI,'(1X,A1,16X,A,16X,A1)')
     &  '!','F MATRIX TRANSFORMATIONS TO COMPUTE:  ','!'
      WRITE(LUPRI,'(1X,A1,70("="),A1)') '+','+'

      WRITE(LUPRI, '(1X,A1,1X,A4,1X,A1,1X,A6,1X,A1,1X,A52,1X,A1)')
     &   '|','TYPE','|','# VEC.','|',
     &   'TRANSFORMED:                                       ','|'

      WRITE(LUPRI,'(1X,A1,70("-"),A1)') '+','+'

      IF (NLRTLBL.GT.0) WRITE(LUPRI,2) '|',' F1 ','|', NLRTLBL,'|',
     &   'first-order amplitude response (R1) vector         ','|'
      IF (NLC1LBL.GT.0) WRITE(LUPRI,2) '|',' FC ','|', NLC1LBL,'|',
     &   'first-order right Cauchy (RC) vector               ','|'
      IF (NLRM.GT.0)    WRITE(LUPRI,2) '|',' FR ','|',    NLRM,'|',
     &   'right eigenvectors (RE)                            ','|'
      IF (NL2LBL.GT.0)  WRITE(LUPRI,2) '|',' F2 ','|',  NL2LBL,'|',
     &   'second-order amplitude response (R2) vector        ','|'
      IF (NCL2LBL.GT.0) WRITE(LUPRI,2) '|',' CF2','|', NCL2LBL,'|',
     &   'second-order right Cauchy (CR2) vector             ','|'
      IF (NL3LBL.GT.0)  WRITE(LUPRI,2) '|',' F3 ','|',  NL3LBL,'|',
     &   'third-order amplitude response (R3) vector         ','|'

      WRITE(LUPRI,'(1X,A1,70("="),A1,///)') '+','+'

C======================================================================
C     Solve first-order CPHF kappa equations.
C     (since these equations are independent of the CC amplitudes,
C      we solve them only the first time we enter this routine)
C======================================================================
      IF ( .NOT. R1HFSKIP ) THEN

         CALL CC_CPHF('R1 ',LRTHFLBL,IDUMMY,IDUMMY,DUMMY,
     &                ISYLRTHF,FRQLRTHF,IDUMMY,NLRTHFLBL,MAXTLBL,
     &                WORK,LWORK)

         R1HFSKIP = .TRUE.

      END IF

C======================================================================
C    Calculate effective Fock matrices from Kappa^(1) one-index
C    transformed integrals:
C======================================================================

C     IF (NSYM.EQ.1) THEN
        I1DXORD = 1
        IOPREL  = 0
        CALL CC_GRAD2(I1DXORD,WORK,LWORK)
C       CALL CCEFFFOCK(I1DXORD,IOPREL,WORK,LWORK)
C     ELSE
C       WRITE (LUPRI,*) 'Warning: no 1-index transformed eff. Fock '
C       WRITE (LUPRI,*) 'matrices calculated... '
C       WRITE (LUPRI,*) 
C    *             'RELAXED second-order properties will be wrong!.'
C     END IF

C======================================================================
C     precalculate right-hand-sides vectors "O1" and the eta vectors 
C     for those first-order cluster amplitude and Lagrangian multiplier
C     equations for which this is neccessary, because
C     they involve orbital relaxaton and/or derivative integrals
C     precalculate first-order chi vectors "X1" for those first-order
C======================================================================
Ctmp
      CALL GPINQ('SKIP_O1','EXIST',O1SKIP)
Ctmp
      IF ( .NOT. O1SKIP ) THEN
 
          CALL CCRHSVEC('O1 ',LBLO1,IDUMMY,IDUMMY,DUMMY,ISYO1,FRQO1,
     &                  LORXO1,IDUMMY,NO1LBL,MAXO1LBL,0,
     &                  WORK,LWORK)
 
      ELSE
        WRITE (LUPRI,*)
     &        'Skipped calculation of O1 right hand side vectors.'
        O1SKIP = .FALSE.
        CALL SYSTEM('rm SKIP_O1')
      END IF
C
C------------------------------------------------------------------
C     Start loop over self-consistent solvent solution of R1/L1
C     CCSLV98, OC
C------------------------------------------------------------------
C
      ICCSLIT = 0
C
  999 CONTINUE
C======================================================================
C     Solve right first-order cluster amplitude response equations.
C======================================================================
      LIST  = 'R1 '
      NSTAT = 0
      ORDER = 1
      ISIDE = +1

      IF ( .NOT. R1SKIP ) THEN
         CALL CC_SOLDRV(LIST,NSTAT,ORDER,ISIDE,APROXR12,
     &                  IDUMMY, IDUMMY, DUMMY, LDUMMY,
     &                  ISYLRT, LRTLBL, FRQLRT, IDUMMY,
     &                  ISYOFT, NLRTLBL, MAXTLBL, 
     &                  WORK, LWORK)

      ELSE
         WRITE (LUPRI,*) 'Skipped first-order response (R1) equations'
         R1SKIP = .FALSE.
      END IF

      IF (DEBUG .OR. LOCDBG. OR. IPRINT.GE.5) THEN
         CALL CC_TSTLST(LIST,NSTAT,ORDER,ISIDE,
     &                  IDUMMY, IDUMMY, DUMMY, LDUMMY,
     &                  ISYLRT, LRTLBL, FRQLRT, IDUMMY,
     &                  ISYOFT, NLRTLBL, MAXTLBL, WORK, LWORK)
      END IF

C======================================================================
C     precalculate F * R1 transformation required in CCLR section
C     and to for the first-order coupled cluster multipliers L1
C======================================================================
      LIST  = 'F1 '
      NSTAT = 0
      ORDER = 1
      ISIDE = +1

      IF ( .NOT. F1SKIP ) THEN

         IF (NLRTLBL .GT. 0) THEN

            IF (NLRTLBL .GT. MAXFDIM) THEN
              WRITE (LUPRI,*)
     &              'ERROR: Dimension MAXFDIM to small in CC_SOLRSP.'
              CALL QUIT('Dimension MAXFDIM to small in CC_SOLRSP.')
            END IF
C
C           set up list of F matrix transformations:
C
            DO IDXR1 = 1, NLRTLBL
              IFTRAN(1,IDXR1) = 0
              IFTRAN(2,IDXR1) = IDXR1
              IFTRAN(3,IDXR1) = IDXR1
            END DO
C
C           call CC_FMATRIX to do the F * R1 transformations:
C
            IOPT = 3
            CALL CC_FMATRIX(IFTRAN, NLRTLBL, 'L0 ', 'R1 ', IOPT, 'F1 ',
     &                      IDUMMY, DUMMY, 0, WORK, LWORK )

         END IF

      ELSE
        WRITE (LUPRI,*) 'Skipped of F x R1 transformations'
        F1SKIP = .FALSE.
      END IF 


      IF (DEBUG .OR. LOCDBG. OR. IPRINT.GE.5) THEN
         CALL CC_TSTLST(LIST,NSTAT,ORDER,ISIDE,
     &                  IDUMMY, IDUMMY, DUMMY, LDUMMY,
     &                  ISYLRT, LRTLBL, FRQLRT, IDUMMY,
     &                  ISYOFT, NLRTLBL, MAXTLBL, WORK, LWORK)
      END IF

C======================================================================
C     For Cholesky CC2, calculate the "right" first-order density
C     matrix.
C======================================================================

      IF (CHOINT .AND. CC2) THEN
         IF (.NOT. D01SKIP) THEN

            ORDER = 1
            ISIDE = 1

            DO IDXR1 = 1, NLRTLBL
               IFTRAN(1,IDXR1) = 0
               IFTRAN(2,IDXR1) = IDXR1
               IFTRAN(3,IDXR1) = IDXR1
            ENDDO

            CALL CC_CHOPDEN(IFTRAN,NLRTLBL,'L0 ','R1 ','d01',ORDER,
     &                      ISIDE,'DENSITY_01',.TRUE.,WORK,LWORK)

         ELSE

            WRITE (LUPRI,*) 'Skipped right first-order density (d01)',
     &                      ' calculation'
            D01SKIP = .FALSE.

         ENDIF
      ENDIF

C======================================================================
C     Solve first-order Lagrange multiplier response equations.
C======================================================================
      LIST  = 'L1 '
      NSTAT = 0
      ORDER = 1
      ISIDE = -1

      IF ( .NOT. L1SKIP ) THEN
        CALL CC_SOLDRV(LIST,NSTAT,ORDER,ISIDE,APROXR12,
     &                 IDUMMY, IDUMMY, DUMMY, LDUMMY,
     &                 ISYLRZ, LRZLBL, FRQLRZ, IDUMMY,
     &                 ISYOFZ, NLRZLBL, MAXZLBL, 
     &                 WORK, LWORK)

      ELSE
        WRITE (LUPRI,*)
     &        'Skipped first-order left response (L1) equations'
        L1SKIP = .FALSE.
      END IF

      IF (DEBUG .OR. LOCDBG. OR. IPRINT.GE.5) THEN
        CALL CC_TSTLST(LIST,NSTAT,ORDER,ISIDE,
     &                 IDUMMY, IDUMMY, DUMMY, LDUMMY,
     &                 ISYLRZ, LRZLBL, FRQLRZ, IDUMMY,
     &                 ISYOFZ, NLRZLBL, MAXZLBL, WORK, LWORK)
      END IF

C
C------------------------------------------------------------------
C     If solvent calculation then check vectors are solved
C     self-consistently otherwise continue loop.
C     CCSLV98, OC
C------------------------------------------------------------------
C
      IF (CCSLV) THEN
C
         IF (CC3) CALL QUIT('No solvent for CC3...')
C
         LSLTCVG = .TRUE.
         LSLLCVG = .TRUE.
C
C        Test if there exist vectors that are not converged.
C
         WRITE(LUPRI,'(//,1X,A)') 'Testing R1 list: '
         WRITE(LUPRI,'(/,1X,A)')
     * 'List nr.  Label     Symmetry      Frequency      T1-norm     '
     * //'   T2-norm'
         DO 997  ILSTNR = 1, NLRTLBL
C
            ISYM = ISYLRT(ILSTNR)
            KT1  = 1
            KT2  = KT1  + NT1AM(ISYM)
            IF (.NOT. CCS) THEN
               KEND1= KT2  + NT2AM(ISYM)
            ELSE
               KEND1= KT2
            ENDIF
            LEND1= LWORK - KEND1
            IF (LEND1 .LE. 0 ) CALL QUIT(
     *          'Insufficient workspace in CC_SOLRSP ')
            IOPT  = 3
            CALL CC_RDRSP('R1',ILSTNR,ISYM,IOPT,MODEL1,WORK(KT1),
     *                    WORK(KT2))
            X1  = DDOT(NT1AM(ISYM),WORK(KT1),1,WORK(KT1),1)
            X2  = 0.0D0
            IF (.NOT.CCS) THEN
               X2  = DDOT(NT2AM(ISYM),WORK(KT2),1,WORK(KT2),1)
            ENDIF
            WRITE(LUPRI,'(1X,I5,5X,A8,2X,I4,6X,F15.10,F15.10,F15.10)')
     *           ILSTNR,LRTLBL(ILSTNR),ISYM,FRQLRT(ILSTNR),X1,X2
C
            XLNCCCU = X1+X2
            XLNCCPR = XLRT(ILSTNR)
            IF (ABS(XLNCCPR-XLNCCCU).GT.CVGTSOL) LSLTCVG = .FALSE.
            IF (IPRINT.GT.16)  WRITE(LUPRI,*) ' LSLTCVG = ',LSLTCVG
            XLRT(ILSTNR) = XLNCCCU
C
 997     CONTINUE
C
         WRITE(LUPRI,'(//,1X,A)') 'Testing L1 list: '
         WRITE(LUPRI,'(/,1X,A)')
     * 'List nr.  Label     Symmetry      Frequency      T1-norm     '
     * //'   T2-norm'
         DO 998  ILSTNR = 1, NLRZLBL
C
            ISYM = ISYLRZ(ILSTNR)
            KT1  = 1
            KT2  = KT1  + NT1AM(ISYM)
            IF (.NOT. CCS) THEN
               KEND1= KT2  + NT2AM(ISYM)
            ELSE
               KEND1= KT2
            ENDIF
            LEND1= LWORK - KEND1
            IF (LEND1 .LE. 0 ) CALL QUIT(
     *          'Insufficient workspace in CC_SOLRSP ')
            IOPT  = 3
            CALL CC_RDRSP('L1',ILSTNR,ISYM,IOPT,MODEL1,WORK(KT1),
     *                    WORK(KT2))
            X1  = DDOT(NT1AM(ISYM),WORK(KT1),1,WORK(KT1),1)
            X2  = 0.0D0
            IF (.NOT.CCS) THEN
               X2  = DDOT(NT2AM(ISYM),WORK(KT2),1,WORK(KT2),1)
            ENDIF
            WRITE(LUPRI,'(1X,I5,5X,A8,2X,I4,6X,F15.10,F15.10,F15.10)')
     *           ILSTNR,LRZLBL(ILSTNR),ISYM,FRQLRZ(ILSTNR),X1,X2
C
            XLNCCCU = X1+X2
            XLNCCPR = XLRZ(ILSTNR)
            IF (ABS(XLNCCPR-XLNCCCU).GT.CVGLSOL) LSLLCVG = .FALSE.
            IF (IPRINT.GT.16)  WRITE(LUPRI,*) ' LSLLCVG = ',LSLLCVG
            XLRZ(ILSTNR) = XLNCCCU
C
 998     CONTINUE
C
C
         IF (LSLTCVG.AND.LSLLCVG) THEN
           WRITE(LUPRI,*)
     *      ' Solvent 1-order response equations converged in'
     *      ,ICCSLIT,' solvent iterations'
         ELSE
           ICCSLIT = ICCSLIT + 1
           IF (ICCSLIT.GT.MXCCSLIT) THEN
             CALL QUIT('Too many solvent iterations in CC_SOLRSP')
           ENDIF
           WRITE(LUPRI,*)' Solvent 1-order response equations not ',
     *               'converged in',ICCSLIT,' solvent iterations'
           WRITE(LUPRI,*)' A New iteration starts '
           GOTO 999
         ENDIF
C
      ENDIF
C======================================================================
C     Solve projected 1st-order Lagrange multiplier response equations.
C======================================================================
      LIST  = 'PL1'
      NSTAT = 0
      ORDER = 1
      ISIDE = -1

      IF ( .NOT. L1SKIP ) THEN
        CALL CC_SOLDRV(LIST,NSTAT,ORDER,ISIDE,APROXR12,
     &                 ISYSPL1, ISTPL1, EIGPL1, LPRPL1,
     &                 ISYPL1,  LBLPL1, FRQPL1, IDUMMY,
     &                 ISYOFPL1, NPL1LBL, MAXPL1LBL,
     &                 WORK, LWORK)

      ELSE
        WRITE (LUPRI,*) 
     &        'Skipped proj. 1-order left response (PL1) equations'
        L1SKIP = .FALSE.
      END IF

      IF (DEBUG .OR. LOCDBG. OR. IPRINT.GE.5) THEN
        CALL CC_TSTLST(LIST,NSTAT,ORDER,ISIDE,
     &                 ISYSPL1, ISTPL1, EIGPL1, LPRPL1,
     &                 ISYPL1,  LBLPL1, FRQPL1, IDUMMY,
     &                 ISYOFPL1, NPL1LBL, MAXPL1LBL, WORK, LWORK)
      END IF

C======================================================================
C     Calculate F-transformed eigenvectors
C======================================================================
      LIST = 'FR '
      NSTAT = 1
      ORDER = 0
      ISIDE = +1


      IF ((.NOT.M1SKIP).AND.(.NOT.FRSKIP)) THEN

         IF (NLRM .GT. 0) THEN

            IF (NLRM .GT. MAXFDIM) THEN
              WRITE (LUPRI,*)
     &              'ERROR: Dimension MAXFDIM to small in CC_SOLRSP.'
              CALL QUIT('Dimension MAXFDIM to small in CC_SOLRSP.')
            END IF
C
C           set up list of F matrix transformations:
C
            DO IDXM = 1, NLRM
              IDXRE = ILRM(IDXM)
              IFTRAN(1,IDXM) = 0
              IFTRAN(2,IDXM) = IDXRE
              IFTRAN(3,IDXM) = IDXM
            END DO
C
C           call CC_FMATRIX to do the F * RE transformations:
C
            IOPT = 3
            CALL CC_FMATRIX(IFTRAN, NLRM, 'L0 ', 'RE ', IOPT, 'FR ',
     &                      IDUMMY, DUMMY, 0, WORK, LWORK )


         END IF

      ELSE
        WRITE (LUPRI,*) 'Skipped of F x RE transformations'
        FRSKIP = .FALSE.
      END IF 

      IF (DEBUG .OR. LOCDBG. OR. IPRINT.GE.5) THEN
        CALL CC_TSTLST(LIST,NSTAT,ORDER,ISIDE,
     &                 ISYLRM, ILRM, FRQLRM, LDUMMY,
     &                 IDUMMY,   IDUMMY, DUMMY, IDUMMY,
     &                 ISYOFM, NLRM, MAXM, WORK, LWORK)
      END IF

C======================================================================
C     Solve equations for excited state M Lagrange multipliers:
C======================================================================
      LIST = 'M1 '
      NSTAT = 1
      ORDER = 0
      ISIDE = -1

      IF ( .NOT. M1SKIP ) THEN
        CALL CC_SOLDRV(LIST,NSTAT,ORDER,ISIDE,APROXR12,
     &                 ISYLRM, ILRM, FRQLRM, LDUMMY,
     &                 IDUMMY,   IDUMMY, DUMMY, IDUMMY,
     &                 ISYOFM, NLRM, MAXM,
     &                 WORK, LWORK)
      ELSE
        WRITE (LUPRI,*) 'Skipped first-order M vector (M1) equations'
        M1SKIP = .FALSE.
      END IF

      IF (DEBUG .OR. LOCDBG. OR. IPRINT.GE.5) THEN
        CALL CC_TSTLST(LIST,NSTAT,ORDER,ISIDE,
     &                 ISYLRM, ILRM, FRQLRM, LDUMMY,
     &                 IDUMMY,   IDUMMY, DUMMY, IDUMMY,
     &                 ISYOFM, NLRM, MAXM, WORK, LWORK)
      END IF

C======================================================================
C     Make rhs for excited state Lagrangian multipliers E0:
C======================================================================

      LIST = 'BE '
      NSTAT = 0
      ORDER = 0
      ISIDE = -1

      IF (.NOT.BESKIP) THEN

         IF (NXGRST .GT. 0) THEN

            IF (NXGRST .GT. MAXFDIM) 
     &        CALL QUIT('Dimension MAXFDIM to small in CC_SOLRSP.')

            ! set up list of 'F matrix' transformations:
            DO IDXF = 1, NXGRST
              IFTRAN(1,IDXF) = IXGRST(IDXF)
              IFTRAN(2,IDXF) = IXGRST(IDXF)
              IFTRAN(3,IDXF) = IDXF
            END DO

            ! call CC_FMATRIX to do the LE * B * RE transformations:
            IOPT = 3
            CALL CC_FMATRIX(IFTRAN, NXGRST, 'LE ', 'RE ', IOPT, LIST,
     &                      IDUMMY, DUMMY, 0, WORK, LWORK )

            ! add eta^(0) to complete rhs vectors:
            NTAMP = NT1AMX
            IF (.NOT.CCS) NTAMP = NTAMP + NT2AMX
            IF (CCR12) NTAMP = NTAMP + NTR12AM(1)
            KETA  = 1
            KEND1 = KETA  + NTAMP
            LWRK1 = LWORK - KEND1

            KBEVC = KEND1
            KEND2 = KBEVC + NTAMP
            LWRK2 = LWORK - KEND2
            IF (LWRK2.LT.0) CALL QUIT('Out of memory in CC_SOLRSP.')

            CALL CC_ETA(WORK(KETA),WORK(KEND1),LWRK1)

            DO IDXF = 1, NXGRST
              IOPT = 3
              CALL CC_RDRSP('BE ',IDXF,ISYM0,IOPT,MODEL,
     &                       WORK(KBEVC),WORK(KBEVC+NT1AMX))
              CALL DAXPY(NTAMP,ONE,WORK(KETA),1,WORK(KBEVC),1) 
              IOPT = 3
              CALL CC_WRRSP('BE ',IDXF,ISYM0,IOPT,MODEL,DUMMY,
     &               WORK(KBEVC),WORK(KBEVC+NT1AMX),WORK(KEND2),LWRK2)
              IF (CC3) THEN
                ! for CC3 also correct the effective RHS vectors
                IOPT = 24
                CALL CC_RDRSP('BE ',IDXF,ISYM0,IOPT,MODEL,
     &                        WORK(KBEVC),WORK(KBEVC+NT1AMX))
                CALL DAXPY(NTAMP,ONE,WORK(KETA),1,WORK(KBEVC),1) 
                IOPT = 24
                CALL CC_WRRSP('BE ',IDXF,ISYM0,IOPT,MODEL,DUMMY,
     &                WORK(KBEVC),WORK(KBEVC+NT1AMX),WORK(KEND2),LWRK2)
              END IF
            END DO
         END IF

      ELSE
        WRITE (LUPRI,*) 'Skipped of LE x B x RE transformations'
        BESKIP = .FALSE.
      END IF

      IF (DEBUG .OR. LOCDBG. OR. IPRINT.GE.5) THEN
        NSTAT = 1
        CALL CC_TSTLST(LIST,NSTAT,ORDER,ISIDE,
     &                 IDUMMY,  IXGRST, EIGVAL, LDUMMY,
     &                 IDUMMY,  IDUMMY, DUMMY, IDUMMY,
     &                 ISYOFXG, NXGRST, MXXGST, WORK, LWORK)
      END IF

C======================================================================
C     Solve equations for excited state lagrange multipliers E0:
C======================================================================

      LIST = 'E0 '
      NSTAT = 0
      ORDER = 0
      ISIDE = -1

      IF ( .NOT. E0SKIP ) THEN
        CALL CC_SOLDRV(LIST,NSTAT,ORDER,ISIDE,APROXR12,
     &                 IDUMMY,  IDUMMY, DUMMY,  LDUMMY,
     &                 IDUMMY,  IDUMMY, DUMMY,  IDUMMY,
     &                 ISYOFXG, NXGRST, MXXGST,
     &                 WORK, LWORK)
      ELSE
        WRITE (LUPRI,*) 'Skipped exc. st. Lagr. mult.  (E0) equations'
        E0SKIP = .FALSE.
      END IF

      IF (DEBUG .OR. LOCDBG. OR. IPRINT.GE.5) THEN
        NSTAT = 1
        CALL CC_TSTLST(LIST,NSTAT,ORDER,ISIDE,
     &                 IDUMMY,  IXGRST, EIGVAL, LDUMMY,
     &                 IDUMMY,  IDUMMY, DUMMY,  IDUMMY,
     &                 ISYOFXG, NXGRST, MXXGST, WORK, LWORK)
      END IF

C======================================================================
C     Make rhs for N2-vector equations.
C======================================================================

      LIST = 'BR '
      NSTAT = 2
      ORDER = 0
      ISIDE = -1

      IF ((.NOT.N2SKIP).AND.(.NOT.BRSKIP)) THEN

         IF (NQRN2 .GT. 0) THEN

            IF (NQRN2 .GT. MAXFDIM) THEN
              WRITE (LUPRI,*)
     &              'ERROR: Dimension MAXFDIM to small in CC_SOLRSP.'
              CALL QUIT('Dimension MAXFDIM to small in CC_SOLRSP.')
            END IF
C
C           set up list of 'F matrix' transformations:
C
            DO IDXF = 1, NQRN2
              IFTRAN(1,IDXF) = IIN2(IDXF)
              IFTRAN(2,IDXF) = IFN2(IDXF)
              IFTRAN(3,IDXF) = IDXF
            END DO
C
C           call CC_FMATRIX to do the LE * B * RE transformations:
C
            IOPT = 3
            CALL CC_FMATRIX(IFTRAN, NQRN2, 'LE ', 'RE ', IOPT, LIST,
     &                      IDUMMY, DUMMY, 0, WORK, LWORK )

         END IF

      ELSE
        WRITE (LUPRI,*) 'Skipped of LE x B x RE transformations'
        BRSKIP = .FALSE.
      END IF

      IF (DEBUG .OR. LOCDBG. OR. IPRINT.GE.5) THEN
        CALL CC_TSTLST(LIST,NSTAT,ORDER,ISIDE,
     &                 ISYSN2, ISTN2, EIGN2, LDUMMY,
     &                 IDUMMY,   IDUMMY, DUMMY, IDUMMY,
     &                 ISYOFN2, NQRN2, MAXQRN2, WORK, LWORK)
      END IF

C======================================================================
C     Solve N2-vector for quadratic response 2. residue calculation.
C======================================================================

      LIST = 'N2 '
      NSTAT = 2
      ORDER = 0
      ISIDE = -1

      IF ( .NOT. N2SKIP ) THEN
        CALL CC_SOLDRV(LIST,NSTAT,ORDER,ISIDE,APROXR12,
     &                 ISYSN2, ISTN2, EIGN2, LDUMMY,
     &                 IDUMMY,   IDUMMY, DUMMY, IDUMMY,
     &                 ISYOFN2, NQRN2, MAXQRN2,
     &                 WORK, LWORK)
      ELSE
        WRITE (LUPRI,*) 'Skipped second-order N vector (N2) equations'
        N2SKIP = .FALSE.
      END IF

      IF (DEBUG .OR. LOCDBG. OR. IPRINT.GE.5) THEN
        CALL CC_TSTLST(LIST,NSTAT,ORDER,ISIDE,
     &                 ISYSN2, ISTN2, EIGN2, LDUMMY,
     &                 IDUMMY,   IDUMMY, DUMMY, IDUMMY,
     &                 ISYOFN2, NQRN2, MAXQRN2, WORK, LWORK)
      END IF

C======================================================================
C     Solve first-order right Cauchy equations.
C======================================================================
      LIST  = 'RC '
      NSTAT = 0   
      ORDER = 1   
      ISIDE = +1  

      IF ( .NOT. RCSKIP ) THEN

       IF (NEWCAU) THEN

        ! make sure that RC files are initialized:
        ! if no restart vector available on file, 
        ! use next lower order RC vect.
        DO IDXRC = 1, NLRCLBL
          ISYM  = ISYLRC(IDXRC)
          KVEC1 = 1
          KVEC2 = KVEC1 + NT1AM(ISYM)
          KEND1 = KVEC2 + NT2AM(ISYM)
          LEND1 = LWORK - KEND1
          IF (LEND1 .LT. 0) THEN
            CALL QUIT('INSUFFICIENT WORKSPACE IN CC_SOLRSP.')
          END IF
          IOPT = 33
          CALL CC_RDRSP('RC',IDXRC,ISYM,IOPT,MODEL,
     &                  WORK(KVEC1),WORK(KVEC2)  )

          IF (IOPT.EQ.33) THEN

            NCAU = ILRCAU(IDXRC)
            IF (NCAU.EQ.1) THEN 
              IDXR1 = IR1TAMP(LRCLBL(IDXRC),.FALSE.,0.0d0,ISYM)
              CALL CC_RDRSP('R1',IDXR1,ISYM,IOPT,MODEL,
     &                      WORK(KVEC1),WORK(KVEC2)  )
            ELSE
              IDXR1 = ILRCAMP(LRCLBL(IDXRC),NCAU-1,ISYM)
              CALL CC_RDRSP('RC',IDXR1,ISYM,IOPT,MODEL,
     &                      WORK(KVEC1),WORK(KVEC2)  )
            END IF

            CALL CC_WRRSP('RC',IDXRC,ISYM,IOPT,MODEL,DUMMY,
     &                    WORK(KVEC1),WORK(KVEC2),WORK(KEND1),LEND1)
          END IF

        END DO
       END IF

       ! Put zeroth-order Cauchy moments at the end of the list,
       ! but without increasing the number of equations that will
       ! solved. The read routine CC_RDRSP will matp these vectors
       ! onto the first-order amplitude response vectors which are
       ! available from the R1 section above.
       NLRCLBLSAVE = NLRCLBL
       LRC1OPN = .TRUE. ! open list
       DO IVEC = 1, NLRCLBLSAVE
         INUM = ILRCAMP(LRCLBL(IVEC),0,ISYLRC(IVEC))
       END DO
       LRC1OPN = .FALSE. ! close list


       CALL CC_SOLDRV(LIST,NSTAT,ORDER,ISIDE,APROXR12,
     &                IDUMMY, IDUMMY, DUMMY, LDUMMY,
     &                ISYLRC, LRCLBL, DUMMY, ILRCAU,
     &                ISYOFC, NLRCLBLSAVE, MAXCLBL,
     &                WORK, LWORK)

      ELSE
        WRITE (LUPRI,*) 'Skipped first-order right Cauchy equations'
        RCSKIP = .FALSE.
      END IF

C======================================================================
C     Add (unrelaxed) static R1 vectors to Cauchy RC list as RC(0)
C======================================================================
      LIST  = 'RC '
      NSTAT = 0   
      ORDER = 1   
      ISIDE = +1  

      LRC1OPN = .TRUE.
      DO IVEC = 1, NLRTLBL
        IF (ABS(FRQLRT(IVEC)).LT.TOL .AND. (.NOT.LORXLRT(IVEC)) )  THEN
          INUM = ILRCAMP(LRTLBL(IVEC),0,ISYLRT(IVEC))
        END IF
      END DO
      LRC1OPN = .FALSE.

      IF (DEBUG .OR. LOCDBG. OR. IPRINT.GE.5) THEN
       CALL CC_TSTLST(LIST,NSTAT,ORDER,ISIDE,
     &                IDUMMY, IDUMMY, DUMMY, LDUMMY,
     &                ISYLRC, LRCLBL, DUMMY, ILRCAU,
     &                ISYOFC, NLRCLBL, MAXCLBL, WORK, LWORK)
      END IF

C======================================================================
C     precalculate F * RC transformations required to solve for
C     the first-order left Cauchy vectors LC
C======================================================================
      LIST  = 'FC '
      NSTAT = 0   
      ORDER = 1   
      ISIDE = +1  

      IF ( .NOT. FCSKIP ) THEN

         IF (NLC1LBL .GT. 0) THEN
C
C          check compatibility with dimension of IFTRAN:
C
           IF (NLC1LBL .GT. MAXFDIM) THEN
             WRITE (LUPRI,*)
     &             'ERROR: Dimension MAXFDIM to small in CC_SOLRSP.'
             CALL QUIT('Dimension MAXFDIM to small in CC_SOLRSP.')
           END IF
C
C          set up list of F matrix transformations:
C
           DO IDXLC = 1, NLC1LBL
             IDXRC=ILRCAMP(LBLLC1(IDXLC),ILC1CAU(IDXLC),ISYLC1(IDXLC))
             IFTRAN(1,IDXLC) = 0
             IFTRAN(2,IDXLC) = IDXRC
             IFTRAN(3,IDXLC) = IDXRC
           END DO
C
C          call CC_FMATRIX to do the F * RC transformations:
C
           IOPT = 3
           CALL CC_FMATRIX(IFTRAN, NLC1LBL, 'L0 ', 'RC ', IOPT, 'FC ',
     &                     IDUMMY, DUMMY, 0, WORK, LWORK )

         ENDIF

      ELSE
        WRITE (LUPRI,*) 'Skipped calculation of F x RC transformations'
        FCSKIP = .FALSE.
      END IF

C==================================================================
C     Solve first-order left Cauchy equations:
C==================================================================
      LIST  = 'LC '
      NSTAT = 0   
      ORDER = 1   
      ISIDE = -1  

      IF ( .NOT. LCSKIP ) THEN
        CALL CC_SOLDRV(LIST,NSTAT,ORDER,ISIDE,APROXR12,
     &                 IDUMMY, IDUMMY, DUMMY, LDUMMY,
     &                 ISYLC1, LBLLC1, DUMMY, ILC1CAU,
     &                 ISYOFLC1, NLC1LBL, MAXLC1LBL,
     &                 WORK, LWORK)
      ELSE
        WRITE (LUPRI,*) 'Skipped first-order left Cauchy equations'
        LCSKIP = .FALSE.
      END IF

C======================================================================
C     Add static L1 vectors to Cauchy LC list as LC(0)
C======================================================================
      LIST  = 'LC '
      NSTAT = 0   
      ORDER = 1   
      ISIDE = -1  

      LLC1OPN = .TRUE.
      DO IVEC = 1, NLRZLBL
        IF (ABS(FRQLRZ(IVEC)).LT.TOL .AND. (.NOT.LORXLRZ(IVEC)) )  THEN
          INUM = ILC1AMP(LRZLBL(IVEC),0,ISYLRZ(IVEC))
        END IF
      END DO
      LLC1OPN = .FALSE.

      IF (DEBUG .OR. LOCDBG. OR. IPRINT.GE.5) THEN
        CALL CC_TSTLST(LIST,NSTAT,ORDER,ISIDE,
     &                 IDUMMY, IDUMMY, DUMMY, LDUMMY,
     &                 ISYLC1, LBLLC1, DUMMY, ILC1CAU,
     &                 ISYOFLC1, NLC1LBL, MAXLC1LBL, WORK, LWORK)
      END IF

C======================================================================
C     precalculate right hand side vectors for second-order
C     cluster amplitude response equations.
C======================================================================
      LIST  = 'O2 '
      NSTAT = 0   
      ORDER = 2   
      ISIDE = +1  

      IF ( .NOT. O2SKIP ) THEN

      CALL CCRHSVEC('O2 ',LBLO2,IDUMMY,IDUMMY,DUMMY,ISYO2,FRQO2,
     &              LORXO2,IDUMMY,NO2LBL,MAXO2LBL,0,
     &              WORK,LWORK)

      ELSE
        WRITE (LUPRI,*)
     &        'Skipped calculation of O2 right hand side vectors.'
        O2SKIP = .FALSE.
      END IF

      IF (DEBUG .OR. LOCDBG. OR. IPRINT.GE.5) THEN
        CALL CC_TSTLST(LIST,NSTAT,ORDER,ISIDE,
     &                 IDUMMY, IDUMMY, DUMMY, LDUMMY,
     &                 ISYO2, LBLO2, FRQO2, IDUMMY,
     &                 ISYOFO2, NO2LBL, MAXO2LBL, WORK, LWORK)
      END IF

C======================================================================
C     Solve second-order right cluster amplitude response equations.
C======================================================================
      LIST  = 'R2 '
      NSTAT = 0   
      ORDER = 2   
      ISIDE = +1  

      IF ( .NOT. R2SKIP ) THEN
        CALL CC_SOLDRV(LIST,NSTAT,ORDER,ISIDE,APROXR12,
     &                 IDUMMY, IDUMMY, DUMMY, LDUMMY,
     &                 ISYR2T, LBLR2T, FRQR2T, IDUMMY,
     &                 ISYOFT2, NR2TLBL, MAXT2LBL,
     &                 WORK, LWORK)
      ELSE
        WRITE (LUPRI,*) 'Skipped second-order response (R2) equations'
        R2SKIP = .FALSE.
      END IF

      IF (DEBUG .OR. LOCDBG. OR. IPRINT.GE.5) THEN
        CALL CC_TSTLST(LIST,NSTAT,ORDER,ISIDE,
     &                 IDUMMY, IDUMMY, DUMMY, LDUMMY,
     &                 ISYR2T, LBLR2T, FRQR2T, IDUMMY,
     &                 ISYOFT2, NR2TLBL, MAXT2LBL, WORK, LWORK)
      END IF

C======================================================================
C     precalculate CO2 vectors which contain all contributions to the
C     rhs vectors for the second-order right Cauchy vectors which are 
C     independent of second-order right Cauchy vectors
C======================================================================
      LIST  = 'CO2'
      NSTAT = 0
      ORDER = 2
      ISIDE = +1

      IF ( .NOT. CO2SKIP ) THEN

        CALL CCRHSVEC('CO2',LBLCO2,IDUMMY,IDUMMY,DUMMY,ISYCO2,DUMMY,
     &                LDUMMY,ICO2CAU,NCO2LBL,MAXCO2LBL,0,
     &                WORK,LWORK)

      ELSE
        WRITE (LUPRI,*)
     &        'Skipped calculation of CO2 right hand side vectors.'
        CO2SKIP = .FALSE.
      END IF

C===================================================================
C     Solve second-order right Cauchy response equations.
C===================================================================
      LIST  = 'CR2'
      NSTAT = 0
      ORDER = 2
      ISIDE = +1

      IF ( .NOT. CR2SKIP ) THEN
        CALL CC_SOLDRV(LIST,NSTAT,ORDER,ISIDE,APROXR12,
     &                 IDUMMY, IDUMMY, DUMMY, LDUMMY,
     &                 ISYCR2, LBLCR2, DUMMY, ICR2CAU,
     &                 ISYOFCR2, NCR2LBL, MAXCR2LBL,
     &                 WORK, LWORK)
      ELSE
        WRITE (LUPRI,*)
     &        'Skipped second-order right Cauchy (CR2) equations'
        CR2SKIP = .FALSE.
      END IF

C======================================================================
C     Add static R2 vectors to second-order Cauchy CR2 list as CR2(0,0)
C======================================================================
      LIST  = 'CR2'
      NSTAT = 0
      ORDER = 2
      ISIDE = +1

      LCR2OPN = .TRUE.
      DO IVEC = 1, NR2TLBL
        IF (ABS(FRQR2T(IVEC,1)).LT.TOL .AND. ABS(FRQR2T(IVEC,2)).LT.TOL
     &    .AND. (.NOT.LORXR2T(IVEC,1)) .AND. (.NOT.LORXR2T(IVEC,2))
     &     )  THEN
          INUM = ICR2AMP(LBLR2T(IVEC,1),0,ISYR2T(IVEC,1),
     &                   LBLR2T(IVEC,2),0,ISYR2T(IVEC,2))
        END IF
      END DO
      LCR2OPN = .FALSE.

      IF (DEBUG .OR. LOCDBG. OR. IPRINT.GE.5) THEN
        CALL CC_TSTLST(LIST,NSTAT,ORDER,ISIDE,
     &                 IDUMMY, IDUMMY, DUMMY, LDUMMY,
     &                 ISYCR2, LBLCR2, DUMMY, ICR2CAU,
     &                 ISYOFCR2, NCR2LBL, MAXCR2LBL, WORK, LWORK)
      END IF

C======================================================================
C     Add static O2 vectors to second-order Cauchy CO2 list as CO2(0,0)
C======================================================================
      LIST  = 'CO2'
      NSTAT = 0
      ORDER = 2
      ISIDE = +1

      LCO2OPN = .TRUE.
      DO IVEC = 1, NO2LBL
        IF (ABS(FRQO2(IVEC,1)).LT.TOL .AND. ABS(FRQO2(IVEC,2)).LT.TOL
     &     )  THEN
          INUM = IRHSCR2(LBLO2(IVEC,1),0,ISYO2(IVEC,1),
     &                   LBLO2(IVEC,2),0,ISYO2(IVEC,2))
        END IF
      END DO
      LCO2OPN = .FALSE.

      IF (DEBUG .OR. LOCDBG. OR. IPRINT.GE.5) THEN
        CALL CC_TSTLST(LIST,NSTAT,ORDER,ISIDE,
     &                 IDUMMY, IDUMMY, DUMMY, LDUMMY,
     &                 ISYCO2, LBLCO2, DUMMY, ICO2CAU,
     &                 ISYOFCO2, NCO2LBL, MAXCO2LBL, WORK, LWORK)
      END IF

C======================================================================
C     precalculate second-order chi vectors required to set up the
C     right hande side vector for the second-order multiplier equations
C======================================================================
      IF ( .NOT. X2SKIP ) THEN

        CALL CC_ETADRV('X2 ',LBLX2,IDUMMY,IDUMMY,DUMMY,
     &                 ISYX2,FRQX2,IDUMMY,NX2LBL,MAXX2LBL,
     &                 WORK,LWORK)

      ELSE
        WRITE (LUPRI,*)
     &        'Skipped calculation of X2 right hand side vectors.'
        X2SKIP = .FALSE.
      END IF

C======================================================================
C     precalculate F * R2 transformations required to solve for
C     the second order coupled cluster multiplier reponses L2
C======================================================================
      IF ( .NOT. F2SKIP ) THEN

         IF (NL2LBL .GT. 0) THEN
C
C          check compatibility with dimension of IFTRAN:
C
           IF (NL2LBL .GT. MAXFDIM) THEN
             WRITE (LUPRI,*)
     &             'ERROR: Dimension MAXFDIM to small in CC_SOLRSP.'
             CALL QUIT('Dimension MAXFDIM to small in CC_SOLRSP.')
           END IF
C
C          set up list of F matrix transformations:
C
           DO IDXL2 = 1, NL2LBL
            IDXR2 = 
     &        IR2TAMP(LBLAL2(IDXL2),.FALSE.,FRQAL2(IDXL2),ISYAL2(IDXL2),
     &                LBLBL2(IDXL2),.FALSE.,FRQBL2(IDXL2),ISYBL2(IDXL2))

            IFTRAN(1,IDXL2) = 0
            IFTRAN(2,IDXL2) = IDXR2
            IFTRAN(3,IDXL2) = IDXR2
           END DO
C
C          call CC_FMATRIX to do the F * R2 transformations:
C
           IOPT = 3
           CALL CC_FMATRIX(IFTRAN, NL2LBL, 'L0 ', 'R2 ', IOPT, 'F2 ',
     &                     IDUMMY, DUMMY, 0, WORK, LWORK )

         ENDIF

      ELSE
        WRITE (LUPRI,*) 'Skipped calculation of F x R2 transformations'
        F2SKIP = .FALSE.
      END IF

C==================================================================
C     Solve second-order lagrange multiplier response equations.
C==================================================================
      LIST  = 'L2 '
      NSTAT = 0   
      ORDER = 2   
      ISIDE = -1  

      IF ( .NOT. L2SKIP ) THEN
        CALL CC_SOLDRV(LIST,NSTAT,ORDER,ISIDE,APROXR12,
     &                 IDUMMY, IDUMMY, DUMMY, LDUMMY,
     &                 ISYL2, LBLL2, FRQL2, IDUMMY,
     &                 ISYOFL2, NL2LBL, MAXL2LBL,
     &                 WORK, LWORK)
      ELSE
        WRITE (LUPRI,*)
     &        'Skipped second-order left response (L2) equations'
        L2SKIP = .FALSE.
      END IF

      IF (DEBUG .OR. LOCDBG. OR. IPRINT.GE.5) THEN
        CALL CC_TSTLST(LIST,NSTAT,ORDER,ISIDE,
     &                 IDUMMY, IDUMMY, DUMMY, LDUMMY,
     &                 ISYL2, LBLL2, FRQL2, IDUMMY,
     &                 ISYOFL2, NL2LBL, MAXL2LBL, WORK, LWORK)
      END IF

C======================================================================
C     precalculate second-order Cauchy eta vectors required to set up 
C     the right hande side vectors for the second-order left Cauchy
C     equations
C======================================================================
      IF ( .NOT. CX2SKIP ) THEN

      CALL CC_ETADRV('CX2',LBLCX2,IDUMMY,IDUMMY,DUMMY,
     &               ISYCX2,DUMMY,ICX2CAU,NCX2LBL,MAXCX2LBL,
     &               WORK,LWORK)

      ELSE
        WRITE (LUPRI,*)
     &        'Skipped calculation of CX2 right hand side vectors.'
        CX2SKIP = .FALSE.
      END IF
C======================================================================
C     precalculate F * CR2 transformations required to solve for
C     the second order left Cauchy vectors CL2
C======================================================================
      IF ( .NOT. CF2SKIP ) THEN

         IF (NCL2LBL .GT. 0) THEN
C
C          check compatibility with dimension of IFTRAN:
C
           IF (NCL2LBL .GT. MAXFDIM) THEN
             WRITE (LUPRI,*)
     &             'ERROR: Dimension MAXFDIM to small in CC_SOLRSP.'
             CALL QUIT('Dimension MAXFDIM to small in CC_SOLRSP.')
           END IF
C
C          set up list of F matrix transformations:
C
           DO IDXL2 = 1, NCL2LBL
             IDXR2 = ICR2AMP(
     &                LBLCL2(IDXL2,1),ICL2CAU(IDXL2,1),ISYCL2(IDXL2,1),
     &                LBLCL2(IDXL2,2),ICL2CAU(IDXL2,2),ISYCL2(IDXL2,2))

             IFTRAN(1,IDXL2) = 0
             IFTRAN(2,IDXL2) = IDXR2
             IFTRAN(3,IDXL2) = IDXR2
           END DO
C
C          call CC_FMATRIX to do the F * CR2 transformations:
C
           IOPT = 3
           CALL CC_FMATRIX(IFTRAN, NCL2LBL, 'L0 ', 'CR2', IOPT, 'CF2',
     &                     IDUMMY, DUMMY, 0, WORK, LWORK )

         ENDIF

      ELSE
        WRITE (LUPRI,*) 'Skipped calculation of F x CR2 transformations'
        CF2SKIP = .FALSE.
      END IF

C==================================================================
C     Solve second-order left Cauchy vector ("CL2") equations.
C==================================================================
      LIST  = 'CL2'
      NSTAT = 0   
      ORDER = 2   
      ISIDE = -1  

      IF ( .NOT. CL2SKIP ) THEN
        CALL CC_SOLDRV(LIST,NSTAT,ORDER,ISIDE,APROXR12,
     &                 IDUMMY, IDUMMY, DUMMY, LDUMMY,
     &                 ISYCL2, LBLCL2, DUMMY, ICL2CAU,
     &                 ISYOFCL2, NCL2LBL, MAXCL2LBL,
     &                 WORK, LWORK)
      ELSE
        WRITE (LUPRI,*)
     &        'Skipped second-order left Cauchy (CL2) equations'
        CL2SKIP = .FALSE.
      END IF

C======================================================================
C     Add static L2 vectors to second-order Cauchy CL2 list as CL2(0,0)
C======================================================================
      LIST  = 'CL2'
      NSTAT = 0   
      ORDER = 2   
      ISIDE = -1  

      LCL2OPN = .TRUE.
      DO IVEC = 1, NL2LBL
        IF (ABS(FRQL2(IVEC,1)).LT.TOL .AND. ABS(FRQL2(IVEC,2)).LT.TOL
     &    .AND. (.NOT.LORXL2(IVEC,1)) .AND. (.NOT.LORXL2(IVEC,2))
     &     )  THEN
          INUM = ICL2AMP(LBLL2(IVEC,1),0,ISYL2(IVEC,1),
     &                   LBLL2(IVEC,2),0,ISYL2(IVEC,2))
        END IF
      END DO
      LCL2OPN = .FALSE.

      IF (DEBUG .OR. LOCDBG. OR. IPRINT.GE.5) THEN
        CALL CC_TSTLST(LIST,NSTAT,ORDER,ISIDE,
     &                 IDUMMY, IDUMMY, DUMMY, LDUMMY,
     &                 ISYCL2, LBLCL2, DUMMY, ICL2CAU,
     &                 ISYOFCL2, NCL2LBL, MAXCL2LBL, WORK, LWORK)
      END IF

C======================================================================
C     Add static X2 vectors to second-order Cauchy CX2 list as CX2(0,0)
C======================================================================
      LCX2OPN = .TRUE.
      DO IVEC = 1, NX2LBL
        IF (ABS(FRQX2(IVEC,1)).LT.TOL .AND. ABS(FRQX2(IVEC,2)).LT.TOL
     &     )  THEN
          INUM = IETACL2(LBLX2(IVEC,1),0,ISYX2(IVEC,1),
     &                   LBLX2(IVEC,2),0,ISYX2(IVEC,2))
        END IF
      END DO
      LCX2OPN = .FALSE.

C===================================================================
C     precalculate right hand side vectors for third-order
C     cluster amplitude response equations.
C===================================================================
      CALL CCRHSVEC('O3 ',LBLO3,IDUMMY,IDUMMY,DUMMY,ISYO3,FRQO3,
     &              LORXO3,IDUMMY,NO3LBL,MAXO3LBL,0,
     &              WORK,LWORK)


C===================================================================
C     Solve third-order right cluster amplitude response equations.
C===================================================================
      LIST  = 'R3 '
      NSTAT = 0   
      ORDER = 3   
      ISIDE = +1  

      CALL CC_SOLDRV(LIST,NSTAT,ORDER,ISIDE,APROXR12,
     &               IDUMMY, IDUMMY, DUMMY, LDUMMY,
     &               ISYR3T, LBLR3T, FRQR3T, IDUMMY,
     &               ISYOFT3, NR3TLBL, MAXT3LBL,
     &               WORK, LWORK)


C======================================================================
C     precalculate third-order chi vectors required to set up the
C     right hande side vector for the third-order multiplier equations
C======================================================================

      CALL CC_ETADRV('X3 ',LBLX3,IDUMMY,IDUMMY,DUMMY,
     &               ISYX3,FRQX3,IDUMMY,NX3LBL,MAXX3LBL,
     &               WORK,LWORK)

C======================================================================
C     precalculate F * R3 transformations required to solve for
C     the second order coupled cluster multiplier reponses L3
C======================================================================
      IF (NL3LBL .GT. 0) THEN
C
C       check compatibility with dimension of IFTRAN:
C
        IF (NL3LBL .GT. MAXFDIM) THEN
          WRITE (LUPRI,*)
     &          'ERROR: Dimension MAXFDIM to small in CC_SOLRSP.'
          CALL QUIT('Dimension MAXFDIM to small in CC_SOLRSP.')
        END IF
C
C       set up list of F matrix transformations:
C
        DO IDXL3 = 1, NL3LBL
          IDXR3 = IR3TAMP(LBLL3(IDXL3,1),FRQL3(IDXL3,1),ISYL3(IDXL3,1),
     &                    LBLL3(IDXL3,2),FRQL3(IDXL3,2),ISYL3(IDXL3,2),
     &                    LBLL3(IDXL3,3),FRQL3(IDXL3,3),ISYL3(IDXL3,3))

          IFTRAN(1,IDXL3) = 0
          IFTRAN(2,IDXL3) = IDXR3
          IFTRAN(3,IDXL3) = IDXR3
        END DO
C
C       call CC_FMATRIX to do the F * R3 transformations:
C
        IOPT = 3
        CALL CC_FMATRIX(IFTRAN, NL3LBL, 'L0 ', 'R3 ', IOPT, 'F3 ',
     &                  IDUMMY, DUMMY, 0, WORK, LWORK )

      ENDIF

C==================================================================
C     Solve second-order lagrange multiplier response equations.
C==================================================================
      LIST  = 'L3 '
      NSTAT = 0   
      ORDER = 3   
      ISIDE = -1  

      CALL CC_SOLDRV(LIST,NSTAT,ORDER,ISIDE,APROXR12,
     &               IDUMMY, IDUMMY, DUMMY, LDUMMY,
     &               ISYL3, LBLL3, FRQL3, IDUMMY,
     &               ISYOFL3, NL3LBL, MAXL3LBL,
     &               WORK, LWORK)


C===================================================================
C     precalculate right hand side vectors for first-order
C     right eigenvector response equations.
C===================================================================
      CALL CCRHSVEC('EO1',LBLER1,ISYSER1,ISTER1,EIGER1,ISYOER1,FRQER1,
     &              LORXER1,IDUMMY,NER1LBL,MAXER1LBL,0,
     &              WORK,LWORK)


C===================================================================
C     Solve first-order right eigenvector response equations.
C===================================================================
      LIST  = 'ER1'
      NSTAT = 1   
      ORDER = 1   
      ISIDE = +1  

      CALL CC_SOLDRV(LIST,NSTAT,ORDER,ISIDE,APROXR12,
     &               ISYSER1, ISTER1, EIGER1, LPRER1,
     &               ISYOER1, LBLER1, FRQER1, IDUMMY,
     &               ISYOFER1, NER1LBL, MAXER1LBL,
     &               WORK, LWORK)

      IF (DEBUG .OR. LOCDBG. OR. IPRINT.GE.5) THEN
        CALL CC_TSTLST(LIST,NSTAT,ORDER,ISIDE,
     &                 ISYSER1, ISTER1, EIGER1, LPRER1,
     &                 ISYOER1, LBLER1, FRQER1, IDUMMY,
     &                 ISYOFER1, NER1LBL, MAXER1LBL, WORK, LWORK)
      END IF

C===================================================================
C     precalculate right hand side vectors for second-order
C     right eigenvector response equations.
C===================================================================
      CALL CCRHSVEC('EO2',LBLER2,ISYSER2,ISTER2,EIGER2,ISYOER2,FRQER2,
     &              LORXER2,IDUMMY,NER2LBL,MAXER2LBL,0,
     &              WORK,LWORK)


C===================================================================
C     Solve second-order right eigenvector response equations.
C===================================================================
      LIST  = 'ER2'
      NSTAT = 1   
      ORDER = 2   
      ISIDE = +1  

      CALL CC_SOLDRV(LIST,NSTAT,ORDER,ISIDE,APROXR12,
     &               ISYSER2, ISTER2, EIGER2, LPRER2,
     &               ISYOER2, LBLER2, FRQER2, IDUMMY,
     &               ISYOFER2, NER2LBL, MAXER2LBL,
     &               WORK, LWORK)

C===================================================================
C     precalculate right hand side vectors (EX1) for first-order
C     left eigenvector response (EL1) equations.
C===================================================================
      if (lold) then
      WRITE(LUPRI,*) 'CC_SOLRSP: use old CC_ETADRV for EX1...'
      CALL CC_ETADRV('EX1',LBLEL1,ISYSEL1,ISTEL1,EIGEL1,
     &             ISYOEL1,FRQEL1,IDUMMY,NEL1LBL,MAXEL1LBL,
     &             WORK,LWORK)
      else
       CALL CC_ETADRV1('EX1',LBLEL1,ISYSEL1,ISTEL1,EIGEL1,
     &               ISYOEL1,FRQEL1,LORXEL1,IDUMMY,NEL1LBL,MAXEL1LBL,
     &               WORK,LWORK)
      end if

C===================================================================
C     Solve first-order left eigenvector response equations.
C===================================================================
      LIST  = 'EL1'
      NSTAT = 1   
      ORDER = 1   
      ISIDE = -1  

      CALL CC_SOLDRV(LIST,NSTAT,ORDER,ISIDE,APROXR12,
     &               ISYSEL1, ISTEL1, EIGEL1, LPREL1,
     &               ISYOEL1, LBLEL1, FRQEL1, IDUMMY,
     &               ISYOFEL1, NEL1LBL, MAXEL1LBL,
     &               WORK, LWORK)


      IF (DEBUG .OR. LOCDBG. OR. IPRINT.GE.5) THEN
        CALL CC_TSTLST(LIST,NSTAT,ORDER,ISIDE,
     &                 ISYSEL1, ISTEL1, EIGEL1, LPREL1,
     &                 ISYOEL1, LBLEL1, FRQEL1, IDUMMY,
     &                 ISYOFEL1, NEL1LBL, MAXEL1LBL, WORK, LWORK)
      END IF

C===================================================================
C     precalculate right hand side vectors for second-order
C     left eigenvector response equations.
C===================================================================
      CALL CC_ETADRV('EX2',LBLEL2,ISYSEL2,ISTEL2,EIGEL2,
     &               ISYOEL2,FRQEL2,IDUMMY,NEL2LBL,MAXEL2LBL,
     &               WORK,LWORK)


C===================================================================
C     Solve first-order left eigenvector response equations.
C===================================================================
      LIST  = 'EL2'
      NSTAT = 1   
      ORDER = 2   
      ISIDE = -1  

      CALL CC_SOLDRV(LIST,NSTAT,ORDER,ISIDE,APROXR12,
     &               ISYSEL2, ISTEL2, EIGEL2, LPRER2,
     &               ISYOEL2, LBLEL2, FRQEL2, IDUMMY,
     &               ISYOFEL2, NEL2LBL, MAXEL2LBL,
     &               WORK, LWORK)

C---------------------------------------------------------------------
C     TEST - Loop over list of response vectors 
C---------------------------------------------------------------------
C
C     Test average values for freq = 0 and tot sym vectors.
C 
      IF (DEBUG .OR. (IPRINT.GE.5)) THEN
         DO ILSTNR = 1, NLRTLBL
C
            ISYM = ISYLRT(ILSTNR)
            FREK = FRQLRT(ILSTNR)
          IF (ISYM.EQ.1 .AND. DABS(FREK).LE.1.0D-10) THEN
C
            KT1  = 1
            KT2  = KT1  + NT1AM(ISYM)
            IF (.NOT. CCS) THEN
               KEND1= KT2  + NT2AM(ISYM)
            ELSE
               KEND1= KT2
            ENDIF
            IF (CCR12) THEN
               KTR12= KEND1 
               KEND1= KTR12 + NTR12AM(ISYM)
            ENDIF
            LEND1= LWORK - KEND1
            IF (LEND1 .LE. 0 ) CALL QUIT(
     *       'Insufficient workspace in cc_solrsp ')
            IOPT  = 3
            CALL CC_RDRSP('R1 ',ILSTNR,ISYM,IOPT,MODEL1,WORK(KT1),
     *                    WORK(KT2))
            IF (CCR12) THEN
              IOPT = 32
              CALL CC_RDRSP('R1 ',ILSTNR,ISYM,IOPT,MODEL1,DUMMY,
     *                    WORK(KTR12))
            ENDIF
            CALL CC_TSTAV(ILSTNR,WORK(KT1),WORK(KEND1),LEND1,0)
C
          END IF

         END DO
      ENDIF

      IF (DEBUG .OR. (IPRINT.GE.5)) THEN
         DO ILSTNR = 1, NO1LBL
C
            ISYM = ISYO1(ILSTNR)
            FREK = FRQO1(ILSTNR)
          IF (ISYM.EQ.1 .AND. DABS(FREK).LE.1.0D-10) THEN
C
            KT1  = 1
            KT2  = KT1  + NT1AM(ISYM)
            IF (.NOT. CCS) THEN
               KEND1= KT2  + NT2AM(ISYM)
            ELSE
               KEND1= KT2
            ENDIF
            IF (CCR12) THEN
               KTR12= KEND1
               KEND1= KTR12 + NTR12AM(ISYM)
            ENDIF
            LEND1= LWORK - KEND1
            IF (LEND1 .LE. 0 ) CALL QUIT(
     *       'Insufficient workspace in cc_solrsp ')
            IOPT  = 3
            IF (CCSDT) IOPT = 24
            CALL CC_RDRSP('O1 ',ILSTNR,ISYM,IOPT,MODEL1,WORK(KT1),
     *                    WORK(KT2))
            IF (CCR12) THEN
              IOPT = 32
              CALL CC_RDRSP('O1 ',ILSTNR,ISYM,IOPT,MODEL1,DUMMY,
     *                    WORK(KTR12))
            ENDIF
            CALL CC_TSTAV(ILSTNR,WORK(KT1),WORK(KEND1),LEND1,1)
C
          END IF

         END DO
      ENDIF
C
      IF (DEBUG .OR. (IPRINT.GE.5)) THEN
         DO ILSTNR = 1, NR2TLBL
            ISYM  = ILSTSYM('R2 ',ILSTNR)
            KT1   = 1
            KT2   = KT1  + NT1AM(ISYM)
            KEND1 = KT2
            IF (.NOT. CCS) KEND1 = KT2  + NT2AM(ISYM)
            IF (CCR12) THEN
               KTR12= KEND1
               KEND1= KTR12 + NTR12AM(ISYM)
            ENDIF
            LEND1 = LWORK - KEND1
            IF (LEND1.LE.0) CALL QUIT('Out of workspace in cc_solrsp')
            IOPT  = 3
            CALL CC_RDRSP('R2 ',ILSTNR,ISYM,IOPT,MODEL1,WORK(KT1),
     *                    WORK(KT2))
            IF (CCR12) THEN
              IOPT = 32
              CALL CC_RDRSP('R2 ',ILSTNR,ISYM,IOPT,MODEL1,DUMMY,
     *                    WORK(KTR12))
            ENDIF
            CALL CC_TSTAV2(ILSTNR,WORK(KT1),WORK(KEND1),LEND1,0)
         END DO
      ENDIF

C
      WRITE (LUPRI,*) 
      WRITE (LUPRI,*) 'Solution of CC response equations completed.'
      WRITE (LUPRI,*) 
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('CC_SOLRSP')
      RETURN
C
   2  FORMAT ((1X,A1,1X,A4,1X,A1,1X,I4,3X,A1,1X,A52,1X,A1))
      END
*=====================================================================*
      SUBROUTINE CC_TSTLST(LIST,NSTAT,ORDER,ISIDE,
     &                     ISYMS, ISTAT, EIGVAL, LPROJ,
     &                     ISYMO, LABEL, FREQS,  ICAU,
     &                     ISYOF, NVEC,  MAXVEC, WORK, LWORK)          
*---------------------------------------------------------------------*
*
*     Purpose: calculate and print norm of all vectors on a list
*              (mainly for testing and debugging purposes)
*
*     Christof Haettig, August 1999
*
*---------------------------------------------------------------------*
      IMPLICIT NONE
C                                                                
#include <priunit.h>
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <dummy.h>
Cholesky
#include <ccdeco.h>
      LOGICAL CCSEFF
Cholesky
C
      LOGICAL TRIPLET
      PARAMETER (TRIPLET= .FALSE.)       
C
      CHARACTER LIST*(3)

      INTEGER NSTAT, ORDER, MAXVEC, NVEC, ISIDE, LWORK
      INTEGER ISYOF(8), KT12, NCCVAR
      INTEGER ISYM, KT1, KT2, IOPT, ILSTNR, KEND1, LEND1, NTAMP

      LOGICAL LPROJ(MAXVEC)
      CHARACTER*8 LABEL(MAXVEC,ORDER)
      CHARACTER*10 MODEL
      INTEGER ISYMS(MAXVEC,NSTAT), ISTAT(MAXVEC,NSTAT)
      INTEGER ISYMO(MAXVEC,ORDER), ICAU(MAXVEC,ORDER)

#if defined (SYS_CRAY)
      REAL FREQS(MAXVEC,ORDER), EIGVAL(MAXVEC,NSTAT)
      REAL WORK(LWORK)
      REAL X1, X2, X12, DDOT
#else
      DOUBLE PRECISION FREQS(MAXVEC,ORDER), EIGVAL(MAXVEC,NSTAT)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION X1, X2, X12, DDOT
#endif
                                                            
      INTEGER ILSTSYM

      CCSEFF = CCS .OR. (CHOINT .AND. CC2)  ! No doubles for Chol. CC2, tbp
      
                                                                     
C
         IF (NVEC.LE.0) THEN
           WRITE (LUPRI,'(//1X,3A)')  'List ',LIST,' is empty.'
           RETURN
         END IF
C
         IF     (LIST(2:2).NE.'C'.AND.ORDER.EQ.1.AND.NSTAT.EQ.0) THEN
        
           WRITE (LUPRI,'(//1X,3A,I3,A,/1X,A)') 'Testing ',LIST,
     &        ' list with',NVEC, ' elements on it:',
     &        'List Nr.  Label     Symmetry      Frequency'
     &      //"      T1-norm        T2-norm     T2'-norm"     

         ELSEIF (LIST(1:1).NE.'C'.AND.ORDER.EQ.2.AND.NSTAT.EQ.0) THEN
        
           WRITE (LUPRI,'(//1X,3A,I3,A,/1X,A)') 'Testing ',LIST,
     &        ' list with', NVEC, ' elements on it:',
     &        'List Nr.  Label 1   Label 2   Symmetry      Freq. 1'
     &      //'     Freq. 2      T1-norm        T2-norm'     

         ELSEIF (LIST(2:2).EQ.'C'.AND.ORDER.EQ.1.AND.NSTAT.EQ.0) THEN

           WRITE (LUPRI,'(//1X,3A,I3,A,/1X,A)') 'Testing ',LIST,
     &        ' list with', NVEC, ' elements on it:',
     &        'List Nr.  Label     Symmetry     Order'
     &      //'    T1-norm             T2-norm'     

         ELSEIF (LIST(1:1).EQ.'C'.AND.ORDER.EQ.2.AND.NSTAT.EQ.0) THEN

           WRITE (LUPRI,'(//1X,3A,I3,A,/1X,A)') 'Testing ',LIST,
     &        ' list with', NVEC, ' elements on it:',
     &        'List Nr.  Label 1   Label 2   Symmetry   Order 1'
     &      //'   Order 2  T1-norm             T2-norm'     

         ELSEIF (ORDER.EQ.0 .AND. NSTAT.EQ.1) THEN

           WRITE (LUPRI,'(//1X,3A,I3,A,/1X,A)') 'Testing ',LIST,
     &        ' list with', NVEC, ' elements on it:',
     &        'List Nr.  State     Symmetry     Eigenvalue'
     &      //'      T1-norm        T2-norm'     

         ELSEIF (ORDER.EQ.1 .AND. NSTAT.EQ.1) THEN

           WRITE (LUPRI,'(//1X,3A,I3,A,/1X,A)') 'Testing ',LIST,
     &        ' list with', NVEC, ' elements on it:',
     &        'List Nr.  Label     State     Symmetry     Eigenval.'
     &      //'   Frequency    T1-norm        T2-norm'     

         ELSEIF (ORDER.EQ.0 .AND. NSTAT.EQ.2) THEN

           WRITE (LUPRI,'(//1X,3A,I3,A,/1X,A)') 'Testing ',LIST,
     &        ' list with', NVEC, ' elements on it:',
     &        'List Nr.  State 1   State 2   Symmetry     Eigval 1'
     &      //'    Eigval 2     T1-norm        T2-norm'     

         ELSE
           WRITE (LUPRI,*) 'Unknown list in CC_TSTLST... ignored.'
           RETURN
         ENDIF
C
         DO ILSTNR = 1, NVEC
 
           ISYM = ILSTSYM(LIST,ILSTNR)
           KT1  = 1 
           KEND1= KT1  + NT1AM(ISYM)
Chol       IF (.NOT. CCS) THEN
           IF (.NOT. CCSEFF) THEN
              KT2   = KEND1
              KEND1 = KT2  + NT2AM(ISYM)
           ENDIF 
           IF (CCR12) THEN
              KT12  = KEND1
              KEND1 = KT12  + NTR12AM(ISYM)
           END IF
           LEND1= LWORK - KEND1
           IF (LEND1 .LE. 0) THEN
               WRITE (LUPRI,*) 'Insufficient workspace in CC_TSTLST '
               CALL QUIT('Insufficient workspace in CC_TSTLST ')
           END IF

           IOPT  = 3
           CALL CC_RDRSP(LIST,ILSTNR,ISYM,IOPT,MODEL,WORK(KT1),
     &                   WORK(KT2))

           X1 = DDOT(NT1AM(ISYM),WORK(KT1),1,WORK(KT1),1)
Chol       IF (CCS) THEN
           IF (CCSEFF) THEN
              X2 = 0.0D0
           ELSE
              X2 = DDOT(NT2AM(ISYM),WORK(KT2),1,WORK(KT2),1)
           ENDIF

           IF (CCR12) THEN
            IOPT = 32
            CALL CC_RDRSP(LIST,ILSTNR,ISYM,IOPT,MODEL,DUMMY,WORK(KT12))
            X12 = DDOT(NTR12AM(ISYM),WORK(KT12),1,WORK(KT12),1)
           ELSE
            X12 = 0.0D0
           END IF

           IF     (LIST(2:2).NE.'C'.AND.ORDER.EQ.1.AND.NSTAT.EQ.0) THEN
        
             ! ground state freq.-dep. first order response
             WRITE (LUPRI,'(1X,I5,5X,A8,2X,I4,4X,4(2X,F15.10))')
     &            ILSTNR,LABEL(ILSTNR,1),ISYM,FREQS(ILSTNR,1),X1,X2,X12

           ELSEIF (LIST(1:1).NE.'C'.AND.ORDER.EQ.2.AND.NSTAT.EQ.0) THEN
        
             ! ground state freq.-dep. second order response
             WRITE (LUPRI,'(1X,I5,5X,2(A8,2X),I4,6X,2F12.7,2X,
     &             F15.10,F15.10)')
     &            ILSTNR,LABEL(ILSTNR,1),LABEL(ILSTNR,2),ISYM,
     &                    FREQS(ILSTNR,1),FREQS(ILSTNR,2),X1,X2

           ELSEIF (LIST(2:2).EQ.'C'.AND.ORDER.EQ.1.AND.NSTAT.EQ.0) THEN

             ! cauchy expansion of ground state first order response
             WRITE (LUPRI,'(1X,I5,5X,A8,2X,I4,6X,I6,2X,E20.10,E20.10)')
     &            ILSTNR,LABEL(ILSTNR,1),ISYM,ICAU(ILSTNR,1),X1,X2

           ELSEIF (LIST(1:1).EQ.'C'.AND.ORDER.EQ.2.AND.NSTAT.EQ.0) THEN

             ! cauchy expansion of ground state second order response
             WRITE (LUPRI,'(1X,I5,5X,2(A8,2X),I4,2I10,2X,E20.10,
     &             E20.10)')
     &            ILSTNR,LABEL(ILSTNR,1),LABEL(ILSTNR,2),ISYM,
     &                     ICAU(ILSTNR,1), ICAU(ILSTNR,2),X1,X2

           ELSEIF (ORDER.EQ.0 .AND. NSTAT.EQ.1) THEN

            ! excited state vectors, no response
            IF (LIST(1:3).EQ.'E0 '.OR.LIST(1:3).EQ.'BE ') THEN
             WRITE (LUPRI,'(1X,I5,5X,I4,4X,I4,6X,F15.10,F15.10,F15.10)')
     &            ILSTNR,ISTAT(ILSTNR,1),ISYM,
     &            EIGVAL(ISTAT(ILSTNR,1),1),X1,X2
            ELSE
             WRITE (LUPRI,'(1X,I5,5X,I4,4X,I4,6X,F15.10,F15.10,F15.10)')
     &            ILSTNR,ISTAT(ILSTNR,1),ISYM,EIGVAL(ILSTNR,1),X1,X2
            END IF

           ELSEIF (ORDER.EQ.1 .AND. NSTAT.EQ.1) THEN

             ! excited state first order response
             WRITE (LUPRI,'(1X,I5,5X,A8,2X,I4,4X,I4,6X,4F15.10)')
     &            ILSTNR,LABEL(ILSTNR,1),ISTAT(ILSTNR,1),ISYM,
     &                  EIGVAL(ILSTNR,1),FREQS(ILSTNR,1),X1,X2

           ELSEIF (ORDER.EQ.0 .AND. NSTAT.EQ.2) THEN

             ! excited state - excited state vectors, no response
             WRITE (LUPRI,'(1X,I5,5X,2(I4,4X),I4,6X,2F15.10,F15.10,
     &             F15.10)')
     &            ILSTNR,ISTAT(ILSTNR,1),ISTAT(ILSTNR,2),ISYM,
     &                   EIGVAL(ILSTNR,1),EIGVAL(ILSTNR,2),X1,X2

           ENDIF

         END DO

         IF (CHOINT .AND. CC2) THEN
            WRITE(LUPRI,'(/,1X,A,A)')
     &      'NOTICE: Cholesky CC2 doubles not stored! Artificial',
     &      ' zero norm printed.'
         ENDIF

         WRITE (LUPRI,'(//)')

         RETURN
         END 
*=====================================================================*
