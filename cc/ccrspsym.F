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
c /* deck ccrspsym */
*=====================================================================*
       SUBROUTINE CCRSPSYM(MOLGRD,WORK,LWRK)
*---------------------------------------------------------------------*
*
*    Purpose: symmetry checks for CC response calculations
*
*    Written by Christof Haettig, October 1996.
*    (Linear response residue setup Ove Christiansen 8-11-1996)
*    (New Linear response residue setup Ove Christiansen 23-4-1997)
*    (PL1 vectors and relaxation in EL1, Sonia Coriani, March 2000)
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <nuclei.h>
#include <symmet.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <ccorb.h>
#include <ccrspprp.h>
#include <cclrinf.h>
#include <ccroper.h>
#include <ccropr2.h>
#include <ccexpfck.h>
#include <cc1dxfck.h>
#include <cclrmrsp.h>
#include <ccer1rsp.h>
#include <ccer2rsp.h>
#include <ccel1rsp.h>
#include <ccel2rsp.h>
#include <ccr1rsp.h>
#include <ccr2rsp.h>
#include <ccr3rsp.h>
#include <ccr4rsp.h>
#include <ccl1rsp.h>
#include <ccl2rsp.h>
#include <ccl3rsp.h>
#include <ccl4rsp.h>
#include <ccx1rsp.h>  
#include <ccx2rsp.h>  
#include <ccx3rsp.h>  
#include <ccx4rsp.h>  
#include <cco1rsp.h>  
#include <cco2rsp.h>  
#include <cco3rsp.h>  
#include <cco4rsp.h>  
#include <ccrc1rsp.h>
#include <cclc1rsp.h>
#include <cccr2rsp.h>
#include <ccco2rsp.h>
#include <cccl2rsp.h>
#include <cccx2rsp.h>
#include <ccexgr.h>  
#include <ccn2rsp.h>  
#include <cclres.h>  
#include <ccpl1rsp.h>
#include <ccexci.h>
Cholesky
#include <cclrcho.h>
Cholesky

* local parameters:
      CHARACTER MSGDBG*(18)
      PARAMETER (MSGDBG='[debug] CCRSPSYM> ')
      CHARACTER SECNAM*(8)
      PARAMETER (SECNAM='CCRSPSYM')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
  
      REAL*8  ZERO, ONE, TWO, EIGHT
      REAL*8  TINY
      REAL*8  CKMXPR
      PARAMETER (ZERO = 0.0d0, ONE = 1.0d0, TWO = 2.0d0, EIGHT = 8.0d0)
      PARAMETER (TINY = 1.0d-5, CKMXPR = 1.0d-12)


* variables:
      LOGICAL LDUM, MOLGRD
      CHARACTER*8 LABEL,CDUM
      INTEGER I, IND, IIND, IIMAX, JIND, KSYMPT, IPRP, IERR, ISYM, IDX
      INTEGER LWRK, IDUM, IL, INUM, J, IVEC, IMATRIX, ICAU,ICAU1,ICAU2
      INTEGER KEND1,LEND1,KPROPAO,JSCOOR,ISCOOR,ICORSY,ICOOR,IATOM
      INTEGER ISYM0

      REAL*8  WORK(LWRK),RDUM

* external functions:
      INTEGER ILSTSYM
      INTEGER INDPRPCC
      INTEGER IROPER
      INTEGER IROPER2
      INTEGER IETA1
      INTEGER ICHI2
      INTEGER ICHI3
      INTEGER ICHI4
      INTEGER IRHSR1
      INTEGER IRHSR2
      INTEGER IRHSR3
      INTEGER IRHSR4
      INTEGER IR1KAPPA
      INTEGER IR1TAMP
      INTEGER IR2TAMP
      INTEGER IR3TAMP
      INTEGER IR4TAMP
      INTEGER IL1ZETA
      INTEGER IL2ZETA
      INTEGER IL3ZETA
      INTEGER IL4ZETA
      INTEGER IER1AMP
      INTEGER IER2AMP
      INTEGER IEL1AMP
      INTEGER IEL2AMP
      INTEGER ILRCAMP
      INTEGER ILC1AMP
      INTEGER ICR2AMP
      INTEGER ICL2AMP
      INTEGER IRHSCR2
      INTEGER IETACL2
      INTEGER IEFFFOCK
      INTEGER IPL1ZETA  

* data:
      LOGICAL FIRSTCALL
      SAVE    FIRSTCALL
      DATA    FIRSTCALL /.TRUE./
      

*---------------------------------------------------------------------*
* print header:
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(3X,A,/)') '  '
      WRITE (LUPRI,'(3X,A)')
     *'*********************************************************'//
     *'**********'
      WRITE (LUPRI,'(3X,A)')
     *'*                                                        '//
     *'         *'
      WRITE (LUPRI,'(3X,A)')
     *'*<<<<<<<<   OUTPUT FROM PROPERTY AND SYMMETRY ANALYSIS   '//
     *'>>>>>>>>>*'
      WRITE (LUPRI,'(3X,A)')
     *'*                                                        '//
     *'         *'
      WRITE (LUPRI,'(3X,A,/)')
     *'*********************************************************'//
     *'**********'

*---------------------------------------------------------------------*
* initializations
*---------------------------------------------------------------------*
* init number of response operators:
      NRSOLBL  = 0
      NRSO2LBL = 0
Cholesky
* number of unique operators for Cholesky CC2-LR
      CALL IZERO(NCHOPLR,NSYM)
Cholesky
      LOPROPN  = .TRUE.  ! open list for extension
      LOPR2OPN = .TRUE.  ! open list for extension

* init number of ground state response equations for t and zeta:
      NLRTLBL = 0
      NR2TLBL = 0
      NR3TLBL = 0
      NR4TLBL = 0
      NLRZLBL = 0
      NL2LBL  = 0
      NL3LBL  = 0
      NL4LBL  = 0

* init number of groud state response equations for kappa:
      NLRTHFLBL = 0

* init number of ground state response equations for projected zeta:
      NPL1LBL = 0

* init number of cauchy equations:
      NLRCLBL = 0
      NLC1LBL = 0
      NCR2LBL = 0
      NCL2LBL = 0

* init number of multipliers for oscillator strengths:
      NLRM    = 0
      NQRN2   = 0

* init number of rhs/chi vectors for ground state 
* response and Cauchy equations:
      NO1LBL  = 0
      NO2LBL  = 0
      NO3LBL  = 0
      NO4LBL  = 0
      NX1LBL  = 0
      NX2LBL  = 0
      NX3LBL  = 0
      NX4LBL  = 0
      NCO2LBL = 0
      NCX2LBL = 0

* init number of left/right excited state vector response equations.
      NER1LBL = 0
      NER2LBL = 0
      NEL1LBL = 0
      NEL2LBL = 0

* init number of effective Fock matrices from one-index transformed
* integrals:
      N1DXFLBL = 0

* make sure that the zeroth-order Hamiltonian 'HAM0    ' is on our
* operator list (might be needed for analytic derivatives):
      INUM = INDPRPCC('HAM0    ')

* allocate work space for one set of property AO integrals
      KPROPAO = 1
      KEND1   = KPROPAO + N2BASX
      LEND1   = LWRK - KEND1

      IF (LEND1 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CCRSPSYM')
      END IF

*---------------------------------------------------------------------*
* loop over property labels in the PRPLBL list:
*---------------------------------------------------------------------*
      DO IPRP = 1, NPRLBL
        LABEL = PRPLBL(IPRP)
        
        IF (LABEL(1:5).EQ.'HAM0 ') THEN 
           IERR    = 0
           KSYMPT  = 1
           IMATRIX = 1
        ELSE
           CALL CCPRPAO(LABEL,.TRUE.,WORK(KPROPAO),KSYMPT,IMATRIX,IERR,
     &                  WORK(KEND1),LEND1)
          
        END IF

        IF (IERR.EQ.0 .AND. KSYMPT.GT.0 .AND. KSYMPT.LE.NSYM) THEN

C          ------------------------------------------------------
C          build a list with labels, symmetries and orbital
C          relaxation flags of requested AND available operators:
C          ------------------------------------------------------
           INUM  = IROPER(LABEL,KSYMPT)

C          ----------------------------------------------------
C          save symmetry of integral matrix: 
C            (symmetric=+1 / antisymmetric=-1)
C          ----------------------------------------------------
           ISYMAT(INUM) = IMATRIX

C          ----------------------------------------------------
C          decide whether basis set depends on the perturbation:
C          (default is no PDBS, at present we have PDBS only 
C           for geometric first derivatives):
C          ----------------------------------------------------
           IF ( LABEL(1:5).EQ.'HAM0 ' ) THEN
              LPDBSOP(INUM) = .TRUE.
           ELSE IF ( LABEL(1:5).EQ.'1DHAM' ) THEN
              LPDBSOP(INUM) = .TRUE.
           ELSE IF ( LABEL(1:5).EQ.'dh/dB' ) THEN
              LPDBSOP(INUM) = .TRUE.
           ELSE
              LPDBSOP(INUM) = .FALSE.
           END IF

C          ----------------------------------------------------
C          determine index of associate 'Atom':
C          ----------------------------------------------------
           IF ( LABEL(1:5).EQ.'1DHAM' .OR. LABEL(4:6).EQ.'DPG') THEN
              IF (LABEL(1:5).EQ.'1DHAM') READ(LABEL,'(5X,I3)') JSCOOR
              IF (LABEL(4:6).EQ.'DPG')   READ(LABEL,'(I3)')    JSCOOR
 
              DO IATOM = 1, NUCIND
                 DO ICORSY = 1, NSYM
                    DO ICOOR = 1, 3
                       ISCOOR = IPTCNT(3*(IATOM-1)+ICOOR,ICORSY-1,1)
                       IF (ISCOOR.EQ.JSCOOR) THEN
                          IATOPR(INUM) = IATOM
                       END IF
                    END DO
                 END DO
              END DO

           ELSE
              IATOPR(INUM) = 0
           END IF

        END IF
      END DO


* close list of operators for extension and sort it:
      LOPROPN  = .FALSE.  ! close list for extension
      LOPR2OPN = .FALSE.  ! close list for extension
      LQUIET   = .FALSE.  ! warn if problems in IROPER2

      IF (LOCDBG .AND. NRSOLBL.GT.0) THEN
        WRITE(LUPRI,'(/A)') ' UNSORTED LIST OF REQUIRED OPERATORS:'
        DO I = 1, NRSOLBL
          WRITE(LUPRI,'(I5,3X,A8,2I5,L3,I5)') I, LBLOPR(I), 
     &             ISYOPR(I), ISYMAT(I), LPDBSOP(I), IATOPR(I)
        END DO
      END IF

      CALL CCLSTSORT('o1 ',IDUM, IDUM, RDUM, ISYOPR,LBLOPR,RDUM,IDUM,
     &                    LDUM, ISYOFO,  NRSOLBL, MAXOLBL, LDUM )

*---------------------------------------------------------------------*
* initialize list of expectation values/effective fock matrices
*---------------------------------------------------------------------*
      CALL CC_EXPFCK_INIT(MOLGRD)

*---------------------------------------------------------------------*
* set up the lists of response/cauchy equations to be solved:
*---------------------------------------------------------------------*
* open lists:
      LR1OPN   = .TRUE.
      LL1OPN   = .TRUE.
      LX1OPN   = .TRUE.
      LO1OPN   = .TRUE.
      LR2OPN   = .TRUE.
      LL2OPN   = .TRUE.
      LX2OPN   = .TRUE.
      LO2OPN   = .TRUE.
      LR3OPN   = .TRUE.
      LL3OPN   = .TRUE.
      LX3OPN   = .TRUE.
      LO3OPN   = .TRUE.
      LR4OPN   = .TRUE.
      LL4OPN   = .TRUE.
      LX4OPN   = .TRUE.
      LO4OPN   = .TRUE.
      LN2OPN   = .TRUE.
      LER1OPN  = .TRUE.
      LER2OPN  = .TRUE.
      LEL1OPN  = .TRUE.
      LEL2OPN  = .TRUE.
      LRC1OPN  = .TRUE.
      LLC1OPN  = .TRUE.
      LCR2OPN  = .TRUE.
      LCO2OPN  = .TRUE.
      LCL2OPN  = .TRUE.
      LCX2OPN  = .TRUE.
      LEXPTOPN = .TRUE.
      LEFCKOPN = .TRUE.
      L1DXFOPN = .TRUE.
      LPL1OPN  = .TRUE.

* linear response equations for R1 or M1 vectors required for linear 
* response residues (one-photon transition moments for ground
* to excited state transitions)
      CALL CC_LRSIND
      IF (CCOPA) CALL CC_OPAIND

* linear response equations for R1 or N2 vectors required for quadratic 
* response second residues (one-photon transition moments for
* excited to excited state transitions)
      IF (CCQR2R) CALL CC_QR2IND
      IF (CCXOPA) CALL CC_XOPAIND

* test input for excited state first-order property calculation.
      IF (CCEXGR) CALL CC_EXGRIND  

* set equations for excited state second-order properties:
      IF (CCEXLR) CALL CC_EXLRIND 

* set equations required for second-order transition moments:
      IF (CCTPA)  CALL CC_TPAIND

* set equations required for third-order transition moments:
      CALL CC_TMIND

* set equations required for MCD section:
      CALL CC_MCDIND(WORK,LWRK)

* set equations required for polarizabilities and Cauchy moments:
      CALL CC_LRIND(WORK,LWRK)

* linear response equations for t amplitudes and zeta multipliers
* required for the first hyperpolarizabilities
      CALL CC_QRIND(WORK,LWRK)

* linear and quadratic response equations for t amplitudes and 
* zeta multipliers required for the second hyperpolarizabilities
      CALL CC_CRIND

* linear and quadratic response equations for t amplitudes and 
* zeta multipliers required for the third hyperpolarizabilities
      CALL CC_4RIND

* first-, second- and third-order response equations for t amplitudes
* and zeta multipliers required for the fourth hyperpolarizabilities
      CALL CC_5RIND


*=====================================================================*
* make response/rhs vector lists consistent:
* (uses a waterfall strategy, so the order is important!)
*=====================================================================*

*---------------------------------------------------------------------*
* request fourth-order chi (X4) vectors for all entries in the
* fourth-order zeta multiplier (L4) list:
*---------------------------------------------------------------------*
      DO IVEC = 1, NL4LBL
        INUM = ICHI4(LBLL4(IVEC,1),FRQL4(IVEC,1),ISYL4(IVEC,1),
     &               LBLL4(IVEC,2),FRQL4(IVEC,2),ISYL4(IVEC,2),
     &               LBLL4(IVEC,3),FRQL4(IVEC,3),ISYL4(IVEC,3),
     &               LBLL4(IVEC,4),FRQL4(IVEC,4),ISYL4(IVEC,4) )
      END DO

*---------------------------------------------------------------------*
* request fourth-order amplitude (R4) vectors for all entries in the
* fourth-order multiplier (L4) vector lists:
*---------------------------------------------------------------------*
      DO IVEC = 1, NL4LBL
        INUM = IR4TAMP(LBLL4(IVEC,1),FRQL4(IVEC,1),ISYL4(IVEC,1),
     &                 LBLL4(IVEC,2),FRQL4(IVEC,2),ISYL4(IVEC,2),
     &                 LBLL4(IVEC,3),FRQL4(IVEC,3),ISYL4(IVEC,3),
     &                 LBLL4(IVEC,4),FRQL4(IVEC,4),ISYL4(IVEC,4) )
      END DO

*---------------------------------------------------------------------*
* request third-order multipliers (L3) vectors for all entries in the 
* fourth-order chi (X4) vector list:
*---------------------------------------------------------------------*
      DO IVEC = 1, NX4LBL
        INUM = IL3ZETA(LBLX4(IVEC,1),FRQX4(IVEC,1),ISYX4(IVEC,1),
     &                 LBLX4(IVEC,2),FRQX4(IVEC,2),ISYX4(IVEC,2),
     &                 LBLX4(IVEC,3),FRQX4(IVEC,3),ISYX4(IVEC,3))

        INUM = IL3ZETA(LBLX4(IVEC,1),FRQX4(IVEC,1),ISYX4(IVEC,1),
     &                 LBLX4(IVEC,2),FRQX4(IVEC,2),ISYX4(IVEC,2),
     &                 LBLX4(IVEC,4),FRQX4(IVEC,4),ISYX4(IVEC,4))

        INUM = IL3ZETA(LBLX4(IVEC,1),FRQX4(IVEC,1),ISYX4(IVEC,1),
     &                 LBLX4(IVEC,3),FRQX4(IVEC,3),ISYX4(IVEC,3),
     &                 LBLX4(IVEC,4),FRQX4(IVEC,4),ISYX4(IVEC,4))

        INUM = IL3ZETA(LBLX4(IVEC,2),FRQX4(IVEC,2),ISYX4(IVEC,2),
     &                 LBLX4(IVEC,3),FRQX4(IVEC,3),ISYX4(IVEC,3),
     &                 LBLX4(IVEC,4),FRQX4(IVEC,4),ISYX4(IVEC,4))
      END DO

*---------------------------------------------------------------------*
* request third-order chi (X3) vectors for all entries in the 
* third-order zeta multiplier (L3) vector list:
*---------------------------------------------------------------------*
      DO IVEC = 1, NL3LBL
        INUM = ICHI3(LBLL3(IVEC,1),FRQL3(IVEC,1),ISYL3(IVEC,1),
     &               LBLL3(IVEC,2),FRQL3(IVEC,2),ISYL3(IVEC,2),
     &               LBLL3(IVEC,3),FRQL3(IVEC,3),ISYL3(IVEC,3))
      END DO

*---------------------------------------------------------------------*
* request third-order amplitude (R3) vectors for all entries in the
* L3 and O4 lists:
*---------------------------------------------------------------------*
      DO IVEC = 1, NL3LBL
        INUM = IR3TAMP(LBLL3(IVEC,1),FRQL3(IVEC,1),ISYL3(IVEC,1),
     &                 LBLL3(IVEC,2),FRQL3(IVEC,2),ISYL3(IVEC,2),
     &                 LBLL3(IVEC,3),FRQL3(IVEC,3),ISYL3(IVEC,3))
      END DO

      DO IVEC = 1, NO4LBL
        INUM = IR3TAMP(LBLO4(IVEC,1),FRQO4(IVEC,1),ISYO4(IVEC,1),
     &                 LBLO4(IVEC,2),FRQO4(IVEC,2),ISYO4(IVEC,2),
     &                 LBLO4(IVEC,3),FRQO4(IVEC,3),ISYO4(IVEC,3))

        INUM = IR3TAMP(LBLO4(IVEC,1),FRQO4(IVEC,1),ISYO4(IVEC,1),
     &                 LBLO4(IVEC,2),FRQO4(IVEC,2),ISYO4(IVEC,2),
     &                 LBLO4(IVEC,4),FRQO4(IVEC,4),ISYO4(IVEC,4))

        INUM = IR3TAMP(LBLO4(IVEC,1),FRQO4(IVEC,1),ISYO4(IVEC,1),
     &                 LBLO4(IVEC,3),FRQO4(IVEC,3),ISYO4(IVEC,3),
     &                 LBLO4(IVEC,4),FRQO4(IVEC,4),ISYO4(IVEC,4))

        INUM = IR3TAMP(LBLO4(IVEC,2),FRQO4(IVEC,2),ISYO4(IVEC,2),
     &                 LBLO4(IVEC,3),FRQO4(IVEC,3),ISYO4(IVEC,3),
     &                 LBLO4(IVEC,4),FRQO4(IVEC,4),ISYO4(IVEC,4))
      END DO

*---------------------------------------------------------------------*
* request third-order amplitude rhs (O3) vectors for all entries in
* the third-order amplitude (R3) list:
*---------------------------------------------------------------------*
      DO IVEC = 1, NR3TLBL
        INUM = IRHSR3(LBLR3T(IVEC,1),FRQR3T(IVEC,1),ISYR3T(IVEC,1),
     &                LBLR3T(IVEC,2),FRQR3T(IVEC,2),ISYR3T(IVEC,2),
     &                LBLR3T(IVEC,3),FRQR3T(IVEC,3),ISYR3T(IVEC,3))
      END DO

*---------------------------------------------------------------------*
* request second-order multiplier (L2) vectors for all entries in the
* third-order chi (X3), static vectors for all entries in the CL2
* list, and for all second-order left excited state (EL2) vectors:
*---------------------------------------------------------------------*
      DO IVEC = 1, NX3LBL
        INUM = IL2ZETA(LBLX3(IVEC,1),FRQX3(IVEC,1),ISYX3(IVEC,1),
     &                 LBLX3(IVEC,2),FRQX3(IVEC,2),ISYX3(IVEC,2))

        INUM = IL2ZETA(LBLX3(IVEC,1),FRQX3(IVEC,1),ISYX3(IVEC,1),
     &                 LBLX3(IVEC,3),FRQX3(IVEC,3),ISYX3(IVEC,3))

        INUM = IL2ZETA(LBLX3(IVEC,2),FRQX3(IVEC,2),ISYX3(IVEC,2),
     &                 LBLX3(IVEC,3),FRQX3(IVEC,3),ISYX3(IVEC,3))
      END DO

      DO IVEC = 1, NCL2LBL
        INUM = IL2ZETA(LBLCL2(IVEC,1),0.0d0,ISYCL2(IVEC,1),
     &                 LBLCL2(IVEC,2),0.0d0,ISYCL2(IVEC,2))
      END DO

      DO IVEC = 1, NEL2LBL
        INUM = IL2ZETA(LBLEL2(IVEC,1),FRQEL2(IVEC,1),ISYOEL2(IVEC,1),
     &                 LBLEL2(IVEC,2),FRQEL2(IVEC,2),ISYOEL2(IVEC,2))
      END DO

*---------------------------------------------------------------------*
* request second-order Cauchy eta (CX2) vectors for all entries in the
* second-order left Cauchy (CL2) vector list:
*---------------------------------------------------------------------*
      DO IVEC = 1, NCL2LBL
        INUM = IETACL2(LBLCL2(IVEC,1),ICL2CAU(IVEC,1),ISYCL2(IVEC,1),
     &                 LBLCL2(IVEC,2),ICL2CAU(IVEC,2),ISYCL2(IVEC,2) )
      END DO

*---------------------------------------------------------------------*
* request second-order eta (X2) vectors for all entries in the
* second-order zeta multiplier (L2) list, and static vectors for
* all entries in the CX2 list:
*---------------------------------------------------------------------*
      DO IVEC = 1, NL2LBL
        INUM = ICHI2(LBLAL2(IVEC),.FALSE.,FRQAL2(IVEC),ISYAL2(IVEC),
     &               LBLBL2(IVEC),.FALSE.,FRQBL2(IVEC),ISYBL2(IVEC) )
      END DO

      DO IVEC = 1, NCX2LBL
        INUM = ICHI2(LBLCX2(IVEC,1),.FALSE.,0.0d0,ISYCX2(IVEC,1),
     &               LBLCX2(IVEC,2),.FALSE.,0.0d0,ISYCX2(IVEC,2) )
      END DO
*---------------------------------------------------------------------*
* request second-order right Cauchy vectors for all entries in the
* second-order left Cauchy vector list and for all entries in the
* second-order right Cauchy vector list with higher cauchy order:
*---------------------------------------------------------------------*
      DO IVEC = 1, NCL2LBL
        INUM = ICR2AMP(LBLCL2(IVEC,1),ICL2CAU(IVEC,1),ISYCL2(IVEC,1),
     &                 LBLCL2(IVEC,2),ICL2CAU(IVEC,2),ISYCL2(IVEC,2) )
      END DO

      DO IVEC = 1, NCR2LBL
        DO ICAU1 = 0, ICR2CAU(IVEC,1)
        DO ICAU2 = 0, ICR2CAU(IVEC,2)
          IF ((ICAU1+ICAU2).GT.0) THEN
            INUM = ICR2AMP(LBLCR2(IVEC,1),ICAU1,ISYCR2(IVEC,1),
     &                     LBLCR2(IVEC,2),ICAU2,ISYCR2(IVEC,2))
          END IF
        END DO
        END DO
      END DO

*---------------------------------------------------------------------*
* request second-order right Cauchy rhs (CO2) vectors for all entries 
* in the second-order right Cauchy (R2) list:
*---------------------------------------------------------------------*
      DO IVEC = 1, NCR2LBL
        INUM = IRHSCR2(LBLCR2(IVEC,1),ICR2CAU(IVEC,1),ISYCR2(IVEC,1),
     &                 LBLCR2(IVEC,2),ICR2CAU(IVEC,2),ISYCR2(IVEC,2) )
      END DO

*---------------------------------------------------------------------*
* request second-order amplitude (R2) vectors for all entries in the
* second-order multiplier (L2), third-order rhs (O3) and second-order
* right cauchy (CR2) and right excited state (ER2) vector lists:
*---------------------------------------------------------------------*
      DO IVEC = 1, NL2LBL
        INUM=IR2TAMP(LBLAL2(IVEC),.FALSE.,FRQAL2(IVEC),ISYAL2(IVEC),
     &               LBLBL2(IVEC),.FALSE.,FRQBL2(IVEC),ISYBL2(IVEC) )
      END DO

      DO IVEC = 1, NO3LBL
        INUM=IR2TAMP(LBLO3(IVEC,1),.FALSE.,FRQO3(IVEC,1),ISYO3(IVEC,1),
     &               LBLO3(IVEC,2),.FALSE.,FRQO3(IVEC,2),ISYO3(IVEC,2))

        INUM=IR2TAMP(LBLO3(IVEC,1),.FALSE.,FRQO3(IVEC,1),ISYO3(IVEC,1),
     &               LBLO3(IVEC,3),.FALSE.,FRQO3(IVEC,3),ISYO3(IVEC,3))

        INUM=IR2TAMP(LBLO3(IVEC,2),.FALSE.,FRQO3(IVEC,2),ISYO3(IVEC,2),
     &               LBLO3(IVEC,3),.FALSE.,FRQO3(IVEC,3),ISYO3(IVEC,3))
      END DO

      DO IVEC = 1, NCR2LBL
        INUM=IR2TAMP(LBLCR2(IVEC,1),.FALSE.,0.0d0,ISYCR2(IVEC,1),
     &               LBLCR2(IVEC,2),.FALSE.,0.0d0,ISYCR2(IVEC,2))
      END DO

      DO IVEC = 1, NER2LBL
        INUM=IR2TAMP(LBLER2(IVEC,1),.FALSE.,FRQER2(IVEC,1),
     &                                              ISYOER2(IVEC,1),
     &               LBLER2(IVEC,2),.FALSE.,FRQER2(IVEC,2),
     &                                              ISYOER2(IVEC,2))
      END DO

*---------------------------------------------------------------------*
* request second-order rhs (O2) vectors for all entries in the
* second-order amplitude (R2) list and static vectors for all entries
* in the CO2 list:
*---------------------------------------------------------------------*
      DO IVEC = 1, NR2TLBL
        INUM = IRHSR2(LBLAR2T(IVEC),.FALSE.,FRQAR2T(IVEC),ISYAR2T(IVEC),
     &                LBLBR2T(IVEC),.FALSE.,FRQBR2T(IVEC),ISYBR2T(IVEC))
      END DO

      DO IVEC = 1, NCR2LBL
        INUM = IRHSR2(LBLCR2(IVEC,1),.FALSE.,0.0d0,ISYCR2(IVEC,1),
     &                LBLCR2(IVEC,2),.FALSE.,0.0d0,ISYCR2(IVEC,2))
      END DO

*---------------------------------------------------------------------*
* request first-order left excited state response vectors (EL1) for 
* all entries in the second-order left excited state (EL2) list:
*---------------------------------------------------------------------*
      DO IVEC = 1, NEL2LBL
        INUM = IEL1AMP(ISTEL2(IVEC),  EIGEL2(IVEC),  ISYSEL2(IVEC),
     &                 LBLEL2(IVEC,1),FRQEL2(IVEC,1),ISYOEL2(IVEC,1),
     &                 .FALSE.,LPREL2(IVEC)                          )
        INUM = IEL1AMP(ISTEL2(IVEC),  EIGEL2(IVEC),  ISYSEL2(IVEC),
     &                 LBLEL2(IVEC,2),FRQEL2(IVEC,2),ISYOEL2(IVEC,2),
     &                 .FALSE.,LPREL2(IVEC)                          )
      END DO

*---------------------------------------------------------------------*
* request first-order right excited state response vectors (ER1) for 
* all entries in the second-order right excited state (ER2) list:
*---------------------------------------------------------------------*
      DO IVEC = 1, NER2LBL
        INUM = IER1AMP(ISTER2(IVEC),  EIGER2(IVEC),  ISYSER2(IVEC),
     &                 LBLER2(IVEC,1),FRQER2(IVEC,1),ISYOER2(IVEC,1),
     &                 LPRER2(IVEC)                                  )
        INUM = IER1AMP(ISTER2(IVEC),  EIGER2(IVEC),  ISYSER2(IVEC),
     &                 LBLER2(IVEC,2),FRQER2(IVEC,2),ISYOER2(IVEC,2),
     &                 LPRER2(IVEC)                                  )
      END DO

*---------------------------------------------------------------------*
* request left first-order cauchy vectors for all entries in the
* second-order Cauchy (CL2) and second-order Cauchy eta (CX2) vectors
* lists and all left Cauchy vectors with higher cauchy order:
*---------------------------------------------------------------------*
      DO IVEC = 1, NCL2LBL
        IF (ICL2CAU(IVEC,1).GT.0) 
     &    INUM = ILC1AMP(LBLCL2(IVEC,1),ICL2CAU(IVEC,1),ISYCL2(IVEC,1))
        IF (ICL2CAU(IVEC,2).GT.0) 
     &    INUM = ILC1AMP(LBLCL2(IVEC,2),ICL2CAU(IVEC,2),ISYCL2(IVEC,2))
      END DO

      DO IVEC = 1, NCX2LBL
        IF (ICX2CAU(IVEC,1).GT.0) 
     &    INUM = ILC1AMP(LBLCX2(IVEC,1),ICX2CAU(IVEC,1),ISYCX2(IVEC,1))
        IF (ICX2CAU(IVEC,2).GT.0) 
     &    INUM = ILC1AMP(LBLCX2(IVEC,2),ICX2CAU(IVEC,2),ISYCX2(IVEC,2))
      END DO

      DO IVEC = 1, NLC1LBL
        DO ICAU = 1, ILC1CAU(IVEC)-1
          INUM = ILC1AMP(LBLLC1(IVEC),ICAU,ISYLC1(IVEC))
        END DO
      END DO

*---------------------------------------------------------------------*
* request first-order right Cauchy vectors for all entries in the
* first-order left Cauchy vector and second-order right Cauchy vector
* list, second-order Cauchy rhs vector list and for all entries in 
* the first-order right Cauchy vector list with higher cauchy order:
*---------------------------------------------------------------------*
      DO IVEC = 1, NLC1LBL
        INUM = ILRCAMP(LBLLC1(IVEC),ILC1CAU(IVEC),ISYLC1(IVEC))
      END DO

      DO IVEC = 1, NCR2LBL
        IF (ICR2CAU(IVEC,1).GT.0) 
     &    INUM = ILRCAMP(LBLCR2(IVEC,1),ICR2CAU(IVEC,1),ISYCR2(IVEC,1))
        IF (ICR2CAU(IVEC,2).GT.0) 
     &    INUM = ILRCAMP(LBLCR2(IVEC,2),ICR2CAU(IVEC,2),ISYCR2(IVEC,2))
      END DO

      DO IVEC = 1, NCO2LBL
        IF (ICO2CAU(IVEC,1).GT.0) 
     &    INUM = ILRCAMP(LBLCO2(IVEC,1),ICO2CAU(IVEC,1),ISYCO2(IVEC,1))
        IF (ICO2CAU(IVEC,2).GT.0) 
     &    INUM = ILRCAMP(LBLCO2(IVEC,2),ICO2CAU(IVEC,2),ISYCO2(IVEC,2))
      END DO

      DO IVEC = 1, NLRCLBL
        DO ICAU = 1, ILRCAU(IVEC)-1
          INUM = ILRCAMP(LRCLBL(IVEC),ICAU,ISYLRC(IVEC))
        END DO
      END DO

*---------------------------------------------------------------------*
* request (unrelaxed) first-order multipliers for all entries in the 
* second-order eta (X2) and for all left cauchy vectors:
*---------------------------------------------------------------------*
      DO IVEC = 1, NX2LBL
        INUM = IL1ZETA(LBLAX2(IVEC),.FALSE.,FRQAX2(IVEC),ISYAX2(IVEC))
        INUM = IL1ZETA(LBLBX2(IVEC),.FALSE.,FRQBX2(IVEC),ISYBX2(IVEC))
      END DO

      DO IVEC = 1, NLC1LBL
        INUM = IL1ZETA(LBLLC1(IVEC),.FALSE.,0.0d0,ISYLC1(IVEC))
      END DO

*---------------------------------------------------------------------*
* request (unrelaxed) first-order amplitude response for all entries in
* the second-order rhs (O2), first-order multiplier (L1), first-order
* left and right excited state (EL1/ER1), projected 1st-order 
* multipliers (PL1) and first-order right cauchy vector lists:
*---------------------------------------------------------------------*
      DO IVEC = 1, NO2LBL
        INUM = IR1TAMP(LBLAO2(IVEC),.FALSE.,FRQAO2(IVEC),ISYAO2(IVEC))
        INUM = IR1TAMP(LBLBO2(IVEC),.FALSE.,FRQBO2(IVEC),ISYBO2(IVEC))
      END DO

      DO IVEC = 1, NLRZLBL
        INUM = IR1TAMP(LRZLBL(IVEC),LORXLRZ(IVEC),
     &                 FRQLRZ(IVEC),ISYLRZ(IVEC))
      END DO

      DO IVEC = 1, NER1LBL
        INUM = IR1TAMP(LBLER1(IVEC),.FALSE.,FRQER1(IVEC),ISYOER1(IVEC))
      END DO

      DO IVEC = 1, NEL1LBL
        INUM = IR1TAMP(LBLEL1(IVEC),LORXEL1(IVEC),FRQEL1(IVEC),
     &                                            ISYOEL1(IVEC))
      END DO

      DO IVEC = 1, NPL1LBL
        INUM = IR1TAMP(LBLPL1(IVEC),LORXPL1(IVEC),
     &                 FRQPL1(IVEC),ISYPL1(IVEC))
      END DO

      DO IVEC = 1, NLRCLBL
        INUM = IR1TAMP(LRCLBL(IVEC),.FALSE.,0.0d0,ISYLRC(IVEC))
      END DO

*---------------------------------------------------------------------*
* request right hand side vector for first-order Lagrangian multiplier
* response equations for all entries in the L1 and PL1 lists:
*---------------------------------------------------------------------*
      DO IVEC = 1, NLRZLBL
        INUM=IETA1(LRZLBL(IVEC),LORXLRZ(IVEC),FRQLRZ(IVEC),ISYLRZ(IVEC))
      END DO
      DO IVEC = 1, NPL1LBL
        INUM=IETA1(LBLPL1(IVEC),LORXPL1(IVEC),FRQPL1(IVEC),ISYPL1(IVEC))
      END DO

*---------------------------------------------------------------------*
* request right hand side vector for first-order amplitude response 
* equations for all entries in the R1 list:
*---------------------------------------------------------------------*
      DO IVEC = 1, NLRTLBL
       INUM=IRHSR1(LRTLBL(IVEC),LORXLRT(IVEC),FRQLRT(IVEC),ISYLRT(IVEC))
      END DO

*---------------------------------------------------------------------*
* for all elements of the O1 and X1 lists request the corresponding
* CPHF response equations:
*---------------------------------------------------------------------*
      DO IVEC = 1, NX1LBL
        IF (LORXX1(IVEC)) THEN
          INUM = IR1KAPPA(LBLX1(IVEC),FRQX1(IVEC),ISYX1(IVEC))
        END IF
      END DO

      DO IVEC = 1, NO1LBL
        IF (LORXO1(IVEC)) THEN
          INUM = IR1KAPPA(LBLO1(IVEC),FRQO1(IVEC),ISYO1(IVEC))
        END IF
      END DO

*---------------------------------------------------------------------*
* for all CPHF equations request RHS vectors:
*---------------------------------------------------------------------*
      IF (NLRTHFLBL.GT.0) THEN
        INUM = IEFFFOCK('HAM0    ',1,1)
      END IF

      DO IVEC = 1, NLRTHFLBL
        INUM = IEFFFOCK(LRTHFLBL(IVEC),ISYLRTHF(IVEC),1)
      END DO

*=====================================================================*
* close lists:
      LR1OPN   = .FALSE.
      LL1OPN   = .FALSE.
      LO1OPN   = .FALSE.
      LX1OPN   = .FALSE.
      LR2OPN   = .FALSE.
      LX2OPN   = .FALSE.
      LL2OPN   = .FALSE.
      LO2OPN   = .FALSE.
      LR3OPN   = .FALSE.
      LX3OPN   = .FALSE.
      LL3OPN   = .FALSE.
      LO3OPN   = .FALSE.
      LR4OPN   = .FALSE.
      LX4OPN   = .FALSE.
      LL4OPN   = .FALSE.
      LO4OPN   = .FALSE.
      LN2OPN   = .FALSE.
      LER1OPN  = .FALSE.
      LER2OPN  = .FALSE.
      LEL1OPN  = .FALSE.
      LEL2OPN  = .FALSE.
      LRC1OPN  = .FALSE.
      LLC1OPN  = .FALSE.
      LCR2OPN  = .FALSE.
      LCO2OPN  = .FALSE.
      LCL2OPN  = .FALSE.
      LCX2OPN  = .FALSE.
      LEXPTOPN = .FALSE.
      LEFCKOPN = .FALSE.
      L1DXFOPN = .FALSE.
      LPL1OPN  = .FALSE.

* sort lists:
      CALL CCLSTSORT('O1 ',IDUM, IDUM, RDUM, ISYO1,LBLO1,FRQO1, IDUM,
     &                    LORXO1, ISYOFO1,  NO1LBL, MAXO1LBL, LDUM    )

      CALL CCLSTSORT('R1 ',IDUM, IDUM, RDUM, ISYLRT,LRTLBL,FRQLRT, IDUM,
     &                    LORXLRT, ISYOFT,  NLRTLBL, MAXTLBL, LDUM    )

      CALL CCLSTSORT('RC ',IDUM, IDUM, RDUM, ISYLRC,LRCLBL,RDUM,ILRCAU,
     &                    LDUM, ISYOFC,  NLRCLBL, MAXCLBL, LDUM    )

      CALL CCLSTSORT('X1 ',IDUM, IDUM, RDUM, ISYX1,LBLX1,FRQX1, IDUM,
     &                    LORXX1, ISYOFX1,  NX1LBL, MAXX1LBL, LDUM    )

      CALL CCLSTSORT('L1 ',IDUM, IDUM, RDUM, ISYLRZ,LRZLBL,FRQLRZ, IDUM,
     &                    LORXLRZ, ISYOFZ,  NLRZLBL, MAXZLBL, LDUM    )

      CALL CCLSTSORT('LC ',IDUM, IDUM, RDUM, ISYLC1,LBLLC1,RDUM,ILC1CAU,
     &                    LDUM, ISYOFLC1, NLC1LBL, MAXLC1LBL, LDUM  )

      CALL CCLSTSORT('O2 ',IDUM, IDUM, RDUM, ISYO2, LBLO2, FRQO2, IDUM,
     &                    LDUM, ISYOFO2, NO2LBL,  MAXO2LBL, LDUM   )

      CALL CCLSTSORT('CO2',IDUM,IDUM, RDUM, ISYCO2,LBLCO2,RDUM,ICO2CAU,
     &                    LDUM, ISYOFCO2,NCO2LBL, MAXCO2LBL,LDUM   )
 
      CALL CCLSTSORT('X2 ',IDUM, IDUM, RDUM, ISYX2, LBLX2, FRQX2, IDUM,
     &                    LDUM, ISYOFX2, NX2LBL,  MAXX2LBL, LDUM   )

      CALL CCLSTSORT('CX2',IDUM,IDUM, RDUM, ISYCX2,LBLCX2,RDUM,ICX2CAU,
     &                    LDUM, ISYOFCX2,NCX2LBL, MAXCX2LBL,LDUM   )
 
      CALL CCLSTSORT('R2 ',IDUM, IDUM, RDUM, ISYR2T,LBLR2T,FRQR2T,IDUM,
     &                    LDUM, ISYOFT2, NR2TLBL, MAXT2LBL, LDUM   )
 
      CALL CCLSTSORT('CR2',IDUM,IDUM, RDUM, ISYCR2,LBLCR2,RDUM,ICR2CAU,
     &                    LDUM, ISYOFCR2,NCR2LBL, MAXCR2LBL,LDUM   )
 
      CALL CCLSTSORT('L2 ',IDUM, IDUM, RDUM, ISYL2, LBLL2, FRQL2, IDUM,
     &                    LDUM, ISYOFL2, NL2LBL,  MAXL2LBL, LDUM   )

      CALL CCLSTSORT('CL2',IDUM,IDUM, RDUM, ISYCL2,LBLCL2,RDUM,ICL2CAU,
     &                    LDUM, ISYOFCL2,NCL2LBL, MAXCL2LBL,LDUM   )
 
      CALL CCLSTSORT('O3 ',IDUM, IDUM, RDUM, ISYO3, LBLO3, FRQO3, IDUM,
     &                    LDUM, ISYOFO3, NO3LBL,  MAXO3LBL, LDUM   )

      CALL CCLSTSORT('X3 ',IDUM, IDUM, RDUM, ISYX3, LBLX3, FRQX3, IDUM,
     &                    LDUM, ISYOFX3, NX3LBL,  MAXX3LBL, LDUM   )

      CALL CCLSTSORT('R3 ',IDUM, IDUM, RDUM, ISYR3T,LBLR3T,FRQR3T, IDUM,
     &                    LDUM, ISYOFT3, NR3TLBL, MAXT3LBL, LDUM   )
 
      CALL CCLSTSORT('L3 ',IDUM, IDUM, RDUM, ISYL3, LBLL3, FRQL3, IDUM,
     &                    LDUM, ISYOFL3, NL3LBL,  MAXL3LBL, LDUM   )

      CALL CCLSTSORT('O4 ',IDUM, IDUM, RDUM, ISYO4, LBLO4, FRQO4, IDUM,
     &                    LDUM, ISYOFO4, NO4LBL,  MAXO4LBL, LDUM   )

      CALL CCLSTSORT('X4 ',IDUM, IDUM, RDUM, ISYX4, LBLX4, FRQX4, IDUM,
     &                    LDUM, ISYOFX4, NX4LBL,  MAXX4LBL, LDUM   )

      CALL CCLSTSORT('R4 ',IDUM, IDUM, RDUM, ISYR4T,LBLR4T,FRQR4T, IDUM,
     &                    LDUM, ISYOFT4, NR4TLBL, MAXT4LBL, LDUM   )
 
      CALL CCLSTSORT('L4 ',IDUM, IDUM, RDUM, ISYL4, LBLL4, FRQL4, IDUM,
     &                    LDUM, ISYOFL4, NL4LBL,  MAXL4LBL, LDUM   )

      CALL CCLSTSORT('M1 ',ISYLRM, ILRM, FRQLRM, IDUM, CDUM, RDUM, IDUM,
     &                    LDUM, ISYOFM,  NLRM,    MAXM,     LDUM   )

      CALL CCLSTSORT('N2 ',ISYSN2, ISTN2, EIGN2, IDUM, CDUM, RDUM, IDUM,
     &                    LDUM, ISYOFN2, NQRN2,   MAXQRN2,  LDUM   )

      CALL CCLSTSORT('ER1',ISYSER1, ISTER1, EIGER1, 
     &                     ISYOER1, LBLER1, FRQER1, IDUM, LDUM, 
     &                     ISYOFER1, NER1LBL, MAXER1LBL, LPRER1 )

      CALL CCLSTSORT('ER2',ISYSER2, ISTER2, EIGER2, 
     &                     ISYOER2, LBLER2, FRQER2, IDUM, LDUM, 
     &                     ISYOFER2, NER2LBL, MAXER2LBL, LPRER2 )

      CALL CCLSTSORT('EL1',ISYSEL1, ISTEL1, EIGEL1, 
     &                     ISYOEL1, LBLEL1, FRQEL1, IDUM, LORXEL1,
     &                     ISYOFEL1, NEL1LBL, MAXEL1LBL, LPREL1 )

      CALL CCLSTSORT('EL2',ISYSEL2, ISTEL2, EIGEL2, 
     &                     ISYOEL2, LBLEL2, FRQEL2, IDUM, LDUM,
     &                     ISYOFEL2, NEL2LBL, MAXEL2LBL, LPREL2 )

      CALL CCLSTSORT('PL1',ISYSPL1, ISTPL1, EIGPL1,
     &                     ISYPL1, LBLPL1, FRQPL1, IDUM, LORXPL1,
     &                     ISYOFPL1, NPL1LBL, MAXPL1LBL, LPRPL1)

* print sorted lists to output:
       IF (NRSOLBL.GT.0) THEN
         CALL AROUND('REQUESTED PROPERTY OPERATORS:')
         WRITE(LUPRI,'(13X,A,/,13X,50("-"))')
     &     'Index   Oper. Label  Symmetry  Transp.  PDBS  Atom'
         DO I = 1, NRSOLBL
           WRITE(LUPRI,'(12X,I5,5X,A8,4X,I5,4X,I5,4X,L3,2X,I5)') 
     &      I, LBLOPR(I), ISYOPR(I), ISYMAT(I), LPDBSOP(I),IATOPR(I)
         END DO
         WRITE(LUPRI,'(13X,50("-"),//)')
       END IF

       IF (NEXPFCKLBL.GT.0) THEN
         CALL AROUND('REQUESTED EXPECTATION VALUES:')
         WRITE(LUPRI,'(23X,A,/,23X,29("-"))')
     &     'Index   Oper. Label  Symmetry'
         DO I = 1, NEXPFCKLBL
            IF (LEXPFCK(1,I)) WRITE(LUPRI,'(22X,I5,5X,A8,4X,I5)')
     &             I, LBLEXPFCK(I), ISYEXPFCK(I)
         END DO
         WRITE(LUPRI,'(23X,29("-"),//)')
       END IF

       IF (NEXPFCKLBL.GT.0) THEN
         CALL AROUND('REQUESTED EFFECTIVE FOCK MATRICES:')
         WRITE(LUPRI,'(23X,A,/,23X,29("-"))')
     &     'Index   Oper. Label  Symmetry'
         DO I = 1, NEXPFCKLBL
            IF (LEXPFCK(2,I)) WRITE(LUPRI,'(22X,I5,5X,A8,4X,I5)')
     &             I, LBLEXPFCK(I), ISYEXPFCK(I)
         END DO
         WRITE(LUPRI,'(23X,29("-"),//)')
       END IF

       IF (N1DXFLBL.GT.0) THEN
         CALL AROUND('REQUESTED 1-IDX-TRAN EFF. FOCK M.:')
         WRITE(LUPRI,'(22X,A,/,22X,32("-"))')
     &     'Index   Oper. Label  Type  Index'
         DO I = 1, N1DXFLBL
            WRITE(LUPRI,'(21X,I5,5X,A8,6X,A3,I5)')
     &       I, LBL1DXFCK(I), LST1DXFCK(I), IRELAX1DX(I)
         END DO
         WRITE(LUPRI,'(22X,32("-"),//)')
       END IF

       IF (NLRTHFLBL.GT.0) THEN
         CALL AROUND('REQUESTED FIRST ORDER KAPPA VECTORS:')
         WRITE(LUPRI,'(18X,A,/,18X,41("-"))')
     &     'Index   Oper. Label   Sym.    Frequency'
         DO I = 1, NLRTHFLBL
           WRITE(LUPRI,'(18X,I4,6X,A8,I6,2X,1P,D15.6)') 
     &            I, LRTHFLBL(I), ISYLRTHF(I), FRQLRTHF(I)
         END DO
         WRITE(LUPRI,'(18X,41("-"),//)')
       END IF

       IF (NO1LBL.GT.0) THEN
         CALL AROUND('REQUESTED FIRST ORDER XI VECTORS:')
         WRITE(LUPRI,'(13X,A,/,13X,50("-"))')
     &     'Index   Oper. Label  relaxed  Sym.     Frequency'
         DO I = 1, NO1LBL
           WRITE(LUPRI,'(13X,I4,6X,A8,5X,L3,I6,3X,1P,D15.6)') 
     &            I, LBLO1(I), LORXO1(I), ISYO1(I), FRQO1(I)
Cholesky
           NCHOPLR(ISYO1(I)) = NCHOPLR(ISYO1(I)) + 1
Cholesky
         END DO
         WRITE(LUPRI,'(13X,50("-"),//)')
       END IF

       IF (NLRTLBL.GT.0) THEN
         CALL AROUND('REQUESTED FIRST ORDER T VECTORS:')
         WRITE(LUPRI,'(13X,A,/,13X,50("-"))')
     &     'Index   Oper. Label  relaxed  Sym.     Frequency'
         DO I = 1, NLRTLBL
           WRITE(LUPRI,'(13X,I4,6X,A8,5X,L3,I6,3X,1P,D15.6)') 
     &            I, LRTLBL(I), LORXLRT(I), ISYLRT(I), FRQLRT(I)
         END DO
         WRITE(LUPRI,'(13X,50("-"),//)')
       END IF

       IF (NX1LBL.GT.0) THEN
         CALL AROUND('REQUESTED FIRST ORDER ETA VECTORS:')
         WRITE(LUPRI,'(13X,A,/,13X,50("-"))')
     &     'Index   Oper. Label  relaxed  Sym.     Frequency'
         DO I = 1, NX1LBL
           WRITE(LUPRI,'(13X,I4,6X,A8,5X,L3,I6,3X,1P,D15.6)') 
     &            I, LBLX1(I), LORXX1(I), ISYX1(I), FRQX1(I)
         END DO
         WRITE(LUPRI,'(13X,50("-"),//)')
       END IF

       IF (NLRZLBL.GT.0) THEN
         CALL AROUND('REQUESTED FIRST ORDER ZETA VECTORS:')
         WRITE(LUPRI,'(13X,A,/,13X,50("-"))')
     &     'Index   Oper. Label  relaxed  Sym.     Frequency'
         DO I = 1, NLRZLBL
           WRITE(LUPRI,'(13X,I4,6X,A8,5X,L3,I6,3X,1P,D15.6)') 
     &            I, LRZLBL(I), LORXLRZ(I), ISYLRZ(I), FRQLRZ(I)
         END DO
         WRITE(LUPRI,'(13X,50("-"),//)')
       END IF

       IF (NLRM.GT.0) THEN
         CALL AROUND('REQUESTED FIRST ORDER M-VECTORS:')
         WRITE(LUPRI,'(15X,A,/,15X,50("-"))')
     &     'Index   State Symmetry  Frequency'
         DO I = 1, NLRM   
            WRITE(LUPRI,'(14X,I5,6X,I8,I5,2X,1P,D15.6)') 
     &            I, ILRM(I), ISYLRM(I), FRQLRM(I)
         END DO
         WRITE(LUPRI,'(15X,50("-"),//)')
       END IF

       IF (NLRCLBL.GT.0) THEN
         CALL AROUND(
     &    'REQUESTED FIRST-ORDER RIGHT CAUCHY VECTORS:')
         WRITE(LUPRI,'(15X,A,/,15X,50("-"))')
     &     'Index   Oper. Label  Symmetry  Cauchy Order'
         DO I = 1, NLRCLBL
            WRITE(LUPRI,'(14X,I5,6X,A8,I5,2X,1P,I5)') 
     &            I, LRCLBL(I), ISYLRC(I), ILRCAU(I)
         END DO
         WRITE(LUPRI,'(15X,50("-"),//)')
       END IF

       IF (NLC1LBL.GT.0) THEN
         CALL AROUND(
     &    'REQUESTED FIRST-ORDER LEFT CAUCHY VECTORS:')
         WRITE(LUPRI,'(15X,A,/,15X,50("-"))')
     &     'Index   Oper. Label  Symmetry  Cauchy Order'
         DO I = 1, NLC1LBL
            WRITE(LUPRI,'(14X,I5,6X,A8,I5,2X,1P,I5)') 
     &            I, LBLLC1(I), ISYLC1(I), ILC1CAU(I)
         END DO
         WRITE(LUPRI,'(15X,50("-"),//)')
       END IF

       IF (NO2LBL.GT.0) THEN
         CALL AROUND('REQUESTED SECOND-ORDER XKSI VECTORS:')
         WRITE(LUPRI,'(5X,2A,/,5X,70("-"))')
     &     'Index    Oper. Label  Symmetry  Frequency',
     &           '   Oper. Label  Symmetry  Frequency  Symmetry'
         DO I = 1, NO2LBL
           WRITE(LUPRI,'(I5,2(3X,A8,I3,2X,1P,D15.6),I3)') 
     &            I, LBLAO2(I), ISYAO2(I), FRQAO2(I),
     &               LBLBO2(I), ISYBO2(I), FRQBO2(I),
     &            MULD2H(ISYAO2(I),ISYBO2(I))
         END DO
         WRITE(LUPRI,'(5X,70("-"),//)')
       END IF 

       IF (NX2LBL.GT.0) THEN
         CALL AROUND('REQUESTED SECOND-ORDER ETA VECTORS:')
         WRITE(LUPRI,'(5X,2A,/,5X,70("-"))')
     &     'Index    Oper. Label  Symmetry  Frequency',
     &           '   Oper. Label  Symmetry  Frequency  Symmetry'
         DO I = 1, NX2LBL
           WRITE(LUPRI,'(I5,2(3X,A8,I3,2X,1P,D15.6),I3)') 
     &            I, LBLAX2(I), ISYAX2(I), FRQAX2(I),
     &               LBLBX2(I), ISYBX2(I), FRQBX2(I),
     &            MULD2H(ISYAX2(I),ISYBX2(I))
         END DO
         WRITE(LUPRI,'(5X,70("-"),//)')
       END IF

       IF (NR2TLBL.GT.0) THEN
         CALL AROUND('REQUESTED SECOND-ORDER T VECTORS:')
         WRITE(LUPRI,'(5X,2A,/,5X,70("-"))')
     &     'Index    Oper. Label  Symmetry  Frequency',
     &           '   Oper. Label  Symmetry  Frequency  Symmetry'
         DO I = 1, NR2TLBL
           WRITE(LUPRI,'(I5,2(3X,A8,I3,2X,1P,D15.6),I3)') 
     &            I, LBLAR2T(I), ISYAR2T(I), FRQAR2T(I),
     &               LBLBR2T(I), ISYBR2T(I), FRQBR2T(I),
     &            MULD2H(ISYAR2T(I),ISYBR2T(I))
         END DO
         WRITE(LUPRI,'(5X,70("-"),//)')
       END IF

       IF (NL2LBL.GT.0) THEN
         CALL AROUND('REQUESTED SECOND-ORDER ZETA VECTORS:')
         WRITE(LUPRI,'(5X,2A,/,5X,70("-"))')
     &     'Index    Oper. Label  Symmetry  Frequency',
     &           '   Oper. Label  Symmetry  Frequency  Symmetry'
         DO I = 1, NL2LBL
           WRITE(LUPRI,'(I5,2(3X,A8,I3,2X,1P,D15.6),I3)') 
     &            I, LBLAL2(I), ISYAL2(I), FRQAL2(I),
     &               LBLBL2(I), ISYBL2(I), FRQBL2(I),
     &            MULD2H(ISYAL2(I),ISYBL2(I))
         END DO
         WRITE(LUPRI,'(5X,70("-"),//)')
       END IF

       IF (NCR2LBL.GT.0) THEN
         CALL AROUND('REQUESTED SECOND-ORDER RIGHT'//
     &                    ' CAUCHY VECTORS:')
         WRITE(LUPRI,'(5X,3A,/,5X,70("-"))')
     &     'Index   Oper. Label  Symmetry  Cauchy Order',
     &          '   Oper. Label  Symmetry  Cauchy Order',
     &                        '  Symmetry  Cauchy Order'
         DO I = 1, NCR2LBL
           WRITE(LUPRI,'(I5,2(3X,A8,I3,2X,1P,I5),I3)') 
     &            I, LBLCR2(I,1), ISYCR2(I,1), ICR2CAU(I,1),
     &               LBLCR2(I,2), ISYCR2(I,2), ICR2CAU(I,2),
     &            MULD2H(ISYCR2(I,1),ISYCR2(I,2))
         END DO
         WRITE(LUPRI,'(5X,70("-"),//)')
       END IF

       IF (NCO2LBL.GT.0) THEN
         CALL AROUND('REQUESTED SECOND-ORDER CAUCHY'//
     &                    ' XKSI VECTORS:')
         WRITE(LUPRI,'(5X,3A,/,5X,70("-"))')
     &     'Index   Oper. Label  Symmetry  Cauchy Order',
     &          '   Oper. Label  Symmetry  Cauchy Order',
     &                        '  Symmetry  Cauchy Order'
         DO I = 1, NCO2LBL
           WRITE(LUPRI,'(I5,2(3X,A8,I3,2X,1P,I5),I3)') 
     &            I, LBLCO2(I,1), ISYCO2(I,1), ICO2CAU(I,1),
     &               LBLCO2(I,2), ISYCO2(I,2), ICO2CAU(I,2),
     &            MULD2H(ISYCO2(I,1),ISYCO2(I,2))
         END DO
         WRITE(LUPRI,'(5X,70("-"),//)')
       END IF

       IF (NCL2LBL.GT.0) THEN
         CALL AROUND('REQUESTED SECOND-ORDER LEFT'//
     &                    ' CAUCHY VECTORS:')
         WRITE(LUPRI,'(5X,3A,/,5X,70("-"))')
     &     'Index   Oper. Label  Symmetry  Cauchy Order',
     &          '   Oper. Label  Symmetry  Cauchy Order',
     &                        '  Symmetry  Cauchy Order'
         DO I = 1, NCL2LBL
           WRITE(LUPRI,'(I5,2(3X,A8,I3,2X,1P,I5),I3)') 
     &            I, LBLCL2(I,1), ISYCL2(I,1), ICL2CAU(I,1),
     &               LBLCL2(I,2), ISYCL2(I,2), ICL2CAU(I,2),
     &            MULD2H(ISYCL2(I,1),ISYCL2(I,2))
         END DO
         WRITE(LUPRI,'(5X,70("-"),//)')
       END IF

       IF (NCX2LBL.GT.0) THEN
         CALL AROUND('REQUESTED SECOND-ORDER CAUCHY'//
     &                    ' ETA VECTORS:')
         WRITE(LUPRI,'(5X,3A,/,5X,70("-"))')
     &     'Index   Oper. Label  Symmetry  Cauchy Order',
     &          '   Oper. Label  Symmetry  Cauchy Order',
     &                        '  Symmetry  Cauchy Order'
         DO I = 1, NCX2LBL
           WRITE(LUPRI,'(I5,2(3X,A8,I3,2X,1P,I5),I3)') 
     &            I, LBLCX2(I,1), ISYCX2(I,1), ICX2CAU(I,1),
     &               LBLCX2(I,2), ISYCX2(I,2), ICX2CAU(I,2),
     &            MULD2H(ISYCX2(I,1),ISYCX2(I,2))
         END DO
         WRITE(LUPRI,'(5X,70("-"),//)')
       END IF

       IF (NO3LBL.GT.0) THEN
         CALL AROUND('REQUESTED THIRD-ORDER XKSI VECTORS:')
         WRITE(LUPRI,'(5X,3A,/,5X,70("-"))')
     &     'Index    Oper. Label  Symmetry  Frequency',
     &           '   Oper. Label  Symmetry  Frequency',
     &           '   Oper. Label  Symmetry  Frequency  Symmetry'
         DO I = 1, NO3LBL
           WRITE(LUPRI,'(I5,3(3X,A8,I3,2X,1P,D15.6),I3)') 
     &            I, (LBLO3(I,J), ISYO3(I,J), FRQO3(I,J), J=1, 3),
     &            ILSTSYM('O3',I)
         END DO
         WRITE(LUPRI,'(5X,70("-"),//)')
       END IF

       IF (NR3TLBL.GT.0) THEN
         CALL AROUND('REQUESTED THIRD-ORDER T VECTORS:')
         WRITE(LUPRI,'(5X,3A,/,5X,70("-"))')
     &     'Index    Oper. Label  Symmetry  Frequency',
     &           '   Oper. Label  Symmetry  Frequency',
     &           '   Oper. Label  Symmetry  Frequency  Symmetry'
         DO I = 1, NR3TLBL
           WRITE(LUPRI,'(I5,3(3X,A8,I3,2X,1P,D15.6),I3)') 
     &            I, (LBLR3T(I,J), ISYR3T(I,J), FRQR3T(I,J), J=1, 3),
     &            ILSTSYM('R3',I)
         END DO
         WRITE(LUPRI,'(5X,70("-"),//)')
       END IF

       IF (NX3LBL.GT.0) THEN
         CALL AROUND('REQUESTED THIRD-ORDER ETA VECTORS:')
         WRITE(LUPRI,'(5X,3A,/,5X,70("-"))')
     &     'Index    Oper. Label  Symmetry  Frequency',
     &           '   Oper. Label  Symmetry  Frequency',
     &           '   Oper. Label  Symmetry  Frequency  Symmetry'
         DO I = 1, NX3LBL
           WRITE(LUPRI,'(I5,3(3X,A8,I3,2X,1P,D15.6),I3)') 
     &            I, (LBLX3(I,J), ISYX3(I,J), FRQX3(I,J), J=1, 3),
     &            ILSTSYM('X3',I)
         END DO
         WRITE(LUPRI,'(5X,70("-"),//)')
       END IF

       IF (NL3LBL.GT.0) THEN
         CALL AROUND('REQUESTED THIRD-ORDER ZETA VECTORS:')
         WRITE(LUPRI,'(5X,3A,/,5X,70("-"))')
     &     'Index    Oper. Label  Symmetry  Frequency',
     &           '   Oper. Label  Symmetry  Frequency',
     &           '   Oper. Label  Symmetry  Frequency  Symmetry'
         DO I = 1, NL3LBL
           WRITE(LUPRI,'(I5,3(3X,A8,I3,2X,1P,D15.6),I3)') 
     &            I, (LBLL3(I,J), ISYL3(I,J), FRQL3(I,J), J=1, 3),
     &            ILSTSYM('L3',I)
         END DO
         WRITE(LUPRI,'(5X,70("-"),//)')
       END IF
       
       IF (NER1LBL.GT.0) THEN
         CALL AROUND('REQUESTED FIRST-ORDER RIGHT'//
     &                 ' EIGENVECTOR RESPONSES:')
         WRITE (LUPRI,'(3x,69("-"),/3x,2A,/3x,69("-"))')
     &      'IDX STATE/SYM   EXC. ENERGY    OPE',
     &                  'RATOR/SYM    FREQUENCY    SYM R  P'
         DO I = 1, NER1LBL
          WRITE (LUPRI,'(I5,3X,I3,I3,2X,1P,D15.6,3X,A8,I3,2X,1P,
     &           D15.6,I3,L3)')
     &          I, ISTER1(I), ISYSER1(I), EIGER1(I),
     &             LBLER1(I), ISYOER1(I), FRQER1(I), ILSTSYM('ER1',I),
     &             LPRER1(I)
         END DO
         WRITE(LUPRI,'(3x,69("-"),//)')
       END IF
       
       IF (NER2LBL.GT.0) THEN
         CALL AROUND('REQUESTED SECOND-ORDER RIGHT'//
     &               ' EIGENVECTOR RESPONSES:')
         WRITE(LUPRI,'(3x,69("-"),/)')
         DO I = 1, NER2LBL
          WRITE (LUPRI,'(I5,3X,2I3,2X,1P,D15.6,2(3X,A8,I3,2X,1P,
     &           D15.6),I3,L3)')
     &         I,ISTER2(I),  ISYSER2(I),  EIGER2(I),
     &           LBLER2(I,1),ISYOER2(I,1),FRQER2(I,1), 
     &           LBLER2(I,2),ISYOER2(I,2),FRQER2(I,2),ILSTSYM('ER2',I),
     &           LPRER2(I)
         END DO
         WRITE(LUPRI,'(3x,69("-"),//)')
       END IF

       IF (NEL1LBL.GT.0) THEN
         CALL AROUND('REQUESTED FIRST-ORDER LEFT'//
     &                ' EIGENVECTOR RESPONSES:')
         WRITE (LUPRI,'(3x,69("-"),/3x,2A,/3x,69("-"))')
     &     'IDX STATE/SYM   EXC. ENERGY    OPE',
     &                  'RATOR/SYM    FREQUENCY    SYM R  P'
         DO I = 1, NEL1LBL
          WRITE (LUPRI,'(I5,3X,I3,I3,2X,1P,D15.6,3X,A8,I3,2X,1P,
     &           D15.6,I3,2L3)')
     &          I, ISTEL1(I), ISYSEL1(I), EIGEL1(I),
     &             LBLEL1(I), ISYOEL1(I), FRQEL1(I), ILSTSYM('EL1',I),
     &             LORXEL1(I),LPREL1(I)
         END DO
         WRITE(LUPRI,'(3x,69("-"),//)')
       END IF

       IF (NEL2LBL.GT.0) THEN
         CALL AROUND('REQUESTED SECOND-ORDER LEFT'//
     &                ' EIGENVECTOR RESPONSES:')
         WRITE(LUPRI,'(3x,69("-"))')
         DO I = 1, NEL2LBL
          WRITE (LUPRI,'(I5,3X,2I3,2X,1P,D15.6,2(3X,A8,I3,2X,1P,
     &           D15.6),I3,L3)')
     &         I,ISTEL2(I),  ISYSEL2(I),  EIGEL2(I),
     &           LBLEL2(I,1),ISYOEL2(I,1),FRQEL2(I,1), 
     &           LBLEL2(I,2),ISYOEL2(I,2),FRQEL2(I,2),ILSTSYM('EL2',I),
     &           LPREL2(I)
         END DO
         WRITE(LUPRI,'(3x,69("-"),//)')
       END IF
       
       IF (NXGRST.GT.0) THEN
         WRITE(LUPRI,'(/A)')' LIST OF REQUIRED ZEROTH-ORDER E0 MULTIP.:'
         DO I = 1, NXGRST   
            WRITE(LUPRI,'(I5,3X,I8,I5,2X,1P,D15.6)') 
     &            I, IXGRST(I), ISYEXC(IXGRST(I)), EIGVAL(IXGRST(I))
         END DO
       END IF

       IF (NQRN2.GT.0) THEN
         CALL AROUND('REQUESTED N(i,f) VECTORS:')
         WRITE(LUPRI,'(/3x,69("-"),/3x,2A,/3x,69("-"))')
     &    'IDX STATE/SYM    EXC. ENERGY',
     &      '  STATE/SYM    EXC. ENERGY SYM'
         DO I = 1, NQRN2
           WRITE(LUPRI,'(I5,2(3X,I3,I3,2X,1P,D15.6),I3)')
     &            I, IIN2(I), ISYIN2(I), FRQIN2(I),
     &               IFN2(I), ISYFN2(I), FRQFN2(I),
     &            MULD2H(ISYIN2(I),ISYFN2(I))
         END DO
         WRITE(LUPRI,'(3x,69("-"),//)')
       END IF

       IF (NPL1LBL.GT.0) THEN
         CALL AROUND('REQUESTED PROJECTED FIRST-ORDER ZETA VECTORS:')
         WRITE(LUPRI,'(3x,69("-"))')
         DO I = 1, NPL1LBL
          WRITE (LUPRI,'(I5,3X,I3,I3,2X,1P,D15.6,3X,A8,I3,2X,1P,D15.6,
     &           I3,L3)')
     &          I, ISTPL1(I), ISYSPL1(I), EIGPL1(I),
     &             LBLPL1(I), ISYPL1(I), FRQPL1(I), ISYPL1(I),
     &             LPRPL1(I)
         END DO
         WRITE(LUPRI,'(3x,69("-"),//)')
       END IF

      RETURN
      END
*---------------------------------------------------------------------*
c /* deck cc_exgrind */
*=====================================================================*
       SUBROUTINE CC_EXGRIND
*---------------------------------------------------------------------*
*
*    Purpose: Control input and equations for calculation of 
*             excited state first order properties.
*
*    OC April 1997
*
*=====================================================================*
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <cclrinf.h>
#include <ccrspprp.h>
#include <ccexci.h>  
#include <ccexgr.h>  
#include <ccsdinp.h>  
#include <ccsdsym.h>  
#include <cclr.h>  
#include <ccroper.h>  
      CHARACTER*8 LABEL,LABELA
      INTEGER ISYMA, INUM, IOP, ISYME

      REAL*8  EIGV
C
      LOGICAL FIRSTCALL
      SAVE    FIRSTCALL
      DATA    FIRSTCALL /.TRUE./
C
*------------------------------------------------------------------
* test if operators are available and translate IAXGRO, 
* arrays from the PRPLBL list to the new list maintained by IROPER.
*------------------------------------------------------------------
      IF (FIRSTCALL) THEN
       IOPER = 1
       DO WHILE (IOPER .LE. NAXGRO)
        LABELA = PRPLBL(IAXGRO(IOPER))
        IF (DEBUG) THEN
          WRITE(LUPRI,'(/2X,A,1X,A)')
     &     'CHECK EXGR OPERATOR:',LABELA
        ENDIF
        IF (IROPER(LABELA,ISYMA) .LT. 0) THEN
          WRITE(LUPRI,'(/2X,3A,/2X,2A)')
     &     ' WARNING: THE OPERATOR WITH THE LABEL "',
     &     LABELA,'" IS NOT AVAILABLE.',
     &     ' LINEAR RESPONSE RESIDUE CALCULATION IS CANCELED FOR THIS',
     &     ' OPERATOR.'
          DO IDX = IOPER, NAXGRO-1
            IAXGRO(IDX) = IAXGRO(IDX+1)
          END DO
          NAXGRO = NAXGRO - 1
        ELSE
          IF (DEBUG) THEN
             WRITE(LUPRI,'(/2X,A,1X,A,A)')
     &        'PUT :',LABELA,' ON THE LIST.'
          ENDIF
          IAXGRO(IOPER) = IROPER(LABELA,ISYMA)
          IOPER = IOPER + 1
        END IF
       END DO
       FIRSTCALL = .FALSE.
      END IF ! (FIRSTCALL)
C
      IF (DEBUG) THEN
         WRITE(LUPRI,'(/,A)') ' Updated list'
         DO IOPER = 1, NAXGRO
            WRITE(LUPRI,*) IOPER,IAXGRO(IOPER),' ',
     *                     LBLOPR(IAXGRO(IOPER)),
     *                     ISYOPR(IAXGRO(IOPER))
         ENDDO
      ENDIF
C
C------------------------------------
C     Fill in equations to be solved.
C------------------------------------
C
      NXGRST = 0
C
      DO 100 ISYME = 1, NSYM
       DO 200 IEX = 1, NCCEXCI(ISYME,1)
C
        IF (SELXGR) THEN
C
C        Check state has been calculated. 
C
           DO 300 I = 1,NSEXGR
             IF ((ISYME.EQ.ISEXGR(I,1))
     *           .AND.(IEX.EQ.ISEXGR(I,2))) THEN
                NXGRST = NXGRST + 1
                IXGRST(NXGRST) = ISYOFE(ISYME)+IEX 
                GO TO 350 
             ENDIF 
  300     CONTINUE
C 
C------------------------------------------------------------------------
C         This state is not requested in oscillator strength calculation.
C         GOTO end of loop.
C------------------------------------------------------------------------
C
          GO TO 200
        ELSE
          NXGRST = NXGRST + 1
          IF (NXGRST.GT.MXXGST)  THEN
            WRITE(LUPRI,*) 'NXGRST =',NXGRST,'MXXGST= ',MXXGST
            CALL QUIT(' Too many states for residue calculation')
          ENDIF
          IXGRST(NXGRST) = ISYOFE(ISYME)+IEX
        ENDIF

  350   CONTINUE
  200  CONTINUE
  100 CONTINUE
C
      ISYOFXG(1) = 0
      DO ISYM = 2, NSYM
        ISYOFXG(ISYM) = NXGRST
      END DO
C
      END
*---------------------------------------------------------------------*
c /* deck cc_opaind */
*=====================================================================*
       SUBROUTINE CC_OPAIND
*---------------------------------------------------------------------*
*
*    Purpose: Control input and equations for calculation of 
*             one-photon absorption strengths for ground to
*             excited state transitions
*
*    Christof Haettig, December 2002, Friedrichstal
*
*=====================================================================*
      IMPLICIT NONE
C#include <implicit.h>
C#include <cclrinf.h>
C#include <cclr.h>  
C#include <ccsdsym.h>  
#include <priunit.h>
#include <ccorb.h>
#include <ccexci.h>  
#include <ccexcinf.h>  
#include <ccrspprp.h>
#include <ccopainf.h>  
#include <ccsdinp.h>  
#include <ccroper.h>  

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*8 LABEL
      LOGICAL TAKE_STATE
      INTEGER ISYMO, INUM, ISTATE, ISYME, IOPER, IDX, IEX

      REAL*8  EIGV
* functions:
      INTEGER IR1TAMP,ILRMAMP,IROPER
C
      LOGICAL FIRSTCALL
      SAVE    FIRSTCALL
      DATA    FIRSTCALL /.TRUE./

*------------------------------------------------------------------
* test if operators are available and translate ILRSOP array
* from the PRPLBL list to the new list maintained by IROPER:
*------------------------------------------------------------------
      IF (FIRSTCALL) THEN
       IOPER = 1
       DO WHILE (IOPER .LE. NLRSOP)
        LABEL = PRPLBL(ILRSOP(IOPER))
        IF (DEBUG) WRITE(LUPRI,'(/2X,2A)') 'CHECK OPERATOR: ',LABEL
        IF (IROPER(LABEL,ISYMO) .LT. 0) THEN
          WRITE(LUPRI,'(/2X,3A,/2X,2A)')
     &     ' WARNING: THE OPERATOR WITH THE LABEL "',
     &         LABEL,'" IS NOT AVAILABLE.',
     &     ' CALCULATION OF TRANSITION MOMENTS IS CANCELED FOR THIS',
     &     ' OPERATOR.'
          DO IDX = IOPER, NLRSOP-1
            ILRSOP(IDX) = ILRSOP(IDX+1)
          END DO
          NLRSOP = NLRSOP - 1
        ELSE
          IF(DEBUG)WRITE(LUPRI,'(/2X,3A)')'PUT: ',LABEL,' ON THE LIST.'
          ILRSOP(IOPER) = IROPER(LABEL,ISYMO)
          IOPER = IOPER + 1
        END IF
       END DO
       FIRSTCALL = .FALSE.
      END IF ! (FIRSTCALL)
C
      IF (LOCDBG) THEN
         WRITE(LUPRI,'(/,A)') ' Updated list in CC_OPAIND:'
         DO IOPER = 1, NLRSOP
            WRITE(LUPRI,*) IOPER,ILRSOP(IOPER),' ',
     *              LBLOPR(ILRSOP(IOPER)),ISYOPR(ILRSOP(IOPER))
         ENDDO
      ENDIF
C
C------------------------------------
C     Fill in equations to be solved.
C------------------------------------
C
      NXLRSST = 0
C
      DO ISYME = 1, NSYM
       DO IEX = 1, NCCEXCI(ISYME,1)
        ISTATE = ISYOFE(ISYME)+IEX
        EIGV   = EIGVAL(ISTATE)
 
        IF (SELLRS) THEN
          ! check, if state has been requested
          TAKE_STATE = .FALSE.
          DO IDX = 1,NSELRS
            IF (ISYME.EQ.ISELRSYM(IDX) .AND. IEX.EQ.ISELRSTA(IDX)) THEN
              TAKE_STATE = .TRUE.
            ENDIF 
          END DO
        ELSE
          TAKE_STATE = .TRUE.
        ENDIF

        IF (TAKE_STATE) THEN
          NXLRSST = NXLRSST + 1
          IF (NXLRSST.GT.MXLRSST)  THEN
            WRITE(LUPRI,*) 'NXLRSST =',NXLRSST,'MXLRSST= ',MXLRSST
            CALL QUIT(' Too many states in CC_OPAIND')
          ENDIF
          ILRSST(NXLRSST) = ISTATE

          IF (.NOT.CIS) THEN
            DO IOPER = 1, NLRSOP
              LABEL = LBLOPR(ILRSOP(IOPER))
              ISYMO = ISYOPR(ILRSOP(IOPER))
              IF (ISYME.EQ.ISYMO) THEN
                IF (.NOT.LRS2N1) THEN
                  INUM = IR1TAMP(LABEL,.FALSE.,-EIGV,ISYMO)
                ELSE
                  INUM = ILRMAMP(ISTATE,EIGV,ISYME)
                ENDIF
              ENDIF
            END DO
          ENDIF 
        END IF ! (TAKE_STATE) 

       END DO
      END DO
C
      RETURN
      END
*---------------------------------------------------------------------*
c /* deck cc_xopaind */
*=====================================================================*
       SUBROUTINE CC_XOPAIND
*---------------------------------------------------------------------*
*
*    Purpose: Control input and equations for calculation of 
*             one-photon absorption strengths for excited to
*             to excited state transitions
*
*    Christof Haettig, October 2003, Friedrichstal
*
*=====================================================================*
      IMPLICIT NONE
#include <priunit.h>
#include <ccorb.h>
#include <ccexci.h>  
#include <ccexcinf.h>  
#include <ccrspprp.h>
#include <ccxopainf.h>  
#include <ccsdinp.h>  
#include <ccroper.h>  

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*8 LABEL
      LOGICAL TAKE_STATE_PAIR
      INTEGER ISYMO, INUM, ISYMI, ISYMF, ISTATEI, ISTATEF, IEXI, IEXF,
     &        IDX, ISYMFI, IOPER, I, NSYMF, NEXF

      REAL*8  EIGVI, EIGVF
* functions:
      INTEGER IROPER, IR1TAMP, IN2AMP

      LOGICAL FIRSTCALL
      SAVE    FIRSTCALL
      DATA    FIRSTCALL /.TRUE./

*------------------------------------------------------------------
* test if operators are available and translate IQR2OP array
* from the PRPLBL list to the new list maintained by IROPER:
*------------------------------------------------------------------
      IF (FIRSTCALL) THEN
       IOPER = 1
       DO WHILE (IOPER .LE. NQR2OP)
        LABEL = PRPLBL(IQR2OP(IOPER))
        IF (DEBUG) WRITE(LUPRI,'(/2X,2A)') 'CHECK OPERATOR: ',LABEL
        IF (IROPER(LABEL,ISYMO) .LT. 0) THEN
          WRITE(LUPRI,'(/2X,3A,/2X,2A)')
     &     ' WARNING: THE OPERATOR WITH THE LABEL "',
     &         LABEL,'" IS NOT AVAILABLE.',
     &     ' CALCULATION OF TRANSITION MOMENTS IS CANCELED FOR THIS',
     &     ' OPERATOR.'
          DO IDX = IOPER, NQR2OP-1
            IQR2OP(IDX) = IQR2OP(IDX+1)
          END DO
          NQR2OP = NQR2OP - 1
        ELSE
          IF(DEBUG)WRITE(LUPRI,'(/2X,3A)')'PUT: ',LABEL,' ON THE LIST.'
          IQR2OP(IOPER) = IROPER(LABEL,ISYMO)
          IOPER = IOPER + 1
        END IF
       END DO
       FIRSTCALL = .FALSE.
      END IF ! (FIRSTCALL)
C
      IF (LOCDBG) THEN
         WRITE(LUPRI,'(/,A)') ' Updated list in CC_XOPAIND:'
         DO IOPER = 1, NQR2OP
            WRITE(LUPRI,*) IOPER,IQR2OP(IOPER),' ',
     *              LBLOPR(IQR2OP(IOPER)),ISYOPR(IQR2OP(IOPER))
         ENDDO
      ENDIF
C
C------------------------------------
C     Fill in equations to be solved.
C------------------------------------
C
      IF (LOCDBG) WRITE(LUPRI,*) 'SELQR2:',SELQR2
C
      NXQR2ST = 0
C
      DO ISYMI = 1, NSYM

       NSYMF = NSYM
       IF (.NOT.SELQR2) NSYMF = ISYMI
       DO ISYMF = 1, NSYMF

        DO IEXI = 1, NCCEXCI(ISYMI,1)
         NEXF = NCCEXCI(ISYMF,1)
         IF ((.NOT.SELQR2).AND.(ISYMI.EQ.ISYMF)) NEXF = IEXI - 1
         DO IEXF = 1, NEXF

          IF (LOCDBG) THEN
            WRITE(LUPRI,*) 'check for:',ISYMI,IEXI,ISYMF,IEXF
          END IF

          IF (SELQR2) THEN
           ! check, if state pair has been selected
           TAKE_STATE_PAIR = .FALSE.
           DO I = 1,NSEQR2
            IF (ISYMI.EQ.ISEQR2SYM(I,1).AND.IEXI.EQ.ISEQR2STA(I,1).AND.
     *          ISYMF.EQ.ISEQR2SYM(I,2).AND.IEXF.EQ.ISEQR2STA(I,2)    )
     *      THEN
             TAKE_STATE_PAIR = .TRUE.
            ENDIF 
           END DO
          ELSE
           TAKE_STATE_PAIR = .TRUE.
          END IF

          IF (LOCDBG) WRITE(LUPRI,*) 'TAKE_STATE_PAIR:',TAKE_STATE_PAIR

          IF (TAKE_STATE_PAIR) THEN
            NXQR2ST = NXQR2ST + 1
            IF (NXQR2ST.GT.MXQR2ST)  THEN
              WRITE(LUPRI,*) 'NXQR2ST =',NXQR2ST,'MXQR2ST= ',MXQR2ST
              CALL QUIT(' Too many states in CC_XOPAIND')
            ENDIF
            ISTATEI = ISYOFE(ISYMI)+IEXI
            ISTATEF = ISYOFE(ISYMF)+IEXF  
            IQR2ST(NXQR2ST,1) = ISTATEI
            IQR2ST(NXQR2ST,2) = ISTATEF

            IF (.NOT.CIS) THEN
              ISYMFI = MULD2H(ISYMI,ISYMF)
              EIGVI = EIGVAL(ISTATEI)
              EIGVF = EIGVAL(ISTATEF)

              DO IOPER = 1, NQR2OP
               LABEL = LBLOPR(IQR2OP(IOPER))
               ISYMO = ISYOPR(IQR2OP(IOPER))
               IF (ISYMO.EQ.ISYMFI) THEN
                IF (.NOT.QR22N1) THEN
                 INUM=IR1TAMP(LABEL,.FALSE.,EIGVI-EIGVF,ISYMO)
                 INUM=IR1TAMP(LABEL,.FALSE.,EIGVF-EIGVI,ISYMO)
                ELSE
                 INUM=IN2AMP(ISTATEI,-EIGVI,ISYMI,ISTATEF,+EIGVF,ISYMF)
                 INUM=IN2AMP(ISTATEF,-EIGVF,ISYMF,ISTATEI,+EIGVI,ISYMI)
                END IF
               END IF
              END DO
            END IF 
          END IF 

         END DO
        END DO
       END DO
      END DO
C
      IF (LOCDBG) THEN
        WRITE(LUPRI,'(a,i3,a)') 
     &    'Transition strengths will be computed for',NXQR2ST,
     &    'transitions'
      END IF
C
      RETURN
      END
*---------------------------------------------------------------------*
c /* deck cc_lrsind */
*=====================================================================*
       SUBROUTINE CC_LRSIND
*---------------------------------------------------------------------*
*
*    Purpose: Control input and equations for calculation of 
*             residues of the linear response function.
*
*    OC 8-11-1996/Modified April 1997
*
*=====================================================================*
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <cclrinf.h>
#include <ccrspprp.h>
#include <ccexci.h>  
#include <cclres.h>  
#include <ccsdinp.h>  
#include <ccsdsym.h>  
#include <cclr.h>  
#include <ccroper.h>  
      CHARACTER*8 LABEL,LABELA,LABELB
      INTEGER ISYMB, ISYMA, IFRB, IFRA, INUM, IOP, ISYME

      REAL*8  EIGV
* functions:
      INTEGER IR1TAMP,ILRMAMP
C
      LOGICAL FIRSTCALL
      SAVE    FIRSTCALL
      DATA    FIRSTCALL /.TRUE./
C
C------------------------------------------------------------------
C     Add residue response equations to list to be solved for CCLR. 
C------------------------------------------------------------------
C
*------------------------------------------------------------------
* test if operators are available and translate IALROP, IBLROP
* arrays from the PRPLBL list to the new list maintained by IROPER.
*------------------------------------------------------------------
      IF (FIRSTCALL) THEN
       IOPER = 1
       DO WHILE (IOPER .LE. NLRSOP)
        LABELA = PRPLBL(IALRSOP(IOPER))
        LABELB = PRPLBL(IBLRSOP(IOPER))
        IF (DEBUG) THEN
          WRITE(LUPRI,'(/2X,A,3(1X,A),A)')
     &     'CHECK LRSD DOUBLE:',LABELA, LABELB
        ENDIF
        IF (      (IROPER(LABELA,ISYMA) .LT. 0)
     &       .OR. (IROPER(LABELB,ISYMB) .LT. 0) ) THEN
          WRITE(LUPRI,'(/2X,5A,/2X,2A)')
     &     ' WARNING: ONE OF THE OPERATORS WITH THE LABELS "',
     &     LABELA,'", "', LABELB,'" IS NOT AVAILABLE.',
     &     ' LINEAR RESPONSE RESIDUE CALCULATION IS CANCELED FOR THIS',
     &     ' OPERATOR DOUBLE.'
          DO IDX = IOPER, NLRSOP-1
            IALRSOP(IDX) = IALRSOP(IDX+1)
            IBLRSOP(IDX) = IBLRSOP(IDX+1)
          END DO
          NLRSOP = NLRSOP - 1
        ELSE
          IF (DEBUG) THEN
             WRITE(LUPRI,'(/2X,A,3(1X,A),A)')
     &        'PUT DOUBLE:',LABELA, LABELB,' ON THE LIST.'
          ENDIF
          IALRSOP(IOPER) = IROPER(LABELA,ISYMA)
          IBLRSOP(IOPER) = IROPER(LABELB,ISYMB)
          IOPER = IOPER + 1
        END IF
       END DO
       FIRSTCALL = .FALSE.
      END IF ! (FIRSTCALL)
C
      IF (DEBUG) THEN
         WRITE(LUPRI,'(/,A)') ' Updated list'
         DO IOPER = 1, NLRSOP
            WRITE(LUPRI,*) IOPER,IALRSOP(IOPER),' ',
     *              LBLOPR(IALRSOP(IOPER)),
     *              ISYOPR(IALRSOP(IOPER)),IBLRSOP(IOPER),
     *              ' ',LBLOPR(IBLRSOP(IOPER)),ISYOPR(IBLRSOP(IOPER))
         ENDDO
      ENDIF
C
C------------------------------------
C     Fill in equations to be solved.
C------------------------------------
C
      NXLRSST = 0
C
      DO 100 ISYME = 1, NSYM
       DO 200 IEX = 1, NCCEXCI(ISYME,1)
C
        IF (SELLRS) THEN
C
C        Check state has been calculated. 
C
           DO 300 I = 1,NSELRS
             IF ((ISYME.EQ.ISELRS(I,1))
     *           .AND.(IEX.EQ.ISELRS(I,2))) THEN
                NXLRSST = NXLRSST + 1
                ILRSST(NXLRSST) = ISYOFE(ISYME)+IEX 
                GO TO 350 
             ENDIF 
  300     CONTINUE
C 
C------------------------------------------------------------------------
C         This state is not requested in oscillator strength calculation.
C         GOTO end of loop.
C------------------------------------------------------------------------
C
          GO TO 200
        ELSE
          NXLRSST = NXLRSST + 1
          IF (NXLRSST.GT.MXLRSST)  THEN
            WRITE(LUPRI,*) 'NXLRSST =',NXLRSST,'MXLRSST= ',MXLRSST
            CALL QUIT(' Too many states for residue calculation')
          ENDIF
          ILRSST(NXLRSST) = ISYOFE(ISYME)+IEX
        ENDIF

  350   CONTINUE
          
        IF (.NOT.CIS) THEN
          DO 400 IOPER = 1, NLRSOP

            ISYMA  = ISYOPR(IALRSOP(IOPER))
            ISYMB  = ISYOPR(IBLRSOP(IOPER))

            IF ((ISYME.EQ.ISYMA).AND.(ISYME.EQ.ISYMB)) THEN

              LABELA = LBLOPR(IALRSOP(IOPER))
              LABELB = LBLOPR(IBLRSOP(IOPER))
              IF (.NOT.LRS2N1) THEN
                 EIGV  = -EIGVAL(ILRSST(NXLRSST)) 
                 INUM  = IR1TAMP(LABELB,.FALSE.,EIGV,ISYMB)
              ELSE
                 EIGV  = EIGVAL(ILRSST(NXLRSST))   
                 INUM  = ILRMAMP(ILRSST(NXLRSST),EIGV,ISYMB)
              ENDIF

            ENDIF
  400     CONTINUE
        ENDIF 
  200  CONTINUE
  100 CONTINUE
C
      END
*---------------------------------------------------------------------*
c /* deck cc_qr2ind */
*=====================================================================*
       SUBROUTINE CC_QR2IND
*---------------------------------------------------------------------*
*
*    Purpose: Control input and equations for calculation of 
*             second residues of the quadratic response function.
*
*    Ove Christiansen April 1997
*
*=====================================================================*
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <cclrinf.h>
#include <ccrspprp.h>
#include <ccexci.h>  
#include <cclres.h>  
#include <ccqr2r.h>  
#include <ccn2rsp.h>
#include <ccsdinp.h>  
#include <ccsdsym.h>  
#include <cclr.h>  
#include <ccroper.h>  
      CHARACTER*8 LABEL,LABELA,LABELB
      INTEGER ISYMB, ISYMA, IFRB, IFRA, INUM, IOP, ISYMI, ISYMF

      REAL*8  EIGV

* functions:
      INTEGER IR1TAMP,IN2AMP
C
      LOGICAL FIRSTCALL
      SAVE    FIRSTCALL
      DATA    FIRSTCALL /.TRUE./
C
C------------------------------------------------------------------
C     Add residue response equations to list to be solved for CCQR2 
C------------------------------------------------------------------
C
*------------------------------------------------------------------
* test if operators are available and translate IALROP, IBLROP
* arrays from the PRPLBL list to the new list maintained by IROPER.
*------------------------------------------------------------------
      IF (FIRSTCALL) THEN
       IOPER = 1
       DO WHILE (IOPER .LE. NQR2OP)
        LABELA = PRPLBL(IAQR2OP(IOPER))
        LABELB = PRPLBL(IBQR2OP(IOPER))
        IF (DEBUG) THEN
          WRITE(LUPRI,'(/2X,A,3(1X,A),A)')
     &     'CHECK QR2R DOUBLE:',LABELA, LABELB
        ENDIF
        IF (      (IROPER(LABELA,ISYMA) .LT. 0)
     &       .OR. (IROPER(LABELB,ISYMB) .LT. 0) ) THEN
          WRITE(LUPRI,'(/2X,5A,/2X,2A)')
     &     ' WARNING: ONE OF THE OPERATORS WITH THE LABELS "',
     &     LABELA,'", "', LABELB,'" IS NOT AVAILABLE.',
     &     ' QUADRATIC RESPONSE 2. RESIDUE CALCULATION IS CANCELED ',
     &     'FOR THIS OPERATOR DOUBLE.'
          DO IDX = IOPER, NQR2OP-1
            IAQR2OP(IDX) = IAQR2OP(IDX+1)
            IBQR2OP(IDX) = IBQR2OP(IDX+1)
          END DO
          NQR2OP = NQR2OP - 1
        ELSE
          IF (DEBUG) THEN
             WRITE(LUPRI,'(/2X,A,3(1X,A),A)')
     &        'PUT DOUBLE:',LABELA, LABELB,' ON THE LIST.'
          ENDIF
          IAQR2OP(IOPER) = IROPER(LABELA,ISYMA)
          IBQR2OP(IOPER) = IROPER(LABELB,ISYMB)
          IOPER = IOPER + 1
        END IF
       END DO
       FIRSTCALL = .FALSE.
      END IF ! (FIRSTCALL)
C
      IF (DEBUG) THEN
         WRITE(LUPRI,'(/,A)') ' Updated list'
         DO IOPER = 1, NQR2OP
            WRITE(LUPRI,*) IOPER,IAQR2OP(IOPER),' ',
     *                     LBLOPR(IAQR2OP(IOPER)),
     *                     ISYOPR(IAQR2OP(IOPER)),
     *                     IBQR2OP(IOPER),' ',
     *                     LBLOPR(IBQR2OP(IOPER)),
     *                     ISYOPR(IBQR2OP(IOPER))
         ENDDO
      ENDIF
C
C------------------------------------
C     Fill in equations to be solved.
C------------------------------------
C
      NXQR2ST = 0
C
      DO 100 ISYMFI = 1, NSYM
       DO 200 ISYMF = 1, NSYM
        ISYMI = MULD2H(ISYMF,ISYMFI)
        IF ((.NOT.SELQR2).AND.(ISYMI.GT.ISYMF)) GOTO 200 
        DO 300 IEXF = 1, NCCEXCI(ISYMF,1)
         NEXI = NCCEXCI(ISYMI,1)
         IF ((.NOT.SELQR2).AND.(ISYMI.EQ.ISYMF)) NEXI = IEXF - 1
         DO 400 IEXI = 1, NEXI
C
          IF (SELQR2) THEN
C
C        Check state set has been selected and calculated.
C
            DO 500 I = 1,NSEQR2
             IF ((ISYMI.EQ.ISEQR2(I,1))
     *           .AND.(IEXI.EQ.ISEQR2(I,2))
     *           .AND.(ISYMF.EQ.ISEQR2(I,3))
     *           .AND.(IEXF.EQ.ISEQR2(I,4))) THEN
                NXQR2ST = NXQR2ST + 1
                IQR2STI(NXQR2ST) = ISYOFE(ISYMI)+IEXI
                IQR2STF(NXQR2ST) = ISYOFE(ISYMF)+IEXF
                GO TO 550 
             ENDIF 
  500       CONTINUE
C 
C--------------------------------------------------------------------------
C           This state is not requested in oscillator strength calculation.
C           GOTO end of loop.
C--------------------------------------------------------------------------
C
            GO TO 400

          ELSE
            NXQR2ST = NXQR2ST + 1
            IF (NXQR2ST.GT.MXQR2ST)  THEN
              WRITE(LUPRI,*) 'NXQR2ST =',NXQR2ST,'MXQR2ST= ',MXQR2ST
              CALL QUIT(' Too many states for residue calculation')
            ENDIF
            IQR2STI(NXQR2ST) = ISYOFE(ISYMI)+IEXI
            IQR2STF(NXQR2ST) = ISYOFE(ISYMF)+IEXF
          ENDIF

  550     CONTINUE
          
          IF (.NOT.CIS) THEN
            DO 600 IOPER = 1, NQR2OP
    
              ISYMA  = ISYOPR(IAQR2OP(IOPER))
              ISYMB  = ISYOPR(IBQR2OP(IOPER))
              ISYMAI = MULD2H(ISYMA,ISYMI)
              ISYMBF = MULD2H(ISYMB,ISYMF)
              IF ((ISYMAI.EQ.ISYMF).AND.(ISYMBF.EQ.ISYMI)) THEN
                LABELA = LBLOPR(IAQR2OP(IOPER))
                LABELB = LBLOPR(IBQR2OP(IOPER))
                IF (.NOT.QR22N1) THEN
                   EIGVI = EIGVAL(IQR2STI(NXQR2ST)) 
                   EIGVF = EIGVAL(IQR2STF(NXQR2ST)) 
                   EIGV  = EIGVI - EIGVF
                   INUM  = IR1TAMP(LABELB,.FALSE.,EIGV,ISYMB)
                   EIGV  = - EIGVI + EIGVF
                   INUM  = IR1TAMP(LABELA,.FALSE.,EIGV,ISYMA)
                ELSE
                   EIGVI = EIGVAL(IQR2STI(NXQR2ST)) 
                   EIGVF = EIGVAL(IQR2STF(NXQR2ST)) 
                   INUM  = IN2AMP(IQR2STI(NXQR2ST),-EIGVI,ISYMI,
     *                            IQR2STF(NXQR2ST),EIGVF,ISYMF)
                   INUM  = IN2AMP(IQR2STF(NXQR2ST),-EIGVF,ISYMF,
     *                            IQR2STI(NXQR2ST),EIGVI,ISYMI)
                ENDIF
  
              ENDIF
  600       CONTINUE
          ENDIF 
  400    CONTINUE
  300   CONTINUE
  200  CONTINUE
  100 CONTINUE
C
      END
*---------------------------------------------------------------------*
c /* deck cc_lrind */
*=====================================================================*
       SUBROUTINE CC_LRIND(WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: determine which the response t amplitudes and zeta
*             multipliers required for the dynamic polarizabilities
*
*    Written by Christof Haettig, October 1996.
*     
*    OC 32-10-1996: ASYMSD option.
*    OC  dec. 1996: Cauchy moment section.
*    CH  oct. 1997: ASYMSD option for Cauchy moment section.
*    CH  nov. 1998: relaxed response.
*    CH  feb. 1999: missing (anti-)symmetrization in +/- w introduced.
*    CH  may. 1999: changes for first-order property gradients
*    CH  apr. 2002: changes for CC3 freq.-dep. polarizabilities
*
*=====================================================================*
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <cclrinf.h>
#include <ccrspprp.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <ccroper.h>
#include <mxcent.h>
#include <cclr.h>

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

* variables:
      CHARACTER*8 LABEL, LABELA, LABELB, LABSOP
      CHARACTER*3 LSTRLX
      LOGICAL SKIP_IT, LRLXA, LRLXB, LPDBSA, LPDBSB
      LOGICAL DIFDIP, SYM1ONLY
      INTEGER ISYMB,ISYMA,IFRB,IFRA,INUM,IOP,ICAUCH,ISYH0,IR1A,IR1B
      INTEGER ISYM0, ISYSOP, ISGNSOP
      
      REAL*8  WORK(LWORK), FREQ

* functions:
      INTEGER IR1TAMP
      INTEGER IL1ZETA
      INTEGER ILRCAMP
      INTEGER ILC1AMP
      INTEGER IROPER
      INTEGER IETA1

      LOGICAL FIRSTCALL
      SAVE    FIRSTCALL
      DATA    FIRSTCALL /.TRUE./

*------------------------------------------------------------------
* test if operators are available and translate IALROP, IBLROP
* arrays from the PRPLBL list to the new list maintained by IROPER.
*------------------------------------------------------------------
      IF (FIRSTCALL) THEN

        IOPER = 1
        DO WHILE (IOPER .LE. NLROP)
          SKIP_IT = .FALSE.
          LABELA  = PRPLBL(IALROP(IOPER))
          LABELB  = PRPLBL(IBLROP(IOPER))
          IOPA    = IROPER(LABELA,ISYMA)
          IOPB    = IROPER(LABELB,ISYMB)

          IF ( IOPA.LT.0 .OR. IOPB.LT.0 ) THEN
            WRITE(LUPRI,'(/2X,5A,/2X,2A)')
     &          ' WARNING: ONE OF THE OPERATORS WITH THE LABELS "',
     &          LABELA,'", "', LABELB,'" IS NOT AVAILABLE.',
     &          ' LINEAR RESPONSE CALCULATION IS CANCELED FOR THIS',
     &          ' OPERATOR PAIR.'
            SKIP_IT = .TRUE.
          END IF

          IF (.NOT.SKIP_IT) THEN
             ! if we have field-dependent basis sets, we need also 
             ! to check, if the second-derivative integrals for this
             ! perturbation pair are available
             IF (LPDBSOP(IOPA) .OR. LPDBSOP(IOPB)) THEN
                CALL CC_FIND_SO_OP(LABELA,LABELB,LABSOP,ISYSOP,
     &                             ISGNSOP,INUM,WORK,LWORK)
                IF (INUM.LT.0) SKIP_IT = .TRUE.
             END IF
          END IF

          IF (SKIP_IT) THEN
            DO IDX = IOPER, NLROP-1
              IALROP(IDX) = IALROP(IDX+1)
              IBLROP(IDX) = IBLROP(IDX+1)
              LALORX(IDX) = LALORX(IDX+1)
              LBLORX(IDX) = LBLORX(IDX+1)
            END DO
            NLROP = NLROP - 1
          ELSE
            IALROP(IOPER) = IROPER(LABELA,ISYMA)
            IBLROP(IOPER) = IROPER(LABELB,ISYMB)
            IOPER = IOPER + 1
          END IF

        END DO
        FIRSTCALL = .FALSE.
      END IF ! (FIRSTCALL)


*---------------------------------------------------------------------*
* set: a) linear response equations to be solved 
*      b) effective Fock operators to be calculated
*      c) nuclear contributions to be calculated
*---------------------------------------------------------------------*
      DIFDIP = .FALSE.

      DO IOPER = 1, NLROP
        LABELA = LBLOPR(IALROP(IOPER))
        LABELB = LBLOPR(IBLROP(IOPER))

        ISYMA  = ISYOPR(IALROP(IOPER))
        ISYMB  = ISYOPR(IBLROP(IOPER))

        LRLXA  = LALORX(IOPER)
        LRLXB  = LBLORX(IOPER)

        LPDBSA = LPDBSOP(IALROP(IOPER))
        LPDBSB = LPDBSOP(IBLROP(IOPER))

        IF ((LRLXA.OR.LRLXB.OR.LPDBSA.OR.LPDBSB) .AND. CC3)
     &    CALL QUIT('Relaxed CC3 linear response no implemented.')

        IF (ISYMA.EQ.ISYMB) THEN

          IF (DEBUG) THEN
            WRITE(LUPRI,'(/2X,A,2(1X,A,2L1))')
     &       'require linear responses for double:',
     &        LABELA, LRLXA, LPDBSA, LABELB, LRLXB, LPDBSB
          ENDIF
          DO IFREQ = 1, NBLRFR 
          DO ISIGN = +1, -1, -2
 
            SIGN   = DBLE(ISIGN)
            FREQA  = SIGN * ALRFR(IFREQ)
            FREQB  = SIGN * BLRFR(IFREQ)

            INUM  = IR1TAMP(LABELB,LRLXB,FREQB,ISYMB)
            IF (LRLXB.OR.LPDBSB) INUM = IETA1(LABELB,LRLXB,FREQB,ISYMB)
            IF (CCSLV)        INUM  = IL1ZETA(LABELB,LRLXB,FREQB,ISYMB)
            IF (.NOT. ASYMSD) THEN
              INUM = IR1TAMP(LABELA,LRLXA,FREQA,ISYMA)
              IF (LRLXA.OR.LPDBSA) INUM=IETA1(LABELA,LRLXA,FREQA,ISYMA)
              IF (CCSLV)         INUM=IL1ZETA(LABELA,LRLXA,FREQA,ISYMA)
            ELSE
               INUM  = IL1ZETA(LABELB,LRLXB,FREQB,ISYMB)
            ENDIF


            IF (LRLXB .OR. LPDBSB) THEN
               INUM = IEFFFOCK(LABELA,ISYMA,1)
               INUM =  IEXPECT(LABELA,ISYMA,1)
            END IF
            IF (LRLXA .OR. LPDBSA) THEN
               INUM = IEFFFOCK(LABELB,ISYMB,1)
               INUM =  IEXPECT(LABELB,ISYMB,1)
            END IF
            IF ((LRLXB.OR.LPDBSB) .AND. (LRLXA.OR.LPDBSA)) THEN
               INUM  = I1DXFCK('HAM0    ','R1 ',LABELA,FREQA,ISYMA)
               INUM  = I1DXFCK('HAM0    ','R1 ',LABELB,FREQB,ISYMB)
               INUM  = IEFFFOCK('HAM0    ',1,1)
               INUM  =  IEXPECT('HAM0    ',1,1)
            END IF

            IF (LPDBSA .OR. LPDBSB) THEN
              CALL CC_FIND_SO_OP(LABELA,LABELB,LABSOP,ISYSOP,
     &                           ISGNSOP,INUM,WORK,LWORK)
              IF (INUM.LT.0) CALL QUIT('Operator error in CC_LRIND.')
              INUM = IEFFFOCK(LABSOP,ISYSOP,2)
              INUM = IEXPECT(LABSOP,ISYSOP,2)

              IF (LABSOP(4:6).EQ.'DPG') THEN
                DIFDIP = .TRUE.
              ELSE IF (LABSOP(3:5).EQ.'QDG') THEN
                CONTINUE
              ELSE IF (LABSOP(3:5).EQ.'OCG') THEN
                CONTINUE
              ELSE IF (LABSOP(2:6).EQ.'-CM1 ') THEN
                CONTINUE
              ELSE IF (LABSOP(4:7).EQ.' NST') THEN
                CONTINUE
              ELSE
                WRITE (LUPRI,*) 
     &                'Illegal or unknown label in CC_LRIND:',LABSOP
                CALL QUIT('Illegal or unknown label in CC_LRIND.')
              END IF
            END IF

          END DO
          END DO
        END IF

      END DO
C
C     Note: this is required to get CAUCHY vectors in correct order.
C
      IF (CAUCHY) THEN

        ! switch off a special treatment of cauchy vectors in the
        ! solver which cannot be used with CC3
        IF (CC3) NEWCAU = .FALSE.

        DO ISYM = 1, NSYM
          DO ICAUCH  = 1, NLRDISP
            DO IOPER = 1, NLROP
              LABELA = LBLOPR(IALROP(IOPER))
              LABELB = LBLOPR(IBLROP(IOPER))
              ISYMA  = ISYOPR(IALROP(IOPER))
              ISYMB  = ISYOPR(IBLROP(IOPER))
              LRLXA  = LALORX(IOPER)
              LRLXB  = LBLORX(IOPER)
              LPDBSA = LPDBSOP(IALROP(IOPER))
              LPDBSB = LPDBSOP(IBLROP(IOPER))

              IF (LRLXA .OR. LRLXB) THEN
                 WRITE (LUPRI,*)
     &                'Warning: orbital relaxation is ignored ',
     &                    'in the calculation of Cauchy moments.'
              END IF

              IF (LPDBSA .OR. LPDBSB) THEN
                 WRITE (LUPRI,*)
     &                'Error: Cauchy moments not implemented',
     &                    'for field-dependent basis sets.'
                 CALL QUIT('No Cauchy moments for '//
     &                'field-dep. basis sets.')
              END IF

              IF ((ISYMA.EQ.ISYMB).AND.(ISYM.EQ.ISYMA)) THEN
                INUM  = ILRCAMP(LABELB,ICAUCH,ISYMB)
                IF (ASYMSD) THEN
                  INUM  = ILC1AMP(LABELB,ICAUCH,ISYMB)
                ELSE
                  INUM  = ILRCAMP(LABELA,ICAUCH,ISYMA)
                END IF
              END IF
            END DO
          END DO
        END DO
      END IF
C
C     let abacus precalculate nuclear contributions:
C
      IF (DIFDIP) THEN
         KCSTRA = 1
         KSCTRA = KCSTRA + MXCOOR*MXCOOR
         KEND   = KSCTRA + MXCOOR*MXCOOR
         LEND   = LWORK  - KEND
 
         IF (LEND.LT.0) THEN
            CALL QUIT('Insufficient memory in CC_LRIND.')
         END IF

         SYM1ONLY = .FALSE.
         CALL CC_SETDORPS('1DHAM   ',SYM1ONLY,0)
         CALL DIPNUC(WORK(KCSTRA),WORK(KSCTRA),IPRINT,DIFDIP)

      END IF

      RETURN
      END
*---------------------------------------------------------------------*
c /* deck cc_qrind */
*=====================================================================*
       SUBROUTINE CC_QRIND(WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: Determine which response t amplitudes and zeta
*             multipliers required for the first hyperpolarizabilities
*             and their dispersion coefficients
*
*    Written by Christof Haettig, October 1996.
*    Dispersion coefficients, October 1997 (Christof Haettig)
*    Relaxed response for one of the operators, June 1999 (Ch. Haettig)
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccqrinf.h>
#include <ccrspprp.h>
#include <ccroper.h>
#include <ccsdinp.h>

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

* variables:
      CHARACTER*8 LABELA, LABELB, LABELC, LABSOP
      LOGICAL LORXA,LORXB,LORXC, LPDBSA,LPDBSB,LPDBSC, SKIP_IT
      LOGICAL LRELAX
      INTEGER ISYMB, ISYMC, ISYMA, IFREQ, IDISP, INUM, IOPER, IDX
      INTEGER ICA,ICB,ICC,ICTOT,ISACAU,ISAMA,ISAMB,ISAMC,ISAPROP
      INTEGER IOPA, IOPB, IOPC, NLORX, ISYSOP, LWORK, ISGNSOP
      
      REAL*8  WORK
      REAL*8  FREQA, FREQB, FREQC

* external functions:
      INTEGER IR1TAMP
      INTEGER IL1ZETA
      INTEGER IROPER
      INTEGER ILRCAMP
      INTEGER ILC1AMP
      INTEGER ICR2AMP
      INTEGER IR2TAMP

* data:
      LOGICAL FIRSTCALL
      SAVE    FIRSTCALL
      DATA    FIRSTCALL /.TRUE./
      CHARACTER*7 CISA(-1:1)
      DATA    CISA /'odd    ','unknown','even   '/


*---------------------------------------------------------------------*
* test if operators are available and translate IAQROP, IBQROP, ICQROP
* arrays from the PRPLBL list to the new list maintained by IROPER.
*---------------------------------------------------------------------*
      IF (FIRSTCALL) THEN

       IOPER = 1
       DO WHILE (IOPER .LE. NQROPER)

        SKIP_IT = .FALSE.
        LABELA  = PRPLBL(IAQROP(IOPER))
        LABELB  = PRPLBL(IBQROP(IOPER))
        LABELC  = PRPLBL(ICQROP(IOPER))
        IOPA    = IROPER(LABELA,ISYMA)
        IOPB    = IROPER(LABELB,ISYMB)
        IOPC    = IROPER(LABELC,ISYMC)
        LORXA   = LAQLRX(IOPER)
        LORXB   = LBQLRX(IOPER)
        LORXC   = LCQLRX(IOPER)

        IF (LOCDBG) THEN
          WRITE(LUPRI,'(/2X,A,3(1X,A),A)')
     &     'CHECK TRIPLE:',LABELA, LABELB, LABELC
        END IF

        IF ( IOPA.LT.0 .OR. IOPB.LT.0 .OR. IOPC.LT.0 ) THEN
          WRITE(LUPRI,'(/2X,7A,/2X,2A)')
     &     ' WARNING: ONE OF THE OPERATORS WITH THE LABELS "',
     &     LABELA,'", "', LABELB,'", "', LABELC,'" IS NOT AVAILABLE.',
     &     ' HYPERPOLARIZABILITY CALCULATION IS CANCELED FOR THIS',
     &     ' OPERATOR TRIPLE.'
           SKIP_IT = .TRUE.
        END IF

        NLORX = 0
        IF (LORXA .OR. LPDBSOP(IOPA)) NLORX = NLORX + 1
        IF (LORXB .OR. LPDBSOP(IOPB)) NLORX = NLORX + 1
        IF (LORXC .OR. LPDBSOP(IOPC)) NLORX = NLORX + 1
        
        IF (NLORX.GT.1) THEN
          WRITE(LUPRI,'(/2X,8A,/2X,A,/2X,A)') 
     &     ' WARNING: OPERATOR TRIPLETT "',
     &     LABELA,'", "', LABELB,'", "', LABELC,'"',
     &     ' WITH MORE THAN ONE FIELD WHICH',
     &     ' INVOKES ORBITAL RELAXATION OR A PERTUR.-DEP. BASIS SET.',
     &     ' CALCULATION IS CANCELED FOR THIS OPERATOR TRIPLE.'
        END IF

        IF (USE_R2 .AND. NLORX.GT.0) THEN
           WRITE (LUPRI,*) 'Second-order response vectors not yet',
     &                ' implemented for fields which invoke'
           WRITE(LUPRI,*)
     &          'orbital relaxation or perturb.-dep. basis sets.'
           WRITE(LUPRI,*) 'USE_R2 option turned off.'
           USE_R2 = .FALSE.
        END IF
              
        IF (.NOT. SKIP_IT) THEN
             ! if we have field-dependent basis sets, we need also
             ! to check, if the second-derivative integrals for this
             ! perturbation pair are available
             IF (LPDBSOP(IOPA) .OR. LPDBSOP(IOPB)) THEN
                CALL CC_FIND_SO_OP(LABELA,LABELB,LABSOP,ISYSOP,
     &                             ISGNSOP,INUM,WORK,LWORK)
                IF (INUM.LT.0) SKIP_IT = .TRUE.
             END IF                                     
             IF (LPDBSOP(IOPA) .OR. LPDBSOP(IOPC)) THEN
                CALL CC_FIND_SO_OP(LABELA,LABELC,LABSOP,ISYSOP,
     &                             ISGNSOP,INUM,WORK,LWORK)
                IF (INUM.LT.0) SKIP_IT = .TRUE.
             END IF                                     
             IF (LPDBSOP(IOPB) .OR. LPDBSOP(IOPC)) THEN
                CALL CC_FIND_SO_OP(LABELB,LABELC,LABSOP,ISYSOP,
     &                             ISGNSOP,INUM,WORK,LWORK)
                IF (INUM.LT.0) SKIP_IT = .TRUE.
             END IF                                   
             iF (SKIP_IT) THEN
               WRITE(LUPRI,'(/2X,7A,/2X,A,/2X,A)') 
     &          ' WARNING: FOR THE OPERATOR TRIPLETT "',
     &          LABELA,'", "', LABELB,'", "', LABELC,'"',
     &         ' A SEC. ORD. OPERATOR IS MISSING.',
     &         ' CALCULATION IS IGNORED.'
             END IF
        END IF

        IF (SKIP_IT) THEN
          DO IDX = IOPER, NQROPER-1
            IAQROP(IDX) = IAQROP(IDX+1)
            IBQROP(IDX) = IBQROP(IDX+1)
            ICQROP(IDX) = ICQROP(IDX+1)
            LAQLRX(IDX) = LAQLRX(IDX+1)
            LBQLRX(IDX) = LBQLRX(IDX+1)
            LCQLRX(IDX) = LCQLRX(IDX+1)
          END DO
          NQROPER = NQROPER - 1
        ELSE 
          WRITE(LUPRI,'(/2X,A,3(1X,A),A)')
     &     'PUT TRIPLE:',LABELA, LABELB, LABELC,' ON THE LIST.'
          IAQROP(IOPER) = IROPER(LABELA,ISYMA) 
          IBQROP(IOPER) = IROPER(LABELB,ISYMB)
          ICQROP(IOPER) = IROPER(LABELC,ISYMC)
          IOPER = IOPER + 1
        END IF

       END DO

       FIRSTCALL = .FALSE.

      END IF ! (FIRSTCALL)

*---------------------------------------------------------------------*
* set list entries for the required response vectors:
*---------------------------------------------------------------------*
      IF (CC3) THEN
        WRITE(LUPRI,'(/5x,A/)')
     &   'Prepare CC3 quadratic response calculation.'
        IF (NQRDISP.GT.0) THEN
          NQRDISP = 0
          WRITE(LUPRI,'(/5x,A//)') 
     &    'Dispersion coefficients (.DISPCF) are switched off for CC3.'
        END IF
        IF (USE_R2) THEN
          WRITE(LUPRI,'(2(/5x,A),/)') 
     &    'Note: .USE R2 option will for CC3 call noddy code routines,',
     &    '      which keep triples amplitudes in memory!!!'
        END IF
      END IF

      DO IOPER = 1, NQROPER
        LABELA = LBLOPR(IAQROP(IOPER))
        LABELB = LBLOPR(IBQROP(IOPER))
        LABELC = LBLOPR(ICQROP(IOPER))

        LPDBSA = LPDBSOP(IAQROP(IOPER))
        LPDBSB = LPDBSOP(IBQROP(IOPER))
        LPDBSC = LPDBSOP(ICQROP(IOPER))

        ISYMA  = ISYOPR(IAQROP(IOPER))
        ISYMB  = ISYOPR(IBQROP(IOPER))
        ISYMC  = ISYOPR(ICQROP(IOPER))

        ISAMA  = ISYMAT(IAQROP(IOPER))
        ISAMB  = ISYMAT(IBQROP(IOPER))
        ISAMC  = ISYMAT(ICQROP(IOPER))
  
        IOPA   = IROPER(LABELA,ISYMA)
        IOPB   = IROPER(LABELB,ISYMB)
        IOPC   = IROPER(LABELC,ISYMC)

        ISAPROP = ISAMA * ISAMB * ISAMC

        LORXA  = LAQLRX(IOPER)
        LORXB  = LBQLRX(IOPER)
        LORXC  = LCQLRX(IOPER)

        LRELAX = LORXA.OR.LORXB.OR.LORXC.OR.LPDBSA.OR.LPDBSB.OR.LPDBSC
        IF (LRELAX.AND.CC3) CALL QUIT('No relaxed CC3 quadratic resp.')

c         WRITE(LUPRI,'(/2X,A,3(1X,A),A)')
c    &     'require responses for triple:',LABELA, LABELB, LABELC
c         WRITE(LUPRI,'(/2X,A,A)')
c    &     'symmetry in the sign of the frequency is ',CISA(ISAPROP)

        
        IF (MULD2H(ISYMA,ISYMB).EQ.ISYMC) THEN
      
*          if we have field-dependent basis sets:
*          --------------------------------------
*          we need to check, if the second-derivative integrals 
*          for the perturbation pair are available
           IF (LPDBSOP(IOPA) .OR. LPDBSOP(IOPB)) THEN
              CALL CC_FIND_SO_OP(LABELA,LABELB,LABSOP,ISYSOP,
     &                           ISGNSOP,INUM,WORK,LWORK)
           END IF                                     
           IF (LPDBSOP(IOPA) .OR. LPDBSOP(IOPC)) THEN
              CALL CC_FIND_SO_OP(LABELA,LABELC,LABSOP,ISYSOP,
     &                           ISGNSOP,INUM,WORK,LWORK)
           END IF                                     
           IF (LPDBSOP(IOPB) .OR. LPDBSOP(IOPC)) THEN
              CALL CC_FIND_SO_OP(LABELB,LABELC,LABSOP,ISYSOP,
     &                           ISGNSOP,INUM,WORK,LWORK)
           END IF                                   

*         for frequency-dependent hyperpolarizabilities:
*         ----------------------------------------------
          DO IFREQ = 1, NQRFREQ
            FREQA  = AQRFR(IFREQ)
            FREQB  = BQRFR(IFREQ)
            FREQC  = CQRFR(IFREQ)


*           request (unrelaxed) first-order t response vectors:

            INUM = IR1TAMP(LABELA,LORXA,+FREQA,ISYMA)
            INUM = IR1TAMP(LABELB,LORXB,+FREQB,ISYMB)
            INUM = IR1TAMP(LABELC,LORXC,+FREQC,ISYMC)
            INUM = IR1TAMP(LABELA,LORXA,-FREQA,ISYMA)
            INUM = IR1TAMP(LABELB,LORXB,-FREQB,ISYMB)
            INUM = IR1TAMP(LABELC,LORXC,-FREQC,ISYMC)


*           request first-order zeta response vectors:

            INUM = IL1ZETA(LABELA,LORXA,+FREQA,ISYMA)
            INUM = IL1ZETA(LABELB,LORXB,+FREQB,ISYMB)
            INUM = IL1ZETA(LABELC,LORXC,+FREQC,ISYMC)
            INUM = IL1ZETA(LABELA,LORXA,-FREQA,ISYMA)
            INUM = IL1ZETA(LABELB,LORXB,-FREQB,ISYMB)
            INUM = IL1ZETA(LABELC,LORXC,-FREQC,ISYMC)

*           second-order amplitude (R2) vectors:
            IF (USE_R2) THEN
              INUM = IR2TAMP(LABELA,.FALSE.,+FREQA,ISYMA,
     &                       LABELB,.FALSE.,+FREQB,ISYMB)
              INUM = IR2TAMP(LABELA,.FALSE.,+FREQA,ISYMA,
     &                       LABELC,.FALSE.,+FREQC,ISYMC)
              INUM = IR2TAMP(LABELB,.FALSE.,+FREQB,ISYMB,
     &                       LABELC,.FALSE.,+FREQC,ISYMC)
              INUM = IR2TAMP(LABELA,.FALSE.,-FREQA,ISYMA,
     &                       LABELB,.FALSE.,-FREQB,ISYMB)
              INUM = IR2TAMP(LABELA,.FALSE.,-FREQA,ISYMA,
     &                       LABELC,.FALSE.,-FREQC,ISYMC)
              INUM = IR2TAMP(LABELB,.FALSE.,-FREQB,ISYMB,
     &                       LABELC,.FALSE.,-FREQC,ISYMC)
            END IF

          END DO


         IF (.NOT. LRELAX) THEN

*         for dispersion coefficients:
*         ----------------------------
*         for T(0) = RC(0) vectors is taken care of seperately,
*         they should not be put to the Cauchy vector list, before
*         the equations for the T vectors have been solved.

          DO IDISP = 1, NQRDISP
            ICA  = IQCAUA(IDISP)
            ICB  = IQCAUB(IDISP)
            ICC  = IQCAUC(IDISP)

            ICTOT  = ICA + ICB + ICC
            ISACAU = 2*( (ICTOT/2)*2 - ICTOT ) + 1

            IF (ISACAU.EQ.ISAPROP .OR. ISAPROP.EQ.0 .OR. ALLDSPCF) THEN

*             request first-order right Cauchy vectors:

              IF (ICA.GT.0) INUM = ILRCAMP(LABELA,ICA,ISYMA)
              IF (ICB.GT.0) INUM = ILRCAMP(LABELB,ICB,ISYMB)
              IF (ICC.GT.0) INUM = ILRCAMP(LABELC,ICC,ISYMC)

*             request first order left Cauchy vectors:

              IF (ICA.GT.0) INUM = ILC1AMP(LABELA,ICA,ISYMA)
              IF (ICB.GT.0) INUM = ILC1AMP(LABELB,ICB,ISYMB)
              IF (ICC.GT.0) INUM = ILC1AMP(LABELC,ICC,ISYMC)

*             second-order right Cauchy (CR2) vectors:
              IF (USE_R2) THEN
                IF ((ICA+ICB).GT.0)
     &              INUM = ICR2AMP(LABELA,ICA,ISYMA,LABELB,ICB,ISYMB)
                IF ((ICA+ICC).GT.0)
     &              INUM = ICR2AMP(LABELA,ICA,ISYMA,LABELC,ICC,ISYMC)
                IF ((ICB+ICC).GT.0)
     &              INUM = ICR2AMP(LABELB,ICB,ISYMB,LABELC,ICC,ISYMC)
              END IF

            END IF

          END DO

         END IF

        END IF

      END DO


      RETURN
      END
*---------------------------------------------------------------------*
c /* deck cc_crind */
*=====================================================================*
       SUBROUTINE CC_CRIND
*---------------------------------------------------------------------*
*
*    Purpose: Determine which response t amplitudes and zeta
*             multipliers required for the second hyperpolarizabilities
*             and their dispersion coefficients
*
*    Written by Christof Haettig, October 1996.
*    Dispersion coefficients Februar 1998 (Christof Haettig).
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <cccrinf.h>
#include <ccrspprp.h>
#include <ccroper.h>
#include <cccperm.h>
#include <ccsdinp.h>

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

* variables:
      CHARACTER*8 LABELA, LABELB, LABELC, LABELD
      CHARACTER*8 LABEL1, LABEL2, LABEL3, LABEL4
      INTEGER ISYMB, ISYMC, ISYMA, ISYMD, IFREQ, INUM, IOPER, IDX
      INTEGER ICAUA, ICAUB, ICAUC, ICAUD, IDISP, ISYM1, ISYM2
      INTEGER ICAU1, ICAU2, ICAU3, ICAU4, ISYM3, ISYM4, P
      
      REAL*8  FREQA, FREQB, FREQC, FREQD

* external functions:
      INTEGER IR2TAMP
      INTEGER IL2ZETA
      INTEGER IR1TAMP
      INTEGER IL1ZETA
      INTEGER IROPER
      INTEGER ICHI2
      INTEGER IRHSR3
      INTEGER IRHSR2
      INTEGER IR3TAMP
      INTEGER ILRCAMP
      INTEGER ILC1AMP
      INTEGER ICR2AMP
      INTEGER ICL2AMP
      INTEGER IETACL2
      INTEGER IRHSCR2

* data:
      LOGICAL FIRSTCALL
      SAVE    FIRSTCALL
      DATA    FIRSTCALL /.TRUE./


      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CC_CRIND> NCROPER = ',NCROPER
      END IF

*---------------------------------------------------------------------*
* test if operators are available and translate IACROP, IBCROP, ICCROP
* and IDCROP arrays from the PRPLBL list to the new list maintained 
* by IROPER.
*---------------------------------------------------------------------*
      IF (FIRSTCALL) THEN

       IOPER = 1
       DO WHILE (IOPER .LE. NCROPER)

        LABELA = PRPLBL(IACROP(IOPER))
        LABELB = PRPLBL(IBCROP(IOPER))
        LABELC = PRPLBL(ICCROP(IOPER))
        LABELD = PRPLBL(IDCROP(IOPER))

        IF (      (IROPER(LABELA,ISYMA) .LT. 0)
     &       .OR. (IROPER(LABELB,ISYMB) .LT. 0)
     &       .OR. (IROPER(LABELC,ISYMC) .LT. 0) 
     &       .OR. (IROPER(LABELD,ISYMD) .LT. 0) ) THEN

          WRITE(LUPRI,'(/2X,9A,/2X,2A)')
     &     ' WARNING: ONE OF THE OPERATORS WITH THE LABELS "',
     &     LABELA,'", "', LABELB,'", "', LABELC,'", "',LABELD,
     &     '" IS NOT AVAILABLE.',
     &     ' HYPERPOLARIZABILITY CALCULATION IS CANCELED FOR THIS',
     &     ' OPERATOR QUADRUPLE.'

          DO IDX = IOPER, NCROPER-1
            IACROP(IDX) = IACROP(IDX+1)
            IBCROP(IDX) = IBCROP(IDX+1)
            ICCROP(IDX) = ICCROP(IDX+1)
            IDCROP(IDX) = IDCROP(IDX+1)
          END DO

          NCROPER = NCROPER - 1

        ELSE 
          IACROP(IOPER) = IROPER(LABELA,ISYMA) 
          IBCROP(IOPER) = IROPER(LABELB,ISYMB)
          ICCROP(IOPER) = IROPER(LABELC,ISYMC)
          IDCROP(IOPER) = IROPER(LABELD,ISYMD)

          IOPER = IOPER + 1
        END IF

       END DO

       FIRSTCALL = .FALSE.

      END IF ! (FIRSTCALL)

*---------------------------------------------------------------------*
* set list entries for the required response vectors:
*---------------------------------------------------------------------*
      IF (CC3) THEN
        WRITE(LUPRI,'(/5x,A/)')'Prepare CC3 cubic response calculation.'
        IF (USE_LBCD) THEN
          USE_LBCD = .FALSE.
          WRITE(LUPRI,'(/5x,A//)') 
     &       'USE_LBCD flag (.L2 BCD) is switched off for CC3.'
        END IF
        IF (USE_L2BC) THEN
          USE_L2BC = .FALSE.
          WRITE(LUPRI,'(/5x,A//)') 
     &       'USE_L2BC flag (.L2 BC ) is switched off for CC3.'
        END IF
        IF (L_USE_CHI2) THEN
          L_USE_CHI2 = .FALSE.
          WRITE(LUPRI,'(/5x,A//)') 
     &       'L_USE_CHI2 flag (.USECHI) is switched off for CC3.'
        END IF
        IF (L_USE_XKS3) THEN
          L_USE_XKS3 = .FALSE.
          WRITE(LUPRI,'(/5x,A//)') 
     &       'L_USE_XKS3 flag (.USEXKS) is switched off for CC3.'
        END IF
        IF (NCRDISP.GT.0) THEN
          NCRDISP = 0
          WRITE(LUPRI,'(/5x,A//)') 
     &    'Dispersion coefficients (.DISPCF) are switched off for CC3.'
        END IF
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'USE_L2BC:',USE_L2BC
        WRITE (LUPRI,*) 'USE_LBCD:',USE_LBCD
        IF (USE_LBCD) THEN
          WRITE (LUPRI,*) 'use L2(BC),L2(BD),L2(CD) to eliminate the'
          WRITE (LUPRI,*) 'R2(AD),R2(AC),R2(AB) vectors...'
        ELSE IF (USE_L2BC) THEN
          WRITE (LUPRI,*) 'use L2(BC) to eliminate R2(AD)...'
        ELSE
          WRITE (LUPRI,*) 'use symmetric 2n+1/2n+2 rule formula...'
        END IF
      END IF

      DO IOPER = 1, NCROPER
        LABELA = LBLOPR(IACROP(IOPER))
        LABELB = LBLOPR(IBCROP(IOPER))
        LABELC = LBLOPR(ICCROP(IOPER))
        LABELD = LBLOPR(IDCROP(IOPER))

        ISYMA  = ISYOPR(IACROP(IOPER))
        ISYMB  = ISYOPR(IBCROP(IOPER))
        ISYMC  = ISYOPR(ICCROP(IOPER))
        ISYMD  = ISYOPR(IDCROP(IOPER))

        
        IF (MULD2H(ISYMA,ISYMB).EQ.MULD2H(ISYMC,ISYMD)) THEN
      
*         for frequency-dependent hyperpolarizabilities:
*         ----------------------------------------------
          DO IFREQ = 1, NCRFREQ
            FREQA  = ACRFR(IFREQ)
            FREQB  = BCRFR(IFREQ)
            FREQC  = CCRFR(IFREQ)
            FREQD  = DCRFR(IFREQ)

            IF (LOCDBG) THEN
              WRITE (LUPRI,*) 'CC_CRIND> put on the list:',
     &          LABELA,'(',FREQA,'),  ', LABELB,'(',FREQB,'),  ',
     &          LABELC,'(',FREQC,'),  ', LABELD,'(',FREQD,')'
            END IF

*           request second-order l and t response vectors:

            IF (USE_LBCD) THEN
              INUM=IL2ZETA(LABELB,+FREQB,ISYMB,LABELC,+FREQC,ISYMC)!B,C
              INUM=IL2ZETA(LABELB,+FREQB,ISYMB,LABELD,+FREQD,ISYMD)!B,D
              INUM=IL2ZETA(LABELC,+FREQC,ISYMC,LABELD,+FREQD,ISYMD)!C,D

              INUM=IL2ZETA(LABELB,-FREQB,ISYMB,LABELC,-FREQC,ISYMC)!B,C
              INUM=IL2ZETA(LABELB,-FREQB,ISYMB,LABELD,-FREQD,ISYMD)!B,D
              INUM=IL2ZETA(LABELC,-FREQC,ISYMC,LABELD,-FREQD,ISYMD)!C,D

              INUM=IRHSR2(LABELA,.FALSE.,+FREQA,ISYMA,
     &                    LABELB,.FALSE.,+FREQB,ISYMB)!A,B
              INUM=IRHSR2(LABELA,.FALSE.,+FREQA,ISYMA,
     &                    LABELC,.FALSE.,+FREQC,ISYMC)!A,C
              INUM=IRHSR2(LABELA,.FALSE.,+FREQA,ISYMA,
     &                    LABELD,.FALSE.,+FREQD,ISYMD)!A,D

              INUM=IRHSR2(LABELA,.FALSE.,-FREQA,ISYMA,
     &                    LABELB,.FALSE.,-FREQB,ISYMB)!A,B
              INUM=IRHSR2(LABELA,.FALSE.,-FREQA,ISYMA,
     &                    LABELC,.FALSE.,-FREQC,ISYMC)!A,C
              INUM=IRHSR2(LABELA,.FALSE.,-FREQA,ISYMA,
     &                    LABELD,.FALSE.,-FREQD,ISYMD)!A,D
            ELSE IF (USE_L2BC) THEN
              INUM=IL2ZETA(LABELB,       +FREQB,ISYMB,
     &                     LABELC,       +FREQC,ISYMC)!B,C
              INUM=IL2ZETA(LABELB,       -FREQB,ISYMB,
     &                     LABELC,       -FREQC,ISYMC)!B,C

              INUM=IR2TAMP(LABELA,.FALSE.,+FREQA,ISYMA,
     &                     LABELB,.FALSE.,+FREQB,ISYMB)!A,B
              INUM=IR2TAMP(LABELA,.FALSE.,+FREQA,ISYMA,
     &                     LABELC,.FALSE.,+FREQC,ISYMC)!A,C

              INUM=IR2TAMP(LABELA,.FALSE.,-FREQA,ISYMA,
     &                     LABELB,.FALSE.,-FREQB,ISYMB)!A,B
              INUM=IR2TAMP(LABELA,.FALSE.,-FREQA,ISYMA,
     &                     LABELC,.FALSE.,-FREQC,ISYMC)!A,C

              INUM=IRHSR2(LABELA,.FALSE.,+FREQA,ISYMA,
     &                    LABELD,.FALSE.,+FREQD,ISYMD)!A,D
              INUM=IRHSR2(LABELA,.FALSE.,-FREQA,ISYMA,
     &                    LABELD,.FALSE.,-FREQD,ISYMD)!A,D
            ELSE 
              INUM=IR2TAMP(LABELA,.FALSE.,+FREQA,ISYMA,
     &                     LABELB,.FALSE.,+FREQB,ISYMB)!A,B
              INUM=IR2TAMP(LABELA,.FALSE.,+FREQA,ISYMA,
     &                     LABELC,.FALSE.,+FREQC,ISYMC)!A,C
              INUM=IR2TAMP(LABELA,.FALSE.,+FREQA,ISYMA,
     &                     LABELD,.FALSE.,+FREQD,ISYMD)!A,D

              INUM=IR2TAMP(LABELA,.FALSE.,-FREQA,ISYMA,
     &                     LABELB,.FALSE.,-FREQB,ISYMB)!A,B
              INUM=IR2TAMP(LABELA,.FALSE.,-FREQA,ISYMA,
     &                     LABELC,.FALSE.,-FREQC,ISYMC)!A,C
              INUM=IR2TAMP(LABELA,.FALSE.,-FREQA,ISYMA,
     &                     LABELD,.FALSE.,-FREQD,ISYMD)!A,D
            END IF

              INUM=IR2TAMP(LABELB,.FALSE.,+FREQB,ISYMB,
     &                     LABELC,.FALSE.,+FREQC,ISYMC)!B,C
              INUM=IR2TAMP(LABELB,.FALSE.,+FREQB,ISYMB,
     &                     LABELD,.FALSE.,+FREQD,ISYMD)!B,D
              INUM=IR2TAMP(LABELC,.FALSE.,+FREQC,ISYMC,
     &                     LABELD,.FALSE.,+FREQD,ISYMD)!C,D

              INUM=IR2TAMP(LABELB,.FALSE.,-FREQB,ISYMB,
     &                     LABELC,.FALSE.,-FREQC,ISYMC)!B,C
              INUM=IR2TAMP(LABELB,.FALSE.,-FREQB,ISYMB,
     &                     LABELD,.FALSE.,-FREQD,ISYMD)!B,D
              INUM=IR2TAMP(LABELC,.FALSE.,-FREQC,ISYMC,
     &                     LABELD,.FALSE.,-FREQD,ISYMD)!C,D

*           request second-order chi vectors:

            IF (L_USE_CHI2) THEN
c             INUM = ICHI2(LABELA,.FALSE.,+FREQA,ISYMA,
c    &                     LABELB,.FALSE.,+FREQB,ISYMB)!A,B
c             INUM = ICHI2(LABELA,.FALSE.,+FREQA,ISYMA,
c    &                     LABELC,.FALSE.,+FREQC,ISYMC)!A,C
c             INUM = ICHI2(LABELA,.FALSE.,+FREQA,ISYMA,
c    &                     LABELD,.FALSE.,+FREQD,ISYMD)!A,D
c             INUM = ICHI2(LABELB,.FALSE.,+FREQB,ISYMB,
c    &                     LABELC,.FALSE.,+FREQC,ISYMC)!B,C
c             INUM = ICHI2(LABELB,.FALSE.,+FREQB,ISYMB,
c    &                     LABELD,.FALSE.,+FREQD,ISYMD)!B,D
c             INUM = ICHI2(LABELC,.FALSE.,+FREQC,ISYMC,
c    &                     LABELD,.FALSE.,+FREQD,ISYMD)!C,D

c             INUM = ICHI2(LABELA,.FALSE.,-FREQA,ISYMA,
c    &                     LABELB,.FALSE.,-FREQB,ISYMB)!A,B
c             INUM = ICHI2(LABELA,.FALSE.,-FREQA,ISYMA,
c    &                     LABELC,.FALSE.,-FREQC,ISYMC)!A,C
c             INUM = ICHI2(LABELA,.FALSE.,-FREQA,ISYMA,
c    &                     LABELD,.FALSE.,-FREQD,ISYMD)!A,D
c             INUM = ICHI2(LABELB,.FALSE.,-FREQB,ISYMB,
c    &                     LABELC,.FALSE.,-FREQC,ISYMC)!B,C
c             INUM = ICHI2(LABELB,.FALSE.,-FREQB,ISYMB,
c    &                     LABELD,.FALSE.,-FREQD,ISYMD)!B,D
c             INUM = ICHI2(LABELC,.FALSE.,-FREQC,ISYMC,
c    &                     LABELD,.FALSE.,-FREQD,ISYMD)!C,D

              INUM = IL2ZETA(LABELA,+FREQA,ISYMA,
     &                       LABELB,+FREQB,ISYMB)!A,B
              INUM = IL2ZETA(LABELA,+FREQA,ISYMA,
     &                       LABELC,+FREQC,ISYMC)!A,C
              INUM = IL2ZETA(LABELA,+FREQA,ISYMA,
     &                       LABELD,+FREQD,ISYMD)!A,D
              INUM = IL2ZETA(LABELB,+FREQB,ISYMB,
     &                       LABELC,+FREQC,ISYMC)!B,C
              INUM = IL2ZETA(LABELB,+FREQB,ISYMB,
     &                       LABELD,+FREQD,ISYMD)!B,D
              INUM = IL2ZETA(LABELC,+FREQC,ISYMC,
     &                       LABELD,+FREQD,ISYMD)!C,D

              INUM = IL2ZETA(LABELA,-FREQA,ISYMA,
     &                       LABELB,-FREQB,ISYMB)!A,B
              INUM = IL2ZETA(LABELA,-FREQA,ISYMA,
     &                       LABELC,-FREQC,ISYMC)!A,C
              INUM = IL2ZETA(LABELA,-FREQA,ISYMA,
     &                       LABELD,-FREQD,ISYMD)!A,D
              INUM = IL2ZETA(LABELB,-FREQB,ISYMB,
     &                       LABELC,-FREQC,ISYMC)!B,C
              INUM = IL2ZETA(LABELB,-FREQB,ISYMB,
     &                       LABELD,-FREQD,ISYMD)!B,D
              INUM = IL2ZETA(LABELC,-FREQC,ISYMC,
     &                       LABELD,-FREQD,ISYMD)!C,D
            END IF

*           request third-order amplitude rhs vectors:

            IF (L_USE_XKS3) THEN
              INUM = IR3TAMP(LABELA,+FREQA,ISYMA, LABELB,+FREQB,ISYMB,
     &                      LABELC,+FREQC,ISYMC) ! A,B,C
              INUM = IR3TAMP(LABELA,+FREQA,ISYMA, LABELB,+FREQB,ISYMB,
     &                      LABELD,+FREQD,ISYMD) ! A,B,D
              INUM = IR3TAMP(LABELA,+FREQA,ISYMA, LABELC,+FREQC,ISYMC,
     &                      LABELD,+FREQD,ISYMD) ! A,C,D
              INUM = IR3TAMP(LABELB,+FREQB,ISYMB, LABELC,+FREQC,ISYMC,
     &                      LABELD,+FREQD,ISYMD) ! B,C,D

              INUM = IR3TAMP(LABELA,-FREQA,ISYMA, LABELB,-FREQB,ISYMB,
     &                      LABELC,-FREQC,ISYMC) ! A,B,C
              INUM = IR3TAMP(LABELA,-FREQA,ISYMA, LABELB,-FREQB,ISYMB,
     &                      LABELD,-FREQD,ISYMD) ! A,B,D
              INUM = IR3TAMP(LABELA,-FREQA,ISYMA, LABELC,-FREQC,ISYMC,
     &                      LABELD,-FREQD,ISYMD) ! A,C,D
              INUM = IR3TAMP(LABELB,-FREQB,ISYMB, LABELC,-FREQC,ISYMC,
     &                      LABELD,-FREQD,ISYMD) ! B,C,D
            END IF

*           request (unrelaxed) first-order t response vectors:

            INUM = IR1TAMP(LABELA,.FALSE.,+FREQA,ISYMA)
            INUM = IR1TAMP(LABELB,.FALSE.,+FREQB,ISYMB)
            INUM = IR1TAMP(LABELC,.FALSE.,+FREQC,ISYMC)
            INUM = IR1TAMP(LABELD,.FALSE.,+FREQD,ISYMD)
            INUM = IR1TAMP(LABELA,.FALSE.,-FREQA,ISYMA)
            INUM = IR1TAMP(LABELB,.FALSE.,-FREQB,ISYMB)
            INUM = IR1TAMP(LABELC,.FALSE.,-FREQC,ISYMC)
            INUM = IR1TAMP(LABELD,.FALSE.,-FREQD,ISYMD)


*           request first order zeta response vectors:

            INUM = IL1ZETA(LABELA,.FALSE.,+FREQA,ISYMA)
            INUM = IL1ZETA(LABELB,.FALSE.,+FREQB,ISYMB)
            INUM = IL1ZETA(LABELC,.FALSE.,+FREQC,ISYMC)
            INUM = IL1ZETA(LABELD,.FALSE.,+FREQD,ISYMD)
            INUM = IL1ZETA(LABELA,.FALSE.,-FREQA,ISYMA)
            INUM = IL1ZETA(LABELB,.FALSE.,-FREQB,ISYMB)
            INUM = IL1ZETA(LABELC,.FALSE.,-FREQC,ISYMC)
            INUM = IL1ZETA(LABELD,.FALSE.,-FREQD,ISYMD)
          END DO

*         for dispersion coefficients :
*         -----------------------------
*         RC(0), CR2(0,0) and CL2(0,0) vectors are calculated
*         as R1(0), R2(0,0) and L2(0,0)... the identification
*         with the cauchy vectors with these response vectors it
*         ensured later by the CC_RDRSP routine... but here we a
*         must no put them on the Cauchy lists but on the response
*         vector lists...

          DO IDISP = 1, NCRDISP
            ICAUA = ICCAUA(IDISP)
            ICAUB = ICCAUB(IDISP)
            ICAUC = ICCAUC(IDISP)
            ICAUD = ICCAUD(IDISP)

*           request first-order right Cauchy vectors:

            IF (ICAUA.GT.0) INUM = ILRCAMP(LABELA,ICAUA,ISYMA)
            IF (ICAUA.EQ.0) INUM = IR1TAMP(LABELA,.FALSE.,0.0d0,ISYMA)

            IF (ICAUB.GT.0) INUM = ILRCAMP(LABELB,ICAUB,ISYMB)
            IF (ICAUB.EQ.0) INUM = IR1TAMP(LABELB,.FALSE.,0.0d0,ISYMB)

            IF (ICAUC.GT.0) INUM = ILRCAMP(LABELC,ICAUC,ISYMC)
            IF (ICAUC.EQ.0) INUM = IR1TAMP(LABELC,.FALSE.,0.0d0,ISYMC)

            IF (ICAUD.GT.0) INUM = ILRCAMP(LABELD,ICAUD,ISYMD)
            IF (ICAUD.EQ.0) INUM = IR1TAMP(LABELD,.FALSE.,0.0d0,ISYMD)

*           request first-order left Cauchy vectors:

            IF (ICAUA.GT.0) INUM = ILC1AMP(LABELA,ICAUA,ISYMA)
            IF (ICAUA.EQ.0) INUM = IL1ZETA(LABELA,.FALSE.,0.0d0,ISYMA)

            IF (ICAUB.GT.0) INUM = ILC1AMP(LABELB,ICAUB,ISYMB)
            IF (ICAUB.EQ.0) INUM = IL1ZETA(LABELB,.FALSE.,0.0d0,ISYMB)

            IF (ICAUC.GT.0) INUM = ILC1AMP(LABELC,ICAUC,ISYMC)
            IF (ICAUC.EQ.0) INUM = IL1ZETA(LABELC,.FALSE.,0.0d0,ISYMC)

            IF (ICAUD.GT.0) INUM = ILC1AMP(LABELD,ICAUD,ISYMD)
            IF (ICAUD.EQ.0) INUM = IL1ZETA(LABELD,.FALSE.,0.0d0,ISYMD)


*           request second-order right Cauchy vectors:

            IF ( NO_2NP1_RULE ) THEN
*             ... if we do not use the 2N+1 rule for the second-order
*                 Cauchy intermediates, we need for all pair of 
*                 operator and accompanied Cauchy order the 
*                 second-order amplitude Cauchy vectors "CR2"

*             .... (A,B) pair ...
              IF (ICAUA.GT.0 .OR. ICAUB.GT.0) THEN 
                INUM = ICR2AMP(LABELA,ICAUA,ISYMA,LABELB,ICAUB,ISYMB)
              ELSE
                INUM = IR2TAMP(LABELA,.FALSE.,0.0d0,ISYMA,
     &                         LABELB,.FALSE.,0.0d0,ISYMB)
              END IF

*             .... (A,C) pair ...
              IF (ICAUA.GT.0 .OR. ICAUC.GT.0) THEN 
                INUM = ICR2AMP(LABELA,ICAUA,ISYMA,LABELC,ICAUC,ISYMC)
              ELSE
                INUM = IR2TAMP(LABELA,.FALSE.,0.0d0,ISYMA,
     &                         LABELC,.FALSE.,0.0d0,ISYMC)
              END IF

*             .... (A,D) pair ...
              IF (ICAUA.GT.0 .OR. ICAUD.GT.0) THEN 
                INUM = ICR2AMP(LABELA,ICAUA,ISYMA,LABELD,ICAUD,ISYMD)
              ELSE
                INUM = IR2TAMP(LABELA,.FALSE.,0.0d0,ISYMA,
     &                         LABELD,.FALSE.,0.0d0,ISYMD)
              END IF

*             .... (B,C) pair ...
              IF (ICAUB.GT.0 .OR. ICAUC.GT.0) THEN
                INUM = ICR2AMP(LABELB,ICAUB,ISYMB,LABELC,ICAUC,ISYMC)
              ELSE
                INUM = IR2TAMP(LABELB,.FALSE.,0.0d0,ISYMB,
     &                         LABELC,.FALSE.,0.0d0,ISYMC)
              END IF

*             .... (B,D) pair ...
              IF (ICAUB.GT.0 .OR. ICAUD.GT.0) THEN
                INUM = ICR2AMP(LABELB,ICAUB,ISYMB,LABELD,ICAUD,ISYMD)
              ELSE
                INUM = IR2TAMP(LABELB,.FALSE.,0.0d0,ISYMB,
     &                         LABELD,.FALSE.,0.0d0,ISYMD)
              END IF

*             .... (C,D) pair ...
              IF (ICAUC.GT.0 .OR. ICAUD.GT.0) THEN
                INUM = ICR2AMP(LABELC,ICAUC,ISYMC,LABELD,ICAUD,ISYMD)
              ELSE
                INUM = IR2TAMP(LABELC,.FALSE.,0.0d0,ISYMC,
     &                         LABELD,.FALSE.,0.0d0,ISYMD)
              END IF

            ELSE
*             ... if we use the 2n+1/2n+2 rules for the second-order
*                 Cauchy intermediates we have more sophisticated
*                 settings with a three-fold case switch for each of
*                 the three different couples of pairs
*                 [(A,B)/(C,D)],  [(A,D)/(B,C)]  and [(A,C)/(D,B)]

              DO P = 1, 3
                LABEL1 = LBLOPR(ICROP(IOPER,I1(P)))
                LABEL2 = LBLOPR(ICROP(IOPER,I2(P)))
                LABEL3 = LBLOPR(ICROP(IOPER,I3(P)))
                LABEL4 = LBLOPR(ICROP(IOPER,I4(P)))
                ICAU1  = ICCAU(IDISP,I1(P))
                ICAU2  = ICCAU(IDISP,I2(P))
                ICAU3  = ICCAU(IDISP,I3(P))
                ICAU4  = ICCAU(IDISP,I4(P))
                ISYM1  = ISYOPR(ICROP(IOPER,I1(P)))
                ISYM2  = ISYOPR(ICROP(IOPER,I2(P)))
                ISYM3  = ISYOPR(ICROP(IOPER,I3(P)))
                ISYM4  = ISYOPR(ICROP(IOPER,I4(P)))

                IF      ( (ICAU1+ICAU2) .GT. (ICAU3+ICAU4) )THEN
                 INUM = IETACL2(LABEL1,ICAU1,ISYM1,LABEL2,ICAU2,ISYM2)
                 INUM = IRHSCR2(LABEL1,ICAU1,ISYM1,LABEL2,ICAU2,ISYM2)
                 IF ( (ICAU3+ICAU4).GT.0 ) THEN
                   INUM=ICL2AMP(LABEL3,ICAU3,ISYM3,LABEL4,ICAU4,ISYM4)
                   INUM=ICR2AMP(LABEL3,ICAU3,ISYM3,LABEL4,ICAU4,ISYM4)
                 ELSE
                   INUM=IL2ZETA(LABEL3,0.0d0,ISYM3,LABEL4,0.0d0,ISYM4)
                   INUM=IR2TAMP(LABEL3,.FALSE.,0.0d0,ISYM3,
     &                          LABEL4,.FALSE.,0.0d0,ISYM4)
                 END IF
                ELSE IF ( (ICAU1+ICAU2) .EQ. (ICAU3+ICAU4) )THEN

                 IF ( (ICAU1+ICAU2).GT.0 ) THEN
                   INUM=IETACL2(LABEL1,ICAU1,ISYM1,LABEL2,ICAU2,ISYM2)
                   INUM=ICR2AMP(LABEL1,ICAU1,ISYM1,LABEL2,ICAU2,ISYM2)
                 ELSE
                   INUM=ICHI2(  LABEL1,.FALSE.,0.0d0,ISYM1,
     &                          LABEL2,.FALSE.,0.0d0,ISYM2)
                   INUM=IR2TAMP(LABEL1,.FALSE.,0.0d0,ISYM1,
     &                          LABEL2,.FALSE.,0.0d0,ISYM2)
                 END IF
                 IF ( (ICAU3+ICAU4).GT.0 ) THEN
                   INUM=IETACL2(LABEL3,ICAU3,ISYM3,LABEL4,ICAU4,ISYM4)
                   INUM=ICR2AMP(LABEL3,ICAU3,ISYM3,LABEL4,ICAU4,ISYM4)
                 ELSE
                   INUM=ICHI2(  LABEL3,.FALSE.,0.0d0,ISYM3,
     &                          LABEL4,.FALSE.,0.0d0,ISYM4)
                   INUM=IR2TAMP(LABEL3,.FALSE.,0.0d0,ISYM3,
     &                          LABEL4,.FALSE.,0.0d0,ISYM4)
                 END IF
                 IF      (ICAU1.EQ.1 .AND. ICAU2.EQ.0) THEN
                  INUM=IL2ZETA(LABEL1,0.0d0,ISYM1,LABEL2,0.0d0,ISYM2)
                 ELSE IF (ICAU1.GT.0                 ) THEN
                  INUM=ICL2AMP(LABEL1,ICAU1-1,ISYM1,LABEL2,ICAU2,ISYM2)
                 END IF
                 IF      (ICAU2.EQ.1 .AND. ICAU1.EQ.0) THEN
                  INUM=IL2ZETA(LABEL1,0.0d0,ISYM1,LABEL2,0.0d0,ISYM2)
                 ELSE IF (ICAU2.GT.0                 ) THEN
                  INUM=ICL2AMP(LABEL1,ICAU1,ISYM1,LABEL2,ICAU2-1,ISYM2)
                 END IF
                 IF      (ICAU3.EQ.1 .AND. ICAU4.EQ.0) THEN
                  INUM=IL2ZETA(LABEL3,0.0d0,ISYM3,LABEL4,0.0d0,ISYM4)
                 ELSE IF (ICAU3.GT.0                 ) THEN
                  INUM=ICL2AMP(LABEL3,ICAU3-1,ISYM3,LABEL4,ICAU4,ISYM4)
                 END IF
                 IF      (ICAU4.EQ.1 .AND. ICAU3.EQ.0) THEN
                  INUM=IL2ZETA(LABEL3,0.0d0,ISYM3,LABEL4,0.0d0,ISYM4)
                 ELSE IF (ICAU4.GT.0                 ) THEN
                  INUM=ICL2AMP(LABEL3,ICAU3,ISYM3,LABEL4,ICAU4-1,ISYM4)
                 END IF

                ELSE IF ( (ICAU1+ICAU2) .LT. (ICAU3+ICAU4) )THEN

                 INUM = IETACL2(LABEL3,ICAU3,ISYM3,LABEL4,ICAU4,ISYM4)
                 INUM = IRHSCR2(LABEL3,ICAU3,ISYM3,LABEL4,ICAU4,ISYM4)
                 IF ( (ICAU1+ICAU2).GT.0 ) THEN
                   INUM=ICL2AMP(LABEL1,ICAU1,ISYM1,LABEL2,ICAU2,ISYM2)
                   INUM=ICR2AMP(LABEL1,ICAU1,ISYM1,LABEL2,ICAU2,ISYM2)
                 ELSE
                   INUM=IL2ZETA(LABEL1,0.0d0,ISYM1,LABEL2,0.0d0,ISYM2)
                   INUM=IR2TAMP(LABEL1,.FALSE.,0.0d0,ISYM1,
     &                          LABEL2,.FALSE.,0.0d0,ISYM2)
                 END IF

                END IF

              END DO ! IPAIRS

            END IF  ! (NO_2NP1_RULE)

          END DO
        END IF

      END DO


      RETURN
      END
*---------------------------------------------------------------------*
c /* deck CC_4RIND */
*=====================================================================*
       SUBROUTINE CC_4RIND
*---------------------------------------------------------------------*
*
*    Purpose: Determine which response t amplitudes and zeta
*             multipliers required for the third hyperpolarizabilities
*
*    Written by Christof Haettig, April 1997.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <cc4rinf.h>
#include <ccrspprp.h>
#include <ccroper.h>
#include <ccl2rsp.h>

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

* variables:
      CHARACTER*8 LABELA, LABELB, LABELC, LABELD, LABELE
      INTEGER ISYMB, ISYMC, ISYMA, ISYMD, ISYME
      INTEGER IFREQ, INUM, IOPER, IDX
      
      REAL*8  FREQA, FREQB, FREQC, FREQD, FREQE

* external functions:
      INTEGER IR2TAMP
      INTEGER IR1TAMP
      INTEGER IL1ZETA
      INTEGER IL2ZETA
      INTEGER IROPER
      INTEGER ICHI3

* data:
      LOGICAL FIRSTCALL
      SAVE    FIRSTCALL
      DATA    FIRSTCALL /.TRUE./


      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CC_4RIND> N4ROPER = ',N4ROPER
        WRITE (LUPRI,*) 'LL2OPN:',LL2OPN
      END IF

*---------------------------------------------------------------------*
* test if operators are available and translate IA4ROP, IB4ROP, etc.
* arrays from the PRPLBL list to the new list maintained by IROPER.
*---------------------------------------------------------------------*
      IF (FIRSTCALL) THEN

       IOPER = 1
       DO WHILE (IOPER .LE. N4ROPER)

        LABELA = PRPLBL(IA4ROP(IOPER))
        LABELB = PRPLBL(IB4ROP(IOPER))
        LABELC = PRPLBL(IC4ROP(IOPER))
        LABELD = PRPLBL(ID4ROP(IOPER))
        LABELE = PRPLBL(IE4ROP(IOPER))

        IF (      (IROPER(LABELA,ISYMA) .LT. 0)
     &       .OR. (IROPER(LABELB,ISYMB) .LT. 0)
     &       .OR. (IROPER(LABELC,ISYMC) .LT. 0) 
     &       .OR. (IROPER(LABELD,ISYMD) .LT. 0)
     &       .OR. (IROPER(LABELE,ISYME) .LT. 0) ) THEN

          WRITE(LUPRI,'(/2X,9A,/2X,2A)')
     &     ' WARNING: ONE OF THE OPERATORS WITH THE LABELS "',
     &     LABELA,'", "', LABELB,'", "', LABELC,'", "',LABELD,
     &     '", "',LABELE,'" IS NOT AVAILABLE.',
     &     ' HYPERPOLARIZABILITY CALCULATION IS CANCELED FOR THIS',
     &     ' OPERATOR QUINTUPLE.'

          DO IDX = IOPER, N4ROPER-1
            IA4ROP(IDX) = IA4ROP(IDX+1)
            IB4ROP(IDX) = IB4ROP(IDX+1)
            IC4ROP(IDX) = IC4ROP(IDX+1)
            ID4ROP(IDX) = ID4ROP(IDX+1)
            IE4ROP(IDX) = IE4ROP(IDX+1)
          END DO

          N4ROPER = N4ROPER - 1

        ELSE 
          IA4ROP(IOPER) = IROPER(LABELA,ISYMA) 
          IB4ROP(IOPER) = IROPER(LABELB,ISYMB)
          IC4ROP(IOPER) = IROPER(LABELC,ISYMC)
          ID4ROP(IOPER) = IROPER(LABELD,ISYMD)
          IE4ROP(IOPER) = IROPER(LABELE,ISYME)

          IOPER = IOPER + 1
        END IF

       END DO

       FIRSTCALL = .FALSE.

      END IF ! (FIRSTCALL)

*---------------------------------------------------------------------*
* set list entries for the required response vectors:
*---------------------------------------------------------------------*
      DO IOPER = 1, N4ROPER
        LABELA = LBLOPR(IA4ROP(IOPER))
        LABELB = LBLOPR(IB4ROP(IOPER))
        LABELC = LBLOPR(IC4ROP(IOPER))
        LABELD = LBLOPR(ID4ROP(IOPER))
        LABELE = LBLOPR(IE4ROP(IOPER))

        ISYMA  = ISYOPR(IA4ROP(IOPER))
        ISYMB  = ISYOPR(IB4ROP(IOPER))
        ISYMC  = ISYOPR(IC4ROP(IOPER))
        ISYMD  = ISYOPR(ID4ROP(IOPER))
        ISYME  = ISYOPR(IE4ROP(IOPER))

        
        IF (MULD2H(ISYMA,ISYMB).EQ.MULD2H(MULD2H(ISYMC,ISYMD),ISYME)
     &    ) THEN
      
          DO IFREQ = 1, N4RFREQ
            FREQA  = A4RFR(IFREQ)
            FREQB  = B4RFR(IFREQ)
            FREQC  = C4RFR(IFREQ)
            FREQD  = D4RFR(IFREQ)
            FREQE  = E4RFR(IFREQ)

            IF (LOCDBG) THEN
              WRITE (LUPRI,*) 'CC_4RIND> put on the list:',
     &          LABELA,'(',FREQA,'),  ', LABELB,'(',FREQB,'),  ',
     &          LABELC,'(',FREQC,'),  ', LABELD,'(',FREQD,'),  ',
     &          LABELE,'(',FREQE,')'
            END IF

*           request second-order t response vectors:

            INUM = IR2TAMP(LABELA,.FALSE.,+FREQA,ISYMA,
     &                     LABELB,.FALSE.,+FREQB,ISYMB)!A,B
            INUM = IR2TAMP(LABELA,.FALSE.,+FREQA,ISYMA,
     &                     LABELC,.FALSE.,+FREQC,ISYMC)!A,C
            INUM = IR2TAMP(LABELA,.FALSE.,+FREQA,ISYMA,
     &                     LABELD,.FALSE.,+FREQD,ISYMD)!A,D
            INUM = IR2TAMP(LABELA,.FALSE.,+FREQA,ISYMA,
     &                     LABELE,.FALSE.,+FREQE,ISYME)!A,E
            INUM = IR2TAMP(LABELB,.FALSE.,+FREQB,ISYMB,
     &                     LABELC,.FALSE.,+FREQC,ISYMC)!B,C
            INUM = IR2TAMP(LABELB,.FALSE.,+FREQB,ISYMB,
     &                     LABELD,.FALSE.,+FREQD,ISYMD)!B,D
            INUM = IR2TAMP(LABELB,.FALSE.,+FREQB,ISYMB,
     &                     LABELE,.FALSE.,+FREQE,ISYME)!B,E
            INUM = IR2TAMP(LABELC,.FALSE.,+FREQC,ISYMC,
     &                     LABELD,.FALSE.,+FREQD,ISYMD)!C,D
            INUM = IR2TAMP(LABELC,.FALSE.,+FREQC,ISYMC,
     &                     LABELE,.FALSE.,+FREQE,ISYME)!C,E
            INUM = IR2TAMP(LABELD,.FALSE.,+FREQD,ISYMD,
     &                     LABELE,.FALSE.,+FREQE,ISYME)!D,E

            INUM = IR2TAMP(LABELA,.FALSE.,-FREQA,ISYMA,
     &                     LABELB,.FALSE.,-FREQB,ISYMB)!A,B
            INUM = IR2TAMP(LABELA,.FALSE.,-FREQA,ISYMA,
     &                     LABELC,.FALSE.,-FREQC,ISYMC)!A,C
            INUM = IR2TAMP(LABELA,.FALSE.,-FREQA,ISYMA,
     &                     LABELD,.FALSE.,-FREQD,ISYMD)!A,D
            INUM = IR2TAMP(LABELA,.FALSE.,-FREQA,ISYMA,
     &                     LABELE,.FALSE.,-FREQE,ISYME)!A,E
            INUM = IR2TAMP(LABELB,.FALSE.,-FREQB,ISYMB,
     &                     LABELC,.FALSE.,-FREQC,ISYMC)!B,C
            INUM = IR2TAMP(LABELB,.FALSE.,-FREQB,ISYMB,
     &                     LABELD,.FALSE.,-FREQD,ISYMD)!B,D
            INUM = IR2TAMP(LABELB,.FALSE.,-FREQB,ISYMB,
     &                     LABELE,.FALSE.,-FREQE,ISYME)!B,E
            INUM = IR2TAMP(LABELC,.FALSE.,-FREQC,ISYMC,
     &                     LABELD,.FALSE.,-FREQD,ISYMD)!C,D
            INUM = IR2TAMP(LABELC,.FALSE.,-FREQC,ISYMC,
     &                     LABELE,.FALSE.,-FREQE,ISYME)!C,E
            INUM = IR2TAMP(LABELD,.FALSE.,-FREQD,ISYMD,
     &                     LABELE,.FALSE.,-FREQE,ISYME)!D,E

*           request second-order zeta response vectors:

            INUM = IL2ZETA(LABELA,+FREQA,ISYMA,LABELB,+FREQB,ISYMB)!A,B
            INUM = IL2ZETA(LABELA,+FREQA,ISYMA,LABELC,+FREQC,ISYMC)!A,C
            INUM = IL2ZETA(LABELA,+FREQA,ISYMA,LABELD,+FREQD,ISYMD)!A,D
            INUM = IL2ZETA(LABELA,+FREQA,ISYMA,LABELE,+FREQE,ISYME)!A,E
            INUM = IL2ZETA(LABELB,+FREQB,ISYMB,LABELC,+FREQC,ISYMC)!B,C
            INUM = IL2ZETA(LABELB,+FREQB,ISYMB,LABELD,+FREQD,ISYMD)!B,D
            INUM = IL2ZETA(LABELB,+FREQB,ISYMB,LABELE,+FREQE,ISYME)!B,E
            INUM = IL2ZETA(LABELC,+FREQC,ISYMC,LABELD,+FREQD,ISYMD)!C,D
            INUM = IL2ZETA(LABELC,+FREQC,ISYMC,LABELE,+FREQE,ISYME)!C,E
            INUM = IL2ZETA(LABELD,+FREQD,ISYMD,LABELE,+FREQE,ISYME)!D,E

            INUM = IL2ZETA(LABELA,-FREQA,ISYMA,LABELB,-FREQB,ISYMB)!A,B
            INUM = IL2ZETA(LABELA,-FREQA,ISYMA,LABELC,-FREQC,ISYMC)!A,C
            INUM = IL2ZETA(LABELA,-FREQA,ISYMA,LABELD,-FREQD,ISYMD)!A,D
            INUM = IL2ZETA(LABELA,-FREQA,ISYMA,LABELE,-FREQE,ISYME)!A,E
            INUM = IL2ZETA(LABELB,-FREQB,ISYMB,LABELC,-FREQC,ISYMC)!B,C
            INUM = IL2ZETA(LABELB,-FREQB,ISYMB,LABELD,-FREQD,ISYMD)!B,D
            INUM = IL2ZETA(LABELB,-FREQB,ISYMB,LABELE,-FREQE,ISYME)!B,E
            INUM = IL2ZETA(LABELC,-FREQC,ISYMC,LABELD,-FREQD,ISYMD)!C,D
            INUM = IL2ZETA(LABELC,-FREQC,ISYMC,LABELE,-FREQE,ISYME)!C,E
            INUM = IL2ZETA(LABELD,-FREQD,ISYMD,LABELE,-FREQE,ISYME)!D,E

*           request third-order chi vectors:
            IF (L_USE_CHI3) THEN

              INUM = ICHI3(LABELC,+FREQC,ISYMC, LABELD,+FREQD,ISYMD,
     &                                          LABELE,+FREQE,ISYME)
              INUM = ICHI3(LABELB,+FREQB,ISYMB, LABELD,+FREQD,ISYMD,
     &                                          LABELE,+FREQE,ISYME)
              INUM = ICHI3(LABELB,+FREQB,ISYMB, LABELC,+FREQC,ISYMC,
     &                                          LABELE,+FREQE,ISYME)
              INUM = ICHI3(LABELB,+FREQB,ISYMB, LABELC,+FREQC,ISYMC,
     &                                          LABELD,+FREQD,ISYMD)
              INUM = ICHI3(LABELA,+FREQA,ISYMA, LABELD,+FREQD,ISYMD,
     &                                          LABELE,+FREQE,ISYME)
              INUM = ICHI3(LABELA,+FREQA,ISYMA, LABELC,+FREQC,ISYMC,
     &                                          LABELE,+FREQE,ISYME)
              INUM = ICHI3(LABELA,+FREQA,ISYMA, LABELC,+FREQC,ISYMC,
     &                                          LABELD,+FREQD,ISYMD)
              INUM = ICHI3(LABELA,+FREQA,ISYMA, LABELB,+FREQB,ISYMB,
     &                                          LABELE,+FREQE,ISYME)
              INUM = ICHI3(LABELA,+FREQA,ISYMA, LABELB,+FREQB,ISYMB,
     &                                          LABELD,+FREQD,ISYMD)
              INUM = ICHI3(LABELA,+FREQA,ISYMA, LABELB,+FREQB,ISYMB,
     &                                          LABELC,+FREQC,ISYMC)

              INUM = ICHI3(LABELC,-FREQC,ISYMC, LABELD,-FREQD,ISYMD,
     &                                          LABELE,-FREQE,ISYME)
              INUM = ICHI3(LABELB,-FREQB,ISYMB, LABELD,-FREQD,ISYMD,
     &                                          LABELE,-FREQE,ISYME)
              INUM = ICHI3(LABELB,-FREQB,ISYMB, LABELC,-FREQC,ISYMC,
     &                                          LABELE,-FREQE,ISYME)
              INUM = ICHI3(LABELB,-FREQB,ISYMB, LABELC,-FREQC,ISYMC,
     &                                          LABELD,-FREQD,ISYMD)
              INUM = ICHI3(LABELA,-FREQA,ISYMA, LABELD,-FREQD,ISYMD,
     &                                          LABELE,-FREQE,ISYME)
              INUM = ICHI3(LABELA,-FREQA,ISYMA, LABELC,-FREQC,ISYMC,
     &                                          LABELE,-FREQE,ISYME)
              INUM = ICHI3(LABELA,-FREQA,ISYMA, LABELC,-FREQC,ISYMC,
     &                                          LABELD,-FREQD,ISYMD)
              INUM = ICHI3(LABELA,-FREQA,ISYMA, LABELB,-FREQB,ISYMB,
     &                                          LABELE,-FREQE,ISYME)
              INUM = ICHI3(LABELA,-FREQA,ISYMA, LABELB,-FREQB,ISYMB,
     &                                          LABELD,-FREQD,ISYMD)
              INUM = ICHI3(LABELA,-FREQA,ISYMA, LABELB,-FREQB,ISYMB,
     &                                          LABELC,-FREQC,ISYMC)
  
            END IF


*           request (unrelaxed) first-order t response vectors:

            INUM = IR1TAMP(LABELA,.FALSE.,+FREQA,ISYMA)
            INUM = IR1TAMP(LABELB,.FALSE.,+FREQB,ISYMB)
            INUM = IR1TAMP(LABELC,.FALSE.,+FREQC,ISYMC)
            INUM = IR1TAMP(LABELD,.FALSE.,+FREQD,ISYMD)
            INUM = IR1TAMP(LABELE,.FALSE.,+FREQD,ISYME)

            INUM = IR1TAMP(LABELA,.FALSE.,-FREQA,ISYMA)
            INUM = IR1TAMP(LABELB,.FALSE.,-FREQB,ISYMB)
            INUM = IR1TAMP(LABELC,.FALSE.,-FREQC,ISYMC)
            INUM = IR1TAMP(LABELD,.FALSE.,-FREQD,ISYMD)
            INUM = IR1TAMP(LABELE,.FALSE.,-FREQE,ISYME)


*           request (unrelaxed) first-order zeta response vectors:

            INUM = IL1ZETA(LABELA,.FALSE.,+FREQA,ISYMA)
            INUM = IL1ZETA(LABELB,.FALSE.,+FREQB,ISYMB)
            INUM = IL1ZETA(LABELC,.FALSE.,+FREQC,ISYMC)
            INUM = IL1ZETA(LABELD,.FALSE.,+FREQD,ISYMD)
            INUM = IL1ZETA(LABELE,.FALSE.,+FREQE,ISYME)

            INUM = IL1ZETA(LABELA,.FALSE.,-FREQA,ISYMA)
            INUM = IL1ZETA(LABELB,.FALSE.,-FREQB,ISYMB)
            INUM = IL1ZETA(LABELC,.FALSE.,-FREQC,ISYMC)
            INUM = IL1ZETA(LABELD,.FALSE.,-FREQD,ISYMD)
            INUM = IL1ZETA(LABELE,.FALSE.,-FREQE,ISYME)
          END DO

        END IF

      END DO


      RETURN
      END
*---------------------------------------------------------------------*
c /* deck CC_5RIND */
*=====================================================================*
       SUBROUTINE CC_5RIND
*---------------------------------------------------------------------*
*
*    Purpose: Determine which response t amplitudes and zeta
*             multipliers required for the fourth hyperpolarizabilities
*             (pentic response function)
*
*    Written by Christof Haettig, Maj 1997.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <cc5rinf.h>
#include <cc5perm.h>
#include <ccrspprp.h>
#include <ccroper.h>

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

* variables:
      CHARACTER*8 LABEL(6)
      INTEGER ISYM(6)
      INTEGER IFREQ, INUM, IOPER, IDX, IDXA, IDXB, IDXC, JDX, ISYMTOT
      
      REAL*8  FREQ(6)

* external functions:
      INTEGER IR3TAMP
      INTEGER ICHI3
      INTEGER IROPER

* data:
      LOGICAL FIRSTCALL
      SAVE    FIRSTCALL
      DATA    FIRSTCALL /.TRUE./


*---------------------------------------------------------------------*
* test if operators are available and translate I5ROP array
* from the PRPLBL list to the new list maintained by IROPER.
*---------------------------------------------------------------------*
      IF (FIRSTCALL) THEN

       IOPER = 1
       DO WHILE (IOPER .LE. N5ROPER)

        DO IDX = 1, 6
          LABEL(IDX) = PRPLBL(I5ROP(IOPER,IDX))
        END DO

        IF (      IROPER(LABEL(A),ISYM(A)) .LT. 0
     &       .OR. IROPER(LABEL(B),ISYM(B)) .LT. 0
     &       .OR. IROPER(LABEL(C),ISYM(C)) .LT. 0 
     &       .OR. IROPER(LABEL(D),ISYM(D)) .LT. 0
     &       .OR. IROPER(LABEL(E),ISYM(E)) .LT. 0 
     &       .OR. IROPER(LABEL(F),ISYM(F)) .LT. 0 ) THEN

          WRITE(LUPRI,'(/2X,9A,/2X,2A)')
     &     ' WARNING: ONE OF THE OPERATORS WITH THE LABELS "',
     &     LABEL(A),'", "', LABEL(B),'", "', LABEL(C),'", "',LABEL(D),
     &     '", "',LABEL(E), '", "',LABEL(F),'" IS NOT AVAILABLE.',
     &     ' HYPERPOLARIZABILITY CALCULATION IS CANCELED FOR THIS',
     &     ' OPERATOR HEXTUPLE.'

C         WRITE (LUPRI,*) 'I5ROP:',(I5ROP(IOPER,IDX),IDX=1,6)

          DO JDX = IOPER, N5ROPER-1
          DO IDX = 1, 6
            I5ROP(JDX,IDX) = I5ROP(JDX+1,IDX)
          END DO
          END DO

          N5ROPER = N5ROPER - 1

        ELSE 
          DO IDX = 1, 6
            I5ROP(IOPER,IDX) = IROPER(LABEL(IDX),ISYM(IDX)) 
          END DO

          IOPER = IOPER + 1
        END IF

       END DO

       FIRSTCALL = .FALSE.

      END IF ! (FIRSTCALL)

*---------------------------------------------------------------------*
* set list entries for the required response vectors:
*---------------------------------------------------------------------*
      DO IOPER = 1, N5ROPER
        ISYMTOT = 1
        DO IDX = 1, 6
          LABEL(IDX) = LBLOPR(I5ROP(IOPER,IDX))
          ISYM(IDX)  = ISYOPR(I5ROP(IOPER,IDX))
          ISYMTOT    = MULD2H(ISYMTOT,ISYM(IDX))
        END DO
        
        IF ( ISYMTOT.EQ.1 ) THEN
      
          DO IFREQ = 1, N5RFREQ
            DO IDX = 1, 6
              FREQ(IDX)  = FREQ5(IFREQ,IDX)
            END DO

            IF (LOCDBG) THEN
             WRITE (LUPRI,*) 'CC_5RIND> put on the list:',
     &        LABEL(A),'(',FREQ(A),'),  ', LABEL(B),'(',FREQ(B),'),  ',
     &        LABEL(C),'(',FREQ(C),'),  ', LABEL(D),'(',FREQ(D),'),  ',
     &        LABEL(E),'(',FREQ(E),'),  ', LABEL(F),'(',FREQ(F),')'
            END IF

*           request third-order t response vectors and third-order
*           chi vectors (which implies, that the second-order
*           Lagrangian multipliers will be computed):
            DO IDXA = 1, 6
            DO IDXB = IDXA+1, 6
            DO IDXC = IDXB+1, 6
               INUM = IR3TAMP(LABEL(IDXA),+FREQ(IDXA),ISYM(IDXA),
     &                        LABEL(IDXB),+FREQ(IDXB),ISYM(IDXB),
     &                        LABEL(IDXC),+FREQ(IDXC),ISYM(IDXC) )

               INUM = IR3TAMP(LABEL(IDXA),-FREQ(IDXA),ISYM(IDXA),
     &                        LABEL(IDXB),-FREQ(IDXB),ISYM(IDXB),
     &                        LABEL(IDXC),-FREQ(IDXC),ISYM(IDXC) )

               INUM = ICHI3(LABEL(IDXA),+FREQ(IDXA),ISYM(IDXA),
     &                      LABEL(IDXB),+FREQ(IDXB),ISYM(IDXB),
     &                      LABEL(IDXC),+FREQ(IDXC),ISYM(IDXC) )

               INUM = ICHI3(LABEL(IDXA),-FREQ(IDXA),ISYM(IDXA),
     &                      LABEL(IDXB),-FREQ(IDXB),ISYM(IDXB),
     &                      LABEL(IDXC),-FREQ(IDXC),ISYM(IDXC) )
            END DO
            END DO
            END DO

          END DO
        END IF
      END DO


      RETURN
      END
*---------------------------------------------------------------------*
c /* deck cc_tpaind */
*=====================================================================*
       SUBROUTINE CC_TPAIND 
*---------------------------------------------------------------------*
*
*    Purpose: Determine which vectors are needed for the calculation
*             of two-photon absorption strength
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <cctpainf.h>
#include <ccrspprp.h>
#include <ccexci.h>
#include <ccexcinf.h>
#include <ccroper.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <ccslvinf.h>

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

* variables:
      CHARACTER*8 LABELA, LABELB
      INTEGER ISYMB, ISYMA, ISYMAB, ISYME
      INTEGER INUM, IOPPAIR, IDX, ISTATE, IEXCI, IRSD
      
      REAL*8  FREQA, FREQB, EIGV 

* external functions:
      INTEGER IROPER
      INTEGER ICHI2
      INTEGER IRHSR2
      INTEGER IR1TAMP
      INTEGER ILRMAMP
      INTEGER IL1ZETA

* data:
      LOGICAL FIRSTCALL
      SAVE    FIRSTCALL
      DATA    FIRSTCALL /.TRUE./


      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CC_TPAIND> NSMOPER = ',NSMOPER
      END IF

*---------------------------------------------------------------------*
* test if operators are available and translate IASMOP, IBSMOP, ICSMOP
* and IDSMOP arrays from the PRPLBL list to the new list maintained 
* by IROPER.
*---------------------------------------------------------------------*
      IF (FIRSTCALL) THEN

       IOPPAIR = 1
       DO WHILE (IOPPAIR .LE. NSMOPER)

        LABELA = PRPLBL(IASMOP(IOPPAIR))
        LABELB = PRPLBL(IBSMOP(IOPPAIR))

        IF (      (IROPER(LABELA,ISYMA) .LT. 0)
     &       .OR. (IROPER(LABELB,ISYMB) .LT. 0) ) THEN
          WRITE(LUPRI,'(/2X,5A,/2X,2A)')
     &     ' WARNING: ONE OF THE OPERATORS WITH THE LABELS "',
     &     LABELA,'", "', LABELB,'" IS NOT AVAILABLE.',
     &     ' SECOND MOMENT CROSS SECTION CALCULATION IS CANCELED ',
     &     ' FOR THIS OPERATOR PAIR.'
          DO IDX = IOPPAIR, NSMOPER-1
            IASMOP(IDX) = IASMOP(IDX+1)
            IBSMOP(IDX) = IBSMOP(IDX+1)
          END DO
          NSMOPER = NSMOPER - 1
        ELSE 
          IASMOP(IOPPAIR) = IROPER(LABELA,ISYMA) 
          IBSMOP(IOPPAIR) = IROPER(LABELB,ISYMB)
          IOPPAIR = IOPPAIR + 1
        END IF

       END DO

       FIRSTCALL = .FALSE.

      END IF ! (FIRSTCALL)

*---------------------------------------------------------------------*
* if no states were selected use by default all states:
*---------------------------------------------------------------------*
      IF ( .NOT. SELSMST ) THEN
         NSMSEL = 0
         DO ISYME = 1, NSYM
           DO IEXCI = 1, NCCEXCI(ISYME,1)
             NSMSEL = NSMSEL + 1
             ISMSEL(NSMSEL,1) = ISYME
             ISMSEL(NSMSEL,2) = IEXCI
           END DO
         END DO
      END IF

*---------------------------------------------------------------------*
* if HALFFR flag is given, set here the laser frequency:
*---------------------------------------------------------------------*
      IF ( HALFFR .OR. (.NOT. SELSMST) ) THEN
        DO IRSD = 1, NSMSEL
          ISYME = ISMSEL(IRSD,1)        
          IEXCI = ISMSEL(IRSD,2)        
          EIGV  = EIGVAL(ISYOFE(ISYME) + IEXCI)        
          
          BSMFR(IRSD) = 0.5d0 * EIGV 
        END DO
      END IF

*---------------------------------------------------------------------*
* for CC3 we can switch off LTPA_USE_O2 & LTPA_USE_X2:
*---------------------------------------------------------------------*
      IF (CC3 .AND. LTPA_USE_O2) THEN
        WRITE(LUPRI,*) 'Info: the .USE O2 option cannot be use for '
        WRITE(LUPRI,*) '      in *CCTPA for CC3... it is turned off' 
        LTPA_USE_O2 = .FALSE.
      END IF

      IF (CC3 .AND. LTPA_USE_X2) THEN
        WRITE(LUPRI,*) 'Info: the .USE X2 option cannot be use for '
        WRITE(LUPRI,*) '      in *CTPA for CC3... it is turned off' 
        LTPA_USE_X2 = .FALSE.
      END IF

*---------------------------------------------------------------------*
* set list entries for the required response vectors:
* note that for S^0f_AB,AB(w_B) we need M^AB_0f(-w_B) and M^AB_f0(w_B)
*---------------------------------------------------------------------*
      DO IOPPAIR = 1, NSMOPER
       LABELA = LBLOPR(IASMOP(IOPPAIR))
       LABELB = LBLOPR(IBSMOP(IOPPAIR))

       ISYMA  = ISYOPR(IASMOP(IOPPAIR))
       ISYMB  = ISYOPR(IBSMOP(IOPPAIR))
       ISYMAB = MULD2H(ISYMA,ISYMB)

       DO IRSD = 1, NSMSEL
 
        ISYME = ISMSEL(IRSD,1) ! irrep

        IF (ISYME.EQ.ISYMAB) THEN
         
         IEXCI  = ISMSEL(IRSD,2)        ! state number within irrep
         ISTATE = ISYOFE(ISYME) + IEXCI ! number over all irreps
         EIGV   = EIGVAL(ISTATE)        ! excitation energie

         FREQB  = BSMFR(IRSD)  ! frequency for field B
         FREQA  = EIGV-FREQB   ! frequency for field A

         IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'CC_TPAIND> put on the list:',
     &       LABELA,'(',FREQA,'),  ', LABELB,'(',FREQB,'),  ',
     &      ISTATE,EIGV
         END IF

         IF (LTPA_USE_X2) THEN
*          request second order chi vectors:
           INUM = ICHI2(LABELA,.FALSE.,-FREQA,ISYMA,
     &                  LABELB,.FALSE.,-FREQB,ISYMB)
         END IF
 
         IF (LTPA_USE_O2) THEN
*          request second-order rhs vectors 
           INUM =IRHSR2(LABELA,.FALSE.,-FREQA,ISYMA,
     &                  LABELB,.FALSE.,-FREQB,ISYMB) 
           INUM =IRHSR2(LABELA,.FALSE.,+FREQA,ISYMA,
     &                  LABELB,.FALSE.,+FREQB,ISYMB) 
         END IF

*        request first order t response vectors:
         INUM = IR1TAMP(LABELA,.FALSE.,-FREQA,ISYMA)
         INUM = IR1TAMP(LABELB,.FALSE.,-FREQB,ISYMB)
         INUM = IR1TAMP(LABELB,.FALSE.,+FREQB,ISYMB)
         INUM = IR1TAMP(LABELA,.FALSE.,+FREQA,ISYMA)


*        request unrelaxed first order zeta response vectors:
         IF (CCSLV) THEN
           INUM = IL1ZETA(LABELA,.FALSE.,-FREQA,ISYMA)
           INUM = IL1ZETA(LABELB,.FALSE.,-FREQB,ISYMB)
           INUM = IL1ZETA(LABELA,.FALSE.,FREQA,ISYMA)
           INUM = IL1ZETA(LABELB,.FALSE.,FREQB,ISYMB)
         ELSE
           INUM = IL1ZETA(LABELA,.FALSE.,-FREQA,ISYMA)
           INUM = IL1ZETA(LABELB,.FALSE.,-FREQB,ISYMB)
         END IF

*        request M1 lagrangian multipliers:
         INUM = ILRMAMP(ISTATE,EIGV,ISYME)

        END IF
       END DO
      END DO

      RETURN
      END
*---------------------------------------------------------------------*
c /* deck cc_tmind */
*=====================================================================*
       SUBROUTINE CC_TMIND 
*---------------------------------------------------------------------*
*
*    Purpose: Determine which vectors are needed in third moment
*             calculations, flags are set for the following :
*             chi vectors , second order rhs vectors,
*             first order t respons vectors, m vectors
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <cctm.h>
#include <cctminf.h>
#include <ccrspprp.h>
#include <ccexci.h>
#include <ccroper.h>

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

* variables:
      CHARACTER*8 LABELA, LABELB, LABELC,
     *            LABELD, LABELE, LABELF 
      INTEGER ISYMB, ISYMC, ISYMA, ISYMD, ISYME, ISYMF, ISYMABC
      INTEGER IFREQ, INUM, IOPER, IDX, IOFFST, I
      
      REAL*8  FREQEX, FREQB, FREQC, EIGV 

* external functions:
      INTEGER IROPER
      INTEGER ICHI3
      INTEGER ILRMAMP
      INTEGER IRHSR3

* data:
      LOGICAL FIRSTCALL
      SAVE    FIRSTCALL
      DATA    FIRSTCALL /.TRUE./


      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CC_TMIND> NTMOPER = ',NTMOPER
      END IF

      IF (FIRSTCALL) THEN

*---------------------------------------------------------------------*
* test if operators are available and translate IATMOP, IBTMOP, ICTMOP
* IDTMOP, IETMOP and IFTMOP arrays from the PRPLBL 
* list to the new list maintained by IROPER.
*---------------------------------------------------------------------*
       IOPER = 1
       DO WHILE (IOPER .LE. NTMOPER)

        LABELA = PRPLBL(IATMOP(IOPER))
        LABELB = PRPLBL(IBTMOP(IOPER))
        LABELC = PRPLBL(ICTMOP(IOPER))
        LABELD = PRPLBL(IDTMOP(IOPER))
        LABELE = PRPLBL(IETMOP(IOPER))
        LABELF = PRPLBL(IFTMOP(IOPER))

        IF (      (IROPER(LABELA,ISYMA) .LT. 0)
     &       .OR. (IROPER(LABELB,ISYMB) .LT. 0)
     &       .OR. (IROPER(LABELC,ISYMC) .LT. 0) 
     &       .OR. (IROPER(LABELD,ISYMD) .LT. 0)
     &       .OR. (IROPER(LABELE,ISYME) .LT. 0) 
     &       .OR. (IROPER(LABELF,ISYMF) .LT. 0) ) THEN

          WRITE(LUPRI,'(/2X,A, /2X,7A/2X,4A)')
     &     ' WARNING: ONE OF THE OPERATORS WITH THE LABELS "',
     &     LABELA,'", "', LABELB,'", "', LABELC,'","',
     &     LABELD,'", "', LABELE,'", "', LABELF,
     &     '" IS NOT AVAILABLE.',
     &     ' THIRD MOMENT CROSS SECTION CALCULATION IS CANCELED ', 
     &     ' FOR THIS OPERATOR SIXTUPLE.'

          DO IDX = IOPER, NTMOPER-1
            IATMOP(IDX) = IATMOP(IDX+1)
            IBTMOP(IDX) = IBTMOP(IDX+1)
            ICTMOP(IDX) = ICTMOP(IDX+1)
            IDTMOP(IDX) = IDTMOP(IDX+1)
            IETMOP(IDX) = IETMOP(IDX+1)
            IFTMOP(IDX) = IFTMOP(IDX+1)
          END DO

          NTMOPER = NTMOPER - 1

        ELSE 
          IATMOP(IOPER) = IROPER(LABELA,ISYMA) 
          IBTMOP(IOPER) = IROPER(LABELB,ISYMB)
          ICTMOP(IOPER) = IROPER(LABELC,ISYMC)
          IDTMOP(IOPER) = IROPER(LABELD,ISYMD)
          IETMOP(IOPER) = IROPER(LABELE,ISYME)
          IFTMOP(IOPER) = IROPER(LABELF,ISYMF)

          IOPER = IOPER + 1
        END IF

       END DO

       FIRSTCALL = .FALSE.

      END IF ! (FIRSTCALL)

*--------------------------------------------------------------------*
* sort list of selected states according to symmetry and canonical
* order within each symmetry
*--------------------------------------------------------------------*

      CALL CC_TMSORT

*---------------------------------------------------------------------*
* set list entries for the required response vectors:
*---------------------------------------------------------------------*
      DO IOPER = 1, NTMOPER
        LABELA = LBLOPR(IATMOP(IOPER))
        LABELB = LBLOPR(IBTMOP(IOPER))
        LABELC = LBLOPR(ICTMOP(IOPER))
        LABELD = LBLOPR(IDTMOP(IOPER))
        LABELE = LBLOPR(IETMOP(IOPER))
        LABELF = LBLOPR(IFTMOP(IOPER))

        ISYMA  = ISYOPR(IATMOP(IOPER))
        ISYMB  = ISYOPR(IBTMOP(IOPER))
        ISYMC  = ISYOPR(ICTMOP(IOPER))
        ISYMD  = ISYOPR(IDTMOP(IOPER))
        ISYME  = ISYOPR(IETMOP(IOPER))
        ISYMF  = ISYOPR(IFTMOP(IOPER))

        ISYMABC = MULD2H(MULD2H(ISYMA,ISYMB),ISYMC) 
        IF (ISYMABC. EQ. MULD2H( MULD2H(ISYMD,ISYMF),ISYME) ) THEN
          DO I = 1, NTMSELX(ISYMABC)  
              IFREQ  = ITMSELX(ISYMABC) + I
              FREQEX  = EXTMFR(IFREQ)
              FREQB  = BTMFR(IFREQ)
              FREQC  = CTMFR(IFREQ)
              IF (LOCDBG) THEN
                WRITE (LUPRI,*) 'CC_TMIND> put on the list:',
     &            LABELA,'(',FREQEX,'),  ', LABELB,'(',FREQB,'),  ',
     &            LABELC,'(',FREQC,'),  ',
     &           IFREQ,FREQEX
              END IF


*           request third order chi vectors:
  
           INUM = ICHI3(LABELA,-FREQEX+FREQB+FREQC,ISYMA,
     &                  LABELB,-FREQB,ISYMB,LABELC,-FREQC,ISYMC)
           INUM = ICHI3(LABELD,-FREQEX+FREQB+FREQC,ISYMD,
     &                  LABELE,-FREQB,ISYME,LABELF,-FREQC,ISYMF)

 
*           request third order rhs vectors 
 
           INUM = IRHSR3(LABELA,-FREQEX+FREQB+FREQC,ISYMA,
     &                  LABELB,-FREQB,ISYMB,LABELC,-FREQC,ISYMC)
           INUM = IRHSR3(LABELA,+FREQEX-FREQB-FREQC,ISYMA,
     &                  LABELB,+FREQB,ISYMB,LABELC,+FREQC,ISYMC)
           INUM = IRHSR3(LABELD,-FREQEX+FREQB+FREQC,ISYMD,
     &                  LABELE,-FREQB,ISYME,LABELF,-FREQC,ISYMF)
           INUM = IRHSR3(LABELD,+FREQEX-FREQB-FREQC,ISYMD,
     &                  LABELE,+FREQB,ISYME,LABELF,+FREQC,ISYMF)

*           request m vectors for different excitation energies


            IOFFST = ISYOFE(ISYMABC) +  ITMSEL(IFREQ,2)
            EIGV  =  EIGVAL(IOFFST)
            INUM   = ILRMAMP(IOFFST,EIGV,ISYMABC)
            CALL FLSHFO(LUPRI)
c           WRITE(LUPRI,*) ' ioffst,eigv,inum,isymabc,ifreq'
c           CALL FLSHFO(LUPRI)
c           WRITE (LUPRI,*) ioffst,eigv,inum,isymabc,ifreq

          END DO

        END IF

      END DO


      RETURN
      END
*---------------------------------------------------------------------*
c /* deck cc_mcdind */
*=====================================================================*
       SUBROUTINE CC_MCDIND(WORK,LWORK) 
*---------------------------------------------------------------------*
*  Purpose: Determine which vectors are needed in magnetic circular
*           dichroism calculations 
*  Flags are set for: 2nd-order rhs vectors for T^AB, 
*                     1st-order T^X (w_X) response amplitudes 
*                     M^f(w_f) lagrangian vectors, 
*                     eigenvectors responses E^fX, Ebar^fX
*                     1st order rhs vectors for Tbar^A (eta part)
*                     projected Tbar^A (PL1)
*
*  Written by Sonia Coriani
*  Version: 04/04-2000
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccmcdinf.h>
#include <ccrspprp.h>
#include <ccexcinf.h>
#include <ccexci.h>
#include <ccroper.h>

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

* variables:
      CHARACTER*8 LABELA, LABELB, LABELC, LABSOP
      INTEGER IOPA,IOPB,IOPC
      LOGICAL LORXA,LORXB,LORXC, LPDBSA,LPDBSB,LPDBSC, SKIP_IT, LRELAX
      INTEGER ISYMA, ISYMB, ISYMC, ISYMAB, ISYMS_F, ISYMS, ISTATE
      INTEGER IEIGV_F, ISTAT_F, IEXCI_F, INUM, IOPER, IDX, IDXS
      INTEGER ISGNSOP,ISYSOP,NLORX,LWORK
      LOGICAL LPROJ
      
      REAL*8  EIGVA_F, ZERO, WORK(LWORK) 

      PARAMETER ( ZERO = 0.0d0 )

* external functions:
      INTEGER IROPER
      INTEGER IRHSR1
      INTEGER IRHSR2
      INTEGER IR1TAMP
      INTEGER IL1ZETA
      INTEGER ILRMAMP
      INTEGER IER1AMP
      INTEGER IEL1AMP
      INTEGER IETA1
      INTEGER IPL1ZETA

* data:
      LOGICAL FIRSTCALL
      SAVE    FIRSTCALL
      DATA    FIRSTCALL /.TRUE./

*----------------------------------------------------------------------*
* Begin
*----------------------------------------------------------------------*

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CC_MCDIND> NMCDOPER = ',NMCDOPER
      END IF

*----------------------------------------------------------------------*
* test if operators are available and translate IAMCDOP,IBMCDOP,ICMCDOP
* arrays from the PRPLBL list to the new list maintained by IROPER.
* Note that NMCDOPER is the number of operator-triples (r,L,r)
*----------------------------------------------------------------------*
  
      LPROJ = .FALSE.

      IF (FIRSTCALL) THEN

        IOPER = 1
        DO WHILE (IOPER .LE. NMCDOPER)

          SKIP_IT = .FALSE.
          LABELA = PRPLBL(IAMCDOP(IOPER))
          LABELB = PRPLBL(IBMCDOP(IOPER))
          LABELC = PRPLBL(ICMCDOP(IOPER))
          LORXA  = LAMCDRX(IOPER)
          LORXB  = LBMCDRX(IOPER)
          LORXC  = LCMCDRX(IOPER)
          IOPA   = IROPER(LABELA,ISYMA)
          IOPB   = IROPER(LABELB,ISYMB)
          IOPC   = IROPER(LABELC,ISYMC)

          WRITE(LUPRI,'(/2X,A,3(1X,A),A)')
     &     'CHECK TRIPLET:',LABELA, LABELB, LABELC
          CALL FLSHFO(LUPRI)


          IF ( (IOPA.LT.0) .OR. (IOPB.LT.0) .OR. (IOPC.LT.0) ) THEN

             WRITE(LUPRI,'(/2X,7A,/2X,2A)')
     &     ' WARNING: ONE OF THE OPERATORS WITH THE LABELS "',
     &       LABELA,'", "', LABELB,'", "', LABELC,'" IS NOT AVAILABLE.',
     &     ' MAGNE.CIRCUL.DICHR. CALCULATION IS CANCELED FOR THIS',
     &     ' OPERATOR TRIPLET.'
 
             SKIP_IT = .TRUE.
          END IF

          NLORX = 0
          IF (LORXA .OR. LPDBSOP(IOPA)) NLORX = NLORX + 1
          IF (LORXB .OR. LPDBSOP(IOPB)) NLORX = NLORX + 1
          IF (LORXC .OR. LPDBSOP(IOPC)) NLORX = NLORX + 1

          IF (NLORX.GT.1) THEN
            WRITE(LUPRI,'(/2X,8A,/2X,A,/2X,A)')
     &       ' WARNING: OPERATOR TRIPLET "',
     &         LABELA,'", "', LABELB,'", "', LABELC,'"',
     &       ' WITH MORE THAN ONE FIELD WHICH',
     &       ' INVOKES ORBITAL RELAXATION OR A PERTUR.-DEP. BASIS SET.',
     &       ' CALCULATION IS CANCELED FOR THIS OPERATOR TRIPLE.'
          END IF
          IF (.NOT. SKIP_IT) THEN
             ! if we have field-dependent basis sets, we need also
             ! to check, if the second-derivative integrals for this
             ! perturbation pair are available
             IF (LPDBSOP(IOPA) .OR. LPDBSOP(IOPB)) THEN
                CALL CC_FIND_SO_OP(LABELA,LABELB,LABSOP,ISYSOP,
     &                             ISGNSOP,INUM,WORK,LWORK)
                IF (INUM.LT.0) SKIP_IT = .TRUE.
             END IF
             IF (LPDBSOP(IOPA) .OR. LPDBSOP(IOPC)) THEN
                CALL CC_FIND_SO_OP(LABELA,LABELC,LABSOP,ISYSOP,
     &                             ISGNSOP,INUM,WORK,LWORK)
                IF (INUM.LT.0) SKIP_IT = .TRUE.
             END IF
             IF (LPDBSOP(IOPB) .OR. LPDBSOP(IOPC)) THEN
                CALL CC_FIND_SO_OP(LABELB,LABELC,LABSOP,ISYSOP,
     &                             ISGNSOP,INUM,WORK,LWORK)
                IF (INUM.LT.0) SKIP_IT = .TRUE.
             END IF
             IF (SKIP_IT) THEN
               WRITE(LUPRI,'(/2X,7A,/2X,A,/2X,A)')
     &          ' WARNING: FOR THE OPERATOR TRIPLET "',
     &            LABELA,'", "', LABELB,'", "', LABELC,'"',
     &         ' A SEC. ORD. OPERATOR IS MISSING.',
     &         ' CALCULATION IS IGNORED.'
             END IF
          END IF



          IF (SKIP_IT) THEN
            DO IDX = IOPER, NMCDOPER-1
              IAMCDOP(IDX) = IAMCDOP(IDX+1)
              IBMCDOP(IDX) = IBMCDOP(IDX+1)
              ICMCDOP(IDX) = ICMCDOP(IDX+1)
              LAMCDRX(IDX) = LAMCDRX(IDX+1)
              LBMCDRX(IDX) = LBMCDRX(IDX+1)
              LCMCDRX(IDX) = LCMCDRX(IDX+1)
            END DO
            NMCDOPER = NMCDOPER - 1        !decrease # of triplets
          ELSE 
            WRITE(LUPRI,'(/2X,A,3(1X,A),A)')
     &      'PUT TRIPLET:',LABELA, LABELB, LABELC,' ONTO THE LIST.'

            IAMCDOP(IOPER) = IROPER(LABELA,ISYMA)
            IBMCDOP(IOPER) = IROPER(LABELB,ISYMB)
            ICMCDOP(IOPER) = IROPER(LABELC,ISYMC)
            IOPER = IOPER + 1
          END IF

        END DO

        FIRSTCALL = .FALSE.

      END IF                             ! end if (FIRSTCALL)
*
*--------------------------------------------------------------------*
* set/check now symmetries and indices of the excited states
*--------------------------------------------------------------------*
*
      IF (SELMCDST) THEN

* check if all required states available, if not remove them from the list:

        IDXS = 1
        DO WHILE (IDXS .LE. NMCDST)
          IF ( IMCDSTNR(IDXS).GT.NCCEXCI(IMCDSTSY(IDXS),1)) THEN
            WRITE(LUPRI,'(/2X,A,I2,A,I2,A,/2X,A)')
     &       ' WARNING: THE STATE WITH SYMMETRY ',IMCDSTSY(IDXS),
     &       ' AND INDEX ',IMCDSTNR(IDXS) ,
     &       ' IS NOT AVAILABLE.',
     &       ' B TERM CALCULATION IS CANCELED FOR THIS STATE.'
            DO IDX = IDXS, NMCDST-1
              IMCDSTNR(IDX) = IMCDSTNR(IDX+1)     !move next index
              IMCDSTSY(IDX) = IMCDSTSY(IDX+1)     !one step back
            END DO
            NMCDST = NMCDST - 1
          ELSE
            IDXS = IDXS + 1
          END IF
        END DO
      
      ELSE

* Use default: MCD for all states specified in *CCEXCI

        DO ISYMS = 1, NSYM
          DO ISTATE = 1, NCCEXCI(ISYMS,1)
            IF (NMCDST.LT.MXMCDST) THEN
              NMCDST = NMCDST + 1
              IMCDSTSY(NMCDST) = ISYMS
              IMCDSTNR(NMCDST) = ISTATE
            END IF
          END DO
        END DO


      END IF

*
*--------------------------------------------------------------------*
* set list entries for all the required response vectors:
* NMCDOPER is # of operator triples (A,B,C)
*--------------------------------------------------------------------*
*
      DO 100 IOPER = 1, NMCDOPER
 
        LPROJ = .FALSE.

        LABELA = LBLOPR(IAMCDOP(IOPER))             !get labels back
        LABELB = LBLOPR(IBMCDOP(IOPER))
        LABELC = LBLOPR(ICMCDOP(IOPER))
        
        LPDBSA = LPDBSOP(IAMCDOP(IOPER))
        LPDBSB = LPDBSOP(IBMCDOP(IOPER))
        LPDBSC = LPDBSOP(ICMCDOP(IOPER))

        LORXA  = LAMCDRX(IOPER)
        LORXB  = LBMCDRX(IOPER)
        LORXC  = LCMCDRX(IOPER)

        ISYMA  = ISYOPR(IAMCDOP(IOPER))             !get symmetries back
        ISYMB  = ISYOPR(IBMCDOP(IOPER))
        ISYMC  = ISYOPR(ICMCDOP(IOPER))
        
        ISYMAB = MULD2H(ISYMA,ISYMB)

        LRELAX = LORXA.OR.LORXB.OR.LORXC.OR.LPDBSA.OR.LPDBSB.OR.LPDBSC

        WRITE(LUPRI,'(/2X,A,3(1X,A),A)')
     &   'require responses for op. triplet:',LABELA, LABELB, LABELC
          call flshfo(6)

        IF (ISYMAB.EQ.ISYMC) THEN

           DO 101 IDX = 1, NMCDST
              ISYMS_F = IMCDSTSY(IDX)        !symmetry of excited state
              ISTAT_F = IMCDSTNR(IDX)        !index of exc.state within symmetry
              !absolute index of the exc. state (pointer)
              IEXCI_F = ISYOFE(ISYMS_F) + ISTAT_F
              EIGVA_F = EIGVAL(IEXCI_F)              !excitation energy

              IF (ISYMS_F.EQ.ISYMC) THEN

                 IF (LOCDBG) THEN
                    WRITE (LUPRI,*) 'CC_MCDIND> put onto the list:',
     &               LABELA,'(',-EIGVA_F,'),  ', LABELB,'(',ZERO,'),  ',
     &               IEXCI_F,EIGVA_F
                 END IF

                 INUM = IR1TAMP(LABELA,LORXA,-EIGVA_F,ISYMA)
                 INUM = IR1TAMP(LABELB,LORXB,ZERO,ISYMB)
                 INUM = ILRMAMP(IEXCI_F,EIGVA_F,ISYMC)
                 INUM = IER1AMP(IEXCI_F,EIGVA_F,ISYMC,
     &                          LABELA,-EIGVA_F,ISYMA,.FALSE.)
                 INUM = IETA1(LABELB,LORXB,ZERO,ISYMB)
                 IF (ISYMB .EQ. 1) LPROJ = .TRUE.
                 INUM = IEL1AMP(IEXCI_F,EIGVA_F,ISYMC,
     &                          LABELB, ZERO,ISYMB,LORXB,LPROJ)
                 IF (.NOT.LUSE2N1) THEN
                   INUM = IR1TAMP(LABELC,LORXC,-EIGVA_F,ISYMC)      
                 END IF
                 INUM = IRHSR1(LABELC,LORXC,EIGVA_F,ISYMC)
                 INUM = IETA1(LABELC,LORXC,EIGVA_F,ISYMC)

                 IF (LUSEPL1) THEN
                    IF (ISYMB .EQ. 1) LPROJ = .TRUE.
                    INUM = IPL1ZETA(LABELA,LORXA,-EIGVA_F,ISYMA,
     &                              LPROJ,IEXCI_F, EIGVA_F,ISYMC)
                 ELSE

                    INUM = IRHSR2(LABELA,LORXA,-EIGVA_F,ISYMA,
     &                            LABELB,LORXB,ZERO,ISYMB) 
                    IF (ISYMB .EQ. 1) LPROJ = .TRUE.
                    INUM = IER1AMP(IEXCI_F,EIGVA_F,ISYMC,
     &                             LABELB,ZERO,ISYMB,LPROJ)
                    INUM = IETA1(LABELA,LORXA,-EIGVA_F,ISYMA)
                 END IF

              END IF
 101       CONTINUE
        END IF
 100  CONTINUE                  
      CALL FLSHFO(LUPRI)

      RETURN
      END
*---------------------------------------------------------------------*
c /* deck cc_exlrind */
*=====================================================================*
       SUBROUTINE CC_EXLRIND
*---------------------------------------------------------------------*
*
*    Purpose: setup of the equations that have to be solved for
*             the excited state linear response properties
*
*    Written by Christof Haettig, July 1997.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccexlrinf.h>
#include <ccrspprp.h>
#include <ccroper.h>
#include <cclr.h>
#include <ccexci.h>
#include <ccsdinp.h>

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

* variables:
      CHARACTER*8 LABELA, LABELB
      LOGICAL LPRJ
      INTEGER ISYMB, ISYMA, IFREQ, INUM, IOPER, ISYMS, ISTATE
      INTEGER IDX, IEXCII, ISTATI, IEXCIF, ISTATF, IDXS, ISYMSI, ISYMSF
      
      REAL*8  HALF, FREQA, FREQB, EIGVI, EIGVF

      PARAMETER ( HALF = 0.5d0 )


* external functions:
      INTEGER IER1AMP
      INTEGER IEL1AMP
      INTEGER IROPER
      INTEGER IRHSR2
      INTEGER IN2AMP

* data:
      LOGICAL FIRSTCALL
      SAVE    FIRSTCALL
      DATA    FIRSTCALL /.TRUE./

*---------------------------------------------------------------------*
* test if operators are available and translate IAQROP, IBQROP, ICQROP
* arrays from the PRPLBL list to the new list maintained by IROPER.
*---------------------------------------------------------------------*
      IF (FIRSTCALL) THEN
       WRITE (LUPRI,*) 'CC_EXLRIND> NEXLROPER = ',NEXLROPER

       IOPER = 1
       DO WHILE (IOPER .LE. NEXLROPER)

          WRITE(LUPRI,'(/2X,A,3I5)') 
     & 'IOPER,IAEXLROP,IBEXLROP:',IOPER,IAEXLROP(IOPER),IBEXLROP(IOPER)
        LABELA = PRPLBL(IAEXLROP(IOPER))
        LABELB = PRPLBL(IBEXLROP(IOPER))
          WRITE(LUPRI,'(/2X,A,2(1X,A),A)') 'CHECK PAIR:',LABELA, LABELB

        IF (      (IROPER(LABELA,ISYMA) .LT. 0)
     &       .OR. (IROPER(LABELB,ISYMB) .LT. 0) ) THEN

          WRITE(LUPRI,'(/2X,5A,/2X,2A)')
     &     ' WARNING: ONE OF THE OPERATORS WITH THE LABELS "',
     &     LABELA,'", "', LABELB,'" IS NOT AVAILABLE.',
     &     ' POLARIZABILITY CALCULATION IS CANCELED FOR THIS',
     &     ' OPERATOR PAIR.'

          DO IDX = IOPER, NEXLROPER-1
            IAEXLROP(IDX) = IAEXLROP(IDX+1)
            IBEXLROP(IDX) = IBEXLROP(IDX+1)
          END DO

          NEXLROPER = NEXLROPER - 1

        ELSE 
          WRITE(LUPRI,'(/2X,A,2(1X,A),A)')
     &     'PUT PAIR:',LABELA, LABELB,' ONT THE LIST.'
          IAEXLROP(IOPER) = IROPER(LABELA,ISYMA) 
          IBEXLROP(IOPER) = IROPER(LABELB,ISYMB)

          IOPER = IOPER + 1
        END IF

       END DO

       FIRSTCALL = .FALSE.

      END IF ! (FIRSTCALL)

*---------------------------------------------------------------------*
* process the excited state information
*---------------------------------------------------------------------*
      IF (ALLSTATES) THEN

* set now symmetries and indeces of the excited states:
* (diagonal cases, i.e., excited state response functions, only)
        DO ISYMS = 1, NSYM
          DO ISTATE = 1, NCCEXCI(ISYMS,1)
            IF (NEXLRST.LT.MXEXLRST) THEN
              NEXLRST = NEXLRST + 1
              IELRSYM(NEXLRST,1) = ISYMS
              IELRSTA(NEXLRST,1) = ISTATE
              IELRSYM(NEXLRST,2) = ISYMS
              IELRSTA(NEXLRST,2) = ISTATE
            END IF
          END DO
        END DO

      ELSE 

* check if all states available, if not remove them from the list:
        IDXS = 1
        DO WHILE (IDXS .LE. NEXLRST)
          IF ( IELRSTA(IDXS,1).GT.NCCEXCI(IELRSYM(IDXS,1),1)
     &        .OR. IELRSTA(IDXS,2).GT.NCCEXCI(IELRSYM(IDXS,2),1) ) THEN
            WRITE(LUPRI,'(2(/2X,A,I2,A,I2),A,/2X,A)')
     &       ' WARNING: THE STATE WITH SYMMETRY ',IELRSYM(IDXS,1),
     &       ' AND INDEX ',IELRSTA(IDXS,1) ,
     &       ' OR THE STATE WITH SYMMETRY ',IELRSYM(IDXS,2),
     &       ' AND INDEX ',IELRSTA(IDXS,2) ,
     &       ' IS NOT AVAILABLE.',
     &       ' POLARIZABILITY CALCULATION IS CANCELED FOR THIS STATE.'
            DO IDX = IDXS, NEXLRST-1
              IELRSTA(IDX,1) = IELRSTA(IDX+1,1)
              IELRSYM(IDX,1) = IELRSYM(IDX+1,1)
              IELRSTA(IDX,2) = IELRSTA(IDX+1,2)
              IELRSYM(IDX,2) = IELRSYM(IDX+1,2)
            END DO
            NEXLRST = NEXLRST - 1
          ELSE
            IDXS = IDXS + 1
          END IF
        END DO

      END IF

*---------------------------------------------------------------------*
* check for HALFFR option:
*---------------------------------------------------------------------*
      IF ( HALFFR .AND. NEXLRFREQ.NE.1 ) THEN
        WRITE (LUPRI,*) 'error in CC_EXLRIND: HALFFR option is',
     &             ' incompatible with a frequency list.' 
        CALL QUIT('error in CC_EXLRIND.')
      END IF

*---------------------------------------------------------------------*
* for CC3 we can switch off USE_O2/USE_EL1 since it can not be used: 
*---------------------------------------------------------------------*
      IF (CC3 .AND. USE_O2) THEN
        WRITE(LUPRI,*) 'Info: the .USE O2 option cannot be use for '
        WRITE(LUPRI,*) '      in *CCEXLR for CC3... it is turned off' 
        USE_O2 = .FALSE.
      END IF

      IF (CC3 .AND. USE_EL1) THEN
        WRITE(LUPRI,*) 'Info: the .USELEF option cannot be use for '
        WRITE(LUPRI,*) '      in *CCEXLR for CC3... it is turned off' 
        USE_EL1 = .FALSE.
      END IF

*---------------------------------------------------------------------*
* set list entries for the required response vectors:
*---------------------------------------------------------------------*
      DO IOPER = 1, NEXLROPER
        LABELA = LBLOPR(IAEXLROP(IOPER))
        LABELB = LBLOPR(IBEXLROP(IOPER))

        ISYMA  = ISYOPR(IAEXLROP(IOPER))
        ISYMB  = ISYOPR(IBEXLROP(IOPER))

C         WRITE(LUPRI,'(/2X,A,3(1X,A),A)')
C    &     'require responses for pair:',LABELA, LABELB

        

      DO IDXS = 1, NEXLRST
        ISYMSI = IELRSYM(IDXS,1) 
        ISTATI = IELRSTA(IDXS,1) 
        ISYMSF = IELRSYM(IDXS,2) 
        ISTATF = IELRSTA(IDXS,2) 
        IEXCII = ISYOFE(ISYMSI) + ISTATI
        EIGVI  = EIGVAL(IEXCII)
        IEXCIF = ISYOFE(ISYMSF) + ISTATF
        EIGVF  = EIGVAL(IEXCIF)

        IF (MULD2H(ISYMA,ISYMB) .EQ. MULD2H(ISYMSI,ISYMSF) ) THEN

          DO IFREQ = 1, NEXLRFREQ
            FREQB  = BEXLRFR(IFREQ)
            IF (IEXCII.EQ.IEXCIF) THEN
              FREQA  = -FREQB
              LPRJ   = .NOT. NOPROJ
            ELSE
              IF ( HALFFR )  FREQB = HALF * (EIGVI-EIGVF)
              FREQA  = EIGVI - EIGVF -FREQB
              LPRJ   = .FALSE.
            END IF

*           request first order right excited state response vectors:
            IF (.NOT. USE_EL1) THEN
             INUM=IER1AMP(IEXCIF,EIGVF,ISYMSF,LABELA,+FREQA,ISYMA,LPRJ)
             INUM=IER1AMP(IEXCIF,EIGVF,ISYMSF,LABELB,+FREQB,ISYMB,LPRJ)
             INUM=IER1AMP(IEXCII,EIGVI,ISYMSI,LABELA,-FREQA,ISYMA,LPRJ)
             INUM=IER1AMP(IEXCII,EIGVI,ISYMSI,LABELB,-FREQB,ISYMB,LPRJ)
            END IF

*           request first order left excited state response vectors:
            IF (USE_EL1) THEN
             INUM=IEL1AMP(IEXCII,EIGVI,ISYMSI,
     &                    LABELA,+FREQA,ISYMA,.FALSE.,LPRJ)
             INUM=IEL1AMP(IEXCII,EIGVI,ISYMSI,
     &                    LABELB,+FREQB,ISYMB,.FALSE.,LPRJ)
             INUM=IEL1AMP(IEXCIF,EIGVF,ISYMSF,
     &                    LABELA,-FREQA,ISYMA,.FALSE.,LPRJ)
             INUM=IEL1AMP(IEXCIF,EIGVF,ISYMSF,
     &                    LABELB,-FREQB,ISYMB,.FALSE.,LPRJ)
            END IF

*           request zeroth-order excited state lagrange vectors:
            INUM = IN2AMP(IEXCII,-EIGVI,ISYMSI,IEXCIF,EIGVF,ISYMSF)
            INUM = IN2AMP(IEXCIF,-EIGVF,ISYMSF,IEXCII,EIGVI,ISYMSI)

*           request right hand side vector for T2:
            IF (USE_O2) THEN
              INUM = IRHSR2(LABELA,.FALSE.,+FREQA,ISYMA,
     &                      LABELB,.FALSE.,+FREQB,ISYMB)
              INUM = IRHSR2(LABELA,.FALSE.,-FREQA,ISYMA,
     &                      LABELB,.FALSE.,-FREQB,ISYMB)
            END IF

          END DO
        END IF
      END DO

      END DO


      RETURN
      END
*---------------------------------------------------------------------*
*=====================================================================*
C  /* Deck iroper */
*=====================================================================*
      INTEGER FUNCTION IROPER(NEWLBL,ISYM)
*---------------------------------------------------------------------*
*
* maintain the list of operators labels for the response calculations
* the operators are specified by a character*8 label (NEWLBL) 
*
* in difference to the list maintained by the INDPRPCC function,
* the list maintained by IROPER is ordered (see routine CCLSTSORT).
*
* Christof Haettig, November 1996, modified Januar 97: 
*
*   if NEWLBL is on the list return list index and set ISYM,
*   if NEWLBL is NOT on the list:
*        LOPROPN=.true.  --> extend list, and return index
*        LOPROPN=.false. --> return -1   
*
*=====================================================================*
      IMPLICIT NONE
#include <ccroper.h>
#include <priunit.h>
C
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*8 NEWLBL
      INTEGER I, ISYM

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'IROPER>',NEWLBL,ISYM
        CALL FLSHFO(LUPRI)
      END IF

      DO I = 1,NRSOLBL
         IF ( NEWLBL .EQ. LBLOPR(I) ) THEN
            IROPER = I
            ISYM   = ISYOPR(IROPER)
            IF (LOCDBG)
     &        WRITE(LUPRI,*) 'IROPER>',IROPER,LBLOPR(IROPER),
     &           ISYOPR(IROPER)
            RETURN
         END IF
      END DO  

      IF (LOPROPN) THEN
        NRSOLBL = NRSOLBL + 1

        IF (NRSOLBL.GT.MAXOLBL) THEN
         WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *    ' NUMBER OF SPECIFIED OPERATORS EXCEED THE MAXIMUM ALLOWED',
     *    ' MAXOLBL =',MAXOLBL,' NRSOLBL= ',NRSOLBL
         CALL QUIT(' IROPER: TOO MANY OPERATORS SPECIFIED')
        END IF

        LBLOPR(NRSOLBL) = NEWLBL
        ISYOPR(NRSOLBL) = ISYM
        IROPER = NRSOLBL

      ELSE
        WRITE(LUPRI,'(/3A)')
     *   ' WARNING: OPERATOR WITH LABEL "',NEWLBL,'" NOT AVAILABLE.'
        IROPER = -1
      END IF

      IF (LOCDBG)
     &  WRITE (LUPRI,*) 
     &      'IROPER>', IROPER, LBLOPR(IROPER), ISYOPR(IROPER)

      RETURN
      END
*=====================================================================*
C  /* Deck ir2tamp */
      INTEGER FUNCTION IR2TAMP(NEWLBLA,LORXA,FRQANEW,ISYMA,
     *                         NEWLBLB,LORXB,FRQBNEW,ISYMB )
*---------------------------------------------------------------------*
C
C maintain the list of second order right response vectors
C
C   if vector is on the list return list index and set ISYMA,ISYMB
C   if vector is NOT on the list:
C        LR2OPN=.true.  --> extend list, and return index
C        LR2OPN=.false. --> return -1   
C
C        NEWLBLA / NEWLBLB -- operator labels
C        LORXA   / LORXB   -- flags for orbital relaxation
C        FRQANEW / FRQBNEW -- frequencies
C        ISYMA   / ISYMB   -- symmetries    
C
C Christof Haettig, Februar 97
C LORXA, LORXB flags introduced in July 1999
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <ccr2rsp.h>
#include <priunit.h>
C
      LOGICAL LORXA, LORXB
      INTEGER ISYMA, ISYMB
      REAL*8  FRQANEW,FRQBNEW,TOL

      PARAMETER(TOL=1.0D-12)

      CHARACTER*8 NEWLBLA, NEWLBLB
      INTEGER I

      DO I = 1,NR2TLBL
         IF ( NEWLBLA.EQ.LBLAR2T(I).AND. NEWLBLB.EQ.LBLBR2T(I)
     *       .AND. (LORXA .EQV. LORXAR2T(I))
     *       .AND. (LORXB .EQV. LORXBR2T(I))  
     *       .AND. (ABS(FRQANEW-FRQAR2T(I)).LT.TOL)
     *       .AND. (ABS(FRQBNEW-FRQBR2T(I)).LT.TOL) 
     *      ) THEN
            IR2TAMP = I
            ISYMA   = ISYAR2T(IR2TAMP)
            ISYMB   = ISYBR2T(IR2TAMP)
            RETURN
         END IF
         IF ( NEWLBLB.EQ.LBLAR2T(I).AND. NEWLBLA.EQ.LBLBR2T(I)
     *       .AND. (LORXB .EQV. LORXAR2T(I))
     *       .AND. (LORXA .EQV. LORXBR2T(I))  
     *       .AND. (ABS(FRQBNEW-FRQAR2T(I)).LT.TOL)
     *       .AND. (ABS(FRQANEW-FRQBR2T(I)).LT.TOL) 
     *      ) THEN
            IR2TAMP = I
            ISYMB   = ISYAR2T(IR2TAMP)
            ISYMA   = ISYBR2T(IR2TAMP)
            RETURN
         END IF
      END DO  

      IF (LR2OPN) THEN
        NR2TLBL = NR2TLBL + 1

        IF (NR2TLBL.GT.MAXT2LBL) THEN
          WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *    '@ NUMBER OF SPECIFIED EQUATIONS EXCEED THE MAXIMUM ALLOWED',
     *    '@ MAXT2LBL =',MAXT2LBL,' NR2TLBL= ',NR2TLBL
          CALL QUIT(' IR2TAMP: TOO MANY EQUATIONS SPECIFIED')
        END IF

        LBLAR2T(NR2TLBL)  = NEWLBLA
        LBLBR2T(NR2TLBL)  = NEWLBLB
        LORXAR2T(NR2TLBL) = LORXA
        LORXBR2T(NR2TLBL) = LORXB
        FRQAR2T(NR2TLBL)  = FRQANEW
        FRQBR2T(NR2TLBL)  = FRQBNEW
        ISYAR2T(NR2TLBL)  = ISYMA
        ISYBR2T(NR2TLBL)  = ISYMB
        IR2TAMP = NR2TLBL

      ELSE
        WRITE(LUPRI,'(3A,L2,A,1P,D12.5,3A,L2,A,1P,D12.5,2A)') 
     *   '@ WARNING: R2 VECTOR FOR ',
     *            NEWLBLA,'(',LORXA,',',FRQANEW,'), ',
     *            NEWLBLB,'(',LORXB,',',FRQBNEW,')',
     *              ' IS NOT AVAILABLE.'
        IR2TAMP = -1
      END IF

      RETURN
      END
*=====================================================================*
C  /* Deck ir1tamp */
      INTEGER FUNCTION IR1TAMP(NEWLBL,LORX,FRQNEW,ISYM)
*---------------------------------------------------------------------*
C
C maintain the list of first order t amplitude responses
C
C   if vector is on the list return list index and set ISYM
C   if vector is NOT on the list:
C        LR1OPN=.true.  --> extend list, and return index
C        LR1OPN=.false. --> return -1
C
C        NEWLBL -- operator label
C        LORX   -- flag for orbital relaxation
C        FRQNEW -- frequency (ignored for unrelaxed orbitals)
C        ISYM   -- symmetry
C
C Christof Haettig, Oktober 1996
C LORX flag introduced and some clean up in Juni 1998
C
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <ccr1rsp.h>
#include <priunit.h>
 
      LOGICAL LORX, LORXI
      INTEGER ISYM
      REAL*8  FRQNEW,TOL

      PARAMETER(TOL=1.0D-12)

      CHARACTER*8 NEWLBL
      INTEGER I

      DO I = 1,NLRTLBL
         IF ( NEWLBL .EQ. LRTLBL(I) .AND. (LORX .EQV. LORXLRT(I)) .AND.
     *       (ABS(FRQNEW-FRQLRT(I)) .LT. TOL) ) THEN
            IR1TAMP = I
            ISYM    = ISYLRT(IR1TAMP)
            RETURN
         END IF
      END DO  

      IF (LR1OPN) THEN
        NLRTLBL = NLRTLBL + 1

        IF (NLRTLBL.GT.MAXTLBL) THEN
          WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *    '@ NUMBER OF SPECIFIED EQUATIONS EXCEED THE MAXIMUM ALLOWED',
     *    '@ MAXTLBL =',MAXTLBL,' NLRTLBL= ',NLRTLBL
          CALL QUIT(' IR1TAMP: TOO MANY EQUATIONS SPECIFIED')
        END IF

        LRTLBL(NLRTLBL)  = NEWLBL
        LORXLRT(NLRTLBL) = LORX
        FRQLRT(NLRTLBL)  = FRQNEW
        ISYLRT(NLRTLBL)  = ISYM
        IR1TAMP          = NLRTLBL

      ELSE
        WRITE(LUPRI,'(/3A,L2,A,1P,D12.5,2A)')
     *   '@ WARNING: R1 VECTOR FOR ',NEWLBL,'(',LORX,',',FRQNEW,')',
     *                  ' IS NOT AVAILABLE.'
        WRITE(LUPRI,'(/A)') ' LIST OF FIRST-ORDER T VECTORS:'
        DO I = 1, NLRTLBL
           WRITE(LUPRI,'(I5,3X,A8,L3,I5,2X,1P,D15.6)') 
     &            I, LRTLBL(I), LORXLRT(I), ISYLRT(I), FRQLRT(I)
           WRITE (LUPRI,*) 
     &            ' NEWLBL .EQ. LRTLBL(I):', NEWLBL .EQ. LRTLBL(I) 
           WRITE (LUPRI,*) 
     &            '(LORX .EQV. LORXLRT(I)):',(LORX .EQV. LORXLRT(I)) 
           WRITE (LUPRI,*) 
     &            'FRQNEW=FRQLRT:',(ABS(FRQNEW-FRQLRT(I)) .LT. TOL) 
        END DO
        IR1TAMP = -1
      END IF

      RETURN
      END
*=====================================================================*
C  /* Deck ir1kappa */
      INTEGER FUNCTION IR1KAPPA(NEWLBL,FRQNEW,ISYM)
*---------------------------------------------------------------------*
C
C maintain the list of first order orbital responses
C
C   if vector is on the list return list index and set ISYM
C   if vector is NOT on the list:
C        LR1OPN=.true.  --> extend list, and return index
C        LR1OPN=.false. --> return -1
C
C        NEWLBL -- operator label
C        FRQNEW -- frequency (ignored for unrelaxed orbitals)
C        ISYM   -- symmetry
C
C   Note that this list shares common block with IR1TAMP list
C
C Christof Haettig, July 2003
C
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <ccr1rsp.h>
#include <priunit.h>
 
      INTEGER ISYM
      REAL*8  FRQNEW,TOL

      PARAMETER(TOL=1.0D-12)

      CHARACTER*8 NEWLBL
      INTEGER I

      DO I = 1,NLRTHFLBL
         IF ( NEWLBL .EQ.  LRTHFLBL(I) .AND. 
     *       (ABS(FRQNEW-FRQLRTHF(I)) .LT. TOL) ) THEN
            IR1KAPPA = I
            ISYM     = ISYLRTHF(IR1KAPPA)
            RETURN
         END IF
      END DO  

      IF (LR1OPN) THEN
        NLRTHFLBL = NLRTHFLBL + 1

        IF (NLRTHFLBL.GT.MAXTLBL) THEN
          WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *    '@ NUMBER OF SPECIFIED EQUATIONS EXCEED THE MAXIMUM ALLOWED',
     *    '@ MAXTLBL =',MAXTLBL,' NLRTHFLBL= ',NLRTHFLBL
          CALL QUIT(' IR1KAPPA: TOO MANY EQUATIONS SPECIFIED')
        END IF

        LRTHFLBL(NLRTHFLBL) = NEWLBL
        FRQLRTHF(NLRTHFLBL) = FRQNEW
        ISYLRTHF(NLRTHFLBL) = ISYM
        IR1KAPPA            = NLRTHFLBL

      ELSE
        WRITE(LUPRI,'(/3A,1P,D12.5,2A)')
     *   '@ WARNING: R1 KAPPA VECTOR FOR ',NEWLBL,'(',FRQNEW,')',
     *                  ' IS NOT AVAILABLE.'
        WRITE(LUPRI,'(/A)') ' LIST OF FIRST-ORDER KAPPA VECTORS:'
        DO I = 1, NLRTHFLBL
           WRITE(LUPRI,'(I5,3X,A8,I5,2X,1P,D15.6)') 
     &            I, LRTHFLBL(I), ISYLRTHF(I), FRQLRTHF(I)
        END DO
        IR1KAPPA = -1
      END IF

      RETURN
      END
*=====================================================================*
*=====================================================================*
C /* Deck ilrcamp */
      INTEGER FUNCTION ILRCAMP(NEWLBL,ICAUCH,ISYM)
*---------------------------------------------------------------------*
C
C maintain the list of right cauchy vectors. 
C
C   if vector is on the list return list index and set symmetry
C   if vector is NOT on the list, then
C      if LRC1OPN = .true.  --> extend list and return index
C      if LRC1OPN = .false. --> return -1
C
C Christof Haettig, october 1997 
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <ccrc1rsp.h>
#include <priunit.h>
C
      INTEGER ISYM,ICAUCH

      CHARACTER*8 NEWLBL
      INTEGER I

      DO I = 1,NLRCLBL
         IF ( NEWLBL .EQ. LRCLBL(I).AND.
     *      (ICAUCH.EQ.ILRCAU(I))) THEN
            ILRCAMP = I
            ISYM    = ISYLRC(ILRCAMP)
            RETURN
         END IF
      END DO  

      IF (LRC1OPN) THEN
        NLRCLBL = NLRCLBL + 1

        IF (NLRCLBL.GT.MAXCLBL) THEN
          WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *    '@ NUMBER OF SPECIFIED EQUATIONS EXCEED THE MAXIMUM ALLOWED',
     *    '@ MAXCLBL =',MAXCLBL,' NLRCLBL= ',NLRCLBL
          CALL QUIT(' ILRCAMP: TOO MANY EQUATIONS SPECIFIED')
        END IF

        LRCLBL(NLRCLBL) = NEWLBL
        ILRCAU(NLRCLBL) = ICAUCH
        ISYLRC(NLRCLBL) = ISYM
        ILRCAMP = NLRCLBL
      ELSE
        WRITE(LUPRI,'(3A,I3,A)') 
     *   '@ WARNING: RC1 VECTOR FOR ',NEWLBL,
     *   ' CAUCHY ORDER',ICAUCH,' IS NOT AVAILABLE.'
        ILRCAMP = -1
      END IF

      RETURN
      END
*---------------------------------------------------------------------*
*=====================================================================*
C /* Deck ILC1AMP */
      INTEGER FUNCTION ILC1AMP(NEWLBL,ICAUCH,ISYM)
*---------------------------------------------------------------------*
C
C maintain the list of left cauchy vectors. 
C
C   if vector is on the list return list index and set symmetry
C   if vector is NOT on the list, then
C      if LLC1OPN = .true.  --> extend list and return index
C      if LLC1OPN = .false. --> return -1
C
C Christof Haettig, october 1997 
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <cclc1rsp.h>
#include <priunit.h>
C
      INTEGER ISYM,ICAUCH

      CHARACTER*8 NEWLBL
      INTEGER I

      DO I = 1,NLC1LBL
         IF ( NEWLBL.EQ.LBLLC1(I) .AND. ICAUCH.EQ.ILC1CAU(I) ) THEN
            ILC1AMP = I
            ISYM    = ISYLC1(ILC1AMP)
            RETURN
         END IF
      END DO  

      IF (LLC1OPN) THEN
        NLC1LBL = NLC1LBL + 1

        IF (NLC1LBL.GT.MAXLC1LBL) THEN
          WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *    '@ NUMBER OF SPECIFIED EQUATIONS EXCEED THE MAXIMUM ALLOWED',
     *    '@ MAXLC1LBL =',MAXLC1LBL,' NLC1LBL= ',NLC1LBL
          CALL QUIT(' ILC1AMP: TOO MANY EQUATIONS SPECIFIED')
        END IF

        LBLLC1(NLC1LBL)  = NEWLBL
        ILC1CAU(NLC1LBL) = ICAUCH
        ISYLC1(NLC1LBL)  = ISYM
        ILC1AMP          = NLC1LBL
      ELSE
        WRITE(LUPRI,'(3A,I3,A)') 
     *   '@ WARNING: LC1 VECTOR FOR ',NEWLBL, 
     *   ' CAUCHY ORDER',ICAUCH,' IS NOT AVAILABLE.'
        ILC1AMP = -1
      END IF

      RETURN
      END
*---------------------------------------------------------------------*
*=====================================================================*
C  /* Deck il1zeta */
      INTEGER FUNCTION IL1ZETA(NEWLBL,LORX,FRQNEW,ISYM)
*---------------------------------------------------------------------*
C
C maintain the list of first order zeta amplitude responses
C
C   if vector is on the list return list index and set ISYMA,ISYMB
C   if vector is NOT on the list:
C        LL1OPN=.true.  --> extend list, and return index
C        LL1OPN=.false. --> return -1
C
C        NEWLBL -- operator label
C        LORX   -- flag for orbital relaxation
C        FRQNEW -- frequency (ignored for unrelaxed orbitals)
C        ISYM   -- symmetry
C
C Christof Haettig, Oktober 1996
C LORX flag introduced and some clean up in Juni 1998
C
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <ccl1rsp.h>
#include <priunit.h>
C

      LOGICAL LORX
      INTEGER ISYM, I
      REAL*8  FRQNEW, TOL

      PARAMETER(TOL=1.0D-12)


      CHARACTER*8 NEWLBL

      DO I = 1,NLRZLBL
         IF ( NEWLBL .EQ. LRZLBL(I).AND. (LORX .EQV. LORXLRZ(I)) .AND.
     *      (ABS(FRQNEW-FRQLRZ(I)).LT.TOL)) THEN
            IL1ZETA = I
            ISYM    = ISYLRZ(IL1ZETA)
            RETURN
         END IF
      END DO  

      IF (LL1OPN) THEN
        NLRZLBL = NLRZLBL + 1

        IF (NLRZLBL.GT.MAXZLBL) THEN
          WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *    '@ NUMBER OF SPECIFIED EQUATIONS EXCEED THE MAXIMUM ALLOWED',
     *    '@ MAXZLBL =',MAXZLBL,' NLRZLBL= ',NLRZLBL
          CALL QUIT(' IL1ZETA: TOO MANY EQUATIONS SPECIFIED')
        END IF

        LRZLBL(NLRZLBL)  = NEWLBL
        LORXLRZ(NLRZLBL) = LORX
        FRQLRZ(NLRZLBL)  = FRQNEW
        ISYLRZ(NLRZLBL)  = ISYM
        IL1ZETA          = NLRZLBL

      ELSE
        WRITE(LUPRI,'(3A,L2,A,1P,D12.5,2A)')
     *   '@ WARNING: L1 VECTOR FOR ',NEWLBL,'(',LORX,',',FRQNEW,')',
     *                  ' IS NOT AVAILABLE.'
        IL1ZETA = -1
      END IF

      RETURN
      END
*=====================================================================*
*---------------------------------------------------------------------*
      INTEGER FUNCTION ILRMAMP(IEXCI,FRQNEW,ISYM)
C
C maintain the list of transition moment lagrangian multipliers
C
C Ove Christiansen April 1997
C
      IMPLICIT NONE
#include <cclrmrsp.h>
#include <priunit.h>
C
      INTEGER ISYM,IEXCI,I
      REAL*8  FRQNEW,TOL

      PARAMETER(TOL=1.0D-12)

      DO I = 1,NLRM
         IF ( IEXCI .EQ. ILRM(I).AND.
     *      (ABS(FRQNEW-FRQLRM(I)).LT.TOL)) THEN
            ILRMAMP = I
            ISYM    = ISYLRM(ILRMAMP)
            RETURN
         END IF
      END DO

      NLRM    = NLRM    + 1

      IF (NLRM   .GT.MAXM   ) THEN
         WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *   '@ NUMBER OF SPECIFIED EQUATIONS EXCEED THE MAXIMUM ALLOWED',
     *   '@ MAXM    =',MAXM   ,' NLRM   = ',NLRM
         CALL QUIT(' ILRMAMP: TOO MANY EQUATIONS SPECIFIED')
      END IF

      ILRM(NLRM)      = IEXCI
      FRQLRM(NLRM)    = FRQNEW
      ISYLRM(NLRM)    = ISYM
      ILRMAMP = NLRM

      RETURN
      END
*---------------------------------------------------------------------*
C  /* Deck cclstsort */
*=====================================================================*
      SUBROUTINE CCLSTSORT (TYPE, 
     &                      ISYMS, ISTAT, EIGVAL,
     &                      ISYMO, LABEL, FREQ, ICAU, LORX,
     &                      ISYOF, NVEC,  MAXVEC, LPROJ )
*---------------------------------------------------------------------*
*
* PURPOSE: sort list of vectors/equations according:
*
*          1.)  over-all symmetry (obtained by calling ILSTSYM)
*          2.)  individual symmetries of the states  (ISYMS)
*          3.)  state indeces  (ISTAT) 
*          4.)  projection flag (LPROJ)
*          5.)  over-all cauchy order (ICAU)
*          6.)  individual symmetries of the first operators (ISYMO)
*          7.)  operator labels (LABEL)
*          8.)  frequencies (FREQ)
*          9.)  individual cauchy orders (ICAU)
*         10.)  orbital relaxation (LORX)
*
*          sets up symmetry offsets ISYOF
*
*          print sorted lists to output
*
*          number of operators (sym., labels, freqs, istat) used
*          depends on TYPE (see subroutines CCLSTCMP and CCLSTSWAP).
*          EIGVAL array is not used for comparison, but is sorted 
*          with the list.
*
*          implemented: o1, 
*                       O1, O2, O3, 
*                       R1, R2, R3, 
*                       X1, X2, X3, 
*                       L1, L2, L3, 
*                           CO2
*                       RC, CR2
*                           CX2
*                       LC, CL2
*                       M1
*                       N2
*                       ER1, ER2
*                       EL1, EL2
*                       PL1
*
*          not tested for RE, LE, E0
*
*
* Christof Haettig, October 1996
* generalized for open ended strategy may 1997
* projection flag 1998
* orbital relaxation flag 1999
* PL1 vectors, Sonia march 2000
*
Cholesky
* Swapped sorting order so that LABEL sort is done
* after FREQ sort by modifying CCLSTCMP
* tbp 2003. Only tested for linear response!
Cholesky
*
*=====================================================================*
      IMPLICIT NONE
#include <priunit.h>
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER NVEC, MAXVEC, NSTAT, ORDER, JSYM, I
      LOGICAL LPROJ(MAXVEC), LORX(MAXVEC,*)
      CHARACTER*(8) LABEL(MAXVEC,*)
      INTEGER ISYMO(MAXVEC,*), ISTAT(MAXVEC,*), ISYOF(8)
      INTEGER ISYMS(MAXVEC,*), ICAU(MAXVEC,*)
      CHARACTER*(*) TYPE

      REAL*8  FREQ(MAXVEC,*)
      REAL*8  EIGVAL(MAXVEC,*)
      
      LOGICAL CHANGES
      INTEGER IVEC

* external functions:
      LOGICAL CCLSTCMP 
      INTEGER ILSTSYM

* check TYPE and determine number of states involved and resp. order:
      IF      (     TYPE(1:2).EQ.'R1' .OR. TYPE(1:2).EQ.'L1' 
     &         .OR. TYPE(1:3).EQ.'O1 '.OR. TYPE(1:3).EQ.'X1 '
     &         .OR. TYPE(1:2).EQ.'RC' .OR. TYPE(1:2).EQ.'LC'
     &         .OR. TYPE(1:2).EQ.'o1'                        ) THEN
        NSTAT = 0
        ORDER = 1
      ELSE IF (     TYPE(1:2).EQ.'R2' .OR. TYPE(1:2).EQ.'L2'
     &         .OR. TYPE(1:2).EQ.'O2' .OR. TYPE(1:2).EQ.'X2' 
     &         .OR. TYPE(1:3).EQ.'CR2'.OR. TYPE(1:3).EQ.'CL2'
     &         .OR. TYPE(1:3).EQ.'CO2'.OR. TYPE(1:3).EQ.'CX2' ) THEN
        NSTAT = 0
        ORDER = 2
      ELSE IF (     TYPE(1:2).EQ.'R3' .OR. TYPE(1:2).EQ.'L3'
     &         .OR. TYPE(1:2).EQ.'O3' .OR. TYPE(1:2).EQ.'X3' ) THEN
        NSTAT = 0
        ORDER = 3
      ELSE IF (     TYPE(1:2).EQ.'R4' .OR. TYPE(1:2).EQ.'L4'
     &         .OR. TYPE(1:2).EQ.'O4' .OR. TYPE(1:2).EQ.'X4' ) THEN
        NSTAT = 0
        ORDER = 4
      ELSE IF (     TYPE(1:2).EQ.'RE' .OR. TYPE(1:2).EQ.'LE' 
     &         .OR. TYPE(1:2).EQ.'E0' .OR. TYPE(1:2).EQ.'M1' ) THEN
        NSTAT = 1
        ORDER = 0
      ELSE IF (     TYPE(1:2).EQ.'N2'                        ) THEN
        NSTAT = 2
        ORDER = 0
      ELSE IF (     TYPE(1:3).EQ.'ER1'.OR. TYPE(1:3).EQ.'EL1') THEN
        NSTAT = 1
        ORDER = 1
      ELSE IF (     TYPE(1:3).EQ.'ER2'.OR. TYPE(1:3).EQ.'EL2') THEN
        NSTAT = 1
        ORDER = 2
      ELSE IF (     TYPE(1:3).EQ.'PL1') THEN
        NSTAT = 1
        ORDER = 1
      ELSE
        WRITE (LUPRI,*) 'unknown list ',TYPE,' in CCLSTSORT.'
        CALL QUIT('unknown list TYPE in CCLSTSORT.')
      END IF

* bubble sort:
      CHANGES = .TRUE.

      DO WHILE (CHANGES)
        CHANGES = .FALSE.

        DO IVEC = 1, NVEC-1
          IF( CCLSTCMP(TYPE,NSTAT,ORDER,IVEC,ISYMS,ISTAT,
     &         ISYMO,LABEL,FREQ,ICAU,LORX,MAXVEC,LPROJ) ) THEN

            CALL CCLSTSWAP(TYPE,NSTAT,ORDER,IVEC, ISYMS,ISTAT,EIGVAL,
     &                     ISYMO,LABEL,FREQ,ICAU,LORX,MAXVEC,LPROJ)

            CHANGES = .TRUE.

          END IF
        END DO

        IF (LOCDBG .AND. (TYPE(2:2).NE.'C'.AND.TYPE(1:1).NE.'C')) THEN
          DO IVEC = 1, NVEC
            WRITE(LUPRI,'(I5,I3,2(3X,A8,I3,2X,1P,D15.6))') 
     &        IVEC, ILSTSYM(TYPE,IVEC),
     &         (ISYMS(IVEC,I),ISTAT(IVEC,I),EIGVAL(IVEC,I),I=1,NSTAT),
     &         (LABEL(IVEC,I),ISYMO(IVEC,I),FREQ(IVEC,I),I=1,ORDER)
          END DO
        ELSE IF ( LOCDBG .AND. (TYPE(2:2).EQ.'C'.OR.TYPE(1:1).EQ.'C')
     &          ) THEN
          WRITE(LUPRI,'(3A)') 'sorted ',TYPE,' list:'
          DO IVEC = 1, NVEC
            JSYM = ILSTSYM(TYPE,IVEC)
            WRITE(LUPRI,'(I5,I3,2(3X,A8,I3,2X,I3))') 
     &        IVEC, JSYM,
     &         (LABEL(IVEC,I),ISYMO(IVEC,I),ICAU(IVEC,I),I=1,ORDER)
          END DO
          CALL FLSHFO(LUPRI)
        END IF

      END DO

      IVEC = 0
      DO JSYM = 1, 8
        ISYOF(JSYM) = IVEC
        IF (NVEC.GT.0) THEN
          DO WHILE(IVEC.LT.NVEC .AND.
     &             ILSTSYM(TYPE,MIN(IVEC+1,NVEC)).EQ.JSYM)
            IVEC = IVEC + 1
          END DO
        END IF
      END DO

      RETURN
      END
*=====================================================================*
*                     END OF SUBROUTINE CCLSTSORT                     *
*=====================================================================*
C  /* Deck cclstcmp */
*=====================================================================*
      LOGICAL FUNCTION CCLSTCMP(TYPE,NSTAT,ORDER,IVEC,ISYMS,ISTAT,
     &                        ISYMO,LABEL,FREQ,ICAU,LORX,MAXVEC,LPROJ) 
*---------------------------------------------------------------------*
* PURPOSE: do comparison for CCLSTSORT according to:
*   
*          1.)  over-all symmetry
*          2.)  individual symmetries of the states  (ISYMS)
*          3.)  state indeces  (ISTAT)
*          4.)  projection flag (LPROJ)
*          5.)  over-all cauchy order
*          6.)  individual symmetries of the first operators (ISYMO)
*          7.)  operator labels (LABEL)
*          8.)  frequencies (FREQ)
*          9.)  individual cauchy orders (ICAU)
*         10.)  orbital relaxation flags (LORX)
*
*          cauchy orders only used for
*                'RC', 'LC', 'CRn', 'COn', 'CLn', 'CXn'
*
*          orbital relaxation flags LORX only used for 
*                'o1 ', 'O1 ', 'R1 ', 'X1 ', 'L1 ' ,'PL1 ', 'EL1 '  
*
*          special treatments: 
*             'o1 '          -- no frequency and no orbital relaxation
*
*
* Christof Haettig, October 1996, 
* generalized for an open ended strategy in may 1997
* PL1 vectors, LORX in EL1 ... Sonia Coriani  2000
*=====================================================================*
      IMPLICIT NONE
#include <priunit.h>
Cholesky
#include <maxorb.h>
#include <ccdeco.h>
Cholesky
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER NSTAT, ORDER, MAXVEC,IVEC, I, IDX, NCAU, NCAU1, NSTAT1
      LOGICAL LPROJ(MAXVEC), LORX(MAXVEC,ORDER)
      CHARACTER*(8) LABEL(MAXVEC,ORDER)
      INTEGER ISYMS(MAXVEC,NSTAT), ISTAT(MAXVEC,NSTAT)
      INTEGER ISYMO(MAXVEC,ORDER), ICAU(MAXVEC,ORDER)
      CHARACTER*(*) TYPE

      REAL*8  FREQ(MAXVEC,ORDER)

* external function:
      INTEGER ILSTSYM
      
*---------------------------------------------------------------------*
* compare over-all symmetry:
*---------------------------------------------------------------------*
      IF      ( ILSTSYM(TYPE,IVEC) .GT. ILSTSYM(TYPE,IVEC+1) ) THEN
         IF (LOCDBG) WRITE (LUPRI,*) 'swap because of overall symmetry.'
         CCLSTCMP = .TRUE.
         RETURN
      ELSE IF ( ILSTSYM(TYPE,IVEC) .LT. ILSTSYM(TYPE,IVEC+1) ) THEN
         CCLSTCMP = .FALSE.
         RETURN
      END IF

*---------------------------------------------------------------------*
* compare the symmetries of the individual states involved:
*---------------------------------------------------------------------*
* we have already sorted according to the over-all symmetry, so for
* zeroth-order vectors we can only sort after NSTAT-1 state symmetries
*
      NSTAT1 = NSTAT
      IF (ORDER.EQ.0) NSTAT1 = NSTAT - 1

      DO IDX = 1, NSTAT1
        IF      ( ISYMS(IVEC,IDX) .GT. ISYMS(IVEC+1,IDX) ) THEN
          IF (LOCDBG) WRITE (LUPRI,*)
     &          'swap because of state symmetries.'
          CCLSTCMP = .TRUE.
          RETURN
        ELSE IF ( ISYMS(IVEC,IDX) .LT. ISYMS(IVEC+1,IDX) ) THEN
          CCLSTCMP = .FALSE.
          RETURN
        END IF
      END DO

*---------------------------------------------------------------------*
* compare the indices of the individual states involved:
*---------------------------------------------------------------------*
      DO IDX = 1, NSTAT
        IF      ( ISTAT(IVEC,IDX) .GT. ISTAT(IVEC+1,IDX) ) THEN
          IF (LOCDBG) WRITE (LUPRI,*) 'swap because of state indices.'
          CCLSTCMP = .TRUE.
          RETURN
        ELSE IF ( ISTAT(IVEC,IDX) .LT. ISTAT(IVEC+1,IDX) ) THEN
          CCLSTCMP = .FALSE.
          RETURN
        END IF
      END DO

*---------------------------------------------------------------------*
* for excited state response vectors or projected response multipliers
* compare projection flag:
*---------------------------------------------------------------------*
      IF (TYPE(1:2).EQ.'ER' .OR. TYPE(1:2).EQ.'EL'
     &                      .OR. TYPE(1:3).EQ.'PL1') THEN
        IF ( (.NOT.LPROJ(IVEC)) .AND. LPROJ(IVEC+1) ) THEN
         IF (LOCDBG) WRITE (LUPRI,*) 'swap because of projection flag.'
         CCLSTCMP = .TRUE.
         RETURN
        ELSE IF ( LPROJ(IVEC) .AND. (.NOT.LPROJ(IVEC+1)) ) THEN
         CCLSTCMP = .FALSE.
         RETURN
        END IF
      END IF
*---------------------------------------------------------------------*
* for cauchy vectors compare over-all cauchy order:
*---------------------------------------------------------------------*
      IF (TYPE(1:2).EQ.'RC' .OR. TYPE(1:2).EQ.'LC' .OR.
     &    TYPE(1:2).EQ.'CR' .OR. TYPE(1:2).EQ.'CL' .OR.
     &    TYPE(1:2).EQ.'CO' .OR. TYPE(1:2).EQ.'CX'      ) THEN

        NCAU  = 0
        NCAU1 = 0
        DO IDX = 1, ORDER
          NCAU  = NCAU  + ICAU(IVEC,IDX)
          NCAU1 = NCAU1 + ICAU(IVEC+1,IDX)
        END DO

        IF      ( NCAU .GT. NCAU1 ) THEN
         IF (LOCDBG) WRITE (LUPRI,*)
     &          'swap because of overall cauchy order.'
         CCLSTCMP = .TRUE.
         RETURN
        ELSE IF ( NCAU .LT. NCAU1 ) THEN
         CCLSTCMP = .FALSE.
         RETURN
        END IF

      END IF

*---------------------------------------------------------------------*
* compare the symmetries of the ORDER-1 first operators
*---------------------------------------------------------------------*
      DO IDX = 1, ORDER-1
        IF      ( ISYMO(IVEC,IDX) .GT. ISYMO(IVEC+1,IDX) ) THEN
          IF (LOCDBG) WRITE (LUPRI,*) 
     &          'swap because of operator symmetries.'
          CCLSTCMP = .TRUE.
          RETURN
        ELSE IF ( ISYMO(IVEC,IDX) .LT. ISYMO(IVEC+1,IDX) ) THEN
          CCLSTCMP = .FALSE.
          RETURN
        END IF
      END DO

* If Cholesky, sort before after frequencies

      IF (.NOT. CHOINT) THEN

*---------------------------------------------------------------------*
* compare the labels
*---------------------------------------------------------------------*
      DO IDX = 1, ORDER 
        DO I = 1, 8
          IF ( LGT(LABEL(IVEC,IDX)(I:I),LABEL(IVEC+1,IDX)(I:I)) ) THEN
            IF (LOCDBG) WRITE (LUPRI,*)
     &            'swap because of operator labels.'
            CCLSTCMP = .TRUE.
            RETURN
          END IF
          IF ( LLT(LABEL(IVEC,IDX)(I:I),LABEL(IVEC+1,IDX)(I:I)) ) THEN
            CCLSTCMP = .FALSE.
            RETURN
          END IF
        END DO
      END DO

      END IF      ! Cholesky

*---------------------------------------------------------------------*
* compare the frequencies
*---------------------------------------------------------------------*
      IF ( TYPE(1:2).NE.'RC' .AND. TYPE(1:2).NE.'LC' .AND.
     &     TYPE(1:2).NE.'CR' .AND. TYPE(1:2).NE.'CL' .AND.
     &     TYPE(1:2).NE.'CO' .AND. TYPE(1:2).NE.'CX' .AND.
     &     TYPE(1:2).NE.'o1'                                ) THEN

        DO IDX = 1, ORDER
          IF      ( FREQ(IVEC,IDX) .GT. FREQ(IVEC+1,IDX) ) THEN
            IF (LOCDBG) WRITE (LUPRI,*) 'swap because of frequencies.'
            CCLSTCMP = .TRUE.
            RETURN
          ELSE IF ( FREQ(IVEC,IDX) .LT. FREQ(IVEC+1,IDX) ) THEN
            CCLSTCMP = .FALSE.
            RETURN
          END IF
        END DO

      END IF

* If Cholesky, sort now after frequencies

      IF (CHOINT) THEN

*---------------------------------------------------------------------*
* compare the labels
*---------------------------------------------------------------------*
      DO IDX = 1, ORDER
        DO I = 1, 8
          IF ( LGT(LABEL(IVEC,IDX)(I:I),LABEL(IVEC+1,IDX)(I:I)) ) THEN
            IF (LOCDBG) WRITE (LUPRI,*)
     &            'swap because of operator labels.'
            CCLSTCMP = .TRUE.
            RETURN
          END IF
          IF ( LLT(LABEL(IVEC,IDX)(I:I),LABEL(IVEC+1,IDX)(I:I)) ) THEN
            CCLSTCMP = .FALSE.
            RETURN
          END IF
        END DO
      END DO

      END IF     ! Cholesky

*---------------------------------------------------------------------*
* compare the cauchy orders:
*---------------------------------------------------------------------*
      IF (TYPE(1:2).EQ.'RC' .OR. TYPE(1:2).EQ.'LC' .OR.
     &    TYPE(1:2).EQ.'CR' .OR. TYPE(1:2).EQ.'CL' .OR.
     &    TYPE(1:2).EQ.'CO' .OR. TYPE(1:2).EQ.'CX'      ) THEN

        DO IDX = 1, ORDER
          IF      ( ICAU(IVEC,IDX) .GT. ICAU(IVEC+1,IDX) ) THEN
            IF (LOCDBG) WRITE (LUPRI,*) 'swap because of cauchy orders.'
            CCLSTCMP = .TRUE.
            RETURN
          ELSE IF ( ICAU(IVEC,IDX) .LT. ICAU(IVEC+1,IDX) ) THEN
            CCLSTCMP = .FALSE.
            RETURN
          END IF
        END DO

      END IF

*---------------------------------------------------------------------*
* compare orbital relaxation flags: 
*---------------------------------------------------------------------*
      IF (TYPE(1:3).EQ.'O1 '.OR. TYPE(1:2).EQ.'R1' .OR.
     &    TYPE(1:3).EQ.'X1 '.OR. TYPE(1:2).EQ.'L1' .OR.
     &    TYPE(1:3).EQ.'PL1'.OR. TYPE(1:3).EQ.'EL1'     ) THEN

        DO IDX = 1, ORDER
          IF ( (.NOT.LORX(IVEC,IDX)) .AND. LORX(IVEC+1,IDX) ) THEN
           IF (LOCDBG) WRITE (LUPRI,*)
     &            'swap because of orb. relax. flag.'
           CCLSTCMP = .TRUE.
           RETURN
          ELSE IF (LORX(IVEC,IDX) .AND. (.NOT.LORX(IVEC+1,IDX))) THEN
           CCLSTCMP = .FALSE.
           RETURN
          END IF
        END DO

      END IF

*---------------------------------------------------------------------*
* both entries are the same???
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(1X,4A)') 'WARNING FROM CCLSTCMP: ',
     &   'The ',TYPE(1:2),' list contains a redundant entry.'
      WRITE (LUPRI,'(1X,A,I2,A,I2,A)') 'Entries ',IVEC,' AND ',IVEC+1,
     &   ' are the same.'

      CCLSTCMP = .FALSE.

      RETURN
      END
*=====================================================================*
*                   END OF SUBROUTINE CCLSTCMP                        *
*=====================================================================*
C  /* Deck cclstswap */
*=====================================================================*
      SUBROUTINE CCLSTSWAP(TYPE,NSTAT,ORDER,IVEC, ISYMS,ISTAT,EIGVAL,
     &                     ISYMO,LABEL,FREQ,ICAU,LORX,MXVEC,LPROJ)
*---------------------------------------------------------------------*
*
* PURPOSE: swap two list elements for CCLSTSORT:
*
*          swaps in general ORDER operators symmetries, labels and 
*          frequencies or cauchy orders, and NSTAT state symmetries,
*          state indeces and eigenvalues
*
*          cauchy orders only used for
*                'RC', 'LC', 'CRn', 'COn', 'CLn', 'CXn'
*
*          orbital relaxation flags LORX only used for 
*                'o1 ', 'O1 ', 'R1 ', 'X1 ', 'L1 ', 'PL1 ', 'EL1 '   
*
*          for 'ELn' and 'ERn' and 'PL1' also the projection 
*                                        flag is swapped
*
*          special treatment:
*              o1 -- no frequency and no orbital relaxation, but 
*                    we swap in addition: ISYMAT, IATOPR, LPDBSOP
*                    
*
* Christof Haettig, October 1996
* generalized for an open ended strategy in may 1997
* Sonia Coriani: PL1 and LORX for EL1
*=====================================================================*
      IMPLICIT NONE
#include <priunit.h>
#include <ccroper.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER IVEC, MXVEC, ORDER, IDX, NSTAT, ISWAP
      LOGICAL LPROJ(MXVEC), LORX(MXVEC,ORDER), LSWAP
      CHARACTER*(8) LABEL(MXVEC,ORDER), LBLSWAP
      INTEGER ISYMO(MXVEC,ORDER)
      INTEGER ICAU(MXVEC,ORDER)
      INTEGER ISYMS(MXVEC,NSTAT)
      INTEGER ISTAT(MXVEC,NSTAT)
      CHARACTER*(*) TYPE

      REAL*8  FREQ(MXVEC,ORDER), EIGVAL(MXVEC,NSTAT), RSWAP

*---------------------------------------------------------------------*
* swap symmetries:
*---------------------------------------------------------------------*
      DO IDX = 1, ORDER
        ISWAP             = ISYMO(IVEC,IDX)
        ISYMO(IVEC,IDX)   = ISYMO(IVEC+1,IDX)
        ISYMO(IVEC+1,IDX) = ISWAP
      END DO

      DO IDX = 1, NSTAT
        ISWAP             = ISYMS(IVEC,IDX)
        ISYMS(IVEC,IDX)   = ISYMS(IVEC+1,IDX)
        ISYMS(IVEC+1,IDX) = ISWAP
      END DO

*---------------------------------------------------------------------*
* swap state indices and eigenvalues:
*---------------------------------------------------------------------*
      DO IDX = 1, NSTAT
        ISWAP             = ISTAT(IVEC,IDX)
        ISTAT(IVEC,IDX)   = ISTAT(IVEC+1,IDX)
        ISTAT(IVEC+1,IDX) = ISWAP
      END DO

      DO IDX = 1, NSTAT
        RSWAP              = EIGVAL(IVEC,IDX)
        EIGVAL(IVEC,IDX)   = EIGVAL(IVEC+1,IDX)
        EIGVAL(IVEC+1,IDX) = RSWAP 
      END DO

*---------------------------------------------------------------------*
* swap projection flag:
*---------------------------------------------------------------------*
      IF (TYPE(1:2).EQ.'ER' .OR. TYPE(1:2).EQ.'EL'
     &                      .OR. TYPE(1:3).EQ.'PL1') THEN
        LSWAP         = LPROJ(IVEC)
        LPROJ(IVEC)   = LPROJ(IVEC+1)
        LPROJ(IVEC+1) = LSWAP
      END IF

*---------------------------------------------------------------------*
* swap labels:
*---------------------------------------------------------------------*
      DO IDX = 1, ORDER
        LBLSWAP           = LABEL(IVEC,IDX)
        LABEL(IVEC,IDX)   = LABEL(IVEC+1,IDX)
        LABEL(IVEC+1,IDX) = LBLSWAP
      END DO

*---------------------------------------------------------------------*
* swap frequencies:
*---------------------------------------------------------------------*
      IF ( TYPE(1:2).NE.'RC'  .AND. TYPE(1:2).NE.'LC' .AND.
     &     TYPE(1:2).NE.'CR'  .AND. TYPE(1:2).NE.'CO' .AND.
     &     TYPE(1:2).NE.'CL'  .AND. TYPE(1:2).NE.'CX' .AND. 
     &     TYPE(1:2).NE.'o1'                                ) THEN
        DO IDX = 1, ORDER
          RSWAP            = FREQ(IVEC,IDX)
          FREQ(IVEC,IDX)   = FREQ(IVEC+1,IDX)
          FREQ(IVEC+1,IDX) = RSWAP
        END DO
      END IF

*---------------------------------------------------------------------*
* swap cauchy orders:
*---------------------------------------------------------------------*
      IF (TYPE(1:2).EQ.'RC' .OR. TYPE(1:2).EQ.'LC' .OR.
     &    TYPE(1:2).EQ.'CR' .OR. TYPE(1:2).EQ.'CL' .OR.
     &    TYPE(1:2).EQ.'CO' .OR. TYPE(1:2).EQ.'CX'      ) THEN
        DO IDX = 1, ORDER
          ISWAP            = ICAU(IVEC,IDX)
          ICAU(IVEC,IDX)   = ICAU(IVEC+1,IDX)
          ICAU(IVEC+1,IDX) = ISWAP
        END DO
      END IF

*---------------------------------------------------------------------*
* swap orbital relaxation flags: 
*---------------------------------------------------------------------*
      IF (TYPE(1:3).EQ.'O1 '.OR. TYPE(1:2).EQ.'R1' .OR.
     &    TYPE(1:3).EQ.'X1 '.OR. TYPE(1:2).EQ.'L1' .OR.
     &    TYPE(1:3).EQ.'PL1'.OR. TYPE(1:3).EQ.'EL1'   ) THEN
        DO IDX = 1, ORDER
          LSWAP            = LORX(IVEC,IDX)
          LORX(IVEC,IDX)   = LORX(IVEC+1,IDX)
          LORX(IVEC+1,IDX) = LSWAP
        END DO
      END IF

*---------------------------------------------------------------------*
* for 'o1' list swap in addition: ISYMAT, IATOPR, LPDBSOP
*---------------------------------------------------------------------*
      IF (TYPE(1:2).EQ.'o1') THEN
        ISWAP           = ISYMAT(IVEC)
        ISYMAT(IVEC)    = ISYMAT(IVEC+1)
        ISYMAT(IVEC+1)  = ISWAP

        ISWAP           = IATOPR(IVEC)
        IATOPR(IVEC)    = IATOPR(IVEC+1)
        IATOPR(IVEC+1)  = ISWAP

        LSWAP           = LPDBSOP(IVEC)
        LPDBSOP(IVEC)   = LPDBSOP(IVEC+1)
        LPDBSOP(IVEC+1) = LSWAP
      END IF

*---------------------------------------------------------------------*
* return:
*---------------------------------------------------------------------*
      RETURN

      END
*=====================================================================*
*                   END OF SUBROUTINE CCLSTSWAP                       *
*=====================================================================*
*=====================================================================*
C  /* Deck ilstsym */
*=====================================================================*
      INTEGER FUNCTION ILSTSYM(LIST, INDEX)
*---------------------------------------------------------------------*
* PURPOSE: get symmetry for vector on list
*
*          LIST : list type
*          INDEX: index of the vector on the list
*
* Christof Haettig, November 1996
* PL1 introduced Sonia
*=====================================================================*
      IMPLICIT NONE
#include <priunit.h>
#include <ccorb.h>
#include <ccroper.h>
#include <cclrmrsp.h>
#include <ccer1rsp.h>
#include <ccer2rsp.h>
#include <ccel1rsp.h>
#include <ccel2rsp.h>
#include <ccr1rsp.h>
#include <ccr2rsp.h>
#include <ccr3rsp.h>
#include <ccr4rsp.h>
#include <ccx1rsp.h>
#include <ccx2rsp.h>
#include <ccx3rsp.h>
#include <ccx4rsp.h>
#include <cco1rsp.h>
#include <cco2rsp.h>
#include <cco3rsp.h>
#include <cco4rsp.h>
#include <ccl1rsp.h>
#include <ccl2rsp.h>
#include <ccl3rsp.h>
#include <ccl4rsp.h>
#include <ccn2rsp.h>
#include <ccrc1rsp.h>
#include <cclc1rsp.h>
#include <cccr2rsp.h>
#include <ccco2rsp.h>
#include <cccx2rsp.h>
#include <cccl2rsp.h>
#include <ccexci.h>
#include <ccpl1rsp.h>

      CHARACTER*(3) LIST
      INTEGER INDEX
      LOGICAL LEOOR

      LEOOR = .FALSE.

* begin:
      IF (LIST(1:2).EQ.'o1') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NRSOLBL) LEOOR = .TRUE.
         ILSTSYM = ISYOPR(INDEX)
      ELSE IF (LIST(1:3).EQ.'O1 '.OR.LIST(1:3).EQ.'O1e') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NO1LBL) LEOOR = .TRUE.
         ILSTSYM = ISYO1(INDEX)
      ELSE IF (LIST(1:2).EQ.'O2' ) THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NO2LBL) LEOOR = .TRUE.
         ILSTSYM = MULD2H(ISYAO2(INDEX),ISYBO2(INDEX))
      ELSE IF (LIST(1:2).EQ.'O3' ) THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NO3LBL) LEOOR = .TRUE.
         ILSTSYM = MULD2H(ISYO3(INDEX,1),ISYO3(INDEX,2))
         ILSTSYM = MULD2H(ILSTSYM,ISYO3(INDEX,3))
      ELSE IF (LIST(1:2).EQ.'O4' ) THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NO4LBL) LEOOR = .TRUE.
         ILSTSYM = MULD2H(ISYO4(INDEX,1),ISYO4(INDEX,2))
         ILSTSYM = MULD2H(ILSTSYM,ISYO4(INDEX,3))
         ILSTSYM = MULD2H(ILSTSYM,ISYO4(INDEX,4))
      ELSE IF (LIST(1:3).EQ.'CO2') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NCO2LBL) LEOOR = .TRUE.
         ILSTSYM = MULD2H(ISYCO2(INDEX,1),ISYCO2(INDEX,2))
      ELSE IF (LIST(1:3).EQ.'X1 '.OR.LIST(1:3).EQ.'X1e') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NX1LBL) LEOOR = .TRUE.
         ILSTSYM = ISYX1(INDEX)
      ELSE IF (LIST(1:2).EQ.'X2' ) THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NX2LBL) LEOOR = .TRUE.
         ILSTSYM = MULD2H(ISYAX2(INDEX),ISYBX2(INDEX))
      ELSE IF (LIST(1:2).EQ.'X3' ) THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NX3LBL) LEOOR = .TRUE.
         ILSTSYM = MULD2H(ISYX3(INDEX,1),ISYX3(INDEX,2))
         ILSTSYM = MULD2H(ILSTSYM,ISYX3(INDEX,3))
      ELSE IF (LIST(1:2).EQ.'X4' ) THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NX4LBL) LEOOR = .TRUE.
         ILSTSYM = MULD2H(ISYX4(INDEX,1),ISYX4(INDEX,2))
         ILSTSYM = MULD2H(ILSTSYM,ISYX4(INDEX,3))
         ILSTSYM = MULD2H(ILSTSYM,ISYX4(INDEX,4))
      ELSE IF (LIST(1:3).EQ.'CX2') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NCX2LBL) LEOOR = .TRUE.
         ILSTSYM = MULD2H(ISYCX2(INDEX,1),ISYCX2(INDEX,2))
C
Cholesky
C
      ELSE IF (LIST(1:3).EQ.'d00') THEN
         ILSTSYM = 1
C
Cholesky
C
      ELSE IF (LIST(1:2).EQ.'L0' .OR. LIST(1:2).EQ.'R0') THEN
         ILSTSYM = 1
      ELSE IF (LIST(1:2).EQ.'D0') THEN
         ILSTSYM = 1
      ELSE IF (LIST(1:2).EQ.'E0' .OR. LIST(1:2).EQ.'BE') THEN
         ILSTSYM = 1
      ELSE IF (LIST(1:2).EQ.'LE' .OR. LIST(1:2).EQ.'RE') THEN
         ILSTSYM = ISYEXC(INDEX)
      ELSE IF (LIST(1:3).EQ.'ER1'.OR.LIST(1:3).EQ.'EO1') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NER1LBL) LEOOR = .TRUE.
         ILSTSYM = MULD2H(ISYSER1(INDEX),ISYOER1(INDEX))
      ELSE IF (LIST(1:3).EQ.'ER2'.OR.LIST(1:3).EQ.'EO2') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NER2LBL) LEOOR = .TRUE.
         ILSTSYM = MULD2H(ISYSER2(INDEX),ISYOER2(INDEX,1))
         ILSTSYM = MULD2H(ILSTSYM,ISYOER2(INDEX,2))
      ELSE IF (LIST(1:3).EQ.'EL1'.OR.LIST(1:3).EQ.'EX1') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NEL1LBL) LEOOR = .TRUE.
         ILSTSYM = MULD2H(ISYSEL1(INDEX),ISYOEL1(INDEX))
      ELSE IF (LIST(1:3).EQ.'EL2'.OR.LIST(1:3).EQ.'EX2') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NEL2LBL) LEOOR = .TRUE.
         ILSTSYM = MULD2H(ISYSEL2(INDEX),ISYOEL2(INDEX,1))
         ILSTSYM = MULD2H(ILSTSYM,ISYOEL2(INDEX,2))
      ELSE IF (LIST(1:2).EQ.'L1'.OR.LIST(1:3).EQ.'X1B') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NLRZLBL) LEOOR = .TRUE.
         ILSTSYM = ISYLRZ(INDEX)
      ELSE IF (LIST(1:2).EQ.'M1'.OR.LIST(1:2).EQ.'FR') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NLRM   ) LEOOR = .TRUE.
         ILSTSYM = ISYLRM(INDEX)
Cholesky
Chol  ELSE IF (LIST(1:2).EQ.'R1' .OR. LIST(1:2).EQ.'F1') THEN
      ELSE IF (LIST(1:2).EQ.'R1' .OR. LIST(1:2).EQ.'F1' .OR.
     &         LIST(1:3).EQ.'XF1' .OR. LIST(1:3).EQ.'d01' .OR.
     &         LIST(1:3).EQ.'eO1') THEN
Cholesky
         IF (INDEX.LE.0 .OR. INDEX.GT.NLRTLBL) LEOOR = .TRUE.
         ILSTSYM = ISYLRT(INDEX)
      ELSE IF ((LIST(1:2).EQ.'RC').OR.(LIST(1:2).EQ.'FC')) THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NLRCLBL) LEOOR = .TRUE.
         ILSTSYM = ISYLRC(INDEX)
      ELSE IF ((LIST(1:2).EQ.'LC').OR.(LIST(1:2).EQ.'XC')) THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NLC1LBL) LEOOR = .TRUE.
         ILSTSYM = ISYLC1(INDEX)
      ELSE IF (LIST(1:3).EQ.'CR2'.OR.LIST(1:3).EQ.'CF2') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NCR2LBL) LEOOR = .TRUE.
         ILSTSYM = MULD2H(ISYCR2(INDEX,1),ISYCR2(INDEX,2))
      ELSE IF (LIST(1:2).EQ.'R2' .OR. LIST(1:2).EQ.'F2') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NR2TLBL) LEOOR = .TRUE.
         ILSTSYM = MULD2H(ISYAR2T(INDEX),ISYBR2T(INDEX))
      ELSE IF (LIST(1:2).EQ.'R3' .OR. LIST(1:2).EQ.'F3') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NR3TLBL) LEOOR = .TRUE.
         ILSTSYM = MULD2H(ISYR3T(INDEX,1),ISYR3T(INDEX,2))
         ILSTSYM = MULD2H(ILSTSYM,ISYR3T(INDEX,3))
      ELSE IF (LIST(1:2).EQ.'R4' .OR. LIST(1:2).EQ.'F4') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NR4TLBL) LEOOR = .TRUE.
         ILSTSYM = MULD2H(ISYR4T(INDEX,1),ISYR4T(INDEX,2))
         ILSTSYM = MULD2H(ILSTSYM,ISYR4T(INDEX,3))
         ILSTSYM = MULD2H(ILSTSYM,ISYR4T(INDEX,4))
      ELSE IF (LIST(1:2).EQ.'L2') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NL2LBL) LEOOR = .TRUE.
         ILSTSYM = MULD2H(ISYAL2(INDEX),ISYBL2(INDEX))
      ELSE IF (LIST(1:2).EQ.'L3') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NL3LBL) LEOOR = .TRUE.
         ILSTSYM = MULD2H(ISYL3(INDEX,1),ISYL3(INDEX,2))
         ILSTSYM = MULD2H(ILSTSYM,ISYL3(INDEX,3))
      ELSE IF (LIST(1:2).EQ.'L4') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NL4LBL) LEOOR = .TRUE.
         ILSTSYM = MULD2H(ISYL4(INDEX,1),ISYL4(INDEX,2))
         ILSTSYM = MULD2H(ILSTSYM,ISYL4(INDEX,3))
         ILSTSYM = MULD2H(ILSTSYM,ISYL4(INDEX,4))
      ELSE IF (LIST(1:3).EQ.'CL2') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NCL2LBL) LEOOR = .TRUE.
         ILSTSYM = MULD2H(ISYCL2(INDEX,1),ISYCL2(INDEX,2))
      ELSE IF (LIST(1:2).EQ.'N2' .OR. LIST(1:2).EQ.'BR') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NQRN2 ) LEOOR = .TRUE.
         ILSTSYM = MULD2H(ISYIN2(INDEX),ISYFN2(INDEX))
      ELSE IF (LIST(1:3).EQ.'PL1') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NPL1LBL) LEOOR = .TRUE.
         ILSTSYM = ISYPL1(INDEX)
      ELSE
         WRITE(LUPRI,*) 'Unknown LIST in ILSTSYM:"',LIST(1:3),'"'
C to force a core dump:
C        WRITE (LUPRI,*) LIST(999999:999999)
         CALL QUIT('Unknown LIST in ILSTSYM.')
      END IF

      IF (LEOOR) THEN  
        WRITE (LUPRI,*) 'INDEX out of range in ILSTSYM:'
        WRITE (LUPRI,*) 'LIST,INDEX:',LIST(1:3),INDEX
C to force a core dump:
C        WRITE (LUPRI,*) LIST(-999999:-999999)
        CALL QUIT('INDEX out of range in ILSTSYM.')
      END IF

      IF (ILSTSYM.LT.1 .OR. ILSTSYM.GT.NSYM) THEN
        NWARN = NWARN + 1
        WRITE (LUPRI,*) 'WARNING from ILSTSYM: symmetry out of range:'
        WRITE (LUPRI,*) 'LIST,INDEX,ILSTSYM:',LIST(1:3),INDEX,ILSTSYM
      END IF

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck ilstsymrlx */
*=====================================================================*
      INTEGER FUNCTION ILSTSYMRLX(LIST,INDEX)
*---------------------------------------------------------------------*
* PURPOSE: get symmetry for orbital relaxation vector on list
*
*          LIST : list type
*          INDEX: index of the vector on the list
*
* Christof Haettig, November 1996
* PL1 introduced Sonia
*=====================================================================*
      IMPLICIT NONE
#include <priunit.h>
#include <ccorb.h>
#include <ccr1rsp.h>

      CHARACTER*(3) LIST
      INTEGER INDEX
      LOGICAL LEOOR

      LEOOR = .FALSE.

* begin:
      IF (LIST(1:2).EQ.'o1') THEN
         CALL QUIT('Illegal list in ILSTSYMRLX.')
      ELSE IF (LIST(1:2).EQ.'R1') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NLRTHFLBL) LEOOR = .TRUE.
         ILSTSYMRLX = ISYLRTHF(INDEX)
      ELSE
         WRITE(LUPRI,*) 'Unknown LIST in ILSTSYM:"',LIST(1:3),'"'
C to force a core dump:
C        WRITE (LUPRI,*) LIST(999999:999999)
         CALL QUIT('Unknown LIST in ILSTSYM.')
      END IF

      IF (LEOOR) THEN  
        WRITE (LUPRI,*) 'INDEX out of range in ILSTSYMRLX:'
        WRITE (LUPRI,*) 'LIST,INDEX:',LIST(1:3),INDEX
C to force a core dump:
C        WRITE (LUPRI,*) LIST(-999999:-999999)
        CALL QUIT('INDEX out of range in ILSTSYMRLX.')
      END IF

      IF (ILSTSYMRLX.LT.1 .OR. ILSTSYMRLX.GT.NSYM) THEN
        NWARN = NWARN + 1
        WRITE (LUPRI,*) 
     &    'WARNING from ILSTSYMRLX: symmetry out of range:'
        WRITE (LUPRI,*) 'LIST,INDEX,ILSTSYMRLX:',
     &    LIST(1:3),INDEX,ILSTSYMRLX
      END IF

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck freqlst */
*=====================================================================*
      REAL*8  FUNCTION FREQLST(LIST, INDEX)
*---------------------------------------------------------------------*
* PURPOSE: return frequency for vector on list
*
*          LIST : list type
*          INDEX: index of the vector on the list
*
* Christof Haettig, April 2002
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit.h>
#include <ccorb.h>
#include <ccroper.h>
#include <cclrmrsp.h>
#include <ccer1rsp.h>
#include <ccer2rsp.h>
#include <ccel1rsp.h>
#include <ccel2rsp.h>
#include <ccr1rsp.h>
#include <ccr2rsp.h>
#include <ccr3rsp.h>
#include <ccr4rsp.h>
#include <ccx1rsp.h>
#include <ccx2rsp.h>
#include <ccx3rsp.h>
#include <ccx4rsp.h>
#include <cco1rsp.h>
#include <cco2rsp.h>
#include <cco3rsp.h>
#include <cco4rsp.h>
#include <ccl1rsp.h>
#include <ccl2rsp.h>
#include <ccl3rsp.h>
#include <ccl4rsp.h>
#include <ccn2rsp.h>
#include <ccrc1rsp.h>
#include <cclc1rsp.h>
#include <cccr2rsp.h>
#include <ccco2rsp.h>
#include <cccx2rsp.h>
#include <cccl2rsp.h>
#include <ccexci.h>
#include <ccpl1rsp.h>

      CHARACTER*(3) LIST
      INTEGER INDEX
      LOGICAL LEOOR

      LEOOR = .FALSE.

* begin:
      IF (LIST(1:2).EQ.'o1') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NRSOLBL) LEOOR = .TRUE.
         CALL QUIT('Illegal list in function FREQLST: '//LIST)
      ELSE IF (LIST(1:3).EQ.'O1 ') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NO1LBL) LEOOR = .TRUE.
         FREQLST = FRQO1(INDEX)
      ELSE IF (LIST(1:2).EQ.'O2' ) THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NO2LBL) LEOOR = .TRUE.
         FREQLST = FRQO2(INDEX,1) + FRQO2(INDEX,2)
      ELSE IF (LIST(1:2).EQ.'O3' ) THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NO3LBL) LEOOR = .TRUE.
         FREQLST = FRQO3(INDEX,1) + FRQO3(INDEX,2) + FRQO3(INDEX,3)
      ELSE IF (LIST(1:2).EQ.'O4' ) THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NO4LBL) LEOOR = .TRUE.
         FREQLST = FRQO4(INDEX,1) + FRQO4(INDEX,2) + 
     &             FRQO4(INDEX,3) + FRQO4(INDEX,4)
      ELSE IF (LIST(1:3).EQ.'CO2') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NCO2LBL) LEOOR = .TRUE.
         FREQLST = 0.0D0
      ELSE IF (LIST(1:3).EQ.'X1 ') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NX1LBL) LEOOR = .TRUE.
         FREQLST = FRQX1(INDEX)
      ELSE IF (LIST(1:2).EQ.'X2' ) THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NX2LBL) LEOOR = .TRUE.
         FREQLST = FRQX2(INDEX,1) + FRQX2(INDEX,2)
      ELSE IF (LIST(1:2).EQ.'X3' ) THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NX3LBL) LEOOR = .TRUE.
         FREQLST = FRQX3(INDEX,1) + FRQX3(INDEX,2) + FRQX3(INDEX,3)
      ELSE IF (LIST(1:2).EQ.'X4' ) THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NX4LBL) LEOOR = .TRUE.
         FREQLST = FRQX4(INDEX,1) + FRQX4(INDEX,2) + 
     &             FRQX4(INDEX,3) + FRQX4(INDEX,4)
      ELSE IF (LIST(1:3).EQ.'CX2') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NCX2LBL) LEOOR = .TRUE.
         FREQLST = 0.0D0
      ELSE IF (LIST(1:2).EQ.'L0' .OR. LIST(1:2).EQ.'R0') THEN
         FREQLST = 0.0D0
      ELSE IF (LIST(1:2).EQ.'D0') THEN
         CALL QUIT('Illegal list in function FREQLST: '//LIST)
      ELSE IF (LIST(1:2).EQ.'E0' .OR. LIST(1:2).EQ.'BE') THEN
         FREQLST = 0.0D0
      ELSE IF (LIST(1:2).EQ.'LE') THEN
         FREQLST = -EIGVAL(INDEX)
      ELSE IF (LIST(1:2).EQ.'RE') THEN
         FREQLST = +EIGVAL(INDEX)
      ELSE IF (LIST(1:3).EQ.'ER1'.OR.LIST(1:3).EQ.'EO1') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NER1LBL) LEOOR = .TRUE.
         FREQLST =  EIGER1(INDEX) + FRQER1(INDEX)
         WRITE(LUPRI,*)'FUNCTION FREQLST NOT DEBUGED FOR "ER1" VECTORS'
      ELSE IF (LIST(1:3).EQ.'ER2'.OR.LIST(1:3).EQ.'EO2') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NER2LBL) LEOOR = .TRUE.
         FREQLST =  EIGER2(INDEX) + FRQER2(INDEX,1) + FRQER2(INDEX,2)
         WRITE(LUPRI,*)'FUNCTION FREQLST NOT DEBUGED FOR "ER2" VECTORS'
      ELSE IF (LIST(1:3).EQ.'EL1'.OR.LIST(1:3).EQ.'EX1') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NEL1LBL) LEOOR = .TRUE.
         FREQLST =  EIGEL1(INDEX) + FRQEL1(INDEX)
         WRITE(LUPRI,*)'FUNCTION FREQLST NOT DEBUGED FOR "EL1" VECTORS'
      ELSE IF (LIST(1:3).EQ.'EL2'.OR.LIST(1:3).EQ.'EX2') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NEL2LBL) LEOOR = .TRUE.
         FREQLST =  EIGEL2(INDEX) + FRQEL2(INDEX,1) + FRQEL2(INDEX,2)
         WRITE(LUPRI,*)'FUNCTION FREQLST NOT DEBUGED FOR "EL2" VECTORS'
      ELSE IF (LIST(1:2).EQ.'L1'.OR.LIST(1:3).EQ.'X1B') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NLRZLBL) LEOOR = .TRUE.
         FREQLST = FRQLRZ(INDEX)
      ELSE IF (LIST(1:2).EQ.'M1'.OR.LIST(1:2).EQ.'FR') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NLRM   ) LEOOR = .TRUE.
         FREQLST = FRQLRM(INDEX)
         WRITE(LUPRI,*)'FUNCTION FREQLST NOT DEBUGED FOR "M1 " VECTORS'
      ELSE IF (LIST(1:2).EQ.'R1' .OR. LIST(1:2).EQ.'F1') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NLRTLBL) LEOOR = .TRUE.
         FREQLST = FRQLRT(INDEX)
      ELSE IF ((LIST(1:2).EQ.'RC').OR.(LIST(1:2).EQ.'FC')) THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NLRCLBL) LEOOR = .TRUE.
         FREQLST = 0.0D0
      ELSE IF ((LIST(1:2).EQ.'LC').OR.(LIST(1:2).EQ.'XC')) THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NLC1LBL) LEOOR = .TRUE.
         FREQLST = 0.0D0
      ELSE IF (LIST(1:3).EQ.'CR2'.OR.LIST(1:3).EQ.'CF2') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NCR2LBL) LEOOR = .TRUE.
         FREQLST = 0.0D0
      ELSE IF (LIST(1:2).EQ.'R2' .OR. LIST(1:2).EQ.'F2') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NR2TLBL) LEOOR = .TRUE.
         FREQLST = FRQR2T(INDEX,1) + FRQR2T(INDEX,2)
      ELSE IF (LIST(1:2).EQ.'R3' .OR. LIST(1:2).EQ.'F3') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NR3TLBL) LEOOR = .TRUE.
         FREQLST = FRQR3T(INDEX,1) + FRQR3T(INDEX,2) + FRQR3T(INDEX,3)
      ELSE IF (LIST(1:2).EQ.'R4' .OR. LIST(1:2).EQ.'F4') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NR4TLBL) LEOOR = .TRUE.
         FREQLST = FRQR4T(INDEX,1) + FRQR4T(INDEX,2) + 
     &             FRQR4T(INDEX,3) + FRQR4T(INDEX,4)
      ELSE IF (LIST(1:2).EQ.'L2') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NL2LBL) LEOOR = .TRUE.
         FREQLST = FRQL2(INDEX,1) + FRQL2(INDEX,2)
      ELSE IF (LIST(1:2).EQ.'L3') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NL3LBL) LEOOR = .TRUE.
         FREQLST = FRQL3(INDEX,1) + FRQL3(INDEX,2) + FRQL3(INDEX,3)
      ELSE IF (LIST(1:2).EQ.'L4') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NL4LBL) LEOOR = .TRUE.
         FREQLST = FRQL4(INDEX,1) + FRQL4(INDEX,2) + 
     &             FRQL4(INDEX,3) + FRQL4(INDEX,4)
      ELSE IF (LIST(1:3).EQ.'CL2') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NCL2LBL) LEOOR = .TRUE.
         FREQLST = 0.0D0
      ELSE IF (LIST(1:2).EQ.'N2' .OR. LIST(1:2).EQ.'BR') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NQRN2 ) LEOOR = .TRUE.
         FREQLST = -EIGN2(INDEX,1) - EIGN2(INDEX,2)
         WRITE(LUPRI,*)'FUNCTION FREQLST NOT DEBUGGED FOR "N2 " VECTORS'
      ELSE IF (LIST(1:3).EQ.'PL1') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NPL1LBL) LEOOR = .TRUE.
         FREQLST = FRQPL1(INDEX)
      ELSE
         WRITE(LUPRI,*) 'Unknown LIST in FREQLST:"',LIST(1:3),'"'
C to force a core dump:
C        WRITE (LUPRI,*) LIST(999999:999999)
         CALL QUIT('Unknown LIST in function FREQLST: '//LIST)
      END IF

      IF (LEOOR) THEN  
        WRITE (LUPRI,*) 'INDEX out of range in FREQLST:'
        WRITE (LUPRI,*) 'LIST,INDEX: ',LIST(1:3),INDEX
C to force a core dump:
C        WRITE (LUPRI,*) LIST(-999999:-999999)
        CALL QUIT('INDEX out of range in FREQLST.')
      END IF

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck idxsym */
*=====================================================================*
      INTEGER FUNCTION IDXSYM(LIST,ISYM,INDEX)
*---------------------------------------------------------------------*
* PURPOSE: Get symmetry for vector on list and calculate the nr. 
*          relative to the offset. Make new list
*
*          LIST : list type
*          INDEX: index of the vector on the list
*          SYM:   Symmetry of vectors 
*          IDXSYM:Index of the vector on the list 
*                 reduced by symmetry offset.
*
* Christof Haettig, November 1996(ILSTSYM), Ove Christiansen Feb. 1997
* PL1 vectors, Sonia 2000
* Cholesky CC2 vectors, tbp 2003
*=====================================================================*
      IMPLICIT NONE
#include <ccorb.h>
#include <ccroper.h>
#include <ccer1rsp.h>
#include <ccer2rsp.h>
#include <ccel1rsp.h>
#include <ccel2rsp.h>
#include <ccr1rsp.h>
#include <ccr2rsp.h>
#include <ccr3rsp.h>
#include <ccr4rsp.h>
#include <ccx1rsp.h>
#include <ccx2rsp.h>
#include <ccx3rsp.h>
#include <ccx4rsp.h>
#include <ccl1rsp.h>
#include <ccl2rsp.h>
#include <ccl3rsp.h>
#include <ccl4rsp.h>
#include <cco1rsp.h>
#include <cco2rsp.h>
#include <cco3rsp.h>
#include <cco4rsp.h>
#include <ccn2rsp.h>
#include <cclrmrsp.h>
#include <ccrc1rsp.h>
#include <cclc1rsp.h>
#include <cccr2rsp.h>
#include <ccco2rsp.h>
#include <cccl2rsp.h>
#include <cccx2rsp.h>
#include <ccexci.h>
#include <ccpl1rsp.h>
#include <priunit.h>

      CHARACTER*(*) LIST
      INTEGER INDEX,ISYM
      LOGICAL LEOOR

      LEOOR = .FALSE.

* begin:
      IF (LIST(1:2).EQ.'L0' .OR. LIST(1:2).EQ.'R0') THEN
         IDXSYM  = 1
      ELSE IF (LIST(1:2).EQ.'D0') THEN
         IDXSYM  = 1
Cholesky
      ELSE IF (LIST(1:3).EQ.'d00') THEN
         IDXSYM = 1
Chol  ELSE IF (LIST(1:2).EQ.'R1' .OR. LIST(1:2).EQ.'F1') THEN
      ELSE IF (LIST(1:2).EQ.'R1' .OR. LIST(1:2).EQ.'F1' .OR.
     &         LIST(1:3).EQ.'XF1' .OR. LIST(1:3).EQ.'d01' .OR.
     &         LIST(1:3).EQ.'eO1') THEN
Cholesky
         IF (INDEX.LT.0 .OR. INDEX.GT.NLRTLBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFT(ISYM) 
      ELSE IF (LIST(1:2).EQ.'R2' .OR. LIST(1:2).EQ.'F2') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NR2TLBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFT2(ISYM)
      ELSE IF (LIST(1:2).EQ.'R3' .OR. LIST(1:2).EQ.'F3') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NR3TLBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFT3(ISYM)
      ELSE IF (LIST(1:2).EQ.'R4' .OR. LIST(1:2).EQ.'F4') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NR4TLBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFT4(ISYM)
      ELSE IF (LIST(1:3).EQ.'O1 '.OR.LIST(1:3).EQ.'O1e') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NO1LBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFO1(ISYM)
      ELSE IF (LIST(1:2).EQ.'O2') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NO2LBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFO2(ISYM)
      ELSE IF (LIST(1:2).EQ.'O3') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NO3LBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFO3(ISYM)
      ELSE IF (LIST(1:2).EQ.'O4') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NO4LBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFO4(ISYM)
      ELSE IF (LIST(1:3).EQ.'CO2') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NCO2LBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFCO2(ISYM) 
      ELSE IF (LIST(1:3).EQ.'X1 '.OR.LIST(1:3).EQ.'X1e') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NX1LBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFX1(ISYM)
      ELSE IF (LIST(1:2).EQ.'X2') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NX2LBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFX2(ISYM)
      ELSE IF (LIST(1:2).EQ.'X3') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NX3LBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFX3(ISYM)
      ELSE IF (LIST(1:2).EQ.'X4') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NX4LBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFX4(ISYM)
      ELSE IF (LIST(1:3).EQ.'CX2') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NCX2LBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFCX2(ISYM) 
      ELSE IF (LIST(1:2).EQ.'L1') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NLRZLBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFZ(ISYM) 
      ELSE IF (LIST(1:2).EQ.'L2') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NL2LBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFL2(ISYM)
      ELSE IF (LIST(1:2).EQ.'L3') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NL3LBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFL3(ISYM)
      ELSE IF (LIST(1:2).EQ.'L4') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NL4LBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFL4(ISYM)
      ELSE IF (LIST(1:3).EQ.'CL2') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NCL2LBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFCL2(ISYM) 
      ELSE IF (LIST(1:2).EQ.'E0' .OR. LIST(1:2).EQ.'BE') THEN
         IDXSYM  = INDEX
      ELSE IF (LIST(1:2).EQ.'LE' .OR. LIST(1:2).EQ.'RE') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NEXCI) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFE(ISYM) 
      ELSE IF (LIST(1:3).EQ.'ER1'.OR.LIST(1:3).EQ.'EO1') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NER1LBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFER1(ISYM) 
      ELSE IF (LIST(1:3).EQ.'ER2'.OR.LIST(1:3).EQ.'EO2') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NER2LBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFER2(ISYM) 
      ELSE IF (LIST(1:3).EQ.'EL1'.OR.LIST(1:3).EQ.'EX1') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NEL1LBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFEL1(ISYM) 
      ELSE IF (LIST(1:3).EQ.'EL2'.OR.LIST(1:3).EQ.'EX2') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NEL2LBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFEL2(ISYM) 
      ELSE IF (LIST(1:2).EQ.'M1' .OR. LIST(1:2).EQ.'FR') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NLRM   ) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFM(ISYM) 
      ELSE IF (LIST(1:2).EQ.'RC' .OR. LIST(1:2).EQ.'FC') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NLRCLBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFC(ISYM) 
      ELSE IF (LIST(1:2).EQ.'LC') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NLC1LBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFLC1(ISYM) 
      ELSE IF (LIST(1:3).EQ.'CR2'.OR.LIST(1:3).EQ.'CF2') THEN
         IF (INDEX.LE.0 .OR. INDEX.GT.NCR2LBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFCR2(ISYM) 
      ELSE IF (LIST(1:2).EQ.'N2'.OR.LIST(1:2).EQ.'BR') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NQRN2 ) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFN2(ISYM)
      !PL1 vectors indices within symmetry class (Sonia)
      ELSE IF (LIST(1:3).EQ.'PL1') THEN
         IF (INDEX.LT.0 .OR. INDEX.GT.NPL1LBL) LEOOR = .TRUE.
         IDXSYM  = INDEX - ISYOFPL1(ISYM)
      ELSE
         WRITE (LUPRI,*) 'Unknown LIST in IDXSYM:"',LIST(1:3),'".'
C to force a core dump:
         WRITE (LUPRI,*) 'core dump:',LIST(9999999:9999999)
         CALL QUIT('Unknown LIST in IDXSYM.')
      END IF

      IF (LEOOR) THEN  
        WRITE (LUPRI,*) 'INDEX out of range in IDXSYM:'
        WRITE (LUPRI,*) 'LIST,INDEX:',LIST,INDEX
        CALL QUIT('INDEX out of range in IDXSYM.')
      END IF
 
C     write(LUPRI,*) 'index,idxsym',index,idxsym

      RETURN
      END
*=====================================================================*
C  /* Deck irhsr2 */
      INTEGER FUNCTION IRHSR2(NEWLBLA,LORXA,FRQANEW,ISYMA,
     *                        NEWLBLB,LORXB,FRQBNEW,ISYMB )
*---------------------------------------------------------------------*
C
C maintain the list right hand side vectors for the 
C second-order coupled cluster amplitude equations
C
C   if vector is on the list return list index and set ISYMA,ISYMB
C   if vector is NOT on the list:
C        LO2OPN=.true.  --> extend list, and return index
C        LO2OPN=.false. --> return -1   
C
C        NEWLBLA / NEWLBLB -- operator labels
C        LORXA   / LORXB   -- flags for orbital relaxation
C        FRQANEW / FRQBNEW -- frequencies 
C        ISYMA   / ISYMB   -- symmetries
C
C Christof Haettig, April 97
C LORXA, LORXB flags introduced in July 1999
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <cco2rsp.h>
#include <priunit.h>
C
      INTEGER ISYMA, ISYMB
      REAL*8  FRQANEW,FRQBNEW,TOL

      PARAMETER(TOL=1.0D-12)

      CHARACTER*8 NEWLBLA, NEWLBLB
      LOGICAL LORXA, LORXB
      INTEGER I

      DO I = 1,NO2LBL
         IF ( NEWLBLA.EQ.LBLAO2(I).AND. NEWLBLB.EQ.LBLBO2(I)
     *       .AND. (LORXA .EQV. LORXAO2(I))
     *       .AND. (LORXB .EQV. LORXBO2(I))
     *       .AND. (ABS(FRQANEW-FRQAO2(I)).LT.TOL)
     *       .AND. (ABS(FRQBNEW-FRQBO2(I)).LT.TOL) 
     *      ) THEN
            IRHSR2 = I
            ISYMA  = ISYAO2(IRHSR2)
            ISYMB  = ISYBO2(IRHSR2)
            RETURN
         END IF
         IF ( NEWLBLB.EQ.LBLAO2(I).AND. NEWLBLA.EQ.LBLBO2(I)
     *       .AND. (LORXB .EQV. LORXAO2(I))
     *       .AND. (LORXA .EQV. LORXBO2(I))
     *       .AND. (ABS(FRQBNEW-FRQAO2(I)).LT.TOL)
     *       .AND. (ABS(FRQANEW-FRQBO2(I)).LT.TOL) 
     *      ) THEN
            IRHSR2 = I
            ISYMB  = ISYAO2(IRHSR2)
            ISYMA  = ISYBO2(IRHSR2)
            RETURN
         END IF
      END DO  

      IF (LO2OPN) THEN
        NO2LBL = NO2LBL + 1

        IF (NO2LBL.GT.MAXO2LBL) THEN
          WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *    '@ NUMBER OF SPECIFIED VECTORS EXCEED THE MAXIMUM ALLOWED',
     *    '@ MAXO2LBL =',MAXO2LBL,' NO2LBL= ',NO2LBL
          CALL QUIT(' IRHSR2: TOO MANY VECTORS SPECIFIED')
        END IF

        LBLAO2(NO2LBL)  = NEWLBLA
        LBLBO2(NO2LBL)  = NEWLBLB
        LORXAO2(NO2LBL) = LORXA
        LORXBO2(NO2LBL) = LORXB
        FRQAO2(NO2LBL)  = FRQANEW
        FRQBO2(NO2LBL)  = FRQBNEW
        ISYAO2(NO2LBL)  = ISYMA
        ISYBO2(NO2LBL)  = ISYMB
        IRHSR2 = NO2LBL

      ELSE
        WRITE(LUPRI,'(3A,L2,A,1P,D12.5,3A,L2,A,1P,D12.5,2A)') 
     *   '@ WARNING: RHSR2 VECTOR FOR ',
     *            NEWLBLA,'(',LORXA,',',FRQANEW,'), ',
     *            NEWLBLB,'(',LORXB,',',FRQBNEW,')',
     *              ' IS NOT AVAILABLE.'
        IRHSR2 = -1
      END IF

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck ichi2 */
      INTEGER FUNCTION ICHI2(NEWLBLA,LORXA,FRQANEW,ISYMA,
     *                       NEWLBLB,LORXB,FRQBNEW,ISYMB )
*---------------------------------------------------------------------*
C
C maintain the list of second-order chi vectors:
C
C   if vector is on the list return list index and set ISYMA,ISYMB
C   if vector is NOT on the list:
C        LX2OPN=.true.  --> extend list, and return index
C        LX2OPN=.false. --> return -1   
C
C Christof Haettig, April 97
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <ccx2rsp.h>
#include <priunit.h>
C
      INTEGER ISYMA, ISYMB
      LOGICAL LORXA, LORXB
      REAL*8  FRQANEW,FRQBNEW,TOL

      PARAMETER(TOL=1.0D-12)

      CHARACTER*8 NEWLBLA, NEWLBLB
      INTEGER I

      DO I = 1,NX2LBL
         IF ( NEWLBLA.EQ.LBLAX2(I).AND. NEWLBLB.EQ.LBLBX2(I)
     *       .AND. (LORXA .EQV. LORXAX2(I))
     *       .AND. (LORXB .EQV. LORXBX2(I))
     *       .AND. (ABS(FRQANEW-FRQAX2(I)).LT.TOL)
     *       .AND. (ABS(FRQBNEW-FRQBX2(I)).LT.TOL) 
     *      ) THEN
            ICHI2  = I
            ISYMA  = ISYAX2(ICHI2)
            ISYMB  = ISYBX2(ICHI2)
            RETURN
         END IF
         IF ( NEWLBLB.EQ.LBLAX2(I).AND. NEWLBLA.EQ.LBLBX2(I)
     *       .AND. (LORXB .EQV. LORXAX2(I))
     *       .AND. (LORXA .EQV. LORXBX2(I))
     *       .AND. (ABS(FRQBNEW-FRQAX2(I)).LT.TOL)
     *       .AND. (ABS(FRQANEW-FRQBX2(I)).LT.TOL) 
     *      ) THEN
            ICHI2  = I
            ISYMB  = ISYAX2(ICHI2)
            ISYMA  = ISYBX2(ICHI2)
            RETURN
         END IF
      END DO  

      IF (LX2OPN) THEN
        NX2LBL = NX2LBL + 1

        IF (NX2LBL.GT.MAXX2LBL) THEN
          WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *    '@ NUMBER OF SPECIFIED VECTORS EXCEED THE MAXIMUM ALLOWED',
     *    '@ MAXX2LBL =',MAXX2LBL,' NX2LBL= ',NX2LBL
          CALL QUIT(' ICHI2: TOO MANY VECTORS SPECIFIED')
        END IF

        LBLAX2(NX2LBL)  = NEWLBLA
        LBLBX2(NX2LBL)  = NEWLBLB
        LORXAX2(NX2LBL) = LORXA
        LORXBX2(NX2LBL) = LORXB
        FRQAX2(NX2LBL)  = FRQANEW
        FRQBX2(NX2LBL)  = FRQBNEW
        ISYAX2(NX2LBL)  = ISYMA
        ISYBX2(NX2LBL)  = ISYMB
        ICHI2 = NX2LBL

      ELSE
        WRITE(LUPRI,'(3A,L2,A,1P,D12.5,3A,L2,A,1P,D12.5,2A)') 
     *   '@ WARNING: X2 VECTOR FOR ',
     *            NEWLBLA,'(',LORXA,',',FRQANEW,'), ',
     *            NEWLBLB,'(',LORXB,',',FRQBNEW,')',
     *              ' IS NOT AVAILABLE.'
        ICHI2 = -1
      END IF

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck il2zeta */
      INTEGER FUNCTION IL2ZETA(NEWLBLA,FRQANEW,ISYMA,
     *                         NEWLBLB,FRQBNEW,ISYMB )
*---------------------------------------------------------------------*
C
C maintain the list of second-order lagrangian multiplier vectors:
C
C   if vector is on the list return list index and set ISYMA,ISYMB
C   if vector is NOT on the list:
C        LL2OPN=.true.  --> extend list, and return index
C        LL2OPN=.false. --> return -1   
C
C Christof Haettig, April 97
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <ccl2rsp.h>
#include <priunit.h>
C
      INTEGER ISYMA, ISYMB
      REAL*8  FRQANEW,FRQBNEW,TOL

      PARAMETER(TOL=1.0D-12)

      CHARACTER*8 NEWLBLA, NEWLBLB
      INTEGER I

      DO I = 1,Nl2LBL
         IF ( NEWLBLA.EQ.LBLAL2(I).AND. NEWLBLB.EQ.LBLBL2(I)
     *       .AND. (ABS(FRQANEW-FRQAL2(I)).LT.TOL)
     *       .AND. (ABS(FRQBNEW-FRQBL2(I)).LT.TOL) 
     *      ) THEN
            IL2ZETA  = I
            ISYMA    = ISYAL2(IL2ZETA)
            ISYMB    = ISYBL2(IL2ZETA)
            RETURN
         END IF
         IF ( NEWLBLB.EQ.LBLAL2(I).AND. NEWLBLA.EQ.LBLBL2(I)
     *       .AND. (ABS(FRQBNEW-FRQAL2(I)).LT.TOL)
     *       .AND. (ABS(FRQANEW-FRQBL2(I)).LT.TOL) 
     *      ) THEN
            IL2ZETA  = I
            ISYMB    = ISYAL2(IL2ZETA)
            ISYMA    = ISYBL2(IL2ZETA)
            RETURN
         END IF
      END DO  

      IF (LL2OPN) THEN
        NL2LBL = NL2LBL + 1

        IF (NL2LBL.GT.MAXL2LBL) THEN
          WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *    '@ NUMBER OF SPECIFIED VECTORS EXCEED THE MAXIMUM ALLOWED',
     *    '@ MAXL2LBL =',MAXL2LBL,' NL2LBL= ',NL2LBL
          CALL QUIT(' IL2ZETA: TOO MANY VECTORS SPECIFIED')
        END IF

        LBLAL2(NL2LBL) = NEWLBLA
        LBLBL2(NL2LBL) = NEWLBLB
        FRQAL2(NL2LBL) = FRQANEW
        FRQBL2(NL2LBL) = FRQBNEW
        ISYAL2(NL2LBL) = ISYMA
        ISYBL2(NL2LBL) = ISYMB
        IL2ZETA = NL2LBL

      ELSE
        WRITE(LUPRI,'(3A,1P,D12.5,3A,1P,D12.5,2A)') 
     *   '@ WARNING: L2 VECTOR FOR ',NEWLBLA,'(',FRQANEW,'), ',
     *                               NEWLBLB,'(',FRQBNEW,')',
     *              ' IS NOT AVAILABLE.'
        IL2ZETA = -1
      END IF

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck iveclist */
      INTEGER FUNCTION IVECLIST(LABELN,FREQN,ISYMN,
     *                          LABELL,FREQL,ISYML,
     *                          TYPE,ORDER,MAXLIST,NLIST,LOPEN)
*---------------------------------------------------------------------*
C maintain a list of response vectors:
C
C   new vector specified by LABELN, FREQN, ISYMN
C   vector list LABELL, FREQL, ISYML
C
C   if vector is on the list return list index
C   if vector is NOT on the list:
C        LOPEN=.true.  --> extend list, and return index
C        LOPEN=.false. --> do not extend the list, but return -1   
C
C Christof Haettig, maj 97
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit.h>
C
      LOGICAL LOPEN, CHANGES, SWAP, NOSWAP, LFOUND
      CHARACTER*(*) TYPE
      INTEGER ORDER, MAXLIST, NLIST, IERR, ILIST, ISY, IOP
      INTEGER ISYMN(ORDER), ISYML(MAXLIST,ORDER)
      REAL*8  FREQN(ORDER), FREQL(MAXLIST,ORDER), FRQ

      CHARACTER*8 LABELN(ORDER), LABELL(MAXLIST,ORDER), LBL
      INTEGER I

*---------------------------------------------------------------------*
* sort after labels, frequencies and symmetries:
*---------------------------------------------------------------------*
      CHANGES = .TRUE.
      DO WHILE (CHANGES)

        CHANGES = .FALSE.

        DO IOP = 1, ORDER-1

          SWAP   = .FALSE.
          NOSWAP = .FALSE.

          DO I = 1, 8
            IF (LGT(LABELN(IOP)(I:I),LABELN(IOP+1)(I:I))
     &          .AND. .NOT. NOSWAP) SWAP = .TRUE.
            IF (LLT(LABELN(IOP)(I:I),LABELN(IOP+1)(I:I))
     &          .AND. .NOT. SWAP) NOSWAP = .TRUE.
          END DO

          IF (FREQN(IOP).GT.FREQN(IOP+1)
     &          .AND. .NOT. NOSWAP) SWAP = .TRUE.
          IF (FREQN(IOP).LT.FREQN(IOP+1)
     &          .AND. .NOT. SWAP) NOSWAP = .TRUE.

          IF (ISYMN(IOP).GT.ISYMN(IOP+1)
     &          .AND. .NOT. NOSWAP) SWAP = .TRUE.
          IF (ISYMN(IOP).LT.ISYMN(IOP+1)
     &          .AND. .NOT. SWAP) NOSWAP = .TRUE.

          IF (SWAP) THEN
            CHANGES = .TRUE.
            LBL           = LABELN(IOP)
            LABELN(IOP)   = LABELN(IOP+1)
            LABELN(IOP+1) = LBL
            FRQ           = FREQN(IOP)
            FREQN(IOP)    = FREQN(IOP+1)
            FREQN(IOP+1)  = FRQ
            ISY           = ISYMN(IOP)
            ISYMN(IOP)    = ISYMN(IOP+1)
            ISYMN(IOP+1)  = ISY
          END IF

        END DO
      END DO

*---------------------------------------------------------------------*
* search list for vector:
*---------------------------------------------------------------------*
      DO ILIST = 1, NLIST

        LFOUND = .TRUE.
        DO IOP = 1, ORDER
          IF ( LABELN(IOP) .NE. LABELL(ILIST,IOP) ) LFOUND = .FALSE.
          IF ( FREQN(IOP)  .NE. FREQL(ILIST,IOP)  ) LFOUND = .FALSE.
        END DO

        IF (LFOUND) THEN
          DO IOP = 1, ORDER
            ISYMN(IOP) = ISYML(ILIST,IOP)
          END DO
          IVECLIST = ILIST
          IERR = 0
          RETURN
        END IF

      END DO  

      IF (LOPEN) THEN
        NLIST = NLIST + 1

        IF (NLIST.GT.MAXLIST) THEN
          WRITE(LUPRI,'(4A,/A,I5,A,I5)')
     *    'NUMBER OF SPECIFIED VECTORS FOR THE ',TYPE,'-VECTOR LIST ',
     *    'EXCEED THE ALLOWED MAXIMUM.',
     *    'MAXIMUM =',MAXLIST,'   ---    SPECIFIED = ',NLIST
          CALL QUIT(' IVECLIST: TOO MANY '//TYPE(1:3)
     *                                         //'-VECTORS SPECIFIED')
        END IF

        DO IOP = 1, ORDER
          ISYML(NLIST,IOP)  = ISYMN(IOP)
          FREQL(NLIST,IOP)  = FREQN(IOP)
          LABELL(NLIST,IOP) = LABELN(IOP)
        END DO

        IVECLIST = NLIST
        IERR     = 0

      ELSE

        WRITE(LUPRI,'(2A,2(3A,1P,D12.5))') 
     *   'WARNING: ',TYPE,'-VECTOR FOR ',
     *     (LABELN(IOP), '(', FREQN(IOP), '), ', IOP=1, ORDER)
        WRITE(LUPRI,'(A)') ' IS NOT AVAILABLE.'
        IVECLIST = -1

      END IF

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck il4zeta */
      INTEGER FUNCTION IL4ZETA(LBLA,FRQA,ISYA,LBLB,FRQB,ISYB,
     *                         LBLC,FRQC,ISYC,LBLD,FRQD,ISYD)
*---------------------------------------------------------------------*
C maintain the list of fourth-order lagrangian multiplier vectors
C Christof Haettig, maj 97
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <ccl4rsp.h>
      CHARACTER*2 TYPE
      INTEGER ORDER
      PARAMETER (TYPE='L4', ORDER=4)
 
      INTEGER ISYA, ISYB, ISYC, ISYD
      INTEGER ISYM(4)
      REAL*8  FRQA,FRQB,FRQC,FRQD
      REAL*8  FREQ(4)
      CHARACTER*8 LBLA, LBLB, LBLC, LBLD
      CHARACTER*8 LABEL(4)
     
* external function:
      INTEGER IVECLIST

      ISYM(1) = ISYA
      ISYM(2) = ISYB
      ISYM(3) = ISYC
      ISYM(4) = ISYD

      FREQ(1) = FRQA
      FREQ(2) = FRQB
      FREQ(3) = FRQC
      FREQ(4) = FRQD

      LABEL(1) = LBLA
      LABEL(2) = LBLB
      LABEL(3) = LBLC
      LABEL(4) = LBLD

      IL4ZETA = IVECLIST(LABEL,FREQ,ISYM,LBLL4,FRQL4,ISYL4,
     &                   TYPE,ORDER,MAXL4LBL,NL4LBL,LL4OPN)

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck ir4tamp */
      INTEGER FUNCTION IR4TAMP(LBLA,FRQA,ISYA,LBLB,FRQB,ISYB,
     *                         LBLC,FRQC,ISYC,LBLD,FRQD,ISYD)
*---------------------------------------------------------------------*
C maintain the list of fourth-order amplitude response vectors
C Christof Haettig, maj 97
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <ccr4rsp.h>
      CHARACTER*2 TYPE
      INTEGER ORDER
      PARAMETER (TYPE='R4',ORDER=4)
 
      INTEGER ISYA, ISYB, ISYC, ISYD
      INTEGER ISYM(4)
      REAL*8  FRQA,FRQB,FRQC,FRQD
      REAL*8  FREQ(4)
      CHARACTER*8 LBLA, LBLB, LBLC, LBLD
      CHARACTER*8 LABEL(4)

* external function:
      INTEGER IVECLIST

      ISYM(1) = ISYA
      ISYM(2) = ISYB
      ISYM(3) = ISYC
      ISYM(4) = ISYD

      FREQ(1) = FRQA
      FREQ(2) = FRQB
      FREQ(3) = FRQC
      FREQ(4) = FRQD

      LABEL(1) = LBLA
      LABEL(2) = LBLB
      LABEL(3) = LBLC
      LABEL(4) = LBLD

      IR4TAMP = IVECLIST(LABEL,FREQ,ISYM,LBLR4T,FRQR4T,ISYR4T,
     &                   TYPE,ORDER,MAXT4LBL,NR4TLBL,LR4OPN)

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck irhsr4 */
      INTEGER FUNCTION IRHSR4(LBLA,FRQA,ISYA,LBLB,FRQB,ISYB,
     *                        LBLC,FRQC,ISYC,LBLD,FRQD,ISYD)
*---------------------------------------------------------------------*
C maintain the list of fourth-order amplitude right hand side vectors
C Christof Haettig, maj 97
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <cco4rsp.h>
      CHARACTER*2 TYPE
      INTEGER ORDER
      PARAMETER (TYPE='O4',ORDER=4)
 
      INTEGER ISYA, ISYB, ISYC, ISYD
      INTEGER ISYM(4)
      REAL*8  FRQA,FRQB,FRQC,FRQD
      REAL*8  FREQ(4)
      CHARACTER*8 LBLA, LBLB, LBLC, LBLD
      CHARACTER*8 LABEL(4)

* external function:
      INTEGER IVECLIST

      ISYM(1) = ISYA
      ISYM(2) = ISYB
      ISYM(3) = ISYC
      ISYM(4) = ISYD

      FREQ(1) = FRQA
      FREQ(2) = FRQB
      FREQ(3) = FRQC
      FREQ(4) = FRQD

      LABEL(1) = LBLA
      LABEL(2) = LBLB
      LABEL(3) = LBLC
      LABEL(4) = LBLD

      IRHSR4 = IVECLIST(LABEL,FREQ,ISYM,LBLO4,FRQO4,ISYO4,
     &                  TYPE,ORDER,MAXO4LBL,NO4LBL,LO4OPN)

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck ichi4 */
      INTEGER FUNCTION ICHI4(LBLA,FRQA,ISYA,LBLB,FRQB,ISYB,
     *                       LBLC,FRQC,ISYC,LBLD,FRQD,ISYD)
*---------------------------------------------------------------------*
C maintain the list of fourth-order chi vectors
C Christof Haettig, maj 97
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <ccx4rsp.h>
      CHARACTER*2 TYPE
      INTEGER ORDER
      PARAMETER (TYPE='X4',ORDER=4)
 
      INTEGER ISYA, ISYB, ISYC, ISYD
      INTEGER ISYM(4)
      REAL*8  FRQA,FRQB,FRQC,FRQD
      REAL*8  FREQ(4)
      CHARACTER*8 LBLA, LBLB, LBLC, LBLD
      CHARACTER*8 LABEL(4)

* external function:
      INTEGER IVECLIST

      ISYM(1) = ISYA
      ISYM(2) = ISYB
      ISYM(3) = ISYC
      ISYM(4) = ISYD

      FREQ(1) = FRQA
      FREQ(2) = FRQB
      FREQ(3) = FRQC
      FREQ(4) = FRQD

      LABEL(1) = LBLA
      LABEL(2) = LBLB
      LABEL(3) = LBLC
      LABEL(4) = LBLD

      ICHI4 = IVECLIST(LABEL,FREQ,ISYM,LBLX4,FRQX4,ISYX4,
     &                 TYPE,ORDER,MAXX4LBL,NX4LBL,LX4OPN)

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck il3zeta */
      INTEGER FUNCTION IL3ZETA(LBLA,FRQA,ISYA,LBLB,FRQB,ISYB,
     *                         LBLC,FRQC,ISYC)
*---------------------------------------------------------------------*
C maintain the list of third-order lagrangian multiplier vectors
C Christof Haettig, maj 97
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <ccl3rsp.h>
      CHARACTER*2 TYPE
      INTEGER ORDER
      PARAMETER (TYPE='L3', ORDER=3)
 
      INTEGER ISYA, ISYB, ISYC
      INTEGER ISYM(ORDER)
      REAL*8  FRQA,FRQB,FRQC
      REAL*8  FREQ(ORDER)
      CHARACTER*8 LBLA, LBLB, LBLC
      CHARACTER*8 LABEL(ORDER)

* external function:
      INTEGER IVECLIST

      ISYM(1) = ISYA
      ISYM(2) = ISYB
      ISYM(3) = ISYC

      FREQ(1) = FRQA
      FREQ(2) = FRQB
      FREQ(3) = FRQC

      LABEL(1) = LBLA
      LABEL(2) = LBLB
      LABEL(3) = LBLC

      IL3ZETA = IVECLIST(LABEL,FREQ,ISYM,LBLL3,FRQL3,ISYL3,
     &                   TYPE,ORDER,MAXL3LBL,NL3LBL,LL3OPN)

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck ir3tamp */
      INTEGER FUNCTION IR3TAMP(LBLA,FRQA,ISYA,LBLB,FRQB,ISYB,
     *                         LBLC,FRQC,ISYC)
*---------------------------------------------------------------------*
C maintain the list of fourth-order amplitude response vectors
C Christof Haettig, maj 97
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <ccr3rsp.h>
      CHARACTER*2 TYPE
      INTEGER ORDER
      PARAMETER (TYPE='R3',ORDER=3)
 
      INTEGER ISYA, ISYB, ISYC
      INTEGER ISYM(ORDER)
      REAL*8  FRQA,FRQB,FRQC
      REAL*8  FREQ(ORDER)
      CHARACTER*8 LBLA, LBLB, LBLC
      CHARACTER*8 LABEL(ORDER)

* external function:
      INTEGER IVECLIST

      ISYM(1) = ISYA
      ISYM(2) = ISYB
      ISYM(3) = ISYC

      FREQ(1) = FRQA
      FREQ(2) = FRQB
      FREQ(3) = FRQC

      LABEL(1) = LBLA
      LABEL(2) = LBLB
      LABEL(3) = LBLC

      IR3TAMP = IVECLIST(LABEL,FREQ,ISYM,LBLR3T,FRQR3T,ISYR3T,
     &                   TYPE,ORDER,MAXT3LBL,NR3TLBL,LR3OPN)

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck irhsr3 */
      INTEGER FUNCTION IRHSR3(LBLA,FRQA,ISYA,LBLB,FRQB,ISYB,
     *                       LBLC,FRQC,ISYC)
*---------------------------------------------------------------------*
C maintain the list of third-order amplitude right hand side vectors
C Christof Haettig, maj 97
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <cco3rsp.h>
      CHARACTER*2 TYPE
      INTEGER ORDER
      PARAMETER (TYPE='O3',ORDER=3)
 
      INTEGER ISYA, ISYB, ISYC
      INTEGER ISYM(ORDER)
      REAL*8  FRQA,FRQB,FRQC
      REAL*8  FREQ(ORDER)
      CHARACTER*8 LBLA, LBLB, LBLC
      CHARACTER*8 LABEL(ORDER)

* external function:
      INTEGER IVECLIST

      ISYM(1) = ISYA
      ISYM(2) = ISYB
      ISYM(3) = ISYC

      FREQ(1) = FRQA
      FREQ(2) = FRQB
      FREQ(3) = FRQC

      LABEL(1) = LBLA
      LABEL(2) = LBLB
      LABEL(3) = LBLC

      IRHSR3 = IVECLIST(LABEL,FREQ,ISYM,LBLO3,FRQO3,ISYO3,
     &                  TYPE,ORDER,MAXO3LBL,NO3LBL,LO3OPN)

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck ichi3 */
      INTEGER FUNCTION ICHI3(LBLA,FRQA,ISYA,LBLB,FRQB,ISYB,
     *                       LBLC,FRQC,ISYC)
*---------------------------------------------------------------------*
C maintain the list of fourth-order chi vectors
C Christof Haettig, maj 97
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <ccx3rsp.h>
      CHARACTER*2 TYPE
      INTEGER ORDER
      PARAMETER (TYPE='X3',ORDER=3)
 
      INTEGER ISYA, ISYB, ISYC
      INTEGER ISYM(ORDER)
      REAL*8  FRQA,FRQB,FRQC
      REAL*8  FREQ(ORDER)
      CHARACTER*8 LBLA, LBLB, LBLC
      CHARACTER*8 LABEL(ORDER)

* external function:
      INTEGER IVECLIST

      ISYM(1) = ISYA
      ISYM(2) = ISYB
      ISYM(3) = ISYC

      FREQ(1) = FRQA
      FREQ(2) = FRQB
      FREQ(3) = FRQC

      LABEL(1) = LBLA
      LABEL(2) = LBLB
      LABEL(3) = LBLC

      ICHI3 = IVECLIST(LABEL,FREQ,ISYM,LBLX3,FRQX3,ISYX3,
     &                 TYPE,ORDER,MAXX3LBL,NX3LBL,LX3OPN)

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck ier1amp */
      INTEGER FUNCTION IER1AMP(IEXCI,  EIGVNEW,ISYMS,
     *                         NEWLBLA,FRQANEW,ISYMA, LPROJ )
*---------------------------------------------------------------------*
C
C maintain the list of first-order right excited state vectors:
C
C   if vector is on the list return list index and set ISYMS,ISYMA
C   if vector is NOT on the list:
C        LER1OPN=.true.  --> extend list, and return index
C        LER1OPN=.false. --> return -1   
C
C Christof Haettig, july 97
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <ccer1rsp.h>
#include <priunit.h>
C
      LOGICAL LPROJ, LPROJ1
      INTEGER ISYMA, ISYMS, IEXCI
      REAL*8  FRQANEW,EIGVNEW,TOL

      PARAMETER(TOL=1.0D-12)

      CHARACTER*8 NEWLBLA
      INTEGER I

      LPROJ1 = LPROJ

*     for non-total symmetric operators we can ignore projection
      IF (ISYMA.NE.1) LPROJ1 = .FALSE.

      DO I = 1,NER1LBL
         IF ( (NEWLBLA.EQ.LBLER1(I)) .AND. (IEXCI.EQ.ISTER1(I))
     *       .AND. (ABS(FRQANEW-FRQER1(I)).LT.TOL) 
     *       .AND. (ABS(EIGVNEW-EIGER1(I)).LT.TOL)
     *       .AND. (LPROJ1.EQV.LPRER1(I))
     *      ) THEN
            IER1AMP  = I
            ISYMS    = ISYSER1(IER1AMP)
            ISYMA    = ISYOER1(IER1AMP)
            RETURN
         END IF
      END DO  

      IF (LER1OPN) THEN
        NER1LBL = NER1LBL + 1

        IF (NER1LBL.GT.MAXER1LBL) THEN
          WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *    '@ NUMBER OF SPECIFIED VECTORS EXCEED THE MAXIMUM ALLOWED',
     *    '@ MAXER1LBL =',MAXER1LBL,' NER1LBL= ',NER1LBL
          CALL QUIT(' IER1AMP: TOO MANY VECTORS SPECIFIED')
        END IF

        ISTER1(NER1LBL) = IEXCI
        EIGER1(NER1LBL) = EIGVNEW
        ISYSER1(NER1LBL)= ISYMS
        LBLER1(NER1LBL) = NEWLBLA
        FRQER1(NER1LBL) = FRQANEW
        ISYOER1(NER1LBL)= ISYMA
        LPRER1(NER1LBL) = LPROJ1
        IER1AMP = NER1LBL

      ELSE
        WRITE(LUPRI,'(A,I3,A,1P,D12.5,3A,1P,D12.5,2A)') 
     *   '@ WARNING: ER1 VECTOR FOR',IEXCI,'(',EIGVNEW,'), ',
     *                                NEWLBLA,'(',FRQANEW,')',
     *              ' IS NOT AVAILABLE.'
        IER1AMP = -1
      END IF

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck ier2amp */
      INTEGER FUNCTION IER2AMP(IEXCI,  EIGVNEW,ISYMS,
     *                         NEWLBLA,FRQANEW,ISYMA,
     *                         NEWLBLB,FRQBNEW,ISYMB, LPROJ )
*---------------------------------------------------------------------*
C
C maintain the list of second-order right excited state vectors:
C
C   if vector is on the list return list index and set symmetries
C   if vector is NOT on the list:
C        LER2OPN=.true.  --> extend list, and return index
C        LER2OPN=.false. --> return -1   
C
C Christof Haettig, july 97
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <ccer2rsp.h>
#include <priunit.h>
C
      LOGICAL LPROJ, LPROJ1
      INTEGER ISYMA, ISYMB, ISYMS, IEXCI
      REAL*8  FRQANEW,FRQBNEW,EIGVNEW,TOL

      PARAMETER(TOL=1.0D-12)

      CHARACTER*8 NEWLBLA, NEWLBLB
      INTEGER I

      LPROJ1 = LPROJ

*     for non-total symmetric operators we ignore projection
      IF (ISYMA.NE.1 .AND. ISYMB.NE.1 .AND. ISYMB.NE.ISYMA) THEN
        LPROJ1 = .FALSE.
      END IF

      DO I = 1,NER2LBL
         IF ( (IEXCI.EQ.ISTER2(I)) .AND. (ABS(EIGVNEW-EIGER2(I)).LT.TOL)
     *       .AND. (NEWLBLA.EQ.LBLER2(I,1))
     *         .AND. (ABS(FRQANEW-FRQER2(I,1)).LT.TOL) 
     *       .AND. (NEWLBLA.EQ.LBLER2(I,2))
     *         .AND. (ABS(FRQBNEW-FRQER2(I,2)).LT.TOL) 
     *       .AND. (LPROJ1.EQV.LPRER2(I))
     *      ) THEN
            IER2AMP  = I
            ISYMS    = ISYSER2(IER2AMP)
            ISYMA    = ISYOER2(IER2AMP,1)
            ISYMB    = ISYOER2(IER2AMP,2)
            RETURN
         END IF
         IF ( (IEXCI.EQ.ISTER2(I)) .AND. (ABS(EIGVNEW-EIGER2(I)).LT.TOL)
     *       .AND. (NEWLBLA.EQ.LBLER2(I,2))
     *         .AND. (ABS(FRQANEW-FRQER2(I,2)).LT.TOL) 
     *       .AND. (NEWLBLA.EQ.LBLER2(I,1))
     *         .AND. (ABS(FRQBNEW-FRQER2(I,1)).LT.TOL) 
     *       .AND. (LPROJ1.EQV.LPRER2(I))
     *      ) THEN
            IER2AMP  = I
            ISYMS    = ISYSER2(IER2AMP)
            ISYMA    = ISYOER2(IER2AMP,2)
            ISYMB    = ISYOER2(IER2AMP,1)
            RETURN
         END IF
      END DO  

      IF (LER2OPN) THEN
        NER2LBL = NER2LBL + 1

        IF (NER2LBL.GT.MAXER2LBL) THEN
          WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *    '@ NUMBER OF SPECIFIED VECTORS EXCEED THE MAXIMUM ALLOWED',
     *    '@ MAXER2LBL =',MAXER2LBL,' NER2LBL= ',NER2LBL
          CALL QUIT(' IER2AMP: TOO MANY VECTORS SPECIFIED')
        END IF

        ISTER2(NER2LBL)   = IEXCI
        EIGER2(NER2LBL)   = EIGVNEW
        ISYSER2(NER2LBL)  = ISYMS
        LBLER2(NER2LBL,1) = NEWLBLA
        FRQER2(NER2LBL,1) = FRQANEW
        ISYOER2(NER2LBL,1)= ISYMA
        LBLER2(NER2LBL,2) = NEWLBLB
        FRQER2(NER2LBL,2) = FRQBNEW
        ISYOER2(NER2LBL,2)= ISYMB
        LPRER2(NER2LBL)   = LPROJ1
        IER2AMP = NER2LBL

      ELSE
        WRITE(LUPRI,'(A,I3,A,1P,D12.5,2(3A,1P,D12.5),2A)') 
     *   '@ WARNING: ER2 VECTOR FOR',IEXCI,'(',EIGVNEW,'), ',
     *        NEWLBLA,'(',FRQANEW,')', NEWLBLB,'(',FRQBNEW,')',
     *              ' IS NOT AVAILABLE.'
        IER2AMP = -1
      END IF

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck iel1amp */
      INTEGER FUNCTION IEL1AMP(IEXCI,  EIGVNEW,ISYMS,
     *                         NEWLBLA,FRQANEW,ISYMA,LORXA,LPROJ )
*---------------------------------------------------------------------*
C
C maintain the list of first-order right excited state vectors:
C
C   if vector is on the list return list index and set ISYMS,ISYMA
C   if vector is NOT on the list:
C        LEL1OPN=.true.  --> extend list, and return index
C        LEL1OPN=.false. --> return -1   
C
C Christof Haettig, july 97
C LORXA flag introduced, Sonia Coriani april 2000
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <ccel1rsp.h>
#include <priunit.h>
C
      LOGICAL LPROJ, LPROJ1, LORXA
      INTEGER ISYMA, ISYMS, IEXCI
      REAL*8  FRQANEW,EIGVNEW,TOL

      PARAMETER(TOL=1.0D-12)

      CHARACTER*8 NEWLBLA
      INTEGER I

      LPROJ1 = LPROJ

*     for non-total symmetric operators we can ignore projection
      IF (ISYMA.NE.1) LPROJ1 = .FALSE.

      DO I = 1,NEL1LBL
         IF ( (NEWLBLA.EQ.LBLEL1(I)) .AND. (IEXCI.EQ.ISTEL1(I))
     *       .AND. (LORXA .EQV. LORXEL1(I))
     *       .AND. (ABS(FRQANEW-FRQEL1(I)).LT.TOL) 
     *       .AND. (ABS(EIGVNEW-EIGEL1(I)).LT.TOL)
     *       .AND. (LPROJ1.EQV.LPREL1(I))
     *      ) THEN
            IEL1AMP  = I
            ISYMS    = ISYSEL1(IEL1AMP)
            ISYMA    = ISYOEL1(IEL1AMP)
            RETURN
         END IF
      END DO  

      IF (LEL1OPN) THEN
        NEL1LBL = NEL1LBL + 1

        IF (NEL1LBL.GT.MAXEL1LBL) THEN
          WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *    '@ NUMBER OF SPECIFIED VECTORS EXCEED THE MAXIMUM ALLOWED',
     *    '@ MAXEL1LBL =',MAXEL1LBL,' NEL1LBL= ',NEL1LBL
          CALL QUIT(' IEL1AMP: TOO MANY VECTORS SPECIFIED')
        END IF

        ISTEL1(NEL1LBL) = IEXCI
        EIGEL1(NEL1LBL) = EIGVNEW
        ISYSEL1(NEL1LBL)= ISYMS
        LBLEL1(NEL1LBL) = NEWLBLA
        LORXEL1(NEL1LBL) = LORXA
        FRQEL1(NEL1LBL) = FRQANEW
        ISYOEL1(NEL1LBL)= ISYMA
        LPREL1(NEL1LBL) = LPROJ1
        IEL1AMP = NEL1LBL

      ELSE
        WRITE(LUPRI,'(A,I3,A,1P,D12.5,3A,L2,A,1P,D12.5,2A)')
     &   '@ WARNING: EL1 VECTOR FOR',IEXCI,'(',EIGVNEW,'), ',
     &                       NEWLBLA,'(',LORXA,',',FRQANEW,')',
     &              ' IS NOT AVAILABLE.'
        IEL1AMP = -1
      END IF

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck iel2amp */
      INTEGER FUNCTION IEL2AMP(IEXCI,  EIGVNEW,ISYMS,
     *                         NEWLBLA,FRQANEW,ISYMA,
     *                         NEWLBLB,FRQBNEW,ISYMB,LPROJ )
*---------------------------------------------------------------------*
C
C maintain the list of second-order left excited state vectors:
C
C   if vector is on the list return list index and set symmetries
C   if vector is NOT on the list:
C        LEL2OPN=.true.  --> extend list, and return index
C        LEL2OPN=.false. --> return -1   
C
C Christof Haettig, july 97
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <ccel2rsp.h>
#include <priunit.h>
C
      LOGICAL LPROJ, LPROJ1
      INTEGER ISYMA, ISYMB, ISYMS, IEXCI
      REAL*8  FRQANEW,FRQBNEW,EIGVNEW,TOL

      PARAMETER(TOL=1.0D-12)

      CHARACTER*8 NEWLBLA, NEWLBLB
      INTEGER I

      LPROJ1 = LPROJ

*     for non-total symmetric operators we ignore projection
      IF (ISYMA.NE.1 .AND. ISYMB.NE.1 .AND. ISYMB.NE.ISYMA) THEN
        LPROJ1 = .FALSE.
      END IF

      DO I = 1,NEL2LBL
         IF ((IEXCI.EQ.ISTEL2(I)) .AND. (ABS(EIGVNEW-EIGEL2(I)).LT.TOL)
     *       .AND. (NEWLBLA.EQ.LBLEL2(I,1))
     *         .AND. (ABS(FRQANEW-FRQEL2(I,1)).LT.TOL) 
     *       .AND. (NEWLBLA.EQ.LBLEL2(I,2))
     *         .AND. (ABS(FRQBNEW-FRQEL2(I,2)).LT.TOL) 
     *       .AND. (LPROJ1.EQV.LPREL2(I))
     *      ) THEN
            IEL2AMP  = I
            ISYMS    = ISYSEL2(IEL2AMP)
            ISYMA    = ISYOEL2(IEL2AMP,1)
            ISYMB    = ISYOEL2(IEL2AMP,2)
            RETURN
         END IF
         IF ((IEXCI.EQ.ISTEL2(I)) .AND. (ABS(EIGVNEW-EIGEL2(I)).LT.TOL)
     *       .AND. (NEWLBLA.EQ.LBLEL2(I,2))
     *         .AND. (ABS(FRQANEW-FRQEL2(I,2)).LT.TOL) 
     *       .AND. (NEWLBLA.EQ.LBLEL2(I,1))
     *         .AND. (ABS(FRQBNEW-FRQEL2(I,1)).LT.TOL) 
     *       .AND. (LPROJ1.EQV.LPREL2(I))
     *      ) THEN
            IEL2AMP  = I
            ISYMS    = ISYSEL2(IEL2AMP)
            ISYMA    = ISYOEL2(IEL2AMP,2)
            ISYMB    = ISYOEL2(IEL2AMP,1)
            RETURN
         END IF
      END DO  

      IF (LEL2OPN) THEN
        NEL2LBL = NEL2LBL + 1

        IF (NEL2LBL.GT.MAXEL2LBL) THEN
          WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *    '@ NUMBER OF SPECIFIED VECTORS EXCEED THE MAXIMUM ALLOWED',
     *    '@ MAXEL2LBL =',MAXEL2LBL,' NEL2LBL= ',NEL2LBL
          CALL QUIT(' IEL2AMP: TOO MANY VECTORS SPECIFIED')
        END IF

        ISTEL2(NEL2LBL)   = IEXCI
        EIGEL2(NEL2LBL)   = EIGVNEW
        ISYSEL2(NEL2LBL)  = ISYMS
        LBLEL2(NEL2LBL,1) = NEWLBLA
        FRQEL2(NEL2LBL,1) = FRQANEW
        ISYOEL2(NEL2LBL,1)= ISYMA
        LBLEL2(NEL2LBL,2) = NEWLBLB
        FRQEL2(NEL2LBL,2) = FRQBNEW
        ISYOEL2(NEL2LBL,2)= ISYMB
        LPREL2(NEL2LBL)   = LPROJ1
        IEL2AMP = NEL2LBL

      ELSE
        WRITE(LUPRI,'(A,I3,A,1P,D12.5,2(3A,1P,D12.5),2A)') 
     *   '@ WARNING: EL2 VECTOR FOR',IEXCI,'(',EIGVNEW,'), ',
     *        NEWLBLA,'(',FRQANEW,')', NEWLBLB,'(',FRQBNEW,')',
     *              ' IS NOT AVAILABLE.'
        IEL2AMP = -1
      END IF

      RETURN
      END
*=====================================================================*
C  /* Deck In2amp */
      INTEGER FUNCTION IN2AMP(IIEX,FRQINEW,ISYMI,
     *                        IFEX,FRQFNEW,ISYMF )
*---------------------------------------------------------------------*
C
C maintain the list of N(if)(omegai,omegaf) multipliers for calculation
C of quadratic response function residues.
C
C   if vector is on the list return list index and set ISYMI,ISYMF
C   if vector is NOT on the list:
C        LN2OPN=.true.  --> extend list, and return index
C        LN2OPN=.false. --> return -1   
C
C Ove Christiansen, April 97
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <ccn2rsp.h>
#include <priunit.h>
C
      INTEGER ISYMI, ISYMF
      REAL*8  FRQINEW,FRQFNEW,TOL

      PARAMETER(TOL=1.0D-12)

      INTEGER I,IIEX,IFEX

      DO I = 1,NQRN2
         IF (IIEX.EQ.IIN2(I).AND. IFEX.EQ.IFN2(I)
     *       .AND. (ABS(FRQINEW-FRQIN2(I)).LT.TOL)
     *       .AND. (ABS(FRQFNEW-FRQFN2(I)).LT.TOL)) THEN
            IN2AMP   = I
            ISYMI    = ISYIN2(IN2AMP)
            ISYMF    = ISYFN2(IN2AMP)
            RETURN
         END IF
         IF (IFEX.EQ.IIN2(I).AND. IIEX.EQ.IFN2(I)
     *       .AND. (ABS(FRQFNEW-FRQIN2(I)).LT.TOL)
     *       .AND. (ABS(FRQINEW-FRQFN2(I)).LT.TOL)) THEN
            IN2AMP   = I
            ISYMF    = ISYIN2(IN2AMP) 
            ISYMI    = ISYFN2(IN2AMP) 
            RETURN
         END IF
      END DO  

      IF (LN2OPN) THEN
        NQRN2  = NQRN2  + 1

        IF (NQRN2 .GT.MAXQRN2 ) THEN
          WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *    '@ NUMBER OF SPECIFIED VECTORS EXCEED THE MAXIMUM ALLOWED',
     *    '@ MAXQRN2  =',MAXQRN2 ,' NQRN2 = ',NQRN2 
          CALL QUIT(' IN2AMP: TOO MANY VECTORS SPECIFIED')
        END IF

        IIN2(NQRN2 )  = IIEX   
        IFN2(NQRN2 )  = IFEX   
        FRQIN2(NQRN2 ) = FRQINEW
        FRQFN2(NQRN2 ) = FRQFNEW
        ISYIN2(NQRN2 ) = ISYMI
        ISYFN2(NQRN2 ) = ISYMF
        IN2AMP  = NQRN2 

      ELSE
        WRITE(LUPRI,'(1A,I3,A,1P,D12.5,A,I3,A,1P,D12.5,2A)') 
     *   '@ WARNING: N2 VECTOR FOR ',IIEX,'(',FRQINEW,'), ',
     *                               IFEX,'(',FRQFNEW,')',
     *              ' IS NOT AVAILABLE.'
        IN2AMP  = -1
      END IF

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck ICL2AMP */
      INTEGER FUNCTION ICL2AMP(NEWLBLA,ICAUA,ISYMA,
     *                         NEWLBLB,ICAUB,ISYMB )
*---------------------------------------------------------------------*
C
C maintain the list of second-order left Cauchy vectors:
C
C   if vector is on the list return list index and set ISYMA,ISYMB
C   if vector is NOT on the list:
C        LCL2OPN=.true.  --> extend list, and return index
C        LCL2OPN=.false. --> return -1   
C
C Christof Haettig, March 98
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <cccl2rsp.h>
#include <priunit.h>
C
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
C
      CHARACTER*8 NEWLBLA, NEWLBLB
      INTEGER ICAUA, ICAUB, I
      INTEGER ISYMA, ISYMB

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) '[DEBUG] ICL2AMP> entered with input:'
        WRITE (LUPRI,*) '[DEBUG] ICL2AMP> LABELS    :',NEWLBLA,NEWLBLB
        WRITE (LUPRI,*) '[DEBUG] ICL2AMP> SYMMETRIES:',ISYMA,ISYMB
        WRITE (LUPRI,*) '[DEBUG] ICL2AMP> CAUCHY ORD:',ICAUA,ICAUB
      END IF

      DO I = 1,NCL2lBL
         IF (     NEWLBLA.EQ.LBLCL2(I,1) .AND. NEWLBLB.EQ.LBLCL2(I,2)
     *       .AND. ICAUA.EQ.ICL2CAU(I,1) .AND. ICAUB.EQ.ICL2CAU(I,2)
     *      ) THEN
            ICL2AMP  = I
            ISYMA    = ISYCL2(ICL2AMP,1)
            ISYMB    = ISYCL2(ICL2AMP,2)
            IF (LOCDBG) THEN
              WRITE (LUPRI,*)
     &              '[DEBUG] ICL2AMP> entry found on the list:'
              WRITE (LUPRI,*) '[DEBUG] ICL2AMP> INDEX    :',ICL2AMP
              WRITE (LUPRI,*) '[DEBUG] ICL2AMP> SYMMETRIES:',ISYMA,ISYMB
            END IF
            RETURN
         END IF
         IF (     NEWLBLA.EQ.LBLCL2(I,2) .AND. NEWLBLB.EQ.LBLCL2(I,1)
     *       .AND. ICAUA.EQ.ICL2CAU(I,2) .AND. ICAUB.EQ.ICL2CAU(I,1)
     *      ) THEN
            ICL2AMP  = I
            ISYMB    = ISYCL2(ICL2AMP,1)
            ISYMA    = ISYCL2(ICL2AMP,2)
            IF (LOCDBG) THEN
              WRITE (LUPRI,*)
     &              '[DEBUG] ICL2AMP> entry found on the list:'
              WRITE (LUPRI,*) '[DEBUG] ICL2AMP> INDEX    :',ICL2AMP
              WRITE (LUPRI,*) '[DEBUG] ICL2AMP> SYMMETRIES:',ISYMA,ISYMB
            END IF
            RETURN
         END IF
      END DO  

      IF (LCL2OPN) THEN
        NCL2lBL = NCL2lBL + 1

        IF (NCL2lBL.GT.MAXCL2LBL) THEN
          WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *    '@ NUMBER OF SPECIFIED VECTORS EXCEED THE MAXIMUM ALLOWED',
     *    '@ MAXCL2LBL =',MAXCL2LBL,' NCL2lBL= ',NCL2lBL
          CALL QUIT(' ICL2AMP: TOO MANY VECTORS SPECIFIED')
        END IF

        LBLCL2(NCL2lBL,1)  = NEWLBLA
        LBLCL2(NCL2lBL,2)  = NEWLBLB
        ICL2CAU(NCL2lBL,1) = ICAUA
        ICL2CAU(NCL2lBL,2) = ICAUB
        ISYCL2(NCL2lBL,1)  = ISYMA
        ISYCL2(NCL2lBL,2)  = ISYMB
        ICL2AMP = NCL2lBL

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) '[DEBUG] ICL2AMP> put entry on the list:'
          WRITE (LUPRI,*) '[DEBUG] ICL2AMP> INDEX    :',ICL2AMP
        END IF
      ELSE
        WRITE(LUPRI,'(3A,I2,3A,I2,2A)') 
     *   '@ WARNING: CL2 VECTOR FOR ',NEWLBLA,'(',ICAUA,'), ',
     *                                NEWLBLB,'(',ICAUB,')',
     *              ' IS NOT AVAILABLE.'
        ICL2AMP = -1
      END IF

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck icr2amp */
      INTEGER FUNCTION ICR2AMP(NEWLBLA,ICAUA,ISYMA,
     *                         NEWLBLB,ICAUB,ISYMB )
*---------------------------------------------------------------------*
C
C maintain the list of second-order right Cauchy vectors:
C
C   if vector is on the list return list index and set ISYMA,ISYMB
C   if vector is NOT on the list:
C        LCR2OPN=.true.  --> extend list, and return index
C        LCR2OPN=.false. --> return -1   
C
C Christof Haettig, March 98
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <cccr2rsp.h>
#include <priunit.h>
C
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
C
      CHARACTER*8 NEWLBLA, NEWLBLB
      INTEGER ICAUA, ICAUB, I
      INTEGER ISYMA, ISYMB

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) '[DEBUG] ICR2AMP> entered with input:'
        WRITE (LUPRI,*) '[DEBUG] ICR2AMP> LABELS    :',NEWLBLA,NEWLBLB
        WRITE (LUPRI,*) '[DEBUG] ICR2AMP> SYMMETRIES:',ISYMA,ISYMB
        WRITE (LUPRI,*) '[DEBUG] ICR2AMP> CAUCHY ORD:',ICAUA,ICAUB
      END IF

      DO I = 1,NCR2LBL
         IF (     NEWLBLA.EQ.LBLCR2(I,1) .AND. NEWLBLB.EQ.LBLCR2(I,2)
     *       .AND. ICAUA.EQ.ICR2CAU(I,1) .AND. ICAUB.EQ.ICR2CAU(I,2)
     *      ) THEN
            ICR2AMP  = I
            ISYMA    = ISYCR2(ICR2AMP,1)
            ISYMB    = ISYCR2(ICR2AMP,2)
            IF (LOCDBG) THEN
              WRITE (LUPRI,*)
     &              '[DEBUG] ICR2AMP> entry found on the list:'
              WRITE (LUPRI,*) '[DEBUG] ICR2AMP> INDEX    :',ICR2AMP
              WRITE (LUPRI,*) '[DEBUG] ICR2AMP> SYMMETRIES:',ISYMA,ISYMB
            END IF
            RETURN
         END IF
         IF (     NEWLBLA.EQ.LBLCR2(I,2) .AND. NEWLBLB.EQ.LBLCR2(I,1)
     *       .AND. ICAUA.EQ.ICR2CAU(I,2) .AND. ICAUB.EQ.ICR2CAU(I,1)
     *      ) THEN
            ICR2AMP  = I
            ISYMB    = ISYCR2(ICR2AMP,1)
            ISYMA    = ISYCR2(ICR2AMP,2)
            IF (LOCDBG) THEN
              WRITE (LUPRI,*)
     &              '[DEBUG] ICR2AMP> entry found on the list:'
              WRITE (LUPRI,*) '[DEBUG] ICR2AMP> INDEX    :',ICR2AMP
              WRITE (LUPRI,*) '[DEBUG] ICR2AMP> SYMMETRIES:',ISYMA,ISYMB
            END IF
            RETURN
         END IF
      END DO  

      IF (LCR2OPN) THEN
        NCR2LBL = NCR2LBL + 1

        IF (NCR2LBL.GT.MAXCR2LBL) THEN
          WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *    '@ NUMBER OF SPECIFIED VECTORS EXCEED THE MAXIMUM ALLOWED',
     *    '@ MAXCR2LBL =',MAXCR2LBL,' NCR2LBL= ',NCR2LBL
          CALL QUIT(' ICR2AMP: TOO MANY VECTORS SPECIFIED')
        END IF

        LBLCR2(NCR2LBL,1)  = NEWLBLA
        LBLCR2(NCR2LBL,2)  = NEWLBLB
        ICR2CAU(NCR2LBL,1) = ICAUA
        ICR2CAU(NCR2LBL,2) = ICAUB
        ISYCR2(NCR2LBL,1)  = ISYMA
        ISYCR2(NCR2LBL,2)  = ISYMB
        ICR2AMP = NCR2LBL

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) '[DEBUG] ICR2AMP> put entry on the list:'
          WRITE (LUPRI,*) '[DEBUG] ICR2AMP> INDEX    :',ICR2AMP
        END IF
      ELSE
        WRITE(LUPRI,'(3A,I2,3A,I2,2A)') 
     *   '@ WARNING: CR2 VECTOR FOR ',NEWLBLA,'(',ICAUA,'), ',
     *                                NEWLBLB,'(',ICAUB,')',
     *              ' IS NOT AVAILABLE.'
        ICR2AMP = -1
      END IF

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck IETACL2 */
      INTEGER FUNCTION IETACL2(NEWLBLA,ICAUA,ISYMA,
     *                         NEWLBLB,ICAUB,ISYMB )
*---------------------------------------------------------------------*
C
C maintain the list of second-order right Cauchy vectors:
C
C   if vector is on the list return list index and set ISYMA,ISYMB
C   if vector is NOT on the list:
C        LCX2OPN=.true.  --> extend list, and return index
C        LCX2OPN=.false. --> return -1   
C
C Christof Haettig, March 98
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <cccx2rsp.h>
#include <priunit.h>
C
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
C
      CHARACTER*8 NEWLBLA, NEWLBLB
      INTEGER ICAUA, ICAUB, I
      INTEGER ISYMA, ISYMB

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) '[DEBUG] IETACL2> entered with input:'
        WRITE (LUPRI,*) '[DEBUG] IETACL2> LABELS    :',NEWLBLA,NEWLBLB
        WRITE (LUPRI,*) '[DEBUG] IETACL2> SYMMETRIES:',ISYMA,ISYMB
        WRITE (LUPRI,*) '[DEBUG] IETACL2> CAUCHY ORD:',ICAUA,ICAUB
      END IF

      DO I = 1,NCX2LBL
         IF (     NEWLBLA.EQ.LBLCX2(I,1) .AND. NEWLBLB.EQ.LBLCX2(I,2)
     *       .AND. ICAUA.EQ.ICX2CAU(I,1) .AND. ICAUB.EQ.ICX2CAU(I,2)
     *      ) THEN
            IETACL2  = I
            ISYMA    = ISYCX2(IETACL2,1)
            ISYMB    = ISYCX2(IETACL2,2)
            IF (LOCDBG) THEN
              WRITE (LUPRI,*)
     &              '[DEBUG] IETACL2> entry found on the list:'
              WRITE (LUPRI,*) '[DEBUG] IETACL2> INDEX    :',IETACL2
              WRITE (LUPRI,*) '[DEBUG] IETACL2> SYMMETRIES:',ISYMA,ISYMB
            END IF
            RETURN
         END IF
         IF (     NEWLBLA.EQ.LBLCX2(I,2) .AND. NEWLBLB.EQ.LBLCX2(I,1)
     *       .AND. ICAUA.EQ.ICX2CAU(I,2) .AND. ICAUB.EQ.ICX2CAU(I,1)
     *      ) THEN
            IETACL2  = I
            ISYMB    = ISYCX2(IETACL2,1)
            ISYMA    = ISYCX2(IETACL2,2)
            IF (LOCDBG) THEN
              WRITE (LUPRI,*)
     &              '[DEBUG] IETACL2> entry found on the list:'
              WRITE (LUPRI,*) '[DEBUG] IETACL2> INDEX    :',IETACL2
              WRITE (LUPRI,*) '[DEBUG] IETACL2> SYMMETRIES:',ISYMA,ISYMB
            END IF
            RETURN
         END IF
      END DO  

      IF (LCX2OPN) THEN
        NCX2LBL = NCX2LBL + 1

        IF (NCX2LBL.GT.MAXCX2LBL) THEN
          WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *    '@ NUMBER OF SPECIFIED VECTORS EXCEED THE MAXIMUM ALLOWED',
     *    '@ MAXCX2LBL =',MAXCX2LBL,' NCX2LBL= ',NCX2LBL
          CALL QUIT(' IETACL2: TOO MANY VECTORS SPECIFIED')
        END IF

        LBLCX2(NCX2LBL,1)  = NEWLBLA
        LBLCX2(NCX2LBL,2)  = NEWLBLB
        ICX2CAU(NCX2LBL,1) = ICAUA
        ICX2CAU(NCX2LBL,2) = ICAUB
        ISYCX2(NCX2LBL,1)  = ISYMA
        ISYCX2(NCX2LBL,2)  = ISYMB
        IETACL2 = NCX2LBL

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) '[DEBUG] IETACL2> put entry on the list:'
          WRITE (LUPRI,*) '[DEBUG] IETACL2> INDEX    :',IETACL2
        END IF
      ELSE
        WRITE(LUPRI,'(3A,I2,3A,I2,2A)') 
     *   '@ WARNING: CX2 VECTOR FOR ',NEWLBLA,'(',ICAUA,'), ',
     *                                NEWLBLB,'(',ICAUB,')',
     *              ' IS NOT AVAILABLE.'
        IETACL2 = -1
      END IF

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck IRHSCR2 */
      INTEGER FUNCTION IRHSCR2(NEWLBLA,ICAUA,ISYMA,
     *                         NEWLBLB,ICAUB,ISYMB )
*---------------------------------------------------------------------*
C
C maintain the list of rhs vectors for second-order 
C right Cauchy vector equations:
C
C   if vector is on the list return list index and set ISYMA,ISYMB
C   if vector is NOT on the list:
C        LCO2OPN=.true.  --> extend list, and return index
C        LCO2OPN=.false. --> return -1   
C
C Christof Haettig, March 98
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <ccco2rsp.h>
#include <priunit.h>
C
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
C
      CHARACTER*8 NEWLBLA, NEWLBLB
      INTEGER ICAUA, ICAUB, I
      INTEGER ISYMA, ISYMB

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) '[DEBUG] IRHSCR2> entered with input:'
        WRITE (LUPRI,*) '[DEBUG] IRHSCR2> LABELS    :',NEWLBLA,NEWLBLB
        WRITE (LUPRI,*) '[DEBUG] IRHSCR2> SYMMETRIES:',ISYMA,ISYMB
        WRITE (LUPRI,*) '[DEBUG] IRHSCR2> CAUCHY ORD:',ICAUA,ICAUB
      END IF

      DO I = 1,NCO2LBL
         IF (     NEWLBLA.EQ.LBLCO2(I,1) .AND. NEWLBLB.EQ.LBLCO2(I,2)
     *       .AND. ICAUA.EQ.ICO2CAU(I,1) .AND. ICAUB.EQ.ICO2CAU(I,2)
     *      ) THEN
            IRHSCR2  = I
            ISYMA    = ISYCO2(IRHSCR2,1)
            ISYMB    = ISYCO2(IRHSCR2,2)
            IF (LOCDBG) THEN
              WRITE (LUPRI,*)
     &              '[DEBUG] IRHSCR2> entry found on the list:'
              WRITE (LUPRI,*) '[DEBUG] IRHSCR2> INDEX    :',IRHSCR2
              WRITE (LUPRI,*) '[DEBUG] IRHSCR2> SYMMETRIES:',ISYMA,ISYMB
            END IF
            RETURN
         END IF
         IF (     NEWLBLA.EQ.LBLCO2(I,2) .AND. NEWLBLB.EQ.LBLCO2(I,1)
     *       .AND. ICAUA.EQ.ICO2CAU(I,2) .AND. ICAUB.EQ.ICO2CAU(I,1)
     *      ) THEN
            IRHSCR2  = I
            ISYMB    = ISYCO2(IRHSCR2,1)
            ISYMA    = ISYCO2(IRHSCR2,2)
            IF (LOCDBG) THEN
              WRITE (LUPRI,*)
     &              '[DEBUG] IRHSCR2> entry found on the list:'
              WRITE (LUPRI,*) '[DEBUG] IRHSCR2> INDEX    :',IRHSCR2
              WRITE (LUPRI,*) '[DEBUG] IRHSCR2> SYMMETRIES:',ISYMA,ISYMB
            END IF
            RETURN
         END IF
      END DO  

      IF (LCO2OPN) THEN
        NCO2LBL = NCO2LBL + 1

        IF (NCO2LBL.GT.MAXCO2LBL) THEN
          WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *    '@ NUMBER OF SPECIFIED VECTORS EXCEED THE MAXIMUM ALLOWED',
     *    '@ MAXCO2LBL =',MAXCO2LBL,' NCO2LBL= ',NCO2LBL
          CALL QUIT(' IRHSCR2: TOO MANY VECTORS SPECIFIED')
        END IF

        LBLCO2(NCO2LBL,1)  = NEWLBLA
        LBLCO2(NCO2LBL,2)  = NEWLBLB
        ICO2CAU(NCO2LBL,1) = ICAUA
        ICO2CAU(NCO2LBL,2) = ICAUB
        ISYCO2(NCO2LBL,1)  = ISYMA
        ISYCO2(NCO2LBL,2)  = ISYMB
        IRHSCR2 = NCO2LBL

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) '[DEBUG] IRHSCR2> put entry on the list:'
          WRITE (LUPRI,*) '[DEBUG] IRHSCR2> INDEX    :',IRHSCR2
        END IF
      ELSE
        WRITE(LUPRI,'(3A,I2,3A,I2,2A)') 
     *   '@ WARNING: CO2 VECTOR FOR ',NEWLBLA,'(',ICAUA,'), ',
     *                                NEWLBLB,'(',ICAUB,')',
     *              ' IS NOT AVAILABLE.'
        IRHSCR2 = -1
      END IF

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck irshr1 */
      INTEGER FUNCTION IRHSR1(NEWLBL,LORX,FRQINP,ISYM)
*---------------------------------------------------------------------*
C
C maintain the list of right hand side vectors for
C first-order t amplitude responses
C
C   if vector is on the list return list index and set ISYM
C   if vector is NOT on the list:
C        LO1OPN=.true.  --> extend list, and return index
C        LO1OPN=.false. --> return -1   
C
C        NEWLBL -- operator label
C        LORX   -- flag for orbital relaxation
C        FRQINP -- frequency (ignored for unrelaxed orbitals)
C        ISYM   -- symmetry
C
C Christof Haettig, Juni 1998
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <cco1rsp.h>
#include <priunit.h>
#include <ccsdinp.h>
 
      LOGICAL LORX
      INTEGER ISYM

      REAL*8  FRQNEW,TOL, FRQINP

      PARAMETER(TOL=1.0D-12)

      CHARACTER*8 NEWLBL
      INTEGER I

* if LORX false ignore frequency (set to zero internally):
      IF (LORX.OR.CCSDT) THEN
        FRQNEW = FRQINP
      ELSE
        FRQNEW = 0.0d0
      END IF

      DO I = 1,NO1LBL
         IF ( (NEWLBL .EQ. LBLO1(I)) .AND. (LORX .EQV. LORXO1(I)) .AND.
     &        (ABS(FRQNEW-FRQO1(I)).LT.TOL)) THEN
            IRHSR1 = I
            ISYM   = ISYO1(IRHSR1)
            RETURN
         END IF
      END DO  

      IF (LO1OPN) THEN
        NO1LBL = NO1LBL + 1

        IF (NO1LBL.GT.MAXO1LBL) THEN
          WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *    '@ NUMBER OF SPECIFIED EQUATIONS EXCEED THE MAXIMUM ALLOWED',
     *    '@ MAXO1LBL =',MAXO1LBL,' NO1LBL= ',NO1LBL
          CALL QUIT(' IRHSR1: TOO MANY EQUATIONS SPECIFIED')
        END IF

        LBLO1(NO1LBL)  = NEWLBL
        ISYO1(NO1LBL)  = ISYM
        LORXO1(NO1LBL) = LORX
        FRQO1(NO1LBL)  = FRQNEW
        IRHSR1         = NO1LBL

      ELSE
        WRITE(LUPRI,'(3A,L2,A,1P,D12.5,2A)')
     *   '@ WARNING: RHSR1 VECTOR FOR ',NEWLBL,'(',LORX,',',FRQNEW,')',
     *                  ' IS NOT AVAILABLE.'
        IRHSR1 = -1
      END IF

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck ieta1 */
      INTEGER FUNCTION IETA1(NEWLBL,LORX,FRQINP,ISYM)
*---------------------------------------------------------------------*
C
C maintain the list of right hand side vectors for
C first-order lagrangian multiplier responses
C
C   if vector is on the list return list index and set ISYM
C   if vector is NOT on the list:
C        LX1OPN=.true.  --> extend list, and return index
C        LX1OPN=.false. --> return -1   
C
C        NEWLBL -- operator label
C        LORX   -- flag for orbital relaxation
C        FRQINP -- frequency (ignored for unrelaxed orbitals)
C        ISYM   -- symmetry 
C
C Christof Haettig, Juni 1998
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <ccx1rsp.h>
#include <priunit.h>
#include <ccsdinp.h>
 
      LOGICAL LORX, LORX1
      INTEGER ISYM

      REAL*8  FRQNEW,TOL, FRQINP

      PARAMETER(TOL=1.0D-12)

      CHARACTER*8 NEWLBL
      INTEGER I

* if LORX false and CCSDT false ignore frequency by setting it
* to zero internally:
      IF (LORX.OR.CCSDT) THEN
        FRQNEW = FRQINP
      ELSE
        FRQNEW = 0.0d0
      END IF


      DO I = 1,NX1LBL
         ! the following crappy comparison of LORX with LORXX1 was
         ! necessary to get it through the XLF compilers
         LORX1 = (LORX.AND.LORXX1(I)) .OR.
     &           ((.NOT.LORX).AND.(.NOT.LORXX1(I)))
         IF ( (NEWLBL .EQ. LBLX1(I)) .AND. LORX1 .AND. 
     &        (ABS(FRQNEW-FRQX1(I)).LT.TOL)) THEN
            IETA1 = I
            ISYM  = ISYX1(IETA1)
            RETURN
         END IF
      END DO  

      IF (LX1OPN) THEN
        NX1LBL = NX1LBL + 1

        IF (NX1LBL.GT.MAXX1LBL) THEN
          WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *    '@ NUMBER OF SPECIFIED EQUATIONS EXCEED THE MAXIMUM ALLOWED',
     *    '@ MAXX1LBL =',MAXX1LBL,' NX1LBL= ',NX1LBL
          CALL QUIT(' IETA1: TOO MANY EQUATIONS SPECIFIED')
        END IF

        LBLX1(NX1LBL)  = NEWLBL
        ISYX1(NX1LBL)  = ISYM
        LORXX1(NX1LBL) = LORX
        FRQX1(NX1LBL)  = FRQNEW
        IETA1          = NX1LBL

      ELSE
        WRITE(LUPRI,'(3A,L2,A,1P,D12.5,2A)') 
     *   '@ WARNING: ETA1 VECTOR FOR ',NEWLBL,'(',LORX,',',FRQNEW,')',
     *                  ' IS NOT AVAILABLE.'
        IETA1 = -1
      END IF

      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck ipl1zeta */
      INTEGER FUNCTION IPL1ZETA(NEWLBLA,LORXA,FRQANEW,ISYMA,LPROJ,
     &                          IEXCI,EIGVNEW,ISYMS)
*---------------------------------------------------------------------*
C
C maintain the list of projected first order zeta amplitude responses
C onto the orthogonal complement of E^f
C
C   if vector is on the list return list index and set ISYMS,ISYMA
C
C   if vector is NOT on the list:
C        LPL1OPN=.true.  --> extend list, and return index IPL1ZETA
C        LPL1OPN=.false. --> return -1
C
C        NEWLBLA -- operator A label
C        LORXA   -- flag for orbital relaxation
C        FRQANEW -- frequency
C        ISYMA   -- symmetry of operator A and of projected PL1
C        LPROJ   -- flag for projection
C
C        IEXCI   -- index for the excited state
C        EIGVNEW -- its eigenvalue (exc. energy)
C        ISYMS   -- its symmetry
C
C Sonia Coriani, March 2000
C based of IL1ZETA and IEL1AMP
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <ccpl1rsp.h>
#include <priunit.h>

      LOGICAL LORXA, LPROJ, LPROJ1
      INTEGER ISYMA, ISYMS, IEXCI, I
      REAL*8  FRQANEW, EIGVNEW, TOL

      PARAMETER(TOL=1.0D-12)

      CHARACTER*8 NEWLBLA


      LPROJ1 = LPROJ
*
* Projection might only be necessary if ISYMA = ISYMS
*
      IF (ISYMA.NE.ISYMS) THEN
         LPROJ1 = .FALSE.
*         WRITE(LUPRI,*) ' Inside IPL1ZETA: LPROJ1 reset to FALSE'
      END IF
*
      DO I = 1,NPL1LBL
         IF ( (NEWLBLA .EQ. LBLPL1(I))         .AND.
     &        (LORXA .EQV. LORXPL1(I))         .AND.
     &        (ABS(FRQANEW-FRQPL1(I)).LT.TOL)  .AND.
     &        (IEXCI .EQ.  ISTPL1(I))          .AND.
     &        (ABS(EIGVNEW-EIGPL1(I)).LT.TOL)  .AND.
     &        (LPROJ1 .EQV. LPRPL1(I))        ) THEN

            IPL1ZETA = I
            ISYMA    = ISYPL1(IPL1ZETA)
            ISYMS    = ISYSPL1(IPL1ZETA)
            RETURN
         END IF
      END DO

      IF (LPL1OPN) THEN
        NPL1LBL = NPL1LBL + 1

        IF (NPL1LBL.GT.MAXPL1LBL) THEN
          WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *    '@ NUMBER OF SPECIFIED EQUATIONS EXCEED THE MAXIMUM ALLOWED',
     *    '@ MAXPL1LBL =',MAXPL1LBL,' NPL1LBL= ',NPL1LBL
          CALL QUIT(' IPL1ZETA: TOO MANY EQUATIONS SPECIFIED')
        END IF

        !the A-operator-for-response infos
        LBLPL1(NPL1LBL)  = NEWLBLA
        LORXPL1(NPL1LBL) = LORXA
        FRQPL1(NPL1LBL)  = FRQANEW
        ISYPL1(NPL1LBL)  = ISYMA
        !the excitated-state-for-projection infos
        ISTPL1(NPL1LBL)  = IEXCI
        ISYSPL1(NPL1LBL) = ISYMS
        EIGPL1(NPL1LBL)  = EIGVNEW
        !the PL^A-vector extra infos
        LPRPL1(NPL1LBL)  = LPROJ1
        IPL1ZETA         = NPL1LBL

      ELSE
        WRITE(LUPRI,'(3A,L2,A,1P,D12.5,2A)')
     *   '@ WARNING: PL1 VECTOR FOR ',NEWLBLA,'(',LORXA,',',FRQANEW,')',
     *                  ' IS NOT AVAILABLE.'
        IPL1ZETA = -1
      END IF

      RETURN
      END
*=====================================================================*
