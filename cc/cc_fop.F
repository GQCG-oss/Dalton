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
C  /* Deck cc_fop */
       SUBROUTINE CC_FOP(WORK,LWORK)
C
C-----------------------------------------------------------------------------
C
C     Purpose: Direct calculation of Coupled Cluster
C              first order properties
C
C              CCS(CIS/HF), MP2, CCD, CCSD, CC3, MCC2
C
C              CCSDT-1a, CCSDT-1b
C
C              and calculates modified triples corrections MCCSD(T), MCC(3)
C
C     Solves for CC t-bar amplitudes = Lagrangian multipliers.
C     For relaxed properties also for orbital multipliers.
C     Calculates various first order one-electron properties.
C
C     Initiated by Ove Christiansen 15 November 1994.
C     CCSD  one electron FOP by Asger Halkier April 1996.
C     MP2   one electron FOP by Asger Halkier September 1996.
C     New CC solvers introduced, Ove Christiansen November 1996.
C     Frozen core contribution to unrelaxed density Ove Christiansen May 1996.
C     MCC2, Ove Christiansen, july 1997
C     Major clean-up of overall structure by Asger Halkier March 1998.
C     New MP2 & CCSD version based on canonical orbitals throughout the whole
C     surface by Asger Halkier Spring 1998. This includes frozen core for
C     the relaxed density.
C
C     Relaxed CC2 FOP by A. Halkier & S. Coriani January 2000.
C     No frozen core possible for Relaxed CC2 initially.
C
C-----------------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
#include <alphac.h>
#include <iratdef.h>
#include <ccfop.h>
#include <cclr.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccrspprp.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <ccroper.h>
#include <ccfield.h>
#include <exeinf.h>
#include <infvar.h>
#include <inftap.h>
#include <dipole.h>
#include <quadru.h>
#include <nqcc.h>
#include <ccfdgeo.h>
#include <ccfro.h>
#include <ccinftap.h>
#include <ccslvinf.h>
C
      LOGICAL CC1BSV,CC1ASV,FASAVE,CICLC,HFCLC, DAR2SA,
     *        TRPCLC,OOTV, EXCLC, RLORBS, LPROJECT, EX, TRIPLET,
     *        LDUM
      DIMENSION WORK(LWORK), ELSEMO(3,3), SKODE(3,3), SKODN(3,3)
      CHARACTER*1 LR, CDUM
      CHARACTER MODEL*10,ETY1*5
      CHARACTER MODELPRI*4, MODELPRI2*14
      CHARACTER*8 LABEL
      CHARACTER*3 LIST
C
      INTEGER ISYOF(8)
C
#include <leinf.h>
C
      CALL QENTER('CC_FOP')
C
C------------------------------------
C     Header of Property calculation.
C------------------------------------
C
      WRITE (LUPRI,'(1X,A,/)') '  '
      WRITE (LUPRI,'(1X,A)')
     *'*********************************************************'//
     *'**********'
      WRITE (LUPRI,'(1X,A)')
     *'*                                                        '//
     *'         *'
      WRITE (LUPRI,'(1X,A)')
     *'*<<<<<<<<<<<<< OUTPUT FROM COUPLED CLUSTER RESPONSE  >>>>'//
     *'>>>>>>>>>*'
      IF ( CCFOP  ) THEN
         WRITE (LUPRI,'(1X,A)')
     *   '*                                                        '//
     *   '         *'
         WRITE (LUPRI,'(1X,A)')
     *   '*<<<<<<<<<<    CALCULATION OF FIRST ORDER PROPERTIES    >'//
     *   '>>>>>>>>>*'
      ENDIF
      WRITE (LUPRI,'(1X,A)')
     *'*                                                        '//
     *'         *'
      WRITE (LUPRI,'(1X,A,/)')
     *'*********************************************************'//
     *'**********'
C
      MODEL = 'CCSD'
      IF (CC2.AND.(.NOT.MCC2)) THEN
         CALL AROUND('Coupled Cluster model is: CC2')
         MODEL = 'CC2'
         MODELPRI = ' CC2'
      ENDIF
      IF (MP2) THEN
         CALL AROUND('Model is second order pert. theory: MP2 ')
         MODEL = 'MP2'
         MODELPRI = ' MP2'
      ENDIF
      IF (CCS.AND.(.NOT.CIS)) THEN
         CALL AROUND('Coupled Cluster model is: CCS')
         MODEL = 'CCS'
         MODELPRI = ' CCS'
      ENDIF
      IF (CCS.AND.CIS) THEN
         CALL AROUND('CIS model in use ')
         MODEL = 'CCS'
         MODELPRI = ' CIS'
      ENDIF
      IF (CCD) THEN
         CALL AROUND('Coupled Cluster model is: CCD')
         MODEL = 'CCD'
         MODELPRI = ' CCD'
      ENDIF
      IF (CC3  ) THEN
         CALL AROUND('Coupled Cluster model is: CC3')
         MODEL = 'CC3'
         MODELPRI = ' CC3'
         CALL QUIT('CC3 first order properties not implemented')
      ENDIF
      IF (CC1A) THEN
         CALL AROUND('Coupled Cluster model is: CCSDT-1a')
         MODEL = 'CCSDT-1a'
         CALL QUIT('CCSDT-1a first order properties not implemented')
      ENDIF
      IF (CC1B) THEN
         CALL AROUND('Coupled Cluster model is: CCSDT-1b')
         MODEL = 'CCSDT-1b'
         CALL QUIT('CCSDT-1b first order properties not implemented')
      ENDIF
      IF (CCSD) THEN
         CALL AROUND('Coupled Cluster model is: CCSD')
         MODEL = 'CCSD'
         MODELPRI = 'CCSD'
      ENDIF
C
      IF (MP2) RLORBS = .TRUE.
      IF (CC2) RLORBS = .FALSE.
      IF (RELORB .AND. CC2) THEN
         IF ((FROIMP) .OR. (FROEXP)) THEN
            WRITE(LUPRI,*) 
     *         'No frozen core for relaxed cc2 implemented yet'
            CALL QUIT('NO FROZEN CORE FOR RELAXED CC2 YET')
         ENDIF
      ENDIF
C
      IF ((.NOT.RELORB) .AND. MP2) THEN
         WRITE(LUPRI,*) 'WARNING: MP2 unrelaxed first order properties '
     *              //'not implemented '
         WRITE(LUPRI,*) 'Orbital relaxation switched on for MP2.'
         RLORBS = RELORB
         RELORB = .TRUE. 
      ENDIF
C
      IF (IPRINT.GT.10) WRITE(LUPRI,*) 'CC_FOP-1: Workspace:',LWORK
C
C-----------------------------
C     Initialize Variables.
C-----------------------------
C
      ISYMTR = ISYMOP
      LIST  = 'L0 '
C
C----------------------------------------------------------------
C     In case of CCS calculation, no equations need to be solved,
C     and we jump directly to calculating the requested first
C     order properties, which are identical to the HF-results.
C     In case of MP2 calculation, we need not solve equations to
C     obtaine the amplitude multipliers, which are evaluated
C     straightforwardly from integrals L(iajb).
C     Otherwise we must set up the right hand side and solve the
C     equations.
C----------------------------------------------------------------
C
      IF (CCS) GOTO 47
      IF (L0SKIP) GOTO 46
C
      IF (MP2) THEN
C
         KMP2LA = 1
         KWRK1  = KMP2LA + NT1AMX + NT2AMX
         LWRK1  = LWORK  - KWRK1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Needed:', KWRK1, 'Available:', LWORK
            CALL QUIT('Insufficient memory for initial allocation in '//
     &                'cc_fop')
         ENDIF
C
         CALL DZERO(WORK(KMP2LA),NT1AMX + NT2AMX)
C
         CALL MP_LAM(WORK(KMP2LA),WORK(KWRK1),LWRK1)
C
         KWRK2 = KMP2LA
         LWRK2 = LWRK1
C
C
         IF ( IPRINT .GT. 10 .OR. DEBUG) THEN
            RHO1N = DDOT(NT1AM(ISYMTR),WORK(KWRK2),1,WORK(KWRK2),1)
            RHO2N = DDOT(NT2AM(ISYMTR),WORK(KWRK2+NT1AMX),1,
     *                   WORK(KWRK2+NT1AMX),1)
            WRITE(LUPRI,*) 'Norm of singles Lambda vector :',RHO1N
            WRITE(LUPRI,*) 'Norm of doubles Lambda vector :',RHO2N
         ENDIF
C
         IF ( IPRINT .GT. 30 ) THEN
            CALL AROUND('CCLR_FOP: Lambda vector in mo basis' )
            CALL OUTPUT(WORK(KWRK2),1,NT1AMX+NT2AMX,1,1,
     *                                NT1AMX+NT2AMX,1,1,LUPRI)
         ENDIF
C
         IF (IPRINT.GT.1) THEN
           DDUMMY  = 0.0D0
           WRITE(LUPRI,'(//1X,A)')
     *       'Analysis of the undifferentiated Lagrangian multipliers:'
           WRITE(LUPRI,'(1X,A)')
     *       '--------------------------------------------------------'
           CALL CC_PRAM(WORK(KWRK2),DDUMMY,ISYMTR)
         END IF
C
         KWRK3  = KWRK2 + NT1AMX + NT2AMX
         LWRK3  = LWORK - KWRK3
C
         IOPT   = 3
         CALL CC_WRRSP('L0',0,1,IOPT,MODEL,DUMMY,
     *                 WORK(KWRK2),WORK(KWRK2+NT1AM(ISYMTR)),
     *                 WORK(KWRK3),LWRK3)

         IF ( IPRINT .GT. 10 .OR. DEBUG) THEN
            RHO1N = DDOT(NT1AM(ISYMTR),WORK(KWRK2),1,WORK(KWRK2),1)
            RHO2N = DDOT(NT2AM(ISYMTR),WORK(KWRK2+NT1AMX),1,
     *                   WORK(KWRK2+NT1AMX),1)
            WRITE(LUPRI,*) 'Norm of singles Lambda vector :',RHO1N
            WRITE(LUPRI,*) 'Norm of doubles Lambda vector :',RHO2N
         ENDIF
      ELSE

         NSTAT = 0
         ORDER = 0
         ISIDE = -1

         ISYOF(1) = 0
         DO I = 2, NSYM
           ISYOF(I) = 1
         END DO

         CALL CC_SOLDRV(LIST,NSTAT,ORDER,ISIDE,
     *                  IDUM,IDUM,RDUM,LDUM,
     *                  IDUM,CDUM,RDUM,IDUM,
     *                  ISYOF,1,1,WORK,LWORK)

      ENDIF
C
C---------------------------------------------------
C SLV98,OC Solvent part 1
C           Calculate norm and test for convergence.
C---------------------------------------------------
C
      IF (CCSLV) THEN
C
              KLAM   = 1
              KLAM2  = 1 + NT1AMX
C
              IOPT = 3
              CALL CC_RDRSP('L0',0,1,IOPT,MODEL,WORK(KLAM),WORK(KLAM2))

              XLNCCCU = DDOT(NT1AMX+NT2AMX,WORK(KLAM),1,WORK(KLAM),1)
              IF (ABS(XLNCCPR-XLNCCCU).LT.CVGLSOL) LSLLCVG = .TRUE.
              IF (IPRINT.GT.2) THEN
                WRITE(LUPRI,*)
     *          'Norm of L-amplitudes in this solvent it.:',XLNCCCU
                WRITE(LUPRI,*)
     *          'Norm of L-amplitudes in prev solvent it.:',XLNCCPR
                WRITE(LUPRI,*) 'LSLLCVG: ',LSLLCVG
              ENDIF
              WRITE(LUPRI,*)
     *        ' Change in norm^2 of L-amplitudes in this solvent it.:',
     *        XLNCCCU-XLNCCPR

              XLNCCPR = XLNCCCU
C
              KWRK3  = KLAM  + NT1AMX + NT2AMX
              KRHO1  = KWRK3
              KRHO2  = KRHO1 + NT1AMX
              KWRK4  = KRHO2 + NT2AMX
              LWRK4  = LWORK - KWRK4
              IF (LWRK4.LE.0) CALL QUIT(' Too little work in cc_fop')
              CALL DZERO(WORK(KRHO1),NT1AMX)
              IF (.NOT.CCS) CALL DZERO(WORK(KRHO2),NT2AMX)
              LR = '0'
              CALL CCSL_LTRB(WORK(KRHO1),WORK(KRHO2),DUM1,DUM2,
     *                       ISYMOP,LR,WORK(KWRK4),LWRK4)
              KOMEG1 = KWRK4
              KOMEG2 = KWRK4 + NT1AMX
              LUOME = -9000
              CALL GPOPEN(LUOME,'CC_OME','UNKNOWN',' ',
     *            'UNFORMATTED',IDUMMY,.FALSE.)
              REWIND (LUOME)
              READ(LUOME) (WORK(KOMEG1+K-1), K = 1,NT1AMX)
              IF (.NOT.CCS) THEN
                 READ(LUOME) (WORK(KOMEG2+K-1), K = 1,NT2AMX)
              ENDIF
              CALL GPCLOSE(LUOME,'KEEP')
C
              CALL DAXPY(NT1AM(ISYMOP),-ONE,WORK(KRHO1),1,
     *                   WORK(KOMEG1),1)
              IF (.NOT. CCS ) THEN
                CALL DAXPY(NT2AM(ISYMOP),-ONE,WORK(KRHO2),1,
     *                     WORK(KOMEG2),1)
              ENDIF
C
              ECCP1 = DDOT(NT1AMX,WORK(KLAM),1,WORK(KOMEG1),1)
              ECCP2 = 0.0D0
              IF (.NOT.CCS) THEN
c                CALL CCLR_DIASCL(WORK(KOMEG2),0.5D0,ISYMTR)
                 ECCP2 = DDOT(NT2AMX,WORK(KLAM2),
     *                        1,WORK(KOMEG2),1)
              ENDIF
              IF (IPRINT .GE. 3) THEN
                 WRITE(LUPRI,*) 'Norm of omega1 in cc_fop:',
     *               DDOT(NT1AM(ISYMOP),WORK(KOMEG1),1,WORK(KOMEG1),1)
                 WRITE(LUPRI,*) 'Norm of omega2 in cc_fop:',
     *               DDOT(NT2AM(ISYMOP),WORK(KOMEG2),1,WORK(KOMEG2),1)
              ENDIF
              ECCL = ECCP1 + ECCP2
              ECCGRS = ECCGRS + ECCL
              WRITE(LUPRI,*)'Total <Lambda|H|CC> energy: ',ECCGRS
              WRITE(LUPRI,'(12X,A,F25.10)')
     *        'The singles contribution is:', ECCP1
              WRITE(LUPRI,'(12X,A,F25.10)')
     *        'The doubles contribution is:', ECCP2
C
      ENDIF
C
C---------------------------------------------------
C SLV98,OC solvent part 1 end
C---------------------------------------------------
C
      CALL FLSHFO(LUPRI)
C
  46  CONTINUE 
C
C-----------------------------------------------------------------
C     Calculate the coupled cluster energy using density matrices,
C     in order to check the unrelaxed CC-density.
C-----------------------------------------------------------------
C
      IF ((TSTDEN) .AND. (CCSD)) THEN
C
         CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND LUSIFC
C
         CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
         READ (LUSIFC) POTNUC
         CALL GPCLOSE(LUSIFC,'KEEP')
C
         KDENS = 1
         KWRK2 = KDENS + N2BST(ISYMOP)
         LWRK2 = LWORK - KWRK2
C
         IF (LWRK2 .LT. 0)
     *      CALL QUIT(' Too little workspace in cc_fop ')
C
         IOPT = 2
         CALL CC_DEN(POTNUC,WORK(KDENS),WORK(KWRK2),WORK(KWRK2),
     *               LWRK2,IOPT)
C
      ENDIF
C
      IF (RELORB) THEN
C
C---------------------------------------------------------
C        Set up diagonal block parts of Zeta-kappa-0, for
C        which no coupled equations needs to be solved,
C        and right hand side for ai-part of the equations.
C---------------------------------------------------------
C
         LENDEN = 2*NT1AMX    + NMATIJ(1)   + NMATAB(1)
     *          + 2*NCOFRO(1) + 2*NT1FRO(1) 
C
         KZKAM  = 1
         KETAAI = KZKAM  + LENDEN
         KEXVAL = KETAAI + NALLAI(1)
         KSOLUT = KEXVAL + 1
         KWRK2  = KSOLUT + NALLAI(1)
         LWRK2  = LWORK  - KWRK2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK2
            CALL QUIT('Insufficient memory for ETA(kappa) in CC_FOP')
         ENDIF
C
         CALL DZERO(WORK(KZKAM),LENDEN)
         CALL DZERO(WORK(KETAAI),NALLAI(1))
         CALL DZERO(WORK(KSOLUT),NALLAI(1))
C
         IF (MP2) THEN
            CALL MP2_ZKDIA(MODEL,WORK(KZKAM),WORK(KWRK2),LWRK2)
            CALL MP2_KANEW(MODEL,WORK(KETAAI),WORK(KZKAM),
     &           WORK(KWRK2),LWRK2)
         ELSE IF (CC2) THEN
            IOPT = 1
            CALL CC2_DEN(WORK(KETAAI),WORK(KZKAM),WORK(KWRK2),LWRK2,
     *                   IOPT)
            IOPT = 2
            CALL DZERO(WORK(KETAAI),NALLAI(1))
            CALL CC2_DEN(WORK(KETAAI),WORK(KZKAM),WORK(KWRK2),LWRK2,
     *                   IOPT)
         ELSE IF (CCSD) THEN
            IOPT = 1
            CALL CC_DEN(DUMMY,WORK(KETAAI),WORK(KZKAM),WORK(KWRK2),
     *                  LWRK2,IOPT)
         ENDIF
C
C------------------------------------------------------------
C        Open files for right hand side and solution vectors,
C        and residual vectors
C------------------------------------------------------------
C
         LUREVE = -2000
         LUSOVE = -2001
         LUGDVE = -2002
         CALL GPOPEN(LUREVE,'ZEKA0RES','UNKNOWN',' ','UNFORMATTED',
     *               IDUMMY,.FALSE.)
C
         CALL GPOPEN(LUSOVE,'ZEKA0SOL','UNKNOWN',' ','UNFORMATTED',
     *               IDUMMY,.FALSE.)
         REWIND(LUSOVE)
C
         CALL GPOPEN(LUGDVE,'ZEKA0RHS','UNKNOWN',' ','UNFORMATTED',
     *               IDUMMY,.FALSE.)
         REWIND(LUGDVE)
         CALL WRITT(LUGDVE,NALLAI(1),WORK(KETAAI))
CCH
         RHSNORM = DDOT(NALLAI(1),WORK(KETAAI),1,WORK(KETAAI),1)
         WRITE(LUPRI,*) 'CC_FOP> Norm of RHS vector:',RHSNORM
CCH
C
C        CALL HEADER('RHS vectors, MP2', -1)
C        CALL OUTPUT(WORK(KETAAI),1,NALLAI(1),1,1,NALLAI(1),1,1,LUPRI)
C
C----------------------------------------------------
C        Solve equations for ai-part of Zeta-kappa-0.
C----------------------------------------------------
C
         FASAVE = FABRHS
         NCOSAV = NCONF
C
         IF (DIRECT) CALL CCDFFOP
C
C     Close the 'AOTWOINT' file before entering the abarsp.
C
         IF (LUINTA .GT. 0) THEN
            CALL GPCLOSE(LUINTA,'KEEP')
            LUINTA = 0
         ENDIF
C
C     Open the 'SIRIFC' file before entering the abarsp.
C
         IF (LUSIFC .LE. 0) THEN
           CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         END IF
C
         CICLC  = .FALSE.
         HFCLC  = .TRUE. 
         TRPCLC = .FALSE.
         OOTV   = .FALSE.
         IOPSYM = 1
         EXCLC  = .FALSE.
         WORK(KEXVAL)= ZERO
         NEXVAL = 1
         NABATY = 1
         NABAOP = 1
         MXRM   = 40
         MXPHP  = 1
C
         FABRHS = .TRUE.
         NCONF  = 1
C
         IF (RHSNORM.GT.1.0D-12) THEN
C
          CALL HEADER('Solving for orbital relaxation vector',-1)
C
          LABEL = 'ETAKAPPA'
C
          CALL ABARSP(CICLC,HFCLC,TRPCLC,OOTV,IOPSYM,EXCLC,WORK(KEXVAL),
     *               NEXVAL,NABATY,NABAOP,LABEL,LUGDVE,LUSOVE,LUREVE,
     *               THRLEQ,MAXITE,IPRINT,MXRM,MXPHP,WORK(KWRK2),LWRK2)
C
          REWIND(LUSOVE)
          CALL READT(LUSOVE,NALLAI(1),WORK(KSOLUT))
C        
         ELSE
          CALL HEADER('Skipped solving for orbital relax. vector',-1)
          CALL DZERO(WORK(KSOLUT),NALLAI(1)) 
         END IF 
C
C        CALL HEADER('After ABARSP, MP2', -1)
C        CALL OUTPUT(WORK(KSOLUT),1,NALLAI(1),1,1,NALLAI(1),1,1,LUPRI)
C
         IF (LUINTA .LE. 0) THEN
           CALL GPOPEN(LUINTA,'AOTWOINT','UNKNOWN',' ','UNFORMATTED',
     *               IDUMMY,.FALSE.)
         END IF
C
C---------------------------------------------------------------
C        Unclosed leftover from respons-solver has to be closed.
C---------------------------------------------------------------
C
         CALL GPCLOSE(LUSOVE,'DELETE')
         CALL GPCLOSE(LUGDVE,'DELETE')
         CALL GPCLOSE(LUREVE,'DELETE')
C
         CALL GPCLOSE(LUSIFC,'KEEP')
         IF (LUPROP .GT. 0) CALL GPCLOSE(LUPROP,'KEEP')
         IF (LUINTM .GT. 0) CALL GPCLOSE(LUINTM,'DELETE')
C
C        save a copy on file CCL0___0
C
         IOPT = 4
         CALL CC_WRRSP('L0',0,1,IOPT,MODEL,WORK(KSOLUT),DUMMY,DUMMY,
     &                 WORK(KWRK2),LWRK2)
C
         FABRHS = FASAVE
         NCONF  = NCOSAV
C
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,48) 'Equations converged to residual less than:',
     *               THRLEQ
  48     FORMAT(3X,A42,2X,F8.6)
C
         CALL FLSHFO(LUPRI)
C
C------------------------------------------------------------------
C        Scale and reorder solution vector according to coupled
C        cluster standards, and write result to disc for later use.
C------------------------------------------------------------------
C
         CALL DSCAL(NALLAI(1),-ONE,WORK(KSOLUT),1)
C
         CALL CC_KABRE(WORK(KSOLUT),WORK(KZKAM),WORK(KWRK2),LWRK2)
C
         IF (IPRINT .GT. 9) THEN
            ZKNOR = DDOT(LENDEN,WORK(KZKAM),1,WORK(KZKAM),1)
            WRITE(LUPRI,*) ' '
            WRITE(LUPRI,*) 'Norm of zeta-kappa-0:', ZKNOR
         ENDIF
C
         LUBAR0 = -516
         CALL GPOPEN(LUBAR0,'CCKABAR0','UNKNOWN',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         REWIND(LUBAR0)
         WRITE(LUBAR0) (WORK(KZKAM+I-1), I = 1,LENDEN)
         CALL GPCLOSE(LUBAR0,'KEEP')
C
C-------------------------------------------------------------
C        Calculate the coupled cluster energy using density
C        matrices, in order to check the effective CC-density.
C-------------------------------------------------------------
C
         IF ((TSTDEN) .AND. (CCSD)) THEN
C
            CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &                  .FALSE.)
            REWIND LUSIFC
C
            CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
            READ (LUSIFC) POTNUC
            CALL GPCLOSE(LUSIFC,'KEEP')
C
            KSCRD = KWRK2
            KENDD = KSCRD + N2BST(ISYMOP)
            LENDD = LWORK - KENDD
C
            IF (LENDD .LT. 0)
     *         CALL QUIT(' Too little workspace in cc_fop ')
C
            IOPT = 3
            CALL CC_DEN(POTNUC,WORK(KSCRD),WORK(KENDD),WORK(KENDD),
     *                  LENDD,IOPT)
C
         ENDIF
C
      ELSE
C 
         KWRK2 = 1
C
      ENDIF
C
      KDENS = KWRK2     
      KWRK3 = KDENS + N2BST(ISYMOP)
      LWRK3 = LWORK - KWRK3
C
      IF (LWRK3 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK3
         CALL QUIT('Insufficient memory for one e- density in CC_FOP')
      ENDIF
C
      IF (.NOT.(CCSLV.OR.DIPMOM.OR.QUADRU.OR.NQCC.OR.
     *    RELCOR.OR.SECMOM.OR.TSTDEN.OR.(NAFOP.GT.0))) GOTO 47
C
C----------------------------------------------------------
C     Calculate one electron AO-density and CC nat.occ.num.
C----------------------------------------------------------
C

      ILSTNR = 1
      CALL CC_D1AO(WORK(KDENS),WORK(KZKAM),WORK(KWRK3),LWRK3,MODELPRI,
     *             LIST,ILSTNR,.TRUE.)
C
      IF ((FROIMP .OR. FROEXP) .AND. (.NOT. MP2)) THEN
C
C
        CALL CC_FCD1AO(WORK(KDENS),WORK(KWRK3),LWRK3,MODELPRI)
C
C
      ENDIF
C
      CALL FLSHFO(LUPRI)
C
      IF (IPRINT .GT. 50) THEN
         CALL AROUND('One electron density with orb.rel in cc_fop')
         CALL CC_PRFCKAO(WORK(KDENS),1)
      ENDIF
      CALL FLSHFO(LUPRI)
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SLV98,OC Solvent section 2
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF (CCSLV) THEN
C
        KETLM  = KWRK3
        KWRK4  = KETLM + 2*NLMCU
        LWRK4  = LWORK - KWRK4
        IF (LWRK4 .LT. 0) THEN
          WRITE(LUPRI,*) 'Needed:', KWRK4, 'Available:', LWORK
          CALL QUIT('Insufficient memory for solvent alloc in cc_fop')
        ENDIF
        CALL CC_SLV(WORK(KDENS),WORK(KETLM),WORK(KWRK4),LWRK4)
C
      ENDIF
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SLV98,OC End of solvent section 2.
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C---------------------------------------------------------------------
C     Calculate the simple one electron AO-density in CCS calculation.
C---------------------------------------------------------------------
C
  47  WRITE(LUPRI,*) ' '
C
      IF (CCS) THEN
C
         KDENS = 1
         KWRK3 = KDENS + N2BST(ISYMOP)
         LWRK3 = LWORK - KWRK3
C
         IF (LWRK3 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK3
            CALL QUIT('Insufficient memory for CCS AO-density in '//
     &                'CC_FOP')
         ENDIF
C
         CALL CCS_D1AO(WORK(KDENS),WORK(KWRK3),LWRK3)
         IF (FROIMP .OR. FROEXP) THEN
           CALL CC_FCD1AO(WORK(KDENS),WORK(KWRK3),LWRK3,MODELPRI)
         ENDIF
         IF (IPRINT .GT. 50) THEN
            CALL AROUND('CCS One electron density in cc_fop')
            CALL CC_PRFCKAO(WORK(KDENS),1)
         ENDIF
C
      ENDIF
C
      MODELPRI2 = '  Relaxed '//MODELPRI
      IF (.NOT. RELORB) MODELPRI2 = 'Unrelaxed '//MODELPRI
      IF (DIPMOM.OR.QUADRU.OR.NQCC.OR.RELCOR.OR.SECMOM.OR.
     *   (NAFOP.GT.0)) THEN
      CALL AROUND(MODELPRI2//' First-order one-electron properties: ')
      ENDIF
C
C=======================================
C     Calculate molecular dipole moment.
C=======================================
C
      IF (DIPMOM) THEN
C
         CALL AROUND(' Electric Dipole Moment ')
C
C-------------------------------------------
C        Calculate the nuclear contribution.
C-------------------------------------------
C
         IASGER = IPRINT - 9
         CALL DIPNUC(WORK(KWRK3),WORK(KWRK3),IASGER,.FALSE.)
C
         DO 100 IDIP = 1,3
C
            IF (IDIP .EQ. 1) LABEL = 'XDIPLEN '
            IF (IDIP .EQ. 2) LABEL = 'YDIPLEN '
            IF (IDIP .EQ. 3) LABEL = 'ZDIPLEN '
C
C----------------------------------
C           get property integrals.
C----------------------------------
C
            KONEP  = KWRK3
            KWRK4  = KONEP  + N2BST(ISYMOP)
            LWRK4  = LWORK  - KWRK4
C
            IF (LWRK4 .LT. 0) THEN
               WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK4
               CALL QUIT('Insufficient memory for DIPLEN-int. in '//
     &                   'CC_FOP')
            ENDIF
C
            CALL DZERO(WORK(KONEP),N2BST(ISYMOP))
            FF = 1.0D0
            ISY = -1
            CALL CC_ONEP(WORK(KONEP),WORK(KWRK4),LWRK4,FF,ISY,LABEL)
C
            IF (IPRINT .GT. 50) THEN
               CALL AROUND('One electron property integrals in cc_fop')
               CALL CC_PRFCKAO(WORK(KONEP),ISYMOP)
            ENDIF
C
C----------------------------------------------
C        Calculate the electronic contribution.
C----------------------------------------------
C
            IF (ISY .EQ. 1 ) THEN
               DIPME(IDIP) = -DDOT(N2BST(ISYMOP),WORK(KONEP),1,
     *                             WORK(KDENS),1)
            ELSE
               DIPME(IDIP) = 0
            ENDIF
            DIPMN(IDIP) = DIPMN(IDIP) + DIPME(IDIP)
C
  100    CONTINUE
C
C---------------------
C        Print result.
C---------------------
C
         IF (IASGER .GT. 0) THEN
            CALL HEADER('Electronic contribution to dipole moment',-1)
            CALL DP0PRI(DIPME)
         ENDIF
         CALL HEADER('Total Molecular Dipole Moment',-1)
         CALL DP0PRI(DIPMN)
C
         CALL FLSHFO(LUPRI)
C
      ENDIF
C
C===========================================
C     Calculate molecular quadrupole moment.
C===========================================
C
      IF (QUADRU) THEN
C
         CALL AROUND(' Electric Quadrupole Moment ')
C
C-------------------------------------------
C        Calculate the nuclear contribution.
C-------------------------------------------
C
         IOPT   = 1
         IASGER = -1
         CALL CCNUCQUA(WORK(KWRK3),LWRK3,IOPT,IASGER)
         CALL DZERO(QDREL,9)
C
         IJ = 0
         DO 110 I = 1,3
            DO 120 J = I,3
               IJ = IJ + 1
C
               IF (IJ .EQ. 1) LABEL = 'XXTHETA '
               IF (IJ .EQ. 2) LABEL = 'XYTHETA '
               IF (IJ .EQ. 3) LABEL = 'XZTHETA '
               IF (IJ .EQ. 4) LABEL = 'YYTHETA '
               IF (IJ .EQ. 5) LABEL = 'YZTHETA '
               IF (IJ .EQ. 6) LABEL = 'ZZTHETA '
C
C-------------------------------------
C              get property integrals.
C-------------------------------------
C
               KONEP  = KWRK3
               KWRK4  = KONEP  + N2BST(ISYMOP)
               LWRK4  = LWORK  - KWRK4
C
               IF (LWRK4 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK4
                  CALL QUIT('Insufficient memory for THETA-int. in '//
     &                      'CC_FOP')
               ENDIF
C
               CALL DZERO(WORK(KONEP),N2BST(ISYMOP))
               FF = 1.0D0
               ISY = -1
               CALL CC_ONEP(WORK(KONEP),WORK(KWRK4),LWRK4,FF,ISY,LABEL)
C
               IF (IPRINT .GT. 50) THEN
                  CALL AROUND('One electron property int. in cc_fop')
                  CALL CC_PRFCKAO(WORK(KONEP),ISYMOP)
               ENDIF
C
C-------------------------------------------------
C           Calculate the electronic contribution.
C-------------------------------------------------
C
               LENGTH = N2BST(ISYMOP)
C
               IF ( ISY .EQ. 1) THEN
                  CALL CCELQUA(WORK(KONEP),WORK(KDENS),LENGTH,I,J,QDREL)
               ENDIF
C
  120       CONTINUE
  110    CONTINUE
C
C------------------------
C        Reorder storing.
C------------------------
C
         CALL CC_QUAREO(QDREL,SKODE)
         CALL CC_QUAREO(QDRNUC,SKODN)
C
C---------------------
C        Print result.
C---------------------
C
         IF (IPRINT .GT. 9) THEN
            CALL HEADER('Nuclear contr. to quadrupole moment',-1)
            WRITE(LUPRI,474) 'X','Y','Z'
            CALL OUTPUT(SKODN,1,3,1,3,3,3,1,LUPRI)
            CALL HEADER('Electronic contr. to quadrupole moment',-1)
            WRITE(LUPRI,474) 'X','Y','Z'
            CALL OUTPUT(SKODE,1,3,1,3,3,3,1,LUPRI)
         ENDIF
         CALL DAXPY(9,-ONE,SKODE,1,SKODN,1)
         CALL HEADER('Total Molecular quadrupole moment',-1)
         WRITE(LUPRI,474) 'X','Y','Z'
         CALL OUTPUT(SKODN,1,3,1,3,3,3,1,LUPRI)
C
         CALL FLSHFO(LUPRI)
C
      ENDIF
C
C==================================================
C     Calculate electronic second moment of charge.
C==================================================
C
      IF (SECMOM) THEN
C
         CALL AROUND(' Electronic second moment of charge ')
C
         CALL DZERO(ELSEMO,9)
C
         IJ = 0
         DO 115 I = 1,3
            DO 125 J = I,3
               IJ = IJ + 1
C
               IF (IJ .EQ. 1) LABEL = 'XXSECMOM'
               IF (IJ .EQ. 2) LABEL = 'XYSECMOM'
               IF (IJ .EQ. 3) LABEL = 'XZSECMOM'
               IF (IJ .EQ. 4) LABEL = 'YYSECMOM'
               IF (IJ .EQ. 5) LABEL = 'YZSECMOM'
               IF (IJ .EQ. 6) LABEL = 'ZZSECMOM'
C
C-------------------------------------
C              get property integrals.
C-------------------------------------
C
               KONEP  = KWRK3
               KWRK4  = KONEP  + N2BST(ISYMOP)
               LWRK4  = LWORK  - KWRK4
C
               IF (LWRK4 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK4
                  CALL QUIT('Insufficient memory for SECMOM-int. in '//
     &                      'CC_FOP')
               ENDIF
C
               CALL DZERO(WORK(KONEP),N2BST(ISYMOP))
               FF = 1.0D0
               ISY = -1
               CALL CC_ONEP(WORK(KONEP),WORK(KWRK4),LWRK4,FF,ISY,LABEL)
C
               IF (IPRINT .GT. 50) THEN
                  CALL AROUND('One electron property int. in cc_fop')
                  CALL CC_PRFCKAO(WORK(KONEP),ISYMOP)
               ENDIF
C
C-------------------------------------------------
C           Calculate the electronic contribution.
C-------------------------------------------------
C
               LENGTH = N2BST(ISYMOP)
C
               IF (ISY.EQ.1) THEN
                 CALL CCELQUA(WORK(KONEP),WORK(KDENS),LENGTH,I,J,ELSEMO)
               ENDIF
C
  125       CONTINUE
  115    CONTINUE
C
C------------------------
C        Reorder storing.
C------------------------
C
         CALL CC_QUAREO(ELSEMO,SKODE)
C
C---------------------
C        Print result.
C---------------------
C
         WRITE(LUPRI,474) 'X','Y','Z'
         CALL OUTPUT(SKODE,1,3,1,3,3,3,1,LUPRI)
         CALL CC_TNSRAN(SKODE,WORK(KWRK3),LWRK3)
C
         CALL FLSHFO(LUPRI)
C
      ENDIF
C
  474 FORMAT(20X,A1,14X,A1,14X,A1)
C
C=======================================
C     Calculate electric field gradient.
C=======================================
C
      IF (NQCC) THEN
C
         CALL AROUND(' Electric Field Gradients ')
C
C-------------------------------------------
C        Calculate the nuclear contribution.
C-------------------------------------------
C
         IOPT   = 2
         IASGER = IPRINT - 5
         CALL CCNUCQUA(WORK(KWRK3),LWRK3,IOPT,IASGER)
C
C----------------------------------------------
C        Calculate the electronic contribution.
C----------------------------------------------
C
         LENGTH = N2BST(ISYMOP)
         CALL CCELEFG(WORK(KDENS),LENGTH,WORK(KWRK3),LWRK3,IASGER)
C
C---------------------
C        Print result.
C---------------------
C
         KDIAG = KWRK3
         KAXIS = KDIAG + 3*MXCENT
         KWRK4 = KAXIS + 9*MXCENT
         LWRK4 = LWORK - KWRK4
C
         IF (LWRK4 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK4
            CALL QUIT('Insufficient memory for EFG-results in CC_FOP')
         ENDIF
C
         IASGER = 2
         ICCPRI = 2
         CALL NQCRES(IASGER,WORK(KDIAG),WORK(KAXIS),ICCPRI)
C
         CALL FLSHFO(LUPRI)
C
      ENDIF
C
C==============================================
C     Calculate first-order relativistic energy
C     corrections within the DPT framework.
C==============================================
C
      IF (DPTECO) THEN
C
         CALL AROUND(' First-order DPT corrections to the ground-state'
     *               //' energy ')
C
         LABEL = 'DERXXPVP'
C
C----------------------------------------------------
C        Calculate the first and simplest correction.
C----------------------------------------------------
C
         KONEP  = KWRK3
         KWRK4  = KONEP  + N2BST(ISYMOP)
         LWRK4  = LWORK  - KWRK4
C
         IF (LWRK4 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK4
            CALL QUIT('Insufficient memory for DPT-integrals in CC_FOP')
         ENDIF
C
         CALL DZERO(WORK(KONEP),N2BST(ISYMOP))
         FF = 1.0D0
         ISY = 1
         CALL CC_ONEP(WORK(KONEP),WORK(KWRK4),LWRK4,FF,ISY,LABEL)
C
         IF (IPRINT .GT. 50) THEN
            CALL AROUND('Relativistic integrals in cc_fop')
            CALL CC_PRFCKAO(WORK(KONEP),ISYMOP)
         ENDIF
C
         DPTONE = ALPHA2*DDOT(N2BST(ISYMOP),WORK(KONEP),1,WORK(KDENS),1)
C
C----------------------------------------------------------
C        Calculate the second "one-electron term" - similar
C        to the reorthonormalization term of the gradient.
C----------------------------------------------------------
C
         RESONE = ZERO
         REORTH = ZERO
         IGROPT = 2
         CALL CC_GRAD(RESONE,REORTH,WORK(KWRK4),LWRK4,IGROPT)
         REORTH = ALPHA2*REORTH
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,136) 'DPTONE:', DPTONE
         WRITE(LUPRI,138) '------ '
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,136) 'DPTREO:', REORTH
         WRITE(LUPRI,138) '------ '
C
  136    FORMAT(9X,A7,F17.9)
  137    FORMAT(9X,A33,F17.9)
  138    FORMAT(9X,A7)
  139    FORMAT(9X,A32)
C
C------------------------------------------------------------
C        Calculate the "ordinary two-electron term" - similar
C        to the "simple" two-electron term of the gradient.
C------------------------------------------------------------
C
         DAR2SA = DAR2EL
         IF (DAR2EL) DAR2EL = .FALSE.
         IOPREL = 2
         CALL CC_2EEXP(WORK(KWRK4),LWRK4,IOPREL)
         DAR2EL = DAR2SA
C
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,136) 'DPTTWO:', WORK(KWRK4)
         WRITE(LUPRI,138) '------ '
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,137) 'Total first-order DPT correction:',
     *                DPTONE+REORTH+WORK(KWRK4)
         WRITE(LUPRI,139) '--------------------------------'
C
      ENDIF
C
C========================================================================
C     Standar scalar relativistic corrections to the ground-state energy.
C========================================================================
C
      IF (RELCOR) THEN
C
         CALL AROUND(' Pauli relativistic corrections to the'
     *               //' ground-state energy ')
C
         DO 130 IRC = 1,2
C
            IF (IRC .EQ. 1) LABEL = 'DARWIN  '
            IF (IRC .EQ. 2) LABEL = 'MASSVELO'
C
C-----------------------------
C           get the integrals.
C-----------------------------
C
            KONEP  = KWRK3
            KWRK4  = KONEP  + N2BST(ISYMOP)
            LWRK4  = LWORK  - KWRK4
C
            IF (LWRK4 .LT. 0) THEN
               WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK4
               CALL QUIT('Insufficient memory for Darwin-int. in '//
     &                   'CC_FOP')
            ENDIF
C
            CALL DZERO(WORK(KONEP),N2BST(ISYMOP))
            FF = 1.0D0
            ISY = 1
            CALL CC_ONEP(WORK(KONEP),WORK(KWRK4),LWRK4,FF,ISY,LABEL)
C
            IF (IPRINT .GT. 50) THEN
               CALL AROUND('Relativistic integrals in cc_fop')
               CALL CC_PRFCKAO(WORK(KONEP),ISYMOP)
            ENDIF
C
C-------------------------------------
C           Calculate the corrections.
C-------------------------------------
C
            IF (IRC .EQ. 1) THEN
               DARW = DDOT(N2BST(ISYMOP),WORK(KONEP),1,WORK(KDENS),1)
            ELSE IF (IRC .EQ. 2) THEN
               VELO = DDOT(N2BST(ISYMOP),WORK(KONEP),1,WORK(KDENS),1)
            ENDIF
C
  130    CONTINUE
C
C----------------------
C     Write out result.
C----------------------
C
      WRITE(LUPRI,*) ' '
      WRITE(LUPRI,131) '1-elec. Darwin term:', DARW
      WRITE(LUPRI,132) '------------------- '
      WRITE(LUPRI,*) ' '
      WRITE(LUPRI,131) 'Mass-Velocity term: ', VELO
      WRITE(LUPRI,132) '------------------  '
      WRITE(LUPRI,*) ' '
      WRITE(LUPRI,133) 'Mass-Velocity + 1-elec. Darwin terms:',DARW+VELO
      WRITE(LUPRI,134) '------------------------------------ '
C
  131 FORMAT(9X,A20,F17.9)
  132 FORMAT(9X,A20)
  133 FORMAT(9X,A37,1X,F17.9)
  134 FORMAT(9X,A37)             
C
      ENDIF
C
C--------------------------------------------------------------------
C     Calculate the relativistic two-electron Darwin term correction.
C--------------------------------------------------------------------
C
      IF (DAR2EL) THEN
         IF (RELCOR) THEN
            IOPREL = 1
            WORK(KWRK3) = DARW + VELO
         ELSE
            IOPREL = 0
         ENDIF
         CALL CC_2EEXP(WORK(KWRK3),LWRK3,IOPREL)
      ENDIF
C
C--------------------------------------------------------------
C     Section for general operator APROP represented by LABEL.
C     Note that only the electronic contribution is calculated.
C--------------------------------------------------------------
C
      DO 140 IOP = 1, NAFOP
C
         LABEL = PRPLBL(IAFOP(IOP))
C
         IF (IOP .EQ. 1) CALL AROUND( 
     *               ' Electronic contribution to operator ')
C
C--------------------------
C        get the integrals.
C--------------------------
C
         KONEP  = KWRK3
         KWRK4  = KONEP  + N2BST(ISYMOP)
         LWRK4  = LWORK  - KWRK4
C
         IF (LWRK4 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK4
            CALL QUIT('Insufficient memory for property integrals '//
     &                'in CC_FOP')
         ENDIF
C
         CALL DZERO(WORK(KONEP),N2BST(ISYMOP))
         FF = 1.0D0
         ISY = -1
         CALL CC_ONEP(WORK(KONEP),WORK(KWRK4),LWRK4,FF,ISY,LABEL)
C
         IF (IPRINT .GT. 50) THEN
            CALL AROUND('APROP integrals in cc_fop')
            CALL CC_PRFCKAO(WORK(KONEP),ISYMOP)
         ENDIF
C
C--------------------------------------------------------------------
C        Calculate the electronic contribution to the given property.
C--------------------------------------------------------------------
C
         IF (ISY.EQ.1) THEN
            PROP = DDOT(N2BST(ISYMOP),WORK(KONEP),1,WORK(KDENS),1)
         ELSE
            PROP = 0.0D0
         ENDIF
C
C-------------------------
C        Write out result.
C-------------------------
C
         WRITE(LUPRI,*) ' '
         IF (ISY.EQ.1) WRITE(LUPRI,141) LABEL//':', PROP
         IF (ISY.NE.1) WRITE(LUPRI,142) LABEL//':','zero by symmetry'
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) ' '
C
  141    FORMAT(20X,A9,1X,F12.8)
  142    FORMAT(20X,A9,1X,A)
C
  140 CONTINUE
C
C-------------------------------------------------------
C        Calculate energy for modifies CCSD(T) or CC(3).
C-------------------------------------------------------
C
      IF (CCPT .OR. CCP3) THEN
C
         CALL AROUND( ' Modified triples corrections ')
         CCSDT = .TRUE.
C
         IF (CCPT) THEN
            CC1BSV = CC1B
            CC1B   = .TRUE.
            CC1ASV = CC1A
            CC1A   = .TRUE.
         ENDIF
C
C---------------------------
C        Dynamic allocation.
C---------------------------
C
         KT1AM   = 1
         KOMEG1  = KT1AM   + NT1AM(ISYMOP)
         KOMEG2  = KOMEG1  + NT1AM(ISYMOP)
         IF (OMEGSQ) THEN
            KT2AM = KOMEG2
     *      + MAX(NT2AMX,NT2AM(ISYMOP),NT2AO(ISYMOP),NT2AOS(ISYMOP))
         ELSE
            KT2AM = KOMEG2
     *      + MAX(NT2AMX,NT2AM(ISYMOP),NT2AO(ISYMOP),2*NT2ORT(ISYMOP))
         ENDIF
         KSCR2   = KT2AM   + NT2AMX
         KEND1   = KSCR2   + NT2AMX + NT1AMX
         LWRK1   = LWORK   - KEND1
C
         IF ( LWRK1 .LT. 0  ) THEN
            CALL QUIT('Insufficient space in CC_FOP ')
         ENDIF
C
         IOPT = 3
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),WORK(KT2AM))
C
         IF ( IPRINT .GT. 50 ) THEN
           CALL AROUND( 'In CC_FOP:  (T1,T2)  vector before ' )
           CALL CC_PRP(WORK(KT1AM),WORK(KT2AM),1,1,1)
         ENDIF
C
         RSPIM = .FALSE.
         CALL CCRHSN(WORK(KOMEG1),WORK(KOMEG2),WORK(KT1AM),
     *               WORK(KT2AM),WORK(KEND1),LWRK1)
C
         RSPIM = .TRUE.
C
         IF (CCPT) THEN
            CC1B   = CC1BSV
            CC1A   = CC1ASV
         ENDIF
C
         KFOCKD = KEND1
         KEND1  = KFOCKD + NORBTS
         LWRK1  = LWORK  - KEND1
C
C----------------------------------------
C        Read canonical orbital energies.
C----------------------------------------
C
         CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND LUSIFC
C
         CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
         READ (LUSIFC) POTNUC,EMY,EACTIV,EMCSCF,ISTATE,ISPIN,NACTEL,
     *                 LSYM
C
         ESCF = EMCSCF
C
         CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
         READ (LUSIFC)
         READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBTS)
C
         CALL GPCLOSE(LUSIFC,'KEEP')
C
C-------------------------------------------------------------
C        Change symmetry-ordering of the Fock-matrix diagonal.
C-------------------------------------------------------------
C
         IF (FROIMP .OR. FROEXP)
     *       CALL CCSD_DELFRO(WORK(KFOCKD),WORK(KEND1),LWRK1)
C
         CALL FOCK_REORDER(WORK(KFOCKD),WORK(KEND1),LWRK1)
C
         ETY1 = 'CCSD'
         IT1 = 1
         ITER = 0
         CALL CCSD_ECCSD(WORK(KT1AM),WORK(KT2AM),WORK(KFOCKD),
     *                   WORK(KEND1),LWRK1,EN2,POTNUC,ESCF,
     *                   ETY1,IALPR,IT1,ITER)
C
         NTAMP = NT1AMX + NT2AMX
C
         KLAM  = KT2AM
         KEND1 = KLAM + NTAMP
         LWRK1 = LWORK   - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Needed:', KEND1, 'Available:', LWORK
            CALL QUIT('Insufficient memory for allocation in cc_fop')
         ENDIF
C
         IOPT   = 3
         CALL CC_RDRSP('L0',0,1,IOPT,MODEL,WORK(KLAM),
     *                 WORK(KLAM+NT1AMX)) 
C
         KLAM2 = KLAM + NT1AMX
C
         IF ( IPRINT .GT. 50 ) THEN
           CALL AROUND( 'In CC_FOP:  (L1,L2)  vector ' )
           CALL CC_PRP(WORK(KLAM),WORK(KLAM2),1,1,1)
         ENDIF
C
         CALL CCLR_DIASCL(WORK(KOMEG2),0.5D0,ISYMTR)
         ECCP1 = DDOT(NT1AMX,WORK(KLAM),1,WORK(KOMEG1),1)
         ECCP2 = DDOT(NT2AMX,WORK(KLAM2),1,WORK(KOMEG2),1)
C
         ETOT = EN2 + ECCP1 + ECCP2
         WRITE(LUPRI,'(//,21X,A)') 'Perturbative triples corrections'
         WRITE(LUPRI,'(21X,A,/)')  '--------------------------------'
         WRITE(LURES,'(//,21X,A)') 'Perturbative triples corrections'
         WRITE(LURES,'(21X,A,/)')  '--------------------------------'
         IF (CCPT) THEN
            WRITE(LUPRI,'(12X,A,F30.10)') 'Total energy MCCSD(T):',ETOT
            WRITE(LURES,'(12X,A,F30.10)') 'Total energy MCCSD(T):',
     *                                        ETOT
         ELSE
            WRITE(LUPRI,'(12X,A,F30.10)') 'Total energy MCC(3):',ETOT
            WRITE(LURES,'(12X,A,F30.10)') 'Total energy MCC(3):',ETOT
         ENDIF
C
         WRITE(LUPRI,'(12X,A,F25.10)')
     *        'The E4 doubles and triples:', ECCP2
         WRITE(LUPRI,'(12X,A,F25.10)')
     *        'The E5 singles and triples:', ECCP1
         WRITE(LURES,'(12X,A,F25.10)')
     *        'The E4 doubles and triples:', ECCP2
         WRITE(LURES,'(12X,A,F25.10)')
     *        'The E5 singles and triples:', ECCP1
         ECCGRS = ETOT
      ENDIF
C
C------------------------------------
C     Restore RELORB for CC2 and MP2.
C------------------------------------
C
      IF ((.NOT.RLORBS).AND.MP2) RELORB = RLORBS 
C
      CALL QEXIT('CC_FOP')
      RETURN
      END
c*DECK CC_ETA
      SUBROUTINE CC_ETA(ETA,WORK,LWORK)
C
C-----------------------------------------------------------------------------
C
C     Purpose: Calculate ETA vector.
C
C              Use F-hat and (iajb) on scratch.
C
C     Written by Ove Christiansen 22 November 1994
C
C-----------------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <ccorb.h>
#include <iratdef.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccinftap.h>
#include <ccsdinp.h>
C
      PARAMETER( TWO = 2.0D00 )
      CHARACTER*10 MODEL
      DIMENSION ETA(*),WORK(LWORK)
C
      CALL QENTER('CC_ETA')
C
      IF ( IPRINT .GT. 10 ) THEN
         CALL AROUND( 'IN CC_ETA: Constructing Eta vector ')
      ENDIF
C
      IF ( CCS ) THEN
         CALL DZERO(ETA,NT1AM(ISYMOP))
         CALL QEXIT('CC_ETA')
         RETURN
      ENDIF
C
C---------------------------------------------------
C     Make eta(ai,bj) from integrals (iajb) on disk.
C---------------------------------------------------
C
      REWIND(LUIAJB)
      CALL READI(LUIAJB,IRAT*NT2AM(ISYMOP),ETA(1+NT1AMX))
C
      IF (IPRINT .GT. 40 ) THEN
         CALL AROUND( 'In CC_ETA:  Integrals (ia|jb) ' )
         CALL CC_PRP(DUM,ETA(1+NT1AMX),1,0,1)
      ENDIF
C
      IOPTTCME = 1
      CALL CCSD_TCMEPK(ETA(1+NT1AMX),1.0D0,ISYMOP,IOPTTCME)
C
      KFOCK  = 1
      KT1AM  = KFOCK  + N2BST(ISYMOP)
      KLAMDP = KT1AM  + NT1AM(ISYMOP)
      KLAMDH = KLAMDP + NLAMDT
      KEND1  = KLAMDH + NLAMDT
      LWRK1  = LWORK  - KEND1
C
C----------------------------------------------------
C     Make eta(ai) from AO fock matrix store on disk.
C----------------------------------------------------
C
      IF ( RSPIM ) THEN
C
         LUFCK = 0
         CALL GPOPEN(LUFCK,'CC_FCKH','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     *               .FALSE.)
         REWIND(LUFCK )
         READ (LUFCK )(WORK(KFOCK + I-1),I = 1,N2BST(ISYMOP))
         CALL GPCLOSE(LUFCK,'KEEP')
C
      ENDIF
C
      IF (IPRINT .GT.140) THEN
         CALL AROUND( 'Usual Fock AO matrix' )
         ISYFAO = 1
         CALL CC_PRFCKAO(WORK(KFOCK),ISYFAO)
      ENDIF
C
      CALL DZERO(WORK(KT1AM),NT1AM(1))
      IF (.NOT.(CCS.OR.CCP2)) THEN
         IOPT = 1
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),DUMMY)
      ENDIF
C
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     *            WORK(KEND1),LWRK1)
C
      ISYFAO = 1
      ISYMPA = 1
      ISYMHO = 1
C
      CALL CC_FCKMO(WORK(KFOCK),WORK(KLAMDP),WORK(KLAMDH),
     *              WORK(KEND1),LWRK1,ISYFAO,ISYMPA,ISYMHO)
C
      IF (IPRINT .GT. 50) THEN
         CALL AROUND( 'In CC_ETA: Fock MO matrix' )
         CALL CC_PRFCKMO(WORK(KFOCK),ISYMOP)
      ENDIF
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMA = MULD2H(ISYMI,ISYMOP)
C
         DO 110 I = 1,NRHF(ISYMI)
C
            DO 120 A = 1,NVIR(ISYMA)
C
               KOFF1 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
               KOFF2 = IFCVIR(ISYMI,ISYMA) + NORB(ISYMI)*(A - 1) + I
C
               ETA(KOFF1) = WORK(KOFF2)
C
  120       CONTINUE
  110    CONTINUE
C
  100 CONTINUE
C
      CALL DSCAL(NT1AMX+NT2AMX,TWO,ETA,1)
      IF (IPRINT .GT. 40 ) THEN
         CALL AROUND( 'In CC_ETA:  Eta vector ' )
         CALL CC_PRP(ETA(1),ETA(1+NT1AMX),1,1,1)
      ENDIF
C
      IF ( IPRINT .GT. 10 ) THEN
         ETA1 = DDOT(NT1AMX,ETA(1),1,ETA(1),1)
         ETA2 = DDOT(NT2AMX,ETA(1+NT1AMX),1,ETA(1+NT1AMX),1)
         WRITE(LUPRI,*) 'Norm of eta1: ',ETA1
         WRITE(LUPRI,*) 'Norm of eta2: ',ETA2
         CALL AROUND( 'END OF CC_ETA ')
      ENDIF
C
      CALL QEXIT('CC_ETA')
      RETURN
      END
C  /* Deck cc_d1orre */
      SUBROUTINE CC_D1ORRE(D1AO,ZKAM,WORK,LWORK)
C
C     Written by Asger Halkier 4/4 - 1996
C
C     Version: 1.0
C
C     Purpose: To add the orbital relaxation term to the
C              CC one electron density in AO basis!
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION D1AO(*), ZKAM(*), WORK(LWORK)
#include <inftap.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      CALL QENTER('CC_D1ORRE')
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      LENGHT = MAX(NLAMDT,NLAMDS)
C
      KCTRAN = 1
      KEND1  = KCTRAN + LENGHT
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for first allocation in '//
     &             'CC_D1ORRE')
      ENDIF
C
C----------------------------------------------------
C     Read MO-coefficient matrix from interface file.
C----------------------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KCTRAN+I-1), I=1,NLAMDS)
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C------------------------------------------------------------
C     Reorder MO-coefficient matrix to lampda matrix storage.
C------------------------------------------------------------
C
      CALL CMO_REORDER(WORK(KCTRAN),WORK(KEND1),LWRK1)
C
      DO 100 ISYM = 1,NSYM
C
C----------------------------------
C        Work space allocation two.
C----------------------------------
C
         KSCR  = KEND1
         KEND2 = KSCR  + NBAS(ISYM)*NRHF(ISYM)
         LWRK2 = LWORK - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient work for second allocation in '//
     &                'CC_D1ORRE')
         ENDIF
C
         CALL DZERO(WORK(KSCR),NBAS(ISYM)*NRHF(ISYM))
C
C------------------------------------
C        Calculate the contributions.
C------------------------------------
C
         KOFF1  = KCTRAN + ILMVIR(ISYM)
         KOFF2  = IT1AM(ISYM,ISYM) + 1
C
         NTOTAL = MAX(NBAS(ISYM),1)
         NTOTA  = MAX(NVIR(ISYM),1)
C
         CALL DGEMM('N','N',NBAS(ISYM),NRHF(ISYM),NVIR(ISYM),ONE,
     *              WORK(KOFF1),NTOTAL,ZKAM(KOFF2),NTOTA,ZERO,
     *              WORK(KSCR),NTOTAL)
C
         KOFF3  = KCTRAN + ILMRHF(ISYM)
         KOFF4  = IAODIS(ISYM,ISYM) + 1
C
         NTOTAL = MAX(NBAS(ISYM),1)
         NTOTBE = MAX(NBAS(ISYM),1)
C
         CALL DGEMM('N','T',NBAS(ISYM),NBAS(ISYM),NRHF(ISYM),TWO,
     *              WORK(KSCR),NTOTAL,WORK(KOFF3),NTOTBE,ONE,
     *              D1AO(KOFF4),NTOTAL)
C
  100 CONTINUE
C
      CALL QEXIT('CC_D1ORRE')
C
      RETURN
      END
C  /* Deck ccdffop */
      subroutine CCDFFOP
C
C     Written by Asger Halkier 5/4 - 1996
C
C     Version: 1.0
C
C     Purpose: Set flags for response solver properly for integral
C              direct calculations!
C
C
#include <implicit.h>
#include <mxcent.h>
#include <abainf.h>
#include <inftra.h>
C
      CALL QENTER('CCDFFOP')
C
      DODRCT = .TRUE.
      USEDRC = .TRUE.
C
      CALL QEXIT('CCDFFOP')
C
      RETURN
      END
C  /* Deck ccnucqau */
      subroutine CCNUCQUA(WORK,LWORK,IOPT,IASGER)
C
C     Written by Asger Halkier 9/4 - 1996
C
C     Version: 1.0
C
C     Purpose: Calculate the nuclear contribution to the
C              molecular quadrupole moment (based on the
C              equivalent ABACUS-routines)!
C
C
#include <implicit.h>
#include <iratdef.h>
#include <priunit.h>
#include <mxcent.h>
      DIMENSION WORK(LWORK)
#include <cbiher.h>
#include <orgcom.h>
#include <nuclei.h>
C
      CALL QENTER('CCNUCQUA')
C
      KGEOM = 1
      KMASS = KGEOM + 3*(NATOMS + NFLOAT)
      KNAT  = KMASS + NATOMS + NFLOAT
      KEND1 = KNAT  + (NATOMS + NFLOAT + 1)/IRAT
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for allocation in CCNUCQUA')
      ENDIF
C
      DO 100 I = 1, MXCENT
         ISOTOP(I) = 1
  100 CONTINUE
C
      CALL CMMASS(WORK(KGEOM),WORK(KMASS),WORK(KNAT),ISOTOP,IASGER)
      CALL DZERO(CMXYZ,3)
C
      IF (IOPT .EQ. 1) THEN
         CALL NUCQDR(WORK(KGEOM),CMXYZ,LUPRI,IASGER)
      ELSE IF (IOPT .EQ. 2) THEN
         CALL NUCNQC(WORK(KGEOM),LUPRI,IASGER)
      ENDIF
C
      CALL QEXIT('CCNUCQUA')
C
      RETURN
      END
C  /* Deck ccelqau */
      subroutine CCELQUA(XONEP,DENS,LENGTH,I,J,RESVEC)
C
C     Written by Asger Halkier 9/4 - 1996
C
C     Version: 1.0
C
C     Purpose: Calculate the electronic contribution to the
C              molecular quadrupole moment (based on the
C              equivalent ABACUS-routines)!
C
C
#include <implicit.h>
#include <mxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
      DIMENSION XONEP(*), DENS(*), RESVEC(3,3)
#include <symmet.h>
#include <quadru.h>
C
      CALL QENTER('CCELQUA')
C
      RESVEC(IPTAX(J,1),IPTAX(I,1)) = DDOT(LENGTH,XONEP,1,DENS,1)
      RESVEC(IPTAX(I,1),IPTAX(J,1)) = DDOT(LENGTH,XONEP,1,DENS,1)
C
      CALL QEXIT('CCELQUA')
C
      RETURN
      END
C  /* Deck cc_quareo */
      subroutine CC_QUAREO(QORI,QNEW)
C
C     Written by Asger Halkier 19/3 - 1998
C
C     Version: 1.0
C
C     Purpose: Reorder quadrupole and second moment tensors to
C              CC storing.
C
C
#include <implicit.h>
#include <mxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
      PARAMETER (ZERO = 0.0D0)
      DIMENSION QORI(3,3), QNEW(3,3)
#include <symmet.h>
#include <quadru.h>
C
      CALL QENTER('CC_QUAREO')
C
      DO 100 I = 1,3
         DO 110 J = 1,3
            QNEW(I,J) = ZERO
  110 CONTINUE
  100 CONTINUE
C
      DO 120 I = 1,3
         DO 130 J = I,3
            QNEW(I,J) = QORI(IPTAX(J,1),IPTAX(I,1))
            QNEW(J,I) = QORI(IPTAX(I,1),IPTAX(J,1))
  130 CONTINUE
  120 CONTINUE
C
      CALL QEXIT('CC_QUAREO')
C
      RETURN
      END
C  /* Deck ccelefg */
      subroutine CCELEFG(DENS,LENGTH,WORK,LWORK,IASGER)
C
C     Written by Asger Halkier 16/4 - 1996
C
C     Version: 1.0
C
C     Purpose: Calculate the electronic contribution to the
C              electric field gradients (based on the
C              equivalent ABACUS-routines)!
C
C
#include <implicit.h>
#include <iratdef.h>
#include <priunit.h>
#include <mxcent.h>
      DIMENSION DENS(*), WORK(LWORK)
#include <nuclei.h>
C
      CALL QENTER('CCELEFG')
C
      NCOMP = 9*NUCDEP
C
C---------------------------
C     Work space allocation.
C---------------------------
C
      KDOTPR = 1
      KCAINT = KDOTPR + NCOMP
      KEND1  = KCAINT + LENGTH
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for initial allocation '//
     &             'in CCELEFG')
      ENDIF
C
      CALL DZERO(WORK(KDOTPR),NCOMP)
C
C------------------------------------------------------------------
C     Calculate contraction of density and cartesian efg-integrals.
C------------------------------------------------------------------
C
      ITYPE = 30
      CALL CCELEFG1(WORK(KDOTPR),DENS,WORK(KCAINT),
     *              WORK(KEND1),LWRK1,NCOMP,LENGTH,IASGER)
C
C-------------------------------------------
C     Calculate the contribution to the EFG.
C-------------------------------------------
C
      KSCR1 = KEND1
      KEND2 = KSCR1 + 9*NUCDEP
      LWRK2 = LWORK - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
         CALL QUIT('Insufficient memory for final allocation '//
     &             'in CCELEFG')
      ENDIF
C
      CALL NQCEL(WORK(KSCR1),WORK(KDOTPR),NCOMP,IASGER)
C
      CALL QEXIT('CCELEFG')
C
      RETURN
      END
C  /* Deck ccelefg1 */
      subroutine CCELEFG1(DOTPRO,DENS,EFGINT,WORK,LWORK,
     *                    NCOMP,LENGTH,IASGER)
C
C     Written by Asger Halkier 16/4 - 1996
C
C     Version: 1.0
C
C     Purpose: To read in appropriate cartesian electric field
C              gradient integrals and contract these with the
C              one electron density matrix (Based on the equivalent
C              ABACUS routines)!
C
C     Merge to Dalton1.0 Ove 16-4-1997
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <mxcent.h>
#include <mxorb.h>
      DIMENSION DENS(*), EFGINT(*), DOTPRO(NCOMP), WORK(LWORK)
      CHARACTER*8 LABEL
#include <nuclei.h>
#include <symmet.h>
#include <chrxyz.h>
#include <chrnos.h>
#include <ibtfun.h>
C
      CALL QENTER('CCELEFG1')
C
C---------------------------
C     Set up loop structure.
C---------------------------
C
      ITYP = 0
C
      DO 100 IATOM = 1,NUCIND
         DO 110 ICOOR1 = 1,3
            DO 120 ICOOR2 = ICOOR1,3
C
               ISYMIJ = IBTXOR(ISYMAX(ICOOR1,1),ISYMAX(ICOOR2,1))
C
               IOFF = 0

               DO 130 IREPC = 0, MAXREP
C
                  IF (IBTAND(ISTBNU(IATOM),IBTXOR(IREPC,ISYMIJ))
     *                .EQ.0) THEN
C
C---------------------------------------------------------------------
C                    Get the integrals and contract with integrals.
C---------------------------------------------------------------------
C
                     IOFF = IOFF + 1 
                     ITYP = ITYP + 1
C
                     LABEL = CHRXYZ(ICOOR1)//CHRXYZ(ICOOR2)//'EFG'//
     *                  CHRNOS(IATOM/10)//CHRNOS(IATOM)//CHRNOS(IOFF)
C
                     CALL DZERO(EFGINT,LENGTH)
                     FF = 1.0D0
                     ISY = -1
                     CALL CC_ONEP(EFGINT,WORK,LWORK,FF,ISY,LABEL)
C
                     IF (IASGER .GT. 45) THEN
                        CALL AROUND('Cartesian EFG-int. in cc_fop')
                        CALL CC_PRFCKAO(EFGINT,ISY)
                     ENDIF
C
                     IF (ISY .EQ. 1) THEN
                        DOTPRO(ITYP) = DDOT(LENGTH,DENS,1,EFGINT,1)
                     ELSE
                        DOTPRO(ITYP) = 0.0D0
                     ENDIF
C
                  ENDIF
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CCELEFG1')
C
      RETURN
      END
C  /* Deck ccs_d1ao */
      SUBROUTINE CCS_D1AO(AODEN,WORK,LWORK)
C
C     Written by Asger Halkier 17/4 - 1996
C
C     Version: 1.0
C
C     Purpose: To set up the one electron AO-density in case
C              of a CCS calculation (equal to HF density)!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION AODEN(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      CALL QENTER('CCS_D1AO')
C
C---------------------------
C     Work space allocation.
C---------------------------
C
      KONEAI = 1
      KONEAB = KONEAI + NT1AMX
      KONEIJ = KONEAB + NMATAB(1)
      KONEIA = KONEIJ + NMATIJ(1)
      KT1AM  = KONEIA + NT1AMX
      KLAMDH = KT1AM  + NT1AMX
      KLAMDP = KLAMDH + NLAMDT
      KEND1  = KLAMDP + NLAMDT
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for work allocation '//
     &             'in CCS_D1AO')
      ENDIF
C
C--------------------------------------------------------------
C     Initialize arrays (note that the t1-amplitudes are zero).
C--------------------------------------------------------------
C
      CALL DZERO(WORK(KONEAI),NT1AMX)
      CALL DZERO(WORK(KONEAB),NMATAB(1))
      CALL DZERO(WORK(KONEIJ),NMATIJ(1))
      CALL DZERO(WORK(KONEIA),NT1AMX)
      CALL DZERO(WORK(KT1AM),NT1AMX)
C
C-----------------------
C     Set up MO-density.
C-----------------------
C
      DO 100 ISYM = 1,NSYM
         DO 110 I = 1,NRHF(ISYM)
C
            NII = IMATIJ(ISYM,ISYM) + NRHF(ISYM)*(I - 1) + I
C
            WORK(KONEIJ + NII - 1) = TWO
C
  110    CONTINUE
  100 CONTINUE
C
C-------------------------------
C     Get MO coefficient matrix.
C-------------------------------
C
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND1),
     *            LWRK1)
C
C-----------------------------------
C     Transform density to AO basis.
C-----------------------------------
C
      CALL DZERO(AODEN,N2BST(1))
C
      ISDEN = 1
      CALL CC_DENAO(AODEN,ISDEN,WORK(KONEAI),WORK(KONEAB),
     *              WORK(KONEIJ),WORK(KONEIA),ISDEN,WORK(KLAMDP),1,
     *              WORK(KLAMDH),1,WORK(KEND1),LWRK1)
C
      CALL QEXIT('CCS_D1AO')
C
      RETURN
      END
C  /* Deck mp_lam */
      SUBROUTINE MP_LAM(TBAM,WORK,LWORK)
C
C     Written by Asger Halkier 6/9 - 1996
C
C     Version: 1.0
C
C     Purpose: To set up the zero'th order Lagrangian multipliers
C              in the MP2 case.
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION TBAM(*), WORK(LWORK)
#include <maxorb.h>
#include <ccorb.h>
#include <iratdef.h>
#include <inftap.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccinftap.h>
#include <ccsdio.h>
#include <ccsdinp.h>
C
      CALL QENTER('MP_LAM')
C
C-----------------------------------------------------------------
C     Read integrals (ia|jb) from disc (file always assumed open).
C-----------------------------------------------------------------
C
      REWIND(LUIAJB)
      READ(LUIAJB) (TBAM(NT1AMX + I), I = 1,NT2AM(ISYMOP))
C
C-----------------------------------------------
C     Take two coulomb minus exchange on vector.
C-----------------------------------------------
C
      IOPTTCME = 1
      CALL CCSD_TCMEPK(TBAM(1+NT1AMX),1.0D0,ISYMOP,IOPTTCME)
C
C---------------------------
C     Work space allocation.
C---------------------------
C
      KFOCKD = 1
      KEND1  = KFOCKD + NORBTS
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need:', KEND1, 'Available:', LWORK
         CALL QUIT('Insufficient memory for allocation in MP_LAM')
      ENDIF
C
C-------------------------------------
C     Read canonical orbital energies.
C-------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND (LUSIFC)
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD + I - 1), I = 1,NORBTS)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C----------------------------------------------------------------
C     Change symmetry ordering of the canonical orbital energies.
C----------------------------------------------------------------
C
      IF (FROIMP .OR. FROEXP)
     *    CALL CCSD_DELFRO(WORK(KFOCKD),WORK(KEND1),LWRK1)
C
      CALL FOCK_REORDER(WORK(KFOCKD),WORK(KEND1),LWRK1)
C
C-------------------------------------
C     Divide with orbital differences.
C-------------------------------------
C
      CALL CCSD_GUESS(TBAM(1),TBAM(1+NT1AMX),WORK(KFOCKD),IPRINT)
C
C-----------------------------------------
C     Final scalings for obtaining result.
C-----------------------------------------
C
      CALL DSCAL(NT2AM(ISYMOP),TWO,TBAM(1+NT1AMX),1)
C
      CALL QEXIT('MP_LAM')
C
      RETURN
      END
C  /* Deck mp2_kari */
      SUBROUTINE MP2_KARI(ETAAI,WORK,LWORK)
C
C     Written by Asger Halkier 7/9 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the right hand side ETAAI for the
C              equations for the zero'th order orbital rotation
C              multipliers in CCPT2 calculations.
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <maxash.h>
#include <mxorb.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION INDEXA(MXCORB)
      DIMENSION ETAAI(*), WORK(LWORK)
#include <ccorb.h>
#include <infind.h>
#include <blocks.h>
#include <ccsdinp.h>
#include <ccinftap.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <distcl.h>
#include <cbieri.h>
#include <eribuf.h>
#include <cclr.h>
C
      CALL QENTER('MP2_KARI')
C
      CALL HEADER('Constructing right-hand-side for MP2-kappa-0(ai)',-1)
C
      TIMETO = ZERO
      TIMETO = SECOND()
C
C----------------------------------------------------------------------
C     Both and t-vectors and tbar-vectors (zeta) are totally symmetric.
C----------------------------------------------------------------------
C
      ISYMTR = 1
      ISYMOP = 1
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C  
      KT2AM  = 1
      KXMAT  = KT2AM  + NT2AMX
      KYMAT  = KXMAT  + NMATIJ(1)
      KXTMAT = KYMAT  + NMATAB(1)
      KYTMAT = KXTMAT + NMATIJ(1)
      KDENSI = KYTMAT + NMATAB(1)
      KFOCK  = KDENSI + N2BAST
      KLAMDP = KFOCK  + N2BST(ISYMOP)
      KLAMDH = KLAMDP + NLAMDT
      KZ2AM  = KLAMDH + NLAMDT
      KT1AM  = KZ2AM  + NT2SQ(1)
      KZ1AM  = KT1AM  + NT1AMX
      KEND1  = KZ1AM  + NT1AMX
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for initial allocation '//
     &             'in MP2_KARI')
      ENDIF
C
C----------------------------------------
C     Read zero'th order zeta amplitudes.
C----------------------------------------
C
      IOPT   = 3
      CALL CC_RDRSP('L0',0,1,IOPT,MODEL,WORK(KZ1AM),WORK(KZ2AM))
C
      KEND1 = KZ1AM
      LWRK1 = LWORK  - KEND1
C
C--------------------------------
C     Square up zeta2 amplitudes.
C--------------------------------
C
      CALL DCOPY(NT2AMX,WORK(KZ2AM),1,WORK(KT2AM),1)
      CALL CC_T2SQ(WORK(KT2AM),WORK(KZ2AM),1)
C
C
C-------------------------------------------
C     Read zero'th order cluster amplitudes.
C-------------------------------------------
C
      IOPT = 3
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),WORK(KT2AM))
C
C----------------------------------
C     Calculate the lambda matrices.
C----------------------------------
C
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND1),
     *            LWRK1)
C
      KEND1 = KT1AM
      LWRK1 = LWORK  - KEND1
C
C
C--------------------------------------------------------
C     Calculate X-intermediate of tbar- and t-amplitudes.
C--------------------------------------------------------
C
      CALL CC_XI(WORK(KXMAT),WORK(KZ2AM),ISYMTR,WORK(KT2AM),ISYMOP,
     *             WORK(KEND1),LWRK1)
C
C--------------------------------------------------------
C     Calculate Y-intermediate of tbar- and t-amplitudes.
C--------------------------------------------------------
C
      CALL CC_YI(WORK(KYMAT),WORK(KZ2AM),ISYMTR,WORK(KT2AM),ISYMOP,
     *           WORK(KEND1),LWRK1)
C
C---------------------------------------
C     Set up 2C-E of cluster amplitudes.
C---------------------------------------
C
      ISYOPE = 1
      IOPTTCME = 1
      CALL CCSD_TCMEPK(WORK(KT2AM),1.0D0,ISYOPE,IOPTTCME)
C
C--------------------------------------------------------------------
C     Set up special modified amplitudes needed in the integral loop.
C     (By doing it this way, we only need one packed vector in core
C     along with the integral distribution in the delta loop.)
C--------------------------------------------------------------------
C
      IOPT   = 3
      CALL CC_RDRSP('L0',0,1,IOPT,MODEL,WORK(KZ1AM),WORK(KZ2AM))
C
      CALL DSCAL(NT2AMX,TWO,WORK(KT2AM),1)
      CALL DAXPY(NT2AMX,ONE,WORK(KZ2AM),1,WORK(KT2AM),1)
C
C----------------------------------
C     Calculate the density matrix.
C----------------------------------
C
      ISYMH = 1
      IC    = 1
      CALL CC_AODENS(WORK(KLAMDP),WORK(KLAMDH),WORK(KDENSI),ISYMH,
     *               IC,WORK(KEND1),LWRK1)
C
      KEND1 = KLAMDH
      LWRK1 = LWORK  - KEND1
C
C------------------------------------------------
C     Read one-electron integrals in Fock-matrix.
C------------------------------------------------
C
      CALL CCRHS_ONEAO(WORK(KFOCK),WORK(KEND1),LWRK1)
C
C------------------------------------------------------- 
C     Calculate special modified X- and Y-intermediates.
C------------------------------------------------------- 
C
      CALL DCOPY(NMATAB(1),WORK(KYMAT),1,WORK(KYTMAT),1)
      CALL DCOPY(NMATIJ(1),WORK(KXMAT),1,WORK(KXTMAT),1)
      CALL CC_EITR(WORK(KYTMAT),WORK(KXTMAT),WORK(KEND1),LWRK1,1)
      CALL DAXPY(NMATAB(1),ONE,WORK(KYMAT),1,WORK(KYTMAT),1)
      CALL DAXPY(NMATIJ(1),ONE,WORK(KXMAT),1,WORK(KXTMAT),1)
C
C-----------------------------------
C     Start the loop over integrals.
C-----------------------------------
C
      KENDS2 = KEND1
      LWRKS2 = LWRK1
C
      IF (DIRECT) THEN
         IF (HERDIR) THEN
           CALL HERDI1(WORK(KEND1),LWRK1,IPRERI)
         ELSE 
           KCCFB1 = KEND1
           KINDXB = KCCFB1 + MXPRIM*MXCONT
           KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
           LWRK1  = LWORK  - KEND1
           CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     *                 KODPP1,KODPP2,KRDPP1,KRDPP2,
     *                 KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     *                 WORK(KEND1),LWRK1,IPRERI)
           KEND1 = KFREE
           LWRK1 = LFREE
         END IF
         NTOSYM = 1
      ELSE
         NTOSYM = NSYM
      ENDIF
C
      KENDSV = KEND1
      LWRKSV = LWRK1
C
      ICDEL1 = 0
      DO 100 ISYMD1 = 1,NTOSYM
C
         IF (DIRECT) THEN
            IF (HERDIR) THEN
              NTOT = MAXSHL
            ELSE
              NTOT = MXCALL
            END IF
         ELSE
            NTOT = NBAS(ISYMD1)
         ENDIF
C
         DO 110 ILLL = 1,NTOT
C
C---------------------------------------------
C           If direct calculate the integrals.
C---------------------------------------------
C
            IF (DIRECT) THEN
C
               KEND1 = KENDSV
               LWRK1 = LWRKSV
C
               DTIME  = SECOND()
               IF (HERDIR) THEN
                 CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,
     &                       IPRERI)
               ELSE
                 CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     *                       WORK(KODCL1),WORK(KODCL2),
     *                       WORK(KODBC1),WORK(KODBC2),
     *                       WORK(KRDBC1),WORK(KRDBC2),
     *                       WORK(KODPP1),WORK(KODPP2),
     *                       WORK(KRDPP1),WORK(KRDPP2),
     *                       WORK(KCCFB1),WORK(KINDXB),
     *                       WORK(KEND1), LWRK1,IPRERI)
               END IF
               DTIME   = SECOND() - DTIME
               TIMHE2 = TIMHE2 + DTIME
C
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
                  CALL QUIT('Insufficient memory for integrals '//
     &                      'in MP2_KARI')
               END IF
C
            ELSE
               NUMDIS = 1
            ENDIF
C
C-----------------------------------------------------
C           Loop over number of distributions in disk.
C-----------------------------------------------------
C
            DO 120 IDEL2 = 1,NUMDIS
C
               IF (DIRECT) THEN
                  IDEL  = INDEXA(IDEL2)
                  ISYMD = ISAO(IDEL)
               ELSE
                  IDEL  = IBAS(ISYMD1) + ILLL
                  ISYMD = ISYMD1
               ENDIF
C
C----------------------------------------
C              Work space allocation two.
C----------------------------------------
C
               ISYDIS = MULD2H(ISYMD,ISYMOP)
C
               KXINT  = KEND1
               KEND2  = KXINT + NDISAO(ISYDIS)
               LWRK2  = LWORK - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
                  CALL QUIT('Insufficient memory for integrals '//
     &                      'in MP2_KARI')
               ENDIF
C
C--------------------------------------------
C              Read AO integral distribution.
C--------------------------------------------
C
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                     WORK(KRECNR),DIRECT)
C
C-------------------------------------------
C              Calculate the AO-Fock matrix.
C-------------------------------------------
C
               ISYDEN = 1
               CALL CC_AOFOCK(WORK(KXINT),WORK(KDENSI),WORK(KFOCK),
     *                        WORK(KEND2),LWRK2,IDEL,ISYMD,ISYDEN)
C
C------------------------------------------
C              Work space allocation three.
C------------------------------------------
C
               KDSRHF = KEND2
               K3OINT = KDSRHF + NDSRHF(ISYMD)
               KSCRTI = K3OINT + NMAIJK(ISYDIS)
               KEND3  = KSCRTI + NT2BCD(ISYDIS)
               LWRK3  = LWORK  - KEND3
C
               IF (LWRK3 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND3,'Available : ',LWORK
                  CALL QUIT('Insufficient memory for integrals '//
     &                      'in MP2_KARI')
               ENDIF
C
C---------------------------------------------------------------------
C              Calculate partially backtransformed modified amplitude.
C---------------------------------------------------------------------
C
               CALL CC_TI(WORK(KSCRTI),ISYMD,WORK(KT2AM),ISYMOP,
     *                    WORK(KLAMDP),1,WORK(KEND3),LWRK3,IDEL,ISYMD)
C
C--------------------------------------------------------
C              Transform one index in the integral batch.
C--------------------------------------------------------
C
               CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KLAMDP),ISYMOP,
     *                     WORK(KEND3),LWRK3,ISYDIS)
C
C------------------------------------------------------------------
C              Calculate contributions involving integrals (vv|ov).
C------------------------------------------------------------------
C
               CALL CCPT_3VT(ETAAI,WORK(KSCRTI),WORK(KDSRHF),
     *                       WORK(KLAMDP),WORK(KEND3),LWRK3,ISYDIS)
C
               CALL CCPT_YTV(ETAAI,WORK(KYTMAT),WORK(KDSRHF),
     *                       WORK(KLAMDP),WORK(KEND3),LWRK3,IDEL,ISYMD)
C
C-------------------------------------------------------------------
C              Calculate integral batch with three occupied indices.
C-------------------------------------------------------------------
C
               CALL CC_INT3O(WORK(K3OINT),WORK(KDSRHF),WORK(KLAMDP),
     *                      ISYMOP,WORK(KLAMDP),WORK(KEND3),LWRK3,
     *                      IDEL,ISYMD,LUDUM,'DUMMY')
C
C------------------------------------------------------------------
C              Calculate contributions involving integrals (oo|ov).
C------------------------------------------------------------------
C
               CALL CCPT_3OT(ETAAI,WORK(KSCRTI),WORK(K3OINT),
     *                       ISYDIS)
C
               CALL CCPT_NXY(ETAAI,WORK(KXMAT),WORK(KYMAT),WORK(K3OINT),
     *                       WORK(KDSRHF),WORK(KLAMDP),WORK(KEND3),
     *                       LWRK3,IDEL,ISYMD)
C
               CALL CCPT_XTO(ETAAI,WORK(KXTMAT),WORK(K3OINT),
     *                       WORK(KLAMDP),WORK(KEND3),LWRK3,IDEL,ISYMD)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C------------------------
C     Recover work space.
C------------------------
C
      KEND1 = KENDS2
      LWRK1 = LWRKS2
C
C------------------------------------------
C     Transform AO Fock matrix to MO basis.
C------------------------------------------
C
      IHELP = 1
      CALL CC_FCKMO(WORK(KFOCK),WORK(KLAMDP),WORK(KLAMDP),
     *                 WORK(KEND1),LWRK1,IHELP,IHELP,IHELP)
C
C-------------------------------------------------------
C     Calculate contributions involving the Fock matrix.
C-------------------------------------------------------
C
      CALL CCPT_FCK(ETAAI,WORK(KFOCK),WORK(KXTMAT),WORK(KYTMAT),
     *              WORK(KEND1),LWRK1)
C
C---------------------------------
C     Write out result and timing.
C---------------------------------
C
      IF (IPRINT .GT. 20) THEN
C
         CALL AROUND('Eta-kappa-0 vector exiting MP2_KARI')
C
         DO 20 ISYM = 1,NSYM
C
            WRITE(LUPRI,*) ' '
            WRITE(LUPRI,444) 'Sub-symmetry block number:', ISYM
            WRITE(LUPRI,555) '--------------------------'
  444       FORMAT(3X,A26,2X,I1)
  555       FORMAT(3X,A25)
C
            KOFF = IT1AM(ISYM,ISYM) + 1
            CALL OUTPUT(ETAAI(KOFF),1,NVIR(ISYM),1,NRHF(ISYM),
     *                  NVIR(ISYM),NRHF(ISYM),1,LUPRI)
C
            IF ((NVIR(ISYM) .EQ. 0) .OR. (NRHF(ISYM) .EQ. 0)) THEN
               WRITE(LUPRI,*) 'This sub-symmetry is empty'
            ENDIF
C
  20     CONTINUE
      ENDIF
C
      IF (IPRINT .GT. 9) THEN
         ETAKAN = DDOT(NT1AMX,ETAAI,1,ETAAI,1)
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'Norm of Eta-kappa-0:', ETAKAN
      ENDIF
C
      TIMETO = SECOND() - TIMETO
C
      IF (IPRINT .GT. 3) THEN
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) ' CCPT2 Eta-0(kappa) calculation completed'
         WRITE(LUPRI,*) 'Total time used in MP2_KARI:', TIMETO
      ENDIF
C
      CALL QEXIT('MP2_KARI')
C
      RETURN
      END
C  /* Deck ccpt_fck */
      SUBROUTINE CCPT_FCK(ETAAI,FCKMO,XTMAT,YTMAT,WORK,LWORK)
C
C     Written by Asger Halkier 9/9 - 1996.
C
C     Version: 1.0
C
C     Purpose: To calculate the Fock matrix contributions to 
C              ETAAI(CCPT2).
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ETAAI(*), FCKMO(*), XTMAT(*), YTMAT(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      CALL QENTER('CCPT_FCK')
C
      IF (LWORK .LT. NT1AMX) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', NT1AMX
         CALL QUIT('Insufficient memory for allocation in CCPT_FCK')
      ENDIF
C
C-----------------------------------------------------------------
C     Copy out needed part of Fock matrix F(ka) and store as T1AM.
C-----------------------------------------------------------------
C
      DO 100 ISYMC = 1,NSYM
C
         ISYMK = MULD2H(ISYMC,ISYMOP)
C
         DO 110 K = 1,NRHF(ISYMK)
C
            DO 120 C = 1,NVIR(ISYMC)
C
               KOFF1 = IFCVIR(ISYMK,ISYMC) + NORB(ISYMK)*(C - 1) + K
               KOFF2 = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1) + C
C
               WORK(KOFF2) = FCKMO(KOFF1)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      DO 130 ISYMA = 1,NSYM
C
         ISYMI = MULD2H(ISYMA,ISYMOP)
         ISYMK = MULD2H(ISYMA,ISYMOP)
         ISYMC = ISYMK
C
C-------------------------------------
C        Calculate XTMAT contribution.
C-------------------------------------
C
         KOFF1 = IT1AM(ISYMA,ISYMK)  + 1
         KOFF2 = IMATIJ(ISYMK,ISYMI) + 1
         KOFF3 = IT1AM(ISYMA,ISYMI)  + 1
C
         NTOTA = MAX(NVIR(ISYMA),1)
         NTOTK = MAX(NRHF(ISYMK),1)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMK),-ONE,
     *              WORK(KOFF1),NTOTA,XTMAT(KOFF2),NTOTK,ONE,
     *              ETAAI(KOFF3),NTOTA)
C
C-------------------------------------
C        Calculate YTMAT contribution.
C-------------------------------------
C
         KOFF4 = IMATAB(ISYMA,ISYMC) + 1
         KOFF5 = IT1AM(ISYMC,ISYMI)  + 1
         KOFF6 = IT1AM(ISYMA,ISYMI)  + 1
C
         NTOTA = MAX(NVIR(ISYMA),1)
         NTOTC = MAX(NVIR(ISYMC),1)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMC),-ONE,
     *              YTMAT(KOFF4),NTOTA,WORK(KOFF5),NTOTC,ONE,
     *              ETAAI(KOFF6),NTOTA)
C
  130 CONTINUE
C
      CALL QEXIT('CCPT_FCK')
C
      RETURN
      END
C  /* Deck ccpt_3ot */
      SUBROUTINE CCPT_3OT(ETAAI,TSCR,X3OINT,ISYDIS)
C
C     Written by Asger Halkier 10/9 - 1996.
C
C     Version: 1.0
C
C     Purpose: To calculate the contributions to ETAAI(CCPT2)
C              originating from amplitudes directly contracted 
C              with integrals (oo|ov).
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ETAAI(*), TSCR(*), X3OINT(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      CALL QENTER('CCPT_3OT')
C
      ISYALK = ISYDIS
      ISYLIK = ISYDIS
C
      DO 100 ISYMK = 1,NSYM
C
         ISYMAL = MULD2H(ISYMK,ISYALK)
         ISYMLI = MULD2H(ISYMK,ISYLIK)
C
         DO 110 K = 1,NRHF(ISYMK)
C
            DO 120 ISYMA = 1,NSYM
C
               ISYMI = ISYMA
               ISYML = MULD2H(ISYMA,ISYMAL)
C
C-----------------------------------------
C              Calculate the contribution.
C-----------------------------------------
C
               KOFF1 = IT2BCD(ISYMAL,ISYMK) + NT1AM(ISYMAL)*(K - 1)
     *               + IT1AM(ISYMA,ISYML)   + 1
               KOFF2 = IMAIJK(ISYMLI,ISYMK) + NMATIJ(ISYMLI)*(K - 1)
     *               + IMATIJ(ISYML,ISYMI)  + 1
               KOFF3 = IT1AM(ISYMA,ISYMI)   + 1
C
               NTOTA = MAX(NVIR(ISYMA),1)
               NTOTL = MAX(NRHF(ISYML),1)
C
               CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYML),
     *                    -ONE,TSCR(KOFF1),NTOTA,X3OINT(KOFF2),NTOTL,
     *                    ONE,ETAAI(KOFF3),NTOTA)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CCPT_3OT')
C
      RETURN
      END
C  /* Deck ccpt_3vt */
      SUBROUTINE CCPT_3VT(ETAAI,TSCR,DSRHF,XLAMDP,WORK,LWORK,ISYDIS)
C
C     Written by Asger Halkier 10/9 - 1996.
C
C     Version: 1.0
C
C     Purpose: To calculate the contributions to ETAAI(CCPT2)
C              originating from amplitudes directly contracted 
C              with integrals (oo|ov).
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ETAAI(*), TSCR(*), DSRHF(*), XLAMDP(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      CALL QENTER('CCPT_3VT')
C
      DO 100 ISYMK = 1,NSYM
C
         ISALBE = MULD2H(ISYMK,ISYDIS)
         ISYMAD = MULD2H(ISYMK,ISYDIS)
         ISYMDI = MULD2H(ISYMK,ISYDIS)
C
C----------------------------------
C        Work space allocation one.
C----------------------------------
C
         KAOINT = 1
         KEND1  = KAOINT + N2BST(ISALBE)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient memory for first allocation '//
     &                'in CCPT_3VT')
         ENDIF
C
         DO 110 K = 1,NRHF(ISYMK)
C
C----------------------------------------
C           Unpack integral distribution.
C----------------------------------------
C
            KOFF1 = IDSRHF(ISALBE,ISYMK) + NNBST(ISALBE)*(K - 1) + 1
C
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KAOINT))
C
            DO 120 ISYMA = 1,NSYM
C
               ISYMAL = ISYMA
               ISYMI  = ISYMA
               ISYMD  = MULD2H(ISYMA,ISYMAD)
               ISYMBE = ISYMD
C
C----------------------------------------
C              Work space allocation two.
C----------------------------------------
C
               KSCRAO = KEND1
               KSCRMO = KSCRAO + NBAS(ISYMAL)*NVIR(ISYMD)
               KEND2  = KSCRMO + NVIR(ISYMA)*NVIR(ISYMD)
               LWRK2  = LWORK  - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
                  CALL QUIT('Insufficient memory for allocation '//
     &                      'in CCPT_B3VT')
               ENDIF
C
               CALL DZERO(WORK(KSCRAO),NBAS(ISYMAL)*NVIR(ISYMD))
               CALL DZERO(WORK(KSCRMO),NVIR(ISYMA)*NVIR(ISYMD))
C
C--------------------------------------------------------------
C              Perform the three contractions to obtain result.
C--------------------------------------------------------------
C
               KOFF2  = KAOINT + IAODIS(ISYMAL,ISYMBE)
               KOFF3  = ILMVIR(ISYMD) + 1
C
               NTOTAL = MAX(NBAS(ISYMAL),1)
               NTOTBE = MAX(NBAS(ISYMBE),1)
C
               CALL DGEMM('N','N',NBAS(ISYMAL),NVIR(ISYMD),NBAS(ISYMBE),
     *                    ONE,WORK(KOFF2),NTOTAL,XLAMDP(KOFF3),NTOTBE,
     *                    ZERO,WORK(KSCRAO),NTOTAL)
C
               KOFF4  = ILMVIR(ISYMA) + 1
C
               NTOTAL = MAX(NBAS(ISYMAL),1)
               NTOTA  = MAX(NVIR(ISYMA),1)
C
               CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMD),NBAS(ISYMAL),
     *                    ONE,XLAMDP(KOFF4),NTOTAL,WORK(KSCRAO),NTOTAL,
     *                    ZERO,WORK(KSCRMO),NTOTA)
C
               KOFF5 = IT2BCD(ISYMDI,ISYMK) + NT1AM(ISYMDI)*(K - 1)
     *               + IT1AM(ISYMD,ISYMI)   + 1
               KOFF6 = IT1AM(ISYMA,ISYMI) + 1 
C
               NTOTA = MAX(NVIR(ISYMA),1)
               NTOTD = MAX(NVIR(ISYMD),1)
C
               CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMD),
     *                    ONE,WORK(KSCRMO),NTOTA,TSCR(KOFF5),NTOTD,
     *                    ONE,ETAAI(KOFF6),NTOTA)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CCPT_3VT')
C
      RETURN
      END
C  /* Deck ccpt_nxy */
      SUBROUTINE CCPT_NXY(ETAAI,XMAT,YMAT,X3OINT,DSRHF,XLAMDP,WORK,
     *                    LWORK,IDEL,ISYDEL)
C
C     Written by Asger Halkier 10/9 - 1996.
C
C     Version: 1.0
C
C     Purpose: To calculate the contributions to ETAAI(CCPT2)
C              containing the original (i.e. nonsymmetrized) 
C              X- and Y-matrices.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, FOUR = 4.0D0)
      DIMENSION ETAAI(*), XMAT(*), YMAT(*), X3OINT(*), DSRHF(*)
      DIMENSION XLAMDP(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      CALL QENTER('CCPT_NXY')
C
      ISYMA = ISYDEL
      ISYMI = ISYMA
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      KAVEC = 1
      KEND1 = KAVEC + NVIR(ISYMA)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation in '//
     &             'CCPT_NXY')
      ENDIF
C
      CALL DZERO(WORK(KAVEC),NVIR(ISYMA))
C
C-------------------------------------
C     Copy vector out of lambda matrix.
C-------------------------------------
C
      KOFF1 = ILMVIR(ISYMA) + IDEL - IBAS(ISYDEL)
C
      CALL DCOPY(NVIR(ISYMA),XLAMDP(KOFF1),NBAS(ISYDEL),WORK(KAVEC),1)
C
C----------------------------------------------
C     X- and Y- matrices are totally symmetric.
C----------------------------------------------
C
      ISYMKL = 1
      ISYMCD = 1
      ISALBE = ISYMCD
C
      DO 100 I = 1,NRHF(ISYMI)
C
C-----------------------------------------
C        Calculate contribution from XMAT.
C-----------------------------------------
C
         KOFF2 = IMAIJK(ISYMKL,ISYMI) + NMATIJ(ISYMKL)*(I - 1) + 1
         KOFF3 = IT1AM(ISYMA,ISYMI)   + NVIR(ISYMA)*(I - 1)    + 1
C
         FACT  = DDOT(NMATIJ(ISYMKL),XMAT,1,X3OINT(KOFF2),1)
C
         CALL DAXPY(NVIR(ISYMA),-FOUR*FACT,WORK(KAVEC),1,ETAAI(KOFF3),1)
C
C----------------------------------
C        Work space allocation two.
C----------------------------------
C
         KAOINT = KEND1
         KMOINT = KAOINT + N2BST(ISALBE)
         KEND2  = KMOINT + NMATAB(ISYMCD)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient memory for allocation in CCPT_NXY')
         ENDIF
C
         CALL DZERO(WORK(KMOINT),NMATAB(ISYMCD))
C
C-------------------------------------
C        Unpack integral distribution.
C-------------------------------------
C
         KOFF4 = IDSRHF(ISALBE,ISYMI) + NNBST(ISALBE)*(I - 1) + 1
C
         CALL CCSD_SYMSQ(DSRHF(KOFF4),ISALBE,WORK(KAOINT))
C
         DO 110 ISYMD = 1,NSYM
C
            ISYMAL = ISYMD
            ISYMC  = MULD2H(ISYMD,ISYMCD)
            ISYMBE = ISYMC
C
C---------------------------------------
C           Work space allocation three.
C---------------------------------------
C
            KSCRAO = KEND2
            KEND3  = KSCRAO + NBAS(ISYMAL)*NVIR(ISYMC)
            LWRK3  = LWORK  - KEND3
C
            IF (LWRK3 .LT. 0) THEN
               WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND3
               CALL QUIT('Insufficient memory for allocation in '//
     &                   'CCPT_NXY')
            ENDIF
C
            CALL DZERO(WORK(KSCRAO),NBAS(ISYMAL)*NVIR(ISYMC))
C
C-------------------------------------------
C           Transform integrals to MO basis.
C-------------------------------------------
C
            KOFF5  = KAOINT + IAODIS(ISYMAL,ISYMBE)
            KOFF6  = ILMVIR(ISYMC) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTBE = MAX(NBAS(ISYMBE),1)
C
            CALL DGEMM('N','N',NBAS(ISYMAL),NVIR(ISYMC),NBAS(ISYMBE),
     *                 ONE,WORK(KOFF5),NTOTAL,XLAMDP(KOFF6),NTOTBE,
     *                 ONE,WORK(KSCRAO),NTOTAL)
C
            KOFF7  = ILMVIR(ISYMD) + 1
            KOFF8  = KMOINT + IMATAB(ISYMD,ISYMC)
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTD  = MAX(NVIR(ISYMD),1)
C
            CALL DGEMM('T','N',NVIR(ISYMD),NVIR(ISYMC),NBAS(ISYMAL),
     *                 ONE,XLAMDP(KOFF7),NTOTAL,WORK(KSCRAO),NTOTAL,
     *                 ONE,WORK(KOFF8),NTOTD)
C
  110    CONTINUE
C
C------------------------------------------
C        Calculate contributions from YMAT.
C------------------------------------------
C
         FACT  = DDOT(NMATAB(ISYMCD),YMAT,1,WORK(KMOINT),1)
C
         KOFF9 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
C
         CALL DAXPY(NVIR(ISYMA),FOUR*FACT,WORK(KAVEC),1,ETAAI(KOFF9),1)
C
  100 CONTINUE
C
      CALL QEXIT('CCPT_NXY')
C
      RETURN
      END
C  /* Deck ccpt_xto */
      SUBROUTINE CCPT_XTO(ETAAI,XTMAT,X3OINT,XLAMDP,WORK,
     *                    LWORK,IDEL,ISYMD)
C
C     Written by Asger Halkier 10/9 - 1996.
C
C     Version: 1.0
C
C     Purpose: To calculate the contribution to ETAAI(CCPT2)
C              involving the symmetrized X-matrix (XTMAT) and the
C              (oo|ov) integrals.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, FOUR = 4.0D0)
      DIMENSION ETAAI(*), XTMAT(*), X3OINT(*), XLAMDP(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      CALL QENTER('CCPT_XTO')
C
      ISYMA  = ISYMD
      ISYMI  = ISYMA
      ISYMKL = 1
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      KAVEC = 1
      KIVEC = KAVEC + NVIR(ISYMA)
      KEND1 = KIVEC + NRHF(ISYMI)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation '//
     &             'in CCPT_XTO')
      ENDIF
C
      CALL DZERO(WORK(KAVEC),NVIR(ISYMA))
      CALL DZERO(WORK(KIVEC),NRHF(ISYMI))
C
C-------------------------------------
C     Copy vector out of lambda matrix.
C-------------------------------------
C
      KOFF1 = ILMVIR(ISYMA) + IDEL - IBAS(ISYMD)
C
      CALL DCOPY(NVIR(ISYMA),XLAMDP(KOFF1),NBAS(ISYMD),WORK(KAVEC),1)
C
      DO 100 ISYML = 1,NSYM
C
         ISYMK  = MULD2H(ISYML,ISYMKL)
         ISYMIK = MULD2H(ISYMI,ISYMK)
C
         DO 110 L = 1,NRHF(ISYML)
C
C--------------------------------------------------------
C           Contract integrals with symmetrized X-matrix.
C--------------------------------------------------------
C
            KOFF2 = IMAIJK(ISYMIK,ISYML) + NMATIJ(ISYMIK)*(L - 1)
     *            + IMATIJ(ISYMI,ISYMK)  + 1
            KOFF3 = IMATIJ(ISYMK,ISYML)  + NRHF(ISYMK)*(L - 1) + 1
C
            NTOTI = MAX(NRHF(ISYMI),1)
C
            CALL DGEMV('N',NRHF(ISYMI),NRHF(ISYMK),ONE,X3OINT(KOFF2),
     *                 NTOTI,XTMAT(KOFF3),1,ONE,WORK(KIVEC),1)
C
  110    CONTINUE
  100 CONTINUE
C
C-----------------------------
C     Final storage in result.
C-----------------------------
C
      DO 120 I = 1,NRHF(ISYMI)
C
         KOFF4 = KIVEC + I - 1
         KOFF5 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
C
         CALL DAXPY(NVIR(ISYMA),WORK(KOFF4),WORK(KAVEC),1,
     *              ETAAI(KOFF5),1)
C
  120 CONTINUE
C
      CALL QEXIT('CCPT_XTO')
C
      RETURN
      END
C  /* Deck ccpt_ytv */
      SUBROUTINE CCPT_YTV(ETAAI,YTMAT,DSRHF,XLAMDP,WORK,
     *                    LWORK,IDEL,ISYDEL)
C
C     Written by Asger Halkier 10/9 - 1996.
C
C     Version: 1.0
C
C     Purpose: To calculate the contribution to ETAAI(CCPT2)
C              involving the symmetrized Y-matrix (YTMAT) and the
C              (vv|ov) integrals.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, FOUR = 4.0D0)
      DIMENSION ETAAI(*), YTMAT(*), DSRHF(*), XLAMDP(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      CALL QENTER('CCPT_YTV')
C
      ISYMC = ISYDEL
      ISYMD = ISYMC
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      KCVEC = 1
      KDVEC = KCVEC + NVIR(ISYMC)
      KEND1 = KDVEC + NVIR(ISYMD)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation in '//
     &             'CCPT_YTV')
      ENDIF
C
      CALL DZERO(WORK(KCVEC),NVIR(ISYMC))
      CALL DZERO(WORK(KDVEC),NVIR(ISYMD))
C
C-------------------------------------
C     Copy vector out of lambda matrix.
C-------------------------------------
C
      KOFF1 = ILMVIR(ISYMC) + IDEL - IBAS(ISYDEL)
C
      CALL DCOPY(NVIR(ISYMC),XLAMDP(KOFF1),NBAS(ISYDEL),WORK(KCVEC),1)
C
C----------------------------------------
C     Contract with symmetrized Y-matrix.
C----------------------------------------
C
      KOFF1 = IMATAB(ISYMD,ISYMC) + 1
C
      NTOTD = MAX(NVIR(ISYMD),1)
C
      CALL DGEMV('N',NVIR(ISYMD),NVIR(ISYMC),ONE,YTMAT(KOFF1),NTOTD,
     *           WORK(KCVEC),1,ZERO,WORK(KDVEC),1)
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMA  = ISYMI
         ISYMAL = ISYMA
         ISYMBE = ISYMD
         ISALBE = MULD2H(ISYMAL,ISYMBE)
C
C----------------------------------
C        Work space allocation two.
C----------------------------------
C
         KAOINT = KEND1
         KSCRAO = KAOINT + N2BST(ISALBE)
         KMOINT = KSCRAO + NBAS(ISYMAL)*NVIR(ISYMD)
         KEND2  = KMOINT + NVIR(ISYMA)*NVIR(ISYMD)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient work space for allocation '//
     &                'in CCPT_YTV')
         ENDIF
C
         CALL DZERO(WORK(KSCRAO),NBAS(ISYMAL)*NVIR(ISYMD))
         CALL DZERO(WORK(KMOINT),NVIR(ISYMA)*NVIR(ISYMD))
C
         DO 110 I = 1,NRHF(ISYMI)
C
C----------------------------------------
C           Unpack integral distribution.
C----------------------------------------
C
            KOFF2 = IDSRHF(ISALBE,ISYMI) + NNBST(ISALBE)*(I - 1) + 1
C
            CALL CCSD_SYMSQ(DSRHF(KOFF2),ISALBE,WORK(KAOINT))
C
C-------------------------------------------
C           Transform integrals to MO basis.
C-------------------------------------------
C
            KOFF3  = KAOINT + IAODIS(ISYMAL,ISYMBE)
            KOFF4  = ILMVIR(ISYMD) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTBE = MAX(NBAS(ISYMBE),1)
C
            CALL DGEMM('N','N',NBAS(ISYMAL),NVIR(ISYMD),NBAS(ISYMBE),
     *                 ONE,WORK(KOFF3),NTOTAL,XLAMDP(KOFF4),NTOTBE,
     *                 ZERO,WORK(KSCRAO),NTOTAL)
C
            KOFF5  = ILMVIR(ISYMA) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTA  = MAX(NVIR(ISYMA),1)
C
            CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMD),NBAS(ISYMAL),
     *                 ONE,XLAMDP(KOFF5),NTOTAL,WORK(KSCRAO),NTOTAL,
     *                 ZERO,WORK(KMOINT),NTOTA)
C
            KOFF6 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
C
            NTOTA = MAX(NVIR(ISYMA),1)
C
            CALL DGEMV('N',NVIR(ISYMA),NVIR(ISYMD),-ONE,WORK(KMOINT),
     *                 NTOTA,WORK(KDVEC),1,ONE,ETAAI(KOFF6),1)
C
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CCPT_YTV')
C
      RETURN
      END
C  /* Deck cc_dedian */
      SUBROUTINE CC_DEDIAN(DENSI,MODEL,WORK,LWORK)
C
C     Written by Asger Halkier 18/3 - 1998
C
C     Version: 1.0
C
C     Purpose: To diagonalize and analyse the correlated
C              one-electron density matrix.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION DENSI(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
#include <ccsdinp.h>
C
      CHARACTER MODEL*4
C
      CALL QENTER('CC_DEDIAN')
C
C---------------------------
C     Work space allocation.
C---------------------------
C
      KNATOC = 1
      KIMANO = KNATOC + NORBT
      KIV1   = KIMANO + NORBT
      KFV1   = KIV1   + NORBT
      KEND1  = KFV1   + NORBT
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for allocation in CC_DEDIAN')
      ENDIF
C
C------------------------------------------------
C     Diagonalize the density in symmetry blocks.
C------------------------------------------------
C
      KOFF1  = 1
      KOFF2  = KNATOC
      KOFF3  = KIMANO
C
      CALL AROUND(MODEL//' Natural Occupations')
C
      DO 100 ISYM = 1,NSYM
C
         CALL DZERO(WORK(KIV1),NORBT)
         CALL DZERO(WORK(KFV1),NORBT)
C
         MATZ  = 0
C
         CALL RG(NORB(ISYM),NORB(ISYM),DENSI(KOFF1),WORK(KOFF2),
     *           WORK(KOFF3),MATZ,DUMMY,WORK(KIV1),WORK(KFV1),IERR)
C
         IF (IERR .NE. 0) THEN
            WRITE(LUPRI,*) 'RG returned non-zero status of IERR'
            WRITE(LUPRI,*) 'Diagonalization of one electron '//
     &           'density failed'
         ENDIF
C
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,444) 'Symmetry block number:', ISYM
         WRITE(LUPRI,555) '---------------------'
         WRITE(LUPRI,*) ' '
         IF (NORB(ISYM) .EQ. 0) THEN
            WRITE(LUPRI,777) 'No orbitals in this symmetry block'
         ELSE
            CALL SORTASH(WORK(KOFF2),WORK(KOFF3),NORB(ISYM))
            WRITE(LUPRI,666) (WORK(KOFF2 + I - 1), I = NORB(ISYM),1,-1)
C
            SUMSYM = ZERO
C
            DO 110 I = 1,NORB(ISYM)
C
               SUMSYM = SUMSYM + WORK(KOFF2 + I - 1)
C
  110       CONTINUE
C
            WRITE(LUPRI,*) ' '
            WRITE(LUPRI,888) 'Sum in this symmetry class:', SUMSYM
C
         ENDIF
C
         IF (IPRINT .GT. 50) THEN
C
            WRITE(LUPRI,*) ' '
            WRITE(LUPRI,555) 'Natocc imaginary part'
            WRITE(LUPRI,*) ' '
            WRITE(LUPRI,666) (WORK(KOFF3 + I - 1), I = NORB(ISYM),1,-1)
C
         ENDIF
C
  444    FORMAT(3X,A22,2X,I1)
  555    FORMAT(3X,A21)
  666    FORMAT(5F13.8)
  777    FORMAT(3X,A34)
  888    FORMAT(3X,A27,2X,F9.6)
C
         KOFF1 = KOFF1 + NORB(ISYM)*NORB(ISYM)
         KOFF2 = KOFF2 + NORB(ISYM)
         KOFF3 = KOFF3 + NORB(ISYM)
C
  100 CONTINUE
C
      CALL SORTASH(WORK(KNATOC),WORK(KIMANO),NORBT)
C
      CALL CCNAOCAN(WORK(KNATOC),WORK(KIMANO))
C
      CALL QEXIT('CC_DEDIAN')
C
      RETURN
      END
C  /* Deck mp_zkdia */
      SUBROUTINE MP2_ZKDIA(MODEL,ZKDIA,WORK,LWORK)
C
C     Written by Asger Halkier 20/3 - 1998
C
C     Version: 1.0
C
C     Purpose: To calculate the pp, ab, & ij parts of kappa-bar-0
C              that do not need the solution of any coupled equations.
C              ZKDIA holds all the blocks pq in the following order:
C              ij, ab, ai, ia; and these are stored full blocks after
C              each other. After these, the blocks containing frozen
C              core indices come: first cJ and then kJ.
C
#include <implicit.h>
#include <dummy.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      CHARACTER MODEL*10
      DIMENSION ZKDIA(*), WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccorb.h>
#include <iratdef.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <ccinftap.h>
#include <ccfro.h>
C
      CALL QENTER('MP2_ZKDIA')
C
      TIMETO = SECOND()
C
      IF (IPRINT .GT. 3) THEN
         CALL HEADER('Calculating diagonal blocks of zeta-kappa-0',-1)
      ENDIF
C
C------------------------------------------------------------------
C     Both t-vectors and tbar-vectors (zeta) are totally symmetric.
C------------------------------------------------------------------
C
      ISYMTR = 1
      ISYMOP = 1
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      KT2AM  = 1
      KXMAT  = KT2AM  + NT2AMX
      KYMAT  = KXMAT  + NMATIJ(1)
      KZ2AM  = KYMAT  + NMATAB(1)
      KT1AM  = KZ2AM  + NT2SQ(1)
      KZ1AM  = KT1AM  + NT1AMX
      KEND1  = KZ1AM  + NT1AMX
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for initial allocation '//
     &             'in MP2_ZKDIA')
      ENDIF
C
C----------------------------------------
C     Read zero'th order zeta amplitudes.
C----------------------------------------
C
      IOPT   = 3
      CALL CC_RDRSP('L0',0,1,IOPT,MODEL,WORK(KZ1AM),WORK(KZ2AM))
C
      KEND1 = KZ1AM
      LWRK1 = LWORK  - KEND1
C
C--------------------------------
C     Square up zeta2 amplitudes.
C--------------------------------
C
      CALL DCOPY(NT2AMX,WORK(KZ2AM),1,WORK(KT2AM),1)
      CALL CC_T2SQ(WORK(KT2AM),WORK(KZ2AM),1)
C
C
C-------------------------------------------
C     Read zero'th order cluster amplitudes.
C-------------------------------------------
C
      IOPT = 3
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),WORK(KT2AM))
C
      KEND1 = KT1AM
C
C
C--------------------------------------------------------
C     Calculate X-intermediate of tbar- and t-amplitudes.
C--------------------------------------------------------
C
      CALL CC_XI(WORK(KXMAT),WORK(KZ2AM),ISYMTR,WORK(KT2AM),ISYMOP,
     *             WORK(KEND1),LWRK1)
C
C--------------------------------------------------------
C     Calculate Y-intermediate of tbar- and t-amplitudes.
C--------------------------------------------------------
C
      CALL CC_YI(WORK(KYMAT),WORK(KZ2AM),ISYMTR,WORK(KT2AM),ISYMOP,
     *           WORK(KEND1),LWRK1)
C
C--------------------------------------------------------------------------
C     Calculate the diagonal elements ZK0(ii) = -X(ii) and ZK0(aa) = Y(aa).
C--------------------------------------------------------------------------
C
      DO 100 ISYMI = 1,NSYM
         DO 110 I = 1,NRHF(ISYMI)
C
            NII = IMATIJ(ISYMI,ISYMI) + NRHF(ISYMI)*(I - 1) + I
C
            ZKDIA(NII) = -WORK(KXMAT + NII - 1)
C
  110    CONTINUE
  100 CONTINUE
C
      DO 120 ISYMA = 1,NSYM
         DO 130 A = 1,NVIR(ISYMA)
C
            NAA = IMATAB(ISYMA,ISYMA) + NVIR(ISYMA)*(A - 1) + A
C
            ZKDIA(NMATIJ(1) + NAA) = WORK(KYMAT + NAA - 1)
C
  130    CONTINUE
  120 CONTINUE
C
C---------------------------------------
C     Set up 2C-E of cluster amplitudes.
C---------------------------------------
C
      ISYOPE = 1
      IOPTTCME = 1
      CALL CCSD_TCMEPK(WORK(KT2AM),1.0D0,ISYOPE,IOPTTCME)
C
C-------------------------------------------------------------
C     Set up special modified amplitudes T(2c-e) + Tbar.
C     Store it squared in KZ2AM to make smart contraction with
C     packed integrals (ai|bj) using the X- and Y-routines.
C-------------------------------------------------------------
C
      IOPT   = 3
      CALL CC_RDRSP('L0',0,1,IOPT,MODEL,WORK(KZ1AM),WORK(KZ2AM))
C
      CALL DSCAL(NT2AMX,TWO,WORK(KT2AM),1)
      CALL DAXPY(NT2AMX,ONE,WORK(KZ2AM),1,WORK(KT2AM),1)
      CALL CC_T2SQ(WORK(KT2AM),WORK(KZ2AM),1)
C
C-----------------------------------------------
C     Read integrals (ai|bj) = (ia|jb) from disc
C     (file always assumed open) into KT2AM.
C-----------------------------------------------
C
      REWIND(LUIAJB)
      READ(LUIAJB) (WORK(KT2AM + I - 1), I = 1,NT2AMX)
C
C-----------------------------------------------
C     Calculate modified X- and Y-intermediates.
C-----------------------------------------------
C
      CALL CC_XI(WORK(KXMAT),WORK(KZ2AM),ISYMTR,WORK(KT2AM),ISYMOP,
     *             WORK(KEND1),LWRK1)
C
      CALL CC_YI(WORK(KYMAT),WORK(KZ2AM),ISYMTR,WORK(KT2AM),ISYMOP,
     *           WORK(KEND1),LWRK1)
C
C---------------------------------------------
C     Calculate the ZK0(ab) and ZK0(ij) blocks
C     from modified X- and Y-intermediates.
C---------------------------------------------
C
      CALL MP2_ZKBLO(ZKDIA,WORK(KXMAT),WORK(KYMAT),WORK(KEND1),LWRK1)
C
C--------------------------------------------------- 
C     Calculate frozen core occupied blocks ZK0(iJ).
C--------------------------------------------------- 
C
      IF (FROIMP) THEN
         KOFRES = NMATIJ(1) + NMATAB(1) + 2*NT1AMX + 2*NT1FRO(1) + 1
         CALL MP2_ZKFCB(ZKDIA(KOFRES),WORK(KZ2AM),WORK(KEND1),LWRK1)
      ENDIF
C
C------------------------------------------------
C     Write out timings and results if requested.
C------------------------------------------------
C
      IF (IPRINT .GT. 3) THEN
         CALL AROUND('Zeta-kappa-0 diagonal blocks')
         ZKAPI1 = DDOT(NMATIJ(1),ZKDIA(1),1,ZKDIA(1),1)
         ZKAPA1 = DDOT(NMATAB(1),ZKDIA(NMATIJ(1)+1),1,
     *                 ZKDIA(NMATIJ(1)+1),1)
         ZKAPIJ = ZKAPI1**0.5
         ZKAPAB = ZKAPA1**0.5
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'Norm of occupied-occupied block:', ZKAPIJ
         WRITE(LUPRI,*) 'Norm of virtual-virtual block:', ZKAPAB
         IF (FROIMP) THEN
            ZKAPF1 = DDOT(NCOFRO(1),ZKDIA(KOFRES),1,
     *                    ZKDIA(KOFRES),1)
            ZKAPFR = ZKAPF1**0.5
         WRITE(LUPRI,*) 'Norm of frozen-core-occupied block:', ZKAPFR
         ENDIF
C
         IF (IPRINT .GT. 50) THEN
            DO 140 ISYM = 1,NSYM
               WRITE(LUPRI,*) ' '
               WRITE(LUPRI,*) 'Symmetry block:', ISYM
               KIJ = IMATIJ(ISYM,ISYM) + 1
               KAB = IMATAB(ISYM,ISYM) + 1 + NMATIJ(1)
               CALL AROUND('occ-occ block')
               CALL OUTPUT(ZKDIA(KIJ),1,NRHF(ISYM),1,NRHF(ISYM),
     *                     NRHF(ISYM),NRHF(ISYM),1,LUPRI)
               CALL AROUND('vir-vir block')
               CALL OUTPUT(ZKDIA(KAB),1,NVIR(ISYM),1,NVIR(ISYM),
     *                     NVIR(ISYM),NVIR(ISYM),1,LUPRI)
  140       CONTINUE
         ENDIF
      ENDIF
C
      TIMETO = SECOND() - TIMETO
C
      IF (IPRINT .GT. 3) THEN
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'Diagonal blocks of Zeta-kappa-0 calculated'
         WRITE(LUPRI,*) 'Total time used in MP2_ZKDIA:', TIMETO
      ENDIF
C
      CALL QEXIT('MP2_ZKDIA')
      RETURN
      END
C  /* Deck mp_zkblo */
      SUBROUTINE MP2_ZKBLO(ZKDIA,XMAT,YMAT,WORK,LWORK)
C
C     Written by Asger Halkier 22/3 - 1998
C
C     Version: 1.0
C
C     Purpose: To calculate the ab & ij parts of kappa-bar-0,
C              from modified X- and Y-intermediates (XMAT & YMAT)
C              and canonical orbital energies.
C
C     If degeneracies occur among the orbitals, the divergent terms
C     with the corresponding orbital energy difference denominators
C     are skipped. This is controlled via the THRDEM parameter.
C
C     Small modifications for CC2 by A. Halkier & S. Coriani
C     14/01-2000. Introduce factor FACT to control antisymmetrization
C     of eta_ij and eta_ab.
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THRDEM = 1.0D-12)
      DIMENSION ZKDIA(*), XMAT(*), YMAT(*), WORK(LWORK)
#include <maxorb.h>
#include <ccorb.h>
#include <iratdef.h>
#include <inftap.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
C
      CALL QENTER('MP2_ZKBLO')
C
      IF (MP2) FACT = ONE
      IF (CC2) FACT = ZERO
C
C---------------------------
C     Work space allocation.
C---------------------------
C
      KFOCKD = 1
      KEND1  = KFOCKD + NORBTS
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need:', KEND1, 'Available:', LWORK
         CALL QUIT('Insufficient memory for allocation in MP2_ZKBLO')
      ENDIF
C
      CALL DZERO(WORK(KFOCKD),NORBTS)
C
C-------------------------------------
C     Read canonical orbital energies.
C-------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND (LUSIFC)
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD + I - 1), I = 1,NORBTS)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C----------------------------------------------------------------
C     Change symmetry ordering of the canonical orbital energies.
C----------------------------------------------------------------
C
      IF (FROIMP)
     *    CALL CCSD_DELFRO(WORK(KFOCKD),WORK(KEND1),LWRK1)
C
      CALL FOCK_REORDER(WORK(KFOCKD),WORK(KEND1),LWRK1)
C
C---------------------------
C     Calculate the results:
C     Occupied block:
C---------------------------
C
      DO 100 ISYMI = 1,NSYM
         ISYMJ = ISYMI
         DO 110 J = 1,NRHF(ISYMJ)
            KOFFJ = KFOCKD + IRHF(ISYMJ) + J - 1
            DO 120 I = J+1,NRHF(ISYMI)
               KOFFI = KFOCKD + IRHF(ISYMI) + I - 1
               DENOM = WORK(KOFFJ) - WORK(KOFFI)
               IF (ABS(DENOM) .GT. THRDEM) THEN
                 NIJ   = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + I
                 NJI   = IMATIJ(ISYMJ,ISYMI) + NRHF(ISYMJ)*(I - 1) + J
C
                 ZKDIA(NIJ) = HALF*(XMAT(NIJ) - FACT*XMAT(NJI))/DENOM
                 ZKDIA(NJI) = ZKDIA(NIJ)
               ENDIF
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C-------------------
C     Virtual block:
C-------------------
C
      DO 130 ISYMA = 1,NSYM
         ISYMB = ISYMA
         DO 140 B = 1,NVIR(ISYMB)
            KOFFB = KFOCKD + IVIR(ISYMB) + B - 1
            DO 150 A = B+1,NVIR(ISYMA)
               KOFFA = KFOCKD + IVIR(ISYMA) + A - 1
               DENOM = WORK(KOFFB) - WORK(KOFFA)
               IF (ABS(DENOM) .GT. THRDEM) THEN
                 NAB   = IMATAB(ISYMA,ISYMB) + NVIR(ISYMA)*(B - 1) + A
                 NBA   = IMATAB(ISYMB,ISYMA) + NVIR(ISYMB)*(A - 1) + B
C
                 ZKDIA(NMATIJ(1) + NAB) =
     *                       HALF*(YMAT(NAB) - FACT*YMAT(NBA))/DENOM
                 ZKDIA(NMATIJ(1) + NBA) = ZKDIA(NMATIJ(1) + NAB)
               ENDIF
C
  150       CONTINUE
  140    CONTINUE
  130 CONTINUE
C
      CALL QEXIT('MP2_ZKBLO')
C
      RETURN
      END
C  /* Deck mp2_kanew */
      SUBROUTINE MP2_KANEW(MODEL,ETAAI,ZKDIA,WORK,LWORK)
C
C     Written by Asger Halkier 23/3 - 1998
C
C     Version: 1.0
C
C     Purpose: To calculate the right hand side ETAAI for the
C              equations for the zero'th order orbital rotation
C              multipliers in MP2 calculations.
C
C     Modifications for inclusion of frozen core contributions
C     by Asger Halkier 28/5 - 1998.
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <maxash.h>
#include <mxorb.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      CHARACTER MODEL*10
      DIMENSION INDEXA(MXCORB)
      DIMENSION ETAAI(*), ZKDIA(*), WORK(LWORK)
#include <ccorb.h>
#include <infind.h>
#include <blocks.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccinftap.h>
#include <ccsdio.h>
#include <distcl.h>
#include <cbieri.h>
#include <eribuf.h>
#include <cclr.h>
#include <ccfro.h>
C
      CALL QENTER('MP2_KANEW')
C
      CALL HEADER('Constructing right-hand-side for MP2-kappa-0(ai)',-1)
C
      TIMETO = ZERO
      TIMETO = SECOND()
C
C------------------------------------------------------------------
C     Both t-vectors and tbar-vectors (zeta) are totally symmetric.
C------------------------------------------------------------------
C
      ISYMTR = 1
      ISYMOP = 1
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C  
      KAFROI = 1
      KT2AM  = KAFROI + NT1FRO(1)
      KLAMDP = KT2AM  + NT2AMX
      KLAMDH = KLAMDP + NLAMDT
      KZ2AM  = KLAMDH + NLAMDT
      KT1AM  = KZ2AM  + NT2AMX
      KEND1  = KT1AM  + NT1AMX
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for initial allocation '//
     &             'in MP2_KANEW')
      ENDIF
C
      CALL DZERO(WORK(KAFROI),NT1FRO(1))
C
C-------------------------------------------
C     Read zero'th order cluster amplitudes.
C-------------------------------------------
C
      IOPT = 3
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),WORK(KT2AM))
C
      CALL DZERO(WORK(KT1AM),NT1AMX)
C
C----------------------------------
C     Calculate the lambda matrices.
C----------------------------------
C
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND1),
     *            LWRK1)
C
C----------------------------------------
C     Read zero'th order zeta amplitudes.
C----------------------------------------
C
      IOPT   = 3
      CALL CC_RDRSP('L0',0,1,IOPT,MODEL,WORK(KT1AM),WORK(KZ2AM))
C
C---------------------------------------
C     Set up 2C-E of cluster amplitudes.
C---------------------------------------
C
      ISYOPE = 1
      IOPTTCME = 1
      CALL CCSD_TCMEPK(WORK(KT2AM),1.0D0,ISYOPE,IOPTTCME)
C
C--------------------------------------------------------------------
C     Set up special modified amplitudes needed in the integral loop.
C     (By doing it this way, we only need one packed vector in core
C     along with the integral distribution in the delta loop.)
C--------------------------------------------------------------------
C
      CALL DSCAL(NT2AMX,TWO,WORK(KT2AM),1)
      CALL DAXPY(NT2AMX,ONE,WORK(KZ2AM),1,WORK(KT2AM),1)
C
      KEND1 = KLAMDH
      LWRK1 = LWORK - KEND1
C
C--------------------------------------------------------------------
C     Calculate the full MO coefficient matrix for frozen core calcs.
C--------------------------------------------------------------------
C
      IF (FROIMP) THEN
C
         KCMO  = KEND1
         KEND1 = KCMO  + NLAMDS
         LWKR1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient memory for allocation '//
     &                'in MP2_KANEW')
         ENDIF
C
         CALL CMO_ALL(WORK(KCMO),WORK(KEND1),LWRK1)
C
      ENDIF
C
C-----------------------------------
C     Start the loop over integrals.
C-----------------------------------
C
      KENDS2 = KEND1
      LWRKS2 = LWRK1
C
      IF (DIRECT) THEN
         IF (HERDIR) THEN
           CALL HERDI1(WORK(KEND1),LWRK1,IPRERI)
         ELSE
           KCCFB1 = KEND1
           KINDXB = KCCFB1 + MXPRIM*MXCONT
           KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
           LWRK1  = LWORK  - KEND1
           CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     *                 KODPP1,KODPP2,KRDPP1,KRDPP2,
     *                 KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     *                 WORK(KEND1),LWRK1,IPRERI)
           KEND1 = KFREE
           LWRK1 = LFREE
         END IF
         NTOSYM = 1
      ELSE
         NTOSYM = NSYM
      ENDIF
C
      KENDSV = KEND1
      LWRKSV = LWRK1
C
      ICDEL1 = 0
      DO 100 ISYMD1 = 1,NTOSYM
C
         IF (DIRECT) THEN
            IF (HERDIR) THEN
              NTOT = MAXSHL
            ELSE
              NTOT = MXCALL
            END IF
         ELSE
            NTOT = NBAS(ISYMD1)
         ENDIF
C
         DO 110 ILLL = 1,NTOT
C
C---------------------------------------------
C           If direct calculate the integrals.
C---------------------------------------------
C
            IF (DIRECT) THEN
C
               KEND1 = KENDSV
               LWRK1 = LWRKSV
C
               DTIME  = SECOND()
               IF (HERDIR) THEN
                 CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,
     &                       IPRERI)
               ELSE
                 CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     *                       WORK(KODCL1),WORK(KODCL2),
     *                       WORK(KODBC1),WORK(KODBC2),
     *                       WORK(KRDBC1),WORK(KRDBC2),
     *                       WORK(KODPP1),WORK(KODPP2),
     *                       WORK(KRDPP1),WORK(KRDPP2),
     *                       WORK(KCCFB1),WORK(KINDXB),
     *                       WORK(KEND1), LWRK1,IPRERI)
               END IF
               DTIME   = SECOND() - DTIME
               TIMHE2 = TIMHE2 + DTIME
C
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
                  CALL QUIT('Insufficient memory for integrals '//
     &                      'in MP2_KANEW')
               END IF
C
            ELSE
               NUMDIS = 1
            ENDIF
C
C-----------------------------------------------------
C           Loop over number of distributions in disk.
C-----------------------------------------------------
C
            DO 120 IDEL2 = 1,NUMDIS
C
               IF (DIRECT) THEN
                  IDEL  = INDEXA(IDEL2)
                  ISYMD = ISAO(IDEL)
               ELSE
                  IDEL  = IBAS(ISYMD1) + ILLL
                  ISYMD = ISYMD1
               ENDIF
C
C----------------------------------------
C              Work space allocation two.
C----------------------------------------
C
               ISYDIS = MULD2H(ISYMD,ISYMOP)
C
               KXINT  = KEND1
               KEND2  = KXINT + NDISAO(ISYDIS)
               LWRK2  = LWORK - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
                  CALL QUIT('Insufficient memory for integrals '//
     &                      'in MP2_KANEW')
               ENDIF
C
C--------------------------------------------
C              Read AO integral distribution.
C--------------------------------------------
C
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                     WORK(KRECNR),DIRECT)
C
C------------------------------------------
C              Work space allocation three.
C------------------------------------------
C
               KDSRHF = KEND2
               K3OINT = KDSRHF + NDSRHF(ISYMD)
               KSCRTI = K3OINT + NMAIJK(ISYDIS)
               IF (FROIMP) THEN
                  KDSFRO = KSCRTI + NT2BCD(ISYDIS)
                  KOFOIN = KDSFRO + NDSFRO(ISYDIS)
                  KEND3  = KOFOIN + NOFROO(ISYDIS) 
               ELSE
                  KEND3  = KSCRTI + NT2BCD(ISYDIS)
               ENDIF
               LWRK3  = LWORK  - KEND3
C
               IF (LWRK3 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND3,'Available : ',LWORK
                  CALL QUIT('Insufficient memory for integrals '//
     &                      'in MP2_KANEW')
               ENDIF
C
C---------------------------------------------------------------------
C              Calculate partially backtransformed modified amplitude.
C---------------------------------------------------------------------
C
               CALL CC_TI(WORK(KSCRTI),ISYMD,WORK(KT2AM),ISYMOP,
     *                    WORK(KLAMDP),1,WORK(KEND3),LWRK3,IDEL,ISYMD)
C
C----------------------------------------------------------------------
C              Transform one index in the integral batch to correlated.
C----------------------------------------------------------------------
C
               CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KLAMDP),ISYMOP,
     *                     WORK(KEND3),LWRK3,ISYDIS)
C
C------------------------------------------------------------------
C              Transform one index in the integral batch to frozen.
C------------------------------------------------------------------
C
               IF (FROIMP) THEN
C
                  CALL CC_GTOFRO(WORK(KXINT),WORK(KDSFRO),WORK(KCMO),
     *                           WORK(KEND3),LWRK3,ISYDIS)
C
C--------------------------------------------------------------
C                 Calculate integral batch (cor fro | cor del).
C--------------------------------------------------------------
C
                  CALL CC_OFROIN(WORK(KDSRHF),WORK(KOFOIN),WORK(KCMO),
     *                           WORK(KEND3),LWRK3,ISYDIS)
C
C---------------------------------------------------------------
C                 Calculate direct contribution to frozen block.
C---------------------------------------------------------------
C
                  CALL MP2_ETFRD(WORK(KAFROI),WORK(KOFOIN),
     *                           WORK(KSCRTI),ISYDIS)
C
C-------------------------------------------------------------------------
C                 Calculate indirect virtual contribution to frozen block.
C-------------------------------------------------------------------------
C
                  CALL MP2_EIDV1(WORK(KAFROI),WORK(KDSFRO),
     *                           ZKDIA(NMATIJ(1)+1),WORK(KCMO),
     *                           WORK(KEND3),LWRK3,IDEL,ISYMD)
C
                  CALL MP2_EIDV2(WORK(KAFROI),WORK(KDSFRO),
     *                           ZKDIA(NMATIJ(1)+1),WORK(KCMO),
     *                           WORK(KEND3),LWRK3,IDEL,ISYMD)
C
C----------------------------------------------------------------------------
C                 Calculate indirect correlated contribution to frozen block.
C----------------------------------------------------------------------------
C
                  CALL MP2_EIDC1(WORK(KAFROI),WORK(KDSFRO),
     *                           ZKDIA(1),WORK(KCMO),WORK(KEND3),
     *                           LWRK3,IDEL,ISYMD)
C
                  CALL MP2_EIDC2(WORK(KAFROI),WORK(KOFOIN),
     *                           ZKDIA(1),WORK(KCMO),WORK(KEND3),
     *                           LWRK3,IDEL,ISYMD)
C
C-----------------------------------------------------------------------------
C                 Calculate indirect frozen contribution to both parts of eta.
C-----------------------------------------------------------------------------
C
                  KOFFJK = NMATIJ(1)   + NMATAB(1) + 2*NT1AMX
     *                   + 2*NT1FRO(1) + 1
C
                  CALL MP2_EIDF1(ETAAI,WORK(KOFOIN),ZKDIA(KOFFJK),
     *                           WORK(KCMO),WORK(KEND3),LWRK3,
     *                           IDEL,ISYMD)
C
                  CALL MP2_EIDF2(WORK(KAFROI),WORK(KDSFRO),
     *                           ZKDIA(KOFFJK),WORK(KCMO),WORK(KEND3),
     *                           LWRK3,IDEL,ISYMD)
C
                  CALL MP2_EIDF3(ETAAI,WORK(KOFOIN),ZKDIA(KOFFJK),
     *                           WORK(KCMO),WORK(KEND3),LWRK3,
     *                           IDEL,ISYMD)
C
                  CALL MP2_EIDF4(ETAAI,WORK(KDSFRO),ZKDIA(KOFFJK),
     *                           WORK(KCMO),WORK(KEND3),LWRK3,
     *                           IDEL,ISYMD)
C
                  CALL MP2_EIDF5(WORK(KAFROI),WORK(KDSRHF),
     *                           ZKDIA(KOFFJK),WORK(KCMO),WORK(KEND3),
     *                           LWRK3,IDEL,ISYMD)
C
                  CALL MP2_EIDF6(WORK(KAFROI),WORK(KDSFRO),
     *                           ZKDIA(KOFFJK),WORK(KCMO),WORK(KEND3),
     *                           LWRK3,IDEL,ISYMD)
C
               ENDIF
C
C------------------------------------------------------------------
C              Calculate contributions involving integrals (vv|ov).
C------------------------------------------------------------------
C
               CALL CCPT_3VT(ETAAI,WORK(KSCRTI),WORK(KDSRHF),
     *                       WORK(KLAMDP),WORK(KEND3),LWRK3,ISYDIS)
C
               CALL MP2_YTV(ETAAI,ZKDIA(NMATIJ(1)+1),WORK(KDSRHF),
     *                      WORK(KLAMDP),WORK(KEND3),LWRK3,IDEL,ISYMD)
C
C-------------------------------------------------------------------
C              Calculate integral batch with three occupied indices.
C-------------------------------------------------------------------
C
               CALL CC_INT3O(WORK(K3OINT),WORK(KDSRHF),WORK(KLAMDP),
     *                      ISYMOP,WORK(KLAMDP),WORK(KEND3),LWRK3,
     *                      IDEL,ISYMD,LUDUM,'DUMMY')
C
C------------------------------------------------------------------
C              Calculate contributions involving integrals (oo|ov).
C------------------------------------------------------------------
C
               CALL CCPT_3OT(ETAAI,WORK(KSCRTI),WORK(K3OINT),
     *                       ISYDIS)
C
               CALL MP2_NXY(ETAAI,ZKDIA(1),ZKDIA(NMATIJ(1)+1),
     *                      WORK(K3OINT),WORK(KDSRHF),WORK(KLAMDP),
     *                      WORK(KEND3),LWRK3,IDEL,ISYMD)
C
               CALL MP2_XTO(ETAAI,ZKDIA(1),WORK(K3OINT),
     *                      WORK(KLAMDP),WORK(KEND3),LWRK3,IDEL,ISYMD)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C---------------------
C     Reorder results.
C---------------------
C
      CALL CC_ETARE(ETAAI,WORK(KAFROI),WORK(KENDS2),LWRKS2)
C
C---------------------------------
C     Write out result and timing.
C---------------------------------
C
      IF (IPRINT .GT. 20) THEN
C
         CALL AROUND('Eta-kappa-0 vector exiting MP2_KANEW')
C
         DO 20 ISYM = 1,NSYM
C
            WRITE(LUPRI,*) ' '
            WRITE(LUPRI,444) 'Sub-symmetry block number:', ISYM
            WRITE(LUPRI,555) '--------------------------'
  444       FORMAT(3X,A26,2X,I1)
  555       FORMAT(3X,A25)
C
            KOFF = IALLAI(ISYM,ISYM) + 1
            CALL OUTPUT(ETAAI(KOFF),1,NVIR(ISYM),1,NRHFS(ISYM),
     *                  NVIR(ISYM),NRHFS(ISYM),1,LUPRI)
C
            IF ((NVIR(ISYM) .EQ. 0) .OR. (NRHFS(ISYM) .EQ. 0)) THEN
               WRITE(LUPRI,*) 'This sub-symmetry is empty'
            ENDIF
C
  20     CONTINUE
      ENDIF
C
      IF (IPRINT .GT. 9) THEN
         ETAKAN = DDOT(NALLAI(1),ETAAI,1,ETAAI,1)
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'Norm of Eta-kappa-0:', ETAKAN
      ENDIF
C
      TIMETO = SECOND() - TIMETO
C
      IF (IPRINT .GT. 3) THEN
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'MP2 Eta-0(kappa) calculation completed'
         WRITE(LUPRI,*) 'Total time used in MP2_KANEW:', TIMETO
      ENDIF
C
      CALL QEXIT('MP2_KANEW')
C
      RETURN
      END
C  /* Deck mp2_nxy */
      SUBROUTINE MP2_NXY(ETAAI,XMAT,YMAT,X3OINT,DSRHF,XLAMDP,WORK,
     *                    LWORK,IDEL,ISYDEL)
C
C     Written by Asger Halkier 23/3 - 1998.
C
C     Version: 1.0
C
C     Purpose: To calculate the contributions to ETAAI(MP2)
C              originating from the coulomb part of the "extra
C              terms" from the diagonal orbital multipliers.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, FOUR = 4.0D0)
      DIMENSION ETAAI(*), XMAT(*), YMAT(*), X3OINT(*), DSRHF(*)
      DIMENSION XLAMDP(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      CALL QENTER('MP2_NXY')
C
      ISYMA = ISYDEL
      ISYMI = ISYMA
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      KAVEC = 1
      KEND1 = KAVEC + NVIR(ISYMA)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation in MP2_NXY')
      ENDIF
C
      CALL DZERO(WORK(KAVEC),NVIR(ISYMA))
C
C-------------------------------------
C     Copy vector out of lambda matrix.
C-------------------------------------
C
      KOFF1 = ILMVIR(ISYMA) + IDEL - IBAS(ISYDEL)
C
      CALL DCOPY(NVIR(ISYMA),XLAMDP(KOFF1),NBAS(ISYDEL),WORK(KAVEC),1)
C
C----------------------------------------------
C     X- and Y- matrices are totally symmetric.
C----------------------------------------------
C
      ISYMKL = 1
      ISYMCD = 1
      ISALBE = ISYMCD
C
      DO 100 I = 1,NRHF(ISYMI)
C
C-----------------------------------------
C        Calculate contribution from XMAT.
C-----------------------------------------
C
         KOFF2 = IMAIJK(ISYMKL,ISYMI) + NMATIJ(ISYMKL)*(I - 1) + 1
         KOFF3 = IT1AM(ISYMA,ISYMI)   + NVIR(ISYMA)*(I - 1)    + 1
C
         FACT  = DDOT(NMATIJ(ISYMKL),XMAT,1,X3OINT(KOFF2),1)
C
         CALL DAXPY(NVIR(ISYMA),FOUR*FACT,WORK(KAVEC),1,ETAAI(KOFF3),1)
C
C----------------------------------
C        Work space allocation two.
C----------------------------------
C
         KAOINT = KEND1
         KMOINT = KAOINT + N2BST(ISALBE)
         KEND2  = KMOINT + NMATAB(ISYMCD)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient memory for allocation in MP2_NXY')
         ENDIF
C
         CALL DZERO(WORK(KMOINT),NMATAB(ISYMCD))
C
C-------------------------------------
C        Unpack integral distribution.
C-------------------------------------
C
         KOFF4 = IDSRHF(ISALBE,ISYMI) + NNBST(ISALBE)*(I - 1) + 1
C
         CALL CCSD_SYMSQ(DSRHF(KOFF4),ISALBE,WORK(KAOINT))
C
         DO 110 ISYMD = 1,NSYM
C
            ISYMAL = ISYMD
            ISYMC  = MULD2H(ISYMD,ISYMCD)
            ISYMBE = ISYMC
C
C---------------------------------------
C           Work space allocation three.
C---------------------------------------
C
            KSCRAO = KEND2
            KEND3  = KSCRAO + NBAS(ISYMAL)*NVIR(ISYMC)
            LWRK3  = LWORK  - KEND3
C
            IF (LWRK3 .LT. 0) THEN
               WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND3
               CALL QUIT('Insufficient memory for allocation '//
     &                   'in MP2_NXY')
            ENDIF
C
            CALL DZERO(WORK(KSCRAO),NBAS(ISYMAL)*NVIR(ISYMC))
C
C-------------------------------------------
C           Transform integrals to MO basis.
C-------------------------------------------
C
            KOFF5  = KAOINT + IAODIS(ISYMAL,ISYMBE)
            KOFF6  = ILMVIR(ISYMC) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTBE = MAX(NBAS(ISYMBE),1)
C
            CALL DGEMM('N','N',NBAS(ISYMAL),NVIR(ISYMC),NBAS(ISYMBE),
     *                 ONE,WORK(KOFF5),NTOTAL,XLAMDP(KOFF6),NTOTBE,
     *                 ONE,WORK(KSCRAO),NTOTAL)
C
            KOFF7  = ILMVIR(ISYMD) + 1
            KOFF8  = KMOINT + IMATAB(ISYMD,ISYMC)
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTD  = MAX(NVIR(ISYMD),1)
C
            CALL DGEMM('T','N',NVIR(ISYMD),NVIR(ISYMC),NBAS(ISYMAL),
     *                 ONE,XLAMDP(KOFF7),NTOTAL,WORK(KSCRAO),NTOTAL,
     *                 ONE,WORK(KOFF8),NTOTD)
C
  110    CONTINUE
C
C------------------------------------------
C        Calculate contributions from YMAT.
C------------------------------------------
C
         FACT  = DDOT(NMATAB(ISYMCD),YMAT,1,WORK(KMOINT),1)
C
         KOFF9 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
C
         CALL DAXPY(NVIR(ISYMA),FOUR*FACT,WORK(KAVEC),1,ETAAI(KOFF9),1)
C
  100 CONTINUE
C
      CALL QEXIT('MP2_NXY')
C
      RETURN
      END
C  /* Deck mp2_xto */
      SUBROUTINE MP2_XTO(ETAAI,XTMAT,X3OINT,XLAMDP,WORK,
     *                    LWORK,IDEL,ISYMD)
C
C     Written by Asger Halkier 23/3 - 1998.
C
C     Version: 1.0
C
C     Purpose: To calculate the (oo|ov) contributions to ETAAI(MP2)
C              originating from the exchange part of the "extra
C              terms" from the diagonal orbital multipliers.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION ETAAI(*), XTMAT(*), X3OINT(*), XLAMDP(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      CALL QENTER('MP2_XTO')
C
      ISYMA  = ISYMD
      ISYMI  = ISYMA
      ISYMKL = 1
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      KAVEC = 1
      KIVEC = KAVEC + NVIR(ISYMA)
      KEND1 = KIVEC + NRHF(ISYMI)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation in MP2_XTO')
      ENDIF
C
      CALL DZERO(WORK(KAVEC),NVIR(ISYMA))
      CALL DZERO(WORK(KIVEC),NRHF(ISYMI))
C
C-------------------------------------
C     Copy vector out of lambda matrix.
C-------------------------------------
C
      KOFF1 = ILMVIR(ISYMA) + IDEL - IBAS(ISYMD)
C
      CALL DCOPY(NVIR(ISYMA),XLAMDP(KOFF1),NBAS(ISYMD),WORK(KAVEC),1)
C
      DO 100 ISYML = 1,NSYM
C
         ISYMK  = MULD2H(ISYML,ISYMKL)
         ISYMIK = MULD2H(ISYMI,ISYMK)
C
         DO 110 L = 1,NRHF(ISYML)
C
C--------------------------------------------------------
C           Contract integrals with symmetrized X-matrix.
C--------------------------------------------------------
C
            KOFF2 = IMAIJK(ISYMIK,ISYML) + NMATIJ(ISYMIK)*(L - 1)
     *            + IMATIJ(ISYMI,ISYMK)  + 1
            KOFF3 = IMATIJ(ISYMK,ISYML)  + NRHF(ISYMK)*(L - 1) + 1
C
            NTOTI = MAX(NRHF(ISYMI),1)
C
            CALL DGEMV('N',NRHF(ISYMI),NRHF(ISYMK),ONE,X3OINT(KOFF2),
     *                 NTOTI,XTMAT(KOFF3),1,ONE,WORK(KIVEC),1)
C
  110    CONTINUE
  100 CONTINUE
C
C-----------------------------
C     Final storage in result.
C-----------------------------
C
      DO 120 I = 1,NRHF(ISYMI)
C
         KOFF4 = KIVEC + I - 1
         KOFF5 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
C
         CALL DAXPY(NVIR(ISYMA),-TWO*WORK(KOFF4),WORK(KAVEC),1,
     *              ETAAI(KOFF5),1)
C
  120 CONTINUE
C
      CALL QEXIT('MP2_XTO')
C
      RETURN
      END
C  /* Deck mp2_ytv */
      SUBROUTINE MP2_YTV(ETAAI,YTMAT,DSRHF,XLAMDP,WORK,
     *                    LWORK,IDEL,ISYDEL)
C
C     Written by Asger Halkier 23/3 - 1998.
C
C     Version: 1.0
C
C     Purpose: To calculate the (vv|ov) contributions to ETAAI(MP2)
C              originating from the exchange part of the "extra
C              terms" from the diagonal orbital multipliers.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION ETAAI(*), YTMAT(*), DSRHF(*), XLAMDP(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      CALL QENTER('MP2_YTV')
C
      ISYMC = ISYDEL
      ISYMD = ISYMC
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      KCVEC = 1
      KDVEC = KCVEC + NVIR(ISYMC)
      KEND1 = KDVEC + NVIR(ISYMD)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation in MP2_YTV')
      ENDIF
C
      CALL DZERO(WORK(KCVEC),NVIR(ISYMC))
      CALL DZERO(WORK(KDVEC),NVIR(ISYMD))
C
C-------------------------------------
C     Copy vector out of lambda matrix.
C-------------------------------------
C
      KOFF1 = ILMVIR(ISYMC) + IDEL - IBAS(ISYDEL)
C
      CALL DCOPY(NVIR(ISYMC),XLAMDP(KOFF1),NBAS(ISYDEL),WORK(KCVEC),1)
C
C----------------------------------------
C     Contract with symmetrized Y-matrix.
C----------------------------------------
C
      KOFF1 = IMATAB(ISYMD,ISYMC) + 1
C
      NTOTD = MAX(NVIR(ISYMD),1)
C
      CALL DGEMV('N',NVIR(ISYMD),NVIR(ISYMC),ONE,YTMAT(KOFF1),NTOTD,
     *           WORK(KCVEC),1,ZERO,WORK(KDVEC),1)
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMA  = ISYMI
         ISYMAL = ISYMA
         ISYMBE = ISYMD
         ISALBE = MULD2H(ISYMAL,ISYMBE)
C
C----------------------------------
C        Work space allocation two.
C----------------------------------
C
         KAOINT = KEND1
         KSCRAO = KAOINT + N2BST(ISALBE)
         KMOINT = KSCRAO + NBAS(ISYMAL)*NVIR(ISYMD)
         KEND2  = KMOINT + NVIR(ISYMA)*NVIR(ISYMD)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient work space for allocation '//
     &                'in MP2_YTV')
         ENDIF
C
         CALL DZERO(WORK(KSCRAO),NBAS(ISYMAL)*NVIR(ISYMD))
         CALL DZERO(WORK(KMOINT),NVIR(ISYMA)*NVIR(ISYMD))
C
         DO 110 I = 1,NRHF(ISYMI)
C
C----------------------------------------
C           Unpack integral distribution.
C----------------------------------------
C
            KOFF2 = IDSRHF(ISALBE,ISYMI) + NNBST(ISALBE)*(I - 1) + 1
C
            CALL CCSD_SYMSQ(DSRHF(KOFF2),ISALBE,WORK(KAOINT))
C
C-------------------------------------------
C           Transform integrals to MO basis.
C-------------------------------------------
C
            KOFF3  = KAOINT + IAODIS(ISYMAL,ISYMBE)
            KOFF4  = ILMVIR(ISYMD) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTBE = MAX(NBAS(ISYMBE),1)
C
            CALL DGEMM('N','N',NBAS(ISYMAL),NVIR(ISYMD),NBAS(ISYMBE),
     *                 ONE,WORK(KOFF3),NTOTAL,XLAMDP(KOFF4),NTOTBE,
     *                 ZERO,WORK(KSCRAO),NTOTAL)
C
            KOFF5  = ILMVIR(ISYMA) + 1
C
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTA  = MAX(NVIR(ISYMA),1)
C
            CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMD),NBAS(ISYMAL),
     *                 ONE,XLAMDP(KOFF5),NTOTAL,WORK(KSCRAO),NTOTAL,
     *                 ZERO,WORK(KMOINT),NTOTA)
C
            KOFF6 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
C
            NTOTA = MAX(NVIR(ISYMA),1)
C
            CALL DGEMV('N',NVIR(ISYMA),NVIR(ISYMD),-TWO,WORK(KMOINT),
     *                 NTOTA,WORK(KDVEC),1,ONE,ETAAI(KOFF6),1)
C
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('MP2_YTV')
C
      RETURN
      END
C  /* Deck cc_kanew */
      SUBROUTINE CC_KANEW(ETAAI,ZKDIA,WORK,LWORK)
C
C     Written by Asger Halkier 10/8 - 1998
C
C     Version: 1.0
C
C     Purpose: To calculate the contributions to the right hand
C              side ETAAI from the diagonal multiplier blocks for
C              the equations for kappa-bar-0. This includes the
C              frozen core contributions.
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <maxash.h>
#include <mxorb.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION INDEXA(MXCORB)
      DIMENSION ETAAI(*), ZKDIA(*), WORK(LWORK)
#include <ccorb.h>
#include <infind.h>
#include <blocks.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <distcl.h>
#include <cbieri.h>
#include <eribuf.h>
#include <cclr.h>
#include <ccfro.h>
C
      CALL QENTER('CC_KANEW')
C
      CALL HEADER('Calculating diagonal contributions to eta-bar-0',-1)
C
      TIMETO = ZERO
      TIMETO = SECOND()
C
      ISYMOP = 1
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C  
      KAFROI = 1
      KLAMDP = KAFROI + NT1FRO(1)
      KLAMDH = KLAMDP + NLAMDT
      KT1AM  = KLAMDH + NLAMDT
      KEND1  = KT1AM  + NT1AMX
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for first allocation '//
     &             'in CC_KANEW')
      ENDIF
C
      KOFFAI = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 1
      CALL DZERO(WORK(KAFROI),NT1FRO(1))
      CALL DCOPY(NT1FRO(1),ZKDIA(KOFFAI),1,WORK(KAFROI),1)
      CALL DZERO(ZKDIA(KOFFAI),2*NT1FRO(1))
      CALL DZERO(WORK(KT1AM),NT1AMX)
C
C----------------------------------
C     Calculate the lambda matrices.
C----------------------------------
C
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND1),
     *            LWRK1)
C
      KEND1 = KLAMDH
      LWRK1 = LWORK - KEND1
C
C--------------------------------------------------------------------
C     Calculate the full MO coefficient matrix for frozen core calcs.
C--------------------------------------------------------------------
C
      IF (FROIMP) THEN
C
         KCMO  = KEND1
         KEND1 = KCMO  + NLAMDS
         LWKR1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient memory for allocation in CC_KANEW')
         ENDIF
C
         CALL CMO_ALL(WORK(KCMO),WORK(KEND1),LWRK1)
C
      ENDIF
C
C-----------------------------------
C     Start the loop over integrals.
C-----------------------------------
C
      KENDS2 = KEND1
      LWRKS2 = LWRK1
C
      IF (DIRECT) THEN
         IF (HERDIR) THEN
           CALL HERDI1(WORK(KEND1),LWRK1,IPRERI)
         ELSE
           KCCFB1 = KEND1
           KINDXB = KCCFB1 + MXPRIM*MXCONT
           KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
           LWRK1  = LWORK  - KEND1
           CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     *                 KODPP1,KODPP2,KRDPP1,KRDPP2,
     *                 KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     *                 WORK(KEND1),LWRK1,IPRERI)
           KEND1 = KFREE
           LWRK1 = LFREE
         END IF
         NTOSYM = 1
      ELSE
         NTOSYM = NSYM
      ENDIF
C
      KENDSV = KEND1
      LWRKSV = LWRK1
C
      ICDEL1 = 0
      DO 100 ISYMD1 = 1,NTOSYM
C
         IF (DIRECT) THEN
            IF (HERDIR) THEN
              NTOT = MAXSHL
            ELSE
              NTOT = MXCALL
            END IF
         ELSE
            NTOT = NBAS(ISYMD1)
         ENDIF
C
         DO 110 ILLL = 1,NTOT
C
C---------------------------------------------
C           If direct calculate the integrals.
C---------------------------------------------
C
            IF (DIRECT) THEN
C
               KEND1 = KENDSV
               LWRK1 = LWRKSV
C
               DTIME  = SECOND()
               IF (HERDIR) THEN
                 CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,
     &                       IPRERI)
               ELSE
                 CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     *                       WORK(KODCL1),WORK(KODCL2),
     *                       WORK(KODBC1),WORK(KODBC2),
     *                       WORK(KRDBC1),WORK(KRDBC2),
     *                       WORK(KODPP1),WORK(KODPP2),
     *                       WORK(KRDPP1),WORK(KRDPP2),
     *                       WORK(KCCFB1),WORK(KINDXB),
     *                       WORK(KEND1), LWRK1,IPRERI)
               END IF
               DTIME   = SECOND() - DTIME
               TIMHE2 = TIMHE2 + DTIME
C
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
                  CALL QUIT('Insufficient memory for integrals '//
     &                      'in CC_KANEW')
               END IF
C
            ELSE
               NUMDIS = 1
            ENDIF
C
C-----------------------------------------------------
C           Loop over number of distributions in disk.
C-----------------------------------------------------
C
            DO 120 IDEL2 = 1,NUMDIS
C
               IF (DIRECT) THEN
                  IDEL  = INDEXA(IDEL2)
                  ISYMD = ISAO(IDEL)
               ELSE
                  IDEL  = IBAS(ISYMD1) + ILLL
                  ISYMD = ISYMD1
               ENDIF
C
C----------------------------------------
C              Work space allocation two.
C----------------------------------------
C
               ISYDIS = MULD2H(ISYMD,ISYMOP)
C
               KXINT  = KEND1
               KEND2  = KXINT + NDISAO(ISYDIS)
               LWRK2  = LWORK - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
                  CALL QUIT('Insufficient memory for integrals '//
     &                      'in CC_KANEW')
               ENDIF
C
C--------------------------------------------
C              Read AO integral distribution.
C--------------------------------------------
C
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                     WORK(KRECNR),DIRECT)
C
C------------------------------------------
C              Work space allocation three.
C------------------------------------------
C
               KDSRHF = KEND2
               K3OINT = KDSRHF + NDSRHF(ISYMD)
               IF (FROIMP) THEN
                  KDSFRO = K3OINT + NMAIJK(ISYDIS)
                  KOFOIN = KDSFRO + NDSFRO(ISYDIS)
                  KEND3  = KOFOIN + NOFROO(ISYDIS) 
               ELSE
                  KEND3  = K3OINT + NMAIJK(ISYDIS)
               ENDIF
               LWRK3  = LWORK  - KEND3
C
               IF (LWRK3 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND3,'Available : ',LWORK
                  CALL QUIT('Insufficient memory for integrals '//
     &                      'in CC_KANEW')
               ENDIF
C
C----------------------------------------------------------------------
C              Transform one index in the integral batch to correlated.
C----------------------------------------------------------------------
C
               CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KLAMDP),ISYMOP,
     *                     WORK(KEND3),LWRK3,ISYDIS)
C
C------------------------------------------------------------------
C              Transform one index in the integral batch to frozen.
C------------------------------------------------------------------
C
               IF (FROIMP) THEN
C
                  CALL CC_GTOFRO(WORK(KXINT),WORK(KDSFRO),WORK(KCMO),
     *                           WORK(KEND3),LWRK3,ISYDIS)
C
C--------------------------------------------------------------
C                 Calculate integral batch (cor fro | cor del).
C--------------------------------------------------------------
C
                  CALL CC_OFROIN(WORK(KDSRHF),WORK(KOFOIN),WORK(KCMO),
     *                           WORK(KEND3),LWRK3,ISYDIS)
C
C-------------------------------------------------------------------------
C                 Calculate indirect virtual contribution to frozen block.
C-------------------------------------------------------------------------
C
                  CALL MP2_EIDV1(WORK(KAFROI),WORK(KDSFRO),
     *                           ZKDIA(NMATIJ(1)+1),WORK(KCMO),
     *                           WORK(KEND3),LWRK3,IDEL,ISYMD)
C
                  CALL MP2_EIDV2(WORK(KAFROI),WORK(KDSFRO),
     *                           ZKDIA(NMATIJ(1)+1),WORK(KCMO),
     *                           WORK(KEND3),LWRK3,IDEL,ISYMD)
C
C----------------------------------------------------------------------------
C                 Calculate indirect correlated contribution to frozen block.
C----------------------------------------------------------------------------
C
                  CALL MP2_EIDC1(WORK(KAFROI),WORK(KDSFRO),
     *                           ZKDIA(1),WORK(KCMO),WORK(KEND3),
     *                           LWRK3,IDEL,ISYMD)
C
                  CALL MP2_EIDC2(WORK(KAFROI),WORK(KOFOIN),
     *                           ZKDIA(1),WORK(KCMO),WORK(KEND3),
     *                           LWRK3,IDEL,ISYMD)
C
C-----------------------------------------------------------------------------
C                 Calculate indirect frozen contribution to both parts of eta.
C-----------------------------------------------------------------------------
C
                  KOFFJK = NMATIJ(1)   + NMATAB(1) + 2*NT1AMX
     *                   + 2*NT1FRO(1) + 1
C
                  CALL MP2_EIDF1(ETAAI,WORK(KOFOIN),ZKDIA(KOFFJK),
     *                           WORK(KCMO),WORK(KEND3),LWRK3,
     *                           IDEL,ISYMD)
C
                  CALL MP2_EIDF2(WORK(KAFROI),WORK(KDSFRO),
     *                           ZKDIA(KOFFJK),WORK(KCMO),WORK(KEND3),
     *                           LWRK3,IDEL,ISYMD)
C
                  CALL MP2_EIDF3(ETAAI,WORK(KOFOIN),ZKDIA(KOFFJK),
     *                           WORK(KCMO),WORK(KEND3),LWRK3,
     *                           IDEL,ISYMD)
C
                  CALL MP2_EIDF4(ETAAI,WORK(KDSFRO),ZKDIA(KOFFJK),
     *                           WORK(KCMO),WORK(KEND3),LWRK3,
     *                           IDEL,ISYMD)
C
                  CALL MP2_EIDF5(WORK(KAFROI),WORK(KDSRHF),
     *                           ZKDIA(KOFFJK),WORK(KCMO),WORK(KEND3),
     *                           LWRK3,IDEL,ISYMD)
C
                  CALL MP2_EIDF6(WORK(KAFROI),WORK(KDSFRO),
     *                           ZKDIA(KOFFJK),WORK(KCMO),WORK(KEND3),
     *                           LWRK3,IDEL,ISYMD)
C
               ENDIF
C
C------------------------------------------------------------------
C              Calculate contributions involving integrals (vv|ov).
C------------------------------------------------------------------
C
               CALL MP2_YTV(ETAAI,ZKDIA(NMATIJ(1)+1),WORK(KDSRHF),
     *                      WORK(KLAMDP),WORK(KEND3),LWRK3,IDEL,ISYMD)
C
C-------------------------------------------------------------------
C              Calculate integral batch with three occupied indices.
C-------------------------------------------------------------------
C
               CALL CC_INT3O(WORK(K3OINT),WORK(KDSRHF),WORK(KLAMDP),
     *                      ISYMOP,WORK(KLAMDP),WORK(KEND3),LWRK3,
     *                      IDEL,ISYMD,LUDUM,'DUMMY')
C
C------------------------------------------------------------------
C              Calculate contributions involving integrals (oo|ov).
C------------------------------------------------------------------
C
               CALL MP2_NXY(ETAAI,ZKDIA(1),ZKDIA(NMATIJ(1)+1),
     *                      WORK(K3OINT),WORK(KDSRHF),WORK(KLAMDP),
     *                      WORK(KEND3),LWRK3,IDEL,ISYMD)
C
               CALL MP2_XTO(ETAAI,ZKDIA(1),WORK(K3OINT),
     *                      WORK(KLAMDP),WORK(KEND3),LWRK3,IDEL,ISYMD)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C---------------------
C     Reorder results.
C---------------------
C
      CALL CC_ETARE(ETAAI,WORK(KAFROI),WORK(KENDS2),LWRKS2)
C
C---------------------------------
C     Write out result and timing.
C---------------------------------
C
      IF (IPRINT .GT. 20) THEN
C
         CALL AROUND('Eta-bar-0-ai vector exiting CC_KANEW')
C
         DO 20 ISYM = 1,NSYM
C
            WRITE(LUPRI,*) ' '
            WRITE(LUPRI,444) 'Sub-symmetry block number:', ISYM
            WRITE(LUPRI,555) '--------------------------'
  444       FORMAT(3X,A26,2X,I1)
  555       FORMAT(3X,A25)
C
            KOFF = IALLAI(ISYM,ISYM) + 1
            CALL OUTPUT(ETAAI(KOFF),1,NVIR(ISYM),1,NRHFS(ISYM),
     *                  NVIR(ISYM),NRHFS(ISYM),1,LUPRI)
C
            IF ((NVIR(ISYM) .EQ. 0) .OR. (NRHFS(ISYM) .EQ. 0)) THEN
               WRITE(LUPRI,*) 'This sub-symmetry is empty'
            ENDIF
C
  20     CONTINUE
      ENDIF
C
      IF (IPRINT .GT. 9) THEN
         ETAKAN = DDOT(NALLAI(1),ETAAI,1,ETAAI,1)
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'Norm of Eta-bar-0:', ETAKAN
      ENDIF
C
      TIMETO = SECOND() - TIMETO
C
      IF (IPRINT .GT. 3) THEN
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'CCSD Eta-bar-0 calculation completed'
         WRITE(LUPRI,*) 'Total time used in CC_KANEW:', TIMETO
      ENDIF
C
      CALL QEXIT('CC_KANEW')
C
      RETURN
      END
C  /* Deck cc_2eexp */
      SUBROUTINE CC_2EEXP(WORK,LWORK,IOPREL)
C
C     Written by Asger Halkier january 1999.
C
C     Version: 1.0
C
C     Purpose: To calculate the contribution to the gradient
C              from the derivative two-electron integrals 
C              using the Coupled Cluster density matrices and
C              the new integral program!
C
C     Current models: CCS, MP2, CCD, CCSD
C
C     CC2 (without frozen core) by A. Halkier & S. Coriani 20/01-2000.
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <maxash.h>
#include <mxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <aovec.h>
#include <iratdef.h>
#include <nuclei.h>
#include <symmet.h>
#include <chrnos.h>
#include <eridst.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (FOUR = 4.0D0)
      LOGICAL SAVDIR, LEX, SAVHER, OLDDX
      DIMENSION INDEXA(MXCORB)
      DIMENSION IADR(MXCORB,MXDIST)
      DIMENSION WORK(LWORK)
      CHARACTER*8 LABEL
      CHARACTER*10 MODEL
#include <ccorb.h>
#include <infind.h>
#include <blocks.h>
#include <ccfield.h>
#include <ccfop.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <distcl.h>
#include <cbieri.h>
#include <eribuf.h>
#include <cclr.h>
#include <ccfro.h>
#include <drw2el.h>
C
      CALL QENTER('CC_2EEXP')
C
C------------------------------
C     Initialization of result.
C------------------------------
C
      IF (IPRINT .GT. 9) CALL AROUND('Entering CC_2EEXP')
      CALL FLSHFO(LUPRI)
      RE2DAR = ZERO
      IF (IOPREL .EQ. 1) RELCO1 = WORK(1)
C
C-----------------------------------------
C     Initialization of timing parameters.
C-----------------------------------------
C
      TIMTOT = ZERO
      TIMTOT = SECOND()
      TIMDEN = ZERO
      TIMDAO = ZERO
      TIRDAO = ZERO
      TIMHE2 = ZERO
      TIMONE = ZERO
      TIMONE = SECOND()
C
C----------------------------------------------------
C     Both zeta- and t-vectors are totally symmetric.
C----------------------------------------------------
C
      ISYMTR = 1
      ISYMOP = 1
C
      IF (CC2) THEN
C
C
C-----------------------------------
C     Initial work space allocation.
C-----------------------------------
C
         N2BSTM = 0
         DO ISYM = 1, NSYM
           N2BSTM = MAX(N2BSTM,N2BST(ISYM))
         END DO

         KFCKEF = 1
         KAODEN = KFCKEF + N2BST(1)
         KCMO   = KAODEN + N2BSTM
         KT2AM  = KCMO   + NLAMDS
         KZ2AM  = KT2AM  + NT2AMX
         KLAMDP = KZ2AM  + NT2SQ(1)
         KLAMDH = KLAMDP + NLAMDT
         KT1AM  = KLAMDH + NLAMDT
         KZ1AM  = KT1AM  + NT1AMX
         KEND1  = KZ1AM  + NT1AMX
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT(
     *      'Insufficient core for initial allocation in CC_2EEXP')
         ENDIF
C
C-------------------------------------------------------------
C        Read MO-coefficients from interface file and reorder.
C-------------------------------------------------------------
C
         LUSIFC = -993
         CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND LUSIFC
         CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
         READ (LUSIFC)
         READ (LUSIFC)
         READ (LUSIFC) (WORK(KCMO+I-1), I=1,NLAMDS)
         CALL GPCLOSE(LUSIFC,'KEEP')
C
         CALL CMO_REORDER(WORK(KCMO),WORK(KEND1),LWRK1)
C
C-------------------------------------------
C        Read zero'th order zeta amplitudes.
C-------------------------------------------
C
         IOPT   = 3
         CALL CC_RDRSP('L0',0,1,IOPT,MODEL,WORK(KZ1AM),WORK(KZ2AM))
C
C-----------------------------------
C        Square up zeta2 amplitudes.
C-----------------------------------
C
         CALL DCOPY(NT2AMX,WORK(KZ2AM),1,WORK(KT2AM),1)
         CALL CC_T2SQ(WORK(KT2AM),WORK(KZ2AM),1)
C
C----------------------------------------------
C        Read zero'th order cluster amplitudes.
C----------------------------------------------
C
         IOPT = 3
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),WORK(KT2AM))
C
C-------------------------------------
C        Calculate the lambda matrices.
C-------------------------------------
C
         CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND1),
     *               LWRK1)
C
C
C-----------------------------------------------
C     Set up 2C-E of cluster amplitudes and save
C     in KT2AM, as we only need T(2c-e) below.
C-----------------------------------------------
C
         ISYOPE = 1
         IOPTTCME = 1
         CALL CCSD_TCMEPK(WORK(KT2AM),1.0D0,ISYOPE,IOPTTCME)
         KT2AMT = KT2AM                  !for safety
C
C-------------------------------
C     Work space allocation one.
C     Note that D(ai) = ZETA(ai)
C     and both D(ia) and h(ia)
C     are stored transposed!
C-------------------------------
C
         LENBAR = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 2*NT1FRO(1)
     *          + 2*NCOFRO(1)
C
         KONEAI = KZ1AM
         KONEAB = KONEAI + NT1AMX
         KONEIJ = KONEAB + NMATAB(1)
         KONEIA = KONEIJ + NMATIJ(1)
         KONINT = KONEIA + NT1AMX
         KKABAR = KONINT + N2BST(ISYMOP)
         KDHFAO = KKABAR + LENBAR
         KKABAO = KDHFAO + N2BST(1)
         KINTIJ = KKABAO + N2BST(1)
         KEND1  = KINTIJ + NMATIJ(1)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient core for allocation 1 in CC_2EEXP')
         ENDIF
C
C
C------------------------------------------------------
C     Initialize remaining one electron density arrays.
C------------------------------------------------------
C
         CALL DZERO(WORK(KONEAB),NMATAB(1))
         CALL DZERO(WORK(KONEIJ),NMATIJ(1))
         CALL DZERO(WORK(KONEIA),NT1AMX)
C
C--------------------------------------------------------
C     Construct remaining blocks of one electron density.
C--------------------------------------------------------
C
         CALL DZERO(WORK(KINTIJ),NMATIJ(1))
         CALL DIJGEN(WORK(KONEIJ),WORK(KINTIJ))
         CALL DIAGEN(WORK(KONEIA),WORK(KT2AMT),WORK(KONEAI))
C
C
C--------------------------------------------------------
C     Backtransform the one electron density to AO-basis.
C--------------------------------------------------------
C
         CALL DZERO(WORK(KAODEN),N2BST(1))
C
         ISDEN = 1
         CALL CC_DENAO(WORK(KAODEN),ISDEN,WORK(KONEAI),WORK(KONEAB),
     *                 WORK(KONEIJ),WORK(KONEIA),ISDEN,WORK(KLAMDP),1,
     *                 WORK(KLAMDH),1,WORK(KEND1),LWRK1)
C
C----------------------------------------------
C     Read orbital relaxation vector from disc.
C----------------------------------------------
C
         CALL DZERO(WORK(KKABAR),LENBAR)
C
         LUCCK = -987
         CALL GPOPEN(LUCCK,'CCKABAR0','UNKNOWN',' ','UNFORMATTED',
     *               IDUMMY,FALSE)
         REWIND(LUCCK)
         READ(LUCCK) (WORK(KKABAR+I-1), I = 1,LENBAR)
         CALL GPCLOSE(LUCCK,'KEEP')
C
C--------------------------------------------------------------
C     Calculate ao-transformed zeta-kappa-bar-0 and HF density.
C--------------------------------------------------------------
C
         KOFDIJ = KKABAR
         KOFDAB = KOFDIJ + NMATIJ(1)
         KOFDAI = KOFDAB + NMATAB(1)
         KOFDIA = KOFDAI + NT1AMX
C
         ISDEN = 1
         CALL DZERO(WORK(KKABAO),N2BST(1))
         CALL CC_DENAO(WORK(KKABAO),ISDEN,WORK(KOFDAI),WORK(KOFDAB),
     *                 WORK(KOFDIJ),WORK(KOFDIA),ISDEN,WORK(KCMO),1,
     *                 WORK(KCMO),1,WORK(KEND1),LWRK1)
C
         CALL CCS_D1AO(WORK(KDHFAO),WORK(KEND1),LWRK1)
         IF (FROIMP .OR. FROEXP) THEN
           MODEL = 'DUMMY'
           CALL CC_FCD1AO(WORK(KDHFAO),WORK(KEND1),LWRK1,MODEL)
         ENDIF
C
C------------------------------------------------------------
C        Add orbital relaxation for effective density matrix.
C------------------------------------------------------------
C
         CALL DAXPY(N2BST(1),ONE,WORK(KKABAO),1,WORK(KAODEN),1)
C
      ELSE IF (CCSD) THEN
C
C-----------------------------------
C     Initial work space allocation.
C-----------------------------------
C
         N2BSTM = 0
         DO ISYM = 1, NSYM
           N2BSTM = MAX(N2BSTM,N2BST(ISYM))
         END DO

         KFCKEF = 1
         KAODSY = KFCKEF + N2BST(1)
         KAODEN = KAODSY + N2BSTM
         KZ2AM  = KAODEN + N2BSTM
         KT2AM  = KZ2AM  + NT2SQ(1)
         KT2AMT = KT2AM  + NT2AMX
         KLAMDP = KT2AMT + NT2AMX
         KLAMDH = KLAMDP + NLAMDT
         KT1AM  = KLAMDH + NLAMDT
         KZ1AM  = KT1AM  + NT1AMX
         KEND1  = KZ1AM  + NT1AMX
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT(
     *      'Insufficient core for first allocation in CC_2EEXP')
         ENDIF
C
C----------------------------------------
C     Read zero'th order zeta amplitudes.
C----------------------------------------
C
         IOPT   = 3
         CALL CC_RDRSP('L0',0,1,IOPT,MODEL,WORK(KZ1AM),WORK(KZ2AM))
C
C--------------------------------
C     Square up zeta2 amplitudes.
C--------------------------------
C
         CALL DCOPY(NT2AMX,WORK(KZ2AM),1,WORK(KT2AM),1)
         CALL CC_T2SQ(WORK(KT2AM),WORK(KZ2AM),1)
C
C-------------------------------------------
C     Read zero'th order cluster amplitudes.
C-------------------------------------------
C
         IOPT = 3
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),WORK(KT2AM))
C
C------------------------------------------------
C     Zero out single vectors in CCD-calculation.
C------------------------------------------------
C
         IF (CCD) THEN
            CALL DZERO(WORK(KT1AM),NT1AMX)
            CALL DZERO(WORK(KZ1AM),NT1AMX)
         ENDIF
C
C----------------------------------
C     Calculate the lambda matrices.
C----------------------------------
C
         CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND1),
     *               LWRK1)
C
C---------------------------------------
C     Set up 2C-E of cluster amplitudes.
C---------------------------------------
C
         ISYOPE = 1
C
         CALL DCOPY(NT2AMX,WORK(KT2AM),1,WORK(KT2AMT),1)
         IOPTTCME = 1
         CALL CCSD_TCMEPK(WORK(KT2AMT),1.0D0,ISYOPE,IOPTTCME)
C
C-------------------------------
C     Work space allocation one.
C     Note that D(ai) = ZETA(ai)
C     and both D(ia) and h(ia) 
C     are stored transposed!
C-------------------------------
C
         LENBAR = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 2*NT1FRO(1)
     *          + 2*NCOFRO(1)
C
         KONEAI = KZ1AM
         KONEAB = KONEAI + NT1AMX
         KONEIJ = KONEAB + NMATAB(1)
         KONEIA = KONEIJ + NMATIJ(1)
         KXMAT  = KONEIA + NT1AMX
         KYMAT  = KXMAT  + NMATIJ(1)
         KMINT  = KYMAT  + NMATAB(1)
         KONINT = KMINT  + N3ORHF(1)
         KMIRES = KONINT + N2BST(ISYMOP)
         KD1ABT = KMIRES + N3ORHF(1)
         KD1IJT = KD1ABT + NMATAB(1)
         KKABAR = KD1IJT + NMATIJ(1)
         KDHFAO = KKABAR + LENBAR
         KKABAO = KDHFAO + N2BST(1)
         KCMO   = KKABAO + N2BST(1)
         KEND1  = KCMO   + NLAMDS
         LWRK1  = LWORK  - KEND1
C
         IF (FROIMP) THEN
            KCMOF = KEND1
            KEND1 = KCMOF + NLAMDS
            LWRK1 = LWORK - KEND1
         ENDIF
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient memory for allocation 1 CC_2EEXP')
         ENDIF
C
         IF (FROIMP) THEN
C
C----------------------------------------------
C           Get the FULL MO coefficient matrix.
C----------------------------------------------
C
            CALL CMO_ALL(WORK(KCMOF),WORK(KEND1),LWRK1)
C
         ENDIF
C
C------------------------------------------------------
C     Initialize remaining one electron density arrays.
C------------------------------------------------------
C
         CALL DZERO(WORK(KONEAB),NMATAB(1))
         CALL DZERO(WORK(KONEIJ),NMATIJ(1))
         CALL DZERO(WORK(KONEIA),NT1AMX)
C
C--------------------------------------------------------
C     Calculate X-intermediate of zeta- and t-amplitudes.
C--------------------------------------------------------
C
         CALL CC_XI(WORK(KXMAT),WORK(KZ2AM),ISYMTR,WORK(KT2AM),ISYMOP,
     *                WORK(KEND1),LWRK1)
C
C--------------------------------------------------------
C     Calculate Y-intermediate of zeta- and t-amplitudes.
C--------------------------------------------------------
C
         CALL CC_YI(WORK(KYMAT),WORK(KZ2AM),ISYMTR,WORK(KT2AM),ISYMOP,
     *              WORK(KEND1),LWRK1)
C
C--------------------------------------------------------------
C     Construct three remaining blocks of one electron density.
C--------------------------------------------------------------
C
         CALL DCOPY(NMATAB(1),WORK(KYMAT),1,WORK(KONEAB),1)
         CALL CC_EITR(WORK(KONEAB),WORK(KONEIJ),WORK(KEND1),LWRK1,1)
         CALL DIJGEN(WORK(KONEIJ),WORK(KXMAT))
         CALL DIAGEN(WORK(KONEIA),WORK(KT2AMT),WORK(KONEAI))
C
C---------------------------------
C     Set up transposed densities.
C---------------------------------
C
         CALL DCOPY(NMATAB(1),WORK(KONEAB),1,WORK(KD1ABT),1)
         CALL DCOPY(NMATIJ(1),WORK(KONEIJ),1,WORK(KD1IJT),1)
         CALL CC_EITR(WORK(KD1ABT),WORK(KD1IJT),WORK(KEND1),LWRK1,1)
C
C----------------------------------------------
C     Read orbital relaxation vector from disc.
C----------------------------------------------
C
         CALL DZERO(WORK(KKABAR),LENBAR)
C
         LUCCK = -678
         CALL GPOPEN(LUCCK,'CCKABAR0','UNKNOWN',' ',
     *               'UNFORMATTED',IDUMMY,.FALSE.)
         REWIND(LUCCK)
         READ(LUCCK) (WORK(KKABAR+I-1), I = 1,LENBAR)
         CALL GPCLOSE(LUCCK,'KEEP')
C
C----------------------------------------------------------
C     Read MO-coefficients from interface file and reorder.
C----------------------------------------------------------
C
         LUSIFC = 0
         CALL GPOPEN(LUSIFC,'SIRIFC','UNKNOWN',' ',
     *               'UNFORMATTED',IDUMMY,.FALSE.)
         REWIND LUSIFC
         CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
         READ (LUSIFC)
         READ (LUSIFC)
         READ (LUSIFC) (WORK(KCMO+I-1), I=1,NLAMDS)
         CALL GPCLOSE (LUSIFC,'KEEP')
C
         CALL CMO_REORDER(WORK(KCMO),WORK(KEND1),LWRK1)
C
C--------------------------------------------------------------
C     Calculate ao-transformed zeta-kappa-bar-0 and HF density.
C--------------------------------------------------------------
C
         KOFDIJ = KKABAR
         KOFDAB = KOFDIJ + NMATIJ(1)
         KOFDAI = KOFDAB + NMATAB(1)
         KOFDIA = KOFDAI + NT1AMX
C
         ISDEN = 1
         CALL DZERO(WORK(KKABAO),N2BST(1))
         CALL CC_DENAO(WORK(KKABAO),ISDEN,WORK(KOFDAI),WORK(KOFDAB),
     *                 WORK(KOFDIJ),WORK(KOFDIA),ISDEN,WORK(KCMO),1,
     *                 WORK(KCMO),1,WORK(KEND1),LWRK1)
C
         CALL CCS_D1AO(WORK(KDHFAO),WORK(KEND1),LWRK1)
         IF (FROIMP .OR. FROEXP) THEN
           MODEL = 'DUMMY'
           CALL CC_FCD1AO(WORK(KDHFAO),WORK(KEND1),LWRK1,MODEL)
         ENDIF
C
C------------------------------------------------------------
C        Add orbital relaxation for effective density matrix.
C------------------------------------------------------------
C
         CALL DCOPY(N2BST(1),WORK(KKABAO),1,WORK(KAODEN),1)
C
C------------------------------------------------------
C        Add frozen core contributions to AO densities.
C------------------------------------------------------
C
         IF (FROIMP) THEN
C
            KOFFAI = KKABAR + NMATIJ(1) + NMATAB(1) + 2*NT1AMX
            KOFFIA = KOFFAI + NT1FRO(1)
            KOFFIJ = KOFFIA + NT1FRO(1)
            KOFFJI = KOFFIJ + NCOFRO(1)
C
            ISDEN = 1
            ICON  = 1
            CALL CC_D1FCB(WORK(KAODEN),WORK(KOFFIJ),WORK(KOFFJI),
     *                    WORK(KOFFAI),WORK(KOFFIA),WORK(KEND1),
     *                    LWRK1,ISDEN,ICON)
C
            ISDEN = 1
            ICON  = 2
            CALL CC_D1FCB(WORK(KKABAO),WORK(KOFFIJ),WORK(KOFFJI),
     *                    WORK(KOFFAI),WORK(KOFFIA),WORK(KEND1),
     *                    LWRK1,ISDEN,ICON)
C
         ENDIF
C
C------------------------------------------------------------
C     Backtransform the one electron density to AO-basis.
C     We thus have the entire effective one-electron density.
C------------------------------------------------------------
C
         ISDEN = 1
         CALL CC_DENAO(WORK(KAODEN),ISDEN,WORK(KONEAI),WORK(KONEAB),
     *                 WORK(KONEIJ),WORK(KONEIA),ISDEN,WORK(KLAMDP),1,
     *                 WORK(KLAMDH),1,WORK(KEND1),LWRK1)
C
C--------------------------------------------------------
C     Calculate M-intermediate of zeta- and t-amplitudes.
C--------------------------------------------------------
C
         CALL CC_MI(WORK(KMINT),WORK(KZ2AM),ISYMTR,WORK(KT2AM),ISYMOP,
     *              WORK(KEND1),LWRK1)
C
C--------------------------------------------------------
C     Calculate resorted M-intermediate M(imjk)->M(mkij). 
C--------------------------------------------------------
C
         CALL CC_MIRS(WORK(KMIRES),WORK(KMINT))
C
      ELSE IF (MP2) THEN
C
C---------------------------------
C     First work space allocation.
C---------------------------------
C
         N2BSTM = 0
         DO ISYM = 1, NSYM
           N2BSTM = MAX(N2BSTM,N2BST(ISYM))
         END DO
C
         LENBAR = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 2*NCOFRO(1)
     *          + 2*NT1FRO(1)
C
         KFCKEF = 1
         KAODSY = KFCKEF + N2BST(1)
         KAODEN = KAODSY + N2BSTM
         KONEAI = KAODEN + N2BSTM
         KONEAB = KONEAI + NT1AMX
         KONEIJ = KONEAB + NMATAB(1)
         KONEIA = KONEIJ + NMATIJ(1)
         KCMO   = KONEIA + NT1AMX
         KKABAR = KCMO   + NLAMDS
         KDHFAO = KKABAR + LENBAR
         KKABAO = KDHFAO + N2BST(1)
         KLAMDH = KKABAO + N2BST(1)
         KLAMDP = KLAMDH + NLAMDT
         KEND1  = KLAMDP + NLAMDT
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT(
     *      'Insufficient memory for work allocation in CC_2EEXP')
         ENDIF
C
C--------------------------
C        Initialize arrays.
C--------------------------
C
         CALL DZERO(WORK(KONEAI),NT1AMX)
         CALL DZERO(WORK(KONEAB),NMATAB(1))
         CALL DZERO(WORK(KONEIJ),NMATIJ(1))
         CALL DZERO(WORK(KONEIA),NT1AMX)
         CALL DZERO(WORK(KKABAR),LENBAR)
C
C-----------------------------------------------------------
C        Calculate correlated part of MO coefficient matrix.
C-----------------------------------------------------------
C
         CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KONEAI),
     *               WORK(KEND1),LWRK1)
         CALL DZERO(WORK(KONEAI),NT1AMX)
C
C-------------------------------------------------
C        Read orbital relaxation vector from disc.
C-------------------------------------------------
C
         LUCCK = -6347
         CALL GPOPEN(LUCCK,'CCKABAR0','UNKNOWN',' ',
     *               'UNFORMATTED',IDUMMY,.FALSE.)
         REWIND(LUCCK)
         READ(LUCCK) (WORK(KKABAR+I-1), I = 1,LENBAR)
         CALL GPCLOSE(LUCCK,'KEEP')
C
C----------------------------------------------------------------
C        Set up the relaxation (correlation) part of the density.
C----------------------------------------------------------------
C
         CALL DCOPY(NMATIJ(1),WORK(KKABAR),1,WORK(KONEIJ),1)
         CALL DCOPY(NMATAB(1),WORK(KKABAR+NMATIJ(1)),1,WORK(KONEAB),1)
         CALL DCOPY(NT1AMX,WORK(KKABAR+NMATIJ(1)+NMATAB(1)),1,
     *              WORK(KONEAI),1)
         CALL DCOPY(NT1AMX,WORK(KONEAI),1,WORK(KONEIA),1)
C
C-------------------------------------
C        Add the Hartree-Fock density.
C-------------------------------------
C
         DO 80 ISYM = 1,NSYM
            DO 85 I = 1,NRHF(ISYM)
               NII = IMATIJ(ISYM,ISYM) + NRHF(ISYM)*(I - 1) + I
               WORK(KONEIJ + NII - 1) = WORK(KONEIJ + NII - 1) + TWO
   85       CONTINUE
   80    CONTINUE
C
C--------------------------------------
C        Transform density to AO basis.
C--------------------------------------
C
         CALL DZERO(WORK(KAODEN),N2BST(1))
C
         ISDEN = 1
         CALL CC_DENAO(WORK(KAODEN),ISDEN,WORK(KONEAI),WORK(KONEAB),
     *                 WORK(KONEIJ),WORK(KONEIA),ISDEN,WORK(KLAMDP),1,
     *                 WORK(KLAMDH),1,WORK(KEND1),LWRK1)
C
C--------------------------------------------------------------
C     Calculate ao-transformed zeta-kappa-bar-0 and HF density.
C--------------------------------------------------------------
C
         KOFDIJ = KKABAR
         KOFDAB = KOFDIJ + NMATIJ(1)
         KOFDAI = KOFDAB + NMATAB(1)
         KOFDIA = KOFDAI + NT1AMX
C
         ISDEN = 1
         CALL DZERO(WORK(KKABAO),N2BST(1))
         CALL CC_DENAO(WORK(KKABAO),ISDEN,WORK(KOFDAI),WORK(KOFDAB),
     *                 WORK(KOFDIJ),WORK(KOFDIA),ISDEN,WORK(KLAMDP),1,
     *                 WORK(KLAMDH),1,WORK(KEND1),LWRK1)
C
         CALL CCS_D1AO(WORK(KDHFAO),WORK(KEND1),LWRK1)
         IF (FROIMP .OR. FROEXP) THEN
           MODEL = 'DUMMY'
           CALL CC_FCD1AO(WORK(KDHFAO),WORK(KEND1),LWRK1,MODEL)
         ENDIF
C
C-------------------------------------------
C        Get the FULL MO coefficient matrix.
C-------------------------------------------
C
         CALL CMO_ALL(WORK(KCMO),WORK(KEND1),LWRK1)
C
C------------------------------------------------------
C        Add frozen core contributions to AO densities.
C------------------------------------------------------
C
         IF (FROIMP) THEN
C
            KOFFAI = KKABAR + NMATIJ(1) + NMATAB(1) + 2*NT1AMX
            KOFFIA = KOFFAI + NT1FRO(1)
            KOFFIJ = KOFFIA + NT1FRO(1)
            KOFFJI = KOFFIJ + NCOFRO(1)
C
            ISDEN = 1
            ICON  = 1
            CALL CC_D1FCB(WORK(KAODEN),WORK(KOFFIJ),WORK(KOFFJI),
     *                    WORK(KOFFAI),WORK(KOFFIA),WORK(KEND1),
     *                    LWRK1,ISDEN,ICON)
C
            ISDEN = 1
            ICON  = 2
            CALL CC_D1FCB(WORK(KKABAO),WORK(KOFFIJ),WORK(KOFFJI),
     *                    WORK(KOFFAI),WORK(KOFFIA),WORK(KEND1),
     *                    LWRK1,ISDEN,ICON)
C
         ENDIF
C
C----------------------------------
C        Work space allocation two.
C----------------------------------
C
         KT2AM = KEND1
         KZ2AM = KT2AM + NT2AMX
         KSKOD = KZ2AM + NT2AMX
         KEND1 = KSKOD + NT1AMX
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT(
     *      'Insufficient memory for work allocation in CC_2EEXP')
         ENDIF
C
C----------------------------------------
C     Read zero'th order zeta amplitudes.
C----------------------------------------
C
         IOPT   = 3
         CALL CC_RDRSP('L0',0,1,IOPT,MODEL,WORK(KSKOD),WORK(KZ2AM))
C
C-------------------------------------------
C     Read zero'th order cluster amplitudes.
C-------------------------------------------
C
         IOPT = 3
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KSKOD),WORK(KT2AM))
C
C-----------------------------------------------------------------------
C        Set up special modified amplitudes needed in the integral loop.
C        (By doing it this way, we only need one packed vector in core
C        along with the integral distribution in the delta loop.)
C-----------------------------------------------------------------------
C
         ISYOPE = 1
         IOPTTCME = 1
         CALL CCSD_TCMEPK(WORK(KT2AM),1.0D0,ISYOPE,IOPTTCME)
         CALL DSCAL(NT2AMX,TWO,WORK(KT2AM),1)
         CALL DAXPY(NT2AMX,ONE,WORK(KZ2AM),1,WORK(KT2AM),1)
C
         KEND1 = KSKOD
         LWRK1 = LWORK - KEND1
C
      ELSE IF (CCS) THEN
C
C---------------------------------
C     First work space allocation.
C---------------------------------
C
         N2BSTM = 0
         DO ISYM = 1, NSYM
           N2BSTM = MAX(N2BSTM,N2BST(ISYM))
         END DO

         KFCKEF = 1
         KAODSY = KFCKEF + N2BST(1)
         KAODEN = KAODSY + N2BSTM
         KCMO   = KAODEN + N2BSTM
         KEND1  = KCMO   + NLAMDS
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT
     *      ('Insufficient memory for work allocation in CC_2EEXP')
         ENDIF
C
         CALL CCS_D1AO(WORK(KAODEN),WORK(KEND1),LWRK1)
         IF (FROIMP .OR. FROEXP) THEN
           MODEL = 'DUMMY'
           CALL CC_FCD1AO(WORK(KAODEN),WORK(KEND1),LWRK1,MODEL)
         ENDIF
C
C-------------------------------------------
C        Get the FULL MO coefficient matrix.
C-------------------------------------------
C
         CALL CMO_ALL(WORK(KCMO),WORK(KEND1),LWRK1)
C
      ENDIF
C
C-----------------------------------------
C     Test: calculate energy contribution.
C-----------------------------------------
C
      IF (.FALSE.) THEN
         KTEST1 = KEND1
         KENDTS = KEND1 + N2BST(1)
         LWRKTS = LWORK - KENDTS
         CALL CCRHS_ONEAO(WORK(KTEST1),WORK(KENDTS),LWRKTS)
         ECCSD1 = DDOT(N2BST(1),WORK(KTEST1),1,WORK(KAODEN),1)
      ENDIF
C
      TIMONE = SECOND() - TIMONE
      CALL FLSHFO(LUPRI)
C
C-----------------------------------
C     Start the loop over integrals.
C-----------------------------------
C
      SAVDIR = DIRECT
      SAVHER = HERDIR
      DIRECT = .TRUE.
      HERDIR = .TRUE.
C
C
      IF (IOPREL .EQ. 2) THEN
         DPTINT = .TRUE.
      ENDIF
      IF (DAR2EL) THEN
         DO2DAR = .TRUE.
         AD2DAR = .FALSE.
         S4CENT = .FALSE.
      ENDIF
C
      KEND1A = KEND1
      LWRK1A = LWRK1
C
      DTIME  = SECOND()
      IF (HERDIR) THEN
         CALL HERDI1(WORK(KEND1),LWRK1,IPRERI)
      ELSE
         KCCFB1 = KEND1
         KINDXB = KCCFB1 + MXPRIM*MXCONT
         KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
         LWRK1  = LWORK  - KEND1
C
         CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     *               KODPP1,KODPP2,KRDPP1,KRDPP2,
     *               KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     *               WORK(KEND1),LWRK1,IPRERI)
         KEND1  = KFREE
         LWRK1  = LFREE
      ENDIF
      DTIME  = SECOND() - DTIME
      TIMHE2 = TIMHE2 + DTIME
      NTOSYM = 1
C
      KENDSV = KEND1
      LWRKSV = LWRK1
C
      ICDEL1 = 0
      IF (HERDIR) THEN
         NTOT = MAXSHL
      ELSE
         NTOT = MXCALL
      ENDIF
C
      DO 100 ILLL = 1,NTOT
C
C---------------------------------------------------------------
C        Determine which delta's to be calculated in this round.
C---------------------------------------------------------------
C
         KEND1 = KENDSV
         LWRK1 = LWRKSV
C
         DTIME  = SECOND()
         IF (HERDIR) THEN
            CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,
     &                  IPRERI)
         ELSE
            CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     *                  WORK(KODCL1),WORK(KODCL2),
     *                  WORK(KODBC1),WORK(KODBC2),
     *                  WORK(KRDBC1),WORK(KRDBC2),
     *                  WORK(KODPP1),WORK(KODPP2),
     *                  WORK(KRDPP1),WORK(KRDPP2),
     *                  WORK(KCCFB1),WORK(KINDXB),
     *                  WORK(KEND1), LWRK1,IPRERI)
         ENDIF
         DTIME  = SECOND() - DTIME
         TIMHE2 = TIMHE2 + DTIME
C
         KRECNR = KEND1
         KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
         LWRK1  = LWORK  - KEND1
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient core in CC_2EEXP')
         END IF
C
C-------------------------------------------------------
C        Open file for effective two electron densities.
C-------------------------------------------------------
C
#if defined (SYS_IRIX) || defined (SYS_NEC) || defined (SYS_DEC)
         NFRL = 2
#endif
#if defined (SYS_LINUX) || defined (SYS_CONVEX) || defined (SYS_AIX) || defined (SYS_PARAGON) || defined (SYS_ALLIANT) || defined (SYS_SUN) || defined (SYS_HPUX) || defined (SYS_HAL)
         NFRL = 8
#endif
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
         NFRL = 16
#endif
C
C         !OLD VERSION
C         LDECH = N2BSTM*NFRL+1
C         OPEN(LUDE,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE='CCTWODEN',
C     *        ACCESS='DIRECT',RECL= LDECH)
C
         LDECH = N2BSTM*NFRL+1
         LUDE = 0
         CALL GPOPEN(LUDE,'CCTWODEN','UNKNOWN','DIRECT','UNFORMATTED',
     *               LDECH,OLDDX)
C
C------------------------------------------------
C        Loop over number of delta distributions.
C------------------------------------------------
C
         DO 110 IDEL2 = 1,NUMDIS
C
            IDEL  = INDEXA(IDEL2)
            ISYMD = ISAO(IDEL)
C
C-------------------------------------
C           Work space allocation two.
C-------------------------------------
C
            ISYDEN = ISYMD
C
            IF (CCSD .OR. CC2) THEN
               KD2IJG = KEND1
               KD2AIG = KD2IJG + ND2IJG(ISYDEN)
               KD2IAG = KD2AIG + ND2AIG(ISYDEN)
               KD2ABG = KD2IAG + ND2AIG(ISYDEN)
               KEND2  = KD2ABG + ND2ABG(ISYDEN)
               LWRK2  = LWORK  - KEND2
            ELSE IF (MP2) THEN
               KD2IJG = KEND1
               KD2IAG = KD2IJG + NF2IJG(ISYDEN)
               KEND2  = KD2IAG + ND2AIG(ISYDEN)
               LWRK2  = LWORK  - KEND2
            ELSE IF (CCS) THEN
               KD2IJG = KEND1
               KEND2  = KD2IJG + NF2IJG(ISYDEN)
               LWRK2  = LWORK  - KEND2
            ENDIF
C
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:',KEND2
               CALL QUIT(
     *              'Insufficient core for allocation 2 in CC_2EEXP')
            ENDIF
C
C--------------------------------------------------
C           Initialize two electron density arrays.
C--------------------------------------------------
C
            AUTIME = SECOND()
C
            CALL DZERO(WORK(KD2IJG),NF2IJG(ISYDEN))
            IF (.NOT. CCS) THEN
               CALL DZERO(WORK(KD2IAG),ND2AIG(ISYDEN))
               IF (CCSD .OR. CC2) THEN
                  CALL DZERO(WORK(KD2AIG),ND2AIG(ISYDEN))
                  CALL DZERO(WORK(KD2ABG),ND2ABG(ISYDEN))
                  CALL DZERO(WORK(KD2IJG),ND2IJG(ISYDEN))
               ENDIF
            ENDIF
C
C----------------------------------------------------------------
C           Calculate the two electron density d(pq,gamma;delta).
C----------------------------------------------------------------
C
            IF (CCSD) THEN
               CALL CC_DEN2(WORK(KD2IJG),WORK(KD2AIG),WORK(KD2IAG),
     *                      WORK(KD2ABG),WORK(KZ2AM),WORK(KT2AM),
     *                      WORK(KT2AMT),WORK(KMINT),WORK(KXMAT),
     *                      WORK(KYMAT),WORK(KONEAB),WORK(KONEAI),
     *                      WORK(KONEIA),WORK(KMIRES),WORK(KLAMDH),1,
     *                      WORK(KLAMDP),1,WORK(KEND2),LWRK2,IDEL,
     *                      ISYMD)
            ELSE IF (CC2) THEN
               CALL CC_DEN2(WORK(KD2IJG),WORK(KD2AIG),WORK(KD2IAG),
     *                      WORK(KD2ABG),WORK(KZ2AM),WORK(KT2AM),
     *                      WORK(KT2AMT),WORK(KEND2),WORK(KEND2),
     *                      WORK(KEND2),WORK(KONEAB),WORK(KONEAI),
     *                      WORK(KONEIA),WORK(KEND2),WORK(KLAMDH),1,
     *                      WORK(KLAMDP),1,WORK(KEND2),LWRK2,IDEL,ISYMD)
            ELSE IF (MP2) THEN
               CALL CCS_DEN2(WORK(KD2IJG),WORK(KCMO),WORK(KEND2),
     *                       LWRK2,IDEL,ISYMD)
               CALL MP2_DEN2(WORK(KD2IAG),WORK(KT2AM),WORK(KLAMDH),
     *                       WORK(KEND2),LWRK2,IDEL,ISYMD)
            ELSE IF (CCS) THEN
               CALL CCS_DEN2(WORK(KD2IJG),WORK(KCMO),WORK(KEND2),
     *                       LWRK2,IDEL,ISYMD)
            ENDIF
            AUTIME = SECOND() - AUTIME
            TIMDEN = TIMDEN + AUTIME
C
C---------------------------------------------------
C           Start loop over second AO-index (gamma).
C---------------------------------------------------
C
            DO 120 ISYMG = 1, NSYM
               DO 130 G  = 1, NBAS(ISYMG)
C
                  IGAM   = G + IBAS(ISYMG)
                  ISYMPQ = MULD2H(ISYMG,ISYDEN)
C
C--------------------------------------------------------
C                 Set addresses for 2-electron densities.
C--------------------------------------------------------
C
                  AUTIME = SECOND()
                  IF (CCSD .OR. CC2) THEN
                     KD2GIJ = KD2IJG + ID2IJG(ISYMPQ,ISYMG)
     *                      + NMATIJ(ISYMPQ)*(G - 1) 
                     KD2GAI = KD2AIG + ID2AIG(ISYMPQ,ISYMG)
     *                      + NT1AM(ISYMPQ)*(G - 1)
                     KD2GAB = KD2ABG + ID2ABG(ISYMPQ,ISYMG)
     *                      + NMATAB(ISYMPQ)*(G - 1)
                     KD2GIA = KD2IAG + ID2AIG(ISYMPQ,ISYMG)
     *                      + NT1AM(ISYMPQ)*(G - 1)
                  ELSE IF (MP2) THEN
                     KD2GIJ = KD2IJG + IF2IJG(ISYMPQ,ISYMG)
     *                      + NFROIJ(ISYMPQ)*(G - 1)
                     KD2GIA = KD2IAG + ID2AIG(ISYMPQ,ISYMG)
     *                      + NT1AM(ISYMPQ)*(G - 1)
                  ELSE IF (CCS) THEN
                     KD2GIJ = KD2IJG + IF2IJG(ISYMPQ,ISYMG)
     *                      + NFROIJ(ISYMPQ)*(G - 1)
                     KLAMDH = KEND3
                  ENDIF
C
C----------------------------------------------------------
C                 Calculate frozen core contributions to d.
C----------------------------------------------------------
C
                  CALL DZERO(WORK(KAODEN),N2BST(ISYMPQ))
C
                  IF ((CCSD) .AND. (FROIMP)) THEN
C
                     KFD2IJ = KEND2
                     KFD2JI = KFD2IJ + NCOFRO(ISYMPQ)
                     KFD2AI = KFD2JI + NCOFRO(ISYMPQ)
                     KFD2IA = KFD2AI + NT1FRO(ISYMPQ)
                     KFD2II = KFD2IA + NT1FRO(ISYMPQ)
                     KEND3  = KFD2II + NFROFR(ISYMPQ)
                     LWRK3  = LWORK  - KEND3
C
                     IF (LWRK3 .LT. 0) THEN
                        WRITE(LUPRI,*) 'Available:', LWORK
                        WRITE(LUPRI,*) 'Needed:', KEND3
                        CALL QUIT('Insufficient work space in CC_2EEXP')
                     ENDIF
C
                     CALL DZERO(WORK(KFD2IJ),NCOFRO(ISYMPQ))
                     CALL DZERO(WORK(KFD2JI),NCOFRO(ISYMPQ))
                     CALL DZERO(WORK(KFD2AI),NT1FRO(ISYMPQ))
                     CALL DZERO(WORK(KFD2IA),NT1FRO(ISYMPQ))
                     CALL DZERO(WORK(KFD2II),NFROFR(ISYMPQ))
C
                     CALL CC_FD2BL(WORK(KFD2II),WORK(KFD2IJ),
     *                             WORK(KFD2JI),WORK(KFD2AI),
     *                             WORK(KFD2IA),WORK(KONEIJ),
     *                             WORK(KONEAB),WORK(KONEAI),
     *                             WORK(KONEIA),WORK(KCMOF),
     *                             WORK(KLAMDH),WORK(KLAMDP),
     *                             WORK(KEND3),LWRK3,IDEL,
     *                             ISYMD,G,ISYMG)
C
                     CALL CC_FD2AO(WORK(KAODEN),WORK(KFD2II),
     *                             WORK(KFD2IJ),WORK(KFD2JI),
     *                             WORK(KFD2AI),WORK(KFD2IA),
     *                             WORK(KCMOF),WORK(KLAMDH),
     *                             WORK(KLAMDP),WORK(KEND3),LWRK3,
     *                             ISYMPQ)
C
                     CALL CC_D2GAF(WORK(KD2GIJ),WORK(KD2GAB),
     *                             WORK(KD2GAI),WORK(KD2GIA),
     *                             WORK(KONEIJ),WORK(KONEAB),
     *                             WORK(KONEAI),WORK(KONEIA),
     *                             WORK(KCMOF),IDEL,ISYMD,G,ISYMG)
C
                     KEND4 = KEND3
                     LWRK4 = LWRK3
C
                  ELSE
C
                     KEND4 = KEND2
                     LWRK4 = LWRK2
                     IF (CCS) KLAMDH = KEND4
C
                  ENDIF
                  AUTIME = SECOND() - AUTIME
                  TIMDEN = TIMDEN + AUTIME
C
C---------------------------------------------------------
C                 Backtransform density fully to AO basis.
C---------------------------------------------------------
C
                  AUTIM1 = SECOND()
                  IF (CCSD .OR. CC2) THEN
                     CALL CC_DENAO(WORK(KAODEN),ISYMPQ,
     *                             WORK(KD2GAI),WORK(KD2GAB),
     *                             WORK(KD2GIJ),WORK(KD2GIA),ISYMPQ,
     *                             WORK(KLAMDP),1,WORK(KLAMDH),1,
     *                             WORK(KEND4),LWRK4)
                  ELSE
                     CALL CCMP_DAO(WORK(KAODEN),WORK(KD2GIJ),
     *                             WORK(KD2GIA),WORK(KCMO),
     *                             WORK(KLAMDH),WORK(KEND4),
     *                             LWRK4,ISYMPQ)
                  ENDIF
C
C-----------------------------------------------------
C                 Add relaxation terms to set up 
C                 effective density. We thus have the
C                 entire effective 2-electron density.
C-----------------------------------------------------
C
                  IF (.NOT. CCS) THEN
                     ICON = 2
                     CALL CC_D2EFF(WORK(KAODEN),G,ISYMG,IDEL,ISYMD,
     *                             WORK(KKABAO),WORK(KDHFAO),ICON)
                     CALL CC_D2EFF(WORK(KAODEN),G,ISYMG,IDEL,ISYMD,
     *                             WORK(KDHFAO),WORK(KKABAO),ICON)
                  ENDIF
                  AUTIM1 = SECOND() - AUTIM1
                  TIMDAO = TIMDAO + AUTIM1
C
C-----------------------------------------------------
C                 Write effective density to disc for 
C                 subsequent use in integral program,
C                 which performs the contraction of
C                 the density with the 2 e- integrals.
C-----------------------------------------------------
C
                  AUTIME = SECOND()
                  NDAD   = NBAST*(IDEL2 - 1) + IGAM
                  NDENEL = N2BST(ISYMPQ)
                  CALL DUMP2DEN(LUDE,WORK(KAODEN),NDENEL,NDAD)
                  AUTIME = SECOND() - AUTIME
                  TIRDAO = TIRDAO + AUTIME
C
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
C
C------------------------------------------------
C        Loop over number of delta distributions.
C------------------------------------------------
C
         DO 140 IDEL2 = 1,NUMDIS
C
            IDEL   = INDEXA(IDEL2)
            ISYMD  = ISAO(IDEL)
            ISYDEN = ISYMD
C
C---------------------------------
C           Work space allocation.
C---------------------------------
C
            ISYDIS = MULD2H(ISYMD,ISYMOP)
C
            KXINT  = KEND1
            KEND2  = KXINT  + NDISAO(ISYDIS)
            LWRK2  = LWORK  - KEND2
C
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:',KEND2
               CALL QUIT('Insufficient core for allocation in CC_2EEXP')
            ENDIF
C
C-----------------------------------------
C           Read AO integral distribution.
C-----------------------------------------
C
            AUTIME = SECOND()
            CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                  WORK(KRECNR),DIRECT)
            AUTIME = SECOND() - AUTIME
            TIRDAO = TIRDAO + AUTIME
C
C---------------------------------------------------
C           Start loop over second AO-index (gamma).
C---------------------------------------------------
C
            DO 150 ISYMG = 1, NSYM
               DO 160 G  = 1, NBAS(ISYMG)
C
                  IGAM   = G + IBAS(ISYMG)
                  ISYMPQ = MULD2H(ISYMG,ISYDEN)
C
C--------------------------------------------
C                 Work space allocation four.
C--------------------------------------------
C
                  KINTAO = KEND2
                  KEND3  = KINTAO + N2BST(ISYMPQ)
                  KCHE3  = KEND3  + N2BST(ISYMPQ)
                  LWRK3  = LWORK  - KCHE3
C
                  IF (LWRK3 .LT. 0) THEN
                     WRITE(LUPRI,*) 'Available:', LWORK
                     WRITE(LUPRI,*) 'Needed:', KCHE3
                     CALL QUIT('Insufficient work space in CC_2EEXP')
                  ENDIF
C
C----------------------------------------------------
C                 Square up AO-integral distribution.
C----------------------------------------------------
C
                  KOFFIN = KXINT + IDSAOG(ISYMG,ISYDIS)
     *                   + NNBST(ISYMPQ)*(G - 1)
C
                  CALL CCSD_SYMSQ(WORK(KOFFIN),ISYMPQ,WORK(KINTAO))
C
C----------------------------------------------
C                 Read density block from disc.
C----------------------------------------------
C
                  AUTIME = SECOND()
                  NDAD   = NBAST*(IDEL2 - 1) + IGAM
                  NDENEL = N2BST(ISYMPQ)
                  CALL RETR2DEN(LUDE,WORK(KEND3),NDENEL,NDAD)
                  AUTIME = SECOND() - AUTIME
                  TIRDAO = TIRDAO + AUTIME
C
C--------------------------------------------------------
C                 calculate the 2 e- density contribution
C                 to the requested property.
C--------------------------------------------------------
C
                  RE2DAR = RE2DAR + HALF*DDOT(N2BST(ISYMPQ),
     *                     WORK(KEND3),1,WORK(KINTAO),1)
C
  160          CONTINUE
  150       CONTINUE
  140    CONTINUE
C
C---------------------------------------------------------
C        Close file with effective two electron densities.
C---------------------------------------------------------
C
         CALL GPCLOSE(LUDE,'DELETE')
C
  100 CONTINUE
C
C------------------------------------------------
C     Restore logical flags for integral program.
C------------------------------------------------
C
      DIRECT = SAVDIR
      HERDIR = SAVHER
      IF (DAR2EL) DO2DAR = .FALSE.
      IF (IOPREL .EQ. 2) THEN
         DPTINT = .FALSE.
      ENDIF
C
C----------------------
C     Print out result.
C----------------------
C
      IF (IOPREL .EQ. 2) THEN
         WORK(1) = RE2DAR
      ELSE IF ((DAR2EL).AND.(IOPREL.NE.2)) THEN
C
         IF (IOPREL .NE. 1) THEN
            CALL AROUND('Relativistic two-electron Darwin correction')
         ENDIF
C
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,131) '2-elec. Darwin term:', RE2DAR
         WRITE(LUPRI,132) '------------------- '
C
         IF (IOPREL .EQ. 1) THEN
            RELCO1 = RELCO1 + RE2DAR
            WRITE(LUPRI,*) ' '
            WRITE(LUPRI,133) 'Total relativistic correction:', RELCO1
            WRITE(LUPRI,134) '----------------------------- '
         ENDIF
C
  131    FORMAT(9X,A20,1X,F17.9)
  132    FORMAT(9X,A20)
  133    FORMAT(9X,A30,1X,F17.9)
  134    FORMAT(9X,A30)
C
      ENDIF
C
      IF (.FALSE.) THEN
C
         LUSIFC = 0
         CALL GPOPEN(LUSIFC,'SIRIFC','UNKNOWN',' ','UNFORMMATED',
     *               IDUMMY,.FALSE.)
         REWIND LUSIFC
C
         CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
         READ (LUSIFC) POTNUC
         CALL GPCLOSE (LUSIFC,'KEEP')
C
         ECCSD = ECCSD1 + RE2DAR + POTNUC
C
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'Coupled Cluster energy constructed'
         WRITE(LUPRI,*) 'from density matrices:'
         WRITE(LUPRI,*) 'CCSD-energy:', ECCSD
         WRITE(LUPRI,*) 'H1 energy, ECCSD1 = ',ECCSD1
c        WRITE(LUPRI,*) 'H2 energy, ECCSD2 = ',RE2DAR
         WRITE(LUPRI,*) 'Two-electron contribution to FODPT:',RE2DAR
         WRITE(LUPRI,*) 'Nuc. Pot. energy  = ',POTNUC
C
      ENDIF
C
C-----------------------
C     Write out timings.
C-----------------------
C
  99  TIMTOT = SECOND() - TIMTOT
C
      IF (IPRINT .GT. 3) THEN
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'Two electron first-order property'//
     *              ' calculation completed'
         WRITE(LUPRI,*) 'Total time used in CC_2EEXP:', TIMTOT
      ENDIF
      IF (IPRINT .GT. 9) THEN
         WRITE(LUPRI,*) 
     *        'Time used for setting up d(pq,ga,de):', TIMDEN
         WRITE(LUPRI,*) 
     *        'Time used for full AO backtransformation:', TIMDAO
         WRITE(LUPRI,*) 
     *        'Time used for reading and writing d and I:',TIRDAO
         WRITE(LUPRI,*) 
     *        'Time used for calculating 2 e- AO-integrals:',TIMHE2
         WRITE(LUPRI,*) 
     *        'Time used for 1 e- density & intermediates:',TIMONE
      ENDIF
C
      CALL QEXIT('CC_2EEXP')
C
      RETURN
  165 CALL QUIT('Error reading CCTWODEN')
      END
C
C/* Deck dump2den */
      SUBROUTINE DUMP2DEN(LUDE,DEN,LENDEN,NDAD)
C
C     Written by Asger Halkier 25/1 - 99.
C
C     Purpose: Write block of effective two electron density matrix
C              (DEN) to disc.
C
C
#include <implicit.h>
      DIMENSION DEN(LENDEN)
C
      CALL QENTER('DUMP2DEN')
C
      WRITE(LUDE,REC=NDAD) (DEN(I), I=1,LENDEN)
C
      CALL QEXIT('DUMP2DEN')
C
      RETURN
      END
C/* Deck retr2den */
      SUBROUTINE RETR2DEN(LUDE,DEN,LENDEN,NDAD)
C
C     Written by Asger Halkier 25/1 - 99.
C
C     Purpose: Read block of effective two electron density matrix
C              (AODEN) from disc.
C
C
#include <implicit.h>
      DIMENSION DEN(LENDEN)
C
      CALL QENTER('RETR2DEN')
C
      READ(LUDE,ERR=1000,REC=NDAD) DEN
C
      CALL QEXIT('RETR2DEN')
C
      RETURN
 1000 CALL QUIT('Error reading CCTWODEN')
      END
