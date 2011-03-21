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
C  /* Deck cc_grad0 */
      SUBROUTINE CC_GRAD0(GRADH1,GRADFS,GRADH2,WORK,LWORK)
C
C     Written by Asger Halkier january 1999, based on CC_GRAD.
C     Modified by Christof Haettig and Sonia Coriani in November 2000
C     for new integral program, re-merged  GRADEE,GRADKE,GRADFS
C
C     Version: 2.0
C
C     Purpose: To calculate the various contribution to the gradient:
C              - from derivative one-electron integrals GRADH1
C              - reorthonomalization part GRADFS
C              - from derivative two-electron integrals GRADH2
C              using the Coupled Cluster density matrices and
C              the new integral program!
C
C     Current models: CCS, MP2, CCD, CC2, CCSD
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxash.h>
#include <maxorb.h>
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
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      DIMENSION INDEXA(MXCORB), INDEXB(MXCORB)
      DIMENSION IPNTAB(MXCORB,2)
      DIMENSION IADR(MXCORB,MXDIST)
      DIMENSION WORK(LWORK)
      DIMENSION GRADH1(*),GRADFS(*),GRADH2(*)
      CHARACTER*8 LABEL
      CHARACTER*10 MODEL
      LOGICAL OLDDX, DIRSAV
#include <ccorb.h>
#include <infind.h>
#include <blocks.h>
#include <ccfield.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccinftap.h>
#include <inftap.h>
#include <distcl.h>
#include <cbieri.h>
#include <eribuf.h>
#include <cclr.h>
#include <ccfro.h>
C
      CALL QENTER('CC_GRAD0')
C
C------------------------------
C     Initialization of result.
C------------------------------
C
      IF (IPRINT .GT. 9 .OR. LOCDBG) THEN
        CALL AROUND('Entering CC_GRAD0')
        CALL FLSHFO(LUPRI)
      END IF
C
      CALL DZERO(GRADH2,MXCOOR)
      CALL DZERO(GRADH1,MXCOOR)
      CALL DZERO(GRADFS,MXCOOR)

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
      N2BSTM = 0
      DO ISYM = 1, NSYM
         N2BSTM = MAX(N2BSTM,N2BST(ISYM))
      END DO
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
CSC   Insert CC2 the dummy way, duplicating all calls
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      IF (CC2) THEN

         IF (FROIMP) THEN
            WRITE(LUPRI,*)
     &     'Frozen-core not (yet) available for CC2 gradient'
            CALL QUIT('No frozen-core for CC2 gradient yet')
         END IF
C
C-----------------------------------
C     Initial work space allocation.
C-----------------------------------
C
         KFCKEF = 1
         KAODSY = KFCKEF + N2BST(1)
         KAODEN = KAODSY + N2BSTM
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
     &       'Insufficient core for initial allocation in CC_GRAD')
         ENDIF
C
C-------------------------------------------------------------
C        Read MO-coefficients from interface file and reorder.
C-------------------------------------------------------------
C
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
            CALL QUIT('Insufficient core for allocation 1 in CC_GRAD')
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
         LUBAR0 = -903
         CALL GPOPEN(LUBAR0,'CCKABAR0','OLD',' ','UNFORMATTED',IDUMMY,
     *               .FALSE.)
         REWIND(LUBAR0)
         READ(LUBAR0) (WORK(KKABAR+I-1), I = 1,LENBAR)
         CALL GPCLOSE(LUBAR0,'KEEP')

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
      ELSE IF (CCSD.OR.CCD) THEN
C
C-----------------------------------
C     Initial work space allocation.
C-----------------------------------
C
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
     *         'Insufficient core for first allocation in CC_GRAD2E')
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
            CALL QUIT(
     *         'Insufficient memory for allocation 1 CC_GRAD2E')
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
         LUBAR0 = -904
         CALL GPOPEN(LUBAR0,'CCKABAR0','OLD',' ','UNFORMATTED',
     *               IDUMMY,.FALSE.)
         REWIND(LUBAR0)
         READ(LUBAR0) (WORK(KKABAR+I-1), I = 1,LENBAR)
         CALL GPCLOSE(LUBAR0,'KEEP')
C
C----------------------------------------------------------
C     Read MO-coefficients from interface file and reorder.
C----------------------------------------------------------
C
         LUSIFC = -1
         CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     *               IDUMMY,.FALSE.)
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
     *          'Insufficient memory for work allocation in CC_GRAD2E')
         ENDIF
C
C--------------------------
C        Initialize arrays:
C--------------------------
C
         CALL DZERO(WORK(KONEAI),NT1AMX)
         CALL DZERO(WORK(KONEAB),NMATAB(1))
         CALL DZERO(WORK(KONEIJ),NMATIJ(1))
         CALL DZERO(WORK(KONEIA),NT1AMX)
         CALL DZERO(WORK(KKABAR),LENBAR)
C
C-----------------------------------------------------------
C        Calculate correlated part of MO coefficient matrix:
C-----------------------------------------------------------
C
         CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KONEAI),
     *               WORK(KEND1),LWRK1)
         CALL DZERO(WORK(KONEAI),NT1AMX)
C
C-------------------------------------------------
C        Read orbital relaxation vector from disc:
C-------------------------------------------------
C
         LUBAR0 = -905
         CALL GPOPEN(LUBAR0,'CCKABAR0','OLD',' ','UNFORMATTED',
     *               IDUMMY,.FALSE.)
         REWIND(LUBAR0)
         READ(LUBAR0) (WORK(KKABAR+I-1), I = 1,LENBAR)
         CALL GPCLOSE(LUBAR0,'KEEP')
         IF (DEBUG .OR. LOCDBG) WRITE(LUPRI,*) 'NORM^2 of CCKABAR0:',
     *      DDOT(LENBAR,WORK(KKABAR),1,WORK(KKABAR),1)
C
C----------------------------------------------------------------
C        Set up the relaxation (correlation) part of the density:
C----------------------------------------------------------------
C
         CALL DCOPY(NMATIJ(1),WORK(KKABAR),1,WORK(KONEIJ),1)
         CALL DCOPY(NMATAB(1),WORK(KKABAR+NMATIJ(1)),1,WORK(KONEAB),1)
         CALL DCOPY(NT1AMX,WORK(KKABAR+NMATIJ(1)+NMATAB(1)),1,
     *              WORK(KONEAI),1)
         CALL DCOPY(NT1AMX,WORK(KONEAI),1,WORK(KONEIA),1)
C
C-------------------------------------
C        Add the Hartree-Fock density:
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
C        Transform density to AO basis:
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
C        Get the FULL MO coefficient matrix:
C-------------------------------------------
C
         CALL CMO_ALL(WORK(KCMO),WORK(KEND1),LWRK1)
C
C------------------------------------------------------
C        Add frozen core contributions to AO densities:
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
C        Work space allocation two:
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
     *         'Insufficient memory for work allocation in CC_GRAD2E')
         ENDIF
C
C-------------------------------------
C     Read zero-order zeta amplitudes:
C-------------------------------------
C
         IOPT   = 3
         CALL CC_RDRSP('L0',0,1,IOPT,MODEL,WORK(KSKOD),WORK(KZ2AM))
         IF (DEBUG.OR.LOCDBG) WRITE(LUPRI,*) 'NORM^2 of L0 (doubles):',
     *      DDOT(NT2AMX,WORK(KZ2AM),1,WORK(KZ2AM),1)
C
C----------------------------------------
C     Read zero-order cluster amplitudes:
C----------------------------------------
C
         IOPT = 3
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KSKOD),WORK(KT2AM))
         IF (DEBUG.OR.LOCDBG) WRITE(LUPRI,*) 'NORM^2 of T0 (doubles):',
     *      DDOT(NT2AMX,WORK(KT2AM),1,WORK(KT2AM),1)
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
         KFCKEF = 1
         KAODSY = KFCKEF + N2BST(1)
         KAODEN = KAODSY + N2BSTM
         KCMO   = KAODEN + N2BSTM
         KEND1  = KCMO   + NLAMDS
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT(
     *         'Insufficient memory for work allocation in CC_GRAD2E')
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
#ifdef asm_debug
! ALFREDO OK ? I think you by mistake left in some debug code? See also
! asm_debug below
casm
      xnorm = ddot(n2bst(1),work(kaoden),1,work(kaoden),1)
      write(lupri,*) 'Norm of density :', xnorm
      koff = kaoden
      do isym = 1,nsym
         write(lupri,*) 'Symmetry block :',isym
         call output(work(koff),1,nbas(isym),1,nbas(isym),
     &               nbas(isym),nbas(isym),1,lupri)
         koff = koff + nbas(isym)*nbas(isym)
      end do
#endif
C
C===========================================================
C Derivative one-electron integral contribution to gradient:
C===========================================================
C
C------------------------------------------------------
C     Loop over symmetry distinct atoms and contract
C     1-electron density with h^(1) integrals (IATOM =
C     0 for zeroth-order effective Fock matrix/energy).
C------------------------------------------------------
C
      DO IATOM = 0, NUCIND
C
        MAXCOMP = 3
        IF (IATOM .EQ. 0) MAXCOMP = 1
C
        DO ICOOR  = 1, MAXCOMP
C
           ICORSY = 1
C
           IF (IATOM .GT. 0) THEN
              ISCOOR = IPTCNT(3*(IATOM-1)+ICOOR,ICORSY-1,1)
           ELSE
              ISCOOR = 1
           ENDIF
C
           IF (ISCOOR .GT. 0) THEN
              K1DHAM = KEND1
              KEND2  = K1DHAM + N2BST(ICORSY)
              LWRK2  = LWORK - KEND2
C
              IF (LWRK2 .LE. 0) THEN
                 CALL QUIT('Insufficient work space in CC_GRAD.')
              END IF
C
              IF (IATOM .EQ. 0) THEN
C
C----------------------------------------------------
C                Test: calculate energy contribution.
C----------------------------------------------------

                 CALL DZERO(WORK(K1DHAM),N2BST(1))
                 CALL CCRHS_ONEAO(WORK(K1DHAM),WORK(KEND2),LWRK2)
                 ECCSD1 = DDOT(N2BST(1),WORK(K1DHAM),1,WORK(KAODEN),1)
                 IF (IPRINT .GT. 10 .OR. LOCDBG) THEN
                    WRITE(LUPRI,*)
     &                     'ECCSD1: 1-e energy contribution', ECCSD1
                 ENDIF
C
C-------------------------------------------------------------------
C                Calculate zeroth-order effective Fock matrix contr.
C-------------------------------------------------------------------
C
                 DO ISYMA = 1, NSYM
                   KOFF1 = IAODIS(ISYMA,ISYMA)
                   CALL TRSREC(NBAS(ISYMA),NBAS(ISYMA),
     &                     WORK(KAODEN+KOFF1),WORK(KAODSY+KOFF1))
                 END DO
                 CALL DAXPY(N2BST(1),ONE,WORK(KAODEN),1,WORK(KAODSY),1)
C
                 DO ISYMA = 1, NSYM
C
                   NBASA = MAX(NBAS(ISYMA),1)
C
                   KOFF1 = KAODSY + IAODIS(ISYMA,ISYMA)
                   KOFF2 = K1DHAM + IAODIS(ISYMA,ISYMA)
                   KOFF3 = KFCKEF + IAODIS(ISYMA,ISYMA)
C
                   CALL DGEMM('N','N',NBAS(ISYMA),NBAS(ISYMA),
     &                         NBAS(ISYMA),HALF,WORK(KOFF1),NBASA,
     &                         WORK(KOFF2),NBASA,ZERO,
     &                         WORK(KOFF3),NBASA)
C
                 END DO
C
C--------------------------------------------------------
C                Test trace of the effective Fock matrix:
C--------------------------------------------------------
C
                 FTRACE = ZERO
                 DO ISYM = 1, NSYM
                   KOFF1 = KFCKEF + IAODIS(ISYM,ISYM)
                   DO I = 1, NBAS(ISYM)
                     FTRACE = FTRACE + WORK(KOFF1+NBAS(ISYM)*(I-1)+I-1)
                   END DO
                 END DO
                 IF (IPRINT .GT. 10 .OR. LOCDBG) THEN
                    WRITE(LUPRI,*)
     &                     'Trace of 1el cont. to eff. Fock:',FTRACE
                 ENDIF
C
              ELSE
C
                 LABEL = '1DHAM'//CHRNOS(ISCOOR/100)
     &                          //CHRNOS(MOD(ISCOOR,100)/10)
     &                          //CHRNOS(MOD((MOD(ISCOOR,100)),10))
C
                 CALL CCPRPAO(LABEL,.TRUE.,WORK(K1DHAM),IRREP,ISYM,IERR,
     &                        WORK(KEND2),LWRK2)
C
                 IF (IERR.NE.0 .OR. IRREP.NE.ICORSY) THEN
                    CALL QUIT('CC_DERIV: error while reading '
     &                           //LABEL//' integrals.')
                 END IF
C
C---------------------------------------------------
C                Calculate contribution to gradient:
C---------------------------------------------------
C
                 GRADH1(ISCOOR) = DDOT(N2BST(1),WORK(K1DHAM),1,
     *                                 WORK(KAODEN),1)
C
              ENDIF
C
              IF (IPRINT .GT. 10 .OR. LOCDBG) THEN
                WRITE (LUPRI,*) 'IATOM :', IATOM
                WRITE (LUPRI,*) 'ICOOR :', ICOOR
                WRITE (LUPRI,*) 'ICORSY:', ICORSY
                WRITE (LUPRI,*) 'ISCOOR:', ISCOOR
                WRITE (LUPRI,*) 'GRADH1:', GRADH1(ISCOOR)
              END IF
           END IF
        END DO
      END DO
C
      TIMONE = SECOND() - TIMONE
      CALL FLSHFO(LUPRI)
C
C==============================================
C Two-electron density dependent contributions:
C==============================================
C----------------------------------------------------
C     Open file for effective two electron densities:
C----------------------------------------------------
C
      LDECH = N2BSTM*2 + 1
      LUDE = -1
      CALL GPOPEN(LUDE,'CCTWODEN','UNKNOWN','DIRECT',
     *            'UNFORMATTED',LDECH,OLDDX)
C
C------------------------------------------------
C     Start the loop over two-electron integrals:
C------------------------------------------------
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C Change this to have non-direct option for the
C undifferentiated integrals!!
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      DIRSAV = DIRECT
      DIRECT = .TRUE.
C
      KEND1A = KEND1
      LWRK1A = LWRK1

      KCCFB1 = KEND1
      KINDXB = KCCFB1 + MXPRIM*MXCONT
      KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
      LWRK1  = LWORK  - KEND1
C
      NTOSYM = 1
      DTIME  = SECOND()
      CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     *            KODPP1,KODPP2,KRDPP1,KRDPP2,
     *            KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     *            WORK(KEND1),LWRK1,IPRERI)
      KEND1  = KFREE
      LWRK1  = LFREE
      DTIME  = SECOND() - DTIME
      TIMHE2 = TIMHE2 + DTIME
C
      KENDSV = KEND1
      LWRKSV = LWRK1
C
      NTOT   = MXCALL
C
C-------------------------------------------
C     Loop over sets of delta distributions:
C-------------------------------------------
C
      IF (LOCDBG) WRITE(LUPRI,*) 'number of sets:',NTOT
C
      ECCSD2 = ZERO
C
      DO 100 ILLD = 1,NTOT
C
         KEND1 = KENDSV
         LWRK1 = LWRKSV
C
C---------------------------------------
C        Get delta indices for this set:
C---------------------------------------
C
         CALL ERIIDX(ILLD,INDEXA,NUMDISD,WORK(KINDXB),IPRERI)
C
C-----------------------------------------------
C        Compute the undifferentiated integrals:
C-----------------------------------------------

         CALL ERIDI2(ILLD,INDEXA,NUMDISD,0,0,
     *               WORK(KODCL1),WORK(KODCL2),
     *               WORK(KODBC1),WORK(KODBC2),
     *               WORK(KRDBC1),WORK(KRDBC2),
     *               WORK(KODPP1),WORK(KODPP2),
     *               WORK(KRDPP1),WORK(KRDPP2),
     *               WORK(KCCFB1),WORK(KINDXB),
     *               WORK(KEND1), LWRK1,IPRERI)

         KRECNR = KEND1
         KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
         LWRK1  = LWORK  - KEND1
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient core in CC_GRAD2E (ERIDI2)')
         END IF

C
C------------------------------------------------
C        Loop over number of delta distributions:
C------------------------------------------------
C
         DO 110 IDEL2 = 1,NUMDISD
C
            IDEL  = INDEXA(IDEL2)
            ISYMD = ISAO(IDEL)
            ISYDEN = ISYMD
C
C-------------------------------------
C           Work space allocation two:
C-------------------------------------
C
            IF (CCSD.OR.CC2) THEN
               KD2IJG = KEND1
               KD2AIG = KD2IJG + ND2IJG(ISYDEN)
               KD2IAG = KD2AIG + ND2AIG(ISYDEN)
               KD2ABG = KD2IAG + ND2AIG(ISYDEN)
               KEND2  = KD2ABG + ND2ABG(ISYDEN)
            ELSE IF (MP2) THEN
               KD2IJG = KEND1
               KD2IAG = KD2IJG + NF2IJG(ISYDEN)
               KEND2  = KD2IAG + ND2AIG(ISYDEN)
            ELSE IF (CCS) THEN
               KD2IJG = KEND1
               KEND2  = KD2IJG + NF2IJG(ISYDEN)
            ENDIF
            LWRK2  = LWORK  - KEND2
C
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:',KEND2
               CALL QUIT(
     *           'Insufficient core for allocation 2 in CC_GRAD2E')
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
               IF (CCSD.OR.CC2) THEN
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
#ifdef asm_debug
casm
            xnorm = ddot(NF2IJG(ISYDEN),WORK(KD2IJG),1,WORK(KD2IJG),1)
            write(lupri,'(a,3i5,d20.10)') 'isyden,illd,idel2,xnorm :',
     &              isyden,illd,idel2,xnorm
#endif
C
C-------------------------------------------------------------------
C        Read in undifferentiated 2e-integrals for Eff. Fock matrix:
C-------------------------------------------------------------------
C
            KXINT = KEND2
            KEND3 = KXINT + NDISAO(ISYMD)
            LWRK3 = LWORK - KEND3

            IF (LWRK3 .LT. 0) THEN
              WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:',KEND3
              CALL QUIT('Insufficient memory in CC_GRAD2E')
            END IF

            CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND3),LWRK3,
     *                    WORK(KRECNR),DIRECT)
#ifdef asm_debug
            write(lupri,'(a,2i8,f16.10)')'idel,idel2,integral norm :',
     &         idel,idel2,
     &         ddot(ndisao(isymd),work(kxint),1,work(kxint),1)
#endif
C
C---------------------------------------------------
C           Start loop over second AO-index (gamma).
C---------------------------------------------------
C
C   for a 2-index approach, but it does not work, IDEL2 and ILLG are
C   interchanged!!!
C
C            DO 120 ILLG = 1, NTOT
C
C               CALL ERIIDX(ILLG,INDEXB,NUMDISG,WORK(KINDXB),IPRERI)
C
C               DO 130 IGAM2  = 1, NUMDISG
C
C                  IGAM = INDEXB(IGAM2)
C                  ISYMG = ISAO(IGAM)
C                  ISYMPQ = MULD2H(ISYMG,ISYDEN)
C                  G = IGAM - IBAS(ISYMG)

            DO 120 ISYMG = 1, NSYM
               DO 130 G = 1, NBAS(ISYMG)
                  IGAM = G + IBAS(ISYMG)
                  ISYMPQ = MULD2H(ISYMG,ISYDEN)
C
C--------------------------------------------------------
C                 Set addresses for 2-electron densities.
C--------------------------------------------------------
C
                  AUTIME = SECOND()
                  IF (CCSD.OR.CC2) THEN
                     KD2GIJ = KD2IJG + ID2IJG(ISYMPQ,ISYMG)
     &                      + NMATIJ(ISYMPQ)*(G - 1) 
                     KD2GAI = KD2AIG + ID2AIG(ISYMPQ,ISYMG)
     &                      + NT1AM(ISYMPQ)*(G - 1)
                     KD2GAB = KD2ABG + ID2ABG(ISYMPQ,ISYMG)
     &                      + NMATAB(ISYMPQ)*(G - 1)
                     KD2GIA = KD2IAG + ID2AIG(ISYMPQ,ISYMG)
     &                      + NT1AM(ISYMPQ)*(G - 1)
                  ELSE IF (MP2) THEN
                     KD2GIJ = KD2IJG + IF2IJG(ISYMPQ,ISYMG)
     &                      + NFROIJ(ISYMPQ)*(G - 1)
                     KD2GIA = KD2IAG + ID2AIG(ISYMPQ,ISYMG)
     &                      + NT1AM(ISYMPQ)*(G - 1)
                  ELSE IF (CCS) THEN
                     KD2GIJ = KD2IJG + IF2IJG(ISYMPQ,ISYMG)
     &                      + NFROIJ(ISYMPQ)*(G - 1)
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
                     KFD2IJ = KEND3
                     KFD2JI = KFD2IJ + NCOFRO(ISYMPQ)
                     KFD2AI = KFD2JI + NCOFRO(ISYMPQ)
                     KFD2IA = KFD2AI + NT1FRO(ISYMPQ)
                     KFD2II = KFD2IA + NT1FRO(ISYMPQ)
                     KEND4  = KFD2II + NFROFR(ISYMPQ)
                     LWRK4  = LWORK  - KEND4
C
                     IF (LWRK4 .LT. 0) THEN
                        WRITE(LUPRI,*) 'Available:', LWORK
                        WRITE(LUPRI,*) 'Needed:', KEND4
                        CALL QUIT( 
     &                    'Insufficient work space in CC_GRAD2E')
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
     *                             WORK(KEND4),LWRK4,IDEL,
     *                             ISYMD,G,ISYMG)
C
                     CALL CC_FD2AO(WORK(KAODEN),WORK(KFD2II),
     *                             WORK(KFD2IJ),WORK(KFD2JI),
     *                             WORK(KFD2AI),WORK(KFD2IA),
     *                             WORK(KCMOF),WORK(KLAMDH),
     *                             WORK(KLAMDP),WORK(KEND4),LWRK4,
     *                             ISYMPQ)
C
                     CALL CC_D2GAF(WORK(KD2GIJ),WORK(KD2GAB),
     *                             WORK(KD2GAI),WORK(KD2GIA),
     *                             WORK(KONEIJ),WORK(KONEAB),
     *                             WORK(KONEAI),WORK(KONEIA),
     *                             WORK(KCMOF),IDEL,ISYMD,G,ISYMG)
C
                     KEND5 = KEND4
                     LWRK5 = LWRK4
C
                  ELSE
C
                     KEND5 = KEND3
                     LWRK5 = LWRK3
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
                  IF (CCSD.OR.CC2) THEN
                     CALL CC_DENAO(WORK(KAODEN),ISYMPQ,
     *                             WORK(KD2GAI),WORK(KD2GAB),
     *                             WORK(KD2GIJ),WORK(KD2GIA),ISYMPQ,
     *                             WORK(KLAMDP),1,WORK(KLAMDH),1,
     *                             WORK(KEND5),LWRK5)
                  ELSE
                     CALL CCMP_DAO(WORK(KAODEN),WORK(KD2GIJ),
     *                             WORK(KD2GIA),WORK(KCMO),
     *                             WORK(KLAMDH),WORK(KEND5),
     *                             LWRK5,ISYMPQ)
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
C----------------------------------------------------------
C                 Calculate 2e- contribution to the Energy: 
C                 using 2 e- density d in memory
C----------------------------------------------------------

                  KINTAO = KEND5
                  KEND6  = KINTAO + N2BST(ISYMPQ)
                  LWRK6  = LWORK  - KEND6

                  IF (LWRK6 .LT. 0) THEN
                     WRITE(LUPRI,*) 'Available:', LWORK
                     WRITE(LUPRI,*) 'Needed:', KEND6
                     CALL QUIT('Insufficient work space in CC_GRAD2E')
                  ENDIF

                  KOFFIN = KXINT + IDSAOG(ISYMG,ISYMD)
     *                   + NNBST(ISYMPQ)*(G - 1)

                  CALL CCSD_SYMSQ(WORK(KOFFIN),ISYMPQ,WORK(KINTAO))

                  ECCSD2 = ECCSD2 + HALF*DDOT(N2BST(ISYMPQ),
     *                     WORK(KAODEN),1,WORK(KINTAO),1)
C--------------------------------------------------------------
C                 calculate the 2 e- density contribution to
C                 the zeroth-order effective Fock matrix.
C--------------------------------------------------------------
                  DO ISYMP = 1, NSYM
C
                     ISYMQ = MULD2H(ISYMP,ISYMPQ)
C
                     KOFF1 = IAODIS(ISYMP,ISYMQ)
                     KOFF2 = IAODIS(ISYMQ,ISYMP)
C
                     CALL TRSREC(NBAS(ISYMP),NBAS(ISYMQ),
     &                        WORK(KAODEN+KOFF1),WORK(KAODSY+KOFF2))
C
                  END DO
C
                  CALL DAXPY(N2BST(ISYMPQ),ONE,WORK(KAODEN),1,
     &                           WORK(KAODSY),1)
C
                  DO ISYMP = 1, NSYM
C
                     ISYMQ = MULD2H(ISYMP,ISYMPQ)
C
                     NBASP = MAX(NBAS(ISYMP),1)
                     NBASQ = MAX(NBAS(ISYMQ),1)
C
                     KOFF1 = KAODSY + IAODIS(ISYMP,ISYMQ)
                     KOFF2 = KINTAO + IAODIS(ISYMQ,ISYMP)
                     KOFF3 = KFCKEF + IAODIS(ISYMP,ISYMP)
C
                     CALL DGEMM('N','N',NBAS(ISYMP),
     &                              NBAS(ISYMP),NBAS(ISYMQ),
     &                              HALF,WORK(KOFF1),NBASP,
     &                              WORK(KOFF2),NBASQ,ONE,
     &                              WORK(KOFF3),NBASP)
C
                  END DO
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

  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
              
C------------------------------------------------------------
C   Derivative-two-electron-integral contribution to gradient
C   from a two-index approach.
C------------------------------------------------------------
C
C-------------------------------------------------
C        Loop over sets of gamma distributions:
C-------------------------------------------------
C
         DO ILLG = 1, NTOT
C
C-----------------------------------------------------
C           Get sets of gammas for this set and set up
C           pointer arrays for ERI:
C-----------------------------------------------------
C
            CALL ERIIDX(ILLG,INDEXB,NUMDISG,WORK(KINDXB),IPRERI) 
C
            CALL IZERO(IPNTAB,2*MXCORB)
C
            DO IDEL2 = 1, NUMDISD
               IDEL   = INDEXA(IDEL2)
               IPNTAB(IDEL,1) = IDEL2
            END DO
C
            DO IGAM2 = 1, NUMDISG
               IGAM   = INDEXB(IGAM2)
               IPNTAB(IGAM,2) = IGAM2
            END DO
C            
C-------------------------------------------------------------
C           Work space allocation
C-------------------------------------------------------------
C           Loop over number of delta and gamma distributions:
            KDENSITY = KEND1 
            LDENSITY = NUMDISD*NUMDISG*NBAST*NBAST
C
            KEND2    = KDENSITY + LDENSITY
            LWRK2    = LWORK - KEND2
C
            KSCRATCH = KEND2
            KEND3    = KSCRATCH + NBAST*NBAST
            LWRK3    = LWORK - KEND3
C
            IF (LWRK3 .LT. 0) THEN
               WRITE(LUPRI,*) 'Available:', LWORK
               WRITE(LUPRI,*) 'Needed:', KEND3
               CALL QUIT('Insufficient work space in CC_GRAD2E')
            ENDIF
C
            CALL DZERO(WORK(KDENSITY),LDENSITY)
C
C-------------------------------------------------------------
C           Loop over number of delta and gamma distributions:
C-------------------------------------------------------------
C
          DO IDEL2 = 1,NUMDISD
C
             IDEL   = INDEXA(IDEL2)
             ISYMD  = ISAO(IDEL)
             ISYDEN = ISYMD
C
             DO IGAM2 = 1, NUMDISG 
C
               IGAM  = INDEXB(IGAM2)
               ISYMG = ISAO(IGAM)
C
               ISYMPQ = MULD2H(ISYMG,ISYDEN)
C
C----------------------------------------------
C              Read density block from disc.
C----------------------------------------------
C
               AUTIME = SECOND()
               NDAD   = NBAST*(IDEL2 - 1) + IGAM
               NDENEL = N2BST(ISYMPQ)
               CALL RETR2DEN(LUDE,WORK(KSCRATCH),NDENEL,NDAD)
               AUTIME = SECOND() - AUTIME
               TIRDAO = TIRDAO + AUTIME
C
C----------------------------------------------
C              expand density matrix:
C----------------------------------------------
C
               IADRDEN = KDENSITY + 
     &            (NUMDISG*(IDEL2-1)+IGAM2-1)*NBAST*NBAST
C
               DO ISYMA = 1, NSYM
                  ISYMB = MULD2H(ISYMPQ,ISYMA)

                  DO A = 1, NBAS(ISYMA)
                  DO B = 1, NBAS(ISYMB)
                     IALP = A + IBAS(ISYMA)
                     IBET = B + IBAS(ISYMB)

                     KOFF1A = KSCRATCH-1 + IAODIS(ISYMA,ISYMB) + 
     &                            NBAS(ISYMA)*(B-1) + A
                     KOFF1B = KSCRATCH-1 + IAODIS(ISYMB,ISYMA) + 
     &                            NBAS(ISYMB)*(A-1) + B
                     KOFF2 = IADRDEN-1 + NBAST*(IBET-1) + IALP

                     WORK(KOFF2) = HALF*( WORK(KOFF1A) + WORK(KOFF1B) )

                  END DO
                  END DO
               END DO
C
C----------------------------------------------
C         close loop over distributions:     
C----------------------------------------------
             END DO ! IGAM2 
          END DO ! IDEL2 
C
C---------------------------------------------------------------
C        Call ERI to compute derivative integrals and to 
C        contract them with the density
C---------------------------------------------------------------
C
         DTIME  = SECOND()
         CALL ERIDID(ILLD,ILLG,WORK(KDENSITY),IPNTAB,
     &               NUMDISD,NUMDISG,
     *               WORK(KODCL1),WORK(KODCL2),
     *               WORK(KODBC1),WORK(KODBC2),
     *               WORK(KRDBC1),WORK(KRDBC2),
     *               WORK(KODPP1),WORK(KODPP2),
     *               WORK(KRDPP1),WORK(KRDPP2),
     *               WORK(KCCFB1),WORK(KINDXB),
     *               WORK(KEND2), LWRK2,IPRERI)
         DTIME  = SECOND() - DTIME
         TIMHE2 = TIMHE2 + DTIME
c 
C--------------------------------------------
C     Close loops over sets of distributions:
C--------------------------------------------
C
         END DO ! ILLG 
C
  100 CONTINUE
C
      CALL GPCLOSE(LUDE,'DELETE')
c 
C---------------------------------------------
C     Test trace of the effective Fock matrix:
C---------------------------------------------
C
      FTRACE = ZERO
      DO ISYM = 1, NSYM
         KOFF1 = KFCKEF + IAODIS(ISYM,ISYM)
         DO I = 1, NBAS(ISYM)
            FTRACE = FTRACE + WORK(KOFF1+NBAS(ISYM)*(I-1)+I-1) 
         END DO
      END DO

      IF (LOCDBG) THEN
        FNORM = DDOT(N2BST(1),WORK(KFCKEF),1,WORK(KFCKEF),1)
        WRITE(LUPRI,*) 
     &    'Norm^2 of the eff. Fock matrix before transformation:',FNORM
        ! CALL OUTPUT(WORK(KFCKEF),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C------------------------------------------------------
C     Transform effective Fock matrix to contravariant.
C
C     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C     NOTE: Change this, so S^(0) is read in from disc.
C     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C------------------------------------------------------
C
      KEFFCO = KEND1A
      KCMODE = KEFFCO + N2BST(1)
      KEND2  = KCMODE + N2BST(1)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
         CALL QUIT( 
     &     'Insufficient memory for initial allocation in CC_GRAD')
      ENDIF
C
      CALL DZERO(WORK(KCMODE),N2BST(1))
      CALL DCOPY(N2BST(1),WORK(KFCKEF),1,WORK(KEFFCO),1)
C
      IF ((CCSD) .AND. (FROIMP)) THEN
         CALL DCOPY(NLAMDS,WORK(KCMOF),1,WORK(KCMO),1)
      ENDIF
C
      DO ISYM = 1,NSYM
C
         NTOT  = MAX(NBAS(ISYM),1)
C
         KOFF1 = KCMO   + ILVISI(ISYM)
         KOFF2 = KCMODE + IAODIS(ISYM,ISYM)
C
         CALL DGEMM('N','T',NBAS(ISYM),NBAS(ISYM),NVIRS(ISYM),ONE,
     *              WORK(KOFF1),NTOT,WORK(KOFF1),NTOT,ONE,
     *              WORK(KOFF2),NTOT)
C
         KOFF3 = KCMO   + ILRHSI(ISYM)
         KOFF4 = KCMODE + IAODIS(ISYM,ISYM)
C
         CALL DGEMM('N','T',NBAS(ISYM),NBAS(ISYM),NRHFS(ISYM),ONE,
     *              WORK(KOFF3),NTOT,WORK(KOFF3),NTOT,ONE,
     *              WORK(KOFF4),NTOT)
C
      END DO
C
      IF (LOCDBG) THEN
        XNORM = DDOT(N2BST(1),WORK(KCMODE),1,WORK(KCMODE),1)
        WRITE(LUPRI,*) 
     &    'Norm^2 of the Overlap matrix:',FNORM
        ! CALL OUTPUT(WORK(KCMODE),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      DO ISYM = 1,NSYM
C
         NTOT  = MAX(NBAS(ISYM),1)
C
         KOFF5 = KEFFCO + IAODIS(ISYM,ISYM)
         KOFF6 = KCMODE + IAODIS(ISYM,ISYM)
         KOFF7 = KFCKEF + IAODIS(ISYM,ISYM)
C
         CALL DGEMM('N','N',NBAS(ISYM),NBAS(ISYM),NBAS(ISYM),ONE,
     *              WORK(KOFF5),NTOT,WORK(KOFF6),NTOT,ZERO,
     *              WORK(KOFF7),NTOT)
C
      END DO
C
      IF (LOCDBG) THEN
        FNORM = DDOT(N2BST(1),WORK(KFCKEF),1,WORK(KFCKEF),1)
        WRITE(LUPRI,*) 
     &    'Norm^2 of the zeroth-order eff. Fock matrix:',FNORM
        ! CALL OUTPUT(WORK(KFCKEF),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C--------------------------------------------------------------
C     calculate reorthonormalization contributions to gradient.
C--------------------------------------------------------------
C
      DO IATOM = 1, NUCIND
        DO ICOOR  = 1, 3
C
           ICORSY = 1
           ISCOOR = IPTCNT(3*(IATOM-1)+ICOOR,ICORSY-1,1)
C
           IF (ISCOOR .GT. 0) THEN
C
              K1DOVL = KEND1A
              KEND2  = K1DOVL + N2BST(ICORSY)
              LWRK2  = LWORK - KEND2
C
              IF (LWRK2 .LE. 0) THEN
                 CALL QUIT('Insufficient work space in CC_GRAD.')
              END IF
C
              LABEL = '1DOVL'//CHRNOS(ISCOOR/100) 
     &                       //CHRNOS(MOD(ISCOOR,100)/10)
     &                       //CHRNOS(MOD((MOD(ISCOOR,100)),10))
C
              CALL CCPRPAO(LABEL,.TRUE.,WORK(K1DOVL),IRREP,ISYM,IERR,
     &                     WORK(KEND2),LWRK2)
C
              IF (IERR.NE.0 .OR. IRREP.NE.ICORSY) THEN
                 CALL QUIT('CC_DERIV: error while reading '
     &                        //LABEL//' integrals.')
              ENDIF
C
              GRADFS(ISCOOR) = DDOT(N2BST(1),WORK(K1DOVL),1,
     *                              WORK(KFCKEF),1)
C
              IF (IPRINT .GT. 10 .OR. LOCDBG) THEN
                WRITE (LUPRI,*) 'IATOM :', IATOM
                WRITE (LUPRI,*) 'ICOOR :', ICOOR
                WRITE (LUPRI,*) 'ICORSY:', ICORSY
                WRITE (LUPRI,*) 'ISCOOR:', ISCOOR
                WRITE (LUPRI,*) 'GRADFS:', GRADFS(ISCOOR)
                WRITE(LUPRI,*)
     &            'Norm^2 of the derivative overlap matrix:',
     &               DDOT(N2BST(1),WORK(K1DOVL),1,WORK(K1DOVL),1)
                WRITE(LUPRI,*) 
     &            'Norm^2 of the zeroth-order eff. Fock matrix:',
     &               DDOT(N2BST(1),WORK(KFCKEF),1,WORK(KFCKEF),1)
c               CALL HEADER('Derivative overlap matrix:',-1)
c               CALL OUTPUT(WORK(K1DOVL),1,NBAST,1,NBAST,NBAST,NBAST,
c    &             1,LUPRI)
              END IF
           END IF
        END DO
      END DO
C
C---------------------------------------------------
C     Read in potential energy and write out energy.
C---------------------------------------------------
C
C      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
C     &            .FALSE.)
C      REWIND LUSIFC
C
C      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
C      READ (LUSIFC) POTNUC
C      CALL GPCLOSE(LUSIFC,'KEEP')
C
C      ECCSD = ECCSD1 + ECCSD2 + POTNUC
C
C      WRITE(LUPRI,*) ' '
C      IF (IPRINT .GT. 10 .OR. LOCDBG) THEN
C         WRITE(LUPRI,*) 'Coupled Cluster energy constructed'
C         WRITE(LUPRI,*) 'from density matrices:'
C         IF (CCS)  WRITE(LUPRI,*) 'CCS-energy:',  ECCSD
C         IF (MP2)  WRITE(LUPRI,*) 'MP2-energy:',  ECCSD
C         IF (CCD)  WRITE(LUPRI,*) 'CCD-energy:',  ECCSD
C         IF (CCSD) WRITE(LUPRI,*) 'CCSD-energy:', ECCSD
C         WRITE(LUPRI,*) 'H1 energy, ECCSD1 = ',ECCSD1
C         WRITE(LUPRI,*) 'H2 energy, ECCSD2 = ',ECCSD2
C         WRITE(LUPRI,*) 'Nuc. Pot. energy  = ',POTNUC
C         WRITE(LUPRI,*) 'FTRACE            = ',FTRACE
C
C         FNORM = DDOT(N2BST(1),WORK(KFCKEF),1,WORK(KFCKEF),1)
C         WRITE(LUPRI,*) 
C     &    'Norm^2 of the zeroth-order eff. Fock matrix:',FNORM
C      ENDIF
C
C      TIMREO = SECOND() - TIMREO
C
C
C===================================================================


      IF (IPRINT.GT.3 .OR. LOCDBG)  THEN
C
C      --------------------------------
C      Write out timings and test info:
C      --------------------------------
C
        CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &              .FALSE.)
        REWIND LUSIFC
        CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
        READ (LUSIFC) POTNUC
        CALL GPCLOSE(LUSIFC,'KEEP')

        WRITE(LUPRI,*) ' '
        WRITE(LUPRI,*) 'NUCLEAR   :',POTNUC 
        WRITE(LUPRI,*) '1e- ENERGY:',ECCSD1 
        WRITE(LUPRI,*) '2e- ENERGY:',ECCSD2 
        WRITE(LUPRI,*) 'TOTAL     :',POTNUC+ECCSD1+ECCSD2 
 
        WRITE(LUPRI,*) ' '
        WRITE(LUPRI,*) 'One electron and reorthonormalization'//
     *              ' gradient calculation completed'
        WRITE(LUPRI,*) 'Two electron derivative gradient'//
     *              ' calculation completed'

        TIMTOT = SECOND() - TIMTOT
        WRITE(LUPRI,'(A,f10.2)') 'Total time used in CC_GRAD2E:',TIMTOT

        IF (IPRINT .GT. 9) THEN
         WRITE(LUPRI,'(/6(/A,f10.2))') 
     &      'Time used for setting up d(pq,ga,de)       :',TIMDEN,
     &      'Time used for full AO backtransformation   :',TIMDAO,
     &      'Time used for reading 2 e- AO-integrals    :',TIRDAO,
     &      'Time used for calculating 2 e- AO-integrals:',TIMHE2,
     &      'Time used for 1 e- part & intermediates    :',TIMONE
C    &     ,'Time used for reorthonormalization part    :',TIMREO
        ENDIF
        CALL FLSHFO(LUPRI)
      END IF
C
C     restore DIRECT flag
C
      DIRECT = DIRSAV
C
      CALL QEXIT('CC_GRAD0')
      RETURN
  165 CALL QUIT('Error reading CCTWODEN')
      END
