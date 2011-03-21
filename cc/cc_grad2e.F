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
C  /* Deck cc_grad2e */
      SUBROUTINE CC_GRAD2E(GRADH2,WORK,LWORK)
C
C     Written by Asger Halkier january 1999.
C     Modified November 2000 for new integral program
C
C     Version: 2.0
C
C     Purpose: To calculate the contribution to the gradient
C              from the derivative two-electron integrals 
C              using the Coupled Cluster density matrices and
C              the new integral program!
C
C     Current models: CCS, MP2, CCD, CCSD
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
      LOGICAL ENERGYTEST
      PARAMETER (ENERGYTEST = .TRUE.)
      DIMENSION INDEXA(MXCORB), INDEXB(MXCORB)
      DIMENSION IPNTAB(MXCORB,2)
      DIMENSION IADR(MXCORB,MXDIST)
      DIMENSION WORK(LWORK)
      DIMENSION GRADH2(*)
      CHARACTER*8 LABEL
      CHARACTER*10 MODEL
      LOGICAL OLDDX
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
      CALL QENTER('CC_GRAD2E')
C
C------------------------------
C     Initialization of result.
C------------------------------
C
      IF (IPRINT .GT. 9) CALL AROUND('Entering CC_GRAD2E')
      CALL FLSHFO(LUPRI)
C
      CALL DZERO(GRADH2,MXCOOR)
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
      IF (CCSD) THEN
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
C     Calculate the lamda matrices.
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
         LUBAR0 = -907
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
     *          'Insufficient memory for work allocation in CC_GRAD2E')
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
         CALL GPOPEN(LUBAR0,'CCKABAR0','OLD',' ','UNFORMATTED',
     *               IDUMMY,.FALSE.)
         REWIND(LUBAR0)
         READ(LUBAR0) (WORK(KKABAR+I-1), I = 1,LENBAR)
         CALL GPCLOSE(LUBAR0,'KEEP')
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
     *         'Insufficient memory for work allocation in CC_GRAD2E')
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
C
C-----------------------------------------
C     Test: calculate energy contribution.
C-----------------------------------------
C
      IF (ENERGYTEST) THEN
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
C-------------------------------------------------------
C     Open file for effective two electron densities.
C-------------------------------------------------------
C
      LDECH = N2BSTM*2 + 1
      LUDE = -1
      CALL GPOPEN(LUDE,'CCTWODEN','UNKNOWN','DIRECT',
     *            'UNFORMATTED',LDECH,OLDDX)
C
C-----------------------------------
C     Start the loop over integrals.
C-----------------------------------
C
      IF (ENERGYTEST) ECCSD2 = ZERO
      DIRECT = .TRUE.
C
      KEND1A = KEND1
      LWRK1A = LWRK1
C
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
      ICDEL1 = 0
      NTOT   = MXCALL
C
C-------------------------------------------------------
C     Loop over sets of delta distributions:
C-------------------------------------------------------
C
      WRITE(LUPRI,*) 'number of sets:',NTOT
C
      DO 100 ILLD = 1,NTOT
C
         KEND1 = KENDSV
         LWRK1 = LWRKSV
C
C------------------------------------------------
C        Get delta indeces for this sets:
C------------------------------------------------
C
         CALL ERIIDX(ILLD,INDEXA,NUMDISD,WORK(KINDXB),IPRERI)
C
C---------------------------------------------------------------
C        if ENERGYTEST compute the undifferentiated integrals:
C---------------------------------------------------------------
C
         IF (ENERGYTEST) THEN
C
           CALL ERIDI2(ILLD,INDEXA,NUMDISD,0,0,
     *                 WORK(KODCL1),WORK(KODCL2),
     *                 WORK(KODBC1),WORK(KODBC2),
     *                 WORK(KRDBC1),WORK(KRDBC2),
     *                 WORK(KODPP1),WORK(KODPP2),
     *                 WORK(KRDPP1),WORK(KRDPP2),
     *                 WORK(KCCFB1),WORK(KINDXB),
     *                 WORK(KEND1), LWRK1,IPRERI)
C
           KRECNR = KEND1
           KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
           LWRK1  = LWORK  - KEND1
           IF (LWRK1 .LT. 0) THEN
              CALL QUIT('Insufficient core in CC_GRAD2E (ERIDI2)')
           END IF
C
         END IF
C
C------------------------------------------------
C        Loop over number of delta distributions.
C------------------------------------------------
C
         DO 110 IDEL2 = 1,NUMDISD
C
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
            IF (CCSD) THEN
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
               IF (CCSD) THEN
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
                  IF (CCSD) THEN
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
                     IF (CCS) KLAMDH = KEND3
C
                     IF (LWRK3 .LT. 0) THEN
                        WRITE(LUPRI,*) 'Available:', LWORK
                        WRITE(LUPRI,*) 'Needed:', KEND3
                        CALL QUIT( 
     *                    'Insufficient work space in CC_GRAD2E')
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
                  IF (CCSD) THEN
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
C------------------------------------------------------------
C        IF (ENERGYTEST) read the undifferentiated integrals:
C------------------------------------------------------------
C
         IF (ENERGYTEST) THEN
            DO IDEL2 = 1, NUMDISD

              IDEL   = INDEXA(IDEL2)
              ISYMD  = ISAO(IDEL)
              ISYDEN = ISYMD
              ISYDIS = ISYDEN

              KXINT = KEND1
              KEND2 = KXINT + NDISAO(ISYDIS)
              LWRK2 = LWORK - KEND2
 
              IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:',KEND2
               CALL QUIT('Insufficient memory in CC_GRAD2E')
              END IF

              CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                    WORK(KRECNR),DIRECT)

              DO ISYMG = 1, NSYM
                 DO G = 1, NBAS(ISYMG)
                  IGAM   = G + IBAS(ISYMG)
                  ISYMPQ = MULD2H(ISYMG,ISYDEN)

                  KINTAO = KEND2
                  KDENAO = KINTAO + N2BST(ISYMPQ)
                  KEND3  = KDENAO + N2BST(ISYMPQ)
                  LWRK3  = LWORK  - KEND3
 
                  IF (LWRK3 .LT. 0) THEN
                     WRITE(LUPRI,*) 'Available:', LWORK
                     WRITE(LUPRI,*) 'Needed:', KEND3
                     CALL QUIT('Insufficient work space in CC_GRAD2E')
                  ENDIF

                  KOFFIN = KXINT + IDSAOG(ISYMG,ISYDIS)
     *                   + NNBST(ISYMPQ)*(G - 1)
 
                  CALL CCSD_SYMSQ(WORK(KOFFIN),ISYMPQ,WORK(KINTAO))

                  NDAD   = NBAST*(IDEL2 - 1) + IGAM
                  NDENEL = N2BST(ISYMPQ)
                  CALL RETR2DEN(LUDE,WORK(KDENAO),NDENEL,NDAD)

                  ECCSD2 = ECCSD2 + HALF*DDOT(N2BST(ISYMPQ),
     *                     WORK(KDENAO),1,WORK(KINTAO),1)

C                 WRITE(LUPRI,*),'CC_GRAD2E> IDEL,IGAM,ECCSD2:',
C    *                                       IDEL,IGAM,ECCSD2

                 END DO
              END DO

            END DO

         END IF ! ENERGYTEST
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
C           Work space allocattion
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
C-----------------------
C     Write out timings.
C-----------------------
C
  99  TIMTOT = SECOND() - TIMTOT
C
      IF (ENERGYTEST) THEN
        CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &              .FALSE.)
        REWIND LUSIFC
        CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
        READ (LUSIFC) POTNUC
        CALL GPCLOSE(LUSIFC,'KEEP')

        WRITE(LUPRI,*) 'NUCLEAR   :',POTNUC 
        WRITE(LUPRI,*) '1e- ENERGY:',ECCSD1 
        WRITE(LUPRI,*) '2e- ENERGY:',ECCSD2 
        WRITE(LUPRI,*) 'TOTAL     :',POTNUC+ECCSD1+ECCSD2 
      END IF
C
      IF (IPRINT .GT. 3) THEN
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'Two electron derivative gradient'//
     *              ' calculation completed'
         WRITE(LUPRI,*) 'Total time used in CC_GRAD2E:', TIMTOT
      ENDIF
C
      IF (IPRINT .GT. 9) THEN
       WRITE(LUPRI,*)'Time used for setting up d(pq,ga,de):', TIMDEN
       WRITE(LUPRI,*)'Time used for full AO backtransformation:',TIMDAO
       WRITE(LUPRI,*)'Time used for reading and writing d and I:',TIRDAO
       WRITE(LUPRI,*)'Time used for calculating 2 e- AO-integrals:',
     *              TIMHE2
       WRITE(LUPRI,*)'Time used for 1 e- density & intermediates:',
     *              TIMONE
      ENDIF
C
      CALL QEXIT('CC_GRAD2E')
      RETURN
  165 CALL QUIT('Error reading CCTWODEN')
      END
