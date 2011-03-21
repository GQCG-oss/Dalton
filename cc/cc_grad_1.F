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
C  /* Deck cc_grad_1 */
      SUBROUTINE CC_GRAD_1(GRADH1,GRADFS,WORK,LWORK,IOPT)
C
C     Written by Asger Halkier & Christof Haettig May/June 1998.
C     Revised January 1999 to interface with new integral program.
C     Re-revised November 2000
C
C     Version: 2.0
C
C     Purpose: To calculate various contributions to the gradient
C              using the Coupled Cluster density matrices!
C
C     January 1999: Compute only the one electron part and the
C                   reorthonormalization parts and do the heavy
C                   derivative two electron part separately later.
C
C     Current models: CCS, MP2, CC2, CCD, CCSD, CCSD(T)
C
C     IOPT = 1, gradient (never used basically)
C     IOPT = 2, used in DPT
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
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      LOGICAL LOCDBG
      PARAMETER ( LOCDBG = .FALSE. )
      PARAMETER (FOUR = 4.0D0)
      DIMENSION INDEXA(MXCORB)
      DIMENSION WORK(LWORK)
      DIMENSION GRADH1(*), GRADFS(*)
      CHARACTER*8 LABEL1
      CHARACTER*10 MODEL
!
      INTEGER LUPTIA
      CHARACTER*5 FNDPTIA
      logical LCCSD_SAV
!
C
#include <ccorb.h>
#include <ccisao.h>
#include <r12int.h>
#include <blocks.h>
#include <ccfield.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccinftap.h>
#include <ccsdio.h>
#include <inftap.h>
#include <distcl.h>
#include <cbieri.h>
#include <eribuf.h>
#include <cclr.h>
#include <ccfro.h>
C
      CALL QENTER('CC_GRAD_1')
C
C------------------------------
C     Trick for CCSD(T)
C------------------------------
C
       IF (CCPT) THEN
          LCCSD_SAV = CCSD
          CCSD = .TRUE.
       END IF
C
C-----------------------------------------
C     Initialization of result
C-----------------------------------------
C
      IF (IPRINT .GT. 9) CALL AROUND('Entering CC_GRAD')
      CALL FLSHFO(LUPRI)
C
      IF (IOPT .EQ. 1) CALL DZERO(GRADH1,MXCOOR)
      IF (IOPT .EQ. 1) CALL DZERO(GRADFS,MXCOOR)
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
      TIMREO = ZERO
      TIMONE = SECOND()
C
C----------------------------------
C     Initialization of file units.
C----------------------------------
C
      LUBAR0 = -324
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

      IF (CC2) THEN
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
         IOPTRW = 3
         CALL CC_RDRSP('L0',0,1,IOPTRW,MODEL,WORK(KZ1AM),WORK(KZ2AM))
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
         IOPTRW = 3
         CALL CC_RDRSP('R0',0,1,IOPTRW,MODEL,WORK(KT1AM),WORK(KT2AM))
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
      ELSE IF (CCD.OR.CCSD.OR.CCPT) THEN
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
     *        'Insufficient memory for first allocation in CC_GRAD')
         ENDIF
C
C----------------------------------------
C     Read zero'th order zeta amplitudes.
C----------------------------------------
C
         IOPTRW = 3
         CALL CC_RDRSP('L0',0,1,IOPTRW,MODEL,WORK(KZ1AM),WORK(KZ2AM))
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
         IOPTRW = 3
         CALL CC_RDRSP('R0',0,1,IOPTRW,MODEL,WORK(KT1AM),WORK(KT2AM))
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
            CALL QUIT('Insufficient memory for allocation 1 CC_GRAD')
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
     &              WORK(KEND1),LWRK1)
C
C--------------------------------------------------------
C     Calculate Y-intermediate of zeta- and t-amplitudes.
C--------------------------------------------------------
C
         CALL CC_YI(WORK(KYMAT),WORK(KZ2AM),ISYMTR,WORK(KT2AM),ISYMOP,
     &              WORK(KEND1),LWRK1)
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
         CALL GPOPEN(LUBAR0,'CCKABAR0','OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LUBAR0)
         READ(LUBAR0) (WORK(KKABAR+I-1), I = 1,LENBAR)
         CALL GPCLOSE(LUBAR0,'KEEP')
C
C----------------------------------------------------------
C     Read MO-coefficients from interface file and reorder.
C----------------------------------------------------------
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
C--------------------------------------------------------------
C     For CCSD(T) include contribution from the D_ia(T) density
C     We thus have the entire effective one-electron density
C     for CCSD(T) in AO basis
C--------------------------------------------------------------
C
 
         IF (CCPT) THEN
           KONEIA_PT = KEND1
           KEND1     = KONEIA_PT + NT1AMX
           LWRK1     = LWORK - KEND1
           IF (LWRK1 .LT. 0) THEN
              WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
              CALL QUIT(
     &       'Insufficient core for allocation in CC_2EEXP (CCSD(T))')
           ENDIF
C
           LUPTIA = -1
           FNDPTIA  = 'DPTIA'
           CALL WOPEN2(LUPTIA,FNDPTIA,64,0)
           IOFF = 1
           CALL GETWA2(LUPTIA,FNDPTIA,WORK(KONEIA_PT),IOFF,NT1AMX)
!           IF (LOCDBG) THEN
!             WRITE(LUPRI,*)'Norm of (T) D_ia:',
!     &       DDOT(NT1AMX,WORK(KONEIA_PT),1,WORK(KONEIA_PT),1)
!           END IF
           CALL WCLOSE2(LUPTIA,FNDPTIA,'KEEP')
 
           KONEAI_PT = KEND1
           KONEIJ_PT = KONEAI_PT + NT1AMX
           KONEAB_PT = KONEIJ_PT + NMATIJ(1)
           KEND11    = KONEAB_PT + NMATAB(1)
           KDTEST    = KEND11
           LWRK11    = LWORK - KEND11
           IF (LWRK11 .LT. 0) THEN
              WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
              CALL QUIT(
     &       'Insufficient core for allocation in CCPT_GRAD0')
           ENDIF
 
           CALL DZERO(WORK(KONEAI_PT),NT1AMX)
           CALL DZERO(WORK(KONEAB_PT),NMATAB(1))
           CALL DZERO(WORK(KONEIJ_PT),NMATIJ(1))
 
           ISDEN = 1
           CALL CC_DENAO(WORK(KAODEN),ISDEN,WORK(KONEAI_PT),
     *                 WORK(KONEAB_PT),
     *                 WORK(KONEIJ_PT),WORK(KONEIA_PT),ISDEN,
     *                 WORK(KCMO),1,
     *                 WORK(KCMO),1,WORK(KEND11),LWRK11)
 
         test=ddot(N2BST(1),WORK(KAODEN),1,WORK(KAODEN),1)
         write(lupri,*) 'Effective 1-e density in AO basis All: ',test
 
         END IF
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
     *         'Insufficient memory for work allocation in CC_GRAD')
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
     *         'Insufficient memory for work allocation in CC_GRAD')
         ENDIF
C
C----------------------------------------
C     Read zero'th order zeta amplitudes.
C----------------------------------------
C
         IOPTRW = 3
         CALL CC_RDRSP('L0',0,1,IOPTRW,MODEL,WORK(KSKOD),WORK(KZ2AM))
C
C-------------------------------------------
C     Read zero'th order cluster amplitudes.
C-------------------------------------------
C
         IOPTRW = 3
         CALL CC_RDRSP('R0',0,1,IOPTRW,MODEL,WORK(KSKOD),WORK(KT2AM))
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
     *         'Insufficient memory for work allocation in CC_GRAD')
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
C------------------------------------------------------
C     Loop over symmetry distinct atoms and contract
C     1-electron density with h^(1) integrals (IATOM =
C     0 for zeroth-order effective Fock matrix/energy).
C------------------------------------------------------
C
      IF (IOPT .EQ. 2) THEN
         NUSAIN = NUCIND
         NUCIND = 0
      ENDIF
      !loop never entered for IOPT=2
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
C
                 CALL DZERO(WORK(K1DHAM),N2BST(1))
                 CALL CCRHS_ONEAO(WORK(K1DHAM),WORK(KEND2),LWRK2)
C
C     Including finite fields (WK/UniKA/11-03-2004).
      DO IF = 1, NFIELD
        IF ( .NOT. NHFFIELD(IF) ) THEN
          CALL CC_ONEP(WORK(K1DHAM),WORK(KEND2),LWRK2,EFIELD(IF),
     *                 1,LFIELD(IF))
        END IF
      END DO
C
                 ECCSD1 = DDOT(N2BST(1),WORK(K1DHAM),1,WORK(KAODEN),1)
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
                 LABEL1 = '1DHAM'//CHRNOS(ISCOOR/100)
     &                          //CHRNOS(MOD(ISCOOR,100)/10)
     &                          //CHRNOS(MOD((MOD(ISCOOR,100)),10))
C
                 CALL CCPRPAO(LABEL1,.TRUE.,WORK(K1DHAM),IRREP,ISYM,
     &                        IERR,WORK(KEND2),LWRK2)
C
                 IF (IERR.NE.0 .OR. IRREP.NE.ICORSY) THEN
                    CALL QUIT('CC_DERIV: error while reading '
     &                           //LABEL1//' integrals.')
                 END IF
C
C------------------------------------------------
C                Calculate gradient contribution:
C------------------------------------------------
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
C-----------------------------------
C     Start the loop over integrals.
C-----------------------------------
C
      ECCSD2 = ZERO
C
      KEND1A = KEND1
      LWRK1A = LWRK1
C
      IF (DIRECT) THEN
         DTIME = SECOND()
C
         IF (HERDIR) THEN
            CALL HERDI1(WORK(KEND1),LWRK1,IPRERI)
         ELSE
            KCCFB1 = KEND1
            KINDXB = KCCFB1 + MXPRIM*MXCONT
            KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
            LWRK1  = LWORK  - KEND1
            CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     *                  KODPP1,KODPP2,KRDPP1,KRDPP2,
     *                  KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     *                  WORK(KEND1),LWRK1,IPRERI)
            KEND1  = KFREE
            LWRK1  = LFREE
         ENDIF
         DTIME  = SECOND() - DTIME
         TIMHE2 = TIMHE2 + DTIME
         NTOSYM = 1
      ELSE
         NTOSYM = NSYM
      ENDIF
C
      KENDSV = KEND1
      LWRKSV = LWRK1
C
      ICDEL1 = 0
      ICDEL2 = 0
      DO 100 ISYMD1 = 1,NTOSYM
C
         IF (DIRECT) THEN
            IF (HERDIR) THEN
               NTOT = MAXSHL
            ELSE
               NTOT = MXCALL
            ENDIF
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
     &                        IPRERI)
               ELSE
                  CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     *                        WORK(KODCL1),WORK(KODCL2),
     *                        WORK(KODBC1),WORK(KODBC2),
     *                        WORK(KRDBC1),WORK(KRDBC2),
     *                        WORK(KODPP1),WORK(KODPP2),
     *                        WORK(KRDPP1),WORK(KRDPP2),
     *                        WORK(KCCFB1),WORK(KINDXB),
     *                        WORK(KEND1), LWRK1,IPRERI)
               ENDIF
               DTIME  = SECOND() - DTIME
               TIMHE2 = TIMHE2 + DTIME
C
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient core in CC_GRAD')
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
                  IF (NOAUXB) THEN
                     IDUM = 1
                     CALL IJKAUX(IDEL,IDUM,IDUM,IDUM)
                  END IF
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
               ISYDEN = ISYMD
C
               IF (CCSD .OR. CC2 .OR. CCPT) THEN
                  KD2IJG = KEND1
                  KD2AIG = KD2IJG + ND2IJG(ISYDEN)
                  KD2IAG = KD2AIG + ND2AIG(ISYDEN)
                  KD2ABG = KD2IAG + ND2AIG(ISYDEN)
                  KEND2  = KD2ABG + ND2ABG(ISYDEN)
                  LWRK2  = LWORK  - KEND2
!!!!
                  IF (CCPT) THEN
                    KD2IJG_PT = KEND2
                    KD2AIG_PT = KD2IJG_PT + ND2IJG(ISYDEN)
                    KD2IAG_PT = KD2AIG_PT + ND2AIG(ISYDEN)
                    KD2ABG_PT = KD2IAG_PT + ND2AIG(ISYDEN)
                    KEND2     = KD2ABG_PT + ND2ABG(ISYDEN)
                    LWRK2  = LWORK  - KEND2
                  END IF
!!!!
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
     &                'Insufficient space for allocation 2 in CC_GRAD')
               ENDIF
C
C-----------------------------------------------------
C              Initialize two electron density arrays.
C-----------------------------------------------------
C
               CALL DZERO(WORK(KD2IJG),NF2IJG(ISYDEN))
               IF (.NOT. CCS) THEN
                  CALL DZERO(WORK(KD2IAG),ND2AIG(ISYDEN))

                  IF (CCSD .OR. CC2 .OR. CCPT) THEN
                     CALL DZERO(WORK(KD2AIG),ND2AIG(ISYDEN))
                     CALL DZERO(WORK(KD2ABG),ND2ABG(ISYDEN))
                     CALL DZERO(WORK(KD2IJG),ND2IJG(ISYDEN))
CSONIA SONIA
                     IF (CCPT) THEN
                       CALL DZERO(WORK(KD2IAG_PT),ND2AIG(ISYDEN))
                       CALL DZERO(WORK(KD2AIG_PT),ND2AIG(ISYDEN))
                       CALL DZERO(WORK(KD2ABG_PT),ND2ABG(ISYDEN))
                       CALL DZERO(WORK(KD2IJG_PT),ND2IJG(ISYDEN))
                     END IF
CSONIA SONIA
                  ENDIF

               ENDIF
C
C-------------------------------------------------------------------
C              Calculate the two electron density d(pq,gamma;delta).
C-------------------------------------------------------------------
C
               AUTIME = SECOND()
C
               IF (CCSD.OR.CCPT) THEN
                  CALL CC_DEN2(WORK(KD2IJG),WORK(KD2AIG),WORK(KD2IAG),
     *                         WORK(KD2ABG),WORK(KZ2AM),WORK(KT2AM),
     *                         WORK(KT2AMT),WORK(KMINT),WORK(KXMAT),
     *                         WORK(KYMAT),WORK(KONEAB),WORK(KONEAI),
     *                         WORK(KONEIA),WORK(KMIRES),WORK(KLAMDH),1,
     *                         WORK(KLAMDP),1,WORK(KEND2),LWRK2,IDEL,
     *                         ISYMD)
CSONIA SONIA
                  IF (CCPT) THEN
                     CALL CC_DEN2_PT(WORK(KD2IJG_PT),WORK(KD2AIG_PT),
     &                            WORK(KD2IAG_PT),WORK(KD2ABG_PT),
     &                            WORK(KCMO),1,WORK(KEND2),LWRK2,
     &                            IDEL,ISYMD)
                  END IF
CSONIA SONIA
               ELSE IF (CC2) THEN
                  CALL CC_DEN2(WORK(KD2IJG),WORK(KD2AIG),WORK(KD2IAG),
     *                      WORK(KD2ABG),WORK(KZ2AM),WORK(KT2AM),
     *                      WORK(KT2AMT),WORK(KEND2),WORK(KEND2),
     *                      WORK(KEND2),WORK(KONEAB),WORK(KONEAI),
     *                      WORK(KONEIA),WORK(KEND2),WORK(KLAMDH),1,
     *                      WORK(KLAMDP),1,WORK(KEND2),LWRK2,IDEL,ISYMD)
               ELSE IF (MP2) THEN
                  CALL CCS_DEN2(WORK(KD2IJG),WORK(KCMO),WORK(KEND2),
     *                          LWRK2,IDEL,ISYMD)
                  CALL MP2_DEN2(WORK(KD2IAG),WORK(KT2AM),WORK(KLAMDH),
     *                          WORK(KEND2),LWRK2,IDEL,ISYMD)
               ELSE IF (CCS) THEN
                  CALL CCS_DEN2(WORK(KD2IJG),WORK(KCMO),WORK(KEND2),
     *                          LWRK2,IDEL,ISYMD)
               ENDIF
C
               AUTIME = SECOND() - AUTIME
               TIMDEN = TIMDEN + AUTIME
C
C------------------------------------------
C              Work space allocation three.
C------------------------------------------
C
               ISYDIS = MULD2H(ISYMD,ISYMOP)
C
               KXINT  = KEND2
               KEND3  = KXINT  + NDISAO(ISYDIS)
               LWRK3  = LWORK  - KEND3
C
               IF (LWRK3 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:',KEND3
                  CALL QUIT(
     &               'Insufficient space for allocation 3 in CC_GRAD')
               ENDIF
C
C--------------------------------------------
C              Read AO integral distribution.
C--------------------------------------------
C
               AUTIME = SECOND()
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND3),LWRK3,
     *                     WORK(KRECNR),DIRECT)
               AUTIME = SECOND() - AUTIME
               TIRDAO = TIRDAO + AUTIME
C
C------------------------------------------------------
C              Start loop over second AO-index (gamma).
C------------------------------------------------------
C
               DO 130 ISYMG = 1, NSYM
                  DO 140 G = 1, NBAS(ISYMG)
C
C-----------------------------------------------------------
C                    Set addresses for 2-electron densities.
C-----------------------------------------------------------
C
                     AUTIME = SECOND()
C
                     ISYMPQ = MULD2H(ISYMG,ISYDEN)
C
                     IF (CCSD .OR. CC2 .OR. CCPT) THEN
                        KD2GIJ = KD2IJG + ID2IJG(ISYMPQ,ISYMG)
     *                         + NMATIJ(ISYMPQ)*(G - 1) 
                        KD2GAI = KD2AIG + ID2AIG(ISYMPQ,ISYMG)
     *                         + NT1AM(ISYMPQ)*(G - 1)
                        KD2GAB = KD2ABG + ID2ABG(ISYMPQ,ISYMG)
     *                         + NMATAB(ISYMPQ)*(G - 1)
                        KD2GIA = KD2IAG + ID2AIG(ISYMPQ,ISYMG)
     *                         + NT1AM(ISYMPQ)*(G - 1)
CSONIA SONIA
                        IF (CCPT) THEN
 
                           KD2GIJ_PT = KD2IJG_PT + ID2IJG(ISYMPQ,ISYMG)
     &                         + NMATIJ(ISYMPQ)*(G - 1)
                           KD2GAI_PT = KD2AIG_PT + ID2AIG(ISYMPQ,ISYMG)
     &                         + NT1AM(ISYMPQ)*(G - 1)
                           KD2GAB_PT = KD2ABG_PT + ID2ABG(ISYMPQ,ISYMG)
     &                         + NMATAB(ISYMPQ)*(G - 1)
                           KD2GIA_PT = KD2IAG_PT + ID2AIG(ISYMPQ,ISYMG)
     &                         + NT1AM(ISYMPQ)*(G - 1)
                        END IF
CSONIA SONIA
                     ELSE IF (MP2) THEN
                        KD2GIJ = KD2IJG + IF2IJG(ISYMPQ,ISYMG)
     *                         + NFROIJ(ISYMPQ)*(G - 1)
                        KD2GIA = KD2IAG + ID2AIG(ISYMPQ,ISYMG)
     *                         + NT1AM(ISYMPQ)*(G - 1)
                     ELSE IF (CCS) THEN
                        KD2GIJ = KD2IJG + IF2IJG(ISYMPQ,ISYMG)
     *                         + NFROIJ(ISYMPQ)*(G - 1)
                        KLAMDH = KEND3
                     ENDIF
C
C-------------------------------------------------------------
C                    Calculate frozen core contributions to d.
C-------------------------------------------------------------
C
                     CALL DZERO(WORK(KAODEN),N2BST(ISYMPQ))
C
                     IF ((CCSD .OR. CCPT) .AND. (FROIMP)) THEN
C
                        KFD2IJ = KEND3
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
                         CALL QUIT('Insufficient work space in CC_GRAD')
                        ENDIF
C
                        CALL DZERO(WORK(KFD2IJ),NCOFRO(ISYMPQ))
                        CALL DZERO(WORK(KFD2JI),NCOFRO(ISYMPQ))
                        CALL DZERO(WORK(KFD2AI),NT1FRO(ISYMPQ))
                        CALL DZERO(WORK(KFD2IA),NT1FRO(ISYMPQ))
                        CALL DZERO(WORK(KFD2II),NFROFR(ISYMPQ))
CSONIA SONIA SONIA
                        IF (CCPT) THEN
                          KFD2IJ_PT = KEND3
                          KFD2JI_PT = KFD2IJ_PT + NCOFRO(ISYMPQ)
                          KFD2AI_PT = KFD2JI_PT + NCOFRO(ISYMPQ)
                          KFD2IA_PT = KFD2AI_PT + NT1FRO(ISYMPQ)
                          KFD2II_PT = KFD2IA_PT + NT1FRO(ISYMPQ)
                          KEND3     = KFD2II_PT + NFROFR(ISYMPQ)
                          LWRK3     = LWORK  - KEND3
                          IF (LWRK3 .LT. 0) THEN
                             WRITE(LUPRI,*) 'Available:', LWORK
                             WRITE(LUPRI,*) 'Needed:', KEND3
                             CALL QUIT(
     &                      'Insufficient work space in CC_GRAD2E')
                          ENDIF
                          CALL DZERO(WORK(KFD2IJ_PT),NCOFRO(ISYMPQ))
                          CALL DZERO(WORK(KFD2JI_PT),NCOFRO(ISYMPQ))
                          CALL DZERO(WORK(KFD2AI_PT),NT1FRO(ISYMPQ))
                          CALL DZERO(WORK(KFD2IA_PT),NT1FRO(ISYMPQ))
                          CALL DZERO(WORK(KFD2II_PT),NFROFR(ISYMPQ))
                        END IF
CSONIA SONIA SONIA
                        CALL CC_FD2BL(WORK(KFD2II),WORK(KFD2IJ),
     *                                WORK(KFD2JI),WORK(KFD2AI),
     *                                WORK(KFD2IA),WORK(KONEIJ),
     *                                WORK(KONEAB),WORK(KONEAI),
     *                                WORK(KONEIA),WORK(KCMOF),
     *                                WORK(KLAMDH),WORK(KLAMDP),
     *                                WORK(KEND3),LWRK3,IDEL,
     *                                ISYMD,G,ISYMG)
C
                        CALL CC_FD2AO(WORK(KAODEN),WORK(KFD2II),
     *                                WORK(KFD2IJ),WORK(KFD2JI),
     *                                WORK(KFD2AI),WORK(KFD2IA),
     *                                WORK(KCMOF),WORK(KLAMDH),
     *                                WORK(KLAMDP),WORK(KEND3),LWRK3,
     *                                ISYMPQ)
C
                        CALL CC_D2GAF(WORK(KD2GIJ),WORK(KD2GAB),
     *                                WORK(KD2GAI),WORK(KD2GIA),
     *                                WORK(KONEIJ),WORK(KONEAB),
     *                                WORK(KONEAI),WORK(KONEIA),
     *                                WORK(KCMOF),IDEL,ISYMD,G,ISYMG)
C
                        IF (CCPT) THEN
                           CALL CCPT_FD2BL(WORK(KFD2II_PT),
     &                                  WORK(KFD2IJ_PT),
     &                                  WORK(KFD2JI_PT),WORK(KFD2AI_PT),
     &                                  WORK(KFD2IA_PT),WORK(KONEIA_PT),
     &                                     WORK(KCMO),WORK(KCMOF),
     &                                     WORK(KEND3),LWRK3,
     &                                     IDEL,ISYMD,G,ISYMG)
                           CALL CC_FD2AO(WORK(KAODEN),WORK(KFD2II_PT),
     &                                WORK(KFD2IJ_PT),WORK(KFD2JI_PT),
     &                                WORK(KFD2AI_PT),WORK(KFD2IA_PT),
     &                                WORK(KCMOF),WORK(KCMO),
     &                                WORK(KCMO),WORK(KEND3),LWRK3,
     &                                ISYMPQ)
                           CALL CCPT_D2GAF(WORK(KD2GIA_PT),
     &                                WORK(KONEIA_PT),WORK(KCMOF),
     &                                IDEL,ISYMD,G,ISYMG)
                        END IF

                     ENDIF
C
C
C-----------------------------------------------
C                    Work space allocation four.
C-----------------------------------------------
C
                     KINTAO = KEND3
                     KEND4  = KINTAO + N2BST(ISYMPQ)
                     LWRK4  = LWORK  - KEND4
C
                     IF (LWRK4 .LT. 0) THEN
                       WRITE(LUPRI,*) 'Available:', LWORK
                       WRITE(LUPRI,*) 'Needed:', KEND4
                       CALL QUIT('Insufficient work space in CC_GRAD')
                     ENDIF
C
C-------------------------------------------------------
C                    Square up AO-integral distribution.
C-------------------------------------------------------
C
                     KOFFIN = KXINT + IDSAOG(ISYMG,ISYDIS)
     *                      + NNBST(ISYMPQ)*(G - 1)
C
                     CALL CCSD_SYMSQ(WORK(KOFFIN),ISYMPQ,WORK(KINTAO))
C
C------------------------------------------------------------
C                    Backtransform density fully to AO basis.
C------------------------------------------------------------
C
                     AUTIM1 = SECOND()
                     IF (CCSD .OR. CC2 .OR. CCPT) THEN
C
                        CALL CC_DENAO(WORK(KAODEN),ISYMPQ,
     *                                WORK(KD2GAI),WORK(KD2GAB),
     *                                WORK(KD2GIJ),WORK(KD2GIA),ISYMPQ,
     *                                WORK(KLAMDP),1,WORK(KLAMDH),1,
     *                                WORK(KEND4),LWRK4)
CSONIA SONIA
                        IF (CCPT) THEN
                           CALL CC_DENAO(WORK(KAODEN),ISYMPQ,
     &                                WORK(KD2GAI_PT),WORK(KD2GAB_PT),
     &                                WORK(KD2GIJ_PT),WORK(KD2GIA_PT),
     &                                ISYMPQ,
     &                                WORK(KCMO),1,WORK(KCMO),1,
     &                                WORK(KEND4),LWRK4)
                        END IF
CSONIA SONIA
C
                     ELSE
C
                        CALL CCMP_DAO(WORK(KAODEN),WORK(KD2GIJ),
     *                                WORK(KD2GIA),WORK(KCMO),
     *                                WORK(KLAMDH),WORK(KEND4),
     *                                LWRK4,ISYMPQ)
C
                     ENDIF
                     AUTIME = SECOND() - AUTIME
                     TIMDAO = TIMDAO + AUTIME
C
C--------------------------------------------------------
C                    Add relaxation terms to set up 
C                    effective density. We thus have the
C                    entire effective 2-electron density.
C--------------------------------------------------------
C
                     IF (.NOT. CCS) THEN
C
                       ICON = 2
C
                       CALL CC_D2EFF(WORK(KAODEN),G,ISYMG,IDEL,ISYMD,
     *                               WORK(KKABAO),WORK(KDHFAO),ICON)
C
                       CALL CC_D2EFF(WORK(KAODEN),G,ISYMG,IDEL,ISYMD,
     *                               WORK(KDHFAO),WORK(KKABAO),ICON)
C
                     ENDIF
C
C-----------------------------------------------------------
C                    calculate the 2 e- density contribution
C                    to the zeroth-order energy.
C-----------------------------------------------------------
C
                     ECCSD2 = ECCSD2 + HALF*DDOT(N2BST(ISYMPQ),
     *                        WORK(KAODEN),1,WORK(KINTAO),1)
C
C--------------------------------------------------------------
C                    calculate the 2 e- density contribution to
C                    the zeroth-order effective Fock matrix.
C--------------------------------------------------------------
C
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
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C------------------------------
C        Remove integral files.
C------------------------------
C
      IF (KEEPAOTWO.LT.2) THEN
         IF (LUINTA.LE.0) THEN
           LUINTA = -1
           CALL GPOPEN(LUINTA,'AOTWOINT','UNKNOWN',' ','UNFORMATTED',
     *                 IDUMMY,.FALSE.)
         END IF
         CALL GPCLOSE(LUINTA,'DELETE')
      END IF
C
      TIMREO = SECOND()
C
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
      IF ((CCSD .or. CCPT) .AND. (FROIMP)) THEN
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
      IF (IOPT .EQ. 1) THEN
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
              LABEL1 = '1DOVL'//CHRNOS(ISCOOR/100)
     &                       //CHRNOS(MOD(ISCOOR,100)/10)
     &                       //CHRNOS(MOD((MOD(ISCOOR,100)),10))
C
              CALL CCPRPAO(LABEL1,.TRUE.,WORK(K1DOVL),IRREP,ISYM,IERR,
     &                     WORK(KEND2),LWRK2)
C
              IF (IERR.NE.0 .OR. IRREP.NE.ICORSY) THEN
                 CALL QUIT('CC_DERIV: error while reading '
     &                        //LABEL1//' integrals.')
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
      ELSE IF (IOPT .EQ. 2) THEN
C
         K1DOVL = KEND1A
         KEND2  = K1DOVL + N2BST(1)
         LWRK2  = LWORK - KEND2
C
         IF (LWRK2 .LE. 0) THEN
            CALL QUIT('Insufficient work space in CC_GRAD.')
         END IF
C
         LABEL1 = 'KINENERG'
         CALL CCPRPAO(LABEL1,.TRUE.,WORK(K1DOVL),IRREP,ISYM,IERR,
     &                WORK(KEND2),LWRK2)
C
         IF (IERR.NE.0 .OR. IRREP.NE.1) THEN
            CALL QUIT('CC_DERIV: error while reading '
     &                   //LABEL1//' integrals.')
         ENDIF
C
         GRADFS(1) = -HALF*DDOT(N2BST(1),WORK(K1DOVL),1,
     *                         WORK(KFCKEF),1)
         NUCIND = NUSAIN
      ENDIF
C
C---------------------------------------------------
C     Read in potential energy and write out energy.
C---------------------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC) POTNUC
      CALL GPCLOSE(LUSIFC,'KEEP')
C
      ECCSD = ECCSD1 + ECCSD2 + POTNUC
C
      WRITE(LUPRI,*) ' '
      IF (IPRINT .GT. 10 .OR. LOCDBG) THEN
         WRITE(LUPRI,*) 'Coupled Cluster energy constructed'
         WRITE(LUPRI,*) 'from density matrices:'
         IF (CCS)  WRITE(LUPRI,*) 'CCS-energy:',  ECCSD
         IF (MP2)  WRITE(LUPRI,*) 'MP2-energy:',  ECCSD
         IF (CCD)  WRITE(LUPRI,*) 'CCD-energy:',  ECCSD
         IF (CCSD) WRITE(LUPRI,*) 'CCSD-energy:', ECCSD
         IF (CCPT) WRITE(LUPRI,*) 'CCSD(T)-energy:', ECCSD
         WRITE(LUPRI,*) 'H1 energy, ECCSD1 = ',ECCSD1
         WRITE(LUPRI,*) 'H2 energy, ECCSD2 = ',ECCSD2
         WRITE(LUPRI,*) 'Nuc. Pot. energy  = ',POTNUC
         WRITE(LUPRI,*) 'FTRACE            = ',FTRACE
C
         FNORM = DDOT(N2BST(1),WORK(KFCKEF),1,WORK(KFCKEF),1)
         WRITE(LUPRI,*) 
     &    'Norm^2 of the zeroth-order eff. Fock matrix:',FNORM
      ENDIF
C
      TIMREO = SECOND() - TIMREO
C
C
C-----------------------
C     Write out timings.
C-----------------------
C
  99  TIMTOT = SECOND() - TIMTOT
C
      IF (IPRINT .GT. 3) THEN
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'One electron and reorthonormalization'//
     *              ' gradient calculation completed'
         WRITE(LUPRI,*) 'Total time used in CC_GRAD:', TIMTOT
      ENDIF
      IF (IPRINT .GT. 9) THEN
       WRITE(LUPRI,*) 'Time used for setting up d(pq,ga,de):', TIMDEN
       WRITE(LUPRI,*) 'Time used for full AO backtransformation:',TIMDAO
       WRITE(LUPRI,*) 'Time used for reading 2 e- AO-integrals:', TIRDAO
       WRITE(LUPRI,*) 'Time used for calculating 2 e- AO-integrals:',
     *            TIMHE2
       WRITE(LUPRI,*) 'Time used for 1 e- part & intermediates:',
     *            TIMONE
       WRITE(LUPRI,*) 'Time used for reorthonormalization part:', TIMREO
      ENDIF
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('CC_GRAD_1')

       IF (CCPT) THEN
          CCSD = LCCSD_SAV
       END IF

      RETURN
      END
