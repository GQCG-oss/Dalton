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
C  /* Deck ccefffock */
*======================================================================*
      SUBROUTINE CCEFFFOCK(I1DXORD,IOPREL,WORK,LWORK)
*----------------------------------------------------------------------*
C
C     Written by Christof Haettig march 2000, based on CC_2EEXP.
C
C     Purpose: To calculate expectation values of one- and two-electron
C              operators and effective Fock matrices using the
C              one- and two-electron density matrices
C
C              The information about which expectation values and
C              which Fock matrices are to be calculated is read
C              from the common blocks in ccexpfck.h and cc1dxfck.h.
C              Expectation values are returned on the common block,
C              the eff. Fock matrices are written to a direct access
C              file and only the start addresses are returned on
C              common blocks in ccexpfck.h and cc1dxfck.h
C
C              I1DXORD = 0  --  calculate usual effective Fock matrices
C                               and expectation values
C
C              I1DXORD = 1  --  calculate effective Fock matrices
C                               from one-index transformed densities
C
C     Current models: CCS, MP2, CCD, CCSD
C
C     CC2 (without frozen core) by A. Halkier & S. Coriani 20/01-2000.
C
*----------------------------------------------------------------------*
      call quit('CCEFFFOCK is not working')
      return
      end
#ifdef HJAAJ_REMOVED
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
      INTEGER  LUBAR0
#else
#include <implicit.h>
#endif
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
#include <ccroper.h>
#include <ccropr2.h>
#include <ccexpfck.h>
#include <cc1dxfck.h>
#include <ccr1rsp.h>
#include <chrxyz.h>
#include <inftap.h>
#include <ccisao.h>
#include <r12int.h>
C
      LOGICAL LOCDBG ! debug flag
      PARAMETER (LOCDBG = .FALSE.)
      LOGICAL DIRINT ! direct derivative integrals
      LOGICAL OLDDX  ! flag to see if file is open for 'direct' read
      PARAMETER (DIRINT = .FALSE.)
C
      INTEGER LUFCK, LUAOIN, LUDE, ISYM0
      PARAMETER (ISYM0 = 1)
      LOGICAL SAVDIR, LEX, SAVHER
      INTEGER INDEXA(MXCORB)
      CHARACTER*8 LABEL, LABELH, FILFCK
      CHARACTER*10 MODEL
      CHARACTER*3 LSTRLX

      INTEGER LWORK

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL RE2DAR, RELCO1, ECCSD
      REAL TIMTOT, TIMDEN, TIMDAO, TIRDAO, TIMHE2, TIMONE
      REAL ECCSD1, DDOT, POTNUC, AUTIM1, SECOND, DTIME, AUTIME
      REAL ZERO, HALF, ONE, TWO
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION RE2DAR, RELCO1, ECCSD
      DOUBLE PRECISION TIMTOT, TIMDEN, TIMDAO, TIRDAO, TIMHE2, TIMONE
      DOUBLE PRECISION ECCSD1, DDOT, POTNUC, AUTIM1, SECOND, DTIME
      DOUBLE PRECISION AUTIME
      DOUBLE PRECISION ZERO, HALF, ONE, TWO
#endif
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)

      LOGICAL LDOTHISATOM, LANYTHING, LEXPEC, LFOCK, LTWO, SQRINT
      LOGICAL LONLYEXP, SYM1ONLY
      LOGICAL DO2DAR, AD2DAR, S4CENT, DPTINT
      LOGICAL LDERINT(8,3)
      INTEGER I1DXORD, IOPREL, NGLMDTM, NALLAIM
      INTEGER IDX, ISYMTR, ISYMOP, ISYOPE, ISYDIS, ISYRLX, ISYPQ1
      INTEGER ISCOOR, ICORSY, ICOOR, ISYMG, ISYMPQ, IGAM, IDLST1
      INTEGER NTOSYM, ISYMD1, ILLL, NUMDIS, KRECNR, LDECH, IDEL2, IDEL
      INTEGER ISYMD, ISYDEN, IDLST, ISYM, LENBAR, KEND0, LWRK0
      INTEGER KTEST1, KENDTS, IATOM, NFOCK1LBL, JATOM, NTOT
      INTEGER N2BSTM, KFCKEF, KAODEN, KCMO, KT2AM, KZ2AM, KLAMDP, KLAMDH
      INTEGER KDNS1D, KRMAT, KKAPPA, KCMOPQ, KCMOHQ, KB1DHFAO, KB2DHFAO
      INTEGER KLAMDPQ, KLAMDHQ, KOVERLP, KQAOS, KB1KABAO, KB2KABAO
      INTEGER KLAMDPQ2, KLAMDHQ2, KT1AM, KZ1AM, KEND1, LWRK1, KT1AMT
      INTEGER KONEAI, KONEAB, KONEIJ, KONEIA, KKABAR, KDHFAO, KINTIJ
      INTEGER KEND1, LWRK1, KOFDIJ, KOFDAB, KOFDAI, KOFDIA, KT2AMT
      INTEGER KXMAT, KYMAT, KMINT, KONINT, KMIRES, KD1ABT, KD1IJT
      INTEGER KCCFB1, KINDXB, KODCL1, KODCL2, KODBC1, KODBC2, KRDBC1
      INTEGER KRDBC2, KODPP1, KODPP2, KRDPP1, KRDPP2, KFREE, LFREE
      INTEGER KD2IJG,   KD2AIG,   KD2IAG,   KD2ABG, KCMOF, KWRKTS
      INTEGER KDB2IJG,  KDB2AIG,  KDB2IAG,  KDB2ABG, LWRKTS
      INTEGER KDB22IJG, KDB22AIG, KDB22IAG, KDB22ABG, KKABAO, KOFFIN
      INTEGER KOFFAI, KOFFIA, KOFFIJ, KOFFJI, KXINT, KINTAO, KCHE3
      INTEGER KFD2IJ, KFD2JI, KFD2AI, KFD2IA, KFD2II, KENDSV, LWRKSV
      INTEGER KD2GIJ,   KD2GAI,   KD2GAB,   KD2GIA, NDAD, NDENEL
      INTEGER KDB2GIJ,  KDB2GAI,  KDB2GAB,  KDB2GIA, KEND3, LWRK3
      INTEGER KDB22GIJ, KDB22GAI, KDB22GAB, KDB22GIA, KEND4, LWRK4
      INTEGER KEND2, LWRK2, NFOCK0LBL, IDLST0, ISYMOPR, ISYFCK, MXCOMP
      INTEGER NVEDUM, MODDUM, NFRL, NII, ICON, ISDEN, IOPER, IORDER
      INTEGER NFOCKLBL, IROPER, IDXRLX

* external functions:
      INTEGER IR1KAPPA, ILSTSYM
C
      CALL QENTER('CCEFFFOCK')
C
C------------------------------
C     check if anything to do:
C------------------------------
C
      LANYTHING = .FALSE.
      LFOCK     = .FALSE.
      LONLYEXP  = .FALSE.

      IF (IOPREL.GT.0) THEN
         LANYTHING = .TRUE.
         LONLYEXP  = .TRUE.
      END IF

      IF (.NOT.LONLYEXP) THEN
        IF (I1DXORD.EQ.0) THEN
          DO IDLST = 1, NEXPFCKLBL
            IF (LEXPFCK(1,IDLST).OR.LEXPFCK(2,IDLST)) LANYTHING=.TRUE.
            IF (LEXPFCK(2,IDLST)) LFOCK = .TRUE.
          END DO
        ELSE IF (I1DXORD.EQ.1) THEN
          IF (N1DXFLBL.GT.0) LANYTHING = .TRUE.
          IF (N1DXFLBL.GT.0) LFOCK     = .TRUE.
        END IF
      END IF

      IF (.NOT. LANYTHING) GOTO 999

      IF (LFOCK .AND. (.NOT.RELORB)) THEN
        WRITE(LUPRI,*) 
     *     'Calculation of eff. Fock mat. needs RELORB = .TRUE.'
        CALL QUIT(
     *     'Calculation of eff. Fock mat. needs RELORB = .TRUE.')
      END IF

      IF (LOCDBG) THEN
         WRITE(LUPRI,*) 'Entering CCEFFFOCK with I1DXORD = ',I1DXORD
      END IF
C
C---------------------------------------------------------
C     Initialization of result for the relativistic stuff.
C---------------------------------------------------------
C
      IF (IPRINT.GT.9 .OR. LOCDBG) CALL AROUND('Entering CCEFFFOCK')
      CALL FLSHFO(LUPRI)
      RE2DAR = ZERO
      IF (IOPREL .EQ. 1) RELCO1 = WORK(1)
C
C------------------------------------------------------------
C     for I1DXORD = 1, initialize list of relaxation vectors:
C------------------------------------------------------------
C
      IF (.NOT.LONLYEXP) THEN
        DO IDX = 1, N1DXFLBL
          IF (LST1DXFCK(IDX)(1:2) .EQ. 'R1') THEN
            IRELAX1DX(IDX) = 
     &          IR1KAPPA(LABRLX(IDX),FRQRLX(IDX),ISYMRLX(IDX))
          ELSE
            CALL QUIT('Unknown List in CCEFFFOCK.')
          END IF
        END DO
      END IF
C
C-----------------------------------------
C     Initialization of timing parameters.
C-----------------------------------------
C
      TIMTOT = SECOND()
      TIMDEN = ZERO
      TIMDAO = ZERO
      TIRDAO = ZERO
      TIMHE2 = ZERO
      TIMONE = SECOND()
C
C----------------------------------------------------
C     Both zeta- and t-vectors are totally symmetric.
C----------------------------------------------------
C
      ISYMTR = 1
      ISYMOP = 1
C
C----------------------------------------------------
C     do some other initializations:
C----------------------------------------------------
C
      N2BSTM  = 0
      NALLAIM = 0
      NGLMDTM = 0
      DO ISYM = 1, NSYM
        N2BSTM  = MAX(N2BSTM,N2BST(ISYM))
        NALLAIM = MAX(NALLAIM,NALLAI(ISYM))
        NGLMDTM = MAX(NGLMDTM,NGLMDT(ISYM))
      END DO

      LENBAR = 2*NT1AMX+NMATIJ(1)+NMATAB(1)+2*NT1FRO(1)+2*NCOFRO(1)

      KEND0 = 1
C
C---------------------------------------------------------------------
C     allocate work space for intermediates needed for the calculation
C     of the effective Fock matrices:
C---------------------------------------------------------------------
C
      IF (.NOT.LONLYEXP) THEN
        KDNS1D   = KEND0                ! 1-index transformed density
        KRMAT    = KDNS1D   + N2BSTM    ! connection matrix
        KKAPPA   = KRMAT    + N2BSTM    ! orbital relaxation vector
        KCMOPQ   = KKAPPA   + 2*NALLAIM ! derivative MO vector (ket)
        KCMOHQ   = KCMOPQ   + NGLMDTM   ! derivative MO vector (bra)
        KB1DHFAO = KCMOHQ   + NGLMDTM   ! alpha-idx transf. HF den
        KB2DHFAO = KB1DHFAO + N2BSTM    ! beta-idx transf. HF den
        KLAMDPQ  = KB2DHFAO + N2BSTM    ! derivative Lambda^p
        KLAMDHQ  = KLAMDPQ  + NGLMDTM   ! derivative Lambda^h
        KEND0    = KLAMDHQ  + NGLMDTM  
C
        IF (I1DXORD.GT.0) THEN
         KOVERLP  = KEND0               ! overlap matrix
         KQAOS    = KOVERLP  + N2BST(1) ! Q^ao x S^AO
         KB1KABAO = KQAOS    + N2BSTM   ! alpha-idx transf. relax. con.
         KB2KABAO = KB1KABAO + N2BSTM   ! beta-idx transf. relax. con.
         KLAMDPQ2 = KB2KABAO + N2BSTM   ! derivative lambda matrices
         KLAMDHQ2 = KLAMDPQ2 + NGLMDTM  ! with left & right Q exchanged
         KEND0    = KLAMDHQ2 + NGLMDTM  
        END IF
      END IF

      LWRK0 = LWORK - KEND0

C
C------------------------------------------------------------------
C     if we are going to calculate one-index transformed effective
C     Fock matrices, we read here the overlap matrix:
C------------------------------------------------------------------
C
      IF (I1DXORD.EQ.1 .AND. (.NOT.LONLYEXP)) THEN
         IF (LWRK0.LT.NBAST) 
     &     CALL QUIT('Insufficient work space in CCEFFFOCK.')
         CALL RDONEL('OVERLAP ',.TRUE.,WORK(KEND1),NBAST)
         CALL CCSD_SYMSQ(WORK(KEND1),ISYM0,WORK(KOVERLP))
      END IF 

*======================================================================*
* do model specific initializations for calculating the densities:
*======================================================================*
C
      IF (CC2) THEN
C
C
C-----------------------------------
C     Initial work space allocation.
C-----------------------------------
C
         KFCKEF = KEND0
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
     *     'Insufficient core for initial allocation in CCEFFFOCK')
         ENDIF
C
C-------------------------------------------------------------
C        Read MO-coefficients from interface file and reorder.
C-------------------------------------------------------------
C
         IF (LUSIFC .LE. 0) THEN
            LUSIFC = -1
            CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     *                  IDUMMY,.FALSE.)
         ENDIF
C
         REWIND LUSIFC
         CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
         READ (LUSIFC)
         READ (LUSIFC)
         READ (LUSIFC) (WORK(KCMO+I-1), I=1,NLAMDS)
         CALL GPCLOSE (LUSIFC,'KEEP')
C
         CALL CMO_REORDER(WORK(KCMO),WORK(KEND1),LWRK1)
C
C-------------------------------------------
C        Read zero'th order zeta amplitudes.
C-------------------------------------------
C
         IOPT   = 3
         CALL CC_RDRSP('L0',0,ISYM0,IOPT,MODEL,WORK(KSKOD),WORK(KZ2AM))
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
         CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KT1AM),WORK(KT2AM))
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
            CALL QUIT(
     *         'Insufficient core for allocation 1 in CCEFFFOCK')
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
         LUBAR0 = -876
         CALL GPOPEN(LUBAR0,'CCKABAR0','OLD',' ','UNFORMATTED',
     *               IDUMMY,.FALSE.)
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
      ELSE IF (CCSD) THEN
C
C-----------------------------------
C     Initial work space allocation.
C-----------------------------------
C
         KFCKEF = KEND0
         KAODEN = KFCKEF + N2BST(1)
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
     *       'Insufficient core for first allocation in CCEFFFOCK')
         ENDIF
C
C----------------------------------------
C     Read zero'th order zeta amplitudes.
C----------------------------------------
C
         IOPT   = 3
         CALL CC_RDRSP('L0',0,ISYM0,IOPT,MODEL,WORK(KSKOD),WORK(KZ2AM))
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
         CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KT1AM),WORK(KT2AM))
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
     *           'Insufficient memory for allocation 1 CCEFFFOCK')
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
         CALL CC_GET_CMO(WORK(KCMO))
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
         KFCKEF = KEND0
         KAODEN = KFCKEF + N2BST(1)
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
     *    'Insufficient memory for work allocation in CCEFFFOCK')
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
         KT1AM = KEND1
         KT2AM = KT1AM + NT1AMX
         KZ2AM = KT2AM + NT2AMX
         KZ1AM = KZ2AM + NT2AMX
         KEND1 = KZ1AM + NT1AMX
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT(
     *    'Insufficient memory for work allocation in CCEFFFOCK')
         ENDIF
C
C----------------------------------------
C     Read zero'th order zeta amplitudes.
C----------------------------------------
C
         IOPT   = 3
         CALL CC_RDRSP('L0',0,ISYM0,IOPT,MODEL,WORK(KSKOD),WORK(KZ2AM))
C
         CALL DZERO(WORK(KZ1AM),NT1AMX)
C
C-------------------------------------------
C     Read zero'th order cluster amplitudes.
C-------------------------------------------
C
         IOPT = 3
         CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KT1AM),WORK(KT2AM))
C
         CALL DZERO(WORK(KT1AM),NT1AMX)
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
         LWRK1 = LWORK - KEND1
C
      ELSE IF (CCS) THEN
C
C---------------------------------
C     First work space allocation.
C---------------------------------
C
         KFCKEF = KEND0
         KAODEN = KFCKEF + N2BST(1)
         KCMO   = KAODEN + N2BSTM
         KDHFAO = KCMO   + NLAMDS
         KT1AM  = KDHFAO + N2BST(1)
         KKABAR = KT1AM  + NT1AM(1)
         KEND1  = KKABAR + LENBAR
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT(
     *    'Insufficient memory for work allocation in CCEFFFOCK')
         ENDIF
C
         CALL DZERO(WORK(KT1AM),NT1AM(1))
         CALL DZERO(WORK(KKABAR),LENBAR)
C
         CALL CCS_D1AO(WORK(KAODEN),WORK(KEND1),LWRK1)
         CALL CCS_D1AO(WORK(KDHFAO),WORK(KEND1),LWRK1)
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
C-----------------------------------------------------------------------
C     loop over list of Fock matrices / expectation values and calculate
C     the 1el contributions to the exp. values and eff. Fock matrices:
C-----------------------------------------------------------------------
C
      IF (.NOT.LONLYEXP) THEN
        CALL CCEFFFOCK1(I1DXORD,WORK(KAODEN),WORK(KDHFAO),
     *                  WORK(KEND1),LWRK1)
      END IF
C
C------------------------------------------------------
C     if required, set flags for relativistic integrals
C------------------------------------------------------
C
      IF (LONLYEXP) THEN
        SAVDIR = DIRECT
        SAVHER = HERDIR
        DIRECT = .TRUE.
        HERDIR = .TRUE.
        INQUIRE(UNIT=93,EXIST=LEX)
        IF (LEX) CLOSE(93,STATUS='DELETE')             
        INQUIRE(UNIT=42,EXIST=LEX)
        IF (LEX) CLOSE(42,STATUS='DELETE')             
        IF (IOPREL .EQ. 2) THEN
C          DPTINT = .TRUE.
           BPH2OO = .TRUE.
        ENDIF
        IF (DAR2EL) THEN
           DO2DAR = .TRUE.
           AD2DAR = .FALSE.
           S4CENT = .FALSE.
        ENDIF
      ELSE
        DIRECT = .FALSE.
      END IF
C
C----------------------------------------------------------------------
C     open file for effective Fock matrices:
C----------------------------------------------------------------------
C
      IF (.NOT.LONLYEXP) THEN
        IF (I1DXORD.EQ.0) THEN
          FILFCK = FILFCKEFF
        ELSE IF (I1DXORD.EQ.1) THEN
          FILFCK = FIL1DXFCK
        ELSE
          CALL QUIT('Illegal value of I1DXORD in CCEFFFOCK.')
        END IF

        LUFCK = -1
        CALL WOPEN2(LUFCK,FILFCK,64,0)
      END IF
C
C----------------------------------------------------------------------
C     start loop over atoms:
C----------------------------------------------------------------------
C
      DO IATOM = 0, NUCIND

        IF (LONLYEXP) THEN
          LDOTHISATOM = (IATOM .EQ. 0) 
        ELSE
          CALL CCEFFFOCKHLP1(IATOM,I1DXORD,LDOTHISATOM,LABELH)
        END IF

        IF (LDOTHISATOM) THEN

          ! since we don't have yet direct derivative integrals,
          ! we precalculate the complete set of integrals here
          ! calling HERMIT via CCDER1.
          IF (LABELH(1:5).EQ.'1DHAM' .OR. LABELH(1:5).EQ.'dh/dB') THEN 
             ! set dorps.h common blocks...
             ! (this needs still to be adapted for symmetry...)
             SYM1ONLY = .TRUE.
             CALL CC_SETDORPS(LABELH,SYM1ONLY,0)
             ! precalculate and sort derivative integrals for this atom
             CALL CCDER1(IATOM,LABELH,LDERINT,WORK(KEND1),LWRK1)
          END IF

        ELSE

          GOTO 90 ! skip everything for this atom

        END IF
C
C----------------------------------------------------------------------
C     for 1-index transformed fock matrices start here a loop over the 
C     fock matrices, since for each matrix we need to compute 1-index
C     transformed densities: 
C----------------------------------------------------------------------
C
      IF (LONLYEXP .OR. I1DXORD.EQ.0) THEN
        NFOCK1LBL = 1
      ELSE IF (I1DXORD.EQ.1) THEN
        NFOCK1LBL = N1DXFLBL
      END IF

      DO IDLST1 = 1, NFOCK1LBL

        IF (I1DXORD.EQ.1 .AND. (.NOT.LONLYEXP)) THEN
           IDLST = IDLST1

           CALL CCEFFFOCKHLP2(IDLST,JATOM,IATOM,ISYMOPR,ISYRLX,ISYFCK,
     &                        ISCOOR,ICOOR,ICORSY,MXCOMP,
     &                        LEXPEC,LFOCK,LTWO,LABEL,IOPER,
     &                        WORK(KKAPPA),WORK(KRMAT),WORK(KQAOS),
     &                        WORK(KT1AM),WORK(KOVERLP),
     &                        WORK(KLAMDPQ),WORK(KLAMDHQ),
     &                        WORK(KCMOPQ),WORK(KCMOHQ),
     &                        WORK(KB1DHFAO),WORK(KB2DHFAO), 
     &                        WORK(KB1KABAO),WORK(KB2KABAO), 
     &                        WORK(KLAMDPQ2),WORK(KLAMDHQ2),
     &                        WORK(KCMO),WORK(KKABAR),
     &                        WORK(KEND1),LWRK1)
        END IF

        IF (LONLYEXP .OR. I1DXORD.EQ.0 .OR. (JATOM.EQ.IATOM .AND. LTWO))
     &  THEN
C
C-----------------------------------
C     Start the loop over integrals.
C-----------------------------------
C
      IF (DIRECT) THEN
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
     *                 KODPP1,KODPP2,KRDPP1,KRDPP2,
     *                 KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     *                 WORK(KEND1),LWRK1,IPRERI)
           KEND1  = KFREE
           LWRK1  = LFREE
        ENDIF
        DTIME  = SECOND() - DTIME
        TIMHE2 = TIMHE2 + DTIME
        NTOSYM = 1
C
      ELSE
        NTOSYM = NSYM
      END IF
C
      KENDSV = KEND1
      LWRKSV = LWRK1
C
      DO ISYMD1 = 1, NTOSYM
C
        IF (DIRECT) THEN
          IF (HERDIR) THEN
             NTOT = MAXSHL
          ELSE
             NTOT = MXCALL
          ENDIF
        ELSE
          NTOT = NBAS(ISYMD1)
        END IF
C
      DO 100 ILLL = 1,NTOT
C
C---------------------------------------------------------------
C        Determine which delta's to be calculated in this round.
C---------------------------------------------------------------
C
         IF (DIRECT) THEN
C
           KEND1 = KENDSV
           LWRK1 = LWRKSV
C
           DTIME  = SECOND()
           IF (HERDIR) THEN
              CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,
     &                    IPRERI)
           ELSE
              CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     *                    WORK(KODCL1),WORK(KODCL2),
     *                    WORK(KODBC1),WORK(KODBC2),
     *                    WORK(KRDBC1),WORK(KRDBC2),
     *                    WORK(KODPP1),WORK(KODPP2),
     *                    WORK(KRDPP1),WORK(KRDPP2),
     *                    WORK(KCCFB1),WORK(KINDXB),
     *                    WORK(KEND1), LWRK1,IPRERI)
           ENDIF
           DTIME  = SECOND() - DTIME
           TIMHE2 = TIMHE2 + DTIME
C
           KRECNR = KEND1
           KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
           LWRK1  = LWORK  - KEND1
           IF (LWRK1 .LT. 0) THEN
              CALL QUIT('Insufficient core in CCEFFFOCK')
           END IF
C
         ELSE
           NUMDIS = 1
         END IF
C
C-------------------------------------------------------
C        Open file for efffective two electron densities.
C-------------------------------------------------------
C
#if defined (SYS_IRIX) || defined (SYS_NEC) || defined (SYS_DEC)
         NFRL = 2
#endif
#if defined (SYS_LINUX) || defined (SYS_AIX) || defined (SYS_PARAGON) || defined (SYS_ALLIANT) \
|| defined (SYS_SUN) || defined (SYS_HPUX) || defined (SYS_HAL)
         NFRL = 8
#endif
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
         NFRL = 16
#endif
C
C
         LDECH = N2BSTM*NFRL+1
         LUDE = -1
         CALL GPOPEN(LUDE,'CCTWODEN','UNKNOWN','DIRECT','UNFORMATTED',
     *               IDUMMY,OLDDX)
C
C------------------------------------------------
C        Loop over number of delta distributions.
C------------------------------------------------
C
         DO 110 IDEL2 = 1,NUMDIS
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
            END IF
C
C-------------------------------------
C           Work space allocation two.
C-------------------------------------
C
            ISYDEN = ISYMD
            ISYDNQ = MULD2H(ISYDEN,ISYRLX)
C
            IF (CCSD .OR. CC2) THEN
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

            IF (I1DXORD.EQ.1) THEN
               IF (CCSD .OR. CC2) THEN
                  ! for lambda^q transformed densities
                  KDB2IJG = KEND2
                  KDB2AIG = KDB2IJG + ND2IJG(ISYDNQ)
                  KDB2IAG = KDB2AIG + ND2AIG(ISYDNQ)
                  KDB2ABG = KDB2IAG + ND2AIG(ISYDNQ)
                  KEND2   = KDB2ABG + ND2ABG(ISYDNQ)

                  ! for lambda^q2 transformed densities
                  KDB22IJG = KEND2
                  KDB22AIG = KDB22IJG + ND2IJG(ISYDNQ)
                  KDB22IAG = KDB22AIG + ND2AIG(ISYDNQ)
                  KDB22ABG = KDB22IAG + ND2AIG(ISYDNQ)
                  KEND2    = KDB22ABG + ND2ABG(ISYDNQ)
               ELSE IF (MP2) THEN
                  ! for lambda^q transformed densities
                  KDB2IJG = KEND2
                  KDB2IAG = KDB2IJG + NF2IJG(ISYDNQ)
                  KEND2   = KDB2IAG + ND2AIG(ISYDNQ)

                  ! for lambda^q2 transformed densities
                  KDB22IJG = KEND2
                  KDB22IAG = KDB22IJG + NF2IJG(ISYDNQ)
                  KEND2    = KDB22IAG + ND2AIG(ISYDNQ)
               ELSE IF (CCS) THEN
                  ! for lambda^q transformed densities
                  KDB2IJG = KEND2
                  KEND2   = KDB2IJG + NF2IJG(ISYDNQ)

                  ! for lambda^q2 transformed densities
                  KDB22IJG = KEND2
                  KEND2    = KDB22IJG + NF2IJG(ISYDNQ)
               ENDIF
            END IF
C
            LWRK2  = LWORK  - KEND2
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:',KEND2
               CALL QUIT(
     *           'Insufficient core for allocation 2 in CCEFFFOCK')
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

            IF (I1DXORD.EQ.1) THEN
              CALL DZERO(WORK(KDB2IJG),NF2IJG(ISYDNQ))
              IF (.NOT. CCS) THEN
                 CALL DZERO(WORK(KDB2IAG),ND2AIG(ISYDNQ))
                 IF (CCSD .OR. CC2) THEN
                    CALL DZERO(WORK(KDB2AIG),ND2AIG(ISYDNQ))
                    CALL DZERO(WORK(KDB2ABG),ND2ABG(ISYDNQ))
                    CALL DZERO(WORK(KDB2IJG),ND2IJG(ISYDNQ))
                 ENDIF
              ENDIF

              CALL DZERO(WORK(KDB22IJG),NF2IJG(ISYDNQ))
              IF (.NOT. CCS) THEN
                 CALL DZERO(WORK(KDB22IAG),ND2AIG(ISYDNQ))
                 IF (CCSD .OR. CC2) THEN
                    CALL DZERO(WORK(KDB22AIG),ND2AIG(ISYDNQ))
                    CALL DZERO(WORK(KDB22ABG),ND2ABG(ISYDNQ))
                    CALL DZERO(WORK(KDB22IJG),ND2IJG(ISYDNQ))
                 ENDIF
              ENDIF
            END IF
C
C----------------------------------------------------------------
C           Calculate the two electron density d(pq,gamma;delta).
C----------------------------------------------------------------
C
            IF (CCSD) THEN
              CALL CC_DEN2(WORK(KD2IJG),WORK(KD2AIG),WORK(KD2IAG),
     *                     WORK(KD2ABG),WORK(KZ2AM),WORK(KT2AM),
     *                     WORK(KT2AMT),WORK(KMINT),WORK(KXMAT),
     *                     WORK(KYMAT),WORK(KONEAB),WORK(KONEAI),
     *                     WORK(KONEIA),WORK(KMIRES),
     *                     WORK(KLAMDH),1,WORK(KLAMDP),1,
     *                     WORK(KEND2),LWRK2,IDEL,ISYMD)
              IF (I1DXORD.EQ.1) THEN
                CALL CC_DEN2(WORK(KDB2IJG),WORK(KDB2AIG),
     *                       WORK(KDB2IAG),WORK(KDB2ABG),
     *                       WORK(KZ2AM),WORK(KT2AM),WORK(KT2AMT),
     *                       WORK(KMINT),WORK(KXMAT),WORK(KYMAT),
     *                       WORK(KONEAB),WORK(KONEAI),WORK(KONEIA),
     *                       WORK(KMIRES),
     *                       WORK(KLAMDHQ),ISYRLX,WORK(KLAMDP),1,
     *                       WORK(KEND2),LWRK2,IDEL,ISYMD)
                CALL CC_DEN2(WORK(KDB2IJG),WORK(KDB2AIG),
     *                       WORK(KDB2IAG),WORK(KDB2ABG),
     *                       WORK(KZ2AM),WORK(KT2AM),WORK(KT2AMT),
     *                       WORK(KMINT),WORK(KXMAT),WORK(KYMAT),
     *                       WORK(KONEAB),WORK(KONEAI),WORK(KONEIA),
     *                       WORK(KMIRES),
     *                       WORK(KLAMDH),1,WORK(KLAMDPQ),ISYRLX,
     *                       WORK(KEND2),LWRK2,IDEL,ISYMD)
                CALL CC_DEN2(WORK(KDB22IJG),WORK(KDB22AIG),
     *                       WORK(KDB22IAG),WORK(KDB22ABG),
     *                       WORK(KZ2AM),WORK(KT2AM),WORK(KT2AMT),
     *                       WORK(KMINT),WORK(KXMAT),WORK(KYMAT),
     *                       WORK(KONEAB),WORK(KONEAI),WORK(KONEIA),
     *                       WORK(KMIRES),
     *                       WORK(KLAMDHQ2),ISYRLX,WORK(KLAMDP),1,
     *                       WORK(KEND2),LWRK2,IDEL,ISYMD)
                CALL CC_DEN2(WORK(KDB22IJG),WORK(KDB22AIG),
     *                       WORK(KDB22IAG),WORK(KDB22ABG),
     *                       WORK(KZ2AM),WORK(KT2AM),WORK(KT2AMT),
     *                       WORK(KMINT),WORK(KXMAT),WORK(KYMAT),
     *                       WORK(KONEAB),WORK(KONEAI),WORK(KONEIA),
     *                       WORK(KMIRES),
     *                       WORK(KLAMDH),1,WORK(KLAMDPQ2),ISYRLX,
     *                       WORK(KEND2),LWRK2,IDEL,ISYMD)
              END IF
            ELSE IF (CC2) THEN
              CALL CC_DEN2(WORK(KD2IJG),WORK(KD2AIG),WORK(KD2IAG),
     *                     WORK(KD2ABG),WORK(KZ2AM),WORK(KT2AM),
     *                     WORK(KT2AMT),WORK(KEND2),WORK(KEND2),
     *                     WORK(KEND2),WORK(KONEAB),WORK(KONEAI),
     *                     WORK(KONEIA),WORK(KEND2),
     *                     WORK(KLAMDH),1,WORK(KLAMDP),1,
     *                     WORK(KEND2),LWRK2,IDEL,ISYMD)
              IF (I1DXORD.EQ.1) THEN
                CALL CC_DEN2(WORK(KDB2IJG),WORK(KDB2AIG),
     *                       WORK(KDB2IAG),WORK(KDB2ABG),
     *                       WORK(KZ2AM),WORK(KT2AM),WORK(KT2AMT),
     *                       WORK(KEND2),WORK(KEND2),WORK(KEND2),
     *                       WORK(KONEAB),WORK(KONEAI),WORK(KONEIA),
     *                       WORK(KEND2),
     *                       WORK(KLAMDHQ),ISYRLX,WORK(KLAMDP),1,
     *                       WORK(KEND2),LWRK2,IDEL,ISYMD)
                CALL CC_DEN2(WORK(KDB2IJG),WORK(KDB2AIG),
     *                       WORK(KDB2IAG),WORK(KDB2ABG),
     *                       WORK(KZ2AM),WORK(KT2AM),WORK(KT2AMT),
     *                       WORK(KEND2),WORK(KEND2),WORK(KEND2),
     *                       WORK(KONEAB),WORK(KONEAI),WORK(KONEIA),
     *                       WORK(KEND2),
     *                       WORK(KLAMDH),1,WORK(KLAMDPQ),ISYRLX,
     *                       WORK(KEND2),LWRK2,IDEL,ISYMD)
                CALL CC_DEN2(WORK(KDB22IJG),WORK(KDB22AIG),
     *                       WORK(KDB22IAG),WORK(KDB22ABG),
     *                       WORK(KZ2AM),WORK(KT2AM),WORK(KT2AMT),
     *                       WORK(KEND2),WORK(KEND2),WORK(KEND2),
     *                       WORK(KONEAB),WORK(KONEAI),WORK(KONEIA),
     *                       WORK(KEND2),
     *                       WORK(KLAMDHQ2),ISYRLX,WORK(KLAMDP),1,
     *                       WORK(KEND2),LWRK2,IDEL,ISYMD)
                CALL CC_DEN2(WORK(KDB22IJG),WORK(KDB22AIG),
     *                       WORK(KDB22IAG),WORK(KDB22ABG),
     *                       WORK(KZ2AM),WORK(KT2AM),WORK(KT2AMT),
     *                       WORK(KEND2),WORK(KEND2),WORK(KEND2),
     *                       WORK(KONEAB),WORK(KONEAI),WORK(KONEIA),
     *                       WORK(KEND2),
     *                       WORK(KLAMDH),1,WORK(KLAMDPQ2),ISYRLX,
     *                       WORK(KEND2),LWRK2,IDEL,ISYMD)
              END IF
            ELSE IF (MP2) THEN
              CALL CCS_DEN2(WORK(KD2IJG),WORK(KCMO),WORK(KEND2),
     *                      LWRK2,IDEL,ISYMD)
              CALL MP2_DEN2(WORK(KD2IAG),WORK(KT2AM),WORK(KLAMDH),
     *                      WORK(KEND2),LWRK2,IDEL,ISYMD)
              IF (I1DXORD.EQ.1) THEN
                CALL QUIT(
     *             '1-index transf. eff. Fock for MP2 unfinished.')
              END IF 
            ELSE IF (CCS) THEN
              CALL CCS_DEN2(WORK(KD2IJG),WORK(KCMO),WORK(KEND2),
     *                      LWRK2,IDEL,ISYMD)
              IF (I1DXORD.EQ.1) THEN
                CALL QUIT(
     *             '1-index transf. eff. Fock for CCS unfinished.')
              END IF 
            ENDIF
            AUTIME = SECOND() - AUTIME
            TIMDEN = TIMDEN + AUTIME
C
C----------------------------------------------------------------------
C           for I1DXORD=0 loop here over the operators and find the 
C           ones which belong to this "atom":
C             - compute one-index connection matrix transf. SCF density
C             - find coordinate index
C----------------------------------------------------------------------
C
            IF ( (.NOT.LONLYEXP) .AND. (I1DXORD.EQ.0) ) THEN
              NFOCK0LBL = NEXPFCKLBL
            ELSE
              NFOCK0LBL = 1
            END IF

            DO IDLST0 = 1, NFOCK0LBL

              IF ( (.NOT.LONLYEXP) .AND. (I1DXORD.EQ.0) ) THEN
                IDLST = IDLST0

                CALL CCEFFFOCKHLP3(IDLST,IATOM,JATOM,ISYMOPR,ISYFCK,
     &                             LABEL,IOPER,IORDER,
     &                             MXCOMP,ISCOOR,ICORSY,ICOOR,
     &                             LEXPEC,LFOCK,LTWO,WORK(KDNS1D),
     &                             WORK(KEND2),LWRK2)
              END IF
              
              IF (JATOM.EQ.IATOM .AND. LTWO .AND.
     &            (LABEL(1:5).EQ.'HAM0    '.OR.LDERINT(ICORSY,ICOOR)) )
     &        THEN
C
C---------------------------------------------------
C           Start loop over second AO-index (gamma).
C---------------------------------------------------
C
            DO 120 ISYMG = 1, NSYM
               DO 130 G  = 1, NBAS(ISYMG)
C
                  IGAM   = G + IBAS(ISYMG)
C
                  ISYMPQ = MULD2H(ISYMG,ISYDEN)
                  ISYPQ1 = MULD2H(ISYMG,ISYDNQ)
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
                     IF (I1DXORD.EQ.1) THEN
                       KDB2GIJ  = KDB2IJG  + ID2IJG(ISYPQ1,ISYMG)
     *                          + NMATIJ(ISYPQ1)*(G - 1) 
                       KDB2GAI  = KDB2AIG  + ID2AIG(ISYPQ1,ISYMG)
     *                          + NT1AM(ISYPQ1)*(G - 1)
                       KDB2GAB  = KDB2ABG  + ID2ABG(ISYPQ1,ISYMG)
     *                          + NMATAB(ISYPQ1)*(G - 1)
                       KDB2GIA  = KDB2IAG  + ID2AIG(ISYPQ1,ISYMG)
     *                          + NT1AM(ISYPQ1)*(G - 1)
                       KDB22GIJ = KDB22IJG + ID2IJG(ISYPQ1,ISYMG)
     *                          + NMATIJ(ISYPQ1)*(G - 1) 
                       KDB22GAI = KDB22AIG + ID2AIG(ISYPQ1,ISYMG)
     *                          + NT1AM(ISYPQ1)*(G - 1)
                       KDB22GAB = KDB22ABG + ID2ABG(ISYPQ1,ISYMG)
     *                          + NMATAB(ISYPQ1)*(G - 1)
                       KDB22GIA = KDB22IAG + ID2AIG(ISYPQ1,ISYMG)
     *                          + NT1AM(ISYPQ1)*(G - 1)
                     END IF
                  ELSE IF (MP2) THEN
                     KD2GIJ = KD2IJG + IF2IJG(ISYMPQ,ISYMG)
     *                      + NFROIJ(ISYMPQ)*(G - 1)
                     KD2GIA = KD2IAG + ID2AIG(ISYMPQ,ISYMG)
     *                      + NT1AM(ISYMPQ)*(G - 1)
                     IF (I1DXORD.EQ.1) THEN
                       KDB2GIJ  = KDB2IJG  + IF2IJG(ISYPQ1,ISYMG)
     *                          + NFROIJ(ISYPQ1)*(G - 1)
                       KDB2GIA  = KDB2IAG  + ID2AIG(ISYPQ1,ISYMG)
     *                          + NT1AM(ISYPQ1)*(G - 1)
                       KDB22GIJ = KDB22IJG + IF2IJG(ISYPQ1,ISYMG)
     *                          + NFROIJ(ISYPQ1)*(G - 1)
                       KDB22GIA = KDB22IAG + ID2AIG(ISYPQ1,ISYMG)
     *                          + NT1AM(ISYPQ1)*(G - 1)
                     END IF
                  ELSE IF (CCS) THEN
                     KD2GIJ = KD2IJG + IF2IJG(ISYMPQ,ISYMG)
     *                      + NFROIJ(ISYMPQ)*(G - 1)
                     IF (I1DXORD.EQ.1) THEN
                       KDB2GIJ  = KDB2IJG  + IF2IJG(ISYPQ1,ISYMG)
     *                          + NFROIJ(ISYPQ1)*(G - 1)
                       KDB22GIJ = KDB22IJG + IF2IJG(ISYPQ1,ISYMG)
     *                          + NFROIJ(ISYPQ1)*(G - 1)
                     END IF
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
                        CALL QUIT(
     *                    'Insufficient work space in CCEFFFOCK')
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
                     IF (I1DXORD.EQ.1) THEN
                        CALL QUIT(
     *                    'frozen core not finished for I1DXORD=1.')
                     END IF
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
                  IF (.NOT. LONLYEXP) THEN
C
C                   ----------------------------------------------
C                   compute the contributions to the expectation
C                   values and effective Fock matrices:
C                   ----------------------------------------------
C
                    D = IDEL - IBAS(ISYMD)
C
                    CALL DSCAL(N2BST(ISYM0),HALF,WORK(KDHFAO),1)

                    CALL CCEFFFOCK2(I1DXORD,LABEL,ICOOR,ICORSY,
     &                              MXCOMP,ISYRLX,
     &                              G,ISYMG,D,ISYMD,
     &                              LUFCK,FILFCK,ISYFCK,IDLST,
     &                              LFOCK,LEXPEC,ISYDEN,ISYMOPR,
     &                              WORK(KLAMDP),WORK(KLAMDH),
     &                              WORK(KQAOS),WORK(KDHFAO),
     &                              WORK(KDNS1D),WORK(KKABAO),
     &                              WORK(KB1KABAO),WORK(KB1DHFAO),
     &                              WORK(KB2KABAO),WORK(KB2DHFAO),
     &                              WORK(KD2GAI),  WORK(KD2GAB),  
     &                              WORK(KD2GIJ),  WORK(KD2GIA),
     &                              WORK(KDB2GAI),WORK(KDB2GAB),
     &                              WORK(KDB2GIJ),WORK(KDB2GIA),
     &                              WORK(KDB22GAI),WORK(KDB22GAB),
     &                              WORK(KDB22GIJ),WORK(KDB22GIA),
     &                              DIRINT,LDERINT,
     &                              WORK(KEND4),LWRK4)

                    CALL DSCAL(N2BST(ISYM0),TWO,WORK(KDHFAO),1)

                  ELSE
C
C                   ----------------------------------------
C                   Backtransform density fully to AO basis.
C                   ----------------------------------------
C
                    AUTIM1 = SECOND()
                    IF (CCSD .OR. CC2) THEN
                       CALL CC_DENAO(WORK(KAODEN),ISYMPQ,
     *                               WORK(KD2GAI),WORK(KD2GAB),
     *                               WORK(KD2GIJ),WORK(KD2GIA),ISYMPQ,
     *                               WORK(KLAMDP),1,WORK(KLAMDH),1,
     *                               WORK(KEND4),LWRK4)
                    ELSE
                       CALL CCMP_DAO(WORK(KAODEN),WORK(KD2GIJ),
     *                               WORK(KD2GIA),WORK(KCMO),
     *                               WORK(KLAMDH),WORK(KEND4),
     *                               LWRK4,ISYMPQ)
                    ENDIF
C
C                   ------------------------------------
C                   Add relaxation terms to set up 
C                   effective density. We thus have the
C                   entire effective 2-electron density.
C                   ------------------------------------
C
                    IF (.NOT. CCS) THEN
                       ICON = 2
                       CALL CC_D2EFF(WORK(KAODEN),G,ISYMG,IDEL,ISYMD,
     *                               WORK(KKABAO),WORK(KDHFAO),ICON)
                       CALL CC_D2EFF(WORK(KAODEN),G,ISYMG,IDEL,ISYMD,
     *                               WORK(KDHFAO),WORK(KKABAO),ICON)
                    ENDIF
                    AUTIM1 = SECOND() - AUTIM1
                    TIMDAO = TIMDAO + AUTIM1
C
C                   ----------------------------------
C                   Write effective density to disc for 
C                   subsequent use in integral program,
C                   which performs the contraction of
C                   the density with the 2 e- integrals.
C                   ----------------------------------
C
                    AUTIME = SECOND()
                    NDAD   = NBAST*(IDEL2 - 1) + IGAM
                    NDENEL = N2BST(ISYMPQ)
                    CALL DUMP2DEN(LUDE,WORK(KAODEN),NDENEL,NDAD)
                    AUTIME = SECOND() - AUTIME
                    TIRDAO = TIRDAO + AUTIME
                  END IF
C
  130          CONTINUE
  120       CONTINUE
            END IF
            END DO ! IDLST0 
  110    CONTINUE
C
C------------------------------------------------
C        Loop over number of delta distributions.
C------------------------------------------------
C
         IF (LONLYEXP) THEN
          DO 140 IDEL2 = 1,NUMDIS
C
            IF (NOAUXB) CALL QUIT('Not adapted for R12...')
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
               CALL QUIT(
     *             'Insufficient core for allocation in CCEFFFOCK')
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
                     CALL QUIT(
     *                  'Insufficient work space in CCEFFFOCK')
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
  140     CONTINUE
C
         END IF
C
C---------------------------------------------------------
C        Close file with effective two electron densities.
C---------------------------------------------------------
C
         CLOSE(LUDE,STATUS='DELETE')
C
  100 CONTINUE
C
C----------------------------------------------
C     close loops over fock matrices and atoms:
C----------------------------------------------
C
      END DO ! NTOSYM
      END IF !(LONLYEXP.OR.I1DXORD.EQ.0.OR.(JATOM.EQ.IATOM .AND. LTWO))
      END DO ! IDLST1

  90  CONTINUE
      END DO ! IATOM
C
C---------------------------------------
C     close file for (.NOT. LONLYEXP)
C---------------------------------------
C
      IF (.NOT.LONLYEXP) THEN
        CALL WCLOSE2(LUFCK,FILFCK,'KEEP')
      END IF
C
C------------------------------------------------
C     Restore logical flags for integral program.
C------------------------------------------------
C
      IF (LONLYEXP) THEN
        DIRECT = SAVDIR
        HERDIR = SAVHER
        IF (DAR2EL) DO2DAR = .FALSE.
        IF (IOPREL .EQ. 2) THEN
           DPTINT = .FALSE.
           BPH2OO = .FALSE.
        ENDIF
      END IF
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
C
C----------------------------------------------------------------
C    some debug print out:
C----------------------------------------------------------------
C
      IF ((.NOT.LONLYEXP).AND. (LOCDBG.OR.DEBUG.OR.IPRINT.GE.2)) THEN
         WRITE(LUPRI,'(//,2X,A))') 'Results of CCEFFFOCK:'
         WRITE(LUPRI,'(2X,A,/)')   '======================'
         WRITE(LUPRI,'(2X,A,I5)')  'one-index transformation order:',
     &            I1DXORD 

         WRITE(LUPRI,'(/2X,3A)') 'Label     Sym. ',
     &    ' <1el> CC    <2el> CC    <1el> HF    <2el> HF',
     &    '    Fock idx '
         WRITE(LUPRI,'(2X,72("-"))')

         IF (I1DXORD.EQ.0) THEN
            NFOCKLBL = NEXPFCKLBL
         ELSE IF (I1DXORD.EQ.1) THEN
            NFOCKLBL = N1DXFLBL
         ELSE
            CALL QUIT('Illegal value of I1DXORD in CC_GRAD2.')
         END IF

         DO IDLST = 1, NFOCKLBL
            IF (I1DXORD.EQ.0) THEN
              LABEL   = LBLEXPFCK(IDLST)
              ISYMOPR = ISYEXPFCK(IDLST)
              LEXPEC  = LEXPFCK(1,IDLST)
              LFOCK   = LEXPFCK(2,IDLST)
              ISYFCK  = ISYMOPR
            ELSE
              LABEL  = LBL1DXFCK(IDLST)
              IOPER  = IROPER(LABEL,ISYMOPR)
              LEXPEC = .FALSE.
              LFOCK  = .TRUE.
              LSTRLX = LST1DXFCK(IDLST)
              IDXRLX = IRELAX1DX(IDLST)
              ISYRLX = ILSTSYM(LSTRLX,IDXRLX)
              ISYFCK = MULD2H(ISYMOPR,ISYRLX)
            END IF
            IF (LEXPEC.AND.LFOCK) THEN
             WRITE(LUPRI,'(2X,A8,2I3,4G12.4,2I5)') LABEL, ISYMOPR, 
     &           ISYFCK,
     &           (EXPVALUE(I,IDLST),I=1,4),
     &           (IADRFCK(I,IDLST),I=1,2)
            ELSE IF (LEXPEC) THEN
             WRITE(LUPRI,'(2X,A8,2I3,4G12.4,A10)') LABEL,ISYMOPR,
     &           ISYFCK,
     &           (EXPVALUE(I,IDLST),I=1,4), 
     &           '  ---  ---'
            ELSE IF (LFOCK)  THEN
             WRITE(LUPRI,'(2X,A8,2I3,A48,2I5)') LABEL,ISYMOPR,
     &           ISYFCK,
     &           '  --.-----    --.-----    --.-----    --.-----  ',
     &             (IADRFCK(I,IDLST),I=1,2)
            ELSE
             WRITE(LUPRI,'(2X,A8,2I3,A48,A10)') LABEL,ISYMOPR,
     &           ISYFCK,
     &           '  --.-----    --.-----    --.-----    --.-----  ',
     &           '  ---  ---'
            END IF
         END DO
C
         WRITE(LUPRI,'(2X,72("-"),/)')
      END IF
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
         WRITE(LUPRI,*) 'Total time used in CCEFFFOCK:', TIMTOT
      ENDIF
      IF (IPRINT .GT. 9) THEN
         WRITE(LUPRI,*) 'Time used for setting up d(pq,ga,de):', TIMDEN
         WRITE(LUPRI,*) 'Time used for full AO backtransformation:', 
     *                   TIMDAO
         WRITE(LUPRI,*) 'Time used for reading and writing d and I:',
     *                   TIRDAO
         WRITE(LUPRI,*) 'Time used for calculating 2 e- AO-integrals:',
     *              TIMHE2
         WRITE(LUPRI,*) 'Time used for 1 e- density & intermediates:',
     *              TIMONE
      ENDIF
C
  999 CONTINUE
      CALL QEXIT('CCEFFFOCK')
      RETURN
  165 CALL QUIT('Error reading CCTWODEN')
      END
#endif /* HJAAJ_REMOVED */
C
