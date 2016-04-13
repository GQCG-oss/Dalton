!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2016 (2015), see http://daltonprogram.org"
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
#ifdef OLD_REVISIONLOG
!===========================================================================
!/* Comdeck revision_log */
!
!Revision 1.7  2000/05/24 19:01:23  hjj
!polish; fixing FNDLB2 parameter lists (ftnchek)
!
!Revision 1.6  2000/05/16 09:47:13  hjj
!s/SETSOPPA2/SET2SOPPA/ (so different from SETSOPPA in first six char.)
!Moved ESRCAL outside linear response loop.
!New GETREF calls with LCREF to determine if CSFs or dets. Read NCDETS from LUSIFC.
!Use TRPFLG for SETCI2 (if any operator is triplet, then dets must be used in response).
!Save disk space and i/o: do not write CI diagonal of S[2] matrix (it is KZCONF 1's).
!Polish of some print.
!
!Revision 1.5  2000/05/01 09:56:14  hjj
!Removed reading of KLWOP (KLWOP is removed from common block).
!
!Revision 1.2 Jan 2000     hjj
!polish: print if TRPLET operator (singlet or triplet was
!  missing in output, which made it difficult to follow at times!)
!  removed unnecessary "STOP"s (which also could cause MPI problems)
!  changed (2) to (*) in DIMENSIONs
!  polishing
!
!960304-spas: initialization of C6IFC moved from C6INP to RSPINF
!960304-spas: SOPPA(CCSD) included in RSPINP, RSPMC, RSPINF
!940708-hjaaj: SOPPA changes (KOFFTY = 0, MAXPHP = 0,
!940602-hjaaj: RSPDRV; set UDV and PVX for TDHF open shell
!940408-hjaaj: RSPMC; revised TDHF def. for open shell; removed EACTIV test
!940405-ekd: Added INFSOP if SOPPA-calculation
!931130-hjaaj: RSPMC: revised for proper treatment of ABACUS calls
!  modified a lot of print, incl. inserting '@' in warning/error output
!931116-mjp+hjaaj: RSPMC: set SOPPA
!931011-hjaaj: C6INP: added C6IFC definition
!931007-hjaaj: added RSPMC and GETREF (moved from rspe2c.u)
!920721-Hinne Hettema
!C6INP : inserted missing GO TO 100; corrected output of FREQUENCIES
!RSPSET: define new RSPSUP variable (use super symmetry)
!RSPDRV: initialize arrays for symmetrizing MCSCF gradient.
!===========================================================================
#endif
C
C  /* Deck rspdrv */
      SUBROUTINE RSPDRV(WRK,LWRK)
C
C Driver for RESPONSE.
C
#include "implicit.h"
C
      DIMENSION WRK(LWRK)
C
#include "iratdef.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
C
C Used from common blocks:
C INFTRA : USEDRC
C
#include "maxorb.h"
#include "priunit.h"
#include "pgroup.h"
#include "inforb.h"
#include "infpri.h"
#include "infmp2.h"
#include "wrkrsp.h"
#include "infrsp.h"
#include "infdim.h"
#include "inftra.h"
#include "inftap.h"
#include "rspprp.h"
#include "inflr.h"
#include "infpp.h"
#include "infc6.h"
#include "infhyp.h"
#include "infsmo.h"
#include "infs0.h"
#include "infesr.h"
#include "infcr.h"
#include "inftmo.h"
#include "inftpa.h"
#include "infsop.h"
#include "gnrinf.h"
#include "infinp.h"
#include "absorp.h"
#include "esg.h"
#include "numder.h"
#include "gtensor.h"
#include "infave.h"
#include "mxcent.h"
#include "esrhfc.h"
#include "abslrs.h"
#include "infvar.h"
C
      LOGICAL SAVLR, DORESP, FLAGSV, EX, NEED_PVX
C
      CALL QENTER('RSPDRV')
C
      TIMRSP = SECOND()
      CALL GETTIM(CSTR,WSTR)
      WRITE (LUPRI,'(/1X,78A1/A/1X,78A1)') ('-',I=1,78),
     &   '  RESPONSE  -  '//
     &   'an MCSCF, MC-srDFT, DFT, and SOPPA response property program',
     &   ('-',I=1,78)
C
C     *************************
C     ***** INPUT SECTION *****
C     *************************
C
C Define wave function variables, also if SIRIUS has not been called
C
      CALL SETSIR(WRK,LWRK)
C
C Define variables that are common for all perturbation symmetries
C
      CALL RSPSET
C
C Read **RESPONS input
C
      CALL RSPREA
C
C  If SOPPA then set up information from MP2 and find number of 2p-2h
C  excitations and offsets to various index arrays.
C
      IF (SOPPA) THEN
         CALL MP2SET
         CALL PHADR1(WRK,1) ! to calculate NPHTOT(1) needed in RSPACT
         CALL SETSOPPA
      END IF
C
C     set *ACT* variables in INFRSP (different settings for SOPPA and
C     MCSCF/DFT)
C
      CALL RSPACT
C
      IF (INPTES) THEN
         WRITE(LUPRI,'(//A/)') ' *** END of input test for RESPONSE ***'
         CALL QUIT('*** END of input test for RESPONSE ***')
      END IF
      CALL FLSHFO(LUPRI)
C
C     ORGANIZE CALCULATION FOR EACH PERTURBATION OPERATOR
C
      IF (DFT_SO) THEN
         IF (RSPCI) THEN
            CALL QUIT('.DFT_SO is not compatible with CI response')
         END IF
         KCMO  = 1
         KLST  = KCMO   + NCMOT
         LFREE = LWRK-KLST
         IF (KLST.GT.LWRK) CALL ERRWRK('RSPDRV 0',KLST-1,LWRK)
         REWIND LUSIFC
         CALL MOLLAB(LBSIFC,LUSIFC,LUPRI)
         READ (LUSIFC)
         READ (LUSIFC)
         CALL READT (LUSIFC,NCMOT,WRK(KCMO))
         CALL DFTSOI(WRK(KCMO),WRK(KLST),LFREE,IPRONE)
      END IF
      KWSYM = 1
      LWSYM = LWRK  - KWSYM
      IF (LWSYM.LT.0) CALL ERRWRK('RSPDRV 2',KWSYM,LWRK)
      CALL RSPSYM(WRK(KWSYM),LWSYM)
C
C     990505-hjaaj:
C     Check if we only want expectation values
C     (we don't need PV or MO integrals in that case)
C
      DORESP = CRCAL.OR.TOMOM.OR.TPAMP
      IF (HYPCAL.OR.SOMOM.OR.EXMOM) THEN
         IF (DORESP) THEN
            NWARN = NWARN + 1
            WRITE (LUPRI,'(/A/A)')
     &      ' WARNING: Quadratic response input ignored',
     &      ' because Cubic response has been requested.'
         END IF
         DORESP = .TRUE.
      END IF
      I = 0
      DO KSYMOP = 1,NSYM
         IF ( NPPCNV(KSYMOP).GT.0 .OR. NGPLR(KSYMOP).GT.0
     *       .OR.  NGPC6(KSYMOP).GT.0 ) I = I + 1
      END DO
      IF (I.GT.0) THEN
         IF (DORESP) THEN
            NWARN = NWARN + 1
            WRITE (LUPRI,'(/A/A)')
     &      ' WARNING: all PP,LR,C6 input ignored because',
     &      ' because cubic or quadratic response has been requested.'
         END IF
      END IF
      DORESP = DORESP .OR. ESRCAL .OR. ABSORP .OR. ABSLRS .OR. GCALC
     &                .OR. I.GT.0
      DORESP = DORESP .OR. ESRHFC
      TRPFLG = TRPFLG .OR. ESRHFC
C
C     PERFORM INTEGRAL TRANSFORMATION
C
      IF (.NOT.NOITRA .AND. DORESP) THEN
         KCMO   = 1
         KWTRA  = KCMO   + NCMOT
         LWTRA  = LWRK   - KWTRA
         IF (LWTRA.LT.0) CALL ERRWRK('RSPDRV 1',KWTRA-1,LWRK)
         REWIND LUSIFC
         IF (RSPCI) THEN
            CALL MOLLAB('CIRESPON',LUSIFC,LUPRI)
         ELSE
            CALL MOLLAB(LBSIFC,LUSIFC,LUPRI)
         END IF
         READ (LUSIFC)
         READ (LUSIFC)
         CALL READT (LUSIFC,NCMOT,WRK(KCMO))
         IF ( .NOT. RSPCI ) USEDRC = .TRUE.
         IF (SOPPA.OR.HYPCAL.OR.SOMOM.OR.EXMOM.OR.
     &       CRCAL.OR.TOMOM.OR.TPAMP) THEN
            JTRLVL = -10
C  TODO hjaaj:  When NEWTRA and USEDRC, all integrals are
C           stored and calculated twice, both as Mulliken
C           and Dirac; but only occ-occ Dirac distributions
C           are needed in rspoli.F (if .not. NEWTRA then
C           only occ-occ Dirac dist. are saved on MODRCINT).
C           hjaaj Nov 2001
         ELSE
            JTRLVL = -4
         END IF
         FLAGSV = DORSP
         DORSP = .TRUE.
         CALL SIR_INTOPEN
         DORSP = FLAGSV
         CALL TRACTL(JTRLVL,WRK(KCMO),WRK(KWTRA),LWTRA)
         CALL FLSHFO(LUPRI)
      ELSE IF (SOPPA) THEN
         FLAGSV = DORSP
         DORSP = .TRUE.
         CALL SIR_INTOPEN
         DORSP = FLAGSV
      END IF
C
C     ALLOCATE WORK SPACE FOR MATRICES THAT WILL BE KEPT DURING THE
C     WHOLE RESPONSE CALCULATION AND READ IN THE MATRICES
C
      NEED_PVX = DORESP .OR. N2AVTOT .GT. 0
      KFREE  = 1
      LFREE  = LWRK
      CALL MEMGET('REAL',KINDX,LCINDX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KCMO ,NCMOT ,WRK,KFREE,LFREE)
      LUDV   = NACTT*NACTT*2
      CALL MEMGET('REAL',KUDV ,LUDV  ,WRK,KFREE,LFREE)
      IF (RSPCI) THEN
         CALL MEMGET('REAL',KPVX ,0     ,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KFOCK,0     ,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KFC  ,0     ,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KFV  ,0     ,WRK,KFREE,LFREE)
      ELSE
         IF ( .NOT. NEED_PVX ) THEN
!        two-elec. dens.mat. not needed for one-el. exp.val.s
            LPVX = 0
         ELSE IF (TRPFLG) THEN
C NEED BOTH TRIPLET AND SINGLET TWO ELECTRON DENSITY MATRICES
            LPVX = 2*LPVMAT
         ELSE
C NEED ONLY SINGLET TWO ELECTRON DENSITY MATRIX
            LPVX = LPVMAT
         END IF
         CALL MEMGET('REAL',KPVX ,LPVX  ,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KFOCK,N2ORBT,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KFC  ,NNORBT,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KFV  ,NNORBT,WRK,KFREE,LFREE)
      END IF
      CALL MEMGET('REAL',KFCAC,NNASHX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2AC,NNASHX*NNASHX,WRK,KFREE,LFREE)
      KTOT  =  KFREE
      KWRK1  = KFREE
      LWRK1  = LFREE
C
      IF (SOPPA) THEN
C
C  Initialize XINDX
C
         CALL DZERO(WRK(KINDX),LCINDX)
C
C  Find address array's for SOPPA calculation
C
         CALL SET2SOPPA(WRK(KINDX+KABSAD-1),WRK(KINDX+KABTAD-1),
     *                  WRK(KINDX+KIJSAD-1),WRK(KINDX+KIJTAD-1),
     *                  WRK(KINDX+KIJ1AD-1),WRK(KINDX+KIJ2AD-1),
     *                  WRK(KINDX+KIJ3AD-1),WRK(KINDX+KIADR1-1))
      ENDIF
C
C
      IF (.NOT. NEED_PVX) KPVX = KWRK1
      CALL RSPMC(WRK(KCMO),WRK(KUDV),WRK(KPVX),WRK(KFOCK),WRK(KFC),
     *           WRK(KFV),WRK(KFCAC),WRK(KH2AC),WRK(KINDX),
     *           WRK(KWRK1),LWRK1)
C     CALL RSPMC(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,WRK,LWRK)
C
C
      IF (DOMCSRDFT) THEN
         CALL TITLER('MC-srDFT response calculation','*',103)
      ELSE IF (DOHFSRDFT) THEN
         CALL TITLER('HF-srDFT response calculation','*',103)
      ELSE IF (DODFT .AND. NASHT.GT.0) THEN
         CALL TITLER('High-spin DFT response calculation','*',103)
      ELSE IF (DODFT) THEN
         CALL TITLER('DFT response calculation (TD-DFT)','*',103)
      ELSE IF (HSROHF) THEN
         CALL TITLER('High-spin ROHF response calculation','*',103)
      ELSE IF (TDHF) THEN
         CALL TITLER('RHF response calculation (TDHF)','*',103)
      ELSE IF (CCPPA) THEN
         CALL TITLER('SOPPA(CCSD) response calculation','*',103)
         IF (PARCAL) CALL PARQUIT('SOPPA(CCSD) RESPONSE')
      ELSE IF (SOPPA) THEN
         CALL TITLER('SOPPA response calculation','*',103)
         IF (PARCAL) CALL PARQUIT('SOPPA RESPONSE')
      ELSE IF (RSPCI) THEN
         CALL TITLER('CI response calculation','*',103)
      ELSE
         CALL TITLER('MCSCF response calculation','*',103)
         ! IF (PARCAL) CALL PARQUIT('MCSCF RESPONSE')
      END IF
C
C
      CALL ASSERT_TRIPLET_RESPONSE_WITH_DETERMINANTS()
C
      IF ( NEED_PVX ) THEN
C
      IF ( .NOT.TDHF ) THEN
         CALL GETCIX(WRK(KINDX),IREFSY,IREFSY,WRK(KFREE),LFREE,0)
      END IF
C
      IF ((.NOT.TDHF) .AND. (.NOT.RSPCI) ) THEN
C
C     ... CALCULATE ONE- AND TWO- BODY DENSITY MATRICES
C         ( THE SYMMETRIC MC TWO BODY DENSITY MATRIX CANNOT BE USED IN
C           RESPONSE CALCULATION )
C
         KCREF  = KWRK1
         KTOT   = KCREF + NCREF
         LFREE  = LWRK  - KTOT
         IF (LFREE.LT.0) CALL ERRWRK('RSPDRV 4',KTOT-1,LWRK)
         KFREE_KTOT  = 1
C
         CALL GETREF(WRK(KCREF),NCREF)
         IF ( IPRRSP.GT.110 ) THEN
            WRITE(LUPRI,'(/A)')
     &           ' ***** ONE BODY DENSITY MATRIX from SIRIFC file'
            CALL OUTPUT(WRK(KUDV),1,NASHT,1,NASHT,NASHT,NASHT,-1,LUPRI)
         END IF
         ISPIN1 = 0
         ISPIN2 = 0
         CALL RSPDM(IREFSY,IREFSY,NCREF,NCREF,WRK(KCREF),WRK(KCREF),
     *              WRK(KUDV),WRK(KPVX),
     *              ISPIN1,ISPIN2,.FALSE.,.FALSE.,
     *              WRK(KINDX),WRK(KTOT),KFREE_KTOT,LFREE)
C        CALL RSPDM(ILSYM,IRSYM,NCLDIM,NCRDIM,CL,CR, RHO1,RHO2,
C    *              ISPIN1,ISPIN2,TDM,NORHO2,XNDXCI,WORK,
C    *              KFREE,LFREE)
C
         IF (TRPFLG) THEN
            ISPIN1 = 1
            ISPIN2 = 1
            CALL RSPDM(IREFSY,IREFSY,NCREF,NCREF,WRK(KCREF),WRK(KCREF),
     *                 WRK(KUDV),WRK(KPVX+LPVMAT),
     *                 ISPIN1,ISPIN2,.FALSE.,.FALSE.,
     *                 WRK(KINDX),WRK(KTOT),KFREE_KTOT,LFREE)
C
         END IF
C
C Include triplet one-electron density
C
            CALL RSPDM(IREFSY,IREFSY,NCREF,NCREF,WRK(KCREF),WRK(KCREF),
     *                 WRK(KUDV+N2ASHX),WRK(KFREE),
     *                 1,0,.FALSE.,.TRUE.,
     *                 WRK(KINDX),WRK(KTOT),KFREE,LFREE)
         IF ( IPRRSP.GT.110 ) THEN
            WRITE(LUPRI,'(/A)')' ** ONE BODY DENSITY MATRIX FROM RSPDM'
            CALL OUTPUT(WRK(KUDV),1,NASHT,1,NASHT,NASHT,NASHT,-1,LUPRI)
            KUDV2=KUDV+N2ASHX
            CALL OUTPUT(WRK(KUDV2),1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
         END IF
      ELSE IF ( TDHF .AND. NASHT .EQ. 1 ) THEN
         WRK(KUDV) = D1
         WRK(KUDV+1) = D1
         WRK(KPVX) = D0
         IF (TRPLET) WRK(KPVX+1) = D0
      ELSE IF ( TDHF .AND. HSROHF ) THEN
         CALL DUNIT(WRK(KUDV),NASHT)
         CALL DZERO(WRK(KPVX),LPVMAT)
      END IF
      END IF ! for IF ( NEED_PVX )
C
C     Calculate expectation value properties (first order)
C     ====================================================
C
      CALL MIXS0(WRK(KCMO),WRK(KUDV),WRK(KTOT),LFREE)
C     CALL MIXS0(CMO,UDV,WRK,LWRK)
C
      CALL AVEPRP(WRK(KCMO),WRK(KUDV),WRK(KPVX),WRK(KTOT),LFREE)
      CALL FLSHFO(LUPRI)
      CALL ZFSDRV(WRK,LWRK)
C
      IF (.NOT.DORESP) GO TO 9999
C
C Open RSPVEC for response vectors, existing file may exist for restart
C
      CALL GPINQ('RSPVEC','EXIST',EX)
      IF (.NOT.EX) THEN
         CALL GPOPEN(LURSP,'RSPVEC','NEW',' ','UNFORMATTED',IDUMMY,
     &        .FALSE.)
         WRITE (LURSP) NISH,NASH,NORB,NBAS,NSYM
         WRITE (LURSP) 'EOFLABEL'
      ELSE
         CALL GPOPEN(LURSP,'RSPVEC','OLD',' ','UNFORMATTED',IDUMMY,
     &        .FALSE.)
      END IF
C Open files for trial and sigma vectors
      CALL GPOPEN(LURSP3,'RSP_BVEC','UNKNOWN',' ',
     &   'UNFORMATTED',IDUMMY,.FALSE.)
      CALL GPOPEN(LURSP4,'RSP_DIAG','UNKNOWN',' ',
     &   'UNFORMATTED',IDUMMY,.FALSE.)
      CALL GPOPEN(LURSP5,'RSP_SVEC','UNKNOWN',' ',
     &   'UNFORMATTED',IDUMMY,.FALSE.)
C
C     Calculate open shell triplet expectation values
C     ====================================================
C
      IF ( ESRCAL ) THEN
         KSYMOP = 1
C
C        DEFINE VARIABLES THAT DEPEND ON SYMMETRY
C
         CALL RSPVAR(WRK(KUDV),WRK(KFOCK),WRK(KFC),WRK(KFV),
     *               WRK(KFCAC),WRK(KH2AC),WRK(KINDX),
     *               WRK(KWRK1),LWRK1)
         CALL RSPESR(WRK(KCMO),WRK(KUDV),WRK(KPVX),
     *               WRK(KFC),WRK(KFV),WRK(KFCAC),WRK(KH2AC),
     *               WRK(KINDX),WRK(KWRK1),LWRK1)
         CALL FLSHFO(LUPRI)
      END IF
C
      IF ( GCALC ) THEN ! calculate g-tensor etc.
         CALL GDRV(WRK,LWRK)
      END IF
C
C     Hyperfine coupling constants driver for RESPONSE
C     ====================================================
C
      IF (ESRHFC) THEN
         CALL HFCDRV(WRK(KCMO),WRK(KUDV),WRK(KPVX),WRK(KFOCK),WRK(KFC),
     *               WRK(KFV),WRK(KFCAC),WRK(KH2AC),WRK(KINDX),
     *               WRK(KWRK1),LWRK1)
      END IF
C
C     Calculate response properties with absorption
C     =============================================
C
      IF (ABSORP .OR. (ABSLRS .AND. (NCONF.GT.1))) THEN
C     Use the old linear response solver
C     for MCSCF (or SCF if OLDCPP is requested by user)
         CALL ABSCALC(WRK(KCMO),WRK(KUDV),WRK(KPVX),WRK(KFOCK),
     *                WRK(KFC),WRK(KFV),WRK(KFCAC),WRK(KH2AC),
     *                WRK(KINDX),WRK(KWRK1),LWRK1)
         GO TO 9999
      END IF
C
C     Calculate response properties with absorption using new cpp
C     ===========================================================
C
      IF (ABSLRS) THEN
         CALL ABSLRSCALC(WRK(KCMO),WRK(KUDV),WRK(KPVX),WRK(KFOCK),
     *                WRK(KFC),WRK(KFV),WRK(KFCAC),WRK(KH2AC),
     *                WRK(KINDX),WRK(KWRK1),LWRK1)
         GO TO 9999
      END IF
C
C     Calculate cubic response properties (fourth order)
C     ==================================================
C
      IF (CRCAL.OR.TOMOM.OR.TPAMP) THEN
         CALL CRCALC(WRK(KCMO),WRK(KUDV),WRK(KPVX),WRK(KFOCK),
     *               WRK(KFC),WRK(KFV),WRK(KFCAC),WRK(KH2AC),
     *               WRK(KINDX),WRK(KWRK1),LWRK1)
         GO TO 9999
      END IF
C
C     Calculate quadratic response properties (third order)
C     =====================================================
C
      IF (HYPCAL.OR.SOMOM.OR.EXMOM) THEN
         CALL QRCALC(WRK(KCMO),WRK(KUDV),WRK(KPVX),WRK(KFOCK),
     *               WRK(KFC),WRK(KFV),WRK(KFCAC),WRK(KH2AC),
     *               WRK(KINDX),WRK(KWRK1),LWRK1)
         GO TO 9999
      END IF
C
C     Calculate linear response properties (second order)
C     ===================================================
C
      DO 100 KSYMOP = 1,NSYM
         IF ( NPPCNV(KSYMOP).EQ.0 .AND. NGPLR(KSYMOP).EQ.0
     *       .AND.  NGPC6(KSYMOP).EQ.0 )
     *      GO TO 100
C     ... skip this symmetry if no operators
         WRITE(LUPRI,'(A/A/A,I5,3A/A)')
     &      ' -------------------------------------------------'//
     &      '---------------',
     &      ' ----- Linear response calculation',
     &      ' ----- Symmetry of excitation/property operator(s)',
     &      KSYMOP,'  ( ',REP(KSYMOP-1),')',
     &      ' -------------------------------------------------'//
     &      '---------------'
         WRITE(LUPRI,'(3(/A,I5))')
     *      ' Number of excitations of this symmetry        ',
     *      NPPCNV(KSYMOP),
     *      ' Number of response properties of this symmetry',
     *      NGPLR(KSYMOP),
     *      ' Number of C6/C8 properties of this symmetry   ',
     *      NGPC6(KSYMOP)
         IF (KSYMOP .EQ. 1 .AND. RSPSUP) THEN
            NINFO = NINFO + 1
            WRITE (LUPRI,'(3(/A))')
     &      ' INFORMATION: Note that only the totally symmetric',
     &      '   representation of the super symmetry is treated',
     &      '   not all of symmetry 1 in the Abelian symmetry'
         END IF
C
C        DEFINE VARIABLES THAT DEPEND ON SYMMETRY
C
         CALL RSPVAR(WRK(KUDV),WRK(KFOCK),WRK(KFC),WRK(KFV),
     *               WRK(KFCAC),WRK(KH2AC),WRK(KINDX),
     *               WRK(KWRK1),LWRK1)
C        CALL RSPVAR(UDV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
         IF ( KZVAR.EQ.0) THEN
            NWARN = NWARN + 1
            WRITE(LUPRI,'(/A,I3)')'@ ****WARNING******'//
     *         ' No variables in symmetry',KSYMOP
            GO TO 100
         END IF
         CALL FLSHFO(LUPRI)
         IF ( NPPCNV(KSYMOP) .GT.0 ) THEN
C
C        ... FIND EXCITATION ENERGIES AND TRANSITION MOMENTS
C
            IF (IPRRSP .GT. 0) TIMPP = SECOND()
            CALL RSPPP(WRK(KCMO),WRK(KUDV),WRK(KPVX),
     *                 WRK(KFC),WRK(KFV),WRK(KFCAC),WRK(KH2AC),
     *                 WRK(KINDX),WRK(KWRK1),LWRK1)
C           CALL RSPPP(CMO,UDV,PV,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
            IF (IPRRSP .GT. 0) THEN
               TIMPP = SECOND() - TIMPP
               WRITE (LUPRI,'(//A,F10.2,A,I2)')
     *         ' Time used in polarization propagator calculation is',
     *         TIMPP,' CPU seconds for symmetry',KSYMOP
            END IF
         END IF
         IF ( NGPLR(KSYMOP) .GT. 0 ) THEN
C
C        ... DETERMINE SECOND ORDER MOLECULAR PROPERTIES
C
            IF (IPRRSP .GT. 0) TIMLR = SECOND()
            CALL RSPLR(WRK(KCMO),WRK(KUDV),WRK(KPVX),
     *                 WRK(KFC),WRK(KFV),WRK(KFCAC),WRK(KH2AC),
     *                 WRK(KINDX),WRK(KWRK1),LWRK1)
            IF (IPRRSP .GT. 0) THEN
               TIMLR = SECOND() - TIMLR
               WRITE (LUPRI,'(//A,F10.2,A,I2)')
     *         ' Time used in linear response calculation is',
     *         TIMLR,' CPU seconds for symmetry',KSYMOP
            END IF
         END IF
         CALL FLSHFO(LUPRI)
         IF ( NGPC6(KSYMOP) .GT. 0 ) THEN
C
C        ... DETERMINE ALPHA(IW)
C
            IF (IPRRSP .GT. 0) TIMC6 = SECOND()
            SAVLR  = RESTLR
            RESTLR = RESTC6
            CALL RSPC6(WRK(KCMO),WRK(KUDV),WRK(KPVX),
     *                 WRK(KFC),WRK(KFV),WRK(KFCAC),WRK(KH2AC),
     *                 WRK(KINDX),WRK(KWRK1),LWRK1)
            RESTLR = SAVLR
            RESTC6 = .FALSE.
            IF (IPRRSP .GT. 0) THEN
               TIMC6 = SECOND() - TIMC6
               WRITE (LUPRI,'(//A,F10.2,A,I2)')
     *         ' Time used in c6/c8 calculation ',
     *         TIMC6,' CPU seconds for symmetry',KSYMOP
            END IF
         END IF
         CALL FLSHFO(LUPRI)
C
C  CALCULATE THE EXCITED STATE GRADIENT
C
         IF ( ESG .AND. ISYME.EQ.KSYMOP .AND.
     $        .NOT. (NMORDR.GT.0 .AND. NARDRP.EQ.0)) THEN
          ESGTIM_XVECS = TIMPP
          CALL RSPESG(WRK(KCMO),WRK(KUDV),WRK(KPVX),WRK(KFOCK),
     *                WRK(KFC),WRK(KFV),WRK(KFCAC),WRK(KH2AC),
     *                WRK(KINDX),WRK(KWRK1),LWRK1)
         END IF
  100 CONTINUE
      CALL GPCLOSE(LUSIFC,'KEEP')
      CALL GPCLOSE(LURSP ,'KEEP')
      CALL GPCLOSE(LURSP3,'DELETE')
      CALL GPCLOSE(LURSP4,'DELETE')
      CALL GPCLOSE(LURSP5,'DELETE')
      IF (LUINTM .GT. 0) CALL GPCLOSE(LUINTM,'KEEP')
      IF (LUMHSO .GT. 0) CALL DARMOV(LUMHSO)
C
 9999 CALL GETTIM(CEND,WEND)
      CTOT = CEND - CSTR
      WTOT = WEND - WSTR
      WRITE (LUPRI,'()')
      IF (QMMM) CALL QMMMTIMES('RESPONSE')
      CALL TIMTXT(' Total CPU  time used in RESPONSE:',CTOT,LUPRI)
      CALL TIMTXT(' Total wall time used in RESPONSE:',WTOT,LUPRI)
C
C     *******************************************
C
      CALL MEMREL('End of RSPDRV',WRK,1,1,KFREE,LFREE)

      CALL QEXIT('RSPDRV')
      RETURN
      END
C  /* Deck rsprea */
      SUBROUTINE RSPREA
C
C GENERAL INPUT FOR RESPONSE
C
#include "implicit.h"
#include "dummy.h"
      PARAMETER ( NTABLE = 11 )
C
#include "priunit.h"
#include "infpri.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "rspprp.h"
#include "inflr.h"
#include "infpp.h"
#include "infc6.h"
#include "infhyp.h"
#include "infsmo.h"
#include "infesr.h"
#include "infcr.h"
#include "inftmo.h"
#include "inftpa.h"
#include "maxorb.h"
#include "infinp.h"
#include "absorp.h"
#include "mxcent.h"
#include "esrhfc.h"
#include "abslrs.h"
! JMHO: pcmlog.h included to check if CUBIC and PCM -> quit
#include "pcmlog.h"
C
      CHARACTER WORD*7, PROMPT*1, TABLE(NTABLE)*7
C
      DATA TABLE /'*END OF', '**RESPO', '*LINEAR', '*QUADRA',
     *            '*CUBIC ', '*C6    ', '*ESR   ', '*ABSORP',
     *            '*SPIN-O', '*ESG   ', '*HFC   '/

      logical, external :: fun_is_ready_for_qr
      logical, external :: fun_is_ready_for_cr
C
C     Include block data RSPINF which initializes /INFRSP/
C
      EXTERNAL RSPINF
C
      WRITE (LUPRI,'(//,A,//)')
     *   ' -------- OUTPUT FROM RESPONSE INPUT PROCESSING --------'
C
C     ***** PROCESS INPUT FOR RESPONSE *****
C
      CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
C
C Variable settings
C
      DO I=1,8
         NPPSTV(I) = 0
         NPPSIM(I) = 0
         NPPCNV(I) = 0
         NTPCN1(I) = 0
         NTPCN2(I) = 0
         NTMCNV(I) = 0
         DO J=1,MAXLPP
            LCMEXC(J,I) = .TRUE.
         END DO
      END DO
C
C Non-equilibrium solvent parameter read in SIRIUS needs to be reset
C
      IF (INERSI) EPSOL = EPPN
C
      REWIND (LUCMD,IOSTAT=IOS)
  100 READ (LUCMD, '(A7)',END=2100) WORD
      CALL UPCASE(WORD)
      IF (WORD .NE. TABLE(2)) GO TO 100
C
  200 CONTINUE
      IF (WORD(1:2) .EQ. '**' .AND. .NOT. WORD .EQ. TABLE(2)) GO TO 1
C     end of this input module
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '*') THEN
         DO I=1, NTABLE
            IF (WORD .EQ. TABLE(I)) THEN
               GO TO (1,2,3,4,5,6,7,8,9,10,11), I
            END IF
         END DO
C
         WRITE (LUPRI,'(/,3A,/)')
     *        ' Keyword ',WORD,' not recognized in RSPREA.'
         CALL PRTAB(NTABLE,TABLE,WORD//' input keywords',LUPRI)

         CALL QUIT(' ILLEGAL "*" OPTION IN RSPREA')
      ELSE
         WRITE (LUPRI,'(/,3A,/)')
     *      '@PROMPT "',PROMPT,'" ILLEGAL OR OUT OF ORDER.'
         CALL QUIT(' PROGRAM STOPPED IN RSPREA, ERROR IN PROMPT')
      END IF
 2    WORD = '*RESPON'
      CALL RSPINP(WORD)
      GO TO 200
 3    CALL LINEAR(WORD)
      GO TO 200

 4    call quadra(word)
#ifndef DISABLE_XC_RESPONSE_SANITY_CHECK
      if (dodft) then
         if (.not. fun_is_ready_for_qr()) then
            write(lupri, *) 'ERROR: functional not fully implemented'
            write(lupri, *) '       or tested for QR'
            write(lupri, *) 'to disable this stop recompile'
            write(lupri, *) 'with -DDISABLE_XC_RESPONSE_SANITY_CHECK'
            write(lupri, *) 'note that GGAKEY functionals are always'
            write(lupri, *) 'stopped although they could be correct'
            call quit('functional not fully implemented/tested for QR')
         end if
      end if
#endif
      go to 200

 5    call cubic(word)
! JMHO: quit if CUBIC and PCM (until bug is found)
      if (pcm) then
        call quit('PCM and cubic response not available')
      end if
#ifndef DISABLE_XC_RESPONSE_SANITY_CHECK
      if (dodft) then
         if (.not. fun_is_ready_for_cr()) then
            write(lupri, *) 'ERROR: functional not fully implemented'
            write(lupri, *) '       or tested for CR'
            write(lupri, *) 'to disable this stop recompile'
            write(lupri, *) 'with -DDISABLE_XC_RESPONSE_SANITY_CHECK'
            write(lupri, *) 'note that GGAKEY functionals are always'
            write(lupri, *) 'stopped although they could be correct'
            call quit('functional not fully implemented/tested for CR')
         end if
      end if
#endif
      go to 200

 6    CALL C6INP(WORD)
      GO TO 200
 7    CALL ESRINP(WORD)
      GO TO 200
 8    CALL ABSORP_INPUT(WORD)
      GO TO 200
 9    CALL HSOINP(WORD)
      GO TO 200
 10   CALL ESGINP(WORD)
      GO TO 200
 11   CALL HFCINP(WORD)
      GO TO 200
C
 1    CONTINUE
      IPRRSP_default = IPRRSP
      IPRRSP = MAX(IPRPP ,IPRLR ,IPRC6)         ! LR
      IPRRSP = MAX(IPRRSP,IPRHYP,IPRSMO)        ! QR
      IPRRSP = MAX(IPRRSP,IPRCR,IPRTMO,IPRTPA)  ! CR
      IPRRSP = MAX(IPRRSP,IPRESR)
      IPRRSP = MAX(IPRRSP,IPRABS)
      IF (IPRRSP .ne. IPRRSP_default) THEN
         WRITE(LUPRI,'(/A,I10)')
     &   ' General print level in RESPONS is IPRRSP =',IPRRSP
      END IF
C
C     We arrive at 2100 if no input file found. We return.
C
 2100 CONTINUE
      CALL GPCLOSE(LUCMD,'KEEP')
C
      RETURN
C
C *** END OR RSPREA
C
      END
C  /* Deck rspinp */
      SUBROUTINE RSPINP(WORD)
C
C     Read general ".xxx" input options for "**RESPONSE"
C
#include "implicit.h"
C
#include "priunit.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "infs0.h"
#include "infpri.h"
#include "infdim.h"
#include "maxorb.h"
#include "infvar.h"
#include "gnrinf.h"
#include "infinp.h"
#include "rspprp.h"
#include "qrinf.h"
#include "infave.h"
#include "inftap.h"
#include "inforb.h"
#include "abslrs.h"
C
      LOGICAL NEWDEF, SET
      PARAMETER ( NTABLE = 25 )
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
      SAVE TABLE, SET
C
      DATA TABLE /'.INPTES', '.NOAVDI', '.SOPPA ',
     *            'XXXXXXX', '.NOITRA', '.ORBSFT',
     *            '.ORBSPC', '.TRPFLG', '.S0MIX ',
     *            '.OPTORB', '.MAXPHP', '.PHPRES',
     *            '.MAXRM ', '.QRREST', '.PROPAV',
     *            '.HIRPA ', '.SOPW4 ', '.SOPPA(',
     *            '.NODOIT', '.DFT-SO', '.PROP2A',
     *            '.TDA   ', '.CIS   ', '.SOPRPA',
     *            '.TRDQF '/
      DATA SET /.FALSE./
C
      IF (.NOT.SET) THEN
         SET = .TRUE.
         IF (.NOT. DIRFCK) DIRFCK = DIRCAL
         IF (DOHFSRDFT .OR. DOMCSRDFT) DIRFCK = .TRUE.
         IF (NASHT .LE. 1 .OR. HSROHF) NOITRA = DIRFCK
         IF (DOMCSRDFT .OR. DOMC) NOITRA = .FALSE.
         IF (DOMC) DIRFCK = .FALSE.
         MAXPHP = 100
         N1AVTOT = 0
         N2AVTOT = 0
         TRDQF = .FALSE.
      END IF
      NEWDEF = (WORD .EQ. '*RESPON')
      ICHANG = 0
      IF (NEWDEF) THEN
      WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') GO TO 100
            IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,
     *                      11,12,13,14,15,16,17,18,19,20,21,22,23,
     *                      24,25), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') '@KEYWORD "',WORD,
     *            '" NOT RECOGNIZED IN RSPINP.'
               CALL QUIT(' ILLEGAL KEYWORD IN RSPINP ')
    1          CONTINUE ! .INPTEST
                  INPTES = .TRUE.
               GO TO 100
    2          CONTINUE ! .NOAVDIA
                  AVDIA  = .FALSE.
               GO TO 100
    3          CONTINUE ! .SOPPA
                  SOPPA = .TRUE.
                  IF (PARCAL) THEN
                    CALL PARQUIT("SOPPA combined with linear response")
                  END IF
                  IF (SOPPA .AND. EMBEDDING) THEN
cjje: change this if full solvent soppa is ever implemented
                    GO TO 24 ! go to .SOPRPA option
                  END IF
               GO TO 100
    4          CONTINUE
               GO TO 100
    5          CONTINUE ! .NOITRA
                  NOITRA = .TRUE.
               GO TO 100
    6          CONTINUE ! .ORBSFT
                  READ (LUCMD,*) ORBSFT
               GO TO 100
    7          CONTINUE ! .ORBSPC
                  ORBSPC = .TRUE.
               GO TO 100
    8          CONTINUE ! TRPFLG
                  TRPLET = .TRUE.
                  TRPFLG = .TRUE.
               GO TO 100
    9          CONTINUE ! .S0MIX
                  NOS0MX = .FALSE.
               GO TO 100
   10          CONTINUE ! .OPTORB
                  OPTORB = .TRUE.
               GO TO 100
   11          CONTINUE ! .MAXPHP
                  READ (LUCMD,*) MAXPHP
               GO TO 100
   12          CONTINUE ! .PHPRES
                  PHPRES = .TRUE.
               GO TO 100
   13          CONTINUE ! .MAXRM
                  READ (LUCMD,*) MAXRM
               GO TO 100
   14          CONTINUE
                  QRREST = .TRUE.
               GO TO 100
   15          CONTINUE ! .PROPAV
                  N1AVTOT = N1AVTOT + 1
                  IF ( N1AVTOT .GT. MAX1AV ) THEN
                     WRITE(LUPRI,'(/,2A,I5)')
     *               '@rspinp: NUMBER OF OPERATORS FOR .PROPAV'
     *               ,'CALCULATION EXCEEDS THE ALLOWED NUMBER :'
     *               ,MAX1AV
                   CALL QUIT('RSPINP, TOO MANY AVE OPERATORS SPECIFIED')
                  END IF
                  READ(LUCMD,'(A)') LBL1AV(N1AVTOT)
               GO TO 100
   16          CONTINUE ! .HIRPA
                  HIRPA = .TRUE.
                  SOPPA = .TRUE.
               GO TO 100
   17          CONTINUE ! .SOPW4
                  SOPW4 = .TRUE.
               GO TO 100
   18          CONTINUE ! .SOPPA(CCSD)
                  SOPPA = .TRUE.
                  CCPPA = .TRUE.
               GO TO 100
   19          CONTINUE ! .NODOIT
                  DIROIT = .FALSE.
               GO TO 100
 20            CONTINUE ! .DFT-SO
c              compute DFT_SO integrals
                  DFT_SO = .TRUE.
               GO TO 100
C              .PROP2Average
 21            CONTINUE
                  N2AVTOT = N2AVTOT + 1
                  IF ( N2AVTOT .GT. MAX2AV ) THEN
                   WRITE(LUPRI,'(//A,I5)')
     &               '@RSPINP: Number of operators for .PROP2AV '//
     &               'calculation exceeds the allowed number :',MAX2AV
                   CALL QUIT('TOO MANY .PROP2AV OPERATORS SPECIFIED')
                  END IF
                  READ(LUCMD,'(A)')LBL2AV(1,N2AVTOT),LBL2AV(2,N2AVTOT)
               GO TO 100
 22            CONTINUE ! .TDA
                  TDA = .TRUE.
               GO TO 100
 23            CONTINUE ! .CIS
                  CISRPA = .TRUE.
               GO TO 100
 24            CONTINUE ! .SOPRPA
                 SOPRPA = .TRUE.
                 SOPPA  = .TRUE.
                 WRITE (LUPRI,'(/A)')
     &           '@ **INFO**'//
     &            ' This is a SOPPA embedding calculation w/ dynamic'
                 WRITE (LUPRI,'(/A)')
     &           '@ **INFO**'//
     &            ' environment response at RPA level only,'
                 WRITE (LUPRI,'(/A)')
     &           '@ **INFO**'//
     &            ' cf. JCC, 2012-2022, 33, (2012) for reference.'
                 IF (SOPRPA .AND. (.NOT. EMBEDDING)) THEN
                  WRITE (LUPRI,'(/A)') '**RESPONSE INPUT ERROR: '//
     &             'Vacuum calculation - use .SOPPA instead of .SOPRPA'
                  CALL QUIT('RSPINP: Vacuum calc. - use .SOPPA instead')
                 END IF
               GO TO 100
 25            CONTINUE ! .TRDQF
                 IF(.NOT. QFIT) THEN
                   CALL QUIT(
     &             'RSPINP: Charge fitting with .TRDQF requires .QFIT')
                 ENDIF
                 TRDQF = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') '@ Prompt "',WORD,
     &            '" not recognized in RSPINP.'
               CALL QUIT('Illegal keyword prompt in RSPINP ')
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('CHANGES OF DEFAULTS FOR RSPINP:',0)
         IF (INPTES) WRITE(LUPRI,'(/A,L1)')
     *      ' INPUT TEST.  INPTES = ',INPTES
         IF (TDA) WRITE(LUPRI,'(/A,L1)')
     *      ' TAMM-DANCOFF APPROXIMATION (B = 0). TDA=',TDA
         IF (DODFT .AND. CISRPA)
     *      CALL QUIT('TDDFT CALC WITH .CIS KEYWORD - USE .TDA INSTEAD')
         IF (CISRPA) WRITE(LUPRI,'(/A,L1)')
     *      ' CI SINGLES CALCULATION (RPA/TDA (B = 0)). CIS=',CISRPA
         IF (HIRPA)  THEN
            WRITE(LUPRI,'(/A)')
     *   ' Higher RPA Polarization Propagator Approximation -- HRPA'
         ELSE IF (SOPPA .AND. .NOT. CCPPA) THEN
            WRITE(LUPRI,'(/A)')
     *   ' Second Order Polarization Propagator Approximation -- SOPPA'
         ELSE IF (SOPPA .AND. CCPPA) THEN
            WRITE(LUPRI,'(/2A/14X,A/A)')
     &   ' SOPPA(CCSD) :',
     &   ' Second Order Polarization Propagator Approximation',
     &   ' with Coupled Cluster Singles and Doubles Amplitudes',
     &   ' Reference : S. P. A. Sauer, '//
     &      'J. Phys. B: At. Mol. Opt. Phys. 30, 3773-3780 (1997)'
         END IF
         IF (DIRFCK) THEN
            WRITE(LUPRI,'(/A)')
     *      ' AO-direct Fock matrix calculations.'
         ELSE IF (NOITRA) THEN
            WRITE(LUPRI,'(/2A,L1)')
     *      ' MO integrals for RESPONS are already available.'
     *      ,'  NOITRA = ',NOITRA
         END IF
         IF (TRPFLG) WRITE(LUPRI,'(/2A,L1,A,L1)')
     *      ' ** CALCULATION OF TRIPLET RESPONSE PROPERTIES  **'
     *      ,'  TRPFLG = ',TRPFLG
     *      ,'  TRPLET = ',TRPLET
         IF (AVDIA) THEN
            WRITE(LUPRI,'(/2A/2A/A,L1)')
     *      ' Default : Using Fock type decoupling of the two-electron',
     *      ' density matrix :',
     *      '    Add DV*(FC+FV) instead of DV*FC to E[2] approximate',
     *      ' orbital diagonal'
         ELSE
            WRITE(LUPRI,'(/A)')
     *      ' .NOAVDIA : Using exact E[2] orbital diagonal. '
         END IF
         IF (ORBSFT.NE.1.0D-4) WRITE(LUPRI,'(/A,D13.6)')
     *      ' Diagonal orbital E[2] shifted with ORBSFT =',ORBSFT
         IF (ORBSPC) WRITE(LUPRI,'(/A)') ' .ORBSPC:'//
     *      ' Calculation with only orbital operators included.'
         IF (.NOT.NOS0MX) WRITE(LUPRI,'(/2A,L1)')
     *      ' S(0) sum rule calculated in mixed representation',
     *      ' NOS0MX = ',NOS0MX
         IF (OPTORB) WRITE(LUPRI,'(/2A,/A,L1)')
     *      ' .OPTORB: Orbital trial vectors calculated with optimal',
     *      ' orbital trial vector algorithm.'
         IF (MAXPHP.NE.100) WRITE(LUPRI,'(/A,I6)')
     *      ' Dimension of H(0) subspace. MAXPHP =',MAXPHP
         IF (PHPRES) WRITE(LUPRI,'(/2A,L1)')
     *      ' H(0) subspace selected from size of residual elements',
     *      ' PHPRES = ',PHPRES
         IF (ABSLRS) THEN
           IF (ABS_MAXRM.NE.600) WRITE(LUPRI,'(/A,I6)')
     *        ' MAX DIMENSION OF REDUCED SPACE. ABS_MAXRM =',ABS_MAXRM
         ELSE
           IF (MAXRM.NE.600) WRITE(LUPRI,'(/A,I6)')
     *        ' Max dimension of reduced space. MAXRM =',MAXRM
         ENDIF
         IF (QRREST) WRITE(LUPRI,'(/A/A)')
     *      ' Quadratic response calculation restarted',
     *      '    from previously calculated solution vectors'
         IF (.NOT. DIROIT) WRITE(LUPRI,'(/A/A,L1)')
     *      ' Direct linear transformations NOT used in RESPONSE',
     *      ' calculation, DIROIT=',DIROIT
         IF (DFT_SO) WRITE (LUPRI,'(/A)')
     &      ' .DFT-SO: DFT spin-orbit integrals calculated.'
      END IF
      IF (CCPPA .AND. RSPCI) THEN
         WRITE (LUPRI,'(/A)') ' Input ERROR:'//
     &      ' CI response and ".SOPPA(CCSD)" are incompatible.'
         CALL QUIT('**RESPONSE: Input error')
      ELSEIF (HIRPA .AND. RSPCI) THEN
         WRITE (LUPRI,'(/A)')
     &      ' Input ERROR: CI response and ".HIRPA " are incompatible.'
         CALL QUIT('**RESPONSE: Input error')
      ELSEIF (SOPPA .AND. RSPCI) THEN
         WRITE (LUPRI,'(/A)')
     &      ' Input ERROR: CI response and ".SOPPA " are incompatible.'
         CALL QUIT('**RESPONSE: Input error')
      END IF
C
C *** END OF RSPINP
C
      RETURN
      END
C  /* Deck bdinfrsp */
      BLOCK DATA RSPINF
C
C Initialize response common blocks
C 25-Jan-1988; l.r. 14-Sep-1990
C
#include "implicit.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "rspprp.h"
#include "inflr.h"
#include "infpp.h"
#include "infc6.h"
#include "infhyp.h"
#include "infsmo.h"
#include "infs0.h"
#include "trhso.h"
#include "maxorb.h"
#include "infvar.h"
#include "qrinf.h"
#include "infesr.h"
#include "infspi.h"
#include "infhso.h"
#include "infslv.h"
#include "infcr.h"
#include "inftmo.h"
#include "inftpa.h"
#include "mxcent.h"
#include "elweak.h"
#include "absorp.h"
#include "esg.h"
#include "channel.h"
#include "esrhfc.h"
#include "abslrs.h"
C
C     Initialization of RSPINP
      DATA INPTES/.FALSE./, AVDIA /.TRUE./
      DATA SOPPA /.FALSE./, RSPCI /.FALSE./, NOITRA/.FALSE./
      DATA CCPPA /.FALSE./
      DATA NOS0MX/.TRUE./
      DATA TRPLET/.FALSE./, TRPFLG/.FALSE./
      DATA ORBSPC/.FALSE./
      DATA ORBSFT/1.0D-4/
      DATA QRREST/.FALSE./
      DATA OPTORB/.FALSE./
      DATA MAXRM/600/
      DATA PHPRES/.FALSE./
      DATA E3TEST/.FALSE./
      DATA A2TEST/.FALSE./
      DATA X2TEST/.FALSE./
      DATA RSPECD/.FALSE./
      DATA RSPOCD/.FALSE./
      DATA TDA/.FALSE./
      DATA CISRPA/.FALSE./
      DATA SOPRPA/.FALSE./
C     Initialization of INFRSP for both PPINP and LRINP
      DATA MAXITO/5/, ANTTES/.FALSE./, IPRRSP/2/
C INITIALIZING ELWEAK
      DATA PVSO/.FALSE./, PVPSO/.FALSE./,ELWEAK/.FALSE./
C     Initialization of INFRSP for LRINP
      DATA ABOCHK/.FALSE./
      DATA ISTOCK/1/, MAXOCK/6/
      DATA NFREQ/1/, FREQ(1)/0.0D0/
C     Initialization of WRKRSP
      DATA RESTLR/.FALSE./, RESTC6/.FALSE./
      DATA RESTPP/.FALSE./, ABCHK/.FALSE./, ABSYM/.FALSE./
      DATA SOPRSY/.FALSE./
C     Initialization of INFPP
      DATA EXMOM/.FALSE./, EXMTES/.FALSE./
      DATA OLSEN/.FALSE./
      DATA MAXITP/60/, IPRPP/2/, THCPP/1.0D-3/, IPREXM/2/
C     Initialization of INFLR
      DATA MAXITL/60/, IPRLR/2/, THCLR/1.0D-3/, THRNRM/1.0D-9/
C     Initialization of INFC6
      DATA THCMOM/1.0D-1/
      DATA MAXITC/60/, IPRC6/2/, THCC6/1.0D-3/, MAXMOM/6/
      DATA NCFREQ/0/ , CFREQ(1)/0.0D0/
      DATA PADE/.FALSE./
      DATA GSLEGN/.FALSE./
      DATA NGRID/10/
      DATA C6IFC/.FALSE./
C     Initialization of TRHSO for HSOINP
      DATA OLDTRA/.FALSE./
C     Initialization of INFHSO for HSOINP
      DATA TESTZY/.FALSE./, DOSO1/.TRUE./, DOSO2/.TRUE./, IPRHSO/0/
      DATA X2MAT/.FALSE./, A2MAT/.FALSE./, X2GRAD/.FALSE./,
     &     PHOSPH/.FALSE./, MNFPHO/.FALSE./, ECPHOS/.FALSE./,           !MK
     &     PHOSPV/.FALSE./, CPPHOL/.FALSE./, CPPHOV/.FALSE./,           !MK
     &     CPPHMF/.FALSE./, CPPHEC/.FALSE./                             !MK
C     Initialization of INFHYP
      DATA HYPCAL/.FALSE./
      DATA QRSPEC/.FALSE./,QRSHG/.FALSE./,QRPOCK/.FALSE./
      DATA QROPRF/.FALSE./,SOCOLL/.FALSE./,SSCOLL/.FALSE./
      DATA IPRHYP/2/
      DATA NBQRFR/1/ , BQRFR(1)/0.0D0/
      DATA NCQRFR/1/ , CQRFR(1)/0.0D0/
      DATA REFCHK/.FALSE./
      DATA IAABB/0/
C     Initialization of INFSMO
      DATA SOMOM/.FALSE./ , TWOPHO/.FALSE./, MCDCAL/.FALSE./,
     &     MCDPRJ/.FALSE./, LTPCD/.FALSE./
      DATA IPRSMO/2/
      DATA NBSMFR/1/ , BSMFR(1)/0.0D0/
C     Initialization of INFESR
      DATA IPRESR/2/ , MAXESR/60/ ,THCESR/1.0D-5/ ,ESRCAL/.FALSE./
C     Initialization of INFSPI
      DATA ISPINA/0/,ISPINB/0/,ISPINC/0/,ISPIND/0/
      DATA MULSP/0,1,1,0/
      DATA DIROIT/.TRUE./
C  Initializing INFSLV
      DATA SLVANT/.FALSE./
      DATA IPRSLV/2/
C  Initializing INFCR
      DATA CRCAL/.FALSE./
      DATA IPRCR/2/
      DATA NBCRFR/1/, NCCRFR/1/, NDCRFR/1/
      DATA BCRFR(1)/0.0D0/, CCRFR(1)/0.0D0/, DCRFR(1)/0.0D0/
      DATA INVEXP/.FALSE./
      DATA CRTHG/.FALSE./,CRSHG/.FALSE./,CRKERR/.FALSE./,CRIDRI/.FALSE./
      DATA CRSPEC/.FALSE./
      DATA GAMALL/.TRUE./
C  Initializing INFTMO
      DATA TOMOM/.FALSE./,CTMOHG/.TRUE./
      DATA NBTMFR/1/, NCTMFR/1/
      DATA BTMFR(1)/0.0D0/, CTMFR(1)/0.0D0/
      DATA IPRTMO/2/
C  Initializing INFTPA
      DATA TPAMP/.FALSE./, TPALP/.FALSE./
      DATA NBTPFR/1/
      DATA BTPFR(1)/0.0D0/
      DATA IPRTPA/2/
C  Initializing ABSORP
      DATA ABSORP/.FALSE./,ABS_INTERVAL/.FALSE./,ABS_REDUCE/.FALSE./
      DATA ABS_ALPHA/.FALSE./,ABS_BETA/.FALSE./,ABS_GAMMA/.FALSE./
      DATA ABS_SHG/.FALSE./,ABS_MCD/.FALSE./,ABS_ANALYZE/.FALSE./
      DATA ABS_NGD/.FALSE./,ABS_GDCOMPLEX/.FALSE./
      DATA NFREQ_BATCH/10/,MAX_MACRO/30/,MAX_MICRO/10/,MAX_ITORB/5/
      DATA NEXCITED_STATES/2/,IPRABS/0/
      DATA NFREQ_ALPHA/1/,NFREQ_BETA_B/1/,NFREQ_BETA_C/1/
      DATA FREQ_INTERVAL/0.0D0,0.0D0,1.0D0/
      DATA THCLR_ABSORP/1.0D-3/,THCPP_ABSORP/1.0D-3/
      DATA FREQ_ALPHA(1)/0.0D0/
      DATA FREQ_BETA_B(1)/0.0D0/,FREQ_BETA_C(1)/0.0D0/
      DATA ABS_OLDCPP/.FALSE./,ABS_REBD/.TRUE./
C  Initializing CHANNEL
      DATA CHANNEL_CALC/.FALSE./,CHANNEL_VCALC/.FALSE./
      DATA CHANNEL_VIRT/-1,-1,-1,-1,-1,-1,-1,-1/
C  Initializing ESG
      DATA ESG/.FALSE./
      DATA IESG/1/
      DATA ISYME/1/
C  Initializing ESRHFC
      DATA ESRHFC/.FALSE./
      DATA HFCFC/.FALSE./,HFCSD/.FALSE./,HFCSO/.FALSE./
      DATA UNGAUSS/.TRUE./
      DATA HSOBRT/.FALSE./,HSOEFF/.FALSE./
      DATA IPRHFC/0/
C  Initializing ABSLRS
       DATA ABSLRS/.FALSE./,ABSLRS_INTERVAL/.FALSE./
       DATA ABSLRS_ALPHA/.FALSE./,ABS_BATCH/.FALSE./
       DATA ABSLRS_ANALYZE/.FALSE./,ABS_IMFREQ/.FALSE./
       DATA ABS_XREL/.FALSE./
       DATA ABS_MAXITER/150/,IPRABSLRS/0/,ABS_NFREQ_ALPHA/1/
       DATA ABS_NFREQ_BETA_B/1/,ABS_NFREQ_BETA_C/1/
       DATA ABS_MAXRM/1000/
       DATA ABS_FREQ_INTERVAL/0.0D0,0.0D0,1.0D0/,ABS_THD_LD/1.0D-8/
       DATA ABS_THCLR/1.0D-3/,ABS_FREQ_ALPHA(1)/0.0D0/
       DATA ABS_FREQ_BETA_B(1)/0.0D0/,ABS_FREQ_BETA_C(1)/0.0D0/
       DATA LUSB/80/,LUAB/81/,LUSS/82/,LUAS/83/,LUABSVECS/84/
       DATA LUE1RED/85/,LUE2RED/86/,LUSRED/87/
      END
C  /* Deck ppinp */
      SUBROUTINE PPINP(WORD)
C
C Read input options for PP = "Polarization Propagator" excitation energies
C                           = "Linear Response Single Residue"
C
C tbp, may 2005: keywords .ECD  (electronic circular dichroism) and
C                         .OECD (oriented electronic circular dichroism)
C                for linear response residues introduced.
C
      use pelib_interface, only: use_pelib
#include "implicit.h"
#include "mxcent.h"
#include "pcmdef.h"
      LOGICAL NEWDEF, ISLINR
      PARAMETER ( NTABLE = 55 )
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER*8 LABEL
      INTEGER NUMBER_ORBS(8)
C
#include "priunit.h"
#include "gnrinf.h"
C
#include "inforb.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "infpri.h"
#include "rspprp.h"
#include "infpp.h"
#include "infspi.h"
#include "infhso.h"
#include "inflr.h"
#include "inflin.h"
#include "maxorb.h"
#include "infinp.h"
#include "inftap.h"
#include "pcm.h"
#include "pcmlog.h"
#include "channel.h"
C
      DATA TABLE /'.DIPLEN', '.DIPLNX', '.DIPLNY', '.DIPLNZ','.DIPVEL',
     *            '.DIPVLX', '.DIPVLY', '.DIPVLZ', '.DIPMAG','.DIPMGX',
     *            '.DIPMGY', '.DIPMGZ', '.QUADMO', '.QUADXX','.QUADXY',
     *            '.QUADXZ', '.QUADYY', '.QUADYZ', '.QUADZZ','.SPIN-O',
     *            '.SPNORX', '.SPNORY', '.SPNORZ', '.PROPRT','.ROOTS ',
     *            '.MAX IT', '.THCPP ', '.OLSEN ', '.MAXITO','.SROOTS',
     *            '.ANTTES', '.RESTPP', '.ABCHK ', '.ABSYM ','.PRINT ',
     *            '.IPREXM', '.E3TEST', '.EXMTES', '.ISPABC','.ECD   ',
     *            '.X2TEST', '.A2TEST', '.SINGLE', '.DOUBLE','.THCLR ',
     *            '.MAXITL', '.TRIPLE', '.NSIMUL', '.NSTART','.OECD  ',
     *            '.ESA   ', '.CHANNE', '.VIRTUA', '.TDA   ','.CIS   '/
      INTEGER IXARR(MAXLPP)
C
C READ IN  INPUT
C
      ISLINR = WORD.EQ.'*LINEAR'
      NEWDEF = (ISLINR .OR. WORD.EQ.'*QUADRA')
      ICHANG = 0
      NPPTOT = 0
      IF (NEWDEF) THEN
         DO I=1,NSYM
            NPPCNV(I) = 1
         END DO
         NPPTOT = NSYM
         WORD1  = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') GO TO 100
            IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO I=1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
     *                      16,17,18,19,20,21,22,23,24,25,26,27,
     *                      28,29,30,31,32,33,34,35,36,37,38,39,
     *                      40,41,42,43,44,45,46,47,48,49,50,51,
     *                      52,53,54,55), I
                  END IF
               END DO
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/4A/)') '@ Keyword "',WORD,
     *            '" is not recognized in RSPPP for section ',WORD1
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT(' Illegal KEYWORD in RSPPP ')
 1             CONTINUE   ! .DIPLEN
               LPPOP( INDPRP('XDIPLEN ')) = .TRUE.
               LPPOP( INDPRP('YDIPLEN ')) = .TRUE.
               LPPOP( INDPRP('ZDIPLEN ')) = .TRUE.
               GO TO 100
 2             CONTINUE   ! .DIPLNX
               LPPOP( INDPRP('XDIPLEN ')) = .TRUE.
               GO TO 100
 3             CONTINUE   ! .DIPLNY
               LPPOP( INDPRP('YDIPLEN ')) = .TRUE.
               GO TO 100
 4             CONTINUE   ! .DIPLNZ
               LPPOP( INDPRP('ZDIPLEN ')) = .TRUE.
               GO TO 100
 5             CONTINUE   ! .DIPVEL
               LPPOP( INDPRP('XDIPVEL ')) = .TRUE.
               LPPOP( INDPRP('YDIPVEL ')) = .TRUE.
               LPPOP( INDPRP('ZDIPVEL ')) = .TRUE.
               GO TO 100
 6             CONTINUE   ! .DIPVLX
                  LPPOP( INDPRP('XDIPVEL ')) = .TRUE.
               GO TO 100
 7             CONTINUE   ! .DIPVLY
               LPPOP( INDPRP('YDIPVEL ')) = .TRUE.
               GO TO 100
 8             CONTINUE   ! .DIPVLZ
               LPPOP( INDPRP('ZDIPVEL ')) = .TRUE.
               GO TO 100
 9             CONTINUE   ! .DIPMAG
               LPPOP( INDPRP('XANGMOM ')) = .TRUE.
               LPPOP( INDPRP('YANGMOM ')) = .TRUE.
               LPPOP( INDPRP('ZANGMOM ')) = .TRUE.
               GO TO 100
 10            CONTINUE   ! .DIPMGX
               LPPOP( INDPRP('XANGMOM ')) = .TRUE.
               GO TO 100
 11            CONTINUE   ! .DIPMGY
               LPPOP( INDPRP('YANGMOM ')) = .TRUE.
               GO TO 100
 12            CONTINUE   ! .DIPMGZ
               LPPOP( INDPRP('ZANGMOM ')) = .TRUE.
               GO TO 100
 13            CONTINUE   ! .QUADMO
               LPPOP( INDPRP('XXQUADRU')) = .TRUE.
               LPPOP( INDPRP('XYQUADRU')) = .TRUE.
               LPPOP( INDPRP('XZQUADRU')) = .TRUE.
               LPPOP( INDPRP('YYQUADRU')) = .TRUE.
               LPPOP( INDPRP('YZQUADRU')) = .TRUE.
               GO TO 100
 14            CONTINUE   ! .QUADXX
               LPPOP( INDPRP('XXQUADRU')) = .TRUE.
               GO TO 100
 15            CONTINUE   ! .QUADXY
               LPPOP( INDPRP('XYQUADRU')) = .TRUE.
               GO TO 100
 16            CONTINUE   ! .QUADXZ
               LPPOP( INDPRP('XZQUADRU')) = .TRUE.
               GO TO 100
 17            CONTINUE   ! .QUADYY
               LPPOP( INDPRP('YYQUADRU')) = .TRUE.
               GO TO 100
 18            CONTINUE   ! .QUADYZ
               LPPOP( INDPRP('YZQUADRU')) = .TRUE.
               GO TO 100
 19            CONTINUE   ! .QUADZZ
               LPPOP( INDPRP('ZZQUADRU')) = .TRUE.
               GO TO 100
 20            CONTINUE   ! .SPIN-ORB
               LPPOP( INDPRP('X SPNORB')) = .TRUE.
               LPPOP( INDPRP('Y SPNORB')) = .TRUE.
               LPPOP( INDPRP('Z SPNORB')) = .TRUE.
               DOSO1 = .TRUE.
               DOSO2 = .TRUE.
               TRPLET = .TRUE.
               TRPFLG = .TRUE.
               GO TO 100
 21            CONTINUE   ! .SPNORX
               LPPOP( INDPRP('X SPNORB')) = .TRUE.
               DOSO1 = .TRUE.
               DOSO2 = .TRUE.
               TRPLET = .TRUE.
               TRPFLG = .TRUE.
               GO TO 100
 22            CONTINUE   ! .SPNORY
               LPPOP( INDPRP('Y SPNORB')) = .TRUE.
               DOSO1 = .TRUE.
               DOSO2 = .TRUE.
               TRPLET = .TRUE.
               TRPFLG = .TRUE.
               GO TO 100
 23            CONTINUE   ! .SPNORZ
               LPPOP( INDPRP('Z SPNORB')) = .TRUE.
               DOSO1 = .TRUE.
               DOSO2 = .TRUE.
               TRPLET = .TRUE.
               TRPFLG = .TRUE.
               GO TO 100
 24            CONTINUE   ! .PROPRT
               READ(LUCMD,'( A )')LABEL
               LPPOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 25            CONTINUE   ! .ROOTS
               NPPTOT = 0
               READ (LUCMD,*) (NPPCNV(MULD2H(J,LSYMRF)),J=1,NSYM)
               DO J=1,NSYM
                  NPPTOT = NPPTOT + NPPCNV(J)
               END DO
               GO TO 100
 26            CONTINUE   ! .MAX IT
               READ (LUCMD,*) MAXITP
               GO TO 100
 27            CONTINUE   ! .THCPP
               READ (LUCMD,*) THCPP
               GO TO 100
 28            CONTINUE   ! .OLSEN
               OLSEN  = .TRUE.
               GO TO 100
 29            CONTINUE   ! .MAXITO
               READ (LUCMD,*) MAXITO
               GO TO 100
 30            CONTINUE   ! .SROOTS
c              SROOTS: compute transition moments between
c              selected excited states
               NPPTOT = 0
               WRITE(LUPRI,*)'Will compute transition moments between ',
     *              ' selected roots'
               WRITE(LUPRI,*) 'Symmetry : Roots '
               DO JSYM=1,NSYM
                  READ (LUCMD,*) IXL, (IXARR(J), J=1,IXL)
                  IMAX = 0
                  DO J = 1, MAXLPP
                     LCMEXC(J,MULD2H(JSYM,LSYMRF)) = .FALSE.
                  END DO
                  DO J = 1, IXL
                     LCMEXC(IXARR(J),MULD2H(JSYM,LSYMRF)) = .TRUE.
                     IMAX = MAX(IMAX,IXARR(J))
                  END DO
                  WRITE(LUPRI,'(I9," : ",20I3)')
     &                 JSYM,(IXARR(J),J=1,IXL)
                  NPPCNV(MULD2H(JSYM,LSYMRF)) = IMAX + 1
               END DO
               DO JSYM=1,NSYM
                  NPPTOT = NPPTOT + NPPCNV(JSYM)
               END DO
               GO TO 100
 31            CONTINUE   ! .ANTTES
               ANTTES = .TRUE.
               GO TO 100
 32            CONTINUE   ! .RESTPP
               RESTPP = .TRUE.
               GO TO 100
 33            CONTINUE   ! .ABCHK
               ABCHK = .TRUE.
               GO TO 100
 34            CONTINUE   ! .ABSYM
               ABSYM = .TRUE.
               GO TO 100
 35            CONTINUE   ! .PRINT
               READ (LUCMD,*) IPRPP
               GO TO 100
 36            CONTINUE   ! .IPREXM
               READ (LUCMD,*) IPREXM
               GO TO 100
 37            CONTINUE   ! .E3TEST
               E3TEST = .TRUE.
               GO TO 100
 38            CONTINUE   ! .EXMTES
               EXMTES = .TRUE.
               GO TO 100
 39            CONTINUE   ! .ISPABC
               READ (LUCMD,*) ISPINA,ISPINB,ISPINC
               GO TO 100
 40            CONTINUE   ! .ECD
               IF (.NOT. RSPECD) THEN
                  RSPECD = .TRUE. ! electronic CD
                  LPPOP( INDPRP('XDIPLEN ')) = .TRUE.
                  LPPOP( INDPRP('YDIPLEN ')) = .TRUE.
                  LPPOP( INDPRP('ZDIPLEN ')) = .TRUE.
                  LPPOP( INDPRP('XDIPVEL ')) = .TRUE.
                  LPPOP( INDPRP('YDIPVEL ')) = .TRUE.
                  LPPOP( INDPRP('ZDIPVEL ')) = .TRUE.
                  LPPOP( INDPRP('XANGMOM ')) = .TRUE.
                  LPPOP( INDPRP('YANGMOM ')) = .TRUE.
                  LPPOP( INDPRP('ZANGMOM ')) = .TRUE.
               END IF
               GO TO 100
 41            CONTINUE   ! .X2TEST
               X2TEST = .TRUE.
               GO TO 100
 42            CONTINUE   ! .A2TEST
               A2TEST = .TRUE.
               GO TO 100
 43            CONTINUE   ! .SINGLE
               IF (EXMOM) THEN
                 CALL QUIT('PPINP implementation ERROR: please report')
               END IF
               GO TO 100
 44            CONTINUE   ! .DOUBLE
               IF (.NOT.EXMOM) THEN
                 CALL QUIT('PPINP implementation ERROR: please report')
               END IF
               GO TO 100
 45            CONTINUE   ! .THCLR
                  READ (LUCMD,*) THCLR
               GO TO 100
 46            CONTINUE   ! .MAXITL
                  READ (LUCMD,*) MAXITL
               GO TO 100
 47            CONTINUE   ! .TRIPLET
                  TRPLET = .TRUE.
                  TRPFLG = .TRUE.
               GO TO 100
 48            CONTINUE   ! .NSUMUL
                  READ (LUCMD,*) (NPPSIM(J),J=1,NSYM)
               GO TO 100
 49            CONTINUE   ! .NSSTART
                  READ (LUCMD,*) (NPPSTV(J),J=1,NSYM)
               GO TO 100
 50            CONTINUE   ! .OECD
               IF (.NOT. RSPOCD) THEN
                  RSPOCD = .TRUE. ! oriented electronic CD
                  LPPOP( INDPRP('XXSECMOM')) = .TRUE.
                  LPPOP( INDPRP('XYSECMOM')) = .TRUE.
                  LPPOP( INDPRP('XZSECMOM')) = .TRUE.
                  LPPOP( INDPRP('YYSECMOM')) = .TRUE.
                  LPPOP( INDPRP('YZSECMOM')) = .TRUE.
                  LPPOP( INDPRP('ZZSECMOM')) = .TRUE.
                  LPPOP( INDPRP('XXROTSTR')) = .TRUE.
                  LPPOP( INDPRP('XYROTSTR')) = .TRUE.
                  LPPOP( INDPRP('XZROTSTR')) = .TRUE.
                  LPPOP( INDPRP('YYROTSTR')) = .TRUE.
                  LPPOP( INDPRP('YZROTSTR')) = .TRUE.
                  LPPOP( INDPRP('ZZROTSTR')) = .TRUE.
                  IF (.NOT.RSPECD) THEN
                     RSPECD = .TRUE. ! electronic CD
                     LPPOP( INDPRP('XDIPLEN ')) = .TRUE.
                     LPPOP( INDPRP('YDIPLEN ')) = .TRUE.
                     LPPOP( INDPRP('ZDIPLEN ')) = .TRUE.
                     LPPOP( INDPRP('XDIPVEL ')) = .TRUE.
                     LPPOP( INDPRP('YDIPVEL ')) = .TRUE.
                     LPPOP( INDPRP('ZDIPVEL ')) = .TRUE.
                     LPPOP( INDPRP('XANGMOM ')) = .TRUE.
                     LPPOP( INDPRP('YANGMOM ')) = .TRUE.
                     LPPOP( INDPRP('ZANGMOM ')) = .TRUE.
                  END IF
               END IF
               GO TO 100
 51            CONTINUE   ! .ESA
C panor:
C     Excited state absorption is calculated from the lowest excited
C     state of symmetry ESASYM to other excited states specified on
C     subsequent line.
               DOESA = .TRUE.
               READ (LUCMD,*) ESASYM
               LPPOP( INDPRP('XDIPLEN ')) = .TRUE.
               LPPOP( INDPRP('YDIPLEN ')) = .TRUE.
               LPPOP( INDPRP('ZDIPLEN ')) = .TRUE.
               NPPTOT = 0
               READ (LUCMD,*) (NPPCNV(J),J=1,NSYM)
               IF (NPPCNV(ESASYM).EQ.0) NPPCNV(ESASYM)=1
               DO J=1,NSYM
                  NPPTOT = NPPTOT + NPPCNV(J)
               END DO
               GO TO 100
C .CHANNEL - restrict to excitations from particular orbitals
 52            CONTINUE   ! .CHANNEL
               CHANNEL_CALC = .TRUE.
               CHANNEL_NORB = 0
               READ(LUCMD,*) (NUMBER_ORBS(I), I=1,NSYM)
               DO I=1,NSYM
                  CHANNEL_NORB = CHANNEL_NORB + NUMBER_ORBS(I)
               END DO
               IF (CHANNEL_NORB .GT. CHANNEL_MAX) CALL QUIT(
     &               'Too many orbitals for .CHANNEL calculation: ' //
     &               'increase CHANNEL_MAX and recompile')
               IOFF=0
               DO ISYM=1,NSYM
                  IF (NUMBER_ORBS(ISYM).GT.0) THEN
                     READ(LUCMD,*) (CHANNEL_ORB(I+IOFF),
     &                  I=1,NUMBER_ORBS(ISYM))
                     DO I=1,NUMBER_ORBS(ISYM)
                        CHANNEL_ORB(I+IOFF) =
     &                  CHANNEL_ORB(I+IOFF) + IORB(ISYM)
                     END DO
                     IOFF = IOFF + NUMBER_ORBS(ISYM)
                  ELSE
                     READ(LUCMD,*)
                  END IF
               END DO
               GOTO 100
C .VIRTUAL - restrict to excitations to lowest virtuals in each symmetry
 53            CONTINUE   ! .VIRTUAL
               CHANNEL_VCALC = .TRUE.
               READ (LUCMD,*) (CHANNEL_VIRT(I),I=1,NSYM)
               GOTO 100
 54            CONTINUE   ! .TDA
               TDA = .TRUE.
               GO TO 100
 55            CONTINUE   ! .CIS
               CISRPA = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') '@PROMPT "',WORD,
     *            '" NOT RECOGNIZED IN PPINP.'
               CALL QUIT(' ILLEGAL PROMPT IN PPINP ')
            END IF
      END IF   ! for IF (NEWDEF)
C ------------------
  300 CONTINUE
      IF (THR_REDFAC .GT. 0.0D0) THEN
         ICHANG = ICHANG + 1
         WRITE (LUPRI,'(3A,1P,D10.2)') '@ INFO ',WORD1,
     &   ' thresholds multiplied with general factor',THR_REDFAC
         THCLR = THCLR*THR_REDFAC
         THCPP = THCPP*THR_REDFAC
      END IF
C
      IF (NPPTOT .EQ. 0) THEN
         NINFO = NINFO + 1
         WRITE (LUPRI,'(/3A)') '@INFO: ',WORD1,
     *      ' input ignored because no roots requested.'
      ELSE
         IF (EXMOM) THEN
         CALL HEADER(' Quadratic Response double residue calculation',0)
         ELSE
            CALL HEADER(' Linear Response single residue calculation',0)
         END IF
         IF (FLAG(16) .AND. .NOT.INERSI) THEN
            WRITE(LUPRI,'(A/A,F8.4)')
     *      ' Equilibrium continuum solvent model requested.',
     *      ' Dielectric constant                            : EPSOL  ='
     *      ,EPSOL
         END IF
         IF (FLAG(16) .AND. INERSI) THEN
            WRITE(LUPRI,'(A,2(/A,F8.4))')
     *      ' Non-equilibrium continuum solvent model requested.',
     *      ' Static dielectric constant                     : EPSTAT ='
     *      ,EPSTAT,
     *      ' Optical dielectric constant                    : EPSOL  ='
     *      ,EPSOL
         END IF

         IF (QM3 .OR. QMMM .OR. QMNPMM) THEN
            WRITE(LUPRI,'(A)') ' QM/MM response calculation'
         ENDIF

         IF (USE_PELIB()) THEN
           WRITE(LUPRI,'(A)') ' Environment is modeled using'//
     &                  ' the polarizable embedding scheme (PE library)'
         END IF

         IF (PCM .AND. .NOT.NEQRSP) THEN
            WRITE(LUPRI,'(A,L1,/)')
     *      ' Equilibrium PCM solvent model requested        : SOLVNT ='
     *      ,PCM
            WRITE(LUPRI,'(A,F8.4)')
     *      ' Dielectric constant                            : EPSOL  ='
     *      ,EPS
         END IF
         IF (PCM .AND. NEQRSP) THEN
            WRITE(LUPRI,'(A,L1,/)')
     *      ' Non-equilibrium PCM solvent model requested    : INERSI ='
     *      ,PCM
            WRITE(LUPRI,'(A,F8.4)')
     *      ' Static dielectric constant                     : EPSTAT ='
     *      ,EPS
            WRITE(LUPRI,'(A,F8.4)')
     *      ' Optical dielectric constant                    : EPSOL  ='
     *      ,EPSINF
         END IF
C
         IF (ICHANG .GT. 0) THEN
            WRITE(LUPRI,'(/I5,A/)') ICHANG,' input options by user.'
         END IF

         IF (DOESA) THEN
            WRITE(LUPRI,'(A,L1,/)')
     *      ' Excited state absorption calculattion requested: DOESA  ='
     *      ,DOESA
         END IF
         IF (CHANNEL_CALC) THEN
            WRITE(LUPRI,'(/A,15(/,10I5))')
     &         ' Restricted channel RPA,   excitations from orbs:',
     &           (CHANNEL_ORB(I),I=1,CHANNEL_NORB)
         END IF
         IF (CHANNEL_VCALC) THEN
          DO I=1,NSYM
             WRITE(LUPRI,'(/A,1I5,/A,1I5)') 'Including only the first',
     &           CHANNEL_VIRT(I),' virtual orbitals in symmetry',I
          ENDDO
         ENDIF
         WRITE(LUPRI,'(A,I4)')
     *      ' Print level                                    : IPRPP  ='
     *       ,IPRPP
         WRITE(LUPRI,'(A,I4)')
     *      ' Maximum number of iterations for eigenval.eqs. : MAXITP ='
     *       ,MAXITP
         WRITE(LUPRI,'(A,1P,D10.3)')
     *      ' Threshold for convergence of eigenvalue eqs.   : THCPP  ='
     *      ,THCPP
         IF (EXMOM) THEN
            WRITE(LUPRI,'(A,I4)')
     *      ' Maximum number of iterations lin.eqs. for A    : MAXITL ='
     *       ,MAXITL
            WRITE(LUPRI,'(A,1P,D10.3)')
     *      ' Threshold for convergence of lin.eqs. for A    : THCLR  ='
     *      ,THCLR
         END IF
         WRITE(LUPRI,'(A,I4)')
     *      ' Maximum iterations in optimal orbital algorithm: MAXITO ='
     *      ,MAXITO
         IF ( ABCHK.AND.ABSYM ) THEN
            ABCHK =.FALSE.
            NWARN = NWARN + 1
            WRITE(LUPRI,'(/A,/A)')'@*******WARNING***********',
     *      ' ABCHK AND ABSYM CANNOT BOTH BE TRUE. ABCHK RESET TO FALSE'
         END IF
         IF (ANTTES) WRITE(LUPRI,'(/A,L1)')
     *      ' ANTISYMMETRY TEST OPTION. ANTTES =',ANTTES
         IF (RESTPP) WRITE(LUPRI,'(/A,L1)')
     *      ' USE TRIAL VECTORS AVAILABLE FOR RESTART. RESTPP =',RESTPP
         IF (ABCHK) WRITE(LUPRI,'(/A,L1)')
     *      ' SET UP E(2) AND S(2) EXPLICITLY. ABCHK =',ABCHK
         IF (ABSYM) WRITE(LUPRI,'(/A,L1)')
     *      ' CHECK SYMMETRY OF E(2) AND S(2) IN REDUCED SPACE. ABSYM ='
     *      ,ABSYM
         IF (OLSEN) WRITE(LUPRI,'(/A,L1)')
     *      ' CI TRIAL VECTORS WITH OLSEN ALGORITHM. OLSEN=',OLSEN
         IF (TDA) WRITE(LUPRI,'(/A,L1)')
     *      ' TAMM-DANCOFF APPROXIMATION (B = 0). TDA=',TDA
         IF (DODFT .AND. CISRPA)
     *      CALL QUIT('TDDFT CALC WITH .CIS KEYWORD - USE .TDA INSTEAD')
         IF (CISRPA) WRITE(LUPRI,'(/A,L1)')
     *      ' CI SINGLES CALCULATION (RPA/TDA (B = 0)). CIS=',CISRPA
       IF (EXMOM) THEN
         IF (IPREXM.NE.2) WRITE(LUPRI,'(/A,I5)')
     *        '    PRINT LEVEL IN EXM ROUTINES. IPREXM =',IPREXM
         IF (EXMTES) WRITE (LUPRI,'(/A)')
     *        '    Test that <i/A/j> = <j/A/i>  (EXMTES true).'
         IF (E3TEST) WRITE (LUPRI,'(/A)')
     *         '    Test of contributions from E3 and S3 terms.'
         WRITE (LUPRI,'(/A,I5)')' Spin of operator A , ISPINA=',ISPINA
         WRITE (LUPRI,'( A,I5)')' Spin of operator B ,'//
     *        ' (Excitation energy) ISPINB=',ISPINB
         WRITE (LUPRI,'( A,I5)') ' Spin of operator C ,'//
     *        ' (Excitation energy) ISPINC=',ISPINC
         IF (ISPINA+ISPINB+ISPINC .GT. 0) THEN
            TRPLET = .TRUE.
            TRPFLG = .TRUE.
         END IF
       END IF
      END IF
C
C  *** END OF PPINP
C
      RETURN
      END
C  /* Deck lrinp */
      SUBROUTINE LRINP(WORD)
C
C Read input options for LR = "Linear Response"
C
C
      use pelib_interface, only: use_pelib
#include "implicit.h"
C
#include "priunit.h"
#include "mxcent.h"
#include "pcmdef.h"
#include "infrsp.h"
#include "inforb.h"
#include "wrkrsp.h"
#include "rspprp.h"
#include "inflr.h"
#include "infpri.h"
#include "infhso.h"
#include "pcm.h"
#include "pcmlog.h"
#include "maxorb.h"
#include "infinp.h"
#include "inftap.h"
#include "nuclei.h"
#include "gnrinf.h"
#include "chrnos.h"
#include "elweak.h"
#include "infrank.h"
#include "channel.h"
C--------------------
C
      LOGICAL NEWDEF, SPNDEF, NXTLAB
      PARAMETER ( NTABLE = 48 )
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER*8 LABEL, LABFND
      INTEGER NUMBER_ORBS(8)
C
!                    1          2          3          4
      DATA TABLE /'.DIPLEN', '.DIPLNX', '.DIPLNY', '.DIPLNZ',
!                    5          6          7          8
     &            '.DIPVEL', '.DIPVLX', '.DIPVLY', '.DIPVLZ',
!                    9         10         11         12
     &            '.DIPMAG', '.DIPMGX', '.DIPMGY', '.DIPMGZ',
!                   13         14         15         16
     &            '.QUADMO', '.QUADXX', '.QUADXY', '.QUADXZ',
!                   17         18         19         20
     &            '.QUADYY', '.QUADYZ', '.QUADZZ', '.SPIN-O',
!                   21         22         23         24
     &            '.SPNORX', '.SPNORY', '.SPNORZ', '.PROPRT',
!                   25         26         27         28
     &            '.MAX IT', '.THCLR ', '.MAXITO', '.PRINT ',
!                   29         30         31         32
     &            '.FREQUE', '.ABOCHK', '.ISTOCK', '.MAXOCK',
!                   33         34         35         36
     &            '.SOPRSY', '.RESTLR', '.DIROIT', '.TRIPLE',
!                   37         38         39         40
     &            '.PV PSO', '.PV SO ', '.PV SO2', '.FERMI ',
!                   41         42         43         44
     &            '.SPIN-D', '.DIPLF ', '.DIPLFX', '.DIPLFY',
!                   45         46         47         48
     &            '.DIPLFZ', '.CHANNE', '.VIRTUA', '.xxxxxx'/
C
C
C
C READ IN  INPUT
C
      NEWDEF = (WORD .EQ. '*LINEAR')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') GO TO 100
            IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
     *                      16,17,18,19,20,21,22,23,24,25,26,27,
     *                      28,29,30,31,32,33,34,35,36,37,38,39,40,41,
     *                      42,43,44,45,46,47,48),I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') '@KEYWORD "',WORD,
     *            '" NOT RECOGNIZED IN LRINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT(' ILLEGAL KEYWORD IN LRINP ')
 1             CONTINUE
                  LLROP( INDPRP('XDIPLEN ')) = .TRUE.
                  LLROP( INDPRP('YDIPLEN ')) = .TRUE.
                  LLROP( INDPRP('ZDIPLEN ')) = .TRUE.
               GO TO 100
 2             CONTINUE
                  LLROP( INDPRP('XDIPLEN ')) = .TRUE.
               GO TO 100
 3             CONTINUE
                  LLROP( INDPRP('YDIPLEN ')) = .TRUE.
               GO TO 100
 4             CONTINUE
                  LLROP( INDPRP('ZDIPLEN ')) = .TRUE.
               GO TO 100
 5             CONTINUE
                  LLROP( INDPRP('XDIPVEL ')) = .TRUE.
                  LLROP( INDPRP('YDIPVEL ')) = .TRUE.
                  LLROP( INDPRP('ZDIPVEL ')) = .TRUE.
               GO TO 100
 6             CONTINUE
                  LLROP( INDPRP('XDIPVEL ')) = .TRUE.
               GO TO 100
 7             CONTINUE
                  LLROP( INDPRP('YDIPVEL ')) = .TRUE.
               GO TO 100
 8             CONTINUE
                  LLROP( INDPRP('ZDIPVEL ')) = .TRUE.
               GO TO 100
 9             CONTINUE
                  LLROP( INDPRP('XANGMOM ')) = .TRUE.
                  LLROP( INDPRP('YANGMOM ')) = .TRUE.
                  LLROP( INDPRP('ZANGMOM ')) = .TRUE.
               GO TO 100
 10            CONTINUE
                  LLROP( INDPRP('XANGMOM ')) = .TRUE.
               GO TO 100
 11            CONTINUE
                  LLROP( INDPRP('YANGMOM ')) = .TRUE.
               GO TO 100
 12            CONTINUE
                  LLROP( INDPRP('ZANGMOM ')) = .TRUE.
               GO TO 100
 13            CONTINUE
                  LLROP( INDPRP('XXQUADRU')) = .TRUE.
                  LLROP( INDPRP('XYQUADRU')) = .TRUE.
                  LLROP( INDPRP('XZQUADRU')) = .TRUE.
                  LLROP( INDPRP('YYQUADRU')) = .TRUE.
                  LLROP( INDPRP('YZQUADRU')) = .TRUE.
               GO TO 100
 14            CONTINUE
                  LLROP( INDPRP('XXQUADRU')) = .TRUE.
               GO TO 100
 15            CONTINUE
                  LLROP( INDPRP('XYQUADRU')) = .TRUE.
               GO TO 100
 16            CONTINUE
                  LLROP( INDPRP('XZQUADRU')) = .TRUE.
               GO TO 100
 17            CONTINUE
                  LLROP( INDPRP('YYQUADRU')) = .TRUE.
               GO TO 100
 18            CONTINUE
                  LLROP( INDPRP('YZQUADRU')) = .TRUE.
               GO TO 100
 19            CONTINUE
                  LLROP( INDPRP('ZZQUADRU')) = .TRUE.
               GO TO 100
 20            CONTINUE
                  LLROP( INDPRP('X SPNORB')) = .TRUE.
                  LLROP( INDPRP('Y SPNORB')) = .TRUE.
                  LLROP( INDPRP('Z SPNORB')) = .TRUE.
                  OPRANK( INDPRP('X SPNORB')) = 1
                  OPRANK( INDPRP('Y SPNORB')) = 1
                  OPRANK( INDPRP('Z SPNORB')) = 1
                  DOSO1 = .TRUE.
                  DOSO2 = .TRUE.
               GO TO 100
 21            CONTINUE
                  LLROP( INDPRP('X SPNORB')) = .TRUE.
                  OPRANK( INDPRP('X SPNORB')) = 1
                  DOSO1 = .TRUE.
                  DOSO2 = .TRUE.
               GO TO 100
 22            CONTINUE
                  LLROP( INDPRP('Y SPNORB')) = .TRUE.
                  OPRANK( INDPRP('Y SPNORB')) = 1
                  DOSO1 = .TRUE.
                  DOSO2 = .TRUE.
               GO TO 100
 23            CONTINUE
                  LLROP( INDPRP('Z SPNORB')) = .TRUE.
                  OPRANK( INDPRP('Z SPNORB')) = 1
                  DOSO1 = .TRUE.
                  DOSO2 = .TRUE.
               GO TO 100
 24            CONTINUE
                  READ(LUCMD,'(BN,A,I8)')LABEL,IRANK
                  LLROP( INDPRP(LABEL)) = .TRUE.
                  OPRANK(INDPRP(LABEL)) = IRANK
               GO TO 100
 25            CONTINUE
                  READ (LUCMD,*) MAXITL
               GO TO 100
 26            CONTINUE
                  READ (LUCMD,*) THCLR
               GO TO 100
 27            CONTINUE
                  READ (LUCMD,*) MAXITO
               GO TO 100
 28            CONTINUE
                  READ (LUCMD,*) IPRLR
               GO TO 100
 29            CONTINUE
                  READ (LUCMD,*) NFREQ
                  IF (NFREQ.LE.MFREQ) THEN
                     READ (LUCMD,*) (FREQ(J),J=1,NFREQ)
                  ELSE
                     WRITE (LUPRI,'(3(/,A,I5),/)')
     *               '@NUMBER OF FREQUENCIES SPECIFIED    :',NFREQ,
     *               '@IS GREATER THAN THE ALLOWED NUMBER :',MFREQ,
     *               '@THE NUMBER IS RESET TO THE MAXIMUM :',MFREQ
                     READ (LUCMD,*) (FREQ(J),J=1,MFREQ),
     *                              (FFFF,J=MFREQ+1,NFREQ)
                     NFREQ = MFREQ
                  END IF
               GO TO 100
 30            CONTINUE
                  ABOCHK =.TRUE.
               GO TO 100
 31            CONTINUE
                  READ (LUCMD,*) ISTOCK
               GO TO 100
 32            CONTINUE
                  READ (LUCMD,*) MAXOCK
               GO TO 100
 33            CONTINUE
                  SOPRSY=.TRUE.
               GO TO 100
 34            CONTINUE
                  RESTLR=.TRUE.
               GO TO 100
 35            CONTINUE
                  DIROIT=.TRUE.
               GO TO 100
 36            CONTINUE
                  TRPLET=.TRUE.
                  TRPFLG=.TRUE.
               GO TO 100
c ach
 37            CONTINUE
                  TRPLET=.TRUE.
                  TRPFLG=.TRUE.
                  ELWEAK=.TRUE.
                  PVPSO =.TRUE.
C PV from pso

                  LLROP( INDPRP('PVIOLA X')) = .TRUE.
                  LLROP( INDPRP('PVIOLA Y')) = .TRUE.
                  LLROP( INDPRP('PVIOLA Z')) = .TRUE.
C for pv  er det aldri symetri ?
C                  if (NSYM .ne 1 ) call QUIT('PV-cannot be calculated'
C#     &                 // ' with symmetry' )
C
                  ILAB=0
                  DO  IATOM = 1, Natoms
                     DO  ICOOR = 1, 3
                        ILAB = ILAB+1
                        IFIRST = ILAB/100
                        ISECND = MOD(ILAB,100)/10
                        ITHIRD = MOD(MOD(ILAB,100),10)
                        LLROP( INDPRP(
     &                       'PSO '//CHRNOS(IFIRST)//CHRNOS(ISECND)
     &                       //CHRNOS(ITHIRD)//' '))=.TRUE.
                     enddo
                  enddo
               GO TO 100
 38            CONTINUE
C PV from spin-orbit
                  PVSO =.TRUE.
                  TRPLET=.TRUE.
                  TRPFLG=.TRUE.
                  ELWEAK = .TRUE.
C
                  LLROP( INDPRP('PVIOLA X')) = .TRUE.
                  LLROP( INDPRP('PVIOLA Y')) = .TRUE.
                  LLROP( INDPRP('PVIOLA Z')) = .TRUE.
C
                  LLROP( INDPRP('X1SPNORB')) = .TRUE.
                  LLROP( INDPRP('Y1SPNORB')) = .TRUE.
                  LLROP( INDPRP('Z1SPNORB')) = .TRUE.
               GO TO 100
C
 39            CONTINUE
C PV from spin-orbit 2 el
               PVSO2 =.TRUE.
               TRPLET=.TRUE.
               TRPFLG=.TRUE.
               ELWEAK = .TRUE.

               LLROP( INDPRP('PVIOLA X')) = .TRUE.
               LLROP( INDPRP('PVIOLA Y')) = .TRUE.
               LLROP( INDPRP('PVIOLA Z')) = .TRUE.
               DOSO1 = .TRUE.
               DOSO2 = .TRUE.
               LLROP( INDPRP('X2SPNORB')) = .TRUE.
               LLROP( INDPRP('Y2SPNORB')) = .TRUE.
               LLROP( INDPRP('Z2SPNORB')) = .TRUE.
               GO TO 100
C
 40            CONTINUE
               CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',
     &                     'UNFORMATTED',IDUMMY,.FALSE.)
 199           CONTINUE
               IF (NXTLAB(LABFND,RTNLBL,LUPROP)) THEN
                  IF (LABFND(1:3) .EQ. 'FC ')
     &                 LLROP(INDPRP(LABFND)) = .TRUE.
                  GOTO 199
               END IF
               CALL GPCLOSE(LUPROP,'KEEP')
               GO TO 100
C
 41            CONTINUE
               CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',
     &                     'UNFORMATTED',IDUMMY,.FALSE.)
 198           CONTINUE
               IF (NXTLAB(LABFND,RTNLBL,LUPROP)) THEN
                  IF (LABFND(1:3) .EQ. 'SD ')
     &                 LLROP(INDPRP(LABFND)) = .TRUE.
                  GOTO 198
               END IF
               CALL GPCLOSE(LUPROP,'KEEP')
               GO TO 100
 42            CONTINUE
                  LLROP( INDPRP('XDIPLOC ')) = .TRUE.
                  LLROP( INDPRP('YDIPLOC ')) = .TRUE.
                  LLROP( INDPRP('ZDIPLOC ')) = .TRUE.
               GO TO 100
 43            CONTINUE
                  LLROP( INDPRP('XDIPLOC ')) = .TRUE.
               GO TO 100
 44            CONTINUE
                  LLROP( INDPRP('YDIPLOC ')) = .TRUE.
               GO TO 100
 45            CONTINUE
                  LLROP( INDPRP('ZDIPLOC ')) = .TRUE.
               GO TO 100
C
C .CHANNEL
 46            CONTINUE
               CHANNEL_CALC = .TRUE.
               CHANNEL_NORB = 0
               READ(LUCMD,*) (NUMBER_ORBS(I), I=1,NSYM)
               DO I=1,NSYM
                  CHANNEL_NORB = CHANNEL_NORB + NUMBER_ORBS(I)
               END DO
               IF (CHANNEL_NORB.GT.CHANNEL_MAX) CALL QUIT(
     &               'Too many orbitals for CHANNEL calculation: ' //
     &               'increase CHANNEL_MAX and recompile'
     &     )
               IOFF=0
               DO ISYM=1,NSYM
                  IF (NUMBER_ORBS(ISYM).GT.0) THEN
                     READ(LUCMD,*)
     &                    (CHANNEL_ORB(I+IOFF), I=1,NUMBER_ORBS(ISYM))
                     DO I=1,NUMBER_ORBS(ISYM)
                        CHANNEL_ORB(I+IOFF) =
     &                       CHANNEL_ORB(I+IOFF) + IORB(ISYM)
                     END DO
                     IOFF = IOFF + NUMBER_ORBS(ISYM)
                  ELSE
                     READ(LUCMD,*)
                  END IF
               END DO
               GO TO 100
C     .VIRTUAL
 47            CONTINUE
               CHANNEL_VCALC = .TRUE.
               READ (LUCMD,*) (CHANNEL_VIRT(I),I=1,NSYM)
               GO TO 100

 48            continue
!                do nothing
                 go to 100

            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') '@PROMPT "',WORD,
     *              '" NOT RECOGNIZED IN LRINP.'
               CALL QUIT(' ILLEGAL PROMPT IN LRINP ')
            END IF
      END IF
C     ---------------------------------------------
  300 CONTINUE
      IF (THR_REDFAC .GT. 0.0D0) THEN
         ICHANG = ICHANG + 1
         WRITE (LUPRI,'(3A,1P,D10.2)') '@ INFO ',WORD1,
     &   ' thresholds multiplied with general factor',THR_REDFAC
         THCLR = THCLR*THR_REDFAC
      END IF
C
      NLRTOT = 0
      IF (ICHANG .GT. 0) THEN
         DO 500 I = 1,NPRLBL
            IF (LLROP(I)) NLRTOT = NLRTOT + 1
 500     CONTINUE
         IF (NLRTOT .EQ. 0) THEN
            NINFO = NINFO + 1
            WRITE (LUPRI,'(/3A)') '@INFO: ',WORD1,
     *      ' input ignored because no operators requested.'
         END IF
      END IF
      IF (NLRTOT .GT. 0) THEN
         CALL HEADER('Linear Response calculation',0)
         IF (FLAG(16) .AND. .NOT.INERSI) THEN
            WRITE(LUPRI,'(A,L1,/)')
     *      ' Equilibrium solvent model requested            : SOLVNT ='
     *      ,FLAG(16)
            WRITE(LUPRI,'(A,F8.4)')
     *      ' Dielectric constant                            : EPSOL  ='
     *      ,EPSOL
         END IF
         IF (FLAG(16) .AND. INERSI) THEN
            WRITE(LUPRI,'(A,L1,/)')
     *      ' Non-equilibrium solvent model requested        : INERSI ='
     *      ,INERSI
            WRITE(LUPRI,'(A,F8.4)')
     *      ' Static dielectric constant                     : EPSTAT ='
     *      ,EPSTAT
            WRITE(LUPRI,'(A,F8.4)')
     *      ' Optical dielectric constant                    : EPSOL  ='
     *      ,EPSOL
         END IF
         IF (PCM .AND. .NOT.NEQRSP) THEN
            WRITE(LUPRI,'(A,L1,/)')
     *      ' Equilibrium PCM solvent model requested        : SOLVNT ='
     *      ,PCM
            WRITE(LUPRI,'(A,F8.4)')
     *      ' Dielectric constant                            : EPSOL  ='
     *      ,EPS
         END IF
         IF (PCM .AND. NEQRSP) THEN
            WRITE(LUPRI,'(A,L1,/)')
     *      ' Non-equilibrium PCM solvent model requested    : INERSI ='
     *      ,PCM
            WRITE(LUPRI,'(A,F8.4)')
     *      ' Static dielectric constant                     : EPSTAT ='
     *      ,EPS
            WRITE(LUPRI,'(A,F8.4)')
     *      ' Optical dielectric constant                    : EPSOL  ='
     *      ,EPSINF
         END IF

         IF (QM3 .OR. QMMM) THEN
           WRITE(LUPRI,'(A)') ' QM/MM response calculation'
         ENDIF

         IF (QMNPMM) THEN
             WRITE(LUPRI,'(A)') ' QM/NP/MM response calculation'
         END IF

         IF (USE_PELIB()) THEN
            WRITE(LUPRI,'(A)') ' Environment is modeled using'//
     &                  ' the polarizable embedding scheme (PE library)'
         END IF

         IF (CHANNEL_CALC) THEN
            WRITE(LUPRI,'(/A,10(/,10I5))')
     &         ' Restricted channel RPA,   excitations from orbs:',
     &           (CHANNEL_ORB(I),I=1,CHANNEL_NORB)
         END IF
         IF (CHANNEL_VCALC) THEN
            DO I=1,NSYM
              WRITE(LUPRI,'(/A,1I5,/A,1I5)') 'Including only the first',
     &              CHANNEL_VIRT(I),' virtual orbitals in symmetry',I
         ENDDO
      ENDIF
         WRITE(LUPRI,'(A,I4)')
     *      ' Print level                                    : IPRLR  ='
     *       ,IPRLR
         WRITE(LUPRI,'(A,I4)')
     *      ' Maximum number of iterations                   : MAXITL ='
     *       ,MAXITL
         WRITE(LUPRI,'(A,1P,D10.3)')
     *      ' Threshold for relative convergence             : THCLR  ='
     *      ,THCLR
         WRITE(LUPRI,'(A,I4)')
     *      ' Maximum iterations in optimal orbital algorithm: MAXITO ='
     *      ,MAXITO
         WRITE(LUPRI,'(/I3,A,(1P,5D14.6))')
     *      NFREQ,' B-frequencies',
     *      (FREQ(I),I=1,NFREQ)
         IF (RESTLR) WRITE(LUPRI,'(/A,L1)')
     *      ' USE TRIAL VECTORS AVAILABLE FOR RESTART. RESTLR =',RESTLR
         IF (ABOCHK) THEN
            WRITE(LUPRI,'(/2A,L1)')
     *      ' test option : SET UP ORBITAL PART OF E(2) AND S(2).',
     *      ' ABOCHK',ABOCHK
            IF ((ISTOCK.NE.1).OR.(MAXOCK.NE.6))
     *         WRITE(LUPRI,'(/A,I5,A,I5)')
     *         ' SUBSPACE SELECTED FROM ROW',ISTOCK,' TO',MAXOCK
         END IF
         IF (SOPRSY) WRITE(LUPRI,'(/A/A/A,L1)')
     *      ' TEST OPTION: BOTH UPPER AND LOWER BLOCK OF',
     *      ' POLARISABILITIES CALCULATED TO CHECK QUADRATIC ACCURACY',
     *      ' SOPRSY =',SOPRSY
      END IF
C
C *** END OF LRINP
C
      RETURN
      END
C  /* Deck c6inp */
      SUBROUTINE C6INP(WORD)
C
#include "implicit.h"
C
#include "priunit.h"
#include "gnrinf.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "rspprp.h"
#include "infc6.h"
#include "infpri.h"
#include "inftap.h"
C
      LOGICAL NEWDEF
      PARAMETER ( NTABLE = 41 )
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER*8 LABEL
C
      DATA TABLE /'.DIPLEN', '.DIPLNX', '.DIPLNY', '.DIPLNZ',
     *            '.DIPVEL', '.DIPVLX', '.DIPVLY', '.DIPVLZ',
     *            '.DIPMAG', '.DIPMGX', '.DIPMGY', '.DIPMGZ',
     *            '.QUADMO', '.QUADXX', '.QUADXY', '.QUADXZ',
     *            '.QUADYY', '.QUADYZ', '.QUADZZ', '.PROPRT',
     *            '.C6SPH ', '.C8SPH ', '.MAX IT', '.THCC6 ',
     *            '.MAXITO', '.PRINT ', '.MAXMOM', '.XXXXXX',
     *            '.FREQUE', '.XXXXXX', '.GSLEGN', '.NGRID ',
     *            '.XXXXXX', '.C10SPH', '.C6ATM ', '.C8ATM ',
     *            '.C10ATM', '.C6LMO ', '.C8LMO ', '.C10LMO',
     *            '.RESTC6'/
C
C
C READ IN  INPUT
C
      NEWDEF = (WORD .EQ. '*C6    ')
      ICHANG = 0
      IF (NEWDEF) THEN
      WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') GO TO 100
            IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,
     *                      18,19,20,21,22,23,24,25,26,27,28,29,30,31,
     *                      32,33,34,35,36,37,38,39,40,41), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') '@KEYWORD "',WORD,
     *            '" NOT RECOGNIZED IN C6INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT(' ILLEGAL KEYWORD IN C6INP ')
    1          CONTINUE
                  LC6OP( INDPRP('XDIPLEN ')) = .TRUE.
                  LC6OP( INDPRP('YDIPLEN ')) = .TRUE.
                  LC6OP( INDPRP('ZDIPLEN ')) = .TRUE.
               GO TO 100
    2          CONTINUE
                  LC6OP( INDPRP('XDIPLEN ')) = .TRUE.
               GO TO 100
    3          CONTINUE
                  LC6OP( INDPRP('YDIPLEN ')) = .TRUE.
               GO TO 100
    4          CONTINUE
                  LC6OP( INDPRP('ZDIPLEN ')) = .TRUE.
               GO TO 100
    5          CONTINUE
                  LC6OP( INDPRP('XDIPVEL ')) = .TRUE.
                  LC6OP( INDPRP('YDIPVEL ')) = .TRUE.
                  LC6OP( INDPRP('ZDIPVEL ')) = .TRUE.
               GO TO 100
    6          CONTINUE
                  LC6OP( INDPRP('XDIPVEL ')) = .TRUE.
               GO TO 100
    7          CONTINUE
                  LC6OP( INDPRP('YDIPVEL ')) = .TRUE.
               GO TO 100
    8          CONTINUE
                  LC6OP( INDPRP('ZDIPVEL ')) = .TRUE.
               GO TO 100
    9          CONTINUE
                  LC6OP( INDPRP('XANGMOM ')) = .TRUE.
                  LC6OP( INDPRP('YANGMOM ')) = .TRUE.
                  LC6OP( INDPRP('ZANGMOM ')) = .TRUE.
               GO TO 100
   10          CONTINUE
                  LC6OP( INDPRP('XANGMOM ')) = .TRUE.
               GO TO 100
   11          CONTINUE
                  LC6OP( INDPRP('YANGMOM ')) = .TRUE.
               GO TO 100
   12          CONTINUE
                  LC6OP( INDPRP('ZANGMOM ')) = .TRUE.
               GO TO 100
   13          CONTINUE
                  LC6OP( INDPRP('XXQUADRU')) = .TRUE.
                  LC6OP( INDPRP('XYQUADRU')) = .TRUE.
                  LC6OP( INDPRP('XZQUADRU')) = .TRUE.
                  LC6OP( INDPRP('YYQUADRU')) = .TRUE.
                  LC6OP( INDPRP('YZQUADRU')) = .TRUE.
               GO TO 100
   14          CONTINUE
                  LC6OP( INDPRP('XXQUADRU')) = .TRUE.
               GO TO 100
   15          CONTINUE
                  LC6OP( INDPRP('XYQUADRU')) = .TRUE.
               GO TO 100
   16          CONTINUE
                  LC6OP( INDPRP('XZQUADRU')) = .TRUE.
               GO TO 100
   17          CONTINUE
                  LC6OP( INDPRP('YYQUADRU')) = .TRUE.
               GO TO 100
   18          CONTINUE
                  LC6OP( INDPRP('YZQUADRU')) = .TRUE.
               GO TO 100
   19          CONTINUE
                  LC6OP( INDPRP('ZZQUADRU')) = .TRUE.
               GO TO 100
   20          CONTINUE
                  READ(LUCMD,'(A)')LABEL
                  LC6OP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
C     .C6SPH
   21          CONTINUE
                  C6IFC = .TRUE.
                  LC6OP( INDPRP('SM01-01 ')) = .TRUE.
                  LC6OP( INDPRP('SM01+00 ')) = .TRUE.
                  LC6OP( INDPRP('SM01+01 ')) = .TRUE.
               GO TO 100
C     .C8SPH
   22          CONTINUE
                  C6IFC = .TRUE.
                  LC6OP( INDPRP('SM01-01 ')) = .TRUE.
                  LC6OP( INDPRP('SM01+00 ')) = .TRUE.
                  LC6OP( INDPRP('SM01+01 ')) = .TRUE.
                  LC6OP( INDPRP('SM02-02 ')) = .TRUE.
                  LC6OP( INDPRP('SM02-01 ')) = .TRUE.
                  LC6OP( INDPRP('SM02+00 ')) = .TRUE.
                  LC6OP( INDPRP('SM02+01 ')) = .TRUE.
                  LC6OP( INDPRP('SM02+02 ')) = .TRUE.
                  LC6OP( INDPRP('SM03-03 ')) = .TRUE.
                  LC6OP( INDPRP('SM03-02 ')) = .TRUE.
                  LC6OP( INDPRP('SM03-01 ')) = .TRUE.
                  LC6OP( INDPRP('SM03+00 ')) = .TRUE.
                  LC6OP( INDPRP('SM03+01 ')) = .TRUE.
                  LC6OP( INDPRP('SM03+02 ')) = .TRUE.
                  LC6OP( INDPRP('SM03+03 ')) = .TRUE.
               GO TO 100
   23          CONTINUE
                  READ (LUCMD,*) MAXITC
               GO TO 100
   24          CONTINUE
                  READ (LUCMD,*) THCC6
               GO TO 100
   25          CONTINUE
                  READ (LUCMD,*) MAXITO
               GO TO 100
   26          CONTINUE
                  READ (LUCMD,*) IPRC6
               GO TO 100
   27          CONTINUE
                  READ (LUCMD,*) MAXMOM
               GO TO 100
   28          CONTINUE
               GO TO 100
   29          CONTINUE
                  READ (LUCMD,*) NCFREQ
                  IF (NCFREQ.LE.MCFREQ) THEN
                     READ (LUCMD,*) (CFREQ(J),J=1,NCFREQ)
                  ELSE
                     WRITE (LUPRI,'(3(/,A,I5),/)')
     *               '@NUMBER OF FREQUENCIES SPECIFIED    :',NCFREQ,
     *               '@IS GREATER THAN THE ALLOWED NUMBER :',MCFREQ,
     *               '@THE NUMBER IS RESET TO THE MAXIMUM :',MCFREQ
                     READ (LUCMD,*) (CFREQ(J),J=1,MCFREQ),
     *                              (FFFF,J=MCFREQ+1,NCFREQ)
                     NCFREQ = MCFREQ
                  END IF
               GO TO 100
   30          CONTINUE
               GO TO 100
   31          CONTINUE
                  GSLEGN  =.TRUE.
               GO TO 100
   32          CONTINUE
                  READ (LUCMD,*) NGRID
               GO TO 100
C     .C6CALC
   33          CONTINUE
                  LC6OP( INDPRP('XDIPLEN ')) = .TRUE.
                  LC6OP( INDPRP('YDIPLEN ')) = .TRUE.
                  LC6OP( INDPRP('ZDIPLEN ')) = .TRUE.
               GO TO 100
C     .C10SPH
   34          CONTINUE
                  C6IFC = .TRUE.
                  LC6OP( INDPRP('SM01-01 ')) = .TRUE.
                  LC6OP( INDPRP('SM01+00 ')) = .TRUE.
                  LC6OP( INDPRP('SM01+01 ')) = .TRUE.
                  LC6OP( INDPRP('SM02-02 ')) = .TRUE.
                  LC6OP( INDPRP('SM02-01 ')) = .TRUE.
                  LC6OP( INDPRP('SM02+00 ')) = .TRUE.
                  LC6OP( INDPRP('SM02+01 ')) = .TRUE.
                  LC6OP( INDPRP('SM02+02 ')) = .TRUE.
                  LC6OP( INDPRP('SM03-03 ')) = .TRUE.
                  LC6OP( INDPRP('SM03-02 ')) = .TRUE.
                  LC6OP( INDPRP('SM03-01 ')) = .TRUE.
                  LC6OP( INDPRP('SM03+00 ')) = .TRUE.
                  LC6OP( INDPRP('SM03+01 ')) = .TRUE.
                  LC6OP( INDPRP('SM03+02 ')) = .TRUE.
                  LC6OP( INDPRP('SM03+03 ')) = .TRUE.
                  LC6OP( INDPRP('SM04-04 ')) = .TRUE.
                  LC6OP( INDPRP('SM04-03 ')) = .TRUE.
                  LC6OP( INDPRP('SM04-02 ')) = .TRUE.
                  LC6OP( INDPRP('SM04-01 ')) = .TRUE.
                  LC6OP( INDPRP('SM04+00 ')) = .TRUE.
                  LC6OP( INDPRP('SM04+01 ')) = .TRUE.
                  LC6OP( INDPRP('SM04+02 ')) = .TRUE.
                  LC6OP( INDPRP('SM04+03 ')) = .TRUE.
                  LC6OP( INDPRP('SM04+04 ')) = .TRUE.
                  LC6OP( INDPRP('SM05-05 ')) = .TRUE.
                  LC6OP( INDPRP('SM05-04 ')) = .TRUE.
                  LC6OP( INDPRP('SM05-03 ')) = .TRUE.
                  LC6OP( INDPRP('SM05-02 ')) = .TRUE.
                  LC6OP( INDPRP('SM05-01 ')) = .TRUE.
                  LC6OP( INDPRP('SM05+00 ')) = .TRUE.
                  LC6OP( INDPRP('SM05+01 ')) = .TRUE.
                  LC6OP( INDPRP('SM05+02 ')) = .TRUE.
                  LC6OP( INDPRP('SM05+03 ')) = .TRUE.
                  LC6OP( INDPRP('SM05+04 ')) = .TRUE.
                  LC6OP( INDPRP('SM05+05 ')) = .TRUE.
               GO TO 100
C     .C6ATM
   35          CONTINUE
                  C6IFC = .TRUE.
                  LC6OP( INDPRP('SM01+00 ')) = .TRUE.
               GO TO 100
C     .C8ATM
   36          CONTINUE
                  C6IFC = .TRUE.
                  LC6OP( INDPRP('SM01+00 ')) = .TRUE.
                  LC6OP( INDPRP('SM02+00 ')) = .TRUE.
                  LC6OP( INDPRP('SM03+00 ')) = .TRUE.
               GO TO 100
C     .C10ATM
   37          CONTINUE
                  C6IFC = .TRUE.
                  LC6OP( INDPRP('SM01+00 ')) = .TRUE.
                  LC6OP( INDPRP('SM02+00 ')) = .TRUE.
                  LC6OP( INDPRP('SM03+00 ')) = .TRUE.
                  LC6OP( INDPRP('SM04+00 ')) = .TRUE.
                  LC6OP( INDPRP('SM05+00 ')) = .TRUE.
               GO TO 100
C     .C6LMO
   38          CONTINUE
                  C6IFC = .TRUE.
                  LC6OP( INDPRP('SM01+00 ')) = .TRUE.
                  LC6OP( INDPRP('SM01+01 ')) = .TRUE.
               GO TO 100
C     .C8LMO
   39          CONTINUE
                  C6IFC = .TRUE.
                  LC6OP( INDPRP('SM01+00 ')) = .TRUE.
                  LC6OP( INDPRP('SM01+01 ')) = .TRUE.
                  LC6OP( INDPRP('SM02+00 ')) = .TRUE.
                  LC6OP( INDPRP('SM02+01 ')) = .TRUE.
                  LC6OP( INDPRP('SM02+02 ')) = .TRUE.
                  LC6OP( INDPRP('SM03+00 ')) = .TRUE.
                  LC6OP( INDPRP('SM03+01 ')) = .TRUE.
                  LC6OP( INDPRP('SM03+02 ')) = .TRUE.
                  LC6OP( INDPRP('SM03+03 ')) = .TRUE.
               GO TO 100
C     .C10LMO
   40          CONTINUE
                  C6IFC = .TRUE.
                  LC6OP( INDPRP('SM01+00 ')) = .TRUE.
                  LC6OP( INDPRP('SM01+01 ')) = .TRUE.
                  LC6OP( INDPRP('SM02+00 ')) = .TRUE.
                  LC6OP( INDPRP('SM02+01 ')) = .TRUE.
                  LC6OP( INDPRP('SM02+02 ')) = .TRUE.
                  LC6OP( INDPRP('SM03+00 ')) = .TRUE.
                  LC6OP( INDPRP('SM03+01 ')) = .TRUE.
                  LC6OP( INDPRP('SM03+02 ')) = .TRUE.
                  LC6OP( INDPRP('SM03+03 ')) = .TRUE.
                  LC6OP( INDPRP('SM04+00 ')) = .TRUE.
                  LC6OP( INDPRP('SM04+01 ')) = .TRUE.
                  LC6OP( INDPRP('SM04+02 ')) = .TRUE.
                  LC6OP( INDPRP('SM04+03 ')) = .TRUE.
                  LC6OP( INDPRP('SM04+04 ')) = .TRUE.
                  LC6OP( INDPRP('SM05+00 ')) = .TRUE.
                  LC6OP( INDPRP('SM05+01 ')) = .TRUE.
                  LC6OP( INDPRP('SM05+02 ')) = .TRUE.
                  LC6OP( INDPRP('SM05+03 ')) = .TRUE.
                  LC6OP( INDPRP('SM05+04 ')) = .TRUE.
                  LC6OP( INDPRP('SM05+05 ')) = .TRUE.
               GO TO 100
   41          CONTINUE
                  RESTC6 = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') '@PROMPT "',WORD,
     *            '" NOT RECOGNIZED IN C6INP.'
               CALL QUIT(' ILLEGAL PROMPT IN C6INP ')
            END IF
      END IF
C
  300 CONTINUE
      IF (THR_REDFAC .GT. 0.0D0) THEN
         ICHANG = ICHANG + 1
         WRITE (LUPRI,'(3A,1P,D10.2)') '@ INFO ',WORD1,
     &   ' thresholds multiplied with general factor',THR_REDFAC
         THCC6 = THCC6*THR_REDFAC
      END IF
C
      NC6TOT = 0
      IF (ICHANG .GT. 0) THEN
         DO 500 I = 1,NPRLBL
            IF (LC6OP(I)) NC6TOT = NC6TOT + 1
 500     CONTINUE
         IF (NC6TOT .EQ. 0) WRITE (LUPRI,'(/A)')
     *      '@ *RSPC6 input ignored because no operators requested.'
      END IF
      IF (NC6TOT .GT. 0) THEN
         CALL HEADER('CHANGES OF DEFAULTS IN C6INP:',0)
         IF (MAXITC.NE.60) WRITE(LUPRI,'(/A,I6)')
     *      ' MAXIMUM NUMBER OF ITERATIONS. MAXITC =',MAXITC
         IF (THCC6.NE.1.0D-3) WRITE(LUPRI,'(/A,1P,D13.6)')
     *      ' THRESHOLD FOR CONVERGENCE. THCC6 =',THCC6
         IF (MAXITO.NE.5) WRITE(LUPRI,'(/A,I5)')
     *      ' MAXIMUM ITERATIONS IN OPTIMAL ORBITAL ALGORITHM. MAXITO ='
     *      ,MAXITO
         IF (MAXMOM.NE.6) WRITE(LUPRI,'(/A,I5)')
     *      ' MAXIMUM NUMBER OF MOMENTS CALCULATED. MAXMOM =',MAXMOM
         IF (RESTC6) WRITE(LUPRI,'(/A)')
     *      ' RESTART OF MOMENT ITERATIONS.'
         IF (GSLEGN) WRITE(LUPRI,'(/2A/A,L1)')
     *      ' Gauss Legendre GRID USED FOR INTEGRATION OF IMAGINARY'
     *      ,'POLARISABILITIES',' GSLEGN =',GSLEGN
         IF (NGRID.NE.10) WRITE(LUPRI,'(/A,I5)')
     *      ' NUMBER OF GRID POINTS. NGRID = ',NGRID
         IF (NCFREQ.NE.0) WRITE(LUPRI,'(/A,I5,A/,(1P,5D14.6))')
     *      ' POLARIZABILITIES CALCULATED AT',NCFREQ,' FREQUENCIES',
     *      (CFREQ(I),I=1,NCFREQ)
         IF (IPRC6.NE.2) WRITE(LUPRI,'(/A,I5)')
     *      ' PRINT LEVEL IN C6 ROUTINES. IPRC6 =',IPRC6
         IF (C6IFC) WRITE(LUPRI,'(/A)')
     &      ' Interface file for C6, C8, or C10 calculation written.'
      END IF
C
C *** END OF C6INP
C
      RETURN
      END
C  /* Deck rspset */
      SUBROUTINE RSPSET
C
C 1) INITIALIZE VARIABLES AND DIMENSIONS WHICH ARE COMMON TO ALL
C    THE RESPONSE CALCULATIONS AND WHICH ARE NOT DEFINED BY SIRSET
C    (INPUT TO THE SIRIUS PROGRAM)
C
C 2) INITIALIZE VARIABLES FOR INPUT
C
C    List of updates:
C    21-July-1992 : Hinne Hettema : Test RSPSUP added.
C
C    Used from common blocks:
C    INFINP : SUPSYM
C    INFORB : NSYM, NSSYM
C
#include "implicit.h"
C
#include "maxorb.h"
#include "priunit.h"
#include "infinp.h"
#include "inforb.h"
#include "infvar.h"
C
#include "infrsp.h"
#include "wrkrsp.h"
#include "inftap.h"
#include "rspprp.h"
#include "infpp.h"
#include "inflr.h"
#include "infc6.h"
#include "infsmo.h"
#include "infhyp.h"
#include "qrinf.h"
#include "infesr.h"
#include "trhso.h"
#include "infcr.h"
C
      LOGICAL FNDLAB
C
C ****** Used unit names: ******
C        (the file numbers are now assigned by GPOPEN)
C
C LUSIFC   FILE WITH CONVERGED MCSCF RESUSTS
C
C LURSP*   RESPONSE ITERATION FILES
C          LURSP1 SAVE NECESSARY INFORMATION FOR RESTART
C          LURSP2 SAVE LR SOLUTION AND RESIDUAL VECTORS
C          LURSP3 SAVE TRIAL VECTORS
C          LURSP4 SAVE INFORMATION FOR DIAGONAL E(2) AND S(2) AND
C                      PHPINFORMATION
C          LURSP5 SAVE E(2) LINEAR TRANSFORMED TRIAL VECTORS
C                      EQUATIONS FOR QUADRATIC RESPONSE
C          LURD   unit for openining ABACUS file with RD vectors.
C          LUAHSO spin-orbit atomic two-electron integrals
C          LUMHSO spin-orbit molecular two-electron integrals
C          LUCRVE response vector information
C          LUCRV1 one-index linear response vectors
C          LUCRV2 two-index linear response vectors
C          LUXYVE XY-vectors
C
      IF (LUSIFC .LE. 0) CALL GPOPEN(LUSIFC,FNSIFC,'OLD',' ',
     &     'UNFORMATTED',IDUMMY,.FALSE.)
      REWIND LUSIFC
      IF (FNDLAB('CIRESPON',LUSIFC)) THEN
         RSPCI = .TRUE.
         WRITE(LUPRI,'(/A/A)')
     *     ' This is a configuration interaction response calculation.',
     *     '   Note: this is equivalent to a CI sum-over-states'//
     *     ' calculation of properties.'
      ELSE
         RSPCI = .FALSE.
      END IF
C
      NPRLBL = 0
      DO 100 I = 1,MAXLBL
         LPPOP(I) = .FALSE.
         LLROP(I) = .FALSE.
         LC6OP(I) = .FALSE.
         AQROP(I) = .FALSE.
         BQROP(I) = .FALSE.
         CQROP(I) = .FALSE.
         ASMOP(I) = .FALSE.
         BSMOP(I) = .FALSE.
         ESROPS(I) = .FALSE.
         ESROPT(I) = .FALSE.
         ACROP(I) = .FALSE.
         BCROP(I) = .FALSE.
         CCROP(I) = .FALSE.
         DCROP(I) = .FALSE.
 100  CONTINUE
C
C  *** Set Supersymetry value if necessary
C      Check if supersymmetry is identical to "D2h" symmetry
C
      RSPSUP = SUPSYM
      IF (SUPSYM .AND. NSSYM .LE. NSYM .AND. MXDGSS .EQ. 1) THEN
         DO 220 ISYM = 1,NSYM
            JSYM = 0
            DO 210 ISSYM = 1,NSSYM
               IF (NINFSS(ISSYM,3) .EQ. ISYM) JSYM = JSYM + 1
  210       CONTINUE
            IF (JSYM .GT. 1) GO TO 300
  220    CONTINUE
         RSPSUP = .FALSE.
      END IF
  300 CONTINUE
C
C *** END OF RSPSET
C
      RETURN
      END
C  /* Deck indprp */
      INTEGER FUNCTION INDPRP(NEWLBL)
C
#include "priunit.h"
#include "rspprp.h"
C
      CHARACTER*8 NEWLBL
      DO 100 I = 1,NPRLBL
         IF ( NEWLBL .EQ. PRPLBL(I) ) THEN
            INDPRP = I
            RETURN
         END IF
 100  CONTINUE
      NPRLBL = NPRLBL + 1
      IF (NPRLBL.GT.MAXLBL) THEN
         WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *   '@ NUMBER OF SPECIFIED PROPERTIES EXCEEDS THE MAXIMUM ALLOWED',
     *   '@ MAXPRP =',MAXLBL,' NPRLBL= ',NPRLBL
         CALL QUIT(' INDPRP: TOO MANY DIFFERENT PROPERTIES SPECIFIED')
      END IF
      PRPLBL(NPRLBL) = NEWLBL
      INDPRP = NPRLBL
      RETURN
      END
C  /* Deck rspsym */
      SUBROUTINE RSPSYM(WRK,LWRK)
C
C DETERMINE HOW MANY OPERATORS OF EACH SYMMETRY THAT HAS TO BE
C CALCULATED IN PP AND LR CALCULATION
C
#include "implicit.h"
#include "iratdef.h"
#include "dummy.h"
C
      DIMENSION  WRK(*)
C
      PARAMETER ( D1 = 1.0D0 , D2 = 2.0D0 , D8=8.0D0 , D1000=1000.D0)
      PARAMETER ( D1INF = 0.99999D0 ,THREQL = 1.0D-8 ,CKMXPR = 1.0D-12)
C
#include "maxorb.h"
#include "maxash.h"
#include "priunit.h"
#include "infrsp.h"
#include "infind.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "rspprp.h"
#include "infpp.h"
#include "inflr.h"
#include "infc6.h"
#include "infhyp.h"
#include "indqr.h"
#include "infsmo.h"
#include "infs0.h"
#include "infesr.h"
#include "infcr.h"
#include "indcr.h"
#include "inftmo.h"
#include "inftpa.h"
#include "orgcom.h"
#include "gtensor.h"
C
      LOGICAL FNDLAB,FNDLB2,TRIMAT,TOFILE
C
      CHARACTER*8 TABOVL, TABLEN(3) ,TABVEL(3) ,RTNLBL(2)
      CHARACTER*8 LABEL
      CHARACTER*7 PRPWRD
      CHARACTER*8 LABINT(3)
C
C     Note that the length of LABINT may change if we allow
C     integrals with more then three components to be calculated by
C     HERMIT "on the fly", K.Ruud-Nov-96
C
      DATA TABLEN/'XDIPLEN ', 'YDIPLEN ', 'ZDIPLEN '/
      DATA TABVEL/'XDIPVEL ', 'YDIPVEL ', 'ZDIPVEL '/
      DATA TABOVL/'OVERLAP '/
C
      CALL QENTER('RSPSYM')
C
C GIVE ARBITRARY EXCITATION ENERGIES TO OBTAIN THE RIGHT
C NUMBER OF LINEAR RESPONSE EQUATIONS THAT HAVE TO BE SOLVED
C IN CALCULATION OF TRANSITION MOMENTS BETWEEN EXCITED STATES
C AND SECOND ORDER TRANSITION MOMENTS
C THESE VALUES WILL LATER BE OVERWRITTEN WITH THE CALCULATED
C VALUES
C
      DO 40 K = 1,2
         DO 40 J = 1,MXEXQR
            DO 40 I = 1,8
               EXCITA(I,J,K) = D1000 + I + J
  40  CONTINUE
C
C Repeat for third moment calculation and two-photon absorption
C
      DO 50 J = 1,MXEXCR
      DO 50 I = 1,8
         EXCIT2(I,J) = D1000 + I + J
  50  CONTINUE
C
C DEFINE DEFAULT VALUES FOR NUMBER OF SIMULTANEOUS VECTORS AND FOR
C NUMBER OF START VECTORS, IF THEY ARE NOT DEFINED IN INPUT.
C INITIALIZE PROPERTY COUNTS TO ZERO.
C
      NLRLBL = 0
      NTRLBL = 0
      NEXLBL = 0
      DO 100 ISYM = 1,NSYM
         NPPSIM(ISYM) = MAX(NPPSIM(ISYM),NPPCNV(ISYM))
         NPPSTV(ISYM) = MAX(NPPSTV(ISYM),NPPSIM(ISYM))
         NGPPP(ISYM)  = 0
         NGPLR(ISYM)  = 0
         NGPC6(ISYM)  = 0
         NGPS0(ISYM)  = 0
         NAQROP(ISYM) = 0
         NBQROP(ISYM) = 0
         NCQROP(ISYM) = 0
         NASMOP(ISYM) = 0
         NBSMOP(ISYM) = 0
         NESRS(ISYM)  = 0
         NESRT(ISYM)  = 0
         NACROP(ISYM) = 0
         NBCROP(ISYM) = 0
         NCCROP(ISYM) = 0
         NDCROP(ISYM) = 0
  100 CONTINUE
C
C DETERMINE THE TOTAL NUMBER OF PERT. OPERATORS OF EACH SYMMETRY
C THAT HAS TO BE CALCULATED IN PP AND IN LR CALCULATION AND
C LABEL EACH OPERATOR
C
      IF (IPRRSP.GE.40) THEN
         WRITE(LUPRI,'(/A)')' **** Output from RSPSYM ****'
         WRITE(LUPRI,'(/A/)')
     &      ' I, LPPOP(I), PRPLBL(I)  (for excitations)'
         DO 224 I = 1,NPRLBL
            WRITE(LUPRI,'(I5,L10,2X,A8)')I,LPPOP(I),PRPLBL(I)
 224     CONTINUE
         WRITE(LUPRI,'(/A/)')
     &      ' I, LLROP(I), PRPLBL(I)  (for linear response)'
         DO 324 I = 1,NPRLBL
            WRITE(LUPRI,'(I5,L10,2X,A8)')I,LLROP(I),PRPLBL(I)
 324     CONTINUE
         WRITE(LUPRI,'(/A/)')' I, LC6OP(I), PRPLBL(I)  (for C6/C8)'
         DO 424 I = 1,NPRLBL
            WRITE(LUPRI,'(I5,L10,2X,A8)')I,LC6OP(I),PRPLBL(I)
 424     CONTINUE
         WRITE(LUPRI,'(/A/)') ' I, AQROP(I), BQROP(I), CQROP(I),'//
     &      ' PRPLBL(I)  (for quadratic response)'
         DO 425 I = 1,NPRLBL
            WRITE(LUPRI,'(I5,3L10,2X,A8)')
     &         I,AQROP(I),BQROP(I),CQROP(I),PRPLBL(I)
 425     CONTINUE
         WRITE(LUPRI,'(/A/)')
     &      ' I, ASMOP(I), BSMOP(I), PRPLBL(I)  (for second mom.)'
         DO 428 I = 1,NPRLBL
            WRITE(LUPRI,'(I5,2L10,2X,A8)')I,ASMOP(I),BSMOP(I),PRPLBL(I)
 428     CONTINUE
C
         WRITE(LUPRI,'(/A/)')
     &      ' I, ESROPS(I), ESROPT(I), PRPLBL(I) (for ESR prop.)'
         DO 431 I = 1,NPRLBL
           WRITE(LUPRI,'(I5,2L10,2X,A8)')I,ESROPS(I),ESROPT(I),PRPLBL(I)
 431     CONTINUE
C
         WRITE(LUPRI,'(/A/)') ' I, ACROP(I), BCROP(I),'//
     &      ' CCROP(I), DCROP(I), PRPLBL(I)  (for cubic response)'
         DO 433 I = 1,NPRLBL
            WRITE(LUPRI,'(I5,4L10,2X,A8)')
     &         I,ACROP(I),BCROP(I),CCROP(I),DCROP(I),PRPLBL(I)
 433     CONTINUE
         WRITE(LUPRI,'(/A/)')
     &   ' I, ATMOP(I), BTMOP(I), CTMOP(I), PRPLBL(I)  (for third mom.)'
         DO 440 I = 1,NPRLBL
            WRITE(LUPRI,'(I5,3L10,2X,A8)')
     &         I,ATMOP(I),BTMOP(I),CTMOP(I),PRPLBL(I)
 440     CONTINUE
         WRITE(LUPRI,'(/A/)')
     &   ' I, ATPOP(I), BTPOP(I), PRPLBL(I) (for cubic double res.)'
         DO 443 I = 1,NPRLBL
            WRITE(LUPRI,'(I5,2L10,2X,A8)')I,ATPOP(I),BTPOP(I),PRPLBL(I)
 443     CONTINUE
      END IF
      IF (.NOT.NOS0MX) THEN
C        ... make sure DIPLEN and DIPVEL labels for S0MIX are in PRPLBL.
C            They will be calculated automatically with GET1INT call below
C            if not already available on AOPROPER.  /hjaaj June 2004
         DO IPRP = 1,3
            J =  INDPRP(TABLEN(IPRP))
            J =  INDPRP(TABVEL(IPRP))
         END DO
      END IF
C
      IF (IPRRSP.GE.30) THEN
         WRITE(LUPRI,'(/A,I5,A//A)')
     & ' NPRLBL=',NPRLBL,' labels in property calculations in RESPONS:',
     & '  I      PRPLBL(I), I=1,NPRLBL '
         DO I = 1,NPRLBL
            WRITE(LUPRI,'(I5,5X,A)')I,PRPLBL(I)
         END DO
      END IF
      CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      DO 250 IPRP = 1,NPRLBL
         JTURN = 0
 248     CONTINUE
         JTURN = JTURN + 1
         REWIND (LUPROP)
         LABEL = PRPLBL(IPRP)
         IF (LABEL(3:8).EQ.'SPNORB') LABEL(2:2) = '1'
         IF (FNDLB2(LABEL,RTNLBL,LUPROP)) THEN
            IF (RTNLBL(2) .EQ. 'SQUARE  ') THEN
               CALL READT(LUPROP,N2BASX,WRK)
               IND   = IDAMAX(N2BASX,WRK,1)
               IF (ABS(WRK(IND)).GT.CKMXPR) THEN
                  IIND = (IND-1)/NBAST + 1
                  JIND = IND - (IIND-1)*NBAST
                  KSYMPT = MULD2H( ISAO(IIND) , ISAO(JIND) )
               ELSE
                  KSYMPT = 0
                  NWARN  = NWARN + 1
                  WRITE(LUPRI,'(/A,/3A,1P,D15.7,/A)')'@ ** WARNING **',
     *               '@ ELEMENTS ON PROPERTY FILE  WITH LABEL: ',
     *               LABEL,       ' ARE SMALLER THAN :',CKMXPR,
     *               '@ PROPERTY CALCULATION NOT CARRIED OUT '
               END IF
            ELSE
               CALL READT(LUPROP,NNBASX,WRK)
               IND   = IDAMAX(NNBASX,WRK,1)
               IF (ABS(WRK(IND)).GT.CKMXPR) THEN
                  IIND  = INT( (SQRT(D1+D8*IND) - D1INF)/ D2 )
                  IIMAX = IIND*(IIND+1)/2
                  IF ( IIMAX .EQ. IND ) THEN
                     JIND = IIND
                  ELSE
                     IIND = IIND + 1
                     JIND = IND  - IIMAX
                  END IF
                  KSYMPT = MULD2H( ISAO(IIND) , ISAO(JIND) )
               ELSE
                  KSYMPT = 0
                  NWARN  = NWARN + 1
                  WRITE(LUPRI,'(/A,/3A,1P,D15.7,/A)')'@ ** WARNING **',
     *               '@ ELEMENTS ON PROPERTY FILE  WITH LABEL: ',
     *               LABEL,' ARE SMALLER THAN :',CKMXPR,
     *               '@ PROPERTY CALCULATION NOT CARRIED OUT '
               END IF
            END IF
         ELSE
            IF (JTURN .GT. 1) THEN
C           ... capture input errors as asking for ".DIPLENX" or "UDIPLEN"
C               which otherwise would cause an infinite loop here /hjaaj June 2004
               WRITE(LUPRI,'(/3A)') '@ RSPSYM ERROR: PROPERTY LABEL "',
     &            LABEL,'" NOT FOUND ON AOPROPER'
               CALL QUIT('RSPSYM ERROR: PROPERTY LABEL NOT FOUND ON '//
     &                   'AOPROPER')
            ELSE IF (LABEL(1:5) .EQ. 'dh/dB') THEN
               TRIMAT = .TRUE.
               PRPWRD = 'MAGMOM '
            ELSE IF (LABEL(1:7) .EQ. 'd|S>/dB') THEN
               TRIMAT = .FALSE.
               PRPWRD = 'S1MAGR '
            ELSE IF (LABEL(2:7) .EQ. 'ANGMOM') THEN
               TRIMAT = .TRUE.
               PRPWRD = 'ANGMOM '
            ELSE IF (LABEL(2:7) .EQ. 'DIPLEN') THEN
               TRIMAT = .TRUE.
               PRPWRD = 'DIPLEN '
            ELSE IF (LABEL(2:7) .EQ. 'DIPVEL') THEN
               TRIMAT = .TRUE.
               PRPWRD = 'DIPVEL '
            ELSE IF (LABEL(3:8) .EQ. 'SPNORB') THEN
               TRIMAT = .TRUE.
               PRPWRD = 'SPNORB'
            ELSE IF (LABEL(3:8) .EQ. 'MNF-SO') THEN
               TRIMAT = .TRUE.
               PRPWRD = 'MNF-SO'
            ELSE IF (LABEL(3:8) .EQ. 'SPNSCA') THEN
               TRIMAT = .TRUE.
               PRPWRD = 'SOSCALE'
            ELSE IF (LABEL(2:7) .EQ. 'ANGECC') THEN
               TRIMAT = .TRUE.
               PRPWRD = 'ANGECC'
            ELSE IF (LABEL(2:7) .EQ. 'DIPLOC') THEN
               TRIMAT = .TRUE.
               PRPWRD = 'DIPLOC '
            ELSE
               WRITE(LUPRI,'(/3A)') '@ RSPSYM ERROR: PROPERTY LABEL "',
     &            LABEL,'" NOT FOUND ON AOPROPER'
               CALL QUIT('RSPSYM ERROR: PROPERTY LABEL NOT FOUND ON '//
     &                   'AOPROPER')
            END IF
C
C     Ask HERMIT to compute the integrals. As we will probably need more
C     components later, we write to disk and read them in later when
C     necessary, K.Ruud,June 10 1996
C
            CALL GPCLOSE(LUPROP,'KEEP')
            NPATOM = 0
            NCOMP  = 0
            TOFILE = .TRUE.
            KINTAD = 1
            KINTRP = KINTAD + 4/IRAT
            KLAST  = KINTRP + 4/IRAT
            LWRK1  = LWRK   - KLAST
            IF (LWRK1.LT.0) CALL ERRWRK('RSPSYM',KLAST,LWRK)
C
C     If you are doing a PCM calculation with symmetry, you may risk
C     that a symmetry-breaking tesserae geometry is left over in the
C     DIPORG array. Reset it to a symmetry-conserving number
C
            DO I = 1, 3
               DIPORG(I) = CMXYZ(I)
            END DO
            CALL GET1IN(DUMMY,PRPWRD,NCOMP,WRK(KLAST),LWRK1,LABINT,
     &                  WRK(KINTRP),WRK(KINTAD),IDUMMY,TOFILE,NPATOM,
     &                  TRIMAT,DUMMY,.FALSE.,DUMMY,IPRRSP - 120)
            CALL GPOPEN(LUPROP,'AOPROPER','OLD',' ','UNFORMATTED',
     &                  IDUMMY,.FALSE.)
C
C     Finished making the integrals, then we read them from disk.
C
            GOTO 248
         END IF
         IF ( (KSYMPT .GT. 0) .AND. (KSYMPT .LE. NSYM) ) THEN
            IF (LPPOP(IPRP)) THEN
               NGPPP(KSYMPT) = NGPPP(KSYMPT) + 1
               LBLPP(KSYMPT,NGPPP(KSYMPT)) = PRPLBL(IPRP)
            END IF
            IF (LLROP(IPRP)) THEN
               NGPLR(KSYMPT) = NGPLR(KSYMPT) + 1
               LBLLR(KSYMPT,NGPLR(KSYMPT)) = PRPLBL(IPRP)
            END IF
            IF (LC6OP(IPRP)) THEN
               NGPC6(KSYMPT) = NGPC6(KSYMPT) + 1
               LBLC6(KSYMPT,NGPC6(KSYMPT)) = PRPLBL(IPRP)
            END IF
            IF (AQROP(IPRP)) THEN
               NAQROP(KSYMPT) = NAQROP(KSYMPT) + 1
               AQRLB(KSYMPT,NAQROP(KSYMPT)) = PRPLBL(IPRP)
            END IF
            IF (BQROP(IPRP)) THEN
               NBQROP(KSYMPT) = NBQROP(KSYMPT) + 1
               BQRLB(KSYMPT,NBQROP(KSYMPT)) = PRPLBL(IPRP)
            END IF
            IF (CQROP(IPRP)) THEN
               NCQROP(KSYMPT) = NCQROP(KSYMPT) + 1
               CQRLB(KSYMPT,NCQROP(KSYMPT)) = PRPLBL(IPRP)
            END IF
            IF (ASMOP(IPRP)) THEN
               NASMOP(KSYMPT) = NASMOP(KSYMPT) + 1
               ASMLB(KSYMPT,NASMOP(KSYMPT)) = PRPLBL(IPRP)
            END IF
            IF (BSMOP(IPRP)) THEN
               NBSMOP(KSYMPT) = NBSMOP(KSYMPT) + 1
               BSMLB(KSYMPT,NBSMOP(KSYMPT)) = PRPLBL(IPRP)
            END IF
            IF (ESROPS(IPRP)) THEN
               NESRS(KSYMPT) = NESRS(KSYMPT) + 1
               LBESRS(KSYMPT,NESRS(KSYMPT)) = PRPLBL(IPRP)
            END IF
            IF (ESROPT(IPRP)) THEN
               NESRT(KSYMPT) = NESRT(KSYMPT) + 1
               LBESRT(KSYMPT,NESRT(KSYMPT)) = PRPLBL(IPRP)
            END IF
            IF (ACROP(IPRP)) THEN
               NACROP(KSYMPT) = NACROP(KSYMPT) + 1
               ACRLB(KSYMPT,NACROP(KSYMPT)) = PRPLBL(IPRP)
            END IF
            IF (BCROP(IPRP)) THEN
               NBCROP(KSYMPT) = NBCROP(KSYMPT) + 1
               BCRLB(KSYMPT,NBCROP(KSYMPT)) = PRPLBL(IPRP)
            END IF
            IF (CCROP(IPRP)) THEN
               NCCROP(KSYMPT) = NCCROP(KSYMPT) + 1
               CCRLB(KSYMPT,NCCROP(KSYMPT)) = PRPLBL(IPRP)
            END IF
            IF (DCROP(IPRP)) THEN
               NDCROP(KSYMPT) = NDCROP(KSYMPT) + 1
               DCRLB(KSYMPT,NDCROP(KSYMPT)) = PRPLBL(IPRP)
            END IF
            IF (ATMOP(IPRP)) THEN
               NATMOP(KSYMPT) = NATMOP(KSYMPT) + 1
               ATMLB(KSYMPT,NATMOP(KSYMPT)) = PRPLBL(IPRP)
            END IF
            IF (BTMOP(IPRP)) THEN
               NBTMOP(KSYMPT) = NBTMOP(KSYMPT) + 1
               BTMLB(KSYMPT,NBTMOP(KSYMPT)) = PRPLBL(IPRP)
            END IF
            IF (CTMOP(IPRP)) THEN
               NCTMOP(KSYMPT) = NCTMOP(KSYMPT) + 1
               CTMLB(KSYMPT,NCTMOP(KSYMPT)) = PRPLBL(IPRP)
            END IF
            IF (ATPOP(IPRP)) THEN
               NATPOP(KSYMPT) = NATPOP(KSYMPT) + 1
               ATPLB(KSYMPT,NATPOP(KSYMPT)) = PRPLBL(IPRP)
            END IF
            IF (BTPOP(IPRP)) THEN
               NBTPOP(KSYMPT) = NBTPOP(KSYMPT) + 1
               BTPLB(KSYMPT,NBTPOP(KSYMPT)) = PRPLBL(IPRP)
            END IF
         END IF
  250 CONTINUE
C
      IF (.NOT.NOS0MX) THEN
         NS0TOT = 0
         DO 220 IPRP = 1,3
            NPRP = 0
            REWIND (LUPROP)
            IF ( FNDLAB(TABLEN(IPRP),LUPROP) ) NPRP = NPRP + 1
            REWIND (LUPROP)
            IF ( FNDLAB(TABVEL(IPRP),LUPROP) ) NPRP = NPRP + 1
            IF (NPRP .EQ. 2) THEN
               CALL READT(LUPROP,NNBASX,WRK)
               IND   = IDAMAX(NNBASX,WRK,1)
               IF (ABS(WRK(IND)).GT.CKMXPR) THEN
                  VAL = D1 + D8*IND
                  IIND  = INT( (SQRT(VAL) - D1INF)/ D2 )
                  IIMAX = IIND*(IIND+1)/2
                  IF ( IIMAX .EQ. IND ) THEN
                     JIND = IIND
                  ELSE
                     IIND = IIND + 1
                     JIND = IND  - IIMAX
                  END IF
                  KSYMPT = MULD2H( ISAO(IIND) , ISAO(JIND) )
               ELSE
                  KSYMPT = 0
                  NINFO  = NINFO + 1
                  WRITE(LUPRI,'(/A,/3A,1P,D15.7,/2A)')'@** INFO **',
     *            '@ELEMENTS ON PROPERTY FILE WITH LABEL: ',
     *            TABVEL(IPRP),' SMALLER THAN :',CKMXPR,
     *            '@PROPERTY CALCULATION NOT CARRIED OUT FOR S(0) IN',
     *            ' MIXED REPRESENTATION'
               END IF
            ELSE
               KSYMPT = 0
               NINFO  = NINFO + 1
               WRITE(LUPRI,'(/5A/4A)')
     *            '@ ** INFO ** LABEL ',TABLEN(IPRP),' OR LABEL ',
     *            TABVEL(IPRP),' NOT FOUND ON LUPROP.',
     *            '@ S(0) IN MIXED REPRESENTATION WILL NOT BE',
     *            ' CALCULATED FOR ',TABLEN(IPRP)(1:1),' OPERATOR.'
            END IF
            IF ( (KSYMPT .GT. 0) .AND. (KSYMPT .LE. NSYM) ) THEN
               NS0TOT = NS0TOT + 1
               NGPS0(KSYMPT) = NGPS0(KSYMPT) + 1
               LBLS0(KSYMPT,NGPS0(KSYMPT),1)  = TABLEN(IPRP)
               LBLS0(KSYMPT,NGPS0(KSYMPT),2)  = TABVEL(IPRP)
            END IF
 220     CONTINUE
         IF (NS0TOT .EQ. 0) NOS0MX = .TRUE.
      END IF
C
      CALL GPCLOSE(LUPROP,'KEEP')
C
      IF (IPRRSP.GT.5) THEN
        WRITE(LUPRI,'(/A/)')  ' RSPSYM: number of property vectors :'
        WRITE(LUPRI,'(A,8I5)')' NGPPP(I) ,I=1,NSYM',(NGPPP(I) ,I=1,NSYM)
        WRITE(LUPRI,'(A,8I5)')' NGPLR(I) ,I=1,NSYM',(NGPLR(I) ,I=1,NSYM)
        WRITE(LUPRI,'(A,8I5)')' NGPC6(I) ,I=1,NSYM',(NGPC6(I) ,I=1,NSYM)
        WRITE(LUPRI,'(A,8I5)')' NAQROP(I),I=1,NSYM',(NAQROP(I),I=1,NSYM)
        WRITE(LUPRI,'(A,8I5)')' NBQROP(I),I=1,NSYM',(NBQROP(I),I=1,NSYM)
        WRITE(LUPRI,'(A,8I5)')' NCQROP(I),I=1,NSYM',(NCQROP(I),I=1,NSYM)
        WRITE(LUPRI,'(A,8I5)')' NASMOP(I),I=1,NSYM',(NASMOP(I),I=1,NSYM)
        WRITE(LUPRI,'(A,8I5)')' NBSMOP(I),I=1,NSYM',(NBSMOP(I),I=1,NSYM)
        WRITE(LUPRI,'(A,8I5)')' NESRS(I) ,I=1,NSYM',(NESRS(I) ,I=1,NSYM)
        WRITE(LUPRI,'(A,8I5)')' NESRT(I) ,I=1,NSYM',(NESRT(I) ,I=1,NSYM)
        WRITE(LUPRI,'(A,8I5)')' NACROP(I),I=1,NSYM',(NACROP(I),I=1,NSYM)
        WRITE(LUPRI,'(A,8I5)')' NBCROP(I),I=1,NSYM',(NBCROP(I),I=1,NSYM)
        WRITE(LUPRI,'(A,8I5)')' NCCROP(I),I=1,NSYM',(NCCROP(I),I=1,NSYM)
        WRITE(LUPRI,'(A,8I5)')' NDCROP(I),I=1,NSYM',(NDCROP(I),I=1,NSYM)
        WRITE(LUPRI,'(A,8I5)')' NATMOP(I),I=1,NSYM',(NATMOP(I),I=1,NSYM)
        WRITE(LUPRI,'(A,8I5)')' NBTMOP(I),I=1,NSYM',(NBTMOP(I),I=1,NSYM)
        WRITE(LUPRI,'(A,8I5)')' NCTMOP(I),I=1,NSYM',(NCTMOP(I),I=1,NSYM)
        WRITE(LUPRI,'(A,8I5)')' NATPOP(I),I=1,NSYM',(NATPOP(I),I=1,NSYM)
        WRITE(LUPRI,'(A,8I5)')' NBTPOP(I),I=1,NSYM',(NBTPOP(I),I=1,NSYM)
      END IF
      NTOTLR = 0
      NTOTC6 = 0
      DO 231 ISYM = 1,NSYM
         NTOTLR = NTOTLR + NGPLR(ISYM)
         NTOTC6 = NTOTC6 + NGPC6(ISYM)
  231 CONTINUE
      NTOTRS = NTOTLR + NTOTC6
      IF ( (NTOTRS.GT.0) .AND. (HYPCAL.OR.SOMOM.OR.EXMOM) ) THEN
         WRITE(LUPRI,'(/A/)')
     *   '@ INPUT SPECIFICATION FOR BOTH LINEAR AND QUADRATIC RESPONSE'
     *   ,'CALCULATIONS NOT ALLOWED'
         IF (NTOTLR.GT.0) WRITE(LUPRI,'(A)')
     *      '@ OPERATORS SPECIFIED IN *RSPLR SECTION'
         IF (NTOTC6.GT.0) WRITE(LUPRI,'(A)')
     *      '@ OPERATORS SPECIFIED IN *RSPC6 SECTION'
         WRITE(LUPRI,'(A)')
     *      '@ THESE OPERATORS MUST NOT BE SPECIFIED TOGETHER WITH '
         IF (HYPCAL) WRITE(LUPRI,'(/A)')'@ .HYPCAL'
         IF (SOMOM)  WRITE(LUPRI,'(/A)')'@ .SOMOM '
         IF (EXMOM)  WRITE(LUPRI,'(/A)')'@ .EXMOM '
         CALL QUIT('RSPSYM: BOTH LR AND QR CALC SPECIFIED')
      END IF
      IF ( (CRCAL.OR.TOMOM.OR.TPAMP).AND.(HYPCAL.OR.SOMOM.OR.EXMOM) )
     *   WRITE(LUPRI,'(/A/)')
     *   '@ INPUT SPECIFICATION FOR BOTH QUADRATIC AND CUBIC RESPONSE'
     *   ,'CALCULATIONS NOT ALLOWED'
      IF ( (NTOTRS.GT.0) .AND. (CRCAL.OR.TOMOM) ) THEN
         WRITE(LUPRI,'(/A/)')
     *   '@ INPUT SPECIFICATION FOR BOTH LINEAR AND CUBIC RESPONSE'
     *   ,'CALCULATIONS NOT ALLOWED'
      END IF
C
      IF (CRCAL) THEN
C
C SET UP POINTERS FOR LINEAR RESPONSE EQUATIONS
C
         CALL CR1IND
         CALL CR2IND
      END IF
C
      IF (TOMOM) THEN
C
C SET UP POINTERS FOR LINEAR RESPONSE EQUATIONS
C AND EXCITATION ENERGIES
C
         CALL TM1IND
         CALL TM2IND
      END IF
C
      IF (TPAMP) THEN
C
C SET UP POINTERS FOR LINEAR RESPONSE EQUATIONS
C AND EXCITATION ENERGIES
C
         CALL TP1IND
         CALL TP2IND
      END IF
C
      IF (EXMOM) THEN
         WRITE(LUPRI,'(/A)')
     * '** CALCULATION OF TRANSITION MOMENTS BETWEEN EXCITED STATES **'
         CALL ANASYM(NGPPP,NPPCNV,NPPCNV,WRK)
C
C SET UP POINTERS FOR EXCITATION ENERGIES
C
         CALL EXIIND(NPPCNV)
C
C SET UP POINTERS FOR LINEAR RESPONSE EQUATIONS
C
         CALL EXMIND
      END IF
      DO 230 ISYM = 1,NSYM
         IF (NPPCNV(ISYM).GT.0) THEN
            WRITE(LUPRI,'(/I7,2A,I5)')
     *      NPPCNV(ISYM),' Excitation energies',
     *      ' are calculated for symmetry no.',ISYM
            IF (NGPPP(ISYM).GT.0) THEN
               WRITE(LUPRI,'(/I7,2A/)')
     *         NGPPP(ISYM),' property residues',
     *         ' are calculated with labels:'
               DO 225 IRESID = 1,NGPPP(ISYM)
                  WRITE(LUPRI,'(15X,A)')LBLPP(ISYM,IRESID)
 225           CONTINUE
            END IF
         END IF
 230  CONTINUE
      IF (IPRRSP.GT.10) THEN
         WRITE(LUPRI,'(/A/)')'  SYMMETRY   NO.CONV  NO.START NO.SIMULT'
         DO 105 ISYM=1,NSYM
            WRITE(LUPRI,'(4I10)')
     *           ISYM,NPPCNV(ISYM),NPPSTV(ISYM),NPPSIM(ISYM)
 105     CONTINUE
      END IF
      IF (HYPCAL) THEN
C
C SET UP POINTERS FOR LINEAR RESPONSE EQUATIONS
C
         CALL HYPIND
      END IF
      IF (SOMOM) THEN
C
C SET UP POINTERS FOR EXCITATION ENERGIES
C
         CALL EXIIND(NSMCNV)
C
C SET UP POINTERS FOR LINEAR RESPONSE EQUATIONS
C
         CALL SOMIND
      END IF
C
      IF (NTOTRS.GT.0) THEN
         DO 330 ISYM = 1,NSYM
            IF (NGPLR(ISYM).GT.0) THEN
               WRITE(LUPRI,'(/I5,A,I5,A/)')NGPLR(ISYM),
     *         ' second order properties calculated with symmetry no.',
     *         ISYM,' and labels:'
               DO 325 IRESID = 1,NGPLR(ISYM)
                  WRITE(LUPRI,'(10X,A)')LBLLR(ISYM,IRESID)
 325           CONTINUE
            END IF
 330     CONTINUE
         DO 430 ISYM = 1,NSYM
            IF (NGPC6(ISYM).GT.0) THEN
               WRITE(LUPRI,'(/I5,A,I5,A/)')NGPC6(ISYM),
     *       ' C6/C8 module calculated for operators with symmetry no.',
     *         ISYM,' and labels:'
               DO 526 IRESID = 1,NGPC6(ISYM)
                  WRITE(LUPRI,'(10X,A)')LBLC6(ISYM,IRESID)
 526           CONTINUE
            END IF
 430     CONTINUE
      ELSE
         IF (CRCAL) THEN
            CALL AN2SYM(NACROP,NBCROP,NCCROP,NDCROP,WRK)
            DO 480 ISYM = 1,NSYM
               IF (NACROP(ISYM).GT.0) THEN
                  WRITE(LUPRI,'(/I5,A,I5,A/)')NACROP(ISYM),
     *            ' A OPERATORS OF SYMMETRY NO:',
     *            ISYM,' AND LABELS:'
                  DO 485 IRESID = 1,NACROP(ISYM)
                     WRITE(LUPRI,'(10X,A)')ACRLB(ISYM,IRESID)
 485              CONTINUE
               END IF
 480        CONTINUE
            DO 481 ISYM = 1,NSYM
               IF (NBCROP(ISYM).GT.0) THEN
                  WRITE(LUPRI,'(/I5,A,I5,A/)')NBCROP(ISYM),
     *            ' B OPERATORS OF SYMMETRY NO:',
     *            ISYM,' AND LABELS:'
                  DO 486 IRESID = 1,NBCROP(ISYM)
                     WRITE(LUPRI,'(10X,A)')BCRLB(ISYM,IRESID)
 486              CONTINUE
               END IF
 481        CONTINUE
            DO 482 ISYM = 1,NSYM
               IF (NCCROP(ISYM).GT.0) THEN
                  WRITE(LUPRI,'(/I5,A,I5,A/)')NCCROP(ISYM),
     *            ' C OPERATORS OF SYMMETRY NO:',
     *            ISYM,' AND LABELS:'
                  DO 487 IRESID = 1,NCCROP(ISYM)
                     WRITE(LUPRI,'(10X,A)')CCRLB(ISYM,IRESID)
 487              CONTINUE
               END IF
 482        CONTINUE
            DO 483 ISYM = 1,NSYM
               IF (NDCROP(ISYM).GT.0) THEN
                  WRITE(LUPRI,'(/I5,A,I5,A/)')NDCROP(ISYM),
     *            ' D OPERATORS OF SYMMETRY NO:',
     *            ISYM,' AND LABELS:'
                  DO 484 IRESID = 1,NDCROP(ISYM)
                     WRITE(LUPRI,'(10X,A)')DCRLB(ISYM,IRESID)
 484              CONTINUE
               END IF
 483        CONTINUE
         END IF
         IF (TOMOM) THEN
            CALL AN2SYM(NATMOP,NBTMOP,NCTMOP,NTMCNV,WRK)
            DO 790 ISYM = 1,NSYM
               IF (NATMOP(ISYM).GT.0) THEN
                  WRITE(LUPRI,'(/I5,A,I5,A/)')NATMOP(ISYM),
     *            ' A OPERATORS OF SYMMETRY NO:',
     *            ISYM,' AND LABELS:'
                  DO 795 IRESID = 1,NATMOP(ISYM)
                     WRITE(LUPRI,'(10X,A)')ATMLB(ISYM,IRESID)
 795              CONTINUE
               END IF
 790        CONTINUE
            DO 791 ISYM = 1,NSYM
               IF (NBTMOP(ISYM).GT.0) THEN
                  WRITE(LUPRI,'(/I5,A,I5,A/)')NBTMOP(ISYM),
     *            ' B OPERATORS OF SYMMETRY NO:',
     *            ISYM,' AND LABELS:'
                  DO 796 IRESID = 1,NBTMOP(ISYM)
                     WRITE(LUPRI,'(10X,A)')BTMLB(ISYM,IRESID)
 796              CONTINUE
               END IF
 791        CONTINUE
            DO 798 ISYM = 1,NSYM
               IF (NCTMOP(ISYM).GT.0) THEN
                  WRITE(LUPRI,'(/I5,A,I5,A/)')NCTMOP(ISYM),
     *            ' C OPERATORS OF SYMMETRY NO:',
     *            ISYM,' AND LABELS:'
                  DO 797 IRESID = 1,NCTMOP(ISYM)
                     WRITE(LUPRI,'(10X,A)')CTMLB(ISYM,IRESID)
 797              CONTINUE
               END IF
 798        CONTINUE
         END IF
         IF (TPAMP) THEN
            CALL AN2SYM(NATPOP,NBTPOP,NTPCN1,NTPCN2,WRK)
            DO 390 ISYM = 1,NSYM
               IF (NATPOP(ISYM).GT.0) THEN
                  WRITE(LUPRI,'(/I5,A,I5,A/)')NATPOP(ISYM),
     *            ' A OPERATORS OF SYMMETRY NO:',
     *            ISYM,' AND LABELS:'
                  DO 395 IRESID = 1,NATPOP(ISYM)
                     WRITE(LUPRI,'(10X,A)')ATPLB(ISYM,IRESID)
 395              CONTINUE
               END IF
 390        CONTINUE
            DO 391 ISYM = 1,NSYM
               IF (NBTPOP(ISYM).GT.0) THEN
                  WRITE(LUPRI,'(/I5,A,I5,A/)')NBTPOP(ISYM),
     *            ' B OPERATORS OF SYMMETRY NO:',
     *            ISYM,' AND LABELS:'
                  DO 396 IRESID = 1,NBTPOP(ISYM)
                     WRITE(LUPRI,'(10X,A)')BTPLB(ISYM,IRESID)
 396              CONTINUE
               END IF
 391        CONTINUE
         END IF
         IF (HYPCAL) THEN
            CALL ANASYM(NAQROP,NBQROP,NCQROP,WRK)
            DO 490 ISYM = 1,NSYM
               IF (NAQROP(ISYM).GT.0) THEN
                  WRITE(LUPRI,'(/I5,A,I5,A/)')NAQROP(ISYM),
     *            ' A OPERATORS OF SYMMETRY NO:',
     *            ISYM,' AND LABELS:'
                  DO 495 IRESID = 1,NAQROP(ISYM)
                     WRITE(LUPRI,'(10X,A)')AQRLB(ISYM,IRESID)
 495              CONTINUE
               END IF
 490        CONTINUE
            DO 491 ISYM = 1,NSYM
               IF (NBQROP(ISYM).GT.0) THEN
                  WRITE(LUPRI,'(/I5,A,I5,A/)')NBQROP(ISYM),
     *            ' B OPERATORS OF SYMMETRY NO:',
     *            ISYM,' AND LABELS:'
                  DO 496 IRESID = 1,NBQROP(ISYM)
                     WRITE(LUPRI,'(10X,A)')BQRLB(ISYM,IRESID)
 496              CONTINUE
               END IF
 491        CONTINUE
            DO 492 ISYM = 1,NSYM
               IF (NCQROP(ISYM).GT.0) THEN
                  WRITE(LUPRI,'(/I5,A,I5,A/)')NCQROP(ISYM),
     *            ' C OPERATORS OF SYMMETRY NO:',
     *            ISYM,' AND LABELS:'
                  DO 497 IRESID = 1,NCQROP(ISYM)
                     WRITE(LUPRI,'(10X,A)')CQRLB(ISYM,IRESID)
 497              CONTINUE
               END IF
 492        CONTINUE
         END IF
         IF (SOMOM) THEN
            CALL ANASYM(NASMOP,NBSMOP,NSMCNV,WRK)
            DO 590 ISYM = 1,NSYM
               IF (NASMOP(ISYM).GT.0) THEN
                  WRITE(LUPRI,'(/I5,A,I5,A/)')NASMOP(ISYM),
     *            ' A OPERATORS OF SYMMETRY NO:',
     *            ISYM,' AND LABELS:'
                  DO 595 IRESID = 1,NASMOP(ISYM)
                     WRITE(LUPRI,'(10X,A)')ASMLB(ISYM,IRESID)
 595              CONTINUE
               END IF
 590        CONTINUE
            DO 591 ISYM = 1,NSYM
               IF (NBSMOP(ISYM).GT.0) THEN
                  WRITE(LUPRI,'(/I5,A,I5,A/)')NBSMOP(ISYM),
     *            ' B OPERATORS OF SYMMETRY NO:',
     *            ISYM,' AND LABELS:'
                  DO 596 IRESID = 1,NBSMOP(ISYM)
                     WRITE(LUPRI,'(10X,A)')BSMLB(ISYM,IRESID)
 596              CONTINUE
               END IF
 591        CONTINUE
         END IF
      END IF
      IF ( CRCAL.OR.TOMOM ) THEN
       IF (IPRRSP.GT.40) THEN
         WRITE(LUPRI,'(/A,I5,A)')
     *   ' LIST OF NLRLBL= ',NLRLBL,
     *   ' SINGLET LINEAR EQUATIONS FOR CUBIC RESPONSE'
         WRITE(LUPRI,'(/A)')'   I  CRLBL(I)  CRFREQ(I)  ISYMCR(I) '
         DO 890 I = 1,NLRLBL
            WRITE(LUPRI,'(I5,5X,A,1P,D14.6,I12)')
     *      I,CRLBL(I),CRFREQ(I),ISYMCR(I)
  890    CONTINUE
       END IF
      END IF
      IF ( EXMOM .OR. SOMOM .OR. HYPCAL ) THEN
       IF (IPRRSP.GT.40) THEN
         WRITE(LUPRI,'(/A,I5,A)')
     *   ' LIST OF NLRLBL= ',NLRLBL,
     *   ' SINGLET LINEAR EQUATIONS FOR QUADRATIC RESPONSE'
         WRITE(LUPRI,'(/A)')'   I  QRLBL(I)  QRFREQ(I)  ISYMQR(I) '
         DO 990 I = 1,NLRLBL
            WRITE(LUPRI,'(I5,5X,A,1P,D14.6,I12)')
     *      I,QRLBL(I),QRFREQ(I),ISYMQR(I)
  990    CONTINUE
         IF ( TRPFLG .AND. NTRLBL.GT.0 ) THEN
            WRITE(LUPRI,'(/A,I5,A)')
     *      ' LIST OF NTRLBL= ',NTRLBL,
     *      ' TRIPLET LINEAR EQUATIONS FOR QUADRATIC RESPONSE'
            WRITE(LUPRI,'(/A)')'   I  TRLBL(I)  TRFREQ(I)  ISYMTR(I) '
            DO 991 I = 1,NTRLBL
               WRITE(LUPRI,'(I5,5X,A,1P,D14.6,I12)')
     *         I,TRLBL(I),TRFREQ(I),ISYMTR(I)
  991       CONTINUE
         END IF
       END IF
      END IF
      IF ( ESRCAL ) THEN
         DO 436 ISYM = 1,NSYM
            IF (NESRS(ISYM).GT.0) THEN
               WRITE(LUPRI,'(/I5,A,I5,A/)')NESRS(ISYM),
     *         ' ESR CALC. : SINGLET OP.  SYMMETRY NO:',
     *         ISYM,' AND LABEL:'
               DO 536 IRESID = 1,NESRS(ISYM)
                  WRITE(LUPRI,'(10X,A)')LBESRS(ISYM,IRESID)
 536           CONTINUE
            END IF
 436     CONTINUE
         DO 437 ISYM = 1,NSYM
            IF (NESRT(ISYM).GT.0) THEN
               WRITE(LUPRI,'(/I5,A,I5,A/)')NESRT(ISYM),
     *         ' ESR CALC. : TRIPLET OP.  SYMMETRY NO:',
     *         ISYM,' AND LABEL:'
               DO 537 IRESID = 1,NESRT(ISYM)
                  WRITE(LUPRI,'(10X,A)')LBESRT(ISYM,IRESID)
 537           CONTINUE
            END IF
 437     CONTINUE
      END IF
C
C *** END OF RSPSYM
C
      CALL QEXIT('RSPSYM')
      RETURN
      END
C  /* Deck anasym */
      SUBROUTINE ANASYM(NAOPNU,NBOPNU,NCOPNU,IWRK)
C
C THIS ROUTINE ANALYZES IF OPERATORS THAT ARE SPECIFIED IN QUADRATIC
C RESPONSE CALCULATIONS HAVE NO MATCHING PARTNERS DUE TO SYMMETRY.
C SUCH OPERATORS ARE REMOVED FROM THE CALCULATION IN THIS ROUTINE.
C
#include "implicit.h"
C
#include "priunit.h"
#include "inforb.h"
#include "infrsp.h"
C
      DIMENSION NAOPNU(*),NBOPNU(*),NCOPNU(*),IWRK(NSYM,*)
      IF (IPRRSP.GT.5) THEN
         WRITE(LUPRI,'(/A)')
     *   ' ANASYM: NUMBER OF QR PROPERTY VECTORS BEFORE REDUCTION'
         WRITE(LUPRI,'(/A)')' NAOPNU(I),I=1,NSYM '
         WRITE(LUPRI,'(8I5)')( NAOPNU(I),I=1,NSYM)
         WRITE(LUPRI,'(/A)')' NBOPNU(I),I=1,NSYM '
         WRITE(LUPRI,'(8I5)')( NBOPNU(I),I=1,NSYM)
         WRITE(LUPRI,'(/A)')' NCOPNU(I),I=1,NSYM '
         WRITE(LUPRI,'(8I5)')( NCOPNU(I),I=1,NSYM)
      END IF
      DO 100 I = 1,NSYM
         IWRK(I,1) = 0
         IWRK(I,2) = 0
         IWRK(I,3) = 0
 100  CONTINUE
      DO 200 IASYM = 1,NSYM
         IF ( NAOPNU(IASYM).LE.0 ) GO TO 200
         DO 300 IBSYM = 1,NSYM
            IF (NBOPNU(IBSYM).LE.0) GO TO 300
            ICSYM = MULD2H(IASYM,IBSYM)
            IF(NCOPNU(ICSYM).GT.0) THEN
               IWRK(IASYM,1) = NAOPNU(IASYM)
               IWRK(IBSYM,2) = NBOPNU(IBSYM)
               IWRK(ICSYM,3) = NCOPNU(ICSYM)
            END IF
 300     CONTINUE
 200  CONTINUE
      NAOPTO = 0
      NBOPTO = 0
      NCOPTO = 0
      DO 400 ISYM = 1,NSYM
         IF (IWRK(ISYM,1).NE.NAOPNU(ISYM)) THEN
           WRITE(LUPRI,'(/I5,A,I5,A)')
     *     NAOPNU(ISYM),' A OPERATORS OF SYMMETRY',ISYM,' NOT INCLUDED'
           NAOPNU(ISYM) = IWRK(ISYM,1)
         END IF
         IF (IWRK(ISYM,2).NE.NBOPNU(ISYM)) THEN
           WRITE(LUPRI,'(/I5,A,I5,A)')
     *     NBOPNU(ISYM),' B OPERATORS OF SYMMETRY',ISYM,' NOT INCLUDED'
           NBOPNU(ISYM) = IWRK(ISYM,2)
         END IF
         IF (IWRK(ISYM,3).NE.NCOPNU(ISYM)) THEN
           WRITE(LUPRI,'(/I5,A,I5,A)')
     *     NCOPNU(ISYM),' C OPERATORS OF SYMMETRY',ISYM,' NOT INCLUDED'
           NCOPNU(ISYM) = IWRK(ISYM,3)
         END IF
         NAOPTO = NAOPTO + NAOPNU(ISYM)
         NBOPTO = NBOPTO + NBOPNU(ISYM)
         NCOPTO = NCOPTO + NCOPNU(ISYM)
 400  CONTINUE
      IF (IPRRSP.GT.5) THEN
         WRITE(LUPRI,'(/A)')
     *   ' ANASYM: NUMBER OF QR PROPERTY VECTORS AFTER REDUCTION'
         WRITE(LUPRI,'(/A)')' NAOPNU(I),I=1,NSYM '
         WRITE(LUPRI,'(8I5)')( NAOPNU(I),I=1,NSYM)
         WRITE(LUPRI,'(/A)')' NBOPNU(I),I=1,NSYM '
         WRITE(LUPRI,'(8I5)')( NBOPNU(I),I=1,NSYM)
         WRITE(LUPRI,'(/A)')' NCOPNU(I),I=1,NSYM '
         WRITE(LUPRI,'(8I5)')( NCOPNU(I),I=1,NSYM)
      END IF
      NOPTOT = NAOPTO + NBOPTO + NCOPTO
      IF (NOPTOT.EQ.0) THEN
         NINFO = NINFO + 1
         WRITE(LUPRI,'(3(/A))')
     *   ' INFO: Symmetry of specified operators does not match.',
     *   ' INFO: Response properties are thus trivially ZERO.',
     *   ' INFO: Welcome back with new input.'
C
      END IF
      RETURN
      END
C  /* Deck rspvar */
      SUBROUTINE RSPVAR(UDV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
C
C Purpose: this subroutine must be called in the RESPONSE module
C each time the operator symmetry KSYMOP changes;
C here variables which depend on KSYMOP are defined. /hjaaj sep 2003
C
#include "implicit.h"
#include "dummy.h"
C
      DIMENSION UDV(NASHT,*),FOCK(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
C
#include "maxorb.h"
#include "maxash.h"
#include "priunit.h"
#include "infinp.h"
#include "infind.h"
#include "infpri.h"
#include "infdim.h"
#include "infvar.h"
#include "inforb.h"
#include "infopt.h"
#include "inftap.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "rspprp.h"
#include "infpp.h"
#include "inflr.h"
#include "infsop.h"
#include "channel.h"
C
      CHARACTER*8 OPLBL(8)
      DATA OPLBL/'EXOPSYM1','EXOPSYM2','EXOPSYM3','EXOPSYM4',
     *           'EXOPSYM5','EXOPSYM6','EXOPSYM7','EXOPSYM8'/
C
      PARAMETER ( D1 = 1.0D0, THR_PRUNE = 1.0D-4 )
C
      CALL QENTER('RSPVAR')
C
      KSYMST = MULD2H(IREFSY,KSYMOP)
C
C READ IN ORBITAL EXITATION OPERATOR
C
      IF (RSPCI) THEN
         NWOPT  = 0
         KZWOPT = 0
         KZYWOP = 0
      ELSE
!        IF ( KSYMOP .NE. JWOPSY ) THEN
            JWOPSY = KSYMOP
            LUINDX = -1
            CALL GPOPEN(LUINDX,'LUINDF','UNKNOWN',' ','UNFORMATTED',
     &                  IDUMMY,.FALSE.)
            REWIND LUINDX
            CALL MOLLAB(OPLBL(KSYMOP),LUINDX,LUPRI)
            READ (LUINDX) NWOPT,NWOPH
            CALL READI(LUINDX,(2*NWOPT),JWOP)
            CALL GPCLOSE(LUINDX,'KEEP')
C
C           Prune some classes of redundant orbital rotations from
C           the orbital rotation vector JWOP:
C              UDV(u,u) = 2
C              UDV(v,v) = 0
C              UDV(u,u) = UDV(v,v) [only possible for RAS]
C           Threshold for redundancy test: THR_PRUNE
C           / hjaaj June 2010
            MWOPT = 0
            DO IWOPT = 1,NWOPT
               K = JWOP(1,IWOPT)
               L = JWOP(2,IWOPT)
               ITYPK = IOBTYP(K)
               ITYPL = IOBTYP(L)
               KL_is_OK = 0
               IF (ITYPK .EQ. JTINAC) THEN
                  IF (ITYPL .EQ. JTACT) THEN
                  ! inact-act rotation
                     NWL = ISW(L) - NISHT
                     IF ( (2.0D0-UDV(NWL,NWL)) .GT. THR_PRUNE) THEN
                        KL_is_OK = 1
                     ELSE IF (IPRRSP .GE. 2) THEN
                        WRITE (LUPRI,*)
     &                  'INFO: Inact-act orbital rotation',K,L,' pruned'
                     END IF
                  ELSE
                     KL_is_OK = 1
                  END IF
               ELSE IF (ITYPL .EQ. JTACT) THEN
               ! act-act rotation (i.e. RAS case)
                  NWK = ISW(K) - NISHT
                  NWL = ISW(L) - NISHT
                  IF (ABS(UDV(NWK,NWK)-UDV(NWL,NWL)) .GT. THR_PRUNE)
     &            THEN
                     KL_is_OK = 1
                  ELSE IF (IPRRSP .GE. 2) THEN
                     WRITE (LUPRI,*)
     &               'INFO:   Act-act orbital rotation',K,L,' pruned'
                  END IF
               ELSE
               ! act-sec rotation
                  NWK = ISW(K) - NISHT
                  IF (UDV(NWK,NWK) .GT. THR_PRUNE) THEN
                     KL_is_OK = 1
                  ELSE IF (IPRRSP .GE. 2) THEN
                     WRITE (LUPRI,*)
     &               'INFO:   Act-sec orbital rotation',K,L,' pruned'
                  END IF
               END IF
               IF (KL_is_OK .EQ. 1) THEN
                  MWOPT = MWOPT + 1
                  JWOP(1,MWOPT) = K
                  JWOP(2,MWOPT) = L
               END IF
            END DO
            IF (MWOPT .NE. NWOPT) THEN
               WRITE (LUPRI,'(/A,I5,A)') 'INFO:',NWOPT-MWOPT,
     &         ' redundant orbital rotations pruned.'
               NWOPT = MWOPT
            END IF
!        END IF
C
         KZWOPT = NWOPT
         KZYWOP = 2*KZWOPT
         NWOPH  = NWOPT ! do not include redundant act-act rotations
                        ! in response, even when we call SIRIUS routines
                        ! /hjaaj June 2010
         IF (CHANNEL_CALC ) CALL CHANNEL_VAR()
         IF (CHANNEL_VCALC) CALL CHANNEL_VIR()
      END IF
C
C READ IN CONFIGURATION SPACE VARIABLES
C
      IF (TDHF) THEN
         IF (SOPPA) THEN
            IF (TRPLET) THEN
               KZCONF = N2P2HT(KSYMOP)
            ELSE
               KZCONF = N2P2HS(KSYMOP)
            ENDIF
            NCONF  = 0
            KZYCON = 2*KZCONF
         ELSE
            NCONF  = 0
            KZCONF = 0
            KZYCON = 0
         ENDIF
      ELSE
         CALL SETCI2(NCONFX,KSYMST,TRPFLG,0)
C...     TRPFLG and not TRPLET because we must use determinants
C        in qr and cr if any operator is triplet ! /hjaaj-may2000
         KZCONF = NCONFX
         KZYCON = 2*KZCONF
      END IF
C
C TOTAL NUMBER OF VARIABLES
C
      KZVAR  = KZWOPT + KZCONF
      KZYVAR = 2*KZVAR
      IF (SOPPA) THEN
         IF (IPRRSP .GE. 2) WRITE(LUPRI,'(/4(/A,I8))')
     *      ' Perturbation symmetry.     KSYMOP:',KSYMOP,
     *      ' p-h variables.             KZWOPT:',KZWOPT,
     *      ' 2p-2h variables.           KZCONF:',KZCONF,
     *      ' Total number of variables. KZVAR :',KZVAR
         IF (TRPLET) THEN
            IF (IPRRSP .GE. 6) WRITE(LUPRI,'(3(/A,I8))')
     *         ' T(1) 2p-2h variables.         :',N12P2H(KSYMOP),
     *         ' T(2) 2p-2h variables.         :',N22P2H(KSYMOP),
     *         ' T(3) 2p-2h variables.         :',N32P2H(KSYMOP)
         ELSE
            IF (IPRRSP .GE. 6) WRITE(LUPRI,'(2(/A,I8))')
     *         ' R(1) 2p-2h variables.         :',NS2P2H(KSYMOP),
     *         ' R(2) 2p-2h variables.         :',NT2P2H(KSYMOP)
         ENDIF
         IF (IPRRSP.GT.45) WRITE(LUPRI,'(/A,3I8)')
     *      ' KZYWOP,KZYCON,KZYVAR:',KZYWOP,KZYCON,KZYVAR
      ELSE
         IF (IPRRSP .GE. 2) WRITE(LUPRI,'(//A,I8/A,L8,3(/A,I8))')
     *      ' Perturbation symmetry.     KSYMOP:',KSYMOP,
     *      ' Perturbation spin symmetry.TRPLET:',TRPLET,
     *      ' Orbital variables.         KZWOPT:',KZWOPT,
     *      ' Configuration variables.   KZCONF:',KZCONF,
     *      ' Total number of variables. KZVAR :',KZVAR
         IF (IPRRSP.GT.45) WRITE(LUPRI,'(/A,3I8)')
     *      ' KZYWOP,KZYCON,KZYVAR:',KZYWOP,KZYCON,KZYVAR
      END IF
C
C DEFINE IF THE REFERENCE STATE FUNCTION HAS TO BE READ IN AS
C THE FIRST TRIAL VECTOR ON LURSP3
C
      IF ( (KZCONF .GT. 0) .AND. (KSYMOP .EQ. 1) .AND.
     *    (.NOT.SOPPA)) THEN
Chj-000406: now always DETERM for TRPLET
Chj  *    ( DETERM .OR. .NOT.TRPLET) .AND. (.NOT.SOPPA)) THEN
         KOFFTY = 1
      ELSE
         KOFFTY = 0
      ENDIF
C
C DEFINE NUMBER OF START ,SIMULTANEOUS AND CONVERGED VECTORS
C
      IF (NPPCNV(KSYMOP) .GT. KZVAR-KOFFTY) THEN
         WRITE (LUPRI,91) NPPCNV(KSYMOP),KZVAR-KOFFTY
 91      FORMAT(/' ** RSPVAR ** Number of PP start vectors',I5,
     &   /'   exceeds the maximum number of variables',I5,
     &   /'   Number of start vectors reset to number of variables')
         NPPCNV(KSYMOP) = KZVAR-KOFFTY
      END IF
      KEXCNV = NPPCNV(KSYMOP)
      KEXSIM = MIN(NPPSIM(KSYMOP), KZVAR)
      KEXSTV = MIN(NPPSTV(KSYMOP), KZVAR)
C
C CALCULATE DIAGONAL CONFIGURATION PART OF E[2]
C
         CALL REWSPL(LURSP4)
         IF (KZCONF.GE.1) THEN
            KECDIA = 1
            KWRK1  = KECDIA + KZCONF
            LWRK1  = LWRK   - KWRK1
            IF (LWRK1.LT.0) CALL ERRWRK('RSPVAR 1',KWRK1-1,LWRK)
            NCONSV = NCONF
            NCONF  = KZCONF
            IF (SOPPA) THEN
               CALL SOPDIAG(KSYMOP,FC,WRK(KECDIA),XINDX(KABSAD),
     *                      XINDX(KABTAD),XINDX(KIJSAD),XINDX(KIJTAD),
     *                      XINDX(KIJ1AD),XINDX(KIJ2AD),XINDX(KIJ3AD),
     *                      WRK(KWRK1),LWRK1,IPRRSP)
            ELSE
C              LUIT2 = 0 means do not write CI diag. to file
               LUIT2S = LUIT2
               LUIT2 = 0
               LSYM_SAVE = LSYM
Chj            dirty trick: LSYM = 0 forces det. diag. for triplet
               IF (TRPLET) LSYM = 0
               CALL CIDIAG(KSYMST,.FALSE.,FCAC,H2AC,XINDX,
     *                     WRK(KECDIA),WRK(KWRK1),LWRK1)
C              CALL CIDIAG(ICSYM,NOH2,FCAC,H2AC,XNDXCI,DIAGC,WRK,LFREE)
               LUIT2 = LUIT2S
               LSYM  = LSYM_SAVE
               NCONF  = NCONSV
               DO 50 I=1,KZCONF
                  WRK(KECDIA-1+I) = WRK(KECDIA-1+I) - EACTIV
 50            CONTINUE
            ENDIF
            CALL WRITT(LURSP4,KZCONF,WRK(KECDIA))
            IF (IPRRSP.GT.110) THEN
               WRITE(LUPRI,'(/A)')
     *         ' DIAGONAL CONFIGURATION PART OF E[2] '
               CALL OUTPUT(WRK(KECDIA),1,KZCONF,1,1,KZCONF,1,-1,LUPRI)
            END IF
         END IF
C
C CALCULATE FOCK MATRIX CONTRIBUTION TO DIAGONAL ORBITAL HESSIAN
C
         KEODIA = 1
         KSODIA = KEODIA + KZWOPT
         KSCDIA = KSODIA + KZWOPT
C
         IF (KZWOPT.GT.0) THEN
             IF (IPRRSP.GT.100) THEN
               IF (NASHT.GT.0) THEN
                  WRITE(LUPRI,'(/A)')
     *            ' ONE ELECTRON ACTIVE DENSITY MATRIX - UNPACKED'
                  CALL OUTPUT(UDV,1,NASHT,1,NASHT,NASHT,NASHT,-1,LUPRI)
                  WRITE(LUPRI,'(/A)') ' FOCK VALENCE MATRIX'
                  CALL OUTPKB(FV,NORB,NSYM,-1,LUPRI)
               END IF
               WRITE(LUPRI,'(/A)') ' FOCK CORE MATRIX'
               CALL OUTPKB(FC,  NORB,NSYM,-1,LUPRI)
               WRITE(LUPRI,'(/A)') ' TOTAL FOCK  MATRIX'
               CALL OUTPTB(FOCK,NORB,NSYM,-1,LUPRI)
             END IF
C
            CALL ORBDIA(FOCK,FC,FV,UDV,WRK(KEODIA),WRK(KSODIA))
C           CALL ORBDIA(FOCK,FC,FV,UDV,EDIA,SDIA)
C
C WRITE diagonals of orbital part of E[2]-matrix needed for RSPLIN
C
            CALL WRITT(LURSP4,KZWOPT,WRK(KEODIA))
C
            IF (IPRRSP.GT.100) THEN
               IF (KZWOPT.GT.0) THEN
                  WRITE(LUPRI,'(/2A)')
     *            ' FOCK MATRIX CONTRIBUTION TO DIAGONAL ORBITAL PART',
     *            ' OF E[2]'
                  IF (AVDIA) WRITE (LUPRI,'(/2A)')
     *            ' (FOCK TYPE DECOUPLING OF TWO ELECTRON DENSITY',
     *            ' MATRIX USED)'
                  CALL PRWOP(WRK(KEODIA),LUPRI)
                  IF (IPRRSP .GT. 120)
     *             CALL OUTPUT(WRK(KEODIA),1,NWOPT,1,1,NWOPT,1,-1,LUPRI)
               END IF
            END IF
         END IF
C
C *** now WRITE orbital diagonal of S[2]-matrix
C
         IF (KZWOPT .GT. 0) THEN
            CALL  WRITT(LURSP4,KZWOPT,WRK(KSODIA))
         END IF
C
         IF (IPRRSP.GT.100) THEN
            IF (KZWOPT.GT.0) THEN
               WRITE(LUPRI,'(/A)') ' DIAGONAL ORBITAL PART OF S[2]'
               CALL OUTPUT(WRK(KSODIA),1,NWOPT,1,1,NWOPT,1,-1,LUPRI)
            END IF
         END IF
C
C WORK SPACE ALLOCATION
C
         IF (KZCONF.GT.0) THEN
            CALL PHPINI(LPHPMX,KZCONF,KZWOPT,MAXPHP)
            KDIAE  = 1
            KTOT   = KDIAE  + LPHPMX
            LTOT   = LWRK   - KTOT
            IF (LTOT.LT.0) CALL ERRWRK('RSPVAR 2',KTOT-1,LWRK)
C DEFINE PARAMETERS FOR PHPGET
            IPWAY  = 1
            ECORE  = -EACTIV
            IF (MAXPHP .GT. 0) THEN
               CALL RSPEDG(WRK(KDIAE))
            END IF
            IPRPHP = IPRRSP / 10
            CALL PHPGET(KSYMST,KZCONF,XINDX,FCAC,H2AC,WRK(KDIAE),
     *            ECORE,IPWAY,IPRPHP,WRK(KTOT),LTOT)
            CALL RSPEDG(WRK(KDIAE))
C WRITE H(0)inverse INFORMATION ON DISK
            CALL PHPDSK(WRK(KDIAE),LPHPMX,.TRUE.)
         ELSE
            LPHPMX = KZVAR
C LPHPMX IS USED LATER AND ASSUMED ALWAYS TO HAVE LENGTH KZVAR
         END IF
C
C *** END OF RSPVAR
C
      CALL QEXIT('RSPVAR')
      RETURN
      END
C  /* Deck rspmc */
      SUBROUTINE RSPMC(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
C
C 28 JAN 87 BASED ON READMC
C Revised May 2001 /hjaaj
C
C Purpose:
C   Read that MC information written by WRSIFC in SIRIUS which
C   is needed for RESPONSE calculation.
C
C   We assume that either the SIRIUS module or the SETSIR
C   subroutine has been called before RSPMC.
C
      use pelib_interface, only: use_pelib
#if defined (HAS_PCMSOLVER)
      use pcm_config, only: pcm_configuration, pcm_cfg
#endif
#include "implicit.h"
#include "infdim.h"
C
      DIMENSION CMO(*),UDV(*),PV(*),XINDX(LCINDX),
     *          FOCK(*),FC(*),FV(*),FCAC(*),H2AC(*),WRK(*)
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0 )
      LOGICAL   HSFOCK
C
C Used from common blocks:
C SCBRHF : IOPRHF
C QRINF  : MZCONF(1)
C
#include "maxorb.h"
#include "priunit.h"
#include "infinp.h"
#include "pcmlog.h"
#include "inftap.h"
#include "infopt.h"
#include "infpri.h"
#include "inforb.h"
#include "infvar.h"
#include "infsop.h"
#include "inflin.h"
#include "qrinf.h"
#include "scbrhf.h"
#include "infrsp.h"
#include "wrkrsp.h"
C--------------------
C CBN+JK 03.01.06
C--------------------
#include "gnrinf.h"
#include "qm3.h"
#include "mxcent.h"
C--------------------
C CBN+JK 03.01.06
C--------------------
C
      CALL QENTER('RSPMC ')
C
C     ====================================================
C     Initial checks based on info in Sirius common blocks
C     (we assume that either the SIRIUS module or the SETSIR
C      subroutine has been called before RSPMC).
C
C     Find out if TDHF or not :
C
      HSFOCK = HSROHF .OR. (DODFT .AND. NASHT .GT. 0)
      ! Open-shell DFT response only works with HSROHF code
      TDHF = NASHT.EQ.0 .OR. (NASHT.EQ.1 .AND. NACTEL.EQ.1) .OR. HSROHF
      ! Closed-shell HF or DFT, or one open shell (.OPEN SHELL keyword in SIRIUS)
      ! or high-spin ROHF (.SINGLY OCCUPIED keyword in SIRIUS)
      HSFOCK = HSFOCK .OR. (DIRFCK .AND. TDHF .AND. NASHT .GT. 0)
      ! Some of the open-shell direct response only works with HSROHF code
C
C     If any triplet operators, then we must use determinants
C     in RESPONSE also for the reference symmetry, even if
C     CSF's were used in SIRIUS. /hjaaj May 2001
C
      MCONRF = NCONRF
      IF (TRPFLG .AND. .NOT.FLAG(27) .AND. .NOT.TDHF) THEN
         FLAG(27) = .TRUE.
C        ... call setsir to reset Sirius common blocks to determinants
         CALL SETSIR(WRK,LWRK)
      END IF
      NCONRF = MAX(1,NCONRF) ! SETSIR(abaset.F) sets NCONRF=0 for Hartree-Fock
                             ! in RSP - especially for RSPESR - we want NCREF=1 for HSROHF
C
C     ====================================================
C     Read and check information on LUSIFC
C
      REWIND LUSIFC
      IF (RSPCI) THEN
         CALL MOLLAB('CIRESPON',LUSIFC,LUPRI)
      ELSE
         CALL MOLLAB(LBSIFC,LUSIFC,LUPRI)
      END IF
      READ (LUSIFC) POTNUC,EMY,EACTIV,EMCSCF,
     &              ISTATE,ISPIN,NACTEL,LSYM,MS2
      NDUM = 1 + 64 + 6*8 + 4
      READ (LUSIFC) MISHT,MASHT, MOCCT, MORBT, MBAST, MCONF,MWOPT,MWOPH,
     &              NCDETS,MCMOT,MNASHX,MNASHY,MNORBT,M2ORBT,
     &              (IDUM,I=1,NDUM),MCTYPE_IFC
      HSROHF = MCTYPE_IFC .EQ. -1
C     READ (LUSIFC) NISHT,NASHT, NOCCT, NORBT, NBAST, NCONF,NWOPT,NWOPH,
C    &              NCDETS,NCMOT,NNASHX,NNASHY,NNORBT,N2ORBT
C    &              NSYM,MULD2H, NRHF,NFRO,
C    &              NISH,NASH,NORB,NBAS,
C    &              NELMN1, NELMX1, NELMN3, NELMX3, MCTYPE_IFC,
C    &              NAS1, NAS2, NAS3
C
C     Transfer information
C
      NCREF  = NCONRF
      IREFSY = LSYMRF
      DETERM = FLAG(27)
C
C     Define MZCONF(1)
C     ... MZCONF(1) is used in qr and cr for cref vector and
C         it will not be defined in qr/cr modules if no
C         perturbation operators of symmetry 1. /hjaaj Mar 2001
C
      IF (TDHF) THEN
         MZCONF(1) = 0
      ELSE
         MZCONF(1) = NCREF
      END IF
C
      IF (IPRRSP .GT. 0) THEN
      IF (NASHT.EQ.0) THEN
         WRITE (LUPRI,'(//,(A,F25.15))')
     *   '   SCF energy         :',EMCSCF,
     *   ' -- inactive part     :',EMY,
     *   ' -- nuclear repulsion :',POTNUC
      ELSE IF (NASHT.EQ.1 .AND. NACTEL.EQ.1) THEN
         WRITE (LUPRI,'(//,(A,F25.15))')
     *   ' open shell SCF energy:',EMCSCF,
     *   ' -- inactive part     :',EMY,
     *   ' --   active part     :',EACTIV,
     *   ' -- nuclear repulsion :',POTNUC
      ELSE IF (HSROHF) THEN
         WRITE (LUPRI,'(//,(A,F25.15))')
     *   ' high-spin open-shell SCF energy:',EMCSCF,
     *   ' -- inactive part               :',EMY,
     *   ' --   active part               :',EACTIV,
     *   ' -- nuclear repulsion           :',POTNUC
      ELSE
         WRITE (LUPRI,'(//,(A,F25.15))')
     *   ' MCSCF energy         :',EMCSCF,
     *   ' -- inactive part     :',EMY,
     *   ' --   active part     :',EACTIV,
     *   ' -- nuclear repulsion :',POTNUC
      END IF
      END IF
C
      IF ( IPRRSP.GT.10 ) WRITE(LUPRI,'(/A,6I6)')
     *      ' NISHT,NASHT,NORBT,NBAST,NCREF,IREFSY',
     *        NISHT,NASHT,NORBT,NBAST,NCREF,IREFSY
      IF (RSPSUP) WRITE (LUPRI,'(/A)')
     &   '@ Super symmetry will be used in perturbation symmetry 1'
C
C CHECK INFORMATION FROM LUSIFC WITH CORRESPONDING INFORMATION FROM
C SIRINP
C
      NUMERR = 0
      IF ( MISHT .NE. NISHT ) THEN
         NUMERR = NUMERR + 1
         WRITE(LUPRI,'(A/A,I5,A,I5)')
     *   '@ *** FATAL ERROR IN RSPMC ***',
     *   '@ TOTAL NUMBER OF INACTIVE ORBITALS FROM SIRINP:',NISHT,
     *   '@ FROM SIRIFC:',MISHT
      END IF
      IF (TDHF .AND. NASHT .EQ. 0 .AND. EACTIV .NE. D0) THEN
         NUMERR = NUMERR + 1
         WRITE(LUPRI,'(A/A,F25.15)')
     *   '@ *** FATAL ERROR IN RSPMC ***',
     *   '@ TDHF calculation, but EACTIV .ne. 0; EACTIV is',EACTIV
      END IF
      IF ( MASHT .NE. NASHT ) THEN
         NUMERR = NUMERR + 1
         WRITE(LUPRI,'(A/A,I5,/A,I5)')
     *   '@ *** FATAL ERROR IN RSPMC ***',
     *   '@ TOTAL NUMBER OF ACTIVE ORBITALS FROM SIRINP:',NASHT,
     *   '@ FROM SIRIFC:',MASHT
      END IF
      IF ( MOCCT .NE. NOCCT ) THEN
         NUMERR = NUMERR + 1
         WRITE(LUPRI,'(A/A,I5,/A,I5)')
     *   '@ *** FATAL ERROR IN RSPMC ***',
     *   '@ TOTAL NUMBER OF OCCUPIED ORBITALS FROM SIRINP:',NOCCT,
     *   '@ FROM SIRIFC:',MOCCT
      END IF
      IF ( MORBT .NE. NORBT ) THEN
         NUMERR = NUMERR + 1
         WRITE(LUPRI,'(A/A,I5,/A,I5)')
     *   '@ *** FATAL ERROR IN RSPMC ***',
     *   '@ TOTAL NUMBER OF  ORBITALS FROM SIRINP:',NORBT,
     *   '@ FROM SIRIFC ',MORBT
      END IF
      IF ( MBAST .NE. NBAST ) THEN
         NUMERR = NUMERR + 1
         WRITE(LUPRI,'(A/A,I5,/A,I5)')
     *   '@ *** FATAL ERROR IN RSPMC ***',
     *   '@ TOTAL NUMBER OF ATOMIC ORBITALS FROM SIRINP:',NBAST,
     *   '@ FROM SIRIFC ',MBAST
      END IF
      IF ( MCONRF .NE. MCONF .AND. MCONF .GT. 1 .AND. MCONRF .GT. 1 )
     &   THEN
         NUMERR = NUMERR + 1
         WRITE(LUPRI,'(A/A,I5/A,I5,A,L5)')
     *  '@ *** FATAL ERROR IN RSPMC ***',
     *  '@ TOTAL NUMBER OF REFERENCE CONFIGURATIONS FROM SIRINP:',
     *  MCONRF,'@ AND FROM SIRIFC',MCONF,';  TDHF flag is',TDHF
      END IF
C     NWOPT test disabled because RSPMC may be called from ABACUS
C     with NCONF,NWOPT not of symmetry 1 /25-Nov-1992 PJ+HJAaJ
C     30-Nov-93 hjaaj: test reactivated with JWOPSY = 1 check
C                      (this should also catch SUPSYM inconsistencies)
      IF ( JWOPSY .EQ. 1 .AND. NWOPT .NE. MWOPT ) THEN
         IF ( NWOPT .LT. MWOPT) THEN
            ! some rotations may have been proned in RSPVAR
            LUINDX = -1
            CALL GPOPEN(LUINDX,'LUINDF','UNKNOWN',' ','UNFORMATTED',
     &                  IDUMMY,.FALSE.)
            REWIND LUINDX
            CALL MOLLAB('EXOPSYM1',LUINDX,LUPRI)
            READ (LUINDX) NWOPT_orig
            CALL GPCLOSE(LUINDX,'KEEP')
         ELSE
            NWOPT_orig = NWOPT
         END IF
         IF (NWOPT_orig .NE. MWOPT) THEN
            NUMERR = NUMERR + 1
            WRITE(LUPRI,'(A/A,I5,/A,I5)')
     *      '@ *** FATAL ERROR IN RSPMC ***',
     *      '@ TOTAL NUMBER OF sym 1 orb. rot. FROM SIRINP:',NWOPT,
     *      '@ FROM SIRIFC ',MWOPT
         END IF
      END IF
      IF (NUMERR .GT. 0)
     *   CALL QUIT('RSPMC ERROR, SIRIFC file inconsistent '//
     &             'with SIRIUS input.')
C
C READ IN MATRICES NEEDED IN RESPONSE CALCULATION
C
      CALL READT (LUSIFC,NCMOT,CMO)
      READ (LUSIFC)
      KDV   = 1 + N2BASX ! DV matrix may be needed later
      KFREE = KDV + NNASHX
      IF (NASHT .GT. 0) THEN
C        read and unpack DV
         CALL READT (LUSIFC,NNASHX,WRK(KDV))
         CALL DSPTSI(NASHT,WRK(KDV),UDV)
         IF (.NOT.RSPCI) THEN
            CALL READT (LUSIFC,N2ORBT,FOCK)
            IF (HSROHF) THEN
               READ (LUSIFC)
            ELSE
               CALL READT (LUSIFC,NNASHY,PV)
            END IF
            CALL READT (LUSIFC,NNORBT,FC)
            CALL READT (LUSIFC,NNORBT,FV)
         END IF
         CALL READT (LUSIFC,NNASHX,FCAC)
         CALL READT (LUSIFC,(NNASHX*NNASHX),H2AC)
      ELSE
         READ (LUSIFC)
         CALL READT (LUSIFC,N2ORBT,FOCK)
         READ (LUSIFC)
         CALL READT (LUSIFC,NNORBT,FC)
         CALL DZERO(FV,NNORBT)
      END IF
C
C     If a solvent-response calculation:
C     (i) read solvent t matrix (ii) subtract it from FC
C     (because the response program requires a FC matrix
C      without solvent contributions).
C
      IF (QM3) THEN
        IF (OLDTG .OR. (.NOT. LOSPC)) THEN
          CALL MOLLAB('SOLVTMAT',LUSIFC,LUPRI)
          CALL READT (LUSIFC,NNORBT,WRK)
          DO 21 I = 1,NNORBT
             FC(I) = FC(I) - WRK(I)
 21       CONTINUE
        END IF
      END IF
      IF (FLAG(16)) THEN
        CALL MOLLAB('SOLVTMAT',LUSIFC,LUPRI)
        CALL READT (LUSIFC,NNORBT,WRK)
        DO 20 I = 1,NNORBT
           FC(I) = FC(I) - WRK(I)
 20     CONTINUE
Cjje
      ELSE IF(PCM .and. (.not. SOPRPA)) then
!      ELSEIF(PCM) THEN
        CALL MOLLAB('PCMJXMAT',LUSIFC,LUPRI)
        CALL READT (LUSIFC,NNORBT,WRK)
        DO I = 1,NNORBT
           FC(I) = FC(I) + WRK(I)
        ENDDO
#if defined (HAS_PCMSOLVER)
      else if (pcm_cfg%do_pcm) then
       CALL MOLLAB('EXTPCMMAT',LUSIFC,LUPRI)
       CALL READT (LUSIFC,NNORBT,WRK)
       DO I = 1,NNORBT
          FC(I) = FC(I) + WRK(I)
       ENDDO
#endif
      END IF
Cjje
C
      IF ((NCONRF .GT. 1) .AND. QMMM ) THEN
         write(lupri,*) 'MCSCF detected (NCONRF > 1) in QMMM calc.'
         write(lupri,*) 'NCONRF...........................', NCONRF
         write(lupri,*) '=> Subtracting QMMM Fock matrix (PCM-like)'
         CALL MOLLAB('QMMMTMAT',LUSIFC,LUPRI)
         CALL READT (LUSIFC,NNORBT,WRK)
         DO I = 1,NNORBT
             FC(I) = FC(I) - WRK(I)
         ENDDO
      END IF
C
      IF (USE_PELIB() .AND. .NOT. TDHF) THEN
         CALL MOLLAB('PEFMAT  ', LUSIFC, LUPRI)
         CALL READT (LUSIFC, NNORBX, WRK)
         FC(1:NNORBX) = FC(1:NNORBX) - WRK(1:NNORBX)
      ENDIF
C
      IF (SOPPA) THEN
C
         IF (CCPPA) THEN
            CALL MOLLAB('CCSDINFO',LUSIFC,LUPRI)
         ELSE
            CALL MOLLAB('MP2INFO ',LUSIFC,LUPRI)
         ENDIF
C
C        reads the MP2 or CCSD correlation coefficients into PV
C
         CALL READT (LUSIFC,LPVMAT,PV)
C
         IF (IPRRSP.GT.20 .AND. .NOT.TRPFLG) THEN
            IF (CCPPA) THEN
              WRITE(LUPRI,'(/A)')
     &              ' RSPMC : CCSD correlation coefficients'
            ELSE
              WRITE(LUPRI,'(/A)')' RSPMC : MP2 correlation coefficients'
            ENDIF
            CALL OUTPUT(PV,1,1,1,LPVMAT,1,LPVMAT,-1,LUPRI)
         END IF
C
         IF (TRPFLG)
     &      CALL TRPKAP(PV(1),PV(LPVMAT+1),XINDX(KIADR1),WRK(KFREE))
C
         IF (IPRRSP.GT.20 .AND. TRPFLG) THEN
            IF (CCPPA) THEN
               WRITE(LUPRI,'(/A)')
     &      ' RSPMC : CCSD singlet and triplet correlation coefficients'
            ELSE
               WRITE(LUPRI,'(/A)')
     &      ' RSPMC : MP2 singlet and triplet correlation coefficients'
            ENDIF
            CALL OUTPUT(PV,1,2,1,LPVMAT,2,LPVMAT,-1,LUPRI)
         END IF
C
         CALL READT (LUSIFC,NACTT*NACTT,UDV)
C
         IF (IPRRSP.GT.20) THEN
            IF (CCPPA) THEN
               WRITE(LUPRI,'(/A)')' RSPMC : CCSD "second order" density'
            ELSE
               WRITE(LUPRI,'(/A)')' RSPMC : MP2 second order density'
            END IF
            CALL OUTPUT(UDV,1,NACTT,1,NACTT,NACTT,NACTT,-1,LUPRI)
         END IF
C
C        UDV contains the MP2 one-density. Remove the diagonal
C        contribution from the zeroth order. (Added in MP2FAC)
C
         CALL SOPUDV(UDV)
      END IF
C
      IF (HSFOCK .AND. .NOT.HSROHF) THEN
         ! We want to use HSROHF in response, but SIRIUS did not use HSROHF
         ! FC, FV are then FC, FV; and not FC+FV+FV(exch), -FV(exch) as expected for HSROHF
         HSROHF = .TRUE.
         WRITE(LUPRI,'(/A)')
     &   'INFO: switching to HSROHF Fock matrices for RESPONS'
         CALL FCKMAT((NASHT.EQ.0),WRK(KDV),CMO,EMY,FC,FV,
     &      WRK(KFREE),LWRK-KFREE)
         CALL DAXPY(NNORBT,-1.0D0,FV,1,FC,1)
      END IF
C
      IF (IPRRSP.GT.20) THEN
         DO 60 ISYM = 1,NSYM
            WRITE(LUPRI,'(/A,I5)') ' **** RSPMC **** '//
     &      'MO coefficients read from SIRIFC for symmetry',ISYM
            CALL OUTPUT(CMO(ICMO(ISYM)+1),1,NBAS(ISYM),1,NORB(ISYM),
     &                  NBAS(ISYM),NORB(ISYM),-1,LUPRI)
 60      CONTINUE
      END IF
C
      IF ( IPRRSP.GT.110 .AND. .NOT. TDHF) THEN
         CALL GETREF(WRK,NCONRF)
         WRITE(LUPRI,'(/A)')
     &      ' **** RSPMC **** Reference configuration coefficients'
         CALL OUTPUT(WRK,1,NCONRF,1,1,NCONRF,1,1,LUPRI)
      END IF
C
C *** End of RSPMC
C
      CALL QEXIT('RSPMC ')
      RETURN
      END
C  /* Deck getref */
      SUBROUTINE GETREF(CREF,LCREF)
C
C 28 JAN 87 BASED ON READMC(RSPMC)
C l.r. Sep 2004 hjaaj
C
C Purpose:
C  GET MC REFERENCE STATE from SIRIFC
C
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION CREF(LCREF)
      LOGICAL   FNDLAB, CIRESPON
C
#include "inftap.h"
C
      IF (LCREF .EQ. 1) THEN
         CREF(1) = 1.0D0
      ELSE IF (LCREF .EQ. 0) THEN
C        ... Hartree-Fock case ....
      ELSE
         REWIND LUSIFC
         IF ( FNDLAB(LBSIFC,LUSIFC) ) THEN
            CIRESPON = .FALSE.
         ELSE
            CIRESPON = .TRUE.
            REWIND LUSIFC
            CALL MOLLAB('CIRESPON',LUSIFC,LUPRI)
         END IF
         READ (LUSIFC)
         READ (LUSIFC) MISHT,MASHT,MOCCT,MORBT,MBAST,
     *                 MCONF,MWOPT,MWOPH,MDETS
         IF (LCREF .EQ. MCONF) THEN
C
            READ (LUSIFC)
         ELSE IF (LCREF .EQ. MDETS) THEN
            CALL MOLLAB('CREFDETS',LUSIFC,LUPRI)
         ELSE
            WRITE(LUPRI,'(/A/3I20)')
     &         'GETREF error, LCREF .ne. NCONF/NCDETS from SIRIFC:',
     &         LCREF,MCONF,MDETS
            CALL QUIT('GETREF error, LCREF .ne. NCONF')
         END IF
         CALL READT (LUSIFC,LCREF,CREF)
      END IF
C
C *** End of GETREF
C
      RETURN
      END
      SUBROUTINE LINEAR(WORD)
C
C Linear response input
C
#include "implicit.h"
#include "priunit.h"
C
#include "infpri.h"
#include "infrsp.h"
C
      PARAMETER ( NTABLE = 1 )
      CHARACTER WORD*7, STRING*7, TABLE(NTABLE)*7
      PARAMETER ( STRING = '*LINEAR' )
      DATA TABLE /'.SINGLE'/
C
  100 READ (LUCMD, '(A7)',END=2100) WORD
      CALL UPCASE(WORD)
      IF (WORD(1:1).EQ.'*') GO TO 11
      DO I = 1, NTABLE
         IF (WORD .EQ. TABLE(I)) THEN
            GO TO (1), I
         END IF
      END DO
      GO TO 100
C
C Without keyword a normal alpha calculation is assumed
C
 11   REWIND (LUCMD,IOSTAT=IOS)
 111  READ (LUCMD, '(A7)',END=2100) WORD
      CALL UPCASE(WORD)
      IF (WORD .NE. STRING) GO TO 111
      CALL LRINP(WORD)
      RETURN
C
C A single residue of the linear response function is calculated
C
 1    REWIND (LUCMD,IOSTAT=IOS)
 110  READ (LUCMD, '(A7)',END=2100) WORD
      CALL UPCASE(WORD)
      IF (WORD .NE. STRING) GO TO 110
      CALL PPINP(WORD)
      RETURN
C
 2100 CALL QUIT('Final "**<something>" not found in DALTON input')
C
      RETURN
      END

      SUBROUTINE QUADRA(WORD)
C
C Quadratic response input
C
#include "implicit.h"
#include "priunit.h"
C
#include "infpri.h"
#include "infrsp.h"
#include "rspprp.h"
#include "infpp.h"
C
      PARAMETER ( NTABLE = 13 )                                         !MK
      CHARACTER WORD*7, STRING*7, TABLE(NTABLE)*7
      PARAMETER ( STRING = '*QUADRA' )
      DATA TABLE /'.SINGLE', '.DOUBLE', '.PHOSPH', '.TWO-PH',
     &            '.MNFPHO', '.TPCD  ', '.ESA   ', '.ECPHOS',           !MK
     &            '.PHOSPV', '.CPPHOL', '.CPPHOV', '.CPPHMF',           !MK
     &            '.CPPHEC'/                                            !MK
C
  100 READ (LUCMD, '(A7)',END=2100) WORD
      CALL UPCASE(WORD)
      IF (WORD(1:1).EQ.'*') GO TO 11
      DO I = 1, NTABLE
         IF (WORD .EQ. TABLE(I)) THEN
            GO TO (1,2,3,4,5,4,6,8,3,3,3,5,8), I                        !MK
         END IF
      END DO
      GO TO 100
C
C Without keyword a normal beta calculation is assumed
C
 11   REWIND (LUCMD,IOSTAT=IOS)
 111  READ (LUCMD, '(A7)',END=2100) WORD
      CALL UPCASE(WORD)
      IF (WORD .NE. STRING) GO TO 111
      CALL HYPINP(WORD)
      RETURN
C
C A single residue of the quadratic response function is calculated
C
 1    REWIND (LUCMD,IOSTAT=IOS)
 110  READ (LUCMD, '(A7)',END=2100) WORD
      CALL UPCASE(WORD)
      IF (WORD .NE. STRING) GO TO 110
      CALL SMOINP(WORD)
      RETURN
C
C A double residue of the quadratic response function is calculated
C
 2    REWIND (LUCMD,IOSTAT=IOS)
 120  READ (LUCMD, '(A7)',END=2100) WORD
      CALL UPCASE(WORD)
      IF (WORD .NE. STRING) GO TO 120
      EXMOM = .TRUE.
      CALL PPINP(WORD)
      RETURN
C
C Spin orbit calculation
C
 3    REWIND (LUCMD,IOSTAT=IOS)
 130  READ (LUCMD, '(A7)',END=2100) WORD
      CALL UPCASE(WORD)
      IF (WORD .NE. STRING) GO TO 130
      CALL SMOINP(WORD)
      RETURN
C
C Two-photon absorption calculation
C
 4    REWIND (LUCMD,IOSTAT=IOS)
 140  READ (LUCMD, '(A7)',END=2100) WORD
      CALL UPCASE(WORD)
      IF (WORD .NE. STRING) GO TO 140
      CALL SMOINP(WORD)
      RETURN
C
C     Phosphorescence calculation using mean-field integrals
C
 5    REWIND (LUCMD,IOSTAT=IOS)
 150  READ (LUCMD, '(A7)',END=2100) WORD
      CALL UPCASE(WORD)
      IF (WORD .NE. STRING) GO TO 150
      CALL SMOINP(WORD)
      RETURN
C
C     Phosphorescence calculation using effective charge integrals
C
 8    REWIND (LUCMD,IOSTAT=IOS)
 180  READ (LUCMD, '(A7)',END=2100) WORD
      CALL UPCASE(WORD)
      IF (WORD .NE. STRING) GO TO 180
      CALL SMOINP(WORD)
      RETURN
C
C Excited state absorption calculation
C
 6    REWIND (LUCMD,IOSTAT=IOS)
 160  READ (LUCMD, '(A7)',END=2100) WORD
      IF (WORD .NE. STRING) GO TO 160
      EXMOM = .TRUE.
      CALL PPINP(WORD)
      RETURN
C
 2100 CALL QUIT('**QUADRA not found in DALTON input')
C
      RETURN
      END
      SUBROUTINE CUBIC(WORD)
C
C Cubic response input
C
#include "implicit.h"
#include "priunit.h"
C
#include "infpri.h"
#include "infrsp.h"
C
      PARAMETER ( NTABLE = 3 )
      CHARACTER WORD*7, STRING*7, TABLE(NTABLE)*7
      PARAMETER ( STRING = '*CUBIC ' )
      DATA TABLE /'.SINGLE', '.DOUBLE','.THREE-'/
C
  100 READ (LUCMD, '(A7)',END=2100) WORD
      CALL UPCASE(WORD)
      IF (WORD(1:1).EQ.'*') GO TO 11
      DO I = 1, NTABLE
         IF (WORD .EQ. TABLE(I)) THEN
            GO TO (1,2,3), I
         END IF
      END DO
      GO TO 100
C
C Without keyword a normal gamma calculation is assumed
C
 11   REWIND (LUCMD,IOSTAT=IOS)
 111  READ (LUCMD, '(A7)',END=2100) WORD
      CALL UPCASE(WORD)
      IF (WORD .NE. STRING) GO TO 111
      CALL CRINP(WORD)
      RETURN
C
C A single residue of the cubic response function is calculated
C
 1    REWIND (LUCMD,IOSTAT=IOS)
 110  READ (LUCMD, '(A7)',END=2100) WORD
      CALL UPCASE(WORD)
      IF (WORD .NE. STRING) GO TO 110
      CALL TMOINP(WORD)
      RETURN
C
C A double residue of the cubic response function is calculated
C
 2    REWIND (LUCMD,IOSTAT=IOS)
 120  READ (LUCMD, '(A7)',END=2100) WORD
      CALL UPCASE(WORD)
      IF (WORD .NE. STRING) GO TO 120
      CALL TPAINP(WORD)
      RETURN
C
C Three-photon absorption calculation
C
 3    REWIND (LUCMD,IOSTAT=IOS)
 130  READ (LUCMD, '(A7)',END=2100) WORD
      CALL UPCASE(WORD)
      IF (WORD .NE. STRING) GO TO 130
      CALL TMOINP(WORD)
      RETURN
C
 2100 CALL QUIT('**END OF not found in DALTON input after *CUBIC ')
C
      RETURN
      END

      SUBROUTINE ESGINP(WORD)
C
C Excited State Gradient input
C
#include "implicit.h"
#include "priunit.h"
C
#include "infpri.h"
#include "rspprp.h"
#include "esg.h"
C
      PARAMETER ( NTABLE = 1 )
      CHARACTER WORD*7, STRING*7
      PARAMETER ( STRING = '*ESG' )
      ESG = .TRUE.
C
C A single residue of the linear response function is calculated
C
      REWIND (LUCMD,IOSTAT=IOS)
 110  READ (LUCMD, '(A7)',END=2100) WORD
      CALL UPCASE(WORD)
      IF (WORD .NE. STRING) GO TO 110
      CALL ESGGRINP(WORD)
      RETURN
C
 2100 CALL QUIT('ESGINP: "*ESG" not found in DALTON input')
C
      RETURN
      END


      SUBROUTINE ESGGRINP(WORD)
      use pelib_interface, only: use_pelib
#include "implicit.h"
      LOGICAL NEWDEF, SPNDEF
      PARAMETER ( NTABLE = 8 )
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER*8 LABEL
C
#include "priunit.h"
#include "inforb.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "infpri.h"
#include "rspprp.h"
#include "infpp.h"
#include "infspi.h"
#include "infhso.h"
#include "inflr.h"
#include "inflin.h"
#include "maxorb.h"
#include "infinp.h"
#include "inftap.h"
#include "esg.h"
#include "gnrinf.h"
      DATA TABLE /'.ROOTS', '.THCPP ', '.MAX IT', '.MAXITO', '.THCLR ',
     *            '.RESTPP','.PRINT ', '.STATE '/
C
      DATA SPNDEF/.FALSE./
C
C READ IN  INPUT
C
      NEWDEF = (WORD.EQ.'*ESG')
      ICHANG = 0
      NPPTOT = 0
      IF (NEWDEF) THEN
         DO I=1,NSYM
            NPPCNV(I) = 0
         END DO
         NPPTOT = NSYM
         WORD1  = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') GO TO 100
            IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO I=1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8), I
                  END IF
               END DO
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') '@KEYWORD "',WORD,
     *            '" NOT RECOGNIZED IN ESG.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT(' ILLEGAL KEYWORD IN ESG INPUT ')
               GO TO 100
C
 1             CONTINUE
               NPPTOT = 0
               READ (LUCMD,*) (NPPCNV(MULD2H(J,LSYMRF)),J=1,NSYM)
               DO J=1,NSYM
                  NPPTOT = NPPTOT + NPPCNV(J)
               END DO
               GO TO 100
 2             CONTINUE
               READ (LUCMD,*) THCPP
               GO TO 100
 3             CONTINUE
               READ (LUCMD,*) MAXITP
               GO TO 100
 4             CONTINUE
               READ (LUCMD,*) MAXITO
               GO TO 100
 5             CONTINUE
               READ (LUCMD,*) THCLR
               GO TO 100
 6             CONTINUE
               RESTPP = .TRUE.
               GO TO 100
 7             CONTINUE
               READ (LUCMD,*) IPRPP
               GO TO 100
 8             CONTINUE
               READ (LUCMD,*) IESG, ISYME
               IF (NPPCNV(MULD2H(ISYME,LSYMRF)).LT.IESG) THEN
                  NPPTOT = NPPTOT - NPPCNV(MULD2H(ISYME,LSYMRF)) + IESG
                  NPPCNV(MULD2H(ISYME,LSYMRF)) = IESG
               END IF
               GOTO 100

            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') '@PROMPT "',WORD,
     *            '" NOT RECOGNIZED IN ESG INPUT.'
               CALL QUIT(' ILLEGAL PROMPT IN ESG INPUT ')
            END IF
      END IF
  300 CONTINUE
      IF (THR_REDFAC .GT. 0.0D0) THEN
         ICHANG = ICHANG + 1
         WRITE (LUPRI,'(3A,1P,D10.2)') '@ INFO ',WORD1,
     &   ' thresholds multiplied with general factor',THR_REDFAC
         THCPP = THCPP*THR_REDFAC
         THCLR = THCLR*THR_REDFAC
      END IF
      IF (ICHANG .GT. 0) THEN
      IF (NPPTOT .EQ. 0) THEN
         NINFO = NINFO + 1
         WRITE (LUPRI,'(/3A)') '@INFO: ',WORD1,
     *      ' input ignored because no roots requested.'
      ELSE
         CALL HEADER(' Excited State Gradient calculation',0)

         WRITE(LUPRI,'(A,I5)')
     *        'Number of roots calculated in response   : ', NPPTOT
         WRITE(LUPRI,'(A,I5)')
     *        'Root selected for excited state gradient : ', IESG
         WRITE(LUPRI,'(A,I5,/)')
     *        'Symmetry of selected root                : ', ISYME

         IF (FLAG(16) .AND. .NOT.INERSI) THEN
            WRITE(LUPRI,'(A,L1,/)')
     *      ' Equilibrium solvent model requested            : SOLVNT ='
     *      ,FLAG(16)
            WRITE(LUPRI,'(A,F8.4)')
     *      ' Dielectric constant                            : EPSOL  ='
     *      ,EPSOL
         END IF
         IF (FLAG(16) .AND. INERSI) THEN
            WRITE(LUPRI,'(A,L1,/)')
     *      ' Non-equilibrium solvent model requested        : INERSI ='
     *      ,INERSI
            WRITE(LUPRI,'(A,F8.4)')
     *      ' Static dielectric constant                     : EPSTAT ='
     *      ,EPSTAT
            WRITE(LUPRI,'(A,F8.4)')
     *      ' Optical dielectric constant                    : EPSOL  ='
     *      ,EPSOL
         END IF
C
         IF (QM3 .OR. QMMM) THEN
           WRITE(LUPRI,'(A)') ' QM/MM response calculation'
         ENDIF

         IF (USE_PELIB()) THEN
            WRITE(LUPRI,'(A)') ' Environment is modeled using'//
     &                  ' the polarizable embedding scheme (PE library)'
         END IF

         WRITE(LUPRI,'(A,I4)')
     *      ' Print level for eigenval.eqs.                  : IPRPP  ='
     *       ,IPRPP
         WRITE(LUPRI,'(A,I4)')
     *      ' Maximum number of iterations for eigenval.eqs. : MAXITP ='
     *       ,MAXITP
         WRITE(LUPRI,'(A,1P,D10.3)')
     *      ' Threshold for convergence of eigenval.eqs.     : THCPP  ='
     *      ,THCPP
         WRITE(LUPRI,'(A,I4)')
     *      ' Maximum number of linear response iterations   : MAXITL ='
     *       ,MAXITL
         WRITE(LUPRI,'(A,1P,D10.3)')
     *      ' Threshold for linear response convergence      : THCLR  ='
     *      ,THCLR
         WRITE(LUPRI,'(A,I4)')
     *      ' Maximum iterations in optimal orbital algorithm: MAXITO ='
     *      ,MAXITO

         IF (RESTPP) WRITE(LUPRI,'(/A,L1)')
     *      ' USE TRIAL VECTORS AVAILABLE FOR RESTART. RESTPP =',RESTPP
         WRITE (LUPRI,'(/A,I5)')' Spin of operator A , ISPINA=',ISPINA
         WRITE (LUPRI,'( A,I5)')' Spin of operator B ,'//
     *        ' (Excitation energy) ISPINB=',ISPINB
         WRITE (LUPRI,'( A,I5)') ' Spin of operator C ,'//
     *        ' (Excitation energy) ISPINC=',ISPINC
      END IF
      END IF
C
C  *** END OF ESG
C
      RETURN
      END

      SUBROUTINE ASSERT_TRIPLET_RESPONSE_WITH_DETERMINANTS()
      IMPLICIT NONE
#include <maxorb.h>
#include <infinp.h>
#include <infrsp.h>
#include <priunit.h>
      CALL QENTER('ASSERT_TRIPLET_RESPONSE_WITH_DETERMINANTS')
      IF (TRPFLG .AND. .NOT.FLAG(27) .AND. .NOT.TDHF) THEN
         WRITE(LUPRI, *)
     &      'Please use **WAVE FUNCTION *OPTIMIZATION .DETERMINANTS'
         CALL QUIT(
     &      'A TRIPLET RESPONSE CALCULATION REQUIRES DETERMINANTS'
     &      )
      END IF
      CALL QEXIT('ASSERT_TRIPLET_RESPONSE_WITH_DETERMINANTS')
      END
