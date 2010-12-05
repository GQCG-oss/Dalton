#include <single.h>
C
C  /* Deck rp_init */
      SUBROUTINE RP_INIT(WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, Januar 1998
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Initialize some pointers and sort AO-integrals
C     into distributions. Based on CC_INP and ???
C
#include <implicit.h>
#include <priunit.h>
C
#include <ccsdinp.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccorb.h>
CSPAS:10/11-2003: replaced with new version from CC_INP
#include <ccsections.h>
#include <maxorb.h>
#include <ccpack.h>
#include <ccroper.h>
CKeinSPASmehr
#include <gnrinf.h>
#include <ccfield.h>
#include <r12int.h>
C
      DIMENSION WORK(LWORK)
      LOGICAL OPENED
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_INIT')
C
C----------------------------------------
C     Check whether SIRIFC is still open.
C----------------------------------------
C
      CALL GPINQ('SIRIFC','OPENE',OPENED) 
      IF (OPENED) THEN
         INQUIRE (FILE='SIRIFC',NUMBER=LUSIFC)
         WRITE(LUPRI,'(2A,I3)') ' RP_INIT: file SIRIFC is already ',
     &                          'opened with unit number ',LUSIFC
         CALL GPCLOSE (LUSIFC,'KEEP')
         WRITE(LUPRI,'(2A,I3)') ' RP_INIT: file SIRIFC was closed again'
      ENDIF
C
C---------------------------------------------------------
C     Initialize common blocks CCSDINP, CCLR, and CCSDSYM.
C---------------------------------------------------------
C
      ISYMOP = 1
CSPAS:10/11-2003: replaced with new version from CC_INP
C     SKIP   = .FALSE.
C     CCRSTR = .FALSE.
C     CCSDT  = .FALSE.
C     CC2    = .FALSE.
C     RESTLR = .FALSE.
C     MP2    = .FALSE.
C     CIS    = .FALSE.
C     CCS    = .FALSE.
C     CCD    = .FALSE.
C     CCP2   = .FALSE.
C     CC1A   = .FALSE.
C     CC1B   = .FALSE.
C     CC2A   = .FALSE.
C     CC2B   = .FALSE.
C     CCPT   = .FALSE.
C     CCP3   = .FALSE.
C     CCRT   = .FALSE.
C     MCC2   = .FALSE.
C     CCR3   = .FALSE.
C     CCR1A  = .FALSE.
C     CCR1B  = .FALSE.
C     RSPIM  = .FALSE.
C     LSEC   = .FALSE.
C     LCOR   = .FALSE.
C     NOCCIT = .FALSE.
C     FROIMP = .FALSE.
C     FROEXP = .FALSE.
C     CCFOP  = .FALSE.
C     CCSTST = .FALSE.
C     T2TCOR = .TRUE.
C     DEBUG  = .FALSE.
C     CCLR   = .FALSE.
C     CCQR   = .FALSE.
C     CCCR   = .FALSE.
C     CC4R   = .FALSE.
C     CC5R   = .FALSE.
C     CCSM   = .FALSE.
C     CCTM   = .FALSE.
C     CCLRSD = .FALSE.
C     CCQR2R = .FALSE.
C     CCEXLR = .FALSE.
C     F1SKIP = .FALSE.
C     R1SKIP = .FALSE.
C     L1SKIP = .FALSE.
C     RESKIP = .FALSE.
C     E0SKIP = .FALSE.
C     L0SKIP = .FALSE.
C     IMSKIP = .FALSE.
C     LISKIP = .FALSE.
C     M1SKIP = .FALSE.
C     FRSKIP = .FALSE.
C     LESKIP = .FALSE.
C     B0SKIP = .FALSE.
C     O2SKIP = .FALSE.
C     R2SKIP = .FALSE.
C     X2SKIP = .FALSE.
C     F2SKIP = .FALSE.
C     L2SKIP = .FALSE.
C     MAXRED = 3
C     MXLRV  = MAXRED - 1
C     MXDIIS = 8
C     MAXITE = 40
C     THRLEQ = 1.0D-05
C     THRENR = 1.0D-08
C     NFIELD = 0
C     NSIMLE = 0
C
C     CCSAVE = .TRUE.
C
C     NPRLBL = 0
C
C     Default values based on gnrinf common block
C
C     DIRECT = DIRCAL
C     IPRINT = 1
C
C     CALL IZERO(NRHFFR,8)
C     CALL IZERO(NVIRFR,8)
C     CALL IZERO(KFRRHF,8*MAXFRO)
C     CALL IZERO(KFRVIR,8*MAXFRO)
C     CALL IZERO(NOMINP,8)
C
C     CALL DZERO(EFIELD,MXFELT)
C
C     ITEST  = 0
C     MINSCR = .TRUE.
C     MINMEM = .FALSE.
C
      LGLO = .FALSE.
      SKIP   = .FALSE.
      CCRSTR = .FALSE.
      CCSDT  = .FALSE.
      CC2    = .FALSE.
      MP2    = .FALSE.
      CIS    = .FALSE.
      CCS    = .FALSE.
      CCD    = .FALSE.
      CCP2   = .FALSE.
      CC1A   = .FALSE.
      CC1B   = .FALSE.
      CCPT   = .FALSE.
      CCP3   = .FALSE.
      CCRT   = .FALSE.
      CCR3   = .FALSE.
      CCR1A  = .FALSE.
      CCR1B  = .FALSE.
      RSPIM  = .FALSE.
      TRIPIM = .FALSE.
      LSEC   = .FALSE.
      LCOR   = .FALSE.
      NOCCIT = .FALSE.
      FROIMP = .FALSE.
      FROEXP = .FALSE.
      CCFOP  = .FALSE.
      CCSTST = .FALSE.
      T2TCOR = .TRUE.
      DEBUG  = .FALSE.
      CCLR   = .FALSE.
      CCQR   = .FALSE.
      CCCR   = .FALSE.
      CC4R   = .FALSE.
      CC5R   = .FALSE.
      CCTM   = .FALSE.
      CCMCD  = .FALSE.
      CCLRSD = .FALSE.
      CCQR2R = .FALSE.
      CCOPA  = .FALSE.
      CCXOPA = .FALSE.
      CCTPA  = .FALSE.
      CCEXLR = .FALSE.
      CCDERI = .FALSE.

      LVVVV  = .FALSE.

      L0SKIP = .FALSE.

      F1SKIP = .FALSE.
      R1SKIP = .FALSE.
      L1SKIP = .FALSE.

      RCSKIP = .FALSE.
      FCSKIP = .FALSE.
      LCSKIP = .FALSE.

      O2SKIP = .FALSE.
      R2SKIP = .FALSE.
      X2SKIP = .FALSE.
      F2SKIP = .FALSE.
      L2SKIP = .FALSE.

      CO2SKIP= .FALSE.
      CR2SKIP= .FALSE.
      CX2SKIP= .FALSE.
      CF2SKIP= .FALSE.
      CL2SKIP= .FALSE.

      RESKIP = .FALSE.
      LESKIP = .FALSE.
      E0SKIP = .FALSE.
      IMSKIP = .FALSE.
      LISKIP = .FALSE.
      M1SKIP = .FALSE.
      FRSKIP = .FALSE.
      BESKIP = .FALSE.
      N2SKIP = .FALSE.
      BRSKIP = .FALSE.
      ETADSC = .FALSE.

      FREEZE = .FALSE.
      CCSLV  = .FALSE.
      NCCSLV = 0

      HERDIR = .FALSE.

      ANAAOD   = .FALSE.
      MAXRED   = 200
      MXLRV    = MAXRED - 1
      MXDIIS   = 8
      MAXITE   = 40
      NOEONL   = .FALSE.
      THRLEQ   = 1.0D-05
      THRENR   = 1.0D-08
      THRLDPHF = 1.0D-13
CWMK  NFIELD = 0
      NSIMLE = 0
      NEWCAU = .FALSE.
      LBUF   = 250000
C
CSPAS:14/11-2003: I have to play with these options later on
      KEEPAOTWO = 0
      KEEPAOIN  = .TRUE.
      LPACKINT  = .FALSE.
CKeinSPASmehr
      THRPCKINT = 1.0D-15
C
CSPAS 22.10.2003 implementing SOPPA(CCSD)
      SIRSOP    = .FALSE.
CKeinSPASmehr
C
CSPAS 15.03.2006 implementing AO-SOPPA
      AOSOPPA   = .FALSE.
CKeinSPASmehr
C
      CONNECTION = 'SYMMETR'
C
      CALL CC_NODINP(WORD,.TRUE.)
C
      MSYMS = MSYM 
C----------------------------------------------------------------
C Initialize additional labels for CCSD(T)-gradient optimization
C and FOP
C----------------------------------------------------------------

      ETACCPT = .FALSE.
      DIRKAPB = .FALSE.
C
C---------------------------------------------------------------------
C     CCSAVE is true for all case, but can be turned of in *ccexci
C     Only used for linear response - mandatory for oscstr and 
C     polarizabilities
C---------------------------------------------------------------------
C
      CCSAVE = .TRUE. 
C
C     Initialize the length of the list of response operator labels
C
      NPRLBL = 0
C
C     Default for evaluation of first derivative integrals
C
      DIRGRD = DIRCAL
C
C     Default values based on gnrinf common block
C
      DIRECT = DIRCAL
      IPRINT = IPRUSR
C
      CALL IZERO(NRHFFR,8)
      CALL IZERO(NVIRFR,8)
      CALL IZERO(KFRRHF,8*MAXFRO)
      CALL IZERO(KFRVIR,8*MAXFRO)
      CALL IZERO(NOMINP,3*8)
C
      CALL DZERO(EFIELD,MXFELT)
C
      ITEST  = 0
C
      MINSCR = .TRUE.
      MINMEM = .FALSE.
      CCPAIR = .FALSE.
CSPAS: 15.11.2009 including AO-SOPPA: MKVAJKL is no longer in COMMON
C     MKVAJKL = .FALSE.
CKeinSPASmehr
ckr      R12CAL = .TRUE.
ckr      R12NOA = .FALSE.
ckr      R12NOP = .FALSE.
ckr      R12NOB = .FALSE.
ckr      R12HYB = .TRUE.
ckr      NORXR  = .FALSE.
cwk      VCLTHR =  0D0
cwk      SVDTHR =  1D-15
ckr      R12XXL = .FALSE.

      BRASCL = 0.5D0
      KETSCL = 2.0D0
CKeinSPASmehr
C
C     Use LABEL (WK/UniKA/04-11-2002).
      LABEL = 'TRCCINT '
C
      CALL CCSD_INIT1(WORK,LWORK)
C
      IF (.NOT. DIRECT) THEN
         CALL CCSD_SORTAO(WORK,LWORK)
      ENDIF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('RP_INIT')
C
      RETURN
      END
