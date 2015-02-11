!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2015 (2015), see http://daltonprogram.org"
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
C
C  /* Deck ccsd_init0 */
      SUBROUTINE CCSD_INIT0(WORD)
C 
C-------------------------------------------------
C
C     08-November 2013 Made by Stephan P. A. Sauer
C
C     Initialize /CCSDINP/ ,/CCLR / and /CCSDSYM/
C     for the CC and AOSOPPA modules 
C
C-------------------------------------------------
C
#include "implicit.h"
#include "ccsdinp.h"      
#include "cclr.h"         
#include "ccsdsym.h"      
C
#include "ccorb.h"        
#include "ccsections.h"   
#include "maxorb.h"      
! maxorb.h needed by ccpack.h
#include "ccpack.h"       
#include "ccroper.h"      
#include "ccrspprp.h"      
#include "gnrinf.h"       
#include "ccfield.h"      
#include "r12int.h"       
#include "mxcent.h"
#include "eribuf.h"        
Cholesky
#include "chodbg.h"       
CSPAS:8/11-13: The following common blocks are probably not needed
C#include "inftap.h"
C#include "ccfop.h"
C#include "leinf.h"
C#include "cbieri.h"
C#include "cch2d.h"
CCholesky
C#include "cc_cho.h"
C#include "ccdeco.h"
C#include "chomp2.h"
C#include "chocc2.h"
CC
C#include "center.h"
CKeinSPASmehr
C
      CHARACTER WORD*7
C
C     Initialize /CCSDINP/ ,/CCLR / and /CCSDSYM/
C
      LGLO = .FALSE.
      SKIP   = .FALSE.
      CCRSTR = .FALSE.
      CCSDT  = .FALSE.
      CC2    = .FALSE.
      MP2    = .FALSE.
      DCPT2  = .FALSE.
      CIS    = .FALSE.
      CCS    = .FALSE.
      CCD    = .FALSE.
!SONIA/FRAN
      rCCD   = .FALSE.
      drCCD  = .FALSE.
      rTCCD  = .FALSE.
      SOSEX  = .FALSE.
!
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
!     Lanczos
      CCLRLCZ = .FALSE.
C     
Cholesky
C 
      CHOPT  = .FALSE.
      CALL CC_CHOPTINIT
      CHODBG = .FALSE.
C     
Cholesky
C 
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
      ONLYMO = .FALSE.
      CCSLV  = .FALSE.
CSPAS:8/11-13: NCCSLV of ccslvinf.h is initialized here, but this common block
C              was not included in CC_INPUT
C     NCCSLV = 0
CKeinSPASmehr

      HERDIR = .FALSE.

      NOSORT = .FALSE.

      ANAAOD   = .FALSE.
      MAXRED   = 200
      MXLRV    = MAXRED - 1
      MXDIIS   = 8
      MAXITE   = 40
      NOEONL   = .FALSE.
      THRLEQ   = 1.0D-05
      THRENR   = 1.0D-08
      THRVEC   = 1.0D9
      MTRIP    = .FALSE.
      THRLDPHF = 1.0D-13
CWMK  NFIELD = 0
      NSIMLE = 0
      NEWCAU = .FALSE.
      LBUF   = 250000
C
      KEEPAOTWO = 0
      KEEPAOIN  = .FALSE.
      LPACKINT  = .FALSE.
      THRPCKINT = 1.0D-15
C
C DRCCD/DRPA related stuff
C
C IT2UPD: decide on doubles amplitude updating scheme for RCCD and DRCCD
C         (not used for other methods)
C   =0 --- use conventional MP2-like update
C   =1 --- use gradient update of Henderson and Scuseria,
C          Mol. Phys. 108, 2511-2517 (2010)
      IT2UPD=1
C
C IT2START: decide on doubles amplitude start guess for RCCD and DRCCD
C         (not used for other methods)
C   =-1 --- use zero amplitudes as initial guess (DEC-CC)
C   =0  --- use conventional MP2 guess
C   =1  --- use gradient start guess along the lines of Henderson and Scuseria,
C           Mol. Phys. 108, 2511-2517 (2010)
      IT2START=1
C
C HURWITZ_CHECK: check that solution is stabilizing for DRCCD
      HURWITZ_CHECK=.FALSE.
C end

C
C     initialize AO-SOPPA variables
C
      SIRSOP    = .FALSE.
      AOSOPPA   = .FALSE.
C
      CONNECTION = 'SYMMETR'
C
      CALL CC_NODINP(WORD,.TRUE.)
C
C----------------------------------------------------------------
C Initialize additional labels for CCSD(T)-gradient optimization
C and FOP
C----------------------------------------------------------------

      ETACCPT = .FALSE.
      DIRKAPB = .FALSE.
C
C-----------------------------
C     Cholesky initializations
C-----------------------------
C
      CALL CC_CHOMP2INIT
      CALL CC_CHOCC2INIT
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
      NPRLBL_CC = 0
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
CSPAS:8/11-13: NOMINP of ccexcinf.h is initialized here but the common block was
C              not included in CC_INPUT
C     CALL IZERO(NOMINP,3*8)
CKeinSPASmehr
C
      ITEST = 0
      DO I=1,8
        NRHFFR(I) = LOCFRO(I)
        ITEST = ITEST + NRHFFR(I)
      END DO
      FROIMP = ITEST .NE. 0
C
      CALL DZERO(EFIELD,MXFELT)
C
      ITEST  = 0
C
      MINSCR = .TRUE.
      MINMEM = .FALSE.
      CCPAIR = .FALSE.
ccn      MKVAJKL = .FALSE.
ckr      R12CAL = .TRUE.
ckr      R12NOA = .FALSE.
ckr      R12NOP = .FALSE.
ckr      R12NOB = .FALSE.
ckr      R12HYB = .TRUE.
ckr      NORXR  = .FALSE.
cwk      VCLTHR =  0D0
cwk      SVDTHR =  1D-15
ckr      R12XXL = .FALSE.

      BRASCL = 1.0D0
      KETSCL = 1.0D0
C
      RETURN
      END
