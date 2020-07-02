!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
C
!===========================================================================
!Revision 1.9  2000/05/24 12:39:05  hjj
!only MAXCIT=3 as default if MC combination with DOSCF or DOMP2
!/hjaaj dec 1999 - apr 2000:
!1) removed obsolete RNALL references
!2) do not activate QCSCF (flag(21)) for solvent or writing interface
!   if MXDIIS.gt.0; solvent and CALL WR_SIRIFC are now also implemented
!   in the DIIS module.
!----------------------------
!revision 1.3 2000/04/14 08:23:47 hjj
!removed "OLDINP", old obsolete input for SIRIUS /hjaaj apr 2000
!----------------------------
!revision 1.2 2000/04/14 07:07:51 hjj
!1) initialize LUSUPM = -9999, to signalize default is use AOSUPMAT
!   (before GPOPEN this was done with LUSUPM=17 in SIRINI)
!2) changed many LUERR to LUPRI in SIRRDI
!=============================================================================
!951017 C*FM
!PRTINP : FM change format nb of determinantsi and CSF  info
!901026-hjaaj:
!-- indfoert common block SCBRHF options i SIRRDI, SCFINP
!-- OPTINP: FLAG(66) = .TRUE. when .DETERMINANTS
!-- indsat MOTECC-90 in-line documentation
!-- SIRRDI: check CORHOL calculations
!-- PRTINP: tell if core hole calculation or open shell HF
!900920-hjaaj:
!-- OPTINP: .COREHOLE keyword + .CORERELAX keyword
!-- SIRRDI: CORHOL and CORRLX implemented, no check if calc. ok yet!
!===========================================================================
C  /* Deck sirinp */
      SUBROUTINE SIRINP(WORK,LFREE)
C
C Revisions:
C   10-Jul-1984 hjaaj
C    7-Jan-1985 hjaaj
C      Oct-1990 hjaaj
C
C ***** PROGRAM SIRIUS, CONTROL ROUTINE FOR INPUT SECTION *****
C
C OBJECTIVE: READ AND CHECK INPUT DATA
C
C MOTECC-90: This module, SIRINP, is described in the input/ouput
C            Documentation of MOTECC-90.
C
C SUBROUTINE CALLS: SIRRDI (READ AND CHECK INPUT DATA)
C                   ORTHO  (REORDERING AND ORTHONORMALIZATION
C                           OF TRIAL MOLECULAR ORBITALS).
C
C
#include "implicit.h"
      DIMENSION WORK(LFREE)
C
C -- local constants
C
#include "dummy.h"
C
C -- common blocks
C    INFINP : FLAG(*),?
C    INFORB : NCMOT,...
C    INFPRI : LUPNCH
C
#include "maxorb.h"
#include "priunit.h"
#include "infinp.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
C
      CALL QENTER('SIRINP')
C
C *****            READ AND CHECK INPUT DATA           *****
C ***** NOTE: IT IS ASSUMED THAT LFREE IS LARGE ENOUGH *****
C ***** ====  TO STORE ONE SET OF MOLECULAR ORBITALS.  *****
C
C     Molecular orbitals are returned in WORK(1)
C
      CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      CALL SIRRDI(WORK(1),LFREE,IREST)
      CALL GPCLOSE(LUCMD,'KEEP')
C
C ***** ORTHONORMALIZE ORBITALS, UNLESS RESTART
C       feb. 2003 hjaaj: now always orthonormalize orbitals,
C       becuase the user may have used .RESTART with
C       SIRIUS.RST from another geometry ...
C
Chj   IF (IREST.EQ.0 .OR. IORTO.EQ.1) THEN
         KW1=1
         KW2=KW1+NCMOT
         KW3=KW2+NNBAST
         IF (KW3.GT.LFREE) CALL ERRWRK('SIRINP.ORTHO',-KW3,LFREE)
         LW3 = LFREE - KW3
         CALL ORTHO(WORK(1),WORK(KW2),WORK(KW3),LW3)
C        CALL ORTHO(CMO,S,SIN,LSIN)
Chj   END IF
C
C     If (FLAG(31)) then punch orbitals to LUPNCH
C
      IF (FLAG(31)) THEN
         CALL PUNMO(0,WORK,DUMMY)
C        CALL PUNMO(IPCTL,CMO,OCC)
      END IF
C
C
C ***** Write molecular orbitals *****
C       Stop if only reorthonormalization.
C
      IF (IORTO.EQ.1) THEN
         CALL NEWORB('REORTHO ',WORK,.TRUE.)
         CALL QUIT(
     &      'SIRIUS NORMAL STOP AFTER ORBITAL ORTHONORMALIZATION.')
Chj   ELSE IF (IREST .EQ. 0) THEN
      ELSE
         CALL NEWORB('REORTHO ',WORK,.FALSE.)
      END IF
C
C end of SIRINP
C
      CALL QEXIT('SIRINP')
      RETURN
      END
C  /* Deck sirrdi */
      SUBROUTINE SIRRDI(CMO,LCMO,IREST)
C
C Revisions
C   14-Oct-1984 hjaaj (included ACTROT option)
C    6-Aug-1984 hjaaj / 15-May-1984 hjaaj
C   23-Mar-1985 hjaaj (new namelist options, common /INFTRA/ for TRACTL)
C   10-Nov-1985 hjaaj
C    1-OCT-1987 HA     Keyword directed input for all parameters
C
C OBJECTIVE: Read input data from file LUCMD and check,
C            get initial molecular orbitals using READMO or from LUCMD.
C
C Input:
C  none
C
C Output:
C  CMO, starting molecular orbitals
C  IREST = 0 not a restart
C        = 1 if old, standard restart for SIRIUS.RST
C        > 1 if start from old SIRIUS.RST or SIRIFC,
C             e.g. from another geometry
C
C Called from SIRINP
C
C Subroutine calls: MOLLAB (FIND RECORD ON DATA FILES)
C                   SIR_PRTINP (PRINT INPUT)
C                   ERRINP (ERROR message OUTPUT ROUTINE)
C
C******************************************
C
      use pelib_interface, only: use_pelib
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "ecpinf.h"
      DIMENSION CMO(LCMO)
C
C -- local variables and constants
C
      CHARACTER WORD*8
      LOGICAL   DELEMO, FNDLAB
      PARAMETER (D0=0.0D0, D1=1.0D0, D2=2.0D0)
      PARAMETER (THR_DEFAULT=1.0D-5, THRRHF_MP2=1.0D-6)
C
C Used from common blocks:
C   GNRINF : WESTA,EMBEDDING,QM3,QMMM, ...
C   CCOM   : THRS ! integral threshold, used for THRP default
C   INFINP : ...,ITRLVL,ITRFIN,..., DO*, DIRFCK,...
C   INFORB : NSYM,?
C   INFVAR : MAXOCC,NWOPT,JWOP(2,*),JWOPSY
C   INFOPT : EPOT,?
C   INFTAP : LUCME,?
C   INFDIM : MAXPHP,?
C   INFTRA : THRP,USEDRC,IPRTRA
C   SCBRHF : IOPRHF, NFRRHF(8), MXHFMA, MXHFMI, RHFCAN, MAXCK, MXDIIS, THRRHF
C   INFMP2 : NFRMP2(*), SRMP2_SELFCONSISTENT
C   INFPOP : IPOPTYP; depends on mxcent.h
C   INFLOC : BOYSEL
C   PCMLOG : PCM
C   HUCKEL : DOHUCKEL
C   CBIHRS : RUNSUP
C   SPINFO : MULTS, MS2
C
#include "gnrinf.h"
#include "maxash.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"

#include "ccom.h"
#include "infinp.h"
#include "inforb.h"
#include "infind.h"
#include "infvar.h"
#include "infopt.h"
#include "inftap.h"
#include "infpri.h"
#include "infdim.h"
#include "inftra.h"
#include "scbrhf.h"
#include "infmp2.h"
#include "infpop.h"
#include "infloc.h"
#include "pcmlog.h"
#include "huckel.h"
#include "cbihrs.h"
#include "spinfo.h"
#include "dftcom.h"
#include "dfterg.h"
#include "stex.h"
#include "center.h"

C
      INTEGER*4   NSYM_i4, NBAS_i4(8),NORB_i4(8), NRHF_i4(8), IOPRHF_i4
      ! hjaaj Apr 2011: such that Dalton with 8 byte integers can read a
      ! SIRIUS.RST from a Dalton with 4 byte integers, and vice versa.
      INTEGER MBAS(8)
      LOGICAL ONELOK

      CALL QENTER('SIRRDI')
C
C     Skip next section if processed earlier
C
      IF (SIR_INPPRC) THEN
         INPERR = 0
         NUMRUN = 1
C        ... only highest calculation specification is done
C
C     Several common block are restored
C     (MUST be called before NEWSYM reset to false)
C
         CALL SRESTR
         IF (RNRESP .OR. RNABAC) THEN
C           ... write SIRIFC and calculate WOP of all symmetries
C               for RESPONS and/or ABACUS
            FLAG(25) = .TRUE.
            FLAG(11) = .TRUE.
         END IF
         IF (WESTA) THEN
C           ... write SIRIFC and CI info for WESTA program
            FLAG(25) = .TRUE.
            FLAG(26) = .TRUE.
         END IF
C
C     We obviously have MO's from previous iteration
C     (.MOSTART NEWORB because abawlk.F(WLKSAV) saves
C     on LUIT1 with label NEWORB)
C
         IRDMO = 9
         FLAG(13) = .TRUE.
C
C     ... unless a new symmetry has been found, vb, april-96
C
         IF (NEWSYM) THEN
            IRDMO = 8
            FLAG(13) = .FALSE.
            IF (DOSCF) THEN
               AUTOCC = NSYM .GT. 1
               INIOCC = -1
               NRHF(1:8) = 0
               IOPRHF = 0
            ELSE
               WRITE (LUPRI,'(/A)')
     &            ' Lowering of symmetry for MCSCF is not'//
     &            ' implemented in SIRIUS yet, sorry!'
               CALL QUIT('Lowering of symmetry for MCSCF not'//
     &            ' implemented in SIRIUS yet, sorry!')
            END IF
            ICI0 = 1
         ELSE
            ICI0 = 6
            AUTOCC = .FALSE.
            INIOCC = 1
C           ... SCF occupation is already in NRHF() and IOPRHF
         END IF
C
C     We want to use Gram-Schmidt orthogonlization in order to change
C     our MO's as little as possible when we read them from file
C
         FLAG(32) = .FALSE.
C
C     We also want to check the geometry prediction (ICI0=6)
C     If excited state force ALWAYS Newton-Raphson
C     because CI vectors for lower roots needed for NEO not available
C
         IF (ISTATE.GT.1 .AND. NMCAVER.EQ.1) FLAG(39) = .TRUE.
         IREST = 0
      ELSE
C
C  ***** READ AND CHECK LABEL *****
C
C
         REWIND (LUCMD,IOSTAT=IOS)
C        ... IOSTAT to avoid program abort on some systems
C         if reading input from a terminal
  900    READ (LUCMD,'(A)',END=910,ERR=920) WORD
         CALL UPCASE(WORD)
         IF (WORD(1:8) .NE. '**SIRIUS' .AND.
     &       WORD(1:8) .NE. '**WAVE F' .AND.
     &       WORD(1:8) .NE. '**WAVEFU') THEN
            GO TO 900
         END IF
         GO TO 930
  910    CONTINUE
            CALL QUIT('ERROR: '//
     &      'End of .dal input file, no wave function input found')
  920    CONTINUE
            CALL QUIT('ERROR reading .dal '//
     &      'input file before finding wave function input')
  930    CONTINUE
      END IF
C
C ***** Read NSYM, NBAS(1:NSYM), POTNUC and CGTO data from LUONEL
C
      ONELOK = .FALSE.
      CALL RDONEL('OPEN    ',ONELOK,CMO,LCMO)
      IF (ONELOK) THEN
         CALL RDONEL('MLCLINFO',.TRUE.,CMO,LCMO)
      ELSE
         NWARN = NWARN + 1
         WRITE (LUPRI,'(/2A/A)')
     *      ' SIRRDI WARNING: could not open LUONEL file ',FNONEL,
     *      ' NBAS(*) set to zero.'
         CALL IZERO(NBAS,8)
         NSYM   = 8
         POTNUC = D0
      END IF
      EPOT = POTNUC
C
C     Initialize local variables and NORB(*)
C
      DO 1670 I = 1,8
         NORB(I) = NBAS(I)
 1670 CONTINUE
C
C     *********************************************
C     ***** SET DEFAULT VALUES FOR INPUT DATA *****
C     *********************************************
C     (most variables are in infinp.h)
C
C     Jump again if processed earlier
C
      IF (SIR_INPPRC) GOTO 124
C
      IREST = 0
C     IREST .ne. 0 for a restart run
C
C********** Default values section SIR_GENINP *****
C
C     NTIT     Number of title cards given
C     IORTO    = 1 stops calculation after orthogonalization
      NTIT  = 0
      IORTO = 0
C
      DO 660 I = 1,NFLAG
         FLAG(I) = .FALSE.
  660 CONTINUE
      DO 670 I = 1,NPFLAG
         P4FLAG(I) = .FALSE.
         P6FLAG(I) = .FALSE.
  670 CONTINUE
C
      DOSCF     = .FALSE.
      DODFT     = .FALSE.
      DOMP2     = .FALSE.
      DOCI      = .FALSE.
      DOMC      = .FALSE.
      DONEVPT   = .FALSE.
      DOLUCITA  = .FALSE.
      DOAPSG    = .FALSE.
      DO_VIRTRUNC = .FALSE.
C  ---> Initialize SR-DFT stuff here !
      DOMCSRDFT = .FALSE.
      DOCISRDFT = .FALSE.
      DOHFSRDFT = .FALSE.
      ADDSRI    = .FALSE.
      SRHYBR    = .FALSE.
      DFT_SPINDNS = .FALSE.
      DFT_LOCALSPIN = .FALSE.
C <---
      DOCINO    = .FALSE.
      DOFCI     = .FALSE. ! Do Full CI
      DOCISD    = .FALSE. ! DO CISD
      DORSP     = RNRESP
      ABAIPH    = RNABAC
C     ... It is followed by an abacus run
      CORHOL    = .FALSE.
      CORRLX    = .FALSE.
C added by Bin Gao for generating cube files
C March 9, 2012
      DO_CUBE = .false.
C
Cholesky orbitals
      ACTSEL = .FALSE.
C
C********** Default values section HAMINP *****
C
CWMK  NFIELD = 0
C
C********** Default values section WAVINP and more *****
C
C
      DO 650 I = 1, 8
         NFRO(I) = 0
         NRHF(I) = 0
         NFRRHF(I) = 0
         NMVO(I) = 0
         NFRMP2(I) = 0
         NISH(I) = 0
         NASH(I) = 0
         NAS1(I) = 0
         NAS2(I) = 0
         NAS3(I) = 0
  650 CONTINUE
      NMVOT = 0
      NASHT = 0
      CALL IZERO(NOROT,MXCORB)
      LNOROT = .FALSE.
      NACTEL = -1
      ISPIN  = -2
      LSYM   = -1
      MCTYPE = 0
      NELMN1 = 0
      NELMX1 = 0
      NELMN3 = 0
      NELMX3 = 0
      MULTS  = -2     ! MULTS .lt. 0 is code for not set yet
      MS2    = -12345 ! -12345 is code for not set by .MS2 option
C
C********** Default values section STEX *****
C
      DOSTEX=.FALSE.
C
C********** Default values section MP2INP *****
C
      MP2_LSHIFT      = D0
      MP2_SCALEFAC(1) =  D2 ! scale factor on (ij|ab)
      MP2_SCALEFAC(2) = -D1 ! scale factor on <ij|ab> = (ia|jb)
      LAMSR           = D0
      MP2_SCALED   = .FALSE.
      MP2_SCS      = .FALSE.
      MP2_SOS      = .FALSE.
      MP2_TDA      = .FALSE.
      MP2_NO_OCCVIR= .FALSE.
      SRMP2_SELFCONSISTENT = .FALSE.
      SRMP2_SRINTS = .FALSE.
      DCPT2        = .FALSE.
      DO_RSDHF     = .FALSE.
      SAVE_MP2WF1  = .FALSE.
C
C********** Default values section POPINP *****
C
      IPOPTYP = -1
C
C
C********** Default values section OPTINP *****
C
C     MAXMAC   Maximum number of macro iterations
C     MAXMIC   Maximum number of micro iterations (reset later)
C     MAXJT    Max number of micro iterations in any macro iteration
      MAXMAC = -1
      MAXMIC = -1
      MAXJT  = -1
      IMCCNO = 0
C
C     THR* : Convergence thresholds.
C            Default set later, now set to -1.0d0 so we
C            can check later if they have been set by user.
C            /hjaaj June09
C     THRGRD   Threshold for convergence on gradient in MC/QC-SCF opt.
C
      THRRHF = -D1
      THRCI  = -D1
      THRMC  = -D1
      THRGRD = -D1
C
C     NROOTS   Number of roots used in micro iterations
C              LROOTS number of start trial vectors,
C              NROOTS number of simultaneous trial vectors
C              in each micro iteration, after the very first.
C     IROOT(*) Sequence numbers for roots to be used in optimization
C     ISTATE   always = IROOT(1), the root you want to converge to.
      NROOTS = 1
      LROOTS = 1
      IROOT(1) = 1
      ISTATE   = 1
      NMCAVER  = 1
      WEIGHT_MCAVER(1) = 1.0D0
C
C     Initialize to zero orbital rotations
      NWOPT = 0
C THRPWF: threshold for print CI-coefficient in PRWF
      THRPWF = -D1
      THRCGR = 0.1D0
C ICHECK: symmetry check of CI vectrors
C         (default 2, choose symmetry of lowest CI root)
      ICHECK = -100
      MAXCIT = -1
      MAXABS = 3
      MAXAPM = 5
C MAXRL : Max dimension of reduced L matrix
      MAXRL  = 120
C Use average (approximate) orbital Hessian diagonal
      FLAG(12) = .FALSE.
      RESPHP = .FALSE.
      JOLSEN = .FALSE.
      MAXPHP = 0
!
!     default CI program used as MCSCF backbone 
      ci_program = 'SIRIUS-CI'

      !> set default level (no spin-densities requested) for spin-density calculation (works only for LUCITA)
      spindens_lvl = 0
C
C********** Default values section CICINP *****
C
C MAXRC : Max dimension of reduced CI matrix
      MAXRC  = 150
      ISTACI = 0
C     ... do not converge to a specific state.
      NROOCI = 1
      MXCIMA = 25
      ICICNO = 0
C
C********** Default values section SCFINP *****
C           and common scbrhf.h
C           (except DIRFCK and HSROHF which are in infinp.h)
C
      IOPRHF = 0
      MXDIIS = 60
      MAXEVC = -1
      RHFCAN = .TRUE.
      MAXFCK = 0
      MXHFMA = 60
      MXHFMI = 12
      HSROHF = .FALSE.
      BCKSTP = .TRUE.
      CALL IZERO(NROHF,8)
C 960625-hjaaj: NISHT and NACTEL used in MOLCHR in
C               ICHARG = ICHRGNU - (2*NISHT+NACTEL)
      NISHT  = 0
      NACTEL = 0
      CALL MOLCHR(ICHARG)
C     ICHARG is now total nuclear charge;
C     KCHARG is input in .mol file for total molecular charge
      NACTEL = -1
      NRHFEL = ICHARG - KCHARG
      INIOCC = -1
      IF (NRHFEL .GT. 0) THEN
         AUTOCC = NSYM .GT. 1
      ELSE
         AUTOCC = .FALSE.
      END IF
      DIRFCK = DIRCAL
C
C********** Default value sections SCFINP and AUXINP *****
C           (default: use AOSUPMAT if available,
C                     LUSUPM=-1 means do not use AOSUPMAT)
C
      IF (.NOT. RUNSUP) THEN
         LUSUPM = -1
      ELSE
         LUSUPM = -9999
      END IF
C
C********** Default values section DFTINP *****
C
      CALL DFTINI
C
C********** Default values section dfterg.h *****
C
C (initialize energy values to zero, so they are defined
C  for non-DFT calcualtions, i.e. standard HF or MCSCF)
C
      EDFTY   = D0
      ESRDFTY = D0
      ESRDFT  = D0
C
C
C********** Default values section ORBINP *****
C
C     LMOORD   True if mo:s have to be reordered
C     BOYORB   True if Boys localization to be used
C     PIPORB   True if Pipek localization to be used
C
      LMOORD = .FALSE.
      BOYORB = .FALSE.
      PIPORB = .FALSE.
C
C     THRSSY   the threshold for identifying supersymmetries and
C              their degeneracy based on matrix elements of the
C              kinetic energy matrix.  Larger values (as 5.D-8) means
C              more approximate symmetries will be enforced (e.g.
C              C3v with 8 digits in the geometries, which will not
C              be recognized with the old value of 1.0D-9). On the
C              other hand, larger values increase the risk for
C              accidental degeneracies, which will cause the super
C              symmetry analysis to fail.
C
      THRSSY = 5.0D-8
      SUPSYM = .FALSE.
C     ... Oct 2003 hjaaj: new default
c     SUPSYM = .NOT. (RNABAC .OR. RNRESP)
C     supsym false forced for DOCCSD set after reading of sirius input
C
C IRDMO : select how to get starting orbitals
C         (default 8, diagonalize one-el. Ham.
C          corresponding to .MOSTART( H1DIAG )
C KDEL  : if .gt. 0, then delete 3s component in d; 4p in f; etc.
C         If IRDMO .eq. 8, i.e. new calculation, call CMODEL
C         to find 3s,4p,... components and delete them.
C         if .lt. 0, orbitals have been deleted in input
C CMAXMO: max allowed MO coefficient (to ensure 2-el integral accuracy)
C THROVL: limit for basis set numerical linear dependence
C DEFLVL: default shift level - must be 0 or negative.
C
      IRDMO  = 8
      FLAG(13) = .FALSE.
      KDEL   = 0
      THROVL = 1.0D-6
      CMAXMO = 1.0D3
      FCVORB   = .FALSE.
      DEFLVL = 1.0D0 ! code value to tell not set by user
C
C********** Default values section CIVINP *****
C
C ICI0  : select how to get starting CI vector in SIROPT.OPTST
C         (default 1, start with MAXCIT CI iterations)
C
      ICI0   = 1
C
C********** Default values section STPINP *****
C
      BETA   = 1.0D0
      BETMIN = 0.2D0
      BETMAX = 1.D6
      RTTOL  = 1.1D0
      RATREJ = 0.25D0
      RATMIN = 0.40D0
      RATGOD = 0.80D0
      STPMAX = 0.70D0
      STPINC = 1.2D0
      STPRED = 0.67D0
      RTRUST = STPMAX
      THQMIN = 0.1D0
      THQLIN = 0.2D0
      THQKVA = -1.0D0
C
C********** Default values section POPINP *****
C
      IPRMUL = 1
C
C********** Default values section TRAINP *****
C
C     ITRLVL   Type of integral transformation level.
C              Zero for first order; two,three, or four for SIRIUS.
C              Five for MP2
C
C     ITRFIN   Level of final integral transformation (when FLAG(25),
C              i.e. call wrgeom in siropt)
C
      THRP   = MAX(THRS,1.0D-15)
      ITRLVL = -6789
      ITRFIN = -6789
      IPRTRA = -1
C     MWORK is defined in block data SBDDIM
C
      USEDRC= .FALSE.
C     ... USEDRC may be reset by transformation module
      AO2INTFILE_LABEL = 'AOTWOINT'
C
C
C********** Default values section UPDINP *****
C
C Default : do not switch to update
      FLAG(37) = .TRUE.
      MAXUIT = 48
C
C********** Default values section SOLINP *****
C
      EPSOL  = -D1
      EPSTAT = -D1
      RSOL(1)= -D1
      RSOL(2)= -D1
      RSOL(3)= -D1
      LSOLMX = -1
      INERSI = .FALSE.
      INERSF = .FALSE.
C
C********** Default values section PRIINP *****
C +++ remember to change reference manual if values are changed!
C
C     IPRI*    Print parameters for output on LUW4 and LUPRI
C     LUW*     Output units for print output (now defined in SIRSET)
      IPRI4  = IPRUSR + 5
      IPRI6  = IPRUSR
      IPRSIR = IPRUSR
      LUW4   = LUPRI
      IPRAVE = IPRUSR
      IPRKAP = IPRUSR
      IPRCIX = IPRUSR
      IPRDIA = IPRUSR - 1
      IPRSIG = IPRUSR
      IPRDNS = IPRUSR
      IPRRHF = IPRUSR - 1
      IPRMP2 = IPRUSR + 4
      IPRSOL = IPRUSR + 4
      IPRCNO = IPRUSR + 1
      IPRFCK = IPRUSR
      CMOPRI = .FALSE.
C     set limit for POPPRI in each DIIS/QCHF/MCSCF iter to 30 nuclei
      LIM_POPPRI = 30
C
C  ***************************************************
C  *****        Read Sirius input                *****
C  ***************************************************
C
C     Extra default: CALL SIROUT for summary
      FLAG(5) = .TRUE.
C
C INPERR = number of input errors. Stop if different from zero.
C NUMRUN = number of wave functions (e.g. RHF + MP2 + MCSCF: NUMRUN=3)
C
      INPERR = 0
      NUMRUN = 0
      CALL SIR_NEWINP(INPERR,NUMRUN,IRDMO,IREST,NSYM,NORB)

      IF (.NOT.DOSCF) MAXFCK = 0
      IF (.NOT.DOCI ) THEN
C     .. ignore .CINO if no CI calculation
         DOCINO = .FALSE.
         DOFCI  = .FALSE.
         DOCISD = .FALSE.
      END IF
      IF (.NOT.DOMC .AND. .NOT.DOCI) THEN
         MCTYPE = 0
         IF (HSROHF) MCTYPE = -1
      END IF
C     ... ignore any *CONFIGURATION input
C         if neither DOMC nor DOCI /hjaaj aug 2005
C
      IF (DOSCF .AND. .NOT.FLAG(21)) THEN
        IF (MXDIIS .LE. 0) THEN
C        ... March 2000/hjaaj: solvent and writing interface are
C            now also defined in DIIS module.
          IF (FLAG(16)) THEN
            NINFO = NINFO + 1
            WRITE (LUPRI,'(/5X,A/5X,A)')
     &      'SIRRDI INFO, ".NOQCSCF" ignored because solvent is only',
     &      '             implemented in DIIS and QCSCF modules.'
            FLAG(21) = .TRUE.
          ELSE IF (FLAG(25) .AND. .NOT.(DOMC .OR. DOCI)) THEN
            NINFO = NINFO + 1
            WRITE (LUPRI,'(/5X,A/5X,A)')
     &     'SIRRDI INFO, ".NOQCSCF" ignored because writing interface',
     &     '             is only implemented in DIIS and QCSCF modules.'
            FLAG(21) = .TRUE.
          END IF
        END IF
        IF (.NOT.FLAG(21)) THEN
C       ... if not QCSCF zero iteration limit (for consistent output in SIR_PRTINP)
          MXHFMA = 0
          MXHFMI = 0
        END IF
      END IF
C
c.hjaaj aug04: user's responsibility now SUPSYM false is default
c.             it should be OK if closed shell and not huckel
c     IF ((DOCCSD .OR. DOMP2) .AND. SUPSYM) THEN
c        SUPSYM = .FALSE.
c        NINFO = NINFO + 1
c        WRITE (LUPRI,'(/A)') ' INFO: SUPSYM disabled'//
c    &      ' because SUPSYM not implemented in MP2 and CC'
c     END IF
      IF (DOCCSD) THEN
         FLAG(11) = .TRUE.
         IF (DOSCF) FLAG(25) = .TRUE.
      END IF
      IF (RNRESP .OR. RNABAC) THEN
         FLAG(25) = .TRUE.
         FLAG(11) = .TRUE.
      END IF
      IF (WESTA) THEN
C        ... write SIRIFC and CI info for WESTA program
         FLAG(25) = .TRUE.
         FLAG(26) = .TRUE.
      END IF
      IF (BOYORB .OR. PIPORB .OR. BOYSEL) THEN
         IF (SUPSYM) THEN
            WRITE(LUPRI,'(/A)')
     &      'WARNING: .SUPSYM ignored because incompatible with .LOCALI'
            SUPSYM = .FALSE.
            NWARN = NWARN + 1
         END IF
      END IF
      IF (RNABAC .AND. SUPSYM .AND. ISPIN .NE. 1) THEN
      ! .SUPSYM does not work ABACUS if not closed shell
         SUPSYM = .FALSE.
         NWARN = NWARN + 1
         WRITE (LUPRI,'(/A)') ' WARNING: .SUPSYM disabled'//
     &      ' because SUPSYM not implemented in **PROPERTIES (ABACUS)'
C...     the orbital gradient may be non-zero without SUPSYM
C        and the WF is then not variational in abacus
      END IF
      IF (FLAG(21)) FLAG(2) = .TRUE.
C     ... qc-scf needs flag(2) for file opening
      IF (SUPSYM .AND. NWOPT.GT.0) THEN
         SUPSYM = .FALSE.
         NWARN = NWARN + 1
         WRITE (LUPRI,'(/A)') ' WARNING: .SUPSYM disabled'//
     &      ' because it is not compatible with .ACTROT'
C        ... active rotations will be deleted when
C            SETWOP is called again because orbitals
C            have been reordered after ORDRSS or ORD2SS /hjaaj aug 04
      END IF
      IF (HFXMU .NE. 0.0D0 .AND. .NOT.DIRCAL) THEN
         NINFO = NINFO + 1
         IF (.NOT.DIRFCK) WRITE (LUPRI,'(/A)')
     &   'INFO: switched to AO-direct Fock matrices because long-range'
     &   //' corrected models as CAM-B3LYP will not work otherwise.'
         DIRCAL = .TRUE.
         DIRFCK = .TRUE.
      END IF
C
C *** Open print unit(s).
C     and initialize units for memory allocation routines
C
      IF (LUW4 .NE. LUPRI) THEN
         LUW4 = -1
         CALL GPOPEN(LUW4,'SIRIUS.RES','UNKNOWN',' ','FORMATTED',IDUMMY,
     &               .FALSE.)
      END IF
C
C     Several common blocks are saved to next iteration
C
      CALL SSAVE
C
C *************************************************************************
C     "124" is where we jump back in if input has been processed earlier,
C     i.e. if SIR_INPPRC is .true.
C
 124  CONTINUE
C
      IF (DIRFCK) THEN
         LUSUPM = -1
      END IF
C
C     set ITRFIN for ABAIPH
!     Apr. 2016: not any more, Abacus calls TRACTL
!     and that is better, then it is known exactly what
!     transformation level is needed.
!
!     IF (ABAIPH .AND. ITRFIN .EQ. -6789) THEN
!        IF (DOMC .OR. DOCI) THEN
!           ITRFIN = 4
!           IF (NEWTRA) ITRFIN = 5 ! to get occ-sec integrals
!        END IF
!     END IF
C
      IF (DOCI .AND. ISTACI .LE. 0) THEN
         IF (DORSP .AND. .NOT.DOMC) THEN
C           MAERKE 921124-hjaaj: STOP/WARN if ISTACI not
C              defined for CIRSP ??? or ISTACI = NROOCI?
            ISTACI = MAX(1,ISTATE)
            ISTATE = ISTACI
         END IF
         IF (DOCINO) ISTACI = MAX(1,ISTATE)
      IF (ISTACI .GT. 0) THEN
         NWARN  = NWARN + 1
         WRITE (LUPRI,'(2(/5X,A),/5X,A,I3)')
     & 'WARNING  Reference state is required for CI response and CI-NO',
     & 'WARNING  ".STATE" was not specified for "*CI INPUT",',
     & 'WARNING  the reference state is set to be state no.',ISTACI
      END IF
      END IF
C
C     Save explicitly defined IPRSIR, IPRI4, IPRI6
C
      IF (IPRSIR .EQ. IPRUSR) THEN
C        .. if user hasn't set IPRSIR with .PRINT but has set
C           either IPRI4 or IPRI6, then use max value /hjaaj Oct09
         IPRSIR = MAX(IPRSIR,IPRI4-5,IPRI6)
      END IF
      MPRSIR = IPRSIR
      MPRI6  = IPRI6
      MPRI4  = IPRI4
C
      IF (THRPWF .LT. D0) THEN
         IF (ISTATE .EQ. 1) THEN
            THRPWF = 0.05D0
         ELSE
            THRPWF = 0.1D0
         END IF
      END IF
      IF (MCTYPE .LE. 0 .AND. (DOMC .OR. DOCI)) THEN
         WRITE (LUPRI,'(/A)')
     &      'INPUT ERROR: MCSCF or CI requested,'//
     &      ' but neither CAS nor RAS space defined.'
         INPERR = INPERR + 1
      END IF
      IF (DONEVPT) THEN
         IF (.NOT.DOMC .AND. .NOT.DOCI) THEN
            WRITE (LUPRI,'(/A)') 'INPUT ERROR: .NEVPT2 requires,'//
     &      ' also .CI or .MCSCF.'
            INPERR = INPERR + 1
         ELSE IF (MCTYPE .NE. 1) THEN
            WRITE (LUPRI,'(/A)') 'INPUT ERROR: NEVPT2 requested,'//
     &   ' NEVPT2 is only implemented for complecte active space (CAS).'
            INPERR = INPERR + 1
         END IF
      END IF
C
C     No need for AUTOCC if C_1 symmetry
C
      IF (NRHFEL .LE. 0) THEN
         NRHFEL = 2*NRHFT
         IF (IOPRHF .GT. 0) NRHFEL = NRHFEL + 1
      END IF
      IF (NSYM .EQ. 1 .AND. INIOCC .LT. 0) THEN
         NRHF(1:8) = 0
         MOCC = NRHFEL/2
         NRHF(1) = MOCC
         IF (2*MOCC .NE. NRHFEL) THEN
            IOPRHF = 1
         ELSE
            IOPRHF = 0
         END IF
         INIOCC = 6
         AUTOCC = .FALSE.
      END IF
C
C     Disable AUTOCC for RHF if RHF not used (HJ aug95)
C
      IF (MCTYPE .GT. 0 .AND. .NOT.DOSCF) THEN
         AUTOCC = .FALSE.
         INIOCC = -1
      END IF
C
C     *** DEFINE WAVE FUNCTION PARAMETERS FOR HARTREE-FOCK ***
C     If Hartree-Fock or MP2 is the top level, then the Hartree-
C     Fock parameters must be transferred to the CI/MCSCF wave
C     function parameters.
C
C     910723-hjaaj: define wave function for Hartree-Fock,
C     even if DOSCF false and DOMP2 false; JWOP information
C     may be needed for RESPONSE
C
      IF (MCTYPE .LE. 0) THEN
         CALL ICOPY(8,NFRRHF,1,NFRO,1)
         CALL ICOPY(8,NRHF,1,NISH,1)
         CALL IZERO(NASH,8)
         IF (IOPRHF .GT. 0) THEN
            LSYM   = IOPRHF
            ISPIN  = 2
            MS2    = 1
            NASH(LSYM) = 1
            NASHT  = 1
            NACTEL = 1
            NROHF(1:8) = NASH(1:8) ! so it will work if HSROHF switched on because open-shell DFT /hjaaj Nov 2016
         ELSE IF (HSROHF) THEN
            LSYM = 1
            DO ISYM = 2, NSYM
               IF (MOD(NROHF(ISYM),2).EQ.1) LSYM = MULD2H(LSYM,ISYM)
            END DO
            NASHT  = ISUM(NSYM,NROHF,1)
            ISPIN  = NASHT + 1
            MS2    = NASHT
            NACTEL = NASHT
            CALL ICOPY(8,NROHF,1,NASH,1)
         ELSE
            LSYM   = 1
            ISPIN  = 1
            MS2    = 0
            NASHT  = 0
            NACTEL = 0
         END IF
C
C        If NASHT .le. 1 then MO integrals not needed in optimization
C        However, if ITRLVL has been set in input, then this overrides.
         IF (ITRLVL .EQ. -6789) THEN
            FLAG(14) = .TRUE.
            FLAG(34) = .TRUE.
         ELSE
            FLAG(34) = .FALSE.
         END IF
C
      END IF
C
      IF (NFIELD .GT. MXFELT) THEN
         WRITE (LUPRI,*) 'Too many field-dependent Hamiltonian terms'
         WRITE (LUPRI,*) 'Current maximum is',MXFELT
         WRITE (LUPRI,*) 'You have specified',NFIELD
         INPERR = INPERR + 1
      END IF
      IF (NFIELD .GT. 0 .AND. SUPSYM) THEN
         SUPSYM = .FALSE.
         NWARN  = NWARN + 1
         WRITE(LUPRI,'(/A)') 'SIRIUS WARNING: .SUPSYM ignored because'//
     &      ' not implemented for finite field calculationss'
      END IF
C
C     New test to avoid conflicting input for subsequent CCSD runs
C
      IF (NASHT .GE. 1 .AND. DOCCSD) THEN
         WRITE (LUPRI,*) 'Only closed-shell RHF wave functions'//
     &        ' can be used for CCSD calculations'
         INPERR = INPERR + 1
      END IF
C
C930406-hjaaj: if (IREST.ne.0) then reset options here and not in
C     sir_geninp as before, where something may not be defined yet
C---------------
      IF (IREST .NE. 0) THEN
         FLAG(13) = .TRUE.
         LMOORD   = .FALSE.
         IF (IREST .EQ. 1) THEN
C           ... restart from SIRIUS.RST
            IRDMO = 9
            ICI0  = 5
C
C           If we restart with DIIS or Roothan Hartree-Fock, then we
C           probably don't have trust radius infomation on file
C           This statement disables the reading of such info
C
            IF ( (MAXFCK .GT. 0 .OR. MXDIIS .GT. 0)
     &           .AND. NASHT .LE. 1 ) ICI0 = 1
         ELSE IF (IREST .EQ. 2) THEN
C           ... start from SIRIUS.RST, perhaps from another geometry
            IRDMO = 9
            IF (.NOT.DOSCF) ICI0 = 4
         ELSE IF (IREST .EQ. 3) THEN
C           ... start from SIRIFC, perhaps from another geometry
            IRDMO = 10
            ICI0 = 10
C           ... TODO hjaaj: implement read CI coeff. from SIRIFC
            CALL QUIT('.IFCSTART is not implemented yet in SIRIUS')
         ELSE
            CALL QUIT('SIRIUS: Unknown start/restart option')
         END IF
      END IF
C
C     Define /INFTRA/
C     Set normal level of integral transformations, unless
C     explicitly set in input.
C
      IF (ITRLVL .EQ. -6789) THEN
         IF (DOCI .AND. .NOT. DOMC) THEN
            ITRLVL = 0
         ELSE IF (FLAG(12)) THEN
C        ... include (ii/aa) and (ia/ia) for orbital diagonal
            ITRLVL = 2
         ELSE
C        ... exclude (ii/aa) and (ia/ia) for orbital diagonal
            ITRLVL = 3
         END IF
C        IF (DONEVPT) ITRLVL=4
C        ... this is instead done just before the KOOPRO4 call in SIRCTL
C            in this way we don't spend CPU time on too big integral
C            transformations in the MCSCF step /hjaaj June 09
      END IF
C
C     Set print parameters in /INFTRA/ for TRACTL
C     (only set IPRTRA if not set explicitly in input)
C
      IF (IPRTRA .EQ. -1) IPRTRA = IPRSIR
C
C     If MAXCIT not specified in input:
      IF (MAXCIT .LT. 0) THEN
         IF (.NOT.(DOSCF .OR. DOMP2 .OR. DOCINO) .AND. FLAG(13)) THEN
C           ... only MC and reading orbitals from file:
C               these may be converged orbitals thus MAXCIT = 20
            MAXCIT = 20
         END IF
      END IF
C     If MAXMAC, MAXJT, MAXMIC not specified in input:
      IF (MAXMAC .LT. 0) THEN
         IF (ISTATE .EQ. 1) THEN
            MAXMAC = 25
         ELSE
            MAXMAC = 40
         END IF
      END IF
      IF (MAXJT  .LT. 0) MAXJT  = 24
      IF (MAXMIC .LT. 0) MAXMIC = MAXJT*MAXMAC
C
C ********  END OF READ OF INPUT DATA  **********
C
C
C *** If requested (flag(13)), get starting orbitals with READMO
C
C hj-aug99: flag(13) is used in sirdiis.F to set initial screening
C   in DIRFCK: loose scr. threshold if H1DIAG or Huckel, tight scr. thr.
C   if orbitals are read in.
C
C
C     if "MOSTART 8" then call cmodel to
C        kdel .eq. 0: define unit matrix
C        kdel .gt. 0: delete 3s in d, etc.
C        kdel .lt. 0: DELETE specified, MOSTART 8 illegal
C        also call DELMO to delete numerically linear dependent
C             orbitals, as defined by THROVL, using canonical
C             orthonormalization. (Flag(33) means .AO DELETE
C             has been specified.)
C
      IF (IRDMO .EQ. 8) THEN
         CALL SETORB
         CALL CMODEL(CMO,.TRUE.)
C        CALL CMODEL(CMO,GETCMO)
C
         CALL SETORB
         KSH1  = 1    + NCMOT
         KSCRA = KSH1 + NNBAST
         LSCRA = LCMO - KSCRA
         CALL DELMO(CMO(1),CMO(KSCRA),LSCRA,THROVL,CMAXMO,DELEMO)
C        ... DELEMO true means orbitals have been deleted.
         IF (DELEMO) THEN
            IF (.NOT.FLAG(33)) THEN
               NWARN = NWARN + 1
               WRITE (LUPRI,'(/A/A)') 'WARNING:'//
     *         ' Orbitals have been deleted because of numerical',
     *         '   linear dependence although not requested in input.'
            END IF
         END IF
C
         CALL SETORB
C
C        Diagonalize H1 for start orbitals
C
         CALL H1MO(CMO,CMO(KSH1),CMO(KSCRA),LSCRA)
C        CALL H1MO(CMO,SH1,SCRA,LSCRA)
C
C        Find Huckel start orbitals if DOHUCKEL and not AUTOCC
C        (If automatic occupation HUCDRV is called later in this routine)
C        Note: H1MO called first anyway, to generate symmetry adapted
C        MOs when supersymmetry exists. /hjaaj
C
         IF ( DOHUCKEL .AND. (.NOT.AUTOCC .OR. NSYM.EQ.1) )
     &      CALL HUCDRV(CMO,CMO(KSCRA),LSCRA)
C
         FLAG(13) = .FALSE.
      ELSE IF (FLAG(13)) THEN
C
C        ... read old orbitals from file
         NERROR = 0
         NDELT = 0
         IF (IRDMO .EQ. 7 .OR. IRDMO .EQ. 9) THEN
C           ... read initial MO's from SIRIUS.RST
C               Find out if any orbitals deleted
            CALL GPOPEN(LUIT1,'SIRIUS.RST','OLD',' ','UNFORMATTED',
     &                  IDUMMY,.FALSE.)
            REWIND LUIT1
            IF (FNDLAB('BASINFO ',LUIT1)) THEN
C           ... 921202: old LUIT1 files do not have BASINFO
C
               IF (INIOCC .LE. 0 .OR. IREST .GT. 0) THEN
                  READ (LUIT1) NSYM_i4,NBAS_i4,NORB_i4,NRHF_i4,IOPRHF_i4
                  IF (INIOCC .LT. 0) THEN
C                 ... 081205 hjaaj: for INIOCC .lt. 0, we
C                     may restart with old orbitals, but new charge.
                     NRHFEL_test = 0
                     DO I = 1, NSYM_i4
                        NRHFEL_test = NRHFEL_test + 2*NRHF_i4(I)
                     END DO
                     IF (IOPRHF_i4 .GT. 0) NRHFEL_test = NRHFEL_test + 1
                     IF (NRHFEL_test .eq. NRHFEL) THEN
                        INIOCC  = 2
                        NRHF(:) = NRHF_i4(:)
                        IOPRHF  = IOPRHF_i4
                     END IF
                  ELSE
                     INIOCC  = 2
                     NRHF(:) = NRHF_i4(:)
                     IOPRHF  = IOPRHF_i4
                  END IF
               ELSE
                  READ (LUIT1) NSYM_i4,NBAS_i4,NORB_i4
               END IF
               MSYM    = NSYM_i4
               MBAS(:) = NBAS_i4(:)
               NORB(:) = NORB_i4(:)
               NERROR  = -1
            END IF
            CALL GPCLOSE(LUIT1,'KEEP')
         ELSE IF (IRDMO .EQ. 10) THEN
C           ... read initial MO's from LUSIFC
C               Find out if any orbitals deleted
            CALL GPOPEN(LUSIFC,FNSIFC,'UNKNOWN',' ','UNFORMATTED',
     &                  IDUMMY,.FALSE.)
            REWIND LUSIFC
            IF (FNDLAB(LBSIFC,LUSIFC)) THEN
               NERROR = -1
               READ (LUSIFC) XPOTNC,EMY,EACTIV,EMCSCF,
     &               JSTATE,JSPIN,MACTEL,LXSYM,MS2X
C              READ (LUSIFC) POTNUC,EMY,EACTIV,EMCSCF,
C    *               ISTATE,ISPIN,NACTEL,LSYM,MS2
             IF (INIOCC .EQ. 0 .AND. MACTEL .LE. 1) THEN
C              ... 081205 hjaaj: not for INIOCC .lt. 0, because we
C                  may restart with old orbitals, but new charge.
               READ (LUSIFC) (IDUM,I=1,13),MSYM,(IDUM,I=1,64),
     &            NRHF,(IDUM,I=1,24),NORB,MBAS
               IF (MACTEL .EQ. 1) THEN
                  IOPRHF = LXSYM
               ELSE
                  IOPRHF = 0
               END IF
               INIOCC = 3
             ELSE
               READ (LUSIFC) (IDUM,I=1,13),MSYM,(IDUM,I=1,96),NORB,MBAS
C              READ (LUSIFC) NISHT,NASHT,NOCCT,NORBT,NBAST,NCONF,NWOPT,
C    *               NWOPH,NCMOT,NNASHX,NNASHY,NNORBT,N2ORBT,
C    *               NSYM,MULD2H, NRHF,NFRO,
C    *               NISH,NASH,NORB,NBAS,
C    *               NELMN1, NELMX1, NELMN3, NELMX3, MCTYPE,
C    *               NAS1, NAS2, NAS3
             END IF
C            IF (XPOTNC .NE. POTNUC) INPERR = INPERR + 1
C            IF (XPOTNC .NE. POTNUC .AND. IREST .NE. 0) IREST = -1
C              ... we cannot use this test because we
C                  may have a new geometry (geometry optimization)
C                  However, code info for SIR_PRTINP in IREST if restarted
C                  calculation /hjaaj feb. 2003
            ELSE
               INPERR = INPERR + 1
               WRITE (LUPRI,'(/2A/3A)')
     &         ' Input ERROR, MO coefficients to be read from'//
     &         ' SIRIUS interface file ',FNSIFC,' but ',
     &         FNSIFC,' with proper information was not provided.'
            END IF
            CALL GPCLOSE(LUSIFC,'KEEP')
         END IF
         IF (NERROR .EQ. -1) THEN
C           test NSYM,NBAS information from LUIT1 or LUSIFC,
C           are they for same basisset and same no. of symmetries ?
            NERROR = 0
            IF (MSYM .NE. NSYM) THEN
               NERROR = NERROR + 1
            ELSE
               DO 40 ISYM = 1,NSYM
                  IF (NBAS(ISYM) .NE. MBAS(ISYM)) NERROR = NERROR + 1
                  NDELT = NDELT + (NBAS(ISYM) - NORB(ISYM))
   40          CONTINUE
            END IF
         END IF
         IF (NERROR .GT. 0) THEN
            INPERR = INPERR + 1
            WRITE (LUPRI,'(/A,I5/5X,A,I3,5X,8I5)')
     &        'Input ERROR in connection with ".MOSTART" of type',IRDMO,
     &        'NSYM, NBAS from input       :',NSYM,(NBAS(I),I=1,NSYM)
            WRITE (LUPRI,'(5X,A,I3,5X,8I5)')
     &        'NSYM, NBAS from MOSTART file:',MSYM,(MBAS(I),I=1,MSYM)
         ELSE IF (NDELT .NE. 0) THEN
            NINFO = NINFO + 1
            WRITE (LUPRI,'(/5X,A/I10,A)')
     &      'INFO: Information regarding orbitals read from file:',
     &      NDELT,' basis functions are deleted'
         END IF
      END IF
C
C *** set up orbital indexing arrays
C
      CALL SETORB
C
C *** Check for input errors for core hole calculation
C
      IF (CORHOL) THEN
         IF (DOSCF .AND. AUTOCC .AND. NSYM .GT. 1) THEN
            INPERR = INPERR + 1
            WRITE (LUPRI,'(/A/A)')
     &     'ERROR: Input error for core-hole calculation:',
     &     'ERROR:  - automatic occupation of SCF orbitals not allowed.'
         END IF
         IF (IRDMO .EQ. 8) THEN
            INPERR = INPERR + 1
            WRITE (LUPRI,'(/A/A)')
     &      'ERROR: Input error for core-hole calculation:',
     &      'ERROR:  - no initial molecular orbitals from'//
     &         ' a previous calculation specified.'
         END IF
         IF (DOSCF .AND. DOMC) THEN
            INPERR = INPERR + 1
            WRITE (LUPRI,'(/A/A)')
     &      'ERROR: Input error for core-hole calculation:',
     &      'ERROR: - both Hartree-Fock and MCSCF specified.'
         END IF
         IF (NASHT .GT. 1) THEN
          IF (NELMN1 .LE. 0 .AND. NELMX1 .LE. 0) THEN
            NELMN1 = 1
            NELMX1 = 1
          END IF
          IF (NELMN1 .NE. NELMX1 .OR. MCTYPE .NE. 2 .OR.
     &        NAS1(JCHSYM) .NE. 1 .OR. NAS1T .NE. 1) THEN
            INPERR = INPERR + 1
            WRITE (LUPRI,'(/A/A)')
     &       ' Input ERROR(s) for core-hole calculation:',
     &       ' -----------------------------------------'
           IF (MCTYPE .NE. 2)
     &      WRITE (LUPRI,*) 'Wave function type is not RAS'
           IF (NELMN1 .NE. NELMX1) THEN
            WRITE (LUPRI,*) 'Number of electrons in RAS1 is not fixed'
            WRITE (LUPRI,*) '-- minimum no. of electrons : ',NELMN1
            WRITE (LUPRI,*) '-- maximum no. of electrons : ',NELMX1
           END IF
           IF (NAS1T .NE. 1)
     &        WRITE (LUPRI,*) 'No. of RAS1 orbitals different from one'
           IF (NAS1(JCHSYM) .NE. 1) WRITE (LUPRI,*)
     & 'Not exactly one RAS1 orbital in core orbital symmetry : ',JCHSYM
           WRITE (LUPRI,*)
          END IF
         ELSE IF (NASHT .EQ. 1) THEN
            IF (JCHSYM .NE. IOPRHF) THEN
               INPERR = INPERR + 1
               WRITE (LUPRI,*) 'Input ERROR: '//
     &         'core-hole sym. and RHF open shell sym. different'
            END IF
         ELSE
            INPERR = INPERR + 1
            WRITE (LUPRI,*) 'Input ERROR: '//
     &      'core-hole not compatible with closed-shell RHF'
         END IF
         IF (CORRLX) THEN
C           ... no reordering of orbitals :
            LMOORD   = .FALSE.
C           ... NEO ALWAYS false :
            FLAG(38) = .FALSE.
C           ... NR ALWAYS true :
            FLAG(39) = .TRUE.
C           ... NO ABSORPTION :
            FLAG(51) = .FALSE.
            FLAG(52) = .FALSE.
            FLAG(53) = .FALSE.
         ELSE
            IF (JCHORB .GT. NISH(JCHSYM)+1) THEN
               INPERR = INPERR + 1
               WRITE (LUPRI,'(/A,I4,A,I2,A)')
     &            'ERROR: Input error for core-hole calculation:'//
     &            'ERROR: Core hole orbital',JCHORB,' in symmetry',
     &               JCHSYM,' is neither inactive nor in RAS1.'
            ELSE IF (JCHORB .LE. NISH(JCHSYM)) THEN
               IF (LMOORD) THEN
                  INPERR = INPERR + 1
                  WRITE (LUPRI,'(/A/A)')
     &            'ERROR: Input error for core-hole calculation:',
     &            'ERROR: .COREHOLE and .REORDER are not compatible'//
     &               ' if core hole orbital is inactive.'
               END IF
               LMOORD = .TRUE.
               DO 120 I = 1,NORBT
  120             IMOORD(I) = I
               DO 130 I = IORB(JCHSYM)+JCHORB,IORB(JCHSYM)+NISH(JCHSYM)
  130             IMOORD(I) = I+1
               IMOORD(IORB(JCHSYM)+NISH(JCHSYM)+1) = IORB(JCHSYM)+JCHORB
            END IF
            LNOROT = .TRUE.
            NOROT(IORB(JCHSYM)+NISH(JCHSYM)+1) = 1
         END IF
      END IF
C
C *** Read initial molecular orbitals (if IRDMO .eq. 8, they
C     are already available).
C
      IF (FLAG(13) .AND. INPERR .EQ. 0) THEN
         K1   = 1    + NCMOT
         LSCR = LCMO - K1
         JRDMO = -IRDMO
         IPRSIR_save = IPRSIR
         IPRSIR = max(1, IPRSIR)
         ! .ge. 1 to print info about orbitals on SIRIUS.RST when
         ! reading initial orbitals from SIRIUS.RST
         CALL READMO(CMO,JRDMO)
         IPRSIR = IPRSIR_save
C        JRDMO negative : return with error status in JRDMO
C        JRDMO positive : abort in READMO if error
         IF (JRDMO .GT. 0) THEN
            WRITE(LUPRI,'(/A)')
     &         ' SIRRDI input ERROR, not able to read mo coefficients.'
            INPERR = INPERR + 1
            GO TO 200
         END IF
      END IF
  200 CONTINUE
C
C *** If requested, reorder mo:s, using reorder array in IMOORD(*).
C
      IF (LMOORD .AND. INPERR .EQ. 0) THEN
         K1 = 1 + NCMOT
         CALL REORD(CMO,CMO(K1),IMOORD)
      END IF
C
C *** AUTOCC : find initial occupation (hjaaj-aug95)
C              and Huckel start orbitals if DOHUCKEL
C
      IF (DOSCF .AND. AUTOCC .AND. IREST .EQ. 0
     &   .AND. INPERR .EQ. 0) THEN
C
C        Initial occupation may have been specified in input or read
C        from LUIT1 or LUSIFC. If not, and more than one symmetry,
C        we find it now from diagonal elements of the H1 matrix.
C        If DOHUCKEL use Huckel guess (INIOCC = 5).
C
         KSCRA = 1    + NCMOT
         LSCRA = LCMO - KSCRA
         IF (INIOCC .LE. 0) THEN
            IF (DOHUCKEL) THEN
               CALL HUCDRV(CMO,CMO(KSCRA),LSCRA)
               FLAG(13) = .FALSE.
               INIOCC = 5
            ELSE
               CALL H1OCC(CMO,CMO(KSCRA),1,LSCRA)
               INIOCC = 4
            END IF
         ELSE IF (INIOCC .EQ. 1 .AND. DOHUCKEL) THEN
            AUTOCC = .FALSE.
C           ... otherwise HUCDRV will destroy user input
C               of initial occupation guess
            CALL HUCDRV(CMO,CMO(KSCRA),LSCRA)
            FLAG(13) = .FALSE.
            AUTOCC = .TRUE.
         END IF
      END IF
      IF (NSYM .EQ. 1) AUTOCC = .FALSE.
C     ... occupation cannot change if only one symmetry
C
C     Check for input errors of initial occupation
C
      NERROR = 0
      DO ISYM = 1,NSYM
         IF (NRHF(ISYM) .GT. NORB(ISYM)) NERROR = NERROR + 1
      END DO
      INPERR = INPERR + NERROR
      IF (NERROR .GT. 0) THEN
         WRITE(LUPRI,'(/A)') 'INPUT ERROR - RHF occupation too big'
      END IF
C
C
C     If INIOCC .ge. 2 and RHF/MP2 highest level (i.e. MCTYPE .LE. 0),
C     we need to retransfer SCF occupation to MCSCF parameters as
C     NRHF(*),IOPRHF may have changed since previous RHF setting above.
C     (INIOCC .ge. 2: NRHF(*) from LUIT1, LUSIFC, H1OCC, or HUCDRV)
C
      IF (INIOCC .GE. 2 .AND. MCTYPE.LE.0) THEN
         CALL ICOPY(8,NRHF,1,NISH,1)
         CALL IZERO(NASH,8)
         IF (IOPRHF .NE. 0) THEN
            LSYM   = IOPRHF
            ISPIN  = 2
            MS2    = 1
            NASH(LSYM) = 1
            NASHT  = 1
            NACTEL = 1
         ELSE IF (HSROHF) THEN
            LSYM = 1
            DO ISYM = 2, NSYM
               IF (MOD(NROHF(ISYM),2).EQ.1) LSYM = MULD2H(LSYM,ISYM)
            END DO
            NASHT  = ISUM(NSYM,NROHF,1)
            ISPIN  = NASHT + 1
            MS2    = NASHT
            NACTEL = NASHT
            CALL ICOPY(8,NROHF,1,NASH,1)
         ELSE
            LSYM   = 1
            ISPIN  = 1
            MS2    = 0
            NASHT  = 0
            NACTEL = 0
         END IF
         CALL SETORB
      END IF

C     To be able to use spin-dependent functionals for closed-shell Hartree-Fock
      IF (DFT_SPINDNS .AND. NASHT .EQ. 0) THEN
         DFT_SPINDNS = .FALSE.
      END IF

C
C *** Check if initial occupation for RHF has been defined.
C
      IF ((DOSCF .OR. DOMP2 .OR. DOCCSD) .AND. INIOCC .LE. 0) THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,'(/A/A)')
     &     'ERROR: Input error for DFT, HF, MP2, or CC calculation:',
     &     'ERROR: Initial SCF occupation has not been defined.'
         INIOCC = 0
C        ... set INIOCC = 0 for SIR_PRTINP
      END IF
C
C CCSD only for closed shell RHF
C
      IF (DOCCSD .AND. NASHT .GE. 1) THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,'(/A)') ' Input error: No open shells '//
     &        'allowed with CC calculation presently'
      END IF
C
C     ***** Set convergence thresholds, if not set by user
C
C If (MP2) force HF convergence to THRRHF_MP2
C If (MP2) force canonical Hartree-Fock orbitals
C hjaaj Aug 04: Also if (CC)
C
      IF (DOMP2 .OR. DOCCSD) THEN
        IF (DOCCSD. and. DOSCF .AND. (IOPRHF .GT. 0 .OR. HSROHF)) THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,'(/A)') ' Input ERROR: closed shell CC'//
     &      ' is not compatible with open shell RHF'
        END IF
        IF (DOHFSRDFT) THEN
C hjaaj May 2008 : THRRHF_MP2 too small for srMP2  NB NOTE BEMAERK TODO OBS
           THRRHF_MAX = 50.0D0*THRRHF_MP2
        ELSE
           THRRHF_MAX = THRRHF_MP2
        END IF
        IF (THRRHF .LE. D0) THEN
           THRRHF = THRRHF_MAX
        ELSE IF (THRRHF .GT. THRRHF_MAX) THEN
          NWARN = NWARN + 1
          WRITE (LUPRI,'(/5X,2(A,1P,D10.2)/5X,A)')
     &   'SIRRDI WARNING: .THRRHF =',THRRHF,
     &   ' is greater than ',THRRHF_MAX,
     &   'SIRRDI WARNING: This may cause problems '//
     &   ' because MP2/CC requires converged orbitals'
        END IF
        IF (.NOT.RHFCAN) THEN
         NINFO = NINFO + 1
         WRITE (LUPRI,'(/5X,A)')
     &  'SIRRDI INFO, ".NONCANONICAL" ignored because MP2/CC requires'//
     &  ' canonical orbitals.'
         RHFCAN = .TRUE.
        END IF
      END IF
      IF (THRRHF .LE. D0) THRRHF = THR_DEFAULT
      IF (THRCI  .LE. D0) THRCI  = THR_DEFAULT
      IF (THRMC  .LE. D0) THRMC  = THR_DEFAULT
C
C     set THRGRD for test of value below
C     set THQKVA if not specified in input
C
      IF (DOMC) THEN
         THRGRD = THRMC
         IF (THQKVA .LE. D0) THQKVA = 8.0D0
      ELSE IF (DOCI) THEN
         THRGRD = THRCI
         IF (NASHT .GT. 0) THEN
            IF (DOCISD) THEN
               INPERR = INPERR + 1
               WRITE(LUPRI,'(/A)')
     &         'INPUT ERROR: both .CISD   and active spaces specifed.'
            END IF
            IF (DOFCI ) THEN
               INPERR = INPERR + 1
               WRITE(LUPRI,'(/A)')
     &         'INPUT ERROR: both .FULLCI and active spaces specifed.'
            END IF
         ELSE IF (DOCISD) THEN
            ! MCTYPE = 2 ! RAS (set in CI INPUT)
            NISHT = 0
            NAS1T = 0
            DO I = 1,NSYM
               NISH(I) = NFRO(I)
               NAS1(I) = NRHF(I) - NISH(I)
               NAS2(I) = 0
               NAS3(I) = NORB(I) - NAS1(I)
               NASH(I) = NAS1(I) + NAS3(I)
               NISHT   = NISHT   + NISH(I)
               NAS1T   = NAS1T   + NAS1(I)
            END DO
            NELMN1 = 2*NAS1T - 2
            NELMX1 = 2*NAS1T
            NELMN3 = 0
            NELMX3 = 2
            NACTEL = ICHARG - KCHARG - 2*NISHT
            IF (LSYM .LT. 0) LSYM = 1
         ELSE IF (DOFCI) THEN
            ! MCTYPE = 1 ! CAS (set in CI INPUT)
            NISHT = 0
            DO I = 1,NSYM
               NISH(I) = NFRO(I)
               NAS1(I) = 0
               NAS2(I) = NORB(I) - NISH(I)
               NAS3(I) = 0
               NASH(I) = NAS2(I)
               NISHT   = NISHT   + NISH(I)
            END DO
            NELMN1 = 0
            NELMX1 = 0
            NELMN3 = 0
            NELMX3 = 0
            NACTEL = ICHARG - KCHARG - 2*NISHT
            IF (LSYM .LT. 0) LSYM = 1
         END IF
      ELSE IF (DOSCF) THEN
         THRGRD = THRRHF
         IF (THQKVA .LE. D0) THQKVA = 0.8D0
      END IF

      IF (NSYM .EQ. 1 .AND. LSYM .LT. 0) LSYM = 1 ! if no symmetry, then do not require .SYMMETRY in input

      IF (LSYM .LT. 0) THEN
         WRITE(LUPRI,'(///A/)') 'FATAL INPUT ERROR: '//
     &      '.SYMMETRY must be specified for MCSCF/CI calculations.'
         CALL QUIT( 'FATAL INPUT ERROR: '//
     &      '.SYMMETRY must be specified for MCSCF/CI calculations.')
      END IF

      IF (MS2 .EQ. -12345) THEN ! MS2 = 2*MS
         IF (MULTS .LE. 0) MULTS = MOD(NACTEL,2) + 1
         IF (FLAG(27) .OR. DOMCSRDFT) THEN ! determinants or open-shell MC-srDFT
            ! MC-srDFT: we need the high-spin MS value for the srDFT part
            MS2 = MULTS - 1 ! MS = max value = S, to minimize risk of converging to wrong spin state
         ELSE
            MS2 = MOD(NACTEL,2) ! MS = 0 or 1/2, for most efficient code
         END IF
      ELSE
         IF (MULTS .LE. 0) MULTS = abs(MS2) + 1
         IF (FLAG(27) .AND. MULTS.ne.abs(MS2)+1) THEN ! determinants
            WRITE(LUPRI,'(//A,I0,A,I0,A)')
     &      ' INPUT ERROR: specified .MS2 and .SPIN MULTIPLICITY (',
     &      MS2,' and ',MULTS,') are not compatible for determinants.'
            CALL QUIT('INPUT ERROR: .MS2 and .SPIN M not compatible.')
         END IF
      END IF
      ISPIN = MULTS
C     .. ISPIN in infinp.h, MULTS in spinfo.h for CI module
C
C     Call AVESET to set up super symmetry information
C     (If SUPSYM false, then super symmetry information
C      will be a copy of the symmetry information)
C
      IF (INPERR .GT. 0 .AND. SUPSYM) THEN
         WRITE(LUPRI,'(/A)') ' NOTE: super symmetry analysis'//
     &   ' not possible because of input errors.'
         SUPSYM = .FALSE.
      END IF
      IF (NFIELD .GT. 0) THEN
         WRITE(LUPRI,'(/A)') ' NOTE: super symmetry analysis'//
     &   ' not possible because of finite field.'
         SUPSYM = .FALSE.
      END IF

      K1    = 1 + NCMOT
      KFREE = 1
      LFREE = LCMO - K1
      CALL AVESET(CMO,CMO(K1),KFREE,LFREE)
C     CALL AVESET(CMO,WRK,KFREE,LFREE)
C
C
C *** Set up various index arrays etc.
C     If SUPSYM check if OK in AVECHK
C
      IF (IORTO .NE. 1) THEN
C     ... only SETORB needed for reorthonormalization and
C         SETORB has already been called.
         K1 = 1 + NCMOT
         JWOPSY = 1
         CALL SIRSET(CMO(K1),LCMO-K1,.FALSE.)
C        CALL SIRSET(WRK,LWRK,OLDWOP)
         CALL AVECHK(INPERR)
      END IF
C
C  Do we use embedding ?
C  (that is, is QM part embedded in environment)
C
      IF (QMMM .AND. ECP) THEN
         CALL QUIT('QMMM is not compatible with the use of ECPs')
      END IF
      IF (QMMM .AND. DOCI) THEN
         CALL QUIT('QMMM is currently not compatible with .CI')
      END IF
      IF (USE_PELIB() .AND. DOCI) THEN
         CALL QUIT('.PEQM is currently not compatible with .CI')
      END IF
      EMBEDDING = FLAG(16) .or. PCM .or. QM3 .or. QMMM .OR. QMNPMM
     &            .OR. USE_PELIB()
C
C ***  PRINT INPUT DATA *****
C
      IF (IPRSIR.GT.2 .OR. IPRI4.GT.7 .OR. IPRI6.GT.2)
     &   SIR_INPPRC = .FALSE.
C     ... force output if print level increased by user
      CALL SIR_PRTINP(LUPRI,NUMRUN,IRDMO)
      IF (IREST.NE.0) WRITE (LUPRI,3050)
      IF (IREST.EQ.-1) WRITE (LUPRI,3052)
      CALL FLSHFO(LUPRI)
      IF (LUW4 .NE. LUPRI .AND. P4FLAG(1)) THEN
         CALL SIR_PRTINP(LUW4,NUMRUN,IRDMO)
         IF (IREST.NE.0) WRITE (LUW4,3050)
         IF (IREST.EQ.-1) WRITE (LUW4,3052)
         CALL FLSHFO(LUW4)
      END IF
      IF (LUCME .GT. 0) THEN
         CALL SIR_PRTINP(LUCME,NUMRUN,IRDMO)
         IF (IREST.NE.0) WRITE (LUCME,3050)
         IF (IREST.EQ.-1) WRITE (LUCME,3052)
         CALL FLSHFO(LUCME)
      END IF
 3050 FORMAT(/5X,'This is a restarted calculation.')
 3052 FORMAT( 5X,'(Probably from another geometry:'//
     &       ' nuclear repulsion energy has changed.)')
C
C
C *****  PROCESS AND CHECK INPUT DATA *****
C
      IF(NSYM  .LE.0 .OR.NSYM  .GT.8)      CALL ERRINP( 4,INPERR)
      IF(NORBT .LT.0 .OR.NORBT .GT.MXCORB) CALL ERRINP( 8,INPERR)
      IF (IORTO .EQ. 1) GO TO 8000
C
      IF (DOSCF .OR. DOCI .OR. DOMC) THEN
         IF(MAXMAC.LE.0 .OR.MAXMAC.GT.200) CALL ERRINP( 5,INPERR)
         IF(THRGRD.LT.D0.OR.THRGRD.GT.D1)  CALL ERRINP( 7,INPERR)
      END IF
      IF(NRHFT .LT.0 .OR.NRHFT .GT.MAXOCC) CALL ERRINP( 9,INPERR)
      IF(NASHT .LT.0 .OR.NASHT .GT.MAXASH) CALL ERRINP(10,INPERR)
!     IF(DOSCF .AND. AUTOCC .AND. MCTYPE .GT. 0 .AND. NSYM .GT. 1) THEN
!        INPERR = INPERR + 1
!        WRITE (LUPRI,'(/A)') 'ERROR: HF occupation must be specified'//
!    *         ' when SCF followed by CI or MCSCF'
!     END IF
      IF(NOCCT .LT.0 .OR. NOCCT .GT.MAXOCC) THEN
         INPERR = INPERR + 1
         WRITE(LUPRI,'(/A,I5/A,I5/A)')
     *      'ERROR: Number of occ. orbitals :',
     *      NOCCT,' is .lt. 0 or .gt. MAXOCC :',MAXOCC,
     *      'Increase MAXOCC in DALTON/include/maxorb.h and recompile'
      END IF
      IF(NSSHT .LT.0 .OR.NSSHT .GT.MXCORB) CALL ERRINP(11,INPERR)
C
      NERROR = 0
      DO 310 ISYM = 1,NSYM
         IF (NFRO(ISYM).GT.NISH(ISYM)) NERROR = NERROR + 1
  310 CONTINUE
      IF (NERROR.GT.0) THEN
         INPERR = INPERR + NERROR
         WRITE (LUPRI,3101) NERROR
      END IF
 3101 FORMAT(/' SIRRDI input ERROR, NFRO is larger than NISH for',
     *        I2,' symmetries.',
     *       /' (Only inactive orbitals can be frozen with ".FROZEN",',
     *        ' use ".FREEZE" instead)')
C
C Population analysis (flag(6)) requires natural/canonical orbitals
C
      IF (FLAG(6)) THEN
         RHFCAN = .TRUE.
         ICICNO = 1
         IMCCNO = 1
      END IF
C
C If (SOLVNT) check epsilon, cavity radius, and lmax
C
      IF (FLAG(16)) THEN
         NERROR = 0
         RSOLMN = MIN(RSOL(1),RSOL(2),RSOL(3))
         IF (EPSOL .LE. D1) NERROR = NERROR + 1
         IF (EPSTAT.LE. D1) NERROR = NERROR + 1
         IF (RSOLMN.LT. D0) NERROR = NERROR + 1
         IF (LSOLMX.LT. 0)  NERROR = NERROR + 1
         IF (RSOL(1) .NE. RSOL(2) .OR. RSOL(1) .NE. RSOL(3))
     &      NERROR = NERROR + 1
         IF (NERROR .GT. 0) THEN
            INPERR = INPERR + 1
            WRITE (LUPRI,'(/A,/A,I12,2(/A,F12.6)/A,3F12.6)')
     &     ' Input ERROR in one or more reaction field specifications:',
     &     ' Lmax                =',LSOLMX,
     &     ' Dielectric constant =',EPSOL,
     &     ' Static diel. const. =',EPSTAT,
     &     ' Cavity radius       =',RSOL(1),RSOL(2),RSOL(3)
         END IF
C
C        orbital absorption not safe
C
         IF (FLAG(51) .OR. FLAG(52) .OR. FLAG(53)) THEN
            WRITE (LUPRI,'(/A/)') ' STOP because ABSORPTION is not'//
     &         ' yet implemented for SOLVENT calculations'
            INPERR = INPERR + 1
         END IF
      END IF
C
C Check IMOORD(*) if orbitals are to be reordered.
C
      IF (LMOORD) THEN
         NP1=0
         NP2=0
         NERROR = 0
         DO 320 I=1,NSYM
            J1 = IORB(I)+1
            J2 = IORB(I)+NORB(I)
            DO 315 J = J1,J2
               NP1 = NP1 + J
               NP2 = NP2 + IMOORD(J)
  315       CONTINUE
            IF (NP1.NE.NP2) NERROR = NERROR + 1
  320    CONTINUE
         IF (NERROR.GT.0) CALL ERRINP(12,INPERR)
      END IF
C
      IF (FLAG(2) .AND. .NOT. FLAG(4)) THEN
C        ... if not ci calculation then
         IF (NROOTS.LT.1.OR.NROOTS.GT.MAXRTS)
     *                                CALL ERRINP(13,INPERR)
         DO 36 I=1,NROOTS
            IF (IROOT(I).LE.0)        CALL ERRINP(14,INPERR)
            DO 36 J=1,(I-1)
               IF (IROOT(J).EQ.IROOT(I))
     *                                CALL ERRINP(14,INPERR)
   36    CONTINUE
      END IF
C
      IF (NWOPT.GT.MAXWOP)            CALL ERRINP(16,INPERR)
C
      IF (IRDMO .EQ. 8 .AND. KDEL .LT. 0) THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,*)
     *      'Sirius input ERROR: DELETE is illegal with MOSTART 8'
      END IF
C
      NERROR = 0
      DO 370 J = 1,NSYM
         IF (NSSH(J).LT.0) NERROR = NERROR + 1
  370 CONTINUE
      INPERR = INPERR + NERROR
      IF (NERROR.NE.0) WRITE (LUPRI,3700) (NSSH(J),J=1,NSYM)
 3700 FORMAT(/' SIRRDI-ERROR, ao integral file not consistent with ',
     *        'input specification'
     *       /T16,'NSSH(*) = NBAS(*) - NOCC(*):',8I5)
C
C
      IF (NSYM.NE.1 .AND. NSYM.NE.2 .AND.
     *    NSYM.NE.4 .AND. NSYM.NE.8) THEN
         WRITE (LUPRI,'(/A)')' SIRRDI, input ERROR.'//
     *      ' Wrong number of symmetries,'//
     *      ' this cannot be a subgroup of D2h.'
         INPERR = INPERR + 1
      END IF
C
      IF (MOD((NACTEL+1),2) .NE. MOD(ISPIN,2) .OR.
     *   (NACTEL+1) .LT. ISPIN) THEN
         WRITE (LUPRI,3710)
         INPERR = INPERR + 1
      END IF
 3710 FORMAT(/' SIRRDI, input ERROR. Number of active electrons',
     *        ' is not consistent with spin multiplicity.')
C
C *****  set some non-specified variables
C
      IF (ICHECK .EQ. -100) THEN
         IF (NROOTS .EQ. 1) THEN
            ICHECK = -1
         ELSE
            ICHECK = 2
         END IF
      END IF
C
      IF (MCTYPE .GT. 0 .AND. PARCAL) THEN
         ! IF (RNRESP) CALL PARQUIT('MCSCF in RESPONSE')
         ! IF (RNABAC) CALL PARQUIT('MCSCF in ABACUS')
      END IF

C
C *****  END OF INPUT PROCESSING *****
C
 8000 IF (INPERR .NE. 0) GO TO 9010
      IF (FLAG(1)) CALL QUIT('SIRIUS "STOP AFTER INPUT" STOP')
C
C
C *** Open other files
C
C
      CALL SIR_OPEN
C
C
      NEWSYM = .FALSE.
      SIR_INPPRC = .TRUE.
      CALL QEXIT('SIRRDI')
      RETURN
C
C *****  ERROR MESSAGES AND ERROR STOP *****
C
 9010 WRITE (LUW4,9000) INPERR
      IF (LUERR.NE.LUW4) WRITE(LUERR,9000) INPERR
      IF (LUPRI.NE.LUW4) WRITE(LUPRI,9000) INPERR
 9000 FORMAT(//20X,I5,' ERROR(s) occured in SIRIUS input processing.'
     *        /20X,'Program cannot continue.'
     *        /20X,'Correct ERROR(s) and be welcome back')
      CALL QTRACE(LUPRI)
      CALL QUIT('*** ERROR *** SIRIUS input errors.')
C
C     END OF SIRRDI.
C
      END
C  /* Deck sir_newinp */
      SUBROUTINE SIR_NEWINP(INPERR,NUMRUN,IRDMO,IREST,NSYM,NORB)
C
C Revisions
C    1-OCT-1987 HA     Keyword directed input for all parameters
C
C OBJECTIVE: Read input data from file LUCMD
C            using new standard Sirius input
C
C
C Modified the routine to include Coupled-Cluster input routines.
C In addition SIR_GENINP is now always called, as it is required.
C *GENINP no longer included in input file
C
C Input:
C  none
C
C Output:
C  IREST .ne. 0 if restart
C
C Called from SIRRDI
C
C
C Calls the following input modules:
C
C  *GENERAL
C  *HAMILTONIAN
C  *OPTIMIZATION
C  *SCF INPUT               used to be *HF INPUT
C  *MP2 INPUT               used to be *MP2 CALCULATION
C  *CONFIGURATION INPUT     used to be *WAVE FUNCTION
C  *CI INPUT                used to be *CI CALCULATION
C  *CC INPUT                used to be *CC
C  *NEVPT2 INPUT
C  *ORBITAL INFORMATION
C  *CI VECTOR
C  *STEP CONTROL
C  *POPULATION ANALYSIS
C  *TRANSFORMATION
C  *RESPONSE CALCULATION
C  *UPDATE
C  *SOLVENT CALCULATION
C  *PRINT LEVELS
C  *AUXILIARY INPUT
C  *LUCITA CALCULATION     ! large-scale CI
C  *APSG INPUT
C  *STEX INPUT
C  *END OF (marks end of keyword input)
C
C******************************************
C List of internal control flags for SIRIUS
C -- saved in FLAG(1:NFLAG) in infinp.h --
C******************************************
C
C Updated: 1-Oct-1987 ha
C
C flag unit    description                                    INPUT
C ---- ----    -----------                                    -----
C                                                           *INP
C 1  SIRCTL  Input test (stop after processing input)      GEN : STOP
C 2  SIRCTL  MCSCF optimization (call SIROPT)              GEN : MCSCF
C 3  SIROUT  Canonical and natural orbitals, if converged  PRI : CANONI
C 4  SIRCTL  CI calculation (CALL SIRCI, no orbital opt.)  CIC : general
C 5  SIRCTL  Optimization result summary (CALL SIROUT)     PRI : SUMMARY
C 6  SIRCTL  Population analysis (CALL SIRPOP)             POP : general
C 7
C 8          (was CALL RSPDRV, now controlled by DORSP variable)
C 9          (was CALL MP2CTL, now controlled by DOMP2 variable)
C 10
C 11 SIRSET  Calculate non-symmetric orbital rotations
C            and write them to LUINDX (labeled)
C 12 GRAD    Exact orbital Hessian diagonal                WAV : FOCKDI
C            (otherwise only Fock matrix contributions)
C 13 various CALL READMO for initial orbitals              ORB : vario
C 14 various Integral transformation not needed.           TRA : OLD TR
C 15 SIRSAV  transformation to NO's each macro iteration.  OPT : TRACI
C 16 SIROPT  Spherical cavity solvent calculation          SOL : genera
C 17 various SUPSYM : find and use 'super symmetry'        OPT : SUPSYM
C 18 SIROUT  Punch MO:s (canonical if flag 3) to LUPNCH    ORB : PUNCHO
C 19 SIROPT  Use UPDATE from the beginning /990821-hjaaj test option/
C 20 SIROPT  No MC optimization, just calculate energy and GEN : STOP
C            gradient.
C 21 various Restricted Hartree-Fock optimization (QC-SCF)
C 22 SRDIIS  use C2-DIIS algorithm                         RHF : C2DIIS
C 23 various Active-active rotations included              OPT : ACTROT
C 24 various No Active-active rotations in RAS             OPT : NO ACT
C 25 var.OPT Write MC information to ABACUS on unit LUSIFC OPT : ABACUS
C                  when MC is converged
C 26 various string information on DISK (for WESTA)        OPT : WESTA
C 27 var.OPT use determinant basis not CSF basis           OPT : DETERM
C 28
C 29
C 30
C 31 SIRINP  Punch orbitals to LUPNCH                      ORB : PUNCHI
C 32 SIRINP  NEG: Gram-Schmidt; POS: sym. orbital orthog.  ORB ORTHGR..
C 33 SIRINP  Delete orbitals based on AO overlap matrix    ORB : AODEL
C 34 intern. MO integrals not needed in MCSCF optimization
C 35 SIRSTP  tight step control for ground state also      STP : TIGHT
C 36 SIRNEO  Dont use additional termination tests.        STP : NO TER
C 37 SIROPT  Dont switch to Broydn update in local region. UPD : genera
C 38 SIROPT  Dont switch from NEO to NR in local region.   OPT : NEO AL
C 39 SIROPT  Only NR local optimization, no NEO steps.     OPT : NR ALW
C 40
C 41 various OPT Optimal orbital trial vector (CALL NEXKAP)OPT : OPTIMA
C 42
C 43
C 44
C 45 SIRCAN  Redo CI with natural orbitals (CALL SIRCI)    OPT : CI CAN
C 46 SIRCAN  Only NO transformation, inactive and          OPT : NATONL
C    OPTST,SIRSAV sec. orbitals unchanged (no can. orbs)
C 47 SIRCAN  temp. - use CICTL, not TRACI)
C 48 SIRCAN,OPTST,SIRSAV CO's: diag. FD=FC+FV blocks       OPT : FOCKON
C 49 SIRDIIS do not use difference density                 DIIS: DIFDEN
C 50
C 51 various Orb. absor. level 1;   inact.-act.,act.-act.  OPT : ABSORB
C 52 various Orb. absor. level 2; + inact.-sec.            OPT : ABSORB
C 53 various Orb. absor. level 3; +   act.-sec. (i.e. all) OPT : ABSORB
C 54 SIROPT  Do not switch absorption off in local region. OPT : ALWAYS
C 55 NEO,NR  include g(orb) H(diag)  ??
C 56 var.UPD Fixed Hessian                                 UPD : FIXED
C 57 UPDINI  Inverse Hess, diag. is initial approx. to     UPD : INVER
C            inverse Hess
C 58 CIST    Take plus combination of degenerate roots     CIV : PLUS C
C 59
C 60
C 61 SIRNEO  Construct E[2] explicitly and quit
C 62
C 63
C 64
C 65
C 66 CISIG1  temp. dont use ms=0 sym. in determinant cisig CIC : EQUALT
C 67 SIRCI   calculate density matrices                    CIC : CIDENS
C 68 SETDET
C 69 SIRCI   force H2AC integrals on disk                  CIC : DISKH2
C 70
C 71 SIRPOP  Do (not) calculate dipole moment              POP : DIPMOM
C 72 SIRPOP  Do (not) calculate quadrupole moments         POP : QUADRP
C 73 SIRPOP  Do (not) perform Mulliken population analysis.POP : MULLIK
C 74 SIRPOP  Do (not) perform Virial analysis.             POP : VIRIAL
C 75
C 76 CINEX   WERESD-energy weighted residuals              CIC : WEIGHT
C 77 CINEX   zero small elements in CI trial vectors       CIC : ZEROEL
C 78 CINEX   symmetrize CI vector                          CIC : CISYMM
C 79
C 80 OTRIAL=TRUE; use orbital trial vectors each macro     OPT : ORB_TR
C
C******************************************
C
      use lucita_cfg
      use lucita_orbital_spaces
#include "implicit.h"
C
C -- local variables and constants
C
      PARAMETER (NTBMOD = 27)
      CHARACTER WORD*7, PROMPT*1, TABMOD(NTBMOD)*7
      LOGICAL   ALLOPT, PTABMO, lucita_input_provided
      LOGICAL   ensemble_input_provided
C needed for lucita: nsym, norb
      DIMENSION NORB(*)
C
#include "priunit.h"
#include "gnrinf.h"
C
      CALL QENTER('SIR_NEWINP')
C
      DATA TABMOD/'*END OF','*CONFIG','*OPTIMI','*CI INP','*ORBITA',   !  1: 5
     *            '*CI VEC','*STEP C','*POPULA','*TRANSF','*SCF IN',   !  6:10
     *            '*UPDATE','*SOLVEN','*PRINT ','*AUXILI','*HF INP',   ! 11:15
     *            '*MP2 IN','*HAMILT','*CC INP','*DFT IN','*NEVPT2',   ! 16:20
     *            '*AUXILL','*LUCITA','*STEX I','*CUBE  ','*CHOACT',   ! 21:25
     *            '*APSG I','*ENSEMB'/
C
C
C  ***************************************************
C  ***** Branch to the *-marked general sections *****
C  ***************************************************
C
      ALLOPT = .FALSE.
      PTABMO = .FALSE.
      lucita_input_provided   = .false.
      ensemble_input_provided = .false.
C
      WORD = '*GENERA'
      CALL SIR_GENINP(WORD,INPERR,IREST,NUMRUN,ALLOPT)
C
  910 PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '*') THEN
         DO 920 I = 1, NTBMOD
            IF (TABMOD(I) .EQ. WORD) THEN
               GO TO (8000,1200,1300,1400,1500,1600,1700,1800,1900,
     *                2500,2100,2200,2300,2400,2500,2600,2700,2800,
     *                2900,3000,2400,3200,3300,3400,3500,3600,3700), I
            END IF
 920     CONTINUE
         IF (WORD(1:2) .EQ. '**') GO TO 8000
         WRITE (LUPRI,'(/,3A,/)') ' Keyword ',WORD,
     *                        ' not recognized in Sirius input.'
         CALL PRTAB(NTBMOD,TABMOD,'General Sirius input keywords',LUPRI)
         CALL QUIT('Illegal keyword in top level SIRIUS input')
      ELSE IF (PROMPT .EQ. '#' .OR. PROMPT .EQ. '!') THEN
         READ (LUCMD,'(A7)') WORD
         GO TO 910
      ELSE
         IF (WORD .EQ. '.OPTION' .OR. WORD .EQ. '.ALLOPT') THEN
            CALL PRTAB(NTBMOD,TABMOD,'General Sirius input keywords',
     &                 LUPRI)
            IF (WORD .EQ. '.ALLOPT') ALLOPT = .TRUE.
            PTABMO = .TRUE.
            READ (LUCMD,'(A7)') WORD
            CALL UPCASE(WORD)
            GO TO 910
         END IF
         WRITE (LUPRI,'(/3A/)') ' Prompter "',PROMPT,'" illegal'
         CALL PRTAB(NTBMOD,TABMOD,'General Sirius input keywords',LUPRI)
         CALL QUIT('ILLEGAL PROMPT IN SIRIUS INPUT')
      END IF
C
 1200  CALL WAVINP(WORD,INPERR,ALLOPT)
       GO TO 910
 1300  CALL OPTINP(WORD,INPERR,ALLOPT)
       GO TO 910
 1400  CALL CICINP(WORD,INPERR,ALLOPT)
       GO TO 910
 1500  CALL ORBINP(WORD,INPERR,IRDMO,ALLOPT)
       GO TO 910
 1600  CALL CIVINP(WORD,INPERR,ALLOPT)
       GO TO 910
 1700  CALL STPINP(WORD,INPERR,ALLOPT)
       GO TO 910
 1800  CALL POPINP(WORD,INPERR,ALLOPT)
       GO TO 910
 1900  CALL TRAINP(WORD,INPERR,ALLOPT)  ! *TRANFORMATION
       GO TO 910
C2000  CONTINUE
C      GO TO 1200
 2100  CALL UPDINP(WORD,INPERR,ALLOPT)
       GO TO 910
 2200  CALL SOLINP(WORD,INPERR,ALLOPT)
       GO TO 910
 2300  CALL PRIINP(WORD,INPERR,ALLOPT)
       GO TO 910
C *AUXILIARY (and misspelled *AUXILLIARY for backwards compatibility)
 2400  CALL AUXINP(WORD,INPERR,ALLOPT)
       GO TO 910
 2500  CALL SCFINP(WORD,INPERR,ALLOPT)
       GO TO 910
 2600  CALL MP2INP(WORD,INPERR,ALLOPT)
       GO TO 910
 2700  CALL HAMINP(WORD,INPERR,ALLOPT)
       GO TO 910
 2800  MSYM = NSYM
       CALL CC_INPUT(WORD,IREST,MSYM)
       IF (NSYM .NE. MSYM) THEN
          CALL QUIT(
     &       'CC_INPUT ERROR: NSYM from SIRIFC .ne. NSYM in this run')
       END IF
       GOTO 910
 2900  CALL DFTINP(WORD,INPERR,ALLOPT)
       GO TO 910
 3000  CALL NEVPT2INP(WORD,INPERR,ALLOPT)
       GO TO 910
!      *LUCITA
!      this is to keep backwards compatibility - lucita input is now
!      processed via the new input reader module.
 3200  lucita_input_provided = .true.    ! *LUCITA
       call move_to_next_star(word,lucmd)
       GO TO 910
 3300  CALL STXINP(WORD,INPERR,ALLOPT)   ! *STEX
       GO TO 910
 3400  CONTINUE                          ! *CUBE
#if defined(BUILD_GEN1INT)
C      added by Bin Gao, to generate cube files to visualize the results (March 9, 2012)
       call gen1int_host_cube_init(LUCMD, LUPRI, WORD)
#else
       call quit('*CUBE not valid because GEN1INT module missing.')
#endif
       goto 910
 3500  CALL CC_ACTINP(WORD,NSYM)         ! *CHOACT
       GO TO 910
!      *APSG INPUT
!3600  CALL APSG_INPUT(WORD,INPERR,ALLOPT)
 3600  CALL QUIT('*APSG INPUT not available in this version.')
       GO TO 910
 3700  ensemble_input_provided = .true.  ! *ENSEMB
!      this is to keep backwards compatibility - ensemble dft input is now
!      processed via the new input reader module.
       call move_to_next_star(word,lucmd)
       GO TO 910
C  Here if "*END OF"
 8000  CONTINUE
       IF (ALLOPT .AND. .NOT.PTABMO) THEN
          CALL PRTAB(NTBMOD,TABMOD,'General Sirius input keywords',
     &               LUPRI)
       END IF

!      read LUCITA input if applicable
       if(lucita_input_provided)then
         call read_menu_input('DALTON.INP', lucmd,
     &                        '*LUCITA', lucita_input_provided)
         if (.not. lucita_input_provided) then
         call quit('Oops, read_menu_input could not find *LUCITA input')
         end if
       end if

!      read ensemble dft input if applicable
       if(ensemble_input_provided)then
         call read_menu_input('DALTON.INP', lucmd,
     &                        '*ENSEMB', ensemble_input_provided)
         if (.not. ensemble_input_provided) then
         call quit('Oops, read_menu_input could not find *ENSEMB input')
         end if
       end if

C
C ********  END OF READ OF INPUT DATA  **********
C
C
C     We anyway need to process CC default settings
C
      JREST = 0
      CALL CC_INPUT('**WAVE F',JREST,NSYM)
C
      CALL QEXIT('SIR_NEWINP')
      RETURN
C
C     END OF SIR_NEWINP.
C
      END
C  /* Deck sir_geninp */
      SUBROUTINE SIR_GENINP(WORD,INPERR,IREST,NUMRUN,ALLOPT)
C
C 1-Oct 1987 Hans Agren
C
#include "implicit.h"
      PARAMETER (NTABLE = 40)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7,
     &          REWORD*20, TITLIN*80, SR_FUNCS*80
      LOGICAL ALLOPT
      PARAMETER (NPDIM = 50)
      INTEGER NPATH(NPDIM)
      character(80) :: functional_line
      logical functional_line_read
C
C Used from common:
C  INFINP : FLAG(*), IORTO, DO*
C  INFORB : NSYM
C  INFPRI : P6FLAG,
C  INFTAP : LUSUPM
C  GNRINF : DOCCSD,
C  CBIHER : SUPMAT,
C
#include "maxorb.h"
#include "priunit.h"
#include "infinp.h"
#include "inforb.h"
#include "infpri.h"
#include "inftap.h"
#include "gnrinf.h"
#include "inftra.h"
#include "dftcom.h"
#include "stex.h"
C
#if !defined (VAR_NONAMELIST)
      NAMELIST /NMLSIR/ NPATH, MAXMIC
#endif
      DATA TABLE /'.CI    ','.FLAGS ','.MCSCF ','.CC    ','.RESTAR',  !  1: 5
     *            '.HF    ','.STOP  ','.PRINT ','.TITLE ','.MP2   ',  !  6:10
     *            '.NSYM  ','.CC ONL','.WESTA ','.INTERF','.BASIS ',  ! 11:15
     *            '.SCF   ','.DFT   ','.NEVPT2','.IFCSTA','.IT1STA',  ! 16:20
     *            '.KSPT2 ','.GASCI ','xXXXXXX','xXXXXXX','xXXXXXX',  ! 21:25
     *            '.CISRDF','.HFSRDF','.MCSRDF','.SRLDA ','.SRFUN ',  ! 26:30
     *            '.SRHYBR','.COPFAC','.STEX  ','.APSG  ','.VIRTRU',  ! 31:35
     *            '.MULOFA','.DSFAC ','.PVALUE','.HFXFAC','.ERF   '/  ! 36:40
C
      functional_line_read = .false.
C
      SR_FUNCS = ' ERROR: .SRFUN was not specified in input.'
C
C     ***** PROCESS INPUT FOR SIR_GENINP *****
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
         CALL UPCASE(WORD)
 1101    PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102,103,104,105,106,107,108,109,110,
     &                   111,112,113,114,115,106,117,118,119,120,
     &                   121,122,123,124,125,126,127,128,129,130,
     &                   131,132,133,134,135,136,137,138,139,140), II
               END IF
 1102       CONTINUE
            IF (WORD .EQ. '.OPTION' .OR. WORD .EQ. '.ALLOPT') THEN
               IF (.NOT. ALLOPT)
     &            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
               IF (WORD .EQ. '.ALLOPT') ALLOPT = .TRUE.
               GO TO 100
            END IF
            WRITE (LUPRI,'(/4A/)') ' Keyword ',WORD,
     *         ' not recognized for ',WORD1
         ELSE IF (PROMPT .EQ. '#' .OR. PROMPT .EQ. '!') THEN
            GO TO 100
         ELSE IF (PROMPT .EQ. '*') THEN
            GO TO 8000
         ELSE
            WRITE (LUPRI,'(/3A/2A/)')
     *         ' Keyword ',WORD,' does not begin with',
     *         ' one of the four characters ".*!#" for ',WORD1
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT(' Non-valid keyword for *WAVE FUNCTION')
C
C *** Option 1 >CI    <  do a CI calculation (before MCSCF)
 101  CONTINUE
         DOCI = .TRUE.
         FLAG(4) = .TRUE.
         NUMRUN = NUMRUN + 1
      GO TO 100
C
C *** Option 2 >FLAGS <  Set Flags directly as in old input
C     set up control flags and print flags
C     read namelist &NMLSIR ... &END
 102  CONTINUE
         CALL IZERO(NPATH,NPDIM)
#if defined (VAR_NONAMELIST)
         CALL QUIT(
     &   'ERROR, .FLAGS namelist not implemented for this machine')
#else
         READ (LUCMD,NML=NMLSIR)
#endif
C
         DO 286 I = 1,NPDIM
            IF (IABS(NPATH(I)) .GT. NFLAG) THEN
               INPERR = INPERR + 1
               WRITE (LUPRI,'(/A,I10)')
     &         ' *** ERROR Sirius .FLAGS, illegal NPATH value:',NPATH(I)
            ELSE IF (NPATH(I).LT.0) THEN
               FLAG(-NPATH(I)) = .FALSE.
               WRITE(LUPRI,'(A,I3,A)') 'INFO: flag',-NPATH(I),
     &         ' deactivated with .FLAGS'
            ELSE IF (NPATH(I).GT.0) THEN
               FLAG(NPATH(I)) = .TRUE.
               WRITE(LUPRI,'(A,I3,A)') 'INFO: flag',NPATH(I),
     &         ' activated with .FLAGS'
            ELSE
               GO TO 299
            END IF
  286    CONTINUE
  299    CONTINUE
      GO TO 100
C
C *** Option 3 >MCSCF <  MCSCF calculation
  103 CONTINUE
         DOMC    = .TRUE.
         FLAG(2) = .TRUE.
         NUMRUN  = NUMRUN + 1
      GO TO 100
C
C *** Option 4 >CC< We want to do a subsequent Coupled Cluster calculation
  104 CONTINUE
         DOCCSD = .TRUE.
         NUMRUN = NUMRUN + 1
         IF (.NOT. DOSCF) THEN
            DOSCF  = .TRUE.
            NUMRUN = NUMRUN + 1
         END IF
         FLAG(21) = .TRUE.
      GO TO 100
C
C *** Option 5 >RESTAR< restart SIRIUS, SIRIUS.RST (LUIT1) must be available
  105 CONTINUE
        IREST = 1
      GO TO 100
C
C *** Option 6 >HF    <  Restricted Hartree-Fock SCF (alias: .SCF)
  106 CONTINUE
         IF (.NOT.DOSCF) THEN
C           ... DOSCF may also be set with .CC
            DOSCF    = .TRUE.
            FLAG(21) = .TRUE.
C           ... flag for QCSCF -- quadratically convergent HF or DFT
            NUMRUN = NUMRUN + 1
         END IF
      GO TO 100
C
C *** Option 7 >STOP  <  Stop acc. to instruction given on next line
  107 CONTINUE
         READ (LUCMD,'(A20)') REWORD
         CALL UPCASE(REWORD)
         IF (INDEX(REWORD,'AFTER ') .LE. 0) THEN
            WRITE(LUPRI,'(/4A)') '"',REWORD,'" is incorrect input,',
     *         ' " AFTER <stop point>" expected after ".STOP" '
            INPERR = INPERR + 1
         ELSE
            IF (INDEX(REWORD,'INPUT') .GT. 0) THEN
               FLAG(1) = .TRUE.
            ELSE IF (INDEX(REWORD,'MO-OR') .GT. 0) THEN
               IORTO = 1
            ELSE IF (INDEX(REWORD,'GRADI') .GT. 0) THEN
               FLAG(20) = .TRUE.
            ELSE
               WRITE(LUPRI,'(/3A/A)') ' Stop "',REWORD,
     *            '" is not recognized,',
     *            ' allowed is stop after "INPUT", "MO-OR", or "GRADI"'
               INPERR = INPERR + 1
            END IF
         END IF
      GO TO 100
C
C *** Option 8 >PRINT <
  108 CONTINUE
         READ (LUCMD,*) IPRSIR
      GO TO 100
C
C *** Option  9 >TITLE <  Title
  109 CONTINUE
    1 NTIT = NTIT + 1
      IF (NTIT .LE. 6) THEN
         READ (LUCMD,'(A)') TITLE(NTIT)
         WORD = TITLE(NTIT)(1:7)
      ELSE
         READ (LUCMD,'(A7)') WORD
      END IF
      CALL UPCASE(WORD)
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '.' .OR. PROMPT .EQ. '*') THEN
         NTIT = MIN(6,NTIT - 1)
         GO TO 1101
      ELSE IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         NTIT = NTIT - 1
         GO TO 1
      ELSE
         GO TO 1
      END IF
C
C *** Option  10 >MP2  <  Moeller-Plesset PT2 calculation
  110 CONTINUE
        DOMP2    = .TRUE.
        NUMRUN   = NUMRUN + 1
        FLAG(3)  = .FALSE.
      GO TO 100
C
C *** Option 11 >NSYM  < Number of symmetries.
C    Specified by integral program, entry used for consistency check.
C    (for example in template inputs which are made for
C     a certain point group symmetry)
  111 CONTINUE
         READ (LUCMD,*) MSYM
         IF (MSYM .NE. NSYM) THEN
            INPERR = INPERR + 1
            WRITE (LUPRI,'(/A,2I5/A)')
     *   ' --- SIR_GENINP-ERROR,'//
     *   ' NSYM specified .ne. number of symmetries:',MSYM,NSYM,
     &   '     NSYM from input used in reading remaining input.'
            NSYM = MSYM
         END IF
      GO TO 100
C
C *** Option 12 >CC ONL<  CC without RHF (restart from RHF)
  112 CONTINUE
         DOCCSD   = .TRUE.
         IREST    = 1
         NUMRUN   = NUMRUN + 1
      GO TO 100
C
C *** Option 13 >WESTA <  Write information to WESTA on unit LUSYMB
  113 CONTINUE
        FLAG(26) = .TRUE.
      GO TO 100
C
C *** Option 14 >INTERF<  Write interface information to LUSIFC
  114 CONTINUE
        FLAG(25) = .TRUE.
      GO TO 100
C
C *** Option 15 >BASIS <  Specify basis set
C    Default: specified by integral program
  115 CONTINUE
        READ (LUCMD,*) (NBAS(I), I = 1,NSYM)
      GO TO 100


C *** Option 17 >DFT   <  Restricted Kohn-Sham DFT
  117 CONTINUE
         IF (DOSCF) THEN
C           ... DOSCF may also be set with .CC
            WRITE (LUPRI,'(//A)')
     &      'Input ERROR: both HF and DFT specified!'
            INPERR = INPERR + 1
         END IF
         DODFT    = .TRUE.
         DOSCF    = .TRUE.
         FLAG(21) = .TRUE.
C        ... flag for QCSCF -- quadratically convergent HF or DFT

c     dodft should be enough but it is sometimes necessary to have also dftrun set.
c     (especially for ABACUS)
         DFTRUN = .TRUE.
         DFTADD = .TRUE.
         LUSUPM   = -1
C        ... no P-supermatrix for DFT
         NUMRUN = NUMRUN + 1
 1171    READ (LUCMD, '(A)') TITLIN
         PROMPT = TITLIN(1:1)
         IF (PROMPT .EQ. '#' .OR. PROMPT .EQ. '!') GO TO 1171
         IF (PROMPT .EQ. '.' .OR. PROMPT .EQ. '*') THEN
            INPERR = INPERR + 1
            WRITE (LUPRI,'(/A/A//A)')
     &      '--> Input ERROR for input line following .DFT',
     &      '    expected functional specification but read:',
     &      TITLIN
            WORD = TITLIN(1:7)
            GO TO 1101
         END IF

         ! we save the functional line
         ! we will parse it later
         ! at this moment we do not know the integrator yet
         functional_line(1:80) = titlin(1:80)
         functional_line_read = .true.

      GO TO 100
C
C *** Option 18: NEVPT2
  118 CONTINUE
         DONEVPT  = .TRUE.
         NUMRUN = NUMRUN + 1
      GO TO 100
C
C *** Option 19: .IFCSTART
  119 CONTINUE
         IREST = 3
      GO TO 100
C
C *** Option 20: .IT1START
  120 CONTINUE
         IREST = 2
      GO TO 100
C *** Option 21: .KSPT2
  121 CONTINUE
        DOMP2    = .TRUE.
        NUMRUN   = NUMRUN + 1
        FLAG(3)  = .FALSE. ! do NOT call sircan for HF canonical orbitals
        READ (LUCMD,*) WDFTMP
      GO TO 100
C *** Option 22: .GASCI
  122 CONTINUE
        DOLUCITA = .TRUE.
        NUMRUN   = NUMRUN + 1
      GO TO 100

! *** Option 23: not in use
  123 CONTINUE
      GO TO 100

! *** Option 24: not in use
  124 CONTINUE
      GO TO 100

C *** Option 25 not used
  125 CONTINUE
      GO TO 100
C
C     ... .CISRDFT hybrid
  126 CONTINUE
         SRDFTRUN  = .TRUE.
         DOCI      = .TRUE.
         FLAG(4)   = .TRUE.
         DOCISRDFT = .TRUE.
         ADDSRI    = .TRUE.
         THRCIDFT  = 1.D-6
         LUSUPM = -1
C        ... Disable use of super-matrix integrals,
C        ... we want to use SIRFCK for generating U_sr integrals.
         NUMRUN = NUMRUN + 1
      GO TO 100
C
C     ... .HFSRDFT  Hartree-Fock srDFT hybrid
  127 CONTINUE
         SRDFTRUN  = .TRUE.
         DOSCF     = .TRUE.
         DOHFSRDFT = .TRUE.
         ADDSRI    = .TRUE.
         FLAG(21)  = .TRUE.
C        ... flag for QCHF -- quadratically convergent Hartree-Fock
         LUSUPM = -1
C        ... Disable use of super-matrix integrals,
C        ... we want to use SIRFCK for generating U_sr integrals.
         NUMRUN = NUMRUN + 1
      GO TO 100
C
C     ... .MCSRDFT hybrid
  128 CONTINUE
         SRDFTRUN  = .TRUE.
         DOMC      = .TRUE.
         DOMCSRDFT = .TRUE.
         FLAG(2)   = .TRUE.
         ADDSRI    = .TRUE.
         LUSUPM = -1
C        ... Disable use of super-matrix integrals,
C        ... we want to use SIRFCK for generating U_sr integrals.
         NUMRUN = NUMRUN + 1
      GO TO 100
C&&& .SRLDA (SR-DFT) : Use SR-LDA functionals for both
C                     exchange and correlation
  129 SR_FUNCS = 'SRXLDA SRCLDA'
      DFTADD = .TRUE.
      IF (.NOT.HFXSET) THEN
         HFXFAC = 0.0D0
         HFXSET = .TRUE.
      END IF
      GO TO 100
C&&& SRFUN : Read Correlation (SRCFUN) and
C&&&         Exchange (SRXFUN) functionals
  130 READ (LUCMD,'(A)') SR_FUNCS
      DFTADD = .TRUE.
      GO TO 100
C&&&  SRHYBR : Perform a DFT MC-SRDFT hybrid calculation.
C&&&           Treat core with regular DFT and valence with MC-SRDFT
  131 SRHYBR = .TRUE.
Chjaaj TODO what about HFXFAC here ????
      GO TO 100
C
C     .COPFAC (for some srDFT testing)
  132 READ (LUCMD,*) COPFAC
      GO TO 100
C *** Option 32: .STEX
  133 CONTINUE
         DOSTEX = .TRUE.
      GO TO 100
C *** Option 34: .APSG
  134 CONTINUE
        DOAPSG = .TRUE.
        MCTYPE = 11
        NUMRUN   = NUMRUN + 1
      GO TO 100
C *** Option 35: .VIRTRUNC
  135 CONTINUE
        READ (LUCMD,*) N_in_RN, THR_VIRTRUNC
        DO_VIRTRUNC = .TRUE.
C
  136 READ (LUCMD,*) XMULFAC_READIN   ! .MULOFA
      SRCMULOFAC = .TRUE. 
      GO TO 100 
C     Multiplication factor for DS in the localspin model 
  137 READ (LUCMD,*) DSFAC            ! .DSFAC
      DSLOCALFAC = .TRUE.
      GO TO 100
  138 READ (LUCMD,*) HEAVISIDE_PVALUE ! .PVALUE
      GO TO 100
  139 CONTINUE  ! .HFXFAC
      IF (HFXSET) WRITE(LUPRI,*)
     &   'Default HFXFAC = ',HFXFAC,' is reset by .HFXFAC input.'
      READ (LUCMD,*) HFXFAC
      IF (HFXSET) WRITE(LUPRI,*)
     &   'New     HFXFAC = ',HFXFAC
      HFXSET = .TRUE.
      GO TO 100
  140 CONTINUE ! .ERF
      READ (LUCMD,*) CHIVAL
      GO TO 100
C
 8000 CONTINUE

!     if long-range WFT short-range DFT, process functional keywords :
#ifdef MOD_SRDFT
      IF (ADDSRI) THEN
         CALL SRFUN_INPUT(SR_FUNCS)
         IF (CHIVAL .LT. 0.0D0) THEN
            CHIVAL = 0.4D0
            WRITE(LUPRI,'(/A)')
     &      ' mu value for srDFT set to 0.4 (default)'
         ELSE
            WRITE(LUPRI,'(/A,F5.2,A)')
     &      ' mu value for srDFT is',CHIVAL,' (user input)'
         END IF
         DOSRIN = .TRUE.
         ! default is erf
      END IF
#else
      IF (ADDSRI) CALL QUIT('srDFT is not available in this version')
#endif

      ! here we parse the functional
      if (functional_line_read) then
         call dftsetfunc(functional_line, inperr)
         IF (WDFTMP.NE.0.0D0 .AND. .NOT.DOMP2) THEN ! double-hybrid DFT requested
            DOMP2    = .TRUE.
            NUMRUN   = NUMRUN + 1
            FLAG(3)  = .FALSE. ! do NOT call sircan for HF canonical orbitals
         END IF
      end if

C
      IF (DOAPSG) THEN
         IF (DOMC) THEN
            WRITE(LUPRI,'(/2A)') WORD1,
     &      ' input error: .MCSCF and .APSG not both allowed.'
            INPERR = INPERR + 1
         END IF
         DOMC = .TRUE.
         call quit('APSG is not available in this version')
      END IF
C
      RETURN
      END

C  /* Deck haminp */
      SUBROUTINE HAMINP(WORD,INPERR,ALLOPT)
C
C Copyright Nov 1990 Hans Joergen Aa. Jensen
C
#include "implicit.h"
      PARAMETER (NTABLE = 2)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7
      LOGICAL ALLOPT
C
C Used from common:
C
C INFINP : NFIELD,EFIELD,LFIELD
C INFPRI : P6FLAG
C
#include "maxorb.h"
#include "priunit.h"
#include "infinp.h"
#include "infpri.h"
C
      DATA TABLE /'.FIELD ','.PRINT '/
C
C     ***** PROCESS INPUT FOR HAMINP *****
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
         CALL UPCASE(WORD)
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102), II
               END IF
 1102       CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               IF (.NOT.ALLOPT)
     &            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/4A/)') ' Keyword ',WORD,
     *         ' not recognized for ',WORD1
         ELSE IF (PROMPT .EQ. '#' .OR. PROMPT .EQ. '!') THEN
            GO TO 100
         ELSE IF (PROMPT .EQ. '*') THEN
            RETURN
         ELSE
            WRITE (LUPRI,'(/3A/2A/)')
     *         ' Keyword ',WORD,' does not begin with',
     *         ' one of the four characters ".*!#" for ',WORD1
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal keyword for "*HAMILTONIAN"')
C
C *** Option 1 >FIELD TERM<   Field dependent Hamiltonian term
  101 CONTINUE
         NFIELD = NFIELD + 1
         IF (NFIELD .LE. MXFELT) THEN
            READ(LUCMD,ERR=103,FMT=*)     EFIELD(NFIELD)
            READ(LUCMD,ERR=104,FMT='(A)') LFIELD(NFIELD)
         ELSE
            INPERR = INPERR + 1
            WRITE (LUPRI,'(/A/A,I5/)')
     *         ' Input ERROR: number of ".FIELD TERM" entries',
     *         '    is greater than current maximum of',MXFELT
            READ (LUCMD,'()')
            READ (LUCMD,'()')
         END IF
      GO TO 100
C
C *** Option 2 >PRINT <  General print level in SIRH1
  102 READ(LUCMD,*) IPRH1
      IF (IPRH1 .GT. 0) P6FLAG(18) = .TRUE.
      GO TO 100
C
 103  WRITE(LUPRI,*) '.FIELD: could not read field strength.'
      CALL QUIT('.FIELD: could not read field strength.')
 104  WRITE(LUPRI,*) '.FIELD: could not read field operator.'
      CALL QUIT('.FIELD: could not read field operator.')
      RETURN
      END
C  /* Deck wavinp */
      SUBROUTINE WAVINP(WORD,INPERR,ALLOPT)
C
C Purpose: process *CONFIGURATION input
C
C 1-Oct 1987 Hans Agren
C
#include "implicit.h"
      PARAMETER (NTABLE = 12)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7
      LOGICAL ALLOPT
C
C Used from common:
C
C INFORB : NASH, NISH, NSYM,...
C INFINP : NACTEL,ISPIN,FLAG(*),MCTYP,..
C SPINFO : MULTS, MS2
C
#include "maxorb.h"
#include "priunit.h"
#include "inforb.h"
#include "infinp.h"
#include "spinfo.h"
C
      LOGICAL RAS1EL
C
      DATA TABLE /'.CAS SP','.INACTI','.ELECTR','.MS2   ','.SPIN M',
     *            '.SYMMET','.RAS1 H','.RAS1 S','.RAS2 S','.RAS3 S',
     *            '.RAS1 E','.RAS3 E'/
C
      NHL1MN = -1
      NHL1MX = -1
      RAS1EL = .FALSE.
C
C     ***** PROCESS INPUT FOR WAVINP *****
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
         CALL UPCASE(WORD)
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102,103,104,105,106,107,108,109,110,
     *                   111,112), II
               END IF
 1102       CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               IF (.NOT.ALLOPT)
     &            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/4A/)') ' Keyword ',WORD,
     *         ' not recognized for ',WORD1
         ELSE IF (PROMPT .EQ. '#' .OR. PROMPT .EQ. '!') THEN
            GO TO 100
         ELSE IF (PROMPT .EQ. '*') THEN
            GO TO 2000
         ELSE
            WRITE (LUPRI,'(/3A/2A/)')
     *         ' Keyword ',WORD,' does not begin with',
     *         ' one of the four characters ".*!#" for ',WORD1
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal keyword for "*CONFIG"')
C
C *** Option 1 >CAS  S<  CAS calculation: Active orbitals
  101 CONTINUE
      IF (MCTYPE .NE. 0 .AND. MCTYPE.NE.1 .AND. MCTYPE.NE.11) THEN
         WRITE (LUPRI,*) '*CONFIG ERROR, both CAS and RAS specified.'
         INPERR = INPERR + 1
      END IF
         IF (MCTYPE .NE. 11) MCTYPE = 1 ! if (not APSG) then CASSCF
         READ(LUCMD,*,ERR=8010) (NAS2(I),I=1,NSYM)
         NAS1(1:NSYM) = 0
         NAS3(1:NSYM) = 0
      GO TO 100
C
C *** Option 2 >INACTI<  Inactive orbitals
  102 CONTINUE
         READ(LUCMD,*,ERR=8010) (NISH(I),I=1,NSYM)
         NRHFT = ISUM(NSYM,NRHF,1)
         IF (NRHFT .EQ. 0) THEN
            NRHF(1:NSYM) = NISH(1:NSYM)
         END IF
      GO TO 100
C
C *** Option 3 >ELECTR< Number of active electrons
  103 READ(LUCMD,*) NACTEL
      GO TO 100
C
C *** Option 4 >MS2   < 2 * M_S value
  104 READ(LUCMD,*) MS2
      IF (MULTS .GT. 0) THEN
         IF (MS2 .LT. -(MULTS-1) .OR. MS2 .GT. MULTS-1) THEN
            WRITE (LUPRI,'(/A,2I10/)')
     &      ' *CONFIG ERROR, .SPIN M and .MS2 not compatible:',MULTS,MS2
            INPERR = INPERR + 1
         END IF
      END IF
      GO TO 100
C
C *** Option 5 >SPIN M<  state spin multiplicity
C                        SPIN QUANTUM NUMBER*2 + 1
  105 READ(LUCMD,*) MULTS
      IF (MS2 .NE. -12345) THEN
         IF (MS2 .LT. -(MULTS-1) .OR. MS2 .GT. MULTS-1) THEN
            WRITE (LUPRI,'(/A,2I10/)')
     &      ' *CONFIG ERROR, .SPIN M and .MS2 not compatible:',MULTS,MS2
            INPERR = INPERR + 1
         END IF
      END IF
      GO TO 100
C
C *** Option 6 >SYMMET<  Space symmetry of MCSCF wave function
  106 READ(LUCMD,*) LSYM
      GO TO 100
C
C *** Option 7 >RAS1 H< number of holes in RAS1
  107 CONTINUE
         IF (RAS1EL) THEN
            WRITE (LUPRI,'(/A/)')
     *         ' *CONFIG ERROR, both .RAS1 H and .RAS1 E specified.'
            INPERR = INPERR + 1
         END IF
         READ(LUCMD,*) NHL1MN,NHL1MX
         RAS1EL = .TRUE.
      GO TO 100
C
C *** Option 8 >RAS1 S<   RAS1 orbital space
  108 CONTINUE
      IF (MCTYPE .NE. 0 .AND. MCTYPE .NE. 2) THEN
         WRITE (LUPRI,*) '*CONFIG ERROR, both CAS and RAS specified.'
         INPERR = INPERR + 1
      END IF
         MCTYPE = 2
         READ(LUCMD,*,ERR=8010) (NAS1(I),I=1,NSYM)
      GO TO 100
C
C *** Option 9 >RAS2 S<   RAS2 orbital space
  109 CONTINUE
      IF (MCTYPE .NE. 0 .AND. MCTYPE .NE. 2) THEN
         WRITE (LUPRI,*) '*CONFIG ERROR, both CAS and RAS specified.'
         INPERR = INPERR + 1
      END IF
         MCTYPE = 2
         READ(LUCMD,*,ERR=8010) (NAS2(I),I=1,NSYM)
         NAS2T = 0
         do i = 1, NSYM
           NAS2T = NAS2T + NAS2(i)
         end do
      GO TO 100
C
C *** Option 10 >RAS3 S<   RAS3 orbital space
  110 CONTINUE
      IF (MCTYPE .NE. 0 .AND. MCTYPE .NE. 2) THEN
         WRITE (LUPRI,*) '*CONFIG ERROR, both CAS and RAS specified.'
         INPERR = INPERR + 1
      END IF
         MCTYPE = 2
         READ(LUCMD,*,ERR=8010) (NAS3(I),I=1,NSYM)
      GO TO 100
C
C *** Option 11 >RAS1 E<   RAS1 electrons
  111 CONTINUE
         IF (RAS1EL) THEN
            WRITE (LUPRI,'(/A/)')
     *         ' *CONFIG ERROR, both .RAS1 H and .RAS1 E specified.'
            INPERR = INPERR + 1
         END IF
         READ(LUCMD,*) NELMN1,NELMX1
         RAS1EL = .TRUE.
      GO TO 100
C
C *** Option 12 >RAS3 E<   RAS3 electrons
  112 CONTINUE
         READ(LUCMD,*) NELMN3, NELMX3
      GO TO 100
C
 2000 CONTINUE
      NASHT = 0
      DO I = 1,NSYM
         NASH(I) = NAS1(I) + NAS2(I) + NAS3(I)
         NASHT   = NASHT   + NASH(I)
      END DO
      IF (NHL1MX .NE. -1 .OR. NHL1MN .NE. -1) THEN
         NAS1T = 0
         DO I = 1,NSYM
            NAS1T = NAS1T + NAS1(I)
         END DO
         NELMN1 = 2*NAS1T - NHL1MX
         NELMX1 = 2*NAS1T - NHL1MN
      END IF
      IF (MCTYPE .GT. 0 .AND. NACTEL .EQ. -1) THEN
         INPERR = INPERR + 1
         WRITE(LUPRI,'(/A)')
     &   ' *CONFIG ERROR: .ELECTRONS not specified'
      END IF
      RETURN
C ----------------------------------------------
 8010 CONTINUE ! NSYM error (probably)
      WRITE (LUPRI,'(///2A/A,I5/)')
     &   ' INPUT ERROR in *OPTIMI for keyword ',WORD,
     &   ' INPUT ERROR is probably wrong number of symmetries, expected'
     &      ,NSYM
      CALL QUIT('Input error in *OPTIMI')
      END
C  /* Deck scfinp */
      SUBROUTINE SCFINP(WORD,INPERR,ALLOPT)
C
C 8-May-1989 Hans Joergen Aa. Jensen
C 13-May-1995 K.Ruud, added keyword .ELECTRONS.
C 1999 Helgaker, DFT input added
C
#include "implicit.h"
      PARAMETER (NTABLE = 26)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7
      LOGICAL ALLOPT
C
C Used from common:
C
C GNRINF : THR_REDFAC
C SCBRHF : NFRRHF(8), IOPRHF, RHFCAN, MXHFMA, MXHFMI, MAXFCK, MXDIIS, THRRHF
C INFINP : FLAG(*),DIRFCK
C INFORB : NRHF(8)
C INFTAP : LUSUPM
C INFPRI : IPRSIR,IPRRHF
C DFTACB : DFTIPT, DFTBR1, DFTBR2
C
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
C
#include "gnrinf.h"
#include "scbrhf.h"
#include "infinp.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "dftcom.h"
#include "dfterg.h"
#include "dftacb.h"
C
      LOGICAL NODIIS,SET,ELECTR
      SAVE SET,ELECTR
      DATA SET,ELECTR /2*.FALSE./
C
      DATA TABLE /'.DOUBLY','.THRESH','.FOCK I','.H1VIRT','.NOQCSC',
     *            '.PRINT ','.OPEN S','.MAX MA','.MAX MI','.NONCAN',
     *            '.FROZEN','.DIRFOC','.NODIIS','.MAX DI','.C2DIIS',
     *            '.MAX ER','.NO BAC','.COREHO','.CORERE','.HF OCC',
     *            '.ELECTR','.AUTOCC','.FC MVO','.SINGLY','.SHIFT ',
     *            '.MEP   '/
C
C     ***** PROCESS INPUT FOR SCFINP *****
C
      IF (WORD .EQ. '*HF INP') THEN
         NWARN = NWARN + 1
         WRITE (LUPRI,'(//A)')
     &   'WARNING: obsolete "*HF INP" specified, please use "*SCF IN".'
         WORD = '*SCF IN'
      END IF
      IF (SET .AND. WORD .EQ. '*SCF IN') THEN
         INPERR = INPERR + 1
         WRITE(LUPRI,*)
     &   ' Input ERROR: *SCF INPUT specified more than once'
      END IF
      NODIIS = .FALSE.
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
         CALL UPCASE(WORD)
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102,103,104,105,106,107,108,109,110,
     &                   111,112,113,114,115,116,117,118,119,120,
     &                   121,122,123,124,125,126),II
               END IF
 1102       CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               IF (.NOT.ALLOPT)
     &            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/4A/)') ' Keyword ',WORD,
     *         ' not recognized for ',WORD1
         ELSE IF (PROMPT .EQ. '#' .OR. PROMPT .EQ. '!') THEN
            GO TO 100
         ELSE IF (PROMPT .EQ. '*') THEN
            GO TO 2000
         ELSE
            WRITE (LUPRI,'(/3A/2A/)')
     *         ' Keyword ',WORD,' does not begin with',
     *         ' one of the four characters ".*!#" for ',WORD1
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal keyword for "*SCF INPUT"')
C
C *** Option 1 >DOUBLY<  doubly occ. orb.s, for RHF and MP2 calc. or DFT
  101 READ(LUCMD,*) (NRHF(I),I=1,NSYM)
      IF (INIOCC .LT. 0) AUTOCC = .FALSE.
C     ... AUTOCC may have been set true by .AUTOCC
      IF (ELECTR) THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,*)
     &' .DOUBLY must not be specified together with .ELECTRONS'
      END IF
      NRHFEL = 0
      INIOCC = 1
      GO TO 100
C
C *** Option 2 >THRESH<  Hartree-Fock convergence threshold
  102 READ (LUCMD,*) THRRHF
      GO TO 100
C
C *** Option 3 >FOCK I<  Maximum number of Fock iterations RHF
  103 READ (LUCMD,*) MAXFCK
      GO TO 100
C
C *** Option 4 >H1VIRT<  Virtual orbitals that diagonalize 1-el. H
  104 CALL IZERO(NMVO,8)
      NMVOT = 0
      IF (FCVORB) THEN
         NINFO = NINFO + 1
         WRITE (LUPRI,*)
     &'SCFINP INFO: more than one specification of .H1VIRT and .FC MVO'
         WRITE (LUPRI,*) '.H1VIRT will be used.'
      END IF
      FCVORB = .TRUE.
      GO TO 100
C
C *** Option 5 >NOQCSCF<  No quadratically convergent SCF
  105 FLAG(21) = .FALSE.
      GO TO 100
C
C *** Option 6 >PRINT <  General print level in Hartree-Fock
  106 READ(LUCMD,*) IPRRHF
      GO TO 100
C
C *** Option 7 >OPEN S< Symmetry of the open shell in a one open
C                       shell calculation.
  107 READ (LUCMD,*) IOPRHF
      IF (INIOCC .LT. 0) AUTOCC = .FALSE.
C     ... AUTOCC may have been set true by .AUTOCC
      IF (ELECTR) THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,*)
     &' .OPEN SHELL is not allowed together with .ELECTRONS'
      END IF
      IF (HSROHF) THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,*)
     &' .OPEN SHELL is not allowed together with .SINGLY'
      END IF
      IF (IOPRHF .LT. 1 .OR. IOPRHF .GT. NSYM) THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,*)
     &' Illegal symmetry of .OPEN SHELL :',IOPRHF
      END IF
      NRHFEL = 0
      INIOCC = 1
      GO TO 100
C
C *** Option 8 >MAX MA<    maximum number of macro iterations
  108 READ (LUCMD,*) MXHFMA
      GO TO 100
C
C *** Option 9 >MAX MI< maximum number of micro iterations per macro
C                       iteration
  109 READ (LUCMD,*) MXHFMI
      GO TO 100
C
C *** Option 10 >NONCAN< No transformation to canonical orbitals
  110 RHFCAN = .FALSE.
      GO TO 100
C
C *** Option 11 >FROZEN<  Frozen RHF orbitals per symmetry
  111 READ (LUCMD,*) (NFRRHF(I),I=1,NSYM)
      GO TO 100
C
C *** Option 12 >DIRFOC< Direct Fock matrix construction
  112 DIRFCK = .TRUE.
      LUSUPM = -1
      GO TO 100
C
C *** Option 13 >NODIIS< Do not use DIIS in Fock iterations
  113 NODIIS = .TRUE.
      GO TO 100
C
C *** Option 14 >MAX DI< Max. number of DIIS iterations
  114 READ (LUCMD,*) MXDIIS
      GO TO 100
C
C *** Option 15 >C2DIIS< Use Sellers C2-DIIS algorithm for DIIS
  115 FLAG(22) = .TRUE.
      GO TO 100
C
C *** Option 16 >MAX ER< Max. number of DIIS error vectors
  116 READ (LUCMD,*) MAXEVC
      GO TO 100
C
C *** Option 17 >NO BAC< No backsteps
  117 CONTINUE
      BCKSTP = .FALSE.
      GO TO 100
C
C *** Option 18 >COREHO< core hole calculation
  118 READ (LUCMD,*) JCHSYM,JCHORB
      CORHOL = .TRUE.
      GO TO 100
C
C *** Option 19 >CORERE< core hole relaxation
  119 CORRLX = .TRUE.
      GO TO 100
C *** Option 20 >HF OCC< (temp. for Dalton 2.0 -- Oct 2003/hjaaj)
  120 CONTINUE
         NWARN = NWARN + 1
         WRITE(LUPRI,'(//A)')
     &   'WARNING: ".HF OCC" is obsolete, please use ".DOUBLY OCCUPIED"'
         GO TO 101
c     GO TO 100
C *** Option 21 >ELECTR< No. of electrons, for automagic det. of occ. orb.s
  121 READ(LUCMD,*) NRHFEL
      IF (INIOCC .EQ. 1) THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,*) ' .ELECTR must not be specified'//
     &   ' together with any of .DOUBLY, .OPEN S and .SINGLY'
      END IF
      AUTOCC = .TRUE.
      ELECTR = .TRUE.
      IF (INIOCC .LT. 0) INIOCC = 0
      GO TO 100
C *** Option 22 >AUTOCC< automatic determination of occ. orbs. in DIIS,
C                        with given inital occ. from user input or file
  122 IF (NSYM .EQ. 1) THEN
         AUTOCC = .FALSE.
         WRITE (LUPRI,'(/A)')
     &      'INFO: .AUTOCC ignored because C_1 symmetry'
      ELSE
         AUTOCC = .TRUE.
         IF (INIOCC .LT. 0) INIOCC = 0
      END IF
      GO TO 100
C
C *** Option 23 >FC MVO<  Modified virtual orbitals from FC from NMVO
C     Reference: C.W. Bauschlicher, JCP 72 (1980) 880.
  123 READ (LUCMD,*) (NMVO(I), I = 1,NSYM)
      NMVOT = ISUM(8,NMVO,1)
      IF (FCVORB) THEN
         NINFO = NINFO + 1
         WRITE (LUPRI,*)
     &'SCFINP INFO: more than one specification of .H1VIRT and .FC MVO'
         WRITE (LUPRI,*) 'Last .FC MVO will be used.'
      END IF
      FCVORB = .TRUE.
      GO TO 100
C
C *** Option 24 >SINGLY< High spin orbitals input
 124  HSROHF = .TRUE.
      READ(LUCMD,*) (NROHF(I),I=1,NSYM)
C
      IF (INIOCC .LT. 0) AUTOCC = .FALSE.
C     ... AUTOCC may have been set true by .AUTOCC
      IF (ELECTR .OR. AUTOCC) THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,*)
     &' .SINGLY must not be specified together with .ELECTR or .AUTOCC'
      END IF
      IF (IOPRHF .GT. 0) THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,*)
     &' .SINGLY must not be specified together with .OPEN SHELL'
      END IF
      NRHFEL = 0
      INIOCC = 1
      NASHT  = ISUM(NSYM,NROHF,1)
      IF (NASHT .LT. 0) THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,*)
     &' Input error: negative occupation for .SINGLY :',NROHF(1:NSYM)
      ELSE IF (NASHT .EQ. 0) THEN
         HSROHF = .FALSE.
         WRITE (LUPRI,*)
     &' INFO: .SINGLY specified with zero occupations;'//
     &' calculation type changed to closed-shell SCF'
      END IF
      GO TO 100
C *** Option 25 >SHIFT < Shift level
 125  READ(LUCMD,*) DEFLVL
      DEFLVL = -ABS(DEFLVL)
      GO TO 100
C     Option 26 >MEP   <
 126  DOMEP = .TRUE.
      GO TO 100
C
 2000 CONTINUE
      IF (THR_REDFAC .GT. 0.0D0) THEN
C         ICHANG = ICHANG + 1
         WRITE (LUPRI,'(3A,1P,D10.2)') '@ INFO ',WORD1,
     &   ' thresholds multiplied with general factor',THR_REDFAC
         THRRHF = THRRHF*THR_REDFAC
      END IF
      IF (NODIIS) MXDIIS = 0
      IF (CORHOL) THEN
         IF (IOPRHF .NE. 0) THEN
            INPERR = INPERR + 1
            WRITE (LUPRI,*) 'SCFINP input ERROR:'//
     &         ' .OPEN SHELL and .COREHOLE not compatible'
         END IF
         IOPRHF = JCHSYM
         NRHF(JCHSYM) = NRHF(JCHSYM) - 1
         IF (NRHF(JCHSYM) .LT. 0) THEN
            INPERR = INPERR + 1
            WRITE (LUPRI,*) 'SCFINP input ERROR:'//
     &         ' no HF OCCUPATION in core hole symmetry'
         END IF
         MXDIIS = 0
         MAXFCK = 0
      END IF
      IF (IOPRHF .NE. 0 .AND. HSROHF) THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,*) 'SCFINP input ERROR:'//
     &      ' .OPEN SHELL and .SINGLY both specified'
      END IF
      SET = .TRUE.
      RETURN
      END
C  /* Deck optinp */
      SUBROUTINE OPTINP(WORD,INPERR,ALLOPT)
C
C 1-Oct 1987 Hans Agren
C
#include "implicit.h"
      PARAMETER (NTABLE = 36)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7,
     *          REWORD*12
      LOGICAL ALLOPT
C
C Used from common:
C
C GNRINF : THR_REDFAC
C INFINP : FLAG(*),ICHECK,*ROOT,ISTATE,MAX*,THR*,SUPSYM,THRSSY
C INFVAR : NWOPT, JWOP(*)
C INFDIM : MAXPHP
C INFOPT : RAT*
C INFOPT : NFRO
C INFPRI : LUPRI
C
#include "maxorb.h"
#include "priunit.h"
#include "gnrinf.h"
#include "infinp.h"
#include "infvar.h"
#include "infdim.h"
#include "inftap.h"
#include "infopt.h"
#include "inforb.h"
C
      DATA TABLE /'.ABSORP','.ALWAYS','.ACTROT','.FROZEN','.MAXAPM',
     *            '.SYM CH','.THRCGR','.STATE ','.COREHO','.CORERE',
     *            '.MAX MA','.MAXABS','.MAX CI','.EXACTD','.MAX MI',
     *            '.NEO AL','.PHPRES','.NO ABS','.NR ALW','.OPTIMA',
     *            '.ORB_TR','.MCAVER','.SIMULT','.THRESH','.OLSEN ',
     *            '.TRACI ','.FOCKDI','.NOTRAC','.CI PHP','.NO ACT',
     *            '.DETERM','.NATONL','.FOCKON','.CI PRO','.SUPSYM',
     &            '.THRSSY'/
C
C Internal flags
      IFNSIM = 0
      ISTATE_SPECIFIED = 0
C
C     ***** PROCESS INPUT FOR OPTINP *****
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
         CALL UPCASE(WORD)
 1001    PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102,103,104,105,106,107,108,109,110,
     *                   111,112,113,114,115,116,117,118,119,120,
     *                   121,122,123,124,125,126,127,128,129,130,
     *                   131,132,133,134,135,136), II
               END IF
 1102       CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               IF (.NOT.ALLOPT)
     &            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/4A/)') ' Keyword ',WORD,
     *         ' not recognized for ',WORD1
         ELSE IF (PROMPT .EQ. '#' .OR. PROMPT .EQ. '!') THEN
            GO TO 100
         ELSE IF (PROMPT .EQ. '*') THEN
            GO TO 9000
         ELSE
            WRITE (LUPRI,'(/3A/2A/)')
     *         ' Keyword ',WORD,' does not begin with',
     *         ' one of the four characters ".*!#" for ',WORD1
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal keyword for "*OPTIMIZATION"')
C
C
C *** Option 1 >ABSORP<  Orbital absorption
 101  CONTINUE
         READ (LUCMD,'(A12)') REWORD
         CALL UPCASE(REWORD)
         IF (INDEX(REWORD,'LEVEL 1') .GT. 0) THEN
            FLAG(51) = .TRUE.
         ELSE IF (INDEX(REWORD,'LEVEL 2') .GT. 0) THEN
            FLAG(52) = .TRUE.
         ELSE IF (INDEX(REWORD,'LEVEL 3') .GT. 0) THEN
            FLAG(53) = .TRUE.
         ELSE
            WRITE(LUPRI,'(4A)') '"',REWORD,'" IS INCORRECT ',
     *         'INPUT, " LEVEL n" (n=1,2,3) EXPECTED AFTER ".ABSORP"'
            INPERR = INPERR + 1
            IF (REWORD(1:1) .EQ. '*' .OR. REWORD(1:1) .EQ. '.') THEN
               WORD = REWORD(1:7)
               GO TO 1001
            END IF
         END IF
      GO TO 100
C
C *** Option 2 >ALWAYS<  Absorption all macroiterations
 102  CONTINUE
        FLAG(54) = .TRUE.
      GO TO 100
C
C *** Option 3 >ACTROT<  include specified active-active rotations
  103 READ (LUCMD,*) NWOPT
      IF (NWOPT .GT. MAXWOP) THEN
          WRITE(LUPRI,'(A,I8)')
     *    ' TOO MANY ACTIVE ROTATIONS SPECIFIED, MAX =',MAXWOP
          CALL QUIT('SIRIUS INPUT: TOO MANY ACTIVE ROTATIONS SPECIFIED')
      END IF
        DO 195 I = 1,NWOPT
          READ (LUCMD,*) JWOP(1,I),JWOP(2,I)
  195   CONTINUE
        FLAG(23)=.TRUE.
      GO TO 100
C
C *** Option 4 >FROZEN< frozen core orbitals
  104 CONTINUE
      READ (LUCMD,*) (NFRO(I),I=1,NSYM)
      GO TO 100
C
C *** Option 5 >MAXAPM< maximum number absorptions within a macro iteration
  105 CONTINUE
      READ (LUCMD,*) MAXAPM
      GO TO 100
C
C *** Option 6 >SYM CH<  Check symmetry of CI vectors and remove
C                        vectors of undesired symmetry.
  106 CONTINUE
         READ (LUCMD,*) ICHECK
      GO TO 100
C
C *** Option 7 >THRCGR<  Threshold for print of CI gradient.
  107 CONTINUE
        READ (LUCMD,*) THRCGR
      GO TO 100
C
C *** Option 8 >STATE <  Root to converge calculation
  108 CONTINUE
         IF ( ISTATE_SPECIFIED .NE. 0) THEN
            WRITE (LUPRI,*) '*OPTIMI ERROR, not both .STATE and .MCAVER'
            CALL QUIT('OPTINP ERROR: not both .STATE and .MCAVER')
         END IF
         ISTATE_SPECIFIED = 1
C
         READ (LUCMD,*) ISTATE
         NMCAVER = 1
         IROOT(1) = ISTATE
C
         IF ( ISTATE .GT. MAXRTS ) THEN
            WRITE (LUPRI,*) '*OPTIMI ERROR: .STATE input is greater'//
     &                      ' than internal maximum of',MAXRTS
            CALL QUIT(
     &         'OPTINP ERROR: .STATE input exceeds internal limit')
         END IF
         IF (ISTATE .GT. 1) THEN
            RATREJ = 0.40D0
            RATMIN = 0.60D0
            RATGOD = 0.80D0
         END IF
      GO TO 100
C
C *** Option 9 >COREHO< core hole calculation
  109 CONTINUE
         READ (LUCMD,*) JCHSYM,JCHORB
         CORHOL = .TRUE.
         FLAG(52) = .TRUE.
      GO TO 100
C
C *** Option 10 >CORERE< core hole relaxation
  110 CONTINUE
         CORHOL = .TRUE.
         CORRLX = .TRUE.
      GO TO 100
C
C *** Option 11 >MAX MA<  max # of macro iterations
  111 READ(LUCMD,*) MAXMAC
      GO TO 100
C
C *** Option 12 >MAXABS<  max # of macro iterations with absorption
  112 READ(LUCMD,*) MAXABS
      GO TO 100
C
C *** Option 13 >MAX CI<  max # of CI iterations before MCSCF
  113 READ(LUCMD,*) MAXCIT
      GO TO 100
C
C *** Option 14 >EXACTD<  Use the exact orbital Hessian diagonal.
  114 CONTINUE
      FLAG(12) = .TRUE.
      GO TO 100
C
C *** Option 15 >MAX MI<  max # of micro iterations per macro it.
  115 READ(LUCMD,*) MAXJT
      GO TO 100
C
C *** Option 16 >NEO AL<  always norm-extended optimization
C                         (.NR ALWAYS takes precedence)
  116 CONTINUE
        FLAG(38) = .TRUE.
      GO TO 100
C
C *** Option 17 >.PHPRES< Select configurations for PHP matrix based on
C                         largest residual rather than lowest diagonal
C                         elements.
  117 RESPHP = .TRUE.
      GO TO 100
C
C *** Option 18 >NO ABS<  No orbital absorption
C (this is default except for .COREHOLE, however, flags disabled if
C  this option is specified after the .ABSORPTION in input)
  118 CONTINUE
         FLAG(51) = .FALSE.
         FLAG(52) = .FALSE.
         FLAG(53) = .FALSE.
         MAXABS   = 0
         MAXAPM   = 0
      GO TO 100
C
C *** Option 19 >NR ALW<  Always Newton-Raphson optimization
  119 CONTINUE
        FLAG(39) = .TRUE.
      GO TO 100
C
C *** Option 20 >OPTIMA<  Optimal orbital trial vectors
  120 CONTINUE
        FLAG(41) = .TRUE.
      GO TO 100
C
C *** Option 21 >ORB_TR<  Use orbital trial vectors as start vectors
  121 CONTINUE
        FLAG(80) = .TRUE.
      GO TO 100
C
C *** Option 22 >MCAVER< State Average MCSCF
C               Convention: ISTATE = max root in average
C               /hjaaj nov 2006
  122 CONTINUE
         IF ( ISTATE_SPECIFIED .NE. 0) THEN
            WRITE (LUPRI,*) '*OPTIMI ERROR, not both .STATE and .MCAVER'
            CALL QUIT('OPTINP ERROR: not both .STATE and .MCAVER')
         END IF
         ISTATE_SPECIFIED = 1
         READ (LUCMD,*) NMCAVER
         READ (LUCMD,*) (IROOT(I),I=1,NMCAVER)
         READ (LUCMD,*) (WEIGHT_MCAVER(I),I=1,NMCAVER)
C
C        Consistency checks, find ISTATE,
C        and ensure MCAVER roots are ordered
         ISTATE = 0
         DO I = 1,NMCAVER
            IF (WEIGHT_MCAVER(I) .LE. 0.0D0) THEN
               CALL QUIT('.MCAVER ERROR: state weight not positive!!!')
            END IF
            ISTATE = MAX(ISTATE,IROOT(I))
            DO J = 1,I-1
               IF (IROOT(I) .EQ. IROOT(J)) THEN
                  CALL QUIT('.MCAVER ERROR: same state duplicated')
               ELSE IF (IROOT(I) .LT. IROOT(J)) THEN
                  K = IROOT(I)
                  IROOT(I) = IROOT(J)
                  IROOT(J) = K
                  WK = WEIGHT_MCAVER(I)
                  WEIGHT_MCAVER(I) = WEIGHT_MCAVER(J)
                  WEIGHT_MCAVER(J) = WK
               END IF
            END DO
         END DO
C
         IF ( ISTATE .GT. MAXRTS ) THEN
            WRITE (LUPRI,*) 'OPTINP ERROR: .STATE input is greater than'
            WRITE (LUPRI,*) 'internal maximum of',MAXRTS
            CALL QUIT(
     &         'OPTINP ERROR: .STATE input exceeds internal limit')
         END IF
         IF (ISTATE .GT. 1) THEN
            RATREJ = 0.40D0
            RATMIN = 0.60D0
            RATGOD = 0.80D0
         END IF
      GO TO 100
C
C *** Option 23 >SIMULT<  # of simultaneous roots
  123 IFNSIM = 1
         READ (LUCMD,*) NROOTS, LROOTS
         IF (LROOTS.LT.NROOTS) LROOTS = NROOTS
      GO TO 100
C
C *** Option 24 >THRESH< Threshold for energy gradient in optimization
  124 READ (LUCMD,*) THRMC
      GO TO 100
C
C *** Option 25 >OLSEN< Use Jeppe Olsen's generalization of the
C                       Davidson algorithm.
  125 JOLSEN = .TRUE.
      GO TO 100
C
C *** Option 26 >TRACI < Call TRACI for Transform of CI-vectors
  126 IMCCNO = 1
      GO TO 100
C
C *** Option 27 >FOCKDI< Use only Fock contribution to orbital hessian
C                        diagonal
  127 FLAG(12) = .FALSE.
      GO TO 100
C
C *** Option 28 >NOTRAC<
C     Disable default transformation of CI-vectors by TRACI
  128 IMCCNO = -1
      GO TO 100
C
C *** Option 29 >CI PHP<  Explicit CI Hamiltonian matrix set up
  129 CONTINUE
        READ (LUCMD,*) MAXPHP
      GO TO 100
C
C *** Option 30 >NO ACT<  No active-active rotations in RAS
  130 CONTINUE
        FLAG(24) = .TRUE.
      GO TO 100
C
C *** Option 31 >DETERM<  use determinant basis instead of CSF basis
  131 CONTINUE
        FLAG(27) = .TRUE.
        FLAG(66) = .TRUE.
C       ... 66: do not use permutation symmetry because
C               this causes numerical problems
      GO TO 100
C
C *** Option 32 >NATONL< Only natural (no co:s) for TRACI and output
 132  CONTINUE
        IMCCNO   = 1
        FLAG(46) = .TRUE.
      GO TO 100
C
C *** Option 33 >FOCKON< Only FOCK (no no:s ) for TRACI and output
 133  CONTINUE
        IMCCNO   = 1
        FLAG(48) = .TRUE.
      GO TO 100
C
C *** Option 34 >CI PRO< set CI program to be used as MCSCF backbone
 134  CONTINUE
        ci_program = '         '
        READ (LUCMD,*) ci_program(1:9)
      GO TO 100
C
C *** Option 35 >SUPSYM< activate super symmetry detection and use
 135  SUPSYM = .TRUE.
      GO TO 100
C
C *** Option 36 >THRSSY< threshold for degeneracy when SUPSYM
 136  READ(LUCMD,*) THRSSY
      GO TO 100
C
 9000 CONTINUE

!     set default CI program
      CALL UPCASE(ci_program)
      if(ci_program .eq. 'SIRIUS-CI' .or. 
     &   ci_program .eq. 'LUCITA   ')then
        write(lupri,'(/a,a9)') ' CI program in use: ',
     &                           ci_program
      else   
        write(lupri,'(/3a)')
     &  ' WARNING: Non-valid CI program entered:', ci_program,
     &  ', will use SIRIUS-CI.'
        NWARN = NWARN + 1
        ci_program = 'SIRIUS-CI'
      end if

      IF (THR_REDFAC .GT. 0.0D0) THEN
C         ICHANG = ICHANG + 1
         WRITE (LUPRI,'(3A,1P,D10.2)') '@ INFO ',WORD1,
     &   ' thresholds multiplied with general factor',THR_REDFAC
         THRMC = THRMC*THR_REDFAC
      END IF
      NROOTS = MAX(ISTATE,NROOTS)
      IF ( ISTATE .GT. MAXRTS ) THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,9008) ISTATE, MAXRTS
      ELSE IF ( NROOTS .GT. MAXRTS ) THEN
         NWARN = NWARN + 1
         WRITE (LUPRI,9009) NROOTS, MAXRTS
         NROOTS = MAXRTS
      END IF
      LROOTS = MAX(LROOTS,NROOTS)
C     Fill out rest of IROOT /hjaaj nov 2006
      J = NMCAVER
      DO I = 1,NROOTS
         I1 = I
         DO K = 1,NMCAVER
            IF (IROOT(K) .EQ. I) I1 = 0
         END DO
         IF (I1 .NE. 0) THEN
            J = J + 1
            IROOT(J) = I
         END IF
      END DO
 9008 FORMAT(/' Input ERROR: .STATE/.MCAVER value of,',I3,
     *  ', is greater than internal maximum of,',I3)
 9009 FORMAT(/' OPTINP WARNING: number of simultaneous vectors,',I3,
     *  ', is greater than maximum,',I3,
     *  T11,'The number is reduced to the maximum.')
      RETURN
C
      END
C  /* Deck cicinp */
      SUBROUTINE CICINP(WORD,INPERR,ALLOPT)
C
C 1-Oct 1987 Hans Agren
C
#include "implicit.h"
      PARAMETER (NTABLE = 13)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7
      LOGICAL ALLOPT
C
C Used from common:
C
C GNRINF : THR_REDFAC
C INFINP : DOCINO, DOFCI, DOCISD,FLAG(*), NROOCI, THRCI, THRPWF
C
#include "maxorb.h"
#include "priunit.h"
#include "gnrinf.h"
#include "infinp.h"
C
      DATA TABLE /'.CINO  ','.CIDENS','.CIROOT','.DISKH2','.THRESH',
     *            '.THRPWF','.WEIGHT','.ZEROEL','.STATE ','.FULLCI',
     *            '.CISD  ','.MAX IT','.SPDENS'/
C
C     ***** PROCESS INPUT FOR CICINP *****
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
         CALL UPCASE(WORD)
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102,103,104,105,106,107,108,109,110,
     *                   111,112,113) II
               END IF
 1102       CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               IF (.NOT.ALLOPT)
     &            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/4A/)') ' Keyword ',WORD,
     *         ' not recognized for ',WORD1
         ELSE IF (PROMPT .EQ. '#' .OR. PROMPT .EQ. '!') THEN
            GO TO 100
         ELSE IF (PROMPT .EQ. '*') THEN
            GO TO 8000
         ELSE
            WRITE (LUPRI,'(/3A/2A/)')
     *         ' Keyword ',WORD,' does not begin with',
     *         ' one of the four characters ".*!#" for ',WORD1
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal keyword for "*CI INPUT"')
C
C *** Option 1 >CINO  < Calculate CI natural orbitals.
  101 CONTINUE
      DOCINO   = .TRUE.
      FLAG(4)  = .TRUE.
      FLAG(67) = .TRUE.
      FLAG(68) = .TRUE.
      GO TO 100
C
C *** Option 2 >CIDENS<  Calculate density matrices
  102 CONTINUE
      FLAG(67) = .TRUE.
      GO TO 100
C
C *** Option 3 >CIROOT<  CI calculations with NROOCI roots
  103 CONTINUE
        READ (LUCMD,*) NROOCI
        NROOCI = MAX(NROOCI,ISTACI)
      GO TO 100
C
C *** Option 4 >DISKH2<  Active 2-el MO integrals on disk
  104 CONTINUE
      FLAG(69) = .TRUE.
      GO TO 100
C
C *** Option 5 >THRESH< Threshold for CI energy gradient
  105 READ (LUCMD,*) THRCI
      GO TO 100
C
C *** Option 6 >THRPWF< CI-coefficient threshold for w.f. printout
  106 READ (LUCMD,*) THRPWF
      GO TO 100
C
C *** Option 7 >WEIGHT<  Use energy weighted residuals
  107 CONTINUE
      FLAG(76) = .TRUE.
      GO TO 100
C
C *** Option 8 >ZEROEL<  Zero small elements in CI trial vectors
  108 CONTINUE
      FLAG(77) = .TRUE.
      GO TO 100
C
C *** Option 9 >STATE <  Converge root number ISTATE to threshold
  109 CONTINUE
      READ (LUCMD,*) ISTACI
      NROOCI = MAX(NROOCI,ISTACI)
      GO TO 100
C *** Option 10 >FULLCI<
  110 CONTINUE
      DOFCI  = .TRUE.
      MCTYPE = 1 ! CAS
      GO TO 100
C *** Option 11 >CISD  <
  111 CONTINUE
      DOCISD = .TRUE.
      MCTYPE = 2 ! RAS
      GO TO 100
C
C *** Option 12 >MAX IT< Maximum number of iterations in
C                        CI diagonalization.
  112 CONTINUE
      READ (LUCMD,*) MXCIMA
      GO TO 100
C     !> >SPDENS< set level (1-/2-particle) spin density calculation
  113 CONTINUE
        READ (LUCMD,*) spindens_lvl
      GO TO 100
C
 8000 CONTINUE
      IF (THR_REDFAC .GT. 0.0D0) THEN
C         ICHANG = ICHANG + 1
         WRITE (LUPRI,'(3A,1P,D10.2)') '@ INFO ',WORD1,
     &   ' thresholds multiplied with general factor',THR_REDFAC
         THRCI = THRCI*THR_REDFAC
      END IF
      IF (DOFCI .AND. DOCISD) THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,'(/A)')
     &      'INPUT ERROR, both .CISD and .FULLCI specified'
      END IF
      RETURN
      END
C  /* Deck orbinp */
      SUBROUTINE ORBINP(WORD,INPERR,IRDMO,ALLOPT)
C
C 1-Oct 1987 Hans Agren
C
#include "implicit.h"
#include "priunit.h"
C
C -- local variables and constants
C
      PARAMETER (NTABLE = 18)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7,
     &          REWORD*12
      PARAMETER (MAXNRT = 40)
      DIMENSION INORO(MAXNRT),NREOR(8),NNOR(8)
      LOGICAL ALLOPT
C
C Used from common:
C
C INFINP : FLAG,IMOORD,CMAXMO
C INFORB : NBAS
C HUCKEL : DOHUCKEL
C
#include "maxorb.h"
#include "infinp.h"
#include "inforb.h"
#include "huckel.h"
#include "r12int.h"
#include "infloc.h"
C
      DATA TABLE /'XXXXXXX','.5D7F9G','.DELETE','.MOSTAR','.FREEZE',
     *            '.GRAM-S','.SYMMET','.PUNCHI','.PUNCHO','.REORDE',
     *            '.FROZEN','.AO DEL','.CMOMAX','xxxxxxx','xxxxxxx',
     *            'xxxxxxx','.LOCALI','.R12ORB'/
C
      INLOC = 0
      INDEL = 0
C
C     ***** PROCESS INPUT FOR ORBINP *****
C
      WORD1 = WORD
      IF(ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                      LUPRI)
  100    READ (LUCMD, '(A7)') WORD
         CALL UPCASE(WORD)
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102,103,104,105,106,107,108,109,110,
     *                   111,112,113,114,115,116,117,118), II
               END IF
 1102       CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               IF (.NOT.ALLOPT)
     &           CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                      LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/4A/)') ' Keyword ',WORD,
     *         ' not recognized for ',WORD1
         ELSE IF (PROMPT .EQ. '#' .OR. PROMPT .EQ. '!') THEN
            GO TO 100
         ELSE IF (PROMPT .EQ. '*') THEN
            GO TO 9000
         ELSE
            WRITE (LUPRI,'(/3A/2A/)')
     *         ' Keyword ',WORD,' does not begin with',
     *         ' one of the four characters ".*!#" for ',WORD1
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal keyword for "*ORBITAL INPUT"')
C
C *** Option 1 >XXXXXX<
  101 CONTINUE
      GO TO 100
C
C *** Option 2 >5D7F9G<  Delete unwanted gaussian components
C                        (s in d, p in f, s and d in g)
  102 CONTINUE
      IF (INDEL .NE. 0) THEN
        WRITE(LUPRI,'(/,A,/)') ' .5D7F9G too late in input stream'
        CALL QUIT(' .5D7F9G too late in input stream')
      ENDIF
C
C     check KDEL option: delete 3d(3s),4f(4p),... components in CMO
C
      IF (KDEL .NE. 0) THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,'(/A,I5,A/A)')
     *      ' FATAL ERROR Sirius ORBINP, KDEL =',KDEL,' at .5D7F9G',
     *      ' KDEL non-zero means inconsistency in deleting orbitals.'
      ELSE
         KDEL  = 1
C        ... to specify this option and that DELETE is now illegal.
         CALL CMODEL(CMO,.FALSE.)
C        CALL CMODEL(CMO,GETCMO)
C        ... call cmodel to get NORB(*); CMO will not be used.
C            CMODEL will be called later to get CMO, if MOSTART 8.
      END IF
      GO TO 100
C
C *** Option 3 >DELETE<  Delete orbitals
 103  CONTINUE
      IF (INDEL .NE. 0) THEN
        WRITE(LUPRI,'(/,A,/)') ' .DELETE too late in input stream'
        CALL QUIT(' .DELETE too late in input stream')
      ENDIF
      READ(LUCMD,*) (NORB(I),I=1,NSYM)
      DO 1205 I = 1,NSYM
         NORB(I) = NBAS(I) - NORB(I)
 1205 CONTINUE
      IF (KDEL .NE. 0) THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,'(/A,I5,A/A)')
     *      ' FATAL ERROR Sirius ORBINP, KDEL =',KDEL,' at .DELETE',
     *      ' KDEL non-zero means inconsistency in deleting orbitals.'
      ELSE
         KDEL  = -1
C        ... to specify that 5D7F9G is now illegal.
      END IF
      GO TO 100
C
C *** Option 4 >MOSTAR<   Molecular orbital input
 104  CONTINUE
        READ (LUCMD,'(A)') REWORD
        CALL UPCASE(REWORD)
          IF (INDEX(REWORD,'FORM18') .GT. 0) THEN
             IRDMO = 3
          ELSE IF (INDEX(REWORD,'FORM12') .GT. 0) THEN
             IRDMO = 4
          ELSE IF (INDEX(REWORD,'OLDORB') .GT. 0) THEN
             IRDMO = 7
          ELSE IF (INDEX(REWORD,'H1DIAG') .GT. 0) THEN
             IRDMO = 8
             DOHUCKEL = .FALSE.
C         ... make sure Huckel guess is disabled,
C             although code in GNRLIN in dalton.F also takes care of that
          ELSE IF (INDEX(REWORD,'NEWORB') .GT. 0) THEN
             IRDMO = 9
          ELSE IF (INDEX(REWORD,'SIRIFC') .GT. 0) THEN
             IRDMO = 10
          ELSE IF (INDEX(REWORD,'LUCIA') .GT. 0) THEN
             IRDMO = 11
          ELSE IF (INDEX(REWORD,'HUCKEL') .GT. 0 .OR.
     &             INDEX(REWORD,'EHT')    .GT. 0) THEN
             IRDMO = 8
             IF (.NOT.DOHUCKEL) THEN
                WRITE (LUPRI,'(/A)') ' INFO: .MOSTART HUCKEL changed'//
     &             ' to H1DIAG because Huckel guess not available.'
             END IF
          ELSE IF (INDEX(REWORD,'EWMO') .GT. 0) THEN
             IRDMO = 8
             IF (.NOT.DOHUCKEL) THEN
                WRITE (LUPRI,'(/A)') ' INFO: .MOSTART HUCKEL changed'//
     &             ' to H1DIAG because Huckel guess not available.'
             END IF
          ELSE
             WRITE(LUPRI,'(3A/A)') '"',REWORD,
     *          '" is incorrect input after ".MOSTART"',
     *          ' Select from .FORM12, .FORM18, .H1DIAG, .HUCKEL,',
     *          ' .EWMO, .SIRIFC, .OLDORB, and .NEWORB'
             CALL QUIT(' INPUT ERROR IN ORBINP')
          END IF
          IF (IRDMO .EQ. 8) THEN
             FLAG(13) = .FALSE.
          ELSE
             FLAG(13) = .TRUE.
          END IF
          IF (IRDMO .NE. 8) DOHUCKEL = .FALSE.
C         ... make sure Huckel guess is disabled,
C             although code in GNRLIN in dalton.F takes care of that
      GO TO 100
C
C *** Option 5 >FREEZE< Orbitals not to be rotated.
  105 CONTINUE
         IF (FLAG(33)) CALL QUIT(
     &      'ORBINP ERROR, .FREEZE and .AO DELETE both specified.')
C        set INDEL=2 to flag that DELETE cannot be specified any more
         INDEL = 2
         READ(LUCMD,*)(NNOR(I),I=1,NSYM)
         IJ=0
         DO 250 I=1,NSYM
            NORI=NNOR(I)
            IF(NORI.EQ.0)GO TO 251
            LNOROT = .TRUE.
            IF (NORI.GT.MAXNRT) THEN
               WRITE(LUPRI,2311) NORI,I,MAXNRT
               INPERR = INPERR + 1
               READ(LUCMD,*) (IDUM,J=1,NORI)
               GO TO 251
            END IF
            READ(LUCMD,*)(INORO(J),J=1,NORI)
            DO 252 J=1,NORI
               II=IJ+INORO(J)
               NOROT(II)=1
  252       CONTINUE
  251       IJ=IJ+NORB(I)
  250    CONTINUE
      GO TO 100
 2311 FORMAT(/' Sirius input ERROR for .FREEZE,',I5,' orbitals speci'//
     *    'fied not to be rotated in symmetry',I2/20X,'Maximum is',I3)
C
C *** Option 6 >GRAM-S<  Gram-Schmidt orthogonalization
  106 CONTINUE
        FLAG(32) = .FALSE.
      GO TO 100
C
C *** Option 7 >SYMMET<  Symmetric orthogonalization
  107 CONTINUE
        FLAG(32) = .TRUE.
      GO TO 100
C
C *** Option 8 >PUNCHI<  Punch input orbitals
  108 CONTINUE
        FLAG(31) = .TRUE.
      GO TO 100
C
C *** Option 9 >PUNCHO<  Punch output orbitals
  109 CONTINUE
        FLAG(18) = .TRUE.
      GO TO 100
C
C *** Option 10 >REORDE<  reordering of molecular orbitals
  110 LMOORD = .TRUE.
      INDEL = 1
C     ... to signal that DELETE cannot be specified any more
      IF (FLAG(33)) THEN
         CALL QUIT(
     &      'ORBINP ERROR, .REORDER and .AO DELETE both specified.')
      END IF
      READ(LUCMD,*) (NREOR(I),I=1,NSYM)
      DO 801 I = 1,MXCORB
         IMOORD(I) = I
  801 CONTINUE
      IORBI = 0
      DO 810 I = 1,NSYM
         IF (NREOR(I) .GT. 0) THEN
            NR2 = 2*NREOR(I)
            IF (NR2.GT.MAXNRT) THEN
               WRITE (LUPRI,*) 'FATAL ERROR for ".REORDER MO"'
               WRITE (LUPRI,*) 'Internal limit of',MAXNRT/2,
     *         ' orbitals to be reordered exceeded in symmetry',I
               CALL QUIT('ERROR in ORBINP, internal limit exceeded')
            END IF
            READ(LUCMD,*) (INORO(J),J=1,NR2)
C           ... this line corresponds to
C               READ(LUCMD,*) (IMONEW(J,I),IMOOLD(J,I),J=1,NREOR(I))
C               in the MOTECC-90 Input/Output documentation
            DO 805 J = 1,NR2,2
              IMOORD(IORBI+INORO(J)) = IORBI+INORO(J+1)
  805       CONTINUE
         END IF
         IORBI = IORBI + NORB(I)
  810 CONTINUE
      GO TO 100
C
C *** Option 11 >FROZEN<
  111 CONTINUE
      READ (LUCMD,*) (LOCFRO(I),I=1,NSYM)
      GO TO 100
C
C *** Option 12 >AO DEL< Delete MO:s based on AO overlap matrix.
  112 IF (INDEL .NE. 0 .OR. KDEL .LT. 0) THEN
         WRITE(LUPRI,*) 'ERROR, ".AO DELETE" cannot be defined together'
         WRITE(LUPRI,*) 'with ".DELETE",".REORDE", or ".FREEZE"'
         CALL QUIT('ERROR in ORBINP, conflict with ".AO DELETE"')
      END IF
      INDEL = 3
      FLAG(33) = .TRUE.
      READ(LUCMD,*) THROVL
      GO TO 100
C
C *** Option 13 >CMOMAX< max allowed absolute value in CMO
  113 CONTINUE
      READ (LUCMD,*) CMAXMO
      GO TO 100
C
C *** Option 14 >
  114 CONTINUE
      GO TO 100
C
C *** Option 15 >
  115 CONTINUE
      GO TO 100
C
C *** Option 16 >
  116 CONTINUE
      GO TO 100
C
*** Option 17 >LOCALI< localize molecular orbitals
  117 INLOC = INLOC + 1
      IF (INLOC .GT. 1) THEN
         INPERR = INPERR + 1
         WRITE(LUPRI,'(/2A)')
     &   WORD1,' INPUT ERROR: .LOCALI specified more than once'
      END IF
      READ (LUCMD,'(A)') REWORD
      CALL UPCASE(REWORD)
      IF (INDEX(REWORD,'BOYLOC') .GT. 0) THEN
         BOYORB = .TRUE.
      ELSE IF (INDEX(REWORD,'SELECT') .GT. 0) THEN
         BOYSEL = .TRUE.
         CALL IZERO(NBOYS,8)
         CALL IZERO(BOYSORB,MXBOYS*8)
         READ(LUCMD,*)(NBOYS(I),I=1,NSYM)
         DO I=1,NSYM
            IF (NBOYS(I).GT.0) THEN
               READ(LUCMD,*)(BOYSORB(J,I),J=1,NBOYS(I))
            END IF
         END DO
      ELSE IF (INDEX(REWORD,'PIPLOC') .GT. 0) THEN
         PIPORB = .TRUE.
         CALL QUIT('Pipek-Mezey localization not yet implemented')
      ELSE
         WRITE(LUPRI,'(3A/A)') '"',REWORD(1:7),'" is incorrect input'//
     *         ' after ".LOCALI"',' Select " BOYLOC" or " PIPLOC"'
         INPERR = INPERR + 1
      END IF
      GO TO 100
C
C *** Option 18 >R12ORB<
  118 READ (LUCMD,*) (NRXR12(I), I=1,NSYM)
      GO TO 100
C
C ---
C
 9000 CONTINUE
      RETURN
      END
C  /* Deck civinp */
      SUBROUTINE CIVINP(WORD,INPERR,ALLOPT)
C
C 1-Oct 1987 Hans Agren
C
#include "implicit.h"
      PARAMETER (NTABLE = 5)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7
      LOGICAL ALLOPT
C
C Used from common:
C
C INFINP : FLAG(*), ICI0,
C
#include "maxorb.h"
#include "priunit.h"
#include "infinp.h"
C
      DATA TABLE /'.PLUS C','.SELECT','.STARTH','.ABACUS','.STARTO'/
C
C     ***** PROCESS INPUT FOR CIVINP *****
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
         CALL UPCASE(WORD)
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102,103,104,105), II
               END IF
 1102       CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               IF (.NOT.ALLOPT)
     &            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/4A/)') ' Keyword ',WORD,
     *         ' not recognized for ',WORD1
         ELSE IF (PROMPT .EQ. '#' .OR. PROMPT .EQ. '!') THEN
            GO TO 100
         ELSE IF (PROMPT .EQ. '*') THEN
            RETURN
         ELSE
            WRITE (LUPRI,'(/3A/2A/)')
     *         ' Keyword ',WORD,' does not begin with',
     *         ' one of the four characters ".*!#" for ',WORD1
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal keyword for "*CI VECTOR"')
C
C *** Option 1 >PLUS C< Choose plus combination of degenerate roots
  101 CONTINUE
        FLAG(58) = .TRUE.
      GO TO 100
C
C *** Option 2 >SELECT<  SELECT -ICI0 as start configuration
  102 CONTINUE
         READ(LUCMD,*) ICI0
         ICI0  = - ICI0
      GO TO 100
C
C *** Option 3 >STARTH<  Conmpute start vector from H-diagonal
  103 CONTINUE
         ICI0  = 1
      GO TO 100
C
C *** Option 4 >ABACUS<  Geometry walk, new geometry ('GEOSAVE')
  104 CONTINUE
         ICI0  = 6
      GO TO 100
C
C *** Option 5 >STARTO<  Start from old CI-vector stored on unit 21
  105 CONTINUE
         ICI0  = 4
      GO TO 100
C
C
      END
C  /* Deck stpinp */
      SUBROUTINE STPINP(WORD,INPERR,ALLOPT)
C
C 1-Oct 1987 Hans Agren
C
#include "implicit.h"
C
C -- local variables and constants
C
      PARAMETER (NTABLE = 16)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7
      LOGICAL ALLOPT
C
C Used from common:
C
C INFINP : FLAG(*)
C INFOPT : BET*,RT*,RAT*,STP*
C
#include "maxorb.h"
#include "priunit.h"
#include "infinp.h"
#include "infopt.h"
C
      DATA TABLE /'.DAMPIN','.MAX DA','.MIN DA','.NO EXT','.GOOD R',
     *            '.MIN RA','.REJECT','.TRUST ','.TOLERA','.INCREM',
     *            '.MAX ST','.DECREM','.TIGHT ','.THQMIN','.THQLIN',
     *            '.THQKVA'/
C
C     ***** PROCESS INPUT FOR STPINP  *****
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
         CALL UPCASE(WORD)
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102,103,104,105,106,107,108,109,110,
     *                   111,112,113,114,115,116), II
               END IF
 1102       CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               IF (.NOT.ALLOPT)
     &            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/4A/)') ' Keyword ',WORD,
     *         ' not recognized for ',WORD1
         ELSE IF (PROMPT .EQ. '#' .OR. PROMPT .EQ. '!') THEN
            GO TO 100
         ELSE IF (PROMPT .EQ. '*') THEN
            RETURN
         ELSE
            WRITE (LUPRI,'(/3A/2A/)')
     *         ' Keyword ',WORD,' does not begin with',
     *         ' one of the four characters ".*!#" for ',WORD1
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal keyword for "*STEP CONTROL"')
C
C *** Option 1 >DAMPIN <  Initial value of damping (BETA)
  101 READ(LUCMD,*) BETA
      GO TO 100
C
C *** Option 2 >MAX DA<  Maximum damping value
  102 READ(LUCMD,*) BETMAX
      GO TO 100
C
C *** Option 3 >MIN DA<  Minimum damping value
  103 READ(LUCMD,*) BETMIN
      GO TO 100
C
C *** Option 4 >NO TER<  Skip extra termination tests
  104 FLAG(36) = .TRUE.
      GO TO 100
C
C *** Option 5 >GOOD R<  "Good" trust radius
  105 READ(LUCMD,*) RATGOD
      GO TO 100
C
C *** Option 6 >MIN RA<  "Bad" trust radius
  106 READ(LUCMD,*) RATMIN
      GO TO 100
C
C *** Option 7 >REJECT<  Ratio predicted/actual step that rejects step
  107 READ(LUCMD,*) RATREJ
      GO TO 100
C
C *** Option 8 >TRUST <  Initial trust radius
  108 READ(LUCMD,*) RTRUST
      GO TO 100
C
C *** Option 9 >TOLERA<  Tolerable trust radius
  109 READ(LUCMD,*) RTTOL
      GO TO 100
C
C *** Option 10 >INCREM<  Increment factor on trust radius
  110 READ (LUCMD,*) STPINC
      GO TO 100
C
C *** Option 11 >MAX ST<  Maximum acceptable step length
  111 READ (LUCMD,*) STPMAX
      GO TO 100
C
C *** Option 12 >DECREM<  Decrement factor on trust radius
  112 READ (LUCMD,*) STPRED
      GO TO 100
C
C *** Option 13 >TIGHT <  Tight step contol for ground state also
  113 FLAG(35) = .TRUE.
      GO TO 100
C
C *** Option 14 >THQMIN<
  114 READ (LUCMD,*) THQMIN
      GO TO 100
C
C *** Option 15 >THQLIN<
  115 READ (LUCMD,*) THQLIN
      GO TO 100
C
C *** Option 16 >THQKVA<
  116 READ (LUCMD,*) THQKVA
      GO TO 100
C
      END
C  /* Deck popinp */
      SUBROUTINE POPINP(WORD,INPERR,ALLOPT)
C
C 1-Oct 1987 Hans Agren
C 10 Aug 89 hjaaj: default to nothing (flag(71:74) = .true.)
C                  any population analysis option sets flag(73) false
C                  new option: .PRINT
C 13-sep-90 hjaaj: IPRMUL = MAX(IPRMUL,xxx)
C
#include "implicit.h"
C
C -- local variables and constants
C
      PARAMETER (NTABLE = 11)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7
      LOGICAL ALLOPT
C
C Used from common:
C
C INFINP : FLAG(*)
C INFPRI : IPRMUL
C
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "infinp.h"
#include "infpri.h"
#include "infpop.h"
C
      DATA TABLE /'.ALL   ','.DIPMOM','.GROSSA','.GROSSM','.MULLIK',
     *            '.NETALL','.NETMO ','.QUADRP','.VIRIAL','.PRINT ',
     &            '.HJAAJP'/
C
C External flags
      FLAG( 6) = .TRUE.
      FLAG(71) = .TRUE.
      FLAG(72) = .TRUE.
      FLAG(73) = .TRUE.
      FLAG(74) = .TRUE.
C
C     ***** PROCESS INPUT FOR POPINP *****
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
         CALL UPCASE(WORD)
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102,103,104,105,106,107,108,109,110,
     &                   111), II
               END IF
 1102       CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               IF (.NOT.ALLOPT)
     &            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/4A/)') ' Keyword ',WORD,
     *         ' not recognized for ',WORD1
         ELSE IF (PROMPT .EQ. '#' .OR. PROMPT .EQ. '!') THEN
            GO TO 100
         ELSE IF(PROMPT .EQ. '*') THEN
            RETURN
         ELSE
            WRITE (LUPRI,'(/3A/2A/)')
     *         ' Keyword ',WORD,' does not begin with',
     *         ' one of the four characters ".*!#" for ',WORD1
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal keyword for "*POPULATION ANALYSIS"')
C
C *** Option 1 >ALL   <  Do all options
  101 FLAG(71) = .FALSE.
      FLAG(72) = .FALSE.
      FLAG(73) = .FALSE.
      FLAG(74) = .FALSE.
      IPRMUL = MAX(IPRMUL,5)
      GO TO 100
C
C *** Option 2 >DIPMOM<  Dipol moments
  102 FLAG(71) = .FALSE.
      GO TO 100
C
C *** Option 3 >GROSSA<  All gross population analysis
  103 FLAG(73) = .FALSE.
      IPRMUL = MAX(IPRMUL,2)
      GO TO 100
C
C *** Option 4 >GROSSM<  Gross MO population analysis
  104 FLAG(73) = .FALSE.
      IPRMUL = MAX(IPRMUL,1)
      GO TO 100
C
C *** Option 5 >MULLIK<  Mulliken population analysis
  105 FLAG(73) = .FALSE.
      IPOPTYP  = 1
      GO TO 100
C
C *** Option 6 >NETALL<  All net population analysis
  106 FLAG(73) = .FALSE.
      IPRMUL = MAX(IPRMUL,5)
      GO TO 100
C
C *** Option 7 >NETMO <  Net MO population analysis
  107 FLAG(73) = .FALSE.
      IPRMUL = MAX(IPRMUL,3)
      GO TO 100
C
C *** Option 8 >QUADRP<  Quadrupole moments
  108 FLAG(72) = .FALSE.
      GO TO 100
C
C *** Option 9 >VIRIAL<  Virial analysis
  109 FLAG(74) = .FALSE.
      GO TO 100
C
C *** Option 10 >PRINT <  Print level
  110 READ(LUCMD,*) IPRXXX
      IPRMUL = MAX(IPRMUL,IPRXXX)
      if (iprxxx .eq. 7) IPOPTYP = 2
      GO TO 100
C
C *** Option 11 >HJAAJP<  HJAaJ population analysis
  111 FLAG(73) = .FALSE.
      IPOPTYP = 3
      GO TO 100
C
      END
C  /* Deck trainp */
      SUBROUTINE TRAINP(WORD,INPERR,ALLOPT)
C
C 1-Oct 1987 Hans Agren
C
#include "implicit.h"
C
C -- local variables and constants
C
      PARAMETER (NTABLE = 8)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7
      LOGICAL ALLOPT
C
C Used from common:
C
C GNRINF : THR_REDFAC
C INFINP : FLAG(*)
C INFDIM : MWORK
C INFTRA : ITRLVL, THRP
C
#include "maxorb.h"
#include "priunit.h"
#include "gnrinf.h"
#include "infinp.h"
#include "infdim.h"
#include "inftra.h"
C
      DATA TABLE /'.RESIDE','.FINAL ','.LEVEL ','.OLD TR','.PRINT ',
     *            '.P THRE','xxxxxxx','xxxxxxx'/
C
C     ***** PROCESS INPUT FOR TRAINP *****
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
         CALL UPCASE(WORD)
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102,103,104,105,106,107,108), II
               END IF
 1102       CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               IF (.NOT.ALLOPT)
     &            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/4A/)') ' Keyword ',WORD,
     *         ' not recognized for ',WORD1
         ELSE IF (PROMPT .EQ. '#' .OR. PROMPT .EQ. '!') THEN
            GO TO 100
         ELSE IF (PROMPT .EQ. '*') THEN
            GO TO 8000
         ELSE
            WRITE (LUPRI,'(/3A/2A/)')
     *         ' Keyword ',WORD,' does not begin with',
     *         ' one of the four characters ".*!#" for ',WORD1
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal keyword for "*TRANSFORMATION INPUT"')
C
C *** Option 1 >RESIDE<  Desired length of work array in transf.
  101 READ(LUCMD,*) MWORK
      GO TO 100
C
C *** Option 2 >FINAL < Final integral transformation level
  102 READ (LUCMD,*) ITRFIN
      GO TO 100
C
C *** Option 3 >LEVEL < Transformation level
  103 READ(LUCMD,*) ITRLVL
      GO TO 100
C
C *** Option 4 >OLD TR<  Use existing transformed integrals
  104 FLAG(14) = .TRUE.
      GO TO 100
C
C *** Option 5 >PRINT <  Print level in transformation
  105 READ(LUCMD,*) IPRTRA
      GO TO 100
C
C *** Option  106 >P THRE<
  106 READ(LUCMD,*) THRP
      GO TO 100
C
C *** Option  107 >      <
  107 CONTINUE
      GO TO 100
C
C *** Option  108 >      <
  108 CONTINUE
      GO TO 100
C
 8000 CONTINUE
      IF (THR_REDFAC .GT. 0.0D0) THEN
C         ICHANG = ICHANG + 1
         WRITE (LUPRI,'(3A,1P,D10.2)') '@ INFO ',WORD1,
     &   ' thresholds multiplied with general factor',THR_REDFAC
         THRP = THRP*THR_REDFAC
      END IF
      END
C  /* Deck updinp */
      SUBROUTINE UPDINP(WORD,INPERR,ALLOPT)
C
C 1-Oct 1987 Hans Agren
C
#include "implicit.h"
C
C -- local variables and constants
C
      PARAMETER (NTABLE = 2)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7
      LOGICAL ALLOPT
C
C Used from common:
C
C INFINP : FLAG(*)
C
#include "maxorb.h"
#include "priunit.h"
#include "infinp.h"
C
      DATA TABLE /'.FIXED ','.INVERS'/
C
      FLAG(37) = .FALSE.
C
C     ***** PROCESS INPUT FOR UPDINP *****
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
         CALL UPCASE(WORD)
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102), II
               END IF
 1102       CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               IF (.NOT.ALLOPT)
     &            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/4A/)') ' Keyword ',WORD,
     *         ' not recognized for ',WORD1
         ELSE IF (PROMPT .EQ. '#' .OR. PROMPT .EQ. '!') THEN
            GO TO 100
         ELSE IF (PROMPT .EQ. '*') THEN
            RETURN
         ELSE
            WRITE (LUPRI,'(/3A/2A/)')
     *         ' Keyword ',WORD,' does not begin with',
     *         ' one of the four characters ".*!#" for ',WORD1
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal keyword for "*UPDATE"')
C
C *** Option 1 >FIXED <  Use fixed hessian
  101 FLAG(56) = .TRUE.
      GO TO 100
C
C *** Option 2 >INVERS<  Use inverse start inverse hessian diagonal
  102 FLAG(57) = .TRUE.
      GO TO 100
C
      END
C  /* Deck solinp */
      SUBROUTINE SOLINP(WORD,INPERR,ALLOPT)
C
C 1-Oct 1987 Hans Agren
C Revised Nov 1989 hjaaj
C
#include "implicit.h"
C
C -- local variables and constants
C
      PARAMETER (NTABLE = 6, D1 = 1.0D0)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7
      LOGICAL ALLOPT
C
C Used from common:
C
C INFINP : FLAG(*), EPSOL, RSOL, NLMSOL, LSOLMX, EPSTAT,
C          INERSI, INERSF, EPPN
C
C EPPN epsilon value in non-equilibrium solvent response calculation
C
C
C INFPRI : IPRSOL
C
#include "maxorb.h"
#include "priunit.h"
#include "infinp.h"
#include "infpri.h"
C
      DATA TABLE /'.CAVITY','.DIELEC','.MAX L ','.PRINT ',
     *            '.INERSI','.INERSF'/
C
C     ***** PROCESS INPUT FOR SOLINP *****
C
      INPEPS = 0
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
         CALL UPCASE(WORD)
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102,103,104,105,106), II
               END IF
 1102       CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               IF (.NOT.ALLOPT)
     &            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/4A/)') ' Keyword ',WORD,
     *         ' not recognized for ',WORD1
         ELSE IF (PROMPT .EQ. '#' .OR. PROMPT .EQ. '!') THEN
               GO TO 100
         ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 9999
         ELSE
            WRITE (LUPRI,'(/3A/2A/)')
     *         ' Keyword ',WORD,' does not begin with',
     *         ' one of the four characters ".*!#" for ',WORD1
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal keyword for "*SOLVENT"')
C
C *** Option 1 >CAVITY<  Cavity radius
  101 READ(LUCMD,*) RSOLAV
      RSOL(1) = RSOLAV
      RSOL(2) = RSOLAV
      RSOL(3) = RSOLAV
      GO TO 100
C
C *** Option 2 >DIELEC<  Dielectric constant
  102 READ(LUCMD,*) EPSOL
      EPSTAT = EPSOL
      INPEPS = INPEPS + 1
      IF (EPSOL.LT.D1) THEN
         INPERR = INPERR + 1
         WRITE(LUPRI,'(/A/A,1P,2D10.2)')
     &      ' INPUT ERROR: Invalid dielectric constant for .DIELEC',
     &      ' Input value :',EPSOL
      END IF
      GO TO 100
C
C *** Option 3 >MAX L <  Maximum L quantum number
  103 READ(LUCMD,*) LSOLMX
      NLMSOL = (LSOLMX+1)*(LSOLMX+1)
      GO TO 100
C
C *** Option 4 >PRINT <  Print level
  104 READ(LUCMD,*) IPRSOL
      GO TO 100
C
C *** Option 5 >INERSI<  Initial state defining inertial pol.
  105 INERSI = .TRUE.
      INPEPS = INPEPS + 1
      READ(LUCMD,*) EPSOL, EPPN
      EPSTAT = EPSOL
      IF (EPSOL.LT.D1 .OR. EPPN.LT.D1) THEN
         INPERR = INPERR + 1
         WRITE(LUPRI,'(/A/A,1P,2D10.2)')
     &      ' INPUT ERROR: Invalid dielectric constants for .INERSI',
     &      ' Static and dynamic input values :',EPSOL,EPPN
      END IF
      GO TO 100
C
C *** Option 6 >INERSF<  Final state with inertial pol. from
C                        previous calc. with ".INERSI"
  106 INERSF = .TRUE.
      INPEPS = INPEPS + 1
      READ(LUCMD,*) EPSTAT, EPSOL
      IF (EPSOL.LT.D1 .OR. EPSTAT.LT.D1) THEN
         INPERR = INPERR + 1
         WRITE(LUPRI,'(/A/A,1P,2D10.2)')
     &      ' INPUT ERROR: Invalid dielectric constants for .INERSF',
     &      ' Static and dynamic input values :',EPSTAT,EPSOL
      END IF
      GO TO 100
C
 9999 CONTINUE
      IF (INPEPS .GT. 1) THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,'(/1X,A/A)') WORD1//' input ERROR',
     *   ' Only one of ".DIELEC", ".INERSI", and ".INERSF"'//
     *   ' may be specified.'
      END IF
C     ... set flag for solvent calc.
      IF (INPEPS .EQ. 1) FLAG(16) = .TRUE.
      RETURN
      END
C  /* Deck priinp */
      SUBROUTINE PRIINP(WORD,INPERR,ALLOPT)
C
C 1-Oct 1987 Hans Agren
C
#include "implicit.h"
C
C -- local variables and constants
C
      PARAMETER (NTABLE = 17)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7
      PARAMETER (NPDIM = 50)
      INTEGER   NPPATH(NPDIM)
      LOGICAL ALLOPT
C
C Used from common:
C
C INFINP : FLAG(*),
C INFPRI : LU*, IPR*,
C INFTRA : IPRTRA
C
#include "maxorb.h"
#include "priunit.h"
#include "infinp.h"
#include "infpri.h"
#include "inftra.h"
C
      DATA TABLE /'.CMOPRI','.IPRKAP','.IPRDIA','.IPRSIG','.IPRDNS',
     *            '.IPRSOL','.RESFIL','.PRINTL','.PRINTF','.IPRCIX',
     *            '.THRPWF','.IPRCNO','.NOSUMM','.CANONI','.IPRAVE',
     *            '.IPRFCK','.POPPRI'/
C
C     ***** PROCESS INPUT FOR PRIINP *****
C
      IF (WORD .NE. '*PRINT ') THEN
         WRITE(LUPRI,*) 'ERROR, PRIINP call with ',WORD,
     *   ' instead of "*PRINT "'
         CALL QUIT('PROGRAM ERROR SIRIUS, illegal call of PRIINP')
      END IF
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
         CALL UPCASE(WORD)
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102,103,104,105,106,107,108,109,110,
     *                   111,112,113,114,115,116,117), II
               END IF
 1102       CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               IF (.NOT.ALLOPT)
     &            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     *                       LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/4A/)') ' Keyword ',WORD,
     *         ' not recognized for ',WORD1
         ELSE IF (PROMPT .EQ. '#' .OR. PROMPT .EQ. '!') THEN
            GO TO 100
         ELSE IF (PROMPT .EQ. '*') THEN
            RETURN
         ELSE
            WRITE (LUPRI,'(/3A/2A/)')
     *         ' Keyword ',WORD,' does not begin with',
     *         ' one of the four characters ".*!#" for ',WORD1
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal keyword for "*PRINT "')
C
C *** Option 1 >CMOPRI<  Print all MO coefficients
C
  101 CONTINUE
      CMOPRI = .TRUE.
      GO TO 100
C
C *** Option 2 >IPRKAP<  Prinout level for NEXKAP section.
C                        (Optimal Orbital Trial Vectors).
C
  102 READ(LUCMD,*) IPRKAP
      GO TO 100
C
C *** Option 3 >IPRDIA<
  103 READ(LUCMD,*) IPRDIA
      GO TO 100
C
C *** Option 4 >IPRSIG<
  104 READ(LUCMD,*) IPRSIG
      GO TO 100
C
C *** Option 5 >IPRDNS<
  105 READ(LUCMD,*) IPRDNS
      GO TO 100
C
C *** Option 6 >IPRSOL<
  106 READ(LUCMD,*) IPRSOL
      GO TO 100
C
C *** Option 7 >RESFIL<  Write sirius result file separate from LUPRI
 107  LUW4 = -8888
C
C     Note, the unit number for LUW4 will be reassigned when the file is opened
C     in the main input-reading routine.
C
      GO TO 100
C
C *** Option 8 >PRINTL<  Print levels on units LUPRI and LUW4
  108 READ(LUCMD,*) IPRI6,IPRI4
      IPRIX = MOD(IPRI4,100)
      IPRIX = MIN(IPRIX,NPFLAG)
      DO 280 I = 1,IPRIX
         P4FLAG(I) = .TRUE.
  280 CONTINUE
      IPRIX = MOD(IPRI6,100)
      IPRIX = MIN(IPRIX,NPFLAG)
      DO 290 I = 1,IPRIX
         P6FLAG(I) = .TRUE.
  290 CONTINUE
      GO TO 100
C
C *** Option 9 >PRINTF<  Individual print flag settings
  109 READ(LUCMD,*) NUM6, NUM4
         IF (NUM6 .GT. 0) THEN
            READ (LUCMD,*) (NPPATH(I), I=1,NUM6)
            DO 298 I = 1,NUM6
               IF (IABS(NPPATH(I)) .GT. NPFLAG) THEN
                  INPERR = INPERR + 1
                  WRITE (LUPRI,'(/A,I10)')
     *               ' ERROR .PRINTF, illegal pri6 value:',NPPATH(I)
               ELSE IF (NPPATH(I).LT.0) THEN
                  P6FLAG(-NPPATH(I)) = .FALSE.
               ELSE IF (NPPATH(I).GT.0) THEN
                  P6FLAG(NPPATH(I))  = .TRUE.
               ELSE
                  GO TO 299
               END IF
  298       CONTINUE
  299       CONTINUE
         END IF
         IF (NUM4 .GT. 0) THEN
            READ (LUCMD,*) (NPPATH(I), I=1,NUM4)
            DO 288 I = 1,NUM4
               IF (IABS(NPPATH(I)) .GT. NPFLAG) THEN
                  INPERR = INPERR + 1
                  WRITE (LUPRI,'(/A,I5)')
     *               ' ERROR .PRINTF, illegal pri4 value:',NPPATH(I)
               ELSE IF (NPPATH(I).LT.0) THEN
                  P4FLAG(-NPPATH(I)) = .FALSE.
               ELSE IF (NPPATH(I).GT.0) THEN
                  P4FLAG(NPPATH(I))  = .TRUE.
               ELSE
                  GO TO 289
               END IF
  288       CONTINUE
  289       CONTINUE
         END IF
C
      GO TO 100
C
C *** Option 10 >IPRCIX<
  110 READ(LUCMD,*) IPRCIX
      GO TO 100
C
C *** Option 11 >THRPWF< Threshold for printout of wave function CI
C                        coefficients.
  111 READ(LUCMD,*) THRPWF
      GO TO 100
C
C *** Option 12 >IPRCNO< Printout level for Canonical/Natural orbital
C                        section.
  112 READ(LUCMD,*) IPRCNO
      GO TO 100
C
C *** Option 13 >NOSUMM<  No final summary
C
  113 FLAG(5) = .FALSE.
      GO TO 100
C
C *** Option 14 >CANONI< in summary, if converged then call sircan for
C                        generation of canonical/natural orbitals
C
  114 FLAG(3) = .TRUE.
      GO TO 100
C
C *** Option 15 >IPRAVE< Printout level for super sym. section
  115 READ(LUCMD,*) IPRAVE
      GO TO 100
C
C *** Option 16 >IPRFCK< Printout level for super sym. section
  116 READ(LUCMD,*) IPRFCK
      GO TO 100
C
C *** Option 17 >POPPRInt< print popana each iter for first LIM_POPPRI atoms
  117 READ(LUCMD,*) LIM_POPPRI
      GO TO 100
C
      END
C  /* Deck auxinp */
      SUBROUTINE AUXINP(WORD,INPERR,ALLOPT)
C
C 1-Oct 1987 Hans Agren
C
#include "implicit.h"
C
C -- local variables and constants
C
      PARAMETER (NTABLE = 2)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7
      LOGICAL ALLOPT
C
C Used from common:
C
C INFINP :
C INFTAP : LUSUPM
C
#include "maxorb.h"
#include "priunit.h"
#include "infinp.h"
#include "inftap.h"
C
      DATA TABLE /'XXXXXXX','.NOSUPM'/
C
C     ***** PROCESS INPUT FOR AUXINP *****
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
         CALL UPCASE(WORD)
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102), II
               END IF
 1102       CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               IF (.NOT.ALLOPT)
     &            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/4A/)') ' Keyword ',WORD,
     *         ' not recognized for ',WORD1
         ELSE IF (PROMPT .EQ. '#' .OR. PROMPT .EQ. '!') THEN
            GO TO 100
         ELSE IF (PROMPT .EQ. '*') THEN
            RETURN
         ELSE
            WRITE (LUPRI,'(/3A/2A/)')
     *         ' Keyword ',WORD,' does not begin with',
     *         ' one of the four characters ".*!#" for ',WORD1
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal keyword for "*AUXILIARY INPUT"')
C
C *** Option 1 >ONESUP<  use same units for SUPM. and ONEL. integrals
C              (obsolete, removed March 2000 /Kenneth Ruud)
  101 GO TO 100
C
C *** Option 2 >NOSUPMAT<  Do not use SUPMAT
  102 LUSUPM = -1
      GO TO 100
C
      END
C  /* Deck errinp */
      SUBROUTINE ERRINP(ITYP,INPERR)
C
C     Last revision 31-Jan-1984 hjaaj
C     (based on CASSCF RELEASE 79 01 01)
C
C     ***** ERROR MESSAGES IN INPUT SECTION *****
C
C
#include "implicit.h"
C
      CHARACTER TAB1*15,TAB2*30
      PARAMETER (TAB1=' ERROR MESSAGE:')
      PARAMETER (TAB2=' IS NOT WITHIN GIVEN LIMITS OF')
C
C
C Used from common blocks:
C   INFVAR : NWOPT
C
#include "maxash.h"
#include "maxorb.h"
#include "priunit.h"
#include "infinp.h"
#include "infvar.h"
#include "inforb.h"
#include "infdim.h"
#include "inftap.h"
#include "infpri.h"
C
C
      INPERR = INPERR + 1
      GO TO ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     *       11,12,13,14,15,16,17),ITYP
C
    1 CONTINUE
      RETURN
C
    2 WRITE(LUPRI,200) TAB1,NACTEL,TAB2
  200 FORMAT(/A/' NUMBER OF ACTIVE ELECTRONS=',I3,A,' OF 0 TO 35')
      RETURN
C
    3 S = ISPIN-1
      S = S / 2
      WRITE(LUPRI,300) TAB1,S,TAB2
  300 FORMAT(/A/' SPIN QUANTUM NUMBER=',F5.1,A,' 0 TO 10')
      RETURN
C
    4 WRITE(LUPRI,400) TAB1,NSYM,TAB2
  400 FORMAT(/A/' NUMBER OF SYMMETRIES=',I3,A,' 1 TO 8')
      RETURN
C
    5 WRITE(LUPRI,500) TAB1,MAXMAC,TAB2
  500 FORMAT(/A/' MAXIMUM NUMBER OF ITERATIONS=',I3,A,' 1 TO 50')
      RETURN
C
    6 CONTINUE
      RETURN
C
    7 WRITE(LUPRI,700) TAB1,THRGRD,TAB2
  700 FORMAT(/A/' THRESHOLD FOR GRADIENT CONVERGENCE=',1P,E10.2,A,
     *       ' 0.0 TO 1.0')
      RETURN
C
    8 WRITE(LUPRI,800) TAB1,NORBT,TAB2,MXCORB
  800 FORMAT(/A/' TOTAL NUMBER OF ORBITALS=',I4,A,' 0 TO',I6)
      RETURN
C
    9 WRITE(LUPRI,900) TAB1,NRHFT,TAB2,MAXOCC
  900 FORMAT(/A/' NUMBER OF RHF ORBITALS=',I4,A,' 0 TO',I6)
      RETURN
C
   10 WRITE(LUPRI,1000) TAB1,NASHT,TAB2,MAXASH
 1000 FORMAT(/A/' NUMBER OF ACTIVE ORBITALS=',I4,A,' 0 TO',I6)
      RETURN
C
   11 WRITE(LUPRI,1100) TAB1,NSSHT,TAB2,MXCORB
 1100 FORMAT(/A/' NUMBER OF SECONDARY ORBITALS=',I4,A,' 0 TO',I6)
      RETURN
C
   12 WRITE(LUPRI,1200) TAB1,(IMOORD(I),I=1,NORBT)
 1200 FORMAT(/A/' REORDERING INDICES DO NOT GIVE CORRECT SUM CHECK',
     *       /(23X,10I4))
      RETURN
C
   13 WRITE(LUPRI,1300) TAB1,NROOTS,TAB2,MAXRTS
 1300 FORMAT(/A/' NUMBER OF NEO ROOTS=',I3,A,' 1 TO',I4)
      RETURN
C
   14 WRITE(LUPRI,1400) TAB1,(IROOT(I),I=1,NROOTS)
 1400 FORMAT(/A/' NEO ROOT INDICES DUPLICATE OR NEGATIVE'/(23X,10I3))
      RETURN
C
   15 RETURN
C
   16 WRITE(LUPRI,1600) TAB1,NWOPT,TAB2,MAXWOP
 1600 FORMAT(/A/' NUMBER OF ORBITAL ROTATIONS =',I8,A,' 0 TO',I8)
      RETURN
C
   17 RETURN
C
      END
C  /* Deck SIR_prtinp */
      SUBROUTINE SIR_PRTINP(IWUNIT,NUMRUN,IRDMO)
C
C Last revisions
C   12-Oct-1984 / 2-Jul-1984 / 4-May-1984 hjaaj
C    7-Jan-1985 hjaaj (write date and time to output if VAX)
C   DFT modifications tuh
C
C  ***** PROGRAM SIRIUS (INPUT SECTION): PRINT ON UNIT IWUNIT *****
C
C  OBJECTIVE: Print input information on unit IWUNIT
C             using at most 72 characters per row
C
C  CALLED FROM SIRRDI
C
      use pelib_interface, only: use_pelib
#include "implicit.h"
C
C Used from common blocks:
C   INFINP : DOMC,DOCI,DIRFCK,FLAG(21),MCTYPE,THR*,NFIELD,...,?
C   INFVAR : NCONF,NWOPT,NVAR
C   INFIND : ISMO(),?
C   INFDIM : MAXPHP
C   SCBRHF : IOPRHF,NFRRHF(*),MAXFCK
C   INFMP2 : NFRMP2(*)
C   DFTCOM : DODFT, DOSRBCK, ...
C   DFTACB : DFTIPT, DFTBR1, DFTBR2
C   gnrinf : NEWSYM,SIR_INPPRC,EMBEDDING,QM3,QMMM,?
C   infpar.h : PARHER
C   pgroup : GROUP, REP
C   pcmlog : PCM
C
#include "priunit.h"
#include "mxcent.h"
#include "maxash.h"
#include "maxorb.h"
C
#include "gnrinf.h"
#include "infinp.h"
#include "infvar.h"
#include "infind.h"
#include "inforb.h"
#include "infdim.h"
#include "scbrhf.h"
#include "infmp2.h"
#include "inftap.h"
#include "infpri.h"
#include "huckel.h"
#include "dftcom.h"
#include "dfterg.h"
#include "dftacb.h"
#include "infpar.h"
#include "pgroup.h"
#include "pcmlog.h"
#include "spinfo.h"
C
      CHARACTER*6 RDMOLB(12)
      DATA RDMOLB/'--1-- ','--2-- ','FORM18','FORM12','--5-- ',
     &            '--6-- ','OLDORB','H1DIAG','NEWORB','SIRIFC',
     &            'EHT   ','EWMO  '/
      CHARACTER*33 TXTOCC(0:6)
      DATA TXTOCC/'NOTHING because of ERROR(s).     ',
     &            'user input.                      ',
     &            'SIRIUS restart file (SIRIUS.RST).',
     &            'SIRIUS interface file (SIRIFC).  ',
     &            'diagonal of H1 matrix.           ',
     &            'extended Huckel guess.           ',
     &            'total charge (only one symmetry!)'/
C
      CHARACTER*70 STARS,RELEAS
      DATA STARS(1:36)  /'************************************'/
      DATA STARS(37:70) /'**********************************'/
      DATA RELEAS(1:36) /'*SIRIUS* a direct, restricted step, '/
      DATA RELEAS(37:70)/'second order MCSCF program       *'/
C
#include "orbtypdef.h"
C
C     ***** WRITE PROGRAM IDENTIFICATION *****
C
      CALL QENTER('SIR_PRTINP')
      WRITE(IWUNIT,200) STARS,RELEAS,STARS
 200  FORMAT(3(/1X,A))

      IF (SIR_INPPRC) THEN
         IF (NEWSYM) THEN
            WRITE(IWUNIT,'(/A)')
     &   ' (Wave function specification NOT suppressed here,'//
     &   ' because symmetry has changed.)'
         ELSE
            WRITE(IWUNIT,'(/A)')
     &   ' (Wave function specification suppressed here,'//
     &   ' see output above for initial geometry.)'
            GO TO 8000
         END IF
      END IF
#ifdef VAR_AUTHORS
         WRITE(IWUNIT,400)
  400    FORMAT(/' Main authors of the SIRIUS module:'
     *   /T5,'Hans Joergen Aa. Jensen, Univ. Southern Denmark, Denmark',
     *       ' (principal author)',
     *   /T5,'Hans Agren,              Royal Inst. Technology, Sweden')
         WRITE(IWUNIT,402)
  402    FORMAT(/' Direct CI and density matrix routines:'
     *   /T5,'Jeppe Olsen,             University of Aarhus,   Denmark')
         CALL TRAUTH(IWUNIT)
c        WRITE(IWUNIT,410)
c410     FORMAT(/,' Questions and correspondence should be directed',
c    *        ' to H.J.Aa. Jensen',/,T5,'E-mail: hjj@sdu.dk',
c    *        //,T2,'Standard reference:',
c    *        /,T5,'H.J.Aa. Jensen, H. Agren, and J. Olsen,',
c    *        ' Chapter 8 in MOTECC-90')
         IF (FLAG(16)) CALL SOLREF(IWUNIT)
#endif
C
C     Stamp date, time and hostname to IWUNIT
C
      CALL TSTAMP(' ',IWUNIT)
C
C  ***** Print title *****
C
C
      IF (NTIT .GT. 0) THEN
         WRITE(IWUNIT,500) TITLE(1)
         DO 5 ITIT=2,NTIT
            WRITE(IWUNIT,600) TITLE(ITIT)
    5    CONTINUE
      END IF
  500 FORMAT(/' Title     :',T14,A60)
  600 FORMAT(T14,A60)
      WRITE (IWUNIT,700) (TITMOL(J), J=1,2)
  700 FORMAT(/' Title lines from ".mol" input file:'/T6,A72/T6,A72)
C
      WRITE(IWUNIT,1700) LUPRI,IPRI6,LUW4,IPRI4
 1700 FORMAT(/' Print level on unit LUPRI =',I4,' is',I4,
     *       /' Print level on unit LUW4  =',I4,' is',I4)
C
C Type of calculation:
C
 1050 FORMAT(/'@',T6,A)
 1051 FORMAT('@',T16,A)
 1052 FORMAT(/'@',T6,A,$)
 1053 FORMAT(A)
 1054 FORMAT('@',T16,A,$)
      IF (IORTO .EQ. 1) THEN
         WRITE (IWUNIT,1050) 'Orthonormalization of molecular orbitals.'
         WRITE (IWUNIT,1100) (I,       I=1,NSYM)
         WRITE (IWUNIT,1110) ('  ---', I=1,NSYM+1)
         WRITE (IWUNIT,1450) (NORBT,NORB(I),I=1,NSYM)
         WRITE (IWUNIT,1500) (NBAST,NBAS(I),I=1,NSYM)
         GO TO 8000
      END IF
      WFTYPE = ' '
      NUMRUN_LOCAL = NUMRUN
      IF (DOLUCITA) THEN
         IF (WFTYPE .EQ. ' ') WFTYPE='LUCITA'
         WRITE (IWUNIT,1050) 'Generalized Active Space CI (with LUCITA)'
         NUMRUN_LOCAL = NUMRUN_LOCAL - 1
         IF (NUMRUN_LOCAL .GT. 0) WRITE (IWUNIT,1051) 'preceded by :'
      END IF
      IF ( DOMC .OR. DOCI ) THEN
         IF ( DONEVPT ) then
            IF (WFTYPE .EQ. ' ') WFTYPE='NEVPT2'
C Oct 2003 hjaaj: from koopro4.F:
C     This program calculates the second order correction to the energy
C     using the n-electron valence state perturbation theory (NEVPT)
C     in the partially contracted (PC) ad strongly contracted (SC)
C     variants using Dyall's Hamiltonian for the definition of the
C     zero order energies. All the equations are based on the spin free
C     formulation of NEVPT.
           WRITE (IWUNIT,1050)
     &     "N-electron valence state 2. order perturbation theory "//
     &     "using Ken Dyall's zeroth-order Hamiltonian"
           NUMRUN_LOCAL = NUMRUN_LOCAL - 1
           IF (NUMRUN_LOCAL .GT. 0) WRITE (IWUNIT,1051) 'preceded by :'
         ENDIF
         IF (DOAPSG) THEN
            IF (WFTYPE .EQ. ' ') WFTYPE='APSG'
            WRITE (IWUNIT,1050)
     &      'Antisymmetrized Product of Strongly Orthogonal Geminals'
            NUMRUN_LOCAL = NUMRUN_LOCAL - 1
         ELSE IF ( DOMCSRDFT ) THEN
            IF (WFTYPE .EQ. ' ') WFTYPE='MC-srDFT'
            WRITE (IWUNIT,1050) 'MC-srDFT optimization.'
            NUMRUN_LOCAL = NUMRUN_LOCAL - 1
            IF ( DORSP ) WRITE (IWUNIT,1050)
     *         'MC-srDFT response calculation.'
         ELSE IF ( DOMC ) THEN
            IF (WFTYPE .EQ. ' ') WFTYPE='MC-SCF'
            WRITE (IWUNIT,1050) 'MCSCF optimization.'
            NUMRUN_LOCAL = NUMRUN_LOCAL - 1
            IF ( DORSP ) WRITE (IWUNIT,1050)
     *         'MCSCF response calculation.'
         ELSE
            IF (DOCINO) THEN
               WRITE (IWUNIT,1050) 'CI-NO calculation.'
            ELSE
               WRITE (IWUNIT,1050)
     *            'Configuration interaction calculation.'
               IF ( DORSP ) WRITE (IWUNIT,1050)
     *            'CI response calculation.'
            END IF
            IF (DOFCI) THEN
               IF (WFTYPE .EQ. ' ') WFTYPE='Full CI'
               WRITE (IWUNIT,1050) 'CI type: Full CI.'
            ELSE IF (DOCISD) THEN
               IF (WFTYPE .EQ. ' ') WFTYPE='CISD'
               WRITE (IWUNIT,1050) 'CI type: CISD.'
            ELSE
               IF (WFTYPE .EQ. ' ') WFTYPE='CI'
               WRITE (IWUNIT,1050)
     &           'CI type: user specified active spaces.'
            END IF
            NUMRUN_LOCAL = NUMRUN_LOCAL - 1
         END IF
         IF (MCTYPE .EQ. 2) THEN
            IF (CORRLX) THEN
               WRITE(IWUNIT,1051)
     *         'Type: Relaxed single core hole calculation using RAS'
            ELSE IF (CORHOL) THEN
               WRITE(IWUNIT,1051)
     *         'Type: Frozen single core hole calculation using RAS'
               WRITE(IWUNIT,'(T16,A,I4,A,I2,A)')
     &         '      (The core hole orbital, no.',JCHORB,
     &         ', in symmetry',JCHSYM,' has been moved to RAS1)'
            ELSE
               WRITE(IWUNIT,1051)
     *         'Type: restricted active space calculation (RAS)'
            END IF
         ELSE IF (MCTYPE .EQ. 1) THEN
            WRITE(IWUNIT,1051)
     *         'Type: complete active space calculation (CAS).'
         ELSE
C          hjaaj Aug 2005 - program validity check, should not be possible
C          to get here unless someone modifies SIRIUS incorrectly
           CALL QUIT('SIR_PRTINP ERROR - MCTYPE invalid for DOMC/DOCI')
         END IF
         IF (FLAG(23)) THEN
            WRITE(IWUNIT,1051)
     &      '      Selected active-active rotations included by user.'
         END IF
      ELSE IF (DOMP2 .AND. DODFT .AND. WDFTMP.NE.0.0D0) THEN
         IF (MP2_SCS) THEN
            IF (WFTYPE .EQ. ' ') WFTYPE='DH-SCS-MP2-DFT'
            WRITE(IWUNIT,1050) 'DH-SCS-MP2-DFT, '//
     &         'double hybrid DFT with SCS-MP2.'
         ELSE IF (MP2_SOS) THEN
            IF (WFTYPE .EQ. ' ') WFTYPE='DH-SOS-MP2-DFT'
            WRITE(IWUNIT,1050) 'DH-SOS-MP2-DFT, '//
     &         'double hybrid DFT with SOS-MP2.'
         ELSE
            IF (WFTYPE .EQ. ' ') WFTYPE='DH-MP2-DFT'
            WRITE(IWUNIT,1050) 'DH-MP2-DFT, '//
     &         'double hybrid DFT with MP2.'
         END IF
         IF (IOPRHF .NE. 0 .OR. HSROHF) WRITE (IWUNIT,1051)
     &      '      (ROHF open shell orbitals are not included.)'
         NUMRUN_LOCAL = NUMRUN_LOCAL - 1
      ELSE IF ( DOMP2 ) THEN
         IF (MP2_SCS) THEN
            IF (WFTYPE .EQ. ' ') WFTYPE='SCS-MP2'
            WRITE(IWUNIT,1050) 'SCS-MP2'
         ELSE IF (MP2_SOS) THEN
            IF (WFTYPE .EQ. ' ') WFTYPE='SOS-MP2'
            WRITE(IWUNIT,1050) 'SOS-MP2'
         ELSE
            IF (WFTYPE .EQ. ' ') WFTYPE='MP2'
            WRITE(IWUNIT,1050) 'MP2, '//
     &         'closed-shell Moeller-Plesset second-order calculation.'
         END IF
         IF (IOPRHF .NE. 0 .OR. HSROHF) WRITE (IWUNIT,1051)
     &      '      (ROHF open shell orbitals are not included.)'
         NUMRUN_LOCAL = NUMRUN_LOCAL - 1
      ELSE IF ( DOCCSD) THEN
         IF (WFTYPE .EQ. ' ') WFTYPE='CC'
         WRITE (IWUNIT,1050) '(Integral direct) CC calculation.'
         NUMRUN_LOCAL = NUMRUN_LOCAL - 1
      ELSE IF ( DOSCF ) THEN
         IF (HSROHF) THEN
            WRITE(IWUNIT,1052) 'Restricted, high spin open shell'
         ELSE IF (IOPRHF .EQ. 0) THEN
            WRITE(IWUNIT,1052) 'Restricted, closed shell'
         ELSE
            WRITE(IWUNIT,1052) 'Restricted, one open shell'
         END IF
         IF (DOHFSRDFT) THEN
            IF (WFTYPE .EQ. ' ') WFTYPE='HF-srDFT'
            WRITE (IWUNIT,1053) ' HF-srDFT hybrid calculation.'
            IF ( DORSP ) WRITE (IWUNIT,1050)
     &      'Time-dependent HF-srDFT calculation (response properties).'
         ELSE IF (DODFT) THEN
             IF (WFTYPE .EQ. ' ') WFTYPE='KS-DFT'
             WRITE (IWUNIT,1053) ' Kohn-Sham DFT calculation.'
             IF ( DORSP ) WRITE (IWUNIT,1050)
     &           'Time-dependent Kohn-Sham DFT calculation (TD-DFT).'
         ELSE
             IF (WFTYPE .EQ. ' ') WFTYPE='HF'
             WRITE (IWUNIT,1053) ' Hartree-Fock calculation.'
             IF ( DORSP ) WRITE (IWUNIT,1050)
     &           'Time-dependent Hartree-Fock calculation '//
     &           '(random phase approximation).'
         END IF
         NUMRUN_LOCAL = NUMRUN_LOCAL - 1
      ELSE
         WRITE (IWUNIT,1050) 'No wave function optimization.'
         IF ( DORSP ) WRITE (IWUNIT,1050)
     &      'Response calculation.'
      END IF
      IF (NUMRUN_LOCAL .GT. 0) THEN
         WRITE(IWUNIT,1050)
     &      'This is a combination run starting with'
         IF (DOSCF) THEN
            IF (HSROHF) THEN
               WRITE(IWUNIT,1054) 'a restricted, high spin open shell'
            ELSE IF (IOPRHF .EQ. 0) THEN
               WRITE(IWUNIT,1054) 'a restricted, closed shell'
            ELSE
               WRITE(IWUNIT,1054) 'a restricted, one open shell'
            END IF
            IF (DOHFSRDFT) THEN
               WRITE (IWUNIT,1053) ' HF-srDFT hybrid calculation'
            ELSE IF (DODFT) THEN
               WRITE (IWUNIT,1053) ' Kohn-Sham DFT calculation'
            ELSE
               WRITE (IWUNIT,1053) ' Hartree-Fock calculation'
            END IF
            NUMRUN_LOCAL = NUMRUN_LOCAL - 1
         END IF
         IF (DOMP2 .AND. NUMRUN_LOCAL .GT. 0) THEN
            IF (MP2_SCS) THEN
               WRITE(IWUNIT,1051) 'an SCS-MP2 calculation'
            ELSE IF (MP2_SOS) THEN
               WRITE(IWUNIT,1051) 'an SOS-MP2 calculation'
            ELSE
               WRITE(IWUNIT,1051) 'an MP2 calculation'
            END IF
            IF (MP2_TDA) WRITE(IWUNIT,1051) '   (saving for TDA-SOPPA)'
            NUMRUN_LOCAL = NUMRUN_LOCAL - 1
         END IF
         IF (FCVORB) THEN
            IF (NMVOT .EQ. 0) THEN
               WRITE(IWUNIT,1051)'diagonal h1 virtual orbitals'
            ELSE
               WRITE(IWUNIT,1051)'diagonal FC(MVO) virtual orbitals'
            END IF
         END IF
         IF (DOCI .AND. NUMRUN_LOCAL .GT. 0) THEN
            IF (DOCINO) THEN
               WRITE(IWUNIT,1051) 'a CI-NO calculation'
            ELSE
               WRITE(IWUNIT,1051) 'a CI calculation'
            END IF
            NUMRUN_LOCAL = NUMRUN_LOCAL - 1
         END IF
         WRITE (IWUNIT,'()')
      END IF
      IF ( DIRFCK ) THEN
      IF ( PARHER ) THEN
         WRITE (IWUNIT,'(A)')
     &      ' Fock matrices are calculated directly and in parallel'//
     &      ' without use of integrals on disk.'
      ELSE
         WRITE (IWUNIT,'(A)')
     &      ' Fock matrices are calculated directly'//
     &      ' without use of integrals on disk.'
      END IF
      END IF
C
C ***** Hamiltonian specifications
C
      IF (NFIELD .GT. 0) THEN
         WRITE (IWUNIT,'(/A/A/A)')
     &      '@  The molecule is placed in a static field.',
     &      '@  Field strength (a.u.)          Field operator',
     &      '@  ---------------------------------------------'
         DO IFIELD = 1,NFIELD
            WRITE (IWUNIT,'(A,1P,G18.8,18X,A8)')
     &      '@  ',EFIELD(IFIELD),LFIELD(IFIELD)
         END DO
         WRITE (IWUNIT,'(A)')
     &      '@  ---------------------------------------------'
      END IF
      IF (IRDMO .GT. 0 .AND. IRDMO .LE. 10) THEN
C
         JRDMO = IRDMO
         IF (IRDMO.EQ.8 .AND. DOHUCKEL) THEN
C        ... if HUCKEL use JRDMO = 11 or 12 to get correct output
            IF (EWMO) THEN
               JRDMO = 12
            ELSE
               JRDMO = 11
            END IF
         END IF
         WRITE (IWUNIT,'(/A/3A)')
     &   ' Initial molecular orbitals are obtained according to',
     &   ' ".MOSTART ',RDMOLB(JRDMO),'" input option'
      END IF
      IF (LMOORD) THEN
         WRITE (IWUNIT,'(A)') ' and reordered as requested '//
     &   'with .REORDER or indirectly with .CORHOL'
      END IF
C
      IF (EMBEDDING) THEN
         WRITE (IWUNIT,'(/A)')
     &   '@    QM part is embedded in an environment :'
       IF (FLAG(16)) THEN
        IF (RSOL(1) .EQ. RSOL(2) .AND. RSOL(1) .EQ. RSOL(3)) THEN
         WRITE (IWUNIT,'(/A,9X,A,2(/A,15X,A,F12.6),/A,15X,A,I12)')
     &      '@','SOLVATION MODEL: molecule is in a spherical'//
     &        ' cavity in a dielectric medium,',
     &      '@',' dielectric constant =',EPSOL,
     &      '@','cavity radius (bohr) =',RSOL(1),
     &      '@','         max l value =',LSOLMX
        ELSE
         WRITE (IWUNIT,'(/A,9X,A,2(/A,15X,A,3F12.6),/A,15X,A,I12)')
     &      '@','SOLVATION MODEL: molecule is in a ellipsoidal'//
     &        ' cavity in a dielectric medium,',
     &      '@','     dielectric constant =',EPSOL,
     &      '@','cavity dimensions (bohr) =',RSOL(1),RSOL(2),RSOL(3),
     &      '@','             max l value =',LSOLMX
        END IF
        IF (INERSI) WRITE (IWUNIT,'(2(/A,9X,A))')
     &      '@','This is the initial state in a calculation with'//
     &      ' inertial polarization.',
     &      ' ','Converged T(lm) values are saved on the Sirius'//
     &      ' interface file.'
        IF (INERSF) WRITE (IWUNIT,'(/A,9X,A/A,14X,A,F12.6)')
     &      '@','This is the final state in a calculation with'//
     &      ' inertial polarization.',
     &      '@','static dielectric constant =',EPSTAT
       END IF
       IF (PCM)  WRITE (IWUNIT,'(/A,9X,A)') '@','Solvation model: PCM'
       IF (QM3)  WRITE (IWUNIT,'(/A,9X,A)') '@','Model: QM3'
       IF (QMMM) WRITE (IWUNIT,'(/A,9X,A)') '@','Model: QMMM'
       IF (USE_PELIB()) WRITE (IWUNIT,'(/A,9X,A)') '@','Model: PE'
      END IF
C
C ***** Wave function specifications
C
      WRITE(IWUNIT,1000)
      IF (DOCCSD) THEN
         WRITE(IWUNIT,'(/5X,A/)')
     &   'For the specification of the Coupled Cluster: see below.'
      END IF
      IF (DOLUCITA) THEN
         WRITE(IWUNIT,'(/5X,A/)')
     &   'For the specification of the GAS CI in LUCITA: see below.'
      END IF
      IF (DOAPSG) THEN
         WRITE(IWUNIT,'(/5X,A/)')
     &   'For the specification of the APSG: see below.'
      END IF
      IF (USE_PELIB()) THEN
         WRITE(IWUNIT,'(/5X,A/5X,A/)')
     &   'Wave function environment is modeled using the',
     &   'polarizable embedding model (PE library).'
        WFTYPE = 'PE-'//WFTYPE
      END IF
C     for SIRIUS (not LUCITA) :
      NINEL = 2*NISHT
      CALL MOLCHR(ICHARG)
C     ... get total charge
      IF (.NOT. (DOLUCITA .AND. NUMRUN.EQ.1) ) WRITE(IWUNIT,1002)
     &   TRIM(WFTYPE),NINEL,NACTEL,ICHARG,ISPIN,MS2,
     &   NSYM,GROUP,LSYM,REP(LSYM-1)
 1000 FORMAT(/5X,'Wave function specification'
     &       /5X,'============================')
 1002 FORMAT( '@    Wave function type        --- ',A,' ---',
     &       /'@    Number of closed shell electrons',I12,
     &       /'@    Number of electrons in active shells',I8,
     &       /'@    Total charge of the molecule',I16,
     &      //'@    Spin multiplicity and 2 M_S',I17,I10,
     &       /'@    Total number of symmetries',I18,
     &        ' (point group: ',A,')',
     &       /'@    Reference state symmetry',I20,
     &        ' (irrep name : ',A,')')
C
      IF (MCTYPE .EQ. 2) THEN
         WRITE (IWUNIT,1075) NAS1T, NAS2T, NAS3T
         IF (NAS1T .GT. 0) THEN
            NHLMN1 = 2*NAS1T - NELMX1
            NHLMX1 = 2*NAS1T - NELMN1
            WRITE (IWUNIT,1071) NELMN1,NHLMX1, NELMX1,NHLMN1
         END IF
         IF (NAS3T .GT. 0) WRITE (IWUNIT,1073) NELMN3, NELMX3
      END IF
 1071 FORMAT(/'@    Minimum number of electrons in RAS1 :',I4,
     &        ' (maximum',I4,' holes)',
     &       /'@    Maximum number of electrons in RAS1 :',I4,
     &        ' (minimum',I4,' holes)')
 1073 FORMAT(/'@    Minimum number of electrons in RAS3 :',I4,
     &       /'@    Maximum number of electrons in RAS3 :',I4)
 1075 FORMAT(/'Number of active orbitals in RAS1, RAS2, and RAS3 :',
     &       3I4)
C
C     DFT run (standard Kohn-Sham -- not short range)
C
      IF (DODFT) THEN
         CALL DFTREPORT
      END IF

      IF (DOHFSRDFT .OR. DOCISRDFT .OR. DOMCSRDFT) THEN
C        ... Check for inconsistencies :

         WRITE (IWUNIT,'(/A,2(/A,6X,2A))')
     &        '@    This is a lrWFT-srDFT calculation using the ',
     &        '@',SRXFUN,' short range exchange functional',
     &        '@',SRCFUN,' short range correlation functional'
         IF (SRLOCALSPIN .NE. ' ') WRITE (IWUNIT,'(A,6X,2A)')
     &        '@',SRLOCALSPIN, ' local spin model'
         IF (DOSRC_PBELO) WRITE(IWUNIT,'(A,1P,D10.2)')
     &      'HEAVISIDE_PVALUE = ', HEAVISIDE_PVALUE

         IF (IWINT.NE.0) WRITE(IWUNIT,'(/5X,A,F10.6)')
     &        'Weighted short range interpolation, COPFAC = ',COPFAC
         WRITE (IWUNIT,'(/A,2F12.6)')
     &        '@    sr-DFT and exact sr-HF exchange weights:',
     &         (1.0D0-HFXFAC),HFXFAC

         IF (DOSRBCK.AND.DOMCSRDFT)
     &      CALL QUIT(' SRBCK functional not implemented for MCSRDFT !')
!hjaaj Oct 09: allow chival=0 for tests
!        IF (DOMCSRDFT .AND. CHIVAL.LE.0.D0)  THEN
!          CALL QUIT(' MCSRDFT not allowed for CHIVAL=0')
!        ENDIF
         IF (DODFT) THEN
           CALL QUIT('DFT and srDFT not both allowed in same run')
         END IF

      ENDIF
C
C ***** Orbital specifications *****
C
      WRITE(IWUNIT,1100) (I,       I=1,NSYM)
      WRITE(IWUNIT,1102) (REP(I-1),I=1,NSYM)
      WRITE(IWUNIT,1110) ('  ---', I=1,NSYM+1)
      IF (MCTYPE .GT. 0) THEN
         WRITE(IWUNIT,1200) NISHT,(NISH(I),I=1,NSYM)
         IF (NASHT .GT. 0) THEN
            WRITE(IWUNIT,1300) NASHT,(NASH(I),I=1,NSYM)
            IF (MCTYPE .GT. 1) THEN
               IF (NAS1T .GT. 0)
     &            WRITE(IWUNIT,1320) NAS1T,(NAS1(I),I=1,NSYM)
               IF (NAS2T .GT. 0)
     &            WRITE(IWUNIT,1330) NAS2T,(NAS2(I),I=1,NSYM)
               IF (NAS3T .GT. 0)
     &            WRITE(IWUNIT,1340) NAS3T,(NAS3(I),I=1,NSYM)
            END IF
         END IF
         WRITE(IWUNIT,1400) NSSHT,(NSSH(I),I=1,NSYM)
      ELSE IF (.NOT.AUTOCC) THEN
         WRITE(IWUNIT,1550) NRHFT,(NRHF(I),I=1,NSYM)
         IF (IOPRHF .NE. 0 .OR. HSROHF)
     *      WRITE (IWUNIT,1555) NASHT,(NASH(I),I=1,NSYM)
         WRITE(IWUNIT,1400) NSSHT,(NSSH(I),I=1,NSYM)
      END IF
      WRITE(IWUNIT,1450) NORBT,(NORB(I),I=1,NSYM)
      WRITE(IWUNIT,1500) NBAST,(NBAS(I),I=1,NSYM)
      IF (DOSCF .AND. (AUTOCC .OR. MCTYPE .GT. 0)) THEN
         IF (AUTOCC) THEN
           IF (DODFT) THEN
             WRITE(IWUNIT,'(/5X,A)')
     &       ' ** Automatic occupation of RKS orbitals **'
           ELSE IF (DOHFSRDFT) THEN
             WRITE(IWUNIT,'(/5X,A)')
     &       ' ** Automatic occupation of RHF-srDFT orbitals **'
           ELSE
             WRITE(IWUNIT,'(/5X,A)')
     &       ' ** Automatic occupation of RHF orbitals **'
           END IF
         END IF
         WRITE(IWUNIT,'(/5X,2A/5X,A)')
     &       ' -- Initial occupation of symmetries is determined from ',
     &       TXTOCC(INIOCC),
     &       ' -- Initial occupation of symmetries is :'
         WRITE(IWUNIT,1550) NRHFT,(NRHF(I),I=1,NSYM)
         IF (HSROHF) THEN
            NROHFT = ISUM(8,NROHF,1)
            WRITE (IWUNIT,1555) NROHFT,(NROHF(I),I=1,NSYM)
         END IF
         IF (IOPRHF .NE. 0) WRITE (IWUNIT,1560) IOPRHF
      END IF
      IF (SUPSYM) CALL AVEPRT(IWUNIT)
      IF (BOYORB .AND. (DOSCF .OR. DOMC)) WRITE(IWUNIT,'(/5X,A)')
     &   'Boys localization of all doubly occupied orbitals'//
     %   ' will be performed after RHF/DFT/MCSCF'
 1100 FORMAT(/5X,'Orbital specifications',
     *       /5X,'======================',
     *       /'@    Abelian symmetry species',T38,'  All |',8I5)
 1102 FORMAT('@',T44,'|',8(2X,A3))
 1110 FORMAT(T38,A5,' |',8A5)
 1200 FORMAT('@    Inactive orbitals',T38,I5,' |',8I5)
 1300 FORMAT('@    Active orbitals',T38,I5,' |',8I5)
 1320 FORMAT('@    - Active orbitals in RAS1',T38,I5,' |',8I5)
 1330 FORMAT('@    - Active orbitals in RAS2',T38,I5,' |',8I5)
 1340 FORMAT('@    - Active orbitals in RAS3',T38,I5,' |',8I5)
 1400 FORMAT('@    Secondary orbitals',T38,I5,' |',8I5)
 1450 FORMAT('@    Total number of orbitals',T38,I5,' |',8I5)
 1500 FORMAT('@    Number of basis functions',T38,I5,' |',8I5)
 1550 FORMAT('@    Occupied SCF orbitals',T38,I5,' |',8I5)
 1555 FORMAT('@    Open shell SCF orbitals',T38,I5,' |',8I5)
 1560 FORMAT(/'@    plus one open shell SCF orbital of symmetry',I4)
C
C ***** Optimization information *****
C       FLAG(27) : use determinants
C
      IF (MCTYPE .GT. 0 .OR. (FLAG(21) .AND. .NOT. AUTOCC)) THEN
         IF (DOMC .OR. MCTYPE .LE. 0) THEN
            NWOPT_print = NWOPT
         ELSE ! no MCSCF, only CI
            NWOPT_print = 0
         END IF
         NVAR_print = NCONF + NWOPT_print
         IF (DOAPSG) THEN
            WRITE (IWUNIT,5020) NCONF,NWOPT_print,NVAR_print
         ELSE IF (FLAG(27)) THEN
            WRITE (IWUNIT,5010) NCONF,NWOPT_print,NVAR_print
         ELSE
            WRITE (IWUNIT,5000) NCONF,NWOPT_print,NVAR_print
         END IF
      END IF
 5000 FORMAT (/,5X,'Optimization information',
     *        /,5X,'========================',
     *        /'@    Number of configurations   ',I15,
     *        /'@    Number of orbital rotations',I15,
     *        /,5X, '------------------------------------------',
     *        /'@    Total number of variables  ',I15)
 5010 FORMAT (/,5X, 'Optimization information',
     *        /,5X, '========================',
     *        /'@    Number of determinants     ',I15,
     *        /'@    Number of orbital rotations',I15,
     *        /,5X, '------------------------------------------',
     *        /'@    Total number of variables  ',I15)
 5020 FORMAT (/,5X, 'Optimization information',
     *        /,5X, '========================',
     *        /'@    Number of geminal coeffients',I14,
     *        /'@    Number of orbital rotations',I15,
     *        /,5X, '------------------------------------------',
     *        /'@    Total number of variables  ',I15)
C
C
C
      IF (DOMC) THEN
         WRITE(IWUNIT,1600) MAXMAC,MAXMIC,THRMC,
     &                      LROOTS,MAXCIT,NROOTS
         IF (ISTATE .GT. 1) THEN
            IF (NMCAVER .EQ. 1) THEN
               WRITE(IWUNIT,'(/A,I4,A)')
     &         '@    This calculation converges to state no.',ISTATE,
     &         ' for the specified symmetry and spin species.'
            ELSE
               WRITE(IWUNIT,'(/A,I4,A//A)')
     &         '@    State-average MCSCF of',NMCAVER,' states'//
     &         ' for the specified symmetry and spin species.',
     &         '@      State no.        Weight'
               WRITE(IWUNIT,'(A,I12,F12.5)')
     &         ('@',IROOT(I),WEIGHT_MCAVER(I), I=1,NMCAVER)
            END IF
         ELSE
            WRITE(IWUNIT,'(/A)')
     &         '@    This calculation converges to the lowest state'//
     &         ' for the specified symmetry and spin species.'
         END IF
         WRITE(IWUNIT,2500) MAXJT
         IF (JOLSEN) WRITE(IWUNIT,2550)
         IF (MAXPHP.GT.1) THEN
            WRITE(IWUNIT,2560) MAXPHP
            IF (RESPHP) WRITE(IWUNIT,2570)
         END IF
      ELSE IF (DOCI) THEN
         WRITE(IWUNIT,1800) MXCIMA, THRCI, NROOCI
         IF (JOLSEN) WRITE(IWUNIT,2550)
         IF (MAXPHP.GT.1) THEN
            WRITE(IWUNIT,2560) MAXPHP
            IF (RESPHP) WRITE(IWUNIT,2570)
         END IF
         IF (DORSP) THEN
            WRITE(IWUNIT,2600) ISTACI
         ELSE IF (ISTACI .GT. 0) THEN
            WRITE(IWUNIT,'(/5X,A,I4)')
     *      'The CI calculation converges only root no.',ISTACI
         ELSE IF (ISTACI .LT. 0 .AND. DOCISRDFT) THEN
            WRITE(IWUNIT,'(/5X,A,I4)')
     *      'The CI-DFT calculation gets DV(ref) from root no.',-ISTACI
         END IF
      ELSE IF (DOSCF) THEN
         WRITE(IWUNIT,1900) MAXFCK,MXDIIS,MXHFMA, THRRHF
      END IF
      IF (DOCCSD) CALL CC_PRTINP(IWUNIT)
 1600 FORMAT(/5X,'Maximum number of macro iterations',I8,
     &       /5X,'Maximum number of micro iterations',I8,
     &       /5X,'Threshold for MCSCF gradient    ',1P,D10.2,
     &       /5X,'Number of initial trial vectors   ',I8,
     &       /5X,'Number of initial CI iterations   ',I8,
     &       /5X,'Number of simultaneous trial vectors',I6)
 1800 FORMAT(/5X,'Maximum number of CI iterations ',I10,
     &       /5X,'Threshold for CI convergence    ',1P,D10.2,
     &       /5X,'Number of CI roots to converge  ',I10)
 1900 FORMAT(/5X,'Maximum number of Fock   iterations',I7,
     &       /5X,'Maximum number of DIIS   iterations',I7,
     &       /5X,'Maximum number of QC-SCF iterations',I7,
     &       /5X,'Threshold for SCF convergence ',1P,D12.2)
 2500 FORMAT(/5X,'Maximum number of NEO/NR iterations',I4)
 2550 FORMAT(/5X,'CI trial vectors are generated with the Olsen',
     &       ' algorithm.')
 2560 FORMAT(/5X,'Dimension of CI PHP subspace' :,I6)
 2570 FORMAT(/5X,
     &   ' CI PHP subspace selected from size of residual elements.')
 2600 FORMAT(/'@    Reference state in CI response calculation :',I8)
C
C ***
C
      NFROT   = ISUM(NSYM,NFRO,1)
      NFRRHFT = ISUM(NSYM,NFRRHF,1)
      NFRMP2T = ISUM(NSYM,NFRMP2,1)
      IF (DOSCF .AND. NFRHFT .GT. 0) THEN
         IF (DODFT) THEN
            WRITE (IWUNIT,6110) NFRRHFT,(NFRRHF(I),I=1,NSYM)
         ELSE IF (DOHFSRDFT) THEN
            WRITE (IWUNIT,6112) NFRRHFT,(NFRRHF(I),I=1,NSYM)
         ELSE
            WRITE (IWUNIT,6111) NFRRHFT,(NFRRHF(I),I=1,NSYM)
         END IF
      END IF
      IF (DOMP2 .AND. NFRMP2T .GT. 0) THEN
         IF (NFRMP2T .EQ. 0 .AND. NFRRHFT .GT. 0) THEN
            WRITE (LUPRI,'(//A/A)')
     &    '@INFO: Frozen SCF orbitals must also be frozen in MP2.',
     &    '@INFO: Therefore frozen MP2 orbitals are set'//
     &    ' equal to frozen SCF orbitals (default).'
            NFRMP2(1:8) = NFRRHF(1:8)
            NFRMP2T     = NVRRHFT
         END IF
         WRITE (IWUNIT,6120) NFRMP2T,(NFRMP2(I),I=1,NSYM)
         IF (NFRMP2T .LT. NFRRHFT) THEN
            WRITE (IWUNIT, '(/A)') '@ FATAL ERROR: '//
     &      'a frozen SCF orbital must also be frozen in MP2'
            CALL QUIT(' FATAL ERROR: '//
     &      'a frozen SCF orbital must also be frozen in MP2')
         END IF
      END IF
      IF (DOMC .AND. NFROT.GT.0) THEN
         WRITE (IWUNIT,6100) NFROT,(NFRO(I),I=1,NSYM)
      END IF
 6110 FORMAT(//'@    Number of frozen SCF orbitals',T38,I5,' |',8I5)
 6111 FORMAT(//'@    Number of frozen KS orbitals',T38,I5,' |',8I5)
 6112 FORMAT(//'@    Number of frozen srHF orbitals',T38,I5,' |',8I5)
 6120 FORMAT(//'@    Number of frozen MP2 orbitals',T38,I5,' |',8I5)
 6100 FORMAT(//'@    Number of frozen MCSCF orbitals',T38,I5,' |',8I5)
C
      NNOROT = 0
      DO 6200 I = 1,NORBT
 6200    IF (NOROT(I).NE.0) NNOROT = NNOROT + 1
      IF (NNOROT.GT.0) THEN
         LNOROT = .TRUE.
         WRITE (IWUNIT,6300) NNOROT
         DO 6250 I = 1,NORBT
         IF (NOROT(I) .NE. 0) THEN
            ISYMI = ISMO(I)
            ITYPI = IOBTYP(I)
            II    = I - IORB(ISYMI)
            WRITE (IWUNIT,6350) I,ISYMI,II,COBTYP(ITYPI)
         END IF
 6250    CONTINUE
      END IF
 6300 FORMAT(/'@    Number of orbitals frozen with ".FREEZE"',I3,
     &      //'     Orbital  Symmetry   No. in   Type',
     &       /'       no.              symmetry')
 6350 FORMAT(3I10,3X,A9)
C
 8000 CONTINUE
C
C *** End of subroutine SIR_PRTINP
C
      CALL QEXIT('SIR_PRTINP')
      RETURN
      END
C  /* Deck nevpt2inp */
      SUBROUTINE NEVPT2INP(WORD,INPERR,ALLOPT)
C
C  2003 -- input for NEVPT2
C
#include "implicit.h"
      PARAMETER (NTABLE = 3)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7
      LOGICAL ALLOPT
C
C Used from common:
C
C GNRINF : THR_REDFAC
C INFORB : NSYM
C INFMP2 : NFRMP2(*)
C INFPRI : IPRMP2
C
#include "maxorb.h"
#include "priunit.h"
#include "gnrinf.h"
#include "inforb.h"
#include "infpri.h"
#include "infpt2.h"
C
      DATA TABLE /'.FROZEN','.THRESH','.STATE'/
C
C     ***** PROCESS INPUT FOR NEVPT2INP *****
C
      do i=1,nsym
         NFRNEVPT2(i) = 0
      enddo
      THRNEVPT = 0.d0
      WORD1    = WORD
      ISTNEVCI = 0
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
         CALL UPCASE(WORD)
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102,103), II
               END IF
 1102       CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               IF (.NOT.ALLOPT)
     &            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
               GO TO 100
            END IF
            WRITE (LUPRI,'(/4A/)') ' Keyword ',WORD,
     *         ' not recognized for ',WORD1
         ELSE IF (PROMPT .EQ. '#' .OR. PROMPT .EQ. '!') THEN
            GO TO 100
         ELSE IF (PROMPT .EQ. '*') THEN
            GO TO 8000
         ELSE
            WRITE (LUPRI,'(/3A/2A/)')
     *         ' Keyword ',WORD,' does not begin with',
     *         ' one of the four characters ".*!#" for ',WORD1
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal keyword for "*NEVPT2"')
C
C *** Option 1 >NEVPT2 FROZEN<   Frozen NEVPT2 orbitals
  101 CONTINUE
         READ(LUCMD,*) (NFRNEVPT2(I),I=1,NSYM)
      GO TO 100
C
C *** Option 2 >THRESHOLD < Threshold for the determinant coefficients for
C                           the calculation of density matrices in koopro4
  102 READ(LUCMD,*) THRNEVPT
      GO TO 100
C
C *** Option 3 >ISTNEVCI < If .CI, state considered for NEVPT2
C
  103 READ(LUCMD,*) ISTNEVCI
      GO TO 100
C
 8000 CONTINUE
      IF (THR_REDFAC .GT. 0.0D0) THEN
C         ICHANG = ICHANG + 1
         WRITE (LUPRI,'(3A,1P,D10.2)') '@ INFO ',WORD1,
     &   ' thresholds multiplied with general factor',THR_REDFAC
         THRNEVPT = THRNEVPT*THR_REDFAC
      END IF
      END
C --- end of sirinp.F ---
