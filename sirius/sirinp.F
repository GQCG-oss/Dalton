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
#include <single.h>
C
#ifdef UNDEF
===========================================================================
Revision 1.9  2000/05/24 12:39:05  hjj
1) only MAXCIT=3 as default if MC combination with DORHF or DOMP2
2) polish of output

Revision 1.8  2000/05/02 12:45:25  hjj
polishing output

Revision 1.7  2000/05/01 14:11:43  hjj
removed obsolete KAVER, KSYM(4); removed .AVERAG from input

revision 1.6 2000/04/14 09:11:03 hjj
removed .IPRERR in PRIINP (IPRERR is now set under **DALTON) /hjaaj
----------------------------
revision 1.5 2000/04/14 08:32:47 hjj
1) CALL QUIT if VAR_NONAMELIST and .FLAGS specified
   (removed special do nothing for T3D, seemed inconsistent ???)
2) removed old not used format statements in input routines
----------------------------
revision 1.4 2000/04/14 08:27:22 hjj
/hjaaj dec 1999 - apr 2000:
1) removed obsolete RNALL references
2) do not activate QCHF (flag(21)) for solvent or writing interface
if MXDIIS.gt.0; solvent and CALL WRSIFC are now also implemented
in the DIIS module.
3) removed absolete input tests on NACTEL and ISPIN
4) if not QCHF then zero iteration limits, for consistent output in PRTINP
5) BUG FIX: hucdrv was not called for start orbitals when
         not automatic occupation (when HF OCC specified in input)
6) fix: corrected output error in max # of QCHF iteration
7) polishing ...
----------------------------
revision 1.3 2000/04/14 08:23:47 hjj
removed "OLDINP", old obsolete input for SIRIUS /hjaaj apr 2000
----------------------------
revision 1.2 2000/04/14 07:07:51 hjj
1) initialize LUSUPM = -9999, to signalize default is use AOSUPMAT
   (before GPOPEN this was done with LUSUPM=17 in SIRINI)
2) changed many LUERR to LUPRI in SIRRDI
=============================================================================
960905-hjaaj
RHFINP : new '.GRAM-S', old '.GRAM S' kept for compatibility
951017 C*FM
PRTINP : FM change format nb of determinantsi and CSF  info 
Til Hans fra Hans Joergen, MANDAG 29. okt. har jeg aendret:
PRIINP : + .NOSUMM og .CANONI (for hvor egen skyld)
GENINP : namlist til NMLSIR i stedet for INPUT (saa det er
         i overensstemmelse med manualen).
901026-hjaaj:
-- indfoert common block SCBRHF options i SIRRDI, RHFINP
-- OPTINP: FLAG(66) = .TRUE. when .DETERMINANTS
-- indsat MOTECC-90 in-line documentation
-- SIRRDI: check CORHOL calculations
-- PRTINP: tell if core hole calculation or open shell HF
900924-hjaaj: GENINP: + .INTERFACE
-- WAVINP: rettet trykfejl s/NHLMN1/NHL1MN/ og s/NHLMX1/NHL1MX/
900920-hjaaj: ORBINP: s/READ M/MOSTAR/
-- OPTINP: space before LEVEL n
-- OPTINP: .COREHOLE keyword + .CORERELAX keyword
-- SIRRDI: CORHOL and CORRLX implemented, no check if calc. ok yet!
changes from 17/9 1990 (ha)
-- *CI INP to *CI CALC
-- *MP2 CALC
-- MP2INP routine included
      END
===========================================================================
#endif
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
#include <implicit.h>
      DIMENSION WORK(LFREE)
C
C -- local constants
C
#include <dummy.h>
C
C -- common blocks
C    INFINP : FLAG(*),?
C    INFORB : NCMOT,...
C    INFPRI : LUPNCH
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
C
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
C
      IF (IREST.EQ.0 .OR. IORTO.EQ.1) THEN
         KW1=1
         KW2=KW1+NCMOT
         KW3=KW2+NNBAST
         IF (KW3.GT.LFREE) CALL ERRWRK('SIRINP.ORTHO',-KW3,LFREE)
         LW3 = LFREE - KW3
         CALL ORTHO(WORK(1),WORK(KW2),WORK(KW3),LW3)
C        CALL ORTHO(CMO,S,SIN,LSIN)
      END IF
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
      ELSE IF (IREST .EQ. 0) THEN
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
C  IREST, = 1 if restart
C
C Called from SIRINP
C
C Subroutine calls: MOLLAB (FIND RECORD ON DATA FILES)
C                   PRTINP (PRINT INPUT)
C                   ERRINP (ERROR message OUTPUT ROUTINE)
C
C******************************************
C
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
      DIMENSION CMO(LCMO)
C
C -- local variables and constants
C
      CHARACTER WORD*8
      LOGICAL   DELEMO, FNDLAB
      PARAMETER (D0=0.0D0, D1=1.0D0, THFMP2=1.0D-6)
C
C Used from common blocks:
C   INFINP : ...,ITRLVL,ITRFIN,..., DO*, DIRFCK,MXDIIS,...
C   INFORB : NSYM,?
C   INFVAR : MAXOCC,NWOPT,JWOP(2,*),JWOPSY
C   INFOPT : EPOT,?
C   INFDIM : MAXPHP,?
C   INFTRA : THRP,THRQ,USEDRC,IPRTRA
C   SCBRHF : IOPRHF, NFRRHF(8), MXHFMA, MXHFMI, RHFCAN
C   INFMP2 : NFRMP2(*)
C   HUCKEL : ADDSTO
C
#include <maxash.h>
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
#include <infopt.h>
#include <inftap.h>
#include <infpri.h>
#include <infdim.h>
#include <inftra.h>
#include <scbrhf.h>
#include <infmp2.h>
#include <gnrinf.h>
#include <siripc.h>
#include <huckel.h>
#include <dftcom.h>
#include <dfterg.h>
C
      DIMENSION MBAS(8)
      LOGICAL ONELOK
C
C -- DATA statements
C
      DATA INPERR/0/, NUMRUN /0/
C
C INPERR = number of input errors. Stop if different from zero.
C
      CALL QENTER('SIRRDI')
C
C  ***** Check dimension
C
      IF ( LCMO .LT. MXCORB*MXCORB )
     *   CALL ERRWRK('SIRRDI.cmo',-MXCORB*MXCORB,LCMO)
C
C     Skip next section if processed earlier
C
      IF (INPPRC) THEN
C
C     Several common block are restored
C     (MUST be called before NEWSYM reset to false)
C
         CALL SRESTR
         IF (RNRESP) THEN
            FLAG(25) = .TRUE.
            FLAG(11) = .TRUE.
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
            IF (DORHF) THEN
               AUTOCC = .TRUE.
               INIOCC = 0
               DO 13 I = 1, 8
                  NRHF(I) = 0
  13           CONTINUE
               IOPRHF = 0
            ELSE
               WRITE (LUPRI,'(/A)')
     &            ' Lowering of symmetry for MCSCF not'//
     &            ' implemented in SIRIUS yet, sorry!'
               CALL QUIT('Lowering of symmetry for MCSCF not'//
     &            ' implemented in SIRIUS yet, sorry!')
            END IF
            NEWSYM = .FALSE.
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
         IF (ISTATE.GT.1) FLAG(39) = .TRUE.
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
         IF (WORD(1:8) .NE. '**SIRIUS' .AND.
     &       WORD(1:8) .NE. '**WAVE F' .AND.
     &       WORD(1:8) .NE. '**WAVEFU') THEN
            GO TO 900
         END IF
         GO TO 930
  910    CONTINUE
            CALL QUIT('End of file on LUCMD, no Sirius input found')
  920    CONTINUE
            CALL QUIT('Error reading LUCMD, no Sirius input found')
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
C
C
C     Jump again if processed earlier
C
      IF (INPPRC) GOTO 124
C
      IREST = 0
C     IREST = 1 for a restart run
C
C********** Default values section GENINP *****
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
      DORHF = .FALSE.
      DOMP2 = .FALSE.
      DOCINO= .FALSE.
      DOCI  = .FALSE.
      DOMC  = .FALSE.
      DODFT = .FALSE.
      DORSP = RNRESP
      ABAIPH= RNABAC
C     ... It's followed by an abacus run
      CORHOL= .FALSE.
      CORRLX= .FALSE.
C
C********** Default values section HAMINP *****
C
      NFIELD = 0
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
C     THRGRD   Threshold for convergence on gradient in MC/QC-HF opt.
      THRRHF = THFMP2
      THRCI  = 1.D-05
      THRMC  = 1.D-05
      THRGRD = 1.D-05
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
C
C     Initialize to zero orbital rotations
      NWOPT = 0
C THRPWF: threshold for print CI-coefficient in PRWF
      THRPWF = -1.0D0
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
C
C********** Default values section CICINP *****
C
C MAXRC : Max dimension of reduced CI matrix
      MAXRC  = 150
      ISTACI = 0
C     ... do not converge to a specific state.
      NROOCI = 1
      MXCIMA = 20
      ICICNO = 0
C
C********** Default values section RHFINP *****
C
      IOPRHF = 0
      DIRFCK = DIRCAL
      MXDIIS = 60
      MAXEVC = -1
      RHFCAN = .TRUE.
      MAXFCK = 0
      MXHFMA = 15
      MXHFMI = 12
C 960625-hjaaj: NISHT and NACTEL used in MOLCHR in
C                       ICHRG = ICHRGNU - (2*NISHT+NACTEL)
      NISHT  = 0
      NACTEL = 0
      CALL MOLCHR(ICHRG)
      NACTEL = -1
      NRHFEL = ICHRG - KCHARG
      INIOCC = -1
      IF (NRHFEL .GT. 0) THEN
         AUTOCC = .TRUE.
      ELSE
         AUTOCC = .FALSE.
      END IF
C
C********** Default value sections RHFINP and AUXINP *****
C           (default: use AOSUPMAT if available,
C                     LUSUPM=-1 means do not use AOSUPMAT)
C
      LUSUPM = -9999
C
C********** Default values section DFTINP *****
C
      CALL DFTINI
C
C
C********** Default values section ORBINP *****
C
C     LMOORD   True if mo:s have to be reordered
C
      LMOORD = .FALSE.
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
      SUPSYM = .NOT. (RNABAC .OR. RNRESP) 
C     supsym false forced for DOCCSD set after reading of sirius input
C
C IRDMO : select how to get starting orbitals
C         (default 8, diagonalize one-el. Ham.
C          corresponding to .MOSTART( H1DIAG )
C KDEL  : if .gt. 0, then delete 3s component in d; 4p in f; etc.
C         If IRDMO .eq. 8, i.e. new calculation, call CMODEL
C         to find 3s,4p,... components and delete them.
C         if .lt. 0, orbitals have been deleted in input
C THROVL: limit for basis set numerical linear dependence
C
      IRDMO  = 8
      FLAG(13) = .FALSE.
      KDEL   = 0
      THROVL = 1.0D-6
      CMAXMO = 1.0D4
      FCVORB   = .FALSE.
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
      THRP   = 1.0D-15
      THRQ   = 1.0D-15
      ITRLVL = -6789
      ITRFIN = -6789
      IPRTRA = -1
C     MWORK is defined in block data SBDDIM
C
      USEDRC= .FALSE.
C     ... USEDRC may be reset by transformation module
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
      IPRDIA = IPRUSR
      IPRSIG = IPRUSR
      IPRDNS = IPRUSR
      IPRRHF = IPRUSR - 1
      IPRMP2 = IPRUSR + 5
      IPRSOL = IPRUSR + 5
      IPRCNO = IPRUSR + 1
      IPRFCK = IPRUSR
C
C  ***************************************************
C  *****        Read Sirius input                *****
C  ***************************************************
C
C     Extra default: CALL SIROUT for summary
      FLAG(5) = .TRUE.
C
      CALL NEWINP(INPERR,NUMRUN,IRDMO,IREST,NSYM)
      IF (.NOT.DORHF) MAXFCK = 0
C
      IF (DORHF .AND. .NOT.FLAG(21)) THEN
        IF (MXDIIS .LE. 0) THEN
C        ... March 2000/hjaaj: solvent and writing interface are
C            now also defined in DIIS module.
          IF (FLAG(16)) THEN
            NINFO = NINFO + 1
            WRITE (LUPRI,'(/5X,A/5X,A)')
     &      'SIRRDI INFO, ".NOQCHF" ignored because solvent is only',
     &      '             implemented in DIIS and QCHF modules.'
            FLAG(21) = .TRUE.
          ELSE IF (FLAG(25) .AND. .NOT.(DOMC .OR. DOCI)) THEN
            NINFO = NINFO + 1
            WRITE (LUPRI,'(/5X,A/5X,A)')
     &      'SIRRDI INFO, ".NOQCHF" ignored because writing interface',
     &      '             is only implemented in DIIS and QCHF modules.'
            FLAG(21) = .TRUE.
          END IF
        END IF
        IF (.NOT.FLAG(21)) THEN
C       ... if not QCHF zero iteration limit (for consistent output in PRTINP)
          MXHFMA = 0
          MXHFMI = 0
        END IF
      END IF
C
      IF (DOCCSD) THEN
         FLAG(11) = .TRUE.
         IF (DORHF) FLAG(25) = .TRUE.
         SUPSYM = .FALSE.
      END IF
      IF (RNRESP .or. RNABAC) THEN
         FLAG(25) = .TRUE.
         FLAG(11) = .TRUE.
      END IF
      IF (FLAG(21)) FLAG(2) = .TRUE.
C     ... qc-scf needs flag(2) for file opening
C
C *** Open print unit(s).
C     and initialize units for memory allocation routines
C
      IF (LUW4 .NE. LUPRI) THEN
         LUW4 = 0
         CALL GPOPEN(LUW4,'SIRIUS.RES','UNKNOWN',' ','FORMATTED',IDUMMY,
     &               .FALSE.)
      END IF
C
C     set THRGRD for test of value below
C     set THQKVA if not specified in input
C
      IF (DOMC) THEN
         THRGRD = THRMC
         IF (THQKVA .LE. D0) THQKVA = 8.0D0
      ELSE IF (DOCI) THEN
         THRGRD = THRCI
      ELSE IF (DORHF) THEN
         THRGRD = THRRHF
         IF (THQKVA .LE. D0) THQKVA = 0.8D0
      END IF
C
C     Several common blocks are saved to next iteration
C
      CALL SSAVE
C
C     This is where we jump back in if input has been processed earlier.
C
      INPPRC = .TRUE.
 124    CONTINUE
      IF (DIRFCK) THEN
         LUSUPM = -1
         IF (DOMP2 .OR. DOCI .OR. DOMC) THEN
            WRITE (LUPRI,'(/A/A)')
     &      'INPUT ERROR: Direct Fock matrix construction',
     &      '             not possible for MP2, CI, and MCSCF.'
            INPERR = INPERR + 1
         END IF
      END IF
C
C     set ITRFIN for ABAIPH
      IF (ABAIPH .AND. ITRFIN .EQ. -6789) THEN
         IF (DOMC .OR. DOCI) ITRFIN = 4
      END IF
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
         WRITE (LUW4,'(2(/5X,A),/A,I3)')
     &  'WARNING, reference state is required for CI response and CINO',
     &  '         ".STATE" was not specified for "*CI INPUT"',
     &  '         reference state is set to be state no.',ISTACI
      END IF
      END IF
C
C     Save explicitly defined IPRSIR, IPRI4, IPRI6
C
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
      IF (MCTYPE .EQ. 0 .AND. (DOMC .OR. DOCI .OR. DOCINO)) THEN
         WRITE (LUPRI,'(/A)')
     &      'INPUT ERROR: MCSCF or CI requested,'//
     &      ' but neither CAS nor RAS space defined.'
         INPERR = INPERR + 1
      END IF
C
C     Disable AUTOCC for RHF if RHF not used (HJ aug95)
C
      IF (MCTYPE .GT. 0 .AND. .NOT.DORHF) THEN
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
C     even if DORHF false and DOMP2 false; JWOP information
C     may be needed for RESPONSE
C
      IF (MCTYPE .EQ. 0) THEN
         CALL ICOPY(8,NFRRHF,1,NFRO,1)
         CALL ICOPY(8,NRHF,1,NISH,1)
         CALL IZERO(NASH,8)
         IF (IOPRHF .NE. 0) THEN
            LSYM   = IOPRHF
            ISPIN  = 2
            NASH(LSYM) = 1
            NASHT  = 1
            NACTEL = 1
         ELSE
            LSYM   = 1
            ISPIN  = 1
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
C
C     New test to avoid conflicting input for subsequent CCSD runs
C
      IF (NASHT .GE. 1 .AND. DOCCSD) THEN
         WRITE (LUPRI,*) 'Only closed-shell RHF wave functions'//
     &        ' can be used for CCSD calculations'
         INPERR = INPERR + 1
      END IF
C
C930406-hjaaj: if (IREST.eq.1) then reset options here and not in
C     geninp as before, where something may not be defined yet
C---------------
      IF (IREST .EQ. 1) THEN
         IRDMO = 9
         ICI0  = 5
         FLAG(13) = .TRUE.
         LMOORD   = .FALSE.
C
C	 If we restart with DIIS or Roothan Hartree-Fock, then we
C        probably don't have trust radius infomation on file
C        This statement disables the reading of such info
C
         IF (MAXFCK .GT. 0 .OR. MXDIIS .GT. 0 .AND. .NOT.
     &       (NASHT .GT. 1)) ICI0 = 1
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
      END IF
C
C     Set print parameters in /INFTRA/ for TRACTL
C     (only set IPRTRA if not set explicitly in input)
C
      IF (IPRTRA .EQ. -1) IPRTRA = IPRSIR
C
C     If MAXCIT, MAXMAC, MAXJT, MAXMIC not specified in input:
      IF (MAXCIT .LT. 0) THEN
         IF (DORHF .OR. DOMP2 .OR. DOCINO .OR. .NOT.FLAG(13)) THEN
            MAXCIT = 3
         ELSE
C           ... only MC and reading orbitals from file:
C               these may be converged orbitals thus MAXCIT = 20
            MAXCIT = 20
         END IF
      END IF
      IF (MAXMAC .LT. 0) THEN
         IF (ISTATE .EQ. 1) THEN
            MAXMAC = 15
         ELSE
            MAXMAC = 25
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
C     if "MOSTRT 8" then call cmodel to
C        kdel .eq. 0: define unit matrix
C        kdel .gt. 0: delete 3s in d, etc.
C        kdel .lt. 0: DELETE specified, MOSTRT 8 illegal
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
         DELEMO = .FALSE.
         CALL DELMO(CMO(1),CMO(KSCRA),LSCRA,THROVL,DELEMO)
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
C        Find Huckel start orbitals if ADDSTO and not AUTOCC
C        (If automatic occupation HUCDRV is called later in this routine)
C
         IF (ADDSTO .AND. .NOT.AUTOCC)
     &      CALL HUCDRV(CMO(KSCRA),LSCRA,CMO)
C
         FLAG(13) = .FALSE.
      ELSE IF (FLAG(13)) THEN
C
C        ... read old orbitals from file
         NERROR = 0
         NDELT = 0
         IF (IRDMO .EQ. 7 .OR. IRDMO .EQ. 9) THEN
C           ... read initial MO's from LUIT1
C               Find out if any orbitals deleted
            CALL GPOPEN(LUIT1,'SIRIUS.RST','OLD',' ','UNFORMATTED',
     &                  IDUMMY,.FALSE.)
            REWIND LUIT1
            IF (FNDLAB('BASINFO ',LUIT1)) THEN
C           ... 921202: old LUIT1 files do not have BASINFO
               IF (INIOCC .LE. 0) THEN
                  READ (LUIT1) MSYM,MBAS,NORB,NRHF,IOPRHF
                  INIOCC = 2
               ELSE
                  READ (LUIT1) MSYM,MBAS,NORB
               END IF
               NERROR = -1
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
     &               JSTATE,JSPIN,MACTEL,LXSYM,KXAVER
            IF (INIOCC .LE. 0 .AND. MACTEL .LE. 1) THEN
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
            END IF
C              READ (LUSIFC) POTNUC,EMY,EACTIV,EMCSCF,
C    *               ISTATE,ISPIN,NACTEL,LSYM
C              READ (LUSIFC) NISHT,NASHT,NOCCT,NORBT,NBAST,NCONF,NWOPT,
C    *               NWOPH,NCMOT,NNASHX,NNASHY,NNORBT,N2ORBT,
C    *               NSYM,MULD2H, NRHF,NFRO,
C    *               NISH,NASH,NORB,NBAS,
C    *               NELMN1, NELMX1, NELMN3, NELMX3, MCTYPE,
C    *               NAS1, NAS2, NAS3
#if defined (VAR_SECSEC)
hvordan ved vi om vi vil have samme beregning,
eller om vi vil bruge orbitalerne ved f.eks. en ny geometri
eller en ny tilstand
               NERR = 0
               IF (XPOTNC .NE. POTNUC) NERR = NERR + 1
               IF (NERR .GT. 0) THEN
                  INPERR = INPERR + 1
               END IF
#endif
            ELSE
               INPERR = INPERR + 1
               WRITE (LUPRI,'(/2A/3A)')
     &         ' Input error, MO coefficients to be read from'//
     &         ' SIRIUS interface file ',FNSIFC,' but ',
     &         FNSIFC,' with proper information was not provided.'
            END IF
            CALL GPCLOSE(LUSIFC,'KEEP')
         END IF
         IF (NERROR .EQ. -1) THEN
C           test NSYM,NBAS information from LUIT1 or LUSIFC
            NERROR = 0
C           IF (MSYM .NE. NSYM) THEN
C              NERROR = NERROR + 1
C           ELSE
C              DO 40 ISYM = 1,NSYM
C                 IF (NBAS(ISYM) .NE. MBAS(ISYM)) NERROR = NERROR + 1
C                 NDELT = NDELT + (NBAS(ISYM) - NORB(ISYM))
C  40          CONTINUE
C           END IF
         END IF
         IF (NERROR .GT. 0) THEN
            INPERR = INPERR + 1
            WRITE (LUPRI,'(/A/5X,A,I3,5X,8I5)')
     &         'Input error in connection with ".MOSTART"',
     &         'NSYM, NBAS from input       :',NSYM,(NBAS(I),I=1,NSYM)
            WRITE (LUPRI,'(5X,A,I3,5X,8I5)')
     &         'NSYM, NBAS from MOSTART file:',MSYM,(MBAS(I),I=1,MSYM)
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
         IF (DORHF .AND. AUTOCC) THEN
            INPERR = INPERR + 1
            WRITE (LUPRI,*) 'Input error for core-hole calculation:'
            WRITE (LUPRI,*)
     &         'automatic occupation of HF orbitals not allowed.'
         END IF
         IF (IRDMO .EQ. 8) THEN
            INPERR = INPERR + 1
            WRITE (LUPRI,*) 'Input error for core-hole calculation:'
            WRITE (LUPRI,*) 'no initial molecular orbitals from'//
     &         ' previous calculation specified'
         END IF
         IF (DORHF .AND. DOMC) THEN
            INPERR = INPERR + 1
            WRITE (LUPRI,*) 'Input error for core-hole calculation:'
            WRITE (LUPRI,*) 'both Hartree-Fock and MCSCF specified.'
         END IF
         IF (NASHT .GT. 1) THEN
          IF (NELMN1 .EQ. 0 .AND. NELMX1 .EQ. 0) THEN
            NELMN1 = 1
            NELMX1 = 1
          END IF
          IF (NELMN1 .NE. NELMX1 .OR. MCTYPE .NE. 2 .OR.
     &        NAS1(JCHSYM) .NE. 1 .OR. NAS1T .NE. 1) THEN
            INPERR = INPERR + 1
            WRITE (LUPRI,*) 'Input error(s) for core-hole calculation:'
            WRITE (LUPRI,*) '-----------------------------------------'
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
               WRITE (LUPRI,*) 'Input error: '//
     &         'core-hole sym. and RHF open shell sym. different'
            END IF
         ELSE
            INPERR = INPERR + 1
            WRITE (LUPRI,*) 'Input error: '//
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
                  WRITE (LUPRI,'(/A/A)')
     &               ' Input error for core-hole calculation',
     &               ' Core hole orbital is not inactive'
            ELSE IF (JCHORB .LE. NISH(JCHSYM)) THEN
               IF (LMOORD) THEN
                  INPERR = INPERR + 1
                  WRITE (LUPRI,'(/A/A)')
     &               ' Input error for core-hole calculation',
     &               ' .COREHOLE and .REORDER are not compatible'//
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
         CALL GPOPEN(LUIT1,'SIRIUS.RST','OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         CALL READMO(CMO,JRDMO,CMO(K1),LSCR)
         CALL GPCLOSE(LUIT1,'KEEP')
C        JRDMO negative : return with error status in JRDMO
C        JRDMO positive : abort in READMO if error
         IF (JRDMO .GT. 0) THEN
            SUPSYM = .FALSE.
            WRITE(LUPRI,'(/A/A)')
     &         ' SIRRDI input error, not able to read mo coefficients.',
     &         ' "SUPSYM" disabled because mo coeff. are not'//
     &         ' available for analysis.'
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
C              and Huckel start orbitals if ADDSTO
C
      IF (NRHFEL .LE. 0) THEN
         NRHFEL = 2*NRHFT
         IF (IOPRHF .GT. 0) NRHFEL = NRHFEL + 1
      END IF
      IF (DORHF .AND. AUTOCC .AND. IREST .EQ. 0
     &   .AND. INPERR .EQ. 0) THEN
C
C        Initial occupation may have been specified in input
C        or read from LUIT1 or LUSIFC. If not, we find it now
C        from diagonal elements of the H1 matrix.
C        If ADDSTO use Huckel guess (INIOCC = 5).
C
         KSCRA = 1    + NCMOT
         LSCRA = LCMO - KSCRA
         IF (INIOCC .LE. 0) THEN
            IF (ADDSTO) THEN
               CALL HUCDRV(CMO(KSCRA),LSCRA,CMO(1))
               FLAG(13) = .FALSE.
               INIOCC = 5
            ELSE
               CALL H1OCC(CMO,CMO(KSCRA),1,LSCRA)
               INIOCC = 4
            END IF
         ELSE IF (INIOCC .EQ. 1 .AND. ADDSTO) THEN
            AUTOCC = .FALSE.
C           ... otherwise HUCDRV will destroy user input
C               of initial occupation guess
            CALL HUCDRV(CMO(KSCRA),LSCRA,CMO(1))
            FLAG(13) = .FALSE.
            AUTOCC = .TRUE.
         END IF
      END IF
C
C       
C     If INIOCC .ge. 2 and RHF/MP2 highest level (i.e. MCTYPE .eq. 0),
C     we need to retransfer HF occupation to MCSCF parameters as
C     NRHF(*),IOPRHF may have changed since previous RHF setting above.
C     (INIOCC .ge. 2: NRHF(*) from LUIT1, LUSIFC, H1OCC, or HUCDRV)
C
      IF (INIOCC .GE. 2 .AND. MCTYPE.EQ.0) THEN
         CALL ICOPY(8,NRHF,1,NISH,1)
         CALL IZERO(NASH,8)
         IF (IOPRHF .NE. 0) THEN
            LSYM   = IOPRHF
            ISPIN  = 2
            NASH(LSYM) = 1
            NASHT  = 1
            NACTEL = 1
         ELSE
            LSYM   = 1
            ISPIN  = 1
            NASHT  = 0
            NACTEL = 0
         END IF
         CALL SETORB
      END IF
C
C
C *** Check if initial occupation for RHF has been defined.
C
      IF ((DORHF .OR. DOMP2 .OR. DOCCSD) .AND. INIOCC .LE. 0) THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,'(/A/A)')
     &      ' Input error for SCF, MP2, or CC calculation:',
     &      ' Initial SCF occupation has not been defined.'
         INIOCC = 0
C        ... set INIOCC = 0 for PRTINP
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
C If (MP2) force HF convergence to THFMP2
C If (MP2) force canonical Hartree-Fock orbitals
C
      IF (DOMP2) THEN
      IF (DORHF .AND. IOPRHF .NE. 0) THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,'(/A)') ' Input error: closed shell MP2'//
     &      ' not compatible with open shell RHF'
      END IF
      IF (THRRHF .GT. THFMP2) THEN
         NINFO = NINFO + 1
         WRITE (LUPRI,'(/5X,2(A,1P,D10.2)/A)')
     &  'SIRRDI INFO, THRRHF reset from',THRRHF,' to',THFMP2,
     &  '             because MP2 requires converged orbitals'
         THRRHF = THFMP2
      END IF
      IF (.NOT.RHFCAN) THEN
         NINFO = NINFO + 1
         WRITE (LUPRI,'(/5X,A)')
     &  'SIRRDI INFO, ".NONCANONICAL" ignored because MP2 requires'//
     &  ' canonical orbitals.'
         RHFCAN = .TRUE.
      END IF
      END IF
C
C     Call AVESET to set up supersymmetry information
C     (If SUPSYM false, then supersymmetry information
C      will be a copy of the symmetry information)
C
      IF (INPERR .GT. 0 .AND. SUPSYM) THEN
         WRITE(LUPRI,'(/A)') ' NOTE: super symmetry analysis'//
     &   ' not possible because of input errors.'
         SUPSYM = .FALSE.
      END IF
      K1 = 1 + NCMOT
      KFREE = 1
      LFREE = LCMO - K1
      CALL AVESET(CMO,CMO(K1),KFREE,LFREE)
C     CALL AVESET(CMO,WRK,KFREE,LFREE)
C
C
C *** Set up various index arrays etc.
C     If AVERAGE check if OK in AVECHK
C
      IF (IORTO .NE. 1) THEN
C     ... only SETORB needed for reorthonormalization and
C         SETORB has already been called.
         K1 = 1 + NCMOT
         JWOPSY = 1
         CALL SIRSET(CMO(K1),LCMO-K1,.FALSE.)
C        CALL SIRSET(WRK,LWRK,OLDWOP)
         CALL AVECHK(CMO(K1),LCMO-K1,INPERR)
      END IF
C
C
C ***  PRINT INPUT DATA *****
C
      CALL PRTINP(LUPRI,NUMRUN,IRDMO)
      IF (IREST.EQ.1) WRITE (LUPRI,3050)
      CALL FLSHFO(LUPRI)
      IF (LUW4 .NE. LUPRI .AND. P4FLAG(1)) THEN
         CALL PRTINP(LUW4,NUMRUN,IRDMO)
         IF (IREST.EQ.1) WRITE (LUW4,3050)
         CALL FLSHFO(LUW4)
      END IF
 3050 FORMAT(/5X,'This is a restarted calculation.')
C
C
C *****  PROCESS AND CHECK INPUT DATA *****
C
      IF(NSYM  .LE.0 .OR.NSYM  .GT.8)      CALL ERRINP( 4,INPERR)
      IF(NORBT .LT.0 .OR.NORBT .GT.MXCORB) CALL ERRINP( 8,INPERR)
      IF (IORTO .EQ. 1) GO TO 8000
C
      IF (NUMRUN .GT. 0) THEN
         IF(MAXMAC.LE.0 .OR.MAXMAC.GT.200) CALL ERRINP( 5,INPERR)
         IF(THRGRD.LT.D0.OR.THRGRD.GT.D1)  CALL ERRINP( 7,INPERR)
      END IF
      IF(NRHFT .LT.0 .OR.NRHFT .GT.MAXOCC) CALL ERRINP( 9,INPERR)
      IF(NASHT .LT.0 .OR.NASHT .GT.MAXASH) CALL ERRINP(10,INPERR)
      IF(DORHF .AND. AUTOCC .AND. MCTYPE .GT. 0) THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,'(/A)') ' HF occupation must be specified when '//
     *         'SCF followed by CI or MCSCF'
      END IF
      IF(NOCCT .LT.1 .OR. NOCCT .GT.MAXOCC) THEN
         INPERR = INPERR + 1
         WRITE(LUPRI,'(/A,I5/A,I5)') 'Number of occ. orbitals  :',NOCCT,
     *      'is .lt. 1 or .gt. MAXOCC :',MAXOCC
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
 3101 FORMAT(/' SIRRDI input error, NFRO is larger than NISH for',
     *        I2,' symmetries.',
     *       /' (Only inactive orbitals can be frozen with ".FROZEN",',
     *        ' use ".FREEZE")')
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
     *      ' Input error in one or more solvent specifications:',
     *      ' Lmax                =',LSOLMX,
     *      ' Dielectric constant =',EPSOL,
     *      ' Static diel. const. =',EPSTAT,
     *      ' Cavity radius       =',RSOL(1),RSOL(2),RSOL(3)
         END IF
C
C        orbital absorption not safe
C
         IF (FLAG(51) .OR. FLAG(52) .OR. FLAG(53)) THEN
            WRITE (LUPRI,'(/2A/)') ' STOP because ABSORPTION is not',
     *         ' yet implemented for SOLVENT calculations'
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
     *      'Sirius input error: DELETE is illegal with MOSTRT 8'
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
         WRITE (LUPRI,'(/A/2A)')' SIRRDI, input error.',
     *      ' Wrong number of symmetries,',
     *      ' this cannot be a subgroup of D2h.'
         INPERR = INPERR + 1
      END IF
C
      IF (MOD((NACTEL+1),2) .NE. MOD(ISPIN,2) .OR.
     *   (NACTEL+1) .LT. ISPIN) THEN
         WRITE (LUPRI,3710)
         INPERR = INPERR + 1
      END IF
 3710 FORMAT(/' SIRRDI, input error.'/' Number of active electrons',
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
C *****  END OF INPUT PROCESSING *****
C
 8000 IF (INPERR .NE. 0) GO TO 9010
      IF (FLAG(1)) CALL QUIT('SIRIUS "STOP AFTER INPUT" STOP')
C
C
C *** Open other files
C
C
      CALL SIROPN
C
C
      CALL QEXIT('SIRRDI')
      RETURN
C
C *****  ERROR MESSAGES AND ERROR STOP *****
C
 9010 WRITE (LUERR,9000) INPERR
      WRITE(LUW4,9000) INPERR
      IF (LUPRI.NE.LUW4) WRITE(LUPRI,9000) INPERR
 9000 FORMAT(//20X,I5,' errors occured in SIRIUS input processing.'
     *        /20X,'Check error file for information about errors.'
     *        /20X,'Program cannot continue.'
     *        /20X,'Correct error(s) and be welcome back')
      CALL QTRACE(LUPRI)
      CALL QUIT('*** ERROR *** SIRIUS input errors.')
C
C     END OF SIRRDI.
C
      END
C  /* Deck newinp */
      SUBROUTINE NEWINP(INPERR,NUMRUN,IRDMO,IREST,NSYM)
C
C Revisions
C    1-OCT-1987 HA     Keyword directed input for all parameters
C
C OBJECTIVE: Read input data from file LUCMD
C            using new standard Sirius input
C
C
C Modified the routine to include Coupled-Cluster input routines.
C In addition GENINP is now always called, as it is required.
C *GENINP no longer included in input file
C
C Input:
C  none
C
C Output:
C  IREST, =1 if restart
C
C Called from SIRRDI
C
C
C Calls the following input modules:
C
C  *GENERAL
C  *HAMILTONIAN
C  *OPTIMIZATION
C  *HF INPUT                used to be *RHF CALCULATION
C  *MP2 INPUT               used to be *MP2 CALCULATION
C  *CONFIGURATION INPUT     used to be *WAVE FUNCTION
C  *CI INPUT                used to be *CI CALCULATION
C  *CC INPUT                used to be *CC
C  *ORBITAL INFORMATION
C  *CI VECTOR
C  *STEP CONTROL
C  *POPULATION ANALYSIS
C  *TRANSFORMATION
C  *RESPONSE CALCULATION
C  *UPDATE
C  *SOLVENT CALCULATION
C  *PRINT LEVELS
C  *AUXILLIARY INPUT
C  *END OF (marks end of keyword input)
C
C******************************************
C List of internal control flags for SIRIUS
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
C 8  SIRCTL  Response calculation (CALL RSPDRV)            GEN : RESPON
C 9  SIRCTL  CALL MP2CTL (OBS FLAG 3 will destroy)         GEN : MP2
C 10
C 11 SIRSET  Calculate non-symmetric orbital rotations
C            and write them to LUINDX (labeled)
C 12 GRAD    Exact orbital Hessian diagonal                WAV : FOCKDI
C                  (otherwise only Fock matrix contributions)
C 13 various CALL READMO for initial orbitals              ORB : vario
C 14 various Integral transformation not needed.           TRA : OLD TR
C 15 SIRSAV  transformation to NO's each macro iteration.  OPT : TRACI
C 16 SIROPT  Solvent calculation                           SOL : genera
C 17 various SUPSYM : find and use 'super symmetry'        ORB : NOSUPS
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
C 41 variuos OPT Optimal orbital trial vector (CALL NEXKAP)OPT : OPTIMA
C 42
C 43
C 44
C 45 SIRCAN  Redo CI with natural orbitals (CALL SIRCI)    OPT : CI CAN
C 46 SIRCAN  Only NO transformation, inactive and          OPT : NATONL
C    OPTST,SIRSAV sec. orbitals unchanged (no can. orbs)
C 47 SIRCAN  temp. - use CICTL, not TRACI)
C 48 SIRCAN,OPTST,SIRSAV CO's: diag. FD=FC+FV blocks       OPT : FOCKON
C 49
C 50
C 51 various Orb. absor. level 1;   inact.-act.,act.-act.  OPT : ABSORB
C 52 various Orb. absor. level 2; + inact.-sec.            OPT : ABSORB
C 53 various Orb. absor. level 3; +   act.-sec. (i.e. all) OPT : ABSORB
C 54 SIROPT  Do not switch absorption off in local region. OPT : ALWAYS
C 55 NEO,NR  include g(orb) H(diag)  ??
C 56 var.UPD Fixed Hessian                                 UPD : FIXED
C 57 UPDINI  Inverse Hess, diag. is initial approx. to     UPD : INVER
C            inverse Hess
C 58 CIST    take plus combination of degenerate roots     CIV : PLUS C
C 59
C 60
C 61
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
#include <implicit.h>
      LOGICAL SKIP
C
C -- local variables and constants
C
      PARAMETER (NTBMOD = 20)
      CHARACTER WORD*7, PROMPT*1, PROMPX*2, TABMOD(NTBMOD)*7,
     *          REWORD*12, RWORD*7, TITLIN*80
      LOGICAL   ALLOPT, PTABMO
C
#include <priunit.h>
#include <gnrinf.h>
C
      CALL QENTER('NEWINP')
C
      DATA TABMOD/'XXXXXXX','*CONFIG','*OPTIMI','*CI INP','*ORBITA',
     *            '*CI VEC','*STEP C','*POPULA','*TRANSF','*SCF IN',
     *            '*UPDATE','*SOLVEN','*PRINT ','*AUXILL','*HF INP',
     *            '*MP2 IN','*HAMILT','*CC INP','*DFT IN','*END OF'/
C
C
C  ***************************************************
C  ***** Branch to the *-marked general sections *****
C  ***************************************************
C
      ALLOPT = .FALSE.
      PTABMO = .FALSE.
C
      WORD = '*GENERA'
      CALL GENINP(WORD,INPERR,IREST,NUMRUN,ALLOPT)
C
  910 PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '*') THEN
         DO 920 I = 1, NTBMOD
            IF (TABMOD(I) .EQ. WORD) THEN
               GO TO (1100,1200,1300,1400,1500,1600,1700,1800,1900,
     *                2500,2100,2200,2300,2400,2500,2600,2700,2800,
     *                2500,8000), I
            END IF
 920     CONTINUE
         IF (WORD(1:2) .EQ. '**') GO TO 8000
         WRITE (LUPRI,'(/,3A,/)') ' Keyword ',WORD,
     *                        ' not recognized in Sirius input.'
         CALL PRTAB(NTBMOD,TABMOD,'General Sirius input keywords',LUPRI)
         CALL QUIT('ILLEGAL KEYWORD IN SIRIUS INPUT')
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
            GO TO 910
         END IF
         WRITE (LUPRI,'(/3A/)') ' Prompter "',PROMPT,'" illegal'
         CALL PRTAB(NTBMOD,TABMOD,'General Sirius input keywords',LUPRI)
         CALL QUIT('ILLEGAL PROMPT IN SIRIUS INPUT')
      END IF
C
 1100  CALL QUIT('Program error in NEWINP reading SIRIUS input.')
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
 1900  CALL TRAINP(WORD,INPERR,ALLOPT)
       GO TO 910
 2000  CONTINUE
       GO TO 1200
 2100  CALL UPDINP(WORD,INPERR,ALLOPT)
       GO TO 910
 2200  CALL SOLINP(WORD,INPERR,ALLOPT)
       GO TO 910
 2300  CALL PRIINP(WORD,INPERR,ALLOPT)
       GO TO 910
 2400  CALL AUXINP(WORD,INPERR,ALLOPT)
       GO TO 910
 2500  CALL RHFINP(WORD,INPERR,ALLOPT)
       GO TO 910
 2600  CALL MP2INP(WORD,INPERR,ALLOPT)
       GO TO 910
 2700  CALL HAMINP(WORD,INPERR,ALLOPT)
       GO TO 910
 2800  MSYM = NSYM
       CALL CC_INP(WORD,IREST,MSYM)
       IF (NSYM .NE. MSYM) THEN
          CALL QUIT(
     &       'CC_INP error: NSYM from SIRIFC .ne. NSYM in this run')
       END IF
       GO TO 910
C  Here if "*END OF"
 8000  CONTINUE
       IF (ALLOPT .AND. .NOT.PTABMO) THEN
          CALL PRTAB(NTBMOD,TABMOD,'General Sirius input keywords',
     &               LUPRI)
       END IF
C
C ********  END OF READ OF INPUT DATA  **********
C
C
C     We anyway need to process CC default settings
C
      JREST = 0
      CALL CC_INP('**WAVE F',JREST,NSYM)
C
      CALL QEXIT('NEWINP')
      RETURN
C
C     END OF NEWINP.
C
      END
C  /* Deck geninp */
      SUBROUTINE GENINP(WORD,INPERR,IREST,NUMRUN,ALLOPT)
C
C 1-Oct 1987 Hans Agren
C
#include <implicit.h>
      PARAMETER (NTABLE = 17)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7,
     *          REWORD*20, TITLIN*80
      LOGICAL ALLOPT
      PARAMETER (NPDIM = 50)
      INTEGER NPATH(NPDIM)
C
C Used from common:
C  INFINP : FLAG(*), IORTO, ICI0, DO*
C  INFORB : NSYM
C  INFPRI : P6FLAG,
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infpri.h>
#include <gnrinf.h>
C
#if !defined (VAR_NONAMELIST)
      NAMELIST /NMLSIR/ NPATH, MAXMIC
#endif
      DATA TABLE /'.CI    ','.FLAGS ','.MCSCF ','.CC    ','.RESTAR',
     *            '.HF    ','.STOP  ','.PRINT ','.TITLE ','.MP2   ',
     *            '.NSYM  ','.CC ONL','.WESTA ','.INTERF','.BASIS ',
     *            '.SCF   ','.DFT   '/
C
C
C
C
C     ***** PROCESS INPUT FOR GENINP *****
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
 1101    PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102,103,104,105,106,107,108,109,110,
     *                   111,112,113,114,115,106,116), II
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
            RETURN
         ELSE
            WRITE (LUPRI,'(/3A/2A/)')
     *         ' Keyword ',WORD,' does not begin with',
     *         ' one of the four characters ".*!#" for ',WORD1
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT(' ILLEGAL KEYWORD IN GENINP ')
C
C *** Option 1 >CI    <  only CI calculation
 101  CONTINUE
         DOCI = .TRUE.
         FLAG(4) = .TRUE.
         NUMRUN = NUMRUN + 1
      GO TO 100
C
C *** Option 2 >FLAGS <  Set Flags directly as in old input
C     set up control flags and print flags
C     read namelist $INPUT ... $END
 102  CONTINUE
         CALL IZERO(NPATH,NPDIM)
#if defined (VAR_NONAMELIST)
         CALL QUIT(
     &   'ERROR, .FLAGS namelist not implemented for this machine')
#else
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_NEC) || defined (SYS_T90)
         READ (LUCMD,NMLSIR)
#else
         READ (LUCMD,NML=NMLSIR)
#endif
#endif
C
         DO 286 I = 1,NPDIM
            IF (IABS(NPATH(I)) .GT. NFLAG) THEN
               INPERR = INPERR + 1
               WRITE (LUERR,'(/A,I5)')
     *            ' *** ERROR Sirius GENINP, NPATH outside range:',NPATH
            ELSE IF (NPATH(I).LT.0) THEN
               FLAG(-NPATH(I)) = .FALSE.
            ELSE IF (NPATH(I).GT.0) THEN
               FLAG(NPATH(I)) = .TRUE.
            ELSE
               GO TO 299
            END IF
  286    CONTINUE
  299    CONTINUE
      GO TO 100
C
C *** Option 3 >MCSCF <  MCSCF calculation
  103 CONTINUE
         DOMC = .TRUE.
         FLAG(2) = .TRUE.
         NUMRUN = NUMRUN + 1
      GO TO 100
C
C *** Option 4 >CC< We want to do a subsequent Coupled Cluster calculation
  104 CONTINUE
         DOCCSD = .TRUE.
         NUMRUN = NUMRUN + 1
         IF (.NOT. DORHF) THEN
            DORHF  = .TRUE.
            NUMRUN = NUMRUN + 1
         END IF
         FLAG(21) = .TRUE.
      GO TO 100
C
C *** Option 5 >RESTAR< restart SIRIUS, LUIT1 must be available
  105 CONTINUE
        IREST = 1
      GO TO 100
C
C *** Option 6 >HARTRE<  Restricted, closed shell Hartree-Fock SCF
  106 CONTINUE
         IF (.NOT.DORHF) THEN
C           ... DORHF may also be set with .CC
            DORHF    = .TRUE.
            FLAG(21) = .TRUE.
C           ... flag for QCHF -- quadratically convergent Hartree-Fock
            NUMRUN = NUMRUN + 1
         END IF
      GO TO 100
C
C *** Option 7 >STOP  <  Stop acc. to instruction given on next line
  107 CONTINUE
         READ (LUCMD,'(A20)') REWORD
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
C *** Option  10 >MP2  <  Moller-Pleset PT2 calculation
  110 CONTINUE
        DOMP2    = .TRUE.
        FLAG(9)  = .TRUE.
        NUMRUN   = NUMRUN + 1
        FLAG(3)  = .FALSE.
      GO TO 100
C
C *** Option 11 >NSYM  < Number of symmetries.
C    Specified by integral program, used for consistency check.
  111 CONTINUE
         READ (LUCMD,*) MSYM
         IF (MSYM .NE. NSYM) THEN
            INPERR = INPERR + 1
            WRITE (LUPRI,'(/A,2I5/A)')
     *   ' >>> GENINP-ERROR, NSYM from LUONEL .ne. NSYM from input:',
     *      NSYM,MSYM,
     &   '     NSYM from input used in reading remaining input.'
            NSYM = MSYM
         END IF
      GO TO 100
C
C *** Option 12 >CC ONL<  Write information to ABACUS on unit LUSIFC
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
C
C *** Option 6 >HARTRE<  Restricted, closed shell Hartree-Fock SCF
  116 CONTINUE
         IF (.NOT.DORHF) THEN
C           ... DORHF may also be set with .CC
            DODFT    = .TRUE.
            DORHF    = .TRUE.
            FLAG(21) = .TRUE.
C           ... flag for QCHF -- quadratically convergent Hartree-Fock
            NUMRUN = NUMRUN + 1
         END IF
      GO TO 100
C
C
      END
C  /* Deck haminp */
      SUBROUTINE HAMINP(WORD,INPERR,ALLOPT)
C
C Copyright Nov 1990 Hans Joergen Aa. Jensen
C
#include <implicit.h>
      PARAMETER (NTABLE = 2)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7,
     *          REWORD*12
      LOGICAL ALLOPT
C
C Used from common:
C
C INFINP : NFIELD,EFIELD,LFIELD
C INFPRI : P6FLAG
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infpri.h>
C
      DATA TABLE /'.FIELD ','.PRINT '/
C
C     ***** PROCESS INPUT FOR HAMINP *****
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
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
         CALL QUIT(' ILLEGAL KEYWORD IN HAMINP ')
C
C *** Option 1 >FIELD TERM<   Field dependent Hamiltonian term
  101 CONTINUE
         NFIELD = NFIELD + 1
         IF (NFIELD .LE. MXFELT) THEN
            READ(LUCMD,*) EFIELD(NFIELD)
            READ(LUCMD,'(1X,A8)') LFIELD(NFIELD)
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
      END
C  /* Deck wavinp */
      SUBROUTINE WAVINP(WORD,INPERR,ALLOPT)
C
C 1-Oct 1987 Hans Agren
C
#include <implicit.h>
      PARAMETER (NTABLE = 12)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7,
     *          REWORD*12,  TITLIN*80
      LOGICAL ALLOPT
C
C Used from common:
C
C INFORB : NASH, NISH, NSYM,...
C INFINP : NACTEL,ISPIN,FLAG(*),MCTYP,..
C
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infinp.h>
C
      LOGICAL RAS1EL
C
      DATA TABLE /'.CAS SP','.INACTI','.ELECTR','xxxxxxx','.SPIN M',
     *            '.SYMMET','.RAS1 H','.RAS1 S','.RAS2 S','.RAS3 S',
     *            '.RAS1 E','.RAS3 E'/
C
      NHL1MN = 0
      NHL1MX = -1
      RAS1EL = .FALSE.
C
C     ***** PROCESS INPUT FOR WAVINP *****
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
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
         CALL QUIT(' ILLEGAL KEYWORD IN WAVINP ')
C
C *** Option 1 >CAS  S<  CAS calculation: Active orbitals
  101 CONTINUE
      IF (MCTYPE .EQ. 2) THEN
         WRITE (LUPRI,*) 'Sirius WAVINP ERROR, both CAS and '//
     &                   'RAS specified.'
         INPERR = INPERR + 1
      END IF
         MCTYPE = 1
         READ(LUCMD,*) (NAS2(I),I=1,NSYM)
         DO 25 ISYM = 1,NSYM
            NAS1(ISYM) = 0
            NAS3(ISYM) = 0
   25    CONTINUE
      GO TO 100
C
C *** Option 2 >INACTI<  Inactive orbitals
  102 CONTINUE
         READ(LUCMD,*) (NISH(I),I=1,NSYM)
         NRHFT = ISUM(NSYM,NRHF,1)
         IF (NRHFT .EQ. 0) THEN
            DO 402 I = 1, NSYM
               NRHF(I) = NISH(I)
  402       CONTINUE
         END IF
      GO TO 100
C
C *** Option 3 >ELECTR< Number of active electrons
  103 READ(LUCMD,*) NACTEL
      GO TO 100
C
C *** Option 4 >xxxxxx<
  104 CONTINUE
      GO TO 100
C
C *** Option 5 >SPIN M<  state spin multiplicity
C                         SPIN QUANTUM NUMBER*2+1
  105 READ(LUCMD,*) ISPIN
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
     *         ' WAVINP error, both .RAS1 H and .RAS1 E specified.'
            INPERR = INPERR + 1
         END IF
         READ(LUCMD,*) NHL1MN,NHL1MX
         RAS1EL = .TRUE.
      GO TO 100
C
C *** Option 8 >RAS1 S<   RAS1 orbital space
  108 CONTINUE
      IF (MCTYPE .EQ. 1) THEN
         WRITE (LUPRI,*) 'Sirius WAVINP ERROR, both CAS '//
     &                   'and RAS specified.'
         INPERR = INPERR + 1
      END IF
         MCTYPE = 2
         READ(LUCMD,*) (NAS1(I),I=1,NSYM)
      GO TO 100
C
C *** Option 9 >RAS2 S<   RAS2 orbital space
  109 CONTINUE
      IF (MCTYPE .EQ. 1) THEN
         WRITE (LUPRI,*) 'Sirius WAVINP ERROR, both CAS '//
     &                   'and RAS specified.'
         INPERR = INPERR + 1
      END IF
         MCTYPE = 2
         READ(LUCMD,*) (NAS2(I),I=1,NSYM)
      GO TO 100
C
C *** Option 10 >RAS3 S<   RAS3 orbital space
  110 CONTINUE
      IF (MCTYPE .EQ. 1) THEN
         WRITE (LUPRI,*) 'Sirius WAVINP ERROR, both CAS '//
     &                   'and RAS specified.'
         INPERR = INPERR + 1
      END IF
         MCTYPE = 2
         READ(LUCMD,*) (NAS3(I),I=1,NSYM)
      GO TO 100
C
C *** Option 11 >RAS1 E<   RAS1 electrons
  111 CONTINUE
         IF (RAS1EL) THEN
            WRITE (LUPRI,'(/A/)')
     *         ' WAVINP error, both .RAS1 H and .RAS1 E specified.'
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
      DO 2100 I = 1,NSYM
         NASH(I) = NAS1(I) + NAS2(I) + NAS3(I)
         NASHT   = NASHT   + NASH(I)
 2100 CONTINUE
      IF (NHL1MX .GT. NHL1MN) THEN
         NAS1T = 0
         DO 2200 I = 1,NSYM
            NAS1T = NAS1T + NAS1(I)
 2200    CONTINUE
         NELMN1 = 2*NAS1T - NHL1MX
         NELMX1 = 2*NAS1T - NHL1MN
      END IF
      RETURN
      END
C  /* Deck rhfinp */
      SUBROUTINE RHFINP(WORD,INPERR,ALLOPT)
C
C 8-May-1989 Hans Joergen Aa. Jensen
C 13-May-1995 K.Ruud, added keyword .ELECTRONS.
C 1999 Helgaker, DFT input added 
C
#include <implicit.h>
      PARAMETER (NTABLE = 38)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7,
     *          REWORD*12
      LOGICAL ALLOPT
C
C Used from common:
C
C SCBRHF : NFRRHF(8), IOPRHF, RHFCAN, MXHFMA, MXHFMI
C INFINP : THRRHF,FLAG(*),DIRFCK,MXDIIS
C INFORB : NRHF(8)
C INFTAP : LUSUPM
C INFPRI : IPRSIR,IPRRHF
C
#include <maxorb.h>
#include <priunit.h>
#include <scbrhf.h>
#include <infinp.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <dftcom.h>
#include <dfterg.h>
C
      LOGICAL NODIIS,SET,ELECTR
      SAVE SET,ELECTR
      DATA SET,ELECTR /2*.FALSE./
C
      DATA TABLE /'.HF OCC','.THRESH','.FOCK I','.H1VIRT','.NOQCHF',
     *            '.PRINT ','.OPEN S','.MAX MA','.MAX MI','.NONCAN',
     *            '.FROZEN','.DIRFOC','.NODIIS','.MAX DI','.C2DIIS',
     *            '.MAX ER','.DDSCF ','.COREHO','.CORERE','XXXXXXX',
     *            '.ELECTR','.AUTOCC','.FC MVO','.HFXFAC','.LDA   ',
     *            '.BLYP  ','.B3LYP ','.FROMVX','.WAH   ','.DFTWEI',
     *            '.DFTAC ','.DFTHES','.DFTTHR','.DFTELS','.RADINT',
     *            '.ANGINT','.ANGMIN','.NOPRUN'/
C
C     ***** PROCESS INPUT FOR RHFINP *****
C
      IF (SET .AND. WORD .EQ. '*HF INP') THEN
         INPERR = INPERR + 1
         WRITE(LUPRI,*)'Input error: *HF INP specified more than once'
      END IF
      NODIIS = .FALSE.
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102,103,104,105,106,107,108,109,110,
     &                   111,112,113,114,115,116,117,118,119,120,
     &                   121,122,123,124,125,126,127,128,129,130,
     &                   131,132,133,134,135,136,137,138),II
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
         CALL QUIT(' ILLEGAL KEYWORD IN RHFINP ')
C
C *** Option 1 >HF OCC<  Hartree-Fock occupancy for RHF and MP2 calc.
  101 READ(LUCMD,*) (NRHF(I),I=1,NSYM)
      IF (INIOCC .LT. 0) AUTOCC = .FALSE.
      IF (ELECTR) THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,*)
     &' .ELECTR must not be specified together with .HF OCC and .OPEN S'
      END IF
      NRHFEL = 0
C     ... AUTOCC may have been set true by .AUTOCC
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
     &'RHFINP INFO: more than one specification of .H1VIRT and .FC MVO'
         WRITE (LUPRI,*) '.H1VIRT will be used.'
      END IF
      FCVORB = .TRUE.
      GO TO 100
C
C *** Option 5 >NOQCHF<  No quadratically convergent Hartree-Fock
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
     &' .ELECTR must not be specified together with .HF OCC and .OPEN S'
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
C *** Option 11 >FROZEN<  Frozen orbitals per symmetry
  111 READ (LUCMD,*) (NFRRHF(I),I=1,NSYM)
      GO TO 100
C
C *** Option 12 >DIRFOC< Direct Fock matrix construction
  117 WRITE (LUPRI,*) 'WARNING: .DDSCF has been replaced by .DIRFOCK'
      WRITE (LUPRI,*) 'WARNING: .DDSCF may be removed in future relase'
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
C 117: see 112
C
C
C *** Option 18 >COREHO< core hole calculation
  118 READ (LUCMD,*) JCHSYM,JCHORB
      CORHOL = .TRUE.
      GO TO 100
C
C *** Option 19 >CORERE< core hole relaxation
  119 CORRLX = .TRUE.
      GO TO 100
C *** Option 20 >XXXXXX<
  120 CONTINUE
      GO TO 100
C *** Option 21 >ELECTR< No. of electrons, for automagic det. of HF occ.
  121 READ(LUCMD,*) NRHFEL
      IF (INIOCC .EQ. 1) THEN
         INPERR = INPERR + 1
         WRITE (LUPRI,*)
     &' .ELECTR must not be specified together with .HF OCC and .OPEN S'
      END IF
      AUTOCC = .TRUE.
      ELECTR = .TRUE.
      IF (INIOCC .LT. 0) INIOCC = 0
      GO TO 100
C *** Option 22 >AUTOCC< automatic determination of  HF occ. in DIIS,
C                        with given inital occ. from user input or file
  122 AUTOCC = .TRUE.
      IF (INIOCC .LT. 0) INIOCC = 0
      GO TO 100
C
C *** Option 23 >FC MVO<  Modified virtual orbitals from FC from NMVO
C     Reference: C.W. Bauschlicher, JCP 72 (1980) 880.
  123 READ (LUCMD,*) (NMVO(I), I = 1,NSYM)
      NMVOT = ISUM(8,NMVO,1)
      IF (FCVORB) THEN
         NINFO = NINFO + 1
         WRITE (LUPRI,*)
     &'RHFINP INFO: more than one specification of .H1VIRT and .FC MVO'
         WRITE (LUPRI,*) 'Last .FC MVO will be used.'
      END IF
      FCVORB = .TRUE.
      GO TO 100
C
C     DFT input
C
  124 CONTINUE
      READ(LUCMD,*) HFXFAC
      DFTADD = .TRUE.
      DFTRUN = .TRUE.
      GO TO 100
C
  125 CONTINUE
      DFTTYP = 'LDA   '
      DFTADD = .TRUE.
      DFTRUN = .TRUE.
      WDFTX  = D1
      WDFTB  = D0
      WDFTC  = D1
      WDFTL  = D0
      HFXFAC = D0
      GO TO 100
  126 CONTINUE
      DFTTYP = 'BLYP  '
      DFTADD = .TRUE.
      DFTRUN = .TRUE.
      WDFTX  = D1
      WDFTB  = D1
      WDFTC  = D0
      WDFTL  = D1
      HFXFAC = D0
      GO TO 100
  127 CONTINUE
      DFTTYP = 'B3LYP '
      DFTADD = .TRUE.
      DFTRUN = .TRUE.
      WDFTX  = 0.80D0
      WDFTB  = 0.72D0
      WDFTL  = 0.81D0
      WDFTC  = D1 - WDFTL
      HFXFAC = D1 - WDFTX
      GO TO 100
  128 CONTINUE
      DFTPOT = .TRUE.
      GO TO 100
  129 CONTINUE
      DFTTYP = 'WAH   '
      DFTORD = .TRUE.
      DFTADD = .TRUE.
      DFTRUN = .TRUE.
      WDFTX  = 0.95D0
      WDFTB  = 0.72D0
      WDFTL  = 0.81D0
      WDFTC  = D1 - WDFTL
      HFXFAC = D1 - WDFTX
      GO TO 100
  130 CONTINUE
      DFTTYP = 'OWN   '
      DFTADD = .TRUE.
      DFTRUN = .TRUE.
      WDFTX  = 0.95D0
      WDFTB  = 0.72D0
      WDFTL  = 0.81D0
      READ(LUCMD,*) WDFTX
      WDFTC  = D1 - WDFTL
      HFXFAC = D1 - WDFTX
      GO TO 100
 131  CONTINUE
      DFTASC = .TRUE.
      DFTPOT = .TRUE.
      READ (LUCMD,*) DFTIPT, DFTBR1, DFTBR2
      GO TO 100
 132  CONTINUE
      DFTHES = .TRUE.
      GO TO 100
 133  CONTINUE
      DFTHRS = .TRUE.
      READ (LUCMD,*) DFTHR0, DFTHRL, DFTHRI
      GO TO 100
 134  CONTINUE
      READ (LUCMD,*) DFTELS 
      GOTO 100
C&&& RADINT: Give relative precision of the radial integration
 135  CONTINUE
      READ (LUCMD,*) RADINT
      GO TO 100
C&&& ANGINT: Give precision of the angular integration.
C            The integration is exact for spherical harmonics up to
C            order ANGINT
 136  CONTINUE
      READ(LUCMD,*) ANGINT
      GO TO 100
C&&& ANGMIN: Set minimum order of Lebedev scheme after pruning
 137  CONTINUE
      READ(LUCMD,*) ANGMIN
      GO TO 100
C&&& NOPRUN: Turn off pruning of angular grid.
 138  CONTINUE
      NOPRUN = .TRUE.
      GO TO 100
C
 2000 CONTINUE
      IF (NODIIS) MXDIIS = 0
      IF (CORHOL) THEN
         IF (IOPRHF .NE. 0) THEN
            INPERR = INPERR + 1
            WRITE (LUPRI,*) 'RHFINP input error:'//
     &         ' .OPEN SHELL and .COREHOLE not compatible'
         END IF
         IOPRHF = JCHSYM
         NRHF(JCHSYM) = NRHF(JCHSYM) - 1
         IF (NRHF(JCHSYM) .LT. 0) THEN
            INPERR = INPERR + 1
            WRITE (LUPRI,*) 'RHFINP input error:'//
     &         ' no HF OCCUPATION in core hole symmetry'
         END IF
         MXDIIS = 0
         MAXFCK = 0
      END IF
      SET = .TRUE.
      RETURN
      END
C  /* Deck mp2inp */
      SUBROUTINE MP2INP(WORD,INPERR,ALLOPT)
C
C 17/9 - 1990 Hans Agren
C
#include <implicit.h>
      PARAMETER (NTABLE = 2)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7,
     *          REWORD*12
      LOGICAL ALLOPT
C
C Used from common:
C
C INFORB : NSYM
C INFMP2 : NFRMP2(*)
C INFPRI : IPRMP2
C
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infmp2.h>
#include <infpri.h>
C
      DATA TABLE /'.MP2 FR','.PRINT '/
CMAERKE 900918: MP2SAV not implemented, should it be?
CMAERKE 900918: always NEWORB in this version.
C
C     ***** PROCESS INPUT FOR MP2INP *****
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
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
         CALL QUIT(' ILLEGAL KEYWORD IN MP2INP ')
C
C *** Option 1 >MP2 FROZEN<   Frozen MP2 orbitals
  101 CONTINUE
         READ(LUCMD,*) (NFRMP2(I),I=1,NSYM)
      GO TO 100
C
C *** Option 2 >PRINT <  General print level in MP2
  102 READ(LUCMD,*) IPRMP2
      GO TO 100
C
      END
C  /* Deck optinp */
      SUBROUTINE OPTINP(WORD,INPERR,ALLOPT)
C
C 1-Oct 1987 Hans Agren
C
#include <implicit.h>
      PARAMETER (NTABLE = 33)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7,
     *          REWORD*12
      LOGICAL ALLOPT
C
C Used from common:
C
C INFINP : FLAG(*),ICHECK,*ROOT,ISTATE,MAX*,THR*
C INFVAR : NWOPT, JWOP(*)
C INFDIM : MAXPHP
C INFOPT : RAT*
C INFPRI : LUERR
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infvar.h>
#include <infdim.h>
#include <inftap.h>
#include <infopt.h>
C
      DATA TABLE /'.ABSORP','.ALWAYS','.ACTROT','xXXXXXX','.MAXAPM',
     *            '.SYM CH','.THRCGR','.STATE ','.COREHO','.CORERE',
     *            '.MAX MA','.MAXABS','.MAX CI','.EXACTD','.MAX MI',
     *            '.NEO AL','.PHPRES','.NO ABS','.NR ALW','.OPTIMA',
     *            '.ORB_TR','xXXXXXX','.SIMULT','.THRESH','.OLSEN ',
     *            '.TRACI ','.FOCKDI','.NOTRAC','.CI PHP','.NO ACT',
     *            '.DETERM','.NATONL','.FOCKON'/
C
C Internal flags
      IFNSIM = 0
C
C     ***** PROCESS INPUT FOR OPTINP *****
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
 1001    PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102,103,104,105,106,107,108,109,110,
     *                   111,112,113,114,115,116,117,118,119,120,
     *                   121,122,123,124,125,126,127,128,129,130,
     *                   131,132,133), II
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
         CALL QUIT(' ILLEGAL KEYWORD IN OPTINP ')
C
C
C *** Option 1 >ABSORP<  Orbital absorption
 101  CONTINUE
         READ (LUCMD,'(A12)') REWORD
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
C *** Option 4 >XXXXXX<  Call SIRCAN for canonical/Fock analysis
  104 CONTINUE
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
         READ (LUCMD,*) ISTATE
C
         IF ( ISTATE .GT. MAXRTS ) THEN
            WRITE (LUPRI,*) 'OPTINP error: .STATE input is greater than'
            WRITE (LUPRI,*) 'internal maximum of',MAXRTS
            CALL QUIT(
     &         'OPTINP error: .STATE input exceeds internal limit')
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
C *** Option 22 >XXXXXX< Call SIROUT for optimization result summary
  122 CONTINUE
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
  128 IMCCNO = -1
      GO TO 100
C
C *** Option 29 >CI PHP<  Explicit CI Hamiltonian matrix set up
  129 CONTINUE
        READ (LUCMD,*) MAXPHP
      GO TO 100
C
C *** Option 30 > NO ACT<  No active-active rotations in RAS
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
 9000 CONTINUE
      NROOTS = MAX(ISTATE,NROOTS)
      IF ( ISTATE .GT. MAXRTS ) THEN
         INPERR = INPERR + 1
         WRITE (LUERR,9008) ISTATE, MAXRTS
      ELSE IF ( NROOTS .GT. MAXRTS ) THEN
         NWARN = NWARN + 1
         WRITE (LUERR,9009) NROOTS, MAXRTS
         NROOTS = MAXRTS
      END IF
      LROOTS = MAX(LROOTS,NROOTS)
      IROOT(1) = ISTATE
      J = 1
      DO 9100 I = 1,NROOTS
         IF (I .NE. ISTATE) THEN
            J = J + 1
            IROOT(J) = I
         END IF
 9100 CONTINUE
 9008 FORMAT(/' Input error: .STATE value of,',I3,
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
#include <implicit.h>
      PARAMETER (NTABLE = 11)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7,
     &          REWORD*12, RWORD*7
      LOGICAL ALLOPT
C
C Used from common:
C
C INFINP : DOCINO, FLAG(*), NROOCI, THRCI, THRPWF
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
C
      DATA TABLE /'.CINO  ','.CIDENS','.CIROOT','.DISKH2','.THRESH',
     *            '.THRPWF','.WEIGHT','.ZEROEL','.STATE ','XXXXXXX',
     *            '.MAX IT'/
C
C     ***** PROCESS INPUT FOR CICINP *****
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102,103,104,105,106,107,108,109,110,
     *                   111) II
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
         CALL QUIT(' ILLEGAL KEYWORD IN CICINP ')
C
C *** Option 1 >CINO  < Calculate CI natural orbitals.
  101 CONTINUE
      DOCINO   = .TRUE.
      FLAG(68) = .TRUE.
      FLAG(67) = .TRUE.
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
C *** Option 10 >XXXXXX<
  110 CONTINUE
      GO TO 100
C
C *** Option 11 >MAX IT< Maximum number of iterations in
C                        CI diagonalization.
  111 CONTINUE
      READ (LUCMD,*) MXCIMA
      GO TO 100
C
      END
C  /* Deck orbinp */
      SUBROUTINE ORBINP(WORD,INPERR,IRDMO,ALLOPT)
C
C 1-Oct 1987 Hans Agren
C
#include <implicit.h>
#include <priunit.h>
C
C -- local variables and constants
C
      PARAMETER (NTABLE = 16)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7,
     &          REWORD*12, RWORD*7
      PARAMETER (MAXNRT = 40)
      DIMENSION NOR(8),INORO(MAXNRT),NREOR(8),NNOR(8)
      LOGICAL ALLOPT
C
C Used from common:
C
C INFINP : FLAG,IMOORD,CMAXMO,SUPSYM,THRSSY
C INFORB : NFRO, NBAS
C HUCKEL : ADDSTO
C
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <huckel.h>
C
      DATA TABLE /'XXXXXXX','.5D7F9G','.DELETE','.MOSTAR','.FREEZE',
     *            '.GRAM-S','.SYMMET','.PUNCHI','.PUNCHO','.REORDE',
     *            '.FROZEN','.AO DEL','.CMOMAX','.NOSUPS','.SUPSYM',
     *            '.THRSSY'/
C
      INDEL = 0
C
C     ***** PROCESS INPUT FOR ORBINP *****
C
      WORD1 = WORD
      IF(ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                      LUPRI)
  100    READ (LUCMD, '(A7)') WORD
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
     &           CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                      LUPRI)
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
         CALL QUIT(' ILLEGAL KEYWORD IN ORBINP ')
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
C            CMODEL will be called later to get CMO, if MOSTRT 8.
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
        READ (LUCMD,'(A12)') REWORD
Chj-aug99: accept both starting in column 2 (used to be required)
C          and in column 1 (as some other places).
          IF (REWORD(1:1) .EQ. ' ' .OR. REWORD(1:1) .EQ. '.') THEN
             RWORD = REWORD(2:7)
          ELSE
             RWORD = REWORD(1:6)
          END IF
          IF (RWORD .EQ. 'FORM18') THEN
             IRDMO = 3
          ELSE IF (RWORD .EQ. 'FORM12') THEN
             IRDMO = 4
          ELSE IF (RWORD .EQ. 'CANORB') THEN
             IRDMO = 1
          ELSE IF (RWORD .EQ. 'MOLORB') THEN
             IRDMO = 2
          ELSE IF (RWORD .EQ. 'OLDORB') THEN
             IRDMO = 7
          ELSE IF (RWORD .EQ. 'H1DIAG') THEN
             IRDMO = 8
             ADDSTO = .FALSE.
C         ... make sure Huckel guess is disabled,
C             although code in GNRLIN in dalton.F takes care of that
          ELSE IF (RWORD .EQ. 'NEWORB') THEN
             IRDMO = 9
          ELSE IF (RWORD .EQ. 'SIRIFC') THEN
             IRDMO = 10
          ELSE IF (RWORD .EQ. 'HUCKEL') THEN
             IRDMO = 8
             IF (.NOT.ADDSTO) THEN
                WRITE (LUPRI,'(/A)') ' INFO: .MOSTART HUCKEL changed'//
     &             ' to H1DIAG because Huckel guess not available.'
             END IF
          ELSE
             WRITE(LUPRI,'(3A/A)') '"',REWORD(1:7),
     *          '" is incorrect input after ".MOSTART"',
     *          ' Select from .FORM12, .FORM18, .H1DIAG, .HUCKEL,',
     *          ' .SIRIFC, .OLDORB, and .NEWORB'
             CALL QUIT(' INPUT ERROR IN ORBINP')
          END IF
          IF (IRDMO .EQ. 8) THEN
             FLAG(13) = .FALSE.
          ELSE
             FLAG(13) = .TRUE.
          END IF
          IF (IRDMO .NE. 8) ADDSTO = .FALSE.
C         ... make sure Huckel guess is disabled,
C             although code in GNRLIN in dalton.F takes care of that
      GO TO 100
C
C *** Option 5 >FREEZE< Orbitals not to be rotated.
  105 CONTINUE
      IF (FLAG(33)) THEN
         CALL QUIT(
     &      'ORBINP error, .FREEZE and .AO DELETE both specified.')
      END IF
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
 2311 FORMAT(/' Sirius ORBINP-ERROR,',I3,' orbitals specified not to ',
     *    'be rotated in symmetry',I2/20X,'Maximum is',I3)
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
     &      'ORBINP error, .REORDE and .AO DELETE both specified.')
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
               CALL QUIT('Error in ORBINP, internal limit exceeded')
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
C *** Option 11 >FROZEN< frozen core orbitals
  111 CONTINUE
      READ (LUCMD,*) (NFRO(I),I=1,NSYM)
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
C *** Option 14 >NOSUPS< do not use super symmetry
  114 SUPSYM = .FALSE.
      GO TO 100
C
C *** Option 15 >SUPSYM< do use super symmetry
  115 SUPSYM = .TRUE.
      GO TO 100
C
C *** Option 16 >THRSSY< read super symmetry threshold
  116 READ(LUCMD,*) THRSSY
      GO TO 100
C
      END
C  /* Deck civinp */
      SUBROUTINE CIVINP(WORD,INPERR,ALLOPT)
C
C 1-Oct 1987 Hans Agren
C
#include <implicit.h>
      PARAMETER (NTABLE = 5)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7
      LOGICAL ALLOPT
C
C Used from common:
C
C INFINP : FLAG(*), ICI0,
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
C
      DATA TABLE /'.PLUS C','.SELECT','.STARTH','.ABACUS','.STARTO'/
C
C     ***** PROCESS INPUT FOR CIVINP *****
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
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
         CALL QUIT(' ILLEGAL KEYWORD IN CIVINP ')
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
#include <implicit.h>
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
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infopt.h>
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
         CALL QUIT(' ILLEGAL KEYWORD IN STPINP ')
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
C *** Option 7 >REJ RA<  Ratio predicted/actual step that rejects step
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
#include <implicit.h>
C
C -- local variables and constants
C
      PARAMETER (NTABLE = 10)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7
      LOGICAL ALLOPT
C
C Used from common:
C
C INFINP : FLAG(*)
C INFPRI : IPRMUL
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infpri.h>
C
      DATA TABLE /'.ALL   ','.DIPMOM','.GROSSA','.GROSSM','.MULLIK',
     *            '.NETALL','.NETMO ','.QUADRP','.VIRIAL','.PRINT '/
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
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102,103,104,105,106,107,108,109,110), II
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
         CALL QUIT(' ILLEGAL KEYWORD IN POPINP ')
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
      GO TO 100
C
C *** Option 6 >NETALL<  All net population analysis
  106 FLAG(73) = .FALSE.
      IPRMUL = MAX(IPRMUL,5)
      GO TO 100
C
C *** Option 7 >NETMO <  Net MO population analysis
  107 FLAG(73) = .FALSE.
      IPRMUL = MAX(IPRMUL,2)
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
      GO TO 100
C
      END
C  /* Deck trainp */
      SUBROUTINE TRAINP(WORD,INPERR,ALLOPT)
C
C 1-Oct 1987 Hans Agren
C
#include <implicit.h>
C
C -- local variables and constants
C
      PARAMETER (NTABLE = 8)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7
      LOGICAL ALLOPT
C
C Used from common:
C
C INFINP : FLAG(*)
C INFDIM : MWORK
C INFTRA : ITRLVL, THRP, THRQ
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infdim.h>
#include <inftra.h>
C
      DATA TABLE /'.RESIDE','.FINAL ','.LEVEL ','.OLD TR','.PRINT ',
     *            '.P THRE','.Q THRE','.PRESOR'/
C
C     ***** PROCESS INPUT FOR TRAINP *****
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
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
            RETURN
         ELSE
            WRITE (LUPRI,'(/3A/2A/)')
     *         ' Keyword ',WORD,' does not begin with',
     *         ' one of the four characters ".*!#" for ',WORD1
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT(' ILLEGAL KEYWORD IN TRAINP ')
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
C *** Option  107 >Q THRE<
  107 READ(LUCMD,*) THRQ
      GO TO 100
C
C *** Option  107 >SORTED<
  108 CONTINUE
      NEWTRA=.TRUE.
      GO TO 100
C
      END
C  /* Deck updinp */
      SUBROUTINE UPDINP(WORD,INPERR,ALLOPT)
C
C 1-Oct 1987 Hans Agren
C
#include <implicit.h>
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
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
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
         CALL QUIT(' ILLEGAL KEYWORD IN UPDINP ')
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
#include <implicit.h>
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
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infpri.h>
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
         CALL QUIT(' ILLEGAL KEYWORD IN SOLINP ')
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
         WRITE (LUERR,'(/1X,A/A)') WORD1//' input error',
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
#include <implicit.h>
C
C -- local variables and constants
C
      PARAMETER (NTABLE = 16)
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
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infpri.h>
#include <inftra.h>
C
      DATA TABLE /'XXXXXXX','.IPRKAP','.IPRDIA','.IPRSIG','.IPRDNS',
     *            '.IPRSOL','.RESFIL','.PRINTL','.PRINTF','.IPRCIX',
     *            '.THRPWF','.IPRCNO','.NOSUMM','.CANONI','.IPRAVE',
     *            '.IPRFCK'/
C
C     ***** PROCESS INPUT FOR PRIINP *****
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
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
         CALL QUIT(' ILLEGAL KEYWORD IN PRIINP ')
C
C *** Option 1 >      < 
C
  101 CONTINUE
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
C           READ (LUCMD,*) (NPP6TH(I), I=1,NUM6)
            DO 298 I = 1,NUM6
               IF (IABS(NPPATH(I)) .GT. NPFLAG) THEN
                  INPERR = INPERR + 1
                  WRITE (LUERR,'(/A,I5)')
     *               ' ERROR PRIINP, NP6PTH outside range',NPPATH(I)
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
C           READ (LUCMD,*) (NP4PTH(I), I=1,NUM4)
            DO 288 I = 1,NUM4
               IF (IABS(NPPATH(I)) .GT. NPFLAG) THEN
                  INPERR = INPERR + 1
                  WRITE (LUERR,'(/A,I5)')
     *               ' ERROR PRIINP, NP4PTH outside range',NPPATH(I)
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
      END
C  /* Deck auxinp */
      SUBROUTINE AUXINP(WORD,INPERR,ALLOPT)
C
C 1-Oct 1987 Hans Agren
C
#include <implicit.h>
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
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inftap.h>
C
      DATA TABLE /'XXXXXXX','.NOSUPM'/
C
C     ***** PROCESS INPUT FOR AUXINP *****
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
  100    READ (LUCMD, '(A7)') WORD
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
         CALL QUIT(' ILLEGAL KEYWORD IN AUXINP ')
C
C *** Option 1 >ONESUP<  use same units for SUPM. and ONEL. integrals
C              (obsolete, removed March 2000 /Kenneth Ruud)
  101 GO TO 100
C
C *** Option 2 >NOSUPM<  Do not use SUPMAT
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
#include <implicit.h>
C
      CHARACTER TAB1*15,TAB2*30
      PARAMETER (TAB1=' ERROR MESSAGE:')
      PARAMETER (TAB2=' IS NOT WITHIN GIVEN LIMITS OF')
C
C
C Used from common blocks:
C   INFVAR : NWOPT
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infvar.h>
#include <inforb.h>
#include <infdim.h>
#include <inftap.h>
#include <infpri.h>
C
C
      INPERR = INPERR + 1
      GO TO ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     *       11,12,13,14,15,16,17),ITYP
C
    1 CONTINUE
      RETURN
C
    2 WRITE(LUERR,200) TAB1,NACTEL,TAB2
  200 FORMAT(/A/' NUMBER OF ACTIVE ELECTRONS=',I3,A,' OF 0 TO 35')
      RETURN
C
    3 S = ISPIN-1
      S = S / 2
      WRITE(LUERR,300) TAB1,S,TAB2
  300 FORMAT(/A/' SPIN QUANTUM NUMBER=',F5.1,A,' 0 TO 10')
      RETURN
C
    4 WRITE(LUERR,400) TAB1,NSYM,TAB2
  400 FORMAT(/A/' NUMBER OF SYMMETRIES=',I3,A,' 1 TO 8')
      RETURN
C
    5 WRITE(LUERR,500) TAB1,MAXMAC,TAB2
  500 FORMAT(/A/' MAXIMUM NUMBER OF ITERATIONS=',I3,A,' 1 TO 50')
      RETURN
C
    6 CONTINUE
      RETURN
C
    7 WRITE(LUERR,700) TAB1,THRGRD,TAB2
  700 FORMAT(/A/' THRESHOLD FOR GRADIENT CONVERGENCE=',1P,E10.2,A,
     *       ' 0.0 TO 1.0')
      RETURN
C
    8 WRITE(LUERR,800) TAB1,NORBT,TAB2,MXCORB
  800 FORMAT(/A/' TOTAL NUMBER OF ORBITALS=',I4,A,' 0 TO',I4)
      RETURN
C
    9 WRITE(LUERR,900) TAB1,NRHFT,TAB2,MAXOCC
  900 FORMAT(/A/' NUMBER OF RHF ORBITALS=',I4,A,' 0 TO',I4)
      RETURN
C
   10 WRITE(LUERR,1000) TAB1,NASHT,TAB2,MAXASH
 1000 FORMAT(/A/' NUMBER OF ACTIVE ORBITALS=',I4,A,' 0 TO',I4)
      RETURN
C
   11 WRITE(LUERR,1100) TAB1,NSSHT,TAB2,MXCORB
 1100 FORMAT(/A/' NUMBER OF SECONDARY ORBITALS=',I4,A,' 0 TO',I4)
      RETURN
C
   12 WRITE(LUERR,1200) TAB1,(IMOORD(I),I=1,NORBT)
 1200 FORMAT(/A/' REORDERING INDICES DO NOT GIVE CORRECT SUM CHECK',
     *       /(23X,10I4))
      RETURN
C
   13 WRITE(LUERR,1300) TAB1,NROOTS,TAB2,MAXRTS
 1300 FORMAT(/A/' NUMBER OF NEO ROOTS=',I3,A,' 1 TO',I4)
      RETURN
C
   14 WRITE(LUERR,1400) TAB1,(IROOT(I),I=1,NROOTS)
 1400 FORMAT(/A/' NEO ROOT INDICES DUPLICATE OR NEGATIVE'/(23X,10I3))
      RETURN
C
   15 RETURN
C
   16 WRITE(LUERR,1600) TAB1,NWOPT,TAB2,MAXWOP
 1600 FORMAT(/A/' NUMBER OF ORBITAL ROTATIONS =',I8,A,' 0 TO',I5)
      RETURN
C
   17 RETURN
C
      END
C  /* Deck prtinp */
      SUBROUTINE PRTINP(IWUNIT,NUMRUN,IRDMO)
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
C  SUBROUTINES CALLED:NONE
C
C  (Based on CASSCF RELEASE 79 01 10)
C
#include <implicit.h>
C
C Used from common blocks:
C   INFINP : DOMC,DOCI,DIRFCK,FLAG(21),MAXFCK,MCTYPE,THR*,NFIELD,...,?
C   INFVAR : NCONF,NWOPT,NVAR
C   INFIND : ISMO(),?
C   INFDIM : MAXPHP
C   SCBRHF : IOPRHF,NFRRHF(*)
C   INFMP2 : NFRMP2(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infvar.h>
#include <infind.h>
#include <inforb.h>
#include <infdim.h>
#include <scbrhf.h>
#include <infmp2.h>
#include <inftap.h>
#include <infpri.h>
#include <huckel.h>
#include <gnrinf.h>
#include <dftcom.h>
#include <dfterg.h>
C
      CHARACTER*6 RDMOLB(11)
      DATA RDMOLB/'CANORB','MOLORB','FORM18','FORM12','--5-- ',
     &            '--6-- ','OLDORB','H1DIAG','NEWORB','SIRIFC',
     &            'HUCKEL'/
      CHARACTER*33 TXTOCC(0:5)
      DATA TXTOCC/'NOTHING because of error(s).     ',
     &            'user input.                      ',
     &            'SIRIUS restart file (SIRIUS.RST).',
     &            'SIRIUS interface file (SIRIFC).  ',
     &            'diagonal of H1 matrix.           ',
     &            'Huckel guess.                    '/
C
      CHARACTER*70 STARS,RELEAS
      DATA STARS(1:36)  /'************************************'/
      DATA STARS(37:70) /'**********************************'/
      DATA RELEAS(1:36) /'*SIRIUS* a direct, restricted step, '/
      DATA RELEAS(37:70)/'second order MCSCF program       *'/
C
#include <orbtypdef.h>
C
C     ***** WRITE PROGRAM IDENTIFICATION *****
C
      WRITE(IWUNIT,200) STARS,RELEAS,STARS
 200  FORMAT(1X,A)
#ifdef VAR_AUTHORS
         WRITE(IWUNIT,400)
 400     FORMAT(/,T2,'Main authors of the SIRIUS module:'
     *    /T5,'Hans Joergen Aa. Jensen, University of Odense,  Denmark',
     *        ' (principal author)',
     *    /T5,'Hans Agren,              University of Uppsala, Sweden')
         CALL CIAUTH(IWUNIT)
         CALL TRAUTH(IWUNIT)
c        WRITE(IWUNIT,410)
c410     FORMAT(/,T2,'Questions and correspondence should be directed',
c    *        ' to H.J.Aa. Jensen',/,T5,'E-mail: hjj@dou.dk',
c    *        //,T2,'Standard reference:',
c    *        /,T5,'H.J.Aa. Jensen, H. Agren, and J. Olsen,',
c    *        ' Chapter 8 in MOTECC-90')
         IF (FLAG(16)) CALL SOLREF(IWUNIT)
#endif
C
C     stamp date, time and hostname to IWUNIT
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
      WRITE (IWUNIT,700) ((TITMOL(I,J),I=1,12), J=1,2)
  700 FORMAT(/' Title lines from integral program:'/T6,12A6/T6,12A6)
C
      WRITE(IWUNIT,1700) LUPRI,IPRI6,LUW4,IPRI4
 1700 FORMAT(/' Print level on unit LUPRI =',I4,' is',I4,
     *       /' Print level on unit LUW4  =',I4,' is',I4)
C
C Type of calculation:
C
 1050 FORMAT(/T6,A/)
 1051 FORMAT(T16,A)
      IF (IORTO .EQ. 1) THEN
         WRITE (IWUNIT,1050) 'Orthonormalization of molecular orbitals.'
         WRITE (IWUNIT,1100) (I,      I=1,NSYM)
         WRITE (IWUNIT,1110) ('  --', I=1,NSYM)
         WRITE (IWUNIT,1450) (NORB(I),I=1,NSYM)
         WRITE (IWUNIT,1500) (NBAS(I),I=1,NSYM)
         GO TO 8000
      ELSE IF ( DOMC .OR. DOCI ) THEN
         IF ( DOMC ) THEN
            WRITE (IWUNIT,1050) 'MC-SCF optimization.'
            IF ( DORSP ) WRITE (IWUNIT,1050)
     *         'Multi-configurational response calculation.'
         ELSE
            IF (DOCINO) THEN
               WRITE (IWUNIT,1050) 'CI-NO calculation.'
            ELSE
               WRITE (IWUNIT,1050)
     *            'Configuration interaction calculation.'
               IF ( DORSP ) WRITE (IWUNIT,1050)
     *            'CI response calculation.'
            END IF
         END IF
         IF (MCTYPE .EQ. 2) THEN
            IF (CORRLX) THEN
               WRITE(IWUNIT,1051)
     *         'Type: Relaxed single core hole calculation using RAS'
            ELSE IF (CORHOL) THEN
               WRITE(IWUNIT,1051)
     *         'Type: Frozen single core hole calculation using RAS'
            ELSE
               WRITE(IWUNIT,1051)
     *         'Type: restricted active space calculation (RAS)'
            END IF
         ELSE IF (FLAG(23)) THEN
            WRITE(IWUNIT,1051)
     *      'Type: reduced configuration space (not CAS).'
         ELSE
            WRITE(IWUNIT,1051)
     *      'Type: complete active space calculation (CAS).'
            IF (LUCME.GT.0) WRITE(LUCME,'(A/)')
     *      'Type: complete active space calculation (CAS).'
         END IF
      ELSE IF ( DOMP2 ) THEN
         WRITE (IWUNIT,1050) ' Moller-Plesset PT2 calculation.'
      ELSE IF ( DOCCSD) THEN
         WRITE (IWUNIT,1050) ' (Integral direct) CC calculation.'
      ELSE IF ( DORHF ) THEN
         IF (.NOT.DODFT) THEN
            IF (IOPRHF .EQ. 0) THEN
               WRITE(IWUNIT,1050)
     *            'Restricted, closed shell Hartree-Fock calculation.'
               WRITE(LUCME,'(A/)')
     *            'Restricted, closed shell Hartree-Fock calculation.'
            ELSE
               WRITE(IWUNIT,1050)
     *            'Restricted, one open shell Hartree-Fock calculation.'
            END IF
            IF ( DORSP ) WRITE (IWUNIT,1050)
     *         'Time-dependent Hartree-Fock calculation '//
     *         '(random phase approximation).'
         ELSE
            IF (IOPRHF .EQ. 0) THEN
               WRITE(IWUNIT,1050)
     *            'Restricted, closed shell Kohn-Sham calculation.'
               WRITE(LUCME,'(A/)')
     *            'Restricted, closed shell Kohn-Sham calculation.'
            ELSE
               WRITE(IWUNIT,1050)
     *            'Restricted, one open shell Kohn-Sham calculation.'
            END IF
            IF ( DORSP ) WRITE (IWUNIT,1050)
     *         'Time-dependent Kohn-Sham calculation '//
     *         '(random phase approximation).'
         END IF
      ELSE
         WRITE (IWUNIT,1050) 'No wave function optimization.'
         IF ( DORSP ) WRITE (IWUNIT,1050)
     *      'Response calculation.'
      END IF
      IF (NUMRUN .GE. 2) THEN
         WRITE(IWUNIT,1050)
     *      'This is a combination run starting with'
         IRUN = 1
         IF (DORHF) THEN
            IF (IOPRHF .EQ. 0) THEN
               WRITE(IWUNIT,1051) 'a RHF calculation'
            ELSE
               WRITE(IWUNIT,1051) 'an open shell RHF calculation'
            END IF
            IRUN = IRUN + 1
         END IF
         IF (DOMP2 .AND. IRUN .LT. NUMRUN) THEN
            WRITE(IWUNIT,1051) 'an MP2 calculation '
            IRUN = IRUN + 1
         END IF
         IF (FCVORB) THEN
         IF (NMVOT .EQ. 0) THEN
            WRITE(IWUNIT,1051)'diagonal h1 virtual orbitals'
         ELSE
            WRITE(IWUNIT,1051)'diagonal FC(MVO) virtual orbitals'
         END IF
         END IF
         IF (DOCI .AND. IRUN .LT. NUMRUN) THEN
            IF (DOCINO) THEN
               WRITE(IWUNIT,1051) 'a  CI-NO  calculation '
            ELSE
               WRITE(IWUNIT,1051) 'a  CI  calculation '
            END IF
            IRUN = IRUN + 1
         END IF
         WRITE (IWUNIT,'()')
      END IF
      IF ( DIRFCK ) THEN
      IF ( PARCAL ) THEN
         WRITE (IWUNIT,1051)
     *      'Fock matrices are calculated directly and in parallel',
     *      'without use of integrals on disk.'
      ELSE
         WRITE (IWUNIT,1051)
     *      'Fock matrices are calculated directly',
     *      'without use of integrals on disk.'
      END IF
      END IF
C
C ***** Hamiltonian specifications
C
      IF (NFIELD .GT. 0) THEN
         WRITE (IWUNIT,'(/A/A/A)')
     *      ' The molecule is placed in a static field.',
     *      ' Field strength (a.u.)          Field operator',
     *      ' ---------------------          --------------'
         DO 100 IFIELD = 1,NFIELD
            WRITE (IWUNIT,'(1P,G18.8,18X,A8)')
     *         EFIELD(IFIELD),LFIELD(IFIELD)
  100    CONTINUE
         WRITE (IWUNIT,'()')
      END IF
      IF (IRDMO .GT. 0 .AND. IRDMO .LE. 10) THEN
C
         JRDMO = IRDMO
         IF (ADDSTO) JRDMO = 11
C        ... if HUCKEL use JRDMO = 11 to get correct output
         WRITE (IWUNIT,'(/A/3A)')
     &   ' Initial molecular orbitals are obtained according to',
     &   ' ".MOSTART ',RDMOLB(JRDMO),'" input option.'
C
      IF (FLAG(16)) THEN
      IF (RSOL(1) .EQ. RSOL(2) .AND. RSOL(1) .EQ. RSOL(3)) THEN
         WRITE (IWUNIT,'(/5X,A,2(/16X,A,F12.6),/16X,A,I12)')
     *      'SOLVATION MODEL: molecule is in a spherical'//
     *        ' cavity in a dielectric medium,',
     *      'dielectric constant =',EPSOL,
     *      '      cavity radius =',RSOL(1),
     *      '        max l value =',LSOLMX
      ELSE
         WRITE (IWUNIT,'(/5X,A,2(/16X,A,3F12.6),/16X,A,I12)')
     *      'SOLVATION MODEL: molecule is in a ellipsoidal'//
     *        ' cavity in a dielectric medium,',
     *      'dielectric constant =',EPSOL,
     *      '  cavity dimensions =',RSOL(1),RSOL(2),RSOL(3),
     *      '        max l value =',LSOLMX
      END IF
         IF (INERSI) WRITE (IWUNIT,'(2(/5X,A))')
     *      'This is the initial state in a calculation with'//
     *      ' inertial polarization.',
     *      'Converged T(lm) values are saved on the Sirius'//
     *      ' interface file.'
         IF (INERSF) WRITE (IWUNIT,'(/5X,A/10X,A,F12.6)')
     *      'This is the final state in a calculation with'//
     *      ' inertial polarization.',
     *      'static dielectric constant =',EPSTAT
      END IF
      END IF
C
C ***** Wave function specifications
C
      NINEL = 2*NISHT
      WRITE(IWUNIT,1000) NINEL,NACTEL,NASHT,NORBT,ISPIN,NSYM,LSYM
      IF (LUCME.GT.0)
     *WRITE(LUCME,1000) NINEL,NACTEL,NASHT,NORBT,ISPIN,NSYM,LSYM
 1000 FORMAT(/5X,'Wave function specification'
     *       /5X,'============================',
     *       /5X,'Number of closed shell electrons',7X,I4,
     *       /5X,'Number of electrons in active shells',4X,I3,
     *       /5X,'Number of active orbitals',14X,I4,
     *       /5X,'Total number of orbitals',I19,
     *      //5X,'Spin multiplicity',22X,I4,
     *       /5X,'Total number of symmetries',I17,
     *       /5X,'Reference state symmetry',I19)
C
      IF (MCTYPE .EQ. 2) THEN
         WRITE (IWUNIT,1075) NAS1T, NAS2T, NAS3T
         IF (NAS1T .GT. 0) WRITE (IWUNIT,1071) NELMN1, NELMX1
         IF (NAS3T .GT. 0) WRITE (IWUNIT,1073) NELMN3, NELMX3
      END IF
 1071 FORMAT(/5X,'Minimum number of electrons in RAS1 :',I4,
     &       /5X,'Maximum number of electrons in RAS1 :',I4)
 1073 FORMAT(/5X,'Minimum number of electrons in RAS3 :',I4,
     &       /5X,'Maximum number of electrons in RAS3 :',I4)
 1075 FORMAT(/5X,'Number of active orbitals in RAS1, RAS2, and RAS3 :',
     &       3I4)
C
C     DFT run
C
      IF (DFTRUN) THEN
         WRITE (IWUNIT,'(/,4X,2A)')
     &        ' This is a DFT calculation of type: ',DFTTYP
         WRITE (IWUNIT,'(/,4X,A,3F12.6)')
     &        ' Dirac, Becke, and HF exchange weights:',
     &          WDFTX,WDFTB,HFXFAC
         WRITE (IWUNIT,'(4X,A,2F12.6)')
     &        ' WVN and LYP correlation weights:      ',
     &          WDFTC, WDFTL
         IF (DFTHRS) THEN
            WRITE (IWUNIT,'(4X,A,17X,3D12.2)')
     &        ' Redefined thresholds:', DFTHR0, DFTHRL, DFTHRI 
         ELSE
            WRITE (IWUNIT,'(4X,A,17X,3D12.2)')
     &        ' Default thresholds:  ', DFTHR0, DFTHRL, DFTHRI 
         END IF
         IF (DFTHRS) THEN
            WRITE (IWUNIT,'(4X,A,17X,3D12.2)')
     &        ' Redefined threshold for number of electrons: ', DFTELS 
         ELSE
            WRITE (IWUNIT,'(4X,A,17X,3D12.2)')
     &        ' Default threshold for number of electrons:   ', DFTELS 
         END IF
         WRITE (IWUNIT,'(4X,A,E10.2,I4)') ' DFT quadrature accuracy:',
     &        RADINT,ANGINT
         WRITE(IWUNIT,*) ANGINT
         IF (DFTPOT) WRITE (IWUNIT,'(/,A)')
     &        ' Kohn-Sham matrix calculated from potential.'
         IF (DFTASC) THEN
            WRITE (IWUNIT,'(/,4X,A)') ' Asymptotic correction applied'
            WRITE (IWUNIT,'(10X,A,F10.5)')
     &        ' Ionization potential used: ', DFTIPT      
            WRITE (IWUNIT,'(10X,A,2F10.5)')
     &        ' Radius factors:            ', DFTBR1, DFTBR2      
         END IF 
         IF (DFTHES) THEN
            WRITE (IWUNIT,'(/,4X,A)') 
     &        ' DFT electronic Hessian constructed.'
         END IF
      END IF
C
C ***** Orbital specifications *****
C
      WRITE(IWUNIT,1100) (I,      I=1,NSYM)
      WRITE(IWUNIT,1110) ('  --', I=1,NSYM)
      IF (MCTYPE .GT. 0) THEN
         WRITE(IWUNIT,1200) (NISH(I),I=1,NSYM)
         IF (NASHT .GT. 0) THEN
            IF (MCTYPE .LE. 1) THEN
               WRITE(IWUNIT,1300) (NASH(I),I=1,NSYM)
            ELSE
               IF (NAS1T .GT. 0) WRITE(IWUNIT,1320) (NAS1(I),I=1,NSYM)
               IF (NAS2T .GT. 0) WRITE(IWUNIT,1330) (NAS2(I),I=1,NSYM)
               IF (NAS3T .GT. 0) WRITE(IWUNIT,1340) (NAS3(I),I=1,NSYM)
            END IF
         END IF
         WRITE(IWUNIT,1400) (NSSH(I),I=1,NSYM)
      ELSE IF (.NOT.AUTOCC) THEN
         WRITE(IWUNIT,1550) (NRHF(I),I=1,NSYM)
         IF (IOPRHF .NE. 0) WRITE (IWUNIT,1555) (NASH(I),I=1,NSYM)
         WRITE(IWUNIT,1400) (NSSH(I),I=1,NSYM)
      END IF
      WRITE(IWUNIT,1450) (NORB(I),I=1,NSYM)
      WRITE(IWUNIT,1500) (NBAS(I),I=1,NSYM)
      IF (DORHF .AND. (AUTOCC .OR. MCTYPE .GT. 0)) THEN
         IF (AUTOCC) THEN
           IF (.NOT.DODFT) THEN
             WRITE(IWUNIT,'(/5X,A/5X,2A/5X,A)')
     &       ' ** Automatic occupation of RHF orbitals **',
     &       ' -- Initial occupation of symmetries is determined from ',
     &       TXTOCC(INIOCC),
     &       ' -- Initial occupation of symmetries is : --'
           ELSE
             WRITE(IWUNIT,'(/5X,A/5X,2A/5X,A)')
     &       ' ** Automatic occupation of RKS orbitals **',
     &       ' -- Initial occupation of symmetries is determined from ',
     &       TXTOCC(INIOCC),
     &       ' -- Initial occupation of symmetries is : --'
           END IF
         END IF
         WRITE(IWUNIT,'()')
         WRITE(IWUNIT,1550) (NRHF(I),I=1,NSYM)
         IF (IOPRHF .NE. 0) WRITE (IWUNIT,1560) IOPRHF
      END IF
      IF (SUPSYM) THEN
         CALL AVEPRT(IWUNIT)
      END IF
 1100 FORMAT(/5X,'Orbital specifications',
     *       /5X,'======================',
     *       /5X,'Abelian symmetry species',T38,8I4)
 1110 FORMAT(T38,8(A4))
 1200 FORMAT(5X,'Inactive orbitals',T38,8I4)
 1300 FORMAT(5X,'Active orbitals',T38,8I4)
 1320 FORMAT(5X,'Active orbitals in RAS1',T38,8I4)
 1330 FORMAT(5X,'Active orbitals in RAS2',T38,8I4)
 1340 FORMAT(5X,'Active orbitals in RAS3',T38,8I4)
 1400 FORMAT(5X,'Secondary orbitals',T38,8I4)
 1450 FORMAT(5X,'Total number of orbitals',T38,8I4)
 1500 FORMAT(5X,'Number of basis functions',T38,8I4)
 1550 FORMAT(5X,'Occupied HF orbitals',T38,8I4)
 1555 FORMAT(5X,'Open shell HF orbitals',T38,8I4)
 1560 FORMAT(/5X,'plus one open shell HF orbital of symmetry',I4)
C
C ***** Optimization information *****
C       FLAG(27) : use determinants
C
      IF (MCTYPE .GT. 0 .OR. (FLAG(21) .AND. .NOT. AUTOCC)) THEN
         IF (FLAG(27)) THEN
            WRITE (IWUNIT,5010) NCONF,NWOPT,NVAR
         ELSE
            WRITE (IWUNIT,5000) NCONF,NWOPT,NVAR
         END IF
      END IF
C*FM  change format nb of determinants and CSFs from I8 to I12
 5000 FORMAT (/,5X,'Optimization information',
     *        /,5X,'========================',
     *        /,5X,'Number of configurations   ',I12,
     *        /,5X,'Number of orbital rotations',I12,
     *        /,5X,'---------------------------------------',
     *        /,5X,'Total number of variables  ',I12)
 5010 FORMAT (/,5X,'Optimization information',
     *        /,5X,'========================',
     *        /,5X,'Number of determinants     ',I12,
     *        /,5X,'Number of orbital rotations',I12,
     *        /,5X,'---------------------------------------',
     *        /,5X,'Total number of variables  ',I12)
C
C
C
      IF (DOMC) THEN
         WRITE(IWUNIT,1600) MAXMAC,MAXMIC,THRMC,
     *                      LROOTS,MAXCIT,NROOTS
         IF (ISTATE .GT. 1) THEN
            WRITE(IWUNIT,'(/5X,A,I4/5X,A)')
     *         'This calculation converges to state no.',ISTATE,
     *         'for the specified symmetry and spin species.'
         ELSE
            WRITE(IWUNIT,'(/5X,A/5X,A)')
     *         'This calculation converges to the lowest state',
     *         'for the specified symmetry and spin species.'
         END IF
         WRITE(IWUNIT,2500) MAXJT
         IF (JOLSEN) WRITE(IWUNIT,2550)
         IF (MAXPHP.GT.1) THEN
            WRITE(IWUNIT,2560) MAXPHP
            IF (RESPHP) WRITE(IWUNIT,2570)
         END IF
      ELSE IF (DOCI) THEN
         WRITE(IWUNIT,1800) MAXMAC, THRCI, NROOCI
         IF (JOLSEN) WRITE(IWUNIT,2550)
         IF (MAXPHP.GT.1) THEN
            WRITE(IWUNIT,2560) MAXPHP
            IF (RESPHP) WRITE(IWUNIT,2570)
         END IF
         IF (DORSP ) THEN
            WRITE(IWUNIT,2600) ISTACI
         ELSE IF (ISTACI .GT. 0) THEN
            WRITE(IWUNIT,'(/5X,A,I4)')
     *      'The CI calculation converges only root no.',ISTACI
         END IF
      ELSE IF (DORHF) THEN
         WRITE(IWUNIT,1900) MAXFCK,MXDIIS,MXHFMA, THRRHF
      END IF
      IF (DOCCSD) CALL CC_PRTI
 1600 FORMAT(/5X,'Maximum number of macro iterations',I8,
     *       /5X,'Maximum number of micro iterations',I8,
     *       /5X,'Threshold for gradient',10X,1P,D10.2,
     *       /5X,'Number of initial trial vectors   ',I8,
     *       /5X,'Number of initial CI iterations   ',I8,
     *       /5X,'Number of simultaneous trial vectors',I6)
 1800 FORMAT(/5X,'Maximum number of CI iterations ',I10,
     *       /5X,'Threshold for CI convergence    ',1P,D10.2,
     *       /5X,'Number of CI roots to converge  ',I10)
 1900 FORMAT(/5X,'Maximum number of Fock  iterations',I8,
     &       /5X,'Maximum number of DIIS  iterations',I8,
     &       /5X,'Maximum number of QC-HF iterations',I8,
     *       /5X,'Threshold for HF convergence ',1P,D13.2)
 2500 FORMAT(/5X,'Maximum number of NEO/NR iterations',I4)
 2550 FORMAT(/5X,'CI trial vectors are generated with the Olsen',
     *       ' algorithm.')
 2560 FORMAT(/5X,'Dimension of CI PHP subspace' :,I6)
 2570 FORMAT(/5X,
     *   ' CI PHP subspace selected from size of residual elements.')
 2600 FORMAT(/5X,'Reference state in CI response calculation :',I8)
C
C ***
C
      NFROT = ISUM(NSYM,NFRO,1)
      NFRHT = ISUM(NSYM,NFRRHF,1)
      NFRMT = ISUM(NSYM,NFRMP2,1)
      IF (DORHF .AND. NFRHT .GT. 0) THEN
         IF (.NOT.DODFT) THEN
            WRITE (IWUNIT,6110) (NFRRHF(I),I=1,NSYM)
         ELSE
            WRITE (IWUNIT,6111) (NFRRHF(I),I=1,NSYM)
         END IF
      END IF
      IF (DOMP2 .AND. NFRMT .GT. 0) THEN
         WRITE (IWUNIT,6120) (NFRMP2(I),I=1,NSYM)
      END IF
      IF (DOMC .AND. NFROT.GT.0) THEN
         WRITE (IWUNIT,6100) (NFRO(I),I=1,NSYM)
      END IF
 6110 FORMAT(//,5X,'Number of frozen HF orbitals',T38,8I4)
 6111 FORMAT(//,5X,'Number of frozen KS orbitals',T38,8I4)
 6120 FORMAT(//,5X,'Number of frozen MP2 orbitals',T38,8I4)
 6100 FORMAT(//,5X,'Number of frozen MCSCF orbitals',T38,8I4)
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
 6300 FORMAT(/5X,'Number of orbitals frozen with ".FREEZE"',I3,
     &      //'     Orbital  Symmetry   No. in   Type',
     &       /'       no.              symmetry')
 6350 FORMAT(3I10,3X,A9)
C
 8000 CONTINUE
C
C *** End of subroutine PRTINP
C
      RETURN
      END
