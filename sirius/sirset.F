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
C File: sirius/sirset.F
C ===========================================================================
C  /* Deck sirini */
      SUBROUTINE SIRINI
C
C  2-Jan-1989 Hans Joergen Aa. Jensen
C
C  input/output unit numbers and buffer sizes
C  .. unit numbers are now defined with GPOPEN
C  initialize IROW
C
#include <implicit.h>
#include <iratdef.h>
C
#include <litinfo.h>
C ... for calculating LBINFO
C
C Used from common blocks:
C   INFIND : IROW(*)
C   INFDIM : LBUFMA
C   INFTAP : LB*
C
#include <maxash.h>
#include <maxorb.h>
#include <infind.h>
#include <infdim.h>
#include <inftap.h>
C#include <infpri.h>
C
      CALL QENTER('SIRINI')
C
C PART 1 : ******* about some I/O units used in SIRIUS
C
C LUCMD    input unit for program control
C LUPRI    output unit for print output
C LUW4     do.
C LUPNCH   Formatted output of final MCSCF or natural/canonical orbitals
C LUERR    output unit for error and warning messages
C LUSTAT   output unit for statistics
C LUINDX   orbital operators of symmetry 1 to NSYM.
C LUSYMB   symbolic matrix elements
C LUINTA   AO integral file
C LUINTM   MO integral file, Mulliken order (1 1 / 2 2)
C LUINTD   MO integral file, Dirac    order <1 2 / 1 2>
C LUONEL   AO one electron integrals
C LUSUPM   AO super matrix file
C LUSOL    Dielectric solvation integrals
C LUPROP   property integral file
C LUIT*    MC iteration files
C          LUIT1 is added for saving old orbitals and CI-vectors
C                for 1) start and 2) back step.
C          LUIT2 is for L-matrix diagonals
C          LUIT3 and LUIT5 are for trial and sigma vectors, resp.
C          LUIT7 is for update information in SIRUPD
C LUINF    file for writing information on convergence etc.
C
C *** Define buffer sizes (LBXXXX is buffer size for file LUXXXX)
C     LBINTM defined in SIR_INTOPEN (891107-hjaaj)
C     LBONEL defined in RDONEL, if LBONEL .lt. 0 (900221)
C     MAERKE 900108-hjaaj: LBXXXX should always be defined in routines
C
      LBINFO=LDINFO+LIINFO/IRAT
      LBONEL=600
      LBUFMA=LBONEL
C
C     loop to set up triangular indexing array ...
C
      IF (MXCORB .GE. LIROW) THEN
         CALL QUIT('SIRINI error: MXCORB exceeds LIROW')
      END IF
      IROW(1) = 0
      DO 400 I = 1,(LIROW-1)
         IROW(I+1) = IROW(I) + I
  400 CONTINUE
C
      CALL QEXIT('SIRINI')
      RETURN
      END
C  /* Deck sirset */
      SUBROUTINE SIRSET(WORK,LFREE,OLDWOP)
C
C
C  2-Jan-1989 Hans Joergen Aa. Jensen
C  890712-hjaaj: OLDWOP parameter
C  900803-hjaaj: JWOPSY must be defined before call of SIRSET.
C  910708-hjaaj: corrected some combination run problems
C
C SIRIUS SETUP ROUTINE FOR COMMON-BLOCKS AND INDEXING ARRAYS
C
C CONVENTIONS for variable names, particularly in INFORB:
C
C  PREFIXES N, NN, N2  FOR SIZES OF STRAIGHT, LOWER TRIANGULAR AND
C                      SQUARED BLOCK ARRAYS WITHIN A SYMMETRY.
C  PREFIXES I, II, I2  FOR SIZES OF CUMULATIVE VALUES OF LOWER
C                      SYMMETRIES OF STRAIGHT, LOWER TRIANGULAR AND
C                      SQUARED BLOCK ARRAYS.
C  SUFFIX  MA          FOR MAXIMUM VALUE OVER SYMMETRIES
C  SUFFIX  T           FOR SUM OF TRIANGLES OR SQUARES OVER SYMMETRIES
C  SUFFIX  X           WHEN SUM OF ARRAY IS TAKEN BEFORE BEING
C                      TRIANGULARIZED OR SQUARED.
C  SUFFIX  Y           EXAMPLE: NASHY = NASHX*(NASHX + 1)/2
C  SUFFIX  (K)         SYMMETRY NUMBER
C
#include <implicit.h>
      DIMENSION WORK(LFREE)
      LOGICAL   OLDWOP
#include <iratdef.h>
#include <litinfo.h>
C -- local constants
      PARAMETER (D2=2.0D0)
#include <dummy.h>
C
C Used from common blocks:
C   INFINP : LSYM,ITRLVL,?
C   INFVAR : NCONF,NWOPT,NVAR,JWOPSY
C   INFORB : MULD2H(8,8),NISH(8),...
C   INFIND : ?
C   INFDIM : NWOPMA,MAXPHP,LPHPMX,...
C   INFOPT : ?
C   INFLIN : most of it (defined here)
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infvar.h>
#include <inforb.h>
#include <infind.h>
#include <infdim.h>
#include <infopt.h>
#include <inflin.h>
#include <inftap.h>
#include <infpri.h>
C
      LOGICAL FIRST, GETWOP
      CHARACTER*8 OPLBL(8)
      SAVE    FIRST
      DATA    FIRST /.TRUE./
      DATA OPLBL /'EXOPSYM1','EXOPSYM2','EXOPSYM3','EXOPSYM4',
     *            'EXOPSYM5','EXOPSYM6','EXOPSYM7','EXOPSYM8'/
C
C
      CALL QENTER('SIRSET')
C
C Initialize
C
      IF (FIRST) THEN
         NWOPMA = 0
      END IF
C
C PART 1 : ***** SET UP ORBITAL DATA ****************************
C
      CALL SETORB
C
C PART 2 : **** CI PARAMETERS ***********************************
C
      SPIN = (ISPIN-1)/D2
      CALL SETCI(NCONF,NCDETS,LSYM,WORK,LFREE,0)
      NCONDI = MAX(1,NCONF)
C
C PART 3 : **** define /INFLIN/ *********************************
C          (NWOPPT and NVARPT are defined in GETWOP)
C
      IPRLIN = IPRSIR
      LSYMRF = LSYM
      LSYMPT = JWOPSY
      LSYMST = MULD2H(LSYMRF,LSYMPT)
      NCONRF = NCONF
      NCONST = NCONF
C
C PART 4 : **** ORBITAL ROTATION PARAMETERS *********************
C
C     Set orbital rotation index vectors JWOP, KLWOP, ? :
C
      IF ( (DOMP2 .OR. DOCI .OR. DOCINO) .AND.
     &     .NOT.(DORHF .OR. DOMC)) THEN
         NWOPT  = 0
         NWOPH  = 0
         NWOPDI = 1
         NWOPPT = 0
         NVARPT = NCONST + NWOPPT
      ELSE
         IF (OLDWOP) THEN
            IF ( .NOT. GETWOP(OPLBL(JWOPSY)) ) THEN
               WRITE(LUERR,5000) OPLBL(JWOPSY)
               WRITE(LUW4,5000) OPLBL(JWOPSY)
               IF (LUPRI.NE.LUW4) WRITE(LUPRI,5000) OPLBL(JWOPSY)
 5000          FORMAT(/' SIRSET FATAL ERROR, OLDWOP true but "',
     *                 A,'" not found on LUINDX.')
               CALL QTRACE(LUERR)
               CALL QUIT('ERROR in SIRSET call with OLDWOP = .TRUE.')
            END IF
         ELSE
            MWOPMA = NWOPMA
            CALL SETWOP(WORK,LFREE)
            NWOPMA = MAX(NWOPMA,MWOPMA)
            NWOPPT = NWOPT
            NVARPT = NCONST + NWOPPT
         END IF
      END IF
C
      NVAR   = NCONF + NWOPT
      NVARH  = NCONF + NWOPH
      NVARMA = NCONMA+ NWOPMA
C
C
C PART 5 : **** MISCELLANEOUS ******************************************
C
C     Define PHP parameters for diagonals
C
      CALL PHPINI(LPHPMX,NCONF,NWOPT,MAXPHP)
C
C     Test for inconsistencies:
C
      NUMERR = 0
      IF (DOMC .AND. ITRLVL.LT.2) THEN
         IF (NWOPT.GT.0 .AND. NASHT.GT.1 .AND. .NOT.FLAG(20)) THEN
C                  ^- not CI        ^- not RHF      ^- not only grad
            WRITE (LUW4,1020) ITRLVL
            IF (LUPRI.NE.LUW4) WRITE (LUPRI,1020) ITRLVL
            WRITE (LUERR,1020) ITRLVL
            NUMERR = NUMERR + 1
         END IF
      END IF
      IF (NCONF.LE.0) THEN
         WRITE (LUW4,1610) NCONF
         IF (LUPRI.NE.LUW4) WRITE (LUPRI,1610) NCONF
         WRITE (LUERR,1610) NCONF
         NUMERR = NUMERR + 1
      END IF
      IF (NASHT.GT.1 .AND. FLAG(21) .AND.
     *    .NOT.(DOCI .OR. DOCINO .OR. DOMC .OR. HSROHF) )THEN
         WRITE (LUW4,1730)
         IF (LUPRI.NE.LUW4) WRITE (LUPRI,1730)
         WRITE (LUERR,1730)
         NUMERR = NUMERR + 1
      END IF
      IF (NASHT.EQ.1 .AND. NACTEL.NE.1) THEN
         WRITE (LUW4,1731)
         IF (LUPRI.NE.LUW4) WRITE (LUPRI,1731)
         WRITE (LUERR,1731)
         NUMERR = NUMERR + 1
      END IF
 1020 FORMAT(/' SIRSET-ERROR, Transformation level must be 2, 3, or 4 ',
     *       /' for MCSCF optimization; specified level is:',I3)
 1610 FORMAT(/,' SIRSET-ERROR, NCONF .le. 0 : NCONF =',I10)
 1730 FORMAT(/,' SIRSET-ERROR, NASHT .gt. 1 inconsistent with ',
     *   'RHF calculation')
 1731 FORMAT(/,' SIRSET-ERROR, NASHT .eq. 1 is only implemented',
     *       /,'               for one active electron (doublet)')
C
C ***
C
      IF (NUMERR .GT. 0) THEN
         IF (FIRST) CALL SIR_PRTINP(LUW4,-1,-1)
         WRITE (LUW4,'(///A/)')
     *      ' SIRSET - FATAL ERRORS DETECTED, SEE ABOVE'
         IF (LUPRI .NE. LUW4) THEN
            IF (FIRST) CALL SIR_PRTINP(LUPRI,-1,-1)
            WRITE (LUPRI,'(///A/)')
     *         ' SIRSET - FATAL ERRORS DETECTED, SEE ABOVE'
         END IF
         WRITE (LUERR,'(//A)')
     *       ' *** ERROR *** input is inconsistent (SIRSET)'
         CALL QTRACE(LUERR)
         CALL QUIT('*** ERROR *** input inconsistent (SIRSET)')
      END IF
C
      FIRST = .FALSE.
C
      CALL QEXIT('SIRSET')
      RETURN
      END
C  /* Deck sir_open */
      SUBROUTINE SIR_OPEN
C
C  5-May-1989 Hans Joergen Aa. Jensen
C
C  Open input/output units required by this Sirius calculation,
C  as specified in command input.
C  Mo integral files are opened in SIR_INTOPEN
C
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
C
C Used from common blocks:
C   INFINP : FLAG(*)
C   INFTAP : LU*
C   INFPRI : NINFO
C
#include <maxorb.h>
#include <infinp.h>
#include <inftap.h>
#include <infpri.h>
C
      LOGICAL FEXIST
C
C     ***** CALL SIR_INTOPEN to open MO integral files
C
      CALL SIR_INTOPEN
C
C     ***** Open other files
C
      IF (FLAG(2) .OR. FLAG(21)) CALL GPOPEN(LUINF,'SIRIUS.ITINFO',
     &     'UNKNOWN',' ','UNFORMATTED',IDUMMY,.FALSE.)
      IF (FLAG(16))
     *   CALL GPOPEN(LUSOL,FNSOL,'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
C
      CALL GPOPEN(LUIT1,'SIRIUS.RST','UNKNOWN',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      CALL GPOPEN(LUIT2,' ','UNKNOWN',' ','UNFORMATTED',IDUMMY,.FALSE.)
      CALL GPOPEN(LUIT3,' ','UNKNOWN',' ','UNFORMATTED',IDUMMY,.FALSE.)
      CALL GPOPEN(LUIT5,' ','UNKNOWN',' ','UNFORMATTED',IDUMMY,.FALSE.)
      RETURN
      END
C  /* Deck setorb */
      SUBROUTINE SETORB
C
C  21-Sep-1988 Hans Joergen Aa. Jensen
C
C  Purpose: set common block /INFORB/ based on
C           NSYM, MULD2H, NISH(*), NASH(*), NORB(*), NBAS(*)
C           and           NRHF(*), NFRO(*)
C
C           and ISW ISX ISMO ISAO  NSM ICH  IOBTYP IACTYP
C           in common block /INFIND/
C
C           and max dimensions in /INFDIM/
C
#include <implicit.h>
#include <priunit.h>
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
#include <infdim.h>
#include <r12int.h>
C
C
C     (1) common /INFORB/
C
      NRHFT  = 0
      NVIRT  = 0
      NISHT  = 0
      NASHT  = 0
      NAS1T  = 0
      NAS2T  = 0
      NAS3T  = 0
      NORBT  = 0
      NBAST  = 0
      NNRHFT = 0
      NNVIRT = 0
      NNOCCT = 0
      NNORBT = 0
      NNBAST = 0
      N2RHFT = 0
      N2VIRT = 0
      N2OCCT = 0
      N2ORBT = 0
      N2BAST = 0
      NCMOT  = 0

      N_ERRORS = 0

      DO 100 I = 1,NSYM
         NOCC(I)  = NISH(I) + NASH(I)
         NSSH(I)  = NORB(I) - NOCC(I)
         IF (NSSH(I) .LT. 0) N_ERRORS = N_ERRORS + 1
         NNORB(I) = (NORB(I)*(NORB(I)+1))/2
         NVIR(I)  = NORB(I) - NRHF(I)
         NNBAS(I) = (NBAS(I)*(NBAS(I)+1))/2
         N2ORB(I) = NORB(I)*NORB(I)
         N2BAS(I) = NBAS(I)*NBAS(I)
         NRHFT    = NRHFT   + NRHF(I)
         NVIRT    = NVIRT   + NVIR(I)
         NISHT    = NISHT   + NISH(I)
         NASHT    = NASHT   + NASH(I)
         NAS1T    = NAS1T   + NAS1(I)
         NAS2T    = NAS2T   + NAS2(I)
         NAS3T    = NAS3T   + NAS3(I)
         NORBT    = NORBT   + NORB(I)
         NBAST    = NBAST   + NBAS(I)
         NNRHFT   = NNRHFT  + (NRHF(I)*(NRHF(I)+1))/2
         NNVIRT   = NNVIRT  + (NVIR(I)*(NVIR(I)+1))/2
         NNOCCT   = NNOCCT  + (NOCC(I)*(NOCC(I)+1))/2
         NNORBT   = NNORBT  + (NORB(I)*(NORB(I)+1))/2
         NNBAST   = NNBAST  + (NBAS(I)*(NBAS(I)+1))/2
         N2RHFT   = N2RHFT  + NRHF(I)*NRHF(I)
         N2VIRT   = N2VIRT  + NVIR(I)*NVIR(I)
         N2OCCT   = N2OCCT  + NOCC(I)*NOCC(I)
         N2ORBT   = N2ORBT  + NORB(I)*NORB(I)
         N2BAST   = N2BAST  + NBAS(I)*NBAS(I)
         NCMOT    = NCMOT   + NORB(I)*NBAS(I)
C        For multiple basis sets (WK/UniKA/04-11-2002).
         NORB1(I) = NORB(I)       
  100 CONTINUE
C
      IF (N_ERRORS .GT. 0) THEN
         WRITE(lupri,'(//A)') 'SETORB FATAL ERROR'//
     &      ' Number of occupied orbitals bigger than total number'//
     &      ' of orbitals !'
         WRITE(lupri,*) 'NISH ',NISH
         WRITE(lupri,*) 'NASH ',NASH
         WRITE(lupri,*) 'NOCC ',NOCC
         WRITE(lupri,*) 'NSSH ',NSSH
         WRITE(lupri,*) 'NORB ',NORB
         WRITE(lupri,*) 'NBAS ',NBAS
         CALL QENTER('SETORB')
         CALL QUIT('No. of occ. orb. > total no. of orb.')
      END IF
C
      IISH(1)  = 0
      IASH(1)  = 0
      ISSH(1)  = 0
      IOCC(1)  = 0
      IORB(1)  = 0
      IBAS(1)  = 0
      IIISH(1) = 0
      IIASH(1) = 0
      IIORB(1) = 0
      IIBAS(1) = 0
      I2ORB(1) = 0
      I2BAS(1) = 0
      ICMO(1)  = 0
      DO 200 I = 2,NSYM
         IISH(I)  = IISH(I-1)  + NISH(I-1)
         IASH(I)  = IASH(I-1)  + NASH(I-1)
         ISSH(I)  = ISSH(I-1)  + NSSH(I-1)
         IOCC(I)  = IOCC(I-1)  + NOCC(I-1)
         IORB(I)  = IORB(I-1)  + NORB(I-1)
         IBAS(I)  = IBAS(I-1)  + NBAS(I-1)
         IIISH(I) = IIISH(I-1) + (NISH(I-1)*(NISH(I-1)+1))/2
         IIASH(I) = IIASH(I-1) + (NASH(I-1)*(NASH(I-1)+1))/2
         IIORB(I) = IIORB(I-1) + NNORB(I-1)
         IIBAS(I) = IIBAS(I-1) + NNBAS(I-1)
         I2ORB(I) = I2ORB(I-1) + N2ORB(I-1)
         I2BAS(I) = I2BAS(I-1) + N2BAS(I-1)
         ICMO(I)  = ICMO(I-1)  + NORB(I-1)*NBAS(I-1)
  200 CONTINUE
C
      NOCCT  = NISHT + NASHT
      NSSHT  = NORBT - NOCCT
C
      NNRHFX = (NRHFT*NRHFT+NRHFT)/2
      N2RHFX = NRHFT*NRHFT
      NNVIRX = (NVIRT*NVIRT+NVIRT)/2
      N2VIRX = NVIRT*NVIRT
      N2ISHX = NISHT*NISHT
      NNASHX = (NASHT*(NASHT+1))/2
      N2ASHX = NASHT*NASHT
      NNASHY = (NNASHX*(NNASHX+1))/2
      NNOCCX = (NOCCT*(NOCCT+1))/2
      N2OCCX = NOCCT*NOCCT
      NNORBX = (NORBT*(NORBT+1))/2
      N2ORBX = NORBT*NORBT
      NNBASX = (NBAST*(NBAST+1))/2
      N2BASX = NBAST*NBAST
C
C     (2) common /INFIND/
C
C     Orbital index arrays:
C     ISW ISX ISMO ISAO JWOP NSM ICH  IOBTYP IACTYP
C
C     setup of IOBTYP and IACTYP
C     (JTFRO=frozen, JTINAC=inactive, JTACT=active, JTSEC=secondary)
C
      INDX = 0
      INDT = 0
      DO 500 ISYM = 1,NSYM
         DO 600 NF = 1,NFRO(ISYM)
            INDX = INDX+1
  600       IOBTYP(INDX) = JTFRO
         DO 700 NI = NFRO(ISYM)+1,NISH(ISYM)
            INDX = INDX+1
  700       IOBTYP(INDX) = JTINAC
         DO 800 NT = NISH(ISYM)+1,NOCC(ISYM)
            INDX = INDX+1
  800       IOBTYP(INDX) = JTACT
         DO 900 NA = NOCC(ISYM)+1,NORB(ISYM)
            INDX = INDX+1
  900       IOBTYP(INDX) = JTSEC
C
         DO 910 NT = 1,NAS1(ISYM)
  910       IACTYP(INDT+NT) = 1
         INDT = INDT + NAS1(ISYM)
         DO 920 NT = 1,NAS2(ISYM)
  920       IACTYP(INDT+NT) = 2
         INDT = INDT + NAS2(ISYM)
         DO 930 NT = 1,NAS3(ISYM)
  930       IACTYP(INDT+NT) = 3
         INDT = INDT + NAS3(ISYM)
  500 CONTINUE
C     setup of ISAO and ISMO, symmetry index of basis fu. and orbitals
      II = 0
      JJ = 0
      DO 1000 ISYM = 1,NSYM
         NBASI = NBAS(ISYM)
         NORBI = NORB(ISYM)
         DO 1010 I = 1,NBASI
            II = II+1
            ISAO(II) = ISYM
            IF (I.LE.NORBI) THEN
               JJ = JJ+1
               ISMO(JJ) = ISYM
            END IF
 1010    CONTINUE
 1000 CONTINUE
C
C     Compute reordering indices for molecular orbitals.
C      ISW reorders to inactive-active-secondary ordering
C      ISX reorders back
C
C     ICH and NSM from old CASSCF CICTL
C
C
      IIN  = 0
      IAC  = NISHT
      ISEC = NOCCT
      IU   = 0
      IR   = 0
      I    = 0
      DO 1030 ISYM = 1,NSYM
         NISHI = NISH(ISYM)
         DO 1040 KISH = 1,NISHI
            I   = I   + 1
            IIN = IIN + 1
            IR  = IR  - 1
            ISW(I)   = IIN
            ISX(IIN) = I
            ICH(I)   = IR
 1040    CONTINUE
         NASHI = NASH(ISYM)
         DO 1050 KASH = 1,NASHI
            I   = I   + 1
            IAC = IAC + 1
            IU  = IU  + 1
            ISW(I)   = IAC
            ISX(IAC) = I
            ICH(I)   = IU
            NSM(IU)  = ISYM
 1050    CONTINUE
         NSSHI = NSSH(ISYM)
         DO 1060 KSEC = 1,NSSHI
            ISEC = ISEC + 1
            I    = I    + 1
            ISW(I)    = ISEC
            ISX(ISEC) = I
            ICH(I)    = 0
 1060    CONTINUE
 1030 CONTINUE
C
C     (3) common /INFDIM/
C
      NISHMA = 0
      NASHMA = 0
      NOCCMA = 0
      NSSHMA = 0
      NORBMA = 0
      NBASMA = 0
      NNBASM = 0
      N2BASM = 0
      DO 2000 ISYM = 1,NSYM
         NISHMA = MAX(NISHMA,NISH(ISYM))
         NASHMA = MAX(NASHMA,NASH(ISYM))
         NSSHMA = MAX(NSSHMA,NSSH(ISYM))
         NOCCMA = MAX(NOCCMA,NOCC(ISYM))
         NORBMA = MAX(NORBMA,NORB(ISYM))
         NBASMA = MAX(NBASMA,NBAS(ISYM))
         NNBASM = MAX(NNBASM,NNBAS(ISYM))
         N2BASM = MAX(N2BASM,N2BAS(ISYM))
 2000 CONTINUE
C
C
      NISHDI = MAX(1,NISHT)
      NASHDI = MAX(1,NASHT)
C     end of SETORB
      RETURN
      END
C  /* Deck setwop */
      SUBROUTINE SETWOP(KLWOP,LKLWOP)
C
C   2-Jan-1989 Hans Joergen Aa. Jensen
C
C  Purpose: set orbital rotation index vectors
C
#include <implicit.h>
#include <iratdef.h>
#include <dummy.h>
C
      DIMENSION KLWOP(NOCCT,NORBT)
C
C Used from common blocks:
C   INFINP : ISTATE,MCTYPE,FLAG(*),?
C   INFVAR : MAXOCC,NWOPT,JWOPSY,JWOP(2,*)
C   INFORB : NISH(8),...
C   INFIND : ICH(*),?
C   INFDIM : NWOPMA,?
C
#include <priunit.h>
#include <maxash.h>
#include <maxorb.h>
#include <infinp.h>
#include <infvar.h>
#include <inforb.h>
#include <infind.h>
#include <infdim.h>
#include <infpri.h>
#include <gnrinf.h>
C -- local variables
C
      LOGICAL ACTALL, MAKABS, FIRST
C
      CHARACTER*8 OPLBL(8), ABSLBL(3), LAB123(3)
      DATA OPLBL /'EXOPSYM1','EXOPSYM2','EXOPSYM3','EXOPSYM4',
     *            'EXOPSYM5','EXOPSYM6','EXOPSYM7','EXOPSYM8'/
      DATA ABSLBL/'ABS1SYM1','ABS2SYM1','ABS3SYM1'/
      DATA LAB123/3*'********'/
      DATA FIRST /.TRUE./
      SAVE FIRST
C
C
      CALL QENTER('SETWOP')
C
C
C   ********Pointer JWOP for orbital operators of symmetry IOPSY
C           (Brillouin-matrix elements ).
C           Arrangement is symmetry first, then within
C           each symmetry order is IT,IA,TA.
C
      IF (WRINDX) NWOPMA = 0
C
C If active-active rotations specified with ACTROT in SIRINP
C check if the specified rotations are OK
C
      IF (FIRST) THEN
         NWOPA = NWOPT
      ELSE
         NWOPA = 0
      END IF
      IACT  = 0
      IF (FLAG(23) .AND. NWOPA.GT.0) THEN
         INPERR = 0
         DO 3210 I = 1,NWOPA
            K = JWOP(1,I)
            L = JWOP(2,I)
            IF (K.GE.L) INPERR = INPERR + 1
            K = ICH(K)
            L = ICH(L)
            IF (K.LE.0 .OR. L.LE.0) THEN
               INPERR = INPERR + 100
            ELSE IF (NSM(K).NE.NSM(L)) THEN
               INPERR = INPERR + 10000
            END IF
 3210    CONTINUE
         IF (INPERR.GT.0) THEN
            WRITE (LUERR,3212) INPERR
            WRITE (LUPRI,3212) INPERR
            WRITE (LUPRI,3214) (I,JWOP(1,I),JWOP(2,I),I=1,NWOPA)
            CALL QTRACE(LUPRI)
            CALL QUIT('*** ERROR-SETWOP *** errors in ACTROT input')
         END IF
         ACTALL = .FALSE.
      ELSE IF (FLAG(23)) THEN
         ACTALL = .TRUE.
      ELSE
         ACTALL = .FALSE.
      END IF
 3212 FORMAT(/'  SETWOP-ERROR, error in specified active rotations',
     *  ' code =',I6)
 3214 FORMAT(/' -SETWOP-INFO,  list of specified rotations:',/,
     *  (I10,':',2I4))
C
      IF (FLAG(11) .AND. NWOPA.GT.0) THEN
         WRITE (LUPRI,'(/A,A)') ' *** ERROR ***',
     *      ' (SETWOP) ACTROT option not implemented for FLAG(11)'
         CALL QTRACE(LUPRI)
         CALL QUIT('*ERROR-SETWOP**ACTROT not implemented for FLAG(11)')
      END IF
C
C     If .NOT.FLAG(11) we only want the totally symmetric operator
C     .. maybe, but now generalized so it should work ok also for
C        jwopsy .ne. 1 /hjaaj aug 2004
C
      IF (FLAG(11) .AND. WRINDX) THEN
         IOPSY1 = 1
         NOPSY = NSYM
      ELSE
         IOPSY1 = JWOPSY
         NOPSY = JWOPSY
      END IF
      MAKABS = (FLAG(51) .OR. FLAG(52) .OR. FLAG(53))
Cold  WRINDX = (FLAG(11) .OR. MAKABS)
C890712-hjaaj: now always wrindx in first call,
C              in order to implement combination runs.
C     used to be: WRINDX = FIRST, now for DALTON:
C     WRINDX moved to common block GNRINF in order to have external control
C
      LUINDX = -1
      IF (WRINDX) THEN
         CALL GETDAT(LAB123(2),LAB123(3))
         CALL GPOPEN(LUINDX,'LUINDF','UNKNOWN',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         REWIND LUINDX
      END IF
C
      DO 1295 IOPSY = IOPSY1,NOPSY
         LVLABS = 0
 1292    CONTINUE
         NWOPX  = NWOPA
         NWOPX1 = 0
         DO 1290 ISYM = 1,NSYM
C
            ISTI  = IORB(ISYM) + NFRO(ISYM) + 1
            IENDI = IORB(ISYM) + NISH(ISYM)
            ISTT  = IENDI + 1
            IENDT = IENDI + NASH(ISYM)
C
            JSYM  = MULD2H(ISYM,IOPSY)
            JSTT  = IORB(JSYM) + NISH(JSYM) + 1
            JENDT = IORB(JSYM) + NOCC(JSYM)
            JSTA  = JENDT + 1
            JENDA = JENDT + NSSH(JSYM)
C
C           *** INACTIVE-ACTIVE ELEMENTS
C
            IF (NISH(ISYM).EQ.0.OR.NASH(JSYM).EQ.0) GO TO 1200
            DO 1210 IN = ISTI,IENDI
               IF (NOROT(IN).EQ.0) THEN
                  DO 1220 JT = JSTT,JENDT
                     IF (NOROT(JT).EQ.0) THEN
                     IF (IOPSY.NE.1 .OR. ISSMO(IN) .EQ. ISSMO(JT)) THEN
                        NWOPX = NWOPX + 1
                        JWOP(1,NWOPX) = IN
                        JWOP(2,NWOPX) = JT
                     END IF
                     END IF
 1220             CONTINUE
               END IF
 1210       CONTINUE
            IF (NWOPX .GT. MAXWOP) GOTO 2500
C
C           *** INACTIVE-SECONDARY ELEMENTS
C
 1200       IF (NISH(ISYM).EQ.0.OR.NSSH(JSYM).EQ.0) GO TO 1230
            IF (LVLABS .EQ. 1) GO TO 1230
            DO 1240 IN = ISTI,IENDI
               IF (NOROT(IN).EQ.0) THEN
                  DO 1250 JA = JSTA,JENDA
                     IF (NOROT(JA).EQ.0) THEN
                     IF (IOPSY.NE.1 .OR. ISSMO(IN) .EQ. ISSMO(JA)) THEN
                        NWOPX = NWOPX + 1
                        JWOP(1,NWOPX) = IN
                        JWOP(2,NWOPX) = JA
                     END IF
                     END IF
 1250             CONTINUE
               END IF
 1240       CONTINUE
            IF (NWOPX .GT. MAXWOP) GOTO 2500
C
C           *** active-active elements, if RAS or
C           if FLAG(23) (i.e. active-active rotations) or "orbital
C           absorption" on ground state (we assume CI frozen in that
C           case, i.e. active-active rotations become non-redundant).
C           For excited states we may loose control by allowing
C           active-active rotations.
C
 1230       IF (NASH(ISYM).EQ.0.OR.NASH(JSYM).EQ.0) GO TO 1165
            IF (ACTALL .OR. (LVLABS .GT. 0 .AND. ISTATE .EQ. 1) ) THEN
               IACT  = 1
               JSTTX = JSTT
               DO 1175 IT = ISTT,IENDT
                  IF (NOROT(IT).NE.0) GO TO 1175
                  IF (ISYM.EQ.JSYM) JSTTX = IT + 1
                  DO 1185 JT = JSTTX,JENDT
                     IF (NOROT(JT).NE.0) GO TO 1185
                     IF (IOPSY.NE.1 .OR. ISSMO(IT) .EQ. ISSMO(JT)) THEN
                        NWOPX = NWOPX + 1
                        JWOP(1,NWOPX) = IT
                        JWOP(2,NWOPX) = JT
                     END IF
 1185             CONTINUE
 1175          CONTINUE
            ELSE IF (MCTYPE .EQ. 2) THEN
C
C              RAS active-active rotations
C
               IF (.NOT.FLAG(24)) CALL RASWOP(ISYM,JSYM,NWOPX,JWOP)
               IACT = 0
            ELSE
               IACT = 0
            END IF
            IF (NWOPX .GT. MAXWOP) GOTO 2500
C           *** ACTIVE-SECONDARY ELEMENTS
 1165       IF (NASH(ISYM).EQ.0.OR.NSSH(JSYM).EQ.0) GO TO 1260
            IF (LVLABS .EQ.1 .OR. LVLABS .EQ. 2) GO TO 1260
            DO 1270 IT = ISTT,IENDT
               IF (NOROT(IT).NE.0) GO TO 1270
               DO 1280 JA = JSTA,JENDA
                  IF (NOROT(JA).NE.0) GO TO 1280
                  IF (IOPSY.NE.1 .OR. ISSMO(IT) .EQ. ISSMO(JA)) THEN
                     NWOPX = NWOPX+1
                     JWOP(1,NWOPX) = IT
                     JWOP(2,NWOPX) = JA
                  END IF
 1280          CONTINUE
 1270       CONTINUE
 1260       CONTINUE
            NWOPX1 = NWOPX
            IF (NWOPX .GT. MAXWOP) GOTO 2500
 1290    CONTINUE
C
C        Construct KLWOP and remove any duplicate active-active
C        rotations
C
         NKLWOP = NOCCT*NORBT
         IF (NKLWOP/IRAT + 1 .GT. LKLWOP)
     &      CALL ERRWRK('SETWOP',(NKLWOP/IRAT),LKLWOP)
         CALL IZERO(KLWOP,NKLWOP)
         NWOPX = 0
         DO 140 IG = 1,NWOPX1
            KW = ISW(JWOP(1,IG))
            LW = ISW(JWOP(2,IG))-NISHT
            IF (KLWOP(KW,LW) .EQ. 0) THEN
               NWOPX = NWOPX + 1
               KLWOP(KW,LW)  = NWOPX
               JWOP(1,NWOPX) = JWOP(1,IG)
               JWOP(2,NWOPX) = JWOP(2,IG)
            END IF
  140    CONTINUE
         IF (NWOPX .GT. MAXWOP) GOTO 2500
         NWOPT = NWOPX
C
C        Fill up with redundant active-active orbital rotations,
C        they are needed for LINTRN.
C
         IF (IOPSY.EQ.1 .AND. IACT .EQ. 0) THEN
            DO 2400 KW = NISHT+1, NOCCT-1
               DO 2200 LW = KW+1, NOCCT
                  IF (KLWOP(KW,LW-NISHT) .EQ. 0) THEN
                     KX = ISX(KW)
                     LX = ISX(LW)
                     IF (ISSMO(KX) .EQ. ISSMO(LX)) THEN
                        NWOPX = NWOPX + 1
                        JWOP(1,NWOPX) = KX
                        JWOP(2,NWOPX) = LX
                     END IF
                  END IF
 2200          CONTINUE
 2400       CONTINUE
         END IF
         NWOPH = NWOPX
         NWOPMA = MAX(NWOPH,NWOPMA)
C
 2500    CONTINUE
         IF (NWOPX .GT. MAXWOP) THEN
            WRITE (LUERR,2410) IOPSY,NWOPX,MAXWOP
            WRITE (LUPRI,2410) IOPSY,NWOPX,MAXWOP
            CALL QUIT('SETWOP, TOO MANY ORBITAL ROTATIONS.')
         END IF
 2410 FORMAT(///' SETWOP - Too many orbital rotations in symmetry',I3
     &  //' NWOPX  =',I8/' MAXWOP =',I8
     &  //' To run this calculation you must modify the Dalton source:'
     &   /'   Increase MAXWOP in include/infvar.h; make depend; make')
C
         IF (IPRSTAT.GE.5) WRITE (LUERR,'(/A,I3,A,I3,/A,2I5/)')
     *      ' (SETWOP) orbital rotation operators, abs level',
     *      LVLABS,', symmetry',IOPSY,' NWOPT,NWOPH =',NWOPT,NWOPH
         IF (IPRSTAT.GE.10) WRITE (LUERR,'(A//,(2I10,I5))')
     *      ' I, JWOP(1,I), JWOP(2,I) for I = 1,NWOPH:',
     *      (I,JWOP(1,I),JWOP(2,I),I=1,NWOPH)
C
C        If .NOT.FLAG(11) we only want the totally symmetric operator
C        and we jump to 1299, otherwise we write the operator on
C        LUINDX and continue to next operator symmetry.
C
C        Also if "absorption",
C          LVLABS = 1: inac-act + act-act rotations
C                 = 2: + inact-sec rotations
C                 = 3: + act-sec rotations
C        write to LUINDX.
C
      IF (.NOT.WRINDX) GO TO 1299
C     ------------------------------v   exit loop
         IF (LVLABS.EQ.0) THEN
            WRITE (LUINDX) LAB123,OPLBL(IOPSY)
         ELSE
            IF (NWOPT .EQ. 0) THEN
C           ... 950517-hjaaj: no absorption at this level
C               if no orbital rotations!
               FLAG(50+LVLABS) = .FALSE.
               GO TO 1293
            END IF
            WRITE (LUINDX) LAB123,ABSLBL(LVLABS)
         END IF
         WRITE (LUINDX) NWOPT,NWOPH,IOPSY,NKLWOP,(IDUMMY,I=1,4)
         LENGTH = MAX(IRAT*4,2*NWOPH)
         CALL WRITI  (LUINDX,LENGTH,JWOP)
         LENGTH = MAX(IRAT*4,NKLWOP)
         CALL WRITI  (LUINDX,LENGTH,KLWOP)
 1293    IF (IOPSY.EQ.1 .AND. MAKABS) THEN
             LVLABS = LVLABS + 1
             IF (LVLABS .LE. 3) GO TO 1292
         END IF
 1295 CONTINUE
C
C     recover operator of symmetry JWOPSY
C
      REWIND LUINDX
      CALL MOLLAB(OPLBL(JWOPSY),LUINDX,LUPRI)
      READ (LUINDX) NWOPT, NWOPH, JWOPSY1, NKLWOP
      IF (JWOPSY1 .NE. JWOPSY) CALL QUIT(
     &   'PROGRAM ERROR: JWOPSY on LUINDX does not correspond to label')
      CALL READI  (LUINDX,(2*NWOPH),JWOP)
      CALL READI  (LUINDX,NKLWOP,KLWOP)
C
      CALL GPCLOSE(LUINDX,'KEEP')
 1299 CONTINUE
C
C     Define NWOPDI for /INFDIM/:
C
      NWOPDI = MAX(1,NWOPT)
      FIRST  = .FALSE.
      WRINDX = .FALSE.
      CALL QEXIT('SETWOP')
      RETURN
      END
C  /* Deck raswop */
      SUBROUTINE RASWOP(ISYM,JSYM,NWOPX,JWOP)
C
C     27-Jul-1988 Hans Joergen Aa. Jensen
C     Revised 25. Oct. 2003, in special circumstances
C     (not CAS in RAS2 but CAS in RAS1 and/or RAS3)
C     some non-redundant rotations were not incuded.
C
C     Set up RAS active to active orbital rotations
C     from symmetry ISYM to JSYM.
C
#include <implicit.h>
      INTEGER JWOP(2,*)
C
C Used from common blocks:
C  INFINP : MCTYPE,NELMN1,NELMX1,NELMN3,NELMX3,NACTEL
C  INFORB : IORB(*),NISH(*),NAS1(*),NAS2(*),NAS3(*),NAS1T,NAS2T,NAS3T
C  INFIND : ISSMO
C
#include <maxorb.h>
#include <maxash.h>
#include <infinp.h>
#include <inforb.h>
#include <infind.h>
C
      LOGICAL CAS1, CAS2, CAS3
C
      IF (MCTYPE .NE. 2) RETURN
C
      MINCAS = MAX( 0 , NACTEL - 2*(NAS2T+NAS3T) )
      MAXCAS = MIN( 2*NAS1T , NACTEL )
      IF ( NELMN1 .LE. MINCAS .AND. NELMX1 .GE. MAXCAS) THEN
         CAS1 = .TRUE.
      ELSE
         CAS1 = .FALSE.
      END IF
C
      MINCAS = MAX( 0 , NACTEL - 2*(NAS1T+NAS3T) )
      MAXCAS = MIN( 2*NAS2T , NACTEL )
      NELMN2 = MAX( 0, NACTEL - NELMX1 - NELMX3)
      NELMX2 =         NACTEL - NELMN1 - NELMN3
      IF ( NELMN2 .LE. MINCAS .AND. NELMX2 .GE. MAXCAS) THEN
         CAS2 = .TRUE.
      ELSE
         CAS2 = .FALSE.
      END IF
C
      MINCAS = MAX( 0 , NACTEL - 2*(NAS1T+NAS2T) )
      MAXCAS = MIN( 2*NAS3T , NACTEL )
      IF ( NELMN3 .LE. MINCAS .AND. NELMX3 .GE. MAXCAS ) THEN
         CAS3 = .TRUE.
      ELSE
         CAS3 = .FALSE.
      END IF
C
      ISTT1 = IORB(ISYM) + NISH(ISYM) + 1
      IEND1 = ISTT1 + NAS1(ISYM) - 1
      ISTT2 = IEND1 + 1
      IEND2 = ISTT2 + NAS2(ISYM) - 1
      JSTT2 = IORB(JSYM) + NISH(JSYM) + NAS1(JSYM) + 1
      JEND2 = JSTT2 + NAS2(JSYM) - 1
      JSTT3 = JEND2 + 1
      JEND3 = JEND2 + NAS3(JSYM)
C
C     RAS1 to RAS2 and RAS3 rotations
C     (rotations which are redundant because RAS1 or RAS3 actually is
C      a CAS space are removed)
C
      IF ( .NOT. CAS1 .OR. .NOT. CAS2) THEN
C     ... RAS1 to RAS2 rotations, when at least one is not a CAS
         DO 200 IT = ISTT1,IEND1
            IF (NOROT(IT).NE.0) GO TO 200
            DO 100 JT = JSTT2,JEND2
               IF (NOROT(JT).NE.0) GO TO 100
               IF (ISYM.NE.JSYM .OR. ISSMO(IT) .EQ. ISSMO(JT)) THEN
                  NWOPX = NWOPX + 1
                  JWOP(1,NWOPX) = IT
                  JWOP(2,NWOPX) = JT
               END IF
  100       CONTINUE
  200    CONTINUE
      END IF
C
      IF ( .NOT. CAS1 .OR. .NOT. CAS3 ) THEN
C     ... RAS1 to RAS3 rotations, when at least one is not a CAS
         DO 400 IT = ISTT1,IEND1
            IF (NOROT(IT).NE.0) GO TO 400
            DO 300 JT = JSTT3,JEND3
               IF (NOROT(JT).NE.0) GO TO 300
               IF (ISYM.NE.JSYM .OR. ISSMO(IT) .EQ. ISSMO(JT)) THEN
                  NWOPX = NWOPX + 1
                  JWOP(1,NWOPX) = IT
                  JWOP(2,NWOPX) = JT
               END IF
  300       CONTINUE
  400    CONTINUE
      END IF
C
      IF ( .NOT. CAS2 .OR. .NOT. CAS3 ) THEN
C     ... RAS2 to RAS3 rotations, when at least one is not a CAS
         DO 700 IT = ISTT2,IEND2
            IF (NOROT(IT).NE.0) GO TO 700
            DO 600 JT = JSTT3,JEND3
               IF (NOROT(JT).NE.0) GO TO 600
               IF (ISYM.NE.JSYM .OR. ISSMO(IT) .EQ. ISSMO(JT)) THEN
                  NWOPX = NWOPX + 1
                  JWOP(1,NWOPX) = IT
                  JWOP(2,NWOPX) = JT
               END IF
  600       CONTINUE
  700    CONTINUE
      END IF
C
      RETURN
      END
