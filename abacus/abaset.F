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
#ifdef OLD_LOG
!===========================================================================
!
!=============================================================================
!900105-hjaaj: flag(27) true se MAERKE, to do: transfer flags from Sirius
!900216-hjaaj: call sir_intopen to open luintm, revised check if luintm ok.
!900402-tuh:   H2MO defined according to info from LUSIFC
!900418-hjaaj: NSETUP: call sirini
!901203-hjaaj: SETSIR: FLAG now read from SIRIUS
!910214-hjaaj: SETSIR: IPRSTAT = MAX(IPRIN4,IPRIN6), var altid = 1
!              SETSIR: THRP/THRQ = 1.0D-12, var = 1.0D-10
!              SETSIR: print RAS parameters after reading LUSIFC
!920113-hjaaj: SETSIR: IPRDNS=IPRCIX=...=IPRSIR=MAX(IPRIN4,IPRIN6)
!920226-hjaaj: NSETUP: check if NCONF from SETCI same as NCONF from LUSIFC
!921014-hjaaj: SETSIR: CBISOL added, get SOLVNT information
!921216-hjaaj: SETSIR: read solvent information; check SUPSYM
!===========================================================================
#endif
C  /* Deck setsir */
      SUBROUTINE SETSIR(WORK,LWORK)
C
C     Written by hjaaj 1985
C     Modified for symmetry tuh 060988
C     Solvent added dec 92
C     l.r. aug. 98 hjaaj
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxmom.h>
      PARAMETER (D0 = 0.0D0)
C
      LOGICAL DIPOK, FNDLAB, MCSSYM
      CHARACTER*8 RTNLBL(2)
      DIMENSION WORK(LWORK)
C
C Used from common blocks:
C   GNRINF : RNRESP
C   ABAINF : DIPDER, POLAR,  DOWALK
C   INFINP : ISTATE,ISPIN,NACTEL,POTNUC,LSOLMX,RSOL(3)
C   INFOPT : EMCSCF,EPOT,EMY,EACTIV,ESOLT
C   PAST   : PASDIP, PASORT, PASRES
C   CBISOL : SOLVNT, ...
C
#include <gnrinf.h>
#include <abainf.h>
#include <infinp.h>
#include <infopt.h>
#include <inforb.h>
#include <infvar.h>
#include <infpri.h>
#include <cbdtra.h>
#include <inftap.h>
#include <past.h>
#include <cbisol.h>
#include <cbirhs.h>
C
      CALL QENTER('SETSIR  ')
      DST = SECOND( )
C
      IPRNTR = IPRTRA
C     ... IPRTRA from cbirhs.h, used for IPRTRD in cbdtra.h for
C         DERTRA, transferred to inftra.h for TRACTL in NSETUP
C         below.
      IPRINT = MAX(IPRIN4,IPRIN6,IPRNTR)
C
C     Set /INFORB/ and /INFVAR/  parameters written on LUSIFC.
C     (The rest will be defined in NSETUP)
C
      IF (LUSIFC .LE. 0) CALL GPOPEN(LUSIFC,'SIRIFC',
     &   'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
      IERROR = -1
      REWIND LUSIFC
      CALL MOLLB2(LBSIFC,RTNLBL,LUSIFC,IERROR)
      IF (IERROR .LT. 0) THEN
        WRITE (LUPRI,'(/A/3A)')
     *      ' >>> FATAL ERROR (ABACUS.SETSIR)',
     *      ' >>> label "', LBSIFC,
     *      '" not found on SIRIFC, the SIRIUS interface file.'
        CALL QUIT('"'//LBSIFC//'" label not found on SIRIFC file.')
      END IF
      READ (LUSIFC) EPOT,EMY,EACTIV,EMCSCF,ISTATE,ISPIN,NACTEL,LSYM
      POTNUC = EPOT
      READ (LUSIFC) MISHT,MASHT,MOCCT,MORBT,MBAST,NCONF,NWOPT,NWOPH,
     *            NCDETS,NCMOT,MNASHX,MNASHY,MNORBT,M2ORBT,
     *            NSYM, MULD2H,NRHF,  NFRO,
     *            NISH,NASH,NORB,NBAS,
     *            NELMN1, NELMX1, NELMN3, NELMX3, MCTYPE,
     *            NAS1, NAS2, NAS3
C
      REWIND LUSIFC
      IF (FNDLAB('SIRFLAGS',LUSIFC)) THEN
         READ (LUSIFC) (FLAG(I),I=1,NFLAG)
      ELSE
         WRITE (LUPRI,'(/A/A)')
     &      ' >>> FATAL ERROR (ABACUS.SETSIR)',
     &      ' >>> label "SIRFLAGS" not found on SIRIUS interface file.'
         CALL QUIT('"SIRFLAGS" label not found on SIRIUS interface.')
      END IF
C     disable SUPSYM for ABACUS /921216-hjaaj
      MCSSYM = FLAG(17)
      FLAG(17) = .FALSE.
C
      IF (LBSIFC .EQ. 'CIRESPON') THEN
         FLAG(4) = .TRUE.
         NWOPT  = 0
         NWOPH  = 0
C        ... set CI flag for NSETUP
      END IF
      MWOPT = NWOPT
      MWOPH = NWOPH
C
C   Read in necessary information if we are using external electric field.
C
      REWIND LUSIFC
      IF (FNDLAB('EXTFIELD',LUSIFC)) THEN
         READ (LUSIFC) NFIELD
         NFIEL4 = MAX(4,NFIELD)
         READ (LUSIFC) (EFIELD(I), I=1,NFIEL4)
         READ (LUSIFC) (LFIELD(I), I=1,NFIEL4)
      END IF
C
C      Set parameters for LOGINP
C
C      Note the following:
C        1. FLAG is read in from LUSIFC
C        2. Most of the remaining parameters are only used
C           in SIRIUS, except that DIRFCK determines whether
C           a direct calculations is caried out or not.
C
      DORHF  = .FALSE.
      DOMP2  = .FALSE.
      DOCINO = .FALSE.
      DOCI   = .FALSE.
      DOMC   = .FALSE.
      DORSP  = .FALSE.
      FCVORB = .FALSE.
      LNOROT = .FALSE.
      LMOORD = .FALSE.
      IF (NASHT .GT. 1) THEN
         DIRFCK = .FALSE.
      ELSE
         DIRFCK = DODRCT
      END IF
      CORHOL = .FALSE.
      CORRLX = .FALSE.
      RESPHP = .FALSE.
      JOLSEN = .FALSE.
      ABAIPH = .FALSE.
      INERSI = .FALSE.
      INERSF = .FALSE.
C     
C   Check if "SOLVNT" calculation
C   Read information if so. From WRSIFC in SIRIUS:
C
C     If (solvent) then
C      *) label "SOLVINFO"
C      *) EPSOL,EPSTAT,RSOL(3),LSOLMX
C      *) GRDNRM,POTNUC,EMY,EACTIV,ESOLT,EMCSCF
C      *) ERLM(LM,1), LM = 1,NLMSOL)
C      *) (TRLM(LM), LM = 1,NLMSOL)    where TRLM(i) = ERLM(i,2)
C      *) NSYM, NBAS
C     end if
C
      REWIND LUSIFC
      IF (FNDLAB('SOLVINFO',LUSIFC)) THEN
         IF (.NOT.FLAG(16)) CALL QUIT(
     &      'ABACUS.SETSIR error - "SOLVINFO" but FLAG(16) false.')
         READ (LUSIFC) EPSOL,EPSTAT,RSOL,LSOLMX
         READ (LUSIFC) GRDNRM,POTNUC,EMY,EACTIV,ESOLT,EMCSCF
C        CALL READT (LUSIFC,NLMSOL,ERLM(1,1))
C        CALL READT (LUSIFC,NLMSOL,ERLM(1,2))
C        READ (LUSIFC) NSYM, NBAS
         NLMSOL = (LSOLMX+1) ** 2
C        transfer information to CBISOL
C        NCNTCV is defined in READIN
         SOLVNT = .TRUE.
         LCAVMX = LSOLMX
         LMTOT  = (LCAVMX+1) ** 2
         LMNTOT = (LCAVMX+1)*(LCAVMX+2)*(LCAVMX+3) / 6
         RCAV(1)   = RSOL(1)
         RCAV(2)   = RSOL(2)
         RCAV(3)   = RSOL(3)
         EPDIEL    = EPSOL
         IF (EPSTAT .NE. EPSOL) CALL QUIT(
     &      'ABACUS.SETSIR error: solvent EPSTAT .ne. EPSOL')
         IF (LCAVMX+1 .GT. MXQNM) CALL QUIT(
     &      'ABACUS.SETSIR error: solvent l_max .gt. MXQNM parameter')
      ELSE
         IF (FLAG(16)) CALL QUIT(
     &      'ABACUS.SETSIR error - no "SOLVINFO" but FLAG(16) true.')
         SOLVNT = .FALSE.
      END IF
C
C        Check that no frozen orbitals are used
C
      DO 50 I = 1, NSYM
         IF (NFRO(I).NE.0) THEN
            WRITE (LUPRI,'(/A/)') ' >>> ERROR in SETSIR '//
     *           ' - Frozen orbitals not allowed in ABACUS.'
            WRITE (LUPRI,'(A,8I5)') ' NFRO:', (NFRO(J),J=1,NSYM)
            CALL QUIT(' ERROR in ABACUS.SETSIR - frozen orbitals used.')
         END IF
 50   CONTINUE
      CALL SETORB
      IF (FLAG(17)) THEN
         CALL QUIT('SETSIR ERROR: .SUPSYM not implemented in ABACUS.')
      ELSE
         CALL AVESE0
      END IF
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Information from SIRIFC file in SETSIR',-1)
         WRITE (LUPRI,'(/A,I5)') ' NISHT  ', NISHT
         WRITE (LUPRI,'(A,8I5)') ' NISH() ', (NISH(I),I=1,NSYM)
         WRITE (LUPRI,'(A,8I5)') ' IISH() ', (IISH(I),I=1,NSYM)
         WRITE (LUPRI,'(A,8I5)') ' IIISH()', (IIISH(I),I=1,NSYM)
         WRITE (LUPRI,'(/A,I5)') ' NASHT  ', NASHT
         WRITE (LUPRI,'(A,8I5)') ' NASH() ', (NASH(I),I=1,NSYM)
         WRITE (LUPRI,'(A,8I5)') ' IASH() ', (IASH(I),I=1,NSYM)
         WRITE (LUPRI,'(A,8I5)') ' IIASH()', (IIASH(I),I=1,NSYM)
         WRITE (LUPRI,'(/A,I5)') ' NOCCT  ', NOCCT
         WRITE (LUPRI,'(A,8I5)') ' NOCC() ', (NOCC(I),I=1,NSYM)
         WRITE (LUPRI,'(A,8I5)') ' IOCC() ', (IOCC(I),I=1,NSYM)
         WRITE (LUPRI,'(/A,I5)') ' NORBT  ', NORBT
         WRITE (LUPRI,'(A,8I5)') ' NORB() ', (NORB(I),I=1,NSYM)
         WRITE (LUPRI,'(A,8I5)') ' IORB() ', (IORB(I),I=1,NSYM)
         WRITE (LUPRI,'(A,8I5)') ' IIORB()', (IIORB(I),I=1,NSYM)
         WRITE (LUPRI,'(/A,I5)') ' NNORBT ', NNORBT
         WRITE (LUPRI,'(A,8I5)') ' NNORB()', (NNORB(I),I=1,NSYM)
         WRITE (LUPRI,'(/A,I5)') ' N2ORBT ', N2ORBT
         WRITE (LUPRI,'(A,8I5)') ' N2ORB()', (N2ORB(I),I=1,NSYM)
         WRITE (LUPRI,'(A,8I5)') ' I2ORB()', (I2ORB(I),I=1,NSYM)
         WRITE (LUPRI,'(/A,I5)') ' NBAST  ', NBAST
         WRITE (LUPRI,'(A,8I5)') ' NBAS() ', (NBAS(I),I=1,NSYM)
         WRITE (LUPRI,'(A,8I5)') ' IBAS() ', (IBAS(I),I=1,NSYM)
         WRITE (LUPRI,'(A,8I5)') ' IIBAS()', (IIBAS(I),I=1,NSYM)
         WRITE (LUPRI,'(A,8I5)') ' N2BAS()', (N2BAS(I),I=1,NSYM)
         WRITE (LUPRI,'(/A,8I5)') ' ICMO() ', (ICMO(I),I=1,NSYM)
         WRITE (LUPRI,'(/A,I5)') ' NWOPT  ', NWOPT
         WRITE (LUPRI,'( A,I5)') ' NWOPH  ', NWOPH
         WRITE (LUPRI,'(/A,I5)') ' MCTYPE ', MCTYPE
         WRITE (LUPRI,'(/A,I5)') ' NELMN1 ', NELMN1
         WRITE (LUPRI,'( A,I5)') ' NELMX1 ', NELMX1
         WRITE (LUPRI,'( A,I5)') ' NELMN3 ', NELMN3
         WRITE (LUPRI,'( A,I5)') ' NELMX3 ', NELMX3
         WRITE (LUPRI,'(/A,8I5)') ' NAS1() ', (NAS1(I),I=1,NSYM)
         WRITE (LUPRI,'( A,8I5)') ' NAS2() ', (NAS2(I),I=1,NSYM)
         WRITE (LUPRI,'( A,8I5)') ' NAS3() ', (NAS3(I),I=1,NSYM)
      END IF
C
C
      H2MO = NASHT .GT. 1
C
C     Set some MC parameters
C
      NWOPT  = 0
      IPRSTAT = 1
      DO 100 I = 1,NPFLAG
         P4FLAG(I) = .FALSE.
         P6FLAG(I) = .FALSE.
  100 CONTINUE
      IPRI4 = IPRIN4
      DO 101 I = 1,MIN(IPRI4,NPFLAG)
         P4FLAG(I) = .TRUE.
  101 CONTINUE
      IPRI6 = IPRIN6
      DO 102 I = 1,MIN(IPRI6,NPFLAG)
         P6FLAG(I) = .TRUE.
  102 CONTINUE
      IPRSIR = MAX(IPRIN4,IPRIN6)
      IPRSTAT = IPRSIR
      IPRCIX = IPRSIR
      IPRSIG = IPRSIR
      IPRDNS = IPRSIR
      IPRDIA = IPRSIR
      IPRKAP = IPRSIR
C     Flags are read from SIRIUS
C     Reset some of them:
      FLAG(11) = .TRUE.
      FLAG(14) = .TRUE.
C     FLAG( 4) = (NASHT .EQ. NORBT)
C     ... Full CI case (flag(4) true means CI case)
C     ... 910118-hjaaj: NO, may be a RAS calc. (e.g. MC-SD)
C
C     Set DERTRA control parameters in COMMON /CBDTRA/
C
      IPRTRD = IPRNTR
C
C     call NSETUP to define CI variables, Sirius unit numbers etc.
C
      CALL NSETUP(WORK,LWORK,IPRNTR)
C
      IF (NWOPT.NE.MWOPT .OR. NWOPH.NE.MWOPH) THEN
      IF (MCSSYM) THEN
         WRITE (LUPRI,'(/,A,2(/A,I5,A,I5))')
     &      ' >>> (SETSIR) ERROR: super symmetry was used in SIRIUS',
     &      ' NWOPT from SETUP:',NWOPT,' and NWOPT from SIRIFC:',MWOPT,
     &      ' NWOPH from SETUP:',NWOPH,' and NWOPH from SIRIFC:',MWOPH
         WRITE (LUPRI,'(/A/A)')
     &      ' ABACUS cannot continue because MC orbital gradient',
     &      ' may be non-zero without super symmetry averaging.'
      ELSE
         WRITE(LUPRI,'(/A,2(/A,I5,A,I5))') ' >>> (ABACUS.SETSIR) ERROR',
     &      ' NWOPT from SETUP:',NWOPT,' and NWOPT from SIRIFC:',MWOPT,
     &      ' NWOPH from SETUP:',NWOPH,' and NWOPH from SIRIFC:',MWOPH
         WRITE(LUPRI,'(/A)') 'The JWOP array (with a zero vector):'
         CALL DZERO(WORK,NWOPT)
         CALL PRWOP(WORK,LUPRI)
      END IF
         CALL QUIT('SETSIR: inconsistency found for orbital rotations')
      END IF
C
C     Set FLAG(23) for calculation of active-active block
C     of QX matrix
C
      FLAG(23) = .TRUE.
C
C     END OF SETSIR
C
      DTIM = SECOND() - DST
      IF (IPRDEF .GT. 3) WRITE (LUPRI,1423) DTIM
 1423 FORMAT (/' *** SETSIR-INFO, time in NSETUP: ',F10.2,' seconds.'/)
C
      CALL QEXIT('SETSIR  ')
      RETURN
C
      END
C  /* Deck nsetup */
      SUBROUTINE NSETUP(WORK,LFREE,IPRNTR)
C
C  NSETUP is modified from SIRSET 1-Sep-1988 tuh
C Written by hjaaj 21-Feb-1985
C (based on SETUP from SIRIUS)
C Revised 931209-hjaaj: set /INFTRA/ param. for abarsp
C
C Setup routine for SIRIUS common-blocks and indexing arrays
C used in nuclear derivative program.
C
C SIRIUS SETUP ROUTINE FOR COMMON-BLOCKS AND INDEXING ARRAYS
C
C  PART 1: input/output unit numbers and def. of IROW in SIRINI
C  PART 2: Orbital and basis info, Offset info, Index arrays
C  PART 3: CI parameters
C  PART 4: Open MO integral files
C  PART 5: Miscellaneous tests for inconsistencies
C
C CONVENTIONS:
C
C  PREFIXES N, NN, N2  FOR SIZES OF STRAIGHT, LOWER TRIANGULAR AND
C                      SQUARED BLOCK ARRAYS WITHIN A SYMMETRY.
C  PREFIXES I, II, I2  FOR SIZES OF CUMULATIVE VALUES OF LOWER
C                      SYMMETRIES OF STRAIGHT, LOWER TRIANGULAR AND
C                      SQUARED BLOCK ARRAYS.
C
C  SUFFIX  MA          FOR MAXIMUM VALUE OVER SYMMETRIES
C  SUFFIX  T           FOR SUM OF TRIANGLES OR SQUARES OVER SYMMETRIES
C  SUFFIX  X           WHEN SUM OF ARRAY IS TAKEN BEFORE BEING
C                      TRIANGULARIZED OR SQUARED.
C  SUFFIX  Y           EXAMPLE: NASHY = NASHX*(NASHX + 1)/2
C  SUFFIX  (K)         SYMMETRY NUMBER
C
#include <implicit.h>
#include <priunit.h>
#include <maxash.h>
#include <maxorb.h>
#include <iratdef.h>
#include <mxcent.h>
C -- local constants
      PARAMETER (D1=1.0D0, D2=2.0D0)
#include <dummy.h>
      DIMENSION WORK(LFREE)
#include <litinfo.h>
C
C Used from common blocks:
C   INFINP : ISTATE,LSYM,MCTYPE,?
C   INFORB : NISH(8),...
C   INFVAR : NCONF,NWOPT,NVAR,JWOPSY,JWOP(2,*)
C   INFDIM : NCONMA,NWOPMA,NVARMA,NWOPDI,MAXRL
C   INFOPT : ?
C   INFTRA : ITRLVL,THRP,THRQ
C   INFTAP : LUSIFC
C   priunit.h : LUERR,...
C
#include <abainf.h>
#include <infinp.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
#include <infdim.h>
#include <infopt.h>
#include <inftra.h>
#include <inftap.h>
#include <infpri.h>
#include <inflin.h>
#include <linaba.h>
C
      CHARACTER*8 RTNLBL(2)
C
C PART 1 : ******* I/O buffer sizes, initialize IROW()
C
      CALL SIRINI
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C
C PART 2 : ***** SET UP ORBITAL DATA *********************************
C PART 2 : *****  OFFSETS ********************************************
C PART 2 : **** INDEX ARRAYS *****************************************
C     They are: ISW ISX LOC ISMO ISAO JWOP NSM ICH  IOBTYP IACTYP
C
C
C     CALL SETORB
C     aug.98: has already been done in SETSIR
C
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C
C Define remaining parameters for /INFDIM/:
C
      MAXRL  = 120
C
C
C
      SPIN = (ISPIN-1)/D2
C
C
C   If not CI, call setwop to set-up JWOP pointer for orbital operators
C   and all related information.
C
      JWOPSY = 1
      IF (FLAG(4)) THEN
C     ... this is a CI calculation
         NWOPT  = 0
         NWOPH  = 0
         NWOPMA = 0
         NWOPDI = MAX(1,NWOPT)
      ELSE
         CALL SETWOP(WORK,LFREE)
      ENDIF
C
C PART 3 : **** CI PARAMETERS **************************************
C
      CALL SETCI(NCONFX,NCDETS,LSYM,WORK,LFREE,0)
C
      NVAR   = NCONF  + NWOPT
      NVARH  = NCONF  + NWOPH
      NVARMA = NCONMA + NWOPMA
      NCONDI = MAX(1,NCONF)
C
C     Parameters needed for /LINABA/
C
      ABAHF  = NASHT .LE. 1
      ABACI  = FLAG(4)
C
C     Parameters needed for /INFLIN/
C
      NCONRF = NCONF
      IF (ABAHF) NCONRF = 0
      LSYMRF = LSYM
C
C
C PART 4 : *** Set INFTRA
C          Open MO integral file(s)
C          (FLAG(34) false tells that mo integrals needed)
C          (FLAG(14) true tells that old mo integral file exists)
C          Set ITRLVL to minimum transformation level required.
C
      AO2INTFILE_LABEL = 'AOTWOINT'
      THRP = 1.D-15
      THRQ = 1.D-15
C
C     USEDRC must be set to false, kr-nov.95
C
      USEDRC = .FALSE.
      IPRTRA = IPRNTR
      IF ((NASHT.GT.1 .AND. .NOT.HSROHF) .OR. H2MO) THEN
         FLAG(14) = .TRUE.
         FLAG(34) = .FALSE.
         IF (NEWTRA) THEN
C            ... we need (ia/ (inactive-secondary) Mulliken distributions;
C            they are only calculated for level 5 in "NEWTRA"
C            transformation  /hjaaj June 09
            KTRLVL = 5
            ITRLVL = KTRLVL
         ELSE
            KTRLVL = 4
            ITRLVL = KTRLVL
         END IF
         CALL SIR_INTOPEN
C        ... SIR_INTOPEN returns FLAG(14) false if no mo integral file
C        or if ITRLVL on file too small.
      ELSE IF (ABASOP) THEN
C        ... SOPPA requires full integral transformation
         KTRLVL = 10
         ITRLVL = KTRLVL
         FLAG(14) = .TRUE.
         FLAG(34) = .FALSE.
         CALL SIR_INTOPEN
      ELSE
C        ... Hartree-Fock, MO integral file will normally
C        not exist, but open to define name 'MOTWOINT' instead
C        of e.g. 'fort.13', and to avoid problems if it actually
C        does exist. /980529-hjaaj
         KTRLVL = 4
         ITRLVL = KTRLVL
         FLAG(14) = .TRUE.
         FLAG(34) = .TRUE.
         CALL SIR_INTOPEN
      END IF

      IF (.NOT. FLAG(14)) THEN
      ! ... needed integrals not available
         WRITE (LUPRI,1020)
         KCMO  = 1
         KWTRA = KCMO + NCMOT
         LWTRA = LFREE - KWTRA
         REWIND (LUSIFC)
         CALL MOLLAB(LBSIFC,LUSIFC,LUERR)
         READ (LUSIFC)
         READ (LUSIFC)
         CALL READT (LUSIFC,NCMOT,WORK(KCMO))
         CALL TRACTL(KTRLVL,WORK(KCMO),WORK(KWTRA),LWTRA)
         FLAG(14) = .TRUE.
      END IF
C
C
C PART 5 : **** MISCELLANEOUS ******************************************
C
C     Test for inconsistencies:
C
      NUMERR = 0
      IF (NCONF.LE.0 .AND. IORTO.NE.1) THEN
         WRITE (LUPRI,1610) NCONF
         WRITE (LUERR,1610) NCONF
         NUMERR = NUMERR + 1
      ENDIF
      IF (NCONF .NE. NCONFX) THEN
         WRITE (LUPRI,1620) NCONF,NCONFX
         WRITE (LUERR,1620) NCONF,NCONFX
         NUMERR = NUMERR + 1
      END IF
      IF (NASHT.GT.1 .AND. (FLAG(21) .AND. .NOT. HSROHF)) THEN
         WRITE (LUPRI,1730)
         WRITE (LUERR,1730)
         NUMERR = NUMERR + 1
      ENDIF
      IF (NASHT.EQ.1 .AND. NACTEL.NE.1) THEN
         WRITE (LUPRI,1731)
         WRITE (LUERR,1731)
         NUMERR = NUMERR + 1
      END IF
 1020 FORMAT(/,' NSETUP: MO transformation level too low'
     &        ,' or no MO integral file found.',
     &       /,' NSETUP: generating MO 2-el. integral file.'/)
 1610 FORMAT(/,' NSETUP-ERROR, NCONF =',I10)
 1620 FORMAT(/,' NSETUP-ERROR, NCONF from interface file',I12,
     *       /,'               NCONF from SETCI         ',I12)
 1730 FORMAT(/,' NSETUP-ERROR, NASHT .gt. 0 inconsistent with ',
     *   'RHF calculation')
 1731 FORMAT(/,' NSETUP-ERROR, NASHT .eq. 1 is only implemented',
     *       /,'               for one active electron (doublet)')
C
C ***
C
      IF (NUMERR .GT. 0) THEN
         WRITE (LUPRI,'(///A/)')
     *      ' NSETUP - FATAL ERRORS DETECTED, SEE :'
         CALL SIR_PRTINP(LUPRI,0,0)
         WRITE (LUPRI,'(///A/)')
     *      ' NSETUP - FATAL ERRORS DETECTED, SEE ABOVE'
         CALL QUIT('*** ERROR *** input inconsistent (NSETUP)')
      ENDIF
C
      RETURN
      END
