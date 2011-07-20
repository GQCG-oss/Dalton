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
C
C FILE    : herrdn.F
C
C
#include <single.h>
C
#ifdef OLD_REVISION_LOG
c===========================================================================
crevision 1.7-1.8 2000/04/28 10:27:31  hjj
c - implemented new normalized s,p,d,f,g,h,i GTO's in SPHINP as default
c - reactivated Huckel for Q .le. 36 for spherical GTO's because
c   they are now normalized so we will not get the splitting of
c   e.g. degenerate d-orbitals any more.
c----------------------------
crevision 1.6
c2000/04/27 10:18:19  hjj
c- inserted BN in format statements for reading DOOWN transformation
c- inserted BN and changed Fxx.y to Fxx.0 in FMT[123] for readin
c  GTO exponents and cont.coeff., now integers will always be read correctly.
c----------------------------
crevision 1.5 2000/04/27 09:41:47  hjj
cChanged input format for DOOWN scheme (CRT .eq. 'X') and inserted some tests.
cImproved output of cartesian transformation matrices.
c----------------------------
crevision 1.4 2000/04/27 07:56:52  hjj
ctell users that d, f, ... atomic GTO's are not normalized
c(I have had several questions from people comparing MO coefficients from
cDalton with other programs, or comparing AO integrals with other programs).
c----------------------------
crevision 1.3 2000/04/14 14:09:48  hjj
c\hjaaj dec 99-jan 00:
c1) corrected printout for Huckel basis when DOHUCKEL
c   (the true basis set name was printed instead of "Huckel")
c2) corrected output format from I3 to I4 for .gt. 999 primitives
c3) changed ZMATL format from A6 to A1 to give room for change
c   of CHARGE format from F4.0 to F9.0; this was done so that
c   point charges can be entered with more decimal places.
c4) polishing.
c----------------------------
crevision 1.2 2000/04/12 10:15:14  hjj
cMAXPRD was 22 on most machines except T3D and LINUX, where it was 14.
cNow LINUX also has MAXPRD = 22, following suggestion by Asger Halkier
ce-mail to dalton-wizards 30 Nov 99 /hjaaj
c=============================================================================
c960705-kr:
c-- We use the small component basis for keeping a possible extra STO-3G
c   basis needed when using Hueckel initial guess.
c951115-hjaaj:
cCNTINP: always read NAMN with '(A4)' and never
c   with free format; if free format doesn't work then try formatted;
cZMAT: Read NAMN with '(A4)', allow free input for IRIX;
c   for CRAY formatted input made last field for CHARGE very long;
c   removed old END=204 branches which were reminiscences from when
c   the input was read directly from file.
cBUILDZ: changed format and output of internal coordinates;
c   removed printing of zero cartesian coords. before printing of
c   the correct cartesian coords.
c950428-kr/hjaaj:
cBASOUT,ORBOUT: removed NPRIM,NPRIMD,NORB,NORBD in parameter list to
c ORBOUT; we think they are related to an old organization of
c small component basis (before IQM(*,2) etc.)
c941219-hjaaj:
cBASINP: s/IQMAX/NHTYP/
c950208-vebjorn
cREAD_MOL: NSYMOP is first read as a text to determine if symmetry is to
c        be added (indicated by blanks).
cSYMINP: Has been moved after BASINP, because symmetry has to be added before
c        SYMINP is executed. The full group string CLASS has been added
c        as a parameter to the subroutine.
cSYMADD: Completely new module that adds symmetry on request from user. The
c        full group of the molecule is returned as well.
cBASINP: Has been moved above SYMADD and SYMINP. As a result the updating
c        of the array ISTBNU had to be moved to BASPRO. As this is processing
c        rather than input, the move makes good sense.
cBASPRO: Now takes care of the generation of the symmetry-dependant atoms
c        through the array ISTBNU.
c941119-tsaue:
cREADIN: added extra index to arrays IQM,NBLCK,JCO,NUC,NRC,SEG,ALPHA,
c        CPRIMU ta accomodat both large and small components
cACPORB: implemented readin of exponents only
cSYMGRP: New routine determining irreps, characters etc.
cDBLGRP: new routine for handling double group symmetry
cBASINP: major modifications - check it carefully out !!!
c        may consider replacing NOORBT with MAXIQM wich has more information
cCNTINP: info on angular momentum moved to BASINP
c940207-hjaaj
cREAD_MOL: 'ANG' was not initialized; errors in 'ANG' message format
c931220-hjaaj
cWRONEL: corrected calc. of INAMN and IGTO for LUONEL
c931022-hjaaj
cREAINP: INPTST = IPREAD .gt. 3 (INPTST was not defined)
cexpanded tabs with 'expand'
cREAD_MOL: If ID3 then reset to blank (for abawalk)
cBASINP: changed format for nuclear charge from F4.1 to F4.0
c (such that an entry '6' is read as 6. and not as 0.6 )
cCNTINP: read in old 3F10 format if error reading coordinates
c        changed format to 3F20.0 instead of 3F20.15
c===========================================================================
#endif
C  /* Deck reainp */
      SUBROUTINE REAINP(WORD,RELCAL,TSTINP)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER ( NTABLE = 12, D1 = 1.0D0)
      LOGICAL NEWDEF, RELCAL, TSTINP
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
!
!     Parameter MAXPRD in cbirea.h is default for MAXPRI
#include "cbirea.h"
#include "r12int.h"
#include "qm3.h"
C
      DATA TABLE /'.PRINT ', 'xXXXXXX', '.MAXPRI', '.BIGVEC',
     *            '.ONLYOV', '.SYMTHR', '.CM FUN', '.ZCMVAL',
     *            '.WRTLIN', '.UNCONT', '.ANGSTR', '.R12AUX'/
C
C     ... in qm3.h, for QM3:
      ONLYOV = .FALSE.
C
      NEWDEF = (WORD .EQ. '*READIN' .OR. WORD .EQ. '*MOLBAS')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)',IOSTAT=IOS) WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized for *MOLBAS.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword for *MOLBAS.')
C                 '.PRINT '
    1          CONTINUE
                  IPREAD_ini = IPREAD
                  READ (LUCMD, *) IPREAD
                  IF (IPREAD .EQ. IPREAD_ini) ICHANG = ICHANG - 1
               GO TO 100
C                 'xXXXXXX'
    2          CONTINUE
               GO TO 100
C                 '.MAXPRI'
    3          CONTINUE
                  READ (LUCMD, *,IOSTAT=IOS) MAXPRI
                  IF (IOS.NE.0) THEN
                    CALL QUIT('Error reading .MAXPRI')
                  ENDIF
               GO TO 100
C                 '.BIGVEC'
    4          CONTINUE
                  BIGVC = .TRUE.
               GO TO 100
C                 '.ONLYOV'
    5          CONTINUE
                  ONLYOV= .TRUE.
               GO TO 100
C                 '.SYMTHR'
    6          CONTINUE
                  READ (LUCMD, *,IOSTAT=IOS) TOL_SYMADD
                  IF (IOS.NE.0) THEN
                    CALL QUIT('Error in reading .SYMTHR')
                  ENDIF
               GO TO 100
C              '.CM FUN'
    7          CONTINUE 
                  READ (LUCMD,*,IOSTAT=IOS) LCMMAX, CMSTR, CMEND
                  IF (IOS.NE.0) THEN
                    CALL QUIT('Error reading .CM FUN')
                  ENDIF
                  NCMSTR = NINT(2*CMSTR + 0.01D0)
                  NCMEND = NINT(2*CMEND + 0.01D0)
               GO TO 100
C              '.ZCMVAL'
    8          CONTINUE 
                  READ (LUCMD,*,IOSTAT=IOS) ZCMVAL
                  IF (IOS.NE.0) THEN
                    CALL QUIT('Error reading .ZCMVAL')
                  ENDIF
               GO TO 100
C              '.WRTLIN'  (ach easier to debug intgral input)
    9          CONTINUE
                  WRTLIN= .TRUE.
               GO TO 100
C              '.UNCONT'
   10          CONTINUE
                  UNCONT = .TRUE.
               GO TO 100
C              '.ANGSTR'
   11          CONTINUE
                  ANGS = .TRUE.
               GO TO 100
   12          CONTINUE
C              '.R12AUX'
C              (moved here June 09, because LMULBS used in READIN
C              (renamed from .AUXBAS) /hjaaj
                  LMULBS = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized for *MOLBAS.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt for *MOLBAS.')
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for *MOLBAS:',0)
         IF (IPREAD .NE. IPREAD_ini) THEN
            WRITE (LUPRI,'(A,I5)')
     &         ' Print level in molecule setup (READIN):',IPREAD
         END IF
         IF (ONLYOV) WRITE (LUPRI,'(A)')
     &  ' Stop after calculation of overlap integrals (.ONLYOV) for QM3'
         IF (TOL_SYMADD .NE. 0.0D0) THEN
            WRITE (LUPRI,'(A,G9.2)') ' Symmetry threshold: ',TOL_SYMADD
         END IF
         IF (UNCONT) THEN
            WRITE (LUPRI,'(/A)')
     &         ' Uncontracted basis forced, irrespective of '//
     &         'basis input file.'
         END IF
         IF (MAXPRI .NE. MAXPRD) THEN
            WRITE (LUPRI,'(/A,I5)')
     *         ' Maximum number of primitives per integral block :',
     *         MAXPRI
         END IF
         IF (LMULBS) WRITE (LUPRI,'(/A)')
     &   ' .R12AUX: Reading auxiliary basis sets from .mol file for R12'
         IF (BIGVC) THEN
            WRITE (LUPRI,'(/A)')
     *         ' Primitives from different centers treated '//
     *         'simultaneously.'
            WRITE (LUPRI,'(/A)')
     &       ' Option .BIGVEC not implemented in present version.'
             CALL QUIT('Option .BIGVEC not allowed.')
         END IF
         IF (LCMMAX .GT. -1) THEN
            WRITE (LUPRI,'(/A,I5)')
     &           ' Center of mass basis functions used with l(max) =',
     &           LCMMAX
            WRITE (LUPRI,'(2(A,F4.1))')
     &           ' Sequence starts at n=',CMSTR,' and ends at n=', CMEND
            IF (ZCMVAL .NE. D1)
     &           WRITE (LUPRI,'(A,F3.1)') ' Net charge of ionic core  ',
     &           ZCMVAL
            WRITE(LUPRI,'(/A/A/A/A/A)')
     &         ' Proper reference for these basis functions:'
     &        ,' K. Kaufmann, W. Baumeister, and M. Jungen'
     &        ,' "Universal Gaussian basis sets for an optimum'
     &        ,' representation of Rydberg and continuum wavefunctions"'
     &        ,' J. Phys. B: At. Mol. Opt. Phys. 22 (1989) 2223-2240'

            IF (LCMMAX .GT. 4) THEN
               WRITE (LUPRI,'(/A)')
     &           ' FATAL ERROR: Center of mass basis functions is'//
     &           ' not defined for l(max) .gt. 4'
               CALL QUIT('Error in *MOLBAS input processing')
            END IF
         END IF
         WRITE (LUPRI,'()')
      END IF
      RETURN
      END
C  /* Deck reaini */
      SUBROUTINE REAINI(IPREAD_ini,RELCAL,TSTINP)
C
C     Define default THRS in CCOM
C
C     Initialize /CBIREA/, control variables in ecpinf.h,
C     and parts of qm3.h
C
#include "implicit.h"
#include "maxorb.h"
      LOGICAL RELCAL, TSTINP
#include "infpar.h"
#include "maxaqn.h"
#include "ccom.h"
#include "cbirea.h"
#include "ecpinf.h"
#include "qm3.h"
#include "r12int.h"
Caspg adding numder.h to set NOMOVE
#include "numder.h"
C
      LOGICAL FIRST_CALL
      SAVE    FIRST_CALL
      DATA    FIRST_CALL/.TRUE./

#ifdef VAR_DEBUG
      write(0,*) 'info, REAINI called by node ',mynum
      write(0,*) 'info, REAINI    FIRST_CALL  ',FIRST_CALL
#endif
C
      IF (.NOT. FIRST_CALL) GO TO 9999
C     ... make sure variables are not reinitialized, maybe after reading input
C
C     ccom.h:
C
      THRS = 1.0D-15
      SPHNRM = .TRUE. ! basis functions normalized, also for d, f, ...
C
C     cbirea.h:
C
      NOMOVE = .FALSE.
      IPREAD = IPREAD_ini
      LUMLCL = -1
      MAXPRI = MAXPRD
C     ... parameter MAXPRD is in cbirea.h
      UNCONT = .FALSE.
      BIGVC  = .FALSE.
      DIRAC  = RELCAL
      GENCON = .TRUE.
      INPTST = TSTINP
      IF (TSTINP) IPREAD = MAX(4,IPREAD)
      BASIS  = .FALSE.
      ATOMBA = .FALSE.
!     ATOMDF = .FALSE. ! done each time MOLECULE.INP is opened
      WRTLIN = .FALSE.
      TOL_SYMADD = 0.0D0
      ZCMVAL = 1.0D0
      LCMMAX = -1
C
C     Initialize effective core potential parameters (ecpinf.h)
C
      ECP    = .FALSE.
      NTYECP = 0
C
C     Initialize ISUBSI , NSISY from qm3.h
C     (ISUBSI(:), NSISY are not defined as for QM3, but making them all
C     zero will pass if tests correctly also when QM3 is .false. /hjaaj May09
C
      CALL IZERO(ISUBSI,MXQM3)
      NSISY(0) = 0
C
C     Initialize LMULBS from r12int.h (done here because used in READIN)
C
      LMULBS = .FALSE.
C
      FIRST_CALL = .FALSE.
 9999 CONTINUE
      RETURN
      END ! SUBROUTINE REAINI
C  /* Deck readin */
      SUBROUTINE READIN(WORK,LWORK,HERMIT)
C     Based on
C       READIN:Input processing routine for  M O L E C U L E.
C         Jan Almlof, Stockholm, Dec. 1971.
C         - adaption for ABACUS (T.Helgaker, University of Oslo)
C         - Occams razor by T.Saue, University of Oslo, March 10 - 1993
C           Major surgery includes:
C               - common block CCOM :
C                       - added KHK(MXQN),KCK(MXQN),NHKOFF(MXQN)
C                       - added GTOTYP(MXQN*(MXQN+1)*(MXQN+2)/6)
C               - common block NUCLEI:
C                       - added GAUNUC, GNUEXP(MXCENT)
C                       - added NAMN(MXCENT)
C               - NEW common block FRAME
C                       - added POTNUC, DPNUC(3)
C               - common block CBIREA:
C                       - IPRINT changed to IPREAD
C                       - added INPTST,FAMILY
C               - common block SHELLS:
C                       - added NLRGSH,NSMLSH,NLARGE,NSMALL,NORBS
C                       - added LCLASS(MXSHEL)
C               - common block MOLINP: deleted !!!
C               - common block MOLINC: deleted !!!
C               - common block HRUNIT: deleted !!!
C******************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
C
#include "abainf.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "cbirea.h"
C
      DIMENSION WORK(LWORK)
      LOGICAL   HERMIT,NEWFILE,SET,CMDOPN
      CHARACTER KEYWRD*6

#include "molinp.h"
#include "ccom.h"
#include "nuclei.h"
#include "frame.h"
#include "primit.h"
#include "shells.h"
#include "symmet.h"
#include "chrsgn.h"
#include "aosotr.h"
#include "infpar.h"
#include "huckel.h"
#include <optinf.h>
C
C     A SET variable has been replaced by RDINPC. This allows external control
C     of READIN processing. The processing should be performed once every
C     geometry iteration.
C
      IF (RDINPC) RETURN
      CALL QENTER('READIN')
      RDINPC = .TRUE.
C
      CMDOPN = .FALSE.
      IF (SLAVE) THEN
         LUINFO = -MYNUM
         IPREAD = -1000
      ELSE
         IF (LUCMD .LT. 0) THEN
            CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ','FORMATTED',IDUMMY,
     &           .FALSE.)
            CMDOPN = .TRUE.
C           ... must close LUCMD again before RETURN, when opened here
         END IF
         REWIND(LUCMD)
         READ (LUCMD,'(A6)') KEYWRD
         CALL UPCASE(KEYWRD)
         IF (KEYWRD .EQ. 'INTGRL' .OR. KEYWRD(1:5) .EQ. 'BASIS' .OR. 
     &       KEYWRD .EQ. 'ATOMBA') LUMLCL = LUCMD
         IF (ITERNR .EQ. 0 .AND. ITRNMR .EQ. 0 .AND. IPRE .EQ. 0) THEN
            LUINFO = LUMLCL
         ELSE
            LUINFO = -1
         END IF
      END IF
      IF (IPREAD .GE. 11) CALL TIMER('START ',TIMSTR,TIMEND)
      IF (IPREAD .GT.  0) THEN
         CALL TITLER(
     &   'Output of molecule and basis set information','*',115)
      ELSE
         CALL TITLER(
     &   'Output of molecule and basis set info suppressed','*',118)
      END IF
C
C Only the master in a parallel calculation reads the input file. The
C slaves do not read the input file in RDLINE but use the input lines
C read by the master and transferred in common block MOLINP. /hjaaj
C
      NEWFILE = LUINFO .NE. LUCMD .AND. .NOT. SLAVE
      IF (NEWFILE) THEN
         CALL GPOPEN (LUINFO,'MOLECULE.INP','OLD',' ','FORMATTED',
     &                IDUMMY,.FALSE.)
      END IF

      ATOMDF  = .FALSE.
      NODDYDF = .FALSE.

      IF (.NOT.SLAVE) REWIND(LUINFO)
      NMLINE = 0

      CALL RDLINE(LUINFO)
      CALL UPCASE(MLINE(NMLINE))
      READ (MLINE(NMLINE),'(A6)',IOSTAT=IOS) KEYWRD
      IF (IOS.NE.0) THEN
         WRITE (LUPRI,'(/A,I5/A/A)') 'Error in input line no. ',NMLINE,
     &   MLINE(NMLINE),'Error in reading ... KEYWRD'
         CALL QUIT('Error in reading...KEYWRD')
      ENDIF

#ifdef PRG_DIRAC
C     In Dirac we do the readin the same way whatever KEYWRD is.
#else
      IF (KEYWRD(1:5) .EQ. 'BASIS') THEN
         BASIS = .TRUE.
      ELSE IF (KEYWRD .EQ. 'ATOMBA') THEN
         ATOMBA = .TRUE.
      ELSE IF (KEYWRD .EQ. 'ATOMDF') THEN
         ATOMDF = .TRUE.
         ATOMBA = .TRUE.
      ELSE IF (KEYWRD .EQ. 'NDDYDF') THEN
         NODDYDF = .TRUE.
         ATOMDF  = .TRUE.
         ATOMBA  = .TRUE.
      ELSE IF (KEYWRD .NE. 'INTGRL') THEN
         GOTO 5010
      END IF
#endif
C
C     ********************************************
C     ***** Allocate memory **********************
C     ********************************************
C
      IF (DIRAC) THEN
C        Max number of basis sets : large, small, fit LL dens, fit SS dens,
C        Huckel
         KSETS = 5
         ISETHUCKEL=5
      ELSE
C        Max number of basis sets : regular, Huckel, fit
         IF (DOHUCKEL) THEN
            KSETS = 3
            ISETHUCKEL=2
         ELSE
            KSETS = 2
            ISETHUCKEL = -100
         END IF
      ENDIF
      IF (KSETS .GT. MXBSETS_TOT) THEN
         WRITE(LUPRI,'(//A/A,I5,A,I5)')
     &   'FATAL ERROR: Parameter MXBSETS_TOT in mxbsets.h is too small',
     &   '             need',KSETS,' but MXBSETS_TOT is',MXBSETS_TOT
         CALL QUIT('Parameter MXBSETS_TOT in mxbsets.h is too small')
      END IF
C
      IF (.NOT. QM3) THEN
        KATOM  = MXCENT
      ELSE
        KATOM  = MXCENT_QM
      END IF
      KANG   = 2*MXQN+1 ! max ang momentum needed is 2 * max(L) + 1
      KAOVEC = MXAOVC
      KBLOCK = KATOM*KAOVEC
      KPRIM  = MAXPRI
      KCMAT  = KPRIM*KPRIM
      KNONT  = 1
      KIQM   = KNONT  +  KATOM
      KBLCK  = KIQM   +  KATOM*KSETS
      KJCO   = KBLCK  +  KATOM*KSETS
      KNUC   = KJCO   +  KATOM*KANG*KSETS
      KNRC   = KNUC   +  KBLOCK*KSETS
      KSEG   = KNRC   +  KBLOCK*KSETS
      KLAST  = KSEG   +  KBLOCK*KSETS
C
C     Allocate 100 000 for WORK in READ_MOL and
C     use rest of memory for blocks
C     /sep 05 hjaaj
C
      LWRK   = LWORK  -  KLAST - 100000
      L_PER_BLCK = KSETS*(KPRIM + 2*KCMAT)
      MXBLCK = MAX( 30, LWRK / L_PER_BLCK )
      MXBLCK = MIN( MXBLCK, KBLOCK)
C
      KA     = KLAST
      KC     = KA     +  MXBLCK*KPRIM*KSETS
      KCC    = KC     +  MXBLCK*KCMAT*KSETS
      KLAST  = KCC    +  MXBLCK*KCMAT*KSETS
C     WRITE(LUPRI,'(A,2I16)') 'KLAST, LWORK',KLAST,LWORK
      IF(KLAST.GT.LWORK) GOTO 5020
      LWRK   = LWORK  -  KLAST + 1
C
C     ********************************************
C     ***** Read input file and process data *****
C     ********************************************
C
      CALL READ_MOL(LUINFO,WORK(KLAST),LWRK,WORK(KNONT),WORK(KIQM),
     &            WORK(KBLCK),WORK(KJCO),WORK(KNUC),
     &            WORK(KNRC),WORK(KSEG),WORK(KA),
     &            WORK(KC),WORK(KCC),
     &            KATOM,KANG,KSETS,KBLOCK,KPRIM,KAOVEC,MXBLCK,HERMIT)
C     CALL RDLINE(LUINFO)
C     READ (MLINE(NMLINE),'(A6)') KEYWRD
C     IF (KEYWRD .NE. 'FINISH') GOTO 5030
      IF (IPREAD .GT. 10) CALL REAPRI
      IF (NEWFILE) CALL GPCLOSE(LUINFO,'KEEP')
C
      IF (NODDYDF) CALL DENFIT_NODDY_INT
C
      IF (CMDOPN) CALL GPCLOSE(LUCMD,'KEEP')
      CALL QEXIT('READIN')
      RETURN
C
C Error messages:
C
 5010 CONTINUE
      CALL QUIT('Wrong keyword. Expecting INTGRL or BASIS')
 5020 CONTINUE
      CALL STOPIT('READIN','READ_MOL',KLAST,LWORK)
C5030 CONTINUE
      CALL QUIT('Wrong keyword. Expecting FINISH')
      END ! SUBROUITNE READIN
C  /* Deck read_mol */
      SUBROUTINE READ_MOL(LUINFO,WORK,LWORK,
     &           NONT,IQM,NBLCK,JCO,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,
     &           KATOM,KANG,KSETS,KBLOCK,KPRIM,KAOVEC,MXBLCK,HERMIT)
C******************************************************************************
C
C     Reads and processes input file of molecular data
C
C       **** Temporary variables used in READIN ****
C
C       KATOM   - max number of atomic types
C       KANG    - max number of different angular momenta
C       KBLOCK  - max number of AO-blocks
C       KPRIM   - max number of primitives
C
C      (Large/small)
C       NONT            (KATOM)         - number of symmetry independent
C                                         centers for this atomic type
C       IQM             (KATOM)         - highest L-value
C       NBLCK           (KATOM)         - number of AO-blocks for a given
C                                         atomic type
C       JCO             (KANG,KATOM)    - number of AO-blocks for a given
C                                         atomic type and L-value
C       NUC             (KBLOCK)        - number of uncontracted functions
C                                         in a given AO-block
C       NRC             (KBLOCK)        - number of contracted functions
C                                         in a given AO-block
C       SEG             (KBLOCK)        - TRUE: segmented contraction
C       ALPHA           (KPRIM,KBLOCK)  - exponents
C       CPRIM           (KPRIM,KPRIM,KBLOCK)  - normalized contraction
C                                               coefficients
C       CPRIMU          (KPRIM,KPRIM,KBLOCK)  - contraction coefficients
C
C     LOCAL VARIABLES:
C       IPRIMA - runs over primitives for a given atomic type
C       ISHELL - runs over all shells
C       IPRIM  - runs over all primitives
C       NUCIND - runs over all symmetry independent centers
C******************************************************************************
#include "implicit.h"
#include "iratdef.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include <codata.h>
      PARAMETER (D0 = 0.0D0, THRMIN = 1.D-15)
#include "clsfmm.h"
#include "cbirea.h"
#include "huckel.h"
      LOGICAL   HERMIT,DOOWN,ZMAT
      LOGICAL   SEG(KBLOCK,KSETS)
      DIMENSION NONT(KATOM),IQM(KATOM,KSETS),NBLCK(KATOM,KSETS),
     &          JCO(KANG,KATOM,KSETS),ISGEN(KATOM),
     &          NUC(KBLOCK,KSETS),NRC(KBLOCK,KSETS),
     &          ALPHA(KPRIM,KBLOCK,KSETS),
     &          CPRIM(KPRIM,KPRIM,KBLOCK,KSETS),
     &          CPRIMU(KPRIM,KPRIM,KBLOCK,KSETS), WORK(LWORK)
C
C       Local variables
C
      CHARACTER*80 BASREF(10,KATOM,3)
      CHARACTER*1 KASYM(3,3),ID3,SPDCAR,CRT
      CHARACTER*72 TTITLE(2)
      CHARACTER*2 SYMTXT
      CHARACTER*15 SYMGROUP
      DIMENSION IFXYZ(3)
#include <cbihr1.h>
      LOGICAL ANG, ADDSYM
#include "molinp.h"
#include "ccom.h"
#include "nuclei.h"
#include "primit.h"
#include "shells.h"
#include "symmet.h"
#include "aosotr.h"
#include "orgcom.h"
#include "infpar.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
      PARAMETER (LUCME = -1)
#else
#include "gnrinf.h"
#include "inftap.h"
#endif
#include "dgroup.h"

#include <numder.h>
#include <r12int.h>
#include <qm3.h>
#include <pcmlog.h>
      CALL QENTER('READ_MOL')
      ZMAT = .FALSE.
C     Initialize auxiliary basis set variables (WK/UniKA/04-11-2002).
      LAUXBS = .FALSE.
      MBSMAX = 4
      NMULBS = 0
C     Initialize ISUBSY for QM3 info to all atoms QM atoms
      CALL IZERO(ISUBSY,MXQM3)
C
C***********************************************************************
C     Card 1: Keyword (INTGRL/BASIS) - read in READIN
C***********************************************************************
C***********************************************************************
C     Card 2: If basis set library, read additional card containing basis set
C***********************************************************************
      DO J=1,MXMULB
         MULNAM(J) = ' '
      END DO
      IF (BASIS) THEN
         CALL RDLINE(LUINFO)
C        Read basis set name(s), set NMULBS (WK/04-11-2002).
         CALL READ_BASMUL(MLINE(NMLINE))
      END IF
C     Enforce computation of Cartesian second moments (WK/UniKA/04-11-2002).
C     Sep 09/hjaaj: moved to HERINI
C     IF (LMULBS) CALL SET_CARMOM_2
C*****************************************************************************
C     Card 2-3: Title cards
C*****************************************************************************
      CALL RDLINE(LUINFO)
      READ (MLINE(NMLINE),'(A72)',IOSTAT=IOS) TTITLE(1)
      CALL RDLINE(LUINFO)
      READ (MLINE(NMLINE),'(A72)',IOSTAT=IOS) TTITLE(2)
      IF (IPREAD .GT. 0) THEN
         WRITE (LUPRI,'(/4X,A/4X,A/2A/2A/4X,A)')
     &   'The two title cards from your ".mol" input:', SEPARATOR(1:72),
     &   ' 1: ',TTITLE(1),' 2: ',TTITLE(2),SEPARATOR(1:72)
         IF (LUCME.GT.0) WRITE (LUCME,'(2X,A)')
     &      SEPARATOR(1:72),TTITLE(1),TTITLE(2),SEPARATOR(1:72)
      END IF
C***********************************************************************
C     Card 4:
C       CRT     - flag for spherical harmonics or "your own scheme"
C       NONTYP  - number of atomic types
C       SYMTXT  - Symmetry operations read as a text to test if
C                 symmetry is to be added automatically
C       KCHARG  - Molecular charge
C       NSYMOP  - number of independent twofold symmetry operations of
C                 point group to be used
C       KASYM   - 3x3 character array specifying the basic symmetry operations;
C                 specify axis that are reversed under the operation
C       THRS    - integral threshold
C       ADDSYM  - flag for automatic determination of symmetry
C***********************************************************************
      CALL RDLINE(LUINFO)
      CALL UPCASE(MLINE(NMLINE))
      CALL LINE4(MLINE(NMLINE),NONTYP,NSYMOP,CRT,KCHARG,THRS,ADDSYM,
     &           KASYM,ID3,DOCART,DOOWN)
C --> Angstrom or Bohr?
      NMLAU = NMLINE
C     .. NMLAU indentifies input line with Angstrom or Bohr indicator.
      IF (ID3 .NE. ' ' .or. ANGS ) THEN
        ANG = .TRUE.
        IF (IPREAD .GT. 0)
     &      WRITE (LUPRI,'(/2X,A/10X,A,F11.8,A2)')
     &          'Coordinates are entered in Angstrom'//
     &          ' and converted to atomic units.',
     &          '- Conversion factor : 1 bohr =',XTANG,' A'
      ELSE
         ANG = .FALSE.
      ENDIF
C
      IF (NONTYP.EQ.0) THEN
        WRITE (LUPRI,'(/A/A//A//I3,1X,A)')
     &    ' You have specified a molecule with zero atoms,',
     &    ' thus all answers to all your input are zero!',
     &    ' (or you made an input error in the .mol file)',
     &    NMLINE,MLINE(NMLINE)
        CALL QUIT('No atoms according to .mol input!')
      ELSE IF (NONTYP.LT.0) THEN
        WRITE (LUPRI,'(/A,I6//I3,1X,A)')
     &    ' >>> READ_MOL error, no. of atomic types negative:',NONTYP,
     &    NMLINE,MLINE(NMLINE)
        CALL QUIT('Negative number of atoms according to .mol input')
      ELSE IF (NONTYP.GT.KATOM) THEN
        WRITE (LUPRI,'(/A,I6/A,I6//I3,1X,A)')
     &   ' >>> READ_MOL error, no. of atomic types    ',NONTYP,
     &   '   is greater than current maximum number ',KATOM,
     &    NMLINE,MLINE(NMLINE)
        CALL QUIT('Too many atomic types in .mol input')
      END IF
C***********************************************************************
C Check integral threshold: THRS must be .gt. 0
C       ( log(thrs*thrs) taken in oneint) (851005-hjaaj)
C***********************************************************************
      IF (THRS .LT. THRMIN) THEN
         IF (THRS .GE. 0.0D0) THEN
C        ... THRS .lt. 0.0D0 is code from LINE4 that
C            user has not specified THRS
            WRITE (LUPRI,'(1P,2(/2X,A,D12.2))')
     &'INFO from READIN: Threshold for discarding integrals was',THRS,
     &'INFO from READIN: Threshold is reset to minimum value   ',THRMIN
         END IF
         THRS = THRMIN
C        note that her2drv.F has explicit THRESH = MAX(THRS,1.0D-15)
C        - hjaaj, aug99
      ENDIF
C***********************************************************************
C Process Cartesian to spherical or your own transformation
C     note: SPHINP will read additional input lines if DOOWN true
C*****************************************************************************
      CALL CRTEST(DOCART)
      CALL SPHINP(LUINFO,WORK,LWORK,DOOWN,MXSPD)
C
C     ************************************************
C     ***** Read orbital and geometry input data *****
C     ************************************************
C
      SPHNRM = .NOT.(DOCART .OR. DOOWN)
C     ... SPHNRM only true if sperical, i.e. the only case
C         where we know the 5 d's have same norm. Used for Huckel.
      CALL BASINP(LUINFO,WORK,LWORK,NONTYP,NSYMOP,NONT,
     &            IQM,NBLCK,JCO,BASREF,
     &            NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,ISGEN,
     &            KATOM,KANG,KSETS,KBLOCK,KPRIM,KAOVEC,
     &            ANG,ZMAT,MXSPD,ADDSYM,MXBLCK)
C
C     *************************************************
C     ***** Process ZMAT coordinates if requested *****
C     *************************************************
C
      IF (ZMAT) CALL BUILDZ(IPREAD,NONTYP,NSYMOP,NONT,IQM,NBLCK,JCO,
     &               NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,ISGEN,KATOM,
     &               KANG,KBLOCK,KPRIM,KAOVEC,DOOWN,MXSPD,BASREF)
C
C     *************************************************
C     ***** Add symmetry to molecule if requested *****
C     ***** Returns SYMGROUP, the full symmetry   *****
C     *************************************************
C

      LINEAR = .FALSE.
      SYMGROUP = 'N/A'
      IF (ADDSYM) THEN
         IF (IPREAD .GE. 3) THEN
            CALL HEADER('Copy of input to READIN before ADDSYM',0)
            WRITE (LUPRI,'(A)') (MLINE(I), I = 1, NMLINE)
            WRITE (LUPRI,'(/)')
         END IF
         CALL SYMADD(WORK,LWORK,NONTYP,NSYMOP,NONT,KATOM,KASYM,
     &               SYMGROUP,TOL_SYMADD,IPREAD)
         IF (SYMGROUP(3:4).EQ.'oo') LINEAR = .TRUE.
      END IF
C
C     ***************************************
C     ***** Process symmetry input data *****
C     ***************************************
C
      CALL SYMINP(NSYMOP,KASYM,IFXYZ, SYMGROUP)
C
C     ***************************************************
C     ***** Process orbital and geometry input data *****
C     ****************************************************
C
      CALL BASPRO(WORK,LWORK,NONTYP,NSYMOP,NONT,IQM,NBLCK,JCO,NUC,NRC,
     &            SEG,ALPHA,CPRIM,CPRIMU,KATOM,KANG,KSETS,
     &            KBLOCK,KPRIM,DOOWN)

#ifdef BUILD_GEN1INT
C...  added by Bin Gao, Feb. 16, 2011
C...  sets the AO blocks for Gen1Int
      call gen1int_block
#endif
C
C     ***********************************
C     ***** Process PCM speres data *****
C     ***********************************
C
      IF(PCM) CALL PCMSPHGEN
C
C     Dec 04 hjaaj: TODO
C     EWMO is only implemented for case with no
C     symmetry adapted orbitals (yet)
C     This is checked with NUCDEP .NE. NUCIND.
C     (the Hamiltonian matrix must be generated BEFORE
C      the symmetry transformation, thus it requires a
C      special output routine in her1int.F, or a backtransformation
C      from SO basis to AO basis in huckel.F).
C
      IF (NUCDEP .NE. NUCIND) THEN
         EWMO = .FALSE.
      END IF
C
C     Determine center of mass
C
      KGEOM = 1
      KMASS = KGEOM + 3*NUCDEP
      KNAT  = KMASS +   NUCDEP
      KNUMIS= KNAT  +  (NUCDEP + 1)/IRAT
      KLAST = KNUMIS+  (NUCDEP + 1)/IRAT
      IF (KLAST .GT. LWORK) CALL STOPIT('READ_MOL','CMMASS',KLAST,
     &                                   LWORK)
      CALL CMMASS(WORK(KGEOM),WORK(KMASS),WORK(KNAT),WORK(KNUMIS),
     &            IPREAD-1)
      IF (ADDSYM) THEN
         CALL DCOPY(3,CMXYZ,1,DIPORG,1)
         CALL DCOPY(3,CMXYZ,1,GAGORG,1)
      ELSE
chj may 09: now done here because READIN now before HERINP
         CALL DZERO(DIPORG,3)
         CALL DZERO(GAGORG,3)
      END IF
      CALL DCOPY(3,CMXYZ,1,CAVORG,1)
C
C     *************************************************
C     ***** Print orbital and geometry input data *****
C     *************************************************
C
      CALL BASOUT(WORK,LWORK,NONTYP,NONT,IQM,NBLCK,JCO,
     &            BASREF,
     &            NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,KATOM,KANG,
     &            KSETS,KBLOCK,KPRIM,KASYM,NSYMOP,HERMIT)
C
      IF (IPREAD .GE. 2) THEN
         CALL HEADER('Copy of input to READIN',0)
         IF (ADDSYM) WRITE(LUPRI,'(A/)')
     &   ' - as modified by symmetry addition module:'
         WRITE (LUPRI,'(80A1 )') ('-',I = 1, 80)
         WRITE (LUPRI,'(A)') (MLINE(I), I = 1, NMLINE)
         WRITE (LUPRI,'(80A1/)') ('-',I = 1, 80)
      END IF
C
      THRCLS = THRS
C
C     ****************************************
C     ***** Output on LUONEL *****************
C     ****************************************
C
      IF (HERMIT) CALL WRONEL(TTITLE,NONTYP,NONT,IQM,IFXYZ,KATOM,JCO,
     &                        KANG)
C
      CALL QEXIT('READ_MOL')
      RETURN
      END ! SUBROUTINE READ_MOL
C  /* Deck baspar */
      SUBROUTINE BASPAR(NSYMOP)
C***********************************************************************
C
C       Set various basis parameters:
C               KHK(J) - number of spherical (cartesian) components for given J
C               KCK(J) - number of Cartesian components for given J
C                      - tabulate incomplete Gamma function
C                      - determine Cartesian powers
C               NHKOFF(J) - offset for components in list of l-functions
C               MAXOPR - maximum number of operations to loop over
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "pi.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "cbirea.h"
C
#include "ccom.h"
#include "nuclei.h"
#include "symmet.h"
#include "aosotr.h"
#include "ibtfun.h"
C
C     ***** Tabulate Incomplete Gamma Function *****
C
      JMAX = MAX(4,4*(NHTYP - 1) + 2)
      CALL GAMTAB(JMAX)
C
C     ***** Cartesian powers *****
C
      CALL CARPOW
C
C     ***** NHKOFF() *****
C     - offset for components in list of l-functions
C
      IOFF = 0
      DO 100 I = 1, NHTYP
         NHKOFF(I) = IOFF
         IOFF = IOFF + KHK(I)
  100 CONTINUE
C
C     ***** MAXOPR *****
C     - maximum number of operations necessary to loop over
C
      MAXOPR = 1
      IF (NSYMOP .GT. 0) THEN
         MAXLO = 7
         DO 200 I = 1, NUCIND
            MAXLO = IBTAND(ISTBNU(I),MAXLO)
  200    CONTINUE
         II = 1
         DO 300 I = 1, NSYMOP
            II = 2*II
            IF(IBTAND(II,MAXLO).EQ.0) MAXOPR = II
  300    CONTINUE
      END IF
      MAXOPR = MAXOPR - 1
C
      END ! SUBROUTINE BASPAR
C  /* Deck basinp */
      SUBROUTINE BASINP(LUINFO,WORK,LWORK,NONTYP,NSYMOP,NONT,
     &                  IQM,NBLCK,JCO,BASREF,
     &                  NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,ISGEN,
     &                  KATOM,KANG,KSETS,KBLOCK,KPRIM,KAOVEC,
     &                  ANG,ZMAT_INPUT,MXSPD,ADDSYM,MXBLCK)
C***********************************************************************
C
C       Read orbital and geometry input data
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0)
      PARAMETER (DSM = 1.0D-30)
#include "cbirea.h"
#include "symmet.h"
      DIMENSION WORK(LWORK)
      LOGICAL   SEG,ANG,NOORBTS,NOFITFS,ZMAT_INPUT,CNTBAS, MM_ATOMTYPE
      LOGICAL   ADDSYM
      INTEGER   INT_TEMP
      DIMENSION NONT(KATOM),IQM(KATOM,KSETS),JCO(KANG,KATOM,KSETS),
     &          NBLCK(KATOM,KSETS),NUC(KBLOCK,KSETS),NRC(KBLOCK,KSETS),
     &          SEG(KBLOCK,KSETS), ALPHA(KPRIM,KBLOCK,KSETS),
     &          CPRIM( KPRIM,KPRIM,KBLOCK,KSETS),
     &          CPRIMU(KPRIM,KPRIM,KBLOCK,KSETS), ISGEN(KBLOCK)
      CHARACTER SPDCAR*1,BSET*5,ZMATL*4, BKPLIN*80, TSTCHA*2
      CHARACTER*80 BASNAM, AUXNAM, BASNAMSAVE
C...  added by Bin Gao, for effective core potentials, Dec. 6, 2010
      character*80 ECPNAME
C...
      CHARACTER*80 BASREF(10,KATOM,KSETS)
      DIMENSION IBLOCK(MXBSETS_TOT)
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "molinp.h"
#include "ccom.h"
#include "nuclei.h"
#include "huckel.h"
#include "cbisol.h"
#include "aosotr.h"
#include "ecpinf.h"
#include "qm3.h"
#include "infpar.h"
C
C     MBSIP = identifier for Previous basis-set identifier MBSI (WK/UniKA/04-11-2002).
      MBSIP     = 0
C
      NUCIND    = 0
      NCTOT     = 0
      NHTYP     = 0
      NTYECP    = 0
      NHTYPF    = 0
      IBLOCK(1:KSETS) = 1
      CALL IZERO(IQM,KATOM*KSETS)
      CNTBAS = .FALSE.
#ifdef PRG_DIRAC
      BSET_is_GENCON(-mxbsets:mxbsets) = .FALSE.
#endif

C
C     Initialize NSETS and BSET for Dalton,
C     if (DIRAC) these will be modified below:
C     NSETS = 1 for Large, = 2 for "Large" and "Small"
C     and BSET will change to "Large" or "Small".
C
      NSETS = 1
      IF (ATOMDF) THEN
         NSETS = 3
      ELSE IF (DOHUCKEL .OR. DIRAC) THEN
         NSETS = 2
      END IF
C     If (DOHUCKEL) use second set for atomic orbitals for Huckel guess
      BSET = 'Basis'
      MBSI = 1
      IF (NSETS .GT. KSETS) THEN
         CALL QUIT('NSETS .gt. KSETS error')
      END IF
      if (ipread .gt. 10) write(lupri,*) 'setting NSETS =',NSETS
C
C     Run over atomic centers
C     =======================
C
      NONTYP_QM = 0
      MM_ATOMTYPE = .FALSE.
      DO 100 I = 1, NONTYP
         NOORBTS = .TRUE. ! NO ORBiTalS
         NOFITFS = .TRUE. ! NO FIT FunctionS
         INUC = NCTOT
C
C***********************************************************************
C        Card 6: ATOMIC DATA
C         Q     - nuclear charge
C         NONT  - number of symmetry independent atoms
C         IQM   - maximum angular quantum number (s=1,p=2,d=3 etc.)
C         NHTYP - maximum angular quantum number for ALL orbitals
C         JCO   - number of AO-blocks for each l-value
C***********************************************************************
C
         CALL RDLINE(LUINFO)
         BKPLIN = MLINE(NMLINE)
         CALL UPCASE(BKPLIN)
C
C        Test if Z-matrix for first line (I.eq.1):
C
         IF (I .EQ. 1) THEN
           READ (BKPLIN,'(A4)',IOSTAT=IOS) ZMATL
           IF (IOS.NE.0) CALL QUIT('Error read READ..ZMATL')
           IF (ZMATL .EQ. 'ZMAT') THEN
             ZMAT_INPUT = .TRUE.
             WRITE (LUPRI,'(/A)')
     &        'Atoms and their geometry are read in Z-matrix format.'

             IF (DIRAC) CALL QUIT
     &          ('READ_MOL: Z-matrix is not implemented for DIRAC')
             IF (.NOT. BASIS) CALL QUIT
     &          ('READ_MOL: Z-matrix is only implemented with BASIS')

             CALL RDLINE(LUINFO)
           END IF
         END IF
C
C        MM atomtype or QM atomtype ?
C
         IPOS = INDEX(BKPLIN,'MM')
         IF (IPOS .GT. 0) THEN
            MM_ATOMTYPE = .TRUE.
         ELSE
C           Arnfinn, nov -08: add 1 to NONTYP_QM if not MM;
C           counting the numbers of QM type atoms
            NONTYP_QM = NONTYP_QM + 1
C           HJAaJ dec 08: NONTYP_QM sent to slaves when parallel,
C           slaves will then only read transferred MLINE(1:NMLINE) correctly
C           if all QM atoms before any MM atom
            IF (MM_ATOMTYPE) THEN
               WRITE (LUPRI,'(/A,I6/A/A)')
     &            'ERROR in .mol input line no.',NMLINE,MLINE(NMLINE),
     &            'All QM atoms must come before any MM atom.'
               CALL QUIT(
     &            'READ_MOL: All QM atoms must come before any MM atom')
            END IF
         END IF
C
C        New or old input style ?
C
         IF (INDEX(BKPLIN,'CHA') .NE. 0) THEN
C           ... New input style
            CALL LINE5R(MLINE(NMLINE),Q,NONT(I),MBSI,IQM(I,1),
     &                  JCO(1,I,1),KANG,BASIS,ATOMBA,LMULBS,
     &                  BASNAMSAVE,RADIUS_PCM, ALPHA_PCM)
            Q_EXPONENT = D0 ! hjaaj: should be introuduced in LINE5R; GNU=
         ELSE
C           ... Old input style
            Q_EXPONENT = D0
C
C        Read/process data on centers/blocks (old input style)
C        ===================================
C
         IF (DIRAC) THEN
            READ (MLINE(NMLINE),'(BN,F10.0,2I5)') Q,NONT(I),NSETS
            IF (LMULBS) THEN
               CALL QUIT('READ_MOL: LMULBS not implemented for Dirac')
            ELSE
               MBSI = 1
            END IF
         ELSE IF (ZMATL .EQ. 'ZMAT') THEN
C           Z-matrix not implemented with LMULBS (WK/UniKA/04-11-2002).
            IF (LMULBS) THEN
               CALL QUIT(
     &            'READ_MOL: Z-matrix not implemented with LMULBS')
            ELSE
               MBSI = 1
            END IF
            CALL ZMAT(NONT,KATOM)
            Q = CHARGE(NCTOT)
            ISOTOP(NCTOT) = 1
C           ... "Isotope=#" not implemented for ZMAT yet,
C               thus always use isotope no. 1 /hjaaj feb 2003
         ELSE IF (ATOMBA .OR. BASIS) THEN
C           Read basis-set identifier from MLINE(NMLINE)
C           Put MBSI = 1 otherwise (WK/UniKA/04-11-2002).
            IF (LMULBS) THEN
               READ (MLINE(NMLINE),'(BN,I3,F7.0,I5)',ERR=299) 
     &               MBSI,Q,NONT(I)
            ELSE
               READ (MLINE(NMLINE),'(BN,F10.0,I5)',ERR=299) Q,NONT(I)
               MBSI = 1
            END IF
         ELSE
C ach 
            IF (WRTLIN) THEN 
               WRITE(LUPRI,*) MLINE(NMLINE)
               CALL FLSHFO(LUPRI)
            ENDIF
C
            IF (MM_ATOMTYPE) THEN
                READ (MLINE(NMLINE),'(BN,1X,F9.0,I5)') Q,NONT(I)
            ELSE
C              Read basis set identifier from MLINE(NMLINE)
C              Put MBSI = 1 otherwise (WK/UniKA/04-11-2002).
              IF (LMULBS) THEN
                 READ (MLINE(NMLINE),'(BN,I3,F7.0,I5,24I5)',ERR=299)
     &           MBSI,Q,NONT(I),IQM(I,1),(JCO(J,I,1), J=1, IQM(I,1))
              ELSE
                 IF (MLINE(NMLINE)(1:3) .ne. '   ' .AND.
     &               MLINE(NMLINE)(4:4) .eq. ' ' ) THEN
                     WRITE(LUPRI,*)
     &            'ERROR, multiple basis set, but no .R12AUX specified.'
                     WRITE(LUPRI,*) MLINE(NMLINE)
                     CALL QUIT('Error when reading .mol file')
                 END IF
                 READ (MLINE(NMLINE),'(BN,F10.0,I5,24I5)',ERR=299)
     &                Q,NONT(I),IQM(I,1),(JCO(J,I,1), J=1, IQM(I,1))
                 MBSI = 1
              END IF
            END IF
         END IF
C
C        End of processing of old input format for line "5"
C         
         END IF
C
C        Check Multiple Basis Set Identifier, MBSI (WK/UniKA/04-11-2002).
C
         IF (LMULBS .AND. MBSI .LE. 0) CALL QUIT(
     &   'READIN: Please provide proper basis set identifier, i.e. > 0')
         IF (LMULBS .AND. MBSI .LT. MBSIP) CALL QUIT('READ_MOL: '//
     &   'Please provide multiple basis sets in increasing order')
         MBSIP = MBSI
         NMULBS = MAX(MBSI,NMULBS)
         IF (NMULBS .GT. MXMULB)
     &      CALL QUIT('READ_MOL: Too many multiple basis sets')
C
C        ---> Check for Huckel :
C
Chj aug99: changed test from 36 to 20, because current Huckel implementation
C       (1) gives poor orbital energies for some d-orbitals
C           because of normalization of x2, y2, z2 to three
C           (instead of one)
C       (2) and breaks symmetry between e.g. x2-y2 and xy because of (1)
C
         IF (.NOT.SPHNRM) THEN
C        ... SPHNRM is only true if spherical b.f., i.e. the only case
C            where we know the 5 d's have same norm.
            MXQHUCK = 20
         ELSE
Chj - apr00: this has been fixed with normalized spherical GTO's from SPHCOM
C Manu MXQHUCK set to 103 instead of 86 
            MXQHUCK = 103 
         ENDIF
         IF (DOHUCKEL .AND. NINT(Q) .GT. MXQHUCK) THEN
            NSETS = 1
            IF (ATOMDF) NSETS = 3
            DOHUCKEL = .FALSE.
            NINFO = NINFO + 1
            WRITE (LUPRI,'(/A,I3)')
     &           '  INFO: No Huckel start guess due'//
     &           ' to nuclear charge larger than',MXQHUCK
C
C     Just in case we have already processed some Huckel-orbitals for 
C     some centers, we zero out IQM in order to prevent errors when sorting
C     integral shells
C
           CALL DZERO(IQM(1,ISETHUCKEL),KATOM)
        ENDIF
        QEFF = Q
C       ... QEFF is effective charge and will be reset if ECP
C
C       Read/process data on symmetry independent centers
C       =================================================
C       (if ZMAT then this is done later)
C
         IF (ZMAT_INPUT) GO TO 100
C
         IF (IPREAD .GT. 0 .OR. NONT(I).LE.0) THEN
            WRITE(LUPRI,'(/A,I5/A/A,F10.5/A,I5)')
     &         '  Atomic type no.',I,
     &         '  --------------------',
     &         '  Nuclear charge:',Q,
     &         '  Number of symmetry independent centers:',NONT(I)
            IF (NONT(I).LE.0) THEN
               WRITE(LUPRI,'(/A)') 'INPUT ERROR: '//
     &         'Number of symmetry independent centers invalid.'
               CALL QUIT('Input error in .mol file, see output.')
            ENDIF
         ENDIF
         CALL CNTINP(LUINFO,NONT(I),Q,Q_EXPONENT,NSYMOP,ANG,MBSI,
     &               RADIUS_PCM, ALPHA_PCM)
C        .. MBSI has been added (WK/UniKA/04-11-2002).
C
         IF (MBSI .GT. 1) THEN
!           .. no Huckel and density fitting bases for duplicate atomic
!              centers for R12 auxiliary basis
            NSETS_NOW = 1
         ELSE
            NSETS_NOW = NSETS
         END IF
         IF (IPREAD .GT. 0) THEN
            IF (MBSI.gt.1) THEN
               WRITE(LUPRI,'(A,I5)')
     &         '  This is an auxiliary basis for R12, MBSI:',MBSI
            ELSE
               WRITE(LUPRI,'(A,I5)')
     &         '  Number of basis sets to read;',NSETS_NOW
            END IF
         END IF
         DO 300 J = 1,NSETS_NOW
          IF (J.EQ.1 .OR. (J.EQ.ISETHUCKEL) 
     &               .OR. (J.EQ.2 .AND. DIRAC) 
     &               .OR. (J.EQ.3 .AND. ATOMDF) ) THEN
C
C         Construct Rydberg/continuum basis set for X center, if requested
C         ================================================================
C         Center of mass centers are designated with 'x' or 'X'
C
            IF ((LCMMAX .GT. -1) .AND. (NONT(I) .EQ. 1) .AND.
     &          (INDEX('   X   ',NAMN(NCTOT)) .GT. 0 .OR.
     &           INDEX('   x   ',NAMN(NCTOT)) .GT. 0)) THEN
               IF (J .EQ. 1) THEN
                  IQM(I,J) = LCMMAX + 1
                  CALL CMBAS(IQM(I,J),JCO(1,I,J),NUC(IBLOCK(J),J),
     &                 NRC(IBLOCK(J),J),SEG(IBLOCK(J),J),
     &                 ALPHA(1,IBLOCK(J),J),
     &                 CPRIM(1,1,IBLOCK(J),J),CPRIMU(1,1,IBLOCK(J),J),
     &                 NBLOCK,KAOVEC,KPRIM)
               ELSE
                  IQM(I,J) = 0
               END IF
C
C         Get basis set from basis set library
C         ===================================
C
            ELSE IF ((BASIS .OR. ATOMBA 
     &               .OR. (J .EQ. ISETHUCKEL )
     &               .OR. (ATOMDF .AND. J .EQ. 3 )) 
     &               .AND. .NOT. DIRAC) THEN
              IF (J.EQ.ISETHUCKEL) THEN
C
C                However, make pseudo-label 'pointcharge' for all cases
C                where we should not add Huckel basis (to avoid calling BASLIB below)
C                (Feb. 2006: exept for QM3 / MM which is also skipped below)
                 IF (Q .EQ. D0) BASNAM = 'pointcharge'
C                ... avoid floating orbitals
                 IF (IQM(1,1) .EQ. 0) BASNAM = 'pointcharge'
C                ... avoid point charges (maybe set with INTGRL so BASNAM undefined)
                 IF (BASNAM(1:2) .EQ. 'MM') BASNAM = 'pointcharge'
C                ... avoid MM atoms
                 IF (BASNAM(1:2) .EQ. 'MM') BASNAM = 'pointcharge'
C                ... avoid MM atoms

                 IF (BASNAM(1:11) .NE. 'pointcharge') BASNAM = 'HUCKEL'
C                ... tell BASLIB that this is for Huckel
              ELSE IF (BASIS) THEN
C                Basis-set name for identifier MBSI (WK/UniKA/04-11-2002).
                 BASNAM = MULNAM(MBSI)
              ELSE IF (ATOMBA) THEN

C...             added by Bin Gao, for effective core potentials, Dec. 6, 2010
C...             sets the ECP
                 IPOS = index(BKPLIN, 'ECP=')
                 if (IPOS.ne.0 .and. J.eq.1) then
C     Read effective core potential information from file for atom. 
C     The charge of atom type will be set to the effective charge below
C     after J=2: reading Huckel basis (if DOHUCKEL)
                   if(ADDSYM) call QUIT('ECP not implemented '//
     &               'with auto-detected symmetry.')
C...               tries to find the name of ECP
                   IPOS = IPOS+index(MLINE(NMLINE-NONT(I))(IPOS:), '=')
C...               removes extra spaces
                   IPOS2 = 0
                   do IPOS3 = IPOS, len(MLINE(NMLINE-NONT(I)))
                     if (MLINE(NMLINE-NONT(I))(IPOS3:IPOS3).ne.' ') then
                       IPOS2 = IPOS3
                       exit
                     end if
                   end do
                   if (IPOS2.eq.0) then
                     write(LUPRI,'(/A/A//A)')
     &                 ' Incorrect input for ECP name in line:',
     &                 MLINE(NMLINE-NONT(I)),
     &                 ' Format is "ECP=?"'
                     call QUIT('Incorrect input for ECP name')
                   end if
C...               gets the name of ECP
                   IPOS3 = index(trim(MLINE(NMLINE-NONT(I))(IPOS2:)),
     &                           ' ')
                   if (IPOS3.eq.0)
     &               IPOS3 = len_trim(MLINE(NMLINE-NONT(I))(IPOS2:))
                   IPOS3 = IPOS3+IPOS2-1
                   ECPNAME = MLINE(NMLINE-NONT(I))(IPOS2:IPOS3)
C...               reads in ECP data
                   call ECPINP(Q, QEFF, I, IPREAD, ECPNAME)
                   INDECP(NTYECP) = NCTOT+1 - NONT(I)
                   NECP(NTYECP)   = NONT(I)
                   if (NECP(NTYECP).gt.MXNONT) then
                     write(LUPRI,*) ' NUMBER OF ECP ATOMS OF TYPE ',
     &                 NTYECP, ' EXCEEDS MXNONT=', MXNONT
                     call QUIT('Too many ECP atoms of this type')
                   end if
                 end if
C...
                 IPOS = INDEX(BKPLIN,'BAS')
                 IF (IPOS .NE. 0) THEN
                    IPOS2 = INDEX(MLINE(NMLINE-NONT(I))(IPOS:),'=')
                    IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 6)) THEN
                       WRITE (LUPRI,'(/A/A//A)')
     &                   ' Incorrect input for basis set name in line:',
     &                   MLINE(NMLINE-NONT(I)),
     &                   ' Format is "Basis=?"'
                       CALL QUIT('Incorrect input for basis set name')
                    ELSE
                       IF (ATOMDF) THEN
                         IPOS3 = INDEX(BKPLIN,'AUX')
                         IF (IPOS3 .NE. 0) THEN
                           IPOS4 = INDEX(MLINE(NMLINE-NONT(I))(IPOS3:), 
     &                                                             '=')
                           IF (IPOS4 .EQ. 0 .OR. (IPOS4 .GT. 11)) THEN
                              WRITE (LUPRI,'(A/A//A)')
     &                          'Incorrect input for auxiliary basis'//
     &                          ' set name in input line:',
     &                          MLINE(NMLINE-NONT(I)),
     &                          'Format is "Aux=?"'
                              CALL QUIT('Incorrect input for '//
     &                             'auxiliary basis set name')
                           ELSE
                             AUXNAM =
     &                       MLINE(NMLINE-NONT(I))((IPOS3+IPOS4):)
                             BASNAM =
     &                       MLINE(NMLINE-NONT(I))((IPOS+IPOS2):IPOS3-1)
                           ENDIF
                         ELSE
                           WRITE (LUPRI,'(A,I4,A/A)') 'No auxiliary ' //
     &                          'basis set given for atom kind # ', I,
     &                          ' - input line:',MLINE(NMLINE-NONT(I))
                           CALL QUIT('Incorrect input for '//
     &                          'auxiliary basis set name')
                         ENDIF
                       ELSE
                         BASNAM = MLINE(NMLINE-NONT(I))((IPOS+IPOS2):)
                       ENDIF
                    END IF
                 ELSE
                    WRITE (LUPRI,'(A,I4,A/A)') 'No basis set given '//
     &                   'for atom kind # ', I,
     &                    ' - input line:',MLINE(NMLINE-NONT(I))
                    CALL QUIT('Incorrect input for '//
     &                   'basis set name')
                 END IF
                 DO ILEN = 1, LEN(BASNAM)
                    IF (BASNAM(ILEN:ILEN) .NE. ' ') GO TO 33
                 END DO
 33              CONTINUE
                 BASNAM = BASNAM(ILEN:)
                 IF (ATOMDF) THEN
                   DO IAUX = 1, LEN(AUXNAM)
                      IF (AUXNAM(IAUX:IAUX) .NE. ' ') GO TO 34
                   END DO
 34                CONTINUE
                   AUXNAM = AUXNAM(IAUX:)
                 ENDIF
                 LBASNAM = LNBLNK(BASNAM)
                 IF (IPREAD .GT. 0) THEN
                    IF (J.EQ.ISETHUCKEL) THEN
                        IF (BASNAM(1:2) .NE. 'MM') WRITE(LUPRI,'(2X,A)')
     &                    'Huckel basis read for this type.'
                    ELSE IF (BASNAM(1:6) .EQ. 'INTGRL') THEN
                       WRITE (LUPRI,'(2X,A)')
     &                    'The basis set will be read from input file.'
                    ELSE IF (BASNAM(1:11) .EQ. 'pointcharge') THEN
                       WRITE (LUPRI,'(2X,A)') 'This is a point charge'//
     &                    ' without basis functions.'
                    ELSE IF (BASNAM(1:2) .EQ. 'MM') THEN
                       WRITE (LUPRI,'(2X,A)')
     &                    'This is an MM atom without basis functions.'
C                    ELSE IF (INDEX(BASNAM,'ecp') .NE. 0) THEN
C                       WRITE (LUPRI,'(2X,3A)')
C     &                 'Basis set : '//
C     &                 'the effective core potential basis set "',
C     &                 BASNAM(1:LBASNAM),'"'
                    ELSE IF (ATOMDF) THEN
                       LAUXNAM = LNBLNK(AUXNAM)
                       WRITE (LUPRI,'(2X,3A)')
     &                 'The auxiliary density fit basis set is "',
     &                 AUXNAM(1:LAUXNAM),'" from the basis set library.'
                    ELSE
                       WRITE (LUPRI,'(2X,3A)') 'The basis set is "',
     &                 BASNAM(1:LBASNAM),'" from the basis set library.'
                    END IF
                 END IF
              END IF
              BASNAMSAVE = BASNAM
              IF (BASNAM(1:6) .EQ. 'INTGRL') THEN
C                This makes it possible to specify basis sets
C                not in basis library for some atoms.
C                Use DIRAC style format for specifying IQM and JCO.
                 CALL RDLINE(LUINFO)
                 READ (MLINE(NMLINE),'(BN,A5,12I5)')
     &                BSET,IQM(I,J),(JCO(K,I,J),K=1,IQM(I,J))
                 IF (BSET .EQ. '     ') BSET = 'Input'
                 CALL GTOINP(LUINFO,IQM(I,J),JCO(1,I,J),
     &                NUC(IBLOCK(J),J),NRC(IBLOCK(J),J),
     &                SEG(IBLOCK(J),J),
     &                ALPHA(1,IBLOCK(J),J),CPRIM(1,1,IBLOCK(J),J),
     &                CPRIMU(1,1,IBLOCK(J),J),ISGEN(IBLOCK(J)),
     &                NBLOCK,KAOVEC,KPRIM)
              ELSE IF (BASNAM(1:11) .EQ. 'pointcharge' .OR.
     &                 BASNAM(1: 2) .EQ. 'MM') THEN
                 IQM(I,J) = 0 ! tell later sections that no basis functions
              ELSE
C...  removed by Bin Gao, Dec. 29, 2010
C...                 IF (INDEX(BASNAM,'ecp') .NE. 0 .AND. J.EQ.1) THEN
C...C     Read effective core potential information from file for atom. 
C...C     The charge of atom type will be set to the effective charge below
C...C     after J=2: reading Huckel basis (if DOHUCKEL)
C...                    IF(ADDSYM) CALL QUIT('ECP not implemented '//
C...     &                   'with auto-detected symmetry.')
C...                    CALL ECPINP(Q,QEFF,I,IPREAD)
C...                    INDECP(NTYECP) = NCTOT+1 - NONT(I)
C...                    NECP(NTYECP)   = NONT(I)
C...                    IF (NECP(NTYECP).GT.MXNONT) THEN
C...                       WRITE(LUPRI,*) ' NUMBER OF ECP ATOMS OF TYPE ', 
C...     &                      NTYECP, ' EXCEEDS MXNONT=',MXNONT
C...                       CALL QUIT('Too many ECP atoms of this type')
C...                    END IF
C...                 END IF
C
                 IF (ATOMDF .AND. (J .EQ. 3)) THEN
                    BASNAMSAVE = BASNAM
                    BASNAM     = AUXNAM
                 END IF
C
                 CALL BASLIB(IQM(I,J),JCO(1,I,J),NUC(IBLOCK(J),J),
     &                  NRC(IBLOCK(J),J),SEG(IBLOCK(J),J),
     &                  ALPHA(1,IBLOCK(J),J),
     &                  CPRIM(1,1,IBLOCK(J),J), CPRIMU(1,1,IBLOCK(J),J),
     &                  NBLOCK,ISGEN(IBLOCK(J)),KAOVEC,KPRIM, Q, QEFF,
     &                  DSM, UNCONT, BASNAM, BASREF(1,I,J),IPREAD)
C
                 IF (ATOMDF .AND. (J .EQ. 3)) THEN
                    BASNAM     = BASNAMSAVE
                 END IF
              END IF
              BASNAM = BASNAMSAVE
C
            ELSE
C
C         Read basis set from .mol file (INTGRL in first line of .mol file)
C         ==============================
C
               IF (DIRAC) THEN
                  CALL RDLINE(LUINFO)
                  READ (MLINE(NMLINE),'(BN,A5,12I5)')
     &                 BSET,IQM(I,J),(JCO(K,I,J),K=1,IQM(I,J))
               ENDIF
               CALL GTOINP(LUINFO,IQM(I,J),JCO(1,I,J),
     &              NUC(IBLOCK(J),J),NRC(IBLOCK(J),J),
     &              SEG(IBLOCK(J),J),
     &              ALPHA(1,IBLOCK(J),J),CPRIM(1,1,IBLOCK(J),J),
     &              CPRIMU(1,1,IBLOCK(J),J),ISGEN(IBLOCK(J)),
     &              NBLOCK,KAOVEC,KPRIM)
            ENDIF
            NBLCK(I,J) = NBLOCK
            IBLOCK(J)  = IBLOCK(J) + NBLCK(I,J)
            NOORBTS    = NOORBTS .AND.IQM(I,J).EQ.0
            IF (IBLOCK(J) .GT. MXBLCK) THEN
               WRITE (LUPRI,*) ' Too many blocks of basis functions '//
     &              'in .mol input:',IBLOCK(J),MXBLCK
               WRITE (LUPRI,*) ' Increase work memory in Dalton '//
     &              '("dalton -M mb" option) and be welcome back'
               CALL QUIT('Too many basis function blocks in input, '//
     &              'see output file')
            END IF
            IF (IQM(I,J).GT.0) THEN
               IF(IPREAD .GT. 3) THEN
                  WRITE(LUPRI,'(2X,A5,A)') BSET,' set:'
                  WRITE(LUPRI,'(6X,A,I5,A,12I5)')
     &                 'Max.ang.quantum no.:',(IQM(I,J)-1),
     &                 '  Blocks:',(JCO(K,I,J), K=1,IQM(I,J))
               ENDIF
C
C     Check angular momentum quantum number
C
               NHTYP     = MAX(NHTYP,IQM(I,J))
               IF(NHTYP.GT.MXSPD) GOTO 5000
            ENDIF
C
          ELSE
            WRITE (LUPRI,'(/A,I3)') 'INFO: skipping basis set no.',J
          ENDIF !The main check at the begining of the 300-loop
 300     CONTINUE  ! DO 300 J = 1,NSETS
C
         DO 400 K = 1,NONT(I)
            NOORBT(K+INUC) = NOORBTS ! no orbitals ?
            IF (MBSI .GT. 1) THEN
               CHARGE(K+INUC) = 0.0D0
               ! .. necessary to get nuclear repulsion, nuclear dipole
               ! moment etc. right - this nucleus is already added
               ! for MBSI.eq.1
               NQ = NINT(Q)
               IZATOM(K+INUC) = -NQ
               ! .. code -Z to tell MBSI.gt.1 (value of 0 is floating orbitals)
               GNUEXP(K+INUC) = D0
            ELSE
               ! .. reset charge to effective charge (different from Q for ECP)
               CHARGE(K+INUC) = QEFF
               IF (IQM(I,1) .EQ. 0) THEN
               ! ... this is a point charge (no basis fu.) /hjaaj Nov 2003
                  IZATOM(K+INUC) = -1234567890
                  ! .. special code to tell this is a point charge (no basis fu.) /hjaaj Nov 2003
                  GNUEXP(K+INUC) = D0
                  ! .. make sure the point charge IS a point !! /hjaaj Feb 2006
               END IF
            END IF
 400     CONTINUE
C
C     Currently the program cannot handle symmetry if there is a chance 
C     that the atoms in the input may be different due to different basis
C     sets
C     
         IF (ADDSYM .AND. .NOT. BASIS .AND. IZATOM(INUC+1).ge.0) THEN
C
C           Check if same nucleus was in prior "type" because
C           then it may have a different basis set
C           and FNDSYM routine can currently not detect that
C           (e.g. if some Carbon have STO-3G and others have
C           6-31G**) /hjaaj Mar 2004
C           We skip point charges/MM atoms (IZATOM = -1234567890) and multiple
C           basis set (IZATOM = -Z) for this test, because they
C           cannot destroy symmetry check (either no basis fu.
C           thus same basis set on all point charges (!)
C           or - for multiple basis set - already checked for
C           first set where IZATOM .ge. 0).
C
            J = 0
            DO K = 1, INUC
               IF (IZATOM(INUC+1) .EQ. IZATOM(K)) J = J + 1
            END DO
            IF (J .GT. 0) THEN
              WRITE (LUPRI,'(/A/A,I4/A/A)')
     &        ' The symmetry of the molecule cannot'//
     &        ' be determined when INTGRL or ATOMBASIS is used,',
     &        ' and different atom types are used for same atom --'//
     &        ' problem discovered for Z =',IZATOM(INUC+1),
     &        ' Please specify the molecular symmetry in the input'//
     &        ' or do not use different atom types for the same atom',
     &        ' and be welcome back !'
              CALL QUIT('Symmetry cannot be determined -- see output')
            END IF
         END IF
  100 CONTINUE
C
      IF (.NOT. ONLYOV) THEN
        IF (QM3) THEN
          CALL QM3MAS(NONT,NONTYP)
C
C-------------------------------------------------------------
C Performing some QM3 test after having read the nuclei
C positions and generated the center of masses of the systems.
C-------------------------------------------------------------
C
C 1: Testing if the number of QM classical charges read from input
C    matches the number of QM nuclei:
C
          IF ( ICHRGS(0) .NE. NSISY(0) ) THEN
           CALL QUIT('Error in no. of classical charges for '//
     &                'the QM system.')
          END IF
C
C 1: Testing if number of QM nuclei exceeds MXCENT_QM in mxcent.h
C
          JK = 0
          DO 655 IK = 1, NCTOT
            IF (ISUBSY(IK) .EQ. 0 .AND.
     &          ISUBSI(IK) .LE. NSISY(0)) THEN
            JK = JK + 1
            END IF
 655      CONTINUE
          IF (JK . GT. MXCENT_QM) THEN
            WRITE(LUPRI,'(//A,I5,A,/A,I5,A,/A/)')
     &      ' FATAL ERROR! The number of QM atoms,',JK,',',
     &      'exceeds the maximum number of QM atoms allowed,',
     &      MXCENT_QM,' !!',
     &      ' SOLUTION: Increase MXCENT_QM in include/mxcent.h'//
     &      ' and re-compile.'
            CALL QUIT('Maximum number of QM atoms exceeded')
          END IF
C
C 2: Testing if the number of atomic van der Waals parameters
C    read for each system type matches the number of nuclei
C    in the system:
C
          IF (LONEPAR) THEN
            DO 777 LVDW = 0, ISYTP
              IF (NSISY(LVDW) .NE. ISIGEPS(LVDW)) THEN
                WRITE(LUPRI,'(//1X,A,I3,A,I3)')
     *          'Number of atoms in system No.',
     *          LVDW,' is: ',NSISY(LVDW)
                WRITE(LUPRI,'(1X,A,I3,A,I3//)')
     *          'Number of sigma and epsilon parameters for sys.',
     *          LVDW,' is: ',ISIGEPS(LVDW)
                CALL QUIT('Number of sigma parameters '//
     *          'and atoms differ')
              END IF
 777        CONTINUE
          END IF
C
C 3: Testing if the number of diatomic sigma and epsilon parameters
C    read matches the number generated from the knowledge of sytem-
C    types and number of nuclei for each systemtype:
C
          IF (LTWOPAR) THEN
            INT_TEMP = 0
            DO 888 LVDW = 0, ISYTP
              DO 889 kVDW = LVDW, ISYTP
                IF (LVDW .EQ. KVDW) THEN
                  IF (LVDW .NE. 0) THEN
                    INT_TEMP = INT_TEMP + IFAC(NSISY(LVDW))
                  END IF
                ELSE
                    INT_TEMP = INT_TEMP + NSISY(LVDW) * NSISY(KVDW)
                END IF
 889          CONTINUE
 888        CONTINUE
            IF (INT_TEMP .NE. NSIGEPS) THEN
              WRITE(LUPRI,'(1X,A,I3,A)')
     *                    'Needs ',INT_TEMP, 'sigma/epsilon parameters!'
              CALL QUIT('Number of diatomic sigma/epsilon parameters '//
     *                  'incorrect!')
            END IF
          END IF
C
        END IF
      END IF
      RETURN
C
C     Error messages:
C
 5000 CONTINUE
         WRITE (LUPRI,'(6X,A,I3,3A/9X,2A)')
     &      '*  Input specifies highest orbital of atomic type ',
     &      I,' as "',SPDCAR(NHTYP - 1),'".',
     &      ' Highest allowed orbital in this version: ',
     &      SPDCAR(MXSPD - 1)
         IF (NHTYP.GT.MXQN) WRITE (LUPRI,'(9X,2(A,I3),A)')
     &      ' Increase MXQN from',MXQN,' to',NHTYP,' and recompile.'
         CALL QUIT('Too high angular quantum no. specified in input.')
 299  CONTINUE
      WRITE (LUPRI,'(A,1X,I4)') 'Cannot determine the charge of atom '//
     &     'type', I
      WRITE (LUPRI,'(A)') "Correct input format is: Charge=???"
      CALL QUIT('Cannot determine the charge of one atom')
      END ! SUBROUTINE BASINP
C  /* Deck cntinp */
      SUBROUTINE CNTINP(LUINFO,NONT,Q,Q_EXPONENT,NSYMOP,ANG,MBSI,
     &                  RADIUS_PCM, ALPHA_PCM)
C*****************************************************************************
C
C     Read and process data about symmetry independent centers
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "dummy.h"
#include "codata.h"
#include "cbirea.h"
C
      PARAMETER(D0 = 0.0D00, CORMAX = 1.D5)
      LOGICAL ANG
      CHARACTER*4 FRMT

C     include character ASYMB(0:NSYMB)*2 with atomic symbols
#include "asymb.h"

#include "gnrinf.h"

#include "molinp.h"
#include "ccom.h"
#include "nuclei.h"
#include "symmet.h"
#include "r12int.h"
#include "pcmdef.h"
#include "pcmlog.h"
#include "pcm.h"
#include "qm3.h"
C
      NQ = NINT(Q)
C
      IF (GAUNUC.AND.(NQ.GT.0)) THEN
         IF(Q_EXPONENT.LE.D0) THEN
            CALL NUCSIZ(NQ,GEXP)
         ELSE
            GEXP = Q_EXPONENT
         ENDIF
         IF (IPREAD.GE.0) WRITE(LUPRI,'(2X,A2,A,1P,D15.5/)') ASYMB(NQ),
     &   '  : Exponent for nuclear gaussian charge distribution : ',GEXP
      ELSE
         GEXP = D0
      END IF
      NUCSTR = NUCIND + 1
      IF (IPREAD .GT. 3) 
     &     WRITE(LUPRI,'(3X,A)') 'Symmetry independent centers:'
      DO 100 N = 1,NONT
        NCTOT = NCTOT + 1
        IF (NCTOT .GT. MXCENT) GOTO 5000
C
C*****************************************************************************
C       Card 7:
C          NAMN - name of symmetry independent center
C          CORR - coordinates of symmetry independent center
C*****************************************************************************
C
        CALL RDLINE(LUINFO)
        IPOS = INDEX(MLINE(NMLINE)(2:80),' ') + 1
C       sep05-hjaaj: allow both e.g. 'H ' and ' H' as names
        IF (IPOS .LE. 2 .AND. MLINE(NMLINE)(1:1) .EQ. ' ') THEN
          WRITE (LUPRI,'(/A,I5/A/A)') 'Error in input line no. ',NMLINE,
     &    MLINE(NMLINE),'Atom name must start in first or second column'
          CALL QUIT('Error in placement of atom name. See output')
        ELSE IF (IPOS .EQ. 2) THEN
          FRMT = '(A1)'
        ELSE IF (IPOS .EQ. 3) THEN
          FRMT = '(A2)'
        ELSE IF (IPOS .EQ. 4) THEN
          FRMT = '(A3)'
        ELSE
          FRMT = '(A4)'
        END IF
        READ (MLINE(NMLINE),FRMT) NAMN(NCTOT)
        IF (NAMN(NCTOT)(1:4) .EQ. '    ' .AND.
     &      NQ .GE. 0 .AND. NQ .LE.NSYMB) THEN
           WRITE(NAMN(NCTOT),'(A2,I2)') ASYMB(NQ), N
        END IF
C
C       951115-hjaaj: NAMN() must never be read with free format
C        (otherwise it will on e.g. AIX be filled with nulls:
C         e.g. ' C  ' would become 'C\NULL\NULL\NULL'; and on
C         e.g. IRIX free format read of char. var. is not defined).
C
#if !defined (SYS_CRAY) && !defined (VAR_NOFREE) && !defined (SYS_T3D) && !defined (SYS_T90)
C
        IF (QM3) THEN
           READ (MLINE(NMLINE)(IPOS:80),*,ERR=101,END=101)
     &         CORD(1,NCTOT),CORD(2,NCTOT),CORD(3,NCTOT),
     &         ISUBSY(NCTOT),ISUBSI(NCTOT)
        ELSE
           READ (MLINE(NMLINE)(IPOS:80),*,ERR=101,END=101)
     &         CORD(1,NCTOT),CORD(2,NCTOT),CORD(3,NCTOT)
        END IF
        GO TO 104
C
 101    READ (MLINE(NMLINE)(IPOS:80),
     &        '(BN,3F20.0,2I6)',ERR=102)
     &        CORD(1,NCTOT),CORD(2,NCTOT),CORD(3,NCTOT),
     &        ISUBSY(NCTOT),ISUBSI(NCTOT)
        GO TO 104
C
 102    READ (MLINE(NMLINE)(IPOS:80),
     &       '(BN,3F10.0,2I6,I5)',ERR=103)
     &       CORD(1,NCTOT),CORD(2,NCTOT),CORD(3,NCTOT),
     &       ISUBSY(NCTOT),ISUBSI(NCTOT)
        GO TO 104
C
  103   CONTINUE
        WRITE(LUPRI,'(/A,I5/A,I5,A)')
     &     ' ERROR: Unable to read Cartesian coordinates of atom no.',
     &     NCTOT,' from line',NMLINE,' in the MOLECULE input file:'
        WRITE(LUPRI,'(A)') MLINE(NMLINE)
        CALL QUIT('ERROR reading atomic coordinates in MOLECULE input')
C
  104   CONTINUE
C
#else
        ISTART = IPOS
        DO J = 1, 3
          CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,CORD(J,NCTOT),'REA',
    &                IERR)
        END DO
        IF (QM3) THEN
          CALL FREFRM(MLINE(NMLINE),ISTART,ISUBSY(NCTOT),DUMMY,
     &                'INT',IERR)
          CALL FREFRM(MLINE(NMLINE),ISTART,ISUBSI(NCTOT),DUMMY,
     &                'INT',IERR)
        END IF
              
#endif
        IF (PCM .AND. (RADIUS_PCM * ALPHA_PCM .GT. 0.1)) THEN
           IDXSPH(NCTOT) = 1
           RSPH(NCTOT)   = RADIUS_PCM * ALPHA_PCM
        ELSE
           IDXSPH(NCTOT) = 0
           RSPH(NCTOT)   = 0.0D0
        END IF
        IF (ISUBSY(NCTOT) .EQ. 0) NUCIND = NUCIND + 1
        IPOS = INDEX(MLINE(NMLINE),'Isotope=')
        IF (IPOS .NE. 0) THEN
          IPOS = IPOS + 8
          READ (MLINE(NMLINE)(IPOS:),'(I3)') MASSNM
          ISOTOP(NCTOT) = ISOMAS(NQ,MASSNM)
        ELSE
          ISOTOP(NCTOT) = 1
        END IF
C
        IF (IPREAD .GT. 3) THEN
           WRITE(LUPRI,'(6X,A4,3F20.15)') NAMN(NCTOT),
     &          (CORD(J,NCTOT), J = 1,3)
        ENDIF
C
C     To avoid problems in later stages, we now rewrite the coordinates
C     of the nuclei in traditional Dalton format, with 4 characters for the
C     name of the atom
C
C     Arnfinn, nov. -08: ISUBSY and ISUBSI included in MLINE if QM3
C
        IF (QM3) THEN
           WRITE (MLINE(NMLINE),'(A4,3F20.10,2X,2I4,A3)') NAMN(NCTOT),
     &        (CORD(J,NCTOT),J = 1, 3),ISUBSY(NCTOT),ISUBSI(NCTOT),'   '
C
        ELSEIF (IPOS .EQ. 0) THEN
#if (defined(GFORTRAN) && GFORTRAN <= 402)
c     a bug in the compiler (seen with gfortran-4.0.2 20051125)
c     allegedly fixed in 4.1.x series.
          WRITE (MLINE(NMLINE),'(A4,3F20.10,2X,A11)') NAMN(NCTOT),
     &          (CORD(J,NCTOT),J = 1, 3), '             '
#else
          WRITE (MLINE(NMLINE),'(A4,3F20.10,2X,A14)') NAMN(NCTOT),
     &          (CORD(J,NCTOT),J = 1, 3), '                '
#endif
        ELSE
          IF (MASSNM .LT. 10) THEN
            WRITE (MLINE(NMLINE),'(A4,3F20.10,A10,I1,A5)') 
     &             NAMN(NCTOT),(CORD(J,NCTOT),J = 1, 3), 
     &             '  Isotope=', MASSNM, '     '
          ELSE IF (MASSNM .LT. 100) THEN
            WRITE (MLINE(NMLINE),'(A4,3F20.10,A10,I2,A4)') 
     &             NAMN(NCTOT),(CORD(J,NCTOT),J = 1, 3), 
     &             '  Isotope=', MASSNM, '    '
          ELSE
            WRITE (MLINE(NMLINE),'(A4,3F20.10,A10,I3,A3)') 
     &             NAMN(NCTOT),(CORD(J,NCTOT),J = 1, 3), 
     &             '  Isotope=', MASSNM, '   '
          END IF
        END IF
C
C
        NCLINE(NCTOT) = NMLINE
        NAMEX(3*NCTOT)     = NAMN(NCTOT)//' z'
        NAMEX(3*NCTOT - 1) = NAMN(NCTOT)//' y'
        NAMEX(3*NCTOT - 2) = NAMN(NCTOT)//' x'
        DO 200 J = 1,3
           IF(ANG) CORD(J,NCTOT) = CORD(J,NCTOT)/XTANG
           IF(ABS(CORD(J,NCTOT)).GT.CORMAX) GOTO 5010
  200   CONTINUE
C
C*****************************************************************************
C       MULBSI  - basis-set identifier (WK/UniKA/04-11-2002).
C       CHARGE  - charge of center
C       NOORBT  - TRUE: no orbitals on this center
C       GNUEXP  - exponent of Gaussian nuclear charge distribution
C*****************************************************************************
C
        MULBSI(NCTOT) = MAX(MBSI,1)
        IF (MBSI .GT. 1) THEN
          CHARGE(NCTOT) = D0
          IZATOM(NCTOT) = -NQ
C         ... code -1 to tell MBSI (value of 0 is floating orbitals)
          GNUEXP(NCTOT) = D0
        ELSE
          CHARGE(NCTOT) = Q
          IZATOM(NCTOT) = NQ
C         ... if point charge this is reset outside to -1234567890
C             to tell this is a point charge /hjaaj Nov 2003
          GNUEXP(NCTOT) = GEXP
        END IF
  100 CONTINUE   !  DO 100 N = 1,NONT
      RETURN
C
C       Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A/A,I5)')
     &    ' Too many atomic centers: MXCENT exceed in CNTINP,',
     &    ' Current limit:',MXCENT
        CALL QUIT('*** ERROR *** MXCENT exceeded in CNTINP')
 5010 CONTINUE
        WRITE (LUPRI,'(A,E12.6,A/A/A,E12.6)')
     &    ' Atomic coordinate ',CORD(J,NCTOT),
     &    ' too large in CNTINP.',
     &    ' Note: Program is unstable for large coordinates.',
     &    ' Maximum coordinate value:',CORMAX
        CALL QUIT('*** ERROR: Atomic coordinate too large in CNTINP')
      END
C  /* Deck baspro */
      SUBROUTINE BASPRO(WORK,LWORK,NONTYP,NSYMOP,NONT,IQM,NBLCK,JCO,NUC,
     &                  NRC,SEG,ALPHA,CPRIM,CPRIMU,KATOM,KANG,KSETS,
     &                  KBLOCK,KPRIM,DOOWN)
C
C*****************************************************************************
C
C       Process orbital and geometry input data
C
C*****************************************************************************
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <cbirea.h>
      DIMENSION WORK(LWORK)
      LOGICAL SEG,DOOWN
      DIMENSION NONT(KATOM),IQM(KATOM,3),JCO(KANG,KATOM,3),
     &          NBLCK(KATOM,3),NUC(KBLOCK,3),NRC(KBLOCK,3),
     &          SEG(KBLOCK,3),
     &          ALPHA(KPRIM,KBLOCK,3),CPRIM(KPRIM,KPRIM,KBLOCK,3),
     &          CPRIMU(KPRIM,KPRIM,KBLOCK,3)

#include "gnrinf.h"

#include "ccom.h"
#include "nuclei.h"
#include "shells.h"
#include "symmet.h"
#include "primit.h"
#include "huckel.h"
#include "aosotr.h"

#include "veclen.h"

#include "ibtfun.h"
C
C*****************************************************************************
C       MULK    - bitstring of basic operations that stabilise a center
C*****************************************************************************
C
      II = 1
      DO 100 I = 1, NONTYP
         DO 110 J = 1, NONT(I)
            MULK = 0
            LL = 1
            DO 140 L = 1, NSYMOP
               DO 150 M = 1,3
                  IF(IBTAND(LL,ISYMAX(M,1)) .NE. 0) THEN
C                 Basic operation LL moves coordinate M
                     IF(ABS(CORD(M,II)).GE.1.0D-06) GOTO 140
                  END IF
 150           CONTINUE
               MULK = MULK + LL
 140        LL = 2*LL
            ISTBNU(II) = MULK
            II = II + 1
 110     CONTINUE
 100  CONTINUE

C*****************************************************************************
C     Set various basis parameters
C*****************************************************************************
      CALL BASPAR(NSYMOP)

C*****************************************************************************
C     Process nuclear data
C*****************************************************************************
      CALL NUCPRO(WORK,LWORK)

C*****************************************************************************
C     Process orbital data for large components
C*****************************************************************************
      ISHELL = 0
      IPRIM  = 0
      IPRIMD = 0
      IORB   = 0
      IORBD  = 0

      LCOMP  = 1
      CALL ORBPRO(NONTYP,NONT,IQM(1,1),NBLCK(1,1),JCO(1,1,1),NUC(1,1),
     &            NRC(1,1),SEG(1,1),
     &            ALPHA(1,1,1),CPRIM(1,1,1,1),CPRIMU(1,1,1,1),
     &            KATOM,KANG,KBLOCK,KPRIM,ISHELL,
     &            IPRIM,IPRIMD,IORB,IORBD,LCOMP)
      NLRGSH = ISHELL
      NLARGE = IORBD
      NORBL  = IORB
      NPLRG  = IPRIMD
      NPLSH  = IPRIM
C
C     Set up a contraction matrix in case of Douglas-Kroll integrals
C
      IF (DKTRAN) THEN
         KLPRIM = 1
         KDKWRK = KLPRIM + KPRIM
         LDKWRK = LWORK - KDKWRK
         IF (LDKWRK .LT. NPLRG*NLARGE)
     &      CALL STOPIT('READ_MOL ','DKPRO',NPLRG*NLARGE,LDKWRK)
         CALL DKPRO(NONTYP,NONT,IQM(1,1),NBLCK(1,1),JCO(1,1,1),
     &              NUC(1,1),NRC(1,1),SEG(1,1),ALPHA(1,1,1),
     &              CPRIM(1,1,1,1),WORK(KLPRIM),WORK(KDKWRK),
     &              KATOM,KANG,KBLOCK,KPRIM,NPLRG,NLARGE)
      END IF
C*****************************************************************************
C       Process orbital data for small components or Huckel
C*****************************************************************************

      LCOMP = 2
      CALL ORBPRO(NONTYP,NONT,IQM(1,2),NBLCK(1,2),JCO(1,1,2),NUC(1,2),
     &            NRC(1,2),SEG(1,2),
     &            ALPHA(1,1,2),CPRIM(1,1,1,2),CPRIMU(1,1,1,2),
     &            KATOM,KANG,KBLOCK,KPRIM,ISHELL,
     &            IPRIM,IPRIMD,IORB,IORBD,LCOMP)
      IF (DOHUCKEL) THEN
         IF (ISHELL .GT. MXSHEL) GOTO 5040
         IF (IORBD  .GT. MXCORB) GOTO 5050
         IF (IPRIMD .GT. MXPRIM) GOTO 5060
         NSMLSH = ISHELL - NLRGSH ! really "NHUCSH", reusing small comp. variable from Dirac
         NORBS  = IORB
         NSMALL = IORBD  - NLARGE ! really "NHUCBA", reusing small comp. variable from Dirac
         NPSML  = IPRIMD - NPLRG  ! really "NPHUC ", reusing small comp. variable from Dirac
C
         KMAX   = NLRGSH
         NBASIS = NLARGE
         NPBAS  = NPLRG
C
      ELSE
         KMAX   = ISHELL
         IF(KMAX.GT.MXSHEL) GOTO 5000
         NSMLSH = KMAX - NLRGSH
         NORBS  = IORB
         NBASIS = IORBD
         IF (IPREAD .GT. 10) THEN
            WRITE(LUPRI,'(A)') 'PRIEXP:'
            WRITE(LUPRI,'(I5,F15.5)') (KK,PRIEXP(KK),KK=1,IPRIM)
         END IF
         IF (NBASIS .GT. MXCORB) GOTO 5010
         NSMALL = NBASIS - NLARGE
         NPBAS  = IPRIMD
         IF (NPBAS .GT. MXPRIM)  GOTO 5020
         NPSML   = IPRIMD - NPLRG
      END IF
C*****************************************************************************
C       Process orbital data for auxiliary basis (density fitting)
C*****************************************************************************
      IF (ATOMDF) THEN
         LCOMP = 3
         CALL ORBPRO(NONTYP,NONT,IQM(1,3),NBLCK(1,3),JCO(1,1,3),
     &               NUC(1,3),NRC(1,3),SEG(1,3),
     &               ALPHA(1,1,3),CPRIM(1,1,1,3),CPRIMU(1,1,1,3),
     &               KATOM,KANG,KBLOCK,KPRIM,
     &               ISHELL,IPRIM,IPRIMD,IORB,IORBD,LCOMP)
         NAUXSH    = ISHELL - NSMLSH - NLRGSH - 1
         NAUX      = IORBD  - NLARGE - NSMALL
         NORBAUX   = IORB
         NPAUX     = IPRIMD - NPLRG  - NPSML  - 1
         KMAXAUX   = NAUXSH
         NBASISAUX = NAUX
         NPBASAUX  = NPAUX
         KMAXTOT   = ISHELL
      ELSE
         NAUXSH    = 0
         NAUX      = 0 
         NORBAUX   = 0
         NPAUX     = 0 
         KMAXAUX   = 0 
         NBASISAUX = 0 
         NPBASAUX  = 0 
         KMAXTOT   = KMAX  
      END IF

      NPSHEL = IPRIM ! total number of primitive shells

CNECgh980808
      NODD = 1 - MOD(NBASIS,2)
C
C     ***************************************
C     ***** Orbital Symmetry Processing *****
C     ***************************************
C
      KKVAL = 1
      KMVAL = KKVAL + MXAQN
      KNVAL = KMVAL + MXAQN
      KIREP = KNVAL + MXAQN
      KLAST = KIREP + MXCORB
      IF (KLAST .GT. LWORK) GOTO 5030
      CALL SYMPRO(WORK(KKVAL),WORK(KMVAL),WORK(KNVAL),WORK(KIREP),
     &            DOOWN)
C
      RETURN
C
C       Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> BASPRO error, number of shells                    ',KMAX,
     *  '                   current maximum number              ',MXSHEL
        CALL QUIT('BASPRO: Too many shells')
 5010 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> BASPRO error, number of contracted basis functions',
     *  NBASIS,
     *  '                   current maximum number              ',MXCORB
        CALL QUIT('BASPRO: Too many contracted basis functions')
 5020 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> BASPRO error, number of primitive basis functions',NPBAS,
     *  '                   current maximum number             ',MXPRIM
        CALL QUIT('BASPRO: Too many primitives')
 5030 CONTINUE
        CALL STOPIT('BASPRO','SYMPRO',KLAST,LWORK)
 5040 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> BASPRO error, # of shells when adding Huckel orbitals',
     &     ISHELL,
     *  '                   current maximum number                 ',
     &     MXSHEL
        CALL QUIT('BASPRO: Too many shells')
 5050 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> BASPRO error, # of contracted bf. after adding Huckel '//
     &     'orbitals',IORBD,
     *  '                   current maximum number                  '//
     &     '        ',MXCORB
        CALL QUIT('BASPRO: Too many contracted basis functions')
 5060 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> BASPRO error, # of primitive bf after adding Huckel '//
     &     'orbitals',IPRIMD,
     *  '                   current maximum number                 '//
     &     '        ',MXPRIM
        CALL QUIT('BASPRO: Too many primitives')
      END
C  /* Deck basout */
      SUBROUTINE BASOUT(WORK,LWORK,NONTYP,NONT,IQM,NBLCK,JCO, BASREF,
     &                  NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,KATOM,KANG,KSETS,
     &                  KBLOCK,KPRIM,KASYM,NSYMOP,HERMIT)
C*****************************************************************************
C
C     Print orbital and geometry input data
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "pi.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "cbirea.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
      PARAMETER (LUCME = -1)
#else
#include "inftap.h"
#endif
#include "symmet.h"
      DIMENSION WORK(LWORK)
      LOGICAL   SEG,HERMIT,NRMPRI
      DIMENSION NONT(KATOM),IQM(KATOM,KSETS),JCO(KANG,KATOM,KSETS),
     &          NBLCK(KATOM,KSETS),NUC(KBLOCK,KSETS),NRC(KBLOCK,KSETS),
     &          SEG(KBLOCK,KSETS),
     &          ALPHA (KPRIM,KBLOCK,KSETS),
     &          CPRIM (KPRIM,KPRIM,KBLOCK,KSETS),
     &          CPRIMU(KPRIM,KPRIM,KBLOCK,KSETS)
      CHARACTER*80 BASREF(10,KATOM,KSETS)
      CHARACTER*1  KASYM(3,3)
      CHARACTER*(6*MXQN+3) TSTRING(MXBSETS_TOT)
      CHARACTER*37 SET_TYPE(MXBSETS_TOT)
      DIMENSION NPRIM(MXBSETS_TOT), NCONT(MXBSETS_TOT)
      DIMENSION IBLOCK(MXBSETS_TOT)
#include <molde.h>
#include "ccom.h"
#include "nuclei.h"
#include <frame.h>
#include <infpar.h>
#include "shells.h"
#include "primit.h"
#include "aosotr.h"
#include "r12int.h"
#include "qm3.h"
C
      IF (DIRAC) THEN
         SET_TYPE(1) = 'large component basis functions'
         SET_TYPE(2) = 'small component basis functions'
         SET_TYPE(3) = 'primary density fit functions  '
         SET_TYPE(4) = 'secondary density fit functions'
         SET_TYPE(5) = 'extended Huckel basis functions'
C
      ELSE
         SET_TYPE(1) = 'basis functions'
         SET_TYPE(2) = 'extended Huckel basis functions'
         SET_TYPE(3) = 'auxiliary density fit basis functions'
      END IF
C
C     ******************************************
C     ***** Printing of basis information  *****
C     ******************************************
C
      IF (IPREAD .GT. 0) THEN
         CALL HEADER('Atoms and basis sets',1)
         WRITE (LUPRI,'(A,I5 )') '  Number of atom types :',NONTYP
         WRITE (LUPRI,'(A,I5/)') '  Total number of atoms:',NUCDEP
         IF (LUCME.GT.0) THEN
            WRITE (LUCME,'(A,I5 )') '  Number of atom types :',NONTYP
            WRITE (LUCME,'(A,I5/)') '  Total number of atoms:',NUCDEP
         END IF
         IF (BASIS) THEN
          IF (LMULBS) THEN
C           Output for multiple basis sets (WK/UniKA/04-11-2002).
           DO IMULBS = 1, NMULBS
            ILEN = LNBLNK(MULNAM(IMULBS))
            WRITE (LUPRI,'(A,I2,A/)') '  Basis set',IMULBS,
     &      ' refers to the basis "'//MULNAM(IMULBS)(1:ILEN)//
     &      '" from the basis set library.'
           END DO
          ELSE
            ILEN = LNBLNK(MULNAM(1))
            WRITE (LUPRI,'(3A/)') '  Basis set used is "',
     &      MULNAM(1)(1:ILEN),'" from the basis set library.'
          END IF
         END IF
 
         IF (LMULBS) THEN
C         Output for multiple basis sets (WK/UniKA/04-11-2002).
          WRITE (LUPRI,'(A)')
     &        '  label  basis   atoms  charge   prim  cont     basis'
          WRITE (LUPRI,'(2X,78A1)') ('-',I=1,78)
          IF (LUCME.GT.0) THEN
             WRITE (LUCME,'(A)')
     &        '  label  basis   atoms  charge   prim  cont     basis'
             WRITE (LUCME,'(2X,78A1)') ('-',I=1,78)
          END IF
         ELSE
          WRITE (LUPRI,'(A)')
     &        '  label    atoms   charge   prim   cont     basis'
          WRITE (LUPRI,'(2X,70A1)') ('-',I=1,70)
          IF (LUCME.GT.0) THEN
            WRITE (LUCME,'(A)')
     &        '  label    atoms   charge   prim   cont     basis'
            WRITE (LUCME,'(2X,70A1)') ('-',I=1,70)
          END IF
         END IF
      END IF
C
      IBLOCK(1:KSETS) = 1
      ICENT  = 0
      NCHTOT = 0
      NQMBAS = 0
      QM3CHT = 0.0D0
      DO 100 I = 1,NONTYP
C
C        Convert information to a string that can be understood by the innocent user
C
         DO ISET = 1, KSETS
            IF (IQM(I,ISET).GT.0) THEN
               TSTRING(ISET) = ' '
               CALL BASTYP(IQM(I,ISET),JCO(1,I,ISET),
     &                     NRC(IBLOCK(ISET),ISET),
     &                     NUC(IBLOCK(ISET),ISET),
     &                     NPRIM(ISET),NCONT(ISET),TSTRING(ISET))
            ELSE
               NPRIM(ISET)   = 0
               NCONT(ISET)   = 0
               TSTRING(ISET) = 'Point Charge'
!              TSTRING(ISET) = 'No '//SET_TYPE(ISET)//' attached'
            END IF
            IBLOCK(ISET) = IBLOCK(ISET) + NBLCK(I,ISET)
         END DO

         NUCTYP = 0
         NCHTYP = 0
         DO N = 1,NONT(I)
            ICENT  = ICENT + 1
            NUCTYP = NUCTYP + NUCDEG(ICENT)
            CHTTYP = CHARGE(ICENT)

            NCHARG = MAX(0,IZATOM(ICENT))
C           hjaaj: do not add the code for point charges (-1234567890)
C           or the code for multiple basis (-NQ) to the total charge !!
            NCHTOT = NCHTOT + NUCDEG(ICENT)*NCHARG

            QM3CHT = QM3CHT + NUCDEG(ICENT)*CHTTYP
C           ... note that charges for multiple basis sets are not added
C           here because CHARGE of these centers has been reset to zero.
C           However, we do add point charges (incl. MM) /hjaaj Feb 2006
C
            IF (ISUBSY(ICENT) .EQ. 0) NQMBAS = NQMBAS + NCONT(ISET)
         END DO
         NCHTOT = NCHTOT + NCHTYP
C
         IF (IPREAD .GT. 0) THEN
            IF(DIRAC) THEN
               WRITE (LUPRI,'(2X,A2,3X,4I8,6X,2A)') 
     &            NAMN(ICENT)(1:2),NUCTYP,NCHARG,
     &            NPRIM(1),NCONT(1),'L - ',TSTRING(1)
               WRITE (LUPRI,'(23X,2I8,6X,2A)')
     &            NPRIM(2),NCONT(2),'S - ',TSTRING(2)
            ELSE
               IF (LMULBS) THEN
C               Output for multiple basis sets (WK/UniKA/04-11-2002).
                  WRITE (LUPRI,'(2X,A2,3X,I5,I8,F10.4,2I6,6X,A)') 
     &                 NAMN(ICENT)(1:2),MULBSI(ICENT),NUCTYP,
     &                 CHTTYP,NPRIM(1),NCONT(1),TSTRING(1)
                  IF (LUCME.GT.0)
     &            WRITE (LUCME,'(2X,A2,3X,I5,I8,F10.4,2I6,6X,A)')
     &                 NAMN(ICENT)(1:2),MULBSI(ICENT),NUCTYP,
     &                 CHTTYP,NPRIM(1),NCONT(1),TSTRING(1)

               ELSE
                  WRITE (LUPRI,'(2X,A2,3X,I8,F10.4,2I6,6X,A)') 
     &                 NAMN(ICENT)(1:2),NUCTYP,CHTTYP,
     &                 NPRIM(1),NCONT(1),TSTRING(1)
                  IF (LUCME.GT.0)
     &            WRITE (LUCME,'(2X,A2,3X,I8,F10.4,2I6,6X,A)')
     &                 NAMN(ICENT)(1:2),NUCTYP,CHTTYP,
     &                 NPRIM(1),NCONT(1),TSTRING(1)
               END IF
            END IF
         END IF
  100 CONTINUE
      IF (IPREAD .GT. 0) THEN
         IF(DIRAC) THEN
            WRITE (LUPRI,'(2X,70A1)') ('-',I=1,70)
            WRITE(LUPRI,'(23X,2I8,3X,A)')
     &            NPLRG,NLARGE,'L - large components'
            WRITE(LUPRI,'(23X,2I8,3X,A)')
     &            NPSML,NSMALL,'S - small components'
         ENDIF
         TSYSCH = QM3CHT
         IF (LMULBS) THEN
C     Output for multiple basis sets (WK/UniKA/04-11-2002).
            WRITE (LUPRI,'(2X,78A1)') ('-',I=1,78)
            WRITE (LUPRI,'(2X,A,I12,F10.4,2I6)')
     &           'total:',NUCDEP,TSYSCH,NPBAS,NBASIS
            IF (LUCME.GT.0) 
     &           WRITE (LUCME,'(2X,A,I12,F10.4,2I6/)')
     &           'total:',NUCDEP,TSYSCH,NPBAS,NBASIS
            WRITE (LUPRI,'(2X,78A1)') ('-',I=1,78)
         ELSE
            WRITE (LUPRI,'(2X,70A1)') ('-',I=1,70)
            WRITE (LUPRI,'(2X,A,I7,F10.4,2I6)')
     &           'total:',NUCDEP,TSYSCH,NPBAS,NBASIS
            IF (LUCME.GT.0) 
     &           WRITE (LUCME,'(2X,A,I7,F10.4,2I6/)')
     &           'total:',NUCDEP,TSYSCH,NPBAS,NBASIS
            WRITE (LUPRI,'(2X,70A1)') ('-',I=1,70)
         END IF
         IF (NHTYP .GT. 2) THEN
            IF(DOCART) THEN
               WRITE (LUPRI,'(A)') '  Cartesian basis used.',
     &     '  (Note that d, f, ... atomic GTOs are not all normalized.)'
               IF (LUCME.GT.0)
     &         WRITE (LUCME,'(A)') '  Cartesian basis used.',
     &     '  (Note that d, f, ... atomic GTOs are not all normalized.)'
            ELSE
               WRITE (LUPRI,'(A)') '  Spherical harmonic basis used.'
               IF (LUCME.GT.0)
     &         WRITE (LUCME,'(A)') '  Spherical harmonic basis used.'
            ENDIF
         END IF
         WRITE (LUPRI,'(/A,1P,D10.2)')
     &      '  Threshold for integrals:',THRS
      END IF
C
! hjaaj: info not availiable here any more - don't really understand
! that it ever has, because AUXNAM may change for each atom type ???
!     IF (ATOMDF) THEN
!        ILEN = LNBLNK(AUXNAM)
!        WRITE (LUPRI,'(/3A/)') '  Density fitting basis is "',
!    &      AUXNAM(1:ILEN),'" from the basis set library.'
!     END IF
C
C     ******************************************
C     ***** Printing of atomic coordinates *****
C     ******************************************
C
      CALL PRICAR(WORK,LWORK)
C
 1010 FORMAT(I14,' large AO-blocks of ',A1,' type')
 1020 FORMAT(I14,' small AO-blocks of ',A1,' type')
      KIPCON = 1
      KLAST  = KIPCON + KPRIM
      IF (KLAST.GT.LWORK) CALL STOPIT('BASOUT','ORBOUT',KLAST,LWORK)
      IF (.NOT. SLAVE) THEN
         CALL MOLOUT(NONTYP,NONT,IQM(1,1),NBLCK(1,1),
     &               JCO(1,1,1),NUC(1,1),NRC(1,1),SEG(1,1),
     &               WORK(KIPCON),KATOM,KANG,KBLOCK,
     &               KPRIM,CPRIMU,NRMPRI,KASYM,NSYMOP)
         CALL MOGTO(NONTYP,NONT,IQM(1,1),NBLCK(1,1),JCO(1,1,1),
     &              NUC(1,1),NRC(1,1),SEG(1,1),KATOM,KANG,KBLOCK,
     &              KPRIM,CPRIMU,NRMPRI)
      END IF
      IF (IPREAD .GE. 2) THEN
C*****************************************************************************
C        Print information about large component basis
C*****************************************************************************
         IF(NLARGE.GT.0) THEN
            IF (DIRAC) CALL HEADER('Large Components Basis',-1)
            NRMPRI = .FALSE.
            CALL ORBOUT(NONTYP,NONT,IQM(1,1),NBLCK(1,1),JCO(1,1,1),
     &                  NUC(1,1),NRC(1,1),SEG(1,1),
     &                  WORK(KIPCON),KATOM,KANG,KBLOCK,KPRIM,
     &                  CPRIMU,NRMPRI)
            IF (IPREAD .GT. 2) THEN
               NRMPRI = .TRUE.
               CALL ORBOUT(NONTYP,NONT,IQM(1,1),NBLCK(1,1),JCO(1,1,1),
     &                     NUC(1,1),NRC(1,1),SEG(1,1),
     &                     WORK(KIPCON),KATOM,KANG,KBLOCK,KPRIM,
     &                     CPRIM,NRMPRI)
            END IF
         END IF
C*****************************************************************************
C        Print information about small component basis
C*****************************************************************************
         IF (NBASIS.GT.NLARGE) THEN
            CALL HEADER('Small Components Basis',-1)
            NRMPRI = .FALSE.
            CALL ORBOUT(NONTYP,NONT,IQM(1,2),NBLCK(1,2),JCO(1,1,2),
     &                  NUC(1,2),NRC(1,2),SEG(1,2),
     &                  WORK(KIPCON),KATOM,KANG,KBLOCK,KPRIM,
     &                  CPRIMU,NRMPRI)
            IF (IPREAD .GT. 2) THEN
               NRMPRI = .TRUE.
               CALL ORBOUT(NONTYP,NONT,IQM(1,2),NBLCK(1,2),JCO(1,1,2),
     &                     NUC(1,2),NRC(1,2),SEG(1,2),
     &                     WORK(KIPCON),KATOM,KANG,KBLOCK,KPRIM,
     &                     CPRIM,NRMPRI)
            END IF
         END IF
      END IF
C
C     Symmetry basis:
C
CHJ   IF (IPREAD .GE. -1) CALL SYMOUT
CHJ april-09 bugfix: always call SYMOUT because now used for MULBS and AMFI!!!
      CALL SYMOUT
C
C     Symmetrized magnetic moments
C
      IF (HERMIT) CALL MAGCOR(IPREAD)
C
C     *************************************************
C     ***** Determine IAOAO and its inverse JAOAO in aosotr.h
C     *************************************************
C
      CALL AOTOAO(IAOAO,JAOAO,IPREAD)
C
C     **************************
C     **** Determine ICNTAO ****
C     **************************
C
      CALL CNTAO(IPREAD)
C
      RETURN
      END
C  /* Deck bastyp */
      SUBROUTINE BASTYP(IQM,JCO,NRC,NUC,NPRIM,NCONT,TSTRING)
C*****************************************************************************
C
C     Generate string : Uncontracted/contracted basis functions
C
C*****************************************************************************
#include "implicit.h"
#include "maxaqn.h"
      CHARACTER TSTRING*(6*MXAQN+3)
      CHARACTER SPDCAR*1 ! type declaration for external function
      DIMENSION JCO(IQM),NRC(*),NUC(*)
#include "ccom.h"
C*****************************************************************************
      IND = 1
      TSTRING(IND:IND) = '['
C*****************************************************************************
C * P R I M I T I V E S:
C*****************************************************************************
      JBLOCK = 0
      NPRIM = 0
      DO 100 I = 1,IQM
         IPRIM  = 0
         ISHELL = 0
         NCOMP = KHK(I)
C
         DO 200 J = 1,JCO(I)
            JBLOCK = JBLOCK + 1
            IPRIM  = IPRIM  + NUC(JBLOCK)
            NPRIM  = NPRIM  + NCOMP*NUC(JBLOCK)
  200    CONTINUE
C
         IF (IPRIM.GT.0) THEN
            NDIG = 3
C           ... i.e. max 999 s-functions per center etc.
C               for the output to be correct, cannot imagine
C               anyone having more basis functions per atom/hjaaj
            ITEN = 10**NDIG
            J = 0
            DO 300 K = NDIG,1,-1
               ITEN = ITEN/10
               ISHDIG = IPRIM/ITEN
               IF (ISHDIG .GT. 0 .OR. J .GT. 0) THEN
                  J = 1
C                 ... we want next digit, also if zero (!)
                  IDIG = ICHAR('0') + ISHDIG
                  IND = IND + 1
                  TSTRING(IND:IND) = CHAR(IDIG)
                  IPRIM=MOD(IPRIM,ITEN)
               END IF
  300       CONTINUE
            IND = IND + 1
            TSTRING(IND:IND) = SPDCAR(I-1)
         END IF
  100 CONTINUE
C*****************************************************************************
      IND = IND + 1
      TSTRING(IND:IND) = '|'
C*****************************************************************************
C * S H E L L S:
C*****************************************************************************
      JBLOCK = 0
      NCONT  = 0
      DO 400 I = 1,IQM
         ISHELL = 0
         NCOMP = KHK(I)
         DO 500 J = 1,JCO(I)
            JBLOCK = JBLOCK + 1
            ISHELL = ISHELL + NRC(JBLOCK)
            NCONT  = NCONT  + NCOMP*NRC(JBLOCK)
  500    CONTINUE
         IF(ISHELL.GT.0) THEN
            NDIG = 3
C           ... i.e. max 999 s-functions per center etc.
C               for the output to be correct, cannot imagine
C               anyone having more basis functions per atom/hjaaj
            ITEN = 10**NDIG
            J = 0
            DO 600 K = NDIG,1,-1
               ITEN = ITEN/10
               ISHDIG = ISHELL/ITEN
               IF (ISHDIG .GT. 0 .OR. J .GT. 0) THEN
                  J = 1
C                 ... we want next digit, also if zero (!)
                  IDIG = ICHAR('0') + ISHDIG
                  IND = IND + 1
                  TSTRING(IND:IND) = CHAR(IDIG)
                  ISHELL=MOD(ISHELL,ITEN)
               END IF
  600       CONTINUE
            IND = IND + 1
            TSTRING(IND:IND) = SPDCAR(I-1)
         ENDIF
  400 CONTINUE
C*****************************************************************************
      IND = IND + 1
      TSTRING(IND:IND) = ']'
      RETURN
      END
C  /* Deck gtoinp */
      SUBROUTINE GTOINP(LUINFO,IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,
     &                  CPRIMU,ISGEN,NBLOCK,KAOVEC,KPRIM)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
      PARAMETER (DSM=1.0D-30)
      CHARACTER*1 FRMT
      CHARACTER*8 FMT1
      CHARACTER*12 FMT2
C
C
      LOGICAL SEG, FREE, SEGIJ
      DIMENSION JCO(IQM),NUC(KAOVEC),NRC(KAOVEC),SEG(KAOVEC),
     &          ALPHA(KPRIM,KAOVEC),ISGEN(KAOVEC),
     &          CPRIM(KPRIM,KPRIM,KAOVEC),CPRIMU(KPRIM,KPRIM,KAOVEC)
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "cbirea.h"
#include "molinp.h"
#include "ccom.h"
C
      IWARN = 0
      NBLOCK = 0
      DO 100 I = 1, IQM
         IF (UNCONT .AND. DKTRAN .AND. JCO(I) .GT. 1) IWARN = IWARN + 1
         DO 200 J = 1, JCO(I)
            NBLOCK = NBLOCK + 1
            IF(NBLOCK.GT.KAOVEC) GOTO 5000
C
            CALL DZERO(ALPHA (1,  NBLOCK),KPRIM)
            CALL DZERO(CPRIMU(1,1,NBLOCK),KPRIM*KPRIM)
            CALL DZERO(CPRIM (1,1,NBLOCK),KPRIM*KPRIM)
C
C           Card 8:
C
C           NUC - number of uncontracted shells in AO-block
C           NRC - number of contracted shells in AO-block
C           ISGEN - bit packed information on generation of small
C                   components using the kinetic balance relation
C             1 - downwards: small L+1 from large L
C             2 - upwards  : small L-1 from large L
C             0 - if not specified same as 3 (both down and up)
C
            CALL RDLINE(LUINFO)
            READ (MLINE(NMLINE),'(BN,A1,I4,2I5)',IOSTAT=IOS)
     &           FRMT,NUCIJ,NRCIJ,ISGEN(NBLOCK)
            IF (IOS.NE.0) THEN
               WRITE(LUPRI,*) 'Error in reading .mol input, line',NMLINE
               WRITE(LUPRI,*) 'Read MLINE(NMLINE)=',MLINE(NMLINE)
               CALL QUIT('Error in reading BN,A,I4,2I5 ... ')
            ENDIF

            IF (IPREAD .GT. 3) THEN
               WRITE(LUPRI,'(3X,A,I2,A,I4)') 'L= ',(I-1),'  Block: ',J
Chj-aug99:     if (uncont) force uncontracted in ACPORB
               MRCIJ = NRCIJ
               IF (UNCONT) MRCIJ = 0
               WRITE(LUPRI,'(6X,A,I5)') 'Primitives:  ',NUCIJ,
     &                                  'Shells    :  ',MRCIJ
               CALL FLSHFO(LUPRI)
            END IF
            SEGIJ = NRCIJ .GE. 0
            NUCIJ = ABS(NUCIJ)
            NRCIJ = ABS(NRCIJ)
            IF (NUCIJ.GT.KPRIM) GOTO 5010
            IF (NRCIJ.GT.KPRIM) GOTO 5020
C
C           Read in exponents and contraction coefficients
C
            CALL ACPORB(LUINFO,FRMT,NUCIJ,NRCIJ,ALPHA(1,NBLOCK),
     &                  CPRIMU(1,1,NBLOCK),KPRIM,KAOVEC)
C
C           Identify segmented contractions
C
            CALL SEGORB(SEGIJ,NUCIJ,NRCIJ,CPRIMU(1,1,NBLOCK),KPRIM,DSM)
            SEG(NBLOCK) = SEGIJ
C
C           Reorder primitive orbitals
C
            CALL PRIORD(ALPHA(1,NBLOCK),CPRIMU(1,1,NBLOCK),NUCIJ,NRCIJ,
     &                  SEG(NBLOCK),KPRIM,DSM)
C
C           Normalize orbitals
C
            CALL NRMORB(I,NRCIJ,NUCIJ,ALPHA(1,NBLOCK),
     &                  CPRIM(1,1,NBLOCK),CPRIMU(1,1,NBLOCK),
     &                  KPRIM,NBLOCK)
C
            NUC(NBLOCK) = NUCIJ
            NRC(NBLOCK) = NRCIJ
  200    CONTINUE
  100 CONTINUE
      IF (IWARN .GT. 0) THEN
C
C     A small test to avoid Douglas-Kroll problems. However, it is not
C     sufficient, as it will not stop duplication of AOs within a single
C     block, as for instance may happen if you use Pople-type basis sets
C     K.Ruud, Nov.-2006
C
         WRITE (LUPRI,*) ' >>> WARNING GTOINP: If primitive exponents'//
     &        ' are duplicated in several blocks, ',
     &        ' >>> Douglas-Kroll will not work correctly'
      END IF
      RETURN
C
C     Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> GTOINP error, number of AO-blocks            ',NBLOCK,
     *  '                   current maximum number MXAOVC =',KAOVEC
        CALL QUIT('Too many AO-blocks')
 5010 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> GTOINP error, number of primitives per block      ',NUCIJ,
     *  '     current maximum number (use .MAXPRI to increase)  ',KPRIM
        CALL QUIT('Too many primitives')
 5020 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> GTOINP error, number of contracted functions      ',NRCIJ,
     *  '     current maximum number (use .MAXPRI to increase)  ',KPRIM
        CALL QUIT('Too many contracted functions')
      END
C  /* Deck cmbas */
      SUBROUTINE CMBAS(IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,
     &                 CPRIMU,NBLOCK,KAOVEC,KPRIM)
C
C     Basis functions according to
C     K. Kaufmann, W. Baumeister, and M. Jungen
C     "Universal Gaussian basis sets for an optimum representation
C     of Rydberg and continuum wavefunctions"
C     J. Phys. B: At. Mol. Opt. Phys. 22 (1989) 2223-2240
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
      LOGICAL SEG, FREE, SEGIJ
      DIMENSION JCO(IQM),NUC(KAOVEC),NRC(KAOVEC),SEG(KAOVEC),
     &          ALPHA(KPRIM,KAOVEC),
     &          CPRIM(KPRIM,KPRIM,KAOVEC),CPRIMU(KPRIM,KPRIM,KAOVEC),
     &          BASA(0:4), BASB(0:4)
#include "cbirea.h"
#include "molinp.h"
#include "ccom.h"
      DATA (BASA(I),I=0,4) /0.584342, 0.452615, 0.382362, 0.337027,
     &                      0.304679/
      DATA (BASB(I),I=0,4) /0.424483, 0.309805, 0.251333, 0.215013,
     &                      0.189944/
C
      NBLOCK = 0
      DO 100 I = 1, IQM
         NBLOCK = NBLOCK + 1
         IF(NBLOCK.GT.KAOVEC) GOTO 5000
C
         CALL DZERO(ALPHA (1,  NBLOCK),KPRIM)
         CALL DZERO(CPRIMU(1,1,NBLOCK),KPRIM*KPRIM)
         CALL DZERO(CPRIM (1,1,NBLOCK),KPRIM*KPRIM)
C
         NUCI = 0
         NRCI = 0
         SEGIJ  = .TRUE.
         DO 110 J = NCMSTR, NCMEND
            NUCI = NUCI + 1
            NRCI = NRCI + 1
            IF (NUCI.GT.KPRIM) GOTO 5010
            IF (NRCI.GT.KPRIM) GOTO 5020
            ALPHA(NUCI,NBLOCK) = ((ZCMVAL/J)/
     &                           ((BASA(I-1)*J)/2+BASB(I-1)))**2
            CPRIMU(NUCI,NRCI,NBLOCK) = 1.0D0
  110    CONTINUE
         NUC(NBLOCK) = NUCI
         NRC(NBLOCK) = NRCI
C
         JCO(I) = 1
         SEG(NBLOCK) = SEGIJ
C
C        Normalize orbitals
C
         CALL NRMORB(I,NRCI,NUCI,ALPHA(1,NBLOCK),CPRIM(1,1,NBLOCK),
     &               CPRIMU(1,1,NBLOCK),KPRIM,NBLOCK)
C
  100 CONTINUE
      RETURN
C
C     Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> CMBAS  error, number of AO-blocks            ',NBLOCK,
     *  '                   current maximum number MXAOVC =',KAOVEC
        CALL QUIT('CMBAS: Too many AO-blocks')
 5010 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> CMBAS  error, number of primitives per block      ',NUCI,
     *  '     current maximum number (use .MAXPRI to increase)  ',KPRIM
        CALL QUIT('CMBAS: Too many primitives')
 5020 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> CMBAS  error, number of contracted functions      ',NRCI,
     *  '     current maximum number (use .MAXPRI to increase)  ',KPRIM
        CALL QUIT('CMBAS: Too many contracted functions')
      END
C  /* Deck acporb */
      SUBROUTINE ACPORB(LUINFO,FRMT,NUCIJ,NRCIJ,ALPHA,CPRIMU,
     &                  KPRIM,KAOVEC)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "dummy.h"
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0)
      CHARACTER*1 FRMT
      CHARACTER*11 FMT1
      CHARACTER*15 FMT2
      CHARACTER*10 FMT3
      LOGICAL FREE
      DIMENSION ALPHA(KPRIM), CPRIMU(KPRIM,KPRIM)
#include "cbirea.h"
#include "molinp.h"
C
C     Card 9:
C
C     ALPHA  - exponent of primitive
C     CPRIMU - unnormalized contraction coefficients
C
      IF (FRMT .EQ. 'F' .OR. FRMT .EQ. 'f') THEN
         FREE = .TRUE.
      ELSE IF (FRMT .EQ. 'H' .OR. FRMT .EQ. 'h') THEN
         FREE = .FALSE.
         NCOL = 4
         FMT1 = '(BN,4F20.0)'
         FMT2 = '(BN,20X,3F20.0)'
         FMT3 = '(BN,F20.0)'
      ELSE
         FREE = .FALSE.
         NCOL = 8
         FMT1 = '(BN,8F10.0)'
         FMT2 = '(BN,10X,7F10.0)'
         FMT3 = '(BN,F10.0)'
      END IF
      IF(NRCIJ.NE.0) THEN
        DO 100 L = 1, NUCIJ
           CALL RDLINE(LUINFO)
           IF (FREE) THEN
#if defined (SYS_CRAY) || defined (VAR_NOFREE) || defined (SYS_T3D) || defined (SYS_T90)
              ISTART = 1
              CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,
     &                    ALPHA(L),'REA',IERR)
              DO M = 1, NRCIJ
                 CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,
     &                       CPRIMU(L,M),'REA',IERR)
              END DO
#else
              IF (WRTLIN) THEN 
                 WRITE(LUPRI,*) MLINE(NMLINE)
                 CALL FLSHFO(LUPRI)
              ENDIF
              READ (MLINE(NMLINE),*) ALPHA(L), (CPRIMU(L,M),M=1,NRCIJ)
#endif
           ELSE
              READ (MLINE(NMLINE),FMT1) ALPHA(L),
     &           (CPRIMU(L,M), M = 1, MIN(NRCIJ,NCOL - 1))
              DO 200 K = 2, (NRCIJ - 1)/(NCOL - 1) + 1
                 CALL RDLINE(LUINFO)
                 READ (MLINE(NMLINE),FMT2) (CPRIMU(L,M),
     &              M = 1 + (NCOL-1)*(K-1), MIN(NRCIJ,(NCOL-1)*K))
  200         CONTINUE
           END IF
  100   CONTINUE
C
Chj-aug99:
C     If (uncont) force uncontracted
C     (this cannot be done earlier, because we must read CPRIMU from file
C      to position LUINFO correctly)
C
          IF (UNCONT) THEN
            NRCIJ = NUCIJ
            IF (NRCIJ .GT. KAOVEC) THEN
C           in her2drv:PAOSET, each segmented is put into its own block
C           thus for UNCONT each primitive is in its own block and we
C           must have at least KAOVEC blocks
               WRITE (LUPRI,'(/A,I6/A,I6/A)')
     *      ' >>> ACPORB error, number of AO-blocks for UNCONT ',NRCIJ,
     *      '                   current maximum number MXAOVC =',KAOVEC,
     &      '  Increase MXAOVC; make depend; make'
               CALL QUIT('Too many AO-blocks for an uncontracted '//
     *                   'basis sets')
            END IF
            DO IJ = 1, NUCIJ
               DO JI = 1, IJ - 1
                  CPRIMU(IJ,JI) = D0
                  CPRIMU(JI,IJ) = D0
               END DO
               CPRIMU(IJ,IJ) = D1
            END DO
          END IF
      ELSE
C
C     Read only exponents
C
        NRCIJ = NUCIJ
        DO 300 L = 1, NUCIJ
           CALL RDLINE(LUINFO)
           IF (FREE) THEN
#if defined (SYS_CRAY) || defined (VAR_NOFREE) || defined (SYS_T3D) || defined (SYS_T90)
              ISTART = 1
              CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,
     &                    ALPHA(L),'REA',IERR)
#else
              READ (MLINE(NMLINE),*,IOSTAT=IOS) ALPHA(L)
              IF (IOS.NE.0) THEN
                WRITE(LUPRI,*)
     &          'Error in reading of ALPHA(L), NMLINE=',NMLINE,':'
                WRITE(LUPRI,*) MLINE(NMLINE)
                CALL QUIT('Error in basis set reading of ALPHA(L)')
              ENDIF
#endif
           ELSE
              READ (MLINE(NMLINE),FMT3,IOSTAT=IOS) ALPHA(L)
              IF (IOS.NE.0) THEN
                WRITE(LUPRI,*)
     &          'Error in reading of ALPHA(L), NMLINE=',NMLINE,':'
                WRITE(LUPRI,*) MLINE(NMLINE)
                CALL QUIT('Error in basis set reading of ALPHA(L)')
              ENDIF
           END IF
           CPRIMU(L,L) = D1
  300   CONTINUE
      ENDIF
C
C     Test print
C
      IF (IPREAD .GT. 3) THEN
         WRITE(LUPRI,'(A,I4,A,I4,A)')
     &   'Contraction matrix (',NUCIJ,'x',NRCIJ,') :' 
         DO 400 L = 1, NUCIJ
            WRITE (LUPRI,'(1P,G20.10)') ALPHA(L)
            WRITE (LUPRI,'(10F12.8)') (CPRIMU(L,M),M=1,NRCIJ)
  400    CONTINUE
      END IF
      RETURN
      END
C  /* Deck segorb */
      SUBROUTINE SEGORB(SEGIJ,NUCIJ,NRCIJ,CPRIMU,KPRIM,DSM)
#include "implicit.h"
#include "priunit.h"
      LOGICAL SEGIJ
      DIMENSION CPRIMU(KPRIM,KPRIM)
C
      DO L = 1, NUCIJ
         NONZER = 0
         DO M = 1, NRCIJ
            IF (ABS(CPRIMU(L,M)).GT.DSM) NONZER = NONZER + 1
         END DO
         SEGIJ = SEGIJ .AND. (NONZER .LE. 1)
      END DO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck nrmorb */
      SUBROUTINE NRMORB(I,NRCIJ,NUCIJ,ALPHA,CPRIM,CPRIMU,KPRIM,NBLOCK)
#include <implicit.h>
#include <priunit.h>
#include <pi.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D4 = 4.0D0,
     &           DP25 = 0.25D0, DP5 = 0.5D0, DP75 = 0.75D0,
     &           THRMIN = 1.D-17, EXPMIN = 1.D-10)
      CHARACTER*1 SPDCAR
      DIMENSION ALPHA(KPRIM), CPRIM(KPRIM,KPRIM), CPRIMU(KPRIM,KPRIM)
C
      PIPPI = (DP5/PI)**DP75
      IF (NUCIJ.EQ.1 .AND. ABS(ALPHA(1)).LE.EXPMIN) THEN
         CPRIM (1,1) = D1 
         CPRIMU(1,1) = D1 
      ELSE 
         DO N = 1, NRCIJ
            SUM = D0
            DO L = 1, NUCIJ
            DO M = 1, NUCIJ
               T = D2*SQRT(ALPHA(L)*ALPHA(M))/(ALPHA(L)+ALPHA(M))
               SUM = SUM + CPRIMU(L,N)*CPRIMU(M,N)*(T**(I + DP5))
            END DO
            END DO
            IF (SQRT(SUM) .LT. THRMIN) GOTO 1000
            SUM=D1/SQRT(SUM)
            DO L=1, NUCIJ
            CPRIM(L,N)=CPRIMU(L,N)*SUM*(D4*ALPHA(L))**(DP5*I+DP25)*PIPPI
            END DO
         END DO
      END IF
      RETURN
 1000 CONTINUE
         WRITE (LUPRI,'(/A,2(I3,A),A1,A/A,2I4)')
     &    ' INPUT ERROR: CGTO no.',N,' for block',NBLOCK,
     &    ' of ',SPDCAR(I-1),' type has zero norm.',
     &    ' Contraction matrix dimensions: ',NUCIJ,NRCIJ
         CALL OUTPUT(CPRIMU,1,NUCIJ,1,NRCIJ,KPRIM,KPRIM,-1,LUPRI)
         CALL QUIT('NRMORB: CGTO with zero norm.')
      END
C  /* Deck orbpro */
      SUBROUTINE ORBPRO(NONTYP,NONT,IQM,NBLCK,JCO,NUC,NRC,SEG,
     &           ALPHA,CPRIM,CPRIMU,KATOM,KANG,KBLOCK,KPRIM,
     &           ISHELL,IPRIM,IPRIMD,IORB,IORBD,LCOMP)
C******************************************************************************
C
C     Process basis data
C
C******************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <pi.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0, D4 = 4.0D0, DP5 = 0.50D0,
     &           DP75 = 0.75D0)
C
#include <cbirea.h>
      LOGICAL SEG,SPHER
      DIMENSION NONT(KATOM),IQM(KATOM),NBLCK(KATOM),
     &          JCO(KANG,KATOM),NUC(KBLOCK),NRC(KBLOCK),
     &          SEG(KBLOCK),ALPHA(KPRIM,KBLOCK),
     &          CPRIM(KPRIM,KPRIM,KBLOCK),CPRIMU(KPRIM,KPRIM,KBLOCK)
#include <ccom.h>
#include <nuclei.h>
#include <primit.h>
#include <shells.h>
#include <symmet.h>
#include <r12int.h>
#include <aosotr.h>
#include <gnrinf.h>
C
C
      IFSTSHAUX = 0
      ICENT  = 0
      JBLOCK = 0
      DO 10 I = 1,NONTYP ! number of atomic types
         DO 20 N = 1,NONT(I) ! number of symmetry independent centers of
                             ! this type
            ICENT = ICENT + 1
            KBCH = JBLOCK
            NDEG = NUCDEG(ICENT)
            DO 30 J = 1,IQM(I) ! ang.mom. 1=s, 2=p, 3=d, etc.
               KKK = 0
               NCOMP = KHK(J)
               NCCMP = KCK(J)
               SPHER = SPH(J)
               DO 40 K = 1, JCO(J,I) ! number of CGTOs in this AO block for this ang.mom.
                  KBCH  = KBCH + 1 ! next block
                  NCONT = NRC(KBCH)
                  IF(NCONT.GT.MXCONT) GOTO 5000
                  DO 50 KK = 1, NCONT ! AO shell conunter for this block
                     KKK    = KKK + 1 ! local AO shell counter for this ang.mom.
                     ISHELL = ISHELL + 1 ! global AO shell couunter
                     !set the pointer to the first shell for density-fitting
                     !auxiliary basis
                     IF ( (LCOMP.EQ.3) .AND. (IFSTSHAUX.EQ.0) ) THEN
                        IFSTSHAUX = ISHELL
                     END IF
C*****************************************************************************
C     Data on the AO-block associated with a given shell(ISHELL):
C       NUCO    - number of uncontracted functions
C       NRCO    - number of contracted functions
C     Data on a given shell in an AO-block:
C       NCENT   - index of symmetry independent center
C       NUMCF   - index of contracted function (shell) in AO-block
C       NUMCFT  - total index of contracted function (shell) in AO-blocks of this ang.mom.
C       NBCH    - index of block in AO-vector
C       ISTBAO  - stabiliser: basic sym. op. that do not move center
C       NHKT    - angular quantum number (s=1,p=2,d=3 etc.)
C       KHKT    - number of spherical (Cartesian) components
C       KCKT    - number of Cartesian components
C       SEGM    - segmented contraction
C       LCLASS  - class: large component (1), small or Huckel(2), density fitting (0)
C       CENT    - coordinates of center
C*****************************************************************************
                     NUCO  (ISHELL) = NUC(KBCH)
                     NRCO  (ISHELL) = NCONT
                     JSTRT (ISHELL) = IPRIM
                     NSTRT (ISHELL) = IORB
                     KSTRT (ISHELL) = IORBD
                     NCENT (ISHELL) = ICENT
                     NUMCF (ISHELL) = KK
                     NUMCFT(ISHELL) = KKK
                     NBCH  (ISHELL) = KBCH
                     SHARE (ISHELL) = .FALSE.
                     ISTBAO(ISHELL) = ISTBNU(ICENT)
                     NHKT  (ISHELL) = J
                     KHKT  (ISHELL) = NCOMP
                     SPHR  (ISHELL) = SPHER
                     KCKT  (ISHELL) = NCCMP
                     SEGM  (ISHELL) = SEG(KBCH)
                     LCLASS(ISHELL) = LCOMP
                     CENT(ISHELL,1,1) = CORD(1,ICENT)
                     CENT(ISHELL,2,1) = CORD(2,ICENT)
                     CENT(ISHELL,3,1) = CORD(3,ICENT)
C                    Basis-set identifier (WK/UniKa/04-11-2002).
                     MBSID(ISHELL) = MULBSI(ICENT)

                     IORB  = IORB  + NCOMP
                     IORBD = IORBD + NCOMP*NDEG
                     
                     IF (.NOT. SEGM(ISHELL) .AND. LCOMP .EQ. 1)
     &                         SEGBAS = .FALSE.
                     IF (.NOT. SEGM(ISHELL) .AND. LCOMP .EQ. 3)
     &                         SEGAUX = .FALSE.
   50             CONTINUE ! DO 50 KK = 1, NCONT
                  DO 60 L = 1, NUC(KBCH)
                     IPRIM = IPRIM + 1
C*****************************************************************************
C       PRIEXP  - exponent of primitive shell
C       PRICCF  - normalized contraction coefficient
C       PRICRX  - x-coordinate of center
C       PRICRY  - y-coordiante of center
C       PRICRZ  - z-coordinate of center
C*****************************************************************************
                     PRIEXP(IPRIM) = ALPHA(L,KBCH)
                     DO 70  M = 1, NCONT
                        PRICCF(IPRIM,M) = CPRIM(L,M,KBCH)
   70                CONTINUE
                     PRICRX(IPRIM) = CORD(1,ICENT)
                     PRICRY(IPRIM) = CORD(2,ICENT)
                     PRICRZ(IPRIM) = CORD(3,ICENT)
                     IPRSHL(IPRIM) = ISHELL
                     IPRIMD = IPRIMD + NCOMP*NDEG
C
   60             CONTINUE ! DO 60 L = 1, NUC(KBCH)
   40          CONTINUE ! DO 40 K = 1, JCO(J,I)
   30       CONTINUE ! DO 30 J = 1,IQM(I)
   20    CONTINUE ! DO 20 N = 1,NONT(I)
         JBLOCK = JBLOCK + NBLCK(I)
   10 CONTINUE ! DO 10 I = 1,NONTYP
C
      !Add exp(0) function (needed in 2- and 3-index integrals)
      !as an extra "0" class for ATOMDF case
      IF (LCOMP .EQ. 3) THEN
         ISHELL = ISHELL + 1
         NUCO  (ISHELL) = 1
         NRCO  (ISHELL) = 1
         JSTRT (ISHELL) = IPRIM
         NSTRT (ISHELL) = IORB
         KSTRT (ISHELL) = IORBD
         NCENT (ISHELL) = 0
         NUMCF (ISHELL) = 0
         NUMCFT(ISHELL) = 0
         NBCH  (ISHELL) = 0
         SHARE (ISHELL) = .FALSE.
!        ISTBAO(ISHELL) = ISTBNU(ICENT)
         NHKT  (ISHELL) = 1
         KHKT  (ISHELL) = 1
         SPHR  (ISHELL) = .FALSE.
         KCKT  (ISHELL) = 1
         SEGM  (ISHELL) = .FALSE.
         LCLASS(ISHELL) = 0
         CENT(ISHELL,1,1) = 0.0D0
         CENT(ISHELL,2,1) = 0.0D0
         CENT(ISHELL,3,1) = 0.0D0
C        Basis-set identifier (WK/UniKa/04-11-2002).
*        MBSID(ISHELL) = MULBSI(ICENT)
*        IORB  = IORB  + NCOMP
*        IORBD = IORBD + NCOMP*NDEG
         IPRIM = IPRIM + 1
C*****************************************************************************
C       PRIEXP  - exponent of primitive shell
C       PRICCF  - normalized contraction coefficient
C       PRICRX  - x-coordinate of center
C       PRICRY  - y-coordiante of center
C       PRICRZ  - z-coordinate of center
C*****************************************************************************
        PRIEXP(IPRIM) = 0.0D0
        PRICCF(IPRIM,1) = 1.0D0
        PRICRX(IPRIM) = 0.0d0
        PRICRY(IPRIM) = 0.0d0
        PRICRZ(IPRIM) = 0.0d0
        IPRSHL(IPRIM) = ISHELL
        IPRIMD = IPRIMD + 1
      END IF

      RETURN
C
C     Error messages:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *   ' >>> ORBPRO error, '//
     *   'no. of contracted functions in an AO-block',NCONT,
     *   '                   '//
     *   'current maximum number           MXCONT = ',MXCONT
      CALL QUIT('ORBPRO: Too many contracted functions per AO-block')
      END
C  /* Deck orbout */
      SUBROUTINE ORBOUT(NONTYP,NONT,IQM,NBLCK,JCO,NUC,NRC,SEG,
     &           IPCON,KATOM,KANG,KBLOCK,KPRIM,CPRIMU,NRMPRI)
#include "implicit.h"
#include "priunit.h"
#include "pi.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0)
C
#include "ccom.h"
#include "cbirea.h"
#include "nuclei.h"
#include "primit.h"
      CHARACTER*10 CHRSEG
      LOGICAL SEG, NRMPRI
      DIMENSION NONT(KATOM),IQM(KATOM),NBLCK(KATOM),
     &          JCO(KANG,KATOM),NUC(KBLOCK),NRC(KBLOCK),
     &          CPRIMU(KPRIM,KPRIM,KBLOCK),
     &          SEG(KBLOCK)
      DIMENSION IPCON(KPRIM)
C
C     P R I M I T I V E S
C     ===================
C
      IF (.NOT.NRMPRI) THEN
        CALL HEADER('Orbital exponents and contraction coefficients',1)
      ELSE
        CALL HEADER
     &  ('Orbital exponents and normalized contraction coefficients',1)
      END IF
      IPRIMD = 0
      IPRIM =  0
      ICENT  = 0
      JBLOCK = 0
      DO 100 I = 1, NONTYP
         DO 110 N = 1, NONT(I)
            ICENT = ICENT + 1
            NDEG  = NUCDEG(ICENT)
            KBCH  = JBLOCK
            DO 200 J = 1, IQM(I)
            DO 200 K = 1, JCO(J,I)
               KBCH = KBCH + 1
               NNUC  = NUC(KBCH)
               NNRC  = NRC(KBCH)
            IF (NNUC .EQ. 0) GO TO 200
               ITYP = NHKOFF(J)
               IPSTRT = IPRIM + 1
               IPRIM =  IPRIM + NNUC
               IF (SEG(KBCH)) THEN
                  CHRSEG = 'seg. cont.'
               ELSE
                  CHRSEG = 'gen. cont.'
               END IF
               IF (NNRC .LE. 5) THEN
                  DO 300 ICOMP = 1, KHK(J)
                     ITYP = ITYP + 1
                     DO 330 L = 1, NDEG
                        IPRIMD = IPRIMD + 1
                        IF (NDEG .GT. 1) THEN
                           WRITE (LUPRI,1000) NAMN(ICENT),'#',L,
     &                          GTOTYP(ITYP),
     &                          IPRIMD,PRIEXP(IPSTRT),
     &                          (CPRIMU(1,MM,KBCH),MM=1,NNRC)
                        ELSE
                           WRITE (LUPRI,1010) NAMN(ICENT),GTOTYP(ITYP),
     &                          IPRIMD,PRIEXP(IPSTRT),
     &                          (CPRIMU(1,MM,KBCH),MM=1,NNRC)
                        END IF
                        IF (NNUC .GT. 1) THEN
                           IPRIMD = IPRIMD + 1
                           WRITE (LUPRI,1020) CHRSEG,
     &                          IPRIMD,PRIEXP(IPSTRT+1),
     &                          (CPRIMU(2,MM,KBCH),MM=1,NNRC)
                        END IF
                        DO 350 M = 3, NNUC
                           IPRIMD = IPRIMD + 1
                           WRITE (LUPRI,1030) IPRIMD,
     &                          PRIEXP(IPSTRT-1+M),
     &                          (CPRIMU(M,MM,KBCH),MM=1,NNRC)
 350                    CONTINUE
 330                  CONTINUE
 300               CONTINUE
               ELSE
                  DO 400 ICOMP = 1, KHK(J)
                     ITYP = ITYP + 1
                     DO 430 L = 1, NDEG
                        IPRIMD = IPRIMD + 1
                        IF (NDEG .GT. 1) THEN
                           WRITE (LUPRI,1005) NAMN(ICENT),'#',L,
     &                          GTOTYP(ITYP),
     &                          IPRIMD,PRIEXP(IPSTRT),
     &                          (CPRIMU(1,MM,KBCH),MM=1,NNRC)
                        ELSE
                           WRITE (LUPRI,1015) NAMN(ICENT),GTOTYP(ITYP),
     &                          IPRIMD,PRIEXP(IPSTRT),
     &                          (CPRIMU(1,MM,KBCH),MM=1,NNRC)
                        END IF
                        IF (NNUC .GT. 1) THEN
                           IPRIMD = IPRIMD + 1
                           WRITE (LUPRI,1025) CHRSEG,
     &                          IPRIMD,PRIEXP(IPSTRT+1),
     &                          (CPRIMU(2,MM,KBCH),MM=1,NNRC)
                        END IF
                        DO 450 M = 3, NNUC
                           IPRIMD = IPRIMD + 1
                           WRITE (LUPRI,1035) IPRIMD,
     &                          PRIEXP(IPSTRT-1+M),
     &                          (CPRIMU(M,MM,KBCH),MM=1,NNRC)
 450                    CONTINUE
 430                 CONTINUE
 400              CONTINUE
               END IF
 200        CONTINUE
 110     CONTINUE
         JBLOCK = JBLOCK + NBLCK(I)
 100  CONTINUE
 1000 FORMAT(/2X,A4,A1,I1,1X,A4,I5,F16.6,2X,5F10.4)
 1010 FORMAT(/2X,A4,3X,      A4,I5,F16.6,2X,5F10.4)
 1020 FORMAT( 3X,A10,           I5,F16.6,2X,5F10.4)
 1030 FORMAT(13X,               I5,F16.6,2X,5F10.4)
 1005 FORMAT(/2X,A4,A1,I1,1X,A4,I5,F16.6,2X,5F10.4/,(36X,5F10.4))
 1015 FORMAT(/2X,A4,3X,      A4,I5,F16.6,2X,5F10.4/,(36X,5F10.4))
 1025 FORMAT( 3X,A10,           I5,F16.6,2X,5F10.4/,(36X,5F10.4))
 1035 FORMAT(13X,               I5,F16.6,2X,5F10.4/,(36X,5F10.4))
C
      IF (NRMPRI) RETURN
C
C       C O N T R A C T I O N
C       =====================
C
C       Looping is over
C                       - atomic type
C                         - symmetry independent center
C                           - shell
C                             - component
C                               - degeneracy of center
C
      CALL HEADER('Contracted Orbitals',1)
      IPRIMD  = 0
      IORBD   = 0
      IPRIM   = 0
      IORB    = 0
      ICENT   = 0
      JBLOCK  = 0
      DO 500 I = 1,NONTYP
         DO 510 N = 1,NONT(I)
            ICENT = ICENT + 1
            NDEG  = NUCDEG(ICENT)
            KBCH  = JBLOCK
            DO 600 J = 1,IQM(I)
               NCOMP = KHK(J)
               DO 610 K = 1,JCO(J,I)
                  KBCH = KBCH + 1
                  NNUC = NUC(KBCH)
                  NNRC  = NRC(KBCH)
                  IPSTRT = IPRIM + 1
                  IPRIM  = IPRIM + NNUC
                  IORB = IORB + NNRC
                  DO 700 L = 1,NNRC
                     JPRIM = 0
                     JPRIMD = IPRIMD
                     DO 750 M = IPSTRT,IPRIM
                        JPRIMD = JPRIMD + 1
                        IF(PRICCF(M,L).NE.0.00D0) THEN
                           JPRIM = JPRIM + 1
                           IPCON(JPRIM) = JPRIMD
                        ENDIF
  750                CONTINUE
                     ITYP = NHKOFF(J)
                     IOFF = 0
                     DO 800 ICOMP = 1, NCOMP
                        ITYP = ITYP + 1
                        IF(NDEG.GT.1) THEN
                           DO 810 LL = 1,NDEG
                              IORBD = IORBD + 1
                              WRITE(LUPRI,1040)
     &                           IORBD,NAMN(ICENT),'#',LL,GTOTYP(ITYP),
     &                           ((IPCON(M)+IOFF),M=1,JPRIM)
                              IOFF = IOFF + NNUC
  810                      CONTINUE
                        ELSE
                           IORBD = IORBD + 1
                           WRITE(LUPRI,1050) IORBD,NAMN(ICENT),
     &                        GTOTYP(ITYP),((IPCON(M)+IOFF),M=1,JPRIM)
                           IOFF = IOFF + NNUC
                        ENDIF
  800                CONTINUE
  700             CONTINUE
                  IPRIMD = IPRIMD + NNUC*NDEG*NCOMP
  610          CONTINUE
  600       CONTINUE
  510    CONTINUE
         JBLOCK = JBLOCK + NBLCK(I)
  500 CONTINUE
      WRITE (LUPRI,'(/)')
 1040 FORMAT(I5,2X,A4,A1,I1,2X,A4,24I5)
 1050 FORMAT(I5,2X,A4,4X,A4,24I5)
C
      RETURN
      END
C  /* Deck nucpro */
      SUBROUTINE NUCPRO(WORK,LWORK)
C******************************************************************************
C
C     Process molecular data
C
C******************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "pi.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0)
      DIMENSION WORK(LWORK)
C
#include "gnrinf.h"
#include "cbirea.h"
#include "ccom.h"
#include "nuclei.h"
#include "frame.h"
#include "primit.h"
#include "shells.h"
#include "symmet.h"
#include "aosotr.h"
#include "pcmlog.h"
#include "qm3.h"
C
C
#include "ibtfun.h"
C
C     Statement function to "convert" 0 - 7 to a printable
C     "binary" form
C
      IFAKBN(I) = 25*IBTAND(I,4) + 5*IBTAND(I,2) + IBTAND(I,1)

      CALL QENTER('NUCPRO')
C
C     *************************************************************
C     ***** Calculate Distances ,Nuclear Potential Energy and *****
C     ***** Nuclear Contributions to Dipole Moment            *****
C     * LOOPS:
C     * 10      - loop over first symmetry independent centre
C     * 20+IF   - loop over dependent centers
C     * 30      - loop over second symmetry independent centre
C     * 40+IF   - loop over dependent centers
C     ************************************************************
C
      DPNUC(1) = D0
      DPNUC(2) = D0
      DPNUC(3) = D0
      NUCDEP    = 0
      NCDTOT    = 0
      NUCPRE(1) = 0
      NFLOAT    = 0

      DO 100 N = 1, NCTOT
         IVARB = ISTBNU(N)
         NDEG  = MULT(IVARB)
         IF (IPREAD .GE. 5)
     &      WRITE (LUPRI,1000) NAMN(N),(CORD(I,N),I=1,3),CHARGE(N),
     &                        IFAKBN(ISTBNU(N)), NDEG - 1
         NUCDEG(N) = NDEG
         NUCPRE(N + 1) = NUCPRE(N) + NDEG
         IF (CHARGE(N) .EQ. D0 .AND. NAMN(N) .NE. 'cav' .AND.
     &       N .LE. NUCIND) THEN
            NFLOAT = NFLOAT + NDEG
         ELSE IF (NFLOAT .GT. 0) THEN
            CALL QUIT(
     &      'All nuclei MUST come before floating orbital centers')
!           So 1:NATOMS run over the real nuclei
         END IF
         IF ( ISUBSY(N) .EQ. 0 ) THEN
C        ... only QM atoms contribute to nuclear dipole moment
            DPNUC(1) = DPNUC(1) + CORD(1,N)*CHARGE(N)*FMULT(IVARB)
            DPNUC(2) = DPNUC(2) + CORD(2,N)*CHARGE(N)*FMULT(IVARB)
            DPNUC(3) = DPNUC(3) + CORD(3,N)*CHARGE(N)*FMULT(IVARB)
         END IF
         II = 0
         DO 200 LA = 0,MAXREP
            IF (IBTAND(IVARB,LA) .EQ. 0) THEN
               II = II + 1
               IF (N .LE. NUCIND) NUCDEP = NUCDEP + 1
               NCDTOT = NCDTOT + 1
               NUCNUM(N,LA+1) = NCDTOT
               IF (NDEG .EQ. 1 ) THEN
                  NAMDEP(NCDTOT)     = NAMN(N)//'  '
                  NAMDPX(3*NCDTOT-2) = NAMN(N)//'   x'
                  NAMDPX(3*NCDTOT-1) = NAMN(N)//'   y'
                  NAMDPX(3*NCDTOT  ) = NAMN(N)//'   z'
                  NDEGNM(NCDTOT)     = 1
               ELSE
                  ICHARD = ICHAR('0') + II
                  NAMDEP(NCDTOT)     = NAMN(N)//'_'//CHAR(ICHARD)
                  NAMDPX(3*NCDTOT-2) = NAMDEP(NCDTOT)//' x'
                  NAMDPX(3*NCDTOT-1) = NAMDEP(NCDTOT)//' y'
                  NAMDPX(3*NCDTOT  ) = NAMDEP(NCDTOT)//' z'
                  NDEGNM(NCDTOT)     = II
               END IF
            ELSE
               NUCNUM(N,LA+1) = 0
            END IF
  200    CONTINUE
  100 CONTINUE
C     Now: NUCDEP = QM atoms + floating centers
C          NCDTOT = NUCDEP + MM centers
      NATOMS = NUCDEP - NFLOAT
      IF (QM3) NUCDEP = NCDTOT
      IF (NCDTOT .GT. MXCENT) GOTO 5010
C
C     Nuclear repulsion energy
C     ========================
C
      POTNUC = D0
      DO 300 N = 1, NUCIND
        DO 310 M = N, NUCIND
          DO 320 KB = 0, MAXREP
            IF (IBTAND(KB,ISTBNU(M)) .EQ. 0) THEN
              IF (M.EQ.N .AND. KB.EQ.0) GO TO 320
              DS = (CORD(1,N)-CORD(1,M)*PT(IBTAND(ISYMAX(1,1),KB)))**2
     &           + (CORD(2,N)-CORD(2,M)*PT(IBTAND(ISYMAX(2,1),KB)))**2
     &           + (CORD(3,N)-CORD(3,M)*PT(IBTAND(ISYMAX(3,1),KB)))**2
              DIST = SQRT(DS)
              IF (DIST.LT.0.1.AND.CHARGE(M)*CHARGE(N).NE.D0) GOTO 5000
              IF (IBTAND(KB,ISTBNU(N)) .EQ. 0) THEN
                IF (DIST.GT.D0) THEN
                  HKAB  = FMULT(IBTAND(ISTBNU(M),ISTBNU(N)))
                  IF (M .EQ. N) HKAB = DP5*HKAB
                  POTNUC = POTNUC + CHARGE(M)*CHARGE(N)*HKAB/DIST
                END IF
              END IF
            END IF
  320     CONTINUE
  310   CONTINUE
  300 CONTINUE
C
C     **********************************************************
C     ***** Calculate symmetry-adapted nuclear coordinates *****
C     **********************************************************
C
      CALL SYMNCO(LWORK,WORK)
C
      CALL QEXIT('NUCPRO')
      RETURN
C
C     Error messages:
C
 5000 CONTINUE
      WRITE (LUPRI,'(3(A,I5)/A,1P,D15.5/)')
     &     ' Nucleus no.',N,' is too close to the',KB,
     &     'th transformation of nucleus no.',M,' DISTANCE =',DIST
      WRITE (LUPRI,'(A,2F30.10)')
     &     ' X_A, X_B :',CORD(1,N),CORD(1,M)*PT(IBTAND(ISYMAX(1,1),KB)),
     &     ' Y_A, Y_B :',CORD(2,N),CORD(2,M)*PT(IBTAND(ISYMAX(2,1),KB)),
     &     ' Z_A, Z_B :',CORD(3,N),CORD(3,M)*PT(IBTAND(ISYMAX(3,1),KB))
      CALL QUIT('NUCPRO: Nuclei too close.')
 5010 CONTINUE
      WRITE (LUPRI,'(///A/2(A,I5/)/A,I5,A)')
     &      ' Allowed number of atoms exceeded.',
     &      ' - Number of atoms in input:',NUCDEP,
     &      ' - Number of atoms allowed: ',MXCENT,
     &      ' Increase MXCENT to at least',NUCDEP,' and recompile.'
      CALL QUIT('NUCPRO: Too many atoms for current value of MXCENT.')
C
 1000 FORMAT(2X,A4,5X,3F10.5,F20.8/,
     &       16X, 'Stabilizer',I3,', with',I2,
     &       ' symmetry equivalent atoms',/)
      END
C  /* Deck sympro */
      SUBROUTINE SYMPRO(LVALUE,MVALUE,NVALUE,IRREP,DOOWN)
C*****************************************************************************
C
C     SYMPRO generates the symmetry orbitals (SOs) from the list of
C     symmetry independent atomic orbitals(AOs).
C     A given SO of irrep L is given as:
C
C         SO(L,a) = SUM{i} [CTRAN(NSORB,i)*G(i)*AO(NSORB)]
C
C     Here {G(i)} is a set of symmetry operators transforming AO(NSORB)
C     into all corresponding orbitals centered on symmetry related sites
C     in the molecule. The loop over i is only allowed to encounter
C     elements that do not contain any basic operations belonging to the
C     stabiliser of the centre; thus any G uniquely defines transformation
C     to a given centre.
C
C     By allowing all symmetry operations to work on the list of
C     symmetry independent AOs a list of symmetry dependent AOs
C     is generated (counted by NAORB). Redundancies are removed by
C     limiting symmetry operations to the set G(i) and the list thus
C     reduced to the list of non-trivial symmetry dependent AOs
C     (counted by NSORB).
C
C******************************************************************************
C
C     This subroutine determines the following parameters:
C
C     ISYMAO  -  indicates how an atomic orbital behaves under the basic
C                symmetry operations. When a basic operation changes the
C                sign of the AO (centered at origo) the corresponding bit
C                gets the value 1.
C
C     CTRAN   -  transformation coefficients (+1 or -1)
C     IPIND   -  packed integer
C     IPTSYM  -  pointer from redundant list of symmetry dependent AOs
C                to list of non-trivial symmetry dependent AOs. The
C                pointer is zero if the AO does not contribute to the
C                irrep.
C     ITRAN   -  given an operation G(i): pointer from list of symmetry
C                dependent AOs to list of symmetry independent AOs
C     JPRX    -  AO label
C     JTRAN   -  index of operation G moving a given AO to another centre
C     MAMN    -  name of AO centre
C     NPARNU  -  offset pointer from non-symmetric operators for given
C                irrep
C     NPARSU  -  offset pointer for symmetry dependent AOs for given irrep
C     NSORB   -  number of symmetry dependent AOs
C*************************************************************************
C
C     The code originally resided inside READIN
C     tuh 120988
C
C     Polished and annotated - tsaue March 10 1993
C*************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      PARAMETER (D0 = 0.0D0)
C
      LOGICAL DOOWN
      DIMENSION LVALUE(MXAQN), MVALUE(MXAQN), NVALUE(MXAQN),
     &          IRREP(MXCORB)
#include "nuclei.h"
#include "ccom.h"
#include "shells.h"
#include "symmet.h"
#include "pincom.h"
#include "aosotr.h"
#include "sphtrm.h"
#include "cbirea.h"
#include "huckel.h"

#ifdef MOD_ECP
#include "argoscom.h"
#endif

C
#include "ibtfun.h"
C
C*****************************************************************************
C Determine ISYMAO:
C       ISYMAO is a bitstring indicating how atomic orbitals centered in
C       origo behaves under the basic symmetry operations
C*****************************************************************************
C
C     I. Cartesian functions
C     ----------------------
C
      IF (DOCART) THEN
         DO 100 I = 1, NHTYP
            CALL LMNVAL(I,I*(I + 1)/2,LVALUE,MVALUE,NVALUE)
            DO 100 J = 1, I*(I + 1)/2
                LVAR = MOD(LVALUE(J),2)*ISYMAX(1,1)
                MVAR = MOD(MVALUE(J),2)*ISYMAX(2,1)
                NVAR = MOD(NVALUE(J),2)*ISYMAX(3,1)
                ISYMAO(I,J) = IBTXOR(LVAR,IBTXOR(MVAR,NVAR))
  100    CONTINUE
C
C     II. Your own scheme
C     -------------------
C
      ELSE IF (DOOWN) THEN
         IJK = 0
         DO 200 I = 1, NHTYP
            DO 210 K = 1, KHK(I)
               ISYMAO(I,K) = -1
  210       CONTINUE
            CALL LMNVAL(I,I*(I + 1)/2,LVALUE,MVALUE,NVALUE)
            DO 220 J = 1, I*(I + 1)/2
                LVAR = MOD(LVALUE(J),2)*ISYMAX(1,1)
                MVAR = MOD(MVALUE(J),2)*ISYMAX(2,1)
                NVAR = MOD(NVALUE(J),2)*ISYMAX(3,1)
                DO 220 K = 1, KHK(I)
                   IJK = IJK + 1
                   IF (ABS(CSP(IJK)).GT.D0) THEN
                      IF (ISYMAO(I,K) .EQ. -1) THEN
                         ISYMAO(I,K) = IBTXOR(LVAR,IBTXOR(MVAR,NVAR))
                      ELSE
                         WRITE (LUPRI,'(/A/A/A,2I5)')
     &                   ' Incorrect Cartesian transformation input:',
     &                   '  Components belonging to different irreps '
     &                   //'have been combined.',
     &                   '  Angular quantum number and component:',I-1,K
                         CALL QUIT('Error in Cartesian->own bf input')
                      END IF
                   END IF
  220       CONTINUE
  200    CONTINUE
      ELSE
C
C     III. Spherical harmonics
C     ------------------------
C
         DO 300 I = 0, NHTYP - 1
            II = I + 1
            IF (I .EQ. 0) THEN
               ISYMAO(II,1) = IREPLM(0,0)
            ELSE IF (I .EQ. 1) THEN
               ISYMAO(II,1) = IREPLM(1, 1)
               ISYMAO(II,2) = IREPLM(1,-1)
               ISYMAO(II,3) = IREPLM(1, 0)
            ELSE
               DO 310 J = -I,I
                  JJ = J + I + 1
                  ISYMAO(II,JJ) = IREPLM(I,J)
  310          CONTINUE
            END IF
  300    CONTINUE
      END IF
C*****************************************************************************
C Determine transformation from AOs to SOs :
C   * DO-loops:
C       400     - run over irreducible representations
C       410     - run over shells
C       420     - run over components, thus (60+70):run over orbitals
C       430     - run over the symmetry operations that have no part in the
C               stabilizer of the center, that is run over a unique set
C               of symmetry operations transforming between symmetry
C               dependent centers
C   * Counting variables:
C       NSORB   - SO orbital number
C       NAORBD  - AO orbital number
C       NAORB   - AO orbital shell number
C       NAA     -
C       JKB     -
C       IRREPN  -
C*****************************************************************************
      IHNMAO    = 0
      NPARSU(1) = 0
      CALL IZERO(NAOS ,8)
      CALL IZERO(NCOS ,8*MXBSETS_TOT)
      CALL DZERO(CTRAN,8*MXCORB)
      CALL IZERO(ITRAN,8*MXCORB)

      CALL IZERO(IRREP,MXCORB)
C
C     Run over Irreducible Representations LA (for LAMBDA)
C
      NSORB     = 0 ! counter for total SO number over all symmetries
      DO 400 LA = 0, MAXREP
         NAORBD = 0 ! counter for symmetry dependent   total AO number
         NAORB  = 0 ! counter for symmetry independent total AO number
         DO 410 IA = 1, KMAX
            MULA   = ISTBAO(IA)
            NHKTA  = NHKT(IA)
            KHKTA  = KHKT(IA)
            IC     = LCLASS(IA)
            IF (DIRAC) THEN
               JC = IC
            ELSE IF (IC .EQ. ISETHUCKEL) THEN
C              ... Huckel basis
               JC = 0
            ELSE IF (IC .EQ. 0) THEN
C              ... Density fitting basis
               JC = -1
            ELSE IF (LMULBS) THEN
C           IC .eq. 1 .and. LMULBS, and
C           JC is then basis-set identifier for multiple basis sets (WK/UniKA/04-11-2002).
               JC = MBSID(IA)
            ELSE 
               JC = 1
            END IF
            IF (JC .GT. MXBSETS .OR. JC .LT. -MXBSETS) THEN
C              Dimensions of COMP(-MXBSETS:MXBSETS) and NCOS(8,-MXBSETS:MXBSETS) and ? must be changed
C              for abs(JC) .gt. MXBSETS /hjaaj Oct 2007
               CALL QUIT('abs(JC) .gt. MXBSETS is not implemented')
            END IF
            DO 420 NA = 1, KHKTA ! loop over m_l components of this orbital shell
               NAORB = NAORB + 1
               IVARB = IBTXOR(LA,ISYMAO(NHKTA,NA))
C
C              If orbital contributes to this representation:
C
c                (IVARB is a bitstring of basic operations with
C                1 in the positions where the basic operations
C                has a different parity for the irrep and the
C                AO(when centered in origo). These basic
C                operations with different parities must then
C                not be part of the stabilizer of the center
C                for the AO).
C
               IF (IBTAND(MULA,IVARB) .NE. 0) THEN
                  IPTSYM(NAORB,LA) = 0
                  NAORBD = NAORBD + MULT(MULA)
               ELSE
                  NSORB = NSORB + 1
                  IPTSYM(NAORB,LA) = NSORB

                  IRREPN        = IRREP(NAORB)
                  IRREP(NAORB)  = IRREPN + 1

                  IPTYP(NSORB)  = NHKOFF(NHKTA) + NA
                  IPCEN(NSORB)  = NCENT(IA)
                  ICLASS(NSORB) = JC
                  IPIND(NSORB)  = IA*2**16 + NA*2**8 + IRREPN ! IA: AO shell, NA: component of this shell
                  NAOS(LA+1)    = NAOS(LA+1) + 1
                  IF (JC .GE. 0) NCOS(LA+1,JC) = NCOS(LA+1,JC) + 1
C
C                 loop over symmetry dependent centers
C
                  JKB = 0
                  DO KB = 0, MAXREP
                     IF (IBTAND(KB,MULA) .EQ. 0) THEN
                        JKB              = JKB    + 1
                        NAORBD           = NAORBD + 1
                        CTRAN(NSORB,JKB) = PT(IBTAND(KB,IVARB))
                        ITRAN(NSORB,JKB) = NAORBD

                        IAOINFO(NAORBD,1) = IA                     !  sym. indep. AO index
                        IAOINFO(NAORBD,2) = NUCNUM(NCENT(IA),KB+1) !  sym. dep. atom index (e.g. to NAMDEP)
                        IAOINFO(NAORBD,3) = NHKOFF(NHKTA) + NA     !  orbital type

                     END IF
                  END DO
                  JTRAN(NSORB) = JKB ! number of AOs contributing to this SO

                  IF (IC .EQ. ISETHUCKEL) THEN
                     NQ     = IZATOM(NCENT(IA))
C                    ... atomic charge of this center
                   IF (NQ .GT. 0) THEN
C                  ... skip point charges (code NQ=-1234567890),
C                      and centers for multiple basis sets (NQ=-Z)
C                      and floating orbital centers (NQ=0)
                     IHNMAO = IHNMAO + 1
                     IHUCPT(IA) = IHNMAO
                     IQCORE = NINT(CHARGE(NCENT(IA)))
C                    ... effective charge of this center
                     IQCORE = NQ - IQCORE
C                    ... .gt. 0 if ECP otherwise .eq. 0
                     NHUCCO = NUMCF(IA)
                     CALL HUCFUN(NQ,IQCORE,NHKTA,NHUCCO,HUCEXC(IHNMAO))
                   END IF
                  END IF

               END IF
  420       CONTINUE
  410    CONTINUE
         NBASI  = NAOS(LA+1)
         NPARLA = NBASI*(NBASI + 1)/2
         IF (LA .LT. MAXREP) NPARSU(LA+2) = NPARSU(LA+1) + NPARLA
  400 CONTINUE
C
C     Compute symmetry offsets for non-symmetric operators
C     (Note case IREPO = 0 is thereby excluded)
C
      DO 500 IREPO = 1,MAXREP
         IBLK = 0
         DO 500 IREPA = 0,MAXREP
            IREPB = IBTXOR(IREPO,IREPA)
            IF (IREPA .GT. IREPB) THEN
               NPARNU(IREPO+1,IREPA+1) = IBLK
               IBLK = IBLK + NAOS(IREPA+1)*NAOS(IREPB+1)
            ENDIF
 500  CONTINUE
      RETURN
      END
C  /* Deck symout */
      SUBROUTINE SYMOUT
C*****************************************************************************
C
C     Print routine for symmetry
C     Extended to write a SYMINFO file for use with mean-field SO integrals
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
C
#include "nuclei.h"
#include "ccom.h"
#include "shells.h"
#include "symmet.h"
#include "pgroup.h"
#include "pincom.h"
#include "cbirea.h"
#include "aosotr.h"
#include "chrsgn.h"
      CHARACTER COMP(-MXBSETS:MXBSETS)
#include <r12int.h>
#include <huckel.h>
casm
#include <suscc.h>
casm

C
C     Write SYMINFO file
C
      LUSYMF = -1
      IDUMMY = -1
      CALL GPOPEN(LUSYMF,'SYMINFO',' ',' ','FORMATTED',IDUMMY,.FALSE.)
      WRITE (LUSYMF,'(A,/A)')
     &     'Symmetry information from Dalton (Hermit module)',
     &     'funct #, unique centre #, L, M , Func.type #, # of '//
     &     'sym.gen.functions , Phases'
C
C*****************************************************************************
C       S Y M M E T R Y      O R B I T A L S
C*****************************************************************************
      IF (MAXREP.GT.0.AND.IPREAD.GT.0) THEN
         CALL HEADER('Symmetry Orbitals',1)
         WRITE (LUPRI,'(A,6X,8I5)')
     &     '  Number of orbitals in each symmetry: ',
     &     (NAOS(I),I=1,MAXREP+1)
      END IF
      IF(DIRAC) THEN
        COMP(1) = 'L'
        COMP(2) = 'S'
        IF (IPREAD .GT. 0) THEN
           WRITE (LUPRI,'(A,8I5)')
     &        '  Number of large orbitals in each symmetry: ',
     &        (NCOS(I,1),I=1,MAXREP+1)
           WRITE (LUPRI,'(A,8I5)')
     &        '  Number of small orbitals in each symmetry: ',
     &        (NCOS(I,2),I=1,MAXREP+1)
        END IF
      ELSE IF (LMULBS) THEN
C       Output for multiple basis sets. Determine total numbers of orbital-basis 
C       (MBAS1) and auxiliary-basis (MBAS2) functions (WK/UniKA/04-11-2002).
        COMP(1) = 'M'
        COMP(2) = 'A'
        IF (IPREAD.GE.0) THEN
           WRITE (LUPRI,'(/A,8I5)')
     *          '  Number of main orbitals in each symmetry:       ',
     &          (NCOS(I,1),I=1,MAXREP+1)
           WRITE (LUPRI,'(A,8I5)')
     *          '  Number of auxiliary functions in each symmetry: ',
     &          (NCOS(I,2),I=1,MAXREP+1)
        END IF
        MBAS1T = 0
        MBAS2T = 0
        DO I=1,MAXREP+1
          MBAS1(I) = NCOS(I,1)
          MBAS1T   = MBAS1T + MBAS1(I)
          MBAS2(I) = NCOS(I,2)
          MBAS2T   = MBAS2T + MBAS1(I)
        END DO
      ELSE
        COMP(1) = ' '
        COMP(2) = ' '
      ENDIF
C
      COMP(-1) = 'D' ! for density fitting
      COMP( 0) = 'H' ! for Huckel
      IF (ATOMDF .AND. IPREAD.GT.0) WRITE (LUPRI,'(A,8I5)')
     &    '  Number of density fitting functions in each symmetry: ',
     &    (NCOS(I,-1),I=1,MAXREP+1)
Chjaaj Nov 07: Huckel orbitals are not counted here, as KMAX
C in DO 410 loop in SYMPRO does not include Huckel orbitals
C (see IF (DOHUCKEL) ... code in BASPRO)
C     IF (DOHUCKEL .AND. IPREAD.GT.0) WRITE (LUPRI,'(A,8I5)')
C    &    '  Number of Huckel orbitals in each symmetry:     ',
C    &    (NCOS(I,0),I=1,MAXREP+1)
Chjaaj Nov 07 - end.
      I=0
      DO 10 LA=1,MAXREP+1
        NBI=NAOS(LA)
        IF (NBI.GT.0) THEN
           IF (MAXREP.GT.0.AND.IPREAD.GT.0) 
     &          WRITE (LUPRI,'(//2X,A,2X,A3,A1,I2,A1/)')
     &          'Symmetry',REP(LA-1),'(',LA,')'
          IPLMAX = 0
          IPNUM  = 1
          IPLBKP = 0
          ICNTBK = 0
          DO 20 L=1,NBI
            I=I+1
            ICENT = IPCEN(I)
            J=NUCDEG(ICENT)
            IF (MAXREP.GT.0.AND.IPREAD.GT.0)
     &         WRITE (LUPRI,1060) I,NAMN(ICENT),COMP(ICLASS(I)),
     +            GTOTYP(IPTYP(I)), ITRAN(I,1),
     *            (CHRSGN(NINT(CTRAN(I,K))),ITRAN(I,K),K=2,J)
C           Write LUSYMF record for amfi:
            IF (DOCART) THEN
C hjaaj oct 2003: amfi only for spherical GTO, signal this problem
C                 to amfi with negative IPL and IPM
               IPL = -1000
               IPM = -1000
            ELSE IF (IPTYP(I) .EQ. 1) THEN
               IPL = 0
               IPM = 0
            ELSE IF (IPTYP(I) .LE. 4) THEN
               IF (IPTYP(I) .EQ. 2) IPM = 1
               IF (IPTYP(I) .EQ. 3) IPM = -1
               IF (IPTYP(I) .EQ. 4) IPM = 0
               IPL = 1
            ELSE IF (IPTYP(I) .LE. 9) THEN
               IPM = IPTYP(I) - 7
               IPL = 2
            ELSE IF (IPTYP(I) .LE. 16) THEN
               IPM = IPTYP(I) - 13
               IPL = 3
            ELSE IF (IPTYP(I) .LE. 25) THEN
               IPM = IPTYP(I) - 21
               IPL = 4
            ELSE
               IPL = -10000
               IPM = IPL
C hjaaj nov 2001: found out that IPL,IPM was only programmed
C  for s,p,d,f,g orbitals here. In order to avoid wrong results
C  in amfi/symtra.F I have inserted negative IPL as code for
C  undefined type such that amfi can stop (hitherto it would
C  erroneously have used the previous IPL and IPM values!).
C  We don't want to stop here, because if not amfi it doesn't matter.
C
            END IF
            IF (IPL .EQ. IPLBKP .AND. .NOT. ICENT .NE. ICNTBK) THEN
               IF (IPTYP(I) .LE. IPLMAX) THEN
                  IPNUM = IPNUM + 1
               END IF
            ELSE
               IPNUM = 1
            END IF
            IPLMAX = IPTYP(I)
            IPLBKP = IPL
            ICNTBK = ICENT
            WRITE (LUSYMF,1070) I, ICENT, IPL, IPM, IPNUM,
     &           J, (NINT(CTRAN(I,K)),K=1,J)
   20     CONTINUE
        ELSE
           IF (MAXREP.GT.0.AND.IPREAD.GT.0)
     &        WRITE (LUPRI,'(//2X,A,2X,A3,A1,I2,A1)')
     &        'No orbitals in symmetry',REP(LA-1),'(',LA,')'
        END IF
  10  CONTINUE
      IF (NFLOAT .GT. 0) WRITE (LUPRI,'(/A,I5)')
     &   '  Number of floating/multiple basis set orbitals :',NFLOAT
      IF (MAXREP .GT. 0 .AND. IPREAD .GT. 0) THEN
         IF (IPREAD .GE. 10) THEN
            WRITE(LUPRI,'(10X,A/)') 'Symmetry pointer indices'
            WRITE(LUPRI,'(4X,8I5)')
     &           ((IPTSYM(I,J),J = 0,MAXREP),I=1,NORBS)
         END IF
         WRITE (LUPRI,'(/A,3(2X,A3,A1,I1,A1))')
     &        '  Symmetries of electric field:',
     &        (REP(ISYMAX(I,1)),'(',ISYMAX(I,1)+1,')',I=1,3)
         WRITE (LUPRI,'(/A,3(2X,A3,A1,I1,A1))')
     &        '  Symmetries of magnetic field:',
     &        (REP(ISYMAX(I,2)),'(',ISYMAX(I,2)+1,')',I=1,3)
      END IF
 1060 FORMAT(I5,3X,A6,1X,A1,3X,A4,I8,7(1X,A,I5))
 1070 FORMAT(I5,13I6)
      WRITE (LUSYMF,'(A)') 'END'
      CALL GPCLOSE(LUSYMF,'KEEP')
C
C     CC paramagnetic susceptibilty stuff
C
      DO I = 1,3
         KSYMAG(I) = ISYMAX(I,2) + 1
         KXYZ(I) = I
      END DO
C
      DO I = 3,1,-1
         DO J = I-1,1,-1
            IF (KSYMAG(J) .GT. KSYMAG(I)) THEN
               ITMP = KXYZ(I)
               KXYZ(I) = KXYZ(J)
               KXYZ(J) = ITMP
            END IF
         END DO
      END DO
C
      RETURN
      END
C  /* Deck symnco */
      SUBROUTINE SYMNCO(LWORK,WORK)
#include <implicit.h>
#include <priunit.h>
#include <aovec.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxmom.h>
#include <maxorb.h>
C
#include <cbirea.h>
#include <nuclei.h>
#include <symmet.h>
      DIMENSION WORK(LWORK)
#include <chrxyz.h>
#include <ibtfun.h>
      LENGTH = 48*NCTOT
      CALL IZERO(IPTCNT(1,0,1),LENGTH)
      IOFFT = 0
      IOFFR = 0
c     IOFFTN = 0
      DO 100 IREP = 0,MAXREP
         ICENTA = 0
         ICENTB = 0
         DO 110 ICENT = 1,NCTOT
            MULC = ISTBNU(ICENT)
            DO 200 IDIRN = 1,3
               ICENTA = ICENTA + 1
               IF (IBTAND(MULC,IBTXOR(IREP,ISYMAX(IDIRN,1))).EQ.0)THEN
                  IOFFT = IOFFT + 1
                  IPTCNT(ICENTA,IREP,1) = IOFFT
               ENDIF
               IF (IBTAND(MULC,IBTXOR(IREP,ISYMAX(IDIRN,2))).EQ.0)THEN
                  IOFFR = IOFFR + 1
                  IPTCNT(ICENTA,IREP,2) = IOFFR
               ENDIF
  200       CONTINUE
c           DO 300 IDIR1 = 1, 3
c           DO 300 IDIR2 = IDIR1, 3
c              ICENTB = ICENTB + 1
c              ISYMIJ = IBTXOR(ISYMAX(IDIR1,1),ISYMAX(IDIR2,1))
c              IF (IBTAND(MULC,IBTXOR(IREP,ISYMIJ)) .EQ. 0) THEN
c                 IOFFTN = IOFFTN + 1
c              ENDIF
c 300       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      IOFF = 0
      DO 400 IREP = 0, MAXREP
         DO 410 ICENT = 1, NCTOT
            IPTNUC(ICENT,IREP) = 0
  410    CONTINUE
         DO 420 ICENT = 1, NCTOT
            MULC = ISTBNU(ICENT)
            IF (IBTAND(MULC,IREP) .EQ. 0) THEN
               IOFF = IOFF + 1
               IPTNUC(ICENT,IREP) = IOFF
            END IF
  420    CONTINUE
  400 CONTINUE
C
      IF (IPREAD .GT. 10) THEN
         CALL HEADER('Symmetry-adapted nuclear coordinates',-1)
         WRITE(LUPRI,'(6X,A//4X,A,3X,8I5)') '( IPTCNT(:,:,1) array )',
     &   'Irrep:',(IREP, IREP = 0,MAXREP)
         ICENTA = 0
         DO 500 ICENT = 1,NUCIND
            WRITE(LUPRI,'(4X,A6,2X,A1,8I5)')
     *      NAMN(ICENT),'X', (IPTCNT(ICENTA+1,IREP,1),IREP = 0,MAXREP)
            WRITE(LUPRI,'(4X,A6,2X,A1,8I5)')
     *      NAMN(ICENT),'Y', (IPTCNT(ICENTA+2,IREP,1),IREP = 0,MAXREP)
            WRITE(LUPRI,'(4X,A6,2X,A1,8I5)')
     *      NAMN(ICENT),'Z', (IPTCNT(ICENTA+3,IREP,1),IREP = 0,MAXREP)
            ICENTA = ICENTA + 3
  500    CONTINUE
         CALL HEADER('Symmetry-adapted nuclear magnetic moments',-1)
         WRITE(LUPRI,'(6X,A//4X,A,3X,8I5)') '( IPTCNT(:,:,2) array )',
     &   'Irrep:',(IREP, IREP = 0,MAXREP)
         ICENTA = 0
         DO 510 ICENT = 1,NUCIND
            WRITE(LUPRI,'(4X,A6,2X,A1,8I5)')
     *      NAMN(ICENT),'X', (IPTCNT(ICENTA+1,IREP,2),IREP = 0,MAXREP)
            WRITE(LUPRI,'(4X,A6,2X,A1,8I5)')
     *      NAMN(ICENT),'Y', (IPTCNT(ICENTA+2,IREP,2),IREP = 0,MAXREP)
            WRITE(LUPRI,'(4X,A6,2X,A1,8I5)')
     *      NAMN(ICENT),'Z', (IPTCNT(ICENTA+3,IREP,2),IREP = 0,MAXREP)
            ICENTA = ICENTA + 3
  510    CONTINUE
         CALL HEADER('Symmetry-adapted nuclei',-1)
         WRITE(LUPRI,'(6X,A//4X,A,3X,8I5)') '( IPTNUC(:,:) array )',
     &   'Irrep:',(IREP, IREP = 0,MAXREP)
         DO 520 ICENT = 1, NUCIND
            WRITE(LUPRI,'(4X,A6,3X,8I5)')
     &      NAMN(ICENT),(IPTNUC(ICENT,IREP),IREP = 0,MAXREP)
  520    CONTINUE
      END IF
C
C     ***************************
C     **** Initialize IPTCOR ****
C     ***************************
C
      NCOOR  = 3*NUCDEP
      KCSTRA = 1
      KSCTRA = KCSTRA + NCOOR*NCOOR
      KLAST  = KSCTRA + NCOOR*NCOOR
      IF (KLAST.GT.LWORK) CALL STOPIT('NUCPRO','TRACOR',KLAST,LWORK)
      CALL TRACOR(WORK(KCSTRA),WORK(KSCTRA),1,NCOOR,IPREAD)
      CALL TRACOR(WORK(KCSTRA),WORK(KSCTRA),2,NCOOR,IPREAD)
      CALL TRACR(NCOOR)
      END
C  /* Deck pricar */
      SUBROUTINE PRICAR(WORK,LWORK)
C
C     This subroutine prints information about atomic coordinates
C
C     tuh 081188 - Bush elected
Clf       031104 - Bush Jr. re-elected
C
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <gnrinf.h>
C
      CHARACTER*4 NAME
      INTEGER NBASE(8), NSIGN(8)
      LOGICAL LINEAR
      DIMENSION WORK(LWORK)
C
#include <nuclei.h>
#include <pcmnuclei.h>
#include <symmet.h>
#include <pgroup.h>
#include <cbirea.h>
#include <frame.h>
C
#include <orgcom.h>
#include <chrxyz.h>
#include <chrsgn.h>
#include <chrnos.h>
#include <ibtfun.h>
C
      NCOOR = 3*NUCDEP
C
C     *********************************
C     ***** Cartesian Coordinates *****
C     *********************************
C
      IPTEST = MIN(NUCDEP/13,2)
C     hjaaj: x,y,z normally only interesting for
C            small molecules, otherwise just a lot of lines.
      IF (IPREAD .GT. IPTEST) THEN
         CALL HEADER('Cartesian Coordinates (a.u.)',1)
         WRITE (LUPRI,'(A,I5)')'  Total number of coordinates:',NCOOR
      END IF

C
C        MIDAS Interface, write out input coordinates for midas.sym
C        analysis
C
      IF (.TRUE.) THEN 
         LUCOOR = -1
         CALL GPOPEN(LUCOOR,'midasifc.cartrot','NEW', 
     &               ' ','FORMATTED',IDUMMY,.FALSE.)
         CALL PRIGEOLU(LUCOOR,CORD)
         CALL GPCLOSE(LUCOOR,'KEEP')
      ENDIF 
C
C
      IPCMN  = 0
      NCENTS = 0
      KMASS  = NCOOR + 1
C
      ICRX   = 1
      DO 100 ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
         NAME   = NAMEX(3*ICENT)(1:4)
         IF (MULT(MULCNT) .EQ. 1) THEN
            IPCMN = IPCMN + 1
            IF (IPREAD .GT. IPTEST) THEN
               WRITE (LUPRI,'(2X,A,3X," : ",3(I5,2X,A,F15.10))')
     *              NAME, ICRX  , CHRXYZ(-1), CORD(1,ICENT),
     *                    ICRX+1, CHRXYZ(-2), CORD(2,ICENT),
     *                    ICRX+2, CHRXYZ(-3), CORD(3,ICENT)
            END IF
            WORK(KMASS + NCENTS) = DISOTP(IZATOM(ICENT),1,'MASS')
            NCENTS = NCENTS + 1
            WORK(ICRX  ) = CORD(1,ICENT)
            WORK(ICRX+1) = CORD(2,ICENT)
            WORK(ICRX+2) = CORD(3,ICENT)
            ICRX = ICRX + 3
         ELSE
            JATOM = 0
            DO 200 ISYMOP = 0, MAXOPR
               IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
                  IPCMN = IPCMN + 1
                  JATOM  = JATOM + 1
                  WORK(KMASS + NCENTS) = 
     *                 DISOTP(IZATOM(ICENT),1,'MASS')
                  NCENTS = NCENTS + 1
                  WORK(ICRX  ) = PT(IBTAND(ISYMAX(1,1),ISYMOP))
     *                        *CORD(1,ICENT)
                  WORK(ICRX+1) = PT(IBTAND(ISYMAX(2,1),ISYMOP))
     *                        *CORD(2,ICENT)
                  WORK(ICRX+2) = PT(IBTAND(ISYMAX(3,1),ISYMOP))
     *                        *CORD(3,ICENT)
                  IF (IPREAD .GT. IPTEST)
     &            WRITE (LUPRI,'(2X,A,"/",I2," : ",3(I5,2X,A,F15.10))')
     &               NAME,JATOM,ICRX,CHRXYZ(-1),WORK(ICRX),
     &                        ICRX+1,CHRXYZ(-2),WORK(ICRX+1),
     &                        ICRX+2,CHRXYZ(-3),WORK(ICRX+2)
                  ICRX = ICRX + 3
               END IF
  200       CONTINUE
         END IF
  100 CONTINUE
      NPCMN = IPCMN
C
C     ********************************
C     ***** Symmetry Coordinates *****
C     ********************************
C
      ISTRT = NCOOR + NCENTS + 1
      IF (LWORK .LT. 2*NCOOR*NCOOR) CALL STOPIT('PRICAR','TRACOR',
     &                                          2*NCOOR*NCOOR,LWORK)
      CALL TRACOR(WORK(ISTRT),WORK(ISTRT + NCOOR*NCOOR),1,NCOOR,IPREAD)
      IF (MAXREP .GT. 0) THEN
         IF (IPREAD .GT. 0) THEN
            CALL HEADER('Symmetry Coordinates',1)
            WRITE (LUPRI,'(A,8I5)')
     *            '  Number of coordinates in each symmetry: ',
     *            (NCRREP(I,1),I=0,MAXREP)
         END IF
         DO 300 ISYM = 0, MAXREP
         IF (NCRREP(ISYM,1) .GT. 0) THEN
            IF (IPREAD .GT. IPTEST) WRITE (LUPRI,'(/A,2X,A3,A2,I2,A1/)')
     &          '  Symmetry',REP(ISYM),' (',ISYM+1,')'
            DO 400 IATOM = 1, NUCIND
               DO 500 ICOOR = 1, 3
                  ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,ISYM,1)
                  IF (ISCOOR .GT. 0) THEN
                     NB = 0
                     DO 600 I = 1, NCOOR
                        IADR = ISTRT + (I - 1)*NCOOR + ISCOOR - 1
                        NSGN = NINT(WORK(IADR))
                        IF (NSGN .NE. 0) THEN
                           NB = NB + 1
                           NBASE(NB) = I
                           NSIGN(NB) = NSGN
                        END IF
  600                CONTINUE
                     IF (IPREAD .GT. IPTEST) THEN
                     IF (NB .EQ. 1) THEN
                        WRITE (LUPRI,'(I5,3X,A,2X,A,I5)') ISCOOR,
     &                     NAMEX(3*IATOM)(1:4), CHRXYZ(-ICOOR), NBASE(1)
                     ELSE
                        WRITE (LUPRI,'(I5,3X,A,2X,A,3X,A,I3,'
     &                   //CHRNOS(NB-1)//'(2X,A,I5),A,I1)')
     &                   ISCOOR,NAMEX(3*IATOM)(1:4),CHRXYZ(-ICOOR),' [',
     &                   NBASE(1),(CHRSGN(NSIGN(I)),NBASE(I),I=2,NB),
     &                   ' ]/',NB
                     END IF
                     END IF
                  END IF
  500          CONTINUE
  400       CONTINUE
         END IF
  300    CONTINUE
      END IF
C
      IF (NUCDEP .EQ. 1) THEN
         IF (IPREAD .GT. 0)
     &      WRITE (LUPRI,'(//A)') '@ This is an atomic calculation.'
      ELSE
C
C     Print internuclear distances (if not an atom)
C
         LUIP = -100
         IF ((IPREAD .GT. 0) .OR. (GEOALL)) 
     &   CALL GEOANA(CORD,.TRUE.,.FALSE.,NBONDS,LUIP,WORK(ISTRT),LWORK)
C
C        As a new default, we (almost) always print a rotational analysis, kr-02
C
         IF (IPREAD .GT. 0) CALL ROTANA(WORK(1),WORK(KMASS),WORK(ISTRT),
     &            CMXYZ,AIMOM,BIMOM,CIMOM,LINEAR,NCENTS,NCOOR,IPREAD)
C
      IF (IPREAD .GT. 0) WRITE(LUPRI,'(//A,F18.12,A)')
     &   '@  Nuclear repulsion energy :',POTNUC,' Hartree'
      END IF
      RETURN
      END
C  /* Deck priord */
      SUBROUTINE PRIORD(ALPHA,CPRIMU,NPRI,NRCI,SEG,KPRIM,DSM)
C
C     Order primitive basis functions
C
#include "implicit.h"
#include "maxorb.h"
      LOGICAL SEG
      DIMENSION ALPHA(KPRIM), CPRIMU(KPRIM,KPRIM)
C
C     *****************************
C     *** Segmented contraction ***
C     *****************************
C
      IDONE = 0
      IF (SEG) THEN
         DO 100 ICONTR = 1, NRCI
            ISTART = IDONE + 1
            NLEFT  = NPRI - IDONE
C
C           Find first primitive
C           ====================
C
            IMXA = ISTART + IDAMAX(NLEFT,ALPHA(ISTART),1) - 1
            CALL DSWAP(1,ALPHA(ISTART),1,ALPHA(IMXA),1)
            CALL DSWAP(NRCI,CPRIMU(ISTART,1),KPRIM,CPRIMU(IMXA,1),
     &                 KPRIM)
C
C           Find corresponding contracted function
C           ======================================
C
            IMXC = IDAMAX(NRCI,CPRIMU(ISTART,1),KPRIM)
            CALL DSWAP(NPRI,CPRIMU(1,ICONTR),1,CPRIMU(1,IMXC),1)
C
C           Collect other primitives contributing to same contracted
C           ========================================================
C
            IPRI = 1
            DO 200 I = ISTART + 1, NPRI
               IF (ABS(CPRIMU(I,ICONTR)) .GT. DSM) THEN
                  CALL DSWAP(1,ALPHA(I),1,ALPHA(ISTART+IPRI),1)
                  CALL DSWAP(NRCI,CPRIMU(I,1),KPRIM,
     &                            CPRIMU(ISTART+IPRI,1),KPRIM)
                  IPRI = IPRI + 1
               END IF
  200       CONTINUE
C
C           Sort primitives
C           ===============
C
            IF (IPRI .GT. 2) THEN
               DO 300 I = ISTART + 1, ISTART + IPRI - 2
                  DO 400 J = I + 1, ISTART + IPRI - 1
                     IF (ALPHA(J) .GT. ALPHA(I)) THEN
                        CALL DSWAP(1,ALPHA(I),1,ALPHA(J),1)
                        CALL DSWAP(NRCI,CPRIMU(I,1),KPRIM,
     &                                  CPRIMU(J,1),KPRIM)
                     END IF
  400             CONTINUE
  300          CONTINUE
            END IF
C
            IDONE = IDONE + IPRI
  100    CONTINUE
C
C     ***************************
C     *** General contraction ***
C     ***************************
C
      ELSE
         DO 500 I = 1, NPRI - 1
            DO 600 J = I + 1, NPRI
            IF(ALPHA(J) .GT. ALPHA(I)) THEN
               CALL DSWAP(1,ALPHA(I),1,ALPHA(J),1)
               CALL DSWAP(NRCI,CPRIMU(I,1),KPRIM,CPRIMU(J,1),KPRIM)
           END IF
  600      CONTINUE
  500   CONTINUE
      END IF
      RETURN
      END
C  /* Deck nucsiz */
      SUBROUTINE NUCSIZ(INUC,DELTA)
C.......................................................................
C
C     For a nucleus of atomic no. "charge", calculate the exponent delta
C     for a gaussian charge distribution approximation to the finite nuc
C.......................................................................
#include "implicit.h"
#include "priunit.h"
C     include character ASYMB(0:NSYNB)*2 with atomic symbols
#include <asymb.h>
      DIMENSION AMASS(103)
      DATA (AMASS(I),I = 1,92)
     1/ 1.0078246D0, 4.002601D0, 7.01600D0, 9.01218D0, 11.009307D0,
     2  12.000000D0, 14.0030738D0, 15.9949141D0,
     *  18.9984022D0, 19.992441D0,
     3  22.9898D0, 23.98504D0, 26.98153D0, 27.976929D0, 30.973764D0,
     4  31.9720727D0, 34.9688531D0, 39.962386D0, 38.96371D0, 39.96259D0,
     5  44.95592D0, 48.D0, 50.9440D0, 51.9405D0, 54.9380D0,
     6  55.9349D0, 58.9332D0, 57.9353D0, 62.9296D0, 63.9291D0,
     7  68.9257D0, 73.9219D0, 74.9216D0, 79.9165D0, 78.91839D0,
     8  83.91151D0, 84.9117D0, 87.9056D0, 88.9059D0, 89.9043D0,
     9  92.9060D0, 97.9055D0, 98.D0, 101.9037D0, 102.9048D0,
     O  107.90389D0, 106.90509D0, 113.9036D0, 114.9041D0, 120.D0,
     1  120.9038D0, 129.9067D0, 126.90466D0, 131.90416D0, 132.9051D0,
     2  137.9050D0, 138.9061D0, 139.9053D0, 140.9074D0, 141.9075D0,
     3  145.D0, 151.9195D0, 152.9209D0, 157.9241D0, 159.9250D0,
     4  163.9288D0, 164.9303D0, 165.9304D0, 168.9344D0, 173.9390D0,
     5  174.9409D0, 179.9468D0, 180.9480D0, 183.9510D0, 186.9560D0,
     6  192.D0, 192.9633D0, 194.9648D0, 196.9666D0, 201.970625D0,
     7  204.9745D0, 207.9766D0, 208.9804D0, 209.D0, 210.D0,
     8  222.D0, 223.D0, 226.D0, 227.D0, 232.D0, 231.D0, 238.D0 /
      DATA (AMASS(I),I = 93,103)
     1/ 237.D0, 244.D0, 243.D0,
     2  247.D0, 247.D0, 251.D0, 252.D0, 257.D0,
     3  258.D0, 259.D0, 260.D0/
      IF (INUC.LE.103) THEN
         A = AMASS(INUC)
      ELSE
         A = INUC*2.556D0
      END IF
      A = A**(1.0D0/3.0D0)
      DELTA = 1.5D10*(.52917706D0/(.836D0*A+0.57D0))**2
      RETURN
      END
C  /* Deck sphlab */
      SUBROUTINE SPHLAB(IORDER,LABINT)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      CHARACTER LABINT(MXQN**2)*4
      CHARACTER SPDCAR*1, SIGNJ*1
#include "chrnos.h"
C
      IOFF = 0
      DO 100 I = 0, IORDER
         IF (I .EQ. 0) THEN
            LABINT(1) = '1s  '
            IOFF = IOFF + 1
         ELSE IF (I .EQ. 1) THEN
            LABINT(2) = '2px '
            LABINT(3) = '2py '
            LABINT(4) = '2pz '
            IOFF = IOFF + 3
         ELSE
            DO 200 J = -I, I
               IADR = IOFF + J + I + 1
               IF (J .LT. 0) THEN
                  SIGNJ = '-'
               ELSE IF (J .EQ. 0) THEN
                  SIGNJ = ' '
               ELSE
                  SIGNJ = '+'
               END IF
               LABINT(IADR) = CHRNOS(I+1)//SPDCAR(I)//
     &                        CHRNOS(ABS(J))//SIGNJ
  200       CONTINUE
            IOFF = IOFF + 2*I + 1
         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck carlab */
      SUBROUTINE CARLAB(IORDER,LABINT)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      DIMENSION IX(MXAQN), IY(MXAQN), IZ(MXAQN)
      CHARACTER LABINT(MXQN*(MXQN+1)*(MXQN+2)/6)*4
      CHARACTER SPDCAR*1
#include "chrnos.h"
      J = 0
      DO 100 I = 0, IORDER
         NFUN = (I + 1)*(I + 2)/2
         IF (I .LE. 3) THEN
            J = J + NFUN
            CALL SETCLB(I,LABINT,MXQN*(MXQN+1)*(MXQN+2)/6)
         ELSE
            CALL LMNVAL(I+1,NFUN,IX,IY,IZ)
            DO 200 K = 1, NFUN
               J  = J + 1
               NX = IX(K)
               NY = IY(K)
               NZ = IZ(K)
               LABINT(J) = SPDCAR(I)//CHRNOS(NX)//CHRNOS(NY)//CHRNOS(NZ)
  200       CONTINUE
         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck aotoao */
      SUBROUTINE AOTOAO(IAOAO,JAOAO,IPREAD)
C
C     This subroutine sets up pointer IAOAO which converts between
C     two different orderings of AO's: from the ordering of
C     MOLECULE to an ordering in which the outer loop is over atoms.
C     (These orderings are identical when no symmetry is used.)
C
C     The purpose of this ordering is to make the AO's appear in the
C     same order as in the corresponding run with no symmetry.
C
C     tuh 120988
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
C
      DIMENSION IAOAO(MXCORB),JAOAO(MXCORB,2)
C
#include "symmet.h"
#include "nuclei.h"
#include "shells.h"
C
#include "ibtfun.h"
C
      IF (IPREAD .GE. 10) CALL TITLER('Output from AOTOAO','*',103)
      DO I = 1,MXCORB
         JAOAO(I,1) = 0
         JAOAO(I,2) = 0
      END DO
C
C     JATOM: Loop over all atoms
C     IAO  : Loop over AO's in atom order (C1 order)
C
      IAO   = 0
      JATOM = 0
      DO 100 IATOM = 1, NUCIND
         ISTABA = ISTBNU(IATOM)
         DO 200 ISYMOP = 0, MAXREP
         IF (IBTAND(ISYMOP,ISTABA) .EQ. 0) THEN
            JATOM = JATOM + 1
C
C           Loop over all orbitals
C           JAO  : AO's in order when symmetry
C
            JAO = 0
            DO 300 ISHELL = 1, KMAX
               ICENT  = NCENT(ISHELL)
               ISTABO = ISTBAO(ISHELL)
               DO 400 ICMP = 1, KHKT(ISHELL)
Clf
Clf Molden needs d, f,.. functions ordered in the following way
Clf 0,+1,-1,+2,-2,+3,-3 etc...
Clf therefore we provide the following index to reorder
Clf the coefficients in the molden.inp file.
Clf The shift index is ISHIFT and it has to be placed
Clf in KAO=IAO-ISHIFT
Clf WARNING: this only works for spherical functions 
Clf and not for cartesians. 
!
!   Example - f-orbital, LMAX=3
!      index         :  1,  2,  3,  4,  5,  6,  7
!      standard order: -3, -2, -1,  0, +1, +2, +3  is LZ + LMAX + 1
!      molden   order:  0, +1, -1, +2, -2, +3, -3  is 2*LZ, for LZ.gt.0; -2*LZ + 1 for LZ.le.0
!      LZ:       -3 -2 -1  0 +1 +2 +3
!      ISHIFT:   -6 -3  0 +3 +3 +2 +1
!
Clf
                  IF(KHKT(ISHELL) .GE. 5) THEN
                     LMAX = (KHKT(ISHELL)-1)/2
                     LZ   = ICMP - 1 - LMAX
                     IF (LZ .EQ. 0) THEN
                        ISHIFT = LMAX
                     ELSE IF (LZ .GT. 0) THEN
                        ISHIFT = LMAX - LZ + 1
                     ELSE
                        ISHIFT = LMAX + 3*LZ
                     END IF
                  ELSE
                     ISHIFT = 0
                  END IF

                  DO 500 JSYMOP = 0, MAXREP
                  IF (IBTAND(JSYMOP,ISTABO) .EQ. 0) THEN
                     JAO = JAO + 1
                     JCENT = NUCNUM(ICENT,JSYMOP+1)
                     IF (JATOM .EQ. JCENT) THEN
                        IAO = IAO + 1
                        IAOAO(JAO) = IAO
                        JAOAO(IAO,1) = JAO

                        KAO = IAO - ISHIFT ! for molden
                        JAOAO(KAO,2) = JAOAO(KAO,1) + ISHIFT ! for molden
                     END IF
                  END IF
  500             CONTINUE
  400          CONTINUE
  300       CONTINUE
C
C           End loop over orbitals
C
         END IF
  200    CONTINUE
  100 CONTINUE
C
C     End loop over atoms
C
      IF (IPREAD .GE. 10) THEN
         CALL HEADER('i -> IAOAO(i), JAOAO(i,1:2)',6)
         DO I = 1, IAO
            WRITE (LUPRI,'(I9,3I6)') I, IAOAO(I), JAOAO(I,1), JAOAO(I,2)
         END DO
      END IF
      RETURN
      END
C/* Deck cntao */
      SUBROUTINE CNTAO(IPRINT)
C 
C     tuh
C     ICNTAO is used in eri.
C     ICNTAO contains index for sym. indep. center of this symmetry
C     orbital (SO, not AO!!!)
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
C
#include "symmet.h"
#include "nuclei.h"
#include "shells.h"
C
#include "ibtfun.h"
C
      IORB = 0
      DO 100 LA = 0, MAXREP
      DO 200 IA = 1, KMAX
         DO 300 NA = 1, KHKT(IA)
            IVARB = IBTXOR(LA,ISYMAO(NHKT(IA),NA))
            IF (IBTAND(ISTBAO(IA),IVARB) .EQ. 0) THEN
               IORB = IORB + 1
               ICNTAO(IORB)  = NCENT(IA)
            END IF
  300    CONTINUE
  200 CONTINUE
  100 CONTINUE
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('ICNTAO',-1)
         DO 400 I = 1, IORB
            WRITE (LUPRI,'(I9,I6)') I, ICNTAO(I)
  400    CONTINUE
      END IF
      RETURN
      END
C  /* Deck magcor */
      SUBROUTINE MAGCOR(IPREAD)
C
C     Printing of symmetrized nuclear magnetic moments
C
C     tuh March 1991
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "chrxyz.h"
#include "chrsgn.h"
C
      INTEGER NBASE(8), NSIGN(8), NRTREP(0:7)
C
#include "cbiher.h"
#include "nuclei.h"
#include "symmet.h"
#include "pgroup.h"
C
#include "ibtfun.h"
C
      IF (NMRISS .AND. MAXREP .GT. 0) THEN
         IF (IPREAD .GT. 0) THEN
            WRITE (LUPRI,'(/)')
            CALL HEADER('Symmetrized nuclear magnetic moments',1)
         END IF
         DO 100 IREP = 0, MAXREP
            NRTREP(IREP) = 0
            DO 110 ICENT = 1, NUCIND
               DO 120 ICOOR = 1, 3
                  ISCOOR = IPTCNT(3*(ICENT - 1) + ICOOR,IREP,2)
                  IF (ISCOOR .GT. 0) NRTREP(IREP) = NRTREP(IREP) + 1
 120           CONTINUE
 110        CONTINUE
 100     CONTINUE
         IF (IPREAD .GT. 0) THEN
            WRITE (LUPRI,'(A,8I3)')
     &           '  Number of components in each symmetry: ',
     &           (NRTREP(I),I=0,MAXREP)
         END IF
C
         DO 200 IREP = 0, MAXREP
         IF (NRTREP(IREP) .GT. 0) THEN
            IF (IPREAD .GT. 0) THEN
              WRITE (LUPRI,'(//2X,A,2X,A3,A1,I2,A1/)')
     &           'Symmetry',REP(IREP),'(',IREP,')'
            END IF
            DO 210 ICENT = 1, NUCIND
               DO 220 ICOOR = 1, 3
                  ISCOOR = IPTCNT(3*(ICENT - 1) + ICOOR,IREP,2)
                  IF (ISCOOR .GT. 0) THEN
                     IVAR = IBTXOR(IREP,ISYMAX(ICOOR,2))
                     NB = 0
                     DO 230 ISYMOP = 0, MAXOPR
                     IF (IBTAND(ISYMOP,ISTBNU(ICENT)) .EQ. 0) THEN
                        NB = NB + 1
                        NSIGN(NB) = NINT(PT(IBTAND(ISYMOP,IVAR)))
                        NBASE(NB) = 3*(NUCNUM(ICENT,ISYMOP+1)-1) + ICOOR
                     END IF
 230                 CONTINUE
                     IF (IPREAD .GT. 0) THEN
                        WRITE (LUPRI,
     &                  '(2X,A,I3,3X,A,2X,A,3X,I2,7(2X,A,I3))')
     &                  'I',ISCOOR, NAMEX(3*ICENT)(1:4), CHRXYZ(-ICOOR),
     &                  NBASE(1), (CHRSGN(NSIGN(I)),NBASE(I),I=2,NB)
                     END IF
                  END IF
  220          CONTINUE
  210       CONTINUE
         END IF
  200    CONTINUE
      END IF
      RETURN
      END
C  /* Deck tracor */
      SUBROUTINE TRACOR(CSTRA,SCTRA,ITYPE,NCOOR,IPRINT)
C
C     Sets up transformation matrices between Cartesian and symmetry
C     nuclear coordinates
C
C     tuh Jun 29 1988
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <thrzer.h>
      PARAMETER (D0 = 0.0D0)
C
#include <nuclei.h>
#include <symmet.h>
C
      DIMENSION CSTRA(NCOOR,NCOOR), SCTRA(NCOOR,NCOOR)
C
#include <ibtfun.h>
C
      CALL DZERO(CSTRA,NCOOR*NCOOR)
      CALL DZERO(SCTRA,NCOOR*NCOOR)
C
C     *****************************************************************
C     **** Transformation matrix Cartesian to symmetry coordinates ****
C     *****************************************************************
C
      DO 100 IREP = 0, MAXREP
         NCRREP(IREP,ITYPE) = 0
         DO 110 ICENT = 1, NUCIND
            MULCNT = ISTBNU(ICENT)
            DO 120 ICOOR = 1, 3
               IF (ITYPE .EQ. 1) THEN
                  ISCOOR = IPTCNT(3*(ICENT-1)+ICOOR,IREP,1)
                  IF (ISCOOR .GT. 0) THEN
                     NCRREP(IREP,1) = NCRREP(IREP,1) + 1
                     IVAR = IBTXOR(IREP,ISYMAX(ICOOR,1))
                     DO 130 ISYMOP = 0, MAXOPR
                        IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
                           ICCOOR = 3*(NUCNUM(ICENT,ISYMOP+1) - 1)+ICOOR
                           CSTRA(ISCOOR,ICCOOR)=PT(IBTAND(ISYMOP,IVAR))
                        END IF
 130                 CONTINUE
                  END IF
               ELSE
                  LSCOOR = IPTCNT(3*(ICENT-1)+ICOOR,IREP,2)
                  IF (LSCOOR .GT. 0) THEN
                     NCRREP(IREP,2) = NCRREP(IREP,2) + 1
                     IVAR = IBTXOR(IREP,ISYMAX(ICOOR,2))
                     DO 140 ISYMOP = 0, MAXOPR
                        IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
                           ICCOOR = 3*(NUCNUM(ICENT,ISYMOP+1) - 1)+ICOOR
                           CSTRA(LSCOOR,ICCOOR)=PT(IBTAND(ISYMOP,IVAR))
                        END IF
 140                 CONTINUE
                  END IF
               END IF
 120        CONTINUE
 110     CONTINUE
 100  CONTINUE
      IF (IPRINT .GE. 20) THEN
         IF (ITYPE .EQ. 1) THEN
            WRITE (LUPRI,'(/A,8I3,/)')
     *           '  Number of coordinates in each symmetry: ',
     *           (NCRREP(I,1),I=0,MAXREP)
            CALL HEADER('Transformation from Cartesian to symmetry '
     *                //'coordinates',-1)
         ELSE
            CALL HEADER('Transformation from rotational cartesian '//
     &                  'to symmetry coordinates',-1)
         END IF
         DO ISYM = 1, NCOOR
            WRITE (LUPRI,'(10F6.2)') (CSTRA(ISYM,J),J=1,NCOOR)
         END DO
      END IF
C
C     Inverted matrix
C
      DO 210 I = 1, NCOOR
         DO 220 J = 1, NCOOR
           SCTRA(I,J) = CSTRA(J,I)
  220    CONTINUE
  210 CONTINUE
      DO 230 ICOL = 1, NCOOR
         SUM1 = D0
         SUM2 = D0
         DO 240 IROW = 1, NCOOR
            SUM1 = SUM1 + ABS(SCTRA(IROW,ICOL))
  240    CONTINUE
         IF (SUM1 .LT. THRZER) GOTO 230
         DO 250 IROW = 1, NCOOR
            SCTRA(IROW,ICOL) = SCTRA(IROW,ICOL)/SUM1
  250    CONTINUE
  230 CONTINUE
C
      IF (IPRINT .GE. 20) THEN
         CALL HEADER(
     &      'The inverted matrix: symmetry to Cartesian',-1)
         DO ISYM = 1, NCOOR
            WRITE (LUPRI,'(10F6.2)') (CSTRA(ISYM,J),J=1,NCOOR)
         END DO
      END IF
C
      RETURN
      END
C  /* Deck tracr */
      SUBROUTINE TRACR(NCOOR)
C
C     Set up IPTCOR - points from symmetry coordinate to generating
C     Cartesian coordinate. Extracted from old TRACOR routine by
C     K.Ruud, Dec-96
C
#include <implicit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
C
#include <symmet.h>
C
      DO 400 ITYPE = 1, 2
         IJ = 0
         DO 410 J = 0, MAXREP
            DO 420 I = 1, NCOOR
               IF (IPTCNT(I,J,ITYPE) .NE. 0) THEN
                  IJ = IJ + 1
                  IPTCOR(IJ,ITYPE) = I
               END IF
  420       CONTINUE
  410    CONTINUE
  400 CONTINUE
      RETURN
      END
C  /* Deck setclb */
      SUBROUTINE SETCLB(I,LABINT,NDIM)
#include "implicit.h"
      CHARACTER LABINT(NDIM)*4
      IF (I .EQ. 0) THEN
         LABINT(1) = 's   '
      ELSE IF (I .EQ. 1) THEN
         LABINT(2) = 'px  '
         LABINT(3) = 'py  '
         LABINT(4) = 'pz  '
      ELSE IF (I .EQ. 2) THEN
         LABINT(5)  = 'dxx '
         LABINT(6)  = 'dxy '
         LABINT(7)  = 'dxz '
         LABINT(8)  = 'dyy '
         LABINT(9)  = 'dyz '
         LABINT(10) = 'dzz '
      ELSE IF (I .EQ. 3) THEN
         LABINT(11)  = 'fxxx'
         LABINT(12)  = 'fxxy'
         LABINT(13)  = 'fxxz'
         LABINT(14)  = 'fxyy'
         LABINT(15)  = 'fxyz'
         LABINT(16)  = 'fxzz'
         LABINT(17)  = 'fyyy'
         LABINT(18)  = 'fyyz'
         LABINT(19)  = 'fyzz'
         LABINT(20)  = 'fzzz'
      END IF
      RETURN
      END
C  /* Deck spdcar */
      FUNCTION SPDCAR(IORDER)
#include "implicit.h"
      CHARACTER*1 SPDCAR
      IF (IORDER .EQ. 0) THEN
         SPDCAR = 's'
      ELSE IF (IORDER .EQ. 1) THEN
         SPDCAR = 'p'
      ELSE IF (IORDER .EQ. 2) THEN
         SPDCAR = 'd'
      ELSE IF (IORDER .LE. 6) THEN
         SPDCAR = CHAR(ICHAR('f') + IORDER - 3)
      ELSE
C
C        Note: j is not used - therefore this special case
C
         SPDCAR = CHAR(ICHAR('k') + IORDER - 7)
      END IF
      RETURN
      END
C  /* Deck sphinp */
      SUBROUTINE SPHINP(LUINFO,WORK,LWORK,DOOWN,MXSPD)
C=======================================================================
C     Calculate spherical harmonic transformation coefficients
C     or, if DOOWN true, read user specified coefficients.
C
C     Modifications to DIRAC by Trond Saue Aug 18 2006
C=======================================================================
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D1 = 1.0D0)
#include "mxcent.h"
#include "maxaqn.h"
#include "sphtrm.h"
#include "molinp.h"
#include "ccom.h"
#include "cbirea.h"
      CHARACTER*1 KEY
      LOGICAL DOOWN
      DIMENSION WORK(LWORK)
C
      IF (DOOWN) THEN
         WRITE(LUPRI,'(/A)') 
     &   '  Reading user specified transformation matrices for GTOs.'
      ELSE IF (IPREAD .GT. 4) THEN
         WRITE(LUPRI,'(/A)') 
     &   '  Calculation of transformation matrices for spherical GTOs.'
      END IF
C
      MXSPD = MXQN
C
      IOFF = 1
      DO I = 1, MXQN
        KCK(I)    = I*(I+1)/2
        ISPADR(I) = IOFF
        IOFF      = IOFF + KCK(I)**2
        IF (DOOWN) THEN
          SPH(I) = .TRUE.
        ELSE
          IF (DOCART) THEN
            KHK(I) = I*(I+1)/2
            SPH(I) = .FALSE.
          ELSE
            KHK(I) = 2*I - 1
            SPH(I) = I .GT. 2
          END IF
        END IF
      END DO
C
      SPHNRM = .TRUE.
C
      IF (DOCART) THEN
         CALL CARLAB(MXQN-1,GTOTYP)
      ELSE
         CALL DZERO(CSP,NCSP)
         IF (.NOT.DOOWN) THEN
            CALL SPHLAB(MXQN-1,GTOTYP)
            KCSPT = 1
            KWRK  = KCSPT + (MXQN*(MXQN + 1)/2)*(2*MXQN - 1)
            IF (KWRK .GT. LWORK) CALL STOPIT('SPHINP',' ',KWRK,LWORK)
            LWRK  = LWORK - KWRK + 1
            DO 200 LVAL = 0, MXQN - 1
               NLM = 2*LVAL + 1
               NXYZ = (LVAL+1)*(LVAL+2)/2
               IF (LVAL .EQ. 0) THEN
                  CSP(1)  = D1
               ELSE IF (LVAL .EQ. 1) THEN
                  CSP (2) = D1
                  CSP (6) = D1
                  CSP(10) = D1
               ELSE
                  IF (.NOT.SPHNRM) THEN
                     MINTEG = 0
                  ELSE
C                    generate normalized s,p,d,f,g,h,i GTO's
                     MINTEG = 2
                  END IF
                  CALL SPHCOM(LVAL,WORK(KCSPT),NLM,NXYZ,1,MINTEG,
     &                        WORK(KWRK),LWRK,IPREAD)
                  CALL MTRSP(NXYZ,NLM,WORK(KCSPT),NXYZ,
     &                       CSP(ISPADR(LVAL+1)),NLM)
               END IF
  200       CONTINUE
         ELSE
            CALL RDLINE(LUINFO)
            READ (MLINE(NMLINE),'(BN,I5)') MXSPD
            IF (MXSPD .GT. MXQN) THEN
              WRITE (LUPRI,'(/A/I10,A,I4//A//A/A)')
     &        ' ERROR: specified MXSPD for GTO transformation matrix',
     &        MXSPD,' is greater than allowed max value of',MXQN,
     &        ' Dump of last two lines read from .mol file:',
     &        MLINE(NMLINE-1),MLINE(NMLINE)
              CALL QUIT('ERROR in reading .mol file, see output')
            END IF
            IOFF = 1
            DO 300 I = 1, MXSPD
               CALL RDLINE(LUINFO)
               READ (MLINE(NMLINE),'(BN,4X,A1,I5)')  KEY, KHK(I)
               IF (KHK(I) .GT. KCK(I)) THEN
                 WRITE (LUPRI,'(/A,I5,A,I3/A,I5,A//A//A)')
     &           ' ERROR: you specified',KHK(I),
     &           ' components for GTO transformation matrix for L =',
     &           I-1,' and that is greater than',KCK(I),
     &           ', the number of Cartesian components.',
     &           ' Dump of last line read from .mol file:',MLINE(NMLINE)
                 CALL QUIT('ERROR in reading .mol file, see output')
               END IF
               CALL RDCSP(LUINFO,CSP(ISPADR(I)),GTOTYP(IOFF),
     &           KHK(I),KCK(I))
               IOFF = IOFF + KHK(I)
  300       CONTINUE
         END IF
         IF (DOOWN .OR. IPREAD .GT. 4) THEN
            CALL HEADER('Cartesian transformation matrices',-1)
            IF (DOOWN) THEN
               WRITE (LUPRI,'(A/)') '  to user specified combinations'
            ELSE
               WRITE (LUPRI,'(A/)') '  to spherical harmonics'
            END IF
            IOFF = 0
            DO 400 I = 1, MXSPD
               KHKI = KHK(I)
               KCKI = KCK(I)
               WRITE (LUPRI,'(/A,I5/A/,(5X,15(A4,1X)))')
     &            '  Coefficients for angular quantum number ',I-1,
     &            '  to GTOs with labels:',(GTOTYP(IOFF+J),J=1,KHKI)
               IOFF = IOFF + KHKI
               CALL OUTPUT(CSP(ISPADR(I)),1,KHKI,1,KCKI,KHKI,KCKI,1,
     &                     LUPRI)
 400        CONTINUE
            CALL FLSHFO(LUPRI)
         END IF
      END IF
      RETURN
      END
C  /* Deck rdcsp */
      SUBROUTINE RDCSP(LUINFO,CSP,GTOTYP,KHKI,KCKI)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "molinp.h"
      CHARACTER*4 GTOTYP(*)
      DIMENSION CSP(KHKI,KCKI)
      DO 100 I = 1, KHKI
         CALL RDLINE(LUINFO)
         READ (MLINE(NMLINE),'(BN,1X,A,5F15.0:/,(5X,5F15.0))')
     &          GTOTYP(I), (CSP(I,J),J=1,KCKI)
  100 CONTINUE
      RETURN
      END
C  /* Deck rdline */
      SUBROUTINE RDLINE(LUINFO)
C
C     Purpose: Read a line from MOLECULE.INP
C
C     We keep track of which line we've read last time (NMLN_LASTREAD).
C     This allows some routines to read a line, check if it want to use
C     that line, and then discard it by decrementing NMLINE.
C     Note that we only allow decrements of one /jth-20000621.
C     (This is used to check if SMALL basis keyword has been specified,
C     or if default - kinetic balance - should be used.)
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "molinp.h"
#include "infpar.h"

      CHARACTER*(120) NEXT_INPUT_LINE
      INTEGER NMLN_LASTREAD
      DATA NMLN_LASTREAD /0/
      SAVE NMLN_LASTREAD, LEN_MLINE
C
C     If NMLINE is zero, reset NMLN_LASTREAD since this is a new call of readin,
C     e.g. during a geometry optimization.
C
      IF ( NMLINE .EQ. 0) THEN
         NMLN_LASTREAD = 0
         LEN_MLINE = LEN( MLINE(1) )
      END IF
C
   10 CONTINUE
      NMLINE = NMLINE + 1

      IF ( NMLINE .LT. NMLN_LASTREAD ) THEN
         WRITE(LUPRI,'(//A,2(/A,I5))')
     &        '*** ERROR ***: NMLINE < NMLN_LASTREAD in RDLINE',
     &        'NMLINE = ',NMLINE,
     &        'NMLN_LASTREAD = ',NMLN_LASTREAD
         CALL QUIT('*** ERROR reading .mol file')
      ENDIF
C
C     Only the master in a parallel calc. reads the input file.
C     The slaves get the input transferred from the master in
C     the common block MOLINP. /hjaaj
C
      IF (.NOT. SLAVE .AND. NMLINE .GT. NMLN_LASTREAD) THEN
          READ (LUINFO,'(A)',END=12,ERR=11) NEXT_INPUT_LINE
          MLINE(NMLINE) = NEXT_INPUT_LINE
      END IF
      IF (NMLINE.GT.KMLINE) THEN
         WRITE (LUPRI,'(//A,/A,I5)')
     &      ' READIN ERROR: Too many lines in input -,',
     &      ' Increase dimension in /MOLINP/ to at least',NMLINE
         CALL QUIT('*** ERROR *** Too many input lines in .mol file')
      END IF

#ifdef VAR_DEBUG
      write (lupri,*) 'debug output from RDLINE, nmline=',nmline
      write (lupri,*) 'on file>',NEXT_INPUT_LINE
      write (lupri,*) 'RDLINE >',MLINE(NMLINE)
#endif

      IF((MLINE(NMLINE)(1:1) .EQ. '#') .OR.
     &   (MLINE(NMLINE)(1:1) .EQ. '!') .OR.
     &   (MLINE(NMLINE)(1:1) .EQ. '$')) GOTO 10

      IF (LNBLNK(NEXT_INPUT_LINE) .GT. LEN_MLINE) THEN
         WRITE(LUPRI,'(//A,I5/2A/2A)')
     &   ' WARNING -- line in ".mol" file truncated, line',NMLINE,
     &   ' Line from file: ',NEXT_INPUT_LINE,
     &   ' Truncated line: ',MLINE(NMLINE)
      END IF

      NMLN_LASTREAD = NMLINE

      RETURN

   11 MLINE(NMLINE)=' '
      WRITE(LUPRI,*) 'ERROR: I/O error in reading of .mol file; '
     &     //'NODE=',MYNUM
      IF (NMLINE>1) THEN
        WRITE(LUPRI,*) 'Line=', NMLINE, 'Prev line:',
     &        MLINE(NMLINE-1)
      ELSE
        WRITE(LUPRI,*) 'Line=', NMLINE
      END IF
      
   12 CONTINUE
      RETURN
      END
C  /* Deck wronel */
      SUBROUTINE WRONEL(TTITLE,NONTYP,NONT,IQM,IFXYZ,KATOM,JCO2,KANG)
C
C     Write interface records
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "dummy.h"
#include "maxorb.h"
#include "maxaqn.h"
C
      CHARACTER*72 TTITLE(2)
      CHARACTER*8 CDATE, CTIME
      CHARACTER*8 TABLE1(4), TABLE2(4), TABLE3(4),
     &            TABLE4(4), TABLE5(4), TABLE6(4)
      DIMENSION QPOL(6), QQ(3), IFXYZ(3), JFXYZ(3)
      DIMENSION CHRN(MXCENT), COOO(MXCENT,3),
     &          INAMN(MXCORB), IGTO(MXCORB)
      DIMENSION NONT(KATOM), IQM(KATOM), JCO2(KANG,KATOM,2)
#include <ccom.h>
#include <shells.h>
#include <symmet.h>
#include <aosotr.h>
#include <frame.h>
#include <inftap.h>
#include <nuclei.h>
#include <sphtrm.h>
      DATA TABLE1 /'********','        ','        ','ISORDK  '/,
     &     TABLE2 /'********','        ','        ','SCFINP  '/,
     &     TABLE3 /'********','        ','        ','SYMINPUT'/,
     &     TABLE4 /'********','        ','        ','TESTDATA'/,
     &     TABLE5 /'********','        ','        ','SPHERICA'/,
     &     TABLE6 /'********','        ','        ','INTEGRAL'/
      DATA QPOL/6*DUMMY/, QQ/3*DUMMY/, JFXYZ/3*IDUMMY/
#include <ibtfun.h>
C
      CALL GETDAT(CDATE,CTIME)
      TABLE1(2) = CDATE
      TABLE2(2) = CDATE
      TABLE3(2) = CDATE
      TABLE4(2) = CDATE
      TABLE5(2) = CDATE
      TABLE6(2) = CDATE
      TABLE1(3) = CTIME
      TABLE2(3) = CTIME
      TABLE3(3) = CTIME
      TABLE4(3) = CTIME
      TABLE5(3) = CTIME
      TABLE6(3) = CTIME
C
      ICENT = 0
      DO 100 N = 1, NUCIND
      DO 100 IREP = 0,MAXREP
         IF (IBTAND(ISTBNU(N),IREP) .EQ. 0) THEN
            ICENT = ICENT + 1
            CHRN(ICENT)   = CHARGE(N)
            COOO(ICENT,1) = PT(IBTAND(ISYMAX(1,1),IREP))*CORD(1,N)
            COOO(ICENT,2) = PT(IBTAND(ISYMAX(2,1),IREP))*CORD(2,N)
            COOO(ICENT,3) = PT(IBTAND(ISYMAX(3,1),IREP))*CORD(3,N)
      END IF
  100 CONTINUE
C
      I = 0
      DO ISYMOP = 0, MAXREP
         NBI = NAOS(ISYMOP+1)
         DO L = 1,NBI
            I = I + 1
            READ(  NAMN(IPCEN(I)),'(A4)',IOSTAT=IOS) INAMN(I)
            READ(GTOTYP(IPTYP(I)),'(A4)',IOSTAT=IOS) IGTO(I)
         END DO
      END DO
C
      IF (LUONEL .LE. 0)
     &   CALL GPOPEN(LUONEL,'AOONEINT',' ',' ',' ',IDUMMY,.FALSE.)
C
C     0 - unlabeled first record
C     ==========================
C
      REWIND LUONEL
      WRITE (LUONEL) TTITLE
      WRITE (LUONEL) MAXREP+1,(NAOS(I),I=1,MAXREP+1), POTNUC,
     &              (0.D0,I=1,2) ! so record is minimum 32 bytes
#ifndef PRG_DIRAC
C
C     1 - ISORDK
C     ==========
C
      WRITE (LUONEL) TABLE1
      WRITE (LUONEL) DUMMY, DUMMY, DUMMY, DUMMY
      WRITE (LUONEL) CHRN, NUCDEP, COOO
C
C     2 - SCFINP
C     ==========
C
      WRITE (LUONEL) TABLE2
      WRITE (LUONEL) TTITLE
      WRITE (LUONEL) MAXREP+1, (NAOS(I), I = 1,MAXREP+1),
     &               POTNUC,
     &               KMAX, (NCENT(I), I = 1,KMAX),
     &               NBASIS, ( JTRAN(I),             I = 1,NBASIS),
     &                       ((ITRAN(I,J), J = 1,8), I = 1,NBASIS),
     &                       ((CTRAN(I,J), J = 1,8), I = 1,NBASIS),
     &               NBASIS, ( INAMN(I),             I = 1,NBASIS),
     &                       ( IPTYP(I),             I = 1,NBASIS),
     &               DPNUC,
     &               NUCDEP, ((COOO(I,J),J=1,3),I=1,NUCDEP),
     &               IFXYZ, DUMMY, QPOL, QQ, JFXYZ
#endif
C
C     3 - SYMINPUT
C     ============
C
      WRITE (LUONEL) TABLE3
      WRITE (LUONEL) NBASIS,(INAMN(I),I=1,NBASIS),(IGTO(I),I=1,NBASIS),
     &               POTNUC,DUMMY,DUMMY,DUMMY
C
C     4 - TESTDATA
C     ============
C
C     Information for testing magnetic integrals
C
      WRITE (LUONEL) TABLE4
      WRITE (LUONEL) NONTYP, (NONT(I), I = 1, NONTYP),
     &               (IQM(I), (JCO2(J,I,1),J=1,IQM(I)), I = 1, NONTYP),
     &               (IDUMMY,I=1,5)
C
C     5 - SPHERICA
C     ============
C
C     Information for testing spherical transformation
C
      IDOCRT = 1
      IF (.NOT.DOCART) IDOCRT = 0
      WRITE (LUONEL) TABLE5
      WRITE (LUONEL) KMAX, IDOCRT, (NHKT(I),I = 1,KMAX), (IDUMMY,I=1,5)
      WRITE (LUONEL) MAXREP, (ISTBAO(I),I=1,KMAX), (IDUMMY,I=1,6)
      WRITE (LUONEL) MXQN, MXAQN, ((ISYMAO(I,J),I=1,MXQN),J=1,MXAQN)
      WRITE (LUONEL) (CSP(I),I=1,NCSP),(ISPADR(I),I=1,MXQN)
C
      WRITE (LUONEL) TABLE6
      CALL GPCLOSE(LUONEL,'KEEP')
      RETURN
      END
C  /* Deck reapri */
      SUBROUTINE REAPRI
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
C
#include "abainf.h"
#include "cbirea.h"
C
#include "molinp.h"
#include "ccom.h"
#include "nuclei.h"
#include "primit.h"
#include "shells.h"
#include "symmet.h"
#include "chrsgn.h"
#include "aosotr.h"
C
      CALL HEADER('Test output from REAPRI',-1)
      WRITE(LUPRI,'(A,1P,D12.2)') '  THRS  ', THRS
      WRITE(LUPRI,'(A,I5)')    '  NHTYP ', NHTYP
      WRITE(LUPRI,'(A,L5)')    '  DOCART', DOCART
      WRITE(LUPRI,*) 'charge ',(CHARGE(i),i=1,nucind)
      WRITE(LUPRI,*) 'cord1',(cord(1,i),i=1,nucind)
      WRITE(LUPRI,*) 'cord2',(cord(2,i),i=1,nucind)
      WRITE(LUPRI,*) 'cord3',(cord(3,i),i=1,nucind)
      WRITE(LUPRI,*) 'noorbt ',(noorbt(i),i=1,nucind)
      WRITE(LUPRI,*) 'nucind,nucdep',NUCIND, NUCDEP
      WRITE(LUPRI,*) 'nucpre ',(nucpre(i),i=1,nucind)
      WRITE(LUPRI,*) 'nucnum ',((NUCNUM(i,j),i=1,nucind),j=1,8)
      WRITE(LUPRI,*) 'nucdeg ',(nucdeg(i),i=1,nucind)
      WRITE(LUPRI,*) 'istbnu ',(istbnu(i),i=1,nucind)
      WRITE(LUPRI,*) 'ntraco,itraco',NTRACO, (ITRACO(i),i=1,3)
      WRITE(LUPRI,*) 'NATOMS, NFLOAT, NBASIS, NPBAS',
     &          NATOMS, NFLOAT, NBASIS, NPBAS
      WRITE(LUPRI,*) 'namex ',(namex(i),i=1,3*nucind)
      WRITE(LUPRI,*) 'namdep',(namdep(i),i=1,nucdep)
      WRITE(LUPRI,*) 'namdpx',(namdpx(i),i=1,nucdep)
      call header('priexp',-1)
      call output(priexp,1,1,1,npshel,1,mxprim,-1,lupri)
      call header('priccf',-1)
      call output(priccf,1,npshel,1,mxcont,mxprim,mxcont,-1,lupri)
      call header('pricrx',-1)
      call output(pricrx,1,npshel,1,3,mxprim,3,1,-lupri)
      WRITE(LUPRI,*) ' kmax ', kmax
      call header('cent',-1)
      call output(cent,1,kmax,1,3,mxshel,3,1,lupri)
      WRITE(LUPRI,*) ' nhkt ',(nhkt(i),i=1,kmax)
      WRITE(LUPRI,*) ' khkt ',(khkt(i),i=1,kmax)
      WRITE(LUPRI,*) ' kckt ',(kckt(i),i=1,kmax)
      WRITE(LUPRI,*) ' istbao ',(istbao(i),i=1,kmax)
      WRITE(LUPRI,*) ' nuco ',(nuco(i),i=1,kmax)
      WRITE(LUPRI,*) ' jstrt ',(jstrt(i),i=1,kmax)
      WRITE(LUPRI,*) ' nstrt ',(nstrt(i),i=1,kmax)
      WRITE(LUPRI,*) ' mst ',(mst(i),i=1,kmax)
      WRITE(LUPRI,*) ' ncent ',(ncent(i),i=1,kmax)
      WRITE(LUPRI,*) ' share ',(share(i),i=1,kmax)
      WRITE(LUPRI,*) ' nrco ',(nrco(i),i=1,kmax)
      WRITE(LUPRI,*) ' numcf ',(numcf(i),i=1,kmax)
      WRITE(LUPRI,*) ' nbch ',(nbch(i),i=1,kmax)
      WRITE(LUPRI,*) ' kstrt ',(kstrt(i),i=1,kmax)
      WRITE(LUPRI,*) ' segm ',(segm(i),i=1,kmax)
      WRITE(LUPRI,*) ' iptshl ',(iptshl(i),i=1,kmax)
      WRITE(LUPRI,*) ' numcft ',(numcft(i),i=1,kmax)
      WRITE(LUPRI,*) ' sphr ',(sphr(i),i=1,kmax)
      WRITE(LUPRI,*) ' fmult ',(fmult(i),i=0,7)
      WRITE(LUPRI,*) ' pt ',(pt(i),i=0,7)
      WRITE(LUPRI,*) ' mult ',(mult(i),i=0,7)
      WRITE(LUPRI,*) ' maxrep, maxopr ',MAXREP, MAXOPR
      WRITE(LUPRI,*) ' isymax ',((isymax(i,j),i=1,3),j=1,2)
      WRITE(LUPRI,*) ' isymao ',((isymao(i,j),i=1,mxqn),j=1,mxaqn)
      WRITE(LUPRI,*) ' nparsu ',(nparsu(i),i=1,8)
      WRITE(LUPRI,*) ' naos ',(naos(i),i=1,8)
      WRITE(LUPRI,*) ' nparnu ',((nparnu(i,j),i=1,8),j=1,8)
      WRITE(LUPRI,*) ' iptsym ',((iptsym(i,j),i=1,norbs),j=0,7)
      WRITE(LUPRI,*) ' iptcnt ',
     &     (((iptcnt(i,j,k),i=1,3*nucdep),j=0,7),k=1,2)
      WRITE(LUPRI,*) ' ncrrep ',((ncrrep(i,j),i=0,7),j=1,2)
      WRITE(LUPRI,*) ' iptcor ',((iptcor(i,j),i=1,3*nucdep),j=1,2)
      WRITE(LUPRI,*) ' naxrep ',((naxrep(i,j),i=0,7),j=1,2)
      WRITE(LUPRI,*) ' iptax ',((iptax(i,j),i=1,3),j=1,2)
      WRITE(LUPRI,*) ' iptxyz ',(((iptxyz(i,j,k),i=1,3),j=0,7),k=1,2)
      WRITE(LUPRI,*) ' iptnuc ',((iptnuc(i,j),i=1,nucdep),j=0,7)
      call header('ctran',-1)
      call output(ctran,1,nbasis,1,8,mxcorb,8,1,lupri)
      do i = 1,nbasis
         write(LUPRI,*) ' itran  ',i,jtran(i),':',(itran(i,j),j=1,8)
      end do
      WRITE(LUPRI,*) ' iaoao  ',(iaoao(i),i=1,nbasis)
      WRITE(LUPRI,*) ' jaoao1 ',(jaoao(i,1),i=1,nbasis)
      WRITE(LUPRI,*) ' jaoao2 ',(jaoao(i,2),i=1,nbasis)
      return
      end
C  /* Deck zmat */
      SUBROUTINE ZMAT(NONT,KATOM)
C
C  Modified for incorporation in HERMIT/SIRIUS/ABACUS, K.Ruud, May-95
C  Based on program ZMAT from CCQC/ Nov. 1988 (H.J.Aa.Jensen)
C**********************************************************
C*** LAST UPDATED ON SEPTEMBER 11, 1986 BY R. REMINGTON ***
C* REASON: TO OBTAIN HIGHER ACCURACY BY INCREASING THE  ***
C* NUMBER OF DIGITS FOR PI AND CHANGE ATAN TO DATAN #143 **
C**********************************************************
C***LAST UPDATED ON FEBRUARY 01, 1985 BY YUKIO YAMAGUCHI***
C**********************************************************
C   THIS PROGRAM CALCULATES CARTESIAN COORDINATES OF A SYSTEM
C   SPECIFIED BY INTERNAL COORDINATES.
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "dummy.h"
      PARAMETER (D0 = 0.0D0)

      DIMENSION NONT(KATOM)

#include "nuclei.h"
#include "molinp.h"
      COMMON/CBZMAT/IZ(MXCENT,4),BL(MXCENT),ALP(MXCENT),BET(MXCENT)

C     include character ASYMB(0:NSYMB)*2 with atomic symbols
#include <asymb.h>

C
      NCTOT = NCTOT + 1
      NUCIND = NUCIND + 1
      IF (NUCIND .GT. MXCENT) GOTO 5000
      IZ(NUCIND,1) = 0
      IZ(NUCIND,2) = 0
      IZ(NUCIND,3) = 0
      IZ(NUCIND,4) = 0
      BL(NUCIND)=D0
      ALP(NUCIND)=D0
      BET(NUCIND)=D0

      READ (MLINE(NMLINE),'(A4)',IOSTAT=IOS) NAMN(NUCIND)
      IF (IOS.NE.0) THEN
        CALL QUIT('ZMAT: Error in reading NAMN(NUCIND)')
      ENDIF
C
      IOS = 0
      IF(NUCIND .EQ. 1) GO TO 200
      IF(NUCIND .EQ. 2) GO TO 201
      IF(NUCIND .EQ. 3) GO TO 202
      GO TO 203
#if defined (SYS_CRAY) || defined (VAR_NOFREE) || defined (SYS_T3D) || defined (SYS_T90)
  200 ISTART = 5
      CALL FREFRM(MLINE(NMLINE),ISTART,INUM,DUMMY,'INT',IERR)
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,CHARGE(NUCIND),'REA',IERR)
      GO TO 204
  201 ISTART = 5
      CALL FREFRM(MLINE(NMLINE),ISTART,INUM,DUMMY,'INT',IERR)
      CALL FREFRM(MLINE(NMLINE),ISTART,IZ(NUCIND,1),DUMMY,'INT',IERR)
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,BL(NUCIND),'REA',IERR)
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,CHARGE(NUCIND),'REA',IERR)
      GO TO 204
  202 ISTART = 5
      CALL FREFRM(MLINE(NMLINE),ISTART,INUM,DUMMY,'INT',IERR)
      CALL FREFRM(MLINE(NMLINE),ISTART,IZ(NUCIND,1),DUMMY,'INT',IERR)
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,BL(NUCIND),'REA',IERR)
      CALL FREFRM(MLINE(NMLINE),ISTART,IZ(NUCIND,2),DUMMY,'INT',IERR)
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,ALP(NUCIND),'REA',IERR)
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,CHARGE(NUCIND),'REA',IERR)
      GO TO 204
  203 ISTART = 5
      CALL FREFRM(MLINE(NMLINE),ISTART,INUM,DUMMY,'INT',IERR)
      CALL FREFRM(MLINE(NMLINE),ISTART,IZ(NUCIND,1),DUMMY,'INT',IERR)
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,BL(NUCIND),'REA',IERR)
      CALL FREFRM(MLINE(NMLINE),ISTART,IZ(NUCIND,2),DUMMY,'INT',IERR)
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,ALP(NUCIND),'REA',IERR)
      CALL FREFRM(MLINE(NMLINE),ISTART,IZ(NUCIND,3),DUMMY,'INT',IERR)
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,BET(NUCIND),'REA',IERR)
      CALL FREFRM(MLINE(NMLINE),ISTART,IZ(NUCIND,4),DUMMY,'INT',IERR)
      CALL FREFRM(MLINE(NMLINE),ISTART,IDUMMY,CHARGE(NUCIND),'REA',IERR)
#else
  200 READ (MLINE(NMLINE)(5:80),*) INUM,CHARGE(NUCIND)
      GO TO 204
  201 READ (MLINE(NMLINE)(5:80),*,IOSTAT=IOS)
     &     INUM,IZ(NUCIND,1),BL(NUCIND),CHARGE(NUCIND)
      GO TO 204
  202 READ (MLINE(NMLINE)(5:80),*,IOSTAT=IOS)
     &     INUM,IZ(NUCIND,1),BL(NUCIND),IZ(NUCIND,2),
     &     ALP(NUCIND),CHARGE(NUCIND)
      GO TO 204
  203 READ (MLINE(NMLINE)(5:80),*,IOSTAT=IOS)
     &     INUM,IZ(NUCIND,1),BL(NUCIND),IZ(NUCIND,2),
     &     ALP(NUCIND),IZ(NUCIND,3),BET(NUCIND),
     &     IZ(NUCIND,4),CHARGE(NUCIND)
#endif
 204  CONTINUE
      IF (IOS.NE.0) THEN
        CALL QUIT('Error in reading...ZMAT')
      ENDIF
      IF (NAMN(NUCIND)(1:4) .EQ. '    ') THEN
         ICHRG = NINT(CHARGE(NUCIND))
         NAMN(NUCIND)(1:2) = ASYMB(ICHRG)
      END IF
      IF (NAMN(NUCIND)(1:2) .EQ. '  ') THEN
         NAMN(NUCIND)(1:3) = NAMN(NUCIND)(2:4) 
         NAMN(NUCIND)(4:4) = ' '
         GO TO 204
      END IF
!hjaaj : automate numbering of atoms, if not numbered by user
      IF (NAMN(NUCIND)(2:4) .EQ. '   ') THEN
         INUM = MOD(INUM,1000)
         WRITE(NAMN(NUCIND)(2:4),'(I3)') INUM
      ELSE IF (NAMN(NUCIND)(3:4) .EQ. '  ') THEN
         INUM = MOD(INUM,100)
         WRITE(NAMN(NUCIND)(3:4),'(I2)') INUM
      END IF
      NMLINE = NMLINE + 1
      NCLINE(NUCIND) = NMLINE
      NONT(NUCIND) = 1
      NAMEX(3*NUCIND)     = NAMN(NUCIND)//' z'
      NAMEX(3*NUCIND - 1) = NAMN(NUCIND)//' y'
      NAMEX(3*NUCIND - 2) = NAMN(NUCIND)//' x'
C
      RETURN
C
C     Error message:
C
 5000 CONTINUE
        WRITE (LUPRI,'(/A/A,I5)')
     &    ' Too many atomic centers: MXCENT exceed in ZMAT,',
     &    ' Current limit:',MXCENT
        CALL QUIT('*** ERROR *** MXCENT exceeded in ZMAT')
      END
C  /* Deck buildz */
      SUBROUTINE BUILDZ(IPRINT,NONTYP,NSYMOP,NONT,IQM,NBLCK,JCO,
     &           NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,ISGEN,KATOM,
     &           KANG,KBLOCK,KPRIM,KAOVEC,DOOWN,MXSPD,BASREF)
C
C     Modified for HERMIT/SIRIUS/ABACUS -> K.Ruud, May-95
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "codata.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0)
      PARAMETER (DSM = 1.0D-30)
      PARAMETER (AR=PI/180.0D+00)
      PARAMETER (CORMAX = 1.D5)
      LOGICAL   NOORBTS, DOOWN
      COMMON/CBZMAT/IZ(MXCENT,4),BL(MXCENT),ALP(MXCENT),BET(MXCENT)
      DIMENSION A(MXCENT),B(MXCENT),CZ(MXCENT,3),D(MXCENT)
      DIMENSION U1(3),U2(3),U3(3),U4(3),VJ(3),VP(3),V3(3)
      DIMENSION NONT(KATOM),IQM(KATOM,3),JCO(KANG,KATOM,3),
     &          NBLCK(KATOM,3),NUC(KBLOCK,3),NRC(KBLOCK,3),
     &          SEG(KBLOCK,3),
     &          ALPHA(KPRIM,KBLOCK,3),CPRIM(KPRIM,KPRIM,KBLOCK,3),
     &          CPRIMU(KPRIM,KPRIM,KBLOCK,3), CHNON(MXCENT)
      DIMENSION ISGEN(*)
      CHARACTER SPDCAR*1, NAVNK*4, CRT*1
      CHARACTER*80 BASNAM, BASREF(10,KATOM,3)
      DIMENSION IBLOCK(3)
#include "nuclei.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#include "ccom.h"
#include "cbirea.h"
#include "molinp.h"
#include "huckel.h"
#include "cbisol.h"
    1 FORMAT(/'   Internal coordinates (Z matrix)'/)
    2 FORMAT(I5,62X,F15.5)
    3 FORMAT(I5,I5,F14.7,43X,F15.5)
    4 FORMAT(I5,I5,F14.7,I5,F14.7,24X,F15.5)
    5 FORMAT(I5,I5,F14.7,I5,F14.7,I5,F14.7,I5,F15.5)
C   6 FORMAT(//25X,' CARTESIAN COORDINATES in Angstrom',
C    &        /14X,'X',24X,'Y',24X,'Z'/)
    7 FORMAT(I5,3(F20.10,5X),F15.5)
    8 FORMAT(//25X,' Cartesian coordinates in bohr',
     &        /14X,'X',24X,'Y',24X,'Z'/)
    9 format(3f20.10)
C
      IF (SOLVNT) NUCIND = NUCIND - 1
C
      IF (IPRINT .GT. 2) THEN
         WRITE(LUPRI,1)
         WRITE(LUPRI,2) 1,CHARGE(1)
         IF (NUCIND .LE. 1) GOTO 201
         WRITE(LUPRI,3) 2,IZ(2,1),BL(2),CHARGE(2)
         IF(NUCIND.LE.2) GO TO 201
         WRITE(LUPRI,4) 3,IZ(3,1),BL(3),IZ(3,2),ALP(3),CHARGE(3)
         IF(NUCIND.LE.3) GO TO 201
         DO I=4,NUCIND
            WRITE(LUPRI,5) I,IZ(I,1),BL(I),IZ(I,2),ALP(I),
     &                     IZ(I,3),BET(I),IZ(I,4),CHARGE(I)
         END DO
 201  CONTINUE
      END IF
      DO J = 1,3
         DO I = 1,NUCIND
            CZ(I,J) = D0
         END DO
      END DO
      DO I = 1,NUCIND
         ALP(I) = ALP(I)*AR
         BET(I) = BET(I)*AR
      END DO
      CZ(2,3) = BL(2)
      IF (NUCIND .GE. 3) THEN
         CZ(3,1) = BL(3)*DSIN(ALP(3))
         IF ((IZ(3,1)-1) .EQ. 0) THEN
            CZ(3,3) = BL(3)*DCOS(ALP(3))
         ELSE
            CZ(3,3) = CZ(2,3)-BL(3)*DCOS(ALP(3))
         END IF
         DO I = 4,NUCIND
            IF((ABS(CZ(I-1,1))-1.0D-6) .LT. D0) THEN
               CZ(I,1) = BL(I)*DSIN(ALP(I))
               ITEMP = IZ(I,1)
               JTEMP = IZ(I,2)
               CZ(I,3) = CZ(ITEMP,3) - BL(I)*DCOS(ALP(I))*
     &                   DSIGN(D1,CZ(ITEMP,3) - CZ(JTEMP,3))
            ELSE
               GOTO 90
            END IF
         END DO
 90      CONTINUE
         K=I
         IF (K .LE. NUCIND) THEN
            DO 250 J = K,NUCIND
               CAJ = COS(ALP(J))
               SAJ = SIN(ALP(J))
               CBJ = COS(BET(J))
               SBJ = SIN(BET(J))
               IF(IZ(J,4) .EQ. 0) THEN
                  CALL V3VEC(U1,CZ,IZ(J,2),IZ(J,3))
                  CALL V3VEC(U2,CZ,IZ(J,1),IZ(J,2))
                  CALL V3PROD(VP,U1,U2)
                  R = SQRT(D1 - (U1(1)*U2(1) + U1(2)*U2(2)
     &                 + U1(3)*U2(3))**2)
                  DO I = 1,3
                     U3(I) = VP(I)/R
                  END DO
                  CALL V3PROD(U4,U3,U2)
                  DO I = 1,3
                     VJ(I) = BL(J)*(-U2(I)*CAJ + U4(I)*SAJ*CBJ
     &                    + U3(I)*SAJ*SBJ)
                     ITEMP = IZ(J,1)
                     CZ(J,I)=VJ(I)+CZ(ITEMP,I)
                  END DO
               ELSE
                  IF (IABS(IZ(J,4)) .EQ. 1) THEN
                     CALL V3VEC(U1,CZ,IZ(J,1),IZ(J,3))
                     CALL V3VEC(U2,CZ,IZ(J,2),IZ(J,1))
                     ZETA = -(U1(1)*U2(1) + U1(2)*U2(2) +U1(3)*U2(3))
                     A(J) = (-CBJ + ZETA*CAJ)/(D1 - ZETA*ZETA)
                     B(J) = (CAJ - ZETA*CBJ)/(D1 - ZETA*ZETA)
                     R=D0
                     GAMMA=PI/D2
                     IF (ZETA .LT. D0) THEN
                        R=PI
                        GAMMA = DATAN(SQRT(D1-ZETA*ZETA)/ZETA)+R
                     ELSE IF (ZETA .GT. D0) THEN
                        GAMMA = DATAN(SQRT(D1-ZETA*ZETA)/ZETA)+R
                     END IF
                     D(J) = D0
                     IF (ABS(GAMMA+ALP(J)+BET(J)-D2*PI) .GT. 1.0D-6)
     &                  THEN
                        D(J) = IZ(J,4)*(SQRT(D1+A(J)*CBJ-B(J)*CAJ))
     &                       /SQRT(D1-ZETA*ZETA)
                     END IF
                     CALL V3PROD(V3,U1,U2)
                     DO I = 1,3
                        U3(I) = A(J)*U1(I)+B(J)*U2(I)+D(J)*V3(I)
                        VJ(I)=BL(J)*U3(I)
                        ITEMP=IZ(J,1)
                        CZ(J,I)=VJ(I)+CZ(ITEMP,I)
                     END DO
                  ELSE
                     CALL V3VEC(U1,CZ,IZ(J,1),IZ(J,3))
                     CALL V3VEC(U2,CZ,IZ(J,2),IZ(J,1))
                     ZETA = -(U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3))
                     CALL V3PROD(V3,U1,U2)
                     V3MAG = SQRT(V3(1)*V3(1)+V3(2)*V3(2)
     &                    + V3(3)*V3(3))
                     A(J) = V3MAG*CBJ/(D1-ZETA*ZETA)
                     B(J) = SQRT((D1-CAJ*CAJ-A(J)*CBJ*V3MAG)
     &                    /(D1-ZETA*ZETA))
                     IF (IZ(J,4) .NE. 2) B(J)=-B(J)
                     D(J) = B(J)*ZETA+CAJ
                     DO I = 1,3
                        U3(I) = B(J)*U1(I)+D(J)*U2(I)+A(J)*V3(I)
                        VJ(I) = BL(J)*U3(I)
                        ITEMP = IZ(J,1)
                        CZ(J,I) = VJ(I)+CZ(ITEMP,I)
                     END DO
                  END IF
               END IF
 250        CONTINUE
         END IF
      END IF
C
C   ELIMINATE IMAGINARY ATOM(S)
C   Hmmm, what do we do here? Remove them or keep them as dummy atoms?
C   Check to see if there is a basis set associated with the center?
C   Is it at all possible to place a dummy atom on top of an already
C   existing center in internal coordinates? K.Ruud, May-95
C
C      NATOMS=0
C      DO 290 I=1,NUCIND
C         IF(ANZ(I) .NE. D0) THEN
C            NATOMS=NATOMS+1
C            CHRG(NATOMS) = DFLOAT( ANZ(I) )
C            DO 280 J=1,3
C               C(NATOMS,J)=CZ(I,J)
C 280        CONTINUE
C         END IF
C 290  CONTINUE
C
C
      IF (IPRINT .GT. 2) THEN
         WRITE(LUPRI,8)
         DO I=1,NUCIND
            WRITE(LUPRI,7) I,(CZ(I,J)/XTANG,J=1,3),CHARGE(I)
         END DO
      END IF
C
C     Move information to HERMIT common blocks
C
      DO I=1,NUCIND
         DO J = 1, 3
            CORD(J,I) = CZ(I,J)/XTANG
            IF (ABS(CORD(J,I)).GT. CORMAX) GOTO 5010
         END DO
      END DO
C
C     Update MOLECULE.INP. In order to work properly with the symmetry
C     detection routines we need to sort atoms with similar charge
C     together in a block.
C
      ITYP = 1
      NONT(ITYP) = 1
      CHNON(ITYP) = CHARGE(1)
      DO 600 I = 2, NUCIND
         DO 601 J = I-1, 1, -1
            IF (CHARGE(I) .EQ. CHARGE(J)) THEN
               NONTYP = NONTYP - 1
               DO L = 1, ITYP
                  IF (CHARGE(I) .EQ. CHNON(L)) NONT(L) = NONT(L) + 1
               END DO
               CHARGK = CHARGE(I)
               XCORD  = CORD(1,I)
               YCORD  = CORD(2,I)
               ZCORD  = CORD(3,I)
               NAVNK  = NAMN(I)
               DO K = I - 1, J + 1, -1
                  CHARGE(K + 1) = CHARGE(K)
                  CORD(1,K + 1) = CORD(1,K)
                  CORD(2,K + 1) = CORD(2,K)
                  CORD(3,K + 1) = CORD(3,K)
                  NAMN(K + 1)   = NAMN(K)
               END DO
               CHARGE(J + 1) = CHARGK
               CORD(1,J + 1) = XCORD
               CORD(2,J + 1) = YCORD
               CORD(3,J + 1) = ZCORD
               NAMN(J + 1)   = NAVNK
               GOTO 600
            END IF
 601     CONTINUE
         ITYP = ITYP + 1
         NONT(ITYP) = 1
         CHNON(ITYP) = CHARGE(I)
 600  CONTINUE
C     define atom number in IZATOM /hjaaj Mar 2004
C     TODO : define multiple basis sets and point charges
      DO I = 1, NUCIND
         IZATOM(I) = NINT(CHARGE(I))
      END DO
      IF (ITYP.NE.NONTYP) CALL QUIT('Problem ITYP.ne.NONTYP in BUILDZ')
      IF (DOCART) THEN
         CRT = 'C'
      ELSE IF (DOOWN) THEN
         CRT = 'X'
      ELSE
         CRT = ' '
      END IF
      WRITE(MLINE(5),'(A1,I4,I3,I2,10X,1P,D10.2)')
     &     CRT,NONTYP,KCHARG,NSYMOP,THRS
      NMLINE = 6
      NHTYP  = 0
      NSETS  = 1
      IF (DOHUCKEL) THEN
         NSETS = NSETS + 1
         ISETHUCKEL = NSETS
      END IF
      IBLOCK(1:NSETS) = 1

      IATOM = 1
      DO 605 I = 1, NONTYP
         NOORBTS = .TRUE.
         WRITE(MLINE(NMLINE),'(6X,F4.0,I5)') CHARGE(IATOM),NONT(I)
         NMLINE = NMLINE + 1
         DO 606 J = 1, NONT(I)
            NCLINE(IATOM) = NMLINE
            WRITE(MLINE(NMLINE)(1:4),'(A4)') NAMN(IATOM)
            IF (ABS(CORD(1,IATOM)).ge.100.0D0 .OR.
     &          ABS(CORD(2,IATOM)).ge.100.0D0 .OR.
     &          ABS(CORD(3,IATOM)).ge.100.0D0) THEN
               WRITE(MLINE(NMLINE)(5:76),'(3F20.10,12X)')
     &            CORD(1,IATOM),CORD(2,IATOM),CORD(3,IATOM)
            ELSE
               WRITE(MLINE(NMLINE)(5:76),'(3F20.15,12X)')
     &            CORD(1,IATOM),CORD(2,IATOM),CORD(3,IATOM)
            END IF
            NMLINE = NMLINE + 1
            IATOM = IATOM + 1
 606     CONTINUE
C
C     Add basis set to newly order atomic set
C
C     TODO : fix code for LMULBS true below
C
         DO J = 1,NSETS
            Q    = CHARGE(IATOM-1)
            QEFF = Q
C           ... QEFF is effective charge for ECP, ECP not implemented
C               for ZMAT (yet?)
            IF (J .EQ. ISETHUCKEL) THEN
               BASNAM = 'HUCKEL'
            ELSE
               BASNAM = MULNAM(J)
            END IF
            CALL BASLIB(IQM(I,J),JCO(1,I,J),NUC(IBLOCK(J),J),
     &                  NRC(IBLOCK(J),J),SEG(IBLOCK(J),J),
     &                  ALPHA(1,IBLOCK(J),J),
     &                  CPRIM(1,1,IBLOCK(J),J),CPRIMU(1,1,IBLOCK(J),J),
     &                  NBLOCK,ISGEN(IBLOCK(J)),KAOVEC,KPRIM,Q,QEFF,
     &                  DSM,UNCONT,BASNAM,BASREF(1,I,J),IPRINT)
            NBLCK(I,J) = NBLOCK
            IBLOCK(J) = IBLOCK(J) + NBLCK(I,J)
            NOORBTS = NOORBTS .AND. IQM(I,J).EQ.0
            IF (IQM(I,1) .GT. 0) THEN
               NHTYP = MAX(NHTYP,IQM(I,J))
               IF (NHTYP .GT. MXSPD) GOTO 5000
            END IF
         END DO
         DO K = 1,NONT(I)
            NOORBT(K + IATOM - 2) = NOORBTS
         END DO
 605  CONTINUE
      NMLINE = NMLINE - 1
      IF (SOLVNT) NUCIND = NUCIND + 1
      RETURN
 5000 CONTINUE
         WRITE (LUPRI,'(6X,A,I3,3A/9X,2A)')
     &      '*  Input specifies highest orbital of atomic type ',
     &      I,' AS "',SPDCAR(NHTYP - 1),'".',
     &      ' Highest allowed orbital in this version: ',
     &      SPDCAR(MXSPD - 1)
         IF (NHTYP.GT.MXQN) WRITE (LUPRI,'(9X,2(A,I3),A)')
     &      ' Increase MXQN from',MXQN,' to',NHTYP,' and recompile.'
         CALL QUIT('Too high angular quantum no. specified in input.')
 5010 CONTINUE
        WRITE (LUPRI,'(A,E12.6,A/A/A,E12.6)')
     &    ' Atomic coordinate ',CORD(J,NUCIND),
     &    ' too large in BUILDZ',
     &    ' Note: Program is unstable for large coordinates.',
     &    ' Maximum coordinate value:',CORMAX
        CALL QUIT('*** ERROR: Atomic coordinate too large in BUILDZ')
      END
C  /* Deck V3vec */
      SUBROUTINE V3VEC(U,C,J,K)
#include "implicit.h"
#include "mxcent.h"
      PARAMETER (D0 = 0.0D0)
      DIMENSION C(MXCENT,3),R(3),U(3)
C
      R2=D0
      DO 101 I=1,3
         R(I)=C(J,I)-C(K,I)
         R2=R2+R(I)*R(I)
 101  CONTINUE
      R2=SQRT(R2)
      DO 102 I=1,3
         U(I)=R(I)/R2
  102 CONTINUE
      RETURN
      END
C  /* Deck v3prod */
      SUBROUTINE V3PROD(VP,X,Y)
#include "implicit.h"
      DIMENSION VP(3),X(3),Y(3)
C
      VP(1)=X(2)*Y(3)-X(3)*Y(2)
      VP(2)=X(3)*Y(1)-X(1)*Y(3)
      VP(3)=X(1)*Y(2)-X(2)*Y(1)
      RETURN
      END
C/*  Deck frefrm */
      SUBROUTINE FREFRM(STRING,ISTART,IVAL,RVAL,TYPE,IERR)
C
C     Simulate free format input for Cray computers when reading from an
C     internal file, K.Ruud-Feb.97
C
#include <implicit.h>
#include <priunit.h>
C      
      CHARACTER*80 STRING
      CHARACTER TYPE*3, FRMT*7
C
      IERR = 0
      IF (TYPE .EQ. 'INT') THEN
         IVAL = 0
      ELSE
         RVAL = 0.0D0
      END IF
      DO 10 IPOS = ISTART, 80
         IF (STRING(IPOS:IPOS) .NE. ' ') GO TO 20
         ! leave cycle when found first non-white character
 10   CONTINUE
C
 20   CONTINUE
      IEND = INDEX(STRING(IPOS:),' ')
      ILEN = IEND - ISTART - 1 + IPOS
      IF (TYPE .EQ. 'INT') THEN
         IF (ILEN .GT. 9) THEN
            WRITE (FRMT,'(2X,A2,I2,A1)') '(I',ILEN,')'
         ELSE
            WRITE (FRMT,'(3X,A2,I1,A1)') '(I',ILEN,')'
         END IF
         READ (STRING(ISTART:80),FRMT,ERR=9000) IVAL
      ELSE IF (TYPE .EQ. 'REA') THEN
         IF (ILEN .GT. 9) THEN
            WRITE (FRMT,'(A2,I2,A3)') '(F',ILEN,'.0)'
         ELSE
            WRITE (FRMT,'(1X,A2,I1,A3)') '(F',ILEN,'.0)'
         END IF
         READ (STRING(ISTART:80),FRMT,ERR=9000) RVAL
      ELSE
         WRITE (LUPRI,'(/,A)') 'Illegal variable type in FREFRM'
         CALL QUIT('Illegal variable type in FREFRM')
      END IF
      ISTART = IEND + IPOS
      GO TO 9010
 9000 CONTINUE
      IERR = 1
 9010 CONTINUE
      RETURN
      END
C  /* Deck molout */
      SUBROUTINE MOLOUT(NONTYP,NONT,IQM,NBLCK,JCO,NUC,NRC,SEG,
     &     IPCON,KATOM,KANG,KBLOCK,KPRIM,CPRIMU,NRMPRI,KASYM,
     &     NSYMOP)
C
C     This routine based on TECORBOUT dumps a MOLECULE.INP file suitable
C     for direct processing by DALTON
C     Now also dumps a file MNF.INP for use for SO mean-field approximation
C     Written by K.Ruud 160297 
C
#include "implicit.h"
#include "priunit.h"
#include <dummy.h>
#include "pi.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0)
C    
#include <ccom.h>
#include <cbirea.h>
#include <nuclei.h>
#include <primit.h>
#include <gnrinf.h>
#include <molinp.h>
      LOGICAL SEG, NRMPRI
      CHARACTER KASYM(3,3)*1, CRT*1
      DIMENSION NONT(KATOM),IQM(KATOM),NBLCK(KATOM),
     &     JCO(KANG,KATOM),NUC(KBLOCK),NRC(KBLOCK),
     &     CPRIMU(KPRIM,KPRIM,KBLOCK),SEG(KBLOCK)
      DIMENSION IPCON(KPRIM)
      CHARACTER*8 LTIME, LDATE
      CHARACTER*1 SPD(8)
      DATA SPD/'S','P','D','F','G','H','I','J'/
C
      LUTMP1 = -1
      LUTMP2 = -1
      LUSHB  = -1
      CALL GPOPEN(
     &          LUTMP1,'DALTON.BAS',' ',' ','FORMATTED',IDUMMY,.FALSE.)
      CALL GPOPEN(LUTMP2, 'MNF.INP',' ',' ','FORMATTED',IDUMMY,.FALSE.)
      CALL GPOPEN(LUSHB, 'LUSH.BAS',' ',' ','FORMATTED',IDUMMY,.FALSE.)
      REWIND (LUSHB)
C    
      IPRIMD = 0
      IPRIM =  0
      ICENT  = 0
      JBLOCK = 0
C hj-aug99: find title lines, first skip initial comments
      ILINE = 1
   10 IF((MLINE(ILINE)(1:1).EQ.'#') .OR. (MLINE(ILINE)(1:1).EQ.'!'))THEN
         ILINE = ILINE + 1
         GO TO 10
      END IF
Chj:  ILINE must now point to line with first effective line (CARD 1)
      CALL GETDAT(LDATE,LTIME)
      WRITE (LUTMP1,'(2A,1X,A/A/2A)')
     &   '# DALTON.BAS with INTGRL directive created ',LDATE,LTIME,
     &   '# Directive in input mol file was:',
     &   '# > ',MLINE(ILINE)(1:76)
      IF (BASIS) THEN
         WRITE (LUTMP1,'(2A)')
     &   '# > ',MLINE(ILINE+1)(1:76)
         ILINE = ILINE + 2
      ELSE
         ILINE = ILINE + 1
      END IF
Chj:  skip any additional comment lines before title lines
   20 IF((MLINE(ILINE)(1:1).EQ.'#') .OR. (MLINE(ILINE)(1:1).EQ.'!'))THEN
         ILINE = ILINE + 1
         GO TO 20
      END IF
      CRT = ' '
      IF (DOCART) CRT = 'C'
      WRITE (LUTMP1,'(A/A/A)')
     &   'INTGRL',MLINE(ILINE),MLINE(ILINE + 1)
      ILINE = ILINE + 2
      WRITE (LUTMP1,'(A1,I4,I3,I2,9A1,1X,D10.2)')
     &     CRT,NONTYP,KCHARG,NSYMOP,((KASYM(I,J),I=1,3),J=1,3),THRS
      ISYCNT = 1
      WRITE(LUSHB,'(A/A)')'$BASIS',MLINE(2)
      DO 100 I = 1, NONTYP
         WRITE (LUTMP1,'(6X,F4.0,24(I5))') CHARGE(ICENT+1),NONT(I),
     &                             IQM(I),(JCO(J,I), J=1, IQM(I))
         DO N = 1, NONT(I)
            WRITE (LUTMP1,'(A4,3F20.13)') NAMN(ISYCNT),
     &           (CORD(J,ISYCNT), J = 1,3)
            WRITE(LUSHB,'(/F16.1,F19.10,2F20.10)')
     &         CHARGE(ICENT+1),(CORD(J,ISYCNT),J=1,3)
            ISYCNT = ISYCNT + 1
         END DO
         MNCENT = ICENT + 1
         DO 110 N = 1, NONT(I)
            ICENT = ICENT + 1
            NDEG  = NUCDEG(ICENT)
            KBCH  = JBLOCK
            WRITE (LUTMP2,*) CHARGE(MNCENT),IQM(I)-1
            II = 0
            DO 200 J = 1, IQM(I)
               DO 210 K = 1, JCO(J,I)
                  KBCH = KBCH + 1
                  NNUC  = NUC(KBCH)
                  NNRC  = NRC(KBCH)
                  IF (NNUC .EQ. 0) GO TO 200
                  WRITE (LUTMP2,*) NNUC,NNRC
                  ITYP = NHKOFF(J)
                  IPSTRT = IPRIM + 1
                  IPRIM =  IPRIM + NNUC
                  ITYP = ITYP + 1
                  IPRIMD = IPRIMD + 1
                     IF (N.EQ.1) THEN
                     WRITE (LUTMP1,'(A1,I4,I5)') 'H', NNUC, NNRC
                     ENDIF
                     DO 410 M = 1, NNUC
                        IPRIMD = IPRIMD + 1
                        IF (NNRC.LE.3) THEN
                        IF (N.EQ.1) THEN
                           WRITE (LUTMP1,1050) PRIEXP(IPSTRT-1+M),
     &                          (CPRIMU(M,MM,KBCH),MM=1,NNRC)
                        ENDIF
                           WRITE (LUTMP2,1050) PRIEXP(IPSTRT-1+M),
     &                          (CPRIMU(M,MM,KBCH),MM=1,NNRC)
                        ELSE
                        IF (N.EQ.1) THEN
                           WRITE (LUTMP1,1050) PRIEXP(IPSTRT - 1 + M),
     &                          (CPRIMU(M,MM,KBCH),MM = 1,3)
                        ENDIF
                           WRITE (LUTMP2,1050) PRIEXP(IPSTRT - 1 + M),
     &                          (CPRIMU(M,MM,KBCH),MM = 1,3)
                           DO NLOOP = 1, INT(NNRC/3) - 1
                           IF (N.EQ.1) THEN
                              WRITE (LUTMP1,1040) (CPRIMU(M,MM,KBCH),
     &                             MM = NLOOP*3 + 1, (NLOOP + 1)*3)
                           ENDIF
                              WRITE (LUTMP2,1040) (CPRIMU(M,MM,KBCH),
     &                             MM = NLOOP*3 + 1, (NLOOP + 1)*3)
                           END DO
                           IF (MOD(NNRC,3) .NE. 0) THEN
                           IF (N.EQ.1) THEN
                              WRITE (LUTMP1,1040) (CPRIMU(M,MM,KBCH),
     &                             MM = NNRC - MOD(NNRC,3) + 1, NNRC)
                           ENDIF
                              WRITE (LUTMP2,1040) (CPRIMU(M,MM,KBCH),
     &                             MM = NNRC - MOD(NNRC,3) + 1, NNRC)
                           END IF
                        ENDIF
 410                 CONTINUE
                  DO MM=1,NNRC
                     II=II+1
                     WRITE(LUSHB,'(I5,A5,I5)')II, SPD(J), NNUC
                     DO M=1,NNUC
                        WRITE(LUSHB,'(I5,D15.9,D20.10)')
     &                     M,PRIEXP(IPSTRT+M-1),CPRIMU(M,MM,KBCH)
                     END DO
                  END DO
 210           CONTINUE
 200        CONTINUE
 110     CONTINUE
         JBLOCK = JBLOCK + NBLCK(I)
 100  CONTINUE
C
      WRITE (LUTMP2,'(A)') 'END'
      CALL GPCLOSE(LUTMP1,'KEEP')
      CALL GPCLOSE(LUTMP2,'KEEP')
      CALL GPCLOSE(LUSHB ,'KEEP')
 1010 FORMAT(/2X,A4,3X,A4,3X,A10)
 1040 FORMAT(20X,3F20.10)
 1050 FORMAT(F20.10,3F20.10)
      RETURN
      END
C  /* Deck crtest */
      SUBROUTINE CRTEST(DOCART)
#include <implicit.h>
#include <priunit.h>
      LOGICAL DOCART
#include <cbieri.h>
C
      IF (RUNERI.AND.DOCART) THEN
         NWARN = NWARN + 1
         RUNERI = .FALSE.
         WRITE (LUPRI,'(A)')
     &' INFO: Only spherical-harmonic basis is implemented in ERI.'//
     &' HERMIT module (TWOINT) is used instead.'
      END IF
      RETURN
      END
C  /* Deck isomas */
      FUNCTION ISOMAS(ICHARG,MASSNM)
C
C     Function to switch from mass number to isotope number sorted
C     according to abundance, K.Ruud-02
C
#include <implicit.h>
#include <priunit.h>
C
      CALL QENTER('ISOMAS')
      IORD = 0
      DO I = 1, 5
         QMASS = DISOTP(ICHARG,I,'MASS')
         IF (ANINT(QMASS) .EQ. MASSNM) IORD = I
      END DO
      IF (IORD .EQ. 0) THEN
         WRITE (LUPRI,'(/A,I4,A,I4)') 'ERROR: unknown mass',MASSNM,
     &        ' for atom with charge',ICHARG
         CALL QUIT('Unknown mass for chosen atomic charge')
      ELSE
         ISOMAS = IORD
      END IF
      CALL QEXIT('ISOMAS')
      RETURN
      END
C  /* Deck set_carmom_2 */
      SUBROUTINE SET_CARMOM_2
C
C     Enforce computation of Cartesian second moments.
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <mxcent.h>
#include <cbiher.h>
      CARMOM = .TRUE.
      ONEPRP = .TRUE.
      IORCAR = MAX(2,IORCAR)
      RETURN
      END
C  /* Deck READ_basmul */
      SUBROUTINE READ_BASMUL(INPUT_LINE)
C
C     Read basis set labels from input line.
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include "implicit.h"
#include "priunit.h"
#include "cbirea.h"
      CHARACTER*80 INPUT_LINE
      LOGICAL NEXTMB

      N = 0
      NEXTMB = .FALSE.
      DO I=1,80
       IF (INPUT_LINE(I:I) .NE. ' ' .AND. .NOT. NEXTMB) THEN
        NEXTMB = .TRUE.
        K = I
       ELSE IF (INPUT_LINE(I:I) .EQ. ' ' .AND. NEXTMB) THEN
        L = I - 1
        NEXTMB = .FALSE.
        N = N + 1
        IF (N .GT. MXMULB) CALL QUIT('Too many multiple basis sets')
        MULNAM(N)(1:I-K) = INPUT_LINE(K:L)
        IF (IPREAD .GT. 1) THEN
          WRITE (LUPRI,'(/A,I2,A)') ' Basis set',N,' is  "'//
     &    MULNAM(N)(1:I-K)//'" from the basis set library.'
        END IF
       END IF
      END DO
      NMULBS = N

      IF (.NOT. LMULBS .AND. N.GT.1) THEN
         WRITE (LUPRI,'(/A/5X,A/A)')
     &      ' INFO: basis set line contains more than one basis set:',
     &      INPUT_LINE(1:L),
     &      ' INFO: ".R12AUX" activated automatically.'
         LMULBS = .TRUE.
      END IF
      RETURN
      END
      SUBROUTINE LINE4(LINE,NONTYP,NSYMOP,CRT,KCHARG,THRS,ADDSYM,KASYM,
     &                 ID3,DOCART,DOOWN)
C
C     Read in the fourth input line using the new input scheme
C
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <maxorb.h>
#include <infpar.h>
      LOGICAL ADDSYM, DOCART, DOOWN
      CHARACTER*80 LINE
      CHARACTER CRT*1, ID3*1, SYMTXT*2
      CHARACTER*1 KASYM(3,3)
C
C     Initialize
C
      NONTYP = 0
      DOCART = .FALSE.
      DOOWN  = .FALSE.
C
C     Number of different atom types
C
#ifdef VAR_DEBUG
C     to test if .mol input transferred OK to slaves
      print *, mytid,' LINE4: ',LINE
#endif

      IPOS = INDEX(LINE,'ATO')
      IF (IPOS .NE. 0) THEN
         IPOS2 = INDEX(LINE(IPOS:),'=')
         IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 10)) THEN
            WRITE (LUPRI,*) 'Incorrect input for number of atom types.'
            WRITE (LUPRI,*) 'Format is "Atomtypes=?"'
            WRITE (LUPRI,*) LINE
            CALL QUIT('Incorrect input for # atomtypes')
         ELSE
            READ (LINE((IPOS+IPOS2):),*) NONTYP
         END IF
         DOCART = .FALSE.
C
C     Kind of basis set
C
         I = 0
         IPOS = INDEX(LINE,'CAR')
         IF (IPOS .NE. 0) THEN
            I = I + 1
            DOCART = .TRUE.
         END IF
         IPOS = INDEX(LINE,'OWN')
         IF (IPOS .NE. 0) THEN
            I = I + 1
            DOOWN = .TRUE.
         END IF
         IPOS = INDEX(LINE,'SPH')
         IF (IPOS .NE. 0) THEN
            I = I + 1
            DOOWN  = .FALSE.
            DOCART = .FALSE.
         END IF
         IF (I .GT. 1) THEN
            WRITE (LUPRI,'(//A/A)')
     &      ' .mol input error: mixed basis types in this line:',LINE
            CALL QUIT('.mol input ERROR: mixed basis types')
         END IF
C
C     Charge of molecule
C     
         IPOS = INDEX(LINE,'CHA')
         IF (IPOS .NE. 0) THEN
            IPOS2 = INDEX(LINE(IPOS:80),'=')
            IF ((IPOS2 .EQ. 0) .OR. (IPOS2 .GT. 7)) THEN
               WRITE (LUPRI,*) 'Incorrect input for molecular charge'
               WRITE (LUPRI,*) 'Format is "Charge=?"'
               WRITE (LUPRI,*) LINE
               CALL QUIT('Incorrect input for molecular charge')
            ELSE
               READ (LINE(IPOS+IPOS2:),*) KCHARG
            END IF
         ELSE
            KCHARG = 0
         END IF
C     
C     Angstrom?
C
         ID3 = ' '
         IPOS = INDEX(LINE,'ANG') 
         IF (IPOS .NE. 0) ID3 = 'X'
C     
C     Symmetry generators
C     
         DO J = 1, 3
            DO I = 1, 3
               KASYM(I,J) = ' '
            END DO
         END DO
         IPOS = INDEX(LINE,'GEN')
         IPSO = INDEX(LINE,'NOS')
         IF (IPOS .NE. 0 .AND. IPSO .NE. 0) THEN
            WRITE (LUPRI,*) 'Incorrect symmetry input. "Nosymmetry" '//
     &           'requested, but generators given.'
            WRITE (LUPRI,*) LINE
            CALL QUIT('Inconsistent symmetry input')
         ELSE IF (IPOS .NE. 0) THEN
            IPOS2 = INDEX(LINE(IPOS:),'=')
            IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 11)) THEN
               WRITE (LUPRI,*) 'Incorrect input for symmetry generators'
               WRITE (LUPRI,*) 'Format is "Generators=?"'
               WRITE (LUPRI,*) LINE
               CALL QUIT('Incorrect input for symmetry generators')
            ELSE
               READ (LINE((IPOS+IPOS2):),*) NSYMOP
            END IF
            IF (NSYMOP .GT. 0 .AND. NSYMOP .LT. 4) THEN
C     
C     For each generator, determine number of axis changing sign
C
               ISTART = IPOS + IPOS2 + 1
               DO I = 1, NSYMOP
                  DO ISTP = ISTART, 80
                     IF (LINE(ISTP:ISTP) .NE. ' ') GO TO 20
                  END DO
 20               CONTINUE
                  IEND = INDEX(LINE(ISTP:),' ')
                  ILEN = IEND - 1
                  DO J = 1, ILEN
                     IWHERE = ISTP + J - 1
                     KASYM(J,I) = LINE(IWHERE:IWHERE)
                  END DO
                  ISTART = ISTP + IEND
               END DO
            ELSE
               IF (NSYMOP .GT. 3) THEN
                  WRITE (LUPRI,*) 'Too many (>3) group generators'
                  WRITE (LUPRI,*) LINE
                  CALL QUIT('Too many (>3) group generators')
               END IF
            END IF
            ADDSYM = .FALSE.
         ELSE IF (IPSO .NE. 0) THEN
            ADDSYM = .FALSE.
            NSYMOP = 0
         ELSE
            ADDSYM = .TRUE.
            NSYMOP = 0
         END IF
C
C     Change of integral threshold?
C
         IPOS = INDEX(LINE,'INT')
         IF (IPOS .NE. 0) THEN
            IPOS2 = INDEX(LINE(IPOS:),'=')
            IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 11)) THEN
               WRITE (LUPRI,*) 'Incorrect input for integral threshold'
               WRITE (LUPRI,*) 'Format is "Integrals=?"'
               WRITE (LUPRI,*) LINE
               CALL QUIT('Incorrect input for integral threshold')
            ELSE
               CALL FREFRM(LINE,IPOS+IPOS2,IDUMMY,THRS,'REA',IERR)
            END IF
         ELSE
C           ... tell READIN that user hasn't specified THRS
            THRS = -1.0D0
         END IF
      ELSE
C
C***********************************************************************
C     
C     Since the number of atom types have not been specified,
C     assume old input format
C
         READ (LINE,'(BN,A1,I4,I3,A2,10A1,D10.2,6I5)',ERR=199) CRT,
     &        NONTYP,KCHARG,SYMTXT,((KASYM(I,J),I=1,3),J=1,3),ID3,THRS
         IF (LINE(21:30) .EQ. '          ') THRS = -1.0D0
C        ... tell READIN that user hasn't specified THRS
         ADDSYM = (SYMTXT(1:2) .EQ. '  ')
C
C        Also automatic symmetry detection for Dinfh or Cinfv
         ADDSYM = (SYMTXT(1:1) .EQ. 'D') .OR. ADDSYM
         ADDSYM = (SYMTXT(1:1) .EQ. 'C') .OR. ADDSYM
         IF (.NOT. ADDSYM) READ(SYMTXT,'(I2)') NSYMOP
         DOCART = CRT .EQ. 'C' .OR. CRT .EQ. 'c'
         DOOWN  = CRT .EQ. 'X'
      END IF
      IF (ADDSYM) NSYMOP = 0
      RETURN
 199  CONTINUE
      WRITE (LUPRI,'(/A)') ' Error in the determination of the number'
     &   //' of different atomic types'
      WRITE (LUPRI,*) "Correct input structure is: Atomtypes=???"
      WRITE (LUPRI,*) LINE
      CALL QUIT('Error in determining the number of different '//
     &     'atom types')
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE LINE4W(TMPLN,NONTYP,NSYMOP,KCHARG,THRS,AUTOSY,NOSYM,
     &                  KASYM,ID3,DOCART,DOOWN)
C
C     Write out molecule specific information in new input format
C
#include <implicit.h>
      CHARACTER*80 TMPLN
      CHARACTER*1 KASYM(3,3), ID3
      LOGICAL AUTOSY, NOSYM, DOCART, DOOWN
C
C     Arnfinn, nov -08: Changed all letters to upper case
C     (ie Atomtypes -> ATOMTYPES).
C
      DO I = 1, 71, 10
         WRITE (TMPLN(I:),'(A10)') '          '
      END DO
      IF (NONTYP .LT. 10) THEN
         WRITE (TMPLN,'(A10,I1)') 'ATOMTYPES=', NONTYP
         ISTART = 13
      ELSE IF (NONTYP .LT. 100) THEN
         WRITE (TMPLN,'(A10,I2)') 'ATOMTYPES=', NONTYP
         ISTART = 14
      ELSE
         WRITE (TMPLN,'(A10,I3)') 'ATOMTYPES=', NONTYP
         ISTART = 15
      END IF
C
C     We assume we will not pass +/-9 as a total charge for the molecule
C
      IF (KCHARG .LT. 0) THEN
         WRITE (TMPLN(ISTART:),'(A7,I2)') 'CHARGE=',KCHARG
         ISTART = ISTART + 10
      ELSE IF (KCHARG .GT. 0) THEN
         WRITE (TMPLN(ISTART:),'(A7,I1)') 'CHARGE=',KCHARG
         ISTART = ISTART + 9
      END IF
C
      IF (NOSYM) THEN
         WRITE (TMPLN(ISTART:),'(A11)') 'NOSYMMETRY'
         ISTART = ISTART + 12
Ckr      END IF
C
      ELSE
         IF (.NOT. AUTOSY) THEN
            WRITE (TMPLN(ISTART:),'(A11,I1,A1)')  'GENERATORS=',
     &             NSYMOP,' '
            ISTART = ISTART + 13
            DO I = 1, 3
               DO J = 1, 3
                  WRITE (TMPLN(ISTART:ISTART),'(A1)') KASYM(J,I)
                  ISTART = ISTART + 1
               END DO
            END DO
         END IF 
      END IF
C
      IF (DOCART) THEN
         WRITE (TMPLN(ISTART:),'(A10)') 'CARTESIAN '
         ISTART = ISTART + 10
      END IF
      IF (DOOWN) THEN
         WRITE (TMPLN(ISTART:),'(A10)') 'OWN '
         ISTART = ISTART + 4
      END IF
C
      IF (ID3 .NE. ' ') THEN
         WRITE (TMPLN(ISTART:),'(A7)') 'ANGSTROM'
      END IF
      RETURN
      END
C
      SUBROUTINE LINE5R(TMPLN,Q,NONT,MBSI,IQM,JCO,KANG,BASIS,ATOMBA,
     &                  LMULBS,BASNAM,RADIUS_PCM,ALPHA_PCM)
C
C     Read the atom-specific information in new input style
C
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <pcmlog.h>
      CHARACTER TMPLN*80, BASNAM*80, BKPLIN*80, FRMT*5
      LOGICAL BASIS, ATOMBA, LMULBS, FOUNDR, FOUNDA
      DIMENSION JCO(KANG)
C
      CALL QENTER('LINE5R')

      RADIUS_PCM = 0.0D0
      ALPHA_PCM = 0.0D0

      BKPLIN=TMPLN
      CALL UPCASE(BKPLIN)
      IPOS = INDEX(BKPLIN,'CHA')
      IPOS2 = INDEX(BKPLIN(IPOS:),'=')
      IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 7)) THEN
         WRITE (LUPRI,*) 'Incorrect input for atomic charge in'
         WRITE (LUPRI,*) TMPLN
         WRITE (LUPRI,*) 'Format is "Charge=?"'
         CALL QUIT('Incorrect input for atomic charge')
      ELSE
         READ (BKPLIN((IPOS+IPOS2):),*) Q
         IPOS = INDEX(BKPLIN,'ATO')
         IF (IPOS .NE. 0) THEN
            IPOS2 = INDEX(BKPLIN(IPOS:),'=')
            IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 6)) THEN
               WRITE (LUPRI,*) 'Incorrect input for # of atoms in'
               WRITE (LUPRI,*) TMPLN
               WRITE (LUPRI,*) 'Format is "Atoms=?"'
               CALL QUIT('Incorrect input for # of atoms')
            ELSE
               READ (BKPLIN((IPOS+IPOS2):),*) NONT
            END IF
         END IF
C     
C     Multiple basis sets used
C         
         IPOS = INDEX(BKPLIN,'SET')
         IF (LMULBS) THEN
            IF (IPOS .NE. 0) THEN
               IPOS2 = INDEX(BKPLIN(IPOS:),'=')
               IF (IPOS2 .EQ. 0 .OR. (IPOS2. GT. 5)) THEN
                  WRITE (LUPRI,*) 'Incorrect input for # of '//
     &                 'basis sets in'
                  WRITE (LUPRI,*) TMPLN
                  WRITE (LUPRI,*) 'Format is "Sets=?"'
                  CALL QUIT('Incorrect input for # of '//
     &                 'basis sets')
               ELSE
                  READ (BKPLIN((IPOS+IPOS2):),*) MBSI
               END IF
            END IF
         ELSE IF (IPOS .GT. 0) THEN
            WRITE(LUPRI,*)
     &         'ERROR, multiple basis set, but no .R12AUX specified.'
            WRITE(LUPRI,*) TMPLN
            CALL QUIT('Error when reading .mol file')
         END IF
C     
C     Read in basis set information
C     
         IF (.NOT. (BASIS .OR. ATOMBA)) THEN
            IPOS = INDEX(BKPLIN,'BLO')
            IF (IPOS .NE. 0) THEN
               IPOS2 = INDEX(BKPLIN(IPOS:),'=')
               IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 7)) THEN
                  WRITE (LUPRI,*) 'Incorrect input for # of '//
     &                 'integral blocks in'
                  WRITE (LUPRI,*) TMPLN
                  WRITE (LUPRI,*) 'Format is "Blocks=? ? ?"'
                  CALL QUIT('Incorrect input for # of '//
     &                 'integral blocks')
               ELSE
                  ISTART = IPOS + IPOS2
                  CALL FREFRM(BKPLIN,ISTART,IQM,DUMMY,'INT',IERR)
                  DO IQMLOP = 1, IQM
                     CALL FREFRM(BKPLIN,ISTART,JCO(IQMLOP),DUMMY,'INT',
     &                           IERR)
                  END DO
               END IF
            END IF
         ELSE IF (ATOMBA) THEN
            IPOS = INDEX(BKPLIN,'BAS')
            IF (IPOS .NE. 0) THEN
               IPOS2 = INDEX(BKPLIN(IPOS:),'=')
               IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 6)) THEN
                  WRITE (LUPRI,*) 'Incorrect input for choice '//
     &                 'of atomic basis set in'
                  WRITE (LUPRI,*) TMPLN
                  WRITE (LUPRI,*) 'Format is "Basis=? ? ?"'
                  WRITE (LUPRI,*) 'Input line was:'
                  WRITE (LUPRI,*) BKPLIN
                  CALL QUIT('Incorrect input for choice of '//
     &                 'atomic basis set')
               ELSE
                  IPOS3 = INDEX(BKPLIN((IPOS+IPOS2):),' ')
                  IF (IPOS3 .LE. 1) IPOS3 = LEN( BKPLIN((IPOS+IPOS2):) )
                  IF (IPOS3 .LE. 1) THEN
                     WRITE (LUPRI,*) 'Incorrect input for choice '//
     &                 'of atomic basis set'
                     WRITE (LUPRI,*) 'Format is "Basis=? ? ?"'
                     WRITE (LUPRI,*) 'Input line was:'
                     WRITE (LUPRI,*) BKPLIN
                     CALL QUIT('Incorrect input for choice of '//
     &                 'atomic basis set')
                  END IF
                  IF (IPOS3 .LT. 10) THEN
                     WRITE (FRMT,'(A2,I1,A1,1X)') '(A',IPOS3 - 1,')'
                  ELSE
                     WRITE (FRMT,'(A2,I2,A1)') '(A',(IPOS3 - 1),')'
                  END IF
                  READ (TMPLN((IPOS + IPOS2):),FRMT) BASNAM
               END IF
            ELSE
               WRITE (LUPRI,*) 'ATOMBASIS selected, but no '//
     &           'atomic basis set specified for one atom type in'
               WRITE (LUPRI,*) TMPLN
               CALL QUIT( 'ATOMBASIS selected, but no '//
     &              'atomic basis set specified for one atom type')
            END IF
         END IF
C
C Read in cavity information
C
         FOUNDR = .FALSE.
         FOUNDA = .FALSE.
         RADIUS_PCM = 0.0D0
         ALPHA_PCM = 0.0D0
         IF (PCM) THEN
            IPOS = INDEX(BKPLIN,'RADIUS')
            IF (IPOS .NE. 0) THEN
               IPOS2 = INDEX(BKPLIN(IPOS:),'=')
               IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 7)) THEN
                  WRITE (LUPRI,*) 'Incorrect input for PCM radius'
                  WRITE (LUPRI,*) 'Format is "Radius=?"'
                  CALL QUIT('Incorrect input for PCM sphere radius')
               ELSE
                  READ (BKPLIN((IPOS+IPOS2):),*) RADIUS_PCM
                  FOUNDR = .TRUE.
               END IF
            END IF
            IPOS = INDEX(BKPLIN,'ALPHA')
            IF ((IPOS .NE. 0) .AND. FOUNDR) THEN
               IPOS2 = INDEX(BKPLIN(IPOS:),'=')
               IF (IPOS2 .EQ. 0 .OR. (IPOS2 .GT. 7)) THEN
                  WRITE (LUPRI,*) 'Incorrect input for PCM alpha'
                  WRITE (LUPRI,*) 'Format is "Alpha=?"'
                  CALL QUIT('Incorrect input for PCM sphere alpha')
               ELSE
                  READ (BKPLIN((IPOS+IPOS2):),*) ALPHA_PCM
                  FOUNDA = .TRUE.
               END IF
            END IF
            IF(FOUNDR .AND. .NOT. FOUNDA) ALPHA_PCM = 1.0
         END IF
      END IF
      CALL QEXIT('LINE5R')
      RETURN
      END
C
      SUBROUTINE LINE5W(LINE,Q,NONT,MBSI,BASIS,ATOMBA,LMULBS,
     &                  BASNAM,IQM,JCO,KANG,RADIUS_PCM,ALPHA_PCM)
C
C     Write out atom-specific input line in new format
C     
#include <implicit.h>
#include <pcmlog.h>
      CHARACTER LINE*80, BASNAM*80
      DIMENSION JCO(KANG)
      LOGICAL   BASIS, ATOMBA, LMULBS
C
      DO I = 1, 71, 10
         WRITE (LINE(I:),'(A10)') '          '
      END DO
C
      ISTART = 17
      IF (ABS(Q) .LT. 1.0D1) THEN
         WRITE (LINE,'(A7,F8.5)') 'Charge=', Q
      ELSE IF (ABS(Q) .LT. 1.0D2) THEN
         WRITE (LINE,'(A7,F8.4)') 'Charge=', Q
      ELSE IF (ABS(Q) .LT. 1.0D4) THEN
         WRITE (LINE,'(A7,F8.2)') 'Charge=', Q
      ELSE IF (ABS(Q) .LT. 1.0D6) THEN
         WRITE (LINE,'(A7,F8.0)') 'Charge=', Q
      ELSE
         WRITE (LINE,'(A7,1P,D13.6)') 'Charge=', Q
         ISTART = 22
      END IF
C
      IF (NONT .LT. 10) THEN
         WRITE (LINE(ISTART:),'(A6,I1)') 'Atoms=', NONT
         ISTART = ISTART + 8
      ELSE IF (NONT .LT. 100) THEN
         WRITE (LINE(ISTART:),'(A6,I2)') 'Atoms=', NONT
         ISTART = ISTART + 9
      ELSE
         WRITE (LINE(ISTART:),'(A6,I3)') 'Atoms=', NONT
         ISTART = ISTART + 10
      END IF
C
      IF (LMULBS) THEN
         WRITE (LINE(ISTART:),'(A5,I1)') 'Sets=', MBSI
         ISTART = ISTART + 7
      END IF
C
      IF (.NOT. (BASIS .OR. ATOMBA)) THEN
         WRITE (LINE(ISTART:),'(A7,I1)') 'Blocks=',IQM
         ISTART = ISTART + 9
         DO I = 1, IQM
            WRITE (LINE(ISTART:),'(I4)') JCO(I)
            ISTART = ISTART + 4
         END DO
      END IF
C
      IF (PCM.AND.(RADIUS_PCM * ALPHA_PCM .GT. 0.1)) THEN  
         WRITE (LINE(ISTART:),'(A7,F6.2)') 'Radius=',RADIUS_PCM
         WRITE (LINE(ISTART:),'(A7,F6.2)') 'Alpha=',ALPHA_PCM
         ISTART = ISTART + 30
      END IF
c
      IF (ATOMBA) THEN
         WRITE (LINE(ISTART:),'(A6,A25)') 'Basis=',BASNAM(1:25)
         ISTART = ISTART + 31
      END IF
C
      RETURN
      END

C --- end of herrdn.F ---
