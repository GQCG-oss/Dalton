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
#ifdef REV_LOG
c970923-vb
cImproved symmetry classification of normal modes through procedure DETIRP.
cAdded VRML visualization of normal modes.
#endif
C  /* Deck vibinp */
      SUBROUTINE VIBINP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (NTABLE = 7)
#include <abainf.h>
#include <cbivib.h>
C
      LOGICAL NEWDEF, DOINT
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, CARD*80, TYPE*4,
     &          WORD1*7
      DIMENSION IA(4)
      DATA TABLE /'xxxxxxx', '.PRINT ', '.ISOTOP', '.INTERN',
     &            '.XXXXXX', '.HESFIL', '.HESPUN'/
      DATA ISOTPD/0/, ISOTHD/0/, DOINT/.FALSE./
C
      NEWDEF = ((WORD .EQ. '*VIBANA').OR.(WORD.EQ.'*HARMON'))
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
  101       PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     &            '" not recognized in VIBINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in VIBINP.')
    1          CONTINUE
               GO TO 100
    2          CONTINUE
                  READ (LUCMD,*) IPRINT
                  IF (IPRINT .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    3          CONTINUE
                  READ (LUCMD,*) NISOTP, NATM
                  IF (NISOTP .GT. MAXSUB) THEN
                     WRITE (LUPRI,*)' Too many isotopic substitutions'//
     &                    ' requested in *VIBANA input'
                     WRITE (LUPRI,*) ' Increase MAXSUB, recompile and'//
     &                    ' be welcome back'
                     CALL QUIT('Input error in *VIBANA')
                  END IF
                  DO 305 ICOUNT = 1, NISOTP
                     READ (LUCMD,*) (ISOTP(ICOUNT,N), N = 1, NATM)
 305              CONTINUE
               GO TO 100
    4             CONTINUE
  400                DOINT = .TRUE.
                     READ (LUCMD,'(A)') CARD
                     READ (CARD,'(1X,A4)') TYPE
                     NCARD = NCARD + 1
                     IF (TYPE .NE. 'STRE' .AND.
     &                   TYPE .NE. 'INVR' .AND.
     &                   TYPE .NE. 'BEND' .AND.
     &                   TYPE .NE. 'OUT ' .AND.
     &                   TYPE .NE. 'TORS' .AND.
     &                   TYPE .NE. 'LIN1' .AND.
     &                   TYPE .NE. 'LIN2' .AND.
     &                   TYPE .NE. '    ' ) THEN
                         KWORD(NCARD) = '*'
                         READ (CARD,'(A7)') WORD
                         GO TO 101
                     END IF
                     READ (CARD,'(1X,A4,4I5,2F10.5)')TYPE,IA,CCOEF,SCALI
                     IF (TYPE .EQ. '    ') THEN
                        KWORD(NCARD) = ' '
                     ELSE
                        KWORD(NCARD) = 'K'
                        NINTCM = NINTCM + 1
                        READ (TYPE,'(A4)') ITYPCM(NINTCM)
                     END IF
                     READ (TYPE,'(A4)') ITYPCD(NCARD)
                     SCALE(NCARD)    = SCALI
                     COEF(NCARD)     = CCOEF
                     IATOMS(1,NCARD) = IA(1)
                     IATOMS(2,NCARD) = IA(2)
                     IATOMS(3,NCARD) = IA(3)
                     IATOMS(4,NCARD) = IA(4)
                  GO TO 400
 5             CONTINUE
               GO TO 100
 6             CONTINUE
                  HESFIL = .TRUE.
c                  MOLHES = .FALSE.
c                  MOLGRD = .FALSE.
               GOTO 100
 7             CONTINUE
                  HESPUN = .TRUE.
               GOTO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     &            '" not recognized in VIBINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in VIBINP.')
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for VIBANA:',0)
         IF (IPRINT .NE. IPRDEF) THEN
            WRITE (LUPRI,'(A,I5)') ' Print level in VIBANA:',IPRINT
         END IF
         IF (DOINT) THEN
            WRITE (LUPRI,'(A)') ' Forces and force constants will '//
     &          'be printed in internal coordinates.'
         END IF
         IF (HESFIL) THEN
            WRITE (LUPRI,'(A)') ' Vibrational analysis will be '//
     &           'based on Hessian read from file'
         END IF
         IF (HESPUN) THEN
            WRITE(LUPRI,'(A)') ' The vibrational Hessian will be '//
     &           'punched to file'
         END IF
         IF (NISOTP .EQ. 0) THEN
            WRITE (LUPRI, '(A)')
     &           ' Vibrational analysis for parent molecule only.'
         ELSE
            WRITE (LUPRI,'(4X,A,I3,A)')
     &           ' Vibrational analysis for ',NISOTP+1,' molecules'
         END IF
      END IF
      RETURN
      END
C  /* Deck vibini */
      SUBROUTINE VIBINI
C
C     Initialize /CBIVIB/
C
#include <implicit.h>
#include <mxcent.h>
#include <abainf.h>
#include <cbivib.h>
C
      NCARD  = 0
      NINTCM = 0
      DIPOL  = DIPDER
      IPRINT = IPRDEF
      NISOTP = 0
      HESFIL = LINCPL
      HESPUN = .FALSE.
      DO 10 I = 1, MAXSUB
         DO 10 J = 1, MXCENT
         ISOTP(I,J) = 1
 10   CONTINUE
      IF (MOLHES) THEN
         MAXDIF = 2
      ELSE IF (MOLGRD) THEN
         MAXDIF = 1
      END IF
      DOVCD = VCD
      RETURN
      END
C  /* Deck vibctl */
      SUBROUTINE VIBCTL(WORK,LWORK)
C
C     June 1985 PJ
C     December 1985 tuh - internal coordinates
C     September 1989 tuh - symmetry
C
C     Driver for vibrational analysis. This subroutine calls
C
C       1) GETNOS for transformation to non-symmetry basis
C       2) INTERN for internal coordinate analysis
C       3) VIBANA for vibrational analysis
C
C     Note: the vibrational analysis is always carried out in
C     non-symmetry Cartesian coordinates, regardless of whether
C     the gradient and Hessian have been calculated in symmetry
C     coordinates. Therefore GETNOS is called first to set up the
C     geometry, gradient and Hessian in non-symmetry coordinates.
C
C     Internal coordinates are used for two purposes:
C
C       1) for printing the gradient and Hessian
C       2) for correcting the Hessian for non-equilibrium forces
C
C     There is no need to define a set of internal coordinates for
C     the vibrational analysis, but at non-equilibrium geometries the
C     calculated frequencies are more meaningful if internal coordinates
C     are defined. At equilibrium the frequencies are the same with or
C     without internal coordinates defined.
C
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <iratdef.h>
      PARAMETER ( THRGEO = 1.D-4, THRSH = 1.0D-6)
      LOGICAL HESEXS
C
      DIMENSION WORK(LWORK)
C
C Used from common blocks:
C   TAYMOL : HESMOL
C
#include <abainf.h>
#include <taymol.h>
#include <nuclei.h>
#include <symmet.h>
#include <dorps.h>
#include <cbivib.h>
C
      CALL QENTER('VIBCTL')
      LUHES = -1
      CALL TIMER('START ',TIMSTR,TIMEND)
      CALL TITLER
     &   ('ABACUS - ROTATIONAL AND VIBRATIONAL ANALYSIS','*',115)
C
      NCORD  = 3*NATOMS
      NINTER = NINTCM
C
C     ********************************************
C     ***** Transform charges, gradient, and *****
C     ***** Hessian to non-symmetry basis    *****
C     ********************************************
C
      KCHRG  = 1
      KGEOM  = KCHRG  + (NATOMS + IRAT - 1)/IRAT
      KGRAD  = KGEOM  + NCORD
      KHESS  = KGRAD  + NCORD
      KGRDN  = KHESS  + NCORD*NCORD
      KDGRAD = KGRDN  + NCORD
      KCAAT  = KDGRAD + 3*NCORD
      KWRK1  = KCAAT  + 3*NCORD
      LWRK1  = LWORK  - KWRK1
      IF (KWRK1 .GE. LWORK) CALL STOPIT('VIBCTL','GETNOS',KWRK1,LWORK)
      CALL GETNOS(WORK(KCHRG),WORK(KGEOM),WORK(KGRAD),WORK(KHESS),
     &            WORK(KGRDN),WORK(KDGRAD),WORK(KCAAT),
     &            WORK(KWRK1),LWRK1,NATOMS,NCORD,DIPDER,DOVCD,IPRINT,
     &            HESPUN)
      IF (NUMHES) THEN
         KCOOR0 = KWRK1
         KGRAD0 = KCOOR0 + 3*NUCIND
         KHESS0 = KGRAD0 + NCORD
         KGRADF = KHESS0 + NCORD*NCORD
C
C     We use the gradient from the first iteration, since we are not in
C     the minimum geometry when we end geometry steps using NUMHES
C
ckr         CALL REANMHES(WORK(KCOOR0),WORK(KGRAD0),WORK(KHESS0),
         CALL REANMHES(WORK(KGEOM),WORK(KGRAD),WORK(KHESS0),
     &                 WORK(KGRADF),WORK(KHESS),NCORD,3*NUCIND)
C
C        Note that the Hessian is not back-transformed to non-symmetry
C        basis, since it is assumed that the VROA and/or Raman calculation
C        is anyway done without the use of symmetry.\kr\
C
      END IF
C
C     If vibrational analysis is based on Hessian read from file,
C     we overwrite the Hessian information, K.Ruud and G.Hangartner, Oct.-96
C
      IF (HESFIL) THEN
         CALL DZERO(WORK(KHESS),NCORD*NCORD)
         CALL GPINQ('DALTON.HES','EXIST',HESEXS)
         IF (.NOT. HESEXS)
     &        CALL QUIT('Unable to open the file DALTON.HES.')
         CALL GPOPEN(LUHES,'DALTON.HES','OLD',' ','FORMATTED',
     &              IDUMMY,.FALSE.)
         READ(LUHES,*) IDIM
         READ(LUHES,*)
         IF (IDIM .NE. NCORD) CALL QUIT('The Hessian in DALTON.HES '//
     &        'has wrong dimensions.')
         IPOS = KHESS
         DO 10 I = 1, NCORD
            DO 15 J = 1, NCORD
               READ(LUHES,*) WORK(IPOS)
               IPOS = IPOS + 1
 15         CONTINUE
            READ(LUHES,*)
 10      CONTINUE
C
C     We try to see if the user has supplied some coordinates on the
C     file that we should use for updating the coordinate information
C     Also check if they match coordinates from MOLECULE.INP
C
         IDIFF = 0
         DO 16 I = 0, NCORD - 1
            READ (LUHES,*,END=18) TCOOR
            IF (ABS(TCOOR - WORK(KGEOM + I)) .GT. THRSH) THEN
               IDIFF = IDIFF + 1
               WORK(KGEOM + I) = TCOOR
            END IF
 16      CONTINUE
         IF (IDIFF .GT. 0) THEN
            WRITE (LUPRI,'(/,A,/,12X,A)') ' WARNING>>> Coordinates on'//
     &           ' DALTON.HES differ from MOLECULE.INP',
     &           'Coordinates from DALTON.HES will be used'
         END IF
 18      CONTINUE
         CALL GPCLOSE(LUHES,'KEEP')
      END IF
C
C     ****************************************
C     ***** INTERNAL COORDINATE ANALYSIS *****
C     ****************************************
C
C     INTERN transforms gradient and Hessian to internal coordinates.
C     Also, a correction for non-equilibrium geometries is calculated
C     for the Hessian. In the remaining part of the program the
C     corrected Hessian in non-symmetry basis is used.
C
      IF (NINTER .GT. 0) THEN
         KHESCR = KWRK1
         KBMAT  = KHESCR + NCORD*NCORD
         KBADJ  = KBMAT  + NINTER*NCORD
         KQMAT  = KBADJ  + NINTER*NCORD
         KQINV  = KQMAT  + NINTER*NINTER
         KBHES  = KQINV  + NINTER*NINTER
         KFINT  = KBHES  + NINTER*NCORD
         KDINT  = KFINT  + NINTER
         KAATIN = KDINT  + NINTER*3
         KHINT  = KAATIN + NINTER*3
         KW2    = KHINT  + NINTER*NINTER
         KIWRK  = KW2    + NINTER
         KCMAT  = KIWRK  + NINTER
         KQQ    = KCMAT  + NINTER*NCORD*NCORD
         KTYPE  = KQQ    + NINTER
         KWRK1  = KTYPE  + NINTER
         LWRK1  = LWORK  - KWRK1
         IF (KWRK1.GE.LWORK) CALL STOPIT('VIBCTL','INTERN',KWRK1,LWORK)
C
C        Only output from this routine is WORK(KHESCR) which contains
C        the Hessian corrected for non-equilibrium effects
C
         CALL INTERN(NCORD,NINTER,WORK(KGRAD),WORK(KHESS),WORK(KHESCR),
     &               WORK(KDGRAD),WORK(KCAAT),WORK(KBMAT),WORK(KBADJ),
     &               WORK(KQMAT),WORK(KQINV),WORK(KBHES),WORK(KFINT),
     &               WORK(KDINT),WORK(KAATIN),WORK(KHINT),WORK(KW2),
     &               WORK(KIWRK),WORK(KCMAT),WORK(KQQ),WORK(KTYPE),
     &               WORK(KWRK1),LWRK1)
      END IF
C
C     ********************************
C     ***** VIBRATIONAL ANALYSIS *****
C     ********************************
C
      IF (MAXDIF .EQ. 2 .OR. HESFIL .OR. NUMHES) THEN
C
C        ***** Check if minimum geometry *****
C
         GRDNRM = DDOT(NCORD,WORK(KGRAD),1,WORK(KGRAD),1)
         GRDNRM = SQRT(GRDNRM)
         IF (DOREPS(0)) THEN
            WRITE(LUPRI,'(/A,F15.8)')
     &        ' Norm of (unprojected) molecular gradient (au) :',GRDNRM
         END IF
         IF (DOREPS(0) .AND. GRDNRM.GT.THRGEO .AND. NFLOAT.GT.0) THEN
               WRITE (LUPRI,'(/2A,2(/A),/)')
     &          ' Vibrational analysis cannot be performed ',
     &          ' at non-equilibrium geometry when',
     &          ' floating orbitals are used.'
         ELSE
            IF (DOREPS(0) .AND. GRDNRM.GT.THRGEO) THEN
               IF (NINTER .GT. 0) THEN
                 WRITE (LUPRI,'(/3(2A/))')
     &           ' Vibrational analysis performed at a',
     &           ' non-equilibrium geometry,',
     &           ' Hessian has been corrected for rotation coordinates',
     &           ' using procedure described by Peter Pulay',
     &           ' in Applications of Electronic Structure',
     &           ' Theory, ed. H.F.Schaefer, Plenum 1977, p. 165-167)'
               ELSE
                 WRITE (LUPRI,'(/2A,2(/A),/)')
     &           ' WARNING: Vibrational analysis performed at',
     &           ' a non-equilibrium geometry,',
     &           ' better results may be obtained using .INTERN',
     &          ' (this will correct Hessian for rotation coordinates).'
               END IF
            END IF
C
            KAMASS = KWRK1  + 3*NCORD
            KEVAL  = KAMASS + NATOMS
            KEVEC  = KEVAL  + NCORD
            KNUMIS = KEVEC  + NCORD*NCORD
            KMAXIS = KNUMIS + (NATOMS + IRAT - 1)/IRAT
            KWRK1  = KMAXIS + (NATOMS + IRAT - 1)/IRAT
            LWRK1  = LWORK  - KWRK1
            IF(KWRK1.GE.LWORK)CALL STOPIT('VIBCTL','ISOMOL',KWRK1,LWORK)
            CALL ISOMOL(WORK(KAMASS),WORK(KGEOM),
     &                  WORK(KGRAD),WORK(KHESS),
     &                  WORK(KGRDN),WORK(KEVAL),WORK(KEVEC),
     &                  WORK(KDGRAD),WORK(KCAAT),WORK(KWRK1),LWRK1,
     &                  WORK(KNUMIS),WORK(KMAXIS),WORK(KCHRG),NATOMS,
     &                  NCORD,DIPDER,POLAR,VROA.OR.RAMAN)
         END IF
      END IF
      CALL TIMER('VIBCTL',TIMSTR,TIMEND)
      CALL QEXIT('VIBCTL')
      RETURN
      END
C  /* Deck isomol */
      SUBROUTINE ISOMOL(AMASS,GEOM,GRAD,HESS,GRDN,EVAL,EVEC,DGRAD,
     &                  CAAT,WRK,LWRK,NUMIS,MAXIS,NATTYP,NATOMS,NCORD,
     &                  DIPDER,POLAR,RAMPRP)
C
C    This subroutine loops over isotopic molecules, calling routines for
C
C      1) rotational analysis (ROTANA)
C      2) vibrational analysis (VIBANA)
C      3) partition functions  (PRTFUN)
C      4) vibrational polarizabilities (VIBPOL)
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <codata.h>
#include <iratdef.h>
C
      LOGICAL DIPDER, POLAR, LINEAR, NOTALL, NONE, RAMPRP
      DIMENSION NUMIS(NATOMS), MAXIS(NATOMS), NATTYP(NATOMS),
     &          AMASS(NATOMS), GEOM(NCORD),
     &          GRAD(NCORD), HESS(NCORD,NCORD),
     &          GRDN(NCORD),
     &          EVAL(NCORD), EVEC(NCORD,NCORD), DGRAD(3,NCORD),
     &          VINRTI(3,3), CAAT(3,NCORD), WRK(LWRK)
C
#include <gnrinf.h>
#include <cbivib.h>
#include <moldip.h>
#include <symmet.h>
#include <dorps.h>
#include <trkoor.h>
#include <cbilnr.h>
#include <cbisol.h>
C
      DIMENSION NMASS(36)
      DATA NMASS /3,2,2,1,2,  2,2,3,1,3,
     &            1,3,1,3,1,  4,2,3,3,5,
     &            1,5,2,4,1,  4,1,5,2,5,
     &            2,5,1,5,2,  5/
C
C     *************************
C     ***** Print Section *****
C     *************************
C
      WRITE (LUPRI,'(//A//,A,F12.5,/A,F12.5,/A,F12.9,
     &                    /A,F12.6)')
     &   ' Conversion factors used:',
     &   '    Hartree to (cm)-1      : ',XTKAYS,
     &   '    a.m.u.  to a.u.        : ',XFAMU,
     &   '    Bohr    to Angstrom    : ',XTANG,
     &   '    a.u.    to A**4amu-1   : ',XTANG**4*XFAMU
      IF (DIPDER) WRITE (LUPRI,'(A,F12.3)')
     &   '    a.u.    to km (mol)-1  : ',XTKMML
      IF (NISOTP .EQ. 0) THEN
         WRITE (LUPRI, '(//A)')
     &    ' Vibrational analysis for parent molecule only.'
      ELSE
         WRITE (LUPRI,'(//3X,A,I3,2X,A)')
     &    ' Vibrational analysis for ',NISOTP + 1,'molecules'
      END IF
C
C     ******************************************************************
C     **** Project out trans-rot coordinates from molecular Hessian ****
C     ******************************************************************
C
      KTRROV = 1
      KWRK1  = KTRROV + 6*NCORD
      IF (SOLVNT) KWRK1 = KTRROV + 9*NCORD
      IF (KWRK1 .GE. LWRK) CALL STOPIT('ISOMOL','VIBHES',KWRK1,LWRK)
      CALL VIBHES(IPRINT,NCORD,GEOM,KTRRO,GRAD,HESS,WRK(KTRROV),
     &            NATTYP,WRK(KWRK1))
      IF (NINTCM .GT. 0 .AND. KTRRO .NE. NCORD - NINTCM) THEN
         WRITE (LUPRI,'(/A,I2,A,/,A,I2)') ' Number of defined'//
     &     ' internal coordinates',NINTCM,' is not consistent',
     &     ' with the number of independent external coordinates',KTRRO
         CALL QUIT
     &       ('ABACUS.VIBANA: Inconsistency in number of coordinates')
      END IF
C
C     *************************************
C     ***** Isotope selection vectors *****
C     *************************************
C
      NTOT = NISOTP
      DO 100 I = 1, NATOMS
         NUMIS(I)  = 1
         MAXIS(I)  = NMASS(NATTYP(I))
 100  CONTINUE
C
C     ****************************************************************
C     ***** Check that all vibrational symmetries are calculated *****
C     ****************************************************************
C
      NOTALL = .FALSE.
      NONE   = .TRUE.
      DO 150 IREP = 0, MAXREP
         IF (NCRREP(IREP,1) .GT. NTRREP(IREP)) THEN
            IF (.NOT.DOREPS(IREP)) THEN
               NOTALL = .TRUE.
            ELSE
               NONE = .FALSE.
            END IF
         END IF
  150 CONTINUE
      IF (OPTWLK .AND. NONE) THEN
         NONE = .FALSE.
         NOTALL = .FALSE.
      END IF
      IF (NONE) THEN
         WRITE (LUPRI,'(2(/A))')
     &      ' Note: No vibrational analysis is carried out since the'//
     &      ' Cartesian derivatives',
     &      ' of the vibrational symmetries are not available.'
      ELSE IF (NOTALL) THEN
         NTOT = 1
         WRITE (LUPRI,'(3(/A))')
     &      ' Note: A full vibrational analysis is not possible since'//
     &      ' the Cartesian',
     &      ' derivatives of all vibrational symmetries are not'//
     &      ' available.',
     &      ' Also no isotope substituted molecules are analyzed.'
      END IF
C
C
C     ****************************************
C     ***** Loop over isotopic molecules *****
C     ****************************************
C
      ISPEC = 0
      KORIG = 1
      KLAST = KORIG + 3
      IF (KLAST .GE. LWRK) CALL STOPIT('ISOMOL','KMOMS',KLAST,LWRK)
      DO 200 I = 0, NISOTP
         DO 210 IATOM = 1, NATOMS
            IF (I .EQ. 0) THEN
               NUMIS(IATOM) = 1
            ELSE
               NUMIS(IATOM) = ISOTP(I,IATOM)
            END IF
 210     CONTINUE
C
         ISPEC = ISPEC + 1
         WRITE (LUPRI,'(//24X,A,I2)') ' Isotopic Molecule No.',ISPEC
         WRITE (LUPRI,'(  24X,A   )') ' ======================= '
C
C           ****************************************
C           ***** Atomic masses and total mass *****
C           ****************************************
C
C           Output: AMASS  - atomic masses
C                   TOTMAS - molecular mass
C
         CALL VIBMAS(AMASS,TOTMAS,NUMIS,NATTYP,NATOMS,GEOM,
     &               WRK(KORIG),MAX(IPRINT,1))
C
C           *******************************
C           ***** Rotational analysis *****
C           *******************************
C
C           Output: VINRTI              - principal axes of inertia
C                   AIMOM, BIMOM, CIMOM - principal moments of inertia
C                   LINEAR              - true for linear molecules
C
         CALL ROTANA(GEOM,AMASS,VINRTI,WRK(KORIG),
     &               AIMOM,BIMOM,CIMOM,LINEAR,
     &               NATOMS,NCORD,IPRINT)
C
         IF (HESFIL) NONE = .FALSE.
         IF (.NOT.NONE) THEN
C
C              ********************************
C              ***** Vibrational analysis *****
C              ********************************
C
C              Output: EVAL   - eigenvalues of mass-weighted Hessian
C                      EVEC   - normal coordinates
C                      NUMMOD - number of vibrational modes
C
            KNSPEC = KLAST  + 3*NTOT
            KIFRQ  = KNSPEC + (NCORD + IRAT - 1)/IRAT
            KEVECS = KIFRQ  + (NCORD + IRAT - 1)/IRAT
            KFREQ  = KEVECS + NCORD*NCORD
            KDINTX = KFREQ  + NCORD
            KDINTY = KDINTX + NCORD
            KDINTZ = KDINTY + NCORD
            KDINTS = KDINTZ + NCORD
            KAATX  = KDINTS + NCORD
            KAATY  = KAATX  + NCORD
            KAATZ  = KAATY  + NCORD
            KAATS  = KAATZ  + NCORD
            KROTST = KAATS  + NCORD
            KIMAG  = KROTST + NCORD
            KIRP   = KIMAG  + NCORD
            KLAST  = KIRP  + (NCORD + IRAT - 1)/IRAT
C
            IF (RAMPRP) THEN
               KROAAF = KLAST
               KROAGN = KROAAF + 9*NFRVAL*NCORD
               KROAGL = KROAGN + 9*NFRVAL*NCORD
               KROAA  = KROAGL + 9*NFRVAL*NCORD
               KLAST  = KROAA  +27*NFRVAL*NCORD
            ELSE
               KROAAF = KLAST
               KROAGN = KLAST
               KROAGL = KLAST
               KROAA  = KLAST
            END IF
            KWRK1  = KLAST
            LWRK1  = LWRK - KWRK1 + 1
            IF (KWRK1 .GE. LWRK) CALL STOPIT('ISOMOL','VIBANA',KWRK1,
     &                                       LWRK)
            CALL VIBANA(AMASS,HESS,GRAD,GRDN,EVAL,EVEC,
     &                  WRK(KEVECS),VINRTI,WRK(KFREQ),DGRAD,CAAT,
     &                  WRK(KDINTX),WRK(KDINTY),WRK(KDINTZ),
     &                  WRK(KDINTS),WRK(KAATX),WRK(KAATY),WRK(KAATZ),
     &                  WRK(KAATS),WRK(KROTST),WRK(KWRK1),LWRK1,
     &                  WRK(KNSPEC),WRK(KIFRQ),NUMMOD,
     &                  WRK(KIMAG),WRK(KIRP),DIPDER,
     &                  WRK(KROAAF),WRK(KROAGN),WRK(KROAGL),
     &                  WRK(KROAA),WRK(KWRK1),LWRK1,NCORD,NOTALL,
     &                  DOVCD,LINEAR,IPRINT)
C
C              *******************************
C              ***** Partition functions *****
C              *******************************
C
            IF (.NOT.NOTALL) THEN
               CALL PRTFUN(TOTMAS,AIMOM,BIMOM,CIMOM,EVAL,
     &                     LINEAR,NCORD,NUMMOD)
            END IF
C
C              ****************************************
C              ***** Vibrational polarizabilities *****
C              ****************************************
C
            IF (DIPDER) THEN
               CALL VIBPOL(DGRAD,EVAL,EVEC,VINRTI,POLAR,NCORD,NUMMOD)
            END IF
         END IF
C
C        Next isotopic molecule
C
 200  CONTINUE
C
C     End of ISOMOL
C
      RETURN
      END
C  /* Deck rotana */
      SUBROUTINE ROTANA(GEOM,AMASS,VINRTI,VCMNUC,
     &                  AIMOM,BIMOM,CIMOM,
     &                  LINEAR,NATOMS,NCORD,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <codata.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (D0 = 0.0D0, TSTLIN =1.D-4, DSMALL = 1.D-10, D2 = 2.0D0)
      PARAMETER (ROTMHZ = (XTHZ*1.0D-6*(XTANG**2))/(D2*XFAMU),
     &           CINVRS = 1.0D4/CCM)
C
      DIMENSION GEOM(NCORD), AMASS(NATOMS),
     &          AINRTI(3,3), AINRTP(6), VINRTI(3,3), VCMNUC(3),
     &          AWRKJC(3), IWRKJC(3)
C
      LOGICAL  LINEAR, PLANAR
#include <cbisol.h>
#include <symmet.h>
#include <inirep.h>
C
C     ***** Moments of inertia *****
C
      CALL DZERO(AINRTI,9)
      KK = 1
      KMAX = NATOMS
      IF (SOLVNT) KMAX = NATOMS - 1
      DO 200 K = 1, KMAX
C
         RMASS       = AMASS(K)
         XK          = XTANG*(GEOM(KK)   - VCMNUC(1))
         YK          = XTANG*(GEOM(KK+1) - VCMNUC(2))
         ZK          = XTANG*(GEOM(KK+2) - VCMNUC(3))
C
         AINRTI(1,1) = AINRTI(1,1) + RMASS*(YK*YK + ZK*ZK)
         AINRTI(2,1) = AINRTI(2,1) - RMASS* XK*YK
         AINRTI(3,1) = AINRTI(3,1) - RMASS* XK*ZK
         AINRTI(2,2) = AINRTI(2,2) + RMASS*(XK*XK + ZK*ZK)
         AINRTI(3,2) = AINRTI(3,2) - RMASS* YK*ZK
         AINRTI(3,3) = AINRTI(3,3) + RMASS*(XK*XK + YK*YK)
C
         KK = KK + 3
  200 CONTINUE
      AINRTI(1,2) = AINRTI(2,1)
      AINRTI(1,3) = AINRTI(3,1)
      AINRTI(2,3) = AINRTI(3,2)
      IF (IPRINT .GE. 4) WRITE (LUPRI,'(//A,3(/5X,3F12.6))')
     &   ' Moments of inertia (u*A**2) :',
     &   ((AINRTI(K,L),L=1,3),K=1,3)
C
C     ***** Principal moments of inertia *****
C
      AINRTP(1) = AINRTI(1,1)
      AINRTP(2) = AINRTI(2,1)
      AINRTP(3) = AINRTI(2,2)
      AINRTP(4) = AINRTI(3,1)
      AINRTP(5) = AINRTI(3,2)
      AINRTP(6) = AINRTI(3,3)
      CALL DUNIT(VINRTI,3)
      CALL JACO(AINRTP,VINRTI,3,3,3,AWRKJC,IWRKJC)
      AINRTP(2) = AINRTP(3)
      AINRTP(3) = AINRTP(6)
      CALL ORDER(VINRTI,AINRTP,3,3)
      AIMOM = AINRTP(1)
      BIMOM = AINRTP(2)
      CIMOM = AINRTP(3)
C
      CALL HEADER
     &    ('Principal moments of inertia (u*A**2) and principal axes',0)
      WRITE (LUPRI,'(3X,A,F15.6,6X,3F12.6)')
     &     'IA',AIMOM,(VINRTI(K,1),K=1,3),
     &     'IB',BIMOM,(VINRTI(K,2),K=1,3),
     &     'IC',CIMOM,(VINRTI(K,3),K=1,3)
C
C        ***** Symmetries of principal axes *****
C
      DO 400 IPAX = 1, 3
         NREPPI(IPAX) = 0
         DO 410 IREP = 0, MAXREP
            DO 420 ICOOR = 1, 3
               IF (ISYMAX(ICOOR,1) .EQ. IREP) THEN
                  IF (ABS(VINRTI(ICOOR,IPAX)) .GT. DSMALL) THEN
                     NREPPI(IPAX) = NREPPI(IPAX) + 1
                     IREPPI(IPAX,NREPPI(IPAX)) = IREP
                     GO TO 410
                  END IF
               END IF
 420        CONTINUE
 410     CONTINUE
 400  CONTINUE
C
C        ***** Rotational constants *****
C
      CALL HEADER('Rotational constants',0)
      IF (CIMOM .LT. TSTLIN) THEN
C        ... this is an atom !
         WRITE (LUPRI,'(A/)')
     &      ' All are zero as the "molecule" is an atom!'
         LINEAR = .TRUE.
         GO TO 9999
      END IF
      IF ( ABS(CIMOM-BIMOM-AIMOM) .LT. TSTLIN) THEN
         PLANAR = .TRUE.
      ELSE
         PLANAR = .FALSE.
      END IF
      IF (AIMOM .LT. TSTLIN) THEN
         LINEAR = .TRUE.
         ROTA   = D0
      ELSE
         LINEAR = .FALSE.
         ROTA   = ROTMHZ/AIMOM
      END IF
      ROTB = ROTMHZ/BIMOM
      ROTC = ROTMHZ/CIMOM
      IF (LINEAR) THEN
         WRITE (LUPRI,'(A/)') ' The molecule is linear.'
         WRITE (LUPRI,'(T16,A,F16.2,A,F12.6,A)')
     &        'B =', ROTB, ' MHz     (', CINVRS*ROTB, ' cm-1)'
      ELSE IF (PLANAR) THEN
         WRITE (LUPRI,'(A/)') ' The molecule is planar.'
      END IF
      IF (.NOT. LINEAR)
     &     WRITE (LUPRI,'(T16,A,T36,A,T56,A,//3F20.4,A,/3F20.6,A)')
     &     'A','B','C', ROTA, ROTB, ROTC, ' MHz',
     &     CINVRS*ROTA, CINVRS*ROTB, CINVRS*ROTC, ' cm-1'
 9999 RETURN
      END
C  /* Deck vibana */
      SUBROUTINE VIBANA(AMASS,HESS,GRAD,GRDN,EVAL,EVEC,EVECS,VINRTI,
     &                  FREQAU,DGRAD,CAAT,DINTX,DINTY,DINTZ,DINTNS,
     &                  AATX,AATY,AATZ,AATXYZ,ROTSTR,WRK,LWRK,
     &                  NSPEC,IFRQCM,NUMMOD,
     &                  IMAGIN,IREPS,INTENS,ROAAFQ,ROAGNQ,ROAGLQ,ROAAQ,
     &                  WORK,LWORK,NCORD,NOTALL,DOVCD,LINEAR,IPRINT)
#include <implicit.h>
#include <dummy.h>
C
C     Carry out a propertry analysis including finding normal
C     frequencies, normal coordinates and IR intensities
C
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <codata.h>
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D2 = 2.0D0, THRSH = 1.0D-8)
      PARAMETER (ESUDIP = 1.D4*DEBYE**2,
     &           ESUROT =  ECHARGE*XTANG*CCM*1.D40*ECHARGE*HBAR/EMASS)
C
      LOGICAL DONE, SAME, NOTALL, DOVCD
      DIMENSION HESS(NCORD,NCORD), GRAD(NCORD), GRDN(NCORD),
     &          AMASS(*), EVAL(NCORD), EVEC(NCORD,NCORD),
     &          EVECS(NCORD,NCORD), FREQAU(NCORD), WRK(LWRK),
     &          DINTX(NCORD), DINTY(NCORD), DINTZ(NCORD), DINTNS(NCORD),
     &          AATX(NCORD), AATY(NCORD), AATZ(NCORD), AATXYZ(NCORD),
     &          ROTSTR(NCORD), NSPEC(NCORD), VINRTI(3,3), IFRQCM(NCORD),
     &          IMAGIN(NCORD), IREPS(NCORD),
     &          DGRAD(3,NCORD), CAAT(3,NCORD), WORK(LWORK),
     &          ROAAFQ(*), ROAGNQ(*), ROAGLQ(*), ROAAQ(*)
C
#include <molde.h>
#include <abainf.h>
#include <cbilnr.h>
#include <nuclei.h>
#include <symmet.h>
#include <pgroup.h>
#include <taymol.h>
#include <cbiwlk.h>
#include <gnrinf.h>
C
      LOGICAL  INTENS, NEWSPC, ACTIVE, CNTAIN, LINEAR
      DIMENSION DMYLAB(3), DMYMOL(3), AATMOL(3), NUMODX(8)
      CHARACTER*1 CHRIMG(0:1)
      CHARACTER*5 IRPTXT
      DATA CHRIMG /' ','i'/
C
C     *********************************************
C     ***** Diagonalize mass weighted Hessian *****
C     *********************************************
C
      KDKIN = 1
      KHESMW = KDKIN  + NCORD
      KWRK1  = KHESMW + NCORD*(NCORD + 1)/2
      LWRK1  = LWRK   - KWRK1 + 1
      CALL VIBNOR(HESS,AMASS,WRK(KDKIN),WRK(KHESMW),
     &            EVAL,EVEC,EVECS,WRK(KWRK1),LWRK1,
     &            NCORD,NUMMOD,NNEG,NOTALL,IPRINT)
      IF (IPRINT .GE. 8) THEN
         WRITE(LUPRI,'(/A)') ' Eigenvalues (Hartrees):'
         WRITE(LUPRI, * ) (II, EVAL(II),II = 1,NCORD)
      END IF
      IF (IPRINT .GE. 10) THEN
         WRITE(LUPRI,'(/A)')' The eigenvectors (normal coordinates):'
         CALL OUTPUT(EVEC,1,NCORD,1,NCORD,NCORD,NCORD,1,LUPRI)
      END IF
C
C     ******************************************
C     ***** Get frequencies and intensities ****
C     ******************************************
C
      IFREQ = 0
      XDMYXP = SQRT(XFAMU*XTKMML)
      DO 110 IMODE = 1, NUMMOD
         IFREQ = IFREQ + 1
         FREQ  = SQRT(ABS(EVAL(IFREQ)))
         FREQAU(IFREQ) = FREQ
         IFRQCM(IFREQ) = NINT(XTKAYS*FREQ)
         IMAGIN(IFREQ) = 0
         IF (EVAL(IFREQ) .LT. D0) IMAGIN(IFREQ) = 1
         IF (INTENS) THEN
            CALL DGEMM('N','N',3,1,NCORD,1.D0,
     &                 DGRAD,3,
     &                 EVEC(1,IFREQ),NCORD,0.D0,
     &                 DMYLAB,3)
c           CALL MPAB(DMYLAB,1,3,1,3, VINRTI,3,3,3,3, DMYMOL,1,3)
            call dcopy(3,dmylab,1,dmymol,1)
            DINTX(IFREQ)  = DMYMOL(1)
            DINTY(IFREQ)  = DMYMOL(2)
            DINTZ(IFREQ)  = DMYMOL(3)
            DINTNS(IFREQ) = V3DOT(DMYMOL,DMYMOL)
         END IF
         IF (DOVCD) THEN
            CALL DGEMM('N','N',3,1,NCORD,1.D0,
     &                 CAAT,3,
     &                 EVEC(1,IFREQ),NCORD,0.D0,
     &                 DMYLAB,3)
c               CALL MPAB(DMYLAB,1,3,1,3, VINRTI,3,3,3,3, AATMOL,1,3)
            call dcopy(3,dmylab,1,aatmol,1)
            AATX(IFREQ)   = AATMOL(1)
            AATY(IFREQ)   = AATMOL(2)
            AATZ(IFREQ)   = AATMOL(3)
            AATXYZ(IFREQ) = V3DOT(AATMOL,AATMOL)
            ROTSTR(IFREQ) = - V3DOT(AATMOL,DMYMOL)
         END IF
 110  CONTINUE
C
C     *********************
C     ***** Get irreps ****
C     *********************
C
      CALL IZERO(IREPS,NCORD)
      IF (7*NCORD .GT. LWORK) CALL STOPIT('VIBANA',' ',6*NCORD,LWORK)
      DO 150 IMODE = 1, NUMMOD
         FREQ = XTKAYS*FREQAU(IMODE)
         CALL DETIRP(EVEC,IMODE,FREQ,NCORD,WORK(1),WORK(NCORD+1),
     &        WORK(2*NCORD+1),WORK(3*NCORD+1),WORK(4*NCORD+1),
     &        WORK(5*NCORD+1),WORK(6*NCORD+1),LWORK-7*NCORD + 1,IRP)
         IREPS(IMODE) = IRP
 150  CONTINUE
C
C     *********************************************
C     ***** Print Frequencies and Intensities *****
C     *********************************************
C
      WRITE(LUPRI,'()')
      IF (INTENS) THEN
         CALL HEADER('Vibrational Frequencies and IR Intensities',1)
         WRITE (LUPRI,'(//1X,A,/,2X,60A,/,1X,A,/,2X,60A)')
     &  ' mode   irrep        frequency             IR intensity',
     &      ('=',I=1,60), '                ' //
     &  'cm-1       hartrees     km/mol   (D/A)**2/amu',
     &      ('-',I=1,60)
      ELSE
         CALL HEADER('Vibrational Frequencies',1)
            WRITE (LUPRI,'(2X,A,/,2X,A)')
     &    ' mode   irrep     cm-1     hartrees ',
     &    '------------------------------------'
      END IF
      WRITE (LUPRI,'()')
      IFREQ = 0
      DO 250 IRP = -1, MAXREP
         CNTAIN = .FALSE.
         DO 300 IMODE = 1, NUMMOD
            IF (IREPS(IMODE) .EQ. IRP) THEN
               IFREQ = IFREQ + 1
               CNTAIN = .TRUE.
               IF (IREPS(IMODE) .GE. 0) THEN
                  IRPTXT = ' ' // REP(IREPS(IMODE)) // ' '
               ELSE
                  IRPTXT = 'mixed'
               END IF
               IF (INTENS .AND. (IMAGIN(IMODE) .EQ. 0)) THEN
                  WRITE (LUPRI,2000) IMODE,IRPTXT,
     &                 XTKAYS*FREQAU(IMODE),CHRIMG(IMAGIN(IMODE)),
     &                 FREQAU(IMODE),CHRIMG(IMAGIN(IMODE)),
     &                 XFAMU*XTKMML*DINTNS(IMODE),
     &                 XFAMU*((DEBYE/XTANG)**2)*DINTNS(IMODE)
               ELSE IF (INTENS .AND. (IMAGIN(IMODE) .NE. 0)) THEN
                  WRITE (LUPRI,2050) IMODE,IRPTXT,
     &                 XTKAYS*FREQAU(IMODE),CHRIMG(IMAGIN(IMODE)),
     &                 FREQAU(IMODE),CHRIMG(IMAGIN(IMODE))
               ELSE
                  WRITE (LUPRI,2100) IMODE,IRPTXT,
     &                 XTKAYS*FREQAU(IMODE),CHRIMG(IMAGIN(IMODE)),
     &                 FREQAU(IMODE),CHRIMG(IMAGIN(IMODE))
               END IF
            END IF
 300     CONTINUE
         IF (CNTAIN) THEN
            WRITE (LUPRI,'()')
         END IF
 250  CONTINUE
C
C     ************************************
C     ***** Print normal coordinates *****
C     ************************************
C
      WRITE (LUPRI,'(//)')
      CALL HEADER('Normal Coordinates (bohrs*amu**(1/2)):',1)
      SXFAMU = SQRT(XFAMU)
      ISTR = 1
      NBATCH = (NUMMOD + 4)/5
      DO 400 IBATCH = 1, NBATCH
         IEND = MIN(ISTR + 4,NUMMOD)
         NUMB = IEND - ISTR + 1
         WRITE (LUPRI,'(/A12,5(I5,A2,I4,A1))') '            ',
     &      (II,'  ',IFRQCM(II),CHRIMG(IMAGIN(II)), II = ISTR,IEND)
         LENH = 10 + NUMB*12
         WRITE (LUPRI,'(2X,70A1)') ('-', II = 1,LENH)
         WRITE (LUPRI,'()')
         DO 410 ICOOR = 1, NCORD
            WRITE (LUPRI,1000) NAMDPX(ICOOR),
     &         (SXFAMU*EVEC(ICOOR,II),II=ISTR,IEND)
            IF (MOD(ICOOR,3) .EQ. 0) WRITE (LUPRI,'()')
 410     CONTINUE
         ISTR = ISTR + 5
 400  CONTINUE
      IF (MOLDEN) CALL MOFREQ(EVEC,NUMMOD,NCORD,FREQAU)
C
C     We punch out the normal coordinates on the file DALTON.NCA for
C     use with the Gamess-US visualization software
C
      LUIP = -1
      CALL GPOPEN(LUIP,'DALTON.NCA','UNKNOWN',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      WRITE(LUIP,8000)
      WRITE(LUIP,8010)
      WRITE(LUIP,8020) (AMASS(IAT),IAT=1,NATOMS)
      DO 910 IMODE = 1, NUMMOD
         WRITE(LUIP,8030) IMODE,FREQAU(IMODE)*XTKAYS
         WRITE(LUIP,8040) (SXFAMU*EVEC(I,IMODE),I=1,NCORD)
 910  CONTINUE
      WRITE(LUIP,8050)
      CALL GPCLOSE(LUIP,'KEEP')

C     We punch out harmonic freqs and normal coordinates on file
C     DALTON.NOR

      LUNOR = -1
      CALL GPOPEN(LUNOR,'DALTON.NOR','UNKNOWN',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      WRITE(LUNOR,'(A)') 'Harmonic Freqs. in cm^-1'
      DO 911 IMODE = 1, NUMMOD
        WRITE(LUNOR,'(E23.16)') FREQAU(IMODE)*XTKAYS
 911  CONTINUE
      WRITE(LUNOR,*)
      WRITE(LUNOR,'(A)') 'Normal Coordinates'
      DO 912 IMODE = 1, NUMMOD
        WRITE(LUNOR,8041) (SXFAMU*EVEC(I,IMODE),I=1,NCORD)
        WRITE(LUNOR,*)
 912  CONTINUE
      WRITE(LUNOR,'(A)') 'Norm of Vectors'
      DO 913 IMODE = 1, NUMMOD
        ENORM2 = 0.0D0
        DO 914 I = 1, NCORD
          ENORM2 = ENORM2 + EVEC(I,IMODE)**2
 914    CONTINUE
      WRITE(LUNOR,'(E23.16)') SXFAMU*sqrt(ENORM2)
 913  CONTINUE
      CALL GPCLOSE(LUNOR,'KEEP')

C
C     We punch out structure, frequencies, and normal coordinates 
C     for MidasCpp interface on file midasifc.coord 
C     

      LUMID  = -1
      CALL GPOPEN(LUMID,'midasifc.coord','UNKNOWN',' ','FORMATTED',
     &            IDUMMY,.FALSE.)

      WRITE(LUMID,'(I4,A)') NATOMS," au"
      CALL PRIGEOLU(LUMID,CORD)
C     DO ICENT=1,NUCDEP 
C        WRITE(LUMID,'(A,3(2X,E22.16))')
C    &         NAMEX(3*ICENT)(1:4),CORD(1,ICENT),
C    &         CORD(2,ICENT),CORD(3,ICENT)
C     ENDDO

      WRITE(LUMID,'(A,I5)') "FREQ",MAXREP+1
      WRITE(LUMID,'(A,A,E23.16)')
     &    '  cm-1                    irrep   ', ' F_au_to_cm-1:',XTKAYS
      DO 915 IMODE = 1, NUMMOD
        IF (IREPS(IMODE) .GE. 0) THEN
           IRPTXT = ' ' // REP(IREPS(IMODE)) // ' '
        ELSE
           IRPTXT = 'mixed'
        END IF
        WRITE(LUMID,'(E24.16,1X,A)') FREQAU(IMODE)*XTKAYS,
     &                                  IRPTXT
 915  CONTINUE
      WRITE(LUMID,'(A)') "COORD"
      WRITE(LUMID,'(A)') "au"
      DO 916 IMODE = 1, NUMMOD
        WRITE(LUMID,'(3(E24.16,1X))')
     &            (SXFAMU*EVEC(I,IMODE),I=1,NCORD)
        WRITE(LUMID,*)
 916  CONTINUE
      CALL GPCLOSE(LUMID,'KEEP')

C
C     ****************************************************************
C     ***** Normal coordinates in redundant internal coordinates *****
C     ****************************************************************
C
C     This analysis is only performed if the OPTIMIZE module is used
C     and the redundant internal coordinates have been set up.
C
      IF (OPTNEW) THEN
C
C     Cannot include optinf.h due to conflicts, need to get number
C     of Cartesian and internal coordinates by calling NUMCRD
C
         CALL NUMCRD(ICRTCR,IINTCR)
         IF (IINTCR .GT. 0) THEN
            MXRCRD = MAX(MAX(IINTCR, ICRTCR), 8)
            MX2CRD = MAX(MXCOOR,MXRCRD)
            KATMAR = 1
            KICRD  = KATMAR + 8*MXCENT
            KWILBM = KICRD  + MXRCRD
            KBMTRA = KWILBM + MXRCRD*MXCOOR
            KTMPMT = KBMTRA + MXRCRD*MXRCRD
            KTMPM2 = KTMPMT + MX2CRD*MX2CRD
            KLAST  = KTMPM2 + MX2CRD*MX2CRD
            KWRK1  = KLAST
            IF (KWRK1 .GE. LWRK) CALL STOPIT('VIBANA','REDVIB',KWRK1,
     &           LWRK)
            CALL REDVIB(NCORD,NUMMOD,MXRCRD,MX2CRD,EVEC,IFRQCM,IMAGIN,
     &           WORK(KATMAR),WORK(KICRD),WORK(KWILBM),WORK(KBMTRA),
     &           WORK(KTMPMT),WORK(KTMPM2))
         END IF
      END IF
C
C
C     *****************************************************************
C     ***** Dipole derivatives with respect to normal coordinates *****
C     *****************************************************************
C
      IF (INTENS) THEN
         WRITE (LUPRI,'()')
         CALL HEADER('Dipole Gradient in '//
     &               'Normal Coordinate Basis (D/(A*amu**(1/2)))',1)
         WRITE (LUPRI,'(2(/,1X,A))')
     &      ' mode           dMA/dQi        dMB/dQi        dMC/dQi',
     &      ' ----------------------------------------------------'
         FAC = SQRT(XFAMU)*(DEBYE/XTANG)
         IFREQ = 0
         DO 700 IMODE = 1, NUMMOD
            IFREQ = IFREQ + 1
            WRITE (LUPRI,'(2X,I2,5X,3F15.6)') IFREQ,
     &           FAC*DINTX(IFREQ), FAC*DINTY(IFREQ), FAC*DINTZ(IFREQ)
 700     CONTINUE
         WRITE (LUPRI,'()')
      END IF
C
C     *******************************************
C     ***** Dipole and rotational strengths *****
C     *******************************************
C
      IF (DOVCD) THEN
         WRITE (LUPRI,'()')
         CALL HEADER('Dipole and rotational strengths',1)
         WRITE (LUPRI,'(2X,A,/2X,A)')
     &      'Units: 10**(-40) (esu**2)*(cm**2) (dipole strength)',
     &      '       10**(-44) (esu**2)*(cm**2) (rotational strength)'
         WRITE (LUPRI,'(2(/,1X,A))')
     &      ' mode       frequency      dip. str.      rot. str.',
     &      ' --------------------------------------------------'
         IFREQ = 0
         DO 740 IMODE = 1, NUMMOD
            IFREQ = IFREQ + 1
            WRITE (LUPRI,'(2X,I2,2X,F15.2,2F15.3)')
     &           IFREQ, XTKAYS*FREQAU(IFREQ),
     &           ESUDIP*DINTNS(IFREQ)/(D2*FREQAU(IFREQ)),
     &           ESUROT*ROTSTR(IFREQ)
 740     CONTINUE
         WRITE (LUPRI,'()')
      END IF
C
      IF (VROA .OR. RAMAN) THEN
         CALL VIBROA(EVEC,FREQAU,ROAAFQ,ROAGNQ,ROAGLQ,ROAAQ,WORK,LWORK,
     &               NCORD,NUMMOD)
      END IF
      IF (V3CAL .AND. .NOT. NMODIF) THEN
         KFMATF = KWRK1
         KFMATT = KFMATF + NCORD*NCORD*NCORD
         KCORR  = KFMATT + NCORD*NCORD*NCORD
         KPOS   = KCORR  + NCORD
         KLAST  = KPOS   + NCORD
         LLEFT  = LWRK   - KLAST
         CALL VIBV3(EVEC,FREQAU,WORK(KFMATF),
     &              WORK(KFMATT),WORK(KCORR),WORK(KPOS),AMASS,
     &              WORK(KLAST),LLEFT,NCORD,NUMMOD)
      END IF
      IF (LINCPL) THEN
C
C     Transform gradient into normal coordinates (gradient obtained
C     for a different electronic state than the Hessian). For use with
C     the linear coupling model of Franck-Condon factors
C
         CALL DGEMM('N','N',1,NUMMOD,NCORD,1.D0,
     &              GRAD,1,
     &              EVEC,NCORD,0.D0,
     &              WORK(KFMATF),1)
C
         CALL TITLER('Linear coupling model gradient analysis','*',118)
         WRITE (LUPRI,'(/15X,A,/15X,A)')
     &        '  Normal mode    Frequency     Gradient',
     &        '  --------------------------------------'
         DO IMOD = 1, NUMMOD
            WRITE (LUPRI,'(19X,I4,10X,F7.2,6X,F8.6)') IMOD,
     &           FREQAU(IMOD)*XTKAYS,WORK(KFMATF + IMOD - 1)*SQRT(XFAMU)
         END DO
      END IF
C
C     *******************************
C     ***** Zero-point energies *****
C     *******************************
C
      IF (.NOT.NOTALL) THEN
         EVIB = D0
         DO 800 IFREQ = 1,NUMMOD
            IF (EVAL(IFREQ).GT.1.D-8) EVIB = EVIB + SQRT(EVAL(IFREQ))
  800    CONTINUE
         EVIB = DP5 * EVIB
         CALL HEADER ('Total Molecular Energy',-1)
         WRITE (LUPRI,'(12X,A/)')
     &      ' electronic     vibrational           total    energy '
         WRITE (LUPRI,'(10X,3(F13.6,3X),A)')
     &      ERGMOL, EVIB, ERGMOL + EVIB, ' Hartrees'
         WRITE (LUPRI,'(10X,3(F13.2,3X),A)')
     &      XTKAYS*ERGMOL, XTKAYS*EVIB, XTKAYS*(ERGMOL + EVIB), ' cm-1'
         WRITE (LUPRI,'(10X,3(F13.2,3X),A)')
     &      XKJMOL*ERGMOL, XKJMOL*EVIB,XKJMOL*(ERGMOL + EVIB), ' kJ/mol'
         WRITE (LUPRI,'(10X,3(F13.2,3X),A)')
     &      XKCMOL*ERGMOL,XKCMOL*EVIB,XKCMOL*(ERGMOL + EVIB),' kcal/mol'
      END IF
C
C     End of VIBANA
C
      RETURN
C
C     Print formats
C
 1000 FORMAT (4X,A,(T13,5F12.6))
 2000 FORMAT (2X,I2,5X,A5,F10.2,A,2X,F9.6,A,2X,F8.3,F9.4,
     &     10X,I1,7(A,I1))
 2050 FORMAT (2X,I2,5X,A5,F10.2,A,2X,F9.6,A,29X,I1,7(A,I1))
 2100 FORMAT (3X,I2,5X,A5,F11.2,A,F10.7,A,5X,I2,7(A,I2))
C
C     Punch formats
C
 8000 FORMAT('----- START OF NORMAL MODES FOR -MOLPLT- PROGRAM -----')
 8010 FORMAT('ATOMIC MASSES')
 8020 FORMAT(5F12.5)
 8030 FORMAT('MODE',I5,'   FREQUENCY=',F10.5,' (CM**-1)')
 8040 FORMAT(1P,3E17.9)
 8041 FORMAT(1P,3E23.16)
 8050 FORMAT('----- END OF NORMAL MODES FOR -MOLPLT- PROGRAM -----')
C
      END
C  /* Deck vibnor */
      SUBROUTINE VIBNOR(HESS,AMASS,DKIN,HESSMW,
     &                  EVAL,EVEC,EVECS,WRK,LWRK,
     &                  NCORD,NINTRN,NNEG,NOTALL,IPRINT)
#include <implicit.h>
C
C     Find normal coordinates in DNORCR and normal frequences in FREQ
C     for projected Hessian eigenvalue problem
C     (HESS - EVAL * DKIN) * EVEC = 0
C
C     Modification KR 29.10.96 - Sort only according to frequency, not symmetry
C     Modificaton TUH 26.12.86 - Sort according to symmetry
C     Modified TUH 21.09.89 - Explicit Abelian symmetry
C
C     Input:   HESS,AMASS
C     Output:  EVAL,EVEC,EVECS,NNEG
C
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <codata.h>
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0 ,
     &            THRSH1 = 1.0D-14, THRSH2 = 1.0D-3, THRSH3 = 1.0D-7)
      LOGICAL NOTALL
C
      DIMENSION HESS(NCORD,NCORD),
     &          AMASS(*), DKIN(NCORD), HESSMW(NCORD*(NCORD+1)/2),
     &          EVAL(NCORD), EVEC(NCORD,NCORD), EVECS(NCORD,NCORD),
     &          WRK(LWRK)
C
#include <abainf.h>
#include <symmet.h>
C
C     *******************************************
C     ***** Construct mass-weighted Hessian *****
C     *******************************************
C
      DO 100 I = 1, NCORD
         DKIN(I) = D1/SQRT(XFAMU*AMASS((I+2)/3))
 100  CONTINUE
      IJ = 0
      DO 200 I = 1, NCORD
         DO 210 J = 1, I
            IJ = IJ + 1
            HESSMW(IJ) = DKIN(J)*HESS(J,I)*DKIN(I)
 210     CONTINUE
 200  CONTINUE
C
C     *********************************************
C     ***** Diagonalize mass-weighted Hessian *****
C     *********************************************
C
      KWRK  = 1
      KIWRK = KWRK + NCORD
      CALL DUNIT(EVEC,NCORD)
      CALL JACO(HESSMW,EVEC,NCORD,NCORD,NCORD,WRK(KWRK),WRK(KIWRK))
C
C     **********************************************
C     ***** Eigenvalues and normal coordinates *****
C     **********************************************
C
      II = 0
      DO 300 I = 1,NCORD
         II = II + I
         EVAL(I) = HESSMW(II)
         DO 310 J = 1,NCORD
            EVEC(J,I) = DKIN(J)*EVEC(J,I)
 310     CONTINUE
 300  CONTINUE
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Eigenvalues of mass-weighted Hessian',-1)
         CALL OUTPUT(EVAL,1,1,1,NCORD,1,NCORD,1,LUPRI)
         CALL HEADER('Normal coordinates in Cartesian basis',-1)
         CALL OUTPUT(EVEC,1,NCORD,1,NCORD,NCORD,NCORD,1,LUPRI)
      END IF
C
C     *******************************************
C     ***** Eliminate non-vibrational modes *****
C     *******************************************
C
      NZER = 0
      NNEG = 0
C/djw/ VALTHR set to 1.0D-10
C     VALTHR = THRSH1
C     IF (NUMHES) VALTHR = 1.0D-10
      VALTHR = 1.0D-10
      DO 400 I = 1, NCORD
         IF (ABS(EVAL(I)) .LE.  VALTHR) NZER = NZER + 1
         IF (    EVAL(I)  .LT. -VALTHR) NNEG = NNEG + 1
 400  CONTINUE
      CALL ORDER2(EVEC,EVAL,NCORD,NCORD)
      IZERO  = NCORD - NNEG - NZER + 1
      DO 500 INEG = NCORD - NNEG + 1, NCORD
         CALL DSWAP(NCORD,EVEC(1,INEG),1,EVEC(1,IZERO),1)
         CALL DSWAP(1,EVAL(INEG),1,EVAL(IZERO),1)
         IZERO = IZERO + 1
 500  CONTINUE
      NINTRN = NCORD - NZER
      IF (.NOT.NOTALL .AND. (NZER .NE. 5 .AND. NZER .NE. 6)) THEN
         WRITE (LUPRI,'(/A,I2,/,A)') ' The number of translational'//
     &      ' and rotational coordinates in VIBNOR',NZER,' is'//
     &      ' not equal to 5 or 6.'
      END IF
      IF (IPRINT .GT. 5) THEN
         CALL HEADER
     &      ('Non-zero eigenvalues of mass-weighted Hessian',-1)
         CALL OUTPUT(EVAL,1,1,1,NINTRN,1,NCORD,1,LUPRI)
         CALL HEADER
     &      ('Non-zero normal coordinates in Cartesian basis',-1)
         CALL OUTPUT(EVEC,1,NCORD,1,NINTRN,NCORD,NCORD,1,LUPRI)
      END IF
C
C     **************************************
C     ***** Transform to symmetry basis ****
C     **************************************
C
      KCSTRA = 1
      KSCTRA = KCSTRA + NCORD*NCORD
      KLAST  = KSCTRA + NCORD*NCORD
      IF (KLAST .GT. LWRK) CALL STOPIT('VIBNOR','TRACOR',KLAST,LWRK)
      CALL TRACOR(WRK(KCSTRA),WRK(KSCTRA),1,NCORD,0)
C
      CALL DGEMM('T','N',NCORD,NCORD,NCORD,1.D0,
     &           WRK(KSCTRA),NCORD,
     &           EVEC,NCORD,0.D0,
     &           EVECS,NCORD)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Normal coordinates in symmetry basis',-1)
         CALL OUTPUT(EVECS,1,NCORD,1,NCORD,NCORD,NCORD,1,LUPRI)
      END IF
C
C     **********************************************************
C     ***** Sort normal modes according to their frequency *****
C     **********************************************************
C
      DO 600 ICOOR = 1, NINTRN
         DO 610 JCOOR = 1, NINTRN
            IF (EVAL(ICOOR) .GT. EVAL(JCOOR)) THEN
               CALL DSWAP(1,EVAL(ICOOR),1,EVAL(JCOOR),1)
               CALL DSWAP(NCORD,EVEC(1,ICOOR),1,EVEC(1,JCOOR),1)
               IF (MAXREP .GT. 1) THEN
                  CALL DSWAP(NCORD,EVECS(1,ICOOR),1,EVECS(1,JCOOR),1)
               END IF
            END IF
 610     CONTINUE
 600  CONTINUE
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Sorted eigenvalues',-1)
         CALL OUTPUT(EVAL,1,1,1,NINTRN,1,NCORD,1,LUPRI)
         CALL HEADER('Sorted normal coordinates in Cartesian basis',-1)
         CALL OUTPUT(EVEC,1,NCORD,1,NINTRN,NCORD,NCORD,1,LUPRI)
         CALL HEADER('Sorted normal coordinates in symmetry basis',-1)
         CALL OUTPUT(EVECS,1,NCORD,1,NINTRN,NCORD,NCORD,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck prtfun */
      SUBROUTINE PRTFUN(TOTMAS,AIMOM,BIMOM,CIMOM,EVAL,
     &                  LINEAR,NCORD,NUMMOD)
#include <implicit.h>
C
C     Calculate translational, rotational, and vibrational
C     partition functions
C
#include <priunit.h>
#include <codata.h>
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, D12 = 12.0D0)
      PARAMETER (D1P5 = 1.5D0, D2P5 = 2.5D0)
      PARAMETER (RPFFAC = 0.041 229 342 D0,
     &           TPFFAC = 0.025 607 3D0)
      PARAMETER ( MAXTMP = 15 )
C
      LOGICAL  LINEAR
      DIMENSION EVAL(NCORD)
C
      DIMENSION TEMP(MAXTMP)
      DATA TEMP/ 50.0 D0, 100.0D0, 200.0D0, 273.15D0, 298.15D0,
     &          300.0 D0, 400.0D0, 500.0D0, 600.0 D0, 700.0 D0,
     &          800.0 D0, 900.0D0,1000.0D0,1500.0 D0,2000.0 D0/
C
C     ***** Print *****
C
      CALL HEADER('Molecular Partition Functions',-1)
      WRITE (LUPRI,'(4(A,/),//A/)')
     & ' Qtran is evaluated per mol at 1 atm. pressure.',
     & ' Qrot does not include symmetry numbers.',
     & ' Qvib does not include zero point energies',
     & ' (i.e., energy scale has vibrational ground state as zero).',
     & ' Temp.(K)  Translational     Rotational   '//
     & ' Vibrational      Total'
C
C     For Non-linear molecules, we need a constant, so we do this
C     before the loop over temperatures
C
      IF (.NOT.LINEAR) THEN
	            QRCOR  = ( D2/AIMOM - AIMOM/(BIMOM*CIMOM)
     &               + D2/BIMOM - BIMOM/(AIMOM*CIMOM)
     &               + D2/CIMOM - CIMOM/(BIMOM*AIMOM) )
      END IF
C
C     Loop over temperatures
C
      DO 100 ITEMP = 1,MAXTMP
         TEMPI  = TEMP(ITEMP)
C
C        ***** Translational partition function *****
C
         QTRANS = TPFFAC*(TOTMAS**D1P5)*(TEMPI**D2P5)
C
C        ***** Rotational partition function *****
C
C        Reference for correction:
C        Stripp and Kirkwood, J.Chem.Phys. 19(1951)1131
C
         IF (LINEAR) THEN
            QBCLAS = SQRT(RPFFAC*BIMOM*TEMPI)
            QRCLAS = QBCLAS * QBCLAS
            QROT   = RPFLIN(D1/QRCLAS)
          ELSE
            QACLAS = SQRT(RPFFAC*AIMOM*TEMPI)
            QBCLAS = SQRT(RPFFAC*BIMOM*TEMPI)
            QCCLAS = SQRT(RPFFAC*CIMOM*TEMPI)
            QRCLAS = SQRT(PI)*QACLAS*QBCLAS*QCCLAS
            QROT   = QRCLAS*(D1 + QRCOR / (D12*RPFFAC*TEMPI) )
         END IF
C
C        ***** Vibrational partition function *****
C
         QVIB = D1
         DO 200 IFREQ = 1,NUMMOD
            IF (EVAL(IFREQ) .GT. 1.D-8) THEN
               EVIBI = SQRT(EVAL(IFREQ))
               UI    = AUTK*(EVIBI/TEMPI)
               QVIB  = QVIB/(D1 - EXP(-UI))
            END IF
  200    CONTINUE
C
C        ***** Total partition function *****
C
         QTOT   = QTRANS * QROT * QVIB
         WRITE (LUPRI,'(1X,F8.2,1P,D15.4,0PF15.2,F15.4,1P,D15.6)')
     &          TEMPI, QTRANS, QROT, QVIB, QTOT
  100 CONTINUE
C
C     End of PRTFUN
C
      RETURN
      END
C  /* Deck vibpol */
      SUBROUTINE VIBPOL(DGRAD,EVAL,EVEC,VINRTI,POLAR,NCORD,NUMMOD)
#include <implicit.h>
C
C     Vibrational polarization, implemented 19-Jul-1986 HJAaJ
C     Ref: Rinaldi et al., Chem.Phys.Lett. 128 (1986) 177-181
C
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
      PARAMETER (D1 = 1.0D0, PVFAC  = 0.124 748 D0)
C
      LOGICAL POLAR, ALLXYZ
      DIMENSION DGRAD(3,NCORD), VINRTI(3,3), EVAL(NCORD), EVEC(NCORD,*)
#include <moldip.h>
#include <symmet.h>
#include <dorps.h>
C
      DIMENSION ALFVIB(3,3), ALFELE(3,3), AMAT(3,3), AVEC(3)
C
      IF (.NOT.(DOREPS(ISYMAX(1,1)) .OR. DOREPS(ISYMAX(2,1))
     &                              .OR. DOREPS(ISYMAX(3,1)))) RETURN
C
      ALLXYZ = DOREPS(ISYMAX(1,1)) .AND. DOREPS(ISYMAX(2,1))
     &                             .AND. DOREPS(ISYMAX(3,1))
C
      CALL DZERO(ALFVIB,9)
C
      NNEG = 0
      DO 100 IFREQ = 1, NUMMOD
         IF (EVAL(IFREQ) .LT. -1.D-8) NNEG = NNEG + 1
 100  CONTINUE
      IF (NNEG .EQ. 0) THEN
         DO 200 IFREQ = 1, NUMMOD
            CALL DGEMM('N','N',3,1,NCORD,1.D0,
     &                 DGRAD,3,
     &                 EVEC(1,IFREQ),NCORD,0.D0,
     &                 AVEC,3)
            CALL DGEMM('N','T',3,3,1,1.D0,
     &                 AVEC,3,
     &                 AVEC,3,0.D0,
     &                 ALFELE,3)
            FAC = D1/EVAL(IFREQ)
            CALL DAXPY(9,FAC,ALFELE,1,ALFVIB,1)
 200     CONTINUE
         CALL AROUND('Polarizabilities')
         WRITE (LUPRI,'(/A/A)')
     &      ' Vibrational polarizabilities calculated acc. to'//
     &      ' Rinaldi et al.,',' Chem.Phys.Lett. 128 (1986) 177-181.'
C
         CALL TRAPOL(ALFVIB,'TOSYM')
         CALL HEADER('Vibrational Polarizabilities (au)',-1)
         CALL POLPRI(ALFVIB,'   ',1)
         CALL HEADER('Vibrational Polarizabilities (angstroms**3)',-1)
         CALL POLPRI(ALFVIB,'EXP',1)
C
         CALL TRAPOL(ALFVIB,'FROMSYM')
         CALL DGEMM('N','N',3,3,3,1.D0,
     &              ALFVIB,3,
     &              VINRTI,3,0.D0,
     &              AMAT,3)
         CALL DGEMM('T','N',3,3,3,1.D0,
     &              VINRTI,3,
     &              AMAT,3,0.D0,
     &              ALFVIB,3)
         CALL HEADER('Vibrational Polarizabilities (angstroms**3)',-1)
         CALL POLPRI(ALFVIB,'PRIEXP',1)
         IF (ALLXYZ) THEN
            PVIB = PVFAC * (ALFVIB(1,1) + ALFVIB(2,2) + ALFVIB(3,3))
            WRITE (LUPRI,'(//A,F10.4,A)')
     &         ' Vibrational mean polarization :',PVIB,' cm**3 mol**-1'
         END IF
         IF (POLAR) THEN
            CALL DCOPY(9,POLFLT,1,ALFELE,1)
            CALL TRAPOL(ALFELE,'FROMSYM')
            CALL DGEMM('N','N',3,3,3,1.D0,
     &                 ALFELE,3,
     &                 VINRTI,3,0.D0,
     &                 AMAT,3)
            CALL DGEMM('T','N',3,3,3,1.D0,
     &                 VINRTI,3,
     &                 AMAT,3,0.D0,
     &                 ALFELE,3)
            CALL HEADER('Electronic Polarizabilities (angstroms**3)',-1)
            CALL POLPRI(ALFELE,'PRIEXP',1)
            CALL DAXPY(9,D1,ALFVIB,1,ALFELE,1)
            CALL HEADER('Total Polarizabilities (angstroms**3)',-1)
            CALL POLPRI(ALFELE,'PRIEXP',1)
         END IF
      END IF
C
C     End of VIBPOL
C
      RETURN
      END
C  /* Deck vibhes */
      SUBROUTINE VIBHES(IPRINT,NCORD,COORD,KTRRO,GRDPRJ,HESPRJ,TRROVE,
     &                  NATTYP,WRK)
C
C     Purpose:
C
C      Project out translational and rotational coordinates from the
C      molecular Hessian
C
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (THRLDP = 1.0D-4, D0 = 0.D0, DM1 = -1.0D0, D1 = 1.0D0)
      DIMENSION COORD(*), GRDPRJ(*), HESPRJ(NCORD,*), TRROVE(NCORD,*),
     &          NATTYP(*), WRK(*)
      LOGICAL   LINDEP(6)
C
#include <taymol.h>
#include <cbisol.h>
      AMASS(I) = DISOTP(NATTYP(I),1,'MASS')
C
C     Construct 6 vectors that span the translational and rotational
C     vector space
C
C     Translation:
C
C     (1,0,0)  ; (0,1,0)  ; (0,0,1)
C
C     Rotation (contain translational components):
C
C     (-Y,X,0) ; (0,-Z,Y) ; (Z,0,-X)
C
      NCAV = 0
      IF (SOLVNT) NCAV = 3
C
      NDIM = 6 + NCAV
      NTOT = NCORD*NDIM
      CALL DZERO(TRROVE,NTOT)
      NATOM = NCORD/3
      IJ = 1
      DO 1100 J = 1,NATOM
         COOX = COORD(IJ)
         COOY = COORD(IJ+1)
         COOZ = COORD(IJ+2)
         TRROVE(IJ,  1) = D1
         TRROVE(IJ+1,2) = D1
         TRROVE(IJ+2,3) = D1
         TRROVE(IJ,  4) = -COOY
         TRROVE(IJ+1,4) =  COOX
         TRROVE(IJ+1,5) = -COOZ
         TRROVE(IJ+2,5) =  COOY
         TRROVE(IJ,  6) =  COOZ
         TRROVE(IJ+2,6) = -COOX
         IJ  = IJ + 3
 1100 CONTINUE
C
C     Center of mass motion
C
      IF (SOLVNT) THEN
C
C        Total mass
C
         TMASS = D0
         DO 100 I = 1, NATOM - 1
            TMASS = TMASS + AMASS(I)
  100    CONTINUE
C
         IJ = 1
         DO 200 J = 1, NATOM - 1
            TRROVE(IJ,  7) = - AMASS(J)/TMASS
            TRROVE(IJ+1,8) = - AMASS(J)/TMASS
            TRROVE(IJ+2,9) = - AMASS(J)/TMASS
            IJ  = IJ + 3
  200    CONTINUE
         TRROVE(3*(NATOM-1) + 1,7) = D1
         TRROVE(3*(NATOM-1) + 2,8) = D1
         TRROVE(3*(NATOM-1) + 3,9) = D1
      END IF
      IF (IPRINT .GE. 15) THEN
         WRITE(LUPRI,*) '  TRROVE'
         CALL OUTPUT(TRROVE,1,NCORD,1,NDIM,NCORD,NDIM,1,LUPRI)
      END IF
C
C     Orthogonalize the trans rot vectors
C
      KTRRO = NDIM
      CALL ORTVEC(0,KTRRO,NCORD,THRLDP,TRROVE,LINDEP)
C
      IF (SOLVNT) NCAV = 3
      NTRRO = KTRRO - NCAV
      IF (IPRINT .GE. 5) THEN
         IF (NTRRO.EQ.6) WRITE(LUPRI,'(/I5,A,/A)')
     &       NTRRO,' Rotational and translational degrees of freedom',
     &             '    **** non-linear geometry ****'
         IF (NTRRO.EQ.5) WRITE(LUPRI,'(/I5,A,/A)')
     &       NTRRO,' Rotational and translational degrees of freedom',
     &             '    **** linear geometry ****'
      END IF
      IF (NTRRO.NE.5.AND.NTRRO.NE.6) THEN
         WRITE (LUPRI,1550) NTRRO
         CALL QUIT('(ABACUS.VIBHES) Incorrect no. of tra/rot coord.')
      END IF
 1550 FORMAT(/,' *** VIBHES ***  ERROR ',
     &       /I5,' translation and rotational coordinates found')
C
C     Project out rotational and translational degrees of freedom
C     from the molecular Hessian
C
C     Construct the projection operator
C
      CALL DGEMM('N','T',NCORD,NCORD,KTRRO,1.D0,
     &           TRROVE,NCORD,
     &           TRROVE,NCORD,0.D0,
     &           WRK,NCORD)
      IF (IPRINT .GE. 15) THEN
         CALL HEADER(' PROJECTION OPERATOR ',-1)
         CALL OUTPUT(WRK,1,NCORD,1,NCORD,NCORD,NCORD,1,LUPRI)
      END IF
      DO 2920 I = 1,NCORD
         IDIAG = (I-1)*NCORD + I
         WRK(IDIAG) = WRK(IDIAG) - D1
 2920 CONTINUE
      NTOT = NCORD*NCORD
      CALL DSCAL(NTOT,DM1,WRK,1)
      IF (IPRINT .GE. 15) THEN
         CALL HEADER('1 - PROJECTION OPERATOR ',-1)
         CALL OUTPUT(WRK,1,NCORD,1,NCORD,NCORD,NCORD,1,LUPRI)
         CALL HEADER('MOLECULAR HESSIAN HESPRJ',-1)
         CALL OUTPUT(HESPRJ,1,NCORD,1,NCORD,NCORD,NCORD,1,LUPRI)
      END IF
C
C     Calculate norm of projected gradient
C
      KSTRT = NCORD*NCORD + 1
      CALL DGEMM('N','N',NCORD,1,NCORD,1.D0,
     &           WRK,NCORD,
     &           GRDPRJ,NCORD,0.D0,
     &           WRK(KSTRT),NCORD)
      GRDNRM = DDOT(NCORD,WRK(KSTRT),1,WRK(KSTRT),1)
      GRDNRM = SQRT(GRDNRM)
      WRITE (LUPRI,'(1X,A,F12.6)')
     &   ' Norm of projected gradient ',GRDNRM
      IF (IPRINT .GE. 15) THEN
         CALL HEADER(' PROJECTED GRADIENT ',-1)
         CALL OUTPUT(WRK(KSTRT),1,1,1,NCORD,1,NCORD,1,LUPRI)
      END IF
C
C     Project trans-rot out of molecular Hessian
C
      CALL DGEMM('N','N',NCORD,NCORD,NCORD,1.D0,
     &           WRK,NCORD,
     &           HESPRJ,NCORD,0.D0,
     &           WRK(KSTRT),NCORD)
      CALL DGEMM('N','N',NCORD,NCORD,NCORD,1.D0,
     &           WRK(KSTRT),NCORD,
     &           WRK,NCORD,0.D0,
     &           HESPRJ,NCORD)
      IF (IPRINT .GE. 15) THEN
         CALL HEADER(' PROJECTED HESSIAN ',-1)
         CALL OUTPUT(HESPRJ,1,NCORD,1,NCORD,NCORD,NCORD,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck hdot */
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (VAR_REAL)
      REAL FUNCTION HDOT(NDIM,VEC1,VEC2,HMAT)
#else
      DOUBLE PRECISION FUNCTION HDOT(NDIM,VEC1,VEC2,HMAT)
#endif
#include <implicit.h>
      PARAMETER (D0 = 0.0D0)
      DIMENSION VEC1(NDIM), VEC2(NDIM), HMAT(NDIM,NDIM)
      HDOT = D0
      DO 100 I = 1, NDIM
         DO 200 J = 1, NDIM
            HDOT = HDOT + VEC1(J)*HMAT(J,I)*VEC2(I)
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck intern */
      SUBROUTINE INTERN(NCORD,NINTER,GRAD,HESS,HESCOR,DGRAD,CAAT,BMAT,
     &                  BADJIN,QMAT,QINV,BHES,FINT,DINT,AATINT,HINT,W2,
     &                  IWORK,CMAT,QQ,CTYPE,WRK,LWORK)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <codata.h>
#include <facang.h>
      PARAMETER (D0 = 0.D0, D1 = 1.0D0, DM1 = - 1.0D0)
C
#include <cbivib.h>
#include <nuclei.h>
C
      CHARACTER UNITC*3, UNITF*9, UNITD*5
      DIMENSION GRAD(NCORD), HESS(NCORD,NCORD), DGRAD(3,NCORD),
     &          CAAT(3,NCORD), BMAT(NINTER,NCORD), BADJIN(NINTER,NCORD),
     &          QMAT(NINTER,NINTER), QINV(NINTER,NINTER),
     &          HESCOR(NCORD,*), BHES(NINTER,NCORD),
     &          FINT(NINTER), DINT(NINTER,3), HINT(NINTER,NINTER),
     &          CMAT(NINTER,NCORD,NCORD), QQ(NINTER), AATINT(NINTER,3),
     &          CTYPE(NINTER), W2(NINTER), IWORK(NINTER), WRK(*)
C
C
      CALL AROUND('Analysis using curvilinear internal coordinates')
      WRITE (LUPRI,'(2(A/),//)')
     & ' (Subroutine MACHB for transformation to internal coordinates',
     & '  has been supplied by Peter Pulay.)'
C
C     ********************
C     ***** B matrix *****
C     ********************
C
      KWKXA = 1
      KWKYA = KWKXA + NCORD
      KWKZA = KWKYA + NCORD
      KWKBD = KWKZA + NCORD
      LGETB = KWKBD + NINTER*NCORD
      CALL GETB(IPRINT,NINTER,NCORD,BMAT,QQ,CTYPE,WRK(KWKXA),WRK(KWKYA),
     &          WRK(KWKZA),WRK(KWKBD))
      IF (IPRINT .GE. 2) THEN
         CALL HEADER('The B matrix (Bij = dqi/dXj)',-1)
         CALL OUTPUT(BMAT,1,NINTER,1,NCORD,NINTER,NCORD,1,LUPRI)
      END IF
C
C     ***********************************************
C     ***** BADJIN - pseudo-inverse of B matrix *****
C     ***********************************************
C
C     Construct QMAT = BMAT(I,K)*BMAT(J,K)
C
      CALL DGEMM('N','T',NINTER,NINTER,NCORD,1.D0,
     &           BMAT,NINTER,
     &           BMAT,NINTER,0.D0,
     &           QMAT,NINTER)
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('QMAT',-1)
         CALL OUTPUT(QMAT,1,NINTER,1,NINTER,NINTER,NINTER,1,LUPRI)
      END IF
C
C     Invert QMAT
C
      CALL DGEINV(NINTER,QMAT,QINV,IWORK,W2,INFO)
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('QINV',-1)
         CALL OUTPUT(QINV,1,NINTER,1,NINTER,NINTER,NINTER,1,LUPRI)
      END IF
      IF (INFO .NE. 0) THEN
         WRITE (LUPRI,'(//,A,I5,A,/)')
     &      ' ERROR (HESINT) INFO =',INFO,' from DGEINV '
         CALL QUIT('(ABACUS.VIBCTL.INTERN) ERROR in DGEINV')
      END IF
C
C     BADJIN
C
      CALL DGEMM('N','N',NINTER,NCORD,NINTER,1.D0,
     &           QINV,NINTER,
     &           BMAT,NINTER,0.D0,
     &           BADJIN,NINTER)
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('BADJIN',-1)
         CALL OUTPUT(BADJIN,1,NINTER,1,NCORD,NINTER,NCORD,1,LUPRI)
      END IF
C
C     *************************************
C     ***** Calculate internal forces *****
C     *************************************
C
      CALL DGEMM('N','N',NINTER,1,NCORD,1.D0,
     &           BADJIN,NINTER,
     &           GRAD,NCORD,0.D0,
     &           FINT,NINTER)
      CALL DSCAL(NINTER,DM1,FINT,1)
C
C     *************************************************
C     ***** Calculate internal dipole derivatives *****
C     *************************************************
C
      IF (DIPOL) THEN
         CALL DGEMM('N','T',NINTER,3,NCORD,1.D0,
     &              BADJIN,NINTER,
     &              DGRAD,3,0.D0,
     &              DINT,NINTER)
      END IF
C
C     ***********************************
C     ***** Calculate internal AATs *****
C     ***********************************
C
      IF (DOVCD) THEN
         CALL DGEMM('N','T',NINTER,3,NCORD,1.D0,
     &              BADJIN,NINTER,
     &              CAAT,3,0.D0,
     &              AATINT,NINTER)
      END IF
C
C     *************************************
C     ***** Calculate force constants *****
C     *************************************
C
      IF (MAXDIF .EQ. 2) THEN
C
C        Differentiated B matrix
C
         KWKXA = 1
         KWKYA = KWKXA + NCORD
         KWKZA = KWKYA + NCORD
         KWKBP = KWKZA + NCORD
         KWKBM = KWKBP + NINTER*NCORD
         LGETC = KWKBM + NINTER*NCORD
         CALL GETCMAT(IPRINT,NINTER,NCORD,CMAT,QQ,CTYPE,WRK(KWKXA),
     &             WRK(KWKYA),WRK(KWKZA),WRK(KWKBP),WRK(KWKBM))
C
C        Correction matrix
C
         CALL DGEMM('N','N',1,NCORD*NCORD,NINTER,1.D0,
     &              FINT,1,
     &              CMAT,NINTER,0.D0,
     &              HESCOR,1)
         IF (IPRINT .GE. 5) THEN
            CALL HEADER('Correction to Hessian',-1)
            CALL OUTPUT(HESCOR,1,NCORD,1,NCORD,NCORD,NCORD,1,LUPRI)
         END IF
         DO 100 J = 1, NCORD
            DO 110 I = 1, NCORD
               HESS(I,J) = HESS(I,J) + HESCOR(I,J)
  110       CONTINUE
  100    CONTINUE
         IF (IPRINT .GE. 5) THEN
            CALL HEADER('Corrected Hessian',-1)
            CALL OUTPUT(HESS,1,NCORD,1,NCORD,NCORD,NCORD,1,LUPRI)
         END IF
C
C        Hessian matrix in internal coordinates
C
         CALL DGEMM('N','N',NINTER,NCORD,NCORD,1.D0,
     &              BADJIN,NINTER,
     &              HESS,NCORD,0.D0,
     &              BHES,NINTER)
         CALL DGEMM('N','T',NINTER,NINTER,NCORD,1.D0,
     &              BHES,NINTER,
     &              BADJIN,NINTER,0.D0,
     &              HINT,NINTER)
      END IF
C
C     *************************
C     ***** Print Section *****
C     *************************
C
C     Coordinates, forces and dipole derivatives
C     ==========================================
C
      IF (DIPOL) THEN
         CALL HEADER('Coordinates, Forces, and Dipole Derivatives',1)
C
C        Print in atomic units
C
         WRITE (LUPRI,'(A/)') '  a) atomic units:'
         WRITE (LUPRI,'(A/)')
     &     '  type     coordinate        force '//
     &     '                dipole derivatives'
         DO 200 I = 1, NINTER
            WRITE (LUPRI,'(2X,A4,1X,2(F12.6,4X),7X,3F8.4)')
     &         ITYPCM(I),QQ(I),FINT(I), DINT(I,1), DINT(I,2),DINT(I,3)
  200    CONTINUE
C
C        Print in experimental units
C
         WRITE (LUPRI,'(//A/)') '  b) experimental units:'
         WRITE (LUPRI,'(A/)')
     &     '  type     coordinate        force '//
     &     '                dipole derivatives'
         CALL DSCAL(NINTER,XAJOUL,FINT,1)
         CALL DSCAL(NINTER*3,DEBYE,DINT,1)
         DO 210 I = 1, NINTER
            IF (CTYPE(I) .EQ. D1) THEN
               QQ(I)     = FACANG*QQ(I)
               UNITC = 'deg'
               UNITF = 'mdynA/rad'
               UNITD = 'D/rad'
            ELSE
               QQ(I)     = QQ(I)/CTYPE(I)
               FINT(I)   = FINT(I)*CTYPE(I)
               DINT(I,1) = DINT(I,1)*CTYPE(I)
               DINT(I,2) = DINT(I,2)*CTYPE(I)
               DINT(I,3) = DINT(I,3)*CTYPE(I)
               UNITC = 'A  '
               UNITF = 'mdyn     '
               UNITD = 'D/A  '
            END IF
            WRITE (LUPRI,'(2X,A4,1X,2(F12.6,1X,A),1X,3F8.4,1X,A)')
     &         ITYPCM(I), QQ(I), UNITC, FINT(I), UNITF,
     &         DINT(I,1), DINT(I,2), DINT(I,3), UNITD
  210    CONTINUE
C
C     Coordinates and forces
C     ======================
C
      ELSE
         CALL HEADER('Coordinates and Forces',1)
C
C        Print in atomic units
C
         WRITE (LUPRI,'(A/)') '  a) atomic units:'
         WRITE (LUPRI,'(A/)') '  type      coordinate         force '
         DO 300 I = 1, NINTER
            WRITE (LUPRI,'(2X,A4,4X,2(F12.6,4X))')
     &                ITYPCM(I),QQ(I),FINT(I)
  300    CONTINUE
C
C        Print in experimental units
C
         WRITE (LUPRI,'(//A/)') '  b) experimental units:'
         WRITE (LUPRI,'(A/)') '  type      coordinate         force '
         CALL DSCAL(NINTER,XAJOUL,FINT,1)
         DO 310 I = 1, NINTER
            IF (CTYPE(I) .EQ. D1) THEN
               QQ(I)   = FACANG*QQ(I)
               UNITC = 'deg'
               UNITF = 'mdynA/rad'
            ELSE
               QQ(I)   =   QQ(I)/CTYPE(I)
               FINT(I) = FINT(I)*CTYPE(I)
               UNITC = 'A  '
               UNITF = 'mdyn     '
            END IF
            WRITE (LUPRI,'(2X,A4,4X,2(F12.6,1X,A))')
     &            ITYPCM(I),QQ(I),UNITC,FINT(I),UNITF
  310    CONTINUE
      END IF
C
C
C     Force constants
C     ===============
C
      IF (MAXDIF .EQ. 2) THEN
         CALL HEADER('Force Constants',1)
         WRITE (LUPRI,'(A/)') '  a) atomic units:'
         CALL PRHINT(HINT,NINTER)
         WRITE (LUPRI,'(A/)')
     &     '  b) experimental units (mdyn/A, mdyn/rad, mdynA/rad**2):'
         DO 400 J = 1,NINTER
            DO 410 I = 1,NINTER
               HINT(I,J) = XAJOUL*CTYPE(I)*HINT(I,J)*CTYPE(J)
  410       CONTINUE
  400    CONTINUE
         CALL PRHINT(HINT,NINTER)
      END IF
C
      CALL HEADER('Conversion factors used:',1)
      WRITE (LUPRI,'((A,F15.8))')
     &   '  Bohr to Angstrom :',XTANG,
     &   '  Hartree to aJ    :',XAJOUL
      IF (DIPOL) WRITE (LUPRI,'(A,F15.8)')
     &   '  au to Debye      :',DEBYE
C
      IF (DOVCD) THEN
         CALL HEADER('Atomic axial tensors in internal coordinates',1)
C
C        Print in atomic units
C
         WRITE (LUPRI,'(A/)') '  a) atomic units:'
         DO 500 I = 1, NINTER
            WRITE (LUPRI,'(2X,A4,1X,F12.6,4X,7X,3F8.4)')
     &         ITYPCM(I),QQ(I), AATINT(I,1), AATINT(I,2), AATINT(I,3)
  500    CONTINUE
C
C        Print in experimental units
C
         WRITE (LUPRI,'(//A/)') '  b) experimental units:'
         DO 510 I = 1, NINTER
            WRITE (LUPRI,'(2X,A4,1X,F12.6,4X,7X,3F8.4)')
     &         ITYPCM(I),QQ(I), AATINT(I,1), AATINT(I,2), AATINT(I,3)
  510    CONTINUE
      END IF
      RETURN
      END
C  /* Deck prhint */
      SUBROUTINE PRHINT (MATRIX,NINTER)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER ( KCOL = 6 )
C
      INTEGER BEGIN
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (VAR_REAL)
      REAL MATRIX(NINTER,NINTER)
#else
      DOUBLE PRECISION MATRIX(NINTER,NINTER)
#endif
#include <cbivib.h>
#include <nuclei.h>
C
      NROW = NINTER
      LAST = MIN(NROW,KCOL)
      BEGIN = 1
 1050 NCOL  = 1
      WRITE (LUPRI,1000) (ITYPCM(I),I = BEGIN,LAST)
      WRITE (LUPRI,'()')
      DO 40 K = BEGIN,NROW
         WRITE (LUPRI,2000) ' ',ITYPCM(K),
     &                      (MATRIX(K,(BEGIN-1)+J),J=1,NCOL)
         IF (K .LT. (BEGIN+KCOL-1)) NCOL = NCOL + 1
   40 CONTINUE
      WRITE (LUPRI,'()')
      LAST  = MIN(LAST+KCOL,NROW)
      BEGIN = BEGIN + NCOL
      IF (BEGIN.LE.NROW) GO TO 1050
      WRITE (LUPRI,'()')
      RETURN
 1000 FORMAT (10X,6(4X,A4,4X),(5X,A4,3X))
 2000 FORMAT (A1,1X,A4,1X,6F12.5)
      END
C  /* Deck getb */
      SUBROUTINE GETB(IPRINT,NINTER,NCORD,BMAT,QQ,CTYPE,XA,YA,ZA,BDAG)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <nuclei.h>
#include <symmet.h>
      DIMENSION XA(NCORD), YA(NCORD), ZA(NCORD), QQ(NINTER),
     &          BMAT(NINTER,NCORD), BDAG(NCORD,NINTER), CTYPE(NINTER)
      LOGICAL PRINT, QONLY
#include <ibtfun.h>
C
      IATOM = 0
      DO 100 ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
         DO 110 ISYMOP = 0, MAXOPR
            IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
               IATOM = IATOM + 1
               XA(IATOM) = PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENT)
               YA(IATOM) = PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENT)
               ZA(IATOM) = PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENT)
            END IF
  110    CONTINUE
  100 CONTINUE
      IF (IATOM .NE. NATOMS) THEN
         WRITE (LUPRI,*) 'ERROR in GETB: IATOM .ne. NATOMS'
         WRITE (LUPRI,*) 'IATOM  ',IATOM
         WRITE (LUPRI,*) 'NATOMS ',NATOMS
         CALL QUIT('ERROR in GETB')
      END IF
C
      PRINT = .TRUE.
      QONLY = .FALSE.
      CALL MACHB(NCORD,NINTER,BDAG,XA,YA,ZA,QQ,CTYPE,PRINT,QONLY)
      DO 200 I = 1, NINTER
         DO 210 J = 1, NCORD
            BMAT(I,J) = BDAG(J,I)
  210    CONTINUE
  200 CONTINUE
      RETURN
      END
C  /* Deck getcmat */
      SUBROUTINE GETCMAT(IPRINT,NINTER,NCORD,CMAT,QQ,CTYPE,XA,YA,ZA,
     &                BPOS,BMIN)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (DELTA = 1.D-6, DM1 = -1.D00, ONE = 1.0D0)
C
      DIMENSION XA(NCORD), YA(NCORD), ZA(NCORD), QQ(NINTER),
     &          CMAT(NINTER,NCORD,NCORD), BPOS(NCORD,NINTER),
     &          BMIN(NCORD,NINTER), CTYPE(NINTER)
#include <nuclei.h>
#include <symmet.h>
      LOGICAL PRINT, QONLY
#include <ibtfun.h>
C
      IATOM = 0
      DO 100 ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
         DO 110 ISYMOP = 0, MAXOPR
            IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
               IATOM = IATOM + 1
               XA(IATOM) = PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENT)
               YA(IATOM) = PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENT)
               ZA(IATOM) = PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENT)
            END IF
  110    CONTINUE
  100 CONTINUE
      IF (IATOM .NE. NATOMS) THEN
         WRITE (LUPRI,*) 'ERROR in GETCMAT: IATOM .ne. NATOMS'
         WRITE (LUPRI,*) 'IATOM  ',IATOM
         WRITE (LUPRI,*) 'NATOMS ',NATOMS
         CALL QUIT('ERROR in GETCMAT')
      END IF
C
      PRINT = (IPRINT .GE. 21)
      QONLY = .FALSE.
      DELTAI = ONE/(DELTA + DELTA)
      DO 200 ICOOR = 1, NCORD
         IDATOM = (ICOOR + 2)/3
         IDCOOR = ICOOR - (IDATOM - 1)*3
         IF (IDCOOR .EQ. 1) XA(IDATOM) = XA(IDATOM) + DELTA
         IF (IDCOOR .EQ. 2) YA(IDATOM) = YA(IDATOM) + DELTA
         IF (IDCOOR .EQ. 3) ZA(IDATOM) = ZA(IDATOM) + DELTA
         CALL MACHB(NCORD,NINTER,BPOS,XA,YA,ZA,QQ,CTYPE,PRINT,QONLY)
         IF (IPRINT .GE. 15) THEN
            CALL HEADER('BPOS',-1)
            CALL OUTPUT(BPOS,1,NCORD,1,NINTER,NCORD,NINTER,1,LUPRI)
         END IF
         IF (IDCOOR .EQ. 1) XA(IDATOM) = XA(IDATOM) - DELTA - DELTA
         IF (IDCOOR .EQ. 2) YA(IDATOM) = YA(IDATOM) - DELTA - DELTA
         IF (IDCOOR .EQ. 3) ZA(IDATOM) = ZA(IDATOM) - DELTA - DELTA
         CALL MACHB(NCORD,NINTER,BMIN,XA,YA,ZA,QQ,CTYPE,PRINT,QONLY)
         IF (IPRINT .GE. 15) THEN
            CALL HEADER('BMIN',-1)
            CALL OUTPUT(BMIN,1,NCORD,1,NINTER,NCORD,NINTER,1,LUPRI)
         END IF
         IF (IDCOOR .EQ. 1) XA(IDATOM) = XA(IDATOM) + DELTA
         IF (IDCOOR .EQ. 2) YA(IDATOM) = YA(IDATOM) + DELTA
         IF (IDCOOR .EQ. 3) ZA(IDATOM) = ZA(IDATOM) + DELTA
         CALL DAXPY(NINTER*NCORD,DM1,BMIN,1,BPOS,1)
         DO 210 I = 1, NCORD
            DO 220 J = 1, NINTER
               CMAT(J,I,ICOOR) = DELTAI * BPOS(I,J)
  220       CONTINUE
  210    CONTINUE
         IF (IPRINT .GE. 12) THEN
             WRITE (LUPRI, '(/A,I5)') ' CMAT COORDINATE',ICOOR
             CALL OUTPUT(CMAT(1,1,ICOOR),1,NINTER,1,NCORD,NINTER,NCORD,
     &                   1,LUPRI)
         END IF
  200 CONTINUE
      RETURN
      END
      FUNCTION RPFLIN(C)
C
C 27-May-1986 hjaaj + tuh
C
C Calculate rotational partition function for linear molecules
C with classical partition function 1/C.
C
C Reference: "Statistical Thermodynamics", B.J.McClelland, Chapman
C            and Hall, 1973, p.73ff.
C
#include <implicit.h>
      PARAMETER ( CASYMP = 0.01 D0 , D1 = 1.0D0 )
C
C     This implementation is accurate to at least 7 digits or more.
C
      IF ( C .LE. CASYMP) THEN
         RPFLIN = (D1/C + (D1/3.0D0) + C/15.D0 + C*C*(4.0D0/315.D0))
      ELSE
         FROT = D1
         DO 100 J = 1,60
            DJ   = J
            FJ   = (2.0D0*DJ+1.0D0) * EXP(-(DJ*(DJ+1.0D0)) * C)
            FROT = FROT + FJ
            IF (FJ .LT. 1.D-10) GO TO 200
  100    CONTINUE
  200    CONTINUE
         RPFLIN = FROT
      END IF
      RETURN
      END
C  /* Deck getnos */
      SUBROUTINE GETNOS(ICHRG,GEOM,GRAD,HESS,GRDN,DGRAD,CAAT,WRK,
     &                  LWRK,NOATOM,NCORD,DIPDER,DOVCD,IPRINT,HESPUN)
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      LOGICAL DIPDER, DOVCD,HESPUN
      DIMENSION ICHRG(NOATOM), GEOM(NCORD), GRAD(NCORD), GRDN(NCORD),
     &          HESS(NCORD,NCORD), WRK(LWRK),
     &          DGRAD(3,NCORD), CAAT(3,NCORD)
#include <taymol.h>
#include <energy.h>
#include <moldip.h>
#include <nuclei.h>
#include <symmet.h>
#include <aatens.h>
C
#include <ibtfun.h>
C
C     Atomic numbers and geometry
C
      IATOM = 0
      DO 100 I = 1, NUCIND
         DO 200 J = 0, MAXOPR
            IF (IBTAND(J,ISTBNU(I)) .EQ. 0) THEN
               ICHRG(IATOM + 1) = IZATOM(I)
               DO 300 K = 1, 3
                  GEOM(3*IATOM+K) = PT(IBTAND(ISYMAX(K,1),J))*CORD(K,I)
  300          CONTINUE
               IATOM = IATOM + 1
            END IF
  200    CONTINUE
  100 CONTINUE
      KCSTRA = 1
      KSCTRA = KCSTRA + NCORD*NCORD
      KLAST  = KSCTRA + NCORD*NCORD
      KNEED  = KLAST  + NCORD*NCORD
      IF (KNEED .GT. LWRK) CALL STOPIT('GETNOS','TRACOR',KLAST,LWRK)
      CALL TRACOR(WRK(KCSTRA),WRK(KSCTRA),1,NCORD,0)
C
C     Gradient
C
      CALL DGEMM('N','N',NCORD,1,NCORD,1.D0,
     &           WRK(KSCTRA),NCORD,
     &           GRDMOL,MXCOOR,0.D0,
     &           GRAD,NCORD)
C
C     Nuclear Repulsion Gradient
C
      CALL DGEMM('N','N',NCORD,1,NCORD,1.D0,
     &           WRK(KSCTRA),NCORD,
     &           GRADNN,MXCOOR,0.D0,
     &           GRDN,NCORD)
C
C     Hessian
C
      CALL DGEMM('N','N',NCORD,NCORD,NCORD,1.D0,
     &           WRK(KSCTRA),NCORD,
     &           HESMOL,MXCOOR,0.D0,
     &           WRK(KLAST),MXCOOR)
      CALL DGEMM('N','T',NCORD,NCORD,NCORD,1.D0,
     &           WRK(KLAST),MXCOOR,
     &           WRK(KSCTRA),NCORD,0.D0,
     &           HESS,NCORD)
C
C     Dipole gradient
C
      IF (DIPDER) CALL TRADIP(DIPFLT,DGRAD,WRK(KCSTRA),WRK(KSCTRA),
     &                        NCORD,1,1)
C
C     AAT
C
      IF (DOVCD) CALL TRADIP(AATTOT,CAAT,WRK(KCSTRA),WRK(KSCTRA),
     &                       NCORD,2,1)
C
C     Print Section
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Charges in non-symmetry coordinates',-1)
         WRITE (LUPRI,'(12I5)') (ICHRG(I),I=1,NOATOM)
         CALL HEADER('Geometry in non-symmetry coordinates',-1)
         CALL OUTPUT(GEOM,1,1,1,NCORD,1,NCORD,1,LUPRI)
         CALL HEADER('Gradient in non-symmetry coordinates',-1)
         CALL OUTPUT(GRAD,1,1,1,NCORD,1,NCORD,1,LUPRI)
         CALL HEADER('Nuc. Rep. gradient in non-sym. coordinates',-1)
         CALL OUTPUT(GRDN,1,1,1,NCORD,1,NCORD,1,LUPRI)
         CALL HEADER('Hessian in non-symmetry coordinates',-1)
         CALL OUTPUT(HESS,1,NCORD,1,NCORD,NCORD,NCORD,1,LUPRI)
         IF (DIPDER) THEN
            CALL HEADER
     &         ('Dipole gradient in non-symmetry basis - GETNOS',-1)
            CALL OUTPUT(DGRAD,1,3,1,NCORD,3,NCORD,1,LUPRI)
         END IF
         IF (DOVCD) THEN
            CALL HEADER('AAT in non-symmetry basis - GETNOS',-1)
            CALL OUTPUT(CAAT,1,3,1,NCORD,3,NCORD,1,LUPRI)
         END IF
      END IF
C
C     If we have asked that the Hessian be punched, so be it
C
      LUHES = -1
      IF (HESPUN) THEN
         CALL GPOPEN(LUHES,'DALTON.HES','UNKNOWN',' ','FORMATTED',
     &               IDUMMY,.FALSE.)
         REWIND (LUHES)
         WRITE (LUHES,*) NCORD
         WRITE (LUHES,*)
         DO 101 ICOOR = 1, NCORD
            DO 102 JCOOR = 1, NCORD
               WRITE(LUHES,*) HESS(JCOOR,ICOOR)
 102        CONTINUE
            WRITE (LUHES,*)
 101     CONTINUE
         DO 103 ICOOR = 1, NCORD
            WRITE (LUHES,*) GEOM(ICOOR)
 103     CONTINUE
         CALL GPCLOSE(LUHES,'KEEP')
      END IF
      RETURN
      END
C  /* Deck vibv3 */
      SUBROUTINE VIBV3(EVEC,FREQAU,FMATF,FMATT,CORR,
     &                 POS,AMASS,WORK,LWORK,NCORD,NUMMOD)
#include <implicit.h>
#include <dummy.h>
#include <thrzer.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <codata.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0,
     &           D4 = 4.0D0, D24 = 24.0D0, D36 = 36.0D0)
      DIMENSION EVEC(NCORD,NCORD), FREQAU(NCORD),
     &          FMATF(NCORD,NCORD,NCORD), CORR(NUMMOD),
     &          FMATT(NCORD,NCORD,NCORD),
     &          POS(3,NCORD/3), AMASS(*), WORK(LWORK)
C
#include <symmet.h>
#include <nuclei.h>
#include <abainf.h>
#include <inftap.h>
#include <cbivib.h>
#include <cbinum.h>
#include <cbiwlk.h>
#include <gnrinf.h>
#include <ibtfun.h>
C
      CALL GPOPEN(LUWLK,ABAWLK,'OLD','SEQUENTIAL','UNFORMATTED',IDUMMY,
     &          .FALSE.)
      REWIND (LUWLK)
      READ (LUWLK)
      READ (LUWLK)
      READ (LUWLK)
      READ (LUWLK)
      READ (LUWLK) IUMMOD, FMATF
      CALL GPCLOSE(LUWLK,'KEEP')
C
C     Transform V3 matrix to normal coordinates
C
      IF (.NOT. NMODIF) THEN
      DO I = 1, NCORD
         DO J = 1, NCORD
            DO K =1, NCORD
               TMP = D0
               DO KP = 1, NCORD
                  TMP = TMP + FMATF(I,J,KP)*EVEC(KP,K)
               END DO
               FMATT(I,J,K) = TMP
            END DO
         END DO
      END DO
      DO I = 1, NCORD
         DO J = 1, NCORD
            DO K =1, NCORD
               TMP = D0
               DO JP = 1, NCORD
                  TMP = TMP + FMATT(I,JP,K)*EVEC(JP,J)
               END DO
               FMATF(I,J,K) = TMP
            END DO
         END DO
      END DO
      DO I = 1, NUMMOD
         DO J = 1, NUMMOD
            DO K = 1, NUMMOD
               TMP = D0
               DO IP = 1, NCORD
                  TMP = TMP + FMATF(IP,J,K)*EVEC(IP,I)
               END DO
               FMATT(I,J,K) = TMP
            END DO
         END DO
      END DO
      ELSE
         CALL DCOPY(NCORD*NCORD*NCORD,FMATF,1,FMATT,1)
         IF (.NOT.NUMVIB) THEN
            FAC = D1/XFAMU
            CALL DSCAL(NCORD*NCORD*NCORD,FAC,FMATT,1)
         END IF
      END IF

C
C     Print in normal-coordinates
C
      CALL TITLER('Output from Anharmonic Force Constant analysis',
     &            '*',103)
      IF (IPRINT .GT. 2) THEN
         CALL HEADER('Normal coordinates',1)
         CALL OUTPUT(EVEC,1,NCORD,1,NUMMOD,NCORD,NCORD,1,LUPRI)
C
         DO I = 1, NUMMOD
            WRITE (LUPRI,'(/A,I5,/A)')
     &      ' Anharmonic force constants for normal mode(*,*,I):',I,
     &      ' -------------------------------------------------------'
            CALL OUTPUT(FMATT(1,1,I),1,NUMMOD,1,NUMMOD,NCORD,NCORD,
     &                  1,LUPRI)
         END DO
      END IF
C
C     Calculate correction to equilibrium geometry in normal coordinates
C
      IF (DOTEMP) THEN
         DO ITMP = 1, NTEMP
C
C     In case centrifugal distortions are to be included, calculate
C     contributions to the effective geometry here
C
            WRITE (LUPRI,'(////A,/1X,A,F12.6,A,/A)')
     &           '----------------------------'//
     &           '------------------------------',
     &           'Vibrational corrections for temperature T = ',
     &           TEMP(ITMP),' K','----------------------------'//
     &           '------------------------------'
            IF (DOCENT) THEN
               KCCOOR = 1
               KCCORN = KCCOOR + MXCOOR
               KQCENT = KCCORN + MXCOOR*MXCOOR
               KLAST  = KQCENT + MXCOOR
               IF (KLAST .GT. LWORK) CALL STOPIT('VIBV3 ','CENTRIF',
     &                                           KLAST,LWORK)
               LWRK   = LWORK - KLAST + 1
               CALL CENTRIF(NCORD/3,NUMMOD,CORD,AMASS,WORK(KCCOOR),
     &           EVEC,WORK(KCCORN),FREQAU,WORK(KQCENT),
     &           WORK(KLAST),LWRK,TEMP(ITMP),IPRINT)
               CALL HEADER('Centrifugal corrections to effective '//
     &                     'geometries in normal coordinates',1)
               DO IMOD = 1, NUMMOD
                  WRITE (LUPRI,'(A,I5,A,F10.4,A,F10.4)') 'Normal mode:',
     &                 IMOD,'   Frequency:', FREQAU(IMOD),
     &                 '   Displacement:',WORK(KQCENT + IMOD - 1)
               END DO
            END IF
C
C     Vibrational effects now being calculated
C
            CALL HEADER('Effective geometry corrections in normal '//
     &                  'coordinates (vibrational) ',1)
            DO IMOD = 1, NUMMOD
               PFAC = FREQAU(IMOD)*FREQAU(IMOD)*D4*SQRT(XFAMU)
               PFACI = D1/PFAC
               TMP = D0
               DO MMOD = 1, NUMMOD
                  IF (TEMP(ITMP) .LT. THRZER) THEN
                     TFAC2 = D1
                  ELSE
                     TFAC  = FREQAU(MMOD)*AUTK/(D2*TEMP(ITMP))
                     TFAC2 = D1/DTANH(TFAC)
                  END IF
                  TMP = TMP + (FMATT(MMOD,MMOD,IMOD)/FREQAU(MMOD))*TFAC2
               END DO
               CORR(IMOD) = - TMP*PFACI
               WRITE (LUPRI,'(A,I5,A,F10.4,A,F10.4)') 'Normal mode:',
     &              IMOD,'   Frequency:', FREQAU(IMOD),
     &              '   Displacement:',CORR(IMOD)
            END DO
            IF (DOCENT) THEN
               CALL HEADER('Total effective geometry corrections '//
     &              'in normal coordinates ',1)
               DO IMOD = 1, NUMMOD
                  CORR(IMOD) = CORR(IMOD) + WORK(KQCENT + IMOD - 1)
                  WRITE (LUPRI,'(A,I5,A,F10.4,A,F10.4)') 'Normal mode:',
     &                 IMOD,'   Frequency:', FREQAU(IMOD),
     &                 '   Displacement:',CORR(IMOD)
               END DO
            END IF
C
C     Convert the suggested displacement into Cartesian coordinates
C
            CALL DZERO(POS,NCORD)
            FAC = SQRT(XFAMU)
            DO IMOD = 1, NUMMOD
               ICOR = 0
               DO IATOM = 1, NCORD/3
                  DO ICOOR = 1, 3
                     ICOR = ICOR + 1
                     POS(ICOOR,IATOM) = POS(ICOOR,IATOM)
     &                                + CORR(IMOD)*EVEC(ICOR,IMOD)*FAC
                  END DO
               END DO
            END DO
CKR
CKR      Print out temperature somewhere
CKR
            CALL AROUND('Change to effective geometry')
            CALL HEADER('Cartesian coordinates (au)',-1)
            CALL PRIGEO(POS)
C
            DO IATOM = 1, NCORD/3
               DO ICOOR = 1, 3
                  POS(ICOOR,IATOM) = CORD(ICOOR,IATOM) +POS(ICOOR,IATOM)
               END DO
            END DO
C
            CALL AROUND('Estimated effective geometry')
            CALL HEADER('Cartesian coordinates (au)',-1)
            CALL PRIGEO(POS)
            CALL GEOANA(POS,.TRUE.,.FALSE.,NBONDS,-1,WORK,LWORK)
         END DO
      ELSE
         CALL HEADER('Effective geometry corrections in normal '//
     &               'coordinates',1)
         DO IMOD = 1, NUMMOD
            PFAC = FREQAU(IMOD)*FREQAU(IMOD)*D4*SQRT(XFAMU)
            PFACI = D1/PFAC
            TMP = D0
            DO MMOD = 1, NUMMOD
               TMP = TMP + FMATT(MMOD,MMOD,IMOD)/FREQAU(MMOD)
            END DO
            CORR(IMOD) = - TMP*PFACI
            WRITE (LUPRI,'(A,I5,A,F10.4,A,F10.4)') 'Normal mode:',IMOD,
     &           '   Frequency:', FREQAU(IMOD), '   Displacement:',
     &           CORR(IMOD)
         END DO
C
C     Convert the suggested displacement into Cartesian coordinates
C
         CALL DZERO(POS,NCORD)
         FAC = SQRT(XFAMU)
         DO IMOD = 1, NUMMOD
            ICOR = 0
            DO IATOM = 1, NCORD/3
               DO ICOOR = 1, 3
                  ICOR = ICOR + 1
                  POS(ICOOR,IATOM) = POS(ICOOR,IATOM)
     &                             + CORR(IMOD)*EVEC(ICOR,IMOD)*FAC
               END DO
            END DO
         END DO
         CALL AROUND('Change to effective geometry')
         CALL HEADER('Cartesian coordinates (au)',-1)
         IF (NMWALK) THEN
            DO IATOM = 1, NUCDEP
               WRITE (LUPRI,'(1X,A,F17.10,2F24.10)') NAMDEP(IATOM),
     *                                  (POS(ICOOR,IATOM),ICOOR=1,3)
            END DO
         ELSE
            CALL PRIGEO(POS)
         END IF
C
         ICOOR = 0
         IATOM = 0
         DO ICENT = 1, NUCIND
            MULCNT = ISTBNU(ICENT)
            DO IOP = 0, MAXOPR
               IF (IBTAND(IOP,MULCNT) .EQ. 0) THEN
                  IATOM = IATOM + 1
                  DO I = 1, 3
                     ICOOR = ICOOR + 1
                     POS(I,IATOM) = POS(I,IATOM)
     &                    + PT(IBTAND(ISYMAX(I,1),IOP))*CORD(I,ICENT)
                  END DO
               END IF
            END DO
         END DO
C
         CALL AROUND('Estimated effective geometry')
         CALL HEADER('Cartesian coordinates (au)',-1)
         IF (NMWALK) THEN
            DO IATOM = 1, NUCDEP
               WRITE (LUPRI,'(1X,A,F17.10,2F24.10)') NAMDEP(IATOM),
     *                                  (POS(ICOOR,IATOM),ICOOR=1,3)
            END DO
         ELSE
            CALL PRIGEO(POS)
            CALL GEOANA(POS,.TRUE.,.FALSE.,NBONDS,-1,WORK,LWORK)
         END IF
      END IF
C
C     Some analysis of the quality of the above predictions
C
      IF (IPRINT .GE. 1) THEN
         CALL AROUND('Perturbation expansion analysis')
         CALL HEADER ('Linear wave function parameters',-1)
         WRITE (LUPRI,'(/A,/A)')' Mode       a^(1)_1          a^(1)_3',
     &                          ' -------------------------------------'
         DO IMOD = 1, NUMMOD
            PFAC = D4*SQRT(D2)*SQRT(FREQAU(IMOD)**3)
            A11 = D0
            DO MMOD = 1, NUMMOD
               A11 = A11 + FMATT(MMOD,MMOD,IMOD)/FREQAU(MMOD)
            END DO
            A11 = - A11/PFAC
            A13 = -SQRT(D3)*FMATT(IMOD,IMOD,IMOD)
     &                    /(D36*SQRT(FREQAU(IMOD)**5))
            WRITE (LUPRI,'(2X,I3,5X,F11.9,8X,F11.9)') IMOD,A11,A13
         END DO
C
         DO IMOD = 1, NUMMOD
            DO JMOD = 1, NUMMOD
               FAC = D4*FREQAU(IMOD)*SQRT(FREQAU(JMOD))
     &                              *(D2*FREQAU(IMOD) + FREQAU(JMOD))
               FMATF(IMOD,JMOD,1) = - FMATT(IMOD,IMOD,JMOD)/FAC
            END DO
            FMATF(IMOD,IMOD,1) = D0
         END DO
         CALL HEADER('Bilinear wave function parameters',-1)
         CALL OUTPUT(FMATF(1,1,1),1,NUMMOD,1,NUMMOD,NCORD,NCORD,1,LUPRI)
      END IF
C
      IF (IPRINT .GT. 1 .AND. .NOT. NMODIF) THEN
         DO IMOD = 1, NUMMOD
            DO JMOD = 1, NUMMOD
               DO KMOD = 1, NUMMOD
                  FAC = D24*SQRT(FREQAU(IMOD)*FREQAU(JMOD)*FREQAU(KMOD))
     &                 *(FREQAU(IMOD) + FREQAU(JMOD) + FREQAU(KMOD))
                  FMATF(IMOD,JMOD,KMOD) = FMATT(IMOD,JMOD,KMOD)
     &                                   *SQRT(D2)/FAC
               END DO
            END DO
         END DO
C
C     Not all matrix elements are meaningful
C
         DO IMOD = 1, NUMMOD
            FMATF(IMOD,IMOD,IMOD) = D0
            DO JMOD = 1, NUMMOD
               FMATF(IMOD,IMOD,JMOD) = D0
               FMATF(IMOD,JMOD,IMOD) = D0
               FMATF(JMOD,IMOD,IMOD) = D0
            END DO
         END DO
         DO I = 1, NUMMOD
            WRITE (LUPRI,'(/A,I5,/A)')
     &      ' Trilinear wave function parameters for mode (*,*,I):',I,
     &      ' -------------------------------------------------------'
            CALL OUTPUT(FMATF(1,1,I),1,NUMMOD,1,NUMMOD,NCORD,NCORD,
     &                  1,LUPRI)
         END DO
      END IF
C
      RETURN
      END
C  /* Deck vibroa */
      SUBROUTINE VIBROA(EVEC,FREQ,ROAAFQ,ROAGNQ,ROAGLQ,ROAAQ,WORK,LWORK,
     &                  NCORD,NUMMOD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
      DIMENSION EVEC(NCORD,NCORD), FREQ(NCORD),
     &          ROAAFQ(3,3,NFRVAL,NCORD), ROAGNQ(3,3,NFRVAL,NCORD),
     &          ROAGLQ(3,3,NFRVAL,NCORD), ROAAQ(3,3,3,NFRVAL,NCORD)
      DIMENSION WORK(LWORK)
#include <cbilnr.h>
      KBAL   = 1
      KBGML  = KBAL   + NCORD
      KBGMNL = KBGML  + NCORD
      KBA    = KBGMNL + NCORD
      KRMINL = KBA    + NCORD
      KDPOLL = KRMINL + NCORD
      KRMINN = KDPOLL + NCORD
      KDPOLN = KRMINN + NCORD
      KDZL   = KDPOLN + NCORD
      KDZNL  = KDZL   + NCORD
      KDXL   = KDZNL  + NCORD
      KDXNL  = KDXL   + NCORD
      KD0L   = KDXNL  + NCORD
      KD0NL  = KD0L   + NCORD
      KDBL   = KD0NL  + NCORD
      KDBNL  = KDBL   + NCORD
      KALMN  = KDBNL  + NCORD
      KGMMNL = KALMN  + NCORD
      KGMMNN = KGMMNL + NCORD
      KCZL   = KGMMNN + NCORD
      KCZNL  = KCZL   + NCORD
      KCXL   = KCZNL  + NCORD
      KCXNL  = KCXL   + NCORD
      KC0L   = KCXNL  + NCORD
      KC0NL  = KC0L   + NCORD
      KCBL   = KC0NL  + NCORD
      KCBNL  = KCBL   + NCORD
      KCONVE = KCBNL  + NCORD
      KBOLTZ = KCONVE + NCORD
      KROAAF = KBOLTZ + NCORD
      KROAGN = KROAAF + 9*MXFR*MXCOOR
      KROAGL = KROAGN + 9*MXFR*MXCOOR
      KROAAD = KROAGL + 9*MXFR*MXCOOR
      KROAAU = KROAAD + 27*MXFR*MXCOOR
      KROAGU = KROAAU + 9*MXFR
      KROGU  = KROAGU + 9*MXFR
      KROAU  = KROGU + 9*MXFR
      KLAST  = KROAU  + 27*MXFR
      IF (KLAST .GT. LWORK) CALL STOPIT('VIBROA',' ',KLAST,LWORK)
      CALL VIBRO1(EVEC,FREQ,ROAAFQ,ROAGNQ,ROAGLQ,ROAAQ,
     &            WORK(KBAL),WORK(KBGML),WORK(KBGMNL),WORK(KBA),
     &            WORK(KRMINL),WORK(KDPOLL),WORK(KDZL),
     &            WORK(KDZNL),WORK(KDXL),WORK(KDXNL),WORK(KD0L),
     &            WORK(KD0NL),WORK(KDBL),WORK(KDBNL),NCORD,NUMMOD,
     &            WORK(KALMN),WORK(KGMMNL),WORK(KGMMNN),
     &            WORK(KRMINN),WORK(KDPOLN),WORK(KCZL),WORK(KCZNL),
     &            WORK(KCXL),WORK(KCXNL),WORK(KC0L),WORK(KC0NL),
     &            WORK(KCBL),WORK(KCBNL),WORK(KCONVE),WORK(KBOLTZ),
     &            WORK(KROAAF),WORK(KROAGN),WORK(KROAGL),WORK(KROAAD),
     &            WORK(KROAAU),WORK(KROAGU),WORK(KROGU),WORK(KROAU))
      RETURN
      END
C  /* Deck vibro1 */
      SUBROUTINE VIBRO1(EVEC,FREQ,ROAAFQ,ROAGNQ,ROAGLQ,ROAAQ,
     &                  BAL,BGML,BGMNL,BA,RMINL,DPOLL,DZL,DZNL,
     &                  DXL,DXNL,D0L,D0NL,DBL,DBNL,NCORD,NUMMOD,
     &                  ALMN,GMMNL,GMMNNL,RMINN,DPOLN,CZL,CZNL,CXL,CXNL,
     &                  C0L,C0NL,CBL,CBNL,CONVER,BOLTZ,ROAAFD,ROAGND,
     &                  ROAGLD,ROAAD,ROAAFU,ROAGNU,ROAGLU,ROAAU)
#include <implicit.h>
#include <dummy.h>
C
C     Transforms tensor gradients from Cartesian to normal coordinates
C     for Raman vibrational intensities and VROA
C     Calculation of the final properties and output modified
C     by G.Hangartner 20.11.1996
C
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <codata.h>
#include <cbilnr.h>
C
      INTEGER F
      PARAMETER (DE3=1.0D3,DE13 = 1.0D13, D2=2.0D0, DE6=1.0D6 )
C Conversion factors from au to derivated cgs
      PARAMETER (XFOAIN = ALPHAC*XTANG**4*XFAMU*DE6,
     &           XFRAIN = XTANG**4*XFAMU,
     &           XFCONV = 1/(XTANG**2*XFAMU)*DE13)
      DIMENSION EVEC(NCORD,NCORD),
     &          ROAAFQ(3,3,NFRVAL,NCORD), ROAGNQ(3,3,NFRVAL,NCORD),
     &          ROAGLQ(3,3,NFRVAL,NCORD), ROAAQ(3,3,3,NFRVAL,NCORD),
     &          BAL(NCORD), BGML(NCORD), BGMNL(NCORD), BA(NCORD),
     &          RMINL(NCORD), ALMN(NCORD), GMMNL(NCORD), GMMNNL(NCORD),
     &          DPOLL(NCORD),RMINN(NCORD),DPOLN(NCORD), DZL(NCORD),
     &          DZNL(NCORD), FREQ(NCORD), DXL(NCORD),
     &          DXNL(NCORD), D0L(NCORD), D0NL(NCORD), DBL(NCORD),
     &          DBNL(NCORD), CZL(NCORD), CZNL(NCORD), CXL(NCORD),
     &          CXNL(NCORD), C0L(NCORD), C0NL(NCORD), CBL(NCORD),
     &          CBNL(NCORD), CONVER(NCORD),BOLTZ(NCORD),
     &          ROAAFD(3,3,MXFR,MXCOOR), ROAGND(3,3,MXFR,MXCOOR),
     &          ROAGLD(3,3,MXFR,MXCOOR), ROAAD(3,3,3,MXFR,MXCOOR),
     &          ROAAFU(3,3,MXFR), ROAGNU(3,3,MXFR),
     &          ROAGLU(3,3,MXFR), ROAAU(3,3,3,MXFR)
C
#include <abainf.h>
#include <inftap.h>
#include <cbivib.h>
C
C     Read in Cartesian gradients
C
      CALL GPOPEN(LUWLK,ABAWLK,'OLD','SEQUENTIAL','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUWLK
      READ (LUWLK)
      READ (LUWLK)
      READ (LUWLK)
      READ (LUWLK) NFRVAL, FRVAL,
     &             ROAAFU,ROAAFD,ROAGNU,ROAGND,
     &             ROAGLU,ROAGLD, ROAAU, ROAAD
      CALL GPCLOSE(LUWLK,'KEEP')
C
C     Transform to normal coordinates
C
      CALL DGEMM('N','N',9*NFRVAL,NCORD,NCORD,1.D0,
     &           ROAAFD,9*MXFR,
     &           EVEC,NCORD,0.D0,
     &           ROAAFQ,9*NFRVAL)
      CALL DGEMM('N','N',9*NFRVAL,NCORD,NCORD,1.D0,
     &           ROAGND,9*MXFR,
     &           EVEC,NCORD,0.D0,
     &           ROAGNQ,9*NFRVAL)
      CALL DGEMM('N','N',9*NFRVAL,NCORD,NCORD,1.D0,
     &           ROAGLD,9*MXFR,
     &           EVEC,NCORD,0.D0,
     &           ROAGLQ,9*NFRVAL)
      CALL DGEMM('N','N',27*NFRVAL,NCORD,NCORD,1.D0,
     &           ROAAD,27*MXFR,
     &           EVEC,NCORD,0.D0,
     &           ROAAQ,27*NFRVAL)
C
      IF (IPRINT .GT. 5) THEN
         CALL TITLER('Output from VIBRO1','*',103)
         CALL HEADER('Polarizability tensor Cartesian gradient)',1)
         DO 100 F = 1, NFRVAL
            WRITE (LUPRI,'(1X,A,F6.4)') ' frequency:', FRVAL(F)
            CALL OUTPUT(ROAAFD(1,1,F,1),1,9,1,NUMMOD,9*MXFR,NCORD,
     &                  1,LUPRI)
  100    CONTINUE
         IF (VROA) THEN
            CALL HEADER('Non-London G tensor Cartesian gradient',1)
            DO 110 F = 1, NFRVAL
               WRITE (LUPRI,'(1X,A,F6.4)') ' frequency:', FRVAL(F)
               CALL OUTPUT(ROAGND(1,1,F,1),1,9,1,NUMMOD,9*MXFR,NCORD,
     &                     1,LUPRI)
  110       CONTINUE
            CALL HEADER('London G tensor Cartesian gradient',1)
            DO 120 F = 1, NFRVAL
               WRITE (LUPRI,'(1X,A,F6.4)') ' frequency:', FRVAL(F)
               CALL OUTPUT(ROAGLD(1,1,F,1),1,9,1,NUMMOD,9*MXFR,NCORD,
     &                     1,LUPRI)
  120       CONTINUE
            CALL HEADER('A tensor Cartesian gradient',1)
            DO 130 F = 1, NFRVAL
               WRITE (LUPRI,'(1X,A,F6.4)') ' frequency:', FRVAL(F)
               CALL OUTPUT(ROAAD(1,1,1,F,1),1,27,1,NUMMOD,27*MXFR,NCORD,
     &                     1,LUPRI)
  130       CONTINUE
C
         END IF
         CALL HEADER('Normal coordinates',1)
         CALL OUTPUT(EVEC,1,NCORD,1,NUMMOD,NCORD,NCORD,1,LUPRI)
C
         CALL HEADER('Polarizability tensor normal coord. gradient',1)
         DO 200 F = 1, NFRVAL
            WRITE (LUPRI,'(1X,A,F6.4)') ' frequency:', FRVAL(F)
            CALL OUTPUT(ROAAFQ(1,1,F,1),1,9,1,NUMMOD,9*NFRVAL,NCORD,
     &                  1,LUPRI)
  200    CONTINUE
         IF (VROA) THEN
           CALL HEADER('Non-London G tensor normal coordinate gradient',
     &                 1)
            DO 210 F = 1, NFRVAL
               WRITE (LUPRI,'(1X,A,F6.4)') ' frequency:', FRVAL(F)
               CALL OUTPUT(ROAGNQ(1,1,F,1),1,9,1,NUMMOD,9*NFRVAL,NCORD,
     &                     1,LUPRI)
  210       CONTINUE
            CALL HEADER('London G tensor normal coordinate gradient',1)
            DO 220 F = 1, NFRVAL
               WRITE (LUPRI,'(1X,A,F6.4)') ' frequency:', FRVAL(F)
               CALL OUTPUT(ROAGLQ(1,1,F,1),1,9,1,NUMMOD,9*NFRVAL,NCORD,
     &                     1,LUPRI)
  220       CONTINUE
            CALL HEADER('A tensor normal coordinate gradient',1)
            DO 230 F = 1, NFRVAL
               WRITE (LUPRI,'(1X,A,F6.4)') ' frequency:', FRVAL(F)
               CALL OUTPUT(ROAAQ(1,1,1,F,1),1,27,1,NUMMOD,27*NFRVAL,
     &                     NCORD,1,LUPRI)
  230       CONTINUE
         END IF
      END IF
C
C Loop over all freq. Calculate properities used in both, Raman and VROA
C
      DO 500 F = 1,NFRVAL
        DO 520 I = 1,NUMMOD
              ALMN  (I) = ALFMN  (ROAAFQ(1,1,F,I))
              BAL   (I) = BETAAL (ROAAFQ(1,1,F,I))
              CONVER(I) = CROSSK (FRVAL(F),FREQ(I))
              BOLTZ (I) = BOLTZK (FREQ(I))
 520    CONTINUE
C
C Calculate specific Raman properities and print them
C
        IF (RAMAN) THEN
           WRITE(LUPRI,'(/,5X,A,F9.6,A,F9.2,A)')
     &          'Raman related properties for freq. ',FRVAL(F),' au  =',
     &           XTNM/FRVAL(F),' nm'
           WRITE(LUPRI,'(5X,A)')
     &     '-----------------------------------------------'//
     &      '----------------'
           WRITE(LUPRI,'(/,1X,A,/)')
     &          'Mode    Freq.     Alpha**2   Beta(a)**2'//
     &          '   Pol.Int.   Depol.Int.  Dep. Ratio '
C
           DO 550 I = 1,NUMMOD
              RMINL (I) = RAMINL (ALMN(I),BAL(I))
              DPOLL (I) = DEPOLL (ALMN(I),BAL(I))
              WRITE (LUPRI,'(1X,I4,2X,F8.2,5F12.6)') I, XTKAYS*FREQ(I),
     &            ALMN(I)*ALMN(I)*XFRAIN,BAL(I)*XFRAIN,
     &            RMINL(I)*XFRAIN,DPOLL(I)*RMINL(I)*XFRAIN,DPOLL(I)
 550       CONTINUE
           WRITE(LUPRI,'(/,1X,A,/)')
     &          'Mode    Freq.   ConverFac(*E13)  BoltzFac'
C
           DO 560 I = 1,NUMMOD
              WRITE (LUPRI,'(1X,I4,2X,F8.2,F12.6,F14.8)')
     &            I,XTKAYS*FREQ(I),CONVER(I)*XFCONV,BOLTZ(I)
 560       CONTINUE
        END IF
C
C End of Raman Output for one frequency, going VROA or next freq.
C
C Calculate VROA properities and print them
C
        IF (VROA) THEN
           WRITE(LUPRI,'(/,2X,A,F9.6,A,F9.2,A)')
     &        'Raman Optical Activity properties for freq.  ',
     &         FRVAL(F),' au  =',XTNM/FRVAL(F),' nm'
           WRITE(LUPRI,'(2X,A)')
     &    '---------------------------------------'//
     &    '---------------------------------'
C
        DO 600 I = 1, NUMMOD
          BGML  (I) = BETAGM(ROAAFQ(1,1,F,I),ROAGLQ(1,1,F,I))
          BGMNL (I) = BETAGM(ROAAFQ(1,1,F,I),ROAGNQ(1,1,F,I))
          BA    (I) = BETAA (ROAAFQ(1,1,F,I),ROAAQ(1,1,1,F,I),FRVAL(F))
          GMMNL (I) = GMMN  (ROAGLQ(1,1,F,I))
          GMMNNL(I) = GMMN  (ROAGNQ(1,1,F,I))
          RMINN (I) = RAMINN (ALMN(I),BAL(I))
          DPOLN (I) = DEPOLN (ALMN(I),BAL(I))
          DZL   (I) = DELTAZ(BGML(I),BA(I))
          DZNL  (I) = DELTAZ(BGMNL(I),BA(I))
          DXL   (I) = DELTAX(BGML(I),BA(I),ALMN(I),GMMNL(I))
          DXNL  (I) = DELTAX(BGMNL(I),BA(I),ALMN(I),GMMNNL(I))
          D0L   (I) = DELTA0(BGML(I),BA(I),ALMN(I),GMMNL(I))
          D0NL  (I) = DELTA0(BGMNL(I),BA(I),ALMN(I),GMMNNL(I))
          DBL   (I) = DELTAB(BGML(I),BA(I))
          DBNL  (I) = DELTAB(BGMNL(I),BA(I))
          CZL   (I) = CID(DZL(I),RMINN(I)*DPOLN(I))
          CZNL  (I) = CID(DZNL(I),RMINN(I)*DPOLN(I))
          CXL   (I) = CID(DXL(I),RMINN(I))
          CXNL  (I) = CID(DXNL(I),RMINN(I))
          C0L   (I) = CID(D0L(I),D2*RMINN(I))
          C0NL  (I) = CID(D0NL(I),D2*RMINN(I))
          CBL   (I) = CID(DBL(I),D2*RMINN(I))
          CBNL  (I) = CID(DBNL(I),D2*RMINN(I))
  600   CONTINUE
C
        WRITE (LUPRI,'(/,3X,A)')
     &      '**** Raman invariants and Intensity Parameters'//
     &      ' (A**4amu**-1) ****'
        WRITE (LUPRI,'(/,1X,A,/)')
     &      'Mode  Freq.     Alpha**2   Beta(a)**2'//
     &       '   Pol.Int. DepolInt. 180/0Int.  DepRatio'
        DO 700 I = 1, NUMMOD
            WRITE(LUPRI,'(I4,F10.2,2F12.6,4F10.4)')
     &            I, XTKAYS*FREQ(I),ALMN(I)**2*XFRAIN,BAL(I)*XFRAIN,
     &            RMINN(I)*XFRAIN,DPOLN(I)*RMINN(I)*XFRAIN,
     &            D2*RMINN(I)*XFRAIN,DPOLN(I)
  700   CONTINUE
        WRITE (LUPRI,'(/,3X,A)')
     &      '**** ROA invariants (*E6, in A**4amu**-1) ****'
        WRITE (LUPRI,'(/,1X,A,/)')
     &      'Mode   Freq.     a*Gm(Lon)  a*Gm(noL)'//
     &      '   Beta(Gm)**2 (L/noL)    Beta(A)**2'
        DO 705 I = 1, NUMMOD
          WRITE (LUPRI,'(I4,1X,F8.2,2F12.4,3F12.2)') I, XTKAYS*FREQ(I),
     &           ALMN(I)*GMMNL(I)*XFOAIN,ALMN(I)*GMMNNL(I)*XFOAIN,
     &           BGML(I)*XFOAIN,BGMNL(I)*XFOAIN,BA(I)*XFOAIN
 705    CONTINUE
        WRITE (LUPRI,'(/,3X,A)')
     &      '**** Conversion Factor to Cross Section (*E13)'//
     &      ' and Boltzmann Factor'
        WRITE (LUPRI,'(/,1X,A,/)')
     &      'Mode  Freq.     ConverFac(*E13)   BoltzFac'
        DO 710 I = 1, NUMMOD
            WRITE(LUPRI,'(I4,F10.2,F14.6,F16.8)')
     &              I, XTKAYS*FREQ(I),CONVER(I)*XFCONV,BOLTZ(I)
  710   CONTINUE
C
C Difference Parameter = 4/c * (45aGm+7beta(Gm)+beta(A)  etc.
C
        WRITE (LUPRI,'(/,3X,A,//,2X,A,/)') '**** Difference '//
     &        'Parameter R-L (A**4amu**-1) (London orbitals) ****',
     &        'Mode   Frequency     DELTApar      DELTAperp  '//
     &         '    DELTA0       DELTA180'
        DO 750 I = 1, NUMMOD
          WRITE(LUPRI,'(1X,I4,2X,F10.2,4F14.6)')I,XTKAYS*FREQ(I),
     &     DZL(I)*XFRAIN,DXL(I)*XFRAIN,D0L(I)*XFRAIN,DBL(I)*XFRAIN
  750   CONTINUE
C
C  Here are printed the CIDs (Definition of Barron):
C  CID = (R-L) / (2*Intensity)
C
           WRITE (LUPRI,'(/,3X,A,//,2X,A,/)') '*** Circular '//
     &       'Intensity difference CID=(R-L)/(2*Int) (*E3) '//
     &          '(London) ***','Mode   Frequency      CID(par)'//
     &          '      CID(perp)     CID(0)       CID(180)'
        DO 800 I = 1, NUMMOD
            WRITE(LUPRI,'(1X,I4,2X,F10.2,4F14.6)') I, XTKAYS*FREQ(I),
     &              CZL(I)*DE3,CXL(I)*DE3,
     &              C0L(I)*DE3,CBL(I)*DE3
  800   CONTINUE
C
C Difference NoLondon
C
        WRITE (LUPRI,'(/,3X,A,//,2X,A,/)') '**** Difference '//
     &         'Parameter R-L (A**4amu**-1) (NoLondon) ****',
     &        'Mode   Frequency     DELTApar      DELTAperp   '//
     &         '   DELTA0       DELTA180'
        DO 760 I = 1, NUMMOD
          WRITE(LUPRI,'(1X,I4,2X,F10.2,4F14.6)')I,XTKAYS*FREQ(I),
     &    DZNL(I)*XFRAIN,DXNL(I)*XFRAIN,D0NL(I)*XFRAIN,DBNL(I)*XFRAIN
 760    CONTINUE
C
C  CIDs NoLondon
C
           WRITE (LUPRI,'(/,3X,A,//,2X,A,/)') '*** Circular '//
     &       'Intensity difference CID=(R-L)/(2*Int) (*E3) '//
     &       ' (NoLon) ***',
     &          'Mode   Frequency      CID(par)    CID(perp)'//
     &          '     CID(0)       CID(180)'
         DO 810 I = 1, NUMMOD
            WRITE(LUPRI,'(1X,I4,2X,F10.2,4F14.6)') I, XTKAYS*FREQ(I),
     &              CZNL(I)*DE3,CXNL(I)*DE3,
     &              C0NL(I)*DE3,CBNL(I)*DE3
  810    CONTINUE
       END IF
  500 CONTINUE
C
C  End of Raman and VROA output for one frequency,taking next
C
      WRITE (LUPRI,'(/,3X,A/)') '**** For definitions, '//
     &    'see the DALTON manual and references therein. ****'
C
      RETURN
      END
C  /* Deck detirp */
      SUBROUTINE DETIRP(EVEC,II,FREQ,NCORD,EVC1,EVC2,GVC1,GVC2,
     &     ITRANS,ATCHRG,WORK,LWORK,IRP)
C
C     Determines the irep of a given normal coordinate, by comparing
C     the results of symmetry operations with the character table.
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <symmet.h>
#include <codata.h>
#include <pgroup.h>
#include <optinf.h>
C
      DIMENSION EVEC(NCORD,NCORD), EVC1(NCORD), EVC2(NCORD)
      DIMENSION GVC1(NCORD), GVC2(NCORD), ATCHRG(NCORD), WORK(LWORK)
      CHARACTER*5 IRPTXT
      DIMENSION ICHRVC(0:7), ITMP(0:7), ITRANS(NCORD/3)
      LOGICAL FOUND, CHANGE
#include <ibtfun.h>
      SXFAMU = SQRT(XFAMU)
C
C     First we have to expand the atomic coordinates. The geometry
C     is placed in GVC1.
C
      IATOM = 0
      DO 10 ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
         DO 15 ISYMOP = 0, MAXOPR
            IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
               GVC1(IATOM*3+1) =
     &              PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENT)
               GVC1(IATOM*3+2) =
     &              PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENT)
               GVC1(IATOM*3+3) =
     &              PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENT)
               IATOM = IATOM + 1
               ATCHRG(IATOM) = CHARGE(ICENT)
            END IF
 15      CONTINUE
 10   CONTINUE
C
C     The normal coordinates are placed in EVC1.
C
      DO 17 I = 1, NCORD
         EVC1(I) = EVEC(I,II)*SXFAMU
 17   CONTINUE
C
C     If visualization of normal modes is requested, we do it here
C
         IF (VRVIBA) THEN
            CALL MKVRVB(NCORD,IATOM,GVC1,EVC1,ATCHRG,II,FREQ,WORK,LWORK)
         END IF
C
C     Then we let each symmetry operation work on both the atomic
C     and the normal coordinates.
C
      CALL IZERO(ICHRVC,8)
      DO 20 ISYMOP = 0, MAXREP
         DO 25 ICENT = 1, NCORD/3
            EVC2((ICENT-1)*3+1) =
     &        PT(IBTAND(ISYMAX(1,1),ISYMOP))*EVC1((ICENT-1)*3+1)
            GVC2((ICENT-1)*3+1) =
     &        PT(IBTAND(ISYMAX(1,1),ISYMOP))*GVC1((ICENT-1)*3+1)
            EVC2((ICENT-1)*3+2) =
     &        PT(IBTAND(ISYMAX(2,1),ISYMOP))*EVC1((ICENT-1)*3+2)
            GVC2((ICENT-1)*3+2) =
     &        PT(IBTAND(ISYMAX(2,1),ISYMOP))*GVC1((ICENT-1)*3+2)
            EVC2((ICENT-1)*3+3) =
     &        PT(IBTAND(ISYMAX(3,1),ISYMOP))*EVC1((ICENT-1)*3+3)
            GVC2((ICENT-1)*3+3) =
     &        PT(IBTAND(ISYMAX(3,1),ISYMOP))*GVC1((ICENT-1)*3+3)
 25      CONTINUE
C
C     If the geometry is affected by the symmetry operation,
C     we have to swap the transformed normal coordinates accordingly.
C     (Otherwise we can't compare the transformed and the
C     non-transformed normal coordinates). First we have to determine
C     how the atomic centers have been transformed.
C
         CHANGE = .FALSE.
         DO 27 I = 1, NCORD
            IF (GVC1(I)*GVC2(I) .LT. -1.0D-10) CHANGE = .TRUE.
 27      CONTINUE
         IF (CHANGE) THEN
            CALL IZERO(ITRANS,IATOM)
            ICENT = 1
 30         CONTINUE
            IF (ITRANS(ICENT) .EQ. 0) THEN
               CX = GVC1((ICENT-1)*3+1)
               CY = GVC1((ICENT-1)*3+2)
               CZ = GVC1((ICENT-1)*3+3)
               ITRANS(ICENT) = ICENT
               I = ICENT + 1
               FOUND = .FALSE.
 35            CONTINUE
               IF ((.NOT. FOUND) .AND. (I .LE. IATOM)) THEN
                  IF ((ABS(CX-GVC2((I-1)*3+1)) .LT. 1.0D-10)  .AND.
     &                 (ABS(CY-GVC2((I-1)*3+2)) .LT. 1.0D-10)  .AND.
     &                 (ABS(CZ-GVC2((I-1)*3+3)) .LT. 1.0D-10)) THEN
                     ITRANS(I) = ICENT
                     ITRANS(ICENT) = I
                     FOUND = .TRUE.
                  ELSE
                     I = I + 1
                     GOTO 35
                  END IF
               END IF
            END IF
            IF (ICENT .LT. IATOM) THEN
               ICENT = ICENT + 1
               GOTO 30
            END IF
C
C     The "swapped" normal coordinates are temporarily stored in GVC2.
C     Then they're moved back to EVC2.
C
            DO 40 I = 1, NCORD/3
               GVC2((I-1)*3+1) = EVC2((ITRANS(I)-1)*3+1)
               GVC2((I-1)*3+2) = EVC2((ITRANS(I)-1)*3+2)
               GVC2((I-1)*3+3) = EVC2((ITRANS(I)-1)*3+3)
 40         CONTINUE
            DO 45 I = 1, NCORD
               EVC2(I) = GVC2(I)
 45         CONTINUE
         END IF
C
C     We construct the scalar product between the two vectors. The sign
C     of this product determines the character.
C
         SCLPRD = DDOT(NCORD,EVC1,1,EVC2,1)
         IF (SCLPRD .GT. 0.0D0) THEN
            ICHRVC(ISYMOP) = 1
         ELSE
            ICHRVC(ISYMOP) = -1
         END IF
 20   CONTINUE
C
C     The character string is compared with the character table, to
C     determine which irep the given normal mode belongs to.
C
      DO 60 I = 0, 7
         ITMP(I) = ICHRVC(I)
         ICHRVC(I) = 0
 60   CONTINUE
      KK = 0
      IND = 0
 65   CONTINUE
      IF (KK .LE. 7) THEN
         DO 67 I = 0, MAXREP
            IF (JSOP(I) .EQ. KK) THEN
               ICHRVC(KK) = ITMP(IND)
               IND = IND + 1
            END IF
 67      CONTINUE
         KK = KK + 1
         GOTO 65
      END IF
C
C     The character string is compared with the character table, to
C     determine which irep the given normal mode belongs to.
C
      IREP = 0
      FOUND = .FALSE.
 50   CONTINUE
      IF ((.NOT. FOUND) .AND. (IREP .LE. MAXREP)) THEN
         FOUND = .TRUE.
         DO 52 I = 0, MAXREP
            IF (ICHRVC(JSOP(I)) .NE. IXVAL(JSOP(I),IREP))
     &           FOUND = .FALSE.
 52      CONTINUE
         IREP = IREP + 1
         GOTO 50
      END IF
      IREP = IREP - 1
      IF (FOUND) THEN
         IRP = IREP
      ELSE
         IRP = -1
      END IF
      RETURN
      END
C /* Deck prpvav */
      SUBROUTINE PRPVAV(EVAL,NUMMOD,DIPM0,DIPMF,SUSTO0,SUSTF,SUSAV,
     &                  GFAC0,GFACF,GFACAV,QUAD0,QUADF,QUADAV,SIGMA0,
     &                  SIGMAF,SIGMAV,SRC0,SRCF,SRCAV,POLAR0,POLARF,
     &                  POLAV,ALFA0,ALFAF,ALFAV,EFG0,EFGF,EFGAV,SSJ0,
     &                  SSJF,SSJAV,ROANL0,ROANLF,ROANLA,ROALO0,ROALOF,
     &                  ROALOA,CSTRA,SCTRA)
C
C     This subroutine calculates vibrationally averaged molecular
C     properties assuming that the calculation is performed at the effective
C     geometry. K.Ruud, San Diego, March 1999
C
#include <implicit.h>
#include <codata.h>
#include <thrzer.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
      PARAMETER (D1 = 1.0D0, D4 = 4.0D0,
     &           THRESH = 1.0D-5, D0 = 0.0D0, D2 = 2.0D0)
      PARAMETER (FAC = ALPHA2/(PMASS*XFAMU*XFAMU)*XTHZ*1.0D-3)
      PARAMETER (FACTOR = (288.0D-30)*(PI**2)*XFMOL*(XTANG**4),
     &           D3 = 3.0D0)
#include <cbilnr.h>
#include <nuclei.h>
      DIMENSION EVAL(NUMMOD), DIPM0(3), DIPMF(3,NUMMOD), SUSTO0(3,3),
     &          SUSTF(3,3,NUMMOD), SUSAV(3,3), GFAC0(3,3), GFACAV(3,3),
     &          GFACF(3,3,NUMMOD), QUAD0(3,3), QUADF(3,3,NUMMOD),
     &          QUADAV(3,3), SIGMA0(3,3,MXCENT),
     &          SIGMAF(3,3,MXCENT,NUMMOD), SIGMAV(3,3,MXCENT),
     &          SRC0(3,3,MXCENT), SRCF(3,3,MXCENT,NUMMOD),
     &          SRCAV(3,3,MXCENT), POLAR0(3,3), POLARF(3,3,NUMMOD),
     &          POLAV(3,3), ALFA0(3,3,MXFR), ALFAF(3,3,MXFR,NUMMOD),
     &          ALFAV(3,3,MXFR), EFG0(3,3,MXCENT),
     &          EFGF(3,3,MXCENT,NUMMOD),
     &          EFGAV(3,3,MXCENT), SSJ0(MXCOOR,MXCOOR),
     &          SSJF(MXCOOR,MXCOOR,NUMMOD), SSJAV(MXCOOR,MXCOOR),
     &          ROANL0(3,3,MXFR), ROANLF(3,3,MXFR,MXCOOR),
     &          ROANLA(3,3,MXFR), ROALO0(3,3,MXFR),
     &          ROALOF(3,3,MXFR,MXCOOR), ROALOA(3,3,MXFR),
     &          CSTRA(*), SCTRA(*)
C
#include <symmet.h>
#include <cbiwlk.h>
#include <abainf.h>
#include <orgcom.h>
#include <ibtfun.h>
C
      PREFAC = D1/(D4*XFAMU)
C
C     We have to check if we accidentaly get some negative frequencies at
C     the effective geometry
C
      INEG = 0
      DO IMOD = 1, NUMMOD
         IF (EVAL(IMOD) .LT. D0) INEG = INEG + 1
      END DO
      IF (INEG .GT. 0) THEN
         NUMMOD = NUMMOD - INEG
         WRITE (LUPRI,'(/A,I3,A,/A)') ' WARNING: There are ',INEG,
     &        ' negative frequencies at effective geometry',
     &        ' These are ignored in the following analysis'
      END IF
C
C     Print output header
C
      CALL TITLER('Vibrationally averaged molecular properties','*',108)
C
      IF (.NOT. DOTEMP) THEN
         NTEMP = 1
         TEMP(NTEMP) = 0.0D0
      END IF
      CALL HEADER('Vibrational Frequencies',1)
      WRITE (LUPRI,'(2X,A,/,2X,A)')
     &     ' mode   cm-1     hartrees ',
     &     '--------------------------'
      DO IMODE = 1, NUMMOD
         WRITE (LUPRI,'(2X,I2,F10.2,3X,F10.6)')
     &        IMODE,XTKAYS*SQRT(EVAL(IMODE)),
     &        SQRT(EVAL(IMODE))
      END DO
      DO ITMP = 1, NTEMP
            WRITE (LUPRI,'(////A,/1X,A,F12.6,A,/A)')
     &           '----------------------------'//
     &           '------------------------------',
     &           'Vibrational corrections for temperature T = ',
     &           TEMP(ITMP),' K','----------------------------'//
     &           '------------------------------'
C
C     Dipole moment
C
      DIPAVX = D0
      DIPAVY = D0
      DIPAVZ = D0
      DO IMOD = 1, NUMMOD
         IF (TEMP(ITMP) .LT. THRZER) THEN
            TFAC2 = D1/SQRT(EVAL(IMOD))
         ELSE
            TFAC  = SQRT(EVAL(IMOD))*AUTK/(D2*TEMP(ITMP))
            TFAC2 = D1/(DTANH(TFAC)*SQRT(EVAL(IMOD)))
         END IF
         DIPAVX = DIPAVX + DIPMF(1,IMOD)*TFAC2
         DIPAVY = DIPAVY + DIPMF(2,IMOD)*TFAC2
         DIPAVZ = DIPAVZ + DIPMF(3,IMOD)*TFAC2
      END DO
      DIPAVX = DIPAVX*PREFAC
      DIPAVY = DIPAVY*PREFAC
      DIPAVZ = DIPAVZ*PREFAC
C
C     Magnetizabilities
C
      IF (MAGSUS) THEN
         CALL DZERO(SUSAV,9)
         DO IMOD = 1, NUMMOD
            IF (TEMP(ITMP) .LT. THRZER) THEN
               TFAC2 = D1/SQRT(EVAL(IMOD))
            ELSE
               TFAC  = SQRT(EVAL(IMOD))*AUTK/(D2*TEMP(ITMP))
               TFAC2 = D1/(DTANH(TFAC)*SQRT(EVAL(IMOD)))
            END IF
            DO ICOOR = 1, 3
               DO JCOOR = 1, 3
                  SUSAV(ICOOR,JCOOR) = SUSAV(ICOOR,JCOOR)
     &                 + SUSTF(ICOOR,JCOOR,IMOD)*TFAC2
               END DO
            END DO
         END DO
         CALL DSCAL(9,PREFAC,SUSAV,1)
      END IF
C
C     Rotational g tensor
C
      IF (MOLGFA) THEN
         CALL DZERO(GFACAV,9)
         DO IMOD = 1, NUMMOD
            IF (TEMP(ITMP) .LT. THRZER) THEN
               TFAC2 = D1/SQRT(EVAL(IMOD))
            ELSE
               TFAC  = SQRT(EVAL(IMOD))*AUTK/(D2*TEMP(ITMP))
               TFAC2 = D1/(DTANH(TFAC)*SQRT(EVAL(IMOD)))
            END IF
            DO ICOOR = 1, 3
               DO JCOOR = 1, 3
                  GFACAV(ICOOR,JCOOR) = GFACAV(ICOOR,JCOOR)
     &                 + GFACF(ICOOR,JCOOR,IMOD)*TFAC2
               END DO
            END DO
         END DO
         CALL DSCAL(9,PREFAC,GFACAV,1)
      END IF
C
C     Molecular quadrupole moment
C
      IF (QUADRU) THEN
         CALL DZERO(QUADAV,9)
         DO IMOD = 1, NUMMOD
            IF (TEMP(ITMP) .LT. THRZER) THEN
               TFAC2 = D1/SQRT(EVAL(IMOD))
            ELSE
               TFAC  = SQRT(EVAL(IMOD))*AUTK/(D2*TEMP(ITMP))
               TFAC2 = D1/(DTANH(TFAC)*SQRT(EVAL(IMOD)))
            END IF
            DO ICOOR = 1, 3
               DO JCOOR = 1, 3
                  QUADAV(ICOOR,JCOOR) = QUADAV(ICOOR,JCOOR)
     &                 + QUADF(ICOOR,JCOOR,IMOD)*TFAC2
               END DO
            END DO
         END DO
         CALL DSCAL(9,PREFAC,QUADAV,1)
      END IF
C
C     Nuclear shieldings
C
      IF (SHIELD) THEN
         CALL DZERO(SIGMAV,9*MXCENT)
         DO IMOD = 1, NUMMOD
            IF (TEMP(ITMP) .LT. THRZER) THEN
               TFAC2 = D1/SQRT(EVAL(IMOD))
            ELSE
               TFAC  = SQRT(EVAL(IMOD))*AUTK/(D2*TEMP(ITMP))
               TFAC2 = D1/(DTANH(TFAC)*SQRT(EVAL(IMOD)))
            END IF
            DO IATOM = 1, NUCDEP
               DO ICOOR = 1, 3
                  DO JCOOR = 1, 3
                     SIGMAV(ICOOR,JCOOR,IATOM)=SIGMAV(ICOOR,JCOOR,IATOM)
     &                 + SIGMAF(ICOOR,JCOOR,IATOM,IMOD)*TFAC2
                  END DO
               END DO
            END DO
         END DO
         CALL DSCAL(9*MXCENT,PREFAC,SIGMAV,1)
      END IF
C
C     Nuclear spin-rotation constants
C
      IF (SPINRO) THEN
         CALL DZERO(SRCAV,9*MXCENT)
         DO IMOD = 1, NUMMOD
            IF (TEMP(ITMP) .LT. THRZER) THEN
               TFAC2 = D1/SQRT(EVAL(IMOD))
            ELSE
               TFAC  = SQRT(EVAL(IMOD))*AUTK/(D2*TEMP(ITMP))
               TFAC2 = D1/(DTANH(TFAC)*SQRT(EVAL(IMOD)))
            END IF
            DO IATOM = 1, NUCDEP
               DO ICOOR = 1, 3
                  DO JCOOR = 1, 3
                     SRCAV(ICOOR,JCOOR,IATOM) = SRCAV(ICOOR,JCOOR,IATOM)
     &                 + SRCF(ICOOR,JCOOR,IATOM,IMOD)*TFAC2
                  END DO
               END DO
            END DO
         END DO
         CALL DSCAL(9*MXCENT,PREFAC,SRCAV,1)
      END IF
C
C     Polarizabilities
C
      IF (POLAR) THEN
         CALL DZERO(POLAV,9)
         DO IMOD = 1, NUMMOD
            IF (TEMP(ITMP) .LT. THRZER) THEN
               TFAC2 = D1/SQRT(EVAL(IMOD))
            ELSE
               TFAC  = SQRT(EVAL(IMOD))*AUTK/(D2*TEMP(ITMP))
               TFAC2 = D1/(DTANH(TFAC)*SQRT(EVAL(IMOD)))
            END IF
            DO ICOOR = 1, 3
               DO JCOOR = 1, 3
                  POLAV(ICOOR,JCOOR) = POLAV(ICOOR,JCOOR)
     &                 + POLARF(ICOOR,JCOOR,IMOD)*TFAC2
               END DO
            END DO
         END DO
         CALL DSCAL(9,PREFAC,POLAV,1)
      END IF
      IF (ALFA) THEN
         CALL DZERO(ALFAV,9*MXFR)
         DO IMOD = 1, NUMMOD
            IF (TEMP(ITMP) .LT. THRZER) THEN
               TFAC2 = D1/SQRT(EVAL(IMOD))
            ELSE
               TFAC  = SQRT(EVAL(IMOD))*AUTK/(D2*TEMP(ITMP))
               TFAC2 = D1/(DTANH(TFAC)*SQRT(EVAL(IMOD)))
            END IF
            DO IFRQR = 1, NFRVAL
               DO ICOOR = 1, 3
                  DO JCOOR = 1, 3
                     ALFAV(ICOOR,JCOOR,IFRQR) = ALFAV(ICOOR,JCOOR,IFRQR)
     &                  + ALFAF(ICOOR,JCOOR,IFRQR,IMOD)*TFAC2
                  END DO
               END DO
            END DO
         END DO
         CALL DSCAL(9*MXFR,PREFAC,ALFAV,1)
      END IF
C
C     Nuclear quadrupole moments
C
      IF (NQCC) THEN
         CALL DZERO(EFGAV,9*MXCENT)
         DO IMOD = 1, NUMMOD
            IF (TEMP(ITMP) .LT. THRZER) THEN
               TFAC2 = D1/SQRT(EVAL(IMOD))
            ELSE
               TFAC  = SQRT(EVAL(IMOD))*AUTK/(D2*TEMP(ITMP))
               TFAC2 = D1/(DTANH(TFAC)*SQRT(EVAL(IMOD)))
            END IF
            DO IATOM = 1, NUCDEP
               DO ICOOR = 1, 3
               DO JCOOR = 1, 3
                  EFGAV(ICOOR,JCOOR,IATOM) = EFGAV(ICOOR,JCOOR,IATOM) +
     &                 EFGF(ICOOR,JCOOR,IATOM,IMOD)*TFAC2
               END DO
               END DO
            END DO
         END DO
         CALL DSCAL(9*MXCENT,PREFAC,EFGAV,1)
      END IF
C
C     Indirect spin-spin coupling constants
C
      IF (SPNSPN) THEN
         CALL DZERO(SSJAV,MXCOOR*MXCOOR)
         DO IMOD = 1, NUMMOD
            IF (TEMP(ITMP) .LT. THRZER) THEN
               TFAC2 = D1/SQRT(EVAL(IMOD))
            ELSE
               TFAC  = SQRT(EVAL(IMOD))*AUTK/(D2*TEMP(ITMP))
               TFAC2 = D1/(DTANH(TFAC)*SQRT(EVAL(IMOD)))
            END IF
            DO ICOOR = 1, 3*NUCDEP
               DO JCOOR = 1, 3*NUCDEP
                  SSJAV(ICOOR,JCOOR) = SSJAV(ICOOR,JCOOR) +
     &                 SSJF(ICOOR,JCOOR,IMOD)*TFAC2
               END DO
            END DO
         END DO
         CALL DSCAL(MXCOOR*MXCOOR,PREFAC,SSJAV,1)
      END IF
C
C     Optical rotation
C
      IF (OPTROT) THEN
         CALL DZERO(ROANLA,9*MXFR)
         CALL DZERO(ROALOA,9*MXFR)
         DO IMOD = 1, NUMMOD
            IF (TEMP(ITMP) .LT. THRZER) THEN
               TFAC2 = D1/SQRT(EVAL(IMOD))
            ELSE
               TFAC  = SQRT(EVAL(IMOD))*AUTK/(D2*TEMP(ITMP))
               TFAC2 = D1/(DTANH(TFAC)*SQRT(EVAL(IMOD)))
            END IF
            DO IFRQ = 1, NFRVAL
            DO ICOOR = 1, 3
               DO JCOOR = 1, 3
                  ROANLA(ICOOR,JCOOR,IFRQ) = ROANLA(ICOOR,JCOOR,IFRQ) +
     &                 ROANLF(ICOOR,JCOOR,IFRQ,IMOD)*TFAC2
                  ROALOA(ICOOR,JCOOR,IFRQ) = ROALOA(ICOOR,JCOOR,IFRQ) +
     &                 ROALOF(ICOOR,JCOOR,IFRQ,IMOD)*TFAC2
               END DO
            END DO
            END DO
         END DO
         CALL DSCAL(9*MXFR,PREFAC,ROANLA,1)
         CALL DSCAL(9*MXFR,PREFAC,ROALOA,1)
      END IF
C
C     Dipole moments
C
      CALL HEADER('Dipole moment (au)',-1)
      WRITE (LUPRI,'(4X,A,7X,A,6X,A)') 'Effective geometry',
     &     '<P^(0)_2>_eff','Vibrationally corrected'
      WRITE (LUPRI,'(A)') ' -------------------------------'//
     &     '----------------------------------------'
      WRITE (LUPRI,'(3(/2X,A1,3(2X,F16.8,4X)))')
     &     'X',DIPM0(1),DIPAVX,DIPM0(1)+DIPAVX,
     &     'Y',DIPM0(2),DIPAVY,DIPM0(2)+DIPAVY,
     &     'Z',DIPM0(3),DIPAVZ,DIPM0(3)+DIPAVZ
C
C     Magnetizabilities
C
      IF (MAGSUS) THEN
         CALL HEADER ('Magnetizabilities (au)',-1)
         WRITE (LUPRI,'(4X,A,7X,A,6X,A)') 'Effective geometry',
     &        '<P^(0)_2>_eff','Vibrationally corrected'
         WRITE (LUPRI,'(A)') ' -------------------------------'//
     &        '----------------------------------------'
         WRITE (LUPRI,'(6(/1X,A2,3(2X,F16.8,4X)))')
     &        'XX',SUSTO0(1,1),SUSAV(1,1),SUSTO0(1,1)+SUSAV(1,1),
     &        'XY',SUSTO0(1,2),SUSAV(1,2),SUSTO0(1,2)+SUSAV(1,2),
     &        'XZ',SUSTO0(1,3),SUSAV(1,3),SUSTO0(1,3)+SUSAV(1,3),
     &        'YY',SUSTO0(2,2),SUSAV(2,2),SUSTO0(2,2)+SUSAV(2,2),
     &        'YZ',SUSTO0(2,3),SUSAV(2,3),SUSTO0(2,3)+SUSAV(2,3),
     &        'ZZ',SUSTO0(3,3),SUSAV(3,3),SUSTO0(3,3)+SUSAV(3,3)
      END IF
C
C     Rotational g tensor
C
      IF (MOLGFA) THEN
         CALL HEADER ('Rotational g tensor',-1)
         WRITE (LUPRI,'(4X,A,7X,A,6X,A)') 'Effective geometry',
     &        '<P^(0)_2>_eff','Vibrationally corrected'
         WRITE (LUPRI,'(A)') ' -------------------------------'//
     &        '----------------------------------------'
         WRITE (LUPRI,'(9(/1X,A2,3(2X,F16.8,4X)))')
     &        'AX',GFAC0(1,1),GFACAV(1,1),GFAC0(1,1)+GFACAV(1,1),
     &        'BX',GFAC0(1,2),GFACAV(1,2),GFAC0(1,2)+GFACAV(1,2),
     &        'CX',GFAC0(1,3),GFACAV(1,3),GFAC0(1,3)+GFACAV(1,3),
     &        'AY',GFAC0(2,1),GFACAV(2,1),GFAC0(2,1)+GFACAV(2,1),
     &        'BY',GFAC0(2,2),GFACAV(2,2),GFAC0(2,2)+GFACAV(2,2),
     &        'CY',GFAC0(2,3),GFACAV(2,3),GFAC0(2,3)+GFACAV(2,3),
     &        'AZ',GFAC0(3,1),GFACAV(3,1),GFAC0(3,1)+GFACAV(3,1),
     &        'BZ',GFAC0(3,2),GFACAV(3,2),GFAC0(3,2)+GFACAV(3,2),
     &        'CZ',GFAC0(3,3),GFACAV(3,3),GFAC0(3,3)+GFACAV(3,3)
      END IF
C
C     Molecular quadrupole moment
C
      IF (QUADRU) THEN
         CALL HEADER ('Molecular quadrupole moment (au)',-1)
         WRITE (LUPRI,'(4X,A,7X,A,6X,A)') 'Effective geometry',
     &        '<P^(0)_2>_eff','Vibrationally corrected'
         WRITE (LUPRI,'(A)') ' -------------------------------'//
     &        '----------------------------------------'
         WRITE (LUPRI,'(6(/1X,A2,3(2X,F16.8,4X)))')
     &        'XX',QUAD0(1,1),QUADAV(1,1),QUAD0(1,1)+QUADAV(1,1),
     &        'XY',QUAD0(1,2),QUADAV(1,2),QUAD0(1,2)+QUADAV(1,2),
     &        'XZ',QUAD0(1,3),QUADAV(1,3),QUAD0(1,3)+QUADAV(1,3),
     &        'YY',QUAD0(2,2),QUADAV(2,2),QUAD0(2,2)+QUADAV(2,2),
     &        'YZ',QUAD0(2,3),QUADAV(2,3),QUAD0(2,3)+QUADAV(2,3),
     &        'ZZ',QUAD0(3,3),QUADAV(3,3),QUAD0(3,3)+QUADAV(3,3)
      END IF
C
C     Nuclear shieldings
C
      IF (SHIELD) THEN
         CALL HEADER ('Nuclear magnetic shielding constants (ppm)',-1)
         DO IATOM = 1, NUCDEP
            WRITE (LUPRI,'(/A,A)') 'Atom : ',NAMDEP(IATOM)
            WRITE (LUPRI,'(4X,A,7X,A,6X,A)') 'Effective geometry',
     &           '<P^(0)_2>_eff','Vibrationally corrected'
            WRITE (LUPRI,'(A)') ' -------------------------------'//
     &           '----------------------------------------'
            WRITE (LUPRI,'(9(/1X,A4,3(2X,F16.8,4X)))')
     &        'Bxmx',SIGMA0(1,1,IATOM),SIGMAV(1,1,IATOM),
     &           SIGMA0(1,1,IATOM)+SIGMAV(1,1,IATOM),
     &        'Bymx',SIGMA0(2,1,IATOM),SIGMAV(2,1,IATOM),
     &           SIGMA0(2,1,IATOM)+SIGMAV(2,1,IATOM),
     &        'Bzmx',SIGMA0(3,1,IATOM),SIGMAV(3,1,IATOM),
     &           SIGMA0(3,1,IATOM)+SIGMAV(3,1,IATOM),
     &        'Bxmy',SIGMA0(1,2,IATOM),SIGMAV(1,2,IATOM),
     &           SIGMA0(1,2,IATOM)+SIGMAV(1,2,IATOM),
     &        'Bymy',SIGMA0(2,2,IATOM),SIGMAV(2,2,IATOM),
     &           SIGMA0(2,2,IATOM)+SIGMAV(2,2,IATOM),
     &        'Bzmy',SIGMA0(3,2,IATOM),SIGMAV(3,2,IATOM),
     &           SIGMA0(3,2,IATOM)+SIGMAV(3,2,IATOM),
     &        'Bxmz',SIGMA0(1,3,IATOM),SIGMAV(1,3,IATOM),
     &           SIGMA0(1,3,IATOM)+SIGMAV(1,3,IATOM),
     &        'Bymz',SIGMA0(2,3,IATOM),SIGMAV(2,3,IATOM),
     &           SIGMA0(2,3,IATOM)+SIGMAV(2,3,IATOM),
     &        'Bzmz',SIGMA0(3,3,IATOM),SIGMAV(3,3,IATOM),
     &           SIGMA0(3,3,IATOM)+SIGMAV(3,3,IATOM)
         END DO
      END IF
C
C     Nuclear spin-rotation constants
C
      IF (SPINRO) THEN
         CALL HEADER ('Nuclear spin-rotation tensors',-1)
         DO IATOM = 1, NUCDEP
            WRITE (LUPRI,'(/A,A)') 'Atom : ',NAMDEP(IATOM)
            WRITE (LUPRI,'(4X,A,7X,A,6X,A)') 'Effective geometry',
     &           '<P^(0)_2>_eff','Vibrationally corrected'
            WRITE (LUPRI,'(A)') ' -------------------------------'//
     &           '----------------------------------------'
            WRITE (LUPRI,'(9(/1X,A4,3(2X,F16.8,4X)))')
     &        'Jxmx',SRC0(1,1,IATOM),SRCAV(1,1,IATOM),
     &           SRC0(1,1,IATOM)+SRCAV(1,1,IATOM),
     &        'Jymx',SRC0(2,1,IATOM),SRCAV(2,1,IATOM),
     &           SRC0(2,1,IATOM)+SRCAV(2,1,IATOM),
     &        'Jzmx',SRC0(3,1,IATOM),SRCAV(3,1,IATOM),
     &           SRC0(3,1,IATOM)+SRCAV(3,1,IATOM),
     &        'Jxmy',SRC0(1,2,IATOM),SRCAV(1,2,IATOM),
     &           SRC0(1,2,IATOM)+SRCAV(1,2,IATOM),
     &        'Jymy',SRC0(2,2,IATOM),SRCAV(2,2,IATOM),
     &           SRC0(2,2,IATOM)+SRCAV(2,2,IATOM),
     &        'Jzmy',SRC0(3,2,IATOM),SRCAV(3,2,IATOM),
     &           SRC0(3,2,IATOM)+SRCAV(3,2,IATOM),
     &        'Jxmz',SRC0(1,3,IATOM),SRCAV(1,3,IATOM),
     &           SRC0(1,3,IATOM)+SRCAV(1,3,IATOM),
     &        'Jymz',SRC0(2,3,IATOM),SRCAV(2,3,IATOM),
     &           SRC0(2,3,IATOM)+SRCAV(2,3,IATOM),
     &        'Jzmz',SRC0(3,3,IATOM),SRCAV(3,3,IATOM),
     &           SRC0(3,3,IATOM)+SRCAV(3,3,IATOM)
         END DO
      END IF
C
C     Polarizabilities
C
      IF (POLAR) THEN
         CALL HEADER ('Static polarizabilities (au)',-1)
         WRITE (LUPRI,'(/4X,A,7X,A,6X,A)') 'Effective geometry',
     &        '<P^(0)_2>_eff','Vibrationally corrected'
         WRITE (LUPRI,'(A)') ' -------------------------------'//
     &        '----------------------------------------'
         WRITE (LUPRI,'(6(/1X,A2,3(2X,F16.8,4X)))')
     &        'XX',POLAR0(1,1),POLAV(1,1),POLAR0(1,1)+POLAV(1,1),
     &        'XY',POLAR0(1,2),POLAV(1,2),POLAR0(1,2)+POLAV(1,2),
     &        'XZ',POLAR0(1,3),POLAV(1,3),POLAR0(1,3)+POLAV(1,3),
     &        'YY',POLAR0(2,2),POLAV(2,2),POLAR0(2,2)+POLAV(2,2),
     &        'YZ',POLAR0(2,3),POLAV(2,3),POLAR0(2,3)+POLAV(2,3),
     &        'ZZ',POLAR0(3,3),POLAV(3,3),POLAR0(3,3)+POLAV(3,3)
      END IF
      IF (ALFA) THEN
         DO I = 1, NFRVAL
            WRITE(LUPRI,'(//10X,A,F10.6)')
     &           'Polarizability (au) for frequency', FRVAL(I)
            WRITE (LUPRI,'(/4X,A,7X,A,6X,A)') 'Effective geometry',
     &           '<P^(0)_2>_eff','Vibrationally corrected'
            WRITE (LUPRI,'(A)') ' -------------------------------'//
     &           '----------------------------------------'
            WRITE (LUPRI,'(6(/1X,A2,3(2X,F16.8,4X)))')
     &        'XX',ALFA0(1,1,I),ALFAV(1,1,I),ALFA0(1,1,I)+ALFAV(1,1,I),
     &        'XY',ALFA0(1,2,I),ALFAV(1,2,I),ALFA0(1,2,I)+ALFAV(1,2,I),
     &        'XZ',ALFA0(1,3,I),ALFAV(1,3,I),ALFA0(1,3,I)+ALFAV(1,3,I),
     &        'YY',ALFA0(2,2,I),ALFAV(2,2,I),ALFA0(2,2,I)+ALFAV(2,2,I),
     &        'YZ',ALFA0(2,3,I),ALFAV(2,3,I),ALFA0(2,3,I)+ALFAV(2,3,I),
     &        'ZZ',ALFA0(3,3,I),ALFAV(3,3,I),ALFA0(3,3,I)+ALFAV(3,3,I)
         END DO
      END IF
C
C     Nuclear quadrupole moments
C
C
      IF (NQCC) THEN
         CALL HEADER('Nuclear quadrupole moments (au)',-1)
         DO IATOM = 1, NUCDEP
            WRITE (LUPRI,'(//A,A)') 'Atom : ',NAMDEP(IATOM)
            WRITE (LUPRI,'(/4X,A,7X,A,6X,A)') 'Effective geometry',
     &           '<P^(0)_2>_eff','Vibrationally corrected'
            WRITE (LUPRI,'(A)') ' -------------------------------'//
     &           '----------------------------------------'
            WRITE (LUPRI,'(6(/1X,A2,2X,3(2X,F16.8,4X)))')
     &        'XX',EFG0(1,1,IATOM),EFGAV(1,1,IATOM),
     &           EFG0(1,1,IATOM)+EFGAV(1,1,IATOM),
     &        'XY',EFG0(1,2,IATOM),EFGAV(1,2,IATOM),
     &           EFG0(1,2,IATOM)+EFGAV(1,2,IATOM),
     &        'XZ',EFG0(1,3,IATOM),EFGAV(1,3,IATOM),
     &           EFG0(1,3,IATOM)+EFGAV(1,3,IATOM),
     &        'YY',EFG0(2,2,IATOM),EFGAV(2,2,IATOM),
     &           EFG0(2,2,IATOM)+EFGAV(2,2,IATOM),
     &        'YZ',EFG0(2,3,IATOM),EFGAV(2,3,IATOM),
     &           EFG0(2,3,IATOM)+EFGAV(2,3,IATOM),
     &        'ZZ',EFG0(3,3,IATOM),EFGAV(3,3,IATOM),
     &           EFG0(3,3,IATOM)+EFGAV(3,3,IATOM)
         END DO
      END IF
C
C     Indirect spin-spin coupling constants
C
      IF (SPNSPN) THEN
         CALL HEADER('Reduced spin-spin coupling constants at '//
     &               'effective geometry (Hz)',-1)
         CALL PRIHES(SSJ0,'SPNSPN',CSTRA,SCTRA)
         CALL HEADER('Reduced spin-spin coupling constants '//
     &        '<P^(0)_2>_eff (Hz)',-1)
         CALL PRIHES(SSJAV,'SPNSPN',CSTRA,SCTRA)
         CALL HEADER('Vibrationally averaged reduced spin-spin '//
     &               'coupling constants (Hz)',-1)
         CALL DAXPY(MXCOOR*MXCOOR,D1,SSJ0,1,SSJAV,1)
         CALL PRIHES(SSJAV,'SPNSPN',CSTRA,SCTRA)
      END IF
C
C     Optical rotation
C
      IF (OPTROT) THEN
         DO I = 1, NFRVAL
            WRITE(LUPRI,'(//10X,A,F10.6)')
     &           'London G matrix (au) for frequency', FRVAL(I)
            WRITE (LUPRI,'(/4X,A,7X,A,6X,A)') 'Effective geometry',
     &           '<P^(0)_2>_eff','Vibrationally corrected'
            WRITE (LUPRI,'(A)') ' -------------------------------'//
     &           '----------------------------------------'
            WRITE (LUPRI,'(6(/1X,A2,3(2X,F16.8,4X)))')
     &     'XX',ROALO0(1,1,I),ROALOA(1,1,I),ROALO0(1,1,I)+ROALOA(1,1,I),
     &     'XY',ROALO0(1,2,I),ROALOA(1,2,I),ROALO0(1,2,I)+ROALOA(1,2,I),
     &     'XZ',ROALO0(1,3,I),ROALOA(1,3,I),ROALO0(1,3,I)+ROALOA(1,3,I),
     &     'YY',ROALO0(2,2,I),ROALOA(2,2,I),ROALO0(2,2,I)+ROALOA(2,2,I),
     &     'YZ',ROALO0(2,3,I),ROALOA(2,3,I),ROALO0(2,3,I)+ROALOA(2,3,I),
     &     'ZZ',ROALO0(3,3,I),ROALOA(3,3,I),ROALO0(3,3,I)+ROALOA(3,3,I)
         END DO
         DO I = 1, NFRVAL
            WRITE(LUPRI,'(//10X,A,F10.6)')
     &           'No-London G matrix (au) for frequency', FRVAL(I)
            WRITE (LUPRI,'(/4X,A,7X,A,6X,A)') 'Effective geometry',
     &           '<P^(0)_2>_eff','Vibrationally corrected'
            WRITE (LUPRI,'(A)') ' -------------------------------'//
     &           '----------------------------------------'
            WRITE (LUPRI,'(6(/1X,A2,3(2X,F16.8,4X)))')
     &     'XX',ROANL0(1,1,I),ROANLA(1,1,I),ROANL0(1,1,I)+ROANLA(1,1,I),
     &     'XY',ROANL0(1,2,I),ROANLA(1,2,I),ROANL0(1,2,I)+ROANLA(1,2,I),
     &     'XZ',ROANL0(1,3,I),ROANLA(1,3,I),ROANL0(1,3,I)+ROANLA(1,3,I),
     &     'YY',ROANL0(2,2,I),ROANLA(2,2,I),ROANL0(2,2,I)+ROANLA(2,2,I),
     &     'YZ',ROANL0(2,3,I),ROANLA(2,3,I),ROANL0(2,3,I)+ROANLA(2,3,I),
     &     'ZZ',ROANL0(3,3,I),ROANLA(3,3,I),ROANL0(3,3,I)+ROANLA(3,3,I)
         END DO
         TMASS = 0.0D0
         JATOM = 0
         DO IATOM = 1, NUCIND
            DO ISYMOP = 0, MAXOPR
               IF (IBTAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
                  JATOM = JATOM + 1
                  NATTYP = IZATOM(IATOM)
                  IF (NATTYP .NE. 0 .AND. .NOT. NOORBT(IATOM)) THEN
                     AMASS = DISOTP(NATTYP,ISOTOP(JATOM),'MASS')
                     TMASS = TMASS + AMASS
                  END IF
               END IF
            END DO
         END DO
         FACTOT = FACTOR*XTKAYS*XTKAYS
         DO IFRVAL = 1, NFRVAL
            BETAL0 = -(ROALO0(1,1,IFRVAL) + ROALO0(2,2,IFRVAL) +
     &                 ROALO0(3,3,IFRVAL))/(D3*FRVAL(IFRVAL))
            BETNL0 = -(ROANL0(1,1,IFRVAL) + ROANL0(2,2,IFRVAL) +
     &                 ROANL0(3,3,IFRVAL))/(D3*FRVAL(IFRVAL))
            ALPHL0 = FACTOT*BETAL0*FRVAL(IFRVAL)*FRVAL(IFRVAL)/TMASS
            ALPHN0 = FACTOT*BETNL0*FRVAL(IFRVAL)*FRVAL(IFRVAL)/TMASS
            BETALA = -(ROALOA(1,1,IFRVAL) + ROALOA(2,2,IFRVAL) +
     &                 ROALOA(3,3,IFRVAL))/(D3*FRVAL(IFRVAL))
            BETNLA = -(ROANLA(1,1,IFRVAL) + ROANLA(2,2,IFRVAL) +
     &                 ROANLA(3,3,IFRVAL))/(D3*FRVAL(IFRVAL))
            ALPHLA = FACTOT*BETALA*FRVAL(IFRVAL)*FRVAL(IFRVAL)/TMASS
            ALPHNA = FACTOT*BETNLA*FRVAL(IFRVAL)*FRVAL(IFRVAL)/TMASS
            WRITE (LUPRI,'(//1X,A,F12.6,A,5X,F12.6,A)')
     &           'Frequency: ',FRVAL(IFRVAL),' au',
     &           XTNM/FRVAL(IFRVAL), ' nm'
            WRITE (LUPRI,'(/1X,A,3F12.6)') 'Beta(London)    ['//
     &           'eff.geom/average/total] :',
     &           BETAL0, BETALA, BETAL0+BETALA
            WRITE (LUPRI,'(1X,A,3F12.6)') 'Beta(No-London) ['//
     &           'eff.geom/average/total] :',
     &           BETNL0, BETNLA, BETNL0+BETNLA
            WRITE (LUPRI,'(/1X,A,3F12.6)') 'Optical rotation '//
     &           '(London)    [eff.geom/average/total] :',
     &           ALPHL0, ALPHLA, ALPHL0+ALPHLA
            WRITE (LUPRI,'(1X,A,3F12.6)') 'Optical rotation '//
     &           '(No-London) [eff.geom/average/total] :',
     &           ALPHN0, ALPHNA, ALPHN0+ALPHNA
         END DO
      END IF
      END DO
      RETURN
      END
C /* Deck centrif */
      SUBROUTINE CENTRIF(NATM,NINCOR,COOR,AMASS,CCOOR,CORN,CCORN,E,
     &                   QCENT,WORK,LWORK,TEMP,IPRINT)
C ***************************************************************************
C
C CALLED FROM THE SUBROUTINE EXOUT
C
C THIS SUBROUTINE COMPUTES THE CENTRIFUGAL DISTORTION OF THE MOLECULE
C Juhani Lounila 21.1.1998
C
C ***************************************************************************
#include <implicit.h>
#include <codata.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0)
      LOGICAL LINEAR
      DIMENSION COOR(3,NATM), AMASS(NATM), CC(3), AINERT(6), U(3,3)
      DIMENSION CCOOR(NATM,3), CCORN(NINCOR,3*NATM),
     &          CORN(3*NATM,3*NATM),QCENT(NINCOR)
      DIMENSION E(NINCOR), WORK(LWORK)
C
C LOCATE THE COORDINATES CC(1), CC(2) AND CC(3)
C OF THE CENTER OF MASS OF THE MOLECULE
C
      CMASS=D0
      CALL DZERO(CC,3)
      DO 10 K=1,NATM
         DO 10 I=1,3
            CC(I)=CC(I)+AMASS(K)*COOR(I,K)
 10   CONTINUE
      DO K=1,NATM
         CMASS=CMASS+AMASS(K)
      END DO
      DO I=1,3
         CC(I)=CC(I)/CMASS
      END DO
C
C FORM THE INERTIA TENSOR AINERT(I) (I=1,...,6)
C IN A COORDINATE SYSTEM WHOSE ORIGIN IS AT THE CENTER OF MASS
C (ELEMENTS: 1=XX, 2=XY, 3=YY, 4=XZ, 5=YZ, 6=ZZ)
C
      CALL DZERO(AINERT,6)
      DO K=1,NATM
         AM=AMASS(K)
         X=COOR(1,K)-CC(1)
         Y=COOR(2,K)-CC(2)
         Z=COOR(3,K)-CC(3)
         AINERT(1)=AINERT(1)+AM*(Y*Y+Z*Z)
         AINERT(2)=AINERT(2)-AM*X*Y
         AINERT(3)=AINERT(3)+AM*(X*X+Z*Z)
         AINERT(4)=AINERT(4)-AM*X*Z
         AINERT(5)=AINERT(5)-AM*Y*Z
         AINERT(6)=AINERT(6)+AM*(X*X+Y*Y)
      END DO
      IF (IPRINT .GT. 6) THEN
         CALL HEADER('Inertia Tensor in CENTRIF',-1)
         CALL OUTPAK(AINERT,3,1,LUPRI)
      END IF
C
C DIAGONALISE THE INERTIA TENSOR
C RESULTANT MATRIX OF EIGENVECTORS: TRANS(I) (I=1,...,9)
C (STORED COLUMNWISE, IN THE SAME SEQUENCE AS EIGENVALUES)
C ELEMENTS: 1=XX, 2=YX, 3=ZX, 4=XY, 5=YY, 6=ZY, 7=XZ, 8=YZ, 9=ZZ
C EIGENVECTORS: U(1,J)=(1,2,3), U(2,J)=(4,5,6), U(3,J)=(7,8,9)
C
      CALL DUNIT(U,3)
      CALL JACO(AINERT,U,3,3,3,WORK,WORK(10))
C     WRITE (*,*) 'NEW INERTIA TENSOR (XX, XY, YY, XZ, YZ, ZZ)'
C     WRITE (*,*) AINERT(1), AINERT(2), AINERT(3), AINERT(4),
C    *AINERT(5), AINERT(6)
C PRINCIPAL MOMENTS OF INERTIA AIXX, AIYY AND AIZZ
C IN A COORDINATE SYSTEM WHOSE ORIGIN IS AT THE CENTER OF MASS
C
      AINERT(2) = AINERT(3)
      AINERT(3) = AINERT(6)
      CALL ORDER(U,AINERT,3,3)
      AIXX=AINERT(1)
      AIYY=AINERT(2)
      AIZZ=AINERT(3)
C
C EIGENVECTORS U(1,J), U(2,J) AND U(3,J) OF THE TENSOR OF INERTIA
C
      CALL DGETRN(U,3,3)
C
C TRANSFORMATION TO THE NEW COORDINATE SYSTEM
C WHOSE ORIGIN IS AT THE CENTER OF MASS OF THE MOLECULE
C AND WHOSE AXES COINCIDE WITH THE PRINCIPAL AXES OF INERTIA:
C NEW CARTESIAN COORDINATES CCOOR(K,1), CCOOR(K,2) AND CCOOR(K,3)
C Note that the ordering of atoms and coordinates differ in COOR and CCOOR
C
      CALL DZERO(CCOOR,3*NATM)
      DO 60 K=1,NATM
         DO 60 I=1,3
            DO 60 J=1,3
               CCOOR(K,I)=CCOOR(K,I)+U(I,J)*(COOR(J,K)-CC(J))
 60   CONTINUE
C
C NEW 'NORMAL COORDINATES' CCORN(I) (I=1,...,ND**2)
C ND = THE NUMBER OF NORMAL COORDINATES
C ND = 3*NATM, NATM = THE NUMBER OF ATOMS
C ELEMENTS: 1=B(1,1,X), 2=B(1,1,Y), 3=B(1,1,Z), 4=B(1,2,X),...
C B(K,M,I) = THE MATRIX ELEMENT FOR THE TRANSFORMATION BETWEEN
C THE NORMAL COORDINATE K AND THE DISPLACEMENT OF THE CARTESIAN
C COORDINATE I OF THE NUCLEUS M FROM ITS EQUILIBRIUM POSITION
C (UNIT: U**(-1/2), U = THE ATOMIC MASS UNIT)
C
      ND=3*NATM
      CALL DZERO(CCORN,3*NATM*NINCOR)
      DO I=1,NINCOR
         M=ND*(I-1)
         DO K=1,NATM
            DO L=1,3
               DO J=1,3
                  NJ = J + 3*(K-1)
                  NL = L + 3*(K-1)
                  CCORN(I,NL)=CCORN(I,NL)+U(L,J)*CORN(NJ,I)
               END DO
            END DO
         END DO
      END DO
C
C CENTRIFUGAL CONTRIBUTIONS TO THE AVERAGE VALUES
C OF THE VIBRATIONAL NORMAL COORDINATES, QCENT(I)
C (I=1,...,ND-6 OR ND-5)
C UNIT: U**(1/2) ANGSTROM
C
      PREFAC = TEMP/(D2*AUTK*SQRT(XFAMU))
      DO 90 I=1, NINCOR
         AXX=D0
         AYY=D0
         AZZ=D0
         DO 100 K=1,NATM
            AM=AMASS(K)
            X=CCOOR(K,1)
            Y=CCOOR(K,2)
            Z=CCOOR(K,3)
            NX=1+3*(K-1)
            NY=2+3*(K-1)
            NZ=3+3*(K-1)
            BX=CCORN(I,NX)
            BY=CCORN(I,NY)
            BZ=CCORN(I,NZ)
            AXX=AXX+D2*AM*(Y*BY+Z*BZ)
            AYY=AYY+D2*AM*(X*BX+Z*BZ)
            AZZ=AZZ+D2*AM*(X*BX+Y*BY)
 100     CONTINUE
         IF (AIXX .LT. 1.0D-4) THEN
            FF=AZZ/AIZZ+AYY/AIYY
         ELSE
            FF=AXX/AIXX+AYY/AIYY+AZZ/AIZZ
         ENDIF
         QCENT(I)=PREFAC*FF/E(I)**2
 90   CONTINUE
      RETURN
      END
C /*  Deck reanmhes */
      SUBROUTINE REANMHES(COOR0,GRAD0,HESS0,GRADF,HESS,NCART,NCRIND)
C
C     Read numerical Hessian from LUWLK
C
#include "implicit.h"
#include "dummy.h"
      DIMENSION COOR0(NCRIND), GRAD0(NCART), HESS0(NCART,NCART),
     &          GRADF(NCART), HESS(NCART,NCART)
#include "inftap.h"
C
      CALL GPOPEN(LUWLK,ABAWLK,'OLD','SEQUENTIAL','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUWLK
      READ (LUWLK)
      READ (LUWLK)
      READ (LUWLK) COOR0, GRAD0, HESS0, GRADF, HESS
      CALL GPCLOSE(LUWLK,'KEEP')
C
      RETURN
      END
C  /* Deck machb */
      SUBROUTINE MACHB (NCORD,NINTER,BMAT,XA,YA,ZA,QQ,CTYPE,WRI,QONLY)
C
C     Subroutine to set up B matrix for transformation from Cartesian to
C     internal coordinates. This routine has been adapted from Pulay's
C     program TEXAS.
C
C
C     NCORD - 3*number of atoms input parameter
C     NINTER - number of internal coordinates
C     BMAT - the transpose of the B matrix, B(NCORD,NINTER)
C     XA,YA,ZA - cartesian coordinates in angstroms
C     QQ contains the values of the internal coordinates (OUTPUT)
C     CTYPE is 0.0 except -1.0 for stretch and 1.0 for invers stretch
C     WRI is .TRUE. if the definition of internal coordinates is
C         to be O printed
C     QONLY is .TRUE. if only coordinates (no B matrix) are to be
C           calculated
C
C
C                      INPUT DATA
C
C     Each elementary valence coordinate on a separate card
C
C                        COL. 1
C
C     'K' or ' ' (blank). If 'K' a new coordinate begins, if blank
C     then the composite internal coordinate begun earlier is
C     continued. Any other character terminates the input.
C
C                        COLS. 2-9
C
C     Scale factor for the total coordinate (only if there is 'K' in
C     column 1). Blank or zero is interpreted as 1.0.
C
C                        COLS. 21-24
C
C     Coordinate type STRE, INVR, BEND, OUT, TORS, LIN1, LIN2
C
C                        COLS. 31-40,41-50,51-60,61-70
C
C     Participating atoms A,B,C,D (FORMAT 4F10.X).
C
C     A and B are given for 'STRE' - order arbitrary
C
C     A and B are given for 'INVR' - order arbitrary
C
C     A,B,C for 'BEND' - A and B are end atoms, C is the apex atom.
C     Atom A out of the BCD plane - C is the central atom -
C     coordinate positive if A is displaced toward the vector
C     product DB*DC
C
C     Torsion A-B-C-D, positive as in the Wilson-Decius-Cross book
C     Note that the value of the coordinate varies between -pi/2 to
C     3pi/2   NOT  between -pi/2 to +PI/2.
C
C     LIN1 L  collinear bending A-B-C distorted in the plane of ABD
C     Positive if A moves toward D.
C
C     LIN2 linear bending. A-C-B distorted perpendicular to the plane.
C     ABD - positive if A moves toward the vector cross product CD*CA.
C     The linear bendings are A-C-B, i.e., the apex atom is third.
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <codata.h>
      LOGICAL WRI, QONLY, FREAD, NEW, OLD, GETB
      PARAMETER (D1 = 1.0D00, FREAD = .FALSE., D0 = 0.0D00)
C
      CHARACTER TIPUS(7)*4, TYP*4, WE*1, TLAST*4
      CHARACTER*6 NAME1, NAME2, NAME3, NAME4
#include <cbivib.h>
#include <nuclei.h>
      DIMENSION QQ(*),
     &          IA(4), U(3), V(3), W(3), Z(3), X(3), UU(3), VV(3),
     &          WW(3), ZZ(3), UV(12),
     &          BMAT(NCORD,NINTER), XA(*), YA(*), ZA(*), CTYPE(NINTER)
      EQUIVALENCE (KA,IA(1)), (KB,IA(2)), (KC,IA(3)), (KD,IA(4)),
     &            (UV(1),UU(1)), (UV(4),VV(1)), (UV(7),WW(1)),
     &            (UV(10),ZZ(1))
      DATA TIPUS/'STRE','INVR','BEND','OUT ','TORS','LIN1','LIN2'/,
     &     ANULL/1.0/
C
      IF (WRI) THEN
         CALL HEADER('DEFINITION OF INTERNAL COORDINATES',-1)
         WRITE (LUPRI,'(5X,A,25X,A,//)')
     *    '#     type     atoms','coefficient  scaling'
      END IF
      I     = 0
      CNORM = D0
      CSCAL = D0
      TLAST = '    '
      ICARD = 1
      GETB = .NOT. QONLY
   20 CONTINUE
      IF (FREAD) THEN
         READ (LUCMD,'(A1)') WE
         BACKSPACE 5
      ELSE
         WE = KWORD(ICARD)
      END IF
      NEW = WE .EQ. 'K'
      OLD = WE .EQ. ' '
C
C     All coordinates have been read - normalize last coordinate
C
      IF (.NOT.NEW .AND. .NOT.OLD) THEN
         CNORM = CSCAL/SQRT(CNORM)
         QQ(I) = QQ(I)*CNORM
         IF (GETB) THEN
            DO 50 K = 1, NCORD
               BMAT(K,I) = BMAT(K,I)*CNORM
   50       CONTINUE
         END IF
         GO TO 370
      END IF
C
C     Not finished yet...
C
      IF (FREAD) THEN
         READ (LUCMD,70)  WE, CC, CCOEF, TYP, IA
         KWORD(ICARD)   = WE
         SCALE(ICARD)   = CC
         COEF(ICARD)    = CCOEF
         ITYPCD(ICARD)   = TYP
         IATOMS(1,ICARD) = IA(1)
         IATOMS(2,ICARD) = IA(2)
         IATOMS(3,ICARD) = IA(3)
         IATOMS(4,ICARD) = IA(4)
      ELSE
         WE    = KWORD(ICARD)
         CC    = SCALE(ICARD)
         CCOEF = COEF(ICARD)
         TYP   = ITYPCD(ICARD)
         IA(1) = IATOMS(1,ICARD)
         IA(2) = IATOMS(2,ICARD)
         IA(3) = IATOMS(3,ICARD)
         IA(4) = IATOMS(4,ICARD)
      END IF
      NAME1 = '      '
      NAME2 = '      '
      NAME3 = '      '
      NAME4 = '      '
      IF (IA(1) .NE. 0) NAME1 = NAMDEP(IA(1))
      IF (IA(2) .NE. 0) NAME2 = NAMDEP(IA(2))
      IF (IA(3) .NE. 0) NAME3 = NAMDEP(IA(3))
      IF (IA(4) .NE. 0) NAME4 = NAMDEP(IA(4))
      ICARD = ICARD + 1
   70 FORMAT (A1,F9.5,F10.4,A4,6X,4F10.4)
      IF (TYP .EQ. '    ') TYP = TLAST
      TLAST = TYP
      IF (CC    .EQ. D0) CC = D1
      IF (CCOEF .EQ. D0) CCOEF  = D1
C
C     Continuation of coordinate
C
      IF (OLD) THEN
         CNORM = CNORM + CCOEF*CCOEF
C
C     New coordinate
C
      ELSE IF (NEW) THEN
C
C        Normalize old coordinate
C
         IF (I.NE.0) THEN
            IF (WRI) WRITE (LUPRI,'()')
            CNORM = CSCAL/SQRT(CNORM)
            IF (GETB) THEN
               DO 90 K = 1, NCORD
                  BMAT(K,I) = BMAT(K,I)*CNORM
   90          CONTINUE
            END IF
            QQ(I) = QQ(I)*CNORM
         END IF
C
C        Initialize new coordinate
C
         I = I + 1
         QQ(I) = D0
         CSCAL = CC
         CNORM = CCOEF*CCOEF
         IF (GETB) THEN
            DO 120 J = 1, NCORD
               BMAT(J,I) = D0
  120       CONTINUE
         END IF
      END IF
C
C     Determine type
C
      DO 150 K = 1, 7
         IF (TYP .EQ. TIPUS(K)) GO TO 170
  150 CONTINUE
C
C     If type unknown return
C
      WRITE (LUPRI,160) I
  160 FORMAT (/,' Undefined int. coordinate type at No.',I3,
     &        /,1X,10('****'))
      GO TO 380
C
C     Print type
C
  170 CONTINUE
      IF (WRI) THEN
         IF (NEW) THEN
            WRITE (LUPRI,180) I,TYP,NAME1,NAME2,NAME3,NAME4,CCOEF,CSCAL
         ELSE
            WRITE (LUPRI,181)       NAME1,NAME2,NAME3,NAME4,CCOEF,CSCAL
         END IF
      END IF
  180 FORMAT (I6,'.',4X,A4,5X,4(A6,2X),F10.5,F12.6)
  181 FORMAT (20X,            4(A6,2X),F10.5,F12.6)
C
C     Test that all specified atoms are well defined
C
      NAB = NCORD/3
      IF (KA.LT.1 .OR. KA.GT.NAB .OR. KB.LT.1 .OR. KB.GT.NAB) GO TO 350
      IF (K.GT.2  .AND. (KC.LT.1 .OR. KC.GT.NAB)) GO TO 350
      IF (K.GT.3  .AND. (KD.LT.1 .OR. KD.GT.NAB)) GO TO 350
C
      GO TO (190,200,210,230,260,280,300), K
C
C..... stretch
C
  190 CALL VEKTOR (UU,R1,KA,KB,XA,YA,ZA)
      UU(1) =   UU(1)*ANULL
      UU(2) =   UU(2)*ANULL
      UU(3) =   UU(3)*ANULL
      VV(1) = - UU(1)
      VV(2) = - UU(2)
      VV(3) = - UU(3)
      IA(3) = 0
      IA(4) = 0
      QQ(I) = QQ(I) + R1*CCOEF
      CTYPE(I) = D1/XTANG
      GO TO 320
C
C.....inverse
C
  200 CALL VEKTOR (UU,R1,KA,KB,XA,YA,ZA)
      RM1   = D1/R1
      RM2   = RM1*RM1
      UU(1) = - RM2*UU(1)*ANULL
      UU(2) = - RM2*UU(2)*ANULL
      UU(3) = - RM2*UU(3)*ANULL
      VV(1) = - UU(1)
      VV(2) = - UU(2)
      VV(3) = - UU(3)
      IA(3) = 0
      IA(4) = 0
      QQ(I) = QQ(I) + RM1*CCOEF
      CTYPE(I) = XTANG
      GO TO 320
C
C.....bending
C
  210 CALL VEKTOR (U,R1,KA,KC,XA,YA,ZA)
      CALL VEKTOR (V,R2,KB,KC,XA,YA,ZA)
      CO = V3DOT(U,V)
      SI = S2(CO)
      DO 220 L = 1, 3
         UU(L) = (CO*U(L) - V(L))/(SI*R1)
         VV(L) = (CO*V(L) - U(L))/(SI*R2)
         WW(L) = - UU(L) - VV(L)
  220 CONTINUE
      IA(4) = 0
      QQ(I) = QQ(I) + CCOEF*DARCOS(CO)
      CTYPE(I) = D1
      GO TO 320
C
C.....out of plane
C
  230 CALL VEKTOR (U,R1,KA,KD,XA,YA,ZA)
      CALL VEKTOR (V,R2,KB,KD,XA,YA,ZA)
      CALL VEKTOR (W,R3,KC,KD,XA,YA,ZA)
      CALL V3NRML (V,W,Z)
      STETA = V3DOT(U,Z)
      CTETA = S2(STETA)
      CFI1  = V3DOT(V,W)
      SFI1  = S2(CFI1)
      CFI2  = V3DOT(W,U)
      CFI3  = V3DOT(V,U)
      DEN   = CTETA*SFI1**2
      ST2   = (CFI1*CFI2 - CFI3)/(R2*DEN)
      ST3   = (CFI1*CFI3 - CFI2)/(R3*DEN)
      DO 240 L = 1, 3
         VV(L) = Z(L)*ST2
         WW(L) = Z(L)*ST3
  240 CONTINUE
      CALL V3NRML (Z,U,X)
      CALL V3NRML (U,X,Z)
      DO 250 L = 1, 3
         UU(L) = Z(L)/R1
         ZZ(L) = - UU(L) - VV(L) - WW(L)
  250 CONTINUE
      CX = - CCOEF
      IF (STETA .LT. D0) CX = CCOEF
      QQ(I) = QQ(I) - CX*DARCOS(CTETA)
      CTYPE(I) = D1
      GO TO 320
C
C.....torsion
C
  260 CALL VEKTOR (U,R1,KA,KB,XA,YA,ZA)
      CALL VEKTOR (V,R2,KC,KB,XA,YA,ZA)
      CALL VEKTOR (W,R3,KC,KD,XA,YA,ZA)
      CALL V3NRML (U,V,Z)
      CALL V3NRML (W,V,X)
      CO  = V3DOT(U,V)
      CO2 = V3DOT(V,W)
      SI  = S2(CO)
      SI2 = S2(CO2)
      DO 270 L = 1, 3
         UU(L) = Z(L)/(R1*SI)
         ZZ(L) = X(L)/(R3*SI2)
         VV(L) = (R1*CO/R2 - D1)*UU(L) - R3*CO2/R2*ZZ(L)
         WW(L) = -UU(L) - VV(L) - ZZ(L)
  270 CONTINUE
      CO   = V3DOT(Z,X)
      U(1) = Z(2)*X(3) - Z(3)*X(2)
      U(2) = Z(3)*X(1) - Z(1)*X(3)
      U(3) = Z(1)*X(2) - Z(2)*X(1)
      SI3  = SQRT(U(1)**2 + U(2)**2 + U(3)**2)
      CO2  = V3DOT(U,V)
      S    = ARC1(-CO,SI3)
      IF (CO2 .LT. D0) S = - S
      IF (S .GT. (PI*0.5D0)) S = S - 2.0D0*PI
      QQ(I) = QQ(I) - CCOEF*S
      CTYPE(I) = D1
C
C     Remember that the range of this coordinate is -pi/2 to 3*pi/2
C     in order to shift the discontinuity off the planar position.
C
      GO TO 320
C
C.....linear coplanar bending
C
  280 CALL VEKTOR (U,R1,KA,KC,XA,YA,ZA)
      CALL VEKTOR (V,R2,KD,KC,XA,YA,ZA)
      CALL VEKTOR (X,R2,KB,KC,XA,YA,ZA)
      CO    = V3DOT(V,U)
      CO2   = V3DOT(X,V)
      QQ(I) = QQ(I) + CCOEF*(PI - DARCOS(CO) - DARCOS(CO2))
      CALL V3NRML (V,U,W)
      CALL V3NRML (U,W,Z)
      CALL V3NRML (X,V,W)
      CALL V3NRML (W,X,U)
C
C     Coordinate positiv if atom A moves towards atom D
C
      DO 290 L = 1, 3
         UU(L) = Z(L)/R1
         VV(L) = U(L)/R2
         WW(L) = - UU(L) - VV(L)
  290 CONTINUE
      IA(4) = 0
      CTYPE(I) = D1
      GO TO 320
C
C.....linear perpendicular bending
C
  300 CALL VEKTOR (U,R1,KA,KC,XA,YA,ZA)
      CALL VEKTOR (V,R2,KD,KC,XA,YA,ZA)
      CALL VEKTOR (Z,R2,KB,KC,XA,YA,ZA)
      CALL V3NRML (V,U,W)
      CALL V3NRML (Z,V,X)
      DO 310 L = 1, 3
         UU(L) = W(L)/R1
         VV(L) = X(L)/R2
         WW(L) = - UU(L) - VV(L)
  310 CONTINUE
      IA(4) = 0
      CO    = V3DOT(U,W)
      CO2   = V3DOT(Z,W)
      QQ(I) = QQ(I) + CCOEF*(PI - DARCOS(CO) - DARCOS(CO2))
      CTYPE(I) = D1
C
  320 CONTINUE
      IF (GETB) THEN
         DO 340 J = 1, 4
            M = IA(J)
            IF (M .GT. 0) THEN
               M = M - 1
               J1 = 3*(J - 1)
               DO 330 L = 1, 3
                  M1 = 3*M + L
                  L1 = J1 + L
                  BMAT(M1,I) = UV(L1)*CCOEF + BMAT(M1,I)
  330          CONTINUE
            END IF
  340    CONTINUE
      END IF
      GO TO 20
C
  350 CONTINUE
      WRITE (LUPRI,360) I
  360 FORMAT (/' Atoms erroneously defined, coordinate No.',I3,
     &        /1X,10('****'))
  370 CONTINUE
  380 RETURN
      END
C  /* Deck arc1 */
      FUNCTION ARC1 (X,Y)
#include <implicit.h>
#include <pi.h>
C
      IF (ABS(X) .LT. 1.0 D-11) THEN
         ARC1 = PI / 2
      ELSE
         S = ATAN(Y/X)
         IF (X .LT. 0.0D0) S = S + PI
         ARC1 = S
      END IF
      RETURN
      END
C  /* Deck darcos */
      FUNCTION DARCOS (X)
#include <implicit.h>
#include <pi.h>
C
      IF (X .GE. 1.0 D0 ) THEN
         DARCOS = 0.0D0
      ELSE IF (X .LE. -1.0 D0) THEN
         DARCOS = PI
      ELSE IF (ABS(X) .LT. 1.0 D-11) THEN
         DARCOS = PI / 2
      ELSE
         X1 = SQRT(1.0 D0 - X*X)
         S  = ATAN(X1/X)
         IF (X .LT. 0.0D0) S = S + PI
         DARCOS = S
      END IF
      RETURN
      END
C  /* Deck vektor */
      SUBROUTINE VEKTOR (U,R,I,J,XA,YA,ZA)
#include <implicit.h>
      DIMENSION U(3), XA(*), YA(*), ZA(*)
C
C     Bildet den normierten Entfernungsvektor vom Kern J nach Kern I
C     und die Entfernung R.
C
      U(1) = XA(I)-XA(J)
      U(2) = YA(I)-YA(J)
      U(3) = ZA(I)-ZA(J)
      R    = SQRT(V3DOT(U,U))
      CALL V3NOM (U)
      RETURN
      END
C  /* Deck v3nom */
      SUBROUTINE V3NOM (U)
#include <implicit.h>
      DIMENSION U(3)
      X = 1.0 D0 / SQRT(V3DOT(U,U))
      DO 10 I = 1,3
         U(I) = U(I)*X
   10 CONTINUE
      RETURN
      END
C  /* Deck s2 */
      FUNCTION S2 (X)
#include <implicit.h>
      S2 = SQRT(1.0 D0 - X*X)
      RETURN
      END
C  /* Deck v3dot */
      FUNCTION V3DOT (U,V)
#include <implicit.h>
      DIMENSION U(3), V(3)
      V3DOT = U(1)*V(1) + U(2)*V(2) + U(3)*V(3)
      RETURN
      END
C  /* Deck v3nrml */
      SUBROUTINE V3NRML (U,V,W)
#include <implicit.h>
      DIMENSION U(3), V(3), W(3)
C
C     99999...  W wird ein senkrecht auf die Ebene(U,V) stehender
C     Einheitor
C
      W(1) = U(2)*V(3)-U(3)*V(2)
      W(2) = U(3)*V(1)-U(1)*V(3)
      W(3) = U(1)*V(2)-U(2)*V(1)
      CALL V3NOM (W)
      RETURN
      END
C  /* Deck rpflin */
