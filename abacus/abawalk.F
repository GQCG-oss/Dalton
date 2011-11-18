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
C  FILE: abacus/abawalk.F
C
C  Geometry "walking" algorithms.
C  Implementation started in 1985 in Aarhus
C  by Trygve Helgaker and Hans Joergen Aa. Jensen.
C
C  Convention: routine names in this module are of the form 'WLK*'
C
C  /* Deck wlkinp */
      SUBROUTINE WLKINP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER ( NTABLE = 42, NWKTYP = 7, NSCTYP = 3,
     *            D0 = 0.0D0, D1 = 1.0D0 )
      PARAMETER ( TRSTI = 1.2D0, TRSTD = 0.7D0, TRSTDR = 0.5D0 )
C
      LOGICAL NEWDEF, SETRAD
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION IDOREP(0:7), IWKVEC(NWKTYP), ISCVEC(NSCTYP)
#include <abainf.h>
#include <gnrinf.h>
#include <nuclei.h>
#include <cbiwlk.h>
C
      CHARACTER*42 WLKTYP(7)
      DATA TABLE /'.RESTAR', '.PRINT ', '.INDEX ', '.TRUST ',
     &            '.ANHARM', '.DISPLA', '.SCALE ', '.MODE  ',
     &            '.TOLERA', '.HESFIL', '.NOPRED', '.MAXNUC',
     &            '.KEEPSY', '.REJECT', '.RATLIM', '.NEWTON',
     &            '.VIBCNV', '.GRDEXT', '.MODFOL', '.MASSES',
     &            '.REPS  ', '.ZERGRD', '.EIGEN ', '.DYNAMI',
     &            '.NOGRAD', '.MOMENT', '.NUMERI', '.HARMON',
     &            '.FRAGME', '.IMAGE ', '.XXXXXX', '.IRC   ',
     &            '.MAXTRU', '.NOORTH', '.NATCON', '.VIBAVE',
     &            '.NORMAL', '.ECKART', '.TEMPER', '.NO CEN',
     &            '.REUSE ', '.ACCURA'/
      DATA WLKTYP /'Level-shifted Newton walk (mode following)',
     *             'Gradient-extremal walk                    ',
     *             'Dynamic walk                              ',
     *             'Newton step                               ',
     *             'Eigenvector step                          ',
     *             'Numerical differentiation                 ',
     *             'Intrinsic reaction path                   '/
C
      CALL QENTER('WLKINP')
C
      SETRAD = .FALSE.
      NIP    = 0
      IRCSGN = 1
      TRUMAX = 0.10D0
      TRUMX1 = TRUMAX
      CALL IZERO(IWKVEC,NWKTYP)
      CALL IZERO(ISCVEC,NSCTYP)
C
      NEWDEF = (WORD .EQ. '*WALK  ')
      DXMXNU = XMXNUC
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &                      11,12,13,14,15,16,17,18,19,20,
     &                      21,22,23,24,25,26,27,28,29,30,
     &                      31,32,33,34,35,36,37,38,39,40,
     &                      41,42), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in *WALK input.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in *WALK input')
    1          CONTINUE
                  START = .FALSE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD,*) IPRWLK
                  IF (IPRWLK .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    3          CONTINUE
                  READ (LUCMD,*) IWKIND
               GO TO 100
    4          CONTINUE
                  READ (LUCMD,*) TRUSTR, TRUSTI, TRUSTD
                  SETRAD = .TRUE.
               GO TO 100
    5          CONTINUE
                  V3CAL = .TRUE.
                  IWKVEC(6) = 1
               GO TO 100
    6          CONTINUE
                  READ (LUCMD,*) DISPLC
               GO TO 100
    7          CONTINUE
                  ISCVEC(2) = 1
                  DO 7000 ICOOR = 1,3
                     DO 7000 INUC = 1,MXCENT
                        SCALCO(ICOOR,INUC) = D1
 7000             CONTINUE
                  READ (LUCMD,*) NUMNUC
                  DO 7100 INUC = 1,NUMNUC
                     READ (LUCMD,*) IATOM,(SCALCO(J,IATOM), J = 1,3)
 7100             CONTINUE
               GO TO 100
    8          CONTINUE
                  READ (LUCMD,*) IMODE
               GO TO 100
    9          CONTINUE
                  READ (LUCMD,*) TOLST
               GO TO 100
   10          CONTINUE
                  HESFLW = .TRUE.
               GO TO 100
   11          CONTINUE
                  WFPRED = .FALSE.
               GO TO 100
   12          CONTINUE
                  READ (LUCMD,*) XMXNUC
                  IF (XMXNUC .EQ. DXMXNU) ICHANG = ICHANG - 1
               GO TO 100
   13          CONTINUE
                  KEEPSY = .TRUE.
               GO TO 100
   14          CONTINUE
                  REJECT = .TRUE.
                  IF (.NOT.SETRAD) TRUSTD = TRSTDR
               GO TO 100
   15          CONTINUE
                  READ (LUCMD,*) RTRMIN, RTRGOD, REJMIN, REJMAX
               GO TO 100
   16          CONTINUE
                  IWKVEC(4) = 1
               GO TO 100
   17          CONTINUE
               GO TO 100
   18          CONTINUE
Cpeo mass-scaled coordinates as default
                  ISCVEC(3) = 1
                  IWKVEC(2) = 1
               GO TO 100
   19          CONTINUE
                  IWKVEC(1) = 1
               GO TO 100
   20          CONTINUE
                  ISCVEC(3) = 1
               GO TO 100
   21          CONTINUE
                  READ (LUCMD,*) NREPS
                  READ (LUCMD,*) (IDOREP(I),I=1,NREPS)
                  DO 21000 IREP = 1, 7
                     DOREPW(IREP) = .FALSE.
21000             CONTINUE
                  DO 21010 IREP = 1, NREPS
                     DOREPW(IDOREP(IREP)) = DOSYM(IREP + 1)
21010             CONTINUE
               GO TO 100
   22          CONTINUE
                  READ (LUCMD,*) ZERGRD
               GO TO 100
   23          CONTINUE
                  IWKVEC(5) = 1
                  GO TO 100
   24          CONTINUE
                  IWKVEC(3) = 1
                  GO TO 100
   25          CONTINUE
                  READ (LUCMD, *) NZEROG
                  READ (LUCMD, *) (IZEROG(I),I=1,NZEROG)
                  GO TO 100
   26          CONTINUE
                  READ (LUCMD,*) NSTMOM
                  DO 265 IP = 1, NSTMOM
                     READ (LUCMD, *) ISTMOM(IP), STRMOM(IP)
  265             CONTINUE
                  GO TO 100
   27          CONTINUE
                  IWKVEC(6) = 1
                  GO TO 100
   28          CONTINUE
                  READ (LUCMD,*) ANHFAC
                  GO TO 100
   29          CONTINUE
                  READ (LUCMD, *) NIP
                  READ (LUCMD, *) (IPART(IP), IP=1,NIP)
                  GO TO 100
   30          CONTINUE
                  IMAGE = .TRUE.
                  GO TO 100
   31          CONTINUE
                  GO TO 100
   32          CONTINUE
                  IWKVEC(7) = 1
                  READ (LUCMD,*) IRCSGN
                  GO TO 100
   33          CONTINUE
                  READ (LUCMD,*) TRUMX1
                  GO TO 100
   34          CONTINUE
                  NOORTH = .TRUE.
                  GO TO 100
   35          CONTINUE
                  NATCON = .TRUE.
                  GO TO 100
 36            CONTINUE
                  VIBAVE = .TRUE.
                  IWKVEC(6) = 1
                  GO TO 100
 37            CONTINUE
                  NMODIF = .TRUE.
                  GO TO 100
 38            CONTINUE
                  ECKART = .TRUE.
                  DO I = 1, NUCDEP
                     READ (LUCMD,'(7X,F17.10,2F24.10)') 
     &                    (ECKGEO(J,I), J = 1, 3)
                  END DO
                  GO TO 100
 39            CONTINUE
                  DOTEMP = .TRUE.
                  READ (LUCMD, *) NTEMP
                  IF (NTEMP .GT. MAXTMP) THEN
                     WRITE (LUPRI,'(/A/A)')
     &               'Too many temperatures requested in *WALK input.',
     &               'Increase MAXTMP in cbiwlk.h and recompile.'
                     CALL QUIT(
     &               'Too many temperatures requested in *WALK input')
                  END IF
                  READ (LUCMD, *) (TEMP(ITMP), ITMP=1,NTEMP)
                  GOTO 100
 40            CONTINUE
                  DOCENT = .FALSE.
                  GO TO 100
 41            CONTINUE
                  REUSED = .TRUE.
                  GO TO 100
 42            CONTINUE
                  ACCURT = .TRUE.
                  GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized for *WALK'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in *WALK input')
            END IF
      END IF
  300 CONTINUE
C
      IF (.NOT.DOWALK) THEN
         IF (ICHANG .GT. 0) WRITE (LUPRI,'(//A/)')
     &      ' Walk not requested, *WALK input ignored.'
         GO TO 9999
      END IF
C
      IF (.NOT.SETRAD .AND. IWKIND .GT. 0) THEN
         TRUSTR = 0.3D0
      END IF
      IF (.NOT.SETRAD .AND. IWKVEC(3) .EQ. 1) THEN
         TRUSTR = 0.005D0
         TRUSTI = 2.0D0
         TRUSTD = 0.8D0
      END IF
      IF (.NOT.SETRAD .AND. IWKVEC(7) .EQ. 1) THEN
         TRUSTR = 0.02D0
      END IF
      IF (IMAGE .AND. IMODE .EQ. 0) THEN
         IMODE = 1
      END IF
      IF (NMODIF .AND. .NOT. V3CAL) NMODIF = .FALSE.
C
      CALL HEADER('Changes of defaults for WALK:',0)
C
C     *******************************
C     ***** Determine walk type *****
C     *******************************
C
      NSET = ISUM(NWKTYP,IWKVEC,1)
C
C     a) No walk type requested, use defaults
C
      IF (NSET .EQ. 0) THEN
         IF (IWKTYP .EQ. 0) THEN
            IF (IWKIND .EQ. 0 .OR. IMAGE) THEN
               IWKTYP = 1
            ELSE
               IWKTYP = 2
            END IF
         END IF
C
C     b) Walk type explicitly requested
C
      ELSE IF (NSET .EQ .1) THEN
         DO 500 I = 1, NWKTYP
            IF (IWKVEC(I).GT.0) IWKTYP = I
  500    CONTINUE
C
C     c) Error: more than one type requested
C
      ELSE
         WRITE (LUPRI,'(//,2(A,/))')
     *      ' Inconsistent specification of walk type.',
     *      ' The following walks have been requested:'
         DO 510 I = 1, NWKTYP
            IF (IWKVEC(I).GT.0) WRITE (LUPRI,'(5X,A)') WLKTYP(I)
  510    CONTINUE
         CALL QUIT('ERROR in *WALK input.')
      END IF
      WRITE (LUPRI,'(/2A)') ' Walk type: ', WLKTYP(IWKTYP)
      IF (V3CAL) THEN
         WRITE (LUPRI,'(A)') ' Anharmonic force field '//
     &        'calculated using numerical differentiation'
         IF (NMODIF) WRITE(LUPRI,'(A)') ' Numerical differentiation '//
     &        'will be performed in normal coordinates'
         IF (DOTEMP) THEN
            WRITE (LUPRI,'(/A)') ' Effective geometries for the '//
     &           'following temperatures will be calculated (in K):'
            DO ITMP = 1, NTEMP
               WRITE (LUPRI,'(3X,F8.2)') TEMP(ITMP)
            END DO
            IF (DOCENT) THEN
               WRITE (LUPRI,'(A)') ' Centrifugal distortion effects '//
     &              'will be included'
            ELSE
               WRITE (LUPRI,'(A)') ' Centrifugal distortion effects '//
     &              'will NOT be included'
            END IF
         END IF
      END IF
      IF (VIBAVE) THEN
         WRITE (LUPRI,'(A)') ' Vibrational averaging of '//
     &        ' molecular properties at an effective geometry'
         IF (ECKART) THEN
            WRITE (LUPRI,'(A)') ' Molecular properties transformed to'
     &           //' the following Eckart frame'
            CALL PRIGEO(ECKGEO)
         END IF
         IF (DOTEMP) THEN
            WRITE (LUPRI,'(/A)') ' Vibrational corrections will be '//
     &           'calculated for the temperatures (in K):'
            DO ITMP = 1, NTEMP
               WRITE (LUPRI,'(3X,F8.2)') TEMP(ITMP)
            END DO
         END IF
         IF (REUSED) WRITE(LUPRI,'(A)') ' Property derivatives from '//
     &        'a previous geometry will be used'
         IF (ACCURT) WRITE(LUPRI,'(A)') ' A five-point derivative '//
     &        'scheme will be used for property derivatives'
         IF (HESFLW) WRITE (LUPRI,'(A)') ' The molecular Hessian '//
     &        ' be read from the file DALTON.HES'
      END IF
      IF (IMAGE) THEN
         WRITE (LUPRI,'(A)') ' Minimization of image surface '
         IF (START) WRITE (LUPRI,'(A,I2)') ' Initial image mode:',IMODE
ckr         IF (STRICT) WRITE (LUPRI,'(A)') ' Strict mode following.'
      END IF
      IF (IWKTYP .EQ.7) TRUSTR = MIN(TRUSTR,TRUMX1)
C
C     *******************
C     ***** Scaling *****
C     *******************
C
      NSET = ISUM(NSCTYP,ISCVEC,1)
C
C     a) No scale type requested, use defaults
C
      IF (NSET .EQ. 0) THEN
         ISCTYP = 1
         IF ((IWKTYP .EQ. 3) .OR. IWKTYP .EQ. 7) ISCTYP = 3
C
C     b) Input scaling
C
      ELSE IF (NSET .EQ .1) THEN
         DO 550 I = 1, NSCTYP
            IF (ISCVEC(I).GT.0) ISCTYP = I
  550    CONTINUE
         IF (IWKTYP .EQ. 3 .AND. ISCTYP .LT. 3) THEN
            WRITE (LUPRI,'(//,2(A,/))')
     *         ' Inconsistent specification of scaling.',
     *         ' Only mass scaling allowed for dynamic walks.'
            CALL QUIT('ERROR in *WALK input')
         END IF
C
C     c) Error: more than one type requested
C
      ELSE
         WRITE (LUPRI,'(//,2(A,/))')
     *      ' Inconsistent specification of scaling.',
     *      ' The following scaling schemes have been requested:'
         IF (ISCVEC(1).GT.0) WRITE (LUPRI,'(/A)') ' No scaling.'
         IF (ISCVEC(2).GT.0) WRITE (LUPRI,'(/A)') ' Scaling from input.'
         IF (ISCVEC(3).GT.0) WRITE (LUPRI,'(/A)') ' Mass scaling.'
         CALL QUIT('ERROR in *WALK input')
      END IF
      IF (WFPRED) GDALL = .TRUE.
C
      IF (ICHANG .GT. 0) THEN
C
         IF (IWKTYP .LT. 3) THEN
            WRITE (LUPRI,'(/A,I2/)')
     *           ' Index of critical point searched:',IWKIND
            WRITE (LUPRI,'(A/A,F12.8/)')
     *        ' Tolerance for stationary point :',
     *        ' Molecular gradient norm less than',TOLST
            IF (IMODE .NE. 0) WRITE (LUPRI,'(A,I10)')
     *           ' Hessian mode selection :',IMODE
            IF (KEEPSY) WRITE (LUPRI,'(/A)')
     *           ' Symmetry will not be broken.'
            IF (ANHFAC .NE. 100.D0) WRITE (LUPRI,'(/A)')
     *         ' Factor for harmonic dominance :', ANHFAC
         ELSE IF (IWKTYP .EQ. 3) THEN
            IF (NIP .GT. 0) WRITE (LUPRI,'(/A,20I3)')
     *         ' Fragment labels for atoms :', (IPART(IP), IP=1,NIP)
            IF (START) THEN
               IF (NSTMOM.EQ.0) THEN
                  WRITE (LUPRI,'(/A,/A)')
     &               ' No initial momenta have been given.',
     &               ' All momenta are set to zero.'
               ELSE
                  WRITE (LUPRI,'(/A,I3,A/)')
     &               ' Initial momenta have been given for ',NSTMOM,
     &               ' modes:'
                  DO 150 I = 1, NSTMOM
                    WRITE(LUPRI,'(10X,I5,5X,F12.6)')ISTMOM(I),STRMOM(I)
  150             CONTINUE
               END IF
            END IF
         ELSE IF (IWKTYP .EQ. 7) THEN
            IF (START) THEN
               IF (IRCSGN .EQ. 0) THEN
                  WRITE (LUPRI,'(//,A)')
     &               ' Error in *WALK input, IRCSGN equal to zero.'
                     CALL QUIT('IRCSGN equal to zero.')
               END IF
               IMODE = SIGN(1,IRCSGN)
               WRITE (LUPRI,'(/A,SP,I2)')
     &            ' Sign of first step: ', IMODE
            END IF
            IF (TRUMAX .NE. TRUMX1) THEN
               TRUMAX = TRUMX1
               WRITE (LUPRI,'(/A,F12.4)')
     &           ' Maximum arc length of each step ',TRUMX1
            END IF
         END IF
C
         IF (IPRWLK .NE. IPRDEF) THEN
            WRITE (LUPRI,'(A,I10)') ' Print level in WALK    :',IPRWLK
         END IF
C
         IF (IWKTYP.NE.5 .AND. IWKTYP.NE.6) THEN
            WRITE (LUPRI,'(/A/A/,(A,F10.4))')
     *        ' Restricted step control parameters',
     *        ' ----------------------------------',
     *        ' Initial trust radius   :',TRUSTR,
     *        ' Trust radius increment :',TRUSTI,
     *        ' Trust radius decrement :',TRUSTD,
     *        ' Bad prediction ratio   :',RTRMIN,
     *        ' Good prediction ratio  :',RTRGOD,
     *        ' Rejection ratio, low   :',REJMIN,
     *        ' Rejection ratio, high  :',REJMAX
            IF (XMXNUC .NE. DXMXNU) WRITE (LUPRI,'(/A,F10.4)')
     *        ' Maximum movement of any atom  :',XMXNUC
         END IF
         IF (REJECT) THEN
            IF (.NOT.START) THEN
               WRITE (LUPRI,'(/A/A)')
     *            ' Previous step has been rejected,',
     *            ' and the trust radius read in will be reduced.'
            ELSE
               WRITE (LUPRI,'(/A)')
     *     ' .REJECT has only meaning when .RESTART is specified.'
            END IF
         END IF
         IF (IWKTYP .NE. 6) THEN
            IF (ISCTYP .EQ. 1) THEN
               WRITE (LUPRI,'(/A/)') ' No scaling in this run.'
            ELSE IF (ISCTYP .EQ. 2) THEN
               WRITE (LUPRI,'(A)') ' Scaling as specified in input:'
               DO 600 INUC = 1,NUMNUC
                  WRITE (LUPRI,'(3F20.8)') (SCALCO(J,INUC),J=1,3)
  600          CONTINUE
            ELSE IF (ISCTYP .EQ. 3) THEN
               WRITE (LUPRI,'(A)') ' Mass scaling in this run.'
            END IF
         END IF
         IF (IWKTYP .EQ. 6) THEN
            IF (ITERMX .EQ. 20) ITERMX = 6*NUCDEP + 1
            WRITE (LUPRI,'(/A,I5)') 'Maximum number of iterations '//
     *           'reset in numerical differentiation to ',ITERMX
         END IF
         IF (NZEROG .GT. 0) THEN
            WRITE (LUPRI,'(/A,16I5)')
     *         ' The following gradient elements are set to zero: ',
     *               (IZEROG(I),I = 1, NZEROG)
         END IF
         IF (.NOT. WFPRED) THEN
            WRITE (LUPRI,'(/A)') ' No prediction of wave function.'
         END IF
         IF (NOORTH) THEN
            WRITE (LUPRI,'(/A)')
     &         ' Predicted orbitals are not orthogonalized.'
         END IF
         IF (NATCON) THEN
            WRITE (LUPRI,'(/A)')
     &            ' Natural connection used for prediction '//
     &                        '(symmetric is default).'
         END IF
         IF (.NOT. START) THEN
            WRITE (LUPRI,'(/A)') ' Information from old LUWLK read.'
         END IF
         WRITE (LUPRI,'(/)')
      END IF
 9999 CALL QEXIT('WLKINP')
      RETURN
      END
C  /* Deck wlkini */
      SUBROUTINE WLKINI
C
C     Initialize /CBIWLK/
C
#include <implicit.h>
#include <mxcent.h>
#include <abainf.h>
#include <gnrinf.h>
#include <cbiwlk.h>
C
      IF (DOWALK) THEN
         IWKTYP = 1
         GDALL = .TRUE.
         ISCTYP = 1
      ELSE
         IWKTYP = 0
         ISCTYP = 0
      END IF
      TOLST  = 1.D-4
      TRUSTR = 0.5D0
      STMOM  = 1.D-4
      TRUSTI = 1.2D0
      TRUSTD = 0.7D0
      XMXNUC = 0.50D0
      REJMIN = 0.1D0
      REJMAX = 1.9D0
      RTRMIN = 0.4D0
      RTRGOD = 0.8D0
      DISPLC = 1.D-4
      THRLDP = 1.D-4
      IPRWLK = IPRDEF
      IMAGE  = .FALSE.
ckr      STRICT = .FALSE.
      V3CAL  = .FALSE.
      VIBAVE = .FALSE.
      NMODIF = .FALSE.
      ECKART = .FALSE.
      DOTEMP = .FALSE.
      DOCENT = .TRUE.
      REUSED = .FALSE.
      CALL DZERO(ECKGEO,MXCOOR)
      NSTMOM = 0
      IWKIND = 0
      IMODE  = 0
      ANHFAC = 100.D0
      CALL IZERO(IPART,MXCENT)
C
C     Commented out because it has to be set externally.
C
C      START  = .FALSE.
C
      DO 400 I = 0, 7
         DOREPW(I) = DOSYM(I + 1)
  400 CONTINUE
      WFPRED = .TRUE.
      REJECT = WLKREJ
      IF (REJECT) TRUSTD = 0.5D0
      KEEPSY = .FALSE.
      NZEROG = 0
      ZERGRD = 1.D-5
  410 CONTINUE
      NTEMP = 0
C
      NOORTH = .FALSE.
      NATCON = .FALSE.
      RETURN
      END
C  /* Deck wlkdrv */
      SUBROUTINE WLKDRV(POLDD,POLDQ,POLDL,POLDA,SPNTOT,WORK,LWORK)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
C
      DIMENSION POLDD(*), POLDQ(*), POLDL(*), POLDA(*), SPNTOT(*),
     &          WORK(LWORK)
C
#include <cbiwlk.h>
#include <cbisol.h>
#include <cbilnr.h>
#include <nuclei.h>
#include <symmet.h>
#include <trkoor.h>
#include <prkoor.h>
C
      CALL QENTER('WLKDRV')
      IF (IPRWLK .GT. 3) CALL TIMER ('START ',TIMEIN,TIMOUT)
      IF (IPRWLK .GT. 2) CALL TITLER('Output from WLKDRV','*',103)
      IF ((V3CAL .OR. VIBAVE) .AND. (MAXREP .GT. 0)) THEN
         WRITE (LUPRI,'(/A)') ' Program cannot use symmetry when '//
     &        'doing vibrational averaging'
         CALL QUIT('Symmetry not implemented for vibrational '//
     &        'averaging')
      END IF
C
C     NPRREP(IREP): Number of coordinates to be projected out
C
      DO 50 IREP = 0, MAXREP
         IF (SOLVNT) THEN
            NPRREP(IREP) = NTRREP(IREP) + NAXREP(IREP,1)
         ELSE
            NPRREP(IREP) = NTRREP(IREP)
         END IF
   50 CONTINUE
C
      N2CRT  = 0
      NCRTOT = 0
      NTMAT  = 0
      DO 100 IREP = 0, MAXREP
         NCR = NCRREP(IREP,1)
         NPR = NPRREP(IREP)
         DOREPW(IREP) = DOREPW(IREP) .AND. (NCR .GT. NPR)
         IF (IPRWLK .GT. 5) THEN
            WRITE (LUPRI,'(/A,I5)') ' Symmetry',IREP + 1
            WRITE (LUPRI,'( A,I5)') ' Cartesian coordinates',NCR
            WRITE (LUPRI,'( A,I5)') ' External coordinates ',NPR
            WRITE (LUPRI,'( A,I5)') ' Internal coordinates ',NCR - NPR
         END IF
         IF (DOREPW(IREP)) THEN
            NCRTOT = NCRTOT + NCR
            N2CRT  = N2CRT  + NCR*NCR
            NTMAT  = NTMAT  + NPR*NCR
         END IF
  100 CONTINUE
      NCART  = NCRREP(0,1)
      NCRIND = 3*NUCIND
      IF (IPRWLK .GT. 5) THEN
         WRITE (LUPRI,'(/A,I5)')
     *         ' Total number of Cartesian coordinates:', NCRTOT
      END IF
      KGRD  = 1
      KHES  = KGRD  + NCRTOT
      KEVAL = KHES  + N2CRT
      KEVEC = KEVAL + NCRTOT
      KGRDI = KEVEC + N2CRT
      KTMAT = KGRDI + NCRTOT
      KCOR  = KTMAT + NTMAT
      KSTPS = KCOR  + NCRIND
      KSTPC = KSTPS + NCART
      KSTPD = KSTPC + NCRIND
      KSCAL = KSTPD + NCRTOT
      KEVLX = KSCAL + NCRTOT
      KEVCX = KEVLX + NCART
      KPMD  = KEVCX + NCART*NCART
      KPMDX = KPMD  + NCRTOT
      KPMCR = KPMDX + NCRTOT
      KPMCX = KPMCR + NCRTOT
      KOVLP = KPMCX + NCRTOT
      KGRD0 = KOVLP + NCART*NCART
      KHES0 = KGRD0 + NCART
      KCOR0 = KHES0 + NCART*NCART
      KGRDN = KCOR0 + NCRIND
      KHESN = KGRDN + NCART
      KSTPX = KHESN + NCART*NCART
      KDPG0 = KSTPX + NCART
      KDPGF = KDPG0 + 3*NCART
      KDIPM = KDPGF + 3*NCART
      KDIP0 = KDIPM + 3
      KSTDX = KDIP0 + 3
      KGRDX = KSTDX + NCART
      KTRDI = KGRDX + NCART
      KANHR = KTRDI + NCART
      KEDIA = KANHR + NCART
      KAFD  = KEDIA + NCRTOT
      KGND  = KAFD  + 9*MXFR*MXCOOR
      KGLD  = KGND  + 9*MXFR*MXCOOR
      KAD   = KGLD  + 9*MXFR*MXCOOR
      KAFU  = KAD   + 27*MXFR*MXCOOR
      KGNU  = KAFU  + 9*MXFR
      KGLU  = KGNU  + 9*MXFR
      KAU   = KGLU  + 9*MXFR
      KFMATF= KAU   + 27*MXFR
      IF (V3CAL) THEN
         KDIPMF = KFMATF+ NCART*NCART*NCART
      ELSE
         KDIPMF = KFMATF
      END IF
      IF (VIBAVE) THEN
         KSUST0 = KDIPMF
         KSUSTF = KSUST0 + 9
         KGFAC0 = KSUSTF + 9*NCART
         KGFACF = KGFAC0 + 9
         KQUAD0 = KGFACF + 9*NCART
         KQUADF = KQUAD0 + 9
         KQUADT = KQUADF + 9*NCART
         KSIGM0 = KQUADT + 9
         KSIGMF = KSIGM0 + 3*MXCOOR
         KSRC0  = KSIGMF + 3*MXCOOR*NCART
         KSRCF  = KSRC0  + 3*MXCOOR
         KSIGMT = KSRCF  + 9*MXCENT*NCART
         KCSTRA = KSIGMT + 9*MXCENT
         KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
         KPOLA0 = KSCTRA + 9*NUCDEP*NUCDEP
         KPOLAF = KPOLA0 + 9
         KEFG0  = KPOLAF + 9*NCART
         KEFGF  = KEFG0  + 9*MXCENT
         KSSJ0  = KEFGF  + 9*MXCENT*NCART
         KSSJF  = KSSJ0  + MXCOOR*MXCOOR
         KWRK   = KSSJF  + MXCOOR*MXCOOR*NCART
      ELSE
         KSUST0 = KDIPMF
         KSUSTF = KDIPMF
         KGFAC0 = KDIPMF
         KGFACF = KDIPMF
         KQUAD0 = KDIPMF
         KQUADF = KDIPMF
         KQUADT = KDIPMF
         KSIGM0 = KDIPMF
         KSIGMF = KDIPMF
         KSIGMT = KDIPMF
         KCSTRA = KDIPMF
         KSCTRA = KDIPMF
         KSRC0  = KDIPMF
         KSRCF  = KDIPMF
         KPOLA0 = KDIPMF
         KPOLAF = KDIPMF
         KEFG0  = KDIPMF
         KEFGF  = KDIPMF
         KSSJ0  = KDIPMF
         KSSJF  = KDIPMF
         KWRK   = KDIPMF
      END IF
      LWRK  = LWORK - KWRK + 1
      IF (KWRK .GT. LWORK) CALL STOPIT('WLKDRV',' ',KWRK,LWORK)
      CALL WLKDR1(POLDD,POLDQ,POLDL,POLDA,WORK(KGRD),WORK(KHES),
     &            WORK(KEVAL),WORK(KEVEC),WORK(KGRDI),WORK(KTMAT),
     &            WORK(KCOR),WORK(KSTPS),WORK(KSTPC),WORK(KSTPD),
     &            WORK(KSCAL),WORK(KEVCX),WORK(KPMD),WORK(KPMDX),
     &            WORK(KPMCR),WORK(KPMCX),WORK(KOVLP),WORK(KGRD0),
     &            WORK(KHES0),WORK(KCOR0),WORK(KGRDN),WORK(KHESN),
     &            WORK(KEVLX),WORK(KSTPX),WORK(KDPG0),WORK(KDPGF),
     &            WORK(KDIPM),WORK(KDIP0),WORK(KSTDX),
     &            WORK(KGRDX),WORK(KTRDI),WORK(KANHR),WORK(KEDIA),
     &            WORK(KAFD),WORK(KGND),WORK(KGLD),WORK(KAD),WORK(KAFU),
     &            WORK(KGNU),WORK(KGLU),WORK(KAU),WORK(KFMATF),
     &            WORK(KSUST0),WORK(KSUSTF),WORK(KGFAC0),WORK(KGFACF),
     &            WORK(KQUAD0),WORK(KQUADF),WORK(KQUADT),WORK(KSIGM0),
     &            WORK(KSIGMF),WORK(KSIGMT),WORK(KCSTRA),WORK(KSCTRA),
     &            WORK(KSRC0),WORK(KSRCF),WORK(KPOLA0),WORK(KPOLAF),
     &            WORK(KEFG0),WORK(KEFGF),SPNTOT,
     &            WORK(KSSJ0),WORK(KSSJF),WORK(KWRK),
     &            LWRK,NCRIND,NCRTOT,N2CRT,NTMAT,NCART)
      IF (IPRWLK.GT.3) CALL TIMER ('WLKDRV',TIMEIN,TIMOUT)
      CALL QEXIT('WLKDRV')
      RETURN
      END
C  /* Deck wlkdr1 */
      SUBROUTINE WLKDR1(POLDD,POLDQ,POLDL,POLDA,GRDCAR,HESCAR,EVAL,EVEC,
     &                  GRDDIA,TMAT,COOR,STPSYM,STPCAR,STPDIA,SCAL,
     &                  EVECX,PMDIA,PMDIAX,PMCAR,PMCARX,OVLPMA,GRAD0,
     &                  HESS0,COOR0,GRADF,HESSF,EVALX,STPSYX,DIPG0,
     &                  DIPGF,DIPM,DIPM0,STPDIX,GRDDIX,TRSDIA,
     &                  ANHARM,ERGDIA,ROAAFD,ROAGND,ROAGLD,ROAAD,ROAAFU,
     &                  ROAGNU,ROAGLU,ROAAU,FMATF,SUSTO0,SUSTF,GFAC0,
     &                  GFACF,QUAD0,QUADF,QUADT,SIGMA0,SIGMAF,SIGMAT,
     &                  CSTRA,SCTRA,SRC0,SRCF,POLAR0,POLARF,
     &                  EFG0,EFGF,SPNTOT,SSJ0,SSJF,WORK,LWORK,NCRIND,
     &                  NCRTOT,N2CRT,NTMAT,NCART)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (D0 = 0.0D00, D100 = 100.0D0, SMALL = 1.0D-8)
C
      LOGICAL INDXOK, OLD, STOP, STATPO, WLKEND, OLDREA, DOHESS, DODIP,
     &        PRJTRO, DONWTN, NOTALL
C
      DIMENSION POLDD(*), POLDQ(*), POLDL(*), POLDA(*)
      DIMENSION GRDCAR(NCRTOT), HESCAR(N2CRT), EVAL(NCRTOT),
     *          EVEC(N2CRT), TMAT(NTMAT), COOR(NCRIND), STPCAR(NCRIND),
     *          WORK(LWORK), CNDHES(0:7), INDHES(0:7), GRDDIA(NCRTOT),
     *          SCAL(NCRTOT), STPSYM(NCART), STPDIA(NCRTOT),
     *          EVECX(NCART,NCART), PMDIA(NCRTOT), PMCAR(NCRTOT),
     *          PMCARX(NCRTOT), OVLPMA(NCART,NCART), GRAD0(NCART),
     *          HESS0(NCART,NCART), COOR0(NCRIND), GRADF(NCART),
     *          HESSF(NCART,NCART), EVALX(NCART), STPSYX(NCRIND),
     *          DIPG0(3,NCART), DIPGF(3,NCART), DIPM(3), DIPM0(3),
     *          PMDIAX(NCRTOT), STPDIX(NCART), GRDDIX(NCART),
     *          TRSDIA(NCART), ANHARM(NCART), ROAAFD(*), ROAGND(*),
     *          ROAGLD(*), ROAAD(*), ROAAFU(*), ROAGNU(*), ROAGLU(*),
     *          ROAAU(*), FMATF(*), SUSTO0(*), SUSTF(*), GFAC0(*),
     *          GFACF(*), QUAD0(*), QUADF(*), QUADT(*), SIGMA0(*),
     *          SIGMAF(*), SIGMAT(*), CSTRA(*), SCTRA(*), SRC0(*),
     *          SRCF(*), POLAR0(*), POLARF(*),
     *          EFG0(*), EFGF(*), SPNTOT(*), SSJ0(*), SSJF(*)
C
#include <cbiwlk.h>
#include <abainf.h>
#include <nuclei.h>
#include <symmet.h>
#include <trkoor.h>
#include <prkoor.h>
#include <taymol.h>
C
C
      IF (IPRWLK .GT. 5) CALL TITLER('Output from WLKDR1','*',103)
C
      OLD = .NOT.START
      ISGNMOD = SIGN(1,IMODE)
      SGNMOD  = ISGNMOD
      IMODE = ABS(IMODE)
C
C     ***************************
C     ***** Initialize walk *****
C     ***************************
C
      IF (START) THEN
         IPOINT = 0
         ESTART = ERGMOL
         IF (IWKTYP .EQ. 3) THEN
            TACCUM = D0
            TIMSTP = D0 ! dynamically reset later
            CALL DZERO(PMCAR,NCRTOT)
         ELSE IF (IWKTYP .EQ. 7) THEN
            AACCUM = D0
         END IF
      END IF
C
C     ***********************************************
C     ***** Read information from previous step *****
C     ***********************************************
C
      OLDREA = (OLD .AND. (IWKTYP.LT.4 .OR. IWKTYP.GE.6))
      IF (VIBAVE .AND. REUSED) OLDREA = .TRUE.
      IF (OLDREA) THEN
         CALL WLKREA(EVALX,EVECX,PMCARX,GRAD0,HESS0,DIPM0,DIPG0,COOR0,
     &               GRADF,HESSF,DIPGF,ERGMLX,GRDNRX,ERGSCX,TRUSTX,
     &               ESTART,TACCUM,DISPLC,STPSYX,STPNRX,STPDIX,GRDDIX,
     &               EKIN0,TRUMAX,NCRTOT,N2CRT,NCRIND,NCART,IWKTYP,
     &               IWKIND,IMODEX,ISCTYP,IPOINT,ISTIND,ISTTPX,NXYZ,
     &               DODIP,IMAGE,AACCUM,ROAAFD,ROAGND,ROAGLD,
     &               ROAAD,ROAAFU,ROAGNU,ROAGLU,ROAAU,FMATF,VIBAVE,
     &               SUSTO0,SUSTF,GFAC0,GFACF,QUAD0,QUADF,SIGMA0,SIGMAF,
     &               SRC0,SRCF,POLAR0,POLARF,EFG0,EFGF,SSJ0,
     &               SSJF,NUMMOD,IPRWLK)
         IMODE  = IMODEX
         TRUSTR = TRUSTX
         IF ((VIBAVE .AND. OLD) .OR. NMODIF) THEN
            CALL DCOPY(NCART*NCART,EVECX,1,EVEC,1)
            CALL DCOPY(NCART,EVALX,1,EVAL,1)
         ELSE IF (VIBAVE) THEN
C
C     Need to reset a couple of variables read in from WLKREA
C
            IPOINT = 0
            ESTART = ERGMOL
            CALL DCOPY(NCART*NCART,EVECX,1,EVEC,1)
            CALL DCOPY(NCART,EVALX,1,EVAL,1)
         END IF
      ELSE
         TRUSTX = TRUSTR
         IF (DIPDER .AND. (IWKTYP .EQ. 6 .AND. .NOT. 
     &                    (V3CAL .OR. VIBAVE))) THEN
            NXYZ = 0
            DO 100 I = 1, 3
               IF (ISYMAX(I,1) .EQ. 0) NXYZ = NXYZ + 1
  100       CONTINUE
            DODIP = NXYZ .GT. 0
         ELSE
            NXYZ = 0
            DODIP = .FALSE.
         END IF
      END IF
      IF (IPRWLK .GT. 2) WRITE (LUPRI,'(A,I5)') ' Walk type: ',IWKTYP
      PRJTRO = (IWKTYP .NE. 3) .AND. (IWKTYP .NE. 7)
      DOHESS = .NOT. (IWKTYP .EQ. 6 .AND. OLD)
      IF (V3CAL .AND. .NOT. NMODIF) DOHESS = .TRUE.
C
C     *******************************************************
C     ***** Cartesian coordinates, gradient and Hessian *****
C     *******************************************************
C
      CALL WLKCGH(GRDCAR,HESCAR,DIPM,DIPG0,WORK,LWORK,DOREPW,DOHESS,
     *            DODIP,NCRTOT,N2CRT,NXYZ,NCART,START,IPRWLK)
C
C     *******************************************
C     ***** Diagonal Hessian representation *****
C     *******************************************
C
      IF (IWKTYP .NE. 6) THEN
         CALL DZERO(STPDIA,NCRTOT)
         THRIND = D0
         IF (IWKTYP .EQ. 7) THRIND = SMALL
         CALL WLKDIA(GRDCAR,HESCAR,EVAL,EVEC,GRDDIA,TMAT,SCAL,WORK,
     &               THRLDP,THRIND,CNDHES,INDHES,LWORK,IZEROG,NZEROG,
     &               DOREPW,NCRTOT,N2CRT,NTMAT,ISCTYP,IWKTYP,PRJTRO,
     &               IPRWLK)
      END IF
      IF (IWKTYP .EQ. 6 .AND. ((VIBAVE .AND. .NOT. REUSED) .OR. NMODIF)
     &     .AND. START) THEN
C
C        Generate normal coordinates for numerical differentiation
C
         KAMASS = 1
         KNATYP = KAMASS + NUCDEP
         KNUMIS = KNATYP + (NUCDEP + 1)/IRAT
         KORIG  = KNUMIS + (NUCDEP + 1)/IRAT
         KGEOM  = KORIG  + 3
         KDKIN  = KGEOM  + 3*NUCDEP
         KHESMW = KDKIN  + NCART
         KEVECS = KHESMW + NCART*(NCART + 1)/2
         KTRROV = KHESMW + NCART*NCART
         KWRK1  = KTRROV + 6*NCART
         LWRK1  = LWORK  - KWRK1 + 1
         IF (LWRK1 .LE. 0) CALL STOPIT('WLKVIB',' ',KWRK1,LWORK)
C
         IF (HESFLW) THEN
            LUHES = -90
            CALL DZERO(HESCAR,NCART*NCART)
            INQUIRE(FILE='DALTON.HES')
            CALL GPOPEN(LUHES,'DALTON.HES','OLD',' ','FORMATTED',
     &                  IDUMMY,.FALSE.)
            READ(LUHES,*) IDIM
            READ(LUHES,*)
            IF (IDIM .NE. NCART) CALL QUIT('The Hessian in DALTON.HES '
     &         //'has wrong dimensions.')
            IPOS = 1
            DO 10 I = 1, NCART
               DO 15 J = 1, NCART
                  READ(LUHES,*) HESCAR(IPOS)
                  IPOS = IPOS + 1
 15            CONTINUE
               READ(LUHES,*)
 10         CONTINUE
            CALL DCOPY(NCART*NCART,HESCAR,1,HESS0,1)
C
C     We try to see if the user has supplied some coordinates on the
C     file that we should use for updating the coordinate information
C     Also check if they match coordinates from MOLECULE.INP
C
            IDIFF = 0
            CALL GPCLOSE(LUHES,'KEEP')
         END IF
         CALL CMMASS(WORK(KGEOM),WORK(KAMASS),WORK(KNATYP),WORK(KNUMIS),
     &               IPRWLK)
C     
C     Project out translation and rotational motion from the Hessian
C
         CALL VIBHES(IPRWLK,NCART,WORK(KGEOM),KTRRO,GRDCAR,HESCAR,
     &               WORK(KTRROV),WORK(KNATYP),WORK(KWRK1))
         NOTALL = .FALSE.
ckr:  ... assume full hessian for all symmetries has been calculated.
         CALL VIBNOR(HESCAR,WORK(KAMASS),WORK(KDKIN),WORK(KHESMW),
     &               EVAL,EVEC,WORK(KEVECS),WORK(KWRK1),LWRK1,NCART,
     &               NUMMOD,NNEG,NOTALL,IPRWLK)
      END IF
C
C     **************************************
C     ***** Kinetic energy and momenta *****
C     **************************************
C
      IF (IWKTYP .EQ. 3) THEN
         IF (START) THEN
C
C           Initial momenta - diagonal and Cartesian basis
C
            CALL DZERO(PMDIAX,NCART)
            DO 150 I = 1, NSTMOM
               PMDIAX(ISTMOM(I)) = STRMOM(I)
  150       CONTINUE
            CALL WLKPCR(PMDIAX,EVEC,PMCARX,SCAL,NCART,IPRWLK)
         ELSE
C
C           Transform old momenta to diagonal basis
C
            CALL WLKPDI(PMCARX,PMDIAX,EVEC,SCAL,NCRTOT,NCART,IPRWLK)
         END IF
C
C        Total kinetic energy
C
         CALL WLKKIN(EKINTO,PMCARX,SCAL,NCART,IPRWLK)
         IF (START) EKIN0 = EKINTO
C
C        Relative translational energy
C
         CALL WLKTRA(PMCARX,EKINPA,VREL,SCAL,CMDIF,NCART,IPRWLK)
      END IF
C
C     *********************************
C     ***** Check for convergence *****
C     *********************************
C
      STOP = WLKEND(GRDCAR,GRDNRM,NCRREP(0,1),INDHES,INDEX,INDXOK,
     *              STATPO,ISTIND,IPOINT,NUMMOD)
C
C     a) Stop walk
C     ============
C
      IF (STOP) THEN
         IPOINT = IPOINT + 1
         IF (IWKTYP .EQ. 6) THEN
            IF (VIBAVE .OR. NMODIF) THEN
               CALL WLKNCI(STPCAR,COOR,DISPLC,NCART,NCRIND,
     &                     EVEC,IPOINT,NUMMOD,ACCURT,IPRWLK)
            ELSE
               CALL WLKINC(STPSYM,STPCAR,COOR,DISPLC,NCART,NCRIND,
     *                     IPOINT,IPRWLK)
            END IF
            CALL WLKMOL(COOR)
            ERGSEC = 0.0D0
            IF (WFPRED) CALL WLKSAV(ERGSEC,COOR,WORK,LWORK)
         END IF
C
C     b) Continue walk
C     ================
C
      ELSE
         IF (PRJTRO) THEN
            NINTER = NCRREP(0,1) - NPRREP(0)
         ELSE
            NINTER = NCRREP(0,1)
         END IF
C
C        ******************************
C        ***** Identify walk mode *****
C        ******************************
C
         IF (OLD.AND.(((IWKIND.GT.0).AND.(IWKTYP.LT.3)).OR.IMAGE)) THEN
            CALL WLKMOD(EVAL,EVALX,EVEC,EVECX,GRDDIA,OVLPMA,OVLPWM,
     &                  COFMOD,ANHRAD,TRUSTX,NCART,NINTER,IMODE,IMAGE,
     &                  IPRWLK)
         ELSE IF (IWKTYP .EQ. 7) THEN
            IMODE = 1
            ANHRAD = D100
         ELSE
            ANHRAD = D100
         END IF
C
C        ************************************
C        ***** Estimate anharmonicities *****
C        ************************************
C
         IF ((IWKTYP .EQ. 3) .OR. (IWKTYP .EQ. 7)) THEN
            CALL WLKANH(ANHFAC,GRDDIA,EVAL,EVALX,EVEC,EVECX,STPDIX,
     &                  OVLPMA,ANHARM,TRSDIA,NCART,OLD,IPRWLK)
         END IF
C
C        *******************************
C        ***** Update trust region *****
C        *******************************
C
         IF (OLD .AND. (IWKTYP .LT. 4 .OR. IWKTYP .EQ. 7)) THEN
            CALL WLKRAD(ERGMOL,ERGMLX,ERGSCX,STPNRX,OVLPWM,COFMOD,
     *                  RATIO,ISTTPX,IPRWLK)
         END IF
C
C        ************************
C        ***** New geometry *****
C        ************************
C
         IF (START) CALL DCOPY(NCRIND,CORD,1,COOR0,1)
         IF (.NOT.(START .OR. REJECT)) IPOINT = IPOINT + 1
         IF (IWKTYP .EQ. 6) THEN
C
C     Normal coordinate or Cartesian displacements
C
            IF (VIBAVE .OR. NMODIF) THEN
               CALL WLKNCI(STPCAR,COOR,DISPLC,NCART,NCRIND,
     &                     EVEC,IPOINT,NUMMOD,ACCURT,IPRWLK)
            ELSE
               CALL WLKINC(STPSYM,STPCAR,COOR,DISPLC,NCART,NCRIND,
     *                     IPOINT,IPRWLK)
            END IF
C
C     The step should be so small that no important change in energy 
C     should occur, and prediction of the new energy should be zero energy
C     change. KR, July-97
C
            ERGSEC = 0.0D0
         ELSE
            STPNRM = TRUSTR
            CALL WLKGEO(STPNRM,GRDDIA,EVAL,EVEC,STPDIA,STPSYM,STPCAR,
     &                  COOR,PMDIAX,PMDIA,PMCAR,SCAL,XMXNUC,TIMSTP,RNU,
     &                  ERGSEC,SGNMOD,ANHRAD,ZERGRD,TRSDIA,NCART,
     &                  NINTER,NCRIND,INDHES(0),IMODE,IWKIND,KEEPSY,
     &                  STATPO,WFPRED,START,INDXOK,IPRWLK,IWKTYP,ISTTYP,
     &                  IMAGE,DONWTN,TRUSTD)
            IF (IWKTYP .EQ. 7) THEN
C
C      Comment: I cannot see why this should not always be done!
C
               TRUSTR = STPNRM
               AACCUM = AACCUM + STPNRM
            END IF
            IF (IWKTYP .EQ. 3) TACCUM = TACCUM + TIMSTP
C           CALL WLKPRD(ERGSEC,GRDCAR,HESCAR,STPSYM,NCART)
         END IF
C
C        *******************************
C        ***** Molecule input file *****
C        *******************************
C
         CALL WLKMOL(COOR)
C
C        *****************************************
C        ***** Save wave function prediction *****
C        *****************************************
C
         IF (WFPRED) CALL WLKSAV(ERGSEC,COOR,WORK,LWORK)
      END IF
C
C     *************************
C     ***** Print results *****
C     *************************
C
      IF (IWKTYP .NE. 6) THEN
         CALL WLKRES(STPDIA,STPCAR,GRDDIA,GRDNRM,EVAL,EVEC,CNDHES,
     &               STPSYM,SCAL,COOR,RATIO,RNU,TRUSTX,TRUSTR,ERGSEC,
     &               OVLPMA,OVLPWM,COFMOD,ANHRAD,STPSYX,ESTART,ERGMOL,
     &               ERGMLX,ERGSCX,TACCUM,TIMSTP,EKINTO,EKINPA,VREL,
     &               GRDCAR,PMCAR,PMCARX,PMDIA,PMDIAX,CMDIF,ANHARM,
     &               TRSDIA,EKIN0,STPNRM,IMODE,IWKIND,DOREPW,NCRTOT,
     &               N2CRT,NCRIND,NCART,INDHES,IPOINT,ISTTYP,ISTIND,
     &               IWKTYP,STOP,STATPO,OLDREA,DOHESS,PRJTRO,IMAGE,
     &               AACCUM,DONWTN,IPRWLK,REJECT,ERGDIA,
     &               WORK,LWORK)
C
C     *************************************
C     ***** Numerical differentiation *****
C     *************************************
C
      ELSE
         IF (V3CAL .AND. NMODIF) THEN
            CALL WLKV3N(GRDCAR,GRAD0,GRADF,EVEC,
     &                 FMATF,WORK,ERGMOL,ESTART,DISPLC,
     &                 LWORK,NCRTOT,NCART,NUMMOD,IPOINT,IPRWLK)
         ELSE IF (V3CAL) THEN
            CALL WLKV3(GRDCAR,HESCAR,GRAD0,HESS0,GRADF,HESSF,
     &                 FMATF,WORK,ERGMOL,ESTART,DISPLC,
     &                 LWORK,NCRTOT,N2CRT,NCART,IPOINT,
     &                 IPRWLK)
         ELSE IF (VIBAVE) THEN
            CALL WLKVIB(GRDCAR,HESCAR,GRAD0,HESS0,HESSF,EVAL,
     &                  DIPM0,DIPGF,SUSTO0,SUSTF,GFAC0,GFACF,QUAD0,
     &                  QUADF,QUADT,SIGMA0,SIGMAF,SIGMAT,CSTRA,SCTRA,
     &                  SRC0,SRCF,POLAR0,POLARF,POLDD,ROAAFD,ROAAD,
     &                  EFG0,EFGF,SPNTOT,SSJ0,SSJF,POLDA,POLDL,ROAGND,
     &                  ROAGLD,ROAGNU,ROAGLU,ERGMOL,ESTART,WORK,LWORK,
     &                  NCART,NCRTOT,NUMMOD,IPOINT,IPRWLK)
         ELSE
            CALL WLKNUM(GRDCAR,HESCAR,DIPG0,GRAD0,HESS0,GRADF,HESSF,
     &                  DIPGF,WORK,ERGMOL,ESTART,DIPM,DIPM0,DISPLC,
     &                  LWORK,DODIP,NXYZ,NCRTOT,N2CRT,NCART,IPOINT,
     &                  NUMHES,IPRWLK)
            IF (VROA .OR. RAMAN) THEN
            CALL WLKROA(ROAAFD,ROAGND,ROAGLD,ROAAD,ROAAFU,ROAGNU,ROAGLU,
     &                  ROAAU,POLDD,POLDQ,POLDL,POLDA,WORK,LWORK,NXYZ,
     &                  NCART,IPOINT,DISPLC,IPRWLK)
            END IF
         END IF
         IF (STOP) THEN
            WRITE (LUPRI,'(/A)')' Numerical differentiation complete.'
         END IF
      END IF
C
C     ********************************
C     ***** Write interface file *****
C     ********************************
C
      CALL WLKWRT(EVAL,EVEC,PMCAR,GRAD0,HESS0,DIPM0,DIPG0,COOR0,GRADF,
     &            HESSF,DIPGF,ERGMOL,GRDNRM,ERGSEC,TRUSTR,ESTART,
     &            TACCUM,TIMSTP,EKINTO,EKINPA,COOR,GRDCAR,HESCAR,STPCAR,
     &            DISPLC,STPSYM,STPNRM,STPDIA,GRDDIA,EKIN0,TRUMAX,
     &            NCRTOT,N2CRT,NCRIND,NCART,IWKTYP,IWKIND,IMODE,ISCTYP,
     &            ISTIND,ISTTYP,NXYZ,DODIP,IMAGE,IPOINT,AACCUM,
     &            ROAAFD,ROAGND,ROAGLD,ROAAD,ROAAFU,ROAGNU,ROAGLU,ROAAU,
     &            FMATF,VIBAVE,SUSTO0,SUSTF,GFAC0,GFACF,QUAD0,QUADF,
     &            SIGMA0,SIGMAF,SRC0,SRCF,POLAR0,POLARF,
     &            EFG0,EFGF,SSJ0,SSJF,NUMMOD,WORK,LWORK)
      IF ((IPOINT .EQ. 2*NUMMOD) .AND. (V3CAL .AND. NMODIF)) THEN
         KFMATF = 1
         KFMATT = KFMATF + NCART*NCART*NCART
         KCORR  = KFMATT + NCART*NCART*NCART
         KPOS   = KCORR  + NUMMOD
         KFREQA = KPOS   + NCART
         KAMASS = KFREQA + NCART
         KLAST  = KAMASS + NUCDEP
         LLEFT  = LWORK  - KLAST
         DO I = 1, NUMMOD
            WORK(KFREQA + I - 1) = SQRT(ABS(EVAL(I)))
         END DO
         CALL DCOPY(NCART,COOR,1,CORD,1)
C
C     WORK(KFMATF) and WORK(KFMATT) is used as temporary arrays here.
C     The most important point is that we would like to initialize 
C     WORK(KAMASS)
C
         CALL CMMASS(WORK(KFMATF),WORK(KAMASS),WORK(KFMATT),WORK(KLAST),
     &               IPRINT)
         CALL VIBV3(EVEC,WORK(KFREQA),WORK(KFMATF),
     &              WORK(KFMATT),WORK(KCORR),WORK(KPOS),WORK(KAMASS),
     &              WORK(KLAST),LLEFT,NCART,NUMMOD)
      END IF
      RETURN
      END
C  /* Deck wlkdia */
      SUBROUTINE WLKDIA(GRDCAR,HESCAR,EVAL,EVEC,GRDDIA,TMAT,SCAL,
     &                  WORK,THRLDP,THRIND,CNDHES,INDHES,LWORK,IZEROG,
     &                  NZEROG,DOREPW,NCRTOT,N2CRT,NTMAT,ISCTYP,IWKTYP,
     &                  PRJTRO,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
      LOGICAL PRJTRO, DOREPW(0:7)
C
      DIMENSION GRDCAR(NCRTOT), HESCAR(N2CRT),
     *          EVAL(NCRTOT), EVEC(N2CRT), TMAT(NTMAT),
     *          SCAL(NCRTOT), WORK(LWORK), CNDHES(0:7),
     *          INDHES(0:7),  GRDDIA(NCRTOT), IZEROG(NZEROG)
C
#include <nuclei.h>
C
      CALL QENTER('WLKDIA')
C
C     ************************
C     ***** Mass scaling *****
C     ************************
C
      CALL WLKSCL(SCAL,GRDCAR,HESCAR,DOREPW,NCRTOT,N2CRT,ISCTYP,
     *            IWKTYP,IPRINT)
C
C     ************************************************
C     ***** Project out translation and rotation *****
C     ************************************************
C
      IF (PRJTRO) THEN
         CALL WLKPRJ(GRDCAR,HESCAR,SCAL,TMAT,WORK,THRLDP,
     *               LWORK,DOREPW,NCRTOT,N2CRT,NTMAT,IPRINT)
      END IF
C
C     ***********************************
C     ***** Diagonal representation *****
C     ***********************************
C
      CALL WLKEIG(GRDCAR,HESCAR,EVAL,EVEC,GRDDIA,TMAT,THRLDP,THRIND,
     &            WORK,CNDHES,INDHES,LWORK,IZEROG,NZEROG,DOREPW,NCRTOT,
     &            N2CRT,NTMAT,PRJTRO,IPRINT)
      CALL QEXIT('WLKDIA')
      RETURN
      END
C  /* Deck wlkcgh */
      SUBROUTINE WLKCGH(GRDCAR,HESCAR,DIPM,DIPG0,WORK,LWORK,DOREPW,
     *                  DOHESS,DODIP,NCRTOT,N2CRT,NXYZ,NCART,START,
     *                  IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (D1 = 1.0D0)
C
      LOGICAL DOREPW(0:7), DOHESS, DODIP, START
      DIMENSION GRDCAR(NCRTOT), HESCAR(N2CRT), DIPG0(3,NCART),
     *          DIPM(3), WORK(LWORK)
C
#include <taymol.h>
#include <moldip.h>
#include <nuclei.h>
#include <symmet.h>
#include <trkoor.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from WLKCGH','*',103)
C
C     ***** Initialize *****
C
      CALL DZERO(GRDCAR,NCRTOT)
C
C     ***** Multiplicities *****
C
      IF (NCRTOT .GT. LWORK) CALL STOPIT('WLKCGH',' ',LWORK,NCRTOT)
      DO 100 IREP = 0, MAXREP
         DO 110 ICENT = 1, NUCIND
#include <scalar.h>
            DO 120 ICOOR = 1, 3
               ISCOOR = IPTCNT(3*(ICENT-1)+ICOOR,IREP,1)
               IF (ISCOOR .GT. 0) THEN
                  WORK(ISCOOR) = D1/SQRT(FMULT(ISTBNU(ICENT)))
               END IF
 120        CONTINUE
 110     CONTINUE
 100  CONTINUE
C
C     ***** Normalized gradient *****
C
      DO 200 I = 1, NCRREP(0,1)
         GRDCAR(I) = WORK(I)*GRDMOL(I)
  200 CONTINUE
C
C     ***** Normalized Hessian *****
C
      IF (DOHESS) THEN
         CALL DZERO(HESCAR,N2CRT)
         IJ  = 1
         ICR = 0
         DO 300 IREP = 0, MAXREP
            IF (DOREPW(IREP)) THEN
               NCR = NCRREP(IREP,1)
               DO 310 I = ICR + 1, ICR + NCR
                  DO 320 J = ICR + 1, ICR + NCR
                     HESCAR(IJ)=WORK(I)*HESMOL(I,J)*WORK(J)
                     IJ = IJ + 1
  320             CONTINUE
  310          CONTINUE
            END IF
            ICR = ICR + NCR
  300    CONTINUE
      END IF
C
C     ***** Normalized dipole gradient *****
C
      IF (DODIP) THEN
         DIPM(IPTAX(1,1)) = DIP0(1)
         DIPM(IPTAX(2,1)) = DIP0(2)
         DIPM(IPTAX(3,1)) = DIP0(3)
         IF (START) THEN
            CALL DZERO(DIPG0,3*NCART)
            DO 400 I = 1, NXYZ
               DO 410 J = 1, NCRREP(0,1)
                  DIPG0(I,J) = DIPFLT(I,J)*WORK(J)
  410          CONTINUE
  400       CONTINUE
         END IF
      END IF
C
C     ***** Print *****
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Normalization factors in WLKCGH',-1)
         CALL OUTPUT(WORK,1,1,1,NCRTOT,1,NCRTOT,1,LUPRI)
         CALL HEADER('Gradient in WLKCGH',-1)
         CALL OUTPUT(GRDCAR,1,1,1,NCRTOT,1,NCRTOT,1,LUPRI)
         IF (DOHESS) THEN
            CALL HEADER('Hessian in WLKCGH',-1)
            IHESS = 1
            DO 500 IREP = 0, MAXREP
               WRITE (LUPRI, '(/,1X,A,I2,/)') ' Symmetry ', IREP + 1
               IF (DOREPW(IREP)) THEN
                  NCR = NCRREP(IREP,1)
                  CALL OUTPUT(HESCAR(IHESS),1,NCR,1,NCR,NCR,NCR,1,LUPRI)
                  IHESS = IHESS + NCR*NCR
               ELSE
                  WRITE (LUPRI, '(1X,A)')
     *                  ' No internal coordinates of this symmetry.'
               END IF
  500       CONTINUE
         END IF
         IF (DODIP) THEN
            CALL HEADER('Dipole gradient in WLKCGH',-1)
            CALL OUTPUT(DIPG0,1,NXYZ,1,NCART,3,NCART,1,LUPRI)
         END IF
      END IF
      RETURN
      END
C  /* Deck wlkscl */
      SUBROUTINE WLKSCL(SCAL,GRDCAR,HESCAR,DOREP,NCRTOT,N2CRT,ISCTYP,
     &                  IWKTYP,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <codata.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (D1 = 1.0D0)
C
      LOGICAL DOREP(0:7), NOCHRG
      DIMENSION SCAL(NCRTOT), GRDCAR(NCRTOT), HESCAR(N2CRT)
C
#include <nuclei.h>
#include <symmet.h>
#include <trkoor.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from WLKSCL','*',103)
C
C     ***** Check for zero charges *****
C
      IF (ISCTYP .EQ. 3) THEN
         NOCHRG = .FALSE.
         DO 100 INUC = 1, NUCIND
            NOCHRG = NOCHRG .OR. IZATOM(INUC).LE.0
 100     CONTINUE
         IF (NOCHRG) THEN
            WRITE (LUPRI, '(/A/A/A)')
     *      ' WARNING: Scaling reset from mass scaling to no scaling',
     *      '          because at least one center has no charge or is',
     *      '          a point charge, or multiple basis sets are used'
            ISCTYP = 1
         END IF
      END IF
C
C     ***** No scaling *****
C
      IF (ISCTYP .EQ. 1) THEN
         DO 200 I = 1, NCRTOT
            SCAL(I) = D1
 200     CONTINUE
C
C     ***** Scaling factors from input *****
C
      ELSE IF (ISCTYP .EQ. 2) THEN
         CALL DCOPY(NCRTOT,SCALCO,1,SCAL,1)
C
C     ***** Mass scaling *****
C
      ELSE
         NOFF = 0
         DO 300 IREP = 0, MAXREP
            IF (DOREP(IREP)) THEN
               DO 310 IATOM = 1, NUCIND
                  NATTYP = IZATOM(IATOM)
                  DMASS  = DISOTP(NATTYP,ISOTOP(IATOM),'MASS')
                  DO 320 ICOOR = 1, 3
                     JCOOR = 3*(IATOM - 1) + ICOOR
                     ISCOOR = IPTCNT(JCOOR,IREP,1)
                     IF (ISCOOR .NE. 0) THEN
                        NOFF = NOFF + 1
                        SCAL(NOFF) = D1/SQRT(DMASS)
                     END IF
 320              CONTINUE
 310           CONTINUE
            END IF
 300     CONTINUE
         IF (IWKTYP .EQ. 3) CALL DSCAL(NCRTOT,D1/SQRT(XFAMU),SCAL,1)
      END IF
C
C     ***** Exit if no scaling *****
C
      IF (ISCTYP .EQ. 1) GO TO 9999
C
C     ***** Print *****
C
      IF (IPRINT .GT. 2 .AND. ISCTYP .GT. 1) THEN
         CALL HEADER('Coordinate scale factors in this calculation:',5)
         CALL OUTPUT(SCAL,1,1,1,NCRTOT,1,NCRTOT,1,LUPRI)
      END IF
C
C     ***** Scale gradient *****
C
      DO 400 I = 1, NCRREP(0,1)
         GRDCAR(I) = SCAL(I)*GRDCAR(I)
  400 CONTINUE
C
C     ***** Scale Hessian *****
C
      IJ   = 1
      ICR = 0
      DO 500 IREP = 0, MAXREP
         IF (DOREP(IREP)) THEN
            NCR = NCRREP(IREP,1)
            DO 510 I = 1, NCR
               DO 520 J = 1, NCR
                  HESCAR(IJ) = SCAL(ICR+I)*SCAL(ICR+J)*HESCAR(IJ)
                  IJ = IJ + 1
  520          CONTINUE
  510       CONTINUE
            ICR = ICR + NCR
         END IF
  500 CONTINUE
C
C     ***** Print *****
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Gradient in WLKSCL',-1)
         CALL OUTPUT(GRDCAR,1,1,1,NCRTOT,1,NCRTOT,1,LUPRI)
         CALL HEADER('Hessian in WLKSCL',-1)
         IHESS = 1
         DO 600 IREP = 0, MAXREP
            WRITE (LUPRI, '(/,1X,A,I2)') ' Symmetry ', IREP + 1
            IF (DOREP(IREP)) THEN
               NCR = NCRREP(IREP,1)
               CALL OUTPUT(HESCAR(IHESS),1,NCR,1,NCR,NCR,NCR,1,LUPRI)
               IHESS = IHESS + NCR*NCR
            ELSE
               WRITE (LUPRI, '(1X,A/)')
     *               ' No internal coordinates of this symmetry.'
            END IF
  600    CONTINUE
      END IF
 9999 CONTINUE
      RETURN
      END
C  /* Deck wlkprj */
      SUBROUTINE WLKPRJ(GRDCAR,HESCAR,SCAL,TMAT,WORK,THRLDP,
     *                  LWORK,DOREPW,NCRTOT,N2CRT,NTMAT,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <iratdef.h>
C
      LOGICAL DOREPW(0:7)
      DIMENSION GRDCAR(NCRTOT), HESCAR(N2CRT), SCAL(NCRTOT),
     *          TMAT(NTMAT), WORK(LWORK)
C
#include <symmet.h>
#include <trkoor.h>
#include <prkoor.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from WLKPRJ','*',103)
C
      IHESS = 1
      ISCAL = 1
      ITMAT = 1
      DO 100 IREP = 0, MAXREP
         IF (DOREPW(IREP)) THEN
            NCR = NCRREP(IREP,1)
            NPR = NPRREP(IREP)
            NTR = NTRREP(IREP)
            IF (IPRINT .GT. 5) THEN
               WRITE (LUPRI,'(1X,A,I2)') ' Symmetry ', IREP + 1
               WRITE (LUPRI,'(1X,A,I2)') ' NCR      ', NCR
               WRITE (LUPRI,'(1X,A,I2)') ' NTR      ', NTR
               WRITE (LUPRI,'(1X,A,I2)') ' NPR      ', NPR
            END IF
            KPRJCT = 1
            KPRJHS = KPRJCT + NCR*NCR
            KPTCOL = KPRJHS + NCR*NCR
            KLAST  = KPTCOL + (NCOOR + IRAT - 1)/IRAT
            IF (KLAST .GT. LWORK) CALL STOPIT('WLKPRJ',' ',KLAST,LWORK)
            CALL WLKPR1(GRDCAR,HESCAR(IHESS),SCAL(ISCAL),TMAT(ITMAT),
     *                  WORK(KPRJCT),WORK(KPRJHS),THRLDP,WORK(KPTCOL),
     *                  NCR,NPR,NTR,NCOOR,IREP,IPRINT)
            IHESS = IHESS + NCR*NCR
            ISCAL = ISCAL + NCR
            ITMAT = ITMAT + NPR*NCR
         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck wlkpr1 */
      SUBROUTINE WLKPR1(GRDCAR,HESCAR,SCAL,TMAT,PROJCT,PRJHES,THRLDP,
     *                  IPTCOL,NCR,NPR,NTR,NCOOR,IREP,IPRINT)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1 = 1.0D0)
C
      DIMENSION TMAT(NCR,NPR), PROJCT(NCR,NCR), HESCAR(NCR,NCR),
     *          SCAL(NCR), PRJHES(NCR,NCR), IPTCOL(NCOOR), GRDCAR(NCR),
     *          LINDEP(9)
#include <cbisol.h>
C
      IF (IPRINT .GT. 5) THEN
         CALL TITLER('Output from WLKPR1','*',103)
         IF (IREP .GT. 0) THEN
            CALL HEADER('Unprojected gradient in WLKPR1',-1)
            CALL OUTPUT(GRDCAR,1,1,1,NCR,1,NCR,1,LUPRI)
         END IF
         CALL HEADER('Unprojected Hessian in WLKPR1',-1)
         CALL OUTPUT(HESCAR,1,NCR,1,NCR,NCR,NCR,1,LUPRI)
      END IF
C
      CALL DZERO(TMAT,NCR*NPR)
C
C     ***** Construct translation and rotation matrix *****
C
      CALL GETTRO(TMAT,IPTCOL,NCR,NTR,'BOTH','TORTHO','CT',IREP,IPRINT)
C
C     ***** Vector representing translation of cavity *****
C
      IF (SOLVNT) CALL WLKCVT(TMAT,IPTCOL,NCR,NPR,NTR,IREP,IPRINT)
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('T matrix in WLKPR1',-1)
         CALL OUTPUT(TMAT,1,NCR,1,NPR,NCR,NPR,1,LUPRI)
      END IF
C
C     ***** Scale T matrix *****
C
      DO 100 ICR = 1, NCR
         FAC = D1/(SCAL(ICR)*SCAL(ICR))
         CALL DSCAL(NPR,FAC,TMAT(ICR,1),NCR)
  100 CONTINUE
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Scaled T matrix in WLKPR1',-1)
         CALL OUTPUT(TMAT,1,NCR,1,NPR,NCR,NPR,1,LUPRI)
      END IF
C
C     ***** Orthogonalize translations and rotations *****
C
      NPR1 = NPR
      CALL ORTVEC(0,NPR1,NCR,THRLDP,TMAT,LINDEP)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Orthogonalized T matrix in WLKPR1',-1)
         CALL OUTPUT(TMAT,1,NCR,1,NPR,NCR,NPR,1,LUPRI)
      END IF
      IF (NPR1 .NE. NPR) THEN
         WRITE (LUPRI,'(//,2(A,I1),A,/A)')
     *      ' Number of trarot vectors reduced in ',NPR,
     *      ' to ',NPR1,' in ORTVEC called from WLKPRJ.',
     *      ' Program cannot proceed .'
         CALL QUIT('Insufficient number of trarot vectors in WLKPRJ.')
      END IF
C
C     ***** Construct projection operator *****
C
      CALL DUNIT(PROJCT,NCR)
      CALL DGEMM('N','T',NCR,NCR,NPR,-1.D0,
     &           TMAT,NCR,
     &           TMAT,NCR,1.D0,
     &           PROJCT,NCR)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Projection operator in WLKPR1',-1)
         CALL OUTPUT(PROJCT,1,NCR,1,NCR,NCR,NCR,1,LUPRI)
      END IF
C
C     ***** Project out translations and rotations *****
C
      IF (IREP .EQ. 0) THEN
         CALL DGEMM('N','N',NCR,1,NCR,1.D0,
     &              PROJCT,NCR,
     &              GRDCAR,NCR,0.D0,
     &              PRJHES,NCR)
         CALL DCOPY(NCR,PRJHES,1,GRDCAR,1)
         IF (IPRINT .GT. 5) THEN
            CALL HEADER('Projected gradient in WLKPR1',-1)
            CALL OUTPUT(GRDCAR,1,1,1,NCR,1,NCR,1,LUPRI)
         END IF
      END IF
      CALL DGEMM('N','N',NCR,NCR,NCR,1.D0,
     &           PROJCT,NCR,
     &           HESCAR,NCR,0.D0,
     &           PRJHES,NCR)
      IF (IPRINT .GT. 20) THEN
         CALL HEADER('Left-projected Hessian in WLKPR1',-1)
         CALL OUTPUT(PRJHES,1,NCR,1,NCR,NCR,NCR,1,LUPRI)
      END IF
      CALL DGEMM('N','N',NCR,NCR,NCR,1.D0,
     &           PRJHES,NCR,
     &           PROJCT,NCR,0.D0,
     &           HESCAR,NCR)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Projected Hessian in WLKPR1',-1)
         CALL OUTPUT(HESCAR,1,NCR,1,NCR,NCR,NCR,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck wlkeig */
      SUBROUTINE WLKEIG(GRDCAR,HESCAR,EVAL,EVEC,GRDDIA,TMAT,THRLDP,
     &                  THRIND,WORK,CNDHES,INDHES,LWORK,IZEROG,NZEROG,
     &                  DOREPW,NCRTOT,N2CRT,NTMAT,PRJTRO,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <iratdef.h>
      LOGICAL PRJTRO, DOREPW(0:7)
      DIMENSION GRDCAR(NCRTOT), HESCAR(N2CRT), TMAT(NTMAT),
     *          EVEC(N2CRT), EVAL(NCRTOT),  WORK(LWORK),
     *          CNDHES(0:7), INDHES(0:7), GRDDIA(NCRTOT),
     *          IZEROG(NZEROG)
C
#include <symmet.h>
#include <trkoor.h>
#include <prkoor.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from WLKEIG','*',103)
C
      IHESS  = 1
      IEVAL  = 1
      ITMAT  = 1
      CALL DZERO(GRDDIA,NCRTOT)
      CALL IZERO(INDHES(0),MAXREP+1)
      DO 100 IREP = 0, MAXREP
         IF (DOREPW(IREP)) THEN
            NCR = NCRREP(IREP,1)
            NPR = NPRREP(IREP)
            IF (IPRINT .GT. 5) THEN
               WRITE (LUPRI,'(1X,A,I2)') ' Symmetry ', IREP + 1
               WRITE (LUPRI,'(1X,A,I2)') ' NCR      ', NCR
               WRITE (LUPRI,'(1X,A,I2)') ' NPR      ', NPR
            END IF
            KPACK  = 1
            KOVRL  = KPACK  + NCR*(NCR + 1)/2
            KWORK  = KOVRL  + NCR
            KIWORK = KWORK  + NCR
            KTRVEC = KIWORK + (NCR + IRAT - 1)/IRAT
            KLAST  = KTRVEC + (NCR + IRAT - 1)/IRAT
            IF (KLAST .GT. LWORK) CALL STOPIT('WLKEIG',' ',KLAST,LWORK)
            CALL WLKEI1(GRDCAR,HESCAR(IHESS),GRDDIA(IEVAL),EVAL(IEVAL),
     &                  EVEC(IHESS),WORK(KOVRL),TMAT(ITMAT),WORK(KPACK),
     &                  WORK(KWORK),THRLDP,THRIND,WORK(KIWORK),
     &                  WORK(KTRVEC),IZEROG,NZEROG,NCR,NPR,INDEX,INDEXM,
     &                  IREP,PRJTRO,IPRINT)
            INDHES(IREP) = INDEX
            IF (PRJTRO) THEN
               CALL WLKCND(CNDNMB,EVAL(IEVAL),NCR-NPR)
               CNDHES(IREP) = CNDNMB
            END IF
            IHESS  = IHESS + NCR*NCR
            IEVAL  = IEVAL + NCR
            ITMAT  = ITMAT + NCR*NPR
         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck wlkei1 */
      SUBROUTINE WLKEI1(GRDCAR,HESCAR,GRDDIA,EVAL,EVEC,OVERLP,TMAT,
     &                  HESPCK,WORK,THRLDP,THRIND,IWORK,TRVEC,IZEROG,
     &                  NZEROG,NCR,NPR,INDEX,INDEXM,IREP,PRJTRO,
     &                  IPRINT)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0, THRZER = 1.0D-04)
C
      LOGICAL PRJTRO, TRVEC(NCR)
      DIMENSION GRDCAR(NCR), HESCAR(NCR,NCR), GRDDIA(NCR), EVAL(NCR),
     *          EVEC(NCR,NCR), HESPCK(NCR*(NCR + 1)/2), TMAT(NCR,NPR),
     *          WORK(NCR), IWORK(NCR), IZEROG(NZEROG), OVERLP(NCR)
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from WLKEI1','*',103)
C
C     ***** Diagonalize Hessian *****
C
      CALL DSITSP(NCR,HESCAR,HESPCK)
      CALL DUNIT(EVEC,NCR)
      CALL JACO(HESPCK,EVEC,NCR,NCR,NCR,WORK,IWORK)
      DO 100 I=1, NCR
         EVAL(I) = HESPCK(I*(I+1)/2)
 100  CONTINUE
      IF (IREP .EQ. 0) THEN
         DO 200 I=1,NCR
            GRDDIA(I) = DDOT(NCR,GRDCAR,1,EVEC(1,I),1)
 200     CONTINUE
      END IF
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Eigenvalues in WLKEI1',-1)
         CALL OUTPUT(EVAL,1,1,1,NCR,1,NCR,1,LUPRI)
         CALL HEADER('Eigenvectors in WLKEI1',-1)
         CALL OUTPUT(EVEC,1,NCR,1,NCR,NCR,NCR,1,LUPRI)
         CALL HEADER('Gradient (diagonal rep.) in WLKEI1',-1)
         CALL OUTPUT(GRDDIA,1,1,1,NCR,1,NCR,1,LUPRI)
      END IF
C
C     ***** Sort in order of increasing trarot component *****
C
      IF (PRJTRO) THEN
         DO 300 ICR = 1, NCR
            OVLP = D0
            DO 310 ITR = 1, NPR
               OVLP = OVLP + ABS(DDOT(NCR,EVEC(1,ICR),1,TMAT(1,ITR),1))
  310       CONTINUE
            OVERLP(ICR) = OVLP
  300    CONTINUE
         CALL DCOPY(NCR,OVERLP,1,WORK,1)
         CALL ORDER(GRDDIA,WORK,NCR,1)
         CALL DCOPY(NCR,OVERLP,1,WORK,1)
         CALL ORDER(EVAL,  WORK,NCR,1)
         CALL DCOPY(NCR,OVERLP,1,WORK,1)
         CALL ORDER(EVEC,  WORK,NCR,NCR)
      END IF
C
      IF (PRJTRO) THEN
         NVEC = NCR - NPR
      ELSE
         NVEC = NCR
      END IF
C
C     Number of negative eigenvalues
C
      INDEX = 0
      DO 700 I = 1, NVEC
         IF (EVAL(I) .LT. -THRIND) INDEX = INDEX + 1
  700 CONTINUE
      IF (IPRINT.GT.5) THEN
         WRITE (LUPRI,'(/A,I2/)')    ' Index of Hessian:      ',INDEX
      END IF
C
C     Order non-trarot eigenvalues in increasing order
C
      DO 500 I = 1, NVEC
         JMIN = I
         EMIN = EVAL(I)
         DO 510 J = (I + 1), NVEC
            IF (EVAL(J) .LT. EMIN) THEN
               EMIN = EVAL(J)
               JMIN = J
            END IF
  510    CONTINUE
         IF (JMIN .NE. I) THEN
            CALL DSWAP(1,  EVAL  (I),1,EVAL  (JMIN),1)
            CALL DSWAP(NCR,EVEC(1,I),1,EVEC(1,JMIN),1)
            CALL DSWAP(1,GRDDIA(I),1,GRDDIA(JMIN),1)
         END IF
  500 CONTINUE
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Sorted eigenvalues in WLKEI1',-1)
         CALL OUTPUT(EVAL,1,1,1,NCR,1,NCR,1,LUPRI)
         CALL HEADER('Sorted eigenvectors in WLKEI1',-1)
         CALL OUTPUT(EVEC,1,NCR,1,NCR,NCR,NCR,1,LUPRI)
         CALL HEADER('Gradient (sorted) in WLKEI1',-1)
         CALL OUTPUT(GRDDIA,1,1,1,NCR,1,NCR,1,LUPRI)
      END IF
C
C     Set gradient elements zero as requested
C
      DO 600 I =1, NZEROG
         GRDDIA(IZEROG(I)) = D0
  600 CONTINUE
      RETURN
      END
C  /* Deck wlkcnd */
      SUBROUTINE WLKCND(CNDNMB,EVAL,NONTRO)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D1000 = 1.0D3, DEXP10 = 1.0D10)
      DIMENSION EVAL(NONTRO)
C
C     Condition number of Hessian
C
      EMAX = D0
      EMIN = D1000
      DO 100 I = 1, NONTRO
         EVALI = ABS(EVAL(I))
         EMAX  = MAX(EMAX,EVALI)
         EMIN  = MIN(EMIN,EVALI)
  100 CONTINUE
      IF (EMIN .LT. D1/DEXP10) THEN
         WRITE (LUPRI,'(//A,1P,D12.5)')
     *     ' Hessian condition number not calculated since smallest'
     *     //' eigenvalue is ',EMIN
         CNDNMB = DEXP10
      ELSE
         CNDNMB = EMAX/EMIN
      END IF
      RETURN
      END
C  /* Deck wlkrad */
      SUBROUTINE WLKRAD(ERGMOL,ERGMLX,ERGSCX,STPNRX,OVLPWM,COFMOD,
     *                  RATIO,ISTTPX,IPRINT)
C
C     Update trust radius
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (D1=1.0D0, D2=2.0D0, DP5 = 0.5D0, DP25=0.25D0,
     *           THDE=1.0D-10)
C
#include <infinp.h>
#include <cbiwlk.h>
C
C     ****************************************
C     ***** Rejection specified in input *****
C     ****************************************
C
      IF (REJECT) THEN
         TRUSTR = TRUSTD*STPNRX
         IF (IPRINT .GT. 2) WRITE (LUPRI,'(/A/A,F12.6/)')
     &      ' Trust radius has been reduced because step was rejected.',
     &      ' Updated trust radius :', TRUSTR
C
C     *********************************************
C     ***** Calculate ratio, reject or update *****
C     *********************************************
C
      ELSE
C
C        1) Caculate ratio between actual and predicted energy change
C        ------------------------------------------------------------
C
         DEACT = ERGMOL - ERGMLX
         IF (ABS(ERGSCX) .GT. THDE) THEN
            RATIO = DEACT/ERGSCX
            IF (IPRINT .GT. 2) THEN
               CALL HEADER('Energy difference to previous geometry:',1)
               WRITE (LUPRI,'(/A,/)')
     *         '     Actual           /  Predicted        =  Ratio '
               WRITE (LUPRI,'(5X,F15.10,2(A,F15.10))')
     *                DEACT,'  /  ',ERGSCX,'  =  ',RATIO
            END IF
         ELSE
            RATIO = D1
            IF (IPRINT .GT. 2) THEN
               WRITE (LUPRI,'(3(/A),/5X,1P,2D15.5)')
     *             ' Close to convergence, ratio set to one.',
     *             ' Energy difference to previous geometry;',
     *             ' actual and predicted:', DEACT,ERGSCX
            END IF
         ENDIF
C
C        2) Update trust radius based on ratio
C        -------------------------------------
C
C        Case 1: Minimization walk
C
         IF (IWKTYP.EQ.1 .AND. IWKIND.EQ.0 .AND. .NOT.IMAGE) THEN
            IF (RATIO .LT. REJMIN .AND. ISTATE .EQ. 1) THEN
               WRITE (LUPRI,'(/A/A,F12.6)')
     *           ' Walk to minimum stopped because ratio between',
     *           ' predicted and actual energy change is:',RATIO
               CALL QUIT
     *            ('*** WLKCTL *** Ratio outside desired interval.')
            ELSE IF (ABS(RATIO) .LT. RTRMIN) THEN
               TRUSTR = TRUSTD*STPNRX
            ELSE IF (ABS(RATIO) .GT. RTRGOD) THEN
               TRUSTR = MAX(TRUSTI*STPNRX,TRUSTR)
            ELSE
               TRUSTR = STPNRX
            END IF
C
C           if (ratio very good) square increase of trust radius
C
            IF (ABS(RATIO-D1) .LT. DP25*(D1-RTRGOD)) THEN
               TRUSTR = MAX(TRUSTI*TRUSTI*STPNRX,TRUSTR)
            END IF
C
C        Case 2: Transition state walk or dynamic path
C
         ELSE
            IF (IMAGE) THEN
               RED = D1
            ELSE
               RED1 = D1
               RED2 = D1
               IF (ISTTPX .NE. 4) THEN
                  IF((IWKTYP.EQ.1.AND.IWKIND.GT.0).OR.(IWKTYP.EQ.2))THEN
                     IF (ABS(OVLPWM) .LE. 0.95D0) RED1 = TRUSTD*RED1
                     IF (ABS(OVLPWM) .LE. 0.90D0) RED1 = TRUSTD*RED1
                     IF (ABS(OVLPWM) .LE. 0.85D0) RED1 = TRUSTD*RED1
                  END IF
                  IF (IWKTYP .EQ. 2) THEN
                     IF (ABS(COFMOD) .LE. 0.90D0) RED2 = TRUSTD*RED2
                  END IF
               END IF
               RED = MIN(RED1,RED2)
            END IF
C
            IF ((RATIO.LT.REJMIN) .OR. (RATIO.GT.REJMAX)) THEN
               WRITE (LUPRI,'(/2A/A,F12.6)')
     *            ' Walk to transition state stopped because',
     *            ' ratio between',
     *            ' predicted and actual energy change is:',RATIO
               CALL QUIT
     *            ('*** WLKCTL *** Ratio outside desired interval.')
            ELSE IF ((RATIO.LT.RTRMIN) .OR. (RATIO.GT.(D2-RTRMIN))) THEN
               TRUSTR = MIN(RED,TRUSTD)*STPNRX
            ELSE IF ((RATIO.GT.RTRGOD) .AND.(RATIO.LT.(D2-RTRGOD))) THEN
               IF (RED .LT. D1) THEN
                  TRUSTR = RED*TRUSTI*STPNRX
               ELSE
                  TRUSTR = MAX(TRUSTI*STPNRX,TRUSTR)
               END IF
            ELSE
               TRUSTR = RED*STPNRX
            END IF
         ENDIF
         IF (IWKTYP .EQ. 7) TRUSTR = MIN(TRUSTR,TRUMAX)
         IF (IPRINT .GT. 5) THEN
            WRITE (LUPRI,'(/A,2F12.6)') ' STPNRX     ', STPNRX
            WRITE (LUPRI,'(/A,2F12.6)') ' RED1, RED2 ', RED1, RED2
            WRITE (LUPRI,'(/A,F12.6)') ' Final TRUSTR ', TRUSTR
         END IF
      END IF
      RETURN
      END
C  /* Deck wlkmod */
      SUBROUTINE WLKMOD(EVAL,EVALX,EVEC,EVECX,GRDDIA,OVLPMA,XMAXOV,
     &                  COFMOD,ANHRAD,TRUSTX,NCART,NINTER,IMODE,IMAGE,
     &                  IPRINT)
C
C     Determine current mode by comparing eigenvectors with those
C     of last iteration
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D3 = 3.0D0, D100 = 100.0D0,
     *           XMINOV = 0.7D0)
C
      LOGICAL IMAGE
      DIMENSION EVEC(NCART,NCART), EVECX(NCART,NCART),
     *          OVLPMA(NCART,NCART), EVAL(NCART), EVALX(NCART),
     *          GRDDIA(NCART)
C
C     **********************************************************
C     ***** Find largest overlap and identify current mode *****
C     **********************************************************
C
      IF (IMAGE .AND. (IMODE.EQ.1)) THEN
         NMODE  = IMODE
         XMAXOV = ABS(DDOT(NCART,EVEC(1,1),1,EVECX(1,1),1))
      ELSE
         IOLD = IMODE
         XMAXOV = D0
         DO 100 K = 1, NINTER
            X = ABS(DDOT(NCART,EVEC(1,K),1,EVECX(1,IOLD),1))
            IF (X.GT.XMAXOV) THEN
               XMAXOV = X
               NMODE  = K
            END IF
  100    CONTINUE
      END IF
C
C     ***** Print *****
C
      IF (NMODE .NE. IMODE) THEN
         WRITE (LUPRI,'(/A/,2(3X,A,I3/))')
     *      ' Crossing of Hessian eigenvalues has occured:',
     *      ' New walk mode number is ', NMODE,
     *      ' Old walk mode number was', IOLD
      END IF
      WRITE (LUPRI,'(/A,I3,A/A,F8.5)')
     *   ' Overlap between normal coordinate',NMODE,
     *   ' and the normal coordinate walk ',
     *   ' direction of the previous iteration is equal to',
     *     XMAXOV
      IF (XMAXOV .LT. XMINOV) THEN
         WRITE (LUPRI,'(/A)')
     *     ' WARNING: This was largest overlap found.'
      END IF
C
      IMODE = NMODE
C
C     **********************************
C     ***** Calculate all overlaps *****
C     **********************************
C
      IF (IPRINT .GE. 3) THEN
         CALL DGEMM('T','N',NCART,NCART,NCART,1.D0,
     &              EVEC,NCART,
     &              EVECX,NCART,0.D0,
     &              OVLPMA,NCART)
         CALL HEADER('Overlap between old and new eigenvectors:',-1)
         CALL OUTPUT(OVLPMA,1,NCART,1,NCART,NCART,NCART,1,LUPRI)
      END IF
C
      GRDNRM = DNRM2(NINTER,GRDDIA,1)
      IF (GRDNRM .GT. 1.0D-7) THEN
         COFMOD = GRDDIA(IMODE)/GRDNRM
      ELSE
         COFMOD = D1
      END IF
      IF (.NOT.IMAGE) THEN
         DELTA = EVAL(IMODE) - EVALX(IOLD)
         IF (ABS(DELTA) .GT. 1.0D-7) THEN
            ANHRAD = D3*TRUSTX*ABS(EVAL(IMODE)/DELTA)
         ELSE
            ANHRAD = D100
         END IF
C
         IF (IPRINT .GT. 5) THEN
            WRITE (LUPRI,'(/A,F20.10)') ' COFMOD ', COFMOD
            WRITE (LUPRI,'( A,F20.10)') ' ANHRAD ', ANHRAD
            WRITE (LUPRI,'( A,F20.10)') ' EVAL   ', EVAL(IMODE)
            WRITE (LUPRI,'( A,F20.10)') ' EVALX  ', EVALX(IMODE)
            WRITE (LUPRI,'( A,F20.10)') ' TRUSTX ', TRUSTX
         END IF
      END IF
      RETURN
      END
C  /* Deck wlkmol */
      SUBROUTINE WLKMOL(COONEW)
C
C 28-Mar-1985 hjaaj
C Modified 26-Nov-1986 tuh, for symmetry 17-Oct-1989
C
C     Write MOLECULE integral input deck on unit LUMOL.
C
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.0D0, D10 = 10.00D00, D100 = 100.00D00)
      PARAMETER (D1000 = 1000.00D0, D10000 = 10000.00D00)
      PARAMETER (D100000 = 100000.00D0, D1000000 = 1000000.00D00)
C
      CHARACTER NAME*4, TMPTXT*11
      DIMENSION COONEW(3,*)
C
C Used from common blocks:
C  /MOLINP/: everything
C  /NUCLEI/: NUCIND
C  /CBISOL/: SOLVNT,NCNTCV
C
#include <molinp.h>
#include <inftap.h>
#include <symmet.h>
#include <infinp.h>
#include <nuclei.h>
#include <cbisol.h>
#include <ibtfun.h>
C
C
C Write new coordinates into MOLECULE input
C
      IF (SOLVNT) THEN
C     ... make a translation of center of cavity to (0,0,0)
         CCAVX = COONEW(1,NCNTCV)
         CCAVY = COONEW(2,NCNTCV)
         CCAVZ = COONEW(3,NCNTCV)
         IF (CCAVX .NE. D0 .OR. CCAVY .NE. D0 .OR. CCAVZ .NE. D0) THEN
            WRITE (LUPRI,'(/A/A,3F20.10)')
     &         'WLKMOL INFO: Center of solvent cavity moved to (0,0,0)',
     &         '   from',CCAVX,CCAVY,CCAVZ
         END IF
      ELSE
         CCAVX = D0
         CCAVY = D0
         CCAVZ = D0
         IF (NFIELD .GT. 0) THEN
            TOTCHR = D0
C
C     If we have a finite field, we must ensure that there is a zero
C     contribution from the nuclear part of the dipole moment to the 
C     total charge. Thus move center-of nuclear charge to the origin
C
            DO IATOM = 1, NUCIND
               DO IREP = 0, MAXREP
                  IF (IBTAND(IREP,ISTBNU(IATOM)) .EQ. 0) THEN
                     CCAVX = CCAVX + PT(IBTAND(ISYMAX(1,1),IREP))
     &                              *COONEW(1,IATOM)*CHARGE(IATOM)
                     CCAVY = CCAVY + PT(IBTAND(ISYMAX(2,1),IREP))
     &                              *COONEW(2,IATOM)*CHARGE(IATOM)
                     CCAVZ = CCAVZ + PT(IBTAND(ISYMAX(3,1),IREP))
     &                              *COONEW(3,IATOM)*CHARGE(IATOM)
                     TOTCHR = TOTCHR + CHARGE(IATOM)
                  END IF
               END DO
            END DO
            CCAVX = CCAVX/TOTCHR
            CCAVY = CCAVY/TOTCHR
            CCAVZ = CCAVZ/TOTCHR
            IF(CCAVX .NE. D0 .OR. CCAVY .NE. D0 .OR. CCAVZ .NE. D0) THEN
               WRITE (LUPRI,'(/A/A,3F20.10)')
     &         'WLKMOL INFO: Center of nuclear charge moved to (0,0,0)',
     &         '   from',CCAVX,CCAVY,CCAVZ
            END IF
         END IF
      END IF
C
C     Coordinates written are always in atomic units
C
      CALL UPCASE(MLINE(NMLAU))
      IPOS = INDEX(MLINE(NMLAU),'ATO')
      IF (IPOS .EQ. 0) THEN
         MLINE(NMLAU)(20:20) = ' '
      ELSE
         IPOS = INDEX(MLINE(NMLAU),'ANG')
         IF (IPOS .NE. 0) THEN
            IPOS2 = INDEX(MLINE(NMLAU)(IPOS:),' ')
            DO ICHARP = IPOS, IPOS + IPOS2 - 1
               WRITE (MLINE(NMLAU)(ICHARP:ICHARP),'(A1)') ' '
            END DO
         END IF
      END IF
      DO 100 N = 1, NUCIND
         CRX = COONEW(1,N) - CCAVX
         CRY = COONEW(2,N) - CCAVY
         CRZ = COONEW(3,N) - CCAVZ
         NC = NCLINE(N)
         IF (NC .NE. 0) THEN
            READ  (MLINE(NC),9100) NAME
            IPOS = INDEX(MLINE(NC),'Isotope=')
            IF (IPOS .NE. 0) THEN
               READ (MLINE(NC)(IPOS:),'(A11)') TMPTXT
            ELSE
               TMPTXT = '           '
            END IF
            CRDLRG = MAX(MAX(ABS(CRX),ABS(CRY)),ABS(CRZ))
            IF (CRDLRG .LT. D10) THEN
               WRITE (MLINE(NC),9060) NAME,CRX,CRY,CRZ,TMPTXT
            ELSE IF (CRDLRG .LT. D100) THEN
               WRITE (MLINE(NC),9050) NAME,CRX,CRY,CRZ,TMPTXT
            ELSE IF (CRDLRG .LT. D1000) THEN
               WRITE (MLINE(NC),9040) NAME,CRX,CRY,CRZ,TMPTXT
            ELSE IF (CRDLRG .LT. D10000) THEN
               WRITE (MLINE(NC),9030) NAME,CRX,CRY,CRZ,TMPTXT
            ELSE IF (CRDLRG .LT. D100000) THEN
               WRITE (MLINE(NC),9020) NAME,CRX,CRY,CRZ,TMPTXT
            ELSE IF (CRDLRG .LT. D1000000) THEN
               WRITE (MLINE(NC),9010) NAME,CRX,CRY,CRZ,TMPTXT
            ELSE
               WRITE (MLINE(NC),9000) NAME,CRX,CRY,CRZ,TMPTXT
            END IF
         END IF
         COONEW(1,N) = CRX
         COONEW(2,N) = CRY
         COONEW(3,N) = CRZ
  100 CONTINUE
 9000 FORMAT (A4,3F20.10,1X,A11)
 9010 FORMAT (A4,3F20.11,1X,A11)
 9020 FORMAT (A4,3F20.12,1X,A11)
 9030 FORMAT (A4,3F20.13,1X,A11)
 9040 FORMAT (A4,3F20.14,1X,A11)
 9050 FORMAT (A4,3F20.15,1X,A11)
 9060 FORMAT (A4,3F20.16,1X,A11)
 9100 FORMAT (A4)
C
C Punch MOLECULE input with updated coordinates to LUMOL
C
      CALL GPOPEN(LUMOL,'ABACUS.MOL','UNKNOWN',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      DO 200 IMLINE = 1,NMLINE
         WRITE (LUMOL,'(A)') MLINE(IMLINE)
  200 CONTINUE
      CALL GPCLOSE(LUMOL,'KEEP')
C
C End of WLKMOL
C
      RETURN
      END
C  /* Deck wlkanh */
      SUBROUTINE WLKANH(ANHFAC,GRDDIA,EVAL,EVALX,EVEC,EVECX,STPDIX,
     &                  OVLPMA,ANHARM,TRSDIA,NCART,OLD,IPRINT)
C
C     Estimate anharmonicities based on information from this and
C     previous geometry, calculate anharmonic thresholds for each
C     mode
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (TRSMAX = 5.D0, D0EFF = 1.D-8)
C
      LOGICAL OLD
      DIMENSION EVEC(NCART,NCART), EVECX(NCART,NCART),
     *          OVLPMA(NCART,NCART), EVAL(NCART), EVALX(NCART),
     *          STPDIX(NCART), ANHARM(NCART), GRDDIA(NCART),
     *          TRSDIA(NCART)
C
      DO 100 I = 1, NCART
         TRSDIA(I) = TRSMAX
  100 CONTINUE
      IF (OLD) THEN
         IF (IPRINT .GT. 30) THEN
            CALL HEADER('EVAL in WLKANH',-1)
            CALL OUTPUT(EVAL,1,1,1,NCART,1,NCART,1,LUPRI)
            CALL HEADER('EVALX in WLKANH',-1)
            CALL OUTPUT(EVALX,1,1,1,NCART,1,NCART,1,LUPRI)
            CALL HEADER('STPDIX in WLKANH',-1)
            CALL OUTPUT(STPDIX,1,1,1,NCART,1,NCART,1,LUPRI)
         END IF
C
C        **********************************
C        ***** Calculate all overlaps *****
C        **********************************
C
         CALL DGEMM('T','N',NCART,NCART,NCART,1.D0,
     &              EVEC,NCART,
     &              EVECX,NCART,0.D0,
     &              OVLPMA,NCART)
         IF (IPRINT .GT. 30) THEN
            CALL HEADER('Overlap between old and new eigenvectors:',-1)
            CALL OUTPUT(OVLPMA,1,NCART,1,NCART,NCART,NCART,1,LUPRI)
         END IF
C
C        ************************************
C        ***** Estimate anharmonicities *****
C        ************************************
C
         DO 200 I = 1, NCART
            IOLD = IDAMAX(NCART,OVLPMA(I,1),NCART)
            IF (ABS(STPDIX(IOLD)) .GT. D0EFF) THEN
               ANHARM(I) = (EVAL(I) - EVALX(IOLD))/STPDIX(IOLD)
               IF (ABS(ANHARM(I)) .GT. D0EFF) THEN
                  HARSTP = WLKTRD(GRDDIA(I),EVAL(I),ANHARM(I),ANHFAC)
                  TRSDIA(I) = MIN(TRSMAX,HARSTP)
               END IF
            END IF
  200    CONTINUE
C
         IF (IPRINT .GT. 5) THEN
            CALL HEADER('ANHARM in WLKANH',-1)
            CALL OUTPUT(ANHARM,1,1,1,NCART,1,NCART,1,LUPRI)
            CALL HEADER('TRSDIA in WLKANH',-1)
            CALL OUTPUT(TRSDIA,1,1,1,NCART,1,NCART,1,LUPRI)
         END IF
      END IF
      RETURN
      END
C  /* Deck wlkrea */
      SUBROUTINE WLKREA(EVAL,EVEC,PMCAR,GRAD0,HESS0,DIPM0,DIPG0,COOR0,
     &                  GRADF,HESSF,DIPGF,ERGMOL,GRDNRM,ERGSEC,TRUSTR,
     &                  ESTART,TACCUM,DISPLC,STPSYM,STPNRM,STPDIX,
     &                  GRDDIX,EKIN0,TRUMAX,NCRTOT,N2CRT,NCRIND,NCART,
     &                  IWKTYP,IWKIND,IMODE,ISCTYP,IPOINT,ISTIND,ISTTYP,
     &                  NXYZ,DODIP,IMAGE,AACCUM,ROAAFD,ROAGND,
     &                  ROAGLD,ROAAD,ROAAFU,ROAGNU,ROAGLU,ROAAU,
     &                  FMATF,VIBAVE,SUSTO0,SUSTF,GFAC0,GFACF,QUAD0,
     &                  QUADF,SIGMA0,SIGMAF,SRC0,SRCF,POLAR0,POLARF,
     &                  EFG0,EFGF,SSJ0,SSJF,NUMMOD,IPRINT)
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <mxcent.h>
C
      LOGICAL DODIP,IMAGE,VIBAVE
      DIMENSION EVAL(NCART), EVEC(NCART,NCART), PMCAR(NCRTOT),
     &          GRAD0(NCART), HESS0(NCART,NCART), COOR0(NCRIND),
     &          GRADF(NCART), HESSF(NCART,NCART), STPSYM(NCART),
     &          DIPM0(3), DIPG0(3,NCART), DIPGF(3,NCART),
     &          STPDIX(NCART), GRDDIX(NCART),
     &          FMATF(NCART,NCART,NCART), SUSTO0(3,3), SUSTF(3,3,NCART),
     &          GFAC0(3,3), GFACF(3,3,NCART), QUAD0(3,3), 
     &          QUADF(3,3,NCART), SIGMA0(3,3,MXCENT), 
     &          SIGMAF(3,3,MXCENT,NCART), SRC0(3,3,MXCENT), 
     &          SRCF(3,3,MXCENT,NCART), POLAR0(3,3),
     &          POLARF(3,3,NCART),
     &          EFG0(3,3,MXCENT), EFGF(3,3,MXCENT,NCART), 
     &          SSJ0(MXCOOR,MXCOOR), SSJF(MXCOOR,MXCOOR,NCART)
C
#include <inftap.h>
#include <nuclei.h>
#include <cbilnr.h>
      DIMENSION ROAAFD(3,3,MXFR,MXCOOR), ROAGND(3,3,MXFR,MXCOOR), 
     &          ROAGLD(3,3,MXFR,MXCOOR), ROAAD(3,3,3,MXFR,MXCOOR),
     &          ROAAFU(3,3,MXFR), ROAGNU(3,3,MXFR), 
     &          ROAGLU(3,3,MXFR), ROAAU(3,3,3,MXFR)
C
      CALL GPOPEN(LUWLK,ABAWLK,'OLD','SEQUENTIAL','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUWLK
      READ (LUWLK) IWKTYP,IWKIND,IMODE,ISCTYP,
     &             ERGMOL,GRDNRM,ERGSEC,TRUSTR,STPNRM,
     &             IPOINT,ESTART,TACCUM,DISPLC,ISTIND,ISTTYP,IMAGE,
     &             AACCUM,TRUMAX
      READ (LUWLK) EVAL, EVEC, STPSYM
      IF (IWKTYP .EQ. 3) READ (LUWLK) PMCAR, STPDIX, GRDDIX, EKIN0
      IF (IWKTYP .EQ. 7) READ (LUWLK) STPDIX, GRDDIX, EKIN0
      IF (IWKTYP .EQ. 6) READ (LUWLK) COOR0, GRAD0, HESS0, GRADF, HESSF,
     &                                DIPM0, DIPG0, DIPGF, NXYZ, DODIP
      IF (IWKTYP .EQ. 6) READ (LUWLK) NFRVAL, FRVAL,
     &                                ROAAFU,ROAAFD,ROAGNU,ROAGND,
     &                                ROAGLU,ROAGLD, ROAAU, ROAAD
      IF (IWKTYP .EQ. 6) READ (LUWLK) NUMMOD, FMATF
      IF (IWKTYP .EQ. 6 .AND. VIBAVE) READ (LUWLK) SUSTO0, SUSTF, GFAC0,
     &                                GFACF, QUAD0, QUADF, SIGMA0, 
     &                                SIGMAF, SRC0, SRCF, POLAR0, 
     &                                POLARF, EFG0, EFGF, SSJ0, SSJF
      CALL GPCLOSE(LUWLK,'KEEP')
C
      IF (IPRINT .GT. 5) THEN
         CALL TITLER('Output from WLKREA','*',103)
         CALL AROUND('The following has been read in from LUWLK')
         WRITE (LUPRI,'(A,I5)')     ' IWKTYP ', IWKTYP
         WRITE (LUPRI,'(A,I5)')     ' IWKIND ', IWKIND
         WRITE (LUPRI,'(A,I5)')     ' IMODE  ', IMODE
         WRITE (LUPRI,'(A,I5)')     ' ISTIND ', ISTIND
         WRITE (LUPRI,'(A,F20.10)') ' ERGMOL ', ERGMOL
         WRITE (LUPRI,'(A,F20.10)') ' GRDNRM ', GRDNRM
         WRITE (LUPRI,'(A,F20.10)') ' ERGSEC ', ERGSEC
         WRITE (LUPRI,'(A,F20.10)') ' TRUSTR ', TRUSTR
         WRITE (LUPRI,'(A,F20.10)') ' STPNRM ', STPNRM
         WRITE (LUPRI,'(A,I5)')     ' IPOINT ', IPOINT
         WRITE (LUPRI,'(A,F20.10)') ' ESTART ', ESTART
         WRITE (LUPRI,'(A,F20.10)') ' TACCUM ', TACCUM
         WRITE (LUPRI,'(A,F20.10)') ' AACCUM ', AACCUM
         WRITE (LUPRI,'(A,F20.10)') ' DISPLC ', DISPLC
         IF (IPRINT .GT. 10) THEN
            CALL HEADER('Hessian eigenvalues',1)
            CALL OUTPUT(EVAL,1,1,1,NCART,1,NCART,1,LUPRI)
            CALL HEADER('Hessian eigenvectors',1)
            CALL OUTPUT(EVEC,1,NCART,1,NCART,NCART,NCART,1,LUPRI)
            CALL HEADER('Previous step vector',1)
            CALL OUTPUT(STPDIX,1,1,1,NCART,1,NCART,1,LUPRI)
            CALL HEADER('Previous gradient vector',1)
            CALL OUTPUT(GRDDIX,1,1,1,NCART,1,NCART,1,LUPRI)
         END IF
      END IF
      RETURN
      END
C  /* Deck wlkwrt */
      SUBROUTINE WLKWRT(EVAL,EVEC,PMCAR,GRAD0,HESS0,DIPM0,DIPG0,COOR0,
     &                  GRADF,HESSF,DIPGF,ERGMOL,GRDNRM,ERGSEC,TRUSTR,
     &                  ESTART,TACCUM,TIMSTP,EKINTO,EKINPA,COOR,GRDCAR,
     &                  HESCAR,STPCAR,DISPLC,STPSYM,STPNRM,STPDIA,
     &                  GRDDIA,EKIN0,TRUMAX,NCRTOT,N2CRT,NCRIND,NCART,
     &                  IWKTYP,IWKIND,IMODE,ISCTYP,ISTIND,ISTTYP,NXYZ,
     &                  DODIP,IMAGE,IPOINT,AACCUM,ROAAFD,ROAGND,
     &                  ROAGLD,ROAAD,ROAAFU,ROAGNU,ROAGLU,ROAAU,
     &                  FMATF,VIBAVE,SUSTO0,SUSTF,GFAC0,GFACF,QUAD0,
     &                  QUADF,SIGMA0,SIGMAF,SRC0,SRCF,POLAR0,POLARF,
     &                  EFG0,EFGF,SSJ0,SSJF,NUMMOD,WORK,LWORK)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <dummy.h>
C
      LOGICAL DODIP, IMAGE, VIBAVE
      DIMENSION EVAL(NCART), EVEC(NCART,NCART), PMCAR(NCRTOT),
     &          COOR(NCRIND), GRDCAR(NCRTOT), HESCAR(N2CRT),
     &          GRAD0(NCART), HESS0(NCART,NCART), COOR0(NCRIND),
     &          GRADF(NCART), HESSF(NCART,NCART), STPSYM(NCART),
     &          DIPM0(3), DIPG0(3,NCART), DIPGF(3,NCART),
     &          STPDIA(NCART),GRDDIA(NCART),
     &          FMATF(NCART,NCART,NCART), SUSTO0(3,3), SUSTF(3,3,NCART),
     &          GFAC0(3,3), GFACF(3,3,NCART), QUAD0(3,3),
     &          QUADF(3,3,NCART), SIGMA0(3,3,MXCENT), 
     &          SIGMAF(3,3,MXCENT,NCART), SRC0(3,3,MXCENT), 
     &          SRCF(3,3,MXCENT,NCART), POLAR0(3,3), POLARF(3,3,NCART),
     &          EFG0(3,3,MXCENT), EFGF(3,3,MXCENT,NCART), 
     &          SSJ0(MXCOOR,MXCOOR), SSJF(MXCOOR,MXCOOR,NCART),
     &          WORK(LWORK)
C
#include <inftap.h>
#include <nuclei.h>
#include <cbilnr.h>
      DIMENSION ROAAFD(3,3,MXFR,MXCOOR), ROAGND(3,3,MXFR,MXCOOR), 
     &          ROAGLD(3,3,MXFR,MXCOOR), ROAAD(3,3,3,MXFR,MXCOOR),
     &          ROAAFU(3,3,MXFR), ROAGNU(3,3,MXFR), 
     &          ROAGLU(3,3,MXFR), ROAAU(3,3,3,MXFR)
C
      CALL GPOPEN(LUWLK,ABAWLK,'UNKNOWN','SEQUENTIAL','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      REWIND LUWLK
      WRITE (LUWLK) IWKTYP,IWKIND,IMODE,ISCTYP,
     &              ERGMOL,GRDNRM,ERGSEC,TRUSTR,STPNRM,
     &              IPOINT,ESTART,TACCUM,DISPLC,ISTIND,ISTTYP,IMAGE,
     &              AACCUM,TRUMAX
      WRITE (LUWLK) EVAL, EVEC, STPSYM
      IF (IWKTYP .EQ. 3) WRITE (LUWLK) PMCAR, STPDIA, GRDDIA, EKIN0
      IF (IWKTYP .EQ. 7) WRITE (LUWLK) STPDIA, GRDDIA, EKIN0
      IF (IWKTYP .EQ. 6) WRITE (LUWLK) COOR0, GRAD0, HESS0, GRADF,HESSF,
     &                                 DIPM0, DIPG0, DIPGF, NXYZ, DODIP
      IF (IWKTYP .EQ. 6) WRITE (LUWLK) NFRVAL, FRVAL,
     &                                 ROAAFU,ROAAFD,ROAGNU,ROAGND,
     &                                 ROAGLU,ROAGLD,ROAAU, ROAAD
      IF (IWKTYP .EQ. 6) WRITE (LUWLK) NUMMOD, FMATF
      IF (IWKTYP .EQ. 6 .AND. VIBAVE)  WRITE (LUWLK) SUSTO0, SUSTF, 
     &                                 GFAC0, GFACF, QUAD0, QUADF, 
     &                                 SIGMA0, SIGMAF, SRC0, SRCF,
     &                                 POLAR0, POLARF, EFG0, EFGF,
     &                                 SSJ0, SSJF
      CALL GPCLOSE(LUWLK,'KEEP')
C
C     Write information on LUTRJ
C
      IF (IWKTYP .EQ. 3) THEN
         CALL GPOPEN(LUWLK,ABATRJ,'UNKNOWN','SEQUENTIAL','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         REWIND LUWLK
         IF (IPOINT .EQ. 0) THEN
            CALL NEWLAB('POINT   ',LUWLK,LUPRI)
         ELSE
            DO 100 ITER = 0, IPOINT
               CALL MOLLAB('POINT   ',LUWLK,LUPRI)
  100       CONTINUE
         END IF
         WRITE (LUWLK) IPOINT, DUMMY, DUMMY, DUMMY, DUMMY
         WRITE (LUWLK) ERGMOL, ERGSEC, EKINTO, EKINPA, TIMSTP
         WRITE (LUWLK) COOR, GRDCAR, HESCAR, STPCAR, PMCAR
         CALL NEWLAB('POINT   ',LUWLK,LUPRI)
         CALL GPCLOSE(LUWLK,'KEEP')
C
C     Write information on LUIRC
C
      ELSE IF (IWKTYP .EQ. 7) THEN
         CALL GPOPEN(LUWLK,ABAIRC,'UNKNOWN','SEQUENTIAL','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         REWIND LUWLK
         IF (IPOINT .EQ. 0) THEN
            CALL NEWLAB('POINT   ',LUWLK,LUPRI)
         ELSE
            DO 200 ITER = 0, IPOINT
               CALL MOLLAB('POINT   ',LUWLK,LUPRI)
  200       CONTINUE
         END IF
         CALL WLKIRW(IPOINT,AACCUM,STPNRM,WORK,LWORK,IPRINT)
         CALL NEWLAB('POINT   ',LUWLK,LUPRI)
         CALL GPCLOSE(LUWLK,'KEEP')
      END IF
      RETURN
      END
C  /* Deck wlkirw */
      SUBROUTINE WLKIRW(IPOINT,AACCUM,STPNRM,WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      DIMENSION WORK(LWORK)
#include <nuclei.h>
      KMASS  = 1
      KCOOR  = KMASS  + NUCDEP
      KGRAD  = KCOOR  + MXCOOR
      KHESS  = KGRAD  + MXCOOR
      KCSTRA = KHESS  + MXCOOR*MXCOOR
      KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
      KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
      IF (KLAST .GT. LWORK) CALL STOPIT('WLKIRW',' ',KLAST,LWORK)
      CALL WLKIR1(IPOINT,AACCUM,STPNRM,WORK(KCOOR),WORK(KGRAD),
     &            WORK(KHESS),WORK(KMASS),WORK(KCSTRA),
     &            WORK(KSCTRA),IPRINT)
      RETURN
      END
C  /* Deck wlkir1 */
      SUBROUTINE WLKIR1(IPOINT,AACCUM,STPNRM,COOR,CGRAD,CHESS,
     &                  DMASS,CSTRA,SCTRA,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <inftap.h>
#include <symmet.h>
#include <nuclei.h>
#include <taymol.h>
      DIMENSION COOR(MXCOOR), CGRAD(MXCOOR), CHESS(MXCOOR,MXCOOR),
     &          DMASS(NUCDEP), CSTRA(*), SCTRA(*)
#include <ibtfun.h>
C
      NCOOR = 3*NUCDEP
      CALL DZERO(COOR, NCOOR)
      CALL DZERO(CGRAD,NCOOR)
      CALL DZERO(CHESS,NCOOR*NCOOR)
C
      ICOOR = 0
      IATOM = 0
      DO 100 ICENT = 1, NUCIND
         DMSS = DISOTP(IZATOM(ICENT),ISOTOP(ICENT),'MASS')
         MULCNT = ISTBNU(ICENT)
         DO 200 IOP = 0, MAXOPR
            IF (IBTAND(IOP,MULCNT) .EQ. 0) THEN
               IATOM = IATOM + 1
               DMASS(IATOM) = DMSS
               DO 300 I = 1, 3
                 ICOOR = ICOOR + 1
                 COOR(ICOOR) = PT(IBTAND(ISYMAX(I,1),IOP))*CORD(I,ICENT)
  300          CONTINUE
            END IF
  200    CONTINUE
  100 CONTINUE
      CALL TRAGRD(GRDMOL,CGRAD,CSTRA,SCTRA,NCRREP(0,1),NCOOR)
      CALL TRAHES(HESMOL,MXCOOR,CHESS,CSTRA,SCTRA,MXCOOR,NCOOR,1)
      ARC = AACCUM - STPNRM
      WRITE (LUWLK) IPOINT, NCOOR, DUMMY, DUMMY, DUMMY,
     &     DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY
      WRITE (LUWLK) ARC, ERGMOL, DUMMY, DUMMY,
     &     DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY
      WRITE (LUWLK) (DMASS (I),I=1,NUCDEP),
     &     DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY
      WRITE (LUWLK) (COOR (I),I=1,NCOOR)
      WRITE (LUWLK) (CGRAD(I),I=1,NCOOR)
      WRITE (LUWLK) ((CHESS(I,J),I=1,NCOOR),J=1,NCOOR)
      RETURN
      END
C  /* Deck wlkres */
      SUBROUTINE WLKRES(STPDIA,STPCAR,GRDDIA,GRDNRM,EVAL,EVEC,CNDHES,
     &                  STPSYM,SCAL,COOR,RATIO,RNU,TRUSTX,TRUSTR,ERGSEC,
     &                  OVLPMA,OVLPWM,COFMOD,ANHRAD,STPSYX,ESTART,
     &                  ERGMOL,ERGMLX,ERGSCX,TACCUM,TIMSTP,EKINTO,
     &                  EKINPA,VREL,GRDCAR,PMCAR,PMCARX,PMDIA,PMDIAX,
     &                  CMDIF,ANHARM,TRSDIA,EKIN0,STPNRM,IMODE,IWKIND,
     &                  DOREPW,NCRTOT,N2CRT,NCRIND,NCART,INDHES,IPOINT,
     &                  ISTTYP,ISTIND,IWKTYP,STOP,STATPO,OLDREA,DOHESS,
     &                  PRJTRO,IMAGE,AACCUM,DONWTN,
     &                  IPRINT,REJECT,ERGDIA,WORK,LWORK)
C
C     Print of results
C
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <codata.h>
#include <facang.h>
      PARAMETER (D1 = 1.0D0, D0 = 0.0D0, D100 =100.D0, D15 = 1.0D15,
     &           DM13 = 1.0D-13)
C
      LOGICAL STATPO, STOP, DOREPW(0:7), OLDREA, DOHESS, NOTALL, PRJTRO,
     &        IMAGE, DONWTN, REJECT
      DIMENSION STPDIA(NCRTOT), EVAL(NCRTOT), EVEC(N2CRT), COOR(NCRIND),
     &          STPCAR(NCRIND), CNDHES(0:7), INDHES(0:7),
     &          GRDDIA(NCRTOT), SCAL(NCRTOT), STPSYM(NCART),
     &          OVLPMA(NCART,NCART), STPSYX(NCART), PMCAR(NCART),
     &          PMDIA(NCART), PMDIAX(NCART), PMCARX(NCART),
     &          ANHARM(NCART), TRSDIA(NCART), GRDCAR(NCART),
     &          ERGDIA(NCRTOT), WORK(LWORK)
C
#include <nuclei.h>
#include <symmet.h>
#include <trkoor.h>
#include <prkoor.h>
C
      CHARACTER*1 STARM, STARE, STARG, STARS, STARH
      CHARACTER*42 WLKTYP(0:8)
      DATA WLKTYP /'trust-region minimization                 ',
     *             'mode-following                            ',
     *             'gradient extremal                         ',
     *             'dynamic walk                              ',
     *             'Newton step                               ',
     *             'eigenvector step                          ',
     *             'numerical differentiation                 ',
     &             'IRC path                                  ',
     &             'TS by minimization on image surface       '/
C
      CALL TITLER('FINAL RESULTS FROM WALK','*',125)
C
      IPRTYP = IWKTYP
      IF (IWKTYP .EQ. 1 .AND. IWKIND .EQ. 0) IPRTYP = 0
      IF (IMAGE) IPRTYP = 8
      WRITE (LUPRI,'(/2A)') ' Walk type: ', WLKTYP(IPRTYP)
      WRITE (LUPRI,'(/A,I3)') ' Iteration number:', IPOINT
      IF (REJECT) WRITE (LUPRI,'(/A)')
     &      '@ Previous step was rejected - attempting shorter step.'
      WRITE (LUPRI,'(/A,F17.10)') '  Norm of molecular gradient:',GRDNRM
      WRITE (LUPRI,'(/A,F17.10)') '@ Energy at current point:', ERGMOL
      IF (IPOINT .GT. 0) THEN
         WRITE (LUPRI,'( A,F17.10)')
     &                      '  Energy at initial point:', ESTART
         WRITE (LUPRI,'(/A,F13.10,A,2(/,39X,F13.2,A))')
     &         '  Change in energy since start of walk:',
     &                   ERGMOL - ESTART,  ' au',
     &                   XKJMOL*(ERGMOL - ESTART), ' kJ/mol',
     &                   XKCMOL*(ERGMOL - ESTART), ' kCal/mol'
         IF (.NOT.REJECT) WRITE (LUPRI,'(/A,F6.2)')
     &      '@ Ratio between actual and predicted energies:',
     &       (ERGMOL - ERGMLX)/ERGSCX
      END IF
C
      IF (IWKTYP .EQ. 3 .OR. IWKTYP .EQ. 7) THEN
C
         CALL HEADER('Atomic masses',1)
         JATOM = 0
         DO 510 IATOM = 1, NUCIND
            DMASS = DISOTP(IZATOM(IATOM),ISOTOP(IATOM),'MASS')
            NDEG  = NUCDEG(IATOM)
            DO 515 IDEG = 1, NDEG
               WRITE (LUPRI,'(2X,A6,2X,F12.6)') NAMDEP(JATOM+IDEG),DMASS
  515       CONTINUE
            JATOM = JATOM + NDEG
  510    CONTINUE
      END IF
C
      IF (IWKTYP .EQ. 3) THEN
         IF (CMDIF .GT. D0) WRITE (LUPRI,'(/A,2(F12.6,A))')
     &      '@ Distance between centers of mass ', CMDIF,' au  ',
     &      XTANG*CMDIF,' Angstrom'
         WRITE (LUPRI,'(/A,F12.2,A,1P,E12.2,A)')
     &      '@ Accumulated time before step: ',
     &        TACCUM - TIMSTP, ' au  ',D15*XFSEC*(TACCUM - TIMSTP),' fs'
         EKINX = EKIN0 + ESTART - ERGMOL
         WRITE (LUPRI,'(/A,1P,E12.6)')
     &     ' Kinetic energy before step - exact       : ', EKINX
         IF (IPOINT .GT. 0) THEN
            WRITE (LUPRI,'(/A,1P,E12.6)')
     &        ' Kinetic energy at initial point          : ', EKIN0
            WRITE (LUPRI,'(/A,1P,E12.6)')
     &        ' Accumulated kinetic energy since start   : ',
     &          EKINX - EKIN0
            ERROR = D100*(EKINX - EKINTO)/EKINX
            WRITE (LUPRI,'( A,1P,E12.6,A,0P,F5.1,A)')
     &        '                            - from momenta: ', EKINTO,
     &        ' (',ERROR,'% error)'
         END IF
         IF (ABS(VREL) .GT. D0) WRITE(LUPRI,'(/A,1P,E12.6,A,0P,F6.2,A)')
     &      '@ Relative velocity before step            : ', ABS(VREL),
     &      ' (',DM13*ABS(VREL)*XTANG/XFSEC, 'km/s)'
         IF (ABS(EKINPA) .GT. D0) THEN
            WRITE (LUPRI,'( A,1P,E12.6,A,0P,F5.1,A)')
     &      '@ Relative translational energy before step: ', EKINPA,
     &      ' (',D100*EKINPA/EKINTO,'% of total)'
         END IF
C
C        Analysis
C
         CALL WLKDAN(PMCARX,NCART,WORK,LWORK,IPRINT)
      END IF
      IF (IWKTYP .EQ. 7) THEN
         WRITE (LUPRI,'(/A,F10.4,A)')
     &      '@ Accumulated arc length (mass weighted) before step: ',
     &        AACCUM - STPNRM, ' au  '
         IF (STOP) THEN
            WRITE (LUPRI,'(/A)')'@ Geometry optimization has converged.'
         END IF
      END IF
      IF (DOHESS) THEN
         INDEX  = ISUM(MAXREP+1,INDHES(0),1)
         IF (MAXREP .GT. 0) THEN
            NOTALL = .FALSE.
            DO 50 IREP = 0, MAXREP
               IF (.NOT.DOREPW(IREP)) THEN
                  IF (NCRREP(IREP,1)-NPRREP(IREP).GT.0) NOTALL=.TRUE.
               END IF
  50        CONTINUE
            WRITE (LUPRI,'(/A,I2)')
     *        '@ Hessian index of totally symmetric representation:',
     *          INDHES(0)
            IF (NOTALL) THEN
               WRITE (LUPRI,'(A)')
     *          '  Hessian index of all representations             : ?'
            ELSE
               WRITE (LUPRI,'(A,I2)')
     *          '@ Hessian index of all representations             :',
     *             INDEX
            END IF
         ELSE
            WRITE (LUPRI,'(/A,I2)') '@ Hessian index:', INDEX
         END IF
         IF (IMAGE .AND. .NOT.STOP) THEN
            IMGIND = INDHES(0)
            IF (EVAL(IMODE).GT.D0) THEN
               IMGIND = IMGIND + 1
            ELSE
               IMGIND = IMGIND - 1
            END IF
            WRITE (LUPRI,'(A,I2)') '@ Image index:  ', IMGIND
         END IF
      END IF
      IF (IWKTYP .LT. 3) THEN
         IF (STOP) THEN
            WRITE (LUPRI,'(/A)')'@ Geometry optimization has converged.'
         ELSE
            IF (IWKTYP .EQ. 1 .AND. (IWKIND .GT. 0 .OR. IMAGE)) THEN
               WRITE (LUPRI,'(A,I2)') '  Initial index:',ISTIND
               IF (IMAGE) THEN
                  IWKINX = 1
               ELSE
                  IWKINX = IWKIND
               END IF
               WRITE (LUPRI,'(A,I2)') '  Target index: ',IWKINX
            ELSE
               WRITE (LUPRI,'(/A,I2)')'  Index of initial point:',ISTIND
            END IF
         END IF
         IF (IMODE .GT. 0) THEN
            IF (IMAGE) THEN
               IF (.NOT.STOP) WRITE (LUPRI,'(/A,I2,/)')
     &                  '@ Image mode:',IMODE
            ELSE
               WRITE (LUPRI,'(/A,I2,/)') '@ Eigenmode followed:',IMODE
            END IF
         END IF
      ELSE IF (IWKTYP .EQ. 5) THEN
         WRITE (LUPRI,'(/A,I2)') '@ Eigenvector selected:',IMODE
      END IF
      IF (.NOT.STOP) THEN
         IF (IWKTYP.EQ.2) THEN
            STPPAR = ABS(STPDIA(IMODE))
            STPPER = SQRT(STPNRM*STPNRM - STPPAR*STPPAR)
            IF (STPNRM .GT. 1.0D-10) WRITE (LUPRI,'(/A,F6.2,A)')
     *          '@ Angle between step and walk mode:    ',
     *            FACANG*ACOS(ABS(STPPAR/STPNRM)), ' deg '
         END IF
         IF ((IWKTYP.EQ.1 .AND. IWKIND.GT.0) .OR. (IWKTYP.EQ.2)
     &                                       .OR. IMAGE) THEN
            IF (OLDREA) THEN
               WRITE (LUPRI,'(A,F6.2,A)')
     *            '@ Angle between gradient and walk mode:',
     *               FACANG*ACOS(ABS(COFMOD)), ' deg '
               WRITE (LUPRI,'(A,F6.2,A)')
     *            '@ Angle between old and new walk modes:',
     *              FACANG*ACOS(ABS(MIN(OVLPWM,D1))), ' deg '
               DOTNO = DDOT(NCART,STPSYM,1,STPSYX,1)
               DNRMN = DNRM2(NCART,STPSYM,1)
               DNRMO = DNRM2(NCART,STPSYX,1)
               WRITE (LUPRI,'(A,F6.2,A)')
     *            '@ Angle between old and new steps:     ',
     *              FACANG*ACOS(ABS(DOTNO/(DNRMN*DNRMO))), ' deg '
               IF (IPRINT .GT. 5) THEN
                  CALL HEADER
     *               ('Overlaps between old and new eigenvectors:',-1)
                  CALL OUTPUT(OVLPMA,1,NCART,1,NCART,NCART,NCART,1,
     *                        LUPRI)
               END IF
            END IF
         END IF
      END IF
C
      IF (DOHESS) THEN
         IF (PRJTRO) THEN
            WRITE (LUPRI,'(/2(/2X,A))')
     *            ' sym        #int        #ext       '
     *            //'gradient     index       cond. num.',
     *            ' ----------------------------------'
     *            //'-----------------------------------'
         ELSE
            WRITE (LUPRI,'(/2(/2X,A))')
     *          ' sym        #int        #ext       gradient     index',
     *          ' -----------------------------------------------------'
         END IF
         DO 100 IREP = 0, MAXREP
            IF (DOREPW(IREP)) THEN
               IF (IREP .EQ. 0) THEN
                  GRADNR = GRDNRM
               ELSE
                  GRADNR = D0
               END IF
               IF (PRJTRO) THEN
                  WRITE (LUPRI,'(3(3X,I2,7X),1P,D12.6,3X,I2,10X,D8.2)')
     *                  IREP + 1, NCRREP(IREP,1) - NPRREP(IREP),
     *                  NPRREP(IREP), GRADNR, INDHES(IREP), CNDHES(IREP)
               ELSE
                  WRITE (LUPRI,'(3(3X,I2,7X),1P,D12.6,3X,I2)')
     *                  IREP + 1, NCRREP(IREP,1) - NPRREP(IREP),
     *                  NPRREP(IREP), GRADNR, INDHES(IREP)
               END IF
            ELSE
               WRITE (LUPRI,'(3(3X,I2,7X))')
     *                IREP + 1, NCRREP(IREP,1) - NPRREP(IREP),
     *                NPRREP(IREP)
            END IF
  100    CONTINUE
      END IF
C
      IF (DOHESS .AND. IWKTYP .NE. 3) THEN
         WRITE (LUPRI,'()')
         IF (STOP) THEN
            CALL AROUND
     *         ('Gradient and Hessian in diagonal representation')
            WRITE (LUPRI,'(2(/10X,A))')
     *      ' @sym        mode          gradient        eigenvalue    ',
     *      ' @-------------------------------------------------------'
         ELSE
            CALL AROUND('Analysis in diagonal representation')
            WRITE (LUPRI,'(2(/A))') '@  iter   sym     mode      '
     *        //'gradient     eigenvalue       step        energy ',
     *        '@  -------------------------'
     *        //'---------------------------------------------------'
         END IF
         DO 150 I = 1, NCRTOT
            ERGDIA(I) = WLKTAY(GRDDIA(I),EVAL(I),STPDIA(I))
  150    CONTINUE
         MAXGRD = IDAMAX(NCRTOT,GRDDIA,1)
         MAXHES = IDAMAX(NCRTOT,EVAL  ,1)
         MAXSTP = IDAMAX(NCRTOT,STPDIA,1)
         MAXERG = IDAMAX(NCRTOT,ERGDIA,1)
         IOFF = 0
         DO 200 IREP = 0, MAXREP
            IF (DOREPW(IREP)) THEN
               IF (PRJTRO) THEN
                  NVEC = NCRREP(IREP,1) - NPRREP(IREP)
               ELSE
                  NVEC = NCRREP(IREP,1)
               END IF
               DO 210 I = 1, NVEC
                  II = IOFF + I
                  IF (STOP) THEN
                     WRITE (LUPRI,'(12X,A,I2,I12,4X,2F17.10)')
     &                  '@',IREP + 1, I, GRDDIA(II), EVAL(II)
                  ELSE
                     STARM = ' '
                     STARG = ' '
                     STARH = ' '
                     STARS = ' '
                     STARE = ' '
                     IF (II .EQ. IMODE)  STARM = '<'
                     IF (II .EQ. MAXGRD) STARG = '<'
                     IF (II .EQ. MAXHES) STARH = '<'
                     IF (II .EQ. MAXSTP) STARS = '<'
                     IF (II .EQ. MAXERG) STARE = '<'
                     WRITE
     &                  (LUPRI,'(A,I5,I6,I8,A,2X,4(F13.8,A))')
     &                  '@',IPOINT, IREP + 1, I,STARM,GRDDIA(II),STARG,
     &                  EVAL(II),STARH, STPDIA(II),STARS,
     &                  ERGDIA(II),STARE
                  END IF
  210          CONTINUE
               IOFF = IOFF +  NCRREP(IREP,1)
            END IF
  200    CONTINUE
         IF ((ISTTYP .EQ. 4) .OR. (IWKTYP .EQ. 7 .AND. DONWTN)) THEN
            WRITE (LUPRI,'(/A)') '  This is Newton step.'
         END IF
         IF (IWKTYP .EQ. 7 .AND. .NOT.STOP .AND. IPOINT.GT.0) THEN
            ANGLE = WINKEL(GRDDIA,STPDIA,NCART)
            WRITE (LUPRI,'(/A,F6.2,A)')
     &         '  Angle between gradient and step: ',ANGLE,' degrees'
         END IF
         IF (.NOT.STOP .AND. IWKTYP.EQ.2) THEN
            WRITE (LUPRI,'(3(/,A,F12.6))')
     *         '  Step size:                   ', STPNRM,
     *         '     - parallel component:     ', STPPAR,
     *         '     - perpendicular component:', STPPER
         END IF
      ELSE IF (IWKTYP .EQ. 3) THEN
         CALL AROUND('Cartesian gradient and momenta')
         DO 40 I = 1, NCRREP(0,1)
            WRITE (LUPRI,'(18X,A6,2F17.10)')
     &            NAMEX(IPTCOR(I,1)), GRDCAR(I), PMCARX(I)
   40    CONTINUE
         WRITE (LUPRI,'()')
         CALL AROUND('Analysis in diagonal representation')
         WRITE (LUPRI,'(2(/A))') '  sym mode    grad'
     *     //'     Hessian    momentum     step     incr mom    energy',
     *       '  ----------------'
     *     //'---------------------------------------------------------'
         IOFF = 0
         DO 220 IREP = 0, MAXREP
            IF (DOREPW(IREP)) THEN
               IF (PRJTRO) THEN
                  NVEC = NCRREP(IREP,1) - NPRREP(IREP)
               ELSE
                  NVEC = NCRREP(IREP,1)
               END IF
               DO 230 I = 1, NVEC
                  II = IOFF + I
                  WRITE (LUPRI,'(2X,I2,2X,I2,1X,6F11.6)')
     &               IREP + 1, I, GRDDIA(II), EVAL(II), PMDIAX(II),
     &               STPDIA(II), PMDIA(II) - PMDIAX(II),
     &               WLKTAY(GRDDIA(II),EVAL(II),STPDIA(II))
  230          CONTINUE
               IOFF = IOFF +  NCRREP(IREP,1)
            END IF
  220    CONTINUE
         WRITE (LUPRI,'(/,A,F12.6)')
     &      ' Step size:', STPNRM
         WRITE (LUPRI,'(/A,F12.2,A,1P,E12.2,A)')
     &      ' Time step                  :',
     &      TIMSTP, ' au ',D15*XFSEC*TIMSTP,' fs '
         WRITE (LUPRI,'( A,F12.2,A,1P,E12.2,A)')
     &      ' Accumulated time after step:',
     &      TACCUM, ' au ',D15*XFSEC*TACCUM,' fs '
      END IF
      IF (IWKTYP .EQ. 7) THEN
         WRITE (LUPRI,'(/,A,F10.6,A)')
     &      '  Arc length (mass-weighted) of step:               ',
     &         STPNRM, ' au '
         WRITE (LUPRI,'(A,F10.6,A)')
     &      '  Accumulated arc length (mass-weighted) after step:',
     &         AACCUM, ' au '
      END IF
      IF (.NOT.STOP .AND. (IWKTYP .EQ. 3 .OR. IWKTYP .EQ. 7)
     &              .AND. (IPOINT .GT. 0)) THEN
         CALL HEADER(' mode   anharmonic estimate   harmonic region',1)
         DO 240 I = 1, NCRREP(0,1)
            WRITE (LUPRI,'(2X,I2,2X,2F24.12)')
     &      I, ANHARM(I), TRSDIA(I)
  240    CONTINUE
      END IF
C
C     ***** Print eigenvectors *****
C
      IF (DOHESS) THEN
         CALL AROUND('Hessian eigenvectors')
         IOFF  = 0
         IEVEC = - 1
         DO 300 IREP = 0, MAXREP
            NCR = NCRREP(IREP,1)
            IF (DOREPW(IREP)) THEN
               IF (PRJTRO) THEN
                  NONTRO = NCR - NPRREP(IREP)
               ELSE
                  NONTRO = NCR
               END IF
               NBATCH = (NONTRO + 4)/5
               ISTR = 1
               IEVC = IEVEC
               DO 310 IBATCH = 1, NBATCH
                  IEND = MIN(ISTR + 4,NONTRO)
                  NUMB = IEND - ISTR + 1
                  NEND = 1 + (NUMB - 1)*NCR
                  WRITE (LUPRI,'(/,1X,A,I1,4X,5(4X,I5,3X))')
     *              '  sym ',IREP + 1, (I,I = ISTR,IEND)
                  WRITE (LUPRI,'(2X,70A1,/)') ('-',I=1,10+NUMB*12)
                  DO 320 ICOOR = 1, NCR
                     WRITE (LUPRI,'(3X,A,(T13,5F12.6))')
     *                   NAMEX(IPTCOR(IOFF + ICOOR,1)),
     *                  (EVEC(IEVC + ICOOR + I),I=1,NEND,NCR)
                     IF (MAXREP+MOD(ICOOR,3) .EQ. 0) WRITE (LUPRI,'()')
 320              CONTINUE
                  ISTR = ISTR + 5
                  IEVC = IEVC + 5*NCR
 310           CONTINUE
               IEVEC = IEVEC + NCR*NCR
               WRITE (LUPRI,'()')
            END IF
            IOFF = IOFF + NCR
 300     CONTINUE
      END IF
C
C     *****************************
C     ***** Geometry analysis *****
C     *****************************
C
      IF (STOP) THEN
         CALL AROUND('Final geometry')
         CALL HEADER('Cartesian coordinates (au)',-1)
         CALL PRIGEO(CORD)
         CALL GEOANA(CORD,.TRUE.,.FALSE.,NBONDS,-1,WORK,LWORK)
      ELSE
         IF (OLDREA) THEN
            WRITE(LUPRI,'(/,A,F12.6)')
     *       ' Old trust radius:          ', TRUSTX
         ELSE
            WRITE(LUPRI,'(/,A,F12.6)')
     *       ' Input trust radius:        ', TRUSTX
         END IF
         WRITE (LUPRI,'(A,F12.6)')
     *       ' Updated trust radius:      ', TRUSTR
         IF (OLDREA.AND.
     &      (((IWKIND.GT.0).AND.(IWKTYP.LT.3)).AND..NOT.IMAGE)) THEN
            IF (ANHRAD .LT. D100) WRITE (LUPRI,'(A,F12.6)')
     *          ' Finite-difference estimate of harmonic radius',ANHRAD
         END IF
         IF (IWKTYP .EQ. 1) THEN
            WRITE (LUPRI,'(A,F15.10)')
     &          ' Level shift employed:       ', RNU
         END IF
         WRITE (LUPRI,'(//A,F15.10)')
     *       ' Predicted energy change ', ERGSEC
         CALL AROUND('Change in geometry')
         CALL HEADER('Cartesian displacements (au)',-1)
         CALL PRIGEO(STPCAR)
         CALL GEOANA(CORD,.FALSE.,.FALSE.,NBONDS,-1,WORK,LWORK)
         CALL GEOANA(COOR,.TRUE. ,.TRUE. ,NBONDS,-1,WORK,LWORK)
         CALL AROUND('Geometry after step')
         CALL HEADER('Cartesian coordinates (au)',-1)
         CALL PRIGEO(COOR)
         CALL GEOANA(COOR,.TRUE.,.FALSE.,NBONDS,-1,WORK,LWORK)
      END IF
C
      RETURN
      END
C  /* Deck wlkend */
      LOGICAL FUNCTION WLKEND(GRDCAR,GRDNRM,NCART,INDHES,INDEX,
     *                        INDXOK,STATPO,ISTIND,IPOINT,NUMMOD)
#include <implicit.h>
#include <priunit.h>
#include <gnrinf.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
C
      LOGICAL INDXOK, STATPO
      DIMENSION GRDCAR(NCART), INDHES(0:7)
C
#include <cbiwlk.h>
#include <symmet.h>
C
      INDEX  = ISUM(MAXREP+1,INDHES(0),1)
      IF (START) ISTIND = INDHES(0)
C
      IF (IWKTYP .EQ. 2) THEN
         INDXOK = INDHES(0) .NE. ISTIND
      ELSE IF (IMAGE) THEN
         INDXOK = INDHES(0) .EQ. 1
      ELSE IF (IWKTYP .EQ. 7) THEN
         INDXOK = INDHES(0) .EQ. 0
      ELSE
         INDXOK = INDHES(0) .EQ. IWKIND
      END IF
      GRDNRM = SQRT(DDOT(NCART, GRDCAR,1,GRDCAR,1))
      STATPO = GRDNRM .LT. TOLST
C
      IF (IWKTYP .EQ. 3)THEN
         WLKEND = IPOINT .EQ. ITERMX
         GEOCNV = WLKEND
      ELSE IF (IWKTYP .EQ. 4) THEN
         WLKEND = .FALSE.
         GEOCNV = WLKEND
      ELSE IF (IWKTYP .EQ. 5) THEN
         WLKEND = .FALSE.
         GEOCNV = WLKEND
      ELSE IF (IWKTYP .EQ. 6) THEN
         IF (VIBAVE .OR. NMODIF) THEN
            IF (ACCURT) THEN
               WLKEND = IPOINT .EQ. 4*NUMMOD - 1
            ELSE
               WLKEND = IPOINT .EQ. 2*NUMMOD - 1
            END IF
         ELSE
            WLKEND = IPOINT .EQ. 2*NCART - 1
         END IF
         GEOCNV = WLKEND
      ELSE IF (IWKTYP .EQ. 7) THEN
         WLKEND = (IPOINT .EQ. ITERMX - 1) .OR. (STATPO .AND. INDXOK)
         GEOCNV = WLKEND
      ELSE
         WLKEND = .NOT.START .AND. STATPO .AND. INDXOK
         GEOCNV = STATPO .AND. INDXOK
      END IF
C
      IF (IPRWLK .GT. 2) THEN
         WRITE (LUPRI,'(/A,1P,D12.6)') ' Norm of gradient   ',GRDNRM
        IF (IWKTYP .NE. 6) THEN
C        ... Hessian index not calculated for numerical differentiation
         WRITE (LUPRI,'(A,I5)')        ' Totally sym.  index',INDHES(0)
         WRITE (LUPRI,'(A,I5)')        ' Total Hessian index',INDEX
        END IF
         WRITE (LUPRI,'(A,I5)')        ' Number of steps    ',IPOINT
         WRITE (LUPRI,'(A,L5)')        ' WLKEND             ',WLKEND
      END IF
C
C     GRADML saves molecular gradient.
C
      GRADML = GRDNRM
      RETURN
      END
C  /* Deck wlkkin */
      SUBROUTINE WLKKIN(EKINTO,PMCAR,SCAL,NCART,IPRINT)
C
C     Total kinetic energy from momenta
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.D0, DP5 = 0.5D0)
      DIMENSION PMCAR(NCART), SCAL(NCART)
      EKINTO = D0
      DO 100 I = 1, NCART
         EKINTO = EKINTO + PMCAR(I)*PMCAR(I)*SCAL(I)*SCAL(I)
  100 CONTINUE
      EKINTO = DP5*EKINTO
C
C     Print
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('PMCAR in WLKKIN',-1)
         CALL OUTPUT(PMCAR,1,1,1,NCART,1,NCART,1,LUPRI)
         CALL HEADER('SCAL in WLKKIN',-1)
         CALL OUTPUT(SCAL,1,1,1,NCART,1,NCART,1,LUPRI)
         CALL HEADER('Output from WLKKIN',-1)
         WRITE (LUPRI,'(1X,A,1P,E24.12)') ' EKINTO ', EKINTO
      END IF
      RETURN
      END
C  /* Deck wlktra */
      SUBROUTINE WLKTRA(PMCAR,EKINPA,VREL,SCAL,CMDIF,NCART,IPRINT)
C
C      Purpose: Perform analysis of kinetic vs. potential energy.
C               Compute kinetic energy released.
C
C     eu, Oslo 16.February 1989
C     symmetrized tuh 13 April 90
C
#include <implicit.h>
#include <priunit.h>
#include <codata.h>
      PARAMETER (D0 = 0.D0, D1 = 1.D0, DP5 = 0.5D0)
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
      LOGICAL DOREL
      DIMENSION PMCAR(NCART), SCAL(NCART), AM(2), QCM(3,2), PCM(3,2),
     *          QCMDIF(3)
#include <cbiwlk.h>
#include <nuclei.h>
#include <symmet.h>
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('PMCAR in WLKTRA',-1)
         CALL OUTPUT(PMCAR,1,1,1,NCART,1,NCART,1,LUPRI)
      END IF
C
      CALL DZERO(AM,2)
      CALL DZERO(QCM,6)
      CALL DZERO(PCM,6)
      IATOM = 0
      DO 100 ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
         FAC    = D1/SQRT(FMULT(MULCNT))
         NATTYP = IZATOM(ICENT)
         NISOTP = ISOTOP(ICENT)
         AMASS  = XFAMU*DISOTP(NATTYP,NISOTP,'MASS')
         DO 200 ISYMOP = 0, MAXOPR
            IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
               IATOM = IATOM + 1
               IFRAG = IPART(IATOM)
               IF (IFRAG .EQ. 1 .OR. IFRAG .EQ. 2) THEN
                  AM(IFRAG) = AM(IFRAG) + AMASS
                  DO 300 ICOOR = 1, 3
                     ISCOOR = IPTCNT(3*(ICENT-1)+ICOOR,0,1)
                     IF (ISCOOR .NE. 0) THEN
                        SGNXYZ = PT(IBTAND(ISYMAX(ICOOR,1),ISYMOP))
                        QCM(ICOOR,IFRAG) = QCM(ICOOR,IFRAG)
     &                                  + SGNXYZ*AMASS*CORD(ICOOR,ICENT)
                        PCM(ICOOR,IFRAG) = PCM(ICOOR,IFRAG)
     &                                  + SGNXYZ*FAC*PMCAR(ISCOOR)
                     END IF
  300             CONTINUE
               END IF
            END IF
  200    CONTINUE
  100 CONTINUE
      DOREL = (AM(1) .GT. D0) .AND. (AM(2) .GT. D0)
      IF (DOREL) THEN
         DO 400 ICOOR = 1, 3
            QCMDIF(ICOOR) = QCM(ICOOR,1)/AM(1) - QCM(ICOOR,2)/AM(2)
  400    CONTINUE
         CMDIF = DNRM2(3,QCMDIF,1)
         CALL DSCAL(3,(D1/CMDIF),QCMDIF,1)
         P1 = DDOT (3,QCMDIF,1,PCM(1,1),1)
         P2 = DDOT (3,QCMDIF,1,PCM(1,2),1)
         VREL  = P1/AM(1) - P2/AM(2)
C
C        Kinetic energy release
C
         EKINPA = DP5*(AM(1)*AM(2)/(AM(1)+AM(2)))*VREL*VREL
      ELSE
         VREL   = D0
         EKINPA = D0
         CMDIF  = D0
      END IF
C
C     Print
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Output from WLKTRA',-1)
         WRITE (LUPRI,'(1X,A,1P,E24.12)') ' EKINPA ', EKINPA
         WRITE (LUPRI,'(1X,A,1P,E24.12)') ' VREL   ', VREL
      END IF
      RETURN
      END
C  /* Deck wlkpdi */
      SUBROUTINE WLKPDI(PMCARX,PMDIA,EVEC,SCAL,NCRTOT,NCART,IPRINT)
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION PMCARX(NCRTOT), PMDIA(NCRTOT), SCAL(NCRTOT),
     *          EVEC(NCART,NCART)
C
      IF (IPRINT .GT. 2) THEN
         CALL HEADER('Unscaled Cartesian momenta',-1)
         CALL OUTPUT(PMCARX,1,1,1,NCART,1,NCART,1,LUPRI)
         CALL HEADER('Scale factors',-1)
         CALL OUTPUT(SCAL,1,1,1,NCART,1,NCART,1,LUPRI)
      END IF
C
C     Momenta in diagonal representation
C
      CALL DZERO(PMDIA,NCART)
      DO 100 I = 1, NCART
         CALL DAXPY(NCART,SCAL(I)*PMCARX(I),EVEC(I,1),NCART,PMDIA,1)
  100 CONTINUE
      IF (IPRINT .GT. 2) THEN
         CALL HEADER('Momenta in diagonal representation',-1)
         CALL OUTPUT(PMDIA,1,1,1,NCART,1,NCART,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck wlkinc */
      SUBROUTINE WLKINC(STPSYM,STPCAR,COOR,DISPLC,NCART,NCRIND,IPOINT,
     *                  IPRINT)
C
C     Atomic displacements for numerical differentiation
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (D2 = 2.0D0)
C
      DIMENSION STPSYM(NCART), STPCAR(NCRIND), COOR(NCRIND)
C
#include <nuclei.h>
#include <symmet.h>
C
C     ***** STPSYM *****
C
      CALL DZERO(STPSYM,NCART)
      IF (IPOINT .EQ. 0) THEN
         STPSYM(1) = DISPLC
      ELSE IF (IPOINT .EQ. 2*NCART) THEN
         STPSYM(NCART) = DISPLC
      ELSE IF (MOD(IPOINT,2) .EQ. 1) THEN
         STPSYM((IPOINT+1)/2) = - D2*DISPLC
      ELSE
         STPSYM((IPOINT+1)/2) = DISPLC
         STPSYM((IPOINT+2)/2) = DISPLC
      END IF
C
C     ***** STPCAR *****
C
      CALL DZERO(STPCAR,NCRIND)
      DO 100 IATOM = 1, NUCIND
         DO 110 ICOOR = 1, 3
            ICCOOR = 3*(IATOM - 1) + ICOOR
            ISCOOR = IPTCNT(ICCOOR,0,1)
            IF (ISCOOR .GT. 0) THEN
              STPCAR(ICCOOR) = STPSYM(ISCOOR)/SQRT(FMULT(ISTBNU(IATOM)))
            END IF
 110     CONTINUE
 100  CONTINUE
C
C     ***** COOR *****
C
      IJ = 1
      DO 200 IATOM = 1, NUCIND
         DO 210 ICOOR = 1, 3
            COOR(IJ) = CORD(ICOOR,IATOM) + STPCAR(IJ)
            IJ = IJ + 1
 210     CONTINUE
 200  CONTINUE
C
C     ***** Print *****
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('STPSYM in WLKINC',1)
         CALL OUTPUT(STPSYM,1,1,1,NCART,1,NCART,1,LUPRI)
         CALL HEADER('STPCAR in WLKINC',1)
         CALL OUTPUT(STPCAR,1,1,1,NCRIND,1,NCRIND,1,LUPRI)
         CALL HEADER('COOR in WLKINC',1)
         CALL OUTPUT(COOR,1,1,1,NCRIND,1,NCRIND,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck wlknci */
      SUBROUTINE WLKNCI(STPCAR,COOR,DISPLC,NCART,NCRIND,EVEC,
     *                  IPOINT,NUMMOD,ACCURT,IPRINT)
C
C     Atomic displacements along normal modes for numerical differentiation
C
#include <implicit.h>
#include <priunit.h>
#include <codata.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (D2 = 2.0D0, D4 = 4.0D0)
      LOGICAL ACCURT
C
      DIMENSION STPCAR(NCRIND), COOR(NCRIND), EVEC(NCART,NCART)
C
#include <nuclei.h>
#include <symmet.h>
C
C     ***** STPCAR *****
C     Note: This only works without symmetry. If possible to include symmetry,
C           compare with structure in WLKINC. K.Ruud-March 1999
C
      DISPNC = DISPLC*SQRT(XFAMU)
      CALL DZERO(STPCAR,NCRIND)
      DO 100 IATOM = 1, NUCIND
         DO 110 ICOOR = 1, 3
            ICCOOR = 3*(IATOM - 1) + ICOOR
            IF (ACCURT) THEN
               IF (IPOINT .EQ. 0) THEN
                  STPCAR(ICCOOR) = EVEC(ICCOOR,1)*DISPNC
               ELSE IF (IPOINT .EQ. 4*NUMMOD) THEN
                  STPCAR(ICCOOR) = EVEC(ICCOOR,NUMMOD)*DISPNC
               ELSE IF ((MOD(IPOINT,4) .EQ. 1) .OR. 
     &                  (MOD(IPOINT,4) .EQ. 3)) THEN
                  STPCAR(ICCOOR) = EVEC(ICCOOR,(IPOINT+3)/4)*DISPNC
               ELSE IF (MOD(IPOINT,4) .EQ. 2) THEN
                  STPCAR(ICCOOR) = -D4*EVEC(ICCOOR,(IPOINT+3)/4)*DISPNC
               ELSE
                  STPCAR(ICCOOR) = EVEC(ICCOOR,(IPOINT + 3)/4)*DISPNC
     &                           + EVEC(ICCOOR,(IPOINT + 4)/4)*DISPNC
               END IF
            ELSE
               IF (IPOINT .EQ. 0) THEN
                  STPCAR(ICCOOR) = EVEC(ICCOOR,1)*DISPNC
               ELSE IF (IPOINT .EQ. 2*NUMMOD) THEN
                  STPCAR(ICCOOR) = EVEC(ICCOOR,NUMMOD)*DISPNC
               ELSE IF (MOD(IPOINT,2) .EQ. 1) THEN
                  STPCAR(ICCOOR) = -D2*EVEC(ICCOOR,(IPOINT+1)/2)*DISPNC
               ELSE
                  STPCAR(ICCOOR) = EVEC(ICCOOR,(IPOINT + 1)/2)*DISPNC
     &                           + EVEC(ICCOOR,(IPOINT + 2)/2)*DISPNC
               END IF
            END IF
 110     CONTINUE
 100  CONTINUE
C
C     ***** COOR *****
C
      IJ = 1
      DO 200 IATOM = 1, NUCIND
         DO 210 ICOOR = 1, 3
            COOR(IJ) = CORD(ICOOR,IATOM) + STPCAR(IJ)
            IJ = IJ + 1
 210     CONTINUE
 200  CONTINUE
C
C     ***** Print *****
C
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('STPCAR in WLKNCI',1)
         CALL OUTPUT(STPCAR,1,1,1,NCRIND,1,NCRIND,1,LUPRI)
         CALL HEADER('COOR in WLKNCI',1)
         CALL OUTPUT(COOR,1,1,1,NCRIND,1,NCRIND,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck wlknum */
      SUBROUTINE WLKNUM(GRDCAR,HESCAR,DIPG0,GRAD0,HESS0,GRADF,HESSF,
     *                  DIPGF,WORK,ERGMOL,ESTART,DIPM,DIPM0,DISPLC,
     *                  LWORK,DODIP,NXYZ,NCRTOT,N2CRT,NCART,IPOINT,
     *                  NUMHES,IPRINT)
C
C     Numerical differentiation
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0)
C
      LOGICAL FINAL, DODIP, NUMHES
      CHARACTER*9 TYPE
      DIMENSION GRDCAR(NCRTOT), HESCAR(N2CRT), GRAD0(NCART),
     *          HESS0(NCART,NCART), GRADF(NCART), HESSF(NCART,NCART),
     *          DIPG0(3,NCART), DIPGF(3,NCART), DIPM(3), DIPM0(3),
     *          WORK(LWORK)
C
#include <nuclei.h>
#include <symmet.h>
C
#include <chrsgn.h>
C
      FINAL = IPOINT .EQ. 2*NCART
      IF (IPOINT .EQ. 0) THEN
         CALL DCOPY (NCART,GRDCAR,1,GRAD0,1)
         CALL DCOPY (NCART*NCART,HESCAR,1,HESS0,1)
         CALL DZERO(GRADF,NCART)
         CALL DZERO(HESSF,NCART*NCART)
         IF (DODIP) THEN
            CALL DCOPY(3,DIPM,1,DIPM0,1)
            CALL DZERO(DIPGF,3*NCART)
         END IF
      ELSE
         I = (IPOINT+1)/2
         IF (MOD(IPOINT,2) .EQ. 1) THEN
            GRADF(I) = (ERGMOL - ESTART)/DISPLC
            DO 100 J = 1, NCART
               HESSF(J,I) = (GRDCAR(J)-GRAD0(J))/DISPLC
  100       CONTINUE
            DO 110 J = 1, NXYZ
               DIPGF(J,I) = (DIPM(J) - DIPM0(J))/DISPLC
  110       CONTINUE
         ELSE
            GRADF(I) = DP5*(GRADF(I) - (ERGMOL - ESTART)/DISPLC)
            DO 200 J = 1, NCART
               HESSF(J,I) = DP5*(HESSF(J,I)-(GRDCAR(J)-GRAD0(J))/DISPLC)
  200       CONTINUE
            DO 210 J = 1, NXYZ
               DIPGF(J,I) = DP5*(DIPGF(J,I)-(DIPM(J)-DIPM0(J))/DISPLC)
  210       CONTINUE
         END IF
      END IF
C
C     ***** Print *****
C
      IF (.NOT. NUMHES) THEN
      CALL TITLER('ABACUS - NUMERICAL DIFFERENTIATION','*',118)
      IF (MAXREP .GT.0) THEN
         WRITE (LUPRI,'(A,/)')
     *   ' Note: Only totally symmetric derivatives are tested.'
      END IF
      IF (IPOINT .EQ. 0) THEN
         WRITE (LUPRI,'(A,F17.10)')
     *         ' Energy at reference geometry:',ERGMOL
         IF (DODIP) WRITE (LUPRI,'(A,3F15.10)')
     *       ' Dipole moment at ref. geometry:',(DIPM0(I),I=1,NXYZ)
         CALL HEADER('Analytical gradient at reference geometry',1)
         CALL OUTPUT(GRDCAR,1,1,1,NCART,1,NCART,1,LUPRI)
         CALL HEADER('Analytical Hessian at reference geometry',1)
         CALL OUTPUT(HESCAR,1,NCART,1,NCART,NCART,NCART,1,LUPRI)
         IF (DODIP) THEN
            CALL HEADER('Analytical dipole gradient at ref. geometry',1)
            CALL OUTPUT(DIPG0,1,NXYZ,1,NCART,3,NCART,1,LUPRI)
         END IF
      ELSE
         ICUR = (IPOINT+1)/2
         WRITE (LUPRI,'(2A,/,A,10X,21X,A,F12.10)')
     *         ' Coordinate displaced in this calculation:    ',
     *         NAMEX(IPTCOR(ICUR,1)),
     *         ' Displacement:', CHRSGN(2*MOD(IPOINT,2)-1), DISPLC
         IF (FINAL) THEN
            WRITE (LUPRI,'(/A)')
     *       ' All displacements are now done'//
     *       ' - numerical differentiation is complete.'
         END IF
         CALL AROUND('Analytical and numerical gradients')
         WRITE (LUPRI,'(/,15X,A,/,1X,75A)')
     *         'analytical          numerical          difference',
     *          ('-',I=1,75)
         DIFMAX = D0
         DO 300 ICOOR = 1, ICUR
            IF (ICOOR .LT. ICUR) THEN
               TYPE = ' central '
            ELSE
               IF (MOD(IPOINT,2) .EQ. 1) THEN
                  TYPE = ' forward '
               ELSE
                  TYPE = ' central '
               END IF
            END IF
            DIFFER = GRAD0(ICOOR) - GRADF(ICOOR)
            DIFMAX = MAX(DIFMAX,ABS(DIFFER))
            WRITE (LUPRI,'(1X,A,1X,3F20.15,A)')
     *             NAMEX(IPTCOR(ICOOR,1)) , GRAD0(ICOOR), GRADF(ICOOR),
     *             DIFFER, TYPE
  300    CONTINUE
         IF (MOD(IPOINT,2) .EQ. 0) THEN
            WRITE (LUPRI,'(/,A,1P,D7.1//)') ' Largest difference: ',
     *             DIFMAX
         END IF
         CALL AROUND('Analytical and numerical Hessians')
         CALL HEADER('Analytical Hessian',1)
         CALL OUTPUT(HESS0,1,NCART,1,NCART,NCART,NCART,1,LUPRI)
C
         CALL HEADER('Numerical Hessian',1)
         CALL OUTPUT(HESSF,1,NCART,1,ICUR,NCART,NCART,1,LUPRI)
         IF (FINAL) THEN
            DIFMAX = D0
            DO 400 I = 1, NCART
               DO 410 J = 1, I
                  DIFFER = HESSF(J,I) - HESSF(I,J)
                  DIFMAX = MAX(DIFMAX,ABS(DIFFER))
  410          CONTINUE
  400       CONTINUE
            WRITE (LUPRI,'(/A,1P,D7.1,/)')
     *         '  Numerical Hessian is symmetric to ',DIFMAX
         END IF
C
         IJ = 0
         DIFMAX = D0
         DO 500 I = 1, ICUR
            DO 510 J = 1, NCART
               IJ = IJ + 1
               DIFFER = HESS0(J,I) - HESSF(J,I)
               DIFMAX = MAX(DIFMAX,ABS(DIFFER))
               WORK(IJ) = DIFFER
  510       CONTINUE
  500    CONTINUE
         IF (IJ .GT. LWORK) CALL STOPIT('WLKNUM',' ',IJ,LWORK)
         CALL HEADER
     *    ('Difference between analytical and numerical Hessians',1)
         CALL OUTPUT(WORK,1,NCART,1,ICUR,NCART,NCART,1,LUPRI)
         IF (MOD(IPOINT,2) .EQ. 0) THEN
            WRITE (LUPRI,'(/,A,1P,D7.1//)') '  Largest difference: ',
     *         DIFMAX
         END IF
         IF (DODIP) THEN
            CALL AROUND('Analytical and numerical dipole gradients')
            CALL HEADER('Analytical dipole gradient',1)
            CALL OUTPUT(DIPG0,1,NXYZ,1,NCART,3,NCART,1,LUPRI)
            CALL HEADER('Numerical dipole gradient',1)
            CALL OUTPUT(DIPGF,1,NXYZ,1,ICUR,3,NCART,1,LUPRI)
            IJ = 0
            DIFMAX = D0
            DO 600 I = 1, ICUR
               DO 610 J = 1, NXYZ
                  IJ = IJ + 1
                  DIFFER = DIPG0(J,I) - DIPGF(J,I)
                  DIFMAX = MAX(DIFMAX,ABS(DIFFER))
                  WORK(IJ) = DIFFER
  610          CONTINUE
  600       CONTINUE
            IF (IJ .GT. LWORK) CALL STOPIT('WLKNUM',' ',IJ,LWORK)
            CALL HEADER('Difference between analytical and numerical '//
     *                  'dipole gradients',1)
            CALL OUTPUT(WORK,1,NXYZ,1,ICUR,NXYZ,NCART,1,LUPRI)
            IF (MOD(IPOINT,2) .EQ. 0) THEN
               WRITE (LUPRI,'(/,A,1P,D7.1//)') '  Largest difference: ',
     *            DIFMAX
            END IF
         END IF
      END IF
      END IF
      RETURN
      END
C  /* Deck wlkv3 */
      SUBROUTINE WLKV3(GRDCAR,HESCAR,GRAD0,HESS0,GRADF,HESSF,
     *                 FMATF,WORK,ERGMOL,ESTART,DISPLC,
     *                 LWORK,NCRTOT,N2CRT,NCART,IPOINT,
     *                 IPRINT)
C
C     Determine anharmonic force field by numerical differentiation
C     Currently the entire anharmonic is determined only, using 
C     numerical differentiation of the Hessian. In the future V_imm should 
C     be determined in internal coordinates, using numerical differentiation
C     of the gradient. This can be used for speeding of the determination 
C     of the effective geometry needed for rovibrational averaging.
C     
C     K.Ruud, San Diego March 1999
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0)
C
      LOGICAL FINAL
      DIMENSION GRDCAR(NCRTOT), HESCAR(NCART,NCART), GRAD0(NCART),
     *          HESS0(NCART,NCART), GRADF(NCART), HESSF(NCART,NCART),
     *          FMATF(NCART,NCART,NCART), WORK(LWORK)
C
#include <nuclei.h>
#include <symmet.h>
C
#include <chrsgn.h>
C
      FINAL = IPOINT .EQ. 2*NCART
      IF (IPOINT .EQ. 0) THEN
         CALL DCOPY (NCART,GRDCAR,1,GRAD0,1)
         CALL DCOPY (NCART*NCART,HESCAR,1,HESS0,1)
         CALL DZERO(FMATF,NCART*NCART*NCART)
      ELSE
         I = (IPOINT+1)/2
         IF (MOD(IPOINT,2) .EQ. 1) THEN
            DO J = 1, NCART
               DO K = 1, NCART
                  FMATF(J,K,I) = (HESCAR(J,K)-HESS0(J,K))/DISPLC
               END DO
            END DO
         ELSE
            DO J = 1, NCART
               DO K = 1, NCART
                  FMATF(J,K,I) = DP5*(FMATF(J,K,I)
     &                         - (HESCAR(J,K)-HESS0(J,K))/DISPLC)
               END DO
            END DO
         END IF
      END IF
C
C     ***** Print *****
C
      CALL TITLER('ABACUS - NUMERICAL DIFFERENTIATION','*',118)
      IF (IPOINT .EQ. 0) THEN
         WRITE (LUPRI,'(A,F17.10)')
     *         ' Energy at reference geometry:',ERGMOL
         CALL HEADER('Analytical gradient at reference geometry',1)
         CALL OUTPUT(GRDCAR,1,1,1,NCART,1,NCART,1,LUPRI)
         CALL HEADER('Analytical Hessian at reference geometry',1)
         CALL OUTPUT(HESCAR,1,NCART,1,NCART,NCART,NCART,1,LUPRI)
      ELSE
         ICUR = (IPOINT+1)/2
         WRITE (LUPRI,'(2A,/,A,10X,21X,A,F12.10)')
     *         ' Coordinate displaced in this calculation:    ',
     *         NAMEX(IPTCOR(ICUR,1)),
     *         ' Displacement:', CHRSGN(2*MOD(IPOINT,2)-1), DISPLC
         IF (IPRINT .GE. 2) THEN
            CALL HEADER('Analytical gradient at displaced geometry',1)
            CALL OUTPUT(GRDCAR,1,1,1,NCART,1,NCART,1,LUPRI)
            CALL HEADER('Analytical Hessian at displaced geometry',1)
            CALL OUTPUT(HESCAR,1,NCART,1,NCART,NCART,NCART,1,LUPRI)
         END IF
         IF (FINAL) THEN
            WRITE (LUPRI,'(/A)')
     *       ' All displacements are now done'//
     *       ' - numerical differentiation is complete.'

            DO I = 1, NCART
               WRITE (LUPRI,'(/A,I5,/A)') 
     &              ' Anharmonic force constants (*,*,I):', I,
     &              ' ----------------------------------------'
               CALL OUTPUT(FMATF(1,1,I),1,NCART,1,NCART,NCART,NCART,
     &                     1,LUPRI)
            END DO
Ckr
Ckr      Print anharmonic force field and analyze in terms of normal
Ckr      coordinates. Need also some consistency check for the accuray of
Ckr      of the numerical differentiation, qua what is done for the symmetric
Ckr      nature of the Hessian.
Ckr            
         END IF
C
      END IF
      RETURN
      END
C  /* Deck wlkv3n */
      SUBROUTINE WLKV3N(GRDCAR,GRAD0,GRADIN,EVEC,
     &                 FMATF,WORK,ERGMOL,ESTART,DISPLC,
     &                 LWORK,NCRTOT,NCART,NUMMOD,IPOINT,IPRINT)
C
C     Determine anharmonic force field by numerical differentiation
C     along normal coordinates. Only the anharmonic force constants
C     needed for determining an effective geometry is considered.
C     
C     K.Ruud, San Diego March 1999
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0)
C
      LOGICAL FINAL
      DIMENSION GRDCAR(NCRTOT), GRAD0(NCART), GRADIN(NCART),
     *          FMATF(NCART,NCART,NCART), EVEC(NCART,NCART), WORK(LWORK)
C
#include <nuclei.h>
#include <symmet.h>
C
#include <chrsgn.h>
C
      FINAL = IPOINT .EQ. 2*NUMMOD
C      
      CALL DGEMM('N','N',1,NUMMOD,NCART,1.D0,
     &           GRDCAR,1,
     &           EVEC,NCART,0.D0,
     &           GRADIN,1)
      IF (IPOINT .EQ. 0) THEN
         CALL DCOPY (NUMMOD,GRADIN,1,GRAD0,1)
         CALL DZERO(FMATF,NCART*NCART*NCART)
      ELSE
         I = (IPOINT+1)/2
         IF (MOD(IPOINT,2) .EQ. 1) THEN
            DO J = 1, NUMMOD
               FMATF(I,I,J) = (GRADIN(J)-GRAD0(J))/DISPLC**2
            END DO
         ELSE
            DO J = 1, NUMMOD
               FMATF(I,I,J) = FMATF(I,I,J)
     &                      + (GRADIN(J)-GRAD0(J))/DISPLC**2
            END DO
         END IF
      END IF
C
C     ***** Print *****
C
      CALL TITLER('ABACUS - NUMERICAL DIFFERENTIATION','*',118)
      IF (IPOINT .EQ. 0) THEN
         WRITE (LUPRI,'(A,F17.10)')
     *         ' Energy at reference geometry:',ERGMOL
         CALL HEADER('Analytical gradient at reference geometry',1)
         CALL OUTPUT(GRDCAR,1,1,1,NCART,1,NCART,1,LUPRI)
         CALL HEADER('Analytical gradient at reference geometry'//
     &        '(normal coordinate basis)',1)
         CALL OUTPUT(GRADIN,1,1,1,NUMMOD,1,NCART,1,LUPRI)
      ELSE
         ICUR = (IPOINT+1)/2
         WRITE (LUPRI,'(2A,/,A,10X,21X,A,F12.10)')
     *         ' Coordinate displaced in this calculation:    ',
     *         NAMEX(IPTCOR(ICUR,1)),
     *         ' Displacement:', CHRSGN(2*MOD(IPOINT,2)-1), DISPLC
         IF (IPRINT .GE. 2) THEN
            CALL HEADER('Analytical gradient at displaced geometry'//
     &           '(normal coordinate basis)',1)
            CALL OUTPUT(GRADIN,1,1,1,NUMMOD,1,NCART,1,LUPRI)
         END IF
         IF (FINAL) THEN
            WRITE (LUPRI,'(/A)')
     *       ' All displacements are now done'//
     *       ' - numerical differentiation is complete.'

            DO I = 1, NUMMOD
               WRITE (LUPRI,'(/A,I5,/A)') 
     &              ' Anharmonic force constants (*,*,I):', I,
     &              ' ----------------------------------------'
               CALL OUTPUT(FMATF(1,1,I),1,NUMMOD,1,NUMMOD,NCART,NCART,
     &                     1,LUPRI)
            END DO
Ckr
Ckr      Print anharmonic force field and analyze in terms of normal
Ckr      coordinates. Need also some consistency check for the accuray of
Ckr      of the numerical differentiation, qua what is done for the symmetric
Ckr      nature of the Hessian.
Ckr            
         END IF
C
      END IF
      RETURN
      END
C  /* Deck wlkvib */
      SUBROUTINE WLKVIB(GRDCAR,HESCAR,GRAD0,HESS0,HESSF,EVAL,
     &                  DIPM0,DIPMF,SUSTO0,SUSTF,GFAC0,GFACF,QUAD0,
     &                  QUADF,QUADT,SIGMA0,SIGMAF,TMAT,CSTRA,SCTRA,SRC0,
     &                  SRCF,POLAR0,POLARF,POLDD,ALFA0,ALFAF,EFG0,EFGF,
     &                  SPNTOT,SSJ0,SSJF,ROANL,ROALO,ROANLF,ROALOF,
     &                  ROANL0,ROALO0,ERGMOL,ESTART,WORK,LWORK,NCART,
     &                  NCRTOT,NUMMOD,IPOINT,IPRINT)
C
C     Determine second-derivative of molecular properties along the normal
C     coordinates of the molecule at the effective molecular geometry.
C
C     K.Ruud, San Diego March 1999
C
#include <implicit.h>
#include <priunit.h>
#include <codata.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, DM1 = -1.0D0)
      PARAMETER (D3 = 3.0D0, D4  = 4.0D0, D12 = 12.0D0)
      PARAMETER (FACTOR = 1.D6*ALPHA2)
      PARAMETER (AUTOHZ = ALPHA2*ALPHA2/(4*XFAMU*XFAMU*PMASS*PMASS)
     &                   *6.5796838999D15)
C
      LOGICAL FINAL, ODD
#include <cbilnr.h>
      DIMENSION GRDCAR(NCRTOT), HESCAR(NCART,NCART), GRAD0(NCRTOT),
     &          HESS0(NCART,NCART), HESSF(NCART,NCART), DIPM0(3), 
     &          EVAL(NCART), DIPMF(3,NCART), SUSTO0(3,3), 
     &          SUSTF(3,3,NCART),
     &          GFAC0(3,3), GFACF(3,3,NCART), QUAD0(3,3),
     &          QUADF(3,3,NCART), QUADT(3,3), SIGMA0(3,3,MXCENT),
     &          SIGMAF(3,3,MXCENT,NCART), SRC0(3,3,MXCENT), 
     &          SRCF(3,3,MXCENT,NCART), TMAT(3,3,MXCENT), CSTRA(*), 
     &          SCTRA(*), POLAR0(3,3), POLARF(3,3,NCART), 
     &          POLDD(2,3,3,MXFR), ALFA0(3,3,MXFR), 
     &          ALFAF(3,3,MXFR,NCART),
     &          EFG0(3,3,MXCENT), EFGF(3,3,MXCENT,NCART), 
     &          SPNTOT(MXCOOR,MXCOOR), SSJ0(MXCOOR,MXCOOR), 
     &          SSJF(MXCOOR,MXCOOR,NCART), ROANL(2,3,3,MXFR), 
     &          ROALO(2,3,3,MXFR), ROANL0(3,3,MXFR), ROALO0(3,3,MXFR),
     &          ROANLF(3,3,MXFR,MXCOOR), ROALOF(3,3,MXFR,MXCOOR),
     &          WORK(LWORK)
C
#include <nuclei.h>
#include <symmet.h>
#include <abainf.h>
#include <dftcom.h>
C
#include <moldip.h>
#include <suscpt.h>
#include <molgfa.h>
#include <quadru.h>
#include <sigma.h>
#include <spinro.h>
#include <nqcc.h>
#include <cbiwlk.h>
C
#include <chrsgn.h>
C
      IF (MAGSUS) THEN
         DO I = 1, 3
            DO J = 1 ,3
               IF (NOLOND) THEN
                  SUSTOT(I,J) = SUSDIA(I,J) + SUSREL(I,J)
               ELSE
                  SUSTOT(I,J) = SUSREL(I,J) + SUSDIA(I,J)
     &                        + SUS2EL(I,J) + SUSFS (I,J)
     &                        + SUSFSY(I,J)
                 IF (DFTRUN) THEN
                  SUSTOT(I,J) = SUSTOT(I,J) + SUSDFT(I,J)
                 ENDIF
               END IF
            END DO
         END DO
      END IF
      IF (QUADRU) THEN
         CALL DCOPY(9,QDRNUC,1,QUADT,1)
         CALL DAXPY(9,DM1,QDREL,1,QUADT,1)
      END IF
      IF (SHIELD) THEN
         CALL TRADIP(SIGMAT,TMAT,CSTRA,SCTRA,3*NUCDEP,2,2)
         CALL DSCAL(9*MXCENT,FACTOR,TMAT,1)
      END IF
      IF (SPNSPN) CALL DSCAL(MXCOOR*MXCOOR,AUTOHZ,SPNTOT,1)
C
C     Transform properties to the Eckart frame
C
      IF (ECKART) THEN
         KGEOMI = 1
         KAMASS = KGEOMI + 3*MXCENT
         KNATTY = KAMASS + MXCENT
         KECKGE = KNATTY + MXCENT
         KGEOM2 = KECKGE + 3*MXCENT
         KLAST  = KGEOM2 + 3*MXCENT
         LWRK   = LWORK - KLAST + 1
         CALL ECKTRN(GRDCAR,HESCAR,DIP0,SUSTOT,GTRAN,QUADT,TMAT,GTRANT,
     &               POLARS,POLDD,ELNQC,SPNTOT,MAGSUS,MOLGFA,QUADRU,
     &               SHIELD,SPINRO,POLAR,ALFA,NQCC,SPNSPN,
     &               NCART,NCRTOT,MXFR,NFRVAL,FRVAL,WORK(KGEOMI),
     &               WORK(KAMASS),WORK(KNATTY),WORK(KECKGE),
     &               WORK(KGEOM2),IPRINT,.FALSE.,WORK(KLAST),LWRK)
      END IF
      IF (IPOINT .EQ. 0) THEN
         CALL DCOPY (NCART,GRDCAR,1,GRAD0,1)
         CALL DCOPY (NCART*NCART,HESCAR,1,HESS0,1)
         CALL DZERO(HESSF,NCART*NCART)
         CALL DCOPY(3,DIP0,1,DIPM0,1)
         IF (MAGSUS) THEN
            CALL DCOPY(9,SUSTOT,1,SUSTO0,1)
            CALL DSCAL(9,DM1,SUSTO0,1)
         END IF
         IF (MOLGFA) CALL DCOPY(9,GTRAN,1,GFAC0,1)
         IF (QUADRU) CALL DCOPY(9,QUADT,1,QUAD0,1)
         IF (SHIELD) CALL DCOPY(9*MXCENT,TMAT,1,SIGMA0,1)
         IF (SPINRO) CALL DCOPY(9*MXCENT,GTRANT,1,SRC0,1)
         IF (POLAR)  CALL DCOPY(9,POLARS,1,POLAR0,1)
         IF (ALFA)   CALL DCOPY(2*9*MXFR,POLDD,2,ALFA0,1)
         IF (NQCC)   CALL DCOPY(9*MXCENT,ELNQC,1,EFG0,1)
         IF (SPNSPN) CALL DCOPY(MXCOOR*MXCOOR,SPNTOT,1,SSJ0,1)
         IF (OPTROT) CALL DCOPY(2*9*MXFR,ROANL,2,ROANL0,1)
         IF (OPTROT) CALL DCOPY(2*9*MXFR,ROALO,2,ROALO0,1)
C
C     Should we transform the normal coordinates to an Eckart frame?
C
      ELSE
         IF (ACCURT) THEN
            FINAL = IPOINT .EQ. 4*NUMMOD
            I = (IPOINT+3)/4
         ELSE
            FINAL = IPOINT .EQ. 2*NUMMOD
            I = (IPOINT+1)/2
         END IF
         IF (ACCURT) THEN
            IWHERE = MOD(IPOINT,4)
            IF (IWHERE .EQ. 1) THEN
                  HESSF(I,1) = (D4/D3)*(ERGMOL - ESTART)/(DISPLC**2)
            ELSE IF (IWHERE .EQ. 2 .OR. IWHERE .EQ. 3) THEN
                  HESSF(I,1) = HESSF(I,1)
     &                 - (ERGMOL - ESTART)/(D12*(DISPLC**2))
            ELSE
                  HESSF(I,1) = HESSF(I,1) + 
     &                         (D4/D3)*(ERGMOL - ESTART)/(DISPLC**2)
            END IF
         ELSE
            IF (MOD(IPOINT,2) .EQ. 1) THEN
               HESSF(I,1) = (ERGMOL - ESTART)/DISPLC**2
               ODD = .TRUE.
            ELSE
               HESSF(I,1) = HESSF(I,1)
     &                    + (ERGMOL - ESTART)/DISPLC**2
               ODD = .FALSE.
            END IF
         END IF
         IF (.NOT. REUSED) THEN
            IF (ACCURT) THEN
               CALL FEMDIF(DIPM0,DIP0,DIPMF(1,I),DISPLC,IWHERE,3,1)
               IF (MAGSUS) THEN
                  CALL DSCAL(9,DM1,SUSTOT,1)
                  CALL FEMDIF(SUSTO0,SUSTOT,SUSTF(1,1,I),
     &                        DISPLC,IWHERE,9,1)
               END IF
               IF (MOLGFA) CALL FEMDIF(GFAC0,GTRAN,GFACF(1,1,I),
     &                              DISPLC,IWHERE,9,1)
               IF (QUADRU) CALL FEMDIF(QUAD0,QUADT,QUADF(1,1,I),
     &                              DISPLC,IWHERE,9,1)
               IF (SHIELD) CALL FEMDIF(SIGMA0,TMAT,SIGMAF(1,1,1,I),
     &                              DISPLC,IWHERE,9*NUCDEP,1)
               IF (SPINRO) CALL FEMDIF(SRC0,GTRANT,SRCF(1,1,1,I),
     &                              DISPLC,IWHERE,9*NUCDEP,1)
               IF (POLAR) CALL FEMDIF(POLAR0,POLARS,POLARF(1,1,I),
     &                             DISPLC,IWHERE,9,1)
               IF (ALFA) CALL FEMDIF(ALFA0,POLDD,ALFAF(1,1,1,I),
     &                            DISPLC,IWHERE,NFRVAL*9,2)
               IF (NQCC) CALL FEMDIF(EFG0,ELNQC,EFGF(1,1,1,I),
     &                            DISPLC,IWHERE,NUCDEP*9,1)
               IF (SPNSPN) CALL FEMDIF(SSJ0,SPNTOT,SSJF(1,1,I),
     &                              DISPLC,IWHERE,NUCDEP*NUCDEP*9,1)
               IF (OPTROT) THEN
                  CALL FEMDIF(ROANL0,ROANL,ROANLF(1,1,1,I),
     &                        DISPLC,IWHERE,NFRVAL*9,2)
                  CALL FEMDIF(ROALO0,ROALO,ROALOF(1,1,1,I),
     &                        DISPLC,IWHERE,NFRVAL*9,2)
               END IF
            ELSE
               CALL TREDIF(DIPM0,DIP0,DIPMF(1,I),DISPLC,ODD,3,1)
               IF (MAGSUS) THEN
                  CALL DSCAL(9,DM1,SUSTOT,1)
                  CALL TREDIF(SUSTO0,SUSTOT,SUSTF(1,1,I),
     &                        DISPLC,ODD,9,1)
               END IF
               IF (MOLGFA) CALL TREDIF(GFAC0,GTRAN,GFACF(1,1,I),
     &                                 DISPLC,ODD,9,1)
               IF (QUADRU) CALL TREDIF(QUAD0,QUADT,QUADF(1,1,I),
     &                                 DISPLC,ODD,9,1)
               IF (SHIELD) CALL TREDIF(SIGMA0,TMAT,SIGMAF(1,1,1,I),
     &                                 DISPLC,ODD,9*NUCDEP,1)
               IF (SPINRO) CALL TREDIF(SRC0,GTRANT,SRCF(1,1,1,I),
     &                                 DISPLC,ODD,9*NUCDEP,1)
               IF (POLAR) CALL TREDIF(POLAR0,POLARS,POLARF(1,1,I),
     &                                DISPLC,ODD,9,1)
               IF (ALFA) CALL TREDIF(ALFA0,POLDD,ALFAF(1,1,1,I),
     &                               DISPLC,ODD,NFRVAL*9,2)
               IF (NQCC) CALL TREDIF(EFG0,ELNQC,EFGF(1,1,1,I),
     &                               DISPLC,ODD,NUCDEP*9,1)
               IF (SPNSPN) CALL TREDIF(SSJ0,SPNTOT,SSJF(1,1,I),
     &                                 DISPLC,ODD,NUCDEP*NUCDEP*9,1)
               IF (OPTROT) THEN
                  CALL TREDIF(ROANL0,ROANL,ROANLF(1,1,1,I),
     &                        DISPLC,ODD,NFRVAL*9,2)
                  CALL TREDIF(ROALO0,ROALO,ROALOF(1,1,1,I),
     &                        DISPLC,ODD,NFRVAL*9,2)
               END IF
            END IF
         END IF

      END IF
C
C     ***** Print *****
C
      CALL TITLER('ABACUS - VIBRATIONAL AVERAGING','*',118)
      IF (IPOINT .EQ. 0) THEN
         WRITE (LUPRI,'(A,F17.10)')
     *         ' Energy at effective geometry:',ERGMOL
         CALL HEADER('Analytical gradient at effective geometry',1)
         CALL OUTPUT(GRDCAR,1,1,1,NCART,1,NCART,1,LUPRI)
         CALL HEADER('Analytical Hessian at effective geometry',1)
         CALL OUTPUT(HESCAR,1,NCART,1,NCART,NCART,NCART,1,LUPRI)
         CALL HEADER('Dipole moment at effective geometry',1)
         CALL OUTPUT(DIP0,1,1,1,3,1,3,1,LUPRI)
         IF (MAGSUS) THEN
            CALL HEADER('Magnetizabilities at effective geometry',1)
            CALL POLPRI(SUSTO0,'   ',-2)
         END IF
         IF (MOLGFA) THEN
            CALL HEADER('Rotational g tensor at effective geometry',1)
            CALL POLPRI(GFAC0,'   ',-2)
         END IF
         IF (QUADRU) THEN
            CALL HEADER('Molecular quadrupole moment at effective '//
     &                  'geometry',1)
            CALL POLPRI(QUAD0,'   ',-2)
         END IF
         IF (SHIELD) THEN
            CALL HEADER('Chemical shielding tensors at effective '//
     &                  'geometry',1)
            DO IATOM = 1, NUCDEP
               WRITE (LUPRI,'(18X,3(A,13X),/)') 'Bx', 'By', 'Bz'
               DO ICOOR = 1, 3
                  WRITE (LUPRI,'(2X,A,3F15.8)')
     &                 NAMDPX(3*(IATOM-1)+ICOOR), 
     &                 (SIGMA0(K,ICOOR,IATOM),K=1,3)
               END DO
            END DO
         END IF
         IF (SPINRO) THEN
            CALL HEADER('Nuclear spin-rotation tensors at effective '//
     &                  'geometry',1)
            DO IATOM = 1, NUCDEP
               WRITE (LUPRI,'(18X,3(A,13X),/)') 'Bx', 'By', 'Bz'
               DO ICOOR = 1, 3
                  WRITE (LUPRI,'(2X,A,3F15.8)')
     &                 NAMDPX(3*(IATOM-1)+ICOOR), 
     &                 (SRC0(K,ICOOR,IATOM),K=1,3)
               END DO
            END DO
         END IF
         IF (POLAR) THEN
            CALL HEADER('Static polarizability at effective geometry',1)
            CALL POLPRI(POLAR0,'   ',-1)
         END IF
         IF (ALFA) THEN
           CALL HEADER('Dynamic polarizability at effective geometry',1)
            DO IFRQR = 1, NFRVAL
               WRITE (LUPRI,'(/A,F8.5)') ' Frequency:   ',FRVAL(IFRQR)
               CALL POLPRI(ALFA0(1,1,IFRQR),'   ',-1)
            END DO
         END IF
         IF (NQCC) THEN
            CALL HEADER('Nuclear quadrupole moments at effective '//
     &                  'geometry',1)
            WRITE (LUPRI,'(18X,3(A,13X),/)') 'XX', 'YY', 'ZZ'
            DO IATOM = 1, NUCDEP
               DO ICOOR = 1, 3
                  WRITE (LUPRI,'(2X,A,3F15.8)')
     &                 NAMDPX(3*(IATOM-1)+ICOOR), 
     &                 (EFG0(K,ICOOR,IATOM),K=1,3)
               END DO
            END DO
         END IF
         IF (SPNSPN) THEN
            CALL HEADER('Spin-spin coupling constants at effective '//
     &                  'geometry',1)
            CALL PRIHES(SSJ0,'SPNSPN',CSTRA,SCTRA)
         END IF
         IF (OPTROT) THEN
            CALL HEADER('London orbital G matrix at effective '//
     &                  'geometry',1)
            DO IFRQR = 1, NFRVAL
               WRITE (LUPRI,'(/1X,A,F8.5)') 'Frequency:   ',FRVAL(IFRQR)
               WRITE(LUPRI,9009) 'Bx','By','Bz'
               WRITE(LUPRI,9008) 'Ex',(ROALO0(1,JDIP,IFRQR), JDIP=1,3)
               WRITE(LUPRI,9008) 'Ey',(ROALO0(2,JDIP,IFRQR), JDIP=1,3)
               WRITE(LUPRI,9008) 'Ez',(ROALO0(3,JDIP,IFRQR), JDIP=1,3)
            END DO
            CALL HEADER('No-London orbital G matrix at effective '//
     &                  'geometry',1)
            DO IFRQR = 1, NFRVAL
               WRITE (LUPRI,'(/1X,A,F8.5)') 'Frequency:   ',FRVAL(IFRQR)
               WRITE(LUPRI,9009) 'Bx','By','Bz'
               WRITE(LUPRI,9008) 'Ex',(ROANL0(1,JDIP,IFRQR), JDIP=1,3)
               WRITE(LUPRI,9008) 'Ey',(ROANL0(2,JDIP,IFRQR), JDIP=1,3)
               WRITE(LUPRI,9008) 'Ez',(ROANL0(3,JDIP,IFRQR), JDIP=1,3)
            END DO
         END IF
      ELSE
         ICUR = (IPOINT+1)/2
         WRITE (LUPRI,'(2A,/,A,10X,21X,A,F12.10)')
     *         ' Coordinate displaced in this calculation:    ',
     *         NAMEX(IPTCOR(ICUR,1)),
     *         ' Displacement:', CHRSGN(2*MOD(IPOINT,2)-1), DISPLC
         IF (IPRINT .GE. 4 .AND. .NOT. REUSED) THEN
            CALL HEADER('Dipole moment at displaced '//
     &                  'geometry',1)
            CALL OUTPUT(DIP0,1,1,1,3,1,3,1,LUPRI)
            IF (MAGSUS) THEN
               CALL HEADER('Magnetizabilities at displaced geometry',1)
               CALL POLPRI(SUSTOT,'   ',-2)            
            END IF
            IF (MOLGFA) THEN
               CALL HEADER('Rotational g tensor at displaced '//
     &                     'geometry',1)
               CALL POLPRI(GTRAN,'   ',-2)            
            END IF
            IF (QUADRU) THEN
               CALL HEADER('Molecular quadrupole moment at displaced '//
     &                     'geometry',1)
               CALL POLPRI(QUADT,'   ',-1)            
            END IF
            IF (SHIELD) THEN
               CALL HEADER('Chemical shielding tensors at displaced '//
     &                     'geometry',1)
               DO IATOM = 1, NUCDEP
                  WRITE (LUPRI,'(18X,3(A,13X),/)') 'Bx', 'By', 'Bz'
                  DO ICOOR = 1, 3
                     WRITE (LUPRI,'(2X,A,3F15.8)')
     &                    NAMDPX(3*(IATOM-1)+ICOOR), 
     &                    (TMAT(K,ICOOR,IATOM),K=1,3)
                  END DO
               END DO
            END IF
            IF (SPINRO) THEN
               CALL HEADER('Nuclear spin-rotation tensors at '//
     &                     'displaced geometry',1)
               DO IATOM = 1, NUCDEP
                  WRITE (LUPRI,'(18X,3(A,13X),/)') 'Bx', 'By', 'Bz'
                  DO ICOOR = 1, 3
                     WRITE (LUPRI,'(2X,A,3F15.8)')
     &                    NAMDPX(3*(IATOM-1)+ICOOR), 
     &                    (GTRANT(K,ICOOR,IATOM),K=1,3)
                  END DO
               END DO
            END IF
            IF (POLAR) THEN
               CALL HEADER('Polarizability at displaced geometry',1)
               CALL POLPRI(POLARS,'   ',-1)            
            END IF
            IF (ALFA) THEN
               CALL HEADER('Polarizability at displaced geometry',1)
               DO IFRQR = 1, NFRVAL
                  WRITE (LUPRI,'(/1X,A,F8.5)') 'Frequency:   ',
     &                 FRVAL(IFRQR)
                  CALL POLPRI(POLDD(1,1,1,IFRQR),'   ',-1)
               END DO
            END IF
            IF (NQCC) THEN
               CALL HEADER('Nuclear quadrupole moments at displaced '//
     &                     'geometry',1)
               WRITE (LUPRI,'(18X,3(A,13X),/)') 'XX', 'YY', 'ZZ'
               DO IATOM = 1, NUCDEP
                  DO ICOOR = 1, 3
                     WRITE (LUPRI,'(2X,A,3F15.8)')
     &                    NAMDPX(3*(IATOM-1)+ICOOR), 
     &                    (ELNQC(K,ICOOR,IATOM),K=1,3)
                  END DO
               END DO
            END IF
            IF (SPNSPN) THEN
               CALL HEADER('Spin-spin coupling constants at effective'//
     &                     ' geometry',1)
               CALL PRIHES(SPNTOT,'SPNSPN',CSTRA,SCTRA)
            END IF
            IF (OPTROT) THEN
               CALL HEADER('London orbital G matrix at effective '//
     &                     'geometry',1)
               DO IFRQR = 1, NFRVAL
               WRITE (LUPRI,'(/1X,A,F8.5)') 'Frequency:   ',FRVAL(IFRQR)
               WRITE(LUPRI,9009) 'Bx','By','Bz'
               WRITE(LUPRI,9008) 'Ex',(ROALO(1,1,JDIP,IFRQR), JDIP=1,3)
               WRITE(LUPRI,9008) 'Ey',(ROALO(1,2,JDIP,IFRQR), JDIP=1,3)
               WRITE(LUPRI,9008) 'Ez',(ROALO(1,3,JDIP,IFRQR), JDIP=1,3)
               END DO
               CALL HEADER('No-London orbital G matrix at effective '//
     &                  'geometry',1)
               DO IFRQR = 1, NFRVAL
               WRITE (LUPRI,'(/1X,A,F8.5)') 'Frequency:   ',FRVAL(IFRQR)
               WRITE(LUPRI,9009) 'Bx','By','Bz'
               WRITE(LUPRI,9008) 'Ex',(ROANL(1,1,JDIP,IFRQR), JDIP=1,3)
               WRITE(LUPRI,9008) 'Ey',(ROANL(1,2,JDIP,IFRQR), JDIP=1,3)
               WRITE(LUPRI,9008) 'Ez',(ROANL(1,3,JDIP,IFRQR), JDIP=1,3)
               END DO
            END IF
         END IF
         IF (FINAL) THEN
            WRITE (LUPRI,'(/A)')
     *       ' All displacements are now done'//
     *       ' - numerical differentiation is complete.'
            IF (REUSED) THEN
               WRITE (LUPRI,'(A)') ' New harmonic force '//
     &              'field, but derivatives read in from file'
C
C     We can either check what the user wants (by reading input file), 
C     since the user need to calculate properties in new effective geometry
C
               CALL ABAINP('**START',WORK,LWORK)
            END IF
            WRITE (LUPRI,'(/A,/A)') 
     &           ' Dipole moment second derivatives:',
     &           ' ---------------------------------'
            WRITE (LUPRI,'(5X,A4,5X,3(5X,A1,11X))') 'Mode','X',
     &           'Y','Z'
            DO IMOD = 1, NUMMOD
               WRITE (LUPRI,'(6X,I2,6X,F12.6,5X,F12.6,5X,F12.6)')
     &              IMOD,(DIPMF(ICOOR,IMOD),ICOOR = 1, 3)
            END DO
            IF (MAGSUS) THEN
               WRITE (LUPRI,'(/A,/A)') 
     &              ' Magnetizability tensor second derivatives:',
     &              ' ------------------------------------------'
               WRITE (LUPRI,'(/1X,A4,6(6X,A2,4X))') 'Mode','XX','XY',
     &              'YY','XZ','YZ','ZZ'
               DO IMOD = 1, NUMMOD
                  WRITE (LUPRI,'(2X,I2,1X,6(2X,F10.6))') IMOD,
     &               ((SUSTF(ICOOR,JCOOR,IMOD), ICOOR = 1, JCOOR),
     &                                          JCOOR = 1, 3)
               END DO
            END IF
            IF (MOLGFA) THEN
               WRITE (LUPRI,'(/A,/A)') 
     &              ' Rotational g tensor second derivatives:',
     &              ' ---------------------------------------'
               WRITE (LUPRI,'(/1X,A4,6(6X,A2,4X))') 'Mode','XX','XY',
     &              'YY','XZ','YZ','ZZ'
               DO IMOD = 1, NUMMOD
                  WRITE (LUPRI,'(2X,I2,1X,6(2X,F10.6))') IMOD,
     &               ((GFACF(ICOOR,JCOOR,IMOD), ICOOR = 1, JCOOR),
     &                                          JCOOR = 1, 3)
               END DO
            END IF
            IF (QUADRU) THEN
               WRITE (LUPRI,'(/A,/A)') 
     &              ' Molecular quadrupole moment second derivatives:',
     &              ' -----------------------------------------------'
               WRITE (LUPRI,'(/1X,A4,6(6X,A2,4X))') 'Mode','XX','XY',
     &              'YY','XZ','YZ','ZZ'
               DO IMOD = 1, NUMMOD
                  WRITE (LUPRI,'(2X,I2,1X,6(2X,F10.6))') IMOD,
     &               ((QUADF(ICOOR,JCOOR,IMOD), ICOOR = 1, JCOOR),
     &                                          JCOOR = 1, 3)
               END DO
            END IF
            IF (SHIELD) THEN
               DO IATOM = 1, NUCDEP
                  WRITE (LUPRI,'(/A,A,/A)')
     &             ' Nuclear shielding tensor second derivatives for: ',
     &              NAMDEP(IATOM),' --------------------------------'//
     &                 '------------------------'
                  WRITE (LUPRI,'(/1X,A4,1X,9(2X,A4,2X))') 'Mode','Bxmx',
     &              'Bxmy','Bxmz','Bymx','Bymy','Bymz','Bzmx','Bzmy',
     &              'Bzmz'
                  DO IMOD = 1, NUMMOD
                     WRITE (LUPRI,'(2X,I2,2X,9F8.3)') IMOD,
     &                    ((SIGMAF(K,I,IATOM,IMOD),I=1,3),K=1,3)
                  END DO
               END DO
            END IF
            IF (SPINRO) THEN
               DO IATOM = 1, NUCDEP
                  WRITE (LUPRI,'(/A,A,/A)')
     &             ' Nuclear spin-rotation tensor second derivatives '//
     &                 'for: ', NAMDEP(IATOM),
     &                 ' --------------------------------'//
     &                 '----------------------------'
                  WRITE (LUPRI,'(/1X,A4,1X,9(2X,A4,2X))') 'Mode','Jxmx',
     &              'Jxmy','Jxmz','Jymx','Jymy','Jymz','Jzmx','Jzmy',
     &              'Jzmz'
                  DO IMOD = 1, NUMMOD
                     WRITE (LUPRI,'(2X,I2,2X,9F8.3)') IMOD,
     &                    ((SRCF(K,I,IATOM,IMOD),I=1,3),K=1,3)
                  END DO
               END DO
            END IF
            IF (POLAR) THEN
               WRITE (LUPRI,'(/A,/A)') 
     &              ' Polarizability second derivatives:',
     &              ' ----------------------------------'
               WRITE (LUPRI,'(/1X,A4,6(6X,A2,4X))') 'Mode','XX','XY',
     &              'YY','XZ','YZ','ZZ'
               DO IMOD = 1, NUMMOD
                  WRITE (LUPRI,'(2X,I2,1X,6(2X,F10.6))') IMOD,
     &               ((POLARF(ICOOR,JCOOR,IMOD), ICOOR = 1, JCOOR),
     &                                           JCOOR = 1, 3)
               END DO
            END IF
            IF (ALFA) THEN
               DO IFRQR = 1, NFRVAL
                  WRITE (LUPRI,'(/A,2X,F8.5,/A)') 
     &              ' Polarizability second derivatives for frequency:',
     &              FRVAL(IFRQR),
     &              ' -----------------------------------------------'//
     &              '-----------'
                  WRITE (LUPRI,'(/1X,A4,6(6X,A2,4X))') 'Mode','XX','XY',
     &                 'YY','XZ','YZ','ZZ'
                  DO IMOD = 1, NUMMOD
                     WRITE (LUPRI,'(2X,I2,1X,6(2X,F10.6))') IMOD,
     &                    ((ALFAF(ICOOR,JCOOR,IFRQR,IMOD), 
     &                    ICOOR = 1, JCOOR), JCOOR = 1, 3)
                  END DO
               END DO
            END IF
            IF (NQCC) THEN
               DO IATOM = 1, NUCDEP
                  WRITE (LUPRI,'(/A,A,/A)')
     &             ' Nuclear quadrupole moment tensor second '//
     &                 'derivatives for: ', NAMDEP(IATOM),
     &                 ' --------------------------------'//
     &                 '----------------------------'
                  WRITE (LUPRI,'(/1X,A4,6(6X,A2,4X))') 'Mode','XX','XY',
     &                 'YY','XZ','YZ','ZZ'
                  DO IMOD = 1, NUMMOD
                     WRITE (LUPRI,'(2X,I2,1X,6(2X,F10.6))') IMOD,
     &               ((EFGF(ICOOR,JCOOR,IATOM,IMOD), ICOOR = 1, JCOOR),
     &                                               JCOOR = 1, 3)
                  END DO
               END DO
            END IF
            IF (SPNSPN) THEN
               DO IMOD = 1, NUMMOD
                  WRITE (LUPRI,'(/A,I4,/A)') 'Spin-spin coupling '//
     &                 'constant second derivatives for mode',IMOD,
     &                 '----------------------------------------'//
     &                 '------------------'
                  CALL PRIHES(SSJF(1,1,IMOD),'SPNSPN',CSTRA,SCTRA)
               END DO
            END IF
            IF (OPTROT) THEN
               DO IFRQR = 1, NFRVAL
               WRITE (LUPRI,'(/A,F8.5,/A)') 
     &           ' London G matrix second derivatives for frequency:',
     &                 FRVAL(IFRQR),
     &           ' -------------------------------------------------'//
     &           '--------'
               WRITE (LUPRI,'(/1X,A4,6(6X,A2,4X))') 'Mode','XX','XY',
     &              'YY','XZ','YZ','ZZ'
               DO IMOD = 1, NUMMOD
                  WRITE (LUPRI,'(2X,I2,1X,6(2X,F10.6))') IMOD,
     &              ((ROALOF(ICOOR,JCOOR,IFRQR,IMOD), ICOOR = 1, JCOOR),
     &                                           JCOOR = 1, 3)
               END DO
               WRITE (LUPRI,'(/A,F8.5,/A)') 
     &          ' No-London G matrix second derivatives for frequency:',
     &                 FRVAL(IFRQR),
     &          ' --------------------------------------------------'//
     &          '----------'
               WRITE (LUPRI,'(/1X,A4,6(6X,A2,4X))') 'Mode','XX','XY',
     &              'YY','XZ','YZ','ZZ'
               DO IMOD = 1, NUMMOD
                  WRITE (LUPRI,'(2X,I2,1X,6(2X,F10.6))') IMOD,
     &              ((ROANLF(ICOOR,JCOOR,IFRQR,IMOD), ICOOR = 1, JCOOR),
     &                                           JCOOR = 1, 3)
               END DO
               END DO
            END IF
C
C     Calculate rovibrational correction to properties.
C     Note: This relies on the calculation being done at the effective
C           geometry.
C
            KSUSAV = 1
            KGFACA = KSUSAV + 9
            KSRCA  = KGFACA + 9
            KPOLAV = KSRCA  + 9*MXCENT
            KALFAV = KPOLAV + 9
            KEFGAV = KALFAV + 9*MXFR
            KSSJAV = KEFGAV + 9*MXCENT
            KROANA = KSSJAV + MXCOOR*MXCOOR
            KROALA = KROANA + 9*MXFR
            KLAST  = KROALA + 9*MXFR
C
C           We need to use the appropriate frequencies for this geometry
C
            CALL DCOPY(NUMMOD,HESSF(1,1),1,EVAL,1)
            CALL DSCAL(NUMMOD,1.0D0/XFAMU,EVAL,1)
            CALL PRPVAV(EVAL,NUMMOD,DIPM0,DIPMF,SUSTO0,SUSTF,
     &                  WORK(KSUSAV),GFAC0,GFACF,WORK(KGFACA),
     &                  QUAD0,QUADF,QUADT,SIGMA0,SIGMAF,TMAT,SRC0,SRCF,
     &                  WORK(KSRCA),POLAR0,POLARF,WORK(KPOLAV),ALFA0,
     &                  ALFAF,WORK(KALFAV),EFG0,EFGF,WORK(KEFGAV),
     &                  SSJ0,SSJF,WORK(KSSJAV),ROALO0,ROALOF,
     &                  WORK(KROALA),ROANL0,ROANLF,WORK(KROANA),
     &                  CSTRA,SCTRA)
         END IF
C
      END IF
      RETURN
9008  FORMAT (16X,A,2X,3F12.7)
9009  FORMAT (/,27X,3(A,10X)/)
      END
C  /* Deck wlkgeo */
      SUBROUTINE WLKGEO(RADIUS,GRDDIA,EVAL,EVEC,STPDIA,STPSYM,STPCAR,
     &                  COOR,PMDIAX,PMDIA,PMCAR,SCAL,XMXNUC,TIMSTP,RNU,
     &                  ERGSEC,SGNMOD,ANHRAD,ZERGRD,TRSDIA,
     &                  NCART,NINTER,NCRIND,NNEG,IMODE,IWKIND,KEEPSY,
     &                  STATPO,WFPRED,START,INDXOK,IPRINT,IWKTYP,ISTTYP,
     &                  IMAGE,DONWTN,TRUSTD)
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (MAXIT = 50, ANGMAX = 30.0D0, D180 = 180.0D0,
     &           REDANG = 0.8D0)
C
      LOGICAL ADPOK, ANHOK, MCOK, START, INDXOK, WFPRED, KEEPSY, STATPO,
     &        IMAGE, DONWTN, ANGOK
      DIMENSION GRDDIA(NCART), EVAL(NCART),  EVEC(NCART,NCART),
     &          STPDIA(NCART), STPSYM(NCART), STPCAR(NCRIND),
     &          PMDIA(NCART),  PMCAR(NCART),   SCAL(NCART),
     &          COOR(NCRIND), TRSDIA(NCART), PMDIAX(NCART)
C
#include <nuclei.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from WLKGEO','*',103)
      IF (IPRINT .GT. 10) THEN
         WRITE (LUPRI,'(A,I5)')       ' IWKTYP ', IWKTYP
         WRITE (LUPRI,'(A,I5)')       ' IWKIND ', IWKIND
         WRITE (LUPRI,'(A,I5)')       ' IMODE  ', IMODE
         WRITE (LUPRI,'(A,I5)')       ' NINTER ', NINTER
         WRITE (LUPRI,'(A,I5)')       ' NCART  ', NCART
         WRITE (LUPRI,'(A,I5)')       ' NCRIND ', NCRIND
         WRITE (LUPRI,'(A,I5)')       ' NNEG   ', NNEG
         WRITE (LUPRI,'(A,L5)')       ' IMAGE  ', IMAGE
         WRITE (LUPRI,'(A,L5)')       ' START  ', START
         WRITE (LUPRI,'(A,L5)')       ' INDXOK ', INDXOK
         WRITE (LUPRI,'(A,L5)')       ' KEEPSY ', KEEPSY
         WRITE (LUPRI,'(A,L5)')       ' WFPRED ', WFPRED
         WRITE (LUPRI,'(A,1P,D12.6)') ' ANHRAD ', ANHRAD
         WRITE (LUPRI,'(A,1P,D12.6)') ' ZERGRD ', ZERGRD
         WRITE (LUPRI,'(A,1P,D12.6)') ' RADIUS ', RADIUS
         WRITE (LUPRI,'(A,1P,D12.6)') ' XMXNUC ', XMXNUC
         WRITE (LUPRI,'(A,1P,D12.6)') ' TIMSTP ', TIMSTP
         CALL HEADER('Gradient in WLKGEO',-1)
         CALL OUTPUT(GRDDIA,1,1,1,NINTER,1,NCART,1,LUPRI)
         CALL HEADER('Eigenvalues in WLKGEO',-1)
         CALL OUTPUT(EVAL,1,1,1,NINTER,1,NCART,1,LUPRI)
         CALL HEADER('Eigenvectors in WLKGEO',-1)
         CALL OUTPUT(EVEC,1,NCART,1,NINTER,NCART,NCART,1,LUPRI)
         CALL HEADER('Scale factors in WLKGEO',-1)
         CALL OUTPUT(SCAL,1,1,1,NCART,1,NCART,1,LUPRI)
         CALL HEADER('TRSDIA in WLKGEO',-1)
         CALL OUTPUT(TRSDIA,1,1,1,NCART,1,NCART,1,LUPRI)
      END IF
C
C     ***********************
C     ***** Newton step *****
C     ***********************
C
      IF (IWKTYP .EQ. 4) THEN
         ISTTYP = 4
         CALL WLKSTP(GRDDIA,EVAL,EVEC,STPDIA,STPSYM,STPCAR,PMDIA,PMCAR,
     *               SCAL,RADIUS,ERGSEC,RNU,STPNRM,STEPU,TIMSTP,SGNMOD,
     *               ANHRAD,ZERGRD,NCART,NINTER,NCRIND,NNEG,IMODE,
     &               IWKIND,KEEPSY,IPRINT,ISTTYP,IMAGE,DONWTN)
C        IF (WFPRED) THEN
C           CALL WLKMC(REDMC,STPCAR,WRK,LWRK,NCART,MCOK,.FALSE.,IPRINT)
C        END IF
         RADIUS = STPNRM
      ELSE
C
C     *****************************
C     ***** Trust region step *****
C     *****************************
C
C        Decrease trust radius until acceptable step has been found
C
         DO 100 ITER = 1, MAXIT
            IF (IWKTYP .EQ. 1) THEN
               IF (INDXOK .AND. ITER.EQ.1) THEN
                  ISTTYP = 4
               ELSE IF (START .AND. STATPO .AND. (IWKIND.GT.0)
     &                  .AND. .NOT.IMAGE) THEN
                  ISTTYP = 5
               ELSE
                  ISTTYP = IWKTYP
               END IF
            ELSE IF (IWKTYP .EQ. 2) THEN
               IF (INDXOK .AND. ITER.EQ.1) THEN
                  ISTTYP = 4
               ELSE IF (START .AND. STATPO) THEN
                  ISTTYP = 5
               ELSE
                  ISTTYP = IWKTYP
               END IF
            ELSE IF (IWKTYP .EQ. 3) THEN
               ISTTYP = 3
               CALL DCOPY(NCART,PMDIAX,1,PMDIA,1)
            ELSE IF (IWKTYP .EQ. 7) THEN
               IF (START) THEN
                  ISTTYP = 5
               ELSE
                  ISTTYP = 6
               END IF
            ELSE
               ISTTYP = IWKTYP
            END IF
            IF (IPRINT .GT. 2) THEN
               WRITE (LUPRI,'(/A,I3)') ' Iteration in WLKGEO:',ITER
               WRITE (LUPRI,'(A,F15.8)') ' Trust radius:',RADIUS
               WRITE (LUPRI,'(A,I3)')   ' ISTTYP ', ISTTYP
            END IF
C
C           ***** Calculate Cartesian step vector *****
C
            CALL WLKSTP(GRDDIA,EVAL,EVEC,STPDIA,STPSYM,STPCAR,PMDIA,
     &                  PMCAR,SCAL,RADIUS,ERGSEC,RNU,STPNRM,STEPU,
     &                  TIMSTP,SGNMOD,ANHRAD,ZERGRD,NCART,NINTER,NCRIND,
     &                  NNEG,IMODE,IWKIND,KEEPSY,IPRINT,ISTTYP,IMAGE,
     &                  DONWTN)
C
C           ***** Test size of Newton step *****
C
            IF (IWKTYP.LT.3 .AND. ISTTYP.EQ.4) THEN
               IF (STPNRM .GT. RADIUS) GO TO 100
            END IF
C
C           ***** Test atomic displacements *****
C
            IF (IWKTYP .EQ. 3 .OR. IWKTYP .EQ. 7) THEN
               ADPOK = .TRUE.
               CALL WLKCAN(REDANH,STPDIA,TRSDIA,NCART,ANHOK,IPRINT)
            ELSE
               ANHOK = .TRUE.
               CALL WLKADP(REDADP,XMXNUC,STPCAR,NCRIND,ADPOK,IPRINT)
            END IF
C
C           ***** Test angle between step and gradient *****
C
            IF (IWKTYP .EQ. 7) THEN
               IF (START) THEN
                  ANGOK = .TRUE.
               ELSE
                  ANGLE = WINKEL(STPDIA,GRDDIA,NCART)
                  ANGOK = ABS(ANGLE-D180) .LT. ANGMAX
               END IF
            ELSE
               ANGOK = .TRUE.
            END IF
C
C           ***** Calculate wave-function prediction and test *****
C
            IF (.NOT.ADPOK) THEN
               WRITE (LUPRI,'(/A)') '@ Step rejected since atomic'//
     *                              ' displacements are too large.'
               RADIUS = REDADP*MIN(STPNRM,RADIUS)
            ELSE IF (.NOT.ANHOK) THEN
               WRITE (LUPRI,'(/A)') '@ Step rejected since'//
     *                              ' anharmonic term is too large.'
               RADIUS = REDANH*MIN(STPNRM,RADIUS)
            ELSE IF ((IWKTYP .EQ. 7) .AND. (.NOT.ANGOK)) THEN
               WRITE (LUPRI,'(/A)') '@ Step rejected since '//
     *                              'step-gradient angle is too large.'
               WRITE (LUPRI,'(A,F12.6)') '@ Angle: ',ANGLE
               RADIUS = REDANG*MIN(STPNRM,RADIUS)
            ELSE
C              IF (WFPRED) THEN
C                 CALL WLKMC(REDMC,STPCAR,WRK,LWRK,NCART,
C    *                       MCOK,.TRUE.,IPRINT)
C              ELSE
                  MCOK = .TRUE.
C              END IF
               IF (MCOK) GO TO 1000
               RADIUS = TRUSTD*MIN(STPNRM,RADIUS)
               IF (IPRINT .GT. 0) THEN
                  WRITE (LUPRI,'(/,A)') '@ Step rejected since '//
     *                              ' orbital rotations are too large.'
               END IF
            END IF
  100    CONTINUE
C
C        A) No step has been found
C        -------------------------
C
         WRITE (LUPRI,'(//A,/A)') ' Acceptable step has not been'//
     *                     ' found in maximum number of iterations.',
     *                     ' Program cannot proceed.'
         CALL QUIT('Maximum iterations exceeded in WLKGEO.')
C
C        B) Step has been found
C        ----------------------
C
 1000    CONTINUE
         RADIUS = MIN(STPNRM,RADIUS)
         IF (IPRINT .GT. 2) THEN
            WRITE (LUPRI,'(/,A)') ' Acceptable step has been obtained.'
         END IF
      END IF
C
C     ************************
C     ***** New geometry *****
C     ************************
C
      IJ = 1
      DO 200 J = 1, NUCIND
         DO 210 I = 1, 3
            COOR(IJ) = CORD(I,J) + STPCAR(IJ)
            IJ = IJ + 1
 210     CONTINUE
 200  CONTINUE
      RETURN
      END
C  /* Deck wlkstp */
      SUBROUTINE WLKSTP(GRDDIA,EVAL,EVEC,STPDIA,STPSYM,STPCAR,PMDIA,
     &                  PMCAR,SCAL,RADIUS,ERGSEC,RNU,STPNRM,STEPU,
     &                  TIMSTP,SGNMOD,ANHRAD,ZERGRD,
     &                  NCART,NINTER,NCRIND,NNEG,
     &                  IMODE,IWKIND, KEEPSY,IPRINT,ISTTYP,
     &                  IMAGE,DONWTN)
C
C     The following kinds of steps are possible
C
C     1 - Level-shifted Newton step to boundary
C     2 - Gradient-extremal step to boundary
C     3 - Dynamic step to boundary
C     4 - Newton step
C     5 - Hessian eigenvector step to boundary
C     6 - Intrinsic reaction path
C
#include <implicit.h>
#include <priunit.h>
#include <codata.h>
      PARAMETER (D0 = 0.0 D0, DP5 = 0.5 D0)
C
      LOGICAL KEEPSY, IMAGE, DONWTN
      DIMENSION STPDIA(NCART), STPSYM(NCART), STPCAR(NCRIND),
     *          SCAL(NCART), EVAL(NCART), EVEC(NCART,NCART),
     *          GRDDIA(NCART), PMDIA(NCART), PMCAR(NCART)
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from WLKSTP','*',103)
      IF (ISTTYP .EQ. 3) THEN
         NVEC = NCART
      ELSE
         NVEC = NINTER
      END IF
      IF (IPRINT .GT. 10) THEN
         WRITE (LUPRI,'(A,I5)')       ' ISTTYP ', ISTTYP
         WRITE (LUPRI,'(A,I5)')       ' IWKIND ', IWKIND
         WRITE (LUPRI,'(A,I5)')       ' IMODE  ', IMODE
         WRITE (LUPRI,'(A,I5)')       ' NINTER ', NINTER
         WRITE (LUPRI,'(A,I5)')       ' NCART  ', NCART
         WRITE (LUPRI,'(A,I5)')       ' NNEG   ', NNEG
         WRITE (LUPRI,'(A,L5)')       ' IMAGE  ', IMAGE
         WRITE (LUPRI,'(A,L5)')       ' KEEPSY ', KEEPSY
         WRITE (LUPRI,'(A,1P,D12.6)') ' RADIUS ', RADIUS
         WRITE (LUPRI,'(A,1P,D12.6)') ' TIMSTP ', TIMSTP
         CALL HEADER('Gradient in WLKSTP',-1)
         CALL OUTPUT(GRDDIA,1,1,1,NVEC,1,NCART,1,LUPRI)
         CALL HEADER('Eigenvalues in WLKSTP',-1)
         CALL OUTPUT(EVAL,1,1,1,NVEC,1,NCART,1,LUPRI)
         CALL HEADER('Eigenvectors in WLKSTP',-1)
         CALL OUTPUT(EVEC,1,NCART,1,NVEC,NCART,NCART,1,LUPRI)
         CALL HEADER('Scale factors in WLKSTP',-1)
         CALL OUTPUT(SCAL,1,1,1,NCART,1,NCART,1,LUPRI)
      END IF
C
C     *******************************************
C     ***** Step in diagonal representation *****
C     *******************************************
C
      CALL DZERO(STPDIA,NCART)
      RNU = D0
C
      GO TO (1,2,3,4,5,6) ISTTYP
C     ========================
C
C     1) Level-shifted Newton step
C     ----------------------------
C
    1 CONTINUE
         IF (IWKIND .EQ. 0 .OR. IMAGE) THEN
           IF (IMAGE) THEN
              CALL WLKNEG(EVAL,GRDDIA,STPDIA,PMDIA,NINTER,IMODE,IPRINT)
           END IF
           CALL WLKFL0(NCART,NINTER,EVAL,GRDDIA,STPDIA,
     &                 RADIUS,RNU,KEEPSY,ZERGRD,IPRINT)
           IF (IMAGE) THEN
              CALL WLKNEG(EVAL,GRDDIA,STPDIA,PMDIA,NINTER,IMODE,IPRINT)
           END IF
         ELSE
           IF (IPRINT .GT. 2) WRITE (LUPRI,'(//1X,A,I5,A/)')
     *        'Level-shift step is taken following eigenmode ',IMODE,'.'
           CALL WLKFL1(NNEG,IMODE,NCART,NINTER,EVAL,GRDDIA,STPDIA,
     *                 RADIUS,RNU,IPRINT)
         END IF
      GO TO 1000
C
C     2) Gradient-extremal step
C     -------------------------
C
    2 CONTINUE
         IF (IPRINT .GT. 2) WRITE (LUPRI,'(//1X,A,I5,A/)')
     &      'Gradient-extremal step is taken following eigenmode ',
     &      IMODE, '.'
         CALL WLKEXT(IMODE,NCART,NINTER,EVAL,GRDDIA,STPDIA,
     &               RADIUS,KEEPSY,ZERGRD,ANHRAD,IWKIND,IPRINT)
      GO TO 1000
C
C     3) Dynamic step
C     ---------------
C
    3 CONTINUE
         IF (IPRINT .GT. 2) WRITE (LUPRI,'(//1X,A,/)') 'Dynamic step.'
         CALL WLKDYN(TIMSTP,STPDIA,PMDIA,GRDDIA,EVAL,NCART,RADIUS,
     &               IPRINT)
      GO TO 1000
C
C     4) Newton step
C     --------------
C
    4 CONTINUE
         IF (IPRINT.GT.2) WRITE (LUPRI,'(//1X,A,/)')'Newton step taken.'
         DO 100 I = 1, NINTER
            STPDIA(I) = - GRDDIA(I)/EVAL(I)
  100    CONTINUE
      GO TO 1000
C
C     5) Hessian eigenvector step
C     --------------------------
C
    5 CONTINUE
         IF (IPRINT .GT. 0) THEN
            IF (IMODE .GT. 0) THEN
               WRITE (LUPRI,'(//1X,A,I5,A/)')' A step is taken '//
     *            'in the POSITIVE direction along mode ',ABS(IMODE),'.'
            ELSE
               WRITE (LUPRI,'(//1X,A,I5,A/)')' A step is taken '//
     *            'in the NEGATIVE direction along mode ',ABS(IMODE),'.'
            END IF
         END IF
         STPDIA(IMODE) = SGNMOD*RADIUS
      GO TO 1000
C
C     6) Intrinsic reaction path step
C     -------------------------------
C
    6 CONTINUE
         IF (IPRINT .GT. 2) WRITE (LUPRI,'(//1X,A,/)') 'IRC step'
         CALL WLKIRC(STPDIA,GRDDIA,EVAL,NCART,RADIUS,DONWTN,IPRINT)
      GO TO 1000
C
 1000 CONTINUE
C     ========
C
C     *****************************
C     ***** Energy prediction *****
C     *****************************
C
      ERGSEC = DDOT(NINTER,GRDDIA,1,STPDIA,1)
     *       + DP5*DV3DOT(NINTER,STPDIA,EVAL,STPDIA)
C
C     ***** Print *****
C
      IF (IPRINT .GT. 2) THEN
         WRITE (LUPRI,'(/A,F25.15)') ' Predicted energy change', ERGSEC
         CALL HEADER('Step in diagonal representation'//
     *               ' (eigenvalue order)',1)
         CALL OUTPUT(STPDIA,1,1,1,NVEC,1,NCART,1,LUPRI)
         IF (IPRINT .GT. 6) THEN
            CALL HEADER('Eigenvector basis:',1)
            CALL OUTPUT(EVEC,1,NCART,1,NVEC,NCART,NCART,1,LUPRI)
         END IF
      END IF
C
C     ****************************************
C     ***** Scaled Cartesian step vector *****
C     ****************************************
C
      CALL DZERO(STPSYM,NCART)
      DO 200 I = 1, NVEC
         CALL DAXPY(NCART,STPDIA(I),EVEC(1,I),1,STPSYM,1)
  200 CONTINUE
      STEPS = DDOT(NCART,STPSYM,1,STPSYM,1)
      STEPS = SQRT(STEPS)
      IF (ISTTYP .EQ. 3) THEN
         STPNRM = STEPS/SQRT(XFAMU)
      ELSE
         STPNRM = STEPS
      END IF
C
C     ***** Print *****
C
      IF (IPRINT .GT. 2) THEN
         IF (IPRINT .GT. 2) THEN
            CALL HEADER('Scaled Cartesian symmetry step vector',-1)
            CALL OUTPUT(STPSYM,1,1,1,NCART,1,NCART,1,LUPRI)
         END IF
         WRITE (LUPRI,'(/A,F15.10/)') ' Norm of scaled step:',STEPS
         IF (IPRINT .GT. 5) THEN
            CALL HEADER('Scale factors',-1)
            CALL OUTPUT(SCAL,1,1,1,NCART,1,NCART,1,LUPRI)
         END IF
      END IF
C
C     ******************************************
C     ***** Unscaled Cartesian step vector *****
C     ******************************************
C
      DO 300 I = 1, NCART
         STPSYM(I) = STPSYM(I)*SCAL(I)
  300 CONTINUE
      STEPU = DDOT(NCART,STPSYM,1,STPSYM,1)
      STEPU = SQRT(STEPU)
C
C     ***** Print *****
C
      IF (IPRINT .GT. 2) THEN
         CALL HEADER('Unscaled Cartesian step vector',-1)
         CALL OUTPUT(STPSYM,1,1,1,NCART,1,NCART,1,LUPRI)
         WRITE (LUPRI,'(/A,F15.10,/)') ' Norm of unscaled step:',STEPU
      END IF
C
C     ************************************
C     ***** Non-symmetry coordinates *****
C     ************************************
C
      CALL WLKCOR(STPSYM,STPCAR,NCART,NCRIND,IPRINT)
C
C     *****************************
C     ***** Cartesian momenta *****
C     *****************************
C
      IF (ISTTYP .EQ. 3) CALL WLKPCR(PMDIA,EVEC,PMCAR,SCAL,NCART,IPRINT)
      RETURN
      END
C  /* Deck wlkneg */
      SUBROUTINE WLKNEG(EVAL,GRDDIA,STPDIA,PMDIA,NINTER,IMODE,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
      LOGICAL FIRST
      DIMENSION EVAL(*), GRDDIA(*), STPDIA(*), PMDIA(*)
      SAVE EMINUS, GMINUS, FIRST, IMINUS
      DATA FIRST /.TRUE./
C
      IF (FIRST) THEN
         EMINUS = - EVAL  (IMODE)
         GMINUS = - GRDDIA(IMODE)
         EVAL  (IMODE) = EMINUS
         GRDDIA(IMODE) = DUMMY
         CALL ORDER(GRDDIA,EVAL,NINTER,1)
         DO 100 I = 1, NINTER
            IF ((EVAL(I).EQ.EMINUS).AND.(GRDDIA(I).EQ.DUMMY)) IMINUS = I
  100    CONTINUE
         GRDDIA(IMINUS) = GMINUS
         FIRST = .FALSE.
      ELSE
         EVAL  (IMINUS) = - EMINUS
         GRDDIA(IMINUS) = - GMINUS
         CALL DCOPY(NINTER,EVAL,1,PMDIA,1)
         CALL ORDER(GRDDIA,EVAL, NINTER,1)
         CALL ORDER(STPDIA,PMDIA,NINTER,1)
         FIRST = .TRUE.
      END IF
      RETURN
      END
C  /* Deck wlkcor */
      SUBROUTINE WLKCOR(STPSYM,STPCAR,NCART,NCRIND,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
C
#include <nuclei.h>
#include <symmet.h>
C
      DIMENSION STPSYM(NCART), STPCAR(NCRIND)
C
#include <ibtfun.h>
C
C     Transforms coordinates to non-symmetry basis
C
      CALL DZERO(STPCAR,NCRIND)
      DO 100 IATOM = 1, NUCIND
         DO 200 ICOOR = 1, 3
            ICCOOR = 3*(IATOM - 1) + ICOOR
            ISCOOR = IPTCNT(ICCOOR,0,1)
            IF (ISCOOR .GT. 0) THEN
              STPCAR(ICCOOR) = STPSYM(ISCOOR)/SQRT(FMULT(ISTBNU(IATOM)))
            END IF
  200    CONTINUE
  100 CONTINUE
      IF (IPRINT .GT. 2) THEN
         CALL HEADER('Cartesian step vector in non-symmetry basis',-1)
         CALL PRIGEO(STPCAR)
      END IF
      RETURN
      END
C  /* Deck wlkpcr */
      SUBROUTINE WLKPCR(PMDIA,EVEC,PMCAR,SCAL,NCART,IPRINT)
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION PMDIA(NCART),PMCAR(NCART),SCAL(NCART),EVEC(NCART,NCART)
C
      CALL DZERO(PMCAR,NCART)
      DO 100 I = 1, NCART
         CALL DAXPY(NCART,PMDIA(I),EVEC(1,I),1,PMCAR,1)
  100 CONTINUE
      IF (IPRINT .GT. 2) THEN
         CALL HEADER('Scaled Cartesian momenta in WLKPCR',-1)
         CALL OUTPUT(PMCAR,1,1,1,NCART,1,NCART,1,LUPRI)
      END IF
      DO 200 I = 1, NCART
         PMCAR(I) = PMCAR(I)/SCAL(I)
  200 CONTINUE
      IF (IPRINT .GT. 0) THEN
         CALL HEADER('Cartesian momenta in WLKPCR',-1)
         CALL OUTPUT(PMCAR,1,1,1,NCART,1,NCART,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck winkel */
      FUNCTION WINKEL(VEC1,VEC2,NDIM)
#include <implicit.h>
#include <facang.h>
      DIMENSION VEC1(NDIM), VEC2(NDIM)
      VEC12  = DDOT(NDIM,VEC1,1,VEC2,1)
      VEC1N  = DNRM2(NDIM,VEC1,1)
      VEC2N  = DNRM2(NDIM,VEC2,1)
      WINKEL = FACANG*ACOS(VEC12/(VEC1N*VEC2N))
      RETURN
      END
C  /* Deck wlkadp */
      SUBROUTINE WLKADP(REDFAC,XMXNUC,STPCAR,NCRIND,ADPOK,IPRINT)
C
C     Largest atomic displacement
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0 D0, D1 = 1.0 D0, DP99 = 0.99 D0)
      LOGICAL ADPOK
      DIMENSION STPCAR(NCRIND)
C
      ADPMAX = D0
      DO 100 IATOM = 1, NCRIND, 3
         XN = DNRM2(3,STPCAR(IATOM),1)
         ADPMAX = MAX(ADPMAX,XN)
  100 CONTINUE
      IF (ADPMAX .LT. XMXNUC) THEN
         ADPOK = .TRUE.
         REDFAC = D1
      ELSE
         ADPOK = .FALSE.
         REDFAC = DP99*XMXNUC/ADPMAX
      END IF
      IF ((IPRINT.GT.2) .OR. (.NOT.ADPOK)) THEN
         WRITE(LUPRI,'(/A,F15.10,/,A,I3,A)')
     *       ' Largest atomic displacement:',ADPMAX,
     *       ' Reduction of trust radius  :',NINT(100*(D1 - REDFAC)),'%'
      END IF
      RETURN
      END
C  /* Deck wlkcan */
      SUBROUTINE WLKCAN(REDANH,STPDIA,TRSDIA,NCART,ANHOK,IPRINT)
C
C     Check if anharmonic contributions dominate
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1 = 1.0 D0, DP95 = 0.95D0)
      LOGICAL ANHOK
      DIMENSION STPDIA(NCART), TRSDIA(NCART)
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('STPDIA in WLKCAN',-1)
         CALL OUTPUT(STPDIA,1,1,1,NCART,1,NCART,1,LUPRI)
         CALL HEADER('TRSDIA in WLKCAN',-1)
         CALL OUTPUT(TRSDIA,1,1,1,NCART,1,NCART,1,LUPRI)
      END IF
C
      ANHOK  = .TRUE.
      REDANH = D1
      DO 100 I = 1, NCART
         IF (ABS(STPDIA(I)) .GT. TRSDIA(I)) THEN
            ANHOK = .FALSE.
            RED = DP95*TRSDIA(I)/ABS(STPDIA(I))
            REDANH = MIN(REDANH,RED)
         END IF
  100 CONTINUE
      IF (.NOT.ANHOK) THEN
         WRITE (LUPRI,'(/A,F12.6)')
     &      ' Reduction factor for trust radius in WLKCAN:' , REDANH
      END IF
      RETURN
      END
C  /* Deck wlkprd */
      SUBROUTINE WLKPRD(ERGSEC,GRDCAR,HESCAR,STPCAR,NCART)
C
C     13-Jun-1985/29-Nov-1985 hjaaj
C     Calculate second order prediction for energy change.
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      DIMENSION GRDCAR(NCART), HESCAR(NCART,NCART), STPCAR(NCART)
      PARAMETER ( DP5 = 0.5D0 )
C
      ERGSCC = DDOT(NCART,GRDCAR,1,STPCAR,1)
      DO 100 J = 1, NCART
         STPCOJ = STPCAR(J)
         DO 200 I = 1, (J-1)
            ERGSCC = ERGSCC + STPCOJ*HESCAR(I,J)*STPCAR(I)
  200    CONTINUE
         ERGSCC = ERGSCC + DP5 * STPCOJ*STPCOJ*HESCAR(J,J)
  100 CONTINUE
      IF (ABS(ERGSCC - ERGSEC) .GT. 1.0D-15) THEN
         WRITE (LUPRI,'(//A,2(/A,F20.15))')
     *      ' ERROR in WLKPRD, program cannot proceed.',
     *      ' Second-order energy in diagonal representation: ',ERGSEC,
     *      ' Second-order energy in Cartesian representation:',ERGSCC
         CALL QUIT('ERROR in WLKPRD.')
      END IF
      RETURN
      END
C  /* Deck wlkfl0 */
      SUBROUTINE WLKFL0(NCORD,NONTRO,EVAL,GRDDIA,STPDIA,
     *                  TRUSTR,RNU,KEEPSY,ZERGRD,IPRINT)
C
C     This subroutine solves the constrained restricted step
C     equations (the level-shifted Newton equations) in the
C     diagonal representation.  We assume that the Newton step
C     is longer than the trust radius.
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER ( D0 = 0.0D0 , DP5 = 0.5D0 )
C
      DIMENSION EVAL(*), GRDDIA(*), STPDIA(*)
      LOGICAL KEEPSY, SPECAS
C
      EXTERNAL WLKSTL
C
      IF (IPRINT .GT. 5) CALL HEADER('OUTPUT FROM WLKFL0',-1)
C
      IF (KEEPSY) THEN
         DO 10 I = 1, NONTRO
            IF (ABS(GRDDIA(I)) .LT. ZERGRD) GRDDIA(I) = D0
   10    CONTINUE
         STPSYM = WLKSTL(GRDDIA(1),EVAL(1),D0,NONTRO,D0)
         WRITE (LUPRI,'(/A,1P,D12.5)')
     *      ' Length of non-symmetry-breaking Newton step: ',STPSYM
         NEGSYM = 0
         DO 50 I = 1, NONTRO
            IF (GRDDIA(I) .NE. D0 .AND. EVAL(I) .LT. D0) NEGSYM=NEGSYM+1
   50    CONTINUE
         WRITE (LUPRI,'(/2A,I5)')
     *      ' Number of negative Hessian eigenvalues corresponding to ',
     *      ' non-symmetry-breaking eigenvectors: ',NEGSYM
         IF (STPSYM .LT. TRUSTR .AND. NEGSYM .EQ. 0) THEN
             WRITE (LUPRI,'(/A)') ' Newton step is taken .'
             DO 60 I = 1, NONTRO
                STPDIA(I) = - GRDDIA(I)/EVAL(I)
   60        CONTINUE
             RETURN
         END IF
         DO 20 I = 1, NONTRO
            IF (GRDDIA(I) .NE. D0) THEN
              GRDMIN = ABS(GRDDIA(I))
              HESMIN = EVAL(I)
              GO TO 30
            END IF
   20    CONTINUE
   30    CONTINUE
         SPECAS = .FALSE.
      ELSE
         HESMIN = EVAL(1)
         GRDMIN = ABS(GRDDIA(1))
C
C        Test whether the lowest Hessian eigenvalue is negative and the
C        corresponding gradient zero. This case is treated separately
C        as described by Fletcher in "Unconstrained Optimization" p.85.
C
         SPECAS = (HESMIN .LT. D0) .AND. (GRDMIN .LT. ZERGRD)
      END IF
C
      GRDNRM = SQRT(DDOT(NONTRO,GRDDIA,1,GRDDIA,1))
C
      IF (IPRINT .GT. 5) THEN
         WRITE (LUPRI,'(A,1P,D12.5)') ' HESMIN: ', HESMIN
         WRITE (LUPRI,'(A,1P,D12.5)') ' GRDMIN: ', GRDMIN
         WRITE (LUPRI,'(A,1P,D12.5)') ' GRDNRM: ', GRDNRM
      END IF
C
C     ************************
C     ***** General case *****
C     ************************
C
      IF (.NOT. SPECAS) THEN
         IF (IPRINT .GT. 3) WRITE (LUPRI,'(/A/)') ' General case.'
C
C        Determine level shift
C
         XMIN = GRDNRM/TRUSTR - MIN(HESMIN,D0)
         XMAX = MAX(D0, - HESMIN + DP5*GRDMIN/TRUSTR)
         IF (IPRINT .GT. 5) THEN
            WRITE (LUPRI,'(/A,2(1P,D12.5,2X))')
     *         ' XMIN and XMAX before WLKBIS (WLKSTL) minimum walk: ',
     *          XMIN,XMAX
         END IF
         CALL WLKBIS(XMAX,XMIN,RNU,GRDDIA,EVAL,TRUSTR,NONTRO,
     *               WLKSTL,IFAIL)
         IF (IPRINT .GT. 5) THEN
            WRITE (LUPRI,'(/A,2(1P,D12.5,2X))')
     *         ' XMIN and XMAX after WLKBIS (WLKSTL) minimum walk: ',
     *          XMIN,XMAX
         END IF
         IF (IPRINT .GT. 3) WRITE (LUPRI,'(/A,1P,D12.5)')
     *      ' Level shift parameter: ', RNU
         IF (IFAIL.EQ.0) THEN
            WRITE (LUPRI,5250) XMAX, XMIN
            CALL QUIT(' *** ERROR, Wrong interval in WLKBIS (WLKSTL)')
         ELSE IF (IFAIL.EQ.1) THEN
            WRITE (LUPRI,5350)
         END IF
C
C        Determine step vector
C
         DO 100 I = 1, NONTRO
            STPDIA(I) = - GRDDIA(I)/(EVAL(I) + RNU)
 100     CONTINUE
C
C     *************************************************
C     ***** Special case: HESMIN < 0 & GRDMIN = 0 *****
C     *************************************************
C
      ELSE
         IF (IPRINT .GT. 3) WRITE (LUPRI,'(/A/)') ' Special case.'
C
C        Set RNU = - HESMIN and determine step length
C
         STPNRM = WLKSTL(GRDDIA(2),EVAL(2),-HESMIN,NONTRO-1,D0)
         IF (IPRINT .GT. 3) THEN
            WRITE (LUPRI,'(/A,F12.6)')
     *      ' Step length with level shift equal to lowest eigenvalue:',
     *      STPNRM
         END IF
         IF (STPNRM .GT. TRUSTR) THEN
            IF (IPRINT .GT. 3) WRITE (LUPRI,'(/A)')
     *         ' Component along lowest eigenvector ignored.'
C
C           Determine step vector in the usual way ignoring the
C           component along the lowest eigenvector. We now know that
C           the level shift must be greater than - HESMIN.
C
            XMIN = GRDNRM/TRUSTR - MIN(HESMIN,D0)
            XMAX = - HESMIN
            IF (IPRINT .GT. 5) THEN
               WRITE (LUPRI,'(/A,2(1P,D12.5,2X))')
     *           ' XMIN and XMAX before WLKBIS (WLKSTL) minimum walk: ',
     *           XMIN,XMAX
            END IF
            CALL WLKBIS(XMAX,XMIN,RNU,GRDDIA(2),EVAL(2),TRUSTR,NONTRO-1,
     *                  WLKSTL,IFAIL)
            IF (IPRINT .GT. 5) THEN
               WRITE (LUPRI,'(/A,2(1P,D12.5,2X))')
     *         ' XMIN and XMAX after WLKBIS (WLKSTL) minimum walk: ',
     *           XMIN,XMAX
            END IF
            IF (IPRINT .GT. 3) WRITE (LUPRI,'(/A,1P,D12.5)')
     *         ' Level shift parameter: ', RNU
            IF (IFAIL.EQ.0) THEN
               WRITE (LUPRI,5250) XMAX, XMIN
               CALL QUIT
     *            (' *** ERROR, Wrong interval in WLKBIS (WLKSTL)')
            ELSE IF (IFAIL.EQ.1) THEN
               WRITE (LUPRI,5350)
            END IF
C
C           Determine step vector
C
            STPDIA(1) = D0
            DO 200 I = 2, NONTRO
               STPDIA(I) = - GRDDIA(I)/(EVAL(I) + RNU)
 200        CONTINUE
         ELSE
C
C           Determine step vector with level shift - HESMIN and add
C           component along the lowest eigenvector to insure that total
C           step length is equal to the trust radius.
C
            DO 300 I = 2, NONTRO
               STPDIA(I) = - GRDDIA(I)/(EVAL(I) - HESMIN)
 300        CONTINUE
            STP2      = DDOT(NONTRO-1,STPDIA(2),1,STPDIA(2),1)
            IF (IPRINT .GT. 3) WRITE (LUPRI,'(/A,F12.6)')
     *         ' Norm of step orthogonal to lowest eigenvector:',
     *         SQRT(STP2)
            STPDIA(1) = SQRT(TRUSTR*TRUSTR - STP2)
         END IF
      END IF
      IF (IPRINT .GT. 2) THEN
         CALL HEADER('STPDIA after min search',-1)
         WRITE (LUPRI,'(5X,3F15.8)') (STPDIA(I),I=1,NCORD)
      END IF
      RETURN
C
C     FORMATS
C
 5250 FORMAT(/' *** Wrong interval specified in WLKBIS (WLKSTL) ***',
     *       /' XMAX= ',F10.6,'   XMIN= ',F10.6)
 5350 FORMAT(/' *** WARNING WLKBIS (WLKSTL) ***',
     *       /' Desired accuracy not obtained in the specified maximum',
     *       /' number of iterations.')
C
C     End of WLKFL0
C
      END
C  /* Deck wlkfl1 */
      SUBROUTINE WLKFL1(NNEG,NMODE,NCORD,NONTRO,EVAL,GRDDIA,STPDIA,
     *                  TRUSTR,RNU,IPRINT)
C
C     Calculate step vector using the level-shift (Fletcher) algorithm
C     for walks towards saddle point (following mode NMODE).
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
      DIMENSION EVAL(*), GRDDIA(*), STPDIA(*)
      PARAMETER ( D0 = 0.0D0 , D1 = 1.0D0 , D2 = 2.0D0 , D4 = 4.0D0)
      PARAMETER ( DP5 = 0.5D0 , DP25 = 0.25D0 )
      PARAMETER ( DADD = 1.0D-6 , RTOL = 1.2D0 )
C
      EXTERNAL WLKDER, WLKSTL
C
C     Determine scale factor.
C
      SCLFAC = D1
      IF (NNEG.EQ.0) THEN
         IF (NMODE.EQ.1) THEN
            IF (EVAL(1).GT.DP5*EVAL(2)) THEN
               SCLFAC = DP25 * (EVAL(2)/EVAL(1))
            END IF
         ELSE
            SCLFAC = DP25 * (EVAL(1)/EVAL(NMODE))
         END IF
      ELSE IF (NMODE.EQ.1) THEN
         IF (EVAL(1).GT.D2*EVAL(2)) THEN
            SCLFAC = D4 * (EVAL(2)/EVAL(1))
         END IF
      ELSE IF (NMODE .LE. NNEG) THEN
         SCLFAC = D4 * (EVAL(1)/EVAL(NMODE))
      ELSE
         WRITE (LUPRI,'(/A//A,I5,A,/A,I5/)')
     *      ' Incorrect Hessian structure for level shift.',
     *      ' There are ',NNEG,' negative eigenvalues.',
     *      ' The mode determined by overlap is:',NMODE
         CALL QUIT('(WLKCTL) Bad Hessian structure.')
      END IF
C
C     Transform to scaled coordinate system
C
      IF (NMODE .NE. 1 .OR. SCLFAC .NE. D1) THEN
         ENMODE = EVAL(NMODE)
         GNMODE = GRDDIA(NMODE)
         DO 270 I = NMODE,2,-1
            EVAL(I)   = EVAL(I-1)
            GRDDIA(I) = GRDDIA(I-1)
  270    CONTINUE
         EVAL(1)   = SCLFAC * ENMODE
         GRDDIA(1) = SCLFAC * GNMODE
C
         IF (IPRINT .GE. 1) WRITE (LUPRI,'(/A,I3,/A,F10.5)')
     *      ' Scale factor to bring streambed mode',NMODE,
     *      ' down as the lowest mode:',SCLFAC
         IF (IPRINT .GT. 5) THEN
            WRITE (LUPRI,'(/A/)')
     *       ' Eigenvalues of projected Hessian after scaling'
            WRITE (LUPRI,'(5X,3F15.6)') (EVAL(I),I=1,NCORD)
            WRITE (LUPRI,'(/A/)')
     *       ' Gradient in diagonal representation after scaling'
            WRITE (LUPRI,'(5X,3F15.6)') (GRDDIA(I),I=1,NCORD)
         END IF
      END IF
C
C     Determine minimum of step length function for WALKS
C     to transition states
C
      XMAX = - EVAL(1) - DADD
      XMIN = - EVAL(2) + DADD
      CALL WLKBIS(XMAX,XMIN,XDET,GRDDIA,EVAL,DUMMY,NONTRO,
     *            WLKDER,IFAIL)
      IF (IPRINT .GT. 3) THEN
         WRITE (LUPRI,'(A,F15.6,A,F15.6)')
     *      ' Level shift at zero derivative is between',
     *      XMIN,' and',XMAX
         WRITE (LUPRI,'(/A,F15.6)')
     *      ' Level shift used:',XDET
      END IF
      IF (IFAIL.EQ.0) THEN
         WRITE (LUPRI,5200)XMAX,XMIN
         CALL QUIT(' *** ERROR, Wrong interval in WLKBIS (WLKDER)')
      ELSE IF (IFAIL.EQ.1) THEN
         WRITE (LUPRI,5300)
      ENDIF
      X = WLKSTL(GRDDIA,EVAL,XDET,NONTRO,D0)
      WRITE (LUPRI,'(/A,F8.4)')
     *   ' Step length at zero derivative',X
C
 5200 FORMAT(/' *** Wrong interval specified in WLKBIS ***',
     *      /' Location of zero value of derivative step failed',
     *      /' XMAX = ',F10.6,'   XMIN = ',F10.6)
 5300 FORMAT(/' *** WARNING WLKBIS (WLKDER) ***',
     *      /' Desired accuracy not obtained in the specified',
     *      /' maximum number of iterations.')
C
C     Compare step at derivative zero with TRUSTR.
C     If no level shift exists such that step length function
C     has a value of TRUSTR, redefine level shift.
C
      TRUSTN = TRUSTR
      IF (TRUSTR .LT. X*RTOL) TRUSTR = RTOL*X
C
C     Determine interval for level shift parameter
C
      XMIN = XDET
      IF ( EVAL(1)*EVAL(2) .LT. D0 ) THEN
C        Newton step fulfills first and second order
C        requirements, select interval including
C        Newton step (i.e. level shift 0).
         IF (XDET .GT. D0) THEN
            XMAX = - EVAL(2) + DADD
         ELSE
            XMAX = - EVAL(1) - DADD
         END IF
      ELSE IF (ABS(EVAL(2)) .LT. ABS(EVAL(1))) THEN
         XMAX = - EVAL(2) + DADD
      ELSE
         XMAX = - EVAL(1) - DADD
      END IF
      IF (IPRINT .GT. 5) THEN
         WRITE (LUPRI,'(/A/5X,2F15.10)')
     *      ' XMIN,XMAX before WLKBIS (WLKSTL)',XMIN,XMAX
         WRITE (LUPRI,'(A,5X,F15.10,I10)')
     *      ' temporary TRUSTR, NONTRO :',TRUSTR,NONTRO
      END IF
      CALL WLKBIS(XMAX,XMIN,RNU,GRDDIA,EVAL,TRUSTR,NONTRO,
     *              WLKSTL,IFAIL)
      IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A/5X,3F15.10)')
     *      ' XMIN,XMAX,RNU after WLKBIS (WLKSTL)',XMIN,XMAX,RNU
      IF (IFAIL.EQ.0) THEN
          WRITE (LUPRI,5250)XMAX,XMIN
          CALL QUIT('*** ERROR, Wrong interval in WLKBIS (WLKSTL)')
      ELSE IF (IFAIL.EQ.1) THEN
         WRITE (LUPRI,5350)
      END IF
      IF (IPRINT .GE. 3) WRITE (LUPRI,'(/A,F15.10)')
     *      ' Restricted-step level shift =',RNU
      DO 830 I=1,NONTRO
         STPDIA(I) = - GRDDIA(I)/(EVAL(I) + RNU)
  830 CONTINUE
      IF (IPRINT .GE. 4) THEN
         WRITE (LUPRI,'(/A)')
     *      ' Step in diagonal Hessian representation :'
         WRITE (LUPRI,'(5X,3F15.8)') (STPDIA(I),I=1,NCORD)
      END IF
      TRUSTR = TRUSTN
      XSTEP  = DDOT(NONTRO,STPDIA,1,STPDIA,1)
      XSTEP  = TRUSTR/SQRT(XSTEP)
      CALL DSCAL(NONTRO,XSTEP,STPDIA,1)
      IF (IPRINT .GT. 2) THEN
         WRITE (LUPRI,'(/A,F8.4)')
     *   ' Step scaled to trust radius =',TRUSTR
         WRITE (LUPRI,'(/A)') ' Step in diagonal Hessian'//
     *      ' representation after scaling :'
         WRITE (LUPRI,'(5X,3F15.8)') (STPDIA(I),I=1,NCORD)
      END IF
      IF (NMODE .NE. 1 .OR. SCLFAC .NE. D1) THEN
         ENMODE = EVAL(1)   / SCLFAC
         GNMODE = GRDDIA(1) / SCLFAC
         SNMODE = STPDIA(1)
         DO 810 I = 2, NMODE
            EVAL(I-1)   = EVAL(I)
            GRDDIA(I-1) = GRDDIA(I)
            STPDIA(I-1) = STPDIA(I)
  810    CONTINUE
         EVAL(NMODE)   = ENMODE
         GRDDIA(NMODE) = GNMODE
         STPDIA(NMODE) = SNMODE
         IF (IPRINT .GT. 5) THEN
            WRITE (LUPRI,'(/A)')
     *           ' Step in original eigenvector basis'
            WRITE (LUPRI,'(5X,3F15.8)') (STPDIA(I),I=1,NCORD)
         END IF
      END IF
      RETURN
C
C     FORMATS
C
 5250 FORMAT(/' *** Wrong interval specified in WLKBIS (WLKSTL) ***',
     *       /' XMAX= ',F10.6,'   XMIN= ',F10.6)
 5350 FORMAT(/' *** WARNING WLKBIS (WLKSTL) ***',
     *       /' Desired accuracy not obtained in the specified maximum',
     *       /' number of iterations.')
C
C     END OF WLKFL1
C
      END
C  /* Deck wlkext */
      SUBROUTINE WLKEXT(NMODE,NCORD,NONTRO,EVAL,GRDDIA,STPDIA,
     *                  TRUSTR,KEEPSY,ZERGRD,ANHRAD,IWKIND,IPRINT)
C
C     Use gradient-extremal walking algorithm
C
C     IWKIND = 0 : Minimization
C            = 1 : Walk to transition state
C            > 2 : Walk to stationary point with unknown index
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
      LOGICAL KEEPSY
      DIMENSION EVAL(*), GRDDIA(*), STPDIA(*)
      PARAMETER ( D0 = 0.0D0, D2 = 2.0D0, D3 = 3.0D0 )
      PARAMETER ( D6 = 6.0D0, DP5 = 0.5D0 )
      PARAMETER ( ZRHESS = 1.0D-2 )
C
C     Check if bifurcation occurs in connection with symmetry breaking
C
      IF (.NOT.KEEPSY) THEN
         DO 50 I = 1, NONTRO
            IF((EVAL(I).LT.ZRHESS).AND.(ABS(GRDDIA(I)).LT.ZERGRD)) THEN
               XSTP = D0
               DO 45 J = 1, NONTRO
                  IF ((J.NE.I) .AND. (J.NE.NMODE)) THEN
                     STPDIA(J) = - GRDDIA(J)/EVAL(J)
                     XSTP      = XSTP + STPDIA(J)*STPDIA(J)
                  END IF
 45            CONTINUE
               XSTP = SQRT(XSTP)
C
C              Be sure you are on the gradient extremal
C
               IF (XSTP*D3 .GT. TRUSTR) THEN
                  STPDIA(I)     = D0
                  STPDIA(NMODE) = D0
                  IF (XSTP .GT. TRUSTR) THEN
                     CALL DSCAL(NONTRO,(TRUSTR/XSTP),STPDIA,1)
                  END IF
                  WRITE(LUPRI,'(/2A/2A,F15.6)')
     *               ' Walk adjusted to be on the gradient ',
     *               ' extremal in a region close to a bifurcation',
     *               ' step in the extremal directions ',
     *               ' are set to zero. Step length:',XSTP
                  GO TO 9999
               END IF
C
C              Predict step length in new gradient-extremal direction by
C              requiring that this step to first order gives a gradient
C              of the size 1/6 times gradient in the old gradient-
C              extremal direction.
C
               STPNMO    = ABS((GRDDIA(NMODE)/D6)/EVAL(I))
               STPNMO    = MIN (TRUSTR / D2,STPNMO)
               STPDIA(I) = STPNMO
               TRUSTR    = STPNMO
               STPDIA(NMODE) = D0
               NMODE = I
               WRITE(LUPRI,'(/A,I5,A//A,F10.6,A,F10.6)')
     *            ' Bifurcation encountered. Mode No.',I,
     *            ' which breaks symmetry is followed.',
     *            ' Hessian eigenvalue:',EVAL(I),
     *            ', Gradient:',GRDDIA(I)
               GO TO 9999
            END IF
 50      CONTINUE
      END IF
C
C     Determine Newton step in all directions orthogonal to the
C     gradient extremal. if walk to transition state (IWKIND.EQ.1) then
C     ensure downhill movement in all directions orthogonal to the
C     gradient extremal.
C
      XSTP = D0
      DO 100 I = 1, NONTRO
         STPDIA(I) = - GRDDIA(I)/EVAL(I)
         IF ((IWKIND.EQ.1) .AND. (I.NE.NMODE) .AND. (EVAL(I).LT.D0))THEN
            STPDIA(I) = -STPDIA(I)
            WRITE(LUPRI,'(/2A,I5)')
     *         'Sign has been changed on Newton step along',
     *         'Hessian eigendirection No.:',I
         END IF
         IF ((EVAL(I).LT.D0).AND.(I.NE.NMODE)) THEN
            WRITE(LUPRI,'(/A,I3,A)')
     *       ' *** WARNING - Hessian eigenvalue No.',I,' is negative'
         END IF
         IF (I.NE.NMODE) XSTP = XSTP + STPDIA(I)*STPDIA(I)
 100  CONTINUE
      XSTP = SQRT(XSTP)
      IF (IPRINT.GT.2) WRITE(LUPRI,'(/A,F10.6)')
     *   ' Step orthogonal to gradient extremal: ',XSTP
      IF (XSTP.GT.TRUSTR) THEN
         STPNMO = D0
         SCALFC = TRUSTR/XSTP
         DO 300 I = 1, NONTRO
            IF (I.NE.NMODE) STPDIA(I) = STPDIA(I)*SCALFC
 300     CONTINUE
         WRITE(LUPRI,'(/A/A,F15.6,2A)')
     *      ' *** WARNING *** ',
     *      ' Step of length',TRUSTR,' is taken exclusively in ',
     *      ' directions orthogonal to the gradient extremal.'
      ELSE
         STPNMO = SQRT(TRUSTR*TRUSTR - XSTP*XSTP)
      END IF
      STPNMO = MIN(STPNMO,DP5*ANHRAD)
      STPSGN = STPDIA(NMODE)/ABS(STPDIA(NMODE))
      IF (IWKIND.GT.0) THEN
         IF (EVAL(NMODE).LT.D0) THEN
            STPDIA(NMODE) = STPNMO*STPSGN
         ELSE
            STPDIA(NMODE) = -STPNMO*STPSGN
         END IF
      ELSE
         IF (EVAL(NMODE).GT.D0) THEN
            STPDIA(NMODE) = STPNMO*STPSGN
         ELSE
            STPDIA(NMODE) = -STPNMO*STPSGN
         END IF
      END IF
 9999 CONTINUE
      RETURN
C
C     END WLKEXT
C
      END
C  /* Deck wlkdyn */
      SUBROUTINE WLKDYN(TIMSTP,Q,P,GRDDIA,EVAL,NCART,RADIUS,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <codata.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D10 = 10.0D0)
      PARAMETER (EFFZER = 1.0D-10)
      PARAMETER (MAXIT = 1000, D100 = 1.0D02)
      DIMENSION GRDDIA(NCART), EVAL(NCART), Q(NCART), P(NCART)
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Gradient in WLKDYN',-1)
         CALL OUTPUT(GRDDIA,1,1,1,NCART,1,NCART,1,LUPRI)
         CALL HEADER('Eigenvalues in WLKDYN',-1)
         CALL OUTPUT(EVAL,1,1,1,NCART,1,NCART,1,LUPRI)
         WRITE (LUPRI,'(1X,A,F12.6)') ' RADIUS ', RADIUS
      END IF
C
      TRADAU = RADIUS*SQRT(XFAMU)
C
C     Find highest eigenvalue
C
      HESMAX = EVAL(1)
      DO 100 I = 2, NCART
         HESMAX = MAX(HESMAX,EVAL(I))
  100 CONTINUE
C
C     Find lowest eigenvalue
C
      HESMIN = EVAL(NCART)
      DO 110 I = 1, NCART - 1
         HESMIN = MIN(HESMIN,EVAL(I))
  110 CONTINUE
      IF (IPRINT .GT. 5) THEN
         WRITE (LUPRI,'(1X,A,F12.6)') ' HESMAX ', HESMAX
         WRITE (LUPRI,'(1X,A,F12.6)') ' HESMIN ', HESMIN
      END IF
C
C     If all non-zero eigenvalues are positive,
C     set time limit for second-order walk
C     to about one vibrational period of lowest mode
C
      IF (HESMIN .LT. EFFZER) THEN
         TLIM = - D1
      ELSE
         HESLOW = EVAL(NCART)
         DO 150 I = 1, NCART - 1
            IF (EVAL(I) .GT. EFFZER) THEN
               COSFAC = ABS(GRDDIA(I)/EVAL(I) + Q(I))
               SINFAC = ABS(P(I)/SQRT(EVAL(I)))
               IF ((COSFAC .GT. EFFZER) .OR. (SINFAC .GT. EFFZER)) THEN
                  HESLOW = MIN(HESLOW,EVAL(I))
               END IF
            END IF
  150    CONTINUE
         IF (HESLOW .GT. EFFZER) THEN
            TLIM = D10/SQRT(HESLOW)
            IF (IPRINT .GT. 5) THEN
               WRITE (LUPRI,'(A,F24.12)')
     *         ' Time limit ', TLIM
            END IF
         END IF
      END IF
C
C     Determine initial time step
C
      DELTAT = D1/(D10*SQRT(HESMAX))
      RAD  = D0
      STEP = WLKSTD(DELTAT,Q,P,GRDDIA,EVAL,NCART,RAD)
      DELTAT = (DELTAT*TRADAU)/(D10*STEP)
      IF (IPRINT .GT. 5) THEN
         WRITE (LUPRI,'(A,F24.12)')
     *   ' Trust radius in mass-weighted atomic units ', TRADAU
         WRITE (LUPRI,'(A,F24.12)')
     *   ' Lowest vibrating eigenvalue in WLKDYN ', HESLOW
         WRITE (LUPRI,'(A,F24.12)')
     *   ' Largest eigenvalue in WLKDYN ', HESMAX
         WRITE (LUPRI,'(A,F24.12)')
     *   ' Initial time step ', DELTAT
      END IF
C
C     Determine TMIN and TMAX bracketing TIMSTP for bisection
C
      TMAX = D0
      STPMAX = WLKSTD(TMAX,Q,P,GRDDIA,EVAL,NCART,TRADAU)
      ITER = 0
  200 CONTINUE
         ITER   = ITER + 1
         TMIN   = TMAX
         TMAX   = TMIN + DELTAT
         IF ((TLIM. GT. D0) .AND. (TMAX .GT. TLIM)) THEN
            WRITE (LUPRI,'(//A,I5,A)')
     *                    ' Time limit exceeded in WLKDYN.'
            WRITE (LUPRI,'(A,1P,D12.4)') ' Time limit ', TLIM
            WRITE (LUPRI,'(A,1P,D12.4)') ' Tmax ',       TMAX
            CALL QUIT('Time limit exceeded in WLKDYN.')
         END IF
         STPMIN = STPMAX
         STPMAX = WLKSTD(TMAX,Q,P,GRDDIA,EVAL,NCART,TRADAU)
         IF (IPRINT .GT. 10) THEN
            WRITE (LUPRI, '(/A,I5)' ) ' Iteration in bisection: ',ITER
            WRITE (LUPRI, '(A,1P,D12.4)') ' TMAX   ', TMAX
            WRITE (LUPRI, '(A,1P,D12.4)') ' STPMAX ', STPMAX
         END IF
         IF (STPMAX .GT. D0) GO TO 300
      IF (ITER .LT. MAXIT) GO TO 200
C
      WRITE (LUPRI,'(//A,I5,A)') ' Bracketing failed after ',
     *                             ITER, ' iterations in WLKDYN.'
      CALL QUIT('Bracketing failed in WLKDYN.')
  300 CONTINUE
      IF (IPRINT .GT. 5) THEN
         WRITE (LUPRI,'(A,1P,2D12.4)')
     *      ' TMIN and TMAX in WLKDYN ', TMIN, TMAX
         WRITE (LUPRI,'(A,1P,2D12.4)')
     *      ' STPMIN and STPMAX in WLKDYN ', STPMIN, STPMAX
         WRITE (LUPRI,'(A,I5)')
     *      ' Number of iterations to find bisection limits',ITER
      END IF
C
C     Bisection
C
      ITER = 0
 400  CONTINUE
         ITER= ITER+1
         TIMSTP= (TMAX+TMIN)/D2
         STP = WLKSTD(TIMSTP,Q,P,GRDDIA,EVAL,NCART,TRADAU)
         IF (STP.GT.D0) TMAX=TIMSTP
         IF (STP.LT.D0) TMIN=TIMSTP
      IF (ABS(STP) .LT. TRADAU/D100) GO TO 401
      IF (ITER .LT. MAXIT) GO TO 400
C
         WRITE (LUPRI,'(//A,I5,A)') ' Bisection failed at after ',ITER,
     &                           ' iterations in WLKDYN.'
         CALL QUIT('Bisection failed in WLKDYN.')
C
 401  CONTINUE
      IF (IPRINT .GT. 2) THEN
         WRITE (LUPRI,'(//A,I5)')  ' Number of bisections   ', ITER
         WRITE (LUPRI,'(A,1P,D12.4)') ' Time to reach boundary ', TIMSTP
         CALL HEADER('Final coordinates and momenta',-1)
         WRITE (LUPRI,'(2A)')
     *    '   #    gradient    Hessian     old coor.  new coor.   ',
     *    '   old mom.    new mom. '
      END IF
      DO 500 K=1,NCART
         Q0   = Q(K)
         P0   = P(K)
         Q(K) = WLKQT(TIMSTP,Q0,P0,GRDDIA(K),EVAL(K))
         P(K) = WLKPT(TIMSTP,Q0,P0,GRDDIA(K),EVAL(K))
         IF (IPRINT .GT. 2) THEN
            WRITE (LUPRI,'(1X,I3,2X,6D12.5)')
     *         K, GRDDIA(K), EVAL(K), Q0, Q(K), P0, P(K)
         END IF
 500  CONTINUE
      STP = WLKSTD(TIMSTP,Q,P,GRDDIA,EVAL,NCART,D0)
      IF (IPRINT .GT. 2) THEN
         WRITE (LUPRI,'(A,1P,D12.4)') ' Norm of step           ', STP
      END IF
      RETURN
      END
C  /* Deck wlkstd */
      FUNCTION WLKSTD(TIMSTP,Q0,P0,GRDDIA,EVAL,NCART,TRUSTR)
C
C     Purpose:
C
C        Calculate step length at time step TIMSTP and
C        subtract TRUSTR
C
C        WLKSTL = //STEP// - TRUSTR
C
#include <implicit.h>
      PARAMETER( D0=0.0D0 )
      DIMENSION Q0(NCART), P0(NCART), GRDDIA(NCART), EVAL(NCART)
C
      STEP = D0
      DO 100 K=1,NCART
         QK = WLKQT(TIMSTP,Q0(K),P0(K),GRDDIA(K),EVAL(K))
         STEP = STEP + (QK - Q0(K))**2
 100  CONTINUE
      WLKSTD = SQRT(STEP) - TRUSTR
      RETURN
      END
C  /* Deck wlktrd */
      FUNCTION WLKTRD(GRAD,HESS,ANHARM,FACTOR)
C
C     Estimate anharmonic threshold
C
#include <implicit.h>
      PARAMETER (D3 = 3.D0, D6 = 6.D0)
      STEP1 = SQRT(D6*ABS(GRAD/(FACTOR*ANHARM)))
      STEP2 =      D3*ABS(HESS/(FACTOR*ANHARM))
      WLKTRD = MAX(STEP1,STEP2)
      RETURN
      END
C  /* Deck wlkqt */
      FUNCTION WLKQT(T,Q0,P0,GDDIA,HESDIA)
C
C     This function gives the coordinate at time T when initial
C     coordinate and momentum are Q0 and P0.
C
C     tuh 890428
C
#include <implicit.h>
      PARAMETER (DP5 = 0.5D0, THRESH = 1.0D-10)
C
      IF (HESDIA .GT. THRESH) THEN
         SQRTHD = SQRT(HESDIA)
         WLKQT = ((GDDIA/HESDIA) + Q0)*COS(SQRTHD*T)
     *         + (P0/SQRTHD)*SIN(SQRTHD*T) - (GDDIA/HESDIA)
      ELSE IF (HESDIA .LT. -THRESH) THEN
         SQRTHD = SQRT(-HESDIA)
         WLKQT = ((GDDIA/HESDIA) + Q0)*COSH(SQRTHD*T)
     *         + (P0/SQRTHD)*SINH(SQRTHD*T) - (GDDIA/HESDIA)
      ELSE
         WLKQT = - DP5*GDDIA*T*T + P0*T + Q0
      END IF
      RETURN
      END
C  /* Deck wlkpt */
      FUNCTION WLKPT(T,Q0,P0,GDDIA,HESDIA)
C
C     This function gives the momentum at time T when initial
C     coordinate and momentum are Q0 and P0.
C
C     tuh 890428
C
#include <implicit.h>
      PARAMETER (THRESH = 1.0D-10)
C
      IF (HESDIA .GT. THRESH) THEN
         SQRTHD = SQRT(HESDIA)
         WLKPT = P0*COS(SQRTHD*T)
     *         - ((GDDIA/SQRTHD) + Q0*SQRTHD)*SIN(SQRTHD*T)
      ELSE IF (HESDIA .LT. -THRESH) THEN
         SQRTHD = SQRT(-HESDIA)
         WLKPT = P0*COSH(SQRTHD*T)
     *         - ((GDDIA/SQRTHD) - Q0*SQRTHD)*SINH(SQRTHD*T)
      ELSE
         WLKPT = - GDDIA*T + P0
      END IF
      RETURN
      END
C  /* Deck wlkirc */
      SUBROUTINE WLKIRC(Q,GRDDIA,EVAL,NCART,RADIUS,DONWTN,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <codata.h>
      PARAMETER (TMAX = 10000.0D0)
      LOGICAL DONWTN
      DIMENSION GRDDIA(NCART), EVAL(NCART), Q(NCART)
#include <abainf.h>
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Gradient in WLKIRC',-1)
         CALL OUTPUT(GRDDIA,1,1,1,NCART,1,NCART,1,LUPRI)
         CALL HEADER('Eigenvalues in WLKIRC',-1)
         CALL OUTPUT(EVAL,1,1,1,NCART,1,NCART,1,LUPRI)
         WRITE (LUPRI,'(1X,A,F12.6)') ' RADIUS ', RADIUS
      END IF
C
C     Take Newton step if:  1) it is smaller than RADIUS
C                           2) all eigenvalues are positive
C
      CALL WLKIRM(GRDDIA,EVAL,Q,NCART,RADIUS,DONWTN,IPRINT)
C
C     Newton step unsatisfactory
C
      IF (.NOT.DONWTN) THEN
C
C        Determine step size using Newton's method such that the
C        integration is to the boundary RADIUS
C
         TIMSTP = WLKIRN(GRDDIA,EVAL,NCART,RADIUS,IPRINT)
         IF (TIMSTP .GT. TMAX) THEN
            WRITE (LUPRI,'(//A,2(/A,1P,E12.6),/A)')
     &         ' >>>>> WARNING: Time step exceeded in WLKIRN!',
     &         ' Time step determined in WLKIRN:', TIMSTP,
     &         ' Maximum time step allowed:     ', TMAX,
     &         ' Calculation proceeds with TIMSTP = TMAX.'
            TIMSTP = TMAX
            NWNABA = NWNABA + 1
         END IF
         IF (IPRINT .GT. 10) THEN
            WRITE (LUPRI,'(/A,F12.6)')
     &         ' Time step after Newton iterations: ', TIMSTP
         END IF
C
C        Get step in diagonal representation
C
         CALL WLKIRQ(TIMSTP,GRDDIA,EVAL,Q,NCART)
      END IF
C
      IF (IPRINT .GT. 10) THEN
         WRITE (LUPRI,'(A,L5)') ' DONWTN ', DONWTN
         IF (.NOT.DONWTN) WRITE (LUPRI,'(A,F12.6)')' Time step ',TIMSTP
         CALL HEADER('Step vector in WLKIRC',-1)
         CALL OUTPUT(Q,1,1,1,NCART,1,NCART,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck wlkirm */
      SUBROUTINE WLKIRM(GRDDIA,EVAL,Q,NCART,RADIUS,DONWTN,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.0D0, SMALL = 1.0D-03, SMALL2 = 1.0D-10)
      LOGICAL DONWTN, INDXOK
      DIMENSION GRDDIA(NCART), EVAL(NCART), Q(NCART)
#include <trkoor.h>
#include <prkoor.h>
C
      DONWTN = .FALSE.
      CALL DZERO(Q,NCART)
C
C     First check index
C
      INDEX = 0
      DO 100 ICOOR = 1, NCART
          IF (EVAL(ICOOR) .LT. - SMALL) INDEX = INDEX + 1
  100 CONTINUE
C
      IF (INDEX .EQ. 0) THEN
C
C        OK index, test norm of Newton step
C
         QNRM  = D0
         DO 200 ICOOR = NPRREP(0) + 1, NCART
             IF (ABS(EVAL(ICOOR)) .GT. SMALL2) THEN
                STP  = - GRDDIA(ICOOR)/EVAL(ICOOR)
                QNRM = QNRM + STP*STP
             END IF
  200    CONTINUE
         QNRM = SQRT(QNRM)
C
C        OK index and norm of step, calculate Newton step
C
         IF (QNRM .LT. RADIUS) THEN
            DONWTN = .TRUE.
            DO 300 ICOOR = NPRREP(0) + 1, NCART
               IF (ABS(EVAL(ICOOR)) .GT. SMALL2) THEN
                  Q(ICOOR) = - GRDDIA(ICOOR)/EVAL(ICOOR)
               END IF
  300       CONTINUE
         END IF
      END IF
      RETURN
      END
C  /* Deck wlkirn */
      FUNCTION WLKIRN(GRDDIA,EVAL,NCART,RADIUS,IPRINT)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D100 = 100.0D0)
      PARAMETER (ITMAX = 1000, SMALL = 10.0D-6, NLOOP = 10)
      DIMENSION GRDDIA(NCART), EVAL(NCART)
C
C     calculates the timestep T corresponding to
C     an arclength DELTA using the Newton-Raphson
C     method:
C      DELTA = DELTA(ARCLENGTH) in [bohr * (amu)**1/2]
C
C     TMIN for bisection
C
      TMIN = D0
C
C     TMAX for bisection
C
      TMAX = D1
      ITER = 0
   10 CONTINUE
      ITER = ITER + 1
      IF (ITER .GT. ITMAX) THEN
         WRITE (LUPRI,'(/A,I5)')
     &      ' Maximum number of steps for TMAX in WLKIRC reached:',
     &        ITMAX
         CALL QUIT('Too many steps to determine TMAX in WLKIRC.')
      END IF
      STP = WLKIRS(TMAX,GRDDIA,EVAL,NCART)
      IF (STP .GT. RADIUS) THEN
         GO TO 50
      ELSE
         TMAX = D2*TMAX
         GO TO 10
      END IF
C
   50 CONTINUE
C
C     Bisection to determine TIME
C
      ITER = 0
 400  CONTINUE
         ITER= ITER+1
         TIME= (TMAX+TMIN)/D2
         STP = WLKIRS(TIME,GRDDIA,EVAL,NCART)
         IF (STP.GT.RADIUS) TMAX=TIME
         IF (STP.LT.RADIUS) TMIN=TIME
      IF (ABS(STP-RADIUS) .LT. RADIUS/D100) GO TO 401
      IF (ITER .LT. ITMAX) GO TO 400
C
         WRITE (LUPRI,'(//A,I5,A)') ' Bisection failed at after ',ITER,
     &                           ' iterations in WLKIRN'
         CALL QUIT('Bisection failed in WLKIRN')
C
 401  CONTINUE
      IF (IPRINT .GT. 10) THEN
         STP = WLKIRS(TIME,GRDDIA,EVAL,NCART)
         WRITE (LUPRI,'(/A,2F12.6)')
     &      ' Time and step size before Newton:', TIME, STP
      END IF
C
C     Loop over Newton iterations
C
      DO 100 ITER = 1, NLOOP
C
C        Current functional value
C
         FCUR = WLKIRS(TIME,GRDDIA,EVAL,NCART)
C
C        Derivative of function
C
         FDER = D0
         DO 200 I = 1, NCART
            IF(ABS(EVAL(I)) .GT. SMALL) THEN
               EXPT = EXP(-EVAL(I)*TIME)
               FDER = FDER - (EXPT - D1)*EXPT*(GRDDIA(I)**2)/EVAL(I)
            END IF
 200     CONTINUE
         FDER = FDER/FCUR
C
C        Newton step
C
         TIME = TIME - (FCUR - RADIUS)/FDER
  100 CONTINUE
      WLKIRN = TIME
      IF (IPRINT .GT. 10) THEN
         STP = WLKIRS(TIME,GRDDIA,EVAL,NCART)
         WRITE (LUPRI,'(/A,2F12.6)')
     &      ' Final time and step size in WLKIRN:', TIME, STP
      END IF
      RETURN
      END
C  /* Deck wlkirs */
      FUNCTION WLKIRS(T,GRDDIA,EVAL,NCART)
#include <implicit.h>
      PARAMETER (SMALL = 10.0D-8, D0 = 0.0D0, D1 = 1.0D0)
      DIMENSION GRDDIA(NCART), EVAL(NCART)
C
      WLKIRS = D0
      DO 10 I=1, NCART
         IF (ABS(EVAL(I)).LT.SMALL) THEN
            QSTEP = - GRDDIA(I)*T
         ELSE
            QSTEP = GRDDIA(I)*(EXP(-EVAL(I)*T) - D1)/EVAL(I)
         ENDIF
         WLKIRS = WLKIRS + QSTEP*QSTEP
 10   CONTINUE
      WLKIRS = SQRT(WLKIRS)
      RETURN
      END
C  /* Deck wlkirq */
      SUBROUTINE WLKIRQ(T,GRDDIA,EVAL,Q,NCART)
#include <implicit.h>
      PARAMETER (SMALL = 10.0D-8, D1 = 1.0D0)
      DIMENSION GRDDIA(NCART), EVAL(NCART), Q(NCART)
C
      DO 10 I=1, NCART
         IF (ABS(EVAL(I)).LT.SMALL) THEN
            Q(I) = - GRDDIA(I)*T
         ELSE
            Q(I) = GRDDIA(I)*(EXP(-EVAL(I)*T) - D1)/EVAL(I)
         ENDIF
 10   CONTINUE
      RETURN
      END
C  /* Deck wlksav */
      SUBROUTINE WLKSAV(ERGSEC,COOR,WORK,LWORK)
C
C     Save second-order MC prediction for found step
C
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (D1 = 1.0D0)
C
      DIMENSION COOR(3,*), WORK(LWORK)
C
C Used from common blocks:
C
C   CBIWLK : STEP, STEPCO(3,*), SCALCO(3,*), IPRWLK, ENRSEC
C   INFTAP : LUSIFC
C   NUCLEI : NUCDEP,CORD(3,*)
C   INFDIM : LBUFMA
C   INFVAR : NCONF,NWOPT,NVAR
C   INFORB : N2BAST,NCMOT
C   INFOPT : EMCSCF,EMCOLD,DEPRED,REJWMI,REJWMA, NREDL
C
#include <cbiwlk.h>
#include <inftap.h>
#include <nuclei.h>
#include <infdim.h>
#include <infvar.h>
#include <inforb.h>
#include <infopt.h>
#include <infpri.h>
C
C     kludge to avoid tranformation to natural orbitals
C     I must talk to HJJ / tuh Jun 92
#include <infinp.h>
      flag(15) = .false.
C
      CALL QENTER('WLKSAV')
C
      CALL ABAVAR(1,.FALSE.,0,WORK,LWORK)
C
C     A) Core allocation
C
      KCSTEP = 1
      KOSTEP = KCSTEP + NCONF
      KCMO   = KOSTEP + NWOPT
      KCREF  = KCMO   + NCMOT
      KWSAV  = KCREF
      KWRK   = KWSAV  + NCONF
      LWSAV  = LWORK  - KWSAV
      LWRK   = LWORK  - KWRK  + 1
      IF (KWRK  .GT. LWORK) CALL STOPIT('WLKSAV',' ',KWRK ,LWORK)
C
C     B) Calculate predicted reference CI vector for new geometry
C        Call abavar to make sure perturbation symmetry is one.
C
      REWIND LUSIFC
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC) EPOT,EMY,EACTIV,EMCSCF
      POTNUC = EPOT
      READ (LUSIFC)
      CALL READT (LUSIFC,NCMOT,WORK(KCMO))
      IF (NCONF .GT. 1) THEN
C
C        No first-order correction has been obtained in this version,
C        therefore set it to zero. Delete this line when WLKMC has been
C        reimplemented.
C
         CALL DZERO(WORK(KCSTEP),NCONF)
         CALL READT (LUSIFC,NCONF,WORK(KCREF))
C
         CALL DAXPY(NCONF,D1,WORK(KCREF),1,WORK(KCSTEP),1)
         CNRM = DDOT(NCONF,WORK(KCSTEP),1,WORK(KCSTEP),1)
         CNRM = D1/SQRT(CNRM)
         CALL DSCAL(NCONF,CNRM,WORK(KCSTEP),1)
         IF (IPRWLK .GT. 20) THEN
            CALL HEADER('Predicted CI vector:',1)
            WRITE (LUPRI,'(10F8.4)') (WORK(KCSTEP-1+I),I=1,NCONF)
         END IF
      ELSE
         WORK(KCSTEP) = D1
      END IF
C
C     C) In SIRSAV:
C     Save predicted reference CI vector and energy change,
C     rotate orbitals and save them
C
      IF (.NOT. REJECT) EMCOLD = EMCSCF
      DEPRED = ERGSEC
      REJWMI = REJMIN
      REJWMA = REJMAX
      NREDL  = 1
      CALL DCOPY(NCONF,WORK(KCSTEP),1,WORK(KWSAV),1)
C
C     No orbital rotation matrix has been generated in this version of
C     the program, therefore set it to zero. Delete this line when WLKMC
C     has been implemented.
C
      CALL DZERO(WORK(KOSTEP),NWOPT)
C
C     Orthogonalize new orbitals
C
      IF (.NOT.NOORTH) THEN
         CALL WLKORT(WORK(KCMO),WORK(KCMO),COOR,WORK(KWRK),LWRK,IPRWLK)
      END IF
C
      CALL GPOPEN(LUIT1,'SIRIUS.RST','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      CALL SIRSAV ('GEOSAVE',WORK(KCMO),DUMMY,DUMMY,DUMMY,
     &             WORK(KOSTEP),DUMMY,WORK(KWSAV),LWSAV)
      CALL GPCLOSE(LUIT1,'KEEP')
      CALL QEXIT('WLKSAV')
      RETURN
      END
C  /* Deck wlkder */
      FUNCTION WLKDER(GDDIA,HESDIA,RNU,NCORD,DUMMY)
C
C     Purpose:
C
C        Calculate derivative of step length function with
C        respect to RNU
C
C     WLKDER = STEP*STEP' / //STEP//
C
C        where
C
C     STEP = -GDDIA / (HESDIA+RNU)
C     STEP'=  GDDIA / ((HESDIA+RNU)**2)
C
#include <implicit.h>
      DIMENSION GDDIA(*),HESDIA(*)
      PARAMETER ( D0=0.0D0 )
      WLKDER=D0
      DEL   =D0
      DO 100 K = 1,NCORD
         DEL    = DEL    + (GDDIA(K)/(HESDIA(K)+RNU))**2
         WLKDER = WLKDER - (GDDIA(K)**2) / ((HESDIA(K)+RNU)**3)
 100  CONTINUE
      DEL    = SQRT(DEL)
      WLKDER = WLKDER/DEL
      RETURN
      END
C  /* Deck wlkstl */
      FUNCTION WLKSTL(GDDIA,HESDIA,RNU,NCORD,RTRUST)
C
C     Purpose:
C
C        Calculate step length at level shift RNU and
C        subtract RTRUST
C
C        WLKSTL = //STEP// - RTRUST
C
C        where
C
C        STEP = - GDDIA/(HESDIA+RNU)
C
#include <implicit.h>
      DIMENSION GDDIA(*),HESDIA(*)
      PARAMETER( D0=0.0D0 )
      STEP = D0
      DO 100 K=1,NCORD
         IF (GDDIA(K) .NE. D0) THEN
            STEPK = GDDIA(K) / (HESDIA(K)+RNU)
            STEP = STEP + STEPK*STEPK
         END IF
 100  CONTINUE
      WLKSTL = SQRT(STEP) - RTRUST
      RETURN
      END
C  /* Deck wlkbis */
      SUBROUTINE WLKBIS(XMAX,XMIN,XDET,GDDIA,HESDIA,TRUSTR,
     *                  NCORD,FUNCT,IFAIL)
C
C Purpose:
C
C  Use bisection to determine zero value of funct(.....)
C
C  Input:  XMAX ; FUNCT(..,XMAX,..) > 0
C          XMIN ; FUNCT(..,XMIN,..) < 0
C  Output: XDET ; FUNCT(..,XDET,..) = 0
C
#include <implicit.h>
      DIMENSION GDDIA(*),HESDIA(*)
      PARAMETER ( MAXIT=100 , D2=2.0D0 , D0=0.0D0 , DTEST=1.0D-7 )
      EXTERNAL FUNCT
      IFAIL = 0
      FMAX  = FUNCT(GDDIA,HESDIA,XMAX,NCORD,TRUSTR)
      FMIN  = FUNCT(GDDIA,HESDIA,XMIN,NCORD,TRUSTR)
      IF ( (FMAX.LT.D0) .OR. (FMIN.GT.D0) )  RETURN
C
      ITER = 0
 100  CONTINUE
         ITER= ITER+1
         XDET= (XMAX+XMIN)/D2
         FDET= FUNCT(GDDIA,HESDIA,XDET,NCORD,TRUSTR)
         IF (FDET.GT.D0) XMAX=XDET
         IF (FDET.LT.D0) XMIN=XDET
      IF (ABS(XMAX-XMIN).LT.DTEST) GO TO 101
      IF (ITER .LT. MAXIT) GO TO 100
C
      IFAIL=1
      RETURN
C
 101  CONTINUE
      IFAIL=5
      RETURN
      END
C  /* Deck wlkmc */
      SUBROUTINE WLKMC(REDFAC,STPCAR,WRK,LWRK,NCORD,MCOK,TEST,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
      PARAMETER (D1=1.0D0, XMXO=0.5D0, DP8=0.8D0)
C
      LOGICAL MCOK, TEST
      DIMENSION STPCAR(NCORD), WRK(LWRK)
C
#include <infvar.h>
C
      KSTPMC = 1
      KWRK1  = KSTPMC + NVAR
      LWRK1  = LWRK - KWRK1 + 1
      IF (KWRK1 .GT. LWRK) CALL STOPIT('WLKMC',' ',KWRK1,LWRK)
      CALL WLKMC1(NCORD,STPCAR,WRK(KSTPMC),XC,XO,X,WRK(KWRK1),LWRK1)
C
      WRITE (LUPRI,'(3(/A,F9.4))') ' Norm of predicted MCSCF step:',X,
     *                             ' configuration part:          ',XC,
     *                             ' orbital part:                ',XO
      IF (IPRINT .GT. 30) THEN
         CALL OUTPUT(WRK(KSTPMC),1,1,1,NVAR,1,NVAR,1,LUPRI)
         CALL HEADER('Predicted MCSCF step vector:',-1)
         IF (IPRINT .GT. 20) THEN
            WRITE (LUPRI,'(/,A,/)') ' CSF part of vector'
            CALL OUTPUT(WRK(KSTPMC),1,1,1,NCONF,1,NCONF,1,LUPRI)
         END IF
         WRITE (LUPRI,'(//,A)') ' Orbital part of vector'
         IF (IPRINT .GE. 20) THEN
            PRFAC = 0.0D0
         ELSE
            PRFAC = 0.1D0
         END IF
         CALL PRKAP(NWOPT,WRK(KSTPMC + NCONF),PRFAC,LUPRI)
      END IF
C
C     Check predicted orbital rotation:
C
      IF (TEST .AND. (XO .GT. XMXO)) THEN
         REDFAC = DP8*(XMXO/XO)
         MCOK   = .FALSE.
C        TRUSTR is undefined
C        WRITE (LUPRI,'(/A,F15.8,/A,F5.2)')
C    *         ' Trust radius is reduced to',TRUSTR,
C    *         ' because predicted orbital step is larger than',XMXO
         WRITE(LUPRI,'(A,F15.10,/,A,I3,A)')
     *      ' Trust radius is reduced since orbital step is too large.',
     *      ' Reduction of trust radius :   ',100*INT(D1 - REDFAC),'%'
      ELSE
         REDFAC = D1
         MCOK   = .TRUE.
      END IF
      RETURN
      END
C  /* Deck wlkmc1 */
      SUBROUTINE WLKMC1(NCORD,STPCOO,STPMC,XC,XO,X,WRK,LWRK)
C
C 13-Jun-1985/29-Nov-1985 hjaaj
C
C  a) calculate second order prediction for energy change.
C  b) if WFPRED, calculate first order prediction for change
C     in MC wave function.
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
      DIMENSION STPCOO(*),STPMC(*),WRK(*)
      LOGICAL OLDGD, OLDRD, GDMOD, GDOK
#include <iratdef.h>
      PARAMETER ( DP5 = 0.5D0 )
      PARAMETER ( THRFAC = 1.D-10 , FCSTEP = 1.0D-2 )
      PARAMETER ( FACLIN = 0.02D0, FACQUA = 0.1D0 )
C
C Used from common blocks:
C  ABAINF: GDALL
C  CBIWLK: WFPRED, ENRSEC, TRUSTR
C  INFTAP: LURDR,LRRDR
C  NUCLEI: DCORD(*,3), DCORGD(*,3)
C  ANRINF: THRNR
C  INFVAR: NCONF,NWOPT,NVAR
C  INFDIM: NVARMA
C
#include <abainf.h>
#include <inflin.h>
#include <cbiwlk.h>
#include <taymol.h>
#include <inftap.h>
#include <nuclei.h>
#include <dorps.h>
#include <anrinf.h>
#include <infvar.h>
#include <infdim.h>
#include <infpri.h>
C
      CHARACTER*1 CXYZ(3)
      DATA CXYZ/'X','Y','Z'/
C
C     Calculate first order prediction for change in MC wave function.
C
C GDMOD=.TRUE.  one linear set of equations is solved where coordinate
C               step multiplied on GD is used as modified gradient
C       .FALSE. RD vectors are used to get prediction about MC step
C
C Determine if RD vectors needed are calculated, otherwise
C set GDMOD .TRUE.
C
      GDMOD = .FALSE.
      GDOK  = .TRUE.
      STPLEN= SQRT( DDOT(NCORD,STPCOO,1,STPCOO,1) )
      IATCO = 0
      ISTART = 1
      DO 1000 IATOM = ISTART, NUCIND
         DO 1010 ICOOR = 1,3
            IATCO = IATCO + 1
            IF ( .NOT.DCORD(IATOM,ICOOR,1) .AND.
     *           ABS(STPCOO(IATCO)) .GT. FCSTEP*STPLEN ) THEN
               GDMOD = .TRUE.
               IF (.NOT. DCORGD(IATOM,ICOOR,1)) GDOK = .FALSE.
            END IF
 1010    CONTINUE
 1000 CONTINUE
      IF (GDALL) GDOK = .TRUE.
      IF (GDMOD .AND. .NOT. GDOK) THEN
         WRITE (LUPRI,'(/A)')
     *      ' GD vector(s) needed in WLKMC1 not available.'
         GO TO 8000
      END IF
C
C     Determine prediction to MC
C
      IF (GDMOD) THEN
C
C        Determine the MCSCF step by solving one set of linear
C        equations
C
         CALL HEADER('Solution of one set of linear equations for'//
     *      ' MC prediction',1)
         CALL GPOPEN(LUGDR,ABAGDR,'OLD','DIRECT',' ',IRAT*NVARMA,OLDGD)
         IF (.NOT. OLDGD) THEN
            WRITE (LUPRI,'(/A)') ' Failed to open LUGDR in WLKMC1.'
            GO TO 8000
         END IF
C
C        Construct modified gradient in WRK
C
         IF (IPRWLK .GE. 30) THEN
            WRITE (LUPRI,'(/A/)')
     *         ' *** Test of WLKMC1 calculating MC step'
         END IF
         CALL DZERO(WRK,NVAR)
         IOFF = 0
         DO 2100 I = 1,NCORD
            FAC = STPCOO(I)
            IF (ABS(FAC) .GT. THRFAC) THEN
               IREC = IOFF + I
               CALL READDX(LUGDR,IREC,IRAT*NVARMA,WRK(1+NVAR))
               CALL DAXPY(NVAR,FAC,WRK(1+NVAR),1,WRK,1)
               IF (IPRWLK .GT. 30) THEN
                  WRITE (LUPRI,'(//A,I3/A,F15.8/A)')
     *               ' Coordinate No. ',I,
     *               ' Step in this coordinate',FAC,
     *               ' GD for this coordinate'
                  WRITE (LUPRI,'(10F8.4)') (WRK(NVAR+J),J=1,NVAR)
               END IF
            END IF
 2100    CONTINUE
         IF (IPRWLK .GE. 30) THEN
            WRITE (LUPRI,'(//A)')
     *         ' GD vector for this step: sum(i) GD(i)*STPCOO(i)'
            WRITE (LUPRI,'(10F8.4)') (WRK(J),J=1,NVAR)
         END IF
C
C        set threshold for MC step vector convergence
C
         THRNRS = THRNR
         GDLEN = SQRT( DDOT(NVAR,WRK,1,WRK,1) )
         IF (STPLEN .GT. DP5*TRUSTR) THEN
            THRNR = MIN(THRNR, FACLIN * GDLEN)
            WRITE (LUPRI,'(//A,2(/A,1P,D10.2))')
     *         ' Global step - ',
     *         '   predicted MC gradient with frozen w.f. :',GDLEN,
     *         '   linear convergence threshold           :',THRNR
         ELSE
            THRNR = MIN(THRNR, FACQUA * GDLEN * GDLEN)
            WRITE (LUPRI,'(//A,2(/A,1P,D10.2))')
     *         ' Local step - ',
     *         '   predicted MC gradient with frozen w.f. :',GDLEN,
     *         '   quadratic convergence threshold        :',THRNR
         END IF
         IF (IPRWLK .GE. 6) THEN
            PTEST = MAX(1.D-10,GDLEN*0.01D0)
            WRITE (LUPRI,3010) PTEST
            DO 301 I = 1,NCONF
               GCI = WRK(I)
               IF (ABS(GCI) .GE. PTEST) WRITE (LUPRI,3011) I,GCI
  301       CONTINUE
            WRITE (LUPRI,3020)
            CALL PRKAP(NWOPT,WRK(1+NCONF),0.1D0,LUPRI)
 3010       FORMAT(//' Predicted configuration gradient after step',
     *              /' -------------------------------------------',
     *              /' (with frozen MC wave function)',
     *             //' Cutoff for print:',1P,D10.2,
     *             //' Configuration no.           value'
     *              /' -----------------           -----')
 3011       FORMAT(I16,F20.10)
 3020       FORMAT(//' Predicted orbital gradient after step',
     *              /' -------------------------------------',
     *              /' (with frozen MC wave function)')
         END IF
C
         CALL RESMC(.TRUE.,WRK,LWRK)
         THRNR  = THRNRS
C
         CALL DCOPY(NVAR,WRK,1,STPMC,1)
      ELSE
C
C        Determine the MCSCF step using the RD vectors
C
         WRITE (LUPRI,'(//A/)')
     *      'MC prediction calculated from available RD vectors.'
C
C        Open LURDR file
C
         CALL GPOPEN(LURDR,ABARDR,'OLD','DIRECT',' ',IRAT*NVARMA,OLDRD)
         CALL DZERO(STPMC,NVAR)
         IREC   = 0
         ISTART = 1
         ISTPCO = 0
         DO 2210 IATOM = ISTART,NUCIND
            DO 2200 ICOOR = 1,3
               IREC   = IREC + 1
               ISTPCO = ISTPCO + 1
               FAC = STPCOO(ISTPCO)
               IF ( ABS(FAC) .GT. THRFAC ) THEN
                  IF ( DCORD(IATOM,ICOOR,1) ) THEN
                     CALL READDX(LURDR,IREC,IRAT*NVAR,WRK)
                     CALL DAXPY(NVAR,FAC,WRK,1,STPMC,1)
                  ELSE
                     WRITE (LUPRI,'(/A,I3,A,A)')
     *               ' Step coordinate',IATOM,CXYZ(ICOOR),
     *               ' is neglected for MC prediction.'
                  END IF
               END IF
 2200       CONTINUE
 2210    CONTINUE
      END IF
      XC = DDOT(NCONF,STPMC,1,STPMC,1)
      XO = DDOT(NWOPT,STPMC(1+NCONF),1,STPMC(1+NCONF),1)
      X  = SQRT(XC+XO)
      XO = SQRT(XO)
      XC = SQRT(XC)
C
C close LUGDR/LURDR
C
      IF (GDMOD) THEN
         CALL GPCLOSE(LUGDR,'DELETE')
      ELSE
         CALL GPCLOSE(LURDR,'DELETE')
      END IF
C
C     End of WLKMC1
C
      RETURN
C
 8000 CONTINUE
      WFPRED = .FALSE.
      WRITE (LUPRI,'(/A/A/)') ' *** WARNING *** insufficient'//
     *   ' information for MC prediction in WALK,',
     *   ' MC prediction not evaluated.'
      RETURN
      END
C  /* Deck wlktay */
      FUNCTION WLKTAY(GRAD,HESS,STEP)
#include <implicit.h>
      PARAMETER (DP5 = 0.5D0)
      WLKTAY = GRAD*STEP + DP5*STEP*HESS*STEP
      RETURN
      END
C  /* Deck wlkdan */
      SUBROUTINE WLKDAN(PMCARX,NCART,WORK,LWORK,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
#include <codata.h>
      PARAMETER (D0 = 0.D0, D1 = 1.D0)
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
      DIMENSION PMCARX(NCART), WORK(LWORK)
#include <cbiwlk.h>
#include <nuclei.h>
#include <symmet.h>
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('PMCARX in WLKDAN',-1)
         CALL OUTPUT(PMCARX,1,1,1,NCART,1,NCART,1,LUPRI)
      END IF
C
      NT = NUCDEP
      NA = 0
      NB = 0
      DO 100 I = 1, NUCDEP
         IF (IPART(I).EQ.1) NA = NA + 1
         IF (IPART(I).EQ.2) NB = NB + 1
  100 CONTINUE
C
      KTMASS = 1
      KTCOR  = KTMASS +   NT
      KTVEL  = KTCOR  + 3*NT
      KTCORR = KTVEL  + 3*NT
      KTVELR = KTCORR + 3*NT
      KTVLRT = KTVELR + 3*NT
C
      KAMASS = KTVLRT + 3*NT
      KACOR  = KAMASS +   NA
      KAVEL  = KACOR  + 3*NA
      KACORR = KAVEL  + 3*NA
      KAVELR = KACORR + 3*NA
      KAVLRT = KAVELR + 3*NA
C
      KBMASS = KAVLRT + 3*NA
      KBCOR  = KBMASS +   NB
      KBVEL  = KBCOR  + 3*NB
      KBCORR = KBVEL  + 3*NB
      KBVELR = KBCORR + 3*NB
      KBVLRT = KBVELR + 3*NB
      KWRK   = KBVLRT + 3*NB
      IF (KWRK .GT. LWORK) CALL STOPIT('WLKDAN',' ',KWRK,LWORK)
      LWRK = LWORK - KWRK + 1
      CALL WLKDA1(PMCARX,NCART,NT,NA,NB,WORK(KWRK),LWRK,IPRINT,
     &            WORK(KTMASS),WORK(KTCOR),WORK(KTVEL),WORK(KTCORR),
     &            WORK(KTVELR),WORK(KTVLRT),
     &            WORK(KAMASS),WORK(KACOR),WORK(KAVEL),WORK(KACORR),
     &            WORK(KAVELR),WORK(KAVLRT),
     &            WORK(KBMASS),WORK(KBCOR),WORK(KBVEL),WORK(KBCORR),
     &            WORK(KBVELR),WORK(KBVLRT))
      RETURN
      END
C  /* Deck wlkda1 */
      SUBROUTINE WLKDA1(PMCARX,NCART,NT,NA,NB,WORK,LWORK,IPRINT,
     &                  TMASS,TCOR,TVEL,TCORR,TVELR,TVLROT,
     &                  AMASS,ACOR,AVEL,ACORR,AVELR,AVLROT,
     &                  BMASS,BCOR,BVEL,BCORR,BVELR,BVLROT)
C
#include <implicit.h>
#include <priunit.h>
#include <codata.h>
      PARAMETER (D0 = 0.D0, DP5 = 0.5D0, D1 = 1.D0, D2 = 2.0D0,
     &           D100 = 100.0D0)
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
      INTEGER X, Y, Z
      LOGICAL TPLAN, APLAN, BPLAN, TLIN, ALIN, BLIN
C
      DIMENSION PMCARX(*),
     &          TMASS(NT), TCOR(NT,3), TVEL(NT,3),
     &          TCORR(NT,3), TVELR(NT,3), TVLROT(NT,3),
     &          TCMXYZ(3), TVLXYZ(3), TVBANG(3),
     &          TLXYZ(3), TOMEGA(3), TINERT(3,3),
     &          AMASS(NA), ACOR(NA,3), AVEL(NA,3),
     &          ACORR(NA,3), AVELR(NA,3), AVLROT(NA,3),
     &          ACMXYZ(3), AVLXYZ(3), AVBANG(3),
     &          ALXYZ(3), AOMEGA(3), AINERT(3,3),
     &          BMASS(NB), BCOR(NB,3), BVEL(NB,3),
     &          BCORR(NB,3), BVELR(NB,3), BVLROT(NB,3),
     &          BCMXYZ(3), BVLXYZ(3), BVBANG(3),
     &          BLXYZ(3), BOMEGA(3), BINERT(3,3),
     &          VECTOR(3)
#include <cbiwlk.h>
#include <nuclei.h>
#include <symmet.h>
C
#include <ibtfun.h>
c     NEXT1(I) = MOD(I,3) + 1
c     NEXT2(I) = NEXT1(NEXT1(I))
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('PMCARX in WLKDA1',-1)
         CALL OUTPUT(PMCARX,1,1,1,NCART,1,NCART,1,LUPRI)
      END IF
C
      IATOM = 0
      DO 100 ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
         FAC    = D1/SQRT(FMULT(MULCNT))
         NATTYP = IZATOM(ICENT)
         NISOTP = ISOTOP(ICENT)
         ATMAS  = XFAMU*DISOTP(NATTYP,NISOTP,'MASS')
         DO 110 ISYMOP = 0, MAXOPR
            IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
               IATOM = IATOM + 1
               TMASS(IATOM) = ATMAS
               DO 120 ICOOR = 1, 3
                  ISCOOR = IPTCNT(3*(ICENT-1)+ICOOR,0,1)
                  IF (ISCOOR .NE. 0) THEN
                     SGN = PT(IBTAND(ISYMAX(ICOOR,1),ISYMOP))
                     TCOR(IATOM,ICOOR) = SGN*CORD(ICOOR,ICENT)
                     TVEL(IATOM,ICOOR) = SGN*FAC*PMCARX(ISCOOR)/ATMAS
                  END IF
  120          CONTINUE
            END IF
  110    CONTINUE
  100 CONTINUE
C
C     Print
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('TMASS in WLKDA1',-1)
         CALL OUTPUT(TMASS,1,1,1,NT,1,1,1,LUPRI)
         CALL HEADER('TCOR in WLKDA1',-1)
         CALL OUTPUT(TCOR,1,NT,1,3,NT,3,1,LUPRI)
         CALL HEADER('TVEL in WLKDA1',-1)
         CALL OUTPUT(TVEL,1,NT,1,3,NT,3,1,LUPRI)
      END IF
C
C     Analysis of total system
C
      CALL WLKDFR('T',NT,TMASS,TCOR,TVEL,TCORR,TVELR,TVLROT,
     &            TCMXYZ,TVLXYZ,TLXYZ,TOMEGA,TVBANG,TINERT,
     &            TTMASS,ETOTT,ETRAT,EROTT,EVIBT,ECORT,
     &            TPLAN,TLIN,IPRINT)
C
C     Fragments A and B
C     =================
C
      IF (NA .GT. 0 .AND. NB .GT. 0) THEN
C
C        Fragment coordinates and velocities
C        ===================================
C
         IA = 0
         IB = 0
         DO 300 I = 1, NT
            IF (IPART(I).EQ.1) THEN
               IA = IA + 1
               AMASS(IA) = TMASS(I)
            ELSE
               IB = IB + 1
               BMASS(IB) = TMASS(I)
            END IF
            DO 310 X = 1, 3
               IF (IPART(I).EQ.1) THEN
                  ACOR(IA,X) = TCORR(I,X)
                  AVEL(IA,X) = TVLROT(I,X)
               ELSE
                  BCOR(IB,X) = TCORR(I,X)
                  BVEL(IB,X) = TVLROT(I,X)
               END IF
  310       CONTINUE
  300    CONTINUE
C
         CALL WLKDFR('A',NA,AMASS,ACOR,AVEL,ACORR,AVELR,AVLROT,
     &               ACMXYZ,AVLXYZ,ALXYZ,AOMEGA,AVBANG,AINERT,
     &               TAMASS,ETOTA,ETRAA,EROTA,EVIBA,ECORA,
     &               APLAN,ALIN,IPRINT)
C
         CALL WLKDFR('B',NB,BMASS,BCOR,BVEL,BCORR,BVELR,BVLROT,
     &               BCMXYZ,BVLXYZ,BLXYZ,BOMEGA,BVBANG,BINERT,
     &               TBMASS,ETOTB,ETRAB,EROTB,EVIBB,ECORB,
     &               BPLAN,BLIN,IPRINT)
C
C        Relative translational energy
C        =============================
C
         VRELX  = AVLXYZ(1) - BVLXYZ(1)
         VRELY  = AVLXYZ(2) - BVLXYZ(2)
         VRELZ  = AVLXYZ(3) - BVLXYZ(3)
         RMASS  = ((TAMASS*TBMASS)/(TAMASS + TBMASS))
         ERLTRA = DP5*RMASS*(VRELX*VRELX + VRELY*VRELY + VRELZ*VRELZ)
         CALL AROUND('Analysis of the system and its fragments')
         IF (TPLAN) WRITE (LUPRI,'(1X,A)') ' Total system is planar.'
         IF (TLIN ) WRITE (LUPRI,'(1X,A)') ' Total system is linear.'
         IF (APLAN) WRITE (LUPRI,'(1X,A)') ' Fragment A is planar.'
         IF (ALIN ) WRITE (LUPRI,'(1X,A)') ' Fragment A is linear.'
         IF (BPLAN) WRITE (LUPRI,'(1X,A)') ' Fragment B is planar.'
         IF (BLIN ) WRITE (LUPRI,'(1X,A)') ' Fragment B is linear.'
         CALL HEADER('Positions of centers of mass',-1)
         WRITE (LUPRI,'(1X,A,3F12.6)')  ' CM of total system:',
     &        TCMXYZ(1), TCMXYZ(2), TCMXYZ(3)
         WRITE (LUPRI,'(1X,A,3F12.6)')  ' CM of fragment A  :',
     &        ACMXYZ(1), ACMXYZ(2), ACMXYZ(3)
         WRITE (LUPRI,'(1X,A,3F12.6)')  ' CM of fragment B  :',
     &        BCMXYZ(1), BCMXYZ(2), BCMXYZ(3)
         CALL HEADER('CM velocities (components and total)',-1)
         WRITE (LUPRI,'(1X,A,4F12.6)') ' CM velocity of total system:',
     &      TVLXYZ(1), TVLXYZ(2), TVLXYZ(3),
     &      SQRT(DDOT(3,TVLXYZ,1,TVLXYZ,1))
         WRITE (LUPRI,'(1X,A,4F12.6)') ' CM velocity of fragment A  :',
     &      AVLXYZ(1), AVLXYZ(2), AVLXYZ(3),
     &      SQRT(DDOT(3,AVLXYZ,1,AVLXYZ,1))
         WRITE (LUPRI,'(1X,A,4F12.6)') ' CM velocity of fragment B  :',
     &      BVLXYZ(1), BVLXYZ(2), BVLXYZ(3),
     &      SQRT(DDOT(3,BVLXYZ,1,BVLXYZ,1))
         CALL HEADER('Angular momenta in CM coordinate systems',-1)
         WRITE (LUPRI,'(1X,A,3F12.6)')
     &       ' Angular momemtum of total system:',
     &        TLXYZ(1), TLXYZ(2), TLXYZ(3)
         WRITE (LUPRI,'(1X,A,3F12.6)')
     &       ' Angular momemtum of fragment A  :',
     &         ALXYZ(1), ALXYZ(2), ALXYZ(3)
         WRITE (LUPRI,'(1X,A,3F12.6)')
     &       ' Angular momentum of fragment B  :',
     &         BLXYZ(1), BLXYZ(2), BLXYZ(3)
         CALL HEADER('Angular velocities:',-1)
         WRITE (LUPRI,'(1X,A,3F12.6)')
     &       ' Angular velocity of total system:',
     &         TOMEGA(1), TOMEGA(2), TOMEGA(3)
         WRITE (LUPRI,'(1X,A,3F12.6)')
     &       ' Angular velocity of fragment A  :',
     &         AOMEGA(1), AOMEGA(2), AOMEGA(3)
         WRITE (LUPRI,'(1X,A,3F12.6)')
     &       ' Angular velocity of fragment B  :',
     &         BOMEGA(1), BOMEGA(2), BOMEGA(3)
         CALL HEADER('Vibrational angular momenta:',-1)
         WRITE (LUPRI,'(1X,A,3F12.6)')
     &       ' Vib. ang. momentum of total system:',
     &         TVBANG(1), TVBANG(2), TVBANG(3)
         WRITE (LUPRI,'(1X,A,3F12.6)')
     &       ' Vib. ang. momentum of fragment A  :',
     &         AVBANG(1), AVBANG(2), AVBANG(3)
         WRITE (LUPRI,'(1X,A,3F12.6)')
     &       ' Vib. ang. momentum of fragment B  :',
     &         BVBANG(1), BVBANG(2), BVBANG(3)
         CALL AROUND('Breakdown of kinetic energy contributions')
         CALL HEADER('Fragment breakdown',-1)
         ETOTCM = ETRAT + EROTT + ECORT
         TOTAL = D0
         CALL WLKDPR('CM            :',ETOTCM,ETOTT,TOTAL)
         CALL WLKDPR('fragment A    :',ETOTA, ETOTT,TOTAL)
         CALL WLKDPR('fragment B    :',ETOTB, ETOTT,TOTAL)
         WRITE (LUPRI,'(1X)')
         CALL WLKDPR('total         :',TOTAL,ETOTT,TOTAL)
C
         CALL HEADER('Modes of motion breakdown',-1)
         TOTAL = D0
         CALL WLKDPR('vibration     :',EVIBA+EVIBB,ETOTT,TOTAL)
         CALL WLKDPR('rotation      :',EROTT+EROTA+EROTB,ETOTT,TOTAL)
         CALL WLKDPR('Coriolis      :',ECORT+ECORA+ECORB,ETOTT,TOTAL)
         CALL WLKDPR('translation   :',ETRAT+ETRAA+ETRAB,ETOTT,TOTAL)
         WRITE (LUPRI,'(1X)')
         CALL WLKDPR('total         :',TOTAL,ETOTT,TOTAL)
C
         CALL HEADER('Complete breakdown',-1)
         TOTAL = D0
         CALL WLKDPR('rotation    CM:',EROTT,ETOTT,TOTAL)
         CALL WLKDPR('Coriolis    CM:',ECORT,ETOTT,TOTAL)
         CALL WLKDPR('translation CM:',ETRAT,ETOTT,TOTAL)
         WRITE (LUPRI,'(1X)')
         CALL WLKDPR('vibration    A:',EVIBA,ETOTT,TOTAL)
         CALL WLKDPR('rotation     A:',EROTA,ETOTT,TOTAL)
         CALL WLKDPR('Coriolis     A:',ECORA,ETOTT,TOTAL)
         CALL WLKDPR('translation  A:',ETRAA,ETOTT,TOTAL)
         WRITE (LUPRI,'(1X)')
         CALL WLKDPR('vibration    B:',EVIBB,ETOTT,TOTAL)
         CALL WLKDPR('rotation     B:',EROTB,ETOTT,TOTAL)
         CALL WLKDPR('Coriolis     B:',ECORB,ETOTT,TOTAL)
         CALL WLKDPR('translation  B:',ETRAB,ETOTT,TOTAL)
         WRITE (LUPRI,'(1X)')
         CALL WLKDPR('total         :',TOTAL,ETOTT,TOTAL)
C
         WRITE (LUPRI,'(/1X,A,F12.6)')
     &      'Relative translational energy between fragments:',ERLTRA
      ELSE
         CALL AROUND('Breakdown of kinetic energy contributions')
         TOTAL = D0
         CALL WLKDPR('translational :',ETRAT,ETOTT,TOTAL)
         CALL WLKDPR('rotational    :',EROTT,ETOTT,TOTAL)
         CALL WLKDPR('vibration     :',EVIBT,ETOTT,TOTAL)
         CALL WLKDPR('Coriolis      :',ECORT,ETOTT,TOTAL)
         WRITE (LUPRI,'(1X)')
         CALL WLKDPR('total         :',TOTAL,ETOTT,TOTAL)
      END IF
C
      RETURN
      END
C  /* Deck wlkdfr */
      SUBROUTINE WLKDFR(FRAG,NATOMS,AMASS,CORABS,VELABS,CORREL,VELREL,
     &                  VELROT,CMXYZ,VELXYZ,ANGXYZ,OMEGA,VIBANG,TINERT,
     &                  TAMASS,ETOT,ETRA,EROT,EVIB,ECOR,PLANAR,LINEAR,
     &                  IPRINT)
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
      PARAMETER (D0 = 0.D0, DP5 = 0.5D0, D1 = 1.D0, D2 = 2.0D0)
#include <mxcent.h>
      INTEGER X, Y, Z
      CHARACTER FRAG*1
      LOGICAL PLANAR, LINEAR
      DIMENSION AMASS(NATOMS), CORABS(NATOMS,3), VELABS(NATOMS,3),
     &          CORREL(NATOMS,3), VELREL(NATOMS,3), VELROT(NATOMS,3),
     &          CMXYZ(3), VELXYZ(3), VIBANG(3),
     &          ANGXYZ(3), OMEGA(3), TINERT(3,3),
     &          VECTOR(3)
C
C     Total kinetic energy from Cartesian velocities
C
      TKIN = D0
      DO 100 X = 1, 3
      DO 100 I = 1, NATOMS
         TKIN = TKIN + DP5*AMASS(I)*VELABS(I,X)*VELABS(I,X)
 100  CONTINUE
C
C     Center of mass
C
      CALL WLKDCM(AMASS,CORABS,NATOMS,TAMASS,CMXYZ)
      CALL WLKDCM(AMASS,VELABS,NATOMS,TAMASS,VELXYZ)
C
C     Cartesian coordinates and velocities relative to CM
C
      DO 200 X = 1, 3
      DO 200 I = 1, NATOMS
         CORREL(I,X) = CORABS(I,X) - CMXYZ(X)
         VELREL(I,X) = VELABS(I,X) - VELXYZ(X)
  200 CONTINUE
C
C     Angular momentum ANGXYZ
C
      CALL WLKDAM(CORREL,VELREL,AMASS,NATOMS,ANGXYZ)
C
C     Angular velocity OMEGA
C
      CALL WLKDIN(CORREL,AMASS,NATOMS,ANGXYZ,TINERT,OMEGA,DUMMY,DUMMY,
     &            .FALSE.,PLANAR,LINEAR)
C
C     Velocities in rotating system
C
      DO 300 I = 1, NATOMS
         CALL WLKCRS(OMEGA,1,CORREL(I,1),NATOMS,VECTOR,1)
         DO 310 X = 1, 3
            VELROT(I,X) = VELREL(I,X) - VECTOR(X)
  310    CONTINUE
  300 CONTINUE
C
C     Vibrational angular momentum
C
      CALL DZERO(VIBANG,3)
      DO 400 I = 1, NATOMS
         CALL WLKCRS(CORREL(I,1),NATOMS,VELROT(I,1),NATOMS,VECTOR,1)
         DO 410 X = 1, 3
            VIBANG(X) = VIBANG(X) + AMASS(I)*VECTOR(X)
  410    CONTINUE
  400 CONTINUE
C
C     Calculation of energy contributions
C     ===================================
C
C     - translational energy
C
      ETRA = DP5*TAMASS*DDOT(3,VELXYZ,1,VELXYZ,1)
C
C     - rotational energy
C
      EROT = D0
      DO 500 I = 1, NATOMS
         CALL WLKCRS(OMEGA,1,CORREL(I,1),NATOMS,VECTOR,1)
         EROT = EROT + DP5*AMASS(I)*DDOT(3,VECTOR,1,VECTOR,1)
  500 CONTINUE
C
C     - vibrational energy
C
      EVIB = D0
      DO 600 I = 1, NATOMS
         EVIB = EVIB + DP5*AMASS(I)*DDOT(3,VELROT(I,1),NATOMS,
     &                                     VELROT(I,1),NATOMS)
  600 CONTINUE
C
C     - Coriolis energy
C
      ECOR = DDOT(3,OMEGA,1,VIBANG,1)
C
C     Total kinetic energy from individual contributions
C
      ETOT = ETRA + EROT + EVIB + ECOR
C
      IF (IPRINT .GT. 10) THEN
         CALL AROUND('Dynamic analysis for fragment '//FRAG)
         WRITE (LUPRI,'(1X,A,I5)')
     &      ' Number of atoms in fragment',NATOMS
         WRITE (LUPRI,'(1X,A,F12.6)')
     &      ' Total mass of fragment ',TAMASS
         CALL HEADER('Masses of fragment '//FRAG,-1)
         CALL OUTPUT(AMASS,1,NATOMS,1,1,NATOMS,1,1,LUPRI)
         CALL HEADER('Coordinates of fragment '//FRAG,-1)
         CALL OUTPUT(CORABS,1,NATOMS,1,3,NATOMS,3,1,LUPRI)
         CALL HEADER('Velocities of fragment '//FRAG,-1)
         CALL OUTPUT(VELABS,1,NATOMS,1,3,NATOMS,3,1,LUPRI)
         CALL HEADER('Energy contributions for fragment'//FRAG,-1)
         WRITE (LUPRI,'(1X,A,2F12.6)')
     &      ' Total kinetic energy of fragment:   ',ETOT,TKIN
         WRITE (LUPRI,'(1X,A,F12.6)')
     &      ' Translational energy of fragment:   ',ETRA
         IF (NATOMS .GT. 1) WRITE (LUPRI,'(1X,A,F12.6)')
     &      ' Rotational    energy of fragment:   ',EROT
         IF (NATOMS .GT. 1) WRITE (LUPRI,'(1X,A,F12.6)')
     &      ' Internal kinetic energy of fragment:',EVIB
         IF (NATOMS .GT. 1) WRITE (LUPRI,'(1X,A,F12.6)')
     &      ' Coriolis energy of fragment:        ',ECOR
         WRITE (LUPRI,'(/,1X,A,3F12.6)')' Position of center of mass ',
     &      CMXYZ(1),CMXYZ(2),CMXYZ(3)
         WRITE (LUPRI,'(1X,A,3F12.6)')' Velocity of center of mass ',
     &      VELXYZ(1),VELXYZ(2),VELXYZ(3)
         WRITE (LUPRI,'(1X,A,3F12.6)')' Angular momentum of fragment',
     &      ANGXYZ(1),ANGXYZ(2),ANGXYZ(3)
         WRITE (LUPRI,'(1X,A,3F12.6)')' Angular velocity of fragment',
     &      OMEGA(1),OMEGA(2),OMEGA(3)
         IF (NATOMS .GT. 1) THEN
            CALL HEADER('Moments of inertia',-1)
            CALL OUTPUT(TINERT,1,3,1,3,3,3,1,LUPRI)
         END IF
         CALL HEADER('Coordinates relative to center of mass',-1)
         CALL OUTPUT(CORREL,1,NATOMS,1,3,NATOMS,3,1,LUPRI)
         CALL HEADER('Velocities relative to center of mass',-1)
         CALL OUTPUT(VELREL,1,NATOMS,1,3,NATOMS,3,1,LUPRI)
         CALL HEADER('Velocities in rotating mass center system',-1)
         CALL OUTPUT(VELROT,1,NATOMS,1,3,NATOMS,3,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck wlkdpr */
      SUBROUTINE WLKDPR(LABEL,PART,TOTAL,ACCUM)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1 = 1.0D0, D100 = 100.0D0)
      CHARACTER LABEL*15
      FRACT = D100*PART/TOTAL
      IF (ABS(FRACT) .GE. D1/(D100*D100)) THEN
         WRITE (LUPRI,'(20X,A,F12.6,2X,F6.2,A)') LABEL, PART, FRACT,' %'
      END IF
      ACCUM = ACCUM + PART
      RETURN
      END
C  /* Deck wlkdcm */
      SUBROUTINE WLKDCM(TMASS,TCOR,N,TOTMAS,CMXYZ)
#include <implicit.h>
      PARAMETER (D0 = 0.0D0)
      INTEGER X
      DIMENSION TMASS(N), TCOR(N,3), CMXYZ(3)
C
C     Total mass
C
      TOTMAS = DSUM(N,TMASS,1)
      DO 100 X = 1, 3
         CMXYZ(X) = D0
         DO 200 I = 1, N
            CMXYZ(X) = CMXYZ(X) + TMASS(I)*TCOR(I,X)
  200    CONTINUE
         CMXYZ(X) = CMXYZ(X)/TOTMAS
  100 CONTINUE
      RETURN
      END
C  /* Deck wlkdam */
      SUBROUTINE WLKDAM(COR,VEL,TMASS,N,ANGMOM)
#include <implicit.h>
      PARAMETER (D0 = 0.0D0)
      INTEGER X
      DIMENSION VEL(N,3), COR(N,3), ANGMOM(3), TMASS(N)
      NEXT1(I) = MOD(I,3) + 1
      NEXT2(I) = NEXT1(NEXT1(I))
C
C     Angular momentum
C
      DO 100 X = 1, 3
         ANGMOM(X) = D0
         DO 200 I = 1, N
            ANGMOM(X) = ANGMOM(X) + TMASS(I)*
     &         (COR(I,NEXT1(X))*VEL(I,NEXT2(X)) -
     &          COR(I,NEXT2(X))*VEL(I,NEXT1(X)))
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck wlkdin */
      SUBROUTINE WLKDIN(COR,TMASS,N,ANGMOM,TINERT,OMEGA,CEPVAL,
     &                  CEPVEC,DOCOPY,PLANAR,LINEAR)
#include <implicit.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, TSTLIN = 1.D-4)
      INTEGER X, Y, Z
      LOGICAL PLANAR, LINEAR, DOCOPY
      DIMENSION COR(N,3), TINERT(3,3), TINVER(3,3), TMASS(N),
     &          ANGMOM(3), OMEGA(3)
      DIMENSION EIGVAL(3), EIGVEC(3,3), AMAT(6), AWRK(3), IWRK(3),
     &          CEPVAL(*), CEPVEC(*)
C
      CALL DZERO(TINERT,9)
      CALL DZERO(OMEGA,3)
C
      IF (N.EQ.1) RETURN
C
C     Moments of inertia tensor
C
      DO 100 I = 1, N
         R2 = TMASS(I)*DDOT(3,COR(I,1),N,COR(I,1),N)
         DO 110 X = 1, 3
            TINERT(X,X) =  TINERT(X,X) + R2
            DO 120 Y = 1, 3
               TINERT(X,Y) = TINERT(X,Y) - TMASS(I)*COR(I,X)*COR(I,Y)
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C     Invert TINERT
C
      AMAT(1) = TINERT(1,1)
      AMAT(2) = TINERT(2,1)
      AMAT(3) = TINERT(2,2)
      AMAT(4) = TINERT(3,1)
      AMAT(5) = TINERT(3,2)
      AMAT(6) = TINERT(3,3)
      CALL DUNIT(EIGVEC,3)
      CALL JACO(AMAT,EIGVEC,3,3,3,AWRK,IWRK)
      EIGVAL(1) = AMAT(1)
      EIGVAL(2) = AMAT(3)
      EIGVAL(3) = AMAT(6)
      CALL ORDER(EIGVEC,EIGVAL,3,3)
C
      IF ( ABS(EIGVAL(3)-EIGVAL(2)-EIGVAL(1)) .LT. TSTLIN) THEN
         PLANAR = .TRUE.
      ELSE
         PLANAR = .FALSE.
      END IF
      IF (EIGVAL(1) .LT. TSTLIN) THEN
         LINEAR = .TRUE.
         PLANAR = .FALSE.
         EIGVAL(1) = D0
      ELSE
         LINEAR = .FALSE.
         EIGVAL(1) = D1/EIGVAL(1)
      END IF
      EIGVAL(2) = D1/EIGVAL(2)
      EIGVAL(3) = D1/EIGVAL(3)
C
C     Inverse moments of inertia
C
      CALL DZERO(TINVER,9)
      DO 200 X = 1, 3
      DO 200 Y = 1, 3
      DO 200 I = 1, 3
         TINVER(X,Y) = TINVER(X,Y) + EIGVEC(X,I)*EIGVEC(Y,I)*EIGVAL(I)
  200 CONTINUE
C
      DO 300 X = 1, 3
      DO 300 Y = 1, 3
         OMEGA(X) = TINVER(X,Y)*ANGMOM(Y)
  300 CONTINUE
C
      IF (DOCOPY) THEN
         CALL DCOPY(3,EIGVAL,1,CEPVAL,1)
         CALL DCOPY(9,EIGVEC,1,CEPVEC,1)
      END IF
      RETURN
      END
C  /* Deck wlkcrs */
      SUBROUTINE WLKCRS(VEC1,N1,VEC2,N2,VEC3,N3)
#include <implicit.h>
      DIMENSION VEC1(N1,3), VEC2(N2,3), VEC3(N3,3)
      VEC3(1,1) = VEC1(1,2)*VEC2(1,3) - VEC1(1,3)*VEC2(1,2)
      VEC3(1,2) = VEC1(1,3)*VEC2(1,1) - VEC1(1,1)*VEC2(1,3)
      VEC3(1,3) = VEC1(1,1)*VEC2(1,2) - VEC1(1,2)*VEC2(1,1)
      RETURN
      END
C  /* Deck wlkroa */
      SUBROUTINE WLKROA(ROAAFD,ROAGND,ROAGLD,ROAAD,ROAAFU,ROAGNU,ROAGLU,
     &                  ROAAU,ROAAFP,ROAAP,ROAGLP,ROAGNP,WORK,LWORK,
     &                  NXYZ,NCART,IPOINT,DISPLC,IPRINT)
C
C     Numerical differentiation for Raman properties
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <codata.h>
C
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D2=2.0D0,DE3=1.0D3)
      PARAMETER (XFRAIN = XTANG**6,
     &           XFOAIN = ALPHAC*XTANG**6*1.0D6)
C
      LOGICAL FINAL, DODIP
      CHARACTER*9 TYPE
      INTEGER X, F
      DIMENSION WORK(LWORK)
C
#include <abainf.h>
#include <cbilnr.h>
#include <nuclei.h>
#include <symmet.h>
C
      DIMENSION ROAAFU(3,3,MXFR),  ROAAFP(2,3,3,MXFR),
     &          ROAAFD(3,3,MXFR,MXCOOR),
     &          ROAGNU(3,3,MXFR),  ROAGNP(2,3,3,MXFR),
     &          ROAGND(3,3,MXFR,MXCOOR),
     &          ROAGLU(3,3,MXFR),  ROAGLP(2,3,3,MXFR),
     &          ROAGLD(3,3,MXFR,MXCOOR),
     &          ROAAU(3,3,3,MXFR), ROAAP(2,3,3,3,MXFR),
     &          ROAAD(3,3,3,MXFR,MXCOOR)
C
#include <chrsgn.h>
C
      FINAL = IPOINT .EQ. 2*NCART
      IF (IPOINT .EQ. 0) THEN
         CALL DCOPY (9*MXFR,ROAAFP,2,ROAAFU,1)
         CALL DCOPY (9*MXFR,ROAGNP,2,ROAGNU,1)
         CALL DCOPY (9*MXFR,ROAGLP,2,ROAGLU,1)
         CALL DCOPY (27*MXFR,ROAAP,2,ROAAU,1)
         CALL DZERO(ROAAFD, 9*MXFR*MXCOOR)
         CALL DZERO(ROAGND, 9*MXFR*MXCOOR)
         CALL DZERO(ROAGLD, 9*MXFR*MXCOOR)
         CALL DZERO(ROAAD ,27*MXFR*MXCOOR)
      ELSE
         X = (IPOINT+1)/2
         IF (MOD(IPOINT,2) .EQ. 1) THEN
            DO 100 I = 1, 3
            DO 100 J = 1, 3
            DO 100 F = 1, NFRVAL
               ROAAFD(I,J,F,X) =(ROAAFP(1,I,J,F) - ROAAFU(I,J,F))/DISPLC
               ROAGND(I,J,F,X) =(ROAGNP(1,I,J,F) - ROAGNU(I,J,F))/DISPLC
               ROAGLD(I,J,F,X) =(ROAGLP(1,I,J,F) - ROAGLU(I,J,F))/DISPLC
            DO 100 K = 1, 3
               ROAAD(I,J,K,F,X)=
     &              (ROAAP(1,I,J,K,F) - ROAAU(I,J,K,F))/DISPLC
  100       CONTINUE
         ELSE
            DO 200 I = 1, 3
            DO 200 J = 1, 3
            DO 200 F = 1, NFRVAL
               ROAAFD(I,J,F,X) = DP5*
     &          (ROAAFD(I,J,F,X)-(ROAAFP(1,I,J,F)-ROAAFU(I,J,F))/DISPLC)
               ROAGND(I,J,F,X) = DP5*
     &          (ROAGND(I,J,F,X)-(ROAGNP(1,I,J,F)-ROAGNU(I,J,F))/DISPLC)
               ROAGLD(I,J,F,X) = DP5*
     &          (ROAGLD(I,J,F,X)-(ROAGLP(1,I,J,F)-ROAGLU(I,J,F))/DISPLC)
            DO 200 K = 1, 3
               ROAAD(I,J,K,F,X) = DP5*
     &       (ROAAD(I,J,K,F,X)-(ROAAP(1,I,J,K,F)-ROAAU(I,J,K,F))/DISPLC)
  200       CONTINUE
         END IF
      END IF
C
C     ***** Print *****
C
      CALL TITLER('ABACUS - Numerical differentiation for VROA','*',114)
      IF (IPOINT .EQ. 0 .AND. IPRINT .GT. 0) THEN
         CALL HEADER('Polarizability tensor at reference geometry',1)
         DO 300 F = 1, NFRVAL
            WRITE (LUPRI,'(1X,A,F6.4)') ' frequency:', FRVAL(F)
            CALL OUTPUT(ROAAFU(1,1,F),1,3,1,3,3,3,1,LUPRI)
  300    CONTINUE
         IF (VROA) THEN
            CALL HEADER('Non-London G tensor at reference geometry',1)
            DO 310 F = 1, NFRVAL
               WRITE (LUPRI,'(1X,A,F6.4)') ' frequency:', FRVAL(F)
               CALL OUTPUT(ROAGNU(1,1,F),1,3,1,3,3,3,1,LUPRI)
  310       CONTINUE
            CALL HEADER('London G tensor at reference geometry',1)
            DO 320 F = 1, NFRVAL
               WRITE (LUPRI,'(1X,A,F6.4)') ' frequency:', FRVAL(F)
               CALL OUTPUT(ROAGLU(1,1,F),1,3,1,3,3,3,1,LUPRI)
  320       CONTINUE
            CALL HEADER('A tensor at reference geometry',1)
            DO 330 F = 1, NFRVAL
               WRITE (LUPRI,'(1X,A,F6.4)') ' frequency:', FRVAL(F)
               CALL OUTPUT(ROAAU(1,1,1,F),1,3,1,9,3,9,1,LUPRI)
  330       CONTINUE
         END IF
C
C     ******* Prints Rayleigh Invariants and scattering parameters
C
         CALL HEADER('Rayleigh properties for parent molecule '//
     &               '(in A**6)',1)
         DO 340 F = 1, NFRVAL
            WRITE (LUPRI,'(2X,A,F8.6,A,F6.2,A,/2X,A)') 
     &         'frequency: ',FRVAL(F),' a.u. = ',XTNM/FRVAL(F),' nm',
     &         '----------------------------------'
            ALMN  = ALFMN  (ROAAFU(1,1,F))
            BAL   = BETAAL (ROAAFU(1,1,F))
            IF (RAMAN) THEN
               RMIN  = RAMINL (ALMN,BAL)
               DPOLR = DEPOLL (ALMN,BAL)
               WRITE (LUPRI,'(/1X,A)') '** Intensities and Depol.'//
     &             ' Ratio for linear polarized incident light **'
               WRITE (LUPRI,'(/5X,A)') 'Alpha**2      '//
     &                'Beta(a)**2    Pol.Int.      Depol.Int.'//
     &                '     Dep.Ratio'
               WRITE (LUPRI,'(5F14.4)') ALMN**2*XFRAIN, BAL*XFRAIN,
     &              RMIN*XFRAIN,RMIN*DPOLR*XFRAIN,DPOLR
            END IF
            IF (VROA) THEN
               RMIN  = RAMINN (ALMN,BAL)
               DPOLR = DEPOLN (ALMN,BAL)
               BGML   = BETAGM (ROAAFU(1,1,F),ROAGLU(1,1,F))
               BGMNL  = BETAGM (ROAAFU(1,1,F),ROAGNU(1,1,F))
               GMMNL  = GMMN   (ROAGLU(1,1,F))
               GMMNNL = GMMN   (ROAGNU(1,1,F))
               BA     = BETAA  (ROAAFU(1,1,F),ROAAU(1,1,1,F),
     &                  FRVAL(F))
               DZL    = DELTAZ (BGML,BA)
               DXL    = DELTAX (BGML,BA,ALMN,GMMNL)
               D0L    = DELTA0 (BGML,BA,ALMN,GMMNL)
               DBL    = DELTAB (BGML,BA)
               CZL    = CID(DZL,RMIN*DPOLR)
               CXL    = CID(DXL,RMIN)
               C0L    = CID(D0L,D2*RMIN)
               CBL    = CID(DBL,D2*RMIN)
               WRITE (LUPRI,'(/1X,A)') '** Intensities and Depol.'//
     &             ' Ratio for circular polarized incident light **'
               WRITE (LUPRI,'(/5X,A)') 'Alpha**2      '//
     &                'Beta(a)**2    Pol.Int.      Depol.Int.'//
     &                '     Dep.Ratio'
               WRITE (LUPRI,'(5F14.4)') ALMN**2*XFRAIN, BAL*XFRAIN,
     &              RMIN*XFRAIN,RMIN*DPOLR*XFRAIN,DPOLR
               WRITE (LUPRI,'(/1X,A)') '** Optical active '//
     &                'Invariants (*E6) **'
               WRITE (LUPRI,'(/7X,A)') ' a*Gm(Lon)     a*Gm(noL)'//
     &                '     Beta(G)**2 (Lon/NoLon)   Beta(A)**2'
               WRITE (LUPRI,'(5F14.4)') ALMN*GMMNL*XFOAIN,
     &             ALMN*GMMNNL*XFOAIN,BGML*XFOAIN,BGMNL*XFOAIN,BA*XFOAIN
               WRITE (LUPRI,'(/1X,A)') '** Difference para'//
     &                'meters L-R (1st line) and '//            
     &                'chirality number (2nd line,*E3) **'          
               WRITE (LUPRI,'(/6X,A,2(/2X,4F14.6),/)') 
     &           'DELTApar(Lon) DELTAperp(Lon) DELTA0(Lon)   '//
     &           'DELTA180(Lon)',
     &           DZL*XFRAIN, DXL*XFRAIN,D0L*XFRAIN,DBL*XFRAIN,
     &           CZL*DE3,CXL*DE3,C0L*DE3,CBL*DE3
            END IF
 340     CONTINUE
C     ***********End Rayleigh (G. Hangartner 24.12.1996)
      ELSE IF (FINAL .OR. (IPRINT.GT.1)) THEN
         ICUR = (IPOINT+1)/2
         WRITE (LUPRI,'(2A,/,A,10X,21X,A,F12.10)')
     &         ' Coordinate displaced in this calculation:    ',
     &         NAMEX(IPTCOR(ICUR,1)),
     &         ' Displacement:', CHRSGN(2*MOD(IPOINT,2)-1), DISPLC
         IF (FINAL) THEN
            WRITE (LUPRI,'(/A)')
     &       ' All displacements are now done'//
     &       ' - numerical differentiation is complete.'
         END IF
C
         CALL HEADER('Polarizability tensor (alpha) gradient',1)
         DO 400 F = 1, NFRVAL
            WRITE (LUPRI,'(1X,A,F6.4)') ' frequency:', FRVAL(F)
            CALL OUTPUT(ROAAFD(1,1,F,1),1,9,1,ICUR,9*MXFR,NCART,1,LUPRI)
  400    CONTINUE
         IF (VROA) THEN
            CALL HEADER('Non-London G tensor gradient',1)
            DO 410 F = 1, NFRVAL
               WRITE (LUPRI,'(1X,A,F6.4)') ' frequency:', FRVAL(F)
               CALL OUTPUT(ROAGND(1,1,F,1),1,9,1,ICUR,9*MXFR,NCART,1,
     &                     LUPRI)
  410       CONTINUE
            CALL HEADER('London G tensor gradient',1)
            DO 420 F = 1, NFRVAL
               WRITE (LUPRI,'(1X,A,F6.4)') ' frequency:', FRVAL(F)
               CALL OUTPUT(ROAGLD(1,1,F,1),1,9,1,ICUR,9*MXFR,NCART,1,
     &                     LUPRI)
  420       CONTINUE
            CALL HEADER('A tensor gradient',1)
            DO 430 F = 1, NFRVAL
                WRITE (LUPRI,'(1X,A,F6.4)') ' frequency:', FRVAL(F)
                CALL OUTPUT(ROAAD(1,1,1,F,1),1,27,1,ICUR,27*MXFR,NCART,
     &                  1,LUPRI)
  430       CONTINUE
         END IF
      END IF
      RETURN
      END
C  /* Deck ircpot */
#if defined (VAR_IRCPOT)
      PROGRAM IRCPOT
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <mxcent.h>
      DIMENSION COOR(MXCOOR), CGRAD(MXCOOR), CHESS(MXCOOR,MXCOOR),
     &          DMASS(MXCENT)
      LUTRJ = -1
      CALL GPOPEN(LUTRJ,'DALTON.IRC','OLD','SEQUENTIAL','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      NPOINT = 2
      REWIND LUTRJ
      DO 100 JPOINT = 1, NPOINT
         READ (LUTRJ)
         READ (LUTRJ) IPOINT, NCOOR
         READ (LUTRJ) AACCUM, ERGMOL
         READ (LUTRJ) (DMASS(I),I=1,NCOOR/3)
         READ (LUTRJ) (COOR (I),I=1,NCOOR)
         READ (LUTRJ) (CGRAD(I),I=1,NCOOR)
         READ (LUTRJ) ((CHESS(I,J),I=1,NCOOR),J=1,NCOOR)
         WRITE (LUPRI,'(//A,I5)') ' Point No. ',IPOINT
         WRITE (LUPRI,'(/A,F12.6)') ' Arc length ', AACCUM
         WRITE (LUPRI,'(/A,F16.8)') ' Energy     ', ERGMOL
         CALL HEADER('Masses',-1)
         CALL OUTPUT(DMASS,1,NCOOR/3,1,1,MXCENT,1,1,LUPRI)
         CALL HEADER('Coordinates',-1)
         CALL OUTPUT(COOR,1,3,1,NCOOR,3,MXCOOR,1,LUPRI)
         CALL HEADER('Gradient',-1)
         CALL OUTPUT(CGRAD,1,3,1,NCOOR,3,MXCOOR,1,LUPRI)
         CALL HEADER('Hessian',-1)
         CALL OUTPUT(CHESS,1,NCOOR,1,NCOOR,MXCOOR,MXCOOR,1,LUPRI)
  100 CONTINUE
      CALL GPCLOSE(LUTRJ,'KEEP')
      STOP
      END
#endif
C  /* Deck wlkort */
      SUBROUTINE WLKORT(CMO,CMOORT,CORNEW,WORK,LWORK,IPRINT)
#include <implicit.h>
#include <maxorb.h>
      DIMENSION CMO(*), CMOORT(*), CORNEW(*), WORK(*)
#include <shells.h>
#include <inforb.h>
C
      KUMO   = 1
      KS11AO = KUMO   + NBAST*NORBT
      KS12AO = KS11AO + N2BASX
      KS22AO = KS12AO + N2BASX
      KS11MO = KS22AO + N2BASX
      KS12MO = KS11MO + N2ORBX
      KS22MO = KS12MO + N2ORBX
      KHALF  = KS22MO + N2ORBX
      KHALFM = KHALF  + NBAST*NORBT
      KS12IN = KHALFM + N2ORBX
      KS22IN = KS12IN + N2ORBX
      KIWRK  = KS22IN + N2ORBX
      KAWRK  = KIWRK  + NORBT
      KWSW   = KAWRK  + NORBT
      KWSWPK = KWSW   + N2ORBX
      KWSWDI = KWSWPK + NNORBX
      KEIGVC = KWSWDI + N2ORBX
      KWSWSQ = KEIGVC + N2ORBX
      KTMAT  = KWSWSQ + N2ORBX
      KOMO   = KTMAT  + N2ORBX
      KSS2OM = KOMO   + NBAST*NORBT
      KCOR1  = KSS2OM + N2ORBX
      KCOR2  = KCOR1  + 3*KMAX
      KLAST  = KCOR2  + 3*KMAX
      LWRK   = LWORK  - KLAST + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('WLKORT',' ',KLAST,LWORK)
      CALL WLKOR1(CMO,CMOORT,CORNEW,WORK(KUMO),WORK(KS11AO),
     &            WORK(KS12AO),WORK(KS22AO),WORK(KS11MO),WORK(KS12MO),
     &            WORK(KS22MO),WORK(KHALF),WORK(KHALFM),WORK(KS12IN),
     &            WORK(KS22IN),WORK(KIWRK),WORK(KAWRK),WORK(KWSW),
     &            WORK(KWSWPK),WORK(KWSWDI),WORK(KEIGVC),WORK(KWSWSQ),
     &            WORK(KTMAT),WORK(KOMO),WORK(KSS2OM),WORK(KCOR1),
     &            WORK(KCOR2),WORK(KLAST),LWRK,IPRINT)
      RETURN
      END
C  /* Deck wlkor1 */
      SUBROUTINE WLKOR1(CMO,CMOORT,CORNEW,UMO,S11AO,S12AO,S22AO,
     &                  S11MO,S12MO,S22MO,HALF,HALFMO,S12INV,S22INV,
     &                  IWRK,AWRK,WSW,WSWPCK,WSWDIA,EIGVEC,WSWSQR,
     &                  TMAT,OMO,S22OMO,COR1,COR2,WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0)
      CHARACTER LABEL(9*MXCENT)*8
      DIMENSION CMO(*), CMOORT(*), CORNEW(3,*),
     &          COR1(KMAX,3), COR2(KMAX,3),
     &          UMO(NBAST,NORBT),
     &          S11AO(NBAST,NBAST), S12AO(NBAST,NBAST),
     &          S22AO(NBAST,NBAST),
     &          S11MO(NORBT,NORBT), S12MO(NORBT,NORBT),
     &          S22MO(NORBT,NORBT),
     &          HALF(NBAST,NORBT), HALFMO(NORBT,NORBT),
     &          S12INV(NORBT,NORBT), S22INV(NORBT,NORBT),
     &          IWRK(NORBT), AWRK(NORBT),
     &          WSW(NORBT,NORBT), WSWPCK(NNORBX),
     &          WSWDIA(NORBT,NORBT), EIGVEC(NORBT,NORBT),
     &          WSWSQR(NORBT,NORBT), TMAT(NORBT,NORBT),
     &          OMO(NBAST,NORBT), S22OMO(NORBT,NORBT),
     &          WORK(LWORK), CMXYZ(3)
#include <cbiwlk.h>
#include <shells.h>
#include <inforb.h>
      INTPRI = MAX(0,IPRINT - 6)
C
C     Get unmodified MO's (UMO's)
C     ===========================
C
      CALL DZERO(UMO,NBAST*NORBT)
      IJ = 0
      DO 100 ISYM = 1, NSYM
      DO 100 J = IORB(ISYM) + 1, IORB(ISYM) + NORB(ISYM)
      DO 100 I = IBAS(ISYM) + 1, IBAS(ISYM) + NBAS(ISYM)
         IJ = IJ + 1
         UMO(I,J) = CMO(IJ)
  100 CONTINUE
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('UMOs in WLKOR1',-1)
         CALL OUTPUT(UMO,1,NBAST,1,NORBT,NBAST,NORBT,1,LUPRI)
      END IF
C
C     Unperturbed and perturbed geometries
C     ====================================
C
      CALL DCOPY(KMAX,CENT(1,1,1),1,COR1(1,1),1)
      CALL DCOPY(KMAX,CENT(1,2,1),1,COR1(1,2),1)
      CALL DCOPY(KMAX,CENT(1,3,1),1,COR1(1,3),1)
      DO 200 I = 1, KMAX
         COR2(I,1) = CORNEW(1,NCENT(I))
         COR2(I,2) = CORNEW(2,NCENT(I))
         COR2(I,3) = CORNEW(3,NCENT(I))
  200 CONTINUE
C
      NCOMP  = 1
      NPATOM = 0
C
C     Get S11AO and transform to MO basis (S12MO)
C     ===========================================
C
      CALL DCOPY(KMAX,COR1(1,1),1,CENT(1,1,1),1)
      CALL DCOPY(KMAX,COR1(1,2),1,CENT(1,2,1),1)
      CALL DCOPY(KMAX,COR1(1,3),1,CENT(1,3,1),1)
      CALL DCOPY(KMAX,COR1(1,1),1,CENT(1,1,2),1)
      CALL DCOPY(KMAX,COR1(1,2),1,CENT(1,2,2),1)
      CALL DCOPY(KMAX,COR1(1,3),1,CENT(1,3,2),1)
      CALL GET1IN(S11AO,'SQOVLAP',NCOMP,WORK,LWORK,LABEL,IREP,DUMMY,
     &            IDUMMY,.FALSE.,NPATOM,.FALSE.,DUMMY,.FALSE.,DUMMY,
     &            INTPRI)
      CALL DCOPY(KMAX,COR1(1,1),1,CENT(1,1,1),1)
      CALL DCOPY(KMAX,COR1(1,2),1,CENT(1,2,1),1)
      CALL DCOPY(KMAX,COR1(1,3),1,CENT(1,3,1),1)
C
      CALL DGEMM('N','N',NBAST,NORBT,NBAST,1.D0,
     &           S11AO,NBAST,
     &           UMO,NBAST,0.D0,
     &           HALF,NBAST)
      CALL DGEMM('T','N',NORBT,NORBT,NBAST,1.D0,
     &           UMO,NBAST,
     &           HALF,NBAST,0.D0,
     &           S11MO,NORBT)
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('S11AO matrix in WLKOR1',-1)
         CALL OUTPUT(S11AO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL HEADER('S11MO matrix in WLKOR1',-1)
         CALL OUTPUT(S11MO,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C     Get S22AO and transform to MO basis (S22MO)
C     ===========================================
C
      CALL DCOPY(KMAX,COR2(1,1),1,CENT(1,1,1),1)
      CALL DCOPY(KMAX,COR2(1,2),1,CENT(1,2,1),1)
      CALL DCOPY(KMAX,COR2(1,3),1,CENT(1,3,1),1)
      CALL DCOPY(KMAX,COR2(1,1),1,CENT(1,1,2),1)
      CALL DCOPY(KMAX,COR2(1,2),1,CENT(1,2,2),1)
      CALL DCOPY(KMAX,COR2(1,3),1,CENT(1,3,2),1)
      CALL GET1IN(S22AO,'SQOVLAP',NCOMP,WORK,LWORK,LABEL,IREP,DUMMY,
     &            IDUMMY,.FALSE.,NPATOM,.FALSE.,DUMMY,.FALSE.,DUMMY,
     &            INTPRI)
      CALL DCOPY(KMAX,COR1(1,1),1,CENT(1,1,1),1)
      CALL DCOPY(KMAX,COR1(1,2),1,CENT(1,2,1),1)
      CALL DCOPY(KMAX,COR1(1,3),1,CENT(1,3,1),1)
C
      CALL DGEMM('N','N',NBAST,NORBT,NBAST,1.D0,
     &           S22AO,NBAST,
     &           UMO,NBAST,0.D0,
     &           HALF,NBAST)
      CALL DGEMM('T','N',NORBT,NORBT,NBAST,1.D0,
     &           UMO,NBAST,
     &           HALF,NBAST,0.D0,
     &           S22MO,NORBT)
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('S22AO matrix in WLKOR1',-1)
         CALL OUTPUT(S22AO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL HEADER('S22MO matrix in WLKOR1',-1)
         CALL OUTPUT(S22MO,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C     Get S12AO and transform to MO basis (S12MO)
C     ===========================================
C
      IF (NATCON) THEN
         CALL DCOPY(KMAX,COR1(1,1),1,CENT(1,1,1),1)
         CALL DCOPY(KMAX,COR1(1,2),1,CENT(1,2,1),1)
         CALL DCOPY(KMAX,COR1(1,3),1,CENT(1,3,1),1)
         CALL DCOPY(KMAX,COR2(1,1),1,CENT(1,1,2),1)
         CALL DCOPY(KMAX,COR2(1,2),1,CENT(1,2,2),1)
         CALL DCOPY(KMAX,COR2(1,3),1,CENT(1,3,2),1)
         CALL GET1IN(S12AO,'SQOVLAP',NCOMP,WORK,LWORK,LABEL,IREP,DUMMY,
     &               IDUMMY,.FALSE.,NPATOM,.FALSE.,DUMMY,.FALSE.,DUMMY,
     &               INTPRI)
         CALL DCOPY(KMAX,COR1(1,1),1,CENT(1,1,1),1)
         CALL DCOPY(KMAX,COR1(1,2),1,CENT(1,2,1),1)
         CALL DCOPY(KMAX,COR1(1,3),1,CENT(1,3,1),1)
C
         CALL DGEMM('N','N',NBAST,NORBT,NBAST,1.D0,
     &              S12AO,NBAST,
     &              UMO,NBAST,0.D0,
     &              HALF,NBAST)
         CALL DGEMM('T','N',NORBT,NORBT,NBAST,1.D0,
     &              UMO,NBAST,
     &              HALF,NBAST,0.D0,
     &              S12MO,NORBT)
C
         IF (IPRINT .GT. 5) THEN
            CALL HEADER('S12AO matrix in WLKOR1',-1)
            CALL OUTPUT(S12AO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            CALL HEADER('S12MO matrix in WLKOR1',-1)
            CALL OUTPUT(S12MO,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
      END IF
C
C     Invert matrices
C     ===============
C
C     S22MO -> S22INV
C
      CALL DGEINV(NORBT,S22MO,S22INV,IWRK,AWRK,INFO)
      IF (INFO .NE. 0) THEN
         WRITE (LUPRI,'(//,A,I5,A,/)')
     &      ' ERROR (WLKTRO) INFO =',INFO,' from DGEINV (S12MO)'
         CALL QUIT('ERROR in DGEINV from WLKTRO')
      END IF
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('S22INV matrix in WLKOR1',-1)
         CALL OUTPUT(S22INV,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C     S12MO -> S12INV
C
      IF (NATCON) THEN
         CALL DGEINV(NORBT,S12MO,S12INV,IWRK,AWRK,INFO)
         IF (INFO .NE. 0) THEN
            WRITE (LUPRI,'(//,A,I5,A,/)')
     &         ' ERROR (WLKTRO) INFO =',INFO,' from DGEINV (S12MO)'
            CALL QUIT('ERROR in DGEINV from WLKTRO')
         END IF
         IF (IPRINT .GT. 5) THEN
            CALL HEADER('S12INV matrix in WLKOR1',-1)
            CALL OUTPUT(S12INV,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
      END IF
C
C     Calculate WSW = S12*S22INV*S21 / S22INV
C     =======================================
C
      IF (NATCON) THEN
         CALL DGEMM('N','N',NORBT,NORBT,NORBT,1.D0,
     &              S12MO,NORBT,
     &              S22INV,NORBT,0.D0,
     &              HALFMO,NORBT)
         CALL DGEMM('N','T',NORBT,NORBT,NORBT,1.D0,
     &              HALFMO,NORBT,
     &              S12MO,NORBT,0.D0,
     &              WSW,NORBT)
         IF (IPRINT .GT. 5) THEN
            CALL HEADER('S12*S22INV*S21 matrix in WLKOR1',-1)
            CALL OUTPUT(WSW,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
         CALL DSITSP(NORBT,WSW,WSWPCK)
      ELSE
         CALL DSITSP(NORBT,S22INV,WSWPCK)
      END IF
C
C     Square root of WSW
C     ==================
C
      CALL DUNIT(EIGVEC,NORBT)
      CALL JACO(WSWPCK,EIGVEC,NORBT,NORBT,NORBT,AWRK,IWRK)
      CALL DZERO(WSWDIA,NORBT*NORBT)
      DO 300 I = 1, NORBT
         WSWDIA(I,I) = SQRT(WSWPCK(I*(I+1)/2))
  300 CONTINUE
      CALL DGEMM('N','N',NORBT,NORBT,NORBT,1.D0,
     &           EIGVEC,NORBT,
     &           WSWDIA,NORBT,0.D0,
     &           HALFMO,NORBT)
      CALL DGEMM('N','T',NORBT,NORBT,NORBT,1.D0,
     &           HALFMO,NORBT,
     &           EIGVEC,NORBT,0.D0,
     &           WSWSQR,NORBT)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('WSWSQR matrix in WLKOR1',-1)
         CALL OUTPUT(WSWSQR,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C     The T matrix
C     ============
C
      IF (NATCON) THEN
         CALL DGEMM('N','N',NORBT,NORBT,NORBT,1.D0,
     &              S12INV,NORBT,
     &              WSWSQR,NORBT,0.D0,
     &              TMAT,NORBT)
      ELSE
         CALL DCOPY(N2ORBX,WSWSQR,1,TMAT,1)
      END IF
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('T matrix in WLKOR1',-1)
         CALL OUTPUT(TMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C     OMOs
C     ====
C
      CALL DGEMM('N','N',NBAST,NORBT,NORBT,1.D0,
     &           UMO,NBAST,
     &           TMAT,NORBT,0.D0,
     &           OMO,NBAST)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('OMOs in WLKOR1',-1)
         CALL OUTPUT(OMO,1,NBAST,1,NORBT,NBAST,NORBT,1,LUPRI)
      END IF
C
C     Test of orthogonality
C     =====================
C
      CALL DGEMM('N','N',NBAST,NORBT,NBAST,1.D0,
     &           S22AO,NBAST,
     &           OMO,NBAST,0.D0,
     &           HALF,NBAST)
      CALL DGEMM('T','N',NORBT,NORBT,NBAST,1.D0,
     &           OMO,NBAST,
     &           HALF,NBAST,0.D0,
     &           S22OMO,NORBT)
C
      DIFMAX = D0
      DO 600 I = 1, NORBT
      DO 600 J = 1, NORBT
         DIFFER = S22OMO(I,J)
         IF (I .EQ. J) DIFFER = DIFFER - D1
         ABSDIF = ABS(DIFFER)
         IF (ABSDIF .GE. DIFMAX) THEN
            DIFMAX = ABSDIF
            IDIFMX = I
            JDIFMX = J
         END IF
  600 CONTINUE
      IF (IPRINT .GT. 5 .OR. DIFMAX .GT. 1.D-7) THEN
         CALL HEADER('OMO overlap matrix in WLKOR1',-1)
         CALL OUTPUT(S22OMO,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
      WRITE (LUPRI,'(/A,E12.6,/A,2I4)')
     &   ' Largest deviation from orthonormality of OMOs:',DIFMAX,
     &   ' - found for element ',IDIFMX,JDIFMX

C
C     Test of similarity
C     ==================
C
      IF (NATCON) THEN
         CALL DGEMM('N','N',NBAST,NORBT,NBAST,1.D0,
     &              S12AO,NBAST,
     &              OMO,NBAST,0.D0,
     &              HALF,NBAST)
         CALL DGEMM('T','N',NORBT,NORBT,NBAST,1.D0,
     &              UMO,NBAST,
     &              HALF,NBAST,0.D0,
     &              WSW,NORBT)
         DO 500 I = 1, NORBT
         DO 500 J = 1, NORBT
            S22OMO(I,J) = (WSW(I,J) - WSW(J,I))/D2
  500    CONTINUE
         DIFMAX = D0
         DO 700 I = 1, NORBT
         DO 700 J = 1, NORBT
            DIFFER = S22OMO(I,J) - S22OMO(J,I)
            ABSDIF = ABS(DIFFER)
            IF (ABSDIF .GE. DIFMAX) THEN
               DIFMAX = ABSDIF
               IDIFMX = I
               JDIFMX = J
            END IF
  700    CONTINUE
         IF (IPRINT .GT. 5 .OR. DIFMAX .GT. 1.D-7) THEN
            CALL HEADER('WT matrix in WLKOR1',-1)
            CALL OUTPUT(WSW,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            CALL HEADER('Antisymmetric part of WT in WLKOR1',-1)
            CALL OUTPUT(S22OMO,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
         WRITE (LUPRI,'(/A,E12.6,/A,2I4)')
     &      ' Largest deviation from Hermiticity of UMO/OMO overlaps:',
     &        DIFMAX,
     &      ' - found for element ',IDIFMX,JDIFMX
      END IF
C
C     Pack OMO coefficients
C     =====================
C
      IJ = 0
      DO 400 ISYM = 1, NSYM
      DO 400 J = IORB(ISYM) + 1, IORB(ISYM) + NORB(ISYM)
      DO 400 I = IBAS(ISYM) + 1, IBAS(ISYM) + NBAS(ISYM)
         IJ = IJ + 1
         CMOORT(IJ) = OMO(I,J)
  400 CONTINUE
      RETURN
      END
C  /* Deck wlkcvt */
      SUBROUTINE WLKCVT(TMAT,IPTCOL,NCR,NPR,NTR,IREP,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
      PARAMETER (D0 = 0.00D0, D1 = 1.00D0)
      DIMENSION TMAT(NCR,NPR), IPTCOL(*)
#include <cbisol.h>
#include <nuclei.h>
#include <symmet.h>
      AMASS(I) = DISOTP(IZATOM(I),ISOTOP(I),'MASS')
C
C     Total mass
C
      TMASS = D0
      DO 100 IATOM = 1, NUCIND
         IF (IATOM .NE. NCNTCV) THEN
            TMASS = TMASS + FMULT(ISTBNU(IATOM))*AMASS(IATOM)
         END IF
  100 CONTINUE
C
      DO 200 IPRJ  = 1, NAXREP(IREP,1)
      DO 200 IATOM = 1, NUCIND
         ISCOOR = IPTCNT(3*(IATOM - 1) + IPTXYZ(IPRJ,IREP,1),IREP,1)
         IF (IATOM .EQ. NCNTCV) THEN
            TMT = D1
         ELSE
            TMT = - FMULT(ISTBNU(IATOM))*AMASS(IATOM)/TMASS
         END IF
         TMAT(IPTCOL(ISCOOR),NTR + IPRJ) = TMT
  200 CONTINUE
C
      RETURN
      END
C /* Deck tredif */
      SUBROUTINE TREDIF(PROP0,PROP,PROPD,DISPLC,ODD,IDIM,IABS)
C
C     Three-point numerical differentiation
C     IABS tells whether the imaginary component of the response function
C     also was calculated. Currently we do not average this component.
C     K.Ruud, June 2001
C
#include <implicit.h>
      LOGICAL ODD
      DIMENSION PROP0(IDIM), PROP(IABS,IDIM), PROPD(IDIM)
C
      IF (ODD) THEN
         DO I = 1, IDIM
            PROPD(I) = (PROP(IABS,I)-PROP0(I))/DISPLC**2
         END DO
      ELSE
         DO I = 1, IDIM
            PROPD(I) = PROPD(I) + (PROP(IABS,I)-PROP0(I))/DISPLC**2
         END DO
      END IF
      RETURN
      END
C /* Deck femdif */
      SUBROUTINE FEMDIF(PROP0,PROP,PROPD,DISPLC,IPOINT,IDIM,IABS)
C
C     Five-point numerical differentiation
C     IABS tells whether the imaginary component of the response function
C     also was calculated. Currently we do not average this component.
C     K.Ruud, June 2001
C
#include <implicit.h>
      PARAMETER (D3 = 3.0D0, D4 = 4.0D0, D12 = D3*D4)
      DIMENSION PROP0(IDIM), PROP(IABS,IDIM), PROPD(IDIM)
C
      IF (IPOINT .EQ. 1) THEN
         DO I = 1, IDIM
            PROPD(I) = (D4/D3)*(PROP(IABS,I)-PROP0(I))/(DISPLC**2)
         END DO
      ELSE IF (IPOINT .EQ. 2 .OR. IPOINT .EQ. 3) THEN
         DO I = 1, IDIM
            PROPD(I) = PROPD(I)
     &               - (PROP(IABS,I)-PROP0(I))/(D12*(DISPLC**2))
         END DO
      ELSE
         DO I = 1, IDIM
            PROPD(I) = PROPD(I)
     &               + (D4/D3)*(PROP(IABS,I)-PROP0(I))/(DISPLC**2)
         END DO
      END IF
      RETURN
      END
C --- end of abawalk.F ---
