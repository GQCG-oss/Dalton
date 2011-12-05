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
C FILE: abaopt.F
#include "single.h"
C
C 950915-vebjornb: New module taking care of optimizations.
C                  It uses both the old cbiwlk.h and the new optinf.h
C
C  /* Deck opinpu */
      SUBROUTINE OPINPU(WORD)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "optinf.h"
#include "infinp.h"
#include "gnrinf.h"
#include "cbiwlk.h"
      PARAMETER (NTABLE = 75, ITRMX = 25, MXREJ = 3, D0 = 0.0D0,
     &     DEFTHE = 1.0D-6, DEFTH2 = 1.0D-4, DEFDSP = 1.0D-3)
      LOGICAL NEWDEF, TRSTCH, TRSTFC, FIRST, REMCRD
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
      SAVE FIRST
      DATA TABLE /'.PRINT ', '.MAX IT', '.TRUSTR', '.TR FAC', !  4
     &            '.TR LIM', '.MAX RE', '.NOTRUS', '.ENERGY', !  8
     &            '.GRADIE', '.STEP T', '.CONDIT', '.NOBREA', ! 12
     &            '.SP BAS', '.PREOPT', '.VISUAL', '.VRML  ', ! 16
     &            '.SYMTHR', '.TRSTRG', '.VR-BON', '.VR-EIG', ! 20
     &            '.INITHE', '.INITEV', '.HESFIL', '.REJINI', ! 24
     &            '.STEEPD', '.RANKON', '.PSB   ', '.DFP   ', ! 28
     &            '.BFGS  ', '.NEWTON', '.QUADSD', '.SCHLEG', ! 32
     &            '.HELLMA', '.BAKER ', '.M-BFGS', '.CARTES', ! 36
     &            '.REDINT', '.INIRED', '.1STORD', '.2NDORD', ! 40
     &            '.GRDINI', '.DISPLA', '.CONSTR', '.MODHES', ! 44
     &            '.REMOVE', '.INIMOD', '.FINDRE', '.CMBMOD', ! 48
     &            '.RF    ', '.GDIIS ', '.DELINT', '.NODIHE', ! 52
     &            '.VR-COR', '.VR-VIB', '.VR-SYM', '.M-PSB ', ! 56
     &            '.LINE S', '.SADDLE', '.MODE  ', '.BOFILL', ! 60
     &            '.NOAUX ', '.BFGSR1', '.STABIL', '.GEOANA', ! 64
     &            '.ADDCRD', '.NOADDA', '.RREDUN', '.NATNRM', ! 68
     &            '.VLOOSE', '.LOOSE ', '.TIGHT ', '.VTIGHT', ! 72
     &            '.NOHSWR', '.FREEZE', '.FRZITR'/
C
      DATA FIRST /.TRUE./
C
      IF (.NOT. FIRST) THEN
         IF ((WORD .NE. '*END OF') .AND. (WORD(1:2) .NE. '**')) THEN
 969        READ(LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .NE. '*') GOTO 969
            CALL UPCASE(WORD)
         END IF
         RETURN
      END IF
      FIRST = .FALSE.
C
      NEWDEF = WORD .EQ. '*OPTIMI'
      ICHANG = 0
      ICHGTH = 0
      ICNLVL = -1
      IPDEF  = 0
      NSPMOD = -1
      SADDLE = .FALSE.
      NOTRST = .FALSE.
      NOBRKS = .FALSE.
      BRKSYM = .FALSE.
      ITRBRK = -1
      NWSYMM = .FALSE.
      NEWSYM = .FALSE.
      NOHSWR = .FALSE.
      DOSPE  = .FALSE.
      DOPRE  = .FALSE.
      FINPRE = .FALSE.
      KEEPHE = .FALSE.
      REJINI = .FALSE.
      GRDINI = .FALSE.
      CHGRDT = .FALSE.
      CONOPT = .FALSE.
      REMCRD = .FALSE.
      ADDCRD = .FALSE.
      REBILD = .FALSE.
      NUMPRE = 0
      IPRE   = 0
      ITOTRJ = 0
      IREDIC = -1
      IINTCR = 0
      ICRTCR = 0
      ISTBLZ = -1
      STEEPD = .FALSE.
      RANKON = .FALSE.
      PSB    = .FALSE.
      DFP    = .FALSE.
      BFGS   = .FALSE.
      BOFILL = .FALSE.
      BFGSR1 = .FALSE.
      MULTI  = .FALSE.
      SCHLEG = .FALSE.
      NEWTON = .FALSE.
      QUADSD = .FALSE.
      FSTORD = .FALSE.
      SNDORD = .FALSE.
      DELINT = .FALSE.
      REDINT = .FALSE.
      CARTCO = .FALSE.
      INRDHS = .FALSE.
      INITHS = .FALSE.
      MODHES = .FALSE.
      CMBMOD = .FALSE.
      INMDHS = .FALSE.
      HSFILE = .FALSE.
      EVLINI = -1.0D0
      FINDRE = .FALSE.
      VISUAL = .FALSE.
      VRML   = .FALSE.
      VRBOND = .FALSE.
      VREIGV = .FALSE.
      VRCORD = .FALSE.
      VRVIBA = .FALSE.
      VRML_SYM = .FALSE.
      TRSTCH = .FALSE.
      TRSTFC = .FALSE.
      IPRINT = 0
      TRSTRA = 0.5D0
      TRSTIN = 1.2D0
      TRSTDE = 0.7D0
      RTENBD = 0.4D0
      RTENGD = 0.8D0
      RTRJMN = 0.1D0
      RTRJMX = 3.0D0
      LNSRCH = .FALSE.
      RATFUN = .FALSE.
      TRSTRG = .FALSE.
      GDIIS  = .FALSE.
      GECONV = .FALSE.
      GETOLL = DEFTHE
      BAKER  = .FALSE.
      NOAUX  = .FALSE.
      NODIHE = .FALSE.
      NOADDA = .FALSE.
      REDRED = .FALSE.
      LINDHD = .TRUE.
C     ... hjaaj Mar 2004: disable generalization of Lindh
C     algorithm for model Hessian as this is not working
C     satisfactorily yet.
      ITRNMR = 0
      ITRMAX = ITRMX
      MAXREJ = MXREJ
      DISPLA = DEFDSP
      GRDTHR = DEFTH2
      THRSTP = DEFTH2
      THRERG = DEFTHE
C     When running linsca version, NATNRM is always used
#ifdef VAR_LINSCA
      NATNRM = .TRUE.
#else
      NATNRM = .FALSE.
#endif
      THRSYM = MIN(5.0D-3,SQRT(DEFTHE))
      ICONDI = 2
      NCRTOT = 0
      NCART  = 0
      ITRFRZ = -1
      ENERGY = D0
      ERGOLD = D0
      ERGPRD = D0
      ERGPRO = D0
      STPNRM = D0
      STPNRO = D0
      GRADNM = D0
      ZERGRD = 1.0D-7
      CALL DZERO(STPDIA,8*MXCENT)
      CALL DZERO(STPSYM,8*MXCENT)
      CALL DZERO(STPINT,8*MXCENT)
      CALL DZERO(GRDDIA,8*MXCENT)
      CALL DZERO(EVAL  ,8*MXCENT)
      CALL DZERO(EVALOL,8*MXCENT)
      CALL DZERO(CRDINT,8*MXCENT)
      CALL DZERO(CRDIN1,8*MXCENT)
      CALL IZERO(ICNSTR,8*MXCENT)
      CALL IZERO(IADDCR,11*4)
      CALL IZERO(IFREEZ,11)
      DO 123 I = 0, 7
         CNDHES(I) = D0
         INDHES(I) = 0
 123  CONTINUE
C
      WORD1 = '*OPTIMI'
      IF (NEWDEF) THEN
         WORD1 = WORD
 100     CONTINUE
         READ (LUCMD, '(A7)') WORD
         CALL UPCASE(WORD)
         PROMPT = WORD(1:1)
         IF ((PROMPT .EQ. '!') .OR. (PROMPT .EQ. '#')) THEN
            GOTO 100
         ELSE IF (PROMPT .EQ. '.') THEN
            ICHANG = ICHANG + 1
            DO 200 I = 1, NTABLE
               IF (TABLE(I) .EQ. WORD) THEN
                  GOTO ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &                  11,12,13,14,15,16,17,18,19,20,
     &                  21,22,23,24,25,26,27,28,29,30,
     &                  31,32,33,34,35,36,37,38,39,40,
     &                  41,42,43,44,45,46,47,48,49,50,
     &                  51,52,53,54,55,56,57,58,59,60,
     &                  61,62,63,64,65,66,67,68,69,70,
     &                  71,72,73,74,75), I
               END IF
 200        CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               GOTO 100
            END IF
            WRITE(LUPRI,'(/4A/)') ' Keyword "',WORD,
     &                         '" not recognized for ',WORD1
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal keyword in '//WORD1//' input.')
 1          READ(LUCMD,*) IPRINT
            GOTO 100
 2          READ(LUCMD,*) ITRMAX
            GOTO 100
 3          READ(LUCMD,*) TRSTRA
            TRSTCH = .TRUE.
            GOTO 100
 4          READ(LUCMD,*) TRSTIN, TRSTDE
            TRSTCH = .TRUE.
            TRSTFC = .TRUE.
            GOTO 100
 5          READ(LUCMD,*) RTENBD, RTENGD, RTRJMN, RTRJMX
            TRSTCH = .TRUE.
            GOTO 100
 6          READ(LUCMD,*) MAXREJ
            GOTO 100
 7          NOTRST = .TRUE.
            GOTO 100
 8          READ(LUCMD,*) THRERG
            CHGRDT = .TRUE.
            GOTO 100
 9          READ(LUCMD,*) GRDTHR
            CHGRDT = .TRUE.
            ICHGTH = ICHGTH + 1
            GOTO 100
 10         READ(LUCMD,*) THRSTP
            CHGRDT = .TRUE.
            ICHGTH = ICHGTH + 1
            GOTO 100
 11         READ(LUCMD,*) ICONDI
            GOTO 100
 12         NOBRKS = .TRUE.
            GOTO 100
 13         READ(LUCMD,*) SPBSTX
            DOSPE = .TRUE.
            GOTO 100
 14         READ(LUCMD,*) NUMPRE
            IF ((NUMPRE .LT. 1) .OR. (NUMPRE .GT. MAXPRE)) THEN
               WRITE(LUPRI,'(/,A,I2,A/)')
     &        ' Number of preoptimization sets must be between 1 and',
     &            MAXPRE,'!'
               CALL QUIT('Illegal number of preoptimization sets.')
            ELSE
               DO 144 I = 1, NUMPRE 
                  READ(LUCMD,'(A60)') PREBTX(I)
 144           CONTINUE
               DOPRE = .TRUE.
            END IF
            GOTO 100
 15         VISUAL = .TRUE.
            GOTO 100
 16         VRML = .TRUE.
            GOTO 100
 17         READ(LUCMD,*) THRSYM
            GOTO 100
 18         TRSTRG = .TRUE.
            GOTO 100
 19         VRBOND = .TRUE.
            GOTO 100
 20         VREIGV = .TRUE.
            GOTO 100
 21         INITHS = .TRUE.
            GOTO 100
 22         READ(LUCMD,*) EVLINI
            GOTO 100
 23         HSFILE = .TRUE.
            GOTO 100
 24         REJINI = .TRUE.
            GOTO 100
 25         STEEPD = .TRUE.
            GOTO 100
 26         RANKON = .TRUE.
            GOTO 100
 27         PSB = .TRUE.
            GOTO 100
 28         DFP    = .TRUE.
            GOTO 100
 29         BFGS   = .TRUE.
            GOTO 100
 30         NEWTON = .TRUE.
            GOTO 100
 31         QUADSD = .TRUE.
            GOTO 100
 32         SCHLEG = .TRUE.
            GOTO 100
 33         HFPROP = .TRUE.
            GOTO 100
 34         BAKER  = .TRUE.
            GOTO 100
 35         MULTI  = .TRUE.
            BFGS   = .TRUE.
            GOTO 100
 36         CARTCO = .TRUE.
            GOTO 100
 37         REDINT = .TRUE.
            GOTO 100
 38         INRDHS = .TRUE.
            GOTO 100
 39         FSTORD = .TRUE.
            GOTO 100
 40         SNDORD = .TRUE.
            GOTO 100
 41         GRDINI = .TRUE.
            GOTO 100
 42         READ (LUCMD,*) DISPLA
            GOTO 100
 43         READ (LUCMD,*) NCON
            IF ((NCON .LE. 0) .OR. (NCON .GE. 8*MXCENT)) THEN
               WRITE(LUPRI,'(/,A,I2,A/)')
     &              ' Number of constrained coordinates must '
     &              // 'be between 1 and', (8*MXCENT-1),'!'
               CALL QUIT('Illegal number of constrained coordinates.')
            END IF
            DO 443 I = 1, NCON
               READ (LUCMD,*) ICON
               IF ((ICON .LE. 0) .OR. (ICON .GT. 8*MXCENT)) THEN
                  WRITE(LUPRI,'(/,A,I2,A/)')
     &                 ' The constrained coordinate should have a '
     &                 // 'value between 1 and', (8*MXCENT),'!'
                  CALL QUIT('Illegal constrained coordinate.')
               END IF
               ICNSTR(ICON) = 1
 443        CONTINUE
            CONOPT = .TRUE.
            GOTO 100
 44         MODHES = .TRUE.
            GOTO 100
 45         READ (LUCMD,*) NREM
            IF ((NREM .LE. 0) .OR. (NREM .GE. 8*MXCENT)) THEN
               WRITE(LUPRI,'(/,A,I2,A/)')
     &              ' Number of coordinates to be removed must '
     &              // 'be between 1 and', (8*MXCENT-1),'!'
               CALL QUIT('Illegal number coordinates to be removed.')
            END IF
            DO 445 I = 1, NREM
               READ (LUCMD,*) IREM
               IF ((IREM .LE. 0) .OR. (IREM .GT. 8*MXCENT)) THEN
                  WRITE(LUPRI,'(/,A,I2,A/)')
     &                 ' The coordinate number should be a '
     &                 // 'value between 1 and', (8*MXCENT),'!'
                  CALL QUIT('Illegal coordinate to be removed.')
               END IF
               ICNSTR(IREM) = 2
 445        CONTINUE
            REMCRD = .TRUE.
            GOTO 100
 46         INMDHS = .TRUE.
            GOTO 100
 47         FINDRE = .TRUE.
            GOTO 100
 48         CMBMOD = .TRUE.
            GOTO 100
 49         RATFUN = .TRUE.
            GOTO 100
 50         GDIIS = .TRUE.
            GOTO 100
 51         DELINT = .TRUE.
            GOTO 100
 52         NODIHE = .TRUE.
            GOTO 100
 53         VRCORD = .TRUE.
            GOTO 100
 54         VRVIBA = .TRUE.
            GOTO 100
C           .VR-SYM
 55         VRML_SYM = .TRUE.
            GOTO 100
 56         MULTI = .TRUE.
            PSB = .TRUE.
            GOTO 100
 57         LNSRCH = .TRUE.
            GOTO 100
 58         SADDLE = .TRUE.
            GOTO 100
 59         READ (LUCMD,*) NSPMOD
            GOTO 100
 60         BOFILL = .TRUE.
            GOTO 100
 61         NOAUX = .TRUE.
            GOTO 100
 62         BFGSR1 = .TRUE.
            GOTO 100
 63         READ (LUCMD,*) ISTBLZ
            GOTO 100
 64         GEOALL = .TRUE.
            GOTO 100
 65         READ (LUCMD,*) NADD
            IF ((NADD .LE. 0) .OR. (NADD .GT. 10)) THEN
               WRITE(LUPRI,'(/,A/)')
     &              ' Number of coordinates to be added must '
     &              // 'be between 1 and 10!'
               CALL QUIT('Illegal number of coordinates to be added.')
            END IF
            IADDCR(0,1) = NADD
            DO 465 I = 1, NADD
               READ (LUCMD,*) IADDCR(I,1), IADDCR(I,2),
     &              IADDCR(I,3), IADDCR(I,4)
 465        CONTINUE
            ADDCRD = .TRUE.
            GOTO 100
 66         NOADDA = .TRUE.
            GOTO 100
 67         REDRED = .TRUE.
            GOTO 100
 68         NATNRM = .TRUE.
            GOTO 100
 69         ICNLVL = 1
            ICHGTH = ICHGTH + 1
            GOTO 100
 70         ICNLVL = 2
            ICHGTH = ICHGTH + 1
            GOTO 100
 71         ICNLVL = 4
            ICHGTH = ICHGTH + 1
            GOTO 100
 72         ICNLVL = 5
            ICHGTH = ICHGTH + 1
            GOTO 100
 73         NOHSWR = .TRUE.
            GOTO 100
 74         READ (LUCMD,*) NFREEZ
            IF ((NFREEZ .LE. 0) .OR. (NFREEZ .GE. MXCENT)) THEN
               WRITE(LUPRI,'(/A,I6,A/)')
     &              ' Number of frozen atoms must '
     &              // 'be between 1 and', MXCENT,'!'
               CALL QUIT('Illegal number of frozen atoms.')
            END IF
            DO 474 I = 1, NFREEZ
               READ (LUCMD,*) NFRZ
               IF ((NFRZ .LE. 0) .OR. (NFRZ .GT. MXCENT)) THEN
                  WRITE(LUPRI,'(/A,I6,A/)')
     &                 ' The constrained coordinate should have a '
     &                 // 'value between 1 and', MXCENT,'!'
                  CALL QUIT('Illegal frozen atom.')
               END IF
               IFREEZ(0) = IFREEZ(0) + 1
               IFREEZ(IFREEZ(0)) = NFRZ
 474        CONTINUE
            GOTO 100
 75         READ (LUCMD,*) ITRFRZ
            GOTO 100
         ELSE IF (PROMPT .EQ. '*') THEN
            GOTO 300
         ELSE
            CALL QUIT('Illegal prompt in OPINPU.')
            WRITE(LUPRI,'(/4A/)')
     &         ' Prompt "',WORD,'" not recognized for ',WORD1
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal prompt in '//WORD1//' input.')
         END IF
      END IF
 300  CONTINUE
C
      IF ((.NOT. OPTNEW) .AND. (ICHANG .GT. 0)) THEN
         WRITE(LUPRI,'(//A/3A/A/)') SEPARATOR,
     &     ' Geometry optimization not requested, ',
     &     WORD1,' input ignored.', SEPARATOR
         RETURN
      END IF
C
C     If a saddle point optimization has been requested,
C     the trust radius modifiers are adjusted.
C
C      IF (SADDLE .AND. (.NOT. TRSTFC)) THEN
C         TRSTIN = 1.5D0
C         TRSTDE = 0.5D0
C      END IF
C
      IF (.NOT. OPTNEW) RETURN
      CALL HEADER('Chosen parameters for '//WORD1//' :',0)
C
C     Check if only redundant internal coordinates should be determinded.
C
      IF (FINDRE) THEN
         WRITE(LUPRI,'(A)') ' Determination of redundant internal '
     &        // 'will be performed.'
         IF (.NOT. VISUAL) THEN
            WRITE(LUPRI,'(/A/A/)') ' *** NOTE! ***',
     &         ' No geometry optimization will be done, ' //
     &         'other keywords will be ignored!!!!!'
         END IF
      END IF
C
C     Check for visualization
C
      IF (VISUAL) THEN
         WRITE(LUPRI,'(A)') ' Visualization has been ' //
     &        'requested. No geometry optimization will be done.',
     &        ' VRML-file of geometry will be created.'
         IF (VRBOND) WRITE(LUPRI,'(A)')
     &        ' Bonds will be drawn between nearby atoms.'
         IF (VRCORD) WRITE(LUPRI,'(A)')
     &        ' Coordinate axes will be drawn.'
         IF (VRML_SYM) THEN
            WRITE(LUPRI,'(A/A)')
     &        ' Symmetry elements will be visualized.',
     &        ' Please note that symmetry should NOT be explicitly' //
     &        ' specified in the input file.'
         END IF
         IF (VREIGV) THEN
            WRITE(LUPRI,'(A)') ' Eigenvectors can ' //
     &       'only be visualized during an optimization.',
     &       ' Keyword will be ignored!'
            VREIGV = .FALSE.
         END IF
         IF (VRVIBA) THEN
            WRITE(LUPRI,'(A)') ' Vibrational modes can ' //
     &       'only be visualized during an optimization.',
     &       ' Keyword will be ignored!'
            VRVIBA = .FALSE.
         END IF
         WRITE(LUPRI,'(A)')
     &      ' Any other keywords in this module are ignored!'
         RETURN
      END IF
C
C     The type of optimization is determined, BFGS is default.
C
      IF (SADDLE) THEN
         WRITE(LUPRI,'(A)')
     &        ' Saddle point optimization has been requested.'
         NOAUX = .TRUE.
      END IF
      ITYP = 0
      IF (STEEPD) THEN
         WRITE(LUPRI,'(A)')
     &        ' 1st order steepest descent method will be used.'
         ITYP = ITYP + 1
      END IF
      IF (RANKON) THEN
         IF (MULTI) THEN
            WRITE(LUPRI,'(A)') ' 1st order method with ' //
     &           '"multiple rank one" updating scheme will be used.'
         ELSE
            WRITE(LUPRI,'(A)')
     &           ' 1st order method with rank one update will be used.'
         END IF
         ITYP = ITYP + 1
      END IF
      IF (BOFILL) THEN
         MULTI = .FALSE.
         WRITE(LUPRI,'(A)')
     &        ' 1st order method with Bofills update will be used.'
         ITYP = ITYP + 1
      END IF
      IF (PSB) THEN
         IF (MULTI) THEN
            WRITE(LUPRI,'(A)') ' 1st order method with ' //
     &           '"multiple PSB" updating scheme will be used.'
         ELSE
            WRITE(LUPRI,'(A)')
     &           ' 1st order method with PSB update will be used.'
         END IF
         ITYP = ITYP + 1
      END IF
      IF (DFP) THEN
         IF (MULTI) THEN
            WRITE(LUPRI,'(A)') ' 1st order method with ' //
     &           '"multiple DFP" updating scheme will be used.'
         ELSE
            WRITE(LUPRI,'(A)')
     &           ' 1st order method with DFP update will be used.'
         END IF
         ITYP = ITYP + 1
      END IF
      IF (BFGS) THEN
         IF (MULTI) THEN
            WRITE(LUPRI,'(A)') ' 1st order method with ' //
     &           '"multiple BFGS" updating scheme will be used.'
         ELSE
            WRITE(LUPRI,'(A)')
     &           ' 1st order method with BFGS update will be used.'
         END IF
         ITYP = ITYP + 1
      END IF
      IF (BFGSR1) THEN
         WRITE(LUPRI,'(A)')
     &           ' 1st order method with BFGS/rank one combination ' // 
     &           'update will be used.'
         ITYP = ITYP + 1
      END IF
      IF (SCHLEG) THEN
         WRITE(LUPRI,'(A)')
     &  ' 1st order method with Schlegels updating scheme will be used.'
         ITYP = ITYP + 1
      END IF
      IF (NEWTON) THEN
         WRITE(LUPRI,'(A)')
     &        ' 2nd order Newton method will be used.'
         ITYP = ITYP + 1
      END IF
      IF (QUADSD) THEN
         WRITE(LUPRI,'(A)')
     &        ' 2nd order quadratic steepest ' //
     &        'descent method will be used.'
         ITYP = ITYP + 1
      END IF
      IF (FSTORD) THEN
         IF (SADDLE) THEN
            WRITE(LUPRI,'(A)')
     &           ' Default 1st order TS-method will be used:' //
     &           '   Bofills update.'
            BOFILL = .TRUE.
         ELSE
            WRITE(LUPRI,'(A)')
     &           ' Default 1st order method will be used:' //
     &           '   BFGS update.'
            BFGS = .TRUE.
         END IF
         REJINI = .TRUE.
         ITYP = ITYP + 1
      END IF
      IF (SNDORD) THEN
         WRITE(LUPRI,'(A)')
     &        ' Default 2nd order method will be used:   Newton method.'
         NEWTON = .TRUE.
         ITYP = ITYP + 1
      END IF
      IF (ITYP .EQ. 0) THEN
         IF (SADDLE) THEN
            WRITE(LUPRI,'(A)')
     &           ' Default 1st order TS-method will be used:' //
     &           '   Bofills update.'
            BOFILL = .TRUE.
         ELSE
            WRITE(LUPRI,'(A)')
     &           ' Default 1st order method will be used:' //
     &           '   BFGS update.'
            BFGS = .TRUE.
         END IF
         FSTORD = .TRUE.
      ELSE IF (ITYP .GT. 1) THEN
         WRITE(LUPRI,'(/A)') ' WARNING! More than one ' //
     &        'optimization method has been selected under '//WORD1
         CALL QUIT
     &        ('More than one optimization method chosen for '//WORD1)
      END IF
      IF (HFPROP) THEN
         WRITE(LUPRI,'(A)') ' The Hellmann-Feynman theorem will ' //
     &        'be utilized to calculate derivatives.'
         WRITE (LUPRI,'(5X,A)') 'This option is currently not working'//
     &        ' correctly, program will stop'
         CALL QUIT('Hellmann-Feynman approximation not working')
      END IF
C
      ITYP = 0
      IF (CARTCO) THEN
         WRITE(LUPRI,'(A)') ' Optimization will be performed ' //
     &        'in Cartesian coordinates.'
         ITYP = ITYP + 1
      END IF
      IF (REDINT) THEN
         WRITE(LUPRI,'(A)') ' Optimization will be performed ' //
     &        'in redundant internal coordinates.'
         ITYP = ITYP + 1
      END IF
      IF (DELINT) THEN
         WRITE(LUPRI,'(A)') ' Optimization will be performed ' //
     &        'in delocalized internal coordinates.'
         ITYP = ITYP + 1
      END IF
      IF (ITYP .EQ. 0) THEN
         IF (FSTORD) THEN
            REDINT = .TRUE.
            WRITE(LUPRI,'(A)') ' Optimization will be performed ' //
     &           'in redundant internal coordinates.'
         ELSE
            CARTCO = .TRUE.
            WRITE(LUPRI,'(A)') ' Optimization will be performed ' //
     &           'in Cartesian coordinates.'
         END IF
      ELSE IF (ITYP .GT. 1) THEN
         WRITE(LUPRI,'(/A)') ' WARNING! More than one ' //
     &        'coordinate system has been selected in '//WORD1
         CALL QUIT
     &        ('More than one coordinate system chosen under '//WORD1)
      END IF
      IF (NOAUX) THEN
            WRITE(LUPRI,'(A)') ' No extra (auxiliary) bonds will ' //
     &        'be added.'
      END IF
      IF (NODIHE) THEN
            WRITE(LUPRI,'(A)') ' No dihedral angles will be used ' //
     &           'as coordinates (just bonds and angles).'
      END IF
      IF (NOADDA) THEN
            WRITE(LUPRI,'(A)') ' No stabilizing additional angles ' //
     &           'will be used.'
      END IF
      IF (REDRED) THEN
            WRITE(LUPRI,'(A)') ' Redundancy will be decreased by ' //
     &           'removing a number of dihedral coordinates.'
      END IF
C
      IF (.NOT. (NEWTON .OR. HSFILE .OR. INITHS .OR. MODHES .OR.
     &     INMDHS .OR. CMBMOD .OR. INRDHS .OR.
     &     (EVLINI .GT. -0.9D0))) THEN
         IF (SADDLE) THEN
            INITHS = .TRUE.
            NOAUX  = .TRUE.
         ELSE
            INMDHS = .TRUE.
         END IF
      END IF
C
      IF (HSFILE) THEN
         IF (NEWTON .OR. QUADSD) THEN
            NINFO = NINFO + 1
            WRITE(LUPRI,'(A)')
     &           'INFO: .HESFIL only has effect when a 1st order ' //
     &           'method has been specified => Keyword ignored.'
         ELSE
            WRITE(LUPRI,'(A)')
     &           ' Initial Hessian will be read from file.'
         END IF
      END IF
      IF (INITHS) THEN
         IF (NEWTON .OR. QUADSD) THEN
            NINFO = NINFO + 1
            WRITE(LUPRI,'(A)')
     &           'INFO: .INITHE only has effect when 1st order ' //
     &           'method has been specified => Keyword ignored.'
         ELSE IF (HSFILE) THEN
            NINFO = NINFO + 1
            WRITE(LUPRI,'(A)')
     &           'INFO: .INITHE has no effect when .HESFIL ' //
     &           'has been specified => Keyword ignored.'
            INITHS = .FALSE.
         ELSE
            WRITE(LUPRI,'(A)')
     &           ' Initial Hessian will be calculated.'
         END IF
      END IF
      IF (MODHES) THEN
         IF (NEWTON .OR. QUADSD) THEN
            NINFO = NINFO + 1
            WRITE(LUPRI,'(A)')
     &           'INFO: .MODHES only has effect when 1st order ' //
     &           'method has been specified => Keyword ignored.'
            MODHES = .FALSE.
         ELSE IF (HSFILE) THEN
            NINFO = NINFO + 1
            WRITE(LUPRI,'(A)')
     &           'INFO: .MODHES has no effect when .HESFIL ' //
     &           'has been specified => Keyword ignored.'
            MODHES = .FALSE.
         ELSE
            WRITE(LUPRI,'(A)')
     &           ' An approximate model Hessian will be used.',
     &           ' The model Hessian parameters ' //
     &           'of Roland Lindh will be used.'
            IF (.NOT. DELINT) REDINT = .TRUE.
         END IF
      END IF
      IF (CMBMOD) THEN
         IF (NEWTON .OR. QUADSD) THEN
            NINFO = NINFO + 1
            WRITE(LUPRI,'(A)')
     &           'INFO: .CMBMOD only has effect when 1st order ' //
     &           'method has been specified => Keyword ignored.'
            CMBMOD = .FALSE.
         ELSE IF (HSFILE) THEN
            NINFO = NINFO + 1
            WRITE(LUPRI,'(A)')
     &           'INFO: .CMBMOD has no effect when .HESFIL ' //
     &           'has been specified => Keyword ignored.'
            CMBMOD = .FALSE.
         ELSE
            WRITE(LUPRI,'(A)')
     &           ' An approximate model Hessian will be used,',
     &           ' with the model Hessian parameters of Roland Lindh.',
     &           ' The Hessian will be updated ' //
     &             'through a combination of a calculated',
     &           ' model Hessian and a BFGS update of the last Hessian.'
            IF (.NOT. DELINT) REDINT = .TRUE.
         END IF
      END IF
      IF (INMDHS) THEN
         IF (NEWTON .OR. QUADSD) THEN
            NINFO = NINFO + 1
            WRITE(LUPRI,'(A)')
     &           'INFO: .INIMOD only has effect when 1st order ' //
     &           'method has been specified => Keyword ignored.'
            INMDHS = .FALSE.
         ELSE IF (HSFILE) THEN
            NINFO = NINFO + 1
            WRITE(LUPRI,'(A)')
     &           'INFO: .INIMOD has no effect when .HESFIL ' //
     &           'has been specified => Keyword ignored.'
            INMDHS = .FALSE.
         ELSE IF (INITHS) THEN
            NINFO = NINFO + 1
            WRITE(LUPRI,'(A)')
     &           'INFO: .INIMOD has no effect when .INITHE ' //
     &           'has been specified => Keyword ignored.'
            INRDHS = .FALSE.
         ELSE
            WRITE(LUPRI,'(A)')
     &           ' Model Hessian will be used as initial Hessian.',
     &           ' The model Hessian ' //
     &             'parameters of Roland Lindh will be used.'
         END IF
      END IF
      IF (INRDHS) THEN
         IF (NEWTON .OR. QUADSD) THEN
            WRITE(LUPRI,'(A)')
     &           ' .INIRED only has effect when 1st order ' //
     &           'method has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
            INRDHS = .FALSE.
         ELSE IF (HSFILE) THEN
            WRITE(LUPRI,'(A)')
     &           ' .INIRED has no effect when .HESFIL ' //
     &           'has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
            INRDHS = .FALSE.
         ELSE IF (INITHS) THEN
            WRITE(LUPRI,'(A)')
     &           ' .INIRED has no effect when .INITHE ' //
     &           'has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
            INRDHS = .FALSE.
         ELSE IF (INMDHS) THEN
            WRITE(LUPRI,'(A)')
     &           ' .INIRED has no effect when .INIMOD ' //
     &           'has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
            INRDHS = .FALSE.
         ELSE
            WRITE(LUPRI,'(A)') ' Initial Hessian will be diagonal in '
     &           // 'internal coordinates.'
         END IF
      END IF
      IF (EVLINI .GT. 0.0D0) THEN
         IF (NEWTON .OR. QUADSD) THEN
            WRITE(LUPRI,'(A)')
     &           ' .INITEV only has effect when 1st order ' //
     &           'method has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
         ELSE IF (HSFILE) THEN
            WRITE(LUPRI,'(A)')
     &           ' .INITEV has no effect when .HESFIL ' //
     &           'has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
         ELSE IF (INITHS) THEN
            WRITE(LUPRI,'(A)')
     &           ' .INITEV has no effect when initial ' //
     &           'Hessian is calculated.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
         ELSE IF (INMDHS) THEN
            WRITE(LUPRI,'(A)')
     &           ' .INITEV has no effect when .INIMOD ' //
     &           'has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
            EVLINI = -1.0D0
         ELSE
            WRITE(LUPRI,'(A,F10.6)') ' Initial diagonal ' //
     &           'Hessian will have elements equal to: ', EVLINI
         END IF
      ELSE IF (.NOT. (DELINT .OR. REDINT .OR. INRDHS .OR. INMDHS)) THEN
         EVLINI = 1.0D0
      END IF
      WRITE(LUPRI,*)
      IF (NOHSWR) THEN
         WRITE(LUPRI,'(A/)')
     &        ' Current Hessian will not be written out to DALTON.HES.'
      END IF
      IF (REJINI) THEN
         IF (NEWTON .OR. QUADSD) THEN
            WRITE(LUPRI,'(A)')
     &           ' .REJINI only has effect when 1st order ' //
     &           'method has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
            WRITE(LUPRI,*)
         ELSE
            WRITE(LUPRI,'(A)') ' Hessian will be reinitialized' //
     &           ' after rejected steps.'
            WRITE(LUPRI,*)
         END IF
      END IF
      IF (GRDINI) THEN
         IF (NEWTON .OR. QUADSD) THEN
            WRITE(LUPRI,'(A)')
     &           ' .GRDINI only has effect when 1st order ' //
     &           'method has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
            WRITE(LUPRI,*)
         ELSE
            WRITE(LUPRI,'(A)') ' Hessian will be reinitialized' //
     &           ' when the norm of the gradient increases.'
            WRITE(LUPRI,*)
         END IF
      END IF
C
      IF (IPRINT .NE. IPDEF) THEN
         WRITE(LUPRI,'(A,I10)') ' Print level in OPTIMI  :',IPRINT
      END IF
      IF (ITRMAX .NE. ITRMX) THEN
         WRITE(LUPRI,'(A,I10)') ' Maximum # of iterations:',ITRMAX
      END IF
      IF (MAXREJ .NE. MXREJ) THEN
         WRITE(LUPRI,'(A,I10)') ' Max # of rejected steps:',MAXREJ
      END IF
      IF (VRML) THEN
         WRITE(LUPRI,'(A)') ' VRML-files for initial and final ' //
     &        'geometries will be created.'
         IF (VRBOND) WRITE(LUPRI,'(A)')
     &        ' Bonds will be drawn between nearby atoms.'
         IF (VREIGV) WRITE(LUPRI,'(A)')
     &        ' Eigenvectors will be visualized.'
      ELSE 
         IF (VRBOND) WRITE(LUPRI,'(A)')
     &        ' .VR-BOND only has effect when .VRML is specified.' //
     &        ' Keyword ignored.'
         IF (VREIGV) WRITE(LUPRI,'(A)')
     &        ' .VR-EIG only has effect when .VRML is specified.' //
     &        ' Keyword ignored.'
      END IF
      IF (DOPRE) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,'(A)') ' Preoptimization will be' //
     &        ' performed with the basis set(s):'
         DO 244 I = 1, NUMPRE
            WRITE(LUPRI,'(A,A60)') '      ',PREBTX(I)
 244     CONTINUE
      END IF
      IF (DOSPE) THEN
         WRITE(LUPRI,'(A)') ' Single point energy will be' //
     &        ' calculated using the basis:'
         WRITE(LUPRI,'(A,A60)') '      ',SPBSTX
      END IF
      IF (TRSTCH) THEN
         WRITE(LUPRI,'(/A/A/,(A,F10.4))')
     &        ' Restricted step control parameters',
     &        ' ----------------------------------',
     &        ' Initial trust radius   :',TRSTRA,
     &        ' Trust radius increment :',TRSTIN,
     &        ' Trust radius decrement :',TRSTDE,
     &        ' Bad prediction ratio   :',RTENBD,
     &        ' Good prediction ratio  :',RTENGD,
     &        ' Rejection ratio, low   :',RTRJMN,
     &        ' Rejection ratio, high  :',RTRJMX
      END IF
      WRITE(LUPRI,*)
      IF (NOTRST) THEN
         WRITE(LUPRI,*) 'No trust region will be used for steps.'
         WRITE(LUPRI,*)
      END IF
      IF (NOBRKS) THEN
         WRITE(LUPRI,*)
     &        'Symmetry will not be broken during optimization.'
         WRITE(LUPRI,*)
      END IF
      IF (THRSYM .NE. MIN(5.0D-3,SQRT(DEFTHE))) THEN
         WRITE(LUPRI,'(A,1P,D13.2)')
     &        ' Threshold for symmetry-break set to       : ', THRSYM
         IF ((THRSYM .LT. D0) .OR. (THRSYM .GT. 0.1D0)) THEN
            WRITE(LUPRI,*) 'Threshold negative or larger than 0.1'
            WRITE(LUPRI,'(A,1P,D13.2/)')
     &           ' Threshold reset to:', MIN(5.0D-3,SQRT(DEFTHE))
            THRSYM = MIN(5.0D-3,SQRT(DEFTHE))
         END IF
      END IF
      IF (.NOT. SADDLE) THEN
         IF (GDIIS) THEN
            WRITE(LUPRI,*) 'Geometrical DIIS will be used ' //
     &           'to control step.'
            RATFUN = .FALSE.
            TRSTRG = .FALSE.
         ELSE IF (RATFUN) THEN
            WRITE(LUPRI,*) 'Rational function method will be used ' //
     &           'to control step.'
            WRITE(LUPRI,*)
            TRSTRG = .FALSE.
         ELSE
            WRITE(LUPRI,*) 'Trust region method will be used ' //
     &           'to control step (default).'
            WRITE(LUPRI,*)
            TRSTRG = .TRUE.
         END IF
      ELSE
         IF (GDIIS) THEN
            WRITE(LUPRI,*) 'Geometrical DIIS not suitable for ' //
     &           'saddle point optimization. Using image function.'
            GDIIS = .FALSE.
            RATFUN = .FALSE.
            TRSTRG = .TRUE.
         ELSE IF (RATFUN) THEN
            WRITE(LUPRI,*) 'Partitioned rational function method ' //
     &           'will be used to control step.'
            TRSTRG = .FALSE.
         ELSE
            WRITE(LUPRI,*) 'Image function method will be used ' //
     &           'to control step (default).'
            TRSTRG = .TRUE.
         END IF
         IF (NSPMOD .LT. 0) THEN
            WRITE(LUPRI,*) 'The eigenvector corresponding to the ' //
     &           'lowest non-zero eigenvalue is chosen'
            WRITE(LUPRI,*) 'as reaction mode (default).'
         ELSE
            WRITE(LUPRI,'(A,I3,A)') 'Eigenvector #',NSPMOD,
     &           ' will be used as reaction mode.'
         END IF
      END IF
      IF (LNSRCH) THEN
         IF (SADDLE) THEN
            WRITE(LUPRI,*)
     &           'Line search disabled because saddle point is sought.'
            WRITE(LUPRI,*)
            LNSRCH = .FALSE.
         ELSE
            WRITE(LUPRI,*) 'Partial line search with bound quartic ' //
     &           'polynomial will be employed.'
            WRITE(LUPRI,*)
         END IF
      END IF
C
      IF (BAKER) THEN
         WRITE(LUPRI,'(/A/)')
     &        ' Baker''s convergence criteria will be used.'
#ifdef VAR_LINSCA
         WRITE(LUPRI,'(/A)')' ERROR! Can not use .BAKER when ' //
     &        'running linsca. Only NATNRM can be used!'
         CALL QUIT('Conflicting specifications of ' //
     &        'convergence setup in '//WORD1)
#endif     
         IF (NATNRM) THEN
            WRITE(LUPRI,'(/A/)')
     &        ' INFO:  .NATNRM is ignored when .BAKER requested!'
            NATNRM = .FALSE.
         END IF
      ELSE IF (NATNRM) THEN
C
C     Tweaked convergence criterias more suitable for large systems (linsca)
C
         WRITE(LUPRI,'(/A)') ' New thresholding scheme for ' //
     &        'geometry convergence will be used'
         IF ((ICNLVL .GT. 0) .AND. ((GRDTHR .NE. DEFTH2)
     &        .OR. (THRSTP .NE. DEFTH2))) THEN
            WRITE(LUPRI,'(/,A)')' WARNING! Conflicting specifications'//
     &           ' of convergence thresholds in '//WORD1
            WRITE(LUPRI,'(A/A/)') ' Please use *either* ''GRADIE'' '//
     &           'and ''STEP T'' *or* one of the keywords',
     &           ' ''.VLOOSE'', ''.LOOSE'', ''.TIGHT'' and ''.VTIGHT.'''
            CALL QUIT('Conflicting specifications of ' //
     &           'convergence threshold in '//WORD1)
         ELSE IF (ICNLVL .GT. 0) THEN
C
C     VLOOSE, LOOSE, TIGHT and VTIGHT are simply predefined threshold for
C     the gradient convergence criteria
C            
            IF (ICNLVL .EQ. 1) THEN
               GRDTHR = 1.D-3
               THRSTP = 1.D-3
            ELSE IF (ICNLVL .EQ. 2) THEN
               GRDTHR = 1.D-4
               THRSTP = 1.D-4
            ELSE IF (ICNLVL .EQ. 4) THEN
               GRDTHR = 1.D-6
               THRSTP = 1.D-6
            ELSE IF (ICNLVL .EQ. 5) THEN
               GRDTHR = 1.D-7
               THRSTP = 1.D-7
            END IF
         ELSE
            IF (GRDTHR .NE. DEFTH2) THEN
               IF ((GRDTHR .LE. D0) .OR. (GRDTHR .GT. 0.1D0))
     &              WRITE(LUPRI,*) 'RMS gradient threshold negative ' //
     &              'or larger than 0.1, value is reset!'
            ELSE
C     Default value for new scheme
               GRDTHR = 1.D-5
            END IF
            IF (THRSTP .NE. DEFTH2) THEN
               IF ((THRSTP .LE. D0) .OR. (THRSTP .GT. 0.1D0))
     &              WRITE(LUPRI,*) 'RMS step threshold negative ' //
     &              'or larger than 0.1, value is reset!'
            ELSE
C     Default value for new scheme
               THRSTP = 1.D-5
            END IF
         END IF
         THGRMX = GRDTHR*1.5D0
         THSTMX = THRSTP*1.5D0
         WRITE(LUPRI,'(A)') ' --------------' //
     &        '-----------------------------------------------'
         WRITE(LUPRI,'(A,1P,D13.2)')
     &     ' Root-mean-square gradient threshold set to    : ', GRDTHR,
     &     ' Maximum gradient element threshold set to     : ', THGRMX,
     &     ' Root-mean-square step threshold set to        : ', THRSTP,
     &     ' Maximum step element threshold set to         : ', THSTMX
      ELSE
C
C     The old (and still default) convergence setup
C
         IF (ICNLVL .GT. 0) THEN
            WRITE(LUPRI,*) 'Keywords VLOOSE, LOOSE, TIGHT and VTIGHT '//
     &           'currently only have effect',
     &           'when .NATNRM is specified.'
            WRITE(LUPRI,*)
         END IF
         IF (GRDTHR .NE. DEFTH2) THEN
            WRITE(LUPRI,'(A,1P,D13.2)')
     &           ' Convergence threshold for gradient set to : ', GRDTHR
            IF ((GRDTHR .LE. D0) .OR. (GRDTHR .GT. 0.1D0)) THEN
               WRITE(LUPRI,'(A,1P,D13.2)')
     &           ' Threshold negative or larger than 0.1,'//
     &           ' threshold reset to:', DEFTH2
               GRDTHR = DEFTH2
            END IF
         END IF
         IF (THRERG .NE. DEFTHE) THEN
            WRITE(LUPRI,'(A,1P,D13.2)')
     &           ' Convergence threshold for energy set to   : ', THRERG
            IF ((THRERG .LE. D0) .OR. (THRERG .GT. 0.1D0)) THEN
               WRITE(LUPRI,'(A,1P,D13.2)')
     &           ' Threshold negative or larger than 0.1,'//
     &           ' threshold reset to:', DEFTHE
               THRERG = DEFTHE
            END IF
         END IF
         IF (THRSTP .NE. DEFTH2) THEN
            WRITE(LUPRI,'(A,1P,D13.2)')
     &           ' Convergence threshold for step set to     : ', THRSTP
            IF ((THRSTP .LE. D0) .OR. (THRSTP .GT. 0.1D0)) THEN
               WRITE(LUPRI,'(A,1P,D13.2)')
     &           ' Threshold negative or larger than 0.1,'//
     &           ' threshold reset to:', DEFTH2
               THRSTP = DEFTH2
            END IF
         END IF
         IF (ICONDI .NE. 2) THEN
            WRITE(LUPRI,'(A,I3)')
     &           ' Number of convergence criteria set to     : ', ICONDI
            IF ((ICONDI .LT. 1) .OR. (ICONDI .GT. 3)) THEN
               WRITE(LUPRI,*) 'Acceptable values are 1, 2 and 3.'
               ICONDI = 2
               WRITE(LUPRI,'(A,I3)')
     &              ' Number reset to:', ICONDI
            END IF
         END IF
      END IF
C
      IF (DISPLA .NE. DEFDSP) THEN
C NB TODO hjaaj Jan09: DISPLA is never used !!!!! Remove ???
         WRITE (LUPRI,'(A,1P,D13.2)')
     &        ' In case of numerical gradients, displacements are : ',
     &        DISPLA
      END IF
      IF (CONOPT) THEN
         WRITE(LUPRI,'(/A)')
     &        ' Constrained optimization has been requested.'
         NOAUX = .TRUE.
         IF (.NOT. REDINT) THEN
            WRITE(LUPRI,'(A)') ' WARNING! Constrained optimizations '
     &           // 'can only be used in conjunction with'
            WRITE(LUPRI,'(A)') ' redundant internal coordinates!. '
     &           // 'Keyword ignored.'
            CONOPT = .FALSE.
         ELSE
            WRITE(LUPRI,'(A)')
     &           ' The following coordinate numbers will be held fixed '
     &           // 'during the optimization:'
            DO 333 I = 1, 8*MXCENT
               IF (ICNSTR(I) .EQ. 1) WRITE(LUPRI,*) '    Coordinate #',I
 333        CONTINUE
            IF (.NOT. NOBRKS) THEN
               NOBRKS = .TRUE.
               WRITE(LUPRI,*)
               WRITE(LUPRI,*) 'Symmetry will not be broken during '
     &              // 'a constrained optimization.'
            END IF
            WRITE(LUPRI,*)
         END IF
      END IF
      IF (REMCRD) THEN
         WRITE(LUPRI,'(A)')
     &        ' Removal of coordinates has been requested.'
         IF (.NOT. REDINT) THEN
            WRITE(LUPRI,'(A)') ' WARNING! Only internal coordinates '
     &           // 'can be removed! Keyword ignored.'
            REMCRD = .FALSE.
         ELSE
            WRITE(LUPRI,'(A)')
     &           ' The following coordinate numbers will be removed: '
            DO 335 I = 1, 8*MXCENT
               IF (ICNSTR(I) .EQ. 2) WRITE(LUPRI,*) '    Coordinate #',I
 335        CONTINUE
            WRITE(LUPRI,*)
         END IF
      END IF
      IF (ADDCRD) THEN
         WRITE(LUPRI,'(A)')
     &        ' Addition of coordinates has been requested.'
         IF (.NOT. REDINT) THEN
            WRITE(LUPRI,'(A)') ' WARNING! Only internal coordinates '
     &           // 'can be added! Keyword ignored.'
            ADDCRD = .FALSE.
         ELSE
            WRITE(LUPRI,'(A)')
     &           ' The following coordinates will be added:'
            DO 345 I = 1, IADDCR(0,1)
               IF (IADDCR(I,4) .GT. 0) THEN
                  WRITE(LUPRI,'(A,4I6)') ' Dihedral between atoms: ',
     &            IADDCR(I,1), IADDCR(I,2), IADDCR(I,3), IADDCR(I,4)
               ELSE IF (IADDCR(I,3) .GT. 0) THEN
                  WRITE(LUPRI,'(A,3I6)') ' Angle between atoms   : ',
     &            IADDCR(I,1), IADDCR(I,2), IADDCR(I,3)
               ELSE IF (IADDCR(I,2) .GT. 0) THEN
                  WRITE(LUPRI,'(A,2I6)') ' Bond between atoms    : ',
     &            IADDCR(I,1), IADDCR(I,2)
               ELSE
                  CALL QUIT(' Error detected in the '
     &                 // 'specification of additional coordinates!')
               END IF
 345        CONTINUE
            WRITE(LUPRI,*)
         END IF
      END IF
C
C     Experimental atom freezing
C
      IF (IFREEZ(0) .GT. 0) THEN
         WRITE(LUPRI,'(A)')
     &        ' Freezing of atoms has been requested.'
         WRITE(LUPRI,'(A)') ' WARNING: .FREEZE is an ' //
     &        'experimental feature, use it at your own risk...'
         IF (.NOT. CONOPT) THEN
            WRITE(LUPRI,'(A)') ' WARNING: Freezing of atoms is really '
     &           // 'intended to be used for constrained optimizations!'
         END IF
         WRITE(LUPRI,'(A,I4,A)') ' The following',IFREEZ(0),
     &        ' atoms will be held frozen: '
C     Do a simple bubble sort
         DO 351 I = 1, IFREEZ(0)-1
            DO 352 J = I+1, IFREEZ(0)
               IF (IFREEZ(I) .GT. IFREEZ(J)) THEN
                  ITMP = IFREEZ(I)
                  IFREEZ(I) = IFREEZ(J)
                  IFREEZ(J) = ITMP
               END IF
 352        CONTINUE
            WRITE(LUPRI,*) '    Atom #      ', IFREEZ(I)
 351     CONTINUE
         WRITE(LUPRI,*) '    Atom #      ', IFREEZ(IFREEZ(0))
         IF (ITRFRZ .GT. 0) WRITE(LUPRI,'(A,I4,A)')
     &        '     Atoms will be kept frozen for',ITRFRZ, ' iterations'
      END IF
C
C     VRML options
C
      IF (VRML) THEN
         WRITE(LUPRI,'(A)')
     &        ' VRML-file of geometry will be created.'
         IF (VRBOND) WRITE(LUPRI,'(A)')
     &        ' Bonds will be drawn between nearby atoms.'
         IF (VRCORD) WRITE(LUPRI,'(A)')
     &        ' Coordinate axes will be drawn.'
         IF (VRML_SYM) THEN
            WRITE(LUPRI,'(A/A)')
     &        ' Symmetry elements will be visualized.',
     &        ' Please note that symmetry should NOT be explicitly' //
     &        ' specified in the input file.'
         END IF
         IF (VREIGV) WRITE(LUPRI,'(A)')
     &        ' Eigenvectors will be visualized.'
         IF (VRVIBA) WRITE(LUPRI,'(A)')
     &        ' Vibrational modes will be visualized.'
      END IF
C
C     Geometry analysis
C
      IF (GEOALL) WRITE(LUPRI,'(/A/)')
     &        ' Geometry analysis will be printed for each iteration.'
C
C     Adjusted accuracy
C
      IF (ISTBLZ .GT. 0) THEN
         WRITE(LUPRI,'(A,I2)')
     &        ' Geometries are stabilized by ignoring all digits '
     &        // 'beyond digit number ', ISTBLZ
         WRITE(LUPRI,'(A/A)')
     &        ' WARNING: This is an experimental feature!', 
     &        ' Use it at your own risk, ' // 
     &        'the result may or may not be meaningful...'
      END IF
      RETURN
      END

C  /* Deck optmin */
      SUBROUTINE OPTMIN(WORK,LWORK,WRKDLM)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "infpar.h"
#include "inftap.h"
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
#endif
#include "gnrinf.h"
#include "cbirea.h"
#include "optinf.h"
      DIMENSION WORK(LWORK)
      LOGICAL TSTIN
      CALL QENTER('OPTMIN')
      IF (NEWTON) THEN
         CALL TITLER('2nd Order Geometry Optimization','>',110)
      ELSE
         CALL TITLER('1st Order Geometry Optimization','>',110)
      END IF
C
C     *** Need to figure out which analytical derivatives ***
C     *** we have. If numerical derivatives, nomove=.true.***
C        
c      CALL FNDANA(NAORDR)
c      .... PROBLEM: WAVE FUNCTION TYPE NOT KNOWN YET !!!
c      IF (NAORDR.EQ.0) THEN
c         NOMOVE = .TRUE.
c         NMWALK = .TRUE.
c      END IF
C
C     The number of Cartesian and redundant internal coordinates
C     are counted to reduce the memory requirement.
C     The input processing in Hermit has to be run to determine
C     these coordinates.
C     The call is also necessary for PREOPT with different basis, see
C     INIPRE.
C     Using "INPTST true" to only run the necessary parts of READIN.
C
      TSTIN  = INPTST
      INPTST = .TRUE.
      IPRUSR_SAV = IPRUSR
      IPREAD_SAV = IPREAD
      IPRUSR = -10
      IPREAD = -10
C     Third parameter .FALSE. in call indicates that LUONEL will NOT be written.
      write (lupri,*) 'Read .mol file to determine number of atoms:'
      CALL READIN(WORK(2),LWORK-2,.FALSE.)
      INPTST = TSTIN
      IPRUSR = IPRUSR_SAV
      IPREAD = IPREAD_SAV
      HRINPC = .FALSE.
C
      IF (FINDRE .OR. DELINT .OR. REDINT .OR. INRDHS .OR. INMDHS) THEN
         KATARR = 2
         KMTMP  = KATARR + 8*MXCENT
         KWRK   = KMTMP  + MXCENT*MXCENT
         LWRK = LWORK  - KWRK + 1
         IF (KWRK+8*MXCENT .GT. LWORK)
     &        CALL STOPIT('OPTMIN',' ',KWRK,LWORK)
         IPSAVE = IPRINT
         IF (.NOT. FINDRE) IPRINT = -1
C     Save ICNSTR vector (coordinate constraints and removals)
         DO 10 I = 1, 8*MXCENT
            WORK(KWRK+I) = ICNSTR(I)
 10      CONTINUE
         CALL FNDRED(WORK(KATARR),WORK(KMTMP))
         DO 15 I = 1, 8*MXCENT
            ICNSTR(I) = WORK(KWRK+I)
 15      CONTINUE
         IPRINT = IPSAVE
         MXRINT = MAX(MAX(IINTCR, ICRTCR), 8)
C
C     If only determination of internal coordinates has been requested,
C     the program stops.
C
         IF (FINDRE .AND. (.NOT. VISUAL)) THEN
            CALL TITLER('Internal coordinates determined','#',103)
            GO TO 9999
         END IF
      ELSE
         KATARR = 2
         KWRK   = KATARR + 8*MXCENT
         LWRK = LWORK  - KWRK + 1
         IF (KWRK .GT. LWORK) CALL STOPIT('OPTMIN',' ',KWRK,LWORK)
         CALL ATMINI(WORK(KATARR),IATOM,.FALSE.)
         MXRINT = MAX(3*IATOM, 8)
      END IF
C
C     If visualization is requested, the program stops afterwards.
C
      IF (VISUAL) THEN
         CALL VISULZ(WORK,LWORK,WRKDLM)
         CALL TITLER('End of visualization','#',103)
         GO TO 9999
      END IF
C
      KGRAD  = 1
      KSTEP  = KGRAD  + MXCOOR
      KSCLVC = KSTEP  + MXCOOR
      KGRDOL = KSCLVC + MXCOOR
      KGRDMT = KGRDOL + MXRINT
      KSTPMT = KGRDMT + MXRINT*25
      KHESOL = KSTPMT + MXRINT*25
      KGRDAR = KHESOL + MXRINT*MXRINT
      KSTPAR = KGRDAR + MXRINT*25
      KINFO  = KSTPAR + MXRINT*25
      KWILBM = KINFO  + (ITRMAX + 1)*6
      KBMTRN = KWILBM + MXRINT*MXCOOR
      KHSINT = KBMTRN + MXRINT*MXRINT
      KVCMOD = KHSINT + MXRINT*MXRINT
      KWRK   = KVCMOD + MXCOOR
      LWRK   = LWORK  - KWRK + 1
      IF (KWRK .GT. LWORK) CALL STOPIT('OPTMIN',' ',KWRK,LWORK)
      CALL RUNOPT(WRKDLM,MXRINT,WORK(KGRAD),WORK(KSTEP),
     &     WORK(KSCLVC),WORK(KGRDOL),WORK(KGRDMT),WORK(KSTPMT),
     &     WORK(KHESOL),WORK(KGRDAR), WORK(KSTPAR),WORK(KINFO),
     &     WORK(KWILBM),WORK(KBMTRN),WORK(KHSINT),WORK(KVCMOD),
     &     WORK(KWRK),LWRK)
C
C     Commented out this part, the finalizing is done in the main
C     program (dalton.F) for all kinds of calculations - vb
C
C#if defined (VAR_MPI)
C      IF (NODTOT .GT. 1 .AND. MYNUM .EQ. 0) THEN
C         IJOB = 0
C         CALL MPI_BCAST(IJOB,1,my_MPI_INTEGER,NODEID(MASTER),
C     &                  MPI_COMM_WORLD,IERR)
C      END IF
C#endif      
 9999 CALL QEXIT('OPTMIN')
      RETURN
      END

C  /* Deck runopt */
      SUBROUTINE RUNOPT(WRKDLM,MXRCRD,EGRAD,CSTEP,SCLVEC,GRDOLD,
     &     GRDMAT,STPMAT,HESOLD,GRDARR,STPARR,GEINFO,WILBMT,
     &     BMTRAN,HESINT,VECMOD,WORK,LWORK)
#include "implicit.h"
#include "dummy.h"
#include "codata.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "optinf.h"
#include "cbirea.h"
#include "molinp.h"
#include "abainf.h"
#include "gnrinf.h"
#include "inftap.h"
#include "maxaqn.h"
#include "symmet.h"
#include "nuclei.h"
#include "taymol.h"
#include "prkoor.h"
#include "trkoor.h"
#include "cbiwlk.h"
#include "molde.h"
#include "numder.h"
      LOGICAL EXHER, EXSIR, EXABA, MINEND, INDXOK, STATPO, ACTIVE
      LOGICAL REJGEO, TRU, FAL, TMPLOG, NEWSTP
      CHARACTER TMPLIN*80, WRDRSP*7
      DIMENSION EGRAD(MXCOOR), CSTEP(MXCOOR), SCLVEC(MXCOOR)
      DIMENSION GRDOLD(MXRCRD), GRDMAT(25,MXRCRD)
      DIMENSION STPMAT(25,MXRCRD), HESOLD(MXRCRD,MXRCRD)
      DIMENSION GRDARR(25,MXRCRD), STPARR(MXRCRD,MXRCRD)
C emg
      CHARACTER LABEL1*8, STHELP*10
C emg
C
C     The array geinfo contains optimization information for each
C     iteration. The first index is the iteration, the second gives
C     the property:   1  -  Energy
C                     2  -  Gradient norm
C                     3  -  Index of Hessian
C                           (a negative index indicates symmetry break)
C                     4  -  Step length
C                     5  -  Trust radius
C                     6  -  # rejected steps
C
      DIMENSION GEINFO(0:ITRMAX,6)
      DIMENSION WILBMT(MXRCRD,MXCOOR), BMTRAN(MXRCRD,MXRCRD)
      DIMENSION HESINT(MXRCRD,MXRCRD), VECMOD(MXCOOR)
      DIMENSION WORK(LWORK)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      CALL QENTER('RUNOPT')
C
C     Initialization of variables.
C
      THRLDP = 1.0D-4
      THRIND = 5.0D-4
      TOLST  = 1.0D-5
      IPRWLK = IPRINT
      EXHER  = .FALSE.
      EXSIR  = .FALSE.
      EXABA  = .FALSE.
      RDINPC = .FALSE.
      RDMLIN = .FALSE.
      CALL DZERO(GEINFO,(ITRMAX+1)*6)
      CALL DZERO(GRDARR,25*MXRCRD)
      CALL DZERO(STPARR,25*MXRCRD)
      RSTARR = .FALSE.
      ACTIVE = .FALSE.
      KEPTIT = 0
      GEINFO(0,5) = TRSTRA
C
C     Perform preoptimization if requested.
C
      IF (DOPRE) CALL INIPRE(WORK,LWORK,WRKDLM)
C
C     Part of the WORK-array is used for storage.
C
      KEHESS = 1
      KALHES = KEHESS + MXCOOR*MXCOOR
      KWRK1  = KALHES + MXCOOR*MXCOOR
      IF (KWRK1 .GT. LWORK) CALL STOPIT('RUNOPT',' ',KWRK1,LWORK)
      LWRK1  = LWORK - KWRK1 + 1
C
C     Calculate energy, gradient and Hessian for second order method and
C     first order method with initial Hessian.
C
      IF (NEWTON .OR. INITHS) THEN
         CALL GTHESS(EGRAD,WORK(KEHESS),WORK(KALHES),EXHER,EXSIR,EXABA,
     &        WORK(KWRK1),LWRK1,WRKDLM)
C
C     First order methods only require the energy and the gradient.
C
      ELSE
         CALL GTGRAD(EGRAD, EXHER,EXSIR,EXABA,
     &        WORK(KWRK1),LWRK1,WRKDLM)
      END IF
C
C     Make VRML-file of initial geometry if requested.
C
      IF (VRML) THEN
         KATARR = KWRK1
         KEVEC  = KATARR + 8*MXCENT
         KEVC1  = KEVEC  + MXCOOR*MXCOOR
         KEVC2  = KEVC1  + MXCOOR
         KWRK2  = KEVC2  + MXCOOR
         LWRK   = LWORK  - KWRK2 + 1
         IF (KWRK2 .GT. LWORK) CALL STOPIT('RUNOPT',' ',
     &        KWRK2,LWORK)
         CALL MKVRML(.FALSE.,WORK(KATARR),MXCOOR,WORK(KEVEC),
     &        WORK(KEVC1),WORK(KEVC2))
      END IF
C
C     Save initial geometry and energy to MOLDEN file
C
      IF (MOLDEN) CALL MOGECON(.FALSE.,ENERGY)
C
C     Count coordinates.
C
      CALL CNTCRD
C
C     We allocate more of the WORK-array
C
      MX2CRD = MAX(MXCOOR,MXRCRD)
      KBMINV = KWRK1
      KPJINM = KBMINV + MXRCRD*MXCOOR
      KEVEC  = KPJINM + MXRCRD*MXRCRD
      KCONMT = KEVEC  + MX2CRD*MX2CRD
      KTEMP1 = KCONMT + MXRCRD*MXRCRD
      KTEMP2 = KTEMP1 + MX2CRD*MX2CRD
      KTEMP3 = KTEMP2 + MX2CRD*MX2CRD
      KTEMP4 = KTEMP3 + MX2CRD*MX2CRD
      KTEMP5 = KTEMP4 + MX2CRD*MX2CRD
      KTEMP6 = KTEMP5 + MX2CRD*MX2CRD
      KTEMP7 = KTEMP6 + MX2CRD*MX2CRD
      KTEMP8 = KTEMP7 + MX2CRD*MX2CRD
      KTEMP9 = KTEMP8 + MX2CRD*MX2CRD
      KWRK2  = KTEMP9 + MX2CRD*MX2CRD
      IF (KWRK2 .GT. LWORK) CALL STOPIT('RUNOPT',' ',KWRK2,LWORK)
      LWRK2  = LWORK - KWRK2 + 1
C
C     Check if redundant internal coordinates should be used.
C
      IF (DELINT .OR. REDINT .OR. INRDHS .OR. INMDHS) THEN
         CALL INIRED(MXRCRD,MX2CRD,WILBMT,BMTRAN,WORK(KBMINV),
     &        WORK(KPJINM),WORK(KTEMP1),WORK(KTEMP2),
     &        WORK(KTEMP3),WORK(KTEMP4),WORK(KTEMP5),WORK(KTEMP6),
     &        WORK(KWRK2),LWRK2)
      END IF
      IF (DELINT .OR. REDINT) THEN
         NCRDHS = IINTCR
      ELSE
         NCRDHS = NCART
      END IF
      IF (RATFUN) NCRDHS = NCRDHS + 1
C
C     aug 99 - hjaaj
C     cut down on hermit and abacus output after initial iteration
C
      IPRUSR_orig = IPRUSR
      IPREAD_orig = IPREAD
      IF (USRIPR) THEN
C        if user has asked for higher print level, no change
         IPRUSR_reduced = IPRUSR
         IPREAD_reduced = IPREAD
      ELSE
         IPRUSR_reduced = -2
         IPREAD_reduced = -2
      END IF
      IPRUSR = IPRUSR_reduced
      IPREAD = IPREAD_reduced
C
C     Initialize Hessian if first order method is used.
C
 7    CONTINUE
      IF (.NOT. NEWTON) CALL INIHES(MXRCRD,MX2CRD,SCLVEC,GRDOLD,HESOLD,
     &     WORK(KTEMP1),WORK(KTEMP2),WORK(KTEMP3),WORK(KTEMP4),
     &     WILBMT,BMTRAN,WORK(KBMINV),HESINT,WORK(KWRK2),LWRK2)
C
      CALL MINCGH(EGRAD,WORK(KEHESS),WORK(KALHES),WORK(KWRK2),LWRK2)
C
C emg ---------------------------------------------------------------
C emg WRIDER: Piece of code to Write Derivatives of Gradient and Hessian
C emg on the optimization procedure (OPTIMIZE) in Cartesian non-mass
C emg weigthed coordinates for further use in MidasCpp.
C emg -----------------------------------------------------------------
      IF (ITRMAX.EQ.0) THEN
         LUHES = -1 
         IF (NEWTON) THEN 
            KCSTRA = KWRK2
            KSCTRA = KCSTRA + (3*NUCDEP)**2
            KH2    = KSCTRA + (3*NUCDEP)**2 
            KWRK3  = KH2    + MXCOOR*MXCOOR
            IF (KWRK3 .GT.LWORK) CALL STOPIT('RUNOPT',' ',KWRK3,LWORK)
            CALL TRAHES(HESMOL,MXCOOR,WORK(KH2),WORK(KCSTRA),
     &                  WORK(KSCTRA),MXCOOR,NUCDEP*3,1)
         ENDIF 
         CALL GTGRAD(EGRAD,EXHER,EXSIR,EXABA,WORK(KWRK1),LWRK1,WRKDLM)
         CALL GPOPEN(LUHES,'midasifc.ederivs','UNKNOWN',' ',
     &            'FORMATTED', IDUMMY,.FALSE.)
         WRITE(LUHES,'(I6)') 3*NUCDEP
         WRITE(LUHES,*)
         DO 413 I = 1, 3*NUCDEP
            WRITE(LUHES,'(E25.16)') EGRAD(I)
413      CONTINUE
         IF (NEWTON) THEN 
            WRITE(LUHES,*)
            DO 411 J = 1, 3*NUCDEP
               DO 412 I = 1, 3*NUCDEP
                  WRITE(LUHES,'(E25.16)') 
     &                  WORK(KH2+(J-1)*MXCOOR+I-1)  
412            CONTINUE
411         CONTINUE
         ENDIF 
         CALL GPCLOSE(LUHES,'KEEP')

         LABEL1 = 'ALL_DONE'
         STHELP = 'THE_END   ' 
         CALL WRIPRO(0.0D0,STHELP,666,
     &               LABEL1,LABEL1,LABEL1,LABEL1,
     &               ZERO,ZERO,ZERO,1,0,0,0)

         CALL GPCLOSE(LUNDPF,'KEEP')
         CALL GPCLOSE(LUNMPF,'KEEP')
         STOP 'Manually stop after punching the Hessian and Gradient'

      ENDIF
C
C emg ---------------------------------------------------------------
C
C     Construct projection operator and use it.
C     Then diagonalize Hessian.
C
      IF (REDINT .OR. DELINT) THEN

         IF (NEWTON) CALL CGHINT(MXRCRD,MX2CRD,SCLVEC,WORK(KTEMP1),
     &        WORK(KTEMP2),WORK(KTEMP3),WORK(KTEMP4),WORK(KTEMP5),
     &        WILBMT,WORK(KBMINV),BMTRAN,HESINT,WORK(KWRK2),LWRK2)
         CALL PRJINT(MXRCRD,IINTCR,WORK(KPJINM),WORK(KCONMT),
     &        HESINT,WORK(KTEMP1),WORK(KTEMP2),WORK(KTEMP3),
     &        WORK(KTEMP4),WORK(KWRK2),LWRK2)
C
C     Note that the contents of WORK(KTEMP7) is passed on
C     from LINSRC to FNSTIN below.
C
         IF (LNSRCH .AND. RATFUN .AND. (ITRNMR .GT. 0))
     &        CALL LINSRC(IINTCR,MXRCRD,GRDINT,GRDARR(1,1),
     &        WORK(KTEMP7),STPARR(1,1),WORK(KTEMP3),WORK(KTEMP4),
     &        ACTIVE,EMOD)
         IF (RATFUN .AND. SADDLE) NCRDHS = NCRDHS - 1
         CALL DIAINT(MXRCRD,MX2CRD,NCRDHS,WORK(KEVEC),WORK(KTEMP1),
     &        WORK(KTEMP2),WORK(KTEMP3),WORK(KTEMP4),THRIND,HESINT,
     &        WORK(KTEMP5),WORK(KWRK2),LWRK2)
         IF (RATFUN .AND. SADDLE) NCRDHS = NCRDHS + 1
      ELSE
         IF (NEWTON) CALL MKSCVC(SCLVEC)
C
C     Note that the contents of WORK(KTEMP1) is passed on
C     from PROJGH to DIAHES below.
C
         CALL PROJGH(EGRAD,WORK(KEHESS),WORK(KALHES),WORK(KTEMP1),
     &        WORK(KTEMP2),WORK(KTEMP3),WORK(KTEMP4))
         IF (LNSRCH .AND. RATFUN .AND. (ITRNMR .GT. 0))
     &        CALL LINSRC(NCART,MXCOOR,EGRAD,GRDARR(1,1),CSTEP,
     &        STPARR(1,1),WORK(KTEMP3),WORK(KTEMP4),ACTIVE,EMOD)
         CALL DIAHES(MXRCRD,MX2CRD,NCRDHS,EGRAD,WORK(KEHESS),
     &        WORK(KALHES),WORK(KTEMP1),THRIND,WORK(KEVEC),WORK(KTEMP2),
     &        WORK(KTEMP3),WORK(KTEMP4),WORK(KWRK2),LWRK2)
      END IF
      GEINFO(0,1) = ENERGY
      GEINFO(0,3) = INDTOT*1.0D0
C
C     Write Hessian to file (for 1st order restarts).
C
      IF (.NOT. NOHSWR)
     &     CALL PNCHES(MXRCRD,MX2CRD,HESINT,WILBMT,BMTRAN,WORK(KTEMP1),
     &     WORK(KTEMP2),WORK(KTEMP3),WORK(KTEMP4),WORK(KWRK2),LWRK2)
C
C     Determine step, check for convergence, print output and
C     and update geometry.
C
      IREJ = 0
 755  CONTINUE
      IF (REDINT .OR. DELINT) THEN
         CALL FNSTIN(MXRCRD,MX2CRD,NCRDHS,SCLVEC,HESINT,WORK(KEVEC),
     &        WORK(KTEMP1),WORK(KTEMP2),WORK(KTEMP3),WORK(KTEMP4),
     &        WORK(KTEMP5),CSTEP,WILBMT,BMTRAN,WORK(KBMINV),GRDARR,
     &        STPARR,ACTIVE,EMOD,VECMOD,WORK(KTEMP7))
      ELSE
         CALL FNDSTP(MXRCRD,MX2CRD,NCRDHS,EGRAD,WORK(KEHESS),
     &        WORK(KEVEC),WORK(KTEMP1),WORK(KTEMP2),WORK(KTEMP3),
     &        WORK(KTEMP4),WORK(KTEMP5),CSTEP,GRDARR,STPARR,
     &        ACTIVE,EMOD,VECMOD)
      END IF
      GECONV = MINEND(MXRCRD,SCLVEC,BMTRAN,WORK(KTEMP1),WORK(KTEMP2))
C
C     If there has been a completely failed step, the geometry has
C     by default not converged.
C
      IF (ABS(GEINFO(0,6)) .GT. 1.0D-3) GECONV = .FALSE.
      CALL PRIALL(CSTEP,WORK(KTEMP1))
C
C     Save this geometry and energy to MOLDEN file
C
      IF (MOLDEN) CALL MOGECON(.FALSE.,ENERGY)
C
      NEWSTP = .FALSE.
C
C     To allow reinitialization
C
      INITHS = .FALSE.
C
      IF (.NOT. GECONV) CALL FNDGEO(CSTEP,EGRAD,WORK(KTEMP1),
     &     WORK(KTEMP2),EXHER,EXSIR,EXABA,WORK(KWRK2),LWRK2,
     &     WRKDLM,IREJ,GEINFO,NEWSTP)
      IF (NEWSTP) GOTO 755
      GEINFO(0,2) = GRADNM
      GEINFO(0,4) = STPNRM
      IF (ITRNMR .LT. ITRMAX) GEINFO(1,5) = TRSTRA
      IF (ABS(GEINFO(0,6)) .LT. 1.0D-3) THEN
         GEINFO(0,6) = IREJ*1.0D0
      ELSE
         GEINFO(0,6) = -(ABS(GEINFO(0,6))+ABS(IREJ)*1.0D0)
      END IF
      ITOTRJ = ITOTRJ + ABS(IREJ)
C
C     Determine value of the various coordinates
C
      IF (REDINT .AND. (IPRINT .GE. 1)) THEN
         CALL ATMINI(WORK(KTEMP1),IATOM,.TRUE.)
         CALL GETINT(IATOM,MXRCRD,WORK(KTEMP1),CRDINT)
         CALL HEADER('New internal coordinates',-1) 
         CALL OUTPUT(CRDINT,1,1,1,IINTCR,1,MXRCRD,1,LUPRI)
         WRITE(LUPRI,'(//)')
      END IF
C
C     If the step has failed
C
      IF (IREJ .LT. 0) THEN
         GOTO 7
C      ELSE IF (REJINI .AND. REDINT .AND. (ITOTRJ .GE. 3)) THEN
C         WRITE(LUPRI,*)'***** NOTE! *****'
C         WRITE(LUPRI,*)
C     &        'The number of dihedral angles will be reduced!'
C         CALL RREDUN
C         ITOTRJ = 0
      END IF
C
C     Make VRML-file of next geometry if requested
C
      IF (VRML) CALL MKVRML(.TRUE.,WORK(KTEMP1),MXCOOR,WORK(KEVEC),
     &        WORK(KTEMP2),WORK(KTEMP3))
C
C     Check if symmetry should be broken.
C
      IF (BRKSYM .AND. (.NOT. NOBRKS)) THEN
         IF (REDINT .OR. DELINT) THEN
            CALL DOBRKI(MXRCRD,MX2CRD,EXHER,EXSIR,EXABA,GEINFO,
     &           WORK(KEVEC),WORK(KBMINV),CSTEP,WORK(KTEMP1),
     &           WORK(KTEMP2),WORK(KTEMP3))
         ELSE
            CALL DOBRK(MXRCRD,MX2CRD,EXHER,EXSIR,EXABA,GEINFO,
     &           WORK(KEVEC),CSTEP,WORK(KTEMP1),WORK(KTEMP2))
         END IF
CChj aug 99: now new basis, reset hermit/abacus print
         IPRUSR = IPRUSR_orig
         IPREAD = IPREAD_orig
C
C     Check if preoptimization is finished
C
      ELSE IF (GECONV .AND. DOPRE .AND. (.NOT. FINPRE)) THEN
         CALL ENDPRE(EXHER,EXSIR,EXABA)
Chj aug 99: now new basis, reset hermit/abacus print
         IPRUSR = IPRUSR_orig
         IPREAD = IPREAD_orig
      END IF
C
C     DO-WHILE loop that runs until geometry has converged or
C     maximum number of iterations is reached.
C
 10   CONTINUE
      IF ((ITRNMR .LT. ITRMAX) .AND. (.NOT. GECONV)) THEN
#ifdef VAR_LINSCA
CANDREAS: reset the SAVEAOPT parameter. There might be a better place to do this, but where?
         CALL DASINI
CANDREAS END.
#endif
         ITRNMR = ITRNMR + 1
         ERGPRO = ERGPRD
         NCRD = NCRTOT
         IF (REDINT .OR. DELINT) NCRD = IINTCR
         DO 20 I = 1, NCRD
            EVALOL(I) = EVAL(I)
 20      CONTINUE
         IF (REDINT .OR. DELINT) THEN
            CALL UPGDST(IINTCR,MXRCRD,GRDARR,STPARR,GRDINT,STPINT)
         ELSE
            CALL UPGDST(NCART,MXRCRD,GRDARR,STPARR,EGRAD,STPSYM)
         END IF
C
C     We go through the same procedure as for the first iteration.
C
         IF (NEWTON) THEN
            CALL GTHESS(EGRAD,WORK(KEHESS),WORK(KALHES),
     &           EXHER,EXSIR,EXABA,WORK(KWRK1),LWRK1,WRKDLM)
         ELSE
            CALL GTGRAD(EGRAD,EXHER,EXSIR,EXABA,
     &           WORK(KWRK1),LWRK1,WRKDLM)
         END IF
C hj aug 99: reset hermit/abacus print to lower level for following iter.
         IPRUSR = IPRUSR_reduced
         IPREAD = IPREAD_reduced
C
C     If redundant internal coordinates are used, Wilson's B matrix,
C     its derivative and its inverse must be updated.
C
         IF (REDINT .OR. DELINT) THEN
            CALL GETWIL(MXRCRD,MX2CRD,WORK(KTEMP1),WILBMT,
     &           BMTRAN,WORK(KTEMP2))
#ifndef VAR_LINSCA
            IF (IPRINT .GE. IPRMAX)
     &           CALL GETDWL(MXRCRD,WORK(KTEMP1),WORK(KTEMP2),
     &           WORK(KTEMP3),WILBMT)
#endif
            CALL GTBINV(MXRCRD,WORK(KTEMP1),WORK(KTEMP2),WORK(KTEMP3),
     &           WORK(KTEMP4),WILBMT,BMTRAN,WORK(KBMINV),WORK(KPJINM),
     &           WORK(KTEMP5),WORK(KTEMP6),WORK(KWRK2),LWRK2)
         END IF
C
C     If new symmetry has been applied, coordinates has to counted
C     again.
C
         IF (NWSYMM) THEN
            CALL CNTCRD
            CALL MINCGH(EGRAD,WORK(KEHESS),WORK(KALHES),
     &           WORK(KWRK2),LWRK2)
            NWSYMM = .FALSE.
         END IF
         IF (REDINT .OR. DELINT) THEN
            NCRDHS = IINTCR
         ELSE
            NCRDHS = NCART
         END IF
         IF (RATFUN) NCRDHS = NCRDHS + 1
         IF (.NOT. NEWTON) THEN
            IF (REBILD) THEN
               CALL INIHES(MXRCRD,MX2CRD,SCLVEC,GRDOLD,HESOLD,
     &              WORK(KTEMP1),WORK(KTEMP2),WORK(KTEMP3),WORK(KTEMP4),
     &              WILBMT,BMTRAN,WORK(KBMINV),HESINT,WORK(KWRK2),LWRK2)
               REBILD = .FALSE.
            ELSE
               CALL UPDHES(MXRCRD,MX2CRD,SCLVEC,GRDOLD,GRDMAT,STPMAT,
     &              HESOLD,WORK(KTEMP1),WORK(KTEMP2),WORK(KTEMP3),
     &              WORK(KTEMP4),WORK(KTEMP5),WORK(KTEMP6),
     &              WORK(KTEMP7),WORK(KTEMP8),WORK(KTEMP9),
     &              WILBMT,BMTRAN,WORK(KBMINV),
     &              HESINT,NINT(ABS(GEINFO(ITRNMR-1,6))),
     &              NINT(ABS(GEINFO(ITRNMR,6))),WORK(KWRK2),LWRK2)
            END IF
            CALL MINCGH(EGRAD,WORK(KEHESS),WORK(KALHES),
     &           WORK(KWRK2),LWRK2)
         END IF
 33      CONTINUE
         IF (REDINT .OR. DELINT) THEN
            IF (NEWTON) CALL CGHINT(MXRCRD,MX2CRD,SCLVEC,WORK(KTEMP1),
     &           WORK(KTEMP2),WORK(KTEMP3),WORK(KTEMP4),WORK(KTEMP5),
     &           WILBMT,WORK(KBMINV),BMTRAN,HESINT,WORK(KWRK2),LWRK2)
            CALL PRJINT(MXRCRD,IINTCR,WORK(KPJINM),WORK(KCONMT),
     &           HESINT,WORK(KTEMP1),WORK(KTEMP2),WORK(KTEMP3),
     &           WORK(KTEMP4),WORK(KWRK2),LWRK2)
            IF (LNSRCH .AND. RATFUN .AND. (ITRNMR .GT. 0))
     &           CALL LINSRC(IINTCR,MXRCRD,GRDINT,GRDARR(1,1),
     &           WORK(KTEMP7),STPARR(1,1),WORK(KTEMP3),WORK(KTEMP4),
     &           ACTIVE,EMOD)
            IF (RATFUN .AND. SADDLE) NCRDHS = NCRDHS - 1
            CALL DIAINT(MXRCRD,MX2CRD,NCRDHS,WORK(KEVEC),WORK(KTEMP1),
     &           WORK(KTEMP2),WORK(KTEMP3),WORK(KTEMP4),THRIND,HESINT,
     &           WORK(KTEMP5),WORK(KWRK2),LWRK2)
            IF (RATFUN .AND. SADDLE) NCRDHS = NCRDHS + 1
         ELSE
            CALL PROJGH(EGRAD,WORK(KEHESS),WORK(KALHES),WORK(KTEMP1),
     &           WORK(KTEMP2),WORK(KTEMP3),WORK(KTEMP4))
            IF (LNSRCH .AND. RATFUN .AND. (ITRNMR .GT. 0))
     &           CALL LINSRC(NCART,MXCOOR,EGRAD,GRDARR(1,1),CSTEP,
     &           STPARR(1,1),WORK(KTEMP3),WORK(KTEMP4),ACTIVE,EMOD)
            CALL DIAHES(MXRCRD,MX2CRD,NCRDHS,EGRAD,WORK(KEHESS),
     &           WORK(KALHES),WORK(KTEMP1),THRIND,WORK(KEVEC),
     &           WORK(KTEMP2),WORK(KTEMP3),WORK(KTEMP4),
     &           WORK(KWRK2),LWRK2)
         END IF
C
C     Update information for this iteration
C
         GEINFO(ITRNMR,1) = ENERGY
         GEINFO(ITRNMR,3) = INDTOT*1.0D0
         IREJ = 0
C
C     Write Hessian to file
C
         IF (.NOT. NOHSWR)
     &      CALL PNCHES(MXRCRD,MX2CRD,HESINT,WILBMT,BMTRAN,WORK(KTEMP1),
     &      WORK(KTEMP2),WORK(KTEMP3),WORK(KTEMP4),WORK(KWRK2),LWRK2)
C
 756     CONTINUE
         IF (REDINT .OR. DELINT) THEN
            CALL FNSTIN(MXRCRD,MX2CRD,NCRDHS,SCLVEC,HESINT,WORK(KEVEC),
     &           WORK(KTEMP1),WORK(KTEMP2),WORK(KTEMP3),WORK(KTEMP4),
     &           WORK(KTEMP5),CSTEP,WILBMT,BMTRAN,WORK(KBMINV),GRDARR,
     &           STPARR,ACTIVE,EMOD,VECMOD,WORK(KTEMP7))
         ELSE
            CALL FNDSTP(MXRCRD,MX2CRD,NCRDHS,EGRAD,WORK(KEHESS),
     &           WORK(KEVEC),WORK(KTEMP1),WORK(KTEMP2),WORK(KTEMP3),
     &           WORK(KTEMP4),WORK(KTEMP5),CSTEP,GRDARR,STPARR,
     &           ACTIVE,EMOD,VECMOD)
         END IF
         GECONV = MINEND(MXRCRD,SCLVEC,BMTRAN,WORK(KTEMP1),WORK(KTEMP2))
         IF (ABS(GEINFO(ITRNMR,6)) .GT. 1.0D-3) GECONV = .FALSE.
         CALL PRIALL(CSTEP,WORK(KTEMP1))
C
C        save this geometry and energy to MOLDEN file
C
         IF (MOLDEN) CALL MOGECON(.FALSE.,ENERGY)
C
         IF (REDINT .AND. (IPRINT .GE. 1)) THEN
            CALL ATMINI(WORK(KTEMP1),IATOM,.TRUE.)
            CALL GETINT(IATOM,MXRCRD,WORK(KTEMP1),CRDINT)
            CALL HEADER('New internal coordinates',-1) 
            CALL OUTPUT(CRDINT,1,1,1,IINTCR,1,MXRCRD,1,LUPRI)
            WRITE(LUPRI,'(//)')
         END IF
         NEWSTP = .FALSE.
         IF ((.NOT. GECONV) .AND. ((.NOT. BRKSYM) .OR. NOBRKS))
     &        CALL FNDGEO(CSTEP,EGRAD,WORK(KTEMP1),WORK(KTEMP2),EXHER,
     &        EXSIR,EXABA,WORK(KWRK2),LWRK2,WRKDLM,IREJ,
     &        GEINFO,NEWSTP)
         IF (NEWSTP) GOTO 756
         GEINFO(ITRNMR,2) = GRADNM
         GEINFO(ITRNMR,4) = STPNRM
         IF (ITRNMR .LT. ITRMAX) GEINFO(ITRNMR+1,5) = TRSTRA
         IF (ABS(GEINFO(ITRNMR,6)) .LT. 1.0D-3) THEN
            GEINFO(ITRNMR,6) = IREJ*1.0D0
         ELSE
            GEINFO(ITRNMR,6) = -(ABS(GEINFO(ITRNMR,6))+ABS(IREJ)*1.0D0)
         END IF
         ITOTRJ = ITOTRJ + ABS(IREJ)
         IF (REBILD) THEN
            CALL DCOPY(IINTCR,STPINT,1,WORK(KTEMP7),1)
            CALL DZERO(STPINT,MXRCRD)
            CALL DZERO(GRDOLD,MXRCRD)
            DO 605 I = 1, IREDIC
               DO 607 J = 1, IINTCR
                  STPINT(I) = STPINT(I) + BMTRAN(I,J)*WORK(KTEMP7+J-1)
                  GRDOLD(I) = GRDOLD(I) + BMTRAN(I,J)*GRDINT(J)
 607           CONTINUE
 605        CONTINUE
         END IF
C
C     If the step has failed
C
         IF (IREJ .LT. 0) THEN
            IF (.NOT. NEWTON) THEN
               CALL INIHES(MXRCRD,MX2CRD,SCLVEC,GRDOLD,HESOLD,
     &              WORK(KTEMP1),WORK(KTEMP2),WORK(KTEMP3),WORK(KTEMP4),
     &              WILBMT,BMTRAN,WORK(KBMINV),HESINT,WORK(KWRK2),LWRK2)
               CALL MINCGH(EGRAD,WORK(KEHESS),WORK(KALHES),
     &              WORK(KWRK2),LWRK2)
            END IF
            GOTO 33
C         ELSE IF (REJINI .AND. REDINT .AND. (ITOTRJ .GE. 5)) THEN
C            WRITE(LUPRI,*)'***** NOTE! *****'
C            WRITE(LUPRI,*)
C     &           'The number of dihedral angles will be reduced!'
C            CALL RREDUN
C            ITOTRJ = 0
         END IF
C
C     Check if rejected steps should cause reinitialization of Hessian.
C
         IF ((.NOT. NEWTON) .AND. (REJINI .AND. (IREJ .GE. 1))) THEN
            WRITE(LUPRI,*)
            WRITE(LUPRI,*)'***** NOTE! *****'
            WRITE(LUPRI,*)
     &           'Due to rejected step, Hessian is reinitialized.'
            WRITE(LUPRI,*)
            CALL INIHES(MXRCRD,MX2CRD,SCLVEC,GRDOLD,HESOLD,
     &           WORK(KTEMP1),WORK(KTEMP2),WORK(KTEMP3),
     &           WORK(KTEMP4),WILBMT,BMTRAN,
     &           WORK(KBMINV),HESINT,WORK(KWRK2),LWRK2)
            TRSTRA = GEINFO(0,5)
            GEINFO(ITRNMR+1,5) = TRSTRA
            RSTARR = .TRUE.
         END IF
C
C     Check if increase of gradient norm should cause reinitialization
C     of Hessian. Reinitialization occurs when the norm of the gradient
C     is larger than the norm of the gradient two iterations earlier.
C
         IF (.NOT.NEWTON .AND. GRDINI .AND. (ITRNMR .GE. 2)) THEN
            IF (GEINFO(ITRNMR,2) .GE. GEINFO(ITRNMR-2,2)) THEN
               WRITE(LUPRI,*)
               WRITE(LUPRI,*)'***** NOTE! *****'
               WRITE(LUPRI,*)'Due to increasing gradient norm, ' //
     &              'Hessian is reinitialized.'
               WRITE(LUPRI,*)
               CALL INIHES(MXRCRD,MX2CRD,SCLVEC,GRDOLD,HESOLD,
     &              WORK(KTEMP1),WORK(KTEMP2),
     &              WORK(KTEMP3),WORK(KTEMP4),WILBMT,BMTRAN,
     &              WORK(KBMINV),HESINT,WORK(KWRK2),LWRK2)
               TRSTRA = GEINFO(0,5)
               GEINFO(ITRNMR+1,5) = TRSTRA
            END IF
         END IF
C
C     Check if symmetry should be broken
C
         IF (BRKSYM .AND. (.NOT. NOBRKS)) THEN
            IF (REDINT .OR. DELINT) THEN
               CALL DOBRKI(MXRCRD,MX2CRD,EXHER,EXSIR,EXABA,GEINFO,
     &              WORK(KEVEC),WORK(KBMINV),CSTEP,WORK(KTEMP1),
     &              WORK(KTEMP2),WORK(KTEMP3))
            ELSE
               CALL DOBRK(MXRCRD,MX2CRD,EXHER,EXSIR,EXABA,GEINFO,
     &              WORK(KEVEC),CSTEP,WORK(KTEMP1),WORK(KTEMP2))
            END IF
Chj aug 99: now new basis, reset hermit/abacus print
            IPRUSR = IPRUSR_orig
            IPREAD = IPREAD_orig
C
C     Check if preoptimization is finished
C
         ELSE IF (GECONV .AND. DOPRE .AND. (.NOT. FINPRE)) THEN
            CALL ENDPRE(EXHER,EXSIR,EXABA)
Chj aug 99: now new basis, reset hermit/abacus print
            IPRUSR = IPRUSR_orig
            IPREAD = IPREAD_orig
         END IF
C
C     Make VRML-file of next geometry if requested
C
         IF (VRML) CALL MKVRML(.TRUE.,WORK(KTEMP1),MXCOOR,WORK(KEVEC),
     &        WORK(KTEMP2),WORK(KTEMP3))
         GOTO 10
C
C     Finished case 1: Geometry has converged.
C
      ELSE IF (GECONV) THEN
C
C     We check to see if there is any properties to be calculated
C
         WORK(1) = WRKDLM
         WRINDX = .TRUE.
Chj aug 99: reset IPRUSR to original level
         IPRUSR = IPRUSR_orig
         IPREAD = IPREAD_orig
         CALL ABAINP('**PROPE',WORK(2),LWORK)
         CALL EXEABA(WORK,LWORK-1,WRKDLM)
C
C     We also check to see if we have requested a RESPONSE calculations
C
         CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ','FORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND (LUCMD)
C
C     Let us first see if we can find any **RESPONSE input
C
 133     CONTINUE
         READ (LUCMD,'(A7)',END=134,ERR=134) WRDRSP
         IF (WRDRSP .EQ. '**RESPO') THEN
            REWIND (LUCMD)
            RNRESP = .TRUE.
            CALL SETSIR(WORK(2),LWORK)
            CALL GPCLOSE(LUSIFC,'KEEP')
            CALL GPCLOSE(LUCMD,'KEEP')
            CALL RSPDRV(WORK(2),LWORK)
         ELSE
            GOTO 133
         END IF
 134     CONTINUE 
         IF (LUCMD .GT. 0) CALL GPCLOSE(LUCMD,'KEEP')
C
C     Single point energy is calculated if requested.
C
         WORK(1) = WRKDLM
         IF (DOSPE) CALL SPNRGY(GEINFO,EXHER,EXSIR,EXABA,WORK,
     &           LWORK,WRKDLM)
C
C     Final results are printed, partially through PRIINF.
C
         IF (MOLDEN) CALL MOGECON(.TRUE.,DUMMY)
         CALL TITLER(' End of Optimization ','<',120)
         CALL PRIINF(GEINFO,WORK(2),LWORK)
         WRITE(LUPRI,*)
         IF (CONOPT) THEN
            WRITE(LUPRI,*) 'Constrained optimization converged in ',
     &           ITRNMR+1, ' iterations!'
            IF (GRADNM .GT. GRDTHR) THEN
               WRITE(LUPRI,*) 'Removing the '
     &           // 'constraint(s) might decrease the energy further.'
            ELSE
               WRITE(LUPRI,*) 'A saddle point might have been reached.'
            END IF
         ELSE
            WRITE(LUPRI,*) 'Geometry converged in ', ITRNMR+1,
     &           ' iterations!'
         END IF
         IF (ITRBRK .GE. 0) THEN
            WRITE(LUPRI,*) 'Please note that symmetry '
     &           // 'was broken during this optimization.'
         ELSE IF (NOBRKS .AND. BRKSYM) THEN
            WRITE(LUPRI,*) 'Symmetry was not broken. '
     &           // 'Minimum was reached within the given symmetry.'
            WRITE(LUPRI,*) 'Please note that breaking this '
     &           // 'symmetry will decrease the energy further.'
         END IF
         IF (NEWTON .AND. SADDLE .AND. (INDTOT .NE. 1)) THEN
            WRITE(LUPRI,'(/A,I3,/A)')
     &         ' WARNING: Please note that the Hessian index',INDTOT,
     &         ' does not correspond to a first order saddle point ' //
     &           '(transition state).'
         END IF
         ENERGY = GEINFO(ITRNMR,1)
         WRITE(LUPRI,'(/A,F14.6,A)')
     &        ' Energy at final geometry is       : ',ENERGY,' a.u.'
         ERGDIF = ENERGY - GEINFO (0,1)
         WRITE(LUPRI,'(A,F14.6,A)')
     &        ' Energy change during optimization : ',ERGDIF,' a.u.'
         ERGDIF = ERGDIF * XKJMOL
         WRITE(LUPRI,'(A,F14.3,A)')
     &        '                                     ',ERGDIF,' kJ/mol'
         IF (DOPRE) THEN
            WRITE(LUPRI,'(/A)') ' Preoptimization was performed using'//
     &           ' the basis set(s):'
            DO 111 I = 1, IPRE-1
               WRITE(LUPRI,'(5X,A60)') PREBTX(I)
 111        CONTINUE
         END IF
         IF (DOSPE) THEN
            ENERGY = GEINFO(ITRNMR+1,1)
            WRITE(LUPRI,'(/A,A60)') ' Using the basis ',SPBSTX
            WRITE(LUPRI,'(A,F14.6,A)')
     &          ' single point energy was calculated: ',ENERGY,' a.u.'
         END IF
         WRITE(LUPRI,*)
C
C     Finished case 2: Exceeded maximum number of iterations.
C
      ELSE
C     No single point energy has been calculated.
         TMPLOG = DOSPE
         DOSPE = .FALSE.
         CALL TITLER('Optimization Control Center','<',120)
         CALL PRIINF(GEINFO,WORK(2),LWORK)
         DOSPE = TMPLOG
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'WARNING: Geometry has NOT converged!'
         WRITE(LUPRI,*) 'Maximum number of iterations (', ITRMAX,
     &                                        ') has been reached and'
         WRITE(LUPRI,*) 'optimization halted. Increase number or ',
     &                                   'restart from last geometry.'
         IF (DOSPE) WRITE(LUPRI,*)
     &                  'No single point energy has been calculated.'
         WRITE(LUPRI,*)
      END IF
      CALL QEXIT('RUNOPT')
      RETURN
      END

C  /* Deck priall */
      SUBROUTINE PRIALL(CSTEP,CORDNW)
C
C     Prints important information for the current iteration.
C
#include "implicit.h"
#include "mxcent.h"
#include "optinf.h"
#include "gnrinf.h"
#include "priunit.h"
#include "trkoor.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "symmet.h"
#include "nuclei.h"
#include "infinp.h"
#ifndef VAR_LINSCA
#include "pcmlog.h"
#endif
      DIMENSION CSTEP(MXCOOR), CORDNW(3,MXCENT)
      CHARACTER*1 MARK
      INTEGER     START
C
C     We have to construct the updated geometry before printing it
C
      IJ = 1
      DO 10 J = 1, NUCIND
         DO 20 I = 1, 3
            CORDNW(I,J) = CORD(I,J) + CSTEP(IJ)
            IJ = IJ + 1
 20      CONTINUE
 10   CONTINUE
C
      CALL TITLER('Optimization Control Center',':',115)
      IF (.NOT. GECONV) THEN
         CALL HEADER('Next geometry (au)',-1)
         CALL PRIGEO(CORDNW)
#ifndef VAR_LINSCA
         IF (PCM) CALL UPDCAV(CORDNW)
#endif
      ELSE
         CALL HEADER('Final geometry (au)',-1)
         CALL PRIGEO(CORD)
      END IF
      CALL HEADER('Optimization information',-1)
      WRITE(LUPRI,'(A,I7)')
     &     ' Iteration number               : ',ITRNMR
      MARK = ' '
      IF (BRKSYM .AND. GECONV) MARK = '*'
      IF (GECONV .AND. DOPRE .AND. (.NOT. FINPRE)) MARK = '*'
      WRITE(LUPRI,'(A,L1,A1)')
     &      ' End of optimization            :       ', GECONV, MARK
      WRITE(LUPRI,678) ' Energy at this geometry is     : ', ENERGY
      IF ((ITRNMR .GT. 0) .AND. (ITRBRK .LT. (ITRNMR - 1))) THEN
         ERGDIF = ENERGY - ERGOLD
         IF (ABS(ERGPRO) .GT. 1.0D-10) THEN
            RATIO = ERGDIF / ERGPRO
         ELSE
            RATIO = 1.0D0
         END IF
 678     FORMAT(A,F14.6)
         WRITE(LUPRI,678)' Energy change from last geom.  : ', ERGDIF
         IF (IPRINT .GE. 3) THEN
            WRITE(LUPRI,678)' Predicted change               : ',ERGPRO
            WRITE(LUPRI,678)' Ratio, actual/predicted change : ',RATIO
         END IF
      END IF
      WRITE(LUPRI,678) ' Norm of gradient               : ', GRADNM
      WRITE(LUPRI,678) ' Norm of step                   : ', STPNRM
      WRITE(LUPRI,678) ' Updated trust radius           : ', TRSTRA
c      if (itrnmr.eq.1) stop ' '
      IF (IPRINT .GE. 3 .OR. INDTOT.NE.INDHES(0)) THEN
         DO 50 I = 0, MAXREP
            IF (INDHES(I) .GT. 0) WRITE(LUPRI,'(A,I1,A,I8)')
     &           ' Hessian index (',I,')              :', INDHES(I)
 50      CONTINUE
      END IF
      WRITE(LUPRI,'(A,I8//)')' Total Hessian index            :',INDTOT
      IF (BRKSYM .AND. GECONV) THEN
         WRITE(LUPRI,'(/A)') ' *) Within given symmetry.'
      ELSE IF (GECONV .AND. DOPRE .AND. (.NOT. FINPRE)) THEN
         WRITE(LUPRI,'(/A)') ' *) End of preoptimization.'
      END IF
      START = 1
      NCRD = NCRTOT
      IF (REDINT .OR. DELINT) NCRD = IINTCR
      IF (RATFUN) THEN
         NCRD = NCRD + 1
         START = 2
      END IF
      IF (IPRINT .GE. 3) THEN
         CALL HEADER('Eigenvalues',-1)
         WRITE(LUPRI,*)
     &        '  #      Current value  Previous value      Change   '
         WRITE(LUPRI,*)
     &        '-----------------------------------------------------'
         DO 100 I = START, NCRD
            EVL   = EVAL(I)
            EVLOL = EVALOL(I)
            IF (EVL .GT. 9.9D3) THEN
               EVL   = 0.0D0
               EVLOL = 0.0D0
            END IF
            NR = I
            IF (RATFUN) NR = NR - 1
            WRITE(LUPRI,'(I4,3F16.6)') NR,EVL,EVLOL,EVL-EVLOL
 100     CONTINUE
      END IF
      RETURN
      END

C  /* Deck priinf */
      SUBROUTINE PRIINF(GEINFO,WORK,LWORK)
C
C     Prints important information at the end of an optimization.
C
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "optinf.h"
      DIMENSION GEINFO(0:ITRMAX,6), WORK(LWORK)
      CHARACTER NWSMRK*1,BRKMRK*1,REDMRK*1
      LOGICAL RED
      TMPNRG = 0.0D0
      RED = .FALSE.
      CALL HEADER('Final geometry',-1)
      CALL PRIGEO(CORD)
      WRITE(LUPRI,*)
     & 'Iter      Energy        Change       GradNorm  Index   StepLen'
     & // '    TrustRad #Rej'
      WRITE(LUPRI,*)
     & '-------------------------------------------------------------'
     & // '------------------'
      DO 10 I = 0, ITRNMR
         IF (I .GT. 0) TMPNRG = GEINFO(I,1) - GEINFO(I-1,1)
C
C     There's three special marks for each iteration, one (*) for a
C     Newton step (that is a step smaller than the trust radius),
C     one (x) for the breaking of symmetry and one (#) for dropping
C     half of dihedral angles.
C
         NWSMRK = ' '
         IF ((GEINFO(I,4)+1.0D-6) .LT. GEINFO(I,5)) NWSMRK = '*'
         BRKMRK = ' '
         IF (NINT(GEINFO(I,3)) .LT. 0) BRKMRK = 'x'
         REDMRK = ' '
         IF (GEINFO(I,6) .LT. -1.0D-3) THEN
            GEINFO(I,6) = ABS(GEINFO(I,6))
            REDMRK = '#'
            RED = .TRUE.
         END IF
         WRITE(LUPRI,'(I4,A1,F15.6,2F13.6,I5,F12.6,A1,F11.6,I4,A1)')
     &        I,BRKMRK,GEINFO(I,1),TMPNRG,GEINFO(I,2),
     &        NINT(ABS(GEINFO(I,3))),GEINFO(I,4),NWSMRK,GEINFO(I,5),
     &        NINT(GEINFO(I,6)),REDMRK
 10   CONTINUE
C
C     We also write the single point energy, if calculated.
C
      IF (DOSPE) WRITE(LUPRI,'(A5,F12.6,F11.6)') '  SP ',
     &          GEINFO(ITRNMR+1,1),GEINFO(ITRNMR+1,1)-GEINFO(ITRNMR,1)
      WRITE(LUPRI,*)
      WRITE(LUPRI,*) '*) Newton step taken.'
      IF (ITRBRK .GE. 0) WRITE(LUPRI,*)
     &     'x) Symmetry was broken after this iteration.'
      IF (RED) THEN
         IF (NEWTON) THEN
            WRITE(LUPRI,*)
     &     '#) Dropped half of the dihedral angles to ' //
     &     'reduce redundancy.'
         ELSE
            WRITE(LUPRI,*) '#) Hessian initialized to unity.'
         END IF
      END IF
      CALL GEOANA(CORD,.TRUE.,.FALSE.,NBONDS,-1,WORK,LWORK)
      RETURN
      END

C  /* Deck gtnrgy */
      SUBROUTINE GTNRGY(EXHER,EXSIR,EXABA,EXESG,WORK,LWORK,WRKDLM)
C
C     Retreives the energy of the current geometry
C     (by running HERMIT and SIRIUS).
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "infopt.h"
#include "infinp.h"
#include "optinf.h"
#include "gnrinf.h"
#include "siripc.h"
#include "dftcom.h"
#ifndef VAR_LINSCA
#include "pcmlog.h"
#endif
#include "molinp.h"
#include "abainf.h"
#include "rspprp.h"
#include "esg.h"
#ifdef VAR_LINSCA
#include "esg2.h"
#endif
Cef begin
#ifndef VAR_LINSCA
#include <incore.h>
#endif
Cef end
      LOGICAL EXHER, EXSIR, EXABA, EXESG, FNDKEY
      DIMENSION WORK(LWORK)
      REAL EXCQQ2
      CALL QENTER('GTNRGY')
      WORK(1) = WRKDLM
      IF (.NOT. EXHER) THEN
         CALL EXEHER(WORK,LWORK-1,WRKDLM)
         EXHER  = .TRUE.
         NEWGEO = .TRUE.
         NPCMIN = .TRUE.
      END IF
#ifndef VAR_LINSCA
      IF (AOSAVE .AND. NEWGEO) THEN
         LINTSV = .FALSE.
         LINTMP = .FALSE.
         INITX = .FALSE.
         MSAVE = .TRUE.
         MMCORE = MMWORK
         LMCORE = MMCORE
         ISCORE = 1
         JSCORE = ISCORE
         N_SHL = 1
         I_SHL = 1
         INDX_SHL1 = 0
         INDX_SHL2 = 0
         INDX_SHL3 = 0
         INDX_SHL4 = 0
C     CALL CLEAR_INCOREMEM()
      END IF
#endif
Cef end
C ach 
      IF (DFTRUN .OR. DOHFSRDFT .OR. DOMCSRDFT) GRDONE = .FALSE.
      IF (.NOT. EXSIR) THEN
C     
C        *** For numerical derivatives we might need to  ***
C        *** reassign the symmetry of the frozen core    ***
C        *** orbitals. tar 01-2001                       ***
C
         IF (NMWALK.AND.DOCCSD) THEN
            KWORD = 2
            KLAST = KWORD + 80*KMLINE
            IF (KLAST.GT.LWORK) CALL QUIT('Memory exceeded ' //
     &                                       'in fndexs.')
            CALL FNDEXS(WORK(KWORD),IPRINT)
         END IF
C
         CALL EXESIR(WORK,LWORK-1,WRKDLM)
C
C        Reduce gradient criterion for DFT optimization
C
         IF ((DFTRUN .OR. DOHFSRDFT .OR. DOMCSRDFT) .AND. .NOT.NATNRM
     &      .AND. (GRDTHR .LT. 5.0D-05 .OR. THRSTP .LT. 5.0D-05)) THEN
            GRDTHR = MAX(GRDTHR,5.0D-05)
            THRSTP = MAX(THRSTP,5.0D-05)
c           IF (CHGRDT) THEN
              NINFO = NINFO + 1
              WRITE (LUPRI,'(/A)') ' INFO>>>> Due to numerical '//
     &        'integration in DFT optimization'
              WRITE (LUPRI,'(A)') ' INFO>>>> thresholds for convergence'
     &        //' of geometry optimization have been reset'
              WRITE (LUPRI,'(/A/,3(/20X,A,F10.8))') ' New thresholds:',
     &        'Gradient norm  ',GRDTHR,'Step norm      ',THRSTP,
     &        'Energy change  ',THRERG
c           END IF
         END IF
         IF (DOCCSD) THEN
C
#ifdef VAR_LINSCA
#ifdef MOD_CC
            CALL CC_GEOFLAG
            CALL EXECC(WORK,LWORK-1,WRKDLM)
#else
            CALL QUIT('CC disabled in this build')
#endif
#else
            CALL CC_GEOFLAG
            CALL EXECC(WORK,LWORK-1,WRKDLM)
#endif
         ENDIF
         EXSIR = .TRUE.
      END IF
      IF (FLAG(16)) THEN
         WRITE (LUPRI,'(/A)') ' Solvent model not allowed with new'//
     &        'optimization module. Use .WALK'
         CALL QUIT('Solvent incompatible with .OPTIMIZE')
      END IF
C
C     ... sep09, kr+hjaaj: a test where we changed the SCF
C     convergence threshold THRGRD from 1.0d-6 to 1.0d-5,
C     the norm of the calculated molecular gradient changed with 1.0d-6.
      IF (.NOT.NATNRM .AND. .NOT. CHGRDT
     &    .AND. (1.99D0*THRGRD) .GT. GRDTHR
     &    .AND. (1.99D0*THRSTP) .GT. GRDTHR ) THEN
         WRITE (LUPRI,'(/A,1P,D10.2/A)') ' INFO>>>> Due to convergence '
     &        //'threshold for the wave function of',THRGRD,
     &        ' INFO>>>> the default thresholds for'
     &        //' convergence of geometry optimization have been reset.'
         WRITE (LUPRI,'(/A/,1P,3(/20X,A,D10.2))') ' Old thresholds:',
     &        'Gradient norm  ',GRDTHR,'Step norm      ',THRSTP,
     &        'Energy change  ',THRERG
         GRDTHR = MAX(GRDTHR,2.0D0*THRGRD)
         THRSTP = MAX(THRSTP,GRDTHR)
         THRERG = GRDTHR
         THRSYM = SQRT(THRERG)
         NINFO = NINFO + 1
         WRITE (LUPRI,'(/A/,1P,3(/20X,A,D10.2))') ' New thresholds:',
     &        'Gradient norm  ',GRDTHR,'Step norm      ',THRSTP,
     &        'Energy change  ',THRERG
      END IF

ckr#ifdef nohjaug99
C hjaaj aug99: no reason to update thresholds as it is
C an equil. optimization in the geometry although excited in MCSCF
C
C     Well, en excited electronic state cannot necessarily be expected to be 
C     as quadratic as the ground state, and thus we should not walk to far 
C     away. This is necessary for instance to get test15 through....
C
#ifdef VAR_LINSCA
      IF ((ISTATE .GT. 1 .OR.(FNDKEY('*LINSCA').AND. FNDKEY('.ESGRSP')))
     &      .AND. ITRNMR .EQ. 0) 
     &   THEN
#else
      IF (ISTATE .GT. 1 .AND. ITRNMR .EQ. 0) THEN
#endif
         WRITE (LUPRI,'(/A)') ' INFO >>>> As this is an excited '//
     &        'optimization some thresholds have been updated'
         IF (TRSTRA .EQ. 0.5D0) TRSTRA = 0.25D0
         IF (TRSTIN .EQ. 1.2D0) TRSTIN = 1.10D0
         IF (RTENBD .EQ. 0.4D0) RTENBD = 0.85D0
         IF (RTRJMN .EQ. 0.1D0) RTRJMN = 0.80D0
         IF (RTENGD .EQ. 0.8D0) RTENGD = 0.90D0
         IF (RTRJMX .EQ. 3.0D0) RTRJMX = 1.20D0
         WRITE(LUPRI,'(/A/A/,(A,F10.4))')
     &        ' Restricted step control parameters',
     &        ' ----------------------------------',
     &        ' Initial trust radius   :',TRSTRA,
     &        ' Trust radius increment :',TRSTIN,
     &        ' Trust radius decrement :',TRSTDE,
     &        ' Bad prediction ratio   :',RTENBD,
     &        ' Good prediction ratio  :',RTENGD,
     &        ' Rejection ratio, low   :',RTRJMN,
     &        ' Rejection ratio, high  :',RTRJMX
      END IF
ckr#endif
      IF ((DOMP2 .AND. .NOT. DOMC) .OR. DOCCSD .OR. ((DOCI .OR. DOCINO)
     &    .AND. .NOT. DOMC)) THEN
         ENERGY = ECORR
      ELSE
         ENERGY = EMCSCF
      END IF

C     Response excited state calculation
C     Maybe reset the step control parameters as above
Cdj moved from gtgrad
      IF (EXESG) THEN
         MOLGRD = .TRUE.
         CALL RSPDRV(WORK,LWORK)
         ENERGY = ENERGY + EXCITA(ISYME,IESG,1)
      END IF
#ifdef VAR_LINSCA
      IF (FNDKEY('*LINSCA').AND. FNDKEY('.ESGRSP')) THEN
         MOLGRD = .TRUE.
         CALL PRIOMA(EXCQQ2)
         ENERGY= ENERGY + EXCQQ2
      END IF
#endif
      CALL QEXIT('GTNRGY')
      RETURN
      END

C  /* Deck gtgrad */
      SUBROUTINE GTGRAD(EGRAD,EXHER,EXSIR,EXABA,WORK,LWORK,WRKDLM)
C
C     Retreives the energy and gradient of the current geometry
C     (by running HERMIT, SIRIUS and ABACUS).
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "taymol.h"
#include "maxaqn.h"
#include "symmet.h"
#include "nuclei.h"
#include "infinp.h"
#include "inftap.h"
#include "infpar.h"
#include "pario.h"
#include "abainf.h"
#include "gnrinf.h"
#include "optinf.h"
#include "pvibav.h"
#include "rspprp.h"
#include "esg.h"
#ifdef VAR_LINSCA
#include "esg2.h"
#endif
      LOGICAL EXHER, EXSIR, EXABA, EXESG, FNDKEY
      DIMENSION WORK(LWORK)
      DIMENSION EGRAD(MXCOOR)
      CALL QENTER('GTGRAD')
      WORK(1) = WRKDLM
      KFREE = 2
      LFREE = LWORK - KFREE + 1
      CALL DZERO(EGRAD,MXCOOR)
      
      EXESG = FNDKEY('*ESG   ')
#if defined (VAR_PARIO)
! WARNING! do not activate VAR_PARIO without fixing
!          the PARIO code !!!! /Feb 2011 hjaaj
      IF ((.NOT.PARIO) .OR. (MYNUM .EQ. 0)) THEN
#endif
         CALL GTNRGY(EXHER,EXSIR,EXABA,EXESG,WORK,LWORK,WRKDLM)
#if defined (VAR_PARIO)
      ELSE
chj-s-090516
chj      HRINPC = .FALSE.
chj      CALL GPOPEN(LUONEL,'AOONEINT',' ',' ',' ',IDUMMY,.FALSE.)
chj      CALL HERINP(WORK,LWORK)
chj      CALL GPCLOSE(LUONEL,'KEEP')
chj      HRINPC = .FALSE.
C        Third parameter .FALSE. in CALL READIN
C        indicates that LUONEL will NOT be written.
         CALL READIN(WORK(KFREE),LFREE,.FALSE.)
chj-e-090516
      END IF
#endif
C
C     Find out if we are running an excited state calculation
#ifdef VAR_LINSCA
      IF (.NOT.( FNDKEY('*ESG   ') .OR. 
     &    (FNDKEY('*LINSCA').AND. FNDKEY('.ESGRSP')) )) THEN
#else
      IF (.NOT. FNDKEY('*ESG   ')) THEN
#endif
C
C     TR: Using FNDANA to find highest order of analytical 
C         derivatives (NAORDR). If this is 0, we use numerical 
C         gradients.
C         NAORDR = 0; Only energies available
C         NAORDR = 1; Analytical gradients available
C         NAORDR = 2; Analytical gradients and hessians available.
C
         CALL FNDANA(NAORDR)
C
         IF (NAORDR.EQ.0) THEN 
            CALL NUMGRD(WORK,LWORK,WRKDLM)
         ELSE
            WRINDX = .TRUE.
            MOLGRD = .TRUE.
            IF (HFPROP) HELFEY = .TRUE.
            IF (ITRNMR .EQ. 0) THEN
               CALL ABAINP('**START',WORK(KFREE),LFREE)
            ELSE
#ifdef VAR_LINSCA
               CALL ABAINP('**PROPE',WORK(KFREE),LWORK)
#else
               CALL ABAINP('**EACH ',WORK(KFREE),LWORK)
#endif
            END IF
C        *** For numerical differentiation (vibrational ***
C        *** average), we don't nescesarily calculate   ***
C        *** all properties.                            ***
C
            IF (.NOT.CNMPRP) THEN
               CALL NVBPIN(IPRINT)
            END IF
C     
            IF ( DOCCSD ) THEN
#ifdef VAR_LINSCA
#ifdef MOD_CC
               CALL CC_GEOFLAG
               CALL EXECC(WORK,LWORK-1,WRKDLM)
#else
               CALL QUIT('CC disabled in this build')
#endif
#else
               CALL CC_GEOFLAG
               CALL EXECC(WORK,LWORK-1,WRKDLM)
#endif
            ELSE
               CALL EXEABA(WORK,LWORK-1,WRKDLM)
            END IF
            EXABA = .TRUE.
            MOLGRD = .FALSE.
            HELFEY = .FALSE.
         END IF
      END IF
      CALL TRAGRD(GRDMOL,EGRAD,WORK(KFREE),
     &            WORK(KFREE + MXCOOR*MXCOOR),NCRREP(0,1),3*NUCDEP)
#ifdef VAR_LINSCA
      IF (FNDKEY('*LINSCA').AND. FNDKEY('.ESGRSP')) THEN
        MOLGRD = .TRUE.
      END IF
#endif
      CALL QEXIT('GTGRAD')
      RETURN
      END

C  /* Deck gthess */
      SUBROUTINE GTHESS(EGRAD,EHESS,ALLHES,EXHER,EXSIR,EXABA,
     &                                 WORK,LWORK,WRKDLM)
C
C     Retreives the energy, gradient and hessian of the current geometry
C     (by running HERMIT, SIRIUS and ABACUS).
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "taymol.h"
#include "maxaqn.h"
#include "symmet.h"
#include "nuclei.h"
#include "abainf.h"
#include "gnrinf.h"
#include "optinf.h"
#include "trkoor.h"
#include "cbiwlk.h"
#include "pvibav.h"
      LOGICAL EXHER, EXSIR, EXABA
      DIMENSION WORK(LWORK)
      DIMENSION EGRAD(MXCOOR), EHESS(MXCOOR,MXCOOR)
      DIMENSION ALLHES(NCRTOT*NCRTOT)
      CALL QENTER('GTHESS')
      WORK(1) = WRKDLM
      KFREE = 2
      LFREE = LWORK - KFREE + 1
      CALL DZERO(EGRAD,MXCOOR)
      CALL DZERO(EHESS,MXCOOR*MXCOOR)
      CALL GTNRGY(EXHER,EXSIR,EXABA,.FALSE.,WORK,LWORK,WRKDLM)
#ifdef VAR_LINSCA
      CLOSE (13)
#endif
      WRINDX = .TRUE.
      MOLGRD = .TRUE.
      MOLHES = .TRUE.
      IF (HFPROP) HELFEY = .TRUE.
      IF (ITRNMR .EQ. 0) THEN
         CALL ABAINP('**START',WORK(KFREE),LFREE)
      ELSE
#ifdef VAR_LINSCA
         CALL ABAINP('**PROPE',WORK(KFREE),LWORK)
#else
         CALL ABAINP('**EACH ',WORK(KFREE),LFREE)
#endif
      END IF
C
C     *** For numerical differentiation (vibrational ***
C     *** average), we don't nescesary calculate all ***
C     *** properties.                                ***
C
      IF (.NOT.CNMPRP) THEN
         CALL NVBPIN(IPRINT)
      END IF
C
      CALL EXEABA(WORK,LWORK-1,WRKDLM)
      EXABA = .TRUE.
      START = .FALSE.
      MOLGRD = .FALSE.
      MOLHES = .FALSE.
      HELFEY = .FALSE.
C
      IF (NMWALK) THEN
C
C        *** For numerical derivatives, the Hessian is ***
C        *** "desymmetrized".                          ***
C
         KCSTRA = KFREE
         KSCTRA = KCSTRA + NCOOR**2
         CALL TRAHES(HESMOL,MXCOOR,EHESS,WORK(KCSTRA),WORK(KSCTRA),
     &               MXCOOR,NCOOR,1)
      ELSE IF ((ITRNMR .GT. 0) .AND. (.NOT. NWSYMM)) THEN
C
C        For all iterations except the first, MINCGH is run here.
C
         CALL MINCGH(EGRAD,EHESS,ALLHES,WORK(KFREE),LFREE)
      END IF
      CALL QEXIT('GTHESS')
      RETURN
      END

C  /* Deck inihes */
      SUBROUTINE INIHES(MXRCRD,MX2CRD,SCLVEC,GRDOLD,HESOLD,
     &     TMPMAT,TMPMT2,TMPMT3,TMPMT4,WILBMT,BMTRAN,BMTINV,
     &     HESINT,WORK,LWORK)
C
C     Initializes Hessian
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "optinf.h"
#include "taymol.h"
#include "nuclei.h"
#include "symmet.h"
      DIMENSION SCLVEC(MXCOOR), GRDOLD(MXRCRD)
      DIMENSION HESOLD(MXRCRD,MXRCRD), TMPMAT(MX2CRD,MX2CRD)
      DIMENSION TMPMT2(MX2CRD,MX2CRD), TMPMT3(MX2CRD,MX2CRD)
      DIMENSION TMPMT4(MX2CRD,MX2CRD), WILBMT(MX2CRD,MX2CRD)
      DIMENSION BMTRAN(MXRCRD,MXRCRD)
      DIMENSION BMTINV(MXRCRD,MXCOOR), HESINT(MXRCRD,MXRCRD)
      DIMENSION WORK(LWORK)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      CALL DZERO(HESOLD,MXRCRD*MXRCRD)
      CALL DZERO(GRDOLD,MXRCRD)
      CALL DZERO(SCLVEC,MXCOOR)
C
C     We initialize the scaling vector (the inverse of the
C     normalization vector in WLKCGH).
C
      DO 1 IREP = 0, MAXREP
         DO 2 ICENT = 1, NUCIND
#include "scalar.h"
            DO 3 ICOOR = 1, 3
               ISCOOR = IPTCNT(3*(ICENT-1)+ICOOR,IREP,1)
               IF (ISCOOR .GT. 0) THEN
                  SCLVEC(ISCOOR) = SQRT(FMULT(ISTBNU(ICENT)))
               END IF
 3          CONTINUE
 2       CONTINUE
 1    CONTINUE
C
      IF (IPRINT .GE. IPRMED) THEN
         CALL TITLER('Output from INIHES','*',103)
         IF (IPRINT .GE. IPRMAX) THEN
            CALL HEADER('Scaling vector',-1)
            CALL OUTPUT(SCLVEC,1,1,1,NCRTOT,1,MXCOOR,1,LUPRI)
         END IF
      END IF
C
C     Check if Hessian should be read from file
C
      IF (HSFILE) THEN
         IF (REDINT .OR. DELINT) THEN
            IF (MAXREP .GT. 0) THEN
               CALL DZERO(TMPMAT,MX2CRD*MX2CRD)
               CALL TRAGRD(GRDMOL,TMPMAT,TMPMT2,TMPMT3,
     &              NCRREP(0,1),3*NUCDEP)
               CALL GX2GQ(MXRCRD,TMPMAT,GRDINT,BMTINV)
               GRADNM = SQRT(DDOT(ICRTCR,TMPMAT,1,TMPMAT,1))
            ELSE
               CALL GX2GQ(MXRCRD,GRDMOL,GRDINT,BMTINV)
               GRADNM = SQRT(DDOT(ICRTCR,GRDMOL,1,GRDMOL,1))
            END IF
         END IF
         CALL REAHES(MXRCRD,MX2CRD,HESINT,HESOLD,TMPMAT,TMPMT2,TMPMT3,
     &        TMPMT4,WILBMT,BMTRAN,BMTINV,WORK,LWORK,IERR)
         CALL DZERO(HESOLD,MXRCRD*MXRCRD)
         IF (IERR .EQ. -1) THEN
            CALL QUIT('Unable to open the file DALTON.HES.')
         ELSE IF (IERR .EQ. -2) THEN
            CALL QUIT('The Hessian in DALTON.HES has ' //
     &                                     'wrong dimensions.')
         END IF
         IF (REDINT .OR. DELINT) THEN
            DO 5 J = 1, IINTCR
               GRDOLD(J) = GRDINT(J)
               DO 7 I = 1, IINTCR
                  HESOLD(I,J) = HESINT(I,J)
 7             CONTINUE
 5          CONTINUE
         ELSE
            DO 12 J = 1, NCRTOT
               GRDOLD(J) = GRDMOL(J)
               DO 14 I = 1, NCRTOT
                  HESOLD(I,J) = HESMOL(I,J)*(SCLVEC(I)*SCLVEC(J))
                  HESMOL(I,J) = HESOLD(I,J)*(SCLVEC(I)*SCLVEC(J))
 14            CONTINUE
 12         CONTINUE
         END IF
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Initial Hessian has been read from file.'
         WRITE(LUPRI,*)
         HSFILE = .FALSE.
C
C     Check if Hessian has been calculated
C
      ELSE IF (INITHS) THEN
         IF (REDINT .OR. DELINT) THEN
            IF (MAXREP .GT. 0) THEN
               CALL DZERO(TMPMAT,MX2CRD*MX2CRD)
               CALL TRAGRD(GRDMOL,TMPMAT,TMPMT2,TMPMT3,
     &                     NCRREP(0,1),3*NUCDEP)
               CALL GX2GQ(MXRCRD,TMPMAT,GRDINT,BMTINV)
               GRADNM = SQRT(DDOT(ICRTCR,TMPMAT,1,TMPMAT,1))
               IF (IPRINT .GE. IPRDBG) THEN
                  CALL HEADER('Cartesian gradient',-1)
                  CALL OUTPUT(TMPMAT,1,1,1,ICRTCR,1,MXCOOR,1,LUPRI)
                  CALL HEADER('Internal gradient',-1)
                  CALL OUTPUT(GRDINT,1,1,1,IINTCR,1,MXRCRD,1,LUPRI)
               END IF
               CALL DZERO(HESOLD,MXRCRD*MXRCRD)
               CALL TRAHES(HESMOL,MXCOOR,HESOLD,TMPMT2,TMPMT3,
     &                     MXCOOR,3*NUCDEP,1)
               CALL HX2HQ(MXRCRD,MX2CRD,TMPMAT,TMPMT2,TMPMT3,TMPMT4,
     &              HESOLD,GRDINT,HESINT,WILBMT,BMTINV,BMTRAN,
     &              WORK,LWORK)
               IF (IPRINT .GE. IPRDBG) THEN
                  CALL HEADER('Cartesian Hessian',-1)
                  CALL OUTPUT(HESOLD,1,ICRTCR,1,ICRTCR,
     &                                      MXCOOR,MXCOOR,1,LUPRI)
               END IF
               CALL DZERO(HESOLD,MXRCRD*MXRCRD)
            ELSE
               CALL GX2GQ(MXRCRD,GRDMOL,GRDINT,BMTINV)
               GRADNM = SQRT(DDOT(ICRTCR,GRDMOL,1,GRDMOL,1))
               IF (IPRINT .GE. IPRDBG) THEN
                  CALL HEADER('Cartesian gradient',-1)
                  CALL OUTPUT(GRDMOL,1,1,1,ICRTCR,1,MXCOOR,1,LUPRI)
                  CALL HEADER('Internal gradient',-1)
                  CALL OUTPUT(GRDINT,1,1,1,IINTCR,1,MXRCRD,1,LUPRI)
               END IF
               CALL HX2HQ(MXRCRD,MX2CRD,TMPMAT,TMPMT2,TMPMT3,TMPMT4,
     &              HESMOL,GRDINT,HESINT,WILBMT,BMTINV,BMTRAN,
     &              WORK,LWORK)
               IF (IPRINT .GE. IPRDBG) THEN
                  CALL HEADER('Cartesian Hessian',-1)
                  CALL OUTPUT(HESMOL,1,ICRTCR,1,ICRTCR,
     &                                      MXCOOR,MXCOOR,1,LUPRI)
               END IF
            END IF
            DO 15 J = 1, IINTCR
               DO 16 I = 1, IINTCR
                  HESOLD(I,J) = HESINT(I,J)
 16            CONTINUE
               GRDOLD(J) = GRDINT(J)
 15         CONTINUE
         ELSE
            DO 18 J = 1, NCRTOT
               DO 19 I = 1, NCRTOT
                  HESOLD(I,J) = HESMOL(I,J)/(SCLVEC(I)*SCLVEC(J))
 19            CONTINUE
               GRDOLD(J) = GRDMOL(J)
 18         CONTINUE
         END IF
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Initial Hessian has been calculated.'
         WRITE(LUPRI,*)
C
C     Next line makes reinitialization possible:
C
         INITHS = .FALSE.
C
C     Check if initial Hessian should be diagonal in redundant internal
C     coordinates.
C
      ELSE IF (REDINT .OR. DELINT .OR. INMDHS .OR. INRDHS) THEN
         NRIC = IINTCR
         IF (DELINT) NRIC = IREDIC
         IF (MAXREP .GT. 0) THEN
            CALL DZERO(TMPMAT,MX2CRD*MX2CRD)
            CALL TRAGRD(GRDMOL,TMPMAT,TMPMT2,TMPMT3,
     &                  NCRREP(0,1),3*NUCDEP)
            CALL GX2GQ(MXRCRD,TMPMAT,GRDINT,BMTINV)
            GRADNM = SQRT(DDOT(ICRTCR,TMPMAT,1,TMPMAT,1))
         ELSE
            CALL GX2GQ(MXRCRD,GRDMOL,GRDINT,BMTINV)
            GRADNM = SQRT(DDOT(ICRTCR,GRDMOL,1,GRDMOL,1))
         END IF
C
         IF (IPRINT .GE. IPRDBG) THEN
            CALL HEADER('Cartesian gradient',-1)
            CALL OUTPUT(GRDMOL,1,1,1,NCART,1,MXCOOR,1,LUPRI)
         END IF
         IF (IPRINT .GE. IPRMAX) THEN
            CALL HEADER('Internal gradient',-1)
            CALL OUTPUT(GRDINT,1,1,1,IINTCR,1,MXRCRD,1,LUPRI)
         END IF
C
C     The default is different values for bonds and other internal
C     coordinates.
C
         CALL DZERO(HESINT,MXRCRD*MXRCRD)
         IF (EVLINI .LE. 0.0D0) THEN
            IF (MODHES .OR. INMDHS .OR. CMBMOD) THEN
               CALL BLDHES(MXRCRD,TMPMAT,HESINT)
            ELSE
               DO 20 I = 1, NRIC
                  IF (INTCRD(I,1) .LT. 10) THEN
                     HESINT(I,I) = 0.5D0
                  ELSE IF (INTCRD(I,1) .LT. 20) THEN
                     HESINT(I,I) = 0.2D0
                  ELSE
                     HESINT(I,I) = 0.1D0
                  END IF
 20            CONTINUE
            END IF
         ELSE
            DO 22 I = 1, IINTCR
               HESINT(I,I) = EVLINI
 22         CONTINUE
         END IF
C
C     If we use delocalized internals, we have to transform from
C     redundant to non-redundant space.
C
         IF (DELINT .AND. (EVLINI .LE. 0.0D0)) THEN
            CALL DZERO(TMPMAT,MX2CRD*MX2CRD)
            DO 200 I = 1, IINTCR
               DO 202 J = 1, NRIC
                  DO 204 K = 1, NRIC
                     TMPMAT(I,J) = TMPMAT(I,J) + BMTRAN(K,I)*HESINT(K,J)
 204              CONTINUE
 202           CONTINUE
 200        CONTINUE
            CALL DZERO(HESINT,MXRCRD*MXRCRD)
            DO 210 I = 1, IINTCR
               DO 212 J = 1, IINTCR
                  DO 214 K = 1, NRIC
                     HESINT(I,J) = HESINT(I,J) + TMPMAT(I,K)*BMTRAN(K,J)
 214              CONTINUE
 212           CONTINUE
 210        CONTINUE
         END IF
         IF (IPRINT .GE. IPRMAX) THEN
            CALL HEADER('Internal Hessian',-1)
            CALL OUTPUT(HESINT,1,IINTCR,1,IINTCR,
     &                                       MXRCRD,MXRCRD,1,LUPRI)
         END IF
         IF ((INRDHS .OR. INMDHS) .AND. ((.NOT. REDINT)
     &        .AND. (.NOT. DELINT))) THEN
            CALL DZERO(HESMOL,MXCOOR*MXCOOR)
            CALL HQ2HX(MXRCRD,MX2CRD,TMPMAT,TMPMT2,TMPMT3,HESINT,
     &           GRDINT,HESMOL,WILBMT,BMTRAN,WORK,LWORK)
            IF (MAXOPR .GT. 0) THEN
               CALL DZERO(TMPMT4,MX2CRD*MX2CRD)
               DO 100 J = 1, ICRTCR
                  DO 102 I = 1, ICRTCR
                     TMPMT4(I,J) = HESMOL(I,J)
 102              CONTINUE
 100           CONTINUE
               IF (IPRINT .GE. IPRDBG) THEN
                  CALL HEADER('Expanded Cartesian Hessian',-1)
                  CALL OUTPUT(HESMOL,1,ICRTCR,1,ICRTCR,MXCOOR,MXCOOR,
     &                                                       1,LUPRI)
               END IF
               DO 715 J = 1, ICRTCR
                  DO 717 I = 1, J-1
                     HESMOL(J,I) = HESMOL(I,J)
 717              CONTINUE
 715           CONTINUE
C
C     Collect the symmetry-to-cartesian transformation matrix in TMPMT3
C
               CALL TRACOR(TMPMT2,TMPMT3,1,ICRTCR,IPRINT)
C
               CALL DGEMM('T','N',ICRTCR,ICRTCR,ICRTCR,1.D0,
     &                    TMPMT3,ICRTCR,
     &                    HESMOL,MXCOOR,0.D0,
     &                    TMPMAT,MXCOOR)
               CALL DZERO(HESMOL,MXCOOR*MXCOOR)
               CALL DGEMM('N','N',ICRTCR,ICRTCR,ICRTCR,1.D0,
     &                    TMPMAT,MXCOOR,
     &                    TMPMT3,ICRTCR,0.D0,
     &                    HESMOL,MXCOOR)
               IF (IPRINT .GE. IPRDBG) THEN
                  CALL HEADER('Totally Symmetric Cartesian Hessian',-1)
                  CALL OUTPUT(HESMOL,1,NCART,1,NCART,MXCOOR,MXCOOR,
     &                                                         1,LUPRI)
               END IF
            END IF
            DO 25 J = 1, ICRTCR
               GRDOLD(J) = GRDMOL(J)
               DO 27 I = 1, ICRTCR
                  HESOLD(I,J) = HESMOL(I,J)*SCLVEC(I)*SCLVEC(J)
                  HESMOL(I,J) = HESMOL(I,J)*SCLVEC(I)*SCLVEC(I)
     &                                             *SCLVEC(J)*SCLVEC(J)
 27            CONTINUE
 25         CONTINUE
         ELSE
            DO 30 J = 1, IINTCR
               GRDOLD(J) = GRDINT(J)
               DO 32 I = 1, IINTCR
                  HESOLD(I,J) = HESINT(I,J)
 32            CONTINUE
 30         CONTINUE
         END IF
C         INMDHS = .FALSE.
C
C     Otherwise the Hessian is set equal to a diagonal matrix. When
C     symmetry is used, the Hessian must be scaled so that it will be
C     correct after normalization.
C
      ELSE
         CALL DZERO(HESMOL,MXCOOR*MXCOOR)
         DO 40 I = 1, NCRTOT
            HESMOL(I,I) = EVLINI*SCLVEC(I)*SCLVEC(I)
            HESOLD(I,I) = EVLINI
            GRDOLD(I) = GRDMOL(I)
 40      CONTINUE
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Initial Hessian is a diagonal matrix.'
         WRITE(LUPRI,*)
      END IF
C
      IF (IPRINT .GE. IPRMED) THEN
         CALL HEADER('Initial Hessian',-1)
         IF (REDINT .OR. DELINT) THEN
            CALL OUTPUT(HESOLD,1,IINTCR,1,IINTCR,MXRCRD,MXRCRD,
     &                                                       1,LUPRI)
         ELSE
            CALL OUTPUT(HESOLD,1,NCART,1,NCART,MXRCRD,MXRCRD,
     &                                                       1,LUPRI)
         END IF
      END IF
      RETURN
      END

C  /* Deck updhes */
      SUBROUTINE UPDHES(MXRCRD,MX2CRD,SCLVEC,GRDOLD,GRDMAT,STPMAT,
     &     HESOLD,GAMMA,TMPMAT,TMPMT2,TMPMT3,TMPMT4,TMPMT5,TMPMT6,
     &     TMPMT7,TMPMT8,WILBMT,BMTRAN,BMTINV,HESINT,IREJOL,IREJNW,
     &     WORK,LWORK)
C
C     Controls Hessian updates
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "optinf.h"
#include "symmet.h"
#include "nuclei.h"
#include "taymol.h"
      DIMENSION SCLVEC(MXCOOR), GRDOLD(MXRCRD)
      DIMENSION GRDMAT(25,MXRCRD), STPMAT(25,MXRCRD)
      DIMENSION HESOLD(MXRCRD,MXRCRD), GAMMA(MXRCRD)
      DIMENSION TMPMAT(MX2CRD*MX2CRD), TMPMT2(MX2CRD*MX2CRD)
      DIMENSION TMPMT3(MX2CRD*MX2CRD), TMPMT4(MX2CRD*MX2CRD)
      DIMENSION TMPMT5(MX2CRD*MX2CRD), TMPMT6(MX2CRD*MX2CRD)
      DIMENSION TMPMT7(MX2CRD*MX2CRD), TMPMT8(MX2CRD*MX2CRD)
      DIMENSION WILBMT(MXRCRD,MXCOOR), BMTRAN(MXRCRD,MXRCRD)
      DIMENSION BMTINV(MXRCRD,MXCOOR), HESINT(MXRCRD,MXRCRD)
      DIMENSION WORK(LWORK)
      LOGICAL RESET, REJLST
C
      CALL QENTER('UPDHES')
C
C     If redundant internal coordinates are used, we have to transform
C     the calculated Cartesian gradient
C
      IF (REDINT .OR. DELINT) THEN
         IF (MAXREP .GT. 0) THEN
            CALL DZERO(TMPMAT,MX2CRD*MX2CRD)
            CALL TRAGRD(GRDMOL,TMPMAT,TMPMT2,TMPMT3,
     &                  NCRREP(0,1),3*NUCDEP)
            CALL GX2GQ(MXRCRD,TMPMAT,GRDINT,BMTINV)
            GRADNM = SQRT(DDOT(ICRTCR,TMPMAT,1,TMPMAT,1))
         ELSE
            CALL GX2GQ(MXRCRD,GRDMOL,GRDINT,BMTINV)
            GRADNM = SQRT(DDOT(ICRTCR,GRDMOL,1,GRDMOL,1))
         END IF
      END IF
C
C     When a new basis set is used in preoptimization, the Hessian is
C     kept as it is, because no step is taken.
C
      IF (KEEPHE) THEN
         CALL UPDOLD(MXRCRD,SCLVEC,GRDOLD,HESOLD,HESINT)
C
C     The requested Hessian update method is used
C
      ELSE
         IF (STEEPD) THEN
            CALL UPOUTP(MXRCRD,SCLVEC,GRDOLD,GAMMA,HESOLD,'UPSTPD')
            CALL INIHES(MXRCRD,MX2CRD,SCLVEC,GRDOLD,HESOLD,TMPMAT,
     &           TMPMT2,TMPMT3,TMPMT4,WILBMT,BMTRAN,BMTINV,HESINT,
     &           WORK,LWORK)
            IF ((.NOT. REDINT) .AND. (.NOT. DELINT)) THEN
               DO 10 J = 1, NCART
                  DO 12 I = 1, NCART
                     HESMOL(I,J) = HESMOL(I,J)/(SCLVEC(I)*SCLVEC(J))
 12               CONTINUE
 10            CONTINUE
            END IF
         ELSE IF (MODHES) THEN
            CALL UPOUTP(MXRCRD,SCLVEC,GRDOLD,GAMMA,HESOLD,'UPMODH')
            ITYPE = 1
            IF (DFP) ITYPE = 2
            IF (PSB) ITYPE = 3
            IF (RANKON) ITYPE = 4
            RESET = .FALSE.
            REJLST = .FALSE.
            IF (ITRBRK .EQ. (ITRNMR-2)) RESET = .TRUE.
C            IF (REJINI .AND. (IREJOL .GT. 0)) RESET = .TRUE.
C            IF (IREJNW .GT. MAXREJ) RESET = .TRUE.
C            IF (IREJOL+IREJNW .GT. 0) REJLST = .TRUE.
            CALL UPMODH(MXRCRD,STPINT,STPMAT,GAMMA,GRDMAT,HESOLD,HESINT,
     &           IINTCR,IREDIC,MXRCRD,TMPMAT(1),TMPMAT(MX2CRD+1),
     &           TMPMT2,TMPMT3,TMPMT4,TMPMT5,TMPMT6,TMPMT7,TMPMT8,
     &           BMTRAN,IINTCR-NPROJ,RESET,REJLST,GRADNM,
     &           ITYPE,DELINT,IPRINT)
         ELSE IF (CMBMOD) THEN
            CALL UPOUTP(MXRCRD,SCLVEC,GRDOLD,GAMMA,HESOLD,'UPCMBM')
            ITYPE = 1
            IF (DFP) ITYPE = 2
            IF (PSB) ITYPE = 3
            IF (RANKON) ITYPE = 4
            RESET = .FALSE.
            REJLST = .FALSE.
            IF (ITRBRK .EQ. (ITRNMR-2)) RESET = .TRUE.
C            IF (REJINI .AND. (IREJOL .GT. 0)) RESET = .TRUE.
            IF (IREJNW .GT. MAXREJ) RESET = .TRUE.
            IF (IREJOL+IREJNW .GT. 0) REJLST = .TRUE.
            CALL UPCMBM(MXRCRD,STPINT,STPMAT,GAMMA,GRDMAT,HESOLD,HESINT,
     &           IINTCR,IREDIC,MXRCRD,TMPMAT(1),TMPMAT(MX2CRD+1),
     &           TMPMT2,TMPMT3,TMPMT4,TMPMT5,TMPMT6,TMPMT7,TMPMT8,
     &           BMTRAN,IINTCR-NPROJ,RESET,REJLST,GRADNM,
     &           ITYPE,DELINT,IPRINT)
         ELSE IF (MULTI) THEN
            CALL UPOUTP(MXRCRD,SCLVEC,GRDOLD,GAMMA,HESOLD,'UPMULT')
            RESET = .FALSE.
            REJLST = .FALSE.
            IF ((ITRNMR .LT. 2) .OR. (ITRBRK .EQ. (ITRNMR-2)))
     &           RESET = .TRUE.
            IF (REJINI .AND. (ABS(IREJOL) .GT. 0)) RESET = .TRUE.
            IF (REJINI .AND. (IREJOL+IREJNW .GT. 0)) RESET = .TRUE.
            ITYPE = 1
            IF (PSB) ITYPE = 2
            IF (DFP) ITYPE = 3
            IF (RANKON) ITYPE = 4
            IF (REDINT .OR. DELINT) THEN
               CALL UPMULT(MXRCRD,STPINT,STPMAT,GAMMA,GRDMAT,HESOLD,
     &              HESINT,IINTCR,MXRCRD,TMPMAT(1),TMPMAT(MX2CRD+1),
     &              TMPMT2,TMPMT3,TMPMT4,TMPMT5,TMPMT6,IINTCR-NPROJ,
     &              RESET,GRADNM,IPRINT,ITYPE,.TRUE.)
            ELSE
               CALL UPMULT(MXRCRD,STPSYM,STPMAT,GAMMA,GRDMAT,HESOLD,
     &              HESMOL,NCART,MXCOOR,TMPMAT(1),TMPMAT(MX2CRD+1),
     &              TMPMT2,TMPMT3,TMPMT4,TMPMT5,TMPMT6,NCART-NPROJ,
     &              RESET,GRADNM,IPRINT,ITYPE,.TRUE.)
            END IF
         ELSE IF (BFGS) THEN 
            CALL UPOUTP(MXRCRD,SCLVEC,GRDOLD,GAMMA,HESOLD,'UPBFGS')
            IF (REDINT .OR. DELINT) THEN
               CALL UPBFGS(MXRCRD,STPINT,GAMMA,HESOLD,HESINT,
     &              IINTCR,MXRCRD,TMPMAT,TMPMT2,TMPMT3,IPRINT)
            ELSE
               CALL UPBFGS(MXRCRD,STPSYM,GAMMA,HESOLD,HESMOL,
     &              NCART,MXCOOR,TMPMAT,TMPMT2,TMPMT3,IPRINT)
            END IF
         ELSE IF (DFP) THEN 
            CALL UPOUTP(MXRCRD,SCLVEC,GRDOLD,GAMMA,HESOLD,'UPDFP ')
            IF (REDINT .OR. DELINT) THEN
               CALL UPDFP(MXRCRD,STPINT,GAMMA,HESOLD,HESINT,
     &              IINTCR,MXRCRD,TMPMAT,TMPMT2,TMPMT3,IPRINT)
            ELSE
               CALL UPDFP(MXRCRD,STPSYM,GAMMA,HESOLD,HESMOL,
     &              NCART,MXCOOR,TMPMAT,TMPMT2,TMPMT3,IPRINT)
            END IF
         ELSE IF (PSB) THEN 
            CALL UPOUTP(MXRCRD,SCLVEC,GRDOLD,GAMMA,HESOLD,'UPPSB ')
            IF (REDINT .OR. DELINT) THEN
               CALL UPPSB(MXRCRD,STPINT,GAMMA,HESOLD,HESINT,
     &              IINTCR,MXRCRD,TMPMAT,TMPMT2,TMPMT3,IPRINT)
            ELSE
               CALL UPPSB(MXRCRD,STPSYM,GAMMA,HESOLD,HESMOL,
     &              NCART,MXCOOR,TMPMAT,TMPMT2,TMPMT3,IPRINT)
            END IF
         ELSE IF (RANKON) THEN 
            CALL UPOUTP(MXRCRD,SCLVEC,GRDOLD,GAMMA,HESOLD,'UPRNKO')
            IF (REDINT .OR. DELINT) THEN
               CALL UPRNKO(MXRCRD,STPINT,GAMMA,HESOLD,HESINT,
     &              IINTCR,MXRCRD,TMPMAT,TMPMT2,TMPMT3,IPRINT)
            ELSE
               CALL UPRNKO(MXRCRD,STPSYM,GAMMA,HESOLD,HESMOL,
     &              NCART,MXCOOR,TMPMAT,TMPMT2,TMPMT3,IPRINT)
            END IF
         ELSE IF (BOFILL) THEN 
            CALL UPOUTP(MXRCRD,SCLVEC,GRDOLD,GAMMA,HESOLD,'UPBOFL')
            IF (REDINT .OR. DELINT) THEN
               CALL UPBOFL(MXRCRD,STPINT,GAMMA,HESOLD,HESINT,IINTCR,
     &              MXRCRD,TMPMAT,TMPMT2,TMPMT3,TMPMT4,TMPMT5,IPRINT)
            ELSE
               CALL UPBOFL(MXRCRD,STPSYM,GAMMA,HESOLD,HESMOL,NCART,
     &              MXCOOR,TMPMAT,TMPMT2,TMPMT3,TMPMT4,TMPMT5,IPRINT)
            END IF
         ELSE IF (BFGSR1) THEN
            CALL UPOUTP(MXRCRD,SCLVEC,GRDOLD,GAMMA,HESOLD,'UPBFR1')
            IF (REDINT .OR. DELINT) THEN
               CALL UPBFR1(MXRCRD,STPINT,GAMMA,HESOLD,HESINT,IINTCR,
     &              MXRCRD,TMPMAT,TMPMT2,TMPMT3,TMPMT4,TMPMT5,IPRINT)
            ELSE
               CALL UPBFR1(MXRCRD,STPSYM,GAMMA,HESOLD,HESMOL,NCART,
     &              MXCOOR,TMPMAT,TMPMT2,TMPMT3,TMPMT4,TMPMT5,IPRINT)
            END IF
         ELSE IF (SCHLEG) THEN
            CALL UPOUTP(MXRCRD,SCLVEC,GRDOLD,GAMMA,HESOLD,'UPSCHL')
            RESET = .FALSE.
            IF ((ITRNMR .LT. 2) .OR. (ITRBRK .EQ. (ITRNMR-2)))
     &           RESET = .TRUE.
            IF (REDINT .OR. DELINT) THEN
               CALL UPSCHL(MXRCRD,STPINT,STPMAT,GAMMA,GRDMAT,HESOLD,
     &              HESINT,IINTCR,MXRCRD,TMPMAT(1),TMPMAT(MX2CRD+1),
     &              TMPMT2,TMPMT3,IINTCR-NPROJ,RESET,GRADNM,IPRINT)
            ELSE
               CALL UPSCHL(MXRCRD,STPSYM,STPMAT,GAMMA,GRDMAT,HESOLD,
     &              HESMOL,NCART,MXCOOR,TMPMAT(1),TMPMAT(MX2CRD+1),
     &              TMPMT2,TMPMT3,NCART-NPROJ,RESET,GRADNM,IPRINT)
            END IF
         END IF
C
C     In all higher symmetries, we set all diagonal elements to 1
C
         IF ((.NOT. DELINT) .AND. (.NOT. REDINT)) THEN
            DO 15 I = (NCART + 1), NCRTOT
               HESMOL(I,I) = SCLVEC(I)*SCLVEC(I)
 15         CONTINUE
C
C     We copy the updated Hessian to HESOLD and the gradient to GRDOLD,
C     making everything ready for the next iteration. The updated
C     Hessian then has to be "off-scaled" before it's sent to MINCGH.
C
            CALL DZERO (HESOLD,MXRCRD*MXRCRD)
            CALL DZERO (GRDOLD,MXRCRD)
            DO 20 J = 1, NCART
               DO 22 I = 1, NCART
                  HESOLD(I,J) = HESMOL(I,J)
                  HESMOL(I,J) = HESMOL(I,J)*SCLVEC(I)*SCLVEC(J)
 22            CONTINUE
               GRDOLD(J) = GRDMOL(J)
 20         CONTINUE
         ELSE
            CALL DZERO (HESOLD,MXRCRD*MXRCRD)
            CALL DZERO (GRDOLD,MXRCRD)
            DO 30 J = 1, IINTCR
               DO 32 I = 1, IINTCR
                  HESOLD(I,J) = HESINT(I,J)
 32            CONTINUE
               GRDOLD(J) = GRDINT(J)
 30         CONTINUE
         END IF
      END IF
      CALL QEXIT('UPDHES')
      RETURN
      END

C  /* Deck upoutp */
      SUBROUTINE UPOUTP(MXRCRD,SCLVEC,GRDOLD,GAMMA,HESOLD,UPDTXT)
C
C     Some common output for the updating subroutines.
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "optinf.h"
#include "taymol.h"
#include "nuclei.h"
#include "symmet.h"
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      DIMENSION SCLVEC(MXCOOR), GRDOLD(MXRCRD) 
      DIMENSION GAMMA(MXRCRD), HESOLD(MXRCRD,MXRCRD) 
      CHARACTER UPDTXT*6, OUTTXT*18
C
      NCOOR = NCART
      IF (REDINT .OR. DELINT) NCOOR = IINTCR
      OUTTXT = 'Output from ' // UPDTXT
      IF (IPRINT .GE. IPRMIN) THEN
         CALL TITLER(OUTTXT,'*',103)
         IF ((IPRINT .GE. IPRDBG) .AND. (.NOT. REDINT)
     &        .AND. (.NOT. DELINT)) THEN
            CALL HEADER('Scaling vector in UPOUTP',-1)
            CALL OUTPUT(SCLVEC,1,1,1,NCART,1,MXCOOR,1,LUPRI)
         END IF
      END IF
      IF (IPRINT .GE. IPRMED) THEN
         CALL HEADER('Step from last geometry',-1)
         IF (REDINT .OR. DELINT) THEN
            CALL OUTPUT(STPINT,1,1,1,NCOOR,1,MXRCRD,1,LUPRI)
         ELSE
            CALL OUTPUT(STPSYM,1,1,1,NCOOR,1,MXCOOR,1,LUPRI)
         END IF
         CALL HEADER('Previous Hessian',-1)
         CALL OUTPUT(HESOLD,1,NCOOR,1,NCOOR,MXRCRD,MXRCRD,1,LUPRI)
         CALL HEADER('Gradient at last geometry',-1)
         CALL OUTPUT(GRDOLD,1,1,1,NCOOR,1,MXRCRD,1,LUPRI)
         CALL HEADER('Gradient at current geometry',-1)
         IF (REDINT .OR. DELINT) THEN
            CALL OUTPUT(GRDINT,1,1,1,NCOOR,1,MXRCRD,1,LUPRI)
         ELSE
            CALL OUTPUT(GRDMOL,1,1,1,NCOOR,1,MXCOOR,1,LUPRI)
         END IF
      END IF
C
C     The gradient difference (gamma) is calculated. The Cartesian
C     vector has to be scaled, to make it "compatible" with the step
C     vector and the Hessian.
C
      CALL DZERO(GAMMA,MXRCRD)
      DO 20 I = 1, NCOOR
         IF (REDINT .OR. DELINT) THEN
            GAMMA(I) = GRDINT(I) - GRDOLD(I)
         ELSE
            GAMMA(I) = (GRDMOL(I) - GRDOLD(I))/SCLVEC(I)
         END IF
 20   CONTINUE
      IF (IPRINT .GE. IPRMED) THEN
         CALL HEADER('Scaled gradient difference',-1)
         CALL OUTPUT(GAMMA,1,1,1,NCOOR,1,MXRCRD,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck updold */
      SUBROUTINE UPDOLD(MXRCRD,SCLVEC,GRDOLD,HESOLD,HESINT)
C
C     The previous Hessian is copied to MOLHES and GRDOLD is updated.
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "optinf.h"
#include "taymol.h"
#include "nuclei.h"
#include "symmet.h"
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      DIMENSION SCLVEC(MXCOOR), GRDOLD(MXRCRD)
      DIMENSION HESOLD(MXRCRD,MXRCRD), HESINT(MXRCRD,MXRCRD)
C
      IF (REDINT .OR. DELINT) THEN
         CALL DZERO(HESINT,MXRCRD*MXRCRD)
         CALL DZERO(GRDOLD,MXRCRD)
         DO 10 J = 1, IINTCR
            DO 15 I = 1, IINTCR
               HESINT(I,J) = HESOLD(I,J)
 15         CONTINUE
            GRDOLD(J) = GRDINT(J)
 10      CONTINUE
      ELSE
         CALL DZERO(HESMOL,MXCOOR*MXCOOR)
         CALL DZERO(GRDOLD,MXRCRD)
         DO 20 J = 1, NCART
            DO 25 I = 1, NCART
               HESMOL(I,J) = HESOLD(I,J)*SCLVEC(I)*SCLVEC(J)
 25         CONTINUE
            GRDOLD(J) = GRDMOL(J)
 20      CONTINUE
         DO 30 I = NCART+1, NCRTOT
            HESMOL(I,I) = SCLVEC(I)*SCLVEC(I)
 30      CONTINUE
      END IF
C
      IF (IPRINT .GE. IPRMED) THEN
         CALL TITLER('Output from UPDOLD','*',103)
         IF (IPRINT .GE. IPRMAX) THEN
            CALL HEADER('Unchanged Hessian',-1)
            CALL OUTPUT(HESOLD,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
            CALL HEADER('Gradient at current geometry',-1)
            CALL OUTPUT(GRDOLD,1,1,1,ICRD,1,MCRD,1,LUPRI)
         ELSE
            WRITE(LUPRI,*)
     &           'Hessian will not be updated in this iteration.'
         END IF
      END IF
      KEEPHE = .FALSE.
      RETURN
      END

C  /* Deck upbfgs */
      SUBROUTINE UPBFGS(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW,ICRD,MCRD,
     &                                 TMPVEC,TMPMAT,TMPMT2,IPRINT)
C
C     Updates the Hessian using BFGS method
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      DIMENSION DELTA(MXRCRD), GAMMA(MXRCRD)
      DIMENSION HESOLD(MXRCRD,MXRCRD), HESNEW(MCRD,MCRD)
      DIMENSION TMPVEC(MCRD), TMPMAT(MCRD,MCRD), TMPMT2(MCRD,MCRD)
C
C     ******************************************************************
C     The BFGS formula is:
C
C                                      T                         T
C                     gamma(n)*gamma(n)    B(n)*delta(n)*delta(n) *B(n)
C     B(n+1) = B(n) + ------------------ - ----------------------------
C                             T                       T
C                     gamma(n) *delta(n)      delta(n) *B(n)*delta(n)
C
C     where
C               gamma(n) = grad(n+1) - grad(n)
C               delta(n) = x(n+1) - x(n)
C
C     The terms in the formula is evaluated one by one below.
C     ******************************************************************
C
C     First we calculate (gamma^T*delta)
C
      FAC = D0
      DO 10 I = 1, ICRD
         FAC = FAC + GAMMA(I)*DELTA(I)
 10   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('(gamma^T*delta)',-1)
         WRITE(LUPRI,'(A,F16.6)') 'Value :    ', FAC
      END IF
C
C     (gamma*gamma^T)/(gamma^T*delta) is calculated
C
      CALL DZERO(TMPMAT,MCRD*MCRD)
      DO 20 J = 1, ICRD
         DO 22 I = 1, ICRD
            TMPMAT(I,J) = GAMMA(I)*GAMMA(J)/FAC
 22      CONTINUE
 20   CONTINUE
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('(gamma*gamma^T)/(gamma^T*delta)',-1)
         CALL OUTPUT(TMPMAT,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     The first two terms of the BFGS formula is placed in TMPMT2.
C
      CALL DZERO(TMPMT2,MCRD*MCRD)
      DO 30 J = 1, ICRD
         DO 32 I = 1, ICRD
            TMPMT2(I,J) = HESOLD(I,J) + TMPMAT(I,J)
 32      CONTINUE
 30   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Sum of first two terms',-1)
         CALL OUTPUT(TMPMT2,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     We place (B*delta) in TMPVEC and calculate (delta^T*B*delta)
C
      CALL DZERO(TMPVEC,MCRD)
      FAC = D0
      DO 40 I = 1, ICRD
         TMP = D0
         DO 42 J = 1, ICRD
            TMP = TMP + HESOLD(I,J)*DELTA(J)
 42      CONTINUE
         TMPVEC(I) = TMP
         FAC = FAC + DELTA(I)*TMP
 40   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('(H*delta)',-1)
         CALL OUTPUT(TMPVEC,1,1,1,ICRD,1,MCRD,1,LUPRI)
         CALL HEADER('(delta^T*B*delta)',-1)
         WRITE(LUPRI,'(A,F16.6)') 'Value :    ', FAC
      END IF
C
C     Since the Hessian is symmetric, the elements in the row vector
C     (delta^T*B) must be equal to the elements in the column vector
C     (B*delta). Calculation of (B*delta*delta^T*B)/(delta^T*B*delta)
C     is therefore quite simple.
C
      CALL DZERO(TMPMAT,MCRD*MCRD)
      DO 50 J = 1, ICRD
         DO 52 I = 1, ICRD
            TMPMAT(I,J) = TMPVEC(I)*TMPVEC(J)/FAC
 52      CONTINUE
 50   CONTINUE
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('(B*delta*delta^T*B)/(delta^T*B*delta)',-1)
         CALL OUTPUT(TMPMAT,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     Finally we obtain the updated Hessian
C
      DO 60 J = 1, ICRD
         DO 62 I = 1, ICRD
            HESNEW(I,J) = TMPMT2(I,J) - TMPMAT(I,J)
 62      CONTINUE
 60   CONTINUE
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Updated Hessian',-1)
         CALL OUTPUT(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck updfp */
      SUBROUTINE UPDFP(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW,ICRD,MCRD,
     &                                 TMPVEC,TMPMAT,TMPMT2,IPRINT)
C
C     Updates the Hessian using DFP method
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      DIMENSION DELTA(MXRCRD), GAMMA(MXRCRD)
      DIMENSION HESOLD(MXRCRD,MXRCRD), HESNEW(MCRD,MCRD)
      DIMENSION TMPVEC(MCRD), TMPMAT(MCRD,MCRD), TMPMT2(MCRD,MCRD)
C
C     ******************************************************************
C     The DFP formula is:
C
C                                  T                                   T
C                          delta(n) *B(n)*delta(n)    gamma(n)*gamma(n)
C     B(n+1) = B(n) + (1 + -----------------------) * ------------------
C                                    T                        T
C                            gamma(n) *delta(n)       gamma(n) *delta(n)
C
C                               T                              T
C              gamma(n)*delta(n) *B(n) + B(n)*delta(n)*gamma(n)
C            - -------------------------------------------------
C                                    T
C                            gamma(n) *delta(n)
C
C     where
C               gamma(n) = grad(n+1) - grad(n)
C               delta(n) = x(n+1) - x(n)
C
C     The terms in the formula is evaluated one by one below.
C     ******************************************************************
C
C     First we calculate (gamma^T*delta)
C
      FAC1 = D0
      DO 10 I = 1, ICRD
         FAC1 = FAC1 + GAMMA(I)*DELTA(I)
 10   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('(gamma^T*delta)',-1)
         WRITE(LUPRI,'(A,F16.6)') 'Value :    ', FAC1
      END IF
C
C     (delta^T*B*delta) is calculated
C
      FAC2 = D0
      DO 20 I = 1, ICRD
         DUMMY = D0
         DO 22 J = 1, ICRD
            DUMMY = DUMMY + HESOLD(I,J)*DELTA(J)
 22      CONTINUE
         FAC2 = FAC2 + DUMMY*DELTA(I)
 20   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('(delta^T*B*delta)',-1)
         WRITE(LUPRI,'(A,F16.6)') 'Value :    ', FAC2
      END IF
      FAC2 = 1.0D0 + FAC2/FAC1
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('(1+(delta^T*B*delta)/(gamma^T*delta))',-1)
         WRITE(LUPRI,'(A,F16.6)') 'Value :    ', FAC2
      END IF
C
C     We calculate (gamma*gamma^T)/(gamma^T*delta)
C
      CALL DZERO(TMPMAT,MCRD*MCRD)
      DO 25 J = 1, ICRD
         DO 27 I = 1, ICRD
            TMPMAT(I,J) = GAMMA(I)*GAMMA(J)/FAC1
 27      CONTINUE
 25   CONTINUE
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('(gamma*gamma^T)/(gamma^T*delta)',-1)
         CALL OUTPUT(TMPMAT,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     The first two terms of the DFP formula is placed in TMPMT2.
C
      CALL DZERO(TMPMT2,MCRD*MCRD)
      DO 30 J = 1, ICRD
         DO 32 I = 1, ICRD
            TMPMT2(I,J) = HESOLD(I,J) + FAC2*TMPMAT(I,J)
 32      CONTINUE
 30   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Sum of first two terms',-1)
         CALL OUTPUT(TMPMT2,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     Next is (delta^T*B)
C
      CALL DZERO(TMPVEC,MCRD)
      DO 40 J = 1, ICRD
         DO 42 I = 1, ICRD
            TMPVEC(J) = TMPVEC(J) + DELTA(I)*HESOLD(I,J)
 42      CONTINUE
 40   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('(delta^T*B)',-1)
         CALL OUTPUT(TMPVEC,1,1,1,ICRD,1,MCRD,1,LUPRI)
      END IF
C
C     Then (gamma*delta^T*B + B*delta*gamma^T)
C
      CALL DZERO(TMPMAT,MCRD*MCRD)
      DO 45 J = 1, ICRD
         DO 47 I = 1, ICRD
            TMPMAT(I,J) = GAMMA(I)*TMPVEC(J) + TMPVEC(I)*GAMMA(J)
 47      CONTINUE
 45   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('(gamma*delta^T*B + B*delta*gamma^T)',-1)
         CALL OUTPUT(TMPMAT,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     And then (gamma*delta^T*B + B*delta*gamma^T)/(gamma^T*delta)
C
      DO 50 J = 1, ICRD
         DO 52 I = 1, ICRD
            TMPMAT(I,J) = TMPMAT(I,J)/FAC1
 52      CONTINUE
 50   CONTINUE
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('(gamma*delta^T*B + B*delta*gamma^T)/' //
     &                                           '(gamma^T*delta)',-1)
         CALL OUTPUT(TMPMAT,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     Finally we obtain the updated Hessian
C
      DO 60 J = 1, ICRD
         DO 62 I = 1, ICRD
            HESNEW(I,J) = TMPMT2(I,J) - TMPMAT(I,J)
 62      CONTINUE
 60   CONTINUE
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Updated Hessian',-1)
         CALL OUTPUT(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck uppsb */
      SUBROUTINE UPPSB(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW,ICRD,MCRD,
     &                                 TMPVEC,TMPMAT,TMPMT2,IPRINT)
C
C     Updates the Hessian using Powell-symmetric-Broyden
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      DIMENSION DELTA(MXRCRD), GAMMA(MXRCRD)
      DIMENSION HESOLD(MXRCRD,MXRCRD), HESNEW(MCRD,MCRD)
      DIMENSION TMPVEC(MCRD), TMPMAT(MCRD,MCRD), TMPMT2(MCRD,MCRD)
C
C     **********************************************************************
C     The Powell-symmetric-Broyden formula is:
C
C                               T                        T
C                      (delta(n) *delta(n))*T(n)*delta(n)
C     B(n+1) = B(n) +  -----------------------------------
C                                     T          2
C                            (delta(n) *delta(n))
C
C                T                        T      T                            T
C       (delta(n) *delta(n))*delta(n)*T(n) -(T(n) *delta(n))*delta(n)*delta(n)
C     + -----------------------------------------------------------------------
C                                        T          2
C                               (delta(n) *delta(n))
C
C     where
C               T(n)     = gamma(n) - B(n)*delta(n)
C               gamma(n) = grad(n+1) - grad(n)
C               delta(n) = x(n+1) - x(n)
C
C     **********************************************************************
C
C     First we calculate (T = gamma-B*delta)
C
      CALL DZERO(TMPVEC,MCRD)
      DO 10 I = 1, ICRD
         DUMMY = D0
         DO 12 J = 1, ICRD
            DUMMY = DUMMY + HESOLD(I,J)*DELTA(J)
 12      CONTINUE
         TMPVEC(I) = GAMMA(I) - DUMMY
 10   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('(T = gamma-B*delta)',-1)
         CALL OUTPUT(TMPVEC,1,1,1,ICRD,1,MCRD,1,LUPRI)
      END IF
C
C     We calculate (delta^T*delta)
C
      FAC1 = D0
      DO 20 I = 1, ICRD
         FAC1 = FAC1 + DELTA(I)*DELTA(I)
 20   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('(delta^T*delta)',-1)
         WRITE(LUPRI,'(A,F16.6)') 'Value :    ', FAC1
      END IF
C
C     Then (T^T*delta)
C
      FAC2 = D0
      DO 22 I = 1, ICRD
         FAC2 = FAC2 + TMPVEC(I)*DELTA(I)
 22   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('(T^T*delta)',-1)
         WRITE(LUPRI,'(A,F16.6)') 'Value :    ', FAC2
      END IF
C
C     and (delta^T*delta)^2
C
      FAC3 = FAC1*FAC1
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('(delta^T*delta)^2',-1)
         WRITE(LUPRI,'(A,G16.6)') 'Value :    ', FAC3
      END IF
C
C     Then we construct the complete term
C
      CALL DZERO(TMPMAT,MCRD*MCRD)
      DO 30 J = 1, ICRD
         DO 32 I = 1, ICRD
            TMPMAT(I,J) = (  FAC1*TMPVEC(I)*DELTA(J)
     &                     + FAC1*DELTA(I)*TMPVEC(J)
     &                     - FAC2*DELTA(I)*DELTA(J) )/FAC3
 32      CONTINUE
 30   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Second term of formula',-1)
         CALL OUTPUT(TMPMAT,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     Finally we obtain the updated Hessian
C
      DO 40 J = 1, ICRD
         DO 45 I = 1, ICRD
            HESNEW(I,J) = HESOLD(I,J) + TMPMAT(I,J)
 45      CONTINUE
 40   CONTINUE
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Updated Hessian',-1)
         CALL OUTPUT(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck upbofl */
      SUBROUTINE UPBOFL(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW,ICRD,MCRD,
     &     TMPVEC,TMPMAT,TMPMT2,TMPMT3,TMPMT4,IPRINT)
C
C     Updates the Hessian using Bofills update
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      DIMENSION DELTA(MXRCRD), GAMMA(MXRCRD)
      DIMENSION HESOLD(MXRCRD,MXRCRD), HESNEW(MCRD,MCRD)
      DIMENSION TMPVEC(MCRD), TMPMAT(MCRD,MCRD), TMPMT2(MCRD,MCRD)
      DIMENSION TMPMT3(MCRD,MCRD), TMPMT4(MCRD,MCRD)
C
C
C     Bofills update is a linear combination of rank one and PSB:
C
C     B(n+1) = phi B_R1 + (1-phi) B_PSB
C
      CALL UPRNKO(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW,
     &     ICRD,MCRD,TMPMAT,TMPMT2,TMPMT3,-1)
      CALL UPPSB(MXRCRD,DELTA,GAMMA,HESOLD,TMPMT4,
     &     ICRD,MCRD,TMPMAT,TMPMT2,TMPMT3,-1)
      CALL DZERO(TMPVEC,MCRD)
      DO 10 I = 1, ICRD
         DO 12 J = 1, ICRD
            TMPVEC(I) = TMPVEC(I) + HESOLD(I,J)*DELTA(J)
 12      CONTINUE
         TMPVEC(I) = GAMMA(I) - TMPVEC(I)
 10   CONTINUE
C
C     We calculate the factor phi
C
      CPHI = DDOT(ICRD,DELTA,1,TMPVEC,1)
      CPHI = 1.0D0-((CPHI*CPHI)/(DDOT(ICRD,DELTA,1,DELTA,1)*
     &     DDOT(ICRD,TMPVEC,1,TMPVEC,1)))
      R1PHI = 1.0D0-CPHI
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Old Hessian',-1)
         CALL OUTPUT(HESOLD,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
         CALL HEADER('Rank one update',-1)
         CALL OUTPUT(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
         CALL HEADER('Powell update',-1)
         CALL OUTPUT(TMPMT4,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
         CALL HEADER('delta',-1)
         CALL OUTPUT(DELTA,1,1,1,ICRD,1,MCRD,1,LUPRI)
         CALL HEADER('gamma',-1)
         CALL OUTPUT(GAMMA,1,1,1,ICRD,1,MCRD,1,LUPRI)
         CALL HEADER('ksi',-1)
         CALL OUTPUT(TMPVEC,1,1,1,ICRD,1,MCRD,1,LUPRI)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'phi:',CPHI
         WRITE(LUPRI,*)
      END IF
      DO 20 I = 1, ICRD
         DO 22 J = 1, ICRD
            HESNEW(I,J) = R1PHI*HESNEW(I,J)+CPHI*TMPMT4(I,J)
 22      CONTINUE
 20   CONTINUE
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Updated Hessian',-1)
         CALL OUTPUT(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck upbfr1 */
      SUBROUTINE UPBFR1(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW,ICRD,MCRD,
     &     TMPVEC,TMPMAT,TMPMT2,TMPMT3,TMPMT4,IPRINT)
C
C     Updates the Hessian using BFGS/rank one combination update (ala Bofill)
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      DIMENSION DELTA(MXRCRD), GAMMA(MXRCRD)
      DIMENSION HESOLD(MXRCRD,MXRCRD), HESNEW(MCRD,MCRD)
      DIMENSION TMPVEC(MCRD), TMPMAT(MCRD,MCRD), TMPMT2(MCRD,MCRD)
      DIMENSION TMPMT3(MCRD,MCRD), TMPMT4(MCRD,MCRD)
C
C     Use a linear combination of rank one and BFGS:
C
C     B(n+1) = phi B_R1 + (1-phi) B_BFGS
C
      CALL UPRNKO(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW,
     &     ICRD,MCRD,TMPMAT,TMPMT2,TMPMT3,-1)
      CALL UPBFGS(MXRCRD,DELTA,GAMMA,HESOLD,TMPMT4,
     &     ICRD,MCRD,TMPMAT,TMPMT2,TMPMT3,-1)
      CALL DZERO(TMPVEC,MCRD)
      DO 10 I = 1, ICRD
         DO 12 J = 1, ICRD
            TMPVEC(I) = TMPVEC(I) + HESOLD(I,J)*DELTA(J)
 12      CONTINUE
         TMPVEC(I) = GAMMA(I) - TMPVEC(I)
 10   CONTINUE
C
C     We calculate the factor phi
C
      CPHI = DDOT(ICRD,DELTA,1,TMPVEC,1)
      CPHI = 1.0D0-((CPHI*CPHI)/(DDOT(ICRD,DELTA,1,DELTA,1)*
     &     DDOT(ICRD,TMPVEC,1,TMPVEC,1)))
      R1PHI = 1.0D0-CPHI
C
C     Output for debugging
C
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Old Hessian',-1)
         CALL OUTPUT(HESOLD,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
         CALL HEADER('Rank one update',-1)
         CALL OUTPUT(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
         CALL HEADER('BFGS update',-1)
         CALL OUTPUT(TMPMT4,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
         CALL HEADER('delta',-1)
         CALL OUTPUT(DELTA,1,1,1,ICRD,1,MCRD,1,LUPRI)
         CALL HEADER('gamma',-1)
         CALL OUTPUT(GAMMA,1,1,1,ICRD,1,MCRD,1,LUPRI)
         CALL HEADER('ksi',-1)
         CALL OUTPUT(TMPVEC,1,1,1,ICRD,1,MCRD,1,LUPRI)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'phi:',CPHI
         WRITE(LUPRI,*)
      END IF
      DO 20 I = 1, ICRD
         DO 22 J = 1, ICRD
            HESNEW(I,J) = R1PHI*HESNEW(I,J)+CPHI*TMPMT4(I,J)
 22      CONTINUE
 20   CONTINUE
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Updated Hessian',-1)
         CALL OUTPUT(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck uprnko */
      SUBROUTINE UPRNKO(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW,ICRD,MCRD,
     &                                 TMPVEC,TMPMAT,TMPMT2,IPRINT)
C
C     Updates the Hessian using rank one method
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      DIMENSION DELTA(MXRCRD), GAMMA(MXRCRD)
      DIMENSION HESOLD(MXRCRD,MXRCRD), HESNEW(MCRD,MCRD)
      DIMENSION TMPVEC(MCRD), TMPMAT(MCRD,MCRD), TMPMT2(MCRD,MCRD)
C
C     ******************************************************************
C     The rank one formula is:
C
C                                                                     T
C                     (gamma(n)-B(n)*delta(n))(gamma(n)-B(n)*delta(n))
C     B(n+1) = B(n) + -------------------------------------------------
C                                                   T
C                           (gamma(n)-B(n)*delta(n)) *delta(n)
C
C     where
C               gamma(n) = grad(n+1) - grad(n)
C               delta(n) = x(n+1) - x(n)
C
C     The terms in the formula is evaluated one by one below.
C     ******************************************************************
C
C     First we calculate (gamma-B*delta)
C
      CALL DZERO(TMPVEC,MCRD)
      DO 10 I = 1, ICRD
         DUMMY = D0
         DO 12 J = 1, ICRD
            DUMMY = DUMMY + HESOLD(I,J)*DELTA(J)
 12      CONTINUE
         TMPVEC(I) = GAMMA(I) - DUMMY
 10   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('(gamma-B*delta)',-1)
         CALL OUTPUT(TMPVEC,1,1,1,ICRD,1,MCRD,1,LUPRI)
      END IF
C
C     (gamma-B*delta)(gamma-B*delta)^T is calculated
C
      CALL DZERO(TMPMAT,MCRD*MCRD)
      DO 20 J = 1, ICRD
         DO 22 I = 1, ICRD
            TMPMAT(I,J) = TMPVEC(I)*TMPVEC(J)
 22      CONTINUE
 20   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('(gamma-B*delta)(gamma-B*delta)^T',-1)
         CALL OUTPUT(TMPMAT,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     ((gamma-B*delta)^T*delta) is calculated
C
      FAC = D0
      DO 25 I = 1, ICRD
         FAC = FAC + TMPVEC(I)*DELTA(I)
 25   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('((gamma-B*delta)^T*delta)',-1)
         WRITE(LUPRI,'(A,F16.6)') 'Value :    ', FAC 
      END IF
C
C     And we complete the term
C     ((gamma-B*delta)(gamma-B*delta)^T)/((gamma-B*delta)^T*delta)
C
      DO 30 J = 1, ICRD
         DO 32 I = 1, ICRD
            TMPMAT(I,J) = TMPMAT(I,J)/FAC
 32      CONTINUE
 30   CONTINUE
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('((gamma-B*delta)(gamma-B*delta)^T)/' //
     &                                 '((gamma-B*delta)^T*delta)',-1)
         CALL OUTPUT(TMPMAT,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     Finally we obtain the updated Hessian
C
      CALL DZERO(HESNEW,MCRD*MCRD)
      DO 40 J = 1, ICRD
         DO 42 I = 1, ICRD
            HESNEW(I,J) = HESOLD(I,J) + TMPMAT(I,J)
 42      CONTINUE
 40   CONTINUE
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Updated Hessian',-1)
         CALL OUTPUT(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck upschl */
      SUBROUTINE UPSCHL(MXRCRD,DELTA,STPMAT,GAMMA,GRDMAT,HESOLD,HESNEW,
     &     ICRD,MCRD,TMPVEC,TMPVC2,RMAT,BMAT,MAXDIM,RESET,GNRM,IPRINT)
C
C     Updates the Hessian using Schlegel's updating scheme
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      DIMENSION DELTA(MXRCRD), STPMAT(25,MCRD)
      DIMENSION GAMMA(MXRCRD), GRDMAT(25,MCRD)
      DIMENSION HESOLD(MXRCRD,MXRCRD), HESNEW(MCRD,MCRD)
      DIMENSION TMPVEC(MCRD), TMPVC2(MCRD)
      DIMENSION RMAT(MCRD,MCRD), BMAT(MCRD,MCRD)
      LOGICAL RESET
C
      CALL DZERO(TMPVEC,MCRD)
      CALL DZERO(TMPVC2,MCRD)
      CALL DZERO(RMAT,MCRD*MCRD)
      CALL DZERO(BMAT,MCRD*MCRD)
C
C     First we have to transfer the last step and gradient difference to
C     STPMAT and GRDMAT respectively.
C
      IF (RESET) THEN
         CALL DZERO(GRDMAT,25*MCRD)
         CALL DZERO(STPMAT,25*MCRD)
         INUM = 1
C
C     After 25 iterations we have to discard the first entries.
C
      ELSE IF (STPMAT(25,1) .GT. 1.0D10) THEN
         DO 10 I = 1, 24
            DO 12 J = 1, ICRD
               STPMAT(I,J) = STPMAT(I+1,J)
               GRDMAT(I,J) = GRDMAT(I+1,J)
 12         CONTINUE
 10      CONTINUE
         STPMAT(25,1) = 0.0D0
         INUM = 24
      ELSE
         INUM = 1
 14      CONTINUE
         IF (STPMAT(INUM,1) .LT. 1.0D10) THEN
            INUM = INUM + 1
            GOTO 14
         END IF
C
C     We also have to check the dimension of our variational space,
C     if the number of displacement vectors exceeds this number, we
C     have to remove some.
C
         IF (INUM .GT. MAXDIM) THEN
            DO 15 I = 1, INUM -1
               DO 16 J = 1, ICRD
                  STPMAT(I,J) = STPMAT(I+1,J)
                  GRDMAT(I,J) = GRDMAT(I+1,J)
 16            CONTINUE
 15         CONTINUE
            INUM = INUM - 1
         END IF
      END IF
C
C     Then we update all vectors
C
      SNMLST = 0.0D0
      DO 17 I = 1, ICRD
         STPMAT(INUM,I) = -DELTA(I)
         GRDMAT(INUM,I) = -GAMMA(I)
         SNMLST = SNMLST + DELTA(I)*DELTA(I)
         DO 19 II = 1, INUM-1
            STPMAT(II,I) = STPMAT(II,I)-DELTA(I)
            GRDMAT(II,I) = GRDMAT(II,I)-GAMMA(I)
 19      CONTINUE
 17   CONTINUE
      SNMLST = SQRT(SNMLST)
      STPMAT(INUM+1,1) = 1.1D10
C
C     **********************************************************************
C     Schlegel's scheme goes as follows:
C
C                    i-1            t
C     r' = (x -x )-  SUM  r ((x -x ) r ) ;   r = r'/|r'| ;   j = i-1, i-2, ...
C      j     j  i   m=j+1  m   j  i   m       j   j   j
C
C                    t      i-1             t               t
C     b   = [ (g -g ) r  -  SUM  b  ((x -x ) r ) ] / (x -x ) r  ;
C      jk       j  i   k   m=j+1  mk   j  i   m        j  i   j
C
C     b   = b   ;   j <= k = i-1, i-2, ...
C      kj    jk
C
C                                       t
C     B  = B   + SUM (b  - r B   r ) r r
C      i    i-1  j k   jk   j i-1 k   j k
C
C     **********************************************************************
C
C     We start by calculating the orthonormal displacement vectors r(j)
C
      DO 20 II = INUM, 1, -1
         TMPNRM = 0.0D0
         DO 22 I = 1, ICRD
            RMAT(II,I) = STPMAT(II,I)
            TMPNRM = TMPNRM + STPMAT(II,I)*STPMAT(II,I)
            TMPVC2(I) = RMAT(II,I)
 22      CONTINUE
         TMPNRM = SQRT(TMPNRM)
         DO 26 I = II+1, INUM
            TMP = 0.0D0
            DO 28 J = 1, ICRD
               TMP = TMP + STPMAT(II,J)*RMAT(I,J)
 28         CONTINUE
            IF(IPRINT .GE. IPRDBG) THEN
               CALL HEADER('(x(j)-x(i))^T * r(m)',-1)
               WRITE(LUPRI,*) TMP
            END IF
            DO 30 J = 1, ICRD
               RMAT(II,J) = RMAT(II,J) - RMAT(I,J)*TMP
               TMPVC2(J) = RMAT(II,J)
 30         CONTINUE
 26      CONTINUE
         IF(IPRINT .GE. IPRDBG) THEN
            CALL HEADER('r(j)''',-1)
            CALL OUTPUT(TMPVC2,1,1,1,ICRD,1,MCRD,1,LUPRI)
         END IF
         TMP = 0.0D0
         DO 32 J = 1, ICRD
            TMP = TMP + RMAT(II,J)*RMAT(II,J)
 32      CONTINUE
         TMP = SQRT(TMP)
C
C     If a vector is in a space already spanned by other vectors
C     (that is less than 15% of the vector is outside this space),
C     all elements are set equal to a very small number.
C
         IF ((TMP .LT. 0.15D0*TMPNRM) .OR.
     &        (TMPNRM .GE. 1.0D2*SNMLST) .OR.
     &        (GNRM .GE. 0.1D0)) THEN
            GRDMAT(II,1) = 1.1D10
            DO 35 J = 1, ICRD
               RMAT(II,J) = 1.0D-25
 35         CONTINUE
            TMP = 1.0D0
         END IF
         DO 40 J = 1, ICRD
            RMAT(II,J) = RMAT(II,J)/TMP
            TMPVC2(J)  = RMAT(II,J)
 40      CONTINUE
         IF (IPRINT .GE. IPRDBG) THEN
            CALL HEADER('Original displacement vector',-1)
            CALL OUTPUT(STPMAT,II,II,1,ICRD,25,MXRCRD,1,LUPRI)
            CALL HEADER('Orthonormalized vector',-1)
            CALL OUTPUT(TMPVC2,1,1,1,ICRD,1,MCRD,1,LUPRI)
         END IF
 20   CONTINUE
C
C     All redundant vectors are removed
C
      DO 45 I = INUM, 1, -1
         IF (GRDMAT(I,1) .GT. 1.0D10) THEN
            DO 47 II = I, INUM
               DO 48 III = 1, ICRD
                  RMAT(II,III)   = RMAT(II+1,III)
                  GRDMAT(II,III) = GRDMAT(II+1,III)
                  STPMAT(II,III) = STPMAT(II+1,III)
 48            CONTINUE
 47         CONTINUE
            INUM = INUM -1
         END IF
 45   CONTINUE
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('Matrix of orthonormalized vectors',-1)
         CALL OUTPUT(RMAT,1,INUM,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     Next task is to calculate the matrix elements b(jk)
C
      DO 50 K = INUM, 1, -1
         DO 52 J = K, 1, -1
            IF(IPRINT .GE. IPRDBG) THEN
               CALL HEADER('g(j)-g(i)',-1)
               CALL OUTPUT(GRDMAT,J,J,1,ICRD,25,MCRD,1,LUPRI)
            END IF
            BMAT(J,K) = 0.0D0
            DO 58 I = 1, ICRD
               BMAT(J,K) = BMAT(J,K) + GRDMAT(J,I)*RMAT(K,I)
 58         CONTINUE
C
            IF(IPRINT .GE. IPRDBG) THEN
               CALL HEADER('x(j)-x(i)',-1)
               CALL OUTPUT(STPMAT,J,J,1,ICRD,25,MCRD,1,LUPRI)
            END IF
C
            DO 70 M = J+1, INUM
               TMP = 0.0D0
               DO 72 I = 1, ICRD
                  TMP = TMP + STPMAT(J,I)*RMAT(M,I)
 72            CONTINUE
               BMAT(J,K) = BMAT(J,K) - BMAT(M,K)*TMP
 70         CONTINUE
            TMP = 0.0D0
            DO 75 I = 1, ICRD
               TMP = TMP + STPMAT(J,I)*RMAT(J,I)
 75         CONTINUE
            BMAT(J,K) = BMAT(J,K)/TMP
            BMAT(K,J) = BMAT(J,K)
 52      CONTINUE
 50   CONTINUE
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('B-matrix of coefficients',-1)
         CALL OUTPUT(BMAT,1,INUM,1,INUM,MCRD,MCRD,1,LUPRI)
      END IF
C
C     And finally we are ready to update the Hessian
C
C     The elements b(jk) are replaced by b(jk)-r(j)^T*B*r(k)
C
      DO 80 K = 1, INUM
         DO 82 I = 1, ICRD
            TMPVEC(I) = 0.0D0
            DO 83 II = 1, ICRD
               TMPVEC(I) = TMPVEC(I) + HESOLD(I,II)*RMAT(K,II)
 83         CONTINUE
 82      CONTINUE
         IF(IPRINT .GE. IPRDBG) THEN
            CALL HEADER('B*r(k)^T',-1)
            CALL OUTPUT(TMPVEC,1,1,1,ICRD,1,MCRD,1,LUPRI)
         END IF
         DO 85 J = 1, INUM
            TMP = 0.0D0
            DO 87 I = 1, ICRD
               TMP = TMP + RMAT(J,I)*TMPVEC(I)
 87         CONTINUE
            BMAT(J,K) = BMAT(J,K) - TMP
 85      CONTINUE
 80   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('b(jk)-r(j)^T*B*r(k)',-1)
         CALL OUTPUT(BMAT,1,INUM,1,INUM,MCRD,MCRD,1,LUPRI)
      END IF
C
C     The old Hessian is moved to HESNEW, so that HESOLD can be used for
C     temporary storage.
C
      CALL DZERO(HESNEW,MCRD*MCRD)
      DO 90 J = 1, ICRD
         DO 92 I = 1, ICRD
            HESNEW(I,J) = HESOLD(I,J)
 92      CONTINUE
 90   CONTINUE
      CALL DZERO(TMPVEC,MCRD)
C
C     Then we construct the matrices r(j)*r(k)^T
C
      DO 100 J = 1, INUM
         DO 110 K = 1, INUM
            CALL DZERO(HESOLD,MXRCRD*MXRCRD)
            DO 114 IJ = 1, ICRD
               DO 116 IK = 1, ICRD
                  HESOLD(IJ,IK) = RMAT(J,IJ)*RMAT(K,IK)
 116           CONTINUE
 114        CONTINUE
            IF (IPRINT .GE. IPRDBG) THEN
               CALL HEADER('r(j)*r(k)^T',-1)
               CALL OUTPUT(HESOLD,1,ICRD,1,ICRD,MXRCRD,MXRCRD,
     &              1,LUPRI)
            END IF
            DO 120 II = 1, ICRD
               DO 122 JJ = 1, ICRD
                  HESNEW(II,JJ) = HESNEW(II,JJ)+BMAT(J,K)*HESOLD(II,JJ)
 122           CONTINUE
 120        CONTINUE
 110     CONTINUE
 100  CONTINUE
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Updated Hessian',-1)
         CALL OUTPUT(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck upmult */
      SUBROUTINE UPMULT(MXRCRD,DELTA,STPMAT,GAMMA,GRDMAT,HESOLD,HESNEW,
     &     ICRD,MCRD,TMPVEC,TMPVC2,TMPMAT,TMPMT2,HESUPD,RMAT,GMAT,
     &     MAXDIM,RESET,GNRM,IPRINT,ITYPE,SMART)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      DIMENSION DELTA(MXRCRD), STPMAT(25,MCRD)
      DIMENSION GAMMA(MXRCRD), GRDMAT(25,MCRD)
      DIMENSION HESOLD(MXRCRD,MXRCRD), HESNEW(MCRD,MCRD)
      DIMENSION TMPVEC(MCRD), TMPVC2(MCRD)
      DIMENSION TMPMAT(MXRCRD*MXRCRD), TMPMT2(MXRCRD*MXRCRD)
      DIMENSION HESUPD(MXRCRD,MXRCRD)
      DIMENSION RMAT(25,MCRD), GMAT(25,MCRD)
      LOGICAL RESET,SMART
C
      ILIM = 25
      IF (.NOT. SMART) ILIM = 5
C
C     ITYPE indicates the update method:
C
C     1 - BFGS (default)
C     2 - PSB
C     3 - DFP
C     4 - Rank one
C
      IF ((ITYPE .LT. 1) .OR. (ITYPE .GT. 4)) ITYPE = 1
C
C     Since STPMAT and GRDMAT are used to store a maximum of 25 vectors,
C     MCRD should be equal to or larger than 25
C
      MXDIM = MAXDIM
      IF (MCRD .LT. 25) MXDIM = MIN(MAXDIM,MCRD)
C
      CALL DZERO(TMPVEC,MCRD)
      CALL DZERO(TMPVC2,MCRD)
      CALL DZERO(GMAT,25*MCRD)
      CALL DZERO(RMAT,25*MCRD)
      CALL DZERO(HESUPD,MXRCRD*MXRCRD)
C
C     First we have to transfer the last step and gradient difference to
C     STPMAT and GRDMAT respectively.
C
      IF (RESET) THEN
         CALL DZERO(GRDMAT,25*MCRD)
         CALL DZERO(STPMAT,25*MCRD)
         INUM = 1
C
C     After ILIM iterations we discard the first entries.
C
      ELSE IF (STPMAT(ILIM,1) .GT. 1.0D10) THEN
         DO 10 I = 1, ILIM-1
            DO 12 J = 1, ICRD
               STPMAT(I,J) = STPMAT(I+1,J)
               GRDMAT(I,J) = GRDMAT(I+1,J)
 12         CONTINUE
 10      CONTINUE
         STPMAT(ILIM,1) = 0.0D0
         INUM = ILIM-1
      ELSE
         INUM = 1
 14      CONTINUE
         IF (STPMAT(INUM,1) .LT. 1.0D10) THEN
            INUM = INUM + 1
            GOTO 14
         END IF
C
C     We also have to check the dimension of our variational space,
C     if the number of displacement vectors exceeds this number, we
C     have to remove some.
C
         IF (INUM .GT. MXDIM) THEN
            DO 15 I = 1, INUM-1
               DO 16 J = 1, ICRD
                  STPMAT(I,J) = STPMAT(I+1,J)
                  GRDMAT(I,J) = GRDMAT(I+1,J)
 16            CONTINUE
 15         CONTINUE
            INUM = INUM - 1
            IF (IPRINT .GE. IPRDBG) THEN
               WRITE(LUPRI,*)
               WRITE(LUPRI,*) 'Removing one ' //
     &              'displacement due to dimension of variational space'
            END IF
         END IF
      END IF
C
C     Then we update all vectors and calculate the norm of the
C     last displacement (SNMLST). The vectors are also copied
C     to RMAT and GMAT.
C
      SNMLST = 0.0D0
      DO 17 I = 1, ICRD
         STPMAT(INUM,I) = -DELTA(I)
         GRDMAT(INUM,I) = -GAMMA(I)
         RMAT(INUM,I) = STPMAT(INUM,I)
         GMAT(INUM,I) = GRDMAT(INUM,I)
         SNMLST = SNMLST + DELTA(I)*DELTA(I)
         DO 19 II = 1, INUM-1
            STPMAT(II,I) = STPMAT(II,I)-DELTA(I)
            GRDMAT(II,I) = GRDMAT(II,I)-GAMMA(I)
            RMAT(II,I) = STPMAT(II,I)
            GMAT(II,I) = GRDMAT(II,I)
 19      CONTINUE
 17   CONTINUE
      SNMLST = SQRT(SNMLST)
      STPMAT(INUM+1,1) = 1.1D10
      RMAT(INUM+1,1) = 1.1D10
      CALL HEADER('Original displacements',-1)
      CALL OUTPUT(STPMAT,1,INUM,1,ICRD,25,MCRD,1,LUPRI)
      CALL HEADER('Original gradient vectors',-1)
      CALL OUTPUT(GRDMAT,1,INUM,1,ICRD,25,MCRD,1,LUPRI)
C
C     Displacements larger than 100 times the last step,
C     and displacements with a gradient difference
C     larger than 0.25D0 are marked (later to be removed).
C
      IF (SMART) THEN
         DO 20 II = 1, INUM
            DNRM = D0
            GNRM = D0
            DO 22 I = 1, ICRD
               DNRM = DNRM + RMAT(II,I)*RMAT(II,I)
               GNRM = GNRM + GMAT(II,I)*GMAT(II,I)
 22         CONTINUE
            DNRM = SQRT(DNRM)
            GNRM = SQRT(GNRM)
            IF (DNRM .GE. 100D0*SNMLST) THEN
C            IF (DNRM .GE. 500D0*SNMLST) THEN
               IF (IPRINT .GE. IPRDBG) WRITE(LUPRI,*)
     &              'Removing one displacement due to distance'
               RMAT(II,1) = -1.1D10
            ELSE IF (GNRM .GE. 0.25D0) THEN
               IF (IPRINT .GE. IPRDBG) WRITE(LUPRI,*)
     &              'Removing one displacement due to gradient norm'
               RMAT(II,1) = -1.1D10
            END IF
 20      CONTINUE
         II = 1
 25      CONTINUE
         IF ((RMAT(II,1) .LE. -1.0D10) .AND. (INUM .GE. 1)) THEN
            DO 27 JJ = II, INUM - 1
               DO 29 J = 1, ICRD
                  RMAT(JJ,J) = RMAT(JJ + 1,J)
                  GMAT(JJ,J) = GMAT(JJ + 1,J)
 29            CONTINUE
 27         CONTINUE
            INUM = INUM - 1
            GOTO 25
         ELSE IF (II .LT. INUM-1) THEN
            II = II + 1
            GOTO 25
         END IF
      END IF
C
C     We check if the last displacement is nearly parallell to an
C     earlier step (dot product larger than 0.75). If that is the case,
C     the older step is removed.
C
      IF ((INUM .GT. 1) .AND. SMART) THEN
         II = 1
 30      CONTINUE
         DOTP = D0
         SNRM1 = D0
         SNRM2 = D0
         DO 32 I = 1, ICRD
            DOTP = DOTP + RMAT(II,I)*RMAT(INUM,I)
            SNRM1 = SNRM1 + RMAT(II,I)*RMAT(II,I)
            SNRM2 = SNRM2 + RMAT(INUM,I)*RMAT(INUM,I)
 32      CONTINUE
         SNRM1 = SQRT(SNRM1)
         SNRM2 = SQRT(SNRM2)
         IF ((SNRM1*SNRM2) .GE. 1.0D-15) THEN
            DOTP = DOTP/(SNRM1*SNRM2)
         ELSE
            DOTP = D0
         END IF
C         IF (DOTP .GE. 0.80D0) THEN
         IF (DOTP .GE. 0.90D0) THEN
            DO 34 J = 1, ICRD
               DO 36 JJ = II, INUM - 1
                  RMAT(JJ,J) = RMAT(JJ+1,J)
                  GMAT(JJ,J) = GMAT(JJ+1,J)
 36            CONTINUE
 34         CONTINUE
            INUM = INUM - 1
            IF (IPRINT .GE. IPRDBG) WRITE(LUPRI,*)
     &           'Removing one nearly parallell and older displacement'
            IF (II .LT. INUM) GOTO 30
         ELSE IF (II .LT. INUM-1) THEN
            II = II + 1
            GOTO 30
         END IF
      END IF
C
C     Some output
C
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('Displacements',-1)
         CALL OUTPUT(RMAT,1,INUM,1,ICRD,25,MCRD,1,LUPRI)
         CALL HEADER('Gradient vectors',-1)
         CALL OUTPUT(GMAT,1,INUM,1,ICRD,25,MCRD,1,LUPRI)
      END IF
C
C     Finally we use do the updating, suppressing output.
C
      DO 60 K = 1, INUM
         CALL DZERO(DELTA,MXRCRD)
         CALL DZERO(GAMMA,MXRCRD)
         DO 62 I = 1, ICRD
             DELTA(I) = -RMAT(K,I)
             GAMMA(I) = -GMAT(K,I)
 62      CONTINUE
         IF (ITYPE .EQ. 1) THEN
            CALL UPBFGS(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW,
     &           ICRD,MCRD,TMPVEC,TMPMAT,TMPMT2,-1)
         ELSE IF (ITYPE .EQ. 2) THEN
            CALL UPPSB(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW,
     &           ICRD,MCRD,TMPVEC,TMPMAT,TMPMT2,-1)
         ELSE IF (ITYPE .EQ. 3) THEN
            CALL UPDFP(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW,
     &           ICRD,MCRD,TMPVEC,TMPMAT,TMPMT2,-1)
         ELSE
            CALL UPRNKO(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW,
     &           ICRD,MCRD,TMPVEC,TMPMAT,TMPMT2,-1)
         END IF
         IF ((IPRINT .GE. IPRDBG) .AND. (K .LT. INUM)) THEN
            CALL HEADER('Updating the Hessian',-1)
            CALL OUTPUT(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)         
         END IF
C
C     The various contributions are collected in UPDHES.
C
         DO 65 J = 1, ICRD
            DO 67 I = 1, ICRD
               HESUPD(I,J) = HESUPD(I,J) + (HESNEW(I,J) - HESOLD(I,J))
 67         CONTINUE
 65      CONTINUE
         IF (IPRINT .GE. IPRDBG) THEN
            CALL HEADER('Update',-1)
            CALL OUTPUT(HESUPD,1,ICRD,1,ICRD,MXRCRD,MXRCRD,1,LUPRI)
         END IF
 60   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Total Update',-1)
         CALL OUTPUT(HESUPD,1,ICRD,1,ICRD,MXRCRD,MXRCRD,1,LUPRI)
      END IF
      IF (INUM .GE. 1) THEN
         FAC = 1.0D0/(1.0D0*INUM)
         DO 70 J = 1, ICRD
            DO 72 I = 1, ICRD
               HESUPD(I,J) = HESUPD(I,J)*FAC
               HESNEW(I,J) = HESOLD(I,J) + HESUPD(I,J)
 72         CONTINUE
 70      CONTINUE
      END IF
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Scaled Update',-1)
         CALL OUTPUT(HESUPD,1,ICRD,1,ICRD,MXRCRD,MXRCRD,1,LUPRI)
      END IF
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Updated Hessian',-1)
         CALL OUTPUT(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck upcmbm */
      SUBROUTINE UPCMBM(MXRCRD,DELTA,STPMAT,GAMMA,GRDMAT,HESOLD,HESNEW,
     &     ICRD,IRCRD,MCRD,TMPVEC,TMPVC2,RMAT,GMAT,TMPMAT,TMPMT2,TMPMT3,
     &     TMPMT4,TMPMT5,BMTRAN,MAXDIM,RESET,REJLST,GNRM,
     &     ITYPE,DELINT,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      DIMENSION DELTA(MXRCRD), STPMAT(25,MCRD)
      DIMENSION GAMMA(MXRCRD), GRDMAT(25,MCRD)
      DIMENSION HESOLD(MXRCRD,MXRCRD), HESNEW(MCRD,MCRD)
      DIMENSION TMPVEC(MCRD), TMPVC2(MCRD)
      DIMENSION TMPMAT(MXRCRD,MXRCRD), TMPMT2(MXRCRD*MXRCRD)
      DIMENSION TMPMT3(MXRCRD,MXRCRD), TMPMT4(MXRCRD,MXRCRD)
      DIMENSION TMPMT5(MXRCRD,MXRCRD), BMTRAN(MXRCRD,MXRCRD)
      DIMENSION RMAT(MCRD,MCRD), GMAT(MCRD,MCRD)
      LOGICAL RESET, REJLST, DELINT
      SAVE IFAC
      DATA IFAC /2/
      IF (RESET) IFAC = 1
      IF (RESET .OR. REJLST) IFAC = 1
C
C     We do the updating, suppressing output.
C
      CALL UPBFGS(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW,
     &     ICRD,MCRD,TMPMAT,TMPMT2,TMPMT3,IPRINT)
C
C     A new model Hessian is calculated
C
      CALL BLDHES(MXRCRD,TMPMAT,TMPMT5)
C
C     If we use delocalized internals, we have to transform from
C     redundant to non-redundant space.
C
      IF (DELINT) THEN
         CALL DZERO(TMPMAT,MXRCRD*MXRCRD)
         DO 200 I = 1, ICRD
            DO 202 J = 1, IRCRD
               DO 204 K = 1, IRCRD
                  TMPMAT(I,J) = TMPMAT(I,J) + BMTRAN(K,I)*TMPMT5(K,J)
 204           CONTINUE
 202        CONTINUE
 200     CONTINUE
         CALL DZERO(TMPMT5,MXRCRD*MXRCRD)
         DO 210 I = 1, ICRD
            DO 212 J = 1, ICRD
               DO 214 K = 1, IRCRD
                  TMPMT5(I,J) = TMPMT5(I,J) + TMPMAT(I,K)*BMTRAN(K,J)
 214           CONTINUE
 212        CONTINUE
 210     CONTINUE
      END IF
      CALL UPBFGS(MXRCRD,DELTA,GAMMA,TMPMT5,TMPMT4,
     &     ICRD,MCRD,TMPMAT,TMPMT2,TMPMT3,IPRINT)
C
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Old Hessian',-1)
         CALL OUTPUT(HESOLD,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)         
         CALL HEADER('BFGS-updated Hessian',-1)
         CALL OUTPUT(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)         
         CALL HEADER('Model Hessian',-1)
         CALL OUTPUT(TMPMT4,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)         
      END IF
      IF (IFAC .GT. 20) THEN
         FAC  = 0.0D0
         FAC2 = 1.0D0
      ELSE
         FAC  = 1.0D0/(1.0D0*IFAC)
         FAC2 = 1.0D0 - FAC
      END IF
      DO 310 I = 1, ICRD
         DO 312 J = 1, ICRD
            HESNEW(I,J) = FAC*TMPMT4(I,J) + FAC2*HESNEW(I,J)
 312     CONTINUE
 310  CONTINUE
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Updated Hessian',-1)
         CALL OUTPUT(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      IFAC = MIN(IFAC*2,25)
      RETURN
      END

C  /* Deck upmodh */
      SUBROUTINE UPMODH(MXRCRD,DELTA,STPMAT,GAMMA,GRDMAT,HESOLD,HESNEW,
     &     ICRD,IRCRD,MCRD,TMPVEC,TMPVC2,RMAT,GMAT,TMPMAT,TMPMT2,TMPMT3,
     &     TMPMT4,TMPMT5,BMTRAN,MAXDIM,RESET,REJLST,GNRM,
     &     ITYPE,DELINT,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      DIMENSION DELTA(MXRCRD), STPMAT(25,MCRD)
      DIMENSION GAMMA(MXRCRD), GRDMAT(25,MCRD)
      DIMENSION HESOLD(MXRCRD,MXRCRD), HESNEW(MCRD,MCRD)
      DIMENSION TMPVEC(MCRD), TMPVC2(MCRD)
      DIMENSION TMPMAT(MXRCRD,MXRCRD), TMPMT2(MXRCRD*MXRCRD)
      DIMENSION TMPMT3(MXRCRD,MXRCRD), TMPMT4(MXRCRD*MXRCRD)
      DIMENSION TMPMT5(MXRCRD*MXRCRD), BMTRAN(MXRCRD,MXRCRD)
      DIMENSION RMAT(MCRD,MCRD), GMAT(MCRD,MCRD)
      LOGICAL RESET, REJLST, DELINT
C
C     A new model Hessian is always calculated as a start
C
      CALL BLDHES(MXRCRD,TMPMAT,HESOLD)
C
C     If we use delocalized internals, we have to transform from
C     redundant to non-redundant space.
C
      IF (DELINT) THEN
         CALL DZERO(TMPMAT,MXRCRD*MXRCRD)
         DO 200 I = 1, ICRD
            DO 202 J = 1, IRCRD
               DO 204 K = 1, IRCRD
                  TMPMAT(I,J) = TMPMAT(I,J) + BMTRAN(K,I)*HESOLD(K,J)
 204           CONTINUE
 202        CONTINUE
 200     CONTINUE
         CALL DZERO(HESOLD,MXRCRD*MXRCRD)
         DO 210 I = 1, ICRD
            DO 212 J = 1, ICRD
               DO 214 K = 1, IRCRD
                  HESOLD(I,J) = HESOLD(I,J) + TMPMAT(I,K)*BMTRAN(K,J)
 214           CONTINUE
 212        CONTINUE
 210     CONTINUE
      END IF
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Model Hessian',-1)
         CALL OUTPUT(HESOLD,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)         
      END IF
C
C     Then we do the updating, suppressing output.
C
      CALL UPMULT(MXRCRD,DELTA,STPMAT,GAMMA,GRDMAT,HESOLD,HESNEW,
     &     ICRD,MCRD,TMPVEC,TMPVC2,TMPMAT,
     &     TMPMT2,TMPMT3,TMPMT4,TMPMT5,MAXDIM,RESET,
     &     GNRM,IPRINT,ITYPE,.FALSE.)
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Updated Hessian',-1)
         CALL OUTPUT(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck upmodo */
      SUBROUTINE UPMODO(MXRCRD,DELTA,STPMAT,GAMMA,GRDMAT,HESOLD,HESNEW,
     &     ICRD,IRCRD,MCRD,TMPVEC,TMPVC2,RMAT,GMAT,TMPMAT,TMPMT2,TMPMT3,
     &     BMTRAN,MAXDIM,RESET,REJLST,GNRM,ITYPE,DELINT,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (D0 = 0.0D0)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      DIMENSION DELTA(MXRCRD), STPMAT(25,MCRD)
      DIMENSION GAMMA(MXRCRD), GRDMAT(25,MCRD)
      DIMENSION HESOLD(MXRCRD,MXRCRD), HESNEW(MCRD,MCRD)
      DIMENSION TMPVEC(MCRD), TMPVC2(MCRD)
      DIMENSION TMPMAT(MXRCRD,MXRCRD), TMPMT2(MXRCRD*MXRCRD)
      DIMENSION TMPMT3(MXRCRD*MXRCRD), BMTRAN(MXRCRD,MXRCRD)
      DIMENSION RMAT(MCRD,MCRD), GMAT(MCRD,MCRD)
      LOGICAL RESET, REJLST, SMART, DELINT
C
C      SMART = .FALSE.
      SMART = .TRUE.
C
C     Since RMAT and GMAT are used to store a maximum of 25 vectors,
C     MCRD should be equal to or larger than 25
C
      MXDIM = MAXDIM
      IF (MCRD .LT. 25) MXDIM = MIN(MAXDIM,MCRD)
C
      CALL DZERO(TMPVEC,MCRD)
      CALL DZERO(TMPVC2,MCRD)
      CALL DZERO(RMAT,MCRD*MCRD)
      CALL DZERO(GMAT,MCRD*MCRD)
      LIMIT = 6
C
C     First we have to transfer the last step and gradient difference to
C     STPMAT and GRDMAT respectively.
C
      IF (RESET) THEN
         CALL DZERO(GRDMAT,25*MCRD)
         CALL DZERO(STPMAT,25*MCRD)
         INUM = 1

C
C     If the last step caused a rejected step, we discard all
C     earlier gradients.
C
C      ELSE IF (REJLST .AND. (INUM .GT. 1)) THEN
C         INUM = 1
C 141     CONTINUE
C         IF (STPMAT(INUM,1) .LT. 1.0D10) THEN
C            INUM = INUM + 1
C            GOTO 141
C         END IF
C         DO 188 I = 1, ICRD
C            STPMAT(1,I) = STPMAT(INUM-1,I)
C            GRDMAT(1,I) = GRDMAT(INUM-1,I)
C            STPMAT(2,I) = 0.0D0
C            GRDMAT(2,I) = 0.0D0
C 188     CONTINUE
C         INUM = 2
C         STPMAT(2,1) = 1.1D10
C         IF (IPRINT .GE. IPRDBG) THEN
C            WRITE(LUPRI,*)
C            WRITE(LUPRI,*) 'Removing all ' //
C     &           'displacement but one due to rejected step.'
C         END IF
C
C     After [LIMIT] iterations we discard the first entries.
C
      ELSE IF (STPMAT(LIMIT,1) .GT. 1.0D10) THEN
         DO 10 I = 1, LIMIT - 1
            DO 12 J = 1, ICRD
               STPMAT(I,J) = STPMAT(I+1,J)
               GRDMAT(I,J) = GRDMAT(I+1,J)
 12         CONTINUE
 10      CONTINUE
         STPMAT(LIMIT,1) = 0.0D0
         INUM = LIMIT - 1
      ELSE
         INUM = 1
 14      CONTINUE
         IF (STPMAT(INUM,1) .LT. 1.0D10) THEN
            INUM = INUM + 1
            GOTO 14
         END IF
C
C     We also have to check the dimension of our variational space,
C     if the number of displacement vectors exceeds this number, we
C     have to remove some.
C
         IF (SMART .AND. (INUM .GT. MXDIM)) THEN
            DO 15 I = 1, INUM-1
               DO 16 J = 1, ICRD
                  STPMAT(I,J) = STPMAT(I+1,J)
                  GRDMAT(I,J) = GRDMAT(I+1,J)
 16            CONTINUE
 15         CONTINUE
            INUM = INUM - 1
            IF (IPRINT .GE. IPRDBG) THEN
               WRITE(LUPRI,*)
               WRITE(LUPRI,*) 'Removing one displacement ' //
     &              'due to dimension of variational space.'
            END IF
         END IF
      END IF
C
C     Then we update all vectors and calculate the norm of the
C     last displacement (SNMLST).
C
      SNMLST = 0.0D0
      DO 17 I = 1, ICRD
         STPMAT(INUM,I) = -DELTA(I)
         GRDMAT(INUM,I) = -GAMMA(I)
         SNMLST = SNMLST + DELTA(I)*DELTA(I)
         DO 19 II = 1, INUM-1
            STPMAT(II,I) = STPMAT(II,I)-DELTA(I)
            GRDMAT(II,I) = GRDMAT(II,I)-GAMMA(I)
 19      CONTINUE
 17   CONTINUE
      SNMLST = SQRT(SNMLST)
      STPMAT(INUM+1,1) = 1.1D10
C
C     We start by calculating the normalized displacement vectors, and
C     we scale the gradient differences. Displacements larger than 100
C     times the last step, and displacements with a gradient difference
C     larger than 0.25D0 are marked (later to be removed).
C
      DO 20 II = 1, INUM
         DNRM = D0
         GNRM = D0
         DO 22 I = 1, ICRD
            DNRM = DNRM + STPMAT(II,I)*STPMAT(II,I)
            GNRM = GNRM + GRDMAT(II,I)*GRDMAT(II,I)
 22      CONTINUE
         DNRM = SQRT(DNRM)
         GNRM = SQRT(GNRM)
         IF (SMART .AND. (DNRM .GE. 100D0*SNMLST)) THEN
            IF (IPRINT .GE. IPRDBG) WRITE(LUPRI,*)
     &           'Removing one displacement due to distance'
            STPMAT(II,1) = 1.1D10
         ELSE IF (SMART .AND. (GNRM .GE. 0.25D0)) THEN
            IF (IPRINT .GE. IPRDBG) WRITE(LUPRI,*)
     &           'Removing one displacement due to gradient norm'
            STPMAT(II,1) = 1.1D10
         ELSE
            DO 24 I = 1, ICRD
C               RMAT(II,I) = STPMAT(II,I)/DNRM
C               GMAT(II,I) = GRDMAT(II,I)/DNRM
               RMAT(II,I) = STPMAT(II,I)
               GMAT(II,I) = GRDMAT(II,I)
 24         CONTINUE
         END IF
 20   CONTINUE
      II = 1
 25   CONTINUE
      IF ((STPMAT(II,1) .GE. 1.0D10) .AND. (INUM .GT. 0)) THEN
         DO 27 JJ = II, INUM - 1
            DO 29 J = 1, ICRD
               STPMAT(JJ,J) = STPMAT(JJ + 1,J)
               GRDMAT(JJ,J) = GRDMAT(JJ + 1,J)
               RMAT(JJ,J)   = RMAT(JJ + 1,J)
               GMAT(JJ,J)   = GMAT(JJ + 1,J)
 29         CONTINUE
 27      CONTINUE
         INUM = INUM - 1
         GOTO 25
      ELSE IF (II .LT. INUM-1) THEN
         II = II + 1
         GOTO 25
      END IF
C
C     We check if the last displacement is nearly parallell to an
C     earlier step (dot product larger than 0.75). If that is the case,
C     the older step is removed.
C
      IF (SMART .AND. (INUM .GT. 1)) THEN
         II = 1
 30      CONTINUE
         DOTP = D0
         DO 32 I = 1, ICRD
            DOTP = DOTP + RMAT(II,I)*RMAT(INUM,I)
 32      CONTINUE
         IF (DOTP .GE. 0.75D0) THEN
            DO 34 J = 1, ICRD
               DO 36 JJ = II, INUM - 1
                  RMAT(JJ,J)   = RMAT(JJ+1,J)
                  STPMAT(JJ,J) = STPMAT(JJ+1,J)
                  GMAT(JJ,J)   = RMAT(JJ+1,J)
                  GRDMAT(JJ,J) = STPMAT(JJ+1,J)
 36            CONTINUE
 34         CONTINUE
            STPMAT(INUM,1) = 1.1D10
            INUM = INUM - 1
            IF (IPRINT .GE. IPRDBG) WRITE(LUPRI,*)
     &           'Removing one nearly parallell and older displacement'
            IF (II .LT. INUM) GOTO 30
         ELSE IF (II .LT. INUM-1) THEN
            II = II + 1
            GOTO 30
         END IF
      END IF
C
C     Some output
C
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('Original displacements',-1)
         CALL OUTPUT(STPMAT,1,INUM,1,ICRD,25,MXRCRD,1,LUPRI)
         CALL HEADER('Normalized displacements',-1)
         CALL OUTPUT(RMAT,1,INUM,1,ICRD,MCRD,MCRD,1,LUPRI)
         CALL HEADER('Original gradient vectors',-1)
         CALL OUTPUT(GRDMAT,1,INUM,1,ICRD,25,MXRCRD,1,LUPRI)
         CALL HEADER('Scaled gradient vectors',-1)
         CALL OUTPUT(GMAT,1,INUM,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
C
C     A new model Hessian is always calculated as a start
C
      CALL BLDHES(MXRCRD,TMPMAT,HESNEW)
C
C     If we use delocalized internals, we have to transform from
C     redundant to non-redundant space.
C
      IF (DELINT) THEN
         CALL DZERO(TMPMAT,MXRCRD*MXRCRD)
         DO 200 I = 1, ICRD
            DO 202 J = 1, IRCRD
               DO 204 K = 1, IRCRD
                  TMPMAT(I,J) = TMPMAT(I,J) + BMTRAN(K,I)*HESNEW(K,J)
 204           CONTINUE
 202        CONTINUE
 200     CONTINUE
         CALL DZERO(HESNEW,MXRCRD*MXRCRD)
         DO 210 I = 1, ICRD
            DO 212 J = 1, ICRD
               DO 214 K = 1, IRCRD
                  HESNEW(I,J) = HESNEW(I,J) + TMPMAT(I,K)*BMTRAN(K,J)
 214           CONTINUE
 212        CONTINUE
 210     CONTINUE
      END IF
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Model Hessian',-1)
         CALL OUTPUT(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)         
      END IF
C
C     Then we use the required subroutine to do the updating,
C     suppressing output.
C
      DO 60 K = 1, INUM
         CALL DZERO(DELTA,MXRCRD)
         CALL DZERO(GAMMA,MXRCRD)
         DO 62 I = 1, ICRD
             DELTA(I) = -RMAT(K,I)
             GAMMA(I) = -GMAT(K,I)
 62      CONTINUE
         IF (ITYPE .EQ. 1) THEN
            CALL UPBFGS(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW,
     &           ICRD,MCRD,TMPVEC,TMPMAT,TMPMT2,-1)
         ELSE IF (ITYPE .EQ. 2) THEN
            CALL UPDFP(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW,
     &           ICRD,MCRD,TMPVEC,TMPMAT,TMPMT2,-1)
         ELSE IF (ITYPE .EQ. 3) THEN
            CALL UPPSB(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW,
     &           ICRD,MCRD,TMPVEC,TMPMAT,TMPMT2,-1)
         ELSE IF (ITYPE .EQ. 4) THEN
            CALL UPRNKO(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW,
     &           ICRD,MCRD,TMPVEC,TMPMAT,TMPMT2,-1)
         END IF
         IF ((IPRINT .GE. IPRDBG) .AND. (K .LT. INUM)) THEN
            CALL HEADER('Updating the Hessian',-1)
            CALL OUTPUT(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)         
         END IF
         CALL DZERO(HESOLD,MXRCRD*MXRCRD)
         DO 65 J = 1, ICRD
            DO 67 I = 1, ICRD
               HESOLD(I,J) = HESNEW(I,J)
 67         CONTINUE
 65      CONTINUE
 60   CONTINUE
      IF (IPRINT .GE. IPRMIN) THEN
         CALL HEADER('Updated Hessian',-1)
         CALL OUTPUT(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck mincgh */
C
C     Calls WLKCGH to do some transformations.
C
      SUBROUTINE MINCGH(EGRAD,EHESS,ALLHES,WORK,LWORK)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "symmet.h"
#include "nuclei.h"
#include "abainf.h"
#include "gnrinf.h"
#include "optinf.h"
#include "cbiwlk.h"
      LOGICAL TRU, FAL
      DIMENSION WORK(LWORK)
      DIMENSION EGRAD(MXCOOR)
      DIMENSION EHESS(MXCOOR,MXCOOR), ALLHES(NCRTOT*NCRTOT)
      IF (REDINT .OR. DELINT) RETURN
      CALL QENTER('MINCGH')
      CALL DZERO(EGRAD,MXCOOR)
      CALL DZERO(EHESS,MXCOOR*MXCOOR)
      CALL DZERO(ALLHES,NCRTOT*NCRTOT)
      TRU = .TRUE.
      FAL = .FALSE.
      CALL WLKCGH(EGRAD,ALLHES,WORK(1),WORK(4),WORK(4+3*NCART),
     &            LWORK-(4+3*NCART),DOREPW,TRU,FAL,NCRTOT,
     &            NCRTOT*NCRTOT,0,NCART,FAL,IPRINT)
      JI = 1
      DO 10 I = 1, NCART
         DO 20 J = 1, NCART
            EHESS(J,I) = ALLHES(JI)
            JI = JI + 1
 20      CONTINUE
 10   CONTINUE
      CALL QEXIT('MINCGH')
      RETURN
      END

C  /* Deck cntcrd */
C
C     Counts the coordinates of all symmetries.
C
      SUBROUTINE CNTCRD
#include "implicit.h"
#include "mxcent.h"
#include "priunit.h"
#include "optinf.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "symmet.h"
#include "prkoor.h"
#include "trkoor.h"
#include "cbiwlk.h"
      NCRTOT = 0
      NTMAT = 0
      NCART = NCRREP(0,1)
      NPROJ = NTRREP(0)
      DO 5 I = 0, 7
         DOREPW(I) = .FALSE.
 5    CONTINUE
      IF (IPRINT .GT. 5) CALL HEADER('Coordinates in CNTCRD',-1)
      DO 10 I  = 0, MAXREP
         NPRREP(I) = NTRREP(I)
         NCR = NCRREP(I,1)
         NPR = NTRREP(I)
         DOREPW(I) = (NCR .GT. NPR)
         IF (DOREPW(I)) THEN
            NCRTOT = NCRTOT + NCR
            NTMAT = NTMAT + NCR*NPR
         END IF
         IF (IPRINT .GT. 5) THEN
            WRITE(LUPRI,'(/A,I2)')' Symmetry', I+1
            WRITE(LUPRI,'(A,I5)')' Cartesian coordinates:', NCR
            WRITE(LUPRI,'(A,I5)')' External coordinates :', NPR
            WRITE(LUPRI,'(A,I5)')' Internal coordinates :', NCR-NPR
            WRITE(LUPRI,'(A,L5)')' Dorepw               :', DOREPW(I)
         END IF
 10   CONTINUE
      IF (IPRINT .GT. 5) WRITE(LUPRI,'(/A,I5)')
     &  ' Total number of Cartesian coordinates:', NCRTOT
      RETURN
      END

C  /* Deck gttmat */
C
C     Construct, scale and orthogonalize T matrix.
C
      SUBROUTINE GTTMAT(TMAT,TMPMAT,NCR,NPR,IREP,THRLDP)
#include "implicit.h"
#include "mxcent.h"
#include "priunit.h"
#include "optinf.h"
      DIMENSION TMAT(NCR,NPR), TMPMAT(MXCOOR)
      CALL DZERO(TMAT,NCR*NPR)
      CALL GETTRO(TMAT,TMPMAT,NCR,NPR,'BOTH','TORTHO',
     &                                         'CT',IREP,IPRINT)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('T matrix in GTTMAT',-1)
         CALL OUTPUT(TMAT,1,NCR,1,NPR,NCR,NPR,1,LUPRI)
      END IF
C
C     No scaling is done yet...
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Scaled T matrix in GTTMAT',-1)
         CALL OUTPUT(TMAT,1,NCR,1,NPR,NCR,NPR,1,LUPRI)
      END IF
      NPR1 = NPR
      CALL ORTVEC(0,NPR1,NCR,THRLDP,TMAT,TMPMAT)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Orthogonalized T matrix in GTTMAT',-1)
         CALL OUTPUT(TMAT,1,NCR,1,NPR,NCR,NPR,1,LUPRI)
      END IF
      IF (NPR1 .NE. NPR) THEN
         WRITE(LUPRI,'(//,2(A,I1),A,/A)')
     *      ' Number of trarot vectors reduced from ', NPR,
     *      ' to ', NPR1, ' in ORTVEC called from WLKPRJ.',
     *      ' Program cannot proceed .'
         CALL QUIT('Insufficient number of trarot vectors in WLKPRJ.')
      END IF
      RETURN
      END

C  /* Deck projgh */
      SUBROUTINE PROJGH(EGRAD,EHESS,ALLHES,TMAT,TMPMAT,TMPMT2,PROJOP)
#include "implicit.h"
#include "mxcent.h"
#include "priunit.h"
#include "optinf.h"
#include "cbiwlk.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "symmet.h"
#include "trkoor.h"
#include "prkoor.h"
      DIMENSION EGRAD(MXCOOR), EHESS(MXCOOR,MXCOOR)
      DIMENSION ALLHES(NCRTOT*NCRTOT), TMAT(NTMAT)
      DIMENSION TMPMAT(MXCOOR,MXCOOR), TMPMT2(MXCOOR)
      DIMENSION PROJOP(MXCOOR,MXCOOR)
      IF (IPRINT .GT. 5) CALL TITLER('Output from PROJGH','*',103)
      IHESS = 1
      ITMAT = 1
      DO 10 IREP = 0, MAXREP
         IF (DOREPW(IREP)) THEN
            NCR = NCRREP(IREP,1)
            NPR = NPRREP(IREP)
            IF (IPRINT .GT. 5) THEN
               WRITE(LUPRI,'(1X,A,I2)') 'Symmetry ', IREP+1
               WRITE(LUPRI,'(1X,A,I2)') 'NCR      ', NCR
               WRITE(LUPRI,'(1X,A,I2)') 'NPR      ', NPR
            END IF
            CALL PRJGH1(EGRAD,ALLHES(IHESS),TMAT(ITMAT),PROJOP,
     &                           TMPMAT,TMPMT2,NCR,NPR,IREP,THRLDP)
            IHESS = IHESS + NCR*NCR
            ITMAT = ITMAT + NPR*NCR
         END IF
 10   CONTINUE
      JI = 1
      DO 20 I = 1, NCART
         DO 30 J = 1, NCART
            EHESS(J,I) = ALLHES(JI)
            JI = JI + 1
 30      CONTINUE
 20   CONTINUE            
      RETURN
      END

C  /* Deck prjgh1 */
C
C     Construct projection operator and use it on gradient and
C     Hessian to remove both rotation and translation.
C
      SUBROUTINE PRJGH1(EGRAD,ALLHES,TMAT,PROJOP,TMPMAT,TMPMT2,
     &                                        NCR,NPR,IREP,THRLDP)
#include "implicit.h"
#include "mxcent.h"
#include "priunit.h"
#include "optinf.h"
      DIMENSION EGRAD(NCR), ALLHES(NCR,NCR)
      DIMENSION TMAT(NCR,NPR), PROJOP(NCR,NCR)
      DIMENSION TMPMAT(MXCOOR,MXCOOR), TMPMT2(MXCOOR)
      CALL DZERO(PROJOP,NCR*NCR)
      CALL DZERO(TMPMAT,MXCOOR*MXCOOR)
C
C     Get translation and rotation matrix.
C
      CALL GTTMAT(TMAT,TMPMT2,NCR,NPR,IREP,THRLDP)
C
C     Construct operator.
C
      CALL DUNIT(PROJOP,NCR)
      CALL DGEMM('N','T',NCR,NCR,NPR,-1.D0,
     &           TMAT,NCR,
     &           TMAT,NCR,1.D0,
     &           PROJOP,NCR)
      IF (IPRINT .GT. 5) THEN
         IF (IREP .EQ. 0) THEN
            CALL HEADER('Unprojected gradient in PROJGH',-1)
            CALL OUTPUT(EGRAD,1,1,1,NCR,1,NCR,1,LUPRI)
         END IF
         CALL HEADER('Unprojected Hessian in PROJGH',-1)
         CALL OUTPUT(ALLHES,1,NCR,1,NCR,NCR,NCR,1,LUPRI)
         CALL HEADER('Projection operator in PROJGH',-1)
         CALL OUTPUT(PROJOP,1,NCR,1,NCR,NCR,NCR,1,LUPRI)
      END IF
C
C     Do projection.
C
      IF (IREP .EQ. 0) THEN
         CALL DGEMM('N','N',NCR,1,NCR,1.D0,
     &              PROJOP,NCR,
     &              EGRAD,NCR,0.D0,
     &              TMPMAT,NCR)
         CALL DCOPY(NCR,TMPMAT,1,EGRAD,1)
         IF (IPRINT .GT. 5) THEN
            CALL HEADER('Projected gradient in PROJGH',-1)
            CALL OUTPUT(EGRAD,1,1,1,NCR,1,NCR,1,LUPRI)
         END IF
      END IF
      CALL DGEMM('N','N',NCR,NCR,NCR,1.D0,
     &           PROJOP,NCR,
     &           ALLHES,NCR,0.D0,
     &           TMPMAT,NCR)
      CALL DGEMM('N','N',NCR,NCR,NCR,1.D0,
     &           TMPMAT,NCR,
     &           PROJOP,NCR,0.D0,
     &           ALLHES,NCR)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Projected Hessian in PROJGH',-1)
         CALL OUTPUT(ALLHES,1,NCR,1,NCR,NCR,NCR,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck diahes */
C
C     Diagonalize Hessian
C
      SUBROUTINE DIAHES(MXRCRD,MX2CRD,NCRDHS,EGRAD,EHESS,ALLHES,
     &     TMAT,THRIND,EVEC,EVCTMP,TMPHES,HESPCK,WORK,LWORK)
#include "implicit.h"
#include "mxcent.h"
#include "priunit.h"
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      PARAMETER (D0 = 0.0D0)
#include "cbiwlk.h"
#include "optinf.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "symmet.h"
      DIMENSION EGRAD(MXCOOR), EHESS(MXCOOR,MXCOOR)
      DIMENSION ALLHES(NCRTOT*NCRTOT), TMAT(NTMAT)
      DIMENSION EVEC(MX2CRD,MX2CRD), EVCTMP(MX2CRD*MX2CRD)
      DIMENSION TMPHES(MX2CRD,MX2CRD)
      DIMENSION HESPCK(NCRDHS*NCRDHS)
      DIMENSION WORK(LWORK)
      LOGICAL PRJTRO
      PRJTRO = .TRUE.
C
      CALL DZERO(GRDDIA,MXRCRD)
      CALL WLKEIG(EGRAD,ALLHES,EVAL,EVCTMP,GRDDIA,TMAT,THRLDP,
     &            THRIND,WORK,CNDHES,INDHES,LWORK,IZEROG,NZEROG,
     &            DOREPW,NCRTOT,NCRTOT*NCRTOT,NTMAT,PRJTRO,IPRINT)
      INDTOT = ISUM(MAXREP+1,INDHES(0),1)
      JI = 1
C
C     Diagonal hessian and eigenvalues are copied.
C
      DO 10 I = 1, NCART
         DO 20 J = 1, NCART
            EHESS(J,I) = ALLHES(JI)
            JI = JI + 1
 20      CONTINUE
 10   CONTINUE
      JI = 1
      II = 0
      DO 30 ISYM = 0, MAXREP
         DO 40 I = 1, NCRREP(ISYM,1)
            DO 50 J = 1, NCRREP(ISYM,1)
               EVEC(II+J,II+I) = EVCTMP(JI)
               JI = JI + 1
 50         CONTINUE
 40      CONTINUE
         II = II + NCRREP(ISYM,1)
 30   CONTINUE
C
C     If we're using the rational function method, we modifiy the
C     totally symmetric part of the Hessian.
C
      IF (RATFUN .AND. (.NOT. SADDLE)) THEN
C
C     Eigenvalues and -vectors must be shifted to make room for one more.
C
         NCR = NCRREP(0,1)
         CALL DZERO(TMPHES,MX2CRD*MX2CRD)
         CALL DZERO(EVCTMP,NCRDHS*NCRDHS)
         IF (MAXREP .GT. 0) THEN
            DO 60 I = NCRTOT, NCR, -1
               DO 62 J = 1, NCRTOT
                  EVEC(I+1,J) = EVEC(I,J)
 62            CONTINUE
 60         CONTINUE
            DO 65 I = NCRTOT, NCR, -1
               DO 67 J = 1, NCRTOT
                  EVEC(J,I+1) = EVEC(J,I)
 67            CONTINUE
               EVAL(I+1) = EVAL(I)
 65         CONTINUE
         END IF
         DO 70 J = 1,  NCR
            DO 72 I = 1,  NCR
               EVCTMP(I+(J-1)*NCRDHS) = EVAL(I)*EVEC(J,I)
 72         CONTINUE
 70      CONTINUE
         DO 75 I = 1,  NCR
            DO 77 J = 1,  NCR
               DO 79 K = 1,  NCR
                  TMPHES(I,J) = TMPHES(I,J)
     &                 + EVEC(I,K)*EVCTMP(K+(J-1)*NCRDHS)
 79            CONTINUE
 77         CONTINUE
 75      CONTINUE
         DO 80 I = 1, NCRTOT
            TMPHES(NCRDHS,I) = EGRAD(I)
            TMPHES(I,NCRDHS) = EGRAD(I)
 80      CONTINUE
         IF (IPRINT .GE. IPRDBG) THEN
            CALL HEADER('Augmented Hessian',-1)
            CALL OUTPUT(TMPHES,1,NCRDHS,1,NCRDHS,MX2CRD,
     &           MX2CRD,1,LUPRI)
         END IF
         TMPHES(NCRDHS,NCRDHS) = D0
         CALL DZERO(EVCTMP,NCRDHS*NCRDHS)
         DO 83 I = 1, NCRDHS
            DO 85 J = 1, NCRDHS
               EVCTMP(J+(I-1)*NCRDHS) = TMPHES(I,J)
 85         CONTINUE
 83      CONTINUE
         CALL DZERO(HESPCK,NCRDHS*NCRDHS)
         CALL DSITSP(NCRDHS,EVCTMP,HESPCK)
         CALL DUNIT(EVCTMP,NCRDHS)
         CALL JACO(HESPCK,EVCTMP,NCRDHS,NCRDHS,NCRDHS,
     &        TMPHES(1,1),TMPHES(1,2))
         DO 90 J = 1, NCRDHS
            EVAL(J) = HESPCK(J*(J+1)/2)
            GRDDIA(J) = DDOT(NCRDHS,EGRAD,1,EVCTMP(1+(J-1)*NCRDHS),1)
            DO 92 I = 1, NCRDHS
               EVEC(I,J) = EVCTMP(I+(J-1)*NCRDHS)
 92         CONTINUE
 90      CONTINUE
         IF ((EVAL(1) .LT. -THRIND) .AND. INDTOT .GT. 0)
     &        INDTOT = INDTOT - 1
         DO 95 I = 1, NCRDHS
            IF (ABS(EVAL(I)) .LE. 1.0D-6) EVAL(I) = EVAL(I) + 1.0D5
 95      CONTINUE
C
C     The eigenvalues are sorted
C
         DO 100 I = 1, NCRDHS
            JMIN = I
            EMIN = EVAL(I)
            DO 105 J = (I + 1), NCRDHS
               IF (EVAL(J) .LT. EMIN) THEN
                  EMIN = EVAL(J)
                  JMIN = J
               END IF
 105        CONTINUE
            IF (JMIN .NE. I) THEN
               CALL DSWAP(1,  EVAL  (I),1,EVAL  (JMIN),1)
               CALL DSWAP(MX2CRD,EVEC(1,I),1,EVEC(1,JMIN),1)
               CALL DSWAP(1,GRDDIA(I),1,GRDDIA(JMIN),1)
            END IF
 100     CONTINUE
         DO 120 I = 1, NCRDHS
            IF (ABS(ABS(EVAL(I))-1.0D5) .LT. 1.0D-3)
     &           EVAL(I) = EVAL(I) - 1.0D5
 120     CONTINUE
         IF (IPRINT .GE. IPRDBG) THEN
            WRITE(LUPRI,*) 'Index of Hessian: ',INDTOT
            CALL HEADER('RF-eigenvalues',-1)
            CALL OUTPUT(EVAL,1,1,1,NCRDHS,1,MXRCRD,1,LUPRI)
            CALL HEADER('RF-eigenvectors',-1)
            CALL OUTPUT(EVEC,1,NCRDHS,1,NCRDHS,MX2CRD,
     &           MX2CRD,1,LUPRI)
         END IF
      END IF
      RETURN
      END

c  /* Deck uptrad */
      SUBROUTINE UPTRAD(REJGEO)
C
C     Updates the trust radius and checks if step should be rejected.
C
#include "implicit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "priunit.h"
#include "infinp.h"
#include "gnrinf.h"
#include "optinf.h"
      PARAMETER (DP25=0.25D0)
      LOGICAL REJGEO
      THDE = 1.0D-10
C
C     CC energies are less precise and we have to decrease the
C     energy threshold
C
      IF (DOCCSD) THDE = 1.0D-8
C
      ERGDIF = ENERGY - ERGOLD
      REJGEO = .FALSE.
      REBILD = .FALSE.
#ifdef VAR_LINSCA
cfilip - print always info about the actual and predicted energy (according to Vebjorn's
c        suggestion). 07.11.2008.
      call header('StepInfo:Energy difference to previous geometry:',-1)
      write(lupri,'(/a),/)')
     &'      Energy       Old_Energy     Actual_diff  Predicted_diff'
     & //'  Ratio'
      write(lupri,'(2F20.10,3F15.10)')
     & ENERGY, ERGOLD, ERGDIF, ERGPRD, ERGDIF/ERGPRD
      write(lupri,*)
cfilip - end.
#endif
C
C     The ratio between actual and predicted energy change is calculated,
C     and this ratio is then used to update the trust radius.
C
      IF (ABS(ERGPRD) .GT. THDE) THEN
         RATIO = ERGDIF/ERGPRD 
         IF (IPRINT .GT. 2) THEN
            CALL HEADER('Energy difference to previous geometry:',-1)
            WRITE(LUPRI,'(/A/5X,F15.10,2(A,F15.10))')
     &      '      Actual           /  Predicted       =    Ratio ',
     &            ERGDIF,'  /  ',ERGPRD,'  =  ',RATIO
         END IF
      ELSE
         RATIO = 1.0D0
         IF (IPRINT .GT. 2) THEN
            WRITE(LUPRI,'(3(/A),/5X,1P,2D16.6)')
     &          ' Close to convergence, ratio set to one.',
     &          ' Energy difference to previous geometry:',
     &          ' actual and predicted:', ERGDIF, ERGPRD
         END IF
      END IF
      IF ((RATIO .LE. RTRJMN) .OR. (RATIO .GE. RTRJMX)) THEN
         IF (SADDLE) THEN
            WRITE(LUPRI,*)
     &           'Trust radius squarely decreased due to bad ratio.'
            TRSTRA = TRSTDE*TRSTDE*STPNRM
         ELSE
            WRITE(LUPRI,'(/A,2(/A,F14.8),/A,F14.2)')
     &           ' Step rejected because ratio between',
     &           '        actual energy change :', ERGDIF,
     &           ' and predicted energy change :', ERGPRD,
     &           '                          is :', RATIO
            REJGEO = .TRUE.
         END IF
      ELSE IF (RATIO .LT. RTENBD) THEN
         WRITE(LUPRI,*) 'Trust radius decreased due to bad ratio.'
         TRSTRA = TRSTDE*STPNRM
      ELSE IF (ABS(RATIO-1.0D0) .LE. DP25*(1.0D0-RTENGD)) THEN
         WRITE(LUPRI,*)
     &       'Trust radius squarely increased due to very good ratio.'
         TRSTRA = MAX(TRSTIN*TRSTIN*STPNRM,TRSTRA)
Chjaug99 IF (ISTATE .GT. 1 .AND. TRSTRA .GT. 0.30D0) TRSTRA = 0.30D0
      ELSE IF (RATIO .GE. RTENGD) THEN
         WRITE(LUPRI,*) 'Trust radius increased due to good ratio.'
         TRSTRA = MAX(TRSTIN*STPNRM,TRSTRA)
Chjaug99 IF (ISTATE .GT. 1 .AND. TRSTRA .GT. 0.30D0) TRSTRA = 0.30D0
      ELSE
         WRITE(LUPRI,*) 'Trust radius set equal to norm of step.'
         TRSTRA = STPNRM
Chjaug99 IF (ISTATE .GT. 1 .AND. TRSTRA .GT. 0.30D0) TRSTRA = 0.30D0
      END IF
C
C     For saddle point optimizations we place both an upper and
C     lower bound on the trust radius.
C
      IF (SADDLE) TRSTRA = MAX(0.025D0, MIN(1.0D0,TRSTRA))
C      IF (SADDLE .AND. DELINT .AND.
C     &     ((RATIO .LE. RTRJMN) .OR. (RATIO .GE. RTRJMX)))
C     &     REBILD = .TRUE.
      IF (DELINT .AND. NEWTON .AND. SADDLE) REBILD = .TRUE.
      IF (.NOT. REJGEO)
     &     WRITE(LUPRI,'(A,F10.5)') ' Updated trust radius', TRSTRA
      RETURN
      END

C  /* Deck rfstp */
      SUBROUTINE RFSTP(MX2CRD,NCHESS,NCRD,ICRD,EVEC,STEP,GRAD,
     &     TMPMAT,HESSMT)
C
C     Determines a rational function (RF) step.
C
#include "implicit.h"
#include "mxcent.h"
#include "priunit.h"
#include "optinf.h"
#include "pi.h"
      DIMENSION EVEC(MX2CRD,MX2CRD)
      DIMENSION STEP(NCRD), GRAD(NCRD)
      DIMENSION TMPMAT(MX2CRD*MX2CRD), HESSMT(NCRD,NCRD)
      LOGICAL STPSCL
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DP5 = 0.5D0)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      STPSCL = .TRUE.
C      STPSCL = .FALSE.
      IF (IPRINT .GE. IPRDBG) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*)
     &        'Lowest eigenvalue (level-shift parameter): ',EVAL(1)
         CALL HEADER('Corresponding eigenvector',-1) 
         CALL OUTPUT(EVEC,1,1,1,NCHESS,1,MX2CRD,1,LUPRI)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Scaling factor: ',EVEC(NCHESS,1)
      END IF
      FAC = EVEC(NCHESS,1)
      IF (ABS(FAC) .GT. 1.0D-8) THEN
         FAC = 1.0D0/FAC
      ELSE
         FAC = 1.0D8
      END IF
      DO 30 I = 1, ICRD
         STEP(I) = EVEC(I,1)*FAC
C
C     For angles and dihedral angles we have to avoid step components
C     giving multiples of 2*pi.
C
         IF (REDINT .AND. INTCRD(I,1) .GT. 10)
     &        STEP(I) = MOD(STEP(I),2.0D0*PI)
C
C     If the step is too large, we simply restrict each element
C     to be below the trust radius.
C
         IF ((ABS(STEP(I)) .GT. TRSTRA) .AND. (.NOT. STPSCL))
     &        STEP(I) = SIGN(TRSTRA,STEP(I))
 30   CONTINUE
      STPNRM = SQRT(DDOT(ICRD,STEP,1,STEP,1))
      IF (IPRINT .GE. IPRDBG) THEN
         WRITE(LUPRI,'(/A,1P,D10.2/)')
     &      'RF-Step length:', STPNRM
         IF (REDINT .OR. DELINT) THEN
            CALL HEADER('RF-Step in internal coordinates',-1) 
         ELSE
            CALL HEADER('RF-step',-1) 
         END IF
         CALL OUTPUT(STEP,1,1,1,ICRD,1,NCRD,1,LUPRI)
      END IF
C
C     Alternatively we restrivt the step norm to be equal or less
C     than the trust radius.
C
      IF ((STPNRM .GT. TRSTRA) .AND. STPSCL) THEN
         FAC = TRSTRA/STPNRM
         DO 32 I = 1, ICRD
            STEP(I) = STEP(I)*FAC
 32      CONTINUE
         IF (IPRINT .GE. IPRDBG) THEN
            WRITE(LUPRI,'(/A,1P,D10.2/)')
     &           'Step too long, step scaled by factor:', FAC
            CALL HEADER('Scaled RF-Step',-1) 
            CALL OUTPUT(STEP,1,1,1,ICRD,1,NCRD,1,LUPRI)
         END IF
         STPNRM = SQRT(DDOT(ICRD,STEP,1,STEP,1))
      END IF
      ICNT = 1
      FAC = 1.0D0
 37   CONTINUE
      IF (ICNT .LE. 10) THEN
         CALL DZERO(TMPMAT,MX2CRD)
         SNDTRM = 0.0D0
         DO 40 I = 1, ICRD
            DO 45 J = 1, ICRD
               TMPMAT(I) = TMPMAT(I) + HESSMT(I,J)*STEP(J)
 45         CONTINUE
            SNDTRM = SNDTRM + TMPMAT(I)*STEP(I)
 40      CONTINUE
         ERGPRD = DDOT(ICRD,GRAD,1,STEP,1)
     &        + 0.5D0*SNDTRM
         IF (ERGPRD .GT. 0.0D0) THEN
            DO 50 I = 1, ICRD
               STEP(I) = EVEC(I,1)*FAC
 50         CONTINUE
            FAC = 0.5D0*FAC
            ICNT = ICNT + 1
            GOTO 37
         END IF
      ELSE
         IF (IPRINT .GE. IPRMIN) THEN
            IF (REDINT .OR. DELINT) THEN
               CALL HEADER('RF-Step in internal coordinates',-1) 
            ELSE
               CALL HEADER('RF-step',-1) 
            END IF
            CALL OUTPUT(STEP,1,1,1,ICRD,1,NCRD,1,LUPRI)
         END IF
         DO 55 I = 1, ICRD
            STEP(I) = -GRAD(I)
 55      CONTINUE
         CALL DZERO(TMPMAT,MX2CRD)
         SNDTRM = 0.0D0
         DO 60 I = 1, ICRD
            DO 65 J = 1, ICRD
               TMPMAT(I) = TMPMAT(I) + HESSMT(I,J)*STEP(J)
 65         CONTINUE
            SNDTRM = SNDTRM + TMPMAT(I)*STEP(I)
 60      CONTINUE
         ERGPRD = DDOT(ICRD,GRAD,1,STEP,1)
     &        + 0.5D0*SNDTRM
      END IF
      IF ((ICNT .GT. 1) .AND. (IPRINT .GE. IPRMIN)) THEN
         IF (REDINT .OR. DELINT) THEN
            CALL HEADER('RF-Step in internal coordinates',-1) 
         ELSE
            CALL HEADER('RF-step',-1) 
         END IF
         CALL OUTPUT(STEP,1,1,1,ICRD,1,NCRD,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck prfstp */
      SUBROUTINE PRFSTP(MX2CRD,NCHESS,NCRD,EVEC,STEP,GRAD,
     &     TMPMAT,HESSMT,IPRF)
C
C     Determines a partitioned rational function (RF) step.
C
#include "implicit.h"
#include "mxcent.h"
#include "priunit.h"
#include "optinf.h"
#include "pi.h"
      DIMENSION EVEC(MX2CRD,MX2CRD)
      DIMENSION STEP(NCRD), GRAD(NCRD)
      DIMENSION TMPMAT(MX2CRD*MX2CRD), HESSMT(NCRD,NCRD)
      LOGICAL STPSCL
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DP5 = 0.5D0)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
      STPSCL = .TRUE.
C      STPSCL = .FALSE.
C
C     This subroutine is called three times. First to minimize one
C     partition, then to maximize another partition, and finally to
C     scale the step and predict the energy change. IPRF keeps track
C     of this (IPRF = 1,2,3).
C
      IF (IPRF .EQ. 3) GOTO 277
C
      IMOD = 1
      IF (IPRF .EQ. 2) IMOD = NCHESS
      IF (IPRINT .GE. IPRDBG) THEN
         WRITE(LUPRI,*)
         IF (IPRF .NE. 2) THEN
            WRITE(LUPRI,*)
     &           'Lowest eigenvalue (level-shift parameter): ',EVAL(1)
         ELSE
            WRITE(LUPRI,*)
     &           'Highest eigenvalue (level-shift parameter): ',
     &           EVAL(NCHESS)
         END IF
         CALL HEADER('Corresponding eigenvector',-1)
         CALL OUTPUT(EVEC,1,NCHESS,IMOD,IMOD,MX2CRD,MX2CRD,1,LUPRI)
      END IF
      FAC = EVEC(NCHESS,IMOD)
      IF (ABS(FAC) .GT. 1.0D-10) THEN
         FAC = 1.0D0/FAC
      ELSE IF ((IMOD .EQ. 1) .AND. (EVAL(IMOD+1) .LT. D0)
     &        .AND. (GRADNM .GT. 1.0D-4)) THEN
 10      CONTINUE
         IMOD = IMOD + 1
         IF (EVAL(IMOD) .LT. D0) THEN
            IF (ABS(EVEC(NCHESS,IMOD)) .GT. 1.0D-10) THEN
               FAC = 1.0D0/EVEC(NCHESS,IMOD)
            ELSE
               GOTO 10
            END IF
         ELSE
            IMOD = 1
            FAC = 1.0D10
         END IF
      ELSE
         FAC = 1.0D10
      END IF
      IF (IPRINT .GE. IPRDBG) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Scaling factor: ',FAC
      END IF
C
C     We only need the preliminary RF-step the first two times the
C     subroutine is called, and we return at this point.
C     We also check the factor in use.
C
      IF ((IPRF .EQ. 1) .OR. (IPRF .EQ. 2)) THEN
         IF ((ABS(FAC) .GT. 1.0D4) .OR. (ABS(FAC) .LT. 1.0D-4)) FAC = D0
         DO 25 I = 1, NCHESS-1
            STEP(I) = EVEC(I,IMOD)*FAC
 25      CONTINUE
         RETURN
      END IF
C
C     The second part of the subroutine scales the step and
C     predicts the energy.
C
 277  CONTINUE
C
      DO 30 I = 1, NCRD
C
C     For angles and dihedral angles we have to avoid step components
C     giving multiples of 2*pi.
C
         IF (REDINT .AND. INTCRD(I,1) .GT. 10)
     &        STEP(I) = MOD(STEP(I),2.0D0*PI)
C
C     If the step is too large, we simply restrict each element
C     to be below the trust radius.
C
         IF ((ABS(STEP(I)) .GT. TRSTRA) .AND. (.NOT. STPSCL))
     &        STEP(I) = SIGN(TRSTRA,STEP(I))
 30   CONTINUE
      STPNRM = SQRT(DDOT(NCRD,STEP,1,STEP,1))
      IF (IPRINT .GE. IPRDBG) THEN
         WRITE(LUPRI,'(/A,1P,D10.2/)')
     &      'RF-Step length:', STPNRM
         IF (REDINT .OR. DELINT) THEN
            CALL HEADER('RF-step in internal coordinates',-1)
         ELSE
            CALL HEADER('RF-step',-1)
         END IF
         CALL OUTPUT(STEP,1,1,1,NCRD,1,NCRD,1,LUPRI)
      END IF
C
C     Alternatively we restrivt the step norm to be equal or less
C     than the trust radius.
C
      IF ((STPNRM .GT. TRSTRA) .AND. STPSCL) THEN
         FAC = TRSTRA/STPNRM
         DO 32 I = 1, NCRD
            STEP(I) = STEP(I)*FAC
 32      CONTINUE
         IF (IPRINT .GE. IPRDBG) THEN
            WRITE(LUPRI,'(/A,1P,D10.2/)')
     &           'Step too long, step scaled by factor:', FAC
            CALL HEADER('Scaled RF-Step',-1)
            CALL OUTPUT(STEP,1,1,1,NCRD,1,NCRD,1,LUPRI)
         END IF
      END IF
      ICNT = 1
      FAC = 1.0D0
 37   CONTINUE
      IF (ICNT .LE. 10) THEN
         CALL DZERO(TMPMAT,MX2CRD)
         SNDTRM = 0.0D0
         DO 40 I = 1, NCRD
            DO 45 J = 1, NCRD
               TMPMAT(I) = TMPMAT(I) + HESSMT(I,J)*STEP(J)
 45         CONTINUE
            SNDTRM = SNDTRM + TMPMAT(I)*STEP(I)
 40      CONTINUE
         ERGPRD = DDOT(NCRD,GRAD,1,STEP,1)
     &        + 0.5D0*SNDTRM
      ELSE
         IF (IPRINT .GE. IPRMIN) THEN
            IF (REDINT .OR. DELINT) THEN
               CALL HEADER('RF-Step in internal coordinates',-1)
            ELSE
               CALL HEADER('RF-step',-1)
            END IF
            CALL OUTPUT(STEP,1,1,1,NCRD,1,NCRD,1,LUPRI)
         END IF
         DO 55 I = 1, NCRD
            STEP(I) = -GRAD(I)
 55      CONTINUE
         CALL DZERO(TMPMAT,MX2CRD)
         SNDTRM = 0.0D0
         DO 60 I = 1, NCRD
            DO 65 J = 1, NCRD
               TMPMAT(I) = TMPMAT(I) + HESSMT(I,J)*STEP(J)
 65         CONTINUE
            SNDTRM = SNDTRM + TMPMAT(I)*STEP(I)
 60      CONTINUE
         ERGPRD = DDOT(NCRD,GRAD,1,STEP,1)
     &        + 0.5D0*SNDTRM
      END IF
      IF ((ICNT .GT. 1) .AND. (IPRINT .GE. IPRMIN)) THEN
         IF (REDINT .OR. DELINT) THEN
            CALL HEADER('RF-Step in internal coordinates',-1)
         ELSE
            CALL HEADER('RF-step',-1)
         END IF
         CALL OUTPUT(STEP,1,1,1,NCRD,1,NCRD,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck prfstc */
      SUBROUTINE PRFSTC(MXRCRD,MX2CRD,NCRDHS,EGRAD,EHESS,EVEC,EVCTMP,
     &     TMPMAT,TMPMT2,TMPMT3,TMPMT4,VECMOD)
C
C     Controls saddle point optimization in Cartesian
C     coordinates using the partitioned rational function approach.
C
#include "implicit.h"
#include "mxcent.h"
#include "priunit.h"
#include "optinf.h"
      DIMENSION EGRAD(MXCOOR), EHESS(MXCOOR,MXCOOR)
      DIMENSION EVEC(MX2CRD,MX2CRD), EVCTMP(NCRDHS*NCRDHS)
      DIMENSION TMPMAT(MX2CRD*MX2CRD),TMPMT2(MX2CRD,MX2CRD)
      DIMENSION TMPMT3(MX2CRD*MX2CRD),TMPMT4(MX2CRD,MX2CRD)
      DIMENSION VECMOD(MXCOOR)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DP5 = 0.5D0)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
C     For saddle point optimizations, we can follow a specific eigenvector.
C     Due to the fact that we are separating one mode for maximization,
C     NCRDHS is temporarily reduced by one.
C
      IMODE = NSPMOD
      NCRDHS = NCRDHS-1
      IF (NSPMOD .GT. 0) THEN
         CALL FNDMOD(.FALSE.,MXRCRD,EVEC,TMPMAT,VECMOD,
     &        TMPMT2,TMPMT3,IMODE)
C
C     If the lowest mode has a gradient element of zero, we have to pick
C     another mode for maximization (or we will end up in a minimum!).
C
      ELSE
         IMODE = 1
 50      CONTINUE
         IF ((ABS(GRDDIA(IMODE)) .LT. 1.0D-10) .AND.
     &        (IMODE .LT. NCRDHS)) THEN
            IMODE = IMODE + 1
            GOTO 50
C     
C     If we find no such mode, we just set IMODE = 1, because we must
C     be at a stationary point.
C     
         ELSE IF (ABS(GRDDIA(IMODE)) .LT. 1.0D-10) THEN
            IMODE = 1
         END IF
         IF (IMODE .NE. 1) THEN
            WRITE (LUPRI,'(A,I3,A)') ' INFO: Mode',IMODE,
     &         ' is selected, because lower modes have zero gradient.'
         END IF
      END IF
      IF (IPRINT .GE. IPRMAX) THEN
         WRITE(LUPRI,'(A,I3,A)') ' Mode',IMODE,
     &      ' will be partitioned out and maximized.'
      END IF
C
C     The selected mode is placed at the very end.
C
      CALL DZERO(TMPMT2,MX2CRD*MX2CRD)
      CALL DZERO(TMPMT4,MX2CRD*MX2CRD)
      DO 400 I = 1, NCRDHS
         DO 402 J = 1, IMODE-1
            TMPMT2(I,J) = EVEC(I,J)
 402     CONTINUE
         DO 403 J = IMODE, NCRDHS-1
            TMPMT2(I,J) = EVEC(I,J+1)
 403     CONTINUE
         TMPMT4(1,I) = EVAL(I)
         TMPMT4(2,I) = GRDDIA(I)
         TMPMT2(I,NCRDHS) = EVEC(I,IMODE)
 400  CONTINUE
      TMPVAL = TMPMT4(1,IMODE)
      DO 406 I = IMODE, NCRDHS-1
         TMPMT4(1,I) = TMPMT4(1,I+1)
         TMPMT4(2,I) = TMPMT4(2,I+1)
 406  CONTINUE
      TMPMT4(1,NCRDHS) = TMPVAL
      TMPMT4(2,NCRDHS) = GRDDIA(IMODE)
C
C     We then make the augmented Hessian that will be minimized.
C
      CALL DZERO(EVCTMP,NCRDHS*NCRDHS)
      DO 410 I = 1, NCRDHS-1
         EVCTMP(I+(I-1)*NCRDHS) = TMPMT4(1,I)
         EVCTMP(I+(NCRDHS-1)*NCRDHS) = TMPMT4(2,I)
         EVCTMP(NCRDHS+(I-1)*NCRDHS) = TMPMT4(2,I)
 410  CONTINUE
      EVCTMP(NCRDHS+(NCRDHS-1)*NCRDHS) = D0
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Augmented Hessian',-1)
         CALL OUTPUT(EVCTMP,1,NCRDHS,1,NCRDHS,NCRDHS,
     &        NCRDHS,1,LUPRI)
      END IF
      CALL DZERO(TMPMT3,MX2CRD*MX2CRD)
      CALL DSITSP(NCRDHS,EVCTMP,TMPMT3)
      CALL DUNIT(EVCTMP,NCRDHS)
      CALL JACO(TMPMT3,EVCTMP,NCRDHS,NCRDHS,NCRDHS,
     &     TMPMAT(1),TMPMAT(1+MX2CRD))
      DO 420 J = 1, NCRDHS
         EVAL(J) = TMPMT3(J*(J+1)/2)
         DO 425 I = 1, NCRDHS
            EVEC(I,J) = EVCTMP(I+(J-1)*NCRDHS)
 425     CONTINUE
 420  CONTINUE
C
C     We add 1.0D5 to all eigenvalues that are essentially zero
C     for the sorting.
C
      DO 427 I = 1, NCRDHS
         IF (ABS(EVAL(I)) .LE. 1.0D-8) EVAL(I) = EVAL(I) + 1.0D5
 427  CONTINUE
      DO 430 I = 1, NCRDHS
         JMIN = I
         EMIN = EVAL(I)
         DO 435 J = (I + 1), NCRDHS
            IF (EVAL(J) .LT. EMIN) THEN
               EMIN = EVAL(J)
               JMIN = J
            END IF
 435     CONTINUE
         IF (JMIN .NE. I) THEN
            CALL DSWAP(1,  EVAL  (I),1,EVAL  (JMIN),1)
            CALL DSWAP(MX2CRD,EVEC(1,I),1,EVEC(1,JMIN),1)
C     CALL DSWAP(1,GRDDIA(I),1,GRDDIA(JMIN),1)
         END IF
 430  CONTINUE
      DO 440 I = 1, NCRDHS
         IF (ABS(ABS(EVAL(I))-1.0D5) .LT. 1.0D-3)
     &        EVAL(I) = EVAL(I) - 1.0D5
 440  CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('RF-eigenvalues',-1)
         CALL OUTPUT(EVAL,1,1,1,NCRDHS,1,MXRCRD,1,LUPRI)
         CALL HEADER('RF-eigenvectors',-1)
         CALL OUTPUT(EVEC,1,NCRDHS,1,NCRDHS,MX2CRD,MX2CRD,1,LUPRI)
      END IF
      CALL PRFSTP(MX2CRD,NCRDHS,MXCOOR,EVEC,STPDIA,EGRAD,
     &     TMPMAT,EHESS,1)
C
C     In the case of saddle point optimization, we also need
C     to take care of the second partition and combine the two.
C
      CMPLIM = MAX(TRSTRA*0.67D0, 0.30D0)
      DO 500 I = 1, NCRDHS-1
         IF (ABS(STPDIA(I)) .GT. CMPLIM)
     &        STPDIA(I) = SIGN(CMPLIM,STPDIA(I))
         TMPMT4(3,I) = STPDIA(I)
 500  CONTINUE
      NCRDHS = NCRDHS + 1
C
C     We then make the augmented Hessian that will be maximized.
C
      CALL DZERO(EVCTMP,NCART*NCART)
      EVCTMP(1) = TMPMT4(1,NCRDHS-1)
      EVCTMP(2) = TMPMT4(2,NCRDHS-1)
      EVCTMP(3) = TMPMT4(2,NCRDHS-1)
      EVCTMP(4) = D0
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Augmented Hessian',-1)
         CALL OUTPUT(EVCTMP,1,2,1,2,2,2,1,LUPRI)
      END IF
      CALL DZERO(TMPMT3,MX2CRD*MX2CRD)
      CALL DSITSP(2,EVCTMP,TMPMT3)
      CALL DUNIT(EVCTMP,2)
      CALL JACO(TMPMT3,EVCTMP,2,2,2,TMPMAT(1),TMPMAT(1+MX2CRD))
      DO 510 J = 1, 2
         EVAL(J) = TMPMT3(J*(J+1)/2)
         DO 515 I = 1, 2
            EVEC(I,J) = EVCTMP(I+(J-1)*2)
 515     CONTINUE
 510  CONTINUE
      DO 517 I = 1, 2
         IF (ABS(EVAL(I)) .LE. 1.0D-8) EVAL(I) = EVAL(I) + 1.0D5
 517  CONTINUE
C
C     The eigenvalues are sorted
C
      IF (EVAL(1) .GT. EVAL(2)) THEN
         CALL DSWAP(1,  EVAL  (1),1,EVAL  (2),1)
         CALL DSWAP(MX2CRD,EVEC(1,1),1,EVEC(1,2),1)
      END IF
      DO 520 I = 1, 2
         IF (ABS(ABS(EVAL(I))-1.0D5) .LT. 1.0D-3)
     &        EVAL(I) = EVAL(I) - 1.0D5
 520  CONTINUE
C
      CALL PRFSTP(MX2CRD,2,MXCOOR,EVEC,STPDIA,EGRAD,
     &     TMPMAT,EHESS,2)
      TMPVAL = STPDIA(1)
      IF (ABS(TMPVAL) .GT. CMPLIM)
     &     TMPVAL = SIGN(CMPLIM,TMPVAL)
      CALL DZERO(STPSYM,MXCOOR)
      DO 530 I = 1, NCRDHS-2
         STPSYM(I) = TMPMT4(3,I)
 530  CONTINUE
      STPSYM(NCRDHS-1) = TMPVAL
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Diagonal RF-step',-1)
         CALL OUTPUT(STPSYM,1,1,1,NCRDHS-1,1,MXRCRD,1,LUPRI)
      END IF
C
C     The symmetry step is constructed from the original eigenvectors
C     (of the normal Hessian) and the diagonal RF-step.
C
      CALL DZERO(STPDIA,MXRCRD)
      CALL DZERO(EVEC,MX2CRD*MX2CRD)
      DO 540 I = 1, NCRDHS-1
         DO 545 J = 1, NCRDHS-1
            EVEC(I,J) = TMPMT2(I,J)
 545     CONTINUE
         EVAL(I) = TMPMT4(1,I)
 540  CONTINUE
      DO 550 I = 1, NCRDHS-1
         DO 555 J = 1, NCRDHS-1
            STPDIA(I) = STPDIA(I) + EVEC(I,J)*STPSYM(J)
 555     CONTINUE
 550  CONTINUE
C
C     The final call to RFSTP to do scaling with respect to the
C     trust radius.
C
      CALL PRFSTP(MX2CRD,NCRDHS,MXCOOR,EVEC,STPDIA,EGRAD,
     &     TMPMAT,EHESS,3)
      RETURN
      END

C  /* Deck fndstp */
      SUBROUTINE FNDSTP(MXRCRD,MX2CRD,NCRDHS,EGRAD,EHESS,EVEC,TMPMAT,
     &     EVCTMP,TMPMT2,TMPMT3,TMPMT4,CSTEP,GRDARR,STPARR,ACTIVE,
     &     EMOD,VECMOD)
C
C     This routine calculates the step that should be taken to obtain
C     the next geometry.
C
#include "implicit.h"
#include "mxcent.h"
#include "priunit.h"
#include "optinf.h"
#include "cbiwlk.h"
#include "trkoor.h"
      PARAMETER (D0 = 0.0D0 , DP5 = 0.5D0)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      DIMENSION EGRAD(MXCOOR), EHESS(MXCOOR,MXCOOR)
      DIMENSION EVEC(MX2CRD,MX2CRD)
      DIMENSION TMPMAT(MX2CRD,MX2CRD),EVCTMP(NCART,NCART)
      DIMENSION TMPMT2(MX2CRD,MX2CRD),TMPMT3(MX2CRD,MX2CRD)
      DIMENSION TMPMT4(MX2CRD,MX2CRD)
      DIMENSION CSTEP(MXCOOR), GRDARR(MXRCRD,25), STPARR(MXRCRD,25)
      DIMENSION VECMOD(MXCOOR)
      LOGICAL INSIDE, ACTIVE, DOSCAL
      CALL QENTER('FNDSTP')
      IF (LNSRCH .AND. (.NOT. RATFUN) .AND. (ITRNMR .GT. 0))
     &     CALL LINSRC(NCART,MXCOOR,EGRAD,GRDARR(1,1),CSTEP,
     &     STPARR(1,1),TMPMAT,TMPMT2,ACTIVE,EMOD)
      IF (ACTIVE) THEN
         DO 5 J = 1, NCART
            DO 7 I = 1, KEPTIT
               STPARR(J,I) = STPARR(J,I) - CSTEP(J)
 7          CONTINUE
            IF (.NOT. RATFUN)
     &           GRDDIA(J) = DDOT(NCART,EGRAD,1,EVEC(1,J),1)
 5       CONTINUE
         IF (.NOT. RATFUN) THEN
            IF (IPRINT .GT. 5) THEN
               CALL HEADER('Diagonal interpolated gradient',1)
               CALL OUTPUT(GRDDIA,1,1,1,NCART,1,MXRCRD,1,LUPRI)
            END IF
         END IF
      END IF
      CALL DZERO(STPDIA,NCART)
      CALL DZERO(STPSYM,NCART)
      NPROJ = NTRREP(0)
      NVEC = NCART-NPROJ
      IF (IPRINT .GT. 2) CALL TITLER('Output from FNDSTP','*',103)
C
C
C     First comes the trust region method
C
      GRADNM = SQRT(DDOT(NVEC,GRDDIA,1,GRDDIA,1))
      IF (TRSTRG .OR. (GDIIS .AND. (KEPTIT .LT. 3))) THEN
C
C     We take a copy of the eigenvectors, to get correct
C     matrix dimensions.
C
         DO 10 I = 1, NCART
            DO 20 J = 1, NCART
               EVCTMP(J,I) = EVEC(J,I)
 20         CONTINUE
 10      CONTINUE
C
C     For saddle point optimizations, we construct the image function.
C
         IF (SADDLE) THEN
            IMODE = NSPMOD
C
C     We can follow a specific eigenvector if needed...
C
            IF (NSPMOD .GT. 0) THEN
               CALL FNDMOD(.FALSE.,MXRCRD,EVEC,TMPMAT,VECMOD,
     &              TMPMT2,TMPMT3,IMODE)
            ELSE
               IMODE = 1
            END IF
            IF (IPRINT .GT. 5) THEN
               WRITE(LUPRI,*)
               WRITE(LUPRI,*) 'Making image function by changing ' //
     &              'the sign of mode ',IMODE
               WRITE(LUPRI,*)
            END IF
            CALL MAKIMG(NCART,NVEC,MXCOOR,EVAL,GRDDIA,
     &           STPDIA,IMODE,.FALSE.)
         END IF
C
C     Newton step is calculated.
C
         DO 30 I = 1, NVEC
            STPDIA(I) = -GRDDIA(I)/EVAL(I)
 30      CONTINUE
         STPNRM = SQRT(DDOT(NVEC,STPDIA,1,STPDIA,1))
C
C     If Newton step is larger that trust radius, we take a step
C     to the boundary. If the Hessian index is larger than zero,
C     the level-shifted step will also be employed, provided the
C     Newton step is larger than 0.5D-3. For saddle points we
C     employ the level-shift when the index is different from 1.
C
         IF (((STPNRM .GT. TRSTRA) .AND. (.NOT. NOTRST)) .OR.
     &        ((.NOT. SADDLE) .AND. (INDTOT .GT. 0) .AND.
     &        (STPNRM .GE. 0.5D-3)) .OR. (SADDLE .AND.
     &        (INDTOT .NE. 1))) THEN
            IF (IPRINT .GT. 5) THEN
               WRITE(LUPRI,'(/A,F15.10)')' Norm of Newton step:', STPNRM
               WRITE(LUPRI,'(A,F15.10/)')' Trust radius       :', TRSTRA
            END IF
            INSIDE = .FALSE.
            IF (STPNRM .LT. TRSTRA) INSIDE = .TRUE.
            CALL LSHFT0(NCART,NVEC,EVAL,GRDDIA,STPDIA,
     &           MIN(TRSTRA,STPNRM),RNU,.FALSE.,ZERGRD,INSIDE,IPRINT)
            STPNRM = SQRT(DDOT(NVEC,STPDIA,1,STPDIA,1))
            IF (IPRINT .GT. 5) THEN
               WRITE(LUPRI,'(/A,F15.10)')' Norm, boundary step:', STPNRM
            END IF
         END IF
C
C     For saddle point optimizations, we check that no single
C     step component is too large, we also restore the
C     original function.
C
         IF (SADDLE) THEN
            DOSCAL = (.NOT. NEWTON)
            IF (INITHS .AND. (ITRNMR .EQ. 0)) DOSCAL = .FALSE.
            IF (DOSCAL) THEN
               CMPLIM = MAX(TRSTRA*0.67D0, 0.3D0)
               DO 35 I = 1, NVEC
                  IF (ABS(STPDIA(I)) .GT. CMPLIM)
     &                 STPDIA(I) = SIGN(CMPLIM,STPDIA(I))
 35            CONTINUE
            END IF
            CALL MAKIMG(NCART,NVEC,MXCOOR,EVAL,GRDDIA,
     &           STPDIA,IMODE,.TRUE.)
         END IF
C
C     Energy is predicted, will be used later to update trust radius.
C
         ERGPRD = DDOT(NVEC,GRDDIA,1,STPDIA,1)
     &        + 0.5D0*DV3DOT(NVEC,STPDIA,EVAL,STPDIA)
         WRITE(LUPRI,'(/A,F25.15)') ' Predicted energy change',ERGPRD
C
C     If the predicted energy is positive, it means the Newton step is
C     towards a maximum/saddle point. We then simply reverse the
C     step direction (a bit dirty, but seems to work).
C
         IF ((.NOT. SADDLE) .AND. (ERGPRD .GT. 0.0D0)) THEN
            WRITE(LUPRI,*) 'Reversing step!'
            DO 40 I = 1, NVEC
               STPDIA(I) = -STPDIA(I)
 40         CONTINUE
            ERGPRD = DDOT(NVEC,GRDDIA,1,STPDIA,1)
     &           + 0.5D0*DV3DOT(NVEC,STPDIA,EVAL,STPDIA)
            WRITE(LUPRI,'(A,F25.15)')
     &           ' New pred. energy change',ERGPRD
         END IF
         CALL HEADER('Step in diagonal representation',1)
         CALL OUTPUT(STPDIA,1,1,1,NVEC,1,MXRCRD,1,LUPRI)
         IF (IPRINT .GT. 5) THEN
            CALL HEADER('Eigenvector basis',1)
            CALL OUTPUT(EVCTMP,1,NCART,1,NVEC,NCART,NCART,1,LUPRI)
         END IF
         DO 150 I = 1, NVEC
            CALL DAXPY(NCART,STPDIA(I),EVEC(1,I),1,STPSYM,1)
 150     CONTINUE
C
C     The rational function method
C
      ELSE IF (RATFUN) THEN
         IF (SADDLE) THEN
            CALL PRFSTC(MXRCRD,MX2CRD,NCRDHS,EGRAD,EHESS,EVEC,EVCTMP,
     &           TMPMAT,TMPMT2,TMPMT3,TMPMT4,VECMOD)
         ELSE
            CALL RFSTP(MX2CRD,NCRDHS,MXCOOR,NCART,EVEC,STPDIA,EGRAD,
     &           TMPMAT,EHESS)
         END IF
         IF (IPRINT .GE. IPRMIN) THEN
            WRITE (LUPRI,'(/A,F25.15)')
     &           ' Predicted energy change', ERGPRD
         END IF
         CALL DZERO(STPSYM,MXCOOR)
         DO 200 I = 1, NCART
            IF (ABS(STPDIA(I)) .GE. 1.0D-6) THEN
               STPSYM(I) = STPDIA(I)
            ELSE
               STPSYM(I) = D0
            END IF
 200     CONTINUE
C
C     The Geometrical DIIS method
C
      ELSE IF (GDIIS) THEN
         CALL DZERO(TMPMAT,MX2CRD*MX2CRD)
         CALL DZERO(EVCTMP,NCART*NCART)
C
C     First we have to construct the inverse Hessian.
C
         DO 210 I = 1, NVEC
            DO 212 J = 1, NCART
               TMPMAT(I,J) = EVEC(J,I)/EVAL(I)
 212        CONTINUE
 210     CONTINUE
         DO 215 I = 1, NCART
            DO 217 J = 1, NCART
               DO 219 K = 1, NVEC
                  EVCTMP(I,J) = EVCTMP(I,J) + EVEC(I,K)*TMPMAT(K,J)
 219           CONTINUE
 217        CONTINUE
 215     CONTINUE
C
C     Then the DIIS-step is determined
C
         CALL GDISTP(MXCOOR,NCART,MXRCRD,MX2CRD,STPDIA,EGRAD,EHESS,
     &        EVCTMP,TMPMAT,TMPMT2,TMPMT3,TMPMT4,GRDARR,STPARR)
         IF (IPRINT .GE. IPRMIN) THEN
            WRITE (LUPRI,'(/A,F25.15)')
     &           ' Predicted energy change', ERGPRD
         END IF
         DO 250 I = 1, NCART
            IF (ABS(STPDIA(I)) .GE. 1.0D-6) THEN
               STPSYM(I) = STPDIA(I)
            ELSE
               STPSYM(I) = D0
            END IF
 250     CONTINUE
      END IF
      IF (ACTIVE) THEN
         DO 300 I = 1, NCART
            STPSYM(I) = STPSYM(I) + CSTEP(I)
 300     CONTINUE
         ERGPRD = ERGPRD + (EMOD-ENERGY)
         WRITE (LUPRI,'(/A,F25.15)')
     &        ' Modified energy prediction due to line search', ERGPRD
      END IF
      STPNRM = SQRT(DDOT(NCART,STPSYM,1,STPSYM,1))
      IF (IPRINT .GT. 2) THEN
         CALL HEADER('Cartesian symmetry step vector',-1)
         CALL OUTPUT(STPSYM,1,1,1,NCART,1,NCART,1,LUPRI)
         WRITE(LUPRI,'(/A,F15.10/)') ' Norm of step:', STPNRM
      END IF
      CALL WLKCOR(STPSYM,CSTEP,NCART,MXCOOR,IPRINT)
      CALL QEXIT('FNDSTP')
      RETURN
      END
      
C  /* Deck fndgeo */
      SUBROUTINE FNDGEO(CSTEP,EGRAD,COONEW,COOOLD,EXHER,EXSIR,EXABA,
     &     WORK,LWORK,WRKDLM,IREJ,GEINFO,NEWSTP)
C
C     If the step is acceptable, the geometry is updated
C     and written to file.
C
#include "implicit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "molinp.h"
#include "optinf.h"
#include "gnrinf.h"
#include "priunit.h"
#include "rspprp.h"
#include "esg.h"
#ifdef VAR_LINSCA
#include "esg2.h"
#endif
      DIMENSION CSTEP(MXCOOR), EGRAD(MXCOOR)
      DIMENSION COONEW(3,MXCENT), COOOLD(3,MXCENT), ICRD(3)
      DIMENSION WORK(LWORK)
      DIMENSION GEINFO(0:ITRMAX,6)
      CHARACTER*10 FILENM
      LOGICAL EXHER,EXSIR,EXABA,REJGEO,NEWSTP
      LOGICAL FAILED
      SAVE FAILED, IFAILD
      DATA FAILED, IFAILD /.FALSE.,0/
      CALL QENTER('FNDGEO')
      NEWSTP = .FALSE.
      EXHER  = .FALSE.
      EXSIR  = .FALSE.
      EXABA  = .FALSE.
      RDINPC = .FALSE.
      RDMLIN = .FALSE.
      HRINPC = .FALSE.
      REJGEO = .TRUE.
      ERGOLD = ENERGY
      IJ = 1
      DO 10 J = 1, NUCIND
         DO 20 I = 1, 3
            COOOLD(I,J) = CORD(I,J)
            COONEW(I,J) = CORD(I,J) + CSTEP(IJ)
            IJ = IJ + 1
 20      CONTINUE
 10   CONTINUE
C
C     Here we start a loop to obtain acceptable step, note that REJGEO
C     is initially set TRUE to enter the loop.
C
 50   CONTINUE
      IF ((IREJ .LE. MAXREJ) .AND. REJGEO) THEN
C
C     If geometry stabilization has been requested (experimental feature!!!),
C     all coordinates with a difference less than the limit are set equal.
C
         IF (ISTBLZ .GT. 0) THEN
            IF (IPRINT .GT. 6) THEN
               CALL HEADER('Non-stabilized geometry',-1)
               CALL PRIGEO(COONEW)
            END IF
            DO 33 J = 1, NUCIND - 1
               ICRD(1) = NINT(COONEW(1,J)*10**ISTBLZ)
               ICRD(2) = NINT(COONEW(2,J)*10**ISTBLZ)
               ICRD(3) = NINT(COONEW(3,J)*10**ISTBLZ)
               DO 34 JJ = J, NUCIND
                  DO 35 I = 1, 3
                     IF (ABS(ICRD(I) - NINT(COONEW(I,JJ)*10**ISTBLZ))
     &                    .LE. 1) COONEW(I,JJ) = COONEW(I,J)
 35               CONTINUE
 34            CONTINUE
 33         CONTINUE
            IF (IPRINT .GT. 6) THEN
               CALL HEADER('Stabilized geometry',-1)
               CALL PRIGEO(COONEW)
            END IF
         END IF
C
         IF (IPRINT .GT. 2) THEN
            CALL HEADER('New geometry',-1)
            CALL PRIGEO(COONEW)
         END IF
         CALL WLKMOL(COONEW)
C
C     Write updated files.
C
         CALL PNCMOL(ITRNMR,IPRINT)
C
C     Memory trap has to be reinitialized.
C
         WORK(1) = WRKDLM
C
C     Calculate energy at new geometry, which is compared to predicted
C     energy(change) in UPTRAD.
C     The temporary update of ITRNMR is in case the molecule input is
C     provided in the DALTON input file
C
         ITRNMR = ITRNMR + 1
         CALL GTNRGY(EXHER,EXSIR,EXABA,.FALSE.,WORK,LWORK,WRKDLM)
         ITRNMR = ITRNMR - 1
         IF (IREJ .EQ. 0 .AND. .NOT. ESG) THEN
            CALL UPTRAD(REJGEO)
            IF (.NOT. REJGEO) THEN
               FAILED = .FALSE.
               IFAILD = 0
            END IF
C
C     After the first failure, we are satisfied if the new energy is below
C     the last. No comparison with predicted energy is done.
C
         ELSE IF ((ENERGY .GT. ERGOLD) .AND. .NOT. ESG) THEN
            IF (FAILED .AND. (IFAILD .LE. MAXREJ) .AND.
     &           (ABS(ENERGY-ERGOLD) .LT. 1.0D-5)) THEN
               WRITE(LUPRI,'(/A)') 'Trouble determining step, ' //
     &              'accepting small energy increase.'
               IFAILD = IFAILD + 1
               REJGEO = .FALSE.
            ELSE
               WRITE(LUPRI,'(/A)')
     &              'Step rejected because energy is increasing.'
               WRITE(LUPRI,'(A,F10.5)')' Updated trust radius', TRSTRA
               REJGEO = .TRUE.
            END IF
         ELSE
            WRITE(LUPRI,'(/A)') 'Acceptable step has been found.'
            REJGEO = .FALSE.
            FAILED = .FALSE.
         END IF
         IF (REJGEO) THEN
            IREJ = IREJ + 1
            EXHER  = .FALSE.
            EXSIR  = .FALSE.
            RDINPC = .FALSE.
            RDMLIN = .FALSE.
            HRINPC = .FALSE.
C
C     Line search based on quadratic model
C
            GRADDI = 0.0D0
            DO 60 I = 1, NCART
               GRADDI = GRADDI + GRDDIA(I)*STPDIA(I)
 60         CONTINUE
            GRADDI = GRADDI/STPNRM
C
            IF (IPRINT .GE. 12) THEN
               CALL HEADER('Line search based on quadratic model',-1)
               WRITE(LUPRI,'(A,F12.6)')
     &              ' Energy at last geometry     : ', ERGOLD
               WRITE(LUPRI,'(A,F12.6)')
     &              ' Energy at rejected geometry : ', ENERGY
               WRITE(LUPRI,'(A,F12.6)')
     &              ' Norm of rejected step       : ', STPNRM
               WRITE(LUPRI,'(A,F12.6)')
     &              ' Norm of gradient            : ', GRADNM
               WRITE(LUPRI,'(A,F12.6)')
     &              ' Gradient along step         : ', GRADDI
            END IF
C
C     The minimum for a quadratic model is calculated with the formula
C                    -f'(0)
C     x     =  -------------------
C      min     2*(f(1)-f(0)-f'(0))
C
            FAC = -0.5D0*GRADDI/(ENERGY-ERGOLD-GRADDI)
C
C     If the factor found is very small or very large, we don't trust
C     it. The factor is replaced by "safer" (but rather atbitrary) numbers.
C
            IF (FAC .LT. 0.1D0) FAC = 0.25D0
            IF (FAC .GT. 0.9D0) FAC = 0.75D0
C
C     If trust region method is used, we find a new level-shift,
C     based on a shorter trust radius. NEWSTP indivates this.
C
C            IF (TRSTRG) THEN
C               TRSTRA = FAC*STPNRM
C               DO 90 J = 1, NUCIND
C                  DO 92 I = 1, 3
C                     CORD(I,J) = COOOLD(I,J)
C 92               CONTINUE
C 90            CONTINUE
C               CALL WLKMOL(COOOLD)
C               NEWSTP = .TRUE.
C               IF (IPRINT .GE. 12) THEN
C                  WRITE(LUPRI,'(A,F12.6)')
C     &                 ' Trust radius decreased to   : ', TRSTRA
C               END IF
C               GO TO 9999
C            END IF
C
C     We have to update both steps and their norm.
C
            DO 70 I = 1, IINTCR
               STPINT(I) = STPINT(I)*FAC
 70         CONTINUE
            DO 75 I = 1, NCART
               STPDIA(I) = STPDIA(I)*FAC
               STPSYM(I) = STPSYM(I)*FAC
 75         CONTINUE
            STPNRM = STPNRM*FAC
C
C     We also set the trust radius equal to the new norm
C
            TRSTRA = STPNRM
C
            WRITE(LUPRI,'(A,F12.6)')
     &           ' Minimum for quadratic model : ', FAC
            WRITE(LUPRI,'(A,F12.6)')
     &           ' Norm of new step            : ', STPNRM
C
C     Finally we construct a new geometry based on the factor found
C
            DO 80 J = 1, NUCIND
               DO 85 I = 1, 3
                  COONEW(I,J)=FAC*COONEW(I,J)+(1.0D0-FAC)*COOOLD(I,J)
 85            CONTINUE
 80         CONTINUE
         END IF
         GOTO 50
      ELSE IF (REJGEO) THEN
C
C     Maximum number of allowed rejections reached
C
         GEINFO(ITRNMR,4) = STPNRM
         IF (ITRNMR .LT. ITRMAX) GEINFO(ITRNMR+1,5) = TRSTRA
         GEINFO(ITRNMR,6) = IREJ*1.0D0
C
C     If redundant internal coordinates are used, we try reduzing the
C     number of dihedral angles to one third the original number (high
C     redundancy might cause problems). We only allow this once before
C     we give up (this should be viewed as an emergency solution!).
C
         IF ((REDINT .AND. (.NOT. FAILED)) .AND. (.NOT. CONOPT)) THEN
            FAILED = .TRUE.
            IREJ = -IREJ
            GEINFO(ITRNMR,6) = 0.0D0
            WRITE(LUPRI,*) 'Maximum number of rejected steps (',MAXREJ,
     &           ') reached.'
            WRITE(LUPRI,'(A)') 'No acceptable step found.'
            WRITE(LUPRI,'(/A)')'***** NOTE! *****'
            WRITE(LUPRI,'(A)')'As an emergency solution, ' //
     &           'the number of dihedral angles will be reduced!'
            CALL RREDUN
            TRSTRA = 0.5D0
            GO TO 9999
         ELSE IF (((.NOT. NEWTON) .AND. (.NOT. FAILED)) .AND.
     &           (.NOT. CONOPT)) THEN
            FAILED = .TRUE.
            IREJ = -IREJ
            GEINFO(ITRNMR,6) = 0.0D0
            WRITE(LUPRI,*) 'Maximum number of rejected steps (',MAXREJ,
     &           ') reached.'
            WRITE(LUPRI,'(A)') 'No acceptable step found.'
            WRITE(LUPRI,'(/A)')'***** NOTE! *****'
            WRITE(LUPRI,'(A)')'As a last resort, ' //
     &           'the Hessian is initialized to unity!'
            EVLINI = 1.0D0
            TRSTRA = 0.5D0
            GO TO 9999
C
C     Otherwise we give up...
C
         ELSE
            CALL PRIINF(GEINFO,WORK(2),LWORK)
            WRITE(LUPRI,*) 'Maximum number of rejected steps (',MAXREJ,
     &        ') reached.'
            WRITE(LUPRI,'(A)') 'No acceptable step found. Aborting.'
            CALL QUIT('*** FNDGEO *** No acceptable step found.')
         END IF
      END IF
 9999 CALL QEXIT('FNDGEO')
      RETURN
      END

C  /* Deck dobrk */
      SUBROUTINE DOBRK(MXRCRD,MX2CRD,EXHER,EXSIR,EXABA,GEINFO,EVEC,
     &                                        CSTEP,TMPMAT,NONT)
C
C     This procedure does the actual breaking of symmetry.
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "molinp.h"
#include "gnrinf.h"
#include "optinf.h"
#include "cbiwlk.h"
#include "cbirea.h"
#include "nuclei.h"
#include "symmet.h"
#include "trkoor.h"
#include "huckel.h"
      LOGICAL EXHER, EXSIR, EXABA, BRKALS, AUTOSY, NOSYM, DOCART, DOOWN,
     &        ADDSYM
      CHARACTER FILENM*10, TMPLN*80, TMPTXT*11, BSNM*80
      CHARACTER*1 CRT, ID3, KASYM(3,3)
      DIMENSION GEINFO(0:ITRMAX,6), NONT(MXCENT), JCO(MXAQN,MXCENT)
      DIMENSION EVEC(MX2CRD,MX2CRD), CSTEP(MXCOOR), IQM(MXCENT)
      DIMENSION TMPMAT(NCRTOT*NCRTOT), Q(MXCENT)
      DIMENSION RADIUS_PCM(MXCENT), ALPHA_PCM(MXCENT)
#include "ibtfun.h"
      CALL DZERO(STPDIA,MXRCRD)
      CALL DZERO(CSTEP,MXCOOR)
      ISYMBR = MAXREP
      START = .TRUE.
      BRKALS = .FALSE.
      CALL IZERO(NUCNUM, MXCENT*8)
C
C     We find the highest symmetry number with a non-zero index.
C     Some test runs indicated this as the most effective way to
C     decrease the energy.
C
 10   CONTINUE
      IF ((INDHES(ISYMBR) .LT. 1) .OR. (.NOT. DOREPW(ISYMBR))) THEN
         ISYMBR = ISYMBR - 1
         IF (ISYMBR .LT. 0) THEN
C
C     If breaking the symmetry causes problems, we remove all symmetry
C
            IF (MAXOPR .GT. 0) THEN
               BRKALS = .TRUE.
               GOTO 77
C
C     If all symmetry has been removed and we're still having trouble,
C     it's time to give up!
C
            ELSE
               CALL QUIT
     &         ('*** DOBRK *** Breaking of symmetry was unsuccesful.')
            END IF
         END IF
         GOTO 10
      END IF
      WRITE(LUPRI,*)
      WRITE(LUPRI,*)'***** NOTE! *****'
      WRITE(LUPRI,*)
     &     'Due to non-zero index of total Hessian, symmetry #',ISYMBR
      WRITE(LUPRI,*)'has to be broken to minimize energy!'
      WRITE(LUPRI,*)
C
C     The eigenvectors of the symmetry to be broken, are copied
C     to TMPMAT.
C
      II = 0
      DO 12 I = 0, ISYMBR - 1
         II = II + NCRREP(I,1)
 12   CONTINUE
      JI = 1
      DO 15 I = 1, NCRREP(ISYMBR,1)
         DO 17 J = 1, NCRREP(ISYMBR,1)
            TMPMAT(JI) = EVEC(II+J,II+I)
            JI = JI + 1
 17      CONTINUE
 15   CONTINUE
      NCR = NCRREP(ISYMBR,1)
      NVC = NCR - NTRREP(ISYMBR)
      DO 19 I = 1, NVC
         STPDIA(I) = 1.0D0/SQRT(NVC*1.0D0)
 19   CONTINUE
      IF (IPRINT .GT. 2) THEN
         CALL HEADER('Geometry before break of symmetry',1)
         CALL PRIGEO(CORD)
         CALL HEADER('Step in diagonal representation',1)
         CALL OUTPUT(STPDIA,1,1,1,NVC,1,MXRCRD,1,LUPRI)
      END IF
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Eigenvector basis',1)
         CALL OUTPUT(TMPMAT,1,NCR,1,NVC,NCRTOT,NCRTOT,1,LUPRI)
      END IF
C
C     We determine the cartesian symmetry-breaking step.
C
      DO 30 I = 1, NVC
         CALL DAXPY(NCR,STPDIA(I),TMPMAT((I-1)*NCR+1),1,CSTEP,1)
 30   CONTINUE
      CALL DZERO(TMPMAT, NCRTOT*NCRTOT)
C
C     Scaling the symmetry-breaking step
C
      DO 32 I = 1, NCR
         TMPMAT(I) = CSTEP(I)*1.00D0
 32   CONTINUE
      STPNRM = SQRT(DDOT(NCR,TMPMAT,1,TMPMAT,1))
      CALL DZERO(CSTEP,MXCOOR)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Cartesian symmetry step vector',1)
         CALL OUTPUT(TMPMAT,1,1,1,NCR,1,NCR,1,LUPRI)
         WRITE(LUPRI,'(/A,F15.10/)') ' Norm of step:', STPNRM
      END IF
      DO 40 IATOM = 1, NUCIND
         DO 42 ICOOR = 1, 3
            ICCOOR = 3*(IATOM - 1) + ICOOR
            ISCOOR = IPTCNT(ICCOOR,ISYMBR,1)
            IF (ISCOOR .GT. 0) THEN
               CSTEP(ICCOOR)=TMPMAT(ICCOOR)/SQRT(FMULT(ISTBNU(IATOM)))
            END IF
 42      CONTINUE
 40   CONTINUE
C
C     Occasionally the symmetry breaking step is zero, and another
C     symmetry must be chosen for breaking.
C
      TMP = 0.0D0
      DO 31 I = 1, NCR
         TMP = TMP + CSTEP(I)*CSTEP(I)
 31   CONTINUE
      IF (TMP .LE. ZERGRD) THEN
         WRITE(LUPRI,*)
     &          'Zero step vector found, breaking another symmetry'
         ISYMBR = ISYMBR - 1
         GOTO 10
      END IF
C
C     The new geometry is calculated.
C
 77   CONTINUE
      IJ = 1
      DO 45 J = 1, NUCIND
         DO 47 I = 1, 3
            CORD(I,J) = CORD(I,J) + CSTEP(IJ)
            IJ = IJ + 1
 47      CONTINUE
 45   CONTINUE
      IF (IPRINT .GT. 2) THEN
         CALL HEADER('Cartesian step vector in non-symmetry basis',1)
         CALL PRIGEO(CSTEP)
         CALL HEADER('Geometry after break of symmetry',1)
         CALL PRIGEO(CORD)
      END IF
      CALL WLKMOL(CORD)
C
C     We remove the number of symmetry operations, so that new symmetry
C     will be added in the next iteration. If we have broken the symmetry
C     earlier but ends up with the same energy, we're having trouble.
C     This is brutally resolved by removing all symmetry, that is the
C     molecule has to be minimized withim the C1 point group.
C
      TMPLN = MLINE(1)
      BASIS = (TMPLN(1:5) .EQ. 'BASIS')
      TMPLN = MLINE(NMLAU)
      II = 0
      DO 49 J = 0, ITRNMR-1
         IF ((GEINFO(J,3) .LT. -0.5D0) .AND.
     &        (ABS(ENERGY-GEINFO(J,1)) .LT. THRERG*10)) II=II+1
 49   CONTINUE
      IF (II .GT. 0) BRKALS = .TRUE.
      CALL UPCASE(TMPLN)
      IF (INDEX(TMPLN,'ATO') .NE. 0) THEN
         CALL LINE4(TMPLN,NONTYP,NSYMOP,CRT,KCHARG,THRS,ADDSYM,KASYM,
     &              ID3,DOCART,DOOWN)
         AUTOSY = .TRUE.
         NOSYM = .FALSE.
         IF (BRKALS) NOSYM = .TRUE.
         ID3 = ' '
         CALL LINE4W(TMPLN,NONTYP,NSYMOP,KCHARG,THRS,AUTOSY,NOSYM,
     &               KASYM,ID3,DOCART,DOOWN)
      ELSE
         READ(TMPLN(1:5),'(I5)') NONTYP
         WRITE(TMPLN(10:20), '(A11)') '           '
         IF (BRKALS) WRITE(TMPLN(10:20), '(A11)') '0          '
      END IF
      MLINE(NMLAU) = TMPLN
C
      IATOM = 1
      I = NCLINE(IATOM)-1
      DO 50 J = 1, NONTYP
         CALL UPCASE(MLINE(I))
         IF (INDEX(MLINE(I),'CHA') .NE. 0) THEN
            CALL LINE5R(MLINE(I),Q(J),NONT(J),MBSI,IQM,JCO(1,J),MXAQN,
     &                 BASIS,ATOMBA,LMULBS,BSNM,
     &                 RADIUS_PCM(J), ALPHA_PCM(J))
         ELSE
            READ(MLINE(I),'(BN,6X,F4.1,I5)') Q, NONT(J)
         END IF
         IATOM = IATOM + NONT(J)
         I = NCLINE(IATOM)-1
 50   CONTINUE
C
C     We run over all atom types and all symmetry independent centres and
C     and expand them to all atoms. The molecule input is modified
C     according to this.
C
      IATOM = 1
      DO 60 ITYP = 1, NONTYP
         IXTRA = 0
         DO 70 I = 1, NONT(ITYP)
            MULCNT = ISTBNU(IATOM)
            DO 80 ISYMOP = 1, MAXOPR
               IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
                  COOX=PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,IATOM)
                  COOY=PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,IATOM)
                  COOZ=PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,IATOM)
                  INSLIN = NCLINE(IATOM)
                  DO 90 J = NMLINE, INSLIN, -1
                     MLINE(J+1) = MLINE(J)
 90               CONTINUE
                  DO 95 J = (IATOM+1), NUCIND
                     NCLINE(J) = NCLINE(J) + 1
 95               CONTINUE
                  TMPLN = MLINE(INSLIN + 1)
                  IPOS = INDEX(TMPLN,'Isotope=')
                  IF (IPOS .EQ. 0) THEN
                     TMPTXT = '           '
                  ELSE
                     TMPTXT = TMPLN(IPOS:IPOS+11)
                  END IF
                  IF (ABS(COOX).ge.100.0D0 .OR.
     &                ABS(COOY).ge.100.0D0 .OR.
     &                ABS(COOZ).ge.100.0D0) THEN
                     WRITE(TMPLN(5:76),'(3F20.10,1X,A11)')
     &                    COOX,COOY,COOZ,TMPTXT
                  ELSE
                     WRITE(TMPLN(5:76),'(3F20.15,1X,A11)')
     &                    COOX,COOY,COOZ,TMPTXT
                  END IF
                  MLINE(INSLIN + 1) = TMPLN
                  NMLINE = NMLINE + 1
                  IXTRA = IXTRA + 1
               END IF
 80         CONTINUE
            IATOM = IATOM + 1
 70      CONTINUE
         TMPLN = MLINE(NCLINE(IATOM - NONT(ITYP)) - 1)
         IF (INDEX(TMPLN,'CHA') .NE. 0) THEN
            CALL LINE5W(TMPLN,Q(ITYP),NONT(ITYP) + IXTRA,MBSI,BASIS,
     &           ATOMBA,LMULBS,BSNM,IQM(ITYP),JCO(1,ITYP),MXAQN,
     &           RADIUS_PCM(ITYP), ALPHA_PCM(ITYP))
         ELSE
            WRITE(TMPLN(11:15),'(I5)') NONT(ITYP) + IXTRA
         END IF
         MLINE(NCLINE(IATOM - NONT(ITYP)) - 1) = TMPLN
 60   CONTINUE
C
C     Write updated geometry to files.
C
      CALL PNCMOL(ITRNMR,IPRINT)
C
C     To mark that symmetry was broken, the index of the hessian is
C     given a negative sign. This is interpreted in PRIINF.
C
      GEINFO(ITRNMR,3) = -ABS(GEINFO(ITRNMR,3))
C
C     Several variables and arrays has to be modified/reset to be able
C     to continue calculation with new symmetry.
C
      BRKSYM = .FALSE.
      ITRBRK = ITRNMR
      INDOLD = INDTOT
      GECONV = .FALSE.
      EXHER  = .FALSE.
      EXSIR  = .FALSE.
      EXABA  = .FALSE.
      RDINPC = .FALSE.
      RDMLIN = .FALSE.
      HRINPC = .FALSE.
      NEWSYM = .TRUE.
      DOHUCKEL = .TRUE.
      KEEPHE = .FALSE.
      RSTARR = .TRUE.
      ERGOLD =  ENERGY
      CALL IZERO(NUCNUM, MXCENT*8)
      CALL IZERO(NCRREP, 16)
      CALL IZERO(IPTCNT, MXCENT*48)
      CALL IZERO(NAXREP, 16)
      CALL IZERO(INDHES, 8)
      RETURN
      END

C  /* Deck visulz */
      SUBROUTINE VISULZ(WORK,LWORK,WRKDLM)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "gnrinf.h"
#include "inftap.h"
      DIMENSION WORK(LWORK)
C
C     HERINP has to be run to process geometry input.
C
C      IPRUSR_SAV = IPRUSR
C      IPRUSR = -2
C      CALL HERINP(WORK,LWORK)
C      CALL GPCLOSE(LUONEL,'KEEP')
C      IPRUSR = IPRUSR_SAV
C      WORK(1) = WRKDLM
C
C     Make VRML-file of geometry
C
      KATARR = 2
      KEVEC  = KATARR + 8*MXCENT
      KEVC1  = KEVEC  + MXCOOR*MXCOOR
      KEVC2  = KEVC1  + MXCOOR
      KWRK   = KEVC2  + MXCOOR
      LWRK = LWORK  - KWRK + 1
      IF (KWRK .GT. LWORK) CALL STOPIT('VISULZ',' ',KWRK,LWORK)
      CALL MKVRML(.FALSE.,WORK(KATARR),MXCOOR,WORK(KEVEC),
     &     WORK(KEVC1),WORK(KEVC2))
      RETURN
      END

C  /* Deck inipre */
      SUBROUTINE INIPRE(WORK,LWORK,WRKDLM)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "optinf.h"
#include "cbirea.h"
#include "molinp.h"
#include "abainf.h"
#include "gnrinf.h"
#include "inftap.h"
      CHARACTER TMPLIN*80
      DIMENSION WORK(LWORK)
C
      WORK(1) = WRKDLM
      KFREE = 2
      LFREE = LWORK - KFREE + 1
C
C     When preoptimization is scheduled, we have to change the basis set.
C     READIN has to be run before we can do the first change.
C
C        Third parameter .FALSE. in CALL READIN
C        indicates that LUONEL will NOT be written.
C
C     CALL READIN(WORK(KFREE),LFREE,.FALSE.)
Chjaaj Juni 2009: READIN has always already been called in OPTMIN
C                 before coming here ...
C
C     Basis set library has to be used for preoptimization to work
C
      IF (BASIS) THEN
         TMPLIN = MLINE(2)
         MLINE(2) = PREBTX(1)
         PREBTX(1) = TMPLIN
         ITRNMR = ITRNMR - 1
         CALL PNCMOL(ITRNMR,IPRINT)
         ITRNMR = ITRNMR + 1
      ELSE
         DOPRE = .FALSE.
         WRITE(LUPRI,'(/A/A/)')
     &      '*** WARNING! *** '//
     &      'Preoptimization can only be done when ' //
     &        'basis set library is used.',
     &      '*** WARNING! *** '//
     &      'Preoptimization has been turned off!'
      END IF
      RDINPC = .FALSE.
      RDMLIN = .FALSE.
      IPRE = IPRE + 1
C
C     Make VRML-file of initial geometry if requested
C
      IF (VRML) THEN
         KATARR = KFREE
         KEVEC  = KATARR + 8*MXCENT
         KEVC1  = KEVEC  + MXCOOR*MXCOOR
         KEVC2  = KEVC1  + MXCOOR
         KWRK   = KEVC2  + MXCOOR
         LWRK = LWORK  - KWRK + 1
         IF (KWRK .GT. LWORK) CALL STOPIT('INIPRE',' ',KWRK,LWORK)
         CALL MKVRML(.FALSE.,WORK(KATARR),MXCOOR,WORK(KEVEC),
     &        WORK(KEVC1),WORK(KEVC2))

      END IF
      RETURN
      END

C  /* Deck endpre */
      SUBROUTINE ENDPRE(EXHER,EXSIR,EXABA)
C
C     This procedure ends preoptimization and starts optimization with
C     the "main" basis set.
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "molinp.h"
#include "gnrinf.h"
#include "optinf.h"
#include "cbiwlk.h"
#include "cbirea.h"
#include "nuclei.h"
#include "symmet.h"
#include "trkoor.h"
      CHARACTER TMPLIN*80
      CHARACTER*1 KASYM(3,3), ID3, CRT
      LOGICAL EXHER, EXSIR, EXABA, DOCART, DOOWN, AUTOSY, NOSYM, ADDSYM
#include "ibtfun.h"
      TMPLIN = MLINE(2)
      IF (IPRE .LT. NUMPRE) THEN
         MLINE(2) = PREBTX(IPRE+1)
         PREBTX(IPRE+1) = PREBTX(IPRE)
         PREBTX(IPRE) = TMPLIN
      ELSE
         MLINE(2) = PREBTX(IPRE)
         PREBTX(IPRE) = TMPLIN
         FINPRE = .TRUE.
      END IF
C
C     If all symmetry has been removed during preoptimization, we try
C     to detect it again after each preoptimization phase.
C
      TMPLIN = MLINE(NMLAU)
      IF ((MAXOPR .EQ. 0) .AND. (ITRBRK .GE. 0)) THEN
         IF (INDEX(TMPLIN,'ATO') .NE. 0) THEN
            CALL LINE4(TMPLIN,NONTYP,NSYMOP,CRT,KCHARG,THRS,ADDSYM,
     &                 KASYM,ID3,DOCART,DOOWN)
            AUTOSY = .TRUE.
            NOSYM = .FALSE.
            ID3 = ' '
            CALL LINE4W(TMPLIN,NONTYP,NSYMOP,KCHARG,THRS,AUTOSY,NOSYM,
     &                  KASYM,ID3,DOCART,DOOWN)
         ELSE
            WRITE(TMPLIN(10:20), '(A11)') '           '
         END IF
         MLINE(NMLAU) = TMPLIN
         NWSYMM = .TRUE.
      END IF
      CALL PNCMOL(ITRNMR,IPRINT)
      IPRE = IPRE + 1
      GECONV = .FALSE.
      EXHER  = .FALSE.
      EXSIR  = .FALSE.
      EXABA  = .FALSE.
      RDINPC = .FALSE.
      RDMLIN = .FALSE.
      HRINPC = .FALSE.
      NEWSYM = .TRUE.
      NEWBAS = .TRUE.
      KEEPHE = .TRUE.
      ERGOLD =  ENERGY
      CALL IZERO(NUCNUM, MXCENT*8)
      CALL IZERO(NCRREP, 16)
      CALL IZERO(IPTCNT, MXCENT*48)
      CALL IZERO(NAXREP, 16)
      CALL IZERO(INDHES, 8)
      RETURN
      END

C  /* Deck mkscvc */
      SUBROUTINE MKSCVC(SCLVEC)
C
C     Makes scaling vector that's necessary for manipulating
C     gradient and step vectors.
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "optinf.h"
#include "nuclei.h"
#include "symmet.h"
      DIMENSION SCLVEC(MXCOOR)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
C
C     We initialize the scaling vector (the inverse of the
C     normalization vector in WLKCGH).
C
      CALL DZERO(SCLVEC,MXCOOR)
      DO 1 IREP = 0, MAXREP
         DO 2 ICENT = 1, NUCIND
#include "scalar.h"
            DO 3 ICOOR = 1, 3
               ISCOOR = IPTCNT(3*(ICENT-1)+ICOOR,IREP,1)
               IF (ISCOOR .GT. 0) THEN
                  SCLVEC(ISCOOR) = SQRT(FMULT(ISTBNU(ICENT)))
               END IF
 3          CONTINUE
 2       CONTINUE
 1    CONTINUE
      IF (IPRINT .GE. IPRMED) THEN
         CALL TITLER('Output from MKSCVC','*',103)
         CALL HEADER('Scaling vector',-1)
         CALL OUTPUT(SCLVEC,1,1,1,NCRTOT,1,MXCOOR,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck minend */
      LOGICAL FUNCTION MINEND(MXRCRD,SCLVEC,BMTRAN,TMPVC1,TMPVC2)
C
C     Determines if the end of the optimization has been reached.
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "taymol.h"
#include "gnrinf.h"
#include "optinf.h"
      DIMENSION SCLVEC(MXCOOR), BMTRAN(MXRCRD,MXRCRD)
      DIMENSION TMPVC1(MXRCRD), TMPVC2(MXRCRD)
      LOGICAL CNVERG, CNVGRD, CNVSTP, INDXOK, CNVRGD(1:4)
      CHARACTER*4 LOGTXT
      CNVERG = .FALSE.
      CNVGRD = .FALSE.
      CNVSTP = .FALSE.
      ICOOR = NCART
      IF (REDINT .OR. DELINT) ICOOR = IINTCR
      IF (CONOPT) GNRM = SQRT(DDOT(IINTCR,GRDINT,1,GRDINT,1))
C
C     First, the convergence criteria of Baker.
C     To make the convergence criteria more fair, we transform
C     the delocalized internals to the primitive space.
C
      IF (BAKER) THEN
         IF (DELINT) THEN
            CALL DZERO(TMPVC1,MXRCRD)
            CALL DZERO(TMPVC2,MXRCRD)
            DO 10 I = 1, IREDIC
               DO 20 J = 1, IINTCR
                  TMPVC1(I) = TMPVC1(I) + BMTRAN(I,J)*GRDINT(J)
                  TMPVC2(I) = TMPVC2(I) + BMTRAN(I,J)*STPINT(J)
 20            CONTINUE
 10         CONTINUE
            CALL MAXELM(TMPVC1,IREDIC,SCLVEC,0,GRDMAX)
            CALL MAXELM(TMPVC2,IREDIC,SCLVEC,0,STPMAX)
         ELSE IF (REDINT) THEN
            CALL MAXELM(GRDINT,IINTCR,SCLVEC,0,GRDMAX)
            CALL MAXELM(STPINT,IINTCR,SCLVEC,0,STPMAX)
         ELSE
            CALL MAXELM(GRDMOL,ICOOR,SCLVEC,2,GRDMAX)
            CALL MAXELM(STPSYM,ICOOR,SCLVEC,1,STPMAX)
         END IF
         EDIFF  = ABS(ERGOLD-ENERGY)
         INDXOK = (INDTOT .EQ. 0)
         IF (SADDLE) INDXOK = (INDTOT .EQ. 1)
C
         CNVGRD = (GRDMAX .LE. 3.0D-4)
         CNVERG = (EDIFF  .LE. 1.0D-6)
         IF ((ITRNMR .LT. 1) .OR. (ITRBRK .EQ. (ITRNMR-1)))
     &          CNVERG = .FALSE.
         CNVSTP = (STPMAX .LE. 3.0D-4)
         MINEND = (INDXOK .AND. CNVGRD) .AND. (CNVERG .OR. CNVSTP)
         IF (NOBRKS .OR. (.NOT. NEWTON) .OR. SADDLE)
     &        MINEND =  (CNVGRD .AND. (CNVERG .OR. CNVSTP))
C
C     Next, the default (old) convergence check
C
      ELSE IF (.NOT. NATNRM) THEN
         ICONV = 0
         EDIFF = ABS(ERGOLD-ENERGY)
         IF (EDIFF .LE. THRERG) THEN
            CNVERG = .TRUE.
            ICONV = ICONV + 1
         END IF
         IF (CONOPT) THEN
            IF (GNRM .LE. GRDTHR) THEN
               CNVGRD = .TRUE.
               ICONV = ICONV + 1
            END IF
         ELSE
            IF (GRADNM .LE. GRDTHR) THEN
               CNVGRD = .TRUE.
               ICONV = ICONV + 1
            END IF
         END IF
         IF (STPNRM .LE. THRSTP) THEN
            CNVSTP = .TRUE.
            ICONV = ICONV + 1
         END IF
         INDXOK = (INDTOT .EQ. 0)
         IF (SADDLE) INDXOK = (INDTOT .EQ. 1)
         MINEND = INDXOK .AND. (ICONV .GE. ICONDI)
         IF (NOBRKS) MINEND = (ICONV .GE. ICONDI)
C
C     Finally, the new convergence scheme (NATNRM)
C
      ELSE
         CALL MAXELM(GRDMOL,ICOOR,SCLVEC,2,GRDMAX)
         CALL MAXELM(STPSYM,ICOOR,SCLVEC,1,STPMAX)
         CNVRGD(1) = .FALSE.
         CNVRGD(2) = .FALSE.
         CNVRGD(3) = .FALSE.
         CNVRGD(4) = .FALSE.
         IF (GRADNM/SQRT(1.D0*ICRTCR) .LT. GRDTHR) CNVRGD(1) = .TRUE.
         IF (GRDMAX .LT. THGRMX) CNVRGD(2)                   = .TRUE.
         IF (STPNRM/SQRT(1.D0*ICRTCR) .LT. THRSTP) CNVRGD(3) = .TRUE.
         IF (STPMAX .LT. THSTMX) CNVRGD(4)                   = .TRUE.
C
         IF (CONOPT) THEN
C
C     For constrained optimizations we examine the change in the RMS gradient
C     and the change in the maximum element. The threshold for these differences
C     are 1/10 of the current regular gradient threshold.
C            
            IF (ITRNMR .EQ. 0) THEN
               GRDTHR = MAX(0.1D0*GRDTHR,1.D-7)
               THGRMX = 1.5D0*GRDTHR
               DIFRMS = 0.D0
               DIFMAX = 0.D0
               PRVRMS = GRADNM/SQRT(1.D0*ICRTCR)
               PRVMAX = GRDMAX
            ELSE
               DIFRMS = ABS(PRVRMS-(GRADNM/SQRT(1.D0*ICRTCR)))
               DIFMAX = ABS(PRVMAX-GRDMAX)
               PRVRMS = GRADNM/SQRT(1.D0*ICRTCR)
               PRVMAX = GRDMAX
            END IF
            IF (DIFRMS .LT. GRDTHR) CNVRGD(1) = .TRUE.
            IF (DIFMAX .LT. THGRMX) CNVRGD(2) = .TRUE.
            WRITE(LUPRI,'(///A)') ' -----------------------' //
     &     '-----------------------------------------------------------'
            WRITE(LUPRI,'(14X,A)') 'RMS grad    Max grad    ' //
     &           'RMS diff*   Max diff*   RMS step    Max step'
            WRITE(LUPRI,'(A)')' -----------------------' //
     &     '-----------------------------------------------------------'
            WRITE(LUPRI,'(A,6F12.8)') ' Curr.val. ',
     &           GRADNM/SQRT(1.D0*ICRTCR), GRDMAX,
     &           DIFRMS, DIFMAX,
     &           STPNRM/SQRT(1.D0*ICRTCR), STPMAX
            WRITE(LUPRI,'(A,24X,4F12.8)') ' Threshold ',
     &           GRDTHR, THGRMX, THRSTP, THSTMX
            WRITE(LUPRI,'(A,14X,4(10X,L1))') ' Converged?      ',
     &           (CNVRGD(I), I = 1,4)
            WRITE(LUPRI,'(A)')' -----------------------' //
     &     '-----------------------------------------------------------'
            WRITE(LUPRI,'(A)') '*) For constrained optimizations' //
     &           ' the change in RMS gradient and the change'
            WRITE(LUPRI,'(A///)') '   in the maximum gradient ' //
     &           'element are used as convergence criteria.'
         ELSE
            WRITE(LUPRI,'(///A)')
     &     ' ----------------------------------------------------------'
            WRITE(LUPRI,'(14X,A)')
     &           'RMS grad    Max grad    RMS step    Max step'
            WRITE(LUPRI,'(A)')
     &     ' ----------------------------------------------------------'
            WRITE(LUPRI,'(A,4F12.8)') ' Curr.val. ',
     &           GRADNM/SQRT(1.D0*ICRTCR), GRDMAX,
     &           STPNRM/SQRT(1.D0*ICRTCR), STPMAX
            WRITE(LUPRI,'(A,4F12.8)') ' Threshold ',
     &           GRDTHR, THGRMX, THRSTP, THSTMX
            WRITE(LUPRI,'(A,4L12)')   ' Converged?',
     &           (CNVRGD(I), I = 1,4)
            WRITE(LUPRI,'(A//)')
     &     ' ----------------------------------------------------------'
         END IF
C
         INDXOK = (INDTOT .EQ. 0)
         IF (SADDLE) INDXOK = (INDTOT .EQ. 1)
         MINEND = INDXOK .AND. CNVRGD(1) .AND. CNVRGD(2) .AND.
     &        CNVRGD(3) .AND. CNVRGD(4)
         IF (NOBRKS) MINEND = CNVRGD(1) .AND. CNVRGD(2) .AND.
     &        CNVRGD(3) .AND. CNVRGD(4)
      END IF
C
C     If the energy difference is below THRSYM and the total Hessian is
C     non-zero, we can be pretty sure that we're approaching a saddle
C     point, and the symmetry should be broken. Only applies to
C     second order methods.
C
      IF ((NEWTON .OR. QUADSD) .AND. (INDTOT .GT. 0) .AND.
     &     (EDIFF .LE. MAX(THRERG,THRSYM)) .AND.
     &     ((ITRNMR-1) .NE. ITRBRK) .AND. (.NOT. SADDLE)) THEN
         BRKSYM = .TRUE.
         RSTARR = .TRUE.
         NWSYMM = .TRUE.
      END IF
C
C     Output from the testing is written.
C
      IF (.NOT. NATNRM .AND. (IPRINT .GT. 2)) THEN
         CALL TITLER('Output from MINEND','*',103)
         LOGTXT = 'no  '
         IF (CNVERG) LOGTXT = 'yes '
         IF ((ITRNMR .LT. 1) .OR. (ITRBRK .EQ. (ITRNMR-1)))
     &        LOGTXT = 'N/A '
         WRITE(LUPRI,'(/A,1P,A5)') ' Energy converged      ',LOGTXT
         LOGTXT = 'no  '
         IF (CNVGRD) LOGTXT = 'yes '
         WRITE(LUPRI,'(A,1P,A5)') ' Gradient converged    ',LOGTXT
         LOGTXT = 'no  '
         IF (CNVSTP) LOGTXT = 'yes '
         WRITE(LUPRI,'(A,1P,A5)') ' Step converged        ',LOGTXT
         IF (.NOT. BAKER) THEN
            WRITE(LUPRI,'(A,1P,I3)') ' Conditions fullfilled ',ICONV
            WRITE(LUPRI,'(A,1P,I3)') ' Required conditions   ',ICONDI
         END IF
         WRITE(LUPRI,'(A,1P,I3)') ' Totally sym. index    ',INDHES(0)
         WRITE(LUPRI,'(A,1P,I3)') ' Hessian index         ',INDTOT
         LOGTXT = 'no  '
         IF (MINEND) LOGTXT = 'yes '
         IF (BRKSYM) LOGTXT = 'yes*'
         WRITE(LUPRI,'(A,1P,A5)') ' End of optimization   ',LOGTXT
         IF (BRKSYM) WRITE(LUPRI,'(/A)') ' *) Within given symmetry.'
      END IF
      RETURN
      END

C  /* Deck spnrgy */
      SUBROUTINE SPNRGY(GEINFO,EXHER,EXSIR,EXABA,WORK,LWORK,WRKDLM)
C
C     Calculates single point energy with another basis if requested.
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "molinp.h"
#include "optinf.h"
#include "gnrinf.h"
      DIMENSION GEINFO(0:ITRMAX,6)
      LOGICAL EXHER,EXSIR,EXABA,REJGEO
      DIMENSION WORK(LWORK)
      EXHER  = .FALSE.
      EXSIR  = .FALSE.
      EXABA  = .FALSE.
      RDINPC = .FALSE.
      RDMLIN = .FALSE.
      HRINPC = .FALSE.
      MLINE(2) = SPBSTX
C
C     New input files are written with desired basis. Iteration number
C     has to be decreased first in order to write to the correct filename.
C
      CALL PNCMOL(ITRNMR,IPRINT)
C
C     Initialization of new memory trap for programs.
C
      WORK(1) = WRKDLM
C
      NEWSYM = .TRUE.
      CALL GTNRGY(EXHER,EXSIR,EXABA,.FALSE.,WORK,LWORK,WRKDLM)
      GEINFO(ITRNMR+1,1) = ENERGY
      RETURN
      END

C  /* Deck lshft0 */
      SUBROUTINE LSHFT0(NCORD,NONTRO,EVAL,GRDDIA,STPDIA,
     &     TRUSTR,RNU,KEEPSY,ZERGRD,INSIDE,IPRINT)
C
C     (Almost identical to WLKFL0 in abawalk.F)
C     This subroutine solves the constrained restricted step
C     equations (the level-shifted Newton equations) in the
C     diagonal representation.  We assume that the Newton step
C     is longer than the trust radius.
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER ( D0 = 0.0D0 , DP5 = 0.5D0 )
C
      DIMENSION EVAL(*), GRDDIA(*), STPDIA(*)
      LOGICAL KEEPSY, INSIDE, SPECAS
C
      EXTERNAL WSTPLN
C
      IF (IPRINT .GT. 5) CALL HEADER('OUTPUT FROM LSHFT0',-1)
C
      IF (KEEPSY) THEN
         DO 10 I = 1, NONTRO
            IF (ABS(GRDDIA(I)) .LT. ZERGRD) GRDDIA(I) = D0
   10    CONTINUE
         STPSYM = WSTPLN(GRDDIA(1),EVAL(1),D0,NONTRO,D0)
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
      GRDNRM = SQRT(DDOT(NONTRO,GRDDIA,1,GRDDIA,1))
      IF (IPRINT .GT. 5) THEN
         WRITE (LUPRI,'(A,1P,D12.5)') ' HESMIN: ', HESMIN
         WRITE (LUPRI,'(A,1P,D12.5)') ' GRDMIN: ', GRDMIN
         WRITE (LUPRI,'(A,1P,D12.5)') ' GRDNRM: ', GRDNRM
         WRITE (LUPRI,'(A,1P,D12.5)') ' ZERGRD: ', ZERGRD
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
     *         ' XMIN and XMAX before WLKBIS (WSTPLN) minimum walk: ',
     *          XMIN,XMAX
         END IF
         CALL WLKBIS(XMAX,XMIN,RNU,GRDDIA,EVAL,TRUSTR,NONTRO,
     *               WSTPLN,IFAIL)
         IF (IPRINT .GT. 5) THEN
            WRITE (LUPRI,'(/A,2(1P,D12.5,2X))')
     *         ' XMIN and XMAX after WLKBIS (WSTPLN) minimum walk: ',
     *          XMIN,XMAX
         END IF
 33      CONTINUE
         IF (IPRINT .GT. 3) WRITE (LUPRI,'(/A,1P,D12.5)')
     *      ' Level shift parameter: ', RNU
C
C     If we're using level-shift for a step less than the trust radius,
C     and we're having trouble with the interval, we simply set
C     the level-shift parameter to zero.
C
         IF ((IFAIL.EQ.0) .AND. INSIDE) THEN
            WRITE (LUPRI,'(/A)')
     *           ' *** ERROR, Wrong interval in WLKBIS (WSTPLN)'
            WRITE (LUPRI,'(A)')
     *           '     Setting level-shift equal to zero.'
            RNU = D0
            IFAIL = -1
            GOTO 33
         ELSE IF (IFAIL.EQ.0) THEN
            WRITE (LUPRI,5250) XMAX, XMIN
            CALL QUIT(' *** ERROR, Wrong interval in WLKBIS (WSTPLN)')
         ELSE IF (IFAIL.EQ.1) THEN
            WRITE (LUPRI,5350)
         END IF
C
C        Determine step vector
C
         DO 100 I = 1, NONTRO
            IF (ABS(EVAL(I) + RNU) .LE. 1.0D-8) THEN
               STPDIA(I) = D0
            ELSE
               STPDIA(I) = - GRDDIA(I)/(EVAL(I) + RNU)
            END IF
 100     CONTINUE
C
C     *************************************************
C     ***** Special case: HESMIN < 0 & GRDMIN = 0 *****
C     *************************************************
C
      ELSE
         IF (IPRINT .GT. 3) WRITE (LUPRI,'(/A/)') ' Special case.'
         IMODE = 2
 150     CONTINUE
         IF (IMODE .LT. NONTRO) THEN
            IF ((EVAL(IMODE) .LT. D0) .AND.
     &           (ABS(GRDDIA(IMODE)) .LT. ZERGRD)) THEN
               IMODE = IMODE + 1
               GOTO 150
            END IF
         END IF
         GRDMIN=GRDDIA(IMODE)
         HESMIN=EVAL(IMODE)
C
C        Set RNU = - HESMIN and determine step length
C
         STPNRM = WSTPLN(GRDDIA(IMODE),EVAL(IMODE),-HESMIN,
     &        NONTRO-IMODE+1,D0)
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
     *           ' XMIN and XMAX before WLKBIS (WSTPLN) minimum walk: ',
     *           XMIN,XMAX
            END IF
            CALL WLKBIS(XMAX,XMIN,RNU,GRDDIA(IMODE),EVAL(IMODE),TRUSTR,
     &           NONTRO-IMODE+1,WSTPLN,IFAIL)
            IF (IPRINT .GT. 5) THEN
               WRITE (LUPRI,'(/A,2(1P,D12.5,2X))')
     *         ' XMIN and XMAX after WLKBIS (WSTPLN) minimum walk: ',
     *           XMIN,XMAX
            END IF
            IF (IPRINT .GT. 3) WRITE (LUPRI,'(/A,1P,D12.5)')
     *         ' Level shift parameter: ', RNU
            IF (IFAIL.EQ.0) THEN
               WRITE (LUPRI,5250) XMAX, XMIN
               CALL QUIT
     *            (' *** ERROR, Wrong interval in WLKBIS (WSTPLN)')
            ELSE IF (IFAIL.EQ.1) THEN
               WRITE (LUPRI,5350)
            END IF
C
C           Determine step vector
C
            DO 199 I = 1, IMODE-1
               STPDIA(I) = D0
 199        CONTINUE
            DO 200 I = IMODE, NONTRO
               STPDIA(I) = - GRDDIA(I)/(EVAL(I) + RNU)
 200        CONTINUE
         ELSE
C
C           Determine step vector with level shift - HESMIN and add
C           component along the lowest eigenvector(s) to insure that total
C           step length is equal to the trust radius.
C
            DO 300 I = IMODE, NONTRO
               STPDIA(I) = - GRDDIA(I)/(EVAL(I) - HESMIN)
 300        CONTINUE
            STP2 = DDOT(NONTRO-IMODE+1,STPDIA(IMODE),1,STPDIA(IMODE),1)
            IF (IPRINT .GT. 3) WRITE (LUPRI,'(/A,F12.6)')
     *         ' Norm of step orthogonal to lowest eigenvector(s):',
     *         SQRT(STP2)
            SXTRA =  SQRT(TRUSTR*TRUSTR - STP2)/SQRT(1.0D0*(IMODE-1))
            IF (IPRINT .GT. 3) WRITE (LUPRI,'(A,F12.6)')
     *         ' Norm of step parallel to lowest eigenvector(s):  ',
     *         SQRT((IMODE-1)*SXTRA*SXTRA)
            DO 400 I = 1, IMODE-1
               STPDIA(I) = SXTRA
 400        CONTINUE
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
 5250 FORMAT(/' *** Wrong interval specified in WLKBIS (WSTPLN) ***',
     *       /' XMAX= ',F10.6,'   XMIN= ',F10.6)
 5350 FORMAT(/' *** WARNING WLKBIS (WSTPLN) ***',
     *       /' Desired accuracy not obtained in the specified maximum',
     *       /' number of iterations.')
      END
 
C  /* Deck pncmol */
      SUBROUTINE PNCMOL(ITER_GEO, IPRINT)
C
C     Punch MOLECULE input with updated coordinates to MOLECULE.INP
C     and XXXX_mol.inp + optional VRML-file
C
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "mxcent.h"
#include "inftap.h"
#include "molinp.h"
C
      CHARACTER*12 FILENM
      PARAMETER(IPRINT_LIM = 2) !hjaaj DEBUG; set IPRINT_LIM to e.g. -1
C
      IF (IPRINT .GT. IPRINT_LIM .OR. NMLINE .LE. 0) THEN
         CALL HEADER('New MOLECULE.INP punched',-1) 
         WRITE(LUPRI,'(A,I5/)') 'Number of lines, NMLINE =',NMLINE
         IF (NMLINE .LE. 0)
     &      CALL QUIT('Fatal error in PNCMOL, NMLINE .le. 0')
      END IF
      CALL GPOPEN(LUMOL,'MOLECULE.INP','UNKNOWN',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND (LUMOL)
      DO IMLINE = 1,NMLINE
         WRITE(LUMOL,'(A)') MLINE(IMLINE)
         IF (IPRINT .GT. IPRINT_LIM)
     &      WRITE (LUPRI,'(I4,2A)') IMLINE,' >',MLINE(IMLINE)
      END DO
      CALL GPCLOSE(LUMOL,'KEEP')
C
      IF (ITER_GEO .GE. 0) THEN
         FILENM = 'XXXX_mol.inp'
         WRITE(FILENM(1:4),'(I4)') (ITER_GEO + 1)
         IF (ITER_GEO .LT. 9)   WRITE(FILENM(3:3),'(A1)') '0'
         IF (ITER_GEO .LT. 99)  WRITE(FILENM(2:2),'(A1)') '0'
         IF (ITER_GEO .LT. 999) WRITE(FILENM(1:1),'(A1)') '0'
         CALL GPOPEN(LUMOL,FILENM,'UNKNOWN',' ','FORMATTED',
     &      IDUMMY,.FALSE.)
         DO IMLINE = 1,NMLINE
            WRITE(LUMOL,'(A)') MLINE(IMLINE)
         END DO
         CALL GPCLOSE(LUMOL,'KEEP')
      END IF
      RETURN
      END

C  /* Deck pnches */
      SUBROUTINE PNCHES(MXRCRD,MX2CRD,HESINT,WILBMT,BMTRAN,
     &     TMPMT1,TMPMT2,TMPMT3,TMPMT4,WORK,LWORK)
C
C     Punch molecular Hessian to the file DALTON.HES, this file can be
C     used to obtain initial Hessian (1st order methods) for restarts
C     or other runs.
C
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "molinp.h"
#include "optinf.h"
#include "gnrinf.h"
#include "taymol.h"
#include "symmet.h"
      DIMENSION HESINT(MXRCRD,MXRCRD)
      DIMENSION WILBMT(MXRCRD,MXCOOR), BMTRAN(MXRCRD,MXRCRD)
      DIMENSION TMPMT1(MX2CRD,MX2CRD), TMPMT2(MX2CRD,MX2CRD)
      DIMENSION TMPMT3(MX2CRD,MX2CRD), TMPMT4(MXCOOR,MXCOOR)
      DIMENSION WORK(LWORK)
C
      LUHES = -1
      CALL GPOPEN(LUHES,'DALTON.HES','UNKNOWN',' ','FORMATTED',IDUMY,
     &     .FALSE.)
      IF (REDINT .OR. DELINT) THEN
         CALL DZERO(TMPMT4,MXCOOR*MXCOOR)
         CALL HQ2HX(MXRCRD,MX2CRD,TMPMT1,TMPMT2,TMPMT3,HESINT,GRDINT,
     &        TMPMT4,WILBMT,BMTRAN,WORK,LWORK)
         WRITE(LUHES,'(I6)') 3*NUCDEP
         WRITE(LUHES,*)
         DO 100 J = 1, 3*NUCDEP
            DO 110 I = 1, 3*NUCDEP
               WRITE(LUHES,'(F20.16)') TMPMT4(I,J)
 110        CONTINUE
            WRITE(LUHES,*)
 100     CONTINUE
      ELSE
         IF (MAXREP .GT. 0) THEN
            CALL DZERO(TMPMT4,MXCOOR*MXCOOR)
            CALL TRAHES(HESMOL,MXCOOR,TMPMT4,TMPMT1,TMPMT2,
     &           MXCOOR,3*NUCDEP,1)
            WRITE(LUHES,*) 3*NUCDEP
            WRITE(LUHES,*)
            DO 200 J = 1, 3*NUCDEP
               DO 210 I = 1, 3*NUCDEP
                  WRITE(LUHES,'(F20.16)') TMPMT4(I,J)
 210           CONTINUE
               WRITE(LUHES,*)
 200        CONTINUE
         ELSE
            WRITE(LUHES,*) NCART
            WRITE(LUHES,*)
            DO 250 J = 1, NCART
               DO 260 I = 1, NCART
                  WRITE(LUHES,'(F20.16)') HESMOL(I,J)
 260           CONTINUE
               WRITE(LUHES,*)
 250        CONTINUE
         END IF
      END IF
      CALL GPCLOSE(LUHES,'KEEP')
      RETURN
      END

C  /* Deck reahes */
      SUBROUTINE REAHES(MXRCRD,MX2CRD,HESINT,ATMARR,TMPMT1,TMPMT2,
     &     TMPMT3,TMPMT4,WILBMT,BMTRAN,BMTINV,WORK,LWORK,IERR)
C
C     Read molecular Hessian from the file DALTON.HES, which is then
C     used as initial Hessian in 1st order methods/restarts. IERR is
C     returned with the value 0 if everything is OK. -1 indicates that
C     the file cannot be opened, -2 that the Hessian in the file has
C     wrong dimensions.
C
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "molinp.h"
#include "optinf.h"
#include "gnrinf.h"
#include "taymol.h"
#include "symmet.h"
      LOGICAL HESEXS
      DIMENSION HESINT(MXRCRD,MXRCRD), ATMARR(MXCENT,8)
      DIMENSION TMPMT1(MX2CRD,MX2CRD), TMPMT2(MX2CRD,MX2CRD)
      DIMENSION TMPMT3(MX2CRD,MX2CRD), TMPMT4(MXCOOR,MXCOOR)
      DIMENSION WILBMT(MXRCRD,MXCOOR), BMTRAN(MXRCRD,MXRCRD)
      DIMENSION BMTINV(MXRCRD,MXCOOR), WORK(LWORK)
C
      LUHES = -1
      CALL DZERO(HESMOL,MXCOOR*MXCOOR)
      INQUIRE(FILE='DALTON.HES',EXIST=HESEXS)
      IF (.NOT. HESEXS) THEN
         IERR = -1
         RETURN
      ELSE
         IERR = 0
      END IF
      CALL GPOPEN(LUHES,'DALTON.HES','OLD',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      READ(LUHES,*) IDIM
      READ(LUHES,*)
      ICRD = 3*NUCDEP
      IF (IDIM .NE. ICRD) THEN
         IERR = -2
         RETURN
      END IF
      CALL DZERO(TMPMT4,MXCOOR*MXCOOR)
      DO 100 J = 1, ICRD
         DO 110 I = 1, ICRD
            READ(LUHES,*) TMPMT4(I,J)
 110     CONTINUE
         READ(LUHES,*)
 100  CONTINUE
      IF (MAXREP .GT. 0) THEN
         CALL TRACOR(TMPMT2,TMPMT3,1,ICRD,0)
         CALL DGEMM('T','N',ICRD,ICRD,ICRD,1.D0,
     &        TMPMT3,ICRD,
     &        TMPMT4,MXCOOR,0.D0,
     &        TMPMT1,MXCOOR)
         CALL DZERO(HESMOL,MXCOOR*MXCOOR)
         CALL DGEMM('N','N',ICRD,ICRD,ICRD,1.D0,
     &        TMPMT1,MXCOOR,
     &        TMPMT3,ICRD,0.D0,
     &        HESMOL,MXCOOR)
      ELSE
         DO 150 J = 1, ICRD
            DO 160 I = 1, ICRD
               HESMOL(I,J) = TMPMT4(I,J)
 160        CONTINUE
 150     CONTINUE
      END IF
      IF (REDINT .OR. DELINT) THEN

         CALL DZERO(HESINT,MXRCRD*MXRCRD)
         CALL HX2HQ(MXRCRD,MX2CRD,ATMARR,TMPMT1,TMPMT2,TMPMT3,
     &        TMPMT4,GRDINT,HESINT,WILBMT,BMTINV,BMTRAN,
     &        WORK,LWORK)

      END IF
      CALL GPCLOSE(LUHES,'KEEP')
      RETURN
      END

C  /* Deck maxelm */
      SUBROUTINE MAXELM(VEC,IDIM,SCLVEC,ISCL,ELMMX)
C
C     Finds the largest elemement (absolute value) of the vector VEC
C     of dimension IDIM. The value is returned through the
C     variable ELMMX.
C
#include "implicit.h"
#include "mxcent.h"
#include "priunit.h"
      DIMENSION VEC(IDIM),SCLVEC(MXCOOR)
      LOGICAL SCALE
      ELMMX = ABS(VEC(1))
      IF (ISCL .GE. 1) ELMMX = ELMMX/SCLVEC(1)
      IF (ISCL .GE. 2) ELMMX = ELMMX/SCLVEC(1)
      DO 10 I = 2, IDIM
         ELM = ABS(VEC(I))
         IF (ISCL .GE. 1) ELM = ELM/SCLVEC(I)
         IF (ISCL .GE. 2) ELM = ELM/SCLVEC(I)
         IF (ELM .GT. ELMMX) ELMMX = ELM
 10   CONTINUE
      RETURN
      END

C  /* Deck wstpln */
      FUNCTION WSTPLN(GDDIA,HESDIA,RNU,NCORD,RTRUST)
C
C     (almost identical to WLKSTL in abawalk.F)
C     Purpose:
C
C        Calculate step length at level shift RNU and
C        subtract RTRUST
C
C        WSTPLN = //STEP// - RTRUST
C
C        where
C
C        STEP = - GDDIA/(HESDIA+RNU)
C
#include "implicit.h"
      DIMENSION GDDIA(*),HESDIA(*)
      PARAMETER(D0=0.0D0, ZERO=1.0D-8 )
      STEP = D0
      DO 100 K=1,NCORD
         IF ((ABS(GDDIA(K)) .GE. ZERO) .AND.
     &        (ABS(HESDIA(K)) .GE. ZERO)) THEN
            STEPK = GDDIA(K) / (HESDIA(K)+RNU)
            STEP = STEP + STEPK*STEPK
         END IF
 100  CONTINUE
      WSTPLN = SQRT(STEP) - RTRUST
      RETURN
      END

C  /* Deck numgrd */
      SUBROUTINE NUMGRD(WORK,LWORK,WRKDLM)
C
C     Performs a numerical differentiation in order to get an
C     finite-difference molecular gradient for use in geometry optimization
C     of MP2, CI or CC wave functions for which there does not exist
C     analytical implementations, K.Ruud, Feb.-6 1997
C
C     Modified in august of 2000 to use the new numerical derivative routines
C
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "inftap.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
C
      LOGICAL FIRST
      PARAMETER (DP5 = 0.50D0, D100=100.0D0)
      CHARACTER*7 WORD
      DIMENSION WORK(LWORK), COOR(3,MXCENT)
C
#include "molinp.h"
#include "taymol.h"
#include "numder.h"
#include "symmet.h"
#include "infpar.h"
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
#endif
#include "nuclei.h"
#include "pgroup.h"
#include "fcsym.h"
#include "gnrinf.h"
#include "optinf.h"
#include "infopt.h"

#if defined (VAR_PARIO)
      SAVE FIRST
      DATA FIRST /.TRUE./
C
! WARNING! do not activate VAR_PARIO without fixing
!          the PARIO code !!!! /Feb 2011 hjaaj
      IF (NODTOT .GT. 1 .AND. MYNUM .EQ. 0) THEN
         IF (FIRST) THEN
            FIRST = .FALSE.
         ELSE
            CALL PARIOT
         END IF
      END IF
#endif

      CALL QENTER('NUMGRD')
      IF (.NOT. CHGRDT) THEN
         GRDTHR = 1.0D-4
         THRSTP = GRDTHR
         THRERG = GRDTHR
         THRSYM = SQRT(THRERG)
         WRITE (LUPRI,'(/A)') ' INFO >>>> Due to limitations '//
     &        'of the accuracy of the numerical gradients'
         WRITE (LUPRI,'(A)') ' INFO >>>> thresholds for convergence'
     &        //' of geometry optimization has been reset'
         WRITE (LUPRI,'(/A,/,3(/,20X,A,F10.8))') ' New thresholds:',
     &        'Gradient norm  ',GRDTHR,'Step norm      ',THRSTP,
     &        'Energy change  ',THRERG
      END IF
C
C     Call new numerical derivative routines.
C
      CALL NMDINI(IPRINT)
      IF (LUCMD .LE. 0) CALL GPOPEN(LUCMD,'DALTON.INP','OLD',
     &   ' ','FORMATTED',IDUMMY,.FALSE.)
      REWIND (LUCMD,IOSTAT=IOS)
 1100 READ (LUCMD,'(A7)',END=1110) WORD
        CALL UPCASE(WORD)
      IF (WORD .NE. '**NMDDR') GOTO 1100
C
      CALL NMDINP(WORD,IPRINT)
 1110 CALL GPCLOSE(LUCMD,'KEEP')
      FCLASS = GROUP
      IF (NMORDR .LT. 1) NMORDR = 1
      CALL NUMDRV(WORK,LWORK,IPRINT,WRKDLM)
C
      CALL QEXIT('NUMGRD')
      RETURN
      END
C --- end of abaopt.F ---
