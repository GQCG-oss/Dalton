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
C  /* Deck hypinp */
      SUBROUTINE HYPINP(WORD)
C
      use pelib_interface, only: use_pelib
#include "implicit.h"
C
#include "priunit.h"
#include "infrsp.h"
#include "rspprp.h"
#include "inflr.h"
#include "infhyp.h"
#include "infpri.h"
#include "infspi.h"
#include "maxorb.h"
#include "infinp.h"
#include "inftap.h"
#include "infrank.h"
#include "inforb.h"
! gnrinf.h: QM3, THR_REDFAC
#include "gnrinf.h"
C
      LOGICAL NEWDEF, BCFREQ, COMFREQ, NXTLAB
      PARAMETER ( NTABLE = 30, ZERO = 0.0D0 , THRFRQ = 1.0D-14 )
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER*8 LABEL, LABFND, RTNLBL(2)
C
      DATA TABLE /'.BPROP ', '.BFREQ ', '.CPROP ', '.CFREQ ',
     *            '.APROP ', '.MAX IT', '.THCLR ', '.MAXITO',
     *            '.PRINT ', '.DIPLEN', '.E3TEST', '.XXXXXX',
     *            '.ISPABC', '.A2TEST', '.X2TEST', '.REFCHK',
     *            '.TSTJEP', '.SHG   ', '.POCKEL', '.OPTREF',
     *            '.FREQUE', '.DIPLNX', '.DIPLNY', '.DIPLNZ',
     *            '.ASPIN ', '.BSPIN ', '.CSPIN ', '.SOSHIE',
     *            '.SOSPIN', '.QLOP '/
C
C READ quadratic response INPUT
C
      NEWDEF = (WORD .EQ. '*QUADRA')
      ICHANG = 0
      IF (NEWDEF) THEN
      HYPCAL = .TRUE.
      BCFREQ = .FALSE.
      COMFREQ= .FALSE.
      QLOP = .FALSE.
      WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') GO TO 100
            IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,
     *                    14,15,16,17,18,19,20,21,22,23,24,
     *                    25,26,27,28,29,30),I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' KEYWORD "',WORD,
     *            '" NOT RECOGNIZED IN RSPQR.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT(' ILLEGAL KEYWORD IN *QUADRA')
    1          CONTINUE
                  READ(LUCMD,'(BN,A,I8)')LABEL, IRANKB
                  IF (LABEL .EQ. 'ANGMOM  ') THEN
                     BQROP(INDPRP('XANGMOM ')) = .TRUE.
                     BQROP(INDPRP('YANGMOM ')) = .TRUE.
                     BQROP(INDPRP('ZANGMOM ')) = .TRUE.
                  ELSE IF (LABEL .EQ. '1SPNORB ') THEN
                     BQROP(INDPRP('X1SPNORB')) = .TRUE.
                     BQROP(INDPRP('Y1SPNORB')) = .TRUE.
                     BQROP(INDPRP('Z1SPNORB')) = .TRUE.
                  ELSE IF (LABEL .EQ. '2SPNORB ') THEN
                     BQROP(INDPRP('X2SPNORB')) = .TRUE.
                     BQROP(INDPRP('Y2SPNORB')) = .TRUE.
                     BQROP(INDPRP('Z2SPNORB')) = .TRUE.
                  ELSE IF (LABEL .EQ. 'MNFSPNOR') THEN
                     BQROP(INDPRP('X1MNF-SO')) = .TRUE.
                     BQROP(INDPRP('Y1MNF-SO')) = .TRUE.
                     BQROP(INDPRP('Z1MNF-SO')) = .TRUE.
                  ELSE IF (LABEL .EQ. 'FERMI CO') THEN
                     CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',
     &                           'UNFORMATTED',IDUMMY,.FALSE.)
 199                 CONTINUE
                     IF (NXTLAB(LABFND,RTNLBL,LUPROP)) THEN
                        IF (LABFND(1:3) .EQ. 'FC ')
     &                       BQROP(INDPRP(LABFND)) = .TRUE.
                        GOTO 199
                     END IF
                     CALL GPCLOSE(LUPROP,'KEEP')
                  ELSE IF (LABEL .EQ. 'SPIN-DIP') THEN
                     CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',
     &                           'UNFORMATTED',IDUMMY,.FALSE.)
 299                 CONTINUE
                     IF (NXTLAB(LABFND,RTNLBL,LUPROP)) THEN
                        IF (LABFND(1:3) .EQ. 'SD ')
     &                       BQROP(INDPRP(LABFND)) = .TRUE.
                        GOTO 299
                     END IF
                     CALL GPCLOSE(LUPROP,'KEEP')
                  ELSE IF (LABEL(1:3) .EQ. 'PSO' .OR.
     &                     LABEL(2:4) .EQ. 'PSO') THEN
                     CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',
     &                           'UNFORMATTED',IDUMMY,.FALSE.)
 399                 CONTINUE
                     IF (NXTLAB(LABFND,RTNLBL,LUPROP)) THEN
                        IF (LABFND(1:3) .EQ. 'PSO')
     &                       BQROP(INDPRP(LABFND)) = .TRUE.
                        GOTO 399
                     END IF
                     CALL GPCLOSE(LUPROP,'KEEP')
                  ELSE
                     BQROP( INDPRP(LABEL)) = .TRUE.
                     OPRANK( INDPRP(LABEL)) = IRANKB
                  END IF
               GO TO 100
    2          CONTINUE
                  BCFREQ = .TRUE.
                  READ (LUCMD,*) NBQRFR
                  IF (NBQRFR.LE.MBQRFR) THEN
                     READ (LUCMD,*) (BQRFR(J),J=1,NBQRFR)
                  ELSE
                     WRITE (LUPRI,'(3(/,A,I5),/)')
     *               ' NUMBER OF FREQUENCIES SPECIFIED    :',NBQRFR,
     *               ' IS GREATER THAN THE ALLOWED NUMBER :',MBQRFR,
     *               ' THE NUMBER IS RESET TO THE MAXIMUM :',MBQRFR
                     READ (LUCMD,*) (BQRFR(J),J=1,MBQRFR),
     *                              (FFFF,J=MBQRFR+1,NBQRFR)
                     NBQRFR = MBQRFR
                  END IF
               GO TO 100
    3          CONTINUE
                  READ(LUCMD,'(BN,A,I8)')LABEL,IRANKC
                  IF (LABEL .EQ. 'ANGMOM  ') THEN
                     CQROP(INDPRP('XANGMOM ')) = .TRUE.
                     CQROP(INDPRP('YANGMOM ')) = .TRUE.
                     CQROP(INDPRP('ZANGMOM ')) = .TRUE.
                  ELSE IF (LABEL .EQ. '1SPNORB ') THEN
                     CQROP(INDPRP('X1SPNORB')) = .TRUE.
                     CQROP(INDPRP('Y1SPNORB')) = .TRUE.
                     CQROP(INDPRP('Z1SPNORB')) = .TRUE.
                  ELSE IF (LABEL .EQ. '2SPNORB ') THEN
                     CQROP(INDPRP('X2SPNORB')) = .TRUE.
                     CQROP(INDPRP('Y2SPNORB')) = .TRUE.
                     CQROP(INDPRP('Z2SPNORB')) = .TRUE.
                  ELSE IF (LABEL .EQ. 'MNFSPNOR') THEN
                     CQROP(INDPRP('X1MNF-SO')) = .TRUE.
                     CQROP(INDPRP('Y1MNF-SO')) = .TRUE.
                     CQROP(INDPRP('Z1MNF-SO')) = .TRUE.
                  ELSE IF (LABEL .EQ. 'FERMI CO') THEN
                     CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',
     &                           'UNFORMATTED',IDUMMY,.FALSE.)
 198                 CONTINUE
                     IF (NXTLAB(LABFND,RTNLBL,LUPROP)) THEN
                        IF (LABFND(1:3) .EQ. 'FC ')
     &                       CQROP(INDPRP(LABFND)) = .TRUE.
                        GOTO 198
                     END IF
                     CALL GPCLOSE(LUPROP,'KEEP')
                  ELSE IF (LABEL .EQ. 'SPIN-DIP') THEN
                     CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',
     &                           'UNFORMATTED',IDUMMY,.FALSE.)
 298                 CONTINUE
                     IF (NXTLAB(LABFND,RTNLBL,LUPROP)) THEN
                        IF (LABFND(1:3) .EQ. 'SD ')
     &                       CQROP(INDPRP(LABFND)) = .TRUE.
                        GOTO 298
                     END IF
                     CALL GPCLOSE(LUPROP,'KEEP')
                  ELSE IF (LABEL(1:3) .EQ. 'PSO' .OR.
     &                     LABEL(2:4) .EQ. 'PSO') THEN
                     CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',
     &                           'UNFORMATTED',IDUMMY,.FALSE.)
 398                 CONTINUE
                     IF (NXTLAB(LABFND,RTNLBL,LUPROP)) THEN
                        IF (LABFND(1:3) .EQ. 'PSO')
     &                       CQROP(INDPRP(LABFND)) = .TRUE.
                        GOTO 398
                     END IF
                     CALL GPCLOSE(LUPROP,'KEEP')
                  ELSE
                     CQROP( INDPRP(LABEL)) = .TRUE.
                     OPRANK( INDPRP(LABEL)) = IRANKC
                  END IF
               GO TO 100
    4          CONTINUE
                  BCFREQ = .TRUE.
                  READ (LUCMD,*) NCQRFR
                  IF (NCQRFR.LE.MCQRFR) THEN
                     READ (LUCMD,*) (CQRFR(J),J=1,NCQRFR)
                  ELSE
                     WRITE (LUPRI,'(3(/,A,I5),/)')
     *               ' NUMBER OF FREQUENCIES SPECIFIED    :',NCQRFR,
     *               ' IS GREATER THAN THE ALLOWED NUMBER :',MCQRFR,
     *               ' THE NUMBER IS RESET TO THE MAXIMUM :',MCQRFR
                     READ (LUCMD,*) (CQRFR(J),J=1,MCQRFR),
     *                              (FFFF,J=MCQRFR+1,NCQRFR)
                     NCQRFR = MCQRFR
                  END IF
               GO TO 100
    5          CONTINUE
                  READ(LUCMD,'(BN,A,I8)')LABEL,IRANKA
                  IF (LABEL .EQ. 'ANGMOM  ') THEN
                     AQROP(INDPRP('XANGMOM ')) = .TRUE.
                     AQROP(INDPRP('YANGMOM ')) = .TRUE.
                     AQROP(INDPRP('ZANGMOM ')) = .TRUE.
                  ELSE IF (LABEL .EQ. '1SPNORB ') THEN
                     AQROP(INDPRP('X1SPNORB')) = .TRUE.
                     AQROP(INDPRP('Y1SPNORB')) = .TRUE.
                     AQROP(INDPRP('Z1SPNORB')) = .TRUE.
                  ELSE IF (LABEL .EQ. '2SPNORB ') THEN
                     AQROP(INDPRP('X2SPNORB')) = .TRUE.
                     AQROP(INDPRP('Y2SPNORB')) = .TRUE.
                     AQROP(INDPRP('Z2SPNORB')) = .TRUE.
                  ELSE IF (LABEL .EQ. 'MNFSPNOR') THEN
                     AQROP(INDPRP('X1MNF-SO')) = .TRUE.
                     AQROP(INDPRP('Y1MNF-SO')) = .TRUE.
                     AQROP(INDPRP('Z1MNF-SO')) = .TRUE.
                  ELSE IF (LABEL .EQ. 'FERMI CO') THEN
                     CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',
     &                           'UNFORMATTED',IDUMMY,.FALSE.)
 197                 CONTINUE
                     IF (NXTLAB(LABFND,RTNLBL,LUPROP)) THEN
                        IF (LABFND(1:3) .EQ. 'FC ')
     &                       AQROP(INDPRP(LABFND)) = .TRUE.
                        GOTO 197
                     END IF
                     CALL GPCLOSE(LUPROP,'KEEP')
                  ELSE IF (LABEL .EQ. 'SPIN-DIP') THEN
                     CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',
     &                           'UNFORMATTED',IDUMMY,.FALSE.)
 297                 CONTINUE
                     IF (NXTLAB(LABFND,RTNLBL,LUPROP)) THEN
                        IF (LABFND(1:3) .EQ. 'SD ')
     &                       AQROP(INDPRP(LABFND)) = .TRUE.
                        GOTO 297
                     END IF
                     CALL GPCLOSE(LUPROP,'KEEP')
                  ELSE IF (LABEL(1:3) .EQ. 'PSO' .OR.
     &                     LABEL(2:4) .EQ. 'PSO') THEN
                     CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',
     &                           'UNFORMATTED',IDUMMY,.FALSE.)
 397                 CONTINUE
                     IF (NXTLAB(LABFND,RTNLBL,LUPROP)) THEN
                        IF (LABFND(1:3) .EQ. 'PSO')
     &                       AQROP(INDPRP(LABFND)) = .TRUE.
                        GOTO 397
                     END IF
                     CALL GPCLOSE(LUPROP,'KEEP')
                  ELSE
                     AQROP( INDPRP(LABEL)) = .TRUE.
                     OPRANK(INDPRP(LABEL)) = IRANKA
                  END IF
               GO TO 100
    6          CONTINUE
                  READ (LUCMD,*) MAXITL
               GO TO 100
    7          CONTINUE
                  READ (LUCMD,*) THCLR
               GO TO 100
    8          CONTINUE
                  READ (LUCMD,*) MAXITO
               GO TO 100
    9          CONTINUE
                  READ (LUCMD,*) IPRHYP
               GO TO 100
   10          CONTINUE
                  LABEL='XDIPLEN'
                  AQROP( INDPRP(LABEL)) = .TRUE.
                  BQROP( INDPRP(LABEL)) = .TRUE.
                  CQROP( INDPRP(LABEL)) = .TRUE.
                  LABEL='YDIPLEN'
                  AQROP( INDPRP(LABEL)) = .TRUE.
                  BQROP( INDPRP(LABEL)) = .TRUE.
                  CQROP( INDPRP(LABEL)) = .TRUE.
                  LABEL='ZDIPLEN'
                  AQROP( INDPRP(LABEL)) = .TRUE.
                  BQROP( INDPRP(LABEL)) = .TRUE.
                  CQROP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
   11          CONTINUE
                  E3TEST = .TRUE.
               GO TO 100
   12          CONTINUE
               GO TO 100
   13          CONTINUE
                  READ (LUCMD,*) ISPINA,ISPINB,ISPINC
               GO TO 100
   14          CONTINUE
                  A2TEST = .TRUE.
               GO TO 100
   15          CONTINUE
                  X2TEST = .TRUE.
               GO TO 100
   16          CONTINUE
                  REFCHK = .TRUE.
               GO TO 100
   17          CONTINUE
                  READ(LUCMD,*) IAABB
               GO TO 100
   18          CONTINUE
                  QRSPEC = .TRUE.
                  QRSHG  = .TRUE.
               GO TO 100
   19          CONTINUE
                  QRSPEC = .TRUE.
                  QRPOCK = .TRUE.
               GO TO 100
   20          CONTINUE
                  QRSPEC = .TRUE.
                  QROPRF = .TRUE.
               GO TO 100
   21          CONTINUE
                  COMFREQ = .TRUE.
                  READ (LUCMD,*) NBQRFR
                  IF (NBQRFR.LE.MBQRFR) THEN
                     READ (LUCMD,*) (BQRFR(J),J=1,NBQRFR)
                  ELSE
                     WRITE (LUPRI,'(3(/,A,I5),/)')
     *               ' NUMBER OF FREQUENCIES SPECIFIED    :',NBQRFR,
     *               ' IS GREATER THAN THE ALLOWED NUMBER :',MBQRFR,
     *               ' THE NUMBER IS RESET TO THE MAXIMUM :',MBQRFR
                     READ (LUCMD,*) (BQRFR(J),J=1,MBQRFR),
     *                              (FFFF,J=MBQRFR+1,NBQRFR)
                     NBQRFR = MBQRFR
                  END IF
                  NCQRFR = NBQRFR
                  DO J = 1,NBQRFR
                     CQRFR(J) = BQRFR(J)
                  END DO
               GO TO 100
   22          CONTINUE
                  LABEL='XDIPLEN'
                  AQROP( INDPRP(LABEL)) = .TRUE.
                  BQROP( INDPRP(LABEL)) = .TRUE.
                  CQROP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
   23          CONTINUE
                  LABEL='YDIPLEN'
                  AQROP( INDPRP(LABEL)) = .TRUE.
                  BQROP( INDPRP(LABEL)) = .TRUE.
                  CQROP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
   24          CONTINUE
                  LABEL='ZDIPLEN'
                  AQROP( INDPRP(LABEL)) = .TRUE.
                  BQROP( INDPRP(LABEL)) = .TRUE.
                  CQROP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
   25          CONTINUE
               READ(LUCMD,*)ISPINA
               GO TO 100
   26          CONTINUE
               READ(LUCMD,*)ISPINB
               GO TO 100
   27          CONTINUE
               READ(LUCMD,*)ISPINC
               GO TO 100
 28            CONTINUE
                  SOCOLL = .TRUE.
               GO TO 100
 29            CONTINUE
                  SSCOLL = .TRUE.
               GO TO 100
 30            CONTINUE
                  QLOP = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' PROMPT "',WORD,
     *            '" NOT RECOGNIZED IN *QUADRA.'
               CALL QUIT(' ILLEGAL PROMPT IN *QUADRA')
            END IF
         GO TO 100
      END IF
  300 CONTINUE
      IF (THR_REDFAC .GT. 0.0D0) THEN
         ICHANG = ICHANG + 1
         WRITE (LUPRI,'(3A,1P,D10.2)') '@ INFO ',WORD1,
     &   ' thresholds multiplied with general factor',THR_REDFAC
         THCLR = THCLR*THR_REDFAC
      END IF
      NAQRTO = 0
      NBQRTO = 0
      NCQRTO = 0
      IF (ICHANG .GT. 0) THEN
         DO 500 I = 1,NPRLBL
            IF (AQROP(I)) THEN
                NAQRTO = NAQRTO + 1
                ! For singlet states couple operator rank to excitation rank
                IF (ISPIN.EQ.1) OPRANK(I) = ISPINA
            END IF
            IF (BQROP(I)) THEN
                NBQRTO = NBQRTO + 1
                IF (ISPIN.EQ.1) OPRANK(I) = ISPINB
            END IF
            IF (CQROP(I)) THEN
                NCQRTO = NCQRTO + 1
                IF (ISPIN.EQ.1) OPRANK(I) = ISPINC
            END IF
  500    CONTINUE
C
C  If special processes specified, we need to put the input freq. into
C  the relevant frequency arrays.
C
         IF (QRSPEC) THEN
            IF (BCFREQ) THEN
               WRITE (LUPRI,'(/A/A/A)')
     &  ' ERROR : inconsistent quadratic response input as',
     &  ' .BFREQ and/or .CFREQ  has been specified together with'//
     &     ' .SHG and/or .POCKEL',
     &  ' Use only .FREQUE with .SHG, .POCKEL or .OPTREF!'
               CALL QUIT('Inconsistent input for quadratic response')
            END IF
            NCQRFR = 0
C
            JFR = 0
            IF (QRPOCK) THEN
               NCQRFR = NCQRFR + 1
               CQRFR(NCQRFR) = ZERO
               DO IFR=1,NBQRFR
                  IF (ABS(BQRFR(IFR)).LE.THRFRQ) JFR = IFR
               ENDDO
            ENDIF
C
            IF (QRSHG) THEN
               DO IFR=1,NBQRFR
                  IF (IFR.NE.JFR) THEN
C                     IFR.EQ.JFR has already been included under QRPOCK
                     NCQRFR = NCQRFR + 1
                     CQRFR(NCQRFR) = BQRFR(IFR)
                  END IF
               END DO
            ENDIF
C
            IF (QROPRF) THEN
               DO IFR=1,NBQRFR
                  IF (IFR.NE.JFR) THEN
C                     IFR.EQ.JFR has already been included under QRPOCK
                     NCQRFR = NCQRFR + 1
                     CQRFR(NCQRFR) = -BQRFR(IFR)
                  END IF
               END DO
            ENDIF
C
         ENDIF
         IF (COMFREQ .AND. BCFREQ) THEN
            WRITE (LUPRI,'(/A/A/A)')
     &  ' *QUADRA ERROR : inconsistent quadratic response input as',
     &  ' .BFREQ and/or .CFREQ  has been specified together with'//
     &     ' .FREQUE',
     &  ' Use either .FREQUE or .BFREQ/.CFREQ !!!!!'
            CALL QUIT('Inconsistent input for quadratic response')
         END IF
      END IF
      NQRCAL = MIN(NAQRTO,NBQRTO,NCQRTO,NBQRFR,NCQRFR)
      IF  (HYPCAL .AND. NQRCAL.LE.0)  THEN
        HYPCAL = .FALSE.
        WRITE(LUPRI,'(/A)') ' Quadratic response input ignored because:'
        IF (NAQRTO .EQ. 0) WRITE (LUPRI,'(A)')
     *     ' - no A operators requested'
        IF (NBQRTO .EQ. 0) WRITE (LUPRI,'(A)')
     *     ' - no B operators requested'
        IF (NBQRFR .LE. 0) WRITE (LUPRI,'(A,I8)')
     &     ' - no B frequencies requested, NBQRFR =',NBQRFR
        IF (NCQRTO .EQ. 0) WRITE (LUPRI,'(A)')
     *     ' - no C operators requested'
        IF (NCQRFR .LE. 0) WRITE (LUPRI,'(A,I8)')
     &     ' - no C frequencies requested, NCQRFR =',NCQRFR
      ELSE IF (HYPCAL) THEN
         CALL HEADER('Quadratic Response calculation',0)
         WRITE (LUPRI,'(A,L2)')
     *      ' First hyperpolarizability calculation : HYPCAL='
     *      ,HYPCAL
         WRITE (LUPRI,'(/A,I5)') ' Spin of operator A , ISPINA=',ISPINA
         WRITE (LUPRI,'( A,I5)') ' Spin of operator B , ISPINB=',ISPINB
         WRITE (LUPRI,'( A,I5)') ' Spin of operator C , ISPINC=',ISPINC
         IF (ISPINA+ISPINB+ISPINC .GT. 0) THEN
            TRPLET = .TRUE.
            TRPFLG = .TRUE.
         END IF
         WRITE(LUPRI,'(/I3,A,(1P,5D14.6))')
     *      NBQRFR,' B-frequencies', (BQRFR(I),I=1,NBQRFR)
         WRITE(LUPRI,'(I3,A,(1P,5D14.6))')
     *      NCQRFR,' C-frequencies', (CQRFR(I),I=1,NCQRFR)
         IF (SOCOLL) WRITE (LUPRI,'(/A)')
     *        ' Quadratic response functions will be collected as '//
     *        'SO-corrections to nuclear shieldings'
         IF (SSCOLL) WRITE (LUPRI,'(/A)')
     *        ' Quadratic response functions will be collected as '//
     *        'SO-corrections to reduced spin-spin couplings'
         IF (FLAG(16) .AND. .NOT.INERSI) THEN
            WRITE(LUPRI,'(/A,L2)')
     &      ' Equilibrium solvent model requested            : SOLVNT ='
     &      ,FLAG(16)
            WRITE(LUPRI,'(A,F8.4)')
     &      '    Dielectric constant                         : EPSOL  ='
     &      ,EPSOL
         END IF
         IF (FLAG(16) .AND. INERSI) THEN
            WRITE(LUPRI,'(/A,L2)')
     &      ' Non-equilibrium solvent model requested        : INERSI ='
     &      ,INERSI
            WRITE(LUPRI,'(A,F8.4)')
     &      '    Static dielectric constant                  : EPSTAT ='
     &      ,EPSTAT
            WRITE(LUPRI,'(A,F8.4)')
     &      '    Optical dielectric constant                 : EPSOL  ='
     &      ,EPSOL
         END IF

         IF (QM3) WRITE(LUPRI,'(/A)')
     &      ' QM/MM response calculation, type .QM3'

         IF (QMMM) WRITE(LUPRI,'(/A)')
     &      ' QM/MM response calculation, type .QMMM'

         IF (PEQM) WRITE(LUPRI,'(/A/A)')
     &      ' Polarizable embedding response calculation, type .PEQM',
     &      ' Environment is modeled using'//
     &         ' the polarizable embedding scheme (PE library)'
         IF (QMNPMM) THEN
            WRITE(LUPRI,'(/A)') ' Heterogeneous environment is '// 
     &                          ' modeled using the QM/NP/MM scheme'
         END IF

         WRITE(LUPRI,'()')

         IF (QRSHG) WRITE (LUPRI,'(A)')
     &      ' Second harmonic generation;        SHG'
         IF (QRPOCK) WRITE (LUPRI,'(A)')
     &      ' Electro-optical Pockels effect;    POCKEL'
         IF (QROPRF) WRITE (LUPRI,'(A)')
     &      ' Optical refractivity calculation;  OPTREF'

         WRITE(LUPRI,'(/A,I4)')
     *      ' Print level                                    : IPRHYP ='
     *       ,IPRHYP
         WRITE(LUPRI,'(A,I4)')
     *      ' Maximum number of iterations in lin.rsp. solver: MAXITL ='
     *       ,MAXITL
         WRITE(LUPRI,'(A,1P,D10.3)')
     *      ' Threshold for convergence of linear resp. eq.s : THCLR  ='
     *      ,THCLR
         WRITE(LUPRI,'(A,I4)')
     *      ' Maximum iterations in optimal orbital algorithm: MAXITO ='
     *      ,MAXITO
         IF (DIROIT) WRITE (LUPRI,'(A,L2)')
     *      ' Direct one-index transformation                : DIROIT ='
     *      ,DIROIT
         IF (E3TEST) WRITE (LUPRI,'(/A)')
     *      ' Test of contributions from E3 and S3 terms.'
         IF (A2TEST) WRITE (LUPRI,'(/A)')
     *      ' Test of contributions from A2 terms.'
         IF (X2TEST) WRITE (LUPRI,'(/A)')
     *      ' Test of contributions from X2 terms.'
         IF (USE_PELIB()) THEN
            WRITE(LUPRI,'(A)') ' Environment is modeled using'//
     &                  ' the polarizable embedding scheme (PE library)'
         END IF
         IF (QMNPMM) THEN
            WRITE(LUPRI,'(/A)') ' Heterogeneous environment is '// 
     &                          ' modeled using the QM/NP/MM scheme'
         END IF
      END IF
C
C *** END OF HYPINP
C
      RETURN
      END
C  /* Deck smoinp */
      SUBROUTINE SMOINP(WORD)
C
C Sonia: added MCDBTE section (1999)
C
      use pelib_interface, only: use_pelib
#include "implicit.h"
C
#include "priunit.h"
#include "infrsp.h"
#include "inforb.h"
#include "rspprp.h"
#include "infsmo.h"
#include "inflr.h"
#include "infpp.h"
#include "infpri.h"
#include "infspi.h"
#include "inflin.h"
#include "infhso.h"
#include "maxorb.h"
#include "infinp.h"
#include "inftap.h"
! gnrinf.h: QM3, THR_REDFAC, PEQM
#include "gnrinf.h"
C
      LOGICAL NEWDEF
      PARAMETER ( NTABLE = 33)                                          !MK
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER*8 LABEL
C
      DATA TABLE /'.BPROP ', '.BFREQ ', '.APROP ', '.ROOTS ',
     *            '.MAXITP', '.THCLR ', '.MAXITL', '.THCPP ',
     *            '.MAXITO', '.PRINT ', '.DIPLEN', '.E3TEST',
     *            '.ISPABC', '.TPCD  ', '.X2TEST', '.A2TEST',
     *            '.DIPLNX', '.DIPLNY', '.DIPLNZ', '.FREQUE',
     *            '.SINGLE', '.PHOSPH', '.MNFPHO', '.TWO-PH',
     *            '.MCDBTE', '.MCDPRJ', '.ECPHOS', '.PHOSPV',           !MK
     *            '.DIPVEL', '.CPPHOL', '.CPPHOV', '.CPPHMF',           !MK
     *            '.CPPHEC'/                                            !MK
C
C READ IN  INPUT
C
      NEWDEF = (WORD .EQ. '*QUADRA')
      ICHANG = 0
      IF (NEWDEF) THEN
         SOMOM = .TRUE.
         JSPABC = 0
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') GO TO 100
            IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,
     &                 17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,
     &                 32,33), I                                        !MK
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' KEYWORD "',WORD,
     *            '" NOT RECOGNIZED IN SMOINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT(' ILLEGAL KEYWORD IN SMOINP ')
 1             CONTINUE
                  READ(LUCMD,'( A )')LABEL
                  BSMOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 2             CONTINUE
                  READ (LUCMD,*) NBSMFR
                  IF (NBSMFR.LE.MBSMFR) THEN
                     READ (LUCMD,*) (BSMFR(J),J=1,NBSMFR)
                  ELSE
                     WRITE (LUPRI,'(3(/,A,I5),/)')
     *               ' NUMBER OF B FREQUENCIES SPECIFIED  :',NBSMFR,
     *               ' IS GREATER THAN THE ALLOWED NUMBER :',MBSMFR,
     *               ' THE NUMBER IS RESET TO THE MAXIMUM :',MBSMFR
                     READ (LUCMD,*) (BSMFR(J),J=1,MBSMFR),
     *                              (FFFF,J=MBSMFR+1,NBSMFR)
                     NBSMFR = MBSMFR
                  END IF
               GO TO 100
 3             CONTINUE
                  READ(LUCMD,'( A )')LABEL
                  ASMOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 4             CONTINUE
                  READ (LUCMD,*) (NSMCNV(MULD2H(J,LSYMRF)),J=1,NSYM)
               GO TO 100
 5             CONTINUE
                  READ (LUCMD,*) MAXITL
               GO TO 100
 6             CONTINUE
                  READ (LUCMD,*) THCLR
               GO TO 100
 7             CONTINUE
                  READ (LUCMD,*) MAXITP
               GO TO 100
 8             CONTINUE
                  READ (LUCMD,*) THCPP
               GO TO 100
 9             CONTINUE
                  READ (LUCMD,*) MAXITO
               GO TO 100
 10            CONTINUE
                  READ (LUCMD,*) IPRSMO
               GO TO 100
 11            CONTINUE
                  LABEL='XDIPLEN'
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  BSMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='YDIPLEN'
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  BSMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='ZDIPLEN'
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  BSMOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 29            CONTINUE                                                 !MK
                  LABEL='XDIPVEL'                                       !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  BSMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='YDIPVEL'                                       !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  BSMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='ZDIPVEL'                                       !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  BSMOP( INDPRP(LABEL)) = .TRUE.                        !MK
               GO TO 100                                                !MK
 12            CONTINUE
                  E3TEST = .TRUE.
               GO TO 100
 13            CONTINUE
                  JSPABC = JSPABC + 1
                  READ (LUCMD,*) ISPINA,ISPINB,ISPINC
               GO TO 100
C              Option 14 .TPACD calculates tensor components for Tinoco
C              original, gauge invariant formulation
 14            CONTINUE
                  TWOPHO = .TRUE. 
                  LTPCD  = .TRUE.
C                 A operators
                  LABEL='XDIPVEL'
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='YDIPVEL'
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='ZDIPVEL'
                  ASMOP( INDPRP(LABEL)) = .TRUE.

                  LABEL='XANGMOM'
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='YANGMOM'
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='ZANGMOM'
                  ASMOP( INDPRP(LABEL)) = .TRUE.

                  LABEL='XXROTSTR'
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='YYROTSTR'
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='ZZROTSTR'
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='XYROTSTR'
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='XZROTSTR'
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='YZROTSTR'
                  ASMOP( INDPRP(LABEL)) = .TRUE.

C                 B operators

                  LABEL='XDIPVEL'
                  BSMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='YDIPVEL'
                  BSMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='ZDIPVEL'
                  BSMOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 15            CONTINUE
                  X2TEST = .TRUE.
               GO TO 100
 16            CONTINUE
                  A2TEST = .TRUE.
               GO TO 100
 17            CONTINUE
                  LABEL='XDIPLEN'
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  BSMOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 18            CONTINUE
                  LABEL='YDIPLEN'
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  BSMOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 19            CONTINUE
                  LABEL='ZDIPLEN'
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  BSMOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 20            CONTINUE
               GO TO 2
 21            CONTINUE
               GO TO 100
 22            CONTINUE
                  PHOSPH = .TRUE.
                  LABEL='XDIPLEN '
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='YDIPLEN '
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='ZDIPLEN '
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='X SPNORB'
                  BSMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='Y SPNORB'
                  BSMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='Z SPNORB'
                  BSMOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 28            CONTINUE                                                 !MK
                  PHOSPV = .TRUE.                                       !MK
                  LABEL='XDIPVEL '                                      !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='YDIPVEL '                                      !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='ZDIPVEL '                                      !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='X SPNORB'                                      !MK
                  BSMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='Y SPNORB'                                      !MK
                  BSMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='Z SPNORB'                                      !MK
                  BSMOP( INDPRP(LABEL)) = .TRUE.                        !MK
               GO TO 100                                                !MK
 30            CONTINUE                                                 !MK
                  CPPHOL = .TRUE.                                       !MK
                  LABEL='XDIPLEN '                                      !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='YDIPLEN '                                      !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='ZDIPLEN '                                      !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='X SPNORB'                                      !MK
                  BSMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='Y SPNORB'                                      !MK
                  BSMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='Z SPNORB'                                      !MK
                  BSMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='XANGMOM'                                       !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='YANGMOM'                                       !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='ZANGMOM'                                       !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
               GO TO 100                                                !MK
 31            CONTINUE                                                 !MK
                  CPPHOV = .TRUE.                                       !MK
                  LABEL='XDIPVEL '                                      !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='YDIPVEL '                                      !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='ZDIPVEL '                                      !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='X SPNORB'                                      !MK
                  BSMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='Y SPNORB'                                      !MK
                  BSMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='Z SPNORB'                                      !MK
                  BSMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='XANGMOM'                                       !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='YANGMOM'                                       !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='ZANGMOM'                                       !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
               GO TO 100                                                !MK
 32            CONTINUE                                                 !MK
                  CPPHMF = .TRUE.                                       !MK
                  LABEL='XDIPVEL '                                      !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='YDIPVEL '                                      !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='ZDIPVEL '                                      !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='XANGMOM'                                       !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='YANGMOM'                                       !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='ZANGMOM'                                       !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='X1MNF-SO'                                      !MK
                  BSMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='Y1MNF-SO'                                      !MK
                  BSMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='Z1MNF-SO'                                      !MK
                  BSMOP( INDPRP(LABEL)) = .TRUE.                        !MK
               GO TO 100                                                !MK
 33            CONTINUE                                                 !MK
                  CPPHEC = .TRUE.                                       !MK
                  LABEL='XDIPVEL '                                      !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='YDIPVEL '                                      !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='ZDIPVEL '                                      !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='X1SPNSCA'                                      !MK
                  BSMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='Y1SPNSCA'                                      !MK
                  BSMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='Z1SPNSCA'                                      !MK
                  BSMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='XANGMOM'                                       !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='YANGMOM'                                       !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
                  LABEL='ZANGMOM'                                       !MK
                  ASMOP( INDPRP(LABEL)) = .TRUE.                        !MK
               GO TO 100                                                !MK
 23            CONTINUE
                  MNFPHO = .TRUE.
                  LABEL='XDIPLEN '
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='YDIPLEN '
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='ZDIPLEN '
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='X1MNF-SO'
                  BSMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='Y1MNF-SO'
                  BSMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='Z1MNF-SO'
                  BSMOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 27            CONTINUE
                  ECPHOS = .TRUE.
                  LABEL='XDIPLEN '
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='YDIPLEN '
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='ZDIPLEN '
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='X1SPNSCA'
                  BSMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='Y1SPNSCA'
                  BSMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='Z1SPNSCA'
                  BSMOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 24            CONTINUE
                  TWOPHO = .TRUE.
                  LABEL='XDIPLEN'
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  BSMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='YDIPLEN'
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  BSMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='ZDIPLEN'
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  BSMOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 25            CONTINUE
                  MCDCAL=.TRUE.
                  LABEL='XDIPLEN'
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='YDIPLEN'
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='ZDIPLEN'
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='XANGMOM'
                  BSMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='YANGMOM'
                  BSMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='ZANGMOM'
                  BSMOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 26            CONTINUE
                  !sonia: mcd B term calculation
                  !with projection of the redundant exc state vector
                  MCDPRJ=.TRUE.
                  MCDCAL=.TRUE.
                  LABEL='XDIPLEN'
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='YDIPLEN'
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='ZDIPLEN'
                  ASMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='XANGMOM'
                  BSMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='YANGMOM'
                  BSMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='ZANGMOM'
                  BSMOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' PROMPT "',WORD,
     *            '" NOT RECOGNIZED IN SMOINP.'
               CALL QUIT(' ILLEGAL PROMPT IN SMOINP ')
            END IF
         GO TO 100
      END IF
  300 CONTINUE
      IF (THR_REDFAC .GT. 0.0D0) THEN
         ICHANG = ICHANG + 1
         WRITE (LUPRI,'(3A,1P,D10.2)') '@ INFO ',WORD1,
     &   ' thresholds multiplied with general factor',THR_REDFAC
         THCLR = THCLR*THR_REDFAC
         THCPP = THCPP*THR_REDFAC
      END IF
C
C If two-photon calculation then we assume that the frequencies of the
C photons match the excitation. At this point we set the frequencies
C to dummy values, later they become half of the excitation energies.
C
      IF (TWOPHO) THEN
         NSMCNV_TOT = ISUM(NSYM,NSMCNV,1)
         IF (NSMCNV_TOT .EQ. 0) THEN
            WRITE (LUPRI,'(/2A)')
     * ' Two-photon input ignored because no excited states specified.',
     * ' Use the .ROOTS keyword.'
            CALL QUIT('***** Input error for two-photon *****')
         END IF
         NBSMFR = NSMCNV_TOT
         DO J=1,NBSMFR
            BSMFR(J)=-1.0D0
         END DO
      END IF
C
      NASMTO = 0
      NBSMTO = 0
      IF (ICHANG .GT. 0) THEN
         DO 500 I = 1,NPRLBL
            IF (ASMOP(I)) NASMTO = NASMTO + 1
            IF (BSMOP(I)) NBSMTO = NBSMTO + 1
 500     CONTINUE
      END IF
      NSMCAL = MIN(NASMTO,NBSMTO,NBSMFR)
      IF  (SOMOM .AND. NSMCAL.LE.0)  THEN
         SOMOM = .FALSE.
         WRITE (LUPRI,'(/A)')
     & ' Quadratic Response single residue calculation ignored because:'
         IF (NASMTO .EQ. 0) WRITE (LUPRI,'(A)')
     &      ' - no A operators requested'
         IF (NBSMTO .EQ. 0) WRITE (LUPRI,'(A)')
     &      ' - no B operators requested'
         IF (NBSMFR .LE. 0) WRITE (LUPRI,'(A,I8)')
     &      ' - no B frequencies requested, NBSMFR =',NBSMFR
      ELSE IF (SOMOM) THEN
         CALL HEADER('Quadratic Response single residue calculation',0)
         WRITE(LUPRI,'(A)')
     &   ' That is, calculation of second order moments.'
         IF (PHOSPH) THEN
           WRITE(LUPRI,'(/A)') ' Phosphorescence calculation requested'
     &     //' with full spin-orbit operator.'
         END IF
         IF (MNFPHO) THEN
           WRITE(LUPRI,'(/A)') ' Phosphorescence calculation requested'
     &     //' with atomic mean-field spin-orbit operator (AMFI).'
         END IF
         IF (ECPHOS) THEN
           WRITE(LUPRI,'(/A)') ' Phosphorescence calculation requested'
     &     //' with effective charge spin-orbit operator.'
         END IF
         IF (PHOSPH .OR. MNFPHO .OR. ECPHOS .OR. PHOSPV .OR.            !MK
     &       CPPHOL .OR. CPPHOV .OR. CPPHMF .OR. CPPHEC) THEN           !MK
                 IF (JSPABC .GT. 0 .AND. ISPINA .NE. 0 .AND.
     &          ISPINB .NE. 1 .AND. ISPINC .NE. 1) THEN
               NWARN = NWARN + 1
               WRITE (LUPRI,'(/A/A,3I4/A)')
     &         ' WARNING : inconsistent phosphorescence input as',
     &         ' .ISPABC has been specified with',ISPINA,ISPINB,ISPINC,
     &         ' .ISPABC is reset to               0   1   1'
            END IF
            ISPINA = 0
            ISPINB = 1
            ISPINC = 1
C
            I = 0
            DO J = 1,NSYM
            IF (NSMCNV(J) .GT. MXPHOS) THEN
               I = I + 1
               NSMCNV(J) = MXPHOS
            END IF
            END DO
            IF (I .GT. 0) THEN
               NWARN = NWARN + 1
               WRITE (LUPRI,'(/A,I3,A/I2,A/)')
     &         ' WARNING: infhso.h is currently dimensioned to max.',
     &            MXPHOS,' excited states of any symmetry.',
     &         I,' specification(s) under .ROOTS have been'//
     &            ' reduced to this value.'
            END IF
         END IF
         IF (TWOPHO) THEN
            WRITE(LUPRI,'(/A/A)')
     &      ' Two-photon transition processes computed (.TWO-PHOTON).',
     &      ' Both photon frequencies are half the excitation energy.'
         ELSE
            WRITE(LUPRI,'(/I3,A,(1P,5D14.6))')
     &      NBSMFR,' B-frequencies', (BSMFR(I),I=1,NBSMFR)
         END IF
         WRITE (LUPRI,'(/A,I5)') ' Spin of operator A , ISPINA=',ISPINA
         WRITE (LUPRI,'( A,I5)') ' Spin of operator B , ISPINB=',ISPINB
         WRITE (LUPRI,'(2A,I5)') ' Spin of operator C ,'
     &    ,' (Excitation energy) ISPINC=',ISPINC
         IF (ISPINA+ISPINB+ISPINC .GT. 0) THEN
            TRPLET = .TRUE.
            TRPFLG = .TRUE.
         END IF
         IF (FLAG(16) .AND. .NOT.INERSI) THEN
            WRITE(LUPRI,'(/A,L2)')
     &      ' Equilibrium solvent model requested            : SOLVNT ='
     &      ,FLAG(16)
            WRITE(LUPRI,'(A,F8.4)')
     &      '    Dielectric constant                         : EPSOL  ='
     &      ,EPSOL
         END IF
         IF (FLAG(16) .AND. INERSI) THEN
            WRITE(LUPRI,'(/A,L2)')
     &      ' Non-equilibrium solvent model requested        : INERSI ='
     &      ,INERSI
            WRITE(LUPRI,'(A,F8.4)')
     &      '    Static dielectric constant                  : EPSTAT ='
     &      ,EPSTAT
            WRITE(LUPRI,'(A,F8.4)')
     &      '    Optical dielectric constant                 : EPSOL  ='
     &      ,EPSOL
         END IF

         IF (QM3) WRITE(LUPRI,'(/A)')
     &      ' QM/MM response calculation, type .QM3'

         IF (QMMM) WRITE(LUPRI,'(/A)')
     &      ' QM/MM response calculation, type .QMMM'

         IF (PEQM) WRITE(LUPRI,'(/A/A)')
     &      ' Polarizable embedding response calculation, type .PEQM',
     &      ' Environment is modeled using'//
     &         ' the polarizable embedding scheme (PE library)'

         IF (MCDCAL) THEN
            WRITE (LUPRI,'(A,L1)')
     *       ' B of Magnetic Circular Dichroism requested    : MCDCAL ='
     *      ,MCDCAL
            !sonia new
            if (MCDPRJ) then
               WRITE (LUPRI,'(A,L1)')
     *       ' Projected B of Magnetic Circular Dichroism    : MCDPRJ='
     *        ,MCDPRJ
            end if
         END IF
         WRITE(LUPRI,'(/A,I5)')
     *      ' Print level                                    : IPRSMO ='
     *      ,IPRSMO
         WRITE(LUPRI,'(A,1P,D10.3)')
     *      ' Threshold for convergence in RSPPP             : THCPP  ='
     *      ,THCPP
         WRITE(LUPRI,'(A,I5)')
     *      ' Maximum number of iterations in RSPPP          : MAXITP ='
     *      ,MAXITP
         WRITE(LUPRI,'(A,1P,D10.3)')
     *      ' Threshold for convergence in RSPLR             : THCLR  ='
     *      ,THCLR
         WRITE(LUPRI,'(A,I5)')
     *      ' Maximum number of iterations in RSPLR          : MAXITL ='
     *      ,MAXITL
         WRITE(LUPRI,'(A,I5)')
     *      ' Maximum iterations in optimal orbital algorithm: MAXITO ='
     *      ,MAXITO
         IF (E3TEST) WRITE (LUPRI,'(/A)')
     *      ' Test of contributions from E3 and S3 terms.'
         IF (USE_PELIB()) THEN
            WRITE(LUPRI,'(A)') ' Environment is modeled using'//
     &                  ' the polarizable embedding scheme (PE library)'
         END IF
      END IF
C
C *** END OF SMOINP
C
      RETURN
      END
C  /* Deck exiind */
      SUBROUTINE EXIIND(NEVAL)
C
C SET UP A LIST WITH THE NUMBER OF EIGENVECTORS THAT HAS
C TO BE SOLVED
C
#include "implicit.h"
C
#include "rspprp.h"
#include "indqr.h"
#include "inforb.h"
#include "infpri.h"
#include "infsmo.h"
C
      DIMENSION NEVAL(*)
C
C If two-photon calculation is specified we only compute a restricted
C number of excitation vectors.
C
      IF (TWOPHO) THEN
         DO ICSYM = 1,NSYM
         DO IBSYM = 1,NSYM
            IASYM = MULD2H(ICSYM,IBSYM)
            IF ( (NSMCNV(ICSYM) .GT. 0) .AND. (NBSMOP(IBSYM) .GT. 0)
     *          .AND. (NASMOP(IASYM).GT.0) ) THEN
               DO IC = 1,NSMCNV(ICSYM)
                  INUM = INQREX(ICSYM,IC)
               END DO
            END IF
         END DO
         END DO
      ELSE
C
C The general case computes all combinations.
C
         DO 100 ISYM = 1,NSYM
            DO 200 IROOT = 1,NEVAL(ISYM)
               INUM = INQREX(ISYM,IROOT)
 200        CONTINUE
 100     CONTINUE
C
      END IF
C
      RETURN
      END
C  /* Deck inqrex */
      INTEGER FUNCTION INQREX(ISYM,IROOT)
C
C FIND THE NUMBER THAT HAS BEEN GIVEN TO THE INUM EXCITATION OF
C SYMMETRY ISYM
C
#include "implicit.h"
C
#include "priunit.h"
#include "infrsp.h"
#include "rspprp.h"
#include "indqr.h"
#include "infpri.h"
C
      DO 100 I=1,NEXLBL
         IF ( (ISYM .EQ. ISEXQR(I) ) .AND. ( IROOT.EQ.JEXQR(I))) THEN
            INQREX = I
            RETURN
         END IF
 100  CONTINUE
      NEXLBL = NEXLBL + 1
      IF ( NEXLBL.GT.MXEXQR ) THEN
         WRITE(LUPRI,'(/A,I6)')
     &      'ERROR INQREX: Number of QR excitation operators exceeds'
     &      //' the allowed number, MXEXQR =',MXEXQR
         CALL QUIT(' INQREX: TOO MANY EXCITATION OPERATORS')
      END IF
      JEXQR(NEXLBL) = IROOT
      ISEXQR(NEXLBL) = ISYM
      INQREX = NEXLBL
      IF (IPRRSP.GT.10) THEN
         WRITE(LUPRI,'(/A)')' INQREX: NEW EXCITATION '
         WRITE(LUPRI,'(A,I5,2X,A,I5,2X,A,I5)')
     *   ' EXCITATION NUMBER:',NEXLBL,'SYMMETRY',ISYM,
     *   ' ROOT NUMBER',IROOT
      END IF
      RETURN
      END
C  /* Deck hypind */
      SUBROUTINE HYPIND
C
C CALCULATE A POINTER TO THE NUMBER OF DIFFERENT
C LINEAR RESPONSE EQUATIONS THAT NEED TO BE SOLVED
C IN A HYPERPOLARIZABILITY CALCULATION
C
#include "implicit.h"
C
#include "rspprp.h"
#include "infhyp.h"
#include "inforb.h"
#include "infrsp.h"
#include "infpri.h"
#include "infspi.h"
C
      DO 300 ICSYM = 1,NSYM
         DO 400 IBSYM = 1,NSYM
            IASYM = MULD2H(ICSYM,IBSYM)
            IF ( (NCQROP(ICSYM) .GT. 0) .AND. (NBQROP(IBSYM) .GT. 0)
     *          .AND. (NAQROP(IASYM).GT.0) ) THEN
               DO 500 ICOP = 1,NCQROP(ICSYM)
                  DO 600 ICFR = 1,NCQRFR
                     IF ( ISPINC.EQ.1 ) THEN
                        INUM = INQRTR(CQRLB(ICSYM,ICOP),
     *                                CQRFR(ICFR),ICSYM)
                     ELSE
                        INUM = INQRLR(CQRLB(ICSYM,ICOP),
     *                                CQRFR(ICFR),ICSYM)
                     END IF
 600              CONTINUE
 500           CONTINUE
               DO 700 IBOP = 1,NBQROP(IBSYM)
                  DO 800 IBFR = 1,NBQRFR
                     IF ( ISPINB.EQ.1 ) THEN
                        INUM = INQRTR(BQRLB(IBSYM,IBOP),
     *                                BQRFR(IBFR),IBSYM)
                     ELSE
                        INUM = INQRLR(BQRLB(IBSYM,IBOP),
     *                                BQRFR(IBFR),IBSYM)
                     END IF
 800              CONTINUE
 700           CONTINUE
C
C  If a special process is specified we only need some of the
C  combinations of the frequencies on the B and C operators
C
               IF (QRSPEC) THEN
                  DO IBFR=1,NBQRFR
C
C  If Pockels effect we need to solve linear response equations for w.
C
                     IF (QRPOCK) THEN
                        DO IAOP=1,NAQROP(IASYM)
                           INUM = INQRLR(AQRLB(IASYM,IAOP),
     &                                   BQRFR(IBFR),IASYM)
                        ENDDO
                     ENDIF
C
C  If SHG we need to solve linear response equations for 2w.
C
                     IF (QRSHG) THEN
                        DO IAOP=1,NAQROP(IASYM)
                           INUM = INQRLR(AQRLB(IASYM,IAOP),
     &                                 2*BQRFR(IBFR),IASYM)
                        ENDDO
                       END IF
                  ENDDO
C
C  If optical refractivity we need to solve linear response equations for 0.
C
                     IF (QROPRF) THEN
                        DO IAOP=1,NAQROP(IASYM)
                           INUM = INQRLR(AQRLB(IASYM,IAOP),
     &                                   0D0,IASYM)
                        ENDDO
                     ENDIF
               ELSE
                  DO 1000 ICFR = 1,NCQRFR
                     DO 1100 IBFR = 1,NBQRFR
                        BPCFR = CQRFR(ICFR) + BQRFR(IBFR)
                        DO 1200 IAOP=1,NAQROP(IASYM)
                           IF ( ISPINA.EQ.1 ) THEN
                              INUM = INQRTR(AQRLB(IASYM,IAOP),
     &                                      BPCFR,IASYM)
                           ELSE
                              INUM = INQRLR(AQRLB(IASYM,IAOP),
     &                                      BPCFR,IASYM)
                           END IF
 1200                   CONTINUE
 1100                CONTINUE
 1000             CONTINUE
               END IF
            END IF
 400     CONTINUE
 300  CONTINUE
      RETURN
      END
C  /* Deck somind */
      SUBROUTINE SOMIND
C
C CALCULATE A POINTER TO THE NUMBER OF DIFFERENT
C LINEAR RESPONSE EQUATIONS THAT ARE USED IN A
C A CALCULATION OF SECOND ORDER TRANSITION MOMENTS
C
#include "implicit.h"
C
#include "rspprp.h"
#include "infsmo.h"
#include "indqr.h"
#include "inforb.h"
#include "infrsp.h"
#include "infpri.h"
#include "infspi.h"
C
C If two-photon calculation and .BFREQ is not specified
C we only compute a restricted number of response vectors
C corresponding to BFREQ=CFREQ=EXCITA/2.
C
      IF (TWOPHO) THEN
         NBSMFR = 0
         DO ICSYM = 1,NSYM
         ICEXC = NSMCNV(ICSYM)
         DO IBSYM = 1,NSYM
            IASYM = MULD2H(ICSYM,IBSYM)
            IF ( (NSMCNV(ICSYM) .GT. 0) .AND. (NBSMOP(IBSYM) .GT. 0)
     *          .AND. (NASMOP(IASYM).GT.0) ) THEN
               DO IC = 1,NSMCNV(ICSYM)
                  IBFR = NBSMFR+IC
                  BSMFR(IBFR) = EXCITA(ICSYM,IC,ISPINC+1)/2
                  DO IAOP=1,NASMOP(IASYM)
                     INUM=INQRLR(ASMLB(IASYM,IAOP),BSMFR(IBFR),IASYM)
                  END DO
                  DO IBOP=1,NBSMOP(IBSYM)
C
C It is the response vector Nb(-wf/2) which is needed, but since if
C
C     Nb(wf/2)  = ( Z Y* ) then Nb(-wf/2) = ( Y* Z )
C
C so we only compute linear response vectors with positive frequencies.
C
                     INUM=INQRLR(BSMLB(IBSYM,IBOP),BSMFR(IBFR),IBSYM)
                  END DO
               END DO
               NBSMFR = NBSMFR + ICEXC
               ICEXC = 0
            END IF
         END DO
         END DO
      ELSE
C
C The general case computes all combinations.
C
      DO 300 ICSYM = 1,NSYM
         DO 400 IBSYM = 1,NSYM
            IASYM = MULD2H(ICSYM,IBSYM)
            IF ( (NSMCNV(ICSYM) .GT. 0) .AND. (NBSMOP(IBSYM) .GT. 0)
     *          .AND. (NASMOP(IASYM).GT.0) ) THEN
               DO 700 IBOP = 1,NBSMOP(IBSYM)
                  DO 800 IBFR = 1,NBSMFR
                     IF ( ISPINB.EQ.1 ) THEN
                        INUM = INQRTR(BSMLB(IBSYM,IBOP),
     *                                -BSMFR(IBFR),IBSYM)
                     ELSE
                        INUM = INQRLR(BSMLB(IBSYM,IBOP),
     *                                -BSMFR(IBFR),IBSYM)
                     END IF
 800              CONTINUE
 700           CONTINUE
               DO 900 IC = 1,NSMCNV(ICSYM)
                  DO 1000 IBFR = 1,NBSMFR
                     CEXMBF = EXCITA(ICSYM,IC,ISPINC+1) -BSMFR(IBFR)
                     DO 1100 IAOP = 1,NASMOP(IASYM)
                        IF ( ISPINA.EQ.0 ) THEN
                           INUM = INQRLR(ASMLB(IASYM,IAOP),CEXMBF,IASYM)
                        ELSE
                           INUM = INQRTR(ASMLB(IASYM,IAOP),CEXMBF,IASYM)
                        END IF
 1100                CONTINUE
 1000             CONTINUE
 900           CONTINUE
            END IF
 400     CONTINUE
 300  CONTINUE
C
      END IF
C
      RETURN
      END
C  /* Deck exmind */
      SUBROUTINE EXMIND
C
C CALCULATE A POINTER TO THE NUMBER OF DIFFERENT LINEAR
C RESPONSE EQUATIONS THAT ARE USED IN A CALCULATION OF
C TRANSITION MOMENTS BETWEEN EXCITED STATES
C
#include "implicit.h"
C
#include "infrsp.h"
#include "rspprp.h"
#include "infpp.h"
#include "infsmo.h"
#include "indqr.h"
#include "inforb.h"
#include "infpri.h"
#include "infspi.h"
#include "priunit.h"
C
C panor:
C     DOESA=.TRUE.
C     For excited state absorption we require the B-state to be the
C     first state of symmetry ESASYM and that the C-state is separate
C     from the B-state.
C
      DO 300 ICSYM = 1,NSYM
         IF (DOESA .OR. EXMTES) THEN
            NBSYME = NSYM
         ELSE
            NBSYME = ICSYM
         END IF
         DO 400 IBSYM = 1,NBSYME
            IF (DOESA.AND.IBSYM.NE.ESASYM) GOTO 400
            IASYM = MULD2H(ICSYM,IBSYM)
            IF ( (NPPCNV(ICSYM) .GT. 0) .AND. (NPPCNV(IBSYM) .GT. 0)
     *          .AND. (NGPPP(IASYM).GT.0) ) THEN
               DO 900 IC = 1,NPPCNV(ICSYM)
                  IF (DOESA .AND. ICSYM.EQ.IBSYM .AND. IC.EQ.1) GOTO 900
                  IF (DOESA) THEN
                     IBE = 1
                  ELSE IF ( EXMTES .OR. ICSYM.NE.IBSYM .OR. 
     *               ISPINB .NE. ISPINC) THEN
                     IBE = NPPCNV(IBSYM)
                  ELSE
                     IBE = IC
                  END IF
                  DO 1000 IB = 1,IBE
                     CEXMBE = EXCITA(ICSYM,IC,ISPINC+1) -
     *                        EXCITA(IBSYM,IB,ISPINB+1)
                     DO 1100 IAOP = 1,NGPPP(IASYM)
                        IF (.NOT.RSPCI.AND.
     *                       LCMEXC(IB,IBSYM).AND.LCMEXC(IC,ICSYM)) THEN
                           IF (ISPINA .EQ. 0) THEN
                              INUM = INQRLR
     *                               (LBLPP(IASYM,IAOP),CEXMBE,IASYM)
                           ELSE
                              INUM = INQRTR
     *                               (LBLPP(IASYM,IAOP),CEXMBE,IASYM)
                           END IF
                        END IF
 1100                CONTINUE
 1000             CONTINUE
 900           CONTINUE
            END IF
 400     CONTINUE
 300  CONTINUE
      RETURN
      END
C  /* Deck inqrlr */
      INTEGER FUNCTION INQRLR(NEWLBL,FRQVAL,ISYM)
C
C DETERMINE IF A POINTER TO LINEAR RESPONSE EQUATION
C CHARACTERIZED BY NEWLBL AND FRQVAL HAS OCCURRED
C IF SO INDQR GIVES ITS NUMBER IN THE LIST
C IF THE POINTER DOES NOT APPEAR IT IS ADDED TO
C LIST AND INDQR GIVES THE ADDED NUMBER
C
#include "implicit.h"
C
#include "priunit.h"
#include "infrsp.h"
#include "rspprp.h"
#include "indqr.h"
C
      PARAMETER ( DIFFR = 1.0D-8 )
      CHARACTER*8 NEWLBL
      DO 100 I=1,NLRLBL
         IF ( (NEWLBL .EQ. QRLBL(I) ) .AND.
     *       ( ABS(FRQVAL- QRFREQ(I)).LE.DIFFR) ) THEN
            INQRLR = I
            RETURN
         END IF
 100  CONTINUE
      NLRLBL = NLRLBL + 1
      IF (NLRLBL.GT.MXLRQR) THEN
         WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *   ' NUMBER OF SPECIFIED PROPERTIES EXCEED THE MAXIMUM ALLOWED',
     *   ' MXLRQR =',MXLRQR,' NLRLBL = ',NLRLBL
         CALL QUIT(' INQRLR: TOO MANY PROPERTIES SPECIFIED')
      END IF
      QRLBL(NLRLBL) = NEWLBL
      QRFREQ(NLRLBL) = FRQVAL
      ISYMQR(NLRLBL) = ISYM
      INQRLR = NLRLBL
      IF (IPRRSP.GT.40) THEN
         WRITE(LUPRI,'(/A)')' INQRLR: NEW OPERATOR'
         WRITE(LUPRI,'(A,I5,2X,A,A,2X,/A,D12.6,2X,A,I5)')
     *   ' NUMBER:',NLRLBL,' LABEL: ',NEWLBL,' FREQUENCY: ',FRQVAL,
     *   'SYMMETRY',ISYM
      END IF
      RETURN
      END
C  /* Deck qrlrsr */
      SUBROUTINE QRLRSR
C
C SORT THE LINEAR RESPONSE EQUATIONS AFTER SYMMETRY
C THE NUMBER OF LINEAR EQUATIONS OF A GIVEN SYMMETRY ARE
C STORED IN NLRQR WITH OFSET IN ILRQR
C
#include "implicit.h"
C
      CHARACTER*8 LABEL
C
#include "priunit.h"
#include "inforb.h"
#include "rspprp.h"
#include "indqr.h"
#include "infpri.h"
#include "infrsp.h"
C
      IF (IPRRSP.GT.40) THEN
         WRITE(LUPRI,'(/A,I5,2X,A)')
     *   ' LIST OF NLRLBL= ',NLRLBL,
     *   ' LINEAR EQUATIONS FOR QUADRATIC RESPONSE'
         WRITE(LUPRI,'(A)')' BEFORE SORTING'
         WRITE(LUPRI,'(/A)')'   I  QRLBL(I)  QRFREQ(I)  ISYMQR(I) '
         DO 90 I = 1,NLRLBL
            WRITE(LUPRI,'(I5,5X,A,D12.6,I12)')
     *      I,QRLBL(I),QRFREQ(I),ISYMQR(I)
  90     CONTINUE
      END IF
      ISMMAX = 1
      DO ISYM = 1,NSYM
         DO ISRT = 1, NLRLBL
            IF (ISRT.GE.ISMMAX .AND. ISYMQR(ISRT).EQ.ISYM) THEN
               IQR   = ISYMQR(ISRT)
               FRVAL = QRFREQ(ISRT)
               LABEL = QRLBL(ISRT)
               ISYMQR(ISRT) = ISYMQR(ISMMAX)
               QRFREQ(ISRT) = QRFREQ(ISMMAX)
               QRLBL(ISRT)  = QRLBL(ISMMAX)
               ISYMQR(ISMMAX) = IQR
               QRFREQ(ISMMAX) = FRVAL
               QRLBL(ISMMAX)  = LABEL
               ISMMAX = ISMMAX + 1
            END IF
         END DO
      END DO
      ITOTOP = 0
      DO 300 ISYM = 1,NSYM
         INUMOP = 0
         DO 400 I = 1,NLRLBL
            IF(ISYMQR(I).EQ.ISYM) THEN
               INUMOP = INUMOP + 1
            END IF
 400     CONTINUE
         ILRQR(ISYM) = ITOTOP
         ITOTOP = ITOTOP + INUMOP
         NLRQR(ISYM) = INUMOP
 300  CONTINUE
      IF (IPRRSP.GT.40) THEN
         WRITE(LUPRI,'(/A,I5,2X,A)')
     *   ' LIST OF NLRLBL= ',NLRLBL,
     *   ' LINEAR EQUATIONS FOR QUADRATIC RESPONSE'
         WRITE(LUPRI,'(A)')'AFTER SORTING'
         DO 301 I = 1,NLRLBL
            WRITE(LUPRI,'(/A)')'   I  QRLBL(I)  QRFREQ(I)  ISYMQR(I) '
            WRITE(LUPRI,'(I5,5X,A,D12.6,I12)')
     *      I,QRLBL(I),QRFREQ(I),ISYMQR(I)
 301     CONTINUE
      END IF
      IF (IPRRSP.GT.35) THEN
         WRITE(LUPRI,'(/A)')
     *   ' NUMBER OF LINEAR RESPONSE EQUATIONS IN VARIOUS SYMMETRIES'
         WRITE(LUPRI,'(/A)')' NLRQR(ISYM),ISYM=1,NSYM '
         WRITE(LUPRI,'(8I5)')( NLRQR(ISYM),ISYM=1,NSYM )
      END IF
      RETURN
      END
C  /* Deck qrexsr */
      SUBROUTINE QREXSR
C
C SORT THE EXCITATION ENERGIES  AFTER SYMMETRY
C THE NUMBER OF EXCITATION ENERGIES OF A GIVEN SYMMETRY ARE
C STORED IN NEXQR WITH OFFSET IN IEXQR
C
#include "implicit.h"
C
#include "priunit.h"
#include "infrsp.h"
#include "inforb.h"
#include "rspprp.h"
#include "indqr.h"
#include "infpri.h"
C
      IF (IPRRSP.GT.40) THEN
         WRITE(LUPRI,'(/A,I5,2X,A)')
     *   ' LIST OF NEXLBL= ',NEXLBL,'EXCITATIONS FOR QUADRATIC RESPONSE'
         WRITE(LUPRI,'(A)')' BEFORE SORTING'
         DO 90 I = 1,NEXLBL
            WRITE(LUPRI,'(/A)')'     I   ISEXQR(I)  JEXQR(I)'
            WRITE(LUPRI,'(/3I10)')I,ISEXQR(I),JEXQR(I)
 90      CONTINUE
      END IF
      ISMMAX = 1
      DO 100 ISYM = 1,NSYM
         I = ISMMAX
         DO 200 ISRT = I,NEXLBL
            IF ( ISEXQR(ISRT).EQ.ISYM) THEN
                IQR   = ISEXQR(ISRT)
                JQR   = JEXQR(ISRT)
                ISEXQR(ISRT) = ISEXQR(ISMMAX)
                JEXQR(ISRT)  = JEXQR(ISMMAX)
                ISEXQR(ISMMAX) = IQR
                JEXQR(ISMMAX)  = JQR
                ISMMAX = ISMMAX + 1
            END IF
 200     CONTINUE
 100  CONTINUE
      ITOTOP = 0
      DO 300 ISYM = 1,NSYM
         INUMOP = 0
         DO 400 I = 1,NEXLBL
            IF(ISEXQR(I).EQ.ISYM) THEN
               INUMOP = INUMOP + 1
            END IF
 400     CONTINUE
         IEXQR(ISYM) = ITOTOP
         ITOTOP = ITOTOP + INUMOP
         NEXQR(ISYM) = INUMOP
 300  CONTINUE
      IF (IPRRSP.GT.40) THEN
         WRITE(LUPRI,'(/A,I5,2X,A)')
     *   ' LIST OF NEXLBL= ',NEXLBL,'EXCITATIONS FOR QUADRATIC RESPONSE'
         WRITE(LUPRI,'(A)')' AFTER SORTING'
         DO 310 I = 1,NEXLBL
            WRITE(LUPRI,'(/A)')'     I   ISEXQR(I)  JEXQR(I)'
            WRITE(LUPRI,'(/3I10)')I,ISEXQR(I),JEXQR(I)
 310     CONTINUE
      END IF
      IF (IPRRSP.GT.35) THEN
         WRITE(LUPRI,'(/A)')
     *   ' NUMBER OF EXCITATIONS IN VARIOUS SYMMETRIES'
         WRITE(LUPRI,'(/A)')' NEXQR(ISYM),ISYM=1,NSYM '
         WRITE(LUPRI,'(8I5)')( NEXQR(ISYM),ISYM=1,NSYM )
      END IF
      RETURN
      END
C  /* Deck inqrtr */
      INTEGER FUNCTION INQRTR(NEWLBL,FRQVAL,ISYM)
C
C DETERMINE IF A POINTER TO TRIPLET LINEAR RESPONSE EQUATION
C CHARACTERIZED BY NEWLBL AND FRQVAL
C HAS OCCURED IF SO INQRTR GIVES ITS NUMBER IN THE
C LIST IF THE POINTER DOES NOT APPEAR IT IS ADDED TO
C LIST AND INQRTR GIVES THE ADDED NUMBER
C
#include "implicit.h"
C
#include "priunit.h"
#include "infrsp.h"
#include "rspprp.h"
#include "indqr.h"
#include "infpri.h"
C
      PARAMETER ( DIFFR = 1.0D-8 )
      CHARACTER*8 NEWLBL
      DO 100 I=1,NTRLBL
         IF ( (NEWLBL .EQ. TRLBL(I) ) .AND.
     *       ( ABS(FRQVAL- TRFREQ(I)).LE.DIFFR) ) THEN
            INQRTR = I
            RETURN
         END IF
 100  CONTINUE
      NTRLBL = NTRLBL + 1
      IF (NTRLBL.GT.MXLRQR) THEN
         WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *   ' NUMBER OF SPECIFIED PROPERTIES EXCEED THE MAXIMUM ALLOWED',
     *   ' MXLRQR =',MXLRQR,' NTRLBL = ',NTRLBL
         CALL QUIT(' INQRTR: TOO MANY PROPERTIES SPECIFIED')
      END IF
      TRLBL(NTRLBL) = NEWLBL
      TRFREQ(NTRLBL) = FRQVAL
      ISYMTR(NTRLBL) = ISYM
      INQRTR = NTRLBL
      IF (IPRRSP.GT.40) THEN
         WRITE(LUPRI,'(/A)')' INQRTR: NEW OPERATOR'
         WRITE(LUPRI,'(A,I5,2X,A,A,2X,/A,D12.6,2X,A,I5)')
     *   ' NUMBER:',NTRLBL,' LABEL: ',NEWLBL,' FREQUENCY: ',FRQVAL,
     *   'SYMMETRY',ISYM
      END IF
      RETURN
      END
C  /* Deck qrtrsr */
      SUBROUTINE QRTRSR
C
C SORT THE TRIPLET LINEAR RESPONSE EQUATIONS AFTER SYMMETRY
C THE NUMBER OF LINEAR EQUATIONS OF A GIVEN SYMMETRY ARE
C STORED IN NTRQR WITH OFSET IN ITRQR
C
#include "implicit.h"
C
      CHARACTER*8 LABEL
C
#include "priunit.h"
#include "inforb.h"
#include "rspprp.h"
#include "indqr.h"
#include "infpri.h"
#include "infrsp.h"
C
      IF (IPRRSP.GT.40) THEN
         WRITE(LUPRI,'(/A,I5,2X,A)')
     *   ' LIST OF NTRLBL= ',NTRLBL,
     *   ' TRIPLET LINEAR EQUATIONS FOR QUADRATIC RESPONSE'
         WRITE(LUPRI,'(A)')' BEFORE SORTING'
         WRITE(LUPRI,'(/A)')'   I  TRLBL(I)  TRFREQ(I)  ISYMTR(I) '
         DO 90 I = 1,NTRLBL
            WRITE(LUPRI,'(I5,5X,A,D12.6,I12)')
     *      I,TRLBL(I),TRFREQ(I),ISYMTR(I)
  90     CONTINUE
      END IF
      ISMMAX = 1
      DO 100 ISYM = 1,NSYM
         I = ISMMAX
         DO 200 ISRT = I,NTRLBL
            IF ( ISYMTR(ISRT).EQ.ISYM) THEN
                IQR   = ISYMTR(ISRT)
                FRVAL  = TRFREQ(ISRT)
                LABEL = TRLBL(ISRT)
                ISYMTR(ISRT) = ISYMTR(ISMMAX)
                TRFREQ(ISRT) = TRFREQ(ISMMAX)
                TRLBL(ISRT)  = TRLBL(ISMMAX)
                ISYMTR(ISMMAX) = IQR
                TRFREQ(ISMMAX) = FRVAL
                TRLBL(ISMMAX)  = LABEL
                ISMMAX = ISMMAX + 1
            END IF
 200     CONTINUE
 100  CONTINUE
      ITOTOP = 0
      DO 300 ISYM = 1,NSYM
         INUMOP = 0
         DO 400 I = 1,NTRLBL
            IF(ISYMTR(I).EQ.ISYM) THEN
               INUMOP = INUMOP + 1
            END IF
 400     CONTINUE
         ITRQR(ISYM) = ITOTOP
         ITOTOP = ITOTOP + INUMOP
         NTRQR(ISYM) = INUMOP
 300  CONTINUE
      IF (IPRRSP.GT.40) THEN
         WRITE(LUPRI,'(/A,I5,2X,A/A)')
     *   ' LIST OF NTRLBL= ',NTRLBL,
     *   ' TRIPLET LINEAR EQUATIONS FOR QUADRATIC RESPONSE',
     *   ' AFTER SORTING'
         DO 301 I = 1,NTRLBL
            WRITE(LUPRI,'(/A)')'   I  TRLBL(I)  TRFREQ(I)  ISYMTR(I) '
            WRITE(LUPRI,'(I5,5X,A,D12.6,I12)')
     *      I,TRLBL(I),TRFREQ(I),ISYMTR(I)
 301     CONTINUE
      END IF
      IF (IPRRSP.GT.35) THEN
         WRITE(LUPRI,'(/2A)')
     *   ' NUMBER OF TRIPLET LINEAR RESPONSE EQUATIONS',
     *   ' IN VARIOUS SYMMETRIES'
         WRITE(LUPRI,'(/A)')' NTRQR(ISYM),ISYM=1,NSYM '
         WRITE(LUPRI,'(8I5)')( NTRQR(ISYM),ISYM=1,NSYM )
      END IF
      RETURN
      END
