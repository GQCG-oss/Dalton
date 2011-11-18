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
      SUBROUTINE TMOINP(WORD)
C
#include <implicit.h>
C
#include <priunit.h>
#include <gnrinf.h>
#include <infrsp.h>
#include <inforb.h>
#include "rspprp.h"
#include <inftmo.h>
#include <inflr.h>
#include <infpp.h>
#include <infpri.h>
#include <infspi.h>
#include <infcr.h>
#include <inflin.h>
#include <inftap.h>
C
      LOGICAL NEWDEF
      PARAMETER ( NTABLE = 20)
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER*8 LABEL
C
      DATA TABLE /'.CPROP ', '.CFREQ ', '.BPROP ', '.BFREQ ',
     *            '.APROP ', '.ROOTS ', '.MAXITP', '.THCPP ', 
     *            '.MAX IT', '.THCLR ', '.MAXITO', '.PRINT ', 
     *            '.NOHG  ', '.DIPLEN', '.DIPLNX', '.DIPLNY', 
     *            '.DIPLNZ', '.SINGLE', '.FREQUE', '.THREE-'/
C
C READ IN  INPUT
C
      NEWDEF = (WORD .EQ. '*CUBIC ')
      ICHANG = 0
      IF (NEWDEF) THEN
         TOMOM = .TRUE.
         WORD1 = WORD
         DO J=1,NSYM
            NTMCNV(J)=1
         END DO
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') GO TO 100
            IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO I=1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,
     &                      15,16,17,18,19,20), I
                  END IF
               END DO
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' KEYWORD "',WORD,
     *            '" NOT RECOGNIZED IN TMOINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT(' ILLEGAL KEYWORD IN TMOINP ')
 1             CONTINUE
                  READ(LUCMD,'( A )')LABEL
                  CTMOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 2             CONTINUE
                  READ (LUCMD,*) NCTMFR
                  IF (NCTMFR.LE.MCTMFR) THEN
                     READ (LUCMD,*) (CTMFR(J),J=1,NCTMFR)
                  ELSE
                     WRITE (LUPRI,'(3(/A,I5),/)')
     *        'WARNING! NUMBER OF C-FREQUENCIES SPECIFIED  :',NCTMFR,
     *        'WARNING! IS GREATER THAN THE ALLOWED NUMBER :',MCTMFR,
     *        'WARNING! THE NUMBER IS RESET TO THE MAXIMUM :',MCTMFR
                     READ (LUCMD,*) (CTMFR(J),J=1,MCTMFR),
     *                              (FFFF,J=MCTMFR+1,NCTMFR)
                     NCTMFR = MCTMFR
                  END IF
               GO TO 100
 3             CONTINUE
                  READ(LUCMD,'( A )')LABEL
                  BTMOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 4             CONTINUE
                  READ (LUCMD,*) NBTMFR
                  IF (NBTMFR.LE.MBTMFR) THEN
                     READ (LUCMD,*) (BTMFR(J),J=1,NBTMFR)
                  ELSE
                     WRITE (LUPRI,'(3(/A,I5),/)')
     *        'WARNING! NUMBER OF B-FREQUENCIES SPECIFIED  :',NBTMFR,
     *        'WARNING! IS GREATER THAN THE ALLOWED NUMBER :',MBTMFR,
     *        'WARNING! THE NUMBER IS RESET TO THE MAXIMUM :',MBTMFR
                     READ (LUCMD,*) (BTMFR(J),J=1,MBTMFR),
     *                              (FFFF,J=MBTMFR+1,NBTMFR)
                     NBTMFR = MBTMFR
                  END IF
               GO TO 100
 5             CONTINUE
                  READ(LUCMD,'( A )')LABEL
                  ATMOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 6             CONTINUE
               READ (LUCMD,*) (NTMCNV(MULD2H(J,LSYMRF)),J=1,NSYM)
               GO TO 100
 7             CONTINUE
                  READ (LUCMD,*) MAXITP
               GO TO 100
 8             CONTINUE
                  READ (LUCMD,*) THCPP
               GO TO 100
 9             CONTINUE
                  READ (LUCMD,*) MAXITL
               GO TO 100
 10            CONTINUE
                  READ (LUCMD,*) THCLR
               GO TO 100
 11            CONTINUE
                  READ (LUCMD,*) MAXITO
               GO TO 100
 12            CONTINUE
                  READ (LUCMD,*) IPRTMO
               GO TO 100
 13            CONTINUE
                  CTMOHG = .FALSE.
               GO TO 100
 14            CONTINUE
                  LABEL='XDIPLEN'
                  ATMOP( INDPRP(LABEL)) = .TRUE.
                  BTMOP( INDPRP(LABEL)) = .TRUE.
                  CTMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='YDIPLEN'
                  ATMOP( INDPRP(LABEL)) = .TRUE.
                  BTMOP( INDPRP(LABEL)) = .TRUE.
                  CTMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='ZDIPLEN'
                  ATMOP( INDPRP(LABEL)) = .TRUE.
                  BTMOP( INDPRP(LABEL)) = .TRUE.
                  CTMOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 15            CONTINUE
                  LABEL='XDIPLEN'
                  ATMOP( INDPRP(LABEL)) = .TRUE.
                  BTMOP( INDPRP(LABEL)) = .TRUE.
                  CTMOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 16            CONTINUE
                  LABEL='YDIPLEN'
                  ATMOP( INDPRP(LABEL)) = .TRUE.
                  BTMOP( INDPRP(LABEL)) = .TRUE.
                  CTMOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 17            CONTINUE
                  LABEL='ZDIPLEN'
                  ATMOP( INDPRP(LABEL)) = .TRUE.
                  BTMOP( INDPRP(LABEL)) = .TRUE.
                  CTMOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 18            CONTINUE
               GO TO 100
 19            CONTINUE
                  READ (LUCMD,*) NBTMFR
                  IF (NBTMFR.LE.MBTMFR) THEN
                     READ (LUCMD,*) (BTMFR(J),J=1,NBTMFR)
                  ELSE
                     WRITE (LUPRI,'(3(/,A,I5),/)')
     *          'WARNING! NUMBER OF B-FREQUENCIES SPECIFIED  :',NBTMFR,
     *          'WARNING! IS GREATER THAN THE ALLOWED NUMBER :',MBTMFR,
     *          'WARNING! THE NUMBER IS RESET TO THE MAXIMUM :',MBTMFR
                     READ (LUCMD,*) (BTMFR(J),J=1,MBTMFR),
     *                              (FFFF,J=MBTMFR+1,NBTMFR)
                     NBTMFR = MBTMFR
                  END IF
                  NCTMFR=NBTMFR
                  DO J=1,NBTMFR
                     CTMFR(J)=BTMFR(J)
                  END DO
               GO TO 100         
 20            CONTINUE
                  THREEPHOTON = .TRUE.
                  LABEL='XDIPLEN'
                  ATMOP( INDPRP(LABEL)) = .TRUE.
                  BTMOP( INDPRP(LABEL)) = .TRUE.
                  CTMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='YDIPLEN'
                  ATMOP( INDPRP(LABEL)) = .TRUE.
                  BTMOP( INDPRP(LABEL)) = .TRUE.
                  CTMOP( INDPRP(LABEL)) = .TRUE.
                  LABEL='ZDIPLEN'
                  ATMOP( INDPRP(LABEL)) = .TRUE.
                  BTMOP( INDPRP(LABEL)) = .TRUE.
                  CTMOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' PROMPT "',WORD,
     *            '" NOT RECOGNIZED IN TMOINP.'
               CALL QUIT(' ILLEGAL PROMPT IN TMOINP ')
            END IF
         GO TO 100
      END IF
  300 CONTINUE
      IF (THR_REDFAC .GT. 0.0D0) THEN
         ICHANG = ICHANG + 1
         WRITE (LUPRI,'(3A,1P,D10.2)') '@ INFO ',WORD1,
     &   ' thresholds multiplied with general factor',THR_REDFAC
         THCPP = THCPP*THR_REDFAC
         THCLR = THCLR*THR_REDFAC
      END IF
      NATMTO = 0
      NBTMTO = 0
      NCTMTO = 0
      IF (ICHANG .GT. 0) THEN
         DO 500 I = 1,NPRLBL
            IF (ATMOP(I)) NATMTO = NATMTO + 1
            IF (BTMOP(I)) NBTMTO = NBTMTO + 1
            IF (CTMOP(I)) NCTMTO = NCTMTO + 1
 500     CONTINUE
         IF (NATMTO .EQ. 0) WRITE (LUPRI,'(/A)')
     * 'WARNING! *CRTMO input ignored because no A operators requested.'
         IF (NBTMTO .EQ. 0) WRITE (LUPRI,'(/A)')
     * 'WARNING! *CRTMO input ignored because no B operators requested.'
         IF (NCTMTO .EQ. 0) WRITE (LUPRI,'(/A)')
     * 'WARNING! *CRTMO input ignored because no C operators requested.'
      END IF
C
C If three-photon calculation then we assume that the frequencies of the 
C photons match the excitation. At this point we set the frequencies
C to dummy values, later they become one third of the excitation energies.
C
      IF (THREEPHOTON) THEN
         NBTMFR = ISUM(NSYM,NTMCNV,1)
         NCTMFR = NBTMFR
         IF (NBTMFR .EQ. 0) THEN
            WRITE (LUPRI,'(/2A)')
     *      ' Input ignored because no excited states specified.',
     *      ' Use the .ROOTS keyword.'
            WRITE(LUPRI,'(/A)') ' <<<<<<< End of DALTON >>>>>>>'
            CALL QUIT('***** End of DALTON *****')
         END IF
         DO J=1,NBTMFR
            BTMFR(J)=1.0D0
            CTMFR(J)=1.0D0
         END DO
      END IF
C
      NTMCAL = MIN(NATMTO,NBTMTO,NCTMTO,NBTMFR,NCTMFR)
      IF  (NTMCAL.GT.0)  THEN
         CALL HEADER('Cubic Response single residue calculation',0)
         WRITE (LUPRI,'(A,L1)')
     *      ' Third moments calculated : TOMOM='
     *      ,TOMOM
         IF (THREEPHOTON) THEN
            WRITE(LUPRI,'(A,L1,/)')
     *      ' Three-photon transition process computed  : THREEPHOTON ='
     *      ,THREEPHOTON
            WRITE(LUPRI,'(/I3,A)')
     *      NBTMFR,' B- and C-frequencies will later be set to 1/3 '//
     *      'of the excitation energies'
         ELSE
            WRITE(LUPRI,'(/I3,A,(T18,1P,5D14.6))')
     *      NBTMFR,' B-frequencies', (BTMFR(I),I=1,NBTMFR)
            WRITE(LUPRI,'(I3,A,(T18,1P,5D14.6))')
     *      NCTMFR,' C-frequencies', (CTMFR(I),I=1,NCTMFR)
         END IF
         WRITE(LUPRI,'(/A,I5)')
     *      ' Print level                                    : IPRTMO ='
     *      ,IPRTMO
         WRITE(LUPRI,'(A,I5)')
     *      ' Maximum number of exc.en. iterations           : MAXITP ='
     *      ,MAXITP
         WRITE(LUPRI,'(A,1P,D10.3)')
     *      ' Threshold for exc.en. convergence              : THCPP  ='
     *      ,THCPP
         WRITE(LUPRI,'(A,I5)')
     *      ' Maximum number of lin.resp. iterations         : MAXITL ='
     *      ,MAXITL
         WRITE(LUPRI,'(A,1P,D10.3)')
     *      ' Threshold for lin.resp. convergence            : THCLR  ='
     *      ,THCLR
         WRITE(LUPRI,'(A,I5)')
     *      ' Maximum iterations in optimal orbital algorithm: MAXITO ='
     *      ,MAXITO
      END IF
C
C *** END OF TMOINP
C
      RETURN
      END
