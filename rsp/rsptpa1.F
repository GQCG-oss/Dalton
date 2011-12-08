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
      SUBROUTINE TPAINP(WORD)
C
#include <implicit.h>
C
#include <priunit.h>
#include <gnrinf.h>
#include <infrsp.h>
#include <inforb.h>
#include <rspprp.h>
#include <inftpa.h>
#include <inflr.h>
#include <infpp.h>
#include <infpri.h>
#include <infspi.h>
#include <infcr.h>
#include <inflin.h>
#include <inftap.h>
C
      LOGICAL NEWDEF
      PARAMETER ( NTABLE = 16)
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER*8 LABEL
C
      DATA TABLE /'.BPROP ', '.BFREQ ', '.APROP ', '.ROOTS ',
     *            '.MAXITP', '.THCPP ', '.MAX IT', '.THCLR ',
     *            '.MAXITO', '.PRINT ', '.FREQUE', '.DIPLEN', 
     *            '.DIPLNX', '.DIPLNY', '.DIPLNZ', '.DOUBLE'/
C
C READ IN  INPUT
C
      NEWDEF = (WORD .EQ. '*CUBIC ')
      ICHANG = 0
      IF (NEWDEF) THEN
         TPAMP = .TRUE.
         WORD1 = WORD
         DO J=1,NSYM
            NTPCN2(J)=1
            NTPCN1(J)=1
         END DO
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') GO TO 100
            IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' KEYWORD "',WORD,
     *            '" NOT RECOGNIZED IN TPAINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT(' ILLEGAL KEYWORD IN TPAINP ')
 1             CONTINUE
                  READ(LUCMD,'( A )')LABEL
                  BTPOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 2             CONTINUE
                  READ (LUCMD,*) NBTPFR
                  IF (NBTPFR.LE.MBTPFR) THEN
                     READ (LUCMD,*) (BTPFR(J),J=1,NBTPFR)
                  ELSE
                     WRITE (LUPRI,'(3(/,A,I5),/)')
     *               ' NUMBER OF FREQUENCIES SPECIFIED    :',NBTPFR,
     *               ' IS GREATER THAN THE ALLOWED NUMBER :',MBTPFR,
     *               ' THE NUMBER IS RESET TO THE MAXIMUM :',MBTPFR
                     READ (LUCMD,*) (BTPFR(J),J=1,MBTPFR),
     *                              (FFFF,J=MBTPFR+1,NBTPFR)
                     NBTPFR = MBTPFR
                  END IF
               GO TO 100
 3             CONTINUE
                  READ(LUCMD,'( A )')LABEL
                  ATPOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 4             CONTINUE
               READ (LUCMD,*) (NTPCN1(MULD2H(J,LSYMRF)),J=1,NSYM)
               DO J=1,NSYM
                  NTPCN2(J)=NTPCN1(J)
               END DO
               GO TO 100
 5             CONTINUE
                  READ (LUCMD,*) MAXITP
               GO TO 100
 6             CONTINUE
                  READ (LUCMD,*) THCPP
               GO TO 100
 7             CONTINUE
                  READ (LUCMD,*) MAXITL
               GO TO 100
 8             CONTINUE
                  READ (LUCMD,*) THCLR
               GO TO 100
 9             CONTINUE
                  READ (LUCMD,*) MAXITO
               GO TO 100
 10            CONTINUE
                  READ (LUCMD,*) IPRTPA
               GO TO 100
 11            CONTINUE
               GO TO 2
 12            CONTINUE
                  TPALP = .TRUE.
                  LABEL='XDIPLEN'
                  ATPOP(INDPRP(LABEL)) = .TRUE.
                  BTPOP(INDPRP(LABEL)) = .TRUE.
                  LABEL='YDIPLEN'
                  ATPOP(INDPRP(LABEL)) = .TRUE.
                  BTPOP(INDPRP(LABEL)) = .TRUE.
                  LABEL='ZDIPLEN'
                  ATPOP(INDPRP(LABEL)) = .TRUE.
                  BTPOP(INDPRP(LABEL)) = .TRUE.
               GO TO 100
 13            CONTINUE
                  TPALP = .TRUE.
                  LABEL='XDIPLEN'
                  ATPOP( INDPRP(LABEL)) = .TRUE.
                  BTPOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 14            CONTINUE
                  TPALP = .TRUE.
                  LABEL='YDIPLEN'
                  ATPOP( INDPRP(LABEL)) = .TRUE.
                  BTPOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 15            CONTINUE
                  TPALP = .TRUE.
                  LABEL='ZDIPLEN'
                  ATPOP( INDPRP(LABEL)) = .TRUE.
                  BTPOP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 16            CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' PROMPT "',WORD,
     *            '" NOT RECOGNIZED IN TPAINP.'
               CALL QUIT(' ILLEGAL PROMPT IN TPAINP ')
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
      NATPTO = 0
      NBTPTO = 0
      NCTPTO = 0
      IF (ICHANG .GT. 0) THEN
         DO 500 I = 1,NPRLBL
            IF (ATPOP(I)) NATPTO = NATPTO + 1
            IF (BTPOP(I)) NBTPTO = NBTPTO + 1
 500     CONTINUE
         IF (NATPTO .EQ. 0) WRITE (LUPRI,'(/A)')
     *      ' *CUBIC input ignored because no A operators requested.'
         IF (NBTPTO .EQ. 0) WRITE (LUPRI,'(/A)')
     *      ' *CUBIC input ignored because no B operators requested.'
      END IF
      NTPCAL = MIN(NATPTO,NBTPTO,NBTPFR)
      IF  (NTPCAL.GT.0)  THEN
         CALL HEADER('Cubic Response double residue calculation',0)
         WRITE (LUPRI,'(A,L1)')
     *      ' Cubic response two-photon absorption TPAMP='
     *      ,TPAMP
         IF (TPALP) WRITE (LUPRI,'(A,L1)')
     *      ' Excited state polarizability requested TPALP='
     *      ,TPALP
         WRITE(LUPRI,'(/A,I5)')
     *      ' Print level                                    : IPRTPA ='
     *      ,IPRTPA
         WRITE(LUPRI,'(A,I5)')
     *      ' Maximum number of iterations                   : MAXITP ='
     *      ,MAXITP
         WRITE(LUPRI,'(A,1P,D10.3)')
     *      ' Threshold for convergence                      : THCPP  ='
     *      ,THCPP
         WRITE(LUPRI,'(A,I5)')
     *      ' Maximum number of iterations                   : MAXITL ='
     *      ,MAXITL
         WRITE(LUPRI,'(A,1P,D10.3)')
     *      ' Threshold for convergence                      : THCLR  ='
     *      ,THCLR
         WRITE(LUPRI,'(A,I5)')
     *      ' Maximum iterations in optimal orbital algorithm: MAXITO ='
     *      ,MAXITO
         WRITE(LUPRI,'(/I3,A,(1P,5D14.6))')
     *      NBTPFR,' B-frequencies',
     *      (BTPFR(I),I=1,NBTPFR)
      END IF
C
C *** END OF TPAINP
C
      RETURN
      END
