C
C...   Copyright (c) 1997 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dalton, an ab initio electronic structure program, Release 1.0
C...   (1997), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, T. Andersen, K. L. Bak, V. Bakken,
C...   O. Christiansen, P. Dahle, E. K. Dalskov, T. Enevoldsen,
C...   H. Heiberg, D. Jonsson, S. Kirpekar, R. Kobayashi, H. Koch,
C...   K. V. Mikkelsen, P. Norman, M. J. Packer, T.Saue,
C...   P. R. Taylor, and O. Vahtras"
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
C
C$Id: rsptmo1.F,v 1.1.1.1 2001-02-08 13:33:14 hjj Exp $
C
#include <single.h>
      SUBROUTINE TMOINP(WORD)
C
#include <implicit.h>
C
#include <priunit.h>
#include <infrsp.h>
#include <inforb.h>
#include <inftmo.h>
#include <inflr.h>
#include <infpp.h>
#include <rspprp.h>
#include <infpri.h>
#include <infspi.h>
#include <infcr.h>
#include <inflin.h>
#include <inftap.h>
C
      LOGICAL NEWDEF
      PARAMETER ( NTABLE = 19)
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER*8 LABEL
C
      DATA TABLE /'.CPROP ', '.CFREQ ', '.BPROP ', '.BFREQ ',
     *            '.APROP ', '.ROOTS ', '.MAXITP', '.THCPP ', 
     *            '.MAXIT ', '.THCLR ', '.MAXITO', '.PRINT ', 
     *            '.NOHG  ', '.DIPLEN', '.DIPLNX', '.DIPLNY', 
     *            '.DIPLNZ', '.SINGLE', '.FREQUE'/
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
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') GO TO 100
            IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO I=1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,
     &                      15,16,17,18,19), I
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
                     WRITE (LUPRI,'(3(/,A,I5),/)')
     *               ' NUMBER OF FREQUENCIES SPECIFIED    :',NCTMFR,
     *               ' IS GREATER THAN THE ALLOWED NUMBER :',MCTMFR,
     *               ' THE NUMBER IS RESET TO THE MAXIMUM :',MCTMFR
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
                     WRITE (LUPRI,'(3(/,A,I5),/)')
     *               ' NUMBER OF FREQUENCIES SPECIFIED    :',NBTMFR,
     *               ' IS GREATER THAN THE ALLOWED NUMBER :',MBTMFR,
     *               ' THE NUMBER IS RESET TO THE MAXIMUM :',MBTMFR
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
     *               ' NUMBER OF FREQUENCIES SPECIFIED    :',NBTMFR,
     *               ' IS GREATER THAN THE ALLOWED NUMBER :',MBTMFR,
     *               ' THE NUMBER IS RESET TO THE MAXIMUM :',MBTMFR
                     READ (LUCMD,*) (BTMFR(J),J=1,MBTMFR),
     *                              (FFFF,J=MBTMFR+1,NBTMFR)
                     NBTMFR = MBTMFR
                  END IF
                  NCTMFR=NBTMFR
                  DO J=1,NBTMFR
                     CTMFR(J)=BTMFR(J)
                  END DO
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
     *      ' *CRTMO input ignored because no A operators requested.'
         IF (NBTMTO .EQ. 0) WRITE (LUPRI,'(/A)')
     *      ' *CRTMO input ignored because no B operators requested.'
         IF (NCTMTO .EQ. 0) WRITE (LUPRI,'(/A)')
     *      ' *CRTMO input ignored because no C operators requested.'
      END IF
      NTMCAL = MIN(NATMTO,NBTMTO,NCTMTO,NBTMFR,NCTMFR)
      IF  (NTMCAL.GT.0)  THEN
         CALL HEADER('Cubic Response single residue calculation',0)
         WRITE (LUPRI,'(A,L1)')
     *      ' Third moments calculated : TOMOM='
     *      ,TOMOM
         WRITE(LUPRI,'(/A,I5)')
     *      ' Print level                                    : IPRTMO ='
     *      ,IPRTMO
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
     *      NBTMFR,' B-frequencies',
     *      (BTMFR(I),I=1,NBTMFR)
         WRITE(LUPRI,'(I3,A,(1P,5D14.6))')
     *      NCTMFR,' C-frequencies',
     *      (CTMFR(I),I=1,NCTMFR)
      END IF
C
C *** END OF TMOINP
C
      RETURN
      END
