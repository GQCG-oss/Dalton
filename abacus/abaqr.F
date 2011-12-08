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
cs--------------------------------------------------------------cs
cs    abaqr.F=insieme di subroutines per la risposta quadratica cs
cs    ottenuta duplicando ?????                                 cs
cs--------------------------------------------------------------cs
C  /* Deck aqrinp */
      SUBROUTINE AQRINP(WORD)
cs
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (NTABLE = 12)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include <inforb.h>
#include <cbiqr.h>
C
#include <abainf.h>
#include <dorps.h>
#include <nuclei.h>
      DATA TABLE /'.SKIP  ', '.IPRINT','.MAXITE','.THRESH',
     *            '.MAXRED', '.MAXPHP','.XXXXXX',
     *            '.OPTORB', '.PRINT ', '.STOP  ', '.B FREQ',
     *            '.C FREQ' /
C
      NEWDEF = (WORD .EQ. '*QRPROP') 
      ICHANG = 0 
      IF (NEWDEF) THEN
         WORD1 = WORD
cs
cs   va a leggere le successive sottokeywords (.xxx) e verifica se
cs   ci sono quelle elencate in TABLE
cs
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
cs
cs legge e assegna le corrispondenti variabili, logiche e non
cs (quelle numeriche le legge da LUCMD)
cs
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
               WRITE (LUPRI,'(/3A/)') ' Keyword "',WORD,
     *            '" not recognized in AQRINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in AQRINP.')
    1          CONTINUE
                  SKIP = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
    2          CONTINUE
                  READ (LUCMD, '(I5)') IPRINT
                  ICHANG = ICHANG + 1
               GO TO 100
    3          CONTINUE
                  READ (LUCMD,'(I5)') MAXITE
                  ICHANG = ICHANG + 1
               GO TO 100
    4          CONTINUE
                  READ (LUCMD,'(D12.6)') THRESH
                  ICHANG = ICHANG + 1
               GO TO 100
    5          CONTINUE
                  READ (LUCMD,'(I5)') MXRM
                  ICHANG = ICHANG + 1
               GO TO 100
    6          CONTINUE
                  READ (LUCMD,'(I5)') MXPHP
                  ICHANG = ICHANG + 1
               GO TO 100
    7          CONTINUE
               GO TO 100
    8          CONTINUE
                  OOTV   = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
    9          CONTINUE
                  READ (LUCMD, '(I5)') IPRQR
                  ICHANG = ICHANG + 1
               GO TO 100
   10          CONTINUE
                  CUT    = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
   11          CONTINUE
cs
cs possiamo avere piu' frequenze per .B FREQ,  fino ad un massimo
cs pari a NFMAX. LBFREQ=label freq, cioe' freq nro 1, 2,...nfmax.
cs
                  READ(LUCMD, *) LBFREQ 
                  IF(LBFREQ.GT.NFMAX) THEN
                    WRITE(LUPRI,'(A)')
     *               ' Too many frequencies in QR input (B OP.)'
                    CALL QUIT('Too many frequencies in QR input (B)')
                  ENDIF
cs
cs le frequenze sono lette nel vettore QBREQ
cs
                  READ(LUCMD, *) (QBFREQ(J),J=1,LBFREQ)
                  ICHANG = ICHANG + 1
               GO TO 100
   12          CONTINUE
cs idem come sopra
                  READ(LUCMD, *) LCFREQ 
                  IF(LCFREQ.GT.NFMAX) THEN
                    WRITE(LUPRI,'(A)')
     *              ' Too many frequencies in QR input (C OP.)'
                    CALL QUIT('Too many frequencies in QR input (C)')
                  ENDIF
                  READ(LUCMD, *) (QCFREQ(J),J=1,LCFREQ)
                  ICHANG = ICHANG + 1
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/3A/)') ' Prompt "',WORD,
     *            '" not recognized in AQRINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in AQRINP.')
            END IF
      END IF
  300 CONTINUE
C
cs    print session
cs
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for AQRDRV:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' AQRDRV skipped in this run.'
         ELSE
            WRITE (LUPRI,'(A,I5)')
     &         ' General print level in AQRDRV        :',IPRQR
            WRITE (LUPRI,'(A,I5)')
     &         ' Integral print level in AQRDRV        :',IPRINT
            WRITE (LUPRI,'(A,1P,E9.2)')
     &           ' Threshold in AQRDRV          :',THRESH
            WRITE(LUPRI,'(A,I5)')' Maximum iterations in AQRDRV :',
     &                           MAXITE
            IF (CUT) THEN
               WRITE (LUPRI,'(/A)') ' Program is stopped after AQRDRV.'
            END IF
cs ricordarsi modificare
            IF (LBFREQ.GT.0) THEN
               DO IBFR=1,LBFREQ
               WRITE (LUPRI,1968) IBFR, QBFREQ(IBFR)
               END DO
            END IF
            IF (LCFREQ.GT.0) THEN
               DO ICFR=1,LCFREQ
               WRITE (LUPRI,1969) ICFR, QCFREQ(ICFR)
               END DO
            END IF
cs ricordarsi modificare
         END IF
      END IF
      RETURN
1968  FORMAT (/1X,'bfreq # =',I4,4X,'bfreq=',F10.5)
1969  FORMAT (/1X,'cfreq # =',I4,4X,'cfreq=',F10.5)
      END
c
c    ---------------- end of subroutine AQRINP ---------------------
c
C  /* Deck aqrini */
      SUBROUTINE AQRINI
C
C     Initialize /QRSRES/
C
#include <implicit.h>
#include <mxcent.h>
#include <cbiqr.h>
#include <abainf.h>
C
cs    iprdef= deve essere il default (applicato per entrambi i print level)
cs
      IPRINT = IPRDEF
      IPRQR  = IPRINT
      SKIP   = .FALSE.
      CUT    = .FALSE.
      OOTV   = .FALSE.
      THRESH = 1.D-04
      MAXITE = 60
      MXRM   = 400
      MXPHP  = 0
cs    LABAPP = 0 ?????????
      NABAPP = 0
      LBFREQ = 1
      LCFREQ = 1
      CALL DZERO (QBFREQ,NFMAX)
      CALL DZERO (QCFREQ,NFMAX)
C
      RETURN
      END
c    ---------------------------------------
C  /* Deck aqrdrv */
      SUBROUTINE AQRDRV(WORK,LWORK)
cs
cs    driver vero e proprio della risposta quadratica in abacus
cs    chiamata quando .MCD, .HYPER e .VERDET = .TRUE.
cs
#include <implicit.h>
#include <iratdef.h>
#include <dummy.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <abainf.h>
#include <exeinf.h>
#include <inforb.h>
#include <cbiqr.h>
#include <infrsp.h>
#include "rspprp.h"
#include <infhyp.h>
#include <infsmo.h>
#include <infdim.h>
#include <inftap.h>
#include <nuclei.h>
#include <symmet.h>
      LOGICAL PASS, HFCLC, CICLC, TRIPLE, EXECLC, FOUND
      CHARACTER*8 LABEL, DLAB(9*MXCENT)
      DIMENSION WORK(LWORK)
C
      IF (SKIP) RETURN
      CALL QENTER('AQRDRV')
      CALL TIMER('START ',TIMEIN,TIMOUT)
      IF (IPRINT .GT. 0) WRITE (LUPRI,'(A,/)')
     *    '  <<<<<<<<<< Output from AQRDRV >>>>>>>>>> '
C
      PASS   = .TRUE.
      IPRRSP = IPRQR
C
C     Get reference state
C     ===================
C
C     1. Work allocations:
C
      KCMO   = 1
      KUDV   = KCMO + NCMOT
      KXINDX = KUDV + N2ASHX
      KWORK1 = KXINDX + LCINDX
      LWORK1 = LWORK - KWORK1
C
      CALL RD_SIRIFC('CMO',FOUND,WORK(KCMO),WORK(KWORK1),LWORK1)
      IF (.NOT.FOUND) CALL QUIT('AQRDRV error: CMO not found on SIRIFC')
      IF (NASHT .GT. 0) THEN
         CALL RD_SIRIFC('DV',FOUND,WORK(KWORK1),WORK(KWORK1),LWORK1)
         IF (.NOT.FOUND)
     &      CALL QUIT('AQRDRV error: DV not found on SIRIFC')
         CALL DSPTSI(NASHT,WORK(KWORK1),WORK(KUDV))
      END IF
C
      LFREE = LWORK1
C
      CALL GETCIX(WORK(KXINDX),IREFSY,IREFSY,WORK(KWORK1),LFREE,0)
C
C     Construct property-integrals and write to LUPROP
C     ================================================
C
C     2. Work allocations:
C
      KIDSYM = KWORK1
      KIDADR = KIDSYM + 9*MXCENT
      KWORK2 = KIDADR + 9*MXCENT
      LWORK2 = LWORK  - KWORK2
      IF (LWORK2.LT.0) CALL QUIT ('No Free Memory left in AQRDRV')
C
      NLBTOT = 0
C
cs
      CALL RSPSET
      IF (HYPER .OR. VERDET .OR. MCD) THEN
         NCOMP  = 0
         NPATOM = 0
         CALL GET1IN(DUMMY,'DIPLEN ',NCOMP,WORK(KWORK2),LWORK2,
     &               DLAB,WORK(KIDSYM),WORK(KIDADR),
     &               IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,DUMMY,
     &               IPRINT)
         NLAB = 3
         CALL LABCOP(NLAB,NLBTOT,DLAB,WORK(KIDSYM),LABAPP,LABSYM)
C
cs
         IF (VERDET .OR. MCD) THEN
            IF (NOLOND) THEN
               CALL GET1IN(DUMMY,'ANGMOM ',NCOMP,WORK(KWORK2),LWORK2,
     &                     DLAB,WORK(KIDSYM),WORK(KIDADR),
     &                     IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,
     &                     DUMMY,IPRINT)
               NLAB = 3
               CALL LABCOP(NLAB,NLBTOT,DLAB,WORK(KIDSYM),
     &                     LABAPP,LABSYM)
            ELSE
              CALL LABCOP(1,NLBTOT,'XLONMAG ',ISYMAX(1,2),LABAPP,LABSYM)
              CALL LABCOP(1,NLBTOT,'YLONMAG ',ISYMAX(2,2),LABAPP,LABSYM)
              CALL LABCOP(1,NLBTOT,'ZLONMAG ',ISYMAX(3,2),LABAPP,LABSYM)
            END IF
C
            IF (NODIFC) THEN
               NCOMP  = 0
               NPATOM = 0
               CALL GET1IN(DUMMY,'HBDO   ',NCOMP,WORK(KWORK2),LWORK2,
     &                     DLAB,WORK(KIDSYM),WORK(KIDADR),
     &                     IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,
     &                     DUMMY,IPRINT)
               NLAB = 3
               CALL LABCOP(NLAB,NLBTOT,DLAB,WORK(KIDSYM),LABAPP,
     &                     LABSYM)
            ENDIF
         END IF
      ENDIF
      DO 10 I=1,NLBTOT
         IF (HYPER .OR. VERDET) THEN 
            HYPCAL = .TRUE.
            IF (LABAPP(I)(2:7) .EQ. 'LONMAG' .OR. 
     &          LABAPP(I)(2:7) .EQ. 'ANGMOM') THEN
               CQROP( INDPRP(LABAPP(I))) = .TRUE. 
            ELSE
               IF (HYPER) CQROP(INDPRP(LABAPP(I))) = .TRUE.
               AQROP( INDPRP(LABAPP(I))) = .TRUE. 
               BQROP( INDPRP(LABAPP(I))) = .TRUE. 
            END IF
            NCQRFR = LCFREQ
            NBQRFR = LBFREQ
            CALL DCOPY (LBFREQ,QBFREQ,1,BQRFR,1)
            CALL DCOPY (LCFREQ,QCFREQ,1,CQRFR,1)
         ENDIF
         IF (MCD) THEN 
            IF (LABAPP(I)(2:7) .EQ. 'LONMAG' .OR. 
     &          LABAPP(I)(2:7) .EQ. 'ANGMOM') THEN
               ASMOP( INDPRP(LABAPP(I))) = .TRUE.
            ELSE
               BSMOP( INDPRP(LABAPP(I))) = .TRUE.
            END IF
            NBSMFR = LBFREQ
cs          NBSMFR = NBFREQ
cs          CALL DCOPY (NBFREQ,BFREQ,1,BSMFR,1)
           CALL DCOPY (LBFREQ,QBFREQ,1,BSMFR,1)
         ENDIF
   10 CONTINUE
C
C     We might need a somewhat more intelligent test for FABRHS
C
      FABRHS = .TRUE.
      CICLC  = .FALSE.
      HFCLC  = NASHT .LE. 1
      TRIPLE = .FALSE.
      EXECLC = .TRUE.
      NABATY = 1
      NABAOP = 1
cs
      NEXENG = 0
C
C     The following two lines are needed in order to pass
C     consistency checks in RSPMC
C
      CALL ABAVAR(1,.FALSE.,0,WORK(KWORK2),LWORK2)
      IF (HFCLC) NCONF = 1
C
cs    Probabilmente con questa chiamata parte il calcolo vero e proprio
cs    ABARSP e' in abarspn.F
cs
      ISYM = 1
      LUSOVE = -1
      LUGDVE = -1
      LUREVE = -1
      CALL GPOPEN(LUSOVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUGDVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUREVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL QUIT('ABARSP call in ABAQR needs to be updated')
C     hjaaj mar 2004: inserted 'LABEL', then number of arguments fits.
C      LABEL needs to be defined.
      CALL ABARSP(CICLC,HFCLC,TRIPLE,OOTV,ISYM,EXECLC,
     *            DUMMY,NEXENG,NABATY,NABAOP,LABEL,LUGDVE,LUSOVE,LUREVE,
     *            THRESH,MAXITE,IPRQR,MXRM,MXPHP,WORK(KWORK2),
     *            LWORK2)
C
C     The following line is needed to reset variables
C
C      CALL ABAVAR(ISYM,TRIPLE,0,WORK(KWORK2),LWORK2)
C
C
      CALL GPCLOSE(LUSOVE,'DELETE')
      CALL GPCLOSE(LUGDVE,'DELETE')
      CALL GPCLOSE(LUREVE,'DELETE')
      CALL TIMER ('AQRDRV',TIMEIN,TIMOUT)
      IF (CUT) THEN
         WRITE (LUPRI,'(/,A)')
     &          ' Program stopped after AQRDRV as required.'
         WRITE (LUPRI,'(A)') ' No restart file has been written.'
         CALL QUIT(' ***** End of DALTON (in AQRDRV) *****')
      END IF
      CALL QEXIT('AQRDRV')
      RETURN
      END
