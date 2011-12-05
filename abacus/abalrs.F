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
C  /* Deck lrsinp */
      SUBROUTINE LRSINP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (NTABLE = 15)
      LOGICAL SET, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include "cbilrs.h"
C
#include <gnrinf.h>
#include <abainf.h>
#include <anrinf.h>
#include <dorps.h>
#include <nuclei.h>
      SAVE SET
      DATA TABLE /'.SKIP  ', '.PRINT ','.MAX IT','.THRESH',
     *            '.MAXRED', '.MAXPHP','.XXXXXX','.XXXXXX',
     *            '.OPTORB', '.XXXXXX','.XXXXXX','.XXXXXX',
     *            '.XXXXXX', '.XXXXXX','.STOP  '/
      DATA SET/.FALSE./
C
      IF (SET) THEN
         IF (WORD .NE. '*END OF') THEN
 969        READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .NE. '*') GO TO 969
         END IF
         RETURN
      END IF
C
      SET = .TRUE.
      CALL LRSINI
C
      NEWDEF = (WORD .EQ. '*LINRES')
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
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in LRSINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in LRSINP.')
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD,*) IPRCLC
               GO TO 100
    3          CONTINUE
                  READ (LUCMD,*) MAXCLC
               GO TO 100
C              .THRESH
    4          CONTINUE
                  READ (LUCMD,*) THRCLC
               GO TO 100
    5          CONTINUE
                  READ (LUCMD,*) MXRM
               GO TO 100
    6          CONTINUE
                  READ (LUCMD,*) MXPHP
               GO TO 100
    7             CONTINUE
               GO TO 100
    8          CONTINUE
                  CONTINUE
               GO TO 100
    9             OOTV   = .TRUE.
               GO TO 100
   10             CONTINUE
               GO TO 100
   11              CONTINUE
               GO TO 100
   12             CONTINUE
               GO TO 100
   13             CONTINUE
               GO TO 100
   14             CONTINUE
               GO TO 100
   15             CUT    = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in LRSINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in LRSINP.')
            END IF
      END IF
  300 CONTINUE
      IF (THR_REDFAC .GT. 0.0D0) THEN
         ICHANG = ICHANG + 1
         WRITE (LUPRI,'(3A,1P,D10.2)') '@ INFO ',WORD1,
     &   ' thresholds multiplied with general factor',THR_REDFAC
         THRCLC = THRCLC*THR_REDFAC
      END IF
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for LINRES:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' LINRES skipped in this run.'
         ELSE
            WRITE (LUPRI,'(A)')
     &         ' Singlet linear response module for properties',
     &         ' as VCD, MAGSUS, SPIN-SPIN, SHIELD, SPINRO, MOLGFA'
            WRITE (LUPRI,'(A,I5)')
     *         ' Print level in LINRES        :',IPRCLC
            WRITE (LUPRI,'(A,1P,D9.2)')
     *         ' Threshold in LINRES          :',THRCLC
            WRITE(LUPRI,'(A,I5)')
     &         ' Maximum iterations in LINRES :',MAXCLC
            IF (CUT) THEN
               WRITE (LUPRI,'(/,A)') ' Program is stopped after LINRES.'
            END IF
         END IF
      END IF
      IPRNR = IPRCLC
      RETURN
      END
C  /* Deck lrsini */
      SUBROUTINE LRSINI
C
C     Initialize /LRSRES/
C
#include <implicit.h>
#include <mxcent.h>
#include <abainf.h>
#include "cbilrs.h"
C
      IPRCLC = IPRDEF
      SKIP   = .FALSE.
      CUT    = .FALSE.
      OOTV   = .FALSE.
      IF (SPINRO .OR. SHIELD .OR. (SHIELD .AND. MAGSUS)
     &          .AND. .NOT. SPNSPN) THEN
         THRCLC = 1.D-04
      ELSE
         THRCLC = 2.D-03
      END IF
      MAXCLC = 60
      MXRM   = 400
      MXPHP  = 0
      RETURN
      END
C  /* Deck lrsdrv */
      SUBROUTINE LRSDRV(WORK,LWORK,PASS)
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <iratdef.h>
#include <priunit.h>
      LOGICAL PASS
      LOGICAL INTTRS, CICLC, HFCLC, TRPCLC
      DIMENSION WORK(LWORK)
C
#include "cbilrs.h"
#include <abainf.h>
#include <gnrinf.h>
#include <spnout.h>
#include <inflin.h>
#include <infvar.h>
#include <infdim.h>
#include <inforb.h>
#include <nuclei.h>
C
      IF (SKIP.OR..NOT.
     &   (MAGSUS .OR.  VCD .OR. SHIELD .OR. SPINRO .OR. MOLGFA
     &           .OR. (SPNSPN.AND.DOPSO))) RETURN
      CALL QENTER('LRSDRV')
      IF (IPRCLC .GT. 0) THEN
         CALL TIMER('START ',TIMEIN,TIMOUT)
         WRITE (LUPRI,'(//A,/)')
     *    '  <<<<<<<<<< Output from LRSDRV >>>>>>>>>> '
      END IF
C
      PASS = .TRUE.
C
C     Allocations
C
      KNABAT = 1
      KZERO  = KNABAT + (3*NUCDEP + 6 + 1)/IRAT
      KLAST  = KZERO  + (3*NUCDEP + 6 + 1)/IRAT
      IF (KLAST .GT. LWORK) CALL STOPIT('LRSDRV',' ',KLAST,LWORK)
      KWRK   = KLAST
      LWRK   = LWORK - KLAST + 1
      CALL LRSDR1(WORK(KNABAT),WORK(KZERO),WORK(KWRK),LWRK)
      IF (IPRCLC .GT. 0) CALL TIMER ('LRSDRV',TIMEIN,TIMOUT)
      PASS = .TRUE.
      IF (CUT) THEN
         WRITE (LUPRI,'(/A)')
     &      ' Program stopped after LRSDRV as requested.'
         CALL QUIT(' ***** End of ABACUS (in LRSDRV) *****')
      END IF
      CALL QEXIT('LRSDRV')
      RETURN
      END
C  /* Deck lrsdr1 */
      SUBROUTINE LRSDR1(NABATY,ZERNRM,WORK,LWORK)
#include <implicit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <iratdef.h>
#include <priunit.h>
      LOGICAL PASS
      LOGICAL INTTRS, CICLC, HFCLC, TRPCLC, ZERNRM(*)
      DIMENSION WORK(LWORK), NABATY(*)
C
#include "cbilrs.h"
#include <abainf.h>
#include <inflin.h>
#include <infvar.h>
#include <infdim.h>
#include <inforb.h>
#include <nuclei.h>
C
      HFCLC  = NASHT .LE. 1
      IF (HFCLC) NCONF = 1
      LUREVE = -1
      LUSOVE = -1
      LUGDVE = -1
      CALL GPOPEN(LUSOVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUGDVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUREVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      IF (VCD .OR. SHIELD .OR. SPNSPN .OR. MAGSUS .OR. SPINRO .OR.
     &    MOLGFA) THEN
         CALL GETMRH(HFCLC,OOTV,NABATY,NABAOP,LUGDVE,LUSOVE,LUREVE,
     &               THRCLC,MAXCLC,IPRCLC,MXRM,MXPHP,ZERNRM,WORK,LWORK)
      ELSE
#ifndef use_getrhs
         CALL QUIT('program error: LRSDR1 not defined for this call')
#else
c is the following old code for NACME which we might want to reactivate?
c -- hjaaj Aug 2004
         TRPCLC = .FALSE.
         INTTRS = .TRUE.
         CICLC  = .FALSE.
         IOPSYM = 1
         NABAOP = 1
         CALL GETRHS(IOPSYM,NEXVAL,LUGDVE,IPRCLC,WORK,LWORK)
         DO 100 I = 1, NABAOP
            NABATY(I) = 1
 100     CONTINUE
         CALL ABARSP(CICLC,HFCLC,TRPCLC,OOTV,IOPSYM,EXCLC,
     &               EXVAL,NEXVAL,NABATY,NABAOP,'DUMMY   ',LUGDVE,
     &               LUSOVE,LUREVE,THRCLC,MAXCLC,IPRCLC,MXRM,MXPHP,
     &               WORK,LWORK)
#endif /* use_getrhs */
C
      END IF
      CALL GPCLOSE(LUSOVE,'DELETE')
      CALL GPCLOSE(LUGDVE,'DELETE')
      CALL GPCLOSE(LUREVE,'DELETE')
      RETURN
      END
C  /* Deck getmrh */
      SUBROUTINE GETMRH(HFCLC,OOTV,NABATY,NABAOP,LUGDVE,LUSOVE,LUREVE,
     &                  THRCLC,MAXCLC,IPRCLC,MXRM,MXPHP,ZERNRM,WORK,
     &                  LWORK)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (D10 = 10.D0, DP5 = 0.5D0)
C
      LOGICAL EXCLC, TRPCLC, CICLC, HFCLC, OOTV, ZERNRM(*), OLDDX
      DIMENSION NABATY(*), WORK(LWORK)
#include <inftap.h>
#include <infvar.h>
#include <infdim.h>
#include <inflin.h>
#include <abainf.h>
#include <spnout.h>
#include <nuclei.h>
#include <inforb.h>
#include <gdvec.h>
#include <abares.h>
#include <chrnos.h>
C
C     Open direct access files
C
      CALL GPOPEN(LUGDI,ABAGDI,'UNKNOWN','DIRECT',' ',IRAT*NVARMA,OLDDX)
      CALL GPOPEN(LURDI,ABARDI,'UNKNOWN','DIRECT',' ',IRAT*NVARMA,OLDDX)
C
      EXCLC  = .FALSE.
      TRPCLC = .FALSE.
      CICLC  = .FALSE.
      EXVAL  = 0.0D0
      NEXVAL = 1
C
      DO 50 ISYM = 1, NSYM
         NREC = NGDVEC(ISYM,2)
         IF (NREC .GT. 0) THEN
            CALL ABAVAR(ISYM,TRPCLC,IPRCLC,WORK,LWORK)
         IF (NVARPT .EQ. 0) GO TO 50
            IF (HFCLC) NCONF = 1
C
C           Collect RHS's
C
            IF (IPRCLC .GT. 10) THEN
               CALL HEADER('RHS vectors in GETMRH for symmetry '//
     &                      CHRNOS(ISYM)//':',-1)
            END IF
            DO 10 I = 1, NREC
               NABAOP = 0
               REWIND LUGDVE
               IREC   = IGDREC(I,ISYM,2)
               ICOOR  = IGDCOR(I,ISYM,2)
               IF (IDORCI(IREC,2) .GT. 1) GO TO 10
               IF (IREC .GT. (3*NUCDEP + 3)) GOTO 10
               IF ((SPNSPN .AND. DOPSO) .OR. (ICOOR .LT. 0)) THEN
                  CALL READDX (LUGDI,IREC,IRAT*NVARPT,WORK)
                  DNORM = DNRM2(NVARPT,WORK,1)
                  IF (DNORM .GT. THRCLC/D10) THEN
                     NABAOP = NABAOP + 1
                     ZERNRM(I) = .FALSE.
                     CALL WRITT(LUGDVE,NVARPT,WORK)
                  ELSE
                     ZERNRM(I) = .TRUE.
                  END IF
                  IF (IPRCLC .GT. 10) THEN
                     WRITE (LUPRI,'(A,I5)') ' Coordinate:',ICOOR
                     WRITE (LUPRI,'(A,I5)') ' Record:    ',IREC
                     CALL OUTPUT(WORK,1,NVARPT,1,1,NVARPT,1,1,LUPRI)
                  END IF
               END IF
C
C           Solve equations
C
               IF (NABAOP .EQ. 1) THEN
                  NABATY(1) = -1
                  IF (.NOT. ZERNRM(I)) THEN
                     NOCONV = .FALSE.
                     FINRES = 0.0D0
                     CALL ABARSP(CICLC,HFCLC,TRPCLC,OOTV,ISYM,EXCLC,
     &                           EXVAL,NEXVAL,NABATY,1,IMGLAB(IREC),
     &                           LUGDVE,LUSOVE,LUREVE,THRCLC,MAXCLC,
     &                           IPRCLC,MXRM,MXPHP,WORK,LWORK)
                     IF (NOCONV) WRITE
     &                 (LUPRI,'(/,3A,2F10.5)') ' No convergence for ',
     &                    IMGLAB(IREC),'. Last residual and threshold:',
     &                    FINRES, THRCLC
                  END IF
               END IF
C
C           Write solutions and residuals
C
               REWIND LUSOVE
               REWIND LUREVE
               IF ((SPNSPN .AND. DOPSO) .OR. (ICOOR .LT. 0)) THEN
C     
C              Solution vector
C
                  IF (ZERNRM(I)) THEN
                     CALL DZERO(WORK,2*NVARPT)
                  ELSE
                     CALL READT(LUSOVE,NVARPT,WORK)
                     CALL READT(LUREVE,NVARPT,WORK(NVARPT + 1))
C
C                 Divide solution by 2 in accordance with ABACUS
C                 solver
C
                     CALL DSCAL(NVARPT,DP5,WORK,1)
                  END IF
                  IREC2 = 2*IREC - 1
                  CALL WRITDX (LURDI,IREC2,IRAT*NVARPT,WORK)
                  IREC2 = 2*IREC
                  CALL WRITDX (LURDI,IREC2,IRAT*NVARPT,
     &                         WORK(NVARPT + 1))
ckr                  IDORCI(IREC,2) = 2
                  CALL ABAWRITRS
                  IF (IPRCLC .GT. 5) THEN
                     CALL HEADER('Response and residuum vectors in '//
     &                   'GETMRH for symmetry '//CHRNOS(ISYM)//':',-1)
                     WRITE (LUPRI,'(A,I5)') ' Coordinate:',ICOOR
                     WRITE (LUPRI,'(A,I5)') ' Record:    ',IREC
                     DN1 = DNRM2(NVARPT,WORK,1)
                     DN2 = DNRM2(NVARPT,WORK(1+NVARPT),1)
                     WRITE (LUPRI,'(A,2(D20.14,2X))') 
     &                      ' Norms:     ',DN1,DN2
                  END IF
                  IF (IPRCLC .GT. 10) THEN
                     CALL OUTPUT(WORK,1,NVARPT,1,2,NVARPT,1,2,LUPRI)
                  END IF
               END IF
 10         CONTINUE
         END IF
 50   CONTINUE
      CALL GPCLOSE(LUGDI,'KEEP')
      CALL GPCLOSE(LURDI,'KEEP')
      RETURN
      END
#ifdef use_getrhs
C  /* Deck getrhs */
      SUBROUTINE GETRHS(IOPSYM,NEXVAL,LUGDVE,IPRINT,WORK,LWORK)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxorb.h>
#include <mxcent.h>
C
      LOGICAL TRIPLE, OLDDX
      DIMENSION WORK(LWORK)
#include <inftap.h>
#include <infvar.h>
#include <infdim.h>
#include <inflin.h>
#include <gdvec.h>
#include <abainf.h>
#include <nuclei.h>
C
      TRIPLE = .FALSE.
      CALL ABAVAR(IOPSYM,TRIPLE,IPRINT,WORK,LWORK)
      CALL GPOPEN(LUGDR,ABAGDR,'UNKNOWN','DIRECT',' ',IRAT*NVARMA,OLDDX)
      REWIND LUGDVE
      IF (NVARPT .GT. LWORK) CALL STOPIT('GETRHS',' ',NVARPT,LWORK)
      DO 100 IOP = 1, NEXVAL
         IREC = IGDREC(IOP,IOPSYM,1)
         CALL READDX(LUGDR,IREC,IRAT*NVARPT,WORK)
         CALL WRITT(LUGDVE,NVARPT,WORK)
 100  CONTINUE
      CALL GPCLOSE(LUGDR,'KEEP')
      RETURN
      END
#endif /* uyy1Gse_getrhs */
C--- end of abalrs.F ---
