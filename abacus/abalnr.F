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
C  /* Deck lnrinp */
      SUBROUTINE LNRINP(WORD)
C
C     Nov. and Dec. 93
C     Written by K.L.Bak and P.Joergensen using EXCITA as template.
C     Purpose: To enable calculations of frequency dependent second
C     order response properties in ABACUS. In particular polariza-
C     bilities and vibrational raman optical activity (VROA).
C
C     Oct - Dev. 99: Extended by K.L.Bak for calculation of
C     vibrational g-factors.
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <cbiexc.h>
#include <gnrinf.h>
#include <cbilnr.h>
#include <abainf.h>
#include <anrinf.h>
#include <dorps.h>
#include <nuclei.h>
#include <absorp.h>
#include <maxaqn.h>
#include <symmet.h>
#include <codata.h>
#include <pcmlog.h>
C
      PARAMETER (NTABLE = 16)
      LOGICAL NEWDEF, LOCFIN
      CHARACTER*8 DIPLEN(3),DIPLOC(3), ANGMOM(3), LONMAG(3)
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
      DATA DIPLEN/'XDIPLEN','YDIPLEN','ZDIPLEN'/
      DATA DIPLOC/'XDIPLOC','YDIPLOC','ZDIPLOC'/
      DATA ANGMOM/'XANGMOM','YANGMOM','ZANGMOM'/
      DATA LONMAG/'XLONMAG','YLONMAG','ZLONMAG'/
C
      DATA TABLE /'.SKIP  ', '.WAVELE','.MAX IT','.THRESH',
     &            '.MAXRED', '.MAXPHP','.STATIC','.VIBGRE',
     *            '.OPTORB', '.FREQUE','.BATCH ','.PRINT ',
     *            '.DAMPIN', '.STOP  ','.FREQ I','.REDUCE'/
C
      NEWDEF = (WORD .EQ. '*ABALNR')
      ICHANG = 0
      NFRVLT = 0
      IF (NEWDEF) THEN
         LOCFIN = .FALSE.
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
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/3A/)') ' Keyword "',WORD,
     *            '" not recognized in LNRINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword under *ABALNR.')
    1          CONTINUE
                  SKIP = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
    2          CONTINUE
                  READ(LUCMD,*) NWAVEL
                  IF (NWAVEL .LT. 0) THEN
                     CALL QUIT('# wavelengths negative under *ABALNR')
                  ELSE IF (NWAVEL .EQ. 0) THEN
                     GO TO 100
                  ELSE
                     IF (.NOT. LOCFIN) THEN
                        NFRVAL = 0
                        LOCFIN = .TRUE.
                     END IF
                  END IF
                  NFRVLT = NFRVLT + NWAVEL
                  NFTOT  = NWAVEL + NFRVAL
                  IF (NFTOT .GT. MXFR) THEN
                     NWAVEL = MXFR - NFRVAL
                  END IF
                  READ(LUCMD,*) (FRVAL(NFRVAL+I),I=1,NWAVEL)
                  DO I = NFRVAL+1,NFRVAL+NWAVEL
                     IF (ABS(FRVAL(I)) .LT. 1.0D-12) THEN
                        WRITE(LUPRI,*) 'Wavelength too small:',
     &                                 FRVAL(I)
                        WRITE(LUPRI,*) 'Use frequency input instead!'
                        CALL QUIT('Wavelength too small under *ABALNR')
                     END IF
                     FRVAL(I) = XTNM/FRVAL(I)
                  END DO
                  NFRVAL = NFRVAL + NWAVEL
                  ICHANG = ICHANG + 1
               GO TO 100
    3          CONTINUE
                  READ (LUCMD,*) INPVAL
                  ICHANG = ICHANG + 1
                  IF (INPVAL .EQ. MAXITE) THEN
                     ICHANG = ICHANG - 1
                  ELSE
                     MAXITE = INPVAL
                  END IF
               GO TO 100
C              .THRESH
    4          CONTINUE
                  READ (LUCMD,*) DTHCLN
                  ICHANG = ICHANG + 1
                  IF (DTHCLN .EQ. THCLNR) THEN
                     ICHANG = ICHANG - 1
                  ELSE
                     THCLNR = DTHCLN
                  END IF
               GO TO 100
    5          CONTINUE
                  READ (LUCMD,*) IMXRM
                  ICHANG = ICHANG + 1
                  IF (IMXRM .EQ. MXRM) THEN
                     ICHANG = ICHANG - 1
                  ELSE
                     MXRM = IMXRM
                  END IF
               GO TO 100
    6          CONTINUE
                  READ (LUCMD,*) IMXPHP
                  ICHANG = ICHANG + 1
                  IF (IMXPHP .EQ. MXPHP) THEN
                     ICHANG = ICHANG - 1
                  ELSE
                     MXPHP = IMXPHP
                  END IF
               GO TO 100
    7          CONTINUE
                  STATIC = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
    8          CONTINUE ! .VIBGREP
CSPAS:20/3-2011: allowing for perturbations in only one irrep
                  VIBGIR = .TRUE.
                  READ(LUCMD,*) IRVIBG
CKeinSPASmehr
               GO TO 100
    9             OOTV   = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
   10             CONTINUE
                  READ(LUCMD,*) NRDFR
                  IF (NRDFR .LT. 0) THEN
                     CALL QUIT('# frequencies negative under *ABALNR')
                  ELSE IF (NRDFR .EQ. 0) THEN
                     GO TO 100
                  ELSE
                     IF (.NOT. LOCFIN) THEN
                        NFRVAL = 0
                        LOCFIN = .TRUE.
                     END IF
                  END IF
                  NFRVLT = NFRVLT + NRDFR
                  NFTOT  = NRDFR  + NFRVAL
                  IF (NFTOT .GT. MXFR) THEN
                     NRDFR = MXFR - NFRVAL
                  END IF
                  READ(LUCMD,*) (FRVAL(NFRVAL+I),I=1,NRDFR)
                  NFRVAL = NFRVAL + NRDFR
                  ICHANG = ICHANG + 1
               GO TO 100
   11             CONTINUE
                  READ(LUCMD,*) NFREQ_BATCH
                  IF (NFREQ_BATCH.GT.MXFR) THEN
                     NFREQ_BATCH = MXFR
                  END IF
               GO TO 100
   12             CONTINUE
                  READ (LUCMD,*) IPRLNR
                  ICHANG = ICHANG + 1
                  IF (IPRLNR .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
   13             CONTINUE
                  ABSORP = .TRUE.
                  ICHANG = ICHANG + 1
                  READ (LUCMD,*) DAMPING
               GO TO 100
   14             CONTINUE
                  CUT    = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
 15            CONTINUE
                  ABS_INTERVAL = .TRUE.
                  ABS_REDUCE = .TRUE.
                  READ(LUCMD,*) (FREQ_INTERVAL(I), I=1,3)
               GO TO 100
 16            CONTINUE
                 ABS_REDUCE=.TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in LNRINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt under *ABALNR.')
            END IF
      END IF
  300 CONTINUE
C
      IF (THR_REDFAC .GT. 0.0D0) THEN
         ICHANG = ICHANG + 1
         WRITE (LUPRI,'(3A,1P,D10.2)') '@ INFO ',WORD1,
     &   ' thresholds multiplied with general factor',THR_REDFAC
         DTHCLN = DTHCLN*THR_REDFAC
      END IF
      CALL IZERO(NOPER,8)
      DO ILABEL=1,3
         ISYM = ISYMAX(ILABEL,1) + 1
         NOPER(ISYM) = NOPER(ISYM) + 1
         LABOP(NOPER(ISYM),ISYM) = DIPLEN(ILABEL)
         IF(PCM.AND.LOCFLD) THEN
            LABOP(NOPER(ISYM),ISYM) = DIPLOC(ILABEL)
         END IF
      END DO
      IF (STATIC) THEN
         DO ILABEL=1,3
            ISYM = ISYMAX(ILABEL,2) + 1
            NOPER(ISYM) = NOPER(ISYM) + 1
            LABOP(NOPER(ISYM),ISYM) = ANGMOM(ILABEL)
         END DO
         DO ILABEL=1,3
            ISYM = ISYMAX(ILABEL,2) + 1
            NOPER(ISYM) = NOPER(ISYM) + 1
            LABOP(NOPER(ISYM),ISYM) = LONMAG(ILABEL)
         END DO
      END IF
C
      IF ((VROA .OR. RAMAN .OR. OPTROT) .AND. NODIFC) THEN
         NWARN = NWARN + 1
         WRITE (LUPRI,'(/2A/A/A/)')
     &   ' WARNING: Raman properties calculation is NOT ',
     &   ' implemented with the NODIFC keyword active. ',
     &   ' WARNING: Raman properties calculation ignord',
     &   ' WARNING: Try calculation again without specifying NODIFC.'
         ROAA = .FALSE.
         ROAG = .FALSE.
      ELSEIF (VROA) THEN
         ROAA = .TRUE.
         ROAG = .TRUE.
      ELSE IF (OPTROT) THEN
         ROAA = .FALSE.
         ROAG = .TRUE.
      ELSE IF (RAMAN) THEN
         ROAA = .FALSE.
         ROAG = .FALSE.
         ALFA = .TRUE.
      ELSE
         ROAA = .FALSE.
         ROAG = .FALSE.
      END IF
C
      IF (VIB_G .AND. NODIFC) THEN
         WRITE(LUPRI,'(/,3A,/)')
     &   ' Vibrational g-factors are NOT',
     &   ' implemented with the NODIFC keyword active. ',
     &   ' Try calculation again without specifying NODIFC.'
         VIB_G = .FALSE.
      END IF
C
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for *ABALNR:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' LNRABA skipped in this run.'
         ELSE
            IF (ABSORP) THEN
               IF (ABS_INTERVAL) THEN
                  NFREQ_INTERVAL = 1 +
     &                 INT((FREQ_INTERVAL(2)-FREQ_INTERVAL(1))/
     &                 FREQ_INTERVAL(3))
                  NFRVAL = NFREQ_INTERVAL
                  DO I=1,NFRVAL
                     FRVAL(I)=FREQ_INTERVAL(1) + (I-1)*FREQ_INTERVAL(3)
                  END DO                  
               ELSE
                  NFREQ_ALPHA = NFRVAL
                  NFREQ_INTERVAL = NFRVAL
                  DO I=1,NFREQ_ALPHA
                     FREQ_ALPHA(I) = FRVAL(I)
                  END DO
               END IF
               THCLR_ABSORP = THCLNR
               IPRABS = IPRLNR
               MAX_MACRO  = MAXITE
               IF(ABS_INTERVAL .AND.
     &              (FREQ_INTERVAL(1).GT.FREQ_INTERVAL(2) .OR.
     &           (FREQ_INTERVAL(2)-FREQ_INTERVAL(1)).LT.FREQ_INTERVAL(3)
     &            .OR. FREQ_INTERVAL(3).LE.0.0D0) ) THEN
                  WRITE(LUPRI,'(/,A,/,A,/,/,A,/,A,/,/,A)')
     & ' Warning: Frequency interval for absorption calculation is not',
     & '          correctly specified according to:',
     & '                 .FREQ_INTERVAL',
     & '                 STARTFREQ ENDFREQ STEPSIZE',
     & '          This input will be ignored. Calculation continues.'
                  ABS_INTERVAL = .FALSE.
                  FREQ_INTERVAL(1)=0.0D0
                  FREQ_INTERVAL(2)=0.0D0
                  FREQ_INTERVAL(3)=1.0D0
               END IF
            ELSE IF (ABS_INTERVAL) THEN
               ABS_INTERVAL = .FALSE.
               WRITE(LUPRI,'(/A,/A)')
     & ' Warning: Frequency interval input only available for'//
     & '          absorption calculation.',
     & '          This input will be ignored. Calculation continues.'
            END IF
            IF (NFRVAL .GT. 0) THEN
               IF (NFRVLT .EQ. 0) NFRVLT = NFRVAL
               J = MIN(MXFR,NFRVAL)
               WRITE (LUPRI,'(A,I4/A,5F10.6:/,(28X,5F10.6))')
     &            ' Number of frequencies     :',NFRVLT,
     &            ' Frequencies               :',
     &            (FRVAL(I), I=1,J)
             IF (NFRVLT .GT. MXFR) THEN
               WRITE(LUPRI,'(/A/A,I5/A)')
     &         'You have asked for too many frequencies under *ABALNR',
     &         'Current maximum is MXFR =',MXFR,
     &         'Reduce number of frequencies or increase'//
     &         ' MXFR in cbilnr.h and recompile.'
               CALL QUIT(
     &         'You have asked for too many frequencies under *ABALNR')
             END IF
            END IF
            IF (STATIC) WRITE(LUPRI,'(A)') 
     &         ' (Also) using static approximation for'//
     &         ' optical activity G tensor.'
            IF (ABSORP) THEN
               WRITE(LUPRI,'(A,F9.6,A)') 
     &              ' Damping parameter equals  :',DAMPING, ' a.u.'
            END IF
            WRITE (LUPRI,'(A,I4)')
     &         ' Print level                    :',IPRLNR
            WRITE (LUPRI,'(A,1P,D10.2)')
     &         ' Lin Resp convergence threshold :',THCLNR
            WRITE (LUPRI,'(A,I4)')
     &         ' Max. Lin Resp iterations       :',MAXITE
            IF (OOTV) WRITE (LUPRI,'(A)')
     &         ' Optimal orbital microiterations used.'
            IF (CUT) THEN
               WRITE (LUPRI,'(/A)') ' Program is stopped after LNRABA.'
            END IF
         END IF
      END IF
C
      RETURN
      END
C  /* Deck lnrini */
      SUBROUTINE LNRINI
C
C     Initialize /LNREXC/
C
#include <implicit.h>
#include <mxcent.h>
#include <cbiexc.h>
#include <cbilnr.h>
#include <abainf.h>
C
      IPR1IN   = IPRDEF
      IPRLNR   = IPRDEF
      SKIP     = .FALSE.
      CUT      = .FALSE.
      OOTV     = .FALSE.
      STATIC   = .FALSE.
      ALFA     = ABA_ALPHA
      NFRVAL   = 1
      FRVAL(1) = 0.0D0
      THCLNR   = 5.D-05
      MAXITE   = 60
      MXRM     = 400
      MXPHP    = 0
CSPAS:20/3-2011: allowing for perturbations in only one irrep
      VIBGIR = .FALSE.
CKeinSPASmehr
      RETURN
      END
C  /* Deck lnraba */
      SUBROUTINE LNRABA(POLDD,POLDQ,POLDL,POLDA,POLVL,POLVV,FOVIBG,
     &                  WORK,LWORK,PASS)
#include <implicit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <iratdef.h>
#include <priunit.h>
#include <cbilnr.h>
      LOGICAL PASS
      LOGICAL CICLC, HFCLC, TRIPLE, EXECLC, FOUND, CONV
      DIMENSION WORK(LWORK),SNDPRP(2)
      DIMENSION POLDD(2,3,3,MXFR), POLDQ(2,3,3,3,MXFR)
      DIMENSION POLDL(2,3,3,MXFR), POLDA(2,3,3,MXFR)
      DIMENSION POLVL(2,3,3,MXFR), POLVV(2,3,3,MXFR)
      DIMENSION POLDLS(3,3), POLDAS(3,3)
CSPAS:6/10-08: trying to add mass independent vibrational g-factor
C     DIMENSION FOVIBG(3*NATOMS,3*NATOMS,NSYM,MXFR)
      DIMENSION FOVIBG(3*NATOMS+3,3*NATOMS+3,NSYM,MXFR)
CKeinSPASmehr
      CHARACTER*8 LABEL1, LABEL2, LABINT(3*MXCOOR), BLANK
      PARAMETER (DP5=0.5D0)
      PARAMETER (LOCDBG = 3)
C
#include <orgcom.h>
#include <cbiexc.h>
#include <inflin.h>
#include <infvar.h>
#include <infdim.h>
#include <inforb.h>
#include <nuclei.h>
#include <inftap.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <maxmom.h>
#include <maxaqn.h>
#include <symmet.h>
#include <abainf.h>
#include <gnrinf.h>
#include <infsop.h>
#include <absorp.h>
#include <pcmlog.h>
C
      IF (SKIP) RETURN
      CALL QENTER('LNRABA')
      CALL TIMER('START ',TIMEIN,TIMOUT)
C
      IF (IPRLNR .GE. 0) 
     &     CALL TITLER('Solving Linear Response Equations','*',118)
C
      IPRRSP = IPRLNR
C
C     Get reference state
C     ===================
C
C     1. Work Allocations:
C
      IF (ABASOP) THEN
         LUDV   = NORBT * NORBT
         LPVX   = LPVMAT
      ELSE
         LUDV   = N2ASHX
         LPVX   = 0
      ENDIF
      KFREE  = 1
      LFREE  = LWORK
C      
      CALL MEMGET('REAL',KCMO  ,NCMOT ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KUDV  ,LUDV  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPVX  ,LPVX  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KXINDX,LCINDX,WORK,KFREE,LFREE)
C
      KWORK1 = KFREE
      LWORK1 = LFREE
C
      CALL RD_SIRIFC('CMO',FOUND,WORK(KCMO),WORK(KFREE),LFREE)
      IF (.NOT.FOUND) CALL QUIT('LNRABA error: CMO not found on SIRIFC')
      IF (NASHT .GT. 0) THEN
         CALL RD_SIRIFC('DV',FOUND,WORK(KWORK1),WORK(KWORK1),LWORK1)
         IF (.NOT.FOUND)
     &      CALL QUIT('LNRABA error: DV not found on SIRIFC')
         CALL DSPTSI(NASHT,WORK(KWORK1),WORK(KUDV))
      END IF
C
      ISYM = 1
      CALL LNRVAR(ISYM,IPRLNR,WORK(KWORK1),LWORK1)
C
      CALL GETCIX(WORK(KXINDX),IREFSY,IREFSY,WORK(KWORK1),LWORK1,0)
C
C
C     SOPPA :
C
      IF (ABASOP) THEN
C
C        Initialize XINDX
C
         CALL DZERO(WORK(KXINDX),LCINDX)
C
C        Find address array's for SOPPA calculation
C
         CALL SET2SOPPA(WORK(KXINDX+KABSAD-1),WORK(KXINDX+KABTAD-1),
     *                  WORK(KXINDX+KIJSAD-1),WORK(KXINDX+KIJTAD-1),
     *                  WORK(KXINDX+KIJ1AD-1),WORK(KXINDX+KIJ2AD-1),
     *                  WORK(KXINDX+KIJ3AD-1),WORK(KXINDX+KIADR1-1))
C
C
         REWIND (LUSIFC)
         IF (CCPPA) THEN
            CALL MOLLAB('CCSDINFO',LUSIFC,LUPRI)
         ELSE
            CALL MOLLAB('MP2INFO ',LUSIFC,LUPRI)
         ENDIF
C
C        reads the MP2 or CCSD correlation coefficients into PV
C
         CALL READT (LUSIFC,LPVMAT,WORK(KPVX))
C
         IF (IPRLNR.GT.10) THEN
            IF (CCPPA) THEN
               WRITE(LUPRI,'(/A)')' EXCIT1 : CCSD correlation ',
     &                           'coefficients'
            ELSE
               WRITE(LUPRI,'(/A,A)')' EXCIT1 :',
     &                              ' MP2 correlation coefficients'
            ENDIF
            CALL OUTPUT(WORK(KPVX),1,LPVMAT,1,1,LPVMAT,1,1,LUPRI)
         END IF
C
C        reads the MP2 or CCSD second order one particle density matrix 
C
         CALL READT (LUSIFC,NORBT*NORBT,WORK(KUDV))
C
C        UDV contains the MP2 one-density. Remove the diagonal
C        contribution from the zeroth order. (Added in MP2FAC)
C
         IF (IPRLNR.GT.10) THEN
            IF (CCPPA) THEN
               WRITE(LUPRI,'(/A)')' RSPMC : CCSD density'
            ELSE
               WRITE(LUPRI,'(/A)')' RSPMC : MP2 density'
            END IF
            CALL OUTPUT(WORK(KUDV),1,NORBT*NORBT,1,1,NORBT*NORBT,1,1,
     &                  LUPRI)
         END IF
C
         CALL SOPUDV(WORK(KUDV))
      END IF
C
C     Construct property-integrals and write to LUPROP
C     ================================================
C
C     2. Work Allocations:
C
      KIDSYM = KWORK1
      KIDADR = KIDSYM + 9*MXCENT
      KWORK2 = KIDADR + 9*MXCENT
      LWORK2 = LWORK  - KWORK2
C
      NLBTOT = 0
C
      IPRINT=IPRLNR
      IF (ALFA .OR. ROAA .OR. ROAG) THEN
         NCOMP  = 0
         NPATOM = 0
Clf
         IF(PCM.AND.LOCFLD) THEN
            CALL GET1IN(DUMMY,'DIPLOC ',NCOMP,WORK(KWORK2),LWORK2,
     &           LABINT,WORK(KIDSYM),WORK(KIDADR),
     &           IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,DUMMY,IPRINT)
            NLAB = 3
            CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
         ELSE
            CALL GET1IN(DUMMY,'DIPLEN ',NCOMP,WORK(KWORK2),LWORK2,
     &           LABINT,WORK(KIDSYM),WORK(KIDADR),
     &           IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,DUMMY,IPRINT)
            NLAB = 3
            CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
         END IF
      ENDIF
      IF (ROAA) THEN
         NCOMP  = 0
         NPATOM = 0
         DIPTMX = DIPORG(1)
         DIPTMY = DIPORG(2)
         DIPTMZ = DIPORG(3)
         DIPORG(1) = 0.0D0
         DIPORG(2) = 0.0D0
         DIPORG(3) = 0.0D0
         CALL GET1IN(DUMMY,'THETA  ',NCOMP,WORK(KWORK2),LWORK2,
     &            LABINT,WORK(KIDSYM),WORK(KIDADR),
     &            IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,DUMMY,
     &            IPRINT)
         DIPORG(1) = DIPTMX
         DIPORG(2) = DIPTMY
         DIPORG(3) = DIPTMZ
         NLAB = 6
         CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
      END IF
C
      IF (ROAG) THEN
C
         CALL LABCOP(1,NLBTOT,'XLONMAG ',ISYMAX(1,2),LABAPP,LABSYM)
         CALL LABCOP(1,NLBTOT,'YLONMAG ',ISYMAX(2,2),LABAPP,LABSYM)
         CALL LABCOP(1,NLBTOT,'ZLONMAG ',ISYMAX(3,2),LABAPP,LABSYM)
C
         NCOMP  = 0
         NPATOM = 0
         CALL GET1IN(DUMMY,'ANGMOM ',NCOMP,WORK(KWORK2),LWORK2,
     &            LABINT,WORK(KIDSYM),WORK(KIDADR),
     &            IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,DUMMY,
     &            IPRINT)
         NLAB = 3
         CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
C
         IF (MVEOR) THEN
            NCOMP  = 0
            NPATOM = 0
            CALL GET1IN(DUMMY,'DIPVEL ',NCOMP,WORK(KWORK2),LWORK2,
     &                  LABINT,WORK(KIDSYM),WORK(KIDADR),
     &                  IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,DUMMY,
     &                  IPRINT)
            NLAB = 3
            CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
         END IF
C
      END IF
C
CSPAS:6/10-08: trying to add mass independent vibrational g-factor
C
      IF (VIB_G) THEN
         NCOMP  = 0
         NPATOM = 0
         CALL GET1IN(DUMMY,'DIPVEL ',NCOMP,WORK(KWORK2),LWORK2,
     &               LABINT,WORK(KIDSYM),WORK(KIDADR),
     &               IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,DUMMY,
     &               IPRINT)
         NLAB = 3
         CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
      END IF
C
CKeinSPASmehr
C
C     Start calculation of polarizabilities and raman optical activity.
C     =================================================================
C
      IF (ALFA .OR. ROAA .OR. ROAG) THEN
C
C     Set variables and logicals
C
      CICLC  = .FALSE.
      HFCLC  = NASHT .LE. 1
      TRIPLE = .FALSE.
      EXECLC = .FALSE.
      NABATY = 1
      TRPLET = .FALSE.
      NABAOP = 1
      BLANK  = '        '
C
      IF (NFRVAL.LE.0) GOTO 999
C
C     Zero the property complex tensors
C
      CALL DZERO(SNDPRP,2)
      CALL DZERO(POLDD,2*9*MXFR)
      CALL DZERO(POLDQ,2*27*MXFR)
      CALL DZERO(POLDL,2*9*MXFR)
      CALL DZERO(POLDA,2*9*MXFR)
      CALL DZERO(POLVL,2*9*MXFR)
      CALL DZERO(POLVV,2*9*MXFR)
      CALL DZERO(POLDLS,9)
      CALL DZERO(POLDAS,9)
C
      IF (NFRVAL.LE.0) GOTO 999
C
C     Start the calculations
C
C     tbp, jan. 2004: pass twice to get DIPVEL as first operator
C                     for velocity gauge optical rotation adjusted
C                     to ensure vanishing rotation at zero frequency
C                     (the modified velocity gauge).
C
      IFRZER = 0
      NFRSAV = NFRVAL
      IF (ROAG .AND. MVEOR .AND. NFRVAL.GT.0) THEN
         NUMRUN = 2
         IFRZER = 0
         XMINFR = 1.0D15
         DO IFRVAL = 1,NFRVAL
            TSTFR  = ABS(FRVAL(IFRVAL))
            IF (TSTFR.LT.1.0D-10 .AND. TSTFR.LT.XMINFR) THEN
               XMINFR = TSTFR
               IFRZER = IFRVAL
            END IF
         END DO
      ELSE
         NUMRUN = 1
      END IF
C
      DO IRUN = 1,NUMRUN
C
         LUSOVE = -1
         LUGDVE = -1
         LUREVE = -1
         CALL GPOPEN(LUSOVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
         CALL GPOPEN(LUGDVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
         CALL GPOPEN(LUREVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
         IF (IRUN .EQ. 2) THEN
            IF (IFRZER .EQ. 0) THEN
               IF (NFRVAL.GE.MXFR) THEN
                  CALL QUIT('Too many frequencies in LNRABA')
               ELSE
                  NFRVAL      = NFRVAL + 1
                  NFREQ_ALPHA = NFRVAL
                  IFRZER      = NFRVAL
                  FRVAL(IFRZER)      = 0.0D0
                  FREQ_ALPHA(IFRZER) = FRVAL(IFRZER)
               END IF
            END IF
         END IF
C
         DO ISYM=1,NSYM
CSPAS:20/3-2011: allowing for VIB_G perturbations in only one irrep
            IF (.NOT.(VIB_G .AND. VIBGIR .AND. (ISYM .NE. IRVIBG))) THEN
CKeinSPASmehr
            KSYMOP = ISYM
            DO IOPER=1,NOPER(KSYMOP)
               KOPER=IOPER
               LABEL1 = LABOP(IOPER,KSYMOP)
               IDIP = 0
               IF ((LABEL1 .EQ.'XDIPLEN').OR.(LABEL1 .EQ.'XDIPLOC')) 
     $              IDIP = 1
               IF ((LABEL1 .EQ.'YDIPLEN').OR.(LABEL1 .EQ.'YDIPLOC')) 
     $              IDIP = 2
               IF ((LABEL1 .EQ.'ZDIPLEN').OR.(LABEL1 .EQ.'ZDIPLOC')) 
     $              IDIP = 3
               IF ((LABEL1(2:7) .EQ.'DIPLEN') .OR. (LABEL1(2:7) .EQ.
     &                             'DIPLOC')) THEN
                  NABATY = 1
               ELSE IF ((LABEL1(2:7) .EQ. 'ANGMOM') .OR.
     &                  (LABEL1(2:7) .EQ. 'LONMAG')) THEN
                  NABATY = -1
               ELSE
                  WRITE (LUPRI,'(2A)') LABEL1, ' not recognised'
                  CALL QUIT('Unrecognised operator in LNRABA')
               END IF
               IF (IRUN .EQ. 2) THEN
                  IF (IDIP .EQ. 1) THEN
                     LABEL1 = 'XDIPVEL '
                  ELSE IF (IDIP .EQ. 2) THEN
                     LABEL1 = 'YDIPVEL '
                  ELSE IF (IDIP .EQ. 3) THEN
                     LABEL1 = 'ZDIPVEL '
                  ELSE
                     GO TO 300  ! cycle operator loop
                  END IF
                  NABATY = -1   ! DIPVEL is anti-Hermitian
               END IF
C     
               CALL LNRVAR(ISYM,IPRLNR,WORK(KWORK2),LWORK2)
C     
               KGD1   = KWORK1
               KWRKG1 = KGD1
               LWRKG1 = LWORK - KWRKG1
               KSLV   = KGD1 + 2*NVARPT
               KLAST  = KSLV + 4*NVARPT
               IF (KLAST.GT.LWORK) CALL STOPIT('LNRABA',' ',KLAST,LWORK)
               KWRK = KLAST
               LWRK = LWORK - KLAST + 1
C     
C     Find right hand side (gradient) of first operator and write to file
C     ===================================================================
C     
               CALL GETGPV(LABEL1,DUMMY,DUMMY,WORK(KCMO),WORK(KUDV),
     &              WORK(KPVX),WORK(KXINDX),ANTSYM,WORK(KWRKG1),LWRKG1)
               REWIND LUGDVE
               CALL WRITT(LUGDVE,2*NVARPT,WORK(KWRKG1))
               IF (IPRLNR.GE.LOCDBG) THEN
                  WRITE (LUPRI,'(A,A,A,D15.6)')
     &            'GP Vector, label: ',LABEL1,'  Norm: ',
     &            DSQRT(DDOT(2*NVARPT,WORK(KWRKG1),1,WORK(KWRKG1),1))
                  IF (IPRLNR.GT.10) THEN
                     CALL OUTPUT(WORK(KGD1),1,NVARPT,1,2,NVARPT,
     &                           2,1,LUPRI)
                  END IF
               END IF
C     
C     Calculate linear response vector and write to file
C     ==================================================
C     
               CALL ABARSP(CICLC,HFCLC,TRIPLE,OOTV,ISYM,EXECLC,
     &              FRVAL,NFRVAL,NABATY,NABAOP,LABEL1,LUGDVE,LUSOVE,
     &              LUREVE,THCLNR,MAXITE,IPRRSP,MXRM,MXPHP,
     &              WORK(KWRK),LWRK)
C     
C     Loop over the second property operators
C     =======================================
C     
               DO 200 IPRLBL = 1, NLBTOT
C     
C     Find label and symmetry of second operator
C     ==========================================
C
                  LABEL2 = LABAPP(IPRLBL)
                  KSYM   = LABSYM(IPRLBL)
C
C     Only ANGMOM in second run (for vel. and mod. vel. opt. rot.).
C     =============================================================
C
                  IF (IRUN.EQ.2 .AND. LABEL2(2:7).NE.'ANGMOM')
     &               GO TO 200
C     
C     If symmetry of first operator equals symmetry of
C     second operator, that is if ISYM = KSYM, then
C     ================================================
C     
                  IF (KSYM.EQ.ISYM) THEN
C
                     IF (IPRLNR.GE.LOCDBG) THEN
                        WRITE(LUPRI,'(/,A,A,A,A,A,/)')
     &                  '***** RESPONSE CALCULATION FOR ',LABEL1,' ',
     &                  LABEL2,' OPERATOR DOUBLE *****'
                        CALL FLSHFO(LUPRI)
                     END IF
C     
C     Find right hand side (gradient) for second operator
C     ===================================================
C     
                     CALL GETGPV(LABEL2,DUMMY,DUMMY,WORK(KCMO),
     &                    WORK(KUDV),WORK(KPVX),WORK(KXINDX),ANTSYM,
     &                    WORK(KWRKG1),LWRKG1)
C
                     IF (IPRLNR.GE.LOCDBG) THEN
                        WRITE (LUPRI,'(A,A,A,D15.6)')
     &                  'GP Vector, label: ',LABEL2,'  Norm: ',
     &                  DSQRT(DDOT(2*NVARPT,WORK(KGD1),1,WORK(KGD1),1))
                        IF (IPRLNR.GT.10) THEN
                           CALL OUTPUT(WORK(KGD1),1,NVARPT,1,2,NVARPT,
     &                                 2,1,LUPRI)
                        END IF
                     ENDIF
C     
C     Form second order properties SNDPRP
C     ===================================
C     
                     IF (.NOT.ABSORP) REWIND LUSOVE
                     DO 100 IFRVAL = 1,NFRVAL
                        IF (ABSORP) THEN
                           CALL GPOPEN(LURSP,'RSPVEC','OLD',' ',
     &                          ' ',IDUMMY,.FALSE.)
                           CALL REARSP(LURSP,4*NVARPT,WORK(KSLV),LABEL1,
     &                          'COMPLEX ',FRVAL(IFRVAL),0.0D0,
     &                          ISYM,0,THCLR_ABSORP,FOUND,CONV,ANTSYM)
                           CALL GPCLOSE(LURSP,'KEEP')
                           SNDPRP(1)=DDOT(2*NVARPT,WORK(KSLV),1,
     &                          WORK(KGD1),1)
                           SNDPRP(2)=DDOT(2*NVARPT,
     &                          WORK(KSLV+2*NVARPT),1,
     &                          WORK(KGD1),1)
                        ELSE
                           CALL READT(LUSOVE,2*NVARPT,WORK(KSLV))
                           SNDPRP(1)=DDOT(2*NVARPT,WORK(KSLV),1,
     &                          WORK(KGD1),1)
                        END IF
C
                        IF (IPRLNR.GE.LOCDBG) THEN
                           WRITE (LUPRI,'(A,I4,A,1P,D15.6)')
     &                     'Solution Vector no. ',IFRVAL,'  Norm: ',
     &                     DSQRT(DDOT(2*NVARPT,WORK(KSLV),1,
     &                                         WORK(KSLV),1))
                           IF (IPRLNR.GT.10) THEN
                              CALL OUTPUT(WORK(KSLV),1,NVARPT,1,2,
     &                                    NVARPT,2,1,LUPRI)
                           END IF
                        ENDIF
                        IF (IPRLNR.GE.(LOCDBG-1) .OR. IDIP.EQ.0) THEN
                           WRITE (LUPRI,'(/,A,F15.8)')
     &                          ' Frequency = ',FRVAL(IFRVAL)
                           WRITE (LUPRI,'(4A,F15.8)')
     &                          ' Second order property for ',
     &                          LABEL1,LABEL2,' = ',SNDPRP(1)
                        ENDIF
C     
C     Write properties into the various property matrices
C     ===================================================
C     Skip this section if first operator is not a dipole
C     operator (tbp, 2004).
C     
                        IF (IDIP.GT.0 .AND. IDIP.LT.4) THEN
                           IF (LABEL2(2:7).EQ.'DIPLEN') THEN
                              IF (LABEL2(1:1).EQ.'X') THEN
                                 POLDD(1,IDIP,1,IFRVAL) = SNDPRP(1)
                                 POLDD(2,IDIP,1,IFRVAL) = SNDPRP(2)
                              END IF
                              IF (LABEL2(1:1).EQ.'Y') THEN
                                 POLDD(1,IDIP,2,IFRVAL) = SNDPRP(1)
                                 POLDD(2,IDIP,2,IFRVAL) = SNDPRP(2)
                              END IF
                              IF (LABEL2(1:1).EQ.'Z') THEN
                                 POLDD(1,IDIP,3,IFRVAL) = SNDPRP(1)
                                 POLDD(2,IDIP,3,IFRVAL) = SNDPRP(2)
                              END IF
Clf
                           ELSE IF (LABEL2(2:7).EQ.'DIPLOC') THEN
                              IF (LABEL2(1:1).EQ.'X') THEN
                                 POLDD(1,IDIP,1,IFRVAL) = SNDPRP(1)
                                 POLDD(2,IDIP,1,IFRVAL) = SNDPRP(2)
                              END IF
                              IF (LABEL2(1:1).EQ.'Y') THEN
                                 POLDD(1,IDIP,2,IFRVAL) = SNDPRP(1)
                                 POLDD(2,IDIP,2,IFRVAL) = SNDPRP(2)
                              END IF
                              IF (LABEL2(1:1).EQ.'Z') THEN
                                 POLDD(1,IDIP,3,IFRVAL) = SNDPRP(1)
                                 POLDD(2,IDIP,3,IFRVAL) = SNDPRP(2)
                              END IF
Clf
                           ELSE IF (LABEL2(3:8).EQ.'THETA ') THEN
                              IF (LABEL2(1:2).EQ.'XX') THEN
                                 POLDQ(1,IDIP,1,1,IFRVAL) = SNDPRP(1)
                                 POLDQ(2,IDIP,1,1,IFRVAL) = SNDPRP(2)
                              END IF
                              IF (LABEL2(1:2).EQ.'XY') THEN
                                 POLDQ(1,IDIP,1,2,IFRVAL) = SNDPRP(1)
                                 POLDQ(2,IDIP,1,2,IFRVAL) = SNDPRP(2)
                              END IF
                              IF (LABEL2(1:2).EQ.'XY') THEN
                                 POLDQ(1,IDIP,2,1,IFRVAL) = SNDPRP(1)
                                 POLDQ(2,IDIP,2,1,IFRVAL) = SNDPRP(2)
                              END IF
                              IF (LABEL2(1:2).EQ.'XZ') THEN
                                 POLDQ(1,IDIP,1,3,IFRVAL) = SNDPRP(1)
                                 POLDQ(2,IDIP,1,3,IFRVAL) = SNDPRP(2)
                              END IF
                              IF (LABEL2(1:2).EQ.'XZ') THEN
                                 POLDQ(1,IDIP,3,1,IFRVAL) = SNDPRP(1)
                                 POLDQ(2,IDIP,3,1,IFRVAL) = SNDPRP(2)
                              END IF
                              IF (LABEL2(1:2).EQ.'YY') THEN
                                 POLDQ(1,IDIP,2,2,IFRVAL) = SNDPRP(1)
                                 POLDQ(2,IDIP,2,2,IFRVAL) = SNDPRP(2)
                              END IF
                              IF (LABEL2(1:2).EQ.'YZ') THEN
                                 POLDQ(1,IDIP,2,3,IFRVAL) = SNDPRP(1)
                                 POLDQ(2,IDIP,2,3,IFRVAL) = SNDPRP(2)
                              END IF
                              IF (LABEL2(1:2).EQ.'YZ') THEN
                                 POLDQ(1,IDIP,3,2,IFRVAL) = SNDPRP(1)
                                 POLDQ(2,IDIP,3,2,IFRVAL) = SNDPRP(2)
                              END IF
                              IF (LABEL2(1:2).EQ.'ZZ') THEN
                                 POLDQ(1,IDIP,3,3,IFRVAL) = SNDPRP(1)
                                 POLDQ(2,IDIP,3,3,IFRVAL) = SNDPRP(2)
                              END IF
                           ELSE IF (LABEL2(2:7).EQ.'LONMAG') THEN
C     
C     Multiply with minus the Bohr-magneton (-0.5) to create the magnetic 
C     dipole operator from the angular momentum operator
C     
                              IF (LABEL2(1:1).EQ.'X') THEN
                                 POLDL(1,IDIP,1,IFRVAL) = -DP5*SNDPRP(1)
                                 POLDL(2,IDIP,1,IFRVAL) = -DP5*SNDPRP(2)
                              END IF
                              IF (LABEL2(1:1).EQ.'Y') THEN
                                 POLDL(1,IDIP,2,IFRVAL) = -DP5*SNDPRP(1)
                                 POLDL(2,IDIP,2,IFRVAL) = -DP5*SNDPRP(2)
                              END IF
                              IF (LABEL2(1:1).EQ.'Z') THEN
                                 POLDL(1,IDIP,3,IFRVAL) = -DP5*SNDPRP(1)
                                 POLDL(2,IDIP,3,IFRVAL) = -DP5*SNDPRP(2)
                              END IF
                           ELSE IF (LABEL2(2:7).EQ.'ANGMOM') THEN
                              IF (IRUN .EQ. 2) THEN
                                 IF (IPRLNR.GT.(LOCDBG+10)) THEN
                                    WRITE(LUPRI,*)
     &                              '...property goes into POLVL'
                                 END IF
                                 IF (LABEL2(1:1).EQ.'X') THEN
                                    POLVL(1,IDIP,1,IFRVAL) =
     &                                                     DP5*SNDPRP(1)
                                    POLVL(2,IDIP,1,IFRVAL) =
     &                                                     DP5*SNDPRP(2)
                                 END IF
                                 IF (LABEL2(1:1).EQ.'Y') THEN
                                    POLVL(1,IDIP,2,IFRVAL) =
     &                                                     DP5*SNDPRP(1)
                                    POLVL(2,IDIP,2,IFRVAL) =
     &                                                     DP5*SNDPRP(2)
                                 END IF
                                 IF (LABEL2(1:1).EQ.'Z') THEN
                                    POLVL(1,IDIP,3,IFRVAL) =
     &                                                     DP5*SNDPRP(1)
                                    POLVL(2,IDIP,3,IFRVAL) =
     &                                                     DP5*SNDPRP(2)
                                 END IF
                              ELSE
                                 IF (LABEL2(1:1).EQ.'X') THEN
                                    POLDA(1,IDIP,1,IFRVAL) =
     &                                                    -DP5*SNDPRP(1)
                                    POLDA(2,IDIP,1,IFRVAL) =
     &                                                    -DP5*SNDPRP(2)
                                 END IF
                                 IF (LABEL2(1:1).EQ.'Y') THEN
                                    POLDA(1,IDIP,2,IFRVAL) =
     &                                                    -DP5*SNDPRP(1)
                                    POLDA(2,IDIP,2,IFRVAL) =
     &                                                    -DP5*SNDPRP(2)
                                 END IF
                                 IF (LABEL2(1:1).EQ.'Z') THEN
                                    POLDA(1,IDIP,3,IFRVAL) =
     &                                                    -DP5*SNDPRP(1)
                                    POLDA(2,IDIP,3,IFRVAL) =
     &                                                    -DP5*SNDPRP(2)
                                 END IF
                              END IF
C-tbp                      ELSE IF (IPRLNR.LE.2) THEN
C                    ... hjaaj mar 2004: this property is not
C                        predefined, should never happen, but
C                        we make sure the value is printed.
C                        (the value was printed above if IPRLNR.gt.2)
C-tbp                         WRITE (LUPRI,'(/A,F15.8,A/4A,F15.8)')
C-tbp&                             ' Frequency = ',FRVAL(IFRVAL),' au',
C-tbp&                             ' Second order property for ',
C-tbp&                             LABEL1,LABEL2,' = ',SNDPRP(1)
                           END IF
                        END IF
 100                 CONTINUE
                  ENDIF
 200           CONTINUE
C
 300           CONTINUE
C     
C     End loop over operator in this symmetry
C     
            END DO
C     
C     End loop over symmetries
C     
CSPAS:20/3-2011: allowing for VIB_G perturbations in only one irrep
            END IF  ! (.NOT.(VIB_G .AND. VIBGIR .AND. (ISYM .NE. IRVIBG)))
CKeinSPASmehr
         END DO
C
C     Close and delete files.
C
         IF (LUSOVE.GT.0) CALL GPCLOSE(LUSOVE,'DELETE')
         IF (LUGDVE.GT.0) CALL GPCLOSE(LUGDVE,'DELETE')
         IF (LUREVE.GT.0) CALL GPCLOSE(LUREVE,'DELETE')
C
C     End repetition loop (IRUN)
C
      END DO
C
      NFRVAL      = NFRSAV
      NFREQ_ALPHA = NFRVAL
      IF (ROAG .AND. MVEOR) THEN
         IF (IFRZER.LT.1 .OR. IFRZER.GT.MXFR) THEN
            CALL QUIT('Logical error in LNRABA')
         END IF
         DO IFRVAL = 1,NFRVAL
            DO IANG = 1,3
               DO IDIP = 1,3
                  DO IR = 1,2
                     POLVV(IR,IDIP,IANG,IFRVAL) =
     &                                        POLVL(IR,IDIP,IANG,IFRVAL)
     &                                      - POLVL(IR,IDIP,IANG,IFRZER)
                  END DO
               END DO
            END DO
         END DO
         IF ((IPRLNR.GE.LOCDBG) .AND. (NFRVAL.GT.0)) THEN
            DO IFRVAL = 1,NFRVAL
               WRITE(LUPRI,'(/,A,1P,D15.6)') 'Frequency: ',FRVAL(IFRVAL)
               WRITE(LUPRI,'(A)') 'POLDL: London'
               WRITE(LUPRI,'(A)') 'POLDA: No-London'
               WRITE(LUPRI,'(A)') 'POLVL: Velocity'
               WRITE(LUPRI,'(A)') 'POLVV: Modified velocity'
               WRITE(LUPRI,'(A,A,/,A,A)')
     &         ' Ei Bj      POLDL           POLDA           POLVL',
     &         '           POLVV',
     &         '-------------------------------------------------',
     &         '---------------------'
               DO IANG = 1,3
                  DO IDIP = 1,3
        WRITE(LUPRI,'(1X,I2,1X,I2,1X,F15.7,1X,F15.7,1X,F15.7,1X,F15.7)')
     &            IDIP,IANG,POLDL(1,IDIP,IANG,IFRVAL),
     &                      POLDA(1,IDIP,IANG,IFRVAL),
     &                      POLVL(1,IDIP,IANG,IFRVAL),
     &                      POLVV(1,IDIP,IANG,IFRVAL)
                     IF (ABSORP) THEN
        WRITE(LUPRI,'(7X,F15.7,1X,F15.7,1X,F15.7,1X,F15.7)')
     &                      POLDL(2,IDIP,IANG,IFRVAL),
     &                      POLDA(2,IDIP,IANG,IFRVAL),
     &                      POLVL(2,IDIP,IANG,IFRVAL),
     &                      POLVV(2,IDIP,IANG,IFRVAL)
                     END IF
                  END DO
               END DO
               WRITE(LUPRI,'(A,A)')
     &         '-------------------------------------------------',
     &         '---------------------'
            END DO
            WRITE(LUPRI,'(/,A,1P,D15.6)') 'Frequency: ',FRVAL(IFRZER)
            WRITE(LUPRI,'(A)') 'POLVL: Velocity'
            WRITE(LUPRI,'(A,/,A)')
     &      ' Ei Bj      POLVL',
     &      '----------------------'
            DO IANG = 1,3
               DO IDIP = 1,3
                  WRITE(LUPRI,'(1X,I2,1X,I2,1X,F15.7)')
     &            IDIP,IANG,POLVL(1,IDIP,IANG,IFRVAL)
                  IF (ABSORP) THEN
                     WRITE(LUPRI,'(7X,F15.7)')
     &               POLVL(2,IDIP,IANG,IFRVAL)
                  END IF
               END DO
            END DO
            WRITE(LUPRI,'(A)')
     &      '----------------------'
         END IF
      END IF
C
C     Skip point if no frequencies specified
C     STATIC = static approximation to optical rotation G tensor
      IF (STATIC)
     &   CALL LNRAST(POLDLS,POLDAS,WORK(KCMO),WORK(KWORK1),LWORK1)
C
 999  CONTINUE
C
C-tbp CALL GPCLOSE(LUSOVE,'DELETE')
C-tbp CALL GPCLOSE(LUGDVE,'DELETE')
C-tbp CALL GPCLOSE(LUREVE,'DELETE')
C
      END IF
C
C     Start calculation of vibrational g-factors.
C     ===========================================
C
      IF (VIB_G) THEN
         CALL HEADER('Calculating terms for vibrational g-factor',-1)
C
C        Set variables and logicals
C
         CICLC  = .FALSE.
         HFCLC  = NASHT .LE. 1
         TRIPLE = .FALSE.
         EXECLC = .FALSE.
         NABATY = -1
         LUSOVE = -1
         LUGDVE = -1
         LUREVE = -1
         NABAOP = 1
         BLANK  = '        '
         IFRVAL = 1 ! SPAS:8/10-08 there is no frequency for VIB_G
C
         CALL GPOPEN(LUSOVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
         CALL GPOPEN(LUGDVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
         CALL GPOPEN(LUREVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
C
CSPAS:6/10-08: trying to add mass independent vibrational g-factor
C        CALL DZERO(FOVIBG,3*NATOMS*3*NATOMS*NSYM*NFRVAL)
         CALL DZERO(FOVIBG,(3*NATOMS+3)*(3*NATOMS+3)*NSYM*NFRVAL)
CKeinSPASmehr
C
C        Loop over the first operators in the linear response function.
C        These are derivatives with respect to the nuclear cartesian 
C        coordinates.
C        ==============================================================
C
         NBR = 0
C
         DO ISYM = 1,NSYM
CSPAS:20/3-2011: allowing for VIB_G perturbations in only one irrep
            IF (.NOT.VIBGIR .OR. (ISYM .EQ. IRVIBG)) THEN
CKeinSPASmehr
C
            CALL LNRVAR(ISYM,IPRLNR,WORK(KWORK1),LWORK1)
C
C           3. Work Allocations:
C
            KGD1   = KWORK1
            KWRKG1 = KGD1
            LWRKG1 = LWORK - KWRKG1
            KSLV   = KGD1 + 2*NVARPT
            KLAST  = KSLV + 2*NVARPT
            IF (KLAST.GT.LWORK) CALL STOPIT('LNRABA',' ',KLAST,LWORK)
            KWRK   = KLAST
            LWRK   = LWORK - KLAST + 1
C
            DO IATOM = 1,NUCIND
               DO ICOOR = 1,3
C
C                 Determine the displacement coordinate number ISCOOR
C                 of symmetry ISYM. If the displacement coordinate
C                 does not have symmetry ISYM, ISCOOR is set to zero.
C                 =================================================== 
C
                  ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,ISYM-1,1)
C
                  IF (ISCOOR.GT.0) THEN
C
                     NBR = NBR + 1
C
                     IF (NBR .NE. ISCOOR) THEN
                        WRITE(LUPRI,'(A/A,2I8)') 
     &                  'Keld, why are NBR and ISCOOR different?',
     &                  'IATOM,NBR,ISCOOR',IATOM,NBR,ISCOOR
                     END IF
C
C                    Find right hand side for first operator 
C                    and write to file
C                    =======================================
C
                     CALL RDS2X(NBR,WORK(KUDV),WORK(KXINDX),
     &                          WORK(KWRKG1),LWRKG1)
                     REWIND LUGDVE
                     CALL WRITT(LUGDVE,2*NVARPT,WORK(KWRKG1))
C
                     IF (IPRLNR.GT.3) THEN
                        WRITE (LUPRI,'(2A)') 'GP Vector, label: ','d/dR'
                        CALL OUTPUT(WORK(KWRKG1),1,NVARPT,1,2,
     &                              NVARPT,2,1,LUPRI)
                     ENDIF
C
C                    Calculate eigenvector for first operator
C                    (i.e. derivatives with respect to the 
C                    nuclear cartesian coordinates) and write to file.
C                    =================================================
C
                     CALL ABARSP(CICLC,HFCLC,TRIPLE,OOTV,ISYM,EXECLC,
     &                           FRVAL,NFRVAL,NABATY,NABAOP,'d/dR    ',
     &                           LUGDVE,LUSOVE,LUREVE,THCLNR,MAXITE,
     &                           IPRRSP,MXRM,MXPHP,WORK(KWRK),LWRK)
C
C                    Loop over the second operators in the linear 
C                    response function.
C                    These are derivatives with respect to the nuclear 
C                    cartesian coordinates.
C                    =================================================
C
                     DO JATOM = 1,NUCIND
                        DO JCOOR = 1,3
C
C                       Determine the displacement coordinate number 
C                       JSCOOR of symmetry ISYM. If the displacement 
C                       coordinate does not have symmetry ISYM, 
C                       JSCOOR is set to zero.
C                       ============================================
C
                           JSCOOR = IPTCNT(3*(JATOM-1)+JCOOR,ISYM-1,1)
C
                           IF (JSCOOR.GT.0) THEN
C
C                             Find right hand side for second operator
C                             ========================================
C
                              CALL RDS2X(JSCOOR,WORK(KUDV),WORK(KXINDX),
     &                                   WORK(KWRKG1),LWRKG1)
C
                              IF (IPRLNR.GT.13) THEN
                                 WRITE (LUPRI,'(2A)') 
     &                               'GP Vector, label: ','d/dR'
                                 CALL OUTPUT(WORK(KWRKG1),1,2*NVARPT,1,
     &                                       1,2*NVARPT,1,1,LUPRI)
                              ENDIF

C                             Form the second order property: 
C                             vibrational g-factor.
C                             ===============================
C
                              REWIND LUSOVE
C
CSPAS:8/10-08 there is no frequency for VIB_G
C                             DO IFRVAL = 1,NFRVAL
CKeinSPASmehr
C
                                 CALL READT(LUSOVE,2*NVARPT,WORK(KSLV))
C
                                 FOVIBG(NBR,JSCOOR,ISYM,IFRVAL) 
     &                               = DDOT(2*NVARPT,WORK(KSLV),1,
     &                                      WORK(KGD1),1)
C
                                 IF (IPRLNR.GT.13) THEN
                                    WRITE (LUPRI,'(2A)') 
     &                                 'Solution Vector, label: ','d/dR'
                                    CALL OUTPUT(WORK(KSLV),1,2*NVARPT,1,
     &                                          1,2*NVARPT,1,1,LUPRI)
                                 ENDIF
C
                                 IF (IPRLNR.GT.2) THEN
CSPAS:8/10-08 there is no frequency for VIB_G
C                                   WRITE (LUPRI,'(A,F15.8)')
C    &                                  'Frequency = ',FRVAL(IFRVAL)
CKeinSPASmehr
                                    WRITE (LUPRI,'(4A,F15.8)')
     &                                  'Second order property for ',
     &                                  'd/dR','d/dR',' = ',
     &                                  FOVIBG(NBR,JSCOOR,ISYM,IFRVAL)
                                 ENDIF
C     
CSPAS:8/10-08 there is no frequency for VIB_G
C                             END DO
CKeinSPASmehr
C     
                           END IF
C     
                        END DO
                     END DO
C
CSPAS:6/10-08: trying to add mass independent vibrational g-factor
C
C                    Loop over dipole velocity operators
C                    ===================================
C     
                     DO JPRLBL = 1, NLBTOT
C     
C                       Find label and symmetry of operator
C                       ===================================
C
                        LABEL2 = LABAPP(JPRLBL)
                        KSYM   = LABSYM(JPRLBL)
C
C                       Check for dipole velocity operator
C                       ==================================
C
                        IF (LABEL2(2:7).EQ.'DIPVEL') THEN
C
                           IF (LABEL2(1:1).EQ.'X') JDIPV = 1
                           IF (LABEL2(1:1).EQ.'Y') JDIPV = 2
                           IF (LABEL2(1:1).EQ.'Z') JDIPV = 3
C     
C                          If symmetry of first operator equals symmetry of
C                          second operator, that is if ISYM = KSYM, then
C                          ================================================
C     
                           IF (KSYM.EQ.ISYM) THEN
C     
C                             Find right hand side (gradient) 
C                             for second operator
C                             ===============================
C     
                              CALL GETGPV(LABEL2,DUMMY,DUMMY,WORK(KCMO),
     &                                    WORK(KUDV),WORK(KPVX),
     &                                    WORK(KXINDX),ANTSYM,
     &                                    WORK(KWRKG1),LWRKG1)
C
                              IF (IPRLNR.GT.3) THEN
                                 WRITE (LUPRI,'(2A)') 
     &                               'GP Vector, label: ',LABEL2
                                 CALL OUTPUT(WORK(KWRKG1),1,2*NVARPT,1,
     &                                       1,2*NVARPT,1,1,LUPRI)
                              ENDIF
C
C                             Form the second order property: 
C                             << X/Y/ZDIPVEL ; d/dR >>
C                             ===============================
C
                              REWIND LUSOVE
C
CSPAS:8/10-08 there is no frequency for VIB_G
C                             DO IFRVAL = 1,NFRVAL
CKeinSPASmehr
C
                                 CALL READT(LUSOVE,2*NVARPT,WORK(KSLV))
C
                                 FOVIBG(NBR,3*NATOMS+JDIPV,ISYM,IFRVAL) 
     &                               = DDOT(2*NVARPT,WORK(KSLV),1,
     &                                      WORK(KGD1),1)
C
                                 IF (IPRLNR.GT.3) THEN
                                    WRITE (LUPRI,'(2A)') 
     &                                 'Solution Vector, label: ','d/dR'
                                    CALL OUTPUT(WORK(KSLV),1,2*NVARPT,1,
     &                                          1,2*NVARPT,1,1,LUPRI)
                                 ENDIF
C
                                 IF (IPRLNR.GT.2) THEN
CSPAS:8/10-08 there is no frequency for VIB_G
C                                   WRITE (LUPRI,'(A,F15.8)')
C    &                                  'Frequency = ',FRVAL(IFRVAL)
CKeinSPASmehr
                                    WRITE (LUPRI,'(4A,F15.8)')
     &                                  'Second order property for ',
     &                                  LABEL2,'d/dR',' = ',
     &                                  FOVIBG(NBR,3*NATOMS+JDIPV,ISYM,
     &                                         IFRVAL)
                                 ENDIF
C     
CSPAS:8/10-08 there is no frequency for VIB_G
C                             END DO
CKeinSPASmehr
C
                           ENDIF
C
                        ENDIF
C
                     END DO
C
CKeinSPASmehr
C
C     
                  END IF
C
               END DO
            END DO
C
CSPAS:6/10-08: trying to add mass independent vibrational g-factor
C
C           Loop over dipole velocity operators
C           ===================================
C     
            DO IPRLBL = 1, NLBTOT
C     
C              Find label and symmetry of operator
C              ===================================
C
               LABEL1 = LABAPP(IPRLBL)
               KSYMOP = LABSYM(IPRLBL)
C
C              Check for dipole velocity operator
C              ==================================
C
               IF (LABEL1(2:7).EQ.'DIPVEL') THEN
C
                  IF (LABEL1(1:1).EQ.'X') IDIPV = 1
                  IF (LABEL1(1:1).EQ.'Y') IDIPV = 2
                  IF (LABEL1(1:1).EQ.'Z') IDIPV = 3
C     
C                 If symmetry of first operator is correct 
C                 that is if ISYM = KSYMOP, then
C                 ========================================
C     
                  IF (KSYMOP.EQ.ISYM) THEN
C     
C                    Find right hand side (gradient) for first operator
C                    ==================================================
C     
                     CALL GETGPV(LABEL1,DUMMY,DUMMY,WORK(KCMO),
     &                           WORK(KUDV),WORK(KPVX),WORK(KXINDX),
     &                           ANTSYM,WORK(KWRKG1),LWRKG1)
                     REWIND LUGDVE
                     CALL WRITT(LUGDVE,2*NVARPT,WORK(KWRKG1))
C
                     IF (IPRLNR.GT.3) THEN
                        WRITE (LUPRI,'(2A)') 
     &                      'GP Vector, label: ',LABEL1
                        CALL OUTPUT(WORK(KWRKG1),1,2*NVARPT,1,
     &                              1,2*NVARPT,1,1,LUPRI)
                     ENDIF
C     
C                    Calculate linear response vector and write to file
C                    ==================================================
C     
                     CALL ABARSP(CICLC,HFCLC,TRIPLE,OOTV,ISYM,EXECLC,
     &                    FRVAL,NFRVAL,NABATY,NABAOP,LABEL1,LUGDVE,
     &                    LUSOVE,LUREVE,THCLNR,MAXITE,IPRRSP,MXRM,MXPHP,
     &                    WORK(KWRK),LWRK)
C
C                    Loop over the second operators in the linear 
C                    response function.
C                    These are derivatives with respect to the nuclear 
C                    cartesian coordinates.
C                    =================================================
C
                     DO JATOM = 1,NUCIND
                        DO JCOOR = 1,3
C
C                       Determine the displacement coordinate number 
C                       JSCOOR of symmetry ISYM. If the displacement 
C                       coordinate does not have symmetry ISYM, 
C                       JSCOOR is set to zero.
C                       ============================================
C
                           JSCOOR = IPTCNT(3*(JATOM-1)+JCOOR,ISYM-1,1)
C
                           IF (JSCOOR.GT.0) THEN
C
C                             Find right hand side for second operator
C                             ========================================
C
                              CALL RDS2X(JSCOOR,WORK(KUDV),WORK(KXINDX),
     &                                   WORK(KWRKG1),LWRKG1)
C
                              IF (IPRLNR.GT.3) THEN
                                 WRITE (LUPRI,'(2A)') 
     &                               'GP Vector, label: ','d/dR'
                                 CALL OUTPUT(WORK(KWRKG1),1,2*NVARPT,1,
     &                                       1,2*NVARPT,1,1,LUPRI)
                              ENDIF

C                             Form the second order property: 
C                             << d/dR ; X/Y/ZDIPVEL >>
C                             ===============================
C
                              REWIND LUSOVE
C
                              CALL READT(LUSOVE,2*NVARPT,WORK(KSLV))
C
                              FOVIBG(3*NATOMS+IDIPV,JSCOOR,ISYM,IFRVAL) 
     &                            = DDOT(2*NVARPT,WORK(KSLV),1,
     &                                   WORK(KGD1),1)
C
                              IF (IPRLNR.GT.3) THEN
                                 WRITE (LUPRI,'(2A)') 
     &                               'Solution Vector, label: ',LABEL1
                                 CALL OUTPUT(WORK(KSLV),1,2*NVARPT,1,
     &                                       1,2*NVARPT,1,1,LUPRI)
                              ENDIF
C
                              IF (IPRLNR.GT.2) THEN
                                 WRITE (LUPRI,'(4A,F15.8)')
     &                               'Second order property for ',
     &                               'd/dR',LABEL1,' = ',
     &                               FOVIBG(3*NATOMS+IDIPV,JSCOOR,ISYM,
     &                                      IFRVAL)
                              ENDIF
C     
                           END IF
C     
                        END DO
                     END DO
C
C                    Loop over dipole velocity operators
C                    ===================================
C     
                     DO JPRLBL = 1, NLBTOT
C     
C                       Find label and symmetry of operator
C                       ===================================
C
                        LABEL2 = LABAPP(JPRLBL)
                        KSYM   = LABSYM(JPRLBL)
C
C                       Check for dipole velocity operator
C                       ==================================
C
                        IF (LABEL2(2:7).EQ.'DIPVEL') THEN
C
                           IF (LABEL2(1:1).EQ.'X') JDIPV = 1
                           IF (LABEL2(1:1).EQ.'Y') JDIPV = 2
                           IF (LABEL2(1:1).EQ.'Z') JDIPV = 3
C     
C                          If symmetry of first operator equals symmetry of
C                          second operator, that is if ISYM = KSYM, then
C                          ================================================
C     
                           IF (KSYM.EQ.ISYM) THEN
C     
C                             Find right hand side (gradient) 
C                             for second operator
C                             ===============================
C     
                              CALL GETGPV(LABEL2,DUMMY,DUMMY,WORK(KCMO),
     &                                    WORK(KUDV),WORK(KPVX),
     &                                    WORK(KXINDX),ANTSYM,
     &                                    WORK(KWRKG1),LWRKG1)
C
                              IF (IPRLNR.GT.3) THEN
                                 WRITE (LUPRI,'(2A)') 
     &                               'GP Vector, label: ',LABEL2
                                 CALL OUTPUT(WORK(KWRKG1),1,2*NVARPT,1,
     &                                       1,2*NVARPT,1,1,LUPRI)
                              ENDIF
C
C                             Form the second order property: 
C                             << X/Y/ZDIPVEL ; X/Y/ZDIPVEL >>
C                             ===============================
C
                              REWIND LUSOVE
C
                              CALL READT(LUSOVE,2*NVARPT,WORK(KSLV))
C
                              FOVIBG(3*NATOMS+IDIPV,3*NATOMS+JDIPV,ISYM,
     &                               IFRVAL) 
     &                            = DDOT(2*NVARPT,WORK(KSLV),1,
     &                                   WORK(KGD1),1)
C
                              IF (IPRLNR.GT.3) THEN
                                 WRITE (LUPRI,'(2A)') 
     &                               'Solution Vector, label: ',LABEL1
                                 CALL OUTPUT(WORK(KSLV),1,2*NVARPT,1,
     &                                       1,2*NVARPT,1,1,LUPRI)
                              ENDIF
C
                              IF (IPRLNR.GT.2) THEN
                                 WRITE (LUPRI,'(4A,F15.8)')
     &                               'Second order property for ',
     &                               LABEL2,LABEL1,' = ',
     &                               FOVIBG(3*NATOMS+IDIPV,
     &                                      3*NATOMS+JDIPV,ISYM,
     &                                      IFRVAL)
                              ENDIF
C
                           ENDIF
C
                        ENDIF
C
                     END DO
C
                  ENDIF
C
               ENDIF
C
            END DO
C
CKeinSPASmehr
C
CSPAS:20/3-2011: allowing for perturbations in only one irrep
            END IF
CKeinSPASmehr
         END DO
C
         CALL GPCLOSE(LUSOVE,'DELETE')
         CALL GPCLOSE(LUGDVE,'DELETE')
         CALL GPCLOSE(LUREVE,'DELETE')
C
      END IF
C
      CALL TIMER ('LNRABA',TIMEIN,TIMOUT)
      PASS = .TRUE.
      IF (CUT) THEN
         WRITE (LUPRI,'(/A/A)')
     &     ' Program stopped after LNRABA as requested.',
     &     ' No restart file has been written.'
         CALL QUIT(' ***** End of ABACUS (.STOP in *ABALNR) *****')
      END IF
C
      CALL QEXIT('LNRABA')
      RETURN
      END
C
C  /* Deck lnrast */ 
      SUBROUTINE LNRAST(POLDLS,POLDAS,CMO,WORK,LWORK)
C
C     _ST_atic approximation to optical rotation G tensors.
C      
#include <implicit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <iratdef.h>
#include <priunit.h>
#include <cbilnr.h>
      LOGICAL FOUND, CONV
      DIMENSION WORK(LWORK), SNDPRP(2), CMO(*)
      DIMENSION POLDLS(3,3), POLDAS(3,3)
      DIMENSION DTOT(N2ORBX), DTOTSV(N2ORBX)
      DIMENSION DIPSLV(N2ORBX), ANGSLV(N2ORBX)
      CHARACTER*8 LABEL1, LABEL2, BLANK
      PARAMETER (DP5=0.5D0)
#include <cbiexc.h>
#include <inflin.h>
#include <infvar.h>
#include <infdim.h>
#include <inforb.h>
#include <nuclei.h>
#include <inftap.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <maxmom.h>
#include <maxaqn.h>
#include <symmet.h>
#include <abainf.h>
#include <gnrinf.h>
#include <infsop.h>
#include <absorp.h>
      BLANK  = '        '
C
C     Construct one-electron density matrix
C
      CALL MKDENS(DTOT,WORK,LWORK,IPRINT)
C
      DO ISYM=1,NSYM
         KSYMOP = ISYM
         DO IOPER=1,NOPER(KSYMOP) 
            KSLV1 = 1
            KSLV2 = KSLV1 + 2*NVARPT
            LABEL1 = LABOP(IOPER,KSYMOP)
            IF (LABEL1(2:7).EQ.'DIPLEN') THEN
               IF (LABEL1 .EQ.'XDIPLEN') IDIP = 1
               IF (LABEL1 .EQ.'YDIPLEN') IDIP = 2
               IF (LABEL1 .EQ.'ZDIPLEN') IDIP = 3
               CALL GPOPEN(LURSP,'RSPVEC','OLD',' ',
     &              ' ',IDUMMY,.FALSE.)
               CALL REARSP(LURSP,2*NVARPT,WORK(KSLV1),LABEL1,
     &              BLANK,0.0D0,0.0D0,ISYM,0,
     &              THCLNR,FOUND,CONV,ANTSYM)
               CALL GPCLOSE(LURSP,'KEEP')
               IF (.NOT.FOUND) THEN
                  WRITE (LUPRI,'(2A)') 'Solution Vector, label: '
     &                 ,LABEL1,' not found on RSPVEC'
                  CALL QUIT('Solution vector 1 not found on RSPVEC')
               ELSE    
                  IF (IPRLNR.GT.3) THEN
                     WRITE(LUPRI,'(2A)') 'Solution Vector, '//
     &                    'label: ',LABEL1
                     CALL OUTPUT(WORK(KSLV1),1,NVARPT,1,2,NVARPT,
     &                    2,1,LUPRI)
                  END IF
                  CALL UPWOP(WORK(KSLV1),DIPSLV,.TRUE.)
                  CALL DGEMM('N','T',NORBT,NORBT,NORBT,1.D0,
     &                 DTOT,NORBT,DIPSLV,NORBT,0.D0,DTOTSV,NORBT)
                  DO IOPER2=1,NOPER(KSYMOP)
                     LABEL2 = LABOP(IOPER2,KSYMOP)
                     IF (LABEL2(2:7).EQ.'ANGMOM') THEN
                        CALL GPOPEN(LURSP,'RSPVEC','OLD',' ',
     &                       ' ',IDUMMY,.FALSE.)
                        CALL REARSP(LURSP,2*NVARPT,WORK(KSLV2),
     &                       LABEL2,BLANK,0.0D0,0.0D0,ISYM,0,
     &                       THCLNR,FOUND,CONV,ANTSYM)
                        CALL GPCLOSE(LURSP,'KEEP')
                        IF (.NOT.FOUND) THEN
                           WRITE (LUPRI,'(2A)') 'Solution Vector, '//
     &                          'label: ',LABEL2,' not found on RSPVEC'
                           CALL QUIT('Solution vector 2 not found '//
     &                          'on RSPVEC')
                        ELSE
                           IF (IPRLNR.GT.3) THEN
                              WRITE(LUPRI,'(2A)') 'Solution '//
     &                             'Vector, label: ', LABEL2
                              CALL OUTPUT(WORK(KSLV2),1,NVARPT,1,2,
     &                             NVARPT,2,1,LUPRI)
                           END IF
                           CALL UPWOP(WORK(KSLV2),ANGSLV,.FALSE.)
                           SNDPRP(1)=DDOT(N2ORBX,DTOTSV,1,
     &                          ANGSLV,1)
                           IF (LABEL2(1:1).EQ.'X')
     &                          POLDAS(IDIP,1) = SNDPRP(1)
                           IF (LABEL2(1:1).EQ.'Y')
     &                          POLDAS(IDIP,2) = SNDPRP(1)  
                           IF (LABEL2(1:1).EQ.'Z')
     &                          POLDAS(IDIP,3) = SNDPRP(1)  
                        END IF   
                     ELSE IF (LABEL2(2:7).EQ.'LONMAG') THEN
                        CALL GPOPEN(LURSP,'RSPVEC','OLD',' ',
     &                       ' ',IDUMMY,.FALSE.)
                        CALL REARSP(LURSP,2*NVARPT,WORK(KSLV2),
     &                       LABEL2,BLANK,0.0D0,0.0D0,ISYM,0,
     &                       THCLNR,FOUND,CONV,ANTSYM)
                        CALL GPCLOSE(LURSP,'KEEP')
                        IF (.NOT.FOUND) THEN
                           WRITE (LUPRI,'(2A)') 'Solution Vector, '//
     &                          'label: ',LABEL2,' not found on RSPVEC'
                           CALL QUIT('Solution vector 2 not found '//
     &                          'on RSPVEC')
                        ELSE 
                           IF (IPRLNR.GT.3) THEN
                              WRITE(LUPRI,'(2A)') 'Solution '//
     &                             'Vector, label: ',LABEL2
                              CALL OUTPUT(WORK(KSLV2),1,NVARPT,1,2,
     &                             NVARPT,2,1,LUPRI)
                           END IF
                           CALL UPWOP(WORK(KSLV2),ANGSLV,.FALSE.)
                           SNDPRP(1)=DDOT(N2ORBX,DTOTSV,1,
     &                          ANGSLV,1)
                           IF (LABEL2(1:1).EQ.'X')
     &                          POLDLS(IDIP,1) = SNDPRP(1)
                           IF (LABEL2(1:1).EQ.'Y')
     &                          POLDLS(IDIP,2) = SNDPRP(1)  
                           IF (LABEL2(1:1).EQ.'Z')
     &                          POLDLS(IDIP,3) = SNDPRP(1)  
                        END IF
                     END IF
                  END DO
               END IF
            END IF
         END DO
      END DO
      WRITE(LUPRI,9003)
      WRITE(LUPRI,9007)
      WRITE(LUPRI,9009) 'Bx','By','Bz'
      WRITE(LUPRI,9008) 'Ex',(POLDLS(1,JDIP), JDIP=1,3)
      WRITE(LUPRI,9008) 'Ey',(POLDLS(2,JDIP), JDIP=1,3)
      WRITE(LUPRI,9008) 'Ez',(POLDLS(3,JDIP), JDIP=1,3)
C
      WRITE(LUPRI,9004)
      WRITE(LUPRI,9007)
      WRITE(LUPRI,9009) 'Bx','By','Bz'
      WRITE(LUPRI,9008) 'Ex',(POLDAS(1,JDIP), JDIP=1,3)
      WRITE(LUPRI,9008) 'Ey',(POLDAS(2,JDIP), JDIP=1,3)
      WRITE(LUPRI,9008) 'Ez',(POLDAS(3,JDIP), JDIP=1,3)         
      DO IFRVAL=1,NFRVAL
         WRITE(LUPRI,9005) FRVAL(IFRVAL)
         WRITE(LUPRI,9007)
         WRITE(LUPRI,9009) 'Bx','By','Bz'
         WRITE(LUPRI,9008) 'Ex',(POLDLS(1,JDIP)*FRVAL(IFRVAL),
     &        JDIP=1,3)
         WRITE(LUPRI,9008) 'Ey',(POLDLS(2,JDIP)*FRVAL(IFRVAL),
     &        JDIP=1,3)
         WRITE(LUPRI,9008) 'Ez',(POLDLS(3,JDIP)*FRVAL(IFRVAL),
     &        JDIP=1,3)   

         WRITE(LUPRI,9006) FRVAL(IFRVAL)
         WRITE(LUPRI,9007)
         WRITE(LUPRI,9009) 'Bx','By','Bz'
         WRITE(LUPRI,9008) 'Ex',(POLDAS(1,JDIP)*FRVAL(IFRVAL),
     &        JDIP=1,3)
         WRITE(LUPRI,9008) 'Ey',(POLDAS(2,JDIP)*FRVAL(IFRVAL),
     &        JDIP=1,3)
         WRITE(LUPRI,9008) 'Ez',(POLDAS(3,JDIP)*FRVAL(IFRVAL),
     &        JDIP=1,3)           
      END DO
 9003 FORMAT (//10X,'London    G tensor using static approximation')
 9004 FORMAT (//10X,'No-London G tensor using static approximation')
 9005 FORMAT (//10X,'London    G tensor stat. appr. for freq ',F12.6,
     &     ' au') 
 9006 FORMAT (/10X,'No-London G tensor stat. appr. for freq ',F12.6,
     &     ' au')
 9007 FORMAT (10X,'------------------------------------------------')
 9008 FORMAT (16X,A,2X,3F12.7)
 9009 FORMAT (27X,3(A,10X)/)
C
      RETURN
      END
C
C  /* Deck lnrvar */
      SUBROUTINE LNRVAR(ISYM,IPRINT,WORK,LWORK)
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
      DIMENSION WORK(LWORK)
#include <inflin.h>
#include <infvar.h>
#include <infrsp.h>
#include <inforb.h>
#include <wrkrsp.h>
#include <abainf.h>
C
      LOGICAL TRIPLE
C
C     Set variables for response modules
C     ==================================
C
      TRIPLE = .FALSE.
      CALL ABAVAR(ISYM,TRIPLE,IPRINT,WORK,LWORK)
      IREFSY = LSYMRF
      NCREF  = NCONRF
      KSYMST = LSYMST
      KSYMOP = LSYMPT
      KZWOPT = NWOPT
      IF ((NASHT .EQ. 0).AND.(.NOT.ABASOP)) THEN
         KZCONF = 0
      ELSE IF (NASHT .EQ. 1) THEN
         KZCONF = 1
      ELSE
         KZCONF = NCONST
      END IF
      KZVAR  = KZWOPT + KZCONF
      KZYWOP = 2*KZWOPT
      KZYCON = 2*KZCONF
      KZYVAR = 2*KZVAR
C
      RETURN
      END
C  /* Deck betagm */
      FUNCTION BETAGM(ALFA,GM)
C
C Calculate Beta(G')**2 =
C BETAGM = 0.5*(3*ALFA(I,J)*GM(I,J) - ALFA(I,I)*GM(J,J))
C
#include <implicit.h>
      PARAMETER ( D0 = 0.0D0  , DP5 = 0.5D0 , D3 = 3.0D0 )
      DIMENSION ALFA(3,3),GM(3,3)
C
      BETAGM = D0
      DO 100 I = 1,3
         DO 200 J = 1,3
            BETAGM = BETAGM +
     &               DP5*(D3*ALFA(I,J)*GM(I,J)-ALFA(I,I)*GM(J,J))
 200     CONTINUE
 100  CONTINUE
C
      RETURN
      END
C  /* Deck betaal */
      FUNCTION BETAAL(ALFA)
C
C Calculate Beta(alfa)**2 =
C BETAAL =  0.5 * (3*ALFA(I,J)*ALFA(I,J) - ALFA(I,I)*ALFA(J,J))
C
#include <implicit.h>
      PARAMETER ( D0 = 0.0D0  , DP5 = 0.5D0 , D3 = 3.0D0 )
      DIMENSION ALFA(3,3)
C
      BETAAL = D0
      DO 100 I = 1,3
         DO 200 J = 1,3
            BETAAL = BETAAL + DP5*(D3*ALFA(I,J)*ALFA(I,J)
     *                             -ALFA(I,I)*ALFA(J,J))
 200     CONTINUE
 100  CONTINUE
C
      RETURN
      END
C  /* Deck betaa */
      FUNCTION BETAA(ALFA,A,OMEGA)
C
C Calculate Beta(A)**2 =
C BETAA =  0.5*OMEGA*ALFA(I,J)*EPSILON(I,K,L)*A(K,L,J)
C
#include <implicit.h>
      PARAMETER ( D0 = 0.0D0 , DP5 = 0.5D0 , D1 = 1.0D0 , DM1 = -1.0D0)
      DIMENSION ALFA(3,3), A(3,3,3)
      DIMENSION EPSILO(3,3,3)
      DATA EPSILO / D0, D0,  D0, D0, D0, DM1, D0,  D1, D0,
     &              D0, D0,  D1, D0, D0, D0,  DM1, D0, D0,
     &              D0, DM1, D0, D1, D0, D0,  D0,  D0, D0/
C
      BETAA = D0
      DO 100 I = 1,3
         DO 200 J = 1,3
            DO 300 K = 1,3
               DO 400 L = 1,3
                  BETAA = BETAA + DP5*OMEGA*ALFA(I,J)*EPSILO(I,K,L)*
     *                                A(K,L,J)
 400           CONTINUE
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
C
      RETURN
      END
C  /* Deck alfmn */
      FUNCTION ALFMN(ALFA)
C
C Calculate AlfaMean =  ALFMN = ( (1/3) * ALFA(I,I) )
C
#include <implicit.h>
      PARAMETER ( D0 = 0.0D0, D3 = 3.0D0 )
      DIMENSION ALFA(3,3)
C
      ALFMN = D0
      DO 100 I = 1,3
         ALFMN = ALFMN + ALFA(I,I)
 100  CONTINUE
C
      ALFMN = ALFMN  / D3
C
      RETURN
      END
C  /* Deck gmmn */
      FUNCTION GMMN(GM)
C
C Calculate G' =  GMMN = ( (1/3) * GM(I,I) )
C
#include <implicit.h>
      PARAMETER ( D0 = 0.0D0, D3 = 3.0D0 )
      DIMENSION GM(3,3)
C
      GMMN = D0
      DO 100 I = 1,3
         GMMN = GMMN + GM(I,I)
 100  CONTINUE
C
      GMMN = GMMN  / D3
C
      RETURN
      END
C  /* Deck raminl */
      FUNCTION RAMINL(ALFMN,BETAAL)
C
C Calculate RamanIntensity =  RAMIN = (45 * ALFMN**2) + (4 * BETA**2)
C for linear pol. incident radiation perpendicular to scattering plane
C 
#include <implicit.h>
      PARAMETER ( D0 = 0.0D0, D4 = 4.0D0, D45 = 45.0D0 )
C
      RAMINL = (D45 * ALFMN**2) + (D4 * BETAAL)
C
      RETURN
      END
C  /* Deck depoll */
      FUNCTION DEPOLL(ALFMN,BETAAL)
C
C Calculate the Depolarization Ratio for rightangle scattering,
C linear polarized incident light, perpendicular (=parall/perpend)
C plane = DEPOLR = (3 * BETA**2) / (45 * ALFMN**2 + 4 * BETA**2)
C
#include <implicit.h>
      PARAMETER ( D0 = 0.0D0, D3 = 3.0D0, D4 = 4.0D0, D45 = 45.0D0 )
C
      A      = D3 * BETAAL
      B      = (D45 * ALFMN**2) + (D4 * BETAAL)
C
      IF ( B .GT. 0.D-6 ) THEN
         DEPOLL = A / B
      ELSE
         DEPOLL = D0
      ENDIF
C
      RETURN
      END
C  /* Deck raminn */
      FUNCTION RAMINN(ALFMN,BETAAL)
C
C Calculate RamanIntensity =  RAMIN = (45 * ALFMN**2) + (7 * BETA**2)
C for natural & circular pol. incident radiation 
C
#include <implicit.h>
      PARAMETER ( D0 = 0.0D0, D7 = 7.0D0, D45 = 45.0D0 )
C
      RAMINN = (D45 * ALFMN**2) + (D7 * BETAAL)
C
      RETURN
      END
C  /* Deck depoln */
      FUNCTION DEPOLN(ALFMN,BETAAL)
C
C Calculate the Depolarization Ratio for rightangle scattering, 
C naturel or circ. pol. incident light (=parall/perpend).
C plane = DEPOLR = (6 * BETA**2) / (45 * ALFMN**2 + 7 * BETA**2)
C
#include <implicit.h>
      PARAMETER ( D0 = 0.0D0, D6 = 6.0D0, D7 = 7.0D0, D45 = 45.0D0 )
C
      A      = D6 * BETAAL
      B      = (D45 * ALFMN**2) + (D7 * BETAAL)
C
      IF ( B .GT. 0.D-6 ) THEN
         DEPOLN = A / B
      ELSE
         DEPOLN = D0
      ENDIF
C
      RETURN
      END
C  /* Deck deltaz */
      FUNCTION DELTAZ(BETAGM,BETAA)
C
C Calculate the difference differential scattering cross section  
C in depolarized rightangle scattering =
C DELTAZ = 4/c ( 6*BETAGM - 2*BETAA)
C
#include <implicit.h>
#include <codata.h>
      PARAMETER ( D2 = 2.0D0, D4 = 4.0D0, D6 = 6.0D0 )
C
      DELTAZ = D4 / CVEL * (D6 * BETAGM - D2 * BETAA)
C
      RETURN
      END
C  /* Deck deltax */
      FUNCTION DELTAX(BETAGM,BETAA,ALFAMN,GMMN)
C
C Calculate the difference differential scattering cross section
C in polarized rightangle scattering =
C DELTAX = 4/c ( 45*ALFAMN*GMMN+7*BETAGM+BETAA)
C
#include <implicit.h>
#include <codata.h>
      PARAMETER ( D4 = 4.0D0, D7 = 7.0D0, D45 = 45.0D0)
C
      DELTAX = D4 / CVEL * (D45*ALFAMN*GMMN+D7*BETAGM+BETAA)
C
      RETURN
      END
C  /* Deck delta0 */
      FUNCTION DELTA0(BETAGM,BETAA,ALFAMN,GMMN)
C
C Calculate the difference differential scattering cross section
C in forward scattering with no analyzer =
C DELTA0 = 4/c (180*ALFAMN*GMMN+4*BETAGM-4*BETAA)
C
#include <implicit.h>
#include <codata.h>
      PARAMETER ( D4 = 4.0D0, D180 = 18.0D1 )
C
      DELTA0 = D4 / CVEL * (D180*ALFAMN*GMMN+D4*BETAGM-D4*BETAA)
C
      RETURN
      END
C  /* Deck deltab */
      FUNCTION DELTAB(BETAGM,BETAA)
C
C Calculate the difference differential scattering cross section
C in backward scattering with no analyzer =
C DELTAB = 4/c ( 24*BETAGM+8*BETAA)
C
#include <implicit.h>
#include <codata.h>
      PARAMETER ( D4 = 4.0D0, D8 = 8.0D0, D24 = 24.0D0 )
C
      DELTAB = D4 / CVEL * (D24*BETAGM + D8*BETAA)
C
      RETURN
      END
C  /* Deck cid */
      FUNCTION CID(DELTA,RMINN)
C
C Calculate the Circular Intensity Difference (CID) for all
C arragements. CID = - DELTA  / (2*RMINN)
C RMIN(0) and (180) with no analyzer = 2*RMIN.
C Rmin (depolarized) = DEPLN * RMINN
#include <implicit.h>
C
      PARAMETER ( D0 = 0.0D0, DP5 = 0.5D0 )
C
      IF ( RMINN .GT. 0.D-6 ) THEN
         CID = DELTA / RMINN * DP5
      ELSE
         CID = D0
      ENDIF
C
      RETURN
      END
C  /* Deck boltzk */
      FUNCTION BOLTZK(FREQ)
C
C Calculate the Boltzmannfactor for Spectra simulation 
C BOLTZK = 1/(1-EXP(-h*freq/k/T)), FREQ is in Hartree = 2Pi*freq
C
#include <implicit.h>
#include <codata.h>
      PARAMETER ( TSTAND = 298.15D0 )
C
      BOLTZK = 1-EXP(-(FREQ*AUTK)/TSTAND)
      RETURN
      END
C  /* Deck crossk */
      FUNCTION CROSSK(FRVAL,FREQ)
C
C Calculate the Constant for the differential scattering cross section
C K = 16/90*Pi**4*((freq0-freq)/C)**3*freq0/C
C Combined with the freq. factor out of PLACZEK approximation
C B = SQRT(1/(4Pi*freq))
C CROSSK = K*B**2
C ****************** FREQ is in Hartree = 2Pi*freq
C
#include <implicit.h>
#include <codata.h>
      PARAMETER ( D180 = 180.D0)
C
      CROSSK = (FRVAL-FREQ)**3*FRVAL/FREQ/D180/CVEL**3
      RETURN
      END
C  /* Deck lnrout */
      SUBROUTINE LNROUT(POLDD,POLDL,POLDA,POLVL,POLVV,FOVIBG,IPRINT,
     &                  WORK,LWORK)
C
#include <implicit.h>
#include <priunit.h>
#include <absorp.h>
#include <inforb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <codata.h>
      PARAMETER (FACTOR = (288.0D-30)*(PI**2)*XFMOL*(XTANG**4),
     &           D3 = 3.0D0)
#include <nuclei.h>
#include <symmet.h>
#include <abainf.h>
#include <cbilnr.h>
      DIMENSION POLDD(2,3,3,MXFR), POLDL(2,3,3,MXFR), POLDA(2,3,3,MXFR),
     &          POLVL(2,3,3,MXFR), POLVV(2,3,3,MXFR),
     &          BETANL(2), BETAL(2), ALPHNL(2), ALPHAL(2),
     &          BETAVL(2), BETAVV(2), ALPHVL(2), ALPHVV(2)
cmbh: helper strings
      character mlab1*8,mlab2*8,mlab3*8
cmbh end
      CHARACTER COORDI*21,COORDJ*21
CSPAS:6/10-08: trying to add mass independent vibrational g-factor
C     DIMENSION FOVIBG(3*NATOMS,3*NATOMS,NSYM,MXFR)
      DIMENSION FOVIBG(3*NATOMS+3,3*NATOMS+3,NSYM,MXFR)
CKeinSPASmehr
      DIMENSION WORK(LWORK)
#include <ibtfun.h>
C
      IF (ALFA .AND. IPRINT .GE. 0) THEN
         CALL TITLER('Frequency dependent polarizabilities','+',118)
         DO 100 IFRVAL = 1,NFRVAL
            WRITE(LUPRI,9001) FRVAL(IFRVAL)
            WRITE(LUPRI,9002)
            IF (ABSORP) WRITE(LUPRI,'(/,10X,A)') 'Real part:'
            WRITE(LUPRI,9003) 'Ex','Ey','Ez'
            WRITE(LUPRI,9004) 'Ex',(POLDD(1,1,JDIP,IFRVAL), JDIP=1,3)
            WRITE(LUPRI,9004) 'Ey',(POLDD(1,2,JDIP,IFRVAL), JDIP=1,3)
            WRITE(LUPRI,9004) 'Ez',(POLDD(1,3,JDIP,IFRVAL), JDIP=1,3)
cmbh: print out polarizabilities for MidasCpp
            mlab1='XDIPLEN'
            mlab2='YDIPLEN'
            mlab3='ZDIPLEN'
            call wripro(POLDD(1,1,1,IFRVAL),"SCF/DFT   ",2,
     *                  mlab1,mlab1,mlab1,mlab1,
     *                  FRVAL(IFRVAL),FRVAL(IFRVAL),FRVAL(IFRVAL),1,0,
     *                  0,0)
            call wripro(POLDD(1,1,2,IFRVAL),"SCF/DFT   ",2,
     *                  mlab1,mlab2,mlab1,mlab2,
     *                  FRVAL(IFRVAL),FRVAL(IFRVAL),FRVAL(IFRVAL),1,0,
     *                  0,0)
            call wripro(POLDD(1,1,3,IFRVAL),"SCF/DFT   ",2,
     *                  mlab1,mlab3,mlab1,mlab3,
     *                  FRVAL(IFRVAL),FRVAL(IFRVAL),FRVAL(IFRVAL),1,0,
     *                  0,0)
            call wripro(POLDD(1,2,2,IFRVAL),"SCF/DFT   ",2,
     *                  mlab2,mlab2,mlab2,mlab2,
     *                  FRVAL(IFRVAL),FRVAL(IFRVAL),FRVAL(IFRVAL),1,0,
     *                  0,0)
            call wripro(POLDD(1,2,3,IFRVAL),"SCF/DFT   ",2,
     *                  mlab2,mlab3,mlab2,mlab3,
     *                  FRVAL(IFRVAL),FRVAL(IFRVAL),FRVAL(IFRVAL),1,0,
     *                  0,0)
            call wripro(POLDD(1,3,3,IFRVAL),"SCF/DFT   ",2,
     *                  mlab3,mlab3,mlab3,mlab3,
     *                  FRVAL(IFRVAL),FRVAL(IFRVAL),FRVAL(IFRVAL),1,0,
     *                  0,0)

cmbh end
            IF (ABSORP) THEN
               WRITE(LUPRI,'(/,10X,A)') 'Imaginary part:'
               WRITE(LUPRI,9003) 'Ex','Ey','Ez'
               WRITE(LUPRI,9004) 'Ex',(POLDD(2,1,JDIP,IFRVAL), JDIP=1,3)
               WRITE(LUPRI,9004) 'Ey',(POLDD(2,2,JDIP,IFRVAL), JDIP=1,3)
               WRITE(LUPRI,9004) 'Ez',(POLDD(2,3,JDIP,IFRVAL), JDIP=1,3)
            END IF
            POLISO = (POLDD(1,1,1,IFRVAL) +
     &                POLDD(1,2,2,IFRVAL) +
     &                POLDD(1,3,3,IFRVAL)) / D3
            WRITE (LUPRI,'(/A,1P,G14.7)')
     &         'Isotropic polarizability: ',POLISO
  100    CONTINUE
      ENDIF
      IF (ROAG) THEN
         IF (ABSORP) THEN
            CALL TITLER('Optical Rotation and Electronic Circular '//
     &           'Dichroism','+',102)
         ELSE
            CALL TITLER('Optical rotation','+',118)
         END IF
         IF (IPRINT .GE. 0) THEN
            DO IFRVAL = 1,NFRVAL
               WRITE(LUPRI,9005) FRVAL(IFRVAL)
               WRITE(LUPRI,9007)
               IF (ABSORP) WRITE(LUPRI,'(/10X,A)') 'Real part:'
               WRITE(LUPRI,9009) 'Bx','By','Bz'
               WRITE(LUPRI,9008) 'Ex',(POLDL(1,1,JDIP,IFRVAL), JDIP=1,3)
               WRITE(LUPRI,9008) 'Ey',(POLDL(1,2,JDIP,IFRVAL), JDIP=1,3)
               WRITE(LUPRI,9008) 'Ez',(POLDL(1,3,JDIP,IFRVAL), JDIP=1,3)
               IF (ABSORP) THEN
                  WRITE(LUPRI,'(/10X,A)') 'Imaginary part:'
                  WRITE(LUPRI,9009) 'Bx','By','Bz'
                  WRITE(LUPRI,9008) 'Ex',(POLDL(2,1,JDIP,IFRVAL), 
     &                 JDIP=1,3)
                  WRITE(LUPRI,9008) 'Ey',(POLDL(2,2,JDIP,IFRVAL), 
     &                 JDIP=1,3)
                  WRITE(LUPRI,9008) 'Ez',(POLDL(2,3,JDIP,IFRVAL), 
     &                 JDIP=1,3)
               END IF
            END DO
            DO IFRVAL = 1,NFRVAL
               WRITE(LUPRI,9006) FRVAL(IFRVAL)
               WRITE(LUPRI,9007)
               IF (ABSORP) WRITE(LUPRI,'(/10X,A)') 'Real part:'
               WRITE(LUPRI,9009) 'Bx','By','Bz'
               WRITE(LUPRI,9008) 'Ex',(POLDA(1,1,JDIP,IFRVAL), JDIP=1,3)
               WRITE(LUPRI,9008) 'Ey',(POLDA(1,2,JDIP,IFRVAL), JDIP=1,3)
               WRITE(LUPRI,9008) 'Ez',(POLDA(1,3,JDIP,IFRVAL), JDIP=1,3)
               IF (ABSORP) THEN
                  WRITE(LUPRI,'(/10X,A)') 'Imaginary part:'
                  WRITE(LUPRI,9009) 'Bx','By','Bz'
                  WRITE(LUPRI,9008) 'Ex',(POLDA(2,1,JDIP,IFRVAL), 
     &                 JDIP=1,3)
                  WRITE(LUPRI,9008) 'Ey',(POLDA(2,2,JDIP,IFRVAL), 
     &                 JDIP=1,3)
                  WRITE(LUPRI,9008) 'Ez',(POLDA(2,3,JDIP,IFRVAL), 
     &                 JDIP=1,3)
               END IF
            END DO
            IF (MVEOR) THEN
               DO IFRVAL = 1,NFRVAL
                  TSTFR = ABS(FRVAL(IFRVAL))
                  IF (TSTFR .LT. 1.0D-8) THEN
                     WRITE(LUPRI,9010) FRVAL(IFRVAL)
                     WRITE(LUPRI,9007)
                     WRITE(LUPRI,'(10X,A)') '--- undefined ---'
                  ELSE
                     FRQINV = 1.0D0/FRVAL(IFRVAL)
                     CALL DSCAL(18,FRQINV,POLVL(1,1,1,IFRVAL),1)
                     WRITE(LUPRI,9010) FRVAL(IFRVAL)
                     WRITE(LUPRI,9007)
                     IF (ABSORP) WRITE(LUPRI,'(/10X,A)') 'Real part:'
                     WRITE(LUPRI,9009) 'Bx','By','Bz'
                     WRITE(LUPRI,9008) 'Ex',(POLVL(1,1,JDIP,IFRVAL),
     &                    JDIP=1,3)
                     WRITE(LUPRI,9008) 'Ey',(POLVL(1,2,JDIP,IFRVAL),
     &                    JDIP=1,3)
                     WRITE(LUPRI,9008) 'Ez',(POLVL(1,3,JDIP,IFRVAL),
     &                    JDIP=1,3)
                     IF (ABSORP) THEN
                        WRITE(LUPRI,'(/10X,A)') 'Imaginary part:'
                        WRITE(LUPRI,9009) 'Bx','By','Bz'
                        WRITE(LUPRI,9008) 'Ex',(POLVL(2,1,JDIP,IFRVAL),
     &                       JDIP=1,3)
                        WRITE(LUPRI,9008) 'Ey',(POLVL(2,2,JDIP,IFRVAL),
     &                       JDIP=1,3)
                        WRITE(LUPRI,9008) 'Ez',(POLVL(2,3,JDIP,IFRVAL),
     &                       JDIP=1,3)
                     END IF
                  END IF
               END DO
               DO IFRVAL = 1,NFRVAL
                  TSTFR = ABS(FRVAL(IFRVAL))
                  IF (TSTFR .LT. 1.0D-8) THEN
                     WRITE(LUPRI,9011) FRVAL(IFRVAL)
                     WRITE(LUPRI,9007)
                     WRITE(LUPRI,'(10X,A)') '--- undefined ---'
                  ELSE
                     FRQINV = 1.0D0/FRVAL(IFRVAL)
                     CALL DSCAL(18,FRQINV,POLVV(1,1,1,IFRVAL),1)
                     WRITE(LUPRI,9011) FRVAL(IFRVAL)
                     WRITE(LUPRI,9007)
                     IF (ABSORP) WRITE(LUPRI,'(/10X,A)') 'Real part:'
                     WRITE(LUPRI,9009) 'Bx','By','Bz'
                     WRITE(LUPRI,9008) 'Ex',(POLVV(1,1,JDIP,IFRVAL),
     &                    JDIP=1,3)
                     WRITE(LUPRI,9008) 'Ey',(POLVV(1,2,JDIP,IFRVAL),
     &                    JDIP=1,3)
                     WRITE(LUPRI,9008) 'Ez',(POLVV(1,3,JDIP,IFRVAL),
     &                    JDIP=1,3)
                     IF (ABSORP) THEN
                        WRITE(LUPRI,'(/10X,A)') 'Imaginary part:'
                        WRITE(LUPRI,9009) 'Bx','By','Bz'
                        WRITE(LUPRI,9008) 'Ex',(POLVV(2,1,JDIP,IFRVAL),
     &                       JDIP=1,3)
                        WRITE(LUPRI,9008) 'Ey',(POLVV(2,2,JDIP,IFRVAL),
     &                       JDIP=1,3)
                        WRITE(LUPRI,9008) 'Ez',(POLVV(2,3,JDIP,IFRVAL),
     &                       JDIP=1,3)
                     END IF
                  END IF
               END DO
            END IF
         END IF
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
         IF (ABSORP) THEN
            WRITE(LUPRI,'(2(/,10X,A))')
     &      'Optical rotation and Electronic Circular Dicroism summary',
     &      '---------------------------------------------------------'
            WRITE(LUPRI,'(/,10X,A)')
     &           'ORD is given in units:'//
     &           ' degrees/(dm g/cm^3)'
            WRITE(LUPRI,'(10X,A)')
     &           'CD extinction coefficient is given in units:'//
     &           ' L/(mol cm)'
         ELSE
            WRITE(LUPRI,'(2(/,10X,A))')
     &           'Optical rotation summary',
     &           '------------------------'
         END IF
         DO IFRVAL = 1, NFRVAL
            TSTFR = ABS(FRVAL(IFRVAL))
            IF (TSTFR .LT. 1.0D-14) THEN
               WRITE (LUPRI,'(/10X,A,F12.6,A)')
     &              'Frequency: ',FRVAL(IFRVAL),' au'
               BETAL(1)  = 0.0D0
               BETANL(1) = 0.0D0
               BETAL(2)  = 0.0D0
               BETANL(2) = 0.0D0
               BETAVL(1) = 0.0D0
               BETAVV(1) = 0.0D0
               BETAVL(2) = 0.0D0
               BETAVV(2) = 0.0D0
            ELSE
               WRITE (LUPRI,'(/10X,A,F12.6,A,F12.6,A)') 
     &              'Frequency: ',FRVAL(IFRVAL),' au  = ',
     &              XTNM/FRVAL(IFRVAL), ' nm'
               XXXDIV   = 1.0D0/(D3*FRVAL(IFRVAL))
               BETAL(1) =-(POLDL(1,1,1,IFRVAL) + POLDL(1,2,2,IFRVAL)
     &                     + POLDL(1,3,3,IFRVAL))*XXXDIV
               BETANL(1)=-(POLDA(1,1,1,IFRVAL) + POLDA(1,2,2,IFRVAL)
     &                     + POLDA(1,3,3,IFRVAL))*XXXDIV
               BETAL(2) =-(POLDL(2,1,1,IFRVAL) + POLDL(2,2,2,IFRVAL)
     &                     + POLDL(2,3,3,IFRVAL))*XXXDIV
               BETANL(2)=-(POLDA(2,1,1,IFRVAL) + POLDA(2,2,2,IFRVAL)
     &                     + POLDA(2,3,3,IFRVAL))*XXXDIV
               BETAVL(1)=-(POLVL(1,1,1,IFRVAL) + POLVL(1,2,2,IFRVAL)
     &                     + POLVL(1,3,3,IFRVAL))*XXXDIV
               BETAVV(1)=-(POLVV(1,1,1,IFRVAL) + POLVV(1,2,2,IFRVAL)
     &                     + POLVV(1,3,3,IFRVAL))*XXXDIV
               BETAVL(2)=-(POLVL(2,1,1,IFRVAL) + POLVL(2,2,2,IFRVAL)
     &                     + POLVL(2,3,3,IFRVAL))*XXXDIV
               BETAVV(2)=-(POLVV(2,1,1,IFRVAL) + POLVV(2,2,2,IFRVAL)
     &                     + POLVV(2,3,3,IFRVAL))*XXXDIV
            END IF
            XXXCON = FACTOT*FRVAL(IFRVAL)*FRVAL(IFRVAL)/TMASS
            ALPHAL(1)=XXXCON*BETAL(1)
            ALPHNL(1)=XXXCON*BETANL(1)
            ALPHAL(2)=XXXCON*BETAL(2)
            ALPHNL(2)=XXXCON*BETANL(2)
            ALPHVL(1)=XXXCON*BETAVL(1)
            ALPHVV(1)=XXXCON*BETAVV(1)
            ALPHVL(2)=XXXCON*BETAVL(2)
            ALPHVV(2)=XXXCON*BETAVV(2)
            IF (ABSORP) THEN
               WRITE(LUPRI,'(10X,A,F12.6,A,F12.2,A)')
     &              'Damping  : ',DAMPING,' au  = ',
     &              DAMPING*XTKAYS,' cm-1'
               IF (MVEOR) THEN
                  WRITE(LUPRI,'(8(/10X,2(A,F15.6)))')
     &              'Beta (London)                :',
     &              BETAL(1), '  + i ',BETAL(2),
     &              'Beta (No-London)             :',
     &              BETANL(1),'  + i ',BETANL(2),
     &              'Beta (Velocity)              :',
     &              BETAVL(1), '  + i ',BETAVL(2),
     &              'Beta (Modified Velocity)     :',
     &              BETAVV(1),'  + i ',BETAVV(2),
     &              'Optical rotation (London)    :',
     &              ALPHAL(1),'  + i ',ALPHAL(2),
     &              'Optical rotation (No-London) :',
     &              ALPHNL(1),'  + i ',ALPHNL(2),
     &              'Optical rotation (Velocity)  :',
     &              ALPHVL(1),'  + i ',ALPHVL(2),
     &              'Optical rotation (Mod. Vel.) :',
     &              ALPHVV(1),'  + i ',ALPHVV(2)
               ELSE
                  WRITE(LUPRI,'(2(/10X,2(A,F12.6)))') 
     &              'Beta (London)                  :',
     &              BETAL(1), '  + i ',BETAL(2),
     &              'Beta (No-London)               :',
     &              BETANL(1),'  + i ',BETANL(2)
                  WRITE(LUPRI,'(4(/10X,A,F16.6))') 
     &              'Optical rotation (London)      :', 
     &              ALPHAL(1),
     &              'Optical rotation (No-London)   :', 
     &              ALPHNL(1),
     &              'Circular dichroism (London)    :', 
     &              ALPHAL(2)*TMASS/3.2988D5,
     &              'Circular dichroism (No-London) :', 
     &              ALPHNL(2)*TMASS/3.2988D5
               END IF
            ELSE
               IF (MVEOR) THEN
                  WRITE(LUPRI,'(8(/10X,A,F15.6))')
     &              'Beta (London)                :',BETAL(1),
     &              'Beta (No-London)             :',BETANL(1),
     &              'Beta (Velocity)              :',BETAVL(1),
     &              'Beta (Modified Velocity)     :',BETAVV(1),
     &              'Optical rotation (London)    :',ALPHAL(1),
     &              'Optical rotation (No-London) :',ALPHNL(1),
     &              'Optical rotation (Velocity)  :',ALPHVL(1),
     &              'Optical rotation (Mod. Vel.) :',ALPHVV(1)
               ELSE
                  WRITE(LUPRI,'(4(/10X,A,F12.6))') 
     &              'Beta (London)                :',BETAL(1),
     &              'Beta (No-London)             :',BETANL(1),
     &              'Optical rotation (London)    :',ALPHAL(1),
     &              'Optical rotation (No-London) :',ALPHNL(1)
               END IF
            END IF
            IF (IPRINT .GE. 1) THEN
               WRITE(LUPRI,'(/A,1P,D16.9)')
     &         'Conversion, beta --> alpha: ',XXXCON
            END IF
         END DO
      ENDIF
C
      IF (ABSORP .AND. IPRINT.GE.0) THEN
         WRITE(LUPRI,'(/10X,A4,A6,A21,/10X,A)') 
     &        'Sym.','State','Excitation energy',
     &        '---------------------------------------------------'
         DO ISYM=1,NSYM
            IF (NOPER(ISYM).GT.0) THEN
               DO ISTATE=1,NEXCITED_STATES
                  WRITE(LUPRI,'(6X,2I6,F12.4,A,F9.3,A,F10.2,A)') 
     &                 ISYM,ISTATE,EXC_ENERGY(ISTATE,ISYM),' a.u.',
     &                 EXC_ENERGY(ISTATE,ISYM)*XTEV,' eV',
     &                 XTNM/EXC_ENERGY(ISTATE,ISYM),' nm'
               END DO
            END IF
         END DO
      END IF
C
      IF (VIB_G) THEN
         NCOOR  = 3*NUCDEP
         KCSTRA = 1
         KSCTRA = KCSTRA + NCOOR*NCOOR
         KLAST  = KSCTRA + NCOOR*NCOOR
         IF (KLAST .GT. LWORK) CALL STOPIT('LNROUT',' ',KLAST,LWORK)
         CALL TRACOR(WORK(KCSTRA),WORK(KSCTRA),1,NCOOR,0)
C
         CALL HEADER('Vibrational g-factor',-1)
CSPAS:8/10-08 there is no frequency for VIB_G
C        DO IFRVAL = 1,NFRVAL
         IFRVAL = 1
CSPAS       WRITE(LUPRI,9110) FRVAL(IFRVAL)
CKeinSPASmehr
            WRITE(LUPRI,9111)
            WRITE(LUPRI,'(16X,A,21X,A,5X,A,8X,A)')
     &         'Coordinates','Sym','(au)','(cm-1 u)'
            WRITE(LUPRI,'(4X,2A)') 
     &          ('--------------------------------------',I=1,2)
            DO ISYM = 1,NSYM
CSPAS:20/3-2011: allowing for VIB_G perturbations in only one irrep
               IF (.NOT.VIBGIR .OR. (ISYM .EQ. IRVIBG)) THEN
CKeinSPASmehr
               DO ISCOOR = 1,3*NATOMS
                  CALL GENCOR(WORK(KCSTRA),NCOOR,ISCOOR,COORDI,ICORSY)
                  IF (ICORSY .EQ. ISYM) THEN
                     DO JSCOOR = 1,3*NATOMS
                        CALL GENCOR(WORK(KCSTRA),NCOOR,JSCOOR,COORDJ,
     &                              JCORSY)
                        IF (JCORSY .EQ. ISYM) THEN
                           VAL = FOVIBG(ISCOOR,JSCOOR,ISYM,IFRVAL)
                           WRITE(LUPRI,'(5X,2A,I3,F12.6,F14.3)') 
     &                           COORDI, COORDJ, ISYM, 
     &                           VAL,VAL*XTKAYS/XFAMU
                        END IF
                     END DO
CSPAS:6/10-08: trying to add mass independent vibrational g-factor
                     DO JDIPV = 1,3
                        IF (JDIPV.EQ.1) THEN
                           COORDJ = "XDIPVEL"
                           JDIPSY = ISYMAX(1,1) + 1
                        ELSE IF (JDIPV.EQ.2) THEN
                           COORDJ = "YDIPVEL"
                           JDIPSY = ISYMAX(2,1) + 1
                        ELSE IF (JDIPV.EQ.3) THEN
                           COORDJ = "ZDIPVEL"
                           JDIPSY = ISYMAX(3,1) + 1
                        END IF
                        IF (JDIPSY .EQ. ISYM) THEN
                           WRITE(LUPRI,'(5X,2A,I3,F12.6)') 
     &                        COORDI, COORDJ, ISYM, 
     &                        FOVIBG(ISCOOR,3*NATOMS+JDIPV,ISYM,IFRVAL)
                        END IF
                     END DO
CKeinSPASmehr
                  END IF
               END DO
CSPAS:6/10-08: trying to add mass independent vibrational g-factor
               DO IDIPV = 1,3
                  IF (IDIPV.EQ.1) THEN
                     COORDI = "XDIPVEL"
                     IDIPSY = ISYMAX(1,1) + 1
                  ELSE IF (IDIPV.EQ.2) THEN
                     COORDI = "YDIPVEL"
                     IDIPSY = ISYMAX(2,1) + 1
                  ELSE IF (IDIPV.EQ.3) THEN
                     COORDI = "ZDIPVEL"
                     IDIPSY = ISYMAX(3,1) + 1
                  END IF
                  IF (IDIPSY .EQ. ISYM) THEN
                     DO JSCOOR = 1,3*NATOMS
                        CALL GENCOR(WORK(KCSTRA),NCOOR,JSCOOR,COORDJ,
     &                              JCORSY)
                        IF (JCORSY .EQ. ISYM) THEN
                           WRITE(LUPRI,'(5X,2A,I3,F12.6)') 
     &                        COORDI, COORDJ, ISYM, 
     &                        FOVIBG(3*NATOMS+IDIPV,JSCOOR,ISYM,IFRVAL)
                        END IF
                     END DO
                     DO JDIPV = 1,3
                        IF (JDIPV.EQ.1) THEN
                           COORDJ = "XDIPVEL"
                           JDIPSY = ISYMAX(1,1) + 1
                        ELSE IF (JDIPV.EQ.2) THEN
                           COORDJ = "YDIPVEL"
                           JDIPSY = ISYMAX(2,1) + 1
                        ELSE IF (JDIPV.EQ.3) THEN
                           COORDJ = "ZDIPVEL"
                           JDIPSY = ISYMAX(3,1) + 1
                        END IF
                        IF (JDIPSY .EQ. ISYM) THEN
                           WRITE(LUPRI,'(5X,2A,I3,F12.6)') 
     &                        COORDI, COORDJ, ISYM, 
     &                 FOVIBG(3*NATOMS+IDIPV,3*NATOMS+JDIPV,ISYM,IFRVAL)
                        END IF
                     END DO
                  END IF
               END DO
CKeinSPASmehr
CSPAS:20/3-2011: allowing for VIB_G perturbations in only one irrep
               END IF
CKeinSPASmehr
            END DO
            WRITE(LUPRI,9111)
CSPAS:8/10-08 there is no frequency for VIB_G
C        END DO
CKeinSPASmehr
      END IF
C
 9001 FORMAT (/10X,'Polarizability tensor for frequency ',F12.6,' au')
 9002 FORMAT (10X,'----------------------------------------------------'
     &       )
 9003 FORMAT (16X,3(A14,1X),/)
 9004 FORMAT (16X,A,2X,1P,3(G14.7,1X))
 9005 FORMAT (/10X,'London    G tensor for frequency ',F12.6,' au')
 9006 FORMAT (/10X,'No-London G tensor for frequency ',F12.6,' au')
 9007 FORMAT (10X,'-------------------------------------------------')
 9008 FORMAT (16X,A,2X,3F12.7)
 9009 FORMAT (27X,3(A,10X)/)
 9010 FORMAT (/10X,'Velocity  G tensor for frequency ',F12.6,' au')
 9011 FORMAT (/10X,'Mod. Vel. G tensor for frequency ',F12.6,' au')
 9110 FORMAT (14X,'Vibrational g-factors for frequency ',F12.6,' au',/)
 9111 FORMAT (4X,76('='))
C
      RETURN
      END
