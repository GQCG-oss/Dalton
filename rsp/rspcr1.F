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
      SUBROUTINE CRINP(WORD)
C
#include <implicit.h>
C
#include <priunit.h>
#include <gnrinf.h>
#include <infrsp.h>
#include "rspprp.h"
#include <inflr.h>
#include <infhyp.h>
#include <infpri.h>
#include <infspi.h>
#include <infcr.h>
#include <inftap.h>
C
      LOGICAL NEWDEF, BCDFRQ, STATIC
      PARAMETER ( NTABLE = 24 , ZERO = 0.0D0 , THRFRQ = 1.0D-14 )
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      CHARACTER*8 LABEL
C
      DATA TABLE /'.BPROP ', '.BFREQ ', '.CPROP ', '.CFREQ ',
     *            '.DPROP ', '.DFREQ ', '.APROP ', '.MAX IT',
     *            '.THCLR ', '.MAXITO', '.PRINT ', '.FREQUE',
     *            '.ISABCD', '.INVEXP', '.THG   ', '.DC-SHG',
     *            '.DC-KER', '.IDRI  ', 'xXXXXXX', '.THRNRM',
     *            '.DIPLEN', '.DIPLNX', '.DIPLNY', '.DIPLNZ'/
C
C  Set BCDFRQ = .TRUE. if frequencies specified for each operator
C  Used to test if special processes (THG,DC-SHG,DC-KERR,IDRI)
C  and indivual frequencies specified then quit.
C
      DATA BCDFRQ/.FALSE./
C
C READ IN  INPUT
C
      NEWDEF = (WORD .EQ. '*CUBIC ')
      ICHANG = 0
      IF (NEWDEF) THEN
         CRCAL = .TRUE.
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
                    GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,
     *                     15,16,17,18,19,20,21,22,23,24),I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' KEYWORD "',WORD,
     *            '" NOT RECOGNIZED IN CRINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT(' ILLEGAL KEYWORD IN CRINP ')
 1             CONTINUE
                  READ(LUCMD,'( A )')LABEL
                  BCROP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 2             CONTINUE
                  BCDFRQ = .TRUE.
                  READ (LUCMD,*) NBCRFR
                  IF (NBCRFR.LE.MBCRFR) THEN
                     READ (LUCMD,*) (BCRFR(J),J=1,NBCRFR)
                  ELSE
                     WRITE (LUPRI,'(3(/,A,I5),/)')
     *               ' NUMBER OF FREQUENCIES SPECIFIED    :',NBCRFR,
     *               ' IS GREATER THAN THE ALLOWED NUMBER :',MBCRFR,
     *               ' THE NUMBER IS RESET TO THE MAXIMUM :',MBCRFR
                     READ (LUCMD,*) (BCRFR(J),J=1,MBCRFR),
     *                              (FFFF,J=MBCRFR+1,NBCRFR)
                     NBCRFR = MBCRFR
                  END IF
               GO TO 100
 3             CONTINUE
                  READ(LUCMD,'( A )')LABEL
                  CCROP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 4             CONTINUE
                  BCDFRQ = .TRUE.
                  READ (LUCMD,*) NCCRFR
                  IF (NCCRFR.LE.MCCRFR) THEN
                     READ (LUCMD,*) (CCRFR(J),J=1,NCCRFR)
                  ELSE
                     WRITE (LUPRI,'(3(/,A,I5),/)')
     *               ' NUMBER OF FREQUENCIES SPECIFIED    :',NCCRFR,
     *               ' IS GREATER THAN THE ALLOWED NUMBER :',MCCRFR,
     *               ' THE NUMBER IS RESET TO THE MAXIMUM :',MCCRFR
                     READ (LUCMD,*) (CCRFR(J),J=1,MCCRFR),
     *                              (FFFF,J=MCCRFR+1,NCCRFR)
                     NCCRFR = MCCRFR
                  END IF
               GO TO 100
 5             CONTINUE
                  READ(LUCMD,'( A )')LABEL
                  DCROP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 6             CONTINUE
                  BCDFRQ = .TRUE.
                  READ (LUCMD,*) NDCRFR
                  IF (NDCRFR.LE.MDCRFR) THEN
                     READ (LUCMD,*) (DCRFR(J),J=1,NDCRFR)
                  ELSE
                     WRITE (LUPRI,'(3(/,A,I5),/)')
     *               ' NUMBER OF FREQUENCIES SPECIFIED    :',NDCRFR,
     *               ' IS GREATER THAN THE ALLOWED NUMBER :',MDCRFR,
     *               ' THE NUMBER IS RESET TO THE MAXIMUM :',MDCRFR
                     READ (LUCMD,*) (DCRFR(J),J=1,MDCRFR),
     *                              (FFFF,J=MDCRFR+1,NDCRFR)
                     NDCRFR = MDCRFR
                  END IF
               GO TO 100
 7             CONTINUE
                  READ(LUCMD,'( A )')LABEL
                  ACROP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 8             CONTINUE
                  READ (LUCMD,*) MAXITL
               GO TO 100
 9             CONTINUE
                  READ (LUCMD,*) THCLR
               GO TO 100
 10            CONTINUE
                  READ (LUCMD,*) MAXITO
                  IF (MAXITO .GT. 0) OPTORB = .TRUE.
               GO TO 100
 11            CONTINUE
                  READ (LUCMD,*) IPRCR
               GO TO 100
 12            CONTINUE
                  READ (LUCMD,*) NBCRFR
                  IF (NBCRFR.LE.MBCRFR) THEN
                     READ (LUCMD,*) (BCRFR(J),J=1,NBCRFR)
                  ELSE
                     WRITE (LUPRI,'(3(/,A,I5),/)')
     *               ' NUMBER OF FREQUENCIES SPECIFIED    :',NBCRFR,
     *               ' IS GREATER THAN THE ALLOWED NUMBER :',MBCRFR,
     *               ' THE NUMBER IS RESET TO THE MAXIMUM :',MBCRFR
                     READ (LUCMD,*) (BCRFR(J),J=1,MBCRFR),
     *                              (FFFF,J=MBCRFR+1,NBCRFR)
                     NBCRFR = MBCRFR
                  END IF
               GO TO 100
 13            CONTINUE
                  READ (LUCMD,*) ISPINA,ISPINB,ISPINC,ISPIND
               GO TO 100
 14            CONTINUE
C                 test option: construct Hessian explicitly
C                 and solve linear equations by matrix inversion
                  INVEXP = .TRUE.
               GO TO 100
 15            CONTINUE
                  CRTHG  = .TRUE.
                  CRSPEC = .TRUE.
               GO TO 100
 16            CONTINUE
                  CRSHG  = .TRUE.
                  CRSPEC = .TRUE.
               GO TO 100
 17            CONTINUE
                  CRKERR = .TRUE.
                  CRSPEC = .TRUE.
               GO TO 100
 18            CONTINUE
                  CRIDRI = .TRUE.
                  CRSPEC = .TRUE.
               GO TO 100
 19            CONTINUE
               GO TO 100
 20            CONTINUE
                  READ (LUCMD,*) THRNRM
               GO TO 100
 21            CONTINUE
                  GAMALL = .FALSE.
                  LABEL='XDIPLEN'
                  ACROP( INDPRP(LABEL)) = .TRUE.
                  BCROP( INDPRP(LABEL)) = .TRUE.
                  CCROP( INDPRP(LABEL)) = .TRUE.
                  DCROP( INDPRP(LABEL)) = .TRUE.
                  LABEL='YDIPLEN'
                  ACROP( INDPRP(LABEL)) = .TRUE.
                  BCROP( INDPRP(LABEL)) = .TRUE.
                  CCROP( INDPRP(LABEL)) = .TRUE.
                  DCROP( INDPRP(LABEL)) = .TRUE.
                  LABEL='ZDIPLEN'
                  ACROP( INDPRP(LABEL)) = .TRUE.
                  BCROP( INDPRP(LABEL)) = .TRUE.
                  CCROP( INDPRP(LABEL)) = .TRUE.
                  DCROP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 22            CONTINUE
                  GAMALL = .FALSE.
                  LABEL='XDIPLEN'
                  ACROP( INDPRP(LABEL)) = .TRUE.
                  BCROP( INDPRP(LABEL)) = .TRUE.
                  CCROP( INDPRP(LABEL)) = .TRUE.
                  DCROP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 23            CONTINUE
                  GAMALL = .FALSE.
                  LABEL='YDIPLEN'
                  ACROP( INDPRP(LABEL)) = .TRUE.
                  BCROP( INDPRP(LABEL)) = .TRUE.
                  CCROP( INDPRP(LABEL)) = .TRUE.
                  DCROP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
 24            CONTINUE
                  GAMALL = .FALSE.
                  LABEL='ZDIPLEN'
                  ACROP( INDPRP(LABEL)) = .TRUE.
                  BCROP( INDPRP(LABEL)) = .TRUE.
                  CCROP( INDPRP(LABEL)) = .TRUE.
                  DCROP( INDPRP(LABEL)) = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' PROMPT "',WORD,
     *            '" NOT RECOGNIZED IN CRINP.'
               CALL QUIT(' ILLEGAL PROMPT IN CRINP ')
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
      NACRTO = 0
      NBCRTO = 0
      NCCRTO = 0
      NDCRTO = 0
      IF (ICHANG .GT. 0) THEN
         DO 500 I = 1,NPRLBL
            IF (ACROP(I)) NACRTO = NACRTO + 1
            IF (BCROP(I)) NBCRTO = NBCRTO + 1
            IF (CCROP(I)) NCCRTO = NCCRTO + 1
            IF (DCROP(I)) NDCRTO = NDCRTO + 1
 500     CONTINUE
         IF (BCDFRQ .AND. CRSPEC) THEN
            WRITE(LUPRI,'(/A)')
     *      ' ***ERROR -- Both b,c,d frequencies and a special process',
     *      ' have been specified'
            CALL QUIT
     *   ('CRINP: Individual frequencies and special process specified')
         ENDIF
      END IF
C
C  If CRSPEC we need to put the read frequencies into the appropriate
C  frequency arrays
C
      IF (CRSPEC) THEN
C        Find the static case, if included in input:
         JFR = 0
         DO IFR=1,NBCRFR
            IF (ABS(BCRFR(IFR)).LE.THRFRQ) JFR = IFR
         ENDDO
         STATIC = (JFR .GT. 0)
C
         NCCRFR = 0
         NDCRFR = 0
         IF (CRKERR .OR. STATIC) THEN
            NCCRFR = NCCRFR + 1
            CCRFR(NCCRFR) = ZERO
         ENDIF
         IF (CRKERR .OR. CRSHG .OR. STATIC) THEN
            NDCRFR = NDCRFR + 1
            DCRFR(NDCRFR) = ZERO
         ENDIF
         IF (CRSHG .OR. CRIDRI .OR. CRTHG) THEN
            DO 50 IFR=1,NBCRFR
               IF (IFR.EQ.JFR) GOTO 50
               NCCRFR = NCCRFR + 1
               CCRFR(NCCRFR) = BCRFR(IFR)
 50         CONTINUE
         ENDIF
         IF (CRIDRI) THEN
            DO 51 IFR=1,NBCRFR
               IF (IFR.EQ.JFR) GOTO 51
               NDCRFR = NDCRFR + 1
               DCRFR(NDCRFR) = - BCRFR(IFR)
 51         CONTINUE
         ENDIF
         IF (CRTHG) THEN
            DO 52 IFR=1,NBCRFR
               IF (IFR.EQ.JFR) GOTO 52
               NDCRFR = NDCRFR + 1
               DCRFR(NDCRFR) = BCRFR(IFR)
 52         CONTINUE
         ENDIF
      ENDIF
C
      NCRCAL = MIN(NACRTO,NBCRTO,NCCRTO,NDCRTO,NBCRFR,NCCRFR,NDCRFR)
      IF  (CRCAL .AND. NCRCAL.LE.0)  THEN
          WRITE (LUPRI,'(/A)') ' Cubic Response input ingnored because:'
         IF (NACRTO .EQ. 0) WRITE (LUPRI,'(A)')
     *      ' - no A operators requested'
         IF (NBCRTO .EQ. 0) WRITE (LUPRI,'(A)')
     *      ' - no B operators requested'
         IF (NBCRFR .EQ. 0) WRITE (LUPRI,'(A)')
     *      ' - no B frequencies requested'
         IF (NCCRTO .EQ. 0) WRITE (LUPRI,'(A)')
     *      ' - no C operators requested'
         IF (NCCRFR .EQ. 0) WRITE (LUPRI,'(A)')
     *      ' - no C frequencies requested'
         IF (NDCRTO .EQ. 0) WRITE (LUPRI,'(A)')
     *      ' - no D operators requested'
         IF (NDCRFR .EQ. 0) WRITE (LUPRI,'(A)')
     *      ' - no D frequencies requested'
          CRCAL = .FALSE.
      ELSE IF (CRCAL) THEN
         CALL HEADER(
     &      'Variable settings for Cubic Response calculation',0)
         WRITE (LUPRI,'(A,L2)')
     *   ' Second hyperpolarizability calculation :  CRCAL='
     *      ,CRCAL
         IF (CRSPEC) WRITE (LUPRI,'(A)')
     *   ' Calculation of a special processes carried out :'
         IF (CRTHG) WRITE (LUPRI,'(A)')
     *   ' - Third harmonic generation calculation'
         IF (CRSHG) WRITE (LUPRI,'(A)')
     *   ' - Dc-second harmonic generation calculation'
         IF (CRKERR) WRITE (LUPRI,'(A)')
     *   ' - Dc-Kerr calculation'
         IF (CRIDRI) WRITE (LUPRI,'(A)')
     *   ' - Intensity dependent refractive index calculation'
         WRITE(LUPRI,'(/I3,A,(1P,5D14.6))')
     *      NBCRFR,' B-frequencies', (BCRFR(I),I=1,NBCRFR)
         WRITE(LUPRI,'(I3,A,(1P,5D14.6))')
     *      NCCRFR,' C-frequencies', (CCRFR(I),I=1,NCCRFR)
         WRITE(LUPRI,'(I3,A,(1P,5D14.6))')
     *      NDCRFR,' D-frequencies', (DCRFR(I),I=1,NDCRFR)
         WRITE(LUPRI,'(/A,I4)')
     *      ' Print level                                    : IPRCR  ='
     *       ,IPRCR
         WRITE(LUPRI,'(A,1P,D10.3)')
     *      ' Threshold for non-zero norm of property vectors: THRNRM ='
     *       ,THRNRM
         WRITE(LUPRI,'(A,1P,D10.3)')
     *      ' Threshold for convergence of lin. resp. eq.s   : THCLR  ='
     *      ,THCLR
         WRITE(LUPRI,'(A,I4)')
     *      ' Maximum number of iterations in LR solver      : MAXITL ='
     *       ,MAXITL
         IF (OPTORB) WRITE(LUPRI,'(A,I4)')
     *      ' Maximum iterations in optimal orbital algorithm: MAXITO ='
     *      ,MAXITO
         IF (DIROIT) WRITE (LUPRI,'(A,L2)')
     *      ' Direct one-index transformation                : DIROIT ='
     *      ,DIROIT
      END IF
C
C *** END OF CRINP
C
      RETURN
      END
      SUBROUTINE CR1IND
C
C CALCULATE A POINTER TO THE NUMBER OF DIFFERENT
C ONE-INDEX LINEAR RESPONSE EQUATIONS THAT NEED TO BE SOLVED
C IN A SECOND HYPERPOLARIZABILITY CALCULATION
C
#include <implicit.h>
C
#include "rspprp.h"
#include <infcr.h>
#include <inforb.h>
#include <infrsp.h>
#include <infpri.h>
#include <infspi.h>
C
      PARAMETER ( ZERO = 0.0D0 )
C
      DO 250 IDSYM = 1,NSYM
      DO 300 ICSYM = 1,NSYM
      DO 400 IBSYM = 1,NSYM
         IASYM = MULD2H(IDSYM,MULD2H(ICSYM,IBSYM))
         IF ( (NDCROP(IDSYM) .GT. 0) .AND. (NCCROP(ICSYM) .GT. 0)
     *       .AND. (NBCROP(IBSYM) .GT. 0)
     *       .AND. (NACROP(IASYM).GT.0) ) THEN
            DO 750 IDOP = 1,NDCROP(IDSYM)
            DO 850 IDFR = 1,NDCRFR
               INUM = INCRLR(DCRLB(IDSYM,IDOP),DCRFR(IDFR),IDSYM)
 850        CONTINUE
 750        CONTINUE
            DO 500 ICOP = 1,NCCROP(ICSYM)
            DO 600 ICFR = 1,NCCRFR
               INUM = INCRLR(CCRLB(ICSYM,ICOP),CCRFR(ICFR),ICSYM)
 600        CONTINUE
 500        CONTINUE
            DO 700 IBOP = 1,NBCROP(IBSYM)
            DO 800 IBFR = 1,NBCRFR
               INUM = INCRLR(BCRLB(IBSYM,IBOP),BCRFR(IBFR),IBSYM)
 800        CONTINUE
 700        CONTINUE
C
C  If a special process is specified we only need some of the
C  combinations of the frequencies on the B,C and D operators
C
            IF (CRSPEC) THEN
               DO IBFR=1,NBCRFR
C
C  If dc-Kerr we need to solve linear response equations for w.
C
                  IF (CRKERR) THEN
                     DO IAOP=1,NACROP(IASYM)
                        INUM =
     *                     INCRLR(ACRLB(IASYM,IAOP),BCRFR(IBFR),IASYM)
                     ENDDO
                  ENDIF
C
C  If dc-SHG we need to solve linear response equations for 2w.
C
                  IF (CRSHG) THEN
                     DO IAOP=1,NACROP(IASYM)
                        INUM =
     *                    INCRLR(ACRLB(IASYM,IAOP),2*BCRFR(IBFR),IASYM)
                     ENDDO
                  ENDIF
C
C  If IDRI we need to solve linear response equations for w only.
C
                  IF (CRIDRI) THEN
                     DO IAOP=1,NACROP(IASYM)
                        INUM =
     *                     INCRLR(ACRLB(IASYM,IAOP),BCRFR(IBFR),IASYM)
                     ENDDO
                  ENDIF
C
C  If THG we need to solve linear response equations for 3w.
C
                  IF (CRTHG) THEN
                     DO IAOP=1,NACROP(IASYM)
                        INUM =
     *                    INCRLR(ACRLB(IASYM,IAOP),3*BCRFR(IBFR),IASYM)
                     ENDDO
                  END IF
               ENDDO
            ELSE
               DO 950 IDFR = 1,NDCRFR
               DO 900 ICFR = 1,NCCRFR
               DO 1000 IBFR = 1,NBCRFR
C
                  BCDFR = DCRFR(IDFR) + CCRFR(ICFR)
     *                  + BCRFR(IBFR)
                  DO 1100 IAOP = 1,NACROP(IASYM)
                     INUM = INCRLR(ACRLB(IASYM,IAOP),BCDFR,IASYM)
 1100             CONTINUE
 1000          CONTINUE
 900           CONTINUE
 950           CONTINUE
            ENDIF
         END IF
 400  CONTINUE
 300  CONTINUE
 250  CONTINUE
      RETURN
      END
      SUBROUTINE CR2IND
C
C CALCULATE A POINTER TO THE NUMBER OF DIFFERENT
C TWO-INDEX LINEAR RESPONSE EQUATIONS THAT NEED TO BE SOLVED
C IN A SECOND HYPERPOLARIZABILITY CALCULATION
C
#include <implicit.h>
C
#include "rspprp.h"
#include <infcr.h>
#include <inforb.h>
#include <infrsp.h>
#include <infpri.h>
#include <infspi.h>
C
      PARAMETER ( ZERO = 0.0D0 )
C
      DO 300 IDSYM = 1,NSYM
      DO 200 ICSYM = 1,NSYM
      DO 100 IBSYM = 1,NSYM
         IASYM = MULD2H(IDSYM,MULD2H(ICSYM,IBSYM))
         IF ( (NDCROP(IDSYM) .GT. 0) .AND. (NCCROP(ICSYM) .GT. 0)
     *       .AND. (NBCROP(IBSYM) .GT. 0)
     *       .AND. (NACROP(IASYM) .GT. 0) ) THEN
            DO 310 IDOP = 1,NDCROP(IDSYM)
            DO 210 ICOP = 1,NCCROP(ICSYM)
            DO 110 IBOP = 1,NBCROP(IBSYM)
C
C  If special processes specified we only need some combinations
C  of frequencies for the two index response vectors
C
               IF (CRSPEC) THEN
                  DO IBFR=1,NBCRFR
C
C  If dc-Kerr we need the combinations w,0 for b,c and b,d
C  and 0,0 for c,d.
C
                     IF (CRKERR) THEN
                        INUM =
     *                     INCR2(BCRLB(IBSYM,IBOP),CCRLB(ICSYM,ICOP),
     *                          BCRFR(IBFR),ZERO,IBSYM,ICSYM)
                        INUM =
     *                     INCR2(BCRLB(IBSYM,IBOP),DCRLB(IDSYM,IDOP),
     *                          BCRFR(IBFR),ZERO,IBSYM,IDSYM)
                        INUM =
     *                     INCR2(CCRLB(ICSYM,ICOP),DCRLB(IDSYM,IDOP),
     *                          ZERO,ZERO,ICSYM,IDSYM)
                     ENDIF
C
C  If dc-SHG we need the combinations w,w for b,c; w,0 for b,d and c,d.
C
                     IF (CRSHG) THEN
                        INUM =
     *                     INCR2(BCRLB(IBSYM,IBOP),CCRLB(ICSYM,ICOP),
     *                          BCRFR(IBFR),BCRFR(IBFR),IBSYM,ICSYM)
                        INUM =
     *                     INCR2(BCRLB(IBSYM,IBOP),DCRLB(IDSYM,IDOP),
     *                          BCRFR(IBFR),ZERO,IBSYM,IDSYM)
                        INUM =
     *                     INCR2(CCRLB(ICSYM,ICOP),DCRLB(IDSYM,IDOP),
     *                          BCRFR(IBFR),ZERO,ICSYM,IDSYM)
                     ENDIF
C
C  If IDRI we need the combinations w,w for b,c; w,-w for b,d and c,d.
C
                     IF (CRIDRI) THEN
                        INUM =
     *                     INCR2(BCRLB(IBSYM,IBOP),CCRLB(ICSYM,ICOP),
     *                          BCRFR(IBFR),BCRFR(IBFR),IBSYM,ICSYM)
                        INUM =
     *                     INCR2(BCRLB(IBSYM,IBOP),DCRLB(IDSYM,IDOP),
     *                          BCRFR(IBFR),-BCRFR(IBFR),IBSYM,IDSYM)
                        INUM =
     *                     INCR2(CCRLB(ICSYM,ICOP),DCRLB(IDSYM,IDOP),
     *                          BCRFR(IBFR),-BCRFR(IBFR),ICSYM,IDSYM)
                     ENDIF
C
C  If THG we need the combinations w,w for all operators
C
                     IF (CRTHG) THEN
                        INUM =
     *                     INCR2(BCRLB(IBSYM,IBOP),CCRLB(ICSYM,ICOP),
     *                          BCRFR(IBFR),BCRFR(IBFR),IBSYM,ICSYM)
                        INUM =
     *                     INCR2(BCRLB(IBSYM,IBOP),DCRLB(IDSYM,IDOP),
     *                          BCRFR(IBFR),BCRFR(IBFR),IBSYM,IDSYM)
                        INUM =
     *                     INCR2(CCRLB(ICSYM,ICOP),DCRLB(IDSYM,IDOP),
     *                          BCRFR(IBFR),BCRFR(IBFR),ICSYM,IDSYM)
                     ENDIF
                  ENDDO
               ELSE
                  DO 320 IDFR = 1,NDCRFR
                  DO 220 ICFR = 1,NCCRFR
                  DO 120 IBFR = 1,NBCRFR
                     INUM =
     *                  INCR2(BCRLB(IBSYM,IBOP),CCRLB(ICSYM,ICOP),
     *                       BCRFR(IBFR),CCRFR(ICFR),IBSYM,ICSYM)
                     INUM =
     *                  INCR2(BCRLB(IBSYM,IBOP),DCRLB(IDSYM,IDOP),
     *                       BCRFR(IBFR),DCRFR(IDFR),IBSYM,IDSYM)
                     INUM =
     *                  INCR2(CCRLB(ICSYM,ICOP),DCRLB(IDSYM,IDOP),
     *                       CCRFR(ICFR),DCRFR(IDFR),ICSYM,IDSYM)
120               CONTINUE
220               CONTINUE
320               CONTINUE
               ENDIF
110         CONTINUE
210         CONTINUE
310         CONTINUE
         END IF
100   CONTINUE
200   CONTINUE
300   CONTINUE
      RETURN
      END
      SUBROUTINE AN2SYM(NAOPNU,NBOPNU,NCOPNU,NDOPNU,IWRK)
C
#include <implicit.h>
C
#include <priunit.h>
#include <inforb.h>
#include <infrsp.h>
#include <infpri.h>
C
C THIS ROUTINE ANALYZE IF OPERATORS THAT ARE SPECIFIED IN CUBIC
C RESPONSE CALCULATIONS HAVE NO MATCHING PARTNERS DUE TO SYMMETRY.
C SUCH OPERATORS ARE REMOVED FROM THE CALCULATION IN THIS ROUTINE.
C
C BASED ON ANASYM
C
      DIMENSION NAOPNU(*),NBOPNU(*),NCOPNU(*),NDOPNU(*),IWRK(*)
C
      IF (IPRRSP.GT.40) THEN
         WRITE(LUPRI,'(/A)')
     *   ' AN2SYM: NUMBER OF PROPERTY VECTORS BEFORE REDUCTION '
         WRITE(LUPRI,'(/A)')' NAOPNU(I),I=1,NSYM '
         WRITE(LUPRI,'(8I5)')( NAOPNU(I),I=1,NSYM)
         WRITE(LUPRI,'(/A)')' NBOPNU(I),I=1,NSYM '
         WRITE(LUPRI,'(8I5)')( NBOPNU(I),I=1,NSYM)
         WRITE(LUPRI,'(/A)')' NCOPNU(I),I=1,NSYM '
         WRITE(LUPRI,'(8I5)')( NCOPNU(I),I=1,NSYM)
         WRITE(LUPRI,'(/A)')' NDOPNU(I),I=1,NSYM '
         WRITE(LUPRI,'(8I5)')( NDOPNU(I),I=1,NSYM)
      END IF
      KASYM = 1
      KBSYM = KASYM + NSYM
      KCSYM = KBSYM + NSYM
      KDSYM = KCSYM + NSYM
      DO 100 I = 1,NSYM
         IWRK(KASYM-1+I) = 0
         IWRK(KBSYM-1+I) = 0
         IWRK(KCSYM-1+I) = 0
         IWRK(KDSYM-1+I) = 0
 100  CONTINUE
      DO 200 IASYM = 1,NSYM
         IF ( NAOPNU(IASYM).LE.0 ) GO TO 200
         DO 300 IBSYM = 1,NSYM
            IF (NBOPNU(IBSYM).LE.0) GO TO 300
            DO 350 ICSYM = 1,NSYM
               IF (NCOPNU(ICSYM).LE.0) GO TO 350
               IDSYM = MULD2H(ICSYM,MULD2H(IASYM,IBSYM))
               IF(NDOPNU(IDSYM).GT.0) THEN
                  IWRK(KASYM-1+IASYM) = NAOPNU(IASYM)
                  IWRK(KBSYM-1+IBSYM) = NBOPNU(IBSYM)
                  IWRK(KCSYM-1+ICSYM) = NCOPNU(ICSYM)
                  IWRK(KDSYM-1+IDSYM) = NDOPNU(IDSYM)
               END IF
 350        CONTINUE
 300     CONTINUE
 200  CONTINUE
      NAOPTO = 0
      NBOPTO = 0
      NCOPTO = 0
      NDOPTO = 0
      DO 400 ISYM = 1,NSYM
         IF (IWRK(KASYM-1+ISYM).NE.NAOPNU(ISYM)) THEN
           WRITE(LUPRI,'(/I5,A,I5,A)')
     *     NAOPNU(ISYM),' A OPERATORS OF SYMMETRY',ISYM,' NOT INCLUDED'
           NAOPNU(ISYM) = IWRK(KASYM-1+ISYM)
         END IF
         IF (IWRK(KBSYM-1+ISYM).NE.NBOPNU(ISYM)) THEN
           WRITE(LUPRI,'(/I5,A,I5,A)')
     *     NBOPNU(ISYM),' B OPERATORS OF SYMMETRY',ISYM,' NOT INCLUDED'
           NBOPNU(ISYM) = IWRK(KBSYM-1+ISYM)
         END IF
         IF (IWRK(KCSYM-1+ISYM).NE.NCOPNU(ISYM)) THEN
           WRITE(LUPRI,'(/I5,A,I5,A)')
     *     NCOPNU(ISYM),' C OPERATORS OF SYMMETRY',ISYM,' NOT INCLUDED'
           NCOPNU(ISYM) = IWRK(KCSYM-1+ISYM)
         END IF
         IF (IWRK(KDSYM-1+ISYM).NE.NDOPNU(ISYM)) THEN
           WRITE(LUPRI,'(/I5,A,I5,A)')
     *     NDOPNU(ISYM),' D OPERATORS OF SYMMETRY',ISYM,' NOT INCLUDED'
           NDOPNU(ISYM) = IWRK(KDSYM-1+ISYM)
         END IF
         NAOPTO = NAOPTO + NAOPNU(ISYM)
         NBOPTO = NBOPTO + NBOPNU(ISYM)
         NCOPTO = NCOPTO + NCOPNU(ISYM)
         NDOPTO = NDOPTO + NDOPNU(ISYM)
 400  CONTINUE
      IF (IPRRSP.GT.40) THEN
         WRITE(LUPRI,'(/A)')
     *   ' AN2SYM: NUMBER OF PROPERTY VECTORS AFTER REDUCTION '
         WRITE(LUPRI,'(/A)')' NAOPNU(I),I=1,NSYM '
         WRITE(LUPRI,'(8I5)')( NAOPNU(I),I=1,NSYM)
         WRITE(LUPRI,'(/A)')' NBOPNU(I),I=1,NSYM '
         WRITE(LUPRI,'(8I5)')( NBOPNU(I),I=1,NSYM)
         WRITE(LUPRI,'(/A)')' NCOPNU(I),I=1,NSYM '
         WRITE(LUPRI,'(8I5)')( NCOPNU(I),I=1,NSYM)
         WRITE(LUPRI,'(/A)')' NDOPNU(I),I=1,NSYM '
         WRITE(LUPRI,'(8I5)')( NDOPNU(I),I=1,NSYM)
      END IF
      NOPTOT = NAOPTO + NBOPTO + NCOPTO + NDOPTO
      IF (NOPTOT.EQ.0) THEN
         WRITE(LUPRI,'(/A,/A)')
     *   ' AN2SYM: SYMMETRY OF SPECIFIED OPERATORS DOES NOT MATCH'
     *   ,' RESPONSE CALCULATION GIVE TRIVIALLY ZERO'
      END IF
      RETURN
      END
      INTEGER FUNCTION INCRLR(NEWLBL,FRQVAL,ISYM)
C
C DETERMINE IF A POINTER TO LINEAR RESPONSE EQUATION
C CHARACTERIZED BY NEWLBL AND FRQVAL HAS OCCURRED
C IF SO INCRLR GIVES ITS NUMBER IN THE LIST
C IF THE POINTER DOES NOT APPEAR IT IS ADDED TO
C LIST AND INCRLR GIVES THE ADDED NUMBER
C
C BASED ON INQRLR
#include <implicit.h>
C
#include <priunit.h>
#include <infrsp.h>
#include "rspprp.h"
#include <indcr.h>
#include <indqr.h>
#include <infpri.h>
C
      PARAMETER ( DIFFR = 1.0D-8 )
      CHARACTER*8 NEWLBL
C
C Label all EXCITLABvectors with positive frequencies
C
      TMPFRQ = ABS(FRQVAL)
C
      IF (NEWLBL.NE.'EXCITLAB') THEN
C
      DO 100 I=1,NLRLBL
         IF ( (NEWLBL .EQ. CRLBL(I) ) .AND.
     *       ( ABS(FRQVAL- CRFREQ(I)).LE.DIFFR) ) THEN
            INCRLR = I + NEXLBL
            RETURN
         END IF
 100  CONTINUE
      NLRLBL = NLRLBL + 1
      IF (NLRLBL.GT.MXLRCR) THEN
         WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *   ' NUMBER OF SPECIFIED PROPERTIES EXCEED THE MAXIMUM ALLOWED',
     *   ' MXLRCR =',MXLRCR,' NLRLBL = ',NLRLBL
         CALL QUIT(' INCRLR: TOO MANY PROPERTIES SPECIFIED')
      END IF
      CRLBL(NLRLBL) = NEWLBL
      CRFREQ(NLRLBL) = FRQVAL
      ISYMCR(NLRLBL) = ISYM
      INCRLR = NLRLBL
      IF (IPRRSP.GT.40) THEN
         WRITE(LUPRI,'(/A)')' INCRLR: NEW OPERATOR'
         WRITE(LUPRI,'(A,I5,2X,A,A,2X,/A,D12.6,2X,A,I5)')
     *   ' NUMBER:',NLRLBL,' LABEL: ',NEWLBL,' FREQUENCY: ',FRQVAL,
     *   'SYMMETRY',ISYM
      END IF
C
      ELSE
C
      IROOT = 1
      DO 200 I=1,NEXLBL
         IF ( ISYM .EQ. ISEXCR(I) )  THEN
            IF ( ABS(TMPFRQ-EXCIT2(ISYM,IROOT)) .LE. DIFFR ) THEN
               INCRLR = I
               RETURN
            END IF
            IROOT = IROOT + 1
         END IF
 200  CONTINUE
      NEXLBL = NEXLBL + 1
      IF ( NEXLBL.GT.MXEXCR ) THEN
         WRITE(LUPRI,'(/A,/A,I5,A,I5)')
     *   'INCRLR: SPECIFIED OPERATORS EXCEED THE ALLOWED NUMBER'
     *   ,' NEXLBL =',NEXLBL,' MXEXCR =',MXEXCR
          CALL QUIT(' INCRLR: TOO MANY EXCITATION OPERATORS')
      END IF
      JEXCR(NEXLBL) = IROOT
      ISEXCR(NEXLBL) = ISYM
      INCRLR = NEXLBL
      IF (IPRRSP.GT.40) THEN
         WRITE(LUPRI,'(/A)')' INCRLR: NEW EXCITATION '
         WRITE(LUPRI,'(A,I5,2X,A,A,2X,/A,D12.6,2(2X,A,I5))')
     *   ' NUMBER:',NEXLBL,' LABEL: ',NEWLBL,' FREQUENCY: ',TMPFRQ,
     *   'SYMMETRY',ISYM,' ROOT NUMBER',IROOT
      END IF
C
      END IF
      RETURN
      END
      INTEGER FUNCTION INCR2(NEWLB1,NEWLB2,FRQ1,FRQ2,ISYM1,ISYM2)
C
C DETERMINE IF A POINTER TO TWO-INDEX LINEAR RESPONSE EQUATION
C CHARACTERIZED BY TWO NEWLBL AND TWO FRQVAL
C IF IN LIST INCR2 GIVES ITS NUMBER IN THE LIST
C IF THE POINTER DOES NOT APPEAR IT IS ADDED TO
C LIST AND INCR2 GIVES THE ADDED NUMBER
C
C BASED ON INCRLR
#include <implicit.h>
C
#include <priunit.h>
#include <infrsp.h>
#include "rspprp.h"
#include <indcr.h>
#include <indqr.h>
#include <infpri.h>
#include <inforb.h>
C
      PARAMETER ( DIFFR = 1.0D-8 )
      CHARACTER*8 NEWLB1,NEWLB2
C
C CHECK LIST IF ALREADY IN, IF SO RETURN
C NOTE THAT Nxy(w1,w2) = Nyx(w2,w1)
C
      DO 100 I=1,NLRLB2
         IF ( (NEWLB1 .EQ. CRLB2(I,1) ) .AND.
     *       ( ABS(FRQ1- CRFRQ2(I,1)).LE.DIFFR) ) THEN
            IF ( (NEWLB2 .EQ. CRLB2(I,2) ) .AND.
     *          ( ABS(FRQ2- CRFRQ2(I,2)).LE.DIFFR) ) THEN
               INCR2 = I
               RETURN
            END IF
         END IF
 100  CONTINUE
C
      DO 200 I=1,NLRLB2
         IF ( (NEWLB2 .EQ. CRLB2(I,1) ) .AND.
     *       ( ABS(FRQ2- CRFRQ2(I,1)).LE.DIFFR) ) THEN
            IF ( (NEWLB1 .EQ. CRLB2(I,2) ) .AND.
     *          ( ABS(FRQ1- CRFRQ2(I,2)).LE.DIFFR) ) THEN
               INCR2 = I
               RETURN
            END IF
         END IF
 200  CONTINUE
C
C ADD ITEM TO LIST
C
      NLRLB2 = NLRLB2 + 1
      IF (NLRLB2.GT.MXLRCR) THEN
         WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *   ' NUMBER OF SPECIFIED PROPERTIES EXCEED THE MAXIMUM ALLOWED',
     *   ' MXLRCR =',MXLRCR,' NLRLB2 = ',NLRLB2
         CALL QUIT(' INCR2: TOO MANY PROPERTIES SPECIFIED')
      END IF
      CRLB2(NLRLB2,1) = NEWLB1
      CRLB2(NLRLB2,2) = NEWLB2
      CRFRQ2(NLRLB2,1) = FRQ1
      CRFRQ2(NLRLB2,2) = FRQ2
      ISMCR2(NLRLB2,1) = ISYM1
      ISMCR2(NLRLB2,2) = ISYM2
      ISMCR2(NLRLB2,3) = MULD2H(ISYM1,ISYM2)
C
C PRINT ITEM
C
      IF (IPRRSP.GT.40) THEN
         WRITE(LUPRI,'(/A,2X,I5)')
     *   ' INCR2: NEW TWO-INDEX VECTOR OF SYMMETRY',ISMCR2(NLRLB2,3)
         WRITE(LUPRI,'(A,I5,2X,A,A,2X,A,A)')' NUMBER:',
     *   NLRLB2,' LABEL1: ',NEWLB1,' LABEL2: ',NEWLB2
         WRITE(LUPRI,'(A,D12.6,2X,A,D12.6)')
     *   ' FREQ1: ',FRQ1,' FREQ2: ',FRQ2
      END IF
C
C RETURN WITH POINTER VALUE
C
      INCR2 = NLRLB2
C
      RETURN
      END
      SUBROUTINE CRCALC(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     *                  XINDX,WRK,LWRK)
C
#include <implicit.h>
#include <iratdef.h>
C
C PURPOSE:
C DRIVER ROUTINE FOR CUBIC RESPONSE
C
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
C
#include <dummy.h>
#include <maxorb.h>
#include <priunit.h>
#include <infvar.h>
#include "rspprp.h"
#include <qrinf.h>
#include <infpp.h>
#include <infpri.h>
#include <inforb.h>
#include <infsmo.h>
#include <infrsp.h>
#include <infcr.h>
#include <inftap.h>
#include <indcr.h>
#include <infdim.h>
#include <inftmo.h>
#include <inftpa.h>
C
      CALL QENTER('CRCALC')
      IF (SOPPA) THEN
         WRITE (LUPRI,*)
     &   'CRCALC ERROR: cubic response not implemented for SOPPA'
         CALL QUIT(
     &   'CRCALC ERROR: cubic response not implemented for SOPPA')
      END IF
C
C Determine solution vectors; first one-index and then two-index
C
      IF (IPRRSP.GT.0) CALL TIMER('START ',TIMSTR,TIMEND)
      KMJWOP = 1
C
C     MJWOP (previously resided in qrinf.h) only needed for MCSCF calculations.
C     
      KTMP = KMJWOP + (16*MAXWOP + 1)/IRAT
      LWRK1  = LWRK - KTMP + 1
C
      CALL CRVEC1(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     *           XINDX,WRK(KMJWOP),WRK(KTMP),LWRK1)
      IF (IPRRSP.GT.0) CALL TIMER('CRVEC1',TIMSTR,TIMEND)
C
      IF (IPRRSP.GT.0) CALL TIMER('START ',TIMSTR,TIMEND)
C
      CALL CRVEC2(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     *           XINDX,WRK(KMJWOP),WRK(KTMP),LWRK1)
      IF (IPRRSP.GT.0) CALL TIMER('CRVEC2',TIMSTR,TIMEND)
C
C Delete files to save disk space
C Keep LURSP and LUXYVE for later use
C
      LURSP1 = -1
      CALL GPOPEN(LURSP1,'RSPRST.E2C','UNKNOWN',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      CALL GPCLOSE(LURSP1,'DELETE')
      CALL GPCLOSE(LURSP3,'DELETE')
      CALL GPCLOSE(LURSP4,'DELETE')
      CALL GPCLOSE(LURSP5,'DELETE')
C
C Layout some workspace
C
      KVECA  = KTMP
      KVECB  = KVECA  + MZYVMX
      KVECC  = KVECB  + MZYVMX
      KVECD  = KVECC  + MZYVMX
      KVECBC = KVECD  + MZYVMX
      KVECBD = KVECBC + MZYVMX
      KVECCD = KVECBD + MZYVMX
      KRESVE = KVECCD + MZYVMX
      KFREE  = KRESVE + MZYVMX
      LFREE  = LWRK   - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('CRCALC',KFREE-1,LWRK)
C
C Calculate the second hyperpolarizability
C
      IF (CRCAL) CALL
     *   CRHYP(WRK(KVECA),WRK(KVECB),WRK(KVECC),WRK(KVECD),
     *         WRK(KVECBC),WRK(KVECBD),WRK(KVECCD),WRK(KRESVE),
     *         CMO,UDV,PV,FOCK,FC,FV,
     *         XINDX,WRK(KMJWOP),WRK(KFREE),LFREE)
      IF (TOMOM) CALL
     *   CRTMO(WRK(KVECA),WRK(KVECB),WRK(KVECC),WRK(KVECD),
     *         WRK(KVECBC),WRK(KVECBD),WRK(KVECCD),WRK(KRESVE),
     *         CMO,UDV,PV,FOCK,FC,FV,
     *         XINDX,WRK(KMJWOP),WRK(KFREE),LFREE)
      IF (TPAMP) CALL
     *   CRTPA(WRK(KVECA),WRK(KVECB),WRK(KVECC),WRK(KVECD),
     *         WRK(KVECBC),WRK(KVECBD),WRK(KVECCD),WRK(KRESVE),
     *         CMO,UDV,PV,FOCK,FC,FV,
     *         XINDX,WRK(KMJWOP),WRK(KFREE),LFREE)
C
      CALL GPCLOSE(LUXYVE,'DELETE')
C
      CALL QEXIT('CRCALC')
      RETURN
      END

      SUBROUTINE CRVEC1(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     *                 XINDX,MJWOP,WRK,LWRK)
C
C PURPOSE:
C  CREATE THE SOLUTION VECTORS THAT ARE REQUIRED TO DO
C  CUBIC RESPONSE
C  SOLUTION OF LINEAR EQUATIONS ARE DETERMINED IN CRLRV1.
C
#include <implicit.h>
#include <iratdef.h>
C
      CHARACTER*8 BLANK
      PARAMETER (BLANK = '        ', D1 = 1.0D0, D0 = 0.0D0)
      LOGICAL FOUND, CONV
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
C
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infvar.h>
#include <infdim.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include "rspprp.h"
#include <indqr.h>
#include <qrinf.h>
      DIMENSION MJWOP(2,MAXWOP,8)
#include <infsmo.h>
#include <infhyp.h>
#include <infpp.h>
#include <inftap.h>
#include <infpri.h>
#include <infspi.h>
#include <infcr.h>
#include <indcr.h>
#include <inftmo.h>
#include <inftpa.h>
C
C INITIALIZE VARIABLES
C
      TRPLET = .FALSE.
      NSOLVC = NLRLBL + NEXLBL
      NVARMA = NCONMA + NWOPMA
      MZYVMX = 2*NVARMA
C
C CHECK NUMBER OF LINEAR RESPONSE EQUATION TO SOLVE
C
      IF (NSOLVC .LE. 0) THEN
         WRITE(LUPRI,'(/A)')
     *     ' CRVEC1: number of solution vectors incorrect'
         WRITE(LUPRI,'(/A)') 'NLRLBL :'
         WRITE(LUPRI,'(/A,4I5)') NLRLBL
         CALL QUIT(' CRVEC1: number of solution vectors incorrect')
      END IF
C
C Solve the eigenvalue problem ( E[2] - w S[2] ) N^F = 0
C
      IF (NEXLBL.GT.0) THEN
C
C Sort eigenvalue eqautions according to symmetry
C
         CALL CREXSR
C
C Solve eigenvalue equation
C
         WRITE(LUPRI,'(//A/)')
     *' <<<<< Solving eigenvalue equations for cubic response >>>>>'
         DO 100 KSYMOP = 1,NSYM
C
C Set up vectors defining the number of eigenvectors that need
C to be solved
C
            NPPCNV(KSYMOP) = NEXCR(KSYMOP)
            NPPSIM(KSYMOP) = NEXCR(KSYMOP)
            NPPSTV(KSYMOP) = NEXCR(KSYMOP)
C
C Define variables that depend on symmetry
C
            CALL RSPVAR(UDV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
            CALL SETZY(MJWOP)
C
C Skip symmetry if no operators
C
            IF ( NEXCR(KSYMOP).EQ.0 ) GOTO 100
C
C     Then we check if there appears to be correct number of eigenvectors
C     of this symmetry on file
C
            CALL REARSP(LURSP,KLEN,WRK,'EXCITLAB',BLANK,D0,D1,KSYMOP,
     &                  NEXCR(KSYMOP),THCPP,FOUND,CONV,ANTSYM)
            IF (FOUND .AND. CONV) THEN
               WRITE (LUPRI,'(A,I3,A)') 'We seem to have all '//
     &              'excitation vectors of symmetry',KSYMOP,' on file'
               GOTO 100
            END IF
C
            IF ( KZVAR.EQ.0 ) THEN
               NWARN = NWARN + 1
               WRITE(LUPRI,'(/2A)')' ****WARNING******',
     *            ' NUMBER OF VARIABLES IN THIS SYMMETRY IS ZERO'
               GO TO 100
            END IF
C
            CALL CRPPVE(CMO,UDV,PV,FC,FV,FCAC,H2AC,XINDX,
     *                  WRK,LWRK)
 100     CONTINUE
C
      END IF
C
C SET UP POINTER FOR THE LINEAR RESPONSE EQUATIONS THAT MUST
C BE SOLVED FOR CUBIC RESPONSE. THE PREVIOUS SET UP GAVE
C THE RIGHT NUMBER OF VECTORS BUT THE FREQUENCIES WERE NOT
C CORRECT FOR FREQUENCIES CONTAINING EXCITATION ENERGIES.
C THE CORRECT EXCITATION ENERGIES ARE STORED IN THE
C VARIABLE EXCIT2
C
      IF (TOMOM) THEN
         NLRLBL = 0
         CALL TM1IND
      ELSE IF (TPAMP) THEN
         NLRLBL = 0
         CALL TP1IND
      END IF
C
C Solve the linear equation ( E[2] - w S[2] ) N^X = X[1]
C
      IF ( NLRLBL.GT.0 ) THEN
C
C SORT SINGLET LINEAR RESPONSE EQUATIONS ACCORDING TO SYMMETRY
C
         CALL CRLRSR
C
C SOLVE SINGLET LINEAR RESPONSE EQUATIONS
C
         WRITE(LUPRI,'(//A/)')
     *   ' <<<<< Solving set of linear response equations'//
     *   ' for cubic response >>>>>'
         KEXSIM = 1
         KEXCNV = 1
         DO 200 KSYMOP = 1,NSYM
C
C Define variables that depend on symmetry
C
            CALL RSPVAR(UDV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
            CALL SETZY(MJWOP)
C
C Skip this symmetry if no operators
C
            IF ( NLRCR(KSYMOP).EQ.0 ) GO TO 200
C
            IF ( KZVAR.EQ.0 ) THEN
               NWARN = NWARN + 1
               WRITE(LUPRI,'(/2A)')' ****WARNING******',
     *            ' NUMBER OF VARIABLES IN THIS SYMMETRY IS ZERO'
               GO TO 200
            END IF
            CALL CRLRV1(CMO,UDV,PV,FC,FV,FCAC,H2AC,XINDX,
     *                  WRK,LWRK)
            CALL FLSHFO(LUPRI)
  200    CONTINUE
      END IF
C
      RETURN
      END
      SUBROUTINE CRLRSR
C
C SORT THE LINEAR RESPONSE EQUATIONS AFTER SYMMETRY
C THE NUMBER OF LINEAR EQUATIONS OF A GIVEN SYMMETRY ARE
C STORED IN NLRCR WITH OFSET IN ILRCR
C
#include <implicit.h>
C
      CHARACTER*8 LABEL
C
#include <priunit.h>
#include <inforb.h>
#include "rspprp.h"
#include <indqr.h>
#include <indcr.h>
#include <infpri.h>
#include <infrsp.h>
C
      IF (IPRRSP.GT.40) THEN
         WRITE(LUPRI,'(/A,I5,2X,A)')
     *   ' LIST OF NLRLBL= ',NLRLBL,
     *   ' LINEAR EQUATIONS FOR CUBIC RESPONSE'
         WRITE(LUPRI,'(A)')' BEFORE SORTING'
         WRITE(LUPRI,'(/A)')'   I  CRLBL(I)  CRFREQ(I)  ISYMCR(I) '
         DO 90 I = 1,NLRLBL
            WRITE(LUPRI,'(I5,5X,A,D12.6,I12)')
     *      I,CRLBL(I),CRFREQ(I),ISYMCR(I)
  90     CONTINUE
      END IF
      ISMMAX = 1
      DO ISYM = 1,NSYM
         DO ISRT = 1,NLRLBL
            IF (ISRT.GE.ISMMAX .AND. ISYMCR(ISRT).EQ.ISYM) THEN
               ICR   = ISYMCR(ISRT)
               FRVAL  = CRFREQ(ISRT)
               LABEL = CRLBL(ISRT)
               ISYMCR(ISRT) = ISYMCR(ISMMAX)
               CRFREQ(ISRT) = CRFREQ(ISMMAX)
               CRLBL(ISRT)  = CRLBL(ISMMAX)
               ISYMCR(ISMMAX) = ICR
               CRFREQ(ISMMAX) = FRVAL
               CRLBL(ISMMAX)  = LABEL
               ISMMAX = ISMMAX + 1
            END IF
         END DO
      END DO
      ITOTOP = 0
      DO 300 ISYM = 1,NSYM
         INUMOP = 0
         DO 400 I = 1,NLRLBL
            IF(ISYMCR(I).EQ.ISYM) THEN
               INUMOP = INUMOP + 1
            END IF
 400     CONTINUE
         ILRCR(ISYM) = ITOTOP
         ITOTOP = ITOTOP + INUMOP
         NLRCR(ISYM) = INUMOP
 300  CONTINUE
      IF (IPRRSP.GT.40) THEN
         WRITE(LUPRI,'(/A,I5,2X,A)')
     *   ' LIST OF NLRLBL= ',NLRLBL,
     *   ' LINEAR EQUATIONS FOR CUBIC RESPONSE'
         WRITE(LUPRI,'(A)')'AFTER SORTING'
         WRITE(LUPRI,'(/A)')'   I  CRLBL(I)  CRFREQ(I)  ISYMCR(I) '
         DO 301 I = 1,NLRLBL
            WRITE(LUPRI,'(I5,5X,A,D12.6,I12)')
     *      I,CRLBL(I),CRFREQ(I),ISYMCR(I)
 301     CONTINUE
      END IF
      IF (IPRRSP.GT.35) THEN
         WRITE(LUPRI,'(/A)')
     *   ' NUMBER OF LINEAR RESPONSE EQUATIONS IN VARIOUS SYMMETRIES'
         WRITE(LUPRI,'(/A)')' NLRCR(ISYM),ISYM=1,NSYM '
         WRITE(LUPRI,'(8I5)')( NLRCR(ISYM),ISYM=1,NSYM )
         WRITE(LUPRI,'(/A)')
     *   ' WITH OFFSET'
         WRITE(LUPRI,'(/A)')' ILRCR(ISYM),ISYM=1,NSYM '
         WRITE(LUPRI,'(8I5)')( ILRCR(ISYM),ISYM=1,NSYM )
      END IF
      RETURN
      END
