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
C ============================================================================
C     newabs.F: Main author Joanna Kauczor
C     This implementation is published in:
C
C     Reference C:
C        J. Kauczor, P. Jorgensen, and P. Norman, 
C        "On the efficiency of algorithms for solving Hartree-Fock and Kohn-Sham
C        response equations",
C        J. Chem. Theory Comput. 7 (2011) 1610.
C ============================================================================
C
      SUBROUTINE ABSORP_INPUT(WORD)
C
C     Purpose:
C     Read in user settings for imaginary polarizabilities describing
C     absorption in the optical processes.
C
C Used from
C  codata.h : XTKAYS
C  gnrinf.h : THR_REDFAC
C  symmet.h : ISYMAX
C  infvar.h : NCONF
C  inftab.h : LUPROP
#include "implicit.h"
#include "priunit.h"
#include "codata.h" 
#include "gnrinf.h"
#include "maxorb.h"
C MXCORB for symmet
#include "maxaqn.h"
C MXQN for symmet
#include "mxcent.h"
C MXCENT for symmet
#include "symmet.h"
C ISYMAX
#include "infvar.h"
C NCONF
#include "abslrs.h"
#include "absorp.h"
#include "inforb.h"
#include "inftap.h"
#include "infinp.h"
#include "infrsp.h"
C
      LOGICAL NEWDEF, ALLCOMP, NXTLAB
      PARAMETER ( NTABLE = 32 )
      PARAMETER ( THRFRQ=1.0D-14, THD=1.0D-6 )
      CHARACTER PROMPT*1, WORD*7, WORD1*7, TABLE(NTABLE)*7
      CHARACTER*8 DIPLEN(3),ANGMOM(3),THETA(6),RTNLBL(2),LABFND
      INTEGER NUMBER_ORBS(8), NUMBER_EXORBS(8)
C
      logical, external :: fun_is_ready_for_qr
c
      DATA DIPLEN/'XDIPLEN','YDIPLEN','ZDIPLEN'/
      DATA ANGMOM/'XANGMOM','YANGMOM','ZANGMOM'/
      DATA THETA/'XXTHETA','XYTHETA','XZTHETA','YYTHETA','YZTHETA',
     &           'ZZTHETA'/
      DATA TABLE /'.ALPHA ','.FREQUE','.THCLR ','.MAXIT ','.MAXRM ',
     &            '.PRINT ','.DAMPIN','.ANALYZ','.FREQ I','.MCD   ',
     &            '.MCHD  ','.NSCD  ','.BETA  ','.BFREQ ','.BFREQI',
     &            '.CFREQ ','.CFREQI','.SHG   ','.XXCOMP','.YYCOMP',
     &            '.ZZCOMP','.IMAG F','.NBATCH','.NOREBD','.OLDCPP',
     &            '.MAX MI','.MAXITO','.EXCITA','.MAX MA','.THCPP ',
     &            '.REDUCE','.BATCH '/
C
      NEWDEF = (WORD .EQ. '*ABSORP')
      ICHANG = 0
C
      ALLCOMP = .TRUE.
      ICOMP=0
      ABSLRS=.TRUE.
      ABSLRS_MCD=.FALSE.
      ABSLRS_MCHD=.FALSE.
      ABSLRS_NSCD=.FALSE.
C
C     Set default values
C
      IF (NEWDEF) THEN
C     Common default value of the damping parameter is set to 
C     be 1000 cm-1 = 4.556333D-3 a.u.
         ABS_DAMP = 1000/XTKAYS
  100    CONTINUE
C
C     Read in input
C
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') GOTO 100
            IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) 
     &                GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
     &                     16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,
     &                     31,32),I
               END DO
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GOTO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' KEYWORD "',WORD,
     *            '" NOT RECOGNIZED IN ABSORP_INPUT.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT(' ILLEGAL KEYWORD IN ABSORP_INPUT ')
C
 1             CONTINUE
               ABSLRS_ALPHA=.TRUE.
               ABSLRS_BETA=.FALSE.
               ABSLRS_GAMMA=.FALSE.
               GOTO 100
 2             CONTINUE
               IF (ABSLRS_ALPHA) THEN
                  READ (LUCMD,*) ABS_NFREQ_ALPHA
                  IF (ABS_NFREQ_ALPHA.LE.NMXFREQ) THEN
                     READ (LUCMD,*) (ABS_FREQ_ALPHA(J),
     &                              J=1,ABS_NFREQ_ALPHA)
                  ELSE
                     WRITE (LUPRI,'(3(/,A,I5),/)')
     &               ' NUMBER OF FREQUENCIES SPECIFIED   :',
     &                                              ABS_NFREQ_ALPHA,
     &               ' IS GREATER THAN THE ALLOWED NUMBER:',NMXFREQ,
     &               ' THE NUMBER IS RESET TO THE MAXIMUM:',NMXFREQ
                     READ (LUCMD,*) (ABS_FREQ_ALPHA(J),J=1,NMXFREQ),
     &                              (FFFF,J=NMXFREQ+1,ABS_NFREQ_ALPHA)
                     ABS_NFREQ_ALPHA = NMXFREQ
                  END IF
                  CALL GPDSRT(ABS_NFREQ_ALPHA,ABS_FREQ_ALPHA,THD)
               ELSE
                  WRITE(LUPRI,1000) WORD
                  READ(LUCMD,*)
                  READ(LUCMD,*)
               END IF
               GOTO 100
 3             CONTINUE
                  READ (LUCMD,*) ABS_THCLR
                  IF (ABS_OLDCPP) THCLR_ABSORP=ABS_THCLR
               GOTO 100
 4             CONTINUE
                  READ (LUCMD,*) ABS_MAXITER
               GOTO 100
 5             CONTINUE
                  READ (LUCMD,*) ABS_MAXRM
               GOTO 100
 6             CONTINUE
                  READ (LUCMD,*) IPRABSLRS
               GOTO 100
 7             CONTINUE
C     User input of damping parameter
                  READ (LUCMD,*) ABS_DAMP
               GOTO 100
 8             CONTINUE
               ABSLRS_ANALYZE=.TRUE.
               GOTO 100
 9             CONTINUE
               ABSLRS_INTERVAL = .TRUE.
c               ABS_REDUCE = .TRUE.
               READ(LUCMD,*) (ABS_FREQ_INTERVAL(I), I=1,3)
               GOTO 100
 10            CONTINUE
               ABSLRS_ALPHA=.FALSE.
               ABSLRS_MCD  =.TRUE.
               GOTO 100
 11            CONTINUE
               ABSLRS_ALPHA=.FALSE.
               ABSLRS_MCHD  =.TRUE.
               GOTO 100
 12            CONTINUE
                ABSLRS_ALPHA=.FALSE.
                ABSLRS_NSCD  =.TRUE.
                GOTO 100
 13            CONTINUE
                ABSLRS_ALPHA=.FALSE.
                ABSLRS_BETA  =.TRUE.
                ABSLRS_GAMMA =.FALSE.
                GOTO 100
 14            CONTINUE
               IF (ABSLRS_BETA .OR. ABSLRS_MCD .OR. ABSLRS_MCHD .OR. 
     &             ABSLRS_NSCD) THEN
                  READ (LUCMD,*) ABS_NFREQ_BETA_B
                  IF (ABS_NFREQ_BETA_B.LE.NMXFREQ) THEN
                     READ (LUCMD,*) (ABS_FREQ_BETA_B(J),
     &               J=1,ABS_NFREQ_BETA_B)
                  ELSE
                     WRITE (LUPRI,'(3(/,A,I5),/)')
     &              ' NUMBER OF FREQUENCIES SPECIFIED   :',
     &                                            ABS_NFREQ_BETA_B,
     &              ' IS GREATER THAN THE ALLOWED NUMBER:',NMXFREQ,
     &              ' THE NUMBER IS RESET TO THE MAXIMUM:',NMXFREQ
                     READ (LUCMD,*) (ABS_FREQ_BETA_B(J),J=1,NMXFREQ),
     &                              (FFFF,J=NMXFREQ+1,ABS_NFREQ_BETA_B)
                     ABS_NFREQ_BETA_B = NMXFREQ
                  END IF
               ELSE
                  WRITE(LUPRI,1010) WORD
                  READ(LUCMD,*)
                  READ(LUCMD,*)
               END IF
               GOTO 100
 15            CONTINUE
               ABSQRS_INTERVAL = .TRUE.
c               ABS_REDUCE = .TRUE.
               READ(LUCMD,*) (ABS_BFREQ_INTERVAL(I), I=1,3)
               GOTO 100
 16            CONTINUE
               IF (ABSLRS_BETA) THEN
                 READ (LUCMD,*) ABS_NFREQ_BETA_C
                 IF (ABS_NFREQ_BETA_C.LE.NMXFREQ) THEN
                   READ (LUCMD,*) (ABS_FREQ_BETA_C(J),
     &                             J=1,ABS_NFREQ_BETA_C)
                 ELSE
                   WRITE (LUPRI,'(3(/,A,I5),/)')
     &          ' NUMBER OF FREQUENCIES SPECIFIED   :',ABS_NFREQ_BETA_C,
     &          ' IS GREATER THAN THE ALLOWED NUMBER:',NMXFREQ,
     &          ' THE NUMBER IS RESET TO THE MAXIMUM:',NMXFREQ
                   READ (LUCMD,*) (ABS_FREQ_BETA_C(J),J=1,NMXFREQ),
     &                              (FFFF,J=NMXFREQ+1,ABS_NFREQ_BETA_C)
                     ABS_NFREQ_BETA_C = NMXFREQ
                  END IF
               ELSE
                  WRITE(LUPRI,1010) WORD
                  READ(LUCMD,*)
                  READ(LUCMD,*)
               END IF
               GOTO 100
 17            CONTINUE
               ABSQRS_CINTERVAL = .TRUE.
c               ABS_REDUCE = .TRUE.
               READ(LUCMD,*) (ABS_CFREQ_INTERVAL(I), I=1,3)
               GOTO 100
 18            CONTINUE
               ABSLRS_ALPHA=.FALSE.
               ABSLRS_BETA =.TRUE.
               ABSLRS_GAMMA=.FALSE.
               ABSLRS_SHG  =.TRUE.
               GOTO 100
 19            CONTINUE      
               ALLCOMP=.FALSE.
               ICOMP=1
               GOTO 100
 20            CONTINUE
               ALLCOMP=.FALSE.
               ICOMP=2
               GOTO 100
 21            CONTINUE
               ALLCOMP=.FALSE.
               ICOMP=3
               GOTO 100
 22            CONTINUE
               ABS_IMFREQ=.TRUE.
               ABS_DAMP=0.0d0
               GOTO 100
 23            CONTINUE
               ABS_BATCH=.TRUE.
               READ(LUCMD,*) ABS_NBATCHMAX
               GOTO 100
 24            CONTINUE
                 ABS_REBD = .FALSE.
               GOTO 100
 25            CONTINUE
                 ABS_OLDCPP = .TRUE.
               GOTO 100
 26            CONTINUE
                 IF (ABS_OLDCPP) THEN
                  READ (LUCMD,*) MAX_MICRO
               ELSE
                  WRITE(LUPRI,1020) WORD
                  READ(LUCMD,*)
               END IF
               GOTO 100
 27            CONTINUE
                 IF (ABS_OLDCPP) THEN
                  READ (LUCMD,*) MAX_ITORB
                  IF (MAX_ITORB .GT. 0) OPTORB = .TRUE.
               ELSE
                  WRITE(LUPRI,1020) WORD
                  READ(LUCMD,*)
               END IF
               GOTO 100
 28             CONTINUE
                 IF (ABS_OLDCPP) THEN
                  READ (LUCMD,*) NEXCITED_STATES
               ELSE
                  WRITE(LUPRI,1020) WORD
                  READ(LUCMD,*)
               END IF
               GOTO 100
 29             CONTINUE
                 IF (ABS_OLDCPP) THEN
                  READ (LUCMD,*) MAX_MACRO
               ELSE
                  WRITE(LUPRI,1020) WORD
                  READ(LUCMD,*)
               END IF
               GOTO 100
 30            CONTINUE
                 IF (ABS_OLDCPP) THEN
                  READ (LUCMD,*) THCPP_ABSORP
                  ABS_THCLR=THCPP_ABSORP
               ELSE
                  WRITE(LUPRI,1020) WORD
                  READ(LUCMD,*)
               END IF
               GOTO 100
 31            CONTINUE
                 IF (ABS_OLDCPP) THEN
                 ABS_REDUCE=.TRUE.
               ELSE
                  WRITE(LUPRI,1020) WORD
               END IF
               GOTO 100
 32            CONTINUE
               READ(LUCMD,*) NFREQ_BATCH
               IF (NFREQ_BATCH.GT.MXFREQ) THEN
                  NFREQ_BATCH = MXFREQ
               END IF
               GOTO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GOTO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' PROMPT "',WORD,
     *            '" NOT RECOGNIZED IN ABSORPTION INPUT.'
               CALL QUIT(' ILLEGAL PROMPT IN ABSORPTION INPUT ')
            END IF
         GOTO 100
      END IF
  300 CONTINUE
      IF (THR_REDFAC .GT. 0.0D0) THEN
         ICHANG = ICHANG + 1
         WRITE (LUPRI,'(3A,1P,D10.2)') '@ INFO ',WORD1,
     &   ' thresholds multiplied with general factor',THR_REDFAC
         ABS_THCLR = ABS_THCLR*THR_REDFAC
         THCPP_ABSORP=ABS_THCLR
      END IF
C
      ABS_MAXRM=MAX(ABS_MAXRM,MAXRM)
C
C     Process user input
C
#ifndef DISABLE_XC_RESPONSE_SANITY_CHECK
!     sanity check for incomplete/untested xc functional implementations
      IF ((ABSLRS_BETA .OR. ABSLRS_GAMMA .OR. ABSLRS_MCD .OR. 
     &     ABSLRS_MCHD .OR. ABSLRS_NSCD) .AND. DODFT) THEN
         IF (.NOT. fun_is_ready_for_qr()) THEN
            WRITE(LUPRI, *) 'ERROR: functional not fully implemented'
            WRITE(LUPRI, *) '       or tested for QR'
            WRITE(LUPRI, *) 'to disable this stop recompile'
            WRITE(LUPRI, *) 'with -DDISABLE_XC_RESPONSE_SANITY_CHECK'
            WRITE(LUPRI, *) 'note that GGAKEY functionals are always'
            WRITE(LUPRI, *) 'stopped although they could be correct'
            CALL QUIT('functional not fully implemented/tested for QR')
         END IF
      END IF
#endif

C
      ABSLRS = ABSLRS_ALPHA .OR. ABSLRS_BETA .OR. ABSLRS_GAMMA .OR.
     &         ABSLRS_MCD .OR. ABSLRS_MCHD .OR. ABSLRS_NSCD
      IF (.NOT.ABSLRS) THEN
         WRITE(LUPRI,'(/A)') ' Absorption input ignored because:'
         WRITE(LUPRI,'(A,L2)')' No process requested:',
     &        ' ABS_ALPHA = ', ABSLRS
      ELSE
C
C     Put operators in lists
C
         CALL IZERO(ABS_NOPER,8)
         CALL IZERO(NOPER,8)
         DO ILABEL=1,3
           IF (ALLCOMP .OR. ICOMP.EQ.ILABEL) THEN
             ISYM = ISYMAX(ILABEL,1)+1
             ABS_NOPER(ISYM) = ABS_NOPER(ISYM) + 1
             ABS_LABOP(ABS_NOPER(ISYM),ISYM) = DIPLEN(ILABEL)
             NOPER(ISYM)=NOPER(ISYM)+1
             LABOP(NOPER(ISYM),ISYM) = DIPLEN(ILABEL)
             IF (ABSLRS_MCD) THEN
                ISYM = ISYMAX(ILABEL,2)+1
                ABS_NOPER(ISYM) = ABS_NOPER(ISYM) + 1
                NOPER(ISYM) = NOPER(ISYM) + 1
                ABS_LABOP(ABS_NOPER(ISYM),ISYM) = ANGMOM(ILABEL)
                LABOP(NOPER(ISYM),ISYM) = ANGMOM(ILABEL)
             ELSE IF (ABSLRS_MCHD) THEN
               ISYM2 = ISYMAX(ILABEL,2)+1
               ABS_NOPER(ISYM2) = ABS_NOPER(ISYM2) + 1
               ABS_LABOP(ABS_NOPER(ISYM2),ISYM2) = ANGMOM(ILABEL)
               NOPER(ISYM2) = NOPER(ISYM2) + 1
               LABOP(NOPER(ISYM2),ISYM2) = ANGMOM(ILABEL)
               DO JLABEL=ILABEL,3
                  ISYM3 = ISYMAX(JLABEL,1)+1
                  ISYM4 = MULD2H(ISYM,ISYM3)
                  ABS_NOPER(ISYM4) = ABS_NOPER(ISYM4) + 1
                  NOPER(ISYM4) = NOPER(ISYM4) + 1
                  IF (ILABEL .EQ.1) THEN
                     KLABEL=JLABEL
                  ELSE
                     KLABEL=JLABEL+ILABEL
                  ENDIF
                  ABS_LABOP(ABS_NOPER(ISYM4),ISYM4) = THETA(KLABEL)
                  LABOP(NOPER(ISYM4),ISYM4) = THETA(KLABEL)
                ENDDO
             END IF
           END IF
         END DO
         IF (ABSLRS_NSCD) THEN
           !stop if we are using symmetry.
           if (NSYM.gt.1) then
               WRITE (LUPRI,*)
     &  'NSCD not yet working symmetry: remove sym and resubmit'
               CALL QUIT('NSCD not fully working with symmetry')
           end if
           CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',
     &                   'UNFORMATTED',IDUMMY,.FALSE.)
 333      CONTINUE
           IF (NXTLAB(LABFND,RTNLBL,LUPROP)) THEN
             IF (LABFND(1:3) .EQ. 'PSO') THEN
               !IND is symmetry of operator, resumed from file. Sonia
               IND = (ICHAR(RTNLBL(1)(1:1))-ICHAR('0'))
               ABS_NOPER(IND) = ABS_NOPER(IND) + 1
               NOPER(IND) = NOPER(IND) + 1
               ABS_LABOP(ABS_NOPER(IND),IND) = LABFND
               LABOP(NOPER(IND),IND) = LABFND
             ENDIF
             GO TO 333
           END IF
           CALL GPCLOSE(LUPROP,'KEEP')
         END IF
C
         IF (ABSLRS_INTERVAL) THEN 
            IF( ABS_FREQ_INTERVAL(1).GT.ABS_FREQ_INTERVAL(2) .OR.
     &         (ABS_FREQ_INTERVAL(2)-ABS_FREQ_INTERVAL(1)).LT.
     &          ABS_FREQ_INTERVAL(3)
     &         .OR. ABS_FREQ_INTERVAL(3).LE.0.0D0 ) THEN
               CALL QUIT('.FREQ_INTERVAL not specify correctly')
            END IF
         END IF
         IF (ABSQRS_INTERVAL) THEN
           IF( ABS_BFREQ_INTERVAL(1).GT.ABS_BFREQ_INTERVAL(2) .OR.
     &        (ABS_BFREQ_INTERVAL(2)-ABS_BFREQ_INTERVAL(1)).LT.
     &         ABS_BFREQ_INTERVAL(3)
     &        .OR. ABS_BFREQ_INTERVAL(3).LE.0.0D0 ) THEN
              CALL QUIT('.BFREQI not specified correctly')
           END IF
           DSMALL=1.0000001
           ABS_NFREQ_BETA_B = 1 +
     &     INT(DSMALL*(ABS_BFREQ_INTERVAL(2)-ABS_BFREQ_INTERVAL(1))/
     &              ABS_BFREQ_INTERVAL(3))
           DO I=1,ABS_NFREQ_BETA_B
             ABS_FREQ_BETA_B(I)=ABS_BFREQ_INTERVAL(1) +
     &              (I-1)*ABS_BFREQ_INTERVAL(3)
           END DO
         ENDIF
         IF (ABSQRS_CINTERVAL) THEN
           IF( ABS_CFREQ_INTERVAL(1).GT.ABS_CFREQ_INTERVAL(2) .OR.
     &        (ABS_CFREQ_INTERVAL(2)-ABS_CFREQ_INTERVAL(1)).LT.
     &         ABS_CFREQ_INTERVAL(3)
     &        .OR. ABS_CFREQ_INTERVAL(3).LE.0.0D0 ) THEN
              CALL QUIT('.CFREQI not specified correctly')
           END IF
           DSMALL=1.0000001
           ABS_NFREQ_BETA_C = 1 +
     &     INT(DSMALL*(ABS_CFREQ_INTERVAL(2)-ABS_CFREQ_INTERVAL(1))/
     &              ABS_CFREQ_INTERVAL(3))
           DO I=1,ABS_NFREQ_BETA_C
             ABS_FREQ_BETA_C(I)=ABS_CFREQ_INTERVAL(1) +
     &              (I-1)*ABS_CFREQ_INTERVAL(3)
           END DO
         ENDIF
C
        IF (ABS_OLDCPP) THEN
          ABSORP=.TRUE.
          ABSLRS=.FALSE.
          WRITE(LUPRI,'(/A)')'.OLDCPP specify' 
          WRITE(LUPRI,'(/A)')'the old CPP solver is used' 
        ELSE 
         CALL AROUND('Variable settings for absorption calculation')
C
         IF (ABSLRS_ALPHA) THEN
            WRITE (LUPRI,'(A,L4)')
     &     ' Linear polarizability calculation requested: 
     &          ABSLRS_ALPHA =',ABSLRS_ALPHA
            IF (.NOT.ABSLRS_INTERVAL) WRITE(LUPRI,'(A,5(4F12.8,/,16X))')
     &        ' at frequencies:',(ABS_FREQ_ALPHA(I),I=1,ABS_NFREQ_ALPHA)
         END IF
C
         IF (ABSLRS_BETA) THEN
            WRITE (LUPRI,'(A,L4)')
     &  ' Nonlinear polarizability calculation requested:ABSLRS_BETA =',
     &           ABSLRS_BETA
            WRITE(LUPRI,'(A,5(4F12.8))')
     &           ' B-FREQ:',(ABS_FREQ_BETA_B(I),I=1,ABS_NFREQ_BETA_B)
            WRITE(LUPRI,'(A,5(4F12.8))')
     &           ' C-FREQ:',(ABS_FREQ_BETA_C(I),I=1,ABS_NFREQ_BETA_C)
         END IF
cC
         IF (ABSLRS_MCD) THEN
            WRITE (LUPRI,'(A,L4)')
     &     ' Magnetic circular dichroism requested : ABSLRS_MCD     = ',
     &           ABSLRS_MCD
            WRITE(LUPRI,'(A,5(4F12.8))')
     &           ' at frequencies:',(ABS_FREQ_BETA_B(I),
     &            I=1,ABS_NFREQ_BETA_B)
         END IF
         IF (ABSLRS_MCHD) THEN
            WRITE (LUPRI,'(A,L4)')
     &     ' Magneto-chiral circular dichroism requested: ABSLRS_MCHD=',
     &           ABSLRS_MCHD
            WRITE(LUPRI,'(A,5(4F12.8))')
     &           ' at frequencies:',(ABS_FREQ_BETA_B(I),
     &            I=1,ABS_NFREQ_BETA_B)
         END IF
cC
         IF (ABSLRS_NSCD) THEN
            WRITE (LUPRI,'(A,L4)')
     &    ' Nuclear-spin circular dichroism/OR requested: ABSLRS_NSCD=',
     &           ABSLRS_NSCD
            WRITE(LUPRI,'(A,5(4F12.8))')
     &           ' at frequencies:',(ABS_FREQ_BETA_B(I),
     &            I=1,ABS_NFREQ_BETA_B)
         END IF
!
         WRITE(LUPRI,'(A,L4)')
     &      ' Absorption calc over frequency interval  : 
     &           ABSLRS_INTERVAL =', ABSLRS_INTERVAL
         IF (ABSLRS_INTERVAL) THEN
c           IF (NCONF.GT.1 ) THEN
c             CALL QUIT('Error occured! Not implemented!')
c           ELSE
             WRITE(LUPRI,'(A,I4)')
     &      ' Number of frequencies              : ABS_NFREQ_INTERVAL='
     &       ,ABS_NFREQ_INTERVAL
           ENDIF
            WRITE(LUPRI,'(3(A,F8.5))')
     &      ' Start:',ABS_FREQ_INTERVAL(1),'Stop:',ABS_FREQ_INTERVAL(2),
     &      '   Step:',ABS_FREQ_INTERVAL(3)
c         END IF
C
         WRITE(LUPRI,'(A,F9.6,F8.1)')
     &      ' Damping parameter (a.u. and cm-1)        : DAMPING      ='
     &      ,ABS_DAMP,ABS_DAMP*XTKAYS
         WRITE(LUPRI,'(A,1P,D8.1)')
     &      ' Threshold of convergence in LR solver    : THCLR        ='
     &      ,ABS_THCLR
         WRITE(LUPRI,'(A,I4)')
     &      ' Maximum iterations in complex LR solver  : MAXIT        ='
     &       ,ABS_MAXITER
         WRITE(LUPRI,'(A,I4)')
     &      ' Maximum size of reduced space            : MAXRM        ='
     &       ,ABS_MAXRM
         WRITE(LUPRI,'(A,I4)')
     &      ' Print level in absorption modules        : PRINT        ='
     &       ,IPRABSLRS
      END IF
      ENDIF
      CALL ABSORP_INTERFACE
C
C Messages
C
 1000 FORMAT(/,
     &  ' Keyword:',A8,' ignored since an ALPHA calc not specified.',
     &  /,' If asked for, the .ALPHA keyword should be placed first',
     &  /,' of *ABSORPTION input keys.',
     &  /,' The program will continue.')

 1010 FORMAT(/,
     &  ' Keyword:',A8,' ignored since a BETA calc not specified.',
     &  /,' If asked for, the .BETA keyword should be placed first',
     &  /,' of *ABSORPTION input keys.',
     &  /,' The program will continue.')


 1020 FORMAT(/,
     &  ' Keyword:',A8,' ignored,calc with old cpp not specified.',
     &  /,' If asked for, the .OLDCPP keyword should be placed first',
     &  /,' of *ABSORPTION input keys.',
     &  /,' The program will continue.')
c
C
C     End of ABSORP_INPUT
C
      RETURN
      END
      SUBROUTINE ABSLRSCALC(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     &                  XINDX,WRK,LWRK)
C
C PURPOSE:
C     Driver routine for the computation of response properties including
C     absorption
C
C Used from:
C   infvar.h : MAXWOP
C   inforb.h : NSYM
C   iradef.h : IRAT
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "abslrs.h"
#include "infvar.h"
#include "inforb.h"
C
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(LWRK)
      LOGICAL EX
C
      KFREE = 1
      LFREE = LWRK
C    
      LUABSVECS = -1
      CALL GPINQ('ABSVECS','EXIST',EX)
      IF (.NOT.EX) THEN
        CALL GPOPEN(LUABSVECS,'ABSVECS','NEW',' ',' ',IDUMMY,.FALSE.)
        WRITE(LUABSVECS) 'EOFLABEL'
      ELSE
        CALL GPOPEN(LUABSVECS,'ABSVECS','OLD',' ',' ',IDUMMY,.FALSE.)
      END IF
C
      CALL MEMGET2('REAL','MJWOP',KMJWOP,(2*8*MAXWOP + 1)/IRAT,WRK,
     &     KFREE,LFREE)
C
      IF (ABSLRS_BETA .OR. ABSLRS_MCD .OR. ABSLRS_MCHD .OR. ABSLRS_NSCD)
     &     THEN
         CALL ABS_BETA_SETUP
      END IF

C
C     Allocate memory for results of linear absorption calc.
C
      DSMALL=1.0000001
      IF (ABSLRS_INTERVAL) THEN
        ABS_NFREQ_INTERVAL = 1 +
     &  INT(DSMALL*(ABS_FREQ_INTERVAL(2)-ABS_FREQ_INTERVAL(1))/
     &              ABS_FREQ_INTERVAL(3))
      ELSE
         ABS_NFREQ_INTERVAL = ABS_NFREQ_ALPHA
      END IF
C
      CALL MEMGET('REAL',KRES,2*ABS_NFREQ_INTERVAL*3*3*3,WRK,
     &                   KFREE,LFREE)
      CALL DZERO(WRK(KRES),2*ABS_NFREQ_INTERVAL*3*3*3)
C
C     Determine linear response vectors
C
      CALL AROUND('Solving Linear Response Equations')
C     default: Use the new SCF solver (see Reference C) 
      CALL ABS_LRS(NSYM,CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     &        XINDX,WRK(KMJWOP),WRK(KRES),WRK(KFREE),LFREE)
C
      IF (ABSLRS_BETA .OR. ABSLRS_MCD .OR. ABSLRS_MCHD .OR. ABSLRS_NSCD)
     &     THEN
         CALL AROUND('Evaluate Quadratic Response Functions')
         CALL ABSQRF(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     &               XINDX,WRK(KMJWOP),WRK(KFREE),LFREE)
      END IF
C
C     Print-out a summary of results
C
      CALL ABSLRSRESULT(WRK(KMJWOP),WRK(KRES),WRK(KFREE),LFREE)
c
      CALL GPCLOSE(LUABSVECS,'KEEP')
C
      RETURN
      END

      SUBROUTINE ABS_LRS(NSYM,CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     &                   XINDX,MJWOP,RESLRF,WRK,LWRK)
C
C PURPOSE:
C     Solve the linear response equations and store response vectors on
C     disk.
C
#include "implicit.h"
#include "priunit.h"
#include "abslrs.h"
c
      LOGICAL CONVERGED
      CHARACTER LABEL*8
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FV(*),FCAC(*),H2AC(*)
c      DIMENSION XINDX(*),MJWOP(2,MAXWOP,8)
      DIMENSION XINDX(*),MJWOP(*)
      DIMENSION RESLRF(2,ABS_NFREQ_INTERVAL,3,3,3),WRK(LWRK)
C     local
      REAL*8  DUMMY, thd
      INTEGER IDUMMY
      PARAMETER ( DUMMY = 1.0D20 , IDUMMY = - 9999999 , thd=1.0D-6  )
      logical negf
C
C
      IF (IPRABSLRS.GT.0) CALL TIMER('START ',TIMSTR,TIMEND)
C 
      negf=.false.
      DO ISYM=1,NSYM
         IF (ABS_NOPER(ISYM).GT.0) THEN
C
            LFREE=LWRK
            KFREE=1
C
            CALL ABSINTER(ISYM,UDV,FOCK,FC,FV,FCAC,H2AC,XINDX,MJWOP,
     &                    KZVAR,WRK(KFREE),LFREE)
C
            KGD    = 1
            kfreqs = KGD + 2*KZVAR
            KXSOL  = kfreqs + abs_nfreq_interval
            KFREE  = KXSOL + 4*ABS_NFREQ_INTERVAL*KZVAR
            LFREE  = LWRK - KFREE 
C
            ABS_KLRED(1)=0
            ABS_KLRED(2)=0
C
            LUSB = -1
            LUAB = -1
            LUSS = -1
            LUAS = -1
C
            CALL GPOPEN(LUSB,'ABS_SB','NEW',' ',' ',IDUMMY,.FALSE.)
            CALL GPOPEN(LUAB,'ABS_AB','NEW',' ',' ',IDUMMY,.FALSE.)
            CALL GPOPEN(LUSS,'ABS_SS','NEW',' ',' ',IDUMMY,.FALSE.)
            CALL GPOPEN(LUAS,'ABS_AS','NEW',' ',' ',IDUMMY,.FALSE.)
C
            LUE1RED = -1
            LUE2RED = -1
            LUSRED  = -1
            CALL GPOPEN(LUE1RED,'ABS_E1RED','NEW',
     &         ' ',' ',IDUMMY,.FALSE.)
            CALL GPOPEN(LUE2RED,'ABS_E2RED','NEW',
     &         ' ',' ',IDUMMY,.FALSE.)
            CALL GPOPEN(LUSRED ,'ABS_SRED' ,'NEW',
     &         ' ',' ',IDUMMY,.FALSE.)
C
            DO IOPER=1,ABS_NOPER(ISYM)
               LABEL=ABS_LABOP(IOPER,ISYM)
               IF (ABSLRS_INTERVAL) THEN
                 DO I=1,ABS_NFREQ_INTERVAL
                   ABS_FREQ_ALPHA(I)=ABS_FREQ_INTERVAL(1) + 
     &              (I-1)*ABS_FREQ_INTERVAL(3)
                 END DO
               END IF
               IF (IPRABSLRS.GE.0) THEN
                  WRITE(LUPRI,'(/2A,/3A,I2,A/,A,5(4F12.8,/,11X))') 
     &                 ' ABSLRS -- Solving linear response equations',
     &                 '( E[2] - {w+iW}*S[2] ) N = B[1]  ',
     &                 ' ABSLRS -- for operator ', LABEL,
     &                 ' of symmetry',ISYM,' at frequencies:',
     &                 ' ABSLRS --',(ABS_FREQ_ALPHA(I),
     &                               I=1,ABS_NFREQ_INTERVAL)
               END IF
C
               CALL GETGPV(LABEL,FC,FV,CMO,UDV,PV,XINDX,ANTSYM,
     &              WRK(KFREE),LFREE)
               CALL DCOPY(2*KZVAR,WRK(KFREE),1,WRK(KGD),1)
               GDNRM=DNRM2(2*KZVAR,WRK(KGD),1)
               RES0=0.0d0
               IF (GDNRM .LE. 1.0d-6) THEN
                 WRK(KXSOL:4*KZVAR-1)=0.0d0
                 DO K=1,ABS_NFREQ_INTERVAL
                   CALL WRITE_XVEC(LUABSVECS,4*KZVAR,WRK(KXSOL),LABEL,
     &                  ABS_FREQ_ALPHA(K),RES0)
                   ENDDO
                 WRITE(LUPRI,*) 'Vectors equal to zero'
               ELSE
                 nfreqs=abs_nfreq_interval
                 do i=1,abs_nfreq_interval
                   if (abs_freq_alpha(i) .lt. 0.d0) negf=.true.
                   wrk(kfreqs-1+i)=abs(abs_freq_alpha(i))
                 enddo
                 if (negf) then
                   call gpdsrt(nfreqs,wrk(kfreqs),thd)
                 endif
                 
                 CALL ABS_CTL(LABEL,KZVAR,WRK(KGD),WRK(KXSOL),
     &                        wrk(kfreqs),nfreqs,
     &                        MJWOP,CMO,UDV,FC,FCAC,FV,PV,XINDX,RESLRF,
     &                        WRK(KFREE),LFREE)
               ENDIF
C
               CALL GET_LRF(LUABSVECS,LABEL,ABS_NFREQ_INTERVAL,
     &                     ABS_FREQ_ALPHA,FC,FV,CMO,UDV,PV,XINDX,
     &                     RESLRF,KZVAR,WRK(KFREE),LFREE)
C
            END DO
C
            CALL GPCLOSE(LUSB,'DELETE')
            CALL GPCLOSE(LUAB,'DELETE')
            CALL GPCLOSE(LUSS,'DELETE')
            CALL GPCLOSE(LUAS,'DELETE')
C
            CALL GPCLOSE(LUE1RED,'DELETE')
            CALL GPCLOSE(LUE2RED,'DELETE')
            CALL GPCLOSE(LUSRED, 'DELETE')
C
         END IF
      END DO
      IF (IPRABSLRS.GT.0) CALL TIMER('ABSVEC',TIMSTR,TIMEND)
C     
      RETURN
      END
      SUBROUTINE ABSINTER(ISYM,UDV,FOCK,FC,FV,FCAC,H2AC,XINDX,MJWOP,
     &                    ABS_KZVAR,WRK,LWRK)
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
C Used from:
C   infvar.h : MAXWOP
C   inforb.h : NSYM
#include "abslrs.h"
#include "absorp.h"
!MAXWOP
#include "infvar.h"
!NSYM
#include "inforb.h"
c
#include "infrsp.h"
#include "wrkrsp.h"
#include "infdim.h"
#include "qrinf.h"
#include "maxorb.h"
#include "maxash.h"
#include "infind.h"
#include "infpri.h"
c
      DIMENSION UDV(*),FOCK(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),MJWOP(2,MAXWOP,8),WRK(LWRK)
      INTEGER ABS_KZVAR,ISYM
C     local
      INTEGER KFREE,LFREE
c     
      KFREE=1
      LFREE=LWRK
C
      KSYMOP      = ISYM
      ABS_GRADSYM = ISYM
      ABS_NSYM    = NSYM
      DO I = 1,NSYM
        ABS_NASH(I) = NASH(I)
        ABS_NISH(I) = NISH(I)
        ABS_NORB(I) = NORB(I)
        DO J = 1,NSYM
          ABS_MULD2H(I,J) = MULD2H(I,J) 
        ENDDO
      ENDDO
      LUABSPRI=LUPRI
      CALL RSPVAR(UDV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK(KFREE),LFREE)
      CALL SETZY(MJWOP)
      ABS_KZVAR=KZVAR
C
      DAMPING=ABS_DAMP
      THCLR_ABSORP=ABS_THCLR
c 
      MAXRM=ABS_MAXRM
      NFREQ_INTERVAL=ABS_NFREQ_INTERVAL
      IPRABS=IPRABSLRS
C
      RETURN
      END
      SUBROUTINE ABS2INTER(MJWOP,ABS_KZVAR)
C
C     Used from:
C       infvar.h: MAXWOP
C       inforb.h: NSYM
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "abslrs.h"
#include "absorp.h"
#include "infvar.h"
#include "inforb.h"
c
#include "infrsp.h"
#include "wrkrsp.h"
#include "infdim.h"
#include "qrinf.h"
#include "maxorb.h"
#include "maxash.h"
#include "infind.h"
#include "infpri.h"
c
      DIMENSION MJWOP(2,MAXWOP,8)
      INTEGER ABS_KZVAR,ISYM
c     
      ABS_GRADSYM = KSYMOP
      ABS_NSYM    = NSYM
      DO I = 1,NSYM
        ABS_NASH(I) = NASH(I)
        ABS_NISH(I) = NISH(I)
        ABS_NORB(I) = NORB(I)
        DO J = 1,NSYM
          ABS_MULD2H(I,J) = MULD2H(I,J)
        ENDDO
      ENDDO
      LUABSPRI=LUPRI
      CALL SETZY(MJWOP)
      ABS_KZVAR=KZVAR
C
      ABS_DAMP=DAMPING
      ABS_THCLR=THCLR_ABSORP
c 
      ABS_MAXRM=MAXRM
      ABS_NFREQ_INTERVAL=NFREQ_INTERVAL
      ABS_NFREQ_ALPHA=NFREQ_ALPHA
      IPRABSLRS=IPRABS
      DO I=1,ABS_NFREQ_INTERVAL
        ABS_FREQ_ALPHA(I)=FREQ_ALPHA(I)
      ENDDO
C
      RETURN
      END
      SUBROUTINE ABS_LINE2(KZVAR,XTMP,ZYMB,ISYMDM,ISYMB,NNMAX,
     &                     CMO,FC,FCAC,FV,UDV,XINDX,MJWOP,WRK,LWRK)
C
#ifdef ENABLE_XCINT_RESPONSE
      use dalton_xcint_interface, only: xcint_is_active,
     &                                  xcint_integrate_rks_lr
#endif
c
C PURPOSE:
C
C
c 
#include "implicit.h"
#include "abslrs.h"
#include "thrzer.h"
#include "priunit.h"
C
C Used from:
C   inforb.h: NORB,NORBT,NASHT,NISHT,N2BASX,NBAST,NBAS,NSYM,MULD2H,ICMO
C   infinp.h: DODFT,HSROHF,DIRFCK
C   dftcom.h: DFT
C   thrzer.h: THZER
C   maxorb.h: MXCORB (needed for inforb)
#include "maxorb.h"
#include "inforb.h"
#include "infinp.h"
#include "dftcom.h"
#include "pcmlog.h"
#include "gnrinf.h"
C
      INTEGER   NNMAX,KZVAR
      DIMENSION XTMP(2*KZVAR,NNMAX)
      DIMENSION ZYMB(NORBT,NORBT,NNMAX)
      DIMENSION CMO(*),FC(*),FCAC(*),UDV(*),FV(*),WRK(LWRK)
      DIMENSION ISYMDM(2*NNMAX),IFCTYP(2*NNMAX),KDMAO(NNMAX),
     &          KFMAO(NNMAX),
     &          KFMMO(NNMAX),KDVAO(NNMAX),KFVAO(NNMAX),KFVMO(NNMAX)
C
      DIMENSION MJWOP(*),XINDX(*)
      PARAMETER (D2 = 2.0D0)
      PARAMETER (D1 = 1.0D0)
C     local
      INTEGER KZYVAR
      REAL*8  DUMMY
      INTEGER IDUMMY
      PARAMETER ( DUMMY = 1.0D20 , IDUMMY = - 9999999 )
C
#ifdef ENABLE_XCINT_RESPONSE
      real(8), allocatable :: dmat_ao(:)
      real(8), allocatable :: xc_mat(:)
      real(8), allocatable :: fmat_mo(:)
      real(8)              :: num_electrons
      logical              :: file_exists
#endif
c
      CALL QENTER('ABS_LINE2')
C
C     Allocate work space for density and Fock matrices
C
      KFREE=1
      LFREE=LWRK
      KZYVAR=2*KZVAR
      IF (NASHT.GT.0) THEN
        NTMAX=2*NNMAX
      ELSE
        NTMAX=NNMAX
      ENDIF
c
      CALL MEMGET('REAL',KDENS,NTMAX*N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KFOCK,NTMAX*N2BASX,WRK,KFREE,LFREE)
      
c
      DO I=1,NNMAX
        KDMAO(I) = KDENS +(I-1)*N2BASX
        KFMAO(I) = KFOCK +(I-1)*N2BASX
      ENDDO
      IF (NASHT.GT.0) THEN
        DO I=1,NNMAX
          KDVAO(I)=KDENS+(NNMAX+I-1)*N2BASX
          KFVAO(I)=KFOCK+(NNMAX+I-1)*N2BASX
        ENDDO
      ENDIF
C
      CALL DZERO(WRK(KDMAO(1)),NTMAX*N2BASX)
c
      NKD=0
C
       DO I=1,NNMAX
          CALL GTZYMT(1,XTMP(1,I),KZYVAR,ISYMB,ZYMB(1,1,I),MJWOP)
C
C     Construct modified density matrices in AO basis
C
          IF (NASHT.GE.0) THEN
            CALL ABS_ZYTRA(ZYMB(1,1,I),CMO,UDV,WRK(KDMAO(I)),
     &                  WRK(KDVAO(I)),WRK(KFREE),LFREE)
          ELSE
            CALL ABS_ZYTRA(ZYMB(1,1,I),CMO,UDV,WRK(KDMAO(I)),
     &                  DUMMY,WRK(KFREE),LFREE)
          ENDIF
        END DO
C     
        IF (NISHT .GT. 0) THEN  
           CALL DSCAL(NNMAX*N2BASX,D2,WRK(KDMAO(1)),1)
        ENDIF
      IF ((NASHT.GT.0 .AND. DODFT).OR.HSROHF) THEN
         CALL DAXPY(NNMAX*N2BASX,D1,WRK(KDVAO(1)),1,WRK(KDMAO(1)),1)
         CALL DSCAL(NNMAX*N2BASX,-D1,WRK(KDVAO(1)),1)
      END IF
        IF (NASHT .GT. 0) THEN  
          DO I=1,NNMAX
            ISYMDM(NNMAX+I)=ISYMB
           ENDDO
        ENDIF
C
C Construct Fock matrices in AO basis
C
C     IFCTYP = XY
C       X indicates symmetry about diagonal
C         X = 0 No symmetry
C         X = 1 Symmetric
C         X = 2 Anti-symmetric
C       Y indicates contributions
C         Y = 0 no contribution !
C         Y = 1 Coulomb
C         Y = 2 Exchange
C         Y = 3 Coulomb + Exchange
C
C     Check if density matrix is unsymmetric (IX=0),
C     symmetric (IX=10), antisymmetric (IX=20), or zero matrix (IX=30)
C     to threshold THRZER
C
      DO I = 1,NTMAX
         IX = 10 * MATSYM(NBAST,NBAST,WRK(KDENS+(I-1)*N2BASX),THRZER)
         IF (IX .EQ. 30) THEN
C     if zero density matrix, do nothing
            IFCTYP(I) = 0
         ELSE IF (IX .EQ. 10) THEN
C     if symmetric, Coulomb+exchange
            IFCTYP(I) = IX + 3
         ELSE IF (IX .EQ. 20) THEN
C     if antisymmetric density matrix, only exchange 
            IFCTYP(I) = IX + 2
         ELSE IF (IX .EQ. 0) THEN
C     if general, Coulomb+exchange
            IFCTYP(I) = IX + 3
         ELSE
            WRITE(LUPRI,'(/A,I3)') 'ABSLIN: error'// 
     &           ' density matrix symmetry:',I
            CALL QUIT('ABSLIN: error density matrix symmetry')
         END IF
      END DO
C
      IF ((NASHT.GT.0.AND.DODFT).OR.HSROHF) THEN
C     Change IFCTYP for all Fa matrices 
        DO I=1,NNMAX  
c          IF (TRPLET) THEN
c             IFCTYP(I)=10*(IFCTYP(I)/10) + 2
c             IFCTYP(NNMAX+I)=10*(IFCTYP(NNMAX+I)/10) + 3
c          ELSE
             IFCTYP(I)=10*(IFCTYP(I)/10) + 3
             IFCTYP(NNMAX+I)=10*(IFCTYP(NNMAX+I)/10) + 2
c          END IF
        END DO
      END IF  
C
      CALL DZERO(WRK(KFMAO(1)),NTMAX*N2BASX)
      CALL SIRFCK(WRK(KFMAO(1)),WRK(KDMAO(1)),NTMAX,ISYMDM,IFCTYP,
     &     DIRFCK,WRK(KFREE),LFREE)
C
      IF ((NASHT.GT.0 .AND. DODFT).OR.HSROHF) THEN
C
C Calculated matrices are (FC+FV,-FV(exch)=FV-Q)
C Input fock matrices (form SIRIFC) are of the form
C (FC+Q,FV-Q), where Q=FV(coul) + 2*FV(exch)
C adapt to this form (which works in RSPORB)
C 
C FC+Q 
         CALL DAXPY(NNMAX*N2BASX,-D1,WRK(KFVAO(1)),1,WRK(KFMAO(1)),1)
      END IF
C Transform Fock matrices to MO basis
C Reuse memory used for density matrices
C
      DO I=1,NNMAX
          KFMMO(I) = KDENS +(I-1)*NORBT*NORBT
      ENDDO
      IF (NASHT.GT.0) THEN
        DO I=1,NNMAX
          KFVMO(I) = KDENS +(NNMAX+I-1)*NORBT*NORBT
        ENDDO
      ENDIF
C     
      NKD=0 
      NASIM=4
      KINT=INT(NNMAX/NASIM)
      KMOD=MOD(NNMAX,NASIM)
      KTOT=KINT
      IF (KMOD.NE.0) THEN
        KTOT=KINT+1
      ENDIF
      DO K=1,KTOT
        IF (K.GT.KINT) NASIM=KMOD
        DO I=1,NASIM
c        DO I=1,NNMAX
          NKD=NKD+1
          CALL FAOMO(ISYMB,WRK(KFMAO(NKD)),WRK(KFMMO(NKD)),CMO,
     &        WRK(KFREE),LFREE)
          IF (NASHT .GT.0) THEN
            CALL DZERO(WRK(KFVMO(NKD)),NORBT*NORBT)
            DO 300 ISYM1=1,NSYM
              ISYM2  = MULD2H(ISYM1,ABS_GRADSYM)
              NORB1  = NORB(ISYM1)
              NORB2  = NORB(ISYM2)
              IF (NORB1.EQ.0 .OR. NORB2.EQ.0) GO TO 300
              CALL AUTPV(ISYM1,ISYM2,CMO(ICMO(ISYM1)+1),
     &                   CMO(ICMO(ISYM2)+1),WRK(KFVAO(NKD)),NBAS,NBAST,
     &                   WRK(KFVMO(NKD)),NORB,NORBT,WRK(KFREE),LFREE)
  300       CONTINUE
          END IF
        ENDDO
c        add DFT contribution
        IF (DODFT) THEN
C         Special case of alpha density only (old code)  
          IF((NISHT.EQ.0) .AND. (NASHT.GT.0)) THEN
            DO IOSIM = 1, NASIM
              call dft_lin_respab(WRK(KFMMO(NKD-NASIM+IOSIM)),
     &                            WRK(KFVMO(NKD-NASIM+IOSIM)),CMO,
     &                            ZYMB(1,1,NKD-NASIM+IOSIM),.FALSE.,
     &                            ISYMB,WRK(KFREE),LFREE,IPRDFT)
            END DO
          ENDIF 
C         Proceed normal way otherwise 
          IF ((NASHT.GT.0) .AND. (NISHT.GT.0)) THEN 
            call dft_lin_respab_b(NASIM,WRK(KFMMO(NKD-NASIM+1)),
     &                            WRK(KFVMO(NKD-NASIM+1)),CMO,
     &                            ZYMB(1,1,NKD-NASIM+1),.FALSE.,
     &                            ISYMB,WRK(KFREE),LFREE,IPRDFT)
          ELSE
#ifdef ENABLE_XCINT_RESPONSE
             if (xcint_is_active()) then

c             if (nasim > 1) then
c                call quit('NASIM > 1 not implemented with ENABLE_XCINT')
c             end if

             inquire(file='num_grid', exist=file_exists)
             if (.not. file_exists) then
                call quit('ERROR: num_grid file not found')
             end if
             allocate(dmat_ao(nbast*nbast*2))
             allocate(xc_mat(nbast*nbast))
             allocate(fmat_mo(norbt*norbt))
             do iosim =1,nasim
               dmat_ao = 0.0d0
               xc_mat = 0.0d0
               fmat_mo = 0.0d0
             ! get unperturbed ao density matrix
               call dft_get_ao_dens_mat(cmo, dmat_ao,
     &                                  wrk(kfree), lfree)
             ! get perturbed ao density matrix
               call deq27(cmo, zymb(1,1,nkd-nasim+iosim), (/0.0d0/), 
     &                    dmat_ao(nbast*nbast + 1),
     &                    (/0.0d0/), wrk(kfree), lfree)

             ! scale density matrices
               call dscal(nbast*nbast*2, 0.5d0, dmat_ao, 1)

               call xcint_integrate_rks_lr(dmat_ao, xc_mat)

             ! transform ao fock matrix contribution to mo
               call lrao2mo(cmo, isymb, xc_mat,
     &                      fmat_mo, wrk(kfree), lfree)

             ! scale and copy on top of fcx
!               call daxpy(norbt*norbt, 0.5d0, fmat_mo, 1, fcx, 1)
               call daxpy(norbt*norbt, 0.5d0, fmat_mo, 1,
     &                    wrk(kfmmo(nkd-nasim+iosim)), 1)
             enddo

             deallocate(dmat_ao)
             deallocate(xc_mat)
             deallocate(fmat_mo)
         else
#endif /* ENABLE_XCINT_RESPONSE */
c         call dft_lin_respf(NNMAX,WRK(KFMMO(1)),CMO,ZYMB,
c     &                 .FALSE.,ISYMB,WRK(KFREE),LFREE,IPRDFT)
            call dft_lin_respf(NASIM,WRK(KFMMO(NKD-NASIM+1)),CMO,
     &                         ZYMB(1,1,NKD-NASIM+1),
     &                        .FALSE.,ISYMB,WRK(KFREE),LFREE,IPRDFT)
#ifdef ENABLE_XCINT_RESPONSE
             endif
#endif
          END IF
        END IF
      END DO
      CALL DZERO(XTMP,KZYVAR*NNMAX)
      IF (NASHT.GT.0) THEN
        CALL MEMGET('REAL',KDUM,N2ORBX,WRK,KFREE,LFREE)
        call DZERO(WRK(KDUM),N2ORBX) 
        CALL FCKOIN(NNMAX,FC,FV,ZYMB,WRK(KFMMO(1)),WRK(KFVMO(1)))
        CALL RSPORB(.TRUE.,NNMAX,FC,WRK(KFMMO(1)),WRK(KFVMO(1)),
     &              WRK(KDUM),WRK(KDUM),UDV,
     &              XTMP,ZYMB)
      ELSE
        CALL FCKOIN(NNMAX,FC,DUMMY,ZYMB,WRK(KFMMO(1)),DUMMY)
C
        CALL RSPORB(.TRUE.,NNMAX,FC,WRK(KFMMO(1)),DUMMY,
     &               DUMMY,DUMMY,UDV,XTMP,ZYMB)
c         CALL ORBEX(ISYMB,ISYMB,KZYVAR,WRK(KFMMO(I)),DUMMY,
c     &        DUMMY,DUMMY,XTMP(1,I),1.0D0,UDV,MJWOP)
      ENDIF
C
      CALL QEXIT('ABS_LINE2')
      RETURN
      END
      SUBROUTINE ABS_LINE2_INTER(KZVAR,VECS,VECB,XTMP,XTEMP,ISYMB,NLT,
     &                 NNMAX,KNVEC,CMO,FC,FCAC,FV,UDV,XINDX,MJWOP,
     &                 WRK,LWRK)
C
C PURPOSE:
C
#include "implicit.h"
#include "abslrs.h"
#include "thrzer.h"
#include "priunit.h"
C Used from
C   inforb.h: NORB,NORBT,NASHT,NISHT,N2BASX,NBAST,NBAS,NSYM,MULD2H,ICMO
C   infinp.h: DODFT,HSROHF,DIRFCK
C   dftcom.h: DFT
C   thrzer.h: THZER
C   maxorb.h: MXCORB (needed for inforb)
#include "inforb.h"
#include "maxorb.h"
#include "infinp.h"
#include "dftcom.h"
c
#include "pcmlog.h"
#include "gnrinf.h"
C
      INTEGER   NLT,NNMAX,KZVAR
      INTEGER   KNVEC(2)
      DIMENSION VECS(KZVAR,NLT)
      DIMENSION VECB(KZVAR,NLT),XTMP(2*KZVAR,NNMAX),XTEMP(2*KZVAR,NNMAX)
      DIMENSION CMO(*),FC(*),FCAC(*),UDV(*),FV(*),WRK(LWRK)
      DIMENSION ISYMDM(2*NNMAX)
C
c      DIMENSION MJWOP(2,MAXWOP,8)
      DIMENSION MJWOP(*),XINDX(*)
      PARAMETER (D2 = 2.0D0)
      PARAMETER (D1 = 1.0D0)
C     local
      INTEGER KZYVAR
      REAL*8  DUMMY
      INTEGER IDUMMY
      PARAMETER ( DUMMY = 1.0D20 , IDUMMY = - 9999999 )
C
      CALL QENTER('ABS_LINE2')
C
C     Allocate work space for density and Fock matrices
C
      KFREE=1
      LFREE=LWRK
      KZYVAR=2*KZVAR
c
      NNS=KNVEC(1)
      NNA=KNVEC(2)
C 
      NKD=0
C
C     Unpack the trial vector of specific parity
c
      DO I=1,MIN(NNS,NNA)
        NKD=NKD+1
        ISYMDM(NKD)=ISYMB
        CALL DCOPY(KZVAR,VECB(1,I),1,XTMP(1,NKD),1)
        CALL DCOPY(KZVAR,VECB(1,I),1,XTMP(KZVAR+1,NKD),1)
        CALL DAXPY(KZVAR,1.0d0,VECB(1,NNS+I),1,XTMP(1,NKD),1)
        CALL DAXPY(KZVAR,-1.0d0,VECB(1,NNS+I),1,XTMP(KZVAR+1,NKD),1)
      END DO
      IF (NNS.NE.NNA) THEN
        IF (NNS.GE.NNA) THEN
          F1=1.0d0
          NNT=0
        ELSE
          F1=-1.0d0
          NNT=NNS
        ENDIF
        DO I=MIN(NNA,NNS)+1,NNMAX
          NKD=NKD+1
          ISYMDM(NKD)=ISYMB
          CALL DCOPY(KZVAR,VECB(1,NNT+I),1,XTMP(1,NKD),1)
          CALL DCOPY(KZVAR,VECB(1,NNT+I),1,XTMP(KZVAR+1,NKD),1)
          CALL DSCAL(KZVAR,F1,XTMP(KZVAR+1,NKD),1)
        ENDDO
      ENDIF
      CALL DCOPY(2*KZVAR*NNMAX,XTMP,1,XTEMP,1)
c
      IF (ABS_BATCH) THEN
        if (nnmax .gt. abs_nbatchmax) then
          KINT=INT(NNMAX/ABS_NBATCHMAX)
          KMOD=MOD(NNMAX,ABS_NBATCHMAX) 
          KTOT=KINT
          KNMAX=ABS_NBATCHMAX
          IF (KMOD.NE.0) THEN
            KTOT=KINT+1
          ENDIF
        else
          kint = 1
          ktot = 1
          knmax = nnmax
        endif
      ELSE
         KINT=1
         KTOT=1
         KNMAX=NNMAX
      ENDIF
      KNTOT=KNMAX
      CALL MEMGET('REAL',KZYMB,KNMAX*NORBT*NORBT,WRK,KFREE,LFREE)
      DO K=1,KTOT
         IF (K.GT.KINT) KNTOT=KMOD
         IF (IPRABSLRS.GE.1) THEN
           WRITE(LUABSPRI,'(A,I4,A,I4)')
     &     'Number of linear transformations',KNTOT,' in batch ',K
         END IF
         CALL ABS_LINE2(KZVAR,XTMP(1,(K-1)*KNMAX+1),WRK(KZYMB),ISYMDM,
     &                 ISYMB,KNTOT,CMO,FC,FCAC,FV,
     &                 UDV,XINDX,MJWOP,WRK(KFREE),LFREE)
      ENDDO
C
C      IF (FLAG(16) .OR. PCM .OR. QM3 .OR. QMNPMM .OR. PEQM) THEN
      IF (EMBEDDING) THEN
         CALL RSPSLV(0,NNMAX,DUMMY,XTEMP,CMO,XINDX,UDV,
     &               XTMP,WRK(KFREE),LFREE)
      ENDIF
      IF (QMMM) THEN
        CALL QUIT('ERROR: CPP cannot be used with QMMM')
      END IF
C
      CALL DZERO(VECS,NLT*KZVAR)
      DO I=1,MIN(NNS,NNA)
         CALL DCOPY(KZVAR,XTMP(1,I),1,VECS(1,I),1)
         CALL DAXPY(KZVAR,1.0d0,XTMP(KZVAR+1,I),1,VECS(1,I),1)
         CALL DCOPY(KZVAR,XTMP(1,I),1,VECS(1,NNS+I),1)
         CALL DAXPY(KZVAR,-1.0d0,XTMP(KZVAR+1,I),1,VECS(1,NNS+I),1)
         CALL DSCAL(KZVAR,0.5d0,VECS(1,I),1)
         CALL DSCAL(KZVAR,0.5d0,VECS(1,NNS+I),1)
      ENDDO
      IF (MIN(NNS,NNA).NE.NNMAX) THEN
        IF (NNS.GE.NNA) THEN
          F1=1.0d0
          NNT=0
        ELSE
          F1=-1.0d0
          NNT=NNS
        ENDIF
        DO I=MIN(NNS,NNA)+1,NNMAX
          CALL DCOPY(KZVAR,XTMP(1,I),1,VECS(1,NNT+I),1)
          CALL DAXPY(KZVAR,F1,XTMP(KZVAR+1,I),1,VECS(1,NNT+I),1)
          CALL DSCAL(KZVAR,0.5d0,VECS(1,NNT+I),1)
         ENDDO
       ENDIF
 
      IF (IPRABSLRS.GE.200) THEN
         WRITE(LUPRI,'(//A)') ' SIGMA VECTORS (IN ORDER REAL-SYM,
     &                REAL-ANTYSYM., IMG-ANTYSYM., IMG-SYM.'
         CALL OUTPUT(VECS,1,KZVAR,1,NLT,KZVAR,NLT,1,LUPRI)
      ENDIF
C
      CALL QEXIT('ABS_LINE2')
      RETURN
      END
      SUBROUTINE ABS_ZYTRA(ZYM,CMO,DV,DXCAO,DXVAO,WRK,LWRK)
C
C
C Purpose: Computes DAO = CMO * ZYM[0,+;-,0] * transp(CMO)
C                         = CMO * transp( CMO * transp(ZYM[]) )
C          
C          ZYM[0,+;-,0](M,N) = +ZYM(M,N), M virtual, N occupied
C                            = -ZYM(M,N), N virtual, M occupied
C                            = 0, otherwhise
C
C Input:
C   ZYM(*)  kappa matrix
C   CMO(*)  MO coefficients
C   TMP(*)  temporary matrix
C
#include "implicit.h"
C
C Used from
C   maxorb.h: MXCORB (needed for inforb)
C   infinp.h: NORBT,NASHT,NISHT,NBAST,N2BASX,NSYM,MULD2H,ICMO,NORB,NBAS,NASH,NISH,IORB,IBAS,ICMO,IASH
#include "priunit.h"
#include "abslrs.h"
#include "inforb.h"
      DIMENSION ZYM(NORBT,*),WRK(LWRK)
      DIMENSION CMO(*),DV(NASHT,*),DXCAO(*),DXVAO(*)
      PARAMETER (HALF = 0.5D0, D2 = 2.0D0, D4 = 4.0D0)
      PARAMETER (HALFM = -0.5D0, D1 = 1.D0, DM1 = -1.D0)
C
      CALL QENTER('ABS_ZYTRA')
C
C loop over symmetries
C
      KFREE=1
      LFREE=LWRK
      LDXAO1 = MAX(NASHT,NISHT)*NBAST
      CALL MEMGET('REAL',KDXAO1,LDXAO1,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDXAO2,N2BASX,WRK,KFREE,LFREE)
C
C     First DXVAO, if nasht.gt.0:
C
      IF (NASHT .GT. 0) THEN
C
C     *************************************************
C     D-I and D-II:  D = D-I - D-II   Active matrices
C     *************************************************
C     DXVAO-I: loop over symmetries then put result in DXVAO
C
         CALL MEMGET('REAL',KDXV ,NASHT*NORBT,WRK,KFREE,LFREE)
         CALL DZERO(WRK(KDXAO2),N2BASX)
         DO 2000 ISYM1 = 1,NSYM
           NASH1 = NASH(ISYM1)
           NISH1 = NISH(ISYM1)
           ISYM2  = MULD2H(ISYM1,ABS_GRADSYM)
           IORB1 = IORB(ISYM1)
           IORB2 = IORB(ISYM2)
           NBAS1 = NBAS(ISYM1)
           NBAS2 = NBAS(ISYM2)
           NORB1 = NORB(ISYM1)
           NORB2 = NORB(ISYM2)
           ICMO1 = ICMO(ISYM1)
           ICMO2 = ICMO(ISYM2)
C **     step 1 of active dm:
C *****  first calculate one-index transformation of second index
C *****  of DV(uv), the active density matrix.
C        DXV(p,u) = sum(v) Bo(v,p) DV(v,u)
           IASH1 = IASH(ISYM1)
           IF (NASH1 .EQ. 0 .OR. NORB2 .EQ. 0) GO TO 2000
           CALL DGEMM('T','N',NORB2,NASH1,NASH1,1.D0,
     &                ZYM(IORB1+NISH1+1,IORB2+1),NORBT,
     &                DV(IASH1+1,IASH1+1),NASHT,0.D0,
     &                WRK(KDXV),NORB2)
C **     step 2 of active dm:
           CALL DGEMM('N','N',NBAS2,NASH1,NORB2,1.D0,CMO(ICMO2+1),
     &                NBAS2,WRK(KDXV),NORB2,0.D0,WRK(KDXAO1),NBAS2)
           IOFMOV = ICMO1 + 1 + NISH1*NBAS1
           IDXAO2 = KDXAO2 + IBAS(ISYM2)*NBAST + IBAS(ISYM1)
           CALL DGEMM('N','T',NBAS1,NBAS2,NASH1,1.D0,CMO(IOFMOV),NBAS1,
     &                WRK(KDXAO1),NBAS2,0.D0,WRK(IDXAO2),NBAST)
C **  this symmetry block finished
 2000    CONTINUE
C
         CALL DCOPY(N2BASX,WRK(KDXAO2),1,DXVAO,1)
C
C     DXVAO-II: loop over symmetries then add results to DXVAO
C               only chnage from previous loop is that first
C               MPATB is here MPAB
C
         CALL DZERO(WRK(KDXAO2),N2BASX)
         DO 3000 ISYM1 = 1,NSYM
           ISYM2  = MULD2H(ISYM1,ABS_GRADSYM)
           NASH1 = NASH(ISYM1)
           NISH1 = NISH(ISYM1)
           NASH2 = NASH(ISYM2)
           NISH2 = NISH(ISYM2)
           NBAS1 = NBAS(ISYM1)
           NBAS2 = NBAS(ISYM2)
           IORB1 = IORB(ISYM1)
           IORB2 = IORB(ISYM2)
           NORB1 = NORB(ISYM1)
           NORB2 = NORB(ISYM2)
           ICMO1 = ICMO(ISYM1)
           ICMO2 = ICMO(ISYM2)
C **     step 1 of active dm:
C *****  first calculate one-index transformation of second index
C *****  of DV(uv), the active density matrix.
C        DXV(p,u) = sum(v) Bo(v,p) DV(v,u)
           IASH1 = IASH(ISYM1)
           IASH2 = IASH(ISYM2)
           IF( NASH2 .EQ. 0 .OR. NORB1 .EQ. 0) GO TO 3000
           CALL DGEMM('N','N',NORB1,NASH2,NASH2,1.D0,
     &                ZYM(IORB1+1,IORB2+NISH2+1),NORBT,
     &                DV(IASH2+1,IASH2+1),NASHT,0.D0,WRK(KDXV),NORB1)
C **     step 2 of active dm:
           CALL DGEMM('N','N',NBAS1,NASH2,NORB1,1.D0,CMO(ICMO1+1),
     &                NBAS1,WRK(KDXV),NORB1,0.D0,WRK(KDXAO1),NBAS1)
           IOFMOV = ICMO2 + 1 + NISH2*NBAS2
           IDXAO2 = KDXAO2 + IBAS(ISYM2)*NBAST + IBAS(ISYM1)
           CALL DGEMM('N','T',NBAS1,NBAS2,NASH2,1.D0,WRK(KDXAO1),NBAS1,
     &                CMO(IOFMOV),NBAS2,0.D0,WRK(IDXAO2),NBAST)
C **  this symmetry block finished
 3000   CONTINUE
C
C subtract D_II from  D-I
        CALL DAXPY(N2BASX,DM1,WRK(KDXAO2),1,DXVAO,1)
      END IF
C
C     DXVAO finished, now DXCAO, if nisht.gt.0:
      IF (NISHT .GT. 0) THEN
        CALL DZERO(WRK(KDXAO2),N2BASX)
        DO 4000 ISYM1 = 1,ABS_NSYM
          NISH1 = NISH(ISYM1)
          IF (NISH1 .EQ. 0) GO TO 4000
          ISYM2  = MULD2H(ISYM1,ABS_GRADSYM)
          NBAS1  = NBAS(ISYM1)
          NBAS2  = NBAS(ISYM2)
          IORB1  = IORB(ISYM1)
          IORB2  = IORB(ISYM2)
          NORB1  = NORB(ISYM1)
          NORB2  = NORB(ISYM2)
          ICMO1  = ICMO(ISYM1)
          ICMO2  = ICMO(ISYM2)
C
          IF (NBAS1*NBAS2*NORB2 .NE. 0) THEN
            CALL DGEMM('N','T',NBAS2,NISH1,NORB2,1.D0,CMO(ICMO2+1),
     &                 NBAS2,ZYM(IORB1+1,IORB2+1),NORBT,0.D0,
     &                 WRK(KDXAO1),NBAS2)
            IDXAO2 = KDXAO2 + IBAS(ISYM2)*NBAST + IBAS(ISYM1)
            CALL DGEMM('N','T',NBAS1,NBAS2,NISH1,1.D0,CMO(ICMO1+1),
     &                 NBAS1,WRK(KDXAO1),NBAS2,0.D0,WRK(IDXAO2),NBAST)
          END IF
 4000   CONTINUE
        CALL DCOPY(N2BASX,WRK(KDXAO2),1,DXCAO,1)
c
        CALL DZERO(WRK(KDXAO2),N2BASX)
        DO 5000 ISYM1 = 1,NSYM      
          NISH1 = NISH(ISYM1)
          IF (NISH1 .EQ. 0) GO TO 5000
          ISYM2  = MULD2H(ISYM1,ABS_GRADSYM)
          NBAS1  = NBAS(ISYM1)
          NBAS2  = NBAS(ISYM2)
          IORB1  = IORB(ISYM1)
          IORB2  = IORB(ISYM2)
          NORB1  = NORB(ISYM1)
          NORB2  = NORB(ISYM2)
          ICMO1  = ICMO(ISYM1)
          ICMO2  = ICMO(ISYM2)
          IF (NBAS1*NBAS2*NORB2 .NE. 0) THEN
            CALL DGEMM('N','N',NBAS2,NISH1,NORB2,1.D0,CMO(ICMO2+1),
     &                 NBAS2,ZYM(IORB2+1,IORB1+1),NORBT,0.D0,
     &                 WRK(KDXAO1),NBAS2)
            IDXAO2 = KDXAO2 + IBAS(ISYM1)*NBAST + IBAS(ISYM2)
            CALL DGEMM('N','T',NBAS2,NBAS1,NISH1,1.D0,WRK(KDXAO1),NBAS2,
     &                 CMO(ICMO1+1),NBAS1,0.D0,WRK(IDXAO2),NBAST)
          END IF
C
 5000   CONTINUE
C
        CALL DAXPY(N2BASX,-1.0d0,WRK(KDXAO2),1,DXCAO,1)
      ENDIF
C
      IF ( IPRABSLRS.GE.200 ) THEN
         WRITE (LUPRI,'(/A)')
     &   'Final result in ABS_ZYTRA'
         CALL OUTPUT(DXCAO,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
      END IF
C 
      CALL QEXIT('ABS_ZYTRA')
C
      RETURN
      END
      SUBROUTINE ABSLRSRESULT(MJWOP,RESLRF,WRK,LWRK)
C
#include "implicit.h"
#include "priunit.h"
C
C Used from
C   inforb.h: NSYM
C   codata.h: XTEV
C    qrinf.h: MZYVAR (kzyvar=mzyvar(isym))
C   infvar.h: MAXWOP
C   infdim.h: NVARMA (when ABS_ANALYZE)
#include "abslrs.h"
#include "inforb.h"
#include "codata.h"
#include "qrinf.h"
#include "infvar.h"
#include "infdim.h"
#include "absorp.h"
!Information on nuclei/nuclear moments required for NSCD/NSOR. Sonia
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "chrnos.h"
#include "nuclei.h"
#include "symmet.h"
#include "cbiher.h"
#include "spnout.h"

C
      LOGICAL FOUND, CONV
C
      DIMENSION MJWOP(2,MAXWOP,8),RESLRF(2,ABS_NFREQ_INTERVAL,3,3,2)
      DIMENSION WRK(LWRK),ISRCORBR(5),ISRCORBI(5)
      !Magneto-chiral dichroism averages. Janusz
      DIMENSION AREAABB(NFREQ_BETA_B),AREABBA(NFREQ_BETA_B)
      DIMENSION AIMAABB(NFREQ_BETA_B),AIMABBA(NFREQ_BETA_B)
      DIMENSION GREABG(NFREQ_BETA_B),GIMABG(NFREQ_BETA_B)
      DIMENSION BIREFR(NFREQ_BETA_B),DICHRO(NFREQ_BETA_B)
      !Nuclear spin CD/OR averages. Sonia
      CHARACTER*8 LABINT, LABX, LABY, LABZ
      DOUBLE PRECISION XNSCD(NPATOM), XNSOR(NPATOM)
C
      KFREE = 1
      LFREE = LWRK
C
      CALL TITLER('FINAL RESULTS FOR ABSORPTION','*',120)
C
      CALL AROUND('Polarizability with damping')
      CALL PRSYMB(LUPRI,'-',66,1)
      WRITE(LUPRI,'(A3,2A10,A12,2A16)') 
     &     ' No','A-oper','B-oper','Frequency','Real part','Imag part'
      CALL PRSYMB(LUPRI,'-',66,1)
      DO IFREQ=1,ABS_NFREQ_INTERVAL
         IF (ABSLRS_INTERVAL) THEN
            FREQ = ABS_FREQ_INTERVAL(1) + (IFREQ-1)*ABS_FREQ_INTERVAL(3)
         ELSE
            FREQ = ABS_FREQ_ALPHA(IFREQ)
         END IF
         WRITE(LUPRI,'(3(/,I3,2A10,F12.6,2F16.6))') 
     &        6*IFREQ-5,'XDIPLEN','XDIPLEN',FREQ,
     &        RESLRF(1,IFREQ,1,1,1),RESLRF(2,IFREQ,1,1,1),
     &        6*IFREQ-4,'YDIPLEN','YDIPLEN',FREQ,
     &        RESLRF(1,IFREQ,2,2,1),RESLRF(2,IFREQ,2,2,1),
     &        6*IFREQ-3,'ZDIPLEN','ZDIPLEN',FREQ,
     &        RESLRF(1,IFREQ,3,3,1),RESLRF(2,IFREQ,3,3,1)
         WRITE(LUPRI,'(3(I3,2A10,F12.6,2F16.6,/))') 
     &        6*IFREQ-2,'XDIPLEN','YDIPLEN',FREQ,
     &        RESLRF(1,IFREQ,1,2,1),RESLRF(2,IFREQ,1,2,1),
     &        6*IFREQ-1,'XDIPLEN','ZDIPLEN',FREQ,
     &        RESLRF(1,IFREQ,1,3,1),RESLRF(2,IFREQ,1,3,1),
     &        6*IFREQ,'YDIPLEN','ZDIPLEN',FREQ,
     &        RESLRF(1,IFREQ,2,3,1),RESLRF(2,IFREQ,2,3,1)
         WRITE(LUPRI,'(6X,A14,3X,F12.6,2F16.6)') 
     &        'Averaged value',FREQ,
     &        (RESLRF(1,IFREQ,1,1,1)+RESLRF(1,IFREQ,2,2,1)+
     &         RESLRF(1,IFREQ,3,3,1))/3,
     &        (RESLRF(2,IFREQ,1,1,1)+RESLRF(2,IFREQ,2,2,1)+
     &         RESLRF(2,IFREQ,3,3,1))/3
      END DO
      WRITE(LUPRI,'(A)')' '
      CALL PRSYMB(LUPRI,'-',66,1)
      WRITE (LUPRI,'(A,F10.6,A,F7.4,A,F7.1,A)')
     &     ' Damping parameter equals',
     &     ABS_DAMP,' au =',
     &     ABS_DAMP*XTEV,' eV =',
     &     ABS_DAMP*XTKAYS,' cm-1'
c
      IF (ABSLRS_BETA) THEN
         CALL AROUND('First-order hyperpolarizability'//
     &        ' with damping')
      ELSE IF (ABSLRS_MCD) THEN
         CALL AROUND('Magnetic circular dichroism'//
     &        ' with damping')
      END IF
c
!-------------------------------------------------------------------!
!     Magneto Optical Activity (CD and OR)                          !
!-------------------------------------------------------------------!
c
      IF (ABSLRS_BETA .OR. ABSLRS_MCD) THEN
         CALL PRSYMB(LUPRI,'-',94,1)
         WRITE(LUPRI,'(A3,6A10,2A16)')
     &        ' No','A-oper','B-oper','C-oper','A-freq',
     &        'B-freq','C-freq','Real part','Imag part'
         CALL PRSYMB(LUPRI,'-',94,1)
         BTMP = 99.9D9
         CTMP = 99.9D9
         BTERM = 0.0
         RMORD = 0.0
         DO IQRF=1,NQRF
           IF (BTMP.NE.QRFFRQ(IQRF,1).OR.CTMP.NE.QRFFRQ(IQRF,2)) THEN
             IF ((ABSLRS_MCD).AND.(.NOT.IQRF.EQ.1)) THEN
               WRITE(LUPRI,'(5X,A21,7X,3F10.6,2F16.6)')
     &               'Orientational average',(QRFFRQ(IQRF-1,I),I=1,3),
     &                BTERM,RMORD
               BTERM = 0.0
               RMORD = 0.0
             END IF
             WRITE(LUPRI,'(A)') ' '
             BTMP = QRFFRQ(IQRF,1)
             CTMP = QRFFRQ(IQRF,2)
           END IF
           IF (ABSLRS_BETA) THEN
              WRITE(LUPRI,'(I5,3A10,3F10.6,2F16.6)') IQRF,
     &             (QRFLAB(IQRF,I), I=1,3),(QRFFRQ(IQRF,I), I=1,3),
     &             (RES_BETA(IQRF,I), I=1,2)
           ELSE IF (ABSLRS_MCD) THEN
             WRITE(LUPRI,'(I3,3A10,3F10.6,2F16.6)') IQRF,
     &            (QRFLAB(IQRF,I), I=1,3),(QRFFRQ(IQRF,I), I=1,3),
     &             RES_BETA(IQRF,2),RES_BETA(IQRF,1)
C
C            Add contribution to orientational average
C
             IF ((QRFLAB(IQRF,1)(:1).EQ.'X')
     &            .AND.(QRFLAB(IQRF,2)(:1).EQ.'Y')
     &            .AND.(QRFLAB(IQRF,3)(:1).EQ.'Z')) THEN
               BTERM = BTERM - RES_BETA(IQRF,2)
               RMORD = RMORD - RES_BETA(IQRF,1)
             ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Y')
     &                .AND.(QRFLAB(IQRF,2)(:1).EQ.'Z')
     &                .AND.(QRFLAB(IQRF,3)(:1).EQ.'X')) THEN
               BTERM = BTERM - RES_BETA(IQRF,2)
               RMORD = RMORD - RES_BETA(IQRF,1)
             ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Z')
     &                .AND.(QRFLAB(IQRF,2)(:1).EQ.'X')
     &                .AND.(QRFLAB(IQRF,3)(:1).EQ.'Y')) THEN
               BTERM = BTERM - RES_BETA(IQRF,2)
               RMORD = RMORD - RES_BETA(IQRF,1)
             ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Z')
     &                .AND.(QRFLAB(IQRF,2)(:1).EQ.'Y')
     &                .AND.(QRFLAB(IQRF,3)(:1).EQ.'X')) THEN
               BTERM = BTERM + RES_BETA(IQRF,2)
               RMORD = RMORD + RES_BETA(IQRF,1)
             ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'X')
     &                .AND.(QRFLAB(IQRF,2)(:1).EQ.'Z')
     &                .AND.(QRFLAB(IQRF,3)(:1).EQ.'Y')) THEN
               BTERM = BTERM + RES_BETA(IQRF,2)
               RMORD = RMORD + RES_BETA(IQRF,1)
             ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Y')
     &                .AND.(QRFLAB(IQRF,2)(:1).EQ.'X')
     &                .AND.(QRFLAB(IQRF,3)(:1).EQ.'Z')) THEN
               BTERM = BTERM + RES_BETA(IQRF,2)
               RMORD = RMORD + RES_BETA(IQRF,1)
             END IF
           END IF
         END DO
         WRITE(LUPRI,'(5X,A21,7X,3F10.6,2F16.6)')
     &                 'Orientational average',(QRFFRQ(NQRF,I), I=1,3),
     &                 BTERM,RMORD
         WRITE(LUPRI,'(A)')' '
         CALL PRSYMB(LUPRI,'-',94,1)
         WRITE (LUPRI,'(A,F10.6,A,F7.4,A,F7.1,A)')
     &        ' Damping parameter equals',
     &        DAMPING,' au =',
     &        DAMPING*XTEV,' eV =',
     &        DAMPING*XTKAYS,' cm-1'
      END IF
!-------------------------------------------------------------------!
!     Nuclear Spin Optical Activity (CD and OR) - S. Coriani, 2013  !
!-------------------------------------------------------------------!
      IF (ABSLRS_NSCD) THEN

        write(lupri,*)'NPATOM=', NPATOM
        write(lupri,*) (IPATOM(I),I = 1, NPATOM)

        CALL AROUND('Nuclear Spin optical rotation/circular dichroism'//
     &        ' with damping')
         CALL PRSYMB(LUPRI,'-',95,1)
         WRITE(LUPRI,'(A3,6A10,2A16)')
     &        ' No','A-oper','B-oper','C-oper','A-freq',
     &        'B-freq','C-freq','Real part','Imag part'
         CALL PRSYMB(LUPRI,'-',95,1)
         BTMP = 99.9D9
         CTMP = 99.9D9
         !NPATOM is the number of nuclei requested in input for which
         !PSO has been computed
         CALL DZERO(XNSCD,NPATOM)
         CALL DZERO(XNSOR,NPATOM)
         !Sonia: note that RES_BETA has opposite sign than std
         !code since scaled by -1. I change the sign again 
         !here to adapt  to old output
         DO IQRF=1,NQRF
            WRITE(LUPRI,'(I3,3A10,3F10.6,2F16.6)') IQRF,
     &            (QRFLAB(IQRF,I), I=1,3),(QRFFRQ(IQRF,I), I=1,3),
     &            -RES_BETA(IQRF,2),-RES_BETA(IQRF,1)
         END DO

         CALL PRSYMB(LUPRI,'=',97,1)
         WRITE(LUPRI,'(A10,2X,3A10,2X,2A14,2X,2A17)')
     &        ' Averages ','  A-freq  ','  B-freq  ','  C-freq  ',
     &        ' 6*Re<<A,B,C>> ',' 6*Im<<A,B,C>> ',' CD(mu rad/M cm) ',
     &        ' OR(mu rad/M cm) ' 
         CALL PRSYMB(LUPRI,'=',97,1)

         DO IDX=1,NPATOM
           IWHICH=IPATOM(IDX)
           !quick fix to resume the info on Gval of atom
           do ii=1,6
              xmmom = DISOTP(IZATOM(IWHICH),ii,'MMOM')
              xspin = DISOTP(IZATOM(IWHICH),ii,'SPIN')
              if (ABS(xmmom).gt.1.0D-4) exit
           end do
           Unit_Factor = xmmom*(1.8687D0)   !(micro rad / (M cm))
           !---------------------------------------------------

         DO IQRF=1,NQRF
            IF (BTMP.NE.QRFFRQ(IQRF,1).OR.CTMP.NE.QRFFRQ(IQRF,2)) THEN
               IF (.NOT.IQRF.EQ.1) THEN
                  WRITE(LUPRI,'(1X,A5,A2,A1,1X,3F10.6,4F16.6)')
     &                 'Atom(',NAMN(IWHICH),')',
     &                 (QRFFRQ(IQRF-1,I),I=1,3),
     &                 XNSCD(IDX),XNSOR(IDX),
     &                 XNSCD(IDX)*Unit_Factor*QRFFRQ(IQRF-1,2),
     &                 XNSOR(IDX)*Unit_Factor*QRFFRQ(IQRF-1,2)
                  XNSCD(IDX) = 0.0D0
                  XNSOR(IDX) = 0.0D0
               END IF
               WRITE(LUPRI,'(A)') ' '
               BTMP = QRFFRQ(IQRF,1)
               CTMP = QRFFRQ(IQRF,2)
            END IF
C           Add contribution to orientational average
            !write (lupri,*)'Name of nucleus is ', NAMN(IDX)
            !WARNING: FOLLOWING ONLY TESTED FOR NO SYMMETRY CASE
             KSYMOP = 1
             !ISCORX = IPTCNT(3*(IDX-1)+1,KSYMOP-1,2)
             !ISCORY = IPTCNT(3*(IDX-1)+2,KSYMOP-1,2)
             !ISCORZ = IPTCNT(3*(IDX-1)+3,KSYMOP-1,2)
             ISCORX = IPTCNT(3*(IWHICH-1)+1,KSYMOP-1,2)
             ISCORY = IPTCNT(3*(IWHICH-1)+2,KSYMOP-1,2)
             ISCORZ = IPTCNT(3*(IWHICH-1)+3,KSYMOP-1,2)
             IF (ISCORX .GT. 0) THEN
               IFIRST = ISCORX/100
               ISECND = MOD(ISCORX,100)/10
               ITHIRD = MOD(MOD(ISCORX,100),10)
               LABX = 'PSO '//CHRNOS(IFIRST)//
     &                          CHRNOS(ISECND)//
     &                          CHRNOS(ITHIRD)//' '
!               write(lupri,*)'Test label PSO X:', LABX
             END IF
             IF (ISCORY .GT. 0) THEN
               IFIRST = ISCORY/100
               ISECND = MOD(ISCORY,100)/10
               ITHIRD = MOD(MOD(ISCORY,100),10)
               LABY = 'PSO '//CHRNOS(IFIRST)//
     &                          CHRNOS(ISECND)//
     &                          CHRNOS(ITHIRD)//' '
               !write(lupri,*)'Test label PSO Y:', LABY
             END IF
             IF (ISCORZ .GT. 0) THEN
               IFIRST = ISCORZ/100
               ISECND = MOD(ISCORZ,100)/10
               ITHIRD = MOD(MOD(ISCORZ,100),10)
               LABZ = 'PSO '//CHRNOS(IFIRST)//
     &                          CHRNOS(ISECND)//
     &                          CHRNOS(ITHIRD)//' '
               !write(lupri,*)'Test label PSO Z:', LABZ
             END IF

!             write(lupri,*)' QRFLAB(IQRF,3)(1:8)',
!     &                       QRFLAB(IQRF,3)(1:8)

             IF ((QRFLAB(IQRF,1)(:1).EQ.'X')
     &              .AND.(QRFLAB(IQRF,2)(:1).EQ.'Y')
     &              .AND.(QRFLAB(IQRF,3)(1:8).EQ.LABZ)) THEN

!                  write(lupri,*)' XYZ '
                  XNSCD(IDX) = XNSCD(IDX) - RES_BETA(IQRF,2)
                  XNSOR(IDX) = XNSOR(IDX) - RES_BETA(IQRF,1)

               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Y')
     &                 .AND.(QRFLAB(IQRF,2)(:1).EQ.'Z')
     &                 .AND.(QRFLAB(IQRF,3)(1:8).EQ.LABX)) THEN
!                  write(lupri,*)' YZX '
                  XNSCD(IDX) = XNSCD(IDX) - RES_BETA(IQRF,2)
                  XNSOR(IDX) = XNSOR(IDX) - RES_BETA(IQRF,1)

               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Z')
     &                 .AND.(QRFLAB(IQRF,2)(:1).EQ.'X')
     &                 .AND.(QRFLAB(IQRF,3)(1:8).EQ.LABY)) THEN
!                  write(lupri,*)' ZXY '
                  XNSCD(IDX) = XNSCD(IDX) - RES_BETA(IQRF,2)
                  XNSOR(IDX) = XNSOR(IDX) - RES_BETA(IQRF,1)

               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Z')
     &                 .AND.(QRFLAB(IQRF,2)(:1).EQ.'Y')
     &                 .AND.(QRFLAB(IQRF,3)(1:8).EQ.LABX)) THEN
!                  write(lupri,*)' ZYX '
                  XNSCD(IDX) = XNSCD(IDX) + RES_BETA(IQRF,2)
                  XNSOR(IDX) = XNSOR(IDX) + RES_BETA(IQRF,1)
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'X')
     &                 .AND.(QRFLAB(IQRF,2)(:1).EQ.'Z')
     &                 .AND.(QRFLAB(IQRF,3)(1:8).EQ.LABY)) THEN
!                  write(lupri,*)' XZY '
                  XNSCD(IDX) = XNSCD(IDX) + RES_BETA(IQRF,2)
                  XNSOR(IDX) = XNSOR(IDX) + RES_BETA(IQRF,1)
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Y')
     &                 .AND.(QRFLAB(IQRF,2)(:1).EQ.'X')
     &                 .AND.(QRFLAB(IQRF,3)(1:8).EQ.LABZ)) THEN
!                  write(lupri,*)' YXZ '
                  XNSCD(IDX) = XNSCD(IDX) + RES_BETA(IQRF,2)
                  XNSOR(IDX) = XNSOR(IDX) + RES_BETA(IQRF,1)
               END IF
         END DO
         WRITE(LUPRI,'(1X,A5,A2,A1,1X,3F10.6,4F16.6)')
     &                 'Atom(',NAMN(IWHICH),')',
     &                  (QRFFRQ(NQRF,I), I=1,3),
     &                 XNSCD(IDX),XNSOR(IDX),
     &                 XNSCD(IDX)*Unit_Factor*QRFFRQ(NQRF,2),
     &                 XNSOR(IDX)*Unit_Factor*QRFFRQ(NQRF,2)
         END DO
         WRITE(LUPRI,'(A)')' '
         CALL PRSYMB(LUPRI,'-',94,1)
         WRITE (LUPRI,'(A,F10.6,A,F7.4,A,F7.1,A)')
     &        ' Damping parameter equals',
     &        DAMPING,' au =',
     &        DAMPING*XTEV,' eV =',
     &        DAMPING*XTKAYS,' cm-1'
      END IF
! -------------------------------------------------------------------
!     Magneto-Chiral Dichroism and Birefringence (J. Cukras, 2013)
! -------------------------------------------------------------------

      IF ( ABSLRS_MCHD ) THEN
         LUFU=-1
         CALL GPOPEN(LUFU,'MCHD_QRF.list','NEW',' ','FORMATTED',
     &   IDUMMY,.FALSE.)
         CALL AROUND('Magneto-chiral circular dichroism'//
     &        ' with damping')
         CALL PRSYMB(LUFU,'-',94,1)
         WRITE(LUFU,'(A3,6A10,2A16)')
     &        ' No','A-oper','B-oper','C-oper','A-freq',
     &        'B-freq','C-freq','G Real part','G Imag part'
         CALL PRSYMB(LUFU,'-',94,1)
         BTMP = 99.9D9
         CTMP = 99.9D9
         GTERMRE = 0.0
         GTERMIM = 0.0
         A1TERMRE = 0.0
         A1TERMIM = 0.0
         A2TERMRE = 0.0
         A2TERMIM = 0.0
         IFREQ=0
C Loop for G         
         DO IQRF=1,NQRF
            IF (BTMP.NE.QRFFRQ(IQRF,1).OR.CTMP.NE.QRFFRQ(IQRF,2)) THEN
               IF (.NOT.IQRF.EQ.1) THEN
                  WRITE(LUFU,'(5X,A21,7X,3F10.6,2F16.6)')
     &          'G_abg orientational average',(QRFFRQ(IQRF-1,I),I=1,3),
     &           GTERMRE,GTERMIM
                IFREQ=IFREQ+1
                GREABG(IFREQ) = GTERMRE
                GIMABG(IFREQ) = GTERMIM
                GTERMRE = 0.0
                GTERMIM = 0.0
               END IF
               WRITE(LUFU,'(A)') ' '
               BTMP = QRFFRQ(IQRF,1)
               CTMP = QRFFRQ(IQRF,2)
            END IF
            WRITE(LUFU,'(I5,3A10,3F10.6,2F16.6)') IQRF,
     &              (QRFLAB(IQRF,I), I=1,3),(QRFFRQ(IQRF,I), I=1,3),
     &              RES_BETA(IQRF,1),RES_BETA(IQRF,2)

C Contribution to orientational avarage for
C G_{\alpha,\beta,\gamma}

               IF ((QRFLAB(IQRF,1)(:1).EQ.'X')
     &              .AND.(QRFLAB(IQRF,2)(:2).EQ.'YA')
     &              .AND.(QRFLAB(IQRF,3)(:1).EQ.'Z')) THEN
                  GTERMRE = GTERMRE - RES_BETA(IQRF,1)
                  GTERMIM = GTERMIM - RES_BETA(IQRF,2)
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Y')
     &                 .AND.(QRFLAB(IQRF,2)(:2).EQ.'ZA')
     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'X')) THEN
                  GTERMRE = GTERMRE - RES_BETA(IQRF,1)
                  GTERMIM = GTERMIM - RES_BETA(IQRF,2)
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Z')
     &                 .AND.(QRFLAB(IQRF,2)(:2).EQ.'XA')
     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'Y')) THEN
                  GTERMRE = GTERMRE - RES_BETA(IQRF,1)
                  GTERMIM = GTERMIM - RES_BETA(IQRF,2)
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Z')
     &                 .AND.(QRFLAB(IQRF,2)(:2).EQ.'YA')
     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'X')) THEN
                  GTERMRE = GTERMRE + RES_BETA(IQRF,1)
                  GTERMIM = GTERMIM + RES_BETA(IQRF,2)
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'X')
     &                 .AND.(QRFLAB(IQRF,2)(:2).EQ.'ZA')
     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'Y')) THEN
                  GTERMRE = GTERMRE + RES_BETA(IQRF,1)
                  GTERMIM = GTERMIM + RES_BETA(IQRF,2)
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Y')
     &                 .AND.(QRFLAB(IQRF,2)(:2).EQ.'XA')
     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'Z')) THEN
                  GTERMRE = GTERMRE + RES_BETA(IQRF,1)
                  GTERMIM = GTERMIM + RES_BETA(IQRF,2)
               END IF

         END DO
         WRITE(LUFU,'(5X,A21,7X,3F10.6,2F16.6)')
     &          'G_abg orientational average',(QRFFRQ(NQRF,I), I=1,3),
     &           GTERMRE,GTERMIM
              IFREQ=IFREQ+1
              GREABG(IFREQ) = GTERMRE
              GIMABG(IFREQ) = GTERMIM

!         CALL PRSYMB(LUPRI,'-',94,1)
!         WRITE(LUPRI,'(A3,6A10,2A16)') 
!     &        ' No','A-oper','B-oper','C-oper','A-freq',
!     &        'B-freq','C-freq','A Real part','A Imag part'
!         CALL PRSYMB(LUPRI,'-',94,1)
         BTMP = 99.9D9
         CTMP = 99.9D9
         IFREQ=0
C Loop for A
         DO IQRF=1,NQRF
            IF (BTMP.NE.QRFFRQ(IQRF,1).OR.CTMP.NE.QRFFRQ(IQRF,2)) THEN
               IF (.NOT.IQRF.EQ.1) THEN
                WRITE(LUFU,'(5X,A21,7X,3F10.6,2F16.6)')
     &        'A_aabb orientational average',(QRFFRQ(IQRF-1,I),I=1,3),
     &         A1TERMRE,A1TERMIM
               WRITE(LUFU,'(5X,A21,7X,3F10.6,2F16.6)')
     &         'A_abba orientational average',(QRFFRQ(IQRF-1,I),I=1,3),
     &         A2TERMRE,A2TERMIM
                  IFREQ=IFREQ+1
                  AREAABB(IFREQ) = A1TERMRE
                  AIMAABB(IFREQ) = A1TERMIM
                  AREABBA(IFREQ) = A2TERMRE
                  AIMABBA(IFREQ) = A2TERMIM
                  A1TERMRE = 0.0
                  A1TERMIM = 0.0
                  A2TERMRE = 0.0
                  A2TERMIM = 0.0
               END IF
!               WRITE(LUPRI,'(A)') ' '
               BTMP = QRFFRQ(IQRF,1)
               CTMP = QRFFRQ(IQRF,2)
            END IF
            WRITE(LUFU,'(I5,3A10,3F10.6,2F16.6)') IQRF,
     &              (QRFLAB(IQRF,I), I=1,3),(QRFFRQ(IQRF,I), I=1,3),
     &              RES_BETA(IQRF,2),RES_BETA(IQRF,1)
! Orientational avarage form
! A_{\alpha,\alpha\beta,\beta}         
               IF ((QRFLAB(IQRF,1)(:1).EQ.'X')
     &              .AND.(QRFLAB(IQRF,2)(:2).EQ.'XX')
     &              .AND.(QRFLAB(IQRF,3)(:1).EQ.'X')) THEN
                  A1TERMRE = A1TERMRE - RES_BETA(IQRF,2)
                  A1TERMIM = A1TERMIM - RES_BETA(IQRF,1)
C Part of the abba                  
                  A2TERMRE = A2TERMRE - RES_BETA(IQRF,2)
                  A2TERMIM = A2TERMIM - RES_BETA(IQRF,1)
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Y')
     &                 .AND.(QRFLAB(IQRF,2)(:2).EQ.'YY')
     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'Y')) THEN
                  A1TERMRE = A1TERMRE - RES_BETA(IQRF,2)
                  A1TERMIM = A1TERMIM - RES_BETA(IQRF,1)
C Part of the abba                                    
                  A2TERMRE = A2TERMRE - RES_BETA(IQRF,2)
                  A2TERMIM = A2TERMIM - RES_BETA(IQRF,1)
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Z')
     &                 .AND.(QRFLAB(IQRF,2)(:2).EQ.'ZZ')
     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'Z')) THEN
                  A1TERMRE = A1TERMRE - RES_BETA(IQRF,2)
                  A1TERMIM = A1TERMIM - RES_BETA(IQRF,1)
C Part of the abba                  
                  A2TERMRE = A2TERMRE - RES_BETA(IQRF,2)
                  A2TERMIM = A2TERMIM - RES_BETA(IQRF,1)
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'X')
     &                 .AND.(QRFLAB(IQRF,2)(:2).EQ.'XY')
     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'Y')) THEN
                  A1TERMRE = A1TERMRE - RES_BETA(IQRF,2)
                  A1TERMIM = A1TERMIM - RES_BETA(IQRF,1)
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Y')
     &                 .AND.(QRFLAB(IQRF,2)(:2).EQ.'YZ')
     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'Z')) THEN
                  A1TERMRE = A1TERMRE - RES_BETA(IQRF,2)
                  A1TERMIM = A1TERMIM - RES_BETA(IQRF,1)
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'X')
     &                 .AND.(QRFLAB(IQRF,2)(:2).EQ.'XZ')
     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'Z')) THEN
                  A1TERMRE = A1TERMRE - RES_BETA(IQRF,2)
                  A1TERMIM = A1TERMIM - RES_BETA(IQRF,1)
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Y')
     &                 .AND.(QRFLAB(IQRF,2)(:2).EQ.'XY')
     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'X')) THEN
                  A1TERMRE = A1TERMRE - RES_BETA(IQRF,2)
                  A1TERMIM = A1TERMIM - RES_BETA(IQRF,1)
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Z')
     &                 .AND.(QRFLAB(IQRF,2)(:2).EQ.'YZ')
     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'Y')) THEN
                  A1TERMRE = A1TERMRE - RES_BETA(IQRF,2)
                  A1TERMIM = A1TERMIM - RES_BETA(IQRF,1)
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Z')
     &                 .AND.(QRFLAB(IQRF,2)(:2).EQ.'XZ')
     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'X')) THEN
                  A1TERMRE = A1TERMRE - RES_BETA(IQRF,2)
                  A1TERMIM = A1TERMIM - RES_BETA(IQRF,1)
C Orientational avarage form
C A_{\alpha,\beta\beta,\alpha}         
C               ELSE IF ((QRFLAB(IQRF,1)(:1).EQ.'X')
C     &              .AND.(QRFLAB(IQRF,2)(:2).EQ.'XX')
C     &              .AND.(QRFLAB(IQRF,3)(:1).EQ.'X')) THEN
C                  A2TERMRE = A2TERMRE - RES_BETA(IQRF,2)
C                  A2TERMIM = A2TERMIM - RES_BETA(IQRF,1)
C               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Y')
C     &                 .AND.(QRFLAB(IQRF,2)(:2).EQ.'YY')
C     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'Y')) THEN
C                  A2TERMRE = A2TERMRE - RES_BETA(IQRF,2)
C                  A2TERMIM = A2TERMIM - RES_BETA(IQRF,1)
C               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Z')
C     &                 .AND.(QRFLAB(IQRF,2)(:2).EQ.'ZZ')
C     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'Z')) THEN
C                  A2TERMRE = A2TERMRE - RES_BETA(IQRF,2)
C                  A2TERMIM = A2TERMIM - RES_BETA(IQRF,1)
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'X')
     &                 .AND.(QRFLAB(IQRF,2)(:2).EQ.'YY')
     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'X')) THEN
                  A2TERMRE = A2TERMRE - RES_BETA(IQRF,2)
                  A2TERMIM = A2TERMIM - RES_BETA(IQRF,1)
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Y')
     &                 .AND.(QRFLAB(IQRF,2)(:2).EQ.'ZZ')
     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'Y')) THEN
                  A2TERMRE = A2TERMRE - RES_BETA(IQRF,2)
                  A2TERMIM = A2TERMIM - RES_BETA(IQRF,1)
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'X')
     &                 .AND.(QRFLAB(IQRF,2)(:2).EQ.'ZZ')
     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'X')) THEN
                  A2TERMRE = A2TERMRE - RES_BETA(IQRF,2)
                  A2TERMIM = A2TERMIM - RES_BETA(IQRF,1)
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Y')
     &                 .AND.(QRFLAB(IQRF,2)(:2).EQ.'XX')
     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'Y')) THEN
                  A2TERMRE = A2TERMRE - RES_BETA(IQRF,2)
                  A2TERMIM = A2TERMIM - RES_BETA(IQRF,1)
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Z')
     &                 .AND.(QRFLAB(IQRF,2)(:2).EQ.'YY')
     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'Z')) THEN
                  A2TERMRE = A2TERMRE - RES_BETA(IQRF,2)
                  A2TERMIM = A2TERMIM - RES_BETA(IQRF,1)
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Z')
     &                 .AND.(QRFLAB(IQRF,2)(:2).EQ.'XX')
     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'Z')) THEN
                  A2TERMRE = A2TERMRE - RES_BETA(IQRF,2)
                  A2TERMIM = A2TERMIM - RES_BETA(IQRF,1)
               END IF
         END DO

            WRITE(LUFU,'(5X,A21,7X,3F10.6,2F16.6)')
     &           'A_aabb orientational average',(QRFFRQ(NQRF,I),I=1,3),
     &            A1TERMRE,A1TERMIM
            WRITE(LUFU,'(5X,A21,7X,3F10.6,2F16.6)')
     &           'A_abba orientational average',(QRFFRQ(NQRF,I),I=1,3),
     &            A2TERMRE,A2TERMIM
                  IFREQ=IFREQ+1
                  AREAABB(IFREQ) = A1TERMRE
                  AIMAABB(IFREQ) = A1TERMIM
                  AREABBA(IFREQ) = A2TERMRE
                  AIMABBA(IFREQ) = A2TERMIM

         WRITE(LUFU,'(A)')' '
         CALL PRSYMB(LUFU,'-',94,1)
         WRITE (LUPRI,'(A,F10.6,A,F7.4,A,F7.1,A)')
     &        ' Damping parameter equals',
     &        DAMPING,' au =',
     &        DAMPING*XTEV,' eV =',
     &        DAMPING*XTKAYS,' cm-1'
         WRITE(LUPRI,*) ''
         CALL PRSYMB(LUPRI,'-',94,1)
!          WRITE(LUPRI,'(A13,A13,A13,A13,A13,A13,A13)') 'FREQ',
!      &   'Re(G_abg)','Im(G_abg)','Re(A_aabb)','Im(A_aabb)',
!      &   "B(A')+B(G)",'Dichro'
         WRITE(LUPRI,'(A5,A10,A13,A13,A13,A13,A13,A13,A13)') '',
     &   'DAMPING','FREQ',
     &   'Re(G_abg)','Im(G_abg)',
     &   'Re(A_aabb)','Im(A_aabb)',
     &   "B(A')+B(G)",'Dichro'
         CALL PRSYMB(LUPRI,'-',94,1)
         DO IFREQ=1,NFREQ_BETA_B
         BIREFR(IFREQ)=(FREQ_BETA_B(IFREQ)/5.0D0)*AIMAABB(IFREQ)*
     &   0.5D0+
     &   GREABG(IFREQ)*0.25D0
         DICHRO(IFREQ)=(FREQ_BETA_B(IFREQ)/5.0D0)*AREAABB(IFREQ)*0.5D0+
     &   GIMABG(IFREQ)*0.25D0
!          WRITE(LUPRI,'(F13.6,F13.6,F13.6,F13.6,F13.6,F13.6,F13.6)')
!      &   FREQ_BETA_B(IFREQ),0.25*GREABG(IFREQ),
!      &   0.25*GIMABG(IFREQ),0.5*AREAABB(IFREQ),0.5*AIMAABB(IFREQ),
!      &   BIREFR(IFREQ),
!      &   DICHRO(IFREQ)
!          WRITE(LUPRI,'(A15,F13.6)') "DBGJC B(A')=",
!      &   0.5*AIMAABB(IFREQ)*
!      &   FREQ_BETA_B(IFREQ)/5
         WRITE(LUPRI,'(A5,F10.6,F13.6,F13.6,
     &   F13.6,F13.6,F13.6,F13.6,F13.6)')
     &   "DBGJC",DAMPING,FREQ_BETA_B(IFREQ),0.25D0*GREABG(IFREQ),
     &   0.25D0*GIMABG(IFREQ),0.5D0*AREAABB(IFREQ),0.5D0*AIMAABB(IFREQ),
     &   BIREFR(IFREQ),
     &   DICHRO(IFREQ)
         WRITE(LUPRI,'(A15,F13.6)') "DBGJC B(A')=",
     &   0.5D0*AIMAABB(IFREQ)*
     &   FREQ_BETA_B(IFREQ)/5.0D0
         END DO
         CALL GPCLOSE(LUFU,'KEEP')
      END IF

C ------------------------------------------------------------------

      IF (ABSLRS_ANALYZE) THEN
         MZYVMX = 2*NVARMA
         CALL MEMGET('REAL',KVEC,2*MZYVMX,WRK,KFREE,LFREE)
         CALL AROUND('Analysis of response vectors')
         DO ISYM=1,NSYM
           IF (ABS_NOPER(ISYM).GT.0) THEN
             KZYVAR=MZYVAR(ISYM)
             DO IOPER=1,ABS_NOPER(ISYM)
               DO IFREQ=1,ABS_NFREQ_INTERVAL
                 IF (ABSLRS_INTERVAL) THEN
                   FREQ = ABS_FREQ_INTERVAL(1) + 
     &                    (IFREQ-1)*ABS_FREQ_INTERVAL(3)
                 ELSE
                   FREQ = ABS_FREQ_ALPHA(IFREQ)
                 END IF
                 WRITE(LUPRI,'(/A11,A10,2(/,A11,F10.6))') 'Property :', 
     &              ABS_LABOP(IOPER,ISYM),
     &              'Frequency:',ABS(FREQ),
     &              'Damping  :',ABS_DAMP
                 CALL READ_XVEC(LUABSVECS,2*KZYVAR,WRK(KVEC),
     &              ABS_LABOP(IOPER,ISYM),ISYM,
     &              ABS(FREQ),ABS_THCLR,FOUND,CONV)
                 IF (.NOT. FOUND) THEN
                   WRITE(LUPRI,'(/A)') 
     &                 ' Response vector not found on file LUABSVECS'
                   CALL QUIT('Response vector not found on file')
                 ELSE IF(.NOT. CONV) THEN
                   WRITE (LUPRI,'(/A)') ' @WARNING: '//
     &                 'Response vector not converged on file LUABSVECS'
                 END IF
C
                 DNORM_RE=DNRM2(KZYVAR,WRK(KVEC),1)
                 DNORM_IM=DNRM2(KZYVAR,WRK(KVEC+KZYVAR),1)
C
                 IF (IPRABSLRS.GT.1) THEN
                   WRITE(LUPRI,'(/A,2F14.8)') 
     &                 ' Norm of vector (real and imag):',
     &                 DNORM_RE, DNORM_IM
C
                    WRITE(LUPRI,'(/7X,2A5,2(2X,2(5X,A2,5X)))')
     &                 'occ','vir','ZR','YR','ZI','YI'
                    DO I=1,KZYVAR/2
                      WRITE(LUPRI,'(I5,2X,2I5,2(2X,2F12.8))') 
     &                    I,MJWOP(1,I,ISYM),MJWOP(2,I,ISYM),
     &                    WRK(KVEC-1+I),WRK(KVEC+KZYVAR/2-1+I),
     &                    WRK(KVEC+KZYVAR-1+I),WRK(KVEC+3*KZYVAR/2-1+I)
                    END DO
C
                 END IF

                 WRITE(LUPRI,'(/A,2F14.8)') 
     &              ' Norm of vector (real and imag):',
     &              DNORM_RE, DNORM_IM
C
C     Sum occupied orbital contributions into aggregates
C     Sonia 2 Joanna: I changed MAXOCC to MAXOCCU because MAXOCC is
C     a constant in one of the new common blocks I had to include for NSCD
C
               MAXOCCU=0
               DO I=1,KZYVAR/2
                  MAXOCCU = MAX(MAXOCCU,MJWOP(1,I,ISYM))
               ENDDO
               CALL MEMGET('REAL',KORBWR,MAXOCCU,WRK,KFREE,LFREE)
               CALL MEMGET('REAL',KORBWI,MAXOCCU,WRK,KFREE,LFREE)
c     real part
               DO I=1,MAXOCCU
                  WRK(KORBWR+I-1) = 0.0D0
               ENDDO
               DO I=1,KZYVAR/2
                  WRK(KORBWR+MJWOP(1,I,ISYM)-1) = 
     &                 WRK(KORBWR+MJWOP(1,I,ISYM)-1) +
     &                 WRK(KVEC-1+I)**2 + WRK(KVEC+KZYVAR/2-1+I)**2
               ENDDO
c     imag part
               DO I=1,MAXOCCU
                  WRK(KORBWI+I-1) = 0.0D0
               ENDDO
               DO I=1,KZYVAR/2
                  WRK(KORBWI+MJWOP(1,I,ISYM)-1) = 
     &                 WRK(KORBWI+MJWOP(1,I,ISYM)-1) +
     &                 WRK(KVEC+KZYVAR-1+I)**2 + 
     &                 WRK(KVEC+3*KZYVAR/2-1+I)**2
               ENDDO
               NSRC=5
               CALL FINDMAXN(WRK(KORBWR),MAXOCCU,ISRCORBR,NSRC)
               NSRC=5
               CALL FINDMAXN(WRK(KORBWI),MAXOCCU,ISRCORBI,NSRC)
               WRITE (LUPRI,*)
     &          'Important occupied orbital contributions (normalized)'
               WRITE (LUPRI,*)
     &          ' occ  real    occ  imag'
               DO I=1,NSRC
                  WRITE (LUPRI,'(I5,F7.3,I6,F7.3)') ISRCORBR(I),
     &                 WRK(KORBWR+ISRCORBR(I)-1)/DNORM_RE**2,
     &                 ISRCORBI(I),
     &                 WRK(KORBWI+ISRCORBI(I)-1)/DNORM_IM**2
               ENDDO
            END DO
         END DO
      END IF
      END DO
      END IF
C     
      RETURN
      END
      SUBROUTINE GET_LRF(LU,LABEL,NFREQ_ABS,FREQ_ABS,
     &     FC,FV,CMO,UDV,PV,XINDX,RESLRF,KZVAR,
     &     WRK,LWRK)
#include "implicit.h"
#include "priunit.h"
#include "abslrs.h"
C
      CHARACTER*8 LABEL,LABGD
      LOGICAL FOUND,CONV
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),XINDX(*),
     &     WRK(LWRK),RESLRF(2,ABS_NFREQ_INTERVAL,3,3,2),
     &     FREQ_ABS(NMXFREQ)
      parameter (d0=-1.D0-16)
C
      KFREE = 1
      LFREE = LWRK
C
      CALL DIPLAB(LABEL,INDEX)
      IF (LABEL(2:8).EQ.'DIPLEN') THEN
         ITYPE = 1
      ELSE IF (LABEL(2:8).EQ.'ANGMOM') THEN
         ITYPE = 2
      ELSE IF (LABEL(3:7).EQ.'THETA') THEN
         ITYPE = 3
      ELSE
         WRITE(LUPRI,'(A)') ' Warning: Unknown operator!'
         RETURN
      END IF
C
      KZYVAR=2*KZVAR
      CALL MEMGET('REAL',KVEC,2*KZYVAR,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KGD,KZYVAR,WRK,KFREE,LFREE)
C
      DO IOPER=1,ABS_NOPER(ABS_GRADSYM)
         LABGD = ABS_LABOP(IOPER,ABS_GRADSYM)
         IF (LABEL(2:8).EQ.LABGD(2:8)) THEN
            CALL DIPLAB(LABGD,INXGD)
            CALL GETGPV(LABGD,FC,FV,CMO,UDV,PV,XINDX,ANTSYM,
     &           WRK(KGD),LFREE)
            DO IFREQ=1,ABS_NFREQ_INTERVAL
                 CALL READ_XVEC(LU,2*KZYVAR,WRK(KVEC),
     &              LABEL,ABS_GRADSYM,
c     &              ABS(FREQ_ABS(IFREQ),ABS_THCLR,FOUND,CONV)
     &              abs(FREQ_ABS(IFREQ)),ABS_THCLR,FOUND,CONV)
               if (freq_abs(ifreq) .lt. d0) then
                 call abs_vec_swap(2*kzyvar,wrk(kvec))
               endif
               RESLRF(1,IFREQ,INXGD,INDEX,ITYPE) =
     &              DDOT(KZYVAR,WRK(KGD),1,WRK(KVEC),1)
               RESLRF(2,IFREQ,INXGD,INDEX,ITYPE) =
     &              DDOT(KZYVAR,WRK(KGD),
     &              1,WRK(KVEC+KZYVAR),1)
            END DO
         END IF
      END DO
C
      RETURN
      END
      SUBROUTINE ABS_QRRDVE(ISYMA,ISYMB,ISYMC,ALAB,BLAB,CLAB,
     &     FREQA,FREQB,FREQC,KZYVA,KZYVB,KZYVC,VECA,VECB,VECC)
C
C PURPOSE: Read in response vectors for quadratic response. 
C
#include "implicit.h"
#include "priunit.h"
#include "absorp.h"
#include "abslrs.h"
#include "inftap.h"
#include "rspprp.h"
#include "inflr.h"
#include "qrinf.h"
C
      LOGICAL FOUND, CONV
      CHARACTER*8 ALAB,BLAB,CLAB,BLANK
      PARAMETER ( D0=0.0D0, DM1=-1.0D0 )
      DIMENSION VECA(*),VECB(*),VECC(*)
C
      KZYVA  = MZYVAR(ISYMA)
      KZYVB  = MZYVAR(ISYMB)
      KZYVC  = MZYVAR(ISYMC)
C
      IF (IPRABS.GE.2) THEN
         WRITE(LUPRI,'(2(/A),2(/A,3(I10)),/A,3A10,/A,3F10.6)')
     &         ' Variables in QRRDVE',
     &         ' ==================================================== ',
     &         ' KZYVA,KZYVB,KZYVC: ',KZYVA,KZYVB,KZYVC,
     &         ' ISYMA,ISYMB,ISYMC: ',ISYMA,ISYMB,ISYMC,
     &         ' ALAB,BLAB,CLAB   : ',ALAB,BLAB,CLAB,
     &         ' FREQA,FREQB,FREQC: ',FREQA,FREQB,FREQC
      END IF
C
C     Read in Na
C
      CALL READ_XVEC(LUABSVECS,2*KZYVA,VECA,ALAB,ISYMA,
     &              ABS(FREQA),ABS_THCLR,FOUND,CONV)
      IF (.NOT. (FOUND .AND. CONV)) THEN
         IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/3A,F7.5,A,I3,/A)') ' Response label ',ALAB,
     &           ' with frequency ',FREQA, ' and symmetry',
     &           ISYMA,' not found on file LUABSVECS'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/3A,F7.5,/A,I3,A)') ' @WARNING>>>>'//
     &           ' Response label ',ALAB,
     &           ' with frequency ',FREQA, ' and symmetry',
     &           ISYMA,' not converged on file LUABSVECS'
         END IF
      END IF
      IF (FREQA .GT. D0) THEN
         CALL DSWAP(KZYVA/2,VECA,1,VECA(1+KZYVA/2),1)
         CALL DSWAP(KZYVA/2,VECA(1+KZYVA),1,VECA(1+KZYVA+KZYVA/2),1)
         CALL DSCAL(KZYVA,DM1,VECA,1)
      END IF
C
      IF (IPRABS.GE.10) THEN
         WRITE(LUPRI,'(/A)') ' Na  vector in ABS_READVEC '
         CALL PRSYMB(LUPRI,'=',72,1)
         WRITE(LUPRI,'(4X,4A15)') 'ZR','YR','ZI','YI'
         CALL OUTPUT(VECA,1,KZYVA/2,1,4,KZYVA/2,4,1,LUPRI)
      END IF
C
C     Read in Nb
C
      CALL READ_XVEC(LUABSVECS,2*KZYVB,VECB,BLAB,ISYMB,
     &              ABS(FREQB),ABS_THCLR,FOUND,CONV)
      IF (.NOT. (FOUND .AND. CONV)) THEN
         IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/3A,F7.5,A,I3,/A)') ' Response label ',BLAB,
     &           ' with frequency ',FREQB, ' and symmetry',
     &           ISYMB,' not found on file LUABSVECS'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/3A,F7.5,/A,I3,A)') ' @WARNING>>>>'//
     &           ' Response label ',BLAB,
     &           ' with frequency ',FREQB, ' and symmetry',
     &           ISYMB,' not converged on file LUABSVECS'
         END IF
      END IF
      IF (FREQB .LT. D0) THEN
         CALL DSWAP(KZYVB/2,VECB,1,VECB(1+KZYVB/2),1)
         CALL DSWAP(KZYVB/2,VECB(1+KZYVB),1,VECB(1+KZYVB+KZYVB/2),1)
         CALL DSCAL(KZYVB,DM1,VECB,1)
      END IF
C
      IF (IPRABS.GE.10) THEN
         WRITE(LUPRI,'(/A)') ' Nb  vector in ABS_READVEC '
         CALL PRSYMB(LUPRI,'=',72,1)
         WRITE(LUPRI,'(4X,4A15)') 'ZR','YR','ZI','YI'
         CALL OUTPUT(VECB,1,KZYVB/2,1,4,KZYVB/2,4,1,LUPRI)
      END IF
C
C     Read in Nc
C
      CALL READ_XVEC(LUABSVECS,2*KZYVC,VECC,CLAB,ISYMC,
     &              ABS(FREQC),ABS_THCLR,FOUND,CONV)
      IF (.NOT. (FOUND .AND. CONV)) THEN
         IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/3A,F7.5,A,I3,/A)') ' Response label ',CLAB,
     &           ' with frequency ',FREQC, ' and symmetry',
     &           ISYMC,' not found on file LUABSVECS'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/3A,F7.5,/A,I3,A)') ' @WARNING>>>>'//
     &           ' Response label ',CLAB,
     &           ' with frequency ',FREQC, ' and symmetry',
     &           ISYMC,' not converged on file LUABSVECS'
         END IF
      END IF
      IF (FREQC .LT. D0) THEN
         CALL DSWAP(KZYVC/2,VECC,1,VECC(1+KZYVC/2),1)
         CALL DSWAP(KZYVC/2,VECC(1+KZYVC),1,VECC(1+KZYVC+KZYVC/2),1)
         CALL DSCAL(KZYVC,DM1,VECC,1)
      END IF
C
      IF (IPRABS.GE.10) THEN
         WRITE(LUPRI,'(/A)') ' Nc  vector in ABS_READVEC '
         CALL PRSYMB(LUPRI,'=',72,1)
         WRITE(LUPRI,'(4X,4A15)') 'ZR','YR','ZI','YI'
         CALL OUTPUT(VECC,1,KZYVC/2,1,4,KZYVC/2,4,1,LUPRI)
      END IF
C
      RETURN
      END
      SUBROUTINE ABS_BETA_SETUP
#include "implicit.h"
#include "priunit.h"
#include "inforb.h"
#include "absorp.h"
#include "abslrs.h"
C
      PARAMETER (THRZERO = 1.0D-6)
C
      LOGICAL DOHYP
      CHARACTER*8 ALAB,BLAB,CLAB
C
      NQRF = 0
      NFREQ_ALPHA = 0
      NFREQ_BETA_B=ABS_NFREQ_BETA_B
      NFREQ_BETA_C=ABS_NFREQ_BETA_C
      DO I=1,ABS_NFREQ_BETA_B
        FREQ_BETA_B(I)=ABS_FREQ_BETA_B(I)
      ENDDO
      DO I=1,ABS_NFREQ_BETA_C
        FREQ_BETA_C(I)=ABS_FREQ_BETA_C(I)
      ENDDO
      ABS_ALPHA=ABSLRS_ALPHA
      ABS_BETA=ABSLRS_BETA
      ABS_MCD=ABSLRS_MCD
      DAMPING=ABS_DAMP
C
      DO ICFR=1,ABS_NFREQ_BETA_C
        DO IBFR=1,ABS_NFREQ_BETA_B
          FREQB = ABS_FREQ_BETA_B(IBFR)
          FREQC = ABS_FREQ_BETA_C(ICFR)
          FREQA = -(FREQB+FREQC)
          IF (ABSLRS_SHG .AND. .NOT.FREQB.EQ.FREQC) GOTO 100
          DO ISYMA=1,NSYM
            DO ISYMB=1,NSYM
              ISYMC = MULD2H(ISYMA,ISYMB)
              IF (ABS_NOPER(ISYMA).GT.0 .AND. ABS_NOPER(ISYMB).GT.0 
     &          .AND. ABS_NOPER(ISYMC).GT.0) THEN
                DO IAOP=1,ABS_NOPER(ISYMA)
                  DO IBOP=1,ABS_NOPER(ISYMB)
                    DO ICOP=1,ABS_NOPER(ISYMC)
                      ALAB = ABS_LABOP(IAOP,ISYMA)
                      BLAB = ABS_LABOP(IBOP,ISYMB)
                      CLAB = ABS_LABOP(ICOP,ISYMC)
                      CALL ABS_QRCHK(DOHYP,ALAB,BLAB,CLAB,
     &                  ISYMA,ISYMB,ISYMC,FREQA,FREQB,FREQC)
                      IF (DOHYP) THEN
                        NQRF = NQRF +1
                        QRFLAB(NQRF,1) = ALAB
                        QRFLAB(NQRF,2) = BLAB
                        QRFLAB(NQRF,3) = CLAB
                        QRFSYM(NQRF,1) = ISYMA
                        QRFSYM(NQRF,2) = ISYMB
                        QRFSYM(NQRF,3) = ISYMC
                        QRFFRQ(NQRF,1) = FREQA
                        QRFFRQ(NQRF,2) = FREQB
                        QRFFRQ(NQRF,3) = FREQC
                      END IF
                    END DO
                  END DO
                END DO
              END IF
            END DO
          END DO
 100     CONTINUE
        END DO
      END DO
C
      CALL GPDSRT(NFREQ_ALPHA,FREQ_ALPHA,THRZERO)
c      ABS_NFREQ_ALPHA = NFREQ_ALPHA
c      DO I=1,ABS_NFREQ_ALPHA
c        ABS_FREQ_ALPHA(I)=FREQ_ALPHA(I)
c      ENDDO
!      IF (IPRABS.GE.0) THEN
         CALL AROUND('Setup of Hyperpolarizability Calculation')
         WRITE (LUPRI,'(2(/A),I4,A)')
     & ' This calculations requires the solution of linear response',
     & ' equations for electric dipole operators at',ABS_NFREQ_ALPHA,
     & ' frequencies:'
         WRITE(LUPRI,'(/A,5(4F12.8,/,9X))')
     &        ' LR FREQ:',(ABS_FREQ_ALPHA(I),I=1,ABS_NFREQ_ALPHA)
         WRITE (LUPRI,'(/A,I5,A)')
     & ' and the evaluation of',NQRF,' quadratic response functions:'
         WRITE(LUPRI,'(/2A,/A3,6A12,/2A)')
     &        '--------------------------------------',
     &        '-------------------------------------',
     &        ' No','A-oper','B-oper','C-oper',
     &        'A-freq','B-freq','C-freq',
     &        '--------------------------------------',
     &        '-------------------------------------'
         DO IQRF=1,NQRF
            WRITE(LUPRI,'(I4,3A12,3F12.8)') IQRF,
     &           (QRFLAB(IQRF,I), I=1,3),(QRFFRQ(IQRF,I), I=1,3)
         END DO
         WRITE(LUPRI,'(2A)')
     &        '--------------------------------------',
     &        '-------------------------------------'
!      END IF
C
C     End of BETA_SETUP
C
      RETURN
      END
      SUBROUTINE ABS_QRCHK(DOHYP,ALAB,BLAB,CLAB,ISYMA,ISYMB,ISYMC,
     &                 FREQA,FREQB,FREQC)
#include "implicit.h"
#include "priunit.h"
#include "absorp.h"
#include "abslrs.h"
C
      PARAMETER (THRZERO = 1.0D-6)
      LOGICAL DOHYP,NEWFRQ
      CHARACTER*8 ALAB,BLAB,CLAB,LAB(3)
      DIMENSION FREQ(3),ISYM(3)
C
      DOHYP = .TRUE.
C
C     Operator requirement for magnetic circular dichroism
C
      IF (ABS_MCD .OR. ABSLRS_MCD) THEN
         IF (ALAB(2:8).NE.'DIPLEN' .OR. BLAB(2:8).NE.'DIPLEN' .OR.
     &       CLAB(2:8).NE.'ANGMOM') DOHYP = .FALSE.
      END IF
C
C     Operator requirement for MCHD
C
      IF (ABSLRS_MCHD) THEN
         IF (ALAB(2:8).NE.'DIPLEN' .OR.
     &       (BLAB(2:8).NE.('ANGMOM') .AND. BLAB(3:8).NE.'THETA') .OR.
     &       CLAB(2:8).NE.'ANGMOM') DOHYP = .FALSE.
      END IF
C
C     Operator requirement for NSCD
C
      IF (ABSLRS_NSCD) THEN
         IF (ALAB(2:8).NE.'DIPLEN' .OR. BLAB(2:8).NE.('DIPLEN') .OR.
     &       CLAB(1:3).NE.'PSO') DOHYP = .FALSE.
      END IF
c
C     Check if equivalent QRF is in the list already.
c
      IF (DAMPING .EQ. 0.0D0) THEN
C     Overall permutational symmetry.
         JCTR=3
         KCTR=1
         LCTR=1
      ELSE
C     Only intrinsic permutational symmetry, so A operator has to match 
C     first operator in list of response functions.
         JCTR=1
         KCTR=2
         LCTR=2
      END IF
      DO IQRF = 1,NQRF
        DO J = 1,JCTR
          DO K = KCTR,3
            IF (K.NE.J) THEN
              DO L = LCTR,3
                IF (L.NE.K .AND. L.NE.J) THEN
C
                  IF ( ALAB.EQ.QRFLAB(IQRF,J) .AND.
     &              BLAB.EQ.QRFLAB(IQRF,K) .AND.
     &              CLAB.EQ.QRFLAB(IQRF,L) .AND.
     &              ISYMA.EQ.QRFSYM(IQRF,J) .AND.
     &              ISYMB.EQ.QRFSYM(IQRF,K) .AND.
     &              ISYMC.EQ.QRFSYM(IQRF,L) .AND.
     &              ABS( FREQA-QRFFRQ(IQRF,J)).LT.THRZERO .AND.
     &              ABS( FREQB-QRFFRQ(IQRF,K)).LT.THRZERO .AND.
     &              ABS( FREQC-QRFFRQ(IQRF,L)).LT.THRZERO ) THEN
                    DOHYP = .FALSE.
                  END IF
C
                END IF
              END DO
            END IF
          END DO
        END DO
      END DO
C
      IF (DOHYP) THEN
C
C     Check if this QRF will inflict new LR solver frequencies.
C
         FREQ(1) = FREQA
         FREQ(2) = FREQB
         FREQ(3) = FREQC
         DO I=1,3
            NEWFRQ = .TRUE.
            DO IFR=1,ABS_NFREQ_ALPHA
              IF (ABS(ABS_FREQ_ALPHA(IFR)-ABS(FREQ(I))).LT.THRZERO) THEN
                  NEWFRQ = .FALSE.
              END IF
            END DO
            IF (NEWFRQ) THEN
               IF (ABS_NFREQ_ALPHA.GE.NMXFREQ) THEN
                  WRITE(LUPRI,'(2(/A),I4,A,/A)')
     & ' The specified calculation requires more than the allowed',
     & ' number of frequencies in the LR solver NMXFREQ=',NMXFREQ,'.',
     & ' The program will stop.'
                  CALL QUIT('Too many frequencies in LR solver.')
               END IF
               ABS_NFREQ_ALPHA = ABS_NFREQ_ALPHA + 1
               ABS_FREQ_ALPHA(ABS_NFREQ_ALPHA) = ABS(FREQ(I))
            END IF
         END DO
C
      END IF
C
      IF (NQRF.GE.MXQRF .AND. DOHYP) THEN
         WRITE(LUPRI,'(2(/A),I4,A,/A)')
     & ' The specified calculation requires more than the allowed',
     & ' number of quadratic response functions MXQRF=',MXQRF,'.',
     & ' The program will stop'
         CALL QUIT('Too many quadratic response functions specified.')
      END IF
C
      RETURN
      END
      SUBROUTINE ABS_VEC_SWAP(NDIM,XVEC)
C
      DIMENSION XVEC(NDIM,4)
      REAL YVEC(NDIM,4)
      parameter (dm1=-1.0d0)

      call dzero(yvec,4*ndim)
      call dcopy(ndim,xvec(1,1),1,yvec(1,2),1)      
      call dcopy(ndim,xvec(1,2),1,yvec(1,1),1)
      call dscal(2*ndim,dm1,yvec(1,1),1)
      call dcopy(ndim,xvec(1,3),1,yvec(1,4),1)      
      call dcopy(ndim,xvec(1,4),1,yvec(1,3),1)
      call dcopy(4*ndim,yvec,1,xvec,1)      

      RETURN
      END
! end of DALTON/rsp/abscomplex.F
