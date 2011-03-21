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
C ============================================================================
C     absorp.F: Main author Patrick Norman
C     This implementation is published in:
C
C        P. Norman, D.M. Bishop, H.J.Aa. Jensen, and J. Oddershede,
C        "Near-resonant absorption in the time-dependent self-consistent
C        field and multiconfigurational self-consistent field approximations",
C        J. Chem. Phys. 115 (2001) 10323-10334.
C
C        and
C
C        P. Norman, D.M. Bishop, H.J.Aa. Jensen, and J. Oddershede,
C        "Nonlinear response theory with relaxation: the first
C        hyperpolarizability",
C        J. Chem. Phys. 123 (2005) 194103.
C ============================================================================
C
      SUBROUTINE ABSORP_INPUT(WORD)
C
C     Purpose:
C     Read in user settings for imaginary polarizabilities describing
C     absorption in the optical processes.
C
#include <implicit.h>
#include <priunit.h>
#include <gnrinf.h>
#include <absorp.h>
#include <infrsp.h>
C
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <symmet.h>
#include <codata.h>
C
      LOGICAL NEWDEF, ALLCOMP
      PARAMETER ( NTABLE = 22 )
      PARAMETER ( D0=0.0D0, D1=1.0D0, THRFRQ=1.0D-14, THD=1.0D-6 )
      CHARACTER PROMPT*1, WORD*7, WORD1*7, TABLE(NTABLE)*7
      CHARACTER*8 DIPLEN(3), ANGMOM(3)
C
      DATA DIPLEN/'XDIPLEN','YDIPLEN','ZDIPLEN'/
      DATA ANGMOM/'XANGMOM','YANGMOM','ZANGMOM'/
      DATA TABLE /'.ALPHA ','.FREQUE','.THCLR ','.MAX MI','.MAXITO',
     &            '.PRINT ','.DAMPIN','.EXCITA','.MAX MA','.THCPP ',
     &            '.BETA  ','.BFREQ ','.CFREQ ','.SHG   ','.ANALYZ',
     &            '.BATCH ','.FREQ I','.MCD   ','.REDUCE','.XXCOMP',
     &            '.YYCOMP','.ZZCOMP'/
C
      NEWDEF = (WORD .EQ. '*ABSORP')
      ICHANG = 0
C
      ALLCOMP = .TRUE.
      ICOMP=0
C
C     Set default values
C
      IF (NEWDEF) THEN
C     Common default value of the damping parameter is set to 
C     be 1000 cm-1 = 4.556333D-3 a.u.
         DAMPING = 1000/XTKAYS
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
     &                GOTO (1,2,3,4,5,6,7,8,9,10,
     &                 11,12,13,14,15,16,17,18,19,20,21,22),I
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
               ABS_ALPHA=.TRUE.
               ABS_BETA =.FALSE.
               ABS_GAMMA=.FALSE.
               GOTO 100
 2             CONTINUE
               IF (ABS_ALPHA) THEN
                  READ (LUCMD,*) NFREQ_ALPHA
                  IF (NFREQ_ALPHA.LE.MXFREQ) THEN
                     READ (LUCMD,*) (FREQ_ALPHA(J),J=1,NFREQ_ALPHA)
                  ELSE
                     WRITE (LUPRI,'(3(/,A,I5),/)')
     &               ' NUMBER OF FREQUENCIES SPECIFIED   :',NFREQ_ALPHA,
     &               ' IS GREATER THAN THE ALLOWED NUMBER:',MXFREQ,
     &               ' THE NUMBER IS RESET TO THE MAXIMUM:',MXFREQ
                     READ (LUCMD,*) (FREQ_ALPHA(J),J=1,MXFREQ),
     &                              (FFFF,J=MXFREQ+1,NFREQ_ALPHA)
                     NFREQ_ALPHA = MXFREQ
                  END IF
                  CALL GPDSRT(NFREQ_ALPHA,FREQ_ALPHA,THD)
               ELSE
                  WRITE(LUPRI,1000) WORD
                  READ(LUCMD,*)
                  READ(LUCMD,*)
               END IF
               GOTO 100
 3             CONTINUE
                  READ (LUCMD,*) THCLR_ABSORP
               GOTO 100
 4             CONTINUE
                  READ (LUCMD,*) MAX_MICRO
               GOTO 100
 5             CONTINUE
                  READ (LUCMD,*) MAX_ITORB
                  IF (MAX_ITORB .GT. 0) OPTORB = .TRUE.
               GOTO 100
 6             CONTINUE
                  READ (LUCMD,*) IPRABS
               GOTO 100
 7             CONTINUE
C     User input of damping parameter
                  READ (LUCMD,*) DAMPING
               GOTO 100
 8             CONTINUE
                  READ (LUCMD,*) NEXCITED_STATES
               GOTO 100
 9             CONTINUE
                  READ (LUCMD,*) MAX_MACRO
               GOTO 100
 10            CONTINUE
                  READ (LUCMD,*) THCPP_ABSORP
               GOTO 100
 11            CONTINUE
               ABS_ALPHA=.FALSE.
               ABS_BETA =.TRUE.
               ABS_GAMMA=.FALSE.
               GOTO 100
 12            CONTINUE
               IF (ABS_BETA .OR. ABS_MCD) THEN
                  READ (LUCMD,*) NFREQ_BETA_B
                  IF (NFREQ_BETA_B.LE.MXFREQ) THEN
                     READ (LUCMD,*) (FREQ_BETA_B(J),J=1,NFREQ_BETA_B)
                  ELSE
                     WRITE (LUPRI,'(3(/,A,I5),/)')
     &              ' NUMBER OF FREQUENCIES SPECIFIED   :',NFREQ_BETA_B,
     &              ' IS GREATER THAN THE ALLOWED NUMBER:',MXFREQ,
     &              ' THE NUMBER IS RESET TO THE MAXIMUM:',MXFREQ
                     READ (LUCMD,*) (FREQ_BETA_B(J),J=1,MXFREQ),
     &                              (FFFF,J=MXFREQ+1,NFREQ_BETA_B)
                     NFREQ_BETA_B = MXFREQ
                  END IF
               ELSE
                  WRITE(LUPRI,1010) WORD
                  READ(LUCMD,*)
                  READ(LUCMD,*)
               END IF
               GOTO 100
 13            CONTINUE
               IF (ABS_BETA) THEN
                  READ (LUCMD,*) NFREQ_BETA_C
                  IF (NFREQ_BETA_C.LE.MXFREQ) THEN
                     READ (LUCMD,*) (FREQ_BETA_C(J),J=1,NFREQ_BETA_C)
                  ELSE
                     WRITE (LUPRI,'(3(/,A,I5),/)')
     &              ' NUMBER OF FREQUENCIES SPECIFIED   :',NFREQ_BETA_C,
     &              ' IS GREATER THAN THE ALLOWED NUMBER:',MXFREQ,
     &              ' THE NUMBER IS RESET TO THE MAXIMUM:',MXFREQ
                     READ (LUCMD,*) (FREQ_BETA_C(J),J=1,MXFREQ),
     &                              (FFFF,J=MXFREQ+1,NFREQ_BETA_C)
                     NFREQ_BETA_C = MXFREQ
                  END IF
               ELSE
                  WRITE(LUPRI,1010) WORD
                  READ(LUCMD,*)
                  READ(LUCMD,*)
               END IF
               GOTO 100
 14            CONTINUE
               ABS_ALPHA=.FALSE.
               ABS_BETA =.TRUE.
               ABS_GAMMA=.FALSE.
               ABS_SHG  =.TRUE.
               GOTO 100
 15            CONTINUE
               ABS_ANALYZE=.TRUE.
               GOTO 100
 16            CONTINUE
               READ(LUCMD,*) NFREQ_BATCH
               IF (NFREQ_BATCH.GT.MXFREQ) THEN
                  NFREQ_BATCH = MXFREQ
               END IF
               GOTO 100
 17            CONTINUE
               ABS_INTERVAL = .TRUE.
               ABS_REDUCE = .TRUE.
               READ(LUCMD,*) (FREQ_INTERVAL(I), I=1,3)
               GOTO 100
 18            CONTINUE
               ABS_ALPHA=.FALSE.
               ABS_BETA =.FALSE.
               ABS_GAMMA=.FALSE.
               ABS_MCD  =.TRUE.
               GOTO 100
 19            CONTINUE
               ABS_REDUCE=.TRUE.
               GOTO 100
 20            CONTINUE
               ALLCOMP=.FALSE.
               ICOMP=1
               GOTO 100
 21            CONTINUE
               ALLCOMP=.FALSE.
               ICOMP=2
               GOTO 100
 22            CONTINUE
               ALLCOMP=.FALSE.
               ICOMP=3
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
         THCLR_ABSORP = THCLR_ABSORP*THR_REDFAC
         THCPP_ABSORP = THCPPPABSORP*THR_REDFAC
      END IF
C
C     Process user input
C
      ABSORP = ABS_ALPHA .OR. ABS_BETA .OR. ABS_GAMMA .OR. ABS_MCD
      IF (.NOT.ABSORP) THEN
         WRITE(LUPRI,'(/A)') ' Absorption input ignored because:'
         WRITE(LUPRI,'(A,L2)')' No process requested:',
     &        ' ABS_ALPHA = ABS_BETA = ABS_GAMMA = ABS_MCD = ', ABSORP
      ELSE 
C
C     Put operators in lists
C
         CALL IZERO(NOPER,8)
         DO ILABEL=1,3
            IF (ALLCOMP .OR. ICOMP.EQ.ILABEL) THEN
               ISYM = ISYMAX(ILABEL,1)+1
               NOPER(ISYM) = NOPER(ISYM) + 1
               LABOP(NOPER(ISYM),ISYM) = DIPLEN(ILABEL)
               IF (ABS_MCD) THEN
                  ISYM = ISYMAX(ILABEL,2)+1
                  NOPER(ISYM) = NOPER(ISYM) + 1
                  LABOP(NOPER(ISYM),ISYM) = ANGMOM(ILABEL)
               END IF
            END IF
         END DO
C
         CALL AROUND('Variable settings for absorption calculation')
C
         IF (ABS_ALPHA) THEN
            WRITE (LUPRI,'(A,L4)')
     &     ' Linear polarizability calculation requested:  ABS_ALPHA =',
     &           ABS_ALPHA
            IF (.NOT.ABS_INTERVAL) WRITE(LUPRI,'(A,5(4F12.8,/,16X))')
     &           ' at frequencies:',(FREQ_ALPHA(I),I=1,NFREQ_ALPHA)
         END IF
C
         IF (ABS_BETA) THEN
            WRITE (LUPRI,'(A,L4)')
     &   ' Nonlinear polarizability calculation requested:  ABS_BETA =',
     &           ABS_BETA
            WRITE(LUPRI,'(A,5(4F12.8))')
     &           ' B-FREQ:',(FREQ_BETA_B(I),I=1,NFREQ_BETA_B)
            WRITE(LUPRI,'(A,5(4F12.8))')
     &           ' C-FREQ:',(FREQ_BETA_C(I),I=1,NFREQ_BETA_C)
         END IF
C
         IF (ABS_MCD) THEN
            WRITE (LUPRI,'(A,L4)')
     &     ' Magnetic circular dichroism requested    : ABS_MCD      =',
     &           ABS_MCD
            WRITE(LUPRI,'(A,5(4F12.8))')
     &           ' at frequencies:',(FREQ_BETA_B(I),I=1,NFREQ_BETA_B)
         END IF
C
         WRITE(LUPRI,'(A,L4)')
     &      ' Absorption calc over frequency interval  : ABS_INTERVAL ='
     &       , ABS_INTERVAL
         IF (ABS_INTERVAL) THEN
            WRITE(LUPRI,'(A,I4)')
     &      ' Number of frequencies per batch          : NFREQ_BATCH  ='
     &       ,NFREQ_BATCH
            WRITE(LUPRI,'(3(A,F8.5))')
     &      ' Start:',FREQ_INTERVAL(1),'   Stop:',FREQ_INTERVAL(2),
     &      '   Step:',FREQ_INTERVAL(3)
         END IF
C
         IF (NEXCITED_STATES .GT. MXSTATES) THEN
            WRITE(LUPRI,'(/2A,/)') ' >>> Warning the number of excited',
     &           ' exceeds maximun, the maximum value is used <<<'
            NEXCITED_STATES = MXSTATES
         END IF
C
         IF (ABS_INTERVAL) THEN 
            IF( FREQ_INTERVAL(1).GT.FREQ_INTERVAL(2) .OR.
     &         (FREQ_INTERVAL(2)-FREQ_INTERVAL(1)).LT.FREQ_INTERVAL(3)
     &         .OR. FREQ_INTERVAL(3).LE.0.0D0 ) THEN
               WRITE(LUPRI,'(/,A,/,A,/,/,A,/,A,/,/,A)')
     & ' Warning: Frequency interval for absorption calculation is not',
     & '          correctly specified according to:',
     & '                 .FREQ_INTERVAL',
     & '                 STARTFREQ ENDFREQ STEPSIZE',
     & '          This input will be ignored. Calculation continues.'
               ABS_INTERVAL = .FALSE.
               FREQ_INTERVAL(1)=D0
               FREQ_INTERVAL(2)=D0
               FREQ_INTERVAL(3)=D1
            END IF
         END IF
C
         WRITE(LUPRI,'(A,F9.6,F8.1)')
     &      ' Damping parameter (a.u. and cm-1)        : DAMPING      ='
     &      ,DAMPING,DAMPING*XTKAYS
         WRITE(LUPRI,'(A,I4)')
     &      ' Number of states used in start iteration : NEXCIT       ='
     &       ,NEXCITED_STATES
         WRITE(LUPRI,'(A,1P,D8.1)')
     &      ' Threshold of convergence in LR solver    : THCLR        ='
     &      ,THCLR_ABSORP
         WRITE(LUPRI,'(A,I4)')
     &      ' Maximum iterations in complex LR solver  : MAX_MACRO    ='
     &       ,MAX_MACRO
         WRITE(LUPRI,'(A,I4)')
     &      ' Maximum iterations in real LR solver     : MAX_MICRO    ='
     &       ,MAX_MICRO
         WRITE(LUPRI,'(A,I4)')
     &      ' Max iter. in optimal orbital algorithm   : MAX_ITORB    ='
     &       ,MAX_ITORB
         WRITE(LUPRI,'(A,I4)')
     &      ' Print level in absorption modules        : IPRABS       ='
     &       ,IPRABS
      END IF
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

C
C     End of ABSORP_INPUT
C
      RETURN
      END
      SUBROUTINE ABSCALC(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     &                  XINDX,WRK,LWRK)
C
C PURPOSE:
C     Driver routine for the computation of response properties including
C     absorption
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <absorp.h>
#include <infvar.h>
#include <infrsp.h>

C
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(LWRK)
C
      IPRRSP = IPRABS - 1
      KFREE = 1
      LFREE = LWRK
C
      CALL MEMGET('REAL',KMJWOP,(2*8*MAXWOP + 1)/IRAT,WRK,
     &     KFREE,LFREE)
C
      IF (ABS_BETA .OR. ABS_MCD) THEN
         CALL BETA_SETUP
      END IF
C
C     Allocate memory for results of linear absorption calc.
C
      IF (ABS_INTERVAL) THEN
         NFREQ_INTERVAL = 1 +
     &        INT((FREQ_INTERVAL(2)-FREQ_INTERVAL(1))/FREQ_INTERVAL(3))
      ELSE
         NFREQ_INTERVAL = NFREQ_ALPHA
      END IF
C
      CALL MEMGET('REAL',KRES,2*NFREQ_INTERVAL*3*3*2,WRK,KFREE,LFREE)
C
C     Determine linear response vectors
C
      CALL AROUND('Solving Linear Response Equations')
      CALL ABSVEC1(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     &             XINDX,WRK(KMJWOP),WRK(KRES),WRK(KFREE),LFREE)
C
C     Determine nonlinear response functions
C
      IF (ABS_BETA.OR.ABS_MCD) THEN
         CALL AROUND('Evaluate Quadratic Response Functions')
         CALL ABSQRF(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     &               XINDX,WRK(KMJWOP),WRK(KFREE),LFREE)
      END IF
C
C     Print-out a summary of results
C
      CALL ABSRESULT(WRK(KMJWOP),WRK(KRES),WRK(KFREE),LFREE)
C
      RETURN
      END
      SUBROUTINE ABSVEC1(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     &                   XINDX,MJWOP,RESLRF,WRK,LWRK)
C
C PURPOSE:
C     Solve the linear response equations and store response vectors on
C     disk.
C
#include <implicit.h>
#include <priunit.h>
#include <absorp.h>
#include <infrsp.h>
#include <inforb.h>
#include <wrkrsp.h>
#include <infvar.h>
C
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),MJWOP(2,MAXWOP,8)
      DIMENSION RESLRF(2,NFREQ_INTERVAL,3,3,2),WRK(LWRK)
C
      KFREE = 1
      LFREE = LWRK
      CALL MEMGET('REAL',KREDE,MAXRM*MAXRM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KREDS,MAXRM*MAXRM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KREDZ,2*MAXRM*MAXRM,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBTYP,MAXRM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIVAL,MAXRM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KRESID,NEXCITED_STATES,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIVEC,2*MXFREQ*MAXRM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KREDGD,MAXRM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KREDZGD,2*MAXRM,WRK,KFREE,LFREE)
C
      IF (IPRABS.GT.0) CALL TIMER('START ',TIMSTR,TIMEND)
C
      DO ISYM=1,NSYM
         IF (NOPER(ISYM).GT.0) THEN
            KSYMOP = ISYM
            CALL RESONANT(WRK(KREDE),WRK(KREDS),WRK(KIBTYP),
     &           WRK(KEIVAL),WRK(KRESID),WRK(KEIVEC),CMO,UDV,PV,
     &           FOCK,FC,FV,FCAC,H2AC,XINDX,MJWOP,WRK(KFREE),LFREE)
C
            DO IOPER=1,NOPER(KSYMOP)
               CALL ABSCTL(IOPER,LABOP(IOPER,KSYMOP),
     &              WRK(KREDE),WRK(KREDS),
     &              WRK(KREDZ),WRK(KREDGD),WRK(KREDZGD),
     &              WRK(KEIVEC),WRK(KIBTYP),
     &              CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,XINDX,
     &              RESLRF,WRK(KFREE),LFREE)
            END DO
         END IF
      END DO
      IF (IPRABS.GT.0) CALL TIMER('ABSVEC',TIMSTR,TIMEND)
      RETURN
      END
      SUBROUTINE RESONANT(REDE,REDS,IBTYP,EIGVAL,RESIDUAL,EIGVEC,
     &     CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,XINDX,MJWOP,WRK,LWRK)
C
C PURPOSE:
C     Solve the eigenvalue equation ( E[2] - w*S[2] )*X = 0 for 
C     the few lowest excited states to be used as startvectors.
C
#include <implicit.h>
#include <priunit.h>
#include <absorp.h>
#include <dummy.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <infvar.h>
C
      DIMENSION REDE(*),REDS(*),IBTYP(*),EIGVAL(*),RESIDUAL(*),EIGVEC(*)
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),MJWOP(2,MAXWOP,8),WRK(LWRK)
      CHARACTER*8 BLANK
C
      PARAMETER ( BLANK='        ' )
C
      IF (IPRABS.GT.0) THEN
         WRITE(LUPRI,'(/,2A,/,2(A,I4),/A)')
     &        ' RESONANT -- Solve the eigenvalue equation',
     &        ' ( E[2] - w*S[2] )*X = 0', 
     &        ' RESONANT -- for the',NEXCITED_STATES,
     &        ' lowest excited states of symmetry',KSYMOP,
     &        ' RESONANT -- to be used as startvectors in LR solver.'
      END IF
C
      KFREE = 1
      LFREE = LWRK
C
C     Initialize variables 
C
      THCRSP = THCPP_ABSORP
      MAXIT  = MAX_MICRO
      CALL RSPVAR(UDV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK(KFREE),LFREE)
      IF (ABS_BETA .OR. ABS_MCD .OR. ABS_ANALYZE) THEN
         CALL SETZY(MJWOP)
      END IF
C
      IF (NEXCITED_STATES.EQ.0) RETURN
C
      KZRED  = 0
      KZYRED = 0
      KEXCNV = NEXCITED_STATES
      KEXSTV = KEXCNV
      KEXSIM = KEXCNV
      CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
     &            .FALSE.,BLANK,BLANK,DUMMY,DUMMY,REDE,REDS,
     &            IBTYP,EIGVAL,RESIDUAL,EIGVEC,XINDX,WRK(KFREE),LFREE)
C
      CALL DCOPY(KEXCNV,EIGVAL,1,EXC_ENERGY(1,KSYMOP),1)
C      
      RETURN
      END
      SUBROUTINE ABSREDUC(LABEL,NFREQ_ABS,FREQ_ABS,
     &     REDE,REDS,REDZ,REDGD,REDZGD,IBTYP,
     &     REDVEC,FC,FV,CMO,UDV,PV,XINDX,NBATCH,RESLRF,WRK,LWRK)
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <absorp.h>
#include <inftap.h>
#include <ibndxdef.h>
#include <qrinf.h>
C
C PURPOSE:
C     Create start vectors from the reduced eigenvalue problem
C     ( E[2] - {w+iW}*S[2] )* (NR + iNI) = B[1]
C
      CHARACTER*8 LABEL
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),XINDX(*),IBTYP(*),
     &     WRK(LWRK),REDE(MAXRM*MAXRM),REDS(MAXRM*MAXRM),
     &     REDZ(2*MAXRM*MAXRM),REDGD(MAXRM),REDZGD(2*MAXRM),
     &     REDVEC(2*MAXRM,MXFREQ),RESLRF(2,NFREQ_INTERVAL,3,3,2),
     &     FREQ_ABS(MXFREQ)
C
      KFREE = 1
      LFREE = LWRK
      CALL MEMGET('INTE',KIPIV,MAXRM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KBVEC,KZYVAR,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KGD,KZYVAR,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KREDE,KZYRED*KZYRED,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KREDS,KZYRED*KZYRED,WRK,KFREE,LFREE)
C
C Construct the reduced gradient with trial vectors from LURSP3
C
      CALL GETGPV(LABEL,FC,FV,CMO,UDV,PV,XINDX,ANTSYM,
     &     WRK(KFREE),LFREE)
      CALL DCOPY(KZYVAR,WRK(KFREE),1,WRK(KGD),1)
      CALL REWSPL(LURSP3)
      KRSP3 = 0
      IF (KOFFTY.EQ.1) THEN
#if defined (VAR_SPLITFILES)
         CALL READST(LURSP3,NRSP5,DUMMY,KRSP3,.FALSE.)
#else        
         READ (LURSP3)
#endif
      END IF
      DO I = 1, KZRED
         IF (IBTYP(KOFFTY+I).EQ.JBCNDX) THEN
#if defined (VAR_SPLITFILES)
            CALL READST(LURSP3,IREC,WRK(KBVEC),KRSP3,.TRUE.)
            IF (IREC.NE.KZCONF) CALL QUIT('Read error in ABSREDUC')
#else        
            CALL READT(LURSP3,KZCONF,WRK(KBVEC))
#endif
            PZY = DDOT(KZCONF,WRK(KBVEC),1,WRK(KGD),1)
            PYZ = DDOT(KZCONF,WRK(KBVEC),1,WRK(KGD+KZVAR),1)
         ELSE
#if defined (VAR_SPLITFILES)
            CALL READST(LURSP3,IREC,WRK(KBVEC),KRSP3,.TRUE.)
            IF (IREC.NE.KZYWOP) CALL QUIT('Read error in ABSREDUC')
#else        
            CALL READT(LURSP3,KZYWOP,WRK(KBVEC))
#endif
            PZY = DDOT(KZWOPT,WRK(KBVEC),1,WRK(KGD+KZCONF),1)
     &          + DDOT(KZWOPT,WRK(KBVEC+KZWOPT),1,
     &            WRK(KGD+KZCONF+KZVAR),1)
            PYZ = DDOT(KZWOPT,WRK(KBVEC),1,WRK(KGD+KZCONF+KZVAR),1)
     &          + DDOT(KZWOPT,WRK(KBVEC+KZWOPT),1,WRK(KGD+KZCONF),1)
         END IF 
         REDGD(2*I-1) = PZY
         REDGD(2*I)   = PYZ
      END DO
C
C     Unpack the triangular packed reduced E[2] and S[2]
C
      CALL DSPTGE(KZYRED,REDE,WRK(KREDE))
      CALL DSPTGE(KZYRED,REDS,WRK(KREDS))
C
      IF (IPRABS.GE.10) THEN
         WRITE(LUPRI,'(2(/,5X,A))') ' Reduced B[1] gradient',
     &        '========================'
         CALL OUTPUT(REDGD,1,2,1,KZRED,2,KZRED,1,LUPRI)
         WRITE(LUPRI,'(2(/,5X,A))') ' Reduced E[2] matrix',
     &        '=====================' 
         CALL OUTPAK(REDE,KZYRED,1,LUPRI)
         WRITE(LUPRI,'(2(/,5X,A))') ' Reduced S[2] matrix',
     &        '=====================' 
         CALL OUTPAK(REDS,KZYRED,1,LUPRI)
      END IF
C
      DO IFREQ=1,NFREQ_ABS
C
C     Put the reduced gradient B[1] in a complex form 
C     (with zero imaginary part) for later call to linear equation 
C     solver ZSYSV, and construct the reduced ( E[2] - {w+iW}*S[2] )
C     matrix.
C
         CALL DZERO(REDZGD,2*KZYRED)
         CALL DCOPY(KZYRED,REDGD,1,REDZGD,2)      
         CALL DZERO(REDZ,2*MAXRM*MAXRM)
         CALL DCOPY(KZYRED*KZYRED,WRK(KREDE),1,REDZ,2)
         CALL DAXPY(KZYRED*KZYRED,-FREQ_ABS(IFREQ),WRK(KREDS),1,
     &              REDZ,2)
         CALL DAXPY(KZYRED*KZYRED,-DAMPING,WRK(KREDS),1,REDZ(2),2)
C
         IF (IPRABS.GT.10) THEN
            WRITE(LUPRI,'(/,5X,A,/,5X,2(A,D10.4),/,5X,A)') 
     &           ' Reduced ( E[2] - {w+iW}*S[2] )  matrix',
     &           ' with w =', FREQ_ABS(IFREQ),' and W =', DAMPING,
     &           '========================================' 
            CALL OUTPUT(REDZ,1,2*KZYRED,1,KZYRED,2*KZYRED,KZYRED,
     &           1,LUPRI)
         END IF
C
         CALL ZSYSV('L',KZYRED,1,REDZ,KZYRED,WRK(KIPIV),REDZGD,KZYRED,
     &              WRK(KFREE),LFREE,INFO)
C
C     Store solution of response vector in reduced space in REDVEC. 
C     Real and imaginary parts of each element in the response vector
C     are stored subsequently. 
C
         CALL DCOPY(2*KZYRED,REDZGD,1,REDVEC(1,IFREQ),1)
C
         IF (IPRABS.GE.5) THEN
            WRITE(LUPRI,'(/,5X,A,/,5X,2(A,D10.4),/,5X,A)') 
     &           ' Reduced ( E[2] - {w+iW}*S[2] )-1 * B[1] vector',
     &           ' with w =', FREQ_ABS(IFREQ),' and W =', DAMPING,
     &           '================================================' 
            CALL OUTPUT(REDVEC(1,IFREQ),1,2,1,KZYRED,2,KZYRED,1,LUPRI)
         END IF
C
         IF (LABEL(2:8).EQ.'DIPLEN') THEN
            ITYPE = 1
         ELSE IF (LABEL(2:8).EQ.'ANGMOM') THEN
            ITYPE = 2
         ELSE
            WRITE(LUPRI,'(A)') ' Warning: Unknown operator!'
            GOTO 99
         END IF
         CALL DIPLAB(LABEL,INDEX)
         RESLRF(1,IFREQ+NBATCH*NFREQ_BATCH,INDEX,INDEX,ITYPE) = 
     &        DDOT(KZYRED,REDGD,1,REDVEC(1,IFREQ),2)
         RESLRF(2,IFREQ+NBATCH*NFREQ_BATCH,INDEX,INDEX,ITYPE) = 
     &        DDOT(KZYRED,REDGD,1,REDVEC(2,IFREQ),2)
 99      CONTINUE
C     
C     End of loop over frequencies
C
      END DO
C
      IF (IPRABS.GE.2) THEN
         WRITE(LUPRI,'(/3A,/A,I2,A,I4,A,/A,I7,A,//6X,A,/6X,A)') 
     &        ' >>> Value of linear polarizability for operator ',
     &        LABEL, 'of <<<',' >>> symmetry',KSYMOP,
     &        ' in the reduced space of dimension',
     &        KZYRED, '. The  <<<',
     &        ' >>> full variational space has dimension',
     &        KZYVAR,'.           <<<',
     &        ' Frequency   Damping         Real   Imaginary',
     &        '----------------------------------------------'
         DO IFREQ=1,NFREQ_ABS
            WRITE(LUPRI,'(6X,F10.4,F11.6,2F12.6)') FREQ_ABS(IFREQ),
     &         DAMPING, (RESLRF(I,IFREQ+NBATCH*NFREQ_BATCH,
     &                          INDEX,INDEX,ITYPE),I=1,2)
         END DO
         CALL FLSHFO(LUPRI)
      END IF
C
      RETURN
      END
      SUBROUTINE FINDMAXN(X,IXLEN,IMAX,NMAX)
C     Put indices to the NMAX elements of X with largest magnitude
C     in IMAX.
      IMPLICIT NONE
      DOUBLE PRECISION X
      INTEGER IXLEN,IMAX,NMAX,I,J,K,N
      DIMENSION X(IXLEN),IMAX(NMAX)
      N = 1
      IMAX(1) = 1
      DO I=2,IXLEN
C     Check if X(I) is larger than elements already stored
         J = N
         DO WHILE (J.GE.1.AND.ABS(X(I)).GT.ABS(X(IMAX(J))))
            J=J-1
         ENDDO
C     now put I after J is there is space
         IF (J.LT.NMAX) THEN
            N = MIN(NMAX,N+1)
            DO K=N,J+2,-1
               IMAX(K) = IMAX(K-1)               
            ENDDO
            IMAX(J+1) = I
         ENDIF           
      ENDDO
      NMAX = N
      END
      SUBROUTINE ABSRESULT(MJWOP,RESLRF,WRK,LWRK)
C
#include <implicit.h>
#include <priunit.h>
#include <absorp.h>
#include <inforb.h>
#include <codata.h>
#include <qrinf.h>
#include <infvar.h>
#include <inftap.h>
#include <infdim.h>
C
      LOGICAL FOUND, CONV
C
      DIMENSION MJWOP(2,MAXWOP,8),RESLRF(2,NFREQ_INTERVAL,3,3,2)
      DIMENSION WRK(LWRK),ISRCORBR(5),ISRCORBI(5)
C
      KFREE = 1
      LFREE = LWRK
C
      CALL TITLER('FINAL RESULTS FOR ABSORPTION','*',120)
      CALL AROUND('Excitation energies for dipole allowed states')
      WRITE(LUPRI,'(16X,A,/17X,A4,A6,A21,/16X,A)') 
     &     '==========================================',
     &     'Sym.','State','Excitation energy',
     &     '=========================================='
      DO ISYM=1,NSYM
         IF (NOPER(ISYM).GT.0) THEN
            DO ISTATE=1,NEXCITED_STATES
               WRITE(LUPRI,'(13X,2I6,F14.6,A,F10.4,A)') 
     &              ISYM,ISTATE,EXC_ENERGY(ISTATE,ISYM),' a.u.',
     &              EXC_ENERGY(ISTATE,ISYM)*XTEV,' eV'
            END DO
         END IF
      END DO
C
      CALL AROUND('Polarizability with damping')
      CALL PRSYMB(LUPRI,'-',66,1)
      WRITE(LUPRI,'(A3,2A10,A12,2A16)') 
     &     ' No','A-oper','B-oper','Frequency','Real part','Imag part'
      CALL PRSYMB(LUPRI,'-',66,1)
      DO IFREQ=1,NFREQ_INTERVAL
         IF (ABS_INTERVAL) THEN
            FREQ = FREQ_INTERVAL(1) + (IFREQ-1)*FREQ_INTERVAL(3)
         ELSE
            FREQ = FREQ_ALPHA(IFREQ)
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
     &     DAMPING,' au =',
     &     DAMPING*XTEV,' eV =',
     &     DAMPING*XTKAYS,' cm-1'
C
      IF (ABS_BETA) THEN
         CALL AROUND('First-order hyperpolarizability'//
     &        ' with damping')
      ELSE IF (ABS_MCD) THEN
         CALL AROUND('Magnetic circular dichroism'//
     &        ' with damping')
      END IF
      IF (ABS_BETA.OR.ABS_MCD) THEN
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
               IF ((ABS_MCD).AND.(.NOT.IQRF.EQ.1)) THEN
                  WRITE(LUPRI,'(5X,A21,7X,3F10.6,2F16.6)') 
     &                 'Orientational average',(QRFFRQ(IQRF-1,I),I=1,3),
     &                 BTERM,RMORD
                  BTERM = 0.0
                  RMORD = 0.0
               END IF
               WRITE(LUPRI,'(A)') ' '
               BTMP = QRFFRQ(IQRF,1)
               CTMP = QRFFRQ(IQRF,2)
            END IF
            IF (ABS_BETA) THEN
               WRITE(LUPRI,'(I3,3A10,3F10.6,2F16.6)') IQRF,
     &              (QRFLAB(IQRF,I), I=1,3),(QRFFRQ(IQRF,I), I=1,3),
     &              (RES_BETA(IQRF,I), I=1,2)
            ELSE IF (ABS_MCD) THEN
               WRITE(LUPRI,'(I3,3A10,3F10.6,2F16.6)') IQRF,
     &              (QRFLAB(IQRF,I), I=1,3),(QRFFRQ(IQRF,I), I=1,3),
     &              RES_BETA(IQRF,2),RES_BETA(IQRF,1)
C
C              Add contribution to orientational average
C
               IF ((QRFLAB(IQRF,1)(:1).EQ.'X')
     &              .AND.(QRFLAB(IQRF,2)(:1).EQ.'Y') 
     &              .AND.(QRFLAB(IQRF,3)(:1).EQ.'Z')) THEN
                  BTERM = BTERM - RES_BETA(IQRF,2) 
                  RMORD = RMORD - RES_BETA(IQRF,1) 
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Y')
     &                 .AND.(QRFLAB(IQRF,2)(:1).EQ.'Z') 
     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'X')) THEN
                  BTERM = BTERM - RES_BETA(IQRF,2)
                  RMORD = RMORD - RES_BETA(IQRF,1)
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Z')
     &                 .AND.(QRFLAB(IQRF,2)(:1).EQ.'X') 
     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'Y')) THEN
                  BTERM = BTERM - RES_BETA(IQRF,2) 
                  RMORD = RMORD - RES_BETA(IQRF,1) 
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Z')
     &                 .AND.(QRFLAB(IQRF,2)(:1).EQ.'Y') 
     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'X')) THEN
                  BTERM = BTERM + RES_BETA(IQRF,2) 
                  RMORD = RMORD + RES_BETA(IQRF,1) 
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'X')
     &                 .AND.(QRFLAB(IQRF,2)(:1).EQ.'Z') 
     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'Y')) THEN
                  BTERM = BTERM + RES_BETA(IQRF,2) 
                  RMORD = RMORD + RES_BETA(IQRF,1) 
               ELSE IF((QRFLAB(IQRF,1)(:1).EQ.'Y')
     &                 .AND.(QRFLAB(IQRF,2)(:1).EQ.'X') 
     &                 .AND.(QRFLAB(IQRF,3)(:1).EQ.'Z')) THEN
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
C
      IF (ABS_ANALYZE) THEN
         MZYVMX = 2*NVARMA
         CALL MEMGET('REAL',KVEC,2*MZYVMX,WRK,KFREE,LFREE)
         CALL AROUND('Analysis of response vectors')
         DO ISYM=1,NSYM
         IF (NOPER(ISYM).GT.0) THEN
            KZYVAR=MZYVAR(ISYM)
            DO IOPER=1,NOPER(ISYM)
            DO IFREQ=1,NFREQ_INTERVAL
               IF (ABS_INTERVAL) THEN
                  FREQ = FREQ_INTERVAL(1) + (IFREQ-1)*FREQ_INTERVAL(3)
               ELSE
                  FREQ = FREQ_ALPHA(IFREQ)
               END IF
               WRITE(LUPRI,'(/A11,A10,2(/,A11,F10.6))') 'Property :', 
     &              LABOP(IOPER,ISYM),
     &              'Frequency:',ABS(FREQ),
     &              'Damping  :',DAMPING
               CALL REARSP(LURSP,2*KZYVAR,WRK(KVEC),
     &              LABOP(IOPER,ISYM),'COMPLEX ',
     &              ABS(FREQ),0.0D0,ISYM,0,THCLR_ABSORP,
     &              FOUND,CONV,ANTSYM)
               IF (.NOT. FOUND) THEN
                  WRITE(LUPRI,'(/A)') 
     &                 ' Response vector not found on file LURSP'
                  CALL QUIT('Response vector not found on file')
               ELSE IF(.NOT. CONV) THEN
                  WRITE (LUPRI,'(/A)') ' @WARNING: '//
     &                 'Response vector not converged on file LURSP'
               END IF
C
               DNORM_RE=DNRM2(KZYVAR,WRK(KVEC),1)
               DNORM_IM=DNRM2(KZYVAR,WRK(KVEC+KZYVAR),1)
C
               IF (IPRABS.GT.1) THEN
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
C     Sum occupied orbital contributions into aggregates
               MAXOCC=0
               DO I=1,KZYVAR/2
                  MAXOCC = MAX(MAXOCC,MJWOP(1,I,ISYM))
               ENDDO
               CALL MEMGET('REAL',KORBWR,MAXOCC,WRK,KFREE,LFREE)
               CALL MEMGET('REAL',KORBWI,MAXOCC,WRK,KFREE,LFREE)
c     real part
               DO I=1,MAXOCC
                  WRK(KORBWR+I-1) = 0.0D0
               ENDDO
               DO I=1,KZYVAR/2
                  WRK(KORBWR+MJWOP(1,I,ISYM)-1) = 
     &                 WRK(KORBWR+MJWOP(1,I,ISYM)-1) +
     &                 WRK(KVEC-1+I)**2 + WRK(KVEC+KZYVAR/2-1+I)**2
               ENDDO
c     imag part
               DO I=1,MAXOCC
                  WRK(KORBWI+I-1) = 0.0D0
               ENDDO
               DO I=1,KZYVAR/2
                  WRK(KORBWI+MJWOP(1,I,ISYM)-1) = 
     &                 WRK(KORBWI+MJWOP(1,I,ISYM)-1) +
     &                 WRK(KVEC+KZYVAR-1+I)**2 + 
     &                 WRK(KVEC+3*KZYVAR/2-1+I)**2
               ENDDO
               NSRC=5
               CALL FINDMAXN(WRK(KORBWR),MAXOCC,ISRCORBR,NSRC)
               NSRC=5
               CALL FINDMAXN(WRK(KORBWI),MAXOCC,ISRCORBI,NSRC)
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
      SUBROUTINE ABSRESID(IOP,LABEL,NFREQ_ABS,FREQ_ABS,
     &     CONVERGED,REDVEC,IBTYP,FC,CMO,UDV,PV,XINDX,WRK,LWRK)
C
C PURPOSE:
C     Compute the complex residual vector to the linear response 
C     equation.
C
C     R = B[1] - ( E[2] - {w+iW}*S[2] )*(NR + iNI)
C
C     or equivalently the real and imaginary parts
C
C     RR = B[1] - E[2]*NR + S[2]*( w*NR - W*NI )
C     RI =      - E[2]*NI + S[2]*( w*NI + W*NR )
C
C     From file we read:
C         LURSP3 contains trial vectors
C         LURSP5 contains sigma vectors
C
C                          / b1 \                    / b2 \
C     NR = sum(i=1,KZRED) |      | * REDVEC(4i-3) + |      | * REDVEC(4i-1)
C                          \ b2 /_i                  \ b1 /_i           
C
C                          / b1 \                    / b2 \
C     NI = sum(i=1,KZRED) |      | * REDVEC(4i-2) + |      | * REDVEC(4i)
C                          \ b2 /_i                  \ b1 /_i           
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <absorp.h>
#include <inftap.h>
#include <ibndxdef.h>
C
      LOGICAL CONVERGED
      CHARACTER*8 LABEL
      DIMENSION REDVEC(2*MAXRM,MXFREQ),IBTYP(MAXRM),FREQ_ABS(MXFREQ),
     &     CMO(*),UDV(*),PV(*),FC(*),XINDX(*),WRK(LWRK)
C
      KFREE = 1
      LFREE = LWRK
      CALL MEMGET('REAL',KRR,KZYVAR,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KRI,KZYVAR,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KBVEC,KZYVAR,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KSLI,KZYVAR,WRK,KFREE,LFREE)
C
      CONVERGED = .TRUE.
C
C     Get the gradient
C
      DO IFREQ=1,NFREQ_ABS
C
C     Initialize residual vector
C
         CALL GETGPV(LABEL,FC,DUMMY,CMO,UDV,PV,XINDX,ANTSYM,
     &        WRK(KFREE),LFREE)
         CALL DCOPY(KZYVAR,WRK(KFREE),1,WRK(KRR),1)
         CALL DZERO(WRK(KRI),KZYVAR)
C
         IF (IPRABS.GE.5) THEN
            WRITE(LUPRI,'(2(/,5X,A))') ' B[1] gradient',
     &           '========================'
            CALL OUTPUT(WRK(KRR),1,KZVAR,1,2,KZVAR,2,1,LUPRI)
         END IF
C
C     Read trial and sigma vectors and add contributions to residual vector
C
         CALL REWSPL(LURSP3)
         KRSP3 = 0 
         CALL REWSPL(LURSP5)
         KRSP5 = 0
         IF (KOFFTY.EQ.1) THEN
#if defined (VAR_SPLITFILES)
            CALL READST(LURSP3,IREC,DUMMY,KRSP3,.FALSE.)
#else        
            READ (LURSP3)
#endif
#if defined (VAR_SPLITFILES)
            CALL READST(LURSP5,IREC,DUMMY,KRSP5,.FALSE.)
#else        
            READ (LURSP5)
#endif
         END IF
         DO I = 1, KZRED
C
            FR1 = REDVEC(4*I-3,IFREQ)
            FR2 = REDVEC(4*I-1,IFREQ)
            FI1 = REDVEC(4*I-2,IFREQ)
            FI2 = REDVEC(4*I,IFREQ)
C
            FR3 = (FREQ_ABS(IFREQ)*FR1-DAMPING*FI1)
            FR4 = (FREQ_ABS(IFREQ)*FR2-DAMPING*FI2)
            FI3 = (FREQ_ABS(IFREQ)*FI1+DAMPING*FR1)
            FI4 = (FREQ_ABS(IFREQ)*FI2+DAMPING*FR2)
C
C     Sigma vectors used that equal E[2]*N
C
#if defined (VAR_SPLITFILES)
            CALL READST(LURSP5,IREC,WRK(KBVEC),KRSP5,.TRUE.)
            IF (IREC.NE.KZYVAR) CALL QUIT('Read error in ABSRESID')
#else        
            CALL READT(LURSP5,KZYVAR,WRK(KBVEC))
#endif
C
            CALL DAXPY(KZYVAR,-FR1,WRK(KBVEC),1,WRK(KRR),1)
            CALL DAXPY(KZVAR,-FR2,WRK(KBVEC+KZVAR),1,WRK(KRR),1)
            CALL DAXPY(KZVAR,-FR2,WRK(KBVEC),1,WRK(KRR+KZVAR),1)
            CALL DAXPY(KZYVAR,-FI1,WRK(KBVEC),1,WRK(KRI),1)
            CALL DAXPY(KZVAR,-FI2,WRK(KBVEC+KZVAR),1,WRK(KRI),1)
            CALL DAXPY(KZVAR,-FI2,WRK(KBVEC),1,WRK(KRI+KZVAR),1)
C
C     Trial vectors used to perform S[2]*N
C
            IF (IBTYP(KOFFTY+I).EQ.JBCNDX) THEN
#if defined (VAR_SPLITFILES)
               CALL READST(LURSP3,IREC,WRK(KBVEC),KRSP3,.TRUE.)
               IF (IREC.NE.KZCONF) CALL QUIT('Read error in ABSRESID')
#else        
               CALL READT(LURSP3,KZCONF,WRK(KBVEC))
#endif
               CALL DZERO(WRK(KSLI),KZYVAR)
               CALL RSPSLI(1,0,WRK(KBVEC),DUMMY,
     &                     UDV,WRK(KSLI),XINDX,WRK(KFREE),LFREE)
            ELSE
#if defined (VAR_SPLITFILES)
               CALL READST(LURSP3,IREC,WRK(KBVEC),KRSP3,.TRUE.)
               IF (IREC.NE.KZYWOP) CALL QUIT('Read error in ABSRESID')
#else        
               CALL READT(LURSP3,KZYWOP,WRK(KBVEC))
#endif
               CALL DZERO(WRK(KSLI),KZYVAR)
               CALL RSPSLI(0,1,DUMMY,WRK(KBVEC),
     &                     UDV,WRK(KSLI),XINDX,WRK(KFREE),LFREE)
            END IF 
C
            CALL DAXPY(KZYVAR,FR3,WRK(KSLI),1,WRK(KRR),1)
            CALL DAXPY(KZVAR,-FR4,WRK(KSLI+KZVAR),1,WRK(KRR),1)
            CALL DAXPY(KZVAR,-FR4,WRK(KSLI),1,WRK(KRR+KZVAR),1)
            CALL DAXPY(KZYVAR,FI3,WRK(KSLI),1,WRK(KRI),1)
            CALL DAXPY(KZVAR,-FI4,WRK(KSLI+KZVAR),1,WRK(KRI),1)
            CALL DAXPY(KZVAR,-FI4,WRK(KSLI),1,WRK(KRI+KZVAR),1)
C
C     End loop over trial and sigma vectors
C
         END DO
C
         DNORM_RR = DNRM2(KZYVAR,WRK(KRR),1)
         DNORM_RI = DNRM2(KZYVAR,WRK(KRI),1)
         DNORM_RT = SQRT(DNORM_RR**2 + DNORM_RI**2)
C
         DNORM_NR = DNRM2(KZYRED,REDVEC(1,IFREQ),2)
         DNORM_NI = DNRM2(KZYRED,REDVEC(2,IFREQ),2)
         DNORM_NT = SQRT( DNORM_NR**2 + DNORM_NI**2 )
C
         RESID(1,IFREQ,IOP,KSYMOP) = DNORM_RR/DNORM_NT
         RESID(2,IFREQ,IOP,KSYMOP) = DNORM_RI/DNORM_NT
         RESID(3,IFREQ,IOP,KSYMOP) = DNORM_RT/DNORM_NT
C
         IF (DNORM_RT/DNORM_NT.GT.THCLR_ABSORP) CONVERGED=.FALSE.
C
C     End loop over frequencies
C
      END DO
C
C     Print norm of residual vector
C
      IF (IPRABS.GE.2) THEN
         WRITE(LUPRI,'(/,1X,A,1P,D9.2,/,1X,A,/,1X,A,I4,A)')
     &        'Convergence of RSP solution vectors, threshold =',
     &        THCLR_ABSORP,
     &  '-------------------------------------------------------------',
     &        '(dimension of reduced space:',KZYRED,')'
         DO IFREQ = 1,NFREQ_ABS
            WRITE(LUPRI,'(1X,A,F7.4,3X,A,F10.6,3X,A,1P,3D9.2)')
     &           'Frequency:',FREQ_ABS(IFREQ),
     &           'Damping:',DAMPING,
     &           'Residual:',(RESID(I,IFREQ,IOP,KSYMOP),I=1,3)
         END DO
      END IF
C
      RETURN
      END
      SUBROUTINE ABSLR(IOP,LABEL,NFREQ_ABS,FREQ_ABS,
     &     REDE,REDS,REDZ,REDGD,REDZGD,IBTYP,
     &     REDVEC,CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <absorp.h>
#include <inftap.h>
#include <ibndxdef.h>
C
C PURPOSE:
C     Solve the coupled equations in this iteration
C
C      (i)    ( E[2] - w*S[2] ) NR = B[1] - W*S[2]*NI
C     (ii)    ( E[2] - w*S[2] ) NI =        W*S[2]*NR
C
      CHARACTER*8 LABEL,BLANK
      PARAMETER (BLANK='        ', D0=0.0D0,DSQRT2=1.4142135623731D0)
      DIMENSION REDE(MAXRM*MAXRM),REDS(MAXRM*MAXRM),REDZ(2*MAXRM*MAXRM),
     &          REDGD(MAXRM),REDZGD(2*MAXRM),REDVEC(2*MAXRM,MXFREQ),
     &          IBTYP(MAXRM),FREQ_ABS(MXFREQ)
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(LWRK)
C
      KFREE = 1
      LFREE = LWRK
      CALL MEMGET('REAL',KGD1,KZYVAR,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KGD2,KZYVAR,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KSLI,KZYVAR,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KBVEC,KZYVAR,WRK,KFREE,LFREE)
C
      DO IFREQ=1,NFREQ_ABS
C
         IF (RESID(3,IFREQ,IOP,KSYMOP).LT.THCLR_ABSORP) THEN
            IF (IPRABS.GE.2) WRITE(LUPRI,'(/,1X,A,F8.4,A)') 
     &           '>>>>> Frequency:',FREQ_ABS(IFREQ),' converged <<<<<'
            GO TO 100
         END IF
C
C     Construct the right-hand-side of equations (i) and (ii)
C
         CALL GETGPV(LABEL,FC,DUMMY,CMO,UDV,PV,XINDX,ANTSYM,
     &        WRK(KFREE),LFREE)
         CALL DCOPY(KZYVAR,WRK(KFREE),1,WRK(KGD1),1)
         CALL DZERO(WRK(KGD2),KZYVAR)
C
C     Read trial vectors and add contributions to gradient
C
         CALL REWSPL(LURSP3)
         KRSP3 = 0
         IF (KOFFTY.EQ.1) THEN
#if defined (VAR_SPLITFILES)
            CALL READST(LURSP3,IREC,DUMMY,KRSP3,.FALSE.)
            IF (IREC.NE.KZCONF) CALL QUIT('Read error in ABSLR')
#else        
            READ (LURSP3)
#endif
         END IF
         DO I = 1, KZRED
            FR1 = DAMPING * REDVEC(4*I-3,IFREQ)
            FR2 = DAMPING * REDVEC(4*I-1,IFREQ)
            FI1 = DAMPING * REDVEC(4*I-2,IFREQ)
            FI2 = DAMPING * REDVEC(4*I,IFREQ)
            IF (IBTYP(KOFFTY+I).EQ.JBCNDX) THEN
#if defined (VAR_SPLITFILES)
               CALL READST(LURSP3,IREC,WRK(KBVEC),KRSP3,.TRUE.)
               IF (IREC.NE.KZCONF) CALL QUIT('Read error in ABSLR')
#else        
               CALL READT(LURSP3,KZCONF,WRK(KBVEC))
#endif
               CALL DZERO(WRK(KSLI),KZYVAR)
               CALL RSPSLI(1,0,WRK(KBVEC),DUMMY,
     &                     UDV,WRK(KSLI),XINDX,WRK(KFREE),LFREE)
               CALL DAXPY(KZCONF,-FI1,WRK(KSLI),1,WRK(KGD1),1)
               CALL DAXPY(KZCONF,FI2,WRK(KSLI),1,WRK(KGD1+KZVAR),1)
               CALL DAXPY(KZCONF,FR1,WRK(KSLI),1,WRK(KGD2),1)
               CALL DAXPY(KZCONF,-FR2,WRK(KSLI),1,WRK(KGD2+KZVAR),1)
            ELSE
#if defined (VAR_SPLITFILES)
               CALL READST(LURSP3,IREC,WRK(KBVEC),KRSP3,.TRUE.)
               IF (IREC.NE.KZYWOP) CALL QUIT('Read error in ABSLR')
#else        
               CALL READT(LURSP3,KZYWOP,WRK(KBVEC))
#endif
               CALL DZERO(WRK(KSLI),KZYVAR)
               CALL RSPSLI(0,1,DUMMY,WRK(KBVEC),
     &                     UDV,WRK(KSLI),XINDX,WRK(KFREE),LFREE)
               CALL DAXPY(KZWOPT,-FI1,WRK(KSLI),1,WRK(KGD1+KZCONF),1)
               CALL DAXPY(KZWOPT,-FI1,WRK(KSLI+KZWOPT),1,
     &                    WRK(KGD1+KZCONF+KZVAR),1)
               CALL DAXPY(KZWOPT,FI2,WRK(KSLI+KZWOPT),1,
     &                    WRK(KGD1+KZCONF),1)
               CALL DAXPY(KZWOPT,FI2,WRK(KSLI),1,
     &                    WRK(KGD1+KZCONF+KZVAR),1)
C
               CALL DAXPY(KZWOPT,FR1,WRK(KSLI),1,WRK(KGD2+KZCONF),1)
               CALL DAXPY(KZWOPT,FR1,WRK(KSLI+KZWOPT),1,
     &                    WRK(KGD2+KZCONF+KZVAR),1)
               CALL DAXPY(KZWOPT,-FR2,WRK(KSLI+KZWOPT),1,
     &                    WRK(KGD2+KZCONF),1)
               CALL DAXPY(KZWOPT,-FR2,WRK(KSLI),1,
     &                    WRK(KGD2+KZCONF+KZVAR),1)
            END IF 
         END DO
C
C     Get norm of solution vectors and gradients
C
         DNORM_NR = DNRM2(KZYRED,REDVEC(1,IFREQ),2)
         DNORM_NI = DNRM2(KZYRED,REDVEC(2,IFREQ),2)
         DNORM_NT = SQRT( DNORM_NR**2 + DNORM_NI**2 )
         DNORM_GD1 = DNRM2(KZYVAR,WRK(KGD1),1)
         DNORM_GD2 = DNRM2(KZYVAR,WRK(KGD2),1)
C
         IF (IPRABS.GE.2) WRITE(LUPRI,'(/,1X,A,F7.4,4(/,1X,A,1P,D9.2))')
     &        'Frequency = ',FREQ_ABS(IFREQ),
     &        'Norm of gradient vector 1 =',DNORM_GD1,
     &        'Norm of solution vector 1 =',DNORM_NR,
     &        'Norm of gradient vector 2 =',DNORM_GD2,
     &        'Norm of solution vector 2 =',DNORM_NI
C
C   Solve Eqs. (i) and (ii)
C
         MAXIT  = MAX_MACRO
         KEXCNV = 1
         KEXSTV = KEXCNV
         KEXSIM = KEXCNV
         RESTLR = .TRUE.
C
         IF (RESID(1,IFREQ,IOP,KSYMOP).GT.RESID(2,IFREQ,IOP,KSYMOP))THEN
C
         IF (RESID(1,IFREQ,IOP,KSYMOP).LT.THCLR_ABSORP/DSQRT2) THEN
            IF (IPRABS.GE.2) WRITE(LUPRI,'(/,1X,A,1P,D9.2,A)') 
     &           '>>>>> Skip real equation because residual is only:',
     &           RESID(1,IFREQ,IOP,KSYMOP), ' <<<<<'
         ELSE
            IF (IPRABS.GE.2) WRITE(LUPRI,'(/,/,1X,A,F8.4,A)') 
     &           '>>>>> Solving real equation for frequency:',
     &           FREQ_ABS(IFREQ),' <<<<<<'
            THCRSP = 0.1D0*THCLR_ABSORP*DNORM_NT/DNORM_NR/DSQRT2
            CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
     &           .TRUE.,BLANK,BLANK,WRK(KGD1),REDGD,REDE,REDS,
     &           IBTYP,FREQ_ABS(IFREQ),RESIDUAL,REDVEC,XINDX,
     &           WRK(KFREE),LFREE)
         END IF
C
         ELSE
C
         IF (RESID(2,IFREQ,IOP,KSYMOP).LT.THCLR_ABSORP/DSQRT2) THEN
            IF (IPRABS.GE.2) WRITE(LUPRI,'(/,1X,A,1P,D9.2,A)') 
     &        '>>>>> Skip imaginary equation because residual is only:',
     &           RESID(2,IFREQ,IOP,KSYMOP), ' <<<<<'
         ELSE
            IF (IPRABS.GE.2) WRITE(LUPRI,'(/,/,1X,A,F8.4,A)') 
     &           '>>>>> Solving imaginary equation for frequency:',
     &           FREQ_ABS(IFREQ),' <<<<<<'
            THCRSP = 0.5D0*THCLR_ABSORP*DNORM_NT/DNORM_NI/DSQRT2
            CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
     &           .TRUE.,BLANK,BLANK,WRK(KGD2),REDGD,REDE,REDS,
     &           IBTYP,FREQ_ABS(IFREQ),RESIDUAL,REDVEC,XINDX,
     &           WRK(KFREE),LFREE)
         END IF
C
         END IF
C
C     End of loop over frequencies
C
 100  CONTINUE
      END DO
C
      RETURN
      END
      SUBROUTINE ABSWRT(LUIN,LUOUT,LABEL,IOP,NFREQ_ABS,FREQ_ABS,
     &     IBTYP,REDVEC,WRK,LWRK)
C
C PURPOSE:
C
C     If LUIN=LURSP3, construct and write response vectors N to file
C     If LUIN=LURSP5, construct and write vectors E[2]*N to file
C
C     From file we read:
C     LUIN contains trial (LUIN=LURSP3) vectors or 
C     sigma vectors (LUIN=LURSP5). Output vectors are written to LUOUT
C     which normally is LURSP, apart from when we perform a reduction of
C     the dimension of the reduced space.
C
C                          / b1 \                    / b2 \
C     NR = sum(i=1,KZRED) |      | * REDVEC(4i-3) + |      | * REDVEC(4i-1)
C                          \ b2 /_i                  \ b1 /_i           
C
C                          / b1 \                    / b2 \
C     NI = sum(i=1,KZRED) |      | * REDVEC(4i-2) + |      | * REDVEC(4i)
C                          \ b2 /_i                  \ b1 /_i           
C
#include <implicit.h>
#include <priunit.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <absorp.h>
#include <ibndxdef.h>
C
      CHARACTER*8 LABEL,BLANK
      PARAMETER (D0 = 0.0D0, D1=1.0D0)
      DIMENSION REDVEC(2*MAXRM,MXFREQ),IBTYP(MAXRM),
     &     FREQ_ABS(MXFREQ),WRK(LWRK)
C
C     Do not allocate with MEMGET since that will give some extra
C     bytes of information in between real and imaginary parts
C
      KBVEC = 1
      KNR   = KBVEC + KZYVAR
      KNI   = KNR   + KZYVAR
      KFREE = KNI   + KZYVAR
      LFREE = LWRK  - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('ABSWRT',KFREE,LWRK)
C
      DO IFREQ=1,NFREQ_ABS
C
C     Read trial vectors and add contributions to response vector
C
         CALL REWSPL(LUIN)
         KRSP3 = 0
         IF (KOFFTY.EQ.1) THEN 
#if defined (VAR_SPLITFILES)
            CALL READST(LUIN,IREC,DUMMY,KRSP3,.FALSE.)
#else        
            READ (LUIN)
#endif
         END IF
C
         CALL DZERO(WRK(KNR),2*KZYVAR)
         DO I = 1, KZRED
            FR1 = REDVEC(4*I-3,IFREQ)
            FR2 = REDVEC(4*I-1,IFREQ)
            FI1 = REDVEC(4*I-2,IFREQ)
            FI2 = REDVEC(4*I,IFREQ)
C
            IF (IBTYP(KOFFTY+I).EQ.JBCNDX) THEN
#if defined (VAR_SPLITFILES)
               CALL READST(LUIN,IREC,WRK(KBVEC),KRSP3,.TRUE.)
               IF (IREC.NE.KZCONF) CALL QUIT('Read error in ABSWRT')
#else        
               CALL READT(LUIN,KZCONF,WRK(KBVEC))
#endif
               CALL DAXPY(KZCONF,FR1,WRK(KBVEC),1,WRK(KNR),1)
               CALL DAXPY(KZCONF,FR2,WRK(KBVEC),1,WRK(KNR+KZVAR),1)
               CALL DAXPY(KZCONF,FI1,WRK(KBVEC),1,WRK(KNI),1)
               CALL DAXPY(KZCONF,FI2,WRK(KBVEC),1,WRK(KNI+KZVAR),1)
            ELSE
#if defined (VAR_SPLITFILES)
               CALL READST(LUIN,IREC,WRK(KBVEC),KRSP3,.TRUE.)
               IF (IREC.NE.KZYWOP) CALL QUIT('Read error in ABSWRT')
#else        
               CALL READT(LUIN,KZYWOP,WRK(KBVEC))
#endif
C
               CALL DAXPY(KZWOPT,FR1,WRK(KBVEC),1,WRK(KNR+KZCONF),1)
               CALL DAXPY(KZWOPT,FR1,WRK(KBVEC+KZWOPT),1,
     &                    WRK(KNR+KZCONF+KZVAR),1)
               CALL DAXPY(KZWOPT,FR2,WRK(KBVEC+KZWOPT),1,
     &                    WRK(KNR+KZCONF),1)
               CALL DAXPY(KZWOPT,FR2,WRK(KBVEC),1,
     &                    WRK(KNR+KZCONF+KZVAR),1)
C
               CALL DAXPY(KZWOPT,FI1,WRK(KBVEC),1,WRK(KNI+KZCONF),1)
               CALL DAXPY(KZWOPT,FI1,WRK(KBVEC+KZWOPT),1,
     &                    WRK(KNI+KZCONF+KZVAR),1)
               CALL DAXPY(KZWOPT,FI2,WRK(KBVEC+KZWOPT),1,
     &                    WRK(KNI+KZCONF),1)
               CALL DAXPY(KZWOPT,FI2,WRK(KBVEC),1,
     &                    WRK(KNI+KZCONF+KZVAR),1)
            END IF 
         END DO
C
C     Write response vectors to file
C
         CALL WRTRSP(LUOUT,2*KZYVAR,WRK(KNR),LABEL,'COMPLEX ',
     &        FREQ_ABS(IFREQ),D0,KSYMOP,0,
     &        RESID(3,IFREQ,IOP,KSYMOP),D1)
C
         IF (IPRABS.GE.10) THEN
            WRITE(LUPRI,'(/A,A10,A,I4,/A,F12.8,A,I4)') 
     &           ' Response vector in ABSWRT for operator',LABEL,
     &           ' of symmetry',
     &           KSYMOP,' and frequency',FREQ_ABS(IFREQ),
     &           ' is written to file LU=',LUOUT
            CALL PRSYMB(LUPRI,'=',72,1)
            WRITE(LUPRI,'(4X,4A15)') 'ZR','YR','ZI','YI'
            CALL OUTPUT(WRK(KNR),1,KZYVAR/2,1,4,KZYVAR/2,4,1,LUPRI)
         END IF
C
      END DO
C
      RETURN
      END
      SUBROUTINE ABSCTL(IOPER,LABEL,
     &     REDE,REDS,REDZ,REDGD,REDZGD,REDVEC,IBTYP,
     &     CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,XINDX,
     &     RESLRF,WRK,LWRK)
C
#include <implicit.h>
#include <priunit.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <absorp.h>
#include <inftap.h>
#include <ibndxdef.h>
C
C PURPOSE:
C
C     Solve the complex LR equation 
C
C     ( E[2] - {w+iW}*S[2] )* (NR + iNI) = B[1]
C
C     or equivalently the pair of real equations
C
C     ( E[2] - w*S[2] )* NR = B[1] - W*S[2]*NI
C     ( E[2] - w*S[2] )* NI = W*S[2]*NR
C
      LOGICAL CONVERGED
      CHARACTER LABEL*8
      DIMENSION FREQ_ABS(MXFREQ)
      DIMENSION REDE(MAXRM*MAXRM),REDS(MAXRM*MAXRM),REDZ(2*MAXRM*MAXRM),
     &     REDGD(MAXRM),REDZGD(2*MAXRM),REDVEC(2*MAXRM,MXFREQ),
     &     IBTYP(MAXRM)
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),RESLRF(2,NFREQ_INTERVAL,3,3,2),WRK(LWRK)
C
      KFREE = 1
      LFREE = LWRK
C
C     Copy the frequencies into local variables to be used in the solver.
C     One can thereby choose to change the frequencies for certain operators.
C
      IF (LABEL(2:8).EQ.'ANGMOM') THEN
         NFREQ_ABS = 1
         FREQ_ABS(1) = 0.0D0
      ELSE
         NFREQ_ABS = NFREQ_ALPHA
         CALL DCOPY(NFREQ_ALPHA,FREQ_ALPHA,1,FREQ_ABS,1)
      END IF
C
C     Check if this an absorption calculation over a specified frequency
C     interval in which case we divide the interval in batches of freqs.
C
C     NFREQ_BATCH = number of frequencies per batch
C
      IBATCH = 0
 10   CONTINUE
C
      IF (ABS_INTERVAL) THEN
         NFREQ_ABS=MIN(NFREQ_BATCH,NFREQ_INTERVAL-IBATCH*NFREQ_BATCH)
         DO I=1,NFREQ_ABS
            FREQ_ABS(I)=FREQ_INTERVAL(1) + 
     &           (I-1+IBATCH*NFREQ_BATCH)*FREQ_INTERVAL(3)
         END DO
      END IF
C
      IF (IPRABS.GE.0) THEN
         WRITE(LUPRI,'(/2A,/3A,I2,A/,A,5(4F12.8,/,11X))') 
     &        ' ABSVEC1 -- Solving linear response equations',
     &        ' ( E[2] - {w-iW}*S[2] ) N = B[1]  ',
     &        ' ABSVEC1 -- for operator ', LABEL,
     &        ' of symmetry',KSYMOP,' at frequencies:',
     &        ' ABSVEC1 --',(FREQ_ABS(I),I=1,NFREQ_ABS)
      END IF
C
C     Check if response solution vectors already exist on file unit
C     LURSP with name RSPVEC.
C
      CALL CHKONFILE(LURSP,CONVERGED,LABEL,KSYMOP,
     &     NFREQ_ABS,FREQ_ABS,THCLR_ABSORP,WRK(KFREE))
      IF (CONVERGED) THEN
         WRITE(LUPRI,'(/A,I3)') 
     &        ' ABSCTL: converged response vectors found on file'//
     &        ' RSPVEC with unit=',LURSP
         CALL GETLRF(LURSP,LABEL,NFREQ_ABS,FREQ_ABS,
     &     FC,FV,CMO,UDV,PV,XINDX,IBATCH,RESLRF,
     &     WRK(KFREE),LFREE)
         GOTO 990
      END IF
C
C     Reduce the dimension of the reduced space
C
      IF (ABS_REDUCE) THEN
         CALL ABSREDUC(LABEL,NFREQ_ABS,FREQ_ABS,
     &        REDE,REDS,REDZ,REDGD,REDZGD,IBTYP,
     &        REDVEC,FC,FV,CMO,UDV,PV,XINDX,IBATCH,RESLRF,
     &        WRK(KFREE),LFREE)
         CALL REDSPACE(LABEL,IOPER,NFREQ_ABS,FREQ_ABS,
     &        REDE,REDS,REDVEC,IBTYP,
     &        FC,FV,CMO,UDV,PV,XINDX,WRK(KFREE),LFREE)
      END IF
C
      ITER = 0
 100  CONTINUE
C     
      IF (IPRABS.GE.2) THEN
         WRITE(LUPRI,'(/A,/A,I3,A,/A)')
     &        ' =======================================',
     &        ' >>>>  Macro iteration number',ITER,'  <<<<<<',
     &        ' ======================================='
C     
      END IF
C     
C     Solve for resonse vector in reduced space
C
      CALL ABSREDUC(LABEL,NFREQ_ABS,FREQ_ABS,
     &     REDE,REDS,REDZ,REDGD,REDZGD,IBTYP,
     &     REDVEC,FC,FV,CMO,UDV,PV,XINDX,IBATCH,RESLRF,
     &     WRK(KFREE),LFREE)
C     
C     Compute the residual in this iteration
C     
      CALL ABSRESID(IOPER,LABEL,NFREQ_ABS,FREQ_ABS,
     &     CONVERGED,REDVEC,IBTYP,FC,CMO,UDV,PV,XINDX,WRK(KFREE),LFREE)
C     
C     Check status
C     
      IF (CONVERGED) GOTO 900
      IF (ITER.GE.MAX_MACRO) GOTO 910
C     
C     Not converged, expand reduced space by solving coupled LR equations
C
      CALL ABSLR(IOPER,LABEL,NFREQ_ABS,FREQ_ABS,
     &     REDE,REDS,REDZ,REDGD,REDZGD,
     &     IBTYP,REDVEC,CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,XINDX,
     &     WRK(KFREE),LFREE)
C     
C     Iteration completed
C     
      ITER = ITER + 1
      GOTO 100
C     
C     Final messages from ABSCTL
C     
 900  CONTINUE
      IF (IPRABS.GE.0) THEN
         WRITE(LUPRI,'(/,1X,A,I4,A)')
     &        '*** ABSCTL: THE REQUESTED',NFREQ_ABS,
     &        ' SOLUTION VECTORS CONVERGED'
      END IF
      GOTO 920
 910  CONTINUE
      IF (IPRABS.GE.0) THEN
         WRITE(LUPRI,'(/,1X,A,I4,A)')
     &        '*** ABSCTL: THE REQUESTED',NFREQ_ABS,
     &        ' SOLUTION VECTORS NOT CONVERGED'
         WRITE(LUPRI,'(/,1X,A)')
     &        '--- MAXIMUM NUMBER OF ITERATIONS REACHED ---'
      END IF
      GOTO 920
C     
C     Construct and write response vectors to file
C     
 920  CONTINUE
      CALL ABSWRT(LURSP3,LURSP,LABEL,IOPER,NFREQ_ABS,FREQ_ABS,
     &     IBTYP,REDVEC,WRK(KFREE),LFREE)
      CALL GETLRF(LURSP,LABEL,NFREQ_ABS,FREQ_ABS,
     &     FC,FV,CMO,UDV,PV,XINDX,IBATCH,RESLRF,
     &     WRK(KFREE),LFREE)
C     
      IF (IPRABS.GE.0) THEN
         WRITE(LUPRI,'(2(/,1X,A),1P,D9.2,/,1X,A,/,1X,A,I4,A)')
     &        'RSP solution vectors written to file.',
     &        'Convergence of RSP solution vectors, threshold =',
     &        THCLR_ABSORP,
     & '-------------------------------------------------------------',
     &        '(dimension of reduced space:',KZYRED,')'
         DO IFREQ = 1,NFREQ_ABS
            WRITE(LUPRI,'(1X,A,F7.4,4X,A,F10.6,4X,A,1P,D9.2)')
     &           'Frequency:',FREQ_ABS(IFREQ),
     &           'Damping:',DAMPING,
     &           'Residual:',RESID(3,IFREQ,IOPER,KSYMOP)
         END DO
         WRITE(LUPRI,'(/2A,3(/A))') 
     &        '@ Value of linear response result for operator: ',
     &        LABEL,
     &        '@ --------------------------------------------',
     &        '@ Frequency   Damping         Real   Imaginary',
     &        '@ --------------------------------------------'
         IF (LABEL(2:8).EQ.'DIPLEN') THEN
            ITYPE = 1
         ELSE IF (LABEL(2:8).EQ.'ANGMOM') THEN
            ITYPE = 2
         ELSE
            WRITE(LUPRI,'(A)') ' Warning: Unknown operator!'
            GOTO 99
         END IF
         CALL DIPLAB(LABEL,INDEX)
         DO IFREQ=1,NFREQ_ABS
            WRITE(LUPRI,'(A,F10.4,F11.6,2F12.6)') '@',
     &         FREQ_ABS(IFREQ), DAMPING, 
     &         (RESLRF(I,IFREQ+IBATCH*NFREQ_BATCH,
     &                 INDEX,INDEX,ITYPE),I=1,2)
         END DO
 99      CONTINUE
      END IF
C
C     If this an absorption calculation over a specified frequency
C     interval then we need to do another batch of freqs. 
C
 990  CONTINUE
      IBATCH = IBATCH + 1
      IF (ABS_INTERVAL .AND. 
     &     NFREQ_INTERVAL-IBATCH*NFREQ_BATCH.GT.0) THEN
         GOTO 10
      END IF
C
      RETURN
      END
      SUBROUTINE GETLRF(LU,LABEL,NFREQ_ABS,FREQ_ABS,
     &     FC,FV,CMO,UDV,PV,XINDX,NBATCH,RESLRF,
     &     WRK,LWRK)
#include <implicit.h>
#include <priunit.h>
#include <wrkrsp.h>
#include <absorp.h>
C
      CHARACTER*8 LABEL,LABGD
      LOGICAL FOUND,CONV
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),XINDX(*),
     &     WRK(LWRK),RESLRF(2,NFREQ_INTERVAL,3,3,2),
     &     FREQ_ABS(MXFREQ)
C
      KFREE = 1
      LFREE = LWRK
C
      CALL DIPLAB(LABEL,INDEX)
      IF (LABEL(2:8).EQ.'DIPLEN') THEN
         ITYPE = 1
      ELSE IF (LABEL(2:8).EQ.'ANGMOM') THEN
         ITYPE = 2
      ELSE
         WRITE(LUPRI,'(A)') ' Warning: Unknown operator!'
         RETURN
      END IF
C
      CALL MEMGET('REAL',KVEC,2*KZYVAR,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KGD,KZYVAR,WRK,KFREE,LFREE)
C
      DO IOPER=1,NOPER(KSYMOP)
         LABGD = LABOP(IOPER,KSYMOP)
         IF (LABEL(2:8).EQ.LABGD(2:8)) THEN
            CALL DIPLAB(LABGD,INXGD)
            CALL GETGPV(LABGD,FC,FV,CMO,UDV,PV,XINDX,ANTSYM,
     &           WRK(KGD),LFREE)
            DO IFREQ=1,NFREQ_ABS
               CALL REARSP(LU,2*KZYVAR,WRK(KVEC),LABEL,'COMPLEX ',
     &              ABS(FREQ_ABS(IFREQ)),0.0D0,KSYMOP,0,THCLR_ABSORP,
     &              FOUND,CONV,ANTSYM)
               RESLRF(1,IFREQ+NBATCH*NFREQ_BATCH,INXGD,INDEX,ITYPE) = 
     &              DDOT(KZYVAR,WRK(KGD),1,WRK(KVEC),1)
               RESLRF(2,IFREQ+NBATCH*NFREQ_BATCH,INXGD,INDEX,ITYPE) = 
     &              DDOT(KZYVAR,WRK(KGD),1,WRK(KVEC+KZYVAR),1)
            END DO
         END IF
      END DO
C
      RETURN
      END
      SUBROUTINE CHKONFILE(LU,FOUND,LABEL,ISYM,NFREQ_ABS,FREQ_ABS,
     &     THD,FLAGS)
#include <priunit.h>
C
      LOGICAL     FOUND,FLAGS(NFREQ_ABS)
      CHARACTER*8 LABEL,LAB1,LAB2
      INTEGER     LU,ISYM,NFREQ_ABS,ISYM1,ISYM2,LEN,NBAS,NORB
      REAL*8      FREQ_ABS(NFREQ_ABS),THD,FREQ1,FREQ2,ANTSYM,RSD,EMCSCF
C
      FOUND = .TRUE.
      DO I=1,NFREQ_ABS
         FLAGS(I)=.FALSE.
      END DO
C
      REWIND(LU)
 100  READ(LU,END=900,ERR=900) LAB1,LAB2,FREQ1,FREQ2,ISYM1,ISYM2,
     &     ANTSYM,RSD,LEN,EMCSCF,NBAS,NORB
C
      IF (LAB1.NE.LABEL .OR. LAB2 .NE.'COMPLEX ' .OR.
     &   ISYM1.NE.ISYM  .OR. ISYM2.NE.0 .OR.
     &   RSD  .GT.THD   .OR. FREQ2.NE.0.0D0) GOTO 100
C
      DO I=1,NFREQ_ABS
         IF (FREQ1.EQ.FREQ_ABS(I)) THEN
            FLAGS(I)=.TRUE.
         END IF
      END DO
C
      READ(LU,END=900,ERR=900)
C      
      GOTO 100
C
 900  CONTINUE
      DO I=1,NFREQ_ABS
         FOUND = FOUND .AND. FLAGS(I)
      END DO
      RETURN
      END
      SUBROUTINE REDSPACE(LABEL,IOPER,NFREQ_ABS,FREQ_ABS,
     &     REDE,REDS,REDVEC,IBTYP,FC,FV,CMO,UDV,PV,XINDX,WRK,LWRK)
#include <implicit.h>
#include <priunit.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <absorp.h>
#include <inftap.h>
#include <ibndxdef.h>
C
      CHARACTER*8 LABEL
      LOGICAL FOUND,CONV,READFLAG
      DIMENSION REDVEC(2*MAXRM,MXFREQ),IBTYP(MAXRM),
     &     REDE(MAXRM*MAXRM),REDS(MAXRM*MAXRM),
     &     FREQ_ABS(MXFREQ),WRK(LWRK)
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),XINDX(*)
C
      IF (IPRABS.GE.2) WRITE(LUPRI,'(2(/A))')
     &     ' INFO: lowering dimension of reduced space'
      CALL FLSHFO(LUPRI)
C
C     Reduction of dimension in reduced space only implented for SCF
C     type wave functions.
C
      IF (KZCONF .GT. 0) RETURN
C
      KBVEC = 1
      KSVEC = KBVEC + KZYVAR
      KB    = KSVEC + KZYVAR
      KS    = KB    + 2*KZYVAR
      KFREE = KS    + 2*KZYVAR
      LFREE = LWRK  - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('ABSWRT',KFREE,LWRK)
C
C     Write response vectors N and sigma vectors E[2]*N to respective
C     temporary files LUBTMP and LUSTMP. Unit numbers are set in GPOPEN.
C
      LUBTMP=-1
      LUSTMP=-1
      CALL GPOPEN(LUBTMP,'RSPBTMP','NEW','SEQUENTIAL','UNFORMATTED',
     &     0,.FALSE.)
      CALL GPOPEN(LUSTMP,'RSPSTMP','NEW','SEQUENTIAL','UNFORMATTED',
     &     0,.FALSE.)
      WRITE(LUBTMP) 'EOFLABEL'
      WRITE(LUSTMP) 'EOFLABEL'
      CALL ABSWRT(LURSP3,LUBTMP,LABEL,IOPER,NFREQ_ABS,FREQ_ABS,
     &     IBTYP,REDVEC,WRK(KFREE),LFREE)
      CALL ABSWRT(LURSP5,LUSTMP,LABEL,IOPER,NFREQ_ABS,FREQ_ABS,
     &     IBTYP,REDVEC,WRK(KFREE),LFREE)
C
C     Read response and sigma vectors from temporary files and perform
C     modified Gram-Schmidt orthogonalization of the respective real and
C     imaginary parts of the response vectors (these then form the new
C     set of trial vectors which are written to LURSP3). The
C     corresponding sigma vectors are contructed and written to LURSP5.
C     We keep NBOLD of the old trial vectors in order to avoid linear
C     dependence.
C
      NBOLD = 0
      NLINDEP = 0
      IFREQ=0
      READFLAG = .TRUE.
      DO I=1+NBOLD,2*NFREQ_ABS+NBOLD
         IF (READFLAG) THEN
            IFREQ=IFREQ+1
            CALL REARSP(LUBTMP,LENGTH,WRK(KB),LABEL,'COMPLEX ',
     &           ABS(FREQ_ABS(IFREQ)),0.0D0,
     &           KSYMOP,0,THCLR,FOUND,CONV,ANTSYM)
            IF (LENGTH.NE.2*KZYVAR) CALL QUIT(' REDSPACE: read error')
            CALL REARSP(LUSTMP,LENGTH,WRK(KS),LABEL,'COMPLEX ',
     &           ABS(FREQ_ABS(IFREQ)),0.0D0,
     &           KSYMOP,0,THCLR,FOUND,CONV,ANTSYM)
            IF (LENGTH.NE.2*KZYVAR) CALL QUIT(' REDSPACE: read error')
            CALL DCOPY(KZYVAR,WRK(KB),1,WRK(KBVEC),1)
            CALL DCOPY(KZYVAR,WRK(KS),1,WRK(KSVEC),1)
            READFLAG = .FALSE.
C
            IF (IPRABS.GE.2) THEN
               CALL GETGPV(LABEL,FC,FV,CMO,UDV,PV,XINDX,ANTSYM,
     &              WRK(KFREE),LFREE)
               VAL1 = DDOT(KZYVAR,WRK(KFREE),1,WRK(KB),1)
               VAL2 = DDOT(KZYVAR,WRK(KFREE),1,WRK(KB+KZYVAR),1)
               WRITE(LUPRI,'(2A,F10.6,2F16.10)') 
     &              ' Value of linear response for ',
     &              LABEL,FREQ_ABS(IFREQ),VAL1,VAL2
            END IF
C
            IF (IPRABS.GE.10) THEN
               WRITE(LUPRI,'(/A,A10,A,I4,/A,F12.8,A,I4)') 
     &              ' Response vector in REDSPACE for operator',LABEL,
     &              ' of symmetry',
     &              KSYMOP,' and frequency',FREQ_ABS(IFREQ),
     &              ' is read from file LUBTMP=',LUBTMP
               CALL PRSYMB(LUPRI,'=',72,1)
               WRITE(LUPRI,'(4X,4A15)') 'ZR','YR','ZI','YI'
               CALL OUTPUT(WRK(KB),1,KZVAR,1,4,KZVAR,4,1,LUPRI)
               WRITE(LUPRI,'(/A,A10,A,I4,/A,F12.8,A,I4)') 
     &              ' Sigma vector in REDSPACE for operator',LABEL,
     &              ' of symmetry',
     &              KSYMOP,' and frequency',FREQ_ABS(IFREQ),
     &              ' is read from file LUSTMP=',LUSTMP
               CALL PRSYMB(LUPRI,'=',72,1)
               WRITE(LUPRI,'(4X,4A15)') 'ZR','YR','ZI','YI'
               CALL OUTPUT(WRK(KS),1,KZVAR,1,4,KZVAR,4,1,LUPRI)
            END IF
         ELSE
            CALL DCOPY(KZYVAR,WRK(KB+KZYVAR),1,WRK(KBVEC),1)
            CALL DCOPY(KZYVAR,WRK(KS+KZYVAR),1,WRK(KSVEC),1)
            READFLAG = .TRUE.
         END IF
C
C     First KZYVAR elements of the space for complex response vectors
C     are now available as scratch.
C
         KBTMP=KB
         KSTMP=KS
C
         KRSP3 = 0
         KRSP5 = 0
         CALL REWSPL(LURSP3)
         CALL REWSPL(LURSP5)
C     
         DO J=1,I-1-NLINDEP
#if defined (VAR_SPLITFILES)
            CALL READST(LURSP3,KZYVAR,WRK(KBTMP),KRSP3,.TRUE.)
            CALL READST(LURSP5,KZYVAR,WRK(KSTMP),KRSP5,.TRUE.)
#else        
            CALL READT(LURSP3,KZYVAR,WRK(KBTMP))
            CALL READT(LURSP5,KZYVAR,WRK(KSTMP))
#endif
C     Orthogonalize against trial vectors (Z,Y)
            FAC = DDOT(KZYVAR,WRK(KBTMP),1,WRK(KBVEC),1)
            CALL DAXPY(KZYVAR,-FAC,WRK(KBTMP),1,WRK(KBVEC),1)
            CALL DAXPY(KZYVAR,-FAC,WRK(KSTMP),1,WRK(KSVEC),1)
C     Orthogonalize against trial vectors (Y,Z)
            FAC = DDOT(KZVAR,WRK(KBTMP+KZVAR),1,WRK(KBVEC),1) +
     &            DDOT(KZVAR,WRK(KBTMP),1,WRK(KBVEC+KZVAR),1)
            CALL DAXPY(KZVAR,-FAC,WRK(KBTMP+KZVAR),1,WRK(KBVEC),1)
            CALL DAXPY(KZVAR,-FAC,WRK(KBTMP),1,WRK(KBVEC+KZVAR),1)
            CALL DAXPY(KZVAR,-FAC,WRK(KSTMP+KZVAR),1,WRK(KSVEC),1)
            CALL DAXPY(KZVAR,-FAC,WRK(KSTMP),1,WRK(KSVEC+KZVAR),1)
         END DO
C
C     Check linear dependence.
C
         BNORM = DNORM2(KZYVAR,WRK(KBVEC),1)
         IF (BNORM.LT.1.0D-6) THEN
            IF (IPRABS.GE.0) WRITE(LUPRI,'(A)') 
     &           ' INFO: vector not used as trial vector'//
     &           ' due to linear dependence'
            NLINDEP = NLINDEP + 1
         ELSE
            FAC = 1.0D0/BNORM
            CALL DSCAL(KZYVAR,FAC,WRK(KBVEC),1)
            CALL DSCAL(KZYVAR,FAC,WRK(KSVEC),1)
C     Orthogonalize trial vector (Z,Y) against (Y,Z) using
C     symemtric orthogonalization to keep paired structure
            OVLP=2.0D0*DDOT(KZVAR,WRK(KBVEC),1,WRK(KBVEC+KZVAR),1)
            C1=0.5D0*(1/SQRT(1+OVLP) + 1/SQRT(1-OVLP))
            C2=0.5D0*(1/SQRT(1+OVLP) - 1/SQRT(1-OVLP))
            CALL DCOPY(KZYVAR,WRK(KBVEC),1,WRK(KBTMP),1)
            CALL DCOPY(KZYVAR,WRK(KSVEC),1,WRK(KSTMP),1)
            CALL DSCAL(KZYVAR,C1,WRK(KBVEC),1)
            CALL DSCAL(KZYVAR,C1,WRK(KSVEC),1)
            CALL DAXPY(KZVAR,C2,WRK(KBTMP),1,WRK(KBVEC+KZVAR),1)
            CALL DAXPY(KZVAR,C2,WRK(KBTMP+KZVAR),1,WRK(KBVEC),1)
            CALL DAXPY(KZVAR,C2,WRK(KSTMP),1,WRK(KSVEC+KZVAR),1)
            CALL DAXPY(KZVAR,C2,WRK(KSTMP+KZVAR),1,WRK(KSVEC),1)
C
#if defined (VAR_SPLITFILES)
            CALL WRITST(LURSP3,KZYVAR,WRK(KBVEC),KRSP3)
            CALL WRITST(LURSP5,KZYVAR,WRK(KSVEC),KRSP5)
#else        
            CALL WRITT(LURSP3,KZYVAR,WRK(KBVEC))
            CALL WRITT(LURSP5,KZYVAR,WRK(KSVEC))
#endif
            IF (IPRABS.GE.10) THEN
               WRITE(LUPRI,'(/A,I4,A,I2)') 
     &              ' INFO: Normalized trial vector',I,
     &              ' is written to file LURSP3=',LURSP3
               IF (IPRABS.GE.0) THEN
                  CALL PRSYMB(LUPRI,'=',72,1)
                  WRITE(LUPRI,'(4X,2A15)') 'Z','Y'
                  CALL OUTPUT(WRK(KBVEC),1,KZVAR,1,2,
     &                 KZVAR,2,1,LUPRI)
                  CALL FLSHFO(LUPRI)
               END IF
            END IF
         END IF
C
      END DO
C
      CALL GPCLOSE(LUBTMP,'DELETE')
      CALL GPCLOSE(LUSTMP,'DELETE')
C
C     Construct the reduced E[2] and S[2] matrices in triangular packed
C     format.
C
      KZRED = 2*NFREQ_ABS+NBOLD-NLINDEP
      KZYRED = 2*KZRED
C
      DO I=1,KZRED
         IROW=2*I-1
         JOFF1 = IROW*(IROW-1)/2
         JOFF2 = JOFF1+IROW
         REWIND(LURSP3)
         DO IDUM=1,I
            CALL READT(LURSP3,KZYVAR,WRK(KBVEC))
         END DO            
         REWIND(LURSP3)
         REWIND(LURSP5)
         DO J=1,I
            CALL READT(LURSP3,KZYVAR,WRK(KB))
            CALL READT(LURSP5,KZYVAR,WRK(KS))
            ICOL=2*J-1
            REDE(JOFF1+ICOL) = 
     &           DDOT(KZYVAR,WRK(KBVEC),1,WRK(KS),1)
            IF (I.NE.J) REDE(JOFF1+ICOL+1) = 
     &              DDOT(KZVAR,WRK(KBVEC),1,WRK(KS+KZVAR),1) +
     &              DDOT(KZVAR,WRK(KBVEC+KZVAR),1,WRK(KS),1)
            REDE(JOFF2+ICOL) = 
     &           DDOT(KZVAR,WRK(KBVEC+KZVAR),1,WRK(KS),1) +
     &           DDOT(KZVAR,WRK(KBVEC),1,WRK(KS+KZVAR),1)
            REDE(JOFF2+ICOL+1) = 
     &           DDOT(KZVAR,WRK(KBVEC+KZVAR),1,WRK(KS+KZVAR),1) +
     &           DDOT(KZVAR,WRK(KBVEC),1,WRK(KS),1)
            CALL DSCAL(KZVAR,-1.0D0,WRK(KB+KZVAR),1)
            REDS(JOFF1+ICOL) = 
     &           DDOT(KZYVAR,WRK(KBVEC),1,WRK(KB),1)
            REDS(JOFF2+ICOL) = 
     &           DDOT(KZVAR,WRK(KBVEC+KZVAR),1,WRK(KB),1) +
     &           DDOT(KZVAR,WRK(KBVEC),1,WRK(KB+KZVAR),1)
            CALL DSCAL(KZYVAR,-1.0D0,WRK(KB),1)
            IF (I.NE.J) REDS(JOFF1+ICOL+1) = 
     &              DDOT(KZVAR,WRK(KBVEC),1,WRK(KB+KZVAR),1) +
     &              DDOT(KZVAR,WRK(KBVEC+KZVAR),1,WRK(KB),1)
            REDS(JOFF2+ICOL+1) = 
     &           DDOT(KZVAR,WRK(KBVEC+KZVAR),1,WRK(KB+KZVAR),1) +
     &           DDOT(KZVAR,WRK(KBVEC),1,WRK(KB),1)
         END DO
      END DO
      CALL DSCAL(KZYRED*(KZYRED+1)/2,2.0D0,REDS,1)
C
      IF (IPRABS.GE.10) THEN
         WRITE(LUPRI,'(/A)') 'Reduced E[2] after reduction = '
         CALL OUTPAK(REDE,KZYRED,1,LUPRI)
         WRITE(LUPRI,'(/A)') 'Reduced S[2] after reduction = '
         CALL OUTPAK(REDS,KZYRED,1,LUPRI)
         CALL FLSHFO(LUPRI)
      END IF
C
      RETURN
      END
      SUBROUTINE BETA_SETUP
#include <implicit.h>
#include <priunit.h>
#include <inforb.h>
#include <absorp.h>
C
      PARAMETER (THRZERO = 1.0D-6)
C
      LOGICAL DOHYP
      CHARACTER*8 ALAB,BLAB,CLAB
C
      NQRF = 0
      NFREQ_ALPHA = 0
C
      DO ICFR=1,NFREQ_BETA_C
      DO IBFR=1,NFREQ_BETA_B
         FREQB = FREQ_BETA_B(IBFR)
         FREQC = FREQ_BETA_C(ICFR)
         FREQA = -(FREQB+FREQC)
         IF (ABS_SHG .AND. .NOT.FREQB.EQ.FREQC) GOTO 100
         DO ISYMA=1,NSYM
         DO ISYMB=1,NSYM
            ISYMC = MULD2H(ISYMA,ISYMB)
            IF (NOPER(ISYMA).GT.0 .AND. NOPER(ISYMB).GT.0 .AND.
     &           NOPER(ISYMC).GT.0) THEN
               DO IAOP=1,NOPER(ISYMA)
               DO IBOP=1,NOPER(ISYMB)
               DO ICOP=1,NOPER(ISYMC)
                  ALAB = LABOP(IAOP,ISYMA)
                  BLAB = LABOP(IBOP,ISYMB)
                  CLAB = LABOP(ICOP,ISYMC)
                  CALL QRCHK(DOHYP,ALAB,BLAB,CLAB,
     &                 ISYMA,ISYMB,ISYMC,FREQA,FREQB,FREQC)
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
C
      IF (IPRABS.GE.0) THEN
         CALL AROUND('Setup of Hyperpolarizability Calculation')
         WRITE (LUPRI,'(2(/A),I4,A)')
     & ' This calculations requires the solution of linear response',
     & ' equations for electric dipole operators at',NFREQ_ALPHA,
     & ' frequencies:'
         WRITE(LUPRI,'(/A,5(4F12.8,/,9X))')
     &        ' LR FREQ:',(FREQ_ALPHA(I),I=1,NFREQ_ALPHA)
         WRITE (LUPRI,'(/A,I4,A)')
     & ' and the evaluation of',NQRF,' quadratic response functions:'
         WRITE(LUPRI,'(/2A,/A3,6A12,/2A)') 
     &        '--------------------------------------',
     &        '-------------------------------------',
     &        ' No','A-oper','B-oper','C-oper',
     &        'A-freq','B-freq','C-freq',
     &        '--------------------------------------',
     &        '-------------------------------------'
         DO IQRF=1,NQRF
            WRITE(LUPRI,'(I3,3A12,3F12.8)') IQRF,
     &           (QRFLAB(IQRF,I), I=1,3),(QRFFRQ(IQRF,I), I=1,3)
         END DO
         WRITE(LUPRI,'(2A)') 
     &        '--------------------------------------',
     &        '-------------------------------------'
      END IF
C
C     End of BETA_SETUP
C
      RETURN
      END
      SUBROUTINE GPDSRT(N,X,THD)
#include <implicit.h>
C
C========================================================================
C
C Purpose: Sort the elements of a vector X containing N double-precision
C          real numbers. Remove repeated occurrences of elements which
C          are seperated by less than THD.
C
C     In: X   - vector of double-precision real numbers
C         N   - number of elements in vector
C         THD - threshold for two elements being considered equal
C
C    Out: X - vector containing unique elements, now sorted
C         N - number of unique elements
C
C Author: Patrick Norman, 2003.
C
C========================================================================
C
      DIMENSION X(N)
C
      IF (N.LE.1) RETURN
C
      M = N
      J = 1
C
 100  CONTINUE
      I = J+1
 200  CONTINUE
C
      IF (ABS(X(J)-X(I)).LT.THD) THEN
         X(I) = X(M)
         M = M-1
         IF (I.LE.M) GOTO 200
      END IF
      IF (X(I).LT.X(J)) THEN
         T = X(J)
         X(J) = X(I)
         X(I) = T
      END IF
C
      I = I+1
      IF (I.LE.M) GOTO 200
C
      J = J+1
      IF (J.LT.M ) GOTO 100
C
      N = M
C
      RETURN
      END
      SUBROUTINE ABSQRF(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     &                  XINDX,MJWOP,WRK,LWRK)
C
C PURPOSE:
C     Compute the quadratice response functions with response vectors
C     found on file.
C
#include <implicit.h>
#include <priunit.h>
#include <absorp.h>
#include <wrkrsp.h>
#include <infvar.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <symmet.h>
#include <cbilrs.h>
#include <infrsp.h>
#include <inforb.h>
#include <qrinf.h>
#include <infdim.h>
C
      PARAMETER ( D0=0.0D0 )
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),MJWOP(2,MAXWOP,8),WRK(LWRK)
C
      DIMENSION HYPVAL(2),TMPVAL(2),VAL(4)
      CHARACTER*8 ALAB,BLAB,CLAB
C
      KFREE = 1
      LFREE = LWRK
C
      ISPINA = 0
      ISPINB = 0
      ISPINC = 0
C
      MZYVMX = 2*NVARMA
      CALL MEMGET('REAL',KVECA,2*MZYVMX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KVECB,2*MZYVMX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KVECC,2*MZYVMX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KGRAD,2*MZYVMX,WRK,KFREE,LFREE)
C
      DO IQRF=1,NQRF
C
         HYPVAL(1) = D0
         HYPVAL(2) = D0
C
         ALAB  = QRFLAB(IQRF,1)
         BLAB  = QRFLAB(IQRF,2)
         CLAB  = QRFLAB(IQRF,3)
         ISYMA = QRFSYM(IQRF,1)
         ISYMB = QRFSYM(IQRF,2)
         ISYMC = QRFSYM(IQRF,3)
         FREQA = QRFFRQ(IQRF,1)
         FREQB = QRFFRQ(IQRF,2)
         FREQC = QRFFRQ(IQRF,3)
C
         WRITE(LUPRI,'(/,2(A,I4),A,3(/A,A10,I4,F10.6))')
     &        ' Quadratic response function no',IQRF,' of',NQRF,'.',
     &        ' A operator, symmetry, frequency: ',ALAB,ISYMA,FREQA,
     &        ' B operator, symmetry, frequency: ',BLAB,ISYMB,FREQB,
     &        ' C operator, symmetry, frequency: ',CLAB,ISYMC,FREQC
C
         CALL QRRDVE(ISYMA,ISYMB,ISYMC,ALAB,BLAB,CLAB,
     &        FREQA,FREQB,FREQC,KZYVA,KZYVB,KZYVC,WRK(KVECA),WRK(KVECB),
     &        WRK(KVECC))
C
C     Calculate Na B[2] Nc type terms (two permutations)
C
         CALL DZERO(WRK(KGRAD),2*MZYVMX)
         CALL X2INIT(1,KZYVA,KZYVC,ISYMA,ISPINA,ISYMC,ISPINC,WRK(KVECA),
     &        WRK(KVECC),WRK(KGRAD),XINDX,UDV,PV,BLAB,ISYMB,ISPINB,
     &        CMO,MJWOP,WRK(KFREE),LFREE)
         VAL(1) = DDOT(KZYVA,WRK(KGRAD),1,WRK(KVECA),1)
         VAL(2) = DDOT(KZYVA,WRK(KGRAD),1,WRK(KVECA+KZYVA),1)
         TMPVAL(1) = VAL(1)
         TMPVAL(2) = VAL(2)
         HYPVAL(1) = HYPVAL(1) + VAL(1)      
         HYPVAL(2) = HYPVAL(2) + VAL(2)
C
         CALL DZERO(WRK(KGRAD),2*MZYVMX)
         CALL X2INIT(1,KZYVA,KZYVB,ISYMA,ISPINA,ISYMB,ISPINB,WRK(KVECA),
     &        WRK(KVECB),WRK(KGRAD),XINDX,UDV,PV,CLAB,ISYMC,ISPINC,
     &        CMO,MJWOP,WRK(KFREE),LFREE)
         VAL(1) = DDOT(KZYVA,WRK(KGRAD),1,WRK(KVECA),1)
         VAL(2) = DDOT(KZYVA,WRK(KGRAD),1,WRK(KVECA+KZYVA),1)
         TMPVAL(1) = TMPVAL(1) + VAL(1)
         TMPVAL(2) = TMPVAL(2) + VAL(2)
         HYPVAL(1) = HYPVAL(1) + VAL(1)      
         HYPVAL(2) = HYPVAL(2) + VAL(2)      
C
         IF (DAMPING.GT.D0) THEN
            CALL DZERO(WRK(KGRAD),2*MZYVMX)
            CALL X2INIT(1,KZYVA,KZYVC,ISYMA,ISPINA,ISYMC,ISPINC,
     &           WRK(KVECA),WRK(KVECC+KZYVC),WRK(KGRAD),XINDX,UDV,PV,
     &           BLAB,ISYMB,ISPINB,CMO,MJWOP,WRK(KFREE),LFREE)
            VAL(1) = -DDOT(KZYVA,WRK(KGRAD),1,WRK(KVECA+KZYVA),1)
            VAL(2) = DDOT(KZYVA,WRK(KGRAD),1,WRK(KVECA),1)
            TMPVAL(1) = TMPVAL(1) + VAL(1)
            TMPVAL(2) = TMPVAL(2) + VAL(2)
            HYPVAL(1) = HYPVAL(1) + VAL(1)      
            HYPVAL(2) = HYPVAL(2) + VAL(2)
C     
            CALL DZERO(WRK(KGRAD),2*MZYVMX)
            CALL X2INIT(1,KZYVA,KZYVB,ISYMA,ISPINA,ISYMB,ISPINB,
     &           WRK(KVECA),WRK(KVECB+KZYVB),WRK(KGRAD),XINDX,UDV,PV,
     &           CLAB,ISYMC,ISPINC,CMO,MJWOP,WRK(KFREE),LFREE)
            VAL(1) = -DDOT(KZYVA,WRK(KGRAD),1,WRK(KVECA+KZYVA),1)
            VAL(2) = DDOT(KZYVA,WRK(KGRAD),1,WRK(KVECA),1)
            TMPVAL(1) = TMPVAL(1) + VAL(1)
            TMPVAL(2) = TMPVAL(2) + VAL(2)
            HYPVAL(1) = HYPVAL(1) + VAL(1)      
            HYPVAL(2) = HYPVAL(2) + VAL(2)
         END IF
C
         IF (IPRABS.GE.0) WRITE(LUPRI,'(A15,4F15.6)') 
     &        ' Na X[2] Ny    ',TMPVAL(1),TMPVAL(2),HYPVAL(1),HYPVAL(2)
C
C     Calculate Nb A[2] Nc type terms (two permutations)
C
         CALL DZERO(WRK(KGRAD),2*MZYVMX)
         CALL A2INIT(1,KZYVB,KZYVC,ISYMB,ISPINB,ISYMC,ISPINC,1,
     &        WRK(KVECC),WRK(KGRAD),XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     &        CMO,MJWOP,WRK(KFREE),LFREE)
         VAL(1) = DDOT(KZYVB,WRK(KGRAD),1,WRK(KVECB),1)
         VAL(2) = DDOT(KZYVB,WRK(KGRAD),1,WRK(KVECB+KZYVB),1)
         TMPVAL(1) = VAL(1)
         TMPVAL(2) = VAL(2)
         HYPVAL(1) = HYPVAL(1) + VAL(1)      
         HYPVAL(2) = HYPVAL(2) + VAL(2)
C
         CALL DZERO(WRK(KGRAD),2*MZYVMX)
         CALL A2INIT(1,KZYVC,KZYVB,ISYMC,ISPINC,ISYMB,ISPINB,1,
     &        WRK(KVECB),WRK(KGRAD),XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     &        CMO,MJWOP,WRK(KFREE),LFREE)
         VAL(1) = DDOT(KZYVC,WRK(KGRAD),1,WRK(KVECC),1)
         VAL(2) = DDOT(KZYVC,WRK(KGRAD),1,WRK(KVECC+KZYVC),1)
         TMPVAL(1) = TMPVAL(1) + VAL(1)
         TMPVAL(2) = TMPVAL(2) + VAL(2)
         HYPVAL(1) = HYPVAL(1) + VAL(1)      
         HYPVAL(2) = HYPVAL(2) + VAL(2)      
C
         IF (DAMPING.GT.D0) THEN
            CALL DZERO(WRK(KGRAD),2*MZYVMX)
            CALL A2INIT(1,KZYVB,KZYVC,ISYMB,ISPINB,ISYMC,ISPINC,1,
     &           WRK(KVECC+KZYVC),WRK(KGRAD),XINDX,UDV,PV,ALAB,ISYMA,
     &           ISPINA,CMO,MJWOP,WRK(KFREE),LFREE)
            VAL(1) = -DDOT(KZYVB,WRK(KGRAD),1,WRK(KVECB+KZYVB),1)
            VAL(2) = DDOT(KZYVB,WRK(KGRAD),1,WRK(KVECB),1)
            TMPVAL(1) = TMPVAL(1) + VAL(1)
            TMPVAL(2) = TMPVAL(2) + VAL(2)
            HYPVAL(1) = HYPVAL(1) + VAL(1)      
            HYPVAL(2) = HYPVAL(2) + VAL(2)
C     
            CALL DZERO(WRK(KGRAD),2*MZYVMX)
            CALL A2INIT(1,KZYVC,KZYVB,ISYMC,ISPINC,ISYMB,ISPINB,1,
     &           WRK(KVECB+KZYVB),WRK(KGRAD),XINDX,UDV,PV,ALAB,ISYMA,
     &           ISPINA,CMO,MJWOP,WRK(KFREE),LFREE)
            VAL(1) = -DDOT(KZYVC,WRK(KGRAD),1,WRK(KVECC+KZYVC),1)*
     &           DAMPING*DAMPING
            VAL(2) = DDOT(KZYVC,WRK(KGRAD),1,WRK(KVECC),1)
            TMPVAL(1) = TMPVAL(1) + VAL(1)
            TMPVAL(2) = TMPVAL(2) + VAL(2)
            HYPVAL(1) = HYPVAL(1) + VAL(1)      
            HYPVAL(2) = HYPVAL(2) + VAL(2)
         END IF      
C
         IF (IPRABS.GE.0) WRITE(LUPRI,'(A15,4F15.6)') 
     &        ' Nx A[2] Ny    ',TMPVAL(1),TMPVAL(2),HYPVAL(1),HYPVAL(2)
C
C     Calculate Na*(E[3]-w*S[3]-i*W*R[3])*Nb*Nc type terms (two permutations)
C
         CALL T3DRV(1,ISYMA,ISYMB,ISYMC,WRK(KVECB),WRK(KVECC),
     &        .FALSE.,WRK(KVECA),-FREQB,-FREQC,XINDX,UDV,PV,MJWOP,
     &        WRK(KFREE),LFREE,CMO,FC,FV)
         IF (IPRABS.GE.10) THEN
            WRITE(LUPRI,'(/A)') ' E[3] vector '
            CALL PRSYMB(LUPRI,'=',72,1)
            WRITE(LUPRI,'(4X,4A15)') 'ZR','YR','ZI','YI'
            CALL OUTPUT(WRK(KFREE),1,KZYVA/2,1,4,KZYVA/2,4,1,LUPRI)
         END IF
         VAL(1) = -DDOT(KZYVA,WRK(KVECA),1,WRK(KFREE),1)
         VAL(2) = -DDOT(KZYVA,WRK(KVECA+KZYVA),1,WRK(KFREE),1)
         VAL(3) = D0
         VAL(4) = D0
         CALL R3DRV(ISYMA,ISYMB,ISYMC,KZYVA,KZYVB,KZYVC,
     &        WRK(KVECA),WRK(KVECB),WRK(KVECC),DAMPING,
     &        XINDX,UDV,MJWOP,WRK(KFREE),LFREE,VAL(3))
         CALL R3DRV(ISYMA,ISYMC,ISYMB,KZYVA,KZYVC,KZYVB,
     &        WRK(KVECA),WRK(KVECC),WRK(KVECB),DAMPING,
     &        XINDX,UDV,MJWOP,WRK(KFREE),LFREE,VAL(3))
         TMPVAL(1) = VAL(1) + VAL(3)
         TMPVAL(2) = VAL(2) + VAL(4)
         HYPVAL(1) = HYPVAL(1) + VAL(1) + VAL(3)      
         HYPVAL(2) = HYPVAL(2) + VAL(2) + VAL(4)      
C
         IF (DAMPING.GT.D0) THEN
            CALL T3DRV(1,ISYMA,ISYMB,ISYMC,WRK(KVECB+KZYVB),WRK(KVECC),
     &           .FALSE.,WRK(KVECA),-FREQB,-FREQC,XINDX,UDV,PV,MJWOP,
     &           WRK(KFREE),LFREE,CMO,FC,FV)
            VAL(1) = -DDOT(KZYVA,WRK(KVECA),1,WRK(KFREE),1)
            VAL(2) = -DDOT(KZYVA,WRK(KVECA+KZYVA),1,WRK(KFREE),1)
            VAL(3) = D0
            VAL(4) = D0
            CALL R3DRV(ISYMA,ISYMB,ISYMC,KZYVA,KZYVB,KZYVC,
     &           WRK(KVECA),WRK(KVECB+KZYVB),WRK(KVECC),DAMPING,
     &           XINDX,UDV,MJWOP,WRK(KFREE),LFREE,VAL(3))
            CALL R3DRV(ISYMA,ISYMC,ISYMB,KZYVA,KZYVC,KZYVB,
     &           WRK(KVECA),WRK(KVECC),WRK(KVECB+KZYVB),DAMPING,
     &           XINDX,UDV,MJWOP,WRK(KFREE),LFREE,VAL(3))
            TMPVAL(1) = VAL(1) + VAL(3)
            TMPVAL(2) = VAL(2) + VAL(4)
            TMPVAL(1) = TMPVAL(1) - VAL(2) - VAL(4)
            TMPVAL(2) = TMPVAL(2) + VAL(1) + VAL(3)
            HYPVAL(1) = HYPVAL(1) - VAL(2) - VAL(4)
            HYPVAL(2) = HYPVAL(2) + VAL(1) + VAL(3)
C
            CALL T3DRV(1,ISYMA,ISYMB,ISYMC,WRK(KVECB),WRK(KVECC+KZYVC),
     &           .FALSE.,WRK(KVECA),-FREQB,-FREQC,XINDX,UDV,PV,MJWOP,
     &           WRK(KFREE),LFREE,CMO,FC,FV)
            VAL(1) = -DDOT(KZYVA,WRK(KVECA),1,WRK(KFREE),1)
            VAL(2) = -DDOT(KZYVA,WRK(KVECA+KZYVA),1,WRK(KFREE),1)
            VAL(3) = D0
            VAL(4) = D0
            CALL R3DRV(ISYMA,ISYMB,ISYMC,KZYVA,KZYVB,KZYVC,
     &           WRK(KVECA),WRK(KVECB),WRK(KVECC+KZYVC),DAMPING,
     &           XINDX,UDV,MJWOP,WRK(KFREE),LFREE,VAL(3))
            CALL R3DRV(ISYMA,ISYMC,ISYMB,KZYVA,KZYVC,KZYVB,
     &           WRK(KVECA),WRK(KVECC+KZYVC),WRK(KVECB),DAMPING,
     &           XINDX,UDV,MJWOP,WRK(KFREE),LFREE,VAL(3))
            TMPVAL(1) = TMPVAL(1) - VAL(2) - VAL(4)
            TMPVAL(2) = TMPVAL(2) + VAL(1) + VAL(3)
            HYPVAL(1) = HYPVAL(1) - VAL(2) - VAL(4)      
            HYPVAL(2) = HYPVAL(2) + VAL(1) + VAL(3)      
C
            CALL T3DRV(1,ISYMA,ISYMB,ISYMC,WRK(KVECB+KZYVB),
     &           WRK(KVECC+KZYVC),
     &           .FALSE.,WRK(KVECA),-FREQB,-FREQC,XINDX,UDV,PV,MJWOP,
     &           WRK(KFREE),LFREE,CMO,FC,FV)
            VAL(1) = -DDOT(KZYVA,WRK(KVECA),1,WRK(KFREE),1)
            VAL(2) = -DDOT(KZYVA,WRK(KVECA+KZYVA),1,WRK(KFREE),1)
            VAL(3) = D0
            VAL(4) = D0
            CALL R3DRV(ISYMA,ISYMB,ISYMC,KZYVA,KZYVB,KZYVC,
     &           WRK(KVECA),WRK(KVECB+KZYVB),WRK(KVECC+KZYVC),DAMPING,
     &           XINDX,UDV,MJWOP,WRK(KFREE),LFREE,VAL(3))
            CALL R3DRV(ISYMA,ISYMC,ISYMB,KZYVA,KZYVC,KZYVB,
     &           WRK(KVECA),WRK(KVECC+KZYVC),WRK(KVECB+KZYVB),DAMPING,
     &           XINDX,UDV,MJWOP,WRK(KFREE),LFREE,VAL(3))
            TMPVAL(1) = TMPVAL(1) - VAL(1) - VAL(3)
            TMPVAL(2) = TMPVAL(2) - VAL(2) - VAL(4)
            HYPVAL(1) = HYPVAL(1) - VAL(1) - VAL(3)      
            HYPVAL(2) = HYPVAL(2) - VAL(2) - VAL(4)      
         END IF
C
         IF (IPRABS.GE.0) WRITE(LUPRI,'(A15,4F15.6)') 
     &        ' Na T[3] NbNc  ',
     &        TMPVAL(1),TMPVAL(2),HYPVAL(1),HYPVAL(2)
C
C     Save value for beta, which equals minus the value of the QRF.
C
         RES_BETA(IQRF,1) = -HYPVAL(1)
         RES_BETA(IQRF,2) = -HYPVAL(2)
C
C     End loop over QRFs
C
      END DO
C
      RETURN
      END
      SUBROUTINE QRCHK(DOHYP,ALAB,BLAB,CLAB,ISYMA,ISYMB,ISYMC,
     &                 FREQA,FREQB,FREQC)
#include <implicit.h>
#include <priunit.h>
#include <absorp.h>
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
      IF (ABS_MCD) THEN
         IF (ALAB(2:8).NE.'DIPLEN' .OR. BLAB(2:8).NE.'DIPLEN' .OR.
     &       CLAB(2:8).NE.'ANGMOM') DOHYP = .FALSE.
      END IF
C
C     Check if equivalent QRF is in the list already.
C
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
            DO IFR=1,NFREQ_ALPHA
               IF (ABS(FREQ_ALPHA(IFR)-ABS(FREQ(I))).LT.THRZERO) THEN
                  NEWFRQ = .FALSE.
               END IF
            END DO
            IF (NEWFRQ) THEN
               IF (NFREQ_ALPHA.GE.MXFREQ) THEN
                  WRITE(LUPRI,'(2(/A),I4,A,/A)')
     & ' The specified calculation requires more than the allowed',
     & ' number of frequencies in the LR solver MXFREQ=',MXFREQ,'.',
     & ' The program will stop.'
                  CALL QUIT('Too many frequencies in LR solver.')
               END IF
               NFREQ_ALPHA = NFREQ_ALPHA + 1
               FREQ_ALPHA(NFREQ_ALPHA) = ABS(FREQ(I))
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
      SUBROUTINE QRRDVE(ISYMA,ISYMB,ISYMC,ALAB,BLAB,CLAB,
     &     FREQA,FREQB,FREQC,KZYVA,KZYVB,KZYVC,VECA,VECB,VECC)
C
C PURPOSE: Read in response vectors for quadratic response. 
C
#include <implicit.h>
#include <priunit.h>
#include <absorp.h>
#include <inftap.h>
#include <rspprp.h>
#include <inflr.h>
#include <qrinf.h>
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
      CALL REARSP(LURSP,2*KZYVA,VECA,ALAB,'COMPLEX ',ABS(FREQA),D0,
     &     ISYMA,0,THCLR,FOUND,CONV,ANTSYM)
      IF (.NOT. (FOUND .AND. CONV)) THEN
         IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/3A,F7.5,A,I3,/A)') ' Response label ',ALAB,
     &           ' with frequency ',FREQA, ' and symmetry',
     &           ISYMA,' not found on file LURSP'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/3A,F7.5,/A,I3,A)') ' @WARNING>>>>'//
     &           ' Response label ',ALAB,
     &           ' with frequency ',FREQA, ' and symmetry',
     &           ISYMA,' not converged on file LURSP'
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
      CALL REARSP(LURSP,2*KZYVB,VECB,BLAB,'COMPLEX ',ABS(FREQB),D0,
     &     ISYMB,0,THCLR,FOUND,CONV,ANTSYM)
      IF (.NOT. (FOUND .AND. CONV)) THEN
         IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/3A,F7.5,A,I3,/A)') ' Response label ',BLAB,
     &           ' with frequency ',FREQB, ' and symmetry',
     &           ISYMB,' not found on file LURSP'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/3A,F7.5,/A,I3,A)') ' @WARNING>>>>'//
     &           ' Response label ',BLAB,
     &           ' with frequency ',FREQB, ' and symmetry',
     &           ISYMB,' not converged on file LURSP'
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
      CALL REARSP(LURSP,2*KZYVC,VECC,CLAB,'COMPLEX ',ABS(FREQC),D0,
     &     ISYMC,0,THCLR,FOUND,CONV,ANTSYM)
      IF (.NOT. (FOUND .AND. CONV)) THEN
         IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/3A,F7.5,A,I3,/A)') ' Response label ',CLAB,
     &           ' with frequency ',FREQC, ' and symmetry',
     &           ISYMC,' not found on file LURSP'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/3A,F7.5,/A,I3,A)') ' @WARNING>>>>'//
     &           ' Response label ',CLAB,
     &           ' with frequency ',FREQC, ' and symmetry',
     &           ISYMC,' not converged on file LURSP'
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
      SUBROUTINE R3DRV(ISYMA,ISYMB,ISYMC,KZYVA,KZYVB,KZYVC,
     &                 VECA,VECB,VECC,FACTOR,
     &                 XINDX,UDV,MJWOP,WRK,LWRK,RESULT)
C
C Compute the third-order relaxation contribution times -i*FACTOR (where
C FACTOR is equal to the common damping parameter). Results for real and
C imaginary parts are added to RESULT(1) and RESULT(2), respectively.
C
C     Na*R[3]*Nb*Nc
C
#include <implicit.h>
#include <priunit.h>
#include <infvar.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <inforb.h>
#include <infdim.h>
#include <qrinf.h>
C
      LOGICAL TDM, NORHO2
c
      PARAMETER ( D0=0.0D0, DH=0.5D0, D1=1.0D0, NTERMS=10 )
      DIMENSION WRK(LWRK),XINDX(*),MJWOP(2,MAXWOP,8),UDV(NASHDI,NASHDI)
      DIMENSION VECA(KZYVA),VECB(KZYVB),VECC(KZYVC),RESULT(2)
      DIMENSION TMP(2,NTERMS)
C
      KFREE=1
      LFREE=LWRK
C
C     There is no third-order relaxation contribution when the 
C     damping parameter is equal to zero. There is also no contribution
C     in Hartree-Fock since the average value of triple exciteations must
C     vanish
C
      IF (FACTOR.EQ.D0 .OR. TDHF) RETURN
C
      CALL DZERO(TMP,NTERMS)
      CALL MEMGET('REAL',KZYMAT,NORBT*NORBT,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDEN1,NASHT*NASHT,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KCREF,MZCONF(1),WRK,KFREE,LFREE)
C
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
      IPRTMP=IPRRSP
      IPRRSP=201
      IPRONE = 0
      CALL DSWAP(KZYVA/2,VECA(1),1,VECA(1+KZYVA/2),1)
      CALL DSWAP(KZYVA/2,VECA(1+KZYVA),1,VECA(1+KZYVA+KZYVA/2),1)
C
C     1/2*<0|[nC,[nB,nA]]|0>
C
      IF (MZWOPT(ISYMA).GT.0 .AND. MZWOPT(ISYMB).GT.0 .AND.
     &    MZWOPT(ISYMC).GT.0) THEN
         CALL TRZYM2(VECA,VECB,VECC,KZYVA,KZYVB,KZYVC,
     &        ISYMA,ISYMB,ISYMC,WRK(KZYMAT),MJWOP,WRK(KFREE),LFREE)
         CALL MELONE(WRK(KZYMAT),1,UDV,D1,TMP(1,1),IPRONE,'R3DRV')
         TMP(1,1) = TMP(1,1)*DH
         CALL TRZYM2(VECA(1+KZYVA),VECB,VECC,KZYVA,KZYVB,KZYVC,
     &        ISYMA,ISYMB,ISYMC,WRK(KZYMAT),MJWOP,WRK(KFREE),LFREE)
         CALL MELONE(WRK(KZYMAT),1,UDV,D1,TMP(2,1),IPRONE,'R3DRV')
         TMP(2,1) = TMP(2,1)*DH
      END IF
C
C     <0|[nA,nC] - [nA,nC](+)|0B>, where (+) is the dagger
C
      IF (MZWOPT(ISYMA).GT.0 .AND. MZCONF(ISYMB).GT.0 .AND.
     &    MZWOPT(ISYMC).GT.0) THEN
         ILSYM  = IREFSY
         IRSYM  = MULD2H(IREFSY,ISYMB)
         NCL    = MZCONF(1)
         NCR    = MZCONF(ISYMB)
         TDM=.TRUE.
         NORHO2=.TRUE.
         CALL DZERO(WRK(KDEN1),NASHT*NASHT)
         CALL RSPDM(ILSYM,IRSYM,NCL,NCR,WRK(KCREF),VECB,
     &              WRK(KDEN1),DUMMY,0,0,TDM,NORHO2,XINDX,WRK(KFREE),
     &              1,LFREE)
         OVLAP = D0
         IF (ILSYM.EQ.IRSYM)
     &         OVLAP = DDOT(NCL,WRK(KCREF),1,VECB,1)
C   
         CALL MEMGET('REAL',KTMP,NORBT*NORBT,WRK,KFREE,LFREE)
         CALL TRZYMT(1,VECA,VECC,KZYVA,KZYVC,ISYMA,ISYMC,WRK(KTMP),
     &               MJWOP,WRK,LWRK)
         CALL DCOPY(NORBT*NORBT,WRK(KTMP),1,WRK(KZYMAT),1)
         CALL DGEMM('T','N',NORBT,NORBT,NORBT,-1.0D0,
     &              WRK(KTMP),NORBT,1.0D0,1,1.0D0,WRK(KZYMAT),NORBT)
         CALL MELONE(ZYMAT,ISYMV1,DEN1,OVLAP,TMP(1,2),IPRONE,'R3DRV')
C
         CALL TRZYMT(1,VECA(1+KZYVA),VECC,KZYVA,KZYVC,ISYMA,ISYMC,
     &               WRK(KTMP),MJWOP,WRK,LWRK)
         CALL DCOPY(NORBT*NORBT,WRK(KTMP),1,WRK(KZYMAT),1)
         CALL DGEMM('T','N',NORBT,NORBT,NORBT,-1.0D0,
     &              WRK(KTMP),NORBT,1.0D0,1,1.0D0,WRK(KZYMAT),NORBT)
         CALL MELONE(ZYMAT,ISYMV1,WRK(KDEN1),OVLAP,TMP(2,2),
     &               IPRONE,'R3DRV')
      END IF
C
C     1/2*<0|nCnB + nC(+)nB(+)|0A>, where (+) is the dagger
C
      IF (MZCONF(ISYMA).GT.0 .AND. MZWOPT(ISYMB).GT.0 .AND.
     &    MZWOPT(ISYMC).GT.0) THEN
         CALL MEMGET('REAL',KDEN2,N2ASHX*N2ASHX,WRK,KFREE,LFREE)
         ILSYM  = IREFSY
         IRSYM  = MULD2H(IREFSY,ISYMA)
         NCL    = MZCONF(1)
         NCR    = MZCONF(ISYMA)
         TDM=.TRUE.
         NORHO2=.FALSE.
         CALL DZERO(WRK(KDEN1),NASHT*NASHT)
         CALL DZERO(WRK(KDEN2),N2ASHX*N2ASHX)
         CALL RSPDM(ILSYM,IRSYM,NCL,NCR,WRK(KCREF),VECA,
     &              WRK(KDEN1),WRK(KDEN2),0,0,TDM,NORHO2,XINDX,
     &              WRK(KFREE),1,LFREE)
         OVLAP = D0
         IF (ILSYM.EQ.IRSYM)
     &         OVLAP = DDOT(NCL,WRK(KCREF),1,VECA,1)
         CALL TWOMOM(ISYMB,ISYMC,KZYVB,KZYVC,VECB,VECC,MJWOP,
     &               WRK(KDEN1),WRK(KDEN2),WRK(KFREE),LFREE,
     &               TMP(1,3))
      END IF
C
      IPRRSP=IPRTMP
C
C     Multiply by -i*W (minus i times the gamma factor).
C
      RESULT(1) =  FACTOR*TMP(1,2)
      RESULT(2) = -FACTOR*TMP(1,1)
C
      RETURN
      END
      SUBROUTINE TWOMOM(ISYMB,ISYMC,KZYVB,KZYVC,VECB,VECC,MJWOP,
     &                 DEN1,DEN2,WRK,LWRK,RESULT)
C
C Calculate <0|nCnB + nC(+)nB(+)|0A> given the one- and two-electron
C transition densities.
C
#include <implicit.h>
#include <priunit.h>
#include <infvar.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <inforb.h>
#include <infdim.h>
#include <qrinf.h>
C
      DIMENSION WRK(LWRK),MJWOP(2,MAXWOP,8)
      DIMENSION DEN1(NASHDI,NASHDI),DEN2(N2ASHX,N2ASHX)
      DIMENSION VECB(KZYVB),VECC(KZYVC),RESULT(2)
C
      KFREE=1
      LFREE=LWRK
      CALL QUIT('TWOMOM: to be implemented....')
C
      RETURN
      END
