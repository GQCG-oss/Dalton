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
C ============================================================================
C     absorp.F: Main author Patrick Norman, written 2001.
C     This implementation is published in:
C        P. Norman, D.M. Bishop, H.J.Aa. Jensen, and J. Oddershede,
C        "Near-resonant absorption in the time-dependent self-consistent
C        field and multiconfigurational self-consistent field approximations",
C        J. Chem. Phys. 115 (2001) 10323-10334.
C ============================================================================
C
      SUBROUTINE ABSORP_INPUT(WORD)
C
C     Purpose:
C     Read in user settings for imaginary polarizabilities describing
C     absorption in the optical processes.
C
#include <implicit.h>
C
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
C
      LOGICAL NEWDEF
      PARAMETER ( CM_TO_AU = 4.556333D-6 )
      PARAMETER ( NTABLE = 10 , D0 = 0.0D0 , THRFRQ = 1.0D-14 )
      CHARACTER PROMPT*1, WORD*7, WORD1*7, TABLE(NTABLE)*7
      CHARACTER*8 DIPLEN(3)
      DATA DIPLEN/'XDIPLEN','YDIPLEN','ZDIPLEN'/
C
      DATA TABLE /'.ALPHA ','.FREQUE','.THCLR ','.MAXIT ','MAXITO',
     &            '.PRINT ','.LIFETI','.EXCITA','.MAX IT','THCPP '/
C
      NEWDEF = (WORD .EQ. '*ABSORP')
      ICHANG = 0
C
C     Set default values
C
      IF (NEWDEF) THEN
         ALPHA_FLG  = .FALSE.
         BETA_FLG   = .FALSE.
         GAMMA_FLG  = .FALSE.
C     Common default value of the excited state lifetimes is set to 
C     be 1000 cm-1.
         LIFETIME = 1000*CM_TO_AU
         NFREQ_ALPHA = 1
         FREQ_ALPHA(1) = D0
         THCLR_ABSORP = 1.0D-3
         THCPP_ABSORP = 1.0D-3
         IPRABS = 0
         NEXCITED_STATES = 5
         MAX_MACRO  = 5
         MAX_MICRO  = 20
         MAX_ITORB = 5
         WORD1 = WORD
  100    CONTINUE
C
C     Read in input
C
            READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') GOTO 100
            IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) GOTO (1,2,3,4,5,6,7,8,9,10), I
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
               ALPHA_FLG=.TRUE.
               GOTO 100
 2             CONTINUE
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
C     User input of lifetime in reciprocal centimeters. Convert to au.
                  READ (LUCMD,*) LIFETIME
                  LIFETIME = LIFETIME*CM_TO_AU
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
C
C     Process user input
C
      ABSORP = ALPHA_FLG .OR. BETA_FLG .OR. GAMMA_FLG
      IF (.NOT.ABSORP) THEN
         WRITE(LUPRI,'(/A)') ' Absorption input ignored because:'
         WRITE(LUPRI,'(A,L2)')' (Hyper)polarizability not requested:',
     &        ' ALPHA_FLG = BETA_FLG = GAMMA_FLG = ', ABSORP
      ELSE 
C
C     Put dipole operators in lists.
C
         CALL IZERO(NOPER,8)
         DO ILABEL=1,3
            ISYM = ISYMAX(ILABEL,1)+1
            NOPER(ISYM) = NOPER(ISYM) + 1
            LABOP(NOPER(ISYM),ISYM) = DIPLEN(ILABEL)
         END DO
C
         CALL HEADER('Variable settings for absorption calculation',0)
C
         IF (ALPHA_FLG) THEN
            WRITE (LUPRI,'(A,L4,/A,5(4F12.8,/,16X))')
     &     ' Linear polarizability calculation requested:  ALPHA_FLG =',
     &           ALPHA_FLG,
     &           ' at frequencies:',(FREQ_ALPHA(I),I=1,NFREQ_ALPHA)
         END IF
C
         IF (NEXCITED_STATES .GT. MXSTATES) THEN
            WRITE(LUPRI,'(/2A,/)') ' >>> Warning the number of excited',
     &           ' exceeds maximun, the maximum value is used <<<'
            NEXCITED_STATES = MXSTATES
         END IF
         WRITE(LUPRI,'(A,F8.1)')
     &      ' Common lifetime for excited states (cm-1)   : LIFETIME  ='
     &      ,LIFETIME/CM_TO_AU
         WRITE(LUPRI,'(A,I4)')
     &      ' Number of excited states used in start iter.: NEXCIT    ='
     &       ,NEXCITED_STATES
         WRITE(LUPRI,'(A,1P,D8.1)')
     &      ' Threshold for convergence of lin. resp. eq.s: THCLR     ='
     &      ,THCLR_ABSORP
         WRITE(LUPRI,'(A,I4)')
     &      ' Maximum iterations in complex LR solver     : MAX_MACRO ='
     &       ,MAX_MACRO
         WRITE(LUPRI,'(A,I4)')
     &      ' Maximum iterations in real LR solver        : MAX_MICRO ='
     &       ,MAX_MICRO
         WRITE(LUPRI,'(A,I4)')
     &      ' Maximum iter. in optimal orbital algorithm  : MAX_ITORB ='
     &       ,MAX_ITORB
         WRITE(LUPRI,'(A,I4)')
     &      ' Print level in absorption modules           : IPRABS    ='
     &       ,IPRABS
C
      END IF
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
#include <iratdef.h>
#include <absorp.h>
#include <dummy.h>
#include <infvar.h>
#include <qrinf.h>
#include <infrsp.h>
#include <inftap.h>
C
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
C
      IPRRSP = IPRABS
      KFREE = 1
      LFREE = LWRK
C
C     Determine linear response vectors
C
      CALL ABSVEC1(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     &             XINDX,WRK(KFREE),LFREE)
C
C     Result section
C
      CALL ABSRESULT(WRK(KFREE),LFREE)
C
      RETURN
      END
      SUBROUTINE ABSVEC1(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     &                   XINDX,WRK,LWRK)
C
C PURPOSE:
C     Solve the linear response equations and store response vectors on
C     disk.
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
C
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
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
     &           FOCK,FC,FV,FCAC,H2AC,XINDX,WRK(KFREE),LFREE)
C
            DO IOPER=1,NOPER(KSYMOP)
               CALL ABSCTL(IOPER,LABOP(IOPER,KSYMOP),
     &              WRK(KREDE),WRK(KREDS),
     &              WRK(KREDZ),WRK(KREDGD),WRK(KREDZGD),
     &              WRK(KEIVEC),WRK(KIBTYP),
     &              CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,XINDX,
     &              WRK(KFREE),LFREE)
            END DO
         END IF
      END DO
      IF (IPRABS.GT.0) CALL TIMER('ABSVEC',TIMSTR,TIMEND)
      RETURN
      END
      SUBROUTINE RESONANT(REDE,REDS,IBTYP,EIGVAL,RESIDUAL,EIGVEC,
     &     CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
C
C PURPOSE:
C     Solve the eigenvalue equation ( E[2] - w*S[2] )*X = 0 for 
C     the few lowest excited states to be used as startvectors.
C
#include <implicit.h>
#include <priunit.h>
#include <absorp.h>
#include <dummy.h>
#include <infvar.h>
#include <infrsp.h>
#include <qrinf.h>
#include <iratdef.h>
#include <infdim.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <symmet.h>
#include <wrkrsp.h>
#include <cbilrs.h>
#include "rspprp.h"
#include <inflr.h>
#include <infpp.h>
C
      DIMENSION REDE(*),REDS(*),IBTYP(*),EIGVAL(*),RESIDUAL(*),EIGVEC(*)
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
      CHARACTER*8 BLANK
C
      PARAMETER ( BLANK='        ', D0=0.0D0 )
C
      KFREE = 1
      LFREE = LWRK
C
C     Initialize variables 
C     SETZY and therby also MJWOP removed 2002-04-19
C
      THCRSP = THCPP_ABSORP
      MAXIT  = MAX_MICRO
      CALL RSPVAR(UDV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK(KFREE),LFREE)
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
      SUBROUTINE ABSREDUC(IOP,LABEL,REDE,REDS,REDZ,REDGD,REDZGD,IBTYP,
     &     REDVEC,FC,CMO,UDV,PV,XINDX,WRK,LWRK)
C
#include <implicit.h>
#include <priunit.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <absorp.h>
#include <inftap.h>
#include <ibndxdef.h>
#include <dummy.h>
      PARAMETER (D0 = 0.0D0)
C
C PURPOSE:
C     Create start vectors from the reduced eigenvalue problem
C     ( E[2] - {w+iW}*S[2] )* (NR + iNI) = B[1]
C
      CHARACTER*8 LABEL
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),XINDX(*),IBTYP(*),WRK(*)
      DIMENSION REDE(MAXRM*MAXRM),REDS(MAXRM*MAXRM),REDZ(2*MAXRM*MAXRM),
     &          REDGD(MAXRM),REDZGD(2*MAXRM),REDVEC(2*MAXRM,MXFREQ)
C
      KFREE = 1
      LFREE = LWRK
      CALL MEMGET('INTE',KIPIV,MAXRM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KBVEC,KZYVAR,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KGD,KZYVAR,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KREDE,KZYRED*KZYRED,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KREDS,KZYRED*KZYRED,WRK,KFREE,LFREE)
C
C Construct the reduced gradient with trial vectors from  LURSP3
C
      CALL GETGPV(LABEL,FC,DUMMY,CMO,UDV,PV,XINDX,ANTSYM,
     &   WRK(KFREE),LFREE)
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
      GRDNRM = DNORM2(KZYRED,REDGD,1)
      IF (GRDNRM.EQ.D0) THEN
         WRITE(LUPRI,'(4(/,5X,A))') 
     &        ' ABSREDUC: The reduced property gradient has',
     &        ' zero norm in the space of trial vectors.',
     &        ' The reason may be that you should employ a higher',
     &        ' point group in your calculation.'
         CALL QUIT('Zero norm of reduced property gradient.')
      END IF
C
      IF (IPRABS.GE.20) THEN
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
      DO IFREQ=1,NFREQ_ALPHA
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
         CALL DAXPY(KZYRED*KZYRED,-FREQ_ALPHA(IFREQ),WRK(KREDS),1,
     &              REDZ,2)
         CALL DAXPY(KZYRED*KZYRED,-LIFETIME,WRK(KREDS),1,REDZ(2),2)
C
         IF (IPRABS.GT.10) THEN
            WRITE(LUPRI,'(/,5X,A,/,5X,2(A,D10.4),/,5X,A)') 
     &           ' Reduced ( E[2] - {w+iW}*S[2] )  matrix',
     &           ' with w =', FREQ_ALPHA(IFREQ),' and W =', LIFETIME,
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
     &           ' with w =', FREQ_ALPHA(IFREQ),' and W =', LIFETIME,
     &           '================================================' 
            CALL OUTPUT(REDVEC(1,IFREQ),1,2,1,KZYRED,2,KZYRED,1,LUPRI)
         END IF
C
         RES_ALPHA(1,IFREQ,IOP,KSYMOP) = 
     &        DDOT(KZYRED,REDGD,1,REDVEC(1,IFREQ),2)
         RES_ALPHA(2,IFREQ,IOP,KSYMOP) = 
     &        DDOT(KZYRED,REDGD,1,REDVEC(2,IFREQ),2)
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
     &        ' Frequency   Lifetime        Real   Imaginary',
     &        '----------------------------------------------'
         DO IFREQ=1,NFREQ_ALPHA
            WRITE(LUPRI,'(6X,F10.4,F11.6,2F12.6)') FREQ_ALPHA(IFREQ),
     &           LIFETIME, (RES_ALPHA(I,IFREQ,IOP,KSYMOP), I=1,2)
         END DO
      END IF
C
      RETURN
      END
      SUBROUTINE ABSRESULT(WRK,LWRK)
C
#include <implicit.h>
#include <priunit.h>
#include <absorp.h>
#include <dummy.h>
#include <infvar.h>
#include <qrinf.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <symmet.h>
#include <cbilrs.h>
#include <infrsp.h>
#include <inforb.h>
C
      PARAMETER ( AU_TO_EV=27.2114, CM_TO_AU = 4.556333D-6 )
      DIMENSION WRK(*)
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
     &              EXC_ENERGY(ISTATE,ISYM)*AU_TO_EV,' eV'
            END DO
         END IF
      END DO
C
      IX = 1
      IY = 1
      IZ = 1
      IF (ISYMAX(1,1).EQ.ISYMAX(2,1)) IY = IY + 1
      IF (ISYMAX(1,1).EQ.ISYMAX(3,1)) IZ = IZ + 1
      IF (ISYMAX(2,1).EQ.ISYMAX(3,1)) IZ = IZ + 1
C
      CALL AROUND('Polarizability with finite lifetime')
      WRITE(LUPRI,'(2A,/A7,A2,4A17,/2A)') 
     &     '=======================================',
     &     '=======================================',
     &     ' Freq.',' ','xx     ','yy     ','zz     ','ave    ',
     &     '=======================================',
     &     '======================================='
      DO IFREQ=1,NFREQ_ALPHA
         AR_XX = RES_ALPHA(1,IFREQ,IX,ISYMAX(1,1)+1)
         AI_XX = RES_ALPHA(2,IFREQ,IX,ISYMAX(1,1)+1)
         AR_YY = RES_ALPHA(1,IFREQ,IY,ISYMAX(2,1)+1)
         AI_YY = RES_ALPHA(2,IFREQ,IY,ISYMAX(2,1)+1)
         AR_ZZ = RES_ALPHA(1,IFREQ,IZ,ISYMAX(3,1)+1)
         AI_ZZ = RES_ALPHA(2,IFREQ,IZ,ISYMAX(3,1)+1)
         AR = ( AR_XX + AR_YY + AR_ZZ )/3
         AI = ( AI_XX + AI_YY + AI_ZZ )/3
         WRITE(LUPRI,'(1X,F6.4,A2,4(F7.2,A2,F7.2,A1))') 
     &        FREQ_ALPHA(IFREQ),' ',
     &        AR_XX,' +',AI_XX,'i',AR_YY,' +',AI_YY,'i',
     &        AR_ZZ,' +',AI_ZZ,'i',AR,' +',AI,'i'
      END DO
      WRITE (LUPRI,'(2A,/,1X,A,F10.6,A,F7.4,A,F7.1,A,/,/)')
     &     '=======================================',
     &     '=======================================',
     &     'Common lifetime equals',
     &     LIFETIME,' au =',
     &     LIFETIME*AU_TO_EV,' eV =',
     &     LIFETIME/CM_TO_AU,' cm-1'
C
      RETURN
      END
      SUBROUTINE ABSRESID(IOP,LABEL,CONVERGED,REDVEC,IBTYP,
     &              FC,CMO,UDV,PV,XINDX,WRK,LWRK)
C
C PURPOSE:
C     Compute the imaginary residual vector to the linear response 
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
#include <infrsp.h>
#include <wrkrsp.h>
#include <absorp.h>
#include <inftap.h>
#include <ibndxdef.h>
#include <dummy.h>
C
      LOGICAL CONVERGED
      CHARACTER*8 LABEL
      DIMENSION REDVEC(2*MAXRM,MXFREQ),IBTYP(MAXRM),
     &     CMO(*),UDV(*),PV(*),FC(*),XINDX(*),WRK(*)
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
      DO IFREQ=1,NFREQ_ALPHA
C
C     Initialize residual vector
C
         CALL GETGPV(LABEL,FC,DUMMY,CMO,UDV,PV,XINDX,ANTSYM,
     &      WRK(KFREE),LFREE)
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
            FR3 = (FREQ_ALPHA(IFREQ)*FR1-LIFETIME*FI1)
            FR4 = (FREQ_ALPHA(IFREQ)*FR2-LIFETIME*FI2)
            FI3 = (FREQ_ALPHA(IFREQ)*FI1+LIFETIME*FR1)
            FI4 = (FREQ_ALPHA(IFREQ)*FI2+LIFETIME*FR2)
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
         DO IFREQ = 1,NFREQ_ALPHA
            WRITE(LUPRI,'(1X,A,F7.4,3X,A,F9.6,3X,A,1P,3D9.2)')
     &           'Frequency:',FREQ_ALPHA(IFREQ),
     &           'Lifetime:',LIFETIME,
     &           'Residual:',(RESID(I,IFREQ,IOP,KSYMOP),I=1,3)
         END DO
      END IF
C
      RETURN
      END
      SUBROUTINE ABSLR(IOP,LABEL,REDE,REDS,REDZ,REDGD,REDZGD,IBTYP,
     &                 REDVEC,CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,XINDX,
     &                 WRK,LWRK)
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
C     Solve the coupled equations in this iteration
C
C      (i)    ( E[2] - w*S[2] ) NR = B[1] - W*S[2]*NI
C     (ii)    ( E[2] - w*S[2] ) NI =        W*S[2]*NR
C
      CHARACTER*8 LABEL,BLANK
      PARAMETER (BLANK='        ', D0=0.0D0)
      DIMENSION REDE(MAXRM*MAXRM),REDS(MAXRM*MAXRM),REDZ(2*MAXRM*MAXRM),
     &          REDGD(MAXRM),REDZGD(2*MAXRM),REDVEC(2*MAXRM,MXFREQ),
     &          IBTYP(MAXRM)
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
C
      KFREE = 1
      LFREE = LWRK
      CALL MEMGET('REAL',KGD1,KZYVAR,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KGD2,KZYVAR,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KSLI,KZYVAR,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KBVEC,KZYVAR,WRK,KFREE,LFREE)
C
      DO IFREQ=1,NFREQ_ALPHA
C
         IF (RESID(3,IFREQ,IOP,KSYMOP).LT.THCLR_ABSORP) THEN
            IF (IPRABS.GE.2) WRITE(LUPRI,'(/,1X,A,F8.4,A)') 
     &           '>>>>> Frequency:',FREQ_ALPHA(IFREQ),' converged <<<<<'
            GO TO 100
         END IF
C
C     Construct the right-hand-side of equations (i) and (ii)
C
         CALL GETGPV(LABEL,FC,FV,CMO,UDV,PV,XINDX,ANTSYM,
     &      WRK(KFREE),LFREE)
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
            FR1 = LIFETIME * REDVEC(4*I-3,IFREQ)
            FR2 = LIFETIME * REDVEC(4*I-1,IFREQ)
            FI1 = LIFETIME * REDVEC(4*I-2,IFREQ)
            FI2 = LIFETIME * REDVEC(4*I,IFREQ)
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
     &        'Frequency = ',FREQ_ALPHA(IFREQ),
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
         IF (RESID(1,IFREQ,IOP,KSYMOP).LT.THCLR_ABSORP) THEN
            IF (IPRABS.GE.2) WRITE(LUPRI,'(/,1X,A,1P,D9.2,A)') 
     &           '>>>>> Skip real equation because residual is only:',
     &           RESID(1,IFREQ,IOP,KSYMOP), ' <<<<<'
         ELSE
            IF (IPRABS.GE.2) WRITE(LUPRI,'(/,/,1X,A,F8.4,A)') 
     &           '>>>>> Solving real equation for frequency:',
     &           FREQ_ALPHA(IFREQ),' <<<<<<'
            THCRSP = THCLR_ABSORP*DNORM_NT/DNORM_NR
            CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
     &           .TRUE.,BLANK,BLANK,WRK(KGD1),REDGD,REDE,REDS,
     &           IBTYP,FREQ_ALPHA(IFREQ),RESIDUAL,REDVEC,XINDX,
     &           WRK(KFREE),LFREE)
         END IF
C
         ELSE
C
         IF (RESID(2,IFREQ,IOP,KSYMOP).LT.THCLR_ABSORP) THEN
            IF (IPRABS.GE.2) WRITE(LUPRI,'(/,1X,A,1P,D9.2,A)') 
     &        '>>>>> Skip imaginary equation because residual is only:',
     &           RESID(2,IFREQ,IOP,KSYMOP), ' <<<<<'
         ELSE
            IF (IPRABS.GE.2) WRITE(LUPRI,'(/,/,1X,A,F8.4,A)') 
     &           '>>>>> Solving imaginary equation for frequency:',
     &           FREQ_ALPHA(IFREQ),' <<<<<<'
            THCRSP = THCLR_ABSORP*DNORM_NT/DNORM_NI
            CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
     &           .TRUE.,BLANK,BLANK,WRK(KGD2),REDGD,REDE,REDS,
     &           IBTYP,FREQ_ALPHA(IFREQ),RESIDUAL,REDVEC,XINDX,
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
      SUBROUTINE ABSWRT(LABEL,IOP,IBTYP,REDVEC,WRK,LWRK)
C
C PURPOSE:
C
C     Construct and write response vectors to file
C
C     From file we read:
C         LURSP3 contains trial vectors
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
#include <inftap.h>
#include <ibndxdef.h>
C
      CHARACTER*8 LABEL,BLANK
      PARAMETER (BLANK='        ', D0 = 0.0D0, D1=1.0D0)
      DIMENSION REDVEC(2*MAXRM,MXFREQ),IBTYP(MAXRM),WRK(*)
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
      DO IFREQ=1,NFREQ_ALPHA
C
C     Read trial vectors and add contributions to response vector
C
         CALL REWSPL(LURSP3)
         KRSP3 = 0
         IF (KOFFTY.EQ.1) THEN 
#if defined (VAR_SPLITFILES)
            CALL READST(LURSP3,IREC,DUMMY,KRSP3,.FALSE.)
#else        
            READ (LURSP3)
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
               CALL READST(LURSP3,IREC,WRK(KBVEC),KRSP3,.TRUE.)
               IF (IREC.NE.KZCONF) CALL QUIT('Read error in ABSWRT')
#else        
               CALL READT(LURSP3,KZCONF,WRK(KBVEC))
#endif
               CALL DAXPY(KZCONF,FR1,WRK(KBVEC),1,WRK(KNR),1)
               CALL DAXPY(KZCONF,FR2,WRK(KBVEC),1,WRK(KNR+KZVAR),1)
               CALL DAXPY(KZCONF,FI1,WRK(KBVEC),1,WRK(KNI),1)
               CALL DAXPY(KZCONF,FI2,WRK(KBVEC),1,WRK(KNI+KZVAR),1)
            ELSE
#if defined (VAR_SPLITFILES)
               CALL READST(LURSP3,IREC,WRK(KBVEC),KRSP3,.TRUE.)
               IF (IREC.NE.KZYWOP) CALL QUIT('Read error in ABSWRT')
#else        
               CALL READT(LURSP3,KZYWOP,WRK(KBVEC))
#endif
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
         CALL WRTRSP(LURSP,2*KZYVAR,WRK(KNR),LABEL,BLANK,
     &        FREQ_ALPHA(IFREQ),D0,KSYMOP,0,
     &        RESID(3,IFREQ,IOP,KSYMOP),D1)
C
      END DO
C
      RETURN
      END
      SUBROUTINE ABSCTL(IOPER,LABEL,
     &     REDE,REDS,REDZ,REDGD,REDZGD,REDVEC,IBTYP,
     &     CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,XINDX,
     &     WRK,LWRK)
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
      DIMENSION REDE(MAXRM*MAXRM),REDS(MAXRM*MAXRM),REDZ(2*MAXRM*MAXRM),
     &     REDGD(MAXRM),REDZGD(2*MAXRM),REDVEC(2*MAXRM,MXFREQ),
     &     IBTYP(MAXRM)
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
C
      KFREE = 1
      LFREE = LWRK 
C
      IF (IPRABS.GE.0) THEN
         WRITE(LUPRI,'(/2A,/3A,I2,A/,A,5(4F12.8,/,11X))') 
     &        ' ABSVEC1 -- Solving linear response equations',
     &        ' ( E[2] - {w-iW}*S[2] ) N = B[1]  ',
     &        ' ABSVEC1 -- for operator ', LABEL,
     &        ' of symmetry',KSYMOP,' at frequencies:',
     &        ' ABSVEC1 --',(FREQ_ALPHA(I),I=1,NFREQ_ALPHA)
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
      CALL ABSREDUC(IOPER,LABEL,
     &     REDE,REDS,REDZ,REDGD,REDZGD,IBTYP,
     &     REDVEC,FC,CMO,UDV,PV,XINDX,WRK(KFREE),LFREE)
C     
C     Compute the residual in this iteration
C     
      CALL ABSRESID(IOPER,LABEL,CONVERGED,REDVEC,IBTYP,
     &     FC,CMO,UDV,PV,XINDX,WRK(KFREE),LFREE)
C     
C     Check status
C     
      IF (CONVERGED) GOTO 900
      IF (ITER.GE.MAX_MACRO) GOTO 910
C     
C     Not converged, expand reduced space by solving coupled LR equations
C     
      CALL ABSLR(IOPER,LABEL,
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
     &        '*** THE REQUESTED',NFREQ_ALPHA,
     &        ' SOLUTION VECTORS CONVERGED'
      END IF
      GOTO 990
 910  CONTINUE
      IF (IPRABS.GE.0) THEN
         WRITE(LUPRI,'(/,1X,A,I4,A)')
     &        '*** THE REQUESTED',NFREQ_ALPHA,
     &        ' SOLUTION VECTORS NOT CONVERGED'
         WRITE(LUPRI,'(/,1X,A)')
     &        '--- MAXIMUM NUMBER OF ITERATIONS REACHED ---'
      END IF
      GOTO 990
C     
C     Construct and write response vectors to file
C     
 990  CONTINUE
      CALL ABSWRT(LABEL,IOPER,IBTYP,REDVEC,
     &     WRK(KFREE),LFREE)
C     
      IF (IPRABS.GE.0) THEN
         WRITE(LUPRI,'(2(/,1X,A),1P,D9.2,/,1X,A,/,1X,A,I4,A)')
     &        'RSP solution vectors written to file.',
     &        'Convergence of RSP solution vectors, threshold =',
     &        THCLR_ABSORP,
     & '-------------------------------------------------------------',
     &        '(dimension of reduced space:',KZYRED,')'
         DO IFREQ = 1,NFREQ_ALPHA
            WRITE(LUPRI,'(1X,A,F7.4,4X,A,F9.6,4X,A,1P,D9.2)')
     &           'Frequency:',FREQ_ALPHA(IFREQ),
     &           'Lifetime:',LIFETIME,
     &           'Residual:',RESID(3,IFREQ,IOPER,KSYMOP)
         END DO
         WRITE(LUPRI,'(/2A,3(/A))') 
     &        '@ Value of linear polarizability for operator: ',
     &        LABEL,
     &        '@ --------------------------------------------',
     &        '@ Frequency   Lifetime        Real   Imaginary',
     &        '@ --------------------------------------------'
         DO IFREQ=1,NFREQ_ALPHA
            WRITE(LUPRI,'(A,F10.4,F11.6,2F12.6)') '@',
     &           FREQ_ALPHA(IFREQ), LIFETIME, 
     &           (RES_ALPHA(I,IFREQ,IOPER,KSYMOP), I=1,2)
         END DO
      END IF
C
      RETURN
      END
