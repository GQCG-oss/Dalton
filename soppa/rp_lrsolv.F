C
C  /* Deck rp_lrsolv */
      SUBROUTINE RP_LRSOLV(ISYMTR,NEXCI,MAXIT,EXVAL,LEXVAL,RESINM,
     &                     LRESINM,CONV,LCONV,FOCKD,LFOCKD,
     &                     REDE,REDS,LMXRED,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, June 1997
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C     Andrea Ligabue, December 2003: linear response functions 
C                                    implemented
C
C     PURPOSE: Solve the RPA linear response equations 
C              using an AO-driven algorithm. 
C
#include <implicit.h>
#include <priunit.h>
C
#include <soppinf.h>
#include <ccsdsym.h>
C
      LOGICAL   NONEWT
C
      CHARACTER*3 CONV(LCONV)
C
      DIMENSION EXVAL(LEXVAL),  RESINM(LRESINM)
      DIMENSION FOCKD(LFOCKD),  WORK(LWORK)
      DIMENSION REDE(LMXRED),   REDS(LMXRED)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_LRSOLV')
C
C==============================================================
C     For checking, calculate E[2] and S[2] matrices explicitly
C     by using unit vectors as trial vectors.
C==============================================================
C
      IF (SOPCHK) THEN
C
         CALL RP_CHECK(FOCKD,LFOCKD,ISYMTR,WORK,LWORK)
C
      END IF
C
C-----------------------------------------------------------------
C     Calculate diagonal parts of E[2] and S[2] and write to disk.
C-----------------------------------------------------------------
C
      DTIME      = SECOND()
      CALL RP_DIAG(FOCKD,LFOCKD,ISYMTR,WORK,LWORK)
      DTIME      = SECOND()   - DTIME
      SOTIME(31) = SOTIME(31) + DTIME
C
C---------------------------------
C     1. allocation of work space.
C---------------------------------
C
      LPOINT = NEXCI
C
      KPOINT = 1
      KEND1  = KPOINT + LPOINT
      LWORK1 = LWORK  - KEND1
C
      CALL SO_MEMMAX ('RP_LRSOLV.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('RP_LRSOLV.1',' ',KEND1,LWORK)
C
C----------------------------------
C     Set up initial trial vectors.
C----------------------------------
C
      DTIME      = SECOND()
      CALL RP_TRIAL1(NNEWTR,WORK(KPOINT),LPOINT,ISYMTR,
     &               NEXCI,WORK(KEND1),LWORK1)
      DTIME      = SECOND()   - DTIME
      SOTIME(32) = SOTIME(32) + DTIME
C
C--------------------------------------------------------------
C     Initialize iteration counter, number of old trialvectors,
C     and logical for 'no new trial vectors.'
C--------------------------------------------------------------
C
      NIT    = 0
C
      NOLDTR = 0
C
      LREDOL = 0
C
      NONEWT = .FALSE.
C
C------------------
C     Write banner.
C------------------
C 
      WRITE(LUPRI,'(/,2X,A,/,16X,A,I1,/,2X,A,/)') 
     &     '*****************************************************'//
     &     '*************',
     &     'RPA iterations, Excitation symmetry ',ISYMTR,
     &     '*****************************************************'//
     &     '*************'
C
C----------------------------------------------------------------------
C     Iteration loop for solving the linear equation/eigenvalue problem
C----------------------------------------------------------------------
C
  100 CONTINUE
C
C-----------------------------------
C        Count number of iterations.
C-----------------------------------
C
         NIT = NIT + 1
C
         IF ( IPRSOP .GE. 2 ) THEN
C
            WRITE(LUPRI,'(/,2X,A,/,16X,I3,A,I1,/,2X,A,/)') 
     &           '==================================================='//
     &           '===============',NIT,
     &           '. RPA iteration, Excitation symmetry ',ISYMTR,
     &           '==================='//
     &           '==============================================='
C
         END IF
C
C--------------------------------------------------------------
C        Make E[2] linear transformation of trialvectors giving 
C        resultvectors.
C--------------------------------------------------------------
C
         DTIME      = SECOND()
         CALL RP_ERES(NOLDTR,NNEWTR,FOCKD,LFOCKD,ISYMTR,
     &                WORK(KEND1),LWORK1)
         DTIME      = SECOND()   - DTIME
         SOTIME(35) = SOTIME(35) + DTIME
C
C
C--------------------------------------------------------------
C        Make S[2] linear transformation of trialvectors giving
C        resultvectors.
C--------------------------------------------------------------
C
         DTIME      = SECOND()
         CALL RP_SRES(NOLDTR,NNEWTR,ISYMTR,WORK(KEND1),LWORK1)
         DTIME      = SECOND()   - DTIME
         SOTIME(40) = SOTIME(40) + DTIME
C
C--------------------------------------------------------
C        Set up and solve the reduced eigenvalue problem.
C--------------------------------------------------------
C
         LREDE  = 2 * ( NOLDTR + NNEWTR )
         LREDS  = 2 * ( NOLDTR + NNEWTR )
C
         DTIME      = SECOND()
         CALL RP_REDEV(NEXCI,NOLDTR,NNEWTR,ISYMTR,REDE,LREDE,REDS,LREDS,
     &                 LREDOL,EXVAL,LEXVAL,WORK(KEND1),LWORK1)
         DTIME      = SECOND()   - DTIME
         SOTIME(33) = SOTIME(33) + DTIME
C
         LREDOL = LREDE
C
C------------------------------------------
C        Reset number of old trial vectors.
C------------------------------------------
C
         NOLDTR = MIN(NOLDTR + NNEWTR, (NSAVMX - 1) * NEXCI )
C
C-------------------------------------------------------------------
C        Determine the residues from the current optimal solution 
C        vectors and decide if convergence is obtained for any of
C        the vectors. For the non-converged vectors create new
C        trial-vectors. These are orthogonalized against the 
C        previous optimal trial-vectors and among themself including
C        the vectors obtained from pairing.
C-------------------------------------------------------------------
C        
         IF (NIT .GE. MAXIT) NONEWT = .TRUE.
C
cLig  <> I have added the parameters needed for the new RP_TRIAL2
         DTIME      = SECOND()
         CALL RP_TRIAL2('EXCITA',NONEWT,NOLDTR,NNEWTR,NLINDP,EXVAL,
     &                  LEXVAL,RESINM,LRESINM,CONV,LCONV,NCONV,ISYMTR,
     &                  NEXCI,DUMMY,DUMMY,WORK(KEND1),LWORK1)
         DTIME      = SECOND()   - DTIME
         SOTIME(34) = SOTIME(34) + DTIME
C
C---------------------------------------------------------------------
C        Write calculated excitation energies and residuals to output.
C---------------------------------------------------------------------
C
         IF ( IPRSOP .GE. 2 ) THEN
C
            WRITE(LUPRI,9010)
            WRITE(LUPRI,9011)
            WRITE(LUPRI,9010)
            DO 200 IEXCI = 1,NEXCI
               WRITE(LUPRI,9012) IEXCI,EXVAL(IEXCI),RESINM(IEXCI),
     &                           CONV(IEXCI)
 200        CONTINUE
            WRITE(LUPRI,9010)
C
         END IF
C
C---------------------------------------
C        Flush the standard output file.
C---------------------------------------
C
         CALL FLSHFO(LUPRI)
C
C---------------------------------------------------------------------
C     Go to next iteration if all eigenvalues are not converged and if
C     the maximum number of iterations have not been reached.
C---------------------------------------------------------------------
C
      IF ( (NNEWTR .GT. 0) .AND. (NIT .LT. MAXIT) ) GO TO 100
C
      IF ( (NNEWTR .EQ. 0) .AND. (NIT.LT.MAXIT) ) THEN
C
         IF ( NLINDP .EQ. 0 ) THEN
C
            WRITE(LUPRI,9001)
            WRITE(LUPRI,9002)
            WRITE(LUPRI,9003)
            WRITE(LUPRI,9006)
            WRITE(LUPRI,9008)
C
         ELSE
C
            WRITE(LUPRI,9001)
            WRITE(LUPRI,9002)
            WRITE(LUPRI,9003)
            WRITE(LUPRI,9004) NLINDP
            WRITE(LUPRI,9005) NEXCI - NLINDP
            WRITE(LUPRI,9008)
C
         END IF
C
      ELSE IF (NIT .EQ. MAXIT) THEN
C
         WRITE(LUPRI,9001)
         WRITE(LUPRI,9002)
         WRITE(LUPRI,9003)
         WRITE(LUPRI,9007) MAXIT
         WRITE(LUPRI,9008)
C
         NNEWTR = 0
C
      ELSE 
C
         CALL QUIT('ERROR occured in RP_LRSOLV')
C
      END IF
C
      CALL FLSHFO(LUPRI)
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('RP_LRSOLV')
C
      RETURN
C
 9001 FORMAT(/,1X,'---------------------------------------------',
     &       '----------------')
 9002 FORMAT(1X,'Iterations stopped since: ')
 9003 FORMAT(1X,'-------------------------')
 9004 FORMAT(1X,I2,' of the excitations have linear dependent new',
     &       ' trial vectors.')
 9005 FORMAT(1X,I2,' of the excitations are converged.')
 9006 FORMAT(1X,'All of the excitations are converged.')
 9007 FORMAT(1X,'Maximum number of ',I3,' iterations is reached.')
 9008 FORMAT(1X,'---------------------------------------------',
     &       '----------------',/)
 9010 FORMAT(12X,'--------------------------------------------------')
 9011 FORMAT(12X,'Excitation   Energy (au)     Residual    Converged')
 9012 FORMAT(13X,I5,2X,F15.8,1X,1P,D14.4,6X,A)
      END
