C
C  /* Deck so_lrsoeq */
      SUBROUTINE SO_LRSOEQ(MODEL,LABEL,ISYMTR,NEXCI,MAXIT,FRVAL,NFRVAL,
     &                     RESINM,LRESINM,CONV,LCONV,DENSIJ,LDENSIJ,
     &                     DENSAB,LDENSAB,DENSAI,LDENSAI,T2MP,LT2MP,
     &                     FOCKD,LFOCKD,REDE,REDS,LMXRED,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Stephan Sauer and Keld Bak, May 1996
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C     Andrea Ligabue, December 2003: linear response functions 
C                                    implemented
C
C     PURPOSE: Solve the SOPPA linear response equations 
C              using an AO-driven algorithm. 
C
C     LABEL             property label
C     ISYMTR            property symmetry
C     FRVAL(NFRVAL)     frequencies
C
#include "implicit.h"
#include "priunit.h"
C
#include "soppinf.h"
#include "ccsdsym.h"
C
      LOGICAL   NONEWT
C
      CHARACTER*3 CONV(LCONV)
      CHARACTER MODEL*5,LABEL*8
C
      DIMENSION FRVAL(NFRVAL),   RESINM(LRESINM)
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB), DENSAI(LDENSAI)
      DIMENSION T2MP(LT2MP),     FOCKD(LFOCKD),   WORK(LWORK)
      DIMENSION REDE(LMXRED),    REDS(LMXRED)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_LRSOEQ')
C
C==============================================================
C     For checking, calculate E[2] and S[2] matrices explicitly
C     by using unit vectors as trial vectors.
C==============================================================
C
      IF (SOPCHK) THEN
C
         CALL SO_CHECK(DENSIJ,LDENSIJ,DENSAB,LDENSAB,T2MP,LT2MP,
     &                 FOCKD,LFOCKD,ISYMTR,WORK,LWORK)
C
      END IF
C
C-----------------------------------------------------------------
C     Calculate diagonal parts of E[2] and S[2] and write to disk.
C-----------------------------------------------------------------
C
      DTIME      = SECOND()
      CALL SO_DIAG(FOCKD,LFOCKD,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &             ISYMTR,WORK,LWORK)
      DTIME      = SECOND()   - DTIME
      SOTIME(31) = SOTIME(31) + DTIME
C
C-------------------------
C     1. Memory allocation
C-------------------------
C
      KREDC   = 1
      KEND1   = (2*NSAVMX)**2
      LWORK1  = LWORK - KEND1
C
      CALL SO_MEMMAX('SO_LRSOEQ.1',LWORK1)
      IF(LWORK1.LT.0) CALL STOPIT('SO_LRSOEQ.1','  ',KEND1,LWORK)
C
C------------------------------
C     Loop over the frequencies
C------------------------------
C
      DO 300 IFREQ = 1,NFRVAL
C
C------------------------------------------------------------------
C     Set up initial trial vectors and write to disk using eq. (19)
C     with (C +-C) as R or the last TV of the previous frequency
C------------------------------------------------------------------
C
      DTIME      = SECOND()
      CALL SO_TRIAL3(NNEWTR,ISYMTR,IFREQ,FRVAL,NFRVAL,
     &               NEXCI,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &               WORK(KEND1),LWORK1)
      DTIME      = SECOND()   - DTIME
cLig  <> we have to define and write properly the index insotime array
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
C     At the moment SO_ERES will allways calculate DENSAI
C     when NIT == 1, it has no sense of "IFREQ". Zeroing it 
C     here will compensate
      IF( MODEL .EQ. 'AOSOP' ) 
     &   CALL DZERO(DENSAI, LDENSAI)
C
C------------------
C     Write banner.
C------------------
C 
      WRITE(LUPRI,'(/,2X,A)') '*********************************'//
     &                        '*********************************'
      IF (MODEL .EQ. 'AOSOP') THEN
         WRITE(LUPRI,'(14X,3A,I1)') 
     &        'SOPPA iterations, Property ',LABEL,', frequency ', IFREQ
      ELSE IF (MODEL .EQ. 'AOSOC') THEN
         WRITE(LUPRI,'(11X,A,I1)') 
     &        'SOPPA(CCSD)  iterations, Excitation symmetry ',ISYMTR
      ELSE
         WRITE(LUPRI,'(A)') 'ERROR: Wrong model in SO_LRSOEQ'
         CALL QUIT('ERROR: Wrong model in SO_LRSOEQ')
      END IF
      WRITE(LUPRI,'(2X,A)') '*********************************'//
     &                        '*********************************'
C
C----------------------------------------------------------------------
C     Iteration loop for solving the linear equation/eigenvalue problem
C----------------------------------------------------------------------
C
  100 CONTINUE
C
C--------------------------------------------------------------
C        Count number of iterations and write header to output.
C--------------------------------------------------------------
C
         NIT = NIT + 1
C
         IF ( IPRSOP .GE. 2 ) THEN
C
            WRITE(LUPRI,'(/,2X,A)') '================================'//
     &                              '=================================='
            IF (MODEL .EQ. 'AOSOP') THEN
               WRITE(LUPRI,'(14X,I3,A,I1)') NIT,
     &              '. SOPPA iteration, Frequency ',IFREQ
            ELSE IF (MODEL .EQ. 'AOSOC') THEN
               WRITE(LUPRI,'(11X,I3,A,I1)') NIT,
     &              '. SOPPA(CCSD)  iteration, Excitation symmetry ',
     &              ISYMTR
            ELSE
               WRITE(LUPRI,'(A)') 'ERROR: Wrong model in SO_LRSOEQ'
               CALL QUIT('ERROR: Wrong model in SO_LRSOEQ')
            END IF
            WRITE(LUPRI,'(2X,A,/)') '================================'//
     &                              '=================================='
C
         END IF
C
C--------------------------------------------------------------
C        Make E[2] linear transformation of trialvectors giving 
C        resultvectors.
C--------------------------------------------------------------
C
         CALL GETTIM (DUMMY,WTIMES)
         DTIME      = SECOND()
         CALL SO_ERES(MODEL,NOLDTR,NNEWTR,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                T2MP,LT2MP,FOCKD,LFOCKD,DENSAI,LDENSAI,NIT,ISYMTR,
     &                WORK(KEND1),LWORK1)
         DTIME      = SECOND()   - DTIME
         SOTIME(35) = SOTIME(35) + DTIME
         CALL GETTIM (DUMMY,WTIMEE)
         SOWTIM(1)  = SOWTIM(1)  + WTIMEE - WTIMES
C
C--------------------------------------------------------------
C        Make S[2] linear transformation of trialvectors giving
C        resultvectors.
C--------------------------------------------------------------
C
         CALL GETTIM (DUMMY,WTIMES)
         DTIME      = SECOND()
         CALL SO_SRES(NOLDTR,NNEWTR,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                ISYMTR,WORK(KEND1),LWORK1)
         DTIME      = SECOND()   - DTIME
         SOTIME(40) = SOTIME(40) + DTIME
         CALL GETTIM (DUMMY,WTIMEE)
         SOWTIM(2)  = SOWTIM(2)  + WTIMEE - WTIMES
C
         IF ( AOTEST ) THEN
C
C--------------------------------------------------------------------
C           Test orthonormality of trial vectors and check the linear
C           transformed trial vectors.
C--------------------------------------------------------------------
C
            CALL SO_TEST(NOLDTR,NNEWTR,ISYMTR,DENSIJ,LDENSIJ,DENSAB,
     &                   LDENSAB,T2MP,LT2MP,FOCKD,LFOCKD,WORK(KEND1),
     &                   LWORK1 )
C
         END IF
C
         IF((NIT.EQ.1).AND.(IFREQ.EQ.1).AND.AOSOP) THEN
C
C-----------------------------------------------------------------------
C          Calculate the gradient property vectors with the right DENSAI
C-----------------------------------------------------------------------
C
           LGPVC1  = 2 * NT1AM(ISYMTR)
CPFP   2009-Jun
           LGPVC2  = 2 * NT2AM(ISYMTR)
C           LGPVC2  = 2 * N2P2HOP(ISYMTR)
Cend-PFP
           KGPVC1  = KEND1
           KGPVC2  = KGPVC1 + LGPVC1 
           KEND2   = KGPVC2 + LGPVC2 
           LWORK2  = LWORK  - KEND2
C
           CALL SO_MEMMAX ('SO_LRSOEQ.2',LWORK2)
           IF (LWORK2 .LT.0) CALL STOPIT('SO_LRSOEQ.2',' ',KEND2,LWORK)
C
           CALL SO_GETGP(WORK(KGPVC1),LGPVC1,WORK(KGPVC2),LGPVC2,LABEL,
     &                  ISYMTR,T2MP,LT2MP,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                  DENSAI,LDENSAI,WORK(KEND2),LWORK2)
C
C---------------------------------------------------------
C          Save the complete property gradients on a file.
C---------------------------------------------------------
C
           CALL SO_OPEN(LUGPV1,FNGPV1,LGPVC1)
           CALL SO_OPEN(LUGPV2,FNGPV2,LGPVC2)
C
           CALL SO_WRITE(WORK(KGPVC1),LGPVC1,LUGPV1,FNGPV1,1)
           CALL SO_WRITE(WORK(KGPVC2),LGPVC2,LUGPV2,FNGPV2,1)
C
           IF(IPRSOP.GE.5) THEN
              LGPVC1H = LGPVC1 /2
              LGPVC2H = LGPVC2 /2
              CALL AROUND('In SO_LRSOEQ:  '//MODEL//
     &                    ' gradient property vector.'//LABEL)
              CALL OUTPUT(WORK(KGPVC1),1,LGPVC1H,1,2,LGPVC1H,2,1,LUPRI)
              CALL OUTPUT(WORK(KGPVC2),1,LGPVC2H,1,2,LGPVC2H,2,1,LUPRI)
           ENDIF
C     
           CALL SO_CLOSE(LUGPV1,FNGPV1,'KEEP')
           CALL SO_CLOSE(LUGPV2,FNGPV2,'KEEP')
C
         ENDIF
C
C----------------------------------------------------
C        Set up and solve the reduced linear problem.
C----------------------------------------------------
C
         LREDE  = 2 * ( NOLDTR + NNEWTR )
         LREDS  = 2 * ( NOLDTR + NNEWTR )
         LREDC  = 2 * ( NOLDTR + NNEWTR )
C
         CALL GETTIM (DUMMY,WTIMES)
         DTIME      = SECOND()
         CALL SO_REDLE(NEXCI,NOLDTR,NNEWTR,LABEL,ISYMTR,REDE,LREDE,
     &                 REDS,LREDS,WORK(KREDC),LREDC,LREDOL,FRVAL,NFRVAL,
     &                 IFREQ,ENORM,PROP,WORK(KEND1),LWORK1)
         DTIME      = SECOND()   - DTIME
         SOTIME(33) = SOTIME(33) + DTIME
         CALL GETTIM (DUMMY,WTIMEE)
         SOWTIM(3)  = SOWTIM(3)  + WTIMEE - WTIMES
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
         CALL GETTIM (DUMMY,WTIMES)
         DTIME      = SECOND()
cLig <> I have to add some parameters to that call and modify so_trial2
         CALL SO_TRIAL2('LINEAR',NONEWT,NOLDTR,NNEWTR,NLINDP,
     &                  FRVAL(IFREQ),1,RESINM,LRESINM,CONV,LCONV,NCONV,
     &                  ISYMTR,NEXCI,DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                  ENORM,LABEL,WORK(KEND1),LWORK1)
         DTIME      = SECOND()   - DTIME
         SOTIME(34) = SOTIME(34) + DTIME
         CALL GETTIM (DUMMY,WTIMEE)
         SOWTIM(4)  = SOWTIM(4)  + WTIMEE - WTIMES
C
C---------------------------------------------------------------------
C        Write calculated excitation energies and residuals to output.
C---------------------------------------------------------------------
C
         IF (IPRSOP .GE. 2) THEN
C
            WRITE(LUPRI,9010)
            WRITE(LUPRI,9011)
            WRITE(LUPRI,9010)
            WRITE(LUPRI,9012) LABEL,FRVAL(IFREQ),PROP,RESINM(1),CONV(1)
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
      IF ((NNEWTR .EQ. 0) .AND. (NIT.LT.MAXIT)) THEN
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
         CALL QUIT('ERROR occured in SO_LRSOEQ')
C
      END IF
C
C----------------------------------------
C    End of the loop over the frequencies
C----------------------------------------
C
  300 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('SO_LRSOEQ')
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
 9006 FORMAT(1X,'All of the properties are converged.')
 9007 FORMAT(1X,'Maximum number of ',I3,' iterations is reached.')
 9008 FORMAT(1X,'---------------------------------------------',
     &       '----------------',/)
 9010 FORMAT(12X,'--------------------------------------------------')
 9011 FORMAT(4X,'Label      Freq        Property (au)   Residual',
     &           '      Converged')
 9012 FORMAT(4X,A8,2X,F8.6,2X,F15.8,1X,1P,D14.4,6X,A)
      END
