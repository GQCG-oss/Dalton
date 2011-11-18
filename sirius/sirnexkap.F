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
C  /* Deck kapbet */
      SUBROUTINE KAPBET (NHORED,RTRUST,BETA,HORED,EVAL,EVEC,WRK,LFREE)
C
C Written 4-Mar-1985 by Hans Joergen Aa. Jensen
C (based on UPDBET from SIRNEO)
C Revised 890809 hjaaj: introduced RTTOL1, RTTOL2 as in UPDBET
C
C  This is the routine for update of beta in the reduced matrix
C  to give the desired step length (trust region control).
C
C  This version is only for ground state optimization.
C
C
#include <implicit.h>
C
      DIMENSION HORED(*),EVEC(NHORED,*),EVAL(*)
      PARAMETER (ITBMAX = 50)
      DIMENSION WRK(ITBMAX,4)
C
#include <iratdef.h>
      PARAMETER (ISTATE = 1,     BETMIN = 0.1D0, RTTOL = 1.1D0)
      PARAMETER (SECFAC = 0.1D0)
      PARAMETER (DHALF=0.5D0, D0=0.0D0, D1=1.0D0, D2=2.0D0)
      PARAMETER (DP1 = 0.1D0, DP9 = 0.9D0)
#include <dummy.h>
C
C
C Used from common blocks:
C  INFPRI: IPRKAP
C
#include <priunit.h>
#include <infpri.h>
C
C Local variable declarations:
C
      LOGICAL SECANT
C
      CALL QENTER('KAPBET')
C
      KEND   = 3*ITBMAX + NHORED*(NHORED+1)/2
      KEND   = KEND + NHORED   + (NHORED-1)/IRAT+1
      IF (KEND.GT.LFREE) CALL ERRWRK('KAPBET',KEND,LFREE)
C
C
      JPRKAP = IPRKAP
      NBOVEC = NHORED - 1
C
C --- Find optimal BETA.
C     Ground state section.
C
      BETL = MAX(BETMIN,ABS(BETA))
      IF (BETL .NE. BETA) THEN
         IF (IPRKAP .LE. 30) THEN
            IPRKAP = 0
         ELSE
            WRITE (LUPRI,'(/A,F12.5)')
     *      ' *** KAPBET calling KAPRED with trial BETA =',BETL
         END IF
         CALL KAPRED (.TRUE.,2,HORED,DUMMY,NBOVEC,0,BETL,
     *                EVEC,EVAL,DUMMY,DUMMY,DUMMY,WRK(1,4))
C        CALL KAPRED (DONEO,ICTL,REDO,REDGO,NBOVEC,N,DAMP,
C    *                EOVEC,EOVAL,GDM,BOVEC,SOVEC,WRK)
         IPRKAP = JPRKAP
      END IF
      SECANT = .FALSE.
      RTTOL1 = RTTOL
      RTTOL2 = D1 / (DP9 + DP1*RTTOL)
      FBI2   = D0
C
C     loop over BETA iterations ...
C
      DO 600 ITB = 1,ITBMAX
         JTB   = ITB
         RBETL = SQRT(EVEC(1,ISTATE) ** (-2) - D1) / BETL
         FBETL = RBETL - RTRUST
         WRK(ITB,1) = BETL
         WRK(ITB,2) = RBETL
         WRK(ITB,3) = EVAL(ISTATE)
C
         IF (RBETL.GT.RTTOL1*RTRUST) THEN
C        (step too long, increase beta)
            IF (SECANT) GO TO 200
            IF (FBI2*FBETL .LT. D0) THEN
               SECANT = .TRUE.
               GO TO 200
            END IF
            FBI2 = FBETL
            BI2  = BETL
            BETL = BETL * MIN(D2,(RBETL/RTRUST))
            IF (ITB .EQ. 1) BETL = MAX(BETL,D1)
            GO TO 300
         ELSE IF (RBETL.LT.RTTOL2*RTRUST) THEN
Chj-890809: if test for RBETL.LT.RTRUST then algorithm may not converge
C           as RBETL may approach RTRUST from below.
C        (step shorter then RTRUST, if NR-like
C        -- i.e. beta = betmin -- exit loop)
CHJ-S-850211: revised to set BETA = D1 when local.
            IF (BETL .EQ. BETMIN) THEN
               IF (BETMIN .LT. D1 .AND. RBETL .LT. DHALF*RTRUST) THEN
                  BETL = D1
                  GO TO 700
               ELSE
                  GO TO 800
               END IF
            END IF
CHJ-E-850211
            IF (SECANT) THEN
               GO TO 200
            ELSE
               IF (FBI2*FBETL .LT. D0) THEN
                  SECANT = .TRUE.
                  GO TO 200
               END IF
               IF (BETL .LE. D1) THEN
                  BETL = BETMIN
               ELSE
                  BETL = BETL * (RBETL/RTRUST)
               END IF
               GO TO 300
            END IF
         ELSE
C        (step within tolerance on RTRUST, exit loop)
            GO TO 800
         END IF
C
C Use the modified "Illinois method", a variation of the secant method,
C  to obtain next guess for BETA (Illinois: SECFAC = one half).
C  (the "Pegasus method", SECFAC = FBI2 / (FBI2+FBETL),
C   did not work as well)
C
  200    CONTINUE
         IF (FBI2*FBETL .LT. D0) THEN
            FBI1 = FBI2
            BI1  = BI2
         ELSE
            FBI1 = SECFAC*FBI1
         END IF
         FBI2 = FBETL
         BI2  = BETL
         BETL = (FBI2*BI1 - FBI1*BI2) / (FBI2 - FBI1)
C
C Solve for eigenvalues and -vectors with new beta
C
  300    CONTINUE
         BETL   = MAX(BETL,BETMIN)
         IPRKAP = 0
         CALL KAPRED (.TRUE.,2,HORED,DUMMY,NBOVEC,0,BETL,
     *                EVEC,EVAL,DUMMY,DUMMY,DUMMY,WRK(1,4))
C        CALL KAPRED (DONEO,ICTL,REDO,REDGO,NBOVEC,N,DAMP,
C    *                EOVEC,EOVAL,GDM,BOVEC,SOVEC,WRK)
         IPRKAP = JPRKAP
C     end loop ITB ...
  600 CONTINUE
      GO TO 900
C
C *** normal exit
C
  700 CONTINUE
      JTB = JTB + 1
      CALL KAPRED (.TRUE.,2,HORED,DUMMY,NBOVEC,0,BETL,
     *             EVEC,EVAL,DUMMY,DUMMY,DUMMY,WRK(1,4))
C     CALL KAPRED (DONEO,ICTL,REDO,REDGO,NBOVEC,N,DAMP,
C    *             EOVEC,EOVAL,GDM,BOVEC,SOVEC,WRK)
      RBETL = SQRT(EVEC(1,ISTATE) ** (-2) - D1) / BETL
      WRK(JTB,1) = BETL
      WRK(JTB,2) = RBETL
      WRK(JTB,3) = EVAL(ISTATE)
  800 CONTINUE
      IF (IPRKAP .GE. 3 .AND. BETA.NE.BETL) THEN
         WRITE (LUPRI,8010) BETL,JTB,NHORED
         IF (IPRKAP .GE. 5) THEN
            WRITE (LUPRI,8012)
            CALL OUTPUT (WRK,1,JTB,1,3,ITBMAX,3,1,LUPRI)
         END IF
      END IF
 8010 FORMAT (/' New beta =',F10.5,' found in',I3,' iterations.',
     *        /' Dimension of reduced augmented orbital Hessian  =',I3)
 8012 FORMAT (/' Found values of beta, step length, and ',
     *         'level shift for this matrix:')
      BETA = BETL
      CALL QEXIT('KAPBET')
      RETURN
C
C ***
C *** ERROR exit
C
  900 CONTINUE
      WRITE (LUPRI,9010) ITBMAX
      CALL OUTPUT (WRK,1,ITBMAX,1,3,ITBMAX,3,1,LUPRI)
      CALL QTRACE(LUPRI)
      CALL QUIT('No new BETA found in KAPBET')
 9010 FORMAT (/' KAPBET-ERROR, new beta not found in',I3,' iterations.')
C
C --  end of subroutine KAPBET
C
      END
C  /* Deck kapctl */
      SUBROUTINE KAPCTL(NRSEQ,MAXKAP,THRKAP,EOVAL,DAMP,NBOFI3,NBOVEC,
     *                  NNEW,ROTRST,REDO,REDGO,EOVEC,BOVEC,SOVEC,GDM,
     *                  DIAOR,CMO,GORB,DV,PV,FC,FV,WRK,LFREE)
C
C Written 29-Dec-1984 by Hans Joergen Aa. Jensen.
C
C Purpose:
C  Do "micro-micro" iterations to find next kappa, that is
C  an orbital rotation trial vector for next micro iteration.
C
C NRSEQ=1 SOLVE LEVEL SHIFTED NEWTON RAPHSON EQUATIONS
C         (LEVEL SHIFT DAMP) TO FIND IMPROVED ORBITAL
C          PARAMETERS
C NRSEQ=0 SOLVE LEVEL SHIFTED NEWTON RAPHSON EQUATIONS
C         AS AN EIGENVALUE PROBLEM. ADJUST LEVEL SHIFT
C         TO OBTAIN STEP LENGTH ROTRST (DAMP IS INITIAL
C         GUESS FOR BETA PARAMETER). Used for "orbital
C         absorption".
C
C Input:
C
C Output:
C
#include <implicit.h>
#include <infvar.h>
      DIMENSION REDO(*),REDGO(*),EOVEC(*),SOVEC(NWOPT,*),BOVEC(NWOPT,*)
      DIMENSION GDM(*),DIAOR(*), CMO(*),GORB(*),DV(*),PV(*),FC(*),FV(*)
      DIMENSION WRK(*), EOVAL(*)
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, MINMIC = 2)
#include <dummy.h>
C
C Used from common blocks:
C   INFVAR : NWOPT
C
#include <maxorb.h>
#include <priunit.h>
#include <infpri.h>
C
      LOGICAL DONEO
C
      CALL QENTER('KAPCTL')
C
      IF (NBOVEC .GE. MAXKAP) THEN
         WRITE (LUPRI,3010) NBOVEC, MAXKAP
         CALL QTRACE(LUPRI)
         CALL QUIT('KAPCTL: INSUFFICIENT SPACE FOR TRIAL VECTORS')
      END IF
      IF (NRSEQ .EQ. 0) THEN
         DONEO = .TRUE.
      ELSE IF (NRSEQ .EQ. 1) THEN
         DONEO = .FALSE.
      ELSE
         WRITE(LUPRI,3020) NRSEQ
         CALL QTRACE(LUPRI)
         CALL QUIT('KAPCTL: NRSEQ ILLEGAL VALUE')
      END IF
 3010 FORMAT(/' ***KAPCTL*** ERROR INSUFFICIENT SPACE FOR TRIAL VECTORS'
     *       /'              NBOVEC =',I5,' MAXKAP =',I5)
 3020 FORMAT(/' ***KAPCTL*** ERROR NRSEQ NOT SET TO 1 OR 0, WAS',I5)
C
C
C Step 1: allocate work space:
C
      KWRK = 1
      LWRK = LFREE
C
C Step 2: construct reduced matrix
C
      IF (DONEO) THEN
         IF (IPRKAP .GE. 10) WRITE(LUPRI,'(/A,F15.5)')
     *      ' (KAPCTL) OBETA in parameter list:',DAMP
         OBETA = MAX(DAMP,D1)
         IF (IPRKAP .GE. 2) WRITE (LUPRI,1000) OBETA,ROTRST,THRKAP
         DAMP  = OBETA
      END IF
 1000 FORMAT(/' (KAPCTL) NEO CALCULATION; OBETA used, ROTRST, THRKAP :',
     *       /(T11,1P,D15.5,0P,F8.4,1P,D10.2))
C
C SET UP INITIAL REDUCED ORBITAL HESSIAN AND FIND SOLUTION VECTOR
C
      CALL KAPRED(DONEO,3,REDO,REDGO,NBOVEC,NNEW,DAMP,
     *            EOVEC,EOVAL,GDM,BOVEC,SOVEC,WRK(KWRK))
C
C --- start of micro-micro iterations ---
C
      IF ( IPRKAP .EQ. 2 ) THEN
         WRITE (LUPRI,'(//A/A)')
     *   ' ORBITAL SUB-ITERATION       ORBITAL RESIDUAL',
     *   ' ---------------------       ----------------'
      END IF
      MICMIC = 0
  100 MICMIC = MICMIC + 1
         IF (IPRKAP .GE. 3) WRITE (LUPRI,'(//A,I5)')
     *      ' *** ORBITAL MICRO-MICRO ITERATION NO.',MICMIC
         IF (MICMIC .LT. MINMIC) THEN
C           do not check convergence in KAPNEX
            JOCONV = 1
         ELSE
            JOCONV = 0
         END IF
         CALL KAPNEX(MICMIC,DONEO,NBOVEC,DAMP,EOVAL,EOVEC,REDGO,
     *               THRKAP,JOCONV,GDM,DIAOR,BOVEC,SOVEC)
      IF (JOCONV.NE.1) GO TO 200
         NBOVEC = NBOVEC + 1
         CALL ORBLIN(1,BOVEC(1,NBOVEC),CMO,GORB,DV,PV,FC,FV,
     &               WRK(KWRK),1,LWRK)
C        CALL ORBLIN(NOSIM,BOVECS,CMO,GORB,DV,PV,FC,FV,
C    &               WRK,KFRSAV,LFRSAV)
         CALL DCOPY(NWOPT,WRK(KWRK),1,SOVEC(1,NBOVEC),1)
         CALL KAPRED(DONEO,3,REDO,REDGO,NBOVEC,1,DAMP,
     *               EOVEC,EOVAL,GDM,BOVEC,SOVEC,WRK(KWRK))
         IF (DONEO) THEN
            CALL KAPBET (NBOVEC+1,ROTRST,DAMP,REDO,EOVAL,EOVEC,
     *                   WRK(KWRK),LWRK)
C           CALL KAPBET (NHORED,RTRUST,BETA,HORED,EVAL,EVEC,WRK,LFREE)
         END IF
         CALL FLSHFO(LUPRI)
      IF ((NBOVEC+1) .LE. MAXKAP) GO TO 100
      WRITE (LUPRI,'(//A,I5/)')
     *   ' Orbital micro-micro iterations stopped'//
     *   ' because max. it. reached:',MICMIC
      GO TO 300
C
C --- end of micro-micro iterations ---
C
  200 CONTINUE
      IF (JOCONV.EQ.0) THEN
C930629-hjaaj: MICMIC is now printed in KAPNEX
C        IF (IPRKAP.GT.0) THEN
C           WRITE (LUPRI,'(/A,I5,A)')
C    *     ' Orbital micro-micro iterations converged in',
C    *     MICMIC,' iterations.'
C        END IF
      ELSE IF (JOCONV.EQ.-1) THEN
         WRITE (LUPRI,'(/A/A,I5)')
     *   ' Orbital micro-micro iterations stopped because of lin. dep.',
     *   ' Number of iterations done =',MICMIC
      END IF
C
C -- PUT KAPPA VECTORS IN GDM
C    921127-hjaaj: Components from IFIL3 are excluded
C    (these has already been added as trial vectors previously)
C
  300 IF (NBOVEC .GT. NBOFI3) THEN
         IF (DONEO) THEN
            IOFF=1
         ELSE
            IOFF=0
         END IF
         CALL DZERO(GDM,NWOPT)
         DO 500 I=NBOFI3+1,NBOVEC
            CALL DAXPY(NWOPT,EOVEC(IOFF+I),BOVEC(1,I),1,GDM,1)
  500    CONTINUE
         IF (DONEO) THEN
            FAC = D1 / (EOVEC(1)*DAMP)
            CALL DSCAL(NWOPT,FAC,GDM,1)
         END IF
         IF (IPRKAP.GE.35) THEN
            WRITE(LUPRI,'(/A)')' ***KAPCTL**** FINAL KAPPA VECTOR'
            CALL OUTPUT(GDM,1,NWOPT,1,1,NWOPT,1,1,LUPRI)
         END IF
      END IF
      CALL FLSHFO(LUPRI)
C
C *** End of subroutine KAPCTL
C
      CALL QEXIT('KAPCTL')
      RETURN
      END
C  /* Deck kapnex */
      SUBROUTINE KAPNEX (MICMIC,DONEO,NBOVEC,DAMP,EOVAL,EOVEC,REDGO,
     *                   THRKAP,JOCONV,GDM,DIAOR,BOVEC,SOVEC)
C
C 28-Dec-1984 PJ+HJAaJ
C (Based on NEXT)
C
C Purpose:
C  Find next guess(es) for kappa for a fixed set of configuration
C  coefficients in micro-micro iterations
C  and orthogonalize to previous b-vectors.
C
C Input:
C  IF (DONEO) THEN
C     SOLVE LEVEL SHIFTED NR ITERATIONS AS AN EIGENVALUE PROBLEM
C     DAMP , current restricted-step factor (BETA)
C     EOVAL, lowest eigenvalue of reduced HESSIAN-matrix
C     EOVEC, lowest eigenvector of reduced HESSIAN-matrix
C  ELSE IF (.NOT. DONEO) THEN
C     SOLVE AS LINEAR SET OF EQUATIONS
C     DAMP , level shift in NR iteration
C     EOVEC, solution to linear set of equations
C  END IF
C  THRKAP, threshold for convergence of the KAPPA vector
C  NBOVEC, number of B-vectors in BOVEC
C  BOVEC(*,1:NBOVEC), previous b-vectors, except first element.
C  SOVEC(*,1:NBOVEC), corresponding sigma vectors
C
C FOR DONEO:
C   We assume implicitly that the B-vectors are extended with
C   the vector (1,0) which are located in the beginning.
C   The corresponding S-vector is (0,OBETA*GDM)
C
C Input :
C  JOCONV =0  check convergence
C         =1  do not check convergence
C
C Output:
C  JOCONV =-1 if the new b-vector for the root we want to
C             converge to is linear dependent
C             with previous b-vectors
C         =1  if micro iterations not converged
C         =0  if micro iterations converged
C  BOVEC(*,NBOVEC+1) the new trial vectors
C
C
C
#include <implicit.h>
#include <infvar.h>
      DIMENSION EOVEC(*),GDM(*),REDGO(*),DIAOR(*)
      DIMENSION BOVEC(NWOPT,*),SOVEC(NWOPT,*)
      LOGICAL   DONEO
C
#include <thrldp.h>
      PARAMETER ( DTEST = 1.D-3, D0 = 0.0D0, D1 = 1.0D0 )
C
C Used from common blocks:
C  INFVAR: NWOPT,?
C
#include <maxorb.h>
#include <priunit.h>
#include <infpri.h>
C
C
      LOGICAL CHKCNV
C
      CALL QENTER('KAPNEX')
C
      CHKCNV = (JOCONV .EQ. 0)
C
      IF (IPRKAP.GE.9) WRITE(LUPRI,'(/A,1P,D15.5,I5,I5)')
     * ' (KAPNEX entered) DAMP,NBOVEC,JOCONV : ',DAMP,NBOVEC,JOCONV
C
      IF (DONEO) THEN
C     --- NEO EIGENVALUE EQUATION
         OBETA  = DAMP
         DSHIFT = EOVAL
C        DSHIFT: SHIFT IN DAVIDSON ALGORITHM
C
C        Find lowest diagonal value: IF NEGATIVE AND LESS THAN EOVAL
C        REDEFINE DSHIFT
C
         DLOW = D0
         IDLOW = 0
         DO 120 I = 1,NWOPT
            IF (DIAOR(I).LT.DLOW) THEN
               DLOW = DIAOR(I)
               IDLOW = I
            END IF
  120    CONTINUE
         IF (IDLOW.NE.0) THEN
            IF ( IPRKAP.GE.3 ) WRITE (LUPRI,1200) IDLOW,DLOW,EOVAL
            IF (EOVAL.GE.DLOW) THEN
               DSHIFT = DLOW*1.1D0
               IF ( IPRKAP.GE.3 ) WRITE (LUPRI,1210) DSHIFT
            END IF
         END IF
 1200    FORMAT(/' (KAPNEX) lowest orbital hessian diagonal element:',
     *           I8,T60,1P,D10.2,
     *          /T11,'lowest eigenvalue of reduced orbital hessian:',
     *           T60,1P,D10.2)
 1210    FORMAT(/T11,'Hence eigenvalues shifted to:',F20.15)
C
C        *******************************************************
C        Step B. B-VECTORS IN BOVEC
C                S-VECTORS IN SOVEC
C                EOVEC are eigenvectors of reduced O-matrix.
C
C        CONSTRUCT ORBITAL PART OF THE RESIDUAL VECTOR
C
C        (first sigma vector is OBETA dependent)
         FAC = EOVEC(1)*OBETA
         CALL DCOPY(NWOPT,GDM,1,BOVEC(1,NBOVEC+1),1)
         CALL DSCAL(NWOPT,FAC,  BOVEC(1,NBOVEC+1),1)
C        (the rest are not)
         DO 300 K = 1,NBOVEC
            FAC = EOVEC(K+1)
            CALL DAXPY(NWOPT,FAC,SOVEC(1,K),1,BOVEC(1,NBOVEC+1),1)
            FAC = -FAC*EOVAL
            CALL DAXPY(NWOPT,FAC,BOVEC(1,K),1,BOVEC(1,NBOVEC+1),1)
  300    CONTINUE
C
C
C        Step C. - convergence test.
C
C        CONTRIBUTION TO RESIDUAL FROM FIRST COMPONENT
C        = sum(i=0,NBOVEC) (S(0,i) - w * B(0,i)) * a(i)
C        = - w * B(0,0) * a(0) + obeta * sum(i=1,NBOVEC) REDGO(i) * a(i)
C        (REDGO(i) = ddot(S(0,*), B(i,*)), i = 1,NBOVEC)
C
C        AND ORBITAL PART CONTRIBUTION
C
         X = D0
         DO 400 K = 1,NBOVEC
            X = X + REDGO(K) * EOVEC(K+1)
  400    CONTINUE
         X = OBETA*X - EOVEC(1)*EOVAL
        IF (IPRKAP.GE.5.OR.ABS(X).GT.DTEST) WRITE(LUPRI,'(/A,1P,D15.8)')
     *     ' COMPONENT 1 CONTRIBUTION TO ORBITAL RESIDUAL (NEO)',X
         X = X*X
         QONORM = DNORM2(NWOPT,BOVEC(1,NBOVEC+1),1)
         QONORM = QONORM**2 + X
         QONORM = SQRT(QONORM) / ABS(EOVEC(1)*OBETA)
C        normalize residual norm so we get same value
C        as if N-R had been used.
C
      ELSE
C     --- N-R LINEAR EQUATIONS
C
         DSHIFT = DAMP
C        DSHIFT: LEVEL SHIFT IN CONJUGATE GRADIENT ALGORITHM
C
C        CONSTRUCT RESIDUAL VECTOR (A-B)*X(I)-DAMP*X(I)-F
C
         CALL DCOPY(NWOPT,GDM,1,BOVEC(1,NBOVEC+1),1)
         DO 1350 K = 1,NBOVEC
            CALL DAXPY(NWOPT,EOVEC(K),SOVEC(1,K),1,
     *                 BOVEC(1,NBOVEC+1),1)
            FAC = -DAMP*EOVEC(K)
            IF (FAC.NE.D0)
     *         CALL DAXPY(NWOPT,FAC,BOVEC(1,K),1,BOVEC(1,NBOVEC+1),1)
 1350    CONTINUE
         QONORM = DNORM2(NWOPT,BOVEC(1,NBOVEC+1),1)
      END IF
C
C     CONVERGENCE TEST (COMMON FOR NEO AND N-R)
C
      IF (CHKCNV .AND. QONORM.LE.THRKAP) GO TO 6000
      JOCONV = 1
      IF ( IPRKAP .GE. 3 ) THEN
         WRITE(LUPRI,'(/A,1P,D10.2)')
     *   ' ORBITAL RESIDUAL IN KAPNEX =',QONORM
      ELSE IF ( IPRKAP .EQ. 2 ) THEN
         WRITE (LUPRI,'(T10,I5,T33,1P,D10.2)') MICMIC,QONORM
      END IF
C
C Step D. CONSTRUCT DAVIDSON OR CONJUGATE GRAD VECTORS
C
      JR = NBOVEC + 1
      DO 1300 I = 1,NWOPT
         D = DIAOR(I) - DSHIFT
         IF (ABS(D) .LE. DTEST) THEN
            BOVEC(I,JR) = BOVEC(I,JR)/SIGN(DTEST,D)
         ELSE
            BOVEC(I,JR) = BOVEC(I,JR)/D
         END IF
 1300 CONTINUE
C
C Step E. Orthogonalization and normalization of new Davidson vectors
C         (AVERAG returns immediately if no averageing requested)
C
      CALL AVERAG(BOVEC,NWOPT,1)
      NNEW = 1
      THRLDV = NWOPT*THRLDP
      CALL ORTVEC(NBOVEC,NNEW,NWOPT,THRLDV,BOVEC,SOVEC(1,NBOVEC+1))
C     CALL ORTVEC(NOLD,NNEW,NVDIM,THRLDP,VEC,LINDEP)
      IF (NNEW.GT.0) GO TO 8000
C
C *** Linear dependence among b-vectors
C
      WRITE (LUPRI,5100) THRLDV
      JOCONV = -1
      GO TO 8000
C
 5100 FORMAT(/' (KAPNEX) diagonalization stops ',
     *        'due to linear dependence,'
     *       /T11,'norm squared of new BOVEC(*,1) is less than: ',
     *        1P,D15.8)
C
C Microiterations converged
C
 6000 JOCONV=0
      IF ( IPRKAP .EQ. 2 ) THEN
         WRITE (LUPRI,'(T10,I5,T33,1P,D10.2,A)')
     &      MICMIC,QONORM,' converged.'
      ELSE IF ( IPRKAP .GT. 0 ) THEN
         WRITE(LUPRI,6010) QONORM, MICMIC
      END IF
 6010 FORMAT(/' Orbital residual in KAPNEX converged to',1P,D10.2,
     &        ' in',I3,' iterations.')
C
C End of subroutine KAPNEX
C
 8000 CONTINUE
      CALL QEXIT('KAPNEX')
      RETURN
      END
C  /* Deck kapred */
      SUBROUTINE KAPRED (DONEO,ICTL,REDO,REDGO,NBOVEC,N,
     *                   DAMP,EOVEC,EOVAL,GDM,BOVEC,SOVEC,WRK)
C
C 28-Dec-1984 PJ+HJAaJ
C  (based on REDMAT)
C
C Input:
C  DONEO, FLOW CONTROL
C        FALSE  SOLVE LEVEL SHIFTED LINEAR SET OF EQUATIONS
C               (LEVEL SHIFT DAMP) TO FIND IMPROVED ORBITAL PARAMETERS
C        TRUE   SOLVE LEVEL SHIFTED LINEAR SET OF EQUATIONS
C               AS AN EIGENVALUE PROBLEM. ADJUST LEVEL SHIFT
C               TO OBTAIN STEP LENGTH RTRUST
C  ICTL,  FLOW CONTROL
C        =1  EXTEND THE REDUCED MATRIX
C        =2  FIND NEW SOLUTION
C        =3  BOTH
C  REDO,  the old reduced matrix
C  BOVEC, the b-vectors
C  SOVEC, the sigma vectors
C  NBOVEC, NUMBER OF B-VECTORS IN BOVEC
C         DONEO: NBOVEC+1 IS DIMENSION OF REDUCED HESSIAN
C                  THE VECTOR (1,0) IS ADDED AS THE FIRST VECTOR
C                  TO OBTAIN THE DIMENSION NBOVEC+1
C     NOT DONEO: NBOVEC IS DIMENSION OF REDUCED HESSIAN
C  N, number of new b-vectors
C  DAMP, LEVEL SHIFT IN LINEAR SET OF EQUATIONS
C        BETA IN NEO EIGENVALUE SYSTEM
C
C Output:
C  REDO,  the new, extended reduced orbital Hessian-matrix
C  REDGO, the new, extended reduced orbital gradient
C  EOVEC :.not.doneo SOLUTION TO LINEAR SET OF EQUATIONS
C        :     doneo LOWEST EIGENVECTOR (NEO)
C
#include <implicit.h>
#include <infvar.h>
      DIMENSION REDO(*),REDGO(*),EOVEC(*)
      DIMENSION GDM(*),BOVEC(NWOPT,*),SOVEC(NWOPT,*)
      DIMENSION WRK(*)
      LOGICAL DONEO
C
      PARAMETER (D0 = 0.0D0)
C
C Used from common blocks:
C   INFVAR : NWOPT
C   INFIND : IROW(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infind.h>
#include <infpri.h>
C
C
      CALL QENTER('KAPRED')
C
      IF (IPRKAP .GE. 9) WRITE (LUPRI,7555) DONEO,ICTL,NBOVEC,N
 7555 FORMAT(/' (KAPRED) CONTROL PARAMETERS DONEO =',L10,
     *       /T11,'ICTL =',I5,' NBOVEC =',I5,' N =',I5)
C
      IF (.NOT. DONEO) GO TO 2000
C
C *******************************************************
C
C SOLVE LINEAR SET OF EQUATIONS AS AN EIGENVALUE EQUATION
C
      IF (ICTL.EQ.2) GO TO 1000
C
C Section 1: extend reduced O-matrix with N new b-vectors
C
      IF (N.LT.1) GO TO 1000
C
C Add N new rows to REDO
C
      IF (NBOVEC+2 .GT. LIROW)
     *   CALL QUIT('KAPRED error: NBOVEC exceeds LIROW')
C        ... IROW(NBOVEC+2) used below
      REDO(1) = D0
      KST = NBOVEC-N+1
      DO 50 K = KST,NBOVEC
         KROW = IROW(K+1)
         REDGO(K) = DDOT(NWOPT,GDM,1,BOVEC(1,K),1)
         REDO(KROW+1) = REDGO(K)
         DO 40 L=1,K
            REDO(KROW+1+L) = DDOT(NWOPT,SOVEC(1,L),1,BOVEC(1,K),1)
   40    CONTINUE
   50 CONTINUE
C
 1000 CONTINUE
C
      IF(ICTL.EQ.1) GO TO 9999
C
C SOLVE THE EIGENVALUE PROBLEM
C
C Insert the only non-zero element in first column
C
      NREDO = NBOVEC + 1
      KTOT  = IROW(NREDO+1)
      KWJ   =   1 + KTOT
      KIWJ  = KWJ + NREDO
C
      CALL DCOPY(KTOT,REDO(1),1,WRK(1),1)
      DO 55 I=1,NREDO
         WRK(IROW(I)+1) = DAMP*WRK(IROW(I)+1)
 55   CONTINUE
      IF ( IPRKAP.GE.5 ) THEN
         WRITE (LUPRI,'(//A)') ' (KAPRED) NEO matrix:'
         CALL OUTPAK(WRK,NREDO,1,LUPRI)
      END IF
      CALL DUNIT(EOVEC,NREDO)
      CALL JACO(WRK,EOVEC,NREDO,NREDO,NREDO,
     *          WRK(KWJ),WRK(KIWJ))
      II = 0
      DO 150 I = 1,NREDO
         II = II + I
         WRK(I) = WRK(II)
 150  CONTINUE
      CALL ORDER (EOVEC,WRK,NREDO,NREDO)
      EOVAL = WRK(1)
      IF ( IPRKAP.GE.3 ) THEN
         WRITE(LUPRI,'(//A)') ' (KAPRED) NEO eigenvalues:'
         CALL OUTPUT(WRK,1,NREDO,1,1,NREDO,1,1,LUPRI)
      END IF
      IF ( IPRKAP.GE.4) THEN
         WRITE(LUPRI,'(/A)') ' - and NEO eigenvectors:'
         CALL OUTPUT(EOVEC,1,NREDO,1,NREDO,NREDO,NREDO,1,LUPRI)
      END IF
C
      GO TO 9999
C
C ********************************************
C
C SOLVE LEVEL SHIFTED NEWTON RAPHSON EQUATIONS
C
 2000 CONTINUE
      IF (ICTL.EQ.2 .OR. N.EQ.0) GO TO 3000
      KST=NBOVEC-N+1
      DO 60 K=KST,NBOVEC
         REDGO(K) = -DDOT(NWOPT,BOVEC(1,K),1,GDM,1)
         KROW = IROW(K)
         DO 70 L=1,K
            REDO(KROW+L) = DDOT(NWOPT,BOVEC(1,K),1,SOVEC(1,L),1)
 70      CONTINUE
 60   CONTINUE
C
 3000 CONTINUE
      IF (ICTL.EQ.1) GO TO 9999
C
C SOLVE LEVEL SHIFTED LINEAR SET OF EQUATIONS
C
C
      KTOT  = IROW(NBOVEC+1)
      KIPVT = 1 + KTOT
C
      CALL DCOPY(KTOT,REDO,1,WRK,1)
      II = 0
      DO 500 I = 1,NBOVEC
         II = II + I
         WRK(II) = WRK(II) - DAMP
 500  CONTINUE
      CALL DCOPY(NBOVEC,REDGO,1,EOVEC,1)
      IF ( IPRKAP.GE.4 ) THEN
         WRITE(LUPRI,'(/A)') ' (KAPRED) Reduced gradient'
         CALL OUTPUT(EOVEC,1,NBOVEC,1,1,NBOVEC,1,1,LUPRI)
         IF ( IPRKAP.GE.5 ) THEN
            WRITE(LUPRI,'(/A)') ' (KAPRED) Reduced orbital Hessian'
            CALL OUTPAK(REDO,NBOVEC,1,LUPRI)
         END IF
         WRITE(LUPRI,'(A,1P,D15.8)')
     *   ' (KAPRED) Reduced orbital Hessian with level shift =',DAMP
      END IF
      CALL DSPSOL(NBOVEC,1,WRK(1),EOVEC,WRK(KIPVT),INFO)
      IF ( IPRKAP.GE.3 ) THEN
         WRITE(LUPRI,'(/A)')
     *      ' (KAPRED) Solutions to N-R set of linear equations'
         CALL OUTPUT(EOVEC,1,NBOVEC,1,1,NBOVEC,1,1,LUPRI)
      END IF
      IF (INFO.NE.0) THEN
         WRITE(LUPRI,8500) INFO
         CALL QTRACE(LUPRI)
         CALL QUIT('KAPRED: no solution to linear equations')
      END IF
 8500 FORMAT(/' (KAPRED) Solution not obtained to linear equations'
     *       /T11,'Check if matrix is singular, LINPACK DSP INFO =',I3)
C
C *** End of subroutine KAPRED
C
 9999 CALL QEXIT('KAPRED')
      RETURN
      END
C  /* Deck kapst */
      SUBROUTINE KAPST(DONEO,NBOVEC,NNEW,REDO,REDGO,NBFIL3,IFIL3,IFIL5,
     *                 IBNDX,IGDM,GDM,CMO,GORB,DIAOR,DSHIFT,
     *                 DV,PV,FC,FV,BOVEC,SOVEC,WRK,LFREE)
C
C WRITTEN JAN-1985 BY PJ,HJJ,AND DLY
C
C  PURPOSE:
C ----------
C DETERMINE THE INITIAL NBOVEC+1 KAPPA VECTORS IN BOVEC(1,*)
C AND THE CORRESPONDING LINEAR TRANSFORMED VECTORS IN SOVEC(1,*)
C THE NORMALIZED MODIFIED GRADIENT VECTOR IS PLACED IN BOVEC(1,1)
C THE RESULTING NBOVEC VECTORS ARE CONSTRUCTED FROM LUFIL3 AND
C LUFIL5 WHEN IGDM=1 AND FROM PREVIOUS ITERATIONS OTHERWISE
C
#include <implicit.h>
#include <infvar.h>
      DIMENSION REDO(*),REDGO(*),IBNDX(*),GDM(*),DIAOR(*)
      DIMENSION CMO(*), GORB(*), DV(*),   PV(*), FC(*),FV(*)
      DIMENSION BOVEC(NWOPT,*),  SOVEC(NWOPT,*), WRK(*)
      LOGICAL   DONEO
#include <thrldp.h>
      PARAMETER ( DTEST = 1.D-3, D0 = 0.0D0, D1=1.D0, DM8=1.0D-8 )
#include <dummy.h>
C
#include <ibndxdef.h>
C
C Used from common blocks:
C
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infpri.h>
C
      CALL QENTER('KAPST ')
      IF (IPRKAP.GE.9) WRITE (LUPRI,'(/A,3I5)')
     *   ' **KAPST; NBOVEC, NBFIL3, IGDM:',NBOVEC,NBFIL3,IGDM
C
      IF (IPRKAP.GE.3) TIMTOT = SECOND()
      IF (IPRKAP.GE.4) TIM    = TIMTOT
      IF (IGDM.EQ.1) THEN
         NBOVEC = 0
      IF (NBFIL3 .GT. 1) THEN
C
C        CONSTRUCT INITIAL ORBITAL VECTORS FROM ORBITAL VECTORS
C        IN MICROITERATIONS
C
         REWIND IFIL5
         IF (NCONF.GT.1) READ(IFIL5)
         READ(IFIL5)
C
         REWIND IFIL3
         READ(IFIL3)
C
C        PLACE B VECTORS (KAPPAS) IN BOVEC AND
C        CORRESPONDING SIGMA VECTORS IN SOVEC
C
         DO 100 I=2,NBFIL3
            IF (IBNDX(I).EQ.JBONDX)THEN
               NBOVEC=NBOVEC+1
               CALL READT(IFIL3,NWOPT,BOVEC(1,NBOVEC))
               IF (NCONF.GT.1) READ(IFIL5)
               CALL READT(IFIL5,NWOPT,SOVEC(1,NBOVEC))
            ELSE
               READ(IFIL3)
               READ(IFIL5)
               READ(IFIL5)
            END IF
  100    CONTINUE
      END IF
      END IF
C
C     CONSTRUCT  GDM DIVIDED BY shifted HESSIAN DIAGONAL
C     IN BOVEC(*,NBOVEC+1)
C
      NNEW = 1
      JR = NBOVEC + 1
      DO 1300 I = 1,NWOPT
         D = DIAOR(I) - DSHIFT
         IF (ABS(D) .LE. DTEST) D = SIGN(DTEST,D)
         BOVEC(I,JR) = GDM(I) / D
 1300 CONTINUE
C     (AVERAG returns immediately if no averageing requested)
      CALL AVERAG(BOVEC(1,NBOVEC+1),NWOPT,NNEW)
      THRLDV = NWOPT*THRLDP
      CALL ORTVEC(NBOVEC,NNEW,NWOPT,THRLDV,BOVEC,SOVEC(1,NBOVEC+1))
C     CALL ORTVEC(NOLD,NNEW,NVDIM,THRLDP,VEC,LINDEP)
      IF (IPRKAP.GE.4) THEN
         TIMNEW = SECOND()
         TIM    = TIMNEW - TIM
         WRITE (LUPRI,'(/A,F10.2)') ' KAPST TIMING; STEP 1 (SETUP) ',TIM
         TIM    = TIMNEW
      END IF
C
C     CONSTRUCT SOVEC(*,NBOVEC+1)
C     NOTE: SPACE ALLOCATION IN NEXKAP ASSUMES NNEW = 1.
C
      IF (NNEW.GE.1) THEN
         CALL ORBLIN(NNEW,BOVEC(1,NBOVEC+1),CMO,GORB,DV,PV,FC,FV,
     &               WRK,1,LFREE)
C        CALL ORBLIN(NOSIM,BOVECS,CMO,GORB,DV,PV,FC,FV,
C    &               WRK,KFRSAV,LFRSAV)
         CALL DCOPY((NNEW*NWOPT),WRK,1,SOVEC(1,NBOVEC+1),1)
      END IF
      IF (IPRKAP.GE.4) THEN
         TIMNEW = SECOND()
         TIM    = TIMNEW - TIM
         WRITE (LUPRI,'(T16,A,F10.2)') 'STEP 2 (ORBLIN)',TIM
         TIM    = TIMNEW
      END IF
C
C     RESET VALUES FOR NUMBER OF NEW VECTORS AND
C     SET UP REDUCED GRADIENT FOR NEW TRIAL VECTORS WHEN NSIM GT 1
C
      IF (IGDM.LE.1) THEN
         NNEW   = NBOVEC + NNEW
         NBOVEC = NNEW
      ELSE
         IF (DONEO) THEN
           REDO(1) = D0
           I1      = 1
           DO 300 I = 1,NBOVEC
              I1 = I1 + I
              REDGO(I) = DDOT(NWOPT,GDM,1,BOVEC(1,I),1)
              REDO(I1) = REDGO(I)
 300       CONTINUE
         ELSE
           DO 200 I = 1,NBOVEC
              REDGO(I) = -DDOT(NWOPT,GDM,1,BOVEC(1,I),1)
 200       CONTINUE
         END IF
         NBOVEC = NBOVEC + NNEW
      END IF
C
C     End of KAPST
C
      CALL QEXIT('KAPST ')
      RETURN
      END
C  /* Deck nexkap */
      SUBROUTINE NEXKAP(NGDM,GDM,NBFIL3,IFIL3,IFIL5,IBNDX,
     *                  NRSEQ,ROTRST,THRKAP,EVAL,DAMP,DIAOR,
     *                  CMO,GORB,DV,PV,FC,FV,WRK,LFREE)
C
C MOTECC-90: The purpose of this module, NEXKAP, and the algorithms
C            used are described in Chapter 8 Section E.2 of MOTECC-90
C            "Optimal Orbital Trial Vectors" and Section E.5
C            "Intermediate Optimization of Orbitals for Fixed
C            Configuration Coefficients"
C
C
C NGDM    NUMBER OF ORBITAL VECTORS
C GDM     MODIFIED GRADIENT VECTORS;RETURN OPTIMAL KAPPA VECTORS
C NBFIL3  NUMBER OF VECTORS ON IT3
C IBNDX   VECTOR SPECIFYING THE CHARACTER OF VECTORS ON IT3
C NRSEQ=1 SOLVE LINEAR SET OF EQUATIONS
C      =0 SOLVE LINEAR SET OF EQUATIONS AS AN EIGENVALUE PROBLEM
C THRKAP  DESIRED CONVERGENCE OF RESIDUAL
C ROTRST  SQRT(RTRUST**2-CISTEP**2) WHERE CISTEP IS THE LENGTH
C         OF THE OPTIMAL CSF VECTOR IN THE REDUCED SPACE.USED
C         ONLY WHEN NRSEQ=0
C EVAL    EIGENVALUES OF REDUCED HESSIAN FOR SIMULTANEOUS ROOT
C         OPTIMIZATION.
C         FIXED LEVEL SHIFTS WHEN LINEAR SET OF EQUATIONS
C         ARE SOLVED.
C DAMP    INITAL DAMPING FACTOR FOR SIMULTANEOUS ROOT OPTIMIZATION
C         (MAX(DAMP,D1) WILL BE USED) .
C
#include <implicit.h>
#include <infvar.h>
      DIMENSION GDM(NWOPT,*),IBNDX(*),THRKAP(*),EVAL(*),DIAOR(*)
      DIMENSION CMO(*),GORB(*),DV(*),PV(*),FC(*),FV(*), 
     &          WRK(*)
C
C Global parameters:
C
#include <ibndxdef.h>
C
C Local parameters:
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
C
C Common blocks:
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <infdim.h>
#include <infpri.h>
C
      LOGICAL DONEO
C
      CALL QENTER('NEXKAP')
C
      DONEO  = (NRSEQ .EQ. 0)
      TIMKAP = SECOND()
      IF (IPRKAP .GE. 2) THEN
         WRITE (LUPRI,'(/A/A)')
     *      ' ORBITAL MICRO-MICRO ITERATION OUTPUT (NEXKAP)',
     *      ' ---------------------------------------------'
         IF (.NOT. DONEO) THEN
            WRITE (LUPRI,'(A)')
     *         ' LINEAR EQUATION FORMULATION USED (NR)'
            IF (IPRKAP.GE.3) WRITE (LUPRI,'(/A/A/,(I10,5X,1P,2D14.4))')
     *         ' VECTOR NO.       LEVEL SHIFT   THRESHOLD',
     *         ' ----------       -----------   ----------',
     *         (I,EVAL(I),THRKAP(I),I=1,NGDM)
         ELSE
            WRITE (LUPRI,'(A)') ' EIGENVALUE FORMULATION USED (NEO)'
            IF (IPRKAP.GE.3) WRITE (LUPRI,'(2(/A,1P,D12.4))')
     *         ' TRUST RADIUS         :',ROTRST,
     *         ' CONVERGENCE THRESHOLD:',THRKAP(1)
         END IF
      END IF
      IF ( IPRKAP.GE.25 ) THEN
         WRITE (LUPRI,'(/A,I5)') ' (NEXKAP TEST OUTPUT) NBFIL3',NBFIL3
         WRITE (LUPRI,'(/A)') ' (NEXKAP) Modified orbital gradient(s)'
         CALL OUTPUT (GDM,1,NWOPT,1,NGDM,NWOPT,NGDM,-1,LUPRI)
         WRITE (LUPRI,'(/A)') ' Orbital gradient and Hessian diagonal'
         WRITE (LUPRI,'(I10,1P,2D20.8)') (I,GORB(I),DIAOR(I),I=1,NWOPT)
      END IF
C
C ALLOCATION OF WORK SPACE
C ------------------------
C
C     DEFINE MAXKAP
C
C     Space in LINTRN for orbital Hessian: LAO + NOSIM*LCO
C        LAO : H2/UGORB, H2X/UDV, BUF,IBUF
C        LCO : SOVECS, DXCAO/FXC,DXVAO/FXV, UTBO,
C              (FXCAO,FXVAO,FSCR) or (FXQ)
C     We only consider NOSIM = 1 in this version. 870508/hjaaj
C
      LAO    = N2ORBX + N2ASHX
      LCO    = NWOPT  + 2*NNBAST + N2ORBT
     *       + MAX(2*NNBAST+NBASMA*(NBASMA+1),NASHT*NORBT)
      LACO   = LAO + LCO
      MAXKAP = 15
 25   CONTINUE
      LSPACE = IROW(MAXKAP+2)+MAXKAP+(MAXKAP+1)**2+2*NWOPT*MAXKAP
     *       + MAX(LACO,IROW(MAXKAP+2)+2*MAXKAP+150)
      IF ((LFREE-LSPACE) .GT. 1000) THEN
C        ... compare to 1000 to leave some extra space
C            for items which where not considered in LSPACE.
         MAXKAP = MAXKAP + 1
         IF (MAXKAP .LE. NWOPT .AND. MAXKAP .LE. (LIROW-2)) GO TO 25
      END IF
      MAXKAP = MAXKAP - 1
      IF ((MAXKAP.LE.15).AND.(MAXKAP.LT.NWOPT)) THEN
         WRITE (LUPRI,30) MAXKAP+1,-LSPACE
         CALL QTRACE(LUPRI)
         CALL QUIT('NEXKAP: insufficient work space')
      END IF
 30   FORMAT(/' ***NOT ENOUGH WORK SPACE IN NEXKAP TO CARRY OUT',
     *        ' MICRO-MICRO ITERATIONS',
     *        /,T12,'MAXKAP=',I5,', ADDITIONAL SPACE NEEDED =',I8)
C
      IF ( IPRKAP.GE.7 ) WRITE(LUPRI,'(A,I6)') ' MAXKAP =',MAXKAP
C
C WORK SPACE CAN NOW BE ALLOCATED
C
      KREDO  = 1
      KREDGO = KREDO  + IROW(MAXKAP + 2)
      KEOVEC = KREDGO + MAXKAP
      KBOVEC = KEOVEC + (MAXKAP + 1)**2
      KSOVEC = KBOVEC + MAXKAP*NWOPT
      KWRK1  = KSOVEC + MAXKAP*NWOPT
      LWRK1  = LFREE  - KWRK1
C
C     FIND NUMBER OF ORBITAL VECTORS TO CHECK WORK SPACE
C
      NBOFI3 = 0
      DO 50 I = 1,NBFIL3
         IF (IBNDX(I).EQ.JBONDX) NBOFI3 = NBOFI3 + 1
 50   CONTINUE
      IF (NBOFI3.GT.(MAXKAP-10) .AND. (MAXKAP.LT.NWOPT)) THEN
         WRITE(LUPRI,55) NBOFI3,MAXKAP,NWOPT
         CALL QTRACE(LUPRI)
         CALL QUIT('NEXKAP: insufficient space for trial vectors')
      END IF
 55   FORMAT(/' ***NEXKAP NOT ENOUGH WORK SPACE FOR 10 TRIAL VECTORS'
     *       /T12,'NBOFI3=',I5,' MAXKAP=',I5,' NWOPT=',I5)
C
C     FIND NGDM OPTIMAL ORBITAL ROTATION PARAMETERS
C
      NBOVEC = 0
      TIMST  = D0
      TIMCTL = D0
      MST    = 0
      MCTL   = 0
      DO 100 IGDM=1,NGDM
         IF (DONEO) THEN
            DSHIFT = D0
         ELSE
            DAMP   = EVAL(IGDM)
            DSHIFT = DAMP
         END IF
C
C        DETERMINE INITIAL SET OF START VECTORS
C
         IF ( IPRKAP.GT.0 ) WRITE(LUPRI,'(/A,I3,A,I3/A)')
     *  ' MICRO-MICRO ITERATIONS ON ORBITAL VECTOR',IGDM,' OUT OF',NGDM,
     *  ' ====================================================='
C
         MST   = MST + 1
         TIMST = TIMST - SECOND()
         CALL KAPST(DONEO,NBOVEC,NNEW,WRK(KREDO),WRK(KREDGO),
     *              NBFIL3,IFIL3,IFIL5,IBNDX,IGDM,
     *              GDM(1,IGDM),CMO,GORB,DIAOR,DSHIFT,DV,PV,FC,FV,
     *              WRK(KBOVEC),WRK(KSOVEC),WRK(KWRK1),LWRK1)
C        CALL KAPST(DONEO,NBOVEC,NNEW,REDO,REDGO,NBFIL3,IFIL3,IFIL5,
C    *              IBNDX,IGDM,GDM,CMO,GORB,DIAOR,DSHIFT,
C    *              DV,PV,FC,FV,BOVEC,SOVEC,WRK,LFREE)
         TIMST = TIMST + SECOND()
C
C
C        DETERMINE OPTIMAL ROTATION PARAMETERS
C
C
         MCTL   = MCTL + 1
         TIMCTL = TIMCTL - SECOND()
         CALL KAPCTL(NRSEQ,MAXKAP,THRKAP(IGDM),EVAL(IGDM),DAMP,
     *               NBOFI3,NBOVEC,NNEW,ROTRST,WRK(KREDO),
     *               WRK(KREDGO),WRK(KEOVEC),
     *               WRK(KBOVEC),WRK(KSOVEC),GDM(1,IGDM),DIAOR,
     *               CMO,GORB,DV,PV,FC,FV,WRK(KWRK1),LWRK1)
         TIMCTL = TIMCTL + SECOND()
C
C
C
      IF ( IPRKAP.GE.15 ) THEN
         STEST = 1.D-14
         HTEST = 1.D-12
         WRITE (LUPRI,'(//,5(A/),2(A,1P,D10.2),A/)')
     *   ' (NEXKAP TEST OUTPUT)',
     *   ' SYMMETRY TEST OF BOVEC AND SOVEC, AND COMPARISON WITH REDO',
     *   ' TYPE 1: NOT SYMMETRIC (TYPE, I, J, HIJ, HJI,      DIFF)',
     *   ' TYPE 2: NOT EQ. REDO  (TYPE, I, J, HIJ, REDO(IJ), DIFF)',
     *   ' TYPE 3: B-MAT NOT UNIT(TYPE, I, J, SIJ                )',
     *   ' test:',HTEST,' for types 1-2 and',STEST,' for type 3'
         IF (DONEO) WRITE (LUPRI,'(A/)') ' (Add 1 to I and J.)'
         IBAD  = 0
         ITEST = 0
         IREDO = KREDO - 1
         IF (DONEO) IREDO = IREDO + 2
         IBV   = KBOVEC
         ISV   = KSOVEC
         DO 65 I = 1,NBOVEC
            JBV = KBOVEC
            JSV = KSOVEC
            DO 60 J = 1,I
               IREDO = IREDO + 1
               ITEST = ITEST + 1
               SIJ = DDOT(NWOPT,WRK(IBV),1,WRK(JBV),1)
               HIJ = DDOT(NWOPT,WRK(IBV),1,WRK(JSV),1)
               IF (J .NE. I) THEN
                  HJI = DDOT(NWOPT,WRK(ISV),1,WRK(JBV),1)
                  DIFF1 = HIJ - HJI
                  IF (ABS(DIFF1) .GT. HTEST) THEN
                     IBAD = IBAD + 1
                     WRITE (LUPRI,'(3I5,1P,3D20.8)') 1,I,J,HIJ,HJI,DIFF1
                  END IF
                  IF (ABS(SIJ) .GT. STEST) THEN
                     IBAD = IBAD + 1
                     WRITE (LUPRI,'(3I5,1P,D20.8)') 3,I,J,SIJ
                  END IF
               ELSE
                  IF (ABS(D1-SIJ) .GT. STEST) THEN
                     IBAD = IBAD + 1
                     WRITE (LUPRI,'(3I5,1P,D20.8)') 3,I,J,SIJ
                  END IF
               END IF
               DIFF2 = WRK(IREDO) - HIJ
               IF (ABS(DIFF2) .GT. HTEST) THEN
                  IBAD = IBAD + 1
                  WRITE (LUPRI,'(3I5,1P,3D20.8)')
     *               2,I,J,HIJ,WRK(IREDO),DIFF2
               END IF
               JBV = JBV + NWOPT
               JSV = JSV + NWOPT
   60       CONTINUE
            IF (DONEO) IREDO = IREDO + 1
            IBV = IBV + NWOPT
            ISV = ISV + NWOPT
   65    CONTINUE
         IF (IBAD.EQ.0) THEN
            WRITE (LUPRI,'(/I5,A)')
     *      ITEST,' elements were tested, all were OK'
         ELSE
            NWARN = NWARN + 1
            WRITE (LUPRI,'(/A/I5,A)')
     *      ' WARNING, errors found in NEXKAP tests',
     *      ITEST,' elements were tested'
            WRITE (LUERR,'(/A)')
     *      ' WARNING, errors found in NEXKAP tests. See output file.'
         END IF
      END IF
C
C     DECIDE NO. OF TRIAL ORBITAL VECTORS FOR NEXT ORBITAL OPTIMIZATION
C
         IF (NBOVEC.GT.(MAXKAP-15)) NBOVEC=MAXKAP-15
 100  CONTINUE
C
C
C
      IF ( IPRKAP.GE.24 ) THEN
         WRITE(LUPRI,'(//A)') ' (NEXKAP) Optimal orbital'//
     *    ' trial vectors from micro-micro iterations'
         IF (NBOFI3 .GT. 0) WRITE(LUPRI,'(A,I5,A)')
     *    '    (orthogonalized against previous',NBOFI3,
     *    ' micro iteration orbital trial vectors)'
         DO 200 IGDM = 1,NGDM
            WRITE(LUPRI,'(/A,I3)') ' OPTIMAL KAPPA VECTOR NO.',IGDM
            WRITE(LUPRI,'(I10,1P,D20.8)') (I,GDM(I,IGDM),I=1,NWOPT)
  200    CONTINUE
      END IF
      IF ( IPRKAP.GE.2) THEN
         TIMKAP = SECOND() - TIMKAP
         WRITE(LUPRI,'(//A,T40,F10.2,2(/T24,A,T40,F10.2,A,I3,A))')
     *      ' >>> TIMING IN NEXKAP, TOTAL TIME',TIMKAP,
     *      'TIME IN KAPST', TIMST, ' FOR',MST, ' CALLS',
     *      'TIME IN KAPCTL',TIMCTL,' FOR',MCTL,' CALLS'
      END IF
C
C END OF NEXKAP
C
      CALL QEXIT('NEXKAP')
      RETURN
      END
C  /* Deck ortvec */
      SUBROUTINE ORTVEC(NOLD,NNEW,NVDIM,THRLDP,VEC,LINDEP)
C
C 15-Mar-1985 hjaaj
C l.r. 4-May-1994 hjaaj (only elim. new vector if norm < THRLDP)
C
#include <implicit.h>
      DIMENSION VEC(NVDIM,*)
      LOGICAL   LINDEP(*)
C
      PARAMETER (THRRND=1.D-2, THRTT=1.D-4, D1=1.0D0)
C
C     Used from common blocks:
C
#include <priunit.h>
#include <infpri.h>
C
      CALL QENTER('ORTVEC')
C
      IF (NNEW.LE.0) GO TO 9999
      TLINDP = SQRT(THRLDP)
C
C     Normalize NNEW new vectors VEC(*,NOLD+1) TO VEC(*,NOLD+NNEW)
C
      IVEC = NOLD
      DO 200 INEW = 1,NNEW
         IVEC = IVEC + 1
         TT   = DNORM2(NVDIM,VEC(1,IVEC),1)
         IF (TT.LE.THRLDP) THEN
            LINDEP(INEW) = .TRUE.
            WRITE (LUPRI,8100) INEW,TT
         ELSE
            LINDEP(INEW) = .FALSE.
            IF (TT.LT.THRTT) THEN
               CALL DSCAL (NVDIM,(D1/TT),VEC(1,IVEC),1)
               TT = DNORM2(NVDIM,VEC(1,IVEC),1)
            END IF
            CALL DSCAL (NVDIM,(D1/TT),VEC(1,IVEC),1)
         END IF
  200 CONTINUE
C
C
C
      IROUND = 0
      ITURN  = 0
 1500 ITURN  = ITURN + 1
C
C     Orthogonalize new vectors against previous vectors
C
      DO 2000 K=1,NOLD
         DO 1900 J = NOLD+1,NOLD+NNEW
            TT = - DDOT(NVDIM,VEC(1,K),1,VEC(1,J),1)
            CALL DAXPY(NVDIM,TT,VEC(1,K),1,VEC(1,J),1)
 1900    CONTINUE
 2000 CONTINUE
C
C     Orthogonalize new vectors against each other
C     and normalization.
C
      DO 2400 INEW = 1,NNEW
         IF (.NOT.LINDEP(INEW)) THEN
C           ... orthogonalize using prev. vectors are normalized
            IVEC = NOLD + INEW
            DO 2300 JNEW = 1,(INEW-1)
               IF (.NOT.LINDEP(JNEW)) THEN
                  JVEC = NOLD + JNEW
                  TT = - DDOT(NVDIM,VEC(1,JVEC),1,VEC(1,IVEC),1)
                  CALL DAXPY(NVDIM,TT,VEC(1,JVEC),1,VEC(1,IVEC),1)
               END IF
 2300       CONTINUE
            TT = DNORM2(NVDIM,VEC(1,IVEC),1)
            IF (TT .LE. THRLDP) THEN
               LINDEP(INEW) = .TRUE.
               WRITE (LUPRI,8100) INEW,TT
            ELSE
               IF (TT .LT. THRRND) IROUND = IROUND+1
               IF (TT .LT. THRTT) THEN
                  CALL DSCAL(NVDIM,(D1/TT),VEC(1,IVEC),1)
                  TT = DNORM2(NVDIM,VEC(1,IVEC),1)
               END IF
               CALL DSCAL(NVDIM,(D1/TT),VEC(1,IVEC),1)
            END IF
         END IF
 2400 CONTINUE
C
C
      IF (IROUND.GT.0 .AND. ITURN.EQ.1) GO TO 1500
C
C
      JNEW = 0
      DO 4400 INEW = 1,NNEW
         IF (.NOT.LINDEP(INEW)) THEN
            JNEW = JNEW + 1
            IF (JNEW .LT. INEW) THEN
               CALL DCOPY (NVDIM,VEC(1,NOLD+INEW),1,VEC(1,NOLD+JNEW),1)
            END IF
         END IF
 4400 CONTINUE
      IF (JNEW .LT. NNEW) THEN
         WRITE (LUPRI,8200) NNEW,JNEW
         NNEW = JNEW
      END IF
C
C
C
 8100 FORMAT(/' ORTVEC: New vector no.',I3,
     *        ' is removed because of linear dependence;'
     *       /' norm of vector after Gram-Schmidt''s orthogonalization',
     *        1P,D9.2)
 8200 FORMAT(/' (ORTVEC) NNEW reduced from',I3,' to',I3)
C
C *** End of subroutine ORTVEC
C
 9999 CALL QEXIT('ORTVEC')
      RETURN
      END
