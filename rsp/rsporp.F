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
#ifdef UNDEF
===========================================================================
/* Comdeck rev_log */
920929-ov    ORPSVE:  work dimension check
920723-hjaaj ORPLIN: inserted missing CMO in CALL RSPOLI
920722-hjaaj ORPNEX,ORPST: new THRLDV definition for lin.dep.
920721-Hinne Hettema ORPORT: inserted averaging code (if RSPSUP)
===========================================================================
#endif
C  /* Deck orpctl */
      SUBROUTINE ORPCTL(MAXORP,KOZRED,AOFREQ,THCORP,A1,REDOE,REDOS,
     *                  CMO,UDV,PVX,FC,FV,FCAC,
     *                  REDA1,OBVEC,OEVEC,OSVEC,SOLEQ,ORBDIA,
     *                  XINDX,WRK,LWRK)
C
C Written 6-JUL-1986 by Poul Joergensen
C Revisions:
C
C THIS SUBROUTINE DIRECT THE SOLUTION THE LINEAR
C SET OF ORBITAL RESPONSE EQUATIONS
C
C ( OE(2)-W*OS(2) )X(I) = A1
C
C WHICH OCCUR WHEN AN OPTIMAL SET OF ORBITAL PARAMETERS
C IS DETERMINED FOR A FIXED SET OF CSF COEFFICIENTS.
C A1 IS THE MODIFIED GRADIENT
C
C THE PAIRED STRUCTURE OF THE LINEAR TRANSFORMATIONS
C
C    OE(2)(X1,X2)T=(U,V)T ; OE(2)(X2,X1)T=(V,U)T
C    OS(2)(X1,X2)T=(Q,R)T ; OS(2)(X2,X1)T=(-R,-Q)T
C
C IS USED TO SET UP THE REDUCED LINEAR SET OF ORBITAL RESPONSE EQUATIONS
C
C KORPST:  NUMBER OF START VECTORS
C MAXPIT: MAXIMUM NUMBER OF MICROITERATIONS
C
#include <implicit.h>
      DIMENSION A1(*),REDOE(*),REDOS(*),REDA1(*)
      DIMENSION OBVEC(KZYWOP,*),OEVEC(KZYWOP,*),OSVEC(KZYWOP,*)
      DIMENSION SOLEQ(*),ORBDIA(*),XINDX(*),WRK(*)
      DIMENSION CMO(*),UDV(*),PVX(*),FC(*),FV(*),FCAC(*)
C
C Used from common blocks:
C  /WRKRSP/: AFREQ,??
C
#include <priunit.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <infpri.h>
C
      ITMIC  = 0
C
C
      CALL ORPST(NEWVE,KOZRED,OBVEC,OEVEC,A1,ORBDIA)
C
C     CALL ORPST(NEWVE,KOZRED,OBVEC,OEVEC,A1,ORBDIN)
C
      IF (NEWVE.NE.0) THEN
         NSIM = 1
      ELSE
         ICTL = 2
         NSIM = 0
         GO TO 700
      ENDIF
C
C     --- MICRO ITERATION LOOP STARTS AT 100 ---
C
 100  CONTINUE
C
         ITMIC = ITMIC + 1
         IF (IPRRSP.GT.5)
     *   WRITE(LUPRI,'(/A,I5)')' ** ORPCTL MICROITERATION NUMBER :',
     &                          ITMIC
C
         NSIM = 1
         CALL ORPLIN(NSIM,OBVEC(1,KOZRED+1),CMO,UDV,PVX,FC,FV,FCAC,
     *               XINDX,WRK,LWRK)
         CALL DCOPY(KZYWOP,WRK(1),1,OEVEC(1,KOZRED+1),1)
         CALL DCOPY(KZYWOP,WRK(1+KZYWOP),1,OSVEC(1,KOZRED+1),1)
C
C        CALL ORPLIN(KOZRED,OBVEC,OEVEC,OSVEC,WRK,LWRK1)
C
         KOZRED = KOZRED+1
         IF (ABSYM) THEN
            CALL E2OSYM(KOZRED,OBVEC,OEVEC,OSVEC,
     *           CMO,UDV,PVX,FC,FV,FCAC,XINDX,WRK,LWRK)
         END IF
C
C        CALL ORPRED(3,..) INCREASE DIMENSION OF REDUCED RPA EQUATION
C        AND SOLVE FOR EIGENVALUES AND EIGENVECTORS
C
         ICTL = 3
         IF (IPRRSP.GE.10)WRITE(LUPRI,*)'AOFREQ,THCORP',AOFREQ,THCORP
 700     CALL ORPRED(ICTL,KOZRED,NSIM,AOFREQ,REDOE,REDOS,REDA1,A1,SOLEQ,
     *               OBVEC,OEVEC,OSVEC,WRK)
C
C        CALL ORPRED(ICTL,KOZRED,N,AOFREQ,REDOE,REDOS,REDA1,A1,SOLEQ,
C    *                   OBVEC,OEVEC,OSVEC,WRK)
C
C        CREATE NEW LINEAR INDEPENDENT TRIAL VECTOR IN ORPNEX
C        FROM THE REDUCED SET OF LINEAR RESPONSE EQUATIONS
C
         CALL ORPNEX(AOFREQ,THCORP,KOZRED,JOCONV,A1,SOLEQ,ORBDIA,
     *                  OBVEC,OEVEC,OSVEC)
C
C        CALL ORPNEX(AOFREQ,THCORP,KOZRED,JOCONV,A1,SOLEQ,ORBDIN,
C    *                  OBVEC,OEVEC,OSVEC)
C
         IF (JOCONV.LT.0) THEN
C           (LINEAR DEPENDENT NEW TRIAL VECTOR )
            WRITE (LUPRI,'(/2A)') ' *** ORPCTL-MICROITERATIONS STOPPED',
     *         ' BECAUSE OF LINEAR DEPENDENT NEW TRIAL VECTOR'
         ELSE  IF(JOCONV.GT.0)THEN
C           (CONVERGED)
            IF (IPRRSP .GE. 5) WRITE(LUPRI,'(/A,I3,A)')
     *         ' *** ORPCTL-MICROITERATIONS CONVERGED IN',
     *         ITMIC,'ITERATIONS.'
         ELSE
C           (NOT CONVERGED)
            IF (KOZRED.GE.MAXORP) THEN
               WRITE(LUPRI,'(/A/A,I4,A,I4)')
     *         ' ORPCTL: DIMENSION OF REDUCED ORBITAL SPACE TOO LARGE'
     *         ,' KOZRED =',KOZRED,' MAXIMUM ALLOWED MAXORP =',MAXORP
               CALL QUIT(' ORPCTL: TOO LARGE REDUCED ORBITAL SPACE')
            END IF
            IF (ITMIC.GE.MAXITO) THEN
C              (MAX NO. OF ITERATIONS REACHED)
               WRITE (LUPRI,'(/A,I4,A)')
     *            ' *** ORPCTL-MAXIMUM NUMBER OF MICROITERATIONS,',
     *            ITMIC,', REACHED'
            ELSE
               GO TO 100
C     ^-----------------
            END IF
         END IF
C
C     --- END OF MICROITERATION LOOP STARTING AT 100 ---
C
C CONSTRUCT SOLUTION VECTOR
C
      CALL DZERO(A1,KZYWOP)
      DO 5000 K=1,KOZRED
         CALL DAXPY(KZYWOP,SOLEQ(2*K-1),OBVEC(1,K),1,A1,1)
         CALL DAXPY(KZWOPT,SOLEQ(2*K),OBVEC(1,K),1,A1(1+KZWOPT),1)
         CALL DAXPY(KZWOPT,SOLEQ(2*K),OBVEC(1+KZWOPT,K),1,A1,1)
 5000 CONTINUE
C
C     END OF ORPCTL
C
      RETURN
      END
C  /* Deck orplin */
      SUBROUTINE ORPLIN(NOSIM,ZYOVEC,CMO,UDV,PVX,FC,FV,FCAC,
     *                  XINDX,WRK,LWRK)
C
C PURPOSE:
C
C  CARRY OUT THE ORBITAL PART OF A LINEAR TRANSFORMATION
C  ON AN ORBITAL TRIAL VECTOR (ZYOVEC)
C
C     EVECS(I) = E[2]*N(I)
C     SVECS(I) = S[2]*N(I)
C
C OUTPUT :
C
C     EVECS(I) AND SVECS(I)
C
C  THE FIRST KZYWOP VARIABLES IN WRK IN EVECS
C  THE NEXT  KZYWOP VARIABLES IN WRK IN SVECS
C
#include <implicit.h>
C
      DIMENSION ZYOVEC(*),CMO(*),UDV(NASHDI,*),PVX(*)
      DIMENSION FC(*),FV(*),FCAC(*),XINDX(*),WRK(*)
C
C  INFDIM : NASHDI
C
#include <inforb.h>
#include <infdim.h>
#include <inftap.h>
#include <wrkrsp.h>
C
      PARAMETER ( DUMMY = 1.0D20 )
C
C ALLOCATE WORK SPACE FOR RSPOLI
C
      KEVECS = 1
      KSVECS = KEVECS +  KZYWOP
      KZYMAT = KSVECS +  KZYWOP
      KFCX   = KZYMAT + NOSIM * NORBT * NORBT
      KFVX   = KFCX   + NOSIM * NORBT * NORBT
      KQAX   = KFVX   + NOSIM * NORBT * NORBT
      KQBX   = KQAX   + NOSIM * NORBT * NASHT
      KWRK1  = KQBX   + NOSIM * NORBT * NASHT
      LWRK1  = LWRK   - KWRK1
      IF (LWRK1.LT.0) CALL ERRWRK('ORPLIN',KWRK1-1,LWRK)
C
C INITIALIZE EVECS AND SVECS
C
      CALL DZERO(WRK,2*KZYWOP)
C
C
C UNPACK ORBITAL VECTORS
C
         CALL RSPZYM(NOSIM,ZYOVEC,WRK(KZYMAT))
C        CALL RSPZYM(NSIM,ZYVEC,ZYMAT)
C
C
C REDEFINE KZCONF IN RSPOLI TO GET ONLY THE ORBITAL PART OF THE
C LINEAR TRANSFORMATIONS . KZCONF IS RESET TO ORIGINAL VALUE AFTER
C RSPOLI
C
      KZCSAV  = KZCONF
      KZTSAV  = KZVAR
      KZCONF  = 0
      KZYCON  = 0
      KZVAR   = KZWOPT
      KZYVAR  = 2*KZVAR
C
C CONSTRUCT ORBITAL PART OF LINEAR TRANSFORMED VECTORS
C
      CALL RSPOLI(0,NOSIM,UDV,DUMMY,0,FC,FV,PVX,WRK(KZYMAT),
     *            WRK(KFCX),WRK(KFVX), WRK(KQAX),WRK(KQBX),
     *            DUMMY,DUMMY,DUMMY,WRK(KEVECS),
     *            XINDX,CMO, WRK(KWRK1),LWRK1)
C     CALL RSPOLI(NCSIM,NOSIM,UDV,ZYCVEC,LZYCVEC,FC,FV,PVX,ZYMAT,
C    *            FCX,FVX, QAX,QBX,
C    *            FVTD,QATD,QBTD,EVECS, XINDX,CMO,WRK,LWRK)
C
      CALL DZERO(WRK(KSVECS),KZYWOP)
      CALL RSPSLI(0,NOSIM,DUMMY,ZYOVEC,UDV,WRK(KSVECS),XINDX,
     *            WRK(KWRK1),LWRK1)
      KZCONF = KZCSAV
      KZYCON = 2*KZCONF
      KZVAR  = KZTSAV
      KZYVAR = 2*KZVAR
C
C *** END OF ORPLIN
C
      RETURN
      END
C  /* Deck orpnex */
      SUBROUTINE ORPNEX(AOFREQ,THCORP,KOZRED,JOCONV,A1,SOLEQ,ORBDIN,
     *                  OBVEC,OEVEC,OSVEC)
C
C PURPOSE: 1)CONSTRUCT (OE(2)-W*OS(2))X(I)
C            WHERE X(I) IS THE SOLUTION OF THE REDUCED
C            SET OF LINEAR ORBITAL RESPONSE EQUATIONS
C          2)TEST FOR CONVERGENCE OF X(I), CONVERGENCE CRITERIUM:
C            //(OE(2)-W*OS(2))X(I)// .LE. THCRSP * //X(I)//
C          3)USE GENERALIZED CONJUGATE GRADIENT ALGORITHM TO DETERMINE
C            NEXT GUESS OF TRIAL VECTORS
C
C PJ JUL-1986
C
#include <implicit.h>
      DIMENSION A1(*),SOLEQ(*),ORBDIN(*)
      DIMENSION OBVEC(KZYWOP,*),OEVEC(KZYWOP,*),OSVEC(KZYWOP,*)
#include <thrldp.h>
C
C Used from common blocks:
C  /INFRSP/:
C  /WRKRSP/: KZWOPT,KZYWOP
C
#include <priunit.h>
#include <infrsp.h>
#include <wrkrsp.h>
C
C
C CONSTRUCT  (OE(2))X(I) WHERE X(I) THE SOLUTION OF THE
C REDUCED SET OF LINEAR RESPONSE EQUATIONS
C
      CALL DZERO(OBVEC(1,KOZRED+1),KZYWOP)
      DO 900 K=1,KOZRED
         SOLZY = SOLEQ(2*K-1)
         SOLYZ = SOLEQ(2*K)
         CALL DAXPY(KZYWOP,SOLZY,OEVEC(1,K),1,OBVEC(1,KOZRED+1),1)
         CALL DAXPY(KZWOPT,SOLYZ,OEVEC(1+KZWOPT,K),1,
     *                                          OBVEC(1,KOZRED+1),1)
         CALL DAXPY(KZWOPT,SOLYZ,OEVEC(1,K),1,
     *                                   OBVEC(1+KZWOPT,KOZRED+1),1)
 900  CONTINUE
      IF (IPRRSP.GT.95) THEN
         WRITE (LUPRI,*) 'FREQUENCY',AOFREQ
         WRITE (LUPRI,*) ' OE(2) X(1) VECTOR'
         WRITE (LUPRI,'(A/,(5D15.8))') ' Z COMPONENT',
     *        (OBVEC(I,KOZRED+1),I=1,KZWOPT)
         WRITE (LUPRI,'(A/,(5D15.8))') ' Y COMPONENT',
     *      (OBVEC(I+KZWOPT,KOZRED+1),I=1,KZWOPT)
      ENDIF
C
C     ADD -W*OS(2)X(I)
C
      IF (ABS(AOFREQ) .GT. THRLDP) THEN
         DO 500 K=1,KOZRED
            SOLZY = -AOFREQ*SOLEQ(2*K-1)
            SOLYZ =  AOFREQ*SOLEQ(2*K)
            CALL DAXPY(KZYWOP,SOLZY,OSVEC(1,K),1,OBVEC(1,KOZRED+1),1)
            CALL DAXPY(KZWOPT,SOLYZ,OSVEC(1+KZWOPT,K),1,
     *                                          OBVEC(1,KOZRED+1),1)
            CALL DAXPY(KZWOPT,SOLYZ,OSVEC(1,K),1,
     *                                   OBVEC(1+KZWOPT,KOZRED+1),1)
 500     CONTINUE
      END IF
C
C TEST FOR CONVERGENCE OF X(I)
C
      IF (IPRRSP.GT.90) THEN
         WRITE (LUPRI,'(//A)')
     *' (OE(2) - AOFREQ * OS(2) )X(1) VECTOR; PROPERTY VECTOR;
     *   DIFFERENCE'
         WRITE (LUPRI,'(/A//,(1P,3D20.8))') ' Z COMPONENT',
     *      (OBVEC(I,KOZRED+1),A1(I),OBVEC(I,KOZRED+1)-A1(I),I=1,KZWOPT)
         WRITE (LUPRI,'(/A//,(1P,3D20.8))') ' Y COMPONENT',
     *      (OBVEC(I,KOZRED+1),A1(I),OBVEC(I,KOZRED+1)-A1(I),
     *                                                I=KZWOPT+1,KZYWOP)
      ENDIF
      DO 1050 I=1,KZYWOP
 1050    OBVEC(I,KOZRED+1) = OBVEC(I,KOZRED+1) - A1(I)
      IF (IPRRSP.GT.150) THEN
         WRITE (LUPRI,*) ' RESIDUAL VECTOR; OE(2) X(1) - A1 VECTOR'
         WRITE (LUPRI,'(A/,(5D15.8))') ' Z COMPONENT',
     *      (OBVEC(I,KOZRED+1),I=1,KZWOPT)
         WRITE (LUPRI,'(A/,(5D15.8))') ' Y COMPONENT',
     *      (OBVEC(I+KZWOPT,KOZRED+1),I=1,KZWOPT)
      ENDIF
      QNORM = DDOT(KZYWOP,OBVEC(1,KOZRED+1),1,OBVEC(1,KOZRED+1),1)
      QNORM = SQRT(QNORM)
      IF (QNORM.LE.THCORP) GO TO 2000
C
C     USE A GENERALIZED CONJUGATE GRADIENT ALGORITHM TO
C     FORM NEW TRIAL VECTORS
C
      IF (IPRRSP.GT.180) THEN
         WRITE (LUPRI,*) ' DIAGONAL, AOFREQ=',AOFREQ
         WRITE (LUPRI,'(A/,(5D15.8))') ' Z COMPONENT',
     &        (ORBDIN(I),I=1,KZWOPT)
         WRITE (LUPRI,'(A/,(5D15.8))') ' Y COMPONENT',
     *      (ORBDIN(I+KZWOPT),I=1,KZWOPT)
      ENDIF
      DO 1400 I=1,KZYWOP
         OBVEC(I,KOZRED+1) = OBVEC(I,KOZRED+1) * ORBDIN(I)
 1400 CONTINUE
      IF (IPRRSP.GT.125) THEN
         WRITE(LUPRI,'(/A)') ' BEFORE ORPORT '
         WRITE(LUPRI,'(/A,I10)')
     *    ' OBVEC(1,KOZRED+1) : KOZRED=',(KOZRED+1)
         CALL OUTPUT(OBVEC,1,KZYWOP,1,(KOZRED+1),
     *               KZYWOP,(KOZRED+1),1,LUPRI)
      END IF
C
C ORTHOGONALIZE TRIAL VECTORS AND EXAMINE FOR LINEAR DEPENDENCE
C
      NSIM = 1
      THRLDV = KZYWOP*THRLDP
      CALL ORPORT(OBVEC,NSIM,KOZRED,THRLDV,OEVEC(1,KOZRED+1))
      IF (IPRRSP.GT.125) THEN
         WRITE(LUPRI,'(/A)') ' AFTER ORPORT '
         WRITE(LUPRI,'(/A,I10)')
     *    ' OBVEC(1,KOZRED+1) : KOZRED=',(KOZRED+1)
         CALL OUTPUT(OBVEC,1,KZYWOP,1,(KOZRED+1),
     *               KZYWOP,(KOZRED+1),1,LUPRI)
      END IF
C
      IF (IPRRSP.GE.5)WRITE(LUPRI,5030) QNORM
 5030 FORMAT(/'*** CONVERGENCE OF NEW TRIAL VECTOR TO:',
     *       /' NORM OF RESIDUAL IS: ',1P,D15.8)
C
      IF (NSIM.EQ.0) THEN
C
C        IF NSIM=0 THEN NEW TRIAL VECTOR IS LINEAR
C        INDEPENDENT OF PREVIOUS VECTORS
C
         WRITE(LUPRI,5010)
 5010    FORMAT(/'*** MICROITERATIONS STOPPED DUE TO LINEAR',
     *          ' DEPENDENCE BETWEEN NEW TRIAL VECTORS')
         IF (IPRRSP.GT.5) WRITE (LUPRI,97) QNORM
 97      FORMAT(/' ORPNEX: Linear  equations NOT converged;',
     *       /' Norm of residual is: ',1P,D15.8)
         JOCONV = -1
      ELSE
         IF (IPRRSP.GT.5) WRITE (LUPRI,97) QNORM
         JOCONV=0
      END IF
C
      RETURN
C
 2000 IF (IPRRSP.GT.5) WRITE (LUPRI,96) QNORM
 96   FORMAT(/' ORPNEX: Linear  equations converged;',
     *       /' Norm of residual is: ',1P,D15.8)
      JOCONV = 1
      RETURN
C
C     END OF ORPNEX
C
      END
C  /* Deck orport */
      SUBROUTINE ORPORT (BVECS,NBX,NBPREV,THRLDP,OLDVEC)
C
C Written 25-Oct-1984 by Hans Jorgen Aa. Jensen
C (based on ORTVEC from MSCSF section)
C REVISED: 22-Jan-1986 hjaaj (Z = Y or Z = -Y check).
C REVISED:  6-JUL-1986 pj (assume bvectors in core)
C REVISED: 21-Jul-1992 hh Averaging if RSPSUP .AND. KSYMOP.EQ.1 added
C Purpose:
C    Orthogonalize the new b-vector against all previous b-vectors
C    Each b-vector is in (Z, Y) form, the (Y, Z) vector to be obtained
C  by switching first and second half. Each of the new b-vectors are
C  orthogonalized both against (Z, Y) and (Y, Z) for each previous
C  b-vector.
C  (Orthogonalization is performed twice if round-off is large,
C   if larger than THRRND).
C
C Input:
C  BVECS(1,1...NBPREV) ,  orthogonal b-vectors
C  BVECS(1,NBPREV+1) , the new b-vector to be orthonormalized
C  NBX=1, the new b-vector in BVECS. NBX MUST BE ONE
C  NBPREV, number of previous b-vectors
C  THRLDP, threshold for linear dependence
C
C Output:
C  BVECS, orthogonal b-vectors
C
C Scratch:
C  OLDVEC(KZYWOP), scratch array for old b-vectors on LU3
C
#include <implicit.h>
      DIMENSION BVECS(KZYWOP,*),OLDVEC(*)
C
      PARAMETER (ZEQLY=1.D-10, THRRND=1.D-4)
      PARAMETER (DP5=0.5D0, D1=1.0D0, D2=2.0D0)
C
C Used from common blocks:
C  /INFRSP/: RSPSUP
C  /WRKRSP/: KZWOPT,KZYWOP,JEX(20),LU3
C
#include <priunit.h>
#include <infrsp.h>
#include <wrkrsp.h>
C
C
C     For zero frequency linear response, the result will
C     be Z = Y for real perturbations and Z = -Y for imaginary
C     perturbations; trial vectors will have the same structue.
C
C     However, if Z = Y or Z = -Y then (Z Y) and (Y Z) are linear
C     dependent and we want instead (Z 0) and (0 Z) as trial vectors.
C
      IF (NBX .NE. 1) CALL QUIT('ORPORT error: NBX .ne. 1')
      NBNEW = NBPREV +NBX
      T1 = DDOT(KZYWOP,BVECS(1,NBNEW),1,BVECS(1,NBNEW),1)
C     test if Z + Y is zero
      DO 1100 I = 1,KZWOPT
         OLDVEC(I) = ( BVECS(I,NBNEW) +
     *                              BVECS(KZWOPT+I,NBNEW) )
 1100 CONTINUE
      TZPY = DDOT(KZWOPT,OLDVEC,1,OLDVEC,1)
      IF (IPRRSP .GE. 50) WRITE (LUPRI,'(/A/,(A,1P,D15.8))')
     *   ' ** TEST OF Z = Y OR Z = -Y:',
     *   ' (Z Y) norm squared =',T1,
     *   ' Z + Y norm squared =',TZPY
      IF (TZPY .LE. ZEQLY*T1) THEN
         IF (IPRRSP .GE. 50) WRITE (LUPRI,'(/A)')
     *      ' Z = -Y IN TRIAL VECTOR. Y COMPONENT REMOVED.'
         CALL DZERO(BVECS(1+KZWOPT,NBNEW),KZWOPT)
      ELSE IF (TZPY .GT. T1) THEN
C        test if Z - Y is zero;
C        only if TZPY = 2 * T1 can Z - Y be zero,
C        so TZPY .gt. T1 is a safe test.
         DO 1200 I = 1,KZWOPT
            OLDVEC(I) = ( BVECS(I,NBNEW) -
     *                             BVECS(KZWOPT+I,NBNEW) )
 1200    CONTINUE
         TZMY = DDOT(KZWOPT,OLDVEC,1,OLDVEC,1)
         IF (IPRRSP .GE. 50) WRITE (LUPRI,'(A,1P,D15.8)')
     *      ' Z - Y norm squared =',TZMY
         IF (TZMY .LE. ZEQLY*T1) THEN
            IF (IPRRSP .GE. 50) WRITE (LUPRI,'(/A)')
     *         ' Z = Y IN TRIAL VECTOR. Y COMPONENT REMOVED.'
            CALL DZERO(BVECS(1+KZWOPT,NBNEW),KZWOPT)
         END IF
      END IF
C
C *** Average the new orbital trial vector Z and Y parts, if necessary
C
      IF ( RSPSUP .AND. (KSYMOP .EQ. 1)) THEN
         CALL RSPAVE(BVECS(1,NBNEW),KZWOPT,2)
      END IF
C
      IROUND=0
      ITURN=0
 1500 ITURN=ITURN+1
C
C        Orthogonalize new b-vectors agains previous b-vectors
C        (both (Z, Y) and (Y, Z))
C
         DO 2000 K=1,NBPREV
            TT = -DDOT(KZYWOP,BVECS(1,K),1,BVECS(1,NBNEW),1)
            CALL DAXPY(KZYWOP,TT,BVECS(1,K),1,BVECS(1,NBNEW),1)
            TT = -DDOT(KZWOPT,BVECS(1,K),1,BVECS(1+KZWOPT,NBNEW),1)
     *           -DDOT(KZWOPT,BVECS(1+KZWOPT,K),1,BVECS(1,NBNEW),1)
            CALL DAXPY(KZWOPT,TT,BVECS(1+KZWOPT,K),1,
     *                                     BVECS(1,NBNEW),1)
            CALL DAXPY(KZWOPT,TT,BVECS(1,K),1,
     *                                     BVECS(1+KZWOPT,NBNEW),1)
 2000    CONTINUE
C
C        Normalization
C
         TT = DDOT(KZYWOP,BVECS(1,NBNEW),1,BVECS(1,NBNEW),1)
         IF (TT .LE. THRLDP) THEN
            IF (IPRRSP.GE.5) WRITE (LUPRI,3250) TT
            NBX= 0
            RETURN
         ELSE
            IF (TT .LT. THRRND) IROUND=IROUND+1
            TT = D1 / SQRT(TT)
            CALL DSCAL(KZYWOP,TT,BVECS(1,NBNEW),1)
         ENDIF
         IF (IROUND.GT.0 .AND. ITURN.EQ.1) GO TO 1500
C     ^----------------------------------------------
 3250 FORMAT(/' ORPORT, b-vector is removed because of linear ',
     *        'dependence; ',1P,D12.5)
C
C Perform symmetric orthonormalization of (Z Y) and (Y Z) pair.
C
C   -1/2       ( C1   C2 )               (  1     OVLPI )
C  S      =   (           ) where S  =  (                )
C              ( C2   C1 )               ( OVLPI     1  )
C
C
      X1 = DDOT(KZYWOP,BVECS(1,NBNEW),1,BVECS(1,NBNEW),1)
      X1 = SQRT(X1)
      OVLPI = D2*DDOT(KZWOPT,BVECS(1,NBNEW),1,
     *                          BVECS(1+KZWOPT,NBNEW),1)
      IF (IPRRSP .GE. 50) WRITE (LUPRI,'(/A,1P,2D20.8)')
     *      ' ORPORT, NORM AND <ZY/YZ> OVERLAP BEFORE S-1/2',X1,OVLPI
      IF (ABS(OVLPI) .GT. THRLDP) THEN
         X1 = DP5 / SQRT(D1+OVLPI)
         X2 = DP5 / SQRT(D1-OVLPI)
         C1 = X1 + X2
         C2 = X1 - X2
         CALL DCOPY(KZYWOP,BVECS(1,NBNEW),1,OLDVEC,1)
         CALL DSCAL(KZYWOP,C1,BVECS(1,NBNEW),1)
         CALL DAXPY(KZWOPT,C2,OLDVEC,1,BVECS(1+KZWOPT,NBNEW),1)
         CALL DAXPY(KZWOPT,C2,OLDVEC(1+KZWOPT),1,
     *                                      BVECS(1,NBNEW),1)
         X1 = DDOT(KZYWOP,BVECS(1,NBNEW),1,BVECS(1,NBNEW),1)
         X1 = SQRT(X1)
         OVLPI = D2*DDOT(KZWOPT,BVECS(1,NBNEW),1,
     *                                      BVECS(1+KZWOPT,NBNEW),1)
         IF (IPRRSP .GE. 50) WRITE (LUPRI,'(A,1P,2D20.8)')
     *         ' ORPORT, NORM AND <ZY/YZ> OVERLAP AFTER  S-1/2',X1,OVLPI
      END IF
C
C *** End of subroutine ORPORT
C
      RETURN
      END
C  /* Deck orppar */
      SUBROUTINE ORPPAR(NOSIM,THCORP,EOVAL,IBTYP,A1,ORBDIE,ORBDIS,
     *                  CMO,UDV,PVX,FC,FV,FCAC,XINDX,WRK,LWRK)
C
C PURPOSE:
C  DETERMINE OPTIMAL SET OF ORBITAL PARAMETERS FOR A FIXED
C  SET OF CSF COEFFICIENTS
C
C INPUT:
C  A1, GRADIENT CONTAINING COUPLING ELEMENTS FOR FIXED CSF COEFFICIENTS.
C OUTPUT:
C  A1, OPTIMAL ORBITAL TRIAL VECTOR
C
#include <implicit.h>
C
      DIMENSION THCORP(*),EOVAL(*),IBTYP(*)
      DIMENSION CMO(*),UDV(*),PVX(*),FC(*),FV(*),FCAC(*)
      DIMENSION A1(KZYWOP,*),ORBDIE(*),ORBDIS(*),XINDX(*),WRK(*)
C
      PARAMETER ( DTEST = 1.0D-4 , DE4 = 1.0D4  )
      PARAMETER ( D1 = 1.0D0 )
C
#include <priunit.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <inftap.h>
#include <infdim.h>
#include <inforb.h>
#include <infpri.h>
C
C ALLOCATE WORK SPACE
C
C MAXORP SET TO HALF MAXRM ASSUMING THAT WITH OPTORB YOU WILL
C NEVER GET MORE THAN HALF OF THE TRIAL VECTORS TO BE ORBITAL
C TRIAL VECTORS
C
C930920-pj+hjaaj: revised next two lines
C  such that MAXORP never can be too small.
      MAXORP = KZRED + NOSIM*(MAXITO+2)
      MXORP  = 2 * MAXORP
C
      KREDOE = 1
      KREDOS = KREDOE + MXORP*(MXORP+1)/2
      KREDA1 = KREDOS + MXORP*(MXORP+1)/2
      KOBVEC = KREDA1 + MXORP
      KOEVEC = KOBVEC + KZYWOP*MAXORP
      KOSVEC = KOEVEC + KZYWOP*MAXORP
      KSOLEQ = KOSVEC + KZYWOP*MAXORP
      KDIAHE = KSOLEQ + MAXORP
      KWRK1  = KDIAHE + KZYWOP
      LWRK1  = LWRK   - KWRK1
      IF (LWRK1.LT.0) CALL ERRWRK('ORPPAR',KWRK1-1,LWRK)
      LWRKRE = MAXORP*(MAXORP+1)/2 + MAXORP
      LWRKLI = 2*KZYWOP + 5*NORBT*NORBT + 4*NASHT*NORBT +
     *         NASHT* NASHT + 2*LBINTM
C     LNDSLI = space needed for RSPSLI in ORPSVE
      LNDSLI = KZYVAR + NCREF + N2ORBX + N2ASHX
      LNEED  = MAX(LWRKRE,LWRKLI,LNDSLI)
      IF (LWRK1.LE.LNEED) THEN
         WRITE(LUPRI,*)' NOT ENOUGH SPACE FOR ORPCTL'
         WRITE(LUPRI,*)' LWRK1,LWRKRE,LWRKLI,LNDSLI'
         WRITE(LUPRI,*)  LWRK1,LWRKRE,LWRKLI,LNDSLI
         CALL QUIT('*** ORPPAR, INSUFFICIENT SPACE. ')
      ENDIF
      IF (ABOCHK) THEN
         CALL E2OCHK(WRK(KOBVEC),WRK(KOEVEC),WRK(KOSVEC),
     *              CMO,UDV,PVX,FC,FV,FCAC,XINDX,WRK(KWRK1),LWRK1)
C        CALL E2OCHK(OBVEC,OEVEC,OSVEC,
C    *              CMO,UDV,PV,FC,FV,FCAC,XINDX,WRK,LWRK)
         WRITE(LUPRI,'(/A)')
     *   ' BLOCK OF ORBITAL SET UP EXPLICITLY : ABOCHK COMPLETED'
         CALL QUIT(' ABOCHK COMPLETED')
      END IF
C
C SET UP INVERSE DIAGONAL ORBITAL HESSIAN
C
      IF (IPRRSP.GE.150) THEN
         WRITE(LUPRI,*)' ORBDIE(I) '
         CALL OUTPUT(ORBDIE,1,KZYWOP,1,1,KZYWOP,1,1,LUPRI)
         WRITE(LUPRI,*)' ORBDIS(I)'
         CALL OUTPUT(ORBDIS,1,KZYWOP,1,1,KZYWOP,1,1,LUPRI)
      ENDIF
      IF (IPRRSP.GE.110) THEN
         WRITE(LUPRI,'(/A)')'  ***** ORPPAR ****'
         WRITE(LUPRI,'(/I5,A)')NOSIM,' MODIFIED GRADIENT VECTORS'
         CALL OUTPUT(A1,1,KZYWOP,1,NOSIM,KZYWOP,NOSIM,1,LUPRI)
      END IF
      KOZRED = 0
      NLOAD = 0
      DO 100 ISIM = 1,NOSIM
         NLOAD = NLOAD + 1
         ESHIFT = -EOVAL(ISIM)
         IF(IPRRSP.GE.70)WRITE(LUPRI,*)' ESHIFT',ESHIFT
         DO 200 K = 1,KZYWOP
            DIAHES = ORBDIE(K) + ESHIFT*ORBDIS(K)
               IF (ABS(DIAHES).LE.DTEST) THEN
                  WRK(KDIAHE-1+K) = SIGN(DE4,DIAHES)
               ELSE
                  WRK(KDIAHE-1+K) = D1 / DIAHES
               ENDIF
 200     CONTINUE
         IF ((NLOAD.EQ.1).AND.(KZRED.GT.0)) THEN
C
            CALL ORPSVE(IBTYP,KOZRED,A1(1,ISIM),WRK(KREDA1),
     *                  WRK(KREDOE),WRK(KREDOS),WRK(KOBVEC),WRK(KOEVEC),
     *                  WRK(KOSVEC),UDV,XINDX,WRK(KWRK1),LWRK1)
C
C           CALL ORPSVE(IBTYP,KOZRED,GD,REDGD,REDOE,REDOS,
C    *                 OBVEC,OEVEC,OSVEC,UDV,XINDX,WRK,LWRK)
C
         ELSE IF (KOZRED.NE.0) THEN
            DO 800 I = 1,KOZRED
               PZY = DDOT(KZYWOP,WRK(KOBVEC+(I-1)*KZYWOP),1,
     *               A1(1,ISIM),1)
               PYZ = DDOT(KZWOPT,WRK(KOBVEC+(I-1)*KZYWOP),1,
     *               A1(1+KZWOPT,ISIM),1)
     *             + DDOT(KZWOPT,WRK(KOBVEC+(I-1)*KZYWOP+KZWOPT),1,
     *               A1(1,ISIM),1)
               WRK(KREDA1+(I-1)*2)     = PZY
               WRK(KREDA1+(I-1)*2+1)   = PYZ
 800        CONTINUE
         ENDIF
      IF (IPRRSP.GE.80) WRITE(LUPRI,*)'EOVAL(ISIM),THCORP(ISIM)'
     *                            ,EOVAL(ISIM),THCORP(ISIM)
      IF (IPRRSP.GE.200) THEN
         WRITE(LUPRI,*)' INVERSE ORBITAL HESSIAN'
         CALL OUTPUT(WRK(KDIAHE),1,KZYWOP,1,1,KZYWOP,1,1,LUPRI)
      ENDIF
         CALL ORPCTL(MAXORP,KOZRED,EOVAL(ISIM),THCORP(ISIM),
     *            A1(1,ISIM),WRK(KREDOE),WRK(KREDOS),
     *            CMO,UDV,PVX,FC,FV,FCAC,
     *            WRK(KREDA1),WRK(KOBVEC),WRK(KOEVEC),WRK(KOSVEC),
     *            WRK(KSOLEQ),WRK(KDIAHE),XINDX,WRK(KWRK1),LWRK1)
C        CALL ORPCTL(MAXORP,KOZRED,AOFREQ,THCORP,
C    *               A1,REDOE,REDOS,
C    *               CMO,UDV,PVX,FC,FV,FCAC,
C    *               REDA1,OBVEC,OEVEC,OSVEC,
C    *               SOLEQ,ORBDIA,XINDX,WRK,LWRK)
 100  CONTINUE
      RETURN
      END
C  /* Deck orpred */
      SUBROUTINE ORPRED(ICTL,KOZRED,N,AOFREQ,REDOE,REDOS,REDA1,A1,SOLEQ,
     *                   OBVEC,OEVEC,OSVEC,WRK)
C
C SOLVE LINEAR SET OF ORBITAL RESPONSE EQUATIONS IN REDUCED SPACE
C
C  (OE(2)-W*OS(2))X = A1
C
C THE STRUCTURE OF OE(2) AND OS(2) IS USED TO FORM A REDUCED
C SET OF LINEAR EQUATIONS OF DIMENSION 2*KOZRED ALTHOUGH
C ONLY  KOZRED LINEAR TRANSFORMATIONS ARE CARRIED OUT.
C
C Input:
C  ICTL, flow control
C       =1, extend reduced set of linear equations
C       =2, solve reduced set of linear equations
C       =3, extend and solve
c
C  REDOE,REDOS, the old reduced matrices (dimension: KOZYRE-2*N)
C  REDA1 the old reduced property vector (dimension:KOZYRE-2*N)
C  SOLEQ, solution to the reduced set of linear equations
C  OBVEC, the N new vectors
C  OEVEC, OE times the N OB-vectors
C  OSVEC  OS times the N OB-vectors
C  N, number of new vectors
C
C
C Output:
C  REDOE, the new, extended reduced OE[2]-matrix (dimension: KOZYRE)
C  REDOS, the new, extended reduced OS[2]-matrix (dimension: KOZYRE)
C  REDA1, the new,extended  reduced A1-vector (dimension: KOZYRE)
C  SOLEQ, the solution to the REDUCed set of linear equations
C
C Scratch:
C  WRK, real scratch array dimension KZYWOP
#include <implicit.h>
      DIMENSION REDOE(*),REDOS(*),REDA1(*),A1(*),SOLEQ(*),
     *          OBVEC(KZYWOP,*),OEVEC(KZYWOP,*),OSVEC(KZYWOP,*), WRK(*)
C
C Used from common blocks:
C  /WRKRSP/: KZWOPT,KOZYRE,KOZRED,KZYWOP,AFREQ,?
C  /INFIND/: IROW(*). ?
C  /INFPRI/: ??
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <infind.h>
#include <infpri.h>
C
      KOZYRE = 2*KOZRED
CSPAS:01.05.2009 bug fix for large number of  excitation energies
      IF (KOZYRE.GT.LIROW) THEN
         WRITE(LUPRI,'(//A/A,A,I5,A,I5,/A)')
     *   ' >>> ERROR IN ORPRED >>>',
     *   ' DIMENSION OF INDEX VECTOR IROW IS TO SMALL; ',
     *   ' NEEDED ',KOZYRE,', GIVEN :',LIROW,
     *   ' CHANGE DIMENSION IN infind.h AND RECOMPILE'
         CALL QUIT('ORPRED: TOO SMALL DIMENSION OF IROW: RECOMPILE')
      END IF
CKeinSPASmehr
C
C Section 1: extend reduced OE(2)-AND OS(2)-matrices
C            and A1 vector
C            with N new b-vectors
C
      IF (ICTL.LT.1 .OR. ICTL.GT.3) GO TO 2000
      IF (ICTL.EQ.2) GO TO 1000
C
C     b-vector times property vector
C
      IOZRED = KOZRED - N
      DO 10 I=1,N
         PZY = DDOT(KZYWOP,OBVEC(1,I+IOZRED),1,A1,1)
         PYZ = DDOT(KZWOPT,OBVEC(1,I+IOZRED),1,A1(1+KZWOPT),1)
     *       + DDOT(KZWOPT,OBVEC(1+KZWOPT,I+IOZRED),1,A1,1)
         REDA1((IOZRED+I)*2-1) = PZY
         REDA1((IOZRED+I)*2)   = PYZ
 10   CONTINUE
C
C     ob-vectors times new S-and E-vectors
C
C
C
      DO 15 J = 1,KOZRED
         IF(J.LE.IOZRED) THEN
            KMIN=1
         ELSE
            KMIN=J-IOZRED
         END IF
         DO 13 K=KMIN,N
            X1 = DDOT(KZYWOP,OBVEC(1,J),1,OEVEC(1,IOZRED+K),1)
            X2 = DDOT(KZWOPT,OBVEC(KZWOPT+1,J),1,OEVEC(1,IOZRED+K),1)
     +         + DDOT(KZWOPT,OBVEC(1,J),1,OEVEC(KZWOPT+1,IOZRED+K),1)
            REDOE(2*J   + IROW(2*(IOZRED+K)))   = X1
            REDOE(2*J-1 + IROW(2*(IOZRED+K)-1)) = X1
            IF(J.NE.(K+IOZRED)) REDOE(2*J+IROW(2*(IOZRED+K)-1))=X2
            REDOE(2*J-1+IROW(2*(IOZRED+K)))=X2
            X1 = DDOT(KZYWOP,OBVEC(1,J),1,OSVEC(1,IOZRED+K),1)
            X2 = DDOT(KZWOPT,OBVEC(KZWOPT+1,J),1,OSVEC(1,IOZRED+K),1)
     +         + DDOT(KZWOPT,OBVEC(1,J),1,OSVEC(KZWOPT+1,IOZRED+K),1)
            REDOS(2*J-1+IROW(2*(IOZRED+K)-1)) =  X1
            REDOS(2*J+IROW(2*(IOZRED+K)))     = -X1
            IF(J.NE.(K+IOZRED)) REDOS(2*J+IROW(2*(IOZRED+K)-1)) = X2
            REDOS(2*J-1+IROW(2*(IOZRED+K)))   = -X2
   13    CONTINUE
   15 CONTINUE
C
C *********************************************************
C Section 2: find solution (SOLEQ) to the reduced set of
C            linear response equations
C
 1000 CONTINUE
      IF (ICTL .EQ. 1) GO TO 2000
C
C     Solve reduced linear orbital response problem in subspace
C
      IJ=0
      DO 70 I=1,KOZYRE
         DO 70 J=1,I
             IJ = IJ + 1
             WRK(IJ) = REDOE(IJ) - AOFREQ*REDOS(IJ)
   70 CONTINUE
      CALL DCOPY(KOZYRE,REDA1,1,SOLEQ,1)
C
C     USE DSPSOL ( WHICH CALLS LINPACK ROUTINES ) TO SOLVE
C     SYSTEM OF LINEAR EQUATIONS
C
      NSIM = 1
      KWRKDS = 1 + KOZYRE*(KOZYRE+1)/2
      CALL DSPSOL(KOZYRE,NSIM,WRK,SOLEQ,WRK(KWRKDS),INFO)
C
C     CALL DSPSOL(N,NSIM,AP,B,KPIVOT,INFO)
C
      IF (INFO .NE. 0) THEN
         WRITE (LUPRI,'(/A,I5)')
     &      ' *** ERROR IN ORPRED.DSPSOL, INFO =',INFO
         CALL QUIT('*** ERROR IN ORPREDO.DSPSOL')
      END IF
C
C
      IF (IPRRSP .GE. 46) THEN
         WRITE (LUPRI,'(/A)')
     *      ' REDUCED SOLUTION:'
         WRITE (LUPRI,'(/I5,1P,D15.6/I5,D15.6)')
     *      (K, SOLEQ(K), K=1,KOZYRE)
      END IF
C
 2000 CONTINUE
C
      IF (IPRRSP .GE. 46) THEN
         WRITE (LUPRI,'(/A)')
     *      ' REDUCED PROPERTY VECTOR:'
         WRITE (LUPRI,'(/I5,1P,D15.6/I5,D15.6)')
     *      (K,REDA1(K), K=1,KOZYRE)
      END IF
      IF (IPRRSP .GE. 46) THEN
         WRITE (LUPRI,'(/A)') ' REDUCED HESSIAN MATRIX:'
         CALL OUTPAK(REDOE,KOZYRE,1,LUPRI)
         WRITE (LUPRI,'(/A)') ' REDUCED METRIC MATRIX:'
         CALL OUTPAK(REDOS,KOZYRE,1,LUPRI)
      END IF
C
C
C *** End of subroutine ORPREDO
C
      RETURN
      END
C  /* Deck orpst */
      SUBROUTINE ORPST(NEWVE,KOZRED,OBVEC,OEVEC,A1,ORBDIN)
C
C PURPOSE: CREATE START VECTOR(S) FOR SOLUTION OF A LINEAR
C          SET OF EQUATIONS
C
C          USE GRADIENT VECTOR MULTIPLIED WITH INVERSE
C          DIAGONAL HESSIAN MATRIX ELEMENTS
C
C
#include <implicit.h>
      DIMENSION OBVEC(KZYWOP,*),OEVEC(KZYWOP,*),A1(*),ORBDIN(*)
#include <thrldp.h>
C
#include <priunit.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <infpri.h>
C
      DO 100 I=1,KZYWOP
            OBVEC(I,KOZRED+1) = A1(I)*ORBDIN(I)
 100  CONTINUE
C
C     NORMALIZE AND SAVE TRIAL VECTOR.
C
      IF (IPRRSP.GT.135) THEN
         WRITE(LUPRI,'(/A)')' ORPST: NON-ORTHOGONALIZED TRIAL VECTOR'
         CALL OUTPUT(OBVEC(1,KOZRED+1),1,KZYWOP,1,1,KZYWOP,1,1,LUPRI)
      END IF
      NEWVE = 1
      THRLDV = KZYWOP*THRLDP
      CALL ORPORT(OBVEC,NEWVE,KOZRED,THRLDV,OEVEC(1,KOZRED+1))
C
C     CALL ORPORT(OBVEC,NBX,NBPREV,THRLDP,OLDVEC)
C
      IF (IPRRSP.GT.135) THEN
         WRITE(LUPRI,'(/A)')' ORPST: ORTHOGONALIZED TRIAL VECTOR'
         CALL OUTPUT(OBVEC(1,KOZRED+1),1,KZYWOP,1,1,KZYWOP,1,1,LUPRI)
      END IF
      IF (NEWVE .EQ. 0) THEN
         IF (IPRRSP.GE.5) WRITE (LUPRI,'(//A,I5)')
     *     ' ORPST, START VECTOR IS LINEAR INDEPENDENT.'
      END IF
C
C     END OF ORPST.
C
      RETURN
      END
C  /* Deck orpsve */
      SUBROUTINE ORPSVE(IBTYP,KOZRED,GD,REDGD,REDOE,REDOS,
     *                 OBVEC,OEVEC,OSVEC,UDV,XINDX,WRK,LWRK)
C
C PURPOSE: CREATE START VECTOR(S) FROM ORBITAL PART OF TRIAL VECTORS
C          FOR SOLUTION OF AN ORBITAL LINEAR SET OF RESPONSE EQUATIONS
C
#include <implicit.h>
#include <dummy.h>
      DIMENSION IBTYP(*),GD(*),REDGD(*),REDOE(*),REDOS(*)
      DIMENSION OBVEC(KZYWOP,*),OEVEC(KZYWOP,*),OSVEC(KZYWOP,*)
      DIMENSION UDV(*),XINDX(*),WRK(*)
#include <ibndxdef.h>
C
#include <inforb.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <inftap.h>
#include <infpri.h>
C
C WORK space check (ov-920923)
C revised 930920-pj+hjaaj, LNDSLI = work space needed for RSPSLI
C
      LNDSLI = KZYVAR + NCREF + N2ORBX + N2ASHX
      IF (LWRK.LT.LNDSLI) CALL ERRWRK('ORPSVE',LNDSLI,LWRK)
C
      KOZRED = 0
C
C READ IN OLD TRIAL VECTORS AND SET UP REDUCED MATRICES
C
      CALL REWSPL(LURSP3)
      CALL REWSPL(LURSP5)
      JRSP3 = 0
      JRSP5 = 0
      IF (KOFFTY.GT.0) THEN
#if !defined (VAR_SPLITFILES)
         READ (LURSP3)
         READ (LURSP5)
#else
         CALL READST(LURSP3,LRSP3,DUMMY,JRSP3,.FALSE.)
         CALL READST(LURSP5,LRSP5,DUMMY,JRSP5,.FALSE.)
#endif
      END IF
      DO 200 ISIM = 1,KZRED
         IF (IBTYP(KOFFTY+ISIM).EQ.JBCNDX) THEN
#if defined (VAR_SPLITFILES)
            CALL READST(LURSP3,LRSP3,DUMMY,JRSP3,.FALSE.)
            CALL READST(LURSP5,LRSP5,DUMMY,JRSP5,.FALSE.)
#else
            READ (LURSP3)
            READ (LURSP5)
#endif
         ELSE
            KOZRED = KOZRED + 1
#if defined (VAR_SPLITFILES)
            CALL READST(LURSP3,LRSP3,OBVEC(1,KOZRED),JRSP3,.TRUE.)
            IF (LRSP3 .NE. KZYWOP) CALL QUIT('Inconsistent record '//
     &           'length on unit LURSP3 in ORPSVE')
            CALL READST(LURSP5,LRSP5,WRK,JRSP5,.TRUE.)
            IF (LRSP5 .NE. KZYVAR) CALL QUIT('Inconsistent record '//
     &           'length on unit LURSP5 in ORPSVE')
#else
            CALL READT(LURSP3,KZYWOP,OBVEC(1,KOZRED))
            CALL READT(LURSP5,KZYVAR,WRK)
#endif
            CALL DCOPY(KZWOPT,WRK(KZCONF+1),1,OEVEC(1,KOZRED),1)
            CALL DCOPY(KZWOPT,WRK(KZVAR+KZCONF+1),1,
     *                                 OEVEC(1+KZWOPT,KOZRED),1)
            CALL DZERO(WRK(1),KZYVAR)
            CALL RSPSLI(0,1,OBVEC(1,KOZRED),OBVEC(1,KOZRED),UDV,
     *                  WRK(1),XINDX,WRK(1+KZYVAR),(LWRK-KZYVAR-1))
            CALL DCOPY(KZWOPT,WRK(KZCONF+1),1,OSVEC(1,KOZRED),1)
            CALL DCOPY(KZWOPT,WRK(KZVAR+KZCONF+1),1,
     *                                 OSVEC(1+KZWOPT,KOZRED),1)
         ENDIF
 200  CONTINUE
C
      CALL ORPRED(1,KOZRED,KOZRED,DUMMY,REDOE,REDOS,REDGD,GD,DUMMY,
     *                   OBVEC,OEVEC,OSVEC,WRK)
C
C     CALL ORPRED(ICTL,KOZRED,N,AOFREQ,REDOE,REDOS,REDA1,A1,SOLEQ,
C    *                   OBVEC,OEVEC,OSVEC,WRK)
C
C
C     END OF ORPSVE.
C
      RETURN
      END
