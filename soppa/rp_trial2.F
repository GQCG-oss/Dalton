C
C  /* Deck rp_trial2 */
      SUBROUTINE RP_TRIAL2(LTYPE,NONEWT,NOLDTR,NNEWTR,NLINDP,
     &                     EIVAL,LEIVAL,RESINM,LRESINM,CONV,LCONV,
     &                     NCONV,ISYMTR,NEXCI,ENORM,LABEL,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, May 1996
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C     Andrea Ligabue, December 2003: linear response function 
C                                    implemented
C
C     PURPOSE: Determine the residues from the current optimal solution
C              vectors and decide if convergence has been obtained for
C              any of the vectors. For the non-converged vectors create
C              new trial-vectors. These are orthonormalized against the
C              privious optimal trial-vectors and among themself 
C              including the ones obtained from pairing.
C
#include <implicit.h>
#include <priunit.h>
C
#include <soppinf.h>
#include <ccsdsym.h>
#include <cbiexc.h>
C
      PARAMETER   (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER   (D100 = 100.0D0)
C
      LOGICAL     NONEWT
C
      CHARACTER*3 YES, NO
      CHARACTER*3 CONV(LCONV)
      CHARACTER*6 LTYPE
C
      DIMENSION   EIVAL(LEIVAL), RESINM(LRESINM)
      DIMENSION   WORK(LWORK)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_TRIAL2')
C
C-----------------------------
C     Set convergence answers.
C-----------------------------
C
      YES = 'yes'
      NO  = ' no'
C
C--------------------------------------------------------------------
C     Work space allocation no. 1.  Notice that the RESULT and TRIAL
C     vectors are of equal length and are chosen to use the same work
C     space as they are not used at the same time. The same apply 
C     to the overlap-RESULT and the DIAGONAL vectors which also
C     use the same work space.
C--------------------------------------------------------------------
C
      LRES1E  = NT1AM(ISYMTR)
      LTR1E   = LRES1E
      LRES1D  = NT1AM(ISYMTR)
      LTR1D   = LRES1D
      LRESI1E = NT1AM(ISYMTR)
      LRESI1D = NT1AM(ISYMTR)
      LRESO1E = NT1AM(ISYMTR)
      LEDIA1  = LRESO1E
      LRESO1D = NT1AM(ISYMTR)
C
      IF (LTYPE .EQ. 'LINEAR') THEN
         LGPVCH = NT1AM(ISYMTR)
         LGPVEC = 2 * LGPVCH 
      ELSE
         LGPVCH = 1
         LGPVEC = 1
      END IF
C
      KRES1E  = 1
      KTR1E   = KRES1E
      KRES1D  = KRES1E + LRES1E
      KTR1D   = KRES1D
      KEND1A  = KRES1D + LRES1D
      LWORK1A = LWORK  - KEND1A
C
      KRESI1E  = KEND1A
      KRESI1D  = KRESI1E + LRESI1E
      KRESO1E  = KRESI1D + LRESI1D
      KEDIA1   = KRESO1E
      KRESO1D  = KRESO1E + LRESO1E
      KEND1B   = KRESO1D + LRESO1D
C
      IF (LTYPE .EQ. 'LINEAR') THEN
         KGPVEC = KEND1B
         KGPV1E = KGPVEC
         KGPV1D = KGPV1E + LGPVCH
         KEND1C = KGPV1D + LGPVCH
      ELSE
         KGPVEC = KEND1B
         KGPV1E = KEND1B
         KGPV1D = KEND1B
         KEND1C = KEND1B
      END IF
C
      LWORK1C = LWORK  - KEND1C
C
      CALL SO_MEMMAX ('RP_TRIAL2.1',LWORK1C)
      IF (LWORK1C .LT. 0) CALL STOPIT('RP_TRIAL2.1',' ',KEND1C,LWORK)
C
      CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
      CALL SO_OPEN(LUTR1D,FNTR1D,LTR1D)
      CALL SO_OPEN(LURS1E,FNRS1E,LRES1E)
      CALL SO_OPEN(LURS1D,FNRS1D,LRES1D)
      CALL SO_OPEN(LURO1E,FNRO1E,LRESO1E)
      CALL SO_OPEN(LURO1D,FNRO1D,LRESO1D)
C
      IF (LTYPE .EQ. 'LINEAR') CALL SO_OPEN(LUGPVE,FNGPVE,LGPVEC)
C
C-------------------------------------------------------------------
C     Initialize number of new trial vectors and number of converged
C     eigenvalues to zero.
C-------------------------------------------------------------------
C
      NNEWTR = 0
C
      NCONV  = 0
C
C------------------------------------------------
C     Loop over number of excitations considered.
C------------------------------------------------
C
      DO 100 IEXCI = 1,NEXCI
C
C----------------------------
C        Read result vectors.
C----------------------------
C
         CALL SO_READ(WORK(KRES1E), LRES1E, LURS1E,FNRS1E,IEXCI)
         CALL SO_READ(WORK(KRES1D), LRES1D, LURS1D,FNRS1D,IEXCI)
         CALL SO_READ(WORK(KRESO1E),LRESO1E,LURO1E,FNRO1E,IEXCI)
         CALL SO_READ(WORK(KRESO1D),LRESO1D,LURO1D,FNRO1D,IEXCI)
C
         IF (LTYPE .EQ. 'LINEAR') 
     &      CALL SO_READ(WORK(KGPVEC),LGPVEC,LUGPVE,FNGPVE,1)
C
C-------------------------------------------------------------
C        Calculate residual vector and norm of result vectors.
C-------------------------------------------------------------
C
cLig  <> it's good to pass EIVAL(IEXCI) becouse is good both for
cLig     linear response (since we have passed only th right value
cLig     of the freq in EIVAL(1) ) and eigenvalue problem
         CALL RP_RESIDUAL(LTYPE,RESINM(IEXCI),WORK(KRESI1E),LRESI1E,
     &                    WORK(KRESI1D),LRESI1D,EIVAL(IEXCI),
     &                    WORK(KRES1E),LRES1E,WORK(KRES1D),LRES1D,
     &                    WORK(KRESO1E),LRESO1E,WORK(KRESO1D),LRESO1D,
     &                    WORK(KGPV1E),LGPVCH,WORK(KGPV1D),LGPVCH,ENORM)
C
C---------------------------------------------------
C        Determine the threshold for the excitation.
C---------------------------------------------------
C
         IF (IEXCI .LE. (NEXCI - NEXCI2(ISYMTR)) ) THEN
            THR = THREXC
         ELSE
            THR = THREX2
         END IF
C
C---------------------
C        If converged.
C---------------------
C
         IF (RESINM(IEXCI) .LE. THR) THEN
C
                 CONV(IEXCI) = YES
C
            NCONV       = NCONV + 1
C
            GO TO 100
C
C--------------------------------------------
C        If not converged but last iteration.
C--------------------------------------------
C
         ELSE IF (NONEWT) THEN
C
                 CONV(IEXCI) = NO
C
            GO TO 100 
C
C------------------------------------------------
C        If not converged and not last iteration.
C------------------------------------------------
C
         ELSE
C
                 CONV(IEXCI) = NO
C
C--------------------------------------------------------------------
C           Check if space is large enough to hold new trial vectors.
C--------------------------------------------------------------------
C
            IF ( (NNEWTR + NOLDTR) .EQ. (LTR1E) ) GO TO 100
C
C-----------------------------------------------
C           Increase number of new trialvectors.
C-----------------------------------------------
C
            NNEWTR = NNEWTR + 1
C
C------------------------------------------------
C           Read diagonal E[2] and S[2] elements.
C------------------------------------------------
C
            CALL GPOPEN(LUDIAG,'RP_DIAG','UNKNOWN',' ','UNFORMATTED',
     &                  IDUMMY,.FALSE.)
            REWIND LUDIAG
C
            READ(LUDIAG) ( WORK(KEDIA1+I-1), I = 1,LEDIA1)
C
            CALL GPCLOSE(LUDIAG,'KEEP')
C
C------------------------------------------------------------
C           Calculate raw new trial vector and write to file.
C------------------------------------------------------------
C
            CALL RP_NEWTRIAL(NNEWTR,NOLDTR,WORK(KTR1E),LTR1E,
     &                       WORK(KTR1D),LTR1D,EIVAL(IEXCI),
     &                       WORK(KEDIA1),LEDIA1,WORK(KRESI1E),LRESI1E,
     &                       WORK(KRESI1D),LRESI1D)
C
         END IF
C
C-----------------------------------------------------------
C     Write residual on file for later use in calculation of
C     linear response functions
C-----------------------------------------------------------
C
        IF (LTYPE .EQ. 'LINEAR') THEN
C
           CALL GPOPEN(LURV1E,FNRV1E,'UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &                 .FALSE.)
           CALL GPOPEN(LURV1D,FNRV1D,'UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &                 .FALSE.)
C          find the right postion
           CALL SO_FNDPOS(LABEL,EIVAL(IEXCI),LURV1E)
           CALL SO_FNDPOS(LABEL,EIVAL(IEXCI),LURV1D)
C      
           CALL SO_WRTVE(WORK(KRESI1E),LRESI1E,ISYMTR,LABEL,
     &                   EIVAL(IEXCI),LURV1E)
           CALL SO_WRTVE(WORK(KRESI1D),LRESI1D,ISYMTR,LABEL,
     &               EIVAL(IEXCI),LURV1D)
C
           CALL GPCLOSE(LURV1E,'KEEP')
           CALL GPCLOSE(LURV1D,'KEEP')
C
        ENDIF
C
  100 CONTINUE
C
C----------------------------------------------------------------------
C     Orthogonalize new trial vector against all previous trial vectors
C     (including the paired ones) and normalize. Make a symmetric
C     orthonormalization of the the new trial vector and its pair trial
C     vector. Finally write the new trial vector to file. 
C----------------------------------------------------------------------
C
      DTIME      = SECOND()
      CALL RP_ORTH_TRN(LTYPE,NOLDTR,NNEWTR,NLINDP,ISYMTR,WORK,LWORK)
      DTIME      = SECOND()   - DTIME
      SOTIME(39) = SOTIME(39) + DTIME
C
C-------------------------------------------------------------------
C     Decrease the number of new trial vectors with the number which
C     have been removed because of linear dependency.
C-------------------------------------------------------------------
C
      NNEWTR = NNEWTR - NLINDP
C
C-----------------
C     Close files.
C-----------------
C
      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
      CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
      CALL SO_CLOSE(LURS1E,FNRS1E,'KEEP')
      CALL SO_CLOSE(LURS1D,FNRS1D,'KEEP')
      CALL SO_CLOSE(LURO1E,FNRO1E,'KEEP')
      CALL SO_CLOSE(LURO1D,FNRO1D,'KEEP')
C
      IF (LTYPE .EQ. 'LINEAR') CALL SO_CLOSE(LUGPVE,FNGPVE,'KEEP')
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('RP_TRIAL2')
C
      RETURN
C
 9001 FORMAT(/,1X,I3,'. excitation, norm of residual is: ',1P,D15.8,/)
C
      END
