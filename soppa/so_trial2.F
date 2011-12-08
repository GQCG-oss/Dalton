#include <single.h>
C
C  /* Deck so_trial2 */
      SUBROUTINE SO_TRIAL2(LTYPE,NONEWT,NOLDTR,NNEWTR,NLINDP,
     &                     EIVAL,LEIVAL,RESINM,LRESINM,CONV,LCONV,
     &                     NCONV,ISYMTR,NEXCI,DENSIJ,LDENSIJ,
     &                     DENSAB,LDENSAB,ENORM,LABEL,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, May 1996
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C     Andrea Ligabue, January 2004: linear response functions included
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
      DIMENSION   DENSIJ(LDENSIJ), DENSAB(LDENSAB)
      DIMENSION   WORK(LWORK)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_TRIAL2')
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
      LRES2E  = NT2AM(ISYMTR)
      LTR2E   = LRES2E
      LRES2D  = NT2AM(ISYMTR)
      LTR2D   = LRES2D
      LRESI1E = NT1AM(ISYMTR)
      LRESI1D = NT1AM(ISYMTR)
      LRESI2E = NT2AM(ISYMTR)
      LRESI2D = NT2AM(ISYMTR)
      LRESO1E = NT1AM(ISYMTR)
      LEDIA1  = LRESO1E
      LRESO1D = NT1AM(ISYMTR)
      LSDIA1  = LRESO1D
      LRESO2E = NT2AM(ISYMTR)
      LEDIA2  = LRESO2E
      LRESO2D = NT2AM(ISYMTR)
C
      IF(LTYPE.EQ.'LINEAR') THEN
              LGPVC1H = NT1AM(ISYMTR)
              LGPVC2H = NT2AM(ISYMTR)
              LGPVC1  = 2 * LGPVC1H
              LGPVC2  = 2 * LGPVC2H
      ELSE
              LGPVC1H = 1
              LGPVC2H = 1
              LGPVC1  = 1
              LGPVC2  = 1
      ENDIF
C
      KRES1E  = 1
      KTR1E   = KRES1E
      KRES1D  = KRES1E + LRES1E
      KTR1D   = KRES1D
      KRES2E  = KRES1D + LRES1D
      KTR2E   = KRES2E
      KRES2D  = KRES2E + LRES2E
      KTR2D   = KRES2D
      KEND1A  = KRES2D + LRES2D
      LWORK1A = LWORK  - KEND1A
C
      KRESI1E  = KEND1A
      KRESI1D  = KRESI1E + LRESI1E
      KRESI2E  = KRESI1D + LRESI1D
      KRESI2D  = KRESI2E + LRESI2E
      KRESO1E  = KRESI2D + LRESI2D
      KEDIA1   = KRESO1E
      KRESO1D  = KRESO1E + LRESO1E
      KSDIA1   = KRESO1D
      KRESO2E  = KRESO1D + LRESO1D
      KEDIA2   = KRESO2E
      KRESO2D  = KRESO2E + LRESO2E
      KEND1B   = KRESO2D + LRESO2D
C
      IF(LTYPE.EQ.'LINEAR') THEN
              KGPVC1  = KEND1B
              KGPVC1E = KGPVC1
              KGPVC1D = KGPVC1E + LGPVC1H
              KGPVC2  = KGPVC1D + LGPVC1H
              KGPVC2E = KGPVC2  
              KGPVC2D = KGPVC2E + LGPVC2H
              KEND1C  = KGPVC2E + LGPVC2H
      ELSE
              KGPVC1  = KEND1B
              KGPVC1E = KEND1B
              KGPVC1D = KEND1B
              KGPVC2  = KEND1B
              KGPVC2E = KEND1B
              KGPVC2D = KEND1B
              KEND1C  = KEND1B
      ENDIF
C
      LWORK1C  = LWORK   - KEND1C
C
      CALL SO_MEMMAX ('SO_TRIAL2.1',LWORK1C)
      IF (LWORK1C .LT. 0) CALL STOPIT('SO_TRIAL2.1',' ',KEND1C,LWORK)
C
      CALL SO_OPEN(LURS1E,FNRS1E,LRES1E)
      CALL SO_OPEN(LURS1D,FNRS1D,LRES1D)
      CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
      CALL SO_OPEN(LUTR1D,FNTR1D,LTR1D)
      CALL SO_OPEN(LURO1E,FNRO1E,LRESO1E)
      CALL SO_OPEN(LURO1D,FNRO1D,LRESO1D)
      CALL SO_OPEN(LURS2E,FNRS2E,LRES2E)
      CALL SO_OPEN(LURS2D,FNRS2D,LRES2D)
      CALL SO_OPEN(LUTR2E,FNTR2E,LTR2E)
      CALL SO_OPEN(LUTR2D,FNTR2D,LTR2D)
C
      IF(LTYPE.EQ.'LINEAR') THEN
              CALL SO_OPEN(LUGPV1,FNGPV1,LGPVC1)
              CALL SO_OPEN(LUGPV2,FNGPV2,LGPVC2)
      ENDIF
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
         CALL SO_READ(WORK(KRES2E), LRES2E, LURS2E,FNRS2E,IEXCI)
         CALL SO_READ(WORK(KRES2D), LRES2D, LURS2D,FNRS2D,IEXCI)
         CALL SO_READ(WORK(KRESO2E),LRESO2E,LUTR2E,FNTR2E,IEXCI)
         CALL SO_READ(WORK(KRESO2D),LRESO2D,LUTR2D,FNTR2D,IEXCI)
         CALL DSCAL(LRESO2D,-ONE,WORK(KRESO2D),1)
C
         IF(LTYPE.EQ.'LINEAR') THEN
                 CALL SO_READ(WORK(KGPVC1),LGPVC1,LUGPV1,FNGPV1,1)
                 CALL SO_READ(WORK(KGPVC2),LGPVC2,LUGPV2,FNGPV2,1)
                 CALL SO_TMLTR(WORK(KGPVC2E),HALF,ISYMTR)
                 CALL SO_TMLTR(WORK(KGPVC2D),HALF,ISYMTR)
C
           IF(IPRSOP.GT.10) THEN
C
              CALL AROUND("Left GP vector before SO_RESIDUAL")
              CALL OUTPUT(WORK(KGPVC1),1,LGPVC1H,1,2,LGPVC1H,2,1,LUPRI)
              CALL OUTPUT(WORK(KGPVC2),1,LGPVC2H,1,2,LGPVC2H,2,1,LUPRI)
C
           ENDIF
C
         ENDIF
C
C-------------------------------------------------------------
C        Calculate residual vector and norm of result vectors.
C-------------------------------------------------------------
C
         CALL SO_RESIDUAL(LTYPE,RESINM(IEXCI),WORK(KRESI1E),
     &                    LRESI1E,WORK(KRESI1D),LRESI1D,
     &                    WORK(KRESI2E),
     &                    LRESI2E,WORK(KRESI2D),LRESI2D,EIVAL(IEXCI),
     &                    WORK(KRES1E),LRES1E,WORK(KRES1D),LRES1D,
     &                    WORK(KRES2E),LRES2E,WORK(KRES2D),LRES2D,
     &                    WORK(KRESO1E),LRESO1E,WORK(KRESO1D),
     &                    LRESO1D,WORK(KRESO2E),LRESO2E,
     &                    WORK(KRESO2D),LRESO2D,
     &                    WORK(KGPVC1E),LGPVC1H,WORK(KGPVC1D),LGPVC1H,
     &                    WORK(KGPVC2E),LGPVC2H,WORK(KGPVC2D),LGPVC2H,
     &                    ENORM)
C
C-----------------------------------------------
C        Determine the threshold the excitation.
C-----------------------------------------------
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
            IF ( (NNEWTR + NOLDTR) .EQ. (LTR1E + LTR2E) ) GO TO 100
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
            CALL GPOPEN(LUDIAG,'SO_DIAG','UNKNOWN',' ','UNFORMATTED',
     &                  IDUMMY,.FALSE.)
            REWIND LUDIAG
C
            READ(LUDIAG) ( WORK(KEDIA1+I-1), I = 1,LEDIA1)
            READ(LUDIAG) ( WORK(KEDIA2+I-1), I = 1,LEDIA2)
            READ(LUDIAG) ( WORK(KSDIA1+I-1), I = 1,LSDIA1)
C
            CALL GPCLOSE (LUDIAG,'KEEP')
C
C------------------------------------------------------------
C           Calculate raw new trial vector and write to file.
C------------------------------------------------------------
C
            CALL SO_NEWTRIAL(NNEWTR,NOLDTR,WORK(KTR1E),LTR1E,
     &                       WORK(KTR1D),LTR1D,WORK(KTR2E),LTR2E,
     &                       WORK(KTR2D),LTR2D,EIVAL(IEXCI),
     &                       WORK(KEDIA1),LEDIA1,WORK(KEDIA2),LEDIA2,
     &                       WORK(KSDIA1),LSDIA1,WORK(KRESI1E),
     &                       LRESI1E,WORK(KRESI1D),LRESI1D,
     &                       WORK(KRESI2E),LRESI2E,WORK(KRESI2D),
     &                       LRESI2D)
C
         END IF
C
         IF(LTYPE.EQ.'LINEAR') THEN
C
C--------------------------------------------------------------
C        Write residual on file for later use in calculation of
C        linear response functions
C--------------------------------------------------------------
C
              CALL GPOPEN(LURV1E,FNRV1E,'UNKNOWN',' ','UNFORMATTED',
     &                    IDUMMY,.FALSE.)
              CALL GPOPEN(LURV1D,FNRV1D,'UNKNOWN',' ','UNFORMATTED',
     &                    IDUMMY,.FALSE.)
              CALL GPOPEN(LURV2E,FNRV2E,'UNKNOWN',' ','UNFORMATTED',
     &                    IDUMMY,.FALSE.)
              CALL GPOPEN(LURV2D,FNRV2D,'UNKNOWN',' ','UNFORMATTED',
     &                    IDUMMY,.FALSE.)
C             find the right positions
              CALL SO_FNDPOS(LABEL,EIVAL(IEXCI),LURV1E)
              CALL SO_FNDPOS(LABEL,EIVAL(IEXCI),LURV1D)
              CALL SO_FNDPOS(LABEL,EIVAL(IEXCI),LURV2E)
              CALL SO_FNDPOS(LABEL,EIVAL(IEXCI),LURV2D)
C
              CALL SO_WRTVE(WORK(KRESI1E),LRESI1E,ISYMTR,LABEL,
     &                      EIVAL(IEXCI),LURV1E)
              CALL SO_WRTVE(WORK(KRESI1D),LRESI1D,ISYMTR,LABEL,
     &                      EIVAL(IEXCI),LURV1D)
              CALL SO_WRTVE(WORK(KRESI2E),LRESI2E,ISYMTR,LABEL,
     &                      EIVAL(IEXCI),LURV2E)
              CALL SO_WRTVE(WORK(KRESI2D),LRESI2D,ISYMTR,LABEL,
     &                      EIVAL(IEXCI),LURV2D)
C
              CALL GPCLOSE(LURV1E,'KEEP')
              CALL GPCLOSE(LURV1D,'KEEP')
              CALL GPCLOSE(LURV2E,'KEEP')
              CALL GPCLOSE(LURV2D,'KEEP')
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
      CALL SO_ORTH_TRN(LTYPE,NOLDTR,NNEWTR,NLINDP,ISYMTR,DENSIJ,LDENSIJ,
     &                 DENSAB,LDENSAB,WORK,LWORK)
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
      CALL SO_CLOSE(LUTR2E,FNTR2E,'KEEP')
      CALL SO_CLOSE(LUTR2D,FNTR2D,'KEEP')
      CALL SO_CLOSE(LURS2E,FNRS2E,'KEEP')
      CALL SO_CLOSE(LURS2D,FNRS2D,'KEEP')
C
      IF(LTYPE.EQ.'LINEAR') THEN
              CALL SO_CLOSE(LUGPV1,FNGPV1,'KEEP')
              CALL SO_CLOSE(LUGPV2,FNGPV2,'KEEP')
      ENDIF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_TRIAL2')
C
      RETURN
C
 9001 FORMAT(/,1X,I3,'. excitation, norm of residual is: ',1P,D15.8,/)
C
      END
