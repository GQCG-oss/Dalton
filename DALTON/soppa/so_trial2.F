C
C  /* Deck so_trial2 */
      SUBROUTINE SO_TRIAL2(MODEL,LTYPE,NONEWT,NOLDTR,NNEWTR,NLINDP,
     &                     EIVAL,LEIVAL,RESINM,LRESINM,CONV,LCONV,
     &                     NCONV,ISYMTR,IMAGPROP,
     &                     NEXCI,DENSIJ,LDENSIJ,
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
      use so_info, only: so_has_doubles, so_singles_second,
     &                   sop_conv_thresh, sop_stat_trh,
     &                   sop_use_seller
C
#include "implicit.h"
#include "priunit.h"
C
#include "soppinf.h"
#include "ccsdsym.h"
#include "iratdef.h"
C
      PARAMETER   (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER   (D100 = 100.0D0)
C
      LOGICAL     NONEWT
C
      CHARACTER*3 YES, NO
      CHARACTER*3 CONV(LCONV)
      CHARACTER*6 LTYPE
      CHARACTER*8 LABEL
      CHARACTER*5 MODEL
C
      DIMENSION   EIVAL(LEIVAL), RESINM(LRESINM)
      DIMENSION   DENSIJ(LDENSIJ), DENSAB(LDENSAB)
      DIMENSION   WORK(LWORK)
      LOGICAL     IMAGPROP, DOUBLES, UNIT_S, STATIC
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
C---------------------------------------------
C     Does this method have doubles aplitudes?
C---------------------------------------------
      DOUBLES = SO_HAS_DOUBLES(MODEL)
C     S matrix is a unit matrix if we don't have second order singles
C     contributions
      UNIT_S = .NOT. SO_SINGLES_SECOND(MODEL)
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
      LSDIA1  = LRESO1D
      IF(DOUBLES)THEN
         LRES2E  = N2P2HOP(ISYMTR)
         LTR2E   = LRES2E
         LRES2D  = N2P2HOP(ISYMTR)
         LTR2D   = LRES2D
         LRESI2E = N2P2HOP(ISYMTR)
         LRESI2D = N2P2HOP(ISYMTR)
         LRESO2E = N2P2HOP(ISYMTR)
         LEDIA2  = LRESO2E
         LRESO2D = N2P2HOP(ISYMTR)
      ELSE
C   Setting all doubles lengths to zero, should allow us to keep
C   the following allocations unmodified
         LRES2E  = 0
         LTR2E   = 0
         LRES2D  = 0
         LTR2D   = 0
         LRESI2E = 0
         LRESI2D = 0
         LRESO2E = 0
         LEDIA2  = 0
         LRESO2D = 0
      ENDIF
C
      IF (LTYPE.EQ.'LINEAR') THEN
         LGPVC1H = NT1AM(ISYMTR)
         IF(DOUBLES)THEN
            LGPVC2H = N2P2HOP(ISYMTR)
         ELSE
            LGPVC2H = 0
         END IF
         LGPVC1  = LGPVC1H
         LGPVC2  = LGPVC2H
      ELSE
         LGPVC1H = 0
         LGPVC2H = 0
         LGPVC1  = 0
         LGPVC2  = 0
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
C     We need to transform trial-vectors with S (S*Tr) , grap some space
      IF (UNIT_S) THEN
C        Only need to change sign on TR1D, done in place
         KTMP1E = KRESO1E
         KTMP1D = KRESO1D
      ELSE
C        Otherwise we borrow the space from the result vectors
         KTMP1E = KRES1E
         KTMP1D = KRES1D
      ENDIF
C
      IF(LTYPE.EQ.'LINEAR') THEN
         KGPVC1  = KEND1B
         KGPVC2  = KGPVC1 + LGPVC1
         KEND1C  = KGPVC2 + LGPVC2
      ELSE
         KGPVC1  = KEND1B
         KGPVC2  = KEND1B
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
      IF(DOUBLES) THEN
         CALL SO_OPEN(LURS2E,FNRS2E,LRES2E)
         CALL SO_OPEN(LURS2D,FNRS2D,LRES2D)
         CALL SO_OPEN(LUTR2E,FNTR2E,LTR2E)
         CALL SO_OPEN(LUTR2D,FNTR2D,LTR2D)
      ENDIF
C
      IF(LTYPE.EQ.'LINEAR') THEN
         CALL SO_OPEN(LUGPV1,FNGPV1,LGPVC1)
         IF(DOUBLES) CALL SO_OPEN(LUGPV2,FNGPV2,LGPVC2)
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
         CALL SO_READ(WORK(KTMP1E),LRESO1E,LUTR1E,FNTR1E,IEXCI)
         CALL SO_READ(WORK(KTMP1D),LRESO1D,LUTR1D,FNTR1D,IEXCI)
         IF (UNIT_S) THEN
C           S(1) == S(0) for RPA, that is simple
            CALL DSCAL(LRESO1D,-ONE,WORK(KRESO1D),1) !KRESO1D = KTMP1D
         ELSE
C           S(2) for doubles, that require some work
            CALL SO_RES_O(WORK(KRESO1E),LRESO1E,WORK(KRESO1D),LRESO1D,
     &                    WORK(KTMP1E), LRESO1E,WORK(KTMP1D) ,LRESO1D,
     &                    DENSIJ,LDENSIJ,DENSAB,LDENSAB,
C           ISYRES == ISYMTR always for abelian symmetry, totally symmetry S..
C           right?
     &                    ISYMTR,ISYMTR)
         ENDIF
C  Warning: For model != aorpa, the space for the result-vectors
C  is used to perform the above transformation.

         CALL SO_READ(WORK(KRES1E), LRES1E, LURS1E,FNRS1E,IEXCI)
         CALL SO_READ(WORK(KRES1D), LRES1D, LURS1D,FNRS1D,IEXCI)
         IF(DOUBLES)THEN
            CALL SO_READ(WORK(KRES2E), LRES2E, LURS2E,FNRS2E,IEXCI)
            CALL SO_READ(WORK(KRES2D), LRES2D, LURS2D,FNRS2D,IEXCI)
C
C           S(0) for doubles, factor of -1 on D part
            CALL SO_READ(WORK(KRESO2E),LRESO2E,LUTR2E,FNTR2E,IEXCI)
            CALL SO_READ(WORK(KRESO2D),LRESO2D,LUTR2D,FNTR2D,IEXCI)
            CALL DSCAL(LRESO2D,-ONE,WORK(KRESO2D),1)
         ENDIF
C
         IF(LTYPE.EQ.'LINEAR') THEN
            CALL SO_READ(WORK(KGPVC1),LGPVC1,LUGPV1,FNGPV1,1)
            IF(DOUBLES)THEN
                 CALL SO_READ(WORK(KGPVC2),LGPVC2,LUGPV2,FNGPV2,1)
C  Property gradient on file is in the wrong basis
                 IF (.NOT.TRIPLET) 
     &              CALL SO_TMLTR(WORK(KGPVC2),HALF,ISYMTR)
            ENDIF
C
            IF(IPRSOP.GT.10) THEN
C
               CALL AROUND("Left GP vector before SO_RESIDUAL")
               CALL OUTPUT(WORK(KGPVC1),1,LGPVC1,1,1,LGPVC1H,1,1,LUPRI)
               IF(DOUBLES)
     &         CALL OUTPUT(WORK(KGPVC2),1,LGPVC2,1,1,LGPVC2H,1,1,LUPRI)
C
            ENDIF
C
         ENDIF
C
C-------------------------------------------------------------
C        Calculate residual vector and norm of result vectors.
C-------------------------------------------------------------
C
         CALL SO_RESIDUAL(DOUBLES,LTYPE,RESINM(IEXCI),ISYMTR,
     &                    WORK(KRESI1E),LRESI1E,WORK(KRESI1D),LRESI1D,
     &                    WORK(KRESI2E),LRESI2E,WORK(KRESI2D),LRESI2D,
     &                    EIVAL(IEXCI),
     &                    WORK(KRES1E),LRES1E,WORK(KRES1D),LRES1D,
     &                    WORK(KRES2E),LRES2E,WORK(KRES2D),LRES2D,
     &                    WORK(KRESO1E),LRESO1E,WORK(KRESO1D),LRESO1D,
     &                    WORK(KRESO2E),LRESO2E,
     &                    WORK(KRESO2D),LRESO2D,
     &                    WORK(KGPVC1),LGPVC1,WORK(KGPVC2),LGPVC2,
     &                    ENORM,IMAGPROP)
C
C-----------------------------------------------
C        Determine the threshold the excitation.
C-----------------------------------------------
C
         THR = SOP_CONV_THRESH
C
         IF ((LTYPE.EQ.'EXCITA') .and.
     &       (IEXCI .GT. (NEXCI - NEXCI2(ISYMTR)))
     &      ) THR = THREX2
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
C            GO TO 100
C
C--------------------------------------------
C        If not converged but last iteration.
C--------------------------------------------
C
         ELSE IF (NONEWT) THEN
C
            CONV(IEXCI) = NO
C
C            GO TO 100
C
C------------------------------------------------
C        If not converged and not last iteration.
C------------------------------------------------
C
         ELSE
C
C------------------------------------------------------------------
C           If imaginary excitation increase NCONV in order to stop
C           iterations if all other excitations are either converged
C           or imaginary, but create a new trial vector in case all
C           excitations haven't converged. The 'i' label in case of
C           imaginary excitation is assigned in SO_ORDEIG
C-------------------------------------------------------------------
CPi 10.08.16
            IF ( CONV(IEXCI) .EQ. '  i') THEN
C
               NCONV = NCONV + 1
C
            ELSE
C
               CONV(IEXCI) = NO
C
            END IF
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
     &                   IDUMMY,.FALSE.)
            REWIND LUDIAG
C
            READ(LUDIAG) ( WORK(KEDIA1+I-1), I = 1,LEDIA1)
            IF (DOUBLES)THEN
               READ(LUDIAG) ( WORK(KEDIA2+I-1), I = 1,LEDIA2)
            ENDIF
            IF(MODEL.NE.'AORPA')THEN
               READ(LUDIAG) ( WORK(KSDIA1+I-1), I = 1,LSDIA1)
            ENDIF
C
            CALL GPCLOSE (LUDIAG,'KEEP')
C
C------------------------------------------------------------
C           Calculate raw new trial vector and write to file.
C------------------------------------------------------------
C
            IF(MODEL.EQ.'AORPA') THEN
C Use special RPA routine, since for RPA "S" is not explicitly calculated
C (and it would be a waste to do so)
              CALL RP_NEWTRIAL(NNEWTR,NOLDTR,WORK(KTR1E),LTR1E,
     &                         WORK(KTR1D),LTR1D,EIVAL(IEXCI),
     &                         WORK(KEDIA1),LEDIA1,
     &                         WORK(KRESI1E),LRESI1E,
     &                         WORK(KRESI1D),LRESI1D)
            ELSE
              CALL SO_NEWTRIAL(DOUBLES,NNEWTR,NOLDTR,
     &                         WORK(KTR1E),LTR1E,WORK(KTR1D),LTR1D,
     &                         WORK(KTR2E),LTR2E,WORK(KTR2D),LTR2D,
     &                         EIVAL(IEXCI),WORK(KEDIA1),LEDIA1,
     &                         WORK(KEDIA2),LEDIA2,WORK(KSDIA1),LSDIA1,
     &                         WORK(KRESI1E),LRESI1E,
     &                         WORK(KRESI1D),LRESI1D,
     &                         WORK(KRESI2E),LRESI2E,WORK(KRESI2D),
     &                         LRESI2D)
            ENDIF
C
         END IF
C
C
C---------------------------------------------------------------------
C        Write residual on file for later use if Seller's formula
C        will be used in the calculation of linear response functions.
C        We will also save the solution vectors. 
C        For methods with non-diagonal S-matrices, we need to re-read
C        these.
C        We must also scale the back to the unnormalized vectors
C---------------------------------------------------------------------
C
         IF (SOP_USE_SELLER.AND.(CONV(IEXCI).EQ.YES)) THEN
            STATIC = ABS(EIVAL(IEXCI)).LT.sop_stat_trh
C
            ENORMINV = 1.0D0/ENORM
            LURV1E = -1
            CALL GPOPEN(LURV1E,FNRV1E,'UNKNOWN','DIRECT','UNFORMATTED',
     &                    IRAT*LRESI1E,OLDDX)
            CALL DSCAL(LRESI1E,ENORMINV,WORK(KRESI1E),1)
            CALL SO_WRTVE(WORK(KRESI1E),LRESI1E,ISYMTR,LABEL,
     &                    EIVAL(IEXCI),LURV1E)
            CALL GPCLOSE(LURV1E,'KEEP')

            IF (.NOT.UNIT_S) THEN
C -- Re-read the trial-vector from file
               CALL SO_READ(WORK(KRESO1E), LRESO1E, LUTR1E, FNTR1E,
     &                      IEXCI)
            END IF
            LUSV1E = -1
            CALL GPOPEN(LUSV1E,FNSV1E,'UNKNOWN','DIRECT','UNFORMATTED',
     &                  IRAT*LRESO1E,OLDDX)
            CALL DSCAL(LRESO1E,ENORMINV,WORK(KRESO1E),1)
            CALL SO_WRTVE(WORK(KRESO1E),LRESO1E,ISYMTR,LABEL,
     &                       EIVAL(IEXCI),LUSV1E)
            CALL GPCLOSE(LUSV1E,'KEEP')

C
            IF (.NOT.STATIC) THEN
               LURV1D = -1
               CALL GPOPEN(LURV1D,FNRV1D,'UNKNOWN','DIRECT',
     &                    'UNFORMATTED',
     &                    IRAT*LRESI1D,OLDDX)
               CALL DSCAL(LRESI1D,ENORMINV,WORK(KRESI1D),1)
               CALL SO_WRTVE(WORK(KRESI1D),LRESI1D,ISYMTR,LABEL,
     &                      EIVAL(IEXCI),LURV1D)
               CALL GPCLOSE(LURV1D,'KEEP')
C
               IF(UNIT_S) THEN
                  DFACT = -ENORMINV
               ELSE 
                  CALL SO_READ(WORK(KRESO1D), LRESO1D, LUTR1D, FNTR1D,
     &                         IEXCI)
                  DFACT = ENORMINV
               END IF
C
               LUSV1D = -1
               CALL GPOPEN(LUSV1D,FNSV1D,'UNKNOWN','DIRECT',
     &                     'UNFORMATTED',
     &                     IRAT*LRESO1D,OLDDX)
               CALL DSCAL(LRESO1D,DFACT,WORK(KRESO1D),1)
               CALL SO_WRTVE(WORK(KRESO1D),LRESO1D,ISYMTR,LABEL,
     &                       EIVAL(IEXCI),LUSV1D)
               CALL GPCLOSE(LUSV1D,'KEEP')

            END IF
C
C---------------------------------------------------------------------
C           For doubles we might as well also save the response vector
C---------------------------------------------------------------------
C
            IF(DOUBLES)THEN
               LURV2E = -1
               CALL GPOPEN(LURV2E,FNRV2E,'UNKNOWN','DIRECT',
     &                    'UNFORMATTED',
     &                    IRAT*LRESI2E,OLDDX)
C              find the right positions
               IF (TRIPLET) THEN
                  CALL DSCAL(LRESI2E,ENORMINV,WORK(KRESI2E),1)
               ELSE
C              Put the residual in the other basis...
C              We need to scale by one half
                  CALL DSCAL(LRESI2E,ENORMINV*0.5D0,WORK(KRESI2E),1)
                  CALL CCSD_TCMEPKX(WORK(KRESI2E),2.0D0,ISYMTR)
               END IF
               CALL SO_WRTVE(WORK(KRESI2E),LRESI2E,ISYMTR,LABEL,
     &                      EIVAL(IEXCI),LURV2E)
               CALL GPCLOSE(LURV2E,'KEEP')
C
               LUSV2E = -1
               CALL GPOPEN(LUSV2E,FNSV2E,'UNKNOWN','DIRECT',
     &                    'UNFORMATTED',
     &                    IRAT*LRESO2E,OLDDX)
C              find the right positions
               CALL DSCAL(LRESO2E,ENORMINV,WORK(KRESO2E),1)
               CALL SO_WRTVE(WORK(KRESO2E),LRESO2E,ISYMTR,LABEL,
     &                       EIVAL(IEXCI),LUSV2E)
               CALL GPCLOSE(LUSV2E,'KEEP')
C
            END IF
            IF(DOUBLES.AND..NOT.STATIC) THEN
               LURV2D = -1
               CALL GPOPEN(LURV2D,FNRV2D,'UNKNOWN','DIRECT',
     &                    'UNFORMATTED',
     &                    IRAT*LRESI2D,OLDDX)
               IF (TRIPLET) THEN
                  CALL DSCAL(LRESI2D,ENORMINV,WORK(KRESI2D),1)
               ELSE
C                 Put the residual in the other basis...
C                 We need to scale by one half               
                  CALL DSCAL(LRESI2D,ENORMINV*0.5D0,WORK(KRESI2D),1)
                  CALL CCSD_TCMEPKX(WORK(KRESI2D),2.0D0,ISYMTR)
               END IF
               CALL SO_WRTVE(WORK(KRESI2D),LRESI2D,ISYMTR,LABEL,
     &                      EIVAL(IEXCI),LURV2D)
               CALL GPCLOSE(LURV2D,'KEEP')
C
               LUSV2D = -1
               CALL GPOPEN(LUSV2D,FNSV2D,'UNKNOWN','DIRECT',
     &                    'UNFORMATTED',
     &                    IRAT*LRESO2D,OLDDX)
C              find the right positions
               CALL DSCAL(LRESO2D,-ENORMINV,WORK(KRESO2D),1)
               CALL SO_WRTVE(WORK(KRESO2D),LRESO2D,ISYMTR,LABEL,
     &                      EIVAL(IEXCI),LUSV2D)
               CALL GPCLOSE(LUSV2D,'KEEP')
            ENDIF
         ENDIF
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
      IF(UNIT_S)THEN
         CALL RP_ORTH_TRN(LTYPE,NOLDTR,NNEWTR,NLINDP,ISYMTR,WORK,LWORK)
      ELSE
         CALL SO_ORTH_TRN(DOUBLES,LTYPE,NOLDTR,NNEWTR,NLINDP,ISYMTR,
     &                    DENSIJ,LDENSIJ,
     &                    DENSAB,LDENSAB,WORK,LWORK)
      ENDIF
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
      IF(DOUBLES)THEN
         CALL SO_CLOSE(LUTR2E,FNTR2E,'KEEP')
         CALL SO_CLOSE(LUTR2D,FNTR2D,'KEEP')
         CALL SO_CLOSE(LURS2E,FNRS2E,'KEEP')
         CALL SO_CLOSE(LURS2D,FNRS2D,'KEEP')
      ENDIF
C
      IF(LTYPE.EQ.'LINEAR') THEN
         CALL SO_CLOSE(LUGPV1,FNGPV1,'KEEP')
         IF(DOUBLES) CALL SO_CLOSE(LUGPV2,FNGPV2,'KEEP')
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
