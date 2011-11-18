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
#include "single.h"
C
C File: abaop2.F
C
C 971020-vebjornb: Just an extension of abaopt.F which is becoming awkwardly large.
C
C  /* Deck linsrc */
      SUBROUTINE LINSRC(ICRD,MCRD,GRAD,GRADOL,STEP,STEPOL,
     &     TMPVEC,TMPVC2,ACTIVE,EMOD)
C
C     This routine calculates the step that should be taken to obtain
C     the next geometry.
C
#include "implicit.h"
#include "mxcent.h"
#include "priunit.h"
#include "optinf.h"
#include "cbiwlk.h"
#include "trkoor.h"
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      PARAMETER (D0 = 0.0D0)
      LOGICAL ACTIVE
      DIMENSION GRAD(MCRD), GRADOL(MCRD), STEP(MCRD), STEPOL(MCRD)
      DIMENSION TMPVEC(MCRD), TMPVC2(MCRD)
      ACTIVE = .TRUE.
      CALL DZERO(TMPVEC,MCRD)
      CALL DZERO(STEP,MCRD)
      STNM = SQRT(DDOT(ICRD,STEPOL,1,STEPOL,1))
      DO 10 I = 1, ICRD
         TMPVEC(I) = -STEPOL(I)/STNM
 10   CONTINUE
      GRAD0 = DDOT(ICRD,GRADOL,1,TMPVEC,1)
      GRAD1 = DDOT(ICRD,GRAD,1,TMPVEC,1)
      IF ((ABS(GRAD0) .LT. 1.0D-6) .OR. (ABS(GRAD1) .LT. 1.0D-6)) THEN
         IF (IPRINT .GE. IPRDBG) THEN
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) 'Too close to minimum, line search skipped.'
         END IF
         ACTIVE = .FALSE.
         RETURN
      END IF
C
C     The energy and gradient from the last and current point, are
C     fitted to a quartic polynomial.
C
      CA = ABS(5.0D0*GRAD0 - 2.0D0*GRAD1 + 3.0D0*ERGOLD - 3.0D0*ENERGY)
      CB = GRAD1 + GRAD0 + 2.0D0*ERGOLD - 2.0D0*ENERGY - 2.0D0*CA
      CC = ENERGY - ERGOLD - GRAD0 - CA - CB
      CD = GRAD0
      CE = ERGOLD
C
C     The line search methdod is only used if the minimum lies
C     between the two points.
C
      IF ((GRAD0 .LT. D0) .AND. (GRAD1 .GT. D0)) THEN
         THRG = 1.0D-5*MIN(ABS(GRAD0),ABS(GRAD1))
         CRDA = D0
         CRDB = 1.0D0
         GRDA = GRAD0
         GRDB = GRAD1
         ISAFE = 0
 15      CONTINUE
         ISAFE = ISAFE + 1
         IF (ISAFE .GE. 200) THEN
            ACTIVE = .FALSE.
            IF (IPRINT .GE. IPRDBG) THEN
               WRITE(LUPRI,*)
               WRITE(LUPRI,*) 'Line search failed, ignoring.'
            END IF
            RETURN
         END IF
         CRDC = CRDA + (CRDB-CRDA)*
     &        MAX(0.1D0,MIN(0.9D0,ABS(GRDA/GRDB)*0.5D0))
         GRDC = ((4.0D0*CA*CRDC + 3.0D0*CB)*CRDC + 2.0D0*CC)*CRDC + CD
         IF (ABS(GRDC) .GT. THRG) THEN
            IF (GRDC .GT. D0) THEN
               CRDB = CRDC
               GRDB = GRDC
            ELSE
               CRDA = CRDC
               GRDA = GRDC
            END IF
            GOTO 15
         END IF
C
C     If the line search ends up almost back at the previous point, we do
C     not trust it and simply discard it.
C
         IF (CRDC .LT. 0.15D0) CRDC = 1.D0
C
         DO 20 I = 1, ICRD
            STEP(I) = STEPOL(I)*(1.0D0-CRDC)
            GRAD(I) = GRADOL(I) + (GRAD(I)-GRADOL(I))*CRDC
 20      CONTINUE
         EMOD = (((CA*CRDC + CB)*CRDC + CC)*CRDC + CD)*CRDC + CE
         IF (IPRINT .GE. IPRDBG) THEN
            CALL HEADER('Interpolated step',-1)
            CALL OUTPUT(STEP,1,1,1,ICRD,1,MCRD,1,LUPRI)
            CALL HEADER('Interpolated gradient',-1)
            CALL OUTPUT(GRAD,1,1,1,ICRD,1,MCRD,1,LUPRI)
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) 'Interpolated energy: ',EMOD
         END IF
      ELSE
         IF (IPRINT .GE. IPRDBG) THEN
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) 'Line search skipped.'
         END IF
         ACTIVE = .FALSE.
         RETURN
      END IF
      RETURN
      END
      
C  /* Deck upgdst */
      SUBROUTINE UPGDST(ICRD,MXRCRD,GRDARR,STPARR,GRAD,STEP)
C
C     Updates arrays containing steps and gradients from
C     previous iterations.
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "optinf.h"
      DIMENSION GRDARR(MXRCRD,25), STPARR(MXRCRD,25)
      DIMENSION GRAD(*), STEP(*)
      IF (RSTARR) THEN
         CALL DZERO(GRDARR,25*MXRCRD)
         CALL DZERO(STPARR,25*MXRCRD)
         KEPTIT = 0
         RSTARR = .FALSE.
      ELSE
         DO 10 I = MIN(24,KEPTIT),1,-1
            DO 20 J = 1, ICRD
               GRDARR(J,I+1) = GRDARR(J,I)
               STPARR(J,I+1) = STPARR(J,I) - STEP(J)
 20         CONTINUE
 10      CONTINUE
         DO 30 I = 1, ICRD
            GRDARR(I,1) = GRAD(I)
            STPARR(I,1) = - STEP(I)
 30      CONTINUE
         KEPTIT = MIN(25,KEPTIT+1)
      END IF
      RETURN
      END

C  /* Deck gdistp */
      SUBROUTINE GDISTP(NCRD,ICRD,MXRCRD,MX2CRD,STEP,GRAD,HESS,HESINV,
     &     TMPMAT,TMPMT2,TMPMT3,TMPMT4,GRDARR,STPARR)
C
C     Updates arrays containing steps and gradients from
C     previous iterations.
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "optinf.h"
      LOGICAL STPSCL
      DIMENSION STEP(NCRD), GRAD(NCRD)
      DIMENSION HESS(NCRD,NCRD), HESINV(ICRD,ICRD)
      DIMENSION TMPMAT(MX2CRD,MX2CRD),TMPMT2(MX2CRD*MX2CRD)
      DIMENSION TMPMT3(MX2CRD*MX2CRD),TMPMT4(MX2CRD*MX2CRD)
      DIMENSION GRDARR(MXRCRD,25), STPARR(MXRCRD,25)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      CALL DZERO(TMPMAT,MX2CRD*MX2CRD)
      CALL DZERO(TMPMT2,MX2CRD*MX2CRD)
      CALL DZERO(TMPMT3,MX2CRD*MX2CRD)
      CALL DZERO(TMPMT4,MX2CRD*MX2CRD)
C
      STPSCL = .TRUE.
C
C     First we have to construct the DIIS matrix from the
C     scalar products of the gradients (we use the gradients as
C     error vectors).
C
      IDIM = KEPTIT + 2
      DO 10 I = 1, KEPTIT
         DO 15 J = I, KEPTIT
            TMPMAT(I+1,J+1) = DDOT(ICRD,GRDARR(1,I),1,GRDARR(1,J),1)
            TMPMAT(J+1,I+1) = TMPMAT(I+1,J+1)
 15      CONTINUE
 10   CONTINUE
      TMPMAT(1,1) = DDOT(ICRD,GRAD(1),1,GRAD(1),1)
      DO 17 I = 1, KEPTIT
         TMPMAT(1,I+1) = DDOT(ICRD,GRAD(1),1,GRDARR(1,I),1)
         TMPMAT(I+1,1) = TMPMAT(1,I+1)
         TMPMAT(I+1,IDIM) = D1
         TMPMAT(IDIM,I+1) = D1
 17   CONTINUE
      TMPMAT(IDIM,IDIM) = D0
      TMPMAT(1,IDIM) = D1
      TMPMAT(IDIM,1) = D1
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('DIIS matrix',-1)
         CALL OUTPUT(TMPMAT,1,IDIM,1,IDIM,
     &        MX2CRD,MX2CRD,1,LUPRI)
      END IF
 25   CONTINUE
      CALL DZERO(TMPMT2,IDIM*IDIM)
      CALL DZERO(TMPMT3,IDIM*IDIM)
      DO 30 J = 1, IDIM
         DO 32 I = 1, IDIM
            TMPMT2(I+(J-1)*IDIM) = TMPMAT(I,J)
 32      CONTINUE
 30   CONTINUE
      CALL DSITSP(IDIM,TMPMT2,TMPMT3)
      CALL DUNIT(TMPMT2,IDIM)
      CALL JACO(TMPMT3,TMPMT2,IDIM,IDIM,IDIM,TMPMT4(1),
     &     TMPMT4(MX2CRD+1))
      CALL DZERO(TMPMT4,MX2CRD*MX2CRD)
      IZER = 0
      EMAX = D0
      DO 34 J = 1, IDIM-1
         TMPMT4(J) = TMPMT3(J*(J+1)/2)
         IF (ABS(TMPMT4(J)) .GT. EMAX) EMAX = ABS(TMPMT4(J))
 34   CONTINUE
      TMPMT4(IDIM) = TMPMT3(IDIM*(IDIM+1)/2)
      DO 35 J = 1, IDIM
         IF (ABS(TMPMT4(J)) .LE. EMAX*1.0D-2) THEN
            IZER = IZER + 1
            LSTZER = J
         END IF
 35   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Initial eigenvalues of DIIS matrix',-1)
         CALL OUTPUT(TMPMT4,1,1,1,IDIM,1,IDIM,1,LUPRI)
      END IF
      IF ((IZER .GT. 0) .AND. (IDIM .GT. 2)) THEN
C
C     We try to remove only the components causing linearity
C     (one at a time)
C
         IFAC = (LSTZER-1)*IDIM
         CMPMAX = D0
         DO 360 K = 1, IDIM-1
            IF (ABS(TMPMT2(K+IFAC)) .GT. CMPMAX) THEN
               CMPMAX = ABS(TMPMT2(K+IFAC))
               MAX1 = K
            END IF
 360     CONTINUE
         CMPMX2 = D0
         DO 361 K = 1, IDIM-1
            IF (K .NE. MAX1) THEN
               IF (ABS(TMPMT2(K+IFAC)) .GT. CMPMX2) THEN
                  CMPMX2 = ABS(TMPMT2(K+IFAC))
                  MAX2 = K
               END IF
            END IF
 361     CONTINUE
         NREMV = MAX(MAX1,MAX2)
         DO 370 I = NREMV, IDIM-1
            DO 371 J = 1, IDIM
               TMPMAT(I,J) = TMPMAT(I+1,J)
 371        CONTINUE
 370     CONTINUE
         DO 380 I = NREMV, IDIM-1
            DO 381 J = 1, IDIM-1
               TMPMAT(J,I) = TMPMAT(J,I+1)
 381        CONTINUE
 380     CONTINUE
         IDIM = IDIM - 1
         GOTO 25
      END IF
      DO 40 J = 1, IDIM
         DO 42 I = 1, IDIM
            TMPMAT(I,J) = TMPMT2(I+(J-1)*IDIM)
 42      CONTINUE
 40   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Eigenvalues of DIIS matrix',-1)
         CALL OUTPUT(TMPMT4,1,1,1,IDIM,1,IDIM,1,LUPRI)
         CALL HEADER('Eigenvectors of DIIS matrix',-1)
         CALL OUTPUT(TMPMAT,1,IDIM,1,IDIM,MX2CRD,MX2CRD,1,LUPRI)
      END IF
C
C     We find the inverse of the DIIS matrix
C
      CALL DZERO(TMPMT2,MX2CRD*MX2CRD)
      CALL DZERO(TMPMT3,MX2CRD*MX2CRD)
      DO 45 I = 1, IDIM
         IF (ABS(TMPMT4(I)) .LE. 1.0D-6) THEN
            TMPMT4(I) = D0
         ELSE
            TMPMT4(I)= 1/TMPMT4(I)
         END IF
 45   CONTINUE
      DO 50 J = 1, IDIM
         DO 52 I = 1, IDIM
            TMPMT2(I+(J-1)*IDIM) = TMPMT4(I)*TMPMAT(J,I)
            TMPMT3(I+(J-1)*IDIM) = TMPMAT(I,J)
 52      CONTINUE
 50   CONTINUE
      CALL DZERO(TMPMAT,MX2CRD*MX2CRD)
      DO 60 J = 1, IDIM
         DO 62 I = 1, IDIM
            DO 64 K = 1, IDIM
               TMPMAT(I,J) = TMPMAT(I,J) + 
     &              TMPMT3(I+(K-1)*IDIM)*TMPMT2(K+(J-1)*IDIM)
 64         CONTINUE
 62      CONTINUE
 60   CONTINUE
      CALL DZERO(TMPMT4,MX2CRD*MX2CRD)
      FAC = D0
      DO 70 I = 1, IDIM-1
         TMPMT4(I) = TMPMAT(I,IDIM)
         FAC = FAC + TMPMT4(I)
 70   CONTINUE
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Inverse of DIIS matrix',-1)
         CALL OUTPUT(TMPMAT,1,IDIM,1,IDIM,MX2CRD,MX2CRD,1,LUPRI)
         CALL HEADER('DIIS coefficients',-1)
         CALL OUTPUT(TMPMT4,1,1,1,IDIM-1,1,MX2CRD,1,LUPRI)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Sum of coefficients: ',FAC
      END IF
      CALL DZERO(TMPMT3,MX2CRD)
      CALL DZERO(STEP,NCRD)
      DO 80 I = 1, ICRD
         TMPMT3(I) = TMPMT4(1)*GRAD(I)
 80   CONTINUE
      DO 82 I = 1, IDIM-2
         DO 84 J = 1, ICRD
            STEP(J) = STEP(J) + TMPMT4(I+1)*STPARR(J,I)
 84      CONTINUE
 82   CONTINUE
C
C     If the step is too large, we simply restrict each element
C     to be below the trust radius.
C     
      IF (.NOT. STPSCL) THEN
         DO 85 I = 1, ICRD
            IF (ABS(STEP(I)) .GT. TRSTRA)
     &           STEP(I) = SIGN(TRSTRA,STEP(I))
 85      CONTINUE
      ELSE
C
C     Alternatively we restrivt the step norm to be equal or less
C     than the trust radius.
C
         STPNRM = SQRT(DDOT(ICRD,STEP,1,STEP,1))
         IF (STPNRM .GT. TRSTRA) THEN
            FAC = TRSTRA/STPNRM
            IF (IPRINT .GE. IPRDBG) THEN
               WRITE(LUPRI,*)
               WRITE(LUPRI,*)
     &              'DIIS-step too long, scaling by factor:', FAC
               WRITE(LUPRI,*)
            END IF
            DO 86 I = 1, ICRD
               STEP(I) = STEP(I)*FAC
 86         CONTINUE
         END IF
      END IF
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Interpolated step',-1)
         CALL OUTPUT(STEP,1,1,1,ICRD,1,NCRD,1,LUPRI)
         CALL HEADER('Interpolated gradient',-1)
         CALL OUTPUT(TMPMT3,1,1,1,ICRD,1,MX2CRD,1,LUPRI)
         CALL HEADER('Inverse of Hessian',-1)
         CALL OUTPUT(HESINV,1,ICRD,1,ICRD,ICRD,ICRD,1,LUPRI)
      END IF
      CALL DZERO(TMPMT2,ICRD)
      DO 100 I = 1, ICRD
         DO 102 J = 1, ICRD
            TMPMT2(I) = TMPMT2(I) + HESINV(I,J)*TMPMT3(J)
 102     CONTINUE
         STEP(I) = STEP(I) - TMPMT2(I)
 100  CONTINUE
      IF (IPRINT .GE. IPRMAX) THEN
         CALL HEADER('Relaxation step',-1)
         CALL OUTPUT(TMPMT2,1,1,1,ICRD,1,NCRD,1,LUPRI)
         CALL HEADER('Total DIIS step',-1)
         CALL OUTPUT(STEP,1,1,1,ICRD,1,NCRD,1,LUPRI)
      END IF
C
C     If the step is too large, we simply restrict each element
C     to be below the trust radius.
C     
      IF (.NOT. STPSCL) THEN
         DO 185 I = 1, ICRD
            IF (ABS(STEP(I)) .GT. TRSTRA)
     &           STEP(I) = SIGN(TRSTRA,STEP(I))
 185     CONTINUE
      ELSE
C     
C     Alternatively we restrivt the step norm to be equal or less
C     than the trust radius.
C     
         STPNRM = SQRT(DDOT(ICRD,STEP,1,STEP,1))
         IF (STPNRM .GT. TRSTRA) THEN
            FAC = TRSTRA/STPNRM
            IF (IPRINT .GE. IPRDBG) THEN
               WRITE(LUPRI,*)
               WRITE(LUPRI,*)
     &              'Relaxation step too long, scaling by factor:', FAC
               WRITE(LUPRI,*)
            END IF
            DO 186 I = 1, ICRD
               STEP(I) = STEP(I)*FAC
 186        CONTINUE
         END IF
      END IF
      ICNT = 1
      FAC = 1.0D0
 200  CONTINUE
      IF (ICNT .LE. 10) THEN
         CALL DZERO(TMPMT2,MX2CRD)
         SNDTRM = 0.0D0
         DO 210 I = 1, ICRD
            DO 215 J = 1, ICRD
               TMPMT2(I) = TMPMT2(I) + HESS(I,J)*STEP(J)
 215        CONTINUE
            SNDTRM = SNDTRM + TMPMT2(I)*STEP(I)
 210     CONTINUE
         ERGPRD = DDOT(ICRD,GRAD,1,STEP,1)
     &        + 0.5D0*SNDTRM
         IF (ERGPRD .GT. 0.0D0) THEN
            DO 220 I = 1, ICRD
               STEP(I) = STEP(I)*FAC
 220        CONTINUE
            FAC = 0.5D0*FAC
            ICNT = ICNT + 1
            GOTO 200
         END IF
      ELSE
         IF (IPRINT .GE. IPRMIN) THEN
            IF (REDINT .OR. DELINT) THEN
               CALL HEADER('DIIS-step in internal coordinates',-1) 
            ELSE
               CALL HEADER('DIIS-step',-1) 
            END IF
            CALL OUTPUT(STEP,1,1,1,ICRD,1,NCRD,1,LUPRI)
         END IF
         DO 230 I = 1, ICRD
            STEP(I) = -GRAD(I)
 230     CONTINUE
         CALL DZERO(TMPMT2,MX2CRD)
         SNDTRM = 0.0D0
         DO 240 I = 1, ICRD
            DO 245 J = 1, ICRD
               TMPMT2(I) = TMPMT2(I) + HESS(I,J)*STEP(J)
 245        CONTINUE
            SNDTRM = SNDTRM + TMPMT2(I)*STEP(I)
 240     CONTINUE
         ERGPRD = DDOT(ICRD,GRAD,1,STEP,1)
     &        + 0.5D0*SNDTRM
      END IF
      STPNRM = SQRT(DDOT(ICRD,STEP,1,STEP,1))
      RETURN
      END

C  /* Deck makimg */
      SUBROUTINE MAKIMG(ICRD,IPRJ,MCRD,EVAL,GRDDIA,STPDIA,IMODE,RESTOR)
C
C     Makes image function and sorts coordinates.
C     It is also used to restore and resort coordinates.
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "dummy.h"
      DIMENSION EVAL(MCRD), GRDDIA(MCRD), STPDIA(MCRD)
      SAVE ENEG, GNEG, ISORT
      LOGICAL RESTOR
      IF (.NOT. RESTOR) THEN
         ENEG = -EVAL(IMODE)
         GNEG = -GRDDIA(IMODE)
         EVAL(IMODE) = ENEG
         GRDDIA(IMODE) = DUMMY
         CALL ORDER(GRDDIA,EVAL,IPRJ,1)
         DO 10 I = 1, IPRJ
            IF ((EVAL(I).EQ.ENEG).AND.(GRDDIA(I).EQ.DUMMY)) ISORT = I
 10      CONTINUE
         GRDDIA(ISORT) = GNEG
      ELSE
         SNEG = STPDIA(ISORT)
         IF (ISORT .LT. IMODE) THEN
            DO 20 I = ISORT, IMODE-1
               EVAL(I)   = EVAL(I+1)
               GRDDIA(I) = GRDDIA(I+1)
               STPDIA(I) = STPDIA(I+1)
 20         CONTINUE
            EVAL(IMODE)   = -ENEG
            GRDDIA(IMODE) = -GNEG
            STPDIA(IMODE) = SNEG
         ELSE IF (IMODE .LT. ISORT) THEN
            DO 30 I = ISORT, IMODE+1,-1
               EVAL(I)   = EVAL(I-1)
               GRDDIA(I) = GRDDIA(I-1)
               STPDIA(I) = STPDIA(I-1)
 30         CONTINUE
            EVAL(IMODE)   = -ENEG
            GRDDIA(IMODE) = -GNEG
            STPDIA(IMODE) = SNEG
         ELSE
            EVAL(IMODE)   = -ENEG
            GRDDIA(IMODE) = -GNEG
         END IF
      END IF
      RETURN
      END

C  /* Deck fndmod */
      SUBROUTINE FNDMOD(INTERN,MXRCRD,EVEC,WILBMT,VECMOD,TMPVEC,
     &     TMPVC2,IMODE)
C
C     Makes sure we follow the same eigenvectormode throughout
C     a saddle point optimization. The mode is selected as the mode
C     (in Cartesian coordinates) with the largest overlap with
C     the previous mode. In redundant coordinates we have to check
C     that the selected mode corresponds to non-zero diagonal element.
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "optinf.h"
      DIMENSION EVEC(MXRCRD,MXRCRD), WILBMT(MXRCRD,MXCOOR)
      DIMENSION VECMOD(MXCOOR), TMPVEC(MXCOOR), TMPVC2(MXCOOR)
      LOGICAL INTERN, REMAIN
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      PARAMETER (D0 = 0.0D0)
      ICRD = NCART
      IF (INTERN) ICRD = IINTCR
      NVEC = ICRD - NPROJ
      IF (IPRINT .GT. IPRDBG) CALL TITLER('Output from FNDMOD','*',103)
      IF (ITRNMR .EQ. 0) THEN
         CALL DZERO(VECMOD,MXCOOR)
         IF (INTERN) THEN
            CALL GQ2GX(MXRCRD,EVEC(1,NSPMOD),VECMOD,WILBMT)
         ELSE
            CALL DCOPY(ICRD,EVEC(1,NSPMOD),1,VECMOD,1)
         END IF
         CALL NRMLVX(NCART,VECMOD)
         IMODE = NSPMOD
      ELSE
         CALL DZERO(TMPVC2,NVEC)
 5       CONTINUE
         REMAIN = .FALSE.
         OVRLAP = D0
         IOVRL  = 0
         DO 10 I = 1, NVEC
            IF (TMPVC2(I) .LT. 1.0D0) THEN
               REMAIN = .TRUE.
               IF (INTERN) THEN
                  CALL GQ2GX(MXRCRD,EVEC(1,I),TMPVEC,WILBMT)
               ELSE
                  CALL DCOPY(ICRD,EVEC(1,I),1,TMPVEC,1)
               END IF
               CALL NRMLVX(NCART,TMPVEC)
               OVR = DDOT(NCART,VECMOD,1,TMPVEC,1)
               IF (ABS(OVR) .GT. OVRLAP) THEN
                  IOVRL = I
                  OVRLAP = ABS(OVR)
               END IF
            END IF
 10      CONTINUE
C
C     If the mode with the largest overlap corresponds to a mode
C     with a gradient component equal to zero, we have to find
C     another mode.
C
         IF ((ABS(GRDDIA(IOVRL)) .LT. 1.0D-10) .AND. REMAIN) THEN
            TMPVC2(IOVRL) = 2.0D0
            GOTO 5
         END IF
         IF (.NOT. REMAIN) IOVRL = IMODE
         CALL DZERO(VECMOD,MXCOOR)
         IF (INTERN) THEN
            CALL GQ2GX(MXRCRD,EVEC(1,IOVRL),VECMOD,WILBMT)
         ELSE
            CALL DCOPY(ICRD,EVEC(1,IOVRL),1,VECMOD,1)
         END IF
         CALL NRMLVX(NCART,VECMOD)
         IMODE = IOVRL
      END IF
      IF (IPRINT .GE. IPRDBG) THEN
         CALL HEADER('Image mode eigenvector',1)
         CALL OUTPUT(EVEC,1,ICRD,IMODE,IMODE,MXRCRD,MXRCRD,1,LUPRI)
      END IF
      RETURN
      END

C  /* Deck numcrd */
      SUBROUTINE NUMCRD(ICART, IINT)
C
C     Simply returns the number of Cartesian (total) coordinates and
C     the number of redundant internal coordinates.
C
#include "implicit.h"
#include "mxcent.h"
#include "optinf.h"
      ICART = ICRTCR
      IINT  = IINTCR
      RETURN
      END
