!
!...   Copyright (c) 2001 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program, Release 1.2
!...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
!...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
!...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
!...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
!...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
!...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
!...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
!...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
!...   and O. Vahtras"
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
!...   For questions concerning this copyright write to:
!...      dalton-admin@kjemi.uio.no
!...
!...   For information on how to get a licence see:
!...      http://www.kjemi.uio.no/software/dalton/dalton.html
!
!
! File: abaop2.F
!
! 971020-vebjornb: Just an extension of abaopt.F which is becoming awkwardly large.
!
!  /* Deck linsrc */
      SUBROUTINE LS_LINSRC(ICRD,MCRD,GRAD,GRADOL,STEP,STEPOL, &
     &     TMPVEC,TMPVC2,ACTIVE,EMOD,lupri,optinfo)
use precision
use optimization_input
use ls_util
use files
use lstiming
!
!     This routine calculates the step that should be taken to obtain
!     the next geometry.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Type(opt_setting) :: optinfo
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      PARAMETER (D0 = 0.0E0_realk)
      LOGICAL ACTIVE
      Real(realk) GRAD(MCRD), GRADOL(MCRD), STEP(MCRD), STEPOL(MCRD)
      Real(realk) TMPVEC(MCRD), TMPVC2(MCRD)
      ACTIVE = .TRUE.
      call ls_dzero(TMPVEC,MCRD)
      call ls_dzero(STEP,MCRD)
      STNM = SQRT(DDOT(ICRD,STEPOL,1,STEPOL,1))
      DO 10 I = 1, ICRD
         TMPVEC(I) = -STEPOL(I)/STNM
 10   CONTINUE
      GRAD0 = DDOT(ICRD,GRADOL,1,TMPVEC,1)
      GRAD1 = DDOT(ICRD,GRAD,1,TMPVEC,1)
      IF ((ABS(GRAD0) .LT. 1.0E-6_realk) .OR. (ABS(GRAD1) .LT. 1.0E-6_realk)) THEN
         IF (optinfo%IPrint .GE. IPRDBG) THEN
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) 'Too close to minimum, line search skipped.'
         END IF
         ACTIVE = .FALSE.
         RETURN
      END IF
!
!     The energy and gradient from the last and current point, are
!     fitted to a quartic polynomial.
!
      CA = ABS(5.0E0_realk*GRAD0 - 2.0E0_realk*GRAD1 + 3.0E0_realk*optinfo%energyOld - 3.0E0_realk*optinfo%energy)
      CB = GRAD1 + GRAD0 + 2.0E0_realk*optinfo%energyOld - 2.0E0_realk*optinfo%energy - 2.0E0_realk*CA
      CC = optinfo%energy - optinfo%energyOld - GRAD0 - CA - CB
      CD = GRAD0
      CE = optinfo%energyOld
!
!     The line search methdod is only used if the minimum lies
!     between the two points.
!
      IF ((GRAD0 .LT. D0) .AND. (GRAD1 .GT. D0)) THEN
         THRG = 1.0E-5_realk*MIN(ABS(GRAD0),ABS(GRAD1))
         CRDA = D0
         CRDB = 1.0E0_realk
         GRDA = GRAD0
         GRDB = GRAD1
         ISAFE = 0
 15      CONTINUE
         ISAFE = ISAFE + 1
         IF (ISAFE .GE. 200) THEN
            ACTIVE = .FALSE.
            IF (optinfo%IPrint .GE. IPRDBG) THEN
               WRITE(LUPRI,*)
               WRITE(LUPRI,*) 'Line search failed, ignoring.'
            END IF
            RETURN
         END IF
         CRDC = CRDA + (CRDB-CRDA)* &
     &        MAX(0.1E0_realk,MIN(0.9E0_realk,ABS(GRDA/GRDB)*0.5E0_realk))
         GRDC = ((4.0E0_realk*CA*CRDC + 3.0E0_realk*CB)*CRDC + 2.0E0_realk*CC)*CRDC + CD
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
!
!     If the line search ends up almost back at the previous point, we do
!     not trust it and simply discard it.
!
         IF (CRDC .LT. 0.15E0_realk) CRDC = 1E0_realk
!
         DO 20 I = 1, ICRD
            STEP(I) = STEPOL(I)*(1.0E0_realk-CRDC)
            GRAD(I) = GRADOL(I) + (GRAD(I)-GRADOL(I))*CRDC
 20      CONTINUE
         EMOD = (((CA*CRDC + CB)*CRDC + CC)*CRDC + CD)*CRDC + CE
         IF (optinfo%IPrint .GE. IPRDBG) THEN
            call lsheader(lupri,'Interpolated step')
            call output(STEP,1,1,1,ICRD,1,MCRD,1,LUPRI)
            call lsheader(lupri,'Interpolated gradient')
            call output(GRAD,1,1,1,ICRD,1,MCRD,1,LUPRI)
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) 'Interpolated energy: ',EMOD
         END IF
      ELSE
         IF (optinfo%IPrint .GE. IPRDBG) THEN
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) 'Line search skipped.'
         END IF
         ACTIVE = .FALSE.
         RETURN
      END IF
      RETURN
      END
      
!  /* Deck upgdst */
      SUBROUTINE LS_UPGDST(ICRD,MXRCRD,GRDARR,STPARR,GRAD,STEP,lupri,optinfo)
use precision
use ls_util
use optimization_input
use files
use lstiming
!
!     Updates arrays containing steps and gradients from
!     previous iterations.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Type(opt_setting) :: optinfo
      Real(realk) GRDARR(MXRCRD,25), STPARR(MXRCRD,25)
      Real(realk) GRAD(*), STEP(*)
      IF (optinfo%Restart) THEN
         call ls_dzero(GRDARR,25*MXRCRD)
         call ls_dzero(STPARR,25*MXRCRD)
         optinfo%KEPTIT = 0
         optinfo%Restart = .FALSE.
      ELSE
         DO 10 I = MIN(24,optinfo%KEPTIT),1,-1
            DO 20 J = 1, ICRD
               GRDARR(J,I+1) = GRDARR(J,I)
               STPARR(J,I+1) = STPARR(J,I) - STEP(J)
 20         CONTINUE
 10      CONTINUE
         DO 30 I = 1, ICRD
            GRDARR(I,1) = GRAD(I)
            STPARR(I,1) = - STEP(I)
 30      CONTINUE
         optinfo%KEPTIT = MIN(25,optinfo%KEPTIT+1)
      END IF
      RETURN
      END

!  /* Deck gdistp */
      SUBROUTINE LS_GDISTP(NCRD,ICRD,MXRCRD,MX2CRD,STEP,GRAD,HESS,HESINV, &
     &     TMPMAT,TMPMT2,TMPMT3,TMPMT4,GRDARR,STPARR,lupri,optinfo)
use precision
use ls_util
use files
use lstiming
use optimization_input
!
!     Updates arrays containing steps and gradients from
!     previous iterations.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Type(opt_setting) :: optinfo
      LOGICAL STPSCL
      Real(realk) STEP(NCRD), GRAD(NCRD)
      Real(realk) HESS(NCRD,NCRD), HESINV(ICRD,ICRD)
      Real(realk) TMPMAT(MX2CRD,MX2CRD),TMPMT2(MX2CRD*MX2CRD)
      Real(realk) TMPMT3(MX2CRD*MX2CRD),TMPMT4(MX2CRD*MX2CRD)
      Real(realk) GRDARR(MXRCRD,25), STPARR(MXRCRD,25)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      PARAMETER (D0 = 0.0E0_realk, D1 = 1.0E0_realk)
      call ls_dzero(TMPMAT,MX2CRD*MX2CRD)
      call ls_dzero(TMPMT2,MX2CRD*MX2CRD)
      call ls_dzero(TMPMT3,MX2CRD*MX2CRD)
      call ls_dzero(TMPMT4,MX2CRD*MX2CRD)
!
      STPSCL = .TRUE.
!
!     First we have to construct the DIIS matrix from the
!     scalar products of the gradients (we use the gradients as
!     error vectors).
!
      IDIM = optinfo%KEPTIT + 2
      DO 10 I = 1, optinfo%KEPTIT
         DO 15 J = I, optinfo%KEPTIT
            TMPMAT(I+1,J+1) = DDOT(ICRD,GRDARR(1,I),1,GRDARR(1,J),1)
            TMPMAT(J+1,I+1) = TMPMAT(I+1,J+1)
 15      CONTINUE
 10   CONTINUE
      TMPMAT(1,1) = DDOT(ICRD,GRAD(1),1,GRAD(1),1)
      DO 17 I = 1, optinfo%KEPTIT
         TMPMAT(1,I+1) = DDOT(ICRD,GRAD(1),1,GRDARR(1,I),1)
         TMPMAT(I+1,1) = TMPMAT(1,I+1)
         TMPMAT(I+1,IDIM) = D1
         TMPMAT(IDIM,I+1) = D1
 17   CONTINUE
      TMPMAT(IDIM,IDIM) = D0
      TMPMAT(1,IDIM) = D1
      TMPMAT(IDIM,1) = D1
      IF (optinfo%IPrint .GE. IPRMAX) THEN
         call lsheader(lupri,'DIIS matrix')
         call output(TMPMAT,1,IDIM,1,IDIM, &
     &        MX2CRD,MX2CRD,1,LUPRI)
      END IF
 25   CONTINUE
      call ls_dzero(TMPMT2,IDIM*IDIM)
      call ls_dzero(TMPMT3,IDIM*IDIM)
      DO 30 J = 1, IDIM
         DO 32 I = 1, IDIM
            TMPMT2(I+(J-1)*IDIM) = TMPMAT(I,J)
 32      CONTINUE
 30   CONTINUE
     call ls_DSITSP(IDIM,TMPMT2,TMPMT3)
     call ls_DUNIT(TMPMT2,IDIM)
     call ls_JACO(TMPMT3,TMPMT2,IDIM,IDIM,IDIM,TMPMT4(1), &
    &     TMPMT4(MX2CRD+1))
      call ls_dzero(TMPMT4,MX2CRD*MX2CRD)
      IZER = 0
      EMAX = D0
      DO 34 J = 1, IDIM-1
         TMPMT4(J) = TMPMT3(J*(J+1)/2)
         IF (ABS(TMPMT4(J)) .GT. EMAX) EMAX = ABS(TMPMT4(J))
 34   CONTINUE
      TMPMT4(IDIM) = TMPMT3(IDIM*(IDIM+1)/2)
      DO 35 J = 1, IDIM
         IF (ABS(TMPMT4(J)) .LE. EMAX*1.0E-2_realk) THEN
            IZER = IZER + 1
            LSTZER = J
         END IF
 35   CONTINUE
      IF (optinfo%IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'Initial eigenvalues of DIIS matrix')
         call output(TMPMT4,1,1,1,IDIM,1,IDIM,1,LUPRI)
      END IF
      IF ((IZER .GT. 0) .AND. (IDIM .GT. 2)) THEN
!
!     We try to remove only the components causing linearity
!     (one at a time)
!
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
         optinfo%NRemove = MAX(MAX1,MAX2)
         DO 370 I = optinfo%NRemove, IDIM-1
            DO 371 J = 1, IDIM
               TMPMAT(I,J) = TMPMAT(I+1,J)
 371        CONTINUE
 370     CONTINUE
         DO 380 I = optinfo%NRemove, IDIM-1
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
      IF (optinfo%IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'Eigenvalues of DIIS matrix')
         call output(TMPMT4,1,1,1,IDIM,1,IDIM,1,LUPRI)
         call lsheader(lupri,'Eigenvectors of DIIS matrix')
         call output(TMPMAT,1,IDIM,1,IDIM,MX2CRD,MX2CRD,1,LUPRI)
      END IF
!
!     We find the inverse of the DIIS matrix
!
      call ls_dzero(TMPMT2,MX2CRD*MX2CRD)
      call ls_dzero(TMPMT3,MX2CRD*MX2CRD)
      DO 45 I = 1, IDIM
         IF (ABS(TMPMT4(I)) .LE. 1.0E-6_realk) THEN
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
      call ls_dzero(TMPMAT,MX2CRD*MX2CRD)
      DO 60 J = 1, IDIM
         DO 62 I = 1, IDIM
            DO 64 K = 1, IDIM
               TMPMAT(I,J) = TMPMAT(I,J) +  &
     &              TMPMT3(I+(K-1)*IDIM)*TMPMT2(K+(J-1)*IDIM)
 64         CONTINUE
 62      CONTINUE
 60   CONTINUE
      call ls_dzero(TMPMT4,MX2CRD*MX2CRD)
      FAC = D0
      DO 70 I = 1, IDIM-1
         TMPMT4(I) = TMPMAT(I,IDIM)
         FAC = FAC + TMPMT4(I)
 70   CONTINUE
      IF (optinfo%IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'Inverse of DIIS matrix')
         call output(TMPMAT,1,IDIM,1,IDIM,MX2CRD,MX2CRD,1,LUPRI)
         call lsheader(lupri,'DIIS coefficients')
         call output(TMPMT4,1,1,1,IDIM-1,1,MX2CRD,1,LUPRI)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Sum of coefficients: ',FAC
      END IF
      call ls_dzero(TMPMT3,MX2CRD)
      call ls_dzero(STEP,NCRD)
      DO 80 I = 1, ICRD
         TMPMT3(I) = TMPMT4(1)*GRAD(I)
 80   CONTINUE
      DO 82 I = 1, IDIM-2
         DO 84 J = 1, ICRD
            STEP(J) = STEP(J) + TMPMT4(I+1)*STPARR(J,I)
 84      CONTINUE
 82   CONTINUE
!
!     If the step is too large, we simply restrict each element
!     to be below the trust radius.
!     
      IF (.NOT. STPSCL) THEN
         DO 85 I = 1, ICRD
            IF (ABS(STEP(I)) .GT. optinfo%TrustRad) &
     &           STEP(I) = SIGN(optinfo%TrustRad,STEP(I))
 85      CONTINUE
      ELSE
!
!     Alternatively we restrivt the step norm to be equal or less
!     than the trust radius.
!
         STPNRM = SQRT(DDOT(ICRD,STEP,1,STEP,1))
         IF (STPNRM .GT. optinfo%TrustRad) THEN
            FAC = optinfo%TrustRad/STPNRM
            IF (optinfo%IPrint .GE. IPRDBG) THEN
               WRITE(LUPRI,*)
               WRITE(LUPRI,*) &
     &              'DIIS-step too long, scaling by factor:', FAC
               WRITE(LUPRI,*)
            END IF
            DO 86 I = 1, ICRD
               STEP(I) = STEP(I)*FAC
 86         CONTINUE
         END IF
      END IF
      IF (optinfo%IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'Interpolated step')
         call output(STEP,1,1,1,ICRD,1,NCRD,1,LUPRI)
         call lsheader(lupri,'Interpolated gradient')
         call output(TMPMT3,1,1,1,ICRD,1,MX2CRD,1,LUPRI)
         call lsheader(lupri,'Inverse of Hessian')
         call output(HESINV,1,ICRD,1,ICRD,ICRD,ICRD,1,LUPRI)
      END IF
      call ls_dzero(TMPMT2,ICRD)
      DO 100 I = 1, ICRD
         DO 102 J = 1, ICRD
            TMPMT2(I) = TMPMT2(I) + HESINV(I,J)*TMPMT3(J)
 102     CONTINUE
         STEP(I) = STEP(I) - TMPMT2(I)
 100  CONTINUE
      IF (optinfo%IPrint .GE. IPRMAX) THEN
         call lsheader(lupri,'Relaxation step')
         call output(TMPMT2,1,1,1,ICRD,1,NCRD,1,LUPRI)
         call lsheader(lupri,'Total DIIS step')
         call output(STEP,1,1,1,ICRD,1,NCRD,1,LUPRI)
      END IF
!
!     If the step is too large, we simply restrict each element
!     to be below the trust radius.
!     
      IF (.NOT. STPSCL) THEN
         DO 185 I = 1, ICRD
            IF (ABS(STEP(I)) .GT. optinfo%TrustRad) &
     &           STEP(I) = SIGN(optinfo%TrustRad,STEP(I))
 185     CONTINUE
      ELSE
!     
!     Alternatively we restrivt the step norm to be equal or less
!     than the trust radius.
!     
         STPNRM = SQRT(DDOT(ICRD,STEP,1,STEP,1))
         IF (STPNRM .GT. optinfo%TrustRad) THEN
            FAC = optinfo%TrustRad/STPNRM
            IF (optinfo%IPrint .GE. IPRDBG) THEN
               WRITE(LUPRI,*)
               WRITE(LUPRI,*) &
     &              'Relaxation step too long, scaling by factor:', FAC
               WRITE(LUPRI,*)
            END IF
            DO 186 I = 1, ICRD
               STEP(I) = STEP(I)*FAC
 186        CONTINUE
         END IF
      END IF
      ICNT = 1
      FAC = 1.0E0_realk
 200  CONTINUE
      IF (ICNT .LE. 10) THEN
         call ls_dzero(TMPMT2,MX2CRD)
         SNDTRM = 0.0E0_realk
         DO 210 I = 1, ICRD
            DO 215 J = 1, ICRD
               TMPMT2(I) = TMPMT2(I) + HESS(I,J)*STEP(J)
 215        CONTINUE
            SNDTRM = SNDTRM + TMPMT2(I)*STEP(I)
 210     CONTINUE
         optinfo%predictedChange = DDOT(ICRD,GRAD,1,STEP,1) &
     &        + 0.5E0_realk*SNDTRM
         IF (optinfo%predictedChange .GT. 0.0E0_realk) THEN
            DO 220 I = 1, ICRD
               STEP(I) = STEP(I)*FAC
 220        CONTINUE
            FAC = 0.5E0_realk*FAC
            ICNT = ICNT + 1
            GOTO 200
         END IF
      ELSE
         IF (optinfo%IPrint .GE. IPRMIN) THEN
            IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
               call lsheader(lupri,'DIIS-step in internal coordinates')
            ELSE
               call lsheader(lupri,'DIIS-step')
            END IF
            call output(STEP,1,1,1,ICRD,1,NCRD,1,LUPRI)
         END IF
         DO 230 I = 1, ICRD
            STEP(I) = -GRAD(I)
 230     CONTINUE
         call ls_dzero(TMPMT2,MX2CRD)
         SNDTRM = 0.0E0_realk
         DO 240 I = 1, ICRD
            DO 245 J = 1, ICRD
               TMPMT2(I) = TMPMT2(I) + HESS(I,J)*STEP(J)
 245        CONTINUE
            SNDTRM = SNDTRM + TMPMT2(I)*STEP(I)
 240     CONTINUE
         optinfo%predictedChange = DDOT(ICRD,GRAD,1,STEP,1) &
     &        + 0.5E0_realk*SNDTRM
      END IF
      STPNRM = SQRT(DDOT(ICRD,STEP,1,STEP,1))
      RETURN
      END
    SUBROUTINE LS_MAKIMG(ICRD,IPRJ,MCRD,IMODE,RESTOR,lupri,optinfo)
!  /* Deck makimg */
use precision
use ls_util
use optimization_input
use files
use lstiming
!
!     Makes image function and sorts coordinates.
!     It is also used to restore and resort coordinates.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Type(opt_setting) :: optinfo
      SAVE ENEG, GNEG, ISORT
      LOGICAL RESTOR
      real(realk),PARAMETER ::DUMMY = 1.0E20_realk
      integer,parameter :: IDUMMY = - 9999999
      IF (.NOT. RESTOR) THEN
         ENEG = -optinfo%EVAL(IMODE)
         GNEG = -optinfo%GRDDIA(IMODE)
         optinfo%EVAL(IMODE) = ENEG
         optinfo%GRDDIA(IMODE) = DUMMY
         call ls_ORDER(optinfo%GRDDIA,optinfo%EVAL,IPRJ,1)
         DO 10 I = 1, IPRJ
            IF ((optinfo%EVAL(I).EQ.ENEG).AND.(optinfo%GRDDIA(I).EQ.DUMMY)) ISORT = I
 10      CONTINUE
         optinfo%GRDDIA(ISORT) = GNEG
      ELSE
         SNEG = optinfo%STPDIA(ISORT)
         IF (ISORT .LT. IMODE) THEN
            DO 20 I = ISORT, IMODE-1
               optinfo%EVAL(I)   = optinfo%EVAL(I+1)
               optinfo%GRDDIA(I) = optinfo%GRDDIA(I+1)
               optinfo%STPDIA(I) = optinfo%STPDIA(I+1)
 20         CONTINUE
            optinfo%EVAL(IMODE)   = -ENEG
            optinfo%GRDDIA(IMODE) = -GNEG
            optinfo%STPDIA(IMODE) = SNEG
         ELSE IF (IMODE .LT. ISORT) THEN
            DO 30 I = ISORT, IMODE+1,-1
               optinfo%EVAL(I)   = optinfo%EVAL(I-1)
               optinfo%GRDINT(I) = optinfo%GRDINT(I-1)
               optinfo%STPDIA(I) = optinfo%STPDIA(I-1)
 30         CONTINUE
            optinfo%EVAL(IMODE)   = -ENEG
            optinfo%GRDINT(IMODE) = -GNEG
            optinfo%STPDIA(IMODE) = SNEG
         ELSE
            optinfo%EVAL(IMODE)   = -ENEG
            optinfo%GRDINT(IMODE) = -GNEG
         END IF
      END IF
      RETURN
      END

!  /* Deck fndmod */
      SUBROUTINE LS_FNDMOD(INTERN,MXRCRD,EVEC,WILBMT,VECMOD,TMPVEC, &
     &     TMPVC2,IMODE,lupri,optinfo)
use precision
use ls_util
use optimization_input
use files
use lstiming
!
!     Makes sure we follow the same eigenvectormode throughout
!     a saddle point optimization. The mode is selected as the mode
!     (in Cartesian coordinates) with the largest overlap with
!     the previous mode. In redundant coordinates we have to check
!     that the selected mode corresponds to non-zero diagonal element.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Type(opt_setting) :: optinfo
      Real(realk) EVEC(MXRCRD,MXRCRD), WILBMT(MXRCRD,MXCOOR)
      Real(realk) VECMOD(MXCOOR), TMPVEC(MXCOOR), TMPVC2(MXCOOR)
      LOGICAL INTERN, REMAIN
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      PARAMETER (D0 = 0.0E0_realk)
      IF (INTERN)THEN
         ICRD = optinfo%NIntCoord
      ELSE
         call lsquit('NCART: no value assigned to this variable',-1)
         ICRD = 0
         !      ICRD = NCART
      ENDIF
      NVEC = ICRD - optinfo%nProjected
      IF (optinfo%IPrint .GT. IPRDBG) call lsheader(lupri,'Output from FNDMOD')
      IF (optinfo%ItrNmr .EQ. 0) THEN
         call ls_dzero(VECMOD,MXCOOR)
         IF (INTERN) THEN
            call ls_GQ2GX(MXRCRD,EVEC(1,optinfo%NSPMod),VECMOD,WILBMT,optinfo)
         ELSE
            call dcopy(ICRD,EVEC(1,optinfo%NSPMod),1,VECMOD,1)
         END IF
         call ls_NRMLVX(NCART,VECMOD)
         IMODE = optinfo%NSPMod
      ELSE
         call ls_dzero(TMPVC2,NVEC)
 5       CONTINUE
         REMAIN = .FALSE.
         OVRLAP = D0
         IOVRL  = 0
         DO 10 I = 1, NVEC
            IF (TMPVC2(I) .LT. 1.0E0_realk) THEN
               REMAIN = .TRUE.
               IF (INTERN) THEN
                  call ls_GQ2GX(MXRCRD,EVEC(1,I),TMPVEC,WILBMT,optinfo)
               ELSE
                  call dcopy(ICRD,EVEC(1,I),1,TMPVEC,1)
               END IF
               call ls_NRMLVX(NCART,TMPVEC)
               OVR = DDOT(NCART,VECMOD,1,TMPVEC,1)
               IF (ABS(OVR) .GT. OVRLAP) THEN
                  IOVRL = I
                  OVRLAP = ABS(OVR)
               END IF
            END IF
 10      CONTINUE
!
!     If the mode with the largest overlap corresponds to a mode
!     with a gradient component equal to zero, we have to find
!     another mode.
!
         IF ((ABS(optinfo%GRDINT(IOVRL)) .LT. 1.0E-10_realk) .AND. REMAIN) THEN
            TMPVC2(IOVRL) = 2.0E0_realk
            GOTO 5
         END IF
         IF (.NOT. REMAIN) IOVRL = IMODE
         call ls_dzero(VECMOD,MXCOOR)
         IF (INTERN) THEN
            call ls_GQ2GX(MXRCRD,EVEC(1,IOVRL),VECMOD,WILBMT,optinfo)
         ELSE
            call dcopy(ICRD,EVEC(1,IOVRL),1,VECMOD,1)
         END IF
         call ls_NRMLVX(NCART,VECMOD)
         IMODE = IOVRL
      END IF
      IF (optinfo%IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'Image mode eigenvector')
         call output(EVEC,1,ICRD,IMODE,IMODE,MXRCRD,MXRCRD,1,LUPRI)
      END IF
      RETURN
      END

!  /* Deck numcrd */
      SUBROUTINE LS_NUMCRD(ICART, IINT,optinfo)
use ls_util
use precision
use optimization_input
use files
use lstiming
!
!     Simply returns the number of Cartesian (total) coordinates and
!     the number of redundant internal coordinates.
!
Implicit Real(realk) (A-H,O-Z)
      Type(opt_setting) :: optinfo
      ICART = optinfo%ICartCoord
      IINT  = optinfo%NIntCoord
      RETURN
      END

!  /* Deck ls_order */
      SUBROUTINE LS_ORDER(EVEC,EVAL,N,NEVEC)
!
!     Copied from gphjj.F
!
! Purpose: order the N values in optinfo%EVAL and their associated vectors
!          in EVEC so optinfo%EVAL(i+1) .ge. optinfo%EVAL(i)
!
! Revisions:
!   29-Jul-1992 hjaaj (only dswap if nevec .gt. 0)
!    2-Nov-1984 hjaaj (new parameter NEVEC, EVEC(1:NEVEC,1:N))
!   27-Oct-1984 hjaaj (reduced number of swaps)
!
use precision
Implicit Real(realk) (A-H,O-Z)
      Real(realk) EVEC(*),EVAL(*)
      IF (N.LE. 1) RETURN
      IN = 1
      DO 10 I=1,N-1
        EMIN = EVAL(I)
        IMIN = I
        DO 20 J=I+1,N
          IF (EVAL(J) .LT. EMIN) THEN
            EMIN = EVAL(J)
            IMIN = J
          ENDIF
   20   CONTINUE
        IF (IMIN.NE.I) THEN
          EVAL(IMIN)=EVAL(I)
          EVAL(I)=EMIN
          IF (NEVEC .GT. 0) THEN
            CALL DSWAP(NEVEC,EVEC(IN),1,EVEC((IMIN-1)*NEVEC+1),1)
          ENDIF
        ENDIF
        IN = IN + NEVEC
   10 CONTINUE
      RETURN
      END
!
!  /* Deck ortvec */
      SUBROUTINE LS_ORTVEC(NOLD,NNEW,NVDIM,THRLDP,VEC,lupri)
use precision
use ls_util
!
! 15-Mar-1985 hjaaj
! l.r. 4-May-1994 hjaaj (only elim. new vector if norm < THRLDP)
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Real(realk) VEC(NVDIM,6)
      LOGICAL   LINDEP(6)
!
      Real(realk), PARAMETER :: THRRND=1E-2_realk, THRTT=1E-4_realk, D1=1.0E0_realk
!
      LINDEP = .FALSE.
      IF (NNEW.LE. 0) CALL LSQUIT('NNEW less than 0',lupri)
      TLINDP = SQRT(THRLDP)
!
!     Normalize NNEW new vectors VEC(*,NOLD+1) TO VEC(*,NOLD+NNEW)
!
      IVEC = NOLD
      DO 200 INEW = 1,NNEW
         IVEC = IVEC + 1
         TT   = dnorm2_ls(NVDIM,VEC(:,IVEC),1)
         IF (TT.LE.THRLDP) THEN
            LINDEP(INEW) = .TRUE.
            WRITE (LUPRI,8100) INEW,TT
         ELSE
            LINDEP(INEW) = .FALSE.
            IF (TT.LT.THRTT) THEN
               CALL DSCAL (NVDIM,(D1/TT),VEC(1,IVEC),1)
               TT = dnorm2_ls(NVDIM,VEC(:,IVEC),1)
            END IF
            CALL DSCAL (NVDIM,(D1/TT),VEC(1,IVEC),1)
         END IF
  200 CONTINUE
!
!
!
      IROUND = 0
      ITURN  = 0
 1500 ITURN  = ITURN + 1
!
!     Orthogonalize new vectors against previous vectors
!
      DO 2000 K=1,NOLD
         DO 1900 J = NOLD+1,NOLD+NNEW
            TT = - DDOT(NVDIM,VEC(1,K),1,VEC(1,J),1)
            CALL DAXPY(NVDIM,TT,VEC(1,K),1,VEC(1,J),1)
 1900    CONTINUE
 2000 CONTINUE
!
!     Orthogonalize new vectors against each other
!     and normalization.
!
      DO 2400 INEW = 1,NNEW
         IF (.NOT.LINDEP(INEW)) THEN
!           ... orthogonalize using prev. vectors are normalized
            IVEC = NOLD + INEW
            DO 2300 JNEW = 1,(INEW-1)
               IF (.NOT.LINDEP(JNEW)) THEN
                  JVEC = NOLD + JNEW
                  TT = - DDOT(NVDIM,VEC(1,JVEC),1,VEC(1,IVEC),1)
                  CALL DAXPY(NVDIM,TT,VEC(1,JVEC),1,VEC(1,IVEC),1)
               END IF
 2300       CONTINUE
            TT = dnorm2_ls(NVDIM,VEC(:,IVEC),1)
            IF (TT .LE. THRLDP) THEN
               LINDEP(INEW) = .TRUE.
               WRITE (LUPRI,8100) INEW,TT
            ELSE
               IF (TT .LT. THRRND) IROUND = IROUND+1
               IF (TT .LT. THRTT) THEN
                  CALL DSCAL(NVDIM,(D1/TT),VEC(1,IVEC),1)
                  TT = dnorm2_ls(NVDIM,VEC(:,IVEC),1)
               END IF
               CALL DSCAL(NVDIM,(D1/TT),VEC(1,IVEC),1)
            END IF
         END IF
 2400 CONTINUE
!
!
      IF (IROUND.GT. 0 .AND. ITURN.EQ. 1) GO TO 1500
!
!
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
!
!
!
 8100 FORMAT(/' ORTVEC: New vector no.',I3, &
      &       ' is removed because of linear dependence; &
      &      norm of vector after Gram-Schmidt''s orthogonalization',&
      &      1P,D9.2)
 8200 FORMAT(/' (ORTVEC) NNEW reduced from',I3,' to',I3)
!
! *** End of subroutine ORTVEC
!
      RETURN
      END
!  /* Deck around */
      SUBROUTINE LS_AROUND(HEAD,lupri)
use precision
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      CHARACTER HEAD*(*)
      LNG = LEN(HEAD) + 2
      IND = (72 - LNG)/2 + 1
      WRITE (LUPRI, '(//,80A)') (' ',I=1,IND), '+', ('-',I=1,LNG), '+'
      WRITE (LUPRI, '(80A)')    (' ',I=1,IND), '! ', HEAD, ' !'
      WRITE (LUPRI, '(80A)')    (' ',I=1,IND), '+', ('-',I=1,LNG), '+'
!     WRITE (LUPRI, '(//,80A)') (' ',I=1,IND), '.', ('-',I=1,LNG), '.'
!     WRITE (LUPRI, '(80A)')    (' ',I=1,IND), '| ', HEAD, ' |'
!     WRITE (LUPRI, '(80A)')    (' ',I=1,IND), '`', ('-',I=1,LNG), ''''
      WRITE (LUPRI, '()')
      RETURN
      END
!  /* Deck dsitsp */
      SUBROUTINE LS_DSITSP(N,ASI,ASP)
!
!  8-Feb-1987 Hans Joergen Aa. Jensen
!
! Purpose: Transform from SI format to SP format, that is:
!          copy upper triangle of symmetric matrix ASI
!          to symmetric, packed matrix ASP.
!
use precision
Implicit Real(realk) (A-H,O-Z)
      Real(realk) ASI(N,*), ASP(*)
!
      DO 200 J = 1,N
         JOFF = (J*J-J)/2
         DO 100 I = 1,J
            ASP(JOFF+I) = ASI(I,J)
  100    CONTINUE
  200 CONTINUE
!
      RETURN
      END
!  /* Deck jaco */
      SUBROUTINE LS_JACO (F,V,NB,NMAX,NROWV,BIG,JBIG)
!
!     F is symmetric packed matrix of Real(realk) NB.
!     The first block of size NMAX will be diagonalized.
!     V is for eigenvectors, only V(NROWV,NMAX) will be referenced.
!     On entry it must correspond the basis vectors corresponding to the
!       F matrix on entry, e.g. unit matrix or AO coefficients for each MO.
!
use ls_util
use precision
Implicit Real(realk) (A-H,O-Z)
      Real(realk) F(*),V(*)
      Real(realk) BIG(*) ,JBIG(*)
      PARAMETER (D0 = 0.0E0_realk, D1 = 1.0E0_realk, ROOT2 = 0.707106781186548E0_realk)
      PARAMETER(DP5 = 0.5E0_realk, D1P5 = 1.5E0_realk, D1P375 = 1.375E0_realk, &
              & D3P875 = 3.875E0_realk, DP25 = 0.25E0_realk)
      PARAMETER (THRZER = 1.0E-14_realk)
      DATA C1,C2,C3,C4,C5,C6/THRZER,THRZER,1E-20_realk,1E-14_realk,1E-9_realk,1E-5_realk/
!     DATA C1,C2,C3,C4,C5,C6/1E-12_realk,1E-12_realk,1E-20_realk,1E-14_realk,1E-9_realk,1E-5_realk/
      IF (NB.LE. 1 .OR. NMAX.LE. 0) RETURN
      IF (NB.EQ. 1) RETURN !900727-hjaaj
      DO 190 I=1,NB
         JBIGI=0
         J=MIN(I-1,NMAX)
         IF (J .GT. 0) THEN
            II = (I*I-I)/2
            ABIGI=D0
            DO 18 K=1,J
            IF (ABIGI .GE. ABS(F(II+K))) GO TO  18
               ABIGI=ABS(F(II+K))
               JBIGI=K
   18       CONTINUE
         END IF
         IF (JBIGI .GT. 0) THEN
            JBIG(I) = JBIGI
            BIG(I)  = F(II+JBIGI)
         ELSE
            JBIG(I) = 0
            BIG(I)  = D0
         END IF
  190 CONTINUE
! 921030-hjaaj: SD = C1 now
      NNB = (NB*NB+NB)/2
!     SD = 1.05E0_realk
!     DO 220 J = 1,NNB
!        SD = MAX(SD, ABS(F(J)) )
! 220 CONTINUE
!     SD=MAX(C1,C2*SD)
      SD=C1
!
      MXITJA = 50*NNB
      ITJACO = 0
  410 ITJACO = ITJACO + 1
      IF (ITJACO .GT. MXITJA) THEN
         CALL LSQUIT('ERROR: JACO did not converge ...',lupri)
      END IF
      T = D0
      DO 230 I=2,NB
      IF (T .GE.  ABS(BIG(I))) GO TO 230
         T = ABS(BIG(I))
         IB= I
  230 CONTINUE
      IF(T.LT.SD) GO TO 420
         IA =JBIG(IB)
         IAA=IA*(IA-1)/2
         IBB=IB*(IB-1)/2
         DIF=F(IAA+IA)-F(IBB+IB)
         IF( ABS(DIF) .GT. C3) GO TO 271
            SX=ROOT2
            CX=ROOT2
         GO TO 270
  271       T2X2 =BIG(IB)/DIF
            T2X25=T2X2*T2X2
         IF(T2X25 .GT. C4) GO TO 240
            CX=1E00_realk
            SX=T2X2
         GO TO 270
  240    IF(T2X25 .GT. C5) GO TO 250
            SX=T2X2*(D1 - D1P5*T2X25)
            CX=D1 - DP5*T2X25
         GO TO 270
  250    IF(T2X25 .GT. C6) GO TO 260
            CX=D1+T2X25*(T2X25*D1P375 - DP5 )
            SX= T2X2*(D1 + T2X25*(T2X25*D3P875 - D1P5))
         GO TO 270
  260       T = DP25  / SQRT(DP25   + T2X25)
            CX= SQRT(DP5   + T)
            SX= SIGN( SQRT(DP5 - T),T2X2)
  270    CONTINUE
         DO 275 IR=1,IA
            T        = F(IAA+IR)*SX
            F(IAA+IR)= F(IAA+IR)*CX+F(IBB+IR)*SX
            F(IBB+IR)=-T           +F(IBB+IR)*CX
  275    CONTINUE
         IEAA=IAA+IA
         IEAB=IBB+IA
         TT  =F(IEAB)
         F(IEAB)=BIG(IB)
         IF (JBIG(IA) .EQ. 0) THEN
            IRST = IA   + 1
            IEAR = IEAA + IA
            IEBR = IEAB + 1
         ELSE
            IRST = IA
            IEAR = IEAA
            IEBR = IEAB
         END IF
         DO 390 IR = IRST,NB
            IF (IR .EQ. IA) THEN
               GO TO 360
!              ... we have checked above that JBIG(IA) .ne. 0
!              IF(JBIG(IR)) 360,380,360
            END IF
            T      = F(IEAR)*SX
            F(IEAR)= F(IEAR)*CX+F(IEBR)*SX
            F(IEBR)=-T         +F(IEBR)*CX
            T   =F(IEAR)
            IT  =IA
            IF(IR-IB) 340,310,320
  310          F(IEAA)=F(IEAA)*CX+F(IEAB)*SX
!              921030+hjaaj: zero f(ab) to avoid round-off errors
!              F(IEAB)=     TT*CX+F(IEBR)*SX
               F(IEAB)=     D0
               F(IEBR)=    -TT*SX+F(IEBR)*CX
            GO TO 360
  320       IF(ABS(T) .GE.  ABS(F(IEBR))) GO TO 340
               T   =F(IEBR)
               IT  =IB
  340       IF(ABS(T) .LT.  ABS(BIG(IR))) GO TO 350
               BIG(IR)  = T
               JBIG(IR) = IT
            GO TO 380
  350       IF(IA .NE. JBIG(IR) .AND. IB .NE. JBIG(IR))  GO TO 380
  360          K= IEAR - IA
               JBIGI = 0
               IR1=MIN (IR-1,NMAX)
               IF (IR1 .GT. 0) THEN
                  ABIGI = D0
                  DO 370 I=1,IR1
                  IF(ABIGI .GE. ABS(F(K+I)))  GO TO 370
                     ABIGI = ABS(F(K+I))
                     JBIGI =I
  370             CONTINUE
               END IF
               IF (JBIGI .GT. 0) THEN
                  JBIG(IR) = JBIGI
                  BIG(IR)  = F(K+JBIGI)
               ELSE
                  JBIG(IR) = 0
                  BIG(IR)  = D0
               END IF
  380       CONTINUE
               IEAR = IEAR + IR
               IF (IR .GE. IB) THEN
                  IEBR = IEBR + IR
               ELSE
                  IEBR = IEBR + 1
               END IF
  390       CONTINUE
         JAA=(IA-1)*NROWV
         JBB=(IB-1)*NROWV
         DO 400 I=1,NROWV
            T=V(JBB+I)*SX
            V(JBB+I)=-V(JAA+I)*SX + V(JBB+I)*CX
  400       V(JAA+I)= V(JAA+I)*CX + T
      GO TO 410
  420 CONTINUE
      RETURN 
      END
!  /* Deck wlkeig */
SUBROUTINE LS_WLKEIG(GRDCAR,HESCAR,EVAL,EVEC,GRDDIA,TMAT,THRLDP,&
     &    THRIND,CNDHES,INDHES,NCRTOT,N2CRT,NTMAT,& 
     &    PRJTRO,IPrint,lupri)
use ls_util
use precision
use memory_handling
use optimization_input
Implicit Real(realk) (A-H,O-Z)
   Integer :: lupri
   LOGICAL PRJTRO
   Real(realk) GRDCAR(NCRTOT), HESCAR(N2CRT), TMAT(NTMAT),EVEC(N2CRT), &
      &      EVAL(NCRTOT),CNDHES(0:7), GRDDIA(NCRTOT)
   Integer  IWORK(NCRTOT)
   Real(realk), pointer :: KPACK(:),KOVRL(:),KWORK(:)
   Logical,pointer :: KTRVEC(:)
   INTEGER :: INDHES
   Integer :: NCRTOT
   INTEGER :: IPRINT
!
      IF (IPrint .GT. 5) CALL LSHEADER(lupri,'Output from WLKEIG')
!
      CALL ls_dzero(GRDDIA,NCRTOT)
      NPR = 6
! Allocate memory fot ls_wlkei1
      Call mem_alloc(KPACK,NCRTOT*(NCRTOT+1)/2) 
      Call mem_alloc(KOVRL,NCRTOT)
      Call mem_alloc(KWORK,NCRTOT)
      Call mem_alloc(KTRVEC,NCRTOT)
!
      CALL LS_WLKEI1(GRDCAR,HESCAR,GRDDIA,EVAL,&
     &               EVEC,KOVRL,TMAT,KPACK,&
     &               KWORK,THRLDP,THRIND,IWORK,&
     &               KTRVEC,NCRTOT,NPR,INDEX,INDEXM,&
     &               PRJTRO,IPrint,lupri)
! Deallocate
      Call mem_dealloc(KPACK)
      Call mem_dealloc(KOVRL)
      Call mem_dealloc(KWORK)
      Call mem_dealloc(KTRVEC)
      INDHESP = INDEX
      IF (PRJTRO) THEN
          CALL LS_WLKCND(CNDNMB,EVAL,NCRTOT-NPR,lupri)
          CNDHES(0) = CNDNMB
      END IF
      RETURN
      END
!  /* Deck wlkei1 */
      SUBROUTINE LS_WLKEI1(GRDCAR,HESCAR,GRDDIA,EVAL,EVEC,OVERLP,TMAT,&
     &                  HESPCK,WORK,THRLDP,THRIND,IWORK,TRVEC, &
     &                  NCR,NPR,INDEX,INDEXM,PRJTRO,&
                        IPrint,lupri)
use precision
use ls_util
use optimization_input
use memory_handling
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      PARAMETER (D0 = 0.0E0_realk, THRZER = 1.0E-04_realk)
!
      LOGICAL PRJTRO, TRVEC(NCR)
      Real(realk) GRDCAR(NCR), HESCAR(NCR,NCR), GRDDIA(NCR), EVAL(NCR),&
      &         EVEC(NCR,NCR), HESPCK(NCR*(NCR + 1)/2), TMAT(NCR,NPR),&
      &         WORK(NCR)
      Integer   IWORK(NCR) 
      Real(realk)  OVERLP(NCR)
      Real(realk), pointer :: Work_mem(:) 
      INTEGER :: IPRINT,INFO
      INFO=0
!
      IF (IPrint .GT. 5) CALL LSHEADER(lupri,'Output from WLKEI1')
!
!     ***** Diagonalize Hessian *****
!
!      CALL LS_DSITSP(NCR,HESCAR,HESPCK)
!Write(*,*)'HESPCK',HESPCK
!      CALL LS_DUNIT(EVEC,NCR)
!         Write(*,*)'Eig',Evec(:,1)
!      CALL LS_JACO(HESPCK,EVEC,NCR,NCR,NCR,WORK,IWORK)
!      DO 100 I=1, NCR
!         EVAL(I) = HESPCK(I*(I+1)/2)
!         Write(*,*)'EVAL',EVAL(I)
! 100  CONTINUE
EVEC = HESCAR
Call mem_alloc(Work_mem,3*NCR)
Call dsyev('V','U',NCR,EVEC,NCR,EVAL,WORK_mem,3*NCR,INFO)
Call mem_dealloc(Work_mem)

!
      DO I=1,NCR
         GRDDIA(I) = DDOT(NCR,GRDCAR,1,EVEC(1:NCR,I),1)
      ENDDO
      IF (IPrint .GT. 5) THEN
         CALL LSHEADER(lupri,'Eigenvalues in WLKEI1')
         CALL OUTPUT(EVAL,1,1,1,NCR,1,NCR,1,LUPRI)
         CALL LSHEADER(lupri,'Eigenvectors in WLKEI1')
         CALL OUTPUT(EVEC,1,NCR,1,NCR,NCR,NCR,1,LUPRI)
         CALL LSHEADER(lupri,'Gradient (diagonal rep.) in WLKEI1')
         CALL OUTPUT(GRDDIA,1,1,1,NCR,1,NCR,1,LUPRI)
      END IF
!
!     ***** Sort in order of increasing trarot component *****
!
      IF (PRJTRO) THEN
         DO 300 ICR = 1, NCR
            OVLP = D0
            DO 310 ITR = 1, NPR
               OVLP = OVLP + ABS(DDOT(NCR,EVEC(1,ICR),1,TMAT(1,ITR),1))
  310       CONTINUE
            OVERLP(ICR) = OVLP
  300    CONTINUE
         CALL DCOPY(NCR,OVERLP,1,WORK,1)
         CALL ls_ORDER(GRDDIA,WORK,NCR,1)
         CALL DCOPY(NCR,OVERLP,1,WORK,1)
         CALL ls_ORDER(EVAL,  WORK,NCR,1)
         CALL DCOPY(NCR,OVERLP,1,WORK,1)
         CALL ls_ORDER(EVEC,  WORK,NCR,NCR)
      END IF
!
      IF (PRJTRO) THEN
         NVEC = NCR - NPR
      ELSE
         NVEC = NCR
      END IF
!
!     Number of negative eigenvalues
!
      INDEX = 0
      DO 700 I = 1, NVEC
         IF (EVAL(I) .LT. -THRIND) INDEX = INDEX + 1
  700 CONTINUE
      IF (IPrint.GT. 5) THEN
         WRITE (LUPRI,'(/A,I2/)')    ' Index of Hessian:      ',INDEX
      END IF
!
!     Order non-trarot eigenvalues in increasing order
!
      DO 500 I = 1, NVEC
         JMIN = I
         EMIN = EVAL(I)
         DO 510 J = (I + 1), NVEC
            IF (EVAL(J) .LT. EMIN) THEN
               EMIN = EVAL(J)
               JMIN = J
            END IF
  510    CONTINUE
         IF (JMIN .NE. I) THEN
            CALL DSWAP(1,  EVAL  (I),1,EVAL  (JMIN),1)
            CALL DSWAP(NCR,EVEC(1,I),1,EVEC(1,JMIN),1)
            CALL DSWAP(1,GRDDIA(I),1,GRDDIA(JMIN),1)
         END IF
  500 CONTINUE
      IF (IPrint .GT. 5) THEN
         CALL LSHEADER(lupri,'Sorted eigenvalues in WLKEI1')
         CALL OUTPUT(EVAL,1,1,1,NCR,1,NCR,1,LUPRI)
         CALL LSHEADER(lupri,'Sorted eigenvectors in WLKEI1')
         CALL OUTPUT(EVEC,1,NCR,1,NCR,NCR,NCR,1,LUPRI)
         CALL LSHEADER(lupri,'Gradient (sorted) in WLKEI1')
         CALL OUTPUT(GRDDIA,1,1,1,NCR,1,NCR,1,LUPRI)
      END IF
      RETURN
      END
!  /* Deck wlkcnd */
      SUBROUTINE LS_WLKCND(CNDNMB,EVAL,NONTRO,lupri)
use precision
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      PARAMETER (D0 = 0.0E0_realk, D1 = 1.0E0_realk, D1000 = 1.0E3_realk, DEXP10 = 1.0E10_realk)
      Real(realk) EVAL(NONTRO)
!
!     Condition number of Hessian
!
      EMAX = D0
      EMIN = D1000
      DO 100 I = 1, NONTRO
         EVALI = ABS(EVAL(I))
         EMAX  = MAX(EMAX,EVALI)
         EMIN  = MIN(EMIN,EVALI)
  100 CONTINUE
      IF (EMIN .LT. D1/DEXP10)  THEN 
      WRITE(LUPRI,'(//A,1P,D12.5)') ' Hessian condition &
      & number not calculated since smallest eigenvalue is ',EMIN
         CNDNMB = DEXP10
      ELSE
         CNDNMB = EMAX/EMIN
      END IF
      RETURN
      END
!  /* Deck wlkbis */
      SUBROUTINE ls_WLKBIS(XMAX,XMIN,XDET,GDDIA,HESDIA,TRUSTR, &
                      & NCORD,FUNCT,IFAIL)
!
! Purpose:
!
!  Use bisection to determine zero value of funct(.....)
!
!  Input:  XMAX ; FUNCT(..,XMAX,..) > 0
!          XMIN ; FUNCT(..,XMIN,..) < 0
!  Output: XDET ; FUNCT(..,XDET,..) = 0
!
use precision
Implicit Real(realk) (A-H,O-Z)
      Real(realk) GDDIA(*),HESDIA(*)
      PARAMETER ( MAXIT=100 , D2=2.0E0_realk , D0=0.0E0_realk , DTEST=1.0E-7_realk )
      EXTERNAL FUNCT
      IFAIL = 0
      FMAX  = wstpln_ls(GDDIA,HESDIA,XMAX,NCORD,TRUSTR)
      IF (abs(FMAX) .le. DTEST) THEN
          XDET = XMAX
          GO TO 8000
      END IF
      FMIN  = wstpln_ls(GDDIA,HESDIA,XMIN,NCORD,TRUSTR)
      IF (abs(FMIN) .le. DTEST) THEN
          XDET = XMIN
          GO TO 8000
       END IF
      IF ( (FMAX.LT. 0E0_realk) .OR. (FMIN.GT. 0E0_realk) )  RETURN
!
      ITER = 0
 100  CONTINUE
         ITER= ITER+1
         XDET= (XMAX+XMIN)/D2
         FDET= wstpln_ls(GDDIA,HESDIA,XDET,NCORD,TRUSTR)
         IF (FDET.GT.D0) XMAX=XDET
         IF (FDET.LT.D0) XMIN=XDET
      IF (ABS(XMAX-XMIN).LT.DTEST) GO TO 8000
      IF (ITER .LT. MAXIT) GO TO 100
!
      IFAIL=1
      RETURN
!
8000  CONTINUE
      IFAIL=5
      RETURN
      END
!========================!
!Calc_TMATRIX            !
!========================!
Subroutine Calc_TMatrix(TMatrix,NCoord,Coord)
! A T matrix from: 
! T.Helgaker, Acta Chemica Scandinavica A 42(1998),515-518,
! equation (16), T here is a transpose of the one in the paper
!
! Input:
! 1) TMatrix - inverse of the T matrix introduced in the reference above
! 2) NCoord - number of Cartesian coordinates
! 3) Coord - cartesian coordinates vector
!
Use precision
Implicit none
Integer :: NCoord
Real(realk) Coord(NCoord) 
Real(realk) TMatrix(NCoord,6)
Integer :: i
! 
TMatrix = 0E0_realk
!
! Upper (translational) part of T
!
Do i = 1, NCoord/3
   TMatrix(i*3-2:i*3,1)=(/1E0_realk,0E0_realk,0E0_realk/)
   TMatrix(i*3-2:i*3,2)=(/0E0_realk,1E0_realk,0E0_realk/)
   TMatrix(i*3-2:i*3,3)=(/0E0_realk,0E0_realk,1E0_realk/)
Enddo
!
! Lower (rotational) part of T
!
Do i = 1,NCoord/3
   TMatrix(3*i-2:3*i,4) = (/0E0_realk,-Coord(3*i),Coord(3*i-1)/)
   TMatrix(3*i-2:3*i,5) = (/Coord(3*i),0E0_realk,-Coord(3*i-2)/)
   TMatrix(3*i-2:3*i,6) = (/-Coord(3*i-1),Coord(3*i-2),0E0_realk/)      
Enddo
!
End Subroutine Calc_TMatrix
!
      FUNCTION dnorm2_ls(N,DX,INCX)
!
!     Forms the two-norm of a vector.
! 19-Sep-1988 -- hjaaj -- based on DNRM2 from LINPACK
!     This version does not use extended precision for intermediates
!     as the LINPACK version does.
!     Equivalent to dnorm2_ls in IBM's ESSL library.
!
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     DNRM2: JACK DONGARRA, LINPACK, 3/11/78.
!
use Precision
Implicit Real(realk) (A-H,O-Z)
!
      Real(realk)  DX(*)
      Real(realk),PARAMETER ::  ZERO = 0.0E0_realk 
!
      dnorm2_ls = ZERO
      IF(N.LE. 0)RETURN
      DTEMP  = ZERO
      IF(INCX.EQ. 1)GO TO 20
!
!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1
!
      IX = 1
      IF(INCX.LT. 0)IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DX(IX)
        IX = IX + INCX
   10 CONTINUE
      dnorm2_ls = SQRT(DTEMP)
      RETURN
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP
!
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DX(I)
   30 CONTINUE
      IF( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DX(I) + DX(I + 1)*DX(I + 1) +&
     &   DX(I + 2)*DX(I + 2) + DX(I + 3)*DX(I + 3) + DX(I + 4)*DX(I + 4)
   50 CONTINUE
   60 dnorm2_ls = SQRT(DTEMP)
      RETURN
      END
