
      SUBROUTINE LS_PRIALL(Molecule,NAtoms,CSTEP,CORDNW,lupri,optinfo)
use ls_util 
use optimization_input
use files
use molecule_type
use molecule_module
use precision
!
!     Prints important information for the current iteration.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri, NAtoms
      Type(MoleculeInfo) :: Molecule
      Type(opt_setting) :: optinfo
      Real(realk) CSTEP(MXCOOR), CORDNW(3,MXCENT)
      CHARACTER*1 MARK
      INTEGER     START
!
!     We have to construct the updated geometry before printing it
!
      IJ = 1
      DO 10 J = 1, NAtoms
         DO 20 I = 1, 3
            CORDNW(I,J) = optinfo%Coordinates(I,J) + CSTEP(IJ)
            IJ = IJ + 1
 20      CONTINUE
 10   CONTINUE
!
      call lsheader(lupri,'Optimization Control Center')
      IF (.NOT. optinfo%GeConv) THEN
         call lsheader(lupri,'Next geometry (au)')
         Do i = 1,NAtoms
                  Molecule%Atom(i)%Center(:)=CORDNW(:,i)
         Enddo
         call Print_Geometry(Molecule,lupri)
      ELSE
         call lsheader(lupri,'Final geometry (au)')
!         Do i = 1,NAtoms
!                  Molecule%Atom(i)%Center(:)=CORDNW(:,i)
!         Enddo
         call Print_Geometry(Molecule,lupri)
 !        Do i = 1,NAtoms
 !                 Molecule%Atom(i)%Center(:)=optinfo%Coordinates(:,i)
 !        Enddo
      END IF
      call lsheader(lupri,'Optimization information')
      WRITE(LUPRI,'(A,I7)') &
     &     ' Iteration number               : ',optinfo%ItrNmr
      MARK = ' '
      IF (optinfo%GeConv .AND. optinfo%DoPre .AND. (.NOT. optinfo%FinPre)) MARK = '*'
      WRITE(LUPRI,'(A,L1,A1)') &
     &      ' End of optimization            :       ', optinfo%GeConv, MARK
      WRITE(LUPRI,678) ' Energy at this geometry is     : ', optinfo%energy
      IF ((optinfo%ItrNmr .GT. 0) .AND. (optinfo%ItBreak .LT. (optinfo%ItrNmr - 1))) THEN
         ERGDIF = optinfo%energy - optinfo%energyOld
         IF (ABS(optinfo%predictedChange) .GT. 1.0E-10_realk) THEN
            RATIO = ERGDIF / optinfo%predictedChange
         ELSE
            RATIO = 1.0E0_realk
         END IF
 678     FORMAT(A,F14.6)
         WRITE(LUPRI,678)' Energy change from last geom.  : ', ERGDIF
         IF (optinfo%IPrint .GE. 3) THEN
            WRITE(LUPRI,678)' Predicted change               : ',optinfo%predictedChange
            WRITE(LUPRI,678)' Ratio, actual/predicted change : ',RATIO
         END IF
      END IF
      WRITE(LUPRI,678) ' Norm of gradient               : ', optinfo%GradNorm
      WRITE(LUPRI,678) ' Norm of step                   : ', optinfo%StepNorm
      WRITE(LUPRI,678) ' Updated trust radius           : ', optinfo%TrustRad
!      if (itrnmr.eq. 1) stop ' '
      IF (optinfo%IPrint .GE. 3) THEN
      IF (optinfo%INDHES .GT. 0) WRITE(LUPRI,'(A,I1,A,I8)') &
     &           ' Hessian index               :', optinfo%INDHES
      END IF
      IF (optinfo%GeConv .AND. optinfo%DoPre .AND. (.NOT. optinfo%FinPre)) THEN
         WRITE(LUPRI,'(/A)') ' *) End of preoptimization.'
      END IF
      START = 1
      NCRD = optinfo%NCoordTot
      IF (optinfo%RedInt .OR. optinfo%DelInt) NCRD = optinfo%NIntCoord
      IF (optinfo%RatFun) THEN
         NCRD = NCRD + 1
         START = 2
      END IF
      IF (optinfo%IPrint .GE. 3) THEN
         call lsheader(lupri,'Eigenvalues')
         WRITE(LUPRI,*) &
     &        '  #      Current value  Previous value      Change   '
         WRITE(LUPRI,*) &
     &        '-----------------------------------------------------'
         DO 100 I = START, NCRD
            EVL   = optinfo%EVAL(I)
            EVLOL = optinfo%EVALOL(I)
            IF (EVL .GT. 9.9E3_realk) THEN
               EVL   = 0.0E0_realk
               EVLOL = 0.0E0_realk
            END IF
            NR = I
            IF (optinfo%RatFun) NR = NR - 1
            WRITE(LUPRI,'(I4,3F16.6)') NR,EVL,EVLOL,EVL-EVLOL
 100     CONTINUE
      END IF
      RETURN
      END

!  /* Deck priinf */
      SUBROUTINE LS_PRIINF(Molecule,NAtoms,GEINFO,lupri,optinfo)
use ls_util 
use precision
use optimization_input
use files
use Molecule_type
use molecule_module

!     Prints important information at the end of an optimization.
!
Implicit Real(realk) (A-H,O-Z)
      Type(MoleculeInfo) :: Molecule
      Integer :: lupri, NAtoms
      TYPE(opt_setting) :: optinfo
      Real(realk) GEINFO(0:optinfo%maxIter+1,6)
      CHARACTER NWSMRK*1,BRKMRK*1,REDMRK*1
      real(realk),PARAMETER :: DUMMY = 1.0E20_realk
      integer,parameter :: IDUMMY = - 9999999
      LOGICAL RED
      TMPNRG = 0.0E0_realk
      RED = .FALSE.
      call lsheader(lupri,'Final geometry')
      Do i = 1,NAtoms
               Molecule%Atom(i)%Center(:)=optinfo%Coordinates(:,i)
      Enddo
      call Print_Geometry(Molecule,lupri)
      WRITE(LUPRI,*) &
     & 'Iter      Energy        Change       GradNorm  Index   StepLen &
     &  TrustRad #Rej'
      WRITE(LUPRI,*) &
     & '-------------------------------------------------------------&
     &------------------'
      DO 10 I = 0, optinfo%ItrNmr
         IF (I .GT. 0) TMPNRG = GEINFO(I,1) - GEINFO(I-1,1)
!
!     There's two special marks for each iteration, one (*) for a
!     Newton step (that is a step smaller than the trust radius)
!     and one (#) for dropping half of dihedral angles.
!
         NWSMRK = ' '
         IF ((GEINFO(I,4)+1.0E-6_realk) .LT. GEINFO(I,5)) NWSMRK = '*'
         BRKMRK = ' '
         REDMRK = ' '
         IF (GEINFO(I,6) .LT. -1.0E-3_realk) THEN
            GEINFO(I,6) = ABS(GEINFO(I,6))
            REDMRK = '#'
            RED = .TRUE.
         END IF
         WRITE(LUPRI,'(I4,A1,F15.6,2F13.6,I5,F12.6,A1,F11.6,I4,A1)') &
     &        I,BRKMRK,GEINFO(I,1),TMPNRG,GEINFO(I,2), &
     &        NINT(ABS(GEINFO(I,3))),GEINFO(I,4),NWSMRK,GEINFO(I,5), &
     &        NINT(GEINFO(I,6)),REDMRK
 10   CONTINUE
      RETURN
      END
!==================!
!  Print_scan      !
!==================!
! Prints final results of scanning PES
Subroutine Print_scan(optinfo,IterMax,lupri)
use ls_util 
use precision
use optimization_input
use files
use Molecule_type
use molecule_module
Implicit none
Integer :: lupri,i,IterMax
Type(opt_setting) :: optinfo
      WRITE(LUPRI,*) 
      call lsheader(lupri,'Scan of potential energy surface along internal coordinate')
      WRITE(LUPRI,*) 'Internal coordinate number'
      WRITE(LUPRI,'(I4)') optinfo%Scancoord
      WRITE(LUPRI,*) 
      WRITE(LUPRI,*) &
     & 'Iter      Internal        Energy  '
      WRITE(LUPRI,*) &
     & '----------------------------------------'
   Do i = 0, IterMax
      Write(*,*) optinfo%scan_info(i+1,:)
      Write(LUPRI,'(I4,F14.4,F20.6)') i, optinfo%Scan_info(i+1,1),optinfo%Scan_info(i+1,2)
   Enddo
!
End subroutine Print_scan
!
!  /* Deck inihes */
      SUBROUTINE LS_INIHES(Molecule,NCOOR,MXRCRD,MX2CRD,GRDOLD,HESOLD, &
     &     TMPMAT,TMPMT2,TMPMT3,TMPMT4,WILBMT,BMTRAN,BMTINV, &
     &     HESINT,lupri,optinfo)
use ls_util 
use optimization_input
use files
use precision
use molecule_type
!
!     Initializes Hessian
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri,NCOOR
      TYPE(MoleculeInfo) :: Molecule
      TYPE(opt_setting) :: optinfo
      Real(realk) GRDOLD(MXCOOR)
      Real(realk) HESOLD(NCOOR,NCOOR), TMPMAT(MX2CRD,MX2CRD)
      Real(realk) TMPMT2(MX2CRD,MX2CRD), TMPMT3(MX2CRD,MX2CRD)
      Real(realk) TMPMT4(MX2CRD,MX2CRD), WILBMT(MX2CRD,MX2CRD)
      Real(realk) BMTRAN(MXRCRD,MXRCRD)
      Real(realk) BMTINV(MXRCRD,MXCOOR), HESINT(MXRCRD,MXRCRD)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
!
      call ls_DZERO(HESOLD,NCOOR*NCOOR)
      call ls_DZERO(GRDOLD,MXCOOR)
!
      IF (optinfo%IPrint .GE. IPRMED) THEN
         call lsheader(lupri,'Output from INIHES')
      END IF
!!!!!! Vladimir: cirrently disabled
!     Check if Hessian should be read from file
!
      IF (optinfo%HessFile) THEN
         !
         Call lsquit('No Hessian file implemented in LSDALTON',lupri)
         !
!         IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
!            call ls_GX2GQ(MXRCRD,optinfo%GradMol,optinfo%GRDINT,BMTINV,optinfo)
!            optinfo%GradNorm = SQRT(DDOT(optinfo%ICartCoord,optinfo%GradMol,1,optinfo%GradMol,1))
!         END IF
!         call ls_REAHES(MXRCRD,MX2CRD,HESINT,HESOLD,TMPMAT,TMPMT2,TMPMT3, &
!     &        TMPMT4,WILBMT,BMTRAN,BMTINV,WORK,LWORK,lupri,IERR)
!         call ls_DZERO(HESOLD,MXRCRD*MXRCRD)
!         IF (IERR .EQ. -1) THEN
!            call lsquit('Unable to open the file DALTON.HES.',lupri)
!         ELSE IF (IERR .EQ. -2) THEN
!            call lsquit('The Hessian in DALTON.HES has  &
!                                       &wrong Real(realk)s.',lupri)
!         END IF
!         IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
!            DO 5 J = 1, optinfo%NIntCoord
!               GRDOLD(J) = optinfo%GRDINT(J)
!               DO 7 I = 1, optinfo%NIntCoord
!                  HESOLD(I,J) = HESINT(I,J)
! 7             CONTINUE
! 5          CONTINUE
!         ELSE
!            DO 12 J = 1, optinfo%NCoordTot
!               GRDOLD(J) = optinfo%GradMol(J)
!               DO 14 I = 1, optinfo%NCoordTot
!                  HESOLD(I,J) = optinfo%HessMol(I,J)
!                  optinfo%HessMol(I,J) = HESOLD(I,J)
! 14            CONTINUE
! 12         CONTINUE
!         END IF
!         WRITE(LUPRI,*)
!         WRITE(LUPRI,*) 'Initial Hessian has been read from file.'
!         WRITE(LUPRI,*)
!         optinfo%HessFile = .FALSE.

!
!     Check if Hessian has been calculated
!
      ELSE IF (optinfo%InitHess) THEN
         IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
            call ls_GX2GQ(optinfo%NIntCoord,optinfo%GradMol,optinfo%GRDINT,BMTINV,optinfo)
            optinfo%GradNorm = SQRT(DDOT(optinfo%ICartCoord,optinfo%GradMol,1,optinfo%GradMol,1))
            IF (optinfo%IPrint .GE. IPRDBG) THEN
               call lsheader(lupri,'Cartesian gradient')
               call ls_output(optinfo%GradMol,1,1,1,optinfo%ICartCoord,1,MXCOOR,1,LUPRI)
               call lsheader(lupri,'Internal gradient')
               call ls_output(optinfo%GRDINT,1,1,1,optinfo%NIntCoord,1,MXRCRD,1,LUPRI)
            END IF
            call ls_HX2HQ(Molecule,MXRCRD,MX2CRD,TMPMAT,TMPMT2,TMPMT3,TMPMT4, &
     &          optinfo%HessMol,optinfo%GRDINT,HESINT,WILBMT,BMTINV,BMTRAN, &
     &           optinfo)
            IF (optinfo%IPrint .GE. IPRDBG) THEN
               call lsheader(lupri,'Cartesian Hessian')
               call ls_output(optinfo%HessMol,1,optinfo%ICartCoord,1,optinfo%ICartCoord, &
     &                                   MXCOOR,MXCOOR,1,LUPRI)
            END IF
            DO 15 J = 1, optinfo%NIntCoord
               DO 16 I = 1, optinfo%NIntCoord
                  HESOLD(I,J) = HESINT(I,J)
 16            CONTINUE
               GRDOLD(J) = optinfo%GRDINT(J)
 15         CONTINUE
         ELSE
            DO 18 J = 1, optinfo%NCoordTot
               DO 19 I = 1, optinfo%NCoordTot
                  HESOLD(I,J) = optinfo%HessMol(I,J)
 19            CONTINUE
               GRDOLD(J) = optinfo%GradMol(J)
 18         CONTINUE
         END IF
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Initial Hessian has been calculated.'
         WRITE(LUPRI,*)
!
!     Next line makes reinitialization possible:
!
         optinfo%InitHess = .FALSE.
!
!     Check if initial Hessian should be diagonal in redundant internal
!     coordinates.
!
      ELSE IF (optinfo%RedInt .OR. optinfo%DelInt .OR. optinfo%InmdHess .OR. optinfo%InrdHess) THEN
         NRIC = optinfo%NIntCoord
         IF (optinfo%DelInt) NRIC = optinfo%NIntCoord
         call ls_GX2GQ(optinfo%NIntCoord,optinfo%GradMol,optinfo%GRDINT,BMTINV,optinfo)
         optinfo%GradNorm = SQRT(DDOT(optinfo%ICartCoord,optinfo%GradMol,1,optinfo%GradMol,1))
!
         IF (optinfo%IPrint .GE. IPRDBG) THEN
            call lsheader(lupri,'Cartesian gradient')
            call ls_output(optinfo%GradMol,1,1,1,optinfo%ICartCoord,1,MXCOOR,1,LUPRI)
         END IF
         IF (optinfo%IPrint .GE. IPRMAX) THEN
            call lsheader(lupri,'Internal gradient')
            call ls_output(optinfo%GRDINT,1,1,1,optinfo%NIntCoord,1,MXRCRD,1,LUPRI)
         END IF
!
!     The default is different values for bonds and other internal
!     coordinates.
!
         call ls_DZERO(HESINT,MXRCRD*MXRCRD)
         IF (optinfo%EvLini .LE. 0.0E0_realk) THEN
            IF (optinfo%ModHes .OR. optinfo%InmdHess .OR. optinfo%CMBMod) THEN
               call ls_BLDHES(Molecule,MXRCRD,TMPMAT,HESINT,lupri,optinfo)
            ELSE
               DO 20 I = 1, NRIC
                  IF (optinfo%INTCRD(I,1) .LT. 10) THEN
                     HESINT(I,I) = 0.5E0_realk
                  ELSE IF (optinfo%INTCRD(I,1) .LT. 20) THEN
                     HESINT(I,I) = 0.2E0_realk
                  ELSE
                     HESINT(I,I) = 0.1E0_realk
                  END IF
 20            CONTINUE
            END IF
         ELSE
            DO 22 I = 1, optinfo%NIntCoord
               HESINT(I,I) = optinfo%EvLini
 22         CONTINUE
         END IF
!
!     If we use delocalized internals, we have to transform from
!     redundant to non-redundant space.
!
         IF (optinfo%DelInt .AND. (optinfo%EvLini .LE. 0.0E0_realk)) THEN
            call ls_DZERO(TMPMAT,MX2CRD*MX2CRD)
            DO 200 I = 1, optinfo%NIntCoord
               DO 202 J = 1, NRIC
                  DO 204 K = 1, NRIC
                     TMPMAT(I,J) = TMPMAT(I,J) + BMTRAN(K,I)*HESINT(K,J)
 204              CONTINUE
 202           CONTINUE
 200        CONTINUE
            call ls_DZERO(HESINT,MXRCRD*MXRCRD)
            DO 210 I = 1, optinfo%NIntCoord
               DO 212 J = 1, optinfo%NIntCoord
                  DO 214 K = 1, NRIC
                     HESINT(I,J) = HESINT(I,J) + TMPMAT(I,K)*BMTRAN(K,J)
 214              CONTINUE
 212           CONTINUE
 210        CONTINUE
         END IF
         IF (optinfo%IPrint .GE. IPRMAX) THEN
            call lsheader(lupri,'Internal Hessian')
            call ls_output(HESINT,1,optinfo%NIntCoord,1,optinfo%NIntCoord, &
     &                                       MXRCRD,MXRCRD,1,LUPRI)
         END IF
!
         IF ((optinfo%InrdHess .OR. optinfo%InmdHess) .AND. ((.NOT. optinfo%RedInt) &
     &        .AND. (.NOT. optinfo%DelInt))) THEN
            call ls_DZERO(optinfo%HessMol,optinfo%ICartCoord*optinfo%ICartCoord)
            call ls_HQ2HX(Molecule,MXRCRD,MX2CRD,TMPMAT,TMPMT2,TMPMT3,HESINT, &
     &           optinfo%GRDINT,optinfo%HessMol,WILBMT,BMTRAN,optinfo)
            DO 25 J = 1, optinfo%ICartCoord
               GRDOLD(J) = optinfo%GradMol(J)
               DO 27 I = 1, optinfo%ICartCoord
                  HESOLD(I,J) = optinfo%HessMol(I,J)
 27            CONTINUE
 25         CONTINUE
         ELSE
            DO 30 J = 1, optinfo%NIntCoord
               GRDOLD(J) = optinfo%GRDINT(J)
               DO 32 I = 1, optinfo%NIntCoord
                  HESOLD(I,J) = HESINT(I,J)
 32            CONTINUE
 30         CONTINUE
         END IF
!         optinfo%InmdHess = .FALSE.
!
!     Otherwise the Hessian is set equal to a diagonal matrix. When
!     symmetry is used, the Hessian must be scaled so that it will be
!     correct after normalization.
!
      ELSE
         call ls_DZERO(optinfo%HessMol,MXCOOR*MXCOOR)
         call ls_DZERO(HESOLD,NCOOR*NCOOR)
         DO 40 I = 1, optinfo%NCoordTot
            optinfo%HessMol(I,I) = optinfo%EvLini
            HESOLD(I,I) = optinfo%EvLini
            GRDOLD(I) = optinfo%GradMol(I)
 40      CONTINUE
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Initial Hessian is a diagonal matrix.'
         WRITE(LUPRI,*)
      END IF
!
      IF (optinfo%IPrint .GE. IPRMED) THEN
         call lsheader(lupri,'Initial Hessian')
         IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
            call ls_output(HESOLD,1,optinfo%NIntCoord,1,optinfo%NIntCoord,MXRCRD,MXRCRD, &
     &                                                       1,LUPRI)
         ELSE
            call ls_output(HESOLD,1,optinfo%ICartCoord,1,optinfo%ICartCoord,MXRCRD,MXRCRD, &
     &                                                       1,LUPRI)
         END IF
      END IF
      RETURN
      END

!  /* Deck updhes */
      SUBROUTINE LS_UPDHES(Molecule,MXRCRD,MX2CRD,GRDOLD,GRDMAT,STPMAT, &
     &     HESOLD,GAMMA,TMPMAT,TMPMT2,TMPMT3,TMPMT4,TMPMT5,TMPMT6, &
     &     TMPMT7,TMPMT8,WILBMT,BMTRAN,BMTINV,HESINT,IREJOL,IREJNW, &
     &     lupri,optinfo)
use ls_util 
use precision
use optimization_input
use files
use molecule_type
!
!     Controls Hessian updates
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      TYPE(MoleculeInfo) :: Molecule
      TYPE(opt_setting) :: optinfo
      Real(realk) GRDOLD(MXRCRD)
      Real(realk) GRDMAT(25,MXRCRD), STPMAT(25,MXRCRD)
      Real(realk) HESOLD(MXRCRD,MXRCRD), GAMMA(MXRCRD)
      Real(realk) TMPMAT(MX2CRD*MX2CRD), TMPMT2(MX2CRD*MX2CRD)
      Real(realk) TMPMT3(MX2CRD*MX2CRD), TMPMT4(MX2CRD*MX2CRD)
      Real(realk) TMPMT5(MX2CRD*MX2CRD), TMPMT6(MX2CRD*MX2CRD)
      Real(realk) TMPMT7(MX2CRD*MX2CRD), TMPMT8(MX2CRD*MX2CRD)
      Real(realk) WILBMT(MXRCRD,MXCOOR), BMTRAN(MXRCRD,MXRCRD)
      Real(realk) BMTINV(MXRCRD,MXCOOR), HESINT(MXRCRD,MXRCRD)
      LOGICAL RESET, REJLST
!
!     If redundant internal coordinates are used, we have to transform
!     the calculated Cartesian gradient
!
      IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
         call ls_GX2GQ(optinfo%NIntCoord,optinfo%GradMol,optinfo%GRDINT,BMTINV,optinfo)
         optinfo%GradNorm = SQRT(DDOT(optinfo%ICartCoord,optinfo%GradMol,1,optinfo%GradMol,1))
      END IF
!
!     When a new basis set is used in preoptimization, the Hessian is
!     kept as it is, because no step is taken.
!
      IF (optinfo%KeepHessian) THEN
         call ls_UPDOLD(MXRCRD,GRDOLD,HESOLD,HESINT,lupri,optinfo)
!
!     The requested Hessian update method is used
!
      ELSE
         IF (optinfo%SteepDescend) THEN
            call ls_UPOUTP(MXRCRD,GRDOLD,GAMMA,HESOLD, &
     &           'UPSTPD',lupri,optinfo)
            call ls_INIHES(Molecule,MXRCRD,MXRCRD,MX2CRD,GRDOLD,HESOLD,TMPMAT, &
     &           TMPMT2,TMPMT3,TMPMT4,WILBMT,BMTRAN,BMTINV,HESINT, &
     &           lupri,optinfo)
            IF ((.NOT. optinfo%RedInt) .AND. (.NOT. optinfo%DelInt)) THEN
               DO 10 J = 1, optinfo%ICartCoord
                  DO 12 I = 1, optinfo%ICartCoord
                     optinfo%HessMol(I,J) = optinfo%HessMol(I,J)
 12               CONTINUE
 10            CONTINUE
            END IF
         ELSE IF (optinfo%ModHes) THEN
            call ls_UPOUTP(MXRCRD,GRDOLD,GAMMA,HESOLD, &
     &           'UPMODH',lupri,optinfo)
            ITYPE = 1
            IF (optinfo%DFP) ITYPE = 2
            IF (optinfo%PSB) ITYPE = 3
            IF (optinfo%RankOn) ITYPE = 4
            RESET = .FALSE.
            REJLST = .FALSE.
            IF (optinfo%ItBreak .EQ. (optinfo%ItrNmr-2)) RESET = .TRUE.
!            IF (optinfo%RejIni .AND. (IREJOL .GT. 0)) RESET = .TRUE.
!            IF (IREJNW .GT. optinfo%MaxRej) RESET = .TRUE.
!            IF (IREJOL+IREJNW .GT. 0) REJLST = .TRUE.
            call ls_UPMODH(Molecule,MXRCRD,optinfo%STPINT,STPMAT,GAMMA, &
     &           GRDMAT,HESOLD,HESINT, &
     &           MXRCRD,TMPMAT(1),TMPMAT(MX2CRD+1), &
     &           TMPMT2,TMPMT3,TMPMT4,TMPMT5,TMPMT6,TMPMT7,TMPMT8, &
     &           BMTRAN,optinfo%NIntCoord-optinfo%nProjected,RESET,REJLST,optinfo%GradNorm, &
     &           ITYPE,lupri,optinfo)
         ELSE IF (optinfo%CMBMod) THEN
            call ls_UPOUTP(MXRCRD,GRDOLD,GAMMA,HESOLD,'UPCMBM',lupri,optinfo)
            ITYPE = 1
            IF (optinfo%DFP) ITYPE = 2
            IF (optinfo%PSB) ITYPE = 3
            IF (optinfo%RankOn) ITYPE = 4
            RESET = .FALSE.
            REJLST = .FALSE.
            IF (optinfo%ItBreak .EQ. (optinfo%ItrNmr-2)) RESET = .TRUE.
!            IF (optinfo%RejIni .AND. (IREJOL .GT. 0)) RESET = .TRUE.
            IF (IREJNW .GT. optinfo%MaxRej) RESET = .TRUE.
            IF (IREJOL+IREJNW .GT. 0) REJLST = .TRUE.
            call ls_UPCMBM(Molecule,MXRCRD,optinfo%STPINT,STPMAT, &
     &           GAMMA,GRDMAT,HESOLD,HESINT,optinfo%NIntCoord,optinfo%NIntCoord, &
     &           MXRCRD,TMPMAT(1:MX2CRD),TMPMAT(MX2CRD+1), &
     &           TMPMT2,TMPMT3,TMPMT4,TMPMT5,TMPMT6,TMPMT7,TMPMT8, &
     &           BMTRAN,optinfo%NIntCoord-optinfo%nProjected,RESET,REJLST,optinfo%GradNorm, &
     &           ITYPE,lupri,optinfo)
         ELSE IF (optinfo%Multi) THEN
            call ls_UPOUTP(MXRCRD,GRDOLD,GAMMA,HESOLD,'UPMULT',lupri,optinfo)
            RESET = .FALSE.
            REJLST = .FALSE.
            IF ((optinfo%ItrNmr .LT. 2) .OR. (optinfo%ItBreak .EQ. (optinfo%ItrNmr-2))) &
     &           RESET = .TRUE.
            IF (optinfo%RejIni .AND. (ABS(IREJOL) .GT. 0)) RESET = .TRUE.
            IF (optinfo%RejIni .AND. (IREJOL+IREJNW .GT. 0)) RESET = .TRUE.
            ITYPE = 1
            IF (optinfo%PSB) ITYPE = 2
            IF (optinfo%DFP) ITYPE = 3
            IF (optinfo%RankOn) ITYPE = 4
            IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
               call ls_UPMULT(MXRCRD,optinfo%STPINT,STPMAT,GAMMA,GRDMAT,HESOLD, &
     &              HESINT,optinfo%NIntCoord,MXRCRD,TMPMAT(1),TMPMAT(MX2CRD+1), &
     &              TMPMT2,TMPMT3,TMPMT4,TMPMT5,TMPMT6,optinfo%NIntCoord-optinfo%nProjected, &
     &              RESET,optinfo%GradNorm,ITYPE,.TRUE.,lupri,optinfo)
            ELSE
               call ls_UPMULT(MXRCRD,optinfo%STPSYM,STPMAT,GAMMA,GRDMAT,HESOLD, &
     &              optinfo%HessMol,optinfo%ICartCoord,MXCOOR,TMPMAT(1),TMPMAT(MX2CRD+1), &
     &              TMPMT2,TMPMT3,TMPMT4,TMPMT5,TMPMT6,optinfo%ICartCoord-optinfo%nProjected, &
     &              RESET,optinfo%GradNorm,ITYPE,.TRUE.,lupri,optinfo)
            END IF
         ELSE IF (optinfo%BFGS) THEN 
            call ls_UPOUTP(MXRCRD,GRDOLD,GAMMA,HESOLD,'UPBFGS',lupri,optinfo)
            IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
               call LS_UPBFGS(MXRCRD,optinfo%STPINT,GAMMA,HESOLD,HESINT, &
     & optinfo%NIntCoord ,MXRCRD,TMPMAT,TMPMT2,TMPMT3,optinfo%IPrint,lupri,optinfo)
            ELSE
               call LS_UPBFGS(MXRCRD,optinfo%STPSYM,GAMMA,HESOLD,optinfo%HessMol, &
     &  optinfo%ICartCoord,MXCOOR,TMPMAT,TMPMT2,TMPMT3,optinfo%IPrint,lupri,optinfo)
            END IF
         ELSE IF (optinfo%DFP) THEN 
            call ls_UPOUTP(MXRCRD,GRDOLD,GAMMA,HESOLD,'UPDFP ',lupri,optinfo)
            IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
               call LS_UPDFP(MXRCRD,optinfo%STPINT,GAMMA,HESOLD,HESINT, &
     & optinfo%NIntCoord,MXRCRD,TMPMAT,TMPMT2,TMPMT3,optinfo%IPrint,lupri,optinfo)
            ELSE
               call LS_UPDFP(MXRCRD,optinfo%STPSYM,GAMMA,HESOLD,optinfo%HessMol, &
     & optinfo%ICartCoord,MXCOOR,TMPMAT,TMPMT2,TMPMT3,optinfo%IPrint,lupri,optinfo)
            END IF
         ELSE IF (optinfo%PSB) THEN 
            call ls_UPOUTP(MXRCRD,GRDOLD,GAMMA,HESOLD,'UPPSB ',lupri,optinfo)
            IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
               call LS_UPPSB(MXRCRD,optinfo%STPINT,GAMMA,HESOLD,HESINT, &
     & optinfo%NIntCoord,MXRCRD,TMPMAT,TMPMT2,TMPMT3,optinfo%IPrint,lupri,optinfo)
            ELSE
               call LS_UPPSB(MXRCRD,optinfo%STPSYM,GAMMA,HESOLD,optinfo%HessMol, &
     & optinfo%ICartCoord,MXCOOR,TMPMAT,TMPMT2,TMPMT3,optinfo%IPrint,lupri,optinfo)
            END IF
         ELSE IF (optinfo%RanKon) THEN 
            call ls_UPOUTP(MXRCRD,GRDOLD,GAMMA,HESOLD,'UPRNKO',lupri,optinfo)
            IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
               call ls_UPRNKO(MXRCRD,optinfo%STPINT,GAMMA,HESOLD,HESINT, &
     &          optinfo%NIntCoord,MXRCRD,TMPMAT,TMPMT2,TMPMT3,-1,lupri,optinfo)
            ELSE
               call ls_UPRNKO(MXRCRD,optinfo%STPSYM,GAMMA,HESOLD,optinfo%HessMol, &
     &          optinfo%ICartCoord,MXCOOR,TMPMAT,TMPMT2,TMPMT3,-1,lupri,optinfo)
            END IF
         ELSE IF (optinfo%Bofill) THEN 
            call ls_UPOUTP(MXRCRD,GRDOLD,GAMMA,HESOLD,'UPBOFL',lupri,optinfo)
            IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
               call ls_UPBOFL(MXRCRD,optinfo%STPINT,GAMMA,HESOLD,HESINT, &
     &  optinfo%NIntCoord,MXRCRD,TMPMAT,TMPMT2,TMPMT3,TMPMT4,TMPMT5, &
     &  optinfo%IPrint,lupri,optinfo)
            ELSE
               call ls_UPBOFL(MXRCRD,optinfo%STPSYM,GAMMA,HESOLD,optinfo%HessMol,& 
     &  optinfo%ICartCoord,MXCOOR,TMPMAT,TMPMT2,TMPMT3,TMPMT4,TMPMT5,&
     &  optinfo%IPrint,lupri,optinfo)
            END IF
         ELSE IF (optinfo%BFGSR1) THEN
            call ls_UPOUTP(MXRCRD,GRDOLD,GAMMA,HESOLD,'UPBFR1',lupri,optinfo)
            IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
               call ls_UPBFR1(MXRCRD,optinfo%STPINT,GAMMA,HESOLD,HESINT, &
     &              optinfo%NIntCoord,MXRCRD,TMPMT5,TMPMAT,TMPMT2,&
     &              TMPMT3,TMPMT4,lupri,optinfo)
            ELSE
               call ls_UPBFR1(MXRCRD,optinfo%STPSYM,GAMMA,HESOLD,optinfo%HessMol,optinfo%ICartCoord, &
     &              MXCOOR,TMPMT5,TMPMAT,TMPMT2,TMPMT3,TMPMT4,lupri,optinfo)
            END IF
         ELSE IF (optinfo%Schleg) THEN
            call ls_UPOUTP(MXRCRD,GRDOLD,GAMMA,HESOLD,'UPSCHL',lupri,optinfo)
            RESET = .FALSE.
            IF ((optinfo%ItrNmr .LT. 2) .OR. (optinfo%ItBreak .EQ. (optinfo%ItrNmr-2))) &
     &           RESET = .TRUE.
            IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
               call ls_UPSCHL(MXRCRD,optinfo%STPINT,STPMAT,GAMMA,GRDMAT,HESOLD, &
     &              HESINT,optinfo%NIntCoord,MXRCRD,TMPMAT(1),TMPMAT(MX2CRD+1), &
     &              TMPMT2,TMPMT3,optinfo%NIntCoord-optinfo%nProjected,RESET,optinfo%GradNorm,lupri,optinfo)
            ELSE
               call ls_UPSCHL(MXRCRD,optinfo%STPSYM,STPMAT,GAMMA,GRDMAT,HESOLD, &
     &              optinfo%HessMol,optinfo%ICartCoord,MXCOOR,TMPMAT(1),TMPMAT(MX2CRD+1), &
     &              TMPMT2,TMPMT3,optinfo%ICartCoord-optinfo%nProjected,RESET,optinfo%GradNorm,lupri,optinfo)
            END IF
         END IF
!
!     In all higher symmetries, we set all diagonal elements to 1
!
         IF ((.NOT. optinfo%DelInt) .AND. (.NOT. optinfo%RedInt)) THEN
            DO 15 I = (optinfo%ICartCoord + 1), optinfo%NCoordTot
               optinfo%HessMol(I,I) = 1E0_realk
 15         CONTINUE
!
!     We copy the updated Hessian to HESOLD and the gradient to GRDOLD,
!     making everything ready for the next iteration. The updated
!     Hessian then has to be "off-scaled" before it's sent to MINCGH.
!
            call ls_DZERO (HESOLD,MXRCRD*MXRCRD)
            call ls_DZERO (GRDOLD,MXRCRD)
            DO 20 J = 1, optinfo%ICartCoord
               DO 22 I = 1, optinfo%ICartCoord
                  HESOLD(I,J) = optinfo%HessMol(I,J)
                  optinfo%HessMol(I,J) = optinfo%HessMol(I,J)
 22            CONTINUE
               GRDOLD(J) = optinfo%GradMol(J)
 20         CONTINUE
         ELSE
            call ls_DZERO (HESOLD,MXRCRD*MXRCRD)
            call ls_DZERO (GRDOLD,MXRCRD)
            DO 30 J = 1, optinfo%NIntCoord
               DO 32 I = 1, optinfo%NIntCoord
                  HESOLD(I,J) = HESINT(I,J)
 32            CONTINUE
               GRDOLD(J) = optinfo%GRDINT(J)
 30         CONTINUE
         END IF
      END IF
!
      RETURN
      END

!  /* Deck upoutp */
      SUBROUTINE LS_UPOUTP(MXRCRD,GRDOLD,GAMMA,HESOLD,UPDTXT,lupri,optinfo)
use ls_util 
use precision
use optimization_input
use files
!
!     Some common output for the updating subroutines.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      TYPE(opt_setting) :: optinfo
      Real(realk) GRDOLD(MXRCRD)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      Real(realk) GAMMA(MXRCRD), HESOLD(MXRCRD,MXRCRD) 
      CHARACTER UPDTXT*6, OUTTXT*18
!
      NCOOR = optinfo%ICartCoord
      IF (optinfo%RedInt .OR. optinfo%DelInt) NCOOR = optinfo%NIntCoord
      OUTTXT = 'Output from ' // UPDTXT
      IF (optinfo%IPrint .GE. IPRMIN) THEN
         call lsheader(lupri,OUTTXT)
      END IF
      IF (optinfo%IPrint .GE. IPRMED) THEN
         call lsheader(lupri,'Step from last geometry')
         IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
            call ls_output(optinfo%STPINT,1,1,1,NCOOR,1,MXRCRD,1,LUPRI)
         ELSE
            call ls_output(optinfo%STPSYM,1,1,1,NCOOR,1,MXCOOR,1,LUPRI)
         END IF
         call lsheader(lupri,'Previous Hessian')
         call ls_output(HESOLD,1,NCOOR,1,NCOOR,MXRCRD,MXRCRD,1,LUPRI)
         call lsheader(lupri,'Gradient at last geometry')
         call ls_output(GRDOLD,1,1,1,NCOOR,1,MXRCRD,1,LUPRI)
         call lsheader(lupri,'Gradient at current geometry')
         IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
            call ls_output(optinfo%GRDINT,1,1,1,NCOOR,1,MXRCRD,1,LUPRI)
         ELSE
            call ls_output(optinfo%GradMol,1,1,1,NCOOR,1,MXCOOR,1,LUPRI)
         END IF
      END IF
!
!     The gradient difference (gamma) is calculated. The Cartesian
!     vector has to be scaled, to make it "compatible" with the step
!     vector and the Hessian.
!
      call ls_DZERO(GAMMA,MXRCRD)
      DO 20 I = 1, NCOOR
         IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
            GAMMA(I) = optinfo%GRDINT(I) - GRDOLD(I)
         ELSE
            GAMMA(I) = (optinfo%GradMol(I) - GRDOLD(I))
         END IF
 20   CONTINUE
      IF (optinfo%IPrint .GE. IPRMED) THEN
         call lsheader(lupri,'Scaled gradient difference')
         call ls_output(GAMMA,1,1,1,NCOOR,1,MXRCRD,1,LUPRI)
      END IF
      RETURN
      END

!  /* Deck updold */
      SUBROUTINE LS_UPDOLD(MXRCRD,GRDOLD,HESOLD,HESINT,lupri,optinfo)
use ls_util 
use precision
use optimization_input
use files
!
!     The previous Hessian is copied to MOLHES and GRDOLD is updated.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      TYPE(opt_setting) :: optinfo
      Real(realk) GRDOLD(MXRCRD)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      Real(realk) HESOLD(MXRCRD,MXRCRD), HESINT(MXRCRD,MXRCRD)
!
      IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
         call ls_DZERO(HESINT,MXRCRD*MXRCRD)
         call ls_DZERO(GRDOLD,MXRCRD)
         DO 10 J = 1, optinfo%NIntCoord
            DO 15 I = 1, optinfo%NIntCoord
               HESINT(I,J) = HESOLD(I,J)
 15         CONTINUE
            GRDOLD(J) = optinfo%GRDINT(J)
 10      CONTINUE
      ELSE
         call ls_DZERO(optinfo%HessMol,MXCOOR*MXCOOR)
         call ls_DZERO(GRDOLD,MXRCRD)
         DO 20 J = 1, optinfo%ICartCoord
            DO 25 I = 1, optinfo%ICartCoord
               optinfo%HessMol(I,J) = HESOLD(I,J)
 25         CONTINUE
            GRDOLD(J) = optinfo%GradMol(J)
 20      CONTINUE
         DO 30 I = optinfo%ICartCoord+1, optinfo%NCoordTot
            optinfo%HessMol(I,I) = 1E0_realk
 30      CONTINUE
      END IF
!
      IF (optinfo%IPrint .GE. IPRMED) THEN
         call lsheader(lupri,'Output from UPDOLD')
         IF (optinfo%IPrint .GE. IPRMAX) THEN
            call lsquit('ICRD argument is not defined',-1)
!            call lsheader(lupri,'Unchanged Hessian')
!            call ls_output(HESOLD,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
!            call lsheader(lupri,'Gradient at current geometry')
!            call ls_output(GRDOLD,1,1,1,ICRD,1,MCRD,1,LUPRI)
         ELSE
            WRITE(LUPRI,*) &
     &           'Hessian will not be updated in this iteration.'
         END IF
      END IF
      optinfo%KeepHessian = .FALSE.
      RETURN
      END

!  /* Deck upbfgs */
      SUBROUTINE LS_UPBFGS(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW,ICRD,MCRD, &
     &                      TMPVEC,TMPMAT,TMPMT2,IPrint,lupri,optinfo)
use ls_util 
use precision
use optimization_input
use files
!
!     Updates the Hessian using optinfo%BFGS method
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      TYPE(opt_setting) :: optinfo
      PARAMETER (D0 = 0.0E0_realk)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
!
      Real(realk) DELTA(MXRCRD), GAMMA(MXRCRD)
      Real(realk) HESOLD(MXRCRD,MXRCRD), HESNEW(MCRD,MCRD)
      Real(realk) TMPVEC(MCRD), TMPMAT(MCRD,MCRD), TMPMT2(MCRD,MCRD)
!
!     ******************************************************************
!     The optinfo%BFGS formula is:
!
!                                      T                         T
!                     gamma(n)*gamma(n)    B(n)*delta(n)*delta(n) *B(n)
!     B(n+1) = B(n) + ------------------ - ----------------------------
!                             T                       T
!                     gamma(n) *delta(n)      delta(n) *B(n)*delta(n)
!
!     where
!               gamma(n) = grad(n+1) - grad(n)
!               delta(n) = x(n+1) - x(n)
!
!     The terms in the formula is evaluated one by one below.
!     ******************************************************************
!
!     First we calculate (gamma^T*delta)
!
      FAC = D0
      DO 10 I = 1, ICRD
         FAC = FAC + GAMMA(I)*DELTA(I)
 10   CONTINUE
      IF (IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'(gamma^T*delta)')
         WRITE(LUPRI,'(A,F16.6)') 'Value :    ', FAC
      END IF
!
!     (gamma*gamma^T)/(gamma^T*delta) is calculated
!
      call ls_DZERO(TMPMAT,MCRD*MCRD)
      DO 20 J = 1, ICRD
         DO 22 I = 1, ICRD
            TMPMAT(I,J) = GAMMA(I)*GAMMA(J)/FAC
 22      CONTINUE
 20   CONTINUE
      IF (IPrint .GE. IPRMAX) THEN
         call lsheader(lupri,'(gamma*gamma^T)/(gamma^T*delta)')
         call ls_output(TMPMAT,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
!
!     The first two terms of the optinfo%BFGS formula are placed in TMPMT2.
!
      call ls_DZERO(TMPMT2,MCRD*MCRD)
      DO 30 J = 1, ICRD
         DO 32 I = 1, ICRD
            TMPMT2(I,J) = HESOLD(I,J) + TMPMAT(I,J)
 32      CONTINUE
 30   CONTINUE
      IF (IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'Sum of first two terms')
         call ls_output(TMPMT2,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
!
!     We place (B*delta) in TMPVEC and calculate (delta^T*B*delta)
!
      call ls_DZERO(TMPVEC,MCRD)
      FAC = D0
      DO 40 I = 1, ICRD
         TMP = D0
         DO 42 J = 1, ICRD
            TMP = TMP + HESOLD(I,J)*DELTA(J)
 42      CONTINUE
         TMPVEC(I) = TMP
         FAC = FAC + DELTA(I)*TMP
 40   CONTINUE
      IF (IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'(H*delta)')
         call ls_output(TMPVEC,1,1,1,ICRD,1,MCRD,1,LUPRI)
         call lsheader(lupri,'(delta^T*B*delta)')
         WRITE(LUPRI,'(A,F16.6)') 'Value :    ', FAC
      END IF
!
!     Since the Hessian is symmetric, the elements in the row vector
!     (delta^T*B) must be equal to the elements in the column vector
!     (B*delta). Calculation of (B*delta*delta^T*B)/(delta^T*B*delta)
!     is therefore quite simple.
!
      call ls_DZERO(TMPMAT,MCRD*MCRD)
      DO 50 J = 1, ICRD
         DO 52 I = 1, ICRD
            TMPMAT(I,J) = TMPVEC(I)*TMPVEC(J)/FAC
 52      CONTINUE
 50   CONTINUE
      IF (IPrint .GE. IPRMAX) THEN
         call lsheader(lupri,'(B*delta*delta^T*B)/(delta^T*B*delta)')
         call ls_output(TMPMAT,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
!
!     Finally we obtain the updated Hessian
!
      DO 60 J = 1, ICRD
         DO 62 I = 1, ICRD
            HESNEW(I,J) = TMPMT2(I,J) - TMPMAT(I,J)
 62      CONTINUE
 60   CONTINUE
      IF (IPrint .GE. IPRMIN) THEN
         call lsheader(lupri,'Updated Hessian')
         call ls_output(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

!  /* Deck updfp */
      SUBROUTINE LS_UPDFP(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW,ICRD,MCRD, &
     &                    TMPVEC,TMPMAT,TMPMT2,IPrint,lupri,optinfo)
use ls_util 
use precision
use optimization_input
use files
!
!     Updates the Hessian using optinfo%DFP method
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      TYPE(opt_setting) :: optinfo
      PARAMETER (D0 = 0.0E0_realk)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
!
      Real(realk) DELTA(MXRCRD), GAMMA(MXRCRD)
      Real(realk) HESOLD(MXRCRD,MXRCRD), HESNEW(MCRD,MCRD)
      Real(realk) TMPVEC(MCRD), TMPMAT(MCRD,MCRD), TMPMT2(MCRD,MCRD)
!
!     ******************************************************************
!     The optinfo%DFP formula is:
!
!                                  T                                   T
!                          delta(n) *B(n)*delta(n)    gamma(n)*gamma(n)
!     B(n+1) = B(n) + (1 + -----------------------) * ------------------
!                                    T                        T
!                            gamma(n) *delta(n)       gamma(n) *delta(n)
!
!                               T                              T
!              gamma(n)*delta(n) *B(n) + B(n)*delta(n)*gamma(n)
!            - -------------------------------------------------
!                                    T
!                            gamma(n) *delta(n)
!
!     where
!               gamma(n) = grad(n+1) - grad(n)
!               delta(n) = x(n+1) - x(n)
!
!     The terms in the formula is evaluated one by one below.
!     ******************************************************************
!
!     First we calculate (gamma^T*delta)
!
      FAC1 = D0
      DO 10 I = 1, ICRD
         FAC1 = FAC1 + GAMMA(I)*DELTA(I)
 10   CONTINUE
      IF (IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'(gamma^T*delta)')
         WRITE(LUPRI,'(A,F16.6)') 'Value :    ', FAC1
      END IF
!
!     (delta^T*B*delta) is calculated
!
      FAC2 = D0
      DO 20 I = 1, ICRD
         DUMMY = D0
         DO 22 J = 1, ICRD
            DUMMY = DUMMY + HESOLD(I,J)*DELTA(J)
 22      CONTINUE
         FAC2 = FAC2 + DUMMY*DELTA(I)
 20   CONTINUE
      IF (IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'(delta^T*B*delta)')
         WRITE(LUPRI,'(A,F16.6)') 'Value :    ', FAC2
      END IF
      FAC2 = 1.0E0_realk + FAC2/FAC1
      IF (IPrint .GE. IPRMAX) THEN
         call lsheader(lupri,'(1+(delta^T*B*delta)/(gamma^T*delta))')
         WRITE(LUPRI,'(A,F16.6)') 'Value :    ', FAC2
      END IF
!
!     We calculate (gamma*gamma^T)/(gamma^T*delta)
!
      call ls_DZERO(TMPMAT,MCRD*MCRD)
      DO 25 J = 1, ICRD
         DO 27 I = 1, ICRD
            TMPMAT(I,J) = GAMMA(I)*GAMMA(J)/FAC1
 27      CONTINUE
 25   CONTINUE
      IF (IPrint .GE. IPRMAX) THEN
         call lsheader(lupri,'(gamma*gamma^T)/(gamma^T*delta)')
         call ls_output(TMPMAT,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
!
!     The first two terms of the optinfo%DFP formula is placed in TMPMT2.
!
      call ls_DZERO(TMPMT2,MCRD*MCRD)
      DO 30 J = 1, ICRD
         DO 32 I = 1, ICRD
            TMPMT2(I,J) = HESOLD(I,J) + FAC2*TMPMAT(I,J)
 32      CONTINUE
 30   CONTINUE
      IF (IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'Sum of first two terms')
         call ls_output(TMPMT2,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
!
!     Next is (delta^T*B)
!
      call ls_DZERO(TMPVEC,MCRD)
      DO 40 J = 1, ICRD
         DO 42 I = 1, ICRD
            TMPVEC(J) = TMPVEC(J) + DELTA(I)*HESOLD(I,J)
 42      CONTINUE
 40   CONTINUE
      IF (IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'(delta^T*B)')
         call ls_output(TMPVEC,1,1,1,ICRD,1,MCRD,1,LUPRI)
      END IF
!
!     Then (gamma*delta^T*B + B*delta*gamma^T)
!
      call ls_DZERO(TMPMAT,MCRD*MCRD)
      DO 45 J = 1, ICRD
         DO 47 I = 1, ICRD
            TMPMAT(I,J) = GAMMA(I)*TMPVEC(J) + TMPVEC(I)*GAMMA(J)
 47      CONTINUE
 45   CONTINUE
      IF (IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'(gamma*delta^T*B + B*delta*gamma^T)')
         call ls_output(TMPMAT,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
!
!     And then (gamma*delta^T*B + B*delta*gamma^T)/(gamma^T*delta)
!
      DO 50 J = 1, ICRD
         DO 52 I = 1, ICRD
            TMPMAT(I,J) = TMPMAT(I,J)/FAC1
 52      CONTINUE
 50   CONTINUE
      IF (IPrint .GE. IPRMAX) THEN
         call lsheader(lupri,'(gamma*delta^T*B + B*delta*gamma^T)/(gamma^T*delta)')
         call ls_output(TMPMAT,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
!
!     Finally we obtain the updated Hessian
!
      DO 60 J = 1, ICRD
         DO 62 I = 1, ICRD
            HESNEW(I,J) = TMPMT2(I,J) - TMPMAT(I,J)
 62      CONTINUE
 60   CONTINUE
      IF (IPrint .GE. IPRMIN) THEN
         call lsheader(lupri,'Updated Hessian')
         call ls_output(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

!  /* Deck uppsb */
      SUBROUTINE LS_UPPSB(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW,ICRD,MCRD, &
     &                  TMPVEC,TMPMAT,TMPMT2,IPrint,lupri,optinfo)
use ls_util 
use precision
use optimization_input
use files
!
!     Updates the Hessian using Powell-symmetric-Broyden
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      TYPE(opt_setting) :: optinfo
      PARAMETER (D0 = 0.0E0_realk)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
!
      Real(realk) DELTA(MXRCRD), GAMMA(MXRCRD)
      Real(realk) HESOLD(MXRCRD,MXRCRD), HESNEW(MCRD,MCRD)
      Real(realk) TMPVEC(MCRD), TMPMAT(MCRD,MCRD), TMPMT2(MCRD,MCRD)
!
!     **********************************************************************
!     The Powell-symmetric-Broyden formula is:
!
!                               T                        T
!                      (delta(n) *delta(n))*T(n)*delta(n)
!     B(n+1) = B(n) +  -----------------------------------
!                                     T          2
!                            (delta(n) *delta(n))
!
!                T                        T      T                            T
!       (delta(n) *delta(n))*delta(n)*T(n) -(T(n) *delta(n))*delta(n)*delta(n)
!     + -----------------------------------------------------------------------
!                                        T          2
!                               (delta(n) *delta(n))
!
!     where
!               T(n)     = gamma(n) - B(n)*delta(n)
!               gamma(n) = grad(n+1) - grad(n)
!               delta(n) = x(n+1) - x(n)
!
!     **********************************************************************
!
!     First we calculate (T = gamma-B*delta)
!
      call ls_DZERO(TMPVEC,MCRD)
      DO 10 I = 1, ICRD
         DUMMY = D0
         DO 12 J = 1, ICRD
            DUMMY = DUMMY + HESOLD(I,J)*DELTA(J)
 12      CONTINUE
         TMPVEC(I) = GAMMA(I) - DUMMY
 10   CONTINUE
      IF (IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'(T = gamma-B*delta)')
         call ls_output(TMPVEC,1,1,1,ICRD,1,MCRD,1,LUPRI)
      END IF
!
!     We calculate (delta^T*delta)
!
      FAC1 = D0
      DO 20 I = 1, ICRD
         FAC1 = FAC1 + DELTA(I)*DELTA(I)
 20   CONTINUE
      IF (IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'(delta^T*delta)')
         WRITE(LUPRI,'(A,F16.6)') 'Value :    ', FAC1
      END IF
!
!     Then (T^T*delta)
!
      FAC2 = D0
      DO 22 I = 1, ICRD
         FAC2 = FAC2 + TMPVEC(I)*DELTA(I)
 22   CONTINUE
      IF (IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'(T^T*delta)')
         WRITE(LUPRI,'(A,F16.6)') 'Value :    ', FAC2
      END IF
!
!     and (delta^T*delta)^2
!
      FAC3 = FAC1*FAC1
      IF (IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'(delta^T*delta)^2')
         WRITE(LUPRI,'(A,G16.6)') 'Value :    ', FAC3
      END IF
!
!     Then we construct the complete term
!
      call ls_DZERO(TMPMAT,MCRD*MCRD)
      DO 30 J = 1, ICRD
         DO 32 I = 1, ICRD
            TMPMAT(I,J) = (  FAC1*TMPVEC(I)*DELTA(J) &
     &                     + FAC1*DELTA(I)*TMPVEC(J) &
     &                     - FAC2*DELTA(I)*DELTA(J) )/FAC3
 32      CONTINUE
 30   CONTINUE
      IF (IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'Second term of formula')
         call ls_output(TMPMAT,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
!
!     Finally we obtain the updated Hessian
!
      DO 40 J = 1, ICRD
         DO 45 I = 1, ICRD
            HESNEW(I,J) = HESOLD(I,J) + TMPMAT(I,J)
 45      CONTINUE
 40   CONTINUE
      IF (IPrint .GE. IPRMIN) THEN
         call lsheader(lupri,'Updated Hessian')
         call ls_output(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

!  /* Deck upbofl */
      SUBROUTINE LS_UPBOFL(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW,ICRD,MCRD, &
     &     TMPVEC,TMPMAT,TMPMT2,TMPMT3,TMPMT4,IPrint,lupri,optinfo)
use ls_util 
use precision
use optimization_input
use files
!
!     Updates the Hessian using Bofills update
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      TYPE(opt_setting) :: optinfo
      PARAMETER (D0 = 0.0E0_realk)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
!
      Real(realk) DELTA(MXRCRD), GAMMA(MXRCRD)
      Real(realk) HESOLD(MXRCRD,MXRCRD), HESNEW(MCRD,MCRD)
      Real(realk) TMPVEC(MCRD), TMPMAT(MCRD,MCRD), TMPMT2(MCRD,MCRD)
      Real(realk) TMPMT3(MCRD,MCRD), TMPMT4(MCRD,MCRD)
!
!
!     Bofills update is a linear combination of rank one and optinfo%PSB:
!
!     B(n+1) = phi B_R1 + (1-phi) B_optinfo%PSB
!
      call ls_UPRNKO(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW, &
     &     ICRD,MCRD,TMPMAT,TMPMT2,TMPMT3,-1,lupri,optinfo)
      call LS_UPPSB(MXRCRD,DELTA,GAMMA,HESOLD,TMPMT4, &
     &     ICRD,MCRD,TMPMAT,TMPMT2,TMPMT3,-1,lupri,optinfo)
      call ls_DZERO(TMPVEC,MCRD)
      DO 10 I = 1, ICRD
         DO 12 J = 1, ICRD
            TMPVEC(I) = TMPVEC(I) + HESOLD(I,J)*DELTA(J)
 12      CONTINUE
         TMPVEC(I) = GAMMA(I) - TMPVEC(I)
 10   CONTINUE
!
!     We calculate the factor phi
!
      CPHI = DDOT(ICRD,DELTA,1,TMPVEC,1)
      CPHI = 1.0E0_realk-((CPHI*CPHI)/(DDOT(ICRD,DELTA,1,DELTA,1)* &
     &     DDOT(ICRD,TMPVEC,1,TMPVEC,1)))
      R1PHI = 1.0E0_realk-CPHI
      IF (IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'Old Hessian')
         call ls_output(HESOLD,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
         call lsheader(lupri,'Rank one update')
         call ls_output(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
         call lsheader(lupri,'Powell update')
         call ls_output(TMPMT4,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
         call lsheader(lupri,'delta')
         call ls_output(DELTA,1,1,1,ICRD,1,MCRD,1,LUPRI)
         call lsheader(lupri,'gamma')
         call ls_output(GAMMA,1,1,1,ICRD,1,MCRD,1,LUPRI)
         call lsheader(lupri,'ksi')
         call ls_output(TMPVEC,1,1,1,ICRD,1,MCRD,1,LUPRI)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'phi:',CPHI
         WRITE(LUPRI,*)
      END IF
      DO 20 I = 1, ICRD
         DO 22 J = 1, ICRD
            HESNEW(I,J) = R1PHI*HESNEW(I,J)+CPHI*TMPMT4(I,J)
 22      CONTINUE
 20   CONTINUE
      IF (IPrint .GE. IPRMIN) THEN
         call lsheader(lupri,'Updated Hessian')
         call ls_output(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

!  /* Deck upbfr1 */
      SUBROUTINE LS_UPBFR1(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW,ICRD,MCRD, &
     &     TMPVEC,TMPMAT,TMPMT2,TMPMT3,TMPMT4,lupri,optinfo)
use ls_util 
use precision
use files
use optimization_input
!
!     Updates the Hessian using optinfo%BFGS/rank one combination update (ala Bofill)
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      TYPE(opt_setting) :: optinfo
      PARAMETER (D0 = 0.0E0_realk)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
!
      Real(realk) DELTA(MXRCRD), GAMMA(MXRCRD)
      Real(realk) HESOLD(MXRCRD,MXRCRD), HESNEW(MCRD,MCRD)
      Real(realk) TMPVEC(MCRD), TMPMAT(MCRD,MCRD), TMPMT2(MCRD,MCRD)
      Real(realk) TMPMT3(MCRD,MCRD), TMPMT4(MCRD,MCRD)
!
!     Use a linear combination of rank one and optinfo%BFGS:
!
!     B(n+1) = phi B_R1 + (1-phi) B_optinfo%BFGS
!
      call ls_UPRNKO(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW, &
     &     ICRD,MCRD,TMPMAT,TMPMT2,TMPMT3,-1,lupri,optinfo)
      call ls_UPBFGS(MXRCRD,DELTA,GAMMA,HESOLD,TMPMT4, &
     &     ICRD,MCRD,TMPMAT,TMPMT2,TMPMT3,-1,lupri,optinfo)
      call ls_DZERO(TMPVEC,MCRD)
      DO 10 I = 1, ICRD
         DO 12 J = 1, ICRD
            TMPVEC(I) = TMPVEC(I) + HESOLD(I,J)*DELTA(J)
 12      CONTINUE
         TMPVEC(I) = GAMMA(I) - TMPVEC(I)
 10   CONTINUE
!
!     We calculate the factor phi
!
      CPHI = DDOT(ICRD,DELTA,1,TMPVEC,1)
      CPHI = 1.0E0_realk-((CPHI*CPHI)/(DDOT(ICRD,DELTA,1,DELTA,1)* &
     &     DDOT(ICRD,TMPVEC,1,TMPVEC,1)))
      R1PHI = 1.0E0_realk-CPHI
!
!     Output for debugging
!
      IF (optinfo%IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'Old Hessian')
         call ls_output(HESOLD,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
         call lsheader(lupri,'Rank one update')
         call ls_output(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
         call lsheader(lupri,'optinfo%BFGS update')
         call ls_output(TMPMT4,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
         call lsheader(lupri,'delta')
         call ls_output(DELTA,1,1,1,ICRD,1,MCRD,1,LUPRI)
         call lsheader(lupri,'gamma')
         call ls_output(GAMMA,1,1,1,ICRD,1,MCRD,1,LUPRI)
         call lsheader(lupri,'ksi')
         call ls_output(TMPVEC,1,1,1,ICRD,1,MCRD,1,LUPRI)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'phi:',CPHI
         WRITE(LUPRI,*)
      END IF
      DO 20 I = 1, ICRD
         DO 22 J = 1, ICRD
            HESNEW(I,J) = R1PHI*HESNEW(I,J)+CPHI*TMPMT4(I,J)
 22      CONTINUE
 20   CONTINUE
      IF (optinfo%IPrint .GE. IPRMIN) THEN
         call lsheader(lupri,'Updated Hessian')
         call ls_output(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

!  /* Deck uprnko */
      SUBROUTINE LS_UPRNKO(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW,ICRD,MCRD, &
     &                   TMPVEC,TMPMAT,TMPMT2,IPrint,lupri,optinfo)
use ls_util 
use precision
use optimization_input
use files
!
!     Updates the Hessian using rank one method
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Integer :: IPrint
      TYPE(opt_setting) :: optinfo
      PARAMETER (D0 = 0.0E0_realk)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
!
      Real(realk) DELTA(MXRCRD), GAMMA(MXRCRD)
      Real(realk) HESOLD(MXRCRD,MXRCRD), HESNEW(MCRD,MCRD)
      Real(realk) TMPVEC(MCRD), TMPMAT(MCRD,MCRD), TMPMT2(MCRD,MCRD)
!
!     ******************************************************************
!     The rank one formula is:
!
!                                                                     T
!                     (gamma(n)-B(n)*delta(n))(gamma(n)-B(n)*delta(n))
!     B(n+1) = B(n) + -------------------------------------------------
!                                                   T
!                           (gamma(n)-B(n)*delta(n)) *delta(n)
!
!     where
!               gamma(n) = grad(n+1) - grad(n)
!               delta(n) = x(n+1) - x(n)
!
!     The terms in the formula is evaluated one by one below.
!     ******************************************************************
!
!     First we calculate (gamma-B*delta)
!
      call ls_DZERO(TMPVEC,MCRD)
      DO 10 I = 1, ICRD
         DUMMY = D0
         DO 12 J = 1, ICRD
            DUMMY = DUMMY + HESOLD(I,J)*DELTA(J)
 12      CONTINUE
         TMPVEC(I) = GAMMA(I) - DUMMY
 10   CONTINUE
      IF (IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'(gamma-B*delta)')
         call ls_output(TMPVEC,1,1,1,ICRD,1,MCRD,1,LUPRI)
      END IF
!
!     (gamma-B*delta)(gamma-B*delta)^T is calculated
!
      call ls_DZERO(TMPMAT,MCRD*MCRD)
      DO 20 J = 1, ICRD
         DO 22 I = 1, ICRD
            TMPMAT(I,J) = TMPVEC(I)*TMPVEC(J)
 22      CONTINUE
 20   CONTINUE
      IF (IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'(gamma-B*delta)(gamma-B*delta)^T')
         call ls_output(TMPMAT,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
!
!     ((gamma-B*delta)^T*delta) is calculated
!
      FAC = D0
      DO 25 I = 1, ICRD
         FAC = FAC + TMPVEC(I)*DELTA(I)
 25   CONTINUE
      IF (IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'((gamma-B*delta)^T*delta)')
         WRITE(LUPRI,'(A,F16.6)') 'Value :    ', FAC 
      END IF
!
!     And we complete the term
!     ((gamma-B*delta)(gamma-B*delta)^T)/((gamma-B*delta)^T*delta)
!
      DO 30 J = 1, ICRD
         DO 32 I = 1, ICRD
            TMPMAT(I,J) = TMPMAT(I,J)/FAC
 32      CONTINUE
 30   CONTINUE
      IF (IPrint .GE. IPRMAX) THEN
         call lsheader(lupri,'((gamma-B*delta)(gamma-B*delta)^T)/((gamma-B*delta)^T*delta)')
         call ls_output(TMPMAT,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
!
!     Finally we obtain the updated Hessian
!
      call ls_DZERO(HESNEW,MCRD*MCRD)
      DO 40 J = 1, ICRD
         DO 42 I = 1, ICRD
            HESNEW(I,J) = HESOLD(I,J) + TMPMAT(I,J)
 42      CONTINUE
 40   CONTINUE
      IF (IPrint .GE. IPRMIN) THEN
         call lsheader(lupri,'Updated Hessian')
         call ls_output(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

!  /* Deck upschl */
      SUBROUTINE LS_UPSCHL(MXRCRD,DELTA,STPMAT,GAMMA,GRDMAT,HESOLD,HESNEW, &
     &     ICRD,MCRD,TMPVEC,TMPVC2,RMAT,BMAT,MAXDIM,RESET,GNRM,lupri,optinfo)
use ls_util 
use precision
use optimization_input
use files
!
!     Updates the Hessian using Schlegel's updating scheme
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      TYPE(opt_setting) :: optinfo
      PARAMETER (D0 = 0.0E0_realk)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
!
      Real(realk) DELTA(MXRCRD), STPMAT(25,MCRD)
      Real(realk) GAMMA(MXRCRD), GRDMAT(25,MCRD)
      Real(realk) HESOLD(MXRCRD,MXRCRD), HESNEW(MCRD,MCRD)
      Real(realk) TMPVEC(MCRD), TMPVC2(MCRD)
      Real(realk) RMAT(MCRD,MCRD), BMAT(MCRD,MCRD)
      LOGICAL RESET
!
      call ls_DZERO(TMPVEC,MCRD)
      call ls_DZERO(TMPVC2,MCRD)
      call ls_DZERO(RMAT,MCRD*MCRD)
      call ls_DZERO(BMAT,MCRD*MCRD)
!
!     First we have to transfer the last step and gradient difference to
!     STPMAT and GRDMAT respectively.
!
      IF (RESET) THEN
         call ls_DZERO(GRDMAT,25*MCRD)
         call ls_DZERO(STPMAT,25*MCRD)
         INUM = 1
!
!     After 25 iterations we have to discard the first entries.
!
      ELSE IF (STPMAT(25,1) .GT. 1.0E10_realk) THEN
         DO 10 I = 1, 24
            DO 12 J = 1, ICRD
               STPMAT(I,J) = STPMAT(I+1,J)
               GRDMAT(I,J) = GRDMAT(I+1,J)
 12         CONTINUE
 10      CONTINUE
         STPMAT(25,1) = 0.0E0_realk
         INUM = 24
      ELSE
         INUM = 1
 14      CONTINUE
         IF (STPMAT(INUM,1) .LT. 1.0E10_realk) THEN
            INUM = INUM + 1
            GOTO 14
         END IF
!
!     We also have to check the Real(realk) of our variational space,
!     if the number of displacement vectors exceeds this number, we
!     have to remove some.
!
         IF (INUM .GT. MAXDIM) THEN
            DO 15 I = 1, INUM -1
               DO 16 J = 1, ICRD
                  STPMAT(I,J) = STPMAT(I+1,J)
                  GRDMAT(I,J) = GRDMAT(I+1,J)
 16            CONTINUE
 15         CONTINUE
            INUM = INUM - 1
         END IF
      END IF
!
!     Then we update all vectors
!
      SNMLST = 0.0E0_realk
      DO 17 I = 1, ICRD
         STPMAT(INUM,I) = -DELTA(I)
         GRDMAT(INUM,I) = -GAMMA(I)
         SNMLST = SNMLST + DELTA(I)*DELTA(I)
         DO 19 II = 1, INUM-1
            STPMAT(II,I) = STPMAT(II,I)-DELTA(I)
            GRDMAT(II,I) = GRDMAT(II,I)-GAMMA(I)
 19      CONTINUE
 17   CONTINUE
      SNMLST = SQRT(SNMLST)
      STPMAT(INUM+1,1) = 1.1E10_realk
!
!     **********************************************************************
!     Schlegel's scheme goes as follows:
!
!                    i-1            t
!     r' = (x -x )-  SUM  r ((x -x ) r ) ;   r = r'/|r'| ;   j = i-1, i-2, ...
!      j     j  i   m=j+1  m   j  i   m       j   j   j
!
!                    t      i-1             t               t
!     b   = [ (g -g ) r  -  SUM  b  ((x -x ) r ) ] / (x -x ) r  ;
!      jk       j  i   k   m=j+1  mk   j  i   m        j  i   j
!
!     b   = b   ;   j <= k = i-1, i-2, ...
!      kj    jk
!
!                                       t
!     B  = B   + SUM (b  - r B   r ) r r
!      i    i-1  j k   jk   j i-1 k   j k
!
!     **********************************************************************
!
!     We start by calculating the orthonormal displacement vectors r(j)
!
      DO 20 II = INUM, 1, -1
         TMPNRM = 0.0E0_realk
         DO 22 I = 1, ICRD
            RMAT(II,I) = STPMAT(II,I)
            TMPNRM = TMPNRM + STPMAT(II,I)*STPMAT(II,I)
            TMPVC2(I) = RMAT(II,I)
 22      CONTINUE
         TMPNRM = SQRT(TMPNRM)
         DO 26 I = II+1, INUM
            TMP = 0.0E0_realk
            DO 28 J = 1, ICRD
               TMP = TMP + STPMAT(II,J)*RMAT(I,J)
 28         CONTINUE
            IF(optinfo%IPrint .GE. IPRDBG) THEN
               call lsheader(lupri,'(x(j)-x(i))^T * r(m)')
               WRITE(LUPRI,*) TMP
            END IF
            DO 30 J = 1, ICRD
               RMAT(II,J) = RMAT(II,J) - RMAT(I,J)*TMP
               TMPVC2(J) = RMAT(II,J)
 30         CONTINUE
 26      CONTINUE
         IF(optinfo%IPrint .GE. IPRDBG) THEN
            call lsheader(lupri,'r(j)''')
            call ls_output(TMPVC2,1,1,1,ICRD,1,MCRD,1,LUPRI)
         END IF
         TMP = 0.0E0_realk
         DO 32 J = 1, ICRD
            TMP = TMP + RMAT(II,J)*RMAT(II,J)
 32      CONTINUE
         TMP = SQRT(TMP)
!
!     If a vector is in a space already spanned by other vectors
!     (that is less than 15% of the vector is outside this space),
!     all elements are set equal to a very small number.
!
         IF ((TMP .LT. 0.15E0_realk*TMPNRM) .OR. &
     &        (TMPNRM .GE. 1.0E2_realk*SNMLST) .OR. &
     &        (GNRM .GE. 0.1E0_realk)) THEN
            GRDMAT(II,1) = 1.1E10_realk
            DO 35 J = 1, ICRD
               RMAT(II,J) = 1.0E-25_realk
 35         CONTINUE
            TMP = 1.0E0_realk
         END IF
         DO 40 J = 1, ICRD
            RMAT(II,J) = RMAT(II,J)/TMP
            TMPVC2(J)  = RMAT(II,J)
 40      CONTINUE
         IF (optinfo%IPrint .GE. IPRDBG) THEN
            call lsheader(lupri,'Original displacement vector')
            call ls_output(STPMAT,II,II,1,ICRD,25,MXRCRD,1,LUPRI)
            call lsheader(lupri,'Orthonormalized vector')
            call ls_output(TMPVC2,1,1,1,ICRD,1,MCRD,1,LUPRI)
         END IF
 20   CONTINUE
!
!     All redundant vectors are removed
!
      DO 45 I = INUM, 1, -1
         IF (GRDMAT(I,1) .GT. 1.0E10_realk) THEN
            DO 47 II = I, INUM
               DO 48 III = 1, ICRD
                  RMAT(II,III)   = RMAT(II+1,III)
                  GRDMAT(II,III) = GRDMAT(II+1,III)
                  STPMAT(II,III) = STPMAT(II+1,III)
 48            CONTINUE
 47         CONTINUE
            INUM = INUM -1
         END IF
 45   CONTINUE
      IF (optinfo%IPrint .GE. IPRMAX) THEN
         call lsheader(lupri,'Matrix of orthonormalized vectors')
         call ls_output(RMAT,1,INUM,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
!
!     Next task is to calculate the matrix elements b(jk)
!
      DO 50 K = INUM, 1, -1
         DO 52 J = K, 1, -1
            IF(optinfo%IPrint .GE. IPRDBG) THEN
               call lsheader(lupri,'g(j)-g(i)')
               call ls_output(GRDMAT,J,J,1,ICRD,25,MCRD,1,LUPRI)
            END IF
            BMAT(J,K) = 0.0E0_realk
            DO 58 I = 1, ICRD
               BMAT(J,K) = BMAT(J,K) + GRDMAT(J,I)*RMAT(K,I)
 58         CONTINUE
!
            IF(optinfo%IPrint .GE. IPRDBG) THEN
               call lsheader(lupri,'x(j)-x(i)')
               call ls_output(STPMAT,J,J,1,ICRD,25,MCRD,1,LUPRI)
            END IF
!
            DO 70 M = J+1, INUM
               TMP = 0.0E0_realk
               DO 72 I = 1, ICRD
                  TMP = TMP + STPMAT(J,I)*RMAT(M,I)
 72            CONTINUE
               BMAT(J,K) = BMAT(J,K) - BMAT(M,K)*TMP
 70         CONTINUE
            TMP = 0.0E0_realk
            DO 75 I = 1, ICRD
               TMP = TMP + STPMAT(J,I)*RMAT(J,I)
 75         CONTINUE
            BMAT(J,K) = BMAT(J,K)/TMP
            BMAT(K,J) = BMAT(J,K)
 52      CONTINUE
 50   CONTINUE
      IF (optinfo%IPrint .GE. IPRMAX) THEN
         call lsheader(lupri,'B-matrix of coefficients')
         call ls_output(BMAT,1,INUM,1,INUM,MCRD,MCRD,1,LUPRI)
      END IF
!
!     And finally we are ready to update the Hessian
!
!     The elements b(jk) are replaced by b(jk)-r(j)^T*B*r(k)
!
      DO 80 K = 1, INUM
         DO 82 I = 1, ICRD
            TMPVEC(I) = 0.0E0_realk
            DO 83 II = 1, ICRD
               TMPVEC(I) = TMPVEC(I) + HESOLD(I,II)*RMAT(K,II)
 83         CONTINUE
 82      CONTINUE
         IF(optinfo%IPrint .GE. IPRDBG) THEN
            call lsheader(lupri,'B*r(k)^T')
            call ls_output(TMPVEC,1,1,1,ICRD,1,MCRD,1,LUPRI)
         END IF
         DO 85 J = 1, INUM
            TMP = 0.0E0_realk
            DO 87 I = 1, ICRD
               TMP = TMP + RMAT(J,I)*TMPVEC(I)
 87         CONTINUE
            BMAT(J,K) = BMAT(J,K) - TMP
 85      CONTINUE
 80   CONTINUE
      IF (optinfo%IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'b(jk)-r(j)^T*B*r(k)')
         call ls_output(BMAT,1,INUM,1,INUM,MCRD,MCRD,1,LUPRI)
      END IF
!
!     The old Hessian is moved to HESNEW, so that HESOLD can be used for
!     temporary storage.
!
      call ls_DZERO(HESNEW,MCRD*MCRD)
      DO 90 J = 1, ICRD
         DO 92 I = 1, ICRD
            HESNEW(I,J) = HESOLD(I,J)
 92      CONTINUE
 90   CONTINUE
      call ls_DZERO(TMPVEC,MCRD)
!
!     Then we construct the matrices r(j)*r(k)^T
!
      DO 100 J = 1, INUM
         DO 110 K = 1, INUM
            call ls_DZERO(HESOLD,MXRCRD*MXRCRD)
            DO 114 IJ = 1, ICRD
               DO 116 IK = 1, ICRD
                  HESOLD(IJ,IK) = RMAT(J,IJ)*RMAT(K,IK)
 116           CONTINUE
 114        CONTINUE
            IF (optinfo%IPrint .GE. IPRDBG) THEN
               call lsheader(lupri,'r(j)*r(k)^T')
               call ls_output(HESOLD,1,ICRD,1,ICRD,MXRCRD,MXRCRD, &
     &              1,LUPRI)
            END IF
            DO 120 II = 1, ICRD
               DO 122 JJ = 1, ICRD
                  HESNEW(II,JJ) = HESNEW(II,JJ)+BMAT(J,K)*HESOLD(II,JJ)
 122           CONTINUE
 120        CONTINUE
 110     CONTINUE
 100  CONTINUE
      IF (optinfo%IPrint .GE. IPRMIN) THEN
         call lsheader(lupri,'Updated Hessian')
         call ls_output(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

!  /* Deck upmult */
      SUBROUTINE LS_UPMULT(MXRCRD,DELTA,STPMAT,GAMMA,GRDMAT,HESOLD,HESNEW, &
     &     ICRD,MCRD,TMPVEC,TMPVC2,TMPMAT,TMPMT2,HESUPD,RMAT,GMAT, &
     &     MAXDIM,RESET,GNRM,ITYPE,SMART,lupri,optinfo)
use ls_util 
use precision
use optimization_input
use files
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      TYPE(opt_setting) :: optinfo
      PARAMETER (D0 = 0.0E0_realk)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      Real(realk) DELTA(MXRCRD), STPMAT(25,MCRD)
      Real(realk) GAMMA(MXRCRD), GRDMAT(25,MCRD)
      Real(realk) HESOLD(MXRCRD,MXRCRD), HESNEW(MCRD,MCRD)
      Real(realk) TMPVEC(MCRD), TMPVC2(MCRD)
      Real(realk) TMPMAT(MXRCRD*MXRCRD), TMPMT2(MXRCRD*MXRCRD)
      Real(realk) HESUPD(MXRCRD,MXRCRD)
      Real(realk) RMAT(25,MCRD), GMAT(25,MCRD)
      LOGICAL RESET,SMART
!
      ILIM = 25
      IF (.NOT. SMART) ILIM = 5
!
!     ITYPE indicates the update method:
!
!     1 - optinfo%BFGS (default)
!     2 - optinfo%PSB
!     3 - optinfo%DFP
!     4 - Rank one
!
      IF ((ITYPE .LT. 1) .OR. (ITYPE .GT. 4)) ITYPE = 1
!
!     Since STPMAT and GRDMAT are used to store a maximum of 25 vectors,
!     MCRD should be equal to or larger than 25
!
      MXDIM = MAXDIM
      IF (MCRD .LT. 25) MXDIM = MIN(MAXDIM,MCRD)
!
      call ls_DZERO(TMPVEC,MCRD)
      call ls_DZERO(TMPVC2,MCRD)
      call ls_DZERO(GMAT,25*MCRD)
      call ls_DZERO(RMAT,25*MCRD)
      call ls_DZERO(HESUPD,MXRCRD*MXRCRD)
!
!     First we have to transfer the last step and gradient difference to
!     STPMAT and GRDMAT respectively.
!
      IF (RESET) THEN
         call ls_DZERO(GRDMAT,25*MCRD)
         call ls_DZERO(STPMAT,25*MCRD)
         INUM = 1
!
!     After ILIM iterations we discard the first entries.
!
      ELSE IF (STPMAT(ILIM,1) .GT. 1.0E10_realk) THEN
         DO 10 I = 1, ILIM-1
            DO 12 J = 1, ICRD
               STPMAT(I,J) = STPMAT(I+1,J)
               GRDMAT(I,J) = GRDMAT(I+1,J)
 12         CONTINUE
 10      CONTINUE
         STPMAT(ILIM,1) = 0.0E0_realk
         INUM = ILIM-1
      ELSE
         INUM = 1
 14      CONTINUE
         IF (STPMAT(INUM,1) .LT. 1.0E10_realk) THEN
            INUM = INUM + 1
            GOTO 14
         END IF
!
!     We also have to check the Real(realk) of our variational space,
!     if the number of displacement vectors exceeds this number, we
!     have to remove some.
!
         IF (INUM .GT. MXDIM) THEN
            DO 15 I = 1, INUM-1
               DO 16 J = 1, ICRD
                  STPMAT(I,J) = STPMAT(I+1,J)
                  GRDMAT(I,J) = GRDMAT(I+1,J)
 16            CONTINUE
 15         CONTINUE
            INUM = INUM - 1
            IF (optinfo%IPrint .GE. IPRDBG) THEN
               WRITE(LUPRI,*)
               WRITE(LUPRI,*) 'Removing one ' // &
     &              'displacement due to Real(realk) of variational space'
            END IF
         END IF
      END IF
!
!     Then we update all vectors and calculate the norm of the
!     last displacement (SNMLST). The vectors are also copied
!     to RMAT and GMAT.
!
      SNMLST = 0.0E0_realk
      DO 17 I = 1, ICRD
         STPMAT(INUM,I) = -DELTA(I)
         GRDMAT(INUM,I) = -GAMMA(I)
         RMAT(INUM,I) = STPMAT(INUM,I)
         GMAT(INUM,I) = GRDMAT(INUM,I)
         SNMLST = SNMLST + DELTA(I)*DELTA(I)
         DO 19 II = 1, INUM-1
            STPMAT(II,I) = STPMAT(II,I)-DELTA(I)
            GRDMAT(II,I) = GRDMAT(II,I)-GAMMA(I)
            RMAT(II,I) = STPMAT(II,I)
            GMAT(II,I) = GRDMAT(II,I)
 19      CONTINUE
 17   CONTINUE
      SNMLST = SQRT(SNMLST)
      STPMAT(INUM+1,1) = 1.1E10_realk
      RMAT(INUM+1,1) = 1.1E10_realk
      call lsheader(lupri,'Original displacements')
      call ls_output(STPMAT,1,INUM,1,ICRD,25,MCRD,1,LUPRI)
      call lsheader(lupri,'Original gradient vectors')
      call ls_output(GRDMAT,1,INUM,1,ICRD,25,MCRD,1,LUPRI)
!
!     Displacements larger than 100 times the last step,
!     and displacements with a gradient difference
!     larger than 0.25E0_realk are marked (later to be removed).
!
      IF (SMART) THEN
         DO 20 II = 1, INUM
            DNRM = D0
            GNRM = D0
            DO 22 I = 1, ICRD
               DNRM = DNRM + RMAT(II,I)*RMAT(II,I)
               GNRM = GNRM + GMAT(II,I)*GMAT(II,I)
 22         CONTINUE
            DNRM = SQRT(DNRM)
            GNRM = SQRT(GNRM)
            IF (DNRM .GE. 100E0_realk*SNMLST) THEN
!            IF (DNRM .GE. 500E0_realk*SNMLST) THEN
               IF (optinfo%IPrint .GE. IPRDBG) WRITE(LUPRI,*) &
     &              'Removing one displacement due to distance'
               RMAT(II,1) = -1.1E10_realk
            ELSE IF (GNRM .GE. 0.25E0_realk) THEN
               IF (optinfo%IPrint .GE. IPRDBG) WRITE(LUPRI,*) &
     &              'Removing one displacement due to gradient norm'
               RMAT(II,1) = -1.1E10_realk
            END IF
 20      CONTINUE
         II = 1
 25      CONTINUE
         IF ((RMAT(II,1) .LE. -1.0E10_realk) .AND. (INUM .GE. 1)) THEN
            DO 27 JJ = II, INUM - 1
               DO 29 J = 1, ICRD
                  RMAT(JJ,J) = RMAT(JJ + 1,J)
                  GMAT(JJ,J) = GMAT(JJ + 1,J)
 29            CONTINUE
 27         CONTINUE
            INUM = INUM - 1
            GOTO 25
         ELSE IF (II .LT. INUM-1) THEN
            II = II + 1
            GOTO 25
         END IF
      END IF
!
!     We check if the last displacement is nearly parallell to an
!     earlier step (dot product larger than 0.75). If that is the case,
!     the older step is removed.
!
      IF ((INUM .GT. 1) .AND. SMART) THEN
         II = 1
 30      CONTINUE
         DOTP = D0
         SNRM1 = D0
         SNRM2 = D0
         DO 32 I = 1, ICRD
            DOTP = DOTP + RMAT(II,I)*RMAT(INUM,I)
            SNRM1 = SNRM1 + RMAT(II,I)*RMAT(II,I)
            SNRM2 = SNRM2 + RMAT(INUM,I)*RMAT(INUM,I)
 32      CONTINUE
         SNRM1 = SQRT(SNRM1)
         SNRM2 = SQRT(SNRM2)
         IF ((SNRM1*SNRM2) .GE. 1.0E-15_realk) THEN
            DOTP = DOTP/(SNRM1*SNRM2)
         ELSE
            DOTP = D0
         END IF
!         IF (DOTP .GE. 0.80E0_realk) THEN
         IF (DOTP .GE. 0.90E0_realk) THEN
            DO 34 J = 1, ICRD
               DO 36 JJ = II, INUM - 1
                  RMAT(JJ,J) = RMAT(JJ+1,J)
                  GMAT(JJ,J) = GMAT(JJ+1,J)
 36            CONTINUE
 34         CONTINUE
            INUM = INUM - 1
            IF (optinfo%IPrint .GE. IPRDBG) WRITE(LUPRI,*) &
     &           'Removing one nearly parallell and older displacement'
            IF (II .LT. INUM) GOTO 30
         ELSE IF (II .LT. INUM-1) THEN
            II = II + 1
            GOTO 30
         END IF
      END IF
!
!     Some output
!
      IF (optinfo%IPrint .GE. IPRMAX) THEN
         call lsheader(lupri,'Displacements')
         call ls_output(RMAT,1,INUM,1,ICRD,25,MCRD,1,LUPRI)
         call lsheader(lupri,'Gradient vectors')
         call ls_output(GMAT,1,INUM,1,ICRD,25,MCRD,1,LUPRI)
      END IF
!
!     Finally we use do the updating, suppressing output.
!
      DO 60 K = 1, INUM
         call ls_DZERO(DELTA,MXRCRD)
         call ls_DZERO(GAMMA,MXRCRD)
         DO 62 I = 1, ICRD
             DELTA(I) = -RMAT(K,I)
             GAMMA(I) = -GMAT(K,I)
 62      CONTINUE
         IF (ITYPE .EQ. 1) THEN
            call ls_UPBFGS(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW, &
     &           ICRD,MCRD,TMPVEC,TMPMAT,TMPMT2,-1,lupri,optinfo)
         ELSE IF (ITYPE .EQ. 2) THEN
            call ls_UPPSB(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW, &
     &           ICRD,MCRD,TMPVEC,TMPMAT,TMPMT2,-1,lupri,optinfo)
         ELSE IF (ITYPE .EQ. 3) THEN
            call ls_UPDFP(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW, &
     &           ICRD,MCRD,TMPVEC,TMPMAT,TMPMT2,-1,lupri,optinfo)
         ELSE
            call ls_UPRNKO(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW, &
     &           ICRD,MCRD,TMPVEC,TMPMAT,TMPMT2,-1,lupri,optinfo)
         END IF
         IF ((optinfo%IPrint .GE. IPRDBG) .AND. (K .LT. INUM)) THEN
            call lsheader(lupri,'Updating the Hessian')
            call ls_output(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)         
         END IF
!
!     The various contributions are collected in UPDHES.
!
         DO 65 J = 1, ICRD
            DO 67 I = 1, ICRD
               HESUPD(I,J) = HESUPD(I,J) + (HESNEW(I,J) - HESOLD(I,J))
 67         CONTINUE
 65      CONTINUE
         IF (optinfo%IPrint .GE. IPRDBG) THEN
            call lsheader(lupri,'Update')
            call ls_output(HESUPD,1,ICRD,1,ICRD,MXRCRD,MXRCRD,1,LUPRI)
         END IF
 60   CONTINUE
      IF (optinfo%IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'Total Update')
         call ls_output(HESUPD,1,ICRD,1,ICRD,MXRCRD,MXRCRD,1,LUPRI)
      END IF
      IF (INUM .GE. 1) THEN
         FAC = 1.0E0_realk/(1.0E0_realk*INUM)
         DO 70 J = 1, ICRD
            DO 72 I = 1, ICRD
               HESUPD(I,J) = HESUPD(I,J)*FAC
               HESNEW(I,J) = HESOLD(I,J) + HESUPD(I,J)
 72         CONTINUE
 70      CONTINUE
      END IF
      IF (optinfo%IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'Scaled Update')
         call ls_output(HESUPD,1,ICRD,1,ICRD,MXRCRD,MXRCRD,1,LUPRI)
      END IF
      IF (optinfo%IPrint .GE. IPRMIN) THEN
         call lsheader(lupri,'Updated Hessian')
         call ls_output(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

!  /* Deck upcmbm */
      SUBROUTINE LS_UPCMBM(Molecule,MXRCRD,DELTA,STPMAT,GAMMA,GRDMAT,HESOLD,HESNEW, &
     &     ICRD,IRCRD,MCRD,TMPVEC,TMPVC2,RMAT,GMAT,TMPMAT,TMPMT2,TMPMT3, &
     &     TMPMT4,TMPMT5,BMTRAN,MAXDIM,RESET,REJLST,GNRM, &
     &     ITYPE,lupri,optinfo)
use ls_util 
use precision
use optimization_input
use files
use Molecule_type
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      TYPE(MoleculeInfo) :: Molecule
      TYPE(opt_setting) :: optinfo
      PARAMETER (D0 = 0.0E0_realk)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
!
      Real(realk) DELTA(MXRCRD), STPMAT(25,MCRD)
      Real(realk) GAMMA(MXRCRD), GRDMAT(25,MCRD)
      Real(realk) HESOLD(MXRCRD,MXRCRD), HESNEW(MCRD,MCRD)
      Real(realk) TMPVEC(MCRD), TMPVC2(MCRD)
      Real(realk) TMPMAT(MXRCRD,MXRCRD), TMPMT2(MXRCRD*MXRCRD)
      Real(realk) TMPMT3(MXRCRD,MXRCRD), TMPMT4(MXRCRD,MXRCRD)
      Real(realk) TMPMT5(MXRCRD,MXRCRD), BMTRAN(MXRCRD,MXRCRD)
      Real(realk) RMAT(MCRD,MCRD), GMAT(MCRD,MCRD)
      LOGICAL RESET, REJLST
      SAVE IFAC
      DATA IFAC /2/
      IF (RESET) IFAC = 1
      IF (RESET .OR. REJLST) IFAC = 1
!
!     We do the updating, suppressing output.
!
      call ls_UPBFGS(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW, &
     &     ICRD,MCRD,TMPMAT,TMPMT2,TMPMT3,optinfo%IPrint,lupri,optinfo)
!
!     A new model Hessian is calculated
!
      call ls_BLDHES(Molecule,MXRCRD,TMPMAT,TMPMT5,lupri,optinfo)
!
!     If we use delocalized internals, we have to transform from
!     redundant to non-redundant space.
!
      IF (optinfo%DelInt) THEN
         call ls_DZERO(TMPMAT,MXRCRD*MXRCRD)
         DO 200 I = 1, ICRD
            DO 202 J = 1, IRCRD
               DO 204 K = 1, IRCRD
                  TMPMAT(I,J) = TMPMAT(I,J) + BMTRAN(K,I)*TMPMT5(K,J)
 204           CONTINUE
 202        CONTINUE
 200     CONTINUE
         call ls_DZERO(TMPMT5,MXRCRD*MXRCRD)
         DO 210 I = 1, ICRD
            DO 212 J = 1, ICRD
               DO 214 K = 1, IRCRD
                  TMPMT5(I,J) = TMPMT5(I,J) + TMPMAT(I,K)*BMTRAN(K,J)
 214           CONTINUE
 212        CONTINUE
 210     CONTINUE
      END IF
      call ls_UPBFGS(MXRCRD,DELTA,GAMMA,TMPMT5,TMPMT4, &
     &     ICRD,MCRD,TMPMAT,TMPMT2,TMPMT3,optinfo%IPrint,lupri,optinfo)
!
      IF (optinfo%IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'Old Hessian')
         call ls_output(HESOLD,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)         
         call lsheader(lupri,'optinfo%BFGS-updated Hessian')
         call ls_output(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)         
         call lsheader(lupri,'Model Hessian')
         call ls_output(TMPMT4,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)         
      END IF
      IF (IFAC .GT. 20) THEN
         FAC  = 0.0E0_realk
         FAC2 = 1.0E0_realk
      ELSE
         FAC  = 1.0E0_realk/(1.0E0_realk*IFAC)
         FAC2 = 1.0E0_realk - FAC
      END IF
      DO 310 I = 1, ICRD
         DO 312 J = 1, ICRD
            HESNEW(I,J) = FAC*TMPMT4(I,J) + FAC2*HESNEW(I,J)
 312     CONTINUE
 310  CONTINUE
      IF (optinfo%IPrint .GE. IPRMIN) THEN
         call lsheader(lupri,'Updated Hessian')
         call ls_output(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      IFAC = MIN(IFAC*2,25)
      RETURN
      END

!  /* Deck upmodh */
      SUBROUTINE LS_UPMODH(Molecule,MXRCRD,DELTA,STPMAT,GAMMA,GRDMAT, &
           HESOLD,HESNEW, &
     &     MCRD,TMPVEC,TMPVC2,RMAT,GMAT,TMPMAT,TMPMT2,TMPMT3, &
     &     TMPMT4,TMPMT5,BMTRAN,MAXDIM,RESET,REJLST,GNRM, &
     &     ITYPE,lupri,optinfo)
use ls_util 
use precision
use optimization_input
use files
use molecule_type
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      TYPE(MoleculeInfo) :: Molecule
      TYPE(opt_setting) :: optinfo
      PARAMETER (D0 = 0.0E0_realk)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
!
      Real(realk) DELTA(MXRCRD), STPMAT(25,MCRD)
      Real(realk) GAMMA(MXRCRD), GRDMAT(25,MCRD)
      Real(realk) HESOLD(MXRCRD,MXRCRD), HESNEW(MCRD,MCRD)
      Real(realk) TMPVEC(MCRD), TMPVC2(MCRD)
      Real(realk) TMPMAT(MXRCRD,MXRCRD), TMPMT2(MXRCRD*MXRCRD)
      Real(realk) TMPMT3(MXRCRD,MXRCRD), TMPMT4(MXRCRD*MXRCRD)
      Real(realk) TMPMT5(MXRCRD*MXRCRD), BMTRAN(MXRCRD,MXRCRD)
      Real(realk) RMAT(MCRD,MCRD), GMAT(MCRD,MCRD)
      LOGICAL RESET, REJLST
!
!     A new model Hessian is always calculated as a start
!
      call ls_BLDHES(Molecule,MXRCRD,TMPMAT,HESOLD,lupri,optinfo)
!
!     If we use delocalized internals, we have to transform from
!     redundant to non-redundant space.
!
      IF (optinfo%DelInt) THEN
         call ls_DZERO(TMPMAT,MXRCRD*MXRCRD)
         DO 200 I = 1, MXRCRD
            DO 202 J = 1,MXRCRD 
               DO 204 K = 1, MXRCRD
                  TMPMAT(I,J) = TMPMAT(I,J) + BMTRAN(K,I)*HESOLD(K,J)
 204           CONTINUE
 202        CONTINUE
 200     CONTINUE
         call ls_DZERO(HESOLD,MXRCRD*MXRCRD)
         DO 210 I = 1, MXRCRD
            DO 212 J = 1, MXRCRD
               DO 214 K = 1, MXRCRD
                  HESOLD(I,J) = HESOLD(I,J) + TMPMAT(I,K)*BMTRAN(K,J)
 214           CONTINUE
 212        CONTINUE
 210     CONTINUE
      END IF
      IF (optinfo%IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'Model Hessian')
         call ls_output(HESOLD,1,MXRCRD,1,MXRCRD,MCRD,MCRD,1,LUPRI)         
      END IF
!
!     Then we do the updating, suppressing output.
!
      call ls_UPMULT(MXRCRD,DELTA,STPMAT,GAMMA,GRDMAT,HESOLD,HESNEW, &
     &     MXRCRD,MCRD,TMPVEC,TMPVC2,TMPMAT, &
     &     TMPMT2,TMPMT3,TMPMT4,TMPMT5,MAXDIM,RESET, &
     &     GNRM,ITYPE,.FALSE.,lupri,optinfo)
      IF (optinfo%IPrint .GE. IPRMIN) THEN
         call lsheader(lupri,'Updated Hessian')
         call ls_output(HESNEW,1,MXRCRD,1,MXRCRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

!  /* Deck upmodo */
      SUBROUTINE LS_UPMODO(Molecule,MXRCRD,DELTA,STPMAT,GAMMA,GRDMAT, &
     &     HESOLD,HESNEW, &
     &     ICRD,IRCRD,MCRD,TMPVEC,TMPVC2,RMAT,GMAT,TMPMAT,TMPMT2,TMPMT3, &
     &     BMTRAN,MAXDIM,RESET,REJLST,GNRM,ITYPE,lupri,optinfo)
use ls_util 
use precision
use optimization_input
use files
use molecule_type
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      TYPE(MoleculeInfo) :: Molecule
      TYPE(opt_setting) :: optinfo
      PARAMETER (D0 = 0.0E0_realk)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
!
      Real(realk) DELTA(MXRCRD), STPMAT(25,MCRD)
      Real(realk) GAMMA(MXRCRD), GRDMAT(25,MCRD)
      Real(realk) HESOLD(MXRCRD,MXRCRD), HESNEW(MCRD,MCRD)
      Real(realk) TMPVEC(MCRD), TMPVC2(MCRD)
      Real(realk) TMPMAT(MXRCRD,MXRCRD), TMPMT2(MXRCRD*MXRCRD)
      Real(realk) TMPMT3(MXRCRD*MXRCRD), BMTRAN(MXRCRD,MXRCRD)
      Real(realk) RMAT(MCRD,MCRD), GMAT(MCRD,MCRD)
      LOGICAL RESET, REJLST, SMART
!
!      SMART = .FALSE.
      SMART = .TRUE.
!
!     Since RMAT and GMAT are used to store a maximum of 25 vectors,
!     MCRD should be equal to or larger than 25
!
      MXDIM = MAXDIM
      IF (MCRD .LT. 25) MXDIM = MIN(MAXDIM,MCRD)
!
      call ls_DZERO(TMPVEC,MCRD)
      call ls_DZERO(TMPVC2,MCRD)
      call ls_DZERO(RMAT,MCRD*MCRD)
      call ls_DZERO(GMAT,MCRD*MCRD)
      LIMIT = 6
!
!     First we have to transfer the last step and gradient difference to
!     STPMAT and GRDMAT respectively.
!
      IF (RESET) THEN
         call ls_DZERO(GRDMAT,25*MCRD)
         call ls_DZERO(STPMAT,25*MCRD)
         INUM = 1

!
!     If the last step caused a rejected step, we discard all
!     earlier gradients.
!
!      ELSE IF (REJLST .AND. (INUM .GT. 1)) THEN
!         INUM = 1
! 141     CONTINUE
!         IF (STPMAT(INUM,1) .LT. 1.0E10_realk) THEN
!            INUM = INUM + 1
!            GOTO 141
!         END IF
!         DO 188 I = 1, ICRD
!            STPMAT(1,I) = STPMAT(INUM-1,I)
!            GRDMAT(1,I) = GRDMAT(INUM-1,I)
!            STPMAT(2,I) = 0.0E0_realk
!            GRDMAT(2,I) = 0.0E0_realk
! 188     CONTINUE
!         INUM = 2
!         STPMAT(2,1) = 1.1E10_realk
!         IF (optinfo%IPrint .GE. IPRDBG) THEN
!            WRITE(LUPRI,*)
!            WRITE(LUPRI,*) 'Removing all ' //
!     &           'displacement but one due to rejected step.'
!         END IF
!
!     After [LIMIT] iterations we discard the first entries.
!
      ELSE IF (STPMAT(LIMIT,1) .GT. 1.0E10_realk) THEN
         DO 10 I = 1, LIMIT - 1
            DO 12 J = 1, ICRD
               STPMAT(I,J) = STPMAT(I+1,J)
               GRDMAT(I,J) = GRDMAT(I+1,J)
 12         CONTINUE
 10      CONTINUE
         STPMAT(LIMIT,1) = 0.0E0_realk
         INUM = LIMIT - 1
      ELSE
         INUM = 1
 14      CONTINUE
         IF (STPMAT(INUM,1) .LT. 1.0E10_realk) THEN
            INUM = INUM + 1
            GOTO 14
         END IF
!
!     We also have to check the Real(realk) of our variational space,
!     if the number of displacement vectors exceeds this number, we
!     have to remove some.
!
         IF (SMART .AND. (INUM .GT. MXDIM)) THEN
            DO 15 I = 1, INUM-1
               DO 16 J = 1, ICRD
                  STPMAT(I,J) = STPMAT(I+1,J)
                  GRDMAT(I,J) = GRDMAT(I+1,J)
 16            CONTINUE
 15         CONTINUE
            INUM = INUM - 1
            IF (optinfo%IPrint .GE. IPRDBG) THEN
               WRITE(LUPRI,*)
               WRITE(LUPRI,*) 'Removing one displacement ' // &
     &              'due to Real(realk) of variational space.'
            END IF
         END IF
      END IF
!
!     Then we update all vectors and calculate the norm of the
!     last displacement (SNMLST).
!
      SNMLST = 0.0E0_realk
      DO 17 I = 1, ICRD
         STPMAT(INUM,I) = -DELTA(I)
         GRDMAT(INUM,I) = -GAMMA(I)
         SNMLST = SNMLST + DELTA(I)*DELTA(I)
         DO 19 II = 1, INUM-1
            STPMAT(II,I) = STPMAT(II,I)-DELTA(I)
            GRDMAT(II,I) = GRDMAT(II,I)-GAMMA(I)
 19      CONTINUE
 17   CONTINUE
      SNMLST = SQRT(SNMLST)
      STPMAT(INUM+1,1) = 1.1E10_realk
!
!     We start by calculating the normalized displacement vectors, and
!     we scale the gradient differences. Displacements larger than 100
!     times the last step, and displacements with a gradient difference
!     larger than 0.25E0_realk are marked (later to be removed).
!
      DO 20 II = 1, INUM
         DNRM = D0
         GNRM = D0
         DO 22 I = 1, ICRD
            DNRM = DNRM + STPMAT(II,I)*STPMAT(II,I)
            GNRM = GNRM + GRDMAT(II,I)*GRDMAT(II,I)
 22      CONTINUE
         DNRM = SQRT(DNRM)
         GNRM = SQRT(GNRM)
         IF (SMART .AND. (DNRM .GE. 100E0_realk*SNMLST)) THEN
            IF (optinfo%IPrint .GE. IPRDBG) WRITE(LUPRI,*) &
     &           'Removing one displacement due to distance'
            STPMAT(II,1) = 1.1E10_realk
         ELSE IF (SMART .AND. (GNRM .GE. 0.25E0_realk)) THEN
            IF (optinfo%IPrint .GE. IPRDBG) WRITE(LUPRI,*) &
     &           'Removing one displacement due to gradient norm'
            STPMAT(II,1) = 1.1E10_realk
         ELSE
            DO 24 I = 1, ICRD
!               RMAT(II,I) = STPMAT(II,I)/DNRM
!               GMAT(II,I) = GRDMAT(II,I)/DNRM
               RMAT(II,I) = STPMAT(II,I)
               GMAT(II,I) = GRDMAT(II,I)
 24         CONTINUE
         END IF
 20   CONTINUE
      II = 1
 25   CONTINUE
      IF ((STPMAT(II,1) .GE. 1.0E10_realk) .AND. (INUM .GT. 0)) THEN
         DO 27 JJ = II, INUM - 1
            DO 29 J = 1, ICRD
               STPMAT(JJ,J) = STPMAT(JJ + 1,J)
               GRDMAT(JJ,J) = GRDMAT(JJ + 1,J)
               RMAT(JJ,J)   = RMAT(JJ + 1,J)
               GMAT(JJ,J)   = GMAT(JJ + 1,J)
 29         CONTINUE
 27      CONTINUE
         INUM = INUM - 1
         GOTO 25
      ELSE IF (II .LT. INUM-1) THEN
         II = II + 1
         GOTO 25
      END IF
!
!     We check if the last displacement is nearly parallell to an
!     earlier step (dot product larger than 0.75). If that is the case,
!     the older step is removed.
!
      IF (SMART .AND. (INUM .GT. 1)) THEN
         II = 1
 30      CONTINUE
         DOTP = D0
         DO 32 I = 1, ICRD
            DOTP = DOTP + RMAT(II,I)*RMAT(INUM,I)
 32      CONTINUE
         IF (DOTP .GE. 0.75E0_realk) THEN
            DO 34 J = 1, ICRD
               DO 36 JJ = II, INUM - 1
                  RMAT(JJ,J)   = RMAT(JJ+1,J)
                  STPMAT(JJ,J) = STPMAT(JJ+1,J)
                  GMAT(JJ,J)   = RMAT(JJ+1,J)
                  GRDMAT(JJ,J) = STPMAT(JJ+1,J)
 36            CONTINUE
 34         CONTINUE
            STPMAT(INUM,1) = 1.1E10_realk
            INUM = INUM - 1
            IF (optinfo%IPrint .GE. IPRDBG) WRITE(LUPRI,*) &
     &           'Removing one nearly parallell and older displacement'
            IF (II .LT. INUM) GOTO 30
         ELSE IF (II .LT. INUM-1) THEN
            II = II + 1
            GOTO 30
         END IF
      END IF
!
!     Some output
!
      IF (optinfo%IPrint .GE. IPRMAX) THEN
         call lsheader(lupri,'Original displacements')
         call ls_output(STPMAT,1,INUM,1,ICRD,25,MXRCRD,1,LUPRI)
         call lsheader(lupri,'Normalized displacements')
         call ls_output(RMAT,1,INUM,1,ICRD,MCRD,MCRD,1,LUPRI)
         call lsheader(lupri,'Original gradient vectors')
         call ls_output(GRDMAT,1,INUM,1,ICRD,25,MXRCRD,1,LUPRI)
         call lsheader(lupri,'Scaled gradient vectors')
         call ls_output(GMAT,1,INUM,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
!
!     A new model Hessian is always calculated as a start
!
      call ls_BLDHES(Molecule,MXRCRD,TMPMAT,HESNEW,lupri,optinfo)
!
!     If we use delocalized internals, we have to transform from
!     redundant to non-redundant space.
!
      IF (optinfo%DelInt) THEN
         call ls_DZERO(TMPMAT,MXRCRD*MXRCRD)
         DO 200 I = 1, ICRD
            DO 202 J = 1, IRCRD
               DO 204 K = 1, IRCRD
                  TMPMAT(I,J) = TMPMAT(I,J) + BMTRAN(K,I)*HESNEW(K,J)
 204           CONTINUE
 202        CONTINUE
 200     CONTINUE
         call ls_DZERO(HESNEW,MXRCRD*MXRCRD)
         DO 210 I = 1, ICRD
            DO 212 J = 1, ICRD
               DO 214 K = 1, IRCRD
                  HESNEW(I,J) = HESNEW(I,J) + TMPMAT(I,K)*BMTRAN(K,J)
 214           CONTINUE
 212        CONTINUE
 210     CONTINUE
      END IF
      IF (optinfo%IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'Model Hessian')
         call ls_output(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)         
      END IF
!
!     Then we use the required subroutine to do the updating,
!     suppressing output.
!
      DO 60 K = 1, INUM
         call ls_DZERO(DELTA,MXRCRD)
         call ls_DZERO(GAMMA,MXRCRD)
         DO 62 I = 1, ICRD
             DELTA(I) = -RMAT(K,I)
             GAMMA(I) = -GMAT(K,I)
 62      CONTINUE
         IF (ITYPE .EQ. 1) THEN
            call ls_UPBFGS(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW, &
     &           ICRD,MCRD,TMPVEC,TMPMAT,TMPMT2,-1,lupri,optinfo)
         ELSE IF (ITYPE .EQ. 2) THEN
            call LS_UPDFP(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW, &
     &           ICRD,MCRD,TMPVEC,TMPMAT,TMPMT2,-1,lupri,optinfo)
         ELSE IF (ITYPE .EQ. 3) THEN
            call LS_UPPSB(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW, &
     &           ICRD,MCRD,TMPVEC,TMPMAT,TMPMT2,-1,lupri,optinfo)
         ELSE IF (ITYPE .EQ. 4) THEN
            call ls_UPRNKO(MXRCRD,DELTA,GAMMA,HESOLD,HESNEW, &
     &           ICRD,MCRD,TMPVEC,TMPMAT,TMPMT2,-1,lupri,optinfo)
         END IF
         IF ((optinfo%IPrint .GE. IPRDBG) .AND. (K .LT. INUM)) THEN
            call lsheader(lupri,'Updating the Hessian')
            call ls_output(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)         
         END IF
         call ls_DZERO(HESOLD,MXRCRD*MXRCRD)
         DO 65 J = 1, ICRD
            DO 67 I = 1, ICRD
               HESOLD(I,J) = HESNEW(I,J)
 67         CONTINUE
 65      CONTINUE
 60   CONTINUE
      IF (optinfo%IPrint .GE. IPRMIN) THEN
         call lsheader(lupri,'Updated Hessian')
         call ls_output(HESNEW,1,ICRD,1,ICRD,MCRD,MCRD,1,LUPRI)
      END IF
      RETURN
      END

!  /* Deck gttmat */
!
!     Construct, scale and orthogonalize T matrix.
!
      SUBROUTINE LS_GTTMAT(TMAT,TMPMAT,NCR,NPR,THRLDP,lupri,optinfo)
use ls_util 
use optimization_input
use files
use precision
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri,NCR
      Real(realk) :: Coord(NCR)
      TYPE(opt_setting) :: optinfo
      Real(realk) TMAT(NCR,NPR), TMPMAT(MXCOOR)
!     Forming Cord array
      Do i = 1, NCR/3
         Coord(3*i-2:3*i) = optinfo%Coordinates(:,i)
      Enddo
! 
      call ls_DZERO(TMAT,NCR*NPR)
      call ls_DZERO(TMPMAT,MXCOOR)
!     Getting the inverse of the T matrix
      Call Calc_TMatrix(TMat,NCR,Coord)
      IF (optinfo%IPrint .GT. 5) THEN
         call lsheader(lupri,'T matrix in GTTMAT')
         call ls_output(TMAT,1,NCR,1,NPR,NCR,NPR,1,LUPRI)
      END IF
      NPR1 = NPR
      call ls_ORTVEC(0,NPR1,NCR,THRLDP,TMAT,lupri)
      IF (optinfo%IPrint .GT. 5) THEN
         call lsheader(lupri,'Orthogonalized T matrix in GTTMAT')
         call ls_output(TMAT,1,NCR,1,NPR,NCR,NPR,1,LUPRI)
      END IF
      IF (NPR1 .NE. NPR) THEN
         WRITE(LUPRI,'(//,2(A,I1),A,/A)')                               &
     &      ' Number of trarot vectors reduced from ', NPR,             &
     &      ' to ', NPR1, ' in ORTVEC called from WLKPRJ.',             &
     &      ' Program cannot proceed .'
          If (NPR1 .EQ. 5) then
              WRITE(LUPRI,*)                          &
     &      ' Currently DALTON can not deal with linear molecules'
          Endif
         call lsquit('Insufficient number of trarot vectors in WLKPRJ',lupri)
      END IF
      RETURN
      END

!  /* Deck projgh */
      SUBROUTINE LS_PROJGH(EGRAD,EHESS,ALLHES,TMAT,TMPMAT,TMPMT2,PROJOP,NCRTOT, &
      & NTempMat,lupri,optinfo)
use ls_util 
use precision
use optimization_input
use files
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Integer :: NCRTOT
      Integer :: NTempMat
      TYPE(opt_setting) :: optinfo
      Real(realk) EGRAD(MXCOOR), EHESS(MXCOOR,MXCOOR)
      Real(realk) ALLHES(NCRTOT*NCRTOT), TMAT(NTempMat)
      Real(realk) TMPMAT(MXCOOR,MXCOOR), TMPMT2(MXCOOR)
      Real(realk) PROJOP(MXCOOR,MXCOOR)
      Real(realk), parameter :: THRLDP = 1.0E-28_realk
      IF (optinfo%IPrint .GT. 5) call lsheader(lupri,'Output from PROJGH')
! We assume that molecule is non-linear, number of coordinates to be projected out is 6
      NPR = 6
      call ls_PRJGH1(EGRAD,ALLHES,TMAT,PROJOP, &
     &      TMPMAT,TMPMT2,NCRTOT,NPR,THRLDP,lupri,NCRTOT,optinfo)
      JI = 1
      DO 20 I = 1, optinfo%ICartCoord
         DO 30 J = 1, optinfo%ICartCoord
            EHESS(J,I) = ALLHES(JI)
            JI = JI + 1
 30      CONTINUE
 20   CONTINUE            
      RETURN
      END

!  /* Deck prjgh1 */
!
!     Construct projection operator and use it on gradient and
!     Hessian to remove both rotation and translation.
!
      SUBROUTINE LS_PRJGH1(EGRAD,ALLHES,TMAT,PROJOP,TMPMAT,TMPMT2, &
     &                     NCR,NPR,THRLDP,lupri,NCRTOT,optinfo)
use ls_util 
use precision
use optimization_input
use files
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      TYPE(opt_setting) :: optinfo
      Real(realk) EGRAD(NCR), ALLHES(NCR,NCR)
      Real(realk) TMAT(NCR,NPR), PROJOP(NCR,NCR)
      Real(realk) TMPMAT(MXCOOR,MXCOOR), TMPMT2(MXCOOR)
      call ls_DZERO(PROJOP,NCR*NCR)
      call ls_DZERO(TMPMAT,MXCOOR*MXCOOR)
!
!     Get translation and rotation matrix.
!
      call ls_GTTMAT(TMAT,TMPMT2,NCR,NPR,THRLDP,lupri,optinfo)
!
!     Construct operator.
!
      call ls_DUNIT(PROJOP,NCR)
      call dgemm('N','T',NCR,NCR,NPR,-1E0_realk,TMAT,NCR, &
     &           TMAT,NCR,1E0_realk,PROJOP,NCR)
      IF (optinfo%IPrint .GT. 5) THEN
         call lsheader(lupri,'Unprojected gradient in PROJGH')
         call ls_output(EGRAD,1,1,1,NCR,1,NCR,1,LUPRI)
         call lsheader(lupri,'Unprojected Hessian in PROJGH')
         call ls_output(ALLHES,1,NCR,1,NCR,NCR,NCR,1,LUPRI)
         call lsheader(lupri,'Projection operator in PROJGH')
         call ls_output(PROJOP,1,NCR,1,NCR,NCR,NCR,1,LUPRI)
      END IF
!
!     Do projection.
!
      call dgemm('N','N',NCR,1,NCR,1E0_realk, &
     &           PROJOP,NCR, &
     &           EGRAD,NCR,0E0_realk, &
     &           TMPMAT,NCR)
      call dcopy(NCR,TMPMAT,1,EGRAD,1)
      IF (optinfo%IPrint .GT. 5) THEN
         call lsheader(lupri,'Projected gradient in PROJGH')
         call ls_output(EGRAD,1,1,1,NCR,1,NCR,1,LUPRI)
      END IF
      call dgemm('N','N',NCR,NCR,NCR,1E0_realk, &
     &           PROJOP,NCR, &
     &           ALLHES,NCR,0E0_realk, &
     &           TMPMAT,NCR)
      call dgemm('N','N',NCR,NCR,NCR,1E0_realk, &
     &           TMPMAT,NCR, &
     &           PROJOP,NCR,0E0_realk, &
     &           ALLHES,NCR)
      IF (optinfo%IPrint .GT. 5) THEN
         call lsheader(lupri,'Projected Hessian in PROJGH')
         call ls_output(ALLHES,1,NCR,1,NCR,NCR,NCR,1,LUPRI)
      END IF
      RETURN
      END

!  /* Deck diahes */
!
!     Diagonalize Hessian
!
      SUBROUTINE LS_DIAHES(MXRCRD,MX2CRD,NCRDHS,EGRAD,EHESS,ALLHES, &
     & TMAT,THRIND,EVEC,EVCTMP,TMPHES,HESPCK,NCRTOT,NTempMat, &
     & lupri,optinfo)
use ls_util 
use precision
use files
use optimization_input
Implicit Real(realk) (A-H,O-Z)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      PARAMETER (D0 = 0.0E0_realk)
      Integer :: lupri
      Integer :: NCRTOT
      Integer :: NTempMAt
      Integer :: NZEROG
! 1  - arbitrary selected, since not used,change when used
      TYPE(opt_setting) :: optinfo
      Real(realk) EGRAD(MXCOOR), EHESS(NCRTOT,NCRTOT)
      Real(realk) ALLHES(MXCOOR*MXCOOR), TMAT(NTempMat)
      Real(realk) EVEC(MX2CRD,MX2CRD), EVCTMP(MX2CRD*MX2CRD)
      Real(realk) TMPHES(MX2CRD,MX2CRD)
      Real(realk) HESPCK(NCRDHS*NCRDHS)
     
      LOGICAL PRJTRO
      PRJTRO = .TRUE.
!
      call ls_DZERO(optinfo%GRDDIA,MXRCRD)
!     Vladimir:
!     NZEROG is number of gradient elements set to zero, feature in WALK
!     Not implemented here, thus, set to 0
      NZEROG = 0
!
      call ls_WLKEIG(EGRAD,ALLHES,optinfo%EVAL,EVCTMP,optinfo%GRDDIA,TMAT,THRLDP, &
     &            THRIND,optinfo%CNDHES,optinfo%INDHES, &
     &   NCRTOT,NCRTOT*NCRTOT,NTempMat,PRJTRO,optinfo%IPrint,lupri)
      JI = 1
!
!     Diagonal hessian and eigenvalues are copied.
!
      DO 10 I = 1, optinfo%ICartCoord
         DO 20 J = 1, optinfo%ICartCoord
            EHESS(J,I) = ALLHES(JI)
            JI = JI + 1
 20      CONTINUE
 10   CONTINUE
      JI = 1
      II = 0
      DO 40 I = 1, NCRTOT
         DO 50 J = 1, NCRTOT
            EVEC(II+J,II+I) = EVCTMP(JI)
            JI = JI + 1
 50      CONTINUE
 40   CONTINUE
!
!     If we're using the rational function method, we modifiy the
!     totally symmetric part of the Hessian.
!
      IF (optinfo%RatFun .AND. (.NOT. optinfo%Saddle)) THEN
!
!     Eigenvalues and -vectors must be shifted to make room for one more.
!
         call ls_DZERO(TMPHES,MX2CRD*MX2CRD)
         call ls_DZERO(EVCTMP,NCRDHS*NCRDHS)
         DO 70 J = 1,  NCRTOT
            DO 72 I = 1,  NCRTOT
               EVCTMP(I+(J-1)*NCRDHS) = optinfo%EVAL(I)*EVEC(J,I)
 72         CONTINUE
 70      CONTINUE
         DO 75 I = 1,  NCRTOT
            DO 77 J = 1,  NCRTOT
               DO 79 K = 1,  NCRTOT
                  TMPHES(I,J) = TMPHES(I,J) &
     &                 + EVEC(I,K)*EVCTMP(K+(J-1)*NCRDHS)
 79            CONTINUE
 77         CONTINUE
 75      CONTINUE
         DO 80 I = 1, NCRTOT
            TMPHES(NCRDHS,I) = EGRAD(I)
            TMPHES(I,NCRDHS) = EGRAD(I)
 80      CONTINUE
         IF (optinfo%IPrint .GE. IPRDBG) THEN
            call lsheader(lupri,'Augmented Hessian')
            call ls_output(TMPHES,1,NCRDHS,1,NCRDHS,MX2CRD, &
     &           MX2CRD,1,LUPRI)
         END IF
         TMPHES(NCRDHS,NCRDHS) = D0
         call ls_DZERO(EVCTMP,NCRDHS*NCRDHS)
         DO 83 I = 1, NCRDHS
            DO 85 J = 1, NCRDHS
               EVCTMP(J+(I-1)*NCRDHS) = TMPHES(I,J)
 85         CONTINUE
 83      CONTINUE
         call ls_DZERO(HESPCK,NCRDHS*NCRDHS)
         call ls_DSITSP(NCRDHS,EVCTMP,HESPCK)
         call ls_DUNIT(EVCTMP,NCRDHS)
         call ls_JACO(HESPCK,EVCTMP,NCRDHS,NCRDHS,NCRDHS, &
     &        TMPHES(1,1),TMPHES(1,2))
         DO 90 J = 1, NCRDHS
            optinfo%EVAL(J) = HESPCK(J*(J+1)/2)
            optinfo%GRDDIA(J) = DDOT(NCRDHS,EGRAD,1,EVCTMP(1+(J-1)*NCRDHS),1)
            DO 92 I = 1, NCRDHS
               EVEC(I,J) = EVCTMP(I+(J-1)*NCRDHS)
 92         CONTINUE
 90      CONTINUE
         IF ((optinfo%EVAL(1) .LT. -THRIND) .AND. optinfo%IndHes .GT. 0) &
     &        optinfo%IndHes = optinfo%IndHes - 1
         DO 95 I = 1, NCRDHS
            IF (ABS(optinfo%EVAL(I)) .LE. 1.0E-6_realk) optinfo%EVAL(I) = optinfo%EVAL(I) + 1.0E5_realk
 95      CONTINUE
!
!     The eigenvalues are sorted
!
         DO 100 I = 1, NCRDHS
            JMIN = I
            EMIN = optinfo%EVAL(I)
            DO 105 J = (I + 1), NCRDHS
               IF (optinfo%EVAL(J) .LT. EMIN) THEN
                  EMIN = optinfo%EVAL(J)
                  JMIN = J
               END IF
 105        CONTINUE
            IF (JMIN .NE. I) THEN
               call dswap(1,  optinfo%EVAL  (I),1,optinfo%EVAL  (JMIN),1)
               call dswap(MX2CRD,EVEC(1,I),1,EVEC(1,JMIN),1)
               call dswap(1,optinfo%GRDDIA(I),1,optinfo%GRDDIA(JMIN),1)
            END IF
 100     CONTINUE
         DO 120 I = 1, NCRDHS
            IF (ABS(ABS(optinfo%EVAL(I))-1.0E5_realk) .LT. 1.0E-3_realk) &
     &           optinfo%EVAL(I) = optinfo%EVAL(I) - 1.0E5_realk
 120     CONTINUE
         IF (optinfo%IPrint .GE. IPRDBG) THEN
            WRITE(LUPRI,*) 'Index of Hessian: ',optinfo%IndHes
            call lsheader(lupri,'RF-eigenvalues')
            call ls_output(optinfo%EVAL,1,1,1,NCRDHS,1,MXRCRD,1,LUPRI)
            call lsheader(lupri,'RF-eigenvectors')
            call ls_output(EVEC,1,NCRDHS,1,NCRDHS,MX2CRD, &
     &           MX2CRD,1,LUPRI)
         END IF
      END IF
      RETURN
      END

!  /* Deck uptrad */
      SUBROUTINE LS_UPTRAD(REJGEO,lupri,optinfo)
use ls_util 
use precision
use optimization_input
use files
!
!     Updates the trust radius and checks if step should be rejected.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      TYPE(opt_setting) :: optinfo
      PARAMETER (DP25=0.25E0_realk)
      LOGICAL REJGEO
      THDE = 1.0E-10_realk
!!!!  Vladimir: desactivated until we have CC in LSDALTON,then should be rewritten
!
!     CC energies are less precise and we have to decrease the
!     energy threshold
!
!      IF (DOCCSD) THDE = 1.0E-8_realk
!
      ERGDIF = optinfo%energy - optinfo%energyOld
      REJGEO = .FALSE.
      optinfo%Rebid = .FALSE.
!filip - print always info about the actual and predicted energy (according to Vebjorn's
!        suggestion). 07.11.2008.
      call lsheader(lupri,'StepInfo:Energy difference to previous geometry:')
      write(lupri,'(/a,/)') &
     &'      Energy       Old_Energy     Actual_diff  Predicted_diff' &
     & //'  Ratio'
      write(lupri,'(2F20.10,3F15.10)') &
     & optinfo%energy, optinfo%energyOld, ERGDIF, optinfo%predictedChange, ERGDIF/optinfo%predictedChange
      write(lupri,*)
!filip - end.
!
!     The ratio between actual and predicted energy change is calculated,
!     and this ratio is then used to update the trust radius.
!
      IF (ABS(optinfo%predictedChange) .GT. THDE) THEN
         RATIO = ERGDIF/optinfo%predictedChange 
         IF (optinfo%IPrint .GT. 2) THEN
            call lsheader(lupri,'Energy difference to previous geometry:')
            WRITE(LUPRI,'(/A/5X,F15.10,2(A,F15.10))') &
     &      '      Actual           /  Predicted       =    Ratio ', &
     &            ERGDIF,'  /  ',optinfo%predictedChange,'  =  ',RATIO
         END IF
      ELSE
         RATIO = 1.0E0_realk
         IF (optinfo%IPrint .GT. 2) THEN
            WRITE(LUPRI,'(3(/A),/5X,1P,2D16.6)') &
     &          ' Close to convergence, ratio set to one.', &
     &          ' Energy difference to previous geometry:', &
     &          ' actual and predicted:', ERGDIF, optinfo%predictedChange
         END IF
      END IF
      IF ((RATIO .LE. optinfo%RTRJMN) .OR. (RATIO .GE. optinfo%RTRJMX)) THEN
         IF (optinfo%Saddle) THEN
            WRITE(LUPRI,*) &
     &           'Trust radius squarely decreased due to bad ratio.'
            optinfo%TrustRad = optinfo%TrustDe*optinfo%TrustDe*optinfo%StepNorm
         ELSE
            WRITE(LUPRI,'(/A/A,F12.6)') &
     &           ' Step rejected because ratio between predicted', &
     &           ' and actual energy change is : ', RATIO
            REJGEO = .TRUE.
         END IF
      ELSE IF (RATIO .LT. optinfo%RTENBD) THEN
         WRITE(LUPRI,*) 'Trust radius decreased due to bad ratio.'
         optinfo%TrustRad = optinfo%TrustDe*optinfo%StepNorm
      ELSE IF (ABS(RATIO-1.0E0_realk) .LE. DP25*(1.0E0_realk-optinfo%RTENGD)) THEN
         WRITE(LUPRI,*) &
     &       'Trust radius squarely increased due to very good ratio.'
         optinfo%TrustRad = MAX(optinfo%TrustIn*optinfo%TrustIn*optinfo%StepNorm,optinfo%TrustRad)
!hjaug99 IF (ISTATE .GT. 1 .AND. optinfo%TrustRad .GT. 0.30E0_realk) optinfo%TrustRad = 0.30E0_realk
      ELSE IF (RATIO .GE. optinfo%RTENGD) THEN
         WRITE(LUPRI,*) 'Trust radius increased due to good ratio.'
         optinfo%TrustRad = MAX(optinfo%TrustIn*optinfo%StepNorm,optinfo%TrustRad)
!hjaug99 IF (ISTATE .GT. 1 .AND. optinfo%TrustRad .GT. 0.30E0_realk) optinfo%TrustRad = 0.30E0_realk
      ELSE
         WRITE(LUPRI,*) 'Trust radius set equal to norm of step.'
         optinfo%TrustRad = optinfo%StepNorm
!hjaug99 IF (ISTATE .GT. 1 .AND. optinfo%TrustRad .GT. 0.30E0_realk) optinfo%TrustRad = 0.30E0_realk
      END IF
!
!     For saddle point optimizations we place both an upper and
!     lower bound on the trust radius.
!
      IF (optinfo%Saddle) optinfo%TrustRad = MAX(0.025E0_realk, MIN(1.0E0_realk,optinfo%TrustRad))
!      IF (optinfo%Saddle .AND. optinfo%DelInt .AND.
!     &     ((RATIO .LE. optinfo%RTRJMN) .OR. (RATIO .GE. optinfo%RTRJMX)))
!     &     optinfo%Rebid = .TRUE.
      IF (optinfo%DelInt .AND. optinfo%Newton .AND. optinfo%Saddle) optinfo%Rebid = .TRUE.
      IF (.NOT. REJGEO) &
     &     WRITE(LUPRI,'(A,F10.5)') ' Updated trust radius', optinfo%TrustRad
      RETURN
      END

!  /* Deck rfstp */
      SUBROUTINE LS_RFSTP(MX2CRD,NCHESS,NCRD,ICRD,EVEC,STEP,GRAD, &
     &     TMPMAT,HESSMT,lupri,optinfo)
use ls_util 
use precision
use optimization_input
use files
use Fundamental
!
!     Determines a rational function (RF) step.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      TYPE(opt_setting) :: optinfo
      Real(realk) EVEC(MX2CRD,MX2CRD)
      Real(realk) STEP(NCRD), GRAD(NCRD)
      Real(realk) TMPMAT(MX2CRD*MX2CRD), HESSMT(NCRD,NCRD)
      LOGICAL STPSCL
      PARAMETER (D0 = 0.0E0_realk, D1 = 1.0E0_realk, DP5 = 0.5E0_realk)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      STPSCL = .TRUE.
!      STPSCL = .FALSE.
      IF (optinfo%IPrint .GE. IPRDBG) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) &
     &        'Lowest eigenvalue (level-shift parameter): ',optinfo%EVAL(1)
         call lsheader(lupri,'Corresponding eigenvector')
         call ls_output(EVEC,1,1,1,NCHESS,1,MX2CRD,1,LUPRI)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Scaling factor: ',EVEC(NCHESS,1)
      END IF
      FAC = EVEC(NCHESS,1)
      IF (ABS(FAC) .GT. 1.0E-8_realk) THEN
         FAC = 1.0E0_realk/FAC
      ELSE
         FAC = 1.0E8_realk
      END IF
      DO 30 I = 1, ICRD
         STEP(I) = EVEC(I,1)*FAC
!
!     For angles and dihedral angles we have to avoid step components
!     giving multiples of 2*pi.
!
         IF (optinfo%RedInt .AND. optinfo%INTCRD(I,1) .GT. 10) &
     &        STEP(I) = MOD(STEP(I),2.0E0_realk*PI)
!
!     If the step is too large, we simply restrict each element
!     to be below the trust radius.
!
         IF ((ABS(STEP(I)) .GT. optinfo%TrustRad) .AND. (.NOT. STPSCL)) &
     &        STEP(I) = SIGN(optinfo%TrustRad,STEP(I))
 30   CONTINUE
      optinfo%StepNorm = SQRT(DDOT(ICRD,STEP,1,STEP,1))
      IF (optinfo%IPrint .GE. IPRDBG) THEN
         WRITE(LUPRI,'(/A,1P,D10.2/)') &
     &      'RF-Step length:', optinfo%StepNorm
         IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
            call lsheader(lupri,'RF-Step in internal coordinates')
         ELSE
            call lsheader(lupri,'RF-step')
         END IF
         call ls_output(STEP,1,1,1,ICRD,1,NCRD,1,LUPRI)
      END IF
!
!     Alternatively we restrivt the step norm to be equal or less
!     than the trust radius.
!
      IF ((optinfo%StepNorm .GT. optinfo%TrustRad) .AND. STPSCL) THEN
         FAC = optinfo%TrustRad/optinfo%StepNorm
         DO 32 I = 1, ICRD
            STEP(I) = STEP(I)*FAC
 32      CONTINUE
         IF (optinfo%IPrint .GE. IPRDBG) THEN
            WRITE(LUPRI,'(/A,1P,D10.2/)') &
     &           'Step too long, step scaled by factor:', FAC
            call lsheader(lupri,'Scaled RF-Step')
            call ls_output(STEP,1,1,1,ICRD,1,NCRD,1,LUPRI)
         END IF
         optinfo%StepNorm = SQRT(DDOT(ICRD,STEP,1,STEP,1))
      END IF
      ICNT = 1
      FAC = 1.0E0_realk
 37   CONTINUE
      IF (ICNT .LE. 10) THEN
         call ls_DZERO(TMPMAT,MX2CRD)
         SNDTRM = 0.0E0_realk
         DO 40 I = 1, ICRD
            DO 45 J = 1, ICRD
               TMPMAT(I) = TMPMAT(I) + HESSMT(I,J)*STEP(J)
 45         CONTINUE
            SNDTRM = SNDTRM + TMPMAT(I)*STEP(I)
 40      CONTINUE
         optinfo%predictedChange = DDOT(ICRD,GRAD,1,STEP,1) &
     &        + 0.5E0_realk*SNDTRM
         IF (optinfo%predictedChange .GT. 0.0E0_realk) THEN
            DO 50 I = 1, ICRD
               STEP(I) = EVEC(I,1)*FAC
 50         CONTINUE
            FAC = 0.5E0_realk*FAC
            ICNT = ICNT + 1
            GOTO 37
         END IF
      ELSE
         IF (optinfo%IPrint .GE. IPRMIN) THEN
            IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
               call lsheader(lupri,'RF-Step in internal coordinates')
            ELSE
               call lsheader(lupri,'RF-step')
            END IF
            call ls_output(STEP,1,1,1,ICRD,1,NCRD,1,LUPRI)
         END IF
         DO 55 I = 1, ICRD
            STEP(I) = -GRAD(I)
 55      CONTINUE
         call ls_DZERO(TMPMAT,MX2CRD)
         SNDTRM = 0.0E0_realk
         DO 60 I = 1, ICRD
            DO 65 J = 1, ICRD
               TMPMAT(I) = TMPMAT(I) + HESSMT(I,J)*STEP(J)
 65         CONTINUE
            SNDTRM = SNDTRM + TMPMAT(I)*STEP(I)
 60      CONTINUE
         optinfo%predictedChange = DDOT(ICRD,GRAD,1,STEP,1) &
     &        + 0.5E0_realk*SNDTRM
      END IF
      IF ((ICNT .GT. 1) .AND. (optinfo%IPrint .GE. IPRMIN)) THEN
         IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
            call lsheader(lupri,'RF-Step in internal coordinates')
         ELSE
            call lsheader(lupri,'RF-step')
         END IF
         call ls_output(STEP,1,1,1,ICRD,1,NCRD,1,LUPRI)
      END IF
      RETURN
      END

!  /* Deck prfstp */
      SUBROUTINE LS_PRFSTP(MX2CRD,NCHESS,NCRD,EVEC,STEP,GRAD, &
     &     TMPMAT,HESSMT,IPRF,lupri,optinfo)
use ls_util 
use precision
use optimization_input
use files
use Fundamental
!
!     Determines a partitioned rational function (RF) step.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      TYPE(opt_setting) :: optinfo
      Real(realk) EVEC(MX2CRD,MX2CRD)
      Real(realk) STEP(NCRD), GRAD(NCRD)
      Real(realk) TMPMAT(MX2CRD*MX2CRD), HESSMT(NCRD,NCRD)
      LOGICAL STPSCL
      PARAMETER (D0 = 0.0E0_realk, D1 = 1.0E0_realk, DP5 = 0.5E0_realk)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
!
      STPSCL = .TRUE.
!      STPSCL = .FALSE.
!
!     This subroutine is called three times. First to minimize one
!     partition, then to maximize another partition, and finally to
!     scale the step and predict the energy change. IPRF keeps track
!     of this (IPRF = 1,2,3).
!
      IF (IPRF .EQ. 3) GOTO 277
!
      IMOD = 1
      IF (IPRF .EQ. 2) IMOD = NCHESS
      IF (optinfo%IPrint .GE. IPRDBG) THEN
         WRITE(LUPRI,*)
         IF (IPRF .NE. 2) THEN
            WRITE(LUPRI,*) &
     &           'Lowest eigenvalue (level-shift parameter): ',optinfo%EVAL(1)
         ELSE
            WRITE(LUPRI,*) &
     &           'Highest eigenvalue (level-shift parameter): ', &
     &           optinfo%EVAL(NCHESS)
         END IF
         call lsheader(lupri,'Corresponding eigenvector')
         call ls_output(EVEC,1,NCHESS,IMOD,IMOD,MX2CRD,MX2CRD,1,LUPRI)
      END IF
      FAC = EVEC(NCHESS,IMOD)
      IF (ABS(FAC) .GT. 1.0E-10_realk) THEN
         FAC = 1.0E0_realk/FAC
      ELSE IF ((IMOD .EQ. 1) .AND. (optinfo%EVAL(IMOD+1) .LT. D0) &
     &        .AND. (optinfo%GradNorm .GT. 1.0E-4_realk)) THEN
 10      CONTINUE
         IMOD = IMOD + 1
         IF (optinfo%EVAL(IMOD) .LT. D0) THEN
            IF (ABS(EVEC(NCHESS,IMOD)) .GT. 1.0E-10_realk) THEN
               FAC = 1.0E0_realk/EVEC(NCHESS,IMOD)
            ELSE
               GOTO 10
            END IF
         ELSE
            IMOD = 1
            FAC = 1.0E10_realk
         END IF
      ELSE
         FAC = 1.0E10_realk
      END IF
      IF (optinfo%IPrint .GE. IPRDBG) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Scaling factor: ',FAC
      END IF
!
!     We only need the preliminary RF-step the first two times the
!     subroutine is called, and we return at this point.
!     We also check the factor in use.
!
      IF ((IPRF .EQ. 1) .OR. (IPRF .EQ. 2)) THEN
         IF ((ABS(FAC) .GT. 1.0E4_realk) .OR. (ABS(FAC) .LT. 1.0E-4_realk)) FAC = D0
         DO 25 I = 1, NCHESS-1
            STEP(I) = EVEC(I,IMOD)*FAC
 25      CONTINUE
         RETURN
      END IF
!
!     The second part of the subroutine scales the step and
!     predicts the energy.
!
 277  CONTINUE
!
      DO 30 I = 1, NCRD
!
!     For angles and dihedral angles we have to avoid step components
!     giving multiples of 2*pi.
!
         IF (optinfo%RedInt .AND. optinfo%INTCRD(I,1) .GT. 10) &
     &        STEP(I) = MOD(STEP(I),2.0E0_realk*PI)
!
!     If the step is too large, we simply restrict each element
!     to be below the trust radius.
!
         IF ((ABS(STEP(I)) .GT. optinfo%TrustRad) .AND. (.NOT. STPSCL)) &
     &        STEP(I) = SIGN(optinfo%TrustRad,STEP(I))
 30   CONTINUE
      optinfo%StepNorm = SQRT(DDOT(NCRD,STEP,1,STEP,1))
      IF (optinfo%IPrint .GE. IPRDBG) THEN
         WRITE(LUPRI,'(/A,1P,D10.2/)') &
     &      'RF-Step length:', optinfo%StepNorm
         IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
            call lsheader(lupri,'RF-step in internal coordinates')
         ELSE
            call lsheader(lupri,'RF-step')
         END IF
         call ls_output(STEP,1,1,1,NCRD,1,NCRD,1,LUPRI)
      END IF
!
!     Alternatively we restrivt the step norm to be equal or less
!     than the trust radius.
!
      IF ((optinfo%StepNorm .GT. optinfo%TrustRad) .AND. STPSCL) THEN
         FAC = optinfo%TrustRad/optinfo%StepNorm
         DO 32 I = 1, NCRD
            STEP(I) = STEP(I)*FAC
 32      CONTINUE
         IF (optinfo%IPrint .GE. IPRDBG) THEN
            WRITE(LUPRI,'(/A,1P,D10.2/)') &
     &           'Step too long, step scaled by factor:', FAC
            call lsheader(lupri,'Scaled RF-Step')
            call ls_output(STEP,1,1,1,NCRD,1,NCRD,1,LUPRI)
         END IF
      END IF
      ICNT = 1
      FAC = 1.0E0_realk
 37   CONTINUE
      IF (ICNT .LE. 10) THEN
         call ls_DZERO(TMPMAT,MX2CRD)
         SNDTRM = 0.0E0_realk
         DO 40 I = 1, NCRD
            DO 45 J = 1, NCRD
               TMPMAT(I) = TMPMAT(I) + HESSMT(I,J)*STEP(J)
 45         CONTINUE
            SNDTRM = SNDTRM + TMPMAT(I)*STEP(I)
 40      CONTINUE
         optinfo%predictedChange = DDOT(NCRD,GRAD,1,STEP,1) &
     &        + 0.5E0_realk*SNDTRM
      ELSE
         IF (optinfo%IPrint .GE. IPRMIN) THEN
            IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
               call lsheader(lupri,'RF-Step in internal coordinates')
            ELSE
               call lsheader(lupri,'RF-step')
            END IF
            call ls_output(STEP,1,1,1,NCRD,1,NCRD,1,LUPRI)
         END IF
         DO 55 I = 1, NCRD
            STEP(I) = -GRAD(I)
 55      CONTINUE
         call ls_DZERO(TMPMAT,MX2CRD)
         SNDTRM = 0.0E0_realk
         DO 60 I = 1, NCRD
            DO 65 J = 1, NCRD
               TMPMAT(I) = TMPMAT(I) + HESSMT(I,J)*STEP(J)
 65         CONTINUE
            SNDTRM = SNDTRM + TMPMAT(I)*STEP(I)
 60      CONTINUE
         optinfo%predictedChange = DDOT(NCRD,GRAD,1,STEP,1) &
     &        + 0.5E0_realk*SNDTRM
      END IF
      IF ((ICNT .GT. 1) .AND. (optinfo%IPrint .GE. IPRMIN)) THEN
         IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
            call lsheader(lupri,'RF-Step in internal coordinates')
         ELSE
            call lsheader(lupri,'RF-step')
         END IF
         call ls_output(STEP,1,1,1,NCRD,1,NCRD,1,LUPRI)
      END IF
      RETURN
      END

!  /* Deck prfstc */
      SUBROUTINE LS_PRFSTC(MXRCRD,MX2CRD,NCRDHS,EGRAD,EHESS,EVEC,EVCTMP, &
     &     TMPMAT,TMPMT2,TMPMT3,TMPMT4,VECMOD,lupri,optinfo)
use ls_util 
use precision
use optimization_input
use files
!
!     Controls saddle point optimization in Cartesian
!     coordinates using the partitioned rational function approach.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      TYPE(opt_setting) :: optinfo
      Real(realk) EGRAD(MXCOOR), EHESS(MXCOOR,MXCOOR)
      Real(realk) EVEC(MX2CRD,MX2CRD), EVCTMP(NCRDHS*NCRDHS)
      Real(realk) TMPMAT(MX2CRD*MX2CRD),TMPMT2(MX2CRD,MX2CRD)
      Real(realk) TMPMT3(MX2CRD*MX2CRD),TMPMT4(MX2CRD,MX2CRD)
      Real(realk) VECMOD(MXCOOR)
      PARAMETER (D0 = 0.0E0_realk, D1 = 1.0E0_realk, DP5 = 0.5E0_realk)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
!
!     For saddle point optimizations, we can follow a specific eigenvector.
!     Due to the fact that we are separating one mode for maximization,
!     NCRDHS is temporarily reduced by one.
!
!      IMODE = optinfo%NSPMod
!      NCRDHS = NCRDHS-1
!      IF (optinfo%NSPMod .GT. 0) THEN
!         call ls_FNDMOD(.FALSE.,MXRCRD,EVEC,TMPMAT,VECMOD, &
!     &        TMPMT2,TMPMT3,IMODE,lupri,optinfo)
!
!     If the lowest mode has a gradient element of zero, we have to pick
!     another mode for maximization (or we will end up in a minimum!).
!
!      ELSE
!         IMODE = 1
! 50      CONTINUE
!         IF ((ABS(optinfo%GRDDIA(IMODE)) .LT. 1.0E-10_realk) .AND. &
!     &        (IMODE .LT. NCRDHS)) THEN
!            IMODE = IMODE + 1
!            GOTO 50
!     
!     If we find no such mode, we just set IMODE = 1, because we must
!     be at a stationary point.
!     
!         ELSE IF (ABS(optinfo%GRDDIA(IMODE)) .LT. 1.0E-10_realk) THEN
!            IMODE = 1
!         END IF
!      END IF
!      IF (optinfo%IPrint .GE. IPRMAX) THEN
!         WRITE(LUPRI,*)
!         WRITE(LUPRI,*) 'Mode ',IMODE,' will be partitioned ' // &
!     &        'out and maximized.'
!         WRITE(LUPRI,*)
!      END IF
!
!     The selected mode is placed at the very end.
!
!      call ls_DZERO(TMPMT2,MX2CRD*MX2CRD)
!      call ls_DZERO(TMPMT4,MX2CRD*MX2CRD)
!      DO 400 I = 1, NCRDHS
!         DO 402 J = 1, IMODE-1
!            TMPMT2(I,J) = EVEC(I,J)
! 402     CONTINUE
!         DO 403 J = IMODE, NCRDHS-1
!            TMPMT2(I,J) = EVEC(I,J+1)
! 403     CONTINUE
!         TMPMT4(1,I) = optinfo%EVAL(I)
!         TMPMT4(2,I) = optinfo%GRDDIA(I)
!         TMPMT2(I,NCRDHS) = EVEC(I,IMODE)
! 400  CONTINUE
!      TMPVAL = TMPMT4(1,IMODE)
!      DO 406 I = IMODE, NCRDHS-1
!         TMPMT4(1,I) = TMPMT4(1,I+1)
!         TMPMT4(2,I) = TMPMT4(2,I+1)
! 406  CONTINUE
!      TMPMT4(1,NCRDHS) = TMPVAL
!      TMPMT4(2,NCRDHS) = optinfo%GRDDIA(IMODE)
!
!     We then make the augmented Hessian that will be minimized.
!
!      call ls_DZERO(EVCTMP,NCRDHS*NCRDHS)
!      DO 410 I = 1, NCRDHS-1
!         EVCTMP(I+(I-1)*NCRDHS) = TMPMT4(1,I)
!         EVCTMP(I+(NCRDHS-1)*NCRDHS) = TMPMT4(2,I)
!         EVCTMP(NCRDHS+(I-1)*NCRDHS) = TMPMT4(2,I)
! 410  CONTINUE
!      EVCTMP(NCRDHS+(NCRDHS-1)*NCRDHS) = D0
!      IF (optinfo%IPrint .GE. IPRDBG) THEN
!         call lsheader(lupri,'Augmented Hessian')
!         call ls_output(EVCTMP,1,NCRDHS,1,NCRDHS,NCRDHS, &
!     &        NCRDHS,1,LUPRI)
!      END IF
!      call ls_DZERO(TMPMT3,MX2CRD*MX2CRD)
!      call ls_DSITSP(NCRDHS,EVCTMP,TMPMT3)
!      call ls_DUNIT(EVCTMP,NCRDHS)
!      call ls_JACO(TMPMT3,EVCTMP,NCRDHS,NCRDHS,NCRDHS, &
!     &     TMPMAT(1),TMPMAT(1+MX2CRD))
!      DO 420 J = 1, NCRDHS
!         optinfo%EVAL(J) = TMPMT3(J*(J+1)/2)
!         DO 425 I = 1, NCRDHS
!            EVEC(I,J) = EVCTMP(I+(J-1)*NCRDHS)
! 425     CONTINUE
! 420  CONTINUE
!
!     We add 1.0E5_realk to all eigenvalues that are essentially zero
!     for the sorting.
!
!      DO 427 I = 1, NCRDHS
!         IF (ABS(optinfo%EVAL(I)) .LE. 1.0E-8_realk) optinfo%EVAL(I) = optinfo%EVAL(I) + 1.0E5_realk
! 427  CONTINUE
!      DO 430 I = 1, NCRDHS
!         JMIN = I
!         EMIN = optinfo%EVAL(I)
!         DO 435 J = (I + 1), NCRDHS
!            IF (optinfo%EVAL(J) .LT. EMIN) THEN
!               EMIN = optinfo%EVAL(J)
!               JMIN = J
!            END IF
! 435     CONTINUE
!         IF (JMIN .NE. I) THEN
!            call dswap(1,  optinfo%EVAL  (I),1,optinfo%EVAL  (JMIN),1)
!            call dswap(MX2CRD,EVEC(1,I),1,EVEC(1,JMIN),1)
!!     call dswap(1,optinfo%GRDDIA(I),1,optinfo%GRDDIA(JMIN),1)
!         END IF
! 430  CONTINUE
!      DO 440 I = 1, NCRDHS
!         IF (ABS(ABS(optinfo%EVAL(I))-1.0E5_realk) .LT. 1.0E-3_realk) &
!     &        optinfo%EVAL(I) = optinfo%EVAL(I) - 1.0E5_realk
! 440  CONTINUE
!      IF (optinfo%IPrint .GE. IPRDBG) THEN
!         call lsheader(lupri,'RF-eigenvalues')
!         call ls_output(optinfo%EVAL,1,1,1,NCRDHS,1,MXRCRD,1,LUPRI)
!         call lsheader(lupri,'RF-eigenvectors')
!         call ls_output(EVEC,1,NCRDHS,1,NCRDHS,MX2CRD,MX2CRD,1,LUPRI)
!      END IF
!      call ls_PRFSTP(MX2CRD,NCRDHS,MXCOOR,EVEC,optinfo%STPDIA,EGRAD, &
!     &     TMPMAT,EHESS,1,lupri,optinfo)
!!
!!     In the case of saddle point optimization, we also need
!!     to take care of the second partition and combine the two.
!!
!      CMPLIM = MAX(optinfo%TrustRad*0.67E0_realk, 0.30E0_realk)
!      DO 500 I = 1, NCRDHS-1
!         IF (ABS(optinfo%STPDIA(I)) .GT. CMPLIM) &
!     &        optinfo%STPDIA(I) = SIGN(CMPLIM,optinfo%STPDIA(I))
!         TMPMT4(3,I) = optinfo%STPDIA(I)
! 500  CONTINUE
!      NCRDHS = NCRDHS + 1
!!
!!     We then make the augmented Hessian that will be maximized.
!!
!      call ls_DZERO(EVCTMP,optinfo%ICartCoord*optinfo%ICartCoord)
!      EVCTMP(1) = TMPMT4(1,NCRDHS-1)
!      EVCTMP(2) = TMPMT4(2,NCRDHS-1)
!      EVCTMP(3) = TMPMT4(2,NCRDHS-1)
!      EVCTMP(4) = D0
!      IF (optinfo%IPrint .GE. IPRDBG) THEN
!         call lsheader(lupri,'Augmented Hessian')
!         call ls_output(EVCTMP,1,2,1,2,2,2,1,LUPRI)
!      END IF
!      call ls_DZERO(TMPMT3,MX2CRD*MX2CRD)
!      call ls_DSITSP(2,EVCTMP,TMPMT3)
!      call ls_DUNIT(EVCTMP,2)
!      call ls_JACO(TMPMT3,EVCTMP,2,2,2,TMPMAT(1),TMPMAT(1+MX2CRD))
!      DO 510 J = 1, 2
!         optinfo%EVAL(J) = TMPMT3(J*(J+1)/2)
!         DO 515 I = 1, 2
!            EVEC(I,J) = EVCTMP(I+(J-1)*2)
! 515     CONTINUE
! 510  CONTINUE
!      DO 517 I = 1, 2
!         IF (ABS(optinfo%EVAL(I)) .LE. 1.0E-8_realk) optinfo%EVAL(I) = optinfo%EVAL(I) + 1.0E5_realk
! 517  CONTINUE
!
!     The eigenvalues are sorted
!
!      IF (optinfo%EVAL(1) .GT. optinfo%EVAL(2)) THEN
!         call dswap(1,  optinfo%EVAL  (1),1,optinfo%EVAL  (2),1)
!         call dswap(MX2CRD,EVEC(1,1),1,EVEC(1,2),1)
!      END IF
!      DO 520 I = 1, 2
!         IF (ABS(ABS(optinfo%EVAL(I))-1.0E5_realk) .LT. 1.0E-3_realk) &
!     &        optinfo%EVAL(I) = optinfo%EVAL(I) - 1.0E5_realk
! 520  CONTINUE
!!
!      call ls_PRFSTP(MX2CRD,2,MXCOOR,EVEC,optinfo%STPDIA,EGRAD, &
!     &     TMPMAT,EHESS,2,lupri,optinfo)
!      TMPVAL = optinfo%STPDIA(1)
!      IF (ABS(TMPVAL) .GT. CMPLIM) &
!     &     TMPVAL = SIGN(CMPLIM,TMPVAL)
!      call ls_DZERO(optinfo%STPSYM,MXCOOR)
!      DO 530 I = 1, NCRDHS-2
!         optinfo%STPSYM(I) = TMPMT4(3,I)
! 530  CONTINUE
!      optinfo%STPSYM(NCRDHS-1) = TMPVAL
!      IF (optinfo%IPrint .GE. IPRDBG) THEN
!         call lsheader(lupri,'Diagonal RF-step')
!         call ls_output(optinfo%STPSYM,1,1,1,NCRDHS-1,1,MXRCRD,1,LUPRI)
!      END IF
!
!     The symmetry step is constructed from the original eigenvectors
!     (of the normal Hessian) and the diagonal RF-step.
!
!      call ls_DZERO(optinfo%STPDIA,MXRCRD)
!      call ls_DZERO(EVEC,MX2CRD*MX2CRD)
!      DO 540 I = 1, NCRDHS-1
!         DO 545 J = 1, NCRDHS-1
!            EVEC(I,J) = TMPMT2(I,J)
! 545     CONTINUE
!         optinfo%EVAL(I) = TMPMT4(1,I)
! 540  CONTINUE
!      DO 550 I = 1, NCRDHS-1
!         DO 555 J = 1, NCRDHS-1
!            optinfo%STPDIA(I) = optinfo%STPDIA(I) + EVEC(I,J)*optinfo%STPSYM(J)
! 555     CONTINUE
! 550  CONTINUE
!
!     The final call ls_to RFSTP to do scaling with respect to the
!     trust radius.
!
!      call ls_PRFSTP(MX2CRD,NCRDHS,MXCOOR,EVEC,optinfo%STPDIA,EGRAD, &
!     &     TMPMAT,EHESS,3,lupri,optinfo)
      RETURN
     END

!  /* Deck fndstp */
      SUBROUTINE LS_FNDSTP(MXRCRD,MX2CRD,NCRDHS,EGRAD,EHESS,EVEC,TMPMAT, &
     &     EVCTMP,TMPMT2,TMPMT3,TMPMT4,CSTEP,GRDARR,STPARR,ACTIVE, &
     &     EMOD,VECMOD,lupri,optinfo)
use ls_util 
use precision
use optimization_input
use files
!
!     This routine calculates the step that should be taken to obtain
!     the next geometry.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      TYPE(opt_setting) :: optinfo
      PARAMETER (ZERGRD = 1.0E-7_realk)
      PARAMETER (D0 = 0.0E0_realk , DP5 = 0.5E0_realk)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      Real(realk) EGRAD(MXCOOR), EHESS(MXCOOR,MXCOOR)
      Real(realk) EVEC(MX2CRD,MX2CRD)
      Real(realk) TMPMAT(MX2CRD,MX2CRD),EVCTMP(optinfo%ICartCoord,optinfo%ICartCoord)
      Real(realk) TMPMT2(MX2CRD,MX2CRD),TMPMT3(MX2CRD,MX2CRD)
      Real(realk) TMPMT4(MX2CRD,MX2CRD)
      Real(realk) CSTEP(MXCOOR), GRDARR(MXRCRD,25), STPARR(MXRCRD,25)
      Real(realk) VECMOD(MXCOOR)
      LOGICAL INSIDE, ACTIVE, DOSCAL
      Real(realk),external::DV3DOT_LS
!
      IF (optinfo%LnSearch .AND. (.NOT. optinfo%RatFun) .AND. (optinfo%ItrNmr .GT. 0)) &
     &     call ls_LINSRC(optinfo%ICartCoord,MXCOOR,EGRAD,GRDARR(1,1),CSTEP, &
     &     STPARR(1,1),TMPMAT,TMPMT2,ACTIVE,EMOD,lupri,optinfo)
      IF (ACTIVE) THEN
         DO 5 J = 1, optinfo%ICartCoord
            DO 7 I = 1, optinfo%KEPTIT
               STPARR(J,I) = STPARR(J,I) - CSTEP(J)
 7          CONTINUE
            IF (.NOT. optinfo%RatFun) &
     &           optinfo%GRDDIA(J) = DDOT(optinfo%ICartCoord,EGRAD,1,EVEC(1,J),1)
 5       CONTINUE
         IF (.NOT. optinfo%RatFun) THEN
            IF (optinfo%IPrint .GT. 5) THEN
               call lsheader(lupri,'Diagonal interpolated gradient')
               call ls_output(optinfo%GRDDIA,1,1,1,optinfo%ICartCoord,1,MXRCRD,1,LUPRI)
            END IF
         END IF
      END IF
      call ls_DZERO(optinfo%STPDIA,optinfo%ICartCoord)
      call ls_DZERO(optinfo%STPSYM,optinfo%ICartCoord)
!     Project out all 6 trans-rot coordinates
      optinfo%nProjected = 6
      NVEC = optinfo%ICartCoord-optinfo%nProjected
      IF (optinfo%IPrint .GT. 2) call lsheader(lupri,'Output from FNDSTP')
!
!     First comes the trust region method
!
      optinfo%GradNorm = SQRT(DDOT(NVEC,optinfo%GRDDIA,1,optinfo%GRDDIA,1))
      IF (optinfo%TrustRg .OR. (optinfo%GDIIS .AND. (optinfo%KEPTIT .LT. 3))) THEN
!
!     We take a copy of the eigenvectors, to get correct
!     matrix Real(realk)s.
!
         DO 10 I = 1, optinfo%ICartCoord
            DO 20 J = 1, optinfo%ICartCoord
               EVCTMP(J,I) = EVEC(J,I)
 20         CONTINUE
 10      CONTINUE
!
!     Newton step is calculated.
!
         DO 30 I = 1, NVEC
            optinfo%STPDIA(I) = -optinfo%GRDDIA(I)/optinfo%EVAL(I)
 30      CONTINUE
       
!!!!!!!!  optinfo%StepNorm = SQRT(DDOT(NVEC,optinfo%STPDIA,1,optinfo%STPDIA,1))
!
!     If Newton step is larger that trust radius, we take a step
!     to the boundary. If the Hessian index is larger than zero,
!     the level-shifted step will also be employed, provided the
!     Newton step is larger than 0.5E-3_realk. For saddle points we
!     employ the level-shift when the index is different from 1.
!
         IF (((optinfo%StepNorm .GT. optinfo%TrustRad) .AND. (.NOT. optinfo%NoTrust)) .OR. &
     &        ((.NOT. optinfo%Saddle) .AND. (optinfo%IndHes .GT. 0) .AND. &
     &        (optinfo%StepNorm .GE. 0.5E-3_realk)) .OR. (optinfo%Saddle .AND. &
     &        (optinfo%IndHes .NE. 1))) THEN
            IF (optinfo%IPrint .GT. 5) THEN
               WRITE(LUPRI,'(/A,F15.10)')' Norm of Newton step:', optinfo%StepNorm
               WRITE(LUPRI,'(A,F15.10/)')' Trust radius       :', optinfo%TrustRad
            END IF
            INSIDE = .FALSE.
            IF (optinfo%StepNorm .LT. optinfo%TrustRad) INSIDE = .TRUE.
            call ls_LSHFT0(optinfo%ICartCoord,NVEC,&
                 &  MIN(optinfo%TrustRad,optinfo%StepNorm),RNU,ZERGRD,INSIDE,lupri,optinfo)
            optinfo%StepNorm = SQRT(DDOT(NVEC,optinfo%STPDIA,1,optinfo%STPDIA,1))
            IF (optinfo%IPrint .GT. 5) THEN
               WRITE(LUPRI,'(/A,F15.10)')' Norm, boundary step:', optinfo%StepNorm
            END IF
         END IF
!
!     Energy is predicted, will be used later to update trust radius.
!
         optinfo%predictedChange = DDOT(NVEC,optinfo%GRDDIA,1,optinfo%STPDIA,1) &
     & + 0.5E0_realk*DV3DOT_LS(NVEC,optinfo%STPDIA,optinfo%EVAL, &
     &optinfo%STPDIA)
         WRITE(LUPRI,'(/A,F25.15)') ' Predicted energy change',optinfo%predictedChange
!
!     If the predicted energy is positive, it means the Newton step is
!     towards a maximum/saddle point. We then simply reverse the
!     step direction (a bit dirty, but seems to work).
!
         IF ((.NOT. optinfo%Saddle) .AND. (optinfo%predictedChange .GT. 0.0E0_realk)) THEN
            WRITE(LUPRI,*) 'Reversing step!'
            DO 40 I = 1, NVEC
               optinfo%STPDIA(I) = -optinfo%STPDIA(I)
 40         CONTINUE
            optinfo%predictedChange = DDOT(NVEC,optinfo%GRDDIA,1,optinfo%STPDIA,1) &
     & + 0.5E0_realk*DV3DOT_LS(NVEC,optinfo%STPDIA,optinfo%EVAL, &
     &  optinfo%STPDIA)
            IF (optinfo%IPrint .GT. 2) THEN
               WRITE(LUPRI,'(A,F25.15)') &
     &              ' New pred. energy change',optinfo%predictedChange
            END IF
         END IF
         call lsheader(lupri,'Step in diagonal representation')
         call ls_output(optinfo%STPDIA,1,1,1,NVEC,1,MXRCRD,1,LUPRI)
         IF (optinfo%IPrint .GT. 5) THEN
            call lsheader(lupri,'Eigenvector basis')
            call ls_output(EVCTMP,1,optinfo%ICartCoord,1,NVEC,optinfo%ICartCoord,optinfo%ICartCoord,1,LUPRI)
         END IF
         DO 150 I = 1, NVEC
            call daxpy(optinfo%ICartCoord,optinfo%STPDIA(I),EVEC(1,I),1,optinfo%STPSYM,1)
 150     CONTINUE
!
!     The rational function method
!
      ELSE IF (optinfo%RatFun) THEN
         IF (optinfo%Saddle) THEN
            call ls_PRFSTC(MXRCRD,MX2CRD,NCRDHS,EGRAD,EHESS,EVEC,EVCTMP, &
     &           TMPMAT,TMPMT2,TMPMT3,TMPMT4,VECMOD,lupri,optinfo)
         ELSE
            call ls_RFSTP(MX2CRD,NCRDHS,MXCOOR,optinfo%ICartCoord,EVEC,optinfo%STPDIA,EGRAD, &
     &           TMPMAT,EHESS,lupri,optinfo)
         END IF
         IF (optinfo%IPrint .GE. IPRMIN) THEN
            WRITE (LUPRI,'(/A,F25.15)') &
     &           ' Predicted energy change', optinfo%predictedChange
         END IF
         call ls_DZERO(optinfo%STPSYM,MXCOOR)
         DO 200 I = 1, optinfo%ICartCoord
            IF (ABS(optinfo%STPDIA(I)) .GE. 1.0E-6_realk) THEN
               optinfo%STPSYM(I) = optinfo%STPDIA(I)
            ELSE
               optinfo%STPSYM(I) = D0
            END IF
 200     CONTINUE
!
!     The Geometrical DIIS method
!
      ELSE IF (optinfo%GDIIS) THEN
         call ls_DZERO(TMPMAT,MX2CRD*MX2CRD)
         call ls_DZERO(EVCTMP,optinfo%ICartCoord*optinfo%ICartCoord)
!
!     First we have to construct the inverse Hessian.
!
         DO 210 I = 1, NVEC
            DO 212 J = 1, optinfo%ICartCoord
               TMPMAT(I,J) = EVEC(J,I)/optinfo%EVAL(I)
 212        CONTINUE
 210     CONTINUE
         DO 215 I = 1, optinfo%ICartCoord
            DO 217 J = 1, optinfo%ICartCoord
               DO 219 K = 1, NVEC
                  EVCTMP(I,J) = EVCTMP(I,J) + EVEC(I,K)*TMPMAT(K,J)
 219           CONTINUE
 217        CONTINUE
 215     CONTINUE
!
!     Then the DIIS-step is determined
!
         call ls_GDISTP(MXCOOR,optinfo%ICartCoord,MXRCRD,MX2CRD,optinfo%STPDIA,EGRAD,EHESS, &
     &        EVCTMP,TMPMAT,TMPMT2,TMPMT3,TMPMT4,GRDARR,STPARR,lupri,optinfo)
         IF (optinfo%IPrint .GE. IPRMIN) THEN
            WRITE (LUPRI,'(/A,F25.15)') &
     &           ' Predicted energy change', optinfo%predictedChange
         END IF
         DO 250 I = 1, optinfo%ICartCoord
            IF (ABS(optinfo%STPDIA(I)) .GE. 1.0E-6_realk) THEN
               optinfo%STPSYM(I) = optinfo%STPDIA(I)
            ELSE
               optinfo%STPSYM(I) = D0
            END IF
 250     CONTINUE
      END IF
      IF (ACTIVE) THEN
         DO 300 I = 1, optinfo%ICartCoord
            optinfo%STPSYM(I) = optinfo%STPSYM(I) + CSTEP(I)
 300     CONTINUE
         optinfo%predictedChange = optinfo%predictedChange + (EMOD-optinfo%energy)
         WRITE (LUPRI,'(/A,F25.15)') &
     &        ' Modified energy prediction due to line search', optinfo%predictedChange
      Write(*,*)'CSTEP',CSTEP
      END IF
      optinfo%StepNorm = SQRT(DDOT(optinfo%ICartCoord,optinfo%STPSYM,1,optinfo%STPSYM,1))
      IF (optinfo%IPrint .GT. 2) THEN
         call lsheader(lupri,'Cartesian step vector')
         call ls_output(optinfo%STPSYM,1,1,1,optinfo%ICartCoord,1,optinfo%ICartCoord,1,LUPRI)
         WRITE(LUPRI,'(/A,F15.10/)') ' Norm of step:', optinfo%StepNorm
      END IF
!
      Do i = 1,MXCOOR
         CSTEP(i) = optinfo%STPSYM(i)
      Enddo
      RETURN
      END
!  /* Deck ls_minend */
      LOGICAL FUNCTION ls_minend(MXRCRD,BMTRAN,TMPVC1,TMPVC2,lupri,optinfo)
use ls_util
use files
use optimization_input
use precision
use dec_typedef_module,only: DECinfo
!
!     Determines if the end of the optimization has been reached.
!
Implicit Real(realk) (A-H,O-Z)
      TYPE(opt_setting) :: optinfo
      Real(realk) BMTRAN(MXRCRD,MXRCRD)
      Real(realk) TMPVC1(MXRCRD), TMPVC2(MXRCRD)
      LOGICAL CNVERG, CNVGRD, CNVSTP, INDXOK, CNVRGD(1:4)
      CHARACTER*4 LOGTXT
      CNVERG = .FALSE.
      CNVGRD = .FALSE.
      CNVSTP = .FALSE.
      ICOOR = optinfo%ICartCoord
      IF (optinfo%RedInt .OR. optinfo%DelInt) ICOOR = optinfo%NIntCoord
      IF (optinfo%ConOpt) GNRM = SQRT(DDOT(optinfo%NIntCoord,optinfo%GRDINT,1,optinfo%GRDINT,1))
!
!     First, the convergence criteria of Baker.
!     To make the convergence criteria more fair, we transform
!     the delocalized internals to the primitive space.
!
      IF (optinfo%Baker) THEN
         IF (optinfo%DelInt) THEN
            call ls_DZERO(TMPVC1,MXRCRD)
            call ls_DZERO(TMPVC2,MXRCRD)
            DO 10 I = 1, optinfo%NIntCoord
               DO 20 J = 1, optinfo%NIntCoord
                  TMPVC1(I) = TMPVC1(I) + BMTRAN(I,J)*optinfo%GRDINT(J)
                  TMPVC2(I) = TMPVC2(I) + BMTRAN(I,J)*optinfo%STPINT(J)
 20            CONTINUE
 10         CONTINUE
            call ls_MAXELM(TMPVC1,optinfo%NIntCoord,0,GRDMAX,lupri,optinfo)
            call ls_MAXELM(TMPVC2,optinfo%NIntCoord,0,STPMAX,lupri,optinfo)
         ELSE IF (optinfo%RedInt) THEN
            call ls_MAXELM(optinfo%GRDINT,optinfo%NIntCoord,0,GRDMAX, &
            & lupri,optinfo)
            call ls_MAXELM(optinfo%STPINT,optinfo%NIntCoord,0,STPMAX, &
            & lupri,optinfo)
         ELSE
            call ls_MAXELM(optinfo%GradMol,ICOOR,2,GRDMAX,lupri,optinfo)
            call ls_MAXELM(optinfo%STPSYM,ICOOR,1,STPMAX,lupri,optinfo)
         END IF
         EDIFF  = ABS(optinfo%energyOld-optinfo%energy)
         INDXOK = (optinfo%IndHes .EQ. 0)
         IF (optinfo%Saddle) INDXOK = (optinfo%IndHes .EQ. 1)
!
         CNVGRD = (GRDMAX .LE. 3.0E-4_realk)
         CNVERG = (EDIFF  .LE. 1.0E-6_realk)
         IF ((optinfo%ItrNmr .LT. 1) .OR. (optinfo%ItBreak .EQ. (optinfo%ItrNmr-1))) &
     &          CNVERG = .FALSE.
         CNVSTP = (STPMAX .LE. 3.0E-4_realk)
         ls_minend = (INDXOK .AND. CNVGRD) .AND. (CNVERG .OR. CNVSTP)
         IF ( (.NOT. optinfo%Newton) .OR. optinfo%Saddle) &
     &        ls_minend =  (CNVGRD .AND. (CNVERG .OR. CNVSTP))
!
!     Next, the default (old) convergence check
!
      ELSE IF (.NOT. optinfo%NatNorm) THEN
         optinfo%IConV = 0
         EDIFF = ABS(optinfo%energyOld-optinfo%energy)
         IF (EDIFF .LE. optinfo%ThrErg) THEN
            CNVERG = .TRUE.
            optinfo%IConV = optinfo%IConV + 1
         END IF
         IF (optinfo%ConOpt) THEN
            IF (GNRM .LE. optinfo%GradThr) THEN
               CNVGRD = .TRUE.
               optinfo%IConV = optinfo%IConV + 1
            END IF
         ELSE
            IF (optinfo%GradNorm .LE. optinfo%GradThr) THEN
               CNVGRD = .TRUE.
               optinfo%IConV = optinfo%IConV + 1
            END IF
         END IF
         IF (optinfo%StepNorm .LE. optinfo%ThrStep) THEN
            CNVSTP = .TRUE.
            optinfo%IConV = optinfo%IConV + 1
         END IF
         INDXOK = (optinfo%IndHes .EQ. 0)
         IF (optinfo%Saddle) INDXOK = (optinfo%IndHes .EQ. 1)
         ls_minend = INDXOK .AND. (optinfo%IConV .GE. optinfo%Condi)
         ls_minend = (optinfo%IConV .GE. optinfo%Condi)
!
!     Finally, the new convergence scheme (optinfo%NatNorm)
!
      ELSE
         IF (optinfo%RedInt) THEN
            call ls_MAXELM(optinfo%GRDINT,optinfo%NIntCoord,0,GRDMAX, &
     &      lupri,optinfo)
            call ls_MAXELM(optinfo%STPINT,optinfo%NIntCoord,0,STPMAX, &
     &      lupri,optinfo)
            SNRM = SQRT(DDOT(optinfo%NIntCoord,optinfo%STPINT,1,optinfo%STPINT,1)/ &
     &           (1E0_realk*optinfo%NIntCoord))
            GNRM = SQRT(DDOT(optinfo%NIntCoord,optinfo%GRDINT,1,optinfo%GRDINT,1)/ &
     &           (1E0_realk*optinfo%NIntCoord))
         ELSE
            call ls_MAXELM(optinfo%GradMol,ICOOR,2,GRDMAX,lupri,optinfo)
            call ls_MAXELM(optinfo%STPSYM,ICOOR,1,STPMAX,lupri,optinfo)
            SNRM = optinfo%StepNorm/SQRT(1E0_realk*optinfo%ICartCoord)
            GNRM = optinfo%GradNorm/SQRT(1E0_realk*optinfo%ICartCoord)
         END IF
         CNVRGD(1) = .FALSE.
         CNVRGD(2) = .FALSE.
         CNVRGD(3) = .FALSE.
         CNVRGD(4) = .FALSE.
         IF (GNRM .LT. optinfo%GradThr) CNVRGD(1)   = .TRUE.
         IF (GRDMAX .LT. optinfo%ThGradMax) CNVRGD(2) = .TRUE.
         IF (SNRM .LT. optinfo%ThrStep) CNVRGD(3)   = .TRUE.
         IF (STPMAX .LT. optinfo%ThStepMax) CNVRGD(4) = .TRUE.
!
         IF (optinfo%ConOpt) THEN
!
!     For constrained optimizations we examine the change in the RMS gradient
!     and the change in the maximum element. The threshold for these differences
!     are 1/10 of the current regular gradient thresholds.
!
            IF (optinfo%ItrNmr .EQ. 0) THEN
               optinfo%GradThr = MAX(0.1E0_realk*optinfo%GradThr,1E-7_realk)
               optinfo%ThGradMax = 5E0_realk*optinfo%GradThr
               DIFRMS = 0E0_realk
               DIFMAX = 0E0_realk
               optinfo%PRVRMS = GNRM
               optinfo%PRVMAX = GRDMAX
            ELSE
               DIFRMS = ABS(optinfo%PRVRMS-GNRM)
               DIFMAX = ABS(optinfo%PRVMAX-GRDMAX)
               optinfo%PRVRMS = GNRM
               optinfo%PRVMAX = GRDMAX
            END IF
            IF (DIFRMS .LT. optinfo%GradThr) CNVRGD(1) = .TRUE.
            IF (DIFMAX .LT. optinfo%ThGradMax) CNVRGD(2) = .TRUE.
            WRITE(LUPRI,'(///A)') ' -----------------------' // &
     &     '-----------------------------------------------------------'
            WRITE(LUPRI,'(14X,A)') 'RMS grad    Max grad    ' // &
     &           'RMS diff*   Max diff*   RMS step    Max step'
            WRITE(LUPRI,'(A)')' -----------------------' // &
     &     '-----------------------------------------------------------'
            WRITE(LUPRI,'(A,6F12.8)') ' Curr.val. ', &
     &           GNRM, GRDMAX, DIFRMS, DIFMAX, SNRM, STPMAX
            WRITE(LUPRI,'(A,24X,4F12.8)') ' Threshold ', &
     &           optinfo%GradThr, optinfo%ThGradMax, optinfo%ThrStep, optinfo%ThStepMax
            WRITE(LUPRI,'(A,24X,4(L1,10X))') ' Converged?      ', &
     &           (CNVRGD(I), I = 1,4)
            WRITE(LUPRI,'(A)')' -----------------------' // &
     &     '-----------------------------------------------------------'
            WRITE(LUPRI,'(A)') '*) For constrained optimizations' // &
     &           ' the change in RMS gradient and the change'
            WRITE(LUPRI,'(A///)') '   in the maximum gradient ' // &
     &           'element are used as convergence criteria.'
         ELSE
            WRITE(LUPRI,'(///A)') &
     &     ' ----------------------------------------------------------' &
     &       //  '--+---------------------------'
            WRITE(LUPRI,'(14X,A)') &
     &           'RMS grad    Max grad    RMS step    Max step' // &
     &           '   |  Tr.rad.        Energy'
            WRITE(LUPRI,'(A)') &
     &     ' ----------------------------------------------------------' &
     &       //  '--+---------------------------'
            WRITE(LUPRI,'(A,4F12.8,A,F10.6,F17.8)') ' Curr.val. ', &
     &           GNRM, GRDMAX, SNRM, STPMAX, '  |', optinfo%TrustRad, optinfo%energy
            WRITE(LUPRI,'(A,4F12.8,A)') ' Threshold ', &
     &           optinfo%GradThr, optinfo%ThGradMax, optinfo%ThrStep, optinfo%ThStepMax, '  |'
            WRITE(LUPRI,'(A,3(L1,10X),L1,6X,A)') ' Converged?      ', &
     &           (CNVRGD(I), I = 1,4), '|'
            WRITE(LUPRI,'(A///)') &
     &     ' ----------------------------------------------------------' &
     &       //  '--+---------------------------'
         END IF
!
         INDXOK = (optinfo%IndHes .EQ. 0)
         IF (optinfo%Saddle) INDXOK = (optinfo%IndHes .EQ. 1)
         ls_minend = INDXOK .AND. CNVRGD(1) .AND. CNVRGD(2) .AND. &
     &        CNVRGD(3) .AND. CNVRGD(4)
         ls_minend = CNVRGD(1) .AND. CNVRGD(2) .AND. &
     &        CNVRGD(3) .AND. CNVRGD(4)
!!$         if(DECinfo%dodec) then ! only check gradient itself for DEC (not step size)
!!$            ls_minend = CNVRGD(1) .AND. CNVRGD(2)
!!$         end if
      END IF
!
!     Output from the testing is written.
!
      IF (.NOT. optinfo%NatNorm .AND. (optinfo%IPrint .GT. 2)) THEN
         call lsheader(lupri,'Output from ls_minend')
         LOGTXT = 'no  '
         IF (CNVERG) LOGTXT = 'yes '
         IF ((optinfo%ItrNmr .LT. 1) .OR. (optinfo%ItBreak .EQ. (optinfo%ItrNmr-1))) &
     &        LOGTXT = 'N/A '
         WRITE(LUPRI,'(/A,1P,A5)') ' Energy converged      ',LOGTXT
         LOGTXT = 'no  '
         IF (CNVGRD) LOGTXT = 'yes '
         WRITE(LUPRI,'(A,1P,A5)') ' Gradient converged    ',LOGTXT
         LOGTXT = 'no  '
         IF (CNVSTP) LOGTXT = 'yes '
         WRITE(LUPRI,'(A,1P,A5)') ' Step converged        ',LOGTXT
         IF (.NOT. optinfo%Baker) THEN
            WRITE(LUPRI,'(A,1P,I3)') ' Conditions fullfilled ',optinfo%IConV
            WRITE(LUPRI,'(A,1P,I3)') ' Required conditions   ',optinfo%Condi
         END IF
         WRITE(LUPRI,'(A,1P,I3)') ' Totally sym. index    ',optinfo%INDHES
         WRITE(LUPRI,'(A,1P,I3)') ' Hessian index         ',optinfo%IndHes
         LOGTXT = 'no  '
         IF (ls_minend) LOGTXT = 'yes '
         WRITE(LUPRI,'(A,1P,A5)') ' End of optimization   ',LOGTXT
      END IF
      RETURN
      END

!  /* Deck lshft0 */
      SUBROUTINE LS_LSHFT0(NCORD,NONTRO,TRUSTR,RNU, &
      & ZERGRD,INSIDE,lupri,optinfo)
use ls_util 
use files
use precision
use optimization_input
!
!     (Almost identical to WLKFL0 in abawalk.F)
!     This subroutine solves the constrained restricted step
!     equations (the level-shifted Newton equations) in the
!     diagonal representation.  We assume that the Newton step
!     is longer than the trust radius.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      TYPE(opt_setting) :: optinfo
      PARAMETER ( D0 = 0.0E0_realk , DP5 = 0.5E0_realk )
!
      LOGICAL KEEPSY, INSIDE, SPECAS
!
      EXTERNAL wstpln_ls
!
      IF (optinfo%IPrint .GT. 5) call lsheader(lupri,'OUTPUT FROM LSHFT0')
!
      HESMIN = optinfo%EVAL(1)
      GRDMIN = ABS(optinfo%GRDDIA(1))

!
!        Test whether the lowest Hessian eigenvalue is negative and the
!        corresponding gradient zero. This case is treated separately
!        as described by Fletcher in "Unconstrained Optimization" p.85.
!
      SPECAS = (HESMIN .LT. D0) .AND. (GRDMIN .LT. ZERGRD)
      optinfo%GradNorm = SQRT(DDOT(NONTRO,optinfo%GRDDIA,1,optinfo%GRDDIA,1))
      IF (optinfo%IPrint .GT. 5) THEN
         WRITE (LUPRI,'(A,1P,D12.5)') ' HESMIN: ', HESMIN
         WRITE (LUPRI,'(A,1P,D12.5)') ' GRDMIN: ', GRDMIN
         WRITE (LUPRI,'(A,1P,D12.5)') ' optinfo%GradNorm: ', optinfo%GradNorm
         WRITE (LUPRI,'(A,1P,D12.5)') ' ZERGRD: ', ZERGRD
      END IF
!
!     ************************
!     ***** General case *****
!     ************************
!
      IF (.NOT. SPECAS) THEN
         IF (optinfo%IPrint .GT. 3) WRITE (LUPRI,'(/A/)') ' General case.'
!
!        Determine level shift
!
         XMIN = optinfo%GradNorm/TRUSTR - MIN(HESMIN,D0)
         XMAX = MAX(D0, - HESMIN + DP5*GRDMIN/TRUSTR)
         IF (optinfo%IPrint .GT. 5) THEN
            WRITE (LUPRI,'(/A,2(1P,D12.5,2X))') &
     &         ' XMIN and XMAX before WLKBIS (wstpln_ls) minimum walk: ', &
     &          XMIN,XMAX
         END IF
        call ls_WLKBIS(XMAX,XMIN,RNU,optinfo%GRDDIA,optinfo%EVAL,TRUSTR,NONTRO, &
    &               wstpln_ls,IFAIL)
         IF (optinfo%IPrint .GT. 5) THEN
            WRITE (LUPRI,'(/A,2(1P,D12.5,2X))') &
     &         ' XMIN and XMAX after WLKBIS (wstpln_ls) minimum walk: ', &
     &          XMIN,XMAX
         END IF
 33      CONTINUE
         IF (optinfo%IPrint .GT. 3) WRITE (LUPRI,'(/A,1P,D12.5)') &
     &      ' Level shift parameter: ', RNU
!     If we're using level-shift for a step less than the trust radius,
!     and we're having trouble with the interval, we simply set
!     the level-shift parameter to zero.
!
         IF ((IFAIL.EQ. 0) .AND. INSIDE) THEN
            WRITE (LUPRI,'(/A)') &
     &           ' *** ERROR, Wrong interval in WLKBIS (wstpln_ls)'
            WRITE (LUPRI,'(A)') &
     &           '     Setting level-shift equal to zero.'
            RNU = 0E0_realk
            IFAIL = -1
            GOTO 33
         ELSE IF (IFAIL.EQ. 0) THEN
            WRITE (LUPRI,5250) XMAX, XMIN
            call lsquit(' *** ERROR, Wrong interval in WLKBIS (wstpln_ls)',lupri)
         ELSE IF (IFAIL.EQ. 1) THEN
            WRITE (LUPRI,5350)
         END IF
!
!        Determine step vector
!
         DO 100 I = 1, NONTRO
            IF (ABS(optinfo%EVAL(I) + RNU) .LE. 1.0E-8_realk) THEN
               optinfo%STPDIA(I) = D0
            ELSE
               optinfo%STPDIA(I) = - optinfo%GRDDIA(I)/(optinfo%EVAL(I) + RNU)
            END IF
 100     CONTINUE
!
!     *************************************************
!     ***** Special case: HESMIN < 0 & GRDMIN = 0 *****
!     *************************************************
!
      ELSE
         IF (optinfo%IPrint .GT. 3) WRITE (LUPRI,'(/A/)') ' Special case.'
         IMODE = 2
 150     CONTINUE
         IF (IMODE .LT. NONTRO) THEN
            IF ((optinfo%EVAL(IMODE) .LT. D0) .AND. &
     &           (ABS(optinfo%GRDDIA(IMODE)) .LT. ZERGRD)) THEN
               IMODE = IMODE + 1
               GOTO 150
            END IF
         END IF
         GRDMIN=optinfo%GRDDIA(IMODE)
         HESMIN=optinfo%EVAL(IMODE)
!
!        Set RNU = - HESMIN and determine step length
!
         optinfo%StepNorm = wstpln_ls(optinfo%GRDDIA(IMODE:),optinfo%EVAL(IMODE:),-HESMIN, &
     &        NONTRO-IMODE+1,D0)
         IF (optinfo%IPrint .GT. 3) THEN
            WRITE (LUPRI,'(/A,F12.6)') &
     &      ' Step length with level shift equal to lowest eigenvalue:', &
     &      optinfo%StepNorm
         END IF
         IF (optinfo%StepNorm .GT. TRUSTR) THEN
            IF (optinfo%IPrint .GT. 3) WRITE (LUPRI,'(/A)') &
     &         ' Component along lowest eigenvector ignored.'
!
!           Determine step vector in the usual way ignoring the
!           component along the lowest eigenvector. We now know that
!           the level shift must be greater than - HESMIN.
!
            XMIN = optinfo%GradNorm/TRUSTR - MIN(HESMIN,D0)
            XMAX = - HESMIN
            IF (optinfo%IPrint .GT. 5) THEN
               WRITE (LUPRI,'(/A,2(1P,D12.5,2X))') &
     &           ' XMIN and XMAX before WLKBIS (wstpln_ls) minimum walk: ', &
     &           XMIN,XMAX
            END IF
            call ls_WLKBIS(XMAX,XMIN,RNU,optinfo%GRDDIA(IMODE:),optinfo%EVAL(IMODE:),TRUSTR, &
    &           NONTRO-IMODE+1,wstpln_ls,IFAIL)
            IF (optinfo%IPrint .GT. 5) THEN
               WRITE (LUPRI,'(/A,2(1P,D12.5,2X))') &
     &         ' XMIN and XMAX after WLKBIS (wstpln_ls) minimum walk: ', &
     &           XMIN,XMAX
            END IF
            IF (optinfo%IPrint .GT. 3) WRITE (LUPRI,'(/A,1P,D12.5)') &
     &         ' Level shift parameter: ', RNU
            IF (IFAIL.EQ. 0) THEN
               WRITE (LUPRI,5250) XMAX, XMIN
               call lsquit &
     &            (' *** ERROR, Wrong interval in WLKBIS (wstpln_ls)',lupri)
            ELSE IF (IFAIL.EQ. 1) THEN
               WRITE (LUPRI,5350)
            END IF
!
!           Determine step vector
!
            DO 199 I = 1, IMODE-1
               optinfo%STPDIA(I) = D0
 199        CONTINUE
            DO 200 I = IMODE, NONTRO
               optinfo%STPDIA(I) = - optinfo%GRDDIA(I)/(optinfo%EVAL(I) + RNU)
 200        CONTINUE
         ELSE
!
!           Determine step vector with level shift - HESMIN and add
!           component along the lowest eigenvector(s) to insure that total
!           step length is equal to the trust radius.
!
            DO 300 I = IMODE, NONTRO
               optinfo%STPDIA(I) = - optinfo%GRDDIA(I)/(optinfo%EVAL(I) - HESMIN)
 300        CONTINUE
            STP2 = DDOT(NONTRO-IMODE+1,optinfo%STPDIA(IMODE:),1,optinfo%STPDIA(IMODE:),1)
            IF (optinfo%IPrint .GT. 3) WRITE (LUPRI,'(/A,F12.6)') &
     &         ' Norm of step orthogonal to lowest eigenvector(s):', &
     &         SQRT(STP2)
            SXTRA =  SQRT(TRUSTR*TRUSTR - STP2)/SQRT(1.0E0_realk*(IMODE-1))
            IF (optinfo%IPrint .GT. 3) WRITE (LUPRI,'(A,F12.6)') &
     &         ' Norm of step parallel to lowest eigenvector(s):  ', &
     &         SQRT((IMODE-1)*SXTRA*SXTRA)
            DO 400 I = 1, IMODE-1
               optinfo%STPDIA(I) = SXTRA
 400        CONTINUE
         END IF
      END IF
      IF (optinfo%IPrint .GT. 2) THEN
         call lsheader(lupri,'optinfo%STPDIA after min search')
         WRITE (LUPRI,'(5X,3F15.8)') (optinfo%STPDIA(I),I=1,NCORD)
      END IF
      RETURN
!
!     FORMATS
!
 5250 FORMAT(/' *** Wrong interval specified in WLKBIS (wstpln_ls) ***', &
     &       /' XMAX= ',F10.6,'   XMIN= ',F10.6)
 5350 FORMAT(/' *** WARNING WLKBIS (wstpln_ls) ***', &
     &       /' Desired accuracy not obtained in the specified maximum', &
     &       /' number of iterations.')
      END
 
!  /* Deck reahes */
!      SUBROUTINE LS_REAHES(MXRCRD,MX2CRD,HESINT,ATMARR,TMPMT1,TMPMT2, &
!     &     TMPMT3,TMPMT4,WILBMT,BMTRAN,BMTINV,WORK,LWORK,IERR,lupri,optinfo)
!use ls_util 
!use optimization_input
!use files
!!
!!     Read molecular Hessian from the file DALTON.HES, which is then
!!     used as initial Hessian in 1st order methods/restarts. IERR is
!!     returned with the value 0 if everything is OK. -1 indicates that
!!     the file cannot be opened, -2 that the Hessian in the file has
!!     wrong Real(realk)s.
!!
!Implicit Real(realk) (A-H,O-Z)
!#include "dummy.h"
!#include "maxorb.h"
!#include "mxcent.h"
!#include "maxaqn.h"
!#include "nuclei.h"
!#include "molinp.h"
!#include "taymol.h"
!#include "symmet.h"
!      Integer :: lupri
!      TYPE(opt_setting) :: optinfo
!      LOGICAL HESEXS
!      Real(realk) HESINT(MXRCRD,MXRCRD), ATMARR(MXCENT,8)
!      Real(realk) TMPMT1(MX2CRD,MX2CRD), TMPMT2(MX2CRD,MX2CRD)
!      Real(realk) TMPMT3(MX2CRD,MX2CRD), TMPMT4(MXCOOR,MXCOOR)
!      Real(realk) WILBMT(MXRCRD,MXCOOR), BMTRAN(MXRCRD,MXRCRD)
!      Real(realk) BMTINV(MXRCRD,MXCOOR), WORK(LWORK)
!!
!      LUHES = -1
!      call ls_DZERO(optinfo%HessMol,MXCOOR*MXCOOR)
!      INQUIRE(FILE='DALTON.HES',EXIST=HESEXS)
!      IF (.NOT. HESEXS) THEN
!         IERR = -1
!         RETURN
!      ELSE
!         IERR = 0
!      END IF
!      call lsopen(LUHES,'DALTON.HES','OLD','FORMATTED')
!      READ(LUHES,*) IDIM
!      READ(LUHES,*)
!      ICRD = 3*NUCDEP
!      IF (IDIM .NE. ICRD) THEN
!         IERR = -2
!         RETURN
!      END IF
!      call ls_DZERO(TMPMT4,MXCOOR*MXCOOR)
!      DO 100 J = 1, ICRD
!         DO 110 I = 1, ICRD
!            READ(LUHES,*) TMPMT4(I,J)
! 110     CONTINUE
!         READ(LUHES,*)
! 100  CONTINUE
!      DO 150 J = 1, ICRD
!         DO 160 I = 1, ICRD
!            optinfo%HessMol(I,J) = TMPMT4(I,J)
! 160     CONTINUE
! 150  CONTINUE
!      IF (optinfo%RedInt .OR. optinfo%DelInt) THEN
!         call ls_DZERO(HESINT,MXRCRD*MXRCRD)
!         call ls_HX2HQ(MXRCRD,MX2CRD,ATMARR,TMPMT1,TMPMT2,TMPMT3, &
!     &        TMPMT4,optinfo%GRDINT,HESINT,WILBMT,BMTINV,BMTRAN, &
!     &        optinfo)
!      END IF
!      call lsclose(LUHES,'KEEP')
!      RETURN
!      END

!  /* Deck maxelm */
      SUBROUTINE LS_MAXELM(VEC,IDIM,ISCL,ELMMX,lupri,optinfo)
use ls_util 
use precision
use optimization_input
use files
!
!     Finds the largest elemement (absolute value) of the vector VEC
!     of Real(realk) IDIM. The value is returned through the
!     variable ELMMX.
!
Implicit Real(realk) (A-H,O-Z)
      Real(realk) VEC(IDIM)
      Integer :: lupri
      Type(opt_setting) :: optinfo
      LOGICAL SCALE
      ELMMX = ABS(VEC(1))
      IF (ISCL .GE. 1) ELMMX = ELMMX
      IF (ISCL .GE. 2) ELMMX = ELMMX
      DO 10 I = 2, IDIM
         ELM = ABS(VEC(I))
         IF (ISCL .GE. 1) ELM = ELM
         IF (ISCL .GE. 2) ELM = ELM
         IF (ELM .GT. ELMMX) ELMMX = ELM
 10   CONTINUE
      RETURN
      END

!  /* Deck wstpln_ls */
      FUNCTION wstpln_ls(GDDIA,HESDIA,RNU,NCORD,RTRUST)
!
!     (almost identical to WLKSTL in abawalk.F)
!     Purpose:
!
!        Calculate step length at level shift RNU and
!        subtract RTRUST
!
!        wstpln_ls = //STEP// - RTRUST
!
!        where
!
!        STEP = - GDDIA/(HESDIA+RNU)
!
use precision
Implicit Real(realk) (A-H,O-Z)
      Real(realk) GDDIA(*),HESDIA(*)
      PARAMETER(D0=0.0E0_realk, ZERO=1.0E-8_realk )
      STEP = 0E0_realk
      DO 100 K=1,NCORD
         IF ((ABS(GDDIA(K)) .GE. ZERO) .AND. &
     &        (ABS(HESDIA(K)) .GE. ZERO)) THEN
            STEPK = GDDIA(K) / (HESDIA(K)+RNU)
            STEP = STEP + STEPK*STEPK
         END IF
 100  CONTINUE
      wstpln_ls = SQRT(STEP) - RTRUST
      RETURN
      END

!  /* Deck dunit */
      SUBROUTINE LS_DUNIT(A,N)
!
!  SUBROUTINE DUNIT SETS THE REAL SQUARE MATRIX A EQUAL
!  TO A UNIT MATRIX.
!  /VER 2/ 14-Sep-1983 hjaaj
!
use precision
Implicit Real(realk) (A-H,O-Z)
      Real(realk) A(*)
      PARAMETER (D1=1.0E00_realk, D0=0.0E00_realk)
!
      NN = N*N
      DO 100 I = 1,NN
         A(I) = D0
  100 CONTINUE
      N1 = N + 1
      DO 200 I = 1,NN,N1
         A(I) = D1
  200 CONTINUE
      RETURN
      END
!================!
! LS_copyGH      !
!================!
Subroutine ls_copyGH(EGRAD,KEHESS,KALHES,optinfo,NCoord)
!
! Copies molecular gradient(optinfo%GradMol) and Hessian(optinfo%HessMol) to vectors
! EGRAD,KALHES,KEHESS
Use precision
use optimization_input
Implicit none
Type(opt_setting) :: optinfo
Integer :: NCoord,i,j,ji
Real(realk) EGRAD(NCoord),KEHESS(NCoord),KALHES(NCoord)
!
      DO i = 1, NCoord
         EGRAD(i) = optinfo%GradMol(i)
      ENDDO  
      ji = 1
      DO i = 1, NCoord
         DO j = 1, NCoord
            KEHESS(ji) = optinfo%HessMol(j,i) 
            KALHES(ji) = optinfo%HessMol(j,i) 
            ji = ji + 1
         ENDDO
      ENDDO
!
End Subroutine LS_copyGH
!==========================!
!  Atom_ini                !
!==========================!
      SUBROUTINE ATOM_INI(ATMARR,Molecule,optinfo,NAtoms,BOHR,lupri)
use precision
use optimization_input
use fundamental
use molecule_type
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri,NAtoms,i
      Type(MOLECULEINFO) :: Molecule
      Type(opt_setting) :: optinfo
      Real(realk) ATMARR(NAtoms,8)
      Real(realk) Coord(3,NAtoms)
      LOGICAL BOHR
!      Real(realk) :: fac
!
!     We initialize the ATMARR array. The first index runs over all
!     atoms, the second marks the following properties:
!
!               1 - Element number
!               2 - X coordinate of atom
!               3 - Y coordinate of atom
!               4 - Z coordinate of atom
!               5 - Covalent radius
!     
      FAC = bohr_to_angstrom
      IF (BOHR) FAC = 1.0E0_realk
!
   Do i = 1, NAtoms
      ! Atomic number/charge
      ATMARR(i,1) = Molecule%Atom(i)%Charge
      ! Coordinates
      ATMARR(i,2) = FAC*optinfo%Coordinates(1,i)
      ATMARR(i,3) = FAC*optinfo%Coordinates(2,i)
      ATMARR(i,4) = FAC*optinfo%Coordinates(3,i)
      ! Covalent radius
      ATMARR(i,5) = CovRad(NINT(ATMARR(i,1)),lupri)
   Enddo     
End subroutine ATOM_INI
!==================!
!    Add_Redspa    !
!==================!
! Combines a model internal Hessian with
! the calculated in reduced space
Subroutine Add_Redspa(optinfo,MaxCoord,Hess_Dim,Red_Hess,HesInt,lupri)
use precision
use optimization_input
use ls_util 
Implicit none
Integer :: Hess_Dim,MaxCoord ! Order of reduced and full Hessians
Real(realk) Red_Hess(Hess_Dim,Hess_Dim) ! Reduced Hessian
Real(realk) HesInt(MaxCoord,MaxCoord) ! Full Hessian
Type(opt_setting) :: optinfo
Integer :: lupri
Integer :: i,j,k
!
Do i = 1, optinfo%NIntCoord
   Do j = 1, Hess_dim
      If (optinfo%Red_space(j) .EQ. i) then
         Do k = 1, Hess_dim
            HesInt(i,optinfo%Red_Space(k)) = Red_Hess(j,k)
         Enddo
      Endif
   Enddo
Enddo
!
call lsheader(lupri,'Internal Hessian: numerical reduced, combined with model')
call ls_output(HESINT,1,optinfo%NIntCoord,1,optinfo%NIntCoord, &
     & MaxCoord,MaxCoord,1,LUPRI)
!
End subroutine Add_Redspa
!=======================!
! Rotate_coordinates    !  
!=======================!
! Translates and rotates Cartesians
! meant for scanning PES
Subroutine Rotate_coordinates(NAtoms,Coordinates,Origin_atom,Active_Atom)
!
Use precision
Implicit none
Integer :: NAtoms, Origin_atom, Active_atom, i
Real(realk), dimension(3,NAtoms) :: Coordinates
Real(realk), dimension(3) :: origin,Axis,normal
Real(realk), dimension(3,3) :: R   ! Rotation matrix
Real(realk) :: distance,theta
! Move molecule so that active atom is the origin
origin = Coordinates(:,Origin_atom)
Call Move_molecule(NAtoms,Coordinates,-origin)
! Find distance between Origin_Atom and Active_Atom
distance = SQRT( dot_product(Coordinates(:,Active_Atom), &
 & Coordinates(:,Active_Atom)) )
! Find the rotation angle
theta = acos(Coordinates(1,Active_Atom)/distance)
! Find the normal to the plane of rotation
normal(1) = 0.00000D0
normal(2) = distance*Coordinates(3,Active_Atom)
normal(3) = -distance*Coordinates(2,Active_Atom)
! Normalize normal
normal(2) = normal(2)/(SQRT(dot_product(normal,normal))) 
normal(3) = normal(3)/(SQRT(dot_product(normal,normal))) 
! Find the rotation matrix
Call Rotation_axis_angle(normal,theta,R)
! Rotate coordinates
Do i = 1,NAtoms
   Coordinates(:,i) = MATMUL(R,Coordinates(:,i))
Enddo
!
End subroutine Rotate_coordinates
!===================
!   Move_Molecule
!===================
Subroutine Move_Molecule(NAtoms, Coordinates, MoveVec)
!
! Translate molecule
!
Use precision
  Implicit None
  Real(realk), Dimension(3,NAtoms) :: Coordinates
  Real(realk), Dimension(3)        :: MoveVec
  Integer                               :: NAtoms, I
!
! Adds specified vector to all the atomic coordinates
!
  Do I = 1, NAtoms
    Coordinates(1,i) = Coordinates(1,i) + MoveVec(1)
    Coordinates(2,i) = Coordinates(2,i) + MoveVec(2)
    Coordinates(3,i) = Coordinates(3,i) + MoveVec(3)
  End Do
End Subroutine Move_Molecule
!======================!
! Rotation_axis_angle  !
!======================!
! Finds a 3 by 3 rotation matrix from
! axis and angle
Subroutine Rotation_axis_angle(normal,angle,Rotation)
Use precision
Implicit none
Real(realk), dimension(3,3) :: Rotation
Real(realk) :: angle
Real(realk), dimension(3) :: normal
!
Rotation(1,1) = cos(angle) + normal(1)**2*(1-cos(angle))
Rotation(1,2) = normal(1)*normal(2)*(1-cos(angle)) - normal(3)*sin(angle)
Rotation(1,3) = normal(1)*normal(3)*(1-cos(angle)) + normal(2)*sin(angle)
Rotation(2,1) = normal(2)*normal(1)*(1-cos(angle)) + normal(3)*sin(angle)
Rotation(2,2) = cos(angle) + normal(2)**2*(1-cos(angle))
Rotation(2,3) = normal(2)*normal(3)*(1-cos(angle)) - normal(1)*sin(angle)
Rotation(3,1) = normal(3)*normal(1)*(1-cos(angle)) - normal(2)*sin(angle)
Rotation(3,2) = normal(3)*normal(2)*(1-cos(angle)) + normal(1)*sin(angle)
Rotation(3,3) = cos(angle) + normal(3)**2*(1-cos(angle))
!
End subroutine Rotation_axis_angle
!=============!
! FM_Energy   !
!=============!
! Adds terms from extermnal force
! to the energy
Subroutine FM_energy(E,optinfo)
Use precision
use optimization_input
Implicit none
Type(opt_setting) :: optinfo
Real(realk) :: R_a(3),R_b(3),direction(3)
Real(realk) :: E
! Define the direction
R_a = optinfo%Coordinates(:,optinfo%Att_atom(1))
R_b = optinfo%Coordinates(:,optinfo%Att_atom(2))
direction = (R_b - R_a)/(sqrt(dot_product(R_b-R_a,R_b-R_a)))
!
Write(*,*)'Non-modified   E=',E
E = E + dot_product(optinfo%Ext_force*direction,R_a) - &
& dot_product(optinfo%Ext_force*direction,R_b)
Write(*,*)'Force-modified E=',E
!
End subroutine FM_energy












