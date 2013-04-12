!==============
!   Accuracy
!==============
Module Accuracy
  Integer, Parameter :: IKind = Selected_Real_Kind(P=15)
End Module Accuracy

!=============
!   Dummies
!=============
Module Dummies
  use precision
  real(realk), Parameter      :: Dummy = -999999E0_realk
  Integer, Parameter          :: IDummy = -999999
End Module Dummies
!====================
!   Dynamics_Param
!====================
Module Dynamics_Param
  use precision
  real(realk)                           :: MaxTime, RotTemp, PathLen, &
                                           StepLen, Temp, CurrTime, VibTemp, &
                                           InitialEnergy, CurrKinetic, &
                                           CurrPotential, CurrEnergy, &
                                           PrevPotential, PrdPotential, &
                                           TimeStep, TotalMass, TrajTime, &
                                           FragDist, CMDist, FragSize, &
                                           IntStepSize, NumIntAccuracy, &
                                           LowerRandom, UpperRandom
  real(realk), Dimension(3)        :: InitialAngMom, CurrAngMom, &
                                           MomInertia
  real(realk), Dimension(3,3)      :: PrincipalMom
  real(realk), Dimension(24)       :: DistLim
  Integer, Dimension(24,0:2)            :: DistCrit
  Integer                               :: FstTraj, MaxIter, &
                                           NumAtoms, NumCoord, NumPaths, &
                                           NUpdates, MaxUpdates, &
                                           HessUpdMethod, NumIntegScheme, &
                                           PrintGnrl, PrintLvl, RotSampMeth, &
                                           StepNum, TotTraj, Trajec, &
                                           VibSampMeth, VibModes, ZeroModes, &
                                           NegModes, MaxNumSteps, TraRotDim, &
                                           RanSeed, MaxRanCnt, NumDistCrit, &
                                           GetRandom, TotCharge
  Logical                               :: ProcessedInput = .False.
  Logical                               :: GradientOnly, InputVeloc, &
                                           TimeL, StepL, PathL, &
                                           VerletAlg, QuadraticMod, &
                                           SecondOrd, FifthOrd, NEW5OR, ThirdOrd,&
                                           OldAlgorithm, BondAnalysis, &
                                           DoProj, SampleOrient, &
                                           NoMove,FOCKMD
End Module Dynamics_Param
!============================
!   Integrators_and_Models
!============================
Module Integrators_and_Models
  Integer, Parameter :: NumSchemes = 5, NumModels = 2
  Character(Len = 24), Dimension(NumSchemes) :: SchemeName = &
      (/ 'Euler''s method          ','Runge-Kutta 2nd order   ', &
         'Runge-Kutta 4th order   ', 'Cash-Karp Runge-Kutta   ', &
         'Bulirsch-Stoer scheme   ' /)
  Character(Len = 24), Dimension(NumModels)  :: ModelName  = &
      (/ 'Quadratic model         ' , &
         'Fifth order model       ' /)
End Module Integrators_and_Models

!==============================
!   Hessian_Updating_Schemes
!==============================
Module Hessian_Updating_Schemes
  Integer, Parameter :: NumMeth = 5
  Character(Len = 30), Dimension(NumMeth) :: UpdateMethod = &
      (/ 'Symmetric Rank One (SR1)      ', &
         'Powell-symmetric-Broyden (PSB)', &
         'Bofill''s update (SR1 & PSB)   ', &
         'Sqrt(Bofill''s update)         ', &
         'Static Hessian (no updates)   ' /)
End Module Hessian_Updating_Schemes
!==============================!
!  Polynomials                 !
!==============================!
! Contains functions for polynomials
! and their derivatives for predictor-corrector
!
Module Polynomials
  use precision
    Implicit None
  Interface y1
    Module Procedure y1
  End Interface

  Interface y2
    Module Procedure y2
  End Interface

  Interface y3
    Module Procedure y3
  End Interface

  Interface y4
    Module Procedure y4
  End Interface

  Interface y5
    Module Procedure y5
  End Interface

  Interface y6
    Module Procedure y6
  End Interface

  Interface DerY1
    Module Procedure DerY1
  End Interface

  Interface DerY2
    Module Procedure DerY2
  End Interface

  Interface DerY3
    Module Procedure DerY3
  End Interface

  Interface DerY4
    Module Procedure DerY4
  End Interface

  Interface DerY5
    Module Procedure DerY5
  End Interface

  Interface DerY6
    Module Procedure DerY6
  End Interface

Contains
!
Function y1(s)
Implicit none
real(realk) :: y1, s
  y1 = 0E0_realk
  y1 = 1E0_realk - 10E0_realk*s**3 + 15E0_realk*s**4 - 6E0_realk*s**5
End function y1
!
Function y2(s,NormStep)
Implicit none
real(realk) :: y2, s, NormStep
  y2 = 0E0_realk
  y2 = NormStep*(s - 6E0_realk*s**3 + 8E0_realk*s**4 - 3E0_realk*s**5)
End function y2
!
Function y3(s,NormStep)
Implicit none
real(realk) :: y3, s, NormStep
  y3 = 0E0_realk
  y3 = (0.5E0_realk*NormStep**2)*(s**2 - 3E0_realk*s**3 + 3E0_realk*s**4 - s**5)
End function y3
!
Function y4(s)
Implicit none
real(realk) :: y4, s
  y4 = 0E0_realk
  y4 = 10E0_realk*s**3 - 15E0_realk*s**4 + 6E0_realk*s**5
End function y4
!
Function y5(s,NormStep)
Implicit none
real(realk) :: y5, s, NormStep
  y5 = 0E0_realk
  y5 = NormStep*(-4E0_realk*s**3 + 7E0_realk*s**4 - 3E0_realk*s**5)
End function y5
!
Function y6(s,NormStep)
Implicit none
real(realk) :: y6, s, NormStep
  y6 = 0E0_realk
  y6 = (0.5E0_realk*NormStep**2)*( s**3 - 2E0_realk*s**4 + s**5)
End function y6
!
Function DerY1(s,NormStep)
Implicit none
real(realk) :: DerY1,s,NormStep
  DerY1 = 0E0_realk
  DerY1 = NormStep**(-1)*(-30E0_realk*s**2 + 60E0_realk*s**3 - 30E0_realk*s**4)
End function DerY1
!
Function DerY2(s)
Implicit none
real(realk) :: DerY2,s
  DerY2 = 0E0_realk
  DerY2 = 1E0_realk - 18E0_realk*s**2 + 32E0_realk*s**3 - 15E0_realk*s**4
End function DerY2
!
Function DerY3(s,NormStep)
Implicit none
real(realk) :: DerY3,s,NormStep
  DerY3 = 0E0_realk
  DerY3 = 0.5E0_realk*NormStep*(2E0_realk*s - 9E0_realk*s**2 + 12E0_realk*s**3 - 5E0_realk*s**4)
End function DerY3
!
Function DerY4(s,NormStep)
Implicit none
real(realk) :: DerY4,s,NormStep
  DerY4 = 0E0_realk
  DerY4 = NormStep**(-1)*(30E0_realk*s**2 - 60E0_realk*s**3 + 30E0_realk*s**4)
End function DerY4
!
Function DerY5(s)
Implicit none
real(realk) :: DerY5,s
  DerY5 = 0E0_realk
  DerY5 = - 12E0_realk*s**2 + 28E0_realk*s**3 - 15E0_realk*s**4
End function DerY5
!
Function DerY6(s,NormStep)
Implicit none
real(realk) :: DerY6,s,NormStep
  DerY6 = 0E0_realk
  DerY6 = 0.5E0_realk*NormStep*(3E0_realk*s**2 - 8E0_realk*s**3 + 5E0_realk*s**4)
End function DerY6
!
End module Polynomials
!==============================!
!  Fock_mat_var                !
!==============================!
Module Fock_mat_var
use precision
! Sets variables to read Fock matrix from SIRIFC file and contains FMD variables
Implicit none
Integer :: FMDim,NPoints,FMADim,PolyOrd,N_FMD 
Logical :: WINDOW, Const_TS
!
Interface Set_Var
   Module procedure Set_Var
End interface !Set_Var
!
Contains
        !
	Subroutine Set_Var(PolyOrd,NPoints,WINDOW,FMADim,Const_TS,lupri)
        ! Sets the size of Fock_matrix_array
        !
        Implicit none
        Integer :: PolyOrd, NPoints, FMADim,lupri
        Logical :: WINDOW, Const_TS
        !
        ! Checking whether NPoints > Polyord. If not - quit!
        !
        If (NPoints .LE. Polyord) then
            Call LSQuit('Number of data-points is less &
            &than polynomial order in Fock matrix dynamics!',lupri)
        Endif
        !
        If ( (WINDOW .EQV. .TRUE.) .AND. (Const_TS .EQV. .FALSE.)) then
           FMADim = NPoints + 4
        Else
           FMADim = NPoints + 1
        Endif        
        End subroutine Set_Var
        !
End module Fock_mat_var
!==============================!
! SCF_convergence              !
!==============================!
Module SCF_convergence
use precision
! Contains variables for calculation of # of everage SCF cycles in trajectory
IMPLICIT NONE
INTEGER :: ITER, TOT_ITER
real(realk) :: AV_ITER
End module SCF_convergence
!=====================!
! Dynamics_time       !
!=====================!
Module Dynamics_time
use precision
! Time of trajectory integration
IMPLICIT NONE
real(realk) :: Start_Time,Stop_Time
End module Dynamics_time
