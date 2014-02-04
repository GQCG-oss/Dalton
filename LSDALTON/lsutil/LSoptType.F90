!==============================================!
! Type containing configuration for optimizer  !
!==============================================!
Module optimization_type
Use precision
Type opt_setting
     !
     Logical :: optimize
     Integer :: MaxIter
     Integer :: scanMaxIter
     Real(realk) :: Change
     Integer :: ChGradThr
     Integer :: IcnLevel
     Integer :: IpDef
     Integer :: NSPMod
     Logical :: Saddle
     Logical :: NoTrust
     Integer :: ItBreak
     Logical :: NoHessianWrite
     Logical :: DoSpE
     Logical :: DoPre
     Logical :: FinPre
     Logical :: KeepHessian
     Logical :: RejIni
     Logical :: GradIni
     Logical :: ChangeGrad
     Logical :: ConOpt
     Logical :: RemCoord
     Logical :: AddCoord
     Logical :: Rebid
     Integer :: NumPre
     Integer :: Pre
     Integer :: TotRj
     Integer :: Redic
     Integer :: NIntCoord
     Integer :: ICartCoord
     Integer :: Stblz
     Logical :: SteepDescend
     Logical :: RanKon
     Logical :: PSB
     Logical :: DFP
     Logical :: BFGS
     Logical :: Bofill
     Logical :: BFGSR1
     Logical :: Multi
     Logical :: Schleg
     Logical :: Newton
     Logical :: QuadSD
     Logical :: FirstOrd
     Logical :: SecondOrd
     Logical :: DelInt
     Logical :: RedInt
     Logical :: CartCoord
     Logical :: InrdHess
     Logical :: InitHess
     Logical :: ModHes
     Logical :: CMBMod
     Logical :: InmdHess
     Logical :: HessFile
     Real(realk) :: EvLini
     Logical :: FindRe
     Logical :: Visual
     Logical :: VRLM
     Logical :: VrBond
     Logical :: Vreigv
     Logical :: VrCord
     Logical :: VrViba
     Logical :: TrustCh
     Logical :: TrustFC
     Integer :: IPrint
     Real(realk) :: TrustRad 
     Real(realk) :: TrustIn
     Real(realk) :: TrustDe
     Real(realk) :: RTENBD
     Real(realk) :: RTENGD
     Real(realk) :: RTRJMN
     Real(realk) :: RTRJMX
     Real(realk) :: PRVRMS
     Real(realk) :: PRVMAX
     Logical :: LnSearch
     Logical :: RatFun
     Logical :: TrustRg
     Logical :: GDIIS
     Logical :: GeConv
     Logical :: Baker
     Logical :: NoAux
     Logical :: NodiHe
     Logical :: NoAdda
     Logical :: RedRed
     Logical :: LINDHD
     Integer :: ItrNmr
     Integer :: MaxRej
     Real(realk) :: Displa
     Real(realk) :: GradThr
     Real(realk) :: ThrStep
     Real(realk) :: ThrErg
     Logical :: NatNorm
     Logical :: GeoAll
     Logical :: Restart
     Integer :: Condi
     Integer :: NCon
     Integer :: NcoordTot
     Integer :: NTempMat ! Number of temporary matrices 
     Integer :: ICon
     Integer :: IConV
     Integer :: ChgThresh
     Integer :: IterFreeze
     Integer :: NAdd
     Integer :: NRemove
     Integer :: IRemove 
     Integer :: nProjected 
     Integer :: KEPTIT 
     Integer, pointer :: AddCoordArray(:,:)
     Integer, pointer :: FreezeArray(:)
     Integer, pointer :: IConstr(:)
     Integer :: NFreeze
     Integer :: NFrAtom
     Character(len=80), pointer :: PreText(:)
     Character(len=80) :: SpBText
     Logical :: doNumGradGeomOpt !To run geom. optim. with a numerical gradient
     Real(realk) :: findif_mesh  !mesh with a numerical gradient 
! Serve as parameters
     Real(realk) :: ZeroGrad ! A threshold for zero gradient 
     Real(realk) :: DefTh2 
     Real(realk) :: DefDsp 
     Real(realk) :: DefThe 
     Integer :: MaxReject 
     Integer :: INDHES  ! Hessian index
! 
     Real(realk), pointer :: STPDIA(:),STPSYM(:),STPINT(:), &
     & GRDDIA(:),EVAL(:),EVALOL(:),GRDINT(:) 
     Real(realk) CNDHES(0:7)
     Real(realk), pointer:: CoordInt(:)
     Real(realk), pointer :: Coordinates(:,:)
     Integer, pointer :: INTCRD(:,:)
     Integer, dimension(0:5) :: ICONF
     Real(realk), pointer :: GradMol(:)
     Real(realk), pointer :: HessMol(:,:)
!
!    Norms of step and molecular gradient
     Real(realk) :: StepNorm
     Real(realk) :: GradNorm
!
     Real(realk) :: ThGradMax
     Real(realk) :: ThStepMax

     Real(realk) :: energy
     Real(realk) :: energyOld
     Real(realk) :: predictedChange
!    Everything about numerical internal Hessian
     Logical :: RedSpa
     Logical :: CartRS
     Integer :: Hess_dim
     Integer :: Red_Atoms
     Integer, pointer :: Cart_Red_Space(:)
     Integer, pointer :: Red_space(:)
     Real(realk) :: displ
     Logical :: ForBac
     INTEGER :: numsteps     !number of geometry steps as input
     INTEGER :: stepnum      !current step number 
     INTEGER :: step_num     ! Needed for standalone 
!    Scanning the internal
     LOGICAL :: simple_scan
     Integer :: ScanCoord
     Real(realk) :: Scan_step  

!    Dynamical threshold (for DEC/FOT)
     Logical :: dynopt
     Logical :: dynamicThreshold
     Logical :: dynamicConvergence
     Logical :: dynamicChange
! 
! Scan_info contains information about scanning an internal coordinate:
!           1 - value of the coordinate
!           2 - energy at this point
    Real(realk), pointer :: Scan_info(:,:)
! Needed for scanning the internals
    Integer, pointer :: Atoms_to_move(:,:)
    Integer, dimension(2) :: N_to_Move
    Logical :: New_stepping
    Logical :: IBT
    Logical :: OLDIBT
    Logical :: Shanks
    Integer :: Deriv_order
End type opt_setting

contains

!Added to avoid "has no symbols" linking warning
subroutine optimization_type_void()
end subroutine optimization_type_void

End module optimization_type
