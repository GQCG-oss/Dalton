!=====================================!
! Geometry optimizer for LSDALTON     !
!=====================================!
!
!==============================================!
! Type containing configuration for optimizer  !
!==============================================!
Module optimization_input
use optimization_type
Use precision
Use ls_util
Use matrix_module
Use memory_handling
!    Serve for size of arrays
     Integer :: MXCENT 
     Integer  MXCOOR 
     real(realk),external :: ddot
     real(realk),external :: dnorm2
     real(realk),external :: vecang_ls
     real(realk),external :: vdwrad_ls
!     real(realk),external :: dv3dot_ls

! This is not needed when you only have 1 sub with same name. TK
!Interface Optimization_set_defaul_config
!  Module procedure Optimization_set_default_config
!End Interface
!
!Interface LS_optimization_input
!  Module procedure LS_optimization_input
!End Interface
!
!Interface Final_opt_config
!  Module procedure Final_opt_config
!End Interface
Contains
!=================================!
! Optimization_set_default_config !
!=================================!
Subroutine Optimization_set_default_config(optinfo)
Implicit none
Type(opt_setting) :: optinfo
! Serve as parameters, do not change
     optinfo%MaxReject = 3
! Setting various thresholds
     optinfo%DefTh2 = 1.0E-5_realk
     optinfo%DefDsp = 1.0E-3_realk
     optinfo%DefThe = 1.0E-6_realk
     optinfo%ZeroGrad = 1.0E-7_realk
!    Norms of step and molecular gradient
     optinfo%StepNorm = 0E0_realk
     optinfo%GradNorm = 0E0_realk
! 
     optinfo%Restart = .FALSE.
! Total number of coordinates
     optinfo%NCoordTot = 0
!
     optinfo%NTempMat = 0
     optinfo%KEPTIT = 0
     optinfo%optimize = .FALSE.
     optinfo%maxIter     = 100
     optinfo%scanMaxIter = 100
     optinfo% Change  = 0E0_realk
     optinfo%ChGradThr    = 0
     optinfo%IcnLevel = -1
     optinfo%IpDef    = 0
     optinfo%NSPMod   = -1
     optinfo%Saddle   = .FALSE. 
     optinfo%NoTrust  = .FALSE.
     optinfo%ItBreak  = -1
     optinfo%NoHessianWrite = .TRUE.
     optinfo%DoSpE = .FALSE.
     optinfo%DoPre = .FALSE.
     optinfo%FinPre = .FALSE.
     optinfo%KeepHessian = .FALSE.
     optinfo%RejIni = .FALSE.
     optinfo%GradIni = .FALSE.
     optinfo%ChangeGrad = .FALSE.
     optinfo%ConOpt = .FALSE.
     optinfo%RemCoord = .FALSE.
     optinfo%AddCoord = .FALSE.
     optinfo%Rebid = .FALSE.
     optinfo%NumPre = 0
     optinfo%Pre = 0
     optinfo%TotRj = 0 
     optinfo%Redic = -1
     optinfo%NIntCoord = 0
     optinfo%Stblz = -1
     optinfo%SteepDescend = .FALSE.
     optinfo%RanKon = .FALSE.
     optinfo%PSB = .FALSE.
     optinfo%DFP = .FALSE.
     optinfo%BFGS = .FALSE.
     optinfo%Bofill = .FALSE.
     optinfo%BFGSR1 = .FALSE.
     optinfo%Multi = .FALSE.
     optinfo%Schleg = .FALSE.
     optinfo%Newton = .FALSE.
     optinfo%QuadSD = .FALSE.
     optinfo%FirstOrd = .FALSE.
     optinfo%SecondOrd = .FALSE.
     optinfo%DelInt = .FALSE.
     optinfo%RedInt = .TRUE.
     optinfo%CartCoord = .FALSE.
     optinfo%InrdHess = .FALSE.
     optinfo%InitHess = .FALSE.
     optinfo%ModHes = .FALSE.
     optinfo%CMBMod = .FALSE.
     optinfo%InmdHess = .FALSE.
     optinfo%HessFile = .FALSE.
     optinfo%EvLini = -1E0_realk
     optinfo%FindRe = .FALSE.
     optinfo%Visual = .FALSE.
     optinfo%VRLM = .FALSE.
     optinfo%VrBond = .FALSE.
     optinfo%Vreigv = .FALSE.
     optinfo%VrCord = .FALSE.
     optinfo%VrViba = .FALSE.
     optinfo%TrustCh = .FALSE.
     optinfo%TrustFC = .FALSE.
     optinfo%IPrint = 0
     optinfo%TrustRad = 0.5E0_realk
     optinfo%TrustIn = 1.2E0_realk
     optinfo%TrustDe = 0.7E0_realk
     optinfo%RTENBD = 0.4E0_realk
     optinfo%RTENGD = 0.8E0_realk
     optinfo%RTRJMN = -0.1E0_realk
     optinfo%RTRJMX = 3.0E0_realk
     optinfo%PRVRMS = 0.0E0_realk
     optinfo%PRVMAX = 0.0E0_realk
     optinfo%LnSearch = .FALSE.
     optinfo%RatFun = .FALSE.
     optinfo%TrustRg = .FALSE.
     optinfo%GDIIS = .FALSE.
     optinfo%GeConv = .FALSE.
     optinfo%Baker = .FALSE.
     optinfo%NoAux = .FALSE.
     optinfo%NodiHe = .FALSE.
     optinfo%NoAdda = .FALSE.
     optinfo%RedRed = .FALSE.
     optinfo%LINDHD = .TRUE.
     optinfo%ItrNmr = 0
     optinfo%MaxRej = 3
     optinfo%GradThr = optinfo%DefTH2
     optinfo%ThrStep = optinfo%DefTH2
     optinfo%ThrErg =  optinfo%DefTHE
     optinfo%NatNorm = .TRUE.
     optinfo%GeoAll = .FALSE.
     optinfo%Condi = 2
     optinfo%NCon = 0
     optinfo%IterFreeze = -1
     optinfo%NAdd = 0
     optinfo%NRemove = 0
     optinfo%IRemove = 0
     optinfo%nProjected  = 0
     optinfo%NFreeze = 0
     optinfo%energy           = 0E0_realk
     optinfo%energyOld        = 0E0_realk
     optinfo%predictedChange  = 0E0_realk
!    Numerical internal Hessian in the reduced space
     optinfo%RedSpa = .FALSE.
     optinfo%CartRS = .FALSE.
     optinfo%ForBac = .FALSE.
     optinfo%displ = 1E-03_realk
     optinfo%numsteps = 0          !number of geometry steps as input
     optinfo%stepnum = 0           !current step number 
     optinfo%IndHes = 0
!    Scanning along coordinates
     optinfo%simple_scan = .FALSE.
     optinfo%ScanCoord = -1
     optinfo%scan_step = 0.1D0
     optinfo%dynopt = .FALSE.
     optinfo%dynamicThreshold   = .FALSE.
     optinfo%dynamicConvergence = .FALSE.
     optinfo%doNumGradGeomOpt = .FALSE.
     optinfo%findif_mesh = 1.0D-5

!    Non-iterative stepping scheme based on perturbation theory (Vladimir R. and Ulf E.)
     optinfo%New_stepping = .TRUE.
     optinfo%IBT = .FALSE.
     optinfo%OLDIBT = .FALSE.
     optinfo%Shanks = .FALSE.
     optinfo%Deriv_order = 20
!    Force_modified PES
     optinfo%FMPES = .FALSE.
     optinfo%Ext_force = 0E0_realk
     optinfo%Att_atom = 0
! 
End subroutine Optimization_set_default_config 
!=========================!
! LS_optimization_input   !
!=========================!
Subroutine LS_optimization_input(optinfo,readword,keyword,lucmd,lupri,NAtoms)
!Use ls_util
!Use memory_handling
! Input processing for linear scaling optimization
Implicit none
Type(opt_setting) :: optinfo
Logical :: readword
Integer :: lucmd, lupri  ! File units
Integer :: NAtoms,i
Character(Len = 70) :: Keyword
Character(Len = 1) :: Prompt
Integer :: FileStatus
Character(Len=7), dimension(88) :: KwordTABLE = &
          (/'.PRINT ', '.MAX IT', '.TRUSTR', '.TR FAC', &
            '.TR LIM', '.MAX RE', '.NOTRUS', '.ENERGY', &
            '.GRADIE', '.STEP T', '.CONDIT', '.NOBREA', &
            '.SP BAS', '.DYNOPT', '.VISUAL', '.VRML  ', & 
            '.SYMTHR', '.TRSTRG', '.VR-BON', '.VR-EIG', &
            '.INITHE', '.INITEV', '.HESFIL', '.REJINI', &
            '.STEEPD', '.RANKON', '.PSB   ', '.DFP   ', &
            '.BFGS  ', '.NEWTON', '.QUADSD', '.SCHLEG', &
            '.HELLMA', '.BAKER ', '.M-BFGS', '.CARTES', &
            '.REDINT', '.INIRED', '.1STORD', '.2NDORD', &
            '.GRDINI', '.DISPLA', '.CONSTR', '.MODHES', &
            '.REMOVE', '.INIMOD', '.FINDRE', '.CMBMOD', &
            '.RF    ', '.GDIIS ', '.DELINT', '.NODIHE', &
            '.VR-COR', '.VR-VIB', '.VR-SYM', '.M-PSB ', &
            '.LINE S', '.SADDLE', '.MODE  ', '.BOFILL', &
            '.NOAUX ', '.BFGSR1', '.STABIL', '.GEOANA', &
            '.ADDCRD', '.NOADDA', '.RREDUN', '.NATNRM', &
            '.VLOOSE', '.LOOSE ', '.TIGHT ', '.VTIGHT', &
            '.NOHSWR', '.FREEZE', '.FRZITR', '.REDSPA', &
            '.CARTRS', '.FORBAC', '.SCANSI', '.SCANST', &
            '.NUMOPT', '.NUMESH', '.NOHOPE', '.ITERBT', &
            '.DERORD', '.OLDIBT', '.SHANKS', '.FMPES '/)
! Number of cartesian coordinates
optinfo%IcartCoord = NAtoms*3
!
! Allocate optinfo%IConstr anyway
!
Call mem_alloc(optinfo%IConstr,30*NAtoms)
optinfo%IConstr = 0
!
! Read keywords
!
Do
      Read(lucmd,'(A70)', IOStat = FileStatus) Keyword
      If (FileStatus > 0) Call LSQuit('Error reading lucmd',lupri)
      If ((FileStatus < 0) ) Exit
      Prompt = Keyword(1:1)
      If ((Prompt == '#') .or. (Prompt == '!')) Cycle
      If (Prompt == '*') then
         ReadWord = .FALSE.
         Exit
      Endif
        If (Any(KWordTable(1:size(KWordTable)) == Keyword(1:7))) Then
          Select Case (Keyword)
                 Case('.PRINT ')
                   Read(lucmd,*) optinfo%IPrint
                 Case('.MAX IT')
                   Read(lucmd,*) optinfo%MaxIter
                 Case('.TRUSTR')
                   Read(lucmd,*) optinfo%TrustRad 
                   optinfo%TrustCh = .TRUE.             
                 Case('.TR FAC') 
                   Read(lucmd,*) optinfo%TrustIn,optinfo%TrustDe
                   optinfo%TrustCh = .TRUE.
                   optinfo%TrustFc = .TRUE.
                 Case('.TR LIM')
                   Read(lucmd,*) optinfo%RTENBD, &
                   optinfo%RTENGD,optinfo%RTRJMN,optinfo%RTRJMX
                   optinfo%TrustCh = .TRUE.
                 Case('.MAX RE')
                   Read (lucmd,*) optinfo%MaxRej
                 Case('.NOTRUS')
                   optinfo%NoTrust = .TRUE.
                 Case('.ENERGY')
                    Call lsquit('.ENERGY not implemented in LSDALTON',lupri)
!                   Read(lucmd,*) optinfo%ThrErg
!                   optinfo%ChangeGrad = .TRUE.
!                 Case('.GRADIE')
!                   Read(lucmd,*) optinfo%GradThr
!                   optinfo%ChangeGrad = .TRUE.
!                   optinfo%ChGradThr = optinfo%ChGradThr + 1
!                 Case('.STEP T')
!                   Read(lucmd,*) optinfo%ThrStep
!                   optinfo%ChangeGrad = .TRUE.
!                   optinfo%ChGradThr = optinfo%ChGradThr + 1
                 Case('.CONDIT')
                    Call lsquit('.CONDIT not implemented in LSDALTON',lupri)
!                   Read(lucmd,*) optinfo%Condi
                 Case('.NOBREA')
                   Write(lupri,*)'.NOBREA is a symmetry keyword'
                   Call LSQuit('Symmetry is not implemented in LSDALTON',lupri)
                 Case('.SP BAS')
                    Call lsquit('.SP BAS not implemented in LSDALTON',lupri)
!                   Read(lucmd,*) optinfo%SpBText
!                   optinfo%DoSpE = .TRUE.
                  Case('.DYNOPT') 
                     optinfo%dynopt = .TRUE.
                 Case('.PREOPT')
                    Call lsquit('.PREOPT not implemented in LSDALTON',lupri)
!                   Read(lucmd,*) optinfo%NumPre
!                   If ((optinfo%NumPre .LT. 1)) then
!                     Write(lupri,'(A,I2,A/)') &
!                     ' Number of preoptimization sets must be greater than 1!'
!                     Call LSQuit('Illegal number of preoptimization sets.',lupri)
!                   Else
!                     Call mem_alloc(optinfo%PreText,optinfo%NumPre)
!                     Do i = 1, optinfo%NumPre
!                        Read(lucmd,'(A60)') optinfo%PreText(i)
!                     Enddo
!                   Endif
                 Case('.VISUAL')
                    Call lsquit('.VISUAL not implemented in LSDALTON',lupri)
!                     optinfo%Visual = .TRUE.
                 Case('.VRLM')
                    Call lsquit('.VRLM not implemented in LSDALTON',lupri)
!                     optinfo%VRLM = .TRUE.
                 Case('.SYMTHR')
                   Write(lupri,*)'.SYMTHR is a symmetry keyword'
                   Call LSQuit('Symmetry is not implemented in LSDALTON',lupri)
                 Case('.TRSTRG')
                   optinfo%TrustRg = .TRUE.
!                 Case('.VR-BON')
!                 Case('.VR-EIG')
                 Case('.INITHE')
                   Write(lupri,*)'INITHE optinon is not available in LSDALTON'
                   Call LSQuit('Hessian not available in LSDALTON',lupri)
                 Case('.INITEV')
                   Read(lucmd,*) optinfo%EvLini
                 Case('.HESFIL')
                    Call lsquit('.HESFIL not implemented in LSDALTON',lupri)
!                   optinfo%HessFile = .TRUE.
                 Case('.REJINI')
                   optinfo%REJINI = .TRUE.
!                 Case('.STEEPD')
                 Case('.RANKON')
                    Call lsquit('.RANKON not implemented in LSDALTON',lupri)
!                   optinfo%RanKon = .TRUE.
                 Case('.PSB   ')
                   optinfo%PSB = .TRUE. 
                 Case('.DFP   ')
                   optinfo%DFP = .TRUE. 
                 Case('.BFGS  ')
                   optinfo%BFGS = .TRUE. 
                 Case('.NEWTON')
                   Write(lupri,*)'NEWTON optinon is not available in LSDALTON'
                   Call LSQuit('Hessian not available in LSDALTON',lupri)
                 Case('.QUADSD')
                    Call lsquit('.QUADSD not implemented in LSDALTON',lupri)
!                   optinfo%QuadSD = .TRUE.
!                 Case('.SCHLEG')
!                 Case('.HELLMA')
                 Case('.NUMOPT')
                    WRITE(LUPRI,*) 'Numerical Gradient geometry optimization calculations are carried out'
                    optinfo%doNumGradGeomOpt = .True.
                 Case('.NUMESH')
                    Read(lucmd,*) optinfo%findif_mesh
                 Case('.BAKER')
                   optinfo%Baker = .TRUE.
!                 Case('.M-BFGS')
                 Case('.CARTES')
                    optinfo%CartCoord = .TRUE. 
                    optinfo%RedInt = .FALSE.
                 Case('.REDINT')
                    If (optinfo%CartCoord) then
                       Call LSQuit('The user must choose&
                                 & between optimization in cartesian or&
                                 & internal coordinates', lupri)
                    Endif       
                    optinfo%RedInt = .TRUE.
                    optinfo%CartCoord = .FALSE.
                 Case('.INIRED')
                    optinfo%InrdHess = .TRUE.
                 Case('.1STORD')
                    optinfo%FirstOrd = .TRUE.
                 Case('.2NDORD')
                    Call lsquit('.2NDORD not implemented in LSDALTON',lupri)
!                    optinfo%SecondOrd = .TRUE.
                 Case('.GRDINI')
                    optinfo%GradIni = .TRUE.
                 Case('.CONSTR')
!                    Call lsquit('.CONSTR not implemented in LSDALTON',lupri)
                    Read(lucmd,*) optinfo%NCon
                    If (optinfo%NCon .LE. 0) then
                       Write(lupri,'(A,I2,A/)') &
                       'Number of constrained coordinates must be at least 1!'
                       Call LSQuit('Illegal number of constrained coordinates.',lupri)
                    Endif
                    Do i = 1,optinfo%NCon
                       Read(lucmd,*) optinfo%ICon
                       If (optinfo%ICon .LE. 0) then
                          Write(lupri,'(/,A,I2,A/)') &
                          'The constrained coordinate should be & 
                          &between 1 and ',(8*NAtoms),'!' 
                          Call LSQuit('Illegal constrained coordinate',lupri)
                       Endif 
                       optinfo%IConstr(optinfo%Icon) = 1
                     Enddo
                     optinfo%ConOpt = .TRUE.  
                 Case('.MODHES')
!                    Call lsquit('.MODHES not implemented in LSDALTON',lupri)
                     optinfo%ModHes = .TRUE.
                 Case('.INIMOD')
!                    Call lsquit('.INIMOD not implemented in LSDALTON',lupri)
                     optinfo%InmdHess = .TRUE.
                 Case('.FINDRE')
!                    Call lsquit('.FINDRE not implemented in LSDALTON',lupri)
                    optinfo%FindRe = .TRUE.
                    optinfo%Redint = .TRUE.
                 Case('.CMBMOD')
                    Call lsquit('.CMBMOD not implemented in LSDALTON',lupri)
!                    optinfo%CMBMod = .TRUE.
                 Case('.RF    ')
                    Call lsquit('.RF not implemented in LSDALTON',lupri)
!                    optinfo%RatFun = .TRUE.
                 Case('.GDIIS')
                    optinfo%GDIIS = .TRUE.
                 Case('.DELINT')
                    Call lsquit('.DELINT not implemented in LSDALTON',lupri)
!                    optinfo%DelInt = .TRUE.
                 Case('.NODIHE')
                    Call lsquit('.NODIHE not implemented in LSDALTON',lupri)
!                    optinfo%NodiHe = .TRUE.
!                 Case('.VR-COR')
!                 Case('.VR-VIB')
!                 Case('.VR-SYM')
                 Case('.M-PSB')
                    Call lsquit('.M-PSB not implemented in LSDALTON',lupri)
!                     optinfo%Multi = .TRUE.
!                     optinfo%PSB = .TRUE.
                 Case('.LINE S')
                    Call lsquit('.LINE S not implemented in LSDALTON',lupri)
!                    optinfo%LnSearch = .TRUE.
                 Case('.SADDLE')
                     optinfo%Saddle = .TRUE.
                 Case('.MODE  ')
                    Call lsquit('.MODE not implemented in LSDALTON',lupri)
!                    Read(lucmd,*)optinfo%NSPMod 
                 Case('.BOFILL')
                    optinfo%Bofill = .TRUE.
                 Case('.NOAUX ')
                    Call lsquit('.NOAUX not implemented in LSDALTON',lupri)
!                    optinfo%NoAux = .TRUE.
                 Case('.BFGSR1')
                    Call lsquit('.BFGSR1 not implemented in LSDALTON',lupri)
!                    optinfo%BFGSR1 = .TRUE.
                 Case('.STABIL')
                    Call lsquit('.STABIL not implemented in LSDALTON',lupri)
!                    Read(lucmd,*) optinfo%Stblz
                 Case('.GEOANA')
                    Call lsquit('.GEOANA not implemented in LSDALTON',lupri)
!                    optinfo%GeoAll = .TRUE.
                 ! Adding coordinates
                 Case('.ADDCRD')
                    Call lsquit('.ADDCRD not implemented in LSDALTON',lupri)
!                    Read(lucmd,*) optinfo%NAdd
!                    If ((optinfo%NAdd .LE. 0) .OR. (optinfo%NAdd .GT. 10)) then
!                       Write(lupri,'(/,A/)') &
!                       ' Number of coordinates to be added must be &
!                       &between 1 and 10'
!                       Call LSQuit(lupri,'Illegal number of coordinates to be added.')
!                    Endif
!                    Call mem_alloc(optinfo%AddCoordArray,optinfo%NAdd,4)
!                    optinfo%AddCoordArray = 0
!                    Do i = 1,optinfo%NAdd
!                       Read(lucmd,*) optinfo%AddCoordArray(i,1), &
!                       optinfo%AddCoordArray(i,2),optinfo%AddCoordArray(i,3),&
!                       optinfo%AddCoordArray(i,4)
!                    Enddo 
!                    optinfo%AddCoord = .TRUE. 
                 ! Removing coordinates
                 Case('.REMOVE')
                    Call lsquit('.REMOVE not implemented in LSDALTON',lupri)
!                    Read(lucmd,*) optinfo%NRemove
!                    If (optinfo%NRemove .LE. 0) then
!                       Write(lupri,'(A,I2,A/)') &
!                       'Number of removed coordinates must be at least 1!'
!                       Call LSQuit('Illegal number of constrained coordinates.',lupri)
!                    Endif
!                    Do i = 1,optinfo%NRemove
!                       Read(lucmd,*) optinfo%IRemove
!                       If (optinfo%IRemove .LE. 0) then
!                          Write(lupri,'(/,A,I2,A/)') &
!                          'The removed coordinate number should be & 
!                          &between 1 and ',(30*NAtoms),'!' 
!                          Call LSQuit('Illegal coordinate number',lupri)
!                       Endif 
!                       optinfo%IConstr(optinfo%IRemove) = 2
!                     Enddo
!                     optinfo%RemCoord = .TRUE.  
                 Case('.NOADDA')
                    Call lsquit('.NOADDA not implemented in LSDALTON',lupri)
!                    optinfo%NoAdda = .TRUE.
                 Case('.PREDUN')
                    Call lsquit('.PREDUN not implemented in LSDALTON',lupri)
!                    optinfo%RedRed = .TRUE.
                 Case('.NATNRM')
                    Call lsquit('.NATNRM not implemented in LSDALTON',lupri)
!                    optinfo%NatNorm = .TRUE.
                 Case('.VLOOSE')
                    optinfo%IcnLevel = 1
                    optinfo%ChgThresh = optinfo%ChgThresh + 1
                 Case('.LOOSE ')
                    optinfo%IcnLevel = 2
                    optinfo%ChgThresh = optinfo%ChgThresh + 1
                 Case('.TIGHT ')
                    optinfo%IcnLevel = 4
                    optinfo%ChgThresh = optinfo%ChgThresh + 1
                 Case('.VTIGHT')
                    optinfo%IcnLevel = 5
                    optinfo%ChgThresh = optinfo%ChgThresh + 1
                 Case('.NOHSWR')
                    Call lsquit('.NOHSWR not implemented in LSDALTON',lupri)
!                    optinfo%NoHessianWrite = .TRUE.
                 Case('.FREEZE')
                    Call lsquit('.FREEZE not implemented in LSDALTON',lupri)
!                    Read(lucmd,*) optinfo%NFreeze !Number of frozen atoms
!                    If ((optinfo%NFreeze .LE. 0) .OR. (optinfo%NAdd .GT. NAtoms)) then
!                       Write(lupri,'(/,A/)') &
!                       'Number of frozen atoms must be &
!                       &between 1 and', NAtoms,'!'
!                       Call LSQuit(lupri,'Illegal number of frozen atoms.')
!                    Endif
!                    Call mem_alloc(optinfo%FreezeArray,optinfo%NFreeze)
!                    optinfo%FreezeArray = 0
!                    Do i=1, optinfo%NFreeze
!                       Read(lucmd,*) optinfo%NFrAtom    ! Number of atom to be frozen
!                       If ((optinfo%NFrAtom .LE. 0) .OR. (optinfo%NFrAtom .GT. NAtoms)) then
!                          Write(lupri,'(/,A/)') &
!                         ' Number of frozen atom must be &
!                          between 1 and', NAtoms,'!'
!                         Call LSQuit(lupri,'Illegal frozen atom.')
!                       Endif
!                      optinfo%FreezeArray(i) = optinfo%NFrAtom 
!                    Enddo
                 Case('.FRZITR')
                    Call lsquit('.FRIZTR not implemented in LSDALTON',lupri)
!                    Read(lucmd,*) optinfo%IterFreeze
                 Case('.REDSPA')
                     optinfo%RedSpa = .TRUE.
                     Read(lucmd,*) optinfo%Hess_dim                              
                     If (optinfo%Hess_dim .LE. 0) then
                        Write(lupri,'(A,I2,A/)') &
                        'Number of  coordinates in reduced space  must be at least 1!'
                        Call LSQuit('Illegal number of coordinates in reduced space.',lupri)
                     Endif
                     Call mem_alloc(optinfo%Red_space,optinfo%Hess_dim)
                     Do i = 1,optinfo%Hess_dim
                        Read(lucmd,*) optinfo%Red_space(i)
                        If (optinfo%Red_Space(i) .LE. 0) then
                           Write(lupri,'(/,A,I2,A/)') &
                           'The coordinate in reduced space number should be & 
                           &between 1 and ',(30*NAtoms),'!' 
                           Call mem_dealloc(optinfo%Red_space)
                           Call LSQuit('Illegal coordinate number',lupri)
                        Endif 
                     Enddo
                 Case('.CARTRS')
                     optinfo%CartRS = .TRUE.
                     Read(lucmd,*) optinfo%Red_Atoms                              
                     If (optinfo%Red_Atoms .LE. 0) then
                        Write(lupri,'(A,I2,A/)') &
                        'Number of atoms in reduced space  must be at least 1!'
                        Call LSQuit('Illegal number of atoms in reduced space.',lupri)
                     Endif
                     Call mem_alloc(optinfo%Cart_Red_space,optinfo%Red_Atoms)
                     Do i = 1,optinfo%Red_Atoms
                        Read(lucmd,*) optinfo%Cart_Red_Space(i)
                        If ((optinfo%Cart_Red_Space(i) .LE. 0) .OR. &
                        &  (optinfo%Cart_Red_Space(i) .GT. NAtoms)) then
                           Write(lupri,'(/,A,I2,A/)') &
                           'The atom in reduced space number should be & 
                           &between 1 and ',(NAtoms),'!' 
                           Call mem_dealloc(optinfo%Cart_Red_Space)
                           Call LSQuit('Illegal coordinate number',lupri)
                        Endif 
                     Enddo
                  Case('.DISPLA')
                      Read(lucmd,*) optinfo%displ
                  Case('.FORBAC')
                      optinfo%ForBac = .TRUE.
                  Case('.SCANSI') 
                     optinfo%simple_scan = .TRUE.
                     Read(lucmd,*) optinfo%ScanCoord
                  Case('.SCANST') 
                     Read(lucmd,*) optinfo%Scan_step
                  Case('.NOHOPE')
                     optinfo%New_stepping = .FALSE.
                  Case('.ITERBT')
                     optinfo%IBT = .TRUE.
                  Case('.OLDIBT')
                     optinfo%OldIBT = .TRUE.
                  Case('.DERORD')
                     Read(lucmd,*) optinfo%Deriv_order
                  Case('.SHANKS')
                     optinfo%Shanks = .TRUE.
                  ! Force-modified PES
                  Case('.FMPES ')
                      optinfo%FMPES = .TRUE.
                      Do i = 1,2
                         Read(lucmd,*) optinfo%Att_atom(i)
                      Enddo
                      Read(lucmd,*) optinfo%Ext_force
                      Write(*,*) optinfo%Att_atom, optinfo%Ext_force

          End select
        Else
           Write(lupri,'(/,3A,/)') ' Keyword "',Keyword, &
                '" not recognized in optimization input'
           Call LSQuit('Illegal keyword in optimization input',lupri)
        End If
Enddo
! The final configuration and printout of input
Call Final_opt_config(optinfo,lupri,NAtoms)
Return
!
End subroutine LS_optimization_input
!=====================!
! Final_opt_config    !
!=====================!
Subroutine Final_opt_config(optinfo,lupri,NAtoms)
! Checking parameters and printing the settings
Implicit none
Type(opt_setting) :: optinfo  ! Information
Integer :: lupri              ! File unit for output 
Integer :: NType,NINFO,i,j,NAtoms,ITmp

NTYPE = 0
!
!     Check if only redundant internal coordinates should be determinded.
!
IF (optinfo%FindRe) THEN
   WRITE(LUPRI,'(A)') ' Determination of redundant internal &
   & coordinates will be performed.'
      WRITE(LUPRI,'(/A/A/)') ' *** NOTE! ***', &
 &         ' No geometry optimization will be done, &
 &         other keywords will be ignored!!!!!'
ELSE
!
!     Check for visualization
!
!      IF (VISUAL) THEN
!         WRITE(LUPRI,'(A)') ' Visualization has been ' //
!     &        'requested. No geometry optimization will be done.',
!     &        ' VRML-file of geometry will be created.'
!         IF (VRBOND) WRITE(LUPRI,'(A)')
!     &        ' Bonds will be drawn between nearby atoms.'
!         IF (VRCORD) WRITE(LUPRI,'(A)')
!     &        ' Coordinate axes will be drawn.'
!         IF (VRSYMM) THEN
!            WRITE(LUPRI,'(A)') ' Symmetry elements will be visualized.',
!     &        ' Please note that symmetry should NOT be specified' //
!     &        ' in the input file.'
!         END IF
!         IF (VREIGV) THEN
!            WRITE(LUPRI,'(A)') ' Eigenvectors can ' //
!     &       'only be visualized during a optimization.',
!     &       ' Keyword will be ignored!'
!            VREIGV = .FALSE.
!         END IF
!         IF (VRVIBA) THEN
!            WRITE(LUPRI,'(A)') ' Vibrational modes can ' //
!     &       'only be visualized during a optimization.',
!     &       ' Keyword will be ignored!'
!            VRVIBA = .FALSE.
!         END IF
!         WRITE(LUPRI,'(A)')
!     &      ' Any other keywords in this module are ignored!'
!         RETURN
!      END IF

!
!     The type of optimization is determined, BFGS is default.
!
      IF (optinfo%Saddle) THEN
         WRITE(LUPRI,'(A)') &
             ' Saddle point optimization has been requested.'
         optinfo%NoAux = .TRUE.
         NType = 0
      END IF
      IF (optinfo%SteepDescend) THEN
         WRITE(LUPRI,'(A)') &
             ' 1st order steepest descent method will be used.'
         NType = NType + 1
      END IF
      IF (optinfo%RANKON) THEN
         IF (optinfo%MULTI) THEN
            WRITE(LUPRI,'(A)') ' 1st order method with &
               &"multiple rank one" updating scheme will be used.'
         ELSE
            WRITE(LUPRI,'(A)') &
                ' 1st order method with rank one update will be used.'
         END IF
         NType = NType + 1
      END IF
      IF (optinfo%BOFILL) THEN
         optinfo%MULTI = .FALSE.
         WRITE(LUPRI,'(A)') &
             ' 1st order method with Bofills update will be used.'
         NType = NType + 1
      END IF
      IF (optinfo%PSB) THEN
         IF (optinfo%MULTI) THEN
            WRITE(LUPRI,'(A)') ' 1st order method with &
               & "multiple PSB" updating scheme will be used.'
         ELSE
            WRITE(LUPRI,'(A)') ' 1st order method with PSB update will be used.'
         END IF
         NType = NType + 1
      END IF
      IF (optinfo%DFP) THEN
         IF (optinfo%MULTI) THEN
            WRITE(LUPRI,'(A)') ' 1st order method with &
            & "multiple DFP" updating scheme will be used.'
         ELSE
            WRITE(LUPRI,'(A)') &
                ' 1st order method with DFP update will be used.'
         END IF
         NType = NType + 1
      END IF
      IF (optinfo%BFGS) THEN
         IF (optinfo%MULTI) THEN
            WRITE(LUPRI,'(A)') ' 1st order method with &
               & "multiple BFGS" updating scheme will be used.'
         ELSE
            WRITE(LUPRI,'(A)') &
                ' 1st order method with BFGS update will be used.'
         END IF
         NType = NType + 1
      END IF
      IF (optinfo%BFGSR1) THEN
         WRITE(LUPRI,'(A)') &
         ' 1st order method with BFGS/rank one combination update will be used.'
         NType = NType + 1
      END IF
!      IF (SCHLEG) THEN
!         WRITE(LUPRI,'(A)')
!     &  ' 1st order method with Schlegels updating scheme will be used.'
!         NType = NType + 1
!      END IF
!      IF (NEWTON) THEN
!         WRITE(LUPRI,'(A)')
!     &        ' 2nd order Newton method will be used.'
!         NType = NType + 1
!      END IF
!      IF (QUADSD) THEN
!         WRITE(LUPRI,'(A)')
!     &        ' 2nd order quadratic steepest ' //
!     &        'descent method will be used.'
!         NType = NType + 1
!      END IF
      IF (optinfo%FirstOrd) THEN
         IF (optinfo%SADDLE) THEN
            WRITE(LUPRI,'(A)') &
               & ' Default 1st order TS-method will be used: Bofills update.'
            optinfo%BOFILL = .TRUE.
         ELSE
            WRITE(LUPRI,'(A)') &
             & ' Default 1st order method will be used: BFGS update.'
            optinfo%BFGS = .TRUE.
         END IF
         optinfo%REJINI = .TRUE.
         NType = NType + 1
      END IF
!      IF (optinfo%SecondOrd) THEN
!         WRITE(LUPRI,'(A)') &
!             ' Default 2nd order method will be used:   Newton method.'
!         optinfo%NEWTON = .TRUE.
!      END IF
      IF (NType .EQ. 0) THEN
         IF (optinfo%SADDLE) THEN
            WRITE(LUPRI,'(A)') &
            ' Default 1st order TS-method will be used: Bofills update.'
            optinfo%BOFILL = .TRUE.
         ELSE
            WRITE(LUPRI,'(A)') &
            ' Default 1st order method will be used: BFGS update.'
            optinfo%BFGS = .TRUE.
         END IF
         optinfo%FirstOrd = .TRUE.
      ELSE IF (NType .GT. 1) THEN
         WRITE(LUPRI,'(/A)') ' WARNING! More than one &
         & optimization method has been selected'
         CALL LSQuit('More than one optimization method chosen' ,lupri)
      END IF
!      IF (HFPROP) THEN
!         WRITE(LUPRI,'(A)') ' The Hellmann-Feynman theorem will ' //
!     &        'be utilized to calculate derivatives.'
!         WRITE (LUPRI,'(5X,A)') 'This option is currently not working'//
!     &        ' correctly, program will stop'
!         CALL QUIT('Hellmann-Feynman approximation not working')
!      END IF
!
      NType = 0
      IF (optinfo%CartCoord) THEN
         WRITE(LUPRI,'(A)') ' Optimization will be performed in Cartesian coordinates.'
         NType = NType + 1
      END IF
      IF (optinfo%REDINT) THEN
         WRITE(LUPRI,'(A)') ' Optimization will be performed in redundant internal coordinates.'
         NType = NType + 1
      END IF
      IF (optinfo%DELINT) THEN
         WRITE(LUPRI,'(A)') ' Optimization will be performed in delocalized internal coordinates.'
         NType = NType + 1
      END IF
      IF (NType .EQ. 0) THEN
         IF (optinfo%FirstOrd) THEN
            optinfo%REDINT = .TRUE.
            WRITE(LUPRI,'(A)') ' Optimization will be performed in redundant internal coordinates.'
         ELSE
            optinfo%CartCoord = .TRUE.
            WRITE(LUPRI,'(A)') ' Optimization will be performed in Cartesian coordinates.'
         END IF
      ELSE IF (NType .GT. 1) THEN
         WRITE(LUPRI,'(/A)') ' WARNING! More than one coordinate system has been selected'
         CALL LSQUIT('More than one coordinate system chosen',lupri)
      END IF
      IF (optinfo%NOAUX) THEN
            WRITE(LUPRI,'(A)') ' No extra (auxiliary) bonds will be added.'
      END IF
      IF (optinfo%NODIHE) THEN
            WRITE(LUPRI,'(A)') ' No dihedral angles will be used &
           & as coordinates (just bonds and angles).'
      END IF
      IF (optinfo%NoAdda) THEN
            WRITE(LUPRI,'(A)') ' No stabilizing additional angles will be used.'
      END IF
      IF (optinfo%RedRed) THEN
            WRITE(LUPRI,'(A)') ' Redundancy will be decreased by&
           & removing a number of dihedral coordinates.'
      END IF
!
      IF (.NOT. (optinfo%HessFile .OR. optinfo%InmdHess .OR. optinfo%ModHes .OR. optinfo%CMBMod &
      .OR. optinfo%InrdHess .OR. (optinfo%EVLINI .GT. -0.9E0_realk))) THEN
         IF (optinfo%SADDLE) THEN
            optinfo%NOAUX  = .TRUE.
         ELSE
            IF (.NOT.optinfo%CartCoord) optinfo%InmdHess = .TRUE.
         END IF
      END IF
!
      IF (optinfo%HessFile) THEN
         IF (optinfo%NEWTON .OR. optinfo%QUADSD) THEN
            call lsquit('NINFO: no value assigned to this variable',-1)
!            NINFO = NINFO + 1
            WRITE(LUPRI,'(A)') &
                'INFO: .HESFIL only has effect when a 1st order&
                & method has been specified => Keyword ignored.'
         ELSE
            WRITE(LUPRI,'(A)') &
                ' Initial Hessian will be read from file.'
         END IF
      END IF
!      IF (INITHS) THEN
!         IF (NEWTON .OR. QUADSD) THEN
!            NINFO = NINFO + 1
!            WRITE(LUPRI,'(A)')
!     &           'INFO: .INITHE only has effect when 1st order ' //
!     &           'method has been specified => Keyword ignored.'
!         ELSE IF (HSFILE) THEN
!            NINFO = NINFO + 1
!            WRITE(LUPRI,'(A)')
!     &           'INFO: .INITHE has no effect when .HESFIL ' //
!     &           'has been specified => Keyword ignored.'
!            INITHS = .FALSE.
!         ELSE
!            WRITE(LUPRI,'(A)')
!     &           ' Initial Hessian will be calculated.'
!         END IF
!      END IF
      IF (optinfo%MODHES) THEN
         WRITE(LUPRI,'(A)') &
  &           ' An approximate model Hessian will be used.',&
  &           ' The model Hessian parameters ', &
  &           'of Roland Lindh will be used.'
           IF (.NOT. optinfo%DELINT) optinfo%REDINT = .TRUE.
      END IF
!
      IF (optinfo%CMBMod) THEN
         IF (optinfo%NEWTON .OR. optinfo%QUADSD) THEN
            call lsquit('NINFO: no value assigned to this variable 838',-1)
!            NINFO = NINFO + 1
            WRITE(LUPRI,'(A)') &
                'INFO: .CMBMOD only has effect when 1st order&
               &  method has been specified => Keyword ignored.'
            optinfo%CMBMOD = .FALSE.
         ELSE IF (optinfo%HessFile) THEN
            call lsquit('NINFO: no value assigned to this variable 838',-1)
!            NINFO = NINFO + 1
            WRITE(LUPRI,'(A)') &
                'INFO: .CMBMOD has no effect when .HESFIL &
              &  has been specified => Keyword ignored.'
           optinfo%CMBMod = .FALSE.
         ELSE
            WRITE(LUPRI,'(A)') &
                ' An approximate model Hessian will be used',&
                &' with the model Hessian parameters of Roland Lindh.',&
                &' The Hessian will be updated &
               & through a combination of a calculated', &
                &' model Hessian and a BFGS update of the last Hessian.'
            IF (.NOT. optinfo%DELINT) optinfo%REDINT = .TRUE.
         END IF
      END IF
      IF (optinfo%InmdHess) THEN
         WRITE(LUPRI,'(A)') &
             ' Model Hessian will be used as initial Hessian.', &
            & ' The model Hessian &
            &   parameters of Roland Lindh will be used.'
      END IF
      IF (optinfo%InrdHess) THEN
         IF (optinfo%NEWTON .OR. optinfo%QUADSD) THEN
            WRITE(LUPRI,'(A)') &
                ' .INIRED only has effect when 1st order & 
             & method has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
            optinfo%InrdHess = .FALSE.
         ELSE IF (optinfo%HessFile) THEN
            WRITE(LUPRI,'(A)') &
               ' .INIRED has no effect when .HESFIL has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
            optinfo%InrdHess = .FALSE.
         ELSE IF (optinfo%InitHess) THEN
            WRITE(LUPRI,'(A)') &
                ' .INIRED has no effect when .INITHE has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
            optinfo%InrdHess = .FALSE.
         ELSE IF (optinfo%InmdHess) THEN
            WRITE(LUPRI,'(A)') &
                ' .INIRED has no effect when .INIMOD has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
            optinfo%InrdHess = .FALSE.
         ELSE
            WRITE(LUPRI,'(A)') ' Initial Hessian will be diagonal in&
           & internal coordinates.'
         END IF
      END IF
      IF (optinfo%EVLINI .GT. 0.0E0_realk) THEN
         IF (optinfo%NEWTON .OR. optinfo%QUADSD) THEN
            WRITE(LUPRI,'(A)') &
            ' .INITEV only has effect when 1st order method has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
         ELSE IF (optinfo%HessFile) THEN
            WRITE(LUPRI,'(A)') ' .INITEV has no effect when .HESFIL has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
         ELSE IF (optinfo%InitHess) THEN
            WRITE(LUPRI,'(A)') &
                ' .INITEV has no effect when initial Hessian is calculated.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
         ELSE IF (optinfo%InmdHess) THEN
            WRITE(LUPRI,'(A)') ' .INITEV has no effect when .INIMOD has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
            optinfo%EVLINI = -1.0E0_realk
         ELSE
            WRITE(LUPRI,'(A,F10.6)') ' Initial diagonal &
           & Hessian will have elements equal to: ', optinfo%EVLINI
         END IF
      ELSE IF (.NOT. (optinfo%DelInt .OR. optinfo%RedInt .OR. optinfo%InrdHess &
             .OR. optinfo%InmdHess)) THEN
         optinfo%EVLINI = 1.0E0_realk
      END IF
      WRITE(LUPRI,*)
      IF (optinfo%NoHessianWrite) THEN
         WRITE(LUPRI,'(A/)') &
            ' Current Hessian will not be written out to DALTON.HES.'
      END IF
      IF (optinfo%REJINI) THEN
         IF (optinfo%NEWTON .OR. optinfo%QUADSD) THEN
            WRITE(LUPRI,'(A)') &
           ' .REJINI only has effect when 1st order method has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
            WRITE(LUPRI,*)
         ELSE
            WRITE(LUPRI,'(A)') ' Hessian will be reinitialized after rejected steps.'
            WRITE(LUPRI,*)
         END IF
      END IF
      IF (optinfo%GradIni) THEN
         IF (optinfo%NEWTON .OR. optinfo%QUADSD) THEN
            WRITE(LUPRI,'(A)') &
            ' .GRDINI only has effect when 1st order method has been specified.'
            WRITE(LUPRI,'(A)') ' Keyword ignored.'
            WRITE(LUPRI,*)
         ELSE
            WRITE(LUPRI,'(A)') &
            ' Hessian will be reinitialized when the norm of the gradient increases.'
            WRITE(LUPRI,*)
         END IF
      END IF
!
      IF (optinfo%IPrint .NE. optinfo%IpDef) THEN
         WRITE(LUPRI,'(A,I10)') ' Print level in OPTIMI  :',optinfo%IPRINT
      END IF
!
      IF (optinfo%DoPre) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,'(A)') ' Preoptimization will be performed with the basis set(s):'
         DO  I = 1, optinfo%NUMPRE
            WRITE(LUPRI,'(A,A60)') '      ',optinfo%PreText(I)
         ENDDO
      END IF
!
      IF (optinfo%DoSpE) THEN
         WRITE(LUPRI,'(A)') &
        ' Single point energy will be calculated using the basis:'
         WRITE(LUPRI,'(A,A60)') '      ',optinfo%SpBText
      END IF
      IF (optinfo%TrustCh) THEN
         WRITE(LUPRI,'(/A/A/,(A,F10.4))') &
             ' Restricted step control parameters', &
             ' ----------------------------------', &
             ' Initial trust radius   :',optinfo%TrustRad, &
             ' Trust radius increment :',optinfo%TrustIn, &
             ' Trust radius decrement :',optinfo%TrustDE, &
             ' Bad prediction ratio   :',optinfo%RTENBD, &
             ' Good prediction ratio  :',optinfo%RTENGD, &
             ' Rejection ratio, low   :',optinfo%RTRJMN, &
             ' Rejection ratio, high  :',optinfo%RTRJMX
      END IF
      WRITE(LUPRI,*)
      IF (optinfo%NoTrust) THEN
         WRITE(LUPRI,*) 'No trust region will be used for steps.'
         WRITE(LUPRI,*)
      END IF
!
      IF (.NOT. optinfo%SADDLE) THEN
         IF (optinfo%GDIIS) THEN
            WRITE(LUPRI,*) &
            'Geometrical DIIS will be used to control step.'
            optinfo%RatFun = .FALSE.
            optinfo%TrustRg = .FALSE.
         ELSE IF (optinfo%RatFun) THEN
            WRITE(LUPRI,*) &
            'Rational function method will be used to control step.'
            WRITE(LUPRI,*)
            optinfo%TrustRg = .FALSE.
         ELSE
            WRITE(LUPRI,*) & 
            'Trust region method will be used to control step (default).'
            WRITE(LUPRI,*)
            optinfo%TrustRg = .TRUE.
         END IF
      ELSE
         IF (optinfo%GDIIS) THEN
            WRITE(LUPRI,*) & 
           'Geometrical DIIS not suitable for saddle point optimization.&
           & Using image function.'
            optinfo%GDIIS = .FALSE.
            optinfo%RatFun = .FALSE.
            optinfo%TrustRg = .TRUE.
         ELSE IF (optinfo%RatFun) THEN
            WRITE(LUPRI,*) & 
           'Partitioned rational function method will be used to control step.'
            optinfo%TrustRg = .FALSE.
         ELSE
            WRITE(LUPRI,*) &
           'Image function method will be used to control step (default).'
            optinfo%TrustRg = .TRUE.
         END IF
         IF (optinfo%NSPMod .LT. 0) THEN
            WRITE(LUPRI,*) & 
            'The eigenvector corresponding to the lowest eigenvalue is chosen'
            WRITE(LUPRI,*) 'as reaction mode (default).'
         ELSE
            WRITE(LUPRI,*) 'Eigenvector #',optinfo%NSPMOD, &
                ' will be used as reaction mode.'
         END IF
      END IF
      IF (optinfo%LnSearch) THEN
         IF (optinfo%SADDLE) THEN
            WRITE(LUPRI,*) &
                'Line search disabled because saddle point is sought.'
            WRITE(LUPRI,*)
            optinfo%LnSearch = .FALSE.
         ELSE
            WRITE(LUPRI,*) &
            'Partial line search with bound quartic polynomial will be employed.'
            WRITE(LUPRI,*)
         END IF
      END IF
! Baker's convergence criteria
IF (optinfo%BAKER) then
    WRITE(LUPRI,'(/A/)')  ' Baker''s convergence criteria will be used.'
ENDIF
!
!
!     Tweaked convergence criterias more suitable for large systems (linsca)
!
IF (optinfo%NatNorm) THEN
 WRITE(LUPRI,'(/A)') &
 ' New thresholding scheme for geometry convergence will be used'
 IF ((optinfo%IcnLevel .GT. 0) .AND. ((optinfo%GradThr .NE. optinfo%DEFTH2) &
    .OR. (optinfo%ThrStep .NE. optinfo%DefTh2))) THEN
    WRITE(LUPRI,'(/,A)') & 
  & ' WARNING! Conflicting specifications of convergence thresholds'
    WRITE(LUPRI,'(A/A/)') & 
   ' Please use *either* ''GRADIE'' and ''STEP T'' *or* one of the keywords',&
  & ' ''.VLOOSE'', ''.LOOSE'', ''.TIGHT'' and ''.VTIGHT.'''
    CALL LSQuit('Conflicting specifications of convergence threshold',lupri)
 ELSE IF (optinfo%IcnLevel .GT. 0) THEN
!
!     VLOOSE, LOOSE, TIGHT and VTIGHT are simply predefined threshold for
!     the gradient convergence criteria
!            
    IF (optinfo%IcnLevel .EQ. 1) THEN
       optinfo%GradThr = 1E-3_realk
       optinfo%ThrStep = 3E-3_realk
    ELSE IF (optinfo%IcnLevel .EQ. 2) THEN
       optinfo%GradThr = 1E-4_realk
       optinfo%ThrStep = 3E-4_realk
    ELSE IF (optinfo%IcnLevel .EQ. 4) THEN
        optinfo%GradThr = 1E-6_realk
        optinfo%ThrStep = 3E-6_realk
    ELSE IF (optinfo%IcnLevel.EQ. 5) THEN
        optinfo%GradThr = 1E-7_realk
        optinfo%ThrStep = 3E-7_realk
    END IF
 ELSE
    IF (optinfo%GradThr .NE. optinfo%DEFTH2) THEN
      IF ((optinfo%GradThr .LE. 0E0_realk) .OR. (optinfo%GradThr .GT. 0.1E0_realk)) &
        &  WRITE(LUPRI,*) 'RMS gradient threshold negative or larger than 0.1, value is reset!'
      ELSE
!     Default value for new scheme
          optinfo%GradThr = 1E-5_realk
      END IF
    ENDIF
    IF (optinfo%ThrStep .NE.optinfo%DefTh2) THEN
      IF ((optinfo%ThrStep .LE. 0E0_realk) .OR. (optinfo%ThrStep .GT. 0.3E0_realk)) &
        &  WRITE(LUPRI,*) 'RMS step threshold negative or larger than 0.3, value is reset!'
      ELSE
!     Default value for new scheme
          optinfo%ThrStep = 3E-5_realk
      END IF
    END IF
    optinfo%ThGradMax = optinfo%GradThr*5E0_realk
    optinfo%ThStepMax = optinfo%ThrStep*5E0_realk
    WRITE(LUPRI,'(A)') ' -------------------------------------------------------------'
    WRITE(LUPRI,'(A,1P,D13.2)') &
   & ' Root-mean-square gradient threshold set to    : ', optinfo%GradThr, &
   & ' Maximum gradient element threshold set to     : ', optinfo%ThGradMax, &
   & ' Root-mean-square step threshold set to        : ', optinfo%ThrStep, &
   & ' Maximum step element threshold set to         : ', optinfo%ThStepMax
! Constrained optimization
 IF (optinfo%ConOpt) THEN
   WRITE(LUPRI,'(/A)') ' Constrained optimization has been requested.'
   optinfo%NoAux = .TRUE.
   IF (.NOT. optinfo%REDINT) THEN
      WRITE(LUPRI,'(A)') ' WARNING! Constrained optimizations&
     & can only be used in conjunction with'
      WRITE(LUPRI,'(A)') ' redundant internal coordinates! Keyword ignored.'
      optinfo%ConOpt = .FALSE.
   ELSE
      WRITE(LUPRI,'(A)') &
      ' The following coordinate numbers will be held fixed during the optimization:'
      DO  I = 1, 8*NAtoms
          IF (optinfo%IConstr(I) .EQ. 1) WRITE(LUPRI,*) '    Coordinate #',I
      ENDDO
      WRITE(LUPRI,*)
   END IF
 END IF
! Remove coordinates
IF (optinfo%RemCoord) THEN
   WRITE(LUPRI,'(A)') ' Removal of coordinates has been requested.'
   IF (.NOT. optinfo%REDINT) THEN
      WRITE(LUPRI,'(A)') ' WARNING! Only internal coordinates can be removed! Keyword ignored.'
      optinfo%RemCoord = .FALSE.
   ELSE
      WRITE(LUPRI,'(A)') ' The following coordinate numbers will be removed: '
      DO I = 1, 8*NAtoms
         IF (optinfo%IConstr(I) .EQ. 2) WRITE(LUPRI,*) '    Coordinate #',I
      ENDDO
      WRITE(LUPRI,*)
   END IF
END IF
! Add coordinates
IF (optinfo%AddCoord) THEN
   WRITE(LUPRI,'(A)') ' Addition of coordinates has been requested.'
   IF (.NOT. optinfo%REDINT) THEN
       WRITE(LUPRI,'(A)') ' WARNING! Only internal coordinates can be added! Keyword ignored.'
       optinfo%AddCoord = .FALSE.
   ELSE
       WRITE(LUPRI,'(A)') ' The following coordinates will be added:'
       DO  I = 1, optinfo%NAdd
           IF (optinfo%AddCoordArray(I,4) .GT. 0) THEN
               WRITE(LUPRI,'(A,4I6)') ' Dihedral between atoms: ', &
              & optinfo%AddCoordArray(I,1), optinfo%AddCoordArray(I,2), &
              & optinfo%AddCoordArray(I,3), optinfo%AddCoordArray(I,4)
           ELSE IF (optinfo%AddCoordArray(I,3) .GT. 0) THEN
               WRITE(LUPRI,'(A,3I6)') ' Angle between atoms   : ', &
              & optinfo%AddCoordArray(I,1), optinfo%AddCoordArray(I,2), optinfo%AddCoordArray(I,3)
           ELSE IF (optinfo%AddCoordArray(I,2) .GT. 0) THEN
               WRITE(LUPRI,'(A,2I6)') ' Bond between atoms    : ', &
              & optinfo%AddCoordArray(I,1), optinfo%AddCoordArray(I,2)
           ELSE
               CALL LSQuit(' Error detected in the specification of additional coordinates!',lupri)
           END IF
           WRITE(LUPRI,*)
       END DO
    ENDIF
END IF
!
!     Experimental atom freezing
!
IF (optinfo%NFreeze .GT. 0) THEN
   WRITE(LUPRI,'(A)') ' Freezing of atoms has been requested.'
   WRITE(LUPRI,'(A)') ' WARNING: .FREEZE is an experimental feature, use it at your own risk...'
   IF (.NOT. optinfo%CONOPT) THEN
       WRITE(LUPRI,'(A)') ' WARNING: Freezing of atoms is really&
      & intended to be used for constrained optimizations!'
   END IF
   WRITE(LUPRI,'(A,I2,A)') ' The following ',optinfo%NFreeze, ' atoms will be held frozen: '
!     Do a simple bubble sort
   DO I = 1, optinfo%NFreeze-1
      DO J = I+1, optinfo%NFreeze
         IF (optinfo%FreezeArray(I) .GT. optinfo%FreezeArray(J)) THEN
            ITMP = optinfo%FreezeArray(I)
            optinfo%FreezeArray(I) = optinfo%FreezeArray(J)
            optinfo%FreezeArray(J) = ITMP
         END IF
      ENDDO
            WRITE(LUPRI,*) '    Atom #      ', optinfo%FreezeArray(I)
   ENDDO
   WRITE(LUPRI,*) '    Atom #      ', optinfo%FreezeArray(optinfo%NFreeze)
   IF (optinfo%IterFreeze .GT. 0) WRITE(LUPRI,'(A,I4,A)') &
  & '     Atoms will be kept frozen for ', optinfo%IterFreeze, ' iterations'
END IF
!
!     VRML options
!
!      IF (VRML) THEN
!         WRITE(LUPRI,'(A)')
!     &        ' VRML-file of geometry will be created.'
!         IF (VRBOND) WRITE(LUPRI,'(A)')
!     &        ' Bonds will be drawn between nearby atoms.'
!         IF (VRCORD) WRITE(LUPRI,'(A)')
!     &        ' Coordinate axes will be drawn.'
!         IF (VRSYMM) THEN
!            WRITE(LUPRI,'(A)') ' Symmetry elements will be visualized.'
!            WRITE(LUPRI,'(A)')
!     &        ' Please note that no symmetry should be specified'
!            WRITE(LUPRI,'(A)')
!     &        ' in the input file.'
!         END IF
!         IF (VREIGV) WRITE(LUPRI,'(A)')
!     &        ' Eigenvectors will be visualized.'
!         IF (VRVIBA) WRITE(LUPRI,'(A)')
!     &        ' Eigenvectors will be visualized.'
!      END IF
!
!     Geometry analysis
!
IF (optinfo%GEOALL) WRITE(LUPRI,'(/A/)') ' Geometry analysis will be printed for each iteration.'
!
!     Adjusted accuracy
!
   IF (optinfo%Stblz .GT. 0) THEN
    WRITE(LUPRI,'(A,I2)') &
   ' Geometries are stabilized by ignoring all &
    &digits beyond digit number ', optinfo%Stblz
    WRITE(LUPRI,'(A/A)') &
             ' WARNING: This is an experimental feature!', & 
            & ' Use it at your own risk the result may or&
             & may not be meaningful...'
   ENDIF
ENDIF
If (optinfo%CartCoord .AND. optinfo%InmdHess) then
   Call lsquit('Model Hessian not implemented for Cartesian optimization.&
   &Try .INIRED keyword instead',lupri)
Endif

IF (optinfo%dynopt) WRITE(LUPRI,'(A)') 'Dynamical optimization is requested.'

!optinfo%ModHes = .FALSE.
!optinfo%InmdHess = .FALSE.
!
End subroutine Final_opt_config
!
End module optimization_input
