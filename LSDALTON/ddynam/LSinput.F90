!======================================!
! LINEAR SCALING MOLECULAR DYNAMICS    !
!======================================!
!
Module ls_dynamics
use ls_dynamicsType
Use precision
Use molecule_type 
Use memory_handling
Use ls_util
Use dynamics_param
Use matrix_module
Use matrix_operations
!
!=======================================!
! Contains information for integration  !
!=======================================! 
type trajtype
     ! Coordinates and masses
     Real(realk), pointer :: Coordinates(:)
     Real(realk), pointer :: Mass(:)
     ! Coordinates  and velocities for Nose-Hoover chain thermostat
     Real(realk), pointer :: eta(:)
     Real(realk), pointer :: v_eta(:)
     ! Charges
     Real(realk), pointer :: Charges(:)
     ! Atom symbols
     Character(Len = 4), pointer :: Labels(:)
     ! Gradient
     Real(realk), pointer :: Gradient(:)
     ! Velocities
     Real(realk), pointer :: Velocities(:)
     ! Accelerations
     Real(realk), pointer :: Accel(:)
     ! Current time
     Real(realk) :: CurrTime
     Real(realk) :: InitialEnergy
     ! Current kinetic,potential and total energies
     Real(realk) :: CurrEnergy
     Real(realk) :: CurrKinetic
     Real(realk) :: CurrPotential
     ! Angular momenta 
     Real(realk), dimension(3) :: InitialAngMom
     Real(realk), dimension(3) :: CurrAngMom
     ! Previous potential energy
     Real(realk) :: PrevPotential
     ! Predicted potential
     Real(realk) :: PrdPotential
     ! Step number
     Integer :: StepNum
     ! Trajectory time
     Real(realk) :: TrajTime
     ! Auxiliary density matrices
     Real(realk), pointer :: Daux(:,:,:)
     ! Array for storing optimized density matrices
     Real(realk), pointer :: Darr(:,:,:)
     ! Array for storing Fock matrices
     Type(matrix), pointer :: Fock_array(:) 
     ! Fock matrix dynamics coefficients
     Real(realk), pointer :: fmd_coef(:) 
     ! Temperature array
     Real(realk), pointer :: T_array(:) 
endtype trajtype
Contains
!==============================!
!  Allocate_traj               !
!==============================!
Subroutine Allocate_traj(NAtoms,Trajectory,nbas,N,Job_type)
! Call mem_allocs space for velocities and other entities for
! trajectory integration
Implicit none
Type(trajtype) :: Trajectory
Integer :: NAtoms, i
Integer, intent(in), optional :: nbas,N 
Character(Len=6), optional :: Job_type 
If (NAtoms .GT. 0) then 
   Call mem_alloc(Trajectory%Velocities,NAtoms*3)
   Call mem_alloc(Trajectory%Coordinates,NAtoms*3)
   Call mem_alloc(Trajectory%Gradient,NAtoms*3)
   Call mem_alloc(Trajectory%Accel,NAtoms*3)
   Call mem_alloc(Trajectory%Mass,NAtoms)
   Call mem_alloc(Trajectory%Charges,NAtoms)
   Call mem_alloc(Trajectory%Labels,NAtoms)
Endif
! Allocate memory for propagation of density matrix of Fock matrix
If ((PRESENT(nbas)).AND. (PRESENT(N))) then
   If (Job_type .EQ. 'TIMREV') then  
   ! Density
      Call mem_alloc(Trajectory%Daux,N,nbas,nbas)
      Call mem_alloc(Trajectory%Darr,N-1,nbas,nbas)
   Endif
   ! Fock
   If (Job_type .EQ. 'FOCKMD') then
      Allocate(Trajectory%Fock_array(N))
      ! Initialize matrices
      Do i = 1, N
         Call mat_init(Trajectory%Fock_array(i),nbas,nbas)
      Enddo
      ! Allocate FMD coefficients
      Call mem_alloc(Trajectory%fmd_coef,N)
   Endif
Endif
!
End subroutine  Allocate_traj
!===========================!
! Deallocate_traj           !
!===========================!
Subroutine Deallocate_traj(Trajectory,TimRev,FMD)
! DeCall mem_allocs space for velocities and other entities for
! trajectory integration
Implicit none
Logical :: TimRev, FMD ! Niklasson, Fock matrix dynamics
Type(trajtype) :: Trajectory
integer :: i
   Call mem_dealloc (Trajectory%Velocities)
   Call mem_dealloc (Trajectory%Coordinates)
   Call mem_dealloc (Trajectory%Gradient)
   Call mem_dealloc (Trajectory%Accel)
   Call mem_dealloc (Trajectory%Mass)
   Call mem_dealloc (Trajectory%Charges)
   Call mem_dealloc (Trajectory%Labels)
   If (TimRev) then
      Call mem_dealloc(Trajectory%Darr)            
      Call mem_dealloc(Trajectory%Daux)
   Endif            
!
   If (FMD) then
      Do i = 1, size(Trajectory%Fock_array)
         Call mat_free(Trajectory%Fock_array(i))
      Enddo
      Deallocate(Trajectory%Fock_array)
      Call mem_dealloc(Trajectory%fmd_coef)
   Endif
!
End subroutine  Deallocate_traj
!=========================!
!   LS_dynamics_init      !
!=========================!
Subroutine LS_dynamics_init(dynamics_input)
  Implicit None
  Type(dyntype)           :: Dynamics_Input
!
! Initialization
!
  Dynamics_Input%do_dynamics = .FALSE.
  Dynamics_Input%PrintLevel = -1
  Dynamics_Input%MaxSam = .False.
  Dynamics_Input%Andersen = .False.
  Dynamics_Input%FockMD   = .False.
  Dynamics_Input%TimRev   = .False. 
  Dynamics_Input%Start_propagation = .False.
  Dynamics_Input%Orb_Con   = .False. 
  Dynamics_Input%MWVel = .False.
  Dynamics_Input%Mass_Weight = .False.
  Dynamics_Input%Filter_order = -1
  Dynamics_Input%PolyOrd  = -1
  Dynamics_Input%NPoints  = -1
  Dynamics_Input%RanSeed   = 398465
  Dynamics_Input% GetRandom = -1
  Dynamics_Input%LowerRandom = -1.1E3_realk
  Dynamics_Input%UpperRandom = -1.1E3_realk
  Dynamics_Input%trajMax = 100
  Dynamics_Input%MaxTime = -1.1E3_realk
  Dynamics_Input%Temp    = -1.1E3_realk
  Dynamics_Input%RotTemp = -1.1E3_realk
  Dynamics_Input%VibTemp = -1.1E3_realk
  Dynamics_Input%PathLen  = -1.1E3_realk
  Dynamics_Input%StepLen  = -1.1E3_realk
  Dynamics_Input%TimeStep = 0.5E0_realk
  Dynamics_Input%PathL    = .False.
  Dynamics_Input%StepL    = .False.
  Dynamics_Input%TimeL    = .False.
  Dynamics_Input%FragDist = -1.1E3_realk
  Dynamics_Input%CMDist   = -1.1E3_realk
  Dynamics_Input%FragSize = -1.1E3_realk
  Dynamics_Input%IntStepSize = 1E-3_realk
  Dynamics_Input%NumIntAccuracy = 1E-12_realk
  Dynamics_Input%Verlet = .TRUE.
  Dynamics_Input%Steered = .FALSE.
  Dynamics_Input%update_step = .FALSE.
  Dynamics_Input%Proj_grad = .FALSE.
  ! Nose-Hoover chain thermostat is turned off by default
  Dynamics_Input%NHChain = .FALSE.
  ! Thermostat 'masses' are set to 0.011 a.u. (2500 cm-1) by default
  Dynamics_Input%omega = 0.0110E0_realk
  ! Nose-Hoover chain length is set to 3 by default
  Dynamics_Input%CLen = 3
  ! Number of Nose-Hoover multisteps is set to 1 by default
  Dynamics_Input%MStep = 1
End subroutine LS_dynamics_init
!=========================!
!   LS_dynamics_input     !
!=========================!
Subroutine LS_dynamics_input(dynamics_input,ReadWord,keyword,lucmd,lupri,NAtoms)
!
! Input processing routine for linear scaling  dynamics parameters
!
  Use LSTiming
  Implicit None
  Type(dyntype)           :: Dynamics_Input
  Integer                 :: lucmd,lupri
  Integer                 :: NAtoms
  Integer                 :: NError, NWarning, NumVel, I, J, NAtm
  Integer, Parameter      :: DefUpdMethod = 3, DefNUpdates = 5, DefInteg = 5
  Integer, Dimension(1:2) :: KWordInd
  Character(Len = 1)      :: Prompt
  Character(Len = 70)     :: Keyword 
  Character(Len = 7), Dimension(1:53) :: KWordTable = &
    (/ '.2NDORD', '.5THORD','.NEW5OR' , '.ACCURA', '.CMDIST', '.DISTCR', &
       '.FRAGDI', '.FRAGSI', '.HESUPD', '.INTSTP', '.INTEGR', '.ANDERS', &
       '.MAX IT', '.MAX TI', '.NOMOVE', '.NOPROJ', '.NUMTRA','.MAXSAM' , &
       '.NUPDAT', '.OLDALG', '.OPTION', '.ORIENT', '.PATH L', '.PATHWA', &
       '.PRINT ', '.QUADRA', '.RANDOM', '.RANMAX', '.RANMIN', '.ROTSAM', &
       '.ROTTEM', '.SEED  ', '.STEP L', '.TEMP  ', '.TIMEST', '.TRAJEC', &
       '.VELOCI', '.VERLET', '.VIBSAM', '.VIBTEM','.FOCKMD' ,'.TIMREV' , &
       '.MASSWE','.MWVEL ' ,'.ORBCON' , '.STEPUP','.STMDYN' ,'.NHCHAI' , &
       '.BATHFR','.CHAINL', '.MULTIS' , '.ININHC' ,'.PROJGR' /)
  Logical :: file_exist,Const_TS
  Logical,intent(inout)   :: ReadWord
  Integer :: FileStatus
!
Call LSHeader(lupri, &
'Input processing for the direct dynamics module')
NError          = 0
NWarning        = 0
!
! Allocating initial velocities
!
Call mem_alloc(Dynamics_input%Initial_velocities,NAtoms*3)

!
! Read keywords
!
Write(*,*)'LS DYNAMICS INPUT'
Do
Write(*,*)'READWORD',ReadWord
      Read(lucmd,'(A70)', IOStat = FileStatus) Keyword
            WRITE(*,*)'READ KEYWORDS',Keyword(1:7)
      If (FileStatus > 0) Call LSQuit('Error reading lucmd',lupri)
      If ((FileStatus < 0)) Exit
      Prompt = Keyword(1:1)
      If ((Prompt == '#') .or. (Prompt == '!')) Cycle
      If (Prompt == '*') then
         ReadWord = .FALSE.
         Exit
      Endif
      If (Prompt == '.') Then
        If (Any(KWordTable(1:size(KWordTable)) .EQ. Keyword(1:7))) Then
          Select Case (Keyword(1:7))
            Case('.OPTION')
            ! Fock matrix dynamics input
            Case('.FOCKMD')
              Dynamics_Input%FockMD = .True.
              Read(lucmd,*) Dynamics_Input%NPoints
              Read(lucmd,*) Dynamics_Input%PolyOrd
            ! Time-reversible propagation
            Case('.TIMREV')
              Dynamics_Input%TimRev = .True.
              Read(lucmd,*) Dynamics_Input%Filter_order
            ! Orbital connection
            Case('.ORBCON')
              Dynamics_Input%Orb_Con = .True.
            !
            Case('.ACCURA')
              Read(lucmd,*) Dynamics_Input%NumIntAccuracy
            Case('.CMDIST')
              Read(lucmd,*) Dynamics_Input%CMDist
            Case('.FRAGDI')
              Read(lucmd,*) Dynamics_Input%FragDist
            Case('.FRAGSI')
              Read(lucmd,*) Dynamics_Input%FragSize
            Case('.INTSTP')
              Read(lucmd,*) Dynamics_Input%IntStepSize
            Case('.NUMTRA')
              Read(lucmd,*) Dynamics_Input%NumTra
              Write(*,*)'NUMTRA',Dynamics_Input%Numtra
            Case('.MAX IT')
              Read(lucmd,*) Dynamics_Input%trajMax
            Case('.MAX TI')
              Read(lucmd,*) Dynamics_Input%MaxTime
            Case('.PATH L')
              Read(lucmd,*) Dynamics_Input%PathLen
            Case('.RANDOM')
              Read(lucmd,*) Dynamics_Input%GetRandom
            Case('.RANMAX')
              Read(lucmd,*) Dynamics_Input%UpperRandom
            Case('.RANMIN')
              Read(lucmd,*) Dynamics_Input%LowerRandom
            Case('.STEP L')
              Read(lucmd,*) Dynamics_Input%StepLen
            Case('.MAXSAM')
              Dynamics_Input%MaxSam = .TRUE.
            Case('.TEMP  ')
              Read(lucmd,*) Dynamics_Input%Temp
            Case('.ANDERS')
              Dynamics_Input%Andersen = .TRUE.
            Case('.NHCHAI')
              Dynamics_Input%NHChain = .TRUE.
            Case('.CHAINL')
              Read(lucmd,*) Dynamics_Input%CLen
            ! If initial conditions for thermostat are
            ! specified (say, to restart)
            Case('.ININHC')
               Dynamics_Input%Init = .TRUE.
               !
               ! Allocate initial conditions for Nose-Hoover chain
               !
               Call mem_alloc(Dynamics_Input%eta,Dynamics_Input%CLen)
               Call mem_alloc(Dynamics_Input%v_eta,Dynamics_Input%CLen)
               Do i = 1, Dynamics_Input%Clen
                  Read(lucmd,*) Dynamics_Input%eta(i)
               Enddo
               Do i = 1, Dynamics_Input%Clen
                  Read(lucmd,*) Dynamics_Input%v_eta(i)
               Enddo
            Case('.MULTIS')
              Read(lucmd,*) Dynamics_Input%MStep
            Case('.BATHFR')
              Read(lucmd,*) Dynamics_Input%omega
            Case('.PROJGR')
              Dynamics_Input%Proj_grad = .TRUE.
            Case('.TIMEST')
              Read(lucmd,*) Dynamics_Input%TimeStep
              Write(*,*)'TIMESTEP=',Dynamics_Input%TimeStep 
            Case('.MWVEL ')
              Dynamics_Input%MWVel = .TRUE. 
            Case('.MASSWE')
              Dynamics_Input%Mass_Weight = .TRUE.
            Case('.STEPUP')
              Dynamics_Input%update_step = .TRUE.
              Read(lucmd,*) Dynamics_Input%NPoints
              Read(lucmd,*) Dynamics_Input%Poly_ord
              Call mem_alloc(Dynamics_Input%Energy_array,Dynamics_Input%NPoints+1)
              Dynamics_Input%Energy_array = 0.0D0
              Call mem_alloc(Dynamics_Input%Time_array,Dynamics_Input%NPoints+1)
              Dynamics_Input%Time_array = 0.0D0
            Case('.VELOCI')
            WRITE(*,*)'READ VELOCITIES'
              Dynamics_Input%InputVeloc = .True.
              Read(lucmd,*) NumVel
              Write(*,*)'NumVel',NumVel
              If (NumVel < 0)  Then
                Write(lupri,'(/,A,I4,A/)') &
                ' Number of velocities to be read in must be at least 1'
                Call LSQuit('Illegal number of velocities.',lupri)
              Else
                 If (NumVel .EQ. 0) Then
                    Dynamics_input%Initial_velocities = 0E0_realk
                    Else
                        If (NumVel .NE. NAtoms) Then
                           Write(lupri,'(/,A,I4,A)') &
                          ' Number of velocities to be read is not equal to number  &
                          & of atoms'
                          Call LSQuit('Illegal number of velocities.',lupri)
                        Else
                            Do I = 1, NumVel
                               Read(lucmd,*) &
                               (Dynamics_input%Initial_velocities(J), &
                               J = (I-1)*3+1, (I-1)*3+3)
                            End Do
                        Endif
                 Endif
              Endif
            Case('.VERLET')
                Dynamics_input%Verlet = .TRUE.
            ! Steered dynamics
            Case('.STMDYN')
                Dynamics_input%Steered = .TRUE.
                Do i = 1,2
                   Read(lucmd,*) Dynamics_input%Att_atom(i)
                Enddo
                Read(lucmd,*) Dynamics_input%Ext_force
            !
            Case('.PRINT ')
              Read(lucmd,*) Dynamics_input%PrintLevel
!            Case('.VIBSAM')
!              Read(lucmd,*) VibSampTxt
!            Case('.VIBTEM')
!              Read(lucmd,*) Vibtemp
!            Case Default
              Write(lupri,'(/,3A,/)') ' Keyword "',Keyword, &
                '" is not yet implemented in Dynamics_Input_Proc.'
          End Select
        Else
          Write(lupri,'(/,3A,/)') ' Keyword "',Keyword, &
            '" not recognized in Dynamics_Input_Proc.'
          Call LSQuit('Illegal keyword in Dynamics_Input_Proc.',lupri)
        End If
      Endif
   Enddo
!
! Printing introduction and starting timer
!
  Call Print_Info(lupri)
!
! The chosen parameters are printed
!
    If (Dynamics_Input%GetRandom > 0) Then
      If (Dynamics_Input%LowerRandom < 0E0_realk) then
         LowerRandom = 0E0_realk
      Endif
      If (Dynamics_Input%UpperRandom < 0E0_realk) then
         Dynamics_Input%UpperRandom = 1E0_realk
      Endif
      Write(lupri,'(A)') ' Generation of random numbers has been requested.'
      Write(lupri,'(A)') ' All other keywords except .SEED, .RANMIN and .RANMAX will ' &
                           // 'thus be ignored!'
      Write(lupri,'(A,I15)') ' Random numbers to be generated: ',&
                               Dynamics_Input%GetRandom
      Write(lupri,'(A,F10.6,A,F10.6,A)') ' Random num. to be sampled in:   [', &
        Dynamics_Input%LowerRandom, ', ', Dynamics_Input%UpperRandom, ']'
      Write(lupri,'(A,I15)') ' Chosen random generator seed  : ', RanSeed
      ProcessedInput = .True.
      Return
    Else 
      If (Dynamics_Input%LowerRandom > 0E0_realk) &
          Write(lupri,'(A)') ' .RANMIN only has effect when sampling random numbers. Keyword ignored.'
      If (Dynamics_Input%UpperRandom > 0E0_realk) &
          Write(lupri,'(A)') ' .RANMAX only has effect when sampling random numbers. Keyword ignored.'
      Dynamics_Input%LowerRandom = -1.1E3_realk
      Dynamics_Input%UpperRandom = -1.1E3_realk
    End If
!
    If (Dynamics_Input%NHchain) then 
      Write(lupri,'(A)') ' Nose-Hoover chain thermostat is used'
      Write(lupri,'(A,F10.6,A)') ' The temperature is ', &
      & Dynamics_Input%Temp, ' K'
    Endif
!   For the first time we run only Verlet!   
    Write(lupri,'(A)') ' Verlet algorithm (1st order) will be used.'
    Write(lupri,'(A,F10.5)') ' Num. integ. step size  : ', Dynamics_Input%TimeStep

    If (Dynamics_Input%trajMax /= -999) Write(lupri,'(A,I10)') &
                         ' Max # of iterations    : ', Dynamics_Input%trajMax
    If (Dynamics_Input%MaxTime > -1E3_realk) Write(lupri,'(A,F10.5)') &
                         ' Max allowed time (fs)  : ', Dynamics_Input%MaxTime
! Some more defaults
    If ((Dynamics_Input%trajMax == -999) .and. (Dynamics_Input%MaxTime < -1E3_realk)) Then
      Dynamics_Input%trajMax = 100
      Dynamics_Input%MaxTime = 100E0_realk
    End If
    If (Dynamics_Input%InputVeloc) &
      Write(lupri,'(A)') ' Velocities have been explicitly specified.'
!
! Sampling
!
    Write(lupri,'(A,I12)') ' Random seed:             ',Dynamics_Input%RanSeed
    If (SampleOrient) &
      Write(lupri,'(A)') ' Orientation will be sampled (random rotation).'
!
  ProcessedInput = .True.
End Subroutine LS_dynamics_input
!===================!
! Pack_coordinates  !
!===================!
Subroutine Pack_Coordinates(Molecule,Coordinates,NAtoms)
! Updates coordinates in config%molecule
Implicit none
Type(Moleculeinfo) :: Molecule
Integer :: NAtoms, i, j
Real(realk) :: Coordinates(NAtoms*3)
!
Do i = 1,NAtoms
   Do j=1,3
      Molecule%Atom(i)%Center(j)=Coordinates(3*(i-1)+j)
   Enddo
Enddo
!
End subroutine Pack_Coordinates
!================
!   Print_Info
!================
Subroutine Print_Info(lupri)
!
! Prints information concerning the dynamics calculation

  Implicit None
  Integer            :: lupri
  Write(lupri,'(//A)') ' ====================================================='
  Write(lupri,'(A)')   '    -----------------------------------------------'
  Write(lupri,'(A)')   '      Calculation of classical trajectories using'
  Write(lupri,'(A)')   '      direct ab initio molecular dynamics (AIMD)'
  Write(lupri,'(A)')   '    -----------------------------------------------'
  Write(lupri,'(A/)')  ' ====================================================='
End Subroutine Print_Info

End module ls_dynamics
