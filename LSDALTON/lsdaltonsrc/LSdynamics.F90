!=========================================!
! Main driver for linear scaling dynamics !
!=========================================!
! Written by Vladimir Rybkin
Module Dynamics_driver
!
use Precision
use fundamental
use ls_util, only: lsheader
use ls_dynamicsType
use ls_dynamics, only: trajtype, Allocate_traj, Deallocate_traj, Pack_coordinates
use TimeRev_propagation, only: propagation
use Fock_MD, only: FMD_run
use dyn_util, only: DoublelinesInt, Underline, print_vector, calc_kinetic_cart,&
     & calc_kinetic, calc_angmom, calc_angmom_cart, write_phasespace, &
     & mass_weight_vector, project_gradient, SinglelinesInt, final_analysis, &
     & print_temp
use configurationType, only: configitem
use matrix_module, only: matrix
use typedefType, only: lsitem
use ks_settings, only: ks_init_incremental_fock, ks_free_incremental_fock
use files, only: lsopen,lsclose
use lstiming, only: lstimer
use lsdalton_rsp_mod, only: lsdalton_response
use temperature, only: maxwell_sampling, andersen_thermostat, NHC_Hamiltonian, &
     & L_c
use energy_and_deriv, only: get_energy, get_gradient
use memory_handling, only: mem_alloc,mem_dealloc
private
public :: LS_dyn_run
Contains
!===================!
! LS_dyn_run        !
!===================!
Subroutine LS_dyn_run(E,config,H1,F,D,S,CMO,ls)
!Use LSTiming
IMPLICIT NONE
Real(realk) :: E(1)
Type(lsitem), intent(inout) :: ls
Type(Matrix), intent(inout) :: F(1),D(1),S,H1
Type(ConfigItem), intent(inout) :: Config
Type(Matrix), intent(inout) :: CMO       ! Orbitals
Type(trajtype) :: Traj
Integer :: lupri, luerr
Integer :: NAtoms,err,i,j
Logical :: Collision  ! Indicates whether a collision has happened
Logical :: Finished   ! Indicates whether integration should be stopped
Real(realk) :: Coll_Freq
Real(realk) :: CPUTime,WallTime
! Bath collision frequency
Coll_Freq = 0.01/config%dynamics%TimeStep 
!
lupri = config%lupri
luerr = config%luerr
!
! Initializing 
!
Finished = .FALSE.
NAtoms = ls%input%Molecule%NAtoms
! Electronic energy is now potential for the nuclei
traj%CurrPotential = E(1)
! Energy and phase space information will be written to a separate file
! called DALTON.PHS. This file can be used to restart trajectory calculations
config%dynamics%Phase = -1
Call LSOpen(config%dynamics%Phase,'DALTON.PHS','NEW','FORMATTED')
Call Prepare_Integration(NAtoms,ls,config,traj,config%dynamics,D(1)%nrow,lupri) 
!
! Get first gradient
!
Call Calc_gradient(lupri,NAtoms,S,F(1),D(1),ls,config,CMO,traj)
! Getting accelerations
If (.NOT. config%dynamics%Mass_Weight) then
   Do i = 1,NAtoms
      Do j = 1,3
         Traj%Accel(3*(i-1)+j)=-Traj%Gradient(3*(i-1)+j)/Traj%Mass(i)
      Enddo
   Enddo
Endif
! Sampling if needed
If (config%dynamics%MaxSam) then
  Call Double_Maxwell(F,D,S,H1,CMO,ls,config,luerr,lupri,NAtoms,traj,config%dynamics)
Endif
!
Do
! Integrate until termination criteria are fulfilled
  Call Initialize_step(NAtoms,ls,traj,config%dynamics,lupri)
  ! Take step
  Call LSTimer('*START',CPUTime,WallTime,lupri)
  !
  If (config%dynamics%NHChain) then ! Canonical ensemble with Nose-Hoover
     Call NH_chain(F,D,S,H1,CMO,ls,config,luerr,lupri,NAtoms,traj,config%dynamics)
  else ! Microcanonical ensemble
     Call Verlet_step(F,D,S,H1,CMO,ls,config,luerr,lupri,NAtoms,traj,config%dynamics)
  Endif
  ! Get properties
  Call lsdalton_response(ls,config,F(1),D(1),S)
  ! Thermostatting if requested
  If (Config%Dynamics%Andersen) then
     Call Andersen_thermostat(NAtoms,Traj%Velocities,Traj%Mass,config%dynamics%Temp,&
     config%dynamics%TimeStep,Coll_Freq,Collision,lupri)
  Endif
  !
  Call Finalize_Step(lupri,NAtoms,ls,config%dynamics,traj,Finished)
  ! 
  Call LSTimer('Integr. step',CPUTime,WallTime,lupri)
  !
  If (Finished) Exit
Enddo
!
Call Final_printout(ls,config%dynamics,traj,NAtoms,lupri)
! A large negative time signals the end of a trajectory
Write(config%dynamics%Phase,'(F12.5)') -999.99999E0_realk
! Close DALTON.PHS
Call LSClose(config%dynamics%Phase,'KEEP')
!
Call Final_Analysis(NAtoms,traj%StepNum,config%Dynamics%NumTra,&
         &lupri,config%dynamics%Phase,config%dynamics%PrintLevel)
Call Deallocate_traj(traj,config%dynamics%TimRev,config%Dynamics%FockMD)
call mem_dealloc(Config%dynamics%Initial_velocities)
! Dellocate NHC thermostat variables
If (config%dynamics%NHChain) then
   call mem_dealloc(traj%eta)
   call mem_dealloc(traj%v_eta)
   call mem_dealloc(traj%T_array)
   call mem_dealloc(config%dynamics%Q)
Endif
!
End subroutine LS_dyn_run
!======================!
! Prepare_Integration  !
!======================!
! Grabs information and sets initial 
! values of variables, allocates memory
!
Subroutine Prepare_Integration(NAtoms,ls,config,traj,dyn,nbast,lupri)
Implicit none
Integer :: NAtoms,i,j,lupri
Type(lsitem), intent(inout) :: ls
Type(configitem), intent(in) :: config
Type(Dyntype) :: dyn
Type(trajtype),intent(inout) :: Traj
Real(realk) :: CM(3) ! Centre of mass
Real(realk), parameter :: fs2au = 41.3413733365613E0_realk
Integer :: nbast
Real(realk), parameter :: kB = 3.166815E0_realk*10E-6 ! [Hartree/(K)]
! Allocating memory
If (dyn%TimRev) then
   ! Density propagation
   Call Allocate_traj(NAtoms,Traj,nbast,dyn%Filter_order,'TIMREV')
Else 
   If (dyn%FockMD) then
   ! Fock matrix dynamics
      Call Allocate_traj(NAtoms,Traj,nbast,dyn%NPoints,'FOCKMD')
   Else
   ! Regular
      Call Allocate_traj(NAtoms,Traj)
   Endif
Endif
! Allocate NHC thermostat variables and temperature array
If (dyn%NHChain) then
   call mem_alloc(traj%eta,dyn%CLen)
   call mem_alloc(traj%v_eta,dyn%CLen)
   call mem_alloc(traj%T_array,dyn%trajMax+1)
Endif      
   
! Extract atom labels
Do i = 1,NAtoms
   traj%Labels(i) = ls%input%Molecule%Atom(i)%Name
Enddo
! Extract charges of the nuclei 
Do i = 1,NAtoms
   traj%Charges(i) = ls%input%Molecule%Atom(i)%Charge
Enddo
! Printing info to DALTON.PHS
Write(dyn%Phase,'(A11,I7,A8,I7)') 'Trajectory ',1,' out of ',dyn%NumTra
! Initializaing integration time
   traj%TrajTime = 0E0_realk
! Initializing step counter
   traj%StepNum = 0
! Grabbing NAtoms
   NAtoms = config%Molecule%NAtoms
!!!!
!  ToDo, Vladimir: Assigning masses to atoms; 
!         may be not the best place to do it
!
Do i = 1, nAtoms
   ls%input%Molecule%Atom(i)%Mass = &
   XFAMU*Isotopes(ls%input%Molecule%Atom(i)%Atomic_number,1,'MASS',lupri)
Enddo
!
! Grabbing information from ls%input%Molecule to Trajectory
!
Do i = 1,nAtoms
   Traj%Mass(i) = ls%input%Molecule%Atom(i)%Mass
   ! Coordinates
   Do j = 1,3
      Traj%Coordinates(3*(i-1)+j)= &
      ls%input%Molecule%Atom(i)%Center(j) 
   Enddo 
Enddo
! Moving origin to centre of mass
!Call Center_of_Mass(NAtoms,traj%Coordinates,Traj%Mass,CM)
!Call Move_Molecule(NAtoms,traj%Coordinates,-CM)
!Write(*,*)'Centre of mass', CM
! Initial velocities sampled from Maxwell distribution if asked
If (dyn%MaxSam) then
   Call Maxwell_sampling(nAtoms,dyn%Initial_velocities,Traj%Mass,dyn%Temp)
Endif
! Velocities
Traj%Velocities = dyn%Initial_Velocities
! NH chain thermostat velocities, coordinates and masses
If (dyn%NHChain) then
   call mem_alloc(dyn%Q,dyn%CLen)
   dyn%Q(1) = NAtoms*3*kB*dyn%Temp/dyn%omega**2
   Do i = 2, dyn%CLen
      dyn%Q(i) = kB*dyn%Temp/dyn%omega**2
   Enddo
   If (dyn%Init) then ! Initial conditions are read
      Do i = 1, dyn%CLen
         traj%eta(i) = dyn%eta(i)
         traj%v_eta(i) = dyn%v_eta(i)
      Enddo 
      Call mem_dealloc(dyn%eta)
      Call mem_dealloc(dyn%v_eta)
   Else ! They are generated
      Do i = 1, dyn%CLen
         Traj%v_eta(i) = sqrt(kB*dyn%Temp/dyn%Q(i))
         Traj%eta(i) = 0.0E0_realk
      Enddo
   Endif
Endif
! Remove mass-weighting if needed
If (dyn%MWVel) then
   Call Mass_weight_vector(nAtoms,Traj%Velocities,Traj%Mass,'REMOVE')
Endif
! Printing initial conditions
Call LSHeader(lupri, 'Initial geometry (au)')
Call Print_Vector(lupri, NAtoms,traj%Labels, Traj%Coordinates)
If (dyn%PrintLevel >= 1) Then
  Call LSHeader(lupri, 'Initial forces (au)')
  Call Print_Vector(lupri, NAtoms,traj%Labels, -Traj%Gradient)
End If
Call LSHeader(lupri, 'Initial velocities (au)')
Call Print_Vector(lupri, NAtoms,traj%Labels, Traj%Velocities)
! Mass-weight if needed
If (dyn%Mass_Weight) then
   Call Mass_weight_vector(nAtoms,traj%Coordinates,Traj%Mass,'WEIGHT')
   Call Mass_weight_vector(nAtoms,traj%Velocities,Traj%Mass,'WEIGHT')
Endif
!
End subroutine Prepare_integration
!==================!
! Double_Maxwell   !
!==================!
! Experimental classical sampling
Subroutine Double_Maxwell(F,D,S,H1,CMO,ls,config,luerr,lupri,NAtoms,Traj,dyn)
Implicit none
Integer :: lupri,luerr,NAtoms,i,j
Type(lsitem), intent(inout) :: ls
Type(ConfigItem), intent(inout) :: Config
Type(Dyntype), intent(inout) :: dyn
Type(trajtype), intent(inout) :: Traj
Type(Matrix), intent(inout) :: F(1),D(1),S,H1  
Type(Matrix), intent(inout) :: CMO       ! Orbitals
Real(realk) :: Initial_Potential
!
Initial_Potential = traj%CurrPotential
Do i = 1,3
   ! Take the Verlet step 
   Do j = 1,3
      Call Verlet_step(F,D,S,H1,CMO,ls,config,luerr,lupri,NAtoms,traj,config%dynamics)
   Enddo
   ! Do the Maxwell sampling of velocities
   Call Maxwell_sampling(nAtoms,Traj%Velocities,Traj%Mass,dyn%Temp)
Enddo
!
End subroutine Double_Maxwell
!===================!
! Initialize_step   !
!===================!
! Starts the step
Subroutine Initialize_step(NAtoms,ls,traj,dyn,lupri)
Implicit none
Integer :: NAtoms,lupri,i
Type(lsitem), intent(in) :: ls
Type(trajtype) :: traj
Type(dyntype)  :: dyn
Real(realk) :: Temperature
Real(realk), pointer :: Cartesian_Coordinates(:)
Real(realk), pointer :: Cartesian_Velocities(:)
Real(realk), parameter :: Boltzmann = 3.166815E0_realk*10E-6 ! [Hartree/(K)]
! Allocate some memory if integration in mass-weighted
If (dyn%Mass_Weight) then
   Call mem_alloc(Cartesian_Coordinates,3*NAtoms)
   Cartesian_Coordinates = Traj%Coordinates
   Call mem_alloc(Cartesian_Velocities,3*NAtoms)
   Cartesian_Velocities = Traj%Velocities
   ! Remove mass-weighting
   Call Mass_weight_vector(nAtoms,Cartesian_Coordinates,Traj%Mass,'REMOVE')
   Call Mass_weight_vector(nAtoms,Cartesian_Velocities,Traj%Mass,'REMOVE')
Endif
!
  Call SinglelinesInt(lupri, 'Trajectory Step Number ', traj%StepNum, 34)
  Write(lupri,'(A,F10.5)') ' Current time       : ', traj%TrajTime
!
  Call LSHeader(lupri, 'Current geometry (au)')
If (.NOT. dyn%Mass_Weight) then   ! Cartesian 
  Call Print_Vector(lupri, NAtoms, traj%Labels, traj%Coordinates)
Else   ! Mass-weighted
  Call Print_Vector(lupri, NAtoms, traj%Labels, Cartesian_Coordinates)
Endif
! Print NHC thermostat variables
If (dyn%NHChain) then
  Call LSHeader(lupri, 'Current Nose-Hoover chain thermostat coordinates (au):')
  Do i = 1, dyn%CLen
     Write(lupri,'(F15.10)') traj%eta(i)
  Enddo
  Call LSHeader(lupri, 'Current Nose-Hoover chain thermostat velocities (au):')
  Do i = 1, dyn%CLen
     Write(lupri,'(F15.10)') traj%v_eta(i)
  Enddo
Endif
! 
  If (dyn%PrintLevel >= 1) Then
    Call LSHeader(lupri, 'Current forces (au)')
    Call Print_Vector(lupri, NAtoms, traj%Labels, -traj%Gradient)
  End If
  If (dyn%PrintLevel >= 3) Then
     Call LSHeader(lupri, 'Current velocities (au)')
     If (.NOT. dyn%Mass_Weight) then   ! Cartesian 
        Call Print_Vector(lupri, NAtoms, traj%Labels, traj%Velocities)
     Else   ! Mass-weighted
        Call Print_Vector(lupri, NAtoms, traj%Labels, Cartesian_Velocities)
     Endif
  End If
!
! Determine kinetic energy, total energy and angular momentum
!
If (.NOT. dyn%Mass_Weight) then   ! Cartesian 
  Call Calc_Kinetic_Cart(NAtoms*3,NAtoms,traj%Mass,traj%Velocities,traj%CurrKinetic)
Else  ! Mass-weighted
  Call Calc_Kinetic(NAtoms*3,NAtoms,traj%Velocities,traj%CurrKinetic)
Endif
If (dyn%NHchain) then
  Call NHC_Hamiltonian(NAtoms,dyn%CLen,traj%CurrPotential+traj%CurrKinetic,&
       &traj%CurrEnergy,dyn%Q,traj%eta,traj%v_eta,dyn%Temp)
Endif
! Setting initial energy if first step
If (traj%StepNum .EQ. 0) then
   If (.NOT. dyn%NHchain) then
      traj%InitialEnergy = traj%CurrPotential + traj%CurrKinetic
   Else 
      traj%InitialEnergy = traj%CurrEnergy
   Endif
Endif
  !
If (.NOT. dyn%NHchain) then
   traj%CurrEnergy = traj%CurrPotential + traj%CurrKinetic
Endif
!
Call LSHeader(lupri, 'Energy conservation (au)')
Write(lupri,'(3(A,F13.6))') ' Total energy: ', traj%CurrEnergy, &
                            '   Potential energy: ',traj%CurrPotential, &
                            '   Kinetic energy: ', traj%CurrKinetic
Write(lupri,'(31X,A,F14.8)') 'Energy conserv.: ',&
&traj%CurrEnergy-traj%InitialEnergy
If (dyn%NHchain) then
   Write(lupri,'(A)') 'Note: For Nose-Hoover chain thermostat total energy is a &
&conserved quantity and not T+V!'
Endif
! Estimating temperature
Temperature = 2*traj%CurrKinetic/(3*NAtoms*Boltzmann) 
Print *, 'Temperature=',Temperature
Write(lupri,'(31X,A,F14.8)') 'Temperature: ',&
&Temperature
If (dyn%NHChain) traj%T_array(traj%StepNum+1) = Temperature
!
If (.NOT. dyn%Mass_Weight) then   ! Cartesian 
  Call Calc_AngMom_Cart(NAtoms*3,NAtoms,traj%Mass,traj%Coordinates,&
       &traj%Velocities,traj%CurrAngMom)
Else  ! Mass-weighted
  Call Calc_AngMom(NAtoms*3,NAtoms,traj%Coordinates,&
       &traj%Velocities,traj%CurrAngMom)
Endif
  ! Setting initial angular momentum if first step
  If (traj%StepNum .EQ. 0) then
      traj%InitialAngMom = traj%CurrAngMom
  Endif
  !
  Call LSHeader(lupri, 'Angular momentum conservation (au)')
  Write(lupri,'(1X,4(A,F16.10))') 'J_tot:', &
                    Sqrt(Sum(traj%CurrAngMom*traj%CurrAngMom)), &
                '   J_x:', traj%CurrAngMom(1), '   J_y:', traj%CurrAngMom(2), &
                '   J_z:', traj%CurrAngMom(3)
  Write(lupri,'(26X,A,F18.14)') 'Angular moment conserv.: ', &
        Sqrt(Sum((traj%CurrAngMom-traj%InitialAngMom)&
        &*(traj%CurrAngMom-traj%InitialAngMom)))
  Write(lupri,*)
!
! Writing the phase space information to DALTON.PHS
!
If (.NOT. dyn%Mass_Weight) then   ! Cartesian 
    Call Write_PhaseSpace(dyn%Phase,NAtoms,&
          &traj%Charges,traj%Coordinates, &
          traj%Velocities,traj%TrajTime,traj%CurrEnergy,traj%CurrPotential, &
          traj%CurrKinetic,traj%CurrEnergy-traj%InitialEnergy, &
          Sqrt(Sum((traj%CurrAngMom-traj%InitialAngMom)*&
          &(traj%CurrAngMom-traj%InitialAngMom))))
Else  ! Mass-weighted
    Call Write_PhaseSpace(dyn%Phase,NAtoms,&
          &traj%Charges,Cartesian_Coordinates, &
          Cartesian_Velocities,traj%TrajTime,traj%CurrEnergy,traj%CurrPotential, &
          traj%CurrKinetic,traj%CurrEnergy-traj%InitialEnergy, &
          Sqrt(Sum((traj%CurrAngMom-traj%InitialAngMom)*&
          &(traj%CurrAngMom-traj%InitialAngMom))))
Endif
! Deallocate some memory if integration in mass-weighted
If (dyn%Mass_Weight) then
   Call mem_dealloc(Cartesian_Coordinates)
   Call mem_dealloc(Cartesian_Velocities)
Endif
!
End subroutine Initialize_step
!===================!
!  Verlet_step      !
!===================!
! Does the velocity-Verlet step,
! or leap-frog step
Subroutine Verlet_step(F,D,S,H1,CMO,ls,config,luerr,lupri,NAtoms,traj,dyn)
Implicit none
Integer :: lupri,luerr,NAtoms,i,j
Type(lsitem), intent(inout) :: ls
Type(ConfigItem), intent(inout) :: Config
Type(Dyntype), intent(inout) :: dyn
Type(trajtype), intent(inout) :: Traj
Type(Matrix), intent(inout) :: F(1),D(1),S,H1
Type(Matrix), intent(inout) :: CMO       ! Orbitals
Real(realk) :: CPUTime,WallTime
Real(realk), parameter :: fs2au = 41.3413733365613E0_realk
Real(realk), pointer :: MW_Gradient(:)
real(realk) :: Eerr,Etmp(1)
! If mass-weighted coordinates are used or projection requested
If (dyn%Mass_Weight .OR. dyn%Proj_grad) then
   Call mem_alloc(MW_Gradient,3*NAtoms)
   MW_Gradient = Traj%Gradient
   Call Mass_weight_vector(nAtoms,MW_Gradient,Traj%Mass,'REMOVE')
   If (.NOT. dyn%Mass_Weight) Call Mass_weight_vector(nAtoms,traj%Coordinates,Traj%Mass,'WEIGHT')
   ! Project gradient
   Call Project_Gradient(NAtoms,traj%Coordinates,traj%Mass,MW_gradient,dyn%PrintLevel,lupri)
   ! Remove mass-weighting if projection in Cartesian
   If (.NOT. dyn%Mass_Weight) then
      Traj%Gradient = MW_Gradient
      Call Mass_weight_vector(nAtoms,traj%Gradient,Traj%Mass,'WEIGHT')
      Call Mass_weight_vector(nAtoms,traj%Coordinates,Traj%Mass,'REMOVE')
   Endif
Endif
!
! Taking first Verlet half-step
!

!
If (.NOT. dyn%Mass_Weight) then  ! In Cartesian
   Traj%Velocities = &
   Traj%Velocities + Traj%Accel*dyn%TimeStep*0.5E0_realk*fs2au
   Traj%Coordinates = &
   Traj%Coordinates + Traj%Velocities*dyn%TimeStep*fs2au
Else   ! In mass-weighted Cartesian
   Traj%Velocities = &
   Traj%Velocities - MW_Gradient*dyn%TimeStep*0.5E0_realk*fs2au
   Traj%Coordinates = &
   Traj%Coordinates + Traj%Velocities*dyn%TimeStep*fs2au
Endif   
!
! Propagate density matrix if requested
!
If (dyn%TimRev) then
  Call Propagation(dyn%Filter_order,D(1)%nrow,traj%StepNum,&
  & D(1),traj%Darr,traj%Daux,dyn%Start_propagation)
Endif
!
! Extrapolate Fock matrix if requested
!
If (dyn%FockMD) then
   Call FMD_run(traj,F(1),traj%StepNum,dyn%NPoints,dyn%PolyOrd,dyn%Start_propagation)
Endif
!
! New energy
!
If (dyn%Mass_Weight) then ! Remove mass-weighting
    Call Mass_weight_vector(nAtoms,Traj%Coordinates,Traj%Mass,'REMOVE')
Endif
Call LSTimer('*START',CPUTime,WallTime,lupri)
!
Call Pack_coordinates(ls%input%Molecule,Traj%Coordinates,NAtoms)
Call Get_Energy(Etmp,Eerr,config,H1,F,D,S,ls,CMO,NAtoms,lupri,luerr)
traj%CurrPotential = Etmp(1)
Call LSTimer('Energy calc.',CPUTime,WallTime,lupri)
!
! New gradient
!
Call LSTimer('*START',CPUTime,WallTime,lupri)
Call Calc_gradient(lupri,NAtoms,S,F(1),D(1),ls,config,CMO,traj)
Call LSTimer('Forces calc.',CPUTime,WallTime,lupri)
! Mass-weight if needed and project if requested
If (dyn%Mass_Weight .OR. dyn%Proj_grad) then
   MW_Gradient = Traj%Gradient
   Call Mass_weight_vector(nAtoms,MW_Gradient,Traj%Mass,'REMOVE')
   Call Mass_weight_vector(nAtoms,Traj%Coordinates,Traj%Mass,'WEIGHT')
   ! Project gradient
   Call Project_Gradient(NAtoms,traj%Coordinates,traj%Mass,MW_gradient,dyn%PrintLevel,lupri)
   ! Remove mass-weighting if projection in Cartesian
   If (.NOT. dyn%Mass_Weight) then
      Traj%Gradient = MW_Gradient 
      Call Mass_weight_vector(nAtoms,traj%Gradient,Traj%Mass,'WEIGHT')
      Call Mass_weight_vector(nAtoms,traj%Coordinates,Traj%Mass,'REMOVE')
   Endif
Endif
!
! New accelerations
!
If (.NOT. dyn%Mass_Weight) then
   Do i = 1,NAtoms
      Do j = 1,3
       Traj%Accel(3*(i-1)+j)=-Traj%Gradient(3*(i-1)+j)/Traj%Mass(i)
      Enddo
   Enddo
Endif
!
! Second velocity half-step
!
If (.NOT. dyn%Mass_Weight) then   ! In Cartesian
   Traj%Velocities = Traj%Velocities + &
   Traj%Accel*dyn%TimeStep*0.5E0_realk*fs2au
Else  ! In mass-weighted Cartesian
   Traj%Velocities = Traj%Velocities - &
   MW_Gradient*dyn%TimeStep*0.5E0_realk*fs2au
Endif

!
! End of Verlet step
!
!
If (dyn%Mass_Weight)  then
   Call Mass_weight_vector(nAtoms,MW_Gradient,Traj%Mass,'WEIGHT')
   Traj%Gradient = MW_Gradient
Endif
! Deallocate mass-weighted gradient
If (dyn%Mass_Weight .OR. dyn%Proj_grad) Call mem_dealloc(MW_Gradient)
!
End subroutine Verlet_Step
!================!
! NH_chain       !
!================!
! Takes a time step for Nose-Hoover chain thermostat.
! Uses Liouville operators and Trotter expansion.
! According to Frenkel and Smit 'Understanding Molecular Simulation".
Subroutine NH_chain(F,D,S,H1,CMO,ls,config,luerr,lupri,NAtoms,traj,dyn)
Implicit none
Integer :: lupri,luerr,NAtoms,i,j
Type(lsitem), intent(inout) :: ls
Type(ConfigItem), intent(inout) :: Config
Type(Dyntype), intent(inout) :: dyn
Type(trajtype), intent(inout) :: Traj
Type(Matrix), intent(inout) :: F(1),D(1)
Type(Matrix), intent(inout) :: S,H1
Type(Matrix), intent(inout) :: CMO       ! Orbitals
Real(realk) :: CPUTime,WallTime
Real(realk), parameter :: fs2au = 41.3413733365613E0_realk
Real(realk), pointer :: mw_velocities(:)
! Allocate some memory if integration in mass-weighted
If (dyn%Mass_Weight) then
   Call mem_alloc(mw_velocities,3*NAtoms)
   mw_Velocities = Traj%Velocities
   ! Remove mass-weighting
   Call Mass_weight_vector(nAtoms,traj%velocities,Traj%Mass,'REMOVE')
Endif
!
  Call SinglelinesInt(lupri, 'Trajectory Step Number ', traj%StepNum, 34)
  Write(lupri,'(A,F10.5)') ' Current time       : ', traj%TrajTime
! First we apply Lc...
Call L_c(traj%eta,traj%v_eta,dyn%TimeStep*fs2au,dyn%Temp,dyn%CLen,&
& dyn%MStep,traj%velocities,dyn%Q,traj%mass,NAtoms)
! Then take a normal Verlet step
Call Verlet_step(F,D,S,H1,CMO,ls,config,luerr,lupri,NAtoms,traj,dyn)
! Finally apply Lc again
Call L_c(traj%eta,traj%v_eta,dyn%TimeStep*fs2au,dyn%Temp,dyn%CLen,&
& dyn%MStep,traj%velocities,dyn%Q,traj%mass,NAtoms)
! Copy mass-weighted velocities back and deallocate
If (dyn%Mass_Weight) then
   Traj%Velocities = mw_Velocities 
   Call mem_dealloc(mw_velocities)
Endif
End subroutine NH_chain
!===============!
! Finalize_Step !
!===============!
! Does the printout and updating after 
! the step is taken
Subroutine Finalize_Step(lupri,NAtoms,ls,dyn,traj,Finished)
Implicit none
Integer :: lupri,NAtoms,i
Type(lsitem), intent(in) :: ls 
Type(dyntype), intent(inout) :: dyn
Type(trajtype), intent(inout) :: traj
Logical :: Finished
Real(realk), pointer :: Cartesian_Coordinates(:)
Real(realk), pointer :: Cartesian_Velocities(:)
! Allocate some memory if integration in mass-weighted
If (dyn%Mass_Weight) then
   Call mem_alloc(Cartesian_Coordinates,3*NAtoms)
   Cartesian_Coordinates = Traj%Coordinates
   Call mem_alloc(Cartesian_Velocities,3*NAtoms)
   Cartesian_Velocities = Traj%Velocities
   ! Remove mass-weighting
   Call Mass_weight_vector(nAtoms,Cartesian_Coordinates,Traj%Mass,'REMOVE')
   Call Mass_weight_vector(nAtoms,Cartesian_Velocities,Traj%Mass,'REMOVE')
Endif
!
Call LSHeader(lupri, 'New geometry (au)')
If (.NOT. dyn%Mass_Weight) then   ! Cartesian 
  Call Print_Vector(lupri, NAtoms, traj%Labels, traj%Coordinates)
Else   ! Mass-weighted
  Call Print_Vector(lupri, NAtoms, traj%Labels, Cartesian_Coordinates)
Endif
!
If (dyn%PrintLevel >= 3) Then
  Call LSHeader(lupri, 'New forces (au)')
  Call Print_Vector(lupri, NAtoms,traj%Labels, -Traj%Gradient)
End If
!
If (dyn%PrintLevel >= 3) Then
   Call LSHeader(lupri, 'New velocities (au)')
   If (.NOT. dyn%Mass_Weight) then   ! Cartesian 
      Call Print_Vector(lupri, NAtoms, traj%Labels, traj%Velocities)
   Else   ! Mass-weighted
      Call Print_Vector(lupri, NAtoms, traj%Labels, Cartesian_Velocities)
   Endif
Endif
! Print NHC thermostat variables
If (dyn%NHChain) then
  Call LSHeader(lupri, 'New Nose-Hoover chain thermostat coordinates (au):')
  Do i = 1, dyn%CLen
     Write(lupri,'(F15.10)') traj%eta(i)
  Enddo
  Call LSHeader(lupri, 'New Nose-Hoover chain thermostat velocities (au):')
  Do i = 1, dyn%CLen
     Write(lupri,'(F15.10)') traj%v_eta(i)
  Enddo
Endif
!
!
traj%StepNum = traj%StepNum + 1
traj%TrajTime = traj%TrajTime + dyn%TimeStep
!
!  Check whether the termination criteria are fulfilled
!
If ((dyn%trajMax >= 0) .and. (traj%StepNum >= dyn%trajMax)) then
   Finished = .TRUE.
Endif
If ((dyn%MaxTime > 0E0_realk) .and. (traj%TrajTime >= dyn%MaxTime-1E-6_realk)) then
   Finished = .TRUE.
Endif
! Deallocate some memory if integration in mass-weighted
If (dyn%Mass_Weight) then
   Call mem_dealloc(Cartesian_Coordinates)
   Call mem_dealloc(Cartesian_Velocities)
Endif
!
End subroutine Finalize_Step
!=================!
! Final_printout  !
!=================!
Subroutine Final_printout(ls,dyn,traj,NAtoms,lupri)
Implicit none
Integer :: lupri,NAtoms,i
Type(lsitem), intent(in) :: ls 
Type(dyntype), intent(inout) :: dyn
Type(trajtype), intent(inout) :: traj
Real(realk), parameter :: Boltzmann = 3.166815E0_realk*10E-6 ! [Hartree/(K)]
!
Call DoublelinesInt(lupri, 'Final information for trajectory ', 1, 44)
Write(lupri,'(A,F10.5)') ' Final time       : ', traj%TrajTime
!
Call Underline(lupri, 'Final geometry (au)', -1)
Call Print_Vector(lupri, NAtoms,traj%Labels, Traj%Coordinates)
If (dyn%PrintLevel >= 3) Then
  Call Underline(lupri, 'Final forces (au)', -1)
  Call Print_Vector(lupri, NAtoms,traj%Labels, -Traj%Gradient)
End If
Call Underline(lupri, 'Final velocities (au)', -1)
Call Print_Vector(lupri, NAtoms,traj%Labels, Traj%Velocities)
!
! Determine kinetic energy, total energy and angular momentum
!
If (.NOT. dyn%Mass_Weight) then   ! Cartesian 
  Call Calc_Kinetic_Cart(NAtoms*3,NAtoms,traj%Mass,traj%Velocities,traj%CurrKinetic)
Else  ! Mass-weighted
  Call Calc_Kinetic(NAtoms*3,NAtoms,traj%Velocities,traj%CurrKinetic)
Endif
  traj%CurrEnergy = traj%CurrPotential + traj%CurrKinetic
  Call Underline(lupri, 'Final energy conservation (au)', -1)
  Write(lupri,'(3(A,F13.6))') ' Total energy: ', traj%CurrEnergy, &
                              '   Potential energy: ',traj%CurrPotential, &
                              '   Kinetic energy: ', traj%CurrKinetic
  Write(lupri,'(31X,A,F14.8)') 'Energy conserv.: ',&
  &traj%CurrEnergy-traj%InitialEnergy
!
If (.NOT. dyn%Mass_Weight) then   ! Cartesian 
  Call Calc_AngMom_Cart(NAtoms*3,NAtoms,traj%Mass,traj%Coordinates,&
       &traj%Velocities,traj%CurrAngMom)
Else  ! Mass-weighted
  Call Calc_AngMom(NAtoms*3,NAtoms,traj%Coordinates,&
       &traj%Velocities,traj%CurrAngMom)
Endif
  !
  Call Underline(lupri, 'Final angular momentum conservation (au)', -1)
  Write(lupri,'(1X,4(A,F16.10))') 'J_tot:', &
                    Sqrt(Sum(traj%CurrAngMom*traj%CurrAngMom)), &
                '   J_x:', traj%CurrAngMom(1), '   J_y:', traj%CurrAngMom(2), &
                '   J_z:', traj%CurrAngMom(3)
  Write(lupri,'(26X,A,F18.14)') 'Angular moment conserv.: ', &
        Sqrt(Sum((traj%CurrAngMom-traj%InitialAngMom)&
        &*(traj%CurrAngMom-traj%InitialAngMom)))
  Write(lupri,'(//)')
! Save final temperature
If (dyn%NHChain) then
   traj%T_array(traj%StepNum+1) = 2*traj%CurrKinetic/(3*NAtoms*Boltzmann) 
   ! Temperature fluctuation in case of Nose-Hoover chain thermostat
   call print_temp(dyn%trajMax,traj%T_array,lupri)
Endif
!
! Writing the phase space information to DALTON.PHS
!
Call Write_PhaseSpace(dyn%Phase,NAtoms,&
          &traj%Charges,traj%Coordinates, &
          traj%Velocities,traj%TrajTime,traj%CurrEnergy,traj%CurrPotential, &
          traj%CurrKinetic,traj%CurrEnergy-traj%InitialEnergy, &
          Sqrt(Sum((traj%CurrAngMom-traj%InitialAngMom)*&
          &(traj%CurrAngMom-traj%InitialAngMom))))
!
End subroutine Final_printout
!===================!
! Calc_gradient     !
!===================!
Subroutine Calc_gradient(lupri,NAtoms,S,F,D,ls,config,C,traj)
! A brief wrapper to call get_gradient for dynamics
Implicit none
Type(trajtype) :: Traj
Integer :: NAtoms,lupri,i
Type(Matrix), intent(inout),target :: S  ! overlap matrices
Type(Matrix), intent(inout) :: F,D   ! Fock and density matrix
Type(Matrix), intent(inout) :: C       ! Orbitals
Type(lsitem) :: ls
Type(ConfigItem), intent(inout) :: Config ! General information
Real(realk), pointer :: Gradient(:,:)
Real(realk) :: Eerr
Real(realk) :: direction(3),R_a(3),R_b(3)
! Allocate gradient
Call mem_alloc(Gradient,3,NAtoms)
! Calculate gradient
Call Get_Gradient(traj%CurrPotential,Eerr,lupri,NAtoms,S,F,D,ls,config,C,Gradient)
! Expand gradient to traj%Gradient
Do i = 1,NAtoms
   traj%Gradient(3*i-2:3*i) = Gradient(:,i)
Enddo
! Deallocate gradient
Call mem_dealloc(Gradient)
! Add external force for steered molecular dynamics
If (config%dynamics%Steered) then
   ! Define the direction
   R_a = traj%Coordinates(config%dynamics%Att_atom(1)*3-2:config%dynamics%Att_atom(1)*3)
   R_b = traj%Coordinates(config%dynamics%Att_atom(2)*3-2:config%dynamics%Att_atom(2)*3)
   direction = (R_b - R_a)/(sqrt(dot_product(R_b-R_a,R_b-R_a)))
   ! Add external force
   traj%Gradient(config%dynamics%Att_atom(1)*3-2:config%dynamics%Att_atom(1)*3) = &
   traj%Gradient(config%dynamics%Att_atom(1)*3-2:config%dynamics%Att_atom(1)*3) &
   & + direction*config%dynamics%Ext_force
   traj%Gradient(config%dynamics%Att_atom(2)*3-2:config%dynamics%Att_atom(2)*3) = &
   traj%Gradient(config%dynamics%Att_atom(2)*3-2:config%dynamics%Att_atom(2)*3) &
   & - direction*config%dynamics%Ext_force
Endif
!
end subroutine Calc_gradient
!
end module Dynamics_driver


