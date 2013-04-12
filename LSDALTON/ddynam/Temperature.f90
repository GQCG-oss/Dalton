!================!
!  Temperature   !
!================!
! Thermostats and connected routines
! Written by Vladimir Rybkin
Module Temperature
Use precision
Use ls_util
Use dyn_util
Use fundamental
  Contains 
!===================!
! Gaussian_Random   !
!===================!  
! Generates random numbers from
! Gaussian distribution: mean=0, variance=1
Subroutine Gaussian_Random(GRand)
Implicit none
Real(realk):: GRand
Real(realk) :: x
! Generate uniformly distributed random [0:1]
Call Random_Number(x)
! Generate Gauss-distributed random [-1:1]
Grand = sqrt(-2E0_realk*log(x))*cos(2*Pi*x) 
write(*,*)'Grand',Grand
!
End subroutine Gaussian_random
!====================!
! Maxell_velocities  !
!====================!
Function Maxwell_velocity(Mass,Temperature,Random)
! Generates velocity from Maxwell distribution
! for a given atom
Implicit none
Real(realk):: Mass,Temperature
Real(realk) Maxwell_Velocity(3), Random(3)
Real(realk), parameter :: Boltzmann = 3.166815E0_realk*10E-9 ! [Hartree/(K)]
!
Maxwell_Velocity = sqrt(Boltzmann*Temperature/mass)*Random
Write(*,*)'Random',Random
!
End function Maxwell_velocity
!=======================!
!  Andersen_thermostat  !
!=======================!
! Andersen thermostat: randomly changes
! the velocities to Maxwell distributed ones 
Subroutine Andersen_thermostat(NAtoms,Velocities,Masses,Temperature,&
& Step,Freq,Collision,lupri)
Implicit none
Integer :: NAtoms, Coll_Atom   
Real(realk) :: Velocities(3*NAtoms)
Real(realk) Random(3) ! Vector of Gauss distr. random numbers 
Real(realk) :: Masses(NAtoms)
Real(realk) :: Temperature, Step, Freq
Real(realk) :: R      ! Uniformly distributed random
Integer(kind=4), parameter :: Seed = 191186
Integer :: lupri ! File unit for output
Real(realk) :: zeta, rand
real :: zeta2
Logical :: Collision
!
Collision = .FALSE.
Write(*,*)'TEMPERATURE MODULE'
! Does the collision take place?
!!#ifdef SYS_AIX
!!zeta = rand()
!!#else
!!zeta = rand(seed)
!!#endif

CALL INIT_RANDOM_SEED()
CALL RANDOM_NUMBER(zeta2)
zeta = zeta2

If (Step*Freq .GE. zeta) Collision = .TRUE.
! If collision happens we change velocities of one atom
If (Collision) then
   Write(*,*)'COLLISION'
   Call Gaussian_Random(Random(1))
   Call Gaussian_Random(Random(2))
   Call Gaussian_Random(Random(3))
   Write(*,*)'RANDOM',Random
   Call Random_Number(R)
   Coll_atom = int(R*NAtoms) + 1
   Velocities(Coll_atom*3-2:3*Coll_atom) = &
   & Maxwell_Velocity(Masses(Coll_atom),Temperature,Random)
   Call LSHeader(lupri, 'Collision occured')
   Call SinglelinesInt(lupri,'Maxwell-distributed velocities assigned for atom ',Coll_Atom, 54)
Endif
!
End subroutine Andersen_thermostat
!====================!
! Maxwell_sampling   !
!====================!
! Randomly samples the velocities from
! Maxwellian distribution
Subroutine Maxwell_sampling(NAtoms,Velocities,Masses,Temperature)
Implicit none
Integer :: NAtoms, i
Real(realk) :: Velocities(3*NAtoms), Masses(NAtoms), Random(3)
Real(realk) :: Temperature
!
Do i = 1,NAtoms
   Call Gaussian_Random(Random(1))
   Call Gaussian_Random(Random(2))
   Call Gaussian_Random(Random(3))
   Velocities(i*3-2:3*i) = &
   & Maxwell_Velocity(Masses(i),Temperature,Random)
Enddo   
!
End subroutine Maxwell_sampling
!
End module Temperature
