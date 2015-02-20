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
!=========!
! L_c     !
!=========!
! Liouville operator for Nose-Hoover chain thermostat
! variables for the half-step
! G. Martyna et al. Mol. Phys. 87, 1117 (1996).
Subroutine L_c(eta,v_eta,dt,T,CLen,MStep,vel,Q,mass,NAtoms)
Implicit none
Integer :: CLen,NAtoms,MStep ! Chain length,number of atoms,multistep
Integer :: step, ys, L, i  
Real(realk), dimension(NAtoms) :: mass
Real(realk), dimension(NAtoms*3) :: vel
Real(realk), dimension(CLen) :: Q(CLen), eta(CLen), v_eta(CLen), G(CLen)
Real(realk) :: T, dt, scal, ts, AA
Integer, parameter :: Nys = 3
Real(realk), dimension(3), parameter :: &
& w = (/1.35120719196E0_realk,-1.70241438392E0_realk,1.35120719196E0_realk/)
Real(realk), parameter :: kB = 3.166815E0_realk*10E-6 ! [Hartree/(K)]
! Update forces/accelerations
! First thermostat accelerations
scal = 1.0E0_realk
G(1) = G_1(NAtoms,mass,Q(1),vel,T,scal)
Do L = 1, CLen - 1
   G(L+1) = G_m(v_eta(L),Q(L),Q(L+1),T)
Enddo
! Multistep procedure
Do step = 1, MStep
   Do ys = 1, Nys
      ts = w(ys)*dt/MStep
      v_eta(CLen) = v_eta(CLen) + G(CLen)*ts/4.0E0_realk
      Do L = 1, CLen-1
      ! Thermostat velocities
         AA = exp(-(ts/8.0E0_realk)*v_eta(CLen-L+1))
         v_eta(CLen-L) = v_eta(CLen-L)*AA**2 + &
              &(ts/4.0E0_realk)*G(CLen-L)*AA
      Enddo ! CLen 
      ! Atomic velocities
      AA = exp(-(ts/2.0E0_realk)*v_eta(1))
      scal = scal*AA
      ! G1
      G(1) = G_1(NAtoms,mass,Q(1),vel,T,scal)
      ! Thermostat positions
      Do L = 1, CLen 
         eta(L) = eta(L) + (ts/2.0E0_realk)*v_eta(L)
      Enddo
      ! Thermostat velocities
      Do L = 1, CLen - 1
         AA = exp(-(ts/8.0E0_realk)*v_eta(L+1))
         v_eta(L) = v_eta(L)*AA**2 + (ts/4.0E0_realk)*G(L)*AA
         G(L+1) = G_m(v_eta(L),Q(L),Q(L+1),T)
      Enddo
      v_eta(CLen) = v_eta(CLen) + G(CLen)*(ts/4.0E0_realk)  
   Enddo ! Nys
Enddo ! Mstep
! Atomic velocities
Do i = 1, NAtoms*3
   vel(i) = vel(i)*scal
Enddo
!
End subroutine L_c
!======!
! G_1  !
!======!
! Calculates G1 factor 
! for NH chain thermostat
Function G_1(NAtoms,mass,Q1,vel,T,scal)
Implicit none
Integer :: NAtoms, i 
Real(realk) mass(NAtoms), vel(NAtoms*3)
Real(realk) :: G_1, Q1, T, scal
Real(realk), parameter :: kB = 3.166815E0_realk*10E-6 ! [Hartree/(K)]
!
G_1 = 0.0E0_realk
Do i = 1, NAtoms
   G_1 = G_1 + mass(i)*dot_product(vel(i*3-2:i*3),vel(i*3-2:i*3))
Enddo
!
G_1 = (G_1*scal**2 - (NAtoms*3)*kB*T)/Q1
! 
End function G_1
!=========!
! G_m     !
!=========!
! Calculates NHC thermostat forces
Function G_m(v,Q,Qm,T)
Implicit none 
Real(realk) :: V,Q,Qm,T,G_m ! V,Q - Vm-1 and Qm-1,Qm  - Qm.
Real(realk), parameter :: kB = 3.166815E0_realk*10E-6 ! [Boltzmann, Hartree/(K)]
!
G_m = (Q*V**2 - kB*T)/Qm
!
End function G_m
!=======================!
! Calc NHC_Hamiltonian  !
!=======================!
! Calculates the conserved
! Nose-Hoover chain "energy"
Subroutine NHC_Hamiltonian(NAtoms,CLen,E,E_NHC,Q,eta,v_eta,T)
Implicit none 
Integer :: NAtoms, CLen,i
Real(realk) :: E, E_NHC, T
Real(realk), dimension(CLen) :: Q(CLen), eta(CLen), v_eta(CLen)
Real(realk), parameter :: kB = 3.166815E0_realk*10E-6 ! [Boltzmann, Hartree/(K)]
!
E_NHC = E
Do i = 1, CLen
   E_NHC = E_NHC + 0.50E0_realk*Q(i)*v_eta(i)**2
   If (i .EQ. 1) then
      E_NHC = E_NHC + 3*NAtoms*kB*T*eta(i)
   Else
      E_NHC = E_NHC + kB*T*eta(i)
   Endif
Enddo 
!
End subroutine NHC_Hamiltonian
!
End module Temperature
