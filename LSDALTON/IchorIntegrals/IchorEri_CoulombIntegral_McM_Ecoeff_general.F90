!> @file
!> Contains the general McM driver 

!> \brief General McMurchie-Davidson Integral scheme
!> \author T. Kjaergaard
!> \date 2013 
MODULE IchorEriCoulombintegralCPUMcMGeneralEcoeffMod
use IchorPrecisionModule
use IchorCommonModule
use IchorEriCoulombintegralCPUMcMGeneralWTUVMod
use IchorEriCoulombintegralCPUMcMspecREcoeffMod
use IchorEriCoulombintegralCPUMcMspecLEcoeffMod
use IchorEriCoulombintegralCPUMcMspecL2EcoeffMod
private 
public :: Ichorbuild_Ecoeff_LHS,Ichorbuild_Ecoeff_RHS, printEcoeff

CONTAINS
Subroutine Ichorbuild_Ecoeff_RHS(nPrimP,nPrimA,nPrimB,maxAngP,maxAng1,maxAng2,nTUV,&
     & ijk,Aexp,Bexp,EcoeffN,Pdistance12,PreExpFac,intprint,lupri,ETIJ,PA,PB,&
     & HPINV,TWOA,TWOB)
  implicit none
  integer,intent(in)        :: nPrimP,nPrimA,nPrimB,lupri,nTUV,ijk
  integer,intent(in)        :: maxAngP,maxAng1,maxAng2,intprint
  real(realk),intent(in)    :: Aexp(nPrimA),Bexp(nPrimB),Pdistance12(3)
  real(realk),intent(in)    :: PreExpFac(nPrimP)
  real(realk),intent(inout) :: Ecoeffn(nPrimP,nTUV,ijk)
  real(realk),intent(inout) :: ETIJ(nPrimP,0:maxAngP,0:maxAng1,0:maxAng2,3)
  real(realk),intent(inout) :: PA(nPrimP,3),PB(nPrimP,3)
  real(realk),intent(inout) :: HPINV(nPrimP),TWOA(nPrimP),TWOB(nPrimP)
  !
  real(realk),parameter :: DHALF=0.5E0_realk,D2=2.0E0_realk
  integer :: iPrimP,I,J,X
  integer,parameter :: nPasses = 1
  real(realk) :: bexpo
  CALL LS_OMP_DZERO(Ecoeffn,nPrimP*nTUV*ijk)
  call build_auxiliary_RHS(nPrimA,nPrimB,Aexp,Bexp,HPINV,TWOA,TWOB,PA,PB,Pdistance12)
  IF(maxAng1+maxAng2.GT.6)THEN
     CALL ICHOR_HERM_ECOEFFS_GENERAL_RHS(ETIJ,nPrimP,maxAng1,maxAng2,&
          & PA,PB,HPINV,TWOA,TWOB,LUPRI)  
  ELSE !special less than 6
     CALL ICHOR_HERM_ECOEFFS(ETIJ,nPrimP,nPasses,maxAng1,maxAng2,&
          & PA,PB,HPINV,TWOA,TWOB,LUPRI)  
  ENDIF
  IF(intprint.GT.30)THEN
!$OMP MASTER
     DO X=1,3
        call PRINT_ETIJ(ETIJ,nPrimP,nPasses,maxAng1,maxAng2,X,LUPRI)
     ENDDO
!$OMP END MASTER
!$OMP BARRIER
  ENDIF
  IF(maxAngP.EQ.0)THEN
     call ICHOR_Ecoeffn_maxAngP0_maxAngA0_RHS(nPrimP,Ecoeffn,ETIJ,PreExpFac)
  ELSEIF(maxAngP.EQ.1)THEN
     IF(maxAng1.EQ.1)THEN
        call ICHOR_Ecoeffn_maxAngP1_maxAngA1_RHS(nPrimP,Ecoeffn,ETIJ,PreExpFac)
     ELSEIF(maxAng1.EQ.0)THEN
        call ICHOR_Ecoeffn_maxAngP1_maxAngA0_RHS(nPrimP,Ecoeffn,ETIJ,PreExpFac)
     ENDIF
  ELSEIF(maxAngP.EQ.2)THEN 
     IF(maxAng1.EQ.2)THEN      !(DS)
        call ICHOR_Ecoeffn_maxAngP2_maxAngA2_RHS(nPrimP,Ecoeffn,ETIJ,PreExpFac)
     ELSEIF(maxAng1.EQ.1)THEN  !(PP)
        call ICHOR_Ecoeffn_maxAngP2_maxAngA1_RHS(nPrimP,Ecoeffn,ETIJ,PreExpFac)
     ELSE !maxAng1.EQ.0        !(SD)
        call ICHOR_Ecoeffn_general_RHS(nPrimP,nTUV,ijk,maxAngP,maxAng1,maxAng2,&
             & Ecoeffn,ETIJ,PreExpFac)
     ENDIF
  ELSEIF(maxAngP.EQ.3)THEN 
     IF(maxAng1.EQ.3)THEN      !(FS)
        call ICHOR_Ecoeffn_maxAngP3_maxAngA3_RHS(nPrimP,Ecoeffn,ETIJ,PreExpFac)
     ELSEIF(maxAng1.EQ.2)THEN  !(DP)
        call ICHOR_Ecoeffn_maxAngP3_maxAngA2_RHS(nPrimP,Ecoeffn,ETIJ,PreExpFac)
     ELSEIF(maxAng1.EQ.1)THEN  !(PD)
        call ICHOR_Ecoeffn_maxAngP3_maxAngA1_RHS(nPrimP,Ecoeffn,ETIJ,PreExpFac)
     ELSE !maxAng1.EQ.0        !(SF)
        call ICHOR_Ecoeffn_maxAngP3_maxAngA0_RHS(nPrimP,Ecoeffn,ETIJ,PreExpFac)
     ENDIF
  ELSEIF(maxAngP.EQ.4)THEN 
     IF(maxAng1.EQ.3)THEN      !(FP)
        call ICHOR_Ecoeffn_maxAngP4_maxAngA3_RHS(nPrimP,Ecoeffn,ETIJ,PreExpFac)
     ELSEIF(maxAng1.EQ.2)THEN  !(DD)
        call ICHOR_Ecoeffn_maxAngP4_maxAngA2_RHS(nPrimP,Ecoeffn,ETIJ,PreExpFac)
     ELSEIF(maxAng1.EQ.1)THEN  !(PF)
        call ICHOR_Ecoeffn_maxAngP4_maxAngA1_RHS(nPrimP,Ecoeffn,ETIJ,PreExpFac)
     ELSE !(GS),(SG)
        call ICHOR_Ecoeffn_general_RHS(nPrimP,nTUV,ijk,maxAngP,maxAng1,maxAng2,&
            & Ecoeffn,ETIJ,PreExpFac)
     ENDIF
  ELSE
     call ICHOR_Ecoeffn_general_RHS(nPrimP,nTUV,ijk,maxAngP,maxAng1,maxAng2,&
          & Ecoeffn,ETIJ,PreExpFac)
  ENDIF
end Subroutine Ichorbuild_Ecoeff_RHS

Subroutine LS_OMP_DZERO(Ecoeffn,N)
implicit none
integer,intent(in) :: N
real(realk),intent(inout) :: Ecoeffn(n)
!local variables
integer :: I
!$OMP DO PRIVATE(I)
DO I=1,N
   Ecoeffn(I) = 0.0E0_realk
ENDDO
!$OMP END DO
end Subroutine LS_OMP_DZERO

Subroutine Ichorbuild_Ecoeff_LHS(nPrimP,nPrimA,nPrimB,maxAngP,maxAng1,maxAng2,nTUV,&
     & ijk,Aexp,Bexp,EcoeffN,Pdistance12,PreExpFac,nPasses,nAtomsA,nAtomsB,&
     & IatomApass,IatomBpass,MaxPasses,intprint,lupri,ETIJ,PA,PB,HPINV,TWOA,TWOB)
  implicit none
  integer,intent(in)        :: nPrimP,nPrimA,nPrimB,lupri,nTUV,ijk,MaxPasses
  integer,intent(in)        :: maxAngP,maxAng1,maxAng2,nPasses,intprint,nAtomsA,nAtomsB
  integer,intent(in)        :: IatomBpass(MaxPasses),IatomApass(MaxPasses)
  real(realk),intent(in)    :: Aexp(nPrimA),Bexp(nPrimB),Pdistance12(3,nAtomsA,nAtomsB)
  real(realk),intent(in)    :: PreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(inout) :: Ecoeffn(nPrimP,nPasses,nTUV,ijk)
  real(realk),intent(inout) :: ETIJ(nPrimP,nPasses,0:maxAngP,0:maxAng1,0:maxAng2,3)
  real(realk),intent(inout) :: PA(nPrimP,nPasses,3),PB(nPrimP,nPasses,3)
  real(realk),intent(inout) :: HPINV(nPrimP),TWOA(nPrimP)
  real(realk),intent(inout) :: TWOB(nPrimP)
  !
  real(realk),parameter :: DHALF=0.5E0_realk,D2=2.0E0_realk
  !
  integer :: K,iPass,iPrimP,iAtomA,iAtomB,X
  CALL LS_OMP_DZERO(Ecoeffn,nPrimP*nPasses*nTUV*ijk)
  call build_auxiliary_LHS(nPrimA,nPrimB,Aexp,Bexp,HPINV,TWOA,TWOB,PA,PB,Pdistance12,&
       & nPasses,MaxPasses,iAtomApass,iAtomBpass,nAtomsA,nAtomsB)
  IF(maxAng1+maxAng2.GT.6 )THEN
     CALL ICHOR_HERM_ECOEFFS_GENERAL_LHS(ETIJ,nPrimP,nPasses,maxAng1,maxAng2,&
          & PA,PB,HPINV,TWOA,TWOB,LUPRI)  
  ELSE !special less than 6
     CALL ICHOR_HERM_ECOEFFS(ETIJ,nPrimP,nPasses,maxAng1,maxAng2,&
          & PA,PB,HPINV,TWOA,TWOB,LUPRI)  
  ENDIF
  IF(intprint.GT.30)THEN
!$OMP MASTER
     do X=1,3
        call PRINT_ETIJ(ETIJ,nPrimP,nPasses,maxAng1,maxAng2,X,LUPRI)
     enddo
!$OMP END MASTER
!$OMP BARRIER
  ENDIF
  IF(maxAngP.EQ.0)THEN
     call ICHOR_Ecoeffn_maxAngP0_maxAngA0_LHS(nPrimP,nPasses,Ecoeffn,ETIJ,PreExpFac,&
          & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
  ELSEIF(maxAngP.EQ.1)THEN
     IF(maxAng1.EQ.1)THEN
        call ICHOR_Ecoeffn_maxAngP1_maxAngA1_LHS(nPrimP,nPasses,Ecoeffn,ETIJ,PreExpFac,&
             & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
     ELSEIF(maxAng1.EQ.0)THEN
        call ICHOR_Ecoeffn_maxAngP1_maxAngA0_LHS(nPrimP,nPasses,Ecoeffn,ETIJ,PreExpFac,&
             & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
     ENDIF
  ELSEIF(maxAngP.EQ.2)THEN 
     IF(maxAng1.EQ.2)THEN      !(DS)
        call ICHOR_Ecoeffn_maxAngP2_maxAngA2_LHS(nPrimP,nPasses,Ecoeffn,ETIJ,PreExpFac,&
             & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
     ELSEIF(maxAng1.EQ.1)THEN  !(PP)
        call ICHOR_Ecoeffn_maxAngP2_maxAngA1_LHS(nPrimP,nPasses,Ecoeffn,ETIJ,PreExpFac,&
             & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
     ELSE !maxAng1.EQ.0        !(SD)
        call ICHOR_Ecoeffn_maxAngP2_maxAngA0_LHS(nPrimP,nPasses,Ecoeffn,ETIJ,PreExpFac,&
             & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
     ENDIF
  ELSEIF(maxAngP.EQ.3)THEN 
     IF(maxAng1.EQ.3)THEN      !(FS)
        call ICHOR_Ecoeffn_maxAngP3_maxAngA3_LHS(nPrimP,nPasses,Ecoeffn,ETIJ,PreExpFac,&
             & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
     ELSEIF(maxAng1.EQ.2)THEN  !(DP)
        call ICHOR_Ecoeffn_maxAngP3_maxAngA2_LHS(nPrimP,nPasses,Ecoeffn,ETIJ,PreExpFac,&
             & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
     ELSEIF(maxAng1.EQ.1)THEN  !(PD)
        call ICHOR_Ecoeffn_maxAngP3_maxAngA1_LHS(nPrimP,nPasses,Ecoeffn,ETIJ,PreExpFac,&
             & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
     ELSE !maxAng1.EQ.0        !(SF)
        call ICHOR_Ecoeffn_maxAngP3_maxAngA0_LHS(nPrimP,nPasses,Ecoeffn,ETIJ,PreExpFac,&
             & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
     ENDIF
  ELSEIF(maxAngP.EQ.4)THEN 
     IF(maxAng1.EQ.3)THEN      !(FP)
        call ICHOR_Ecoeffn_maxAngP4_maxAngA3_LHS(nPrimP,nPasses,Ecoeffn,ETIJ,PreExpFac,&
             & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
     ELSEIF(maxAng1.EQ.2)THEN  !(DD)
        call ICHOR_Ecoeffn_maxAngP4_maxAngA2_LHS(nPrimP,nPasses,Ecoeffn,ETIJ,PreExpFac,&
             & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
     ELSEIF(maxAng1.EQ.1)THEN  !(PF)
        call ICHOR_Ecoeffn_maxAngP4_maxAngA1_LHS(nPrimP,nPasses,Ecoeffn,ETIJ,PreExpFac,&
             & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
     ELSE !(GS),(SG)
        call ICHOR_Ecoeffn_general_LHS(nPrimP,nPasses,nTUV,ijk,maxAngP,maxAng1,maxAng2,&
             & Ecoeffn,ETIJ,PreExpFac,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
     ENDIF
  ELSE
     call ICHOR_Ecoeffn_general_LHS(nPrimP,nPasses,nTUV,ijk,maxAngP,maxAng1,maxAng2,&
          & Ecoeffn,ETIJ,PreExpFac,MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
  ENDIF
end Subroutine Ichorbuild_Ecoeff_LHS
  
SUBROUTINE PRINT_ETIJ(ETIJ,nPrimP,nPasses,MAXI,MAXJ,X,LUPRI)
implicit none
INTEGER          ::  MAXI,MAXJ,X,nPrimP,nPasses,LUPRI
CHARACTER(len=4) ::  WORD
REAL(REALK)      ::  ETIJ(nPrimP,nPasses,0:MAXI+MAXJ,0:MAXI,0:MAXJ,3)
!
INTEGER          ::  I,J,T,K,IP
WRITE(LUPRI,*)'Output from Icor ETIJ'
WRITE(LUPRI,*)'----------------------'
WRITE (LUPRI,'(2X,A,2I5,I7,I7)') 'MAXI,MAXJ,nPrim,nPasses:', MAXI, MAXJ,nPrimP,nPasses
IF (X .EQ. 1) WORD = 'EX00'
IF (X .EQ. 2) WORD = 'EY00'
IF (X .EQ. 3) WORD = 'EZ00'
DO IP = 1, nPasses
   DO I = 0, MAXI 
      DO J = 0, MAXJ 
         DO T = 0, I + J
            WRITE (LUPRI,'(/,2X,A4,A1,I1,A1,I1,A1,I1,A7,I2,A1/)')&
                 &              WORD, '(', I, ',', J, ', ',  T, ',iPass=',iP,')'
            WRITE (LUPRI,'(1X,6ES16.8)') (ETIJ(K,IP,T,I,J,X),K=1,nPrimP)
         END DO
      END DO
   END DO
ENDDO
END SUBROUTINE PRINT_ETIJ

subroutine build_auxiliary_RHS(nPrimA,nPrimB,Aexp,Bexp,HPINV,TWOA,TWOB,PA,PB,Pdistance12)
  implicit none
  integer,intent(in) :: nPrimA,nPrimB
  real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB),Pdistance12(3)
  real(realk),intent(inout) :: HPINV(nPrimA,nPrimB),TWOA(nPrimA,nPrimB)
  real(realk),intent(inout) :: TWOB(nPrimA,nPrimB),PA(nPrimA,nPrimB,3)
  real(realk),intent(inout) :: PB(nPrimA,nPrimB,3)
  real(realk),parameter :: DHALF=0.5E0_realk,D2=2.0E0_realk,D1=1.0E0_realk
  !
  integer :: J,I,X
  !$OMP DO COLLAPSE(2) PRIVATE(I,J) 
  DO J=1,nPrimB
     DO I=1,nPrimA
        HPINV(I,J) = DHALF/(Aexp(I)+Bexp(J)) ! 1/2p
        TWOB(I,J)  = D2*Bexp(J)              ! 2b
        TWOA(I,J)  = D2*Aexp(I)              ! 2a
     ENDDO
  ENDDO
  !$OMP END DO NOWAIT
  !$OMP DO COLLAPSE(3) PRIVATE(X,I,J)
  DO X=1,3
     DO J=1,nPrimB
        DO I=1,nPrimA
           PA(I,J,X) = -Bexp(J)*Pdistance12(X)/(Aexp(I)+Bexp(J))
           PB(I,J,X) =  Aexp(I)*Pdistance12(X)/(Aexp(I)+Bexp(J))
        ENDDO
     ENDDO
  ENDDO
  !$OMP END DO
end subroutine build_auxiliary_RHS

subroutine build_auxiliary_LHS(nPrimA,nPrimB,Aexp,Bexp,HPINV,TWOA,TWOB,PA,PB,Pdistance12,&
     & nPasses,MaxPasses,iAtomApass,iAtomBpass,nAtomsA,nAtomsB)
  implicit none
  integer,intent(in) :: nPrimA,nPrimB,MaxPasses,nPasses,nAtomsA,nAtomsB
  real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB),Pdistance12(3,nAtomsA,nAtomsB)
  real(realk),intent(inout) :: HPINV(nPrimA,nPrimB),TWOA(nPrimA,nPrimB)
  real(realk),intent(inout) :: TWOB(nPrimA,nPrimB),PA(nPrimA,nPrimB,nPasses,3)
  real(realk),intent(inout) :: PB(nPrimA,nPrimB,nPasses,3)
  real(realk),parameter :: DHALF=0.5E0_realk,D2=2.0E0_realk,D1=1.0E0_realk
  integer,intent(in) :: iAtomApass(MaxPasses),iAtomBpass(MaxPasses)
  !
  integer :: J,I,X,IP,iAtomA,iAtomB
  !$OMP DO COLLAPSE(2) PRIVATE(I,J) 
  DO J=1,nPrimB
     DO I=1,nPrimA
        HPINV(I,J) = DHALF/(Aexp(I)+Bexp(J)) ! 1/2p
        TWOB(I,J)  = D2*Bexp(J)              ! 2b
        TWOA(I,J)  = D2*Aexp(I)              ! 2a
     ENDDO
  ENDDO
  !$OMP END DO NOWAIT
  !$OMP DO COLLAPSE(3) PRIVATE(X,I,J,iAtomA,iAtomB)
  DO X=1,3
     DO IP = 1,nPasses
        DO J=1,nPrimB
           iAtomA = iAtomApass(iP)
           iAtomB = iAtomBpass(iP)
           DO I=1,nPrimA
              PA(I,J,IP,X) = -Bexp(J)*Pdistance12(X,iAtomA,iAtomB)/(Aexp(I)+Bexp(J))
              PB(I,J,IP,X) =  Aexp(I)*Pdistance12(X,iAtomA,iAtomB)/(Aexp(I)+Bexp(J))
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !$OMP END DO
end subroutine build_auxiliary_LHS

!> \brief SUBROUTINE TO CALCULATE E's USING THE HERMITE RECURANNCE RELATIONS
!> \author A. Teale 
!> \date 2009
!>
!> Original code was written by Andrew Teale Andrew.Teale@nottingham.ac.uk
!> Copied and modified to Ichor code 2014.
!>
!> VARIABLES
!> =========      
!> nPrimP       NUMBER OF PRIMITIVES IN OVERLAP DISTRIBUTION
!> nPass        NUMBER OF PASSES OF OVERLAP DISTRIBUTIONS
!> MAXI,MAXJ    MAXIMUM I,J VALUES IN ETIJ's (MAX ANGMOM OF G_i and G_j RESPECTIVELY)
!> I,J,K,IJ     COUNTERS     
!> PA           DISTANCE BETWEEN CENTRE OF CHARGE P AND A (THE CENTRE OF PRIMITIVE GAUSSIAN G_i)           
!> PB           DISTANCE BETWEEN CENTRE OF CHARGE P AND B (THE CENTRE OF PRIMITIVE GAUSSIAN G_j)
!> HPINV        1/2p (p = a + b I.E. SUM OF EXPONENTS OF PRIM G_i and G_j)           
!> PINV         1/p
!> APINV, BPINV a/p b/p
!> TWOA, TWOB   2a, 2b
!> ETIJ         E COEFFICIENTS 
!> X            INTEGER SPECIFYING X,Y,Z (1,2,3) COMPONENT OF OVERALL OVERLAP DISTRIBUTION
!> T            COUNTER INDEX USED TO COUNT UPTO SUM OF ANGULAR MOMENTA OF G_i AND G_j 
!>
!> \param ETIJ the single cartesian Ecoeff to be built 
!> \param nprimP the number of primitives
!> \param nPass the number of npasses
!> \param MAXI the maximum angular momentum for orbital 1
!> \param MAXJ the maximum angular momentum for orbital 2
!> \param PA the argument -b/p*Rpq (b is exponent on orbital 2, p=a+b)
!> \param PB the argument  a/p*Rpq (a is exponent on orbital 1, p=a+b)
!> \param HPINV the argument 1/2p
!> \param TWOA 2a (a is exponent on orbital 1)
!> \param TWOB 2b (b is exponent on orbital 2)
!> \param X INTEGER SPECIFYING X,Y,Z (1,2,3) COMPONENT OF OVERALL OVERLAP DISTRIBUTION
!> \param LUPRI the logical unit number for the output file
SUBROUTINE ICHOR_HERM_ECOEFFS_GENERAL_LHS(ETIJ,nPrimP,nPass,MAXI,MAXJ,PA,PB,&
     & HPINV,TWOA,TWOB,LUPRI)
  implicit none
  INTEGER,intent(in)        :: nPrimP,nPass,MAXI,MAXJ,LUPRI
  REAL(REALK),intent(in)    :: PA(nPrimP,nPass,3), PB(nPrimP,nPass,3)
  REAL(REALK),intent(inout) :: ETIJ(nPrimP,nPass,0:MAXI+MAXJ,0:MAXI,0:MAXJ,3)
  REAL(REALK),intent(in)    :: HPINV(nPrimP)
  REAL(REALK)               :: TWOA(nPrimP),TWOB(nPrimP)
  !
  REAL(REALK),PARAMETER     :: D1 = 1.0E0_realk, D2 = 2.0E0_realk
  INTEGER                   :: MAXIJ,iPass,K,IJ,I,J,TJM,TIM,T,T1,X
!$OMP DO COLLAPSE(3) PRIVATE(X,Ipass,K,I,TIM,T,T1,J,TJM,IJ)
  DO X=1,3
   DO Ipass = 1,nPass
    DO K = 1, nPrimP
     !
     ! Run over I (J = 0)                                     
     ! ==================
     !
     DO I = 0, MAXI
        TIM = I-1
        IF (I .LE. 2) THEN
           !           E(0,0,0)                                              
           IF (I .EQ. 0) THEN
              ETIJ(K,ipass,0,0,0,X) = D1
           ELSE IF (I .EQ. 1) THEN
              !           E(T,1,0)                                            
              ETIJ(K,ipass,0,1,0,X) = PA(K,ipass,X)
              ETIJ(K,ipass,1,1,0,X) = HPINV(K)
           ELSE IF (I .EQ. 2) THEN
              !           E(T,2,0)                                              
              ETIJ(K,ipass,0,2,0,X) = PA(K,ipass,X)*PA(K,ipass,X) + HPINV(K) - D1/TWOA(K)
              ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass,X)*HPINV(K)
              ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
           ENDIF
        ELSE
           ETIJ(K,ipass,0,I,0,X) = PA(K,ipass,X)*ETIJ(K,ipass,0,I-1,0,X) + ETIJ(K,ipass,1,I-1,0,X) -&
                & TIM*ETIJ(K,ipass,0,I-2,0,X)/TWOA(K)
           ETIJ(K,ipass,I-1,I,0,X) = HPINV(K)*ETIJ(K,ipass,I-2,I-1,0,X) + &
                & PA(K,ipass,X)*ETIJ(K,ipass,I-1,I-1,0,X)
           ETIJ(K,ipass,I,I,0,X) = HPINV(K)*ETIJ(K,ipass,I-1,I-1,0,X)
           DO T = 1, I - 2                                    
              T1 = T + 1
              ETIJ(K,ipass,T,I,0,X) = HPINV(K)*ETIJ(K,ipass,T-1,I-1,0,X) &
                   & + PA(K,ipass,X)*ETIJ(K,ipass,  T,I-1,0,X)&
                   & + T1*ETIJ(K,ipass,T+1,I-1,0,X) - TIM*ETIJ(K,ipass, T,I-2,0,X)/TWOA(K)
           END DO
        END IF
        !   Run over J                                          
        !   ==========                                          
        DO J = 1, MAXJ
           IJ = I + J
           TJM = J-1
           IF (IJ .LE. 2) THEN
              IF (IJ .EQ. 1) THEN                                                                 
                 !              E(0,1)                                           
                 ETIJ(K,ipass,0,0,1,X) = PB(K,ipass,X)    
                 ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
              ELSE IF (IJ .EQ. 2) THEN     
                 !                 E(0,2)                        
                 IF (I .EQ. 0) THEN
                    ETIJ(K,ipass,0,0,2,X) = PB(K,ipass,X)*PB(K,ipass,X) + HPINV(K) - D1/TWOB(K)
                    ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass,X)*HPINV(K)
                    ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
                 ELSE
                    !                    E(1,1)                                    
                    ETIJ(K,ipass,0,1,1,X) = PA(K,ipass,X)*PB(K,ipass,X) + HPINV(K)
                    ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass,X) + PB(K,ipass,X))*HPINV(K)
                    ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
                 END IF
              ENDIF
           ELSE
              !             E(I,J)                                            
              IF(J.LT. 2)THEN
                 ETIJ(K,ipass,0,I,J,X) = PB(K,ipass,X)*ETIJ(K,ipass,0,I,J-1,X)+ETIJ(K,ipass,1,I,J-1,X)
                 ETIJ(K,ipass,IJ-1,I,J,X) = HPINV(K)*ETIJ(K,ipass,IJ-2,I,J-1,X)&
                      & + PB(K,ipass,X)*ETIJ(K,ipass,IJ-1,I,J-1,X)
                 ETIJ(K,ipass,IJ,I,J,X) = HPINV(K)*ETIJ(K,ipass,IJ-1,I,J-1,X)
                 DO T = 1, IJ - 2
                    T1 = T + 1
                    ETIJ(K,ipass,T,I,J,X)= HPINV(K)*ETIJ(K,ipass,T-1,I,J-1,X)&
                         & + PB(K,ipass,X)*ETIJ(K,ipass,  T,I,J-1,X)&
                         & + T1*ETIJ(K,ipass,T+1,I,J-1,X)! - TJM*ETIJ(K,ipass,  T,I,J-2)/TWOB(K)
                 END DO
              ELSE
                 ETIJ(K,ipass,   0,I,J,X) = PB(K,ipass,X)*ETIJ(K,ipass,   0,I,J-1,X) &
                      & + ETIJ(K,ipass,   1,I,J-1,X)&
                      & -TJM*ETIJ(K,ipass,  0,I,J-2,X)/TWOB(K)
                 ETIJ(K,ipass,IJ-1,I,J,X) = HPINV(K)*ETIJ(K,ipass,IJ-2,I,J-1,X) &
                      & + PB(K,ipass,X)*ETIJ(K,ipass,IJ-1,I,J-1,X)
                 ETIJ(K,ipass,  IJ,I,J,X) = HPINV(K)*ETIJ(K,ipass,IJ-1,I,J-1,X)
                 DO T = 1, IJ - 2
                    T1 = T + 1
                    ETIJ(K,ipass,T,I,J,X)= HPINV(K)*ETIJ(K,ipass,T-1,I,J-1,X) &
                         & + PB(K,ipass,X)*ETIJ(K,ipass,  T,I,J-1,X)&
                         & + T1*ETIJ(K,ipass,T+1,I,J-1,X) - TJM*ETIJ(K,ipass,  T,I,J-2,X)/TWOB(K)
                 END DO
              ENDIF
           END IF
        END DO
     END DO
    ENDDO
   ENDDO
  ENDDO
!$OMP END DO
END SUBROUTINE ICHOR_HERM_ECOEFFS_GENERAL_LHS

!> \brief SUBROUTINE TO CALCULATE E's USING THE HERMITE RECURANNCE RELATIONS
!> \author A. Teale 
!> \date 2009
!>
!> Original code was written by Andrew Teale Andrew.Teale@nottingham.ac.uk
!> Copied and modified to Ichor code 2014.
!>
!> VARIABLES
!> =========      
!> nPrimP       NUMBER OF PRIMITIVES IN OVERLAP DISTRIBUTION
!> MAXI,MAXJ    MAXIMUM I,J VALUES IN ETIJ's (MAX ANGMOM OF G_i and G_j RESPECTIVELY)
!> I,J,K,IJ     COUNTERS     
!> PA           DISTANCE BETWEEN CENTRE OF CHARGE P AND A (THE CENTRE OF PRIMITIVE GAUSSIAN G_i)           
!> PB           DISTANCE BETWEEN CENTRE OF CHARGE P AND B (THE CENTRE OF PRIMITIVE GAUSSIAN G_j)
!> HPINV        1/2p (p = a + b I.E. SUM OF EXPONENTS OF PRIM G_i and G_j)           
!> PINV         1/p
!> APINV, BPINV a/p b/p
!> TWOA, TWOB   2a, 2b
!> ETIJ         E COEFFICIENTS 
!> X            INTEGER SPECIFYING X,Y,Z (1,2,3) COMPONENT OF OVERALL OVERLAP DISTRIBUTION
!> T            COUNTER INDEX USED TO COUNT UPTO SUM OF ANGULAR MOMENTA OF G_i AND G_j 
!>
!> \param ETIJ the single cartesian Ecoeff to be built 
!> \param nprimP the number of primitives
!> \param MAXI the maximum angular momentum for orbital 1
!> \param MAXJ the maximum angular momentum for orbital 2
!> \param PA the argument -b/p*Rpq (b is exponent on orbital 2, p=a+b)
!> \param PB the argument  a/p*Rpq (a is exponent on orbital 1, p=a+b)
!> \param HPINV the argument 1/2p
!> \param TWOA 2a (a is exponent on orbital 1)
!> \param TWOB 2b (b is exponent on orbital 2)
!> \param X INTEGER SPECIFYING X,Y,Z (1,2,3) COMPONENT OF OVERALL OVERLAP DISTRIBUTION
!> \param LUPRI the logical unit number for the output file
SUBROUTINE ICHOR_HERM_ECOEFFS_GENERAL_RHS(ETIJ,nPrimP,MAXI,MAXJ,PA,PB,&
     & HPINV,TWOA,TWOB,LUPRI)
implicit none
INTEGER,intent(in)        :: nPrimP,MAXI,MAXJ,LUPRI
REAL(REALK),intent(in)    :: PA(nPrimP,3), PB(nPrimP,3)
REAL(REALK),intent(inout) :: ETIJ(nPrimP,0:MAXI+MAXJ,0:MAXI,0:MAXJ,3)
REAL(REALK),intent(in)    :: HPINV(nPrimP)
REAL(REALK)               :: TWOA(nPrimP),TWOB(nPrimP)
!
REAL(REALK),PARAMETER     :: D1 = 1.0E0_realk, D2 = 2.0E0_realk
INTEGER                   :: MAXIJ,K,IJ,I,J,TJM,TIM,T,T1,X
!
! Run over I (J = 0)                                     
! ==================
!
!$OMP DO COLLAPSE(2) PRIVATE(X,K,I,TIM,T,T1,J,TJM,IJ)
DO X=1,3
 DO K = 1, nPrimP
  DO I = 0, MAXI
   TIM = I-1
   IF (I .LE. 2) THEN
      !           E(0,0,0)                                              
      IF (I .EQ. 0) THEN
         ETIJ(K,0,0,0,X) = D1
      ELSE IF (I .EQ. 1) THEN
         !           E(T,1,0)                                            
         ETIJ(K,0,1,0,X) = PA(K,X)
         ETIJ(K,1,1,0,X) = HPINV(K)
      ELSE IF (I .EQ. 2) THEN
         !           E(T,2,0)                                              
         ETIJ(K,0,2,0,X) = PA(K,X)*PA(K,X) + HPINV(K) - D1/TWOA(K)
         ETIJ(K,1,2,0,X) = D2*PA(K,X)*HPINV(K)
         ETIJ(K,2,2,0,X) = HPINV(K)*HPINV(K)
      ENDIF
   ELSE
      ETIJ(K,0,I,0,X) = PA(K,X)*ETIJ(K,0,I-1,0,X) + ETIJ(K,1,I-1,0,X) - TIM*ETIJ(K,0,I-2,0,X)/TWOA(K)
      ETIJ(K,I-1,I,0,X) = HPINV(K)*ETIJ(K,I-2,I-1,0,X) + PA(K,X)*ETIJ(K,I-1,I-1,0,X)
      ETIJ(K,I,I,0,X) = HPINV(K)*ETIJ(K,I-1,I-1,0,X)
      DO T = 1, I - 2
         T1 = T + 1                                    
         ETIJ(K,T,I,0,X) = HPINV(K)*ETIJ(K,T-1,I-1,0,X)+ PA(K,X)*ETIJ(K,  T,I-1,0,X)&
              & + T1*ETIJ(K,T+1,I-1,0,X) - TIM*ETIJ(K, T,I-2,0,X)/TWOA(K)
      END DO
   END IF
   !   Run over J                                          
   !   ==========                                          
   DO J = 1, MAXJ
      IJ = I + J
      TJM = J-1     
      IF (IJ .LE. 2) THEN
         IF (IJ .EQ. 1) THEN                                                                 
            !              E(0,1)                                           
            ETIJ(K,0,0,1,X) = PB(K,X)    
            ETIJ(K,1,0,1,X) = HPINV(K) 
         ELSE IF (IJ .EQ. 2) THEN     
            !                 E(0,2)                        
            IF (I .EQ. 0) THEN
               ETIJ(K,0,0,2,X) = PB(K,X)*PB(K,X) + HPINV(K) - D1/TWOB(K)
               ETIJ(K,1,0,2,X) = D2*PB(K,X)*HPINV(K)
               ETIJ(K,2,0,2,X) = HPINV(K)*HPINV(K)
            ELSE
               !                    E(1,1)                                    
               ETIJ(K,0,1,1,X) = PA(K,X)*PB(K,X) + HPINV(K)
               ETIJ(K,1,1,1,X) = (PA(K,X) + PB(K,X))*HPINV(K)
               ETIJ(K,2,1,1,X) = HPINV(K)*HPINV(K)
            END IF
         ENDIF
      ELSE
         !             E(I,J)                                            
         IF(J.LT. 2)THEN
            ETIJ(K,0,I,J,X) = PB(K,X)*ETIJ(K,0,I,J-1,X)+ ETIJ(K,1,I,J-1,X)
            ETIJ(K,IJ-1,I,J,X) = HPINV(K)*ETIJ(K,IJ-2,I,J-1,X)+ PB(K,X)*ETIJ(K,IJ-1,I,J-1,X)
            ETIJ(K,IJ,I,J,X) = HPINV(K)*ETIJ(K,IJ-1,I,J-1,X)
            DO T = 1, IJ - 2
               T1 = T + 1                                    
               ETIJ(K,T,I,J,X)= HPINV(K)*ETIJ(K,T-1,I,J-1,X)+ PB(K,X)*ETIJ(K,  T,I,J-1,X)&
                    & + T1*ETIJ(K,T+1,I,J-1,X)! - TJM*ETIJ(K,  T,I,J-2)/TWOB(K)
            END DO
         ELSE
            ETIJ(K,   0,I,J,X) = PB(K,X)*ETIJ(K,   0,I,J-1,X)+ ETIJ(K,   1,I,J-1,X)&
                 & -TJM*ETIJ(K,  0,I,J-2,X)/TWOB(K)
            ETIJ(K,IJ-1,I,J,X) = HPINV(K)*ETIJ(K,IJ-2,I,J-1,X)+ PB(K,X)*ETIJ(K,IJ-1,I,J-1,X)
            ETIJ(K,  IJ,I,J,X) = HPINV(K)*ETIJ(K,IJ-1,I,J-1,X)
            DO T = 1, IJ - 2
               T1 = T + 1
               ETIJ(K,T,I,J,X)= HPINV(K)*ETIJ(K,T-1,I,J-1,X)+ PB(K,X)*ETIJ(K,  T,I,J-1,X)&
                    & + T1*ETIJ(K,T+1,I,J-1,X) - TJM*ETIJ(K,  T,I,J-2,X)/TWOB(K)
            END DO
         ENDIF
      END IF
   END DO
  END DO
 ENDDO
ENDDO
!$OMP END DO
END SUBROUTINE ICHOR_HERM_ECOEFFS_GENERAL_RHS

SUBROUTINE ICHOR_HERM_ECOEFFS(ETIJ,nPrimP,nPass,MAXI,MAXJ,PA,PB,&
     & HPINV,TWOA,TWOB,LUPRI)
  implicit none
  INTEGER,intent(in)        :: nPrimP,nPass,MAXI,MAXJ,LUPRI
  REAL(REALK),intent(in)    :: PA(nPrimP,nPass,3), PB(nPrimP,nPass,3)
  REAL(REALK),intent(inout) :: ETIJ(nPrimP,nPass,0:MAXI+MAXJ,0:MAXI,0:MAXJ,3)
  REAL(REALK),intent(in)    :: HPINV(nPrimP)
  REAL(REALK)               :: TWOA(nPrimP),TWOB(nPrimP)
  !
  REAL(REALK),PARAMETER     :: D1 = 1.0E0_realk, D2 = 2.0E0_realk
  INTEGER                   :: MAXIJ,iPass,K,X
!$OMP MASTER
  MAXIJ = MAXI + MAXJ
  IF(maxIJ.EQ.0 )THEN
     IF(maxI.EQ.0 )THEN
        !maxI.EQ.0.AND.maxJ.EQ.0
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ENDIF
  ELSEIF(maxIJ.EQ.1 )THEN
     IF(maxI.EQ.1 )THEN
        !maxI.EQ.1.AND.maxJ.EQ.0
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,1,0,X) = PA(K,ipass,X)
                 ETIJ(K,ipass,1,1,0,X) = HPINV(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ELSEIF(maxI.EQ.0)THEN
        !maxI.EQ.0.AND.maxJ.EQ.1
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,0,1,X) = PB(K,ipass,X)    
                 ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ENDIF
  ELSEIF(maxIJ.EQ.2 )THEN
     IF(maxI.EQ.2 )THEN
        !maxI.EQ.2.AND.maxJ.EQ.0
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,1,0,X) = PA(K,ipass,X)
                 ETIJ(K,ipass,1,1,0,X) = HPINV(K)
                 ETIJ(K,ipass,0,2,0,X) = PA(K,ipass,X)*PA(K,ipass,X) + HPINV(K) - D1/TWOA(K)
                 ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ELSEIF(maxI.EQ.1)THEN
        !maxI.EQ.1.AND.maxJ.EQ.1
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,0,1,X) = PB(K,ipass,X)    
                 ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
                 ETIJ(K,ipass,0,1,0,X) = PA(K,ipass,X)
                 ETIJ(K,ipass,1,1,0,X) = HPINV(K)
                 ETIJ(K,ipass,0,1,1,X) = PA(K,ipass,X)*PB(K,ipass,X) + HPINV(K)
                 ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass,X) + PB(K,ipass,X))*HPINV(K)
                 ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ELSEIF(maxI.EQ.0)THEN
        !maxI.EQ.0.AND.maxJ.EQ.2
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,0,1,X) = PB(K,ipass,X)    
                 ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
                 ETIJ(K,ipass,0,0,2,X) = PB(K,ipass,X)*PB(K,ipass,X) + HPINV(K) - D1/TWOB(K)
                 ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ENDIF
  ELSEIF(maxIJ.EQ.3 )THEN
     IF(maxI.EQ.3 )THEN
        !maxI.EQ.3.AND.maxJ.EQ.0
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,1,0,X) = PA(K,ipass,X)
                 ETIJ(K,ipass,1,1,0,X) = HPINV(K)
                 ETIJ(K,ipass,0,2,0,X) = PA(K,ipass,X)*PA(K,ipass,X) + HPINV(K) - D1/TWOA(K)
                 ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,3,0,X) = PA(K,ipass,X)*ETIJ(K,ipass,0,2,0,X)+ ETIJ(K,ipass,1,2,0,X)&
                      & - 2*ETIJ(K,ipass, 0,1,0,X)/TWOA(K)
                 ETIJ(K,ipass,2,3,0,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,3,3,0,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,1,3,0,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,1,2,0,X)&
                      & + 2*ETIJ(K,ipass,2,2,0,X) - 2*ETIJ(K,ipass, 1,1,0,X)/TWOA(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ELSEIF(maxI.EQ.2)THEN
        !maxI.EQ.2.AND.maxJ.EQ.1
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,0,1,X) = PB(K,ipass,X)    
                 ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
                 ETIJ(K,ipass,0,1,0,X) = PA(K,ipass,X)
                 ETIJ(K,ipass,1,1,0,X) = HPINV(K)
                 ETIJ(K,ipass,0,1,1,X) = PA(K,ipass,X)*PB(K,ipass,X) + HPINV(K)
                 ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass,X) + PB(K,ipass,X))*HPINV(K)
                 ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,2,0,X) = PA(K,ipass,X)*PA(K,ipass,X) + HPINV(K) - D1/TWOA(K)
                 ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,2,1,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,2,0,X)+ ETIJ(K,ipass, 1,2,0,X)
                 ETIJ(K,ipass,2,2,1,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,3,2,1,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,1,2,1,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,2,0,X)&
                      & + 2*ETIJ(K,ipass,2,2,0,X)! - 0*ETIJ(K,ipass,1,2,*,X)/TWOB(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ELSEIF(maxI.EQ.1)THEN
        !maxI.EQ.1.AND.maxJ.EQ.2
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,0,1,X) = PB(K,ipass,X)    
                 ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
                 ETIJ(K,ipass,0,0,2,X) = PB(K,ipass,X)*PB(K,ipass,X) + HPINV(K) - D1/TWOB(K)
                 ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,1,0,X) = PA(K,ipass,X)
                 ETIJ(K,ipass,1,1,0,X) = HPINV(K)
                 ETIJ(K,ipass,0,1,1,X) = PA(K,ipass,X)*PB(K,ipass,X) + HPINV(K)
                 ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass,X) + PB(K,ipass,X))*HPINV(K)
                 ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,1,2,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,1,1,X)+ ETIJ(K,ipass, 1,1,1,X)&
                      & -1*ETIJ(K,ipass,0,1,0,X)/TWOB(K)
                 ETIJ(K,ipass,2,1,2,X) = HPINV(K)*ETIJ(K,ipass,1,1,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,1,1,X)
                 ETIJ(K,ipass,3,1,2,X) = HPINV(K)*ETIJ(K,ipass,2,1,1,X)
                 ETIJ(K,ipass,1,1,2,X) = HPINV(K)*ETIJ(K,ipass,0,1,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,1,1,X)&
                      & + 2*ETIJ(K,ipass,2,1,1,X) - 1*ETIJ(K,ipass,1,1,0,X)/TWOB(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ELSEIF(maxI.EQ.0)THEN
        !maxI.EQ.0.AND.maxJ.EQ.3
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,0,1,X) = PB(K,ipass,X)    
                 ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
                 ETIJ(K,ipass,0,0,2,X) = PB(K,ipass,X)*PB(K,ipass,X) + HPINV(K) - D1/TWOB(K)
                 ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,0,3,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,0,2,X)+ ETIJ(K,ipass, 1,0,2,X)&
                      & -2*ETIJ(K,ipass,0,0,1,X)/TWOB(K)
                 ETIJ(K,ipass,2,0,3,X) = HPINV(K)*ETIJ(K,ipass,1,0,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,0,2,X)
                 ETIJ(K,ipass,3,0,3,X) = HPINV(K)*ETIJ(K,ipass,2,0,2,X)
                 ETIJ(K,ipass,1,0,3,X) = HPINV(K)*ETIJ(K,ipass,0,0,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,0,2,X)&
                      & + 2*ETIJ(K,ipass,2,0,2,X) - 2*ETIJ(K,ipass,1,0,1,X)/TWOB(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ENDIF
  ELSEIF(maxIJ.EQ.4 )THEN
     IF(maxI.EQ.4 )THEN
        !maxI.EQ.4.AND.maxJ.EQ.0
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,1,0,X) = PA(K,ipass,X)
                 ETIJ(K,ipass,1,1,0,X) = HPINV(K)
                 ETIJ(K,ipass,0,2,0,X) = PA(K,ipass,X)*PA(K,ipass,X) + HPINV(K) - D1/TWOA(K)
                 ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,3,0,X) = PA(K,ipass,X)*ETIJ(K,ipass,0,2,0,X)+ ETIJ(K,ipass,1,2,0,X)&
                      & - 2*ETIJ(K,ipass, 0,1,0,X)/TWOA(K)
                 ETIJ(K,ipass,2,3,0,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,3,3,0,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,1,3,0,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,1,2,0,X)&
                      & + 2*ETIJ(K,ipass,2,2,0,X) - 2*ETIJ(K,ipass, 1,1,0,X)/TWOA(K)
                 ETIJ(K,ipass,0,4,0,X) = PA(K,ipass,X)*ETIJ(K,ipass,0,3,0,X)+ ETIJ(K,ipass,1,3,0,X)&
                      & - 3*ETIJ(K,ipass, 0,2,0,X)/TWOA(K)
                 ETIJ(K,ipass,3,4,0,X) = HPINV(K)*ETIJ(K,ipass,2,3,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,3,3,0,X)
                 ETIJ(K,ipass,4,4,0,X) = HPINV(K)*ETIJ(K,ipass,3,3,0,X)
                 ETIJ(K,ipass,1,4,0,X) = HPINV(K)*ETIJ(K,ipass,0,3,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,1,3,0,X)&
                      & + 2*ETIJ(K,ipass,2,3,0,X) - 3*ETIJ(K,ipass, 1,2,0,X)/TWOA(K)
                 ETIJ(K,ipass,2,4,0,X) = HPINV(K)*ETIJ(K,ipass,1,3,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,2,3,0,X)&
                      & + 3*ETIJ(K,ipass,3,3,0,X) - 3*ETIJ(K,ipass, 2,2,0,X)/TWOA(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ELSEIF(maxI.EQ.3)THEN
        !maxI.EQ.3.AND.maxJ.EQ.1
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,0,1,X) = PB(K,ipass,X)    
                 ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
                 ETIJ(K,ipass,0,1,0,X) = PA(K,ipass,X)
                 ETIJ(K,ipass,1,1,0,X) = HPINV(K)
                 ETIJ(K,ipass,0,1,1,X) = PA(K,ipass,X)*PB(K,ipass,X) + HPINV(K)
                 ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass,X) + PB(K,ipass,X))*HPINV(K)
                 ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,2,0,X) = PA(K,ipass,X)*PA(K,ipass,X) + HPINV(K) - D1/TWOA(K)
                 ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,2,1,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,2,0,X)+ ETIJ(K,ipass, 1,2,0,X)
                 ETIJ(K,ipass,2,2,1,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,3,2,1,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,1,2,1,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,2,0,X)&
                      & + 2*ETIJ(K,ipass,2,2,0,X)! - 0*ETIJ(K,ipass,1,2,*,X)/TWOB(K)
                 ETIJ(K,ipass,0,3,0,X) = PA(K,ipass,X)*ETIJ(K,ipass,0,2,0,X)+ ETIJ(K,ipass,1,2,0,X)&
                      & - 2*ETIJ(K,ipass, 0,1,0,X)/TWOA(K)
                 ETIJ(K,ipass,2,3,0,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,3,3,0,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,1,3,0,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,1,2,0,X)&
                      & + 2*ETIJ(K,ipass,2,2,0,X) - 2*ETIJ(K,ipass, 1,1,0,X)/TWOA(K)
                 ETIJ(K,ipass,0,3,1,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,3,0,X)+ ETIJ(K,ipass, 1,3,0,X)
                 ETIJ(K,ipass,3,3,1,X) = HPINV(K)*ETIJ(K,ipass,2,3,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,3,0,X)
                 ETIJ(K,ipass,4,3,1,X) = HPINV(K)*ETIJ(K,ipass,3,3,0,X)
                 ETIJ(K,ipass,1,3,1,X) = HPINV(K)*ETIJ(K,ipass,0,3,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,3,0,X)&
                      & + 2*ETIJ(K,ipass,2,3,0,X)! - 0*ETIJ(K,ipass,1,3,*,X)/TWOB(K)
                 ETIJ(K,ipass,2,3,1,X) = HPINV(K)*ETIJ(K,ipass,1,3,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,3,0,X)&
                      & + 3*ETIJ(K,ipass,3,3,0,X)! - 0*ETIJ(K,ipass,2,3,*,X)/TWOB(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ELSEIF(maxI.EQ.2)THEN
        !maxI.EQ.2.AND.maxJ.EQ.2
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,0,1,X) = PB(K,ipass,X)    
                 ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
                 ETIJ(K,ipass,0,0,2,X) = PB(K,ipass,X)*PB(K,ipass,X) + HPINV(K) - D1/TWOB(K)
                 ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,1,0,X) = PA(K,ipass,X)
                 ETIJ(K,ipass,1,1,0,X) = HPINV(K)
                 ETIJ(K,ipass,0,1,1,X) = PA(K,ipass,X)*PB(K,ipass,X) + HPINV(K)
                 ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass,X) + PB(K,ipass,X))*HPINV(K)
                 ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)

                 ETIJ(K,ipass,0,1,2,X) = PB(K,ipass,X)*ETIJ(K,ipass,0,1,1,X)+ ETIJ(K,ipass,1,1,1,X)&
                      & -1*ETIJ(K,ipass,0,1,0,X)/TWOB(K)
                 ETIJ(K,ipass,2,1,2,X) = HPINV(K)*ETIJ(K,ipass,1,1,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,1,1,X)
                 ETIJ(K,ipass,3,1,2,X) = HPINV(K)*ETIJ(K,ipass,2,1,1,X)
                 ETIJ(K,ipass,1,1,2,X) = HPINV(K)*ETIJ(K,ipass,0,1,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,1,1,X)&
                      & + 2*ETIJ(K,ipass,2,1,1,X) - 1*ETIJ(K,ipass,1,1,0,X)/TWOB(K)
                 ETIJ(K,ipass,0,2,0,X) = PA(K,ipass,X)*PA(K,ipass,X) + HPINV(K) - D1/TWOA(K)
                 ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,2,1,X) = PB(K,ipass,X)*ETIJ(K,ipass,0,2,0,X)+ ETIJ(K,ipass,1,2,0,X)
                 ETIJ(K,ipass,2,2,1,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,3,2,1,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,1,2,1,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,2,0,X)&
                      & + 2*ETIJ(K,ipass,2,2,0,X)! - 0*ETIJ(K,ipass,1,2,*,X)/TWOB(K)
                 ETIJ(K,ipass,0,2,2,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,2,1,X)+ ETIJ(K,ipass, 1,2,1,X)&
                      & -1*ETIJ(K,ipass,0,2,0,X)/TWOB(K)
                 ETIJ(K,ipass,3,2,2,X) = HPINV(K)*ETIJ(K,ipass,2,2,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,2,1,X)
                 ETIJ(K,ipass,4,2,2,X) = HPINV(K)*ETIJ(K,ipass,3,2,1,X)
                 ETIJ(K,ipass,1,2,2,X) = HPINV(K)*ETIJ(K,ipass,0,2,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,2,1,X)&
                      & + 2*ETIJ(K,ipass,2,2,1,X) - 1*ETIJ(K,ipass,1,2,0,X)/TWOB(K)
                 ETIJ(K,ipass,2,2,2,X) = HPINV(K)*ETIJ(K,ipass,1,2,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,2,1,X)&
                      & + 3*ETIJ(K,ipass,3,2,1,X) - 1*ETIJ(K,ipass,2,2,0,X)/TWOB(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ELSEIF(maxI.EQ.1)THEN
        !maxI.EQ.1.AND.maxJ.EQ.3
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,0,1,X) = PB(K,ipass,X)    
                 ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
                 ETIJ(K,ipass,0,0,2,X) = PB(K,ipass,X)*PB(K,ipass,X) + HPINV(K) - D1/TWOB(K)
                 ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,0,3,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,0,2,X)+ ETIJ(K,ipass, 1,0,2,X)&
                      & -2*ETIJ(K,ipass,0,0,1,X)/TWOB(K)
                 ETIJ(K,ipass,2,0,3,X) = HPINV(K)*ETIJ(K,ipass,1,0,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,0,2,X)
                 ETIJ(K,ipass,3,0,3,X) = HPINV(K)*ETIJ(K,ipass,2,0,2,X)
                 ETIJ(K,ipass,1,0,3,X) = HPINV(K)*ETIJ(K,ipass,0,0,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,0,2,X)&
                      & + 2*ETIJ(K,ipass,2,0,2,X) - 2*ETIJ(K,ipass,1,0,1,X)/TWOB(K)
                 ETIJ(K,ipass,0,1,0,X) = PA(K,ipass,X)
                 ETIJ(K,ipass,1,1,0,X) = HPINV(K)
                 ETIJ(K,ipass,0,1,1,X) = PA(K,ipass,X)*PB(K,ipass,X) + HPINV(K)
                 ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass,X) + PB(K,ipass,X))*HPINV(K)
                 ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,1,2,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,1,1,X)+ ETIJ(K,ipass, 1,1,1,X)&
                      & -1*ETIJ(K,ipass,0,1,0,X)/TWOB(K)
                 ETIJ(K,ipass,2,1,2,X) = HPINV(K)*ETIJ(K,ipass,1,1,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,1,1,X)
                 ETIJ(K,ipass,3,1,2,X) = HPINV(K)*ETIJ(K,ipass,2,1,1,X)
                 ETIJ(K,ipass,1,1,2,X) = HPINV(K)*ETIJ(K,ipass,0,1,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,1,1,X)&
                      & + 2*ETIJ(K,ipass,2,1,1,X) - 1*ETIJ(K,ipass,1,1,0,X)/TWOB(K)
                 ETIJ(K,ipass,0,1,3,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,1,2,X)+ ETIJ(K,ipass, 1,1,2,X)&
                      & -2*ETIJ(K,ipass,0,1,1,X)/TWOB(K)
                 ETIJ(K,ipass,3,1,3,X) = HPINV(K)*ETIJ(K,ipass,2,1,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,1,2,X)
                 ETIJ(K,ipass,4,1,3,X) = HPINV(K)*ETIJ(K,ipass,3,1,2,X)
                 ETIJ(K,ipass,1,1,3,X) = HPINV(K)*ETIJ(K,ipass,0,1,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,1,2,X)&
                      & + 2*ETIJ(K,ipass,2,1,2,X) - 2*ETIJ(K,ipass,1,1,1,X)/TWOB(K)
                 ETIJ(K,ipass,2,1,3,X) = HPINV(K)*ETIJ(K,ipass,1,1,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,1,2,X)&
                      & + 3*ETIJ(K,ipass,3,1,2,X) - 2*ETIJ(K,ipass,2,1,1,X)/TWOB(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ELSEIF(maxI.EQ.0)THEN
        !maxI.EQ.0.AND.maxJ.EQ.4
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,0,1,X) = PB(K,ipass,X)    
                 ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
                 ETIJ(K,ipass,0,0,2,X) = PB(K,ipass,X)*PB(K,ipass,X) + HPINV(K) - D1/TWOB(K)
                 ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,0,3,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,0,2,X)+ ETIJ(K,ipass, 1,0,2,X)&
                      & -2*ETIJ(K,ipass,0,0,1,X)/TWOB(K)
                 ETIJ(K,ipass,2,0,3,X) = HPINV(K)*ETIJ(K,ipass,1,0,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,0,2,X)
                 ETIJ(K,ipass,3,0,3,X) = HPINV(K)*ETIJ(K,ipass,2,0,2,X)
                 ETIJ(K,ipass,1,0,3,X) = HPINV(K)*ETIJ(K,ipass,0,0,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,0,2,X)&
                      & + 2*ETIJ(K,ipass,2,0,2,X) - 2*ETIJ(K,ipass,1,0,1,X)/TWOB(K)
                 ETIJ(K,ipass,0,0,4,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,0,3,X)+ ETIJ(K,ipass, 1,0,3,X)&
                      & -3*ETIJ(K,ipass,0,0,2,X)/TWOB(K)
                 ETIJ(K,ipass,3,0,4,X) = HPINV(K)*ETIJ(K,ipass,2,0,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,0,3,X)
                 ETIJ(K,ipass,4,0,4,X) = HPINV(K)*ETIJ(K,ipass,3,0,3,X)
                 ETIJ(K,ipass,1,0,4,X) = HPINV(K)*ETIJ(K,ipass,0,0,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,0,3,X)&
                      & + 2*ETIJ(K,ipass,2,0,3,X) - 3*ETIJ(K,ipass,1,0,2,X)/TWOB(K)
                 ETIJ(K,ipass,2,0,4,X) = HPINV(K)*ETIJ(K,ipass,1,0,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,0,3,X)&
                      & + 3*ETIJ(K,ipass,3,0,3,X) - 3*ETIJ(K,ipass,2,0,2,X)/TWOB(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ENDIF
  ELSEIF(maxIJ.EQ.5 )THEN
     IF(maxI.EQ.5 )THEN
        !maxI.EQ.5.AND.maxJ.EQ.0
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,1,0,X) = PA(K,ipass,X)
                 ETIJ(K,ipass,1,1,0,X) = HPINV(K)
                 ETIJ(K,ipass,0,2,0,X) = PA(K,ipass,X)*PA(K,ipass,X) + HPINV(K) - D1/TWOA(K)
                 ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,3,0,X) = PA(K,ipass,X)*ETIJ(K,ipass,0,2,0,X)+ ETIJ(K,ipass,1,2,0,X)&
                      & - 2*ETIJ(K,ipass, 0,1,0,X)/TWOA(K)
                 ETIJ(K,ipass,2,3,0,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,3,3,0,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,1,3,0,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,1,2,0,X)&
                      & + 2*ETIJ(K,ipass,2,2,0,X) - 2*ETIJ(K,ipass, 1,1,0,X)/TWOA(K)
                 ETIJ(K,ipass,0,4,0,X) = PA(K,ipass,X)*ETIJ(K,ipass,0,3,0,X)+ ETIJ(K,ipass,1,3,0,X)&
                      & - 3*ETIJ(K,ipass, 0,2,0,X)/TWOA(K)
                 ETIJ(K,ipass,3,4,0,X) = HPINV(K)*ETIJ(K,ipass,2,3,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,3,3,0,X)
                 ETIJ(K,ipass,4,4,0,X) = HPINV(K)*ETIJ(K,ipass,3,3,0,X)
                 ETIJ(K,ipass,1,4,0,X) = HPINV(K)*ETIJ(K,ipass,0,3,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,1,3,0,X)&
                      & + 2*ETIJ(K,ipass,2,3,0,X) - 3*ETIJ(K,ipass, 1,2,0,X)/TWOA(K)
                 ETIJ(K,ipass,2,4,0,X) = HPINV(K)*ETIJ(K,ipass,1,3,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,2,3,0,X)&
                      & + 3*ETIJ(K,ipass,3,3,0,X) - 3*ETIJ(K,ipass, 2,2,0,X)/TWOA(K)
                 ETIJ(K,ipass,0,5,0,X) = PA(K,ipass,X)*ETIJ(K,ipass,0,4,0,X)+ ETIJ(K,ipass,1,4,0,X)&
                      & - 4*ETIJ(K,ipass, 0,3,0,X)/TWOA(K)
                 ETIJ(K,ipass,4,5,0,X) = HPINV(K)*ETIJ(K,ipass,3,4,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,4,4,0,X)
                 ETIJ(K,ipass,5,5,0,X) = HPINV(K)*ETIJ(K,ipass,4,4,0,X)
                 ETIJ(K,ipass,1,5,0,X) = HPINV(K)*ETIJ(K,ipass,0,4,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,1,4,0,X)&
                      & + 2*ETIJ(K,ipass,2,4,0,X) - 4*ETIJ(K,ipass, 1,3,0,X)/TWOA(K)
                 ETIJ(K,ipass,2,5,0,X) = HPINV(K)*ETIJ(K,ipass,1,4,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,2,4,0,X)&
                      & + 3*ETIJ(K,ipass,3,4,0,X) - 4*ETIJ(K,ipass, 2,3,0,X)/TWOA(K)
                 ETIJ(K,ipass,3,5,0,X) = HPINV(K)*ETIJ(K,ipass,2,4,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,3,4,0,X)&
                      & + 4*ETIJ(K,ipass,4,4,0,X) - 4*ETIJ(K,ipass, 3,3,0,X)/TWOA(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ELSEIF(maxI.EQ.4)THEN
        !maxI.EQ.4.AND.maxJ.EQ.1
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,0,1,X) = PB(K,ipass,X)    
                 ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
                 ETIJ(K,ipass,0,1,0,X) = PA(K,ipass,X)
                 ETIJ(K,ipass,1,1,0,X) = HPINV(K)
                 ETIJ(K,ipass,0,1,1,X) = PA(K,ipass,X)*PB(K,ipass,X) + HPINV(K)
                 ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass,X) + PB(K,ipass,X))*HPINV(K)
                 ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,2,0,X) = PA(K,ipass,X)*PA(K,ipass,X) + HPINV(K) - D1/TWOA(K)
                 ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,2,1,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,2,0,X)+ ETIJ(K,ipass, 1,2,0,X)
                 ETIJ(K,ipass,2,2,1,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,3,2,1,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,1,2,1,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,2,0,X)&
                      & + 2*ETIJ(K,ipass,2,2,0,X)! - 0*ETIJ(K,ipass,1,2,*,X)/TWOB(K)
                 ETIJ(K,ipass,0,3,0,X) = PA(K,ipass,X)*ETIJ(K,ipass,0,2,0,X)+ ETIJ(K,ipass,1,2,0,X)&
                      & - 2*ETIJ(K,ipass, 0,1,0,X)/TWOA(K)
                 ETIJ(K,ipass,2,3,0,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,3,3,0,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,1,3,0,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,1,2,0,X)&
                      & + 2*ETIJ(K,ipass,2,2,0,X) - 2*ETIJ(K,ipass, 1,1,0,X)/TWOA(K)
                 ETIJ(K,ipass,0,3,1,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,3,0,X)+ ETIJ(K,ipass, 1,3,0,X)
                 ETIJ(K,ipass,3,3,1,X) = HPINV(K)*ETIJ(K,ipass,2,3,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,3,0,X)
                 ETIJ(K,ipass,4,3,1,X) = HPINV(K)*ETIJ(K,ipass,3,3,0,X)
                 ETIJ(K,ipass,1,3,1,X)= HPINV(K)*ETIJ(K,ipass,0,3,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,3,0,X)&
                      & + 2*ETIJ(K,ipass,2,3,0,X)! - 0*ETIJ(K,ipass,1,3,*,X)/TWOB(K)
                 ETIJ(K,ipass,2,3,1,X)= HPINV(K)*ETIJ(K,ipass,1,3,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,3,0,X)&
                      & + 3*ETIJ(K,ipass,3,3,0,X)! - 0*ETIJ(K,ipass,2,3,*,X)/TWOB(K)
                 ETIJ(K,ipass,0,4,0,X) = PA(K,ipass,X)*ETIJ(K,ipass,0,3,0,X)+ ETIJ(K,ipass,1,3,0,X)&
                      & - 3*ETIJ(K,ipass, 0,2,0,X)/TWOA(K)
                 ETIJ(K,ipass,3,4,0,X) = HPINV(K)*ETIJ(K,ipass,2,3,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,3,3,0,X)
                 ETIJ(K,ipass,4,4,0,X) = HPINV(K)*ETIJ(K,ipass,3,3,0,X)
                 ETIJ(K,ipass,1,4,0,X) = HPINV(K)*ETIJ(K,ipass,0,3,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,1,3,0,X)&
                      & + 2*ETIJ(K,ipass,2,3,0,X) - 3*ETIJ(K,ipass, 1,2,0,X)/TWOA(K)
                 ETIJ(K,ipass,2,4,0,X) = HPINV(K)*ETIJ(K,ipass,1,3,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,2,3,0,X)&
                      & + 3*ETIJ(K,ipass,3,3,0,X) - 3*ETIJ(K,ipass, 2,2,0,X)/TWOA(K)
                 ETIJ(K,ipass,0,4,1,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,4,0,X)+ ETIJ(K,ipass, 1,4,0,X)
                 ETIJ(K,ipass,4,4,1,X) = HPINV(K)*ETIJ(K,ipass,3,4,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,4,4,0,X)
                 ETIJ(K,ipass,5,4,1,X) = HPINV(K)*ETIJ(K,ipass,4,4,0,X)
                 ETIJ(K,ipass,1,4,1,X) = HPINV(K)*ETIJ(K,ipass,0,4,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,4,0,X)&
                      & + 2*ETIJ(K,ipass,2,4,0,X)! - 0*ETIJ(K,ipass,1,4,*,X)/TWOB(K)
                 ETIJ(K,ipass,2,4,1,X) = HPINV(K)*ETIJ(K,ipass,1,4,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,4,0,X)&
                      & + 3*ETIJ(K,ipass,3,4,0,X)! - 0*ETIJ(K,ipass,2,4,*,X)/TWOB(K)
                 ETIJ(K,ipass,3,4,1,X) = HPINV(K)*ETIJ(K,ipass,2,4,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,4,0,X)&
                      & + 4*ETIJ(K,ipass,4,4,0,X)! - 0*ETIJ(K,ipass,3,4,*,X)/TWOB(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ELSEIF(maxI.EQ.3)THEN
        !maxI.EQ.3.AND.maxJ.EQ.2
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,0,1,X) = PB(K,ipass,X)    
                 ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
                 ETIJ(K,ipass,0,0,2,X) = PB(K,ipass,X)*PB(K,ipass,X) + HPINV(K) - D1/TWOB(K)
                 ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,1,0,X) = PA(K,ipass,X)
                 ETIJ(K,ipass,1,1,0,X) = HPINV(K)
                 ETIJ(K,ipass,0,1,1,X) = PA(K,ipass,X)*PB(K,ipass,X) + HPINV(K)
                 ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass,X) + PB(K,ipass,X))*HPINV(K)
                 ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,1,2,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,1,1,X)+ ETIJ(K,ipass, 1,1,1,X)&
                      & -1*ETIJ(K,ipass,0,1,0,X)/TWOB(K)
                 ETIJ(K,ipass,2,1,2,X) = HPINV(K)*ETIJ(K,ipass,1,1,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,1,1,X)
                 ETIJ(K,ipass,3,1,2,X) = HPINV(K)*ETIJ(K,ipass,2,1,1,X)
                 ETIJ(K,ipass,1,1,2,X)= HPINV(K)*ETIJ(K,ipass,0,1,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,1,1,X)&
                      & + 2*ETIJ(K,ipass,2,1,1,X) - 1*ETIJ(K,ipass,1,1,0,X)/TWOB(K)
                 ETIJ(K,ipass,0,2,0,X) = PA(K,ipass,X)*PA(K,ipass,X) + HPINV(K) - D1/TWOA(K)
                 ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,2,1,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,2,0,X)+ ETIJ(K,ipass, 1,2,0,X)
                 ETIJ(K,ipass,2,2,1,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,3,2,1,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,1,2,1,X)= HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,2,0,X)&
                      & + 2*ETIJ(K,ipass,2,2,0,X)! - 0*ETIJ(K,ipass,1,2,*,X)/TWOB(K)
                 ETIJ(K,ipass,0,2,2,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,2,1,X)+ ETIJ(K,ipass, 1,2,1,X)&
                      & -1*ETIJ(K,ipass,0,2,0,X)/TWOB(K)
                 ETIJ(K,ipass,3,2,2,X) = HPINV(K)*ETIJ(K,ipass,2,2,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,2,1,X)
                 ETIJ(K,ipass,4,2,2,X) = HPINV(K)*ETIJ(K,ipass,3,2,1,X)
                 ETIJ(K,ipass,1,2,2,X)= HPINV(K)*ETIJ(K,ipass,0,2,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,2,1,X)&
                      & + 2*ETIJ(K,ipass,2,2,1,X) - 1*ETIJ(K,ipass,1,2,0,X)/TWOB(K)
                 ETIJ(K,ipass,2,2,2,X)= HPINV(K)*ETIJ(K,ipass,1,2,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,2,1,X)&
                      & + 3*ETIJ(K,ipass,3,2,1,X) - 1*ETIJ(K,ipass,2,2,0,X)/TWOB(K)
                 ETIJ(K,ipass,0,3,0,X) = PA(K,ipass,X)*ETIJ(K,ipass,0,2,0,X)+ ETIJ(K,ipass,1,2,0,X)&
                      & - 2*ETIJ(K,ipass, 0,1,0,X)/TWOA(K)
                 ETIJ(K,ipass,2,3,0,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,3,3,0,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,1,3,0,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,1,2,0,X)&
                      & + 2*ETIJ(K,ipass,2,2,0,X) - 2*ETIJ(K,ipass, 1,1,0,X)/TWOA(K)
                 ETIJ(K,ipass,0,3,1,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,3,0,X)+ ETIJ(K,ipass, 1,3,0,X)
                 ETIJ(K,ipass,3,3,1,X) = HPINV(K)*ETIJ(K,ipass,2,3,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,3,0,X)
                 ETIJ(K,ipass,4,3,1,X) = HPINV(K)*ETIJ(K,ipass,3,3,0,X)
                 ETIJ(K,ipass,1,3,1,X)= HPINV(K)*ETIJ(K,ipass,0,3,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,3,0,X)&
                      & + 2*ETIJ(K,ipass,2,3,0,X)! - 0*ETIJ(K,ipass,1,3,*,X)/TWOB(K)
                 ETIJ(K,ipass,2,3,1,X)= HPINV(K)*ETIJ(K,ipass,1,3,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,3,0,X)&
                      & + 3*ETIJ(K,ipass,3,3,0,X)! - 0*ETIJ(K,ipass,2,3,*,X)/TWOB(K)
                 ETIJ(K,ipass,0,3,2,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,3,1,X)+ ETIJ(K,ipass, 1,3,1,X)&
                      & -1*ETIJ(K,ipass,0,3,0,X)/TWOB(K)
                 ETIJ(K,ipass,4,3,2,X) = HPINV(K)*ETIJ(K,ipass,3,3,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,4,3,1,X)
                 ETIJ(K,ipass,5,3,2,X) = HPINV(K)*ETIJ(K,ipass,4,3,1,X)
                 ETIJ(K,ipass,1,3,2,X)= HPINV(K)*ETIJ(K,ipass,0,3,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,3,1,X)&
                      & + 2*ETIJ(K,ipass,2,3,1,X) - 1*ETIJ(K,ipass,1,3,0,X)/TWOB(K)
                 ETIJ(K,ipass,2,3,2,X)= HPINV(K)*ETIJ(K,ipass,1,3,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,3,1,X)&
                      & + 3*ETIJ(K,ipass,3,3,1,X) - 1*ETIJ(K,ipass,2,3,0,X)/TWOB(K)
                 ETIJ(K,ipass,3,3,2,X)= HPINV(K)*ETIJ(K,ipass,2,3,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,3,1,X)&
                      & + 4*ETIJ(K,ipass,4,3,1,X) - 1*ETIJ(K,ipass,3,3,0,X)/TWOB(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ELSEIF(maxI.EQ.2)THEN
        !maxI.EQ.2.AND.maxJ.EQ.3
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,0,1,X) = PB(K,ipass,X)    
                 ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
                 ETIJ(K,ipass,0,0,2,X) = PB(K,ipass,X)*PB(K,ipass,X) + HPINV(K) - D1/TWOB(K)
                 ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,0,3,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,0,2,X)+ ETIJ(K,ipass, 1,0,2,X)&
                      & -2*ETIJ(K,ipass,0,0,1,X)/TWOB(K)
                 ETIJ(K,ipass,2,0,3,X) = HPINV(K)*ETIJ(K,ipass,1,0,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,0,2,X)
                 ETIJ(K,ipass,3,0,3,X) = HPINV(K)*ETIJ(K,ipass,2,0,2,X)
                 ETIJ(K,ipass,1,0,3,X)= HPINV(K)*ETIJ(K,ipass,0,0,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,0,2,X)&
                      & + 2*ETIJ(K,ipass,2,0,2,X) - 2*ETIJ(K,ipass,1,0,1,X)/TWOB(K)
                 ETIJ(K,ipass,0,1,0,X) = PA(K,ipass,X)
                 ETIJ(K,ipass,1,1,0,X) = HPINV(K)
                 ETIJ(K,ipass,0,1,1,X) = PA(K,ipass,X)*PB(K,ipass,X) + HPINV(K)
                 ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass,X) + PB(K,ipass,X))*HPINV(K)
                 ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,1,2,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,1,1,X)+ ETIJ(K,ipass, 1,1,1,X)&
                      & -1*ETIJ(K,ipass,0,1,0,X)/TWOB(K)
                 ETIJ(K,ipass,2,1,2,X) = HPINV(K)*ETIJ(K,ipass,1,1,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,1,1,X)
                 ETIJ(K,ipass,3,1,2,X) = HPINV(K)*ETIJ(K,ipass,2,1,1,X)
                 ETIJ(K,ipass,1,1,2,X)= HPINV(K)*ETIJ(K,ipass,0,1,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,1,1,X)&
                      & + 2*ETIJ(K,ipass,2,1,1,X) - 1*ETIJ(K,ipass,1,1,0,X)/TWOB(K)
                 ETIJ(K,ipass,0,1,3,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,1,2,X)+ ETIJ(K,ipass, 1,1,2,X)&
                      & -2*ETIJ(K,ipass,0,1,1,X)/TWOB(K)
                 ETIJ(K,ipass,3,1,3,X) = HPINV(K)*ETIJ(K,ipass,2,1,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,1,2,X)
                 ETIJ(K,ipass,4,1,3,X) = HPINV(K)*ETIJ(K,ipass,3,1,2,X)
                 ETIJ(K,ipass,1,1,3,X)= HPINV(K)*ETIJ(K,ipass,0,1,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,1,2,X)&
                      & + 2*ETIJ(K,ipass,2,1,2,X) - 2*ETIJ(K,ipass,1,1,1,X)/TWOB(K)
                 ETIJ(K,ipass,2,1,3,X)= HPINV(K)*ETIJ(K,ipass,1,1,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,1,2,X)&
                      & + 3*ETIJ(K,ipass,3,1,2,X) - 2*ETIJ(K,ipass,2,1,1,X)/TWOB(K)
                 ETIJ(K,ipass,0,2,0,X) = PA(K,ipass,X)*PA(K,ipass,X) + HPINV(K) - D1/TWOA(K)
                 ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,2,1,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,2,0,X)+ ETIJ(K,ipass, 1,2,0,X)
                 ETIJ(K,ipass,2,2,1,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,3,2,1,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,1,2,1,X)= HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,2,0,X)&
                      & + 2*ETIJ(K,ipass,2,2,0,X)! - 0*ETIJ(K,ipass,1,2,*,X)/TWOB(K)
                 ETIJ(K,ipass,0,2,2,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,2,1,X)+ ETIJ(K,ipass, 1,2,1,X)&
                      & -1*ETIJ(K,ipass,0,2,0,X)/TWOB(K)
                 ETIJ(K,ipass,3,2,2,X) = HPINV(K)*ETIJ(K,ipass,2,2,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,2,1,X)
                 ETIJ(K,ipass,4,2,2,X) = HPINV(K)*ETIJ(K,ipass,3,2,1,X)
                 ETIJ(K,ipass,1,2,2,X)= HPINV(K)*ETIJ(K,ipass,0,2,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,2,1,X)&
                      & + 2*ETIJ(K,ipass,2,2,1,X) - 1*ETIJ(K,ipass,1,2,0,X)/TWOB(K)
                 ETIJ(K,ipass,2,2,2,X)= HPINV(K)*ETIJ(K,ipass,1,2,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,2,1,X)&
                      & + 3*ETIJ(K,ipass,3,2,1,X) - 1*ETIJ(K,ipass,2,2,0,X)/TWOB(K)
                 ETIJ(K,ipass,0,2,3,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,2,2,X)+ ETIJ(K,ipass, 1,2,2,X)&
                      & -2*ETIJ(K,ipass,0,2,1,X)/TWOB(K)
                 ETIJ(K,ipass,4,2,3,X) = HPINV(K)*ETIJ(K,ipass,3,2,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,4,2,2,X)
                 ETIJ(K,ipass,5,2,3,X) = HPINV(K)*ETIJ(K,ipass,4,2,2,X)
                 ETIJ(K,ipass,1,2,3,X)= HPINV(K)*ETIJ(K,ipass,0,2,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,2,2,X)&
                      & + 2*ETIJ(K,ipass,2,2,2,X) - 2*ETIJ(K,ipass,1,2,1,X)/TWOB(K)
                 ETIJ(K,ipass,2,2,3,X)= HPINV(K)*ETIJ(K,ipass,1,2,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,2,2,X)&
                      & + 3*ETIJ(K,ipass,3,2,2,X) - 2*ETIJ(K,ipass,2,2,1,X)/TWOB(K)
                 ETIJ(K,ipass,3,2,3,X)= HPINV(K)*ETIJ(K,ipass,2,2,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,2,2,X)&
                      & + 4*ETIJ(K,ipass,4,2,2,X) - 2*ETIJ(K,ipass,3,2,1,X)/TWOB(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ELSEIF(maxI.EQ.1)THEN
        !maxI.EQ.1.AND.maxJ.EQ.4
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,0,1,X) = PB(K,ipass,X)    
                 ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
                 ETIJ(K,ipass,0,0,2,X) = PB(K,ipass,X)*PB(K,ipass,X) + HPINV(K) - D1/TWOB(K)
                 ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,0,3,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,0,2,X)+ ETIJ(K,ipass, 1,0,2,X)&
                      & -2*ETIJ(K,ipass,0,0,1,X)/TWOB(K)
                 ETIJ(K,ipass,2,0,3,X) = HPINV(K)*ETIJ(K,ipass,1,0,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,0,2,X)
                 ETIJ(K,ipass,3,0,3,X) = HPINV(K)*ETIJ(K,ipass,2,0,2,X)
                 ETIJ(K,ipass,1,0,3,X)= HPINV(K)*ETIJ(K,ipass,0,0,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,0,2,X)&
                      & + 2*ETIJ(K,ipass,2,0,2,X) - 2*ETIJ(K,ipass,1,0,1,X)/TWOB(K)
                 ETIJ(K,ipass,0,0,4,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,0,3,X)+ ETIJ(K,ipass, 1,0,3,X)&
                      & -3*ETIJ(K,ipass,0,0,2,X)/TWOB(K)
                 ETIJ(K,ipass,3,0,4,X) = HPINV(K)*ETIJ(K,ipass,2,0,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,0,3,X)
                 ETIJ(K,ipass,4,0,4,X) = HPINV(K)*ETIJ(K,ipass,3,0,3,X)
                 ETIJ(K,ipass,1,0,4,X)= HPINV(K)*ETIJ(K,ipass,0,0,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,0,3,X)&
                      & + 2*ETIJ(K,ipass,2,0,3,X) - 3*ETIJ(K,ipass,1,0,2,X)/TWOB(K)
                 ETIJ(K,ipass,2,0,4,X)= HPINV(K)*ETIJ(K,ipass,1,0,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,0,3,X)&
                      & + 3*ETIJ(K,ipass,3,0,3,X) - 3*ETIJ(K,ipass,2,0,2,X)/TWOB(K)
                 ETIJ(K,ipass,0,1,0,X) = PA(K,ipass,X)
                 ETIJ(K,ipass,1,1,0,X) = HPINV(K)
                 ETIJ(K,ipass,0,1,1,X) = PA(K,ipass,X)*PB(K,ipass,X) + HPINV(K)
                 ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass,X) + PB(K,ipass,X))*HPINV(K)
                 ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,1,2,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,1,1,X)+ ETIJ(K,ipass, 1,1,1,X)&
                      & -1*ETIJ(K,ipass,0,1,0,X)/TWOB(K)
                 ETIJ(K,ipass,2,1,2,X) = HPINV(K)*ETIJ(K,ipass,1,1,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,1,1,X)
                 ETIJ(K,ipass,3,1,2,X) = HPINV(K)*ETIJ(K,ipass,2,1,1,X)
                 ETIJ(K,ipass,1,1,2,X)= HPINV(K)*ETIJ(K,ipass,0,1,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,1,1,X)&
                      & + 2*ETIJ(K,ipass,2,1,1,X) - 1*ETIJ(K,ipass,1,1,0,X)/TWOB(K)
                 ETIJ(K,ipass,0,1,3,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,1,2,X)+ ETIJ(K,ipass, 1,1,2,X)&
                      & -2*ETIJ(K,ipass,0,1,1,X)/TWOB(K)
                 ETIJ(K,ipass,3,1,3,X) = HPINV(K)*ETIJ(K,ipass,2,1,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,1,2,X)
                 ETIJ(K,ipass,4,1,3,X) = HPINV(K)*ETIJ(K,ipass,3,1,2,X)
                 ETIJ(K,ipass,1,1,3,X)= HPINV(K)*ETIJ(K,ipass,0,1,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,1,2,X)&
                      & + 2*ETIJ(K,ipass,2,1,2,X) - 2*ETIJ(K,ipass,1,1,1,X)/TWOB(K)
                 ETIJ(K,ipass,2,1,3,X)= HPINV(K)*ETIJ(K,ipass,1,1,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,1,2,X)&
                      & + 3*ETIJ(K,ipass,3,1,2,X) - 2*ETIJ(K,ipass,2,1,1,X)/TWOB(K)
                 ETIJ(K,ipass,0,1,4,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,1,3,X)+ ETIJ(K,ipass, 1,1,3,X)&
                      & -3*ETIJ(K,ipass,0,1,2,X)/TWOB(K)
                 ETIJ(K,ipass,4,1,4,X) = HPINV(K)*ETIJ(K,ipass,3,1,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,4,1,3,X)
                 ETIJ(K,ipass,5,1,4,X) = HPINV(K)*ETIJ(K,ipass,4,1,3,X)
                 ETIJ(K,ipass,1,1,4,X)= HPINV(K)*ETIJ(K,ipass,0,1,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,1,3,X)&
                      & + 2*ETIJ(K,ipass,2,1,3,X) - 3*ETIJ(K,ipass,1,1,2,X)/TWOB(K)
                 ETIJ(K,ipass,2,1,4,X)= HPINV(K)*ETIJ(K,ipass,1,1,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,1,3,X)&
                      & + 3*ETIJ(K,ipass,3,1,3,X) - 3*ETIJ(K,ipass,2,1,2,X)/TWOB(K)
                 ETIJ(K,ipass,3,1,4,X)= HPINV(K)*ETIJ(K,ipass,2,1,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,1,3,X)&
                      & + 4*ETIJ(K,ipass,4,1,3,X) - 3*ETIJ(K,ipass,3,1,2,X)/TWOB(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ELSEIF(maxI.EQ.0)THEN
        !maxI.EQ.0.AND.maxJ.EQ.5
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,0,1,X) = PB(K,ipass,X)    
                 ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
                 ETIJ(K,ipass,0,0,2,X) = PB(K,ipass,X)*PB(K,ipass,X) + HPINV(K) - D1/TWOB(K)
                 ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,0,3,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,0,2,X)+ ETIJ(K,ipass, 1,0,2,X)&
                      & -2*ETIJ(K,ipass,0,0,1,X)/TWOB(K)
                 ETIJ(K,ipass,2,0,3,X) = HPINV(K)*ETIJ(K,ipass,1,0,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,0,2,X)
                 ETIJ(K,ipass,3,0,3,X) = HPINV(K)*ETIJ(K,ipass,2,0,2,X)
                 ETIJ(K,ipass,1,0,3,X)= HPINV(K)*ETIJ(K,ipass,0,0,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,0,2,X)&
                      & + 2*ETIJ(K,ipass,2,0,2,X) - 2*ETIJ(K,ipass,1,0,1,X)/TWOB(K)
                 ETIJ(K,ipass,0,0,4,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,0,3,X)+ ETIJ(K,ipass, 1,0,3,X)&
                      & -3*ETIJ(K,ipass,0,0,2,X)/TWOB(K)
                 ETIJ(K,ipass,3,0,4,X) = HPINV(K)*ETIJ(K,ipass,2,0,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,0,3,X)
                 ETIJ(K,ipass,4,0,4,X) = HPINV(K)*ETIJ(K,ipass,3,0,3,X)
                 ETIJ(K,ipass,1,0,4,X)= HPINV(K)*ETIJ(K,ipass,0,0,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,0,3,X)&
                      & + 2*ETIJ(K,ipass,2,0,3,X) - 3*ETIJ(K,ipass,1,0,2,X)/TWOB(K)
                 ETIJ(K,ipass,2,0,4,X)= HPINV(K)*ETIJ(K,ipass,1,0,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,0,3,X)&
                      & + 3*ETIJ(K,ipass,3,0,3,X) - 3*ETIJ(K,ipass,2,0,2,X)/TWOB(K)
                 ETIJ(K,ipass,0,0,5,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,0,4,X)+ ETIJ(K,ipass, 1,0,4,X)&
                      & -4*ETIJ(K,ipass,0,0,3,X)/TWOB(K)
                 ETIJ(K,ipass,4,0,5,X) = HPINV(K)*ETIJ(K,ipass,3,0,4,X)+ PB(K,ipass,X)*ETIJ(K,ipass,4,0,4,X)
                 ETIJ(K,ipass,5,0,5,X) = HPINV(K)*ETIJ(K,ipass,4,0,4,X)
                 ETIJ(K,ipass,1,0,5,X)= HPINV(K)*ETIJ(K,ipass,0,0,4,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,0,4,X)&
                      & + 2*ETIJ(K,ipass,2,0,4,X) - 4*ETIJ(K,ipass,1,0,3,X)/TWOB(K)
                 ETIJ(K,ipass,2,0,5,X)= HPINV(K)*ETIJ(K,ipass,1,0,4,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,0,4,X)&
                      & + 3*ETIJ(K,ipass,3,0,4,X) - 4*ETIJ(K,ipass,2,0,3,X)/TWOB(K)
                 ETIJ(K,ipass,3,0,5,X)= HPINV(K)*ETIJ(K,ipass,2,0,4,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,0,4,X)&
                      & + 4*ETIJ(K,ipass,4,0,4,X) - 4*ETIJ(K,ipass,3,0,3,X)/TWOB(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ENDIF
  ELSEIF(maxIJ.EQ.6 )THEN !FF
     IF(maxI.EQ.6 )THEN
        !maxI.EQ.6.AND.maxJ.EQ.0
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,1,0,X) = PA(K,ipass,X)
                 ETIJ(K,ipass,1,1,0,X) = HPINV(K)
                 ETIJ(K,ipass,0,2,0,X) = PA(K,ipass,X)*PA(K,ipass,X) + HPINV(K) - D1/TWOA(K)
                 ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,3,0,X) = PA(K,ipass,X)*ETIJ(K,ipass,0,2,0,X)+ ETIJ(K,ipass,1,2,0,X)&
                      & - 2*ETIJ(K,ipass, 0,1,0,X)/TWOA(K)
                 ETIJ(K,ipass,2,3,0,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,3,3,0,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,1,3,0,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,1,2,0,X)&
                      & + 2*ETIJ(K,ipass,2,2,0,X) - 2*ETIJ(K,ipass, 1,1,0,X)/TWOA(K)
                 ETIJ(K,ipass,0,4,0,X) = PA(K,ipass,X)*ETIJ(K,ipass,0,3,0,X)+ ETIJ(K,ipass,1,3,0,X)&
                      & - 3*ETIJ(K,ipass, 0,2,0,X)/TWOA(K)
                 ETIJ(K,ipass,3,4,0,X) = HPINV(K)*ETIJ(K,ipass,2,3,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,3,3,0,X)
                 ETIJ(K,ipass,4,4,0,X) = HPINV(K)*ETIJ(K,ipass,3,3,0,X)
                 ETIJ(K,ipass,1,4,0,X) = HPINV(K)*ETIJ(K,ipass,0,3,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,1,3,0,X)&
                      & + 2*ETIJ(K,ipass,2,3,0,X) - 3*ETIJ(K,ipass, 1,2,0,X)/TWOA(K)
                 ETIJ(K,ipass,2,4,0,X) = HPINV(K)*ETIJ(K,ipass,1,3,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,2,3,0,X)&
                      & + 3*ETIJ(K,ipass,3,3,0,X) - 3*ETIJ(K,ipass, 2,2,0,X)/TWOA(K)
                 ETIJ(K,ipass,0,5,0,X) = PA(K,ipass,X)*ETIJ(K,ipass,0,4,0,X)+ ETIJ(K,ipass,1,4,0,X)&
                      & - 4*ETIJ(K,ipass, 0,3,0,X)/TWOA(K)
                 ETIJ(K,ipass,4,5,0,X) = HPINV(K)*ETIJ(K,ipass,3,4,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,4,4,0,X)
                 ETIJ(K,ipass,5,5,0,X) = HPINV(K)*ETIJ(K,ipass,4,4,0,X)
                 ETIJ(K,ipass,1,5,0,X) = HPINV(K)*ETIJ(K,ipass,0,4,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,1,4,0,X)&
                      & + 2*ETIJ(K,ipass,2,4,0,X) - 4*ETIJ(K,ipass, 1,3,0,X)/TWOA(K)
                 ETIJ(K,ipass,2,5,0,X) = HPINV(K)*ETIJ(K,ipass,1,4,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,2,4,0,X)&
                      & + 3*ETIJ(K,ipass,3,4,0,X) - 4*ETIJ(K,ipass, 2,3,0,X)/TWOA(K)
                 ETIJ(K,ipass,3,5,0,X) = HPINV(K)*ETIJ(K,ipass,2,4,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,3,4,0,X)&
                      & + 4*ETIJ(K,ipass,4,4,0,X) - 4*ETIJ(K,ipass, 3,3,0,X)/TWOA(K)
                 ETIJ(K,ipass,0,6,0,X) = PA(K,ipass,X)*ETIJ(K,ipass,0,5,0,X)+ ETIJ(K,ipass,1,5,0,X)&
                      & - 5*ETIJ(K,ipass, 0,4,0,X)/TWOA(K)
                 ETIJ(K,ipass,5,6,0,X) = HPINV(K)*ETIJ(K,ipass,4,5,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,5,5,0,X)
                 ETIJ(K,ipass,6,6,0,X) = HPINV(K)*ETIJ(K,ipass,5,5,0,X)
                 ETIJ(K,ipass,1,6,0,X) = HPINV(K)*ETIJ(K,ipass,0,5,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,1,5,0,X)&
                      & + 2*ETIJ(K,ipass,2,5,0,X) - 5*ETIJ(K,ipass, 1,4,0,X)/TWOA(K)
                 ETIJ(K,ipass,2,6,0,X) = HPINV(K)*ETIJ(K,ipass,1,5,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,2,5,0,X)&
                      & + 3*ETIJ(K,ipass,3,5,0,X) - 5*ETIJ(K,ipass, 2,4,0,X)/TWOA(K)
                 ETIJ(K,ipass,3,6,0,X) = HPINV(K)*ETIJ(K,ipass,2,5,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,3,5,0,X)&
                      & + 4*ETIJ(K,ipass,4,5,0,X) - 5*ETIJ(K,ipass, 3,4,0,X)/TWOA(K)
                 ETIJ(K,ipass,4,6,0,X) = HPINV(K)*ETIJ(K,ipass,3,5,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,4,5,0,X)&
                      & + 5*ETIJ(K,ipass,5,5,0,X) - 5*ETIJ(K,ipass, 4,4,0,X)/TWOA(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ELSEIF(maxI.EQ.5)THEN
        !maxI.EQ.5.AND.maxJ.EQ.1
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,0,1,X) = PB(K,ipass,X)    
                 ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
                 ETIJ(K,ipass,0,1,0,X) = PA(K,ipass,X)
                 ETIJ(K,ipass,1,1,0,X) = HPINV(K)
                 ETIJ(K,ipass,0,1,1,X) = PA(K,ipass,X)*PB(K,ipass,X) + HPINV(K)
                 ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass,X) + PB(K,ipass,X))*HPINV(K)
                 ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,2,0,X) = PA(K,ipass,X)*PA(K,ipass,X) + HPINV(K) - D1/TWOA(K)
                 ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,2,1,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,2,0,X)+ ETIJ(K,ipass, 1,2,0,X)
                 ETIJ(K,ipass,2,2,1,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,3,2,1,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,1,2,1,X)= HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,2,0,X)&
                      & + 2*ETIJ(K,ipass,2,2,0,X)! - 0*ETIJ(K,ipass,1,2,*,X)/TWOB(K)
                 ETIJ(K,ipass,0,3,0,X) = PA(K,ipass,X)*ETIJ(K,ipass,0,2,0,X)+ ETIJ(K,ipass,1,2,0,X)&
                      & - 2*ETIJ(K,ipass, 0,1,0,X)/TWOA(K)
                 ETIJ(K,ipass,2,3,0,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,3,3,0,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,1,3,0,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,1,2,0,X)&
                      & + 2*ETIJ(K,ipass,2,2,0,X) - 2*ETIJ(K,ipass, 1,1,0,X)/TWOA(K)
                 ETIJ(K,ipass,0,3,1,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,3,0,X)+ ETIJ(K,ipass, 1,3,0,X)
                 ETIJ(K,ipass,3,3,1,X) = HPINV(K)*ETIJ(K,ipass,2,3,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,3,0,X)
                 ETIJ(K,ipass,4,3,1,X) = HPINV(K)*ETIJ(K,ipass,3,3,0,X)
                 ETIJ(K,ipass,1,3,1,X)= HPINV(K)*ETIJ(K,ipass,0,3,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,3,0,X)&
                      & + 2*ETIJ(K,ipass,2,3,0,X)! - 0*ETIJ(K,ipass,1,3,*,X)/TWOB(K)
                 ETIJ(K,ipass,2,3,1,X)= HPINV(K)*ETIJ(K,ipass,1,3,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,3,0,X)&
                      & + 3*ETIJ(K,ipass,3,3,0,X)! - 0*ETIJ(K,ipass,2,3,*,X)/TWOB(K)
                 ETIJ(K,ipass,0,4,0,X) = PA(K,ipass,X)*ETIJ(K,ipass,0,3,0,X)+ ETIJ(K,ipass,1,3,0,X)&
                      & - 3*ETIJ(K,ipass, 0,2,0,X)/TWOA(K)
                 ETIJ(K,ipass,3,4,0,X) = HPINV(K)*ETIJ(K,ipass,2,3,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,3,3,0,X)
                 ETIJ(K,ipass,4,4,0,X) = HPINV(K)*ETIJ(K,ipass,3,3,0,X)
                 ETIJ(K,ipass,1,4,0,X) = HPINV(K)*ETIJ(K,ipass,0,3,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,1,3,0,X)&
                      & + 2*ETIJ(K,ipass,2,3,0,X) - 3*ETIJ(K,ipass, 1,2,0,X)/TWOA(K)
                 ETIJ(K,ipass,2,4,0,X) = HPINV(K)*ETIJ(K,ipass,1,3,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,2,3,0,X)&
                      & + 3*ETIJ(K,ipass,3,3,0,X) - 3*ETIJ(K,ipass, 2,2,0,X)/TWOA(K)
                 ETIJ(K,ipass,0,4,1,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,4,0,X)+ ETIJ(K,ipass, 1,4,0,X)
                 ETIJ(K,ipass,4,4,1,X) = HPINV(K)*ETIJ(K,ipass,3,4,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,4,4,0,X)
                 ETIJ(K,ipass,5,4,1,X) = HPINV(K)*ETIJ(K,ipass,4,4,0,X)
                 ETIJ(K,ipass,1,4,1,X)= HPINV(K)*ETIJ(K,ipass,0,4,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,4,0,X)&
                      & + 2*ETIJ(K,ipass,2,4,0,X)! - 0*ETIJ(K,ipass,1,4,*,X)/TWOB(K)
                 ETIJ(K,ipass,2,4,1,X)= HPINV(K)*ETIJ(K,ipass,1,4,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,4,0,X)&
                      & + 3*ETIJ(K,ipass,3,4,0,X)! - 0*ETIJ(K,ipass,2,4,*,X)/TWOB(K)
                 ETIJ(K,ipass,3,4,1,X)= HPINV(K)*ETIJ(K,ipass,2,4,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,4,0,X)&
                      & + 4*ETIJ(K,ipass,4,4,0,X)! - 0*ETIJ(K,ipass,3,4,*,X)/TWOB(K)
                 ETIJ(K,ipass,0,5,0,X) = PA(K,ipass,X)*ETIJ(K,ipass,0,4,0,X)+ ETIJ(K,ipass,1,4,0,X)&
                      & - 4*ETIJ(K,ipass, 0,3,0,X)/TWOA(K)
                 ETIJ(K,ipass,4,5,0,X) = HPINV(K)*ETIJ(K,ipass,3,4,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,4,4,0,X)
                 ETIJ(K,ipass,5,5,0,X) = HPINV(K)*ETIJ(K,ipass,4,4,0,X)
                 ETIJ(K,ipass,1,5,0,X) = HPINV(K)*ETIJ(K,ipass,0,4,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,1,4,0,X)&
                      & + 2*ETIJ(K,ipass,2,4,0,X) - 4*ETIJ(K,ipass, 1,3,0,X)/TWOA(K)
                 ETIJ(K,ipass,2,5,0,X) = HPINV(K)*ETIJ(K,ipass,1,4,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,2,4,0,X)&
                      & + 3*ETIJ(K,ipass,3,4,0,X) - 4*ETIJ(K,ipass, 2,3,0,X)/TWOA(K)
                 ETIJ(K,ipass,3,5,0,X) = HPINV(K)*ETIJ(K,ipass,2,4,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,3,4,0,X)&
                      & + 4*ETIJ(K,ipass,4,4,0,X) - 4*ETIJ(K,ipass, 3,3,0,X)/TWOA(K)
                 ETIJ(K,ipass,0,5,1,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,5,0,X)+ ETIJ(K,ipass, 1,5,0,X)
                 ETIJ(K,ipass,5,5,1,X) = HPINV(K)*ETIJ(K,ipass,4,5,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,5,5,0,X)
                 ETIJ(K,ipass,6,5,1,X) = HPINV(K)*ETIJ(K,ipass,5,5,0,X)
                 ETIJ(K,ipass,1,5,1,X)= HPINV(K)*ETIJ(K,ipass,0,5,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,5,0,X)&
                      & + 2*ETIJ(K,ipass,2,5,0,X)! - 0*ETIJ(K,ipass,1,5,*,X)/TWOB(K)
                 ETIJ(K,ipass,2,5,1,X)= HPINV(K)*ETIJ(K,ipass,1,5,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,5,0,X)&
                      & + 3*ETIJ(K,ipass,3,5,0,X)! - 0*ETIJ(K,ipass,2,5,*,X)/TWOB(K)
                 ETIJ(K,ipass,3,5,1,X)= HPINV(K)*ETIJ(K,ipass,2,5,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,5,0,X)&
                      & + 4*ETIJ(K,ipass,4,5,0,X)! - 0*ETIJ(K,ipass,3,5,*,X)/TWOB(K)
                 ETIJ(K,ipass,4,5,1,X)= HPINV(K)*ETIJ(K,ipass,3,5,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,4,5,0,X)&
                      & + 5*ETIJ(K,ipass,5,5,0,X)! - 0*ETIJ(K,ipass,4,5,*,X)/TWOB(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ELSEIF(maxI.EQ.4)THEN
        !maxI.EQ.4.AND.maxJ.EQ.2
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,0,1,X) = PB(K,ipass,X)    
                 ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
                 ETIJ(K,ipass,0,0,2,X) = PB(K,ipass,X)*PB(K,ipass,X) + HPINV(K) - D1/TWOB(K)
                 ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,1,0,X) = PA(K,ipass,X)
                 ETIJ(K,ipass,1,1,0,X) = HPINV(K)
                 ETIJ(K,ipass,0,1,1,X) = PA(K,ipass,X)*PB(K,ipass,X) + HPINV(K)
                 ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass,X) + PB(K,ipass,X))*HPINV(K)
                 ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,1,2,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,1,1,X)+ ETIJ(K,ipass, 1,1,1,X)&
                      & -1*ETIJ(K,ipass,0,1,0,X)/TWOB(K)
                 ETIJ(K,ipass,2,1,2,X) = HPINV(K)*ETIJ(K,ipass,1,1,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,1,1,X)
                 ETIJ(K,ipass,3,1,2,X) = HPINV(K)*ETIJ(K,ipass,2,1,1,X)
                 ETIJ(K,ipass,1,1,2,X)= HPINV(K)*ETIJ(K,ipass,0,1,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,1,1,X)&
                      & + 2*ETIJ(K,ipass,2,1,1,X) - 1*ETIJ(K,ipass,1,1,0,X)/TWOB(K)
                 ETIJ(K,ipass,0,2,0,X) = PA(K,ipass,X)*PA(K,ipass,X) + HPINV(K) - D1/TWOA(K)
                 ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,2,1,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,2,0,X)+ ETIJ(K,ipass, 1,2,0,X)
                 ETIJ(K,ipass,2,2,1,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,3,2,1,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,1,2,1,X)= HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,2,0,X)&
                      & + 2*ETIJ(K,ipass,2,2,0,X)! - 0*ETIJ(K,ipass,1,2,*,X)/TWOB(K)
                 ETIJ(K,ipass,0,2,2,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,2,1,X)+ ETIJ(K,ipass, 1,2,1,X)&
                      & -1*ETIJ(K,ipass,0,2,0,X)/TWOB(K)
                 ETIJ(K,ipass,3,2,2,X) = HPINV(K)*ETIJ(K,ipass,2,2,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,2,1,X)
                 ETIJ(K,ipass,4,2,2,X) = HPINV(K)*ETIJ(K,ipass,3,2,1,X)
                 ETIJ(K,ipass,1,2,2,X)= HPINV(K)*ETIJ(K,ipass,0,2,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,2,1,X)&
                      & + 2*ETIJ(K,ipass,2,2,1,X) - 1*ETIJ(K,ipass,1,2,0,X)/TWOB(K)
                 ETIJ(K,ipass,2,2,2,X)= HPINV(K)*ETIJ(K,ipass,1,2,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,2,1,X)&
                      & + 3*ETIJ(K,ipass,3,2,1,X) - 1*ETIJ(K,ipass,2,2,0,X)/TWOB(K)
                 ETIJ(K,ipass,0,3,0,X) = PA(K,ipass,X)*ETIJ(K,ipass,0,2,0,X)+ ETIJ(K,ipass,1,2,0,X)&
                      & - 2*ETIJ(K,ipass, 0,1,0,X)/TWOA(K)
                 ETIJ(K,ipass,2,3,0,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,3,3,0,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,1,3,0,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,1,2,0,X)&
                      & + 2*ETIJ(K,ipass,2,2,0,X) - 2*ETIJ(K,ipass, 1,1,0,X)/TWOA(K)
                 ETIJ(K,ipass,0,3,1,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,3,0,X)+ ETIJ(K,ipass, 1,3,0,X)
                 ETIJ(K,ipass,3,3,1,X) = HPINV(K)*ETIJ(K,ipass,2,3,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,3,0,X)
                 ETIJ(K,ipass,4,3,1,X) = HPINV(K)*ETIJ(K,ipass,3,3,0,X)
                 ETIJ(K,ipass,1,3,1,X)= HPINV(K)*ETIJ(K,ipass,0,3,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,3,0,X)&
                      & + 2*ETIJ(K,ipass,2,3,0,X)! - 0*ETIJ(K,ipass,1,3,*,X)/TWOB(K)
                 ETIJ(K,ipass,2,3,1,X)= HPINV(K)*ETIJ(K,ipass,1,3,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,3,0,X)&
                      & + 3*ETIJ(K,ipass,3,3,0,X)! - 0*ETIJ(K,ipass,2,3,*,X)/TWOB(K)
                 ETIJ(K,ipass,0,3,2,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,3,1,X)+ ETIJ(K,ipass, 1,3,1,X)&
                      & -1*ETIJ(K,ipass,0,3,0,X)/TWOB(K)
                 ETIJ(K,ipass,4,3,2,X) = HPINV(K)*ETIJ(K,ipass,3,3,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,4,3,1,X)
                 ETIJ(K,ipass,5,3,2,X) = HPINV(K)*ETIJ(K,ipass,4,3,1,X)
                 ETIJ(K,ipass,1,3,2,X)= HPINV(K)*ETIJ(K,ipass,0,3,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,3,1,X)&
                      & + 2*ETIJ(K,ipass,2,3,1,X) - 1*ETIJ(K,ipass,1,3,0,X)/TWOB(K)
                 ETIJ(K,ipass,2,3,2,X)= HPINV(K)*ETIJ(K,ipass,1,3,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,3,1,X)&
                      & + 3*ETIJ(K,ipass,3,3,1,X) - 1*ETIJ(K,ipass,2,3,0,X)/TWOB(K)
                 ETIJ(K,ipass,3,3,2,X)= HPINV(K)*ETIJ(K,ipass,2,3,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,3,1,X)&
                      & + 4*ETIJ(K,ipass,4,3,1,X) - 1*ETIJ(K,ipass,3,3,0,X)/TWOB(K)
                 ETIJ(K,ipass,0,4,0,X) = PA(K,ipass,X)*ETIJ(K,ipass,0,3,0,X)+ ETIJ(K,ipass,1,3,0,X)&
                      & - 3*ETIJ(K,ipass, 0,2,0,X)/TWOA(K)
                 ETIJ(K,ipass,3,4,0,X) = HPINV(K)*ETIJ(K,ipass,2,3,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,3,3,0,X)
                 ETIJ(K,ipass,4,4,0,X) = HPINV(K)*ETIJ(K,ipass,3,3,0,X)
                 ETIJ(K,ipass,1,4,0,X) = HPINV(K)*ETIJ(K,ipass,0,3,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,1,3,0,X)&
                      & + 2*ETIJ(K,ipass,2,3,0,X) - 3*ETIJ(K,ipass, 1,2,0,X)/TWOA(K)
                 ETIJ(K,ipass,2,4,0,X) = HPINV(K)*ETIJ(K,ipass,1,3,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,2,3,0,X)&
                      & + 3*ETIJ(K,ipass,3,3,0,X) - 3*ETIJ(K,ipass, 2,2,0,X)/TWOA(K)
                 ETIJ(K,ipass,0,4,1,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,4,0,X)+ ETIJ(K,ipass, 1,4,0,X)
                 ETIJ(K,ipass,4,4,1,X) = HPINV(K)*ETIJ(K,ipass,3,4,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,4,4,0,X)
                 ETIJ(K,ipass,5,4,1,X) = HPINV(K)*ETIJ(K,ipass,4,4,0,X)
                 ETIJ(K,ipass,1,4,1,X)= HPINV(K)*ETIJ(K,ipass,0,4,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,4,0,X)&
                      & + 2*ETIJ(K,ipass,2,4,0,X)! - 0*ETIJ(K,ipass,1,4,*,X)/TWOB(K)
                 ETIJ(K,ipass,2,4,1,X)= HPINV(K)*ETIJ(K,ipass,1,4,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,4,0,X)&
                      & + 3*ETIJ(K,ipass,3,4,0,X)! - 0*ETIJ(K,ipass,2,4,*,X)/TWOB(K)
                 ETIJ(K,ipass,3,4,1,X)= HPINV(K)*ETIJ(K,ipass,2,4,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,4,0,X)&
                      & + 4*ETIJ(K,ipass,4,4,0,X)! - 0*ETIJ(K,ipass,3,4,*,X)/TWOB(K)
                 ETIJ(K,ipass,0,4,2,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,4,1,X)+ ETIJ(K,ipass, 1,4,1,X)&
                      & -1*ETIJ(K,ipass,0,4,0,X)/TWOB(K)
                 ETIJ(K,ipass,5,4,2,X) = HPINV(K)*ETIJ(K,ipass,4,4,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,5,4,1,X)
                 ETIJ(K,ipass,6,4,2,X) = HPINV(K)*ETIJ(K,ipass,5,4,1,X)
                 ETIJ(K,ipass,1,4,2,X)= HPINV(K)*ETIJ(K,ipass,0,4,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,4,1,X)&
                      & + 2*ETIJ(K,ipass,2,4,1,X) - 1*ETIJ(K,ipass,1,4,0,X)/TWOB(K)
                 ETIJ(K,ipass,2,4,2,X)= HPINV(K)*ETIJ(K,ipass,1,4,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,4,1,X)&
                      & + 3*ETIJ(K,ipass,3,4,1,X) - 1*ETIJ(K,ipass,2,4,0,X)/TWOB(K)
                 ETIJ(K,ipass,3,4,2,X)= HPINV(K)*ETIJ(K,ipass,2,4,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,4,1,X)&
                      & + 4*ETIJ(K,ipass,4,4,1,X) - 1*ETIJ(K,ipass,3,4,0,X)/TWOB(K)
                 ETIJ(K,ipass,4,4,2,X)= HPINV(K)*ETIJ(K,ipass,3,4,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,4,4,1,X)&
                      & + 5*ETIJ(K,ipass,5,4,1,X) - 1*ETIJ(K,ipass,4,4,0,X)/TWOB(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ELSEIF(maxI.EQ.3)THEN
        !maxI.EQ.3.AND.maxJ.EQ.3
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,0,1,X) = PB(K,ipass,X)    
                 ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
                 ETIJ(K,ipass,0,0,2,X) = PB(K,ipass,X)*PB(K,ipass,X) + HPINV(K) - D1/TWOB(K)
                 ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,0,3,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,0,2,X)+ ETIJ(K,ipass, 1,0,2,X)&
                      & -2*ETIJ(K,ipass,0,0,1,X)/TWOB(K)
                 ETIJ(K,ipass,2,0,3,X) = HPINV(K)*ETIJ(K,ipass,1,0,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,0,2,X)
                 ETIJ(K,ipass,3,0,3,X) = HPINV(K)*ETIJ(K,ipass,2,0,2,X)
                 ETIJ(K,ipass,1,0,3,X)= HPINV(K)*ETIJ(K,ipass,0,0,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,0,2,X)&
                      & + 2*ETIJ(K,ipass,2,0,2,X) - 2*ETIJ(K,ipass,1,0,1,X)/TWOB(K)
                 ETIJ(K,ipass,0,1,0,X) = PA(K,ipass,X)
                 ETIJ(K,ipass,1,1,0,X) = HPINV(K)
                 ETIJ(K,ipass,0,1,1,X) = PA(K,ipass,X)*PB(K,ipass,X) + HPINV(K)
                 ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass,X) + PB(K,ipass,X))*HPINV(K)
                 ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,1,2,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,1,1,X)+ ETIJ(K,ipass, 1,1,1,X)&
                      & -1*ETIJ(K,ipass,0,1,0,X)/TWOB(K)
                 ETIJ(K,ipass,2,1,2,X) = HPINV(K)*ETIJ(K,ipass,1,1,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,1,1,X)
                 ETIJ(K,ipass,3,1,2,X) = HPINV(K)*ETIJ(K,ipass,2,1,1,X)
                 ETIJ(K,ipass,1,1,2,X)= HPINV(K)*ETIJ(K,ipass,0,1,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,1,1,X)&
                      & + 2*ETIJ(K,ipass,2,1,1,X) - 1*ETIJ(K,ipass,1,1,0,X)/TWOB(K)
                 ETIJ(K,ipass,0,1,3,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,1,2,X)+ ETIJ(K,ipass, 1,1,2,X)&
                      & -2*ETIJ(K,ipass,0,1,1,X)/TWOB(K)
                 ETIJ(K,ipass,3,1,3,X) = HPINV(K)*ETIJ(K,ipass,2,1,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,1,2,X)
                 ETIJ(K,ipass,4,1,3,X) = HPINV(K)*ETIJ(K,ipass,3,1,2,X)
                 ETIJ(K,ipass,1,1,3,X)= HPINV(K)*ETIJ(K,ipass,0,1,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,1,2,X)&
                      & + 2*ETIJ(K,ipass,2,1,2,X) - 2*ETIJ(K,ipass,1,1,1,X)/TWOB(K)
                 ETIJ(K,ipass,2,1,3,X)= HPINV(K)*ETIJ(K,ipass,1,1,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,1,2,X)&
                      & + 3*ETIJ(K,ipass,3,1,2,X) - 2*ETIJ(K,ipass,2,1,1,X)/TWOB(K)
                 ETIJ(K,ipass,0,2,0,X) = PA(K,ipass,X)*PA(K,ipass,X) + HPINV(K) - D1/TWOA(K)
                 ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,2,1,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,2,0,X)+ ETIJ(K,ipass, 1,2,0,X)
                 ETIJ(K,ipass,2,2,1,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,3,2,1,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,1,2,1,X)= HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,2,0,X)&
                      & + 2*ETIJ(K,ipass,2,2,0,X)! - 0*ETIJ(K,ipass,1,2,*,X)/TWOB(K)
                 ETIJ(K,ipass,0,2,2,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,2,1,X)+ ETIJ(K,ipass, 1,2,1,X)&
                      & -1*ETIJ(K,ipass,0,2,0,X)/TWOB(K)
                 ETIJ(K,ipass,3,2,2,X) = HPINV(K)*ETIJ(K,ipass,2,2,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,2,1,X)
                 ETIJ(K,ipass,4,2,2,X) = HPINV(K)*ETIJ(K,ipass,3,2,1,X)
                 ETIJ(K,ipass,1,2,2,X)= HPINV(K)*ETIJ(K,ipass,0,2,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,2,1,X)&
                      & + 2*ETIJ(K,ipass,2,2,1,X) - 1*ETIJ(K,ipass,1,2,0,X)/TWOB(K)
                 ETIJ(K,ipass,2,2,2,X)= HPINV(K)*ETIJ(K,ipass,1,2,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,2,1,X)&
                      & + 3*ETIJ(K,ipass,3,2,1,X) - 1*ETIJ(K,ipass,2,2,0,X)/TWOB(K)
                 ETIJ(K,ipass,0,2,3,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,2,2,X)+ ETIJ(K,ipass, 1,2,2,X)&
                      & -2*ETIJ(K,ipass,0,2,1,X)/TWOB(K)
                 ETIJ(K,ipass,4,2,3,X) = HPINV(K)*ETIJ(K,ipass,3,2,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,4,2,2,X)
                 ETIJ(K,ipass,5,2,3,X) = HPINV(K)*ETIJ(K,ipass,4,2,2,X)
                 ETIJ(K,ipass,1,2,3,X)= HPINV(K)*ETIJ(K,ipass,0,2,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,2,2,X)&
                      & + 2*ETIJ(K,ipass,2,2,2,X) - 2*ETIJ(K,ipass,1,2,1,X)/TWOB(K)
                 ETIJ(K,ipass,2,2,3,X)= HPINV(K)*ETIJ(K,ipass,1,2,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,2,2,X)&
                      & + 3*ETIJ(K,ipass,3,2,2,X) - 2*ETIJ(K,ipass,2,2,1,X)/TWOB(K)
                 ETIJ(K,ipass,3,2,3,X)= HPINV(K)*ETIJ(K,ipass,2,2,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,2,2,X)&
                      & + 4*ETIJ(K,ipass,4,2,2,X) - 2*ETIJ(K,ipass,3,2,1,X)/TWOB(K)
                 ETIJ(K,ipass,0,3,0,X) = PA(K,ipass,X)*ETIJ(K,ipass,0,2,0,X)+ ETIJ(K,ipass,1,2,0,X)&
                      & - 2*ETIJ(K,ipass, 0,1,0,X)/TWOA(K)
                 ETIJ(K,ipass,2,3,0,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,3,3,0,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,1,3,0,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PA(K,ipass,X)*ETIJ(K,ipass,1,2,0,X)&
                      & + 2*ETIJ(K,ipass,2,2,0,X) - 2*ETIJ(K,ipass, 1,1,0,X)/TWOA(K)
                 ETIJ(K,ipass,0,3,1,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,3,0,X)+ ETIJ(K,ipass, 1,3,0,X)
                 ETIJ(K,ipass,3,3,1,X) = HPINV(K)*ETIJ(K,ipass,2,3,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,3,0,X)
                 ETIJ(K,ipass,4,3,1,X) = HPINV(K)*ETIJ(K,ipass,3,3,0,X)
                 ETIJ(K,ipass,1,3,1,X)= HPINV(K)*ETIJ(K,ipass,0,3,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,3,0,X)&
                      & + 2*ETIJ(K,ipass,2,3,0,X)! - 0*ETIJ(K,ipass,1,3,*,X)/TWOB(K)
                 ETIJ(K,ipass,2,3,1,X)= HPINV(K)*ETIJ(K,ipass,1,3,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,3,0,X)&
                      & + 3*ETIJ(K,ipass,3,3,0,X)! - 0*ETIJ(K,ipass,2,3,*,X)/TWOB(K)
                 ETIJ(K,ipass,0,3,2,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,3,1,X)+ ETIJ(K,ipass, 1,3,1,X)&
                      & -1*ETIJ(K,ipass,0,3,0,X)/TWOB(K)
                 ETIJ(K,ipass,4,3,2,X) = HPINV(K)*ETIJ(K,ipass,3,3,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,4,3,1,X)
                 ETIJ(K,ipass,5,3,2,X) = HPINV(K)*ETIJ(K,ipass,4,3,1,X)
                 ETIJ(K,ipass,1,3,2,X)= HPINV(K)*ETIJ(K,ipass,0,3,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,3,1,X)&
                      & + 2*ETIJ(K,ipass,2,3,1,X) - 1*ETIJ(K,ipass,1,3,0,X)/TWOB(K)
                 ETIJ(K,ipass,2,3,2,X)= HPINV(K)*ETIJ(K,ipass,1,3,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,3,1,X)&
                      & + 3*ETIJ(K,ipass,3,3,1,X) - 1*ETIJ(K,ipass,2,3,0,X)/TWOB(K)
                 ETIJ(K,ipass,3,3,2,X)= HPINV(K)*ETIJ(K,ipass,2,3,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,3,1,X)&
                      & + 4*ETIJ(K,ipass,4,3,1,X) - 1*ETIJ(K,ipass,3,3,0,X)/TWOB(K)
                 ETIJ(K,ipass,0,3,3,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,3,2,X)+ ETIJ(K,ipass, 1,3,2,X)&
                      & -2*ETIJ(K,ipass,0,3,1,X)/TWOB(K)
                 ETIJ(K,ipass,5,3,3,X) = HPINV(K)*ETIJ(K,ipass,4,3,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,5,3,2,X)
                 ETIJ(K,ipass,6,3,3,X) = HPINV(K)*ETIJ(K,ipass,5,3,2,X)
                 ETIJ(K,ipass,1,3,3,X)= HPINV(K)*ETIJ(K,ipass,0,3,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,3,2,X)&
                      & + 2*ETIJ(K,ipass,2,3,2,X) - 2*ETIJ(K,ipass,1,3,1,X)/TWOB(K)
                 ETIJ(K,ipass,2,3,3,X)= HPINV(K)*ETIJ(K,ipass,1,3,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,3,2,X)&
                      & + 3*ETIJ(K,ipass,3,3,2,X) - 2*ETIJ(K,ipass,2,3,1,X)/TWOB(K)
                 ETIJ(K,ipass,3,3,3,X)= HPINV(K)*ETIJ(K,ipass,2,3,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,3,2,X)&
                      & + 4*ETIJ(K,ipass,4,3,2,X) - 2*ETIJ(K,ipass,3,3,1,X)/TWOB(K)
                 ETIJ(K,ipass,4,3,3,X)= HPINV(K)*ETIJ(K,ipass,3,3,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,4,3,2,X)&
                      & + 5*ETIJ(K,ipass,5,3,2,X) - 2*ETIJ(K,ipass,4,3,1,X)/TWOB(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO 
     ELSEIF(maxI.EQ.2)THEN
        !maxI.EQ.2.AND.maxJ.EQ.4
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,0,1,X) = PB(K,ipass,X)    
                 ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
                 ETIJ(K,ipass,0,0,2,X) = PB(K,ipass,X)*PB(K,ipass,X) + HPINV(K) - D1/TWOB(K)
                 ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,0,3,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,0,2,X)+ ETIJ(K,ipass, 1,0,2,X)&
                      & -2*ETIJ(K,ipass,0,0,1,X)/TWOB(K)
                 ETIJ(K,ipass,2,0,3,X) = HPINV(K)*ETIJ(K,ipass,1,0,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,0,2,X)
                 ETIJ(K,ipass,3,0,3,X) = HPINV(K)*ETIJ(K,ipass,2,0,2,X)
                 ETIJ(K,ipass,1,0,3,X)= HPINV(K)*ETIJ(K,ipass,0,0,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,0,2,X)&
                      & + 2*ETIJ(K,ipass,2,0,2,X) - 2*ETIJ(K,ipass,1,0,1,X)/TWOB(K)
                 ETIJ(K,ipass,0,0,4,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,0,3,X)+ ETIJ(K,ipass, 1,0,3,X)&
                      & -3*ETIJ(K,ipass,0,0,2,X)/TWOB(K)
                 ETIJ(K,ipass,3,0,4,X) = HPINV(K)*ETIJ(K,ipass,2,0,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,0,3,X)
                 ETIJ(K,ipass,4,0,4,X) = HPINV(K)*ETIJ(K,ipass,3,0,3,X)
                 ETIJ(K,ipass,1,0,4,X)= HPINV(K)*ETIJ(K,ipass,0,0,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,0,3,X)&
                      & + 2*ETIJ(K,ipass,2,0,3,X) - 3*ETIJ(K,ipass,1,0,2,X)/TWOB(K)
                 ETIJ(K,ipass,2,0,4,X)= HPINV(K)*ETIJ(K,ipass,1,0,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,0,3,X)&
                      & + 3*ETIJ(K,ipass,3,0,3,X) - 3*ETIJ(K,ipass,2,0,2,X)/TWOB(K)
                 ETIJ(K,ipass,0,1,0,X) = PA(K,ipass,X)
                 ETIJ(K,ipass,1,1,0,X) = HPINV(K)
                 ETIJ(K,ipass,0,1,1,X) = PA(K,ipass,X)*PB(K,ipass,X) + HPINV(K)
                 ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass,X) + PB(K,ipass,X))*HPINV(K)
                 ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,1,2,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,1,1,X)+ ETIJ(K,ipass, 1,1,1,X)&
                      & -1*ETIJ(K,ipass,0,1,0,X)/TWOB(K)
                 ETIJ(K,ipass,2,1,2,X) = HPINV(K)*ETIJ(K,ipass,1,1,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,1,1,X)
                 ETIJ(K,ipass,3,1,2,X) = HPINV(K)*ETIJ(K,ipass,2,1,1,X)
                 ETIJ(K,ipass,1,1,2,X)= HPINV(K)*ETIJ(K,ipass,0,1,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,1,1,X)&
                      & + 2*ETIJ(K,ipass,2,1,1,X) - 1*ETIJ(K,ipass,1,1,0,X)/TWOB(K)
                 ETIJ(K,ipass,0,1,3,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,1,2,X)+ ETIJ(K,ipass, 1,1,2,X)&
                      & -2*ETIJ(K,ipass,0,1,1,X)/TWOB(K)
                 ETIJ(K,ipass,3,1,3,X) = HPINV(K)*ETIJ(K,ipass,2,1,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,1,2,X)
                 ETIJ(K,ipass,4,1,3,X) = HPINV(K)*ETIJ(K,ipass,3,1,2,X)
                 ETIJ(K,ipass,1,1,3,X)= HPINV(K)*ETIJ(K,ipass,0,1,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,1,2,X)&
                      & + 2*ETIJ(K,ipass,2,1,2,X) - 2*ETIJ(K,ipass,1,1,1,X)/TWOB(K)
                 ETIJ(K,ipass,2,1,3,X)= HPINV(K)*ETIJ(K,ipass,1,1,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,1,2,X)&
                      & + 3*ETIJ(K,ipass,3,1,2,X) - 2*ETIJ(K,ipass,2,1,1,X)/TWOB(K)
                 ETIJ(K,ipass,0,1,4,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,1,3,X)+ ETIJ(K,ipass, 1,1,3,X)&
                      & -3*ETIJ(K,ipass,0,1,2,X)/TWOB(K)
                 ETIJ(K,ipass,4,1,4,X) = HPINV(K)*ETIJ(K,ipass,3,1,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,4,1,3,X)
                 ETIJ(K,ipass,5,1,4,X) = HPINV(K)*ETIJ(K,ipass,4,1,3,X)
                 ETIJ(K,ipass,1,1,4,X)= HPINV(K)*ETIJ(K,ipass,0,1,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,1,3,X)&
                      & + 2*ETIJ(K,ipass,2,1,3,X) - 3*ETIJ(K,ipass,1,1,2,X)/TWOB(K)
                 ETIJ(K,ipass,2,1,4,X)= HPINV(K)*ETIJ(K,ipass,1,1,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,1,3,X)&
                      & + 3*ETIJ(K,ipass,3,1,3,X) - 3*ETIJ(K,ipass,2,1,2,X)/TWOB(K)
                 ETIJ(K,ipass,3,1,4,X)= HPINV(K)*ETIJ(K,ipass,2,1,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,1,3,X)&
                      & + 4*ETIJ(K,ipass,4,1,3,X) - 3*ETIJ(K,ipass,3,1,2,X)/TWOB(K)
                 ETIJ(K,ipass,0,2,0,X) = PA(K,ipass,X)*PA(K,ipass,X) + HPINV(K) - D1/TWOA(K)
                 ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,2,1,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,2,0,X)+ ETIJ(K,ipass, 1,2,0,X)
                 ETIJ(K,ipass,2,2,1,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,3,2,1,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
                 ETIJ(K,ipass,1,2,1,X)= HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,2,0,X)&
                      & + 2*ETIJ(K,ipass,2,2,0,X)! - 0*ETIJ(K,ipass,1,2,*,X)/TWOB(K)
                 ETIJ(K,ipass,0,2,2,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,2,1,X)+ ETIJ(K,ipass, 1,2,1,X)&
                      & -1*ETIJ(K,ipass,0,2,0,X)/TWOB(K)
                 ETIJ(K,ipass,3,2,2,X) = HPINV(K)*ETIJ(K,ipass,2,2,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,2,1,X)
                 ETIJ(K,ipass,4,2,2,X) = HPINV(K)*ETIJ(K,ipass,3,2,1,X)
                 ETIJ(K,ipass,1,2,2,X)= HPINV(K)*ETIJ(K,ipass,0,2,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,2,1,X)&
                      & + 2*ETIJ(K,ipass,2,2,1,X) - 1*ETIJ(K,ipass,1,2,0,X)/TWOB(K)
                 ETIJ(K,ipass,2,2,2,X)= HPINV(K)*ETIJ(K,ipass,1,2,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,2,1,X)&
                      & + 3*ETIJ(K,ipass,3,2,1,X) - 1*ETIJ(K,ipass,2,2,0,X)/TWOB(K)
                 ETIJ(K,ipass,0,2,3,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,2,2,X)+ ETIJ(K,ipass, 1,2,2,X)&
                      & -2*ETIJ(K,ipass,0,2,1,X)/TWOB(K)
                 ETIJ(K,ipass,4,2,3,X) = HPINV(K)*ETIJ(K,ipass,3,2,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,4,2,2,X)
                 ETIJ(K,ipass,5,2,3,X) = HPINV(K)*ETIJ(K,ipass,4,2,2,X)
                 ETIJ(K,ipass,1,2,3,X)= HPINV(K)*ETIJ(K,ipass,0,2,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,2,2,X)&
                      & + 2*ETIJ(K,ipass,2,2,2,X) - 2*ETIJ(K,ipass,1,2,1,X)/TWOB(K)
                 ETIJ(K,ipass,2,2,3,X)= HPINV(K)*ETIJ(K,ipass,1,2,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,2,2,X)&
                      & + 3*ETIJ(K,ipass,3,2,2,X) - 2*ETIJ(K,ipass,2,2,1,X)/TWOB(K)
                 ETIJ(K,ipass,3,2,3,X)= HPINV(K)*ETIJ(K,ipass,2,2,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,2,2,X)&
                      & + 4*ETIJ(K,ipass,4,2,2,X) - 2*ETIJ(K,ipass,3,2,1,X)/TWOB(K)
                 ETIJ(K,ipass,0,2,4,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,2,3,X)+ ETIJ(K,ipass, 1,2,3,X)&
                      & -3*ETIJ(K,ipass,0,2,2,X)/TWOB(K)
                 ETIJ(K,ipass,5,2,4,X) = HPINV(K)*ETIJ(K,ipass,4,2,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,5,2,3,X)
                 ETIJ(K,ipass,6,2,4,X) = HPINV(K)*ETIJ(K,ipass,5,2,3,X)
                 ETIJ(K,ipass,1,2,4,X)= HPINV(K)*ETIJ(K,ipass,0,2,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,2,3,X)&
                      & + 2*ETIJ(K,ipass,2,2,3,X) - 3*ETIJ(K,ipass,1,2,2,X)/TWOB(K)
                 ETIJ(K,ipass,2,2,4,X)= HPINV(K)*ETIJ(K,ipass,1,2,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,2,3,X)&
                      & + 3*ETIJ(K,ipass,3,2,3,X) - 3*ETIJ(K,ipass,2,2,2,X)/TWOB(K)
                 ETIJ(K,ipass,3,2,4,X)= HPINV(K)*ETIJ(K,ipass,2,2,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,2,3,X)&
                      & + 4*ETIJ(K,ipass,4,2,3,X) - 3*ETIJ(K,ipass,3,2,2,X)/TWOB(K)
                 ETIJ(K,ipass,4,2,4,X)= HPINV(K)*ETIJ(K,ipass,3,2,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,4,2,3,X)&
                      & + 5*ETIJ(K,ipass,5,2,3,X) - 3*ETIJ(K,ipass,4,2,2,X)/TWOB(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO 
     ELSEIF(maxI.EQ.1)THEN
        !maxI.EQ.1.AND.maxJ.EQ.5
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,0,1,X) = PB(K,ipass,X)    
                 ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
                 ETIJ(K,ipass,0,0,2,X) = PB(K,ipass,X)*PB(K,ipass,X) + HPINV(K) - D1/TWOB(K)
                 ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,0,3,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,0,2,X)+ ETIJ(K,ipass, 1,0,2,X)&
                      & -2*ETIJ(K,ipass,0,0,1,X)/TWOB(K)
                 ETIJ(K,ipass,2,0,3,X) = HPINV(K)*ETIJ(K,ipass,1,0,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,0,2,X)
                 ETIJ(K,ipass,3,0,3,X) = HPINV(K)*ETIJ(K,ipass,2,0,2,X)
                 ETIJ(K,ipass,1,0,3,X)= HPINV(K)*ETIJ(K,ipass,0,0,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,0,2,X)&
                      & + 2*ETIJ(K,ipass,2,0,2,X) - 2*ETIJ(K,ipass,1,0,1,X)/TWOB(K)
                 ETIJ(K,ipass,0,0,4,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,0,3,X)+ ETIJ(K,ipass, 1,0,3,X)&
                      & -3*ETIJ(K,ipass,0,0,2,X)/TWOB(K)
                 ETIJ(K,ipass,3,0,4,X) = HPINV(K)*ETIJ(K,ipass,2,0,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,0,3,X)
                 ETIJ(K,ipass,4,0,4,X) = HPINV(K)*ETIJ(K,ipass,3,0,3,X)
                 ETIJ(K,ipass,1,0,4,X)= HPINV(K)*ETIJ(K,ipass,0,0,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,0,3,X)&
                      & + 2*ETIJ(K,ipass,2,0,3,X) - 3*ETIJ(K,ipass,1,0,2,X)/TWOB(K)
                 ETIJ(K,ipass,2,0,4,X)= HPINV(K)*ETIJ(K,ipass,1,0,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,0,3,X)&
                      & + 3*ETIJ(K,ipass,3,0,3,X) - 3*ETIJ(K,ipass,2,0,2,X)/TWOB(K)
                 ETIJ(K,ipass,0,0,5,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,0,4,X)+ ETIJ(K,ipass, 1,0,4,X)&
                      & -4*ETIJ(K,ipass,0,0,3,X)/TWOB(K)
                 ETIJ(K,ipass,4,0,5,X) = HPINV(K)*ETIJ(K,ipass,3,0,4,X)+ PB(K,ipass,X)*ETIJ(K,ipass,4,0,4,X)
                 ETIJ(K,ipass,5,0,5,X) = HPINV(K)*ETIJ(K,ipass,4,0,4,X)
                 ETIJ(K,ipass,1,0,5,X)= HPINV(K)*ETIJ(K,ipass,0,0,4,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,0,4,X)&
                      & + 2*ETIJ(K,ipass,2,0,4,X) - 4*ETIJ(K,ipass,1,0,3,X)/TWOB(K)
                 ETIJ(K,ipass,2,0,5,X)= HPINV(K)*ETIJ(K,ipass,1,0,4,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,0,4,X)&
                      & + 3*ETIJ(K,ipass,3,0,4,X) - 4*ETIJ(K,ipass,2,0,3,X)/TWOB(K)
                 ETIJ(K,ipass,3,0,5,X)= HPINV(K)*ETIJ(K,ipass,2,0,4,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,0,4,X)&
                      & + 4*ETIJ(K,ipass,4,0,4,X) - 4*ETIJ(K,ipass,3,0,3,X)/TWOB(K)
                 ETIJ(K,ipass,0,1,0,X) = PA(K,ipass,X)
                 ETIJ(K,ipass,1,1,0,X) = HPINV(K)
                 ETIJ(K,ipass,0,1,1,X) = PA(K,ipass,X)*PB(K,ipass,X) + HPINV(K)
                 ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass,X) + PB(K,ipass,X))*HPINV(K)
                 ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,1,2,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,1,1,X)+ ETIJ(K,ipass, 1,1,1,X)&
                      & -1*ETIJ(K,ipass,0,1,0,X)/TWOB(K)
                 ETIJ(K,ipass,2,1,2,X) = HPINV(K)*ETIJ(K,ipass,1,1,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,1,1,X)
                 ETIJ(K,ipass,3,1,2,X) = HPINV(K)*ETIJ(K,ipass,2,1,1,X)
                 ETIJ(K,ipass,1,1,2,X)= HPINV(K)*ETIJ(K,ipass,0,1,1,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,1,1,X)&
                      & + 2*ETIJ(K,ipass,2,1,1,X) - 1*ETIJ(K,ipass,1,1,0,X)/TWOB(K)
                 ETIJ(K,ipass,0,1,3,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,1,2,X)+ ETIJ(K,ipass, 1,1,2,X)&
                      & -2*ETIJ(K,ipass,0,1,1,X)/TWOB(K)
                 ETIJ(K,ipass,3,1,3,X) = HPINV(K)*ETIJ(K,ipass,2,1,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,1,2,X)
                 ETIJ(K,ipass,4,1,3,X) = HPINV(K)*ETIJ(K,ipass,3,1,2,X)
                 ETIJ(K,ipass,1,1,3,X)= HPINV(K)*ETIJ(K,ipass,0,1,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,1,2,X)&
                      & + 2*ETIJ(K,ipass,2,1,2,X) - 2*ETIJ(K,ipass,1,1,1,X)/TWOB(K)
                 ETIJ(K,ipass,2,1,3,X)= HPINV(K)*ETIJ(K,ipass,1,1,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,1,2,X)&
                      & + 3*ETIJ(K,ipass,3,1,2,X) - 2*ETIJ(K,ipass,2,1,1,X)/TWOB(K)
                 ETIJ(K,ipass,0,1,4,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,1,3,X)+ ETIJ(K,ipass, 1,1,3,X)&
                      & -3*ETIJ(K,ipass,0,1,2,X)/TWOB(K)
                 ETIJ(K,ipass,4,1,4,X) = HPINV(K)*ETIJ(K,ipass,3,1,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,4,1,3,X)
                 ETIJ(K,ipass,5,1,4,X) = HPINV(K)*ETIJ(K,ipass,4,1,3,X)
                 ETIJ(K,ipass,1,1,4,X)= HPINV(K)*ETIJ(K,ipass,0,1,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,1,3,X)&
                      & + 2*ETIJ(K,ipass,2,1,3,X) - 3*ETIJ(K,ipass,1,1,2,X)/TWOB(K)
                 ETIJ(K,ipass,2,1,4,X)= HPINV(K)*ETIJ(K,ipass,1,1,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,1,3,X)&
                      & + 3*ETIJ(K,ipass,3,1,3,X) - 3*ETIJ(K,ipass,2,1,2,X)/TWOB(K)
                 ETIJ(K,ipass,3,1,4,X)= HPINV(K)*ETIJ(K,ipass,2,1,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,1,3,X)&
                      & + 4*ETIJ(K,ipass,4,1,3,X) - 3*ETIJ(K,ipass,3,1,2,X)/TWOB(K)
                 ETIJ(K,ipass,0,1,5,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,1,4,X)+ ETIJ(K,ipass, 1,1,4,X)&
                      & -4*ETIJ(K,ipass,0,1,3,X)/TWOB(K)
                 ETIJ(K,ipass,5,1,5,X) = HPINV(K)*ETIJ(K,ipass,4,1,4,X)+ PB(K,ipass,X)*ETIJ(K,ipass,5,1,4,X)
                 ETIJ(K,ipass,6,1,5,X) = HPINV(K)*ETIJ(K,ipass,5,1,4,X)
                 ETIJ(K,ipass,1,1,5,X)= HPINV(K)*ETIJ(K,ipass,0,1,4,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,1,4,X)&
                      & + 2*ETIJ(K,ipass,2,1,4,X) - 4*ETIJ(K,ipass,1,1,3,X)/TWOB(K)
                 ETIJ(K,ipass,2,1,5,X)= HPINV(K)*ETIJ(K,ipass,1,1,4,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,1,4,X)&
                      & + 3*ETIJ(K,ipass,3,1,4,X) - 4*ETIJ(K,ipass,2,1,3,X)/TWOB(K)
                 ETIJ(K,ipass,3,1,5,X)= HPINV(K)*ETIJ(K,ipass,2,1,4,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,1,4,X)&
                      & + 4*ETIJ(K,ipass,4,1,4,X) - 4*ETIJ(K,ipass,3,1,3,X)/TWOB(K)
                 ETIJ(K,ipass,4,1,5,X)= HPINV(K)*ETIJ(K,ipass,3,1,4,X)+ PB(K,ipass,X)*ETIJ(K,ipass,4,1,4,X)&
                      & + 5*ETIJ(K,ipass,5,1,4,X) - 4*ETIJ(K,ipass,4,1,3,X)/TWOB(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ELSEIF(maxI.EQ.0)THEN
        !maxI.EQ.0.AND.maxJ.EQ.6
        !!$OMP DO  PRIVATE(X,Ipass,K)
        DO X=1,3
           DO Ipass = 1,nPass
              DO K = 1, nPrimP
                 ETIJ(K,ipass,0,0,0,X) = D1
                 ETIJ(K,ipass,0,0,1,X) = PB(K,ipass,X)    
                 ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
                 ETIJ(K,ipass,0,0,2,X) = PB(K,ipass,X)*PB(K,ipass,X) + HPINV(K) - D1/TWOB(K)
                 ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass,X)*HPINV(K)
                 ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
                 ETIJ(K,ipass,0,0,3,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,0,2,X)+ ETIJ(K,ipass, 1,0,2,X)&
                      & -2*ETIJ(K,ipass,0,0,1,X)/TWOB(K)
                 ETIJ(K,ipass,2,0,3,X) = HPINV(K)*ETIJ(K,ipass,1,0,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,0,2,X)
                 ETIJ(K,ipass,3,0,3,X) = HPINV(K)*ETIJ(K,ipass,2,0,2,X)
                 ETIJ(K,ipass,1,0,3,X)= HPINV(K)*ETIJ(K,ipass,0,0,2,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,0,2,X)&
                      & + 2*ETIJ(K,ipass,2,0,2,X) - 2*ETIJ(K,ipass,1,0,1,X)/TWOB(K)
                 ETIJ(K,ipass,0,0,4,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,0,3,X)+ ETIJ(K,ipass, 1,0,3,X)&
                      & -3*ETIJ(K,ipass,0,0,2,X)/TWOB(K)
                 ETIJ(K,ipass,3,0,4,X) = HPINV(K)*ETIJ(K,ipass,2,0,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,0,3,X)
                 ETIJ(K,ipass,4,0,4,X) = HPINV(K)*ETIJ(K,ipass,3,0,3,X)
                 ETIJ(K,ipass,1,0,4,X)= HPINV(K)*ETIJ(K,ipass,0,0,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,0,3,X)&
                      & + 2*ETIJ(K,ipass,2,0,3,X) - 3*ETIJ(K,ipass,1,0,2,X)/TWOB(K)
                 ETIJ(K,ipass,2,0,4,X)= HPINV(K)*ETIJ(K,ipass,1,0,3,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,0,3,X)&
                      & + 3*ETIJ(K,ipass,3,0,3,X) - 3*ETIJ(K,ipass,2,0,2,X)/TWOB(K)
                 ETIJ(K,ipass,0,0,5,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,0,4,X)+ ETIJ(K,ipass, 1,0,4,X)&
                      & -4*ETIJ(K,ipass,0,0,3,X)/TWOB(K)
                 ETIJ(K,ipass,4,0,5,X) = HPINV(K)*ETIJ(K,ipass,3,0,4,X)+ PB(K,ipass,X)*ETIJ(K,ipass,4,0,4,X)
                 ETIJ(K,ipass,5,0,5,X) = HPINV(K)*ETIJ(K,ipass,4,0,4,X)
                 ETIJ(K,ipass,1,0,5,X)= HPINV(K)*ETIJ(K,ipass,0,0,4,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,0,4,X)&
                      & + 2*ETIJ(K,ipass,2,0,4,X) - 4*ETIJ(K,ipass,1,0,3,X)/TWOB(K)
                 ETIJ(K,ipass,2,0,5,X)= HPINV(K)*ETIJ(K,ipass,1,0,4,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,0,4,X)&
                      & + 3*ETIJ(K,ipass,3,0,4,X) - 4*ETIJ(K,ipass,2,0,3,X)/TWOB(K)
                 ETIJ(K,ipass,3,0,5,X)= HPINV(K)*ETIJ(K,ipass,2,0,4,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,0,4,X)&
                      & + 4*ETIJ(K,ipass,4,0,4,X) - 4*ETIJ(K,ipass,3,0,3,X)/TWOB(K)
                 ETIJ(K,ipass,0,0,6,X) = PB(K,ipass,X)*ETIJ(K,ipass, 0,0,5,X)+ ETIJ(K,ipass, 1,0,5,X)&
                      & -5*ETIJ(K,ipass,0,0,4,X)/TWOB(K)
                 ETIJ(K,ipass,5,0,6,X) = HPINV(K)*ETIJ(K,ipass,4,0,5,X)+ PB(K,ipass,X)*ETIJ(K,ipass,5,0,5,X)
                 ETIJ(K,ipass,6,0,6,X) = HPINV(K)*ETIJ(K,ipass,5,0,5,X)
                 ETIJ(K,ipass,1,0,6,X)= HPINV(K)*ETIJ(K,ipass,0,0,5,X)+ PB(K,ipass,X)*ETIJ(K,ipass,1,0,5,X)&
                      & + 2*ETIJ(K,ipass,2,0,5,X) - 5*ETIJ(K,ipass,1,0,4,X)/TWOB(K)
                 ETIJ(K,ipass,2,0,6,X)= HPINV(K)*ETIJ(K,ipass,1,0,5,X)+ PB(K,ipass,X)*ETIJ(K,ipass,2,0,5,X)&
                      & + 3*ETIJ(K,ipass,3,0,5,X) - 5*ETIJ(K,ipass,2,0,4,X)/TWOB(K)
                 ETIJ(K,ipass,3,0,6,X)= HPINV(K)*ETIJ(K,ipass,2,0,5,X)+ PB(K,ipass,X)*ETIJ(K,ipass,3,0,5,X)&
                      & + 4*ETIJ(K,ipass,4,0,5,X) - 5*ETIJ(K,ipass,3,0,4,X)/TWOB(K)
                 ETIJ(K,ipass,4,0,6,X)= HPINV(K)*ETIJ(K,ipass,3,0,5,X)+ PB(K,ipass,X)*ETIJ(K,ipass,4,0,5,X)&
                      & + 5*ETIJ(K,ipass,5,0,5,X) - 5*ETIJ(K,ipass,4,0,4,X)/TWOB(K)
              ENDDO
           ENDDO
        ENDDO
        !!$OMP END DO
     ENDIF
  ELSE
     CALL ichorquit('EXTEND ICHOR_HERM_ECOEFFS (MAIN2 IN etij)',-1)
  ENDIF
!$OMP END MASTER
!$OMP BARRIER
END SUBROUTINE ICHOR_HERM_ECOEFFS

subroutine printEcoeff(EcoeffN,nTUVQ,ijkQcart,nPrimP,nPass,lupri)
implicit none
integer     :: nTUVQ,ijkQcart,nPrimP,nPass,lupri
real(realk) :: EcoeffN(nPrimP,nPass,nTUVQ,ijkQcart)
!
integer :: i,j,k,ip
k=0
WRITE(LUPRI,'(A)') 'EcoeffN'
WRITE(lupri,'(A,I7)')'nTUVQ:    ',nTUVQ
WRITE(lupri,'(A,I7)')'ijkQcart: ',ijkQcart
WRITE(lupri,'(A,I7)')'nPrimP:   ',nPrimP
WRITE(lupri,'(A,I7)')'nPass:    ',nPass
WRITE(LUPRI,'(1X,A6,1X,A6,1X,A5,12X,A14 )') 'prim  ','tuv   ','iPass ','ijk=1,ijkQcart'
DO ip=1,nPass
 DO i=1,nPrimP
   DO j=1,nTUVQ
      WRITE(LUPRI,'(1X,I4,2X,I4,1X,I3,5ES18.9/,(15X,5ES18.9))') i,j,ip,(EcoeffN(i,ip,j,k),k=1,ijkQcart)
   ENDDO
 ENDDO
ENDDO

end subroutine printEcoeff

end MODULE IchorEriCoulombintegralCPUMcMGeneralEcoeffMod
