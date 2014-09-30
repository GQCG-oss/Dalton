!> @file
!> Contains the general McM driver 

!> \brief General McMurchie-Davidson Integral scheme
!> \author T. Kjaergaard
!> \date 2013 
MODULE IchorEriCoulombintegralCPUMcMGeneralEcoeffMod
use IchorPrecisionModule
use IchorEriCoulombintegralCPUMcMGeneralWTUVMod
private 
public :: Ichorbuild_Ecoeff_LHS,Ichorbuild_Ecoeff_RHS, printEcoeff

CONTAINS
Subroutine Ichorbuild_Ecoeff_RHS(nPrimP,nPrimA,nPrimB,maxAngP,maxAng1,maxAng2,nTUV,&
     & ijk,Aexp,Bexp,EcoeffN,Pdistance12,PreExpFac,intprint,lupri,ETIJ)
  implicit none
  integer,intent(in)        :: nPrimP,nPrimA,nPrimB,lupri,nTUV,ijk
  integer,intent(in)        :: maxAngP,maxAng1,maxAng2,intprint
  real(realk),intent(in)    :: Aexp(nPrimA),Bexp(nPrimB),Pdistance12(3)
  real(realk),intent(in)    :: PreExpFac(nPrimP)
  real(realk),intent(inout) :: Ecoeffn(nPrimP,nTUV,ijk)
  real(realk),intent(inout) :: ETIJ(nPrimP,0:maxAngP,0:maxAng1,0:maxAng2,3)
  !
  real(realk) :: PINV(nPrimP),PA(nPrimP),PB(nPrimP)
  real(realk) :: APINV(nPrimP),BPINV(nPrimP)
  real(realk) :: HPINV(nPrimP),TWOA(nPrimP)
  real(realk) :: TWOB(nPrimP)
  real(realk),parameter :: DHALF=0.5E0_realk,D2=2.0E0_realk
  !
  integer :: iPrimP
  integer,parameter :: nPasses = 1
  call build_auxiliary(nPrimA,nPrimB,Aexp,Bexp,PINV,APINV,BPINV,HPINV,TWOA,TWOB)
  !X DIST
  DO iPrimP=1,nPrimP
     PA(iPrimP) = -BPINV(iPrimP)*Pdistance12(1)
     PB(iPrimP) =  APINV(iPrimP)*Pdistance12(1)
  ENDDO  
  IF(maxAng1+maxAng2.GT.-1 )THEN !6
     CALL ICHOR_HERM_ECOEFFS_GENERAL_RHS(ETIJ,nPrimP,maxAng1,maxAng2,&
          & PA,PB,HPINV,TWOA,TWOB,1,LUPRI)  
  ELSE !special less than 6
     CALL ICHOR_HERM_ECOEFFS(ETIJ,nPrimP,nPasses,maxAng1,maxAng2,&
          & PA,PB,HPINV,TWOA,TWOB,1,LUPRI)  
  ENDIF
  IF(intprint.GT.30)THEN
     call PRINT_ETIJ(ETIJ,nPrimP,nPasses,maxAng1,maxAng2,1,LUPRI)
  ENDIF
  !Y DIST
  DO iPrimP=1,nPrimP
     PA(iPrimP) = -BPINV(iPrimP)*Pdistance12(2)
     PB(iPrimP) =  APINV(iPrimP)*Pdistance12(2)
  ENDDO  
  IF(maxAng1+maxAng2.GT.-1 )THEN !6
     CALL ICHOR_HERM_ECOEFFS_GENERAL_RHS(ETIJ,nPrimP,maxAng1,maxAng2,&
          & PA,PB,HPINV,TWOA,TWOB,2,LUPRI)
  ELSE !special less than 6
     CALL ICHOR_HERM_ECOEFFS(ETIJ,nPrimP,nPasses,maxAng1,maxAng2,&
          & PA,PB,HPINV,TWOA,TWOB,2,LUPRI)
  ENDIF
  IF(intprint.GT.30)THEN
     call PRINT_ETIJ(ETIJ,nPrimP,nPasses,maxAng1,maxAng2,2,LUPRI)
  ENDIF
  !Z DIST
  DO iPrimP=1,nPrimP
     PA(iPrimP) = -BPINV(iPrimP)*Pdistance12(3)
     PB(iPrimP) =  APINV(iPrimP)*Pdistance12(3)
  ENDDO  
  IF(maxAng1+maxAng2.GT.-1 )THEN !6
     CALL ICHOR_HERM_ECOEFFS_GENERAL_RHS(ETIJ,nPrimP,maxAng1,maxAng2,&
          & PA,PB,HPINV,TWOA,TWOB,3,LUPRI)
  ELSE !NORMALY call special for less than 6
     CALL ICHOR_HERM_ECOEFFS(ETIJ,nPrimP,nPasses,maxAng1,maxAng2,&
          & PA,PB,HPINV,TWOA,TWOB,3,LUPRI)
  ENDIF
  IF(intprint.GT.30)THEN
     call PRINT_ETIJ(ETIJ,nPrimP,nPasses,maxAng1,maxAng2,3,LUPRI)
  ENDIF
  IF(maxAngP.EQ.0)THEN
     IF(maxAng1.EQ.0)THEN
        call ICHOR_Ecoeffn_maxAngP0_maxAngA0_RHS(nPrimP,Ecoeffn,ETIJ,PreExpFac)
     ENDIF
  ELSE
     call ICHOR_Ecoeffn_general_RHS(nPrimP,nTUV,ijk,maxAngP,maxAng1,maxAng2,Ecoeffn,ETIJ,PreExpFac)
  ENDIF
end Subroutine Ichorbuild_Ecoeff_RHS

Subroutine Ichorbuild_Ecoeff_LHS(nPrimP,nPrimA,nPrimB,maxAngP,maxAng1,maxAng2,nTUV,&
     & ijk,Aexp,Bexp,EcoeffN,Pdistance12,PreExpFac,nPasses,nAtomsA,nAtomsB,&
     & IatomApass,IatomBpass,MaxPasses,intprint,lupri,ETIJ)
  implicit none
  integer,intent(in)        :: nPrimP,nPrimA,nPrimB,lupri,nTUV,ijk,MaxPasses
  integer,intent(in)        :: maxAngP,maxAng1,maxAng2,nPasses,intprint,nAtomsA,nAtomsB
  integer,intent(in)        :: IatomBpass(MaxPasses),IatomApass(MaxPasses)
  real(realk),intent(in)    :: Aexp(nPrimA),Bexp(nPrimB),Pdistance12(3,nAtomsA,nAtomsB)
  real(realk),intent(in)    :: PreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(inout) :: Ecoeffn(nPrimP,nPasses,nTUV,ijk)
  real(realk),intent(inout) :: ETIJ(nPrimP,nPasses,0:maxAngP,0:maxAng1,0:maxAng2,3)
  !
  real(realk) :: PINV(nPrimP),PA(nPrimP,nPasses),PB(nPrimP,nPasses)
  real(realk) :: APINV(nPrimP),BPINV(nPrimP)
  real(realk) :: HPINV(nPrimP),TWOA(nPrimP)
  real(realk) :: TWOB(nPrimP),dist
  real(realk),parameter :: DHALF=0.5E0_realk,D2=2.0E0_realk
  !
  integer :: K,iPass,iPrimP,iAtomA,iAtomB
  call build_auxiliary(nPrimA,nPrimB,Aexp,Bexp,PINV,APINV,BPINV,HPINV,TWOA,TWOB)
  !X DIST
  do iPass=1,nPasses
     iAtomA = iAtomApass(iPass)
     iAtomB = iAtomBpass(iPass)
     dist=-Pdistance12(1,iAtomA,iAtomB)
     DO iPrimP=1,nPrimP
        PA(iPrimP,iPass) = BPINV(iPrimP)*dist
     ENDDO
     dist=Pdistance12(1,iAtomA,iAtomB)
     DO iPrimP=1,nPrimP
        PB(iPrimP,iPass) =  APINV(iPrimP)*dist
     ENDDO
  enddo
  IF(maxAng1+maxAng2.GT.-1 )THEN !6
     CALL ICHOR_HERM_ECOEFFS_GENERAL_LHS(ETIJ,nPrimP,nPasses,maxAng1,maxAng2,&
          & PA,PB,HPINV,TWOA,TWOB,1,LUPRI)  
  ELSE !special less than 6
     CALL ICHOR_HERM_ECOEFFS(ETIJ,nPrimP,nPasses,maxAng1,maxAng2,&
          & PA,PB,HPINV,TWOA,TWOB,1,LUPRI)  
  ENDIF
  IF(intprint.GT.30)THEN
     call PRINT_ETIJ(ETIJ,nPrimP,nPasses,maxAng1,maxAng2,1,LUPRI)
  ENDIF
  !Y DIST
  do iPass=1,nPasses
     iAtomA = iAtomApass(iPass)
     iAtomB = iAtomBpass(iPass)
     dist=-Pdistance12(2,iAtomA,iAtomB)
     DO iPrimP=1,nPrimP
        PA(iPrimP,iPass) = BPINV(iPrimP)*dist
     ENDDO
     dist=Pdistance12(2,iAtomA,iAtomB)
     DO iPrimP=1,nPrimP
        PB(iPrimP,iPass) =  APINV(iPrimP)*dist
     ENDDO
  enddo
  IF(maxAng1+maxAng2.GT.-1 )THEN !6
     CALL ICHOR_HERM_ECOEFFS_GENERAL_LHS(ETIJ,nPrimP,nPasses,maxAng1,maxAng2,&
          & PA,PB,HPINV,TWOA,TWOB,2,LUPRI)
  ELSE !special less than 6
     CALL ICHOR_HERM_ECOEFFS(ETIJ,nPrimP,nPasses,maxAng1,maxAng2,&
          & PA,PB,HPINV,TWOA,TWOB,2,LUPRI)
  ENDIF
  IF(intprint.GT.30)THEN
     call PRINT_ETIJ(ETIJ,nPrimP,nPasses,maxAng1,maxAng2,2,LUPRI)
  ENDIF
  !Z DIST
  do iPass=1,nPasses
     iAtomA = iAtomApass(iPass)
     iAtomB = iAtomBpass(iPass)
     dist=-Pdistance12(3,iAtomA,iAtomB)
     DO iPrimP=1,nPrimP
        PA(iPrimP,iPass) = BPINV(iPrimP)*dist
     ENDDO
     dist=Pdistance12(3,iAtomA,iAtomB)
     DO iPrimP=1,nPrimP
        PB(iPrimP,iPass) =  APINV(iPrimP)*dist
     ENDDO
  enddo
  IF(maxAng1+maxAng2.GT.-1 )THEN !6
     CALL ICHOR_HERM_ECOEFFS_GENERAL_LHS(ETIJ,nPrimP,nPasses,maxAng1,maxAng2,&
          & PA,PB,HPINV,TWOA,TWOB,3,LUPRI)
  ELSE !NORMALY call special for less than 6
     CALL ICHOR_HERM_ECOEFFS(ETIJ,nPrimP,nPasses,maxAng1,maxAng2,&
          & PA,PB,HPINV,TWOA,TWOB,3,LUPRI)
  ENDIF
  IF(intprint.GT.30)THEN
     call PRINT_ETIJ(ETIJ,nPrimP,nPasses,maxAng1,maxAng2,3,LUPRI)
  ENDIF
  IF(maxAngP.EQ.0)THEN
     IF(maxAng1.EQ.0)THEN
        call ICHOR_Ecoeffn_maxAngP0_maxAngA0_LHS(nPrimP,nPasses,Ecoeffn,ETIJ,PreExpFac,&
             & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
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

subroutine build_auxiliary(nPrimA,nPrimB,Aexp,Bexp,PINV,APINV,BPINV,HPINV,TWOA,TWOB)
    implicit none
    integer,intent(in) :: nPrimA,nPrimB
    real(realk),intent(in) :: Aexp(nPrimA),Bexp(nPrimB)
    real(realk),intent(inout) :: APINV(nPrimA,nPrimB),BPINV(nPrimA,nPrimB)
    real(realk),intent(inout) :: HPINV(nPrimA,nPrimB),TWOA(nPrimA,nPrimB)
    real(realk),intent(inout) :: TWOB(nPrimA,nPrimB),PINV(nPrimA,nPrimB)
    real(realk),parameter :: DHALF=0.5E0_realk,D2=2.0E0_realk,D1=1.0E0_realk
    !
    real(realk) :: bexpo
    integer :: J,I
    DO J=1,nPrimB
       bexpo = Bexp(J)
       DO I=1,nPrimA
          PINV(I,J)  = D1/(Aexp(I)+bexpo)   ! 1/p
          HPINV(I,J) = DHALF   * PINV(I,J)  ! 1/2p
          BPINV(I,J) = bexpo   * PINV(I,J)  ! b/p
          APINV(I,J) = Aexp(I) * PINV(I,J)  ! a/p
       ENDDO
       bexpo = D2*bexpo
       DO I=1,nPrimA
          TWOB(I,J)  = bexpo                ! 2b
       ENDDO
    ENDDO
    DO J=1,nPrimB
       DO I=1,nPrimA
          TWOA(I,J)  = D2*Aexp(I)           ! 2a
       ENDDO
    ENDDO
end subroutine build_auxiliary

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
     & HPINV,TWOA,TWOB,X,LUPRI)
implicit none
INTEGER,intent(in)        :: nPrimP,nPass,MAXI,MAXJ,LUPRI,X
REAL(REALK),intent(in)    :: PA(nPrimP,nPass), PB(nPrimP,nPass)
REAL(REALK),intent(inout) :: ETIJ(nPrimP,nPass,0:MAXI+MAXJ,0:MAXI,0:MAXJ,3)
REAL(REALK),intent(in)    :: HPINV(nPrimP)
REAL(REALK)               :: TWOA(nPrimP),TWOB(nPrimP)
!
REAL(REALK),PARAMETER     :: D1 = 1.0E0_realk, D2 = 2.0E0_realk
INTEGER                   :: MAXIJ,iPass,K,IJ,I,J,TJM,TIM,T,T1
!
! Run over I (J = 0)                                     
! ==================
!
DO I = 0, MAXI
   TIM = I-1
   IF (I .LE. 2) THEN
      !           E(0,0,0)                                              
      IF (I .EQ. 0) THEN
         DO Ipass = 1,nPass
            DO K = 1, nPrimP
               ETIJ(K,ipass,0,0,0,X) = D1
            ENDDO
         ENDDO
      ELSE IF (I .EQ. 1) THEN
         !           E(T,1,0)                                            
         DO Ipass = 1,nPass
            DO K = 1, nPrimP
               ETIJ(K,ipass,0,1,0,X) = PA(K,ipass)
               ETIJ(K,ipass,1,1,0,X) = HPINV(K)
            ENDDO
         ENDDO
      ELSE IF (I .EQ. 2) THEN
         !           E(T,2,0)                                              
         DO Ipass = 1,nPass
            DO K = 1, nPrimP
               ETIJ(K,ipass,0,2,0,X) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
               ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass)*HPINV(K)
               ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
            ENDDO
         ENDDO
      ENDIF
   ELSE
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,I,0,X) = PA(K,ipass)*ETIJ(K,ipass,0,I-1,0,X) + ETIJ(K,ipass,1,I-1,0,X) - TIM*ETIJ(K,ipass,0,I-2,0,X)/TWOA(K)
            ETIJ(K,ipass,I-1,I,0,X) = HPINV(K)*ETIJ(K,ipass,I-2,I-1,0,X) + PA(K,ipass)*ETIJ(K,ipass,I-1,I-1,0,X)
            ETIJ(K,ipass,I,I,0,X) = HPINV(K)*ETIJ(K,ipass,I-1,I-1,0,X)
         END DO
      ENDDO
      DO T = 1, I - 2                                    
         T1 = T + 1
         DO Ipass = 1,nPass
            DO K = 1, nPrimP
               ETIJ(K,ipass,T,I,0,X) = HPINV(K)*ETIJ(K,ipass,T-1,I-1,0,X)+ PA(K,ipass)*ETIJ(K,ipass,  T,I-1,0,X)&
                    & + T1*ETIJ(K,ipass,T+1,I-1,0,X) - TIM*ETIJ(K,ipass, T,I-2,0,X)/TWOA(K)
            ENDDO
         ENDDO
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
            DO Ipass = 1,nPass
               DO K = 1, nPrimP
                  ETIJ(K,ipass,0,0,1,X) = PB(K,ipass)    
                  ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
               ENDDO
            ENDDO
         ELSE IF (IJ .EQ. 2) THEN     
            !                 E(0,2)                        
            IF (I .EQ. 0) THEN
               DO Ipass = 1,nPass
                  DO K = 1, nPrimP
                     ETIJ(K,ipass,0,0,2,X) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
                     ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass)*HPINV(K)
                     ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
                  ENDDO
               ENDDO
            ELSE
               !                    E(1,1)                                    
               DO Ipass = 1,nPass
                  DO K = 1, nPrimP
                     ETIJ(K,ipass,0,1,1,X) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
                     ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
                     ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
                  ENDDO
               ENDDO
            END IF
         ENDIF
      ELSE
         !             E(I,J)                                            
         IF(J.LT. 2)THEN
            DO Ipass = 1,nPass
               DO K = 1, nPrimP
                  ETIJ(K,ipass,0,I,J,X) = PB(K,ipass)*ETIJ(K,ipass,0,I,J-1,X)+ ETIJ(K,ipass,1,I,J-1,X)
                  ETIJ(K,ipass,IJ-1,I,J,X) = HPINV(K)*ETIJ(K,ipass,IJ-2,I,J-1,X)+ PB(K,ipass)*ETIJ(K,ipass,IJ-1,I,J-1,X)
                  ETIJ(K,ipass,IJ,I,J,X) = HPINV(K)*ETIJ(K,ipass,IJ-1,I,J-1,X)
               ENDDO
            ENDDO
            DO T = 1, IJ - 2
               T1 = T + 1
               DO Ipass = 1,nPass
                  DO K = 1, nPrimP
                     ETIJ(K,ipass,T,I,J,X)= HPINV(K)*ETIJ(K,ipass,T-1,I,J-1,X)+ PB(K,ipass)*ETIJ(K,ipass,  T,I,J-1,X)&
                          & + T1*ETIJ(K,ipass,T+1,I,J-1,X)! - TJM*ETIJ(K,ipass,  T,I,J-2)/TWOB(K)
                  ENDDO
               ENDDO
            END DO
         ELSE
            DO Ipass = 1,nPass
               DO K = 1, nPrimP
                  ETIJ(K,ipass,   0,I,J,X) = PB(K,ipass)*ETIJ(K,ipass,   0,I,J-1,X)+ ETIJ(K,ipass,   1,I,J-1,X)&
                       & -TJM*ETIJ(K,ipass,  0,I,J-2,X)/TWOB(K)
                  ETIJ(K,ipass,IJ-1,I,J,X) = HPINV(K)*ETIJ(K,ipass,IJ-2,I,J-1,X)+ PB(K,ipass)*ETIJ(K,ipass,IJ-1,I,J-1,X)
                  ETIJ(K,ipass,  IJ,I,J,X) = HPINV(K)*ETIJ(K,ipass,IJ-1,I,J-1,X)
               ENDDO
            ENDDO
            DO T = 1, IJ - 2
               T1 = T + 1
               DO Ipass = 1,nPass
                  DO K = 1, nPrimP
                     ETIJ(K,ipass,T,I,J,X)= HPINV(K)*ETIJ(K,ipass,T-1,I,J-1,X)+ PB(K,ipass)*ETIJ(K,ipass,  T,I,J-1,X)&
                          & + T1*ETIJ(K,ipass,T+1,I,J-1,X) - TJM*ETIJ(K,ipass,  T,I,J-2,X)/TWOB(K)
                  ENDDO
               ENDDO
            END DO
         ENDIF
      END IF
   END DO
END DO
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
     & HPINV,TWOA,TWOB,X,LUPRI)
implicit none
INTEGER,intent(in)        :: nPrimP,MAXI,MAXJ,LUPRI,X
REAL(REALK),intent(in)    :: PA(nPrimP), PB(nPrimP)
REAL(REALK),intent(inout) :: ETIJ(nPrimP,0:MAXI+MAXJ,0:MAXI,0:MAXJ,3)
REAL(REALK),intent(in)    :: HPINV(nPrimP)
REAL(REALK)               :: TWOA(nPrimP),TWOB(nPrimP)
!
REAL(REALK),PARAMETER     :: D1 = 1.0E0_realk, D2 = 2.0E0_realk
INTEGER                   :: MAXIJ,K,IJ,I,J,TJM,TIM,T,T1
!
! Run over I (J = 0)                                     
! ==================
!
DO I = 0, MAXI
   TIM = I-1
   IF (I .LE. 2) THEN
      !           E(0,0,0)                                              
      IF (I .EQ. 0) THEN
         DO K = 1, nPrimP
            ETIJ(K,0,0,0,X) = D1
         ENDDO
      ELSE IF (I .EQ. 1) THEN
         !           E(T,1,0)                                            
         DO K = 1, nPrimP
            ETIJ(K,0,1,0,X) = PA(K)
            ETIJ(K,1,1,0,X) = HPINV(K)
         ENDDO
      ELSE IF (I .EQ. 2) THEN
         !           E(T,2,0)                                              
         DO K = 1, nPrimP
            ETIJ(K,0,2,0,X) = PA(K)*PA(K) + HPINV(K) - D1/TWOA(K)
            ETIJ(K,1,2,0,X) = D2*PA(K)*HPINV(K)
            ETIJ(K,2,2,0,X) = HPINV(K)*HPINV(K)
         ENDDO
      ENDIF
   ELSE
      DO K = 1, nPrimP
         ETIJ(K,0,I,0,X) = PA(K)*ETIJ(K,0,I-1,0,X) + ETIJ(K,1,I-1,0,X) - TIM*ETIJ(K,0,I-2,0,X)/TWOA(K)
         ETIJ(K,I-1,I,0,X) = HPINV(K)*ETIJ(K,I-2,I-1,0,X) + PA(K)*ETIJ(K,I-1,I-1,0,X)
         ETIJ(K,I,I,0,X) = HPINV(K)*ETIJ(K,I-1,I-1,0,X)
      END DO
      DO T = 1, I - 2                                    
         T1 = T + 1
         DO K = 1, nPrimP
            ETIJ(K,T,I,0,X) = HPINV(K)*ETIJ(K,T-1,I-1,0,X)+ PA(K)*ETIJ(K,  T,I-1,0,X)&
                 & + T1*ETIJ(K,T+1,I-1,0,X) - TIM*ETIJ(K, T,I-2,0,X)/TWOA(K)
         ENDDO
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
            DO K = 1, nPrimP
               ETIJ(K,0,0,1,X) = PB(K)    
               ETIJ(K,1,0,1,X) = HPINV(K) 
            ENDDO
         ELSE IF (IJ .EQ. 2) THEN     
            !                 E(0,2)                        
            IF (I .EQ. 0) THEN
               DO K = 1, nPrimP
                  ETIJ(K,0,0,2,X) = PB(K)*PB(K) + HPINV(K) - D1/TWOB(K)
                  ETIJ(K,1,0,2,X) = D2*PB(K)*HPINV(K)
                  ETIJ(K,2,0,2,X) = HPINV(K)*HPINV(K)
               ENDDO
            ELSE
               !                    E(1,1)                                    
               DO K = 1, nPrimP
                  ETIJ(K,0,1,1,X) = PA(K)*PB(K) + HPINV(K)
                  ETIJ(K,1,1,1,X) = (PA(K) + PB(K))*HPINV(K)
                  ETIJ(K,2,1,1,X) = HPINV(K)*HPINV(K)
               ENDDO
            END IF
         ENDIF
      ELSE
         !             E(I,J)                                            
         IF(J.LT. 2)THEN
            DO K = 1, nPrimP
               ETIJ(K,0,I,J,X) = PB(K)*ETIJ(K,0,I,J-1,X)+ ETIJ(K,1,I,J-1,X)
               ETIJ(K,IJ-1,I,J,X) = HPINV(K)*ETIJ(K,IJ-2,I,J-1,X)+ PB(K)*ETIJ(K,IJ-1,I,J-1,X)
               ETIJ(K,IJ,I,J,X) = HPINV(K)*ETIJ(K,IJ-1,I,J-1,X)
            ENDDO
            DO T = 1, IJ - 2
               T1 = T + 1
               DO K = 1, nPrimP
                  ETIJ(K,T,I,J,X)= HPINV(K)*ETIJ(K,T-1,I,J-1,X)+ PB(K)*ETIJ(K,  T,I,J-1,X)&
                       & + T1*ETIJ(K,T+1,I,J-1,X)! - TJM*ETIJ(K,  T,I,J-2)/TWOB(K)
               ENDDO
            END DO
         ELSE
            DO K = 1, nPrimP
               ETIJ(K,   0,I,J,X) = PB(K)*ETIJ(K,   0,I,J-1,X)+ ETIJ(K,   1,I,J-1,X)&
                    & -TJM*ETIJ(K,  0,I,J-2,X)/TWOB(K)
               ETIJ(K,IJ-1,I,J,X) = HPINV(K)*ETIJ(K,IJ-2,I,J-1,X)+ PB(K)*ETIJ(K,IJ-1,I,J-1,X)
               ETIJ(K,  IJ,I,J,X) = HPINV(K)*ETIJ(K,IJ-1,I,J-1,X)
            ENDDO
            DO T = 1, IJ - 2
               T1 = T + 1
               DO K = 1, nPrimP
                  ETIJ(K,T,I,J,X)= HPINV(K)*ETIJ(K,T-1,I,J-1,X)+ PB(K)*ETIJ(K,  T,I,J-1,X)&
                       & + T1*ETIJ(K,T+1,I,J-1,X) - TJM*ETIJ(K,  T,I,J-2,X)/TWOB(K)
               ENDDO
            END DO
         ENDIF
      END IF
   END DO
END DO
END SUBROUTINE ICHOR_HERM_ECOEFFS_GENERAL_RHS

SUBROUTINE ICHOR_HERM_ECOEFFS(ETIJ,nPrimP,nPass,MAXI,MAXJ,PA,PB,&
     & HPINV,TWOA,TWOB,X,LUPRI)
implicit none
INTEGER,intent(in)        :: nPrimP,nPass,MAXI,MAXJ,LUPRI,X
REAL(REALK),intent(in)    :: PA(nPrimP,nPass), PB(nPrimP,nPass)
REAL(REALK),intent(inout) :: ETIJ(nPrimP,nPass,0:MAXI+MAXJ,0:MAXI,0:MAXJ,3)
REAL(REALK),intent(in)    :: HPINV(nPrimP)
REAL(REALK)               :: TWOA(nPrimP),TWOB(nPrimP)
!
REAL(REALK),PARAMETER     :: D1 = 1.0E0_realk, D2 = 2.0E0_realk
INTEGER                   :: MAXIJ,iPass,K
MAXIJ = MAXI + MAXJ
IF(maxIJ.EQ.0 )THEN
   IF(maxI.EQ.0 )THEN
      !maxI.EQ.0.AND.maxJ.EQ.0
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
         ENDDO
      ENDDO
   ENDIF
ELSEIF(maxIJ.EQ.1 )THEN
   IF(maxI.EQ.1 )THEN
      !maxI.EQ.1.AND.maxJ.EQ.0
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,1,0,X) = PA(K,ipass)
            ETIJ(K,ipass,1,1,0,X) = HPINV(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.0)THEN
      !maxI.EQ.0.AND.maxJ.EQ.1
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,0,1,X) = PB(K,ipass)    
            ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
         ENDDO
      ENDDO
   ENDIF
ELSEIF(maxIJ.EQ.2 )THEN
   IF(maxI.EQ.2 )THEN
      !maxI.EQ.2.AND.maxJ.EQ.0
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,1,0,X) = PA(K,ipass)
            ETIJ(K,ipass,1,1,0,X) = HPINV(K)
            ETIJ(K,ipass,0,2,0,X) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.1)THEN
      !maxI.EQ.1.AND.maxJ.EQ.1
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,0,1,X) = PB(K,ipass)    
            ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
            ETIJ(K,ipass,0,1,0,X) = PA(K,ipass)
            ETIJ(K,ipass,1,1,0,X) = HPINV(K)
            ETIJ(K,ipass,0,1,1,X) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.0)THEN
      !maxI.EQ.0.AND.maxJ.EQ.2
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,0,1,X) = PB(K,ipass)    
            ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
            ETIJ(K,ipass,0,0,2,X) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
         ENDDO
      ENDDO
   ENDIF
ELSEIF(maxIJ.EQ.3 )THEN
   IF(maxI.EQ.3 )THEN
      !maxI.EQ.3.AND.maxJ.EQ.0
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,1,0,X) = PA(K,ipass)
            ETIJ(K,ipass,1,1,0,X) = HPINV(K)
            ETIJ(K,ipass,0,2,0,X) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,3,0,X) = PA(K,ipass)*ETIJ(K,ipass,0,2,0,X)+ ETIJ(K,ipass,1,2,0,X)&
                 & - 2*ETIJ(K,ipass, 0,1,0,X)/TWOA(K)
            ETIJ(K,ipass,2,3,0,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PA(K,ipass)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,3,3,0,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,1,3,0,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PA(K,ipass)*ETIJ(K,ipass,1,2,0,X)&
                 & + 2*ETIJ(K,ipass,2,2,0,X) - 2*ETIJ(K,ipass, 1,1,0,X)/TWOA(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.2)THEN
      !maxI.EQ.2.AND.maxJ.EQ.1
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,0,1,X) = PB(K,ipass)    
            ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
            ETIJ(K,ipass,0,1,0,X) = PA(K,ipass)
            ETIJ(K,ipass,1,1,0,X) = HPINV(K)
            ETIJ(K,ipass,0,1,1,X) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,2,0,X) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,2,1,X) = PB(K,ipass)*ETIJ(K,ipass, 0,2,0,X)+ ETIJ(K,ipass, 1,2,0,X)
            ETIJ(K,ipass,2,2,1,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PB(K,ipass)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,3,2,1,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,1,2,1,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PB(K,ipass)*ETIJ(K,ipass,1,2,0,X)&
                 & + 2*ETIJ(K,ipass,2,2,0,X)! - 0*ETIJ(K,ipass,1,2,*,X)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.1)THEN
      !maxI.EQ.1.AND.maxJ.EQ.2
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,0,1,X) = PB(K,ipass)    
            ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
            ETIJ(K,ipass,0,0,2,X) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,1,0,X) = PA(K,ipass)
            ETIJ(K,ipass,1,1,0,X) = HPINV(K)
            ETIJ(K,ipass,0,1,1,X) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,1,2,X) = PB(K,ipass)*ETIJ(K,ipass, 0,1,1,X)+ ETIJ(K,ipass, 1,1,1,X)&
                 & -1*ETIJ(K,ipass,0,1,0,X)/TWOB(K)
            ETIJ(K,ipass,2,1,2,X) = HPINV(K)*ETIJ(K,ipass,1,1,1,X)+ PB(K,ipass)*ETIJ(K,ipass,2,1,1,X)
            ETIJ(K,ipass,3,1,2,X) = HPINV(K)*ETIJ(K,ipass,2,1,1,X)
            ETIJ(K,ipass,1,1,2,X) = HPINV(K)*ETIJ(K,ipass,0,1,1,X)+ PB(K,ipass)*ETIJ(K,ipass,1,1,1,X)&
                 & + 2*ETIJ(K,ipass,2,1,1,X) - 1*ETIJ(K,ipass,1,1,0,X)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.0)THEN
      !maxI.EQ.0.AND.maxJ.EQ.3
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,0,1,X) = PB(K,ipass)    
            ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
            ETIJ(K,ipass,0,0,2,X) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,0,3,X) = PB(K,ipass)*ETIJ(K,ipass, 0,0,2,X)+ ETIJ(K,ipass, 1,0,2,X)&
                 & -2*ETIJ(K,ipass,0,0,1,X)/TWOB(K)
            ETIJ(K,ipass,2,0,3,X) = HPINV(K)*ETIJ(K,ipass,1,0,2,X)+ PB(K,ipass)*ETIJ(K,ipass,2,0,2,X)
            ETIJ(K,ipass,3,0,3,X) = HPINV(K)*ETIJ(K,ipass,2,0,2,X)
            ETIJ(K,ipass,1,0,3,X) = HPINV(K)*ETIJ(K,ipass,0,0,2,X)+ PB(K,ipass)*ETIJ(K,ipass,1,0,2,X)&
                 & + 2*ETIJ(K,ipass,2,0,2,X) - 2*ETIJ(K,ipass,1,0,1,X)/TWOB(K)
         ENDDO
      ENDDO
   ENDIF
ELSEIF(maxIJ.EQ.4 )THEN
   IF(maxI.EQ.4 )THEN
      !maxI.EQ.4.AND.maxJ.EQ.0
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,1,0,X) = PA(K,ipass)
            ETIJ(K,ipass,1,1,0,X) = HPINV(K)
            ETIJ(K,ipass,0,2,0,X) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,3,0,X) = PA(K,ipass)*ETIJ(K,ipass,0,2,0,X)+ ETIJ(K,ipass,1,2,0,X)&
                 & - 2*ETIJ(K,ipass, 0,1,0,X)/TWOA(K)
            ETIJ(K,ipass,2,3,0,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PA(K,ipass)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,3,3,0,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,1,3,0,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PA(K,ipass)*ETIJ(K,ipass,1,2,0,X)&
                 & + 2*ETIJ(K,ipass,2,2,0,X) - 2*ETIJ(K,ipass, 1,1,0,X)/TWOA(K)
            ETIJ(K,ipass,0,4,0,X) = PA(K,ipass)*ETIJ(K,ipass,0,3,0,X)+ ETIJ(K,ipass,1,3,0,X)&
                 & - 3*ETIJ(K,ipass, 0,2,0,X)/TWOA(K)
            ETIJ(K,ipass,3,4,0,X) = HPINV(K)*ETIJ(K,ipass,2,3,0,X)+ PA(K,ipass)*ETIJ(K,ipass,3,3,0,X)
            ETIJ(K,ipass,4,4,0,X) = HPINV(K)*ETIJ(K,ipass,3,3,0,X)
            ETIJ(K,ipass,1,4,0,X) = HPINV(K)*ETIJ(K,ipass,0,3,0,X)+ PA(K,ipass)*ETIJ(K,ipass,1,3,0,X)&
                 & + 2*ETIJ(K,ipass,2,3,0,X) - 3*ETIJ(K,ipass, 1,2,0,X)/TWOA(K)
            ETIJ(K,ipass,2,4,0,X) = HPINV(K)*ETIJ(K,ipass,1,3,0,X)+ PA(K,ipass)*ETIJ(K,ipass,2,3,0,X)&
                 & + 3*ETIJ(K,ipass,3,3,0,X) - 3*ETIJ(K,ipass, 2,2,0,X)/TWOA(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.3)THEN
      !maxI.EQ.3.AND.maxJ.EQ.1
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,0,1,X) = PB(K,ipass)    
            ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
            ETIJ(K,ipass,0,1,0,X) = PA(K,ipass)
            ETIJ(K,ipass,1,1,0,X) = HPINV(K)
            ETIJ(K,ipass,0,1,1,X) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,2,0,X) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,2,1,X) = PB(K,ipass)*ETIJ(K,ipass, 0,2,0,X)+ ETIJ(K,ipass, 1,2,0,X)
            ETIJ(K,ipass,2,2,1,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PB(K,ipass)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,3,2,1,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,1,2,1,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PB(K,ipass)*ETIJ(K,ipass,1,2,0,X)&
                 & + 2*ETIJ(K,ipass,2,2,0,X)! - 0*ETIJ(K,ipass,1,2,*,X)/TWOB(K)
            ETIJ(K,ipass,0,3,0,X) = PA(K,ipass)*ETIJ(K,ipass,0,2,0,X)+ ETIJ(K,ipass,1,2,0,X)&
                 & - 2*ETIJ(K,ipass, 0,1,0,X)/TWOA(K)
            ETIJ(K,ipass,2,3,0,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PA(K,ipass)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,3,3,0,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,1,3,0,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PA(K,ipass)*ETIJ(K,ipass,1,2,0,X)&
                 & + 2*ETIJ(K,ipass,2,2,0,X) - 2*ETIJ(K,ipass, 1,1,0,X)/TWOA(K)
            ETIJ(K,ipass,0,3,1,X) = PB(K,ipass)*ETIJ(K,ipass, 0,3,0,X)+ ETIJ(K,ipass, 1,3,0,X)
            ETIJ(K,ipass,3,3,1,X) = HPINV(K)*ETIJ(K,ipass,2,3,0,X)+ PB(K,ipass)*ETIJ(K,ipass,3,3,0,X)
            ETIJ(K,ipass,4,3,1,X) = HPINV(K)*ETIJ(K,ipass,3,3,0,X)
            ETIJ(K,ipass,1,3,1,X) = HPINV(K)*ETIJ(K,ipass,0,3,0,X)+ PB(K,ipass)*ETIJ(K,ipass,1,3,0,X)&
                 & + 2*ETIJ(K,ipass,2,3,0,X)! - 0*ETIJ(K,ipass,1,3,*,X)/TWOB(K)
            ETIJ(K,ipass,2,3,1,X) = HPINV(K)*ETIJ(K,ipass,1,3,0,X)+ PB(K,ipass)*ETIJ(K,ipass,2,3,0,X)&
                 & + 3*ETIJ(K,ipass,3,3,0,X)! - 0*ETIJ(K,ipass,2,3,*,X)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.2)THEN
      !maxI.EQ.2.AND.maxJ.EQ.2
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,0,1,X) = PB(K,ipass)    
            ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
            ETIJ(K,ipass,0,0,2,X) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,1,0,X) = PA(K,ipass)
            ETIJ(K,ipass,1,1,0,X) = HPINV(K)
            ETIJ(K,ipass,0,1,1,X) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)

            ETIJ(K,ipass,0,1,2,X) = PB(K,ipass)*ETIJ(K,ipass,0,1,1,X)+ ETIJ(K,ipass,1,1,1,X)&
                                & -1*ETIJ(K,ipass,0,1,0,X)/TWOB(K)
            ETIJ(K,ipass,2,1,2,X) = HPINV(K)*ETIJ(K,ipass,1,1,1,X)+ PB(K,ipass)*ETIJ(K,ipass,2,1,1,X)
            ETIJ(K,ipass,3,1,2,X) = HPINV(K)*ETIJ(K,ipass,2,1,1,X)
            ETIJ(K,ipass,1,1,2,X) = HPINV(K)*ETIJ(K,ipass,0,1,1,X)+ PB(K,ipass)*ETIJ(K,ipass,1,1,1,X)&
                                & + 2*ETIJ(K,ipass,2,1,1,X) - 1*ETIJ(K,ipass,1,1,0,X)/TWOB(K)
            ETIJ(K,ipass,0,2,0,X) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,2,1,X) = PB(K,ipass)*ETIJ(K,ipass,0,2,0,X)+ ETIJ(K,ipass,1,2,0,X)
            ETIJ(K,ipass,2,2,1,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PB(K,ipass)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,3,2,1,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,1,2,1,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PB(K,ipass)*ETIJ(K,ipass,1,2,0,X)&
                                & + 2*ETIJ(K,ipass,2,2,0,X)! - 0*ETIJ(K,ipass,1,2,*,X)/TWOB(K)
            ETIJ(K,ipass,0,2,2,X) = PB(K,ipass)*ETIJ(K,ipass, 0,2,1,X)+ ETIJ(K,ipass, 1,2,1,X)&
                                & -1*ETIJ(K,ipass,0,2,0,X)/TWOB(K)
            ETIJ(K,ipass,3,2,2,X) = HPINV(K)*ETIJ(K,ipass,2,2,1,X)+ PB(K,ipass)*ETIJ(K,ipass,3,2,1,X)
            ETIJ(K,ipass,4,2,2,X) = HPINV(K)*ETIJ(K,ipass,3,2,1,X)
            ETIJ(K,ipass,1,2,2,X) = HPINV(K)*ETIJ(K,ipass,0,2,1,X)+ PB(K,ipass)*ETIJ(K,ipass,1,2,1,X)&
                                & + 2*ETIJ(K,ipass,2,2,1,X) - 1*ETIJ(K,ipass,1,2,0,X)/TWOB(K)
            ETIJ(K,ipass,2,2,2,X) = HPINV(K)*ETIJ(K,ipass,1,2,1,X)+ PB(K,ipass)*ETIJ(K,ipass,2,2,1,X)&
                                & + 3*ETIJ(K,ipass,3,2,1,X) - 1*ETIJ(K,ipass,2,2,0,X)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.1)THEN
      !maxI.EQ.1.AND.maxJ.EQ.3
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,0,1,X) = PB(K,ipass)    
            ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
            ETIJ(K,ipass,0,0,2,X) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,0,3,X) = PB(K,ipass)*ETIJ(K,ipass, 0,0,2,X)+ ETIJ(K,ipass, 1,0,2,X)&
                 & -2*ETIJ(K,ipass,0,0,1,X)/TWOB(K)
            ETIJ(K,ipass,2,0,3,X) = HPINV(K)*ETIJ(K,ipass,1,0,2,X)+ PB(K,ipass)*ETIJ(K,ipass,2,0,2,X)
            ETIJ(K,ipass,3,0,3,X) = HPINV(K)*ETIJ(K,ipass,2,0,2,X)
            ETIJ(K,ipass,1,0,3,X) = HPINV(K)*ETIJ(K,ipass,0,0,2,X)+ PB(K,ipass)*ETIJ(K,ipass,1,0,2,X)&
                 & + 2*ETIJ(K,ipass,2,0,2,X) - 2*ETIJ(K,ipass,1,0,1,X)/TWOB(K)
            ETIJ(K,ipass,0,1,0,X) = PA(K,ipass)
            ETIJ(K,ipass,1,1,0,X) = HPINV(K)
            ETIJ(K,ipass,0,1,1,X) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,1,2,X) = PB(K,ipass)*ETIJ(K,ipass, 0,1,1,X)+ ETIJ(K,ipass, 1,1,1,X)&
                 & -1*ETIJ(K,ipass,0,1,0,X)/TWOB(K)
            ETIJ(K,ipass,2,1,2,X) = HPINV(K)*ETIJ(K,ipass,1,1,1,X)+ PB(K,ipass)*ETIJ(K,ipass,2,1,1,X)
            ETIJ(K,ipass,3,1,2,X) = HPINV(K)*ETIJ(K,ipass,2,1,1,X)
            ETIJ(K,ipass,1,1,2,X) = HPINV(K)*ETIJ(K,ipass,0,1,1,X)+ PB(K,ipass)*ETIJ(K,ipass,1,1,1,X)&
                 & + 2*ETIJ(K,ipass,2,1,1,X) - 1*ETIJ(K,ipass,1,1,0,X)/TWOB(K)
            ETIJ(K,ipass,0,1,3,X) = PB(K,ipass)*ETIJ(K,ipass, 0,1,2,X)+ ETIJ(K,ipass, 1,1,2,X)&
                 & -2*ETIJ(K,ipass,0,1,1,X)/TWOB(K)
            ETIJ(K,ipass,3,1,3,X) = HPINV(K)*ETIJ(K,ipass,2,1,2,X)+ PB(K,ipass)*ETIJ(K,ipass,3,1,2,X)
            ETIJ(K,ipass,4,1,3,X) = HPINV(K)*ETIJ(K,ipass,3,1,2,X)
            ETIJ(K,ipass,1,1,3,X) = HPINV(K)*ETIJ(K,ipass,0,1,2,X)+ PB(K,ipass)*ETIJ(K,ipass,1,1,2,X)&
                 & + 2*ETIJ(K,ipass,2,1,2,X) - 2*ETIJ(K,ipass,1,1,1,X)/TWOB(K)
            ETIJ(K,ipass,2,1,3,X) = HPINV(K)*ETIJ(K,ipass,1,1,2,X)+ PB(K,ipass)*ETIJ(K,ipass,2,1,2,X)&
                 & + 3*ETIJ(K,ipass,3,1,2,X) - 2*ETIJ(K,ipass,2,1,1,X)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.0)THEN
      !maxI.EQ.0.AND.maxJ.EQ.4
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,0,1,X) = PB(K,ipass)    
            ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
            ETIJ(K,ipass,0,0,2,X) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,0,3,X) = PB(K,ipass)*ETIJ(K,ipass, 0,0,2,X)+ ETIJ(K,ipass, 1,0,2,X)&
                 & -2*ETIJ(K,ipass,0,0,1,X)/TWOB(K)
            ETIJ(K,ipass,2,0,3,X) = HPINV(K)*ETIJ(K,ipass,1,0,2,X)+ PB(K,ipass)*ETIJ(K,ipass,2,0,2,X)
            ETIJ(K,ipass,3,0,3,X) = HPINV(K)*ETIJ(K,ipass,2,0,2,X)
            ETIJ(K,ipass,1,0,3,X) = HPINV(K)*ETIJ(K,ipass,0,0,2,X)+ PB(K,ipass)*ETIJ(K,ipass,1,0,2,X)&
                 & + 2*ETIJ(K,ipass,2,0,2,X) - 2*ETIJ(K,ipass,1,0,1,X)/TWOB(K)
            ETIJ(K,ipass,0,0,4,X) = PB(K,ipass)*ETIJ(K,ipass, 0,0,3,X)+ ETIJ(K,ipass, 1,0,3,X)&
                 & -3*ETIJ(K,ipass,0,0,2,X)/TWOB(K)
            ETIJ(K,ipass,3,0,4,X) = HPINV(K)*ETIJ(K,ipass,2,0,3,X)+ PB(K,ipass)*ETIJ(K,ipass,3,0,3,X)
            ETIJ(K,ipass,4,0,4,X) = HPINV(K)*ETIJ(K,ipass,3,0,3,X)
            ETIJ(K,ipass,1,0,4,X) = HPINV(K)*ETIJ(K,ipass,0,0,3,X)+ PB(K,ipass)*ETIJ(K,ipass,1,0,3,X)&
                 & + 2*ETIJ(K,ipass,2,0,3,X) - 3*ETIJ(K,ipass,1,0,2,X)/TWOB(K)
            ETIJ(K,ipass,2,0,4,X) = HPINV(K)*ETIJ(K,ipass,1,0,3,X)+ PB(K,ipass)*ETIJ(K,ipass,2,0,3,X)&
                 & + 3*ETIJ(K,ipass,3,0,3,X) - 3*ETIJ(K,ipass,2,0,2,X)/TWOB(K)
         ENDDO
      ENDDO
   ENDIF
ELSEIF(maxIJ.EQ.5 )THEN
   IF(maxI.EQ.5 )THEN
      !maxI.EQ.5.AND.maxJ.EQ.0
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,1,0,X) = PA(K,ipass)
            ETIJ(K,ipass,1,1,0,X) = HPINV(K)
            ETIJ(K,ipass,0,2,0,X) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,3,0,X) = PA(K,ipass)*ETIJ(K,ipass,0,2,0,X)+ ETIJ(K,ipass,1,2,0,X)&
                 & - 2*ETIJ(K,ipass, 0,1,0,X)/TWOA(K)
            ETIJ(K,ipass,2,3,0,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PA(K,ipass)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,3,3,0,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,1,3,0,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PA(K,ipass)*ETIJ(K,ipass,1,2,0,X)&
                 & + 2*ETIJ(K,ipass,2,2,0,X) - 2*ETIJ(K,ipass, 1,1,0,X)/TWOA(K)
            ETIJ(K,ipass,0,4,0,X) = PA(K,ipass)*ETIJ(K,ipass,0,3,0,X)+ ETIJ(K,ipass,1,3,0,X)&
                 & - 3*ETIJ(K,ipass, 0,2,0,X)/TWOA(K)
            ETIJ(K,ipass,3,4,0,X) = HPINV(K)*ETIJ(K,ipass,2,3,0,X)+ PA(K,ipass)*ETIJ(K,ipass,3,3,0,X)
            ETIJ(K,ipass,4,4,0,X) = HPINV(K)*ETIJ(K,ipass,3,3,0,X)
            ETIJ(K,ipass,1,4,0,X) = HPINV(K)*ETIJ(K,ipass,0,3,0,X)+ PA(K,ipass)*ETIJ(K,ipass,1,3,0,X)&
                 & + 2*ETIJ(K,ipass,2,3,0,X) - 3*ETIJ(K,ipass, 1,2,0,X)/TWOA(K)
            ETIJ(K,ipass,2,4,0,X) = HPINV(K)*ETIJ(K,ipass,1,3,0,X)+ PA(K,ipass)*ETIJ(K,ipass,2,3,0,X)&
                 & + 3*ETIJ(K,ipass,3,3,0,X) - 3*ETIJ(K,ipass, 2,2,0,X)/TWOA(K)
            ETIJ(K,ipass,0,5,0,X) = PA(K,ipass)*ETIJ(K,ipass,0,4,0,X)+ ETIJ(K,ipass,1,4,0,X)&
                 & - 4*ETIJ(K,ipass, 0,3,0,X)/TWOA(K)
            ETIJ(K,ipass,4,5,0,X) = HPINV(K)*ETIJ(K,ipass,3,4,0,X)+ PA(K,ipass)*ETIJ(K,ipass,4,4,0,X)
            ETIJ(K,ipass,5,5,0,X) = HPINV(K)*ETIJ(K,ipass,4,4,0,X)
            ETIJ(K,ipass,1,5,0,X) = HPINV(K)*ETIJ(K,ipass,0,4,0,X)+ PA(K,ipass)*ETIJ(K,ipass,1,4,0,X)&
                 & + 2*ETIJ(K,ipass,2,4,0,X) - 4*ETIJ(K,ipass, 1,3,0,X)/TWOA(K)
            ETIJ(K,ipass,2,5,0,X) = HPINV(K)*ETIJ(K,ipass,1,4,0,X)+ PA(K,ipass)*ETIJ(K,ipass,2,4,0,X)&
                 & + 3*ETIJ(K,ipass,3,4,0,X) - 4*ETIJ(K,ipass, 2,3,0,X)/TWOA(K)
            ETIJ(K,ipass,3,5,0,X) = HPINV(K)*ETIJ(K,ipass,2,4,0,X)+ PA(K,ipass)*ETIJ(K,ipass,3,4,0,X)&
                 & + 4*ETIJ(K,ipass,4,4,0,X) - 4*ETIJ(K,ipass, 3,3,0,X)/TWOA(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.4)THEN
      !maxI.EQ.4.AND.maxJ.EQ.1
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,0,1,X) = PB(K,ipass)    
            ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
            ETIJ(K,ipass,0,1,0,X) = PA(K,ipass)
            ETIJ(K,ipass,1,1,0,X) = HPINV(K)
            ETIJ(K,ipass,0,1,1,X) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,2,0,X) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,2,1,X) = PB(K,ipass)*ETIJ(K,ipass, 0,2,0,X)+ ETIJ(K,ipass, 1,2,0,X)
            ETIJ(K,ipass,2,2,1,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PB(K,ipass)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,3,2,1,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,1,2,1,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PB(K,ipass)*ETIJ(K,ipass,1,2,0,X)&
                 & + 2*ETIJ(K,ipass,2,2,0,X)! - 0*ETIJ(K,ipass,1,2,*,X)/TWOB(K)
            ETIJ(K,ipass,0,3,0,X) = PA(K,ipass)*ETIJ(K,ipass,0,2,0,X)+ ETIJ(K,ipass,1,2,0,X)&
                 & - 2*ETIJ(K,ipass, 0,1,0,X)/TWOA(K)
            ETIJ(K,ipass,2,3,0,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PA(K,ipass)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,3,3,0,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,1,3,0,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PA(K,ipass)*ETIJ(K,ipass,1,2,0,X)&
                 & + 2*ETIJ(K,ipass,2,2,0,X) - 2*ETIJ(K,ipass, 1,1,0,X)/TWOA(K)
            ETIJ(K,ipass,0,3,1,X) = PB(K,ipass)*ETIJ(K,ipass, 0,3,0,X)+ ETIJ(K,ipass, 1,3,0,X)
            ETIJ(K,ipass,3,3,1,X) = HPINV(K)*ETIJ(K,ipass,2,3,0,X)+ PB(K,ipass)*ETIJ(K,ipass,3,3,0,X)
            ETIJ(K,ipass,4,3,1,X) = HPINV(K)*ETIJ(K,ipass,3,3,0,X)
            ETIJ(K,ipass,1,3,1,X)= HPINV(K)*ETIJ(K,ipass,0,3,0,X)+ PB(K,ipass)*ETIJ(K,ipass,1,3,0,X)&
                 & + 2*ETIJ(K,ipass,2,3,0,X)! - 0*ETIJ(K,ipass,1,3,*,X)/TWOB(K)
            ETIJ(K,ipass,2,3,1,X)= HPINV(K)*ETIJ(K,ipass,1,3,0,X)+ PB(K,ipass)*ETIJ(K,ipass,2,3,0,X)&
                 & + 3*ETIJ(K,ipass,3,3,0,X)! - 0*ETIJ(K,ipass,2,3,*,X)/TWOB(K)
            ETIJ(K,ipass,0,4,0,X) = PA(K,ipass)*ETIJ(K,ipass,0,3,0,X)+ ETIJ(K,ipass,1,3,0,X)&
                 & - 3*ETIJ(K,ipass, 0,2,0,X)/TWOA(K)
            ETIJ(K,ipass,3,4,0,X) = HPINV(K)*ETIJ(K,ipass,2,3,0,X)+ PA(K,ipass)*ETIJ(K,ipass,3,3,0,X)
            ETIJ(K,ipass,4,4,0,X) = HPINV(K)*ETIJ(K,ipass,3,3,0,X)
            ETIJ(K,ipass,1,4,0,X) = HPINV(K)*ETIJ(K,ipass,0,3,0,X)+ PA(K,ipass)*ETIJ(K,ipass,1,3,0,X)&
                 & + 2*ETIJ(K,ipass,2,3,0,X) - 3*ETIJ(K,ipass, 1,2,0,X)/TWOA(K)
            ETIJ(K,ipass,2,4,0,X) = HPINV(K)*ETIJ(K,ipass,1,3,0,X)+ PA(K,ipass)*ETIJ(K,ipass,2,3,0,X)&
                 & + 3*ETIJ(K,ipass,3,3,0,X) - 3*ETIJ(K,ipass, 2,2,0,X)/TWOA(K)
            ETIJ(K,ipass,0,4,1,X) = PB(K,ipass)*ETIJ(K,ipass, 0,4,0,X)+ ETIJ(K,ipass, 1,4,0,X)
            ETIJ(K,ipass,4,4,1,X) = HPINV(K)*ETIJ(K,ipass,3,4,0,X)+ PB(K,ipass)*ETIJ(K,ipass,4,4,0,X)
            ETIJ(K,ipass,5,4,1,X) = HPINV(K)*ETIJ(K,ipass,4,4,0,X)
            ETIJ(K,ipass,1,4,1,X) = HPINV(K)*ETIJ(K,ipass,0,4,0,X)+ PB(K,ipass)*ETIJ(K,ipass,1,4,0,X)&
                 & + 2*ETIJ(K,ipass,2,4,0,X)! - 0*ETIJ(K,ipass,1,4,*,X)/TWOB(K)
            ETIJ(K,ipass,2,4,1,X) = HPINV(K)*ETIJ(K,ipass,1,4,0,X)+ PB(K,ipass)*ETIJ(K,ipass,2,4,0,X)&
                 & + 3*ETIJ(K,ipass,3,4,0,X)! - 0*ETIJ(K,ipass,2,4,*,X)/TWOB(K)
            ETIJ(K,ipass,3,4,1,X) = HPINV(K)*ETIJ(K,ipass,2,4,0,X)+ PB(K,ipass)*ETIJ(K,ipass,3,4,0,X)&
                 & + 4*ETIJ(K,ipass,4,4,0,X)! - 0*ETIJ(K,ipass,3,4,*,X)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.3)THEN
      !maxI.EQ.3.AND.maxJ.EQ.2
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,0,1,X) = PB(K,ipass)    
            ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
            ETIJ(K,ipass,0,0,2,X) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,1,0,X) = PA(K,ipass)
            ETIJ(K,ipass,1,1,0,X) = HPINV(K)
            ETIJ(K,ipass,0,1,1,X) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,1,2,X) = PB(K,ipass)*ETIJ(K,ipass, 0,1,1,X)+ ETIJ(K,ipass, 1,1,1,X)&
                 & -1*ETIJ(K,ipass,0,1,0,X)/TWOB(K)
            ETIJ(K,ipass,2,1,2,X) = HPINV(K)*ETIJ(K,ipass,1,1,1,X)+ PB(K,ipass)*ETIJ(K,ipass,2,1,1,X)
            ETIJ(K,ipass,3,1,2,X) = HPINV(K)*ETIJ(K,ipass,2,1,1,X)
            ETIJ(K,ipass,1,1,2,X)= HPINV(K)*ETIJ(K,ipass,0,1,1,X)+ PB(K,ipass)*ETIJ(K,ipass,1,1,1,X)&
                 & + 2*ETIJ(K,ipass,2,1,1,X) - 1*ETIJ(K,ipass,1,1,0,X)/TWOB(K)
            ETIJ(K,ipass,0,2,0,X) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,2,1,X) = PB(K,ipass)*ETIJ(K,ipass, 0,2,0,X)+ ETIJ(K,ipass, 1,2,0,X)
            ETIJ(K,ipass,2,2,1,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PB(K,ipass)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,3,2,1,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,1,2,1,X)= HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PB(K,ipass)*ETIJ(K,ipass,1,2,0,X)&
                 & + 2*ETIJ(K,ipass,2,2,0,X)! - 0*ETIJ(K,ipass,1,2,*,X)/TWOB(K)
            ETIJ(K,ipass,0,2,2,X) = PB(K,ipass)*ETIJ(K,ipass, 0,2,1,X)+ ETIJ(K,ipass, 1,2,1,X)&
                 & -1*ETIJ(K,ipass,0,2,0,X)/TWOB(K)
            ETIJ(K,ipass,3,2,2,X) = HPINV(K)*ETIJ(K,ipass,2,2,1,X)+ PB(K,ipass)*ETIJ(K,ipass,3,2,1,X)
            ETIJ(K,ipass,4,2,2,X) = HPINV(K)*ETIJ(K,ipass,3,2,1,X)
            ETIJ(K,ipass,1,2,2,X)= HPINV(K)*ETIJ(K,ipass,0,2,1,X)+ PB(K,ipass)*ETIJ(K,ipass,1,2,1,X)&
                 & + 2*ETIJ(K,ipass,2,2,1,X) - 1*ETIJ(K,ipass,1,2,0,X)/TWOB(K)
            ETIJ(K,ipass,2,2,2,X)= HPINV(K)*ETIJ(K,ipass,1,2,1,X)+ PB(K,ipass)*ETIJ(K,ipass,2,2,1,X)&
                 & + 3*ETIJ(K,ipass,3,2,1,X) - 1*ETIJ(K,ipass,2,2,0,X)/TWOB(K)
            ETIJ(K,ipass,0,3,0,X) = PA(K,ipass)*ETIJ(K,ipass,0,2,0,X)+ ETIJ(K,ipass,1,2,0,X)&
                 & - 2*ETIJ(K,ipass, 0,1,0,X)/TWOA(K)
            ETIJ(K,ipass,2,3,0,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PA(K,ipass)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,3,3,0,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,1,3,0,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PA(K,ipass)*ETIJ(K,ipass,1,2,0,X)&
                 & + 2*ETIJ(K,ipass,2,2,0,X) - 2*ETIJ(K,ipass, 1,1,0,X)/TWOA(K)
            ETIJ(K,ipass,0,3,1,X) = PB(K,ipass)*ETIJ(K,ipass, 0,3,0,X)+ ETIJ(K,ipass, 1,3,0,X)
            ETIJ(K,ipass,3,3,1,X) = HPINV(K)*ETIJ(K,ipass,2,3,0,X)+ PB(K,ipass)*ETIJ(K,ipass,3,3,0,X)
            ETIJ(K,ipass,4,3,1,X) = HPINV(K)*ETIJ(K,ipass,3,3,0,X)
            ETIJ(K,ipass,1,3,1,X)= HPINV(K)*ETIJ(K,ipass,0,3,0,X)+ PB(K,ipass)*ETIJ(K,ipass,1,3,0,X)&
                 & + 2*ETIJ(K,ipass,2,3,0,X)! - 0*ETIJ(K,ipass,1,3,*,X)/TWOB(K)
            ETIJ(K,ipass,2,3,1,X)= HPINV(K)*ETIJ(K,ipass,1,3,0,X)+ PB(K,ipass)*ETIJ(K,ipass,2,3,0,X)&
                 & + 3*ETIJ(K,ipass,3,3,0,X)! - 0*ETIJ(K,ipass,2,3,*,X)/TWOB(K)
            ETIJ(K,ipass,0,3,2,X) = PB(K,ipass)*ETIJ(K,ipass, 0,3,1,X)+ ETIJ(K,ipass, 1,3,1,X)&
                 & -1*ETIJ(K,ipass,0,3,0,X)/TWOB(K)
            ETIJ(K,ipass,4,3,2,X) = HPINV(K)*ETIJ(K,ipass,3,3,1,X)+ PB(K,ipass)*ETIJ(K,ipass,4,3,1,X)
            ETIJ(K,ipass,5,3,2,X) = HPINV(K)*ETIJ(K,ipass,4,3,1,X)
            ETIJ(K,ipass,1,3,2,X)= HPINV(K)*ETIJ(K,ipass,0,3,1,X)+ PB(K,ipass)*ETIJ(K,ipass,1,3,1,X)&
                 & + 2*ETIJ(K,ipass,2,3,1,X) - 1*ETIJ(K,ipass,1,3,0,X)/TWOB(K)
            ETIJ(K,ipass,2,3,2,X)= HPINV(K)*ETIJ(K,ipass,1,3,1,X)+ PB(K,ipass)*ETIJ(K,ipass,2,3,1,X)&
                 & + 3*ETIJ(K,ipass,3,3,1,X) - 1*ETIJ(K,ipass,2,3,0,X)/TWOB(K)
            ETIJ(K,ipass,3,3,2,X)= HPINV(K)*ETIJ(K,ipass,2,3,1,X)+ PB(K,ipass)*ETIJ(K,ipass,3,3,1,X)&
                 & + 4*ETIJ(K,ipass,4,3,1,X) - 1*ETIJ(K,ipass,3,3,0,X)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.2)THEN
      !maxI.EQ.2.AND.maxJ.EQ.3
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,0,1,X) = PB(K,ipass)    
            ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
            ETIJ(K,ipass,0,0,2,X) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,0,3,X) = PB(K,ipass)*ETIJ(K,ipass, 0,0,2,X)+ ETIJ(K,ipass, 1,0,2,X)&
                 & -2*ETIJ(K,ipass,0,0,1,X)/TWOB(K)
            ETIJ(K,ipass,2,0,3,X) = HPINV(K)*ETIJ(K,ipass,1,0,2,X)+ PB(K,ipass)*ETIJ(K,ipass,2,0,2,X)
            ETIJ(K,ipass,3,0,3,X) = HPINV(K)*ETIJ(K,ipass,2,0,2,X)
            ETIJ(K,ipass,1,0,3,X)= HPINV(K)*ETIJ(K,ipass,0,0,2,X)+ PB(K,ipass)*ETIJ(K,ipass,1,0,2,X)&
                 & + 2*ETIJ(K,ipass,2,0,2,X) - 2*ETIJ(K,ipass,1,0,1,X)/TWOB(K)
            ETIJ(K,ipass,0,1,0,X) = PA(K,ipass)
            ETIJ(K,ipass,1,1,0,X) = HPINV(K)
            ETIJ(K,ipass,0,1,1,X) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,1,2,X) = PB(K,ipass)*ETIJ(K,ipass, 0,1,1,X)+ ETIJ(K,ipass, 1,1,1,X)&
                 & -1*ETIJ(K,ipass,0,1,0,X)/TWOB(K)
            ETIJ(K,ipass,2,1,2,X) = HPINV(K)*ETIJ(K,ipass,1,1,1,X)+ PB(K,ipass)*ETIJ(K,ipass,2,1,1,X)
            ETIJ(K,ipass,3,1,2,X) = HPINV(K)*ETIJ(K,ipass,2,1,1,X)
            ETIJ(K,ipass,1,1,2,X)= HPINV(K)*ETIJ(K,ipass,0,1,1,X)+ PB(K,ipass)*ETIJ(K,ipass,1,1,1,X)&
                 & + 2*ETIJ(K,ipass,2,1,1,X) - 1*ETIJ(K,ipass,1,1,0,X)/TWOB(K)
            ETIJ(K,ipass,0,1,3,X) = PB(K,ipass)*ETIJ(K,ipass, 0,1,2,X)+ ETIJ(K,ipass, 1,1,2,X)&
                 & -2*ETIJ(K,ipass,0,1,1,X)/TWOB(K)
            ETIJ(K,ipass,3,1,3,X) = HPINV(K)*ETIJ(K,ipass,2,1,2,X)+ PB(K,ipass)*ETIJ(K,ipass,3,1,2,X)
            ETIJ(K,ipass,4,1,3,X) = HPINV(K)*ETIJ(K,ipass,3,1,2,X)
            ETIJ(K,ipass,1,1,3,X)= HPINV(K)*ETIJ(K,ipass,0,1,2,X)+ PB(K,ipass)*ETIJ(K,ipass,1,1,2,X)&
                 & + 2*ETIJ(K,ipass,2,1,2,X) - 2*ETIJ(K,ipass,1,1,1,X)/TWOB(K)
            ETIJ(K,ipass,2,1,3,X)= HPINV(K)*ETIJ(K,ipass,1,1,2,X)+ PB(K,ipass)*ETIJ(K,ipass,2,1,2,X)&
                 & + 3*ETIJ(K,ipass,3,1,2,X) - 2*ETIJ(K,ipass,2,1,1,X)/TWOB(K)
            ETIJ(K,ipass,0,2,0,X) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,2,1,X) = PB(K,ipass)*ETIJ(K,ipass, 0,2,0,X)+ ETIJ(K,ipass, 1,2,0,X)
            ETIJ(K,ipass,2,2,1,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PB(K,ipass)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,3,2,1,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,1,2,1,X)= HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PB(K,ipass)*ETIJ(K,ipass,1,2,0,X)&
                 & + 2*ETIJ(K,ipass,2,2,0,X)! - 0*ETIJ(K,ipass,1,2,*,X)/TWOB(K)
            ETIJ(K,ipass,0,2,2,X) = PB(K,ipass)*ETIJ(K,ipass, 0,2,1,X)+ ETIJ(K,ipass, 1,2,1,X)&
                 & -1*ETIJ(K,ipass,0,2,0,X)/TWOB(K)
            ETIJ(K,ipass,3,2,2,X) = HPINV(K)*ETIJ(K,ipass,2,2,1,X)+ PB(K,ipass)*ETIJ(K,ipass,3,2,1,X)
            ETIJ(K,ipass,4,2,2,X) = HPINV(K)*ETIJ(K,ipass,3,2,1,X)
            ETIJ(K,ipass,1,2,2,X)= HPINV(K)*ETIJ(K,ipass,0,2,1,X)+ PB(K,ipass)*ETIJ(K,ipass,1,2,1,X)&
                 & + 2*ETIJ(K,ipass,2,2,1,X) - 1*ETIJ(K,ipass,1,2,0,X)/TWOB(K)
            ETIJ(K,ipass,2,2,2,X)= HPINV(K)*ETIJ(K,ipass,1,2,1,X)+ PB(K,ipass)*ETIJ(K,ipass,2,2,1,X)&
                 & + 3*ETIJ(K,ipass,3,2,1,X) - 1*ETIJ(K,ipass,2,2,0,X)/TWOB(K)
            ETIJ(K,ipass,0,2,3,X) = PB(K,ipass)*ETIJ(K,ipass, 0,2,2,X)+ ETIJ(K,ipass, 1,2,2,X)&
                 & -2*ETIJ(K,ipass,0,2,1,X)/TWOB(K)
            ETIJ(K,ipass,4,2,3,X) = HPINV(K)*ETIJ(K,ipass,3,2,2,X)+ PB(K,ipass)*ETIJ(K,ipass,4,2,2,X)
            ETIJ(K,ipass,5,2,3,X) = HPINV(K)*ETIJ(K,ipass,4,2,2,X)
            ETIJ(K,ipass,1,2,3,X)= HPINV(K)*ETIJ(K,ipass,0,2,2,X)+ PB(K,ipass)*ETIJ(K,ipass,1,2,2,X)&
                 & + 2*ETIJ(K,ipass,2,2,2,X) - 2*ETIJ(K,ipass,1,2,1,X)/TWOB(K)
            ETIJ(K,ipass,2,2,3,X)= HPINV(K)*ETIJ(K,ipass,1,2,2,X)+ PB(K,ipass)*ETIJ(K,ipass,2,2,2,X)&
                 & + 3*ETIJ(K,ipass,3,2,2,X) - 2*ETIJ(K,ipass,2,2,1,X)/TWOB(K)
            ETIJ(K,ipass,3,2,3,X)= HPINV(K)*ETIJ(K,ipass,2,2,2,X)+ PB(K,ipass)*ETIJ(K,ipass,3,2,2,X)&
                 & + 4*ETIJ(K,ipass,4,2,2,X) - 2*ETIJ(K,ipass,3,2,1,X)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.1)THEN
      !maxI.EQ.1.AND.maxJ.EQ.4
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,0,1,X) = PB(K,ipass)    
            ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
            ETIJ(K,ipass,0,0,2,X) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,0,3,X) = PB(K,ipass)*ETIJ(K,ipass, 0,0,2,X)+ ETIJ(K,ipass, 1,0,2,X)&
                 & -2*ETIJ(K,ipass,0,0,1,X)/TWOB(K)
            ETIJ(K,ipass,2,0,3,X) = HPINV(K)*ETIJ(K,ipass,1,0,2,X)+ PB(K,ipass)*ETIJ(K,ipass,2,0,2,X)
            ETIJ(K,ipass,3,0,3,X) = HPINV(K)*ETIJ(K,ipass,2,0,2,X)
            ETIJ(K,ipass,1,0,3,X)= HPINV(K)*ETIJ(K,ipass,0,0,2,X)+ PB(K,ipass)*ETIJ(K,ipass,1,0,2,X)&
                 & + 2*ETIJ(K,ipass,2,0,2,X) - 2*ETIJ(K,ipass,1,0,1,X)/TWOB(K)
            ETIJ(K,ipass,0,0,4,X) = PB(K,ipass)*ETIJ(K,ipass, 0,0,3,X)+ ETIJ(K,ipass, 1,0,3,X)&
                 & -3*ETIJ(K,ipass,0,0,2,X)/TWOB(K)
            ETIJ(K,ipass,3,0,4,X) = HPINV(K)*ETIJ(K,ipass,2,0,3,X)+ PB(K,ipass)*ETIJ(K,ipass,3,0,3,X)
            ETIJ(K,ipass,4,0,4,X) = HPINV(K)*ETIJ(K,ipass,3,0,3,X)
            ETIJ(K,ipass,1,0,4,X)= HPINV(K)*ETIJ(K,ipass,0,0,3,X)+ PB(K,ipass)*ETIJ(K,ipass,1,0,3,X)&
                 & + 2*ETIJ(K,ipass,2,0,3,X) - 3*ETIJ(K,ipass,1,0,2,X)/TWOB(K)
            ETIJ(K,ipass,2,0,4,X)= HPINV(K)*ETIJ(K,ipass,1,0,3,X)+ PB(K,ipass)*ETIJ(K,ipass,2,0,3,X)&
                 & + 3*ETIJ(K,ipass,3,0,3,X) - 3*ETIJ(K,ipass,2,0,2,X)/TWOB(K)
            ETIJ(K,ipass,0,1,0,X) = PA(K,ipass)
            ETIJ(K,ipass,1,1,0,X) = HPINV(K)
            ETIJ(K,ipass,0,1,1,X) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,1,2,X) = PB(K,ipass)*ETIJ(K,ipass, 0,1,1,X)+ ETIJ(K,ipass, 1,1,1,X)&
                 & -1*ETIJ(K,ipass,0,1,0,X)/TWOB(K)
            ETIJ(K,ipass,2,1,2,X) = HPINV(K)*ETIJ(K,ipass,1,1,1,X)+ PB(K,ipass)*ETIJ(K,ipass,2,1,1,X)
            ETIJ(K,ipass,3,1,2,X) = HPINV(K)*ETIJ(K,ipass,2,1,1,X)
            ETIJ(K,ipass,1,1,2,X)= HPINV(K)*ETIJ(K,ipass,0,1,1,X)+ PB(K,ipass)*ETIJ(K,ipass,1,1,1,X)&
                 & + 2*ETIJ(K,ipass,2,1,1,X) - 1*ETIJ(K,ipass,1,1,0,X)/TWOB(K)
            ETIJ(K,ipass,0,1,3,X) = PB(K,ipass)*ETIJ(K,ipass, 0,1,2,X)+ ETIJ(K,ipass, 1,1,2,X)&
                 & -2*ETIJ(K,ipass,0,1,1,X)/TWOB(K)
            ETIJ(K,ipass,3,1,3,X) = HPINV(K)*ETIJ(K,ipass,2,1,2,X)+ PB(K,ipass)*ETIJ(K,ipass,3,1,2,X)
            ETIJ(K,ipass,4,1,3,X) = HPINV(K)*ETIJ(K,ipass,3,1,2,X)
            ETIJ(K,ipass,1,1,3,X)= HPINV(K)*ETIJ(K,ipass,0,1,2,X)+ PB(K,ipass)*ETIJ(K,ipass,1,1,2,X)&
                 & + 2*ETIJ(K,ipass,2,1,2,X) - 2*ETIJ(K,ipass,1,1,1,X)/TWOB(K)
            ETIJ(K,ipass,2,1,3,X)= HPINV(K)*ETIJ(K,ipass,1,1,2,X)+ PB(K,ipass)*ETIJ(K,ipass,2,1,2,X)&
                 & + 3*ETIJ(K,ipass,3,1,2,X) - 2*ETIJ(K,ipass,2,1,1,X)/TWOB(K)
            ETIJ(K,ipass,0,1,4,X) = PB(K,ipass)*ETIJ(K,ipass, 0,1,3,X)+ ETIJ(K,ipass, 1,1,3,X)&
                 & -3*ETIJ(K,ipass,0,1,2,X)/TWOB(K)
            ETIJ(K,ipass,4,1,4,X) = HPINV(K)*ETIJ(K,ipass,3,1,3,X)+ PB(K,ipass)*ETIJ(K,ipass,4,1,3,X)
            ETIJ(K,ipass,5,1,4,X) = HPINV(K)*ETIJ(K,ipass,4,1,3,X)
            ETIJ(K,ipass,1,1,4,X)= HPINV(K)*ETIJ(K,ipass,0,1,3,X)+ PB(K,ipass)*ETIJ(K,ipass,1,1,3,X)&
                 & + 2*ETIJ(K,ipass,2,1,3,X) - 3*ETIJ(K,ipass,1,1,2,X)/TWOB(K)
            ETIJ(K,ipass,2,1,4,X)= HPINV(K)*ETIJ(K,ipass,1,1,3,X)+ PB(K,ipass)*ETIJ(K,ipass,2,1,3,X)&
                 & + 3*ETIJ(K,ipass,3,1,3,X) - 3*ETIJ(K,ipass,2,1,2,X)/TWOB(K)
            ETIJ(K,ipass,3,1,4,X)= HPINV(K)*ETIJ(K,ipass,2,1,3,X)+ PB(K,ipass)*ETIJ(K,ipass,3,1,3,X)&
                 & + 4*ETIJ(K,ipass,4,1,3,X) - 3*ETIJ(K,ipass,3,1,2,X)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.0)THEN
      !maxI.EQ.0.AND.maxJ.EQ.5
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,0,1,X) = PB(K,ipass)    
            ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
            ETIJ(K,ipass,0,0,2,X) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,0,3,X) = PB(K,ipass)*ETIJ(K,ipass, 0,0,2,X)+ ETIJ(K,ipass, 1,0,2,X)&
                 & -2*ETIJ(K,ipass,0,0,1,X)/TWOB(K)
            ETIJ(K,ipass,2,0,3,X) = HPINV(K)*ETIJ(K,ipass,1,0,2,X)+ PB(K,ipass)*ETIJ(K,ipass,2,0,2,X)
            ETIJ(K,ipass,3,0,3,X) = HPINV(K)*ETIJ(K,ipass,2,0,2,X)
            ETIJ(K,ipass,1,0,3,X)= HPINV(K)*ETIJ(K,ipass,0,0,2,X)+ PB(K,ipass)*ETIJ(K,ipass,1,0,2,X)&
                 & + 2*ETIJ(K,ipass,2,0,2,X) - 2*ETIJ(K,ipass,1,0,1,X)/TWOB(K)
            ETIJ(K,ipass,0,0,4,X) = PB(K,ipass)*ETIJ(K,ipass, 0,0,3,X)+ ETIJ(K,ipass, 1,0,3,X)&
                 & -3*ETIJ(K,ipass,0,0,2,X)/TWOB(K)
            ETIJ(K,ipass,3,0,4,X) = HPINV(K)*ETIJ(K,ipass,2,0,3,X)+ PB(K,ipass)*ETIJ(K,ipass,3,0,3,X)
            ETIJ(K,ipass,4,0,4,X) = HPINV(K)*ETIJ(K,ipass,3,0,3,X)
            ETIJ(K,ipass,1,0,4,X)= HPINV(K)*ETIJ(K,ipass,0,0,3,X)+ PB(K,ipass)*ETIJ(K,ipass,1,0,3,X)&
                 & + 2*ETIJ(K,ipass,2,0,3,X) - 3*ETIJ(K,ipass,1,0,2,X)/TWOB(K)
            ETIJ(K,ipass,2,0,4,X)= HPINV(K)*ETIJ(K,ipass,1,0,3,X)+ PB(K,ipass)*ETIJ(K,ipass,2,0,3,X)&
                 & + 3*ETIJ(K,ipass,3,0,3,X) - 3*ETIJ(K,ipass,2,0,2,X)/TWOB(K)
            ETIJ(K,ipass,0,0,5,X) = PB(K,ipass)*ETIJ(K,ipass, 0,0,4,X)+ ETIJ(K,ipass, 1,0,4,X)&
                 & -4*ETIJ(K,ipass,0,0,3,X)/TWOB(K)
            ETIJ(K,ipass,4,0,5,X) = HPINV(K)*ETIJ(K,ipass,3,0,4,X)+ PB(K,ipass)*ETIJ(K,ipass,4,0,4,X)
            ETIJ(K,ipass,5,0,5,X) = HPINV(K)*ETIJ(K,ipass,4,0,4,X)
            ETIJ(K,ipass,1,0,5,X)= HPINV(K)*ETIJ(K,ipass,0,0,4,X)+ PB(K,ipass)*ETIJ(K,ipass,1,0,4,X)&
                 & + 2*ETIJ(K,ipass,2,0,4,X) - 4*ETIJ(K,ipass,1,0,3,X)/TWOB(K)
            ETIJ(K,ipass,2,0,5,X)= HPINV(K)*ETIJ(K,ipass,1,0,4,X)+ PB(K,ipass)*ETIJ(K,ipass,2,0,4,X)&
                 & + 3*ETIJ(K,ipass,3,0,4,X) - 4*ETIJ(K,ipass,2,0,3,X)/TWOB(K)
            ETIJ(K,ipass,3,0,5,X)= HPINV(K)*ETIJ(K,ipass,2,0,4,X)+ PB(K,ipass)*ETIJ(K,ipass,3,0,4,X)&
                 & + 4*ETIJ(K,ipass,4,0,4,X) - 4*ETIJ(K,ipass,3,0,3,X)/TWOB(K)
         ENDDO
      ENDDO
   ENDIF
ELSEIF(maxIJ.EQ.6 )THEN !FF
   IF(maxI.EQ.6 )THEN
      !maxI.EQ.6.AND.maxJ.EQ.0
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,1,0,X) = PA(K,ipass)
            ETIJ(K,ipass,1,1,0,X) = HPINV(K)
            ETIJ(K,ipass,0,2,0,X) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,3,0,X) = PA(K,ipass)*ETIJ(K,ipass,0,2,0,X)+ ETIJ(K,ipass,1,2,0,X)&
                 & - 2*ETIJ(K,ipass, 0,1,0,X)/TWOA(K)
            ETIJ(K,ipass,2,3,0,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PA(K,ipass)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,3,3,0,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,1,3,0,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PA(K,ipass)*ETIJ(K,ipass,1,2,0,X)&
                 & + 2*ETIJ(K,ipass,2,2,0,X) - 2*ETIJ(K,ipass, 1,1,0,X)/TWOA(K)
            ETIJ(K,ipass,0,4,0,X) = PA(K,ipass)*ETIJ(K,ipass,0,3,0,X)+ ETIJ(K,ipass,1,3,0,X)&
                 & - 3*ETIJ(K,ipass, 0,2,0,X)/TWOA(K)
            ETIJ(K,ipass,3,4,0,X) = HPINV(K)*ETIJ(K,ipass,2,3,0,X)+ PA(K,ipass)*ETIJ(K,ipass,3,3,0,X)
            ETIJ(K,ipass,4,4,0,X) = HPINV(K)*ETIJ(K,ipass,3,3,0,X)
            ETIJ(K,ipass,1,4,0,X) = HPINV(K)*ETIJ(K,ipass,0,3,0,X)+ PA(K,ipass)*ETIJ(K,ipass,1,3,0,X)&
                 & + 2*ETIJ(K,ipass,2,3,0,X) - 3*ETIJ(K,ipass, 1,2,0,X)/TWOA(K)
            ETIJ(K,ipass,2,4,0,X) = HPINV(K)*ETIJ(K,ipass,1,3,0,X)+ PA(K,ipass)*ETIJ(K,ipass,2,3,0,X)&
                 & + 3*ETIJ(K,ipass,3,3,0,X) - 3*ETIJ(K,ipass, 2,2,0,X)/TWOA(K)
            ETIJ(K,ipass,0,5,0,X) = PA(K,ipass)*ETIJ(K,ipass,0,4,0,X)+ ETIJ(K,ipass,1,4,0,X)&
                 & - 4*ETIJ(K,ipass, 0,3,0,X)/TWOA(K)
            ETIJ(K,ipass,4,5,0,X) = HPINV(K)*ETIJ(K,ipass,3,4,0,X)+ PA(K,ipass)*ETIJ(K,ipass,4,4,0,X)
            ETIJ(K,ipass,5,5,0,X) = HPINV(K)*ETIJ(K,ipass,4,4,0,X)
            ETIJ(K,ipass,1,5,0,X) = HPINV(K)*ETIJ(K,ipass,0,4,0,X)+ PA(K,ipass)*ETIJ(K,ipass,1,4,0,X)&
                 & + 2*ETIJ(K,ipass,2,4,0,X) - 4*ETIJ(K,ipass, 1,3,0,X)/TWOA(K)
            ETIJ(K,ipass,2,5,0,X) = HPINV(K)*ETIJ(K,ipass,1,4,0,X)+ PA(K,ipass)*ETIJ(K,ipass,2,4,0,X)&
                 & + 3*ETIJ(K,ipass,3,4,0,X) - 4*ETIJ(K,ipass, 2,3,0,X)/TWOA(K)
            ETIJ(K,ipass,3,5,0,X) = HPINV(K)*ETIJ(K,ipass,2,4,0,X)+ PA(K,ipass)*ETIJ(K,ipass,3,4,0,X)&
                 & + 4*ETIJ(K,ipass,4,4,0,X) - 4*ETIJ(K,ipass, 3,3,0,X)/TWOA(K)
            ETIJ(K,ipass,0,6,0,X) = PA(K,ipass)*ETIJ(K,ipass,0,5,0,X)+ ETIJ(K,ipass,1,5,0,X)&
                 & - 5*ETIJ(K,ipass, 0,4,0,X)/TWOA(K)
            ETIJ(K,ipass,5,6,0,X) = HPINV(K)*ETIJ(K,ipass,4,5,0,X)+ PA(K,ipass)*ETIJ(K,ipass,5,5,0,X)
            ETIJ(K,ipass,6,6,0,X) = HPINV(K)*ETIJ(K,ipass,5,5,0,X)
            ETIJ(K,ipass,1,6,0,X) = HPINV(K)*ETIJ(K,ipass,0,5,0,X)+ PA(K,ipass)*ETIJ(K,ipass,1,5,0,X)&
                 & + 2*ETIJ(K,ipass,2,5,0,X) - 5*ETIJ(K,ipass, 1,4,0,X)/TWOA(K)
            ETIJ(K,ipass,2,6,0,X) = HPINV(K)*ETIJ(K,ipass,1,5,0,X)+ PA(K,ipass)*ETIJ(K,ipass,2,5,0,X)&
                 & + 3*ETIJ(K,ipass,3,5,0,X) - 5*ETIJ(K,ipass, 2,4,0,X)/TWOA(K)
            ETIJ(K,ipass,3,6,0,X) = HPINV(K)*ETIJ(K,ipass,2,5,0,X)+ PA(K,ipass)*ETIJ(K,ipass,3,5,0,X)&
                 & + 4*ETIJ(K,ipass,4,5,0,X) - 5*ETIJ(K,ipass, 3,4,0,X)/TWOA(K)
            ETIJ(K,ipass,4,6,0,X) = HPINV(K)*ETIJ(K,ipass,3,5,0,X)+ PA(K,ipass)*ETIJ(K,ipass,4,5,0,X)&
                 & + 5*ETIJ(K,ipass,5,5,0,X) - 5*ETIJ(K,ipass, 4,4,0,X)/TWOA(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.5)THEN
      !maxI.EQ.5.AND.maxJ.EQ.1
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,0,1,X) = PB(K,ipass)    
            ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
            ETIJ(K,ipass,0,1,0,X) = PA(K,ipass)
            ETIJ(K,ipass,1,1,0,X) = HPINV(K)
            ETIJ(K,ipass,0,1,1,X) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,2,0,X) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,2,1,X) = PB(K,ipass)*ETIJ(K,ipass, 0,2,0,X)+ ETIJ(K,ipass, 1,2,0,X)
            ETIJ(K,ipass,2,2,1,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PB(K,ipass)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,3,2,1,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,1,2,1,X)= HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PB(K,ipass)*ETIJ(K,ipass,1,2,0,X)&
                 & + 2*ETIJ(K,ipass,2,2,0,X)! - 0*ETIJ(K,ipass,1,2,*,X)/TWOB(K)
            ETIJ(K,ipass,0,3,0,X) = PA(K,ipass)*ETIJ(K,ipass,0,2,0,X)+ ETIJ(K,ipass,1,2,0,X)&
                 & - 2*ETIJ(K,ipass, 0,1,0,X)/TWOA(K)
            ETIJ(K,ipass,2,3,0,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PA(K,ipass)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,3,3,0,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,1,3,0,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PA(K,ipass)*ETIJ(K,ipass,1,2,0,X)&
                 & + 2*ETIJ(K,ipass,2,2,0,X) - 2*ETIJ(K,ipass, 1,1,0,X)/TWOA(K)
            ETIJ(K,ipass,0,3,1,X) = PB(K,ipass)*ETIJ(K,ipass, 0,3,0,X)+ ETIJ(K,ipass, 1,3,0,X)
            ETIJ(K,ipass,3,3,1,X) = HPINV(K)*ETIJ(K,ipass,2,3,0,X)+ PB(K,ipass)*ETIJ(K,ipass,3,3,0,X)
            ETIJ(K,ipass,4,3,1,X) = HPINV(K)*ETIJ(K,ipass,3,3,0,X)
            ETIJ(K,ipass,1,3,1,X)= HPINV(K)*ETIJ(K,ipass,0,3,0,X)+ PB(K,ipass)*ETIJ(K,ipass,1,3,0,X)&
                 & + 2*ETIJ(K,ipass,2,3,0,X)! - 0*ETIJ(K,ipass,1,3,*,X)/TWOB(K)
            ETIJ(K,ipass,2,3,1,X)= HPINV(K)*ETIJ(K,ipass,1,3,0,X)+ PB(K,ipass)*ETIJ(K,ipass,2,3,0,X)&
                 & + 3*ETIJ(K,ipass,3,3,0,X)! - 0*ETIJ(K,ipass,2,3,*,X)/TWOB(K)
            ETIJ(K,ipass,0,4,0,X) = PA(K,ipass)*ETIJ(K,ipass,0,3,0,X)+ ETIJ(K,ipass,1,3,0,X)&
                 & - 3*ETIJ(K,ipass, 0,2,0,X)/TWOA(K)
            ETIJ(K,ipass,3,4,0,X) = HPINV(K)*ETIJ(K,ipass,2,3,0,X)+ PA(K,ipass)*ETIJ(K,ipass,3,3,0,X)
            ETIJ(K,ipass,4,4,0,X) = HPINV(K)*ETIJ(K,ipass,3,3,0,X)
            ETIJ(K,ipass,1,4,0,X) = HPINV(K)*ETIJ(K,ipass,0,3,0,X)+ PA(K,ipass)*ETIJ(K,ipass,1,3,0,X)&
                 & + 2*ETIJ(K,ipass,2,3,0,X) - 3*ETIJ(K,ipass, 1,2,0,X)/TWOA(K)
            ETIJ(K,ipass,2,4,0,X) = HPINV(K)*ETIJ(K,ipass,1,3,0,X)+ PA(K,ipass)*ETIJ(K,ipass,2,3,0,X)&
                 & + 3*ETIJ(K,ipass,3,3,0,X) - 3*ETIJ(K,ipass, 2,2,0,X)/TWOA(K)
            ETIJ(K,ipass,0,4,1,X) = PB(K,ipass)*ETIJ(K,ipass, 0,4,0,X)+ ETIJ(K,ipass, 1,4,0,X)
            ETIJ(K,ipass,4,4,1,X) = HPINV(K)*ETIJ(K,ipass,3,4,0,X)+ PB(K,ipass)*ETIJ(K,ipass,4,4,0,X)
            ETIJ(K,ipass,5,4,1,X) = HPINV(K)*ETIJ(K,ipass,4,4,0,X)
            ETIJ(K,ipass,1,4,1,X)= HPINV(K)*ETIJ(K,ipass,0,4,0,X)+ PB(K,ipass)*ETIJ(K,ipass,1,4,0,X)&
                 & + 2*ETIJ(K,ipass,2,4,0,X)! - 0*ETIJ(K,ipass,1,4,*,X)/TWOB(K)
            ETIJ(K,ipass,2,4,1,X)= HPINV(K)*ETIJ(K,ipass,1,4,0,X)+ PB(K,ipass)*ETIJ(K,ipass,2,4,0,X)&
                 & + 3*ETIJ(K,ipass,3,4,0,X)! - 0*ETIJ(K,ipass,2,4,*,X)/TWOB(K)
            ETIJ(K,ipass,3,4,1,X)= HPINV(K)*ETIJ(K,ipass,2,4,0,X)+ PB(K,ipass)*ETIJ(K,ipass,3,4,0,X)&
                 & + 4*ETIJ(K,ipass,4,4,0,X)! - 0*ETIJ(K,ipass,3,4,*,X)/TWOB(K)
            ETIJ(K,ipass,0,5,0,X) = PA(K,ipass)*ETIJ(K,ipass,0,4,0,X)+ ETIJ(K,ipass,1,4,0,X)&
                 & - 4*ETIJ(K,ipass, 0,3,0,X)/TWOA(K)
            ETIJ(K,ipass,4,5,0,X) = HPINV(K)*ETIJ(K,ipass,3,4,0,X)+ PA(K,ipass)*ETIJ(K,ipass,4,4,0,X)
            ETIJ(K,ipass,5,5,0,X) = HPINV(K)*ETIJ(K,ipass,4,4,0,X)
            ETIJ(K,ipass,1,5,0,X) = HPINV(K)*ETIJ(K,ipass,0,4,0,X)+ PA(K,ipass)*ETIJ(K,ipass,1,4,0,X)&
                 & + 2*ETIJ(K,ipass,2,4,0,X) - 4*ETIJ(K,ipass, 1,3,0,X)/TWOA(K)
            ETIJ(K,ipass,2,5,0,X) = HPINV(K)*ETIJ(K,ipass,1,4,0,X)+ PA(K,ipass)*ETIJ(K,ipass,2,4,0,X)&
                 & + 3*ETIJ(K,ipass,3,4,0,X) - 4*ETIJ(K,ipass, 2,3,0,X)/TWOA(K)
            ETIJ(K,ipass,3,5,0,X) = HPINV(K)*ETIJ(K,ipass,2,4,0,X)+ PA(K,ipass)*ETIJ(K,ipass,3,4,0,X)&
                 & + 4*ETIJ(K,ipass,4,4,0,X) - 4*ETIJ(K,ipass, 3,3,0,X)/TWOA(K)
            ETIJ(K,ipass,0,5,1,X) = PB(K,ipass)*ETIJ(K,ipass, 0,5,0,X)+ ETIJ(K,ipass, 1,5,0,X)
            ETIJ(K,ipass,5,5,1,X) = HPINV(K)*ETIJ(K,ipass,4,5,0,X)+ PB(K,ipass)*ETIJ(K,ipass,5,5,0,X)
            ETIJ(K,ipass,6,5,1,X) = HPINV(K)*ETIJ(K,ipass,5,5,0,X)
            ETIJ(K,ipass,1,5,1,X)= HPINV(K)*ETIJ(K,ipass,0,5,0,X)+ PB(K,ipass)*ETIJ(K,ipass,1,5,0,X)&
                 & + 2*ETIJ(K,ipass,2,5,0,X)! - 0*ETIJ(K,ipass,1,5,*,X)/TWOB(K)
            ETIJ(K,ipass,2,5,1,X)= HPINV(K)*ETIJ(K,ipass,1,5,0,X)+ PB(K,ipass)*ETIJ(K,ipass,2,5,0,X)&
                 & + 3*ETIJ(K,ipass,3,5,0,X)! - 0*ETIJ(K,ipass,2,5,*,X)/TWOB(K)
            ETIJ(K,ipass,3,5,1,X)= HPINV(K)*ETIJ(K,ipass,2,5,0,X)+ PB(K,ipass)*ETIJ(K,ipass,3,5,0,X)&
                 & + 4*ETIJ(K,ipass,4,5,0,X)! - 0*ETIJ(K,ipass,3,5,*,X)/TWOB(K)
            ETIJ(K,ipass,4,5,1,X)= HPINV(K)*ETIJ(K,ipass,3,5,0,X)+ PB(K,ipass)*ETIJ(K,ipass,4,5,0,X)&
                 & + 5*ETIJ(K,ipass,5,5,0,X)! - 0*ETIJ(K,ipass,4,5,*,X)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.4)THEN
      !maxI.EQ.4.AND.maxJ.EQ.2
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,0,1,X) = PB(K,ipass)    
            ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
            ETIJ(K,ipass,0,0,2,X) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,1,0,X) = PA(K,ipass)
            ETIJ(K,ipass,1,1,0,X) = HPINV(K)
            ETIJ(K,ipass,0,1,1,X) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,1,2,X) = PB(K,ipass)*ETIJ(K,ipass, 0,1,1,X)+ ETIJ(K,ipass, 1,1,1,X)&
                 & -1*ETIJ(K,ipass,0,1,0,X)/TWOB(K)
            ETIJ(K,ipass,2,1,2,X) = HPINV(K)*ETIJ(K,ipass,1,1,1,X)+ PB(K,ipass)*ETIJ(K,ipass,2,1,1,X)
            ETIJ(K,ipass,3,1,2,X) = HPINV(K)*ETIJ(K,ipass,2,1,1,X)
            ETIJ(K,ipass,1,1,2,X)= HPINV(K)*ETIJ(K,ipass,0,1,1,X)+ PB(K,ipass)*ETIJ(K,ipass,1,1,1,X)&
                 & + 2*ETIJ(K,ipass,2,1,1,X) - 1*ETIJ(K,ipass,1,1,0,X)/TWOB(K)
            ETIJ(K,ipass,0,2,0,X) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,2,1,X) = PB(K,ipass)*ETIJ(K,ipass, 0,2,0,X)+ ETIJ(K,ipass, 1,2,0,X)
            ETIJ(K,ipass,2,2,1,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PB(K,ipass)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,3,2,1,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,1,2,1,X)= HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PB(K,ipass)*ETIJ(K,ipass,1,2,0,X)&
                 & + 2*ETIJ(K,ipass,2,2,0,X)! - 0*ETIJ(K,ipass,1,2,*,X)/TWOB(K)
            ETIJ(K,ipass,0,2,2,X) = PB(K,ipass)*ETIJ(K,ipass, 0,2,1,X)+ ETIJ(K,ipass, 1,2,1,X)&
                 & -1*ETIJ(K,ipass,0,2,0,X)/TWOB(K)
            ETIJ(K,ipass,3,2,2,X) = HPINV(K)*ETIJ(K,ipass,2,2,1,X)+ PB(K,ipass)*ETIJ(K,ipass,3,2,1,X)
            ETIJ(K,ipass,4,2,2,X) = HPINV(K)*ETIJ(K,ipass,3,2,1,X)
            ETIJ(K,ipass,1,2,2,X)= HPINV(K)*ETIJ(K,ipass,0,2,1,X)+ PB(K,ipass)*ETIJ(K,ipass,1,2,1,X)&
                 & + 2*ETIJ(K,ipass,2,2,1,X) - 1*ETIJ(K,ipass,1,2,0,X)/TWOB(K)
            ETIJ(K,ipass,2,2,2,X)= HPINV(K)*ETIJ(K,ipass,1,2,1,X)+ PB(K,ipass)*ETIJ(K,ipass,2,2,1,X)&
                 & + 3*ETIJ(K,ipass,3,2,1,X) - 1*ETIJ(K,ipass,2,2,0,X)/TWOB(K)
            ETIJ(K,ipass,0,3,0,X) = PA(K,ipass)*ETIJ(K,ipass,0,2,0,X)+ ETIJ(K,ipass,1,2,0,X)&
                 & - 2*ETIJ(K,ipass, 0,1,0,X)/TWOA(K)
            ETIJ(K,ipass,2,3,0,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PA(K,ipass)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,3,3,0,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,1,3,0,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PA(K,ipass)*ETIJ(K,ipass,1,2,0,X)&
                 & + 2*ETIJ(K,ipass,2,2,0,X) - 2*ETIJ(K,ipass, 1,1,0,X)/TWOA(K)
            ETIJ(K,ipass,0,3,1,X) = PB(K,ipass)*ETIJ(K,ipass, 0,3,0,X)+ ETIJ(K,ipass, 1,3,0,X)
            ETIJ(K,ipass,3,3,1,X) = HPINV(K)*ETIJ(K,ipass,2,3,0,X)+ PB(K,ipass)*ETIJ(K,ipass,3,3,0,X)
            ETIJ(K,ipass,4,3,1,X) = HPINV(K)*ETIJ(K,ipass,3,3,0,X)
            ETIJ(K,ipass,1,3,1,X)= HPINV(K)*ETIJ(K,ipass,0,3,0,X)+ PB(K,ipass)*ETIJ(K,ipass,1,3,0,X)&
                 & + 2*ETIJ(K,ipass,2,3,0,X)! - 0*ETIJ(K,ipass,1,3,*,X)/TWOB(K)
            ETIJ(K,ipass,2,3,1,X)= HPINV(K)*ETIJ(K,ipass,1,3,0,X)+ PB(K,ipass)*ETIJ(K,ipass,2,3,0,X)&
                 & + 3*ETIJ(K,ipass,3,3,0,X)! - 0*ETIJ(K,ipass,2,3,*,X)/TWOB(K)
            ETIJ(K,ipass,0,3,2,X) = PB(K,ipass)*ETIJ(K,ipass, 0,3,1,X)+ ETIJ(K,ipass, 1,3,1,X)&
                 & -1*ETIJ(K,ipass,0,3,0,X)/TWOB(K)
            ETIJ(K,ipass,4,3,2,X) = HPINV(K)*ETIJ(K,ipass,3,3,1,X)+ PB(K,ipass)*ETIJ(K,ipass,4,3,1,X)
            ETIJ(K,ipass,5,3,2,X) = HPINV(K)*ETIJ(K,ipass,4,3,1,X)
            ETIJ(K,ipass,1,3,2,X)= HPINV(K)*ETIJ(K,ipass,0,3,1,X)+ PB(K,ipass)*ETIJ(K,ipass,1,3,1,X)&
                 & + 2*ETIJ(K,ipass,2,3,1,X) - 1*ETIJ(K,ipass,1,3,0,X)/TWOB(K)
            ETIJ(K,ipass,2,3,2,X)= HPINV(K)*ETIJ(K,ipass,1,3,1,X)+ PB(K,ipass)*ETIJ(K,ipass,2,3,1,X)&
                 & + 3*ETIJ(K,ipass,3,3,1,X) - 1*ETIJ(K,ipass,2,3,0,X)/TWOB(K)
            ETIJ(K,ipass,3,3,2,X)= HPINV(K)*ETIJ(K,ipass,2,3,1,X)+ PB(K,ipass)*ETIJ(K,ipass,3,3,1,X)&
                 & + 4*ETIJ(K,ipass,4,3,1,X) - 1*ETIJ(K,ipass,3,3,0,X)/TWOB(K)
            ETIJ(K,ipass,0,4,0,X) = PA(K,ipass)*ETIJ(K,ipass,0,3,0,X)+ ETIJ(K,ipass,1,3,0,X)&
                 & - 3*ETIJ(K,ipass, 0,2,0,X)/TWOA(K)
            ETIJ(K,ipass,3,4,0,X) = HPINV(K)*ETIJ(K,ipass,2,3,0,X)+ PA(K,ipass)*ETIJ(K,ipass,3,3,0,X)
            ETIJ(K,ipass,4,4,0,X) = HPINV(K)*ETIJ(K,ipass,3,3,0,X)
            ETIJ(K,ipass,1,4,0,X) = HPINV(K)*ETIJ(K,ipass,0,3,0,X)+ PA(K,ipass)*ETIJ(K,ipass,1,3,0,X)&
                 & + 2*ETIJ(K,ipass,2,3,0,X) - 3*ETIJ(K,ipass, 1,2,0,X)/TWOA(K)
            ETIJ(K,ipass,2,4,0,X) = HPINV(K)*ETIJ(K,ipass,1,3,0,X)+ PA(K,ipass)*ETIJ(K,ipass,2,3,0,X)&
                 & + 3*ETIJ(K,ipass,3,3,0,X) - 3*ETIJ(K,ipass, 2,2,0,X)/TWOA(K)
            ETIJ(K,ipass,0,4,1,X) = PB(K,ipass)*ETIJ(K,ipass, 0,4,0,X)+ ETIJ(K,ipass, 1,4,0,X)
            ETIJ(K,ipass,4,4,1,X) = HPINV(K)*ETIJ(K,ipass,3,4,0,X)+ PB(K,ipass)*ETIJ(K,ipass,4,4,0,X)
            ETIJ(K,ipass,5,4,1,X) = HPINV(K)*ETIJ(K,ipass,4,4,0,X)
            ETIJ(K,ipass,1,4,1,X)= HPINV(K)*ETIJ(K,ipass,0,4,0,X)+ PB(K,ipass)*ETIJ(K,ipass,1,4,0,X)&
                 & + 2*ETIJ(K,ipass,2,4,0,X)! - 0*ETIJ(K,ipass,1,4,*,X)/TWOB(K)
            ETIJ(K,ipass,2,4,1,X)= HPINV(K)*ETIJ(K,ipass,1,4,0,X)+ PB(K,ipass)*ETIJ(K,ipass,2,4,0,X)&
                 & + 3*ETIJ(K,ipass,3,4,0,X)! - 0*ETIJ(K,ipass,2,4,*,X)/TWOB(K)
            ETIJ(K,ipass,3,4,1,X)= HPINV(K)*ETIJ(K,ipass,2,4,0,X)+ PB(K,ipass)*ETIJ(K,ipass,3,4,0,X)&
                 & + 4*ETIJ(K,ipass,4,4,0,X)! - 0*ETIJ(K,ipass,3,4,*,X)/TWOB(K)
            ETIJ(K,ipass,0,4,2,X) = PB(K,ipass)*ETIJ(K,ipass, 0,4,1,X)+ ETIJ(K,ipass, 1,4,1,X)&
                 & -1*ETIJ(K,ipass,0,4,0,X)/TWOB(K)
            ETIJ(K,ipass,5,4,2,X) = HPINV(K)*ETIJ(K,ipass,4,4,1,X)+ PB(K,ipass)*ETIJ(K,ipass,5,4,1,X)
            ETIJ(K,ipass,6,4,2,X) = HPINV(K)*ETIJ(K,ipass,5,4,1,X)
            ETIJ(K,ipass,1,4,2,X)= HPINV(K)*ETIJ(K,ipass,0,4,1,X)+ PB(K,ipass)*ETIJ(K,ipass,1,4,1,X)&
                 & + 2*ETIJ(K,ipass,2,4,1,X) - 1*ETIJ(K,ipass,1,4,0,X)/TWOB(K)
            ETIJ(K,ipass,2,4,2,X)= HPINV(K)*ETIJ(K,ipass,1,4,1,X)+ PB(K,ipass)*ETIJ(K,ipass,2,4,1,X)&
                 & + 3*ETIJ(K,ipass,3,4,1,X) - 1*ETIJ(K,ipass,2,4,0,X)/TWOB(K)
            ETIJ(K,ipass,3,4,2,X)= HPINV(K)*ETIJ(K,ipass,2,4,1,X)+ PB(K,ipass)*ETIJ(K,ipass,3,4,1,X)&
                 & + 4*ETIJ(K,ipass,4,4,1,X) - 1*ETIJ(K,ipass,3,4,0,X)/TWOB(K)
            ETIJ(K,ipass,4,4,2,X)= HPINV(K)*ETIJ(K,ipass,3,4,1,X)+ PB(K,ipass)*ETIJ(K,ipass,4,4,1,X)&
                 & + 5*ETIJ(K,ipass,5,4,1,X) - 1*ETIJ(K,ipass,4,4,0,X)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.3)THEN
      !maxI.EQ.3.AND.maxJ.EQ.3
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,0,1,X) = PB(K,ipass)    
            ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
            ETIJ(K,ipass,0,0,2,X) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,0,3,X) = PB(K,ipass)*ETIJ(K,ipass, 0,0,2,X)+ ETIJ(K,ipass, 1,0,2,X)&
                 & -2*ETIJ(K,ipass,0,0,1,X)/TWOB(K)
            ETIJ(K,ipass,2,0,3,X) = HPINV(K)*ETIJ(K,ipass,1,0,2,X)+ PB(K,ipass)*ETIJ(K,ipass,2,0,2,X)
            ETIJ(K,ipass,3,0,3,X) = HPINV(K)*ETIJ(K,ipass,2,0,2,X)
            ETIJ(K,ipass,1,0,3,X)= HPINV(K)*ETIJ(K,ipass,0,0,2,X)+ PB(K,ipass)*ETIJ(K,ipass,1,0,2,X)&
                 & + 2*ETIJ(K,ipass,2,0,2,X) - 2*ETIJ(K,ipass,1,0,1,X)/TWOB(K)
            ETIJ(K,ipass,0,1,0,X) = PA(K,ipass)
            ETIJ(K,ipass,1,1,0,X) = HPINV(K)
            ETIJ(K,ipass,0,1,1,X) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,1,2,X) = PB(K,ipass)*ETIJ(K,ipass, 0,1,1,X)+ ETIJ(K,ipass, 1,1,1,X)&
                 & -1*ETIJ(K,ipass,0,1,0,X)/TWOB(K)
            ETIJ(K,ipass,2,1,2,X) = HPINV(K)*ETIJ(K,ipass,1,1,1,X)+ PB(K,ipass)*ETIJ(K,ipass,2,1,1,X)
            ETIJ(K,ipass,3,1,2,X) = HPINV(K)*ETIJ(K,ipass,2,1,1,X)
            ETIJ(K,ipass,1,1,2,X)= HPINV(K)*ETIJ(K,ipass,0,1,1,X)+ PB(K,ipass)*ETIJ(K,ipass,1,1,1,X)&
                 & + 2*ETIJ(K,ipass,2,1,1,X) - 1*ETIJ(K,ipass,1,1,0,X)/TWOB(K)
            ETIJ(K,ipass,0,1,3,X) = PB(K,ipass)*ETIJ(K,ipass, 0,1,2,X)+ ETIJ(K,ipass, 1,1,2,X)&
                 & -2*ETIJ(K,ipass,0,1,1,X)/TWOB(K)
            ETIJ(K,ipass,3,1,3,X) = HPINV(K)*ETIJ(K,ipass,2,1,2,X)+ PB(K,ipass)*ETIJ(K,ipass,3,1,2,X)
            ETIJ(K,ipass,4,1,3,X) = HPINV(K)*ETIJ(K,ipass,3,1,2,X)
            ETIJ(K,ipass,1,1,3,X)= HPINV(K)*ETIJ(K,ipass,0,1,2,X)+ PB(K,ipass)*ETIJ(K,ipass,1,1,2,X)&
                 & + 2*ETIJ(K,ipass,2,1,2,X) - 2*ETIJ(K,ipass,1,1,1,X)/TWOB(K)
            ETIJ(K,ipass,2,1,3,X)= HPINV(K)*ETIJ(K,ipass,1,1,2,X)+ PB(K,ipass)*ETIJ(K,ipass,2,1,2,X)&
                 & + 3*ETIJ(K,ipass,3,1,2,X) - 2*ETIJ(K,ipass,2,1,1,X)/TWOB(K)
            ETIJ(K,ipass,0,2,0,X) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,2,1,X) = PB(K,ipass)*ETIJ(K,ipass, 0,2,0,X)+ ETIJ(K,ipass, 1,2,0,X)
            ETIJ(K,ipass,2,2,1,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PB(K,ipass)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,3,2,1,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,1,2,1,X)= HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PB(K,ipass)*ETIJ(K,ipass,1,2,0,X)&
                 & + 2*ETIJ(K,ipass,2,2,0,X)! - 0*ETIJ(K,ipass,1,2,*,X)/TWOB(K)
            ETIJ(K,ipass,0,2,2,X) = PB(K,ipass)*ETIJ(K,ipass, 0,2,1,X)+ ETIJ(K,ipass, 1,2,1,X)&
                 & -1*ETIJ(K,ipass,0,2,0,X)/TWOB(K)
            ETIJ(K,ipass,3,2,2,X) = HPINV(K)*ETIJ(K,ipass,2,2,1,X)+ PB(K,ipass)*ETIJ(K,ipass,3,2,1,X)
            ETIJ(K,ipass,4,2,2,X) = HPINV(K)*ETIJ(K,ipass,3,2,1,X)
            ETIJ(K,ipass,1,2,2,X)= HPINV(K)*ETIJ(K,ipass,0,2,1,X)+ PB(K,ipass)*ETIJ(K,ipass,1,2,1,X)&
                 & + 2*ETIJ(K,ipass,2,2,1,X) - 1*ETIJ(K,ipass,1,2,0,X)/TWOB(K)
            ETIJ(K,ipass,2,2,2,X)= HPINV(K)*ETIJ(K,ipass,1,2,1,X)+ PB(K,ipass)*ETIJ(K,ipass,2,2,1,X)&
                 & + 3*ETIJ(K,ipass,3,2,1,X) - 1*ETIJ(K,ipass,2,2,0,X)/TWOB(K)
            ETIJ(K,ipass,0,2,3,X) = PB(K,ipass)*ETIJ(K,ipass, 0,2,2,X)+ ETIJ(K,ipass, 1,2,2,X)&
                 & -2*ETIJ(K,ipass,0,2,1,X)/TWOB(K)
            ETIJ(K,ipass,4,2,3,X) = HPINV(K)*ETIJ(K,ipass,3,2,2,X)+ PB(K,ipass)*ETIJ(K,ipass,4,2,2,X)
            ETIJ(K,ipass,5,2,3,X) = HPINV(K)*ETIJ(K,ipass,4,2,2,X)
            ETIJ(K,ipass,1,2,3,X)= HPINV(K)*ETIJ(K,ipass,0,2,2,X)+ PB(K,ipass)*ETIJ(K,ipass,1,2,2,X)&
                 & + 2*ETIJ(K,ipass,2,2,2,X) - 2*ETIJ(K,ipass,1,2,1,X)/TWOB(K)
            ETIJ(K,ipass,2,2,3,X)= HPINV(K)*ETIJ(K,ipass,1,2,2,X)+ PB(K,ipass)*ETIJ(K,ipass,2,2,2,X)&
                 & + 3*ETIJ(K,ipass,3,2,2,X) - 2*ETIJ(K,ipass,2,2,1,X)/TWOB(K)
            ETIJ(K,ipass,3,2,3,X)= HPINV(K)*ETIJ(K,ipass,2,2,2,X)+ PB(K,ipass)*ETIJ(K,ipass,3,2,2,X)&
                 & + 4*ETIJ(K,ipass,4,2,2,X) - 2*ETIJ(K,ipass,3,2,1,X)/TWOB(K)
            ETIJ(K,ipass,0,3,0,X) = PA(K,ipass)*ETIJ(K,ipass,0,2,0,X)+ ETIJ(K,ipass,1,2,0,X)&
                 & - 2*ETIJ(K,ipass, 0,1,0,X)/TWOA(K)
            ETIJ(K,ipass,2,3,0,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PA(K,ipass)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,3,3,0,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,1,3,0,X) = HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PA(K,ipass)*ETIJ(K,ipass,1,2,0,X)&
                 & + 2*ETIJ(K,ipass,2,2,0,X) - 2*ETIJ(K,ipass, 1,1,0,X)/TWOA(K)
            ETIJ(K,ipass,0,3,1,X) = PB(K,ipass)*ETIJ(K,ipass, 0,3,0,X)+ ETIJ(K,ipass, 1,3,0,X)
            ETIJ(K,ipass,3,3,1,X) = HPINV(K)*ETIJ(K,ipass,2,3,0,X)+ PB(K,ipass)*ETIJ(K,ipass,3,3,0,X)
            ETIJ(K,ipass,4,3,1,X) = HPINV(K)*ETIJ(K,ipass,3,3,0,X)
            ETIJ(K,ipass,1,3,1,X)= HPINV(K)*ETIJ(K,ipass,0,3,0,X)+ PB(K,ipass)*ETIJ(K,ipass,1,3,0,X)&
                 & + 2*ETIJ(K,ipass,2,3,0,X)! - 0*ETIJ(K,ipass,1,3,*,X)/TWOB(K)
            ETIJ(K,ipass,2,3,1,X)= HPINV(K)*ETIJ(K,ipass,1,3,0,X)+ PB(K,ipass)*ETIJ(K,ipass,2,3,0,X)&
                 & + 3*ETIJ(K,ipass,3,3,0,X)! - 0*ETIJ(K,ipass,2,3,*,X)/TWOB(K)
            ETIJ(K,ipass,0,3,2,X) = PB(K,ipass)*ETIJ(K,ipass, 0,3,1,X)+ ETIJ(K,ipass, 1,3,1,X)&
                 & -1*ETIJ(K,ipass,0,3,0,X)/TWOB(K)
            ETIJ(K,ipass,4,3,2,X) = HPINV(K)*ETIJ(K,ipass,3,3,1,X)+ PB(K,ipass)*ETIJ(K,ipass,4,3,1,X)
            ETIJ(K,ipass,5,3,2,X) = HPINV(K)*ETIJ(K,ipass,4,3,1,X)
            ETIJ(K,ipass,1,3,2,X)= HPINV(K)*ETIJ(K,ipass,0,3,1,X)+ PB(K,ipass)*ETIJ(K,ipass,1,3,1,X)&
                 & + 2*ETIJ(K,ipass,2,3,1,X) - 1*ETIJ(K,ipass,1,3,0,X)/TWOB(K)
            ETIJ(K,ipass,2,3,2,X)= HPINV(K)*ETIJ(K,ipass,1,3,1,X)+ PB(K,ipass)*ETIJ(K,ipass,2,3,1,X)&
                 & + 3*ETIJ(K,ipass,3,3,1,X) - 1*ETIJ(K,ipass,2,3,0,X)/TWOB(K)
            ETIJ(K,ipass,3,3,2,X)= HPINV(K)*ETIJ(K,ipass,2,3,1,X)+ PB(K,ipass)*ETIJ(K,ipass,3,3,1,X)&
                 & + 4*ETIJ(K,ipass,4,3,1,X) - 1*ETIJ(K,ipass,3,3,0,X)/TWOB(K)
            ETIJ(K,ipass,0,3,3,X) = PB(K,ipass)*ETIJ(K,ipass, 0,3,2,X)+ ETIJ(K,ipass, 1,3,2,X)&
                 & -2*ETIJ(K,ipass,0,3,1,X)/TWOB(K)
            ETIJ(K,ipass,5,3,3,X) = HPINV(K)*ETIJ(K,ipass,4,3,2,X)+ PB(K,ipass)*ETIJ(K,ipass,5,3,2,X)
            ETIJ(K,ipass,6,3,3,X) = HPINV(K)*ETIJ(K,ipass,5,3,2,X)
            ETIJ(K,ipass,1,3,3,X)= HPINV(K)*ETIJ(K,ipass,0,3,2,X)+ PB(K,ipass)*ETIJ(K,ipass,1,3,2,X)&
                 & + 2*ETIJ(K,ipass,2,3,2,X) - 2*ETIJ(K,ipass,1,3,1,X)/TWOB(K)
            ETIJ(K,ipass,2,3,3,X)= HPINV(K)*ETIJ(K,ipass,1,3,2,X)+ PB(K,ipass)*ETIJ(K,ipass,2,3,2,X)&
                 & + 3*ETIJ(K,ipass,3,3,2,X) - 2*ETIJ(K,ipass,2,3,1,X)/TWOB(K)
            ETIJ(K,ipass,3,3,3,X)= HPINV(K)*ETIJ(K,ipass,2,3,2,X)+ PB(K,ipass)*ETIJ(K,ipass,3,3,2,X)&
                 & + 4*ETIJ(K,ipass,4,3,2,X) - 2*ETIJ(K,ipass,3,3,1,X)/TWOB(K)
            ETIJ(K,ipass,4,3,3,X)= HPINV(K)*ETIJ(K,ipass,3,3,2,X)+ PB(K,ipass)*ETIJ(K,ipass,4,3,2,X)&
                 & + 5*ETIJ(K,ipass,5,3,2,X) - 2*ETIJ(K,ipass,4,3,1,X)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.2)THEN
      !maxI.EQ.2.AND.maxJ.EQ.4
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,0,1,X) = PB(K,ipass)    
            ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
            ETIJ(K,ipass,0,0,2,X) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,0,3,X) = PB(K,ipass)*ETIJ(K,ipass, 0,0,2,X)+ ETIJ(K,ipass, 1,0,2,X)&
                 & -2*ETIJ(K,ipass,0,0,1,X)/TWOB(K)
            ETIJ(K,ipass,2,0,3,X) = HPINV(K)*ETIJ(K,ipass,1,0,2,X)+ PB(K,ipass)*ETIJ(K,ipass,2,0,2,X)
            ETIJ(K,ipass,3,0,3,X) = HPINV(K)*ETIJ(K,ipass,2,0,2,X)
            ETIJ(K,ipass,1,0,3,X)= HPINV(K)*ETIJ(K,ipass,0,0,2,X)+ PB(K,ipass)*ETIJ(K,ipass,1,0,2,X)&
                 & + 2*ETIJ(K,ipass,2,0,2,X) - 2*ETIJ(K,ipass,1,0,1,X)/TWOB(K)
            ETIJ(K,ipass,0,0,4,X) = PB(K,ipass)*ETIJ(K,ipass, 0,0,3,X)+ ETIJ(K,ipass, 1,0,3,X)&
                 & -3*ETIJ(K,ipass,0,0,2,X)/TWOB(K)
            ETIJ(K,ipass,3,0,4,X) = HPINV(K)*ETIJ(K,ipass,2,0,3,X)+ PB(K,ipass)*ETIJ(K,ipass,3,0,3,X)
            ETIJ(K,ipass,4,0,4,X) = HPINV(K)*ETIJ(K,ipass,3,0,3,X)
            ETIJ(K,ipass,1,0,4,X)= HPINV(K)*ETIJ(K,ipass,0,0,3,X)+ PB(K,ipass)*ETIJ(K,ipass,1,0,3,X)&
                 & + 2*ETIJ(K,ipass,2,0,3,X) - 3*ETIJ(K,ipass,1,0,2,X)/TWOB(K)
            ETIJ(K,ipass,2,0,4,X)= HPINV(K)*ETIJ(K,ipass,1,0,3,X)+ PB(K,ipass)*ETIJ(K,ipass,2,0,3,X)&
                 & + 3*ETIJ(K,ipass,3,0,3,X) - 3*ETIJ(K,ipass,2,0,2,X)/TWOB(K)
            ETIJ(K,ipass,0,1,0,X) = PA(K,ipass)
            ETIJ(K,ipass,1,1,0,X) = HPINV(K)
            ETIJ(K,ipass,0,1,1,X) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,1,2,X) = PB(K,ipass)*ETIJ(K,ipass, 0,1,1,X)+ ETIJ(K,ipass, 1,1,1,X)&
                 & -1*ETIJ(K,ipass,0,1,0,X)/TWOB(K)
            ETIJ(K,ipass,2,1,2,X) = HPINV(K)*ETIJ(K,ipass,1,1,1,X)+ PB(K,ipass)*ETIJ(K,ipass,2,1,1,X)
            ETIJ(K,ipass,3,1,2,X) = HPINV(K)*ETIJ(K,ipass,2,1,1,X)
            ETIJ(K,ipass,1,1,2,X)= HPINV(K)*ETIJ(K,ipass,0,1,1,X)+ PB(K,ipass)*ETIJ(K,ipass,1,1,1,X)&
                 & + 2*ETIJ(K,ipass,2,1,1,X) - 1*ETIJ(K,ipass,1,1,0,X)/TWOB(K)
            ETIJ(K,ipass,0,1,3,X) = PB(K,ipass)*ETIJ(K,ipass, 0,1,2,X)+ ETIJ(K,ipass, 1,1,2,X)&
                 & -2*ETIJ(K,ipass,0,1,1,X)/TWOB(K)
            ETIJ(K,ipass,3,1,3,X) = HPINV(K)*ETIJ(K,ipass,2,1,2,X)+ PB(K,ipass)*ETIJ(K,ipass,3,1,2,X)
            ETIJ(K,ipass,4,1,3,X) = HPINV(K)*ETIJ(K,ipass,3,1,2,X)
            ETIJ(K,ipass,1,1,3,X)= HPINV(K)*ETIJ(K,ipass,0,1,2,X)+ PB(K,ipass)*ETIJ(K,ipass,1,1,2,X)&
                 & + 2*ETIJ(K,ipass,2,1,2,X) - 2*ETIJ(K,ipass,1,1,1,X)/TWOB(K)
            ETIJ(K,ipass,2,1,3,X)= HPINV(K)*ETIJ(K,ipass,1,1,2,X)+ PB(K,ipass)*ETIJ(K,ipass,2,1,2,X)&
                 & + 3*ETIJ(K,ipass,3,1,2,X) - 2*ETIJ(K,ipass,2,1,1,X)/TWOB(K)
            ETIJ(K,ipass,0,1,4,X) = PB(K,ipass)*ETIJ(K,ipass, 0,1,3,X)+ ETIJ(K,ipass, 1,1,3,X)&
                 & -3*ETIJ(K,ipass,0,1,2,X)/TWOB(K)
            ETIJ(K,ipass,4,1,4,X) = HPINV(K)*ETIJ(K,ipass,3,1,3,X)+ PB(K,ipass)*ETIJ(K,ipass,4,1,3,X)
            ETIJ(K,ipass,5,1,4,X) = HPINV(K)*ETIJ(K,ipass,4,1,3,X)
            ETIJ(K,ipass,1,1,4,X)= HPINV(K)*ETIJ(K,ipass,0,1,3,X)+ PB(K,ipass)*ETIJ(K,ipass,1,1,3,X)&
                 & + 2*ETIJ(K,ipass,2,1,3,X) - 3*ETIJ(K,ipass,1,1,2,X)/TWOB(K)
            ETIJ(K,ipass,2,1,4,X)= HPINV(K)*ETIJ(K,ipass,1,1,3,X)+ PB(K,ipass)*ETIJ(K,ipass,2,1,3,X)&
                 & + 3*ETIJ(K,ipass,3,1,3,X) - 3*ETIJ(K,ipass,2,1,2,X)/TWOB(K)
            ETIJ(K,ipass,3,1,4,X)= HPINV(K)*ETIJ(K,ipass,2,1,3,X)+ PB(K,ipass)*ETIJ(K,ipass,3,1,3,X)&
                 & + 4*ETIJ(K,ipass,4,1,3,X) - 3*ETIJ(K,ipass,3,1,2,X)/TWOB(K)
            ETIJ(K,ipass,0,2,0,X) = PA(K,ipass)*PA(K,ipass) + HPINV(K) - D1/TWOA(K)
            ETIJ(K,ipass,1,2,0,X) = D2*PA(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,2,0,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,2,1,X) = PB(K,ipass)*ETIJ(K,ipass, 0,2,0,X)+ ETIJ(K,ipass, 1,2,0,X)
            ETIJ(K,ipass,2,2,1,X) = HPINV(K)*ETIJ(K,ipass,1,2,0,X)+ PB(K,ipass)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,3,2,1,X) = HPINV(K)*ETIJ(K,ipass,2,2,0,X)
            ETIJ(K,ipass,1,2,1,X)= HPINV(K)*ETIJ(K,ipass,0,2,0,X)+ PB(K,ipass)*ETIJ(K,ipass,1,2,0,X)&
                 & + 2*ETIJ(K,ipass,2,2,0,X)! - 0*ETIJ(K,ipass,1,2,*,X)/TWOB(K)
            ETIJ(K,ipass,0,2,2,X) = PB(K,ipass)*ETIJ(K,ipass, 0,2,1,X)+ ETIJ(K,ipass, 1,2,1,X)&
                 & -1*ETIJ(K,ipass,0,2,0,X)/TWOB(K)
            ETIJ(K,ipass,3,2,2,X) = HPINV(K)*ETIJ(K,ipass,2,2,1,X)+ PB(K,ipass)*ETIJ(K,ipass,3,2,1,X)
            ETIJ(K,ipass,4,2,2,X) = HPINV(K)*ETIJ(K,ipass,3,2,1,X)
            ETIJ(K,ipass,1,2,2,X)= HPINV(K)*ETIJ(K,ipass,0,2,1,X)+ PB(K,ipass)*ETIJ(K,ipass,1,2,1,X)&
                 & + 2*ETIJ(K,ipass,2,2,1,X) - 1*ETIJ(K,ipass,1,2,0,X)/TWOB(K)
            ETIJ(K,ipass,2,2,2,X)= HPINV(K)*ETIJ(K,ipass,1,2,1,X)+ PB(K,ipass)*ETIJ(K,ipass,2,2,1,X)&
                 & + 3*ETIJ(K,ipass,3,2,1,X) - 1*ETIJ(K,ipass,2,2,0,X)/TWOB(K)
            ETIJ(K,ipass,0,2,3,X) = PB(K,ipass)*ETIJ(K,ipass, 0,2,2,X)+ ETIJ(K,ipass, 1,2,2,X)&
                 & -2*ETIJ(K,ipass,0,2,1,X)/TWOB(K)
            ETIJ(K,ipass,4,2,3,X) = HPINV(K)*ETIJ(K,ipass,3,2,2,X)+ PB(K,ipass)*ETIJ(K,ipass,4,2,2,X)
            ETIJ(K,ipass,5,2,3,X) = HPINV(K)*ETIJ(K,ipass,4,2,2,X)
            ETIJ(K,ipass,1,2,3,X)= HPINV(K)*ETIJ(K,ipass,0,2,2,X)+ PB(K,ipass)*ETIJ(K,ipass,1,2,2,X)&
                 & + 2*ETIJ(K,ipass,2,2,2,X) - 2*ETIJ(K,ipass,1,2,1,X)/TWOB(K)
            ETIJ(K,ipass,2,2,3,X)= HPINV(K)*ETIJ(K,ipass,1,2,2,X)+ PB(K,ipass)*ETIJ(K,ipass,2,2,2,X)&
                 & + 3*ETIJ(K,ipass,3,2,2,X) - 2*ETIJ(K,ipass,2,2,1,X)/TWOB(K)
            ETIJ(K,ipass,3,2,3,X)= HPINV(K)*ETIJ(K,ipass,2,2,2,X)+ PB(K,ipass)*ETIJ(K,ipass,3,2,2,X)&
                 & + 4*ETIJ(K,ipass,4,2,2,X) - 2*ETIJ(K,ipass,3,2,1,X)/TWOB(K)
            ETIJ(K,ipass,0,2,4,X) = PB(K,ipass)*ETIJ(K,ipass, 0,2,3,X)+ ETIJ(K,ipass, 1,2,3,X)&
                 & -3*ETIJ(K,ipass,0,2,2,X)/TWOB(K)
            ETIJ(K,ipass,5,2,4,X) = HPINV(K)*ETIJ(K,ipass,4,2,3,X)+ PB(K,ipass)*ETIJ(K,ipass,5,2,3,X)
            ETIJ(K,ipass,6,2,4,X) = HPINV(K)*ETIJ(K,ipass,5,2,3,X)
            ETIJ(K,ipass,1,2,4,X)= HPINV(K)*ETIJ(K,ipass,0,2,3,X)+ PB(K,ipass)*ETIJ(K,ipass,1,2,3,X)&
                 & + 2*ETIJ(K,ipass,2,2,3,X) - 3*ETIJ(K,ipass,1,2,2,X)/TWOB(K)
            ETIJ(K,ipass,2,2,4,X)= HPINV(K)*ETIJ(K,ipass,1,2,3,X)+ PB(K,ipass)*ETIJ(K,ipass,2,2,3,X)&
                 & + 3*ETIJ(K,ipass,3,2,3,X) - 3*ETIJ(K,ipass,2,2,2,X)/TWOB(K)
            ETIJ(K,ipass,3,2,4,X)= HPINV(K)*ETIJ(K,ipass,2,2,3,X)+ PB(K,ipass)*ETIJ(K,ipass,3,2,3,X)&
                 & + 4*ETIJ(K,ipass,4,2,3,X) - 3*ETIJ(K,ipass,3,2,2,X)/TWOB(K)
            ETIJ(K,ipass,4,2,4,X)= HPINV(K)*ETIJ(K,ipass,3,2,3,X)+ PB(K,ipass)*ETIJ(K,ipass,4,2,3,X)&
                 & + 5*ETIJ(K,ipass,5,2,3,X) - 3*ETIJ(K,ipass,4,2,2,X)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.1)THEN
      !maxI.EQ.1.AND.maxJ.EQ.5
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,0,1,X) = PB(K,ipass)    
            ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
            ETIJ(K,ipass,0,0,2,X) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,0,3,X) = PB(K,ipass)*ETIJ(K,ipass, 0,0,2,X)+ ETIJ(K,ipass, 1,0,2,X)&
                 & -2*ETIJ(K,ipass,0,0,1,X)/TWOB(K)
            ETIJ(K,ipass,2,0,3,X) = HPINV(K)*ETIJ(K,ipass,1,0,2,X)+ PB(K,ipass)*ETIJ(K,ipass,2,0,2,X)
            ETIJ(K,ipass,3,0,3,X) = HPINV(K)*ETIJ(K,ipass,2,0,2,X)
            ETIJ(K,ipass,1,0,3,X)= HPINV(K)*ETIJ(K,ipass,0,0,2,X)+ PB(K,ipass)*ETIJ(K,ipass,1,0,2,X)&
                 & + 2*ETIJ(K,ipass,2,0,2,X) - 2*ETIJ(K,ipass,1,0,1,X)/TWOB(K)
            ETIJ(K,ipass,0,0,4,X) = PB(K,ipass)*ETIJ(K,ipass, 0,0,3,X)+ ETIJ(K,ipass, 1,0,3,X)&
                 & -3*ETIJ(K,ipass,0,0,2,X)/TWOB(K)
            ETIJ(K,ipass,3,0,4,X) = HPINV(K)*ETIJ(K,ipass,2,0,3,X)+ PB(K,ipass)*ETIJ(K,ipass,3,0,3,X)
            ETIJ(K,ipass,4,0,4,X) = HPINV(K)*ETIJ(K,ipass,3,0,3,X)
            ETIJ(K,ipass,1,0,4,X)= HPINV(K)*ETIJ(K,ipass,0,0,3,X)+ PB(K,ipass)*ETIJ(K,ipass,1,0,3,X)&
                 & + 2*ETIJ(K,ipass,2,0,3,X) - 3*ETIJ(K,ipass,1,0,2,X)/TWOB(K)
            ETIJ(K,ipass,2,0,4,X)= HPINV(K)*ETIJ(K,ipass,1,0,3,X)+ PB(K,ipass)*ETIJ(K,ipass,2,0,3,X)&
                 & + 3*ETIJ(K,ipass,3,0,3,X) - 3*ETIJ(K,ipass,2,0,2,X)/TWOB(K)
            ETIJ(K,ipass,0,0,5,X) = PB(K,ipass)*ETIJ(K,ipass, 0,0,4,X)+ ETIJ(K,ipass, 1,0,4,X)&
                 & -4*ETIJ(K,ipass,0,0,3,X)/TWOB(K)
            ETIJ(K,ipass,4,0,5,X) = HPINV(K)*ETIJ(K,ipass,3,0,4,X)+ PB(K,ipass)*ETIJ(K,ipass,4,0,4,X)
            ETIJ(K,ipass,5,0,5,X) = HPINV(K)*ETIJ(K,ipass,4,0,4,X)
            ETIJ(K,ipass,1,0,5,X)= HPINV(K)*ETIJ(K,ipass,0,0,4,X)+ PB(K,ipass)*ETIJ(K,ipass,1,0,4,X)&
                 & + 2*ETIJ(K,ipass,2,0,4,X) - 4*ETIJ(K,ipass,1,0,3,X)/TWOB(K)
            ETIJ(K,ipass,2,0,5,X)= HPINV(K)*ETIJ(K,ipass,1,0,4,X)+ PB(K,ipass)*ETIJ(K,ipass,2,0,4,X)&
                 & + 3*ETIJ(K,ipass,3,0,4,X) - 4*ETIJ(K,ipass,2,0,3,X)/TWOB(K)
            ETIJ(K,ipass,3,0,5,X)= HPINV(K)*ETIJ(K,ipass,2,0,4,X)+ PB(K,ipass)*ETIJ(K,ipass,3,0,4,X)&
                 & + 4*ETIJ(K,ipass,4,0,4,X) - 4*ETIJ(K,ipass,3,0,3,X)/TWOB(K)
            ETIJ(K,ipass,0,1,0,X) = PA(K,ipass)
            ETIJ(K,ipass,1,1,0,X) = HPINV(K)
            ETIJ(K,ipass,0,1,1,X) = PA(K,ipass)*PB(K,ipass) + HPINV(K)
            ETIJ(K,ipass,1,1,1,X) = (PA(K,ipass) + PB(K,ipass))*HPINV(K)
            ETIJ(K,ipass,2,1,1,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,1,2,X) = PB(K,ipass)*ETIJ(K,ipass, 0,1,1,X)+ ETIJ(K,ipass, 1,1,1,X)&
                 & -1*ETIJ(K,ipass,0,1,0,X)/TWOB(K)
            ETIJ(K,ipass,2,1,2,X) = HPINV(K)*ETIJ(K,ipass,1,1,1,X)+ PB(K,ipass)*ETIJ(K,ipass,2,1,1,X)
            ETIJ(K,ipass,3,1,2,X) = HPINV(K)*ETIJ(K,ipass,2,1,1,X)
            ETIJ(K,ipass,1,1,2,X)= HPINV(K)*ETIJ(K,ipass,0,1,1,X)+ PB(K,ipass)*ETIJ(K,ipass,1,1,1,X)&
                 & + 2*ETIJ(K,ipass,2,1,1,X) - 1*ETIJ(K,ipass,1,1,0,X)/TWOB(K)
            ETIJ(K,ipass,0,1,3,X) = PB(K,ipass)*ETIJ(K,ipass, 0,1,2,X)+ ETIJ(K,ipass, 1,1,2,X)&
                 & -2*ETIJ(K,ipass,0,1,1,X)/TWOB(K)
            ETIJ(K,ipass,3,1,3,X) = HPINV(K)*ETIJ(K,ipass,2,1,2,X)+ PB(K,ipass)*ETIJ(K,ipass,3,1,2,X)
            ETIJ(K,ipass,4,1,3,X) = HPINV(K)*ETIJ(K,ipass,3,1,2,X)
            ETIJ(K,ipass,1,1,3,X)= HPINV(K)*ETIJ(K,ipass,0,1,2,X)+ PB(K,ipass)*ETIJ(K,ipass,1,1,2,X)&
                 & + 2*ETIJ(K,ipass,2,1,2,X) - 2*ETIJ(K,ipass,1,1,1,X)/TWOB(K)
            ETIJ(K,ipass,2,1,3,X)= HPINV(K)*ETIJ(K,ipass,1,1,2,X)+ PB(K,ipass)*ETIJ(K,ipass,2,1,2,X)&
                 & + 3*ETIJ(K,ipass,3,1,2,X) - 2*ETIJ(K,ipass,2,1,1,X)/TWOB(K)
            ETIJ(K,ipass,0,1,4,X) = PB(K,ipass)*ETIJ(K,ipass, 0,1,3,X)+ ETIJ(K,ipass, 1,1,3,X)&
                 & -3*ETIJ(K,ipass,0,1,2,X)/TWOB(K)
            ETIJ(K,ipass,4,1,4,X) = HPINV(K)*ETIJ(K,ipass,3,1,3,X)+ PB(K,ipass)*ETIJ(K,ipass,4,1,3,X)
            ETIJ(K,ipass,5,1,4,X) = HPINV(K)*ETIJ(K,ipass,4,1,3,X)
            ETIJ(K,ipass,1,1,4,X)= HPINV(K)*ETIJ(K,ipass,0,1,3,X)+ PB(K,ipass)*ETIJ(K,ipass,1,1,3,X)&
                 & + 2*ETIJ(K,ipass,2,1,3,X) - 3*ETIJ(K,ipass,1,1,2,X)/TWOB(K)
            ETIJ(K,ipass,2,1,4,X)= HPINV(K)*ETIJ(K,ipass,1,1,3,X)+ PB(K,ipass)*ETIJ(K,ipass,2,1,3,X)&
                 & + 3*ETIJ(K,ipass,3,1,3,X) - 3*ETIJ(K,ipass,2,1,2,X)/TWOB(K)
            ETIJ(K,ipass,3,1,4,X)= HPINV(K)*ETIJ(K,ipass,2,1,3,X)+ PB(K,ipass)*ETIJ(K,ipass,3,1,3,X)&
                 & + 4*ETIJ(K,ipass,4,1,3,X) - 3*ETIJ(K,ipass,3,1,2,X)/TWOB(K)
            ETIJ(K,ipass,0,1,5,X) = PB(K,ipass)*ETIJ(K,ipass, 0,1,4,X)+ ETIJ(K,ipass, 1,1,4,X)&
                 & -4*ETIJ(K,ipass,0,1,3,X)/TWOB(K)
            ETIJ(K,ipass,5,1,5,X) = HPINV(K)*ETIJ(K,ipass,4,1,4,X)+ PB(K,ipass)*ETIJ(K,ipass,5,1,4,X)
            ETIJ(K,ipass,6,1,5,X) = HPINV(K)*ETIJ(K,ipass,5,1,4,X)
            ETIJ(K,ipass,1,1,5,X)= HPINV(K)*ETIJ(K,ipass,0,1,4,X)+ PB(K,ipass)*ETIJ(K,ipass,1,1,4,X)&
                 & + 2*ETIJ(K,ipass,2,1,4,X) - 4*ETIJ(K,ipass,1,1,3,X)/TWOB(K)
            ETIJ(K,ipass,2,1,5,X)= HPINV(K)*ETIJ(K,ipass,1,1,4,X)+ PB(K,ipass)*ETIJ(K,ipass,2,1,4,X)&
                 & + 3*ETIJ(K,ipass,3,1,4,X) - 4*ETIJ(K,ipass,2,1,3,X)/TWOB(K)
            ETIJ(K,ipass,3,1,5,X)= HPINV(K)*ETIJ(K,ipass,2,1,4,X)+ PB(K,ipass)*ETIJ(K,ipass,3,1,4,X)&
                 & + 4*ETIJ(K,ipass,4,1,4,X) - 4*ETIJ(K,ipass,3,1,3,X)/TWOB(K)
            ETIJ(K,ipass,4,1,5,X)= HPINV(K)*ETIJ(K,ipass,3,1,4,X)+ PB(K,ipass)*ETIJ(K,ipass,4,1,4,X)&
                 & + 5*ETIJ(K,ipass,5,1,4,X) - 4*ETIJ(K,ipass,4,1,3,X)/TWOB(K)
         ENDDO
      ENDDO
   ELSEIF(maxI.EQ.0)THEN
      !maxI.EQ.0.AND.maxJ.EQ.6
      DO Ipass = 1,nPass
         DO K = 1, nPrimP
            ETIJ(K,ipass,0,0,0,X) = D1
            ETIJ(K,ipass,0,0,1,X) = PB(K,ipass)    
            ETIJ(K,ipass,1,0,1,X) = HPINV(K) 
            ETIJ(K,ipass,0,0,2,X) = PB(K,ipass)*PB(K,ipass) + HPINV(K) - D1/TWOB(K)
            ETIJ(K,ipass,1,0,2,X) = D2*PB(K,ipass)*HPINV(K)
            ETIJ(K,ipass,2,0,2,X) = HPINV(K)*HPINV(K)
            ETIJ(K,ipass,0,0,3,X) = PB(K,ipass)*ETIJ(K,ipass, 0,0,2,X)+ ETIJ(K,ipass, 1,0,2,X)&
                 & -2*ETIJ(K,ipass,0,0,1,X)/TWOB(K)
            ETIJ(K,ipass,2,0,3,X) = HPINV(K)*ETIJ(K,ipass,1,0,2,X)+ PB(K,ipass)*ETIJ(K,ipass,2,0,2,X)
            ETIJ(K,ipass,3,0,3,X) = HPINV(K)*ETIJ(K,ipass,2,0,2,X)
            ETIJ(K,ipass,1,0,3,X)= HPINV(K)*ETIJ(K,ipass,0,0,2,X)+ PB(K,ipass)*ETIJ(K,ipass,1,0,2,X)&
                 & + 2*ETIJ(K,ipass,2,0,2,X) - 2*ETIJ(K,ipass,1,0,1,X)/TWOB(K)
            ETIJ(K,ipass,0,0,4,X) = PB(K,ipass)*ETIJ(K,ipass, 0,0,3,X)+ ETIJ(K,ipass, 1,0,3,X)&
                 & -3*ETIJ(K,ipass,0,0,2,X)/TWOB(K)
            ETIJ(K,ipass,3,0,4,X) = HPINV(K)*ETIJ(K,ipass,2,0,3,X)+ PB(K,ipass)*ETIJ(K,ipass,3,0,3,X)
            ETIJ(K,ipass,4,0,4,X) = HPINV(K)*ETIJ(K,ipass,3,0,3,X)
            ETIJ(K,ipass,1,0,4,X)= HPINV(K)*ETIJ(K,ipass,0,0,3,X)+ PB(K,ipass)*ETIJ(K,ipass,1,0,3,X)&
                 & + 2*ETIJ(K,ipass,2,0,3,X) - 3*ETIJ(K,ipass,1,0,2,X)/TWOB(K)
            ETIJ(K,ipass,2,0,4,X)= HPINV(K)*ETIJ(K,ipass,1,0,3,X)+ PB(K,ipass)*ETIJ(K,ipass,2,0,3,X)&
                 & + 3*ETIJ(K,ipass,3,0,3,X) - 3*ETIJ(K,ipass,2,0,2,X)/TWOB(K)
            ETIJ(K,ipass,0,0,5,X) = PB(K,ipass)*ETIJ(K,ipass, 0,0,4,X)+ ETIJ(K,ipass, 1,0,4,X)&
                 & -4*ETIJ(K,ipass,0,0,3,X)/TWOB(K)
            ETIJ(K,ipass,4,0,5,X) = HPINV(K)*ETIJ(K,ipass,3,0,4,X)+ PB(K,ipass)*ETIJ(K,ipass,4,0,4,X)
            ETIJ(K,ipass,5,0,5,X) = HPINV(K)*ETIJ(K,ipass,4,0,4,X)
            ETIJ(K,ipass,1,0,5,X)= HPINV(K)*ETIJ(K,ipass,0,0,4,X)+ PB(K,ipass)*ETIJ(K,ipass,1,0,4,X)&
                 & + 2*ETIJ(K,ipass,2,0,4,X) - 4*ETIJ(K,ipass,1,0,3,X)/TWOB(K)
            ETIJ(K,ipass,2,0,5,X)= HPINV(K)*ETIJ(K,ipass,1,0,4,X)+ PB(K,ipass)*ETIJ(K,ipass,2,0,4,X)&
                 & + 3*ETIJ(K,ipass,3,0,4,X) - 4*ETIJ(K,ipass,2,0,3,X)/TWOB(K)
            ETIJ(K,ipass,3,0,5,X)= HPINV(K)*ETIJ(K,ipass,2,0,4,X)+ PB(K,ipass)*ETIJ(K,ipass,3,0,4,X)&
                 & + 4*ETIJ(K,ipass,4,0,4,X) - 4*ETIJ(K,ipass,3,0,3,X)/TWOB(K)
            ETIJ(K,ipass,0,0,6,X) = PB(K,ipass)*ETIJ(K,ipass, 0,0,5,X)+ ETIJ(K,ipass, 1,0,5,X)&
                 & -5*ETIJ(K,ipass,0,0,4,X)/TWOB(K)
            ETIJ(K,ipass,5,0,6,X) = HPINV(K)*ETIJ(K,ipass,4,0,5,X)+ PB(K,ipass)*ETIJ(K,ipass,5,0,5,X)
            ETIJ(K,ipass,6,0,6,X) = HPINV(K)*ETIJ(K,ipass,5,0,5,X)
            ETIJ(K,ipass,1,0,6,X)= HPINV(K)*ETIJ(K,ipass,0,0,5,X)+ PB(K,ipass)*ETIJ(K,ipass,1,0,5,X)&
                 & + 2*ETIJ(K,ipass,2,0,5,X) - 5*ETIJ(K,ipass,1,0,4,X)/TWOB(K)
            ETIJ(K,ipass,2,0,6,X)= HPINV(K)*ETIJ(K,ipass,1,0,5,X)+ PB(K,ipass)*ETIJ(K,ipass,2,0,5,X)&
                 & + 3*ETIJ(K,ipass,3,0,5,X) - 5*ETIJ(K,ipass,2,0,4,X)/TWOB(K)
            ETIJ(K,ipass,3,0,6,X)= HPINV(K)*ETIJ(K,ipass,2,0,5,X)+ PB(K,ipass)*ETIJ(K,ipass,3,0,5,X)&
                 & + 4*ETIJ(K,ipass,4,0,5,X) - 5*ETIJ(K,ipass,3,0,4,X)/TWOB(K)
            ETIJ(K,ipass,4,0,6,X)= HPINV(K)*ETIJ(K,ipass,3,0,5,X)+ PB(K,ipass)*ETIJ(K,ipass,4,0,5,X)&
                 & + 5*ETIJ(K,ipass,5,0,5,X) - 5*ETIJ(K,ipass,4,0,4,X)/TWOB(K)
         ENDDO
      ENDDO
   ENDIF
ELSE
   CALL lsquit('EXTEND ICHOR_HERM_ECOEFFS (MAIN2 IN etij)',-1)
ENDIF
END SUBROUTINE ICHOR_HERM_ECOEFFS

subroutine ICHOR_Ecoeffn_maxAngP0_maxAngA0_RHS(nPrimP,Ecoeffn,ETIJ,PreExpFac)
  implicit none
  integer,intent(in) :: nPrimP
  real(realk),intent(in) :: ETIJ(nprimP,3),PreExpFac(nPrimP)
  real(realk),intent(inout) :: Ecoeffn(nPrimP)
!  
  integer :: i
  DO i = 1, nPrimP
     EcoeffN(i) = ETIJ(i,1)*ETIJ(i,2)*ETIJ(i,3)*PreExpFac(i)
  ENDDO
End subroutine ICHOR_Ecoeffn_maxAngP0_maxAngA0_RHS

subroutine ICHOR_Ecoeffn_maxAngP0_maxAngA0_LHS(nPrimP,nPass,Ecoeffn,ETIJ,PreExpFac,&
     & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
  implicit none
  integer,intent(in) :: nPrimP,nPass,MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: ETIJ(nprimP,nPass,3),PreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(inout) :: Ecoeffn(nPrimP,nPass,1,  1)
!  
  integer :: i,ii,iAtomA,iAtomB
  DO ii = 1, nPass
     iAtomA = IatomAPass(ii)
     iAtomB = IatomBPass(ii)
     DO i = 1, nPrimP
        EcoeffN(i,ii,1,1) = ETIJ(i,ii,1)*ETIJ(i,ii,2)*ETIJ(i,ii,3)*PreExpFac(i,iAtomA,iAtomB)
     ENDDO
  ENDDO
End subroutine ICHOR_Ecoeffn_maxAngP0_maxAngA0_LHS
  
subroutine ICHOR_Ecoeffn_general_LHS(nPrimP,nPass,nTUV,ijk,JMAX,l1,l2,EcoeffN,ETIJ,PreExpFac,&
     & MaxPasses,nAtomsA,nAtomsB,IatomApass,IatomBpass)
  implicit none
  integer,intent(in) :: nPrimP,nPass,nTUV,ijk,JMAX,MaxPasses,nAtomsA,nAtomsB
  integer,intent(in) :: IatomApass(MaxPasses),IatomBpass(MaxPasses)
  real(realk),intent(in) :: ETIJ(nPrimP,nPass,0:JMAX,0:l1,0:l2,3)
  real(realk),intent(in) :: PreExpFac(nPrimP,nAtomsA,nAtomsB)
  real(realk),intent(inout) :: EcoeffN(nPrimP,nPass,nTUV,ijk)
  !
  integer :: l2,l1,P2,iP2,JP2,kp2,P1,jp1,ip1,kp1,tp,up,vp,ituvp,J,ijkQ,i
  integer :: ii,iAtomA,iAtomB
  ijkQ=0
  DO P2 = 0,l2
   DO iP2=P2,0,-1
    DO jP2=P2-iP2,0,-1
     kP2=P2-iP2-jP2
     DO P1 = 0,l1
      DO iP1=P1,0,-1
       DO jP1=P1-iP1,0,-1
        kP1=P1-iP1-jP1
        IF(iP1+jP1+kP1+iP2+jP2+kP2 .EQ. l1+l2) then
         ijkQ=ijkQ+1
         DO tP=0,iP1+iP2
          DO uP=0,jP1+jP2
           DO vP=0,kP1+kP2
            ituvP = IchorTUVindexFuncFull(tp,up,vp)
            DO ii = 1, nPass
             iAtomA = IatomAPass(ii)
             iAtomB = IatomBPass(ii)
             DO i = 1, nPrimP
              EcoeffN(i,ii,ituvP,ijkQ) = ETIJ(i,ii,tp,iP1,iP2,1) &
                   & *ETIJ(i,ii,up,jP1,jP2,2)*ETIJ(i,ii,vp,kP1,kP2,3)*PreExpFac(i,iAtomA,iAtomB)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDIF
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
end subroutine ICHOR_Ecoeffn_general_LHS

subroutine ICHOR_Ecoeffn_general_RHS(nPrimP,nTUV,ijk,JMAX,l1,l2,EcoeffN,ETIJ,PreExpFac)
  implicit none
  integer,intent(in) :: nPrimP,nTUV,ijk,JMAX
  real(realk),intent(in) :: ETIJ(nPrimP,0:JMAX,0:l1,0:l2,3)
  real(realk),intent(in) :: PreExpFac(nPrimP)
  real(realk),intent(inout) :: EcoeffN(nPrimP,nTUV,ijk)
  !
  integer :: l2,l1,P2,iP2,JP2,kp2,P1,jp1,ip1,kp1,tp,up,vp,ituvp,J,ijkQ,i
  integer :: ii,iAtomA,iAtomB
  real(realk),parameter :: D1=1.0E0_realk,signQ=-1.0E0_realk
  real(realk) :: sign,signijk
  sign = signQ 
  ijkQ=0
  DO P2 = 0,l2
   DO iP2=P2,0,-1
    DO jP2=P2-iP2,0,-1
     kP2=P2-iP2-jP2
     DO P1 = 0,l1
      DO iP1=P1,0,-1
       DO jP1=P1-iP1,0,-1
        kP1=P1-iP1-jP1
        IF(iP1+jP1+kP1+iP2+jP2+kP2 .EQ. l1+l2) then
         ijkQ=ijkQ+1
         DO tP=0,iP1+iP2
          DO uP=0,jP1+jP2
           DO vP=0,kP1+kP2
            IF(MOD(tP+uP+vP,2).EQ. 0)THEN
               signijk = D1
            ELSE
               signijk = signQ
            ENDIF
            ituvP = IchorTUVindexFuncFull(tp,up,vp)
            DO i = 1, nPrimP
             EcoeffN(i,ituvP,ijkQ) = signijk*ETIJ(i,tp,iP1,iP2,1) &
                   & *ETIJ(i,up,jP1,jP2,2)*ETIJ(i,vp,kP1,kP2,3)*PreExpFac(i)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDIF
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
  ENDDO
end subroutine ICHOR_Ecoeffn_general_RHS

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
