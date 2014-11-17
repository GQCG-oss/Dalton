MODULE AGC_CPU_OBS_Sphcontract2Mod
!Automatic Generated Code (AGC) by runSphContractOBS2.f90 in tools directory
use IchorPrecisionMod
  
 CONTAINS
  
 subroutine SphericalContractOBS2_CPU_maxAngQ2_maxAngC2(nlmP,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: nlmP,nContPasses
  real(realk),intent(in)    :: IN(nlmP,  6,nContPasses)
  real(realk),intent(inout) :: OUT(nlmP,  5,nContPasses)
  integer :: iPass,ijkP,iP
  real(realk),parameter :: SPHMAT1_3      =   -2.8867513459481292E-01_realk
  real(realk),parameter :: SPHMAT1_5      =    5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT2_1      =    1.0000000000000000E+00_realk
  real(realk),parameter :: SPHMAT4_5      =   -5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT6_3      =    5.7735026918962584E-01_realk
!$OMP DO PRIVATE(iPass,ijkP)
  DO iPass=1,nContPasses
   DO ijkP=1,nlmP
    OUT(ijkP,1,iPass) = IN(ijkP,2,iPass)
    OUT(ijkP,2,iPass) = IN(ijkP,5,iPass)
    OUT(ijkP,3,iPass) = IN(ijkP,1,iPass)*SPHMAT1_3 + IN(ijkP,4,iPass)*SPHMAT1_3 + IN(ijkP,6,iPass)*SPHMAT6_3
    OUT(ijkP,4,iPass) = IN(ijkP,3,iPass)
    OUT(ijkP,5,iPass) = IN(ijkP,1,iPass)*SPHMAT1_5 + IN(ijkP,4,iPass)*SPHMAT4_5
   ENDDO
  ENDDO
!$OMP END DO
end subroutine SphericalContractOBS2_CPU_maxAngQ2_maxAngC2 
  
  
 subroutine SphericalContractOBS2_CPU_maxAngQ2_maxAngC0(nlmP,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: nlmP,nContPasses
  real(realk),intent(in)    :: IN(nlmP,  6,nContPasses)
  real(realk),intent(inout) :: OUT(nlmP,  5,nContPasses)
  integer :: iPass,ijkP,iP
  real(realk),parameter :: SPHMAT1_3      =   -2.8867513459481292E-01_realk
  real(realk),parameter :: SPHMAT1_5      =    5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT2_1      =    1.0000000000000000E+00_realk
  real(realk),parameter :: SPHMAT4_5      =   -5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT6_3      =    5.7735026918962584E-01_realk
!$OMP DO PRIVATE(iPass,ijkP)
  DO iPass=1,nContPasses
   DO ijkP=1,nlmP
    OUT(ijkP,1,iPass) = IN(ijkP,2,iPass)
    OUT(ijkP,2,iPass) = IN(ijkP,5,iPass)
    OUT(ijkP,3,iPass) = IN(ijkP,1,iPass)*SPHMAT1_3 + IN(ijkP,4,iPass)*SPHMAT1_3 + IN(ijkP,6,iPass)*SPHMAT6_3
    OUT(ijkP,4,iPass) = IN(ijkP,3,iPass)
    OUT(ijkP,5,iPass) = IN(ijkP,1,iPass)*SPHMAT1_5 + IN(ijkP,4,iPass)*SPHMAT4_5
   ENDDO
  ENDDO
!$OMP END DO
end subroutine SphericalContractOBS2_CPU_maxAngQ2_maxAngC0 
  
  
 subroutine SphericalContractOBS2_CPU_maxAngQ3_maxAngC2(nlmP,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: nlmP,nContPasses
  real(realk),intent(in)    :: IN(nlmP, 18,nContPasses)
  real(realk),intent(inout) :: OUT(nlmP, 15,nContPasses)
  integer :: iPass,ijkP,iP
  real(realk),parameter :: SPHMAT1_3      =   -2.8867513459481292E-01_realk
  real(realk),parameter :: SPHMAT1_5      =    5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT2_1      =    1.0000000000000000E+00_realk
  real(realk),parameter :: SPHMAT4_5      =   -5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT6_3      =    5.7735026918962584E-01_realk
!$OMP DO PRIVATE(iPass,ijkP)
  DO iPass=1,nContPasses
   DO ijkP=1,nlmP
    OUT(ijkP,1,iPass) = IN(ijkP,2,iPass)
    OUT(ijkP,2,iPass) = IN(ijkP,5,iPass)
    OUT(ijkP,3,iPass) = IN(ijkP,1,iPass)*SPHMAT1_3 + IN(ijkP,4,iPass)*SPHMAT1_3 + IN(ijkP,6,iPass)*SPHMAT6_3
    OUT(ijkP,4,iPass) = IN(ijkP,3,iPass)
    OUT(ijkP,5,iPass) = IN(ijkP,1,iPass)*SPHMAT1_5 + IN(ijkP,4,iPass)*SPHMAT4_5
    OUT(ijkP,6,iPass) = IN(ijkP,8,iPass)
    OUT(ijkP,7,iPass) = IN(ijkP,11,iPass)
    OUT(ijkP,8,iPass) = IN(ijkP,7,iPass)*SPHMAT1_3 + IN(ijkP,10,iPass)*SPHMAT1_3 + IN(ijkP,12,iPass)*SPHMAT6_3
    OUT(ijkP,9,iPass) = IN(ijkP,9,iPass)
    OUT(ijkP,10,iPass) = IN(ijkP,7,iPass)*SPHMAT1_5 + IN(ijkP,10,iPass)*SPHMAT4_5
    OUT(ijkP,11,iPass) = IN(ijkP,14,iPass)
    OUT(ijkP,12,iPass) = IN(ijkP,17,iPass)
    OUT(ijkP,13,iPass) = IN(ijkP,13,iPass)*SPHMAT1_3 + IN(ijkP,16,iPass)*SPHMAT1_3 + IN(ijkP,18,iPass)*SPHMAT6_3
    OUT(ijkP,14,iPass) = IN(ijkP,15,iPass)
    OUT(ijkP,15,iPass) = IN(ijkP,13,iPass)*SPHMAT1_5 + IN(ijkP,16,iPass)*SPHMAT4_5
   ENDDO
  ENDDO
!$OMP END DO
end subroutine SphericalContractOBS2_CPU_maxAngQ3_maxAngC2 
  
  
 subroutine SphericalContractOBS2_CPU_maxAngQ3_maxAngC1(nlmP,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: nlmP,nContPasses
  real(realk),intent(in)    :: IN(nlmP, 18,nContPasses)
  real(realk),intent(inout) :: OUT(nlmP, 15,nContPasses)
  integer :: iPass,ijkP,iP
  real(realk),parameter :: SPHMAT1_7      =   -2.8867513459481292E-01_realk
  real(realk),parameter :: SPHMAT1_13     =    5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT4_1      =    1.0000000000000000E+00_realk
  real(realk),parameter :: SPHMAT10_13    =   -5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT16_7     =    5.7735026918962584E-01_realk
!$OMP DO PRIVATE(iPass,ijkP)
  DO iPass=1,nContPasses
   DO ijkP=1,nlmP
    OUT(ijkP,1,iPass) = IN(ijkP,4,iPass)
    OUT(ijkP,2,iPass) = IN(ijkP,5,iPass)
    OUT(ijkP,3,iPass) = IN(ijkP,6,iPass)
    OUT(ijkP,4,iPass) = IN(ijkP,13,iPass)
    OUT(ijkP,5,iPass) = IN(ijkP,14,iPass)
    OUT(ijkP,6,iPass) = IN(ijkP,15,iPass)
    OUT(ijkP,7,iPass) = IN(ijkP,1,iPass)*SPHMAT1_7 + IN(ijkP,10,iPass)*SPHMAT1_7 + IN(ijkP,16,iPass)*SPHMAT16_7
    OUT(ijkP,8,iPass) = IN(ijkP,2,iPass)*SPHMAT1_7 + IN(ijkP,11,iPass)*SPHMAT1_7 + IN(ijkP,17,iPass)*SPHMAT16_7
    OUT(ijkP,9,iPass) = IN(ijkP,3,iPass)*SPHMAT1_7 + IN(ijkP,12,iPass)*SPHMAT1_7 + IN(ijkP,18,iPass)*SPHMAT16_7
    OUT(ijkP,10,iPass) = IN(ijkP,7,iPass)
    OUT(ijkP,11,iPass) = IN(ijkP,8,iPass)
    OUT(ijkP,12,iPass) = IN(ijkP,9,iPass)
    OUT(ijkP,13,iPass) = IN(ijkP,1,iPass)*SPHMAT1_13 + IN(ijkP,10,iPass)*SPHMAT10_13
    OUT(ijkP,14,iPass) = IN(ijkP,2,iPass)*SPHMAT1_13 + IN(ijkP,11,iPass)*SPHMAT10_13
    OUT(ijkP,15,iPass) = IN(ijkP,3,iPass)*SPHMAT1_13 + IN(ijkP,12,iPass)*SPHMAT10_13
   ENDDO
  ENDDO
!$OMP END DO
end subroutine SphericalContractOBS2_CPU_maxAngQ3_maxAngC1 
  
  
 subroutine SphericalContractOBS2_CPU_maxAngQ4_maxAngC2(nlmP,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: nlmP,nContPasses
  real(realk),intent(in)    :: IN(nlmP, 36,nContPasses)
  real(realk),intent(inout) :: OUT(nlmP, 25,nContPasses)
  integer :: iPass,ijkP,iP
  real(realk),parameter :: SPHMAT1_13     =    8.3333333333333356E-02_realk
  real(realk),parameter :: SPHMAT1_15     =   -1.4433756729740646E-01_realk
  real(realk),parameter :: SPHMAT1_25     =    2.5000000000000000E-01_realk
  real(realk),parameter :: SPHMAT2_11     =   -2.8867513459481292E-01_realk
  real(realk),parameter :: SPHMAT2_21     =    5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT4_15     =    1.4433756729740646E-01_realk
  real(realk),parameter :: SPHMAT4_25     =   -2.5000000000000000E-01_realk
  real(realk),parameter :: SPHMAT6_13     =   -1.6666666666666671E-01_realk
  real(realk),parameter :: SPHMAT6_23     =    2.8867513459481292E-01_realk
  real(realk),parameter :: SPHMAT8_1      =    1.0000000000000000E+00_realk
  real(realk),parameter :: SPHMAT10_5     =   -5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT12_3     =    5.7735026918962584E-01_realk
  real(realk),parameter :: SPHMAT36_13    =    3.3333333333333343E-01_realk
!$OMP DO PRIVATE(iPass,ijkP)
  DO iPass=1,nContPasses
   DO ijkP=1,nlmP
    OUT(ijkP,1,iPass) = IN(ijkP,8,iPass)
    OUT(ijkP,2,iPass) = IN(ijkP,11,iPass)
    OUT(ijkP,3,iPass) = IN(ijkP,7,iPass)*SPHMAT2_11 + IN(ijkP,10,iPass)*SPHMAT2_11 + IN(ijkP,12,iPass)*SPHMAT12_3
    OUT(ijkP,4,iPass) = IN(ijkP,9,iPass)
    OUT(ijkP,5,iPass) = IN(ijkP,7,iPass)*SPHMAT2_21 + IN(ijkP,10,iPass)*SPHMAT10_5
    OUT(ijkP,6,iPass) = IN(ijkP,26,iPass)
    OUT(ijkP,7,iPass) = IN(ijkP,29,iPass)
    OUT(ijkP,8,iPass) = IN(ijkP,25,iPass)*SPHMAT2_11 + IN(ijkP,28,iPass)*SPHMAT2_11 + IN(ijkP,30,iPass)*SPHMAT12_3
    OUT(ijkP,9,iPass) = IN(ijkP,27,iPass)
    OUT(ijkP,10,iPass) = IN(ijkP,25,iPass)*SPHMAT2_21 + IN(ijkP,28,iPass)*SPHMAT10_5
    OUT(ijkP,11,iPass) = IN(ijkP,2,iPass)*SPHMAT2_11 + IN(ijkP,20,iPass)*SPHMAT2_11 + IN(ijkP,32,iPass)*SPHMAT12_3
    OUT(ijkP,12,iPass) = IN(ijkP,5,iPass)*SPHMAT2_11 + IN(ijkP,23,iPass)*SPHMAT2_11 + IN(ijkP,35,iPass)*SPHMAT12_3
    OUT(ijkP,13,iPass) = IN(ijkP,1,iPass)*SPHMAT1_13 + IN(ijkP,4,iPass)*SPHMAT1_13 + IN(ijkP,6,iPass)*SPHMAT6_13 &
               & + IN(ijkP,19,iPass)*SPHMAT1_13 + IN(ijkP,22,iPass)*SPHMAT1_13 + IN(ijkP,24,iPass)*SPHMAT6_13 &
               & + IN(ijkP,31,iPass)*SPHMAT6_13 + IN(ijkP,34,iPass)*SPHMAT6_13 + IN(ijkP,36,iPass)*SPHMAT36_13
    OUT(ijkP,14,iPass) = IN(ijkP,3,iPass)*SPHMAT2_11 + IN(ijkP,21,iPass)*SPHMAT2_11 + IN(ijkP,33,iPass)*SPHMAT12_3
    OUT(ijkP,15,iPass) = IN(ijkP,1,iPass)*SPHMAT1_15 + IN(ijkP,4,iPass)*SPHMAT4_15 + IN(ijkP,19,iPass)*SPHMAT1_15 &
               & + IN(ijkP,22,iPass)*SPHMAT4_15 + IN(ijkP,31,iPass)*SPHMAT6_23 + IN(ijkP,34,iPass)*SPHMAT2_11
    OUT(ijkP,16,iPass) = IN(ijkP,14,iPass)
    OUT(ijkP,17,iPass) = IN(ijkP,17,iPass)
    OUT(ijkP,18,iPass) = IN(ijkP,13,iPass)*SPHMAT2_11 + IN(ijkP,16,iPass)*SPHMAT2_11 + IN(ijkP,18,iPass)*SPHMAT12_3
    OUT(ijkP,19,iPass) = IN(ijkP,15,iPass)
    OUT(ijkP,20,iPass) = IN(ijkP,13,iPass)*SPHMAT2_21 + IN(ijkP,16,iPass)*SPHMAT10_5
    OUT(ijkP,21,iPass) = IN(ijkP,2,iPass)*SPHMAT2_21 + IN(ijkP,20,iPass)*SPHMAT10_5
    OUT(ijkP,22,iPass) = IN(ijkP,5,iPass)*SPHMAT2_21 + IN(ijkP,23,iPass)*SPHMAT10_5
    OUT(ijkP,23,iPass) = IN(ijkP,1,iPass)*SPHMAT1_15 + IN(ijkP,4,iPass)*SPHMAT1_15 + IN(ijkP,6,iPass)*SPHMAT6_23 &
               & + IN(ijkP,19,iPass)*SPHMAT4_15 + IN(ijkP,22,iPass)*SPHMAT4_15 + IN(ijkP,24,iPass)*SPHMAT2_11
    OUT(ijkP,24,iPass) = IN(ijkP,3,iPass)*SPHMAT2_21 + IN(ijkP,21,iPass)*SPHMAT10_5
    OUT(ijkP,25,iPass) = IN(ijkP,1,iPass)*SPHMAT1_25 + IN(ijkP,4,iPass)*SPHMAT4_25 + IN(ijkP,19,iPass)*SPHMAT4_25 &
               & + IN(ijkP,22,iPass)*SPHMAT1_25
   ENDDO
  ENDDO
!$OMP END DO
end subroutine SphericalContractOBS2_CPU_maxAngQ4_maxAngC2 
  
  
END MODULE AGC_CPU_OBS_Sphcontract2Mod
