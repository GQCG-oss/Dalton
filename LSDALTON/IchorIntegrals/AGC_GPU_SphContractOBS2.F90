MODULE AGC_GPU_OBS_Sphcontract2Mod
!Automatic Generated Code (AGC) by runSphContractOBS2.f90 in tools directory
use IchorPrecisionModule  
  
 CONTAINS
  
subroutine SphericalContractOBS2_GPU_maxAngQ2_maxAngC2(nlmP,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: nlmP,nContPasses
  real(realk),intent(in)    :: IN(nContPasses*nlmP,  6)
  real(realk),intent(inout) :: OUT(nContPasses*nlmP,  5)
  integer :: iPass,ijkP,iP
  real(realk),parameter :: SPHMAT1_3      =   -2.8867513459481292E-01_realk
  real(realk),parameter :: SPHMAT1_5      =    5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT2_1      =    1.0000000000000000E+00_realk
  real(realk),parameter :: SPHMAT4_5      =   -5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT6_3      =    5.7735026918962584E-01_realk
!$ACC PARALLEL LOOP PRIVATE(iP) PRESENT(nlmP,nContPasses,IN,OUT)
  DO iP=1,nContPasses*nlmP
    OUT(iP,1) = IN(iP,2)
    OUT(iP,2) = IN(iP,5)
    OUT(iP,3) = IN(iP,1)*SPHMAT1_3 + IN(iP,4)*SPHMAT1_3 + IN(iP,6)*SPHMAT6_3
    OUT(iP,4) = IN(iP,3)
    OUT(iP,5) = IN(iP,1)*SPHMAT1_5 + IN(iP,4)*SPHMAT4_5
  ENDDO
end subroutine SphericalContractOBS2_GPU_maxAngQ2_maxAngC2 
  
  
subroutine SphericalContractOBS2_GPU_maxAngQ2_maxAngC0(nlmP,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: nlmP,nContPasses
  real(realk),intent(in)    :: IN(nContPasses*nlmP,  6)
  real(realk),intent(inout) :: OUT(nContPasses*nlmP,  5)
  integer :: iPass,ijkP,iP
  real(realk),parameter :: SPHMAT1_3      =   -2.8867513459481292E-01_realk
  real(realk),parameter :: SPHMAT1_5      =    5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT2_1      =    1.0000000000000000E+00_realk
  real(realk),parameter :: SPHMAT4_5      =   -5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT6_3      =    5.7735026918962584E-01_realk
!$ACC PARALLEL LOOP PRIVATE(iP) PRESENT(nlmP,nContPasses,IN,OUT)
  DO iP=1,nContPasses*nlmP
    OUT(iP,1) = IN(iP,2)
    OUT(iP,2) = IN(iP,5)
    OUT(iP,3) = IN(iP,1)*SPHMAT1_3 + IN(iP,4)*SPHMAT1_3 + IN(iP,6)*SPHMAT6_3
    OUT(iP,4) = IN(iP,3)
    OUT(iP,5) = IN(iP,1)*SPHMAT1_5 + IN(iP,4)*SPHMAT4_5
  ENDDO
end subroutine SphericalContractOBS2_GPU_maxAngQ2_maxAngC0 
  
  
subroutine SphericalContractOBS2_GPU_maxAngQ3_maxAngC2(nlmP,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: nlmP,nContPasses
  real(realk),intent(in)    :: IN(nContPasses*nlmP, 18)
  real(realk),intent(inout) :: OUT(nContPasses*nlmP, 15)
  integer :: iPass,ijkP,iP
  real(realk),parameter :: SPHMAT1_3      =   -2.8867513459481292E-01_realk
  real(realk),parameter :: SPHMAT1_5      =    5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT2_1      =    1.0000000000000000E+00_realk
  real(realk),parameter :: SPHMAT4_5      =   -5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT6_3      =    5.7735026918962584E-01_realk
!$ACC PARALLEL LOOP PRIVATE(iP) PRESENT(nlmP,nContPasses,IN,OUT)
  DO iP=1,nContPasses*nlmP
    OUT(iP,1) = IN(iP,2)
    OUT(iP,2) = IN(iP,5)
    OUT(iP,3) = IN(iP,1)*SPHMAT1_3 + IN(iP,4)*SPHMAT1_3 + IN(iP,6)*SPHMAT6_3
    OUT(iP,4) = IN(iP,3)
    OUT(iP,5) = IN(iP,1)*SPHMAT1_5 + IN(iP,4)*SPHMAT4_5
    OUT(iP,6) = IN(iP,8)
    OUT(iP,7) = IN(iP,11)
    OUT(iP,8) = IN(iP,7)*SPHMAT1_3 + IN(iP,10)*SPHMAT1_3 + IN(iP,12)*SPHMAT6_3
    OUT(iP,9) = IN(iP,9)
    OUT(iP,10) = IN(iP,7)*SPHMAT1_5 + IN(iP,10)*SPHMAT4_5
    OUT(iP,11) = IN(iP,14)
    OUT(iP,12) = IN(iP,17)
    OUT(iP,13) = IN(iP,13)*SPHMAT1_3 + IN(iP,16)*SPHMAT1_3 + IN(iP,18)*SPHMAT6_3
    OUT(iP,14) = IN(iP,15)
    OUT(iP,15) = IN(iP,13)*SPHMAT1_5 + IN(iP,16)*SPHMAT4_5
  ENDDO
end subroutine SphericalContractOBS2_GPU_maxAngQ3_maxAngC2 
  
  
subroutine SphericalContractOBS2_GPU_maxAngQ3_maxAngC1(nlmP,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: nlmP,nContPasses
  real(realk),intent(in)    :: IN(nContPasses*nlmP, 18)
  real(realk),intent(inout) :: OUT(nContPasses*nlmP, 15)
  integer :: iPass,ijkP,iP
  real(realk),parameter :: SPHMAT1_7      =   -2.8867513459481292E-01_realk
  real(realk),parameter :: SPHMAT1_13     =    5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT4_1      =    1.0000000000000000E+00_realk
  real(realk),parameter :: SPHMAT10_13    =   -5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT16_7     =    5.7735026918962584E-01_realk
!$ACC PARALLEL LOOP PRIVATE(iP) PRESENT(nlmP,nContPasses,IN,OUT)
  DO iP=1,nContPasses*nlmP
    OUT(iP,1) = IN(iP,4)
    OUT(iP,2) = IN(iP,5)
    OUT(iP,3) = IN(iP,6)
    OUT(iP,4) = IN(iP,13)
    OUT(iP,5) = IN(iP,14)
    OUT(iP,6) = IN(iP,15)
    OUT(iP,7) = IN(iP,1)*SPHMAT1_7 + IN(iP,10)*SPHMAT1_7 + IN(iP,16)*SPHMAT16_7
    OUT(iP,8) = IN(iP,2)*SPHMAT1_7 + IN(iP,11)*SPHMAT1_7 + IN(iP,17)*SPHMAT16_7
    OUT(iP,9) = IN(iP,3)*SPHMAT1_7 + IN(iP,12)*SPHMAT1_7 + IN(iP,18)*SPHMAT16_7
    OUT(iP,10) = IN(iP,7)
    OUT(iP,11) = IN(iP,8)
    OUT(iP,12) = IN(iP,9)
    OUT(iP,13) = IN(iP,1)*SPHMAT1_13 + IN(iP,10)*SPHMAT10_13
    OUT(iP,14) = IN(iP,2)*SPHMAT1_13 + IN(iP,11)*SPHMAT10_13
    OUT(iP,15) = IN(iP,3)*SPHMAT1_13 + IN(iP,12)*SPHMAT10_13
  ENDDO
end subroutine SphericalContractOBS2_GPU_maxAngQ3_maxAngC1 
  
  
subroutine SphericalContractOBS2_GPU_maxAngQ4_maxAngC2(nlmP,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: nlmP,nContPasses
  real(realk),intent(in)    :: IN(nContPasses*nlmP, 36)
  real(realk),intent(inout) :: OUT(nContPasses*nlmP, 25)
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
!$ACC PARALLEL LOOP PRIVATE(iP) PRESENT(nlmP,nContPasses,IN,OUT)
  DO iP=1,nContPasses*nlmP
    OUT(iP,1) = IN(iP,8)
    OUT(iP,2) = IN(iP,11)
    OUT(iP,3) = IN(iP,7)*SPHMAT2_11 + IN(iP,10)*SPHMAT2_11 + IN(iP,12)*SPHMAT12_3
    OUT(iP,4) = IN(iP,9)
    OUT(iP,5) = IN(iP,7)*SPHMAT2_21 + IN(iP,10)*SPHMAT10_5
    OUT(iP,6) = IN(iP,26)
    OUT(iP,7) = IN(iP,29)
    OUT(iP,8) = IN(iP,25)*SPHMAT2_11 + IN(iP,28)*SPHMAT2_11 + IN(iP,30)*SPHMAT12_3
    OUT(iP,9) = IN(iP,27)
    OUT(iP,10) = IN(iP,25)*SPHMAT2_21 + IN(iP,28)*SPHMAT10_5
    OUT(iP,11) = IN(iP,2)*SPHMAT2_11 + IN(iP,20)*SPHMAT2_11 + IN(iP,32)*SPHMAT12_3
    OUT(iP,12) = IN(iP,5)*SPHMAT2_11 + IN(iP,23)*SPHMAT2_11 + IN(iP,35)*SPHMAT12_3
    OUT(iP,13) = IN(iP,1)*SPHMAT1_13 + IN(iP,4)*SPHMAT1_13 + IN(iP,6)*SPHMAT6_13 &
               & + IN(iP,19)*SPHMAT1_13 + IN(iP,22)*SPHMAT1_13 + IN(iP,24)*SPHMAT6_13 &
               & + IN(iP,31)*SPHMAT6_13 + IN(iP,34)*SPHMAT6_13 + IN(iP,36)*SPHMAT36_13
    OUT(iP,14) = IN(iP,3)*SPHMAT2_11 + IN(iP,21)*SPHMAT2_11 + IN(iP,33)*SPHMAT12_3
    OUT(iP,15) = IN(iP,1)*SPHMAT1_15 + IN(iP,4)*SPHMAT4_15 + IN(iP,19)*SPHMAT1_15 &
               & + IN(iP,22)*SPHMAT4_15 + IN(iP,31)*SPHMAT6_23 + IN(iP,34)*SPHMAT2_11
    OUT(iP,16) = IN(iP,14)
    OUT(iP,17) = IN(iP,17)
    OUT(iP,18) = IN(iP,13)*SPHMAT2_11 + IN(iP,16)*SPHMAT2_11 + IN(iP,18)*SPHMAT12_3
    OUT(iP,19) = IN(iP,15)
    OUT(iP,20) = IN(iP,13)*SPHMAT2_21 + IN(iP,16)*SPHMAT10_5
    OUT(iP,21) = IN(iP,2)*SPHMAT2_21 + IN(iP,20)*SPHMAT10_5
    OUT(iP,22) = IN(iP,5)*SPHMAT2_21 + IN(iP,23)*SPHMAT10_5
    OUT(iP,23) = IN(iP,1)*SPHMAT1_15 + IN(iP,4)*SPHMAT1_15 + IN(iP,6)*SPHMAT6_23 &
               & + IN(iP,19)*SPHMAT4_15 + IN(iP,22)*SPHMAT4_15 + IN(iP,24)*SPHMAT2_11
    OUT(iP,24) = IN(iP,3)*SPHMAT2_21 + IN(iP,21)*SPHMAT10_5
    OUT(iP,25) = IN(iP,1)*SPHMAT1_25 + IN(iP,4)*SPHMAT4_25 + IN(iP,19)*SPHMAT4_25 &
               & + IN(iP,22)*SPHMAT1_25
  ENDDO
end subroutine SphericalContractOBS2_GPU_maxAngQ4_maxAngC2 
  
  
END MODULE
