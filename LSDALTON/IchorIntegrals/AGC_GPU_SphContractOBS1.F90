MODULE AGC_GPU_OBS_Sphcontract1Mod
!Automatic Generated Code (AGC) by runSphContractOBS1.f90 in tools directory
use IchorPrecisionModule  
  
 CONTAINS
  
subroutine SphericalContractOBS1_GPU_maxAngP2_maxAngA2(ijkQcart,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: ijkQcart,nContPasses
  real(realk),intent(in)    :: IN(nContPasses,  6,ijkQcart)
  real(realk),intent(inout) :: OUT(nContPasses,  5,ijkQcart)
  integer :: iP,ijkQ
  real(realk),parameter :: SPHMAT1_3      =   -2.8867513459481292E-01_realk
  real(realk),parameter :: SPHMAT1_5      =    5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT2_1      =    1.0000000000000000E+00_realk
  real(realk),parameter :: SPHMAT4_5      =   -5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT6_3      =    5.7735026918962584E-01_realk
!$ACC PARALLEL LOOP PRIVATE(iP,ijkQ) PRESENT(nContPasses,ijkQcart,IN,OUT)
  DO iP=1,nContPasses
   DO ijkQ=1,ijkQcart
    OUT(iP,1,ijkQ) = IN(iP,2,ijkQ)
    OUT(iP,2,ijkQ) = IN(iP,5,ijkQ)
    OUT(iP,3,ijkQ) = IN(iP,1,ijkQ)*SPHMAT1_3 + IN(iP,4,ijkQ)*SPHMAT1_3 + IN(iP,6,ijkQ)*SPHMAT6_3
    OUT(iP,4,ijkQ) = IN(iP,3,ijkQ)
    OUT(iP,5,ijkQ) = IN(iP,1,ijkQ)*SPHMAT1_5 + IN(iP,4,ijkQ)*SPHMAT4_5
   ENDDO
  ENDDO
end subroutine SphericalContractOBS1_GPU_maxAngP2_maxAngA2 
  
  
subroutine SphericalContractOBS1_GPU_maxAngP2_maxAngA0(ijkQcart,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: ijkQcart,nContPasses
  real(realk),intent(in)    :: IN(nContPasses,  6,ijkQcart)
  real(realk),intent(inout) :: OUT(nContPasses,  5,ijkQcart)
  integer :: iP,ijkQ
  real(realk),parameter :: SPHMAT1_3      =   -2.8867513459481292E-01_realk
  real(realk),parameter :: SPHMAT1_5      =    5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT2_1      =    1.0000000000000000E+00_realk
  real(realk),parameter :: SPHMAT4_5      =   -5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT6_3      =    5.7735026918962584E-01_realk
!$ACC PARALLEL LOOP PRIVATE(iP,ijkQ) PRESENT(nContPasses,ijkQcart,IN,OUT)
  DO iP=1,nContPasses
   DO ijkQ=1,ijkQcart
    OUT(iP,1,ijkQ) = IN(iP,2,ijkQ)
    OUT(iP,2,ijkQ) = IN(iP,5,ijkQ)
    OUT(iP,3,ijkQ) = IN(iP,1,ijkQ)*SPHMAT1_3 + IN(iP,4,ijkQ)*SPHMAT1_3 + IN(iP,6,ijkQ)*SPHMAT6_3
    OUT(iP,4,ijkQ) = IN(iP,3,ijkQ)
    OUT(iP,5,ijkQ) = IN(iP,1,ijkQ)*SPHMAT1_5 + IN(iP,4,ijkQ)*SPHMAT4_5
   ENDDO
  ENDDO
end subroutine SphericalContractOBS1_GPU_maxAngP2_maxAngA0 
  
  
subroutine SphericalContractOBS1_GPU_maxAngP3_maxAngA2(ijkQcart,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: ijkQcart,nContPasses
  real(realk),intent(in)    :: IN(nContPasses, 18,ijkQcart)
  real(realk),intent(inout) :: OUT(nContPasses, 15,ijkQcart)
  integer :: iP,ijkQ
  real(realk),parameter :: SPHMAT1_3      =   -2.8867513459481292E-01_realk
  real(realk),parameter :: SPHMAT1_5      =    5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT2_1      =    1.0000000000000000E+00_realk
  real(realk),parameter :: SPHMAT4_5      =   -5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT6_3      =    5.7735026918962584E-01_realk
!$ACC PARALLEL LOOP PRIVATE(iP,ijkQ) PRESENT(nContPasses,ijkQcart,IN,OUT)
  DO iP=1,nContPasses
   DO ijkQ=1,ijkQcart
    OUT(iP,1,ijkQ) = IN(iP,2,ijkQ)
    OUT(iP,2,ijkQ) = IN(iP,5,ijkQ)
    OUT(iP,3,ijkQ) = IN(iP,1,ijkQ)*SPHMAT1_3 + IN(iP,4,ijkQ)*SPHMAT1_3 + IN(iP,6,ijkQ)*SPHMAT6_3
    OUT(iP,4,ijkQ) = IN(iP,3,ijkQ)
    OUT(iP,5,ijkQ) = IN(iP,1,ijkQ)*SPHMAT1_5 + IN(iP,4,ijkQ)*SPHMAT4_5
    OUT(iP,6,ijkQ) = IN(iP,8,ijkQ)
    OUT(iP,7,ijkQ) = IN(iP,11,ijkQ)
    OUT(iP,8,ijkQ) = IN(iP,7,ijkQ)*SPHMAT1_3 + IN(iP,10,ijkQ)*SPHMAT1_3 + IN(iP,12,ijkQ)*SPHMAT6_3
    OUT(iP,9,ijkQ) = IN(iP,9,ijkQ)
    OUT(iP,10,ijkQ) = IN(iP,7,ijkQ)*SPHMAT1_5 + IN(iP,10,ijkQ)*SPHMAT4_5
    OUT(iP,11,ijkQ) = IN(iP,14,ijkQ)
    OUT(iP,12,ijkQ) = IN(iP,17,ijkQ)
    OUT(iP,13,ijkQ) = IN(iP,13,ijkQ)*SPHMAT1_3 + IN(iP,16,ijkQ)*SPHMAT1_3 + IN(iP,18,ijkQ)*SPHMAT6_3
    OUT(iP,14,ijkQ) = IN(iP,15,ijkQ)
    OUT(iP,15,ijkQ) = IN(iP,13,ijkQ)*SPHMAT1_5 + IN(iP,16,ijkQ)*SPHMAT4_5
   ENDDO
  ENDDO
end subroutine SphericalContractOBS1_GPU_maxAngP3_maxAngA2 
  
  
subroutine SphericalContractOBS1_GPU_maxAngP3_maxAngA1(ijkQcart,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: ijkQcart,nContPasses
  real(realk),intent(in)    :: IN(nContPasses, 18,ijkQcart)
  real(realk),intent(inout) :: OUT(nContPasses, 15,ijkQcart)
  integer :: iP,ijkQ
  real(realk),parameter :: SPHMAT1_7      =   -2.8867513459481292E-01_realk
  real(realk),parameter :: SPHMAT1_13     =    5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT4_1      =    1.0000000000000000E+00_realk
  real(realk),parameter :: SPHMAT10_13    =   -5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT16_7     =    5.7735026918962584E-01_realk
!$ACC PARALLEL LOOP PRIVATE(iP,ijkQ) PRESENT(nContPasses,ijkQcart,IN,OUT)
  DO iP=1,nContPasses
   DO ijkQ=1,ijkQcart
    OUT(iP,1,ijkQ) = IN(iP,4,ijkQ)
    OUT(iP,2,ijkQ) = IN(iP,5,ijkQ)
    OUT(iP,3,ijkQ) = IN(iP,6,ijkQ)
    OUT(iP,4,ijkQ) = IN(iP,13,ijkQ)
    OUT(iP,5,ijkQ) = IN(iP,14,ijkQ)
    OUT(iP,6,ijkQ) = IN(iP,15,ijkQ)
    OUT(iP,7,ijkQ) = IN(iP,1,ijkQ)*SPHMAT1_7 + IN(iP,10,ijkQ)*SPHMAT1_7 + IN(iP,16,ijkQ)*SPHMAT16_7
    OUT(iP,8,ijkQ) = IN(iP,2,ijkQ)*SPHMAT1_7 + IN(iP,11,ijkQ)*SPHMAT1_7 + IN(iP,17,ijkQ)*SPHMAT16_7
    OUT(iP,9,ijkQ) = IN(iP,3,ijkQ)*SPHMAT1_7 + IN(iP,12,ijkQ)*SPHMAT1_7 + IN(iP,18,ijkQ)*SPHMAT16_7
    OUT(iP,10,ijkQ) = IN(iP,7,ijkQ)
    OUT(iP,11,ijkQ) = IN(iP,8,ijkQ)
    OUT(iP,12,ijkQ) = IN(iP,9,ijkQ)
    OUT(iP,13,ijkQ) = IN(iP,1,ijkQ)*SPHMAT1_13 + IN(iP,10,ijkQ)*SPHMAT10_13
    OUT(iP,14,ijkQ) = IN(iP,2,ijkQ)*SPHMAT1_13 + IN(iP,11,ijkQ)*SPHMAT10_13
    OUT(iP,15,ijkQ) = IN(iP,3,ijkQ)*SPHMAT1_13 + IN(iP,12,ijkQ)*SPHMAT10_13
   ENDDO
  ENDDO
end subroutine SphericalContractOBS1_GPU_maxAngP3_maxAngA1 
  
  
subroutine SphericalContractOBS1_GPU_maxAngP4_maxAngA2(ijkQcart,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: ijkQcart,nContPasses
  real(realk),intent(in)    :: IN(nContPasses, 36,ijkQcart)
  real(realk),intent(inout) :: OUT(nContPasses, 25,ijkQcart)
  integer :: iP,ijkQ
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
!$ACC PARALLEL LOOP PRIVATE(iP,ijkQ) PRESENT(nContPasses,ijkQcart,IN,OUT)
  DO iP=1,nContPasses
   DO ijkQ=1,ijkQcart
    OUT(iP,1,ijkQ) = IN(iP,8,ijkQ)
    OUT(iP,2,ijkQ) = IN(iP,11,ijkQ)
    OUT(iP,3,ijkQ) = IN(iP,7,ijkQ)*SPHMAT2_11 + IN(iP,10,ijkQ)*SPHMAT2_11 + IN(iP,12,ijkQ)*SPHMAT12_3
    OUT(iP,4,ijkQ) = IN(iP,9,ijkQ)
    OUT(iP,5,ijkQ) = IN(iP,7,ijkQ)*SPHMAT2_21 + IN(iP,10,ijkQ)*SPHMAT10_5
    OUT(iP,6,ijkQ) = IN(iP,26,ijkQ)
    OUT(iP,7,ijkQ) = IN(iP,29,ijkQ)
    OUT(iP,8,ijkQ) = IN(iP,25,ijkQ)*SPHMAT2_11 + IN(iP,28,ijkQ)*SPHMAT2_11 + IN(iP,30,ijkQ)*SPHMAT12_3
    OUT(iP,9,ijkQ) = IN(iP,27,ijkQ)
    OUT(iP,10,ijkQ) = IN(iP,25,ijkQ)*SPHMAT2_21 + IN(iP,28,ijkQ)*SPHMAT10_5
    OUT(iP,11,ijkQ) = IN(iP,2,ijkQ)*SPHMAT2_11 + IN(iP,20,ijkQ)*SPHMAT2_11 + IN(iP,32,ijkQ)*SPHMAT12_3
    OUT(iP,12,ijkQ) = IN(iP,5,ijkQ)*SPHMAT2_11 + IN(iP,23,ijkQ)*SPHMAT2_11 + IN(iP,35,ijkQ)*SPHMAT12_3
    OUT(iP,13,ijkQ) = IN(iP,1,ijkQ)*SPHMAT1_13 + IN(iP,4,ijkQ)*SPHMAT1_13 + IN(iP,6,ijkQ)*SPHMAT6_13 + IN(iP,19,ijkQ)*SPHMAT1_13 &
               & + IN(iP,22,ijkQ)*SPHMAT1_13 + IN(iP,24,ijkQ)*SPHMAT6_13 + IN(iP,31,ijkQ)*SPHMAT6_13 + IN(iP,34,ijkQ)*SPHMAT6_13 &
               & + IN(iP,36,ijkQ)*SPHMAT36_13
    OUT(iP,14,ijkQ) = IN(iP,3,ijkQ)*SPHMAT2_11 + IN(iP,21,ijkQ)*SPHMAT2_11 + IN(iP,33,ijkQ)*SPHMAT12_3
    OUT(iP,15,ijkQ) = IN(iP,1,ijkQ)*SPHMAT1_15 + IN(iP,4,ijkQ)*SPHMAT4_15 + IN(iP,19,ijkQ)*SPHMAT1_15 + IN(iP,22,ijkQ)*SPHMAT4_15 &
               & + IN(iP,31,ijkQ)*SPHMAT6_23 + IN(iP,34,ijkQ)*SPHMAT2_11
    OUT(iP,16,ijkQ) = IN(iP,14,ijkQ)
    OUT(iP,17,ijkQ) = IN(iP,17,ijkQ)
    OUT(iP,18,ijkQ) = IN(iP,13,ijkQ)*SPHMAT2_11 + IN(iP,16,ijkQ)*SPHMAT2_11 + IN(iP,18,ijkQ)*SPHMAT12_3
    OUT(iP,19,ijkQ) = IN(iP,15,ijkQ)
    OUT(iP,20,ijkQ) = IN(iP,13,ijkQ)*SPHMAT2_21 + IN(iP,16,ijkQ)*SPHMAT10_5
    OUT(iP,21,ijkQ) = IN(iP,2,ijkQ)*SPHMAT2_21 + IN(iP,20,ijkQ)*SPHMAT10_5
    OUT(iP,22,ijkQ) = IN(iP,5,ijkQ)*SPHMAT2_21 + IN(iP,23,ijkQ)*SPHMAT10_5
    OUT(iP,23,ijkQ) = IN(iP,1,ijkQ)*SPHMAT1_15 + IN(iP,4,ijkQ)*SPHMAT1_15 + IN(iP,6,ijkQ)*SPHMAT6_23 + IN(iP,19,ijkQ)*SPHMAT4_15 &
               & + IN(iP,22,ijkQ)*SPHMAT4_15 + IN(iP,24,ijkQ)*SPHMAT2_11
    OUT(iP,24,ijkQ) = IN(iP,3,ijkQ)*SPHMAT2_21 + IN(iP,21,ijkQ)*SPHMAT10_5
    OUT(iP,25,ijkQ) = IN(iP,1,ijkQ)*SPHMAT1_25 + IN(iP,4,ijkQ)*SPHMAT4_25 + IN(iP,19,ijkQ)*SPHMAT4_25 + IN(iP,22,ijkQ)*SPHMAT1_25
   ENDDO
  ENDDO
end subroutine SphericalContractOBS1_GPU_maxAngP4_maxAngA2 
  
  
END MODULE 
