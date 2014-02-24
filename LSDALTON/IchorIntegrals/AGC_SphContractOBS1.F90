MODULE AGC_OBS_Sphcontract1Mod
!Automatic Generated Code (AGC) by runSphContractOBS1.f90 in tools directory
use IchorPrecisionModule  
  
 CONTAINS
  
subroutine SphericalContractOBS1_maxAngP2_maxAngA2(ijkQcart,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: ijkQcart,nContPasses
  real(realk),intent(in)    :: IN(  6,ijkQcart,nContPasses)
  real(realk),intent(inout) :: OUT(  5,ijkQcart,nContPasses)
  integer :: iPass,ijkQ
  real(realk),parameter :: SPHMAT1_3      =   -2.8867513459481292E-01_realk
  real(realk),parameter :: SPHMAT1_5      =    5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT2_1      =    1.0000000000000000E+00_realk
  real(realk),parameter :: SPHMAT4_5      =   -5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT6_3      =    5.7735026918962584E-01_realk
  DO iPass=1,nContPasses
   DO ijkQ=1,ijkQcart
    OUT(  1,ijkQ,iPass) = IN(  2,ijkQ,iPass)
    OUT(  2,ijkQ,iPass) = IN(  5,ijkQ,iPass)
    OUT(  3,ijkQ,iPass) = IN(  1,ijkQ,iPass)*SPHMAT1_3      
    OUT(  3,ijkQ,iPass) = OUT(  3,ijkQ,iPass) + IN(  4,ijkQ,iPass)*SPHMAT1_3      
    OUT(  3,ijkQ,iPass) = OUT(  3,ijkQ,iPass) + IN(  6,ijkQ,iPass)*SPHMAT6_3      
    OUT(  4,ijkQ,iPass) = IN(  3,ijkQ,iPass)
    OUT(  5,ijkQ,iPass) = IN(  1,ijkQ,iPass)*SPHMAT1_5      
    OUT(  5,ijkQ,iPass) = OUT(  5,ijkQ,iPass) + IN(  4,ijkQ,iPass)*SPHMAT4_5      
   ENDDO
  ENDDO
end subroutine SphericalContractOBS1_maxAngP2_maxAngA2 
  
  
subroutine SphericalContractOBS1_maxAngP2_maxAngA0(ijkQcart,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: ijkQcart,nContPasses
  real(realk),intent(in)    :: IN(  6,ijkQcart,nContPasses)
  real(realk),intent(inout) :: OUT(  5,ijkQcart,nContPasses)
  integer :: iPass,ijkQ
  real(realk),parameter :: SPHMAT1_3      =   -2.8867513459481292E-01_realk
  real(realk),parameter :: SPHMAT1_5      =    5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT2_1      =    1.0000000000000000E+00_realk
  real(realk),parameter :: SPHMAT4_5      =   -5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT6_3      =    5.7735026918962584E-01_realk
  DO iPass=1,nContPasses
   DO ijkQ=1,ijkQcart
    OUT(  1,ijkQ,iPass) = IN(  2,ijkQ,iPass)
    OUT(  2,ijkQ,iPass) = IN(  5,ijkQ,iPass)
    OUT(  3,ijkQ,iPass) = IN(  1,ijkQ,iPass)*SPHMAT1_3      
    OUT(  3,ijkQ,iPass) = OUT(  3,ijkQ,iPass) + IN(  4,ijkQ,iPass)*SPHMAT1_3      
    OUT(  3,ijkQ,iPass) = OUT(  3,ijkQ,iPass) + IN(  6,ijkQ,iPass)*SPHMAT6_3      
    OUT(  4,ijkQ,iPass) = IN(  3,ijkQ,iPass)
    OUT(  5,ijkQ,iPass) = IN(  1,ijkQ,iPass)*SPHMAT1_5      
    OUT(  5,ijkQ,iPass) = OUT(  5,ijkQ,iPass) + IN(  4,ijkQ,iPass)*SPHMAT4_5      
   ENDDO
  ENDDO
end subroutine SphericalContractOBS1_maxAngP2_maxAngA0 
  
  
subroutine SphericalContractOBS1_maxAngP3_maxAngA2(ijkQcart,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: ijkQcart,nContPasses
  real(realk),intent(in)    :: IN( 18,ijkQcart,nContPasses)
  real(realk),intent(inout) :: OUT( 15,ijkQcart,nContPasses)
  integer :: iPass,ijkQ
  real(realk),parameter :: SPHMAT1_3      =   -2.8867513459481292E-01_realk
  real(realk),parameter :: SPHMAT1_5      =    5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT2_1      =    1.0000000000000000E+00_realk
  real(realk),parameter :: SPHMAT4_5      =   -5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT6_3      =    5.7735026918962584E-01_realk
  DO iPass=1,nContPasses
   DO ijkQ=1,ijkQcart
    OUT(  1,ijkQ,iPass) = IN(  2,ijkQ,iPass)
    OUT(  2,ijkQ,iPass) = IN(  5,ijkQ,iPass)
    OUT(  3,ijkQ,iPass) = IN(  1,ijkQ,iPass)*SPHMAT1_3      
    OUT(  3,ijkQ,iPass) = OUT(  3,ijkQ,iPass) + IN(  4,ijkQ,iPass)*SPHMAT1_3      
    OUT(  3,ijkQ,iPass) = OUT(  3,ijkQ,iPass) + IN(  6,ijkQ,iPass)*SPHMAT6_3      
    OUT(  4,ijkQ,iPass) = IN(  3,ijkQ,iPass)
    OUT(  5,ijkQ,iPass) = IN(  1,ijkQ,iPass)*SPHMAT1_5      
    OUT(  5,ijkQ,iPass) = OUT(  5,ijkQ,iPass) + IN(  4,ijkQ,iPass)*SPHMAT4_5      
    OUT(  6,ijkQ,iPass) = IN(  8,ijkQ,iPass)
    OUT(  7,ijkQ,iPass) = IN( 11,ijkQ,iPass)
    OUT(  8,ijkQ,iPass) = IN(  7,ijkQ,iPass)*SPHMAT1_3      
    OUT(  8,ijkQ,iPass) = OUT(  8,ijkQ,iPass) + IN( 10,ijkQ,iPass)*SPHMAT1_3      
    OUT(  8,ijkQ,iPass) = OUT(  8,ijkQ,iPass) + IN( 12,ijkQ,iPass)*SPHMAT6_3      
    OUT(  9,ijkQ,iPass) = IN(  9,ijkQ,iPass)
    OUT( 10,ijkQ,iPass) = IN(  7,ijkQ,iPass)*SPHMAT1_5      
    OUT( 10,ijkQ,iPass) = OUT( 10,ijkQ,iPass) + IN( 10,ijkQ,iPass)*SPHMAT4_5      
    OUT( 11,ijkQ,iPass) = IN( 14,ijkQ,iPass)
    OUT( 12,ijkQ,iPass) = IN( 17,ijkQ,iPass)
    OUT( 13,ijkQ,iPass) = IN( 13,ijkQ,iPass)*SPHMAT1_3      
    OUT( 13,ijkQ,iPass) = OUT( 13,ijkQ,iPass) + IN( 16,ijkQ,iPass)*SPHMAT1_3      
    OUT( 13,ijkQ,iPass) = OUT( 13,ijkQ,iPass) + IN( 18,ijkQ,iPass)*SPHMAT6_3      
    OUT( 14,ijkQ,iPass) = IN( 15,ijkQ,iPass)
    OUT( 15,ijkQ,iPass) = IN( 13,ijkQ,iPass)*SPHMAT1_5      
    OUT( 15,ijkQ,iPass) = OUT( 15,ijkQ,iPass) + IN( 16,ijkQ,iPass)*SPHMAT4_5      
   ENDDO
  ENDDO
end subroutine SphericalContractOBS1_maxAngP3_maxAngA2 
  
  
subroutine SphericalContractOBS1_maxAngP3_maxAngA1(ijkQcart,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: ijkQcart,nContPasses
  real(realk),intent(in)    :: IN( 18,ijkQcart,nContPasses)
  real(realk),intent(inout) :: OUT( 15,ijkQcart,nContPasses)
  integer :: iPass,ijkQ
  real(realk),parameter :: SPHMAT1_7      =   -2.8867513459481292E-01_realk
  real(realk),parameter :: SPHMAT1_13     =    5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT4_1      =    1.0000000000000000E+00_realk
  real(realk),parameter :: SPHMAT10_13    =   -5.0000000000000000E-01_realk
  real(realk),parameter :: SPHMAT16_7     =    5.7735026918962584E-01_realk
  DO iPass=1,nContPasses
   DO ijkQ=1,ijkQcart
    OUT(  1,ijkQ,iPass) = IN(  4,ijkQ,iPass)
    OUT(  2,ijkQ,iPass) = IN(  5,ijkQ,iPass)
    OUT(  3,ijkQ,iPass) = IN(  6,ijkQ,iPass)
    OUT(  4,ijkQ,iPass) = IN( 13,ijkQ,iPass)
    OUT(  5,ijkQ,iPass) = IN( 14,ijkQ,iPass)
    OUT(  6,ijkQ,iPass) = IN( 15,ijkQ,iPass)
    OUT(  7,ijkQ,iPass) = IN(  1,ijkQ,iPass)*SPHMAT1_7      
    OUT(  7,ijkQ,iPass) = OUT(  7,ijkQ,iPass) + IN( 10,ijkQ,iPass)*SPHMAT1_7      
    OUT(  7,ijkQ,iPass) = OUT(  7,ijkQ,iPass) + IN( 16,ijkQ,iPass)*SPHMAT16_7     
    OUT(  8,ijkQ,iPass) = IN(  2,ijkQ,iPass)*SPHMAT1_7      
    OUT(  8,ijkQ,iPass) = OUT(  8,ijkQ,iPass) + IN( 11,ijkQ,iPass)*SPHMAT1_7      
    OUT(  8,ijkQ,iPass) = OUT(  8,ijkQ,iPass) + IN( 17,ijkQ,iPass)*SPHMAT16_7     
    OUT(  9,ijkQ,iPass) = IN(  3,ijkQ,iPass)*SPHMAT1_7      
    OUT(  9,ijkQ,iPass) = OUT(  9,ijkQ,iPass) + IN( 12,ijkQ,iPass)*SPHMAT1_7      
    OUT(  9,ijkQ,iPass) = OUT(  9,ijkQ,iPass) + IN( 18,ijkQ,iPass)*SPHMAT16_7     
    OUT( 10,ijkQ,iPass) = IN(  7,ijkQ,iPass)
    OUT( 11,ijkQ,iPass) = IN(  8,ijkQ,iPass)
    OUT( 12,ijkQ,iPass) = IN(  9,ijkQ,iPass)
    OUT( 13,ijkQ,iPass) = IN(  1,ijkQ,iPass)*SPHMAT1_13     
    OUT( 13,ijkQ,iPass) = OUT( 13,ijkQ,iPass) + IN( 10,ijkQ,iPass)*SPHMAT10_13    
    OUT( 14,ijkQ,iPass) = IN(  2,ijkQ,iPass)*SPHMAT1_13     
    OUT( 14,ijkQ,iPass) = OUT( 14,ijkQ,iPass) + IN( 11,ijkQ,iPass)*SPHMAT10_13    
    OUT( 15,ijkQ,iPass) = IN(  3,ijkQ,iPass)*SPHMAT1_13     
    OUT( 15,ijkQ,iPass) = OUT( 15,ijkQ,iPass) + IN( 12,ijkQ,iPass)*SPHMAT10_13    
   ENDDO
  ENDDO
end subroutine SphericalContractOBS1_maxAngP3_maxAngA1 
  
  
subroutine SphericalContractOBS1_maxAngP4_maxAngA2(ijkQcart,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: ijkQcart,nContPasses
  real(realk),intent(in)    :: IN( 36,ijkQcart,nContPasses)
  real(realk),intent(inout) :: OUT( 25,ijkQcart,nContPasses)
  integer :: iPass,ijkQ
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
  DO iPass=1,nContPasses
   DO ijkQ=1,ijkQcart
    OUT(  1,ijkQ,iPass) = IN(  8,ijkQ,iPass)
    OUT(  2,ijkQ,iPass) = IN( 11,ijkQ,iPass)
    OUT(  3,ijkQ,iPass) = IN(  7,ijkQ,iPass)*SPHMAT2_11     
    OUT(  3,ijkQ,iPass) = OUT(  3,ijkQ,iPass) + IN( 10,ijkQ,iPass)*SPHMAT2_11     
    OUT(  3,ijkQ,iPass) = OUT(  3,ijkQ,iPass) + IN( 12,ijkQ,iPass)*SPHMAT12_3     
    OUT(  4,ijkQ,iPass) = IN(  9,ijkQ,iPass)
    OUT(  5,ijkQ,iPass) = IN(  7,ijkQ,iPass)*SPHMAT2_21     
    OUT(  5,ijkQ,iPass) = OUT(  5,ijkQ,iPass) + IN( 10,ijkQ,iPass)*SPHMAT10_5     
    OUT(  6,ijkQ,iPass) = IN( 26,ijkQ,iPass)
    OUT(  7,ijkQ,iPass) = IN( 29,ijkQ,iPass)
    OUT(  8,ijkQ,iPass) = IN( 25,ijkQ,iPass)*SPHMAT2_11     
    OUT(  8,ijkQ,iPass) = OUT(  8,ijkQ,iPass) + IN( 28,ijkQ,iPass)*SPHMAT2_11     
    OUT(  8,ijkQ,iPass) = OUT(  8,ijkQ,iPass) + IN( 30,ijkQ,iPass)*SPHMAT12_3     
    OUT(  9,ijkQ,iPass) = IN( 27,ijkQ,iPass)
    OUT( 10,ijkQ,iPass) = IN( 25,ijkQ,iPass)*SPHMAT2_21     
    OUT( 10,ijkQ,iPass) = OUT( 10,ijkQ,iPass) + IN( 28,ijkQ,iPass)*SPHMAT10_5     
    OUT( 11,ijkQ,iPass) = IN(  2,ijkQ,iPass)*SPHMAT2_11     
    OUT( 11,ijkQ,iPass) = OUT( 11,ijkQ,iPass) + IN( 20,ijkQ,iPass)*SPHMAT2_11     
    OUT( 11,ijkQ,iPass) = OUT( 11,ijkQ,iPass) + IN( 32,ijkQ,iPass)*SPHMAT12_3     
    OUT( 12,ijkQ,iPass) = IN(  5,ijkQ,iPass)*SPHMAT2_11     
    OUT( 12,ijkQ,iPass) = OUT( 12,ijkQ,iPass) + IN( 23,ijkQ,iPass)*SPHMAT2_11     
    OUT( 12,ijkQ,iPass) = OUT( 12,ijkQ,iPass) + IN( 35,ijkQ,iPass)*SPHMAT12_3     
    OUT( 13,ijkQ,iPass) = IN(  1,ijkQ,iPass)*SPHMAT1_13     
    OUT( 13,ijkQ,iPass) = OUT( 13,ijkQ,iPass) + IN(  4,ijkQ,iPass)*SPHMAT1_13     
    OUT( 13,ijkQ,iPass) = OUT( 13,ijkQ,iPass) + IN(  6,ijkQ,iPass)*SPHMAT6_13     
    OUT( 13,ijkQ,iPass) = OUT( 13,ijkQ,iPass) + IN( 19,ijkQ,iPass)*SPHMAT1_13     
    OUT( 13,ijkQ,iPass) = OUT( 13,ijkQ,iPass) + IN( 22,ijkQ,iPass)*SPHMAT1_13     
    OUT( 13,ijkQ,iPass) = OUT( 13,ijkQ,iPass) + IN( 24,ijkQ,iPass)*SPHMAT6_13     
    OUT( 13,ijkQ,iPass) = OUT( 13,ijkQ,iPass) + IN( 31,ijkQ,iPass)*SPHMAT6_13     
    OUT( 13,ijkQ,iPass) = OUT( 13,ijkQ,iPass) + IN( 34,ijkQ,iPass)*SPHMAT6_13     
    OUT( 13,ijkQ,iPass) = OUT( 13,ijkQ,iPass) + IN( 36,ijkQ,iPass)*SPHMAT36_13    
    OUT( 14,ijkQ,iPass) = IN(  3,ijkQ,iPass)*SPHMAT2_11     
    OUT( 14,ijkQ,iPass) = OUT( 14,ijkQ,iPass) + IN( 21,ijkQ,iPass)*SPHMAT2_11     
    OUT( 14,ijkQ,iPass) = OUT( 14,ijkQ,iPass) + IN( 33,ijkQ,iPass)*SPHMAT12_3     
    OUT( 15,ijkQ,iPass) = IN(  1,ijkQ,iPass)*SPHMAT1_15     
    OUT( 15,ijkQ,iPass) = OUT( 15,ijkQ,iPass) + IN(  4,ijkQ,iPass)*SPHMAT4_15     
    OUT( 15,ijkQ,iPass) = OUT( 15,ijkQ,iPass) + IN( 19,ijkQ,iPass)*SPHMAT1_15     
    OUT( 15,ijkQ,iPass) = OUT( 15,ijkQ,iPass) + IN( 22,ijkQ,iPass)*SPHMAT4_15     
    OUT( 15,ijkQ,iPass) = OUT( 15,ijkQ,iPass) + IN( 31,ijkQ,iPass)*SPHMAT6_23     
    OUT( 15,ijkQ,iPass) = OUT( 15,ijkQ,iPass) + IN( 34,ijkQ,iPass)*SPHMAT2_11     
    OUT( 16,ijkQ,iPass) = IN( 14,ijkQ,iPass)
    OUT( 17,ijkQ,iPass) = IN( 17,ijkQ,iPass)
    OUT( 18,ijkQ,iPass) = IN( 13,ijkQ,iPass)*SPHMAT2_11     
    OUT( 18,ijkQ,iPass) = OUT( 18,ijkQ,iPass) + IN( 16,ijkQ,iPass)*SPHMAT2_11     
    OUT( 18,ijkQ,iPass) = OUT( 18,ijkQ,iPass) + IN( 18,ijkQ,iPass)*SPHMAT12_3     
    OUT( 19,ijkQ,iPass) = IN( 15,ijkQ,iPass)
    OUT( 20,ijkQ,iPass) = IN( 13,ijkQ,iPass)*SPHMAT2_21     
    OUT( 20,ijkQ,iPass) = OUT( 20,ijkQ,iPass) + IN( 16,ijkQ,iPass)*SPHMAT10_5     
    OUT( 21,ijkQ,iPass) = IN(  2,ijkQ,iPass)*SPHMAT2_21     
    OUT( 21,ijkQ,iPass) = OUT( 21,ijkQ,iPass) + IN( 20,ijkQ,iPass)*SPHMAT10_5     
    OUT( 22,ijkQ,iPass) = IN(  5,ijkQ,iPass)*SPHMAT2_21     
    OUT( 22,ijkQ,iPass) = OUT( 22,ijkQ,iPass) + IN( 23,ijkQ,iPass)*SPHMAT10_5     
    OUT( 23,ijkQ,iPass) = IN(  1,ijkQ,iPass)*SPHMAT1_15     
    OUT( 23,ijkQ,iPass) = OUT( 23,ijkQ,iPass) + IN(  4,ijkQ,iPass)*SPHMAT1_15     
    OUT( 23,ijkQ,iPass) = OUT( 23,ijkQ,iPass) + IN(  6,ijkQ,iPass)*SPHMAT6_23     
    OUT( 23,ijkQ,iPass) = OUT( 23,ijkQ,iPass) + IN( 19,ijkQ,iPass)*SPHMAT4_15     
    OUT( 23,ijkQ,iPass) = OUT( 23,ijkQ,iPass) + IN( 22,ijkQ,iPass)*SPHMAT4_15     
    OUT( 23,ijkQ,iPass) = OUT( 23,ijkQ,iPass) + IN( 24,ijkQ,iPass)*SPHMAT2_11     
    OUT( 24,ijkQ,iPass) = IN(  3,ijkQ,iPass)*SPHMAT2_21     
    OUT( 24,ijkQ,iPass) = OUT( 24,ijkQ,iPass) + IN( 21,ijkQ,iPass)*SPHMAT10_5     
    OUT( 25,ijkQ,iPass) = IN(  1,ijkQ,iPass)*SPHMAT1_25     
    OUT( 25,ijkQ,iPass) = OUT( 25,ijkQ,iPass) + IN(  4,ijkQ,iPass)*SPHMAT4_25     
    OUT( 25,ijkQ,iPass) = OUT( 25,ijkQ,iPass) + IN( 19,ijkQ,iPass)*SPHMAT4_25     
    OUT( 25,ijkQ,iPass) = OUT( 25,ijkQ,iPass) + IN( 22,ijkQ,iPass)*SPHMAT1_25     
   ENDDO
  ENDDO
end subroutine SphericalContractOBS1_maxAngP4_maxAngA2 
  
  
END MODULE AGC_OBS_Sphcontract1Mod
