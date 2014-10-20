module SPAGC_CPU_OBS_Sphcontract1Mod
!Automatic Generated Code (AGC) by runSphContractOBS1.f90 in tools directory
use IchorPrecisionMod
  
 CONTAINS
  
 subroutine SPSphericalContractOBS1_CPU_maxAngP2_maxAngA2(ijkQcart,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: ijkQcart,nContPasses
  real(reals),intent(in)    :: IN(  6,ijkQcart*nContPasses)
  real(reals),intent(inout) :: OUT(  5,ijkQcart*nContPasses)
  integer :: iP,ijkQ
  real(reals),parameter :: SPHMAT1_3      =   -2.8867513459481292E-01_reals
  real(reals),parameter :: SPHMAT1_5      =    5.0000000000000000E-01_reals
  real(reals),parameter :: SPHMAT2_1      =    1.0000000000000000E+00_reals
  real(reals),parameter :: SPHMAT4_5      =   -5.0000000000000000E-01_reals
  real(reals),parameter :: SPHMAT6_3      =    5.7735026918962584E-01_reals
!$OMP DO PRIVATE(iP)
  DO iP=1,ijkQcart*nContPasses
    OUT(1,iP) = IN(2,iP)
    OUT(2,iP) = IN(5,iP)
    OUT(3,iP) = IN(1,iP)*SPHMAT1_3 + IN(4,iP)*SPHMAT1_3 + IN(6,iP)*SPHMAT6_3
    OUT(4,iP) = IN(3,iP)
    OUT(5,iP) = IN(1,iP)*SPHMAT1_5 + IN(4,iP)*SPHMAT4_5
  ENDDO
!$OMP END DO
end subroutine SPSphericalContractOBS1_CPU_maxAngP2_maxAngA2 
  
  
 subroutine SPSphericalContractOBS1_CPU_maxAngP2_maxAngA0(ijkQcart,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: ijkQcart,nContPasses
  real(reals),intent(in)    :: IN(  6,ijkQcart*nContPasses)
  real(reals),intent(inout) :: OUT(  5,ijkQcart*nContPasses)
  integer :: iP,ijkQ
  real(reals),parameter :: SPHMAT1_3      =   -2.8867513459481292E-01_reals
  real(reals),parameter :: SPHMAT1_5      =    5.0000000000000000E-01_reals
  real(reals),parameter :: SPHMAT2_1      =    1.0000000000000000E+00_reals
  real(reals),parameter :: SPHMAT4_5      =   -5.0000000000000000E-01_reals
  real(reals),parameter :: SPHMAT6_3      =    5.7735026918962584E-01_reals
!$OMP DO PRIVATE(iP)
  DO iP=1,ijkQcart*nContPasses
    OUT(1,iP) = IN(2,iP)
    OUT(2,iP) = IN(5,iP)
    OUT(3,iP) = IN(1,iP)*SPHMAT1_3 + IN(4,iP)*SPHMAT1_3 + IN(6,iP)*SPHMAT6_3
    OUT(4,iP) = IN(3,iP)
    OUT(5,iP) = IN(1,iP)*SPHMAT1_5 + IN(4,iP)*SPHMAT4_5
  ENDDO
!$OMP END DO
end subroutine SPSphericalContractOBS1_CPU_maxAngP2_maxAngA0 
  
  
 subroutine SPSphericalContractOBS1_CPU_maxAngP3_maxAngA2(ijkQcart,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: ijkQcart,nContPasses
  real(reals),intent(in)    :: IN( 18,ijkQcart*nContPasses)
  real(reals),intent(inout) :: OUT( 15,ijkQcart*nContPasses)
  integer :: iP,ijkQ
  real(reals),parameter :: SPHMAT1_3      =   -2.8867513459481292E-01_reals
  real(reals),parameter :: SPHMAT1_5      =    5.0000000000000000E-01_reals
  real(reals),parameter :: SPHMAT2_1      =    1.0000000000000000E+00_reals
  real(reals),parameter :: SPHMAT4_5      =   -5.0000000000000000E-01_reals
  real(reals),parameter :: SPHMAT6_3      =    5.7735026918962584E-01_reals
!$OMP DO PRIVATE(iP)
  DO iP=1,ijkQcart*nContPasses
    OUT(1,iP) = IN(2,iP)
    OUT(2,iP) = IN(5,iP)
    OUT(3,iP) = IN(1,iP)*SPHMAT1_3 + IN(4,iP)*SPHMAT1_3 + IN(6,iP)*SPHMAT6_3
    OUT(4,iP) = IN(3,iP)
    OUT(5,iP) = IN(1,iP)*SPHMAT1_5 + IN(4,iP)*SPHMAT4_5
    OUT(6,iP) = IN(8,iP)
    OUT(7,iP) = IN(11,iP)
    OUT(8,iP) = IN(7,iP)*SPHMAT1_3 + IN(10,iP)*SPHMAT1_3 + IN(12,iP)*SPHMAT6_3
    OUT(9,iP) = IN(9,iP)
    OUT(10,iP) = IN(7,iP)*SPHMAT1_5 + IN(10,iP)*SPHMAT4_5
    OUT(11,iP) = IN(14,iP)
    OUT(12,iP) = IN(17,iP)
    OUT(13,iP) = IN(13,iP)*SPHMAT1_3 + IN(16,iP)*SPHMAT1_3 + IN(18,iP)*SPHMAT6_3
    OUT(14,iP) = IN(15,iP)
    OUT(15,iP) = IN(13,iP)*SPHMAT1_5 + IN(16,iP)*SPHMAT4_5
  ENDDO
!$OMP END DO
end subroutine SPSphericalContractOBS1_CPU_maxAngP3_maxAngA2 
  
  
 subroutine SPSphericalContractOBS1_CPU_maxAngP3_maxAngA1(ijkQcart,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: ijkQcart,nContPasses
  real(reals),intent(in)    :: IN( 18,ijkQcart*nContPasses)
  real(reals),intent(inout) :: OUT( 15,ijkQcart*nContPasses)
  integer :: iP,ijkQ
  real(reals),parameter :: SPHMAT1_7      =   -2.8867513459481292E-01_reals
  real(reals),parameter :: SPHMAT1_13     =    5.0000000000000000E-01_reals
  real(reals),parameter :: SPHMAT4_1      =    1.0000000000000000E+00_reals
  real(reals),parameter :: SPHMAT10_13    =   -5.0000000000000000E-01_reals
  real(reals),parameter :: SPHMAT16_7     =    5.7735026918962584E-01_reals
!$OMP DO PRIVATE(iP)
  DO iP=1,ijkQcart*nContPasses
    OUT(1,iP) = IN(4,iP)
    OUT(2,iP) = IN(5,iP)
    OUT(3,iP) = IN(6,iP)
    OUT(4,iP) = IN(13,iP)
    OUT(5,iP) = IN(14,iP)
    OUT(6,iP) = IN(15,iP)
    OUT(7,iP) = IN(1,iP)*SPHMAT1_7 + IN(10,iP)*SPHMAT1_7 + IN(16,iP)*SPHMAT16_7
    OUT(8,iP) = IN(2,iP)*SPHMAT1_7 + IN(11,iP)*SPHMAT1_7 + IN(17,iP)*SPHMAT16_7
    OUT(9,iP) = IN(3,iP)*SPHMAT1_7 + IN(12,iP)*SPHMAT1_7 + IN(18,iP)*SPHMAT16_7
    OUT(10,iP) = IN(7,iP)
    OUT(11,iP) = IN(8,iP)
    OUT(12,iP) = IN(9,iP)
    OUT(13,iP) = IN(1,iP)*SPHMAT1_13 + IN(10,iP)*SPHMAT10_13
    OUT(14,iP) = IN(2,iP)*SPHMAT1_13 + IN(11,iP)*SPHMAT10_13
    OUT(15,iP) = IN(3,iP)*SPHMAT1_13 + IN(12,iP)*SPHMAT10_13
  ENDDO
!$OMP END DO
end subroutine SPSphericalContractOBS1_CPU_maxAngP3_maxAngA1 
  
  
 subroutine SPSphericalContractOBS1_CPU_maxAngP4_maxAngA2(ijkQcart,nContPasses,IN,OUT)
  implicit none
  integer,intent(in)        :: ijkQcart,nContPasses
  real(reals),intent(in)    :: IN( 36,ijkQcart*nContPasses)
  real(reals),intent(inout) :: OUT( 25,ijkQcart*nContPasses)
  integer :: iP,ijkQ
  real(reals),parameter :: SPHMAT1_13     =    8.3333333333333356E-02_reals
  real(reals),parameter :: SPHMAT1_15     =   -1.4433756729740646E-01_reals
  real(reals),parameter :: SPHMAT1_25     =    2.5000000000000000E-01_reals
  real(reals),parameter :: SPHMAT2_11     =   -2.8867513459481292E-01_reals
  real(reals),parameter :: SPHMAT2_21     =    5.0000000000000000E-01_reals
  real(reals),parameter :: SPHMAT4_15     =    1.4433756729740646E-01_reals
  real(reals),parameter :: SPHMAT4_25     =   -2.5000000000000000E-01_reals
  real(reals),parameter :: SPHMAT6_13     =   -1.6666666666666671E-01_reals
  real(reals),parameter :: SPHMAT6_23     =    2.8867513459481292E-01_reals
  real(reals),parameter :: SPHMAT8_1      =    1.0000000000000000E+00_reals
  real(reals),parameter :: SPHMAT10_5     =   -5.0000000000000000E-01_reals
  real(reals),parameter :: SPHMAT12_3     =    5.7735026918962584E-01_reals
  real(reals),parameter :: SPHMAT36_13    =    3.3333333333333343E-01_reals
!$OMP DO PRIVATE(iP)
  DO iP=1,ijkQcart*nContPasses
    OUT(1,iP) = IN(8,iP)
    OUT(2,iP) = IN(11,iP)
    OUT(3,iP) = IN(7,iP)*SPHMAT2_11 + IN(10,iP)*SPHMAT2_11 + IN(12,iP)*SPHMAT12_3
    OUT(4,iP) = IN(9,iP)
    OUT(5,iP) = IN(7,iP)*SPHMAT2_21 + IN(10,iP)*SPHMAT10_5
    OUT(6,iP) = IN(26,iP)
    OUT(7,iP) = IN(29,iP)
    OUT(8,iP) = IN(25,iP)*SPHMAT2_11 + IN(28,iP)*SPHMAT2_11 + IN(30,iP)*SPHMAT12_3
    OUT(9,iP) = IN(27,iP)
    OUT(10,iP) = IN(25,iP)*SPHMAT2_21 + IN(28,iP)*SPHMAT10_5
    OUT(11,iP) = IN(2,iP)*SPHMAT2_11 + IN(20,iP)*SPHMAT2_11 + IN(32,iP)*SPHMAT12_3
    OUT(12,iP) = IN(5,iP)*SPHMAT2_11 + IN(23,iP)*SPHMAT2_11 + IN(35,iP)*SPHMAT12_3
    OUT(13,iP) = IN(1,iP)*SPHMAT1_13 + IN(4,iP)*SPHMAT1_13 + IN(6,iP)*SPHMAT6_13 + IN(19,iP)*SPHMAT1_13 &
               & + IN(22,iP)*SPHMAT1_13 + IN(24,iP)*SPHMAT6_13 + IN(31,iP)*SPHMAT6_13 + IN(34,iP)*SPHMAT6_13 &
               & + IN(36,iP)*SPHMAT36_13
    OUT(14,iP) = IN(3,iP)*SPHMAT2_11 + IN(21,iP)*SPHMAT2_11 + IN(33,iP)*SPHMAT12_3
    OUT(15,iP) = IN(1,iP)*SPHMAT1_15 + IN(4,iP)*SPHMAT4_15 + IN(19,iP)*SPHMAT1_15 + IN(22,iP)*SPHMAT4_15 &
               & + IN(31,iP)*SPHMAT6_23 + IN(34,iP)*SPHMAT2_11
    OUT(16,iP) = IN(14,iP)
    OUT(17,iP) = IN(17,iP)
    OUT(18,iP) = IN(13,iP)*SPHMAT2_11 + IN(16,iP)*SPHMAT2_11 + IN(18,iP)*SPHMAT12_3
    OUT(19,iP) = IN(15,iP)
    OUT(20,iP) = IN(13,iP)*SPHMAT2_21 + IN(16,iP)*SPHMAT10_5
    OUT(21,iP) = IN(2,iP)*SPHMAT2_21 + IN(20,iP)*SPHMAT10_5
    OUT(22,iP) = IN(5,iP)*SPHMAT2_21 + IN(23,iP)*SPHMAT10_5
    OUT(23,iP) = IN(1,iP)*SPHMAT1_15 + IN(4,iP)*SPHMAT1_15 + IN(6,iP)*SPHMAT6_23 + IN(19,iP)*SPHMAT4_15 &
               & + IN(22,iP)*SPHMAT4_15 + IN(24,iP)*SPHMAT2_11
    OUT(24,iP) = IN(3,iP)*SPHMAT2_21 + IN(21,iP)*SPHMAT10_5
    OUT(25,iP) = IN(1,iP)*SPHMAT1_25 + IN(4,iP)*SPHMAT4_25 + IN(19,iP)*SPHMAT4_25 + IN(22,iP)*SPHMAT1_25
  ENDDO
!$OMP END DO
end subroutine SPSphericalContractOBS1_CPU_maxAngP4_maxAngA2 
  
  
END module SPAGC_CPU_OBS_Sphcontract1Mod
