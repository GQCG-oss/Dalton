MODULE IchorParametersModule
!> Spherical Specification
Integer,parameter :: SphericalParam = 1 !spherical harmonic basis set
!Job Specification
integer,parameter :: IchorJobEri = 1     !4 center 2 electronic repulsion integrals
integer,parameter :: IchorJobLink = 2    !Linear scaling exchange (LinK) 
!Input Spec
integer,parameter :: IchorInputNoInput = 1 !no input in inputstorage (no Density matrix)
integer,parameter :: IchorInputDmat = 2    !RHS Density matrix stored in inputStorage
!Parallelization specification 
integer,parameter :: IchorParNone = 1     !no parallelization
!> Screening specification 
integer,parameter :: IchorScreen = 1      !default screening including Cauchy-Schwarz screening, OD and QQR
!> Screening specification 
integer,parameter :: IchorScreenCS = 2    !Only Cauchy-Schwarz (CS) screening
integer,parameter :: IchorScreenOD = 3    !only Overlap Density (OD) screening
integer,parameter :: IchorScreenQQR = 4   !only QQR screening
integer,parameter :: IchorScreenCSOD = 5  !CS and OD screening
integer,parameter :: IchorScreenNone = 6  !no screening
!> Debug info specification - The print Level or IchorDebugNone
integer,parameter :: IchorDebugNone = 0
!> Integral Algorithm specification 
!> IchorAlgoSpec = IchorAlgoOS = 1 means Obara-Saika (Head-Gordon Pople)
Integer,parameter :: IchorAlgoOS = 1
!> filestorageIdentifier = Logical Unit number .OR. IchorNofilestorage
Integer,parameter :: IchorNofilestorage = 0
!> IchorPermuteSpec
Integer,parameter :: IchorPermuteTTT = 1 !(SameLHSaos=.TRUE. , SameRHSaos=.TRUE. , SameODs=.TRUE. ) 
Integer,parameter :: IchorPermuteFFT = 2 !(SameLHSaos=.FALSE., SameRHSaos=.FALSE., SameODs=.TRUE. ) 
Integer,parameter :: IchorPermuteTTF = 3 !(SameLHSaos=.TRUE. , SameRHSaos=.TRUE. , SameODs=.FALSE.) 
Integer,parameter :: IchorPermuteTFF = 4 !(SameLHSaos=.TRUE. , SameRHSaos=.FALSE., SameODs=.FALSE.) 
Integer,parameter :: IchorPermuteFTF = 5 !(SameLHSaos=.FALSE., SameRHSaos=.TRUE. , SameODs=.FALSE.) 
Integer,parameter :: IchorPermuteFFF = 6 !(SameLHSaos=.FALSE., SameRHSaos=.FALSE., SameODs=.FALSE.) 

!public:: SphericalParam, IcorJobEri, IcorInputNoInput, IchorParNone,&
!     & IchorScreen, IchorScreenNone, IchorDebugNone, IchorAlgoOS, IchorPermuteTTT,&
!     & IchorPermuteFFT, IchorPermuteTTF, IchorPermuteTFF, IchorPermuteFTF, Ichor!PermuteFFF, IchorNofilestorage
!private
CONTAINS
subroutine determineScreening(IchorScreenSpec,CSscreen,ODscreen,QQRscreen)
  implicit none
  integer,intent(in) :: IchorScreenSpec
  logical,intent(inout) :: CSscreen,ODscreen,QQRscreen
  SELECT CASE(IchorScreenSpec)          
  CASE(IchorScreen)
     CSscreen = .TRUE. 
     ODscreen = .TRUE.
     QQRscreen = .TRUE.
  CASE(IchorScreenCS)
     CSscreen = .TRUE. 
     ODscreen = .FALSE.
     QQRscreen = .FALSE.
  CASE(IchorScreenOD)
     CSscreen = .FALSE.
     ODscreen = .TRUE.
     QQRscreen = .FALSE.
  CASE(IchorScreenQQR)
     CSscreen = .FALSE. 
     ODscreen = .FALSE.
     QQRscreen = .TRUE.
  CASE(IchorScreenCSOD)
     CSscreen = .TRUE. 
     ODscreen = .TRUE.
     QQRscreen = .FALSE.
  CASE(IchorScreenNone)
     CSscreen = .FALSE.
     ODscreen = .FALSE.
     QQRscreen = .FALSE.          
  CASE DEFAULT
     print*,'IchorScreenSpec = ',IchorScreenSpec,' not recognized in determineScreening'
     print*,'assume IchorScreen = ',IchorScreen
     CSscreen = .TRUE. 
     ODscreen = .TRUE.
     QQRscreen = .TRUE.
  END SELECT
END subroutine determineScreening

END MODULE IchorParametersModule
