MODULE IchorParametersMod
  use IchorCommonMod
!> Spherical Specification
Integer,parameter :: SphericalParam = 1 !spherical harmonic basis set
!Job Specification
integer,parameter :: IchorJobEri = 1     !4 center 2 electronic repulsion integrals
integer,parameter :: IchorJobLink = 2    !Linear scaling exchange (LinK) 
integer,parameter :: IchorJobMoTrans = 3 !Calc Mo Integrals directly
                                         !MO transform the AO integral block 
                                         !directly to MO (nMO1,nMO2,nMO3,nMO4)
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
Integer,parameter :: MaxSpecialAngmom = 2
!> IchorOperatorSpec
integer,parameter :: CoulombOperator=1
integer,parameter :: GGemOperator   =2    ! The Gaussian geminal operator g
integer,parameter :: GGemCouOperator=3    ! The Gaussian geminal divided by the Coulomb operator g/r12
integer,parameter :: GGemGrdOperator=4    ! The double commutator [[T,g],g]
integer,parameter :: GGemSqOperator =5    ! The Gaussian geminal operator squared g^2
real(realk) :: IchorGPUmaxmem

!public:: SphericalParam, IcorJobEri, IcorInputNoInput, IchorParNone,&
!     & IchorScreen, IchorScreenNone, IchorDebugNone, IchorAlgoOS, IchorPermuteTTT,&
!     & IchorPermuteFFT, IchorPermuteTTF, IchorPermuteTFF, IchorPermuteFTF, Ichor!PermuteFFF, IchorNofilestorage
!private
#ifdef VAR_OPENACC
  integer,parameter :: maxnAsyncHandles=16
#else
  integer,parameter :: maxnAsyncHandles=1
#endif
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

subroutine determine_OpereratorSpec(Oper,Operparam)
  implicit none
  Character(len=7),intent(in)     :: Oper
  integer,intent(inout)     :: Operparam

  SELECT CASE(Oper)
  CASE('Coulomb') 
     !standard Coulomb Operator 
     operparam = CoulombOperator
  CASE('GGem   ') 
     ! The Gaussian geminal operator g
     operparam = GGemOperator
  CASE('GGemCou') 
     ! The Gaussian geminal divided by the Coulomb operator g/r12
     operparam = GGemCouOperator
  CASE('GGemGrd') 
     ! The double commutator [[T,g],g]
     operparam = GGemGrdOperator
  CASE('GGemSq ') 
     ! The Gaussian geminal operator squared g^2
     operparam = GGemSqOperator
  CASE DEFAULT
     WRITE(6,'(1X,2A)') 'Unknown Operator =',Oper
     CALL ICHORQUIT('Param_oper_paramfromString',-1)
  END SELECT

end subroutine Determine_OpereratorSpec

subroutine Determine_Operator_From_OpereratorSpec(Oper,Operparam)
  implicit none
  Character(len=7),intent(inout)     :: Oper
  integer,intent(in)     :: Operparam

  SELECT CASE(Operparam)
  CASE(CoulombOperator) 
     oper = 'Coulomb'
  CASE(GGemOperator) 
     oper = 'GGem   '
  CASE(GGemCouOperator) 
     oper = 'GGemCou'
  CASE(GGemGrdOperator) 
     oper = 'GGemGrd'
  CASE(GGemSqOperator) 
     oper = 'GGemSq '
  CASE DEFAULT
     WRITE(6,'(1X,2A)') 'Unknown Operator parameter =',Operparam
     CALL ICHORQUIT('Unknown Operator Param_oper_Stringfromparam',-1)
  END SELECT

end subroutine Determine_Operator_From_OpereratorSpec

END MODULE IchorParametersMod
