!> @file 
!> Contains the precision specifications
MODULE lsparameters
  use precision
! THESE ARE STRING SPECIFIERS FOR THE AOs
  integer,parameter :: AORegular = 1
  integer,parameter :: AOEmpty = 2
  integer,parameter :: AOdfAux = 3
  integer,parameter :: AONuclear = 4
  integer,parameter :: AOpCharge = 5
  integer,parameter :: AOdfCABS = 6
  integer,parameter :: AOdfJK = 7
  integer,parameter :: AOVAL = 8
  integer,parameter :: AOelField = 9
  integer,parameter :: AOadmm = 10  !ADMM basis
  integer,parameter :: AONuclearSpec = 11 !single Nuclei
! THESE ARE STRING SPECIFIERS FOR THE Operator
  integer,parameter :: CoulombOperator = 1
  integer,parameter :: OverlapOperator = 2
  integer,parameter :: KineticOperator = 3
  integer,parameter :: NucpotOperator = 4
  integer,parameter :: CarmomOperator = 5
  integer,parameter :: ANGLONOperator = 6
  integer,parameter :: ANGMOMOperator = 7
  integer,parameter :: DIPVELOperator = 8
  integer,parameter :: DIPLENOperator = 9
  integer,parameter :: ROTSTROperator = 10
  integer,parameter :: THETAOperator = 11
  integer,parameter :: LONMOM1Operator = 12
  integer,parameter :: LONMOM2Operator = 13
  integer,parameter :: PSOOperator = 14
  integer,parameter :: NSTNOLOperator = 15
  integer,parameter :: NSTLONOperator = 16
  integer,parameter :: ELPOTOperator = 17
  integer,parameter :: NucleiOperator = 18
  integer,parameter :: ErfcOperator = 19
  integer,parameter :: ErfOperator = 20
  integer,parameter :: CAMOperator = 21
  integer,parameter :: MulmomOperator = 22
  integer,parameter :: GGemOperator = 23
  integer,parameter :: GGemCouOperator = 24
  integer,parameter :: GGemGrdOperator = 25
  integer,parameter :: MAGMOMOperator = 26
  integer,parameter :: NSTOperator = 27
  integer,parameter :: LONMOMOperator = 28
  integer,parameter :: DCM1Operator = 29
  integer,parameter :: DCM2Operator = 30
! THESE ARE STRING SPECIFIERS FOR THE integralType
  integer,parameter :: ContractedInttype = 1
  integer,parameter :: PrimitiveInttype = 2
! THESE ARE STRING SPECIFIERS FOR THE SPECIFICATION
  integer,parameter :: RegularSpec         = 1 
  integer,parameter :: GradientSpec        = 2 
  integer,parameter :: MagDerivSpec        = 3 
  integer,parameter :: MagGradSpec         = 4
  integer,parameter :: MagDerivRSpec       = 5 
  integer,parameter :: MagDerivLSpec       = 6
  integer,parameter :: pChargeSpec         = 7
  integer,parameter :: EcontribSpec        = 8
  integer,parameter :: GeoDerivSpec        = 9
  integer,parameter :: GeoDerivCoulombSpec = 10
  integer,parameter :: GeoDerivLHSSpec     = 11
  integer,parameter :: GeoDerivRHSSpec     = 12

! THESE ARE MPI JOB SPECIFIERS 
  integer,parameter :: MATRIXTY                     =  1
  integer,parameter :: LSGETINT                     =  2
  integer,parameter :: LSJENGIN                     =  3
  integer,parameter :: LSLINK                       =  4
  integer,parameter :: DFTSETFU                     =  5
  integer,parameter :: LSMPI_IIDFTKSM               =  6
  integer,parameter :: IIDFTGEO                     =  7
  integer,parameter :: IIDFTLIN                     =  8
  integer,parameter :: IIDFTQRS                     =  9
  integer,parameter :: IIDFTMAG                     = 10
  integer,parameter :: IIDFTMAL                     = 11
  integer,parameter :: IIDFTGKS                     = 12
  integer,parameter :: IIDFTGLR                     = 13
  integer,parameter :: MP2INAMP                     = 14
  integer,parameter :: GRIDINIT                     = 15
  integer,parameter :: GRIDEXIT                     = 16
  integer,parameter :: PDMSLAVE                     = 17
  integer,parameter :: LSMPIQUIT                    = 18    
  integer,parameter :: LSMPIPRINTINFO               = 19
  integer,parameter :: CCSDDATA                     = 20
  integer,parameter :: GROUPINIT                    = 21
  integer,parameter :: DECDRIVER                    = 22
  integer,parameter :: DEFAULTGROUPS                = 23
  integer,parameter :: IISCREEN                     = 24
  integer,parameter :: IISCREENFREE                 = 25
  integer,parameter :: IISCREENINIT                 = 26
  integer,parameter :: QUITMOREJOBS                 = 27
  integer,parameter :: QUITNOMOREJOBS               = 28
  integer,parameter :: LSMPITEST                    = 29
  integer,parameter :: ARRAYTEST                    = 30
  integer,parameter :: PDMA4SLV                     = 31
  integer,parameter :: LSMPI_IIDFTKSME              = 32
#ifdef MOD_UNRELEASED
  integer,parameter :: CCSDPTSLAVE                  = 33
#endif
  integer,parameter :: CCSDSLV4E2                   = 34
  integer,parameter :: DFTADDFU                     = 35
  integer,parameter :: LSMPI_IIDFTABSVALOVERLAP     = 36
  integer,parameter :: GIVE_BIRTH                   = 37
  integer,parameter :: LSPDM_GIVE_BIRTH             = 38
  integer,parameter :: SLAVES_SHUT_DOWN_CHILD       = 39
  integer,parameter :: LSPDM_SLAVES_SHUT_DOWN_CHILD = 40
  integer,parameter :: CHILD_SHUT_DOWN              = 41
  integer,parameter :: JOB_LSPDM_INIT_GLOBAL_BUFFER = 42
  integer,parameter :: JOB_LSPDM_FREE_GLOBAL_BUFFER = 43
  integer,parameter :: CCGETGMO                     = 44
  integer,parameter :: RPAGETRESIDUAL               = 45
  integer,parameter :: MOCCSDDATA                   = 46
  integer,parameter :: MO_INTEGRAL_SIMPLE           = 47
  integer,parameter :: DEC_SETTING_TO_SLAVES        = 48
  integer,parameter :: INITSLAVETIME                = 49
  integer,parameter :: GETSLAVETIME                 = 50
  integer,parameter :: RIMP2INAMP                   = 51
  integer,parameter :: SIMPLE_MP2_PAR               = 52
  integer,parameter :: RPAGETFOCK                   = 53
  integer,parameter :: SET_SPLIT_MPI_MSG            = 54
  integer,parameter :: SET_MAX_SIZE_ONE_SIDED       = 55
  integer,parameter :: RIMP2FULL                    = 56
  integer,parameter :: SET_GPUMAXMEM                = 57
  integer,parameter :: SET_TENSOR_BACKEND_TRUE      = 58
  integer,parameter :: SET_TENSOR_DEBUG_TRUE        = 59
  integer,parameter :: F12_INTEGRAL_CALCULATION     = 60
  integer,parameter :: CANONMP2FULL                 = 61
  integer,parameter :: PDMMGRIDINIT                 = 62
  integer,parameter :: PDMMGRIDEXIT                 = 63
  integer,parameter :: PDMMSLAVE                    = 64
  integer,parameter :: MATRIXTY2                    = 65
  integer,parameter :: SET_TENSOR_ALWAYS_SYNC_TRUE  = 66
  integer,parameter :: INIT_BG_BUF                  = 67
  integer,parameter :: FREE_BG_BUF                  = 68
  integer,parameter :: CHANGE_BG_BUF                = 69
  integer,parameter :: LSTHCRIMP2INAMP              = 70
  integer,parameter :: LSTHCRIMP2FULL               = 71

! s
  integer,parameter :: SymFromTriangularPostprocess=1
  integer,parameter :: SymmetricPostprocess=2
  integer,parameter :: AntiSymmetricPostprocess=3

  real(realk) :: GPUMAXMEM
save
INTEGER :: AORdefault
INTEGER :: AODFdefault
CONTAINS
SUBROUTINE init_AO_parameters()
  implicit none
  AORdefault = AORegular
  AODFdefault = AOdfAux
END SUBROUTINE init_AO_parameters

SUBROUTINE set_default_AOs(newAORegular,newAOdfAux)
  implicit none
  integer :: newAORegular,newAOdfAux
  AORdefault = newAORegular
  AODFdefault = newAOdfAux
END SUBROUTINE set_default_AOs

SUBROUTINE get_default_AOs(oldAORegular,oldAOdfAux)
  implicit none
  integer :: oldAORegular,oldAOdfAux
  oldAORegular = AORdefault 
  oldAOdfAux = AODFdefault  
END SUBROUTINE get_default_AOs

subroutine param_oper_paramfromString(Oper,Operparam)
  implicit none
  Character(len=7),intent(in)     :: Oper
  integer,intent(inout)     :: Operparam

  SELECT CASE(Oper)
  CASE('Coulomb') 
     operparam = CoulombOperator
  CASE('Overlap') 
     operparam = OverlapOperator
  CASE('Kinetic') 
     operparam = KineticOperator
  CASE('Nucpot ') 
     operparam = NucpotOperator
  CASE('CARMOM ') 
     operparam = CARMOMOperator
  CASE('ANGLON ') 
     operparam = ANGLONOperator
  CASE('ANGMOM ') 
     operparam = ANGMOMOperator
  CASE('DIPVEL ') 
     operparam = DIPVELOperator
  CASE('DIPLEN ') 
     operparam = DIPLENOperator
  CASE('ROTSTR ') 
     operparam = ROTSTROperator
  CASE('THETA  ') 
     operparam = THETAOperator
  CASE('LONMOM1') 
     operparam = LONMOM1Operator
  CASE('LONMOM2') 
     operparam = LONMOM2Operator
  CASE('PSO    ') 
     operparam = PSOOperator
  CASE('NSTNOL ') 
     operparam = NSTNOLOperator
  CASE('NSTLON ') 
     operparam = NSTLONOperator
  CASE('1ELPOT ') 
     operparam = ELPOTOperator
  CASE('Nuclei ') 
     operparam = NucleiOperator
  CASE('Erfc   ') 
     operparam = ErfcOperator
  CASE('Erf    ') 
     operparam = ErfOperator
  CASE('CAM    ') 
     operparam = CAMOperator
  CASE('GGem   ') 
     operparam = GGemOperator
  CASE('GGemCou') 
     operparam = GGemCouOperator
  CASE('GGemGrd') 
     operparam = GGemGrdOperator
  CASE('MAGMOM ') 
     operparam = MAGMOMOperator
  CASE('NST    ') 
     operparam = NSTOperator
  CASE('LONMOM ') 
     operparam = LONMOMOperator
  CASE('D-CM1  ') 
     operparam = DCM1Operator
  CASE('D-CM2  ') 
     operparam = DCM2Operator
  CASE DEFAULT
     WRITE(6,'(1X,2A)') 'Unknown Operator =',Oper
     CALL LSQUIT('Param_oper_paramfromString',-1)
  END SELECT

end subroutine Param_oper_paramfromString

subroutine param_oper_Stringfromparam(Oper,Operparam)
  implicit none
  Character(len=7),intent(inout)     :: Oper
  integer,intent(in)     :: Operparam

  SELECT CASE(Operparam)
  CASE(CoulombOperator) 
     oper = 'Coulomb'
  CASE(OverlapOperator) 
     oper = 'Overlap'
  CASE(KineticOperator) 
     oper = 'Kinetic'
  CASE(NucpotOperator)  
     oper = 'Nucpot '
  CASE(CARMOMOperator)  
     oper = 'CARMOM '
  CASE(ANGLONOperator)  
     oper = 'ANGLON '
  CASE(ANGMOMOperator)  
     oper = 'ANGMOM '
  CASE(DIPVELOperator)  
     oper = 'DIPVEL '
  CASE(DIPLENOperator)  
     oper = 'DIPLEN '
  CASE(ROTSTROperator)  
     oper = 'ROTSTR '
  CASE(THETAOperator)   
     oper = 'THETA  '
  CASE(LONMOM1Operator) 
     oper = 'LONMOM1'
  CASE(LONMOM2Operator) 
     oper = 'LONMOM2'
  CASE(PSOOperator)     
     oper = 'PSO    '
  CASE(NSTNOLOperator)  
     oper = 'NSTNOL '
  CASE(NSTLONOperator)  
     oper = 'NSTLON '
  CASE(ELPOTOperator)   
     oper = '1ELPOT '
  CASE(NucleiOperator)  
     oper = 'Nuclei '
  CASE(ErfcOperator)    
     oper = 'Erfc   '
  CASE(ErfOperator)     
     oper = 'Erf    '
  CASE(CAMOperator)     
     oper = 'CAM    '
  CASE(GGemOperator) 
     oper = 'GGem   '
  CASE(GGemCouOperator) 
     oper = 'GGemCou'
  CASE(GGemGrdOperator) 
     oper = 'GGemGrd'
  CASE(MAGMOMOperator) 
     oper = 'MAGMOM '
  CASE(NSTOperator) 
     oper = 'NST    '
  CASE(LONMOMOperator) 
     oper = 'LONMOM '
  CASE(DCM1Operator) 
     oper = 'D-CM1  '
  CASE(DCM2Operator) 
     oper= 'D-CM2  '
  CASE DEFAULT
     WRITE(6,'(1X,2A)') 'Unknown Operator parameter =',Operparam
     CALL LSQUIT('Unknown Operator Param_oper_Stringfromparam',-1)
  END SELECT

end subroutine Param_oper_Stringfromparam

subroutine param_AO_Stringfromparam(AO1,AO1param)
  implicit none
  Character(len=8),intent(inout)     :: AO1
  integer,intent(in)     :: AO1param
  SELECT CASE(AO1param)
  CASE(AORegular) 
     AO1 = 'Regular '
  CASE(AOEmpty)   
     AO1 = 'DF-Aux  '
  CASE(AOdfAux)   
     AO1 = 'Empty   '
  CASE(AONuclear) 
     AO1 = 'Nuclear '
  CASE(AOpCharge) 
     AO1 = 'pCharge '
  CASE(AOdfCABS) 
     AO1 = 'CABS    '
  CASE(AOdfJK) 
     AO1 = 'JKAUX   '
  CASE(AOVAL) 
     AO1 = 'VALENCE '
  CASE(AOelField) 
     AO1 = 'elField '
  CASE(AOadmm) 
     AO1 = 'ADMM    '
  CASE DEFAULT
     WRITE(6,'(1X,2A)') 'Unknown AO string =',AO1param
     CALL LSQUIT('Unknown Operator II_determineOperatorparameter',-1)
  END SELECT

end subroutine Param_AO_Stringfromparam

subroutine param_inttype_Stringfromparam(inttype,inttypeparam)
  implicit none
  Character(len=10),intent(inout)     :: inttype
  integer,intent(in)               :: inttypeparam

  SELECT CASE(inttypeparam)
  CASE(ContractedInttype) 
     inttype='Contracted'
  CASE(PrimitiveInttype)  
     inttype='Primitive '
  CASE DEFAULT
     WRITE(6,'(1X,2A)') 'Unknown inttype parameter =',inttypeparam
     CALL LSQUIT('Unknown Operator Param_oper_Stringfromparam',-1)
  END SELECT

end subroutine Param_inttype_Stringfromparam

END MODULE Lsparameters

