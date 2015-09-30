!> @file 
!> Contains integralinput structure which contain all info required to do an integral evaluation.  
MODULE integral_type
use precision
use AO_typetype
use lstensor_typetype
use f12_module
!*****************************************
!*
!* OBJECT CONTAINING INFORMATION ABOUT
!* THE MAIN-INTEGRAL SETTINGS OBJECT
!*   - used in ThermiteDriver.f90
!*
!*****************************************
TYPE INTEGRALINPUT
REAL(REALK) :: Integralthreshold
TYPE(AOITEMPOINTER)   :: AO(4)
Integer               :: AOdim(4)
Integer               :: currentFragment(4)
type(lstensor),pointer :: LST_DLHS
type(lstensor),pointer :: LST_DRHS
type(lstensor),pointer :: LST_GAB_LHS
type(lstensor),pointer :: LST_GAB_RHS
logical               :: noOMP
logical               :: GAB_LHSusePointer
logical               :: GAB_RHSusePointer
logical               :: uselst_DRHS
logical               :: uselst_DLHS
!Real(realk),pointer   :: GAB_LHS(:,:)
!Real(realk),pointer   :: GAB_RHS(:,:)
!TYPE(LSMATRIX),pointer  :: DMAT_LHS
!TYPE(LSMATRIX),pointer  :: DMAT_RHS
!LOGICAL       :: LHS_DMAT,RHS_DMAT
INTEGER       :: FTUVmaxprim
INTEGER       :: LU_MMDATA,MMINDEX
INTEGER       :: LU_MMDATR
INTEGER       :: LU_MMDADA, LU_MMDADR
INTEGER       :: NDMAT_LHS,NDMAT_RHS
LOGICAL       :: DRHS_SYM,DLHS_SYM
!INTEGER       :: NDIM_LHS(2),NDIM_RHS(2)
LOGICAL       :: sameLHSaos,sameRHSaos,sameODs
INTEGER       :: CENTERS,ndim
LOGICAL       :: DO_PASSES
integer       :: maxPasses
LOGICAL       :: HIGH_RJ000_ACCURACY
integer       :: MMunique_ID1
integer       :: MMunique_ID2
integer       :: MMstartA
integer       :: MMstartB
LOGICAL       :: DO_COULOMB,DO_EXCHANGE,DO_FMM,DO_MULMOM,DO_PROP
LOGICAL       :: DO_FOCK,DO_ENERGY,DO_JENGINE,DO_DAJENGINE,DO_DACOULOMB,DO_LINK,DO_DALINK
INTEGER       :: DASCREEN_THRLOG
LOGICAL       :: CS_int !Cauchy-Schwarz integrals
LOGICAL       :: PS_int !Primitive Cauchy-Schwarz integrals
LOGICAL       :: setETUVoutside !Sets up ETUV-tensor before integral-loop
LOGICAL       :: sphericalEcoeff !ETUV-tensor in solid-harmonical form
LOGICAL       :: hermiteEcoeff
!Cauchy-Schwarz screening
LOGICAL       :: CS_SCREEN 
LOGICAL       :: ContAng !Indicates the orbital ordering to be contracted first angular second.
                         !The deafult is angular components first then contracted (only applicable 
                         !to generally contraced basis sets)
Real(realk)   :: CS_THRESHOLD
integer(kind=short)  :: CS_THRLOG
!Cauchy-Schwarz screening
LOGICAL       :: PS_SCREEN 
Real(realk)   :: PS_THRESHOLD
integer(kind=short)  :: PS_THRLOG
!Cauchy-Schwarz screening
LOGICAL       :: PARI_SCREEN 
Real(realk)   :: PARI_THRESHOLD
!Overlap-extent screening (for overlap integrals)
LOGICAL       :: OE_SCREEN 
Real(realk)   :: OE_THRESHOLD
!Screen OD-batches by AO-batch extent
LOGICAL       :: OD_SCREEN 
Real(realk)   :: OD_THRESHOLD
!Multipole Based integral estimate (MBIE) screening
LOGICAL       :: MBIE_SCREEN 
LOGICAL       :: MBIE_INT
!Screen integrals by their non-classical extent
LOGICAL       :: NonClassical_SCREEN 
Real(realk)   :: MM_SCREENTHR != 1E-10_realk
Integer       :: MM_TLMAX     != 8
Integer       :: MM_LMAX      != 20
Logical       :: MM_NOONE     != .false.
Logical       :: MM_NOSCREEN  != .false.
Logical       :: DO_MMGRD     != .true.
!Primitive-Screening maximum element of prim GAB matrix
integer(kind=short) :: CS_MAXELM_LHS
integer(kind=short) :: CS_MAXELM_RHS
integer(kind=short) :: PS_MAXELM_LHS
integer(kind=short) :: PS_MAXELM_RHS
Real(realk)   :: CoulombFactor
Real(realk)   :: exchangeFactor
Logical       :: OD_MOM
Real(realk)   :: MOM_CENTER(3)
!Derivative info
! First the number of components 
INTEGER       :: nGeoderivComp  !total nr. of geometrical derivative components
INTEGER       :: nGeoderivCompP !nr. of geometrical derivative components on P
INTEGER       :: nGeoderivCompQ !nr. of geometrical derivative components on Q
INTEGER       :: nMagderivCompP !nr. of magnetic derivative components on P
INTEGER       :: nMagderivCompQ !nr. of magnetic derivative components on Q
Integer       :: nMultipoleMomentComp !nr. of multipolmomentcomponents     
Integer       :: nCartesianMomentComp !nr. of cartesianmomentcomponents     
! Second the order of derivative 
INTEGER       :: GeoderOrderP ! order of geometrical derivative on P
INTEGER       :: GeoderOrderQ ! order of geometrical derivative on Q
INTEGER       :: MagderOrderP ! order of magnetic derivative on P
INTEGER       :: MagderOrderQ ! order of magnetic derivative on Q
Logical       :: DoMagScreen  ! Exploit Screening effect due to Qmn matrix in magnetic derivative integrals.
Integer       :: MMorder  !order of multipolmoments
Integer       :: CMorder !order of cartesianmoments
Integer       :: CMiMat  !specific cartesianmoment
INTEGER       :: geoderivOrder !the total order of geometrical derivative
INTEGER       :: MagderivOrder !the total order of magnetic derivative
Logical       :: DO_GRADIENT
Logical       :: AC_center

LOGICAL       :: AddToIntegral
!Attenuation parameters 
Real(realk)   :: ATTomega 
Real(realk)   :: ATTalpha
Real(realk)   :: ATTbeta
LOGICAL       :: ATTFACTOR
TYPE(GaussianGeminal) :: GGem
INTEGER       :: iPQxyz
integer       :: LinComCarmomType
Real(realk)   :: PROP_ORIGIN(3)
Logical       :: PropDerivEcoeff
Logical       :: PropKineticEcoeff
Logical       :: PropMomentEcoeff
Logical       :: PropAnti
integer       :: PropType
integer       :: PropMaxD
integer       :: PropMaxM
integer       :: PropRequireBoys
logical       :: fullcontraction
integer       :: operator
integer       :: node
integer       :: numnodes
logical       :: LHSsameAsRHSDmat
logical       :: RealGabMatrix
END TYPE INTEGRALINPUT

public :: INTEGRALINPUT,print_INTEGRALINPUT
private 

CONTAINS
SUBROUTINE print_INTEGRALINPUT(INT_INPUT,iunit)
implicit none
TYPE(INTEGRALINPUT),intent(IN) :: INT_INPUT
INTEGER,intent(IN)             :: IUNIT
!
integer :: i

WRITE(IUNIT,'(1X,A)') '**** Integral input ****'
WRITE(IUNIT,'(3X,A,F15.4)') 'Integralthreshold   = ',INT_INPUT%Integralthreshold
!TYPE(AOITEMPOINTER)   :: AO(4)
WRITE(IUNIT,'(3X,A,4I5)')   'AOdim               = ',(INT_INPUT%AOdim(i),i=1,4)
WRITE(IUNIT,'(3X,A,4I5)')   'currentFragment     = ',(INT_INPUT%currentFragment(i),i=1,4)
!type(lstensor),pointer :: LST_DLHS
!type(lstensor),pointer :: LST_DRHS
!type(lstensor),pointer :: LST_GAB_LHS
!type(lstensor),pointer :: LST_GAB_RHS
WRITE(IUNIT,'(3X,A,L2)')   'noOMP = ',INT_INPUT%noOMP
WRITE(IUNIT,'(3X,A,L2)')   'GAB_LHSusePointer    = ',INT_INPUT%GAB_LHSusePointer
WRITE(IUNIT,'(3X,A,L2)')   'GAB_RHSusePointer    = ',INT_INPUT%GAB_RHSusePointer
WRITE(IUNIT,'(3X,A,L2)')   'uselst_DRHS          = ',INT_INPUT%uselst_DRHS
WRITE(IUNIT,'(3X,A,L2)')   'uselst_DLHS          = ',INT_INPUT%uselst_DLHS
WRITE(IUNIT,'(3X,A,I5)')   'FTUVmaxprim          = ',INT_INPUT%FTUVmaxprim
WRITE(IUNIT,'(3X,A,I5)')   'LU_MMDATA            = ',INT_INPUT%LU_MMDATA
WRITE(IUNIT,'(3X,A,I5)')   'MMINDEX              = ',INT_INPUT%MMINDEX
WRITE(IUNIT,'(3X,A,I5)')   'LU_MMDATR            = ',INT_INPUT%LU_MMDATR
WRITE(IUNIT,'(3X,A,I5)')   'LU_MMDADA            = ',INT_INPUT%LU_MMDADA 
WRITE(IUNIT,'(3X,A,I5)')   'LU_MMDADR            = ',INT_INPUT%LU_MMDADR
WRITE(IUNIT,'(3X,A,I5)')   'NDMAT_LHS            = ',INT_INPUT%NDMAT_LHS
WRITE(IUNIT,'(3X,A,I5)')   'NDMAT_RHS            = ',INT_INPUT%NDMAT_RHS
WRITE(IUNIT,'(3X,A,L2)')   'DRHS_SYM             = ',INT_INPUT%DRHS_SYM
WRITE(IUNIT,'(3X,A,L2)')   'DLHS_SYM             = ',INT_INPUT%DLHS_SYM
WRITE(IUNIT,'(3X,A,L2)')   'sameLHSaos           = ',INT_INPUT%sameLHSaos
WRITE(IUNIT,'(3X,A,L2)')   'sameRHSaos           = ',INT_INPUT%sameRHSaos
WRITE(IUNIT,'(3X,A,L2)')   'sameODs              = ',INT_INPUT%sameODs
WRITE(IUNIT,'(3X,A,I5)')   'CENTERS              = ',INT_INPUT%CENTERS
WRITE(IUNIT,'(3X,A,I5)')   'ndim                 = ',INT_INPUT%ndim
WRITE(IUNIT,'(3X,A,L2)')   'DO_PASSES            = ',INT_INPUT%DO_PASSES
WRITE(IUNIT,'(3X,A,L2)')   'ContAng              = ',INT_INPUT%ContAng 
WRITE(IUNIT,'(3X,A,I5)')   'maxPasses            = ',INT_INPUT%maxPasses
WRITE(IUNIT,'(3X,A,L2)')   'HIGH_RJ000_ACCURACY  = ',INT_INPUT%HIGH_RJ000_ACCURACY
WRITE(IUNIT,'(3X,A,I5)')   'MMunique_ID1         = ',INT_INPUT%MMunique_ID1
WRITE(IUNIT,'(3X,A,I5)')   'MMunique_ID2         = ',INT_INPUT%MMunique_ID2
WRITE(IUNIT,'(3X,A,I5)')   'MMstartA             = ',INT_INPUT%MMstartA
WRITE(IUNIT,'(3X,A,I5)')   'MMstartB             = ',INT_INPUT%MMstartB
WRITE(IUNIT,'(3X,A,L2)')   'DO_COULOMB           = ',INT_INPUT%DO_COULOMB
WRITE(IUNIT,'(3X,A,L2)')   'DO_EXCHANGE          = ',INT_INPUT%DO_EXCHANGE
WRITE(IUNIT,'(3X,A,L2)')   'DO_FMM               = ',INT_INPUT%DO_FMM
WRITE(IUNIT,'(3X,A,L2)')   'DO_MULMOM            = ',INT_INPUT%DO_MULMOM
WRITE(IUNIT,'(3X,A,L2)')   'DO_PROP              = ',INT_INPUT%DO_PROP
WRITE(IUNIT,'(3X,A,L2)')   'DO_FOCK              = ',INT_INPUT%DO_FOCK
WRITE(IUNIT,'(3X,A,L2)')   'DO_ENERGY            = ',INT_INPUT%DO_ENERGY
WRITE(IUNIT,'(3X,A,L2)')   'DO_JENGINE           = ',INT_INPUT%DO_JENGINE
WRITE(IUNIT,'(3X,A,L2)')   'DO_DAJENGINE         = ',INT_INPUT%DO_DAJENGINE
WRITE(IUNIT,'(3X,A,L2)')   'DO_LINK              = ',INT_INPUT%DO_LINK
WRITE(IUNIT,'(3X,A,L2)')   'DO_DALINK            = ',INT_INPUT%DO_DALINK
WRITE(IUNIT,'(3X,A,I5)')   'DASCREEN_THRLOG      = ',INT_INPUT%DASCREEN_THRLOG
WRITE(IUNIT,'(3X,A,L2)')   'CS_int               = ',INT_INPUT%CS_int
WRITE(IUNIT,'(3X,A,L2)')   'PS_int               = ',INT_INPUT%PS_int
WRITE(IUNIT,'(3X,A,L2)')   'setETUVoutside       = ',INT_INPUT%setETUVoutside
WRITE(IUNIT,'(3X,A,L2)')   'sphericalEcoeff      = ',INT_INPUT%sphericalEcoeff
WRITE(IUNIT,'(3X,A,L2)')   'hermiteEcoeff        = ',INT_INPUT%hermiteEcoeff
WRITE(IUNIT,'(3X,A,L2)')   'CS_SCREEN            = ',INT_INPUT%CS_SCREEN 
WRITE(IUNIT,'(3X,A,F15.5)')'CS_THRESHOLD         = ',INT_INPUT%CS_THRESHOLD
WRITE(IUNIT,'(3X,A,I5)')   'CS_THRLOG            = ',INT_INPUT%CS_THRLOG
WRITE(IUNIT,'(3X,A,L2)')   'PS_SCREEN            = ',INT_INPUT%PS_SCREEN 
WRITE(IUNIT,'(3X,A,F15.5)')'PS_THRESHOLD         = ',INT_INPUT%PS_THRESHOLD
WRITE(IUNIT,'(3X,A,I5)')   'PS_THRLOG            = ',INT_INPUT%PS_THRLOG
WRITE(IUNIT,'(3X,A,L2)')   'PARI_SCREEN          = ',INT_INPUT%PARI_SCREEN 
WRITE(IUNIT,'(3X,A,F15.5)')'PARI_THRESHOLD       = ',INT_INPUT%PARI_THRESHOLD
WRITE(IUNIT,'(3X,A,L2)')   'OE_SCREEN            = ',INT_INPUT%OE_SCREEN 
WRITE(IUNIT,'(3X,A,F15.5)')'OE_THRESHOLD         = ',INT_INPUT%OE_THRESHOLD
WRITE(IUNIT,'(3X,A,L2)')   'OD_SCREEN            = ',INT_INPUT%OD_SCREEN 
WRITE(IUNIT,'(3X,A,F15.5)')'OD_THRESHOLD         = ',INT_INPUT%OD_THRESHOLD
WRITE(IUNIT,'(3X,A,L2)')   'MBIE_SCREEN          = ',INT_INPUT%MBIE_SCREEN 
WRITE(IUNIT,'(3X,A,L2)')   'MBIE_INT             = ',INT_INPUT%MBIE_INT
WRITE(IUNIT,'(3X,A,L2)')   'NonClassical_SCREEN  = ',INT_INPUT%NonClassical_SCREEN 
WRITE(IUNIT,'(3X,A,F15.5)')'MM_SCREENTHR         = ',INT_INPUT%MM_SCREENTHR
WRITE(IUNIT,'(3X,A,I5)')   'MM_TLMAX             = ',INT_INPUT%MM_TLMAX    
WRITE(IUNIT,'(3X,A,I5)')   'MM_LMAX              = ',INT_INPUT%MM_LMAX     
WRITE(IUNIT,'(3X,A,L2)')   'MM_NOONE             = ',INT_INPUT%MM_NOONE    
WRITE(IUNIT,'(3X,A,L2)')   'MM_NOSCREEN          = ',INT_INPUT%MM_NOSCREEN 
WRITE(IUNIT,'(3X,A,L2)')   'DO_MMGRD             = ',INT_INPUT%DO_MMGRD    
WRITE(IUNIT,'(3X,A,I5)')   'CS_MAXELM_LHS        = ',INT_INPUT%CS_MAXELM_LHS
WRITE(IUNIT,'(3X,A,I5)')   'CS_MAXELM_RHS        = ',INT_INPUT%CS_MAXELM_RHS
WRITE(IUNIT,'(3X,A,I5)')   'PS_MAXELM_LHS        = ',INT_INPUT%PS_MAXELM_LHS
WRITE(IUNIT,'(3X,A,I5)')   'PS_MAXELM_RHS        = ',INT_INPUT%PS_MAXELM_RHS
WRITE(IUNIT,'(3X,A,F15.5)')'CoulombFactor        = ',INT_INPUT%CoulombFactor
WRITE(IUNIT,'(3X,A,F15.5)')'exchangeFactor       = ',INT_INPUT%exchangeFactor
WRITE(IUNIT,'(3X,A,L2)')   'OD_MOM               = ',INT_INPUT%OD_MOM
WRITE(IUNIT,'(3X,A,F15.5)')'MOM_CENTER(3)        = ',INT_INPUT%MOM_CENTER(3)
WRITE(IUNIT,'(3X,A,I5)')   'nGeoderivCompP       = ',INT_INPUT%nGeoderivCompP
WRITE(IUNIT,'(3X,A,I5)')   'nGeoderivCompQ       = ',INT_INPUT%nGeoderivCompQ
WRITE(IUNIT,'(3X,A,I5)')   'nMagderivCompP       = ',INT_INPUT%nMagderivCompP
WRITE(IUNIT,'(3X,A,I5)')   'nMagderivCompQ       = ',INT_INPUT%nMagderivCompQ
WRITE(IUNIT,'(3X,A,I5)')   'nMultipoleMomentComp = ',INT_INPUT%nMultipoleMomentComp
WRITE(IUNIT,'(3X,A,I5)')   'nCartesianMomentComp = ',INT_INPUT%nCartesianMomentComp
WRITE(IUNIT,'(3X,A,I5)')   'GeoderOrderP         = ',INT_INPUT%GeoderOrderP
WRITE(IUNIT,'(3X,A,I5)')   'GeoderOrderQ         = ',INT_INPUT%GeoderOrderQ
WRITE(IUNIT,'(3X,A,I5)')   'MagderOrderP         = ',INT_INPUT%MagderOrderP
WRITE(IUNIT,'(3X,A,I5)')   'MagderOrderQ         = ',INT_INPUT%MagderOrderQ
WRITE(IUNIT,'(3X,A,L2)')   'DoMagScreen          = ',INT_INPUT%DoMagScreen 
WRITE(IUNIT,'(3X,A,I5)')   'MMorder              = ',INT_INPUT%MMorder
WRITE(IUNIT,'(3X,A,I5)')   'CMorder              = ',INT_INPUT%CMorder
WRITE(IUNIT,'(3X,A,I5)')   'geoderivOrder        = ',INT_INPUT%geoderivOrder
WRITE(IUNIT,'(3X,A,I5)')   'MagderivOrder        = ',INT_INPUT%MagderivOrder
WRITE(IUNIT,'(3X,A,L2)')   'DO_GRADIENT          = ',INT_INPUT%DO_GRADIENT
WRITE(IUNIT,'(3X,A,L2)')   'AC_center            = ',INT_INPUT%AC_center
WRITE(IUNIT,'(3X,A,L2)')   'AddToIntegral        = ',INT_INPUT%AddToIntegral
WRITE(IUNIT,'(3X,A,F15.5)')'ATTomega             = ',INT_INPUT%ATTomega 
WRITE(IUNIT,'(3X,A,F15.5)')'ATTalpha             = ',INT_INPUT%ATTalpha
WRITE(IUNIT,'(3X,A,F15.5)')'ATTbeta              = ',INT_INPUT%ATTbeta
WRITE(IUNIT,'(3X,A,L2)')   'ATTFACTOR            = ',INT_INPUT%ATTFACTOR
!TYPE(GaussianGeminal) ::   'GGem                 = ',INT_INPUT%GGem
WRITE(IUNIT,'(3X,A,I5)')   'iPQxyz               = ',INT_INPUT%iPQxyz
WRITE(IUNIT,'(3X,A,I5)')   'LinComCarmomType     = ',INT_INPUT%LinComCarmomType
WRITE(IUNIT,'(3X,A,F15.5)')'PROP_ORIGIN(3)       = ',INT_INPUT%PROP_ORIGIN(3)
WRITE(IUNIT,'(3X,A,L2)')   'PropDerivEcoeff      = ',INT_INPUT%PropDerivEcoeff
WRITE(IUNIT,'(3X,A,L2)')   'PropKineticEcoeff    = ',INT_INPUT%PropKineticEcoeff
WRITE(IUNIT,'(3X,A,L2)')   'PropMomentEcoeff     = ',INT_INPUT%PropMomentEcoeff
WRITE(IUNIT,'(3X,A,L2)')   'PropAnti             = ',INT_INPUT%PropAnti
WRITE(IUNIT,'(3X,A,I5)')   'PropType             = ',INT_INPUT%PropType
WRITE(IUNIT,'(3X,A,I5)')   'PropMaxD             = ',INT_INPUT%PropMaxD
WRITE(IUNIT,'(3X,A,I5)')   'PropMaxM             = ',INT_INPUT%PropMaxM
WRITE(IUNIT,'(3X,A,I5)')   'PropRequireBoys      = ',INT_INPUT%PropRequireBoys
WRITE(IUNIT,'(3X,A,L2)')   'fullcontraction      = ',INT_INPUT%fullcontraction
WRITE(IUNIT,'(3X,A,I5)')   'operator             = ',INT_INPUT%operator

END SUBROUTINE print_INTEGRALINPUT


END MODULE integral_type
