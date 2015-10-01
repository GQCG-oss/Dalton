!> @file
!> Contains configuration module. 

!> \brief Type definitions for configuration structures and routines for reading input.
!> \author S. Host and T. Kjaergaard
!> \date March 2010
module configurationType
use TYPEDEFTYPE
use profile_type, only: profileinput
use av_utilities, only: AVitem
use diagonalization, only: DiagItem
use typedeftype
use basis_typetype
use molecule_typetype, only: moleculeinfo
use opttype, only: OptItem
use response_wrapper_type_module, only: MCDinputitem, ALPHAinputitem, BETAinputitem, &
     & GAMMAinputitem, TPAinputitem, DTPAinputitem, ESGinputitem, ESDinputitem, &
     & RSPSOLVERinputitem,NMRinputitem
use lsdalton_response_type_mod, only: rsp_tasksitem
use optimization_type, only: opt_setting
use ls_dynamicsType, only: dyntype
use soeo_typedef, only: soeoItem_input
use decompMod, only: DecompItem
use lattice_type, only: lvec_list_t
use davidson_settings, only: RedSpaceItem
use arhDensity, only: solveritem
#ifdef HAS_PCMSOLVER
use ls_pcm_config, only: pcmtype
#endif
use precision

private
public :: responseitem,ConfigItem,LowAccuracyStartType,&
     & Set_Low_accuracy_start_settings, Revert_Low_accuracy_start_settings
!> \brief Contains info, settings and data for response part (defaults or read from input file).
!> \author T. Kjaergaard
!> \date April 2010
type responseitem
   !> Used to store info about MCD calculation
   type(MCDinputitem)  :: MCDinput
   !> Used to store info about polarizability (alpha) calculation
   type(ALPHAinputitem)  :: ALPHAinput
   !> Used to store info about 1st hyperpolarizability (beta) calculation
   type(BETAinputitem)  :: BETAinput
   !> Used to store info about 2nd hyperpolarizability (gamma) calculation
   type(GAMMAinputitem)  :: GAMMAinput
   !> Used to store info about standard TPA calculation.
   type(TPAinputitem)  :: TPAinput
   !> Used to store info about damped TPA calculation.
   type(DTPAinputitem)  :: DTPAinput
   !> Used to store info about excited state gradient calculation.
   type(ESGinputitem)  :: ESGinput
   !> Used to store info about excited state dipole moment calculation.
   type(ESDinputitem)  :: ESDinput
   !> Used to store info about solver that is used for calculation.
   type(RSPSOLVERinputitem) :: RSPSOLVERinput
   !> Used to store info about NMR 
   type(NMRinputitem)  :: NMRinput
   type(rsp_tasksitem) :: tasks
   logical :: noOpenRSP
end type responseitem

!> \brief Contains info, settings and data for entire calculation (defaults or read from input file).
!> \author S. Host
!> \date March 2010
type ConfigItem
   !> Logical unit number for LSDALTON.OUT
   integer              :: lupri
   !> Logical unit number for DALTON.ERR
   integer              :: luerr
   !> Should a Dec calculation be performed
   logical              :: doDEC
   !> Turns off DEC energy contribution for get_energy calls
   logical              :: noDecEnergy
   !> Should Memory Information be printed
   logical              :: PrintMemory
   !> Perform interaction energy calculation using Counter Poise Correction
   logical              :: InteractionEnergy
   !> Same SubSystems in Interaction energies
   logical              :: SameSubSystems
   !> Construct SubSystems Density matrix in Interaction energies
   logical              :: SubSystemDensity
   !> Used for Augmented Roothaan-Hall, direct density optimization etc.
   type(SolverItem),pointer     :: solver
   !> Used for davidson solver in SCF opt
   type(RedSpaceItem)   :: davidSCF
   !> Used for davidson solver in orbital localization
   type(RedSpaceItem)   :: davidOrbLoc
   !> Used to store OAO decomposition of overlap matrix + OAO transformed matrices
   type(DecompItem),pointer     :: decomp
   !> Used to store info about integrals
   type(integralconfig) :: integral
   !> Used to store info about density optimization type
   type(OptItem)        :: opt
   !> Used to store info about SCF averaging
   type(AvItem)         :: av
   !> Used to store info about diagonalization
   type(DiagItem)       :: diag
   !> Used to store info about molecule
   type(moleculeinfo)   :: molecule
   !> Used to store info about which atoms have which basisset
   type(basissetlibraryitem):: lib(nBasisBasParam)
   !> Used to store info about response calculation
   type(responseitem) :: response
   !> Used to store info about geometry optimization
   type(opt_setting)  :: optinfo
   !> Used to store info about dynamics
   type(dyntype) :: dynamics
   !> Used to store info about options in a soeo-calculation
   type(soeoItem_input) :: soeoinp
   !> For use in pbc
   type(lvec_list_t) :: latt_config
   !Only for testing new sparse matrix library, should be removed afterwards!
   logical            :: sparsetest
   !> Used to store info about profile options
   type(profileinput) :: prof
   !> Memory monitor for MPI calculations
   logical            :: mpi_mem_monitor
#ifdef HAS_PCMSOLVER
   !> Used to store info about Polarizable Continuum Model calculation
   type(pcmtype)      :: pcm
#endif
#ifdef MOD_UNRELEASED
   !> Used to store info about geometrical Hessian
   type(geoHessianconfig) :: geoHessian
#endif
   !> Max memory available on gpu measured in GB. By default set to 2 GB
   real(realk) :: GPUMAXMEM
   !> Should a excited state geometry optimization be performed 
   logical              :: doESGopt
   !> Skip LSDALTON calculation and calculate PLT file from existing density
   !> or orbital file (see type pltinfo)
   logical  :: DoPLT
   !> Information about PLT calculation (only used if doPLT=true)
   type(pltinfo) :: PLT
   !> Should we do an F12 calc which requires a CABS basis
   logical              :: doF12
   !> Should we do an RIMP2 calc which requires a AUX basis
   logical              :: doRIMP2
   !> do MPI testing of mpicopy_setting and mpicopy_screen
   logical              :: doTestMPIcopy
   !> do testing of the high-order derivative integrals (HODI)
   logical              :: doTestHodi
   !> maximum order of the high-order derivative integrals (HODI)
   integer              :: testHodiOrder
   !> skip SCF calculations
   logical              :: skipscfloop
   !> test papi
   logical              :: papitest
   !> Memory monitor for Matrices
   logical            :: mat_mem_monitor
end type ConfigItem

type LowAccuracyStartType
logical :: DaLink
INTEGER :: Dascreen_thrlog
REAL(REALK) :: THRESHOLD
logical :: HIGH_RJ000_ACCURACY
end type LowAccuracyStartType

contains
subroutine set_Low_accuracy_start_settings(lupri,ls,config,LAStype)
implicit none
!> Contains info, settings and data for entire calculation
type(configItem), intent(inout) :: config
!> Logical unit number for LSDALTON.OUT
integer, intent(in)             :: lupri
!> Object containing integral settings and molecule
type(lsitem), intent(inout)     :: ls
!> Object containing change
type(LowAccuracyStartType), intent(inout)     :: LAStype
real(realk) :: conv_factor
write(config%lupri,*) '======================================'
write(config%lupri,*) ' use Low accuracy start settings      '
write(config%lupri,*) '======================================'
write(*,*) '======================================'
write(*,*) ' use Low accuracy start settings      '
LAStype%dalink = ls%Setting%scheme%dalink
LAStype%DASCREEN_THRLOG = ls%Setting%scheme%Dascreen_thrlog
LAStype%THRESHOLD = ls%Setting%scheme%THRESHOLD
LAStype%HIGH_RJ000_ACCURACY = ls%Setting%scheme%HIGH_RJ000_ACCURACY

ls%Setting%scheme%dalink = .TRUE.
config%integral%DALINK = .TRUE.
ls%input%dalton%DALINK = .TRUE.
conv_factor = 1.0E-2_realk

if (config%decomp%cfg_unres) then
   config%opt%set_convergence_threshold = conv_factor*sqrt((config%decomp%nocca+config%decomp%noccb)*1.0E0_realk)
else
   config%opt%set_convergence_threshold = conv_factor*sqrt(config%decomp%nocc*2.0E0_realk)
endif
WRITE(lupri,*)
WRITE(config%LUPRI,"('Low accuracy start dynamic convergence threshold for gradient: ',E10.2)") &
     & config%opt%set_convergence_threshold

config%INTEGRAL%DASCREEN_THRLOG = -1
ls%input%dalton%DASCREEN_THRLOG = -1
ls%setting%scheme%DASCREEN_THRLOG = -1

ls%input%dalton%THRESHOLD = 1.0E-6_realk
config%integral%THRESHOLD = 1.0E-6_realk
ls%setting%scheme%THRESHOLD = 1.0E-6_realk

ls%input%dalton%HIGH_RJ000_ACCURACY = .FALSE.
config%integral%HIGH_RJ000_ACCURACY = .FALSE.
ls%setting%scheme%HIGH_RJ000_ACCURACY = .FALSE.

WRITE(config%LUPRI,'(A)')' '
WRITE(config%LUPRI,'(A60,ES12.4)')'The Overall Screening threshold is set to              :',config%integral%THRESHOLD
WRITE(config%LUPRI,'(A60,ES12.4)')'The Screening threshold used for Coulomb               :',&
     & config%integral%THRESHOLD*config%integral%J_THR
WRITE(config%LUPRI,'(A60,ES12.4)')'The Screening threshold used for Exchange              :',&
     &config%integral%THRESHOLD*config%integral%K_THR
WRITE(config%LUPRI,'(A60,ES12.4)')'The Screening threshold used for One-electron operators:',&
     &config%integral%THRESHOLD*config%integral%ONEEL_THR
if(config%integral%DALINK)THEN
   WRITE(config%LUPRI,'(A)')' '
   WRITE(config%LUPRI,'(A,ES12.4)')'   DaLink have been activated, so in addition to using ',&
        & config%integral%THRESHOLD*config%integral%K_THR      
   WRITE(config%LUPRI,'(A)')'   as a screening threshold on the integrals contribution to'
   WRITE(config%LUPRI,'(A)')'   the Fock matrix, we also use a screening threshold'
   WRITE(config%LUPRI,'(A,ES12.4)')'   on the integrals contribution to the Energy:       ',&
        &config%integral%THRESHOLD*config%integral%K_THR*(1.0E+1_realk**(-config%INTEGRAL%DASCREEN_THRLOG))
endif

end subroutine Set_Low_accuracy_start_settings

subroutine revert_Low_accuracy_start_settings(lupri,ls,config,LAStype)
implicit none
!> Contains info, settings and data for entire calculation
type(configItem), intent(inout) :: config
!> Logical unit number for LSDALTON.OUT
integer, intent(in)             :: lupri
!> Object containing integral settings and molecule
type(lsitem), intent(inout)     :: ls
!> Object containing change
type(LowAccuracyStartType), intent(inout)     :: LAStype
write(config%lupri,*) '======================================'
write(config%lupri,*) ' done with Low accuracy start         '
write(config%lupri,*) '======================================'
write(*,*) '======================================'

ls%Setting%scheme%dalink = LAStype%dalink
config%integral%DALINK = LAStype%dalink
ls%input%dalton%DALINK = LAStype%dalink

ls%Setting%scheme%Dascreen_thrlog = LAStype%DASCREEN_THRLOG
ls%input%dalton%DASCREEN_THRLOG = LAStype%DASCREEN_THRLOG
config%INTEGRAL%DASCREEN_THRLOG = LAStype%DASCREEN_THRLOG

ls%Setting%scheme%THRESHOLD = LAStype%THRESHOLD
ls%input%dalton%THRESHOLD = LAStype%THRESHOLD
config%integral%THRESHOLD = LAStype%THRESHOLD

ls%Setting%scheme%HIGH_RJ000_ACCURACY = LAStype%HIGH_RJ000_ACCURACY
ls%input%dalton%HIGH_RJ000_ACCURACY = LAStype%HIGH_RJ000_ACCURACY
config%integral%HIGH_RJ000_ACCURACY = LAStype%HIGH_RJ000_ACCURACY

config%opt%set_convergence_threshold = config%opt%cfg_convergence_threshold

WRITE(config%LUPRI,'(A)')' '
WRITE(config%LUPRI,'(A60,ES12.4)')'The Overall Screening threshold is set to              :',config%integral%THRESHOLD
WRITE(config%LUPRI,'(A60,ES12.4)')'The Screening threshold used for Coulomb               :',&
     & config%integral%THRESHOLD*config%integral%J_THR
WRITE(config%LUPRI,'(A60,ES12.4)')'The Screening threshold used for Exchange              :',&
     &config%integral%THRESHOLD*config%integral%K_THR
WRITE(config%LUPRI,'(A60,ES12.4)')'The Screening threshold used for One-electron operators:',&
     &config%integral%THRESHOLD*config%integral%ONEEL_THR
if(config%integral%DALINK)THEN
   WRITE(config%LUPRI,'(A)')' '
   WRITE(config%LUPRI,'(A,ES12.4)')'   DaLink have been activated, so in addition to using ',&
        & config%integral%THRESHOLD*config%integral%K_THR      
   WRITE(config%LUPRI,'(A)')'   as a screening threshold on the integrals contribution to'
   WRITE(config%LUPRI,'(A)')'   the Fock matrix, we also use a screening threshold'
   WRITE(config%LUPRI,'(A,ES12.4)')'   on the integrals contribution to the Energy:       ',&
        &config%integral%THRESHOLD*config%integral%K_THR*(1.0E+1_realk**(-config%INTEGRAL%DASCREEN_THRLOG))
endif

end subroutine Revert_Low_accuracy_start_settings


end module configurationType
