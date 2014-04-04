MODULE IntegralInterfaceModuleDF
  use precision
  use TYPEDEFTYPE
  use Typedef  
  use Matrix_module
  use Matrix_Operations
  use Integralparameters
  use LStiming
  use ls_Integral_Interface
  use IO
  use molecule_module
  use pari_mod
  use mat3d_mod
  use lstensor_operationsmod
  use lstensor_typetype
  use memory_handling
  use AtomSparse
  use linsolvdf
  use GCtransMod
  use screen_mod
  public :: II_get_df_coulomb_mat,II_get_df_J_gradient, &
       & II_get_df_exchange_mat, II_get_pari_df_exchange_mat
  private
contains
!> \brief Calculates the coulomb matrix using density fitting
!> \author S. Reine
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param F the coulomb matrix
SUBROUTINE II_get_df_coulomb_mat(LUPRI,LUERR,SETTING,D,F,ndmat)
IMPLICIT NONE
INTEGER               :: LUPRI,LUERR,ndmat
TYPE(MATRIX),target   :: D(ndmat),F(ndmat)
TYPE(LSSETTING)       :: SETTING
!
integer               :: usemat
TYPE(MATRIX),target   :: galpha,calpha
Integer             :: nbasis,naux,info
Real(realk)         :: TSTART,TEND
Character(80)       :: Filename
type(matrixp)       :: Jmat(1),Dmat(1)
type(matrixp)       :: Intmat(1)

IF (SETTING%SCHEME%PARI_J) THEN
   IF (SETTING%SCHEME%SIMPLE_PARI) THEN
      IF(ndmat.NE. 1) THEN
           WRITE(*,*)     "The PARI approximation isn't implemented for unrestricted cases yet."
           WRITE(LUPRI,*) "The PARI approximation isn't implemented for unrestricted cases yet."
           CALL LSQUIT('For the time being the pari code is &
           &not implemted for more than 1 density matrix - spam Simen Reine and Patrick Merlot',lupri)
      ENDIF
      CALL II_get_pari_df_coulomb_mat_simple(LUPRI,LUERR,SETTING,D(1),F(1))
   ELSE 
      CALL II_get_pari_df_coulomb_mat(LUPRI,LUERR,SETTING,D,F,ndmat)
      ! CALL II_get_NRpari_df_coulomb_mat(LUPRI,LUERR,SETTING,D,F)
   ENDIF
ELSE IF (SETTING%SCHEME%OVERLAP_DF_J) THEN
   IF (SETTING%SCHEME%DENSFIT.OR.SETTING%SCHEME%PARI_J) THEN
     IF(matrix_type .EQ. mtype_unres_dense)&
          &CALL LSQUIT('Density fitting not implemented for unrestricted - spam Simen Reine',lupri)
   ENDIF
   IF(ndmat.NE. 1) THEN
        WRITE(*,*)     "The PARI approximation isn't implemented for unrestricted cases yet."
        WRITE(LUPRI,*) "The PARI approximation isn't implemented for unrestricted cases yet."
        CALL LSQUIT('For the time being II_get_overlap_df_coulomb_mat is &
         &not implemted for more than 1 density matrix - spam Simen Reine and Patrik',lupri)
   ENDIF
  CALL II_get_overlap_df_coulomb_mat(LUPRI,LUERR,SETTING,D(1),F(1))
ELSE
  CALL II_get_regular_df_coulomb_mat(LUPRI,LUERR,SETTING,D,F,ndmat)
ENDIF

END SUBROUTINE II_get_df_coulomb_mat

!> \brief Calculates the density-fitted electron-electron repulsion contribution to the molecular gradient
!> \author S. Reine
!> \date 2010-03-01
!> \param eeGrad The electron-electron-repulsion gradient
!> \param DmatLHS The left-hand-side (or first electron) density matrix
!> \param DmatRHS The reft-hand-side (or second electron) density matrix
!> \param ndlhs The number of LHS density matrices
!> \param ndrhs The number of RHS density matrices
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_df_J_gradient(eeGrad,DmatLHS,DmatRHS,ndlhs,ndrhs,setting,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Real(realk),intent(INOUT)     :: eeGrad(3,setting%molecule(1)%p%nAtoms)
Integer,intent(IN)            :: lupri,luerr,ndlhs,ndrhs
Type(matrixp),intent(IN)      :: DmatLHS(ndlhs),DmatRHS(ndrhs)
!
!
integer                   :: nAtoms,iDmat,nlhs,nrhs
type(matrixp)             :: eeGradMat(1)
type(matrix),target       :: eeGradTarget
Type(matrixp),pointer     :: DLHS(:),DRHS(:)
logical                   :: same, save_screen
TYPE(matrix),target       :: calpha(ndrhs)
Character(80)             :: Filename
Real(realk)               :: eeGradtmp(3,setting%molecule(1)%p%nAtoms)

integer :: nbast,naux
logical :: ReCalcGab,saveNOSEGMENT

!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%J_THR

nbast = DmatLHS(1)%p%nrow
naux  = setting%molecule(1)%p%nbastAUX
same = ls_same_mats(DmatLHS,DmatRHS,ndlhs,ndrhs)
IF (.NOT.same) CALL lsQUIT('Error in II_get_df_J_gradient. Only working for same lhs and rhs dmat',-1)
!IF ((ndlhs.GT. 1).OR.(ndrhs.GT. 1)) CALL lsQUIT('Error in II_get_df_J_gradient. Only working for ndmat = 1!',-1)
DO idmat=1,ndlhs
  write(Filename,'(A8,I3)') 'LSCALPHA',idmat
  IF (.not.io_file_exist(Filename,SETTING%IO)) call lsquit('Error in II_get_df_J_gradient. CALPHA does not exsit!',-1)
  CALL mat_init(calpha(idmat),naux,1)
  call io_read_mat(calpha(idmat),Filename,Setting%IO,lupri,luerr)
ENDDO


nlhs = ndlhs
nrhs = ndrhs
call mem_alloc(DLHS,nlhs)
call mem_alloc(DRHS,nrhs)

eeGrad = 0E0_realk

!****************** Calculate (rho^e|tilde rho) ***********
DO idmat = 1,nlhs
   DLHS(idmat)%p => DmatLHS(idmat)%p
ENDDO
DO idmat = 1,nrhs
  DRHS(idmat)%p => calpha(idmat)
ENDDO

CALL ls_attachDmatToSetting(DLHS,nlhs,setting,'LHS',1,2,.TRUE.,lupri)
CALL ls_attachDmatToSetting(DRHS,nrhs,setting,'RHS',3,4,.TRUE.,lupri)

nAtoms = setting%molecule(1)%p%nAtoms

eeGradtmp = 0E0_realk
call initIntegralOutputDims(setting%Output,3,nAtoms,1,1,1)
IF(SETTING%SCHEME%FMM) THEN
   ! recalculate the moments
   ! turning off the screening (meaning the screening on the final values in the printmm routines !!)
   ! the reason for turning screening off it that the overlap index in the moments and derivative moments 
   ! has to match (problem if moment is written but not the derivative moment or the other way around)
   ! this should be fixed, the easiest (?) way could be to calculate both moments and derivative moments at the
   ! same time (as done in the FCK3 code)
   SAVE_SCREEN = SETTING%SCHEME%MM_NOSCREEN
   SETTING%SCHEME%MM_NOSCREEN = .TRUE.
   SETTING%SCHEME%CREATED_MMFILES=.false.
   SETTING%SCHEME%DO_MMGRD = .TRUE.
   !we turn off family type basis sets because it does not work
   !for FMM-GRADIENTS - when calculating both 1 and 2 electron 
   !contributions together
   saveNOSEGMENT = SETTING%SCHEME%NOSEGMENT
   SETTING%SCHEME%NOSEGMENT = .TRUE.
   !recalc primscreening matrix
   ReCalcGab = SETTING%SCHEME%ReCalcGab
   SETTING%SCHEME%recalcGab = .TRUE.
   CALL ls_multipolemoment(LUPRI,LUERR,SETTING,nbast,naux,  &
   & DLHS(1)%p%nrow,DLHS(1)%p%ncol,naux,1,&
   & AORdefault,AORdefault,AODFdefault,AOempty,RegularSpec,ContractedInttype,.TRUE.)
   ! now derivative moments
   CALL ls_multipolemoment(LUPRI,LUERR,SETTING,nbast,naux,  &
   & DLHS(1)%p%nrow,DLHS(1)%p%ncol,naux,1,&
   & AORdefault,AORdefault,AODFdefault,AOempty,GradientSpec,ContractedInttype,.TRUE.)
   SETTING%SCHEME%MM_NOSCREEN = SAVE_SCREEN
END IF
CALL ls_jengine(AORdefault,AORdefault,AODFdefault,AOempty,&
     &          CoulombOperator,GradientSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,eeGRADtmp,.FALSE.)
IF(SETTING%SCHEME%FMM)THEN
   call ls_jengineClassicalGRAD(eeGRADtmp,AORdefault,AORdefault,AODFdefault,AOempty,&
        & CoulombOperator,GradientSpec,ContractedInttype,SETTING,LUPRI,LUERR,nAtoms)
ENDIF
CALL ls_freeDmatFromSetting(setting)
eeGrad = eeGrad + eeGradtmp
IF(SETTING%SCHEME%FMM) SETTING%SCHEME%DO_MMGRD = .FALSE.
!**************** End calculate (rho^e|tilde rho) *********

!****************** Calculate (tilde rho^e|rho) ***********
DO idmat = 1,nlhs
   DLHS(idmat)%p => calpha(idmat)
ENDDO
DO idmat = 1,nrhs
  DRHS(idmat)%p => DmatRHS(idmat)%p
ENDDO

CALL ls_attachDmatToSetting(DLHS,nlhs,setting,'LHS',1,2,.TRUE.,lupri)
CALL ls_attachDmatToSetting(DRHS,nrhs,setting,'RHS',3,4,.TRUE.,lupri)

nAtoms = setting%molecule(1)%p%nAtoms

eeGradtmp = 0E0_realk
call initIntegralOutputDims(setting%Output,3,nAtoms,1,1,1)
CALL ls_jengine(AODFdefault,AOempty,AORdefault,AORdefault,&
     &          CoulombOperator,GradientSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,eeGRADtmp,.FALSE.)
IF(SETTING%SCHEME%FMM)THEN
   call ls_jengineClassicalGRAD(eeGRADtmp,AODFdefault,AOempty,AORdefault,AORdefault,&
     & CoulombOperator,GradientSpec,ContractedInttype,SETTING,LUPRI,LUERR,nAtoms)
ENDIF
CALL ls_freeDmatFromSetting(setting)
eeGrad = eeGrad + 0.5E0_realk*eeGradtmp
!**************** End calculate (tilde rho^e|rho) *********

!**************** Calculate (tilde rho^e|tilde rho) ***********
DO idmat = 1,nlhs
DLHS(idmat)%p => calpha(idmat)
ENDDO
DO idmat = 1,nrhs
  DRHS(idmat)%p => calpha(idmat)
ENDDO

CALL ls_attachDmatToSetting(DLHS,nlhs,setting,'LHS',1,2,.TRUE.,lupri)
CALL ls_attachDmatToSetting(DRHS,nrhs,setting,'RHS',3,4,.TRUE.,lupri)

nAtoms = setting%molecule(1)%p%nAtoms

eeGRADtmp = 0E0_realk
call initIntegralOutputDims(setting%Output,3,nAtoms,1,1,1)
CALL ls_jengine(AODFdefault,AOempty,AODFdefault,AOempty,&
     &          CoulombOperator,GradientSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,eeGRADtmp,.FALSE.)
IF(SETTING%SCHEME%FMM)THEN
   call ls_jengineClassicalGRAD(eeGRADtmp,AODFdefault,AOempty,AODFdefault,AOempty,&
     & CoulombOperator,GradientSpec,ContractedInttype,SETTING,LUPRI,LUERR,nAtoms)
ENDIF
eeGrad = eeGrad - 0.5E0_realk*eeGradtmp
CALL ls_freeDmatFromSetting(setting)
!************** End calculate (tilde rho^e|tilde rho) *********

call mem_dealloc(DLHS)
call mem_dealloc(DRHS)

DO idmat=1,ndrhs
  CALL mat_free(calpha(idmat))
ENDDO
IF(SETTING%SCHEME%FMM)THEN
   SETTING%SCHEME%NOSEGMENT = saveNOSEGMENT
   SETTING%SCHEME%ReCalcGab = ReCalcGab 
ENDIF

END SUBROUTINE II_get_df_J_gradient

!> \brief This subroutine calculates the NON-ROBUST Pari-Atomic density-fitted (PARI) Coulomb matrix
!> \author P. Merlot, S. Reine
!> \date 2010-02-03
!> \param lupri Default print-unit for output
!> \param luerr Default print-unit for termination
!> \param setting Contains information about the integral settings
!> \param D Density-matrix
!> \param F Coulomb contribution to the Fock- or KS-matrix
SUBROUTINE II_get_NRpari_df_coulomb_mat(LUPRI,LUERR,SETTING,D,F)
implicit none
TYPE(MATRIX),target,intent(in)    :: D
TYPE(MATRIX),target,intent(inout) :: F
TYPE(LSSETTING),intent(inout)     :: SETTING
INTEGER,intent(in)                :: LUPRI,LUERR
!
TYPE(MATRIX),target   :: cbeta,galpha,galphaFit
Integer               :: nAtoms,temp1,temp2
Integer               :: iAtom,jAtom,nOrbReg,nOrbAux,nBastAux,nBastReg
TYPE(MOLECULARORBITALINFO) :: orbitalInfo
INTEGER,pointer       :: numAtomicOrbitalsReg(:)
INTEGER,pointer       :: numAtomicOrbitalsAux(:)
INTEGER,pointer       :: startAtomicOrbitalsReg(:)
INTEGER,pointer       :: startAtomicOrbitalsAux(:)
type(matrixp)         :: intmat(1)
Integer               :: usemat
Integer               :: ndmat = 1
Integer               :: iAtomA,iAtomB,iAtomC,iAtomD,nRegA,nRegB,nRegC,nRegD,nAuxA,nAuxB,nAuxC,nAuxD
Integer               :: iAlpha,iBeta,iGamma,iRegA,iRegB,iRegC,iRegD
Integer               :: startRegA,startRegB,startRegC,startRegD,startAuxA,startAuxB,startAuxC,startAuxD
Integer               :: endRegA,endRegB,endRegC,endRegD,endAuxA,endAuxB,endAuxC,endAuxD
Integer               :: nAux
Integer               :: info
Real(realk),pointer   :: cbetafull(:,:,:)
Real(realk),pointer   :: extracted_alpha_ab(:,:,:)
Real(realk),pointer   :: extracted_alpha_beta(:,:)
Real(realk),pointer   :: galphafull(:,:,:),galphaFitfull(:,:,:)
Real(realk),pointer   :: Dfull(:,:,:)
Real(realk),pointer   :: Jfull(:,:,:,:,:)
TYPE(matrixp)         :: Jmat(1),Dmat(1)
TYPE(BLOCKINFO)       :: pairAtomic,auxAtomic
TYPE(FragmentInfo)    :: atomicFragments
Real(realk)           :: dfac,factor,TSTART,TEND,tefull,tsfull
Integer               :: iAO,iRegAfull,iRegBfull
TYPE(MAT3D),pointer   :: extracted_calpha_ab(:)
logical               :: saveRecalcGab
TYPE(LSTENSOR),pointer:: regCSfull,auxCSfull

saveRecalcGab = setting%scheme%recalcgab
setting%scheme%recalcgab = .TRUE.
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%J_THR
!Then get the full screening matrices
CALL II_getScreenMatFull(regCSfull,auxCSfull,Setting,lupri,luerr)
!
CALL LSTIMER('START ',TSTART,TEND,LUPRI)
CALL LSTIMER('START ',tsfull,tefull,lupri)
!
CALL getMolecularDimensions(SETTING%MOLECULE(1)%p,nAtoms,nBastReg,nBastAux)

IF(SETTING%SCHEME%FMM) CALL LSQUIT('FMM not yet implemented for PARI!',lupri)

! Alloc full density-matrix and Coulomb-matrix
call mem_alloc(Dfull,D%nrow,D%ncol,1)
call mat_to_full(D,1E0_realk,Dfull(:,:,1))
call mat_zero(F)              
Jmat(1)%p => F 
call mem_alloc(Jfull,nBastReg,nBastReg,1,1,ndmat)
call ls_dzero(Jfull,nBastReg*nBastReg*ndmat)

! --- Build an array with the list of atoms with their respective nb. of orbitals/auxiliary functions
! --- Build another array to know where each atom start in the complete matrices
CALL setMolecularOrbitalInfo(SETTING%MOLECULE(1)%p,orbitalInfo)

numAtomicOrbitalsReg   => orbitalInfo%numAtomicOrbitalsReg
startAtomicOrbitalsReg => orbitalInfo%startAtomicOrbitalsReg
numAtomicOrbitalsAux   => orbitalInfo%numAtomicOrbitalsAux
startAtomicOrbitalsAux => orbitalInfo%startAtomicOrbitalsAux
 
!WRITE(*,*) "nAtoms= ",nAtoms
!WRITE(*,*) "nBastReg= ",nBastReg
!WRITE(*,*) "nBastAux= ",nBastAux
!WRITE(*,*) "numAtomicOrbtialsReg= ",numAtomicOrbitalsReg
!WRITE(*,*) "startAtomicOrbitalsReg= ",startAtomicOrbitalsReg
!WRITE(*,*) "numAtomicOrbtialsAux= ",numAtomicOrbitalsAux
!WRITE(*,*) "startAtomicOrbitalsAux= ",startAtomicOrbitalsAux

! --- Calculate the pair-atomic fitting coefficients
allocate(extracted_calpha_ab(nAtoms))
CALL getPariCoefficients(LUPRI,LUERR,SETTING,extracted_calpha_ab,orbitalInfo,regCSfull,auxCSfull)

!re-set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%J_THR


CALL LSTIMER('coeffs',TSTART,TEND,LUPRI)


! --- compute Cbeta 
call mem_alloc(cbetafull,nBastAux,1,ndmat)
cbetafull = 0.0E0_realk
DO iAtomA=1,nAtoms
  call getAtomicOrbitalInfo(orbitalInfo,iAtomA,nRegA,startRegA,endRegA,nAuxA,startAuxA,endAuxA)
  DO iAtomB=1,nAtoms
    call getAtomicOrbitalInfo(orbitalInfo,iAtomB,nRegB,startRegB,endRegB,nAuxB,startAuxB,endAuxB)
    DO iRegB=1,nRegB
      iRegBfull = iRegB+startRegB-1
      DO iRegA=1,nRegA
        iRegAfull = iRegA+startRegA-1
        dfac = 2E0_realk*Dfull(iRegAfull,iRegBfull,ndmat)
        IF (iAtomA.NE.iAtomB) dfac = dfac + 2E0_realk * Dfull(iRegBfull,iRegAfull,ndmat)
        DO iAlpha=1,nAuxA
          cbetafull(iAlpha+startAuxA-1,1,ndmat) = cbetafull(iAlpha+startAuxA-1,1,ndmat) &
     &       + dfac * extracted_calpha_ab(iAtomA)%elements(iAlpha,iRegA,iRegBfull)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO

CALL LSTIMER('c-beta',TSTART,TEND,LUPRI)
! ----------------------------------------------------------------------------------
! --- COMPUTING THE 3 CONTRIBUTIONS TO THE COULOMB PART (Jab)

! --- add first contribution to Jab
call mat_init(cbeta,nBastAux,1)
call mat_set_from_full(cbetafull(:,:,1),1E0_realk,cbeta)
Dmat(1)%p => cbeta
CALL ls_attachDmatToSetting(Dmat,ndmat,setting,'RHS',3,4,.TRUE.,lupri)
!Jmat,ndmat,
call initIntegralOutputDims(setting%Output,F%nrow,F%ncol,1,1,1)
call ls_jengine(AORdefault,AORdefault,AODFdefault,AOempty,CoulombOperator,RegularSpec,ContractedInttype,&
     &          SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,F,.FALSE.)
IF(SETTING%SCHEME%FMM)THEN
   call ls_jengineClassicalMat(F,AORdefault,AORdefault,AODFdefault,AOempty,&
        & CoulombOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
ENDIF

CALL ls_freeDmatFromSetting(setting)

call mat_to_full(F,1E0_realk,Jfull(:,:,1,1,1))
CALL LSTIMER('J-1   ',TSTART,TEND,LUPRI)

! --- compute gAlpha
call mem_alloc(galphafull,nBastAux,1,ndmat)
call ls_DZERO(galphafull,nBastAux*ndmat)
call mat_init(galpha,nBastAux,1)
call mat_zero(galpha)
Jmat(1)%p => galpha  
Dmat(1)%p => D  
CALL ls_attachDmatToSetting(Dmat,ndmat,setting,'RHS',3,4,.TRUE.,lupri)
!Jmat,ndmat,
call initIntegralOutputDims(setting%Output,galpha%nrow,galpha%ncol,1,1,1)
call ls_jengine(AODFdefault,AOempty,AORdefault,AORdefault,CoulombOperator,RegularSpec,ContractedInttype,&
     &          SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,galpha,.FALSE.)
IF(SETTING%SCHEME%FMM)THEN
   call ls_jengineClassicalMat(galpha,AODFdefault,AOempty,AORdefault,AORdefault&
        & ,CoulombOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
ENDIF

CALL ls_freeDmatFromSetting(setting)
call mat_to_full(galpha,1E0_realk,galphafull(:,:,1))       
call mat_free(galpha)
CALL LSTIMER('galph1',TSTART,TEND,LUPRI)

DO iAtomA=1,nAtoms
  call getAtomicOrbitalInfo(orbitalInfo,iAtomA,nRegA,startRegA,endRegA,nAuxA,startAuxA,endAuxA)
  DO iAtomB=1,nAtoms
    call getAtomicOrbitalInfo(orbitalInfo,iAtomB,nRegB,startRegB,endRegB,nAuxB,startAuxB,endAuxB)
    DO iRegB=1,nRegB
      iRegBfull = iRegB+startRegB-1
      DO iRegA=1,nRegA
        iRegAfull = iRegA+startRegA-1
        DO iAlpha=1,nAuxA
          Jfull(iRegAfull,iRegBfull,1,1,1) &
     &      =  Jfull(iRegAfull,iRegBfull,1,1,1) &
     &        + extracted_calpha_ab(iAtomA)%elements(iAlpha,iRegA,iRegBfull)&
     &          * galphafull(iAlpha+startAuxA-1,1,1)
        ENDDO
      ENDDO
    ENDDO
    IF (iAtomA.NE.iAtomB) THEN
      DO iRegA=1,nRegA
        iRegAfull = iRegA+startRegA-1
        DO iRegB=1,nRegB
          iRegBfull = iRegB+startRegB-1
          DO iBeta=1,nAuxB
            Jfull(iRegAfull,iRegBfull,1,1,1) &
     &        =  Jfull(iRegAfull,iRegBfull,1,1,1) &
     &         + extracted_calpha_ab(iAtomB)%elements(iBeta,iRegB,iRegAfull)&
     &           * galphafull(iBeta+startAuxB-1,1,1)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDDO
ENDDO
CALL LSTIMER('Jfull',TSTART,TEND,LUPRI)

CALL mat_set_from_full(Jfull(:,:,1,1,1),0.5E0_realk,F)    

call mem_dealloc(galphafull)
! call mem_dealloc(galphaFitfull)
call mem_dealloc(cbetafull)
call mat_free(cbeta)

call mem_dealloc(Dfull)
CALL mem_dealloc(Jfull)
call freePariCoefficients(extracted_calpha_ab,orbitalInfo%nAtoms)
deallocate(extracted_calpha_ab)
CALL freeMolecularOrbitalInfo(orbitalInfo)
!CALL freeFragmentInfo(atomicFragments)
IF(setting%IntegralTransformGC)THEN
   call AO2GCAO_transform_matrixF(F,setting,lupri)
ENDIF
CALL LSTIMER('PARI-J',tsfull,tefull,lupri)
setting%scheme%recalcgab = saveRecalcGab

!WRITE(*,*) "End of II_get_pari_df_coulomb_mat "
!WRITE(*,*)
!WRITE(*,*)


CONTAINS
  SUBROUTINE checkConsistency(MOLECULE,nAtoms,nBastAux,nBastReg)
  implicit none
  TYPE(MOLECULE_PT),intent(in) :: MOLECULE(4)
  Integer,intent(in)           :: nAtoms,nBastAux,nBastReg
  !
  integer :: iAO
  ! Consistency checking
  DO iAO=2,4
    IF (nAtoms.NE.MOLECULE(iAO)%p%nAtoms) THEN
       CALL LSQUIT('Error in PARI! Different number of atoms for different AOs',-1)
    ENDIF
    IF (nBastAux.NE.MOLECULE(iAO)%p%nBastAux) THEN
       CALL LSQUIT('Error in PARI! Different number of basis funcitons for different AOs',-1)
    ENDIF
    IF ((nBastReg.NE.MOLECULE(iAO)%p%nBastReg).OR.(nBastReg.NE.MOLECULE(iAO)%p%nBastVAL)) THEN
       CALL LSQUIT('Error in PARI! Different number of auxiliary basis functions for different AOs',-1)
    ENDIF
  ENDDO
  END SUBROUTINE checkConsistency
END SUBROUTINE II_get_NRpari_df_coulomb_mat

SUBROUTINE II_getScreenMatFull(regCSfull,auxCSfull,Setting,lupri,luerr)
implicit none
TYPE(LSTENSOR),pointer        :: regCSfull,auxCSfull
TYPE(LSSETTING),intent(inout) :: SETTING
INTEGER,intent(in)            :: LUPRI,LUERR
!
type(moleculeinfo),pointer :: molecule
Character(80) :: Filename
Character(53) :: identifier
Logical       :: FoundInMem,dofit
Integer       :: molID,THR,THR2
Integer       :: i

IF (setting%scheme%CS_SCREEN) THEN
  THR = ABS(NINT(LOG10(SETTING%SCHEME%intTHRESHOLD)))
  molecule => SETTING%MOLECULE(1)%p
  molID = SETTING%molID(1)
  !Find the precalculated screening matrices (up to 10 orders of magnitude smaller 
  !than the current intTHRESHOLD
  DO I=0,10
    THR2 = THR + I
    CALL io_get_CSidentifier(identifier,THR2,molecule,molecule,Setting%Scheme%CS_SCREEN,Setting%Scheme%PS_SCREEN)
    CALL io_get_filename(Filename,identifier,AORdefault,AORdefault,AORdefault,AORdefault,molID,molID,&
        &              CoulombOperator,ContractedInttype,.FALSE.,LUPRI,LUERR)
    call determine_lst_in_screenlist(Filename,FoundInMem,SETTING%IO)
    IF (FoundInMem) THEN
      call screen_associate(regCSfull,Filename,FoundInMem)
      EXIT
    ENDIF
  ENDDO
  IF (.NOT.FoundInMem) THEN
    CALL io_get_CSidentifier(identifier,THR,molecule,molecule,Setting%Scheme%CS_SCREEN,Setting%Scheme%PS_SCREEN)
    CALL io_get_filename(Filename,identifier,AORdefault,AORdefault,AORdefault,AORdefault,molID,molID,&
      &                  CoulombOperator,ContractedInttype,.FALSE.,LUPRI,LUERR)
    write(lupri,'(1X,1A)') 'Error in II_getScreenMatFull no regular screening matrix found'
    write(lupri,'(1X,2A)') 'Identifier =',identifier
    write(lupri,'(1X,2A)') 'Filename ='  ,Filename
    CALL LSQUIT('Error in II_getScreenMatFull no regular screening matrix found',-1)
  ENDIF
  
  dofit = setting%scheme%densfit.OR.setting%scheme%pari_j .OR. setting%scheme%pari_k
  IF (dofit) THEN
    DO I=0,10
      THR2 = THR + I
      CALL io_get_CSidentifier(identifier,THR2,molecule,molecule,Setting%Scheme%CS_SCREEN,Setting%Scheme%PS_SCREEN)
      CALL io_get_filename(Filename,identifier,AODFdefault,AOempty,AODFdefault,AOempty,molID,molID,&
            &              CoulombOperator,ContractedInttype,.FALSE.,LUPRI,LUERR)
      call determine_lst_in_screenlist(Filename,FoundInMem,SETTING%IO)
      IF (FoundInMem) THEN
        call screen_associate(auxCSfull,Filename,FoundInMem)
        EXIT
      ENDIF
    ENDDO
    IF (.NOT.FoundInMem) THEN
      CALL io_get_CSidentifier(identifier,THR,molecule,molecule,Setting%Scheme%CS_SCREEN,Setting%Scheme%PS_SCREEN)
      CALL io_get_filename(Filename,identifier,AODFdefault,AOempty,AODFdefault,AOempty,molID,molID,&
            &              CoulombOperator,ContractedInttype,.FALSE.,LUPRI,LUERR)
      write(lupri,'(1X,1A)') 'Error in II_getScreenMatFull no auxiliary screening matrix found'
      write(lupri,'(1X,2A)') 'Identifier =',identifier
      write(lupri,'(1X,2A)') 'Filename ='  ,Filename
      CALL LSQUIT('Error in II_getScreenMatFull no auxiliary screening matrix found',-1)
    ENDIF
  ELSE
    NULLIFY(auxCSfull)
  ENDIF
ELSE
  NULLIFY(regCSfull)
  NULLIFY(auxCSfull)
ENDIF

END SUBROUTINE II_getScreenMatFull

!> \brief This subroutine calculates the pari-atomic density-fitted (PARI) Coulomb matrix
!> \author P. Merlot, S. Reine
!> \date 2010-02-03
!> \param lupri Default print-unit for output
!> \param luerr Default print-unit for termination
!> \param setting Contains information about the integral settings
!> \param D Density-matrix
!> \param F Coulomb contribution to the Fock- or KS-matrix
SUBROUTINE II_get_pari_df_coulomb_mat(LUPRI,LUERR,SETTING,D,F,ndmat)
implicit none
TYPE(MATRIX),target,intent(in)    :: D(ndmat)
TYPE(MATRIX),target,intent(inout) :: F(ndmat)
TYPE(LSSETTING),intent(inout)     :: SETTING
INTEGER,intent(in)                :: LUPRI,LUERR,ndmat
!
TYPE(MATRIX),target        :: galpha,galphaFit
Integer                    :: nAtoms,temp1,temp2
Integer                    :: iAtom,jAtom,nOrbReg,nOrbAux,nBastAux,nBastReg
TYPE(MOLECULARORBITALINFO) :: orbitalInfo
INTEGER,pointer            :: numAtomicOrbitalsReg(:)
INTEGER,pointer            :: numAtomicOrbitalsAux(:)
INTEGER,pointer            :: startAtomicOrbitalsReg(:)
INTEGER,pointer            :: startAtomicOrbitalsAux(:)
type(matrixp)              :: intmat(1)
Integer                    :: usemat
Integer                    :: iAtomA,iAtomB,iAtomC,iAtomD,nRegA,nRegB,nRegC,nRegD,nAuxA,nAuxB,nAuxC,nAuxD
Integer                    :: iAlpha,iBeta,iGamma,iRegA,iRegB,iRegC,iRegD
Integer                    :: startRegA,startRegB,startRegC,startRegD,startAuxA,startAuxB,startAuxC,startAuxD
Integer                    :: endRegA,endRegB,endRegC,endRegD,endAuxA,endAuxB,endAuxC,endAuxD
Integer                    :: nAux
Integer                    :: info
Real(realk),pointer        :: cbetafull(:,:,:)
Real(realk),pointer        :: extracted_alpha_ab(:,:,:)
Real(realk),pointer        :: extracted_alpha_beta(:,:)
Real(realk),pointer        :: galphafull(:,:,:),galphaFitfull(:,:,:)
Real(realk),pointer        :: Dfull(:,:,:),DfullAO(:,:,:)
Real(realk),pointer        :: Jfull(:,:,:)
TYPE(matrixp)              :: Jmat(1),Dmat(1)
TYPE(BLOCKINFO)            :: pairAtomic,auxAtomic
TYPE(FragmentInfo)         :: atomicFragments
Real(realk)                :: dfac,factor,TSTART,TEND,tefull,tsfull
Integer                    :: iAO,iRegAfull,iRegBfull
TYPE(MAT3D),pointer        :: extracted_calpha_ab(:)
Integer                    :: ireg,iaux
CHARACTER(LEN=3)           :: nline
logical                    :: saveRecalcGab
integer                    :: idmat,nmat,nrow,ncol
TYPE(LSTENSOR),pointer     :: regCSfull,auxCSfull

IF (ndmat.GT.1) THEN
   WRITE(*,*)     "The PARI approximation isn't implemented for unrestricted cases yet."
   WRITE(LUPRI,*) "The PARI approximation isn't implemented for unrestricted cases yet."
   CALL LSQUIT('Error in II_get_pari_df_coulomb_mat: ndmat>1 not tested',-1)
ENDIF
nrow = D(1)%nrow
ncol = D(1)%ncol
IF ((F(1)%nrow.NE.nrow).OR.(F(1)%ncol.NE.ncol)) CALL LSQUIT('Error in II_get_pari_df_coulomb_mat F/D',-1)

IF (matrix_type .EQ. mtype_unres_dense)THEN
  IF (SETTING%SCHEME%NON_ROBUST_PARI) THEN 
     WRITE(*,*)     "The PARI approximation isn't implemented for unrestricted cases yet."
     WRITE(LUPRI,*) "The PARI approximation isn't implemented for unrestricted cases yet."
     CALL LSQUIT('Error in II_get_pari_df_coulomb_mat. NR and unrestricted',-1)
  ENDIF
  IF(setting%IntegralTransformGC) THEN
    WRITE(*,*)     "The PARI approximation isn't implemented for unrestricted cases yet."
    WRITE(LUPRI,*) "The PARI approximation isn't implemented for unrestricted cases yet."
    CALL LSQUIT('Error in II_get_pari_df_coulomb_mat. GC and unrestricted',-1)
  ENDIF
  IF (SETTING%SCHEME%FMM) THEN
    WRITE(*,*)     "The PARI approximation isn't implemented for unrestricted cases yet."
    WRITE(LUPRI,*) "The PARI approximation isn't implemented for unrestricted cases yet."
    call lsquit('Not allowed combination in II_get_pari_df_coulomb_mat. FMM and unrestricted',-1)
  ENDIF
  nmat = 2*ndmat
  call mem_alloc(Dfull,nrow,ncol,nmat)
  DO Idmat=1,ndmat
    CALL DCOPY(nrow*ncol,D(idmat)%elms,1,Dfull(:,:,2*idmat-1),1)
    CALL DCOPY(nrow*ncol,D(idmat)%elmsb,1,Dfull(:,:,2*idmat),1)
  ENDDO
ELSE
  nmat = ndmat
  ! Alloc full density-matrix and Coulomb-matrix
  call mem_alloc(Dfull,nrow,ncol,ndmat)
  DO idmat=1,ndmat
    call mat_to_full(D(idmat),1E0_realk,Dfull(:,:,idmat))
  ENDDO
ENDIF

IF(setting%IntegralTransformGC)THEN
   call mem_alloc(DfullAO,nrow,ncol,nmat)
   call GCAO2AO_transform_fullD(Dfull,DfullAO,nrow,nmat,setting,lupri)
   call mem_dealloc(Dfull)
ELSE
   DfullAO => Dfull
ENDIF

call mem_alloc(Jfull,nrow,ncol,nmat)
call ls_dzero(Jfull,nrow*ncol*nmat)

CALL LSTIMER('START ',TSTART,TEND,LUPRI)
CALL LSTIMER('START ',tsfull,tefull,lupri)
!
CALL getMolecularDimensions(SETTING%MOLECULE(1)%p,nAtoms,nBastReg,nBastAux)

IF(SETTING%SCHEME%FMM) CALL LSQUIT('FMM not yet implemented for PARI!',lupri)

! --- Build an array with the list of atoms with their respective nb. of orbitals/auxiliary functions
! --- Build another array to know where each atom start in the complete matrices
CALL setMolecularOrbitalInfo(SETTING%MOLECULE(1)%p,orbitalInfo)

numAtomicOrbitalsReg   => orbitalInfo%numAtomicOrbitalsReg
startAtomicOrbitalsReg => orbitalInfo%startAtomicOrbitalsReg
numAtomicOrbitalsAux   => orbitalInfo%numAtomicOrbitalsAux
startAtomicOrbitalsAux => orbitalInfo%startAtomicOrbitalsAux
 
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%J_THR
!Then get the full screening matrices
CALL II_getScreenMatFull(regCSfull,auxCSfull,Setting,lupri,luerr)

! --- Calculate the pair-atomic fitting coefficients
allocate(extracted_calpha_ab(nAtoms))
CALL getPariCoefficients(LUPRI,LUERR,SETTING,extracted_calpha_ab,orbitalInfo,regCSfull,auxCSfull)

!reset threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%J_THR

CALL LSTIMER('coeffs',TSTART,TEND,LUPRI)
  
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%J_THR

DO idmat=1,nmat
  
  ! --- compute Cbeta 
  call mem_alloc(cbetafull,nBastAux,1,1)
  cbetafull = 0.0E0_realk
  DO iAtomA=1,nAtoms
    call getAtomicOrbitalInfo(orbitalInfo,iAtomA,nRegA,startRegA,endRegA,nAuxA,startAuxA,endAuxA)
    DO iAtomB=1,nAtoms
      call getAtomicOrbitalInfo(orbitalInfo,iAtomB,nRegB,startRegB,endRegB,nAuxB,startAuxB,endAuxB)
      DO iRegB=1,nRegB
        iRegBfull = iRegB+startRegB-1
        DO iRegA=1,nRegA
          iRegAfull = iRegA+startRegA-1
          dfac = 2E0_realk*DfullAO(iRegAfull,iRegBfull,idmat)
          IF (iAtomA.NE.iAtomB) dfac = dfac + 2E0_realk * DfullAO(iRegBfull,iRegAfull,idmat)
          DO iAlpha=1,nAuxA
            cbetafull(iAlpha+startAuxA-1,1,1) = cbetafull(iAlpha+startAuxA-1,1,1) &
       &       + dfac * extracted_calpha_ab(iAtomA)%elements(iAlpha,iRegA,iRegBfull)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  
  
  CALL LSTIMER('c-beta',TSTART,TEND,LUPRI)
  ! ----------------------------------------------------------------------------------
  ! --- COMPUTING THE 3 CONTRIBUTIONS TO THE COULOMB PART (Jab)
  
  CALL ls_attachDmatToSetting(cbetafull,nBastAux,1,1,setting,'RHS',3,4,lupri)
  CALL ls_attach_gab_to_setting(setting,regCSfull,auxCSfull)
  call initIntegralOutputDims(setting%Output,nBastReg,nBastReg,1,1,1)
  call ls_jengine(AORdefault,AORdefault,AODFdefault,AOempty,CoulombOperator,RegularSpec,ContractedInttype,&
       &          SETTING,LUPRI,LUERR)
  CALL retrieve_Output(lupri,setting,Jfull(:,:,idmat),.FALSE.)
  IF(SETTING%SCHEME%FMM)THEN
     IF (ndmat.GT.1) CALL LSQUIT('not testet',-1)
     CALL LSQUIT('not testet',-1)
     !call ls_jengineClassicalFull
  ENDIF

  CALL ls_freeDmatFromSetting(setting)
  CALL ls_free_gab_from_setting(setting,lupri)
  
  CALL LSTIMER('J-1   ',TSTART,TEND,LUPRI)
  
  ! --- compute gAlpha
  call mem_alloc(galphafull,nBastAux,1,1)
  call ls_DZERO(galphafull,nBastAux)
  CALL ls_attachDmatToSetting(DfullAO(:,:,idmat:idmat),nBastReg,nBastReg,1,setting,'RHS',3,4,lupri)
  CALL ls_attach_gab_to_setting(setting,auxCSfull,regCSfull)
  call initIntegralOutputDims(setting%Output,nBastAux,1,1,1,1)
  call ls_jengine(AODFdefault,AOempty,AORdefault,AORdefault,CoulombOperator,RegularSpec,ContractedInttype,&
       &          SETTING,LUPRI,LUERR)
  CALL retrieve_Output(lupri,setting,galphafull,.FALSE.)
  IF(SETTING%SCHEME%FMM)THEN
     IF (ndmat.GT.1) CALL LSQUIT('not testet',-1)
     CALL LSQUIT('not testet',-1)
     !call ls_jengineClassicalFull
  ENDIF

  CALL ls_freeDmatFromSetting(setting)
  CALL ls_free_gab_from_setting(setting,lupri)
  CALL LSTIMER('galph3',TSTART,TEND,LUPRI)
  
  ! --- compute gAlphaFit
  call mem_alloc(galphaFitfull,nBastAux,1,1)
  call ls_DZERO(galphaFitfull,nBastAux)
  CALL ls_attachDmatToSetting(cbetafull,nBastAux,1,1,setting,'RHS',3,4,lupri)
  CALL ls_attach_gab_to_setting(setting,auxCSfull,auxCSfull)
  !Jmat,ndmat,
  call initIntegralOutputDims(setting%Output,nBastAux,1,1,1,1)
  call ls_jengine(AODFdefault,AOempty,AODFdefault,AOempty,CoulombOperator,RegularSpec,ContractedInttype,&
       &          SETTING,LUPRI,LUERR)
  CALL retrieve_Output(lupri,setting,galphaFitFull,.FALSE.)
  IF(SETTING%SCHEME%FMM)THEN
     IF (ndmat.GT.1) CALL LSQUIT('not testet',-1)
     CALL LSQUIT('not testet',-1)
     !call ls_jengineClassicalFull
  ENDIF
  CALL ls_freeDmatFromSetting(setting)
  CALL ls_free_gab_from_setting(setting,lupri)
  
  CALL LSTIMER('galph2',TSTART,TEND,LUPRI)
  
  ! --- contract gAlpha and gAlphaFit, then add the second/third contributions to Jab
  ! --- TO BE OPTIMIZED using the sparcity of C_alpha_ab_full
  IF (.NOT.(SETTING%SCHEME%NON_ROBUST_PARI)) THEN
     DO iAlpha=1,nBastAux
        galphafull(iAlpha,1,1)=galphafull(iAlpha,1,1)-galphaFitfull(iAlpha,1,1)
     ENDDO
  ENDIF
  !
  DO iAtomA=1,nAtoms
    call getAtomicOrbitalInfo(orbitalInfo,iAtomA,nRegA,startRegA,endRegA,nAuxA,startAuxA,endAuxA)
    DO iAtomB=1,nAtoms
      call getAtomicOrbitalInfo(orbitalInfo,iAtomB,nRegB,startRegB,endRegB,nAuxB,startAuxB,endAuxB)
      DO iRegB=1,nRegB
        iRegBfull = iRegB+startRegB-1
        DO iRegA=1,nRegA
          iRegAfull = iRegA+startRegA-1
          DO iAlpha=1,nAuxA
            Jfull(iRegAfull,iRegBfull,idmat) &
       &      =  Jfull(iRegAfull,iRegBfull,idmat) &
       &        + extracted_calpha_ab(iAtomA)%elements(iAlpha,iRegA,iRegBfull)&
       &          * galphafull(iAlpha+startAuxA-1,1,1)
          ENDDO
        ENDDO
      ENDDO
      IF (iAtomA.NE.iAtomB) THEN
        DO iRegA=1,nRegA
          iRegAfull = iRegA+startRegA-1
          DO iRegB=1,nRegB
            iRegBfull = iRegB+startRegB-1
            DO iBeta=1,nAuxB
              Jfull(iRegAfull,iRegBfull,idmat) &
       &        =  Jfull(iRegAfull,iRegBfull,idmat) &
       &         + extracted_calpha_ab(iAtomB)%elements(iBeta,iRegB,iRegAfull)&
       &           * galphafull(iBeta+startAuxB-1,1,1)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  ENDDO

  call mem_dealloc(galphafull)
  call mem_dealloc(galphaFitfull)
  call mem_dealloc(cbetafull)

ENDDO
  
CALL LSTIMER('Jfull',TSTART,TEND,LUPRI)


IF (matrix_type .EQ. mtype_unres_dense)THEN
  DO Idmat=1,ndmat
    CALL DCOPY(nrow*ncol,Jfull(:,:,2*idmat-1),1,F(idmat)%elms, 1)
    CALL DCOPY(nrow*ncol,Jfull(:,:,2*idmat  ),1,F(idmat)%elmsb,1)
  ENDDO
ELSE
  factor = 1E0_realk
  IF (SETTING%SCHEME%NON_ROBUST_PARI) factor = 0.5_realk
  DO Idmat=1,ndmat
    call mat_set_from_full(Jfull(:,:,idmat),factor,F(idmat))
    IF(setting%IntegralTransformGC)THEN
       call AO2GCAO_transform_matrixF(F(idmat),setting,lupri)
    ENDIF
  ENDDO
ENDIF

call mem_dealloc(DfullAO)
CALL mem_dealloc(Jfull)
call freePariCoefficients(extracted_calpha_ab,orbitalInfo%nAtoms)
deallocate(extracted_calpha_ab)
CALL freeMolecularOrbitalInfo(orbitalInfo)

CALL LSTIMER('PARI-J',tsfull,tefull,lupri)

CONTAINS
  SUBROUTINE checkConsistency(MOLECULE,nAtoms,nBastAux,nBastReg)
  implicit none
  TYPE(MOLECULE_PT),intent(in) :: MOLECULE(4)
  Integer,intent(in)           :: nAtoms,nBastAux,nBastReg
  !
  integer :: iAO
  ! Consistency checking
  DO iAO=2,4
    IF (nAtoms.NE.MOLECULE(iAO)%p%nAtoms) THEN
       CALL LSQUIT('Error in PARI! Different number of atoms for different AOs',-1)
    ENDIF
    IF (nBastAux.NE.MOLECULE(iAO)%p%nBastAux) THEN
       CALL LSQUIT('Error in PARI! Different number of basis funcitons for different AOs',-1)
    ENDIF
    IF ((nBastReg.NE.MOLECULE(iAO)%p%nBastReg).OR.(nBastReg.NE.MOLECULE(iAO)%p%nBastVAL)) THEN
       CALL LSQUIT('Error in PARI! Different number of auxiliary basis functions for different AOs',-1)
    ENDIF
  ENDDO
  END SUBROUTINE checkConsistency
END SUBROUTINE II_get_pari_df_coulomb_mat



!> \brief This subroutine calculates the pair-atomic density-fitted (PARI) Coulomb matrix
!> \author P. Merlot, S. Reine
!> \date 2011-04-19
!> \param lupri Default print-unit for output
!> \param luerr Default print-unit for termination
!> \param setting Contains information about the integral settings
!> \param D Density-matrix
!> \param F Coulomb contribution to the Fock- or KS-matrix
SUBROUTINE II_get_pari_df_coulomb_mat_simple(LUPRI,LUERR,SETTING,D,F)
  implicit none
  TYPE(MATRIX),target,intent(in)    :: D
  TYPE(MATRIX),target,intent(inout) :: F
  TYPE(LSSETTING),intent(inout)     :: SETTING
  INTEGER,intent(in)                :: LUPRI,LUERR
  !                                    
  Integer               :: nAtoms,nBastAux,nBastReg
  Real(realk),pointer   :: ab_beta_full(:,:,:)
  Real(realk),pointer   :: ab_beta_AB(:,:,:)

  Real(realk),pointer   :: alpha_beta_full(:,:,:,:,:)   ! in retrieve_output_2dim, not implemented with just 2 dimensions !!!
  Real(realk),pointer   :: alpha_beta_AB(:,:) 
  Real(realk),pointer   :: Calpha_ab_full(:,:,:)
!  Real(realk),pointer   :: Calpha_ab_AB(:,:,:)
  Real(realk),pointer   :: Calpha_full(:)
  Real(realk),pointer   :: galpha_full(:)
  Real(realk),pointer   :: galphaFit_full(:)
  Real(realk),pointer   :: Dfull(:,:,:),DfullAO(:,:,:)
  Real(realk),pointer   :: Jfull(:,:,:,:,:)
  Integer               :: ndmat = 1
  TYPE(MOLECULARORBITALINFO) :: orbitalInfo
  Integer               :: nRegA,nRegB,nAuxA,nAuxB,nAuxAB
  Integer               :: iAtomA,iAtomB,iAuxA,iAuxB,iReg,iRegA,iRegB,iRegC,iRegD
  Integer               :: iAux,iAlpha,iBeta
  Integer               :: startRegA,startRegB,startAuxA,startAuxB
  Integer               :: endRegA,  endRegB,  endAuxA,  endAuxB
  Integer               :: info
  Real(realk)           :: tmp
  CHARACTER(LEN=3) :: nline
  logical               :: saveRecalcGab
!   TYPE(MATRIX)          :: F_copy  ! copy of the Fock matrix
!   TYPE(MATRIX)          :: S   !overlap matrix
!   TYPE(MATRIX)          :: Cmo ! eigenvector of the Fock matrix
!   Real(realk), pointer  :: eival(:)  ! eigenvalues of the fock matrix
  !
  !set threshold 
  SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%J_THR
  saveRecalcGab = setting%scheme%recalcgab
  setting%scheme%recalcgab = .TRUE.
  !
  ! Read Molecule infos
  call getMolecularDimensions(SETTING%MOLECULE(1)%p,nAtoms,nBastReg,nBastAux)
  ! Alloc full density-matrix Dfull and Coulomb-matrix Jfull
  call mem_alloc(Dfull,D%nrow,D%ncol,1)
  call mat_to_full(D,1E0_realk,Dfull(:,:,1))

  IF(setting%IntegralTransformGC)THEN
     call mem_alloc(DfullAO,D%nrow,D%ncol,1)
     call GCAO2AO_transform_fullD(Dfull,DfullAO,D%nrow,1,setting,lupri)
     call mem_dealloc(Dfull)
  ELSE
     DfullAO => Dfull
  ENDIF

  call mat_zero(F)
  call mem_alloc(Jfull,nBastReg,nBastReg,1,1,ndmat)
  call ls_dzero(Jfull,nBastReg*nBastReg*ndmat)

  ! Calculate (ab|beta) full
  call mem_alloc(ab_beta_full,nBastAux,nBastReg,nBastReg)
  ab_beta_full = 0E0_realk ! call ls_dzero(ab_beta_full,nBastAux*nBastReg*nBastReg)
  call initIntegralOutputDims(setting%Output,nBastAux,1,nBastReg,nBastReg,1)
  call ls_getIntegrals(AODFdefault,AOempty,AORdefault,AORdefault,CoulombOperator,RegularSpec,&
       &               ContractedInttype,SETTING,LUPRI,LUERR)
  call retrieve_Output(lupri,setting,ab_beta_full,.FALSE.)

  ! Calculate (alpha|beta) full
  call mem_alloc(alpha_beta_full,nBastAux,1,nBastAux,1,1)
  alpha_beta_full = 0E0_realk ! call ls_dzero(alpha_beta_full,nBastAux*nBastReg*nBastReg)
  call initIntegralOutputDims(setting%Output,nBastAux,1,nBastAux,1,1)
  call ls_getIntegrals(AODFdefault,AOempty,AODFdefault,AOempty,CoulombOperator,RegularSpec,&
       &               ContractedInttype,SETTING,LUPRI,LUERR)
  call retrieve_Output(lupri,setting,alpha_beta_full,.FALSE.)

  ! Alloc memory for Calpha_ab_full
  call mem_alloc(Calpha_ab_full,nBastAux,nBastReg,nBastReg)

#if 0
  DO iRegB=1,nBastReg
    DO iRegA=1,nBastReg
      DO iAux=1,nBastAux
        Calpha_ab_full(iAux,iRegA,iRegB) = ab_beta_full(iAux,iRegA,iRegB) 
      ENDDO
    ENDDO
  ENDDO
! Calpha_ab_full = ab_beta_full
  call mem_alloc(alpha_beta_AB,nBastAux,nBastAux)
  DO iAuxB=1,nBastAux
    DO iAuxA=1,nBastAux
      alpha_beta_AB(iAuxA,iAuxB) = alpha_beta_full(iAuxA,1,iAuxB,1,1)
    ENDDO
  ENDDO
! alpha_beta_AB = alpha_beta_full(:,1,:,1,1)
  ! --- solve the linear system: (alpha|beta) Calpha_ab_AB = (ab|beta) e.g. A X = B
  call DPOSV('U',nBastAux,nBastReg*nBastReg,alpha_beta_AB,nBastAux,Calpha_ab_full,nBastAux,info)
  If (info.NE. 0) THEN
     WRITE(LUPRI,'(1X,A,I5)') 'DPOSV error in II_get_pari_df_coulomb_mat_simple. Info =',info
     call LSQUIT('DPOSV error in II_get_pari_df_coulomb_mat_simple',lupri)
  ENDIF
  call mem_dealloc(alpha_beta_AB)
#else
  Calpha_ab_full = 0E0_realk
  ! Calculate Calpha_ab_AB for each {A,B} atom pairs and insert them into Calpha_ab_full
  call setMolecularOrbitalInfo(SETTING%MOLECULE(1)%p,orbitalInfo)
  DO iAtomA=1,nAtoms
     call getAtomicOrbitalInfo(orbitalInfo,iAtomA,nRegA,startRegA,endRegA,nAuxA,startAuxA,endAuxA)
     DO iAtomB=1,nAtoms
        call getAtomicOrbitalInfo(orbitalInfo,iAtomB,nRegB,startRegB,endRegB,nAuxB,startAuxB,endAuxB)
        nAuxAB = nAuxA+nAuxB
        IF (iAtomA.EQ.iAtomB) nAuxAB = nAuxA
        
        ! --- extract (alpha|beta) for alpha,beta in {AUB}
        call mem_alloc(alpha_beta_AB,nAuxAB,nAuxAB)
        alpha_beta_AB = 0E0_realk
!         alpha_beta_AB(1:nAuxA,1:nAuxA) &
!              & = alpha_beta_full(startAuxA:endAuxA,1,startAuxA:endAuxA,1,1)
!         IF (iAtomA.NE.iAtomB) THEN
!           alpha_beta_AB(1:nAuxA,nAuxA+1:nAuxAB) &
!                & = alpha_beta_full(startAuxA:endAuxA,1,startAuxB:endAuxB,1,1)
!           alpha_beta_AB(nAuxA+1:nAuxAB,1:nAuxA) &
!                & = alpha_beta_full(startAuxB:endAuxB,1,startAuxA:endAuxA,1,1)
!           alpha_beta_AB(nAuxA+1:nAuxAB,nAuxA+1:nAuxAB) &
!                & = alpha_beta_full(startAuxB:endAuxB,1,startAuxB:endAuxB,1,1)
!         ENDIF

        DO iAux=1,nAuxA
           DO iAuxA=1,nAuxA
              alpha_beta_AB(iAuxA,iAux) &
                   & = alpha_beta_full(iAuxA+startAuxA-1,1,iAux+startAuxA-1,1,1)
           ENDDO
        ENDDO

        IF (iAtomA.NE.iAtomB) THEN
           DO iAuxB=1,nAuxB
              DO iAuxA=1,nAuxA
                 alpha_beta_AB(iAuxA,nAuxA+iAuxB) &
                      & = alpha_beta_full(iAuxA+startAuxA-1,1,iAuxB+startAuxB-1,1,1)
              ENDDO
           ENDDO
           DO iAuxA=1,nAuxA
              DO iAuxB=1,nAuxB
                 alpha_beta_AB(nAuxA+iAuxB,iAuxA) &
                      & = alpha_beta_full(iAuxB+startAuxB-1,1,iAuxA+startAuxA-1,1,1)
              ENDDO
           ENDDO
           DO iAux=1,nAuxB
              DO iAuxB=1,nAuxB
                 alpha_beta_AB(nAuxA+iAuxB,nAuxA+iAux) &
                      & = alpha_beta_full(iAuxB+startAuxB-1,1,iAux+startAuxB-1,1,1)
              ENDDO
           ENDDO
        ENDIF

        ! --- extract (ab|beta) for a in A, b in B and beta in {AUB}
        call mem_alloc(ab_beta_AB,nAuxAB,nRegA,nRegB)
        ab_beta_AB = 0E0_realk
!         ab_beta_AB(1:nAuxA,1:nRegA,1:nRegB) &
!              & = ab_beta_full(startAuxA:endAuxA,startRegA:endRegA,startRegB:endRegB)
!         IF (iAtomA.NE.iAtomB) THEN
!            ab_beta_AB(nAuxA+1:nAuxAB,1:nRegA,1:nRegB) &
!                 & = ab_beta_full(startAuxB:endAuxB,startRegA:endRegA,startRegB:endRegB)
!         ENDIF

        DO iRegB=1,nRegB
           DO iRegA=1,nRegA
              DO iAuxA=1,nAuxA
                 ab_beta_AB(iAuxA,iRegA,iRegB) &
                      & = ab_beta_full(startAuxA-1+iAuxA,startRegA-1+iRegA,startRegB-1+iRegB)
              ENDDO
              IF (iAtomA.NE.iAtomB) THEN
                 DO iAuxB=1,nAuxB
                    ab_beta_AB(nauxA+iAuxB,iRegA,iRegB) &
                         & = ab_beta_full(startAuxB-1+iAuxB,startRegA-1+iRegA,startRegB-1+iRegB)
                 ENDDO
              ENDIF
           ENDDO
        ENDDO


        ! --- alloc memory for Calpha_ab_AB
!        call mem_alloc(Calpha_ab_AB,nAuxAB,nRegA,nRegB)
!        Calpha_ab_AB = 0E0_realk
        
        ! --- check if (alpha|beta) is a NxN symmetric matrix with N=nAuxAB
        DO iAuxA=1,nAuxAB
           DO iAuxB=iauxA,nAuxAB
              If (ABS(alpha_beta_AB(iAuxA,iAuxB) - alpha_beta_AB(iAuxB,iAuxA)) .GT. 1.0E-10) THEN
                 WRITE(LUPRI,'(1X,A,I5)') 'Non symmetric matrix in II_get_pari_df_coulomb_mat_simple. Info =',info
                 call LSQUIT('Non symmetric matrix in II_get_pari_df_coulomb_mat_simple',lupri)
              ENDIF
           ENDDO
        ENDDO





        !write(*,*) '(alpha|beta) on AB'
        iaux = 0
        DO iauxa=1,nauxab
           DO iauxb=1,nauxab
              iaux=iaux+1
              nline='no'
              if(iaux .eq. nauxab) THEN
                 iaux=0
                 nline='yes'
              ENDIF
              !write(*,'(E12.4)',advance=nline)  alpha_beta_AB(iauxa,iauxb)
           ENDDO
        ENDDO


        !write(*,'(A,I2,A,I2,A,I2,A,I2,A,I2)') 'iA=',iAtomA,', iB=',iAtomB,', nRegAB=',nrega*nregB,', nRegA/nRegB=',nRegA,'/',nRegB


        ireg = 0
        DO iauxa=1,nauxab
           !write(*,'(A,I2)') '(ab|beta) on AB with beta=',iauxa
           DO irega=1,nrega
              DO iregb=1,nregb
                 ireg=ireg+1
                 nline='no'
                 if(ireg .eq. nregb) THEN
                    ireg=0
                    nline='yes'
                 ENDIF
                 !write(*,'(E12.4)',advance=nline)  ab_beta_AB(iauxa,irega,iregb)
              ENDDO
           ENDDO
        ENDDO
        
        !write(*,'(A,I2,A,I2)') 'DPOSV: nAuxAB=',nAuxAB,' nRegAB=',nrega*nregB

        ! --- solve the linear system: (alpha|beta) Calpha_ab_AB = (ab|beta) e.g. A X = B
        call DPOSV('U',nAuxAB,nRegA*nRegB,alpha_beta_AB,nAuxAB,ab_beta_AB,nAuxAB,info)
        If (info.NE. 0) THEN
           WRITE(LUPRI,'(1X,A,I5)') 'DPOSV error in II_get_pari_df_coulomb_mat_simple. Info =',info
           call LSQUIT('DPOSV error in II_get_pari_df_coulomb_mat_simple',lupri)
        ENDIF


        ireg = 0                                                                              
        DO iauxa=1,nauxab
           !write(*,'(A20,I3)') 'Calpha_ab_AB iaux=',iauxa
           DO irega=1,nrega                                                                
              DO iregb=1,nregb                                                                
                 ireg=ireg+1                                                                     
                 nline='no'                                                                      
                 if(ireg .eq. nregb) THEN                                                     
                    ireg=0                                                                       
                    nline='yes'                                                                  
                 ENDIF
                 !write(*,'(E12.4)',advance=nline)  ab_beta_AB(iauxa,irega,iregb)
              ENDDO
           ENDDO
        ENDDO

        ! --- map Calpha_ab_AB (stored into ab_beta_AB) into Calpha_ab_full
!         Calpha_ab_full(startAuxA:endAuxA,startRegA:endRegA,startRegB:endRegB) &
!              & = ab_beta_AB(1:nAuxA,1:nRegA,1:nRegB)
!         IF (iAtomA.NE.iAtomB) THEN
!            Calpha_ab_full(startAuxB:endAuxB,startRegA:endRegA,startRegB:endRegB) &
!                 & = ab_beta_AB(nAuxA+1:nAuxAB,1:nRegA,1:nRegB)
!         ENDIF

        Do iRegB=1,nRegB
           Do iRegA=1,nRegA
              DO iAuxA=1,nAuxA
                 Calpha_ab_full(startAuxA-1+iAuxA,startRegA-1+iRegA,startRegB-1+iRegB) &
                      & = ab_beta_AB(iAuxA,iRegA,iRegB)
              ENDDO
              IF (iAtomA.NE.iAtomB) THEN
                 DO iAuxB=1,nAuxB
                    Calpha_ab_full(startAuxB-1+iAuxB,startRegA-1+iRegA,startRegB-1+iRegB) &
                         & = ab_beta_AB(nAuxA+iAuxB,iRegA,iRegB)
              ENDDO
              ENDIF
           ENDDO
        ENDDO


        ireg = 0
        DO iauxa=1,nBastAux
           !write(*,'(A20,I3)') 'Calpha_ab_full iaux=',iauxa
           DO irega=1,nbastreg
              DO iregb=1,nbastreg
                 ireg=ireg+1
                 nline='no'
                 if(ireg .eq. nregb) THEN
                    ireg=0
                    nline='yes'
                 ENDIF
                 !write(*,'(E12.4)',advance=nline)  Calpha_ab_full(iauxa,irega,iregb)
              ENDDO
           ENDDO
        ENDDO


        ! --- dealloc memory
        call mem_dealloc(alpha_beta_AB) 
        call mem_dealloc(ab_beta_AB)    
     ENDDO ! iAtomB
  ENDDO ! iAtomA
  call freeMolecularOrbitalInfo(orbitalInfo)

#endif

  
  ! Calculate Calpha_full
  call mem_alloc(Calpha_full,nBastAux)
  Calpha_full = 0E0_realk
  DO iAux=1,nBastAux
     tmp = 0E0_realk
     DO iRegD=1,nBastReg
        DO iRegC=1,nBastReg
           tmp = tmp + Calpha_ab_full(iAux,iRegC,iRegD) * 2E0_realk * DfullAO(iRegC,iRegD,ndmat)
!          Calpha_full(iAux) = Calpha_full(iAux) + Calpha_ab_full(iAux,iRegC,iRegD) * DfullAO(iRegC,iRegD,ndmat)
        ENDDO
     ENDDO
     Calpha_full(iAux) = tmp
  ENDDO


  !write(*,*) 'CbetaFull:'
  DO iauxa=1,nbastaux
     !write(*,'(E12.4)')  Calpha_full(iauxa)
  ENDDO
  ! Calculate first contribution to Jab
  DO iRegB=1,nBastReg
     DO iRegA=1,nBastReg
        tmp = 0E0_realk
        DO iAux=1,nBastAux
           tmp = tmp + Calpha_full(iAux) * ab_beta_full(iAux,iRegA,iRegB)
        ENDDO
        Jfull(iRegA,iRegB,1,1,ndmat) = tmp
     ENDDO
  ENDDO


  !write(*,*) 'Jfull before correction terms'
  ireg = 0
  DO irega=1,nbastreg
     DO iregb=1,nbastreg
        ireg=ireg+1
        nline='no'
        if(ireg .eq. nbastreg) THEN
           ireg=0
           nline='yes'
        ENDIF
        !write(*,'(E12.4)',advance=nline)  Jfull(irega,iregb,1,1,1)
     ENDDO
  ENDDO

  ! Calculate galpha (part of 2nd contribution to Jab)
  call mem_alloc(galpha_full,nBastAux)
  DO iAux=1,nBastAux
     tmp = 0E0_realk
     DO iRegD=1,nBastReg
        DO iRegC=1,nBastReg
!            tmp = tmp + ab_beta_full(iAux,iRegC,iRegD) * 2E0_realk * Dfull(iRegC,iRegD,ndmat)
           tmp = tmp + ab_beta_full(iAux,iRegC,iRegD) * 2E0_realk * DfullAO(iRegC,iRegD,ndmat)
        ENDDO
     ENDDO
     galpha_full(iAux) = tmp
  ENDDO

  ! Calculate galphaFit (part of 3rd contribution to Jab)
  call mem_alloc(galphaFit_full,nBastAux)
  DO iAlpha=1,nBastAux
     tmp = 0E0_realk
     DO iBeta=1,nBastAux
        tmp = tmp + alpha_beta_full(iAlpha,1,iBeta,1,1) * Calpha_full(iBeta)
     ENDDO
     galphaFit_full(iAlpha) = tmp
  ENDDO

  ! Add up all contribution of Jab
  IF (SETTING%SCHEME%NON_ROBUST_PARI) THEN
     DO iRegB=1,nBastReg
        DO iRegA=1,nBastReg
           tmp = 0E0_realk
           DO iAlpha=1,nBastAux
              tmp = tmp + Calpha_ab_full(iAlpha,iRegA,iRegB) * galpha_full(iAlpha) ! NON-ROBUST version, still variational
           ENDDO
           Jfull(iRegA,iRegB,1,1,ndmat) = Jfull(iRegA,iRegB,1,1,ndmat) + tmp
        ENDDO
     ENDDO
  ELSE
     DO iRegB=1,nBastReg
        DO iRegA=1,nBastReg
           tmp = 0E0_realk
           DO iAlpha=1,nBastAux
              tmp = tmp + Calpha_ab_full(iAlpha,iRegA,iRegB) * (galpha_full(iAlpha) - galphaFit_full(iAlpha)) ! ROBUST version
           ENDDO
           Jfull(iRegA,iRegB,1,1,ndmat) = Jfull(iRegA,iRegB,1,1,ndmat) + tmp
        ENDDO
     ENDDO
  ENDIF

  call mem_dealloc(galphaFit_full)
  call mem_dealloc(galpha_full)
  
  ! store results into F matrix
  IF (SETTING%SCHEME%NON_ROBUST_PARI) THEN
     call mat_set_from_full(Jfull(:,:,1,1,1),0.5_realk,F) !NON-ROBUST PARI
  ELSE
     call mat_set_from_full(Jfull(:,:,1,1,1),1E0_realk,F) ! ROBUST PARI
  ENDIF

!   ! get eigenvalues of the Fock matrix
!   call mat_init(Cmo,nBastReg,nBastReg)
!   call mat_zero(Cmo)
!   call mem_alloc(eival,nBastReg)
!   call mat_init(F_copy,F%nrow,F%ncol)
!   call mat_zero(F_copy)
!   call mat_copy(1.0_realk,F_copy,F)
!   ! Overlap matrix                                                                                       
!   call mat_init(S,nBastReg,nBastReg)
!   call mat_zero(S)
!   call II_get_overlap(LUPRI,LUERR,SETTING,S)
!   call mat_diag_f (F_copy, S, eival, Cmo)   
!   ! check for negative eigenvalues
!   DO iRegA=1,nBastReg
!      IF (eival(iRegA) <= 0) THEN
!         write (*,*) '     negative eigenvalue(',iRegA,' : ',eival(iRegA)
!      ENDIF
!   ENDDO
!    call mat_free(Cmo)
!    call mat_free(F_copy)
!    call mat_free(S)
!    call mem_dealloc(eival)


   
  IF(setting%IntegralTransformGC)THEN
     call AO2GCAO_transform_matrixF(F,setting,lupri)
  ENDIF
  ! Free Memory
  call mem_dealloc(ab_beta_full)
!  call mem_dealloc(ab_beta_AB)

  call mem_dealloc(alpha_beta_full)
!  call mem_dealloc(alpha_beta_AB)

  call mem_dealloc(Calpha_ab_full)
!  call mem_dealloc(Calpha_ab_AB)
  call mem_dealloc(Calpha_full)


!  call mem_dealloc(galpha_full)
!  call mem_dealloc(galphaFit_full)

!  call mem_dealloc(Dfull)
  call mem_dealloc(DfullAO)
  call mem_dealloc(Jfull)
  !call freeMolecularOrbitalInfo(orbitalInfo)
  setting%scheme%recalcgab = saveRecalcGab
END SUBROUTINE II_get_pari_df_coulomb_mat_simple

!> \brief This subroutine calculates the pair-atomic density-fitted (PARI) Exchange matrix
!> \author P. Merlot, S. Reine
!> \date 2010-03-29
!> \param lupri Default print-unit for output
!> \param luerr Default print-unit for termination
!> \param setting Contains information about the integral settings
!> \param D Density-matrix
!> \param F Exchange contribution to the Fock- or KS-matrix
SUBROUTINE II_get_pari_df_exchange_mat(LUPRI,LUERR,SETTING,D,F,ndmat)
IMPLICIT NONE
TYPE(MATRIX),target,intent(in)    :: D(ndmat)
TYPE(MATRIX),target,intent(inout) :: F(ndmat)
TYPE(LSSETTING),intent(inout)     :: SETTING
Integer,intent(in)                :: LUPRI,LUERR,ndmat
!
Real(realk)                :: TSTART,TEND
Integer                    :: nAtoms,nBastAux,nBastReg
Real(realk),pointer        :: Dfull(:,:,:),DfullAO(:,:,:)
Real(realk),pointer        :: Kfull_3cContrib(:,:,:,:,:)
Real(realk),pointer        :: Kfull_2cContrib(:,:,:,:,:)
Real(realk),pointer        :: Kfull(:,:,:,:,:)
TYPE(matrixp)              :: Kmat(1),Dmat(1)
TYPE(MOLECULARORBITALINFO) :: orbitalInfo
Real(realk),pointer        :: alpha_beta(:,:,:,:,:)
type(matrixp)              :: intmat(1)
Integer                    :: usemat
Integer                    :: iAtomA,iAtomB,iAtomC,iAtomD
TYPE(MAT3D),pointer        :: calpha_ab(:)
Real(realk),pointer        :: dalpha_ad(:,:,:)

TYPE(MATRIX)               :: matrixK
Integer                    :: nRegA,nRegB,nRegC,nRegD,nAuxA,nAuxB,nAuxC,nAuxD
Integer                    :: startRegA,startRegB,startRegC,startRegD,startAuxA,startAuxB,startAuxC,startAuxD
Integer                    :: endRegA,endRegB,endRegC,endRegD,endAuxA,endAuxB,endAuxC,endAuxD
Integer                    :: iAlpha,iRegA,iRegB,iRegC,iRegD
!
TYPE(AtomSparseMat)        :: alphaBeta
TYPE(MoleculeInfo),pointer      :: molecule
Real(realk) :: ts,te,tsfull,tefull
logical               :: saveRecalcGab
integer :: idmat,nmat,nrow,ncol
TYPE(LSTENSOR),pointer :: regCSfull,auxCSfull

IF (ndmat.GT.1) THEN
   WRITE(*,*)     "The PARI approximation isn't implemented for unrestricted cases yet."
   WRITE(LUPRI,*) "The PARI approximation isn't implemented for unrestricted cases yet."
   CALL LSQUIT('Error in II_get_pari_df_exchange_mat: ndmat>1 not tested',-1)
ENDIF
nrow = D(1)%nrow
ncol = D(1)%ncol
IF ((F(1)%nrow.NE.nrow).OR.(F(1)%ncol.NE.ncol)) CALL LSQUIT('Error in II_get_pari_df_exchange_mat F/D',-1)

IF (matrix_type .EQ. mtype_unres_dense)THEN
  IF (SETTING%SCHEME%NON_ROBUST_PARI) THEN 
   WRITE(*,*)     "The PARI approximation isn't implemented for unrestricted cases yet."
   WRITE(LUPRI,*) "The PARI approximation isn't implemented for unrestricted cases yet."
   CALL LSQUIT('Error in II_get_pari_df_exchange_mat. NR and unrestricted',-1)
  ENDIF
  IF(setting%IntegralTransformGC) THEN
     WRITE(*,*)     "The PARI approximation isn't implemented for unrestricted cases yet."
     WRITE(LUPRI,*) "The PARI approximation isn't implemented for unrestricted cases yet."
     CALL LSQUIT('Error in II_get_pari_df_exchange_mat. GC and unrestricted',-1)
  ENDIF
  IF (SETTING%SCHEME%FMM) THEN
     WRITE(*,*)     "The PARI approximation isn't implemented for unrestricted cases yet."
     WRITE(LUPRI,*) "The PARI approximation isn't implemented for unrestricted cases yet."
     call lsquit('Not allowed combination in II_get_pari_df_exchange_mat. FMM and unrestricted',-1)
  ENDIF
  nmat = 2*ndmat
  call mem_alloc(Dfull,nrow,ncol,nmat)
  DO Idmat=1,ndmat
    CALL DCOPY(nrow*ncol,D(idmat)%elms,1,Dfull(:,:,2*idmat-1),1)
    CALL DCOPY(nrow*ncol,D(idmat)%elmsb,1,Dfull(:,:,2*idmat),1)
  ENDDO
ELSE
  nmat = ndmat
  ! Alloc full density-matrix and Coulomb-matrix
  call mem_alloc(Dfull,nrow,ncol,ndmat)
  DO idmat=1,ndmat
    call mat_to_full(D(idmat),1E0_realk,Dfull(:,:,idmat))
  ENDDO
ENDIF

!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%K_THR
!Then get the full screening matrices
CALL II_getScreenMatFull(regCSfull,auxCSfull,Setting,lupri,luerr)
!
molecule => SETTING%MOLECULE(1)%p
!
call getMolecularDimensions(molecule,nAtoms,nBastReg,nBastAux)

!call init_AtomSparseMat(alphaBeta,molecule,molecule,AODFdefault,AODFdefault,ContractedInttype,lupri)
CALL LSTIMER('START ',te,ts,lupri)
CALL LSTIMER('START ',tefull,tsfull,lupri)

IF(SETTING%SCHEME%FMM) call LSQUIT('FMM not yet implemented for PARI-K!',lupri)

call mem_alloc(Kfull_3cContrib,nBastReg,nBastReg,1,1,nmat)
call mem_alloc(Kfull_2cContrib,nBastReg,nBastReg,1,1,nmat)
call ls_DZERO(Kfull_3cContrib,nBastReg*nBastReg*nmat)
call ls_DZERO(Kfull_2cContrib,nBastReg*nBastReg*nmat)

! --- Build an array with the list of atoms with their respective nb. of orbitals/auxiliary functions
! --- Build another array to know where each atom start in the complete matrices
call setMolecularOrbitalInfo(molecule,orbitalInfo)

! --- Get screening integrals G_ab and G_alpha
!> \todo (Patrick) Get screening integrals G_ab and G_alpha
 
! --- Set up domains of atoms
!> \todo (Patrick) Set up domains of atoms

! --- Build the full huge matrices (alpha | beta)
usemat = 0 
call mem_alloc(alpha_beta,nBastAux,1,nBastAux,1,1)
call ls_DZERO(alpha_beta,nBastAux*nBastAux)
call ls_attach_gab_to_setting(setting,auxCSfull,auxCSfull)
call initIntegralOutputDims(setting%output,nBastaux,1,nBastaux,1,1)
call ls_getIntegrals(AODFdefault,AOempty,&
     &AODFdefault,AOempty,CoulombOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
call retrieve_Output(lupri,setting,alpha_beta,.FALSE.)
call ls_free_gab_from_setting(setting,lupri)

! --- Memory allocation (all Density coefficients stored in an INefficient way)

CALL LSTIMER('(al|be)',te,ts,lupri)

! --- Calculate the pair-atomic fitting coefficients 
allocate(calpha_ab(nAtoms))
call getPariCoefficients(LUPRI,LUERR,SETTING,calpha_ab,orbitalInfo,regCSfull,auxCSfull)

CALL LSTIMER('Coeffs ',te,ts,lupri)

!re-set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%K_THR

DO idmat=1,nmat
  ! 
  ! Insert call to screening_matrices over the full molecule for both AORdefault,AORdefault and AODFdefault
  ! interactions
  !
  ! --- Loop over atoms A
  DO iAtomA=1,nAtoms
     call getAtomicOrbitalInfo(orbitalInfo,iAtomA,nRegA,startRegA,endRegA,nAuxA,startAuxA,endAuxA)
  
    ! --- Construct density coefficients (dalpha_ad), alpha,a in A and d in Mol.
    call mem_alloc(dalpha_ad,nAuxA,nRegA,nBastReg)
  
    !> \todo (Patrick) limit the B atoms of the calpha_ab coefficients to the vicinity of A
    call getDensityCoefficients(iAtomA,orbitalInfo,calpha_ab,Dfull(:,:,idmat:idmat),dalpha_ad)
  
    call lstimer('den-cof',te,ts,lupri)
  ! Three-center contributions
    Kfull => Kfull_3cContrib(:,:,:,:,idmat:idmat)
    DfullAO => Dfull(:,:,idmat:idmat)
    call pariK_ThreeCenter(Kfull,DfullAO,setting,calpha_ab,dalpha_ad,&
       &                   iAtomA,orbitalInfo,nBastReg,nAtoms,regCSfull,auxCSfull,lupri,luerr)
    call lstimer('3-CENT ',te,ts,lupri)
  
  ! Two-center contribution
    Kfull => Kfull_2cContrib(:,:,:,:,idmat:idmat)
    call pariK_TwoCenter(Kfull,DfullAO,calpha_ab,dalpha_ad,alpha_beta,&
       &                 iAtomA,orbitalInfo,nBastReg,nAtoms,lupri)
    call lstimer('2-CENT ',te,ts,lupri)
  
    call mem_dealloc(dalpha_ad)
  ENDDO ! --- end loop A
  
  ! --- SYMMETRIZE or ANTI-SYMMETRIZE K, and MULTIPLY BY (1/2)
  ! --- add the contributions of the K matrix (K_ac = K_ac + K_ca = K_ca )
  Kfull_3cContrib(:,:,1,1,idmat) = Kfull_3cContrib(:,:,1,1,idmat) + Kfull_2cContrib(:,:,1,1,idmat)
  Kfull_3cContrib(:,:,1,1,idmat) = Kfull_3cContrib(:,:,1,1,idmat) + transpose(Kfull_3cContrib(:,:,1,1,idmat))
ENDDO
  
call freePariCoefficients(calpha_ab,orbitalInfo%nAtoms)
deallocate(calpha_ab)

CALL freeMolecularOrbitalInfo(orbitalInfo)

! --- ADD THE EXCHANGE CONTRIBUTION (K) TO THE FOCK MATRIX (F)

IF(matrix_type .EQ. mtype_unres_dense)THEN
   DO idmat=1,ndmat
     call DAXPY(F(1)%nrow*F(1)%ncol,-Setting%Scheme%exchangeFactor,Kfull_3cContrib(:,:,1,1,2*idmat-1),1,F(idmat)%elms,1)
     call DAXPY(F(1)%nrow*F(1)%ncol,-Setting%Scheme%exchangeFactor,Kfull_3cContrib(:,:,1,1,2*idmat  ),1,F(idmat)%elmsb,1)
   ENDDO
ELSE !CLOSED_SHELL
   DO idmat=1,ndmat
     call mat_set_from_full(Kfull_3cContrib(:,:,1,1,idmat),-Setting%Scheme%exchangeFactor,F(idmat),'exchange')
   ENDDO
ENDIF

! --- Free temporary allocd memory
!call ls_freeDfull(setting)
call mem_dealloc(Kfull_3cContrib)
call mem_dealloc(Kfull_2cContrib)
call mem_dealloc(Dfull)
call mem_dealloc(alpha_beta)

!call free_AtomSparseMat(alphaBeta)
CALL LSTIMER('PARI-K',tefull,tsfull,lupri)

END SUBROUTINE II_get_pari_df_exchange_mat

!> \brief Calculates the coulomb matrix using regular density fitting
!> \author S. Reine
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param F the coulomb matrix
SUBROUTINE II_get_regular_df_coulomb_mat(LUPRI,LUERR,SETTING,D,F,ndmat)
!use linsolvdf
IMPLICIT NONE
TYPE(MATRIX),target   :: D(ndmat),F(ndmat)
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,ndmat
!
Integer              :: ndmattmp,nbasis,naux,info,nbast,idmat
Type(Matrix)         :: alphabeta
Type(Matrix)         :: galpha
Type(Matrix),target  :: calpha
Type(Matrix),target  :: calpha_old
Real(realk),pointer :: DfullRHS(:,:,:)
Real(realk),pointer :: Ffull(:,:,:)
Real(realk)         :: TSTART,TEND
Character(80)       :: Filename
logical :: inc_scheme,do_inc
integer :: nmat,nrow,ncol
integer :: i,j,natoms
real(realk) :: tmp_sum
Real(realk),pointer   :: eigValphaBeta(:), copy_alpBeta(:,:)
Real(realk)         :: minEigv,maxEigv,conditionNum

IF (matrix_type .EQ. mtype_unres_dense)THEN
  IF (SETTING%SCHEME%FMM) call lsquit('Not allowed combination in II_get_regular_df_coulomb_mat. FMM and unrestricted',-1)
  nmat = 2*ndmat
  nrow = D(1)%nrow
  ncol = D(1)%ncol
  IF ((F(1)%nrow.NE.nrow).OR.(F(1)%ncol.NE.ncol)) CALL LSQUIT('Error in II_get_regular_df_coulomb_mat F/D',-1)
  call mem_alloc(DfullRHS,nrow,ncol,nmat)
  call mem_alloc(Ffull,nrow,ncol,nmat)
  DO Idmat=1,ndmat
    CALL DCOPY(nrow*ncol,D(idmat)%elms,1,DfullRHS(:,:,2*idmat-1),1)
    CALL DCOPY(nrow*ncol,D(idmat)%elmsb,1,DfullRHS(:,:,2*idmat),1)
  ENDDO
ELSE
  nmat = ndmat
ENDIF

!call LSHEADER(lupri,'II_get_regular_df_coulomb_mat')
call LSTIMER('START ',TSTART,TEND,LUPRI)
!set threshold 

nbast = D(1)%nrow


do Idmat=1,nmat
   ndmattmp = 1
   SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%J_THR
   call getMolecularDimensions(SETTING%MOLECULE(1)%p,nAtoms,nBasis,nAux)
   IF(nbast.NE.nbasis)THEN
      call lsquit('ERROR in dimensions in II_get_regular_df_coulomb_mat',-1)
   ENDIF
   !(alpha|rho)
  
   IF(SETTING%SCHEME%FMM)THEN
      CALL ls_attachDmatToSetting(D(idmat:idmat),ndmattmp,setting,'LHS',1,2,.TRUE.,lupri)
      CALL ls_attachDmatToSetting(D(idmat:idmat),ndmattmp,setting,'RHS',3,4,.TRUE.,lupri)
      call ls_multipolemoment(LUPRI,LUERR,SETTING,nbasis,naux,&
        & nbast,nbast,nbast,nbast,&
        &AODFdefault,AOempty,AORdefault,AORdefault,RegularSpec,ContractedInttype,.TRUE.)
      CALL ls_freeDmatFromSetting(setting)
   ENDIF

   IF (matrix_type .EQ. mtype_unres_dense)THEN
     CALL ls_attachDmatToSetting(DfullRHS(:,:,idmat:idmat),nrow,ncol,ndmattmp,setting,'RHS',3,4,lupri)
   ELSE
     CALL ls_attachDmatToSetting(D(idmat:idmat),ndmattmp,setting,'RHS',3,4,.TRUE.,lupri)
   ENDIF
   call mat_init(galpha,naux,ndmattmp)
   call initIntegralOutputDims(setting%Output,naux,1,1,1,ndmattmp)
   call ls_jengine(AODFdefault,AOempty,AORdefault,AORdefault,CoulombOperator,RegularSpec,&
        &          ContractedInttype,SETTING,LUPRI,LUERR)
   CALL retrieve_Output(lupri,setting,galpha,.FALSE.)
   IF(SETTING%SCHEME%FMM)THEN
      IF (ndmattmp.GT.1) CALL LSQUIT('not testet',-1)
      call ls_jengineClassicalMAT(galpha,AODFdefault,AOempty,AORdefault,AORdefault,&
           &CoulombOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
   ENDIF

   CALL ls_freeDmatFromSetting(setting)
   call LSTIMER('GALPHA',TSTART,TEND,LUPRI)
   !(alpha|beta)
   call mat_init(alphabeta,naux,naux)
   call io_get_filename(Filename,'ALBE',AODFdefault,AOEmpty,AODFdefault,AOEmpty,0,0,&
        &CoulombOperator,Contractedinttype,.FALSE.,LUPRI,LUERR)

!  Either calculate the (alpha|beta) matrix and its cholesky-factorization and write
!  OR read cholesky-factorization from file
   IF (io_file_exist(Filename,SETTING%IO)) THEN
      call io_read_mat(alphabeta,Filename,SETTING%IO,LUPRI,LUERR)
   ELSE
      call mat_zero(alphabeta)
      call initIntegralOutputDims(setting%output,naux,1,naux,1,1)
      call ls_getIntegrals(AODFdefault,AOempty,&
           &AODFdefault,AOempty,CoulombOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
      call retrieve_Output(lupri,setting,alphabeta,.FALSE.)

#if 0
     ! checking eigenvalues of the full (alpha|beta) matrix
     call mem_alloc(copy_alpBeta,nAux,nAux)
!      DO j=i,nAux
!         DO i=1,nAux
!            copy_alpBeta(i,j) = alphabeta(i,j)
!         ENDDO
!      ENDDO
     call mat_to_full(alphabeta,1E0_realk,copy_alpBeta)
     call mem_alloc(eigValphaBeta,nAux)
     call II_get_eigv_square_mat(lupri,luerr,copy_alpBeta,eigValphaBeta,nAux) 
     call check_min_max_Array_elem(minEigV,maxEigV,conditionNum,eigValphaBeta,nAux,lupri,luerr)
     call mem_dealloc(eigValphaBeta)
     call mem_dealloc(copy_alpBeta)
     write(lupri,*) "(alpha|beta) full: minEigV of all (alpha|beta) matrix: ",minEigV
     write(*,*)     "(alpha|beta) full: minEigV of all (alpha|beta) matrix: ",minEigV
     write(lupri,*) "(alpha|beta) full: maxEigV of all (alpha|beta) matrix: ",maxEigV
     write(*,*)     "(alpha|beta) full: maxEigV of all (alpha|beta) matrix: ",maxEigV
     write(lupri,*) "(alpha|beta) full: Condition Number (abs(max)/abs(min): ",conditionNum
     write(*,*)     "(alpha|beta) full: Condition Number (abs(max)/abs(min): ",conditionNum
     !call LSQUIT('Testing eigenvalues of the FULL (alpha|beta) matrix - quitting II_get_regular_df_coulomb_mat()',-1)
#endif


!     Make Choleksy-factorization
      CALL mat_dpotrf(alphabeta)
!     Save Cholesky-factors to file
      call io_add_filename(SETTING%IO,Filename,LUPRI)
      call io_write_mat(alphabeta,Filename,SETTING%IO,LUPRI,LUERR)

!TODO in case of SCALAPACK save alphabeta in memory instead of file. TK
!     add a SaveAB_DF logical and so on like CABS.      
   ENDIF
   call LSTIMER('ALBE  ',TSTART,TEND,LUPRI)

   !c_alpha = (alpha|beta)^-1 (beta|rho)
   call mat_init(calpha,naux,ndmattmp)
   CALL mat_assign(calpha,galpha)

   !  Solve the system A*X = B, using the Cholesky-factorization if A and overwriting B with X.
   CALL mat_dpotrs(alphabeta,calpha)

   call mat_free(alphabeta)
   call mat_free(galpha)

   
!  Hack - turn off for unrestricted case (which cannot handle non-square matrices)
   IF (matrix_type .NE. mtype_unres_dense)THEN
      ! **** Write calpha to file in matrix format for storage (used for gradients)
      write(Filename,'(A8,I3)') 'LSCALPHA',idmat
      IF (.NOT.io_file_exist(Filename,SETTING%IO)) THEN
         call io_add_filename(SETTING%IO,Filename,LUPRI)
      ELSE
         IF(setting%scheme%incremental) THEN
            call mat_init(calpha_old,naux,1)
            call io_read_mat(calpha_old,Filename,Setting%IO,lupri,luerr)
            call mat_daxpy(1E0_realk,calpha_old,calpha)
            call mat_free(calpha_old)
         ENDIF
      ENDIF
      call io_write_mat(calpha,Filename,Setting%IO,lupri,luerr)
      ! **** End write
   ENDIF
   
   call LSTIMER('LINSOL',TSTART,TEND,LUPRI)
   !Jfit_ab = (ab|rho_fit)
   IF(SETTING%SCHEME%FMM) THEN
      CALL ls_attachDmatToSetting(D(idmat:idmat),ndmattmp,setting,'LHS',1,2,.TRUE.,lupri)
      CALL ls_attachDmatToSetting(calpha,ndmattmp,setting,'RHS',3,4,.TRUE.,lupri)
      call ls_multipolemoment(LUPRI,LUERR,SETTING,nbasis,naux,&
        & nbast,nbast,naux,1,&
        &AORdefault,AORdefault,AODFdefault,AOempty,RegularSpec,ContractedInttype,.TRUE.)
      CALL ls_freeDmatFromSetting(setting)
   ENDIF
   CALL ls_attachDmatToSetting(calpha,ndmattmp,setting,'RHS',3,4,.TRUE.,lupri)
   call initIntegralOutputDims(setting%Output,nbast,nbast,1,1,1)
   call ls_jengine(AORdefault,AORdefault,AODFdefault,AOempty,CoulombOperator,RegularSpec,ContractedInttype,&
        &          SETTING,LUPRI,LUERR)
   IF (matrix_type .EQ. mtype_unres_dense)THEN
     CALL retrieve_Output(lupri,setting,Ffull(:,:,idmat),setting%IntegralTransformGC)
   ELSE
     CALL retrieve_Output(lupri,setting,F(idmat),setting%IntegralTransformGC)
   ENDIF
   IF(SETTING%SCHEME%FMM)THEN
      IF (ndmat.GT.1) CALL LSQUIT('not testet',-1)
      call ls_jengineClassicalMAT(F(1),AORdefault,AORdefault,AODFdefault,AOempty,&
           & CoulombOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
   ENDIF
   CALL ls_freeDmatFromSetting(setting)
   call mat_free(calpha)
ENDDO

IF (matrix_type .EQ. mtype_unres_dense)THEN
  DO Idmat=1,ndmat
    CALL DCOPY(nrow*ncol,Ffull(:,:,2*idmat-1),1,F(idmat)%elms, 1)
    CALL DCOPY(nrow*ncol,Ffull(:,:,2*idmat  ),1,F(idmat)%elmsb,1)
  ENDDO
  call mem_dealloc(DfullRHS)
  call mem_dealloc(Ffull)
ENDIF

call LSTIMER('FIT-J ',TSTART,TEND,LUPRI)
END SUBROUTINE II_get_regular_df_coulomb_mat

!> \brief Calculates the exchange matrix using regular density fitting
!> \author S. Reine and T. Kjaergaard
!> \date 2012-03-16
!> \param lupri Default print-unit for output
!> \param luerr Default print-unit for termination
!> \param setting Contains information about the integral settings
!> \param D Density-matrix
!> \param F Exchange contribution to the Fock- or KS-matrix
!> \param ndmat number of Density matrices
SUBROUTINE II_get_df_exchange_mat(LUPRI,LUERR,SETTING,Dmat,F,ndmat)
  IMPLICIT NONE
  TYPE(MATRIX),target,intent(in)    :: Dmat(ndmat)
  TYPE(MATRIX),target,intent(inout) :: F(ndmat)
  TYPE(LSSETTING),intent(inout)     :: SETTING
  Integer,intent(in)                :: LUPRI,LUERR,ndmat
  !
  Integer                    :: nAtoms,nBastAux,nBast,nbasis
  Real(realk),pointer        :: Dfull(:,:,:)
  Real(realk),pointer        :: Kfull(:,:,:)
  TYPE(MOLECULARORBITALINFO) :: orbitalInfo
  Integer                    :: iAtomA,iAtomB,iAtomC,iAtomD
  real(realk),pointer        :: Galpha(:,:,:),Calpha(:,:,:),Dalpha(:,:,:),alphabeta(:,:)
  real(realk),pointer        :: alphaAB(:,:,:),alphaCD(:,:,:)
  Integer                    :: nBastLocA,nBastLocC,nAuxA,nAuxC
  Integer                    :: startRegA,startRegC,startAuxA,startAuxC
  Integer                    :: endRegA,endRegC,endAuxA,endAuxC
  Integer                    :: iRegA,iRegB,iRegC,iRegD,INFO
  TYPE(MoleculeInfo),pointer      :: molecule
  TYPE(MoleculeInfo),pointer      :: ATOMS(:)
  Character(80)             :: Filename
  Real(realk) :: TSTART,TEND,TSTARTFULL,TENDFULL,tmp
  integer :: i,idmat,nmat,a,c,D,ialpha
  integer(kind=long) :: nsize,maxsize
  integer(kind=short),pointer :: AtomicGab(:,:)

  call LSTIMER('START ',TSTARTFULL,TENDFULL,LUPRI)
  call LSTIMER('START ',TSTART,TEND,LUPRI)
  nbasis = Dmat(1)%nrow
  maxsize = 1000*1000*1000
  IF(ndmat.GT.1.AND.setting%scheme%incremental)call lsquit('increm and df ndmat>1 not working',-1)
  call getMolecularDimensions(SETTING%MOLECULE(1)%p,nAtoms,nBast,nBastAux)
  IF(nbasis.NE.nBast)THEN
     print*,'nbasis',nbasis
     print*,'nbast',nbast
     Call lsquit('dim mismatch in II_get_df_exchange2',-1)
  ENDIF
  IF (ndmat.GT.1) THEN
      WRITE(*,*)     "The PARI approximation isn't implemented for unrestricted cases yet."
      WRITE(LUPRI,*) "The PARI approximation isn't implemented for unrestricted cases yet."
     CALL LSQUIT('Error in II_get_pari_df_exchange_mat: ndmat>1 not tested',-1)
  ENDIF
  !set threshold 
  SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%K_THR
  
  IF (matrix_type .EQ. mtype_unres_dense)THEN
     nmat = 2*ndmat
     call mem_alloc(Dfull,nbast,nbast,nmat)
     DO Idmat=1,ndmat
        CALL DCOPY(nbast*nbast,Dmat(idmat)%elms,1,Dfull(:,:,2*idmat-1),1)
        CALL DCOPY(nbast*nbast,Dmat(idmat)%elmsb,1,Dfull(:,:,2*idmat),1)
     ENDDO
  ELSE
     nmat = ndmat
     call mem_alloc(Dfull,nbast,nbast,nmat)
     DO Idmat=1,ndmat
        call mat_to_full(Dmat(idmat),1.0d0,Dfull(:,:,idmat))
     ENDDO
  ENDIF
  call mem_alloc(Kfull,nbast,nbast,nmat)
  call ls_dzero(Kfull,nbast*nbast*nmat)
  molecule => SETTING%MOLECULE(1)%p
  nsize = 3*nBast*nBast
  nsize = nsize*nBastAux
  nsize = nsize*8
  IF(nsize.LT.maxsize)THEN !.OR.DEBUGDFKATOMLOOP
     !(alpha|cd)
     !call typedef_setMolecules(setting,molecule,1,2,3,4)
     call initIntegralOutputDims(setting%Output,nBastAux,1,nbast,nbast,1)
     call ls_getIntegrals(AODFdefault,AOempty,AORdefault,AORdefault,CoulombOperator,RegularSpec,&
          &               ContractedInttype,SETTING,LUPRI,LUERR)
     call mem_alloc(Galpha,nbastaux,nBast,nBast)
     CALL retrieve_Output(lupri,setting,Galpha,.FALSE.)
     call LSTIMER('DF-K-INT',TSTART,TEND,LUPRI)

     !(alpha|beta)
     call io_get_filename(Filename,'ALBE',AODFdefault,AOEmpty,AODFdefault,AOEmpty,0,0,&
          &CoulombOperator,Contractedinttype,.FALSE.,LUPRI,LUERR)
     IF (io_file_exist(Filename,SETTING%IO)) THEN
        call mem_alloc(alphabeta,nbastaux,nbastaux)
        call io_read(alphabeta,nbastaux,nbastaux,1,1,1,Filename,SETTING%IO,LUPRI,LUERR)
     ELSE
        call initIntegralOutputDims(setting%output,nbastaux,1,nbastaux,1,1)
        call ls_getIntegrals(AODFdefault,AOempty,AODFdefault,AOempty,CoulombOperator,RegularSpec,&
             &                  ContractedInttype,SETTING,LUPRI,LUERR)
        call mem_alloc(alphabeta,nbastaux,nbastaux)
        call retrieve_Output(lupri,setting,alphabeta,.FALSE.)

        !  Make Choleksy-factorization
        call DPOTRF('U',nbastaux,alphabeta,nbastaux,INFO)
        IF (info.ne. 0) THEN
           WRITE(LUPRI,'(1X,A,I5)') 'DPOTRF error in II_get_df_exchange_mat. Info =',info
           call LSQUIT('DPOTRF error in II_get_df_exchange_mat',lupri)
        ENDIF
        !  Save Cholesky-factors to file
        call io_add_filename(SETTING%IO,Filename,LUPRI)
        call io_write(alphabeta,nbastaux,nbastaux,1,1,1,Filename,SETTING%IO,LUPRI,LUERR)
     ENDIF
     call LSTIMER('ALBE-K',TSTART,TEND,LUPRI)

     !c_alpha = (alpha|beta)^-1 (beta|cd)
     call mem_alloc(Calpha,nBastaux,nBast,nBast)
     call DCOPY(nBastaux*nBast*nBast,Galpha,1,Calpha,1)
     !  c_(alpha,aB) = (alpha|beta)^-1 (beta|aB)
     !  Solve the system A*X = B, overwriting B with X.
     !  Solve  A=(beta|alpha)  X=c_alpha   B = (beta|aB)
     CALL DPOTRS('U',nBastAux,nBast*nBast,alphabeta,nbastaux,Calpha,nbastaux,info)
     IF (info.ne. 0) THEN
        WRITE(LUPRI,'(1X,A,I5)') 'DPOTRS error in II_get_df_exchange_mat. Info =',info
        call LSQUIT('DPOTRS error in II_get_df_exchange_mat',lupri)
     ENDIF
     call mem_dealloc(alphabeta)

     call LSTIMER('LINSOL-K',TSTART,TEND,LUPRI)
     call mem_alloc(Dalpha,nBastaux,nBast,nBast)
     DO idmat=1,nmat
        ! dalpha,ad = sum_d c_alpha,ab D_bd
        CALL DGEMM('N','N',nBast*nbastaux,nBast,nBast,1E0_realk,calpha,nBast*nbastaux,&
             &         Dfull(:,:,idmat),nBast,0E0_realk,dalpha,nBast*nbastaux)

        !Kac = sum_alpha,d dalpha,ad * (alpha|cd) 
        DO i=1,nBast
           CALL DGEMM('T','N',nBast,nBast,nBastaux,-setting%scheme%exchangeFactor,Dalpha(:,:,i),nBastaux,&
                &           galpha(:,:,i),nBastaux,1E0_realk,Kfull(:,:,idmat),nBast)
        ENDDO
     ENDDO

     call mem_dealloc(galpha)
     call mem_dealloc(calpha)
     call mem_dealloc(Dalpha)
     call LSTIMER('FULL FIT-K ',TSTART,TEND,LUPRI)
  ELSE
     !MEMORY OPTIMIZED VERSION
     call getMolecularDimensions(molecule,nAtoms,nBast,nBastAux)
     ! --- Build an array with the list of atoms with their respective nb. of orbitals/auxiliary functions
     ! --- Build another array to know where each atom start in the complete matrices
     call setMolecularOrbitalInfo(molecule,orbitalInfo)

     !find nMaxLocBast
!     IF(nMaxLocBast*nBast*nBastAux*2*8+nBastAux*nBastAux*8.GT.1000*1000*1000)THEN
!        CALL LSQUIT('MEMORY REQUIREMENTS TOO LARGE',-1)
!     ENDIF     
     
!     Constrict Atomic Gab matrix
!     call mem_alloc(AtomicGab,nAtoms,nAtoms)
!     Call ls_BuildAtomicGab(AtomicGab,nAtoms,nAtoms,AORdefault,AORdefault,&
!          & CoulombOperator,SETTING,LUPRI,LUERR)
!     WRITE(lupri,*)'AtomicGab'
!     call shortint_output(AtomicGab,nAtoms,nAtoms,lupri)
!     call mem_dealloc(AtomicGab)

     !(alpha|beta)
     call io_get_filename(Filename,'ALBE',AODFdefault,AOEmpty,AODFdefault,AOEmpty,0,0,&
          &CoulombOperator,Contractedinttype,.FALSE.,LUPRI,LUERR)
     IF (io_file_exist(Filename,SETTING%IO)) THEN
        call mem_alloc(alphabeta,nbastaux,nbastaux)
        call io_read(alphabeta,nbastaux,nbastaux,1,1,1,Filename,SETTING%IO,LUPRI,LUERR)
     ELSE
        call initIntegralOutputDims(setting%output,nbastaux,1,nbastaux,1,1)
        call ls_getIntegrals(AODFdefault,AOempty,AODFdefault,AOempty,CoulombOperator,RegularSpec,&
             &                  ContractedInttype,SETTING,LUPRI,LUERR)
        call mem_alloc(alphabeta,nbastaux,nbastaux)
        call retrieve_Output(lupri,setting,alphabeta,.FALSE.)
        !  Make Choleksy-factorization
        call DPOTRF('U',nbastaux,alphabeta,nbastaux,INFO)
        IF (info.ne. 0) THEN
           WRITE(LUPRI,'(1X,A,I5)') 'DPOTRF error in II_get_df_exchange_mat. Info =',info
           call LSQUIT('DPOTRF error in II_get_df_exchange_mat',lupri)
        ENDIF
        !  Save Cholesky-factors to file
        call io_add_filename(SETTING%IO,Filename,LUPRI)
        call io_write(alphabeta,nbastaux,nbastaux,1,1,1,Filename,SETTING%IO,LUPRI,LUERR)
     ENDIF
     call LSTIMER('ALBE-K',TSTART,TEND,LUPRI)

     allocate(ATOMS(nAtoms))
     CALL pari_set_atomic_fragments(molecule,ATOMS,nAtoms,lupri)
     DO idmat=1,nmat
        !Modify loop to be over possible batch of atoms
        DO iAtomA=1,nAtoms
           call getAtomicOrbitalInfo(orbitalInfo,iAtomA,nBastLocA,startRegA,endRegA,nauxA,startAuxA,endAuxA)
           ! List B should be from Gab list
           ! create molecule with these atoms
           ! (alpha| a fullB)
           call typedef_setMolecules(setting,molecule,1,4,ATOMS(iAtomA),3)
           call initIntegralOutputDims(setting%Output,nBastAux,1,nBastLocA,nBast,1)
           call ls_getIntegrals(AODFdefault,AOEmpty,AORdefault,AORdefault,CoulombOperator,RegularSpec,&
                &               Contractedinttype,SETTING,LUPRI,LUERR)
           call mem_alloc(alphaAB,nBastAux,nBastLocA,nBast)
           CALL retrieve_Output(lupri,setting,alphaAB,.FALSE.)

           !  c_(alpha,aB) = (alpha|beta)^-1 (beta|aB)
           !  Solve the system A*X = B, overwriting B with X.
           !  Solve  A= (beta|alpha)  X=c_alpha   B = (beta|aB)

           call mem_alloc(Calpha,nBastAux,nBastLocA,nBast)
           call DCOPY(nBastaux*nBast*nBastLocA,alphaAB,1,Calpha,1)
           CALL DPOTRS('U',nBastAux,nbast*nBastLocA,alphabeta,nBastAux,Calpha,nBastAux,info)

           ! d_(alpha,aD) = sum_b c_(alpha,aB) * D_bd = [alpha*s,B]*[B,D] 
           call mem_alloc(Dalpha,nbastAux,nBastLocA,nBast)  
           CALL DGEMM('N','N',nBastAux*nBastLocA,nBast,nBast,1E0_realk,Calpha,nBastAux*nBastLocA,&
                &         Dfull(:,:,idmat),nbast,0E0_realk,dalpha,nBastAux*nBastLocA)
           call mem_dealloc(Calpha)

           ! --- Loop over atoms C
           ! BATCHES OF ATOMS so that nBastA*nBast*nBastAux less than 1 GB
           iAtomC=iAtomA
           !   Generate K(A,C) += d_(alpha,aD)*(alpha|cd)
           DO a=1,nBastLocA
              DO c=1,nBastLocA
                 tmp = 0.0E0_realk
                 do D=1,nBast 
                    DO ialpha=1,nBastAux
                       tmp = tmp + dalpha(ialpha,a,D)*alphaAB(ialpha,c,D)    
                    ENDDO
                 ENDDO
                 Kfull(startRegA-1+a,startRegA-1+c,idmat) = &
 & Kfull(startRegA-1+a,startRegA-1+c,idmat) -setting%scheme%exchangeFactor*tmp
              ENDDO
           ENDDO
           call mem_dealloc(alphaAB)

           DO iAtomC=iAtomA+1,nAtoms 
              call getAtomicOrbitalInfo(orbitalInfo,iAtomC,nBastLocC,startRegC,endRegC,nauxC,startAuxC,endAuxC)
              !ThreeCenter Contributions
              ! (alpha| c D)
              call typedef_setMolecules(setting,molecule,1,4,ATOMS(iAtomC),3)
              call initIntegralOutputDims(setting%Output,nBastAux,1,nBastLocC,nBast,1)
              call ls_getIntegrals(AODFdefault,AOEmpty,AORdefault,AORdefault,CoulombOperator,&
                   & RegularSpec,Contractedinttype,SETTING,LUPRI,LUERR)
              call mem_alloc(alphaCD,nBastAux,nBastLocC,nBast)
              CALL retrieve_Output(lupri,setting,alphacD,.FALSE.)

              !   Generate K(A,C) += d_(alpha,aD)*(alpha|cd)
              DO a=1,nBastLocA
                 DO c=1,nBastLocC
                    tmp = 0.0E0_realk
                    do D=1,nBast 
                       DO ialpha=1,nBastAux
                          tmp = tmp + dalpha(ialpha,a,D)*alphacD(ialpha,c,D)    
                       ENDDO
                    ENDDO
                    !WARNING Assumes D symmetric
                    Kfull(startRegA-1+a,startRegC-1+c,idmat) = &
 & Kfull(startRegA-1+a,startRegC-1+c,idmat) -setting%scheme%exchangeFactor*tmp
                    Kfull(startRegC-1+c,startRegA-1+a,idmat) = &
 & Kfull(startRegC-1+c,startRegA-1+a,idmat) -setting%scheme%exchangeFactor*tmp
                 ENDDO
              ENDDO
              call mem_dealloc(alphaCD)
           ENDDO
        ENDDO
        call mem_dealloc(dalpha)  
     ENDDO
     CALL freeMolecularOrbitalInfo(orbitalInfo)
     call mem_dealloc(alphabeta)
     call typedef_setMolecules(setting,molecule,1,2,3,4)
     deallocate(ATOMS)
     call LSTIMER('FRAG FIT-K ',TSTART,TEND,LUPRI)
  ENDIF
  call mem_dealloc(Dfull)  
  IF (matrix_type .EQ. mtype_unres_dense)THEN
     DO Idmat=1,ndmat
        CALL DCOPY(nbast*nbast,Kfull(:,:,2*idmat-1),1,F(idmat)%elms, 1)
        CALL DCOPY(nbast*nbast,Kfull(:,:,2*idmat  ),1,F(idmat)%elmsb,1)
     ENDDO
  ELSE
     DO Idmat=1,ndmat
        call mat_set_from_full(Kfull(:,:,idmat),1E0_realk,F(idmat))
     ENDDO
  ENDIF
  call mem_dealloc(Kfull)
  call LSTIMER('FIT-K ',TSTARTFULL,TENDFULL,LUPRI)

END SUBROUTINE II_get_df_exchange_mat

!> \brief Calculates the coulomb matrix using overlap density fitting
!> \author S. Reine
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param F the coulomb matrix
SUBROUTINE II_get_overlap_df_coulomb_mat(LUPRI,LUERR,SETTING,D,F)
IMPLICIT NONE
TYPE(MATRIX),target   :: D,F
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
!
Integer             :: ndmat,nbasis,naux
Real(realk),pointer :: g1alphafull(:,:,:,:,:)
!Real(realk),pointer :: g2alpha(:,:,:,:,:)
Real(realk),pointer :: calphafull(:,:,:)
TYPE(matrix),target :: g1alpha,calpha,g2alpha,F2
Real(realk)         :: TSTART,TEND
integer             :: oper,natoms
Logical             :: Coulomb
Logical             :: Frag
type(matrixp)       :: Jmat(1),Dmat(1)
type(matrix),target :: tmpF,tmpD

!CALL LSHEADER(lupri,'II_get_overlap_df_coulomb_mat')
CALL LSTIMER('START ',TSTART,TEND,LUPRI)
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%J_THR
oper=OverlapOperator
Coulomb = oper.EQ.CoulombOperator
ndmat = 1
call getMolecularDimensions(SETTING%MOLECULE(1)%p,nAtoms,nBasis,nAux)

!<alpha|w|rho>
CALL ls_attachDmatToSetting(D,ndmat,setting,'LHS',1,2,.TRUE.,lupri)
CALL ls_attachDmatToSetting(D,ndmat,setting,'RHS',3,4,.TRUE.,lupri)
IF(SETTING%SCHEME%FMM.AND.Coulomb) call ls_multipolemoment(LUPRI,LUERR,SETTING,nbasis,naux,&
     & D%nrow,D%ncol,D%nrow,D%ncol,&
     & AODFdefault,AOempty,AORdefault,AORdefault,RegularSpec,ContractedInttype,.TRUE.)
CALL ls_freeDmatFromSetting(setting)
CALL ls_attachDmatToSetting(D,ndmat,setting,'RHS',3,4,.TRUE.,lupri)
!CALL ls_attachDmatToSetting(Dmat,ndmat,setting,'RHS',3,4,lupri)
!Jmat,ndmat,
call initIntegralOutputDims(setting%Output,naux,1,1,1,1)
call ls_jengine(AODFdefault,AOempty,AORdefault,AORdefault,oper,RegularSpec,ContractedInttype,&
     &          SETTING,LUPRI,LUERR)
call mat_init(g1alpha,naux,1) 
CALL mat_ZERO(g1alpha)
CALL retrieve_Output(lupri,setting,g1alpha,.FALSE.)
IF(SETTING%SCHEME%FMM)THEN
   call ls_jengineClassicalMAT(g1alpha,AODFdefault,AOempty,AORdefault,AORdefault,&
        & oper,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
ENDIF
   
CALL ls_freeDmatFromSetting(setting)

call mem_alloc(g1alphafull,naux,1,1,1,ndmat)
call mat_to_full(g1alpha,1E0_realk,g1alphafull(:,:,1,1,1))

call mat_free(g1alpha)
call mem_alloc(calphafull,naux,1,ndmat)
!|rho_fit) = <alpha|w|beta>^-1 <beta|w|rho>
call linsolv_df(calphafull,g1alphafull,AODFdefault,oper,naux,ndmat,SETTING,LUPRI,LUERR)
call mem_dealloc(g1alphafull)
call mat_init(calpha,naux,1)
CALL mat_set_from_full(calphafull(:,:,1),1E0_realk,calpha)
call mem_dealloc(calphafull)
!Calculate multipole moments (lor read from file) if using FMM
CALL ls_attachDmatToSetting(D,ndmat,setting,'LHS',1,2,.TRUE.,lupri)
CALL ls_attachDmatToSetting(calpha,ndmat,setting,'RHS',3,4,.TRUE.,lupri)
IF(SETTING%SCHEME%FMM)call ls_multipolemoment(LUPRI,LUERR,SETTING,nbasis,naux,&
     & D%nrow,D%ncol,calpha%nrow,calpha%ncol,&
     &AORdefault,AORdefault,AODFdefault,AOempty,RegularSpec,ContractedInttype,.TRUE.)
CALL ls_freeDmatFromSetting(setting)
CALL ls_attachDmatToSetting(calpha,ndmat,setting,'RHS',3,4,.TRUE.,lupri)
!Jfit_ab(1) = (ab|rho_fit)
CALL MAT_ZERO(F)
!Jmat(1)%p => F
!Dmat(1)%p => calpha
!CALL ls_attachDmatToSetting(Dmat,ndmat,setting,'RHS',3,4,lupri)
!Jmat,ndmat,
call initIntegralOutputDims(setting%Output,F%nrow,F%ncol,1,1,1)
call ls_jengine(AORdefault,AORdefault,AODFdefault,AOempty,CoulombOperator,RegularSpec,ContractedInttype,&
     &          SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,F,.FALSE.)
IF(SETTING%SCHEME%FMM)THEN
   call ls_jengineClassicalMat(F,AORdefault,AORdefault,AODFdefault,AOempty,&
        & CoulombOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
ENDIF
CALL ls_freeDmatFromSetting(setting)

if(.NOT.Coulomb) then
   !(alpha|rho)
   CALL ls_attachDmatToSetting(D,ndmat,setting,'LHS',1,2,.TRUE.,lupri)
   CALL ls_attachDmatToSetting(calpha,ndmat,setting,'RHS',3,4,.TRUE.,lupri)
   IF(SETTING%SCHEME%FMM) call ls_multipolemoment(LUPRI,LUERR,SETTING,nbasis,naux,&
     & D%nrow,D%ncol,calpha%nrow,calpha%ncol,&
     &AODFdefault,AOempty,AORdefault,AORdefault,RegularSpec,ContractedInttype,.TRUE.)
   CALL ls_freeDmatFromSetting(setting)
   CALL ls_attachDmatToSetting(D,ndmat,setting,'RHS',3,4,.TRUE.,lupri)
   !Jmat,ndmat,
   call initIntegralOutputDims(setting%Output,naux,1,1,1,1)
   call ls_jengine(AODFdefault,AOempty,AORdefault,AORdefault,CoulombOperator,RegularSpec,ContractedInttype,&
        &          SETTING,LUPRI,LUERR)
   call mat_init(g1alpha,naux,1) 
   call mat_zero(g1alpha) 
   CALL retrieve_Output(lupri,setting,g1alpha,.FALSE.)
   IF(SETTING%SCHEME%FMM)THEN
      call ls_jengineClassicalMat(g1alpha,AODFdefault,AOempty,AORdefault,AORdefault,&
           & CoulombOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
   ENDIF

   CALL ls_freeDmatFromSetting(setting)
   !(alpha|rho_fit)
   CALL ls_attachDmatToSetting(D,ndmat,setting,'LHS',1,2,.TRUE.,lupri)
   CALL ls_attachDmatToSetting(calpha,ndmat,setting,'RHS',3,4,.TRUE.,lupri)
   IF(SETTING%SCHEME%FMM)call ls_multipolemoment(LUPRI,LUERR,SETTING,nbasis,naux,&
        & D%nrow,D%ncol,calpha%nrow,calpha%ncol,&
        & AODFdefault,AOempty,AODFdefault,AOempty,RegularSpec,ContractedInttype,.TRUE.)
   CALL ls_freeDmatFromSetting(setting)
   CALL ls_attachDmatToSetting(calpha,ndmat,setting,'RHS',3,4,.TRUE.,lupri)
   call initIntegralOutputDims(setting%Output,naux,1,1,1,1)
   call ls_jengine(AODFdefault,AOempty,AODFdefault,AOempty,CoulombOperator,&
        & RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
   call mat_init(g2alpha,naux,1)
   CALL mat_ZERO(g2alpha)
   CALL retrieve_Output(lupri,setting,g2alpha,.FALSE.)
   IF(SETTING%SCHEME%FMM)THEN
      call ls_jengineClassicalMat(g2alpha,AODFdefault,AOempty,AODFdefault,AOempty&
           & ,CoulombOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
   ENDIF
   CALL ls_freeDmatFromSetting(setting)
   !(alpha|delta_rho)
   call mat_daxpy(-1E0_realk,g2alpha,g1alpha)
   call mat_free(g2alpha)
   !c_alpha = <alpha|w|beta>^-1 (beta|delta_rho)
   call mem_alloc(g1alphafull,naux,1,1,1,ndmat)
   call mat_to_full(g1alpha,1E0_realk,g1alphafull(:,:,1,1,1))
   call mat_free(g1alpha)
   call mem_alloc(calphafull,naux,1,ndmat)
   call linsolv_df(calphafull,g1alphafull,AODFdefault,oper,naux,ndmat,SETTING,LUPRI,LUERR)
   call mem_dealloc(g1alphafull)
   CALL mat_set_from_full(calphafull(:,:,1),1E0_realk,calpha)
   call mem_dealloc(calphafull)

   !Jfit_ab(2,3) = <ab|w|delta_rho>
   CALL ls_attachDmatToSetting(D,ndmat,setting,'LHS',1,2,.TRUE.,lupri)
   CALL ls_attachDmatToSetting(calpha,ndmat,setting,'RHS',3,4,.TRUE.,lupri)
   IF(SETTING%SCHEME%FMM.AND.Coulomb)call ls_multipolemoment(LUPRI,LUERR,SETTING,&
        & nbasis,naux,D%nrow,D%ncol,calpha%nrow,calpha%ncol,&
        &AORdefault,AORdefault,AODFdefault,AOempty,RegularSpec,ContractedInttype,.TRUE.)
   CALL ls_freeDmatFromSetting(setting)
   CALL ls_attachDmatToSetting(calpha,ndmat,setting,'RHS',3,4,.TRUE.,lupri)
   call initIntegralOutputDims(setting%Output,nbasis,nbasis,1,1,1)
   call ls_jengine(AORdefault,AORdefault,AODFdefault,AOempty,oper,RegularSpec,ContractedInttype,&
        &          SETTING,LUPRI,LUERR)
   call mat_init(F2,nbasis,nbasis)
   CALL MAT_ZERO(F2)
   CALL retrieve_Output(lupri,setting,F2,.FALSE.)
   IF(SETTING%SCHEME%FMM)THEN
      call ls_jengineClassicalMat(F2,AORdefault,AORdefault,AODFdefault,AOempty,oper,RegularSpec,&
      &ContractedInttype,SETTING,LUPRI,LUERR)
   ENDIF
   CALL ls_freeDmatFromSetting(setting)
   call mat_daxpy(1E0_realk,F2,F)
   call mat_free(F2)
endif
call mat_free(calpha)

IF(setting%IntegralTransformGC)THEN
   call AO2GCAO_transform_matrixF(F,setting,lupri)
ENDIF

CALL LSTIMER('FIT-JO',TSTART,TEND,LUPRI)

END SUBROUTINE II_get_overlap_df_coulomb_mat

END MODULE INTEGRALINTERFACEMODULEDF
