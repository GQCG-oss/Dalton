MODULE IntegralInterfaceModuleDF
  use precision
  use TYPEDEFTYPE
  use Typedef  
  use Matrix_module
  use Matrix_Operations
  use LSparameters
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
  use dec_typedef_module, only: batchTOorb
  use screen_mod
  use BUILDAOBATCH
  use,intrinsic :: iso_c_binding,only:c_f_pointer, c_loc,C_PTR
  use ThermiteIntTransform_module
  public :: II_get_df_coulomb_mat,II_get_df_J_gradient, &
       & II_get_df_exchange_mat, II_get_pari_df_exchange_mat,&
       & init_IIDF_matrix,free_IIDF_matrix,&
       & II_get_RI_alphaCD_3CenterInt, II_get_RI_alphabeta_2CenterInt, &
       & II_get_RI_alphaCD_3CenterInt2,getRIbasisMPI,getMaxAtomicnAux,&
       & II_get_RI_AlphaCD_3CenterIntFullOnAllNN, GetOperatorFromCharacter,&
       & II_get_RI_AlphaCD_3CenterIntFullOnAllNNdim,II_get_RIMP2_grad,&
       & II_get_RI_AlphaBeta_geoderiv2CenterInt
  private

  SAVE
  logical :: SavealphaBeta
  type(matrix) :: AlphaBetaSave
contains
  subroutine init_IIDF_matrix()
    SavealphaBeta = .FALSE.
  end subroutine init_IIDF_matrix

  subroutine free_IIDF_matrix()
    IF(SavealphaBeta)THEN
       call mat_free(AlphaBetaSave)
    ENDIF
  end subroutine free_IIDF_matrix

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


!> \brief Find the precalculated screening matrices for Reg and Aux basis
!> 
  SUBROUTINE II_getScreenMatFull(regCSfull,auxCSfull,Setting,lupri,luerr)
    implicit none
    TYPE(LSTENSOR),pointer        :: regCSfull,auxCSfull
    TYPE(LSSETTING),intent(inout) :: SETTING
    INTEGER,intent(in)            :: LUPRI,LUERR
    !
    Logical       :: dofit

    IF (setting%scheme%CS_SCREEN) THEN
       call II_getScreenMat(regCSfull,AORdefault,AORdefault,&
            AORdefault,AORdefault,Setting,lupri,luerr)
       
       ! For density fiiting, get the screening matrix for the auxiliary basis
       dofit = setting%scheme%densfit.OR.setting%scheme%pari_j.OR. &
            setting%scheme%pari_k.OR.setting%scheme%mopari_k
       
       IF (dofit) THEN
          call II_getScreenMat(auxCSfull,AODFdefault,AOempty,&
               AODFdefault,AOempty,Setting,lupri,luerr)
       ELSE
          NULLIFY(auxCSfull)
       ENDIF
    ELSE
       NULLIFY(regCSfull)
       NULLIFY(auxCSfull)
    ENDIF

  END SUBROUTINE II_getScreenMatFull

  !> \brief Find the precalculated screening matrix
  !> (up to 10 orders of magnitude smaller than the current intTHRESHOLD)
  !
  !> \param CSfull Precalculated screening matrix
  !> \param AO1 Character string usually 'Regular' or 'Empty' or 'DF-Aux' for center 1
  !> \param AO2 Character string usually 'Regular' or 'Empty' or 'DF-Aux' for center 2
  !> \param AO3 Character string usually 'Regular' or 'Empty' or 'DF-Aux' for center 3
  !> \param AO4 Character string usually 'Regular' or 'Empty' or 'DF-Aux' for center 4
  !> \param SETTING Integral evalualtion settings
  !> \param LUPRI the logical unit number for the output file
  !> \param LUERR the logical unit number for the error file
  subroutine II_getScreenMat(CSfull,AO1,AO2,AO3,AO4,Setting,lupri,luerr)
    implicit none
    TYPE(LSTENSOR),pointer        :: CSfull
    TYPE(LSSETTING),intent(inout) :: SETTING
    INTEGER,intent(in)            :: LUPRI,LUERR,AO1,AO2,AO3,AO4

    type(moleculeinfo),pointer :: molecule
    Character(80) :: Filename
    Character(53) :: identifier
    Logical       :: FoundInMem
    Integer       :: molID,THR,THR2
    Integer       :: i
    THR = ABS(NINT(LOG10(SETTING%SCHEME%intTHRESHOLD)))
       molecule => SETTING%MOLECULE(1)%p
       molID = SETTING%molID(1)

    DO I=0,10
       THR2 = THR + I
       CALL io_get_CSidentifier(identifier,THR2,molecule,molecule,&
            Setting%Scheme%CS_SCREEN,Setting%Scheme%PS_SCREEN)
       CALL io_get_filename(Filename,identifier,AO1,AO2,&
            AO3,AO4,molID,molID,&
            CoulombOperator,ContractedInttype,.FALSE.,LUPRI,LUERR)
       call determine_lst_in_screenlist(Filename,FoundInMem,SETTING%IO)
       IF (FoundInMem) THEN
          call screen_associate(CSfull,Filename,FoundInMem)
          EXIT
       ENDIF
    ENDDO
    IF (.NOT.FoundInMem) THEN
       CALL io_get_CSidentifier(identifier,THR,molecule,molecule,&
            Setting%Scheme%CS_SCREEN,Setting%Scheme%PS_SCREEN)
       CALL io_get_filename(Filename,identifier,AO1,AO2,&
            AO3,AO4,molID,molID,&
            CoulombOperator,ContractedInttype,.FALSE.,LUPRI,LUERR)
       write(lupri,'(1X,1A)') &
            'Error in II_getScreenMatFull no screening matrix found'
       write(lupri,'(1X,2A)') 'Identifier =',identifier
       write(lupri,'(1X,2A)') 'Filename ='  ,Filename
       CALL LSQUIT('Error in II_getScreenMatFull no screening matrix found',-1)
    ENDIF

  end subroutine II_getScreenMat

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
    call Test_if_64bit_integer_required(nBastAux,nBastReg,nBastReg)
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
          call Test_if_64bit_integer_required(nAuxAB,nRegA*nRegB)
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

  !> \brief Calculates the pair-atomic density-fitted (PARI) Exchange matrix
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
    Integer                    :: iAtomA,iAtomB,iAtomC,iAtomD,iAtom
    TYPE(MAT3D),pointer        :: calpha_ab(:)
    Real(realk),pointer        :: dalpha_ad(:,:,:)

    TYPE(MATRIX)               :: matrixK
    Integer                    :: nRegA,nRegB,nRegC,nRegD,nAuxA,nAuxB,nAuxC,nAuxD
    Integer                    :: startRegA,startRegB,startRegC,startRegD
    Integer                    :: startAuxA,startAuxB,startAuxC,startAuxD
    Integer                    :: endRegA,endRegB,endRegC,endRegD
    Integer                    :: endAuxA,endAuxB,endAuxC,endAuxD
    Integer                    :: iAlpha,iRegA,iRegB,iRegC,iRegD,iAuxA
    Integer                    :: nElec,nOcc,iOcc,i,j
    !
    TYPE(AtomSparseMat)        :: alphaBeta
    TYPE(MoleculeInfo),pointer :: molecule
    Real(realk)                :: ts,te,tsfull,tefull
    logical                    :: saveRecalcGab
    integer                    :: idmat,nmat,nrow,ncol
    TYPE(LSTENSOR),pointer     :: regCSfull,auxCSfull
    !
    Real(realk)                :: threshold,norm,dasum
    Logical,pointer            :: neighbours(:,:),MOcontrib(:,:)
    Real(realk)                :: minEigv,maxEigv,conditionNum
    Type(moleculeinfo),pointer :: atoms_A(:)
    Type(mat3d),pointer        :: calpha_ab_mo(:,:)
    Type(mat2d),pointer        :: alpha_beta_mo(:,:)
    Real(realk),pointer        :: MOcoeff(:,:),calpha_ab_block(:,:,:)
    Real(realk),pointer        :: MOmatK(:,:)
    Real(realk),pointer        :: matH_Q(:,:,:),matD_Q(:,:,:)
    Real(realk),pointer        :: tmp(:,:,:)
    Logical                    :: compute_coeff
    Character(80)              :: C_filename

    Integer                    :: INFO
    Integer,pointer            :: PIV(:)
    Real(realk)                :: chol_tol,sum
    Real(realk),pointer        :: matA(:,:),Work(:)

    IF (ndmat.GT.1) THEN
       WRITE(*,*)     &
            "The PARI approximation isn't implemented for unrestricted cases yet."
       WRITE(LUPRI,*) &
            "The PARI approximation isn't implemented for unrestricted cases yet."
       CALL LSQUIT('Error in II_get_pari_df_exchange_mat: ndmat>1 not tested',-1)
    ENDIF

    nrow = D(1)%nrow
    ncol = D(1)%ncol

    IF ((F(1)%nrow.NE.nrow).OR.(F(1)%ncol.NE.ncol)) &
         CALL LSQUIT('Error in II_get_pari_df_exchange_mat F/D',-1)

    IF (matrix_type .EQ. mtype_unres_dense)THEN
       IF (SETTING%SCHEME%NON_ROBUST_PARI) THEN 
          WRITE(*,*)     &
               "The PARI approximation isn't implemented for unrestricted cases yet."
          WRITE(LUPRI,*) &
               "The PARI approximation isn't implemented for unrestricted cases yet."
          CALL LSQUIT('Error in II_get_pari_df_exchange_mat. NR and unrestricted',-1)
       ENDIF
       IF(setting%IntegralTransformGC) THEN
          WRITE(*,*)     &
               "The PARI approximation isn't implemented for unrestricted cases yet."
          WRITE(LUPRI,*) &
               "The PARI approximation isn't implemented for unrestricted cases yet."
          CALL LSQUIT('Error in II_get_pari_df_exchange_mat. GC and unrestricted',-1)
       ENDIF
       IF (SETTING%SCHEME%FMM) THEN
          WRITE(*,*)     &
               "The PARI approximation isn't implemented for unrestricted cases yet."
          WRITE(LUPRI,*) &
               "The PARI approximation isn't implemented for unrestricted cases yet."
          call lsquit('Not allowed combination in II_get_pari_df_exchange_mat.'//&
               'FMM and unrestricted',-1)
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
    threshold=SETTING%SCHEME%intTHRESHOLD

    !Then get the full precalculated screening matrices G_ab and G_alpha
    !G_ab = integer part of log( max( sqrt (ab|ab) ))
    !where a and b runs on batches of atomic orbitals
    CALL II_getScreenMatFull(regCSfull,auxCSfull,Setting,lupri,luerr)
    !
    molecule => SETTING%MOLECULE(1)%p
    !
    call getMolecularDimensions(molecule,nAtoms,nBastReg,nBastAux)

    CALL LSTIMER('START ',te,ts,lupri)
    CALL LSTIMER('START ',tefull,tsfull,lupri)

    ! --- Build an array with the list of atoms with their respective 
    ! --- nb. of orbitals/auxiliary functions
    ! --- Build another array to know where each atom start in the complete matrices
    call setMolecularOrbitalInfo(molecule,orbitalInfo)
    
    !--- MO-based implementation from Manzer et al., JCTC, 2015, 11(2), pp 518-527 
    if (setting%scheme%MOPARI_K) then
       
       call mem_alloc(MOmatK,nBastReg,nBastReg)
       MOmatK=0E0_realk

       !Get the MO coeff by Cholesky Decomposition of the Density Matrix
       !The rank of the density matrix should be equal to the nb of Occ. orbitals
       !except for the iteration 0 where the density can be non valid
       chol_tol=1.0E-12_realk
       call mem_alloc(matA,nBastReg,nBastReg)
       call mem_alloc(PIV,nBastReg)
       call mem_alloc(work,2*nBastReg)
       INFO=-1
       nOcc=0
       PIV(:)=0
       
       do i=1,nBastReg
          do j=1,nBastReg
             matA(i,j) = Dfull(i,j,1)
          enddo
       enddo
           
       call dpstrf('L',nBastReg,matA,nBastReg,PIV,nOcc,chol_tol,Work,INFO)
       
       call mem_alloc(MOcoeff,nBastReg,nOcc)
       MOcoeff=0E0_realk
       do j=1,nOcc
          do i=j,nBastReg
             MOcoeff(PIV(i),j)=MatA(i,j)
          enddo
       enddo
       
       call mem_dealloc(PIV)
       call mem_dealloc(matA)
       call mem_dealloc(work)

       !write(lupri,*) 'Rank of the Cholesky Decomposition:',nOcc

       !write(lupri,'(/A/)') 'Atomic Contributions to the MO'
       !call mem_alloc(MOcontrib,nAtoms,nOcc)
       !MOcontrib=.false.
       !do iOcc=1,nOcc
       !   do iAtomA=1,nAtoms
       !      call getAtomicOrbitalInfo(orbitalInfo,iAtomA,nRegA,startRegA,endRegA,&
       !        nAuxA,startAuxA,endAuxA) 
       !      norm= dasum(nRegA,MOcoeff(startRegA:endRegA,iOcc),1)
       !      if (CEILING(LOG10(norm)).gt.-4) then
       !         MOContrib(iAtomA,iOcc)=.true.
       !      endif
       !   enddo
       !   write(lupri,*) MOcontrib(:,iOcc)
       !enddo
       !call mem_dealloc(MOcontrib)

       ! --- each molecule made of only one atom 
       allocate(atoms_A(nAtoms))  
 
       ! --- Initialize matrix of PARI coefficients
       allocate(calpha_ab_mo(nAtoms,nAtoms))
       do iAtomA=1,nAtoms
          do iAtomB=1,nAtoms
             nullify(calpha_ab_mo(iAtomA,iAtomB)%elements)
          enddo
       enddo

       ! --- Check if Pari-coeff have been computed and stored in previous 
       ! --- SCF iteration
       compute_coeff = .false.
       C_filename = 'CALPHA_AB'
       if (io_file_exist(C_filename,setting%IO)) then
          call io_read_mat3d_mo(calpha_ab_mo,nAtoms,nAtoms,C_filename,&
               setting%io,lupri,luerr)
       else 
          compute_coeff = .true.
       endif

       ! --- Initialize matrix of 2-center integrals
       allocate(alpha_beta_mo(nAtoms,nAtoms))
       do iAtomA=1,nAtoms
          do iAtomB=1,nAtoms
             nullify(alpha_beta_mo(iAtomA,iAtomB)%elements)
          enddo
       enddo

       call pari_set_atomic_fragments(molecule,atoms_A,nAtoms,lupri)

       ! --- Set up domains of atoms
       ! --- For B>=A neighbour(iAtomA,iAtomB= .true. if
       ! --- max(G_ab) > threshold with a in Reg(A) and b in Reg(B)
       call mem_alloc(neighbours,nAtoms,nAtoms)
       neighbours(:,:) = .false.
       call get_neighbours(neighbours,orbitalInfo,regCSfull,threshold,molecule,&
            atoms_A,lupri,luerr,setting)
      
       !write(lupri,'(/A/)') 'Matrix of atomic neighbours'
       !do iAtomA=1,nAtoms
       !   write(lupri,*) neighbours(iAtomA,:)
       !enddo

       ! Loop over A
       do iAtomA=1,nAtoms
          call getAtomicOrbitalInfo(orbitalInfo,iAtomA,nRegA,startRegA,endRegA,&
               nAuxA,startAuxA,endAuxA) 

          ! --- Construction of matrix H_i^nuQ for AtomQ=AtomA
          call mem_alloc(matH_Q,nBastReg,nAuxA,nOcc)
          matH_Q=0E0_realk
          call getHQcoeff(matH_Q,calpha_ab_mo,alpha_beta_mo,iAtomA,neighbours,&
               MOcoeff,orbitalInfo,setting,molecule,atoms_A,regCSfull,auxCSfull,&
               nOcc,lupri,luerr)
          
          ! --- Construction of matrix D_i^muQ for AtomQ=AtomA
          call mem_alloc(matD_Q,nAuxA,nOcc,nBastReg)
          matD_Q=0E0_realk
          call getDQcoeff(matD_Q,calpha_ab_mo,iAtomA,neighbours,MOcoeff,&
               orbitalInfo,setting,molecule,atoms_A,regCSfull,auxCSfull,&
               nOcc,lupri,luerr)
          
          ! --- Addition of the AtomQ=AtomA contribution to the matrix L^munu 
          call dgemm('N','N',nBastReg,nBastReg,nAuxA*nOcc,1.0E0_realk,matH_Q,&
               nBastReg,matD_Q,nAuxA*nOcc,1.0E0_realk,MOmatK,nBastReg)
          
          call mem_dealloc(matD_Q)
          call mem_dealloc(matH_Q)
       Enddo !Loop A
       
       ! --- If first iteration, store the PARI coeff on disk
       if (compute_coeff) then
          CALL io_add_filename(setting%io,C_filename,lupri)
          CALL io_write_mat3d_mo(calpha_ab_mo,nAtoms,nAtoms,C_filename,&
               setting%io,lupri,luerr) 
       endif

       ! --- Deallocations
       do iAtomA=1,nAtoms
          do iAtomB=1,nAtoms
             if (associated(calpha_ab_mo(iAtomA,iAtomB)%elements)) then
                call free_MAT3D(calpha_ab_mo(iAtomA,iAtomB))
             endif
          enddo
       enddo
       do iAtomA=1,nAtoms
          do iAtomB=1,nAtoms
             if (associated(alpha_beta_mo(iAtomA,iAtomB)%elements)) then
                call free_mat2d(alpha_beta_mo(iAtomA,iAtomB))
             endif
          enddo
       enddo
       call mem_dealloc(MOcoeff)
       call mem_dealloc(neighbours)
       call pari_free_atomic_fragments(atoms_A,nAtoms)
       deallocate(atoms_A) 
       deallocate(calpha_ab_mo)
       deallocate(alpha_beta_mo)
      
       ! --- Construct K_munu from L_munu
       do iRegA=1,nBastReg
          do iRegB=iRegA,nBastReg
             MOMatK(iRegA,iRegB) =  MOMatK(iRegA,iRegB) +  MOMatK(iRegB,iRegA)
             MOMatK(iRegB,iRegA) =  MOMatK(iRegA,iRegB)
          enddo
       enddo
       
       ! --- Add the exchange contribution (K) to the Fock matrix (F)
       call mat_set_from_full(MOMatK(:,:),&
            -Setting%Scheme%exchangeFactor,F(1),'exchange')
       
       ! --- Deallocation, resetting molecule to default
       call mem_dealloc(MOmatK)
       call typedef_setMolecules(setting,molecule,1,2,3,4)
       
    else !PARI_K

       call mem_alloc(Kfull_3cContrib,nBastReg,nBastReg,1,1,nmat)
       call mem_alloc(Kfull_2cContrib,nBastReg,nBastReg,1,1,nmat)
       call ls_DZERO(Kfull_3cContrib,nBastReg*nBastReg*nmat)
       call ls_DZERO(Kfull_2cContrib,nBastReg*nBastReg*nmat)

       ! --- Build the full huge matrices (alpha | beta)
       usemat = 0 
       call mem_alloc(alpha_beta,nBastAux,1,nBastAux,1,1)
       call ls_DZERO(alpha_beta,nBastAux*nBastAux)
       call ls_attach_gab_to_setting(setting,auxCSfull,auxCSfull)
       call initIntegralOutputDims(setting%output,nBastaux,1,nBastaux,1,1)
       call ls_getIntegrals(AODFdefault,AOempty,AODFdefault,AOempty,&
            CoulombOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
       call retrieve_Output(lupri,setting,alpha_beta,.FALSE.)
       call ls_free_gab_from_setting(setting,lupri)

       ! --- Memory allocation (all Density coefficients stored in an INefficient way)

       CALL LSTIMER('(al|be)',te,ts,lupri)

       ! --- Calculate the pair-atomic fitting coefficients 
       allocate(calpha_ab(nAtoms))
       call getPariCoefficients(LUPRI,LUERR,SETTING,calpha_ab,orbitalInfo,&
            regCSfull,auxCSfull)

       CALL LSTIMER('Coeffs ',te,ts,lupri)

       !re-set threshold 
       SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%K_THR

       DO idmat=1,nmat
          ! 
          ! Insert call to screening_matrices over the full molecule for both 
          ! AORdefault,AORdefault and AODFdefault
          ! interactions
          !
          ! --- Loop over atoms A
          DO iAtomA=1,nAtoms
             call getAtomicOrbitalInfo(orbitalInfo,iAtomA,nRegA,startRegA,endRegA,&
                  nAuxA,startAuxA,endAuxA)

             ! --- Construct density coefficients (dalpha_ad), alpha,a in A and d in Mol.
             call mem_alloc(dalpha_ad,nAuxA,nRegA,nBastReg)

             !> \todo (Patrick) limit the B atoms of the calpha_ab coefficients to 
             ! the vicinity of A
             call getDensityCoefficients(iAtomA,orbitalInfo,calpha_ab,&
                  Dfull(:,:,idmat:idmat),dalpha_ad)

             call lstimer('den-cof',te,ts,lupri)
             ! Three-center contributions
             Kfull => Kfull_3cContrib(:,:,:,:,idmat:idmat)
             DfullAO => Dfull(:,:,idmat:idmat)
             call pariK_ThreeCenter(Kfull,DfullAO,setting,calpha_ab,dalpha_ad,&
                  iAtomA,orbitalInfo,nBastReg,nAtoms,regCSfull,auxCSfull,lupri,luerr)
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
          Kfull_3cContrib(:,:,1,1,idmat) = Kfull_3cContrib(:,:,1,1,idmat) &
               + Kfull_2cContrib(:,:,1,1,idmat)
          
          call symmetrizeKfull_3cContrib(Kfull_3cContrib(:,:,1,1,idmat),nBastReg)
          !write(lupri,'(/A/)') 'Final MOmatK'
          !do iRegA=1,nBastReg
          !   write(lupri,*) Kfull_3cContrib(iRegA,:,1,1,idmat)
          !enddo
       ENDDO

       call freePariCoefficients(calpha_ab,orbitalInfo%nAtoms)
       deallocate(calpha_ab)
       call mem_dealloc(alpha_beta)
       
       ! --- ADD THE EXCHANGE CONTRIBUTION (K) TO THE FOCK MATRIX (F)
       
       IF(matrix_type .EQ. mtype_unres_dense)THEN
          DO idmat=1,ndmat
             call DAXPY(F(1)%nrow*F(1)%ncol,-Setting%Scheme%exchangeFactor,&
                  Kfull_3cContrib(:,:,1,1,2*idmat-1),1,F(idmat)%elms,1)
             call DAXPY(F(1)%nrow*F(1)%ncol,-Setting%Scheme%exchangeFactor,&
                  Kfull_3cContrib(:,:,1,1,2*idmat  ),1,F(idmat)%elmsb,1)
          ENDDO
       ELSE !CLOSED_SHELL
          DO idmat=1,ndmat
             call mat_set_from_full(Kfull_3cContrib(:,:,1,1,idmat),&
                  -Setting%Scheme%exchangeFactor,F(idmat),'exchange')
          ENDDO
       ENDIF
       
       ! --- Free temporary allocd memory
       !call ls_freeDfull(setting)
       call mem_dealloc(Kfull_3cContrib)
       call mem_dealloc(Kfull_2cContrib)
    endif
    
    CALL freeMolecularOrbitalInfo(orbitalInfo)
    call mem_dealloc(Dfull)

    !call free_AtomSparseMat(alphaBeta)
    CALL LSTIMER('PARI-K',tefull,tsfull,lupri)
    !call lsquit('End of PARI-k',-1)
  END SUBROUTINE II_get_pari_df_exchange_mat

  !> \brief Set up domains of atoms such that max(G_ab) > threshold
  !> \brief with a in Reg(A) and b in Reg(B)
  !> \author E. Rebolini
  !> \date 2015-02
  !> \param neighbours (iAtomA,1:M) contains the indices of the M neighbours B 
  !> \param orbitalInfo Orbital information
  !> \param regCSfull Screening matrix G_ab = log( sqrt( (ab|ab) ) )
  !> \param threshold Minimum value for G_ab
  !> \param molecule Description of the molecule
  !> \param atoms Each molecule made of only one atom 
  !> \param lupri Default print-unit for output
  !> \param luerr Default print-unit for termination
  !> \param setting Contains information about the integral settings
subroutine get_neighbours(neighbours,orbitalInfo,regCSfull,threshold,molecule,&
     atoms,lupri,luerr,setting)
  implicit none
  Logical,pointer                       :: neighbours(:,:)
  Type(molecularorbitalinfo),intent(in) :: orbitalInfo
  Type(lstensor),pointer                :: regCSfull
  Real(realk)                           :: threshold
  Type(moleculeinfo),pointer            :: molecule
  Type(moleculeinfo),pointer            :: atoms(:)
  Integer,intent(in)                    :: lupri,luerr
  Type(lssetting),intent(inout)         :: setting
    
  Integer                               :: nAtoms
  Integer                               :: iAtomA,nRegA,startRegA,endRegA
  Integer                               :: nAuxA,startAuxA,endAuxA
  Integer                               :: iAtomB,nRegB,startRegB,endRegB
  Integer                               :: nAuxB,startAuxB,endAuxB
  TYPE(moleculeinfo),pointer            :: AB
  TYPE(moleculeinfo),target             :: ABtarget
  Integer                               :: nAux
  TYPE(LSTENSOR),pointer                :: auxCSab,regCSab
  INTEGER                               :: atomsAB(2),dummyAtoms(1)
  Integer(kind=short)                   :: maxgab
 
  nAtoms=orbitalInfo%nAtoms
  !write(lupri,*) 'threshold',threshold
  do iAtomA=1,nAtoms
       call getAtomicOrbitalInfo(orbitalInfo,iAtomA,nRegA,startRegA,endRegA,&
            nAuxA,startAuxA,endAuxA)
       !Case A=B
       neighbours(iAtomA,iAtomA)=.true.
       !Case A<>B
       do iAtomB=iAtomA+1,nAtoms
          call getAtomicOrbitalInfo(orbitalInfo,iAtomB,nRegB,startRegB,endRegB,&
               nAuxB,startAuxB,endAuxB)          
          nAux=nAuxA+nAuxB
          CALL pariSetPairFragment(AB,ABtarget,setting%basis(1)%p,molecule,atoms,&
               molecule%nAtoms,iAtomA,iAtomB,nAuxA,nAuxB,nAux,lupri)
          
          NULLIFY(regCSab)
          ALLOCATE(regCSab)
          call ls_subScreenAtomic(regCSab,regCSfull,iAtomA,iAtomB,&
               nRegA,nRegB,.FALSE.)
          maxgab=regCSab%maxgabelm
          !write(lupri,*) iAtomA,iAtomB,maxgab
          if (2*maxgab.GE.(CEILING(LOG10(threshold)))) then
             neighbours(iAtomA,iAtomB)=.true.
          endif
          call lstensor_free(regCSab)
          DEALLOCATE(regCSab)
          CALL pariFreePairFragment(AB,iAtomA,iAtomB)
       enddo      
    enddo
    
end subroutine

subroutine symmetrizeKfull_3cContrib(Kfull_3cContrib,nBastReg)
  implicit none
  integer,intent(in) :: nBastReg
  real(realk),intent(inout) :: Kfull_3cContrib(nBastReg,nBastReg)
  
  Kfull_3cContrib = Kfull_3cContrib + transpose(Kfull_3cContrib)
  
end subroutine symmetrizeKfull_3cContrib

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
      !build Matrix 
      call mat_zero(alphabeta)
      call initIntegralOutputDims(setting%output,naux,1,naux,1,1)
      call ls_getIntegrals(AODFdefault,AOempty,&
           &AODFdefault,AOempty,CoulombOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
      call retrieve_Output(lupri,setting,alphabeta,.FALSE.)
      !Make Choleksy-factorization
      CALL mat_dpotrf(alphabeta)
      !Save Cholesky-factors to file
      call io_add_filename(SETTING%IO,Filename,LUPRI)
      call io_write_mat(alphabeta,Filename,SETTING%IO,LUPRI,LUERR)
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
        call Test_if_64bit_integer_required(nbastaux,nbastaux)
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
     call Test_if_64bit_integer_required(nBastAux,nBast,nBast)
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
           call Test_if_64bit_integer_required(nBastAux,nbast,nBastLocA)
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

SUBROUTINE II_get_RI_AlphaCD_3CenterInt(LUPRI,LUERR,FullAlphaCD,SETTING,&
     & nbasisAux,nbasis)
  IMPLICIT NONE
  Integer,intent(in)            :: LUPRI,LUERR,nbasis,nbasisAux
  REAL(REALK),intent(inout)     :: FullAlphaCD(nbasisAux,nbasis,nbasis)
  TYPE(LSSETTING),intent(inout) :: SETTING
  !
  Integer                    :: nAtoms,nBastAux,nBast,N,K,M,ialpha,v,a
  Integer                    :: BDIAG,IDIAG,ILOC,JLOC,ALPHAAUX,GAMMA,DELTA
  Real(realk)                :: TSTART,TEND,tmp
  logical :: MasterWakeSlaves

  call LSTIMER('START ',TSTART,TEND,LUPRI)  
  call getMolecularDimensions(SETTING%MOLECULE(1)%p,nAtoms,nBast,nBastAux)
  IF(nbasis.NE.nBast)THEN
     print*,'nbasis',nbasis
     print*,'nbast',nbast
     Call lsquit('dim mismatch in II_get_df_exchange2',-1)
  ENDIF
  IF(nbasisAux.NE.nBastAux)THEN
     print*,'nbasisAux',nbasisAux
     print*,'nbastAux',nbastAux
     Call lsquit('dim mismatch in II_get_df_exchange2',-1)
  ENDIF
  !set threshold 
  SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%K_THR
  !(alpha|cd)
  !call typedef_setMolecules(setting,molecule,1,2,3,4)
  call initIntegralOutputDims(setting%Output,nBastAux,1,nbast,nbast,1)
  !MPI is used inside this routine 
  !but since this is called by both master and slaves
  !the master should not wake up the slaves
  MasterWakeSlaves = SETTING%SCHEME%MasterWakeSlaves
  SETTING%SCHEME%MasterWakeSlaves = .FALSE.
  call ls_getIntegrals(AODFdefault,AOempty,AORdefault,AORdefault,&
       & CoulombOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
  SETTING%SCHEME%MasterWakeSlaves = MasterWakeSlaves
  CALL retrieve_Output(lupri,setting,FullAlphaCD,.FALSE.)  
  call LSTIMER('AlphaCD',TSTART,TEND,LUPRI)
END SUBROUTINE II_get_RI_AlphaCD_3CenterInt

SUBROUTINE II_get_RI_AlphaCD_3CenterInt2(LUPRI,LUERR,FullAlphaCD,SETTING,nbasisAux,&
     & nbasis,nvirt,nocc,Cvirt,Cocc,maxsize,mynum,numnodes,InOper)
  IMPLICIT NONE
  Integer,intent(in)            :: LUPRI,LUERR,nbasis,nbasisAux
  Integer,intent(in)            :: nocc,nvirt,mynum,numnodes  
  REAL(REALK),pointer           :: FullAlphaCD(:,:,:)!(nbasisAuxMPI,nvirt,nocc)
  REAL(REALK),intent(in)        :: Cocc(nbasis,nocc)
  REAL(REALK),intent(in)        :: Cvirt(nbasis,nvirt)
  TYPE(LSSETTING),intent(inout) :: SETTING
  integer(kind=long),intent(in) :: maxsize
  integer,optional              :: InOper 
  !
  integer(kind=long)         :: nsize
  Integer                    :: nAtoms,nAtomsAux,nBastAux,nBast,N,K,M,ialpha,v,a
  Integer                    :: BDIAG,IDIAG,ILOC,JLOC,ALPHAAUX,GAMMA,DELTA
  Integer                    :: ALPHA,BETA,I,Oper 
  Real(realk) :: TSTART,TEND,tmp,TMP1
  real(realk),pointer :: AlphaCD(:,:,:),AlphaCD2(:,:,:)
  TYPE(MoleculeInfo),pointer      :: molecule1,molecule2,molecule3,molecule4
  TYPE(MoleculeInfo),pointer      :: ATOMS(:)
  TYPE(MOLECULARORBITALINFO) :: orbitalInfo
  Integer :: iAtomA,nAuxA,B
  integer :: J,mynum2,startF,iatomampi,MynbasisAuxMPI,nBastAuxT
  logical :: doMPI,MasterWakeSlaves
  integer,pointer :: nbasisAuxMPI(:),startAuxMPI(:,:),AtomsMPI(:,:),nAtomsMPI(:),nAuxMPI(:,:)
  !
  integer :: iShell, nAuxShellA,nbatches
  integer, pointer :: batchdim(:)

  IF(present(InOper))THEN
     Oper = InOper
  ELSE
     Oper = CoulombOperator
  ENDIF

  SETTING%scheme%CS_SCREEN = .FALSE.
  SETTING%scheme%PS_SCREEN = .FALSE.
  SETTING%scheme%OD_SCREEN = .FALSE.
  SETTING%scheme%OE_SCREEN = .FALSE.

  nsize = (nbasis*nbasis*nbasisAux+nocc*nbasis*nbasisAux)*mem_realsize
  !in case of setting%molBuild = .TRUE. then the setting%molecule(1)%p can 
  !not be used as pointers
  molecule1 => SETTING%MOLECULE(1)%p
  molecule2 => SETTING%MOLECULE(2)%p
  molecule3 => SETTING%MOLECULE(3)%p
  molecule4 => SETTING%MOLECULE(4)%p
!  print*,'nsize',nsize
!  print*,'maxsize',maxsize
  call LSTIMER('START ',TSTART,TEND,LUPRI)
  IF((numnodes.EQ.1.AND.maxsize.GT.nsize).AND.&
       & (.NOT.setting%scheme%ForceRIMP2memReduced))THEN     
     !serial version
     call getMolecularDimensions(SETTING%MOLECULE(1)%p,nAtomsAux,nBast,nBastAux) !Aux
     call getMolecularDimensions(SETTING%MOLECULE(3)%p,nAtoms,nBast,nBastAuxT)!Reg 
     IF(nbasis.NE.nBast)THEN
        print*,'nbasis',nbasis
        print*,'nbast',nbast
        Call lsquit('dim mismatch in II_get_RI_AlphaCD_3CenterInt2',-1)
     ENDIF
     IF(nbasisAux.NE.nBastAux)THEN
        print*,'nbasisAux',nbasisAux
        print*,'nbastAux',nbastAux
        Call lsquit('dim mismatch in II_get_RI_AlphaCD_3CenterInt2',-1)
     ENDIF
     IF(nbasisAUX.NE.nBastAux)Call lsquit('dim mismatch in ',-1)
    !set threshold 
     SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%K_THR
     !(alpha|cd)
     !call typedef_setMolecules(setting,molecule,1,2,3,4)
     call initIntegralOutputDims(setting%Output,nBastAux,1,nbast,nbast,1)

     !both Master and non-master calls this
     call ls_getIntegrals(AODFdefault,AOempty,AORdefault,AORdefault,&
          & Oper,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
     
     call mem_alloc(AlphaCD,nbastAux,nbast,nbast)
     CALL retrieve_Output(lupri,setting,AlphaCD,.FALSE.)
     
!     print*,'II_get_RI_AlphaCD_3CenterInt2: FullAlphaCD(1:nA,1:4) AO'
!     call ls_output(AlphaCD,1,nbastAux,1,4,nbastAux,nbast*nbast,1,6)

     ! Transform index delta to diagonal occupied index 
     !(alphaAux;gamma,J) = (alphaAux;gamma,delta)*C(delta,J)
     M = nbastAux*nbast !rows of Output Matrix
     N = nocc           !columns of Output Matrix
     K = nbast          !summation dimension
     call mem_alloc(AlphaCD2,nbastAux,nbast,nocc)
     call Test_if_64bit_integer_required(N,K)
     call Test_if_64bit_integer_required(M,K)
     call Test_if_64bit_integer_required(M,N)
     call dgemm('N','N',M,N,K,1.0E0_realk,AlphaCD,M,Cocc,nbast,&
          & 0.0E0_realk,AlphaCD2,M)

!     print*,'II_get_RI_AlphaCD_3CenterInt2: FullAlphaCD(1:nA,1:4) HALF'
!     call ls_output(AlphaCD2,1,nbastAux,1,4,nbastAux,nbast*nocc,1,6)

     call mem_dealloc(AlphaCD)
     call mem_alloc(FullAlphaCD,nbasisAux,nvirt,nocc)
     !(alphaAux,B,J) = (alphaAux,gamma,delta)*C(gamma,B)
!$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) &
!$OMP PRIVATE(BDIAG,IDIAG,ALPHAAUX,GAMMA,TMP) &
!$OMP SHARED(Cvirt,AlphaCD2,FullAlphaCD,nvirt,nocc,nbasisAUX,nbasis)
     do BDIAG = 1,nvirt
        do IDIAG = 1,nocc
           do ALPHAAUX = 1,nbasisAUX
              TMP = 0.0E0_realk
              do GAMMA = 1,nbasis
                 TMP = TMP + Cvirt(GAMMA,BDIAG)*AlphaCD2(ALPHAAUX,GAMMA,IDIAG)
              enddo
              FullAlphaCD(ALPHAAUX,BDIAG,IDIAG) = TMP
           enddo
        enddo
     enddo
!$OMP END PARALLEL DO
     call mem_dealloc(AlphaCD2)
!     print*,'II_get_RI_AlphaCD_3CenterInt2: FullAlphaCD(1:nA,1:4) MO'
!     call ls_output(FullAlphaCD,1,nbastAux,1,4,nbastAux,nvirt*nocc,1,6)

  ELSE
     !need to deactivate MPI inside ls_getIntegrals. 
     !both master and slaves call this routine 
     doMPI = SETTING%SCHEME%doMPI
     MasterWakeSlaves = SETTING%SCHEME%MasterWakeSlaves
     SETTING%SCHEME%doMPI = .FALSE.
     SETTING%SCHEME%MasterWakeSlaves = .FALSE.

     !MPI version and memory reduced version
     !     print*,'MEMORY OPTIMIZED VERSION'
     !     WRITE(lupri,*)'MEMORY OPTIMIZED VERSION'
     !MEMORY OPTIMIZED VERSION
     call getMolecularDimensions(molecule1,nAtomsAux,nBast,nBastAux) !Aux
     call getMolecularDimensions(molecule3,nAtoms,nBast,nBastAuxT)!Reg 
     allocate(ATOMS(nAtomsAux))
     CALL pari_set_atomic_fragments(molecule1,ATOMS,nAtomsAux,lupri)

     ! Split only on of the atomic loops over nodes
     call mem_alloc(nbasisAuxMPI,numnodes)
     call mem_alloc(nAtomsMPI,numnodes)    
     call mem_alloc(startAuxMPI,nAtomsAux,numnodes)
     call mem_alloc(AtomsMPI,nAtomsAux,numnodes)
     call mem_alloc(nAuxMPI,nAtomsAux,numnodes)
     call getRIbasisMPI(molecule1,nAtomsAux,numnodes,nbasisAuxMPI,startAuxMPI,&
          & AtomsMPI,nAtomsMPI,nAuxMPI)
     call mem_dealloc(startAuxMPI)
     MynbasisAuxMPI = nbasisAuxMPI(mynum+1)
     call mem_dealloc(nbasisAuxMPI)
     IF(MynbasisAuxMPI.GT.0)THEN
        call mem_alloc(FullAlphaCD,MynbasisAuxMPI,nvirt,nocc)
        FullalphaCD = 0.0E0_realk
        startF = 0
        DO iAtomAMPI=1,nAtomsMPI(mynum+1)
           iAtomA = AtomsMPI(iAtomAMPI,mynum+1)
           nAuxA = nAuxMPI(iAtomAMPI,mynum+1)
           call typedef_setMolecules(setting,molecule3,3,4,ATOMS(iAtomA),1)             
           nsize = (MynbasisAuxMPI*nvirt*nocc+2*nAuxA*nBast*nBast)*mem_realsize
           IF(nsize.GT.maxsize.OR.setting%scheme%ForceRIMP2memReduced)THEN     
              !Memory reduced version:   We split up the atom into batches
              nullify(batchdim)
              call build_minimalbatchesOfAOs2(lupri,setting,batchdim,nbatches,'D')
              BatchD: do iShell = 1,nbatches
                 nAuxShellA = batchdim(iShell)
                 call initIntegralOutputDims(setting%Output,nAuxShellA,1,nBast,nBast,1)
                 setting%batchindex(1)=iShell
                 setting%batchdim(1)=nAuxShellA
                 call ls_getIntegrals(AODFdefault,AOEmpty,AORdefault,AORdefault,&
                      & Oper,RegularSpec,Contractedinttype,SETTING,LUPRI,LUERR)
                 setting%batchindex(1)=0
                 setting%batchdim(1)=nAuxShellA
                 call mem_alloc(alphaCD,nAuxShellA,nBast,nBast)
                 CALL retrieve_Output(lupri,setting,alphaCD,.FALSE.)
                 ! Transform index delta to diagonal occupied index 
                 !(alphaAux;gamma,J) = (alphaAux;gamma,delta)*C(delta,J)
                 M = nAuxShellA*nbast   !rows of Output Matrix
                 N = nocc               !columns of Output Matrix
                 K = nbast              !summation dimension
                 call mem_alloc(AlphaCD2,nAuxShellA,nBast,nocc)
                 call Test_if_64bit_integer_required(N,K)
                 call Test_if_64bit_integer_required(M,K)
                 call Test_if_64bit_integer_required(M,N)
                 call dgemm('N','N',M,N,K,1.0E0_realk,AlphaCD,M,Cocc,nbast,0.0E0_realk,AlphaCD2,M)
                 call mem_dealloc(AlphaCD)
                 !(alphaAux,B,J) = (alphaAux,gamma,delta)*C(gamma,B)
                 call DF3centerTrans(nvirt,nocc,nbast,nAuxShellA,MynbasisAuxMPI,Cvirt,AlphaCD2,FullAlphaCD,startF)
                 call mem_dealloc(AlphaCD2)
                 startF = startF + nAuxShellA !nAUXA
              ENDDO BatchD
              call mem_dealloc(batchdim)
           ELSE
              !Original code
              call initIntegralOutputDims(setting%Output,nAuxA,1,nBast,nBast,1)
              call ls_getIntegrals(AODFdefault,AOEmpty,AORdefault,AORdefault,&
                   & Oper,RegularSpec,Contractedinttype,SETTING,LUPRI,LUERR)
              call mem_alloc(alphaCD,nAuxA,nBast,nBast)
              CALL retrieve_Output(lupri,setting,alphaCD,.FALSE.)
              ! Transform index delta to diagonal occupied index 
              !(alphaAux;gamma,J) = (alphaAux;gamma,delta)*C(delta,J)
              M = nAuxA*nbast        !rows of Output Matrix
              N = nocc               !columns of Output Matrix
              K = nbast              !summation dimension
              call mem_alloc(AlphaCD2,nAuxA,nBast,nocc)
              call Test_if_64bit_integer_required(N,K)
              call Test_if_64bit_integer_required(M,K)
              call Test_if_64bit_integer_required(M,N)
              call dgemm('N','N',M,N,K,1.0E0_realk,AlphaCD,M,Cocc,nbast,0.0E0_realk,AlphaCD2,M)
              call mem_dealloc(AlphaCD)
              !(alphaAux,B,J) = (alphaAux,gamma,delta)*C(gamma,B)
              call DF3centerTrans(nvirt,nocc,nbast,nAuxA,MynbasisAuxMPI,Cvirt,AlphaCD2,FullAlphaCD,startF)
              call mem_dealloc(AlphaCD2)
              startF = startF + nAUXA
           ENDIF
        ENDDO
        !restore 
        call typedef_setMolecules(setting,molecule1,1,molecule2,2,&
             & molecule3,3,molecule4,4)
        call pari_free_atomic_fragments(ATOMS,nAtomsAux)        
     ENDIF
     deallocate(ATOMS)
     call mem_dealloc(AtomsMPI)
     call mem_dealloc(nAtomsMPI)
     call mem_dealloc(nAuxMPI)
     SETTING%SCHEME%doMPI = doMPI
     SETTING%SCHEME%MasterWakeSlaves = MasterWakeSlaves
  ENDIF
  !restore 
  SETTING%MOLECULE(1)%p => molecule1
  SETTING%MOLECULE(2)%p => molecule2
  SETTING%MOLECULE(3)%p => molecule3
  SETTING%MOLECULE(4)%p => molecule4
  call ls_setDefaultFragments(setting)
  call LSTIMER('AlphaCD',TSTART,TEND,LUPRI)
  SETTING%scheme%CS_SCREEN = .TRUE.
  SETTING%scheme%PS_SCREEN = .TRUE.
  SETTING%scheme%OD_SCREEN = .TRUE.
  SETTING%scheme%OE_SCREEN = .TRUE.
END SUBROUTINE II_get_RI_AlphaCD_3CenterInt2

!> \brief Calculates the 3 center 2 electron repulsion integral assuming all nodes have full Integral
!> This includes a allreduce! This also assumes that all slaves are alive and all call this subroutine        
!> WARNING: in case of MPI the FullAlphaCD(:,:,:) will only contain partial contributions calculated
!>          by this rank so you may need a lsmpi_allreduce or something, depending on your needs.
!> \author Thomas Kjaergaard
!> \date 2015
SUBROUTINE II_get_RI_AlphaCD_3CenterIntFullOnAllNN(LUPRI,LUERR,FullAlphaCD,&
     & SETTING,nAux,n1,n2,intspec,MaxnAux,nMO1,nMO2,AOtoMO,C1,C2,nthreads,dim1,&
     & GindexToLocal,DECPRINTLEVEL,use_bg_bufInput)
  IMPLICIT NONE
  Integer,intent(in)     :: LUPRI,LUERR,n1,n2,nAux,MaxnAux,nthreads
  integer,intent(in)     :: nMO1,nMO2,dim1,DECPRINTLEVEL
  REAL(REALK),pointer    :: FullAlphaCD(:) !dim1,nMO1,nMO2
  TYPE(LSSETTING),intent(inout) :: SETTING
  character,intent(in)  :: intspec(4)
  logical :: AOtoMO
  real(realk) :: C1(n1,nMO1),C2(n2,nMO2)  
  integer :: GindexToLocal(nAux)
  logical,optional :: use_bg_bufInput
  !
  real(realk) :: TSTART,TEND
  real(realk),pointer :: TmpAlphaCD(:),w0(:)
  logical :: MasterWakeSlaves,RestricedSize,use_bg_buf
  integer :: Oper,MaxN,i,startF,N,M,K,TID
  integer :: ao(3),dims(3),GlobalToLocal(nAux)
  integer(kind=8) :: n8,w0size,w1size
  type(C_PTR) :: cpointer
  use_bg_buf = .FALSE.
  IF(present(use_bg_bufInput)) use_bg_buf = use_bg_bufInput
  call nullThermiteIntTransform()
  call InitThermiteIntTransform1(GindexToLocal,nAux) 

  MaxN = MaxnAux
  RestricedSize = MaxN.LT.dim1
  IF(MaxnAux.GT.dim1) MaxN = dim1
  
  IF(.NOT.AOtoMO)THEN
     IF(MaxN.NE.nAux)CALL LSQUIT('dim mismatch in II_get_RI_AlphaCD_3CenterIntFullOnAllNN',-1)
  ENDIF
  call LSTIMER('START ',TSTART,TEND,LUPRI)  
  call GetOperatorFromCharacter(Oper,intspec(4),Setting)
  call GetAOtypes(ao,3,intspec(1:3))
  SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%K_THR
!  WRITE(LUPRI,*)'II_get_RI_AlphaCD_3CenterIntFullOnAll Integral Threshold:',SETTING%SCHEME%intTHRESHOLD
  !(alpha|cd)  
  call nullifyIntegralOutput(setting%Output)
  nullify(setting%output%resulttensor)
  IF(RestricedSize)THEN
     !Build AO integrals in setting%output%ResultMat
     !When AO integrals becomes bigger than (MaxN,n1,n2)
     !the AO integrals are transformed to
     !MO basis with accumulation in setting%output%Result3D

     w1size = dim1*nMO1*nMO2
     IF(DECPRINTLEVEL.GT.2)WRITE(LUPRI,'(A,F13.5,A)')&
          & '3 center RI: Allocating MO (alpha|AI)',dim1*nMO1*nMO2*8E-9_realk,' GB'
     IF(use_bg_buf)THEN
        call mem_pseudo_alloc(FullAlphaCD, w1size)
     ELSE
        call mem_alloc(FullAlphaCD,w1size)
     ENDIF

     call InitThermiteIntTransform(n1,nMO1,n2,nMO2,C1,C2,dim1,MaxN,nthreads) !full
     call initIntegralOutputDims(setting%Output,MaxN,1,n1,n2,nthreads)
     IF(DECPRINTLEVEL.GT.2)WRITE(LUPRI,'(A,F13.5,A,F13.5,A)')&
          & '3 center RI: Allocating 5dim Buffer',MaxN*n1*n2*8E-9_realk,&
          & ' GB for each thread = ',MaxN*n1*n2*nthreads*8E-9_realk,' GB'
     IF(use_bg_buf)THEN
        w0size = MaxN*n1*n2*nthreads
        call mem_pseudo_alloc(w0, w0size)
!        setting%output%ResultMat(1:MaxN,1,1:n1,1:n2,1:nthreads) => w0(1:w0size)
        cpointer = c_loc(w0(1))
        call c_f_pointer(cpointer,setting%output%ResultMat,[MaxN,1,n1,n2,nthreads])
     ELSE
        call mem_alloc(setting%output%ResultMat,MaxN,1,n1,n2,nthreads)
     ENDIF
     n8 = MaxN*n1*n2*nthreads
     call ls_dzero8(setting%output%ResultMat,n8) !due to screening 
     IF(DECPRINTLEVEL.GT.2)WRITE(LUPRI,'(A,F13.5,A)')&
          & '3 center RI: Allocating TmpArray of ',MaxN*n1*nMO2*8E-9_realk,' GB'
     call ThermiteIntTransform_alloc_TmpArray(use_bg_buf) !MaxN,n1,nMO2
     call ls_dzero8(FullAlphaCD,w1size)
!     setting%output%Result3D => FullAlphaCD
     cpointer = c_loc(FullAlphaCD(1))
     call c_f_pointer(cpointer,setting%output%Result3D,[dim1,nMO1,nMO2])
     setting%Output%ndim3D(1) = dim1
     setting%Output%ndim3D(2) = nMO1
     setting%Output%ndim3D(3) = nMO2
     setting%Output%ndim3D(4) = nAux
     !memory requirements:
     !nAux*nMO1*nMO2+MaxN*(n1*nMO2+n1*n2)+n1*nMO1+n2*nMO2 
     !choose
     !MaxN = FLOOR((MaxSize/8.0E-9_realk-nAux*nMO1*nMO2-n1*nMO1-n2*nMO2)/(n1*nMO2+n1*n2))
  ELSE 
     !Full AO integral Fits so the AO integrals are build 
     !directly in setting%output%Result3D (AO to MO happens outside Integral code)
     DoThermiteIntTransform = .FALSE.
     call initIntegralOutputDims(setting%Output,1,1,1,1,1)
     call mem_alloc(setting%output%ResultMat,1,1,1,1,1)
     n8 = dim1*n1*n2
     IF(DECPRINTLEVEL.GT.2)WRITE(LUPRI,'(A,F13.5,A)')&
          & '3 center RI: Allocating AO (alpha|CD)',dim1*n1*n2*8E-9_realk,' GB'
     w1size = dim1*n1*n2
     IF(use_bg_buf)THEN
        call mem_pseudo_alloc(FullAlphaCD, w1size)
     ELSE
        call mem_alloc(FullAlphaCD,w1size)
     ENDIF
     call ls_dzero8(FullAlphaCD,w1size)
     cpointer = c_loc(FullAlphaCD(1))
     call c_f_pointer(cpointer,setting%output%Result3D,[dim1,n1,n2])
!     setting%output%Result3D => FullAlphaCD
     setting%Output%ndim3D(1) = dim1
     setting%Output%ndim3D(2) = n1
     setting%Output%ndim3D(3) = n2
     setting%Output%ndim3D(4) = nAux
     !memory requirements:
     !MAX(nAux*nMO1*nMO2+nAux*n1*nMO2,nAux*n1*n2+nAux*n1*nMO2)
  ENDIF

  !MPI is used inside this routine but since this is called by both master and slaves
  !the master should not wake up the slaves
  MasterWakeSlaves = SETTING%SCHEME%MasterWakeSlaves
  SETTING%SCHEME%MasterWakeSlaves = .FALSE.
  setting%output%FullAlphaCD = .TRUE.  
  call ls_getIntegrals1(ao(1),AOempty,ao(2),ao(3),Oper,RegularSpec,ContractedInttype,0,SETTING,LUPRI,LUERR)
  setting%output%FullAlphaCD = .FALSE.
  SETTING%SCHEME%MasterWakeSlaves = MasterWakeSlaves
  call LSTIMER('AlphaCD',TSTART,TEND,LUPRI)

  IF(AOtoMO)THEN
     IF(DoThermiteIntTransform)THEN
        DO TID=1,nthreads
           IF(iLocalTIT2(TID).NE.0)THEN
              call InitThermiteIntThreadID(TID-1)
              !transform last contributions to MO
              call ThermiteIntTransform_AOtoMOInternalFinal(FullAlphaCD,&
                   & dim1,nMO1,nMO2,setting%output%ResultMat,MaxN,n1,n2)
           ENDIF
        ENDDO
     ELSE
!        print*,'II_get_RI_AlphaCD_3CenterIntFullOnAll:    FullAlphaCD(1:nA,1:4) AO'
!        call ls_output(FullAlphaCD,1,dim1,1,4,dim1,n1*n2,1,6)

        !Perform AO to MO
        M = dim1*n1   !rows of Output Matrix
        N = nMO2      !columns of Output Matrix
        K = n2        !summation dimension
        IF(DECPRINTLEVEL.GT.2)WRITE(LUPRI,'(A,F13.5,A)')&
             & '3 center RI: Allocating Full TmpArray ',M*N*8E-9_realk,' GB'
        n8 = M*N
        IF(use_bg_buf)THEN
           call mem_pseudo_alloc(TmpAlphaCD,n8)
        ELSE
           call mem_alloc(TmpAlphaCD,n8)
        ENDIF
        call dgemm('N','N',M,N,K,1.0E0_realk,FullAlphaCD,M,C2,K,0.0E0_realk,TmpAlphaCD,M)

!        print*,'II_get_RI_AlphaCD_3CenterIntFullOnAll:    FullAlphaCD(1:nA,1:4) HALF'
!        call ls_output(TmpAlphaCD,1,dim1,1,4,dim1,n1*nMO2,1,6)
        w1size = dim1*nMO1*nMO2
        IF(use_bg_buf)THEN
           !do nothing
        ELSE
           call mem_dealloc(FullAlphaCD)
           call mem_alloc(FullAlphaCD,w1size)
        ENDIF

        call DF3centerTrans4(nMO1,nMO2,n1,dim1,C1,TmpAlphaCD,FullAlphaCD(1:w1size))

        IF(use_bg_buf)THEN
           call mem_pseudo_dealloc(TmpAlphaCD)
           !This looks weird but FullAlphaCD is currently pointing to the first 1:dim1*n1*n2 elements
           !of a "permanent" memory array. I only need the first 1:dim1*nMO1*nMO2 elements
           !so I shrink the array dimension by deassociating (NOT deallocating) and reassociate
           call mem_pseudo_dealloc(FullAlphaCD)
           call mem_pseudo_alloc(FullAlphaCD, w1size)
        ELSE
           call mem_dealloc(TmpAlphaCD)
        ENDIF
     ENDIF
!     print*,'II_get_RI_AlphaCD_3CenterIntFullOnAll: FullAlphaCD(1:nA,1:4) MO dim1=',dim1
!     call ls_output(FullAlphaCD,1,dim1,1,4,dim1,n1*n2,1,6)
  ENDIF
  call FreeThermiteIntTransform(use_bg_buf)
  IF(RestricedSize)THEN
     IF(use_bg_buf)THEN
        call mem_pseudo_dealloc(w0)
        nullify(setting%output%ResultMat)
     ELSE
        call mem_dealloc(setting%output%ResultMat)
     ENDIF
  ELSE
     call mem_dealloc(setting%output%ResultMat)
  ENDIF
END SUBROUTINE II_get_RI_AlphaCD_3CenterIntFullOnAllNN

subroutine II_get_RI_AlphaCD_3CenterIntFullOnAllNNdim(setting,nAuxMPI,&
     & IndexToGlobal,numnodes,MaxnAuxMPI,AODF,GindexToLocal,nbasisAux)
implicit none
TYPE(LSSETTING),intent(inout) :: SETTING
integer,intent(in) :: numnodes,nbasisAux
integer,intent(inout) :: nAuxMPI(numnodes),MaxnAuxMPI
integer :: GindexToLocal(nbasisAux)
integer,pointer :: IndexToGlobal(:,:) !intent(inout)
character,intent(in)  :: AODF
call MPIdistributeAOs(setting,AODF,nAuxMPI,numnodes,IndexToGlobal,MaxnAuxMPI,&
     & GindexToLocal,nbasisAux)
end subroutine II_get_RI_AlphaCD_3CenterIntFullOnAllNNdim

subroutine DF3centerTrans4(nvirt,nocc,nbast,nAuxA,Cvirt,AlphaCD2,FullAlphaCD)
  implicit none
  integer,intent(in) :: nvirt,nocc,nbast,nAuxA
  real(realk),intent(in) :: Cvirt(nbast,nvirt),AlphaCD2(nAuxA,nBast,nocc)
  real(realk),intent(inout) :: FullAlphaCD(nAuxA,nvirt,nocc)
  !
  integer ::ADIAG,IDIAG,ALPHAAUX,ALPHA
  real(realk) :: TMP
  !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) &
  !$OMP PRIVATE(ADIAG,IDIAG,ALPHAAUX,ALPHA,TMP) &
  !$OMP SHARED(nvirt,nocc,nbast,nAuxA,Cvirt,AlphaCD2,FullAlphaCD) 
  do IDIAG = 1,nocc
     do ADIAG = 1,nvirt
        do ALPHAAUX = 1,nAUXA           
           FullAlphaCD(ALPHAAUX,ADIAG,IDIAG) = 0.0E0_realk
        enddo
        do ALPHA = 1,nbast
           TMP = Cvirt(ALPHA,ADIAG)
           do ALPHAAUX = 1,nAUXA           
              FullAlphaCD(ALPHAAUX,ADIAG,IDIAG) = FullAlphaCD(ALPHAAUX,ADIAG,IDIAG) + TMP*AlphaCD2(ALPHAAUX,ALPHA,IDIAG)
           enddo
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO
end subroutine DF3centerTrans4

!subroutine ass_D5to3TK(elms,ptr,dims)
!    implicit none
!    integer,intent(in) :: dims(3)
!    real(realk),target,intent(in) :: elms(dims(1),1,dims(2),dims(3),1)
!    real(realk),pointer:: ptr(:,:,:,:,:)
!    ptr => elms
!end subroutine ass_D5to3TK

subroutine GetAOtypes(ao,nAO,intspec)
  implicit none
  integer,intent(in) :: nAO
  integer,intent(inout) :: ao(nAO)
  character,intent(in)  :: intspec(nAO)
  integer :: i
  ! ***** SELECT AO TYPES *****
  DO i=1,nAO
     IF (intSpec(i).EQ.'R') THEN
        !   The regular AO-basis
        ao(i) = AORegular
     ELSEIF(intSpec(i).EQ.'D')THEN
        !   The Aux AO-type basis
        ao(i) = AOdfAux
     ELSE IF (intSpec(i).EQ.'C') THEN
        !   The CABS AO-type basis
        ao(i) = AOdfCABS
     ELSE
        call lsquit('Error in specification of ao1 in GetAOtypes',-1)
     ENDIF
  ENDDO
end subroutine GetAOtypes

subroutine GetOperatorFromCharacter(Oper,intspec,Setting)
implicit none
TYPE(LSSETTING)       :: SETTING
character,intent(in)  :: intspec
integer,intent(inout)  :: Oper

IF (intSpec.EQ.'C') THEN
   ! Regular Coulomb operator 1/r12
   oper = CoulombOperator
ELSE
   IF (intSpec.EQ.'E') THEN
      ! Long-Range operator erf/r12
      oper = ErfOperator
   ELSE
      IF (intSpec.EQ.'G') THEN
         oper = GGemOperator
      ELSE IF (intSpec.EQ.'F') THEN
         oper = GGemCouOperator
      ELSE IF (intSpec.EQ.'D') THEN
         oper = GGemGrdOperator
      ELSE IF (intSpec.EQ.'2') THEN
         oper = GGemQuaOperator
      ELSE
         print*,'intspec:',intspec
         call lsquit('Error in specification of operator in GetOperatorFromCharacter',-1)
      ENDIF
      call setF12Operator(Oper,Setting) 
   ENDIF
ENDIF
end subroutine GetOperatorFromCharacter

subroutine setF12Operator(Oper,Setting) 
  implicit none
  TYPE(LSSETTING)       :: SETTING
  integer,intent(in)  :: Oper
!
  integer             :: i,j,k,l
  real(realk)         :: coeff(6),exponent(6),tmp
  real(realk)         :: coeff2(21),sumexponent(21),prodexponent(21)
  integer             :: IJ,nGaussian,nG2,ao(4),dummy
        
  IF ((oper.EQ.GGemOperator.OR.oper .EQ. GGemCouOperator).OR.&
       & (oper .EQ. GGemGrdOperator.OR.oper .EQ. GGemQuaOperator))THEN
     nGaussian = 6
     nG2 = nGaussian*(nGaussian+1)/2
     call stgfit(1E0_realk,nGaussian,exponent,coeff)
     IJ=0
     DO I=1,nGaussian
        DO J=1,I
           IJ = IJ + 1
           coeff2(IJ) = 2E0_realk * coeff(I) * coeff(J)
           prodexponent(IJ) = exponent(I) * exponent(J)
           sumexponent(IJ) = exponent(I) + exponent(J)
        ENDDO
        coeff2(IJ) = 0.5E0_realk*coeff2(IJ)
     ENDDO
     IF (oper.EQ.GGemOperator) THEN
        ! The Gaussian geminal operatorg 
        call set_GGem(Setting%GGem,coeff,exponent,nGaussian)
     ELSE IF (oper .EQ. GGemCouOperator) THEN
        ! The Gaussian geminal divided by the Coulomb operator g/r12
        call set_GGem(Setting%GGem,coeff,exponent,nGaussian)
     ELSE IF (oper .EQ. GGemGrdOperator) THEN
        ! The double commutator [[T,g],g]
        call set_GGem(Setting%GGem,coeff2,sumexponent,prodexponent,nG2)
     ELSE IF (oper .EQ. GGemQuaOperator) THEN
        ! The Gaussian geminal operator squared g^2
        call set_GGem(Setting%GGem,coeff2,sumexponent,prodexponent,nG2)
     ELSE
        call lsquit('Error in specification of operator in setF12Operator',-1)
     ENDIF
  ENDIF
end subroutine setF12Operator

subroutine DF3centerTrans(nvirt,nocc,nbast,nAuxA,MynbasisAuxMPI,Cvirt,AlphaCD2,FullAlphaCD,startF)
  implicit none
  integer,intent(in) :: nvirt,nocc,nbast,nAuxA,MynbasisAuxMPI,startF
  real(realk),intent(in) :: Cvirt(nbast,nvirt),AlphaCD2(nAuxA,nBast,nocc)
  real(realk),intent(inout) :: FullAlphaCD(MynbasisAuxMPI,nvirt,nocc)
  !
  integer ::BDIAG,IDIAG,ALPHAAUX,ALPHA
  real(realk) :: TMP 
  !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) &
  !$OMP PRIVATE(BDIAG,IDIAG,ALPHAAUX,ALPHA,TMP) &
  !$OMP SHARED(nvirt,nocc,nbast,nAuxA,Cvirt,AlphaCD2,FullAlphaCD,startF) 
  do BDIAG = 1,nvirt
     do IDIAG = 1,nocc
        do ALPHAAUX = 1,nAUXA
           TMP = 0.0E0_realk
           do ALPHA = 1,nbast
              TMP = TMP + Cvirt(ALPHA,BDIAG)*AlphaCD2(ALPHAAUX,ALPHA,IDIAG)
           enddo
           FullAlphaCD(startF + ALPHAAUX,BDIAG,IDIAG) = FullAlphaCD(startF + ALPHAAUX,BDIAG,IDIAG) + TMP
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO
end subroutine DF3centerTrans

subroutine getRIbasisMPI(molecule1,nAtoms,numnodes,nbasisAuxMPI,&
     & startAuxMPI,AtomsMPI,nAtomsMPI,nAuxMPI)
  implicit none
  TYPE(MoleculeInfo),intent(in)    :: molecule1
  integer,intent(in)    :: numnodes,nAtoms
  integer,intent(inout) :: nbasisAuxMPI(numnodes),nAtomsMPI(numnodes)
  integer,intent(inout) :: startAuxMPI(nAtoms,numnodes)
  integer,intent(inout) :: AtomsMPI(nAtoms,numnodes),nAuxMPI(nAtoms,numnodes)
  !
  TYPE(MOLECULARORBITALINFO) :: orbitalInfo
  integer :: mynum2,iAtomA,nBast,nBastAux,idx(1),nb
  integer :: nBastLocA,startRegA,endRegA,nauxA,startAuxA,endAuxA     
  call setMolecularOrbitalInfo(molecule1,orbitalInfo)     
  DO mynum2 = 1,numnodes        
     nbasisAuxMPI(mynum2) = 0
     nAtomsMPI(mynum2) = 0
  ENDDO
  nb = 0
  DO iAtomA=1,nAtoms
     call getAtomicOrbitalInfo(orbitalInfo,iAtomA,nBastLocA,startRegA,endRegA,nauxA,startAuxA,endAuxA)
     idx = MINLOC(nbasisAuxMPI) !the rank this atom should be assign to. 
     nbasisAuxMPI(idx(1)) = nbasisAuxMPI(idx(1)) + nAuxA
     nAtomsMPI(idx(1)) = nAtomsMPI(idx(1)) + 1
     AtomsMPI(nAtomsMPI(idx(1)),idx(1)) = iAtomA
     nAuxMPI(nAtomsMPI(idx(1)),idx(1)) = nAuxA
     startAuxMPI(nAtomsMPI(idx(1)),idx(1)) = nb
     nb = nb + nauxA
  ENDDO
  CALL freeMolecularOrbitalInfo(orbitalInfo)
end subroutine getRIbasisMPI

subroutine getMaxAtomicnAux(molecule1,MaxAtomicnAux,nAtoms)
  implicit none
  TYPE(MoleculeInfo),intent(in)    :: molecule1
  integer,intent(in) :: nAtoms
  integer,intent(inout) :: MaxAtomicnAux
  !
  TYPE(MOLECULARORBITALINFO) :: orbitalInfo
  integer :: iAtomA,nBastLocA,startRegA,endRegA,nauxA,startAuxA,endAuxA
  call setMolecularOrbitalInfo(molecule1,orbitalInfo)     
  MaxAtomicnAux = 0
  DO iAtomA=1,nAtoms
     call getAtomicOrbitalInfo(orbitalInfo,iAtomA,nBastLocA,startRegA,endRegA,nauxA,startAuxA,endAuxA)
     MaxAtomicnAux = MAX(MaxAtomicnAux,nAuxA)
  ENDDO
  CALL freeMolecularOrbitalInfo(orbitalInfo)
end subroutine getMaxAtomicnAux

SUBROUTINE II_get_RI_AlphaBeta_2CenterInt(LUPRI,LUERR,AlphaBeta,SETTING,nbasisAux,InOper)
  IMPLICIT NONE
  Integer,intent(in)                :: LUPRI,LUERR,nbasisAux
  REAL(REALK)                       :: AlphaBeta(nbasisAux,nbasisAux)
  TYPE(LSSETTING),intent(inout)     :: SETTING
  integer,optional                  :: InOper
  !
  Integer                    :: nAtoms,nBastAux,nBast,Oper
  Real(realk) :: TSTART,TEND,tmp
  logical :: MasterWakeSlaves,doMPI
  call LSTIMER('START ',TSTART,TEND,LUPRI)
  call getMolecularDimensions(SETTING%MOLECULE(1)%p,nAtoms,nBast,nBastAux)
  IF(nbasisAux.NE.nBastAux)THEN
     print*,'nbasisAux',nbasisAux
     print*,'nbastAux',nbastAux
     Call lsquit('dim mismatch in II_get_df_exchange2',-1)
  ENDIF
  !set threshold 
  SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%K_THR
  !(alpha|cd)
  !call typedef_setMolecules(setting,molecule,1,2,3,4)
  call initIntegralOutputDims(setting%output,nbastaux,1,nbastaux,1,1)
  !MPI is used inside this routine 
  !but since this is called by both master and slaves
  !the master should not wake up the slaves
  doMPI = SETTING%SCHEME%doMPI
  SETTING%SCHEME%doMPI = .FALSE.
  MasterWakeSlaves = SETTING%SCHEME%MasterWakeSlaves
  SETTING%SCHEME%MasterWakeSlaves = .FALSE.
  IF(present(InOper))THEN
     Oper = InOper
  ELSE
     Oper = CoulombOperator
  ENDIF
  call setF12Operator(Oper,Setting) 
  call ls_getIntegrals(AODFdefault,AOempty,AODFdefault,AOempty,&
       & Oper,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
  SETTING%SCHEME%doMPI = doMPI
  SETTING%SCHEME%MasterWakeSlaves = MasterWakeSlaves
  call retrieve_Output(lupri,setting,AlphaBeta,.FALSE.)
  call LSTIMER('AlphaBeta',TSTART,TEND,LUPRI)
END SUBROUTINE II_get_RI_AlphaBeta_2CenterInt

SUBROUTINE II_get_RI_AlphaBeta_geoderiv2CenterInt(LUPRI,LUERR,AlphaBetaDeriv,&
     & SETTING,nbasisAux,natoms)
  IMPLICIT NONE
  Integer,intent(in)                :: LUPRI,LUERR,nbasisAux,natoms
  REAL(REALK)                       :: AlphaBetaDeriv(nbasisAux,nbasisAux,3*natoms)
  TYPE(LSSETTING),intent(inout)     :: SETTING
  !
  Integer     :: nBastAux,nBast,natoms2,I
  Real(realk) :: TSTART,TEND,tmp
  logical :: MasterWakeSlaves,doMPI
  REAL(REALK),pointer :: AlphaBetaDeriv2(:,:,:)
  integer,pointer :: BACKUPmolindex1(:),BACKUPmolindex2(:)
  
  call LSTIMER('START ',TSTART,TEND,LUPRI)
  call getMolecularDimensions(SETTING%MOLECULE(1)%p,nAtoms2,nBast,nBastAux)
  call mem_alloc(BACKUPmolindex1,nAtoms2)
  call mem_alloc(BACKUPmolindex2,nAtoms2)
  !The molecule contains molecule%p%ATOM(i)%molecularIndex pointing to the
  !correct atom in full molecule but in order to only allocate (nbastAux,nbastAux,3*natoms_frag)
  !instead of (nbastAux,nbastAux,3*natoms) we modify the molecularindex
  do i = 1,size(SETTING%MOLECULE(1)%p%ATOM)
     BACKUPmolindex1(i) = SETTING%MOLECULE(1)%p%ATOM(i)%molecularIndex 
     BACKUPmolindex2(i) = SETTING%MOLECULE(2)%p%ATOM(i)%molecularIndex 
     SETTING%MOLECULE(1)%p%ATOM(i)%molecularIndex = i
     SETTING%MOLECULE(2)%p%ATOM(i)%molecularIndex = i
  enddo

  IF(nbasisAux.NE.nBastAux)THEN
     print*,'nbasisAux',nbasisAux
     print*,'nbastAux',nbastAux
     Call lsquit('dim mismatch in II_get_RI_AlphaBeta_geoderiv2CenterInt',-1)
  ENDIF
  IF(natoms.NE.natoms2)Call lsquit('atoms mismatch in II_get_RI_AlphaBeta_geoderiv2CenterInt',-1)
  !set threshold 
  SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%K_THR
  !(alpha|cd)
  !call typedef_setMolecules(setting,molecule,1,2,3,4)
  call initIntegralOutputDims(setting%output,nbastaux,1,nbastaux,1,3*natoms)
  !MPI is used inside this routine 
  !but since this is called by both master and slaves
  !the master should not wake up the slaves
  doMPI = SETTING%SCHEME%doMPI
  SETTING%SCHEME%doMPI = .FALSE.
  MasterWakeSlaves = SETTING%SCHEME%MasterWakeSlaves
  SETTING%SCHEME%MasterWakeSlaves = .FALSE.

  !(P^x|Q)
  call ls_getIntegrals(AODFdefault,AOempty,AODFdefault,AOempty,&
       & CoulombOperator,GeoDerivLHSSpec,ContractedInttype,SETTING,LUPRI,LUERR)
  call retrieve_Output(lupri,setting,AlphaBetaDeriv,.FALSE.)
  !(P|Q^x)
  call initIntegralOutputDims(setting%output,nbastaux,1,nbastaux,1,3*natoms)
  call ls_getIntegrals(AODFdefault,AOempty,AODFdefault,AOempty,&
       & CoulombOperator,GeoDerivRHSSpec,ContractedInttype,SETTING,LUPRI,LUERR)
  call mem_alloc(AlphaBetaDeriv2,nbasisAux,nbasisAux,3*natoms)
  call retrieve_Output(lupri,setting,AlphaBetaDeriv2,.FALSE.)
  call dcopy(nbastaux*nbastaux*3*natoms,AlphaBetaDeriv2,1,AlphaBetaDeriv,1)
  call mem_dealloc(AlphaBetaDeriv2)
  SETTING%SCHEME%doMPI = doMPI
  SETTING%SCHEME%MasterWakeSlaves = MasterWakeSlaves
  do i = 1,size(SETTING%MOLECULE(1)%p%ATOM)
     SETTING%MOLECULE(1)%p%ATOM(i)%molecularIndex = BACKUPmolindex1(i)
     SETTING%MOLECULE(2)%p%ATOM(i)%molecularIndex = BACKUPmolindex2(i)
  enddo
  call mem_dealloc(BACKUPmolindex1)
  call mem_dealloc(BACKUPmolindex2)

  call LSTIMER('AlphaBetaDeriv',TSTART,TEND,LUPRI)
END SUBROUTINE II_get_RI_AlphaBeta_geoderiv2CenterInt

!3. Calculate 2 of 3 RIMP2 gradient contributions: 
! A. (P^x|bj)*CalphaTheta(P,b,j)
! B. (P|(beta nu)^x)*Cvirt(beta,b)*Cocc(nu,J)*CalphaTheta(P,b,j)
! Loop over P shell (batches) belonging to this node
!      for each P shell:
!      A. 
!      Calculate (P^x|bj)
!      Contract grad(x) = (P^x|bj)*CalphaTheta(P,b,j)
!      B.
!      Contract  CalphaThetaAO(Pshell,beta,nu) = Cvirt(beta,b)*Cocc(nu,J)*CalphaTheta(P,b,j)
!      Calculate (P|(beta nu)^x)
!      Contract grad(x) = (P|(beta nu)^x)*CalphaThetaAO(Pshell,beta,nu)
! End Loop over P shell
SUBROUTINE II_get_RIMP2_grad(LUPRI,LUERR,Gradient,SETTING,nbasisAux,&
     & nbasis,nvirt,nocc,Cvirt,Cocc,maxsize,mynum,numnodes,natoms,&
     & nAuxLoc,CfitLHS,InOper)
  IMPLICIT NONE
  Integer,intent(in)            :: LUPRI,LUERR,nbasis,nbasisAux,nAuxLoc
  Integer,intent(in)            :: nocc,nvirt,mynum,numnodes,natoms
  REAL(REALK),intent(inout)     :: Gradient(3,natoms)
  REAL(REALK),intent(in)        :: Cocc(nbasis,nocc)
  REAL(REALK),intent(in)        :: Cvirt(nbasis,nvirt)
  REAL(REALK),intent(in)        :: CfitLHS(nAuxLoc,nocc,nvirt)
  TYPE(LSSETTING),intent(inout) :: SETTING
  integer(kind=long),intent(in) :: maxsize
  integer,optional              :: InOper 
  !
  integer(kind=long)         :: nsize
  Integer                    :: nAtomsAux,nBastAux,nBast,N,K,M,ialpha,v,a,natoms2
  Integer                    :: BDIAG,IDIAG,ILOC,JLOC,ALPHAAUX,GAMMA,DELTA
  Integer                    :: ALPHA,BETA,I,Oper 
  Real(realk) :: TSTART,TEND,tmp,TMP1
  real(realk),pointer :: AlphaCD(:,:,:,:,:),AlphaCDmo1(:,:,:,:)
  real(realk),pointer :: AlphaCDmo(:,:,:,:)
  TYPE(MoleculeInfo),pointer      :: molecule1,molecule2,molecule3,molecule4
  TYPE(MoleculeInfo),pointer      :: ATOMS(:)
  TYPE(MOLECULARORBITALINFO) :: orbitalInfo
  Integer :: iAtomA,nAuxA,B,offsetLoc,offsetFull
  integer :: J,mynum2,startF,iatomampi,MynbasisAuxMPI,nBastAuxT
  logical :: doMPI,MasterWakeSlaves
  integer,pointer :: nbasisAuxMPI(:),startAuxMPI(:,:),AtomsMPI(:,:),nAtomsMPI(:),nAuxMPI(:,:)
  !
  integer :: iShell, nAuxShellA,nbatches,X
  integer, pointer :: batchdim(:),BACKUPmolindex(:)
!  REAL(REALK)     :: GradientA(3,natoms)
!  REAL(REALK)     :: GradientB(3,natoms)
!  GradientA = 0.0E0_realk
!  GradientB = 0.0E0_realk
  Gradient = 0.0E0_realk

  IF(present(InOper))THEN
     Oper = InOper
  ELSE
     Oper = CoulombOperator
  ENDIF
  SETTING%scheme%CS_SCREEN = .FALSE.
  SETTING%scheme%PS_SCREEN = .FALSE.
  SETTING%scheme%OD_SCREEN = .FALSE.
  SETTING%scheme%OE_SCREEN = .FALSE.

  nsize = (nbasis*nbasis*nbasisAux+nocc*nbasis*nbasisAux)*mem_realsize
  !in case of setting%molBuild = .TRUE. then the setting%molecule(1)%p can 
  !not be used as pointers

  molecule1 => SETTING%MOLECULE(1)%p
  molecule2 => SETTING%MOLECULE(2)%p
  molecule3 => SETTING%MOLECULE(3)%p
  molecule4 => SETTING%MOLECULE(4)%p
  call mem_alloc(BACKUPmolindex,size(molecule3%ATOM))
  call LSTIMER('START ',TSTART,TEND,LUPRI)

  !need to deactivate MPI inside ls_getIntegrals. 
  !both master and slaves call this routine 
  doMPI = SETTING%SCHEME%doMPI
  MasterWakeSlaves = SETTING%SCHEME%MasterWakeSlaves
  SETTING%SCHEME%doMPI = .FALSE.
  SETTING%SCHEME%MasterWakeSlaves = .FALSE.
  
  !MPI version and memory reduced version
  call getMolecularDimensions(molecule1,nAtomsAux,nBast,nBastAux) !Aux
  call getMolecularDimensions(molecule3,nAtoms2,nBast,nBastAuxT)!Reg 
  IF(nAtomsAux.NE.nAtoms)call lsquit('natoms mismatch in II_get_RIMP2_grad',-1)
  allocate(ATOMS(nAtomsAux))
  CALL pari_set_atomic_fragments(molecule1,ATOMS,nAtomsAux,lupri)
  
  ! Split only on of the atomic loops over nodes
  call mem_alloc(nbasisAuxMPI,numnodes)
  call mem_alloc(nAtomsMPI,numnodes)    
  call mem_alloc(startAuxMPI,nAtomsAux,numnodes)
  call mem_alloc(AtomsMPI,nAtomsAux,numnodes)
  call mem_alloc(nAuxMPI,nAtomsAux,numnodes)
  call getRIbasisMPI(molecule1,nAtomsAux,numnodes,nbasisAuxMPI,startAuxMPI,&
       & AtomsMPI,nAtomsMPI,nAuxMPI)
  MynbasisAuxMPI = nbasisAuxMPI(mynum+1)
  IF(MynbasisAuxMPI.NE.nAuxLoc) call lsquit('dim mismatch in II_get_RIMP2_grad',-1)
  call mem_dealloc(nbasisAuxMPI)
  IF(MynbasisAuxMPI.GT.0)THEN
   DO iAtomAMPI=1,nAtomsMPI(mynum+1)
      iAtomA = AtomsMPI(iAtomAMPI,mynum+1)
      offsetFull = startAuxMPI(iAtomA,mynum+1)
      nAuxA = nAuxMPI(iAtomAMPI,mynum+1)
      !While the Atom ATOMS(iAtomA)%ATOM(1) is iAtomA we set it to 1 
      !so that we can build Integral(1:nAuxShellA,1,1:nbast,1:nbast,1:3)
      !instead of Integral(1:nAuxShellA,1,1:nbast,1:nbast,1+(iAtomA-1)*3:3+(iAtomA-1)*3)
      !Naturally we place it the correct place in the gradient in the subroutine
      !RIMP2_Grad_Contract
      ATOMS(iAtomA)%ATOM(1)%molecularIndex=1
      !molecule3 also contain molecule3%ATOM(iatom)%molecularIndex pointing to the
      !correct atom but in order to only allocate (nAuxShellA,nbast,nbast,3*natoms_frag)
      !instead of (nAuxShellA,nbast,nbast,3*natoms) we modify the molecularindex
      do iShell = 1,size(molecule3%ATOM)
         BACKUPmolindex(iShell) = molecule3%ATOM(iShell)%molecularIndex 
         molecule3%ATOM(iShell)%molecularIndex = iShell
      enddo
      call typedef_setMolecules(setting,molecule3,3,4,ATOMS(iAtomA),1)             
      
      nullify(batchdim)
      call build_minimalbatchesOfAOs2(lupri,setting,batchdim,nbatches,'D')
      offsetLoc = 0
      BatchD: do iShell = 1,nbatches
         nAuxShellA = batchdim(iShell)
         !================================================================
         !A.    grad(x) = (P^x|bj)*CalphaTheta(P,b,j)
         !================================================================
         call initIntegralOutputDims(setting%Output,nAuxShellA,1,nbast,nbast,3)
         setting%batchindex(1)=iShell
         setting%batchdim(1)=nAuxShellA
         !Calculat (P^x|beta nu)
         call ls_getIntegrals(AODFdefault,AOEmpty,AORdefault,AORdefault,Oper,&
              & GeoDerivLHSSpec,ContractedInttype,SETTING,LUPRI,LUERR)
         setting%batchindex(1)=0
         setting%batchdim(1)=nAuxShellA
         call mem_alloc(alphaCD,nAuxShellA,1,nbast,nbast,3)
         alphaCD = 0.0E0_realk !obsolete?
         CALL retrieve_Output(lupri,setting,alphaCD,.FALSE.)
         call mem_alloc(alphaCDmo1,nAuxShellA,nbast,nocc,3)
         !AO to MO (P^x|beta nu) -> (P^x|beta j)  
         call AOtoMO_Pgrad4cent(alphaCD,alphaCDmo1,nAuxShellA,nbast,Cocc,nocc,3)
         call mem_dealloc(AlphaCD)
         call mem_alloc(alphaCDmo,nAuxShellA,nvirt,nocc,3)
         !AO to MO (P^x|beta j) -> (P^x|b j)  
         call AOtoMO_Pgrad3cent(alphaCDmo1,alphaCDmo,nAuxShellA,nbast,nocc,Cvirt,nvirt,3)
         call mem_dealloc(AlphaCDmo1)
         startF = offsetFull+offsetLoc
         !print*,'startF',startF
         !Contract grad(x) = (P^x|bj)*CalphaTheta(P,j,b)
         call RIMP2_Grad_ContractA(nvirt,nocc,nAuxShellA,nBasisAux,&
              & CfitLHS,alphaCDmo,Gradient,startF,iAtomA,natoms)
!         call RIMP2_Grad_ContractA(nvirt,nocc,nAuxShellA,nBasisAux,&
!              & CfitLHS,alphaCDmo,GradientA,startF,iAtomA,natoms)

         call mem_dealloc(AlphaCDmo)

         !================================================================
         !B.    grad(x) = (P|(bj)^x)*CalphaTheta(P,b,j)
         !================================================================
         !Calculate (P^x|bj)
         !WARNING: May require Batching over one of the AOs! Or Atomic Batching
         !         Due to the 3*natoms geoderiv components.  
         call initIntegralOutputDims(setting%Output,nAuxShellA,1,nbast,nbast,3*natoms)
         setting%batchindex(1)=iShell
         setting%batchdim(1)=nAuxShellA
         !Calculat (P|(beta nu)^x)
         call ls_getIntegrals(AODFdefault,AOEmpty,AORdefault,AORdefault,Oper,&
              & GeoDerivRHSSpec,ContractedInttype,SETTING,LUPRI,LUERR)
         setting%batchindex(1)=0
         setting%batchdim(1)=nAuxShellA
         call mem_alloc(alphaCD,nAuxShellA,1,nbast,nbast,3*natoms)
         alphaCD = 0.0E0_realk !obsolete?
         CALL retrieve_Output(lupri,setting,alphaCD,.FALSE.)

         call mem_alloc(alphaCDmo1,nAuxShellA,nbast,nocc,3*natoms)
         !AO to MO (P|(beta nu)^x) -> (P^x|(beta j)^x)  
         call AOtoMO_Pgrad4cent(alphaCD,alphaCDmo1,nAuxShellA,nbast,Cocc,nocc,3*natoms)
         call mem_dealloc(AlphaCD)
         call mem_alloc(alphaCDmo,nAuxShellA,nvirt,nocc,3*natoms)
         !AO to MO (P|(beta j)^x) -> (P|(b j)^x)  
         call AOtoMO_Pgrad3cent(alphaCDmo1,alphaCDmo,nAuxShellA,nbast,nocc,Cvirt,nvirt,3*natoms)
         call mem_dealloc(AlphaCDmo1)
         !Contract grad(x) = (P|(bj)^x)*CalphaTheta(P,j,b)
         call RIMP2_Grad_ContractB(nvirt,nocc,nAuxShellA,nBasisAux,&
              & CfitLHS,alphaCDmo,Gradient,startF,natoms)
!         call RIMP2_Grad_ContractB(nvirt,nocc,nAuxShellA,nBasisAux,&
!              & CfitLHS,alphaCDmo,GradientB,startF,natoms)
         call mem_dealloc(AlphaCDmo)

         offsetLoc = offsetLoc + nAuxShellA !offset in Local Aux
      ENDDO BatchD
      call mem_dealloc(batchdim)

!      print*,'PARTIAL A'
!      call ls_output(GradientA,1,3,1,natoms,3,natoms,1,6)      
!      print*,'PARTIAL B'
!      call ls_output(GradientB,1,3,1,natoms,3,natoms,1,6)      
!      print*,'PARTIAL A+B'
!      call ls_output(Gradient,1,3,1,natoms,3,natoms,1,6)
      do iShell = 1,size(molecule3%ATOM)
         molecule3%ATOM(iShell)%molecularIndex = BACKUPmolindex(iShell) 
      enddo
   ENDDO
   !restore 
   call typedef_setMolecules(setting,molecule1,1,molecule2,2,&
        & molecule3,3,molecule4,4)
   call pari_free_atomic_fragments(ATOMS,nAtomsAux)        
  ENDIF
  call mem_dealloc(BACKUPmolindex)
  deallocate(ATOMS)
  call mem_dealloc(AtomsMPI)
  call mem_dealloc(nAtomsMPI)
  call mem_dealloc(nAuxMPI)
  call mem_dealloc(startAuxMPI)
  SETTING%SCHEME%doMPI = doMPI
  SETTING%SCHEME%MasterWakeSlaves = MasterWakeSlaves
  !restore 
  SETTING%MOLECULE(1)%p => molecule1
  SETTING%MOLECULE(2)%p => molecule2
  SETTING%MOLECULE(3)%p => molecule3
  SETTING%MOLECULE(4)%p => molecule4
  call ls_setDefaultFragments(setting)
  call LSTIMER('AlphaCD',TSTART,TEND,LUPRI)
  SETTING%scheme%CS_SCREEN = .TRUE.
  SETTING%scheme%PS_SCREEN = .TRUE.
  SETTING%scheme%OD_SCREEN = .TRUE.
  SETTING%scheme%OE_SCREEN = .TRUE.
END SUBROUTINE II_get_RIMP2_grad

subroutine AOtoMO_Pgrad4cent(alphaCD,alphaCDmo,nAuxShellA,nbast,Cocc,nocc,ngradcomp)
  implicit none
integer,intent(in)        :: nAuxShellA,nbast,nocc,ngradcomp
real(realk),intent(in)    :: AlphaCD(nAuxShellA*nbast,nbast,ngradcomp)
real(realk),intent(in)    :: Cocc(nbast,nocc)
real(realk),intent(inout) :: alphaCDmo(nAuxShellA*nbast,nocc,ngradcomp)
!
integer :: I,PBETA,NU,X
real(realk) :: TMP
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) &
!$OMP PRIVATE(I,PBETA,NU,X,TMP) &
!$OMP SHARED(nocc,nbast,nAuxShellA,Cocc,AlphaCD,AlphaCDmo,ngradcomp) 
DO X=1,ngradcomp
   DO I=1,nocc
      TMP = Cocc(1,I)
      DO PBETA=1,nAuxShellA*nbast
         alphaCDmo(PBETA,I,X) = AlphaCD(PBETA,1,X)*TMP
      ENDDO
      DO NU=2,nbast
         TMP = Cocc(NU,I)
         DO PBETA=1,nAuxShellA*nbast
            alphaCDmo(PBETA,I,X) = alphaCDmo(PBETA,I,X) + AlphaCD(PBETA,NU,X)*TMP
         ENDDO
      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO
end subroutine AOtoMO_Pgrad4cent

subroutine AOtoMO_Pgrad3cent(alphaCD,alphaCDmo,nAuxShellA,nbast,nocc,Cvirt,nvirt,ngradcomp)
  implicit none
integer,intent(in)        :: nAuxShellA,nbast,nvirt,nocc,ngradcomp
real(realk),intent(in)    :: AlphaCD(nAuxShellA,nbast,nocc*ngradcomp)
real(realk),intent(in)    :: Cvirt(nbast,nvirt)
real(realk),intent(inout) :: alphaCDmo(nAuxShellA,nvirt,nocc*ngradcomp)
!
integer :: B,P,BETA,X
real(realk) :: TMP
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) &
!$OMP PRIVATE(B,P,BETA,X,TMP) &
!$OMP SHARED(nvirt,nocc,nbast,nAuxShellA,Cvirt,AlphaCD,AlphaCDmo,ngradcomp) 
DO X=1,ngradcomp*nocc
   DO B=1,nvirt
      TMP = Cvirt(1,B)
      DO P=1,nAuxShellA
         alphaCDmo(P,B,X) = AlphaCD(P,1,X)*TMP
      ENDDO
      DO BETA=2,nbast
         TMP = Cvirt(BETA,B)
         DO P=1,nAuxShellA
            alphaCDmo(P,B,X) = alphaCDmo(P,B,X) + AlphaCD(P,BETA,X)*TMP
         ENDDO
      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO
end subroutine AOtoMO_Pgrad3cent


!grad(1:3,iAtomAux) =+ CfitLHS(nAux,a,i)*alphaCDmo(nAuxShellA,a,i,3)
subroutine RIMP2_Grad_ContractA(nvirt,nocc,nAuxShellA,nAux,&
     & CfitLHS,alphaCDmo,Gradient,startF,iAtomA,natoms)
  implicit none
  Integer,intent(in)            :: nocc,nvirt,natoms,startF,nAux,iAtomA,nAuxShellA
  REAL(REALK),intent(inout)     :: Gradient(3,natoms)
  REAL(REALK),intent(in)        :: CfitLHS(nAux,nocc,nvirt)
  REAL(REALK),intent(in)        :: alphaCDmo(nAuxShellA,nvirt,nocc,3)
  !
  integer :: X,A,I,ALPHA
  real(realk) :: TMP
  !$OMP PARALLEL DEFAULT(none) PRIVATE(X,A,I,ALPHA,TMP) SHARED(nvirt,&
  !$OMP nocc,nAuxShellA,startF,CfitLHS,alphaCDmo,iAtomA,&
  !$OMP Gradient)
  do X=1,3
     TMP = 0.0E0_realk
     !$OMP DO COLLAPSE(2)
     do A=1,nvirt
        do I=1,nocc
           do ALPHA=1,nAuxShellA
              TMP = TMP + CfitLHS(startF + ALPHA,I,A)*alphaCDmo(ALPHA,A,I,X)
           enddo
        enddo
     enddo
     !$OMP END DO
     !$OMP ATOMIC
     Gradient(X,iAtomA) = Gradient(X,iAtomA) + TMP
  enddo
  !$OMP END PARALLEL

end subroutine RIMP2_Grad_ContractA

!grad(1:3,iAtomAux) =+ CfitLHS(nAux,a,i)*alphaCDmo(nAuxShellA,a,i,3)
subroutine RIMP2_Grad_ContractB(nvirt,nocc,nAuxShellA,nAux,&
     & CfitLHS,alphaCDmo,Gradient,startF,natoms)
  implicit none
  Integer,intent(in)            :: nocc,nvirt,natoms,startF,nAux,nAuxShellA
  REAL(REALK),intent(inout)     :: Gradient(3*natoms)
  REAL(REALK),intent(in)        :: CfitLHS(nAux,nocc,nvirt)
  REAL(REALK),intent(in)        :: alphaCDmo(nAuxShellA,nvirt,nocc,3*natoms)
  !
  integer :: X,A,I,ALPHA
  real(realk) :: TMP
  !$OMP PARALLEL DEFAULT(none) PRIVATE(X,A,I,ALPHA,TMP) SHARED(nvirt,&
  !$OMP nocc,nAuxShellA,startF,CfitLHS,alphaCDmo,&
  !$OMP Gradient,natoms)
  do X=1,3*natoms
     TMP = 0.0E0_realk
     !$OMP DO COLLAPSE(2)
     do A=1,nvirt
        do I=1,nocc
           do ALPHA=1,nAuxShellA
              TMP = TMP + CfitLHS(startF + ALPHA,I,A)*alphaCDmo(ALPHA,A,I,X)
           enddo
        enddo
     enddo
     !$OMP END DO
     !$OMP ATOMIC
     Gradient(X) = Gradient(X) + TMP
  enddo
  !$OMP END PARALLEL

end subroutine RIMP2_Grad_ContractB

END MODULE INTEGRALINTERFACEMODULEDF
