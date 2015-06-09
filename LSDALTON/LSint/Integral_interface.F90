MODULE IntegralInterfaceMOD
  use precision
  use TYPEDEFTYPE, only: LSSETTING, LSINTSCHEME, LSITEM, integralconfig
  use Matrix_module, only: MATRIX, MATRIXP
  use LSparameters
  use LSTIMING
  use molecule_typetype, only: MOLECULE_PT, ATOMITEM
  use molecule_type, only: build_pointmolecule, DETERMINE_MAXCOOR, &
       & free_moleculeinfo
  use integral_type, only: INTEGRALINPUT
  use integraloutput_typetype, only: INTEGRALOUTPUT
  use integraloutput_type, only: initintegraloutputdims
  use lstensor_operationsmod, only: lstensor, lstensor_nullify
  use ao_typetype, only: aoitem, BATCHORBITALINFO
  use ao_type, only: free_aoitem, freebatchorbitalinfo, initbatchorbitalinfo, &
       & setbatchorbitalinfo
  use TYPEDEF, only: getNbasis, retrieve_output, gcao2ao_transform_matrixd2, &
       & retrieve_screen_output, ao2gcao_transform_matrixf, &
       & gcao2ao_transform_fulld, ao2gcao_transform_fullf, &
       & ao2gcao_half_transform_matrix,gcao2ao_half_transform_matrix,&
       & GCAO2AO_transform_matrixD2,typedef_setMolecules
  use KS_settings, only: SaveF0andD0,incrD0,incrF0,incremental_scheme,&
       & incrDdiff,do_increment,activate_incremental
  use ls_Integral_Interface, only: ls_same_mats, ls_getintegrals, &
       & ls_freedmatfromsetting, ls_attachdmattosetting, ls_multipolemoment, &
       & ls_get_exchange_mat, ls_jengineclassicalgrad, &
       & ls_jengine, ls_getnucscreenintegrals, setaobatch, &
       & ls_getscreenintegrals1,ls_setDefaultFragments, &
       & ls_get_coulomb_mat,ls_LHSSameAsRHSDmatToSetting,&
       & ls_LHSSameAsRHSDmatToSetting_deactivate,&
       & ls_get_coulomb_and_exchange_mat, ls_jengineclassicalmat
  use molecule_module, only: getMolecularDimensions
  use matrix_operations, only: mat_dotproduct, matrix_type, mtype_unres_dense,&
       & mat_daxpy, mat_init, mat_free, mat_write_to_disk, mat_print, mat_zero,&
       & mat_scal, mat_mul, mat_assign, mat_trans, mat_copy, mat_add, mat_trAB,&
       & mat_sqnorm2,mat_tr,mat_max_elm
  use matrix_util, only: mat_get_isym, util_get_symm_part,util_get_antisymm_part, matfull_get_isym, mcweeney_purify, &
                         util_get_symm_and_antisymm_part_full
  use memory_handling, only: mem_alloc, mem_dealloc
#ifdef VAR_MPI
  use lsmpi_op, only: mpicopy_screen
  use lsmpi_type!, only: ls_mpiFinalizeBuffer, ls_mpiInitBuffer, &
!       & LSMPIBROADCAST,  ls_mpibcast, get_rank_for_comm
  use infpar_module
#endif
  use IIDFTINT
  use IntegralInterfaceModuleDF
  use files, only: lsopen,lsclose
  use io!, only: io_get_filename
  use screen_mod,only: screen_add_associate_item, screen_init, screen_free
  use f12_module, only: free_ggem, set_ggem, stgfit
  use linsolvdf, only: linsolv_df
!  use IntegralInterfaceModuleDF, only: ii_get_df_j_gradient, ii_get_df_coulomb_mat, ii_get_df_exchange_mat, ii_get_pari_df_exchange_mat
  use SphCart_Matrices, only: BUILD_CART_TO_SPH_MAT
  use gridgenerationmodule
  use II_XC_interfaceModule
  use GCtransMod
  use dft_typetype
  use ls_util, only: ls_print_gradient
  public::  II_get_overlap, II_get_mixed_overlap,II_get_mixed_overlap_full,&
       & II_get_h1,II_get_h1_mixed,II_get_h1_mixed_full,II_get_kinetic,&
       & II_get_kinetic_mixed,II_get_kinetic_mixed_Full,II_get_nucel_mat,&
       & II_get_nucel_mat_mixed,II_get_nucel_mat_mixed_full,&
       & II_get_magderivOverlap,II_get_maggradOverlap,II_get_magderivOverlapR,&
       & II_get_magderivOverlapL,II_get_ep_integrals,II_get_ep_integrals2,&
       & II_get_ep_integrals3,II_get_ep_ab,II_GET_MOLECULAR_GRADIENT,&
       & II_get_twoElectron_gradient,II_get_K_gradient,II_get_K_gradientfull,&
       & II_get_regular_K_gradient,II_get_J_gradient,II_get_pso_spec_expval,&
       & II_get_J_gradient_regular,II_get_oneElectron_gradient,&
       & II_get_ne_gradient,II_get_kinetic_gradient,II_get_nucpot,&
       & II_get_nn_gradient,II_get_nuc_Quad,II_get_carmom,&
       & II_get_single_carmom,II_get_nucdip,II_get_prop,II_get_pso_spec,&
       & II_get_prop_expval,II_get_integral,II_get_integral_full,&
       & II_get_sphmom,II_carmom_to_shermom,II_get_3center_overlap,&
       & II_get_2center_eri,II_get_2center_mixed_eri,II_get_4center_eri,&
       & II_get_4center_eri_diff,II_get_1el_diff,II_precalc_ScreenMat, &
#ifdef VAR_MPI
       & II_bcast_screen, II_screeninit, II_screenfree,&
#endif
       & II_get_2int_ScreenMat,II_get_maxGabelm_ScreenMat,II_setIncremental,&
       & II_getBatchOrbitalInfo,II_get_reorthoNormalization,&
       & II_get_reorthoNormalization2,II_get_geoderivKinetic,&
       & II_get_geoderivnucel,II_get_geoderivOverlap,&
       & II_get_geoderivOverlapL,II_get_geoderivExchange,&
       & II_get_geoderivCoulomb,II_get_GaussianGeminalFourCenter,&
       & II_get_magderiv_4center_eri,II_get_magderivF,&
       & II_get_magderivK,II_get_magderivJ, II_get_Econt,II_get_exchangeEcont,&
       & II_get_CoulombEcont,II_get_ABres_4CenterEri,II_get_Fock_mat_full,&
       & II_get_coulomb_mat_full, II_get_coulomb_mat_mixed_full,&
       & II_get_jengine_mat_full, II_get_exchange_mat_full,&
       & ii_get_exchange_mat_mixed_full, II_get_exchange_mat1_full,&
       & II_get_exchange_mat_regular_full, II_get_admm_exchange_mat,get_T23,&
       & II_get_ADMM_K_gradient, II_get_coulomb_mat,ii_get_exchange_mat_mixed,&
       & II_get_exchange_mat,II_get_coulomb_and_exchange_mat, II_get_Fock_mat,&
       & II_get_coulomb_mat_mixed, II_GET_DISTANCEPLOT_4CENTERERI,&
       & II_get_2int_ScreenRealMat,transformed_f2_to_f3,transform_D3_to_D2,&
       & ii_get_2int_batchscreenmat, II_get_nst_spec_expval
  private

INTERFACE II_get_coulomb_mat
   MODULE PROCEDURE II_get_coulomb_mat_array,II_get_coulomb_mat_single
end INTERFACE II_get_coulomb_mat

INTERFACE II_get_exchange_mat
   MODULE PROCEDURE II_get_exchange_mat_array,II_get_exchange_mat_single
end INTERFACE II_get_exchange_mat

INTERFACE II_get_coulomb_and_exchange_mat
   MODULE PROCEDURE II_get_coulomb_and_exchange_mat_array,&
        &II_get_coulomb_and_exchange_mat_single
END INTERFACE II_get_coulomb_and_exchange_mat

INTERFACE II_get_Fock_mat
   MODULE PROCEDURE II_get_Fock_mat_array,II_get_Fock_mat_single
END INTERFACE II_get_Fock_mat

INTERFACE II_get_rho2
   MODULE PROCEDURE II_get_rho2_mat
END INTERFACE II_get_rho2

CONTAINS
!> \brief Calculates overlap integral matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param S the overlap matrix
SUBROUTINE II_get_overlap(LUPRI,LUERR,SETTING,S)
IMPLICIT NONE
TYPE(MATRIX)          :: S
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
!
Integer             :: nbast
real(realk)         :: TS,TE
logical             :: sameMolecule
call time_II_operations1()
CALL LSTIMER('START ',TS,TE,LUPRI)
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
nbast = S%nrow
IF(SETTING%SCHEME%DEBUGOVERLAP)THEN
  call initIntegralOutputDims(setting%Output,nbast,nbast,1,1,1)
  CALL ls_getIntegrals(AORdefault,AORdefault,AOempty,AOempty,&
       &OverlapOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
  IF(ASSOCIATED(setting%molecule(1)%p,setting%molecule(2)%p))THEN
     sameMolecule = .TRUE.
  ELSE
     sameMolecule = .FALSE.
  ENDIF
ELSE
   call initIntegralOutputDims(setting%Output,nbast,1,nbast,1,1)
   CALL ls_getIntegrals(AORdefault,AOempty,AORdefault,AOempty,&
        &OverlapOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
  IF(ASSOCIATED(setting%molecule(1)%p,setting%molecule(3)%p))THEN
     sameMolecule = .TRUE.
  ELSE
     sameMolecule = .FALSE.
  ENDIF
ENDIF
CALL retrieve_Output(lupri,setting,S,setting%IntegralTransformGC)
!it is very important that the Overlap matrix is symmetric due to 
!decomposition 
IF(sameMolecule)THEN
   call util_get_symm_part(S)
ENDIF
IF (setting%scheme%intprint.GE.2) write(lupri,'(A46,F18.8)') 'Overlap-matrix dot product:',mat_dotproduct(S,S)
CALL LSTIMER('OVERLAP',TS,TE,LUPRI)
call time_II_operations2(JOB_II_get_overlap)
END SUBROUTINE II_get_overlap

!> \brief Calculates overlap integral matrix between 2 different AO basis's
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
SUBROUTINE II_get_mixed_overlap(LUPRI,LUERR,SETTING,S,AO1,AO2,GCAO1,GCAO2)
IMPLICIT NONE
TYPE(MATRIX)          :: S
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,IPRINT,AO1,AO2
LOGICAL               :: GCAO1,GCAO2
!
Integer             :: i,j,LU,nbast2,nbast1

!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR

IPRINT=SETTING%SCHEME%INTPRINT
nbast1 = getNbasis(AO1,ContractedintType,SETTING%MOLECULE(1)%p,LUPRI)
nbast2 = getNbasis(AO2,ContractedintType,SETTING%MOLECULE(2)%p,LUPRI)
call Test_if_64bit_integer_required(nbast1,nbast2)

IF(nBast1.NE.S%nrow)CALL LSQUIT('dim1 mismatch in II_get_mixed_overlap',-1)
IF(nBast2.NE.S%ncol)CALL LSQUIT('dim2 mismatch in II_get_mixed_overlap',-1)
call initIntegralOutputDims(setting%output,nbast1,nbast2,1,1,1)
CALL ls_getIntegrals(AO1,AO2,AOempty,AOempty,&
     &OverlapOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,S,.FALSE.)

IF (GCAO1) call AO2GCAO_half_transform_matrix(S,SETTING,LUPRI,1)
IF (GCAO2) call AO2GCAO_half_transform_matrix(S,SETTING,LUPRI,2)

IF(IPRINT.GT. 1000)THEN
   WRITE(LUPRI,'(A,2X,F16.8)')'Mixed Overlap',mat_dotproduct(S,S)
   call mat_print(S,1,nbast1,1,nbast2,lupri)
ENDIF

END SUBROUTINE II_get_mixed_overlap

!> \brief Calculates overlap integral matrix between 2 different AO basis's
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
SUBROUTINE II_get_mixed_overlap_full(LUPRI,LUERR,SETTING,S,nbast1,nbast2,AO1,AO2)
IMPLICIT NONE
INTEGER               :: LUPRI,LUERR,AO1,AO2,nbast1,nbast2
real(realk)           :: S(nbast1,nbast2)
TYPE(LSSETTING)       :: SETTING
!
INTEGER               :: nbast1TEST,nbast2TEST,IPRINT
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
IPRINT=SETTING%SCHEME%INTPRINT
nbast1TEST = getNbasis(AO1,ContractedintType,SETTING%MOLECULE(1)%p,LUPRI)
nbast2TEST = getNbasis(AO2,ContractedintType,SETTING%MOLECULE(1)%p,LUPRI)
IF(nBast1TEST.NE.nbast1)CALL LSQUIT('dim1 mismatch in II_get_mixed_overlap_full',-1)
IF(nBast2TEST.NE.nbast2)CALL LSQUIT('dim2 mismatch in II_get_mixed_overlap_full',-1)
call initIntegralOutputDims(setting%output,nbast1,nbast2,1,1,1)
CALL ls_getIntegrals(AO1,AO2,AOempty,AOempty,&
     &OverlapOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,S,.FALSE.)
END SUBROUTINE II_get_mixed_overlap_full

!> \brief Calculates one electron fock matrix contribution
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param h the one electron fock matrix contribution
SUBROUTINE II_get_h1(LUPRI,LUERR,SETTING,h)
IMPLICIT NONE
TYPE(MATRIX),target   :: h
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
!
Integer             :: nbast
TYPE(MATRIX),target :: tmp
Real(realk)         :: OLDTHRESH

CALL II_get_nucel_mat(LUPRI,LUERR,SETTING,h)
IF (setting%scheme%intprint.GE.2) THEN
  write(lupri,'(A46,F18.8)') 'Nuclear-electron-attraction-matrix dot product::',mat_dotproduct(h,h)
  IF (setting%scheme%intprint.GE.10) THEN
    write(lupri,'(A46)') 'Nuclear-electron-attraction matrix'
    call mat_print(h,1,h%nrow,1,h%ncol,lupri)
  ENDIF
ENDIF

nbast = h%nrow
CALL mat_init(tmp,nbast,nbast)
CALL II_get_kinetic(LUPRI,LUERR,SETTING,tmp)
IF (setting%scheme%intprint.GE.2) write(lupri,'(A46,F18.8)') 'Kinetic-energy-matrix dot product:',mat_dotproduct(tmp,tmp)

call mat_daxpy(1E0_realk,tmp,h)
CALL mat_free(tmp)

END SUBROUTINE II_get_h1

!> \brief Calculates one electron fock matrix contribution
!> \author S. Reine
!> \date May 2012
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param h the one electron fock matrix contribution
SUBROUTINE II_get_h1_mixed(LUPRI,LUERR,SETTING,h,AO1,AO2)
IMPLICIT NONE
TYPE(MATRIX),target   :: h
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,AO1,AO2
!
TYPE(MATRIX),target :: tmp
Real(realk)         :: OLDTHRESH

CALL II_get_nucel_mat_mixed(LUPRI,LUERR,SETTING,h,AO1,AO2)
!write(lupri,*) 'QQQ New mixed h:',mat_dotproduct(h,h)

CALL mat_init(tmp,h%nrow,h%ncol)
CALL II_get_kinetic_mixed(LUPRI,LUERR,SETTING,tmp,AO1,AO2)
!write(lupri,*) 'QQQ New mixed K:',mat_dotproduct(tmp,tmp)

call mat_daxpy(1E0_realk,tmp,h)
CALL mat_free(tmp)

END SUBROUTINE II_get_h1_mixed

!> \brief Calculates one electron fock matrix contribution
!> \author S. Reine
!> \date May 2012
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param h the one electron fock matrix contribution
SUBROUTINE II_get_h1_mixed_full(LUPRI,LUERR,SETTING,h,nbast1,nbast2,AO1,AO2)
IMPLICIT NONE
INTEGER               :: LUPRI,LUERR,AO1,AO2,nbast1,nbast2
real(realk)           :: h(nbast1,nbast2)
TYPE(LSSETTING)       :: SETTING
!
real(realk),pointer   :: tmp(:,:)
INTEGER               :: nbast1TEST,nbast2TEST
nbast1TEST = getNbasis(AO1,ContractedintType,SETTING%MOLECULE(1)%p,LUPRI)
nbast2TEST = getNbasis(AO2,ContractedintType,SETTING%MOLECULE(1)%p,LUPRI)
IF(nBast1TEST.NE.nbast1)CALL LSQUIT('dim1 mismatch in II_get_h1_mixed_full',-1)
IF(nBast2TEST.NE.nbast2)CALL LSQUIT('dim2 mismatch in II_get_h1_mixed_full',-1)
CALL II_get_nucel_mat_mixed_full(LUPRI,LUERR,SETTING,h,nbast1,nbast2,AO1,AO2)
CALL mem_alloc(tmp,nbast1,nbast2)
CALL II_get_kinetic_mixed_full(LUPRI,LUERR,SETTING,tmp,nbast1,nbast2,AO1,AO2)
call daxpy(nbast1*nbast2,1E0_realk,tmp,1,h,1)
CALL mem_dealloc(tmp)
END SUBROUTINE II_get_h1_mixed_full

!> \brief Calculates the kinetic energy fock matrix contribution
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param T the kinetic energy fock matrix contribution
SUBROUTINE II_get_kinetic(LUPRI,LUERR,SETTING,T)
IMPLICIT NONE
TYPE(MATRIX),target   :: T
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
!
Integer               :: nbast
real(realk)           :: OLDTHRESH
real(realk)         :: TS,TE
call time_II_operations1()
CALL LSTIMER('START ',TS,TE,LUPRI)
CALL mat_zero(T)
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
nbast = T%nrow
call initIntegralOutputDims(setting%output,nbast,1,nbast,1,1)
CALL ls_getIntegrals(AORdefault,AOempty,AORdefault,AOempty,&
     &KineticOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,T,setting%IntegralTransformGC)
CALL LSTIMER('Kinetic',TS,TE,LUPRI)

call time_II_operations2(JOB_II_get_kinetic)
END SUBROUTINE II_get_kinetic

!> \brief Calculates the kinetic energy fock matrix contribution
!> \author S Reine
!> \date May 2012
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param T the kinetic energy fock matrix contribution
SUBROUTINE II_get_kinetic_mixed(LUPRI,LUERR,SETTING,T,AO1,AO3)
IMPLICIT NONE
TYPE(MATRIX),target   :: T
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,AO1,AO3
!
Integer               :: nbast1,nbast3
real(realk)           :: OLDTHRESH
real(realk)         :: TS,TE
call time_II_operations1()
CALL LSTIMER('START ',TS,TE,LUPRI)
CALL mat_zero(T)
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR


nbast1 = getNbasis(AO1,ContractedintType,SETTING%MOLECULE(1)%p,LUPRI)
nbast3 = getNbasis(AO3,ContractedintType,SETTING%MOLECULE(3)%p,LUPRI)

IF(nBast1.NE.T%nrow)CALL LSQUIT('dim1 mismatch in II_get_kinetic_mixed',-1)
IF(nBast3.NE.T%ncol)CALL LSQUIT('dim3 mismatch in II_get_kinetic_mixed',-1)

IF (setting%IntegralTransformGC) call lsquit('Error in II_get_kinetic_mixed. GC transform not implemented',-1)
call initIntegralOutputDims(setting%output,nbast1,1,nbast3,1,1)
CALL ls_getIntegrals(AO1,AOempty,AO3,AOempty,&
     &KineticOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,T,setting%IntegralTransformGC)
CALL LSTIMER('Kinetic',TS,TE,LUPRI)

call time_II_operations2(JOB_II_get_kinetic)
END SUBROUTINE II_get_kinetic_mixed

!> \brief Calculates the kinetic energy fock matrix contribution
!> \author S Reine
!> \date May 2012
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param T the kinetic energy fock matrix contribution
SUBROUTINE II_get_kinetic_mixed_full(LUPRI,LUERR,SETTING,T,nbast1,nbast3,AO1,AO3)
IMPLICIT NONE
INTEGER               :: LUPRI,LUERR,nbast1,nbast3,AO1,AO3
real(realk)           :: T(nbast1,nbast3)
TYPE(LSSETTING)       :: SETTING
!
real(realk)           :: OLDTHRESH
real(realk)         :: TS,TE
call time_II_operations1()
CALL LSTIMER('START ',TS,TE,LUPRI)
CALL ls_dzero(T,nbast1*nbast3)
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
IF (setting%IntegralTransformGC) call lsquit('Error in II_get_kinetic_mixed. GC transform not implemented',-1)
call initIntegralOutputDims(setting%output,nbast1,1,nbast3,1,1)
CALL ls_getIntegrals(AO1,AOempty,AO3,AOempty,&
     &KineticOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,T,setting%IntegralTransformGC)
CALL LSTIMER('Kinetic',TS,TE,LUPRI)
call time_II_operations2(JOB_II_get_kinetic)
END SUBROUTINE II_get_kinetic_mixed_Full

!> \brief Calculates the nuclear attraction fock matrix contribution
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param h the nuclear attraction fock matrix contribution
SUBROUTINE II_get_nucel_mat(LUPRI,LUERR,SETTING,h)
IMPLICIT NONE
TYPE(MATRIX)          :: h
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
!
Integer             :: nbast
TYPE(MATRIX)        :: temp
real(realk)         :: OLDTHRESH
real(realk)         :: TS,TE
call time_II_operations1()
CALL LSTIMER('START ',TS,TE,LUPRI)
CALL mat_zero(h)
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
nbast = h%nrow
! Calculate multipole moments when using FMM
IF(SETTING%SCHEME%FMM) THEN
  call mat_init(temp,nbast,nbast)
  CALL mat_zero(temp)
  CALL ls_attachDmatToSetting(temp,1,setting,'LHS',1,2,.TRUE.,lupri)
  CALL ls_attachDmatToSetting(temp,1,setting,'RHS',3,4,.TRUE.,lupri)
  call ls_multipolemoment(LUPRI,LUERR,SETTING,nbast,0,nbast,nbast,nbast,nbast,&
     &AORdefault,AORdefault,AORdefault,AORdefault,RegularSpec,ContractedInttype,.FALSE.)
  CALL ls_freeDmatFromSetting(setting)
  call mat_free(temp)
ENDIF 
call initIntegralOutputDims(setting%output,nbast,nbast,1,1,1)
CALL ls_getIntegrals(AORdefault,AORdefault,AONuclear,AOempty,&
     &NucpotOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,h,setting%IntegralTransformGC)
CALL LSTIMER('NucElec',TS,TE,LUPRI)
call time_II_operations2(JOB_II_get_nucel_mat)

END SUBROUTINE II_get_nucel_mat

!> \brief Calculates the electrostatic potential integrals
!> \author S. Reine
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param ep the electrostatic potential integrals v_ab = sum_C (ab|C) Z_C
!> \param N the number of nuclei
!> \param R the center of the nuclei
!> \param charge the charge of the nuclei
SUBROUTINE II_get_ep_ab(LUPRI,LUERR,SETTING,ep,N,R,charge)
IMPLICIT NONE
TYPE(MATRIX),target :: ep
INTEGER             :: LUPRI,LUERR,N
REAL(realk)         :: R(3,N) !vector R={x,y,z}
REAL(realk)         :: charge(N)
TYPE(LSSETTING)     :: SETTING
!
type(MOLECULEINFO),pointer :: molecule,Point
real(realk)         :: TS,TE
CALL LSTIMER('START ',TS,TE,LUPRI)

!Build the point charges from the centers and charges
allocate(Point)
call build_pointMolecule(Point,R,N,lupri,charge)
molecule => setting%MOLECULE(1)%p

!Specify that we use the regular molecule for AO 1 and 2, and the point charges for AO 3
call typedef_setMolecules(setting,molecule,1,2,point,3)

!Calculate the integrals and contract with the charges: sum_N (ab|R(N)) charge(N)
call II_get_nucel_mat(LUPRI,LUERR,SETTING,ep)

!Free the point charges
call free_Moleculeinfo(Point)
deallocate(Point)

!Reset molecule 
call typedef_setMolecules(setting,molecule,1,2,3,4)

CALL LSTIMER('ep_ab  ',TS,TE,LUPRI)

END SUBROUTINE II_get_ep_ab

!> \brief Calculates the nuclear attraction fock matrix contribution
!> \author S. Reine
!> \date  May 2012
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param h the nuclear attraction fock matrix contribution
SUBROUTINE II_get_nucel_mat_mixed(LUPRI,LUERR,SETTING,h,AO1,AO2)
IMPLICIT NONE
TYPE(MATRIX)   :: h
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,AO1,AO2
!
Integer             :: nbast1,nbast2
TYPE(MATRIX)        :: temp
real(realk)         :: OLDTHRESH
real(realk)         :: TS,TE
call time_II_operations1()
CALL LSTIMER('START ',TS,TE,LUPRI)
CALL mat_zero(h)
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
nbast1 = getNbasis(AO1,ContractedintType,SETTING%MOLECULE(1)%p,LUPRI)
nbast2 = getNbasis(AO2,ContractedintType,SETTING%MOLECULE(2)%p,LUPRI)

!print *,"nbast1:",nbast1
!print *,"nbast2:",nbast2
!print *,"AO1:",AO1
!print *,"AO2:",AO2
!print *,"h%nrow:",h%nrow
!print *,"h€ncol:",h%ncol

IF(nBast1.NE.h%nrow)CALL LSQUIT('dim1 mismatch in II_get_nucel_mat_mixed',-1)
IF(nBast2.NE.h%ncol)CALL LSQUIT('dim2 mismatch in II_get_nucel_mat_mixed',-1)

IF (setting%IntegralTransformGC) CALL LSQUIT('Error in II_get_nucel_mat_mixed. GC transform not implemented',-1)
! Calculate multipole moments when using FMM
IF(SETTING%SCHEME%FMM) THEN
  CALL LSQUIT('Error in II_get_nucel_mat_mixed. FMM not implemented',-1)
! call mat_init(temp,nbast,nbast)
! CALL mat_zero(temp)
! CALL ls_attachDmatToSetting(temp,1,setting,'LHS',1,2,.TRUE.,lupri)
! CALL ls_attachDmatToSetting(temp,1,setting,'RHS',3,4,.TRUE.,lupri)
! call ls_multipolemoment(LUPRI,LUERR,SETTING,nbast,0,nbast,nbast,nbast,nbast,&
!    &AORdefault,AORdefault,AORdefault,AORdefault,RegularSpec,ContractedInttype,.FALSE.)
! CALL ls_freeDmatFromSetting(setting)
! call mat_free(temp)
ENDIF 
call initIntegralOutputDims(setting%output,nbast1,nbast2,1,1,1)
CALL ls_getIntegrals(AO1,AO2,AONuclear,AOempty,&
     &NucpotOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,h,setting%IntegralTransformGC)
CALL LSTIMER('NucElec',TS,TE,LUPRI)
call time_II_operations2(JOB_II_get_nucel_mat)

END SUBROUTINE II_get_nucel_mat_mixed

!> \brief Calculates the nuclear attraction fock matrix contribution
!> \author S. Reine
!> \date  May 2012
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param h the nuclear attraction fock matrix contribution
SUBROUTINE II_get_nucel_mat_mixed_full(LUPRI,LUERR,SETTING,h,nbast1,nbast2,AO1,AO2)
IMPLICIT NONE
INTEGER               :: LUPRI,LUERR,AO1,AO2,nbast1,nbast2
real(realk)           :: h(nbast1,nbast2)
TYPE(LSSETTING)       :: SETTING
!
real(realk),pointer :: temp(:,:)
real(realk)         :: OLDTHRESH
real(realk)         :: TS,TE
call time_II_operations1()
CALL LSTIMER('START ',TS,TE,LUPRI)
CALL ls_dzero(h,nbast1*nbast2)
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
IF (setting%IntegralTransformGC) CALL LSQUIT('Error in II_get_nucel_mat_mixed. GC transform not implemented',-1)
IF(SETTING%SCHEME%FMM) THEN
  CALL LSQUIT('Error in II_get_nucel_mat_mixed. FMM not implemented',-1)
ENDIF 
call initIntegralOutputDims(setting%output,nbast1,nbast2,1,1,1)
CALL ls_getIntegrals(AO1,AO2,AONuclear,AOempty,&
     &NucpotOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,h,setting%IntegralTransformGC)
CALL LSTIMER('NucElec',TS,TE,LUPRI)
call time_II_operations2(JOB_II_get_nucel_mat)
END SUBROUTINE II_get_nucel_mat_mixed_full

!> \brief Calculates the magnetic derivative of the overlap matrix
!> \author T. Kjaergaard
!> \date 2011-08-22
!> \param Sx The 3 magnetic derivative components of the overlap matrix
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_magderivOverlap(Sx,setting,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Integer,intent(IN)            :: lupri,luerr
Type(matrix),intent(INOUT)    :: Sx(3) !Sx,Sy,Sz
!
integer             :: nbast,I
real(realk)         :: OLDTHRESH
real(realk),parameter ::  D05=0.5E0_realk
call time_II_operations1()
!set threshold
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
nbast = Sx(1)%nrow
call initIntegralOutputDims(setting%output,nbast,nbast,1,1,3)
CALL ls_getIntegrals(AORdefault,AORdefault,AOempty,AOempty,&
     &OverlapOperator,MagDerivSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,Sx,setting%IntegralTransformGC)
DO I=1,3
   call mat_scal(0.5_realk,Sx(I))
ENDDO
call time_II_operations2(JOB_II_get_magderivOverlap)
END SUBROUTINE II_get_magderivOverlap

!> \brief Calculates the magnetic derivative of the overlap matrix with D
!> \author T. Kjaergaard
!> \date 2011-08-22
!> \param Sx The 3 magnetic derivative components of the overlap matrix
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_maggradOverlap(maggrad,D,ndmat,setting,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Real(realk),intent(INOUT)       :: maggrad(3)
Integer,intent(IN)            :: lupri,luerr,ndmat
Type(matrix),intent(IN)      :: D(ndmat)
!
Type(matrix)       :: Dmat_AO(ndmat)
integer             :: I
real(realk)         :: OLDTHRESH
!set threshold 
call time_II_operations1()
IF(ndmat.NE.1)call lsquit('option not verified in II_get_maggradOverlap',-1)
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
IF(setting%IntegralTransformGC)THEN
   DO I=1,ndmat
      CALL mat_init(Dmat_AO(I),D(1)%nrow,D(1)%ncol)
      call GCAO2AO_transform_matrixD2(D(I),Dmat_AO(I),setting,lupri)
   ENDDO
   CALL ls_attachDmatToSetting(Dmat_AO,ndmat,setting,'LHS',1,2,.TRUE.,lupri)
ELSE
   CALL ls_attachDmatToSetting(D,ndmat,setting,'LHS',1,2,.TRUE.,lupri)
ENDIF
call initIntegralOutputDims(setting%output,1,1,1,1,3)
CALL ls_getIntegrals(AORdefault,AORdefault,AOempty,AOempty,&
     &OverlapOperator,MagGradSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL ls_freeDmatFromSetting(setting)
CALL retrieve_Output(lupri,setting,maggrad,setting%IntegralTransformGC)
IF(setting%IntegralTransformGC)THEN
   DO I=1,ndmat
      CALL mat_free(Dmat_AO(I))
   ENDDO
ENDIF
maggrad(1) = - 0.5E0_realk * maggrad(1)
maggrad(2) = - 0.5E0_realk * maggrad(2)
maggrad(3) = - 0.5E0_realk * maggrad(3)
call time_II_operations2(JOB_II_get_maggradOverlap)
END SUBROUTINE II_get_maggradOverlap

!> \brief Calculates the magnetic derivative of the overlap matrix. Differentiated on the ket vector 
!> \author T. Kjaergaard
!> \date 2011-08-22
!> \param Sx The 3 magnetic derivative components of the overlap matrix
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_magderivOverlapR(Sx,setting,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Integer,intent(IN)            :: lupri,luerr
Type(matrix),intent(INOUT)    :: Sx(3) !Sx,Sy,Sz
!
integer             :: nbast,I
real(realk)         :: OLDTHRESH
real(realk),parameter ::  D05=0.5E0_realk
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
nbast = Sx(1)%nrow
call initIntegralOutputDims(setting%output,nbast,nbast,1,1,3)
CALL ls_getIntegrals(AORdefault,AORdefault,AOempty,AOempty,&
     &OverlapOperator,MagDerivRSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,Sx,setting%IntegralTransformGC)
DO I=1,3
   call mat_scal(0.5_realk,Sx(I))
ENDDO
END SUBROUTINE II_get_magderivOverlapR

!> \brief Calculates the magnetic derivative of the overlap matrix. Differentiated on the bra vector (LHS) 
!> \author T. Kjaergaard
!> \date 2011-08-22
!> \param Sx The 3 magnetic derivative components of the overlap matrix
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_magderivOverlapL(Sx,setting,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Integer,intent(IN)            :: lupri,luerr
Type(matrix),intent(INOUT)    :: Sx(3) !Sx,Sy,Sz
!
integer             :: nbast,I
real(realk)         :: OLDTHRESH
real(realk),parameter ::  D05=0.5E0_realk
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
nbast = Sx(1)%nrow
call initIntegralOutputDims(setting%output,nbast,nbast,1,1,3)
CALL ls_getIntegrals(AORdefault,AORdefault,AOempty,AOempty,&
     &OverlapOperator,MagDerivLSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,Sx,setting%IntegralTransformGC)
DO I=1,3
   call mat_scal(0.5_realk,Sx(I))
ENDDO
END SUBROUTINE II_get_magderivOverlapL

!> \brief Calculates the nuclear attraction fock matrix contribution
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param h the nuclear attraction fock matrix contribution
SUBROUTINE II_get_ep_integrals(LUPRI,LUERR,SETTING,Integral,R)
IMPLICIT NONE
TYPE(MATRIX)        :: Integral
REAL(realk)         :: R(3,1) !vector R={x,y,z}
TYPE(LSSETTING)     :: SETTING
INTEGER             :: LUPRI,LUERR
!
Integer             :: nbast
real(realk)         :: OLDTHRESH
type(MOLECULE_PT)   :: temp,Point

CALL mat_zero(Integral)
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
nbast = Integral%nrow
call initIntegralOutputDims(setting%output,nbast,nbast,1,1,1)
allocate(Point%p)
call build_pointMolecule(Point%p,R,1,lupri)
temp%p  => setting%MOLECULE(3)%p
setting%MOLECULE(3)%p => Point%p
CALL ls_getIntegrals(AORdefault,AORdefault,AONuclear,AOempty,&
     &NucpotOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,Integral,setting%IntegralTransformGC)
call free_Moleculeinfo(Point%p)
setting%MOLECULE(3)%p => temp%p
END SUBROUTINE II_get_ep_integrals

!> \brief Calculates the nuclear attraction fock matrix contribution
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param h the nuclear attraction fock matrix contribution
SUBROUTINE II_get_ep_integrals2(LUPRI,LUERR,SETTING,R,D,output)
IMPLICIT NONE
TYPE(MATRIX),target :: D
!TYPE(MATRIX)        :: Integral
REAL(realk)         :: R(3,1) !vector R={x,y,z}
REAL(realk)         :: output(1)
TYPE(LSSETTING)     :: SETTING
INTEGER             :: LUPRI,LUERR
!
TYPE(MATRIX)        :: D_AO
Integer             :: nbast
real(realk)         :: OLDTHRESH
type(MOLECULE_PT)   :: temp,Point
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
nbast = D%nrow
call initIntegralOutputDims(setting%output,1,1,1,1,1)
allocate(Point%p)
call build_pointMolecule(Point%p,R,1,lupri)
temp%p  => setting%MOLECULE(1)%p
setting%MOLECULE(1)%p => Point%p
IF(setting%IntegralTransformGC)THEN
   CALL mat_init(D_AO,D%nrow,D%ncol)
   call GCAO2AO_transform_matrixD2(D,D_AO,setting,lupri)
   CALL ls_attachDmatToSetting(D_AO,1,setting,'RHS',3,4,.TRUE.,lupri)
ELSE
   CALL ls_attachDmatToSetting(D,1,setting,'RHS',3,4,.TRUE.,lupri)
ENDIF
CALL ls_jengine(AONuclear,AOempty,AORdefault,AORdefault,&
     &          NucpotOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
!no fmm contribution
CALL ls_freeDmatFromSetting(setting)
IF(setting%IntegralTransformGC)THEN
   CALL mat_free(D_AO)
ENDIF
CALL retrieve_Output(lupri,setting,output,setting%IntegralTransformGC)
call free_Moleculeinfo(Point%p)
setting%MOLECULE(1)%p => temp%p
!call print_mol(setting%MOLECULE(3)%p,lupri)
END SUBROUTINE II_get_ep_integrals2

!> \brief Calculates the nuclear attraction fock matrix contribution
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param h the nuclear attraction fock matrix contribution
SUBROUTINE II_get_ep_integrals3(LUPRI,LUERR,SETTING,R,N,D,output)
IMPLICIT NONE
TYPE(MATRIX),target :: D
!TYPE(MATRIX)        :: Integral
INTEGER             :: LUPRI,LUERR,N
REAL(realk)         :: R(3,N) !vector R={x,y,z}
REAL(realk)         :: output(N)
TYPE(LSSETTING)     :: SETTING
!
TYPE(matrix)        :: D_AO
Integer             :: nbast
real(realk)         :: OLDTHRESH
type(MOLECULE_PT)   :: temp,Point
real(realk)         :: TS,TE
CALL LSTIMER('START ',TS,TE,LUPRI)
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
nbast = D%nrow!Integral%nrow
call initIntegralOutputDims(setting%output,N,1,1,1,1)
!critical important that ls_attachDmatToSetting is called before 
!change of molecule, because GCAO transform requires full molecule
IF(setting%IntegralTransformGC)THEN
   CALL mat_init(D_AO,D%nrow,D%ncol)
   call GCAO2AO_transform_matrixD2(D,D_AO,setting,lupri)
   CALL ls_attachDmatToSetting(D_AO,1,setting,'RHS',3,4,.TRUE.,lupri)
ELSE
   CALL ls_attachDmatToSetting(D,1,setting,'RHS',3,4,.TRUE.,lupri)
ENDIF
allocate(Point%p)
call build_pointMolecule(Point%p,R,N,lupri)
temp%p  => setting%MOLECULE(1)%p
setting%MOLECULE(1)%p => Point%p
!call print_mol(setting%MOLECULE(3)%p,lupri)
CALL ls_jengine(AOpCharge,AOempty,AORdefault,AORdefault,&
     &          NucpotOperator,pChargeSpec,ContractedInttype,SETTING,LUPRI,LUERR)
!no fmm contribution
CALL ls_freeDmatFromSetting(setting)
IF(setting%IntegralTransformGC)THEN
   CALL mat_free(D_AO)
ENDIF
call free_Moleculeinfo(Point%p)
setting%MOLECULE(1)%p => temp%p
CALL retrieve_Output(lupri,setting,output,setting%IntegralTransformGC)
CALL LSTIMER('ElPot ',TS,TE,LUPRI)
!call print_mol(setting%MOLECULE(3)%p,lupri)
END SUBROUTINE II_get_ep_integrals3

!> \brief Calculates the full molecular gradient
!> \author T. Kjaergaard
!> \date 2010-04-26
!> \param lupri Default print unit
!> \param F The Fock/Kohn-Sham matrix
!> \param D density matrix
!> \param setting Integral evalualtion settings
!> \param dodft is it a DFT or HF calc 
!> \param doprint should we print the gradient
subroutine II_get_molecular_gradient(GRAD,lupri,F,D,setting,dodft,doprint)
  implicit none
  integer,intent(in)        :: lupri
  type(Matrix),intent(in)   :: F,D
  TYPE(LSSETTING)           :: SETTING
  logical                   :: dodft,doprint
!
  type(matrixp)               :: Dmat(1)
  type(Matrix),target         :: tempm1,tempm3
  type(Matrix)                :: tempm2
  integer   :: nbast,natom,ndmat,ix,iatom,luerr
  real(realk),pointer  :: tmpGRAD(:,:)
  real(realk), intent(inout) :: GRAD(3,SETTING%MOLECULE(1)%p%Natoms)
  real(realk) :: nrm,ts,te
  call lstimer('START ',ts,te,lupri)

 if(doprint) then
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,'(8X,A)') '*************************************************************'
    write(lupri,'(8X,A)') '*            MOLECULAR GRADIENT RESULTS (in a.u.)           *'
    write(lupri,'(8X,A)') '*************************************************************'
    write(lupri,*) 
    write(lupri,*) 
 end if

  nbast = D%nrow    
  natom = setting%molecule(1)%p%natoms
  ndmat = 1
  luerr = 0
  call mat_init(tempm1,nbast,nbast)
  call mat_assign(tempm1,D)
  call II_get_nn_gradient(Grad,setting,lupri,luerr)
  if(doprint)CALL LS_PRINT_GRADIENT(lupri,setting%molecule(1)%p,GRAD,natom,'nuclear rep')
  call lstimer('nn-grad ',ts,te,lupri)
  Dmat(1)%p => tempm1
  call mem_alloc(tmpGrad,3,nAtom)
  CALL II_get_twoElectron_gradient(tmpGrad,natom,Dmat,Dmat,ndmat,ndmat,setting,lupri,luerr)
  CALL DSCAL(3*nAtom,4E0_realk,tmpGrad,1)
  if(doprint)CALL LS_PRINT_GRADIENT(lupri,setting%molecule(1)%p,tmpGRAD,natom,'twoElectron')
  call lstimer('two-grad',ts,te,lupri)
  CALL DAXPY(3*natom,1E0_realk,tmpGrad,1,Grad,1)
  CALL II_get_ne_gradient(tmpGrad,Dmat,ndmat,setting,lupri,luerr)
  CALL DSCAL(3*nAtom,2E0_realk,tmpGrad,1)
  if(doprint)CALL LS_PRINT_GRADIENT(lupri,setting%molecule(1)%p,tmpGRAD,natom,'nuc-el attr')
  call lstimer('ne-grad ',ts,te,lupri)
  CALL DAXPY(3*natom,1E0_realk,tmpGrad,1,Grad,1)
  if (dodft) THEN
     CALL II_get_xc_geoderiv_molgrad(lupri,luerr,setting,nbast,D,tmpGrad,natom)
     if(doprint)CALL LS_PRINT_GRADIENT(lupri,setting%molecule(1)%p,tmpGRAD,natom,'exchange-corr')
     CALL DAXPY(3*natom,1E0_realk,tmpGrad,1,Grad,1)
     call lstimer('xc-grad ',ts,te,lupri)
  endif
  CALL II_get_kinetic_gradient(tmpGrad,Dmat,ndmat,setting,lupri,luerr)
  CALL DSCAL(3*nAtom,2E0_realk,tmpGrad,1)
  if(doprint)CALL LS_PRINT_GRADIENT(lupri,setting%molecule(1)%p,tmpGRAD,natom,'kinetic')
  CALL DAXPY(3*natom,1E0_realk,tmpGrad,1,Grad,1)
  call lstimer('kin-grad',ts,te,lupri)
  call mat_init(tempm2,nbast,nbast)
  call mat_init(tempm3,nbast,nbast)
  call mat_mul(tempm1,F,'n','n',1E0_realk,0E0_realk,tempm2)
  call mat_mul(tempm2,D,'n','n',-1E0_realk,0E0_realk,tempm3)
  call mat_free(tempm1)
  call mat_free(tempm2)

  Dmat(1)%p => tempm3 !the DFD mat
  CALL II_get_reorthoNormalization(tmpGrad,Dmat,ndmat,setting,lupri,luerr)
  call mat_free(tempm3)
  CALL DSCAL(3*nAtom,2E0_realk,tmpGrad,1)
  if(doprint)CALL LS_PRINT_GRADIENT(lupri,setting%molecule(1)%p,tmpGRAD,natom,'reorthonomal')
  call lstimer('reorth-g',ts,te,lupri)
  CALL DAXPY(3*natom,1E0_realk,tmpGrad,1,Grad,1)
  CALL LS_PRINT_GRADIENT(lupri,setting%molecule(1)%p,GRAD,natom,'TOTAL')

 if(doprint) then
    nrm = 0E0_realk
    DO iatom = 1,natom
      DO ix=1,3
        nrm = nrm + Grad(ix,iatom)*Grad(ix,iatom)
      ENDDO
    ENDDO
    nrm = sqrt(nrm/3/natom)

    write(lupri,*) 
    write(lupri,*)
    write(lupri,'(1X,A24,F23.16)') 'RMS gradient norm (au): ',nrm
    write(lupri,*) 
    write(lupri,*) 
 end if
  call mem_dealloc(tmpGrad)

end subroutine II_GET_MOLECULAR_GRADIENT

!> \brief Calculates the two-electron contribution to the molecular gradient
!> \author S. Reine
!> \date 2010-03-01
!> \param twoGrad The two-electron gradient
!> \param natoms The number of atoms
!> \param DmatLHS The left-hand-side density matrix
!> \param DmatRHS The reft-hand-side density matrix
!> \param ndlhs The number of LHS density matrices
!> \param ndrhs The number of RHS density matrices
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_twoElectron_gradient(twoGrad,natoms,DmatLHS,DmatRHS,ndlhs,ndrhs,setting,lupri,luerr)
IMPLICIT NONE
Integer,intent(IN)            :: lupri,luerr,ndlhs,ndrhs,natoms
TYPE(LSSETTING) :: SETTING
Real(realk),intent(INOUT)     :: twoGrad(3,natoms)
Type(matrixp),intent(IN)      :: DmatLHS(ndlhs),DmatRHS(ndrhs)
!
Real(realk),pointer :: exchangeGrad(:,:)
!nAtoms = setting%molecule(1)%p%nAtoms
CALL mem_alloc(exchangeGrad,3,nAtoms)

CALL II_get_J_gradient(twoGrad,DmatLHS,DmatRHS,ndlhs,ndrhs,setting,lupri,luerr)
CALL DSCAL(3*nAtoms,4.0E0_realk,twoGRAD,1)
CALL LS_PRINT_GRADIENT(lupri,setting%molecule(1)%p,twoGRAD,natoms,'Coulomb')

exchangeGrad = 0E0_realk

CALL II_get_K_gradient(exchangeGrad,DmatLHS,DmatRHS,ndlhs,ndrhs,setting,lupri,luerr)
CALL DSCAL(3*nAtoms,4.0E0_realk,exchangeGrad,1)
CALL LS_PRINT_GRADIENT(lupri,setting%molecule(1)%p,exchangeGRAD,natoms,'exchange')

CALL DAXPY(3*nAtoms,1E0_realk,exchangeGrad,1,twoGrad,1)
CALL DSCAL(3*nAtoms,0.25E0_realk,twoGRAD,1)
CALL mem_dealloc(exchangeGrad)

END SUBROUTINE II_get_twoElectron_gradient

!> \brief Calculates the exchange contribution to the molecular gradient
!> \author S. Reine
!> \date 2010-03-01
!> \param kGrad The exchange gradient
!> \param DmatLHS The left-hand-side (or first electron) density matrix
!> \param DmatRHS The reft-hand-side (or second electron) density matrix
!> \param ndlhs The number of LHS density matrices
!> \param ndrhs The number of RHS density matrices
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_K_gradient(kGrad,DmatLHS,DmatRHS,ndlhs,ndrhs,setting,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING)               :: setting
Integer,intent(IN)            :: ndlhs,ndrhs,lupri,luerr
Real(realk),intent(INOUT)     :: kGrad(3,setting%molecule(1)%p%nAtoms)
Type(matrixp),intent(IN)      :: DmatLHS(ndlhs),DmatRHS(ndrhs)
!
logical                       :: ADMMexchange
!
call time_II_operations1()


! *********************************************************************************
! *                      ADMM exchange
! *********************************************************************************
ADMMexchange = setting%scheme%ADMM_EXCHANGE
IF (ADMMexchange) THEN
   !FixMe Should also work for incremental scheme
   IF(incremental_scheme)THEN
      call lsquit('Auxiliary Density Matrix Calculation requires NOINCREM',-1)
   ENDIF
   CALL II_get_ADMM_K_gradient(kGrad,DmatLHS,DmatRHS,ndlhs,ndrhs,setting,lupri,luerr)
 ELSE
! *********************************************************************************
! *                      regular exchange
! *********************************************************************************
   CALL II_get_regular_K_gradient(kGrad,DmatLHS,DmatRHS,ndlhs,ndrhs,setting,lupri,luerr)
ENDIF
call time_II_operations2(JOB_II_GET_K_GRADIENT)
END SUBROUTINE II_get_K_gradient


!> \brief Calculates the regular exchange contribution to the molecular gradient
!> \author S. Reine
!> \date 2010-03-01
!> \param kGrad The exchange gradient
!> \param DmatLHS The left-hand-side (or first electron) density matrix
!> \param DmatRHS The reft-hand-side (or second electron) density matrix
!> \param ndlhs The number of LHS density matrices
!> \param ndrhs The number of RHS density matrices
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_regular_K_gradient(kGrad,DmatLHS,DmatRHS,ndlhs,ndrhs,setting,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING) :: SETTING
Real(realk),intent(INOUT)       :: kGrad(3,setting%molecule(1)%p%nAtoms)
Integer,intent(IN)            :: lupri,luerr,ndlhs,ndrhs
Type(matrixp),intent(IN)      :: DmatLHS(ndlhs),DmatRHS(ndrhs)
!
integer               :: nAtoms,iDmat,I,J
real(realk)           :: Factor, OLDTHRESH
integer               :: Oper
Integer                :: symLHS(ndlhs),symRHS(ndrhs)
Integer                :: nlhs,nrhs,Dascreen_thrlog
Type(matrixp)          :: DLHS(2*ndlhs),DRHS(2*ndrhs)
Type(matrix)           :: DLHS_AO(2*ndlhs),DRHS_AO(2*ndrhs)
logical                :: Dalink
real(realk),pointer :: grad(:,:)
! Check first if the exchange-contribution should be calculated. If not exit 
! this subroutine
IF (SETTING%SCHEME%exchangeFactor.EQ. 0.0E0_realk) RETURN

! we screenin based on the non differentiated integrals
! so we loosen the screening threshold with a factor 10
! and use DaLink to speed up the calculation. 
OLDTHRESH = SETTING%SCHEME%CS_THRESHOLD
SETTING%SCHEME%CS_THRESHOLD = SETTING%SCHEME%CS_THRESHOLD*1.0E-1_realk
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%K_THR

! Check symetry. Split non-symmetric matrices to symmetric and anti symmetric parts. 
! Symmetric, anti-symmetric and anti-symmetric, symmetric paris vanish.
call II_split_dmats(DmatLHS,DmatRHS,ndlhs,ndrhs,DLHS,DRHS,symLHS,symRHS,nlhs,nrhs)


IF (nlhs.EQ. 0) RETURN
! Actual calculation

nAtoms = setting%molecule(1)%p%nAtoms
IF(setting%IntegralTransformGC)THEN
   DO I=1,nlhs
      CALL mat_init(DLHS_AO(I),DLHS(1)%p%nrow,DLHS(1)%p%ncol)
      call GCAO2AO_transform_matrixD2(DLHS(I)%p,DLHS_AO(I),setting,lupri)
   ENDDO
   DO I=1,nrhs
      CALL mat_init(DRHS_AO(I),DRHS(1)%p%nrow,DRHS(1)%p%ncol)
      call GCAO2AO_transform_matrixD2(DRHS(I)%p,DRHS_AO(I),setting,lupri)
   ENDDO
   CALL ls_attachDmatToSetting(DLHS_AO,nlhs,setting,'LHS',1,3,.TRUE.,lupri)
   CALL ls_attachDmatToSetting(DRHS_AO,nrhs,setting,'RHS',2,4,.TRUE.,lupri)
ELSE
   CALL ls_attachDmatToSetting(DLHS,nlhs,setting,'LHS',1,3,.TRUE.,lupri)
   CALL ls_attachDmatToSetting(DRHS,nrhs,setting,'RHS',2,4,.TRUE.,lupri)
ENDIF

IF (SETTING%SCHEME%CAM) THEN
  Oper = CAMOperator       !Coulomb attenuated method
ELSEIF (SETTING%SCHEME%SR_EXCHANGE) THEN
  Oper = ErfcOperator      !Short-Range Coulomb screened exchange
ELSE
  Oper = CoulombOperator   !Regular Coulomb metric 
ENDIF

Dalink = setting%scheme%daLinK
Dascreen_thrlog = setting%scheme%daScreen_THRLOG
setting%scheme%daLinK = .TRUE.
setting%scheme%daScreen_THRLOG = 0
!Calculates the HF-exchange contribution

call initIntegralOutputDims(setting%Output,3,nAtoms,1,1,1)
call ls_get_exchange_mat(AORdefault,AORdefault,AORdefault,AORdefault,&
     &                   Oper,GradientSpec,ContractedInttype,SETTING,LUPRI,LUERR)
#ifdef VAR_MPI
! Hack - since symmetry is currently turned off in exchange the resultTensor 
!        currently has dimension 3,nAtom*4
  call mem_alloc(grad,3,nAtoms*4)
  CALL retrieve_Output(lupri,setting,grad,.FALSE.)
  Kgrad(:,1:nAtoms) = grad(:,1:nAtoms)
  call mem_dealloc(grad)
#else
CALL retrieve_Output(lupri,setting,kGrad,.FALSE.)
#endif

Factor = SETTING%SCHEME%exchangeFactor
do I=1,3
   do J=1,natoms
      kgrad(I,J)=Factor*kGrad(I,J)
   enddo
enddo

CALL ls_freeDmatFromSetting(setting)
IF(setting%IntegralTransformGC)THEN
   DO I=1,nlhs
      CALL mat_free(DLHS_AO(I))
   ENDDO
   DO I=1,nrhs
      CALL mat_free(DRHS_AO(I))
   ENDDO
ENDIF

! Free allocated memory
call II_free_split_dmats(DLHS,DRHS,symLHS,symRHS,ndlhs,ndrhs)
setting%scheme%daLinK=Dalink
setting%scheme%daScreen_THRLOG = Dascreen_thrlog

SETTING%SCHEME%CS_THRESHOLD = OLDTHRESH
END SUBROUTINE II_get_regular_K_gradient

!> \brief Calculates the exchange contribution to the molecular gradient
!> \author S. Reine
!> \date 2010-03-01
!> \param kGrad The exchange gradient
!> \param DmatLHS The left-hand-side (or first electron) density matrix
!> \param DmatRHS The reft-hand-side (or second electron) density matrix
!> \param ndlhs The number of LHS density matrices
!> \param ndrhs The number of RHS density matrices
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_K_gradientfull(kGrad,DmatLHS,DmatRHS,nbast,ndlhs,ndrhs,setting,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING)               :: setting
Integer,intent(IN)            :: ndlhs,ndrhs,lupri,luerr,nbast
Real(realk),intent(INOUT)     :: kGrad(3,setting%molecule(1)%p%nAtoms)
Real(realk),intent(IN)         :: DmatLHS(nbast,nbast,ndlhs)
Real(realk),intent(IN)         :: DmatRHS(nbast,nbast,ndrhs)
!
logical                       :: ADMMexchange
call time_II_operations1()
ADMMexchange = setting%scheme%ADMM_EXCHANGE
IF (ADMMexchange) THEN
   WRITE(*,*)     "The ADMM approximation isn't implemented for unrestricted cases yet."
   WRITE(LUPRI,*) "The ADMM approximation isn't implemented for unrestricted cases yet."
   call lsquit('ADMM exchange not implemented for full Dmat.',-1)
ELSE
   CALL II_get_regular_K_gradientfull(kGrad,DmatLHS,DmatRHS,nbast,ndlhs,ndrhs,setting,lupri,luerr)
ENDIF
call time_II_operations2(JOB_II_GET_K_GRADIENT)
END SUBROUTINE II_get_K_gradientfull

!> \brief Calculates the regular exchange contribution to the molecular gradient
!> \author S. Reine
!> \date 2010-03-01
!> \param kGrad The exchange gradient
!> \param DmatLHS The left-hand-side (or first electron) density matrix
!> \param DmatRHS The reft-hand-side (or second electron) density matrix
!> \param ndlhs The number of LHS density matrices
!> \param ndrhs The number of RHS density matrices
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_regular_K_gradientfull(kGrad,DmatLHS,DmatRHS,nbast,ndlhs,ndrhs,setting,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING) :: SETTING
Real(realk),intent(INOUT)       :: kGrad(3,setting%molecule(1)%p%nAtoms)
Integer,intent(IN)            :: lupri,luerr,ndlhs,ndrhs,nbast
Real(realk),intent(IN)         :: DmatLHS(nbast,nbast,ndlhs)
Real(realk),intent(IN)         :: DmatRHS(nbast,nbast,ndrhs)
integer :: Oper
IF (SETTING%SCHEME%CAM) THEN
  Oper = CAMOperator       !Coulomb attenuated method
ELSEIF (SETTING%SCHEME%SR_EXCHANGE) THEN
  Oper = ErfcOperator      !Short-Range Coulomb screened exchange
ELSE
  Oper = CoulombOperator   !Regular Coulomb metric 
ENDIF
call II_get_regular_K_gradientmixedfull(kGrad,DmatLHS,DmatRHS,nbast,ndlhs,ndrhs,setting,&
     & AORdefault,AORdefault,AORdefault,AORdefault,Oper,lupri,luerr)

end SUBROUTINE II_get_regular_K_gradientfull


!> \brief Calculates the regular exchange contribution to the molecular gradient
!> \author S. Reine
!> \date 2010-03-01
!> \param kGrad The exchange gradient
!> \param DmatLHS The left-hand-side (or first electron) density matrix
!> \param DmatRHS The reft-hand-side (or second electron) density matrix
!> \param ndlhs The number of LHS density matrices
!> \param ndrhs The number of RHS density matrices
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_regular_K_gradientmixedfull(kGrad,DmatLHS,DmatRHS,nbast,ndlhs,ndrhs,setting,&
     & AO1,AO2,AO3,AO4,Oper,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING) :: SETTING
Real(realk),intent(INOUT)       :: kGrad(3,setting%molecule(1)%p%nAtoms)
Integer,intent(IN)            :: lupri,luerr,ndlhs,ndrhs,AO1,AO2,AO3,AO4,Oper,nbast
Real(realk),intent(IN)         :: DmatLHS(nbast,nbast,ndlhs)
Real(realk),intent(IN)         :: DmatRHS(nbast,nbast,ndrhs)
!Type(matrixp),intent(IN)      :: DmatLHS(ndlhs),DmatRHS(ndrhs)
!
integer               :: nAtoms,iDmat,I,J
real(realk)           :: Factor, OLDTHRESH
Integer                :: symLHS(ndlhs),symRHS(ndrhs)
Integer                :: nlhs,nrhs,Dascreen_thrlog
!Type(matrixp)          :: DLHS(2*ndlhs),DRHS(2*ndrhs)
!Type(matrix)           :: DLHS_AO(2*ndlhs),DRHS_AO(2*ndrhs)
logical                :: Dalink,IntegralTransformGC
Real(realk),pointer   :: DRHS(:,:,:),DRHStmp(:,:,:)
Real(realk),pointer   :: DLHS(:,:,:),DLHStmp(:,:,:)
real(realk),pointer :: grad(:,:)
! Check first if the exchange-contribution should be calculated. If not exit 
! this subroutine
IF (SETTING%SCHEME%exchangeFactor.EQ. 0.0E0_realk) RETURN

! we screenin based on the non differentiated integrals
! so we loosen the screening threshold with a factor 10
! and use DaLink to speed up the calculation. 
OLDTHRESH = SETTING%SCHEME%CS_THRESHOLD
SETTING%SCHEME%CS_THRESHOLD = SETTING%SCHEME%CS_THRESHOLD*1.0E-1_realk
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%K_THR

! Check symetry. Split non-symmetric matrices to symmetric and anti symmetric parts. 
! Symmetric, anti-symmetric and anti-symmetric, symmetric paris vanish.
call II_split_dmatsFULL(DmatLHS,DmatRHS,nbast,ndlhs,ndrhs,DLHStmp,DRHStmp,symLHS,symRHS,nlhs,nrhs)

IF (nlhs.EQ. 0) RETURN

! Actual calculation

nAtoms = setting%molecule(1)%p%nAtoms
IF(setting%IntegralTransformGC)THEN
   call mem_alloc(DRHS,nbast,nbast,nrhs)
   call GCAO2AO_transform_fullD(DRHStmp,DRHS,nbast,nrhs,setting,lupri)
   call mem_alloc(DLHS,nbast,nbast,nlhs)
   call GCAO2AO_transform_fullD(DLHStmp,DLHS,nbast,nlhs,setting,lupri)
   IntegralTransformGC = .TRUE.
   CALL ls_attachDmatToSetting(DLHS,nbast,nbast,nlhs,setting,'LHS',1,3,lupri)
   CALL ls_attachDmatToSetting(DRHS,nbast,nbast,nrhs,setting,'RHS',2,4,lupri)
ELSE
   CALL ls_attachDmatToSetting(DLHStmp,nbast,nbast,nlhs,setting,'LHS',1,3,lupri)
   CALL ls_attachDmatToSetting(DRHStmp,nbast,nbast,nrhs,setting,'RHS',2,4,lupri)
   IntegralTransformGC=.FALSE.
ENDIF
setting%IntegralTransformGC = .FALSE.

Dalink = setting%scheme%daLinK
Dascreen_thrlog = setting%scheme%daScreen_THRLOG
setting%scheme%daLinK = .TRUE.
setting%scheme%daScreen_THRLOG = 0
!Calculates the HF-exchange contribution

call initIntegralOutputDims(setting%Output,3,nAtoms,1,1,1)
call ls_get_exchange_mat(AORdefault,AORdefault,AORdefault,AORdefault,&
     &                   Oper,GradientSpec,ContractedInttype,SETTING,LUPRI,LUERR)
#ifdef VAR_MPI
! Hack - since symmetry is currently turned off in exchange the resultTensor 
!        currently has dimension 3,nAtom*4
  call mem_alloc(grad,3,nAtoms*4)
  CALL retrieve_Output(lupri,setting,grad,.FALSE.)
  Kgrad(:,1:nAtoms) = grad(:,1:nAtoms)
  call mem_dealloc(grad)
#else
CALL retrieve_Output(lupri,setting,kGrad,.FALSE.)
#endif

Factor = SETTING%SCHEME%exchangeFactor
do I=1,3
   do J=1,natoms
      kgrad(I,J)=Factor*kGrad(I,J)
   enddo
enddo

CALL ls_freeDmatFromSetting(setting)
IF(IntegralTransformGC)THEN
   call mem_dealloc(DRHS)
   call mem_dealloc(DLHS)
ENDIF
setting%IntegralTransformGC = IntegralTransformGC

! Free allocated memory
call II_free_split_dmatsFull(DLHStmp,DRHStmp)
setting%scheme%daLinK=Dalink
setting%scheme%daScreen_THRLOG = Dascreen_thrlog

SETTING%SCHEME%CS_THRESHOLD = OLDTHRESH
END SUBROUTINE II_get_regular_K_gradientmixedfull

!> \brief Driver for calculating the electron-electron repulsion contribution to the molecular gradient
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
SUBROUTINE II_get_J_gradient(eeGrad,DmatLHS,DmatRHS,ndlhs,ndrhs,setting,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING) :: SETTING
Real(realk),intent(INOUT)     :: eeGrad(3,setting%molecule(1)%p%nAtoms)
Integer,intent(IN)            :: lupri,luerr,ndlhs,ndrhs
Type(matrixp),intent(IN)      :: DmatLHS(ndlhs),DmatRHS(ndrhs)
!
Integer                :: symLHS(ndlhs),symRHS(ndrhs)
Integer                :: nlhs,nrhs,I
Type(matrixp)          :: DLHS(ndlhs),DRHS(ndrhs)
Type(matrixp)          :: DLHS_AO(ndlhs),DRHS_AO(ndrhs)
! Check symetry and remove all anti-symmetric components (which will not contribute) 
! to the Coulomb ERI gradient
call time_II_operations1()
call II_symmetrize_dmats(DmatLHS,DmatRHS,ndlhs,ndrhs,DLHS,DRHS,symLHS,symRHS,nlhs,nrhs)

IF (nlhs.EQ. 0) RETURN

IF(setting%IntegralTransformGC)THEN
   DO I=1,nlhs
      allocate(DLHS_AO(I)%p)
      CALL mat_init(DLHS_AO(I)%p,DLHS(1)%p%nrow,DLHS(1)%p%ncol)
      call GCAO2AO_transform_matrixD2(DLHS(I)%p,DLHS_AO(I)%p,setting,lupri)
   ENDDO
   DO I=1,nrhs
      allocate(DRHS_AO(I)%p)
      CALL mat_init(DRHS_AO(I)%p,DRHS(1)%p%nrow,DRHS(1)%p%ncol)
      call GCAO2AO_transform_matrixD2(DRHS(I)%p,DRHS_AO(I)%p,setting,lupri)
   ENDDO
ELSE
   DO I=1,nlhs
      DLHS_AO(I)%p => DLHS(I)%p
   ENDDO
   DO I=1,nrhs
      DRHS_AO(I)%p => DRHS(I)%p
   ENDDO
ENDIF

! Actual calculation
IF (setting%scheme%densfit) THEN
  CALL II_get_df_J_gradient(eeGrad,DLHS_AO(1:nlhs),DRHS_AO(1:nrhs),nlhs,nrhs,setting,lupri,luerr)
ELSE
  CALL II_get_J_gradient_regular(eeGrad,DLHS_AO(1:nlhs),DRHS_AO(1:nrhs),nlhs,nrhs,setting,lupri,luerr)
ENDIF

! Free allocated memory
call II_free_symmetrized_dmats(DLHS,DRHS,symLHS,symRHS,ndlhs,ndrhs)
IF(setting%IntegralTransformGC)THEN
   DO I=1,nlhs
      CALL mat_free(DLHS_AO(I)%p)
      deallocate(DLHS_AO(I)%p)
   ENDDO
   DO I=1,nrhs
      CALL mat_free(DRHS_AO(I)%p)
      deallocate(DRHS_AO(I)%p)
   ENDDO
ENDIF
call time_II_operations2(JOB_II_get_J_gradient)

END SUBROUTINE II_get_J_gradient

!> \brief Calculates the (regular) electron-electron repulsion contribution to the molecular gradient
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
SUBROUTINE II_get_J_gradient_regular(eeGrad,DmatLHS,DmatRHS,ndlhs,ndrhs,setting,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Real(realk),intent(INOUT)     :: eeGrad(3,setting%molecule(1)%p%nAtoms)
Integer,intent(IN)            :: lupri,luerr,ndlhs,ndrhs
Type(matrixp),intent(IN)      :: DmatLHS(ndlhs),DmatRHS(ndrhs)
!
integer                   :: nAtoms,iDmat,nlhs,nrhs
type(matrixp)             :: eeGradMat(1)
type(matrix),target       :: eeGradTarget
Type(matrixp),pointer     :: DLHS(:),DRHS(:)
logical                   :: same, save_screen,saveNOSEGMENT
logical :: ReCalcGab

integer :: nbast

!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%J_THR
nbast = DmatLHS(1)%p%nrow

same = ls_same_mats(DmatLHS,DmatRHS,ndlhs,ndrhs)
NULLIFY(DLHS)
NULLIFY(DRHS)
IF (.NOT.same.AND.ndlhs.EQ.ndrhs) THEN
! If the density-matrices are not the same we make a dual set of matrices -
! for two matrices DLHS = D and DRHS = P we have the two contributions
! D_ab(ab'|cd)P_cd and P_ab(ab'|cd)D_cd  (since we only differentiate on the LHS)
  nlhs = 2*ndlhs
  nrhs = 2*ndrhs
  call mem_alloc(DLHS,nlhs)
  call mem_alloc(DRHS,nrhs)
  DO idmat = 1,ndrhs
    DLHS(idmat)%p       => DmatLHS(idmat)%p
    DLHS(ndlhs+idmat)%p => DmatRHS(idmat)%p
    DRHS(idmat)%p       => DmatRHS(idmat)%p
    DRHS(ndlhs+idmat)%p => DmatLHS(idmat)%p
  ENDDO
ELSE
  nlhs = ndlhs
  nrhs = ndrhs
  call mem_alloc(DLHS,nlhs)
  call mem_alloc(DRHS,nrhs)
  DO idmat = 1,nlhs
    DLHS(idmat)%p => DmatLHS(idmat)%p
  ENDDO
  DO idmat = 1,nrhs
    DRHS(idmat)%p => DmatRHS(idmat)%p
  ENDDO
ENDIF
CALL ls_attachDmatToSetting(DLHS,nlhs,setting,'LHS',1,2,.TRUE.,lupri)
CALL ls_attachDmatToSetting(DRHS,nrhs,setting,'RHS',3,4,.TRUE.,lupri)

nAtoms = setting%molecule(1)%p%nAtoms

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
   reCalcGab = SETTING%SCHEME%reCalcGab
   SETTING%SCHEME%reCalcGab = .TRUE.
   CALL ls_multipolemoment(LUPRI,LUERR,SETTING,nbast,0,  &
   & DLHS(1)%p%nrow,DLHS(1)%p%ncol,DRHS(1)%p%nrow,DRHS(1)%p%ncol,&
   & AORdefault,AORdefault,AORdefault,AORdefault,RegularSpec,ContractedInttype,.TRUE.)
   ! now derivative moments
   CALL ls_multipolemoment(LUPRI,LUERR,SETTING,nbast,0,  &
   & DLHS(1)%p%nrow,DLHS(1)%p%ncol,DRHS(1)%p%nrow,DRHS(1)%p%ncol,&
   & AORdefault,AORdefault,AORdefault,AORdefault,GradientSpec,ContractedInttype,.TRUE.)
   SETTING%SCHEME%MM_NOSCREEN = SAVE_SCREEN
END IF
CALL ls_jengine(AORdefault,AORdefault,AORdefault,AORdefault,&
     &          CoulombOperator,GradientSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,eeGRAD,.FALSE.)
IF(SETTING%SCHEME%FMM) THEN
   call ls_jengineClassicalGrad(eeGRAD,AORdefault,AORdefault,AORdefault,AORdefault,&
        & CoulombOperator,GradientSpec,ContractedInttype,SETTING,LUPRI,LUERR,nAtoms)
ENDIF
CALL ls_freeDmatFromSetting(setting)

IF(SETTING%SCHEME%FMM)THEN
   SETTING%SCHEME%DO_MMGRD = .FALSE.
   SETTING%SCHEME%NOSEGMENT = saveNOSEGMENT
   SETTING%SCHEME%ReCalcGab = ReCalcGab
ENDIF
call mem_dealloc(DLHS)
call mem_dealloc(DRHS)

IF (.NOT.same) call dscal(3*nAtoms,0.5E0_realk,eeGrad,1)

END SUBROUTINE II_get_J_gradient_regular


!> \brief Calculates the one-electron contribution to the molecular gradient
!> \author S. Reine
!> \date 2010-03-01
!> \param oneGrad The one-electron gradient contribution
!> \param Dmat The density matrices
!> \param ndmat The number of density matrices
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_oneElectron_gradient(oneGrad,Dmat,ndmat,setting,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Real(realk),intent(INOUT)       :: oneGrad(3,setting%molecule(1)%p%nAtoms)
Integer,intent(IN)            :: lupri,luerr,ndmat
Type(matrixp),intent(IN)      :: Dmat(ndmat)
!
Real(realk),pointer :: gradCont(:,:)
Integer :: nAtoms
!
! Nuclear-electron attraction gradient
CALL II_get_ne_gradient(oneGrad,Dmat,ndmat,setting,lupri,luerr)

nAtoms = setting%molecule(1)%p%nAtoms
CALL mem_alloc(gradCont,3,nAtoms)
! Kinetic energy gradient
CALL II_get_kinetic_gradient(gradCont,Dmat,ndmat,setting,lupri,luerr)
!
CALL DAXPY(3*nAtoms,1E0_realk,gradCont,1,oneGrad,1)
CALL mem_dealloc(gradCont)
!
END SUBROUTINE II_get_oneElectron_gradient


!> \brief Calculates the nuclear-electron attraction contribution to the molecular gradient
!> \author S. Reine
!> \date 2010-03-01
!> \param neGrad The nuclear-electron attraction gradient
!> \param Dmat The density matrices
!> \param ndmat The number of density matrices
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_ne_gradient(neGrad,Dmat,ndmat,setting,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Real(realk),intent(INOUT)       :: neGrad(3,setting%molecule(1)%p%nAtoms)
Integer,intent(IN)            :: lupri,luerr,ndmat
Type(matrixp),intent(IN)      :: Dmat(ndmat)
!
Type(matrixp)       :: Dmat_AO(ndmat)
integer             :: nAtoms,I
real(realk)         :: OLDTHRESH
logical :: saveNOSEGMENT,reCalcGab
call time_II_operations1()

!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
IF(ndmat.NE.1)call lsquit('This option in II_get_ne_gradient have not been verified',-1)
IF(SETTING%SCHEME%FMM) THEN
! the classical contribution to the N-E gradient is calculated together with the 
! j_gradient 
! here we calculated only the non-classical contribution


!we turn off family type basis sets because it does not work
!for FMM-GRADIENTS - when calculating both 1 and 2 electron 
!contributions together
saveNOSEGMENT = SETTING%SCHEME%NOSEGMENT
SETTING%SCHEME%NOSEGMENT = .TRUE.
!recalc primscreening matrix
reCalcGab = SETTING%SCHEME%reCalcGab
SETTING%SCHEME%reCalcGab = .TRUE.
ENDIF 

IF(setting%IntegralTransformGC)THEN
   DO I=1,ndmat
      allocate(Dmat_AO(I)%p)
      CALL mat_init(Dmat_AO(I)%p,Dmat(1)%p%nrow,Dmat(1)%p%ncol)
      call GCAO2AO_transform_matrixD2(Dmat(I)%p,Dmat_AO(I)%p,setting,lupri)
   ENDDO
   CALL ls_attachDmatToSetting(Dmat_AO,ndmat,setting,'LHS',1,2,.TRUE.,lupri)
ELSE
   CALL ls_attachDmatToSetting(Dmat,ndmat,setting,'LHS',1,2,.TRUE.,lupri)
ENDIF

nAtoms = setting%molecule(1)%p%nAtoms

call initIntegralOutputDims(setting%output,3,nAtoms,1,1,1)
CALL ls_getIntegrals(AORdefault,AORdefault,AONuclear,AOempty,&
     &NucpotOperator,GradientSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,neGrad,.FALSE.)

IF(SETTING%SCHEME%FMM) THEN
   !restore values
   SETTING%SCHEME%NOSEGMENT = saveNOSEGMENT
   SETTING%SCHEME%reCalcGab = reCalcGab
ENDIF

CALL ls_freeDmatFromSetting(setting)
IF(setting%IntegralTransformGC)THEN
   DO I=1,ndmat
      CALL mat_free(Dmat_AO(I)%p)
      deallocate(Dmat_AO(I)%p)
   ENDDO
ENDIF
call time_II_operations2(JOB_II_get_ne_gradient)

END SUBROUTINE II_get_ne_gradient

!> \brief Calculates the kinetic energy contribution to the molecular gradient
!> \author S. Reine
!> \date 2010-03-01
!> \param kinGrad The kinetic energy gradient
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_kinetic_gradient(kinGrad,Dmat,ndmat,setting,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Real(realk),intent(INOUT)       :: kinGrad(3,setting%molecule(1)%p%nAtoms)
Integer,intent(IN)            :: lupri,luerr,ndmat
Type(matrixp),intent(IN)      :: Dmat(ndmat)
Type(matrixp)                 :: Dmat_AO(ndmat)
!
integer             :: nAtoms,I
real(realk)         :: OLDTHRESH
call time_II_operations1()
IF(ndmat.NE.1)call lsquit('This option in II_get_kinetic_gradient have not been verified',-1)
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR

IF(setting%IntegralTransformGC)THEN
   DO I=1,ndmat
      allocate(Dmat_AO(I)%p)
      CALL mat_init(Dmat_AO(I)%p,Dmat(1)%p%nrow,Dmat(1)%p%ncol)
      call GCAO2AO_transform_matrixD2(Dmat(I)%p,Dmat_AO(I)%p,setting,lupri)
   ENDDO
   CALL ls_attachDmatToSetting(Dmat_AO,ndmat,setting,'LHS',1,3,.TRUE.,lupri)
ELSE
   CALL ls_attachDmatToSetting(Dmat,ndmat,setting,'LHS',1,3,.TRUE.,lupri)
ENDIF

nAtoms = setting%molecule(1)%p%nAtoms
call initIntegralOutputDims(setting%output,3,nAtoms,1,1,1)
CALL ls_getIntegrals(AORdefault,AOempty,AORdefault,AOempty,&
     &KineticOperator,GradientSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,kinGrad,.FALSE.)
CALL ls_freeDmatFromSetting(setting)
IF(setting%IntegralTransformGC)THEN
   DO I=1,ndmat
      CALL mat_free(Dmat_AO(I)%p)
      deallocate(Dmat_AO(I)%p)
   ENDDO
ENDIF
call time_II_operations2(JOB_II_get_kinetic_gradient)
END SUBROUTINE II_get_kinetic_gradient

!> \brief Calculates the nuclear repulsion energy contribution
!> \author S. Reine and T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param nucpot the nuclear repulsion energy contribution
SUBROUTINE II_get_nucpot(LUPRI,LUERR,SETTING,NUCPOT)
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
integer             :: usemat
INTEGER               :: LUPRI,LUERR
REAL(realk)           :: nucpot
Integer               :: I,J
real(realk)           :: pq(3),distance
logical               :: NOBQBQ
call time_II_operations1()
NOBQBQ = SETTING%SCHEME%NOBQBQ
NUCPOT=0.0E0_realk
DO I=1,SETTING%MOLECULE(1)%p%Natoms
 IF(SETTING%MOLECULE(1)%p%ATOM(I)%phantom)CYCLE
 DO J=I+1,SETTING%MOLECULE(1)%p%Natoms
  IF(SETTING%MOLECULE(1)%p%ATOM(J)%phantom)CYCLE
  
  if(setting%molecule(1)%p%ATOM(I)%Pointcharge .and. &
    &setting%molecule(1)%p%ATOM(J)%Pointcharge .and. NOBQBQ) cycle

  pq(1) = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(1)-SETTING%MOLECULE(1)%p%ATOM(J)%CENTER(1)
  pq(2) = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(2)-SETTING%MOLECULE(1)%p%ATOM(J)%CENTER(2)
  pq(3) = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(3)-SETTING%MOLECULE(1)%p%ATOM(J)%CENTER(3)
  Distance = sqrt(pq(1)*pq(1)+pq(2)*pq(2)+pq(3)*pq(3))
  NUCPOT=NUCPOT+SETTING%MOLECULE(1)%p%ATOM(I)%Charge*SETTING%MOLECULE(1)%p%ATOM(J)%Charge/Distance
 ENDDO
ENDDO
call time_II_operations2(JOB_II_get_nucpot)
END SUBROUTINE II_get_nucpot

!> \brief Calculates the nuclear-nuclear contribution to the molecular gradient
!> \author S. Reine
!> \date 2010-02-29
!> \param nucGrad The nuclear-nuclear molecular gradient vector
!> \param settings The settings for integral evaluation
SUBROUTINE II_get_nn_gradient(nucGrad,SETTING,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING),intent(IN) :: SETTING
REAL(realk),intent(INOUT)    :: nucGrad(3,SETTING%MOLECULE(1)%p%Natoms)
Integer,intent(IN)         :: lupri,luerr
!
Integer     :: I,J
real(realk) :: pq(3),distance,temp,extent2
logical     :: screen
logical               :: NOBQBQ
call time_II_operations1()
NOBQBQ = SETTING%SCHEME%NOBQBQ

if( SETTING%SCHEME%FMM .AND. (.NOT. SETTING%SCHEME%MM_NO_ONE)) then
   write(*,*) 'Attention: the printed nuclear repulsion gradient will only '
   write(*,*) '           include the non-classical part!'
   write(*,*) '           The classical part will be given as part of the '
   write(*,*) '           electron repulsion gradient'
end if
screen = .true.
extent2 = 2E0_realk
!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!
! extent2 = extent(nuc1) + extent(nuc2) = 2.0 
! since,  the artificial extent of the nuclei is set to 1.0E0_realk
! this has to match the settings in the routines 
!    ls_mm_read_in_raw_data (mm_interface.f90) 
!    getODcenter (LSint/ODbatch.f90) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

nucGrad=0E0_realk
DO I=1,SETTING%MOLECULE(1)%p%Natoms
 IF(SETTING%MOLECULE(1)%p%ATOM(I)%phantom)CYCLE
 DO J=I+1,SETTING%MOLECULE(1)%p%Natoms
  IF(SETTING%MOLECULE(1)%p%ATOM(J)%phantom)CYCLE
  IF(setting%molecule(1)%p%ATOM(I)%Pointcharge .and. &
       &setting%molecule(1)%p%ATOM(J)%Pointcharge .and. NOBQBQ) cycle

  pq(1) = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(1)-SETTING%MOLECULE(1)%p%ATOM(J)%CENTER(1)
  pq(2) = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(2)-SETTING%MOLECULE(1)%p%ATOM(J)%CENTER(2)
  pq(3) = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(3)-SETTING%MOLECULE(1)%p%ATOM(J)%CENTER(3)
  Distance = sqrt(pq(1)*pq(1)+pq(2)*pq(2)+pq(3)*pq(3))
  if( SETTING%SCHEME%FMM .AND. (.NOT. SETTING%SCHEME%MM_NO_ONE)) then
     if((Distance .le. extent2)) then
        screen = .true.
     else
        screen = .false.
     endif
  endif
  if(screen) then
     temp = SETTING%MOLECULE(1)%p%ATOM(I)%Charge*&
          &     SETTING%MOLECULE(1)%p%ATOM(J)%Charge/Distance**3
     nucGrad(1,I) = nucGrad(1,I) - pq(1)*temp
     nucGrad(2,I) = nucGrad(2,I) - pq(2)*temp
     nucGrad(3,I) = nucGrad(3,I) - pq(3)*temp
     nucGrad(3,J) = nucGrad(3,J) + pq(3)*temp
     nucGrad(1,J) = nucGrad(1,J) + pq(1)*temp
     nucGrad(2,J) = nucGrad(2,J) + pq(2)*temp
  end if
 ENDDO
ENDDO
call time_II_operations2(JOB_II_get_nn_gradient)
END SUBROUTINE II_get_nn_gradient

!> \brief Calculates the nuclear contribution to the quadrupole moment
!> (FIXME add reference)
!> \author T. Kjaergaard
!> \date 2010-02-29
!> \param Q The nuclear contribution to the quadrupole moment
!> \param settings The settings for integral evaluation
SUBROUTINE II_get_nuc_Quad(Q,SETTING,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING),intent(IN) :: SETTING
REAL(realk),intent(INOUT)    :: Q(3,3)
Integer,intent(IN)         :: lupri,luerr
!
Integer     :: I
real(realk) :: X,Y,Z,Charge,R
real(realk),parameter :: D3=3.0E0_realk
logical               :: NOBQBQ
NOBQBQ = SETTING%SCHEME%NOBQBQ
Q=0E0_realk
DO I=1,SETTING%MOLECULE(1)%p%Natoms
   IF(SETTING%MOLECULE(1)%p%ATOM(I)%phantom)CYCLE
   IF(setting%molecule(1)%p%ATOM(I)%Pointcharge.and. NOBQBQ) cycle

   X = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(1)
   Y = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(2)
   Z = SETTING%MOLECULE(1)%p%ATOM(I)%CENTER(3)
   R = X*X+Y*Y+Z*Z
   Charge = SETTING%MOLECULE(1)%p%ATOM(I)%Charge
   Q(1,1) = Q(1,1) + Charge * (D3 * X*X - R)
   Q(2,2) = Q(2,2) + Charge * (D3 * Y*Y - R)
   Q(3,3) = Q(3,3) + Charge * (D3 * Z*Z - R)
   Q(2,1) = Q(2,1) + Charge * (D3 * Y*X )
   Q(3,1) = Q(3,1) + Charge * (D3 * Z*X )
   Q(3,2) = Q(3,2) + Charge * (D3 * Z*Y )
ENDDO
Q(1,2) = Q(2,1)
Q(1,3) = Q(3,1)
Q(2,3) = Q(3,2)

END SUBROUTINE II_get_nuc_Quad

!> \brief Calculates the cartesian moments 
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param carmom the cartesian moments 
!> \param nmat = (nderiv+1)*(nderiv+2)*(nderiv+3)/6
!> \param nderiv The order of cartesian moments
!> \param X the x-component of the multipole expansion
!> \param Y the y-component of the multipole expansion
!> \param Z the z-component of the multipole expansion
SUBROUTINE II_get_carmom(LUPRI,LUERR,SETTING,carmom,nmat,nderiv,X,Y,Z)
IMPLICIT NONE
INTEGER,intent(IN)              :: LUPRI,LUERR,nderiv,nmat
TYPE(MATRIX),intent(INOUT)      :: carmom(nmat)
TYPE(LSSETTING),intent(INOUT)   :: SETTING
Real(realk),intent(IN)          :: X,Y,Z
!
INTEGER :: nbast

!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
SETTING%SCHEME%OD_MOM = .FALSE.
SETTING%SCHEME%MOM_CENTER = (/ X,Y,Z /)
nbast = carmom(1)%nrow
SETTING%SCHEME%CMorder = nderiv
call initIntegralOutputDims(setting%output,nbast,nbast,1,1,nmat)
CALL ls_getIntegrals(AORdefault,AORdefault,AOempty,AOempty,&
     & CarmomOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
SETTING%SCHEME%CMorder = 0
CALL retrieve_Output(lupri,setting,carmom,setting%IntegralTransformGC)

END SUBROUTINE II_get_carmom

!> \brief Calculates the cartesian moments 
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param carmom the cartesian moments 
!> \param nmat the number of matrices 
!> \param nderiv The order of cartesian moments
!> \param X the x-component of the multipole expansion
!> \param Y the y-component of the multipole expansion
!> \param Z the z-component of the multipole expansion
SUBROUTINE II_get_single_carmom(LUPRI,LUERR,SETTING,carmom,imat,nderiv,X,Y,Z)
IMPLICIT NONE
INTEGER,intent(IN)              :: LUPRI,LUERR,nderiv,imat
TYPE(MATRIX),intent(INOUT)      :: carmom
TYPE(LSSETTING),intent(INOUT)   :: SETTING
Real(realk),intent(IN)          :: X,Y,Z
!
INTEGER :: nbast

!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
SETTING%SCHEME%OD_MOM = .FALSE.
SETTING%SCHEME%MOM_CENTER = (/ X,Y,Z /)
nbast = carmom%nrow
SETTING%SCHEME%CMorder = nderiv
SETTING%SCHEME%CMimat = imat
call initIntegralOutputDims(setting%output,nbast,nbast,1,1,1)
CALL ls_getIntegrals(AORdefault,AORdefault,AOempty,AOempty,&
     & CarmomOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
SETTING%SCHEME%CMorder = 0
SETTING%SCHEME%CMimat = 0
CALL retrieve_Output(lupri,setting,carmom,setting%IntegralTransformGC)
END SUBROUTINE II_get_single_carmom

!> \brief Calculates the dipole moment contribution from the nuclei
!> \param setting Integral evalualtion settings
!> \param nucdip Nuclear contrib to dipole moment
SUBROUTINE II_get_nucdip(setting,nucdip)
IMPLICIT NONE
TYPE(LSSETTING),intent(in) :: setting
REAL(realk),intent(inout)    :: nucdip(3)
INTEGER            :: i
TYPE(ATOMITEM),pointer :: mol(:)
nucdip = 0
mol => setting%MOLECULE(1)%p%ATOM
DO i=1,setting%MOLECULE(1)%p%Natoms
  IF(mol(i)%phantom)CYCLE
  nucdip = nucdip + mol(i)%Charge * mol(i)%CENTER
ENDDO
END SUBROUTINE II_get_nucdip

!> \brief General property integrals
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param matArray the matrices of property integrals
!> \param nOperatorComp the number of matrices 
!> \param Oper the label of property integral
!>
!> IMPLEMENTED SO FAR
!> Oper   :nOperatorComp    Ordering  
!> DIPLEN :3                X,Y,Z comp
!> ANGMOM :3                X,Y,Z comp
!> ANGLON :3                X,Y,Z comp
!> 1ELPOT :1  
!> LONMOM1:3                X,Y,Z comp
!> LONMOM2:3                X,Y,Z comp
!> LONMOM :3                X,Y,Z comp
!> MAGMOM :3                X,Y,Z comp
!> PSO    :3*NATOMS   (3,NATOM) X,Y,Z comp for each atom
!> NSTLON :9*NATOMS   (3,3,NATOM) (magnetic X,Y,Z, atomic X,Y,Z, for each Atom
!> NSTNOL :9*NATOMS   (3,3,NATOM) (magnetic X,Y,Z, atomic X,Y,Z, for each Atom
!> NST    :9*NATOMS   (3,3,NATOM) (magnetic X,Y,Z, atomic X,Y,Z, for each Atom
!> DIPVEL :3                X,Y,Z comp
!> ROTSTR :6                XX,XY,XZ,YY,YZ,ZZ
!> THETA  :6                XX,XY,XZ,YY,YZ,ZZ
!> OVERLAP:1
!> D-CM1  :9          (3,3)       (Electric X,Y,Z, Magnetic X,Y,Z)
!> D-CM2  :18         (3,6)       (Electric X,Y,Z, Magnetic XX,XY,XZ,YY,YZ,ZZ)
!> I will implement more upon request. T. Kjaergaard
recursive SUBROUTINE II_get_prop(LUPRI,LUERR,SETTING,MatArray,nOperatorComp,Oper)
IMPLICIT NONE
Character(len=7)     :: Oper
INTEGER              :: LUPRI,LUERR,nOperatorComp
TYPE(LSSETTING)      :: SETTING
TYPE(MATRIX)         :: MatArray(nOperatorComp)
!
integer :: nbast,I,Operparam
Character(len=7)     :: Oper2
TYPE(MATRIX),pointer :: TMPMatArray(:)
real(realk) :: TS,TE
logical :: CS_screenSAVE, PS_screenSAVE
call time_II_operations1()
CALL LSTIMER('START ',TS,TE,LUPRI)
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
call param_oper_paramfromString(Oper,Operparam)

setting%scheme%do_prop = .TRUE.
setting%scheme%propoper = Operparam

nbast = MatArray(1)%nrow
IF(Oper.EQ.'MAGMOM ')THEN
   !MAGMOM is a sum of 2 contributions ANGLON and LONMOM
   !kinetic energy part of LONMOM
   call II_get_prop(LUPRI,LUERR,SETTING,MatArray,nOperatorComp,'LONMOM1')
   call mem_alloc(TMPMatArray,nOperatorComp)
   do I=1,nOperatorComp
      call mat_init(TMPMatArray(I),nbast,nbast)
   enddo
   !nuclear potential energy part of LONMOM
   call II_get_prop(LUPRI,LUERR,SETTING,TMPMatArray,nOperatorComp,'LONMOM2')
   do I=1,nOperatorComp
      call mat_daxpy(1.0E0_realk,TMPMatArray(I),MatArray(I))
   enddo   
!  ANGLON
   call II_get_prop(LUPRI,LUERR,SETTING,TMPMatArray,nOperatorComp,'ANGLON ')
   do I=1,nOperatorComp
      call mat_daxpy(0.5E0_realk,TMPMatArray(I),MatArray(I))
      call mat_free(TMPMatArray(I))
   enddo
   call mem_dealloc(TMPMatArray)
ELSEIF(Oper.EQ.'LONMOM2')THEN
   call initIntegralOutputDims(setting%output,nbast,nbast,1,1,nOperatorComp)
   CALL ls_getIntegrals(AORdefault,AORdefault,AONuclear,AOempty,&
        & Operparam,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
   CALL retrieve_Output(lupri,setting,MatArray,setting%IntegralTransformGC)
ELSEIF(Oper.EQ.'LONMOM ')THEN
   ! 2 contributions kinetic energy part + nuclear potential energy part 
   call II_get_prop(LUPRI,LUERR,SETTING,MatArray,nOperatorComp,'LONMOM1')
   call mem_alloc(TMPMatArray,nOperatorComp)
   do I=1,nOperatorComp
      call mat_init(TMPMatArray(I),nbast,nbast)
   enddo
   call II_get_prop(LUPRI,LUERR,SETTING,TMPMatArray,nOperatorComp,'LONMOM2')
   do I=1,nOperatorComp
      call mat_daxpy(1.0E0_realk,TMPMatArray(I),MatArray(I))
      call mat_free(TMPMatArray(I))
   enddo
   call mem_dealloc(TMPMatArray)
ELSEIF(Oper.EQ.'NST    ')THEN
   ! 2 contributions lonodon(NSTLON) + nonlondon(NSTNOL)
   call II_get_prop(LUPRI,LUERR,SETTING,MatArray,nOperatorComp,'NSTLON ')
   call mem_alloc(TMPMatArray,nOperatorComp)
   do I=1,nOperatorComp
      call mat_init(TMPMatArray(I),nbast,nbast)
   enddo
   call II_get_prop(LUPRI,LUERR,SETTING,TMPMatArray,nOperatorComp,'NSTNOL ')
   do I=1,nOperatorComp
      call mat_daxpy(1.0E0_realk,TMPMatArray(I),MatArray(I))
      call mat_free(TMPMatArray(I))
   enddo
   call mem_dealloc(TMPMatArray)
ELSEIF(Oper.EQ.'1ELPOT ')THEN
   call initIntegralOutputDims(setting%output,nbast,nbast,1,1,nOperatorComp)
   CALL ls_getIntegrals(AORdefault,AORdefault,AONuclear,AOempty,&
        & Operparam,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
   CALL retrieve_Output(lupri,setting,MatArray,setting%IntegralTransformGC)   
ELSEIF(Oper.EQ.'PSO    '.OR.(Oper.EQ.'NSTLON '.OR.Oper.EQ.'NSTNOL '))THEN
   call initIntegralOutputDims(setting%output,nbast,nbast,1,1,nOperatorComp)
   CALL ls_getIntegrals(AORdefault,AORdefault,AONuclear,AOempty,&
        & Operparam,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
   CALL retrieve_Output(lupri,setting,MatArray,setting%IntegralTransformGC)   
ELSE
   call initIntegralOutputDims(setting%output,nbast,nbast,1,1,nOperatorComp)
   CALL ls_getIntegrals(AORdefault,AORdefault,AOempty,AOempty,&
        & Operparam,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
   CALL retrieve_Output(lupri,setting,MatArray,setting%IntegralTransformGC)
ENDIF

setting%scheme%do_prop = .FALSE.
CALL LSTIMER('PropInt:'//Oper,TS,TE,LUPRI)
call time_II_operations2(JOB_II_get_prop)

END SUBROUTINE II_get_prop

!> \brief specific nuclei moment derivative 
!> \author T. Kjaergaard
!> \date 2015
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param matArray the matrices of property integrals
!> \param iAtom chosen atom to calculate the 3 PSO X,Y,Z components
recursive SUBROUTINE II_get_PSO_spec(LUPRI,LUERR,SETTING,MatArray,iAtom)
IMPLICIT NONE
INTEGER              :: LUPRI,LUERR,IATOM
TYPE(LSSETTING)      :: SETTING
TYPE(MATRIX)         :: MatArray(3)
!
integer :: nbast,I,Operparam
Character(len=7)     :: Oper2
TYPE(MATRIX),pointer :: TMPMatArray(:)
real(realk) :: TS,TE
logical :: CS_screenSAVE, PS_screenSAVE
call time_II_operations1()
CALL LSTIMER('START ',TS,TE,LUPRI)
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR

call param_oper_paramfromString('PSO    ',Operparam)
setting%scheme%do_prop = .TRUE.
setting%scheme%propoper = Operparam
nbast = MatArray(1)%nrow
setting%scheme%AONuclearSpecID = iAtom
call initIntegralOutputDims(setting%output,nbast,nbast,1,1,3)
CALL ls_getIntegrals(AORdefault,AORdefault,AONuclearSpec,AOempty,&
     & Operparam,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,MatArray,setting%IntegralTransformGC)   

setting%scheme%do_prop = .FALSE.
CALL LSTIMER('PSOspec',TS,TE,LUPRI)
call time_II_operations2(JOB_II_get_prop)

END SUBROUTINE II_get_PSO_spec

!> \brief General routine for calculation of expectation value of property integrals
!> So in Other words does dotproduct(PropIntegral,D)
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param expval the matrices of property integrals
!> \param Dmat the matrices the should be contracted with the property integral
!> \param nOperatorComp the number of matrices 
!> \param Oper the label of property integral
!>
!> IMPLEMENTED SO FAR
!> Oper   :nOperatorComp    Ordering  
!> PSO    :3*NATOMS   (3,NATOM) X,Y,Z comp for each atom
!> NSTLON :9*NATOMS   (3,3,NATOM) (magnetic X,Y,Z, atomic X,Y,Z, for each Atom
!> NSTNOL :9*NATOMS   (3,3,NATOM) (magnetic X,Y,Z, atomic X,Y,Z, for each Atom
!> NST    :9*NATOMS   (3,3,NATOM) (magnetic X,Y,Z, atomic X,Y,Z, for each Atom
!> I will implement more upon request. T. Kjaergaard
recursive SUBROUTINE II_get_prop_expval(LUPRI,LUERR,SETTING,expval,Dmat,ndmat,nOperatorComp,Oper)
IMPLICIT NONE
Character(len=7)     :: Oper
INTEGER              :: LUPRI,LUERR,nOperatorComp
TYPE(LSSETTING)      :: SETTING
TYPE(MATRIX)         :: DMat(ndmat)
real(realk)          :: expval(nOperatorComp*ndmat)
!
TYPE(MATRIX)         :: DMat_AO(ndmat)
integer :: I,ndmat,J,operparam
real(realk) :: TS,TE
logical :: CS_screenSAVE, PS_screenSAVE
real(realk),pointer :: TMPexpval(:)
call time_II_operations1()
CALL LSTIMER('START ',TS,TE,LUPRI)
!set threshold 
call param_oper_paramfromString(Oper,Operparam)
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
setting%scheme%do_prop = .TRUE.
setting%scheme%propoper = Operparam
IF(Oper.EQ.'NST    ')THEN
   ! 2 contributions lonodon(NSTLON) + nonlondon(NSTNOL)
   call II_get_prop_expval(LUPRI,LUERR,SETTING,expval,Dmat,ndmat,nOperatorComp,'NSTLON ')
   call mem_alloc(TMPexpval,nOperatorComp)
   call II_get_prop_expval(LUPRI,LUERR,SETTING,TMPexpval,Dmat,ndmat,nOperatorComp,'NSTNOL ')
   do I=1,nOperatorComp
      expval(I) = expval(I) + TMPexpval(I) 
   enddo
   call mem_dealloc(TMPexpval)
ELSE
   IF(setting%IntegralTransformGC)THEN
      DO I=1,ndmat
         CALL mat_init(Dmat_AO(I),Dmat(1)%nrow,Dmat(1)%ncol)
         call GCAO2AO_transform_matrixD2(Dmat(I),Dmat_AO(I),setting,lupri)
      ENDDO
      CALL ls_attachDmatToSetting(Dmat_AO,ndmat,setting,'LHS',1,2,.TRUE.,lupri)
   ELSE
      CALL ls_attachDmatToSetting(Dmat,ndmat,setting,'LHS',1,2,.TRUE.,lupri)
   ENDIF
   call initIntegralOutputDims(setting%output,1,1,1,1,nOperatorComp*ndmat)
   CALL ls_getIntegrals(AORdefault,AORdefault,AONuclear,AOempty,&
        & Operparam,EcontribSpec,ContractedInttype,SETTING,LUPRI,LUERR)
   CALL retrieve_Output(lupri,setting,expval,setting%IntegralTransformGC)   
   CALL ls_freeDmatFromSetting(setting)
   IF(setting%IntegralTransformGC)THEN
      DO I=1,ndmat
         CALL mat_free(Dmat_AO(I))
      ENDDO
   ENDIF
ENDIF
setting%scheme%do_prop = .FALSE.
CALL LSTIMER('Propval:'//Oper,TS,TE,LUPRI)
call time_II_operations2(JOB_II_get_prop_expval)

END SUBROUTINE II_get_prop_expval

!> \brief routine for calculation of expectation value of Specific nuclei PSO integrals
!> So in Other words does dotproduct(PSOIntegral,D)
!> \author T. Kjaergaard
!> \date 2015
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param expval the matrices of property integrals Order: (3 PSO components, ndmat)
!> \param Dmat the matrices the should be contracted with the property integral
!> \param nDmat
!> \param iAtom Chosen atom to calculate the PSO integral for 
recursive SUBROUTINE II_get_pso_spec_expval(LUPRI,LUERR,SETTING,expval,Dmat,ndmat,iAtom)
IMPLICIT NONE
INTEGER              :: LUPRI,LUERR,iAtom
TYPE(LSSETTING)      :: SETTING
TYPE(MATRIX)         :: DMat(ndmat)
real(realk)          :: expval(3*ndmat) !Order: 3,ndmat
!
TYPE(MATRIX)         :: DMat_AO(ndmat)
integer :: I,ndmat,J,operparam
real(realk) :: TS,TE
logical :: CS_screenSAVE, PS_screenSAVE
real(realk),pointer :: TMPexpval(:)
call time_II_operations1()
CALL LSTIMER('START ',TS,TE,LUPRI)
!set threshold 
call param_oper_paramfromString('PSO    ',Operparam)
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
setting%scheme%do_prop = .TRUE.
setting%scheme%propoper = Operparam

IF(setting%IntegralTransformGC)THEN
   DO I=1,ndmat
      CALL mat_init(Dmat_AO(I),Dmat(1)%nrow,Dmat(1)%ncol)
      call GCAO2AO_transform_matrixD2(Dmat(I),Dmat_AO(I),setting,lupri)
   ENDDO
   CALL ls_attachDmatToSetting(Dmat_AO,ndmat,setting,'LHS',1,2,.TRUE.,lupri)
ELSE
   CALL ls_attachDmatToSetting(Dmat,ndmat,setting,'LHS',1,2,.TRUE.,lupri)
ENDIF
setting%scheme%AONuclearSpecID = iAtom
call initIntegralOutputDims(setting%output,1,1,1,1,3*ndmat)
CALL ls_getIntegrals(AORdefault,AORdefault,AONuclearSpec,AOempty,&
     & Operparam,EcontribSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,expval,setting%IntegralTransformGC)   
CALL ls_freeDmatFromSetting(setting)
IF(setting%IntegralTransformGC)THEN
   DO I=1,ndmat
      CALL mat_free(Dmat_AO(I))
   ENDDO
ENDIF
setting%scheme%do_prop = .FALSE.
CALL LSTIMER('PSOSpecExpVal:',TS,TE,LUPRI)
call time_II_operations2(JOB_II_get_prop_expval)

END SUBROUTINE II_get_pso_spec_expval

!> \brief routine for calculation of expectation value of Specific nuclei PSO integrals
!> So in Other words does dotproduct(PSOIntegral,D)
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param expval the matrices of property integrals (3,3) order: Magnetic comp, Magnetic Moment comp
!> \param Dmat the matrices the should be contracted with the property integral
!> \param iAtom the chosen atom for which this should be calculated
recursive SUBROUTINE II_get_nst_spec_expval(LUPRI,LUERR,SETTING,expval,Dmat,iAtom)
IMPLICIT NONE
INTEGER              :: LUPRI,LUERR,iAtom
TYPE(LSSETTING)      :: SETTING
TYPE(MATRIX)         :: DMat(1)
real(realk)          :: expval(9)
!
TYPE(MATRIX)         :: DMat_AO(1)
integer :: I,ndmat,J,operparam
real(realk) :: TS,TE
logical :: CS_screenSAVE, PS_screenSAVE
real(realk) :: TMPexpval(9)
call time_II_operations1()
CALL LSTIMER('START ',TS,TE,LUPRI)
ndmat = 1
IF(setting%IntegralTransformGC)THEN
   CALL mat_init(Dmat_AO(1),Dmat(1)%nrow,Dmat(1)%ncol)
   call GCAO2AO_transform_matrixD2(Dmat(1),Dmat_AO(1),setting,lupri)
   CALL ls_attachDmatToSetting(Dmat_AO,ndmat,setting,'LHS',1,2,.TRUE.,lupri)
ELSE
   CALL ls_attachDmatToSetting(Dmat,ndmat,setting,'LHS',1,2,.TRUE.,lupri)
ENDIF
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR

!set threshold 
call param_oper_paramfromString('NSTLON ',Operparam)
setting%scheme%do_prop = .TRUE.
setting%scheme%propoper = Operparam
setting%scheme%AONuclearSpecID = iAtom
call initIntegralOutputDims(setting%output,1,1,1,1,9)
CALL ls_getIntegrals(AORdefault,AORdefault,AONuclearSpec,AOempty,&
     & Operparam,EcontribSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,TMPexpval,setting%IntegralTransformGC)   

call param_oper_paramfromString('NSTNOL ',Operparam)
setting%scheme%do_prop = .TRUE.
setting%scheme%propoper = Operparam
setting%scheme%AONuclearSpecID = iAtom
call initIntegralOutputDims(setting%output,1,1,1,1,9)
CALL ls_getIntegrals(AORdefault,AORdefault,AONuclearSpec,AOempty,&
     & Operparam,EcontribSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,expval,setting%IntegralTransformGC)   
do I=1,9
   expval(I) = expval(I) + TMPexpval(I) 
enddo
CALL ls_freeDmatFromSetting(setting)
IF(setting%IntegralTransformGC)THEN   
   CALL mat_free(Dmat_AO(1))
ENDIF
setting%scheme%do_prop = .FALSE.
CALL LSTIMER('NSTSpecExpVal:',TS,TE,LUPRI)
call time_II_operations2(JOB_II_get_prop_expval)

END SUBROUTINE II_get_nst_spec_expval

!> \brief Calculates property integrals
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param mat the matrix of property integrals
!> \param nmat the number of matrices 
!> \param STRING the label of property integral
!>
!> IMPLEMENTED SO FAR 
!> DIPLEN :3  XDIPLEN ,YDIPLEN ,ZDIPLEN
!> SECMOM :6  XXSECMOM,XYSECMOM,XZSECMOM,YYSECMOM,YZSECMOM,ZZSECMOM
!> WILL IMPLEMT AT SOME POINT
!> ANGMOM :3  XANGMOM,YANGMOM,ZANGMOM
!>
SUBROUTINE II_get_integral(LUPRI,LUERR,SETTING,mat,nmat,STRING)
IMPLICIT NONE
character(len=7)     :: STRING
TYPE(LSSETTING)      :: SETTING
INTEGER              :: LUPRI,LUERR,nderiv,nmat
TYPE(MATRIX),target  :: mat(nmat)
!
integer              :: Oper
TYPE(MATRIX),pointer :: TMP(:)
INTEGER              :: nbast,I,nmat2

!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
nbast = mat(1)%nrow
Oper = CarmomOperator
SELECT CASE(STRING)
CASE('DIPLEN ')
   IF(nmat.NE. 3)call LSQUIT('II_get_integral:DIPLEN INPUT REQUIRES A DIM = 3',lupri)
   nmat2=4
   nderiv = 1
CASE('SECMOM ')
   IF(nmat.NE. 6)call LSQUIT('II_get_integral:SECMOM INPUT REQUIRES A DIM = 6',lupri)
   nmat2=10
   nderiv = 2
CASE DEFAULT
   WRITE(LUPRI,'(1X,2A)') 'Wrong case in II_get_integral =',STRING
   CALL LSQUIT('Wrong case in II_get_integral',lupri)
END SELECT

call initIntegralOutputDims(setting%output,nbast,nbast,1,1,nmat2)
SETTING%SCHEME%CMORDER = nderiv
CALL ls_getIntegrals(AORdefault,AORdefault,AOempty,AOempty,&
     & Oper,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
SETTING%SCHEME%CMORDER = 0
call mem_ALLOC(TMP,nmat2)
DO I=1,nmat2
   call mat_init(TMP(I),nbast,nbast)
ENDDO
CALL retrieve_Output(lupri,setting,TMP,setting%IntegralTransformGC)

SELECT CASE(STRING)
CASE('DIPLEN ')
   DO I=1,3
      call mat_assign(mat(I),tmp(I+1))
   ENDDO
CASE('SECMOM ')
   DO I=1,6
      call mat_assign(mat(I),tmp(I+4))
   ENDDO
END SELECT

DO I=1,nmat2
   call mat_free(TMP(I))
ENDDO
call mem_DEALLOC(TMP)

END SUBROUTINE II_get_integral

SUBROUTINE II_get_integral_full(LUPRI,LUERR,SETTING,mat,nbast,nmat,STRING)
IMPLICIT NONE
character(len=7)     :: STRING
TYPE(LSSETTING)      :: SETTING
INTEGER              :: LUPRI,LUERR,nderiv,nmat,nbast
real(realk)          :: mat(nbast,nbast,nmat)
!
integer              :: Oper
real(realk),pointer  :: TMP(:,:,:)
INTEGER              :: I,nmat2
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
Oper = CarmomOperator
SELECT CASE(STRING)
CASE('DIPLEN ')
   IF(nmat.NE. 3)call LSQUIT('II_get_integral:DIPLEN INPUT REQUIRES A DIM = 3',lupri)
   nmat2=4
   nderiv = 1
CASE('SECMOM ')
   IF(nmat.NE. 6)call LSQUIT('II_get_integral:SECMOM INPUT REQUIRES A DIM = 6',lupri)
   nmat2=10
   nderiv = 2
CASE DEFAULT
   WRITE(LUPRI,'(1X,2A)') 'Wrong case in II_get_integral =',STRING
   CALL LSQUIT('Wrong case in II_get_integral',lupri)
END SELECT

call initIntegralOutputDims(setting%output,nbast,nbast,1,1,nmat2)
SETTING%SCHEME%CMORDER = nderiv
CALL ls_getIntegrals(AORdefault,AORdefault,AOempty,AOempty,&
     & Oper,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
SETTING%SCHEME%CMORDER = 0
call mem_ALLOC(TMP,nbast,nbast,nmat2)
CALL retrieve_Output(lupri,setting,TMP,setting%IntegralTransformGC)
SELECT CASE(STRING)
CASE('DIPLEN ')
   CALL DAXPY(3*nbast*nbast,1E0_realk,tmp(:,:,2:4),1,MAT,1)
CASE('SECMOM ')
   CALL DAXPY(3*nbast*nbast,1E0_realk,tmp(:,:,5:10),1,MAT,1)
END SELECT
call mem_DEALLOC(TMP)
END SUBROUTINE II_get_integral_full

!> \brief Calculates the spherical moments
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param sphmom the matrices of the spherical moments
!> \param nSpherMom the number of matrices 
!> \param mulMomOrder the order of moments
!> \param X the x-component of the multipole expansion
!> \param Y the y-component of the multipole expansion
!> \param Z the z-component of the multipole expansion
SUBROUTINE II_get_sphmom(LUPRI,LUERR,SETTING,sphmom,nSpherMom,mulMomOrder,X,Y,Z)
IMPLICIT NONE
INTEGER                    :: LUPRI,LUERR,mulMomOrder,nSpherMom,usemat,nbast,I
TYPE(MATRIX),intent(INOUT) :: sphmom(nSpherMom)
TYPE(LSSETTING)            :: SETTING
Real(realk),intent(IN)     :: X,Y,Z
!
TYPE(MATRIX),pointer :: carmom(:)
Integer              :: ncarmom,icarmom

nbast = sphmom(1)%nrow
ncarmom = (mulMomOrder+1)*(mulMomOrder+2)*(mulMomOrder+3)/6
call mem_ALLOC(carmom,ncarmom)
DO icarmom=1,ncarmom
  call mat_init(carmom(icarmom),nbast,nbast)
  call mat_zero(carmom(icarmom))
ENDDO

call II_get_carmom(LUPRI,LUERR,SETTING,carmom,ncarmom,mulMomOrder,X,Y,Z)

call II_carmom_to_shermom(sphmom,carmom,nSpherMom,ncarmom,mulMomOrder,lupri,setting%scheme%intprint)

DO icarmom=1,ncarmom
  call mat_free(carmom(icarmom))
ENDDO
call mem_DEALLOC(carmom)

END SUBROUTINE II_get_sphmom

SUBROUTINE II_carmom_to_shermom(sphmom,carmom,nSpherMom,ncarmom,mulMomOrder,lupri,iprint)
IMPLICIT NONE
INTEGER,INTENT(IN)         :: LUPRI,iprint,mulMomOrder,nSpherMom,nCarMom
TYPE(MATRIX),intent(INOUT) :: sphmom(nSpherMom)
TYPE(MATRIX),intent(IN)    :: carmom(nCarMom)
!
REAL(REALK),pointer    :: TRANSMAT(:)
REAL(REALK),PARAMETER :: D0 = 0E0_realk
!
Integer :: nSpher,nCart,angMom,I,J,startI,startJ
Real(realk) :: coef

! Alloc with the maximum possible size nCarMom*nSpherMom
call mem_alloc(TRANSMAT,nCarMom*nSpherMom)

startJ=0
startI=0
DO ANGMOM=0,mulMomOrder
   nSpher  = 2*ANGMOM + 1
   nCart = (ANGMOM + 1)*(ANGMOM + 2)/2
   CALL BUILD_CART_TO_SPH_MAT(ANGMOM,TRANSMAT,nCarMom*nSpherMom,nSpher,nCart,lupri,iprint)
   DO I = 1, nSpher
      DO J = 1, nCart
         COEF = TRANSMAT(J+(I-1)*nCart)
         IF (ABS(COEF) .GT. D0) THEN
           call mat_daxpy(coef,carmom(startJ+J),sphmom(startI+I))
         END IF
      ENDDO
   ENDDO
   startI=startI+nSpher
   startJ=startJ+nCart
ENDDO

call mem_dealloc(TRANSMAT)

END SUBROUTINE II_carmom_to_shermom

!> \brief Calculates the 3 center overlap integrals
!> \author T. Kjaergaard
!> \date 2010
!> \param lupriold old Default print unit
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param nbast The number of basis functions
!> \param auxoption If the auxiliary basis should be used
SUBROUTINE II_get_3center_overlap(LUPRIOLD,LUPRI,LUERR,SETTING,nbast,auxoption)
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,nbast,lupriold,i,j,k,nbastaux,iprint
TYPE(MATRIX),target   :: TMP
integer             :: usemat
LOGICAL               :: auxoption
Real(realk),pointer   :: integrals(:,:,:,:,:)
Real(realk)           :: SUM
type(matrixp)         :: intmat(1)
logical :: IntegralTransformGC

!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
IntegralTransformGC=.FALSE.
iprint=SETTING%SCHEME%INTPRINT
IF(auxoption)THEN
   nbastaux = SETTING%MOLECULE(1)%p%nbastAUX
   call initIntegralOutputDims(setting%output,nbast,nbast,1,nbastaux,1)
   CALL ls_getIntegrals(AORdefault,AORdefault,AOempty,AODFdefault,&
        &OverlapOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
   call mem_alloc(integrals,nbast,nbast,1,nbastaux,1)
   CALL retrieve_Output(lupri,setting,integrals,IntegralTransformGC)
ELSE
   nbastaux = nbast
   call initIntegralOutputDims(setting%output,nbast,nbast,1,nbast,1)
   CALL ls_getIntegrals(AORdefault,AORdefault,AOempty,AORdefault,&
        &OverlapOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
   call mem_alloc(integrals,nbast,nbast,1,nbast,1)
   CALL retrieve_Output(lupri,setting,integrals,IntegralTransformGC)
ENDIF

SUM=0
IF(IPRINT.GT. 100)THEN
   WRITE(LUPRIOLD,'(2X,A)')'3 Center Overlap (a|b|c)'
   WRITE(LUPRIOLD,'(2X,A)')'FORMAT IS (a,b) c=1, c=2 , ..'
ENDIF
DO i=1,nbast
   DO j=1,nbast
      IF(IPRINT.GT. 100)THEN
         WRITE(LUPRIOLD,'(2X,A1,I2,A1,I2,A2,5F16.6/,(10X,5F16.6))')&
              &'(',i,',',j,')=',&
              &(Integrals(i,j,1,k,1),k=1,nbastaux)
      ENDIF
   ENDDO
ENDDO
WRITE(LUPRIOLD,'(2X,A,F20.12)')'SUM OF ALL 3 CENTER OVERLAPS '
WRITE(LUPRIOLD,'(2X,A,F20.12)')'FROM THE STANDARD DRIVER ',SUM
call mem_dealloc(integrals)

END SUBROUTINE II_get_3center_overlap

!> \brief Calculates the 2 center electron repulsion integrals (eri)
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param F the matrix containing the 2 center eri
SUBROUTINE II_get_2center_eri(LUPRI,LUERR,SETTING,F)
IMPLICIT NONE
TYPE(MATRIX),target   :: F
integer               :: usemat
TYPE(LSSETTING)       :: SETTING
TYPE(INTEGRALINPUT)   :: INTINPUT
TYPE(INTEGRALOUTPUT)  :: INTOUTPUT
INTEGER               :: LUPRI,LUERR
!
Real(realk),pointer :: integrals(:,:,:,:,:)
Integer             :: nbast
type(matrixp)         :: intmat(1)

!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%J_THR
nbast = F%nrow
call initIntegralOutputDims(setting%output,nbast,1,nbast,1,1)
CALL ls_getIntegrals(AORdefault,AOempty,AORdefault,AOempty,CoulombOperator,RegularSpec,ContractedInttype,&
      &              SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,F,setting%IntegralTransformGC)

END SUBROUTINE II_get_2center_eri

!> \brief 2 center eris from two different basis sets
SUBROUTINE II_get_2center_mixed_eri(LUPRI,LUERR,SETTING,S,AO1,AO2,GCAO1,GCAO2)
IMPLICIT NONE
TYPE(MATRIX)          :: S
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,IPRINT,AO1,AO2
LOGICAL               :: GCAO1,GCAO2
!
Integer             :: i,j,LU,nbast2,nbast1

!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR

IPRINT=SETTING%SCHEME%INTPRINT
nbast1 = getNbasis(AO1,ContractedintType,SETTING%MOLECULE(1)%p,LUPRI)
nbast2 = getNbasis(AO2,ContractedintType,SETTING%MOLECULE(2)%p,LUPRI)

IF(nBast1.NE.S%nrow)CALL LSQUIT('dim1 mismatch in II_get_mixed_overlap',-1)
IF(nBast2.NE.S%ncol)CALL LSQUIT('dim2 mismatch in II_get_mixed_overlap',-1)
call initIntegralOutputDims(setting%output,nbast1,1,nbast2,1,1)
CALL ls_getIntegrals(AO1,AOempty,AO2,AOempty,&
     &CoulombOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,S,.FALSE.)

IF (GCAO1) call AO2GCAO_half_transform_matrix(S,SETTING,LUPRI,1)
IF (GCAO2) call AO2GCAO_half_transform_matrix(S,SETTING,LUPRI,2)

IF(IPRINT.GT. 1000)THEN
   WRITE(LUPRI,'(A,2X,F16.8)')'Mixed 2centerERI',mat_dotproduct(S,S)
   call mat_print(S,1,nbast1,1,nbast2,lupri)
ENDIF

END SUBROUTINE II_get_2center_mixed_eri



!> \brief Calculates the explicit 4 center eri tensor in Mulliken (ab|cd) or Dirac noation <a(1)c(2)|r_12^-1|b(1)d(2)>
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param outputintegral the explicit 4 center eri
!> \param DIM1 the dimension of orbital alpha
!> \param DIM2 the dimension of orbital beta
!> \param DIM3 the dimension of orbital delta
!> \param DIM4 the dimension of orbital gamma
!> \param Dirac Specifies Dirac or Mulliken (default) notation
SUBROUTINE II_get_4center_eri(LUPRI,LUERR,SETTING,outputintegral,dim1,dim2,dim3,dim4,intspec,dirac)
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,dim1,dim2,dim3,dim4
REAL(REALK),target,intent(inout) :: outputintegral(dim1,dim2,dim3,dim4)
Logical,optional,intent(in)   :: dirac
character,intent(in)  :: intspec(5)
!
Logical             :: dirac_format
REAL(REALK),pointer :: integrals(:,:,:,:)
integer             :: i,j,k,l,Oper
real(realk)         :: coeff(6),exponent(6),tmp
real(realk)         :: coeff2(21),sumexponent(21),prodexponent(21)
integer             :: IJ,nGaussian,nG2,ao(4),dummy
real(realk)         :: GGem

DO i=1,4
   IF (intSpec(i).NE.'R')call lsquit('Only R allowed in II_get_4center_eri',-1)
ENDDO

IF (intSpec(5).NE.'C') THEN
   nGaussian = 6
   nG2 = nGaussian*(nGaussian+1)/2
   GGem = 0E0_realk
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
ENDIF

! ***** SELECT OPERATOR TYPE *****
IF (intSpec(5).EQ.'C') THEN
   ! Regular Coulomb operator 1/r12
   oper = CoulombOperator
ELSE IF (intSpec(5).EQ.'G') THEN
   ! The Gaussian geminal operator g
   oper = GGemOperator
   call set_GGem(Setting%GGem,coeff,exponent,nGaussian)
ELSE IF (intSpec(5).EQ.'F') THEN
   ! The Gaussian geminal divided by the Coulomb operator g/r12
   oper = GGemCouOperator
   call set_GGem(Setting%GGem,coeff,exponent,nGaussian)
ELSE IF (intSpec(5).EQ.'D') THEN
   ! The double commutator [[T,g],g]
   oper = GGemGrdOperator
   call set_GGem(Setting%GGem,coeff2,sumexponent,prodexponent,nG2)
ELSE IF (intSpec(5).EQ.'2') THEN
   ! The Gaussian geminal operator squared g^2
   oper = GGemOperator
   call set_GGem(Setting%GGem,coeff2,sumexponent,prodexponent,nG2)
ELSE
   call lsquit('Error in specification of operator in ',-1)
ENDIF

!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%J_THR

!Specify Dirac or Mulliken integral format
dirac_format = .FALSE.
IF (present(dirac)) THEN
  dirac_format = dirac
ENDIF

IF (dirac_format) THEN
  call initIntegralOutputDims(setting%output,dim1,dim3,dim2,dim4,1)
  call mem_alloc(integrals,dim1,dim3,dim2,dim4)
ELSE
  call initIntegralOutputDims(setting%output,dim1,dim2,dim3,dim4,1)
  integrals => outputintegral
ENDIF
 
CALL ls_getIntegrals(AORdefault,AORdefault,AORdefault,AORdefault,&
     &Oper,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)

IF(setting%IntegralTransformGC)THEN
  call lsquit('Error in II_get_4center_eri - IntegralTransformGC not implemented',lupri)
ELSE
  CALL retrieve_Output(lupri,setting,integrals,setting%IntegralTransformGC)
  IF (dirac_format) THEN
    DO l=1,dim4
     DO k=1,dim3
      DO j=1,dim2
       DO i=1,dim1
         outputintegral(i,j,k,l) = integrals(i,k,j,l)
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   call mem_dealloc(integrals)
  ENDIF
ENDIF

END SUBROUTINE II_get_4center_eri

!> \brief Calculates the differentiated 4 center eri tensor in Mulliken (ab|cd) or Dirac noation <a(1)c(2)|r_12^-1|b(1)d(2)>
!> \author S. Reine
!> \date 2013-02-05
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param outputintegral the explicit 4 center eri
!> \param DIM1 the dimension of orbital alpha
!> \param DIM2 the dimension of orbital beta
!> \param DIM3 the dimension of orbital delta
!> \param DIM4 the dimension of orbital gamma
!> \param DIM5 the number of differential components
!> \param Dirac Specifies Dirac or Mulliken (default) notation
!> \param geoderiv Specifies the geoemetrical derivative order
SUBROUTINE II_get_4center_eri_diff(LUPRI,LUERR,SETTING,outputintegral,dim1,dim2,dim3,dim4,dim5,dirac,geoderiv)
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,dim1,dim2,dim3,dim4,dim5
REAL(REALK),target,intent(inout) :: outputintegral(:,:,:,:,:) !dim1,dim2,dim3,dim4,dim5
Logical,optional      :: dirac
Integer,optional      :: geoderiv
!
Logical             :: dirac_format,dogeoderiv,nofamily
REAL(REALK),pointer :: integrals(:,:,:,:,:)
integer             :: i,j,k,l,n,intSpec,geoOrder

!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%J_THR

!Specify geometrical derivative order if any
dogeoderiv = .FALSE.
geoOrder   = 0
IF (present(geoderiv)) THEN
  dogeoderiv = geoderiv.NE.0
  geoOrder   = geoderiv
ENDIF

IF (dogeoderiv) THEN
  IF (geoderiv.GE.1) THEN
!   intSpec = GeoDerivLHSSpec
!   intSpec = GeoDerivRHSSpec
    intSpec = GeoDerivSpec
  ELSE
     call lsquit('Error in II_get_4center_eri_diff - only first order geometrical derivative integrals implemented',lupri)
  ENDIF

  IF (.NOT.setting%scheme%nofamily) CALL LSQUIT('II_get_4center_eri_diff only working with .NOFAMILY keyword option',-1)
ELSE
  intSpec = RegularSpec
ENDIF


!Specify Dirac or Mulliken integral format
dirac_format = .FALSE.
IF (present(dirac)) THEN
  dirac_format = dirac
ENDIF

IF (dirac_format) THEN
  call initIntegralOutputDims(setting%output,dim1,dim3,dim2,dim4,dim5)
  call mem_alloc(integrals,dim1,dim3,dim2,dim4,dim5)
  call ls_dzero(integrals,dim1*dim2*dim3*dim4*dim5)
ELSE
  call initIntegralOutputDims(setting%output,dim1,dim2,dim3,dim4,dim5)
  IF (dogeoderiv.AND.(intSpec.NE.GeoDerivSpec)) THEN
    call mem_alloc(integrals,dim1,dim2,dim3,dim4,dim5)
    call ls_dzero(integrals,dim1*dim2*dim3*dim4*dim5)
  ELSE
    integrals => outputintegral
  ENDIF
ENDIF
 
CALL ls_getIntegrals(AORdefault,AORdefault,AORdefault,AORdefault,&
     &CoulombOperator,intSpec,ContractedInttype,SETTING,LUPRI,LUERR,geoOrder)

IF(setting%IntegralTransformGC)THEN
  call lsquit('Error in II_get_4center_eri_diff - IntegralTransformGC not implemented',lupri)
ELSE
  CALL retrieve_Output(lupri,setting,integrals,setting%IntegralTransformGC)
  IF (dogeoderiv) THEN
    IF ((geoderiv.EQ.1).AND.(intSpec.NE.GeoDerivSpec)) THEN
     DO n=1,dim5
      DO l=1,dim4
       DO k=1,dim3
        DO j=1,dim2
         DO i=1,dim1
          outputintegral(i,j,k,l,n) = integrals(i,j,k,l,n) + integrals(k,l,i,j,n)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
     IF (dirac_format) THEN
       call dcopy(dim1*dim2*dim3*dim4*dim5,outputintegral,1,integrals,1)
     ELSE
       call mem_dealloc(integrals)
     ENDIF
    ENDIF
  ENDIF
  IF (dirac_format) THEN
   DO n=1,dim5
    DO l=1,dim4
     DO k=1,dim3
      DO j=1,dim2
       DO i=1,dim1
         outputintegral(i,j,k,l,n) = integrals(i,k,j,l,n)
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
   call mem_dealloc(integrals)
  ENDIF
ENDIF

END SUBROUTINE II_get_4center_eri_diff

!> \brief Calculates the differentiated nuclear-electron attraction integral tensor: d/dR_e sum_C (ab|C)Z_C 
!> \author S. Reine
!> \date 2013-02-05
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param outputintegral the differentiated 1-electron integral - for nucel the
!>        nuclei are summed over according to  d/dR_e sum_C (ab|C)Z_C
!> \param DIM1 the dimension of orbital alpha
!> \param DIM2 the dimension of orbital beta
!> \param DIM5 the number of differential components
!> \param geoderiv Specifies the geoemetrical derivative order
SUBROUTINE II_get_1el_diff(LUPRI,LUERR,SETTING,outputintegral,oneElType,dim1,dim2,dim5,geoderiv)
IMPLICIT NONE
TYPE(LSSETTING)         :: SETTING
INTEGER                 :: LUPRI,LUERR,dim1,dim2,dim5
REAL(REALK),target,intent(inout) :: outputintegral(:,:,:,:,:)
Character(*),intent(IN) :: oneElType
Integer,optional        :: geoderiv
!
Logical             :: dogeoderiv,nofamily
integer             :: i,j,k,l,n,intSpec,geoOrder
integer             :: AO1,AO2,AO3,AO4,Oper,n1,n2,n3,n4,n5

!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%J_THR

!Specify geometrical derivative order if any
dogeoderiv = .FALSE.
geoOrder   = 0
IF (present(geoderiv)) THEN
  dogeoderiv = geoderiv.NE.0
  geoOrder   = geoderiv
ENDIF

IF (dogeoderiv) THEN
  IF (geoderiv.GE.1) THEN
    intSpec = GeoDerivSpec
  ENDIF

  !Simen Hack - GeoDerivSpec currently not working with family-type basis sets
  nofamily = setting%scheme%nofamily
  setting%scheme%nofamily = .TRUE.
ELSE
  intSpec = RegularSpec
ENDIF

n1   = dim1
n2   = dim2
n3   = 1
n4   = 1
n5   = dim5
AO1  = AORdefault
AO2  = AORdefault
AO3  = AOEmpty
AO4  = AOEmpty
select case(oneElType)
case ('nucel')
  AO3  = AONuclear
  Oper = NucpotOperator
case ('overlap')
  Oper = OverlapOperator
! intSpec = GeoDerivLHSSpec
case ('kinetic')
  n2   = 1
  n3   = dim2
  AO2  = AOEmpty
  AO3  = AORdefault
  Oper = KineticOperator
case default
  call lsquit('Error in II_get_1el_diff - Unknown 1-electron case',lupri)
end select

call initIntegralOutputDims(setting%output,n1,n2,n3,n4,n5)
call ls_dzero(outputintegral,n1*n2*n3*n4*n5)
 
CALL ls_getIntegrals(AO1,AO2,AO3,AO4,Oper,intSpec,ContractedInttype,SETTING,LUPRI,LUERR,geoOrder)

IF(setting%IntegralTransformGC)THEN
  call lsquit('Error in II_get_1el_diff - IntegralTransformGC not implemented',lupri)
ELSE
  CALL retrieve_Output(lupri,setting,outputintegral,setting%IntegralTransformGC)
  IF (dogeoderiv) THEN
    setting%scheme%nofamily = nofamily
  ENDIF
ENDIF

END SUBROUTINE II_get_1el_diff

!> \brief Calculates and stores the screening integrals
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
SUBROUTINE II_precalc_ScreenMat(LUPRI,LUERR,SETTING)
IMPLICIT NONE
TYPE(LSSETTING)      :: SETTING
INTEGER,intent(in)   :: LUPRI,LUERR
!
TYPE(lstensor),pointer :: GAB
LOGICAL :: IntegralTransformGC,CSintsave,PSintsave,FoundOnDisk
LOGICAL :: FoundInMem,CS_INT,PS_INT
INTEGER :: THR,I,ilst,J,K,dummy
CHARACTER(len=10) :: INTTYPE(2)
integer :: Oper(8)
Character(80)       :: Filename
Character(53)       :: identifier
Character(22)       :: label1,label2
real(realk)         :: OLDTHRESH
real(realk),pointer :: screenMat(:,:)
Integer   :: natoms,n1,n2
integer(kind=short) :: GABelm
integer(kind=short),pointer :: MAT2(:,:)
logical :: dofit
integer :: AO1,AO2,AO3,AO4
dummy=1
call time_II_operations1()

dofit = SETTING%SCHEME%DENSFIT .OR. SETTING%SCHEME%PARI_J .OR. &
     SETTING%SCHEME%PARI_K .OR. SETTING%SCHEME%MOPARI_K

IF(SETTING%SCHEME%saveGABtoMem)THEN
 IF(SETTING%SCHEME%CS_SCREEN.OR.SETTING%SCHEME%PS_SCREEN &
      & .OR.SETTING%SCHEME%MBIE_SCREEN)THEN
    !set threshold 
    SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*&
         &MIN(SETTING%SCHEME%ONEEL_THR,SETTING%SCHEME%J_THR,SETTING%SCHEME%K_THR)
    Oper(1)=CoulombOperator
    Oper(2)=NucleiOperator
    Oper(3)=ErfcOperator
    Oper(4)=CAMOperator
    Oper(5)=GGemOperator
    Oper(6)=GGemCouOperator
    Oper(7)=GGemGrdOperator
    Oper(8)=CoulombOperator !Density-fitting screening integrals
    CS_INT = SETTING%SCHEME%CS_SCREEN
    PS_INT = SETTING%SCHEME%PS_SCREEN
    DO J=1,size(Oper)
       IF(J.EQ.3.AND.(.NOT.SETTING%SCHEME%SR_EXCHANGE) )CYCLE
       IF(J.EQ.4.AND.(.NOT.SETTING%SCHEME%CAM) )CYCLE
       IF(J.EQ.5.AND.(.NOT.SETTING%GGEM%is_set) )CYCLE
       IF(J.EQ.6.AND.(.NOT.SETTING%GGEM%is_set) )CYCLE
       IF(J.EQ.7.AND.(.NOT.SETTING%GGEM%is_set) )CYCLE
       IF(J.EQ.8.AND.(.NOT.dofit) )CYCLE
       nullify(GAB)
       THR = ABS(NINT(LOG10(SETTING%SCHEME%intTHRESHOLD)))
       IF (ASSOCIATED(Setting%MOLECULE(1)%p)) THEN
          label1 = Setting%MOLECULE(1)%p%label
       ELSE
          label1 = 'Empty_________________'
       ENDIF
       IF (ASSOCIATED(Setting%MOLECULE(2)%p)) THEN
          label2 = Setting%MOLECULE(2)%p%label
       ELSE
          label2 = 'Empty_________________'
       ENDIF
       IF(THR.GT.9)THEN
          write(identifier,'(A2,I2,A1,A22,A1,A22,A1,L1,L1)')'CS',THR,'_',label1,'_',label2,'_',CS_INT,PS_INT
       ELSE
          write(identifier,'(A3,I1,A1,A22,A1,A22,A1,L1,L1)')'CS0',THR,'_',label1,'_',label2,'_',CS_INT,PS_INT
       ENDIF
       IF(Oper(J).EQ.NucleiOperator) THEN
          AO1 = AONuclear
          AO2 = AOEmpty
          AO3 = AONuclear
          AO4 = AOEmpty
       ELSE IF (J.EQ.8) THEN
          AO1 = AODFdefault
          AO2 = AOEmpty
          AO3 = AODFdefault
          AO4 = AOEmpty
       ELSE
          AO1 = AORdefault
          AO2 = AORdefault
          AO3 = AORdefault
          AO4 = AORdefault
       ENDIF
       CALL io_get_filename(Filename,identifier,AO1,AO2,AO3,AO4,&
           &0,0,oper(J),Contractedinttype,.FALSE.,LUPRI,LUERR)
       IF (.NOT.SETTING%SCHEME%saveGABtoMem)THEN
          allocate(GAB)
          call lstensor_nullify(GAB)
       ENDIF
       IntegralTransformGC=.FALSE.
       nAtoms = SETTING%MOLECULE(1)%p%natoms
       CSintsave = setting%scheme%CS_int
       PSintsave = setting%scheme%PS_int
       setting%scheme%CS_int = CS_INT
       setting%scheme%PS_int = PS_INT
       IF(SETTING%SCHEME%saveGABtoMem)THEN
          call screen_add_associate_item(GAB,FILENAME)
       ENDIF
       IF(Oper(J).EQ.NucleiOperator) THEN
          call mem_alloc(screenMat,Natoms,1)
          CALL ls_getNucScreenIntegrals(GAB,ScreenMat,Filename,AONuclear,AOEmpty,&
               & Oper(J),SETTING,LUPRI,LUERR,.TRUE.)
          call mem_dealloc(screenMat)
       ELSE
          IF(Oper(J).EQ.ErfcOperator)setting%scheme%CS_int = .TRUE.
          IF(Oper(J).EQ.CAMOperator)setting%scheme%CS_int = .TRUE.
          call initIntegralOutputDims(setting%Output,dummy,dummy,1,1,1)
          CALL ls_getIntegrals(AO1,AO2,AO3,AO4,&
               &Oper(J),RegularSpec,Contractedinttype,SETTING,LUPRI,LUERR)
          CALL retrieve_screen_Output(lupri,setting,GAB,IntegralTransformGC)   
       ENDIF
       setting%scheme%CS_int = CSintsave
       setting%scheme%PS_int = PSintsave
    ENDDO
 ENDIF
ENDIF
#ifdef VAR_MPI
call ls_mpibcast(IISCREEN,infpar%master,setting%comm)
call II_bcast_screen(setting%comm)
#endif
call time_II_operations2(JOB_II_precalc_ScreenMat)

END SUBROUTINE II_precalc_ScreenMat

#ifdef VAR_MPI
subroutine II_bcast_screen(comm)
implicit none
integer(kind=ls_mpik) :: comm,Master,mynum
logical :: Slave
call get_rank_for_comm(comm,mynum)
Master = infpar%master
Slave = mynum.NE.Master
call ls_mpiInitBuffer(Master,LSMPIBROADCAST,comm)
call mpicopy_screen(slave,master)
call ls_mpiFinalizeBuffer(Master,LSMPIBROADCAST,comm)
end subroutine II_bcast_screen

subroutine II_screeninit(comm)
implicit none
integer(kind=ls_mpik) :: comm
call screen_init
end subroutine II_screeninit

subroutine II_screenfree(comm)
implicit none
integer(kind=ls_mpik) :: comm
call screen_free
end subroutine II_screenfree
#endif

!> \brief Calculates the 4 center 2 eri screening mat
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param Gab the output matrix
SUBROUTINE II_get_2int_BatchScreenMat(LUPRI,LUERR,SETTING,nbatches,batchGAB,nbast)
IMPLICIT NONE
INTEGER,intent(in)                :: LUPRI,LUERR,nbatches,nbast
REAL(realk),intent(inout)         :: batchGAB(nbatches,nbatches)
TYPE(LSSETTING),intent(inout)     :: SETTING
!
REAL(realk),pointer        :: bastGAB(:,:)
IF(setting%IntegralTransformGC)THEN
   !I do not think it makes sense to transform afterwards 
   !so here the basis needs to be transformed
   call lsquit('II_get_2int_ScreenMat and IntegralTransformGC do not work',-1)
ENDIF
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
call initIntegralOutputDims(setting%Output,nbast,nbast,1,1,1)
setting%Output%RealGabMatrix = .TRUE.
CALL ls_getScreenIntegrals1(AORdefault,AORdefault,&
     &CoulombOperator,.TRUE.,.FALSE.,.FALSE.,SETTING,LUPRI,LUERR,.TRUE.)
setting%Output%RealGabMatrix = .FALSE.
call mem_alloc(bastGAB,nbast,nbast)
call ls_dzero(bastGAB,nbast*nbast)
CALL retrieve_Output(lupri,setting,bastGAB,setting%IntegralTransformGC)
call ConvertBASTGabToBatchesGab(nbast,nbatches,setting,bastGAB,batchGAB,&
     & lupri,luerr)
call mem_dealloc(bastGAB)
END SUBROUTINE II_get_2int_BatchScreenMat

subroutine ConvertBASTGabToBatchesGab(nbast,nbatches,setting,BastGAB,BatchGAB,&
     & lupri,luerr)
implicit none
TYPE(LSSETTING),intent(inout)  :: SETTING
real(realk),intent(in) :: BastGab(nbast,nbast)
integer,intent(in) ::  nbast,nbatches,lupri,luerr
real(realk),intent(inout) :: BatchGab(nbatches,nbatches)
!
type(AOITEM) :: AO
integer,pointer :: batchdim(:)
integer :: I,A,dimI,dimJ,J,offsetI,offsetJ,JJ,II,nbastout
real(realk) :: TMP

CALL setAObatch(AO,0,1,nbastout,AORdefault,Contractedinttype,Setting%scheme,&
     & Setting%fragment(1)%p,setting%basis(1)%p,lupri,luerr)
IF(nbastout.NE.nbast)call lsquit('Error in ConvertBASTGabToBatchesGab',-1)
call mem_alloc(batchdim,AO%nbatches)
do I=1,AO%nbatches
   dimI = 0
   DO A=1,AO%BATCH(I)%nAngmom
      dimI = dimI + AO%BATCH(I)%norbitals(A)
   ENDDO
   batchdim(I) = dimI
enddo

offsetJ = 0
do J=1,AO%nbatches
   dimJ = batchdim(J)
   offsetI = 0
   do I=1,AO%nbatches
      dimI = batchdim(I)
      TMP = 0.0E0_realk
      do JJ=1,dimJ
         do II=1,dimI            
            TMP = MAX(TMP,ABS(BastGAB(II+offsetI,JJ+offsetJ)))
         enddo
      enddo
      BatchGAB(I,J) = TMP
      offsetI = offsetI + dimI
   enddo
   offsetJ = offsetJ + dimJ
enddo
CALL free_aoitem(lupri,AO)
call mem_dealloc(batchdim)
end subroutine ConvertBASTGabToBatchesGab

!> \brief Calculates the 4 center 2 eri screening mat
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param Gab the output matrix
SUBROUTINE II_get_2int_ScreenMat(LUPRI,LUERR,SETTING,GAB)
IMPLICIT NONE
TYPE(MATRIX),target   :: GAB
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
!
Integer             :: nbast1,nbast2
!call lsquit('II_get_2int_ScreenMat not implemented ',-LUPRI)
IF(setting%IntegralTransformGC)THEN
   !I do not think it makes sense to transform afterwards 
   !so here the basis needs to be transformed
   call lsquit('II_get_2int_ScreenMat and IntegralTransformGC do not work',-1)
ENDIF
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
nbast1 = GAB%nrow
nbast2 = GAB%ncol
call mat_zero(GAB)
call initIntegralOutputDims(setting%Output,nbast1,nbast2,1,1,1)
setting%Output%RealGabMatrix = .TRUE.
CALL ls_getScreenIntegrals1(AORdefault,AORdefault,&
     &CoulombOperator,.TRUE.,.FALSE.,.FALSE.,SETTING,LUPRI,LUERR,.TRUE.)
setting%Output%RealGabMatrix = .FALSE.
CALL retrieve_Output(lupri,setting,GAB,setting%IntegralTransformGC)

END SUBROUTINE II_get_2int_ScreenMat

!> \brief Calculates the 4 center 2 eri screening mat
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param Gab the output matrix
SUBROUTINE II_get_2int_ScreenRealMat(LUPRI,LUERR,SETTING,nbast,GAB)
IMPLICIT NONE
INTEGER               :: LUPRI,LUERR,nbast
real(realk)           :: GAB(nbast,nbast)
TYPE(LSSETTING)       :: SETTING
IF(setting%IntegralTransformGC)THEN
   !I do not think it makes sense to transform afterwards 
   !so here the basis needs to be transformed
   call lsquit('II_get_2int_ScreenMat and IntegralTransformGC do not work',-1)
ENDIF
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
call ls_dzero(GAB,nbast*nbast)
call initIntegralOutputDims(setting%Output,nbast,nbast,1,1,1)
setting%Output%RealGabMatrix = .TRUE.
CALL ls_getScreenIntegrals1(AORdefault,AORdefault,&
     &CoulombOperator,.TRUE.,.FALSE.,.FALSE.,SETTING,LUPRI,LUERR,.TRUE.)
setting%Output%RealGabMatrix = .FALSE.
CALL retrieve_Output(lupri,setting,GAB,setting%IntegralTransformGC)
END SUBROUTINE II_get_2int_ScreenRealMat

!> \brief Calculates get the maxGabelm eri screening mat
!> \author J. Rekkedal
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evaluation settings
!> \param maxGabelm the output integer(short)
SUBROUTINE II_get_maxGabelm_ScreenMat(LUPRI,LUERR,SETTING,nbast,maxGABelm)
IMPLICIT NONE
INTEGER(short),intent(out)   :: maxGABelm
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,nbast
!
!call lsquit('II_get_2int_ScreenMat not implemented ',-LUPRI)
CALL ls_setDefaultFragments(setting)
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
call initIntegralOutputDims(setting%Output,nbast,nbast,1,1,1)
CALL ls_getScreenIntegrals1(AORdefault,AORdefault,&
     &CoulombOperator,.TRUE.,.FALSE.,.FALSE.,SETTING,LUPRI,LUERR,.TRUE.)

     setting%Output%resulttensor => setting%Output%ScreenTensor

CALL retrieve_Output(lupri,setting,maxGABelm,setting%IntegralTransformGC)

END SUBROUTINE II_get_maxGabelm_ScreenMat

!> \brief set the Incremental scheme
!> \author T. Kjaergaard
!> \date 2010
!> \param scheme the integral scheme
!> \param increm the increm logical
SUBROUTINE II_setIncremental(scheme,increm)
implicit none
TYPE(LSINTSCHEME),INTENT(INOUT) :: scheme
LOGICAL,INTENT(IN)              :: increm
scheme%incremental = increm
END SUBROUTINE II_setIncremental

!> \brief build BatchOrbitalInfo
!> \author T. Kjaergaard
!> \date 2010
!> \param setting Integral evalualtion settings
!> \param AO options AORdefault,AODFdefault or AOempty
!> \param intType options contracted or primitive
!> \param nbast the number of basis functions
!> \param orbtobast for a given basis function it provides the batch index
!> \param nBatches the number of batches
!> \param lupri Default print unit
!> \param luerr Default error print unit
SUBROUTINE II_getBatchOrbitalInfo(Setting,nBast,orbToBatch,nBatches,lupri,luerr)
implicit none
TYPE(LSSETTING),intent(IN) :: Setting
Integer,intent(IN)         :: nBast,lupri,luerr
Integer,intent(INOUT)      :: orbToBatch(nBast)
Integer,intent(INOUT)      :: nBatches
!
integer   :: AO,intType
TYPE(BATCHORBITALINFO) :: BO
TYPE(AOITEM)           :: AObuild
Integer                :: nDim
!
AO = AORdefault
intType = Contractedinttype
CALL initBatchOrbitalInfo(BO,nBast)
CALL setAObatch(AObuild,0,1,nDim,AO,intType,Setting%scheme,Setting%fragment(1)%p,&
     &          setting%basis(1)%p,lupri,luerr)
CALL setBatchOrbitalInfo(BO,AObuild,lupri)
orbToBatch = BO%orbToBatch
nBatches   = BO%nBatches
CALL freeBatchOrbitalInfo(BO)
CALL free_aoitem(lupri,AObuild)
END SUBROUTINE II_getBatchOrbitalInfo

!> \brief Calculates the reorthonormalization gradient term
!> \author S. Reine
!> \date 2010-03-22
!> \param reOrtho The reorthonormalization gradient term
!> \param DFD The matrix product DFD = D F D
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_reorthoNormalization(reOrtho,DFDmat,ndmat,setting,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Real(realk),intent(INOUT)       :: reOrtho(3,setting%molecule(1)%p%nAtoms)
Integer,intent(IN)            :: lupri,luerr,ndmat
Type(matrixp),intent(IN)      :: DFDmat(ndmat)
Logical :: GC
GC = setting%integralTransformGC
!
call II_get_reorthoNormalization_mixed(reOrtho,DFDmat,ndmat,AORdefault,AORdefault,GC,GC,setting,lupri,luerr)
!
END SUBROUTINE II_get_reorthoNormalization


!> \brief Calculates the reorthonormalization gradient term
!> \author S. Reine
!> \date 2010-03-22
!> \param reOrtho The reorthonormalization gradient term
!> \param DFD The matrix product DFD = D F D
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_reorthoNormalization_mixed(reOrtho,DFDmat,ndmat,AO1,AO2,GC1,GC2,setting,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Real(realk),intent(INOUT)     :: reOrtho(3,setting%molecule(1)%p%nAtoms)
Integer,intent(IN)            :: lupri,luerr,ndmat,AO1,AO2
Type(matrixp),intent(IN)      :: DFDmat(ndmat)
Logical,intent(IN)            :: GC1,GC2
!
Type(matrixp)       :: DFDmat_AO(ndmat)
integer             :: nAtoms,I
real(realk)         :: OLDTHRESH
call time_II_operations1()
IF(ndmat.NE.1)call lsquit('option not verified in II_get_reorthoNormalization',-1)
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR

IF (GC1.OR.GC2) THEN
   DO I=1,ndmat
      allocate(DFDmat_AO(I)%p)
      CALL mat_init(DFDmat_AO(I)%p,DFDmat(1)%p%nrow,DFDmat(1)%p%ncol)
      CALL mat_assign(DFDmat_AO(I)%p,DFDmat(I)%p)
      IF (GC1) call GCAO2AO_half_transform_matrix(DFDmat_AO(I)%p,SETTING,LUPRI,1)
      IF (GC2) call GCAO2AO_half_transform_matrix(DFDmat_AO(I)%p,SETTING,LUPRI,2)
   ENDDO
   CALL ls_attachDmatToSetting(DFDmat_AO,ndmat,setting,'LHS',1,2,.TRUE.,lupri)
ELSE
   CALL ls_attachDmatToSetting(DFDmat,ndmat,setting,'LHS',1,2,.TRUE.,lupri)
ENDIF

nAtoms = setting%molecule(1)%p%nAtoms

call initIntegralOutputDims(setting%output,3,nAtoms,1,1,1)
CALL ls_getIntegrals(AO1,AO2,AOempty,AOempty,&
     &OverlapOperator,GradientSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,reOrtho,setting%IntegralTransformGC)
CALL ls_freeDmatFromSetting(setting)
IF(GC1.OR.GC2)THEN
   DO I=1,ndmat
      CALL mat_free(DFDmat_AO(I)%p)
      deallocate(DFDmat_AO(I)%p)
   ENDDO
ENDIF
call time_II_operations2(JOB_II_get_reorthoNormalization)

END SUBROUTINE II_get_reorthoNormalization_mixed

!> \brief Calculates the reorthonormalization gradient term (LHS)
!> \author S. Reine
!> \date 2010-03-22
!> \param reOrtho The reorthonormalization gradient term
!> \param DFD The matrix product DFD = D F D
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_reorthoNormalization2(reOrtho2,DFDmat,ndmat,setting,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Real(realk),intent(INOUT)       :: reOrtho2(3,setting%molecule(1)%p%nAtoms)
Integer,intent(IN)            :: lupri,luerr,ndmat
Type(matrixp),intent(IN)      :: DFDmat(ndmat)
!
Type(matrixp)       :: DFDmat_AO(ndmat)
integer             :: nAtoms,I
real(realk)         :: OLDTHRESH
call time_II_operations1()
IF(ndmat.NE.1)call lsquit('option not verified in II_get_reorthoNormalization',-1)
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR

IF(setting%IntegralTransformGC)THEN
   DO I=1,ndmat
      allocate(DFDmat_AO(I)%p)
      CALL mat_init(DFDmat_AO(I)%p,DFDmat(1)%p%nrow,DFDmat(1)%p%ncol)
      call GCAO2AO_transform_matrixD2(DFDmat(I)%p,DFDmat_AO(I)%p,setting,lupri)
   ENDDO
   CALL ls_attachDmatToSetting(DFDmat_AO,ndmat,setting,'LHS',1,3,.TRUE.,lupri)
ELSE
   CALL ls_attachDmatToSetting(DFDmat,ndmat,setting,'LHS',1,3,.TRUE.,lupri)
ENDIF

nAtoms = setting%molecule(1)%p%nAtoms

call initIntegralOutputDims(setting%output,3,nAtoms,1,1,1)
CALL ls_getIntegrals(AORdefault,AOempty,AORdefault,AOempty,&
     &OverlapOperator,GradientSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,reOrtho2,setting%IntegralTransformGC)
CALL ls_freeDmatFromSetting(setting)
IF(setting%IntegralTransformGC)THEN
   DO I=1,ndmat
      CALL mat_free(DFDmat_AO(I)%p)
      deallocate(DFDmat_AO(I)%p)
   ENDDO
ENDIF
call time_II_operations2(JOB_II_get_reorthoNormalization2)

END SUBROUTINE II_get_reorthoNormalization2

!> \brief Calculates the geometric derivative of the kinetic matrix
!> \author T. Kjaergaard
!> \date 2011-08-22
!> \param hx The 3*Natoms geometric derivative components of the overlap matrix
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_geoderivKinetic(hx,natoms,setting,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Integer,intent(IN)            :: lupri,luerr,natoms
Type(matrix),intent(INOUT)    :: Hx(3*natoms) !hx,hy,hz for each atom
!
integer             :: nbast
real(realk)         :: OLDTHRESH
call time_II_operations1()
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
nbast = Hx(1)%nrow
call initIntegralOutputDims(setting%output,nbast,1,nbast,1,3*natoms)
CALL ls_getIntegrals(AORdefault,AOempty,AORdefault,AOempty,&
     &KineticOperator,GeoDerivLHSSpec,ContractedInttype,SETTING,LUPRI,LUERR)
!doing GeoDerivLHSSpec we only differentiate on the LHS so we need to symmetrize
setting%output%postprocess = SymmetricPostprocess
CALL retrieve_Output(lupri,setting,Hx,setting%IntegralTransformGC)
call time_II_operations2(JOB_II_get_geoderivOverlap)

END SUBROUTINE II_get_geoderivKinetic

!> \brief Calculates the geometric derivative of the Nuclear Attraction matrix
!> \author T. Kjaergaard
!> \date 2011-08-22
!> \param hx The 3*Natoms geometric derivative components of the overlap matrix
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_geoderivnucel(hx,natoms,setting,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Integer,intent(IN)            :: lupri,luerr,natoms
Type(matrix),intent(INOUT)    :: Hx(3*natoms) !hx,hy,hz for each atom
!
integer             :: nbast,I
real(realk)         :: OLDTHRESH
!FIX ME is this necessary 
do I=1,3*natoms
   CALL mat_zero(hx(I))
enddo
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
nbast = Hx(1)%nrow
IF(SETTING%SCHEME%FMM) THEN
   call lsquit('combi II_get_geoderivnucel and FMM not tested',-1)
ENDIF
call initIntegralOutputDims(setting%output,nbast,nbast,1,1,3*natoms)
CALL ls_getIntegrals(AORdefault,AORdefault,AONuclear,AOempty,&
     &NucpotOperator,GeoDerivLHSSpec,ContractedInttype,SETTING,LUPRI,LUERR)
!doing GeoDerivLHSSpec we only differentiate on the LHS so we need to symmetrize
setting%output%postprocess = SymmetricPostprocess
CALL retrieve_Output(lupri,setting,hx,setting%IntegralTransformGC)
END SUBROUTINE II_get_geoderivnucel 

!> \brief Calculates the geometric derivative of the overlap matrix
!> \author T. Kjaergaard
!> \date 2011-08-22
!> \param Sx The 3*Natoms geometric derivative components of the overlap matrix
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_geoderivOverlap(Sx,natoms,setting,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Integer,intent(IN)            :: lupri,luerr,natoms
Type(matrix),intent(INOUT)    :: Sx(3*natoms) !Sx,Sy,Sz for each atom
!
integer             :: nbast
real(realk)         :: OLDTHRESH
call time_II_operations1()
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
nbast = Sx(1)%nrow
call initIntegralOutputDims(setting%output,nbast,nbast,1,1,3*natoms)
CALL ls_getIntegrals(AORdefault,AORdefault,AOempty,AOempty,&
     &OverlapOperator,GeoDerivLHSSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,Sx,setting%IntegralTransformGC)
CALL ls_freeDmatFromSetting(setting)
call time_II_operations2(JOB_II_get_geoderivOverlap)

END SUBROUTINE II_get_geoderivOverlap

!> \brief Calculates the geometric derivative of a rectangular overlap matrix
!>        for DEBUG purposes only.
!> \author Patrick Merlot
!> \date 2013-09-13
!> \param Sx The 3*Natoms geometric derivative components of the overlap matrix
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_geoderivOverlap_mixed(Sx,natoms,AO1,AO2,setting,lupri,luerr)
   IMPLICIT NONE
   TYPE(LSSETTING),intent(INOUT) :: SETTING
   Integer,intent(IN)            :: lupri,luerr,AO1,AO2
   Type(matrix),intent(INOUT)    :: Sx(3*natoms)
   !
   integer             :: nAtoms,n1,n2

   SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
   n1 = Sx(1)%nrow
   n2 = Sx(1)%ncol
   call initIntegralOutputDims(setting%output,n1,n2,1,1,3*nAtoms)
   CALL ls_getIntegrals(AO1,AO2,AOempty,AOempty,&
        &OverlapOperator,GeoDerivLHSSpec,ContractedInttype,SETTING,LUPRI,LUERR)
   CALL retrieve_Output(lupri,setting,Sx,setting%IntegralTransformGC)
   CALL ls_freeDmatFromSetting(setting)
END SUBROUTINE II_get_geoderivOverlap_mixed

!> \brief Calculates the LHS half geometric derivative of the overlap matrix
!> \author T. Kjaergaard
!> \date 2011-08-22
!> \param Sx The 3*Natoms geometric derivative components of the overlap matrix
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_geoderivOverlapL(Sx,natoms,setting,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Integer,intent(IN)            :: lupri,luerr,natoms
Type(matrix),intent(INOUT)    :: Sx(3*natoms) !Sx,Sy,Sz for each atom
!
integer             :: nbast
real(realk)         :: OLDTHRESH
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
nbast = Sx(1)%nrow
call initIntegralOutputDims(setting%output,nbast,1,nbast,1,3*natoms)
CALL ls_getIntegrals(AORdefault,AOempty,AORdefault,AOempty,&
     &OverlapOperator,GeoDerivLHSSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,Sx,setting%IntegralTransformGC)
CALL ls_freeDmatFromSetting(setting)
END SUBROUTINE II_get_geoderivOverlapL

!> \brief Calculates the geometric derivative of the overlap matrix
!> \author T. Kjaergaard
!> \date 2011-08-22
!> \param Sx The 3*Natoms geometric derivative components of the overlap matrix
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_geoderivExchange(Kx,Dmat,natoms,setting,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Integer,intent(IN)            :: lupri,luerr,natoms
Type(matrix),intent(IN)       :: Dmat
Type(matrix),intent(INOUT)    :: Kx(3*natoms) !Kx,Ky,Kz for each atom
!
Type(matrix)        :: Dmat_AO
integer   :: Oper
integer             :: nbast,ndmats
real(realk)         :: OLDTHRESH,TS,TE
real(realk),pointer :: DfullRHS(:,:,:)
call time_II_operations1()
CALL LSTIMER('START ',TS,TE,LUPRI)
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%K_THR
nbast = Dmat%nrow
ndmats=1
IF(matrix_type .EQ. mtype_unres_dense)call lsquit('unres geoderiv K not testet- not implemented',-1)

IF(setting%IntegralTransformGC)THEN
   CALL mat_init(Dmat_AO,Dmat%nrow,Dmat%ncol)
   call GCAO2AO_transform_matrixD2(Dmat,Dmat_AO,setting,lupri)
   CALL ls_attachDmatToSetting(Dmat_AO,ndmats,setting,'RHS',2,4,.FALSE.,lupri)
ELSE
   CALL ls_attachDmatToSetting(Dmat,ndmats,setting,'RHS',2,4,.FALSE.,lupri)
ENDIF

IF (SETTING%SCHEME%CAM) THEN
  Oper = CAMOperator       !Coulomb attenuated method
ELSEIF (SETTING%SCHEME%SR_EXCHANGE) THEN
  Oper = ErfcOperator      !Short-Range Coulomb screened exchange
ELSE
  Oper = CoulombOperator   !Regular Coulomb metric 
ENDIF
!Calculates the HF-exchange contribution
call initIntegralOutputDims(setting%Output,nbast,nbast,1,1,3*Natoms)
!call ls_get_exchange_mat_serial(AORdefault,AORdefault,AORdefault,AORdefault,&
!     &                   Oper,GeoDerivLHSSpec,ContractedInttype,SETTING,LUPRI,LUERR)
call ls_get_exchange_mat(AORdefault,AORdefault,AORdefault,AORdefault,&
     &                   Oper,GeoDerivLHSSpec,ContractedInttype,SETTING,LUPRI,LUERR)
IF(setting%IntegralTransformGC) CALL mat_free(Dmat_AO)
CALL retrieve_Output(lupri,setting,Kx,setting%IntegralTransformGC)
CALL ls_freeDmatFromSetting(setting)
CALL LSTIMER('geoKbuild',TS,TE,LUPRI)
call time_II_operations2(JOB_II_get_geoderivExchange)

END SUBROUTINE II_get_geoderivExchange

!> \brief Calculates the geometric derivative of the overlap matrix
!> \author T. Kjaergaard
!> \date 2011-08-22
!> \param Sx The 3*Natoms geometric derivative components of the overlap matrix
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_geoderivCoulomb(Jx,Dmat,natoms,setting,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING),intent(INOUT) :: SETTING
Integer,intent(IN)            :: lupri,luerr,natoms
Type(matrix),intent(IN)       :: Dmat
Type(matrix),intent(INOUT)    :: Jx(3*natoms) !Kx,Ky,Kz for each atom
!
Type(matrix)        :: Dmat_AO
integer   :: Oper
integer             :: nbast,ndmats
real(realk)         :: OLDTHRESH,TS,TE
real(realk),pointer :: DfullRHS(:,:,:)
call time_II_operations1()
CALL LSTIMER('START ',TS,TE,LUPRI)
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%J_THR
nbast = Dmat%nrow
ndmats=1
IF(matrix_type .EQ. mtype_unres_dense)call lsquit('unres geoderiv K not testet- not implemented',-1)

IF(setting%IntegralTransformGC)THEN
   CALL mat_init(Dmat_AO,Dmat%nrow,Dmat%ncol)
   call GCAO2AO_transform_matrixD2(Dmat,Dmat_AO,setting,lupri)
   CALL ls_attachDmatToSetting(Dmat_AO,ndmats,setting,'RHS',3,4,.FALSE.,lupri)
ELSE
   CALL ls_attachDmatToSetting(Dmat,ndmats,setting,'RHS',3,4,.FALSE.,lupri)
ENDIF

call initIntegralOutputDims(setting%Output,nbast,nbast,1,1,3*Natoms)

CALL ls_getIntegrals(AORdefault,AORdefault,AORdefault,AORdefault,&
     &CoulombOperator,GeoDerivCoulombSpec,ContractedInttype,SETTING,LUPRI,LUERR)

IF(setting%IntegralTransformGC) CALL mat_free(Dmat_AO)
CALL retrieve_Output(lupri,setting,Jx,setting%IntegralTransformGC)
CALL ls_freeDmatFromSetting(setting)
CALL LSTIMER('geoJbuild',TS,TE,LUPRI)
call time_II_operations2(JOB_II_get_geoderivCoulomb)

END SUBROUTINE II_get_geoderivCoulomb

!> \brief Check symetry and remove all anti-symmetric components (which will not contribute) to the Coulomb ERI gradient
!> \author S. Reine
!> \date 2010
!> \param DmatLHS The left-hand-side density matrix
!> \param DmatRHS The right-hand-side density matrix
!> \param ndlhs The number of LHS density matrices
!> \param ndrhs The number of RHS density matrices
!> \param DLHS The output left-hand-side density matrix
!> \param DRHS The output right-hand-side density matrix
!> \param symLHS The symmetry of LHS density matrices
!> \param symRHS The symmetry of RHS density matrices
!> \param nlhs The output number of LHS density matrices
!> \param nrhs The output number of RHS density matrices
SUBROUTINE II_symmetrize_dmats(DmatLHS,DmatRHS,ndlhs,ndrhs,DLHS,DRHS,symLHS,symRHS,nlhs,nrhs)
implicit none
Integer,intent(in)          :: ndlhs,ndrhs
Type(matrixp),intent(IN)    :: DmatLHS(ndlhs),DmatRHS(ndrhs)
Integer,intent(inout)       :: symLHS(ndlhs),symRHS(ndrhs)
Integer,intent(inout)       :: nlhs,nrhs
Type(matrixp),intent(inout) :: DLHS(ndlhs),DRHS(ndrhs)
!
Integer                :: idmat,isym

IF (ndlhs.NE.ndrhs) CALL lsquit('Erorr in II_symmetrize_dmats: ndlhs different from ndrhs.',-1)

nlhs = 0
DO idmat=1,ndlhs
  isym = mat_get_isym(DmatLHS(idmat)%p)
  symLHS(idmat) = isym
  If (isym.EQ. 1) THEN
    nlhs = nlhs + 1
    DLHS(nlhs)%p => DmatLHS(idmat)%p
  ELSE IF  (isym.EQ. 3) THEN
    nlhs = nlhs + 1
    ALLOCATE (DLHS(nlhs)%p)
    call mat_init(DLHS(nlhs)%p,DmatLHS(idmat)%p%nrow,DmatLHS(idmat)%p%ncol)
    call mat_assign(DLHS(nlhs)%p,DmatLHS(idmat)%p)
    call util_get_symm_part(DLHS(nlhs)%p)
  ELSE IF ((isym.EQ. 2).OR.(isym.EQ. 4)) THEN
    ! Do nothing (i.e. remove the component from the calculation)
  ELSE
    call lsquit('Error in II_symmetrize_dmats. isym wrong in lhs-loop!',-1)
  ENDIF
ENDDO

nrhs = 0
DO idmat=1,ndrhs
  isym = mat_get_isym(DmatRHS(idmat)%p)
  symRHS(idmat) = isym
  If (isym.EQ. 1) THEN
    nrhs = nrhs + 1
    DRHS(nrhs)%p => DmatRHS(idmat)%p
  ELSE IF  (isym.EQ. 3) THEN
    nrhs = nrhs + 1
    ALLOCATE (DRHS(nrhs)%p)
    call mat_init(DRHS(nrhs)%p,DmatRHS(idmat)%p%nrow,DmatRHS(idmat)%p%ncol)
    call mat_assign(DRHS(nrhs)%p,DmatRHS(idmat)%p)
    call util_get_symm_part(DRHS(nrhs)%p)
  ELSE IF ((isym.EQ. 2).OR.(isym.EQ. 4)) THEN
    ! Do nothing (i.e. remove the component from the calculation)
  ELSE
    call lsquit('Error in II_symmetrize_dmats. isym wrong in rhs-loop!',-1)
  ENDIF
ENDDO

IF (nlhs.NE.nrhs) THEN
 call lsquit('Error in II_symmetrize_dmats. Different number of lhs and rhs density matrices.',-1)
ENDIF
!DO idmat=1,ndrhs
!  IF (symLHS(idmat).NE.symRHS(idmat)) &
! &  call lsquit('Erorr in II_symmetrize_dmats. Different matrix symmetries not implemented.',-1)
!ENDDO

END SUBROUTINE II_symmetrize_dmats

!> \brief free the symmetrized dmats
!> \author S. Reine
!> \date 2010
!> \param DLHS The left-hand-side density matrix
!> \param DRHS The right-hand-side density matrix
!> \param symLHS The symmetry of LHS density matrices
!> \param symRHS The symmetry of RHS density matrices
!> \param nlhs The output number of LHS density matrices
!> \param nrhs The output number of RHS density matrices
SUBROUTINE II_free_symmetrized_dmats(DLHS,DRHS,symLHS,symRHS,ndlhs,ndrhs)
implicit none
Integer,intent(IN)          :: ndlhs,ndrhs
Integer,intent(IN)          :: symLHS(ndlhs),symRHS(ndrhs)
Type(matrixp),intent(inout) :: DLHS(ndlhs),DRHS(ndrhs)
!
Integer                :: nlhs,nrhs,idmat,isym

nlhs = 0
DO idmat=1,ndlhs
  isym = symLHS(idmat)
  If (isym.EQ. 1) THEN
    nlhs = nlhs + 1
  ELSE IF  (isym.EQ. 3) THEN
    nlhs = nlhs + 1
    call mat_free(DLHS(nlhs)%p)
    DEALLOCATE(DLHS(nlhs)%p)
  ENDIF
ENDDO

nrhs = 0
DO idmat=1,ndrhs
  isym = symRHS(idmat) 
  If (isym.EQ. 1) THEN
    nrhs = nrhs + 1
  ELSE IF  (isym.EQ. 3) THEN
    nrhs = nrhs + 1
    call mat_free(DRHS(nrhs)%p)
    DEALLOCATE (DRHS(nrhs)%p)
  ENDIF
ENDDO

END SUBROUTINE II_free_symmetrized_dmats

!> \brief split the dmats
!> \author S. Reine
!> \date 2010
!> \param DmatLHS The left-hand-side density matrix
!> \param DmatRHS The right-hand-side density matrix
!> \param ndlhs The number of LHS density matrices
!> \param ndrhs The number of RHS density matrices
!> \param DLHS The output left-hand-side density matrix
!> \param DRHS The output right-hand-side density matrix
!> \param symLHS The symmetry of LHS density matrices
!> \param symRHS The symmetry of RHS density matrices
!> \param nlhs The output number of LHS density matrices
!> \param nrhs The output number of RHS density matrices
SUBROUTINE II_split_dmats(DmatLHS,DmatRHS,ndlhs,ndrhs,DLHS,DRHS,symLHS,symRHS,nlhs,nrhs)
  implicit none
  Integer,intent(in)          :: ndlhs,ndrhs
  Type(matrixp),intent(IN)    :: DmatLHS(ndlhs),DmatRHS(ndrhs)
  Integer,intent(inout)       :: symLHS(ndlhs),symRHS(ndrhs)
  Integer,intent(inout)       :: nlhs,nrhs
  Type(matrixp),intent(inout) :: DLHS(2*ndlhs),DRHS(2*ndrhs)
  !
  Integer                :: idmat,isym, ILHS, IRHS

  IF (ndlhs.NE.ndrhs) CALL lsquit('Error in II_split_dmats: ndlhs different from ndrhs.',-1)

  ! Make symmetry check
  do idmat=1,ndlhs
     ILHS = mat_get_isym(DmatLHS(idmat)%p)
     IRHS = mat_get_isym(DmatRHS(idmat)%p)

     if(ILHS==IRHS) then ! Set symmetries to those found by mat_get_isym
        symLHS(idmat) = ILHS
        symRHS(idmat) = IRHS
     else
        ! If symmetries are different, set them both to have no symmetry (=3)
        ! to make all special cases run through.
        print '(a,2i4)', 'WARNING: Different symmetries in II_split_dmats:', ILHS,IRHS
        symLHS(idmat) = 3
        symRHS(idmat) = 3
     end if

  end do

  nlhs = 0
  DO idmat=1,ndlhs
     isym = symLHS(idmat)
     If (isym.EQ. 1) THEN
        nlhs = nlhs + 1
        DLHS(nlhs)%p => DmatLHS(idmat)%p
     ELSE IF  (isym.EQ. 3) THEN
        nlhs = nlhs + 2
        ALLOCATE (DLHS(nlhs-1)%p)
        ALLOCATE (DLHS(nlhs)%p)
        call mat_init(DLHS(nlhs-1)%p,DmatLHS(idmat)%p%nrow,DmatLHS(idmat)%p%ncol)
        call mat_init(DLHS(nlhs)%p,DmatLHS(idmat)%p%nrow,DmatLHS(idmat)%p%ncol)
        call mat_assign(DLHS(nlhs-1)%p,DmatLHS(idmat)%p)
        call util_get_symm_part(DLHS(nlhs-1)%p)
        call util_get_antisymm_part(DmatLHS(idmat)%p,DLHS(nlhs)%p)
     ELSE IF ((isym.EQ. 2).OR.(isym.EQ. 4)) THEN
        ! Do nothing (i.e. remove the component from the calculation)
     ELSE
        call lsquit('Error in II_split_dmats. isym wrong in lhs-loop!',-1)
     ENDIF
  ENDDO

  nrhs = 0
  DO idmat=1,ndrhs
     isym = symRHS(idmat)
     If (isym.EQ. 1) THEN
        nrhs = nrhs + 1
        DRHS(nrhs)%p => DmatRHS(idmat)%p
     ELSE IF  (isym.EQ. 3) THEN
        nrhs = nrhs + 2
        ALLOCATE (DRHS(nrhs-1)%p)
        ALLOCATE (DRHS(nrhs)%p)
        call mat_init(DRHS(nrhs-1)%p,DmatRHS(idmat)%p%nrow,DmatRHS(idmat)%p%ncol)
        call mat_init(DRHS(nrhs)%p,DmatRHS(idmat)%p%nrow,DmatRHS(idmat)%p%ncol)
        call mat_assign(DRHS(nrhs-1)%p,DmatRHS(idmat)%p)
        call util_get_symm_part(DRHS(nrhs-1)%p)
        call util_get_antisymm_part(DmatRHS(idmat)%p,DRHS(nrhs)%p)
     ELSE IF ((isym.EQ. 2).OR.(isym.EQ. 4)) THEN
        ! Do nothing (i.e. remove the component from the calculation)
     ELSE
        call lsquit('II_split_dmats in II_symmetrize_dmats. isym wrong in rhs-loop!',-1)
     ENDIF
  ENDDO


END SUBROUTINE II_split_dmats

!> \brief free the split dmats
!> \author S. Reine
!> \date 2010
!> \param DLHS The left-hand-side density matrix
!> \param DRHS The right-hand-side density matrix
!> \param symLHS The symmetry of LHS density matrices
!> \param symRHS The symmetry of RHS density matrices
!> \param nlhs The output number of LHS density matrices
!> \param nrhs The output number of RHS density matrices
SUBROUTINE II_free_split_dmats(DLHS,DRHS,symLHS,symRHS,ndlhs,ndrhs)
implicit none
Integer,intent(IN)          :: ndlhs,ndrhs
Integer,intent(IN)          :: symLHS(ndlhs),symRHS(ndrhs)
Type(matrixp),intent(inout) :: DLHS(2*ndlhs),DRHS(2*ndrhs)
!
Integer                :: nlhs,nrhs,idmat,isym

nlhs = 0
DO idmat=1,ndlhs
  isym = symLHS(idmat)
  If (isym.EQ. 1) THEN
    nlhs = nlhs + 1
  ELSE IF  (isym.EQ. 3) THEN
    nlhs = nlhs + 2
    call mat_free(DLHS(nlhs-1)%p)
    call mat_free(DLHS(nlhs)%p)
    DEALLOCATE(DLHS(nlhs-1)%p)
    DEALLOCATE(DLHS(nlhs)%p)
  ENDIF
ENDDO

nrhs = 0
DO idmat=1,ndrhs
  isym = symRHS(idmat) 
  If (isym.EQ. 1) THEN
    nrhs = nrhs + 1
  ELSE IF  (isym.EQ. 3) THEN
    nrhs = nrhs + 2
    call mat_free(DRHS(nrhs-1)%p)
    call mat_free(DRHS(nrhs)%p)
    DEALLOCATE (DRHS(nrhs-1)%p)
    DEALLOCATE (DRHS(nrhs)%p)
  ENDIF
ENDDO

END SUBROUTINE II_free_split_dmats

!> \brief split the dmats
!> \author S. Reine
!> \date 2010
!> \param DmatLHS The left-hand-side density matrix
!> \param DmatRHS The right-hand-side density matrix
!> \param ndlhs The number of LHS density matrices
!> \param ndrhs The number of RHS density matrices
!> \param DLHS The output left-hand-side density matrix
!> \param DRHS The output right-hand-side density matrix
!> \param symLHS The symmetry of LHS density matrices
!> \param symRHS The symmetry of RHS density matrices
!> \param nlhs The output number of LHS density matrices
!> \param nrhs The output number of RHS density matrices
SUBROUTINE II_split_dmatsfull(DmatLHS,DmatRHS,nbast,ndlhs,ndrhs,DLHS,DRHS,&
     & symLHS,symRHS,nlhs,nrhs)
  implicit none
  Integer,intent(in)          :: ndlhs,ndrhs,nbast
  real(realk),intent(IN)      :: DmatLHS(nbast,nbast,ndlhs),DmatRHS(nbast,nbast,ndrhs)
  Integer,intent(inout)       :: symLHS(ndlhs),symRHS(ndrhs)
  Integer,intent(inout)       :: nlhs,nrhs
  real(realk),pointer         :: DLHS(:,:,:),DRHS(:,:,:)
  !
  Integer                :: idmat,isym, ILHS, IRHS

  IF (ndlhs.NE.ndrhs) CALL lsquit('Error in II_split_dmats: ndlhs different from ndrhs.',-1)

  ! Make symmetry check
  do idmat=1,ndlhs
     ILHS = matfull_get_isym(DmatLHS(:,:,idmat),nbast,nbast)
     IRHS = matfull_get_isym(DmatRHS(:,:,idmat),nbast,nbast)
     if(ILHS==IRHS) then ! Set symmetries to those found by mat_get_isym
        symLHS(idmat) = ILHS
        symRHS(idmat) = IRHS
     else
        ! If symmetries are different, set them both to have no symmetry (=3)
        ! to make all special cases run through.
        print '(a,2i4)', 'WARNING: Different symmetries in II_split_dmats:', ILHS,IRHS
        symLHS(idmat) = 3
        symRHS(idmat) = 3
     end if
  end do

  nlhs = 0
  DO idmat=1,ndlhs
     isym = symLHS(idmat)
     If (isym.EQ. 1) THEN
        nlhs = nlhs + 1
     ELSE IF  (isym.EQ. 3) THEN
        nlhs = nlhs + 2
     ELSE IF ((isym.EQ. 2).OR.(isym.EQ. 4)) THEN
        ! Do nothing (i.e. remove the component from the calculation)
     ELSE
        call lsquit('Error in II_split_dmats. isym wrong in lhs-loop!',-1)
     ENDIF
  ENDDO
  IF (nlhs.EQ. 0) RETURN

  call mem_alloc(DLHS,nbast,nbast,nlhs)
  nlhs = 0
  DO idmat=1,ndlhs
     isym = symLHS(idmat)
     If (isym.EQ. 1) THEN
        nlhs = nlhs + 1
        call dcopy(nbast*nbast,DmatLHS(:,:,idmat),1,DLHS(:,:,nlhs),1)
     ELSE IF  (isym.EQ. 3) THEN
        nlhs = nlhs + 2
        ! Asym  DLHS(:,:,nlhs-1)
        ! AntiSym  DLHS(:,:,nlhs)
        call util_get_symm_And_antisymm_part_full(DmatLHS(:,:,idmat),DLHS(:,:,nlhs-1),DLHS(:,:,nlhs),nbast)
     ELSE IF ((isym.EQ. 2).OR.(isym.EQ. 4)) THEN
        ! Do nothing (i.e. remove the component from the calculation)
     ELSE
        call lsquit('Error in II_split_dmats. isym wrong in lhs-loop!',-1)
     ENDIF
  ENDDO

  nrhs = 0
  DO idmat=1,ndrhs
     isym = symRHS(idmat)
     If (isym.EQ. 1) THEN
        nrhs = nrhs + 1
     ELSE IF  (isym.EQ. 3) THEN
        nrhs = nrhs + 2
     ELSE IF ((isym.EQ. 2).OR.(isym.EQ. 4)) THEN
        ! Do nothing (i.e. remove the component from the calculation)
     ELSE
        call lsquit('II_split_dmats in II_symmetrize_dmats. isym wrong in rhs-loop!',-1)
     ENDIF
  ENDDO
  call mem_alloc(DRHS,nbast,nbast,nrhs)

  nrhs = 0
  DO idmat=1,ndrhs
     isym = symRHS(idmat)
     If (isym.EQ. 1) THEN
        nrhs = nrhs + 1
        call dcopy(nbast*nbast,DmatRHS(:,:,idmat),1,DRHS(:,:,nrhs),1)
     ELSE IF  (isym.EQ. 3) THEN
        nrhs = nrhs + 2
        ! Asym  DRHS(:,:,nrhs-1)
        ! AntiSym  DRHS(:,:,nrhs)
        call util_get_symm_And_antisymm_part_full(DmatRHS(:,:,idmat),DRHS(:,:,nrhs-1),DRHS(:,:,nrhs),nbast)
     ELSE IF ((isym.EQ. 2).OR.(isym.EQ. 4)) THEN
        ! Do nothing (i.e. remove the component from the calculation)
     ELSE
        call lsquit('II_split_dmats in II_symmetrize_dmats. isym wrong in rhs-loop!',-1)
     ENDIF
  ENDDO

END SUBROUTINE II_split_dmatsfull

!> \brief free the split dmats
!> \author T. Kjaergaard
!> \date 2013
!> \param DLHS The left-hand-side density matrix
!> \param DRHS The right-hand-side density matrix
SUBROUTINE II_free_split_dmatsfull(DLHS,DRHS)
implicit none
real(realk),pointer         :: DLHS(:,:,:),DRHS(:,:,:)
call mem_dealloc(DLHS)
call mem_dealloc(DRHS)
END SUBROUTINE II_free_split_dmatsfull

!> \brief calculate ...
!> \author S. Reine
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param GGem 
!> \param nbast the number of basis functions
!> \param filedump 
SUBROUTINE II_get_GaussianGeminalFourCenter(LUPRI,LUERR,SETTING,GGem,nbast,filedump)
IMPLICIT NONE
Integer,intent(in)    :: nbast
Real(realk)           :: GGem(nbast,nbast,nbast,nbast,1)
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
logical               :: filedump
!
real         :: TS,TE
real(realk)         :: TS2,TE2
real(realk)         :: coeff(6),exponent(6),tmp
real(realk)         :: coeff2(21),sumexponent(21),prodexponent(21)
integer             :: iunit,i,j,k,l,IJ
integer             :: nGaussian,nG2

nGaussian = 6
nG2 = nGaussian*(nGaussian+1)/2

CALL LSTIMER('START ',TS2,TE2,LUPRI)
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR

GGem = 0E0_realk

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

call set_GGem(Setting%GGem,coeff,exponent,nGaussian)
call initIntegralOutputDims(setting%Output,nbast,nbast,nbast,nbast,1)
call cpu_time(TS)
CALL ls_getIntegrals(AORdefault,AORdefault,AORdefault,AORdefault,&
     & CoulombOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
call cpu_time(TE)
print*, ' cpu_time (Coulomb)', TE-TS
CALL retrieve_Output(lupri,setting,GGem,setting%IntegralTransformGC)
IF (filedump) THEN
  iunit = -1
  call LSOPEN(iunit,'coulomb.dat','UNKNOWN','FORMATTED')
  tmp = 0E0_realk
  DO I=1,nbast
    DO J=1,nbast
      DO K=1,nbast
        DO L=1,K
         if (abs(GGem(I,J,K,L,1)).GT. 1E-12_realk) write(iunit,'(4I5,1G24.10)') I,J,K,L,GGem(I,J,K,L,1)
         tmp = tmp + GGem(I,J,K,L,1)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  write(lupri,'(1X,A,G24.10)') 'QQQ:GGem sum:',tmp
  call LSCLOSE(iunit,'KEEP')
ENDIF

call initIntegralOutputDims(setting%Output,nbast,nbast,nbast,nbast,1)
call cpu_time(TS)
CALL ls_getIntegrals(AORdefault,AORdefault,AORdefault,AORdefault,&
     &               GGemOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
call cpu_time(TE)
print*, ' cpu_time (GGem000)', TE-TS
CALL retrieve_Output(lupri,setting,GGem,setting%IntegralTransformGC)
IF (filedump) THEN
  iunit = -1
  call LSOPEN(iunit,'ggem000.dat','UNKNOWN','FORMATTED')
  tmp = 0E0_realk
  DO I=1,nbast
    DO J=1,nbast
      DO K=1,nbast
        DO L=1,K
         if (abs(GGem(I,J,K,L,1)).GT. 1E-12_realk) write(iunit,'(4I5,1G24.10)') I,J,K,L,GGem(I,J,K,L,1)
         tmp = tmp + GGem(I,J,K,L,1)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  write(lupri,'(1X,A,G24.10)') 'QQQ:GGem sum:',tmp
  call LSCLOSE(iunit,'KEEP')
ENDIF

call initIntegralOutputDims(setting%Output,nbast,nbast,nbast,nbast,1)
call cpu_time(TS)
CALL ls_getIntegrals(AORdefault,AORdefault,AORdefault,AORdefault,&
     &               GGemCouOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
call cpu_time(TE)
print*, ' cpu_time (GGemCou)', TE-TS
CALL retrieve_Output(lupri,setting,GGem,setting%IntegralTransformGC)
IF (filedump) THEN
  iunit = -1
  call LSOPEN(iunit,'ggemcou.dat','UNKNOWN','FORMATTED')
  tmp = 0E0_realk
  DO I=1,nbast
    DO J=1,nbast
      DO K=1,nbast
        DO L=1,K
         if (abs(GGem(I,J,K,L,1)).GT. 1E-12_realk) write(iunit,'(4I5,1G24.10)') I,J,K,L,GGem(I,J,K,L,1)
         tmp = tmp + GGem(I,J,K,L,1)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  write(lupri,'(1X,A,G24.10)') 'QQQ:GGem sum:',tmp
  call LSCLOSE(iunit,'KEEP')
ENDIF

call set_GGem(Setting%GGem,coeff2,sumexponent,prodexponent,nG2)

call initIntegralOutputDims(setting%Output,nbast,nbast,nbast,nbast,1)
call cpu_time(TS)
CALL ls_getIntegrals(AORdefault,AORdefault,AORdefault,AORdefault,&
     &               GGemOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
call cpu_time(TE)
print*, ' cpu_time (GGemQua)', TE-TS
CALL retrieve_Output(lupri,setting,GGem,setting%IntegralTransformGC)
IF (filedump) THEN
  iunit = -1
  call LSOPEN(iunit,'ggemqua.dat','UNKNOWN','FORMATTED')
  tmp = 0E0_realk
  DO I=1,nbast
    DO J=1,nbast
      DO K=1,nbast
        DO L=1,K
         if (abs(GGem(I,J,K,L,1)).GT. 1E-12_realk) write(iunit,'(4I5,1G24.10)') I,J,K,L,GGem(I,J,K,L,1)
         tmp = tmp + GGem(I,J,K,L,1)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  write(lupri,'(1X,A,G24.10)') 'QQQ:GGem sum:',tmp
  call LSCLOSE(iunit,'KEEP')
ENDIF

call initIntegralOutputDims(setting%Output,nbast,nbast,nbast,nbast,1)
call cpu_time(TS)
CALL ls_getIntegrals(AORdefault,AORdefault,AORdefault,AORdefault,&
     &               GGemGrdOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
call cpu_time(TE)
print*, ' cpu_time (GGemGrd)', TE-TS
CALL retrieve_Output(lupri,setting,GGem,setting%IntegralTransformGC)
IF (filedump) THEN
  iunit = -1
  call LSOPEN(iunit,'ggemgrd.dat','UNKNOWN','FORMATTED')
  tmp = 0E0_realk
  DO I=1,nbast
    DO J=1,nbast
      DO K=1,nbast
        DO L=1,K
         if (abs(GGem(I,J,K,L,1)).GT. 1E-12_realk) write(iunit,'(4I5,1G24.10)') I,J,K,L,GGem(I,J,K,L,1)
         tmp = tmp + GGem(I,J,K,L,1)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  write(lupri,'(1X,A,G24.10)') 'QQQ:GGem sum:',tmp
  call LSCLOSE(iunit,'KEEP')
ENDIF
call free_GGem(Setting%GGem)

CALL LSTIMER('GGem   ',TS2,TE2,LUPRI)
END SUBROUTINE II_get_GaussianGeminalFourCenter

!> \brief Calculates the magnetic derivative 4 center 2 electron repulsion integrals
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param nbast The number of basis functions
SUBROUTINE II_get_magderiv_4center_eri(LUPRI,LUERR,SETTING,nbast,Dmat,Fx)
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,nbast
Real(realk),pointer   :: integrals(:,:,:,:,:)
INTEGER               :: a,b,c,d,X
type(matrix)          :: Dmat,Fx(3),Dmat_AO
REAL(REALK),pointer   :: Fxfull(:,:,:),Dfull(:,:) 
real(realk) :: KFAC,JFAC,maxCoor,OLDTHRESH
Real(realk),pointer   :: integrals2(:,:,:,:,:)
call time_II_operations1()

!The size of the magnetic derivative integral will increase with 
!the size of the system so we modify the threshold with the 
!largest X,Y or Z distance in the molecule. 
call determine_maxCoor(SETTING%MOLECULE(1)%p,maxCoor)
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%J_THR*(1.0E0_realk/maxCoor)
!derivative on the LHS 
call initIntegralOutputDims(setting%output,nbast,nbast,nbast,nbast,3)
CALL ls_getIntegrals(AORdefault,AORdefault,AORdefault,AORdefault,&
     &CoulombOperator,MagDerivSpec,ContractedInttype,SETTING,LUPRI,LUERR)
KFAC = setting%scheme%exchangeFactor * 0.5E0_realk
JFAC = 2.0E0_realk * 0.5E0_realk
call mem_alloc(integrals,nbast,nbast,nbast,nbast,3)
CALL retrieve_Output(lupri,setting,integrals,.FALSE.)

!WRITE(lupri,*)'ls_getIntegrals  MagDerivRSpec   RHS' 
!!derivative on the RHS 
!call initIntegralOutputDims(setting%output,nbast,nbast,nbast,nbast,3)
!CALL ls_getIntegrals(AORdefault,AORdefault,AORdefault,AORdefault,&
!     &CoulombOperator,MagDerivRSpec,ContractedInttype,SETTING,LUPRI,LUERR)
!call mem_alloc(integrals2,nbast,nbast,nbast,nbast,3)
!CALL retrieve_Output(lupri,setting,integrals2,.FALSE.)
!DO X=1,3
! DO b=1,nbast
!  DO a=1,nbast
!   DO d=1,nbast
!    DO c=1,nbast
!     WRITE(lupri,'(A,I2,A,I2,A,I2,A,I2,A,I2,A,ES20.11)')'Integrals2(',a,',',b,',',c,',',d,',',X,')=',Integrals2(a,b,c,d,X)
!    ENDDO
!   ENDDO
!  ENDDO
! ENDDO
!ENDDO
!DO X=1,3
! DO b=1,nbast
!  DO a=1,nbast
!   DO d=1,nbast
!    DO c=1,nbast
!     IF(ABS(Integrals2(a,b,c,d,X)-Integrals(c,d,a,b,X)).GT.1.0E-10)THEN
!        print*,'a,b,c,d,X',a,b,c,d,X
!        print*,'Integrals2(a,b,c,d,X)',Integrals2(a,b,c,d,X)
!        print*,'Integrals(c,d,a,b,X) ',Integrals(c,d,a,b,X)
!        call lsquit('test done',-1)
!     ENDIF
!    ENDDO
!   ENDDO
!  ENDDO
! ENDDO
!ENDDO
!call mem_dealloc(integrals2)

call mem_alloc(integrals2,nbast,nbast,nbast,nbast,3)
DO X=1,3
   DO b=1,nbast
      DO a=1,nbast
         DO d=1,nbast
            DO c=1,nbast
               Integrals2(a,b,c,d,X) = Integrals(c,d,a,b,X) 
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
DO X=1,3
   DO b=1,nbast
      DO a=1,nbast
         DO d=1,nbast
            DO c=1,nbast
               Integrals(a,b,c,d,X) = Integrals2(a,b,c,d,X)+Integrals(a,b,c,d,X)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

call mem_alloc(Dfull,nbast,nbast)
IF(setting%IntegralTransformGC)THEN
   CALL mat_init(Dmat_AO,Dmat%nrow,Dmat%ncol)
   call GCAO2AO_transform_matrixD2(Dmat,Dmat_AO,setting,lupri)
   CALL DCOPY(nbast*nbast,Dmat_AO%elms,1,Dfull,1)
   CALL mat_free(Dmat_AO)
ELSE
   CALL DCOPY(nbast*nbast,Dmat%elms,1,Dfull,1)
ENDIF
call mem_alloc(Fxfull,nbast,nbast,3)
call ls_dzero(Fxfull,nbast*nbast*3)
do X=1,3
DO a=1,nbast 
   DO b=1,nbast 
      DO c=1,nbast 
         DO d=1,nbast 
            Fxfull(a,b,X)=Fxfull(a,b,X)+Dfull(d,c)*& 
                 & (JFAC*Integrals(a,b,c,d,X) & 
                 & -KFAC*Integrals(a,d,c,b,X) ) 
         ENDDO
      ENDDO
   ENDDO
ENDDO
enddo
do X=1,3
   CALL DCOPY(nbast*nbast,Fxfull(:,:,X),1,Fx(X)%elms,1)
   IF(setting%IntegralTransformGC)THEN
      call AO2GCAO_transform_matrixF(Fx(X),setting,lupri)
   ENDIF
enddo
call mem_dealloc(integrals)
call mem_dealloc(integrals2)
call mem_dealloc(Dfull)
call mem_dealloc(Fxfull)
call time_II_operations2(JOB_II_get_magderiv_4center_eri)

END SUBROUTINE II_get_magderiv_4center_eri

!> \brief Calculates the magnetic derivative Exchange matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param nbast The number of basis functions
SUBROUTINE II_get_magderivF(LUPRI,LUERR,SETTING,nbast,Dmat,Fx)
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,nbast
type(matrix)          :: Dmat(1),Fx(3),temp(3)
!
integer :: I

call II_get_magderivJ(LUPRI,LUERR,SETTING,nbast,Dmat,Fx)

IF (SETTING%SCHEME%exchangeFactor.GT. 0.0E0_realk)THEN
   call mat_init(temp(1),nbast,nbast)
   call mat_init(temp(2),nbast,nbast)
   call mat_init(temp(3),nbast,nbast)
   call II_get_magderivK(LUPRI,LUERR,SETTING,nbast,Dmat,temp)
   DO I =1,3
      call mat_daxpy(1.0E0_realk,temp(I),Fx(I))
      call mat_free(temp(I))
   ENDDO
ENDIF
END SUBROUTINE II_get_magderivF

!> \brief Calculates the magnetic derivative Exchange matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param nbast The number of basis functions
SUBROUTINE II_get_magderivK(LUPRI,LUERR,SETTING,nbast,Dmat,Kx)
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,nbast
type(matrix)          :: Dmat(1),Kx(3)
!
type(matrix)          :: tempm1,DMAT_AO
type(matrix),pointer  :: D2(:),Kx2(:)
real(realk) :: KFAC,maxCoor,OLDTHRESH,TS,TE
Real(realk),pointer   :: DFULLRHS(:,:,:)
integer      :: Oper,isym,ndmat2,idmat,i
logical,pointer :: symmetricD(:)
IF (SETTING%SCHEME%exchangeFactor.EQ. 0.0E0_realk)RETURN
call time_II_operations1()
CALL LSTIMER('START ',TS,TE,LUPRI)

ISYM = mat_get_isym(Dmat(1))
IF(ISYM.EQ.1)THEN !sym
   ndmat2 = 1
   call mem_alloc(symmetricD,1)
   symmetricD(1)=.TRUE.
ELSEIF(ISYM.EQ.2)THEN !antisym
   ndmat2 = 1
   call mem_alloc(symmetricD,1)
   symmetricD(1)=.FALSE.
ELSEIF(ISYM.EQ.3)THEN !nonsym
   ndmat2 = 2 !we split up into sym and anti sym part
   call mem_alloc(symmetricD,2)
   symmetricD(1)=.TRUE.
   symmetricD(2)=.FALSE.
ELSE
   call mat_zero(Kx(1))
   call mat_zero(Kx(2))
   call mat_zero(Kx(3))
   RETURN
ENDIF

call mem_alloc(D2,ndmat2)
do idmat=1,ndmat2
   call mat_init(D2(idmat),Dmat(1)%nrow,Dmat(1)%ncol)
enddo
call mat_assign(D2(1),Dmat(1))
if(ndmat2.EQ.2)THEN
   call util_get_symm_part(D2(1))
   call util_get_antisymm_part(Dmat(1),D2(2))
endif
call mem_alloc(Kx2,ndmat2*3)
do idmat=1,ndmat2
   do i=1,3
      call mat_init(Kx2(i+(idmat-1)*3),Kx(1)%nrow,Kx(1)%ncol)
      call mat_zero(Kx2(i+(idmat-1)*3)) 
   enddo
enddo
call II_get_magderivK1(LUPRI,LUERR,SETTING,nbast,D2,Kx2,symmetricD,ndmat2)
do idmat=1,ndmat2
   IF(symmetricD(idmat))THEN
      call mat_init(tempm1,nbast,nbast)
      DO i=1,3
         call mat_assign(Kx(i),Kx2(i+(idmat-1)*3))
         call mat_trans(Kx(i),tempm1)
         call mat_daxpy(-1.0E0_realk,tempm1,Kx(i))
         call mat_scal(-0.5E0_realk,Kx(I))
      ENDDO
      call mat_free(tempm1)
   ELSE
      call mat_init(tempm1,nbast,nbast)
      DO i=1,3
         call mat_trans(Kx2(i+(idmat-1)*3),tempm1)
         call mat_daxpy(1.0E0_realk,tempm1,Kx2(i+(idmat-1)*3))
         call mat_daxpy(-0.5E0_realk,Kx2(i+(idmat-1)*3),Kx(i))
      ENDDO
      call mat_free(tempm1)
   ENDIF
ENDDO

do idmat=1,ndmat2
   call mat_free(D2(idmat))
   do i=1,3
      call mat_free(Kx2(i+(idmat-1)*3))
   enddo
enddo
call mem_dealloc(Kx2)
call mem_dealloc(D2)
call mem_dealloc(symmetricD)
CALL LSTIMER('magKbuild',TS,TE,LUPRI)
call time_II_operations2(JOB_II_get_magderivK_4center_eri)
END SUBROUTINE II_get_magderivK

!> \brief Calculates the magnetic derivative Exchange matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param nbast The number of basis functions
SUBROUTINE II_get_magderivK1(LUPRI,LUERR,SETTING,nbast,Dmat,Kx,symmetricD,ndmat)
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,nbast
type(matrix)          :: Dmat(ndmat),Kx(3*ndmat)
logical               :: symmetricD(nDmat)
!
type(matrix)          :: tempm1,DMAT_AO(ndmat)
INTEGER               :: ndmat,idmat
real(realk) :: KFAC,maxCoor,OLDTHRESH,TS,TE
Real(realk),pointer   :: DFULLRHS(:,:,:)
integer      :: Oper
IF (ABS(SETTING%SCHEME%exchangeFactor).LT.1.0E0-15_realk)RETURN
!The size of the magnetic derivative integral will increase with 
!the size of the system so we modify the threshold with the 
!largest X,Y or Z distance in the molecule. 
call determine_maxCoor(SETTING%MOLECULE(1)%p,maxCoor)
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%K_THR*(2.0E0_realk/maxCoor)
IF(matrix_type .EQ. mtype_unres_dense)call lsquit('unres magderiv K not testet- not implemented',-1)
IF(setting%IntegralTransformGC)THEN
   do idmat=1,ndmat
      CALL mat_init(Dmat_AO(idmat),Dmat(1)%nrow,Dmat(1)%ncol)
      call GCAO2AO_transform_matrixD2(Dmat(idmat),Dmat_AO(idmat),setting,lupri)
   enddo
   CALL ls_attachDmatToSetting(Dmat_AO,ndmat,setting,'RHS',2,4,.FALSE.,lupri)
ELSE
   CALL ls_attachDmatToSetting(Dmat,ndmat,setting,'RHS',2,4,.FALSE.,lupri)
ENDIF

IF (SETTING%SCHEME%CAM) THEN
  Oper = CAMOperator       !Coulomb attenuated method
ELSEIF (SETTING%SCHEME%SR_EXCHANGE) THEN
  Oper = ErfcOperator      !Short-Range Coulomb screened exchange
ELSE
  Oper = CoulombOperator   !Regular Coulomb metric 
ENDIF
!Calculates the HF-exchange contribution
call initIntegralOutputDims(setting%Output,nbast,nbast,1,1,3*ndmat)
call ls_get_exchange_mat(AORdefault,AORdefault,AORdefault,AORdefault,&
     &                   Oper,MagDerivSpec,ContractedInttype,SETTING,LUPRI,LUERR)
IF(setting%IntegralTransformGC)THEN
   do idmat=1,ndmat
      CALL mat_free(Dmat_AO(idmat))
   enddo
ENDIF
CALL retrieve_Output(lupri,setting,Kx,setting%IntegralTransformGC)
CALL ls_freeDmatFromSetting(setting)
END SUBROUTINE II_get_magderivK1

!> \brief Calculates the magnetic derivative Coulomb integrals
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param nbast The number of basis functions
SUBROUTINE II_get_magderivJ(LUPRI,LUERR,SETTING,nbast,Dmat,Jx)
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,nbast
type(matrix)          :: Dmat(1),Jx(3)
!
INTEGER               :: idmat,isym,ndmat
real(realk) :: TS,TE
type(matrix),pointer  :: Jx2(:)
type(matrix) :: D2(1)
call time_II_operations1()
IF (setting%scheme%densfit) THEN
   CALL II_get_df_magderivJ(LUPRI,LUERR,SETTING,nbast,Dmat,Jx)
ELSE
   CALL LSTIMER('START ',TS,TE,LUPRI)
   ndmat = 1
   ISYM = mat_get_isym(Dmat(1))
   IF(ISYM.EQ.1)THEN !sym
      !we can use the jengine where we only differentiate on the LHS
      !as the magderiv on the RHS is zero for sym D mat
      call II_get_magderivJ1(LUPRI,LUERR,SETTING,nbast,Dmat,Jx,ndmat)
      DO idmat=1,3
         call mat_scal(0.5E0_realk,Jx(Idmat))
      ENDDO
   ELSEIF(ISYM.EQ.2)THEN !antisym
      call II_get_magderivJ1asym(LUPRI,LUERR,SETTING,nbast,Dmat(1),Jx,ndmat)
   ELSEIF(ISYM.EQ.3)THEN !nonsym
      !we split up into sym and anti sym part
      !SYM PART
      call mat_init(D2(1),Dmat(1)%nrow,Dmat(1)%ncol)
      call mat_assign(D2(1),Dmat(1))
      call util_get_symm_part(D2(1))
      call II_get_magderivJ1(LUPRI,LUERR,SETTING,nbast,D2,Jx,ndmat)
      DO idmat=1,3
         call mat_scal(0.5E0_realk,Jx(Idmat))
      ENDDO
      !ASYM PART
      call util_get_antisymm_part(Dmat(1),D2(1))
      call mem_alloc(Jx2,3)
      do idmat=1,3
         call mat_init(Jx2(idmat),Jx(1)%nrow,Jx(1)%ncol)
      enddo
      call II_get_magderivJ1asym(LUPRI,LUERR,SETTING,nbast,D2,Jx2,ndmat)
      call mat_free(D2(1))
      DO idmat=1,3
         call mat_daxpy(-1.0E0_realk,Jx2(Idmat),Jx(idmat))
         call mat_free(Jx2(idmat))
      ENDDO
      call mem_dealloc(Jx2)
   ELSE
      !zero Dmat
      call mat_zero(Jx(1))
      call mat_zero(Jx(2))
      call mat_zero(Jx(3))
   ENDIF
   CALL LSTIMER('magJengin',TS,TE,LUPRI)
ENDIF
call time_II_operations2(JOB_II_get_magderivJ_4center_eri)
end SUBROUTINE II_get_magderivJ

!> \brief Calculates the magnetic derivative Coulomb integrals
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param nbast The number of basis functions
SUBROUTINE II_get_magderivJ1(LUPRI,LUERR,SETTING,nbast,Dmat,Jx,ndmat)
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,nbast,ndmat
type(matrix)          :: Dmat(ndmat),Jx(3*ndmat)
!
INTEGER               :: idmat
real(realk)           :: maxCoor,OLDTHRESH
type(matrix)          :: Dmat_AO(ndmat)

!The size of the magnetic derivative integral will increase with 
!the size of the system so we modify the threshold with the 
!largest X,Y or Z distance in the molecule. 
call determine_maxCoor(SETTING%MOLECULE(1)%p,maxCoor)
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%J_THR*(1.0E0_realk/maxCoor)

IF(setting%IntegralTransformGC)THEN
   do idmat=1,ndmat
      CALL mat_init(Dmat_AO(idmat),nbast,nbast)
      call GCAO2AO_transform_matrixD2(Dmat(idmat),Dmat_AO(idmat),setting,lupri)
   enddo
   CALL ls_attachDmatToSetting(Dmat_AO,ndmat,setting,'RHS',3,4,.TRUE.,lupri)
ELSE
   CALL ls_attachDmatToSetting(Dmat,ndmat,setting,'RHS',3,4,.TRUE.,lupri)
ENDIF
call initIntegralOutputDims(setting%Output,nbast,nbast,1,1,3*ndmat)
call ls_jengine(AORdefault,AORdefault,AORdefault,AORdefault,CoulombOperator,MagDerivSpec,&
     & ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,Jx,setting%IntegralTransformGC)
IF(SETTING%SCHEME%FMM) THEN
   call lsquit('not tested',-1)
!   call ls_jengineClassicalMAT
ENDIF

CALL ls_freeDmatFromSetting(setting)
IF(setting%IntegralTransformGC)THEN
   do idmat=1,ndmat
      CALL mat_free(Dmat_AO(idmat))
   enddo
ENDIF
END SUBROUTINE II_get_magderivJ1

!> \brief Calculates the magnetic derivative Coulomb integrals for asym matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param nbast The number of basis functions
SUBROUTINE II_get_magderivJ1asym(LUPRI,LUERR,SETTING,nbast,Dmat,Jx,ndmat)
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,nbast,ndmat
type(matrix)          :: Dmat(ndmat),Jx(3*ndmat)
!
INTEGER               :: idmat
real(realk)           :: maxCoor,OLDTHRESH
type(matrix)          :: Dmat_AO(ndmat)

IF(SETTING%SCHEME%FMM) call lsquit('not tested',-1)
!The size of the magnetic derivative integral will increase with 
!the size of the system so we modify the threshold with the 
!largest X,Y or Z distance in the molecule. 
call determine_maxCoor(SETTING%MOLECULE(1)%p,maxCoor)
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%J_THR*(1.0E0_realk/maxCoor)

IF(setting%IntegralTransformGC)THEN
   do idmat=1,ndmat
      CALL mat_init(Dmat_AO(idmat),nbast,nbast)
      call GCAO2AO_transform_matrixD2(Dmat(idmat),Dmat_AO(idmat),setting,lupri)
   enddo
   CALL ls_attachDmatToSetting(Dmat_AO,ndmat,setting,'RHS',3,4,.TRUE.,lupri)
ELSE
   CALL ls_attachDmatToSetting(Dmat,ndmat,setting,'RHS',3,4,.TRUE.,lupri)
ENDIF
call initIntegralOutputDims(setting%Output,nbast,nbast,1,1,3*ndmat)
call ls_get_coulomb_mat(AORdefault,AORdefault,AORdefault,AORdefault,CoulombOperator,MagDerivRSpec,&
     & ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,Jx,setting%IntegralTransformGC)

CALL ls_freeDmatFromSetting(setting)
IF(setting%IntegralTransformGC)THEN
   do idmat=1,ndmat
      CALL mat_free(Dmat_AO(idmat))
   enddo
ENDIF
END SUBROUTINE II_get_magderivJ1asym

!> \brief Calculates the Fock matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param F the exchange matrix
!> \param ndmat the number of density matrices
SUBROUTINE II_get_Econt(LUPRI,LUERR,SETTING,D,Econt,ndmat,modThresh)
IMPLICIT NONE
integer               :: ndmat,LUPRI,LUERR
TYPE(MATRIX)          :: D(ndmat)
real(realk)           :: Econt(ndmat)
TYPE(LSSETTING)       :: SETTING
real(realk)           :: modThresh
!
integer               :: idmat,nbast
type(matrix)          :: D_AO(ndmat)
type(matrix)          :: Ddiff2(ndmat)
real(realk) :: DeltaEcontdiff,twoElectronEnergyCont0
real(realk),parameter :: D2=2.0E0_realk
nbast = D(1)%nrow
IF(setting%IntegralTransformGC)THEN
   setting%IntegralTransformGC = .FALSE.
   do idmat=1,ndmat
      CALL mat_init(D_AO(idmat),nbast,nbast)
      call GCAO2AO_transform_matrixD2(D(idmat),D_AO(idmat),setting,lupri)
   enddo
   IF(SaveF0andD0)THEN
      do idmat=1,ndmat
         CALL mat_init(Ddiff2(idmat),nbast,nbast)
         call mat_assign(Ddiff2(idmat),D_AO(idmat))
         call mat_free(D_AO(idmat))
         call mat_daxpy(-1E0_realk,incrD0(1),Ddiff2(idmat)) !Ddiff is now the difference density 
      enddo
      call II_get_Econt1(LUPRI,LUERR,SETTING,Ddiff2,Econt,ndmat,modThresh)
      twoElectronEnergyCont0 = mat_dotproduct(incrD0(1),incrF0(1))
      write(lupri,'(A30,ES26.16)')'Tr(D0,F0) =',twoElectronEnergyCont0 
      do idmat=1,ndmat
         DeltaEcontdiff =  D2*mat_dotproduct(Ddiff2(idmat),incrF0(1))
         write(lupri,'(A30,I3,A2,ES26.16)')'Tr(Ddiff,F0), (idmat=',idmat,')=',DeltaEcontdiff
         write(lupri,'(A30,I3,A2,ES26.16)')'Tr(Ddiff,F(Ddiff), (idmat=',idmat,')=',Econt(idmat)
         CALL mat_free(Ddiff2(idmat))
         Econt(idmat) = Econt(idmat) + DeltaEcontdiff + twoElectronEnergyCont0
         write(lupri,'(A30,I3,A2,ES26.16)')'Econt(idmat=',idmat,')=',Econt(idmat)
      enddo
   ELSE
      twoElectronEnergyCont0 = 0.0E0_realk
      call II_get_Econt1(LUPRI,LUERR,SETTING,D_AO,Econt,ndmat,modThresh)
      do idmat=1,ndmat
         write(lupri,'(A30,I3,A2,ES26.16)')'Econt(idmat=',idmat,')=',Econt(idmat)
      enddo
      do idmat=1,ndmat
         CALL mat_free(D_AO(idmat))
      enddo
   ENDIF
   setting%IntegralTransformGC = .TRUE.
ELSE
   !already in AO
   IF(SaveF0andD0)THEN
      do idmat=1,ndmat
         CALL mat_init(Ddiff2(idmat),nbast,nbast)
         call mat_assign(Ddiff2(idmat),D(idmat))
         call mat_daxpy(-1E0_realk,incrD0(1),Ddiff2(idmat)) !Ddiff is now the difference density 
      enddo
      call II_get_Econt1(LUPRI,LUERR,SETTING,Ddiff2,Econt,ndmat,modThresh)   
      twoElectronEnergyCont0 = mat_dotproduct(incrD0(1),incrF0(1))
      write(lupri,'(A30,ES26.16)')'Tr(D0,F0) =',twoElectronEnergyCont0 
      do idmat=1,ndmat
         DeltaEcontdiff =  D2*mat_dotproduct(Ddiff2(idmat),incrF0(1))
         write(lupri,'(A30,I3,A2,ES26.16)')'Tr(Ddiff,F0), (idmat=',idmat,')=',DeltaEcontdiff
         write(lupri,'(A30,I3,A2,ES26.16)')'Tr(Ddiff,F(Ddiff), (idmat=',idmat,')=',Econt(idmat)
         CALL mat_free(Ddiff2(idmat))
         Econt(idmat) = Econt(idmat) + DeltaEcontdiff + twoElectronEnergyCont0
         write(lupri,'(A30,I3,A2,ES26.16)')'Econt(idmat=',idmat,')=',Econt(idmat)
      enddo
   ELSE
      twoElectronEnergyCont0 = 0.0E0_realk
      call II_get_Econt1(LUPRI,LUERR,SETTING,D,Econt,ndmat,modThresh)   
      do idmat=1,ndmat
         write(lupri,'(A30,I3,A2,ES26.16)')'Econt(idmat=',idmat,')=',Econt(idmat)
      enddo
   ENDIF
ENDIF

END SUBROUTINE II_get_Econt

!> \brief Calculates the Fock matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param F the exchange matrix
!> \param ndmat the number of density matrices
SUBROUTINE II_get_Econt1(LUPRI,LUERR,SETTING,DAO,Econt,ndmat,modThresh)
IMPLICIT NONE
integer               :: ndmat,LUPRI,LUERR
TYPE(MATRIX)          :: DAO(ndmat)
real(realk)           :: Econt(ndmat)
TYPE(LSSETTING)       :: SETTING
real(realk)           :: modThresh
!
Real(realk)  :: Econt2(ndmat)
integer :: idmat
real(realk)         :: OLDTHRESH,Kfac

IF(setting%IntegralTransformGC)THEN
   call lsquit('IntegralTransformGC in II_get_Econt1 use .NOGCBASIS',-1)
ENDIF
Kfac = SETTING%SCHEME%exchangeFactor
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%J_THR*modThresh

WRITE(lupri,'(A,ES20.10)')'Screening Threshold for Coulomb Energy Cont=',SETTING%SCHEME%intTHRESHOLD
WRITE(lupri,'(A,ES20.10)')'Screening Threshold for Coulomb Matrix Cont=',SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%J_THR
call II_get_CoulombEcont(LUPRI,LUERR,SETTING,DAO,Econt,ndmat)
IF (SETTING%SCHEME%exchangeFactor.GT. 0.0E0_realk) THEN
   SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%K_THR*(1.0E0_realk/Kfac)*modThresh
   WRITE(lupri,'(A,ES20.10)')'Screening Threshold for Exchange Energy Cont=',SETTING%SCHEME%intTHRESHOLD
   WRITE(lupri,'(A,ES20.10)')'Screening Threshold for Exchange Matrix Cont='&
   &,SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%K_THR*(1.0E0_realk/Kfac)
   call II_get_exchangeEcont(LUPRI,LUERR,SETTING,DAO,Econt2,ndmat)
   do idmat=1,ndmat
      WRITE(lupri,'(A30,I3,A3,ES20.10)')'Coulomb Energy cont, Idmat=',idmat,' is',Econt(idmat)
      WRITE(lupri,'(A30,I3,A3,ES20.10)')'Exchange Energy cont, Idmat=',idmat,' is',Kfac*Econt2(idmat)
      Econt(idmat) = Econt(idmat) + Kfac*Econt2(idmat)
      WRITE(lupri,'(A30,I3,A3,ES20.10)')'Energy cont, Idmat=',idmat,' is',Econt(idmat)
   enddo
ENDIF
END SUBROUTINE II_get_Econt1

!> \brief 
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param ndmat the number of density matrices
!> \param F the exchange matrix
SUBROUTINE II_get_exchangeEcont(LUPRI,LUERR,SETTING,D,Econt,ndmat)
IMPLICIT NONE
Integer               :: ndmat
TYPE(MATRIX),target   :: D(ndmat)
real(realk)           :: Econt(ndmat)
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
!
Integer             :: idmat,incdmat,nrow,ncol
Real(realk),pointer :: DfullLHS(:,:,:)
Real(realk),pointer :: DfullRHS(:,:,:)
Real(realk)         :: TS,TE,fac
integer    :: Oper,dascreen_thrlog
integer :: nCalcInt,nCalcIntZero,nCalcIntZeroContrib
Logical :: DaLink

IF(setting%IntegralTransformGC)THEN
   call lsquit('IntegralTransformGC in II_get_exchangeEcont use .NOGCBASIS',-1)
ENDIF
DaLink = setting%scheme%daLinK
Dascreen_thrlog = setting%scheme%Dascreen_thrlog
CALL LSTIMER('START ',TS,TE,LUPRI)
!attach matrices
CALL ls_attachDmatToSetting(D,ndmat,setting,'RHS',1,3,.FALSE.,lupri)
CALL ls_attachDmatToSetting(D,ndmat,setting,'LHS',2,4,.FALSE.,lupri)
CALL ls_LHSSameAsRHSDmatToSetting(setting)

setting%scheme%daLinK=setting%scheme%LSdaLinK
IF(.NOT.setting%scheme%LSDASCREEN)setting%scheme%daLinK=.FALSE.
setting%scheme%dascreen_thrlog = setting%scheme%lsdascreen_thrlog
IF (SETTING%SCHEME%CAM) THEN
  Oper = CAMOperator       !Coulomb attenuated method
ELSEIF (SETTING%SCHEME%SR_EXCHANGE) THEN
  Oper = ErfcOperator      !Short-Range Coulomb screened exchange
ELSE
  Oper = CoulombOperator   !Regular Coulomb metric 
ENDIF
!Calculates the HF-exchange contribution
call initIntegralOutputDims(setting%Output,1,1,1,1,ndmat)
call ls_get_exchange_mat(AORdefault,AORdefault,AORdefault,AORdefault,&
     &                   Oper,EcontribSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,Econt,setting%IntegralTransformGC)

CALL ls_freeDmatFromSetting(setting)
CALL ls_LHSSameAsRHSDmatToSetting_deactivate(setting)
CALL LSTIMER('st-KEcont',TS,TE,LUPRI)
setting%scheme%daLinK = DaLink
setting%scheme%Dascreen_thrlog = Dascreen_thrlog
END SUBROUTINE II_get_exchangeEcont

!> \brief Calculates the coulomb matrix using the jengine method (default)
!> \author S. Reine
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the array of density matrix
!> \param E the coulomb energies 
!> \param ndmat the number of density matrix
SUBROUTINE II_get_CoulombEcont(LUPRI,LUERR,SETTING,D,Econt,ndmat)
IMPLICIT NONE
Integer             :: ndmat
TYPE(MATRIX),target   :: D(ndmat)
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,nbasis
real(realk)           :: Econt(ndmat)
!
Real(realk)         :: TS,TE
logical :: FMM,DaJengine,DaCoulomb
integer :: naux,Dascreen_thrlog
integer      :: oper,natoms
Real(realk),pointer :: g1alphafull(:,:,:,:,:)
Real(realk),pointer :: calphafull(:,:,:)

IF(setting%IntegralTransformGC)THEN
   call lsquit('IntegralTransformGC in II_get_CoulombEcont use .NOGCBASIS',-1)
ENDIF
CALL LSTIMER('START ',TS,TE,LUPRI)
FMM = SETTING%SCHEME%FMM
SETTING%SCHEME%FMM = .FALSE.
call getMolecularDimensions(SETTING%MOLECULE(1)%p,nAtoms,nBasis,nAux)

DaJengine = setting%scheme%DaJengine
DaCoulomb = setting%scheme%DaCoulomb
Dascreen_thrlog = setting%scheme%dascreen_thrlog
setting%scheme%dascreen_thrlog = setting%scheme%lsdascreen_thrlog

IF (SETTING%SCHEME%DENSFIT) THEN
   Oper = CoulombOperator
   IF(SETTING%SCHEME%OVERLAP_DF_J) Oper = OverlapOperator
   CALL ls_attachDmatToSetting(D,ndmat,setting,'RHS',3,4,.TRUE.,lupri)
   call initIntegralOutputDims(setting%Output,naux,1,1,1,1)
   call ls_jengine(AODFdefault,AOempty,AORdefault,AORdefault,Oper,RegularSpec,ContractedInttype,&
        &          SETTING,LUPRI,LUERR)
   call mem_alloc(g1alphafull,naux,1,1,1,ndmat)
   CALL retrieve_Output(lupri,setting,g1alphafull,.FALSE.)
   IF(SETTING%SCHEME%FMM.AND.(Oper.EQ.CoulombOperator)) THEN
      CALL LSQUIT('Ecotn FMM not tested',-1)
!      call ls_jengineClassicalfull(g1alphafull,AOempty,AORdefault,AORdefault,&
!           &CoulombOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
   ENDIF
   CALL ls_freeDmatFromSetting(setting)
   call LSTIMER('GALPHA',TS,TE,LUPRI)

   call mem_alloc(calphafull,naux,1,ndmat)
   !|rho_fit) = <alpha|w|beta>^-1 <beta|w|rho>
   call linsolv_df(calphafull,g1alphafull,AODFdefault,Oper,naux,ndmat,SETTING,LUPRI,LUERR)
   call mem_dealloc(g1alphafull)
   call LSTIMER('CALPHA',TS,TE,LUPRI)
   !Econt = C_alpha (alpha|beta) C_beta
   CALL ls_attachDmatToSetting(calphafull,naux,1,ndmat,setting,'LHS',1,2,lupri)
   CALL ls_attachDmatToSetting(calphafull,naux,1,ndmat,setting,'RHS',3,4,lupri)
   DaJengine = setting%scheme%DaJengine
   setting%scheme%DaJengine=setting%scheme%LSDaJengine
   IF(.NOT.setting%scheme%LSDASCREEN)setting%scheme%DaJengine=.FALSE.
   call initIntegralOutputDims(setting%Output,1,1,1,1,ndmat)
   call ls_jengine(AODFdefault,AOempty,AODFdefault,AOempty,CoulombOperator,EcontribSpec,&
        & ContractedInttype,SETTING,LUPRI,LUERR)
   CALL retrieve_Output(lupri,setting,Econt,.FALSE.)   
   IF(SETTING%SCHEME%FMM.AND.(Oper.EQ.CoulombOperator)) THEN
      call lsquit('not tested',-1)
!      call ls_jengineClassicalfull
   ENDIF
   call mem_dealloc(calphafull)
   CALL ls_freeDmatFromSetting(setting)
   CALL LSTIMER('df-J Econt',TS,TE,LUPRI)
ELSE
   CALL ls_attachDmatToSetting(D,ndmat,setting,'LHS',1,2,.TRUE.,lupri)
   CALL ls_attachDmatToSetting(D,ndmat,setting,'RHS',3,4,.TRUE.,lupri)
   CALL ls_LHSSameAsRHSDmatToSetting(setting)
   call initIntegralOutputDims(setting%Output,1,1,1,1,ndmat)
   IF(SETTING%SCHEME%LSDAJENGINE)THEN
      setting%scheme%DaJengine=setting%scheme%LSDaJengine
      IF(.NOT.setting%scheme%LSDASCREEN)setting%scheme%DaJengine=.FALSE.
      call ls_jengine(AORdefault,AORdefault,AORdefault,AORdefault,&
           & CoulombOperator,EcontribSpec,ContractedInttype,SETTING,LUPRI,LUERR)
   ELSEIF(SETTING%SCHEME%LSDACoulomb)THEN
      setting%scheme%DaCoulomb=setting%scheme%LSDaCoulomb
      IF(.NOT.setting%scheme%LSDASCREEN)setting%scheme%DaCoulomb=.FALSE.
      call ls_get_coulomb_mat(AORdefault,AORdefault,AORdefault,AORdefault,&
           & CoulombOperator,EcontribSpec,ContractedInttype,SETTING,LUPRI,LUERR)
   ELSE
      call lsquit('EcontDaJENGINE or EcontDaCoulomb must be true',-1)
   ENDIF
   CALL retrieve_Output(lupri,setting,Econt,setting%IntegralTransformGC)
   IF(SETTING%SCHEME%FMM) THEN
      call lsquit('not tested',-1)
   ENDIF
   CALL ls_freeDmatFromSetting(setting)
   CALL ls_LHSSameAsRHSDmatToSetting_deactivate(setting)
   IF(SETTING%SCHEME%JENGINE)THEN
      CALL LSTIMER('JengineEcont',TS,TE,LUPRI)
   ELSE
      CALL LSTIMER('CoulombEcont',TS,TE,LUPRI)
   ENDIF
ENDIF
SETTING%SCHEME%FMM = FMM
setting%scheme%DaJengine = DaJengine 
setting%scheme%DaCoulomb = DaCoulomb
setting%scheme%dascreen_thrlog = Dascreen_thrlog

END SUBROUTINE II_get_CoulombEcont

!> \brief Calculates the (ab|cd) with fixed a and b batchindexes so that the output would be a 4dim tensor with dim (dimAbatch,dimBbatch,fulldimC,fulldimD)
!> \author T. Kjaergaard
!> \date 2010-03-17
!> \param integrals the (ab|cd) tensor
!> \param batchA the requested A batchindex
!> \param batchB the requested B batchindex
!> \param dimA dimension of the requested  A batchindex
!> \param dimB dimension of the requested  B batchindex
!> \param nbast The number of oribtals (or basis functions) for c and d
!> \param lupri The default print-unit. If -1 on input it opens default LSDALTON.OUT
!> \param luerr The dafault error-unit. If -1 on input it opens default LSDALTON.ERR
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE II_get_ABres_4CenterEri(setting,integrals,batchA,batchB,dimA,dimB,nbast,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
Integer,intent(in)      :: nbast,dimA,dimB,lupri,luerr,batchA,batchB
Real(realk),intent(out) :: integrals(dimA,dimB,nbast,nbast,1)
!
TYPE(MATRIX),target :: tmp_gab
type(lsitem)        :: ls
Integer             :: mtype_save, nbas,I,J
!type(configItem)    :: config
TYPE(integralconfig)   :: integral
logical,pointer       :: OLDsameMOLE(:,:)
character :: intspec(5)
intspec(1) = 'R'
intspec(2) = 'R'
intspec(3) = 'R'
intspec(4) = 'R'
intspec(5) = 'C'
setting%batchindex(1)=batchA
setting%batchindex(2)=batchB
setting%batchdim(1)=dimA
setting%batchdim(2)=dimB
call mem_alloc(OLDsameMOLE,4,4)
OLDsameMOLE = setting%sameMol
setting%sameMol(1,2)=batchA .EQ.batchB
setting%sameMol(2,1)=batchA .EQ.batchB
DO I=1,2
   DO J=3,4
      setting%sameMol(I,J)=.FALSE.
      setting%sameMol(J,I)=.FALSE.
   ENDDO
ENDDO
CALL II_get_4center_eri(LUPRI,LUERR,ls%SETTING,integrals,dimA,dimB,&
     & nbast,nbast,intspec)

setting%batchindex(1)=0
setting%batchindex(2)=0
setting%batchdim(1)=0
setting%batchdim(2)=0
setting%sameMol = OLDsameMOLE
call mem_dealloc(OLDsameMOLE) 

END SUBROUTINE II_get_ABres_4CenterEri

!> \brief Calculates the (ab|cd) with fixed a and b batchindexes so that the output would be a 4dim tensor with dim (dimAbatch,dimBbatch,fulldimC,fulldimD)
!> \author T. Kjaergaard
!> \date 2010-03-17
!> \param integrals the (ab|cd) tensor
!> \param batchA the requested A batchindex
!> \param batchB the requested B batchindex
!> \param dimA dimension of the requested  A batchindex
!> \param dimB dimension of the requested  B batchindex
!> \param nbast The number of oribtals (or basis functions) for c and d
!> \param lupri The default print-unit. If -1 on input it opens default LSDALTON.OUT
!> \param luerr The dafault error-unit. If -1 on input it opens default LSDALTON.ERR
!> Note that if both LSDALTON.OUT and LSDALTON.ERR are already open, and the -1 
!> unit-numbers are provided, the code will crash (when attemting to reopen
!> a file that is already open).
SUBROUTINE II_get_DistancePlot_4CenterEri(setting,luINT,lupri,luerr)
IMPLICIT NONE
TYPE(LSSETTING)         :: SETTING
Integer,intent(in)      :: luint,lupri,luerr
!   
type(moleculeinfo),target :: atomicmoleculeP
type(moleculeinfo),target :: atomicmoleculeQ
integer :: I,J,iAO,nAtoms,nbastP,nbastQ,ATOMP,ATOMQ
integer :: IA,IB,IC,ID
real(realk) :: Xp,Yp,Zp,Xq,Yq,Zq,DISTX,DISTY,DISTZ,Rpq
logical,pointer       :: OLDsameMOLE(:,:)
real(realk),pointer   :: integrals(:,:,:,:)
type(MOLECULE_PT)   :: temp,Point
character :: intspec(5)
intspec(1) = 'R'
intspec(2) = 'R'
intspec(3) = 'R'
intspec(4) = 'R'
intspec(5) = 'C'

temp%p  => setting%MOLECULE(1)%p
call mem_alloc(OLDsameMOLE,4,4)
OLDsameMOLE = setting%sameMol

setting%sameMol(1,2)=.TRUE.
setting%sameMol(2,1)=.TRUE.
setting%sameMol(3,4)=.TRUE.
setting%sameMol(4,3)=.TRUE.
DO I=1,2
   DO J=3,4
      setting%sameMol(I,J)=.FALSE.
      setting%sameMol(J,I)=.FALSE.
   ENDDO
ENDDO
nAtoms = setting%molecule(1)%p%nAtoms
WRITE(luint,'(A16,3X,A16)')'Distance','IntegralValue'
do ATOMP = 1,nAtoms 
   CALL build_atomicmolecule(temp%p,atomicmoleculeP,ATOMP,lupri)
   do iAO=1,2
      Setting%molecule(iAO)%p => atomicmoleculeP
      Setting%fragment(iAO)%p => atomicmoleculeP
   enddo
   nbastP = getNbasis(AORdefault,Contractedinttype,setting%MOLECULE(1)%p,lupri)

   Xp = - Setting%molecule(1)%p%ATOM(1)%CENTER(1) 
   Yp = - Setting%molecule(1)%p%ATOM(1)%CENTER(1) 
   Zp = - Setting%molecule(1)%p%ATOM(1)%CENTER(1) 
   do ATOMQ = ATOMP,nAtoms 
      CALL build_atomicmolecule(temp%p,atomicmoleculeQ,ATOMQ,lupri)
      do iAO=3,4
         Setting%molecule(iAO)%p => atomicmoleculeQ
         Setting%fragment(iAO)%p => atomicmoleculeQ
      enddo
      nbastQ = getNbasis(AORdefault,Contractedinttype,setting%MOLECULE(3)%p,lupri)
      call mem_alloc(integrals,nbastP,nbastP,nbastQ,nbastQ)
      CALL II_get_4center_eri(LUPRI,LUERR,SETTING,integrals,nbastP,nbastP,nbastQ,nbastQ,intspec)

      Xq = Setting%molecule(3)%p%ATOM(1)%CENTER(1) 
      Yq = Setting%molecule(3)%p%ATOM(1)%CENTER(1) 
      Zq = Setting%molecule(3)%p%ATOM(1)%CENTER(1) 
      DISTX = Xp + Xq
      DISTY = Yp + Yq
      DISTZ = Zp + Zq
      Rpq = SQRT(DISTX*DISTX + DISTY*DISTY + DISTZ*DISTZ)
      DO ID=1,nbastQ
         DO IC=ID,nbastQ
            DO IB=1,nbastP
               DO IA=IB,nbastP
                  IF(ABS(integrals(IA,IB,IC,ID)).GT.1.0E-15_realk)THEN
                     WRITE(luint,'(ES16.8,3X,ES16.8)')Rpq,integrals(IA,IB,IC,ID)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      call mem_dealloc(integrals)
      call free_Moleculeinfo(atomicmoleculeQ)
   enddo
   call free_Moleculeinfo(atomicmoleculeP)
enddo
setting%sameMol = OLDsameMOLE
call mem_dealloc(OLDsameMOLE) 
do iAO=1,4
   setting%MOLECULE(iAO)%p => temp%p
   setting%Fragment(iAO)%p => temp%p
enddo

END SUBROUTINE II_get_DistancePlot_4CenterEri

!==========================================================================
!
! Full interfaces using full fortran arrays instead of type matrix     
!
!==========================================================================

!> \brief Calculates the Fock matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param Dsym is the density symmetric
!> \param F the exchange matrix
!> \param ndmat the number of density matrices
SUBROUTINE II_get_Fock_mat_full(LUPRI,LUERR,SETTING,nbast,D,Dsym,F)
IMPLICIT NONE
INTEGER,intent(in)             :: LUPRI,LUERR,nbast
real(realk),target,intent(in)  :: D(nbast,nbast,1)
real(realk),intent(inout)      :: F(nbast,nbast,1)
TYPE(LSSETTING),intent(inout)  :: SETTING
LOGICAL,intent(in)    :: Dsym
!
real(realk),pointer   :: DAO(:,:,:),KAO(:,:,:)
logical               :: default,IntegralTransformGC
real(realk)           :: TS,TE,fac,maxelm
integer               :: i
real(realk), external :: ddot
default = SETTING%SCHEME%DENSFIT.OR.SETTING%SCHEME%PARI_J.OR.SETTING%SCHEME%JENGINE &
     &      .OR. SETTING%SCHEME%LinK .OR. (matrix_type .EQ. mtype_unres_dense).OR. &
     &      SETTING%SCHEME%CAM

fac = 2E0_realk
IF(matrix_type .EQ. mtype_unres_dense)fac = 1E0_realk
IF(setting%IntegralTransformGC)THEN
   setting%IntegralTransformGC = .FALSE. 
   !change D to AO basis (currently in GCAO basis)
   call mem_alloc(DAO,nbast,nbast,1)
   call GCAO2AO_transform_fullD(D,DAO,nbast,1,setting,lupri)
   IntegralTransformGC = .TRUE.
ELSE
   IntegralTransformGC = .FALSE.
   DAO => D
ENDIF
IF (default) THEN
   call II_get_coulomb_mat_full(LUPRI,LUERR,SETTING,nbast,DAO,F)
   WRITE(lupri,*)'The Coulomb energy contribution ',&
        & fac*0.5E0_realk*ddot(nbast*nbast,DAO,1,F,1)

   call mem_alloc(KAO,nbast,nbast,1)
   call ls_dzero(KAO,nbast*nbast)
   IF (setting%scheme%daLinK) THEN
      CALL ls_attachDmatToSetting(DAO,nbast,nbast,1,setting,'LHS',2,4,lupri)
   ENDIF
   call II_get_exchange_mat_full(LUPRI,LUERR,SETTING,nbast,DAO,Dsym,KAO)

   WRITE(lupri,*)'The Exchange energy contribution ',&
        &-ddot(nbast*nbast,DAO,1,KAO,1)
   call daxpy(nbast*nbast,1E0_realk,KAO,1,F,1)

   call mem_dealloc(KAO)
ELSE 
   call lsquit('II_get_coulomb_and_exchange_mat_full not implemented',-1)
ENDIF
IF(IntegralTransformGC)THEN
   !transform back to GCAO basis 
   call AO2GCAO_transform_fullF(F,nbast,setting,lupri)
   CALL mem_dealloc(DAO)
   setting%IntegralTransformGC = .TRUE. !back to original value
ENDIF
WRITE(lupri,*)'The Fock energy contribution ',-ddot(nbast*nbast,D,1,F,1)
END SUBROUTINE II_get_Fock_mat_full

SUBROUTINE II_get_coulomb_mat_full(LUPRI,LUERR,SETTING,nbast,D,F)
IMPLICIT NONE
INTEGER               :: LUPRI,LUERR,nbast
real(realk)           :: D(nbast,nbast,1),F(nbast,nbast,1)
TYPE(LSSETTING)       :: SETTING
CALL II_get_coulomb_mat_mixed_full(LUPRI,LUERR,SETTING,nbast,nbast,nbast,nbast,D,F,&
     & AORdefault,AORdefault,AORdefault,AORdefault,coulombOperator)
END SUBROUTINE II_get_coulomb_mat_full

!> \brief Calculates the coulomb matrix
!> \author S. Reine
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param F the coulomb matrix
!> \param ndmat the number of density matrix
SUBROUTINE II_get_coulomb_mat_mixed_full(LUPRI,LUERR,SETTING,n1,n2,n3,n4,D,F,&
     & AO1,AO2,AO3,AO4,Oper)
IMPLICIT NONE
Integer               :: n1,n2,n3,n4
Real(realk),target    :: D(n3,n4,1)
Real(realk)           :: F(n1,n2,1)
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,AO1,AO2,AO3,AO4,Oper
!
Real(realk)           :: TS,TE
integer               :: I
real(realk),pointer   :: DAO(:,:,:)
real(realk),pointer   :: Dmat(:,:,:)
logical :: IntegralTransformGC
!
call time_II_operations1()
CALL LSTIMER('START ',TS,TE,LUPRI)

IF(setting%IntegralTransformGC)THEN
   setting%IntegralTransformGC = .FALSE. 
   !change D to AO basis (currently in GCAO basis)
   IF (n3.NE.n4) CALL LSQUIT('Error in II_get_coulomb_mat_mixed_full. n3.NE.n4 and TransformGC',lupri)
   call mem_alloc(DAO,n3,n4,1)
   call GCAO2AO_transform_fullD(D,DAO,n3,1,setting,lupri)
   Dmat => DAO
   IntegralTransformGC = .TRUE.
ELSE
   Dmat => D
   IntegralTransformGC=.FALSE.
ENDIF

IF (SETTING%SCHEME%DENSFIT.OR.SETTING%SCHEME%PARI_J) THEN
  !Consistency testing
  IF ((AO1.NE.AORdefault).OR.(AO2.NE.AORdefault).OR.(AO3.NE.AORdefault).OR.(AO4.NE.AORdefault)&
     &.OR.(Oper.NE.CoulombOperator)) call lsquit('Error in II_get_coulomb_mat_mixed',-1)
  call lsquit('II_get_df_coulomb_mat_full not imple',-1)
!  CALL II_get_df_coulomb_mat(LUPRI,LUERR,SETTING,Dmat,F,ndmat)
ELSE
   IF(SETTING%SCHEME%JENGINE)THEN
      call II_get_jengine_mat_full(LUPRI,LUERR,SETTING,n1,n2,n3,n4,D,F,AO1,AO2,AO3,AO4,Oper)
   ELSE
      call lsquit('II_get_4center_coulomb_mat_full not impl. ',-1)
   ENDIF
ENDIF

IF(IntegralTransformGC)THEN
   !transform back to GCAO basis 
   IF (n1.NE.n2) CALL LSQUIT('Error in II_get_coulomb_mat_mixed_full. n1.NE.n2 and TransformGC',lupri)
   call AO2GCAO_transform_fullF(F,n1,setting,lupri)
   CALL mem_dealloc(DAO)
   setting%IntegralTransformGC = .TRUE. 
ENDIF

IF (SETTING%SCHEME%DENSFIT.OR.SETTING%SCHEME%PARI_J) THEN
   IF (SETTING%SCHEME%PARI_J) THEN
      IF (SETTING%SCHEME%SIMPLE_PARI) THEN
         CALL LSTIMER('dfJF-PariSimple',TS,TE,LUPRI)
      ELSE 
         CALL LSTIMER('dfJF-Pari',TS,TE,LUPRI)
      ENDIF
   ELSE IF (SETTING%SCHEME%OVERLAP_DF_J) THEN
      CALL LSTIMER('dfJF-Overlap',TS,TE,LUPRI)
   ELSE
      CALL LSTIMER('dfJF-Jengine',TS,TE,LUPRI)
   ENDIF
ELSE
   IF(SETTING%SCHEME%JENGINE)THEN
      CALL LSTIMER('regF-Jengine',TS,TE,LUPRI)
   ELSE
      CALL LSTIMER('regF-Jbuild',TS,TE,LUPRI)
   ENDIF
ENDIF
call time_II_operations2(JOB_II_get_coulomb_mat)

END SUBROUTINE II_get_coulomb_mat_mixed_full

!> \brief Calculates the overlap integral <rho1 rho2>
!> \author S. Reine
!> \date 2014-11-19
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D1 density matrix number 1
!> \param D2 density matrix number 2
!> \param AO1 spcifies orbitals on density number 1
!> \param AO2 spcifies orbitals on density number 2
!> \param GCAO1 spcifies if AO1 is a GC basis
!> \param GCAO2 spcifies if AO2 is a GC basis
!> \param n1 dimension of AO1
!> \param n2 dimension of AO2
!> \param ndmat the number of density matrix
Function II_get_rho2_mat(LUPRI,LUERR,SETTING,D1mat,D2mat,AO1,AO2,GCAO1,GCAO2,n1,n2)
IMPLICIT NONE
INTEGER               :: LUPRI,LUERR,n1,n2,AO1,AO2
Logical               :: GCAO1,GCAO2
TYPE(matrix),target   :: D1mat
TYPE(matrix),target   :: D2mat
real(realk)           :: II_get_rho2_mat
TYPE(LSSETTING)       :: SETTING
!
logical               :: IntegralTransformGC
TYPE(matrix),pointer  :: D1ao,D2ao
TYPE(matrix)          :: F1ao(1)

IntegralTransformGC = setting%IntegralTransformGC
D1ao => D1mat
D2ao => D2mat
IF (IntegralTransformGC) THEN
   setting%IntegralTransformGC = .FALSE.
   IF (GCAO1) THEN
     allocate(D1ao)
     call mat_init(D1ao,n1,n1)
     call GCAO2AO_transform_matrixD2(D1mat,D1ao,setting,lupri)
   ENDIF
   IF (GCAO2) THEN
     allocate(D2ao)
     call mat_init(D2ao,n2,n2)
     call GCAO2AO_transform_matrixD2(D2mat,D2ao,setting,lupri)
   ENDIF
ENDIF
call mat_init(F1ao(1),n1,n1)

CALL II_get_coulomb_mat_mixed(LUPRI,LUERR,SETTING,(/D2ao/),F1ao,1,AO1,AO1,AO2,AO2,overlapOperator)

II_get_rho2_mat = mat_trAB(F1ao(1),D1ao)

call mat_free(F1ao(1))
IF (IntegralTransformGC) THEN
  setting%IntegralTransformGC = IntegralTransformGC 
  IF (GCAO1) THEN
    call mat_free(D1ao)
    deallocate(D1ao)
  ENDIF
  IF (GCAO2) THEN
    call mat_free(D2ao)
    deallocate(D2ao)
  ENDIF
ENDIF

END Function II_get_rho2_mat

!> \brief Calculates the coulomb matrix using the jengine method (default)
!> \author S. Reine
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param F the coulomb matrix
!> \param ndmat the number of density matrix
SUBROUTINE II_get_jengine_mat_full(LUPRI,LUERR,SETTING,n1,n2,n3,n4,D,F,&
     & AO1,AO2,AO3,AO4,Oper)
IMPLICIT NONE
Integer             :: n1,n2,n3,n4,LUPRI,LUERR,AO1,AO2,AO3,AO4,Oper
Real(realk)         :: D(n3,n4,1)
Real(realk)         :: F(n1,n2,1)
TYPE(LSSETTING)     :: SETTING
!
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%J_THR
CALL ls_dzero(F,n1*n2)
IF(SETTING%SCHEME%FMM)THEN
   call lsquit('Error in II_get_jengine_mat_full: FMM',-1)
ENDIF
CALL ls_attachDmatToSetting(D,n3,n4,1,setting,'RHS',3,4,lupri)
call initIntegralOutputDims(setting%Output,n1,n2,1,1,1)
call ls_jengine(AO1,AO2,AO3,AO4,Oper,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,F,setting%IntegralTransformGC)
CALL ls_freeDmatFromSetting(setting)
END SUBROUTINE II_get_jengine_mat_full

!> \brief Calculates the exchange matrix
!> \author S. Reine
!> \date   May 2012
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param ndmat the number of density matrices
!> \param Dsym is the density symmetric
!> \param F the exchange matrix
SUBROUTINE II_get_exchange_mat_full(LUPRI,LUERR,SETTING,nbast,D,Dsym,F)
IMPLICIT NONE
INTEGER,intent(in)    :: LUPRI,LUERR,nbast
real(realk)           :: D(nbast,nbast,1),F(nbast,nbast,1)
TYPE(LSSETTING)       :: SETTING
LOGICAL,intent(in)    :: Dsym

CALL II_get_exchange_mat_mixed_full(LUPRI,LUERR,SETTING,nbast,D,Dsym,F,&
     &    AORdefault,AORdefault,AORdefault,AORdefault,coulombOperator)

END SUBROUTINE II_get_exchange_mat_full


!> \brief Calculates the exchange matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param ndmat the number of density matrices
!> \param Dsym is the density symmetric
!> \param F the exchange matrix
SUBROUTINE ii_get_exchange_mat_mixed_full(LUPRI,LUERR,SETTING,nbast,D,&
     & Dsym,F,AO1,AO3,AO2,AO4,Oper)
IMPLICIT NONE
INTEGER,intent(in)    :: LUPRI,LUERR,nbast,AO1,AO2,AO3,AO4,Oper
real(realk)           :: D(nbast,nbast,1),F(nbast,nbast,1)
TYPE(LSSETTING)       :: SETTING
LOGICAL,intent(in)    :: Dsym
!
real(realk),pointer :: D2(:,:,:),TMP(:,:,:)
integer :: isym(1)
IF (SETTING%SCHEME%exchangeFactor.EQ. 0.0E0_realk) RETURN
call time_II_operations1()

! Check symetry. Split non-symmetric matrices to symmetric and anti symmetric parts. 
! Make symmetry check
IF(.NOT.DSYM)THEN
   ISYM(1) = matfull_get_isym(D,nbast,nbast)
   IF(ISYM(1).EQ.1.OR.ISYM(1).EQ.2)THEN !sym or antisym
      call II_get_exchange_mat1_full(LUPRI,LUERR,SETTING,nbast,D,F,1,AO1,AO3,AO2,AO4,Oper)
   ELSEIF(ISYM(1).EQ.3)THEN !nonsym
      call mem_alloc(D2,nbast,nbast,2)
      call mem_alloc(TMP,nbast,nbast,2)
      call ls_transpose(D(:,:,1),TMP(:,:,1),nbast)
!      call ls_transpose(D(:,:,1),TMP(:,:,2),nbast)
      call dcopy(nbast*nbast,D,1,D2(:,:,1),1)
      call dcopy(nbast*nbast,D,1,D2(:,:,2),1)
      call dscal(nbast*nbast*2,0.5E0_realk,D2,1)
      call daxpy(nbast*nbast,0.5E0_realk,TMP,1,D2,1)   !Sym
      call daxpy(nbast*nbast,-0.5E0_realk,TMP,1,D2(:,:,2),1)  !Asym
      call ls_dzero(TMP,nbast*nbast*2)
      call II_get_exchange_mat1_full(LUPRI,LUERR,SETTING,nbast,D2,TMP,2,AO1,AO3,AO2,AO4,Oper)
      call daxpy(nbast*nbast,1.0E0_realk,TMP(:,:,1),1,F,1) !Sym 
      call daxpy(nbast*nbast,1.0E0_realk,TMP(:,:,2),1,F,1) !Asym
      call mem_dealloc(TMP)
      call mem_dealloc(D2)
   ELSE !zero 
      call ls_dzero(F,nbast*nbast)
   ENDIF
ELSE
   call II_get_exchange_mat1_full(LUPRI,LUERR,SETTING,nbast,D,F,1,AO1,AO3,AO2,AO4,Oper)
ENDIF
call time_II_operations2(JOB_II_GET_EXCHANGE_MAT)

END SUBROUTINE ii_get_exchange_mat_mixed_full

!> \brief Calculates the exchange matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param ndmat the number of density matrices
!> \param F the exchange matrix
SUBROUTINE II_get_exchange_mat1_full(LUPRI,LUERR,SETTING,nbast,D,F,ndmat,AO1,AO3,AO2,AO4,Oper)
IMPLICIT NONE
INTEGER               :: LUPRI,LUERR,AO1,AO2,AO3,AO4,Oper,nbast,ndmat
real(realk),target    :: D(nbast,nbast,ndmat)
real(realk)           :: F(nbast,nbast,ndmat)
TYPE(LSSETTING)       :: SETTING
!
integer :: idmat
real(realk),pointer   :: DAO(:,:,:)
real(realk),pointer   :: Dmat(:,:,:)
logical :: IntegralTransformGC
!
call time_II_operations1()

IF(setting%IntegralTransformGC)THEN
   setting%IntegralTransformGC = .FALSE. 
   !change D to AO basis (currently in GCAO basis)
   call mem_alloc(DAO,nbast,nbast,ndmat)
   call GCAO2AO_transform_fullD(D,DAO,nbast,ndmat,setting,lupri)
   Dmat => DAO
   IntegralTransformGC = .TRUE.
ELSE
   Dmat => D
   IntegralTransformGC=.FALSE.
ENDIF

IF (SETTING%SCHEME%DF_K) THEN
   CALL LSQUIT('Error in II_get_exchange_mat1. DF_K is only implemented for MAT ',-1)

ELSE IF (SETTING%SCHEME%PARI_K.OR.SETTING%SCHEME%MOPARI_K) THEN
   CALL LSQUIT('Error in II_get_exchange_mat1. PARI_K is only implemented for MAT',-1)
ELSE
   CALL II_get_exchange_mat_regular_full(LUPRI,LUERR,SETTING,nbast,Dmat,F,ndmat,AO1,AO3,AO2,AO4,Oper)
ENDIF

IF(IntegralTransformGC)THEN
   !transform back to GCAO basis 
   do idmat=1,ndmat
      call AO2GCAO_transform_fullF(F(:,:,idmat),nbast,setting,lupri)
   enddo
   CALL mem_dealloc(DAO)
   setting%IntegralTransformGC = .TRUE. 
ENDIF
END SUBROUTINE II_get_exchange_mat1_full

!> \brief Calculates the exchange matrix using explicit 4 center 
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param ndmat the number of density matrices
!> \param F the exchange matrix
SUBROUTINE II_get_exchange_mat_regular_full(LUPRI,LUERR,SETTING,nbast,D,F,ndmat,AO1,AO3,AO2,AO4,OperIn)
IMPLICIT NONE
INTEGER              :: LUPRI,LUERR,AO1,AO2,AO3,AO4,OperIn,nbast,ndmat
real(realk),target   :: D(nbast,nbast,ndmat)
real(realk)          :: F(nbast,nbast,ndmat)
TYPE(LSSETTING)      :: SETTING
!
real(realk),pointer  :: K(:,:,:)
Integer             :: idmat,incdmat,nrow,ncol
Real(realk)         :: TS,TE,fac
integer    :: Oper
integer :: nCalcInt,nCalcIntZero,nCalcIntZeroContrib,i

CALL LSTIMER('START ',TS,TE,LUPRI)
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%K_THR
IF(matrix_type .EQ. mtype_unres_dense) THEN
   call lsquit('error unres II_get_exchange_mat_regular_full',-1)
ENDIF

CALL ls_attachDmatToSetting(D,nbast,nbast,ndmat,setting,'RHS',2,4,lupri)
IF (SETTING%SCHEME%CAM) THEN
  Oper = CAMOperator       !Coulomb attenuated method
ELSEIF (SETTING%SCHEME%SR_EXCHANGE) THEN
  Oper = ErfcOperator      !Short-Range Coulomb screened exchange
ELSE
  Oper = OperIn   !Deafult is the Coulomb metric 
ENDIF
!Calculates the HF-exchange contribution
call initIntegralOutputDims(setting%Output,nbast,nbast,1,1,ndmat)
call ls_get_exchange_mat(AO1,AO3,AO2,AO4,Oper,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
call mem_alloc(K,nbast,nbast,ndmat)
CALL retrieve_Output(lupri,setting,K,setting%IntegralTransformGC)
call daxpy(nbast*nbast*ndmat,-1E0_realk,K,1,F,1)
call mem_dealloc(K)
CALL ls_freeDmatFromSetting(setting)
IF(SETTING%SCHEME%DALINK .AND.SETTING%SCHEME%LINK)THEN
   CALL LSTIMER('DaLINK-KFbuild',TS,TE,LUPRI)
ELSEIF(SETTING%SCHEME%LINK)THEN
   CALL LSTIMER('LINK-KFbuild',TS,TE,LUPRI)
ELSE
   CALL LSTIMER('st-KFbuild',TS,TE,LUPRI)
ENDIF
END SUBROUTINE II_get_exchange_mat_regular_full

SUBROUTINE II_get_admm_exchange_mat(LUPRI,LUERR,SETTING,optlevel,D,F,dXC,ndmat,&
     & EdXC,dsym,Econt,ADMMBASISFILE)
implicit none
Integer,intent(in)            :: lupri,luerr,optlevel,ndmat
real(realk),intent(out)       :: EdXC
real(realk),intent(inout)     :: Econt(5)
TYPE(Matrix),intent(in)       :: D
TYPE(Matrix),intent(inout)    :: F,dXC
type(lssetting),intent(inout) :: setting
logical,intent(in)            :: dsym,ADMMBASISFILE
!
logical             :: GC3
TYPE(Matrix)        :: D2(1),TMP,TMPF,k2(1),x2(1),F3(1),R33,S33,Ftmp(1)
TYPE(Matrix)        :: S22
real(realk)         :: var
logical             :: ADMMexchange,testNelectrons,unres,grid_done
real(realk)         :: ex2(1),ex3(1),Edft_corr,ts,te,tsfull,tefull,hfweight
integer             :: nbast,nbast2,AORold,AO3,nelectrons
character(21)       :: L2file,L3file
real(realk)         :: GGAXfactor,fac,fac_unrest
real(realk)         :: constrain_factor, largeLambda,smallK,largeK,mixK
logical             :: isADMMQ,separateX,DODISP
logical             :: isADMMS, isADMMP,PRINT_EK3,saveDF,printExchangeMetric,unrest
 !
unrest = matrix_type.EQ.mtype_unres_dense
fac_unrest = 2E0_realk
IF (unrest) fac_unrest = 1E0_realk
 !
nelectrons = setting%molecule(1)%p%nelectrons 
isADMMQ = setting%scheme%ADMMQ
isADMMS = setting%scheme%ADMMS
isADMMP = setting%scheme%ADMMP
PRINT_EK3 = setting%scheme%PRINT_EK3
separateX = setting%scheme%admm_separateX
separateX = separateX.OR.(.NOT.setting%do_dft) !Hack for HF - for now SR
separateX = separateX.OR.setting%scheme%cam    !Hack for camb3lyp - for now SR
printExchangeMetric = setting%scheme%ADMMexchangeMetric

IF (setting%scheme%cam) THEN
  GGAXfactor = 1.0E0_realk
ELSE
  GGAXfactor = setting%scheme%exchangeFactor
ENDIF
IF(ndmat.GT.1) THEN
   WRITE(*,*)     "The ADMM approximation isn't implemented for ndmat .GT.0 ."
   WRITE(LUPRI,*) "The ADMM approximation isn't implemented for ndmat .GT.0 ."
   call lsquit('II_get_admm_exchange_mat ndmat.GT.1',-1)
ENDIF
nbast = F%nrow
unres = matrix_type .EQ. mtype_unres_dense

CALL lstimer('START',tsfull,tefull,lupri)
CALL lstimer('START',ts,te,lupri)

nbast2 = getNbasis(AOadmm,Contractedinttype,setting%MOLECULE(1)%p,6)

call mat_init(D2(1),nbast2,nbast2)
call mat_init(k2(1),nbast2,nbast2)
call mat_init(x2(1),nbast2,nbast2)
call mat_init(TMPF,nbast,nbast)         

call mat_zero(dXC)

GC3 = setting%IntegralTransformGC
setting%IntegralTransformGC =.FALSE.

!Select regular basis to either optlevel 2 or 3
IF (optlevel.EQ.3) THEN
  AO3 = AORdefault
ELSEIF (optlevel.EQ.2) THEN
  AO3 = AOVAL
ELSE
  CALL LSQUIT('II_get_admm_exchange_mat:Error in ADMM, unknown optlevel',-1)
ENDIF
       
constrain_factor = 1.0E0_realk
IF (isADMMQ .OR. isADMMS .OR. isADMMP) THEN   
   call get_small_lambda(constrain_factor,D,setting,lupri,luerr,nbast2,nbast,AOadmm,AO3,.FALSE.,GC3)
ENDIF
setting%scheme%ADMM_CONSTRAIN_FACTOR = constrain_factor

!We transform the full Density to a level 2 density D2
call transform_D3_to_D2(D,D2(1),setting,lupri,luerr,nbast2,&
                  & nbast,AOadmm,AO3,setting%scheme%ADMM1,&
                  & .FALSE.,GC3,constrain_factor)

!Store original AO-index
AORold  = AORdefault

!****Calculation of Level 2 exchange matrix from level 2 Density matrix starts here

!ADMM (level 2) AO settings 
call set_default_AOs(AOadmm,AODFdefault)

CALL lstimer('AUX-IN',ts,te,lupri)
call mat_zero(k2(1))
call II_get_exchange_mat(LUPRI,LUERR,SETTING,D2,1,Dsym,k2)


call Transformed_F2_to_F3(TMPF,k2(1),setting,lupri,luerr,nbast2,nbast,&
                        & AOadmm,AO3,.FALSE.,GC3,constrain_factor)
fac = 1E0_realk
IF (isADMMP) fac = constrain_factor**(4.E0_realk)
call mat_daxpy(fac,TMPF,F)
IF(ADMMBASISFILE)Econt(3) = mat_dotproduct(TMPF,D)*fac_unrest/2E0_realk

!Prints the exchange-metric error between the regular and ADMM density
!   int (rho(1,2)-fit(1,2))^2/r_12 dr_1 dr_2
IF (printExchangeMetric) THEN
  smallK = fac*mat_trAB(k2(1),d2(1))
  call mat_init(Ftmp(1),nbast,nbast)
  call mat_zero(Ftmp(1))
  call set_default_AOs(AO3,AODFdefault)
  setting%IntegralTransformGC = GC3
  call II_get_exchange_mat(lupri,luerr,setting,(/D/),1,Dsym,Ftmp)
  setting%IntegralTransformGC =.FALSE.
  call set_default_AOs(AOadmm,AODFdefault)
  largeK = mat_trAB(Ftmp(1),D)
  call mat_zero(Ftmp(1))
  call II_get_exchange_mat_mixed(lupri,luerr,setting,D2,1,Dsym,Ftmp,AO3,AOadmm,AO3,AOadmm,coulombOperator)
  IF (GC3) call AO2GCAO_transform_matrixF(Ftmp(1),setting,lupri)
  mixK   = sqrt(fac)*2.d0*mat_trAB(Ftmp(1),D)
  write(lupri,'(X,A,F18.12)') 'Exchange-metric error:',-(largeK+smallK-mixK)*fac_unrest/2E0_realk
  write(*,'(X,A,F18.12)') 'Exchange-metric error:',-(largeK+smallK-mixK)*fac_unrest/2E0_realk
  call mat_free(Ftmp(1))
ENDIF

!Subtract XC-correction
CALL lstimer('AUX-EX',ts,te,lupri)
!****Calculation of Level 2 XC matrix from level 2 Density matrix starts here
call II_DFTsetFunc(setting%scheme%dft%DFTfuncObject(dftfunc_ADMML2),GGAXfactor,lupri)

dodisp = setting%scheme%dft%dodisp
if(setting%scheme%dft%dodisp) setting%scheme%dft%dodisp = .false.

!choose the ADMM Level 2 grid
setting%scheme%dft%igrid = Grid_ADMML2
  
!Only test electrons if the D2 density matrix is McWeeny purified
testNelectrons = setting%scheme%dft%testNelectrons
setting%scheme%dft%testNelectrons = setting%scheme%ADMM1

!Level 2 XC matrix
call mat_zero(x2(1))
call II_get_xc_Fock_mat(LUPRI,LUERR,SETTING,nbast2,D2,x2,EX2,1)

IF (isADMMS) THEN
   EX2 = constrain_factor**(4./3.)*EX2            ! RE-SCALING EXC2 TO FIT k2
   call mat_scal(constrain_factor**(4./3.),x2(1)) ! RE-SCALING XC2 TO FIT k2  
ENDIF

IF(ADMMBASISFILE)Econt(4) = EX2(1)

!Transform to level 3
call transformed_F2_to_F3(TMPF,x2(1),setting,lupri,luerr,nbast2,nbast,&
                          & AOadmm,AO3,.FALSE.,GC3,constrain_factor)

call mat_daxpy(-fac,TMPF,dXC)
setting%scheme%dft%testNelectrons = testNelectrons

!Re-set to default (level 3) grid
setting%scheme%dft%igrid = Grid_Default

!****Calculation of Level 3 XC matrix from level 2 Density matrix starts here
call set_default_AOs(AORold,AODFdefault)  !Revert back to original settings and free stuff 
setting%IntegralTransformGC = GC3     !Restore GC transformation to level 3

IF (separateX) THEN
  CALL mat_init(F3(1),nbast,nbast)
  CALL mat_zero(F3(1))

  call II_get_xc_Fock_mat(LUPRI,LUERR,SETTING,nbast,(/D/),F3,EX3,1)
  IF(ADMMBASISFILE)Econt(2) = EX3(1)
  
  CALL mat_daxpy(1E0_realk,F3(1),dXC)
  CALL mat_free(F3(1))
  EdXC = (EX3(1)- fac*EX2(1))
ELSE
  EdXC = -fac*EX2(1)
ENDIF

!Restore dft functional to original
IF (setting%do_dft) THEN
  call II_DFTsetFunc(setting%scheme%dft%DFTfuncObject(dftfunc_Default),GGAXfactor,lupri)
  !Augment the functional with the admm gga exchange contribution X
   IF (.NOT.separateX) THEN
     call II_DFTaddFunc(setting%scheme%dft%DFTfuncObject(dftfunc_ADMML2),GGAXfactor,lupri)
   ENDIF
ENDIF

IF (PRINT_EK3) THEN
   write(*,*)     "E(x2)= ", fac*EX2(1)
   write(lupri,*) "E(x2)= ", fac*EX2(1)
   IF (separateX) THEN
      write(*,*)     "E(X3)= ", EX3(1)
      write(lupri,*) "E(X3)= ", EX3(1)
      write(*,*)     "E(X3)-E(x2)= ",EdXC
      write(lupri,*) "E(X3)-E(x2)= ",EdXC
   ENDIF
!***  Simen: 2014-11-19, Added for the basis-set optimization
   !Factor 2 included here because the Coulomb factor used in the FTUVs are set to two for default AOs, and 
   !one for other.
   saveDF = setting%scheme%densfit
   setting%scheme%densfit = .FALSE.
   setting%IntegralTransformGC = .FALSE.
   var = 2E0_realk * II_get_rho2(LUPRI,LUERR,SETTING,D2(1),D2(1),AOadmm,AOadmm,.FALSE.,.FALSE.,nbast2,nbast2)
   setting%IntegralTransformGC = GC3     !Restore GC transformation to level 3
   var = var + II_get_rho2(LUPRI,LUERR,SETTING,D,D,AO3,AO3,GC3,GC3,nbast,nbast)
   var = var - 2E0_realk * II_get_rho2(LUPRI,LUERR,SETTING,D2(1),D,AOadmm,AO3,.FALSE.,GC3,nbast2,nbast)
   write(*,*) "Fitting error = ", fac_unrest * var
   write(lupri,*) "Fitting error = ", fac_unrest * var
   setting%scheme%densfit = saveDF
ENDIF
     


if(dodisp) setting%scheme%dft%dodisp = dodisp


! term of Kadmm coming from dependence of K on lambda: K=[D,lambda(D)]
IF (isADMMQ .OR. isADMMS .OR. isADMMP) THEN

  call mat_init(R33,nbast,nbast) 
  CALL mat_init(S33,nbast,nbast)

  !R = T^T s^-1 T
  CALL get_R33(R33,AOadmm,AO3,GC3,nbast2,nbast,constrain_factor,setting,lupri,luerr)
  CALL get_S33(S33,AO3,GC3,setting,lupri,luerr)
  call mat_daxpy(-1E0_realk,R33,S33)

  call get_large_Lambda(largeLambda,k2(1),x2(1),D2(1),EX2(1),constrain_factor,setting)
  setting%scheme%ADMM_LARGE_LAMBDA = largeLambda

  call mat_daxpy(largeLambda,S33,dXC)

  CALL mat_free(S33)
  CALL mat_free(R33)
ENDIF
call mat_free(TMPF)
call mat_free(k2(1))
call mat_free(x2(1))
call mat_free(D2(1))
CALL lstimer('ADMM',tsfull,tefull,lupri)
END SUBROUTINE II_get_admm_exchange_mat

!> \brief Calculates the ADMM exchange contribution to the molecular gradient
!> \author S. Reine and P. Merlot
!> \date 2013-01-29
!> \param kGrad The ADMM exchange gradient
!> \param DmatLHS The left-hand-side (or first electron) density matrix
!> \param DmatRHS The reft-hand-side (or second electron) density matrix
!> \param ndlhs The number of LHS density matrices
!> \param ndrhs The number of RHS density matrices
!> \param setting Integral evalualtion settings
!> \param lupri Default print unit
!> \param luerr Unit for error printing
SUBROUTINE II_get_ADMM_K_gradient(admm_Kgrad,DmatLHS,DmatRHS,ndlhs,ndrhs,setting,lupri,luerr)
IMPLICIT NONE
!
type(lssetting),intent(inout) :: setting
real(realk),intent(inout)     :: admm_Kgrad(3,setting%molecule(1)%p%nAtoms)
integer,intent(in)            :: lupri,luerr,ndlhs,ndrhs
type(MatrixP),intent(in)      :: DmatLHS(ndlhs),DmatRHS(ndrhs)
!
real(realk),pointer :: grad_k2(:,:),grad_xc2(:,:),grad_XC3(:,:)
real(realk),pointer :: ADMM_proj(:,:),ADMM_charge_term(:,:)
integer             :: nbast,nbast2,AO3,idmat,nAtoms
real(realk)         :: ts,te,hfweight
real(realk)         :: Exc2(1), E_x2
type(Matrix),target :: D2
type(Matrix)        :: k2,xc2,zeromat
type(matrixp)       :: D2p(1)
logical             :: GC3,GC2,testNelectrons,grid_done,DSym,unres
integer             :: AORold
character(len=80)   :: WORD
character(21)       :: L2file,L3file
real(realk)         :: GGAXfactor
real(realk)         :: constrain_factor,nrm
logical             :: PRINT_EK3
logical             :: isADMMQ,isADMMS,isADMMP,DODISP
integer             :: iAtom,iX
!
isADMMQ = setting%scheme%ADMMQ
isADMMS = setting%scheme%ADMMS
isADMMP = setting%scheme%ADMMP
PRINT_EK3 = setting%scheme%PRINT_EK3
call lstimer('START',ts,te,lupri)
IF (setting%scheme%cam) THEN
  GGAXfactor = 1.0E0_realk
ELSE
  GGAXfactor = setting%scheme%exchangeFactor
ENDIF

IF (ndrhs.NE.ndlhs) call lsquit('II_get_ADMM_K_gradient:Different LHS/RHS density matrices not implemented',-1)

nAtoms = setting%molecule(1)%p%nAtoms
call ls_dzero(admm_Kgrad,3*nAtoms)
nbast  = DmatLHS(1)%p%nrow
unres  = matrix_type .EQ. mtype_unres_dense

  
DO idmat=1,ndrhs
   
   nbast2 = getNbasis(AOadmm,Contractedinttype,setting%MOLECULE(1)%p,6)
   
   GC3 = setting%IntegralTransformGC
   setting%IntegralTransformGC = .FALSE.
   
   AO3 = AORdefault

   ! Get the scaling factor derived from constraining the total charge
   ! for the projection, obtained through the minimization
   !
   !      min [ <(rho-rho')^2> + lambda <rho-rho'> ]
   !
   ! and where
   !
   !      constrain_factor = (1-lambda)^-1
   !
   ! Assumes the energy has been calculated first - which stores the constrain_factor
   constrain_factor = setting%scheme%ADMM_CONSTRAIN_FACTOR

   !!We transform the full Density to a level 2 density D2
   call mat_init(D2,nbast2,nbast2)
   call mat_zero(D2)
   call transform_D3_to_D2(DmatLHS(idmat)%p,D2,setting,lupri,luerr,&
                           & nbast2,nbast,AOadmm,AO3,&
                           & setting%scheme%ADMM1,.FALSE.,GC3,&
                           & constrain_factor)
 
   !Store original AO-index
   AORold  = AORdefault
   
   !!****Calculation of Level 2 exchange gradient from level 2 Density matrix starts here
   !ADMM (level 2) AO settings 
   call set_default_AOs(AOadmm,AODFdefault)
   
   D2p(1)%p => D2
   
   ! LEVEL 2 exact exchange matrix
   call mat_init(k2,nbast2,nbast2)
   call mat_zero(k2)
   Dsym = .TRUE. !symmetric Density matrix !!!!!!!!!!!!!!!!!!!    ASSUMPTION WITHOUT LSQUIT HERE      !!!!!!!!!!!!!!!!!!!!!!!!
   call II_get_exchange_mat(lupri,luerr,setting,D2,1,Dsym,k2)
   
   call mem_alloc(grad_k2,3,nAtoms)
   call ls_dzero(grad_k2,3*nAtoms)
   call II_get_regular_K_gradient(grad_k2,D2p,D2p,1,1,setting,lupri,luerr)
   
   call DSCAL(3*nAtoms,4E0_realk,grad_k2,1) !Include factor 4 to use D instead of 2D
   IF (isADMMP) THEN   
      call DSCAL(3*nAtoms,constrain_factor**(4.E0_realk),grad_k2,1)
   ENDIF
   call DAXPY(3*nAtoms,1E0_realk,grad_k2,1,admm_Kgrad,1) 
   CALL LS_PRINT_GRADIENT(lupri,setting%molecule(1)%p,grad_k2,nAtoms,'grad_k2')
   call mem_dealloc(grad_k2)
   
   ! XC-correction
   !****Calculation of Level 2 XC gradient from level 2 Density matrix starts here

   call II_DFTsetFunc(setting%scheme%dft%DFTfuncObject(dftfunc_ADMML2),GGAXfactor,lupri)

   dodisp = setting%scheme%dft%dodisp
   if(setting%scheme%dft%dodisp) setting%scheme%dft%dodisp = .false.
   
   !choose the ADMM Level 2 grid
   setting%scheme%dft%igrid = Grid_ADMML2

   !Only test electrons if the D2 density matrix is McWeeny purified

   testNelectrons = setting%scheme%dft%testNelectrons
   setting%scheme%dft%testNelectrons = .FALSE. !setting%scheme%ADMM1
   
   !Level 2 XC matrix
   call mat_init(xc2,nbast2,nbast2)
   call mat_zero(xc2)
   Exc2(1) = 0.0E0_realk
   E_x2    = 0.0E0_realk
   call II_get_xc_Fock_mat(lupri,luerr,setting,nbast2,D2,xc2,Exc2,1)
   IF (isADMMS) THEN   
      call mat_scal(constrain_factor**(4./3.),xc2)
      E_x2 = Exc2(1)*constrain_factor**(4./3.)
   ELSE
      E_x2 = Exc2(1)
   ENDIF
   
   
   !Level 2 XC gradient
   call mem_alloc(grad_xc2,3,nAtoms)
   call ls_dzero(grad_xc2,3*nAtoms)
   call II_get_xc_geoderiv_molgrad(lupri,luerr,setting,nbast2,D2,grad_xc2,nAtoms)
   
   IF (isADMMS) THEN   
      call DSCAL(3*nAtoms,constrain_factor**(4.E0_realk/3.E0_realk),grad_xc2,1)
   ELSEIF (isADMMP) THEN   
      call DSCAL(3*nAtoms,constrain_factor**(4.E0_realk),grad_xc2,1)
   ENDIF

   call DSCAL(3*nAtoms,-1E0_realk,grad_xc2,1)
   call DAXPY(3*nAtoms,1E0_realk,grad_xc2,1,admm_Kgrad,1)
   CALL LS_PRINT_GRADIENT(lupri,setting%molecule(1)%p,grad_xc2,nAtoms,'grad_xc2')
   call mem_dealloc(grad_xc2)
   
   !Re-set to level 3 grid
   setting%scheme%dft%igrid = Grid_Default
   
   !****Calculation of Level 3 XC gradient
   call set_default_AOs(AORold,AODFdefault)  !Revert back to original settings and free stuff 
   setting%IntegralTransformGC = GC3     !Restore GC transformation to level 3

   call mem_alloc(grad_XC3,3,nAtoms)
   call ls_dzero(grad_XC3,3*nAtoms)
   call II_get_xc_geoderiv_molgrad(lupri,luerr,setting,nbast,DmatLHS(idmat)%p,grad_XC3,nAtoms)
   
   call DAXPY(3*nAtoms,1E0_realk,grad_XC3,1,admm_Kgrad,1)
   
   CALL LS_PRINT_GRADIENT(lupri,setting%molecule(1)%p,grad_XC3,nAtoms,'grad_XC3')
   call mem_dealloc(grad_XC3)

   
   ! set back the default choice for testing the nb. of electrons
   setting%scheme%dft%testNelectrons = testNelectrons
   

   call mem_alloc(ADMM_proj,3,nAtoms)
   call ls_dzero(ADMM_proj,3*nAtoms)
   call mem_alloc(ADMM_charge_term,3,nAtoms)
   call ls_dzero(ADMM_charge_term,3*nAtoms)
   
   IF (isADMMQ.OR.isADMMS.OR.isADMMP) THEN
      call get_ADMM_K_gradient_constrained_charge_term(ADMM_charge_term,k2,xc2,&
                  & E_x2,DmatLHS(idmat)%p,D2,nbast2,nbast,nAtoms,&
                  & AOadmm,AO3,.FALSE.,GC3,setting,lupri,luerr)
      call DSCAL(3*nAtoms,2E0_realk,ADMM_charge_term,1)
      call LS_PRINT_GRADIENT(lupri,setting%molecule(1)%p,ADMM_charge_term,nAtoms,'ADMM-Chrg')  
      call DAXPY(3*nAtoms,1E0_realk,ADMM_charge_term,1,admm_Kgrad,1)
   ENDIF
   
   ! Additional (reorthonormalisation like) projection terms coming from the 
   ! derivative of the small d2 Density matrix
   ! calculating Tr(T^x D3 trans(T) [2 k22(D2) - xc2(D2)]))
   call get_ADMM_K_gradient_projection_term(ADMM_proj,k2,xc2,E_x2,&
               & DmatLHS(idmat)%p,D2,nbast2,nbast,nAtoms,&
               & AOadmm,AO3,.FALSE.,GC3,constrain_factor,setting,lupri,luerr)
   call DSCAL(3*nAtoms,2E0_realk,ADMM_proj,1)
   call LS_PRINT_GRADIENT(lupri,setting%molecule(1)%p,ADMM_proj,nAtoms,'ADMM_proj')  
   call DAXPY(3*nAtoms,1E0_realk,ADMM_proj,1,admm_Kgrad,1)   

   !FREE MEMORY
   call mem_dealloc(ADMM_proj)
   call mem_dealloc(ADMM_charge_term)
   call mat_free(k2)
   call mat_free(xc2)
   call mat_free(D2)
ENDDO !idmat

call DSCAL(3*nAtoms,0.25_realk,admm_Kgrad,1)

!Restore dft functional to original
IF (setting%do_dft) call II_DFTsetFunc(setting%scheme%dft%DFTfuncObject(dftfunc_Default),hfweight,lupri)
if(dodisp) setting%scheme%dft%dodisp = dodisp
!
CONTAINS

   SUBROUTINE get_ADMM_K_gradient_constrained_charge_term(ADMM_charge_term,&
                                 & k2,xc2,E_x2,D3,D2,n2,n3,&
                                 & nAtoms,AO2,AO3,GCAO2,GCAO3,&
                                 & setting,lupri,luerr)
      implicit none
      type(lssetting),intent(inout) :: setting
      real(realk),intent(inout)  :: ADMM_charge_term(3,nAtoms)
      type(matrix),intent(in),target    :: D3 !Regular density matrix
      type(matrix),intent(in),target    :: D2 !ADMM projected density matrix
      type(matrix),intent(in)    :: k2 !ADMM exchange matrix k2(d2)
      type(matrix),intent(in)    :: xc2 !ADMM XC matrix xc2(d2)
      integer,intent(in)         :: nAtoms
      integer,intent(in)         :: n2,n3,lupri,luerr
      integer,intent(in)         :: AO2,AO3
      logical,intent(in)         :: GCAO2,GCAO3
      real(realk),intent(IN)     :: E_x2
      !
      real(realk)                :: LambdaEnergy
      type(matrixp)              :: tmpDFD(1)
      real(realk),pointer        :: reOrtho_D3(:,:),reOrtho_d2(:,:)
      logical                    :: isADMMP
      !
      isADMMP = setting%scheme%ADMMP
      call ls_dzero(ADMM_charge_term,3*nAtoms)

      tmpDFD(1)%p => D3
      call mem_alloc(reOrtho_D3,3,nAtoms)
      call ls_dzero(reOrtho_D3,3*nAtoms)
      call II_get_reorthoNormalization(reOrtho_D3,tmpDFD,1,setting,lupri,luerr)
      call DAXPY(3*nAtoms, 1E0_realk,reOrtho_D3,1,ADMM_charge_term,1)
      
      tmpDFD(1)%p => D2
      call mem_alloc(reOrtho_d2,3,nAtoms)
      call ls_dzero(reOrtho_d2,3*nAtoms)
      call II_get_reorthoNormalization_mixed(reOrtho_d2,tmpDFD,1,AO2,AO2,&
                                       & GCAO2,GCAO2,setting,lupri,luerr)
      IF (isADMMP) THEN   
         call DSCAL(3*nAtoms,constrain_factor**(2.E0_realk),reOrtho_d2,1)
      ENDIF
      call DAXPY(3*nAtoms,-1E0_realk,reOrtho_d2,1,ADMM_charge_term,1)

      ! Assumes the energy has been calculated first - which stores the Lambda
      LambdaEnergy = setting%scheme%ADMM_LARGE_LAMBDA
      call DSCAL(3*nAtoms,LambdaEnergy,ADMM_charge_term,1)   
      
      ! free memory                                 
      call mem_dealloc(reOrtho_D3)
      call mem_dealloc(reOrtho_d2)
   END SUBROUTINE get_ADMM_K_gradient_constrained_charge_term

   
   SUBROUTINE get_ADMM_K_gradient_projection_term(ADMM_proj,k2,xc2,E_x2,D3,D2,&
                                 & n2,n3,nAtoms,AO2,AO3,GCAO2,GCAO3,&
                                 & constrain_factor,setting,lupri,luerr)
      implicit none
      type(lssetting),intent(inout) :: setting
      real(realk),intent(inout)  :: ADMM_proj(3,nAtoms)
      type(matrix),intent(in)    :: D3 !Regular density matrix
      type(matrix),intent(in)    :: D2 !ADMM projected density matrix
      type(matrix),intent(in)    :: k2 !ADMM exchange matrix k2(d2)
      type(matrix),intent(in)    :: xc2 !ADMM XC matrix xc2(d2)
      integer,intent(in)         :: nAtoms
      integer,intent(in)         :: n2,n3,lupri,luerr
      integer,intent(in)         :: AO2,AO3
      logical,intent(in)         :: GCAO2,GCAO3
      real(realk),intent(in)     :: E_x2
      real(realk),intent(in)     :: constrain_factor! Lagrange Mult. for coeff.
      !
      type(matrix),target        :: A22,B32,C22,tmp22
      real(realk)                :: LambdaE ! Lagrange Mult. for energy chg. const.
      type(matrix)               :: S22,S22inv,T23,tmp32
      type(matrixp)              :: tmpDFD(1)
      real(realk),pointer        :: reOrtho1(:,:),reOrtho2(:,:)
      integer                    :: NbEl ! nb. electrons
      Logical                    :: isADMMQ,isADMMS,isADMMP
      !
      isADMMP = setting%scheme%ADMMP
      isADMMS = setting%scheme%ADMMS
      isADMMQ = setting%scheme%ADMMQ
      NbEl = setting%molecule(1)%p%nelectrons
      call mat_init(T23,n2,n3)
      call mat_init(S22,n2,n2)
      call mat_init(S22inv,n2,n2)       

      call get_T23(setting,lupri,luerr,T23,n2,n3,AO2,AO3,GCAO2,GCAO3,&
                  & constrain_factor)
      call get_S22(S22,AO2,GCAO2,setting,lupri,luerr)
      call get_S22inv(S22,S22inv,setting,lupri,luerr)
      
      !    A22 = 2*[ k2(d2) - xc2(d2) ]
      ! or A22 = 2*[ k2(d2) - xc2(d2) ] - LAMBDA_energy S22
      call mat_init(A22,n2,n2)
      call mat_zero(A22)
      call mat_add(2E0_realk,k2,-2E0_realk, xc2, A22)
      IF (isADMMQ.OR.isADMMS.OR.isADMMP) THEN
         ! Assumes the energy has been calculated first - which stores the Lambda
         LambdaE = setting%scheme%ADMM_LARGE_LAMBDA
         IF (isADMMP) THEN
            ! scaling here to avoid scalign too much this Lambda_P s2 contribution
            ! since we scale the whole (k2-x2-lbd_P s2 / xi) by xi**2 afterward
            call mat_daxpy(-2E0_realk*LambdaE/(constrain_factor**(2E0_realk)),S22,A22)
         ELSE
             call mat_daxpy(-2E0_realk*LambdaE,S22,A22)   
         ENDIF   
      ENDIF

      ! B32 = D33 T32 A22 S22inv
      call mat_init(B32,n3,n2)
      call mat_init(tmp32,n3,n2)
      call mat_mul(D3 ,T23,'n','t',1E0_realk,0E0_realk,B32)
      call mat_mul(B32,A22,'n','n',1E0_realk,0E0_realk,tmp32)
      call mat_mul(tmp32,S22inv,'n','n',1E0_realk,0E0_realk,B32)

      ! C22 = D22 A22 S22inv
      call mat_init(C22,n2,n2)
      call mat_init(tmp22,n2,n2)
      call mat_mul(D2 ,A22,'n','n',1E0_realk,0E0_realk,tmp22)
      call mat_mul(tmp22,S22inv,'n','n',1E0_realk,0E0_realk,C22)

      ! 2 correction terms similar to the reorthonormalisation gradient term 
      tmpDFD(1)%p => C22
      call mem_alloc(reOrtho2,3,nAtoms)
      call ls_dzero(reOrtho2,3*nAtoms)
      call II_get_reorthoNormalization_mixed(reOrtho2,tmpDFD,1,AO2,AO2,&
                                       & GCAO2,GCAO2,setting,lupri,luerr)
      tmpDFD(1)%p => B32
      call mem_alloc(reOrtho1,3,nAtoms)
      call ls_dzero(reOrtho1,3*nAtoms)
      call II_get_reorthoNormalization_mixed(reOrtho1,tmpDFD,1,AO3,AO2,&
                                       & GCAO3,GCAO2,setting,lupri,luerr)

      call ls_dzero(ADMM_proj,3*nAtoms)
      
      IF (isADMMQ .OR. isADMMS) THEN
         ! for ADMMQ and ADMMS, transpose(T') s2 T' replaced by S32 T'
         ! so need to compensate for the missing scaling of transpose(T23)
         ! in B32 only, not C22
         call DSCAL(3*nAtoms,constrain_factor,reOrtho1,1)
      ENDIF
      call DAXPY(3*nAtoms, 1E0_realk,reOrtho1,1,ADMM_proj,1) !  Tr(B32 S23e)
      call DAXPY(3*nAtoms,-1E0_realk,reOrtho2,1,ADMM_proj,1) ! -Tr(C22 s2e)
      
      IF (isADMMP) THEN   
         ! scale k2 and x2 by xi^2, BUT (lambda_P s2) scaled by xi *ONLY*
         ! the latter is taken care of earlier in the subroutine
         call DSCAL(3*nAtoms,constrain_factor**(4.E0_realk),ADMM_proj,1)
      ENDIF
      ! -- free memory --
      call mem_dealloc(reOrtho1)
      call mem_dealloc(reOrtho2)
      call mat_free(S22inv)
      call mat_free(S22)
      call mat_free(T23)
      call mat_free(A22)
      call mat_free(B32)
      call mat_free(C22)
      call mat_free(tmp22)
     
      call mat_free(tmp32)
      !
   END SUBROUTINE get_ADMM_K_gradient_projection_term
END SUBROUTINE II_get_ADMM_K_gradient

!> @file
!> Library-like integral-interface routines 

!> \brief Calculates the coulomb matrix
!> \author S. Reine
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param F the coulomb matrix
!> \param ndmat the number of density matrix
SUBROUTINE II_get_coulomb_mat_array(LUPRI,LUERR,SETTING,D,F,ndmat)
IMPLICIT NONE
Integer,intent(in)    :: ndmat
TYPE(MATRIX),intent(in) :: D(ndmat)
TYPE(MATRIX),intent(inout) :: F(ndmat)
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
CALL II_get_coulomb_mat_mixed(LUPRI,LUERR,SETTING,D,F,ndmat,&
     & AORdefault,AORdefault,AORdefault,AORdefault,coulombOperator)
END SUBROUTINE II_get_coulomb_mat_array

!> @file
!> Library-like integral-interface routines 

!> \brief Calculates the coulomb matrix
!> \author S. Reine
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param F the coulomb matrix
!> \param ndmat the number of density matrix
SUBROUTINE II_get_coulomb_mat_single(LUPRI,LUERR,SETTING,D,F,ndmat)
IMPLICIT NONE
Integer,intent(in)      :: ndmat
TYPE(MATRIX),intent(in) :: D
TYPE(MATRIX),intent(inout)  :: F
TYPE(LSSETTING)         :: SETTING
INTEGER                 :: LUPRI,LUERR
!
TYPE(MATRIX)  :: Farray(ndmat)
!NOT OPTIMAL USE OF MEMORY OR ANYTHING
call mat_init(Farray(1),F%nrow,F%ncol)
CALL II_get_coulomb_mat_mixed(LUPRI,LUERR,SETTING,(/D/),Farray(1:1),ndmat,&
     & AORdefault,AORdefault,AORdefault,AORdefault,coulombOperator)
call mat_assign(F,Farray(1))
call mat_free(Farray(1))
END SUBROUTINE II_get_coulomb_mat_single

!> \brief Calculates the coulomb matrix
!> \author S. Reine
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param F the coulomb matrix
!> \param ndmat the number of density matrix
SUBROUTINE II_get_coulomb_mat_mixed(LUPRI,LUERR,SETTING,D,F,ndmat,AO1,AO2,AO3,AO4,Oper)
IMPLICIT NONE
Integer,intent(in)             :: ndmat
TYPE(MATRIX),target,intent(in) :: D(ndmat)
TYPE(MATRIX),intent(inout)     :: F(ndmat)
TYPE(LSSETTING)                :: SETTING
INTEGER,intent(in)             :: LUPRI,LUERR,AO1,AO2,AO3,AO4,Oper
!
Real(realk)           :: TS,TE
integer               :: I
TYPE(MATRIX),target   :: D_AO(ndmat)
TYPE(MATRIX),pointer  :: Dmat(:)
logical :: IntegralTransformGC
!
call time_II_operations1()
CALL LSTIMER('START ',TS,TE,LUPRI)

IF(setting%IntegralTransformGC)THEN
   setting%IntegralTransformGC = .FALSE. 
   !change D to AO basis (currently in GCAO basis)
   DO I=1,ndmat
      CALL mat_init(D_AO(I),D(I)%nrow,D(I)%ncol)
      call GCAO2AO_transform_matrixD2(D(I),D_AO(I),setting,lupri)
   ENDDO
   Dmat => D_AO
   IntegralTransformGC=.TRUE.
ELSE
   Dmat => D
   IntegralTransformGC =.FALSE.
ENDIF

IF (setting%scheme%intprint.GT.10) THEN
  write(lupri,'(A)') 'DMAT in II_get_coulomb_mat_mixed'
  DO I=1,ndmat
    write(lupri,'(1X,A,I3)')  'IDMAT',I
    CALL mat_print(Dmat(I),1,Dmat(1)%ncol,1,Dmat(1)%nrow,lupri)
  ENDDO
ENDIF

IF (SETTING%SCHEME%DENSFIT.OR.SETTING%SCHEME%PARI_J) THEN
  !Consistency testing
  IF ((AO1.NE.AORdefault).OR.(AO2.NE.AORdefault).OR.(AO3.NE.AORdefault).OR.(AO4.NE.AORdefault)&
     &.OR.(Oper.NE.CoulombOperator)) call lsquit('Error in II_get_coulomb_mat_mixed',-1)

  CALL II_get_df_coulomb_mat(LUPRI,LUERR,SETTING,Dmat,F,ndmat)
ELSE
  CALL II_get_regular_coulomb_mat(LUPRI,LUERR,SETTING,Dmat,F,ndmat,AO1,AO2,AO3,AO4,Oper)
ENDIF

IF(IntegralTransformGC)THEN
   DO I=1,ndmat
      CALL mat_free(D_AO(I))
   ENDDO
   !transform back to GCAO basis 
   DO I=1,ndmat
      call AO2GCAO_transform_matrixF(F(I),setting,lupri)
   ENDDO
   setting%IntegralTransformGC = .TRUE. 
ENDIF

IF (setting%scheme%intprint.GT.10) THEN
  write(lupri,'(A)') 'Fmat in II_get_coulomb_mat_mixed'
  DO I=1,ndmat
    write(lupri,'(1X,A,I3)')  'IDMAT',I
    CALL mat_print(F(I),1,F(1)%ncol,1,F(1)%nrow,lupri)
  ENDDO
ENDIF

IF (SETTING%SCHEME%DENSFIT.OR.SETTING%SCHEME%PARI_J) THEN
   IF (SETTING%SCHEME%PARI_J) THEN
      IF (SETTING%SCHEME%SIMPLE_PARI) THEN
         CALL LSTIMER('dfJ-PariSimple',TS,TE,LUPRI)
      ELSE 
         CALL LSTIMER('dfJ-Pari',TS,TE,LUPRI)
      ENDIF
   ELSE IF (SETTING%SCHEME%OVERLAP_DF_J) THEN
      CALL LSTIMER('dfJ-Overlap',TS,TE,LUPRI)
   ELSE
      CALL LSTIMER('dfJ-Jengine',TS,TE,LUPRI)
   ENDIF
ELSE
   IF(SETTING%SCHEME%JENGINE)THEN
      CALL LSTIMER('reg-Jengine',TS,TE,LUPRI)
   ELSE
      CALL LSTIMER('reg-Jbuild',TS,TE,LUPRI)
   ENDIF
ENDIF
call time_II_operations2(JOB_II_get_coulomb_mat)

END SUBROUTINE II_get_coulomb_mat_mixed

!> \brief Calculates the coulomb matrix
!> \author S. Reine
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param F the coulomb matrix
!> \param ndmat the number of density matrix
SUBROUTINE II_get_regular_coulomb_mat(LUPRI,LUERR,SETTING,D,F,ndmat,AO1,AO2,AO3,AO4,Oper)
IMPLICIT NONE
Integer,intent(in)    :: ndmat
TYPE(MATRIX),intent(in)    :: D(ndmat)
TYPE(MATRIX),intent(inout) :: F(ndmat)
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,AO1,AO2,AO3,AO4,Oper
!
Real(realk)         :: fac
integer :: idmat

IF(SETTING%SCHEME%JENGINE)THEN
!   IF(SETTING%SCHEME%FMM)THEN
!     call II_calc_and_write_MMfile_for_FMM(LUPRI,LUERR,SETTING,D)
!     call write_MM_DENS_file(LUPRI,LUERR,SETTING,D)
!     CALL LSTIMER('CalcMM',TS,TE,LUPRI)
!   ENDIF
   call II_get_jengine_mat(LUPRI,LUERR,SETTING,D,F,ndmat,AO1,AO2,AO3,AO4,Oper)
ELSE
  !Consistency testing
  IF ((AO1.NE.AORdefault).OR.(AO2.NE.AORdefault).OR.(AO3.NE.AORdefault).OR.(AO4.NE.AORdefault)&
     &.OR.(Oper.NE.CoulombOperator)) call lsquit('Error in II_get_coulomb_mat_mixed',-1)
  call II_get_4center_coulomb_mat(LUPRI,LUERR,SETTING,D,F,ndmat)
ENDIF

END SUBROUTINE II_get_regular_coulomb_mat

!> \brief Calculates the coulomb matrix using the jengine method (default)
!> \author S. Reine
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param F the coulomb matrix
!> \param ndmat the number of density matrix
SUBROUTINE II_get_jengine_mat(LUPRI,LUERR,SETTING,D,F,ndmat,AO1,AO2,AO3,AO4,Oper)
IMPLICIT NONE
Integer,intent(in)    :: ndmat
TYPE(MATRIX),intent(in)    :: D(ndmat)
TYPE(MATRIX),intent(inout) :: F(ndmat)
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,AO1,AO2,AO3,AO4,Oper
!
Real(realk)         :: TS,TE
integer             :: I

!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%J_THR
!CALL LSHEADER(lupri,'II_get_coulomb_mat')
!CALL LSTIMER('START ',TS,TE,LUPRI)
!ndmat = 1
DO I=1,ndmat
   CALL MAT_ZERO(F(I))
ENDDO
IF(SETTING%SCHEME%FMM)THEN
  IF ((AO1.NE.AORdefault).OR.(AO2.NE.AORdefault).OR.(AO3.NE.AORdefault).OR.(AO4.NE.AORdefault)&
     &.OR.(Oper.NE.CoulombOperator)) call lsquit('Error in II_get_jengine_mat: FMM',-1)
   CALL ls_attachDmatToSetting(D,ndmat,setting,'LHS',1,2,.TRUE.,lupri)
   CALL ls_attachDmatToSetting(D,ndmat,setting,'RHS',3,4,.TRUE.,lupri)
   call ls_multipolemoment(LUPRI,LUERR,SETTING,F(1)%nrow,0,&
        & D(1)%nrow,D(1)%ncol,D(1)%nrow,D(1)%ncol,&
        & AORdefault,AORdefault,AORdefault,AORdefault,RegularSpec,ContractedInttype,.TRUE.)
   CALL ls_freeDmatFromSetting(setting)
ENDIF
CALL ls_attachDmatToSetting(D,ndmat,setting,'RHS',3,4,.TRUE.,lupri)
!Jmat,ndmat,
call initIntegralOutputDims(setting%Output,F(1)%nrow,F(1)%ncol,1,1,ndmat)
call ls_jengine(AO1,AO2,AO3,AO4,Oper,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,F,setting%IntegralTransformGC)
IF(SETTING%SCHEME%FMM.AND.(Oper.EQ.CoulombOperator)) THEN
   IF(ndmat.GT.1)call lsquit('not testet ndmat gt 1 ls_jengineClassicalMAT',-1)
   call ls_jengineClassicalMAT(F(1),AORdefault,AORdefault,AORdefault,AORdefault,&
       & CoulombOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
ENDIF

CALL ls_freeDmatFromSetting(setting)
!CALL LSTIMER('Jengin',TS,TE,LUPRI)

END SUBROUTINE II_get_jengine_mat

!> \brief Calculates the coulomb matrix using explicit 4 center int
!> \author S. Reine
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param F the coulomb matrix
!> \param ndmat the number of density matrix
SUBROUTINE II_get_4center_coulomb_mat(LUPRI,LUERR,SETTING,D,F,ndmat)
IMPLICIT NONE
Integer,intent(in)    :: ndmat
TYPE(MATRIX),intent(in)    :: D(ndmat)
TYPE(MATRIX),intent(inout) :: F(ndmat)
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
!
Real(realk),pointer :: Dfull(:,:,:)
Real(realk),pointer :: Ffull(:,:,:,:,:)
Real(realk)         :: TS,TE
integer             :: I
type(matrixp)       :: Jmat(1),Dmat(1)

!CALL LSHEADER(lupri,'II_get_coulomb_mat')

!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%J_THR
!ndmat = 1
CALL LSTIMER('START ',TS,TE,LUPRI)
call initIntegralOutputDims(setting%Output,F(1)%nrow,F(1)%ncol,1,1,ndmat)
CALL ls_attachDmatToSetting(D,ndmat,setting,'RHS',3,4,.TRUE.,lupri)
call ls_get_coulomb_mat(AORdefault,AORdefault,AORdefault,AORdefault,&
     &                   CoulombOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
!call mem_dealloc(Dfull)
CALL retrieve_Output(lupri,setting,F,setting%IntegralTransformGC)
DO I=1,ndmat
   call mat_scal(2E0_realk,F(I))
ENDDO
CALL ls_freeDmatFromSetting(setting)

CALL LSTIMER('reg-J ',TS,TE,LUPRI)

END SUBROUTINE II_get_4center_coulomb_mat

!> \brief Calculates the exchange matrix
!> \author S. Reine
!> \date   May 2012
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param ndmat the number of density matrices
!> \param Dsym is the density symmetric
!> \param F the exchange matrix
SUBROUTINE II_get_exchange_mat_array(LUPRI,LUERR,SETTING,D,ndmat,Dsym,F)
IMPLICIT NONE
Integer,intent(in)    :: ndmat
TYPE(MATRIX),intent(in)    :: D(ndmat)
TYPE(MATRIX),intent(inout) :: F(ndmat)
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
LOGICAL,intent(in)    :: Dsym

CALL II_get_exchange_mat_mixed(LUPRI,LUERR,SETTING,D,ndmat,Dsym,F,&
     &    AORdefault,AORdefault,AORdefault,AORdefault,coulombOperator)

END SUBROUTINE II_get_exchange_mat_array

!> \brief Calculates the exchange matrix
!> \author S. Reine
!> \date   May 2012
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param ndmat the number of density matrices
!> \param Dsym is the density symmetric
!> \param F the exchange matrix
SUBROUTINE II_get_exchange_mat_single(LUPRI,LUERR,SETTING,D,ndmat,Dsym,F)
IMPLICIT NONE
Integer,intent(in)    :: ndmat
TYPE(MATRIX),intent(in)    :: D
TYPE(MATRIX),intent(inout) :: F
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
LOGICAL,intent(in)    :: Dsym
!
TYPE(MATRIX)  :: Farray(ndmat)
!NOT OPTIMAL USE OF MEMORY OR ANYTHING
call mat_init(Farray(1),F%nrow,F%ncol)
call mat_zero(Farray(1))
CALL II_get_exchange_mat_mixed(LUPRI,LUERR,SETTING,(/D/),ndmat,Dsym,Farray,&
     &    AORdefault,AORdefault,AORdefault,AORdefault,coulombOperator)
call mat_assign(F,Farray(1))
call mat_free(Farray(1))
END SUBROUTINE II_get_exchange_mat_single


!> \brief Calculates the exchange matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param ndmat the number of density matrices
!> \param Dsym is the density symmetric
!> \param F the exchange matrix
SUBROUTINE ii_get_exchange_mat_mixed(LUPRI,LUERR,SETTING,D,ndmat,Dsym,F,AO1,AO3,AO2,AO4,Oper)
IMPLICIT NONE
Integer,intent(in)    :: ndmat
TYPE(MATRIX),intent(in)    :: D(ndmat)
TYPE(MATRIX),intent(inout) :: F(ndmat)
TYPE(LSSETTING)       :: SETTING
INTEGER,intent(in)    :: LUPRI,LUERR,AO1,AO2,AO3,AO4,Oper
LOGICAL               :: Dsym
!
integer :: isym(ndmat),idmat,idmat2,ndmat2
type(matrix),pointer :: D2(:),F2(:)
IF (SETTING%SCHEME%exchangeFactor.EQ. 0.0E0_realk) RETURN
call time_II_operations1()

! Check symetry. Split non-symmetric matrices to symmetric and anti symmetric parts. 
! Make symmetry check
IF(.NOT.DSYM)THEN
   idmat2 = 0
   do idmat=1,ndmat
      ISYM(idmat) = mat_get_isym(D(idmat))
      IF(ISYM(idmat).EQ.1.OR.ISYM(idmat).EQ.2)THEN !sym or antisym
         idmat2 = idmat2 + 1
      ELSEIF(ISYM(idmat).EQ.3)THEN !nonsym
         idmat2 = idmat2 + 2 !we split up into sym and anti sym part
      ELSE
         !do nothing
      ENDIF
   enddo
   ndmat2 = idmat2
   if(ndmat2.GT.0)THEN
      call mem_alloc(D2,ndmat2)
      do idmat2=1,ndmat2
         call mat_init(D2(idmat2),D(1)%nrow,D(1)%ncol)
      enddo
      idmat2 = 0
      do idmat=1,ndmat
         IF(ISYM(idmat).EQ.1.OR.ISYM(idmat).EQ.2)THEN !sym or antisym
            idmat2 = idmat2 + 1
            call mat_copy(1.0E0_realk,D(idmat), D2(idmat2))
         ELSEIF(ISYM(idmat).EQ.3)THEN !nonsym
            call mat_assign(D2(idmat2+1),D(idmat))
            call util_get_symm_part(D2(idmat2+1))
            call util_get_antisymm_part(D(idmat),D2(idmat2+2))
            idmat2 = idmat2 + 2 !we split up into sym and anti sym part
         ELSE
            !do nothing
         ENDIF
      enddo
      call mem_alloc(F2,ndmat2)
      do idmat2=1,ndmat2
         call mat_init(F2(idmat2),F(1)%nrow,F(1)%ncol)
         call mat_zero(F2(idmat2)) 
      enddo
      call II_get_exchange_mat1(LUPRI,LUERR,SETTING,D2,ndmat2,F2,AO1,AO3,AO2,AO4,Oper)
      do idmat2=1,ndmat2
         call mat_free(D2(idmat2))
      enddo
      call mem_dealloc(D2)
      idmat2 = 0
      do idmat=1,ndmat
         IF(ISYM(idmat).EQ.1.OR.ISYM(idmat).EQ.2)THEN !sym or antisym
            idmat2 = idmat2 + 1
            call mat_copy(1.0E0_realk,F2(idmat2), F(idmat))
         ELSEIF(ISYM(idmat).EQ.3)THEN !nonsym
            call mat_add(1E0_realk,F2(idmat2+1),1E0_realk,F2(idmat2+2),F(idmat))
            idmat2 = idmat2 + 2
         ELSE
            call mat_zero(F(idmat))
         ENDIF
      enddo
      do idmat2=1,ndmat2
         call mat_free(F2(idmat2))
      enddo
      call mem_dealloc(F2)
   else
      do idmat=1,ndmat
         call mat_zero(F(idmat))
      enddo
   endif
ELSE
   call II_get_exchange_mat1(LUPRI,LUERR,SETTING,D,ndmat,F,AO1,AO3,AO2,AO4,Oper)
ENDIF
call time_II_operations2(JOB_II_GET_EXCHANGE_MAT)

END SUBROUTINE ii_get_exchange_mat_mixed

!> \brief Calculates the exchange matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param ndmat the number of density matrices
!> \param F the exchange matrix
SUBROUTINE II_get_exchange_mat1(LUPRI,LUERR,SETTING,D,ndmat,F,AO1,AO3,AO2,AO4,Oper)
IMPLICIT NONE
Integer,intent(in)    :: ndmat
TYPE(MATRIX),target,intent(in)    :: D(ndmat)
TYPE(MATRIX),intent(inout) :: F(ndmat)
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,AO1,AO2,AO3,AO4,Oper
!
TYPE(MATRIX),target   :: D_AO(ndmat)
TYPE(MATRIX),pointer  :: Dmat(:)
logical :: IntegralTransformGC
integer :: I

IF(setting%IntegralTransformGC)THEN
   setting%IntegralTransformGC = .FALSE. 
   !change D to AO basis (currently in GCAO basis)
   DO I=1,ndmat
      CALL mat_init(D_AO(I),D(I)%nrow,D(I)%ncol)
      call GCAO2AO_transform_matrixD2(D(I),D_AO(I),setting,lupri)
   ENDDO
   Dmat => D_AO
   IntegralTransformGC = .TRUE.
ELSE
   Dmat => D
   IntegralTransformGC = .FALSE.
ENDIF

IF (SETTING%SCHEME%DF_K) THEN
   IF ((AO1.NE.AORdefault).OR.(AO2.NE.AORdefault).OR.&
        (AO3.NE.AORdefault).OR.(AO4.NE.AORdefault)) &
        CALL LSQUIT('Error in II_get_exchange_mat1. DF_K is only implemented for regular AOs',-1)

   CALL II_get_df_exchange_mat(LUPRI,LUERR,SETTING,Dmat,F,ndmat)
ELSE IF (SETTING%SCHEME%PARI_K.OR.SETTING%SCHEME%MOPARI_K) THEN
   IF ((AO1.NE.AORdefault).OR.(AO2.NE.AORdefault).OR.&
        (AO3.NE.AORdefault).OR.(AO4.NE.AORdefault)) &
        & CALL LSQUIT('Error in II_get_exchange_mat1. PARI_K is only implemented for regular AOs',-1)
   CALL II_get_pari_df_exchange_mat(LUPRI,LUERR,SETTING,Dmat,F,ndmat)
ELSE
   IF (((AO1.NE.AORdefault).OR.(AO2.NE.AORdefault).OR.&
        (AO3.NE.AORdefault).OR.(AO4.NE.AORdefault) &
        & .OR.(Oper.NE.coulombOperator)).AND.setting%IntegralTransformGC) &
        & CALL LSQUIT('Error in II_get_exchange_mat1. Mixed exchange does not work with GC basis',-1)
   CALL II_get_exchange_mat_regular(LUPRI,LUERR,SETTING,Dmat,ndmat,F,AO1,AO3,AO2,AO4,Oper)
ENDIF

IF(IntegralTransformGC)THEN
   DO I=1,ndmat
      CALL mat_free(D_AO(I))
   ENDDO
   !transform back to GCAO basis 
   DO I=1,ndmat
      call AO2GCAO_transform_matrixF(F(I),setting,lupri)
   ENDDO
   setting%IntegralTransformGC = .TRUE. 
ENDIF

END SUBROUTINE II_get_exchange_mat1

!> \brief Calculates the exchange matrix using explicit 4 center 
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param ndmat the number of density matrices
!> \param F the exchange matrix
SUBROUTINE II_get_exchange_mat_regular(LUPRI,LUERR,SETTING,D,ndmat,F,AO1,AO3,AO2,AO4,OperIn)
IMPLICIT NONE
Integer,intent(in)    :: ndmat
TYPE(MATRIX),target,intent(in)    :: D(ndmat)
TYPE(MATRIX),intent(inout) :: F(ndmat)
TYPE(MATRIX),pointer  :: K(:)
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR,AO1,AO2,AO3,AO4,OperIn
!
Integer             :: idmat,incdmat,nrow,ncol,ndmats
Real(realk),pointer :: Kfull(:,:,:,:,:)
Real(realk),pointer :: DfullLHS(:,:,:)
Real(realk),pointer :: DfullRHS(:,:,:)
Real(realk)         :: TS,TE,fac
integer    :: Oper
integer :: nCalcInt,nCalcIntZero,nCalcIntZeroContrib,i

CALL LSTIMER('START ',TS,TE,LUPRI)
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%K_THR
IF(matrix_type .EQ. mtype_unres_dense) THEN
  ndmats = 2*ndmat
ELSE
  ndmats = ndmat
ENDIF

!CALL LSHEADER(lupri,'II_get_exchange_mat')
CALL ls_attachDmatToSetting(D,ndmat,setting,'RHS',2,4,.FALSE.,lupri)

IF (SETTING%SCHEME%CAM) THEN
  Oper = CAMOperator       !Coulomb attenuated method
ELSEIF (SETTING%SCHEME%SR_EXCHANGE) THEN
  Oper = ErfcOperator      !Short-Range Coulomb screened exchange
ELSE
  Oper = OperIn   !Deafult is the Coulomb metric 
ENDIF
!Calculates the HF-exchange contribution
call initIntegralOutputDims(setting%Output,F(1)%nrow,F(1)%ncol,1,1,ndmats)
call ls_get_exchange_mat(AO1,AO3,AO2,AO4,Oper,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
!IF (setting%scheme%daLinK)call mem_dealloc(DfullLHS)
!call mem_dealloc(DfullRHS)
call mem_alloc(K,ndmat)
DO idmat=1,ndmat
   call mat_init(K(idmat),F(1)%nrow,F(1)%ncol)
ENDDO
CALL retrieve_Output(lupri,setting,K,setting%IntegralTransformGC)

DO idmat=1,ndmat
  call mat_daxpy(-1E0_realk,K(idmat),F(idmat))
  call mat_free(K(idmat))
ENDDO
call mem_dealloc(K)

CALL ls_freeDmatFromSetting(setting)

IF(SETTING%SCHEME%DALINK .AND.SETTING%SCHEME%LINK)THEN
   CALL LSTIMER('DaLINK-Kbuild',TS,TE,LUPRI)
ELSEIF(SETTING%SCHEME%LINK)THEN
   CALL LSTIMER('LINK-Kbuild',TS,TE,LUPRI)
ELSE
   CALL LSTIMER('st-Kbuild',TS,TE,LUPRI)
ENDIF
END SUBROUTINE II_get_exchange_mat_regular

!> \brief Calculates the coulomb and exchange matrix using explicit 4 center 
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param F the exchange matrix
!> \param ndmat the number of density matrices
SUBROUTINE II_get_coulomb_and_exchange_mat_array(LUPRI,LUERR,SETTING,D,F,ndmat)
IMPLICIT NONE
Integer,intent(in)    :: ndmat
TYPE(MATRIX),intent(in)    :: D(ndmat)
TYPE(MATRIX),intent(inout) :: F(ndmat)
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
!
Real(realk) :: TS,TE
Integer     :: idmat
call time_II_operations1()

CALL LSTIMER('START ',TS,TE,LUPRI)
!set threshold 
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%K_THR

!CALL LSHEADER(lupri,'II_get_coulomb_and_exchange_mat')
!ndmat = 1
CALL ls_attachDmatToSetting(D,ndmat,setting,'RHS',3,4,.TRUE.,lupri)
call initIntegralOutputDims(setting%Output,F(1)%nrow,F(1)%ncol,1,1,ndmat)
call ls_get_coulomb_and_exchange_mat(AORdefault,AORdefault,AORdefault,AORdefault,&
     &                   CoulombOperator,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,F,setting%IntegralTransformGC)
CALL ls_freeDmatFromSetting(setting)
DO idmat=1,ndmat
  CALL mat_scal(2E0_realk,F(idmat))
ENDDO
CALL LSTIMER('Joint-JK  ',TS,TE,LUPRI)
call time_II_operations2(JOB_II_get_coulomb_and_exchange_mat)

END SUBROUTINE II_get_coulomb_and_exchange_mat_array

!> \brief Calculates the coulomb and exchange matrix using explicit 4 center 
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param F the exchange matrix
!> \param ndmat the number of density matrices
SUBROUTINE II_get_coulomb_and_exchange_mat_single(LUPRI,LUERR,SETTING,D,F,ndmat)
IMPLICIT NONE
Integer,intent(in)    :: ndmat
TYPE(MATRIX),intent(in)    :: D
TYPE(MATRIX),intent(inout) :: F
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
TYPE(MATRIX)  :: Farray(ndmat)
!NOT OPTIMAL USE OF MEMORY OR ANYTHING
call mat_init(Farray(1),F%nrow,F%ncol)
call II_get_coulomb_and_exchange_mat_array(LUPRI,LUERR,SETTING,(/D/),Farray,ndmat)
call mat_assign(F,Farray(1))
call mat_free(Farray(1))
END SUBROUTINE II_get_coulomb_and_exchange_mat_single

!> \brief Calculates the Fock matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param Dsym is the density symmetric
!> \param F the exchange matrix
!> \param ndmat the number of density matrices
SUBROUTINE II_get_Fock_mat_array(LUPRI,LUERR,SETTING,D,Dsym,F,ndmat,setting_incremental_scheme)
IMPLICIT NONE
integer,intent(in)          :: ndmat
TYPE(MATRIX),intent(in)     :: D(ndmat)
TYPE(MATRIX),intent(inout)  :: F(ndmat)
TYPE(LSSETTING),intent(inout)  :: SETTING
INTEGER,intent(in)    :: LUPRI,LUERR
LOGICAL,intent(in)    :: Dsym,setting_incremental_scheme
!
TYPE(MATRIX)          :: D_AO(ndmat),K(ndmat)
logical               :: default,inc_scheme,do_inc
real(realk)           :: TS,TE,fac,maxelm
integer               :: i
default = SETTING%SCHEME%DENSFIT.OR.SETTING%SCHEME%PARI_J.OR.SETTING%SCHEME%JENGINE &
     &      .OR. SETTING%SCHEME%LinK .OR. (matrix_type .EQ. mtype_unres_dense).OR. &
     &      SETTING%SCHEME%CAM

fac = 2E0_realk
IF(matrix_type .EQ. mtype_unres_dense)fac = 1E0_realk
call get_incremental_settings(inc_scheme,do_inc)
!setting_incremental_scheme   from setting 
!incremental_scheme = inc_scheme = input keyword .INCREM activated
!do_increment       = do_inc     = increm activated 
IF(setting_incremental_scheme.AND..NOT.inc_scheme)THEN
   call lsquit('II_get_Fock_mat incremental scheme inconsistensy',lupri)
ENDIF
IF(setting_incremental_scheme.AND.(ndmat.GT.1.AND.(matrix_type .NE. mtype_unres_dense)))THEN
   call lsquit('II_get_Fock_mat incremental scheme inconsistensy2',lupri)
ENDIF

IF(setting%IntegralTransformGC)THEN
   setting%IntegralTransformGC = .FALSE. !we do this manually in order to
                                         !get incremental correct
   !change D to AO basis (currently in GCAO basis)
   DO I=1,ndmat
      CALL mat_init(D_AO(I),D(1)%nrow,D(1)%ncol)
      call GCAO2AO_transform_matrixD2(D(I),D_AO(I),setting,lupri)
   ENDDO
   IF(setting_incremental_scheme)THEN
      call mat_assign(incrDdiff(1),D_AO(1))
      IF(do_inc)THEN !do increment
         call mat_daxpy(-1E0_realk,incrD0(1),incrDdiff(1)) !incrDdiff is now the difference density
      ELSE
         !test if we should activate increm
         call mat_daxpy(-1E0_realk,incrD0(1),incrDdiff(1)) !incrDdiff is now the difference density 
         call activate_incremental(lupri,do_inc)
         IF(.NOT.do_inc)THEN
            !we should not activate increm so we set incrDdiff = D_AO
            call mat_assign(incrDdiff(1),D_AO(1))
         ENDIF
      ENDIF
      IF (default) THEN
         call II_get_coulomb_mat(LUPRI,LUERR,SETTING,incrDdiff,F,ndmat)   
         DO I=1,ndmat
            IF(do_inc.AND.setting_incremental_scheme)THEN ! dot(D_AO * J_AO(incrDdiff))
               WRITE(lupri,*)'The increment Coulomb energy contribution ',&
                    & fac*0.5E0_realk*mat_dotproduct(D_AO(I),F(I))
            ELSE           ! dot(D_AO * J_AO(D_AO))
               WRITE(lupri,*)'The Coulomb energy contribution ',&
                    & fac*0.5E0_realk*mat_dotproduct(D_AO(I),F(I))
            ENDIF
         ENDDO
         DO i=1,ndmat
            call mat_init(K(i),F(1)%nrow,F(1)%ncol)
            call mat_zero(K(i))
         ENDDO
         IF (setting%scheme%daLinK) THEN
            !attach full (NOT diff D) to LHS
            CALL ls_attachDmatToSetting(D_AO,ndmat,setting,'LHS',2,4,.FALSE.,lupri)
         ENDIF
         call II_get_exchange_mat(LUPRI,LUERR,SETTING,incrDdiff,ndmat,Dsym,K)
         DO i=1,ndmat
            IF(do_inc.AND.setting_incremental_scheme)THEN ! dot(D_AO * K_AO(incrDdiff))
               WRITE(lupri,*)'The increment Exchange energy contribution ',&
                    & mat_dotproduct(D_AO(i),K(i))
            ELSE           ! dot(D_AO * K_AO(D_AO)) 
               WRITE(lupri,*)'The Exchange energy contribution ',&
                    &mat_dotproduct(D_AO(i),K(i))
            ENDIF
            call mat_daxpy(1E0_realk,K(i),F(i))
            call mat_free(K(i))
         ENDDO
      ELSE 
         call II_get_coulomb_and_exchange_mat(LUPRI,LUERR,SETTING,incrDdiff,F,ndmat)
         DO i=1,ndmat
            IF(do_inc.AND.setting_incremental_scheme)THEN ! dot(D_AO * F_AO(incrDdiff))
               WRITE(lupri,*)'The increment J+K energy contribution ',&
                    &-mat_dotproduct(D_AO(i),F(i))
            ELSE           ! dot(D_AO * F_AO(D_AO))
               WRITE(lupri,*)'The J+K energy contribution ',&
                    &-mat_dotproduct(D_AO(i),F(i))
            ENDIF
         ENDDO
      ENDIF
      IF(do_inc)THEN !obtain full AO fock matric F = incrF0 + incremental F 
         call mat_daxpy(1E0_realk,incrF0(1),F(1))
      ENDIF
   ELSE
      IF (default) THEN
         call II_get_coulomb_mat(LUPRI,LUERR,SETTING,D_AO,F,ndmat)   
         DO I=1,ndmat
            WRITE(lupri,*)'The Coulomb energy contribution ',&
                 & fac*0.5E0_realk*mat_dotproduct(D_AO(I),F(I))
         ENDDO
         DO i=1,ndmat
            call mat_init(K(i),F(1)%nrow,F(1)%ncol)
            call mat_zero(K(i))
         ENDDO
         IF (setting%scheme%daLinK) THEN
            CALL ls_attachDmatToSetting(D_AO,ndmat,setting,'LHS',2,4,.FALSE.,lupri)
         ENDIF
         call II_get_exchange_mat(LUPRI,LUERR,SETTING,D_AO,ndmat,Dsym,K)
         DO i=1,ndmat
            WRITE(lupri,*)'The Exchange energy contribution ',&
                 &mat_dotproduct(D_AO(i),K(i))
            call mat_daxpy(1E0_realk,K(i),F(i))
            call mat_free(K(i))
         ENDDO
      ELSE 
         call II_get_coulomb_and_exchange_mat(LUPRI,LUERR,SETTING,D_AO,F,ndmat)
         DO i=1,ndmat
               WRITE(lupri,*)'The J+K energy contribution ',&
                    &-mat_dotproduct(D_AO(i),F(i))
         ENDDO
      ENDIF
   ENDIF
   IF(SaveF0andD0)THEN
      call mat_assign(incrD0(1),D_AO(1))
      call mat_assign(incrF0(1),F(1))
   ENDIF
   !transform back to GCAO basis 
   DO i=1,ndmat
      call AO2GCAO_transform_matrixF(F(I),setting,lupri)
   ENDDO
   DO I=1,ndmat
      CALL mat_free(D_AO(I))
   ENDDO
   setting%IntegralTransformGC = .TRUE. !back to original value
ELSE
   IF(setting_incremental_scheme)THEN
      call mat_assign(incrDdiff(1),D(1))
      IF(do_inc)THEN !do increment
         call mat_daxpy(-1E0_realk,incrD0(1),incrDdiff(1)) !incrDdiff is now the difference density 
      ELSE
         !test is we should activate increm
         call mat_daxpy(-1E0_realk,incrD0(1),incrDdiff(1)) !incrDdiff is now the difference density 
         call activate_incremental(lupri,do_inc)
         IF(.NOT.do_inc)THEN
            !we should not activate increm so we set incrDdiff = D_AO
            call mat_assign(incrDdiff(1),D(1))
         ENDIF
      ENDIF
      IF (default) THEN
         call II_get_coulomb_mat(LUPRI,LUERR,SETTING,incrDdiff,F,ndmat)   
         DO I=1,ndmat
            IF(do_inc.AND.setting_incremental_scheme)THEN ! dot(D * J(incrDdiff))
               WRITE(lupri,*)'The increment Coulomb energy contribution ',&
                    &fac*0.5E0_realk*mat_dotproduct(D(I),F(I))
            ELSE           ! dot(D * J(D))
               WRITE(lupri,*)'The Coulomb energy contribution ',fac*0.5E0_realk*mat_dotproduct(D(I),F(I))
            ENDIF
         ENDDO
         DO i=1,ndmat
            call mat_init(K(i),F(I)%nrow,F(I)%ncol)
            call mat_zero(K(i))
         ENDDO
         IF (setting%scheme%daLinK) THEN
            !attach full D (NOT the diff D)
            CALL ls_attachDmatToSetting(D,ndmat,setting,'LHS',2,4,.FALSE.,lupri)
         ENDIF
         call II_get_exchange_mat(LUPRI,LUERR,SETTING,incrDdiff,ndmat,Dsym,K)
         DO i=1,ndmat
            IF(do_inc.AND.setting_incremental_scheme)THEN ! dot(D * K(incrDdiff))
               WRITE(lupri,*)'The increment Exchange energy contribution ',-mat_dotproduct(D(i),K(i))
            ELSE           ! dot(D * K(D))
               WRITE(lupri,*)'The Exchange energy contribution ',-mat_dotproduct(D(i),K(i))
            ENDIF
            call mat_daxpy(1E0_realk,K(i),F(i))
            call mat_free(K(i))
         ENDDO
      ELSE 
         call II_get_coulomb_and_exchange_mat(LUPRI,LUERR,SETTING,D,F,ndmat)
         DO i=1,ndmat
            IF(do_inc.AND.setting_incremental_scheme)THEN ! dot(D * F(incrDdiff))
               WRITE(lupri,*)'The increment J+K energy contribution ',-mat_dotproduct(D(i),F(i))
            ELSE           ! dot(D * F(D))
               WRITE(lupri,*)'The J+K energy contribution ',-mat_dotproduct(D(i),F(i))
            ENDIF
         ENDDO
      ENDIF
      IF(do_inc)THEN !obtain full AO fock matric F = incrF0 + incremental F 
         call mat_daxpy(1E0_realk,incrF0(1),F(1))
      ENDIF
   ELSE
      IF (default) THEN
         call II_get_coulomb_mat(LUPRI,LUERR,SETTING,D,F,ndmat)   
         DO I=1,ndmat
            WRITE(lupri,*)'The Coulomb energy contribution ',fac*0.5E0_realk*mat_dotproduct(D(I),F(I))
         ENDDO
         DO i=1,ndmat
            call mat_init(K(i),F(I)%nrow,F(I)%ncol)
            call mat_zero(K(i))
         ENDDO
         IF (setting%scheme%daLinK) THEN
            CALL ls_attachDmatToSetting(D,ndmat,setting,'LHS',2,4,.FALSE.,lupri)
         ENDIF
         call II_get_exchange_mat(LUPRI,LUERR,SETTING,D,ndmat,Dsym,K)
         DO i=1,ndmat
            WRITE(lupri,*)'The Exchange energy contribution ',-mat_dotproduct(D(i),K(i))
            call mat_daxpy(1E0_realk,K(i),F(i))
            call mat_free(K(i))
         ENDDO
      ELSE 
         call II_get_coulomb_and_exchange_mat(LUPRI,LUERR,SETTING,D,F,ndmat)
         DO i=1,ndmat
            WRITE(lupri,*)'The J+K energy contribution ',-mat_dotproduct(D(i),F(i))
         ENDDO
      ENDIF
   ENDIF
   IF(SaveF0andD0)THEN
      call mat_assign(incrD0(1),D(1))
      call mat_assign(incrF0(1),F(1))
   ENDIF
ENDIF
DO i=1,ndmat
   WRITE(lupri,*)'The Fock energy contribution ',-mat_dotproduct(D(i),F(i))
ENDDO

END SUBROUTINE II_get_Fock_mat_array

!> \brief Calculates the Fock matrix
!> \author T. Kjaergaard
!> \date 2010
!> \param lupri Default print unit
!> \param luerr Default error print unit
!> \param setting Integral evalualtion settings
!> \param D the density matrix
!> \param Dsym is the density symmetric
!> \param F the exchange matrix
!> \param ndmat the number of density matrices
SUBROUTINE II_get_Fock_mat_single(LUPRI,LUERR,SETTING,D,Dsym,F,ndmat,setting_incremental_scheme)
IMPLICIT NONE
integer,intent(in)          :: ndmat
TYPE(MATRIX),intent(in)     :: D
TYPE(MATRIX),intent(inout)  :: F
TYPE(LSSETTING),intent(inout)  :: SETTING
INTEGER,intent(in)    :: LUPRI,LUERR
LOGICAL,intent(in)    :: Dsym,setting_incremental_scheme
TYPE(MATRIX)  :: Farray(ndmat)
!NOT OPTIMAL USE OF MEMORY OR ANYTHING
call mat_init(Farray(1),F%nrow,F%ncol)
call II_get_Fock_mat_array(LUPRI,LUERR,SETTING,(/D/),Dsym,Farray,ndmat,setting_incremental_scheme)
call mat_assign(F,Farray(1))
call mat_free(Farray(1))
END SUBROUTINE II_get_Fock_mat_single

!****ADMM specific routines starts here
!Due to cross dependence on II_routines they need to stay in this module
SUBROUTINE Transformed_F2_to_F3(F,F2,setting,lupri,luerr,n2,n3,AO2,AO3,&
                              & GCAO2,GCAO3,constrain_factor)
  implicit none
  type(matrix),intent(inout) :: F  !level 3 matrix output 
  type(matrix),intent(in)    :: F2 !level 2 matrix input 
  type(lssetting) :: setting
  Integer :: n2,n3,AO2,AO3,lupri,luerr
  Logical :: GCAO2,GCAO3
  real(realk),intent(IN)        :: constrain_factor     
  !
  TYPE(MATRIX) :: S23,T23

  CALL mat_init(T23,n2,n3)
  CALL mat_init(S23,n2,n3)

  CALL get_T23(setting,lupri,luerr,T23,n2,n3,&
               & AO2,AO3,GCAO2,GCAO3,constrain_factor)

  CALL mat_mul(F2,T23,'n','n',1E0_realk,0E0_realk,S23)
  CALL mat_mul(T23,S23,'t','n',1E0_realk,0E0_realk,F)
 
  CALL mat_free(T23)
  CALL mat_free(S23)
END SUBROUTINE TRANSFORMED_F2_TO_F3

SUBROUTINE get_small_lambda(constrain_factor,D3,setting,lupri,luerr,n2,n3,&
                  & AO2,AO3,GCAO2,GCAO3)
   implicit none
   type(matrix),intent(in)    :: D3     !level 3 matrix input 
   real(realk),intent(inout)  :: constrain_factor
   type(lssetting)            :: setting
   integer                    :: n2,n3,AO3,AO2,lupri,luerr
   logical                    :: GCAO2,GCAO3
   !
   TYPE(MATRIX) :: R33
   real(realk)  :: lambda,trace,fac
   integer      :: nelectrons
   logical      :: unrest

   unrest = matrix_type.EQ.mtype_unres_dense

   CALL mat_init(R33,n3,n3)

   CALL get_R33(R33,AO2,AO3,GCAO3,n2,n3,1E0_realk,setting,lupri,luerr)
   trace = mat_trAB(D3,R33)

   nelectrons = setting%molecule(1)%p%nelectrons

   ! The lagrangian multiplier
   ! lambda = 1 - sqrt[ 2/N Tr(D3 S32 T23) ] 
   fac = 2.0E0_realk
   IF (unrest) fac = 1.0E0_realk
   lambda = 1E0_realk - sqrt(fac*trace/nelectrons)
   constrain_factor = 1.0E0_realk / (1E0_realk - lambda)

   CALL mat_free(R33)

END SUBROUTINE get_small_lambda

SUBROUTINE get_T23(setting,lupri,luerr,T23,n2,n3,&
                  & AO2,AO3,GCAO2,GCAO3,constrain_factor)
use io
implicit none
TYPE(lssetting),intent(inout) :: setting
TYPE(MATRIX),intent(inout)    :: T23
Integer,intent(IN)            :: n2,n3,AO2,AO3,lupri,luerr
Logical,intent(IN)            :: GCAO2,GCAO3
real(realk),intent(IN)        :: constrain_factor
!
TYPE(MATRIX) :: S23,S22,S22inv
Character(80) :: Filename
Logical :: McWeeny,ERI2C
real(realk) :: lambda
Logical     :: isADMMQ,isADMMS,isADMMP
!
isADMMQ = setting%scheme%ADMMQ
isADMMS = setting%scheme%ADMMS
isADMMP = setting%scheme%ADMMP
!these options are for the ERI metric
!with McWeeny ADMM1 is assumed, without ADMM2
McWeeny = setting%scheme%ADMM1
ERI2C   = setting%scheme%ADMM_2ERI

write(Filename,'(A8,2L1)') 'ADMM_T23',GCAO2,GCAO3

IF (io_file_exist(Filename,setting%IO)) THEN
  call io_read_mat(T23,Filename,setting%IO,LUPRI,LUERR)
ELSE
  CALL mat_init(S22,n2,n2)
  CALL mat_init(S22inv,n2,n2)
  CALL mat_init(S23,n2,n3)
 
  call get_S22(S22,AO2,GCAO2,setting,lupri,luerr)
  call get_S23(S23,AO2,AO3,GCAO2,GCAO3,setting,lupri,luerr)
  call get_S22inv(S22,S22inv,setting,lupri,luerr)
  CALL mat_mul(S22inv,S23,'n','n',1E0_realk,0E0_realk,T23)

  CALL mat_free(S22inv)
  CALL mat_free(S23)
  CALL mat_free(S22)
  call io_add_filename(setting%IO,Filename,LUPRI)
  call io_write_mat(T23,Filename,setting%IO,LUPRI,LUERR)
ENDIF
! IF constraining the total charge
! Lagrangian multiplier for conservation of the total nb. of electrons
! constrain_factor = 1 / (1-lambda)
IF (isADMMQ .OR. isADMMS) THEN
   call mat_scal(constrain_factor,T23)
ENDIF
END SUBROUTINE get_T23

SUBROUTINE get_S22inv(S22,S22inv,setting,lupri,luerr)
implicit none
TYPE(matrix),intent(INOUT)    :: S22
TYPE(matrix),intent(INOUT)    :: S22inv
TYPE(lssetting),intent(inout) :: setting
INTEGER,intent(IN)            :: lupri,luerr
!
Character(80) :: Filename

Filename='ADMM_S22inv'
IF (io_file_exist(Filename,setting%IO)) THEN
  call io_read_mat(S22inv,Filename,setting%IO,LUPRI,LUERR)
ELSE
  CALL mat_inv(S22,S22inv)
ENDIF
END SUBROUTINE get_S22inv

SUBROUTINE get_S22(S22,AO2,GCAO2,setting,lupri,luerr)
implicit none
TYPE(matrix),intent(INOUT)    :: S22
TYPE(lssetting),intent(inout) :: setting
Integer,intent(IN)            :: AO2,lupri,luerr
Logical,intent(IN)            :: GCAO2
!
Character(80) :: Filename
Logical :: McWeeny,ERI2C
McWeeny = setting%scheme%ADMM1
ERI2C   = setting%scheme%ADMM_2ERI
Filename='ADMM_S22'
IF (io_file_exist(Filename,setting%IO)) THEN
  call io_read_mat(S22,Filename,setting%IO,LUPRI,LUERR)
ELSE
  IF (ERI2C.AND..NOT.McWeeny) THEN
    CALL II_get_2center_mixed_eri(lupri,luerr,setting,S22,AO2,AO2,GCAO2,GCAO2)
  ELSE
    CALL II_get_mixed_overlap(lupri,luerr,setting,S22,AO2,AO2,GCAO2,GCAO2)
  ENDIF
ENDIF
END SUBROUTINE get_S22

SUBROUTINE get_S33(S33,AO3,GCAO3,setting,lupri,luerr)
implicit none
TYPE(matrix),intent(INOUT)    :: S33
TYPE(lssetting),intent(inout) :: setting
Integer,intent(IN)            :: AO3,lupri,luerr
Logical,intent(IN)            :: GCAO3
!
Character(80) :: Filename
write(Filename,'(A8,L1)') 'ADMM_S33',GCAO3
IF (io_file_exist(Filename,setting%IO)) THEN
  call io_read_mat(S33,Filename,setting%IO,LUPRI,LUERR)
ELSE
  CALL II_get_mixed_overlap(lupri,luerr,setting,S33,AO3,AO3,GCAO3,GCAO3)
ENDIF
END SUBROUTINE get_S33

SUBROUTINE get_S23(S23,AO2,AO3,GCAO2,GCAO3,setting,lupri,luerr)
implicit none
TYPE(matrix),intent(INOUT)    :: S23
TYPE(lssetting),intent(inout) :: setting
Integer,intent(IN)            :: AO2,AO3,lupri,luerr
Logical,intent(IN)            :: GCAO2,GCAO3
!
Character(80) :: Filename
Logical :: McWeeny,ERI2C
McWeeny = setting%scheme%ADMM1
ERI2C   = setting%scheme%ADMM_2ERI
write(Filename,'(A8,L1)') 'ADMM_S23',GCAO3
IF (io_file_exist(Filename,setting%IO)) THEN
  call io_read_mat(S23,Filename,setting%IO,LUPRI,LUERR)
ELSE
  IF (ERI2C) THEN
    CALL II_get_2center_mixed_eri(lupri,luerr,setting,S23,AO2,AO3,GCAO2,GCAO3)
  ELSE
    CALL II_get_mixed_overlap(lupri,luerr,setting,S23,AO2,AO3,GCAO2,GCAO3)
  ENDIF
ENDIF
END SUBROUTINE get_S23

SUBROUTINE get_R33(R33,AO2,AO3,GC3,nbast2,nbast3,constrain_factor,setting,lupri,luerr)
implicit none
TYPE(matrix),intent(INOUT)    :: R33
TYPE(lssetting),intent(inout) :: setting
Integer,intent(IN)            :: AO2,AO3,nbast2,nbast3,lupri,luerr
Logical,intent(IN)            :: GC3
real(realk),intent(IN)        :: constrain_factor
!
Character(80) :: Filename
TYPE(matrix)  :: S23,T23

write(Filename,'(A8,L1)') 'ADMM_R33',GC3
IF (io_file_exist(Filename,setting%IO)) THEN
  call io_read_mat(R33,Filename,setting%IO,LUPRI,LUERR)
ELSE
  CALL mat_init(S23,nbast2,nbast3)
  CALL mat_init(T23,nbast2,nbast3)
  
  CALL get_S23(S23,AO2,AO3,.FALSE.,GC3,setting,lupri,luerr)
  CALL get_T23(setting,lupri,luerr,T23,nbast2,nbast3,AO2,AO3,.FALSE.,GC3,1E0_realk)
  CALL mat_mul(S23,T23,'t','n',1E0_realk,0E0_realk,R33)

  CALL mat_free(S23)
  CALL mat_free(T23)
ENDIF

call mat_scal(constrain_factor*constrain_factor,R33)

END SUBROUTINE get_R33


SUBROUTINE transform_D3_to_D2(D,D2,setting,lupri,luerr,n2,n3,AO2,AO3,&
                              & McWeeny,GCAO2,GCAO3,constrain_factor)
  implicit none
  type(matrix),intent(in)    :: D     !level 3 matrix input 
  type(matrix),intent(inout) :: D2 !level 2 matrix input 
  type(lssetting) :: setting
  integer :: n2,n3,AO3,AO2,lupri,luerr
  logical :: McWeeny,GCAO2,GCAO3
  real(realk),intent(IN)        :: constrain_factor
  !
  TYPE(MATRIX)       :: S22,tmp23,T23
  Logical            :: purify_failed

  CALL mat_init(T23,n2,n3)
  CALL mat_init(tmp23,n2,n3)

  CALL get_T23(setting,lupri,luerr,T23,n2,n3,&
               & AO2,AO3,GCAO2,GCAO3,constrain_factor)

  CALL mat_mul(T23,D,'n','n',1E0_realk,0E0_realk,tmp23)
  CALL mat_mul(tmp23,T23,'n','t',1E0_realk,0E0_realk,D2)
 
  IF (McWeeny) THEN
    CALL mat_init(S22,n2,n2)
    CALL get_S22(S22,AO2,GCAO2,setting,lupri,luerr)
    CALL McWeeney_purify(S22,D2,purify_failed)
    IF (purify_failed) THEN
      write(lupri,'(1X,A)') 'McWeeny purification failed for ADMM D2 matrix- reverting to the non-purified D2'
    ELSE
      write(lupri,'(1X,A)') 'McWeeny purified ADMM D2 matrix'
    ENDIF
    CALL mat_free(S22)
  ENDIF
   CALL mat_free(T23)
   CALL mat_free(tmp23)
END SUBROUTINE TRANSFORM_D3_TO_D2

! LAMBDA_Q = (2/N) Tr([k2-xc2] d')
! LAMBDA_S = (2/N) [ Tr(k2 d') - xi**(2/3) Tr(x2 d')] - 2/(3N) xi**(2/3) Ex2(d')
! LAMBDA_P = (2/N) xi**2 Tr([k2-xc2] d)
SUBROUTINE get_large_Lambda(Lambda,k2,x2,D2,E_x2,chi,setting)
implicit none
real(realk),intent(out)    :: Lambda
type(matrix),intent(in)    :: k2          !Exact exchange matrix in the small basis
type(matrix),intent(in)    :: x2          !GGA exchange matrix in the small basis
type(matrix),intent(in)    :: D2          !Projected density matrix
real(realk),intent(in)     :: E_x2        !GGA exchange energy in the small basis
real(realk),intent(in)     :: chi !the small lambda
type(lssetting),intent(in) :: setting
!
integer     :: NbEl
logical     :: unrest
real(realk) :: fac
!
unrest = matrix_type.EQ.mtype_unres_dense
fac = 1E0_realk
IF (unrest) fac=2E0_realk
!
NbEl = setting%molecule(1)%p%nelectrons
IF (setting%scheme%ADMMS) THEN
   Lambda = 2E0_realk/NbEl * ( mat_trAB(k2,D2)/fac-mat_trAB(x2,D2)/fac - E_x2/3.E0_realk )
ELSEIF (setting%scheme%ADMMP) THEN
   Lambda = 2E0_realk/NbEl * chi**(4.E0_realk) * (mat_trAB(k2,D2)/fac - E_x2)
ELSEIF (setting%scheme%ADMMQ) THEN
   Lambda = 2E0_realk/NbEl*(mat_trAB(k2,D2)-mat_trAB(x2,D2))/fac
ELSE
   Lambda = 0E0_realk
ENDIF
END SUBROUTINE get_large_Lambda
!****ADMM specific routines ends here

End MODULE IntegralInterfaceMOD

