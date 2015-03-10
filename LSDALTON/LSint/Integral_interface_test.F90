SUBROUTINE II_test_uncontAObatch(LUPRI,LUERR,SETTING)
  use precision
  use TYPEDEFTYPE  
  use TYPEDEF
  use Matrix_module
  use Matrix_Operations, only: mat_init,mat_add,mat_free,mat_trab
  use ls_Integral_Interface
  use LSparameters
IMPLICIT NONE
TYPE(LSSETTING)       :: SETTING
INTEGER               :: LUPRI,LUERR
!
Integer             :: nbast
TYPE(MATRIX)        :: S1,S2,tmp
logical :: OLD_UNCONT
SETTING%SCHEME%intTHRESHOLD=SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
nbast = getNbasis(AOdfAux,ContractedintType,SETTING%MOLECULE(1)%p,LUPRI)
call mat_init(S1,nbast,nbast)
call mat_init(S2,nbast,nbast)
OLD_UNCONT = SETTING%SCHEME%UNCONT
SETTING%SCHEME%UNCONT = .TRUE.
call initIntegralOutputDims(setting%output,nbast,nbast,1,1,1)
CALL ls_getIntegrals(AORegular,AORegular,AOempty,AOempty,&
     &OverlapOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,S1,.FALSE.)
SETTING%SCHEME%UNCONT = OLD_UNCONT 
call initIntegralOutputDims(setting%output,nbast,nbast,1,1,1)
CALL ls_getIntegrals(AOdfAux,AOdfAux,AOempty,AOempty,&
     &OverlapOperator,RegularSpec,ContractedInttype,SETTING,LUPRI,LUERR)
CALL retrieve_Output(lupri,setting,S2,.FALSE.)
call mat_init(tmp,nbast,nbast)
call mat_add(1E0_realk,S1,-1E0_realk,S2,tmp)
IF(ABS(mat_trab(tmp,tmp)).GT.1.0E-10_realk)THEN
   call lsquit('Error in build_AO with uncont param',-1)
ELSE
   WRITE(lupri,'(A)') 'II_test_uncontAObatch successful'
ENDIF
call mat_free(S1)
call mat_free(S2)
call mat_free(tmp)
END SUBROUTINE II_test_uncontAObatch


