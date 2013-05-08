MODULE density_fitting
  USE mm_global_paras_mod
  IMPLICIT NONE
  PUBLIC :: init_fit_param 
  ! Included for density fitting
  TYPE densfitpar_t
     LOGICAL :: AUX, LHS_AUX, RHS_AUX
     INTEGER :: NAUXMOM
  END TYPE densfitpar_t
!  type(densfitpar_t), SAVE :: fit
  type(densfitpar_t):: fit
CONTAINS
   SUBROUTINE init_fit_param
      IMPLICIT NONE
      fit%AUX     = .FALSE.
      fit%LHS_AUX = .FALSE.
      fit%RHS_AUX = .FALSE.
      fit%NAUXMOM = 0
    END SUBROUTINE init_fit_param
END MODULE density_fitting
