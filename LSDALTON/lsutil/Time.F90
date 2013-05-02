!> @file 
!> timing module
MODULE LSTIMING
use precision
!use typedef
save
Integer :: LSIUNIT = 6
logical :: LEVEL1ACTIVE,LEVEL2ACTIVE
Integer :: matoptlevel
LOGICAL :: LSTIME_PRINT = .FALSE.
!job specifiers for matrix operations 
Integer,parameter :: JOB_mat_set_from_full = 1
Integer,parameter :: JOB_mat_to_full = 2
Integer,parameter :: JOB_mat_trans = 3
Integer,parameter :: JOB_mat_assign = 4
Integer,parameter :: JOB_mat_copy = 5
Integer,parameter :: JOB_mat_tr = 6
Integer,parameter :: JOB_mat_trAB = 7 
Integer,parameter :: JOB_mat_mul = 8
Integer,parameter :: JOB_mat_add = 9
Integer,parameter :: JOB_mat_daxpy = 10
Integer,parameter :: JOB_mat_dotproduct = 11
Integer,parameter :: JOB_mat_sqnorm2 = 12
Integer,parameter :: JOB_mat_diag_f = 13
Integer,parameter :: JOB_mat_create_block = 14
Integer,parameter :: JOB_mat_add_block = 15
Integer,parameter :: JOB_mat_retrieve_block = 16
Integer,parameter :: JOB_mat_scal = 17
Integer,parameter :: JOB_mat_zero = 18
Integer,parameter :: JOB_mat_write_to_disk = 19
Integer,parameter :: JOB_mat_read_from_disk = 20
!job specifiers for integral interface routines. 
Integer,parameter :: JOB_II_get_overlap = 21
Integer,parameter :: JOB_II_get_magderivOverlap = 22
Integer,parameter :: JOB_II_get_maggradOverlap = 23
Integer,parameter :: JOB_II_get_kinetic = 24
Integer,parameter :: JOB_II_get_nucel_mat = 25
Integer,parameter :: JOB_II_GET_K_GRADIENT = 26
Integer,parameter :: JOB_II_get_J_gradient = 27
Integer,parameter :: JOB_II_get_ne_gradient = 28
Integer,parameter :: JOB_II_GET_KINETIC_GRADIENT = 29
Integer,parameter :: JOB_II_get_nucpot = 30
Integer,parameter :: JOB_II_GET_NN_GRADIENT = 31
Integer,parameter :: JOB_II_get_prop = 32
Integer,parameter :: JOB_II_get_prop_expval = 33
Integer,parameter :: JOB_II_get_coulomb_mat = 34
Integer,parameter :: JOB_II_GET_EXCHANGE_MAT = 35
Integer,parameter :: JOB_II_get_coulomb_and_exchange_mat = 36
Integer,parameter :: JOB_II_precalc_ScreenMat = 38
Integer,parameter :: JOB_II_GET_REORTHONORMALIZATION = 39
Integer,parameter :: JOB_II_GET_REORTHONORMALIZATION2 = 40
Integer,parameter :: JOB_II_get_geoderivOverlap = 41
Integer,parameter :: JOB_II_get_geoderivExchange = 42
Integer,parameter :: JOB_II_get_magderiv_4center_eri = 43
Integer,parameter :: JOB_II_get_magderivK_4center_eri = 44
Integer,parameter :: JOB_II_GET_MAGDERIVJ_4CENTER_ERI = 45
Integer,parameter :: JOB_II_GET_DECPACKED4CENTER_J_ERI = 46
!XC
Integer,parameter :: JOB_II_get_xc_Fock_mat = 47
Integer,parameter :: JOB_II_get_xc_geoderiv_molgrad = 48
Integer,parameter :: JOB_II_get_xc_linrsp = 49
Integer,parameter :: JOB_II_get_xc_quadrsp = 50
Integer,parameter :: JOB_II_get_xc_magderiv_kohnsham_mat = 51
Integer,parameter :: JOB_II_get_xc_magderiv_linrsp = 52
Integer,parameter :: JOB_II_get_xc_geoderiv_FxDgrad = 53
Integer,parameter :: JOB_II_get_xc_geoderiv_GxDgrad = 54
Integer,parameter :: JOB_homolumo = 55
Integer,parameter :: JOB_mat_density_from_orbs = 56
Integer,parameter :: JOB_mat_chol = 57
Integer,parameter :: JOB_mat_inv = 58
Integer,parameter :: JOB_mat_dsyev = 59
Integer,parameter :: JOB_mat_setlowertriangular_zero = 60
Integer,parameter :: JOB_mat_dsyevx = 61
Integer,parameter :: JOB_mat_identity = 62
Integer,parameter :: JOB_mat_scal_dia = 63
Integer,parameter :: JOB_mat_scal_dia_vec = 64
Integer,parameter :: JOB_II_GET_DECPACKED4CENTER_K_ERI = 65
Integer,parameter :: JOB_II_get_geoderivCoulomb = 66
!counter for matrix operations
Integer :: nmat_set_from_full(3) 
Integer :: nmat_to_full(3)
Integer :: nmat_trans(3) 
Integer :: nmat_assign(3)
Integer :: nmat_copy(3) 
Integer :: nmat_tr(3) 
Integer :: nmat_trAB(3)
Integer :: nmat_mul(3) 
Integer :: nmat_add(3) 
Integer :: nmat_daxpy(3) 
Integer :: nmat_dotproduct(3)
Integer :: nmat_sqnorm2(3) 
Integer :: nmat_diag_f(3) 
Integer :: nmat_create_block(3)
Integer :: nmat_add_block(3) 
Integer :: nmat_retrieve_block(3) 
Integer :: nmat_scal(3) 
Integer :: nmat_zero(3) 
Integer :: nmat_density_from_orbs(3) 
Integer :: nmat_scal_dia_vec(3) 
Integer :: nmat_scal_dia(3) 
Integer :: nmat_chol(3) 
Integer :: nmat_inv(3)
Integer :: nmat_dsyev(3)
Integer :: nmat_dsyevx(3)
Integer :: nmat_identity(3)
Integer :: nmat_setlowertriangular_zero(3)
Integer :: nmat_write_to_disk(3) 
Integer :: nmat_read_from_disk(3)
!counters for integral interface routines. 
Integer :: nII_get_overlap(3)
Integer :: nII_get_magderivOverlap(3)
Integer :: nII_get_maggradOverlap(3)
Integer :: nII_get_kinetic(3)
Integer :: nII_get_nucel_mat(3)
Integer :: nII_GET_NN_GRADIENT(3)
Integer :: nII_precalc_ScreenMat(3)
Integer :: nII_get_ne_gradient(3)
Integer :: nII_get_geoderivOverlap(3)
Integer :: nII_get_geoderivExchange(3)
Integer :: nII_get_geoderivCoulomb(3)
Integer :: nII_GET_KINETIC_GRADIENT(3)
Integer :: nII_GET_REORTHONORMALIZATION(3)
Integer :: nII_GET_REORTHONORMALIZATION2(3)
Integer :: nII_get_J_gradient(3)
Integer :: nII_GET_K_GRADIENT(3)
Integer :: nII_get_nucpot(3)
Integer :: nII_get_prop(3)
Integer :: nII_get_prop_expval(3)
Integer :: nII_GET_EXCHANGE_MAT(3)
Integer :: nII_get_coulomb_mat(3)
Integer :: nII_get_coulomb_and_exchange_mat(3)
Integer :: nII_get_magderiv_4center_eri(3)
Integer :: nII_get_magderivK_4center_eri(3)
Integer :: nII_GET_MAGDERIVJ_4CENTER_ERI(3)
Integer :: nII_GET_DECPACKED4CENTER_J_ERI(3)
Integer :: nII_GET_DECPACKED4CENTER_K_ERI(3)
!counters for XC
Integer :: nII_get_xc_Fock_mat(3)
Integer :: nII_get_xc_geoderiv_molgrad(3)
Integer :: nII_get_xc_linrsp(3)
Integer :: nII_get_xc_quadrsp(3)
Integer :: nII_get_xc_magderiv_kohnsham_mat(3)
Integer :: nII_get_xc_magderiv_linrsp(3)
Integer :: nII_get_xc_geoderiv_FxDgrad(3)
Integer :: nII_get_xc_geoderiv_GxDgrad(3)
Integer :: nhomolumo(3)

real(realk) :: CPUTIME_TOTAL_mat(3),CPUTIME_mat_set_from_full(3),CPUTIME_mat_to_full(3)
real(realk) :: CPUTIME_mat_trans(3),CPUTIME_mat_assign(3),CPUTIME_mat_copy(3),CPUTIME_mat_tr(3)
real(realk) :: CPUTIME_mat_trAB(3),CPUTIME_mat_mul(3),CPUTIME_mat_add(3),CPUTIME_mat_daxpy(3)
real(realk) :: CPUTIME_mat_dotproduct(3),CPUTIME_mat_sqnorm2(3),CPUTIME_mat_diag_f(3)
real(realk) :: CPUTIME_mat_create_block(3),CPUTIME_mat_add_block(3),CPUTIME_mat_retrieve_block(3)
real(realk) :: CPUTIME_mat_scal(3),CPUTIME_mat_zero(3),CPUTIME_mat_write_to_disk(3)
real(realk) :: CPUTIME_mat_read_from_disk(3),CPUTIME_mat_density_from_orbs(3)
real(realk) :: CPUTIME_mat_scal_dia(3),CPUTIME_mat_scal_dia_vec(3)
real(realk) :: CPUTIME_mat_chol(3),CPUTIME_mat_inv(3),CPUTIME_mat_dsyev(3)
real(realk) :: CPUTIME_mat_setlowertriangular_zero(3),CPUTIME_mat_dsyevx(3),CPUTIME_mat_identity(3)
real(realk) :: WALLTIME_TOTAL_mat(3),WALLTIME_mat_set_from_full(3),WALLTIME_mat_to_full(3)
real(realk) :: WALLTIME_mat_trans(3),WALLTIME_mat_assign(3),WALLTIME_mat_copy(3),WALLTIME_mat_tr(3)
real(realk) :: WALLTIME_mat_trAB(3),WALLTIME_mat_mul(3),WALLTIME_mat_add(3),WALLTIME_mat_daxpy(3)
real(realk) :: WALLTIME_mat_dotproduct(3),WALLTIME_mat_sqnorm2(3),WALLTIME_mat_diag_f(3)
real(realk) :: WALLTIME_mat_create_block(3),WALLTIME_mat_add_block(3),WALLTIME_mat_retrieve_block(3)
real(realk) :: WALLTIME_mat_scal(3),WALLTIME_mat_zero(3),WALLTIME_mat_write_to_disk(3)
real(realk) :: WALLTIME_mat_read_from_disk(3),WALLTIME_mat_density_from_orbs(3)
real(realk) :: WALLTIME_mat_scal_dia(3),WALLTIME_mat_scal_dia_vec(3)
real(realk) :: WALLTIME_mat_chol(3),WALLTIME_mat_inv(3),WALLTIME_mat_dsyev(3)
real(realk) :: WALLTIME_mat_setlowertriangular_zero(3),WALLTIME_mat_dsyevx(3),WALLTIME_mat_identity(3)
!Timers for Integral interface routines
real(realk) :: CPUTIME_II_get_overlap(3),CPUTIME_II_get_magderivOverlap(3)
real(realk) :: CPUTIME_II_get_maggradOverlap(3),CPUTIME_II_get_kinetic(3),CPUTIME_II_get_nucel_mat(3)
real(realk) :: CPUTIME_II_GET_NN_GRADIENT(3),CPUTIME_II_precalc_ScreenMat(3)
real(realk) :: CPUTIME_II_get_ne_gradient(3),CPUTIME_II_get_geoderivOverlap(3),CPUTIME_II_get_geoderivExchange(3)
real(realk) :: CPUTIME_II_get_geoderivCoulomb(3)
real(realk) :: CPUTIME_II_GET_KINETIC_GRADIENT(3),CPUTIME_II_GET_REORTHONORMALIZATION(3),CPUTIME_II_GET_REORTHONORMALIZATION2(3)
real(realk) :: CPUTIME_II_get_J_gradient(3),CPUTIME_II_GET_K_GRADIENT(3)
real(realk) :: CPUTIME_II_get_nucpot(3),CPUTIME_II_get_prop(3),CPUTIME_II_get_prop_expval(3)
real(realk) :: CPUTIME_II_GET_EXCHANGE_MAT(3),CPUTIME_II_get_coulomb_mat(3)
real(realk) :: CPUTIME_II_get_coulomb_and_exchange_mat(3)
real(realk) :: CPUTIME_II_get_magderiv_4center_eri(3),CPUTIME_II_get_magderivK_4center_eri(3)
real(realk) :: CPUTIME_II_GET_MAGDERIVJ_4CENTER_ERI(3),CPUTIME_II_GET_DECPACKED4CENTER_J_ERI(3)
real(realk) :: CPUTIME_II_GET_DECPACKED4CENTER_K_ERI(3)
!Timers for XC
real(realk) :: CPUTIME_II_get_xc_Fock_mat(3),CPUTIME_II_get_xc_geoderiv_molgrad(3)
real(realk) :: CPUTIME_II_get_xc_linrsp(3),CPUTIME_II_get_xc_quadrsp(3)
real(realk) :: CPUTIME_II_get_xc_magderiv_kohnsham_mat(3),CPUTIME_II_get_xc_magderiv_linrsp(3)
real(realk) :: CPUTIME_II_get_xc_geoderiv_FxDgrad(3),CPUTIME_II_get_xc_geoderiv_GxDgrad(3)
real(realk) :: CPUTIME_homolumo(3)
!Timers for Integral interface routines
real(realk) :: WALLTIME_II_get_overlap(3),WALLTIME_II_get_magderivOverlap(3)
real(realk) :: WALLTIME_II_get_maggradOverlap(3),WALLTIME_II_get_kinetic(3),WALLTIME_II_get_nucel_mat(3)
real(realk) :: WALLTIME_II_GET_NN_GRADIENT(3),WALLTIME_II_precalc_ScreenMat(3)
real(realk) :: WALLTIME_II_get_ne_gradient(3),WALLTIME_II_get_geoderivOverlap(3),WALLTIME_II_get_geoderivExchange(3)
real(realk) :: WALLTIME_II_get_geoderivCoulomb(3)
real(realk) :: WALLTIME_II_GET_KINETIC_GRADIENT(3),WALLTIME_II_GET_REORTHONORMALIZATION(3),WALLTIME_II_GET_REORTHONORMALIZATION2(3)
real(realk) :: WALLTIME_II_get_J_gradient(3),WALLTIME_II_GET_K_GRADIENT(3)
real(realk) :: WALLTIME_II_get_nucpot(3),WALLTIME_II_get_prop(3),WALLTIME_II_get_prop_expval(3)
real(realk) :: WALLTIME_II_GET_EXCHANGE_MAT(3),WALLTIME_II_get_coulomb_mat(3)
real(realk) :: WALLTIME_II_get_coulomb_and_exchange_mat(3)
real(realk) :: WALLTIME_II_get_magderiv_4center_eri(3),WALLTIME_II_get_magderivK_4center_eri(3)
real(realk) :: WALLTIME_II_GET_MAGDERIVJ_4CENTER_ERI(3),WALLTIME_II_GET_DECPACKED4CENTER_J_ERI(3)
real(realk) :: WALLTIME_II_GET_DECPACKED4CENTER_K_ERI(3)
!Timers for XC
real(realk) :: WALLTIME_II_get_xc_Fock_mat(3),WALLTIME_II_get_xc_geoderiv_molgrad(3)
real(realk) :: WALLTIME_II_get_xc_linrsp(3),WALLTIME_II_get_xc_quadrsp(3)
real(realk) :: WALLTIME_II_get_xc_magderiv_kohnsham_mat(3),WALLTIME_II_get_xc_magderiv_linrsp(3)
real(realk) :: WALLTIME_II_get_xc_geoderiv_FxDgrad(3),WALLTIME_II_get_xc_geoderiv_GxDgrad(3)
real(realk) :: WALLTIME_homolumo(3)
real(realk) :: CPUTIME_TOTAL_II(3),WALLTIME_TOTAL_II(3)
real(realk) :: matcputime1,matwalltime1,matcputime2,matwalltime2
real(realk) :: IIcputime1,IIwalltime1,IIcputime2,IIwalltime2
!$OMP THREADPRIVATE(cputime1,walltime1,cputime2,walltime2)
contains
subroutine set_matop_timer_optlevel(optl)
implicit none
integer :: optl
matoptlevel = optl
IF(optl.EQ.1)THEN
   LEVEL1ACTIVE = .TRUE.
ENDIF
IF(optl.EQ.2)THEN
   LEVEL2ACTIVE = .TRUE.
ENDIF
end subroutine set_matop_timer_optlevel

subroutine init_timers()
implicit none
LEVEL1ACTIVE = .FALSE.
LEVEL2ACTIVE = .FALSE.
matoptlevel = 3
!mat
CPUTIME_TOTAL_mat(1:3) = 0.0E0_realk
CPUTIME_mat_set_from_full(1:3) = 0.0E0_realk
CPUTIME_mat_to_full(1:3) = 0.0E0_realk
CPUTIME_mat_trans(1:3) = 0.0E0_realk
CPUTIME_mat_assign(1:3) = 0.0E0_realk
CPUTIME_mat_copy(1:3) = 0.0E0_realk
CPUTIME_mat_tr(1:3) = 0.0E0_realk
CPUTIME_mat_trAB(1:3) = 0.0E0_realk
CPUTIME_mat_mul(1:3) = 0.0E0_realk
CPUTIME_mat_add(1:3) = 0.0E0_realk
CPUTIME_mat_daxpy(1:3) = 0.0E0_realk
CPUTIME_mat_dotproduct(1:3) = 0.0E0_realk
CPUTIME_mat_sqnorm2(1:3) = 0.0E0_realk
CPUTIME_mat_diag_f(1:3) = 0.0E0_realk
CPUTIME_mat_create_block(1:3) = 0.0E0_realk
CPUTIME_mat_add_block(1:3) = 0.0E0_realk
CPUTIME_mat_retrieve_block(1:3) = 0.0E0_realk
CPUTIME_mat_scal(1:3) = 0.0E0_realk
CPUTIME_mat_zero(1:3) = 0.0E0_realk
CPUTIME_mat_density_from_orbs(1:3) = 0.0E0_realk
CPUTIME_mat_scal_dia_vec(1:3) = 0.0E0_realk
CPUTIME_mat_scal_dia(1:3) = 0.0E0_realk
CPUTIME_mat_chol(1:3) = 0.0E0_realk
CPUTIME_mat_inv(1:3) = 0.0E0_realk
CPUTIME_mat_dsyev(1:3) = 0.0E0_realk
CPUTIME_mat_dsyevx(1:3) = 0.0E0_realk
CPUTIME_mat_identity(1:3) = 0.0E0_realk
CPUTIME_mat_setlowertriangular_zero(1:3) = 0.0E0_realk
CPUTIME_mat_write_to_disk(1:3) = 0.0E0_realk
CPUTIME_mat_read_from_disk(1:3) = 0.0E0_realk

CPUTIME_II_get_overlap(1:3) = 0.0E0_realk  
CPUTIME_II_get_magderivOverlap(1:3) = 0.0E0_realk  
CPUTIME_II_get_maggradOverlap(1:3) = 0.0E0_realk  
CPUTIME_II_get_kinetic(1:3) = 0.0E0_realk  
CPUTIME_II_get_nucel_mat(1:3) = 0.0E0_realk  
CPUTIME_II_GET_NN_GRADIENT(1:3) = 0.0E0_realk  
CPUTIME_II_precalc_ScreenMat(1:3) = 0.0E0_realk  
CPUTIME_II_get_ne_gradient(1:3) = 0.0E0_realk  
CPUTIME_II_get_geoderivOverlap(1:3) = 0.0E0_realk  
CPUTIME_II_get_geoderivExchange(1:3) = 0.0E0_realk  
CPUTIME_II_get_geoderivCoulomb(1:3) = 0.0E0_realk  
CPUTIME_II_GET_KINETIC_GRADIENT(1:3) = 0.0E0_realk  
CPUTIME_II_GET_REORTHONORMALIZATION(1:3) = 0.0E0_realk  
CPUTIME_II_GET_REORTHONORMALIZATION2(1:3) = 0.0E0_realk  
CPUTIME_II_get_J_gradient(1:3) = 0.0E0_realk  
CPUTIME_II_GET_K_GRADIENT(1:3) = 0.0E0_realk  
CPUTIME_II_get_nucpot(1:3) = 0.0E0_realk  
CPUTIME_II_get_prop(1:3) = 0.0E0_realk  
CPUTIME_II_get_prop_expval(1:3) = 0.0E0_realk  
CPUTIME_II_GET_EXCHANGE_MAT(1:3) = 0.0E0_realk  
CPUTIME_II_get_coulomb_mat(1:3) = 0.0E0_realk  
CPUTIME_II_get_coulomb_and_exchange_mat(1:3) = 0.0E0_realk  
CPUTIME_II_get_magderiv_4center_eri(1:3) = 0.0E0_realk  
CPUTIME_II_get_magderivK_4center_eri(1:3) = 0.0E0_realk  
CPUTIME_II_GET_MAGDERIVJ_4CENTER_ERI(1:3) = 0.0E0_realk  
CPUTIME_II_GET_DECPACKED4CENTER_J_ERI(1:3) = 0.0E0_realk  
CPUTIME_II_GET_DECPACKED4CENTER_K_ERI(1:3) = 0.0E0_realk  
CPUTIME_II_get_xc_Fock_mat(1:3) = 0.0E0_realk  
CPUTIME_II_get_xc_geoderiv_molgrad(1:3) = 0.0E0_realk  
CPUTIME_II_get_xc_linrsp(1:3) = 0.0E0_realk  
CPUTIME_II_get_xc_quadrsp(1:3) = 0.0E0_realk  
CPUTIME_II_get_xc_magderiv_kohnsham_mat(1:3) = 0.0E0_realk  
CPUTIME_II_get_xc_magderiv_linrsp(1:3) = 0.0E0_realk  
CPUTIME_II_get_xc_geoderiv_FxDgrad(1:3) = 0.0E0_realk  
CPUTIME_II_get_xc_geoderiv_GxDgrad(1:3) = 0.0E0_realk  
CPUTIME_homolumo(1:3) = 0.0E0_realk  

WALLTIME_TOTAL_mat(1:3) = 0.0E0_realk
WALLTIME_mat_set_from_full(1:3) = 0.0E0_realk
WALLTIME_mat_to_full(1:3) = 0.0E0_realk
WALLTIME_mat_trans(1:3) = 0.0E0_realk
WALLTIME_mat_assign(1:3) = 0.0E0_realk
WALLTIME_mat_copy(1:3) = 0.0E0_realk
WALLTIME_mat_tr(1:3) = 0.0E0_realk
WALLTIME_mat_trAB(1:3) = 0.0E0_realk
WALLTIME_mat_mul(1:3) = 0.0E0_realk
WALLTIME_mat_add(1:3) = 0.0E0_realk
WALLTIME_mat_daxpy(1:3) = 0.0E0_realk
WALLTIME_mat_dotproduct(1:3) = 0.0E0_realk
WALLTIME_mat_sqnorm2(1:3) = 0.0E0_realk
WALLTIME_mat_diag_f(1:3) = 0.0E0_realk
WALLTIME_mat_create_block(1:3) = 0.0E0_realk
WALLTIME_mat_add_block(1:3) = 0.0E0_realk
WALLTIME_mat_retrieve_block(1:3) = 0.0E0_realk
WALLTIME_mat_scal(1:3) = 0.0E0_realk
WALLTIME_mat_zero(1:3) = 0.0E0_realk
WALLTIME_mat_density_from_orbs(1:3) = 0.0E0_realk
WALLTIME_mat_scal_dia_vec(1:3) = 0.0E0_realk
WALLTIME_mat_scal_dia(1:3) = 0.0E0_realk
WALLTIME_mat_chol(1:3) = 0.0E0_realk
WALLTIME_mat_inv(1:3) = 0.0E0_realk
WALLTIME_mat_dsyev(1:3) = 0.0E0_realk
WALLTIME_mat_dsyevx(1:3) = 0.0E0_realk
WALLTIME_mat_identity(1:3) = 0.0E0_realk
WALLTIME_mat_setlowertriangular_zero(1:3) = 0.0E0_realk
WALLTIME_mat_write_to_disk(1:3) = 0.0E0_realk
WALLTIME_mat_read_from_disk(1:3) = 0.0E0_realk

WALLTIME_II_get_overlap(1:3) = 0.0E0_realk  
WALLTIME_II_get_magderivOverlap(1:3) = 0.0E0_realk  
WALLTIME_II_get_maggradOverlap(1:3) = 0.0E0_realk  
WALLTIME_II_get_kinetic(1:3) = 0.0E0_realk  
WALLTIME_II_get_nucel_mat(1:3) = 0.0E0_realk  
WALLTIME_II_GET_NN_GRADIENT(1:3) = 0.0E0_realk  
WALLTIME_II_precalc_ScreenMat(1:3) = 0.0E0_realk  
WALLTIME_II_get_ne_gradient(1:3) = 0.0E0_realk  
WALLTIME_II_get_geoderivOverlap(1:3) = 0.0E0_realk  
WALLTIME_II_get_geoderivExchange(1:3) = 0.0E0_realk  
WALLTIME_II_get_geoderivCoulomb(1:3) = 0.0E0_realk  
WALLTIME_II_GET_KINETIC_GRADIENT(1:3) = 0.0E0_realk  
WALLTIME_II_GET_REORTHONORMALIZATION(1:3) = 0.0E0_realk  
WALLTIME_II_GET_REORTHONORMALIZATION2(1:3) = 0.0E0_realk  
WALLTIME_II_get_J_gradient(1:3) = 0.0E0_realk  
WALLTIME_II_GET_K_GRADIENT(1:3) = 0.0E0_realk  
WALLTIME_II_get_nucpot(1:3) = 0.0E0_realk  
WALLTIME_II_get_prop(1:3) = 0.0E0_realk  
WALLTIME_II_get_prop_expval(1:3) = 0.0E0_realk  
WALLTIME_II_GET_EXCHANGE_MAT(1:3) = 0.0E0_realk  
WALLTIME_II_get_coulomb_mat(1:3) = 0.0E0_realk  
WALLTIME_II_get_coulomb_and_exchange_mat(1:3) = 0.0E0_realk  
WALLTIME_II_get_magderiv_4center_eri(1:3) = 0.0E0_realk  
WALLTIME_II_get_magderivK_4center_eri(1:3) = 0.0E0_realk  
WALLTIME_II_GET_MAGDERIVJ_4CENTER_ERI(1:3) = 0.0E0_realk  
WALLTIME_II_GET_DECPACKED4CENTER_J_ERI(1:3) = 0.0E0_realk  
WALLTIME_II_GET_DECPACKED4CENTER_K_ERI(1:3) = 0.0E0_realk  
WALLTIME_II_get_xc_Fock_mat(1:3) = 0.0E0_realk  
WALLTIME_II_get_xc_geoderiv_molgrad(1:3) = 0.0E0_realk  
WALLTIME_II_get_xc_linrsp(1:3) = 0.0E0_realk  
WALLTIME_II_get_xc_quadrsp(1:3) = 0.0E0_realk  
WALLTIME_II_get_xc_magderiv_kohnsham_mat(1:3) = 0.0E0_realk  
WALLTIME_II_get_xc_magderiv_linrsp(1:3) = 0.0E0_realk  
WALLTIME_II_get_xc_geoderiv_FxDgrad(1:3) = 0.0E0_realk  
WALLTIME_II_get_xc_geoderiv_GxDgrad(1:3) = 0.0E0_realk  
WALLTIME_homolumo(1:3) = 0.0E0_realk  

nmat_set_from_full = 0
nmat_to_full = 0
nmat_trans = 0
nmat_assign = 0
nmat_copy = 0
nmat_tr = 0
nmat_trAB = 0
nmat_mul = 0
nmat_add = 0
nmat_daxpy = 0
nmat_dotproduct = 0
nmat_sqnorm2 = 0
nmat_diag_f = 0
nmat_create_block = 0
nmat_add_block = 0
nmat_retrieve_block = 0
nmat_scal = 0
nmat_zero = 0
nmat_density_from_orbs = 0
nmat_scal_dia_vec = 0
nmat_scal_dia = 0
nmat_chol = 0
nmat_inv = 0
nmat_dsyev = 0
nmat_dsyevx = 0
nmat_identity = 0
nmat_setlowertriangular_zero = 0
nmat_write_to_disk = 0
nmat_read_from_disk = 0
!
nII_get_overlap(1:3) = 0  
nII_get_magderivOverlap(1:3) = 0  
nII_get_maggradOverlap(1:3) = 0  
nII_get_kinetic(1:3) = 0  
nII_get_nucel_mat(1:3) = 0  
nII_GET_NN_GRADIENT(1:3) = 0  
nII_precalc_ScreenMat(1:3) = 0  
nII_get_ne_gradient(1:3) = 0  
nII_get_geoderivOverlap(1:3) = 0  
nII_get_geoderivExchange(1:3) = 0  
nII_get_geoderivCoulomb(1:3) = 0  
nII_GET_KINETIC_GRADIENT(1:3) = 0  
nII_GET_REORTHONORMALIZATION(1:3) = 0  
nII_GET_REORTHONORMALIZATION2(1:3) = 0  
nII_get_J_gradient(1:3) = 0  
nII_GET_K_GRADIENT(1:3) = 0  
nII_get_nucpot(1:3) = 0  
nII_get_prop(1:3) = 0  
nII_get_prop_expval(1:3) = 0  
nII_GET_EXCHANGE_MAT(1:3) = 0  
nII_get_coulomb_mat(1:3) = 0  
nII_get_coulomb_and_exchange_mat(1:3) = 0  
nII_get_magderiv_4center_eri(1:3) = 0  
nII_get_magderivK_4center_eri(1:3) = 0  
nII_GET_MAGDERIVJ_4CENTER_ERI(1:3) = 0  
nII_GET_DECPACKED4CENTER_J_ERI(1:3) = 0  
nII_GET_DECPACKED4CENTER_K_ERI(1:3) = 0  
nII_get_xc_Fock_mat(1:3) = 0  
nII_get_xc_geoderiv_molgrad(1:3) = 0  
nII_get_xc_linrsp(1:3) = 0  
nII_get_xc_quadrsp(1:3) = 0  
nII_get_xc_magderiv_kohnsham_mat(1:3) = 0  
nII_get_xc_magderiv_linrsp(1:3) = 0  
nII_get_xc_geoderiv_FxDgrad(1:3) = 0  
nII_get_xc_geoderiv_GxDgrad(1:3) = 0  
nhomolumo(1:3) = 0  

end subroutine init_timers

subroutine time_mat_operations1()
implicit none
#ifdef VAR_TIME
!$OMP CRITICAL (timematop)
CALL LS_GETTIM(matcputime1,matwalltime1)
!$OMP END CRITICAL (timematop)
#endif

end subroutine time_mat_operations1

subroutine time_mat_operations2(job)
implicit none
integer :: job
real(realk) :: deltacpu,deltawall
#ifdef VAR_TIME
!$OMP CRITICAL (timematop)
CALL LS_GETTIM(matcputime2,matwalltime2)
DeltaCPU =  matcputime2 - matcputime1 !Node CPU time  += task CPU time
DeltaWall = matwalltime2 - matwalltime1 !Node wall time += task wall time
call add_time_to_job(job,deltacpu,deltawall)
!$OMP END CRITICAL (timematop)
#endif
end subroutine time_mat_operations2

subroutine time_II_operations1()
implicit none
#ifdef VAR_TIME
!$OMP CRITICAL (timematop)
CALL LS_GETTIM(IIcputime1,IIwalltime1)
!$OMP END CRITICAL (timematop)
#endif

end subroutine time_II_operations1

subroutine time_II_operations2(job)
implicit none
integer :: job
real(realk) :: deltacpu,deltawall
#ifdef VAR_TIME
!$OMP CRITICAL (timematop)
CALL LS_GETTIM(IIcputime2,IIwalltime2)
DeltaCPU =  IIcputime2 - IIcputime1 !Node CPU time  += task CPU time
DeltaWall = IIwalltime2 - IIwalltime1 !Node wall time += task wall time
call add_time_to_job(job,deltacpu,deltawall)
!$OMP END CRITICAL (timematop)
#endif
end subroutine time_II_operations2

subroutine add_time_to_job(job,deltacpu,deltawall)
implicit none
integer :: job
real(realk) :: deltacpu,deltawall

SELECT CASE(job)
CASE(JOB_mat_set_from_full)
   nmat_set_from_full(matoptlevel) = nmat_set_from_full(matoptlevel)+1 
   CPUTIME_mat_set_from_full(matoptlevel) = CPUTIME_mat_set_from_full(matoptlevel) + deltacpu
   WALLTIME_mat_set_from_full(matoptlevel) = WALLTIME_mat_set_from_full(matoptlevel) + deltawall
CASE(JOB_mat_to_full)
   nmat_to_full(matoptlevel) = nmat_to_full(matoptlevel)+1 
   CPUTIME_mat_to_full(matoptlevel) = CPUTIME_mat_to_full(matoptlevel) + deltacpu
   WALLTIME_mat_to_full(matoptlevel) = WALLTIME_mat_to_full(matoptlevel) + deltawall
CASE(JOB_mat_trans)
   nmat_trans(matoptlevel) = nmat_trans(matoptlevel)+1 
   CPUTIME_mat_trans(matoptlevel) = CPUTIME_mat_trans(matoptlevel) + deltacpu
   WALLTIME_mat_trans(matoptlevel) = WALLTIME_mat_trans(matoptlevel) + deltawall
CASE(JOB_mat_assign)
   nmat_assign(matoptlevel) = nmat_assign(matoptlevel)+1 
   CPUTIME_mat_assign(matoptlevel) = CPUTIME_mat_assign(matoptlevel) + deltacpu
   WALLTIME_mat_assign(matoptlevel) = WALLTIME_mat_assign(matoptlevel) + deltawall
CASE(JOB_mat_copy)
   nmat_copy(matoptlevel) = nmat_copy(matoptlevel)+1 
   CPUTIME_mat_copy(matoptlevel) = CPUTIME_mat_copy(matoptlevel) + deltacpu
   WALLTIME_mat_copy(matoptlevel) = WALLTIME_mat_copy(matoptlevel) + deltawall
CASE(JOB_mat_tr)
   nmat_tr(matoptlevel) = nmat_tr(matoptlevel)+1 
   CPUTIME_mat_tr(matoptlevel) = CPUTIME_mat_tr(matoptlevel) + deltacpu
   WALLTIME_mat_tr(matoptlevel) = WALLTIME_mat_tr(matoptlevel) + deltawall
CASE(JOB_mat_trAB)
   nmat_trAB(matoptlevel) = nmat_trAB(matoptlevel)+1 
   CPUTIME_mat_trAB(matoptlevel) = CPUTIME_mat_trAB(matoptlevel) + deltacpu
   WALLTIME_mat_trAB(matoptlevel) = WALLTIME_mat_trAB(matoptlevel) + deltawall
CASE(JOB_mat_mul)
   nmat_mul(matoptlevel) = nmat_mul(matoptlevel)+1 
   CPUTIME_mat_mul(matoptlevel) = CPUTIME_mat_mul(matoptlevel) + deltacpu
   WALLTIME_mat_mul(matoptlevel) = WALLTIME_mat_mul(matoptlevel) + deltawall
CASE(JOB_mat_add)
   nmat_add(matoptlevel) = nmat_add(matoptlevel)+1 
   CPUTIME_mat_add(matoptlevel) = CPUTIME_mat_add(matoptlevel) + deltacpu
   WALLTIME_mat_add(matoptlevel) = WALLTIME_mat_add(matoptlevel) + deltawall
CASE(JOB_mat_daxpy)
   nmat_daxpy(matoptlevel) = nmat_daxpy(matoptlevel)+1 
   CPUTIME_mat_daxpy(matoptlevel) = CPUTIME_mat_daxpy(matoptlevel) + deltacpu
   WALLTIME_mat_daxpy(matoptlevel) = WALLTIME_mat_daxpy(matoptlevel) + deltawall
CASE(JOB_mat_dotproduct)
   nmat_dotproduct(matoptlevel) = nmat_dotproduct(matoptlevel)+1 
   CPUTIME_mat_dotproduct(matoptlevel) = CPUTIME_mat_dotproduct(matoptlevel) + deltacpu
   WALLTIME_mat_dotproduct(matoptlevel) = WALLTIME_mat_dotproduct(matoptlevel) + deltawall
CASE(JOB_mat_sqnorm2)
   nmat_sqnorm2(matoptlevel) = nmat_sqnorm2(matoptlevel)+1 
   CPUTIME_mat_sqnorm2(matoptlevel) = CPUTIME_mat_sqnorm2(matoptlevel) + deltacpu
   WALLTIME_mat_sqnorm2(matoptlevel) = WALLTIME_mat_sqnorm2(matoptlevel) + deltawall
CASE(JOB_mat_diag_f)
   nmat_diag_f(matoptlevel) = nmat_diag_f(matoptlevel)+1 
   CPUTIME_mat_diag_f(matoptlevel) = CPUTIME_mat_diag_f(matoptlevel) + deltacpu
   WALLTIME_mat_diag_f(matoptlevel) = WALLTIME_mat_diag_f(matoptlevel) + deltawall
CASE(JOB_mat_create_block)
   nmat_create_block(matoptlevel) = nmat_create_block(matoptlevel)+1 
   CPUTIME_mat_create_block(matoptlevel) = CPUTIME_mat_create_block(matoptlevel) + deltacpu
   WALLTIME_mat_create_block(matoptlevel) = WALLTIME_mat_create_block(matoptlevel) + deltawall
CASE(JOB_mat_add_block)
   nmat_add_block(matoptlevel) = nmat_add_block(matoptlevel)+1 
   CPUTIME_mat_add_block(matoptlevel) = CPUTIME_mat_add_block(matoptlevel) + deltacpu
   WALLTIME_mat_add_block(matoptlevel) = WALLTIME_mat_add_block(matoptlevel) + deltawall
CASE(JOB_mat_retrieve_block)
   nmat_retrieve_block(matoptlevel) = nmat_retrieve_block(matoptlevel)+1 
   CPUTIME_mat_retrieve_block(matoptlevel) = CPUTIME_mat_retrieve_block(matoptlevel) + deltacpu
   WALLTIME_mat_retrieve_block(matoptlevel) = WALLTIME_mat_retrieve_block(matoptlevel) + deltawall
CASE(JOB_mat_scal)
   nmat_scal(matoptlevel) = nmat_scal(matoptlevel)+1 
   CPUTIME_mat_scal(matoptlevel) = CPUTIME_mat_scal(matoptlevel) + deltacpu
   WALLTIME_mat_scal(matoptlevel) = WALLTIME_mat_scal(matoptlevel) + deltawall
CASE(JOB_mat_zero)
   nmat_zero(matoptlevel) = nmat_zero(matoptlevel)+1 
   CPUTIME_mat_zero(matoptlevel) = CPUTIME_mat_zero(matoptlevel) + deltacpu
   WALLTIME_mat_zero(matoptlevel) = WALLTIME_mat_zero(matoptlevel) + deltawall
CASE(JOB_mat_density_from_orbs)
   nmat_density_from_orbs(matoptlevel) = nmat_density_from_orbs(matoptlevel)+1 
   CPUTIME_mat_density_from_orbs(matoptlevel) = CPUTIME_mat_density_from_orbs(matoptlevel) + deltacpu
   WALLTIME_mat_density_from_orbs(matoptlevel) = WALLTIME_mat_density_from_orbs(matoptlevel) + deltawall
CASE(JOB_mat_scal_dia_vec)
   nmat_scal_dia_vec(matoptlevel) = nmat_scal_dia_vec(matoptlevel)+1 
   CPUTIME_mat_scal_dia_vec(matoptlevel) = CPUTIME_mat_scal_dia_vec(matoptlevel) + deltacpu
   WALLTIME_mat_scal_dia_vec(matoptlevel) = WALLTIME_mat_scal_dia_vec(matoptlevel) + deltawall
CASE(JOB_mat_scal_dia)
   nmat_scal_dia(matoptlevel) = nmat_scal_dia(matoptlevel)+1 
   CPUTIME_mat_scal_dia(matoptlevel) = CPUTIME_mat_scal_dia(matoptlevel) + deltacpu
   WALLTIME_mat_scal_dia(matoptlevel) = WALLTIME_mat_scal_dia(matoptlevel) + deltawall
CASE(JOB_mat_chol)
   nmat_chol(matoptlevel) = nmat_chol(matoptlevel)+1 
   CPUTIME_mat_chol(matoptlevel) = CPUTIME_mat_chol(matoptlevel) + deltacpu
   WALLTIME_mat_chol(matoptlevel) = WALLTIME_mat_chol(matoptlevel) + deltawall
CASE(JOB_mat_inv)
   nmat_inv(matoptlevel) = nmat_inv(matoptlevel)+1 
   CPUTIME_mat_inv(matoptlevel) = CPUTIME_mat_inv(matoptlevel) + deltacpu
   WALLTIME_mat_inv(matoptlevel) = WALLTIME_mat_inv(matoptlevel) + deltawall
CASE(JOB_mat_dsyev)
   nmat_dsyev(matoptlevel) = nmat_dsyev(matoptlevel)+1 
   CPUTIME_mat_dsyev(matoptlevel) = CPUTIME_mat_dsyev(matoptlevel) + deltacpu
   WALLTIME_mat_dsyev(matoptlevel) = WALLTIME_mat_dsyev(matoptlevel) + deltawall
CASE(JOB_mat_dsyevx)
   nmat_dsyevx(matoptlevel) = nmat_dsyevx(matoptlevel)+1 
   CPUTIME_mat_dsyevx(matoptlevel) = CPUTIME_mat_dsyevx(matoptlevel) + deltacpu
   WALLTIME_mat_dsyevx(matoptlevel) = WALLTIME_mat_dsyevx(matoptlevel) + deltawall
CASE(JOB_mat_identity)
   nmat_identity(matoptlevel) = nmat_identity(matoptlevel)+1 
   CPUTIME_mat_identity(matoptlevel) = CPUTIME_mat_identity(matoptlevel) + deltacpu
   WALLTIME_mat_identity(matoptlevel) = WALLTIME_mat_identity(matoptlevel) + deltawall
CASE(JOB_mat_setlowertriangular_zero)
   nmat_setlowertriangular_zero(matoptlevel) = nmat_setlowertriangular_zero(matoptlevel)+1 
   CPUTIME_mat_setlowertriangular_zero(matoptlevel) = CPUTIME_mat_setlowertriangular_zero(matoptlevel) + deltacpu
   WALLTIME_mat_setlowertriangular_zero(matoptlevel) = WALLTIME_mat_setlowertriangular_zero(matoptlevel) + deltawall
CASE(JOB_mat_write_to_disk)
   nmat_write_to_disk(matoptlevel) = nmat_write_to_disk(matoptlevel)+1 
   CPUTIME_mat_write_to_disk(matoptlevel) = CPUTIME_mat_write_to_disk(matoptlevel) + deltacpu
   WALLTIME_mat_write_to_disk(matoptlevel) = WALLTIME_mat_write_to_disk(matoptlevel) + deltawall
CASE(JOB_mat_read_from_disk)
   nmat_read_from_disk(matoptlevel) = nmat_read_from_disk(matoptlevel)+1 
   CPUTIME_mat_read_from_disk(matoptlevel) = CPUTIME_mat_read_from_disk(matoptlevel) + deltacpu
   WALLTIME_mat_read_from_disk(matoptlevel) = WALLTIME_mat_read_from_disk(matoptlevel) + deltawall
CASE(JOB_II_get_overlap)
   nII_get_overlap(matoptlevel) = nII_get_overlap(matoptlevel)+1 
   CPUTIME_II_get_overlap(matoptlevel) = CPUTIME_II_get_overlap(matoptlevel) + deltacpu
   WALLTIME_II_get_overlap(matoptlevel) = WALLTIME_II_get_overlap(matoptlevel) + deltawall
CASE(JOB_II_get_magderivOverlap)
   nII_get_magderivOverlap(matoptlevel) = nII_get_magderivOverlap(matoptlevel)+1 
   CPUTIME_II_get_magderivOverlap(matoptlevel) = CPUTIME_II_get_magderivOverlap(matoptlevel) + deltacpu
   WALLTIME_II_get_magderivOverlap(matoptlevel) = WALLTIME_II_get_magderivOverlap(matoptlevel) + deltawall
CASE(JOB_II_get_maggradOverlap)
   nII_get_maggradOverlap(matoptlevel) = nII_get_maggradOverlap(matoptlevel)+1 
   CPUTIME_II_get_maggradOverlap(matoptlevel) = CPUTIME_II_get_maggradOverlap(matoptlevel) + deltacpu
   WALLTIME_II_get_maggradOverlap(matoptlevel) = WALLTIME_II_get_maggradOverlap(matoptlevel) + deltawall
CASE(JOB_II_get_kinetic)
   nII_get_kinetic(matoptlevel) = nII_get_kinetic(matoptlevel)+1 
   CPUTIME_II_get_kinetic(matoptlevel) = CPUTIME_II_get_kinetic(matoptlevel) + deltacpu
   WALLTIME_II_get_kinetic(matoptlevel) = WALLTIME_II_get_kinetic(matoptlevel) + deltawall
CASE(JOB_II_get_nucel_mat)
   nII_get_nucel_mat(matoptlevel) = nII_get_nucel_mat(matoptlevel)+1 
   CPUTIME_II_get_nucel_mat(matoptlevel) = CPUTIME_II_get_nucel_mat(matoptlevel) + deltacpu
   WALLTIME_II_get_nucel_mat(matoptlevel) = WALLTIME_II_get_nucel_mat(matoptlevel) + deltawall
CASE(JOB_II_GET_NN_GRADIENT)
   nII_GET_NN_GRADIENT(matoptlevel) = nII_GET_NN_GRADIENT(matoptlevel)+1 
   CPUTIME_II_GET_NN_GRADIENT(matoptlevel) = CPUTIME_II_GET_NN_GRADIENT(matoptlevel) + deltacpu
   WALLTIME_II_GET_NN_GRADIENT(matoptlevel) = WALLTIME_II_GET_NN_GRADIENT(matoptlevel) + deltawall
CASE(JOB_II_precalc_ScreenMat)
   nII_precalc_ScreenMat(matoptlevel) = nII_precalc_ScreenMat(matoptlevel)+1 
   CPUTIME_II_precalc_ScreenMat(matoptlevel) = CPUTIME_II_precalc_ScreenMat(matoptlevel) + deltacpu
   WALLTIME_II_precalc_ScreenMat(matoptlevel) = WALLTIME_II_precalc_ScreenMat(matoptlevel) + deltawall
CASE(JOB_II_get_ne_gradient)
   nII_get_ne_gradient(matoptlevel) = nII_get_ne_gradient(matoptlevel)+1 
   CPUTIME_II_get_ne_gradient(matoptlevel) = CPUTIME_II_get_ne_gradient(matoptlevel) + deltacpu
   WALLTIME_II_get_ne_gradient(matoptlevel) = WALLTIME_II_get_ne_gradient(matoptlevel) + deltawall
CASE(JOB_II_get_geoderivOverlap)
   nII_get_geoderivOverlap(matoptlevel) = nII_get_geoderivOverlap(matoptlevel)+1 
   CPUTIME_II_get_geoderivOverlap(matoptlevel) = CPUTIME_II_get_geoderivOverlap(matoptlevel) + deltacpu
   WALLTIME_II_get_geoderivOverlap(matoptlevel) = WALLTIME_II_get_geoderivOverlap(matoptlevel) + deltawall
CASE(JOB_II_get_geoderivExchange)
   nII_get_geoderivExchange(matoptlevel) = nII_get_geoderivExchange(matoptlevel)+1 
   CPUTIME_II_get_geoderivExchange(matoptlevel) = CPUTIME_II_get_geoderivExchange(matoptlevel) + deltacpu
   WALLTIME_II_get_geoderivExchange(matoptlevel) = WALLTIME_II_get_geoderivExchange(matoptlevel) + deltawall
CASE(JOB_II_get_geoderivCoulomb)
   nII_get_geoderivCoulomb(matoptlevel) = nII_get_geoderivCoulomb(matoptlevel)+1 
   CPUTIME_II_get_geoderivCoulomb(matoptlevel) = CPUTIME_II_get_geoderivCoulomb(matoptlevel) + deltacpu
   WALLTIME_II_get_geoderivCoulomb(matoptlevel) = WALLTIME_II_get_geoderivCoulomb(matoptlevel) + deltawall
CASE(JOB_II_GET_KINETIC_GRADIENT)
   nII_GET_KINETIC_GRADIENT(matoptlevel) = nII_GET_KINETIC_GRADIENT(matoptlevel)+1 
   CPUTIME_II_GET_KINETIC_GRADIENT(matoptlevel) = CPUTIME_II_GET_KINETIC_GRADIENT(matoptlevel) + deltacpu
   WALLTIME_II_GET_KINETIC_GRADIENT(matoptlevel) = WALLTIME_II_GET_KINETIC_GRADIENT(matoptlevel) + deltawall
CASE(JOB_II_GET_REORTHONORMALIZATION)
   nII_GET_REORTHONORMALIZATION(matoptlevel) = nII_GET_REORTHONORMALIZATION(matoptlevel)+1 
   CPUTIME_II_GET_REORTHONORMALIZATION(matoptlevel) = CPUTIME_II_GET_REORTHONORMALIZATION(matoptlevel) + deltacpu
   WALLTIME_II_GET_REORTHONORMALIZATION(matoptlevel) = WALLTIME_II_GET_REORTHONORMALIZATION(matoptlevel) + deltawall
CASE(JOB_II_GET_REORTHONORMALIZATION2)
   nII_GET_REORTHONORMALIZATION2(matoptlevel) = nII_GET_REORTHONORMALIZATION2(matoptlevel)+1 
   CPUTIME_II_GET_REORTHONORMALIZATION2(matoptlevel) = CPUTIME_II_GET_REORTHONORMALIZATION2(matoptlevel) + deltacpu
   WALLTIME_II_GET_REORTHONORMALIZATION2(matoptlevel) = WALLTIME_II_GET_REORTHONORMALIZATION2(matoptlevel) + deltawall
CASE(JOB_II_get_J_gradient)
   nII_get_J_gradient(matoptlevel) = nII_get_J_gradient(matoptlevel)+1 
   CPUTIME_II_get_J_gradient(matoptlevel) = CPUTIME_II_get_J_gradient(matoptlevel) + deltacpu
   WALLTIME_II_get_J_gradient(matoptlevel) = WALLTIME_II_get_J_gradient(matoptlevel) + deltawall
CASE(JOB_II_GET_K_GRADIENT)
   nII_GET_K_GRADIENT(matoptlevel) = nII_GET_K_GRADIENT(matoptlevel)+1 
   CPUTIME_II_GET_K_GRADIENT(matoptlevel) = CPUTIME_II_GET_K_GRADIENT(matoptlevel) + deltacpu
   WALLTIME_II_GET_K_GRADIENT(matoptlevel) = WALLTIME_II_GET_K_GRADIENT(matoptlevel) + deltawall
CASE(JOB_II_get_nucpot)
   nII_get_nucpot(matoptlevel) = nII_get_nucpot(matoptlevel)+1 
   CPUTIME_II_get_nucpot(matoptlevel) = CPUTIME_II_get_nucpot(matoptlevel) + deltacpu
   WALLTIME_II_get_nucpot(matoptlevel) = WALLTIME_II_get_nucpot(matoptlevel) + deltawall
CASE(JOB_II_get_prop)
   nII_get_prop(matoptlevel) = nII_get_prop(matoptlevel)+1 
   CPUTIME_II_get_prop(matoptlevel) = CPUTIME_II_get_prop(matoptlevel) + deltacpu
   WALLTIME_II_get_prop(matoptlevel) = WALLTIME_II_get_prop(matoptlevel) + deltawall
CASE(JOB_II_get_prop_expval)
   nII_get_prop_expval(matoptlevel) = nII_get_prop_expval(matoptlevel)+1 
   CPUTIME_II_get_prop_expval(matoptlevel) = CPUTIME_II_get_prop_expval(matoptlevel) + deltacpu
   WALLTIME_II_get_prop_expval(matoptlevel) = WALLTIME_II_get_prop_expval(matoptlevel) + deltawall
CASE(JOB_II_GET_EXCHANGE_MAT)
   nII_GET_EXCHANGE_MAT(matoptlevel) = nII_GET_EXCHANGE_MAT(matoptlevel)+1 
   CPUTIME_II_GET_EXCHANGE_MAT(matoptlevel) = CPUTIME_II_GET_EXCHANGE_MAT(matoptlevel) + deltacpu
   WALLTIME_II_GET_EXCHANGE_MAT(matoptlevel) = WALLTIME_II_GET_EXCHANGE_MAT(matoptlevel) + deltawall
CASE(JOB_II_get_coulomb_mat)
   nII_get_coulomb_mat(matoptlevel) = nII_get_coulomb_mat(matoptlevel)+1 
   CPUTIME_II_get_coulomb_mat(matoptlevel) = CPUTIME_II_get_coulomb_mat(matoptlevel) + deltacpu
   WALLTIME_II_get_coulomb_mat(matoptlevel) = WALLTIME_II_get_coulomb_mat(matoptlevel) + deltawall
CASE(JOB_II_get_coulomb_and_exchange_mat)
   nII_get_coulomb_and_exchange_mat(matoptlevel) = nII_get_coulomb_and_exchange_mat(matoptlevel)+1 
   CPUTIME_II_get_coulomb_and_exchange_mat(matoptlevel) = CPUTIME_II_get_coulomb_and_exchange_mat(matoptlevel) + deltacpu
   WALLTIME_II_get_coulomb_and_exchange_mat(matoptlevel) = WALLTIME_II_get_coulomb_and_exchange_mat(matoptlevel) + deltawall
CASE(JOB_II_get_magderiv_4center_eri)
   nII_get_magderiv_4center_eri(matoptlevel) = nII_get_magderiv_4center_eri(matoptlevel)+1 
   CPUTIME_II_get_magderiv_4center_eri(matoptlevel) = CPUTIME_II_get_magderiv_4center_eri(matoptlevel) + deltacpu
   WALLTIME_II_get_magderiv_4center_eri(matoptlevel) = WALLTIME_II_get_magderiv_4center_eri(matoptlevel) + deltawall
CASE(JOB_II_get_magderivK_4center_eri)
   nII_get_magderivK_4center_eri(matoptlevel) = nII_get_magderivK_4center_eri(matoptlevel)+1 
   CPUTIME_II_get_magderivK_4center_eri(matoptlevel) = CPUTIME_II_get_magderivK_4center_eri(matoptlevel) + deltacpu
   WALLTIME_II_get_magderivK_4center_eri(matoptlevel) = WALLTIME_II_get_magderivK_4center_eri(matoptlevel) + deltawall
CASE(JOB_II_GET_MAGDERIVJ_4CENTER_ERI)
   nII_GET_MAGDERIVJ_4CENTER_ERI(matoptlevel) = nII_GET_MAGDERIVJ_4CENTER_ERI(matoptlevel)+1 
   CPUTIME_II_GET_MAGDERIVJ_4CENTER_ERI(matoptlevel) = CPUTIME_II_GET_MAGDERIVJ_4CENTER_ERI(matoptlevel) + deltacpu
   WALLTIME_II_GET_MAGDERIVJ_4CENTER_ERI(matoptlevel) = WALLTIME_II_GET_MAGDERIVJ_4CENTER_ERI(matoptlevel) + deltawall
CASE(JOB_II_GET_DECPACKED4CENTER_J_ERI)
   nII_GET_DECPACKED4CENTER_J_ERI(matoptlevel) = nII_GET_DECPACKED4CENTER_J_ERI(matoptlevel)+1 
   CPUTIME_II_GET_DECPACKED4CENTER_J_ERI(matoptlevel) = CPUTIME_II_GET_DECPACKED4CENTER_J_ERI(matoptlevel) + deltacpu
   WALLTIME_II_GET_DECPACKED4CENTER_J_ERI(matoptlevel) = WALLTIME_II_GET_DECPACKED4CENTER_J_ERI(matoptlevel) + deltawall
CASE(JOB_II_GET_DECPACKED4CENTER_K_ERI)
   nII_GET_DECPACKED4CENTER_K_ERI(matoptlevel) = nII_GET_DECPACKED4CENTER_K_ERI(matoptlevel)+1 
   CPUTIME_II_GET_DECPACKED4CENTER_K_ERI(matoptlevel) = CPUTIME_II_GET_DECPACKED4CENTER_K_ERI(matoptlevel) + deltacpu
   WALLTIME_II_GET_DECPACKED4CENTER_K_ERI(matoptlevel) = WALLTIME_II_GET_DECPACKED4CENTER_K_ERI(matoptlevel) + deltawall
CASE(JOB_II_get_xc_Fock_mat)
   nII_get_xc_Fock_mat(matoptlevel) = nII_get_xc_Fock_mat(matoptlevel)+1 
   CPUTIME_II_get_xc_Fock_mat(matoptlevel) = CPUTIME_II_get_xc_Fock_mat(matoptlevel) + deltacpu
   WALLTIME_II_get_xc_Fock_mat(matoptlevel) = WALLTIME_II_get_xc_Fock_mat(matoptlevel) + deltawall
CASE(JOB_II_get_xc_geoderiv_molgrad)
   nII_get_xc_geoderiv_molgrad(matoptlevel) = nII_get_xc_geoderiv_molgrad(matoptlevel)+1 
   CPUTIME_II_get_xc_geoderiv_molgrad(matoptlevel) = CPUTIME_II_get_xc_geoderiv_molgrad(matoptlevel) + deltacpu
   WALLTIME_II_get_xc_geoderiv_molgrad(matoptlevel) = WALLTIME_II_get_xc_geoderiv_molgrad(matoptlevel) + deltawall
CASE(JOB_II_get_xc_linrsp)
   nII_get_xc_linrsp(matoptlevel) = nII_get_xc_linrsp(matoptlevel)+1 
   CPUTIME_II_get_xc_linrsp(matoptlevel) = CPUTIME_II_get_xc_linrsp(matoptlevel) + deltacpu
   WALLTIME_II_get_xc_linrsp(matoptlevel) = WALLTIME_II_get_xc_linrsp(matoptlevel) + deltawall
CASE(JOB_II_get_xc_quadrsp)
   nII_get_xc_quadrsp(matoptlevel) = nII_get_xc_quadrsp(matoptlevel)+1 
   CPUTIME_II_get_xc_quadrsp(matoptlevel) = CPUTIME_II_get_xc_quadrsp(matoptlevel) + deltacpu
   WALLTIME_II_get_xc_quadrsp(matoptlevel) = WALLTIME_II_get_xc_quadrsp(matoptlevel) + deltawall
CASE(JOB_II_get_xc_magderiv_kohnsham_mat)
   nII_get_xc_magderiv_kohnsham_mat(matoptlevel) = nII_get_xc_magderiv_kohnsham_mat(matoptlevel)+1 
   CPUTIME_II_get_xc_magderiv_kohnsham_mat(matoptlevel) = CPUTIME_II_get_xc_magderiv_kohnsham_mat(matoptlevel) + deltacpu
   WALLTIME_II_get_xc_magderiv_kohnsham_mat(matoptlevel) = WALLTIME_II_get_xc_magderiv_kohnsham_mat(matoptlevel) + deltawall
CASE(JOB_II_get_xc_magderiv_linrsp)
   nII_get_xc_magderiv_linrsp(matoptlevel) = nII_get_xc_magderiv_linrsp(matoptlevel)+1 
   CPUTIME_II_get_xc_magderiv_linrsp(matoptlevel) = CPUTIME_II_get_xc_magderiv_linrsp(matoptlevel) + deltacpu
   WALLTIME_II_get_xc_magderiv_linrsp(matoptlevel) = WALLTIME_II_get_xc_magderiv_linrsp(matoptlevel) + deltawall
CASE(JOB_II_get_xc_geoderiv_FxDgrad)
   nII_get_xc_geoderiv_FxDgrad(matoptlevel) = nII_get_xc_geoderiv_FxDgrad(matoptlevel)+1 
   CPUTIME_II_get_xc_geoderiv_FxDgrad(matoptlevel) = CPUTIME_II_get_xc_geoderiv_FxDgrad(matoptlevel) + deltacpu
   WALLTIME_II_get_xc_geoderiv_FxDgrad(matoptlevel) = WALLTIME_II_get_xc_geoderiv_FxDgrad(matoptlevel) + deltawall
CASE(JOB_II_get_xc_geoderiv_GxDgrad)
   nII_get_xc_geoderiv_GxDgrad(matoptlevel) = nII_get_xc_geoderiv_GxDgrad(matoptlevel)+1 
   CPUTIME_II_get_xc_geoderiv_GxDgrad(matoptlevel) = CPUTIME_II_get_xc_geoderiv_GxDgrad(matoptlevel) + deltacpu
   WALLTIME_II_get_xc_geoderiv_GxDgrad(matoptlevel) = WALLTIME_II_get_xc_geoderiv_GxDgrad(matoptlevel) + deltawall
CASE(JOB_homolumo)
   nhomolumo(matoptlevel) = nhomolumo(matoptlevel)+1 
   CPUTIME_homolumo(matoptlevel) = CPUTIME_homolumo(matoptlevel) + deltacpu
   WALLTIME_homolumo(matoptlevel) = WALLTIME_homolumo(matoptlevel) + deltawall
CASE DEFAULT
   WRITE (LSIUNIT,'(/,3A,/)') ' JobID "',JOB,&
        & '" not recognized in add_time_to_matjob'
   CALL lsQUIT('Illegal JobID in add_time_to_matjob',-1)
END SELECT
IF(job.GT.0.AND.job.LT.21)THEN
   CPUTIME_TOTAL_mat(matoptlevel) = CPUTIME_TOTAL_mat(matoptlevel) + deltacpu
   WALLTIME_TOTAL_mat(matoptlevel) = WALLTIME_TOTAL_mat(matoptlevel) + deltawall
ELSE
   CPUTIME_TOTAL_II(matoptlevel) = CPUTIME_TOTAL_II(matoptlevel) + deltacpu
   WALLTIME_TOTAL_II(matoptlevel) = WALLTIME_TOTAL_II(matoptlevel) + deltawall
ENDIF
end subroutine add_time_to_job

subroutine print_timers(lupri)
  implicit none
  integer :: lupri
#ifdef VAR_TIME
  integer :: I
  character(len=3) :: ID
  DO I=1,3
     IF((I.EQ.1).AND.(.NOT.LEVEL1ACTIVE))CYCLE
     IF((I.EQ.2).AND.(.NOT.LEVEL2ACTIVE))CYCLE
     write(lupri,'(6X,A)')     '==============================================================='
     write(lupri,'(6X,A)')     'Overall Timings for the Matrix Operations and Integral routines'
     IF(I.EQ.1)write(lupri,'(6X,A)')'         in Level 1 The Atomic Calculation                     '
     IF(I.EQ.2)write(lupri,'(6X,A)')'         in Level 2 The Valence Calculation                    '
     IF(I.EQ.3)write(lupri,'(6X,A)')'         in Level 3 The Full molecular Calculation             '
     write(lupri,'(6X,A18,A8,2A18)')   'Matrix Operation  ','# calls',' CPU Time', 'Wall Time'
     write(lupri,'(6X,A)')     '==============================================================='
     IF(I.EQ.1) ID = 'L1 '
     IF(I.EQ.2) ID = 'L2 '
     IF(I.EQ.3) ID = 'L3 '
     IF(nmat_set_from_full(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_set_from_full ',&
          &nmat_set_from_full(I),CPUTIME_mat_set_from_full(I),WALLTIME_mat_set_from_full(I)
     IF(nmat_to_full(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_to_full       ',&
          &nmat_to_full(I),CPUTIME_mat_to_full(I),WALLTIME_mat_to_full(I)
     IF(nmat_trans(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_trans         ',&
          &nmat_trans(I),CPUTIME_mat_trans(I),WALLTIME_mat_trans(I)
     IF(nmat_assign(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_assign        ',&
          &nmat_assign(I),CPUTIME_mat_assign(I),WALLTIME_mat_assign(I)
     IF(nmat_copy(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_copy          ',&
          &nmat_copy(I),CPUTIME_mat_copy(I),WALLTIME_mat_copy(I)
     IF(nmat_tr(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_tr            ',&
          &nmat_tr(I),CPUTIME_mat_tr(I),WALLTIME_mat_tr(I)
     IF(nmat_trAB(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_trAB          ',&
          &nmat_trAB(I),CPUTIME_mat_trAB(I),WALLTIME_mat_trAB(I)
     IF(nmat_mul(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_mul           ',&
          &nmat_mul(I),CPUTIME_mat_mul(I),WALLTIME_mat_mul(I)
     IF(nmat_add(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_add           ',&
          &nmat_add(I),CPUTIME_mat_add(I),WALLTIME_mat_add(I)
     IF(nmat_daxpy(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_daxpy         ',&
          &nmat_daxpy(I),CPUTIME_mat_daxpy(I),WALLTIME_mat_daxpy(I)
     IF(nmat_dotproduct(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_dotproduct    ',&
          &nmat_dotproduct(I),CPUTIME_mat_dotproduct(I),WALLTIME_mat_dotproduct(I)
     IF(nmat_sqnorm2(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_sqnorm2       ',&
          &nmat_sqnorm2(I),CPUTIME_mat_sqnorm2(I),WALLTIME_mat_sqnorm2(I)
     IF(nmat_diag_f(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_diag_f        ',&
          &nmat_diag_f(I),CPUTIME_mat_diag_f(I),WALLTIME_mat_diag_f(I)
     IF(nmat_create_block(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_create_block  ',&
          &nmat_create_block(I),CPUTIME_mat_create_block(I),WALLTIME_mat_create_block(I)
     IF(nmat_add_block(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_add_block     ',&
          &nmat_add_block(I),CPUTIME_mat_add_block(I),WALLTIME_mat_add_block(I)
     IF(nmat_retrieve_block(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_retrieve_block',&
          &nmat_retrieve_block(I),CPUTIME_mat_retrieve_block(I),WALLTIME_mat_retrieve_block(I)
     IF(nmat_scal(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_scal          ',&
          &nmat_scal(I),CPUTIME_mat_scal(I),WALLTIME_mat_scal(I)
     IF(nmat_zero(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_zero          ',&
          &nmat_zero(I),CPUTIME_mat_zero(I),WALLTIME_mat_zero(I)
     IF(nmat_density_from_orbs(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_density_from_orbs          ',&
          &nmat_density_from_orbs(I),CPUTIME_mat_density_from_orbs(I),WALLTIME_mat_density_from_orbs(I)
     IF(nmat_scal_dia_vec(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_scal_dia_vec          ',&
          &nmat_scal_dia_vec(I),CPUTIME_mat_scal_dia_vec(I),WALLTIME_mat_scal_dia_vec(I)
     IF(nmat_scal_dia(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_scal_dia          ',&
          &nmat_scal_dia(I),CPUTIME_mat_scal_dia(I),WALLTIME_mat_scal_dia(I)
     IF(nmat_inv(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_inv          ',&
          &nmat_inv(I),CPUTIME_mat_inv(I),WALLTIME_mat_inv(I)
     IF(nmat_dsyev(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_dsyev          ',&
          &nmat_dsyev(I),CPUTIME_mat_dsyev(I),WALLTIME_mat_dsyev(I)
     IF(nmat_dsyevx(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_dsyevx          ',&
          &nmat_dsyevx(I),CPUTIME_mat_dsyevx(I),WALLTIME_mat_dsyevx(I)
     IF(nmat_identity(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_identity          ',&
          &nmat_identity(I),CPUTIME_mat_identity(I),WALLTIME_mat_identity(I)
     IF(nmat_setlowertriangular_zero(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_lowertriang_z',&
          &nmat_setlowertriangular_zero(I),CPUTIME_mat_setlowertriangular_zero(I),&
          &WALLTIME_mat_setlowertriangular_zero(I)
     IF(nmat_chol(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_chol          ',&
          &nmat_chol(I),CPUTIME_mat_chol(I),WALLTIME_mat_chol(I)
     IF(nmat_write_to_disk(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_write_to_disk ',&
          &nmat_write_to_disk(I),CPUTIME_mat_write_to_disk(I),WALLTIME_mat_write_to_disk(I)
     IF(nmat_read_from_disk(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'mat_read_from_disk',&
          &nmat_read_from_disk(I),CPUTIME_mat_read_from_disk(I),WALLTIME_mat_read_from_disk(I)
     write(lupri,'(6X,A)')     '=============================================================='
     write(lupri,'(3X,A3,A18,8X,2F18.4)') ID,'TOTAL MAT         ',CPUTIME_TOTAL_mat(I),WALLTIME_TOTAL_mat(I)
     write(lupri,'(6X,A)')     '=============================================================='
     IF(nII_get_overlap(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_overlap',&
          &nII_get_overlap(I),CPUTIME_II_get_overlap(I),WALLTIME_II_get_overlap(I)
     IF(nII_get_magderivOverlap(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_magderivOverlap'&
          &,nII_get_magderivOverlap(I),CPUTIME_II_get_magderivOverlap(I),WALLTIME_II_get_magderivOverlap(I)
     IF(nII_get_maggradOverlap(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_maggradOverlap'&
          &,nII_get_maggradOverlap(I),CPUTIME_II_get_maggradOverlap(I),WALLTIME_II_get_maggradOverlap(I)
     IF(nII_get_kinetic(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_kinetic'&
          &,nII_get_kinetic(I),CPUTIME_II_get_kinetic(I),WALLTIME_II_get_kinetic(I)
     IF(nII_get_nucel_mat(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_nucel_mat'&
          &,nII_get_nucel_mat(I),CPUTIME_II_get_nucel_mat(I),WALLTIME_II_get_nucel_mat(I)
     IF(nII_GET_NN_GRADIENT(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_GET_NN_GRADIENT'&
          &,nII_GET_NN_GRADIENT(I),CPUTIME_II_GET_NN_GRADIENT(I),WALLTIME_II_GET_NN_GRADIENT(I)
     IF(nII_precalc_ScreenMat(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_precalc_ScreenMat'&
          &,nII_precalc_ScreenMat(I),CPUTIME_II_precalc_ScreenMat(I),WALLTIME_II_precalc_ScreenMat(I)
     IF(nII_get_ne_gradient(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_ne_gradient'&
          &,nII_get_ne_gradient(I),CPUTIME_II_get_ne_gradient(I),WALLTIME_II_get_ne_gradient(I)
     IF(nII_get_geoderivOverlap(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_geoderivOverlap'&
          &,nII_get_geoderivOverlap(I),CPUTIME_II_get_geoderivOverlap(I),WALLTIME_II_get_geoderivOverlap(I)
     IF(nII_get_geoderivExchange(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_geoderivExchange'&
          &,nII_get_geoderivExchange(I),CPUTIME_II_get_geoderivExchange(I),WALLTIME_II_get_geoderivExchange(I)
     IF(nII_get_geoderivCoulomb(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_geoderivCoulomb'&
          &,nII_get_geoderivCoulomb(I),CPUTIME_II_get_geoderivCoulomb(I),WALLTIME_II_get_geoderivCoulomb(I)
     IF(nII_GET_KINETIC_GRADIENT(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_GET_KINETIC_GRADIENT'&
          &,nII_GET_KINETIC_GRADIENT(I),CPUTIME_II_GET_KINETIC_GRADIENT(I),WALLTIME_II_GET_KINETIC_GRADIENT(I)
     IF(nII_GET_REORTHONORMALIZATION(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_GET_REORTHONORMALIZATION'&
          &,nII_GET_REORTHONORMALIZATION(I),CPUTIME_II_GET_REORTHONORMALIZATION(I),WALLTIME_II_GET_REORTHONORMALIZATION(I)
     IF(nII_GET_REORTHONORMALIZATION2(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_GET_REORTHONORMALIZATION2'&
          &,nII_GET_REORTHONORMALIZATION2(I),CPUTIME_II_GET_REORTHONORMALIZATION2(I),WALLTIME_II_GET_REORTHONORMALIZATION2(I)
     IF(nII_get_J_gradient(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_J_gradient'&
          &,nII_get_J_gradient(I),CPUTIME_II_get_J_gradient(I),WALLTIME_II_get_J_gradient(I)
     IF(nII_get_K_gradient(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_K_gradient'&
          &,nII_get_K_gradient(I),CPUTIME_II_get_K_gradient(I),WALLTIME_II_get_K_gradient(I)
     IF(nII_GET_K_GRADIENT(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_GET_K_GRADIENT'&
          &,nII_GET_K_GRADIENT(I),CPUTIME_II_GET_K_GRADIENT(I),WALLTIME_II_GET_K_GRADIENT(I)
     IF(nII_get_nucpot(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_nucpot'&
          &,nII_get_nucpot(I),CPUTIME_II_get_nucpot(I),WALLTIME_II_get_nucpot(I)
     IF(nII_get_prop(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_prop'&
          &,nII_get_prop(I),CPUTIME_II_get_prop(I),WALLTIME_II_get_prop(I)
     IF(nII_get_prop_expval(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_prop_expval'&
          &,nII_get_prop_expval(I),CPUTIME_II_get_prop_expval(I),WALLTIME_II_get_prop_expval(I)
     IF(nII_GET_EXCHANGE_MAT(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_GET_EXCHANGE_MAT'&
          &,nII_GET_EXCHANGE_MAT(I),CPUTIME_II_GET_EXCHANGE_MAT(I),WALLTIME_II_GET_EXCHANGE_MAT(I)
     IF(nII_get_coulomb_mat(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_coulomb_mat'&
          &,nII_get_coulomb_mat(I),CPUTIME_II_get_coulomb_mat(I),WALLTIME_II_get_coulomb_mat(I)
     IF(nII_get_coulomb_and_exchange_mat(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_coulomb_and_exchange_mat'&
          &,nII_get_coulomb_and_exchange_mat(I),CPUTIME_II_get_coulomb_and_exchange_mat(I),&
          & WALLTIME_II_get_coulomb_and_exchange_mat(I)
     IF(nII_get_magderiv_4center_eri(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_magderiv_4center_eri'&
          &,nII_get_magderiv_4center_eri(I),CPUTIME_II_get_magderiv_4center_eri(I),WALLTIME_II_get_magderiv_4center_eri(I)
     IF(nII_get_magderivK_4center_eri(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_magderivK_4center_eri'&
          &,nII_get_magderivK_4center_eri(I),CPUTIME_II_get_magderivK_4center_eri(I),WALLTIME_II_get_magderivK_4center_eri(I)
     IF(nII_GET_MAGDERIVJ_4CENTER_ERI(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_GET_MAGDERIVJ_4CENTER_ERI'&
          &,nII_GET_MAGDERIVJ_4CENTER_ERI(I),CPUTIME_II_GET_MAGDERIVJ_4CENTER_ERI(I),WALLTIME_II_GET_MAGDERIVJ_4CENTER_ERI(I)
     IF(nII_GET_DECPACKED4CENTER_J_ERI(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_GET_DECPACKED4CENTER_J_ERI'&
          &,nII_GET_DECPACKED4CENTER_J_ERI(I),CPUTIME_II_GET_DECPACKED4CENTER_J_ERI(I),WALLTIME_II_GET_DECPACKED4CENTER_J_ERI(I)
     IF(nII_GET_DECPACKED4CENTER_K_ERI(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_GET_DECPACKED4CENTER_K_ERI'&
          &,nII_GET_DECPACKED4CENTER_K_ERI(I),CPUTIME_II_GET_DECPACKED4CENTER_K_ERI(I),WALLTIME_II_GET_DECPACKED4CENTER_K_ERI(I)
     IF(nII_get_xc_Fock_mat(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_xc_Fock_mat'&
          &,nII_get_xc_Fock_mat(I),CPUTIME_II_get_xc_Fock_mat(I),WALLTIME_II_get_xc_Fock_mat(I)
     IF(nII_get_xc_geoderiv_molgrad(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_xc_geoderiv_molgrad'&
          &,nII_get_xc_geoderiv_molgrad(I),CPUTIME_II_get_xc_geoderiv_molgrad(I),WALLTIME_II_get_xc_geoderiv_molgrad(I)
     IF(nII_get_xc_linrsp(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_xc_linrsp'&
          &,nII_get_xc_linrsp(I),CPUTIME_II_get_xc_linrsp(I),WALLTIME_II_get_xc_linrsp(I)
     IF(nII_get_xc_quadrsp(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_xc_quadrsp'&
          &,nII_get_xc_quadrsp(I),CPUTIME_II_get_xc_quadrsp(I),WALLTIME_II_get_xc_quadrsp(I)
     IF(nII_get_xc_magderiv_kohnsham_mat(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_xc_magderiv_kohnsham_mat'&
          &,nII_get_xc_magderiv_kohnsham_mat(I),CPUTIME_II_get_xc_magderiv_kohnsham_mat(I),&
          & WALLTIME_II_get_xc_magderiv_kohnsham_mat(I)
     IF(nII_get_xc_magderiv_linrsp(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_xc_magderiv_linrsp'&
          &,nII_get_xc_magderiv_linrsp(I),CPUTIME_II_get_xc_magderiv_linrsp(I),WALLTIME_II_get_xc_magderiv_linrsp(I)
     IF(nII_get_xc_geoderiv_FxDgrad(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_xc_geoderiv_FxDgrad'&
          &,nII_get_xc_geoderiv_FxDgrad(I),CPUTIME_II_get_xc_geoderiv_FxDgrad(I),WALLTIME_II_get_xc_geoderiv_FxDgrad(I)
     IF(nII_get_xc_geoderiv_GxDgrad(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'II_get_xc_geoderiv_GxDgrad'&
          &,nII_get_xc_geoderiv_GxDgrad(I),CPUTIME_II_get_xc_geoderiv_GxDgrad(I),WALLTIME_II_get_xc_geoderiv_GxDgrad(I)
     IF(nhomolumo(I).GT.0)write(lupri,'(3X,A3,A18,I8,2F18.4)') ID,'homolumo'&
          &,nhomolumo(I),CPUTIME_homolumo(I),WALLTIME_homolumo(I)
     write(lupri,'(6X,A)')     '=============================================================='
     write(lupri,'(3X,A3,A18,8X,2F18.4)') ID,'TOTAL INTEGRAL    ',CPUTIME_TOTAL_II(I),WALLTIME_TOTAL_II(I)
     write(lupri,'(6X,A)')     '=============================================================='
  ENDDO
#endif
end subroutine print_timers

!> \brief take time
!> \author T. Kjaergaard
!> \date 2010
!> \param text label to print along with timings
!> \param CPUTIME the cpu time
!> \param WALLTIME the wall time
!> \param lupri the logical unit number 
SUBROUTINE LSTIMER(TEXT,CPUTIME,WALLTIME,LUPRI)
implicit none
INTEGER           :: LUPRI,length
CHARACTER*(*)     :: TEXT
CHARACTER(len=15) :: PRINTTEXT
REAL(REALK)       :: TIME1,TIME2,DELTAWALL,CPUTIME,WALLTIME,DELTACPU

INTEGER :: IUNIT

IF (LUPRI.EQ.-1) THEN
   IUNIT = LSIUNIT
ELSE
   IUNIT = LUPRI
ENDIF

length = LEN(TEXT)
IF(length .GT. 15) CALL LSQUIT('TEXTLENGTH PROVIDED TO LSTIMER IS LIMITED TO 15',lupri)
PRINTTEXT='               '
PRINTTEXT(1:length) =  TEXT(1:length) 
IF (PRINTTEXT(1:5) .EQ. 'START') THEN
   CALL LS_GETTIM(CPUTIME,WALLTIME)
ELSE
   CALL LS_GETTIM(TIME1,TIME2)
   DELTACPU=TIME1-CPUTIME
   DELTAWALL=TIME2-WALLTIME
   IF (LSTIME_PRINT) THEN
!$OMP CRITICAL
      CALL LS_TIMTXT('>>>  CPU Time used in '//PRINTTEXT//' is',DELTACPU,IUNIT)
      CALL LS_TIMTXT('>>> wall Time used in '//PRINTTEXT//' is',DELTAWALL,IUNIT)
!$OMP END CRITICAL
   ENDIF
   WALLTIME=TIME2
   CPUTIME=TIME1
   CALL LS_FLSHFO(IUNIT)
ENDIF

END SUBROUTINE LSTIMER

SUBROUTINE SET_LSTIME_PRINT(TIME_PRINT)
implicit none
LOGICAL, INTENT(IN) :: TIME_PRINT
LSTIME_PRINT = TIME_PRINT
END SUBROUTINE SET_LSTIME_PRINT

!> \brief set the lsiunit integer
!> \author T. Kjaergaard
!> \date 2010
!> \param iunit
SUBROUTINE SET_LSIUNIT(IUNIT)
implicit none
INTEGER :: IUNIT
  LSIUNIT = IUNIT
END SUBROUTINE SET_LSIUNIT

END MODULE LSTIMING
