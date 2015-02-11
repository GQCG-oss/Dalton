MODULE KS_settings
use precision
use matrix_module
use matrix_operations
use LS_UTIL, only : capitalize_string
! This module evaluates the Fock/Kohn-Sham matrix using chosen
! algorithm, matrix representation, etc.
!
SAVE
  LOGICAL :: incremental_scheme
  LOGICAL :: do_increment
  LOGICAL :: SaveF0andD0
  TYPE(Matrix) :: incrD0(1), incrF0(1), incrDdiff(1)
CONTAINS

subroutine ks_free_incremental_fock
implicit none
  incremental_scheme = .false.
  SaveF0andD0 = .false.
  do_increment       = .false.
  call mat_free(incrD0(1))
  call mat_free(incrF0(1))
  call mat_free(incrDdiff(1))
end subroutine ks_free_incremental_fock

subroutine ks_init_incremental_fock(nbast)
implicit none
integer :: nbast
  incremental_scheme = .true.
  SaveF0andD0 = .true.
  do_increment       = .false.
  call mat_init(incrD0(1),nbast,nbast)
  call mat_zero(incrD0(1))
  call mat_init(incrF0(1),nbast,nbast)
  call mat_init(incrDdiff(1),nbast,nbast)
end subroutine ks_init_incremental_fock

subroutine ks_init_linesearch_fock(nbast)
implicit none
integer :: nbast
SaveF0andD0 = .true.
IF(.NOT.incremental_scheme)THEN
   call mat_init(incrD0(1),nbast,nbast)
   call mat_zero(incrD0(1))
   call mat_init(incrF0(1),nbast,nbast)
   call mat_init(incrDdiff(1),nbast,nbast)
ENDIF
end subroutine ks_init_linesearch_fock

subroutine ks_free_linesearch_fock()
implicit none
integer :: nbast
  SaveF0andD0 = .false.
  IF(.NOT.incremental_scheme)THEN
     call mat_free(incrD0(1))
     call mat_free(incrF0(1))
     call mat_free(incrDdiff(1))
  ENDIF
end subroutine ks_free_linesearch_fock

subroutine activate_incremental(lupri,do_inc)
implicit none
integer,intent(in) :: lupri
logical :: do_inc
real(realk) :: maxelm

IF(.NOT.do_increment)THEN
   !The incremental scheme is not yet activated. 
   !We determine if we should activate the incremantal scheme
   !this only have an effect on timings in the local region (small differences)
   !if we activate it in the global region there is an unnecessary 
   !accumualtion of error which means that we need to calc 
   !integral with an higher accuracy.  
!   call mat_daxpy(-1E0_realk,D0(1),Ddiff(1)) done outside
   call mat_abs_max_elm(incrDdiff(1),maxelm) 
   IF(maxelm.LT. 0.10E0_realk)THEN
      print*,'activate incremental'
      WRITE(lupri,'(A,F18.12)')'The maximum difference in Density matrix elements are:',maxelm
      WRITE(lupri,'(A)')'The Incremental scheme have therefore been activated.'
      do_increment = .true. 
   ENDIF
ENDIF
do_inc = do_increment
end subroutine activate_incremental

! check if for the functional used the empirical disp. corr. is defined
SUBROUTINE II_dftdispcheck(func,lupri)
  implicit none
  Character(len=80),intent(IN) :: func
  Integer,intent(IN)           :: lupri
  !
  Logical :: disp

  disp = .FALSE.
  
  call capitalize_string(func)
  IF (func == "BP86")  disp = .TRUE.
  IF (func == "BLYP")  disp = .TRUE.
  IF (func == "PBE")   disp = .TRUE.
  IF (func == "B3LYP") disp = .TRUE.
!  IF (insensitveEQUIV(func,"TPSS"))  disp = .TRUE. !!!TPSS not yet implemented in LSDALTON

IF (disp .EQV. .FALSE.) THEN
  write(lupri,'(A)') "###########################################################################################################"
  write(lupri,'(A)') "The empirical dispersion correction is only defined for the following functionals:"
  write(lupri,'(A)') "   BP86"
  write(lupri,'(A)') "   BLYP"
  write(lupri,'(A)') "   PBE"
  write(lupri,'(A)') "   B3LYP"
!  write(lupri,'(A)') "   TPSS"
  write(lupri,'(A)') ""
  write(lupri,'(A)') "At the moment these functionals have to be specified in the LSDALTON.INP with these names explicitly,"
  write(lupri,'(A)') "defining the functionals by specifying the exchange and correlation parts explicitly is not yet possible."
!  write(lupri,'(A)') "TPSS is not implemented yet."
  write(lupri,'(A)') ""
  write(lupri,'(A)') "###########################################################################################################"
  call lsquit("Empirical dispersion correction not defined for the choosen functional",lupri)
ENDIF
END SUBROUTINE II_dftdispcheck

END MODULE KS_settings

SUBROUTINE get_incremental_settings(inc_scheme,do_inc)
use KS_settings
LOGICAL :: inc_scheme,do_inc
inc_scheme = incremental_scheme
do_inc = do_increment
END SUBROUTINE get_incremental_settings
