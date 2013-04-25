MODULE KS_settings
use matrix_module
use matrix_operations
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

END MODULE KS_settings

SUBROUTINE get_incremental_settings(inc_scheme,do_inc)
use KS_settings
LOGICAL :: inc_scheme,do_inc
inc_scheme = incremental_scheme
do_inc = do_increment
END SUBROUTINE get_incremental_settings
