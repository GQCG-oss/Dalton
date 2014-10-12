!> @file
!> Contains definition of realk - normally double precision
!> \brief Contains definition of realk - normally double precision
!> \author T. Kjaergaard
!> \date 2013 
MODULE IchorPrecisionModule
#ifdef VAR_OPENACC
  !OpenACC libary routines 
  use openacc, only: acc_handle_kind
#endif
! Long integer is equivalent to 64 bit integers (8 bytes) and ranges 
! from -9223372036854775808 to 9223372036854775807 (i.e., -2**63 to 2**63-1)
  integer, parameter :: long = 8
! maximum value of a 32 bit integer. 
  INTEGER, PARAMETER :: MAXINT = 2147483640

  INTEGER, PARAMETER :: reals = 4
  INTEGER, PARAMETER :: reald = 8
#ifdef SYS_REAL
  INTEGER, PARAMETER :: realk = 4
#else
  INTEGER, PARAMETER :: realk = 8
#endif
#ifdef VAR_OPENACC
  !OpenACC libary kind specification
  INTEGER,parameter :: acckind = acc_handle_kind
#else
#ifdef VAR_INT64
  INTEGER,parameter :: acckind = 8
#else
  INTEGER,parameter :: acckind = 4 
#endif
#endif
contains
!Added to avoid "has no symbols" linking warning
subroutine Ichorprecision_void()
end subroutine Ichorprecision_void

END MODULE IchorPrecisionModule
