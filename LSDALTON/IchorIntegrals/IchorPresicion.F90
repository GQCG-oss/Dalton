!> @file
!> Contains definition of realk - normally double precision
!> \brief Contains definition of realk - normally double precision
!> \author T. Kjaergaard
!> \date 2013 
MODULE IchorPrecisionModule
#ifdef SYS_REAL
  INTEGER, PARAMETER :: realk = 4
#else
  INTEGER, PARAMETER :: realk = 8
#endif
! Long integer is equivalent to 64 bit integers (8 bytes) and ranges 
! from -9223372036854775808 to 9223372036854775807 (i.e., -2**63 to 2**63-1)
  integer, parameter :: long = 8
!integer defining the size of the MPI input (depending on wheter the mpi library
!is compiled with 64 or 32 bit integers)
contains
!Added to avoid "has no symbols" linking warning
subroutine Ichorprecision_void()
end subroutine Ichorprecision_void

END MODULE IchorPrecisionModule
