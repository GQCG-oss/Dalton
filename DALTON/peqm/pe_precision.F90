module pe_precision

!    integer, parameter :: sp = selected_real_kind(6, 37)
!    integer, parameter :: dp = selected_real_kind(15, 307)
!    integer, parameter :: qp = selected_real_kind(33, 4931)

    integer, parameter :: sp = kind(1.0)
    integer, parameter :: dp = selected_real_kind(2*precision(1.0_sp))
    integer, parameter :: qp = selected_real_kind(2*precision(1.0_dp))

!    FORTRAN 2008
!    use, intrinsic :: iso_fortran_env
!    integer, parameter :: sp = REAL32
!    integer, parameter :: dp = REAL64
!    integer, parameter :: qp = REAL128

end module pe_precision
