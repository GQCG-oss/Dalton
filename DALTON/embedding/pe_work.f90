module pe_work

    use pe_precision

    implicit none

    real(dp), dimension(:), pointer :: pewrk

contains

subroutine dalwrk2pewrk(dalwrk)

    real(dp), dimension(:), target :: dalwrk

    pewrk => dalwrk

end subroutine dalwrk2pewrk

subroutine nullify_pewrk()

    nullify(pewrk)

end subroutine nullify_pewrk

end module pe_work
