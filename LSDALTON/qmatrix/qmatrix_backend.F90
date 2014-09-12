!> \brief a backend matrix module for QMatrix library, which will invoke
!>        the matrix module of LSDALTON
!> \author Bin Gao
!> \date 2014-03-20
module qmatrix_backend
    use precision, only: realk
    use files, only: lsopen, &
                     lsclose
    use Matrix_module, only: Matrix
    use matrix_operations, only: mat_free,           &
                                 mat_init,           &
                                 mat_create_block,   &
                                 mat_add_block,      &
                                 mat_retrieve_block, &
                                 mat_assign,         &
                                 mat_zero,           &
                                 mat_tr,             &
                                 mat_trAB,           &
                                 mat_dotproduct,     &
                                 mat_write_to_disk,  &
                                 mat_print,          &
                                 mat_read_from_disk, &
                                 mat_scal,           &
                                 mat_daxpy,          &
                                 mat_trans,          &
                                 mat_mul

    implicit none

    ! type of real matrix
    type, public :: real_mat_t
        private
        integer :: nrow = -1             !number of rows
        integer :: ncol = -1             !number of columns
        integer :: sym_type = 0          !symmetry of the matrix, not used
        logical :: assembled = .false.   !indicates if the matrix is assembled or not
        type(Matrix), pointer :: ls_mat  !matrix type implemented in LSDALTON
    end type real_mat_t

    public :: Matrix_Create
    public :: Matrix_Destroy
    public :: Matrix_SetSymType
    public :: Matrix_SetDimMat
    public :: Matrix_Assemble
    public :: Matrix_GetSymType
    public :: Matrix_GetDimMat
    public :: Matrix_IsAssembled
    public :: Matrix_SetValues
    public :: Matrix_AddValues
    public :: Matrix_GetValues
    public :: Matrix_Duplicate
    public :: Matrix_ZeroEntries
    public :: Matrix_GetTrace
    public :: Matrix_GetMatProdTrace
    public :: Matrix_Write
    public :: Matrix_Read
    public :: Matrix_Scale
    public :: Matrix_AXPY
    public :: Matrix_Transpose
    public :: Matrix_GEMM

    contains

    ! all matrices should be created by this subroutine before their use
    subroutine Matrix_Create(A)
        type(real_mat_t), intent(inout) :: A
        A%nrow = -1
        A%ncol = -1
        A%sym_type = 0
        A%assembled = .false.
        allocate(A%ls_mat)
        return
    end subroutine Matrix_Create

    ! all matrices should be destroyed by this subroutine after their use
    subroutine Matrix_Destroy(A)
        type(real_mat_t), intent(inout) :: A
        if (A%assembled) then
            call mat_free(A%ls_mat)
        end if
        A%nrow = -1
        A%ncol = -1
        A%sym_type = 0
        A%assembled = .false.
        deallocate(A%ls_mat)
        nullify(A%ls_mat)
        return
    end subroutine Matrix_Destroy

    subroutine Matrix_SetSymType(A, sym_type)
        type(real_mat_t), intent(inout) :: A
        integer, intent(in) :: sym_type
        A%sym_type = sym_type
        return
    end subroutine Matrix_SetSymType

    subroutine Matrix_SetDimMat(A, dim_mat)
        type(real_mat_t), intent(inout) :: A
        integer, intent(in) :: dim_mat
        A%nrow = dim_mat
        A%ncol = dim_mat
        return
    end subroutine Matrix_SetDimMat

    ! in general, all matrices should be assembled before used in the matrix operations
    subroutine Matrix_Assemble(A)
        type(real_mat_t), intent(inout) :: A
        if (A%nrow<1) then
            call lsquit("Matrix_Assemble>> invalid number of rows", A%nrow)
        end if
        if (A%ncol<1) then
            call lsquit("Matrix_Assemble>> invalid number of columns", A%ncol)
        end if
        call mat_init(A%ls_mat, A%nrow, A%ncol)
        A%assembled = .true.
        return
    end subroutine Matrix_Assemble

    subroutine Matrix_GetSymType(A, sym_type)
        type(real_mat_t), intent(in) :: A
        integer, intent(out) :: sym_type
        sym_type = A%sym_type
        return
    end subroutine Matrix_GetSymType

    subroutine Matrix_GetDimMat(A, dim_mat)
        type(real_mat_t), intent(in) :: A
        integer, intent(out) :: dim_mat
        dim_mat = A%nrow
        return
    end subroutine Matrix_GetDimMat

    subroutine Matrix_IsAssembled(A, assembled)
        type(real_mat_t), intent(in) :: A
        logical, intent(out) :: assembled
        assembled = A%assembled
        return
    end subroutine Matrix_IsAssembled

    subroutine Matrix_SetValues(A,             &
                                idx_first_row, &
                                num_row_set,   &
                                idx_first_col, &
                                num_col_set,   &
                                values)
        type(real_mat_t), intent(inout) :: A
        integer, intent(in) :: idx_first_row
        integer, intent(in) :: num_row_set
        integer, intent(in) :: idx_first_col
        integer, intent(in) :: num_col_set
        real(realk), intent(in) :: values(num_row_set*num_col_set)
        if (A%assembled) then
            call mat_create_block(A%ls_mat,      &
                                  values,        &
                                  num_row_set,   &
                                  num_col_set,   &
                                  idx_first_row, &
                                  idx_first_col)
        else
            call lsquit("Matrix_SetValues>> A is not assembled", -1)
        end if
        return
    end subroutine Matrix_SetValues

    subroutine Matrix_AddValues(A,             &
                                idx_first_row, &
                                num_row_add,   &
                                idx_first_col, &
                                num_col_add,   &
                                values)
        type(real_mat_t), intent(inout) :: A
        integer, intent(in) :: idx_first_row
        integer, intent(in) :: num_row_add
        integer, intent(in) :: idx_first_col
        integer, intent(in) :: num_col_add
        real(realk), intent(in) :: values(num_row_add*num_col_add)
        if (A%assembled) then
            call mat_add_block(A%ls_mat,      &
                               values,        &
                               num_row_add,   &
                               num_col_add,   &
                               idx_first_row, &
                               idx_first_col)
        else
            call lsquit("Matrix_AddValues>> A is not assembled", -1)
        end if
        return
    end subroutine Matrix_AddValues

    subroutine Matrix_GetValues(A,             &
                                idx_first_row, &
                                num_row_get,   &
                                idx_first_col, &
                                num_col_get,   &
                                values)
        type(real_mat_t), intent(inout) :: A
        integer, intent(in) :: idx_first_row
        integer, intent(in) :: num_row_get
        integer, intent(in) :: idx_first_col
        integer, intent(in) :: num_col_get
        real(realk), intent(out) :: values(num_row_get*num_col_get)
        if (A%assembled) then
            call mat_retrieve_block(A%ls_mat,      &
                                    values,        &
                                    num_row_get,   &
                                    num_col_get,   &
                                    idx_first_row, &
                                    idx_first_col)
        else
            call lsquit("Matrix_GetValues>> A is not assembled", -1)
        end if
        return
    end subroutine Matrix_GetValues

    ! B should be created, but may not be assembled
    subroutine Matrix_Duplicate(A, duplicate_option, B)
        type(real_mat_t), intent(in) :: A
        integer, intent(in) :: duplicate_option
        type(real_mat_t), intent(inout) :: B
#include "api/qmatrix_f_mat_duplicate.h90"
        ! B has been assembled before
        if (B%assembled) then
            call mat_free(B%ls_mat)
        end if
        B%nrow = A%nrow
        B%ncol = A%ncol
        B%sym_type = A%sym_type
        B%assembled = A%assembled
        if (A%assembled) then
            call mat_init(B%ls_mat, B%nrow, B%ncol)
            ! copies the values
            if (duplicate_option==COPY_PATTERN_AND_VALUE) then
                call mat_assign(B%ls_mat, A%ls_mat)
            end if
        end if
        return
    end subroutine Matrix_Duplicate

    subroutine Matrix_ZeroEntries(A)
        type(real_mat_t), intent(inout) :: A
        if (A%assembled) then
            call mat_zero(A%ls_mat)
        else
            call lsquit("Matrix_ZeroEntries>> A is not assembled", -1)
        end if
        return
    end subroutine Matrix_ZeroEntries

    subroutine Matrix_GetTrace(A, trace)
        type(real_mat_t), intent(in) :: A
        real(realk), intent(out) :: trace
        if (A%assembled) then
            trace = mat_tr(A%ls_mat)
        else
            call lsquit("Matrix_GetTrace>> A is not assembled", -1)
        end if
        return
    end subroutine Matrix_GetTrace

    subroutine Matrix_GetMatProdTrace(A, B, op_B, trace)
        type(real_mat_t), intent(in) :: A
        type(real_mat_t), intent(in) :: B
        integer, intent(in) :: op_B
        real(realk), intent(out) :: trace
#include "api/qmatrix_f_mat_operations.h90"
        if (.not.A%assembled) then
            call lsquit("Matrix_GEMM>> A is not assembled", -1)
        end if
        if (.not.B%assembled) then
            call lsquit("Matrix_GEMM>> B is not assembled", -1)
        end if
        select case (op_B)
        case (MAT_NO_OPERATION)
            if (A%ncol/=B%nrow) then
                call lsquit("Matrix_GetMatProdTrace>> invalid dimensions (1)", &
                            A%ncol-B%nrow)
            end if
            if (A%nrow/=B%ncol) then
                call lsquit("Matrix_GetMatProdTrace>> non-square product (1)", &
                            A%nrow-B%ncol)
            end if
            ! trace = \sum_{ij}A_{ij}*B_{ji}
            trace = mat_trAB(A%ls_mat, B%ls_mat)
        case (MAT_TRANSPOSE)
            if (A%ncol/=B%ncol) then
                call lsquit("Matrix_GetMatProdTrace>> invalid dimensions (2)", &
                            A%ncol-B%ncol)
            end if
            if (A%nrow/=B%nrow) then
                call lsquit("Matrix_GetMatProdTrace>> non-square product (2)", &
                            A%nrow-B%nrow)
            end if
            ! trace = \sum_{ij}A_{ij}*B_{ij}
            trace = mat_dotproduct(A%ls_mat, B%ls_mat)
        case default
            call lsquit("Matrix_GetMatProdTrace>> invalid operation on B", op_B)
        end select
        return
    end subroutine Matrix_GetMatProdTrace

    subroutine Matrix_Write(A, mat_label, view_option)
        type(real_mat_t), intent(in) :: A
        character*(*), intent(in) :: mat_label
        integer, intent(in) :: view_option
#include "api/qmatrix_f_mat_view.h90"
        integer io_matrix
        if (A%assembled) then
            select case (view_option)
            ! BINARY_VIEW
            case (BINARY_VIEW)
                io_matrix = -1
                call lsopen(io_matrix,                   &
                            "qmatrix_"//trim(mat_label), &
                            "unknown",                   &
                            "unformatted")
                call mat_write_to_disk(io_matrix, A%ls_mat, .true.)
                call lsclose(io_matrix, "KEEP")
                io_matrix = -1
                call lsopen(io_matrix,                           &
                            "qmatrix_"//trim(mat_label)//"_dim", &
                            "unknown",                           &
                            "unformatted")
                write(io_matrix) A%nrow, A%ncol, A%sym_type
                call lsclose(io_matrix, "KEEP")
            ! ASCII_VIEW
            case (ASCII_VIEW)
                io_matrix = -1
                call lsopen(io_matrix,                   &
                            "qmatrix_"//trim(mat_label), &
                            "unknown",                   &
                            "formatted")
                call mat_print(A%ls_mat, 1, A%nrow, 1, A%ncol, io_matrix)
                call lsclose(io_matrix, "KEEP")
            case default
                call lsquit("Matrix_Write>> invalid view option", view_option)
            end select
        else
            call lsquit("Matrix_Write>> A is not assembled", -1)
        end if
        return
    end subroutine Matrix_Write

    ! A should be created, but may not be assembled
    subroutine Matrix_Read(A, mat_label, view_option)
        type(real_mat_t), intent(inout) :: A
        character*(*), intent(in) :: mat_label
        integer, intent(in) :: view_option
#include "api/qmatrix_f_mat_view.h90"
        integer io_matrix
        ! A has been assembled before
        if (A%assembled) then
            call mat_free(A%ls_mat)
        end if
        select case (view_option)
        ! BINARY_VIEW
        case (BINARY_VIEW)
            io_matrix = -1
            call lsopen(io_matrix,                           &
                        "qmatrix_"//trim(mat_label)//"_dim", &
                        "old",                               &
                        "unformatted")
            read(io_matrix) A%nrow, A%ncol, A%sym_type
            call lsclose(io_matrix, "KEEP")
            call mat_init(A%ls_mat, A%nrow, A%ncol)
            io_matrix = -1
            call lsopen(io_matrix, "qmatrix_"//trim(mat_label), "old", "unformatted")
            call mat_read_from_disk(io_matrix, A%ls_mat, .true.)
            call lsclose(io_matrix, "KEEP")
        ! ASCII_VIEW
        case (ASCII_VIEW)
            call lsquit("Matrix_Read>> ASCII_VIEW not supported", view_option)
        case default
            call lsquit("Matrix_Read>> invalid view option", view_option)
        end select
        A%assembled = .true.
        return
    end subroutine Matrix_Read

    subroutine Matrix_Scale(scal_number, A)
        real(realk), intent(in) :: scal_number
        type(real_mat_t), intent(inout) :: A
        if (A%assembled) then
            call mat_scal(scal_number, A%ls_mat)
        else
            call lsquit("Matrix_Scale>> A is not assembled", -1)
        end if
        return
    end subroutine Matrix_Scale

    subroutine Matrix_AXPY(multiplier, X, Y)
        real(realk), intent(in) :: multiplier
        type(real_mat_t), intent(in) :: X
        type(real_mat_t), intent(inout) :: Y
        if (.not.X%assembled) then
            call lsquit("Matrix_AXPY>> X is not assembled", -1)
        end if
        if (.not.Y%assembled) then
            call lsquit("Matrix_AXPY>> Y is not assembled", -1)
        end if
        call mat_daxpy(multiplier, X%ls_mat, Y%ls_mat)
        Y%sym_type = X%sym_type
        return
    end subroutine Matrix_AXPY

    ! B should be created, but may not be assembled
    subroutine Matrix_Transpose(op_A, A, B)
        integer, intent(in) :: op_A
        type(real_mat_t), intent(inout) :: A
        type(real_mat_t), intent(inout) :: B
#include "api/qmatrix_f_mat_operations.h90"
#include "api/qmatrix_f_mat_duplicate.h90"
        type(real_mat_t) C
        if (.not.A%assembled) then
            call lsquit("Matrix_Transpose>> A is not assembled", -1)
        end if
        select case (op_A)
        case (MAT_NO_OPERATION)
            ! A=A
            if (associated(B%ls_mat, target=A%ls_mat)) then
                return
            ! B=A
            else
                call Matrix_Duplicate(A, COPY_PATTERN_AND_VALUE, B)
            end if
        case (MAT_TRANSPOSE)
            ! in-place transpose, we copy A to another temporary matrix C
            ! and then perform out-of-place transpose
            if (associated(B%ls_mat, target=A%ls_mat)) then
!FIXME: not work for non-square matrix, also costing extra memory here using C
                call Matrix_Create(C)
                call Matrix_Duplicate(A, COPY_PATTERN_AND_VALUE, C)
                call mat_trans(C%ls_mat, B%ls_mat)
                call Matrix_Destroy(C)
            ! out-of-place transpose
            else
                ! B is not the same shape as A^{T}, or B is not assembled before
                if (B%nrow/=A%ncol .or. B%ncol/=A%nrow) then
                    ! B is assembled before
                    if (B%assembled) then
                        call mat_free(B%ls_mat)
                    end if
                    B%nrow = A%ncol
                    B%ncol = A%nrow
                    B%sym_type = A%sym_type
                    B%assembled = .true.
                    call mat_init(B%ls_mat, B%nrow, B%ncol)
                else if (.not.B%assembled) then
                    B%sym_type = A%sym_type
                    B%assembled = .true.
                    call mat_init(B%ls_mat, B%nrow, B%ncol)
                end if
                call mat_trans(A%ls_mat, B%ls_mat)
            end if
        case default
            call lsquit("Matrix_Transpose>> invalid operation on A", op_A)
        end select
        return
    end subroutine Matrix_Transpose

    ! C should be created, but may not be assembled
    subroutine Matrix_GEMM(op_A, op_B, alpha, A, B, beta, C)
        integer, intent(in) :: op_A
        integer, intent(in) :: op_B
        real(realk), intent(in) :: alpha
        type(real_mat_t), intent(in) :: A
        type(real_mat_t), intent(in) :: B
        real(realk), intent(in) :: beta
        type(real_mat_t), intent(inout) :: C
#include "api/qmatrix_f_mat_operations.h90"
        character transa, transb
        if (.not.A%assembled) then
            call lsquit("Matrix_GEMM>> A is not assembled", -1)
        end if
        if (.not.B%assembled) then
            call lsquit("Matrix_GEMM>> B is not assembled", -1)
        end if
        select case (op_A)
        case (MAT_NO_OPERATION)
            transa = "N"
        case (MAT_TRANSPOSE)
            transa = "T"
        case default
            call lsquit("Matrix_GEMM>> invalid operation on A", op_A)
        end select
        select case (op_B)
        case (MAT_NO_OPERATION)
            transb = "N"
        case (MAT_TRANSPOSE)
            transb = "T"
        case default
            call lsquit("Matrix_GEMM>> invalid operation on B", op_B)
        end select
        ! C is not assembled
        if (.not.C%assembled) then
            C%nrow = A%nrow
            C%ncol = B%ncol
            C%assembled = .true.
            call mat_init(C%ls_mat, C%nrow, C%ncol)
        else
            if (C%nrow/=A%nrow) then
                call lsquit("Matrix_GEMM>> invalid number of rows", C%nrow)
            end if
            if (C%ncol/=B%ncol) then
                call lsquit("Matrix_GEMM>> invalid number of columns", C%ncol)
            end if
        end if
        call mat_mul(A%ls_mat, B%ls_mat, transa, transb, alpha, beta, C%ls_mat)
        C%sym_type = 0
        return
    end subroutine Matrix_GEMM

end module qmatrix_backend
