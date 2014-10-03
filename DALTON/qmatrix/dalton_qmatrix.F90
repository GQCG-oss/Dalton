!> \brief QMatrix interface
!> \author Bin Gao
!> \date 2014-06-15
module dalton_qmatrix
    use qmatrix

    implicit none

    ! type of QMatrix interface information
    type, public :: DalQMat
        private
        integer :: stdin = 5
        integer :: stdout = 6
        logical :: do_test = .false.
        logical :: do_scf = .false.
        integer :: num_alpha_occ = 0
        integer, pointer :: alpha_occ_num(:,:)
        integer :: num_beta_occ = 0
        integer, pointer :: beta_occ_num(:,:)
    end type DalQMat

    public :: dalton_qmatrix_init
    public :: dalton_qmatrix_input
    public :: dalton_qmatrix_test
    !public :: dalton_qmatrix_scf
    public :: dalton_qmatrix_finalize

    contains

    !> \brief initializes the QMatrix interface
    !> \author Bin Gao
    !> \date 2014-06-15
    subroutine dalton_qmatrix_init(dal_qmat)
        use matrix_operations
        type(DalQMat), intent(inout) :: dal_qmat
        dal_qmat%stdin = 5
        dal_qmat%stdout = 6
        dal_qmat%do_test = .false.
        dal_qmat%do_scf = .false.
        dal_qmat%num_alpha_occ = 0
        nullify(dal_qmat%alpha_occ_num)
        dal_qmat%num_beta_occ = 0
        nullify(dal_qmat%beta_occ_num)
        return
    end subroutine dalton_qmatrix_init

    !> \brief processes input
    !> \author Bin Gao
    !> \date 2014-06-15
    subroutine dalton_qmatrix_input(dal_qmat, io_input, io_output)
        type(DalQMat), intent(inout) :: dal_qmat
        integer, intent(in) :: io_input
        integer, intent(in) :: io_output
        character(80) word
        character(2) prompt
        logical readword
        dal_qmat%stdin = io_input
        dal_qmat%stdout = io_output
        ! processes input after **QMATRIX
        readword = .true.
        do while (readword)
            read(io_input,100,err=200,end=200) word
            prompt = word(1:2)
            ! end of **QMATRIX
            if (prompt=='**') then
                exit
            ! comments
            else if ((prompt(1:1)=='!') .or. (prompt(1:1)=='#')) then
                cycle
            else if (prompt(1:1)=='*') then
                select case (trim(word))
                case ('*SCF')
                    dal_qmat%do_scf = .true.
                    do
                        read(io_input,100,err=200,end=200) word
                        if (word(1:2)=='**' .or. word=='*END OF INPUT') then
                            readword = .false.
                            exit
                        else if (word(1:1)=='!' .or. word(1:1)=='#') then
                            cycle
                        else
                            select case (trim(word))
                            case ('.ALPHA')
                                read(io_input,*) dal_qmat%num_alpha_occ
                                if (dal_qmat%num_alpha_occ>0) then
                                    allocate(dal_qmat%alpha_occ_num(2,dal_qmat%num_alpha_occ))
                                    read(io_input,*) dal_qmat%alpha_occ_num
                                end if
                            case ('.BETA')
                                read(io_input,*) dal_qmat%num_beta_occ
                                if (dal_qmat%num_beta_occ>0) then
                                    allocate(dal_qmat%beta_occ_num(2,dal_qmat%num_beta_occ))
                                    read(io_input,*) dal_qmat%beta_occ_num
                                end if
                            case default
                                write(io_output,100) ' Keyword "'//trim(word)// &
                                                     '" not recognized in **QMATRIX.'
                                call QUIT('Illegal keyword in **QMATRIX.')
                            end select
                        end if
                    end do
                case default
                    write(io_output,100) ' Keyword "'//trim(word)// &
                                         '" not recognized in **QMATRIX.'
                    call QUIT('Illegal keyword in **QMATRIX.')
                end select
            else if (prompt(1:1)=='.') then
                select case (trim(word))
                case ('.TEST')
                    dal_qmat%do_test = .true.
                case default
                    write(io_output,100) ' Keyword "'//trim(word)// &
                                         '" not recognized in **QMATRIX.'
                    call QUIT('Illegal keyword in **QMATRIX.')
                end select
            end if
            if (word=='*END OF INPUT') exit
        end do
200     return
100     format(A)
    end subroutine dalton_qmatrix_input

    !> \brief performs the Fortran test of the QMatrix library
    !> \author Bin Gao
    !> \date 2014-06-15
    subroutine dalton_qmatrix_test(dal_qmat)
        type(DalQMat), intent(inout) :: dal_qmat
        if (dal_qmat%do_test) then
            write(dal_qmat%stdout,100) ""
            write(dal_qmat%stdout,100) "Fortran test of the QMatrix library"
            write(dal_qmat%stdout,100) "==================================="
            write(dal_qmat%stdout,100) ""
            call test_f_QMatrix(dal_qmat%stdout)
            write(dal_qmat%stdout,100) ""
            write(dal_qmat%stdout,100) "End of the test of the QMatrix library"
            write(dal_qmat%stdout,100) "======================================"
            write(dal_qmat%stdout,100) ""
        end if
        return
100     format(A)
    end subroutine dalton_qmatrix_test

    !> \brief finalizes the QMatrix interface
    !> \author Bin Gao
    !> \date 2014-06-15
    subroutine dalton_qmatrix_finalize(dal_qmat)
        type(DalQMat), intent(inout) :: dal_qmat
        dal_qmat%stdin = 5
        dal_qmat%stdout = 6
        dal_qmat%do_test = .false.
        dal_qmat%do_scf = .false.
        dal_qmat%num_alpha_occ = 0
        if (associated(dal_qmat%alpha_occ_num)) then
            deallocate(dal_qmat%alpha_occ_num)
            nullify(dal_qmat%alpha_occ_num)
        end if
        dal_qmat%num_beta_occ = 0
        if (associated(dal_qmat%beta_occ_num)) then
            deallocate(dal_qmat%beta_occ_num)
            nullify(dal_qmat%beta_occ_num)
        end if
        return
    end subroutine dalton_qmatrix_finalize

end module dalton_qmatrix
