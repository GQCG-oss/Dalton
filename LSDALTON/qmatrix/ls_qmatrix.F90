!> \brief QMatrix interface
!> \author Bin Gao
!> \date 2014-06-15
module ls_qmatrix
    use precision, only: realk
    use qmatrix_backend
    use qmatrix

    implicit none

    ! type of QMatrix interface information
    type, public :: LSQMat
        private
        integer :: stdin = 5
        integer :: stdout = 6
        logical :: do_test = .false.
        logical :: do_scf = .false.
        integer :: num_alpha_occ = 0
        integer, pointer :: alpha_occ_num(:,:)
        integer :: num_beta_occ = 0
        integer, pointer :: beta_occ_num(:,:)
    end type LSQMat

    public :: ls_qmatrix_init
    public :: ls_qmatrix_input
    public :: ls_qmatrix_test
    !public :: ls_qmatrix_scf
    public :: ls_qmatrix_finalize

    contains

    !> \brief initializes the QMatrix interface
    !> \author Bin Gao
    !> \date 2014-06-15
    subroutine ls_qmatrix_init(ls_qmat)
        type(LSQMat), intent(inout) :: ls_qmat
        ls_qmat%stdin = 5
        ls_qmat%stdout = 6
        ls_qmat%do_test = .false.
        ls_qmat%do_scf = .false.
        ls_qmat%num_alpha_occ = 0
        nullify(ls_qmat%alpha_occ_num)
        ls_qmat%num_beta_occ = 0
        nullify(ls_qmat%beta_occ_num)
        return
    end subroutine ls_qmatrix_init

    !> \brief processes input
    !> \author Bin Gao
    !> \date 2014-06-15
    subroutine ls_qmatrix_input(ls_qmat, io_input, io_output, readword, word)
        type(LSQMat), intent(inout) :: ls_qmat
        integer, intent(in) :: io_input
        integer, intent(in) :: io_output
        logical, intent(inout) :: readword
        character(80), intent(inout) :: word
        character(2) prompt
        ls_qmat%stdin = io_input
        ls_qmat%stdout = io_output
        do while (readword)
            read(io_input,100,err=200,end=200) word
            prompt = word(1:2)
            ! end of **QMATRIX
            if (prompt=='**') then
                readword = .false.
                exit
            ! comments
            else if ((prompt(1:1)=='!') .or. (prompt(1:1)=='#')) then
                cycle
            else if (prompt(1:1)=='*') then
                select case (trim(word))
                case ('*SCF')
                    ls_qmat%do_scf = .true.
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
                                read(io_input,*) ls_qmat%num_alpha_occ
                                if (ls_qmat%num_alpha_occ>0) then
                                    allocate(ls_qmat%alpha_occ_num(2,ls_qmat%num_alpha_occ))
                                    read(io_input,*) ls_qmat%alpha_occ_num
                                end if
                            case ('.BETA')
                                read(io_input,*) ls_qmat%num_beta_occ
                                if (ls_qmat%num_beta_occ>0) then
                                    allocate(ls_qmat%beta_occ_num(2,ls_qmat%num_beta_occ))
                                    read(io_input,*) ls_qmat%beta_occ_num
                                end if
                            case default
                                write(io_output,100) ' Keyword "'//trim(word)// &
                                                     '" not recognized in **QMATRIX.'
                                call lsquit('Illegal keyword in **QMATRIX.', io_output)
                            end select
                        end if
                    end do
                case default
                    write(io_output,100) ' Keyword "'//trim(word)// &
                                         '" not recognized in **QMATRIX.'
                    call lsquit('Illegal keyword in **QMATRIX.', io_output)
                end select
            else if (prompt(1:1)=='.') then
                select case (trim(word))
                case ('.TEST')
                    ls_qmat%do_test = .true.
                case default
                    write(io_output,100) ' Keyword "'//trim(word)// &
                                         '" not recognized in **QMATRIX.'
                    call lsquit('Illegal keyword in **QMATRIX.', io_output)
                end select
            end if
            if (word=='*END OF INPUT') then
                readword = .false.
                exit
            end if
        end do
        return
200     readword = .false.
        return
100     format(A)
    end subroutine ls_qmatrix_input

    !> \brief performs the Fortran test of the QMatrix library
    !> \author Bin Gao
    !> \date 2014-06-15
    subroutine ls_qmatrix_test(ls_qmat)
        type(LSQMat), intent(inout) :: ls_qmat
        if (ls_qmat%do_test) then
            write(ls_qmat%stdout,100) ""
            write(ls_qmat%stdout,100) "Fortran test of the QMatrix library"
            write(ls_qmat%stdout,100) "==================================="
            write(ls_qmat%stdout,100) ""
            call test_f_QMatrix(ls_qmat%stdout)
            write(ls_qmat%stdout,100) ""
            write(ls_qmat%stdout,100) "End of the test of the QMatrix library"
            write(ls_qmat%stdout,100) "======================================"
            write(ls_qmat%stdout,100) ""
        end if
        return
100     format(A)
    end subroutine ls_qmatrix_test

    !> \brief finalizes the QMatrix interface
    !> \author Bin Gao
    !> \date 2014-06-15
    subroutine ls_qmatrix_finalize(ls_qmat)
        type(LSQMat), intent(inout) :: ls_qmat
        ls_qmat%stdin = 5
        ls_qmat%stdout = 6
        ls_qmat%do_test = .false.
        ls_qmat%do_scf = .false.
        ls_qmat%num_alpha_occ = 0
        if (associated(ls_qmat%alpha_occ_num)) then
            deallocate(ls_qmat%alpha_occ_num)
            nullify(ls_qmat%alpha_occ_num)
        end if
        ls_qmat%num_beta_occ = 0
        if (associated(ls_qmat%beta_occ_num)) then
            deallocate(ls_qmat%beta_occ_num)
            nullify(ls_qmat%beta_occ_num)
        end if
        return
    end subroutine ls_qmatrix_finalize

end module ls_qmatrix
