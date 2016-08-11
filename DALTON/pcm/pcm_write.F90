!> @file
!> write functionality for PCM
module pcm_write

implicit none

public init_host_writer
public host_writer
public pcm_write_file
public pcm_write_file_separate

private

character(11), parameter :: pcm_file_name = 'PCM_mep_asc'
integer,       parameter :: pcm_file_unit = 800
logical                  :: pcm_file_exists = .false.
logical                  :: pcm_file_open = .false.
integer                  :: global_print_unit = -1

contains

!> \brief initializes the global_print_unit for host_writer
!> \author R. Di Remigio
!> \date 2014
!> \param print_unit DALTON output print unit
subroutine init_host_writer(print_unit)

   integer :: print_unit

   if (global_print_unit .eq. -1) then
     global_print_unit = print_unit
   else
     write(*, *) "Printing to stdout"
   end if

end subroutine init_host_writer

!> \brief prints message to LSDALTON output
!> \author R. Di Remigio
!> \date 2014
!> \param message the message to be printed out
!>
!> This subroutine is called by PCMSolver to print its own output to
!> LSDALTON output.
subroutine host_writer(message, message_length) bind(c, name='host_writer')

    use, intrinsic :: iso_c_binding, only: c_char, c_size_t

    character(kind=c_char)        :: message(*)
    integer(c_size_t), intent(in), value :: message_length

    character(len=message_length) :: f_message

    call pcmsolver_c2f_string(message, f_message, message_length)
    if (global_print_unit .eq. -1) then
        write(*, '(1000A)') f_message
    else
        write(global_print_unit, '(1000A)') f_message
    end if

end subroutine host_writer

!> \brief print converged charges and potentials to file
!> \author R. Di Remigio
!> \date 2014
!> \param nr_points number of cavity points
!> \param potentials converged potential at cavity points
!> \param charges converged charge at cavity points
subroutine pcm_write_file(nr_points, potentials, charges)

! Passed variables
   integer(8), intent(in) :: nr_points
   real(8), intent(in) :: potentials(nr_points)
   real(8), intent(in) :: charges(nr_points)

! Local variables
   real(8)       :: tot_chg
   integer(4)       :: ipoint
   character(8)  :: for_title  = '(20X, A)'
   character(19) :: for_header = '(A, T27, A, T62, A)'
   character(20) :: for_data   = '(I6, 2(20X, F20.12))'

   inquire(file = pcm_file_name, exist = pcm_file_exists)
   if (pcm_file_exists) then
      open(pcm_file_unit, &
           file = pcm_file_name, &
           status = 'old', &
           form = 'formatted', &
           access = 'sequential')
      close(pcm_file_unit, status = 'delete')
   end if
   open(pcm_file_unit, &
        file = pcm_file_name, &
        status = 'new', &
        form = 'formatted', &
        access = 'sequential')
   rewind(pcm_file_unit)

   write(pcm_file_unit, for_title) "Converged MEP and ASC"
   write(pcm_file_unit, for_header) "Finite element #", "Total MEP", "Total ASC"
   tot_chg = 0.0d0
   do ipoint = 1, nr_points
      tot_chg = tot_chg + charges(ipoint)
      write(pcm_file_unit, for_data) ipoint, potentials(ipoint), charges(ipoint)
   end do
   write(pcm_file_unit, '(A, F20.12)') 'Sum of apparent surface charges ', tot_chg

end subroutine pcm_write_file

!> \brief print converged charges and potentials to file, separate case
!> \author R. Di Remigio
!> \date 2014
!> \param nr_points number of cavity points
!> \param nuc_pot nuclear potential at cavity points
!> \param nuc_chg nuclear charge at cavity points
!> \param ele_pot converged electronic potential at cavity points
!> \param ele_chg converged electronic charge at cavity points
subroutine pcm_write_file_separate(nr_points, nuc_pot, nuc_chg, ele_pot, ele_chg)

! Passed variables
   integer(8), intent(in) :: nr_points
   real(8), intent(in) :: nuc_pot(nr_points), ele_pot(nr_points)
   real(8), intent(in) :: nuc_chg(nr_points), ele_chg(nr_points)

! Local variables
   real(8) :: tot_nuc_chg, tot_ele_chg
   integer(4) :: ipoint
   character(8)  :: for_title  = '(60X, A)'
   character(36) :: for_header = '(A, T27, A, T62, A, T97, A, T132, A)'
   character(20) :: for_data   = '(I6, 4(20X, F20.12))'

   inquire(file = pcm_file_name, exist = pcm_file_exists)
   if (pcm_file_exists) then
      open(pcm_file_unit, &
           file = pcm_file_name, &
           status = 'old', &
           form = 'formatted', &
           access = 'sequential')
      close(pcm_file_unit, status = 'delete')
   end if
   open(pcm_file_unit, &
        file = pcm_file_name, &
        status = 'new', &
        form = 'formatted', &
        access = 'sequential')
   rewind(pcm_file_unit)

   write(pcm_file_unit, for_title) "Converged MEP and ASC"
   write(pcm_file_unit, for_header) "Finite element #", "Nuclear MEP", "Nuclear ASC", "Electronic MEP", "Electronic ASC"
   tot_nuc_chg = 0.0d0
   tot_ele_chg = 0.0d0
   do ipoint = 1, nr_points
      tot_nuc_chg = tot_nuc_chg + nuc_chg(ipoint)
      tot_ele_chg = tot_ele_chg + ele_chg(ipoint)
      write(pcm_file_unit, for_data) ipoint, nuc_pot(ipoint), nuc_chg(ipoint), ele_pot(ipoint), ele_chg(ipoint)
   end do
   write(pcm_file_unit, '(A, F20.12)') 'Sum of nuclear apparent charges ', tot_nuc_chg
   write(pcm_file_unit, '(A, F20.12)') 'Sum of electronic apparent charges ', tot_ele_chg
   write(pcm_file_unit, '(A, F20.12)') 'Sum of apparent surface charges ', tot_nuc_chg + tot_ele_chg

end subroutine pcm_write_file_separate

end module
