module pcm_write
   
   implicit none
      
   character(11), parameter, public :: pcm_file_name = 'PCM_mep_asc'
   integer,       parameter, public :: pcm_file_unit = 800 
   logical                 , public :: pcm_file_exists = .false.
   logical                 , public :: pcm_file_open = .false.

   public report_after_pcm_input
   public pcm_write_file
   public pcm_write_file_separate
   public host_writer

   private
   
   integer                          :: global_print_unit = -1
  
   contains
   
   subroutine report_after_pcm_input(print_unit)
!
! A subroutine to print the relevant settings of both DALTON and PCMSolver.
!
      use pcmmod_cfg

      integer, optional, intent(in) :: print_unit

      if (present(print_unit)) then
              global_print_unit = print_unit
              write(global_print_unit, *) ' ===== Polarizable Continuum Model calculation set-up ====='
              write(global_print_unit, *) '* Polarizable Continuum Model using PCMSolver external module:'                              
              write(global_print_unit, *) '  1: Converged potentials and charges at tesserae representative points written on file.'

              if (pcmmod_old_integration) then
                 write(global_print_unit, *) '  2: Use point-by-point charge-attraction integrals evaluation subroutines.'
              else
                 write(global_print_unit, *) '  2: Use vectorized charge-attraction integrals evaluation subroutines.'
              end if
              
              if (pcmmod_separate) then
                 write(global_print_unit, *) '  3: Separate potentials and apparent charges in nuclear and electronic.'
              else
                 write(global_print_unit, *) '  3: Use total potentials and apparent charges.'
              end if
              
              if (pcmmod_print > 5 .and. pcmmod_print < 10) then
                 write(global_print_unit, *) '  4: Print potentials at tesserae representative points.'
              else if (pcmmod_print > 10) then
                 write(global_print_unit, *) '  4: Print potentials and charges at tesserae representative points.'
              else
                 write(global_print_unit, *) '  4: Do not print potentials and charges.'
              end if     
      ! Should print info for the cavity and the solver here.
              call print_citation
      else
              write(*, *) ' ===== Polarizable Continuum Model calculation set-up ====='                                         
              write(*, *) '* Polarizable Continuum Model using PCMSolver external module:'
              write(*, *) '  1: Converged potentials and charges at tesserae representative points written on file.'
              
              if (pcmmod_separate) then
                 write(*, *) '  2: Separate potentials and apparent charges in nuclear and electronic.'
              else
                 write(*, *) '  2: Use total potentials and apparent charges.'
              end if
              
              if (pcmmod_print > 5 .and. pcmmod_print < 10) then
                 write(*, *) '  3: Print potentials at tesserae representative points.'
              else if (pcmmod_print > 10) then
                 write(*, *) '  3: Print potentials and charges at tesserae representative points.'
              else
                 write(*, *) '  3: Do not print potentials and charges.'
              end if     
      ! Should print info for the cavity and the solver here.
              call print_citation
      end if

   end subroutine

   subroutine pcm_write_file(nr_points, potentials, charges)

! Passed variables
      integer(4), intent(in) :: nr_points
      real(8), intent(in) :: potentials(nr_points)
      real(8), intent(in) :: charges(nr_points)
 
! Local variables
      real(8)       :: tot_chg
      integer       :: ipoint
      character(8)  :: for_title  = '(20X, A)'
      character(19) :: for_header = '(A, T27, A, T62, A)'
      character(20) :: for_data   = '(I6, 2(20X, F15.12))'

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
      write(pcm_file_unit, '(A, F15.12)') 'Sum of apparent surface charges ', tot_chg
    
   end subroutine 
   
   subroutine pcm_write_file_separate(nr_points, nuc_pot, nuc_chg, ele_pot, ele_chg)

! Passed variables
      integer(4), intent(in) :: nr_points
      real(8), intent(in) :: nuc_pot(nr_points), ele_pot(nr_points)
      real(8), intent(in) :: nuc_chg(nr_points), ele_chg(nr_points)
 
! Local variables
      real(8) :: tot_nuc_chg, tot_ele_chg
      integer :: ipoint
      character(8)  :: for_title  = '(60X, A)'
      character(36) :: for_header = '(A, T27, A, T62, A, T97, A, T132, A)'
      character(20) :: for_data   = '(I6, 4(20X, F15.12))'

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
      write(pcm_file_unit, '(A, F15.12)') 'Sum of nuclear apparent charges ', tot_nuc_chg
      write(pcm_file_unit, '(A, F15.12)') 'Sum of electronic apparent charges ', tot_ele_chg
      write(pcm_file_unit, '(A, F15.12)') 'Sum of apparent surface charges ', tot_nuc_chg + tot_ele_chg
    
   end subroutine 
   
   subroutine host_writer(message, message_length) bind(c, name='host_writer')
       
      use, intrinsic :: iso_c_binding

      integer(c_size_t), intent(in) :: message_length
      character(kind=c_char)        :: message(message_length)

      if (global_print_unit .eq. -1) then
        write(*, '(1000A)') message
      else
        write(global_print_unit, '(1000A)') message
      end if
      flush(global_print_unit)

   end subroutine host_writer
   
end module
