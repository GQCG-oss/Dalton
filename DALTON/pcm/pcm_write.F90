module pcm_write
   
   implicit none
      
   character(*), parameter, public :: pcm_file_name = 'PCM_mep_asc'
   integer,      parameter, public :: pcm_file_unit = 800 
   logical                , public :: pcm_file_exists = .false.
   logical                , public :: pcm_file_open = .false.

   public report_after_pcm_input
   public pcm_write_file
   public pcm_write_file_separate

   private
  
   contains
   
   subroutine report_after_pcm_input(print_unit)
!
! A subroutine to print the relevant settings of both DALTON and PCMSolver.
!
      use pcmmod_cfg

      integer, optional :: print_unit
      integer           :: print_here

      if (present(print_unit)) then
              print_here = print_unit
              write(print_here, *) ' ===== Polarizable Continuum Model calculation set-up ====='
              write(print_here, *) '* Polarizable Continuum Model using PCMSolver external module:'                              
              write(print_here, *) '  1: Converged potentials and charges at tesserae representative points written on file.'
              
              if (pcmmod_separate) then
                 write(print_here, *) '  2: Separate potentials and apparent charges in nuclear and electronic.'
              else
                 write(print_here, *) '  2: Use total potentials and apparent charges.'
              end if
              
              if (pcmmod_print > 5 .and. pcmmod_print < 10) then
                 write(print_here, *) '  3: Print potentials at tesserae representative points.'
              else if (pcmmod_print > 10) then
                 write(print_here, *) '  3: Print potentials and charges at tesserae representative points.'
              else
                 write(print_here, *) '  3: Do not print potentials and charges.'
              end if     
      ! Should print info for the cavity and the solver here.
              write(print_here, *) '* References: '                                                               
              write(print_here, *) '  - J. Tomasi, B. Mennucci and R. Cammi:',                &
              '   "Quantum Mechanical Continuum Solvation Models", Chem. Rev., 105 (2005) 2999' 
              write(print_here, *) '  - PCMSolver, an API for the Polarizable Continuum Model electrostatic problem',   &    
              '    L. Frediani et al. '
      else
              print_here = 0
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
              write(*, *) '* References: '                                                               
              write(*, *) '  - J. Tomasi, B. Mennucci and R. Cammi:',                &
              '   "Quantum Mechanical Continuum Solvation Models", Chem. Rev., 105 (2005) 2999' 
              write(*, *) '  - PCMSolver, an API for the Polarizable Continuum Model electrostatic problem',   &    
              '    L. Frediani et al. '
      end if

   end subroutine

   subroutine pcm_write_file(nr_points, potentials, charges)

! Passed variables
      integer, intent(in) :: nr_points
      real(8), intent(in) :: potentials(nr_points)
      real(8), intent(in) :: charges(nr_points)
 
! Local variables
      real(8) :: tot_chg
      integer :: ipoint

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

      tot_chg = 0.0d0
      do ipoint = 1, nr_points
         tot_chg = tot_chg + charges(ipoint)
         write(pcm_file_unit, *) '@point', ipoint, "MEP =", potentials(ipoint), "ASC =", charges(ipoint)
      end do
      write(pcm_file_unit, *) 'Number of tesserae', nr_points, 'Sum of apparent charges', tot_chg
    
   end subroutine 
   
   subroutine pcm_write_file_separate(nr_points, nuc_pot, nuc_chg, ele_pot, ele_chg)

! Passed variables
      integer, intent(in) :: nr_points
      real(8), intent(in) :: nuc_pot(nr_points), ele_pot(nr_points)
      real(8), intent(in) :: nuc_chg(nr_points), ele_chg(nr_points)
 
! Local variables
      real(8) :: tot_nuc_chg, tot_ele_chg
      integer :: ipoint

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

      tot_nuc_chg = 0.0d0
      tot_ele_chg = 0.0d0
      do ipoint = 1, nr_points
         tot_nuc_chg = tot_nuc_chg + nuc_chg(ipoint)
         tot_ele_chg = tot_ele_chg + ele_chg(ipoint)
         write(pcm_file_unit, *) '@point', ipoint
         write(pcm_file_unit, *) "NUC_MEP =", nuc_pot(ipoint), "NUC_ASC =", nuc_chg(ipoint)
         write(pcm_file_unit, *) "ELE_MEP =", ele_pot(ipoint), "ELE_ASC =", ele_chg(ipoint)
      end do
      write(pcm_file_unit, *) 'Number of tesserae', nr_points 
      write(pcm_file_unit, *) 'Sum of nuclear apparent charges', tot_nuc_chg, 'Sum of electronic apparent charges', tot_ele_chg
      write(pcm_file_unit, *) 'Sum of apparent charges', tot_nuc_chg + tot_ele_chg
    
   end subroutine 
   
end module
