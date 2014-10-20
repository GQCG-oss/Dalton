!> @file
!> configure PCM input LSDALTON-side 
module ls_pcm_config

use, intrinsic :: iso_c_binding 

implicit none

private

public pcmtype
public ls_pcm_init
public ls_pcm_input
public host_input

type pcmtype
     logical :: do_pcm
     logical :: separate
     integer :: print_level
     logical :: host_provides_input
end type pcmtype

type(pcmtype), public    :: pcm_config
integer                  :: global_print_unit = -1

! *PCMSOL section
! cavity specification *PCMSOL section
character(len=8), public, save :: pcmmod_cavity_type = 'gepol  '//c_null_char
integer, public, save  :: pcmmod_patch_level = 2
real(8), public, save  :: pcmmod_coarsity = 0.5
real(8), public, save  :: pcmmod_cavity_area = 0.3
real(8), public, save  :: pcmmod_min_distance = 0.1
integer, public, save  :: pcmmod_der_order = 4
logical, public, save :: pcmmod_scaling = .true. 
character(len=8), public, save :: pcmmod_radii_set = 'bondi  '//c_null_char
character(len=20), public, save :: pcmmod_restart_name = 'cavity.npz         '//c_null_char
real(8), public, save  :: pcmmod_min_radius = 100.0
! solver specification *PCMSOL section
character(len=7), public, save :: pcmmod_solver_type = 'iefpcm'//c_null_char
character(len=16), public, save :: pcmmod_solvent = '               '//c_null_char
character(len=11), public, save :: pcmmod_equation_type = 'secondkind'//c_null_char
real(8), public, save  :: pcmmod_correction = 0.0
real(8), public, save  :: pcmmod_probe_radius = 1.0
! green specification *PCMSOL section
character(len=7), public, save :: pcmmod_inside_type = 'vacuum'//c_null_char
character(len=22), public, save :: pcmmod_outside_type = 'uniformdielectric    '//c_null_char
real(8), public, save :: pcmmod_outside_epsilon = 1.0 

type, bind(c) :: cavityInput
        character(kind=c_char) :: cavity_type(8)
        integer(c_int) :: patch_level
        real(c_double) :: coarsity
        real(c_double) :: area
        real(c_double) :: min_distance
        integer(c_int) :: der_order
        logical(c_bool) :: scaling
        character(kind=c_char) :: restart_name(20)
        character(kind=c_char) :: radii_set(8)
        real(c_double) :: min_radius
end type

type, bind(c) :: solverInput
        character(kind=c_char) :: solver_type(7)
        character(kind=c_char) :: solvent(16)
        character(kind=c_char) :: equation_type(11) 
        real(c_double)         :: correction
        real(c_double)         :: probe_radius
end type

type, bind(c) :: greenInput
        character(kind=c_char) :: inside_type(7)
        real(c_double)         :: outside_epsilon
        character(kind=c_char) :: outside_type(22)
end type

contains

!> \brief initializes PCM input section to its default values
!> \author R. Di Remigio
!> \date 2014
!> \param pcm_input PCM input section   
subroutine ls_pcm_init(pcm_input)
  
  type(pcmtype)           :: pcm_input
  ! Polarizable continuum model calculation is turned off by default
  pcm_input%do_pcm = .false.
  ! Use of separate charges and potentials is turned off by default 
  pcm_input%separate = .false.
  ! Print level is set to 0 by default
  pcm_input%print_level = 0
  pcm_config = pcm_input

end subroutine ls_pcm_init

!> \brief processes PCM input section to its default values
!> \author R. Di Remigio
!> \date 2014
!> \param pcm_input PCM input section   
!> \param readword
!> \param keyword
!> \param lucmd
!> \param lupri
subroutine ls_pcm_input(pcm_input, readword, keyword, lucmd, lupri)
  
   use lstiming                                                         
   use ls_util, only: lsheader
   type(pcmtype)           :: pcm_input 
   integer                 :: lucmd, lupri
   integer                 :: NError, NWarning
   integer, dimension(1:2) :: KWordInd
   character(len = 1)      :: Prompt
   character(len = 70)     :: Keyword 
   character(len = 7), dimension(1:2) :: kwordtable = &
     (/ '.SEPARA', '.PRINT ' /)
   logical :: file_exist, const_ts
   logical, intent(inout)   :: readword
   integer                  :: filestatus
   
   pcm_input%host_provides_input = .false.
   call lsheader(lupri, &
   'Input processing for the Polarizable Continuum Model module')
   nerror          = 0
   nwarning        = 0
   ReadKeywords: do 
      Read(lucmd,'(A70)', IOStat = FileStatus) Keyword
      !      WRITE(*,*)'READ KEYWORDS',Keyword(1:7)
      If (FileStatus > 0) Call LSQuit('Error reading lucmd', lupri)
      If ((FileStatus < 0)) Exit
      Prompt = Keyword(1:1)
      If ((Prompt == '#') .or. (Prompt == '!')) Cycle
      If (Prompt == '*') then
         ReadWord = .FALSE.
         Exit
      Endif
      If (Prompt == '.') Then
        If (Any(KWordTable(1:size(KWordTable)) .EQ. Keyword(1:7))) Then
          Select Case (Keyword(1:7))
            ! Use separate or total potential and charges
            case('.SEPARA')
              pcm_input%separate = .true.
            ! Print level 
            case('.PRINT ')
              read(lucmd,*) pcm_input%print_level
            ! Case Default
              Write(lupri,'(/,3A,/)') ' Keyword "',Keyword, &
                '" is not yet implemented in PCM_Input_Proc.'
          End Select
        Else
          Write(lupri,'(/,3A,/)') ' Keyword "',Keyword, &
            '" not recognized in PCM_Input_Proc.'
          Call LSQuit('Illegal keyword in PCM_Input_Proc.',lupri)
        End If
      Endif
   enddo ReadKeywords
   pcm_config = pcm_input
   call report_after_pcm_input(lupri, pcm_input)

end subroutine ls_pcm_input
   
!> \brief print relevant setting of both LSDALTON and PCMSolver 
!> \author R. Di Remigio
!> \date 2014
!> \param print_unit the printing unit to be used
!> \param pcm_input PCM input section   
subroutine report_after_pcm_input(print_unit, pcm_cfg)
   
   use ls_pcm_write, only: init_host_writer
   integer, optional, intent(in) :: print_unit
   type(pcmtype), intent(in)     :: pcm_cfg

   if (present(print_unit)) then
      global_print_unit = print_unit                                                                                                 
      ! Initialize host writer                                                                                                       
      call init_host_writer(global_print_unit)                                                                                       
      write(global_print_unit, *) &
      ' ===== Polarizable Continuum Model calculation set-up ====='                                      
      write(global_print_unit, *) &
      '* Polarizable Continuum Model using PCMSolver external module:'                                   
      write(global_print_unit, *) &
      '  1: Converged potentials and charges at tesserae representative points written on file.'         
                                                                                                                                 
      if (pcm_cfg%separate) then                                                                                                     
         write(global_print_unit, *) &
         '  2: Separate potentials and apparent charges in nuclear and electronic.'                      
      else                                                                                                                           
         write(global_print_unit, *) &
         '  2: Use total potentials and apparent charges.'                                               
      end if                                                                                                                         
                                                                                                                                     
      if (pcm_cfg%print_level > 5 .and. pcm_cfg%print_level < 10) then                                                               
         write(global_print_unit, *) &
         '  3: Print potentials at tesserae representative points.'                                      
      else if (pcm_cfg%print_level > 10) then                                                                                        
         write(global_print_unit, *) &
         '  3: Print potentials and charges at tesserae representative points.'                          
      else                                                                                                                           
         write(global_print_unit, *) &
         '  3: Do not print potentials and charges.'                                                     
      end if                                                                                                                         
      call print_citation                                                                                                            
   else
      write(*, *) &
      ' ===== Polarizable Continuum Model calculation set-up ====='                                              
      write(*, *) &
      '* Polarizable Continuum Model using PCMSolver external module:'                                           
      write(*, *) &
      '  1: Converged potentials and charges at tesserae representative points written on file.'                 
                                                                                                                             
      if (pcm_cfg%separate) then                                                                                             
         write(*, *) &
         '  2: Separate potentials and apparent charges in nuclear and electronic.'                              
      else                                                                                                                   
         write(*, *) &
         '  2: Use total potentials and apparent charges.'                                                       
      end if                                                                                                                 
                                                                                                                         
      if (pcm_cfg%print_level > 5 .and. pcm_cfg%print_level < 10) then                                                       
         write(*, *) &
         '  3: Print potentials at tesserae representative points.'                                              
      else if (pcm_cfg%print_level > 10) then                                                                                
         write(*, *) &
         '  3: Print potentials and charges at tesserae representative points.'                                  
      else                                                                                                                   
         write(*, *) &
         '  3: Do not print potentials and charges.'                                                             
      end if                                                                                                                 
      call print_citation                                                                                                    
   end if

end subroutine report_after_pcm_input
     
!> \brief sets PCMSolver input parameters from LSDALTON input 
!> \author R. Di Remigio
!> \date 2014
!> \param cavity struct holding cavity parameters
!> \param solver struct holding solver parameters
!> \param green  struct holding Green's functions parameters
subroutine host_input(cavity, solver, green) bind(c, name='host_input')
! Performs syntactic checks on PCMSolver input and fills the data structures
! holding input data
! WARNING ! Do NOT change the order in which strings are passed since the
! module relies on that
                                                      
   type(cavityInput) :: cavity
   type(solverInput) :: solver
   type(greenInput)  :: green
   
   cavity%patch_level  = pcmmod_patch_level
   cavity%coarsity     = pcmmod_coarsity
   cavity%area         = pcmmod_cavity_area
   cavity%min_distance = pcmmod_min_distance
   cavity%der_order    = pcmmod_der_order
   cavity%scaling      = pcmmod_scaling
   cavity%min_radius   = pcmmod_min_radius
   ! Pass the strings relevant to the cavity input section
   call push_input_string(pcmmod_cavity_type//c_null_char)
   call push_input_string(pcmmod_radii_set//c_null_char)
   call push_input_string(pcmmod_restart_name//c_null_char)
   
   solver%correction        = pcmmod_correction
   solver%probe_radius      = pcmmod_probe_radius
   ! Pass the strings relevant to the solver input section
   call push_input_string(pcmmod_solver_type//c_null_char)
   call push_input_string(pcmmod_solvent//c_null_char)
   call push_input_string(pcmmod_equation_type//c_null_char)
   
   green%outside_epsilon = pcmmod_outside_epsilon
   ! Pass the strings relevant to the green input section
   call push_input_string(pcmmod_inside_type//c_null_char)
   call push_input_string(pcmmod_outside_type//c_null_char)
                                                      
end subroutine host_input

end module ls_pcm_config
