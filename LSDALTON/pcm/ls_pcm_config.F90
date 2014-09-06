!> @file
!> configure PCM input LSDALTON-side 
module ls_pcm_config

implicit none

private

public pcmtype
public ls_pcm_init
public ls_pcm_input

type pcmtype
     logical :: do_pcm
     logical :: separate
     integer :: print_level
end type pcmtype

type(pcmtype), public    :: pcm_config
integer                  :: global_print_unit = -1

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
      write(global_print_unit, *) ' ===== Polarizable Continuum Model calculation set-up ====='                                      
      write(global_print_unit, *) '* Polarizable Continuum Model using PCMSolver external module:'                                   
      write(global_print_unit, *) '  1: Converged potentials and charges at tesserae representative points written on file.'         
                                                                                                                                 
      if (pcm_cfg%separate) then                                                                                                     
         write(global_print_unit, *) '  2: Separate potentials and apparent charges in nuclear and electronic.'                      
      else                                                                                                                           
         write(global_print_unit, *) '  2: Use total potentials and apparent charges.'                                               
      end if                                                                                                                         
                                                                                                                                     
      if (pcm_cfg%print_level > 5 .and. pcm_cfg%print_level < 10) then                                                               
         write(global_print_unit, *) '  3: Print potentials at tesserae representative points.'                                      
      else if (pcm_cfg%print_level > 10) then                                                                                        
         write(global_print_unit, *) '  3: Print potentials and charges at tesserae representative points.'                          
      else                                                                                                                           
         write(global_print_unit, *) '  3: Do not print potentials and charges.'                                                     
      end if                                                                                                                         
      call print_citation                                                                                                            
   else
      write(*, *) ' ===== Polarizable Continuum Model calculation set-up ====='                                              
      write(*, *) '* Polarizable Continuum Model using PCMSolver external module:'                                           
      write(*, *) '  1: Converged potentials and charges at tesserae representative points written on file.'                 
                                                                                                                             
      if (pcm_cfg%separate) then                                                                                             
         write(*, *) '  2: Separate potentials and apparent charges in nuclear and electronic.'                              
      else                                                                                                                   
         write(*, *) '  2: Use total potentials and apparent charges.'                                                       
      end if                                                                                                                 
                                                                                                                         
      if (pcm_cfg%print_level > 5 .and. pcm_cfg%print_level < 10) then                                                       
         write(*, *) '  3: Print potentials at tesserae representative points.'                                              
      else if (pcm_cfg%print_level > 10) then                                                                                
         write(*, *) '  3: Print potentials and charges at tesserae representative points.'                                  
      else                                                                                                                   
         write(*, *) '  3: Do not print potentials and charges.'                                                             
      end if                                                                                                                 
      call print_citation                                                                                                    
   end if

end subroutine report_after_pcm_input

end module ls_pcm_config
