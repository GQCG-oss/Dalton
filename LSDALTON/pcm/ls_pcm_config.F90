module ls_pcm_config

  implicit none
  
  private

  public pcmtype
  public LS_pcm_init
  public LS_pcm_input
  public ls_pcm_write_file
  public ls_pcm_write_file_separate

  type pcmtype
       logical :: do_pcm
       logical :: separate
       integer :: print_level
  end type pcmtype

  type(pcmtype), public    :: pcm_config
  character(11), parameter :: pcm_file_name = 'PCM_mep_asc'
  integer,       parameter :: pcm_file_unit = 800 
  logical                  :: pcm_file_exists = .false.
  logical                  :: pcm_file_open = .false.
  integer                  :: global_print_unit = -1

contains

subroutine LS_pcm_init(pcm_input)
  type(pcmtype)           :: pcm_input
!
! Initialization
!
  ! Polarizable continuum model calculation is turned off by default
  pcm_input%do_pcm = .false.
  ! Use of separate charges and potentials is turned off by default 
  pcm_input%separate = .false.
  ! Print level is set to 0 by default
  pcm_input%print_level = 0
end subroutine LS_pcm_init

Subroutine LS_pcm_input(pcm_input, ReadWord, keyword, lucmd, lupri)
!
! Input processing routine for PCM 
!
  Use LSTiming
  use ls_util , only: lsheader
  Type(pcmtype)           :: pcm_input 
  Integer                 :: lucmd, lupri
  Integer                 :: NError, NWarning
  Integer, Dimension(1:2) :: KWordInd
  Character(Len = 1)      :: Prompt
  Character(Len = 70)     :: Keyword 
  Character(Len = 7), Dimension(1:2) :: KWordTable = &
    (/ '.SEPARA', '.PRINT ' /)
  Logical :: file_exist,Const_TS
  Logical, intent(inout)   :: ReadWord
  Integer :: FileStatus
!
Call LSHeader(lupri, &
'Input processing for the direct dynamics module')
NError          = 0
NWarning        = 0
!
! Read keywords
!
Write(*,*)'LS PCM INPUT'
Do
Write(*,*)'READWORD',ReadWord
      Read(lucmd,'(A70)', IOStat = FileStatus) Keyword
            WRITE(*,*)'READ KEYWORDS',Keyword(1:7)
      If (FileStatus > 0) Call LSQuit('Error reading lucmd',lupri)
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
!            Case Default
              Write(lupri,'(/,3A,/)') ' Keyword "',Keyword, &
                '" is not yet implemented in PCM_Input_Proc.'
          End Select
        Else
          Write(lupri,'(/,3A,/)') ' Keyword "',Keyword, &
            '" not recognized in PCM_Input_Proc.'
          Call LSQuit('Illegal keyword in PCM_Input_Proc.',lupri)
        End If
      Endif
   Enddo
   pcm_config = pcm_input
!
! The chosen parameters are printed
!
  call report_after_pcm_input(lupri, pcm_input)
end subroutine LS_pcm_input
   
subroutine report_after_pcm_input(print_unit, pcm_cfg)
!
! A subroutine to print the relevant settings of both LSDALTON and PCMSolver.
!
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
   ! Should print info for the cavity and the solver here.
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
           
   ! Should print info for the cavity and the solver here.
           call print_citation
   end if
end subroutine report_after_pcm_input

subroutine ls_pcm_write_file(nr_points, potentials, charges)

! Passed variables
   integer, intent(in) :: nr_points
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
 
end subroutine ls_pcm_write_file
   
subroutine ls_pcm_write_file_separate(nr_points, nuc_pot, nuc_chg, ele_pot, ele_chg)

! Passed variables
   integer, intent(in) :: nr_points
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
 
end subroutine ls_pcm_write_file_separate 

end module ls_pcm_config
