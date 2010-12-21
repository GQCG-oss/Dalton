!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

module input_reader

! radovan: - this module should not initialize anything
!            it is a reader
!          - when we have only modules and no common blocks the initialization
!            problem goes away since variables can be initialized at definition

!radovan: please feel free to copy/use/improve this inside DIRAC/Dalton
!         without asking

  use keyword
  use character_processing
#ifdef PRG_DIRAC
  use dirac_input_processing
#else
  use dalton_input_processing
#endif

  implicit none

  public read_menu_input
  public move_to_next_star

  private

contains

  subroutine read_menu_input(input_section_in,input_found)

!   radovan: - this routine should not include any common blocks
!              this should be done in called subroutines
!            - keep ifndef traps inside subroutines for better readability
!              of this subroutine

!   ----------------------------------------------------------------------------
    character(*), intent(in)  :: input_section_in
    logical,      intent(out) :: input_found
!   ----------------------------------------------------------------------------
    integer                   :: menu_fh
    character(kw_length)      :: input_section
    character(kw_length)      :: word
    character(kw_length)      :: kw_section
    character(len=10)         :: menu_file
    logical                   :: read_all_input
    logical                   :: file_is_open
!   ----------------------------------------------------------------------------

!   this is to catch keywords that appear before some section starts
    kw_section = '       '

!   check for the input section we would like to read
    input_section = uppercase(input_section_in) ! truncate or extend to kw_length

!   if applicable process the complete input (all modules)  
    if (input_section(1:3) == 'ALL') then
       read_all_input = .true.
    else
       read_all_input = .false.
    end if

#ifdef PRG_DIRAC
    menu_file = 'DIRAC.INP'
#else
    menu_file = 'DALTON.INP'
#endif

    inquire(file=menu_file,opened=file_is_open,number=menu_fh)

    if(file_is_open)then
!     use existing file handle
      unit_in = menu_fh
    else
!     open menu file with unique file unit (defined in module keyword)
      unit_in = unit_in_default
      open(unit_in, file = menu_file)
    end if

    rewind(unit_in)

!   initialize 
    input_found = .false.
    do while (.true.)

      read(unit_in, '(a7)', end=1) word

      if (word == '       ') then
!         blank line
      else
         select case (word(1:1))
         
           case ('!', '#')
!            comment
         
           case ('*')
!            section
             kw_section = uppercase(word)
             if (word == '*END OF' .or. word == '**END O') go to 1
         
           case default
!            keyword
             if (read_all_input .or. kw_section == input_section) then
               input_found = .true.
#ifdef PRG_DIRAC
               call read_dirac_input(word, kw_section)
#else
               call read_dalton_input(word, kw_section)
#endif
             end if
         
         end select
      end if

    end do

 1  if(.not.file_is_open)then 
      close(unit_in,status='keep')
    end if

    unit_in = unit_in_default

  end subroutine

  subroutine move_to_next_star(word_io, file_unit)

!   ----------------------------------------------------------------------------
    character(kw_length), intent(inout) :: word_io
    integer,              intent(in)    :: file_unit
!   ----------------------------------------------------------------------------
    character(kw_length)                :: word
!   ----------------------------------------------------------------------------

    do while (.true.)
      read(file_unit, '(a7)') word
      if (word(1:1) == '*') then
        word_io = word
        return
      end if
    end do

  end subroutine

 end module
