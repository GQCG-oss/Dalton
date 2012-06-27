
! (c) Radovan Bast and Stefan Knecht
! licensed under the GNU Lesser General Public License

module input_reader

   use keyword
   use character_processing
   use input_reader_sections

   implicit none

   public read_menu_input
   public move_to_next_star

   private

contains

   subroutine read_menu_input(input_section_in, input_found)

!    ---------------------------------------------------------------------------
     character(*), intent(in)  :: input_section_in
     logical,      intent(out) :: input_found
!    ---------------------------------------------------------------------------
     integer                   :: menu_fh
     character(kw_length)      :: input_section
     character(kw_length)      :: word
     character(kw_length)      :: kw_section
     character(len=10)         :: menu_file
     logical                   :: read_all_input
     logical                   :: file_is_open
!    ---------------------------------------------------------------------------

!    this is to catch keywords that appear before some section starts
     kw_section = '       '

!    check for the input section we would like to read
     input_section = uppercase(input_section_in) ! truncate or extend to kw_length

!    if applicable process the complete input (all modules)
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

     inquire(file=menu_file, opened=file_is_open, number=menu_fh)

     if(file_is_open)then
!       use existing file handle
        unit_in = menu_fh
     else
!       open menu file with unique file unit (defined in module keyword)
        unit_in = unit_in_default
        open(unit_in, file=menu_file)
     end if

     rewind(unit_in)

!    initialize
     input_found = .false.

     do while (.true.)

       read(unit_in, '(a7)', end=1) word

       if (word == '       ') then
!          blank line
       else
          select case (word(1:1))

            case ('!', '#')
!              comment

            case ('*')
!              section
               kw_section = uppercase(word)
               if (word == '*END OF' .or. word == '**END O') go to 1

            case default
!              keyword
               if (read_all_input .or. kw_section == input_section) then
                  input_found = .true.
                  call read_input_sections(word, kw_section)
               end if

          end select
       end if

     end do

1    if(.not. file_is_open)then
        close(unit_in, status='keep')
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
