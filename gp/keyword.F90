!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!

module keyword

!radovan: please feel free to copy/use/improve this inside DIRAC/Dalton
!         without asking

  use character_processing

  implicit none

  public reset_available_kw_list
  public check_whether_kw_found
  public kw_matches
  public kw_read

  integer, parameter, public :: kw_length       = 7
  integer, parameter, public :: unit_in_default = 5
  integer,            public :: unit_in         = unit_in_default

  private

  interface kw_read
    module procedure kw_read_c
    module procedure kw_read_i1
    module procedure kw_read_i3
    module procedure kw_read_r1
    module procedure kw_read_r2
    module procedure kw_read_r3
    module procedure kw_read_r4
  end interface

  integer, parameter         :: max_nr_kw   = 200
  integer, parameter         :: line_length = 80

  character(kw_length)       :: available_kw_list(max_nr_kw)
  integer                    :: nr_available_kw
  logical                    :: kw_found

contains

  subroutine kw_read_c(kw_input, c)

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in)  :: kw_input
    character(*),         intent(out) :: c
!   ----------------------------------------------------------------------------

    read(unit_in, *, err=1) c
    return

 1  call kw_read_error(kw_input)

  end subroutine

  subroutine kw_read_i1(kw_input, i)

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in)  :: kw_input
    integer,              intent(out) :: i
!   ----------------------------------------------------------------------------

    read(unit_in, *, err=1) i
    return

 1  call kw_read_error(kw_input)

  end subroutine

  subroutine kw_read_i3(kw_input, i1, i2, i3)

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in)  :: kw_input
    integer,              intent(out) :: i1, i2, i3
!   ----------------------------------------------------------------------------

    read(unit_in, *, err=1) i1, i2, i3
    return

 1  call kw_read_error(kw_input)

  end subroutine

  subroutine kw_read_r1(kw_input, r)

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in)  :: kw_input
    real(8),              intent(out) :: r
!   ----------------------------------------------------------------------------

    read(unit_in, *, err=1) r
    return

 1  call kw_read_error(kw_input)

  end subroutine

  subroutine kw_read_r2(kw_input, r1, r2)

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in)  :: kw_input
    real(8),              intent(out) :: r1, r2
!   ----------------------------------------------------------------------------

    read(unit_in, *, err=1) r1, r2
    return

 1  call kw_read_error(kw_input)

  end subroutine

  subroutine kw_read_r3(kw_input, r1, r2, r3)

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in)  :: kw_input
    real(8),              intent(out) :: r1, r2, r3
!   ----------------------------------------------------------------------------

    read(unit_in, *, err=1) r1, r2, r3
    return

 1  call kw_read_error(kw_input)

  end subroutine

  subroutine kw_read_r4(kw_input, r1, r2, r3, r4)

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in)  :: kw_input
    real(8),              intent(out) :: r1, r2, r3, r4
!   ----------------------------------------------------------------------------

    read(unit_in, *, err=1) r1, r2, r3, r4
    return

 1  call kw_read_error(kw_input)

  end subroutine

  subroutine kw_read_error(kw_input)

!   ----------------------------------------------------------------------------
    character(kw_length),   intent(in)  :: kw_input
!   ----------------------------------------------------------------------------
    character(line_length)              :: line
!   ----------------------------------------------------------------------------

    backspace unit_in
    read(unit_in, *) line

    write(*, *) 'error in input line:'
    write(*, *) line
    write(*, *) 'following keyword '//kw_input

    call quit('error in line following keyword '//kw_input)

  end subroutine

  subroutine reset_available_kw_list()

    nr_available_kw = 0
    kw_found        = .false.

  end subroutine

  subroutine check_whether_kw_found(kw_input, kw_section)

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in) :: kw_input
    character(kw_length), intent(in) :: kw_section
!   ----------------------------------------------------------------------------
    integer                          :: i
!   ----------------------------------------------------------------------------

    if (.not. kw_found) then

      write(*, *) 'illegal keyword '//kw_input//' in section '//kw_section

      write(*, *) 'list of available keywords in section '//kw_section//':'
      do i = 1, nr_available_kw
        write(*, *) available_kw_list(i)
      end do

      call quit('illegal keyword '//kw_input//' in section '//kw_section)

    end if

  end subroutine

  function kw_matches(kw_input, kw_option)

!   ----------------------------------------------------------------------------
    character(kw_length), intent(in) :: kw_input
    character(kw_length), intent(in) :: kw_option
    logical                          :: kw_matches
!   ----------------------------------------------------------------------------

    if (lowercase(kw_input) == lowercase(kw_option)) then
      kw_matches = .true.
      kw_found   = .true.
    else
      kw_matches = .false.
      nr_available_kw = nr_available_kw + 1
      available_kw_list(nr_available_kw) = kw_option
    end if

  end function

end module
