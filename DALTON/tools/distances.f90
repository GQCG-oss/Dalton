program distances
!
!.... Convert a typical PDB XYZ file to a list of internuclear distances
!
!.... Program runs as a filter (old-school UNIX), reading from stdin
!.... and writing to stdout, and is invoked as
!....     distances [--min lower_bound] [--max upper_bound] < XYZinputfile > Distancesfile
!.... where the optional arguments screen out internuclear
!.... distances greater than upper_bound and less than lower_bound.
!
!.... If --min is omitted, lower_bound is taken as zero.
!.... If --max is omitted, upper_bound is taken as (effectively) infinity.
!
!.... The less explicit, old-school form
!....     distances [-b lower_bound] [-t upper_bound] < XYZinputfile > Distancesfile
!.... will also work.
!
!.... Assumptions: 
!....    All atoms of a given atomic number may not be grouped together.
!....    There are at most 99 atoms of a given type, except for
!....    those with a one-letter atomic symbol.
!....    No more than 999 atoms total.
!....    Highest atomic number handled is 112 (Cn)
!

  use io_channels
  use periodic_table

  implicit none

  integer, parameter        :: MAX_LENGTH=80
  character(LEN=MAX_LENGTH) :: input_line
  character(LEN=MAX_LENGTH) :: output_line(999,MAX_ELTS)
  character(LEN=MAX_LENGTH) :: temp_line(999,MAX_ELTS)
  integer                   :: atoms_of_this_type(MAX_ELTS)

  integer :: num_atoms, num_types

  integer :: arg_count

  character(LEN=2)   :: symbol
  integer            :: atomic_number
  integer            :: total_atoms
  real (kind(0.d0))  :: x, y, z
  real (kind(0.d0))  :: screening_threshold, baseline_threshold
  real (kind(0.d0))  :: x_vec(999), y_vec(999), z_vec(999)
  real (kind(0.d0))  :: dist

  character(LEN=4)   :: name_vec(999)

  logical            :: com_transform



  real (kind(0.d0)), parameter  :: angstrom_to_a0=0.5291772083D0

  character(LEN=2)   :: current_symbol
  integer            :: current_count
  logical            :: one_symbol
  character(LEN=5)   :: charge(MAX_ELTS)
  character(LEN=5)   :: count(MAX_ELTS)

  integer :: i, j, k

! Check for command-line argument for screening threshold

  screening_threshold = 1.d8
  baseline_threshold  = -0.1d0
  if (command_argument_count() .eq. 2) then
     call get_command_argument(1,input_line)
     if (input_line(1:2) .eq. '-t' &
         .or. input_line(1:5) .eq. '--max') then
        call get_command_argument(2,input_line)
        read(input_line,*) screening_threshold
     elseif (input_line(1:2) .eq. '-b' &
         .or. input_line(1:5) .eq. '--min') then
        call get_command_argument(2,input_line)
        read(input_line,*) baseline_threshold
     endif
  endif

  if (command_argument_count() .eq. 4) then
    call get_command_argument(1,input_line)
     if (input_line(1:2) .eq. '-t' &
         .or. input_line(1:5) .eq. '--max') then
        call get_command_argument(2,input_line)
        read(input_line,*) screening_threshold
     elseif (input_line(1:2) .eq. '-b' &
         .or. input_line(1:5) .eq. '--min') then
        call get_command_argument(2,input_line)
        read(input_line,*) baseline_threshold
     endif

    call get_command_argument(3,input_line)
     if (input_line(1:2) .eq. '-t' &
         .or. input_line(1:5) .eq. '--max') then
        call get_command_argument(4,input_line)
        read(input_line,*) screening_threshold
     elseif (input_line(1:2) .eq. '-b' &
         .or. input_line(1:5) .eq. '--min') then
        call get_command_argument(4,input_line)
        read(input_line,*) baseline_threshold
     endif

  endif

  atoms_of_this_type = 0

  read(luinp,*) num_atoms
  read(luinp,'(a)') input_line

  num_types = 0

  current_symbol = '  '
  current_count = 0

  do i = 1,num_atoms
     read(luinp,'(a)') input_line
!     write(luerr,'(a)') input_line

     symbol = input_line(1:2)

     if (symbol .eq. current_symbol) then
        current_count = current_count + 1
        atoms_of_this_type(atomic_number) = current_count

        read(input_line(3:MAX_LENGTH),*) x, y, z

        x = x/angstrom_to_a0
        y = y/angstrom_to_a0
        z = z/angstrom_to_a0

        if (one_symbol) then
           write(output_line(current_count,atomic_number),'(a1,i3.3,3(4x,f14.6))') &
                trim(current_symbol), current_count, x, y, z
        else
           write(output_line(current_count,atomic_number),'(a2,i2.2,3(4x,f14.6))') &
                current_symbol, current_count, x, y, z
        endif

     else

!....  Different element: write out any preceding block before locating it.

        if (current_count .ne. 0) then
           select case(atomic_number)
           case(1:9)
                write(charge(atomic_number),'(f3.1,a)') real(atomic_number)+0.001, '  '
           case(10:99)
                write(charge(atomic_number),'(f4.1,a)') real(atomic_number)+0.001, ' '
           case(100:)
                write(charge(atomic_number),'(f5.1,a)') real(atomic_number)+0.001
           end select

           select case(current_count)
           case(1:9)
              write(count(atomic_number),'(i1)') current_count
           case(10:99)
              write(count(atomic_number),'(i2)') current_count 
           case(100:)
              write(count(atomic_number),'(i3)') current_count
           end select

         endif
           
!....  Locate new element
         do j = 1,MAX_ELTS
            if (symbol .ne. elements(j)) cycle

            current_symbol = symbol
            one_symbol = len(trim(symbol)) .eq. 1
            atomic_number = j
            if (atoms_of_this_type(atomic_number) .eq. 0) num_types = num_types + 1
            current_count = atoms_of_this_type(atomic_number) + 1
            atoms_of_this_type(atomic_number) = current_count

            read(input_line(3:MAX_LENGTH),*) x, y, z

            x = x/angstrom_to_a0
            y = y/angstrom_to_a0
            z = z/angstrom_to_a0

            if (one_symbol) then
               write(output_line(current_count,atomic_number),'(a1,i3.3,3(4x,f14.6))') &
                    trim(current_symbol), current_count, x, y, z
            else
               write(output_line(current_count,atomic_number),'(a2,i2.2,3(4x,f14.6))') &
                    current_symbol, current_count, x, y, z
            endif

            exit
         enddo

      endif

   enddo

!....  Write last block

   if (current_count .ne. 0) then
      select case(atomic_number)
      case(1:9)
         write(charge(atomic_number),'(f3.1,a)') real(atomic_number)+0.001, '  '
      case(10:99)
         write(charge(atomic_number),'(f4.1,a)') real(atomic_number)+0.001, ' '
      case(100:)
         write(charge(atomic_number),'(f5.1,a)') real(atomic_number)+0.001
      end select

      select case(current_count)
      case(1:9)
         write(count(atomic_number),'(i1)') current_count
      case(10:99)
         write(count(atomic_number),'(i2)') current_count 
      case(100:)
         write(count(atomic_number),'(i3)') current_count
      end select

   endif

! Now generate internuclear distances and write out

   total_atoms = 0

   do atomic_number = MAX_ELTS,1,-1

      if (atoms_of_this_type(atomic_number) .gt. 0) then


         do j = 1,atoms_of_this_type(atomic_number)

            total_atoms = total_atoms + 1

            read(output_line(j,atomic_number),'(a4,3(4x,f14.6))') &
                 name_vec(total_atoms), x_vec(total_atoms), &
                 y_vec(total_atoms), z_vec(total_atoms)

         enddo

      endif

   enddo

   write(lupri,'(/,a,/)') 'Table of internuclear distances in Angstrom'
   write(lupri,'(4x,a4,4x,a4,14x,a8,/)') 'Atom', 'Atom', 'Distance'

   do j = 2,total_atoms

      do k = 1,j-1

         dist = angstrom_to_a0 &
                  *sqrt((x_vec(k) - x_vec(j))**2 & 
                      + (y_vec(k) - y_vec(j))**2 & 
                      + (z_vec(k) - z_vec(j))**2)

         if (dist .lt. screening_threshold .and. dist .gt. baseline_threshold) &
              write(lupri,'(4x,a4,4x,a4,8x,f14.6)') &
                name_vec(j), name_vec(k), dist

      enddo

   enddo

end program distances



