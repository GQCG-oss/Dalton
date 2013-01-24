program xyz2dalton
!
!.... Convert a typical PDB XYZ file to a dalton-type input format
!
!.... Program runs as a filter (old-school UNIX), reading from stdin
!.... and writing to stdout, and is invoked as
!....     xyz2dalton [-c] < XYZinputfile > DaltonMolfile
!.... where the optional argument -c if present causes the
!.... output geometry to be translated to a centre-of-mass
!.... coordinate system origin.
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
  integer                   :: atoms_of_this_type(MAX_ELTS)

  integer :: num_atoms, num_types

  integer :: arg_count

  character(LEN=2)   :: symbol
  integer            :: atomic_number
  real (kind(0.d0))  :: x, y, z
  real (kind(0.d0))  :: total_mass
  real (kind(0.d0))  :: x_com, y_com, z_com

  logical            :: com_transform



  real (kind(0.d0)), parameter  :: angstrom_to_a0=0.5291772083D0

  character(LEN=2)   :: current_symbol
  integer            :: current_count
  logical            :: one_symbol
  character(LEN=5)   :: charge(MAX_ELTS)
  character(LEN=5)   :: count(MAX_ELTS)

  integer :: i, j

! Check for command-line argument for transformation to centre of mass

  com_transform = command_argument_count() .gt. 0

  x_com = 0.d0
  y_com = 0.d0
  z_com = 0.d0

  atoms_of_this_type = 0

  read(luinp,*) num_atoms
  read(luinp,'(a)') input_line

  write(lupri,'(a)') 'BASIS'
  write(lupri,'(a)') 'Choose basis here!'
  write(lupri,'(a)') input_line
  write(lupri,'(a)') 'PDB structure'
!  write(lupri,'(a)') 'Move Atomtypes line here!'

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

        x_com = x_com + atomic_weights(atomic_number)*x
        y_com = y_com + atomic_weights(atomic_number)*y
        z_com = z_com + atomic_weights(atomic_number)*z
        total_mass = total_mass + atomic_weights(atomic_number)

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

            x_com = x_com + atomic_weights(atomic_number)*x
            y_com = y_com + atomic_weights(atomic_number)*y
            z_com = z_com + atomic_weights(atomic_number)*z
            total_mass = total_mass + atomic_weights(atomic_number)

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

   x_com = x_com/total_mass
   y_com = y_com/total_mass
   z_com = z_com/total_mass

! Now write output file...

   select case(num_types)
   case(1:9)
      write(lupri,'(a,i1)') 'Atomtypes=', num_types
   case(10:99)
      write(lupri,'(a,i2)') 'Atomtypes=', num_types
   case(100:)
      write(lupri,'(a,i3)') 'Atomtypes=', num_types
   end select


   do atomic_number = MAX_ELTS,1,-1

      if (atoms_of_this_type(atomic_number) .gt. 0) then
         write(lupri,'(a,a,1x,a,a)') &
           'Charge=', charge(atomic_number), &
           'Atoms=', count(atomic_number)

         do j = 1,atoms_of_this_type(atomic_number)

            if (com_transform) then

               read(output_line(j,atomic_number),'(4x,3(4x,f14.6))') x, y, z

               x = x - x_com
               y = y - y_com
               z = z - z_com

               write(output_line(j,atomic_number),'(4x,3(4x,f14.6))') x, y, z

            endif

            write(lupri,'(a)') trim(output_line(j,atomic_number))

         enddo
      endif

   end do

end program xyz2dalton



