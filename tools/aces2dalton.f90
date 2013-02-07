program aces2dalton
!
!  Reformat a basis set from the PNNL library.
!  The "Dalton" format from the latter is so utterly
!  useless from a Dalton library context that is is
!  preferable/easier to export from the library in
!  AcesII (C4) format and reformat from there.
!
!  Remember to ensure "Optimize general contractions" is checked,
!  save the file as a text file from within the browser window, 
!  and run the program as a filter
!    aces2dalton <inputfile >outputfile
!  A file can contain basis sets for as many elements as you like.
!  Move the outputfile to a directory where Dalton will find it
!  (see the Dalton documentation).
!
! Limits and limitations:
!   See the modules "use"d below for dimensioning of tables
!   of element names and angular momenta.
!   It is assumed that no basis set will contain more than 100 exponents
!   of a given angular type.
!   No more than 7 decimal places are possible for any exponent, and
!   no more than 6 for a contraction coefficient.  For exponents larger
!   than 1,000,000 only 3 decimal places are available.
!

  use io_channels
  use periodic_table
  use angular_momenta

  implicit none

  character(LEN=1000)  :: line, basis

  integer              :: length, line_length, max_length, status

  character(LEN=2)     :: atom

  integer              :: atomic_number, l, l_max, chars_out
  integer              :: i, j, k
  integer              :: ang_mom_values(MAX_ANG_MOM)
  integer              :: num_prim(MAX_ANG_MOM), num_contr(MAX_ANG_MOM)

  real (kind(0.d0))    :: exponents(100), contraction_coeffts(100,100)

  character(LEN=1)     :: upcase, downcase

  logical              :: alphabetic

! Read file, echoing lines that begin with a "!" but replaced with a "$"

  write(lupri,'(a)') "$ Basis converted using PRT's splendid code..."
  write(lupri,'(a)') '$'


!
! Ignore pointless first line...
!
  read(luinp,'(a)',IOSTAT=status) line

!  write(lupri,'(a,a)') 'Echo: ', line

  do

     read(luinp,'(a)',IOSTAT=status) line

!     write(lupri,'(a,a)') 'Echo: ', line

     if (status .ne. 0) exit

     length = line_length(line)
     max_length = max(max_length,line_length(line))

     if (length .eq. 0) cycle

     if (line(1:1) .eq. '!') then

        line(1:1) = '$'
        write(lupri,'(a)') line(1:length)
        cycle

     else
!
! New basis
!
     atom(1:1) = upcase(line(1:1))
     atom(2:2) = ' '
     if (line(2:2) .ne. ':') atom(2:2) = downcase(line(2:2))
     endif

!
! Assign atomic number
!
     atomic_number = 0
     do i = 1,MAX_ELTS
        if (elements(i) .ne. atom) cycle
        atomic_number = i
        exit
     enddo

!     write(lupri,'(a,i2)') 'Atomic number ', atomic_number

     if (atomic_number .eq. 0) then
        write(luerr,'(a,a)') &
             ' Error determining atomic number for ', atom
        write(luerr,'(a)') 'Cannot identify atom in basis set file'
        call abort
     endif

     if (atomic_number .lt. 10) then
        write(lupri,'(a2,i1)') 'a ', atomic_number
     else
        write(lupri,'(a2,i2)') 'a ', atomic_number
     endif

     write(lupri,'(a,a)') '$ ', trim(element_names(atomic_number))

! Should be pointless comment line followed by blank line

     read(luinp,'(a)') line
!     write(lupri,'(a,a)') 'Echo: ', line

     read(luinp,'(a)') line
!     write(lupri,'(a,a)') 'Echo: ', line


     if (line_length(line) .ne. 0) then
        write(luerr,'(a)') 'Expected blank line after atom identifier line!'
        call abort
     endif

     read(luinp,*) l_max

!     write(lupri,'(a,i2)') 'Max (real) l', l_max

! In practice this is l+1...

     if (l_max .gt. MAX_ANG_MOM) then
          write(luerr,'(a)') 'Angular momentum too high for this dimensioning!'
          call abort
     endif

     read(luinp,'(10i5)') ang_mom_values
!     write(lupri,'(a,10i5)') 'Ang mom values', (ang_mom_values(i),i=1,l_max)
     read(luinp,'(10i5)') num_contr
!     write(lupri,'(a,10i5)') 'Contracted functions', (num_contr(i),i=1,l_max)
     read(luinp,'(10i5)') num_prim
!     write(lupri,'(a,10i5)') 'Primitive functions ', (num_prim(i),i=1,l_max)

     basis = '$ Primitive basis   ('
     chars_out = 21

     do i = 1,l_max
        if (num_prim(i) .gt. 9) then
           write(basis(chars_out+1:chars_out+4),'(i2,a,a)') num_prim(i), downcase(l_labels(i)), ','
           chars_out = chars_out+4
        else
           write(basis(chars_out+1:chars_out+3),'(i1,a,a)') num_prim(i), downcase(l_labels(i)), ','
           chars_out = chars_out+3
        endif
     enddo
!  Back up over last comma
     write(basis(chars_out:chars_out),'(a)') ')'

     write(lupri,'(a)') basis(1:chars_out)

     basis = '$ Contracted basis  ['
     chars_out = 21

     do i = 1,l_max
        if (num_contr(i) .gt. 9) then
           write(basis(chars_out+1:chars_out+4),'(i2,a,a)') num_contr(i), downcase(l_labels(i)), ','
           chars_out = chars_out+4
        else
           write(basis(chars_out+1:chars_out+3),'(i1,a,a)') num_contr(i), downcase(l_labels(i)), ','
           chars_out = chars_out+3
        endif
     enddo
!  Back up over last comma
     write(basis(chars_out:chars_out),'(a)') ']'

     write(lupri,'(a)') basis(1:chars_out)

!     write(lupri,'(a,10i3)') '$ Primitives by l value  ', (num_prim(i),i=1,l_max)
!     write(lupri,'(a,10i3)') '$ Contractions by l value', (num_contr(i),i=1,l_max)

     read(luinp,'(a)') line

     if (line_length(line) .ne. 0) then
        write(luerr,'(a)') 'Expected blank line after angmom/prim/contraction counts!'
        call abort
     endif

     do i = 1, l_max

        if (ang_mom_values(i)+1 .ne. i) then
             write(luerr,'(a)') 'Angular momenta must be increasing from S, no gaps!'
             call abort
        endif

        if (num_contr(i) .eq. 0) cycle

        read(luinp,'(5f14.0)') (exponents(j),j = 1,num_prim(i))

!        write(lupri,'(a,/,(5f14.6))') 'Exponents', &
!             (exponents(j),j = 1,num_prim(i))

        read(luinp,'(a)') line

        if (line_length(line) .ne. 0) then
             write(luerr,'(a)') 'Expected blank line after exponents!'
             call abort
        endif

        do j = 1,num_prim(i)
           read(luinp,'(7f11.0)') &
                (contraction_coeffts(k,j),k = 1,num_contr(i))
!           write(lupri,'(a,i2,/,(5f10.4))') ' Contractions for prim', j,&
!                (contraction_coeffts(k,j),k = 1,num_contr(i))


        enddo

        write(lupri,'(a,a,a)') '$ ', l_labels(i), '-TYPE FUNCTIONS'

        write(lupri,'(3i5)') num_prim(i), num_contr(i), 0

        do j = 1,num_prim(i)

           if (exponents(j) .lt. 1.d6) then

              write(lupri,'(f15.7,6f12.8)') exponents(j), &
                        (contraction_coeffts(k,j),k = 1,num_contr(i))

           else

              if (num_contr(i) .lt. 7) then 
                 write(lupri,'(f15.3,6f12.8)') exponents(j), &
                          (contraction_coeffts(k,j),k = 1,num_contr(i))
              else
                 write(lupri,'(f15.3,6f12.8)') exponents(j), &
                          (contraction_coeffts(k,j),k = 1,6)
                 write(lupri,'(f15.7,6f12.8)') &
                          (contraction_coeffts(k,j),k = 7,num_contr(i))
              endif

           endif

        enddo

        read(luinp,'(a)') line

        if (line_length(line) .ne. 0) then
             write(luerr,'(a)') 'Expected blank line after contraction coefficients!'
             call abort
        endif

     enddo

  enddo

end program aces2dalton

function upcase(chr)

  implicit none

  character(LEN=1) :: chr, upcase

  upcase = chr

  if ((ichar(chr) - ichar('a') + 1) .gt. 0 .and. &
      (ichar(chr) - ichar('a') + 1) .le. 26) &
     upcase = char(ichar(chr) - ichar('a') + ichar('A'))

  return

end function upcase

function downcase(chr)

  implicit none

  character(LEN=1) :: chr, downcase

  downcase = chr

  if ((ichar(chr) - ichar('A') + 1) .gt. 0 .and. &
      (ichar(chr) - ichar('A') + 1) .le. 26) &
     downcase = char(ichar(chr) - ichar('A') + ichar('a'))

  return

end function downcase

function line_length(line)

  implicit none

  integer       :: line_length

  character(LEN=*) :: line

  line_length = len_trim(line)

  return

end function line_length
