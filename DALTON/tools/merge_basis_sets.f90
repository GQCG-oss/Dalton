program merge_basis_sets
!
!  This program is a (quick) modification of Peter Taylor's
!  aces2dalton script. This program merges 2 basis set files
!  in AcesII format from the EMSL library to one file in the 
!  (LS)Dalton fornat.
!
!  Author: Peter Taylor, modify by Pablo Baudin
!  Date of modification: July 2015
!
!  Remember to ensure "Optimize general contractions" is checked,
!  save the files as a text file from within the browser window, 
!  and run the program as,
!    ./merge_basis_sets input1 input2 > outputfile
!  The files can contain basis sets for as many elements as you like.
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

  character(LEN=1000)  :: line1, line2, basis
  character(LEN=50)    :: filename1, filename2
  character(LEN=2)     :: atom, atom2
  character(LEN=1)     :: upcase, downcase

  integer              :: length1, length2, line_length, stat
  integer              :: atomic_number, chars_out
  integer              :: i, j, k, l, l_max1, l_max2, l_max, iter, inp1, inp2

  integer              :: ang_mom_values1(MAX_ANG_MOM)
  integer              :: ang_mom_values2(MAX_ANG_MOM)
  integer              :: num_prim(MAX_ANG_MOM), num_contr(MAX_ANG_MOM)
  integer              :: num_prim1(MAX_ANG_MOM), num_contr1(MAX_ANG_MOM)
  integer              :: num_prim2(MAX_ANG_MOM), num_contr2(MAX_ANG_MOM)

  real (kind(0.d0))    :: exponents(100), contraction_coeffts(100,100)


  ! get input filenames:
  call getarg(1,filename1)
  call getarg(2,filename2)

  ! Open files to be merged
  inp1=21
  inp2=22
  open(unit=inp1, file=filename1, status="old", FORM="formatted")
  open(unit=inp2, file=filename2, status="old", FORM="formatted")


  ! Loop over all lines in both input files
  ! ***************************************
  iter = 0
  loop_filelines: do

     iter = iter + 1

     ! Read and write initial comments and blancks
     ! *******************************************
     read(inp1,'(a)',IOSTAT=stat) line1
     if (stat .ne. 0) exit loop_filelines
     read(inp2,'(a)',IOSTAT=stat) line2
     if (stat .ne. 0) exit loop_filelines

     length1 = line_length(line1)
     do while ((length1==0).or.(line1(1:1) .eq. '!'))

       if (line1(1:1) .eq. '!') then
          line1(1:1) = '$'
          write(lupri,'(a)') line1(1:length1)
       endif

       read(inp1,'(a)',IOSTAT=stat) line1
       if (stat .ne. 0) exit loop_filelines
       length1 = line_length(line1)

     end do

     if (iter==1) then
        write(lupri,'(a)') '$'
        write(lupri,'(a)') '$     Augmented with'
        write(lupri,'(a)') '$'
     end if

      
     length2 = line_length(line2)
     do while ((length2==0).or.(line2(1:1) .eq. '!'))

       if (line2(1:1) .eq. '!') then
          line2(1:1) = '$'
          write(lupri,'(a)') line2(1:length2)
       endif

       read(inp2,'(a)',IOSTAT=stat) line2
       if (stat .ne. 0) exit loop_filelines
       length2 = line_length(line2)

     end do


     ! A new basis is found and will be treated
     ! ****************************************

     ! Get Atomic symbol
     atom(1:1) = upcase(line1(1:1))
     atom(2:2) = ' '
     if (line1(2:2) .ne. ':') atom(2:2) = downcase(line1(2:2))

     atom2(1:1) = upcase(line2(1:1))
     atom2(2:2) = ' '
     if (line2(2:2) .ne. ':') atom2(2:2) = downcase(line2(2:2))

     if (atom .ne. atom2) then
        write(luerr,*) atom, atom2
        write(luerr,'(a)') &
           ' files input1 and input2 must contain the same atoms in the same order'
        call abort
     end if


     ! Assign atomic number
     atomic_number = 0
     do i = 1,MAX_ELTS
        if (elements(i) .ne. atom) cycle
        atomic_number = i
        exit
     enddo

     !write(lupri,'(a,i2)') 'Atomic number ', atomic_number

     if (atomic_number .eq. 0) then
        write(luerr,'(2a)') ' Error determining atomic number for ', atom
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
     read(inp1,'(a)') line1
     read(inp1,'(a)') line1
     read(inp2,'(a)') line2
     read(inp2,'(a)') line2

     if (line_length(line1) .ne. 0) then
        write(luerr,*) line1, atom
        write(luerr,'(a)') 'Expected blank line after atom identifier line!'
        call abort
     endif
     if (line_length(line2) .ne. 0) then
        write(luerr,*) line2, atom
        write(luerr,'(a)') 'Expected blank line after atom identifier line!'
        call abort
     endif

     ! Read max angular momentum, in practice this is l+1...
     read(inp1,*) l_max1
     read(inp2,*) l_max2

     !write(lupri,'(a,i2)') 'Max (real) l', l_max

     if (l_max1 .gt. MAX_ANG_MOM) then
          write(luerr,'(a)') 'Angular momentum too high for this dimensioning!'
          call abort
     endif
     if (l_max2 .gt. MAX_ANG_MOM) then
          write(luerr,'(a)') 'Angular momentum too high for this dimensioning!'
          call abort
     endif

     ang_mom_values1 = 0
     ang_mom_values2 = 0
     read(inp1,'(10i5)') ang_mom_values1
     read(inp2,'(10i5)') ang_mom_values2
!     write(lupri,'(a,10i5)') 'Ang mom values', (ang_mom_values(i),i=1,l_max)
     num_contr = 0
     num_contr1 = 0
     num_contr1 = 0
     read(inp1,'(10i5)') num_contr1
     read(inp2,'(10i5)') num_contr2
!     write(lupri,'(a,10i5)') 'Contracted functions', (num_contr(i),i=1,l_max)
     num_prim = 0
     num_prim1 = 0
     num_prim1 = 0
     read(inp1,'(10i5)') num_prim1
     read(inp2,'(10i5)') num_prim2
!     write(lupri,'(a,10i5)') 'Primitive functions ', (num_prim(i),i=1,l_max)

     basis = '$ Primitive basis   ('
     chars_out = 21

     ! merge:
     l_max=max(l_max1,l_max2)
     do i=1,l_max
        num_contr(i) = num_contr1(i) + num_contr2(i)
        num_prim(i) = num_prim1(i) + num_prim2(i)
     end do

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

     read(inp1,'(a)') line1
     read(inp2,'(a)') line2

     if (line_length(line1) .ne. 0) then
        write(luerr,*) line1, atom
        write(luerr,'(a)') 'Expected blank line after angmom/prim/contraction counts!'
        call abort
     endif
     if (line_length(line2) .ne. 0) then
        write(luerr,*) line2, atom
        write(luerr,'(a)') 'Expected blank line after angmom/prim/contraction counts!'
        call abort
     endif

     do i = 1, l_max

        !if (ang_mom_values(i)+1 .ne. i) then
        !     write(luerr,'(a)') 'Angular momenta must be increasing from S, no gaps!'
        !     call abort
        !endif

        if (num_contr(i) .eq. 0) cycle

        read(inp1,'(5f14.0)') (exponents(j),j = 1,num_prim1(i))
        if (i<=l_max2) read(inp2,'(5f14.0)') (exponents(j+num_prim1(i)),j = 1,num_prim2(i))

!        write(lupri,'(a,/,(5f14.6))') 'Exponents', &
!             (exponents(j),j = 1,num_prim(i))

        read(inp1,'(a)') line1
        if (i<=l_max2) read(inp2,'(a)') line2

        if (line_length(line1) .ne. 0) then
             write(luerr,*) line1, atom
             write(luerr,'(a)') 'Expected blank line after exponents! 0'
             call abort
        endif
        if (line_length(line2) .ne. 0 .and. (i<=l_max2)) then
             write(luerr,*) line2, atom
             write(luerr,'(a)') 'Expected blank line after exponents! 1'
             call abort
        endif

        do j = 1,num_prim1(i)
           read(inp1,'(7f11.0)') (contraction_coeffts(k,j),k = 1,num_contr1(i))
!           write(lupri,'(a,i2,/,(5f10.4))') ' Contractions for prim', j,&
!                (contraction_coeffts(k,j),k = 1,num_contr(i))

        enddo

        do j = (num_prim1(i)+1),num_prim(i)
           if (i<=l_max2) read(inp2,'(7f11.0)') (contraction_coeffts(k,j),k = (num_contr1(i)+1),num_contr(i))
!           write(lupri,'(a,i2,/,(5f10.4))') ' Contractions for prim', j,&
!                (contraction_coeffts(k,j),k = 1,num_contr(i))

        enddo

        ! put zero in the rest
        do j = 1,num_prim1(i)
           do k=(num_contr1(i)+1),num_contr(i)
              contraction_coeffts(k,j) = 0.0
           end do
        enddo
        do j = (num_prim1(i)+1),num_prim(i)
           do k=1,num_contr1(i)
              contraction_coeffts(k,j) = 0.0
           end do
        enddo

        write(lupri,'(a,a,a)') '$ ', l_labels(i), '-TYPE FUNCTIONS'

        write(lupri,'(3i5)') num_prim(i), num_contr(i), 0

        do j = 1,num_prim(i)

           if (exponents(j) .lt. 1.d6) then

              write(lupri,'(f15.7,60f12.8)') exponents(j), &
                        (contraction_coeffts(k,j),k = 1,num_contr(i))

           else

              if (num_contr(i) .lt. 7) then 
                 write(lupri,'(f15.3,60f12.8)') exponents(j), &
                          (contraction_coeffts(k,j),k = 1,num_contr(i))
              else
                 write(lupri,'(f15.3,60f12.8)') exponents(j), &
                          (contraction_coeffts(k,j),k = 1,6)
                 write(lupri,'(f15.7,60f12.8)') &
                          (contraction_coeffts(k,j),k = 7,num_contr(i))
              endif

           endif

        enddo

        read(inp1,'(a)') line1
        if (i<=l_max2) read(inp2,'(a)') line2

        if (line_length(line1) .ne. 0) then
             write(luerr,*) line1, atom
             write(luerr,'(a)') 'Expected blank line after contraction coefficients! 2'
             call abort
        endif

        if (line_length(line2) .ne. 0 .and. (i<=l_max2) ) then
             write(luerr,*) line2, atom
             write(luerr,'(a)') 'Expected blank line after contraction coefficients! 3'
             call abort
        endif

     enddo

  enddo loop_filelines

  close(unit=inp1)
  close(unit=inp2)

end program merge_basis_sets

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
