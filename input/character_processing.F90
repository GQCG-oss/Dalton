
! (c) Radovan Bast and Stefan Knecht
! licensed under the GNU Lesser General Public License

module character_processing

   implicit none

   public lowercase
   public uppercase
   public prefix_zeros
   public word_count
   public word_contains

   private

contains

   function lowercase(s)

!     --------------------------------------------------------------------------
      character(*), intent(in) :: s
      character(len(s))        :: lowercase
!     --------------------------------------------------------------------------
      integer                  :: off, i, ia
!     --------------------------------------------------------------------------

      lowercase = s

      off = iachar('a') - iachar('A')

      do i = 1, len(s)
         ia = iachar(s(i:i))
         if (ia >= iachar('A') .and. ia <= iachar('Z')) then
            lowercase(i:i) = achar(ia + off)
         end if
      enddo

   end function

   function uppercase(s)

!     --------------------------------------------------------------------------
      character(*), intent(in) :: s
      character(len(s))        :: uppercase
!     --------------------------------------------------------------------------
      integer                  :: off, i, ia
!     --------------------------------------------------------------------------

      uppercase = s

      off = iachar('A') - iachar('a')

      do i = 1, len(s)
         ia = iachar(s(i:i))
         if (ia >= iachar('a') .and. ia <= iachar('Z')) then
            uppercase(i:i) = achar(ia + off)
         end if
      enddo

   end function

   function prefix_zeros(i, n)

!     prefix_zeros(137, 6) returns '000137'

!     --------------------------------------------------------------------------
      integer,      intent(in) :: i
      integer,      intent(in) :: n
      character(n)             :: prefix_zeros
!     --------------------------------------------------------------------------
      integer                  :: k
      character(1)             :: c09(0:9) = (/'0','1','2','3','4','5','6','7','8','9'/)
!     --------------------------------------------------------------------------

      do k = 1, n
         prefix_zeros(n-k+1:n-k+1) = c09(mod(i, 10**k)/10**(k-1))
      end do

   end function

   function word_count(s)

!     --------------------------------------------------------------------------
      character(*), intent(in) :: s
      integer                  :: word_count
!     --------------------------------------------------------------------------
      integer                  :: i
      logical                  :: is_blank
!     --------------------------------------------------------------------------

      word_count = 0

      if (len(s) <= 0) return

      is_blank = .true.

      do i = 1, len(s)
         if (s(i:i) == ' ') then
            is_blank = .true.
         else if (is_blank) then
            word_count = word_count + 1
            is_blank = .false.
         end if
      end do

   end function

   function word_contains(word, substring)

!     --------------------------------------------------------------------------
      character(*), intent(in) :: word
      character(*), intent(in) :: substring
!     --------------------------------------------------------------------------
      logical                  :: word_contains
!     --------------------------------------------------------------------------

      word_contains = .false.
      if (index(word, substring) > 0) then
         word_contains = .true.
      end if

   end function

end module
