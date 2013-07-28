! Copyright 2012 Bin Gao
!           2012 Andreas J. Thorvaldsen
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.

! module matrix_lowlevel 
!
! intended to provide low-level functionaly for type(matrix)
module matrix_lowlevel

  use matrix_backend

  implicit none
  
  ! passed on from matrix_backend
  public matrix
  public mat_init          !initialize and allocate
  public mat_remove        !free and nullify
  public mat_acquire_block !asssociate pointer with block of matrix
  public mat_unacquire_block !disassociate pointer
  public mat_get_nrow      !get nrow
  public mat_get_ncol      !get ncol
  public mat_mpi_bcast     !broadcast matrix to all slaves

  ! implemented here
  public mat_set_block     !copy array into block of matrix
  public mat_get_block     !copy block of matrix to array
  public mat_dot_block     !dot block of matrix with array
  public mat_print         !print formatted
  public mat_zero_like     !allocate in same shape and zero elms

contains

  !> \brief inserts a block of values into a matrix
  !> \author Bin Gao
  !> \date 2012-01-17
  !> \param minrow is the minimum row index
  !> \param maxrow is the maximum row index
  !> \param mincol is the minimum column index
  !> \param maxcol is the maximum column index
  !> \param values contains the block of values
  !> \param trans indicates if transposing the values
  !> \return A is the matrix
  subroutine mat_set_block(A, minrow, maxrow, mincol, maxcol, values, trans)
    type(matrix), intent(inout) :: A
    integer,      intent(in)    :: minrow, maxrow, mincol, maxcol
    real(8),      intent(in)    :: values(minrow : maxrow, mincol : maxcol)
    logical, optional, intent(in) :: trans
    logical p_trans
    real(8), pointer :: block(:,:)
    p_trans = .false.
    if (present(trans)) p_trans = trans
    if (p_trans) then
       block => mat_acquire_block(A, mincol, maxcol, minrow, maxrow)
       block(:,:) = transpose(values)
    else
       block => mat_acquire_block(A, minrow, maxrow, mincol, maxcol)
       block(:,:) = values
    end if
    call mat_unacquire_block(A, block)
  end subroutine

  !> \brief gets (copies) a block of values from (the real part of) a matrix
  !> \author Bin Gao
  !> \date 2012-01-17
  !> \param A is the matrix
  !> \param minrow is the minimum row index
  !> \param maxrow is the maximum row index
  !> \param mincol is the minimum column index
  !> \param maxcol is the maximum column index
  !> \return values contains the block of values
  subroutine mat_get_block(A, minrow, maxrow, mincol, maxcol, values)
    type(matrix)                :: A !intent in, but changed underway
    integer,      intent(in)    :: minrow, maxrow, mincol, maxcol
    real(8),      intent(out)   :: values(minrow : maxrow, mincol : maxcol)
    real(8), pointer :: block(:,:)
    block => mat_acquire_block(A, minrow, maxrow, mincol, maxcol)
    values(:,:) = block
    call mat_unacquire_block(A, block)
  end subroutine

  !> \brief gets the dot product between a block of values and
  !>        a corresponding block of (the real part of) a matrix
  function mat_dot_block(A, minrow, maxrow, mincol, maxcol, values, trans)
    type(matrix)                     :: A !intent in, but changed underway
    integer,           intent(in)    :: minrow, maxrow, mincol, maxcol
    real(8),           intent(in)    :: values(minrow : maxrow, mincol : maxcol)
    logical, optional, intent(in)    :: trans
    real(8) mat_dot_block
    logical p_trans !indicates if transposing the values (private)
    real(8), pointer :: block(:,:)
    call mat_assert_def(A, 'mat_dot_block called with matrix A undefined')
    p_trans = .false.
    if (present(trans)) p_trans = trans
    if (p_trans) then
       block => mat_acquire_block(A, mincol, maxcol, minrow, maxrow)
       ! factor two here due to closed-shell
       mat_dot_block = sum(block * transpose(values))
    else
       block => mat_acquire_block(A, minrow, maxrow, mincol, maxcol)
       ! factor two here due to closed-shell
       mat_dot_block = sum(block * values)
    end if
    call mat_unacquire_block(A, block)
    ! douple of closed-shell
    if (mat_is_closed_shell(A)) mat_dot_block = mat_dot_block * 2
  end function


  ! print matrix A to unit (default 6=stdout), optionally starting with 'label',
  ! optionally making each column 'width' wide, using optional separators
  ! (left brace, column separator, right brace, row separator)
  subroutine mat_print(A, label, unit, width, separators)
    type(matrix)                       :: A !intent in, but changed under way
    character(*), optional, intent(in) :: label
    integer,      optional, intent(in) :: unit, width
    character(4), optional, intent(in) :: separators
    integer      nrow, ncol, uni, wid, pre, dec, i, j, siz
    character(8) fmt
    character(4) spr
    type(matrix) Aa
    real(8), pointer :: ptr(:,:)
    call mat_assert_def(A, 'mat_print called with undefined matrix A')
    ! process optional argument unit, which defaults to stdout
    uni = 6
    if (present(unit)) uni = unit
    ! set pre to the largest number of digits to be printed
    ! before the decimal point (including any minus sign)
    nrow = mat_get_nrow(A)
    ncol = mat_get_ncol(A)
    if (mat_is_zero(A)) then
       pre = -huge(1)
    else if (mat_is_complex(A)) then
       pre =          pre_decimals_mat(mat_get_part(A, imag=.false.))
       pre = max(pre, pre_decimals_mat(mat_get_part(A, imag=.true.)))
    else
       pre = pre_decimals_mat(A)
    end if
    ! process optional width. Ensure the pre-decimals will fit
    wid = 11
    if (present(width)) wid = max(width, max(pre,2)+2) !max, to avoid *****
    ! set dec to the number of decimals to be printed
    dec = wid - max(pre,2) - 1
    ! argument label is optional. If present, print that
    if (present(label)) write (uni,'(a)') label
    ! process optional argument sep(arators), defaulting to spaces
    spr = '    '
    if (present(separators)) spr = separators
    ! create the format string to be used for each element
    write (fmt,'(a2,i2,a1,i2,a1)') '(f', wid, '.', dec, ')'
    ! call the printing routine
    if (mat_is_complex(A)) then
       write (uni,'(a)') '(real part)'
       call print_part(mat_get_part(A, imag=.false.))
       write (uni,'(a)') '(imaginary part)'
       call print_part(mat_get_part(A, imag=.true.))
    else
       call print_part(A)
    end if
    if (present(label)) write (uni,'()') !final blank line'
  contains
    ! number of pre-decimals in the printing of number
    integer function pre_decimals(r)
      real(8), intent(in) :: r
      ! treat nan and inf as zero
      if (.not.(r==r) .or. r > huge(r) .or. r < -huge(r)) then
         pre_decimals = -huge(1)
      else !take logarithm, add one for minus sign
         pre_decimals = ceiling(log10(abs(r))) + merge(1,0,r<0)
      end if
    end function
    ! number of pre-decimals for a (real) matrix
    integer function pre_decimals_mat(A)
      type(matrix) :: A
      real(8), pointer :: elms(:,:)
      integer i, j
      elms => mat_acquire_block(A, 1, nrow, 1, ncol)
      pre_decimals_mat = maxval((/((pre_decimals(elms(i,j)), i=1,nrow), j=1,ncol)/))
      call mat_unacquire_block(A, elms)
    end function
    ! per-part printing routine. 'mat' must be real
    subroutine print_part(mat)
      type(matrix) :: mat !intent in but changed underway
      character(ncol*(wid+1)+3) line
      integer i, j, l
      real(8), pointer :: elms(:,:)
      if (mat_is_zero(mat)) then
         write (uni,'(a)') 'zero matrix'
         return
      end if
      elms => mat_acquire_block(mat, 1, nrow, 1, ncol)
      do i = 1, nrow
         line(1:1) = merge(spr(1:1),' ',i==1)
         line(2:2) = spr(1:1)
         l = 3
         do j = 1, ncol
            write (line(l:l+wid-1), fmt) elms(i,j)
            l = l + wid
            if (j/=ncol) line(l:l) = spr(2:2)
            if (j/=ncol) l = l+1
         end do
         line(l:l) = spr(3:3)
         line(l+1:l+1) = merge(spr(3:3), spr(4:4), i==nrow)
         l = l + 2
         if (present(separators)) then
            write (uni,'(a)') line
         else
            write (uni,'(a)') line(2:len(line)-2)
         end if
      end do
      call mat_unacquire_block(mat, elms)
    end subroutine
  end subroutine


  ! B = 0*A, allocated and filled with zeros
  subroutine mat_zero_like(A, B)
    type(matrix), intent(in)    :: A
    type(matrix), intent(inout) :: B
    ! first set B to a non-allocated zero matrix in the shape of A
    call mat_axpy(0d0, A, .false., .true., B)
    ! get B allocated and filled with zeros
    call mat_ensure_alloc(B)
  end subroutine

end module
