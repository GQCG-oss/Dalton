! Copyright 2009-2012 Andreas J. Thorvaldsen
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.
!
! Module matrix_backend contains a slightly generalized interface
! to LSdalton's module matrix_operations.
!
! Only intended as hidden back-end for matrix_defop.
module matrix_backend

  use matrix_module
  use matrix_operations,    &
      only: mat_nullify,    &
            lsdalton_mat_init => mat_init, &
            lsdalton_mat_free => mat_free, &
            mat_mul,        &
            mat_zero,       &
            mat_daxpy,      &
            mat_scal,       &
            mat_trans,      &
            mat_trab,       &
            mat_dotproduct, &
            mat_init_magic_value
  implicit none

  public matrix
  public realk
  public quit
  public mat_nullify
  public mat_init
  public mat_dup
  public mat_move
  public mat_remove
  public mat_is_defined
  public mat_assert_def
  public mat_is_zero
  public mat_is_complex
  public mat_is_quaternion
  public mat_is_closed_shell
  public mat_is_open_shell
  public mat_is_alias
  public mat_get_nrow
  public mat_get_ncol
  public mat_get_part
  public mat_axpy
  public mat_gemm
  public mat_dot
  public mat_trace
  public mat_hide_temp
  public mat_unhide_temp
  public mat_ensure_alloc
  public mat_imag_times_mat
  public mat_acquire_block
  public mat_unacquire_block
  public mat_mpi_bcast
  public matrix_backend_debug
  private
  
  logical :: matrix_backend_debug = .false.
  integer, parameter :: magic_val_zero      =  2035728172
  integer, parameter :: magic_val_temp      = -1776576656
  !integer, parameter :: magic_val_cplx      =  1474967875
  !integer, parameter :: magic_val_zero_cplx = -2279863826
  !integer, parameter :: magic_val_temp_cplx = -983746554

contains


  ! Scale-add or scale-copy X to Y: Y (+)= a*X^tx.  If += and both X and Y are zero,
  ! ie 0 += a*0^tx, this is a no-op, except that shapes are verified and that Y
  ! will receive any extra component flags from X (complex, open_shell).
  subroutine mat_axpy(a, X, tx, ey, Y)
    real(realk),  intent(in)    :: a
    type(matrix), intent(in)    :: X
    logical,      intent(in)    :: tx, ey
    type(matrix), intent(inout) :: Y
    logical      equals, zero
    type(matrix) Yy
    equals = ey !modifiable copy
    zero = (a==0 .or. mat_is_zero(X))
    ! verify and prepare Y
    if (.not.prepare_inout_mat(merge(X%ncol, X%nrow, tx), &
                               merge(X%nrow, X%ncol, tx), &
                               zero, X%complex, equals, Y)) &
       call quit('mat_axpy called with matrices X^tx and Y with different shapes')
    ! skip, or relay to complex or real handler
    if (zero .and. .not.(X%complex .and. .not.Y%complex)) then
       ! all done
    else if (zero .and. (X%complex .and. .not.Y%complex)) then
       call mat_init(Yy, Y%nrow, Y%ncol, is_zero=.true.) !zero, real
       call mat_set_part(Yy, Y, imag=.true.)
    else if (X%complex .or. Y%complex) then
       ! Yr (+)= a Xr
       Yy = mat_get_part(Y, imag=.false.)
       call axpy_real(a, mat_get_part(X, imag=.false.), tx, equals, Yy)
       call mat_set_part(Yy, Y)
       ! alloc im part of Y if zero
       Yy = mat_get_part(Y, imag=.true.)
       if (mat_is_zero(Yy)) equals = .true.
       if (mat_is_zero(Yy)) &
          call mat_init(Yy, Yy%nrow, Yy%ncol)
       ! Yi (+)= a Xi
       call axpy_real(a, mat_get_part(X, imag=.true.), tx, equals, Yy)
       call mat_set_part(Yy, Y, imag=.true.)
       call mat_nullify(Yy)
    else
       call axpy_real(a, X, tx, equals, Y)
    end if
  end subroutine


  ! private per-part axpy
  subroutine axpy_real(a, X, tx, ey, Y)
    use matrix_operations, only: matrix_type, mtype_dense
    real(realk),  intent(in)    :: a
    type(matrix), intent(in)    :: X
    logical,      intent(in)    :: tx, ey
    type(matrix), intent(inout) :: Y
    type(matrix) Xt
    logical XisY, useXt
    ! determine whether X aliases some part of Y
    XisY = associated(X%init_self_ptr, Y%init_self_ptr)
    if (matrix_type == mtype_dense .and. XisY) &
       XisY = associated(X%elms, Y%elms)
    ! whether we must compute X^T prior to adding it to Y
    ! Y=a*Y^T Y+=a*X^t
    useXt = (tx .and. a/=0 .and. .not.mat_is_zero(X) &
             .and. (XisY .or. a/=1 .or. .not.ey))
    if (useXt) call mat_init(Xt, X%ncol, X%nrow)
    if (useXt) call mat_trans(X, Xt)
    ! Y=0, Y=a*Y, Y+=a*Y, Y=a*X^T
    if (.not.useXt .and. (mat_is_zero(X) .or. (tx .and. ey) .or. XisY)) then
       if (tx .and. ey) call mat_trans(X, Y)
       if ((a==0 .or. mat_is_zero(X)) .and. ey) then
          call mat_zero(Y)
       else if (.not.mat_is_zero(X)) then
          call mat_scal(merge(a,a+1,ey), Y)
       end if
    else !Y=a*X, Y+=a*X, Y+=a*X^T
       if (ey) call mat_zero(Y) !since no out-of-place scal exists
       call mat_daxpy(a, merge(Xt,X,useXt), Y)
       if (useXt) call mat_remove(Xt)
    end if
  end subroutine

    
  ! matrix multiply: C (+)= f * A^ta * B^tb. (Re-)Allocate C if it is
  ! an alias or has fewer parts than A and B
  subroutine mat_gemm(f, A, ta, B, tb, ec, C)
    real(realk),          intent(in)    :: f
    type(matrix), target, intent(in)    :: A, B
    logical,              intent(in)    :: ta, tb, ec
    type(matrix), target, intent(inout) :: C
    logical equals, zero, cplx
    equals = ec !modifiable copy
    zero = (f==0 .or. mat_is_zero(A) .or. mat_is_zero(B))
    cplx = (A%complex .or. B%complex)
    ! verify that A and B can be multiplied
    if (merge(A%nrow, A%ncol, ta) /= merge(B%ncol, B%nrow, tb)) &
       call quit('mat_gemm called with matrices A^ta and B^tb which cannot be multiplied')
    ! verify and prepare C
    if (.not.prepare_inout_mat(merge(A%ncol, A%nrow, ta), &
                               merge(B%nrow, B%ncol, tb), &
                               zero, cplx, equals, C)) &
       call quit('mat_gemm called with A^ta * B^tb which does not fit in C')
    ! skip, or relay to complex or real handler
    if (zero .and. .not.(cplx .and. .not.C%complex)) then
       ! all done
    else if (cplx .or. C%complex) then
       call gemm_complex
    else
       call gemm_real(f, A, B, equals, C)
    end if
  contains
    ! split into parts, multiply parts with do_real, join
    subroutine gemm_complex
      type(matrix) Ar, Ai, Br, Bi, Cc
      Ar = mat_get_part(A, imag=.false.)
      Ai = mat_get_part(A, imag=.true.)
      Br = mat_get_part(B, imag=.false.)
      Bi = mat_get_part(B, imag=.true.)
      ! Cr (+)= f Ar Br
      Cc = mat_get_part(C, imag=.false.)
      call gemm_real(f, Ar, Br, equals, Cc)
      ! Cr   -= f Ai Bi
      if (.not.(mat_is_zero(Ai) .or. mat_is_zero(Bi))) &
         call gemm_real(-f, Ai, Bi, .false., Cc)
      call mat_set_part(Cc, C, imag=.false.)
      Cc = mat_get_part(C, imag=.true.)
      ! if Cc was real, the imaginary parts needs to be allocated
      if (mat_is_zero(Cc)) then
         call mat_init(Cc, Cc%nrow, Cc%ncol)
         equals = .true.
      end if
      ! Ci (+)= f Ai Br
      if (.not.mat_is_zero(Bi)) &
         call gemm_real(f, Ar, Bi, equals, Cc)
      ! Ci   += f Ar Bi
      if (.not.mat_is_zero(Ai)) &
         call gemm_real(f, Ai, Br, equals .and. mat_is_zero(Bi), Cc)
      call mat_set_part(Cc, C, imag=.true.)
      call mat_nullify(Ar)
      call mat_nullify(Ai)
      call mat_nullify(Br)
      call mat_nullify(Bi)
    end subroutine
    subroutine gemm_real(f, A, B, ec, C)
      real(realk),  intent(in)    :: f
      type(matrix), intent(in)    :: A, B
      logical,      intent(in)    :: ec
      type(matrix), intent(inout) :: C
      call mat_mul(A, B, merge('T','N',ta), merge('T','N',tb), &
                   f, merge(0E0_realk, 1E0_realk, ec), C)
    end subroutine
  end subroutine


  subroutine mat_init(A, nrow, ncol, is_zero, is_complex, is_quaternion, is_open_shell)
    type(matrix), target, intent(out) :: A
    integer,              intent(in)  :: nrow, ncol
    logical,    optional, intent(in)  :: is_zero, is_complex, is_open_shell
    logical,    optional, intent(in)  :: is_quaternion
    logical zero, complex
    if (present(is_quaternion)) call quit( &
          "mat_init called with 'is_quaternion' present, but this isn't implemented")
    if (present(is_open_shell)) call quit( &
          "mat_init called with 'is_open_shell' present, but this isn't implemented")
    zero = .false.
    if (present(is_zero)) zero = is_zero
    complex = .false.
    if (present(is_complex)) complex = is_complex
    call mat_nullify(A)
    if (zero) then
       A%complex = complex
       A%nrow = nrow
       A%ncol = ncol
       A%init_magic_tag = magic_val_zero
    else
       call lsdalton_mat_init(A, nrow, ncol, complex=complex)
       !A%elms=0.0/0.0 !NaN to detect use of uninitalized elms (ifort only)
    end if
  end subroutine

  ! Deallocate A if non-alias, nullify if alias
  subroutine mat_remove(A)
    type(matrix), target,intent(inout) :: A
    call mat_assert_def(A, 'mat_remove(A) called with matrix A undefined')
    if (mat_is_zero(A) .or. .not.mat_is_alias(A, A)) then
       call mat_nullify(A)
    else
       call lsdalton_mat_free(A)
    end if
  end subroutine

  ! 'safely' move matrix A to B. If B is allocated, and A isn't its alias,
  ! free B first. Then overwrite B by A (duplicate)
  subroutine mat_move(A, B)
    type(matrix), target, intent(inout) :: A, B
    call mat_assert_def(A, 'mat_move(A,B) called with matrix A undefined')
    if (mat_is_defined(B) .and. .not.mat_is_alias(A, B)) &
       call mat_remove(B)
    call mat_dup(A, B)
    if (mat_is_alias(A, A)) B%init_self_ptr => B
    call mat_nullify(A)
  end subroutine

  ! fields-wise copy A to B. No book-keeping involved
  subroutine mat_dup(A, B)
    type(matrix), target, intent(in)  :: A
    type(matrix), target, intent(out) :: B
    B = A
  end subroutine

  ! whether A is defined
  logical function mat_is_defined(A)
    type(matrix), intent(in) :: A
    mat_is_defined = (A%init_magic_tag == mat_init_magic_value &
                 .or. A%init_magic_tag == magic_val_zero)
  end function

  ! quit if A isn't defined
  subroutine mat_assert_def(A, errmsg)
    type(matrix), target, intent(in) :: A
    character(*),         intent(in) :: errmsg
    if (.not.mat_is_defined(A)) call quit(errmsg)
  end subroutine

  ! whether A is an unallocated zero matrix
  logical function mat_is_zero(A)
    type(matrix), intent(in) :: A
    mat_is_zero = (A%init_magic_tag == magic_val_zero)
  end function

  ! whether A is a complex matrix (Ar Ai)
  logical function mat_is_complex(A)
    type(matrix), intent(in) :: A
    mat_is_complex = A%complex
  end function

  ! whether A is a quaternion matrix (A0 Ai Aj Ak)
  logical function mat_is_quaternion(A)
    type(matrix), intent(in) :: A
    mat_is_quaternion = .false.
  end function

  ! whether A is a closed shell matrix, for which a dot product
  ! carries an addional factor 2
  logical function mat_is_closed_shell(A)
    type(matrix), intent(in) :: A
    mat_is_closed_shell = .true. !FIXME implement non-shell and open-shell
  end function

  ! whether A is an open-shell matrix (Aa Ab). The alpha component will
  ! carries an addional factor 2
  logical function mat_is_open_shell(A)
    type(matrix), intent(in) :: A
    mat_is_open_shell = .false. !FIXME implement open-shell
  end function

  ! whether A aliases B, or B aliases A, or A and B alias the same matrix
  logical function mat_is_alias(A, B)
    type(matrix), target, intent(in) :: A, B
    mat_is_alias = associated(A%init_self_ptr, B)
  end function

  integer function mat_get_nrow(A)
    type(matrix), intent(in) :: A
    mat_get_nrow = A%nrow
  end function

  integer function mat_get_ncol(A)
    type(matrix), intent(in) :: A
    mat_get_ncol = A%ncol
  end function


  ! matrix dot product: sum_ij Aij Bij (ta=F)
  !  or  product trace: sum_ij Aji Bij (ta=F)
  ! A and B are assumed initialized and with equal/transpose shapes
  function mat_dot(A, B, t)
    type(matrix), intent(in) :: A, B
    logical,      intent(in) :: t
    complex(realk) mat_dot
    call mat_assert_def(A, 'mat_dot(A,B) or mat_trace(A,B) called with matrix A undefined')
    call mat_assert_def(B, 'mat_dot(A,B) or mat_trace(A,B) called with matrix B undefined')
    mat_dot = blk_dot(1,1) &
            + (0E0_realk,1E0_realk) * merge(1,-1,t) * blk_dot(2,1) &
            + (0E0_realk,1E0_realk) * blk_dot(1,2) &
            + merge(-1,1,t) * blk_dot(2,2)
    if (mat_is_closed_shell(A) .or. mat_is_closed_shell(B)) &
       mat_dot = 2*mat_dot
  contains
    function blk_dot(ia, ib)
      integer      :: ia, ib
      real(realk)  :: blk_dot
      type(matrix) :: dupA, dupB
      logical      :: Xis0
      integer      :: siz
      blk_dot = 0
      if (ia/=1 .and. .not.A%complex) return
      if (ib/=1 .and. .not.B%complex) return
      if (A%complex) then
         dupA = A !copy all fields
         dupA%complex = .false.
         siz = size(A%elms)/2
         dupA%elms => A%elms(1 + (ia-1)*siz : ia*siz)
      end if
      if (B%complex) then
         dupB = B !copy all fields
         dupB%complex = .false.
         siz = size(B%elms)/2
         dupB%elms => B%elms(1 + (ib-1)*siz : ib*siz)
      end if
      if (t) then
         blk_dot = mat_trab(merge(dupA, A, A%complex), &
                            merge(dupB, B, B%complex))
      else
         blk_dot = mat_dotproduct(merge(dupA, A, A%complex), &
                                  merge(dupB, B, B%complex))
      end if
    end function
  end function


  ! matrix trace = sum Aii, for A square
  function mat_trace(A)
    implicit none
    type(matrix), intent(in) :: A
    real(realk) :: mat_trace
    integer :: i
    i=0
    call mat_assert_def(A, 'mat_trace(A) called with matrix A undefined')
    if (A%nrow/=A%ncol) call quit('mat_trace called with non-square matrix A')
    mat_trace = 0
    if (.not.mat_is_zero(A)) &
       mat_trace = sum((/(A%elms(1+(i-1)*(A%nrow+1)), i=1,A%nrow)/))
    if (mat_is_closed_shell(A)) mat_trace = 2*mat_trace
  end function

  ! return alias to either real or imaginary part of A
  function mat_get_part(A, imag, quat_i, quat_j, beta) result(B)
    type(matrix),      intent(in) :: A
    logical, optional, intent(in) :: imag, quat_i, quat_j, beta
    type(matrix) B
    logical im
    if (present(quat_i)) call quit( &
          "mat_get_part called with 'quat_i' present, which isn't implemented yet")
    if (present(quat_j)) call quit( &
          "mat_get_part called with 'quat_j' present, which isn't implemented yet")
    if (present(beta)) call quit( &
          "mat_get_part called with 'beta' present, which isn't implemented yet")
    im = .false.        
    if (present(imag)) im = imag
    if (mat_is_zero(A) .or. (im .and. .not.A%complex)) then
       call mat_init(B, A%nrow, A%ncol, is_zero=.true.)
    else
       B = A
       if (A%complex .and. .not.im) then
          B%elms => A%elms(:A%nrow*A%ncol)
       else if (A%complex .and. im) then
          B%elms => A%elms(A%nrow*A%ncol+1 : 2*A%nrow*A%ncol)
       end if
    end if
    B%complex = .false.
  end function

  ! private used by mat_axpy and mat_gemm to verify and prepare
  ! their inout args Y and C, respectively
  function prepare_inout_mat(nrow, ncol, zero, cplx, equals, A)
    integer,      intent(in)    :: nrow, ncol
    logical,      intent(in)    :: zero, cplx
    logical,      intent(inout) :: equals
    type(matrix), intent(inout) :: A
    logical      prepare_inout_mat
    type(matrix) reA
    ! verify current contents, if any
    if (mat_is_defined(A)) then
       prepare_inout_mat = (nrow == A%nrow .and. ncol == A%ncol)
       if (.not.equals .and. .not.prepare_inout_mat) return
       ! free if it is A=..., but current A has wrong shape, is non-zero when the result should
       ! be zero, is complex when the result is real, or is an alias
       if (equals .and. (.not.prepare_inout_mat &
                    .or. (zero .and. .not.mat_is_zero(A)) &
                    .or. (A%complex .and. .not.cplx) &
                    .or. (.not.zero .and. .not.mat_is_alias(A, A)))) &
          call mat_remove(A)
    end if
    ! initialize if undefined
    if (.not.mat_is_defined(A)) then
       call mat_init(A, nrow, ncol, is_zero=zero, is_complex=cplx)
    ! otherwise, no need to reallocate if 0(+)=0, A non-alias, or Az+=0z
    else if ((zero .and. mat_is_zero(A)) .or. mat_is_alias(A, A) &
        .or. (zero .and. .not.(cplx .and. .not.A%complex))) then
        if (mat_is_zero(A)) A%complex = (A%complex .or. cplx)
    ! reallocate if A alias or real when result complex
    else 
       if (mat_is_zero(A)) equals = .true.
       call mat_ensure_alloc(A, only_alloc=.true.)
    end if
    prepare_inout_mat = .true.
  end function

  ! If the matrix is zero, allocate and zero elements. If it is an alias,
  ! allocate and copy. Otherwise it is allocated, so do nothing
  ! only_alloc=T
  subroutine mat_ensure_alloc(A, only_alloc)
    type(matrix),      intent(inout) :: A
    logical, optional, intent(in)    :: only_alloc
    type(matrix) B, Bb
    logical      no
    real(realk)  fac
    call mat_assert_def(A, 'mat_ensure_alloc(A) called with matrix A undefined')
    ! if non-alias, return
    if (mat_is_alias(A, A)) return
    ! allocate B
    call mat_init(B, A%nrow, A%ncol, is_complex=A%complex)
    ! decide whether or not to set the data in the newly allocated B (no)
    no = .false.
    if (present(only_alloc)) no = (only_alloc .and. mat_is_zero(A))
    ! decide which factor to use
    fac = merge(0E0_realk, 1E0_realk, mat_is_zero(A))
    ! make the copy (or zero). For complex do the parts separately
    if (.not.no .and. A%complex) then
       Bb = mat_get_part(B, imag=.false.)
       call axpy_real(fac, mat_get_part(A, imag=.false.), .false., .true., Bb)
       call mat_set_part(Bb, B, imag=.false.)
       Bb = mat_get_part(B, imag=.true.)
       call axpy_real(fac, mat_get_part(A, imag=.true.), .false., .true., Bb)
       call mat_set_part(Bb, B, imag=.true.)
    else if (.not.no) then
       call axpy_real(fac, A, .false., .true., B)
    end if
    ! finally move the copy in B over A
    call mat_move(B, A)
  end subroutine

  ! place i*f*A in B. This is the only way to construct a complex
  ! matrix when you only have real ones.
  subroutine mat_imag_times_mat(f, A, B)
    real(realk),          intent(in)    :: f
    type(matrix), target, intent(in)    :: A
    type(matrix), target, intent(inout) :: B
    logical      ignored
    type(matrix) Bb
    ! prepare B so that it can contain i*f*A, return if zero
    ignored = .true.
    ignored = prepare_inout_mat(A%nrow, A%ncol, &
                    f==0 .or. mat_is_zero(A), .true., ignored, B)
    if (mat_is_zero(B)) return
    ! re(B) = -f*im(A)
    Bb = mat_get_part(B, imag=.false.)
    if (A%complex) then
       call axpy_real(-f, mat_get_part(A, imag=.true.), .false., .true., Bb)
    else
       call axpy_real(0E0_realk, A, .false., .true., Bb)
    end if
    call mat_set_part(Bb, B, imag=.false.)
    ! im(B) = f*re*(A)
    Bb = mat_get_part(B, imag=.true.)
    call axpy_real(f, mat_get_part(A, imag=.false.), .false., .true., Bb)
    call mat_set_part(Bb, B, imag=.true.)
  end subroutine


  ! re(B)=A or im(B)=A. A must be real
  subroutine mat_set_part(A, B, imag, quat_i, quat_j, beta)
    use matrix_operations, only: matrix_type, mtype_dense
    type(matrix), target, intent(inout) :: A, B
    logical,    optional, intent(in)    :: imag, quat_i, quat_j, beta
    type(matrix), target :: C
    integer n, i
    logical im, useC, Bzero, Bcplx
    if (present(quat_i)) call quit( &
          "mat_set_part called with 'quat_i' present, which isn't implemented yet")
    if (present(quat_j)) call quit( &
          "mat_set_part called with 'quat_j' present, which isn't implemented yet")
    if (present(beta)) call quit( &
          "mat_set_part called with 'beta' present, which isn't implemented yet")
    im = .false.
    if (present(imag)) im = imag
    ! verify that we don't attempt to set the imaginary part of
    ! a matrix for which complex isn't even implemented
    if (im .and. matrix_type /= mtype_dense) &
       call quit("mat_set_part called with a matrix type for which complex isn't implemented")
    ! A must be real
    if (A%complex) call quit('mat_set_part(A,B) called with non-real matrix A')
    ! if B will become complex, reallocate it, and copy/zero the other part
    n = A%nrow * A%ncol
    Bzero = mat_is_zero(B)
    Bcplx = B%complex
    if ((mat_is_zero(A) .and. Bzero) .or. (im .and. .not.B%complex)) then
       if (.not.Bzero) call mat_move(B, C)
       call mat_init(B, A%nrow, A%ncol, is_complex=(im.or.Bcplx))
       i = merge(0,n,im)
       if (.not.Bzero) B%elms(i+1:i+n) = C%elms(:)
       if (Bzero .and. B%complex) B%elms(i+1:i+n) = 0
       if (.not.Bzero) call mat_remove(C)
    end if
    ! either copy from A, or zero
    if (mat_is_zero(A) .or. B%complex) then
       i = merge(n,0,im)
       if (.not.mat_is_zero(A)) B%elms(i+1:i+n) = A%elms(:)
       if (     mat_is_zero(A)) B%elms(i+1:i+n) = 0
       call mat_remove(A)
    else !or move the entire A
       call mat_move(A, B)
    end if
  end subroutine


  subroutine mat_hide_temp(A)
    type(matrix), intent(inout) :: A
    if (mat_is_alias(A, A)) &
       A%init_magic_tag = magic_val_temp
  end subroutine

  subroutine mat_unhide_temp(A)
    type(matrix), target, intent(inout) :: A
    if (A%init_magic_tag == magic_val_temp) then
       A%init_magic_tag = mat_init_magic_value
       A%init_self_ptr => A
    end if
  end subroutine

  function mat_acquire_block(A, minrow, maxrow, mincol, maxcol)
    type(matrix), target :: A !no intent to allow use in intent(in) routines
    integer,  intent(in) :: minrow, maxrow, mincol, maxcol
    real(realk), pointer :: mat_acquire_block(:,:)
    if (A%init_magic_tag == magic_val_temp) then
       call quit('mat_acquire_block called with matrix A already acquired')
    else if (A%init_magic_tag == magic_val_zero) then
       call quit('mat_acquire_block called with matrix A zero')
    else if (A%init_magic_tag /= mat_init_magic_value) then
       call quit('mat_acquire_block called with matrix A undefined')
    end if
    mat_acquire_block => pointer_to_block(A%elms(:A%nrow*A%ncol))
    A%init_magic_tag = magic_val_temp
  contains
    function pointer_to_block(elms)
      real(realk), target, intent(in) :: elms(A%nrow, A%ncol)
      real(realk), pointer ::  pointer_to_block(:,:)
      pointer_to_block => elms(minrow : maxrow, mincol : maxcol)
    end function
  end function

  subroutine mat_unacquire_block(A, block)
    type(matrix), target :: A !no intent to allow use in intent(in) routines
    real(realk),  pointer :: block(:,:) !intent(inout), but pointer can't have intent
    if (A%init_magic_tag /= magic_val_temp) &
       call quit('mat_unacquire_block called with matrix A not already acquired')
    nullify(block)
    A%init_magic_tag = mat_init_magic_value
  end subroutine

  !> \brief broadcasts matrix
  !> \author Bin Gao
  !> \date 2012-05-13
  !> \param A is the matrix
  !> \param root is the root processor which broadcasts the matrix
  !> \param mat_comm is the MPI communicator
  subroutine mat_mpi_bcast(A, root, mat_comm)
    type(matrix), intent(inout) :: A
    integer,      intent(in)    :: root
    integer,      intent(in)    :: mat_comm
  end subroutine

  ! LSdalton uses lsquit(msg,lupri) instead of quit, so make this proxy.
  ! Only for use here, in matrix_defop, and in matrix_lowlevel
  subroutine quit(msg)
    character(*), intent(in) :: msg
    call lsquit(msg,-1) !-1 means lupri unavailable
  end subroutine

end module
