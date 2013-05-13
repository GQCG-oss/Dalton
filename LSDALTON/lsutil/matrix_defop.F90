! Copyright 2009-2012 Andreas J. Thorvaldsen
! This file is made available under the terms of the
! GNU Lesser General Public License version 3.

! This file is generic for both DALTON and DIRAC.
! DALTON or DIRAC specific things are in matrix_backend
! which is the backend to matrix_defop.
!
! Defined operators for type(matrix)
!
! The operators are as follows:
!
!    equals                       A = B
!    plus                         A + B
!    minus                        A - B
!    matrix multiply              A * B
!    scale by integer             n * A
!    scale by real(realk)         r * A
!    scale by complex(realk)      z * A
!    transpose                    trps(A)
!    trace                        trace(A)
!    scalar product               dot(A,B)
!    product trace                trace(A,B)
!    short-hand freeing           A=0, A(:)=0, A(:,:)=0 ...
!    norm ie. sqrt(dot(A,A))      norm(A)
!
! To increase performance and lower memory needs, formulas are not
! evaluated immediately as they stand, but the operators are
! first collected into 'proxy' matrices of either of the forms
!
!    axpy:      f * A(^T)
!    gemm:      f * A(^T) * B(^T)
!
! To check the state of matrices: Whether defined: isdef(A),
! whether zero: iszero(A) (zero matrices do not allocate memory)
!
! NB: matrix_defop strives to follow fortran 90 standard, but if you run a rigorous
!     checker (eg. 'forcheck'), you will find numerous instances of the following errors:
!
!    **[539 E] procedure must have private accessibility (because arg type private)
!    **[620 E] dummy input argument must not be defined  (because proxy arg intent in)
!
! This is because we have made type(proxy), which is used during formula evaluation,
! private, to prevent users from declaring and programming with type(proxy)
! in their code. We also destroy (clear fields of) each type(proxy) after use, to
! avoid littering memory with the type(matrix) contained within them. If your compiler
! chokes on these even with 'strict' mode turned off, it might be because it wishes
! to optimize by reusing type(proxy), which would anyway break matrix_defop.
! In that case you probably won't be able to use matrix_defop.
!
module matrix_defop

  use matrix_backend, &
           isdef => mat_is_defined, &
           iszero => mat_is_zero, &
           isalias => mat_is_alias

  implicit none

  public matrix
  public isdef
  public iszero
  public assignment(=)
  public operator(+)
  public operator(-)
  public operator(*)
  public trans
  public dot
  public norm
  public trace
  public mat_move
  public mat_ensure_alloc
  public mat_get_part
  public matrix_defop_debug
  
  ! uncomment if compiler complains about type(proxy) being private
  ! public proxy

  ! for switching debugging on or off
  logical :: matrix_defop_debug = .false.

  interface assignment(=)
     module procedure mat_eq_mat
     module procedure mat_eq_prx
     module procedure mat_eq_zero
     module procedure mat_eq_zero_1D
     module procedure mat_eq_zero_2D
     module procedure mat_eq_zero_3D
     module procedure mat_eq_zero_4D
  end interface

  interface operator(+)
     module procedure mat_plus_mat
     module procedure mat_plus_prx
     module procedure prx_plus_mat
     module procedure prx_plus_prx
     module procedure plus_mat
     module procedure plus_prx
  end interface

  interface operator(-)
     module procedure mat_minus_mat
     module procedure mat_minus_prx
     module procedure prx_minus_mat
     module procedure prx_minus_prx
     module procedure minus_mat
     module procedure minus_prx
  end interface

  interface operator(*)
     module procedure real_times_mat
     module procedure real_times_prx
     module procedure complex_times_mat
     module procedure complex_times_prx
     module procedure integer_times_mat
     module procedure integer_times_prx
     module procedure mat_times_mat
     module procedure mat_times_prx
     module procedure prx_times_mat
     module procedure prx_times_prx
  end interface

  interface trans
     module procedure transpose_mat
     module procedure transpose_prx
  end interface

  interface dot
     module procedure mat_dot_mat
     module procedure mat_dot_prx
     !ajt No idea why gfortran complains about this
     ! module procedure prx_dot_mat
     module procedure prx_dot_prx
  end interface

  interface trace
     module procedure mat_prodtr_mat
     module procedure mat_prodtr_prx
     !ajt No idea why gfortran complains about this
     ! module procedure prx_prodtr_mat
     module procedure prx_prodtr_prx
     module procedure mat_trace
     module procedure prx_trace
  end interface

  interface norm
     module procedure mat_norm
     module procedure prx_norm
  end interface

  ! private type used to contain intermediates during matrix
  ! algebra evaluation
  type proxy
     ! private
     real(realk)  f          !scale factor. A must be zero when f=0
     logical      ta, hb, tb !A(^T), have B, B(^T)
     type(matrix) A, B       !matrices in f*A*B
     logical      is_delayed !whether this proxy is delayed, and misses
                             !a contribution kept in 'delayed_matrix'
  end type

  ! to avoid allocating a temporary during evaluation of expressions of the
  ! form "C = f*A(*B) + C", the evaluation of mat_plus_mat, prx_plus_mat and mat_plus_prx
  ! is delayed to see whether a mat_eq_prx comes after it, in which case
  ! no temporary copy of C is needed
  logical               :: have_delayed_matrix = .false.
  !logical               :: delayed_minus
  type(matrix), pointer :: delayed_matrix

  private

contains

  ! copy A to B. If B's current allocation cannot fit A, allocate B
  ! in the shape of A, then copy the data. If A is zero, deallocate B's
  ! current allocation (if any)
  subroutine mat_eq_mat(B, A)
    type(matrix), target, intent(inout) :: B
    type(matrix), target, intent(in)    :: A
    call mat_axpy(1E0_realk, A, .false., .true., B)
  end subroutine

  ! evaluate matrix proxy P and place the result in M
  subroutine mat_eq_prx(M, P)
    type(matrix), target, intent(inout) :: M
    type(proxy),  target, intent(in)    :: P
    call mat_eq_prx_subr(M, P)
  end subroutine

  ! evaluate matrix proxy P and place the result in M. If P is
  ! delayed, and M is the delayed_matrix, call mat_revplusmin_prx
  ! to have it evaluated first
  subroutine mat_eq_prx_subr(M, P)
    type(matrix), target, intent(inout) :: M
    type(proxy),  target                :: P
    logical      Mbusy, equals
    type(matrix) newM
    ! if P has M as delayed contribution: P = f*A(*B) + M, this is really
    ! M += f*A(B), so no need for scratch
    equals = .true.
    if (have_delayed_matrix .and. associated(delayed_matrix, M)) then
       if (.not.P%is_delayed) call quit( &
             'have_delayed_matrix .and. .not.P%is_delayed should be impossible')
       nullify(delayed_matrix)
       have_delayed_matrix = .false.
       P%is_delayed = .false.
       equals = .false.
    end if
    call rectify_input_proxy(P)
    ! if P=1*A or 0, move P%A to M, and don't evaluate
    if (equals .and. (P%f==0 .or. P%f==1) .and. .not.P%ta &
               .and. .not.P%hb .and. isalias(P%A, P%A)) then
       call mat_move(P%A, M)
       return
    end if
    ! determine whether M is involved in the operation, thus can't be
    ! overwritten during it
    if (.not.P%hb) Mbusy = (isalias(P%A, M) .and. P%ta)
    if (     P%hb) Mbusy = (isalias(P%A, M) .or. isalias(P%B, M))
    if (.not.equals .and. Mbusy) call mat_dup(M, newM)
    ! call axpy or gemm
    if (.not.P%hb .and. Mbusy) then
       call mat_axpy(P%f, P%A, P%ta, equals, newM)
    else if (.not.P%hb .and. .not.Mbusy) then
       call mat_axpy(P%f, P%A, P%ta, equals, M)
    else if (P%hb .and. Mbusy) then
       call mat_gemm(P%f, P%A, P%ta, P%B, P%tb, equals, newM)
    else if (P%hb .and. .not.Mbusy) then
       call mat_gemm(P%f, P%A, P%ta, P%B, P%tb, equals, M)
    end if
    if (Mbusy) call mat_move(newM, M)
    call mat_remove(P%A)
    if (P%hb) call mat_remove(P%B)
  end subroutine


  !------------------------------------------
  !   add and subtract operators A+B, A-B
  !------------------------------------------
  
  function mat_plus_mat(A, B) result(P)
    type(matrix), target, intent(in) :: A, B
    type(proxy),  target             :: P
    if (.not.have_delayed_matrix .and. &
              .not.isalias(A, A) .and. isalias(B, B)) then
       call mat_revplusmin_mat(.false., .false., B, A, P)
    else
       call mat_revplusmin_mat(.false., .false., A, B, P)
    end if
  end function

  function mat_minus_mat(A, B) result(P)
    type(matrix), target, intent(in) :: A, B
    type(proxy),  target             :: P
    if (.not.have_delayed_matrix .and. &
              .not.isalias(A, A) .and. isalias(B, B)) then
       call mat_revplusmin_mat(.true., .true., B, A, P)
    else
       call mat_revplusmin_mat(.false., .true., A, B, P)
    end if
  end function

  function mat_plus_prx(A, P) result(R)
    type(matrix), target, intent(in) :: A
    type(proxy),  target, intent(in) :: P
    type(proxy),  target             :: R
    call mat_revplusmin_prx(.false., .true., .false., A, P, R)
  end function

  function mat_minus_prx(A, P) result(R)
    type(matrix), target, intent(in) :: A
    type(proxy),  target, intent(in) :: P
    type(proxy),  target             :: R
    call mat_revplusmin_prx(.false., .false., .true., A, P, R)
  end function

  function prx_plus_mat(P, A) result(R)
    type(proxy),  target, intent(in) :: P
    type(matrix), target, intent(in) :: A
    type(proxy),  target             :: R
    call mat_revplusmin_prx(.false., .false., .false., A, P, R)
  end function

  function prx_minus_mat(P, A) result(R)
    type(proxy),  target, intent(in) :: P
    type(matrix), target, intent(in) :: A
    type(proxy),  target             :: R
    call mat_revplusmin_prx(.false., .true., .true., A, P, R)
  end function

  function prx_plus_prx(P, Q) result(R)
    type(proxy), target, intent(in) :: P, Q
    type(proxy), target             :: R
    call prx_revplusmin_prx(.false., .false., P, Q, R)
  end function

  function prx_minus_prx(P, Q) result(R)
    type(proxy), target, intent(in) :: P, Q
    type(proxy), target             :: R
    logical pq, qp, havePA, haveQA
    call prx_revplusmin_prx(.false., .true., P, Q, R)
  end function

  ! A+B, A-B or B-A. If there is no 'delayed' type(matrix),
  ! and A isn't alias, place that in delayed_matrix, and only return P=+-B.
  ! If not, add/subtract A and B using axpy.
  subroutine mat_revplusmin_mat(rev, minus, A, B, P)
    logical,              intent(in)  :: rev, minus
    type(matrix), target, intent(in)  :: A, B
    type(proxy),  target, intent(out) :: P
    call clear_proxy(P)
    call mat_dup(B, P%A)
    P%f = merge(-1, 1, minus .and. .not.rev)
    if (.not.have_delayed_matrix .and. isalias(A, A) &
                                 .and. .not.(rev .and. minus)) then
       delayed_matrix => A
       !delayed_minus  = (minus .and. rev)
       have_delayed_matrix = .true.
       P%is_delayed = .true.
    else
       call mat_axpy(1E0_realk * merge(-1, 1, minus), A, .false., .false., P%A)
    end if
    call prepare_result_proxy(P)
  end subroutine


  ! R=M+P, R=M-P or R=P-M.
  subroutine mat_revplusmin_prx(Mtemp, rev, minus, M, P, R)
    logical,              intent(in)  :: Mtemp, rev, minus
    type(matrix), target              :: M
    type(proxy),  target              :: P
    type(proxy),  target, intent(out) :: R
    logical minM, Atemp
    ! start by moving P to R
    R = P
    call mat_nullify(P%A)
    if (P%hb) call mat_nullify(P%B)
    call rectify_input_proxy(R)
    ! apply and forget minus signs, if possible
    minM = (minus .and. rev .and. .not.iszero(M))
    if (minus .and. .not.rev) R%f = -R%f
    if (iszero(M)) minM = .false.
    ! 
    Atemp = (R%f==1 .and. .not.R%ta &
                    .and. .not.R%hb .and. isalias(R%A, R%A))
    if (.not.have_delayed_matrix .and. .not.Mtemp .and. isalias(M, M) &
                                 .and. .not.Atemp .and. .not.minM) then
       delayed_matrix => M
       have_delayed_matrix = .true.
       !delayed_minus  =  minM
       R%is_delayed = .true.
       call prepare_result_proxy(R)
       return
    end if
    !
    if (Mtemp .and. .not.(iszero(M) .and. .not.R%hb)) then
       if (minM) R%f = -R%f
       if (.not.R%hb) call mat_axpy(R%f, R%A, R%ta, .false., M)
       if (     R%hb) call mat_gemm(R%f, R%A, R%ta, R%B, R%tb, .false., M)
       R%f = merge(-1, 1, minM)
       call mat_move(M, R%A)
       R%ta = .false.
       if (R%hb) call mat_remove(R%B)
       R%hb = .false.
    else if (.not.R%hb .and. (Mtemp .or. isalias(R%A, R%A))) then
       ! apply f on A first if small, to avoid potential overflow
       if (.not.iszero(M) .and. abs(R%f) < 1/2E0_realk**16) then
          call mat_axpy(R%f, R%A, .false., .true., R%A)
          R%f = 1
       end if
       ! merge: (-1 if minM else 1) / (1 if M=0 else R%f)
       call mat_axpy(merge(-1,1,minM) / merge(1E0_realk,R%f,iszero(M)), &
                     M, R%ta, .false., R%A)
    else
       if (.not.Mtemp) call mat_dup(M, P%A)
       if (     Mtemp) call mat_move(M, P%A)
       if (minM) R%f = -R%f
       if (.not.R%hb) call mat_axpy(R%f, R%A, R%ta, .false., P%A)
       if (     R%hb) call mat_gemm(R%f, R%A, R%ta, R%B, R%tb, .false., P%A)
       R%f = merge(-1, 1, minM)
       call mat_move(P%A, R%A)
       R%ta = .false.
       if (R%hb) call mat_remove(R%B)
       R%hb = .false.
    end if
    call prepare_result_proxy(R)
  end subroutine


  ! R = P+Q, P-Q or -P+Q, where P, Q, R are matrix proxies.
  subroutine prx_revplusmin_prx(rev, minus, P, Q, R)
    logical,             intent(in)  :: rev, minus
    type(proxy), target              :: P, Q
    type(proxy), target, intent(out) :: R
    call rectify_input_proxy(P)
    call rectify_input_proxy(Q)
    if (minus .and. .not.rev) Q%f = -Q%f
    if (minus .and.      rev) P%f = -P%f
    if ((P%ta .or. P%hb) .and. (Q%ta .or. Q%hb)) then
        if (.not.P%hb) call mat_axpy(P%f, P%A, P%ta, .true., R%A)
        if (     P%hb) call mat_gemm(P%f, P%A, P%ta, P%B, P%tb, .true., R%A)
        call mat_move(R%A, P%A)
        if (P%hb) call mat_remove(P%B)
        P%f = 1; P%ta = .false.; P%hb = .false.
    end if
    ! ajt FIXME Factor transfer shouldn't be done with very small factors.
    !           Then the factor should be applied first
    if (.not.(P%ta .or. P%hb) .or. (Q%ta .or. Q%hb)) then
       if (P%f/=-1 .and. P%f/=0) Q%f = Q%f / P%f
       call prepare_result_proxy(Q)
       call mat_revplusmin_prx(.true., .true., P%f==-1, P%A, Q, R)
       if (P%f/=-1 .and. P%f/=0) R%f = R%f * P%f
    else !.not.Q%ta .and. .not.Q%hb
       if (Q%f/=-1 .and. Q%f/=0) P%f = P%f / Q%f
       call prepare_result_proxy(P)
       call mat_revplusmin_prx(.true., .true., Q%f==-1, Q%A, P, R)
       if (Q%f/=-1 .and. Q%f/=0) R%f = R%f * Q%f
    end if
  end subroutine


  !----------------------------------
  !   scale operators f*A, +A, -A
  !----------------------------------

  function plus_mat(A) result(P)
    type(matrix), target, intent(in) :: A
    type(proxy),  target             :: P
    call real_times_mat_subr(1E0_realk, A, P)
  end function

  function minus_mat(A) result(P)
    type(matrix), target, intent(in) :: A
    type(proxy),  target             :: P
    call real_times_mat_subr(-1E0_realk, A, P)
  end function

  function integer_times_mat(r, A) result(P)
    integer,              intent(in) :: r
    type(matrix), target, intent(in) :: A
    type(proxy),  target             :: P
    call real_times_mat_subr(r*1E0_realk, A, P)
  end function

  function real_times_mat(r, A) result(P)
    real(realk),          intent(in) :: r
    type(matrix), target, intent(in) :: A
    type(proxy),  target             :: P
    call real_times_mat_subr(r, A, P)
  end function

  function complex_times_mat(z, A) result(P)
    complex(realk),       intent(in) :: z
    type(matrix), target, intent(in) :: A
    type(proxy),  target             :: P
    ! early return if z real
    if (aimag(z)==0) call real_times_mat_subr(real(z), A, P)
    if (aimag(z)==0) return
    ! imag part, then add real part
    call clear_proxy(P)
    call mat_imag_times_mat(aimag(z), A, P%A)
    call mat_axpy(real(z), A, .false., .false., P%A)
    call prepare_result_proxy(P)
  end function

  subroutine real_times_mat_subr(r, A, P)
    real(realk),          intent(in) :: r
    type(matrix), target, intent(in) :: A
    type(proxy),  target             :: P
    call clear_proxy(P)
    P%f = r * merge(0, 1, iszero(A))
    call mat_dup(A, P%A)
    call prepare_result_proxy(P)
  end subroutine

  function plus_prx(P) result(R)
    type(proxy), target, intent(in) :: P
    type(proxy), target             :: R
    call real_times_prx_subr(1E0_realk, P, R)
  end function

  function minus_prx(P) result(R)
    type(proxy), target, intent(in) :: P
    type(proxy), target             :: R
    call real_times_prx_subr(-1E0_realk, P, R)
  end function

  function integer_times_prx(f, P) result(R)
    integer,             intent(in) :: f
    type(proxy), target, intent(in) :: P
    type(proxy), target             :: R
    call real_times_prx_subr(f*1E0_realk, P, R)
  end function

  function real_times_prx(f, P) result(R)
    real(realk),         intent(in) :: f
    type(proxy), target, intent(in) :: P
    type(proxy), target             :: R
    call real_times_prx_subr(f, P, R)
  end function

  function complex_times_prx(z, P) result(R)
    complex(realk),      intent(in) :: z
    type(proxy), target, intent(in) :: P
    type(proxy), target             :: R
    call complex_times_prx_subr(z, P, R)
  end function

  subroutine real_times_prx_subr(f, P, R)
    real(realk),         intent(in)  :: f
    type(proxy), target              :: P
    type(proxy), target, intent(out) :: R
    call rectify_input_proxy(P)
    call clear_proxy(R)
    R%f = f*P%f
    if (R%f==0 .and. P%hb) then
       call mat_gemm(R%f, P%A, P%ta, P%B, P%tb, .true., R%A)
       call mat_remove(P%A)
       call mat_remove(P%B)
    else if (R%f==0 .and. .not.P%hb) then
       call mat_axpy(R%f, P%A, P%ta, .true., R%A)
       call mat_remove(P%A)
    else
       call mat_move(P%A, R%A)
       R%ta = P%ta
       R%hb = P%hb
       if (P%hb) call mat_move(P%B, R%B)
       if (P%hb) R%tb = P%tb
    end if
    call prepare_result_proxy(R)
  end subroutine

  subroutine complex_times_prx_subr(z, P, R)
    complex(realk),      intent(in)  :: z
    type(proxy), target              :: P
    type(proxy), target, intent(out) :: R
    ! early return if z is real
    if (aimag(z)==0) call real_times_prx_subr(real(z), P, R)
    if (aimag(z)==0) return
    call rectify_input_proxy(P)
    call clear_proxy(R)
    ! reduce P of the form P=f*A*B to f*A
    if (P%hb) then
       call mat_gemm(P%f, P%A, P%ta, P%B, P%tb, .true., R%A)
       call mat_move(R%A, P%A)
       call mat_remove(P%B)
       call clear_proxy(P)
    end if
    ! place i*zi*P%f*P%A in R%A, and add zr*P%f*P%A to it
    call mat_imag_times_mat(aimag(z)*P%f, P%A, R%A)
    call mat_axpy(real(z)*P%f, P%A, .false., .false., R%A)
    call mat_remove(P%A)
    R%ta = P%ta
    call prepare_result_proxy(R)
  end subroutine


  !----------------------------------
  !   transpose operator trps(A)
  !----------------------------------

  ! transpose of matrix
  function transpose_mat(A) result(P)
    type(matrix), target, intent(in) :: A
    type(proxy),  target             :: P
    call clear_proxy(P)
    ! if zero, use axpy, otherwise
    if (iszero(A)) then
       call mat_axpy(1E0_realk, A, .true., .true., P%A)
    else
       call mat_dup(A, P%A)
       P%ta = .true.
    end if
  end function

  ! transpose of matrix proxy
  function transpose_prx(P) result(R)
    type(proxy), target :: P
    type(proxy), target :: R
    call rectify_input_proxy(P)
    R%f = P%f
    if (P%hb) then
       P%hb = .true.
       call mat_move(P%A, R%B)
       R%tb = .not.P%ta
       call mat_move(P%B, R%A)
       R%ta = .not.P%tb
    else if (P%f/=0) then
       call mat_move(P%A, R%A)
       R%ta = .not.P%ta
    else !zero
       call mat_axpy(0E0_realk, P%A, .true., .true., R%A)
       call mat_nullify(P%A)
    end if
  end function


  !-----------------------------------------
  !   matrix-matrix multiply operator A*B
  !-----------------------------------------

  function mat_times_mat(A, B) result(P)
    type(matrix), target, intent(in) :: A, B
    type(proxy),  target             :: P
    call clear_proxy(P)
    if (iszero(A) .or. iszero(B)) then
       P%f = 0
       call mat_gemm(P%f, A, .false., B, .false., .true., P%A)
    else
       call mat_dup(A, P%A)
       call mat_dup(B, P%B)
       P%hb = .true.
    end if
    call prepare_result_proxy(P)
  end function

  function mat_times_prx(A, P) result(R)
    type(matrix), target, intent(in) :: A
    type(proxy),  target, intent(in) :: P
    type(proxy),  target             :: R
    call mat_revtimes_prx(.false., A, .false., .false., P, R)
  end function

  function prx_times_mat(P, A) result(R)
    type(proxy),  target, intent(in) :: P
    type(matrix), target, intent(in) :: A
    type(proxy),  target             :: R
    call mat_revtimes_prx(.false., A, .false., .true., P, R)
  end function

  function prx_times_prx(P, Q) result(R)
    type(proxy), target, intent(in) :: P, Q
    type(proxy), target             :: R
    call prx_times_prx_subr(P, Q, R)
  end function


  ! R = A^ta * P or R = P * A^ta (rev=true)
  subroutine mat_revtimes_prx(atmp, A, ta, rev, P, R)
    logical,              intent(in)  :: atmp, ta, rev
    type(matrix), target              :: A
    type(proxy),  target              :: P
    type(proxy),  target, intent(out) :: R
    call clear_proxy(R)
    ! first reduce P from the form f*A^ta*B^tb to 1*A
    call rectify_input_proxy(P)
    if (.not.iszero(A) .and. P%hb) then
       call mat_gemm(P%f, P%A, P%ta, P%B, P%tb, .true., R%A)
       call mat_move(R%A, P%A)
       P%f  = 1
       P%ta = .false.
       call mat_remove(P%B)
       P%hb = .false.
       P%tb = .false.
    end if
    ! deal with zero separately
    R%f  = P%f * merge(0, 1, iszero(A))
    R%hb = (R%f/=0)
    if (R%f==0) then
       if (.not.P%hb .and. .not.rev) then
          call mat_gemm(R%f, A, ta, P%A, P%ta, .true., R%A)
       else if (.not.P%hb .and. rev) then
          call mat_gemm(R%f, P%A, P%ta, A, ta, .true., R%A)
       else if (P%hb .and. .not.rev) then
          call mat_gemm(R%f, A,    ta,     P%A, P%ta, .true., R%B)
          call mat_gemm(R%f, R%B, .false., P%B, P%tb, .true., R%A)
          call mat_nullify(R%B)
       else if (P%hb .and. rev) then
          call mat_gemm(R%f, P%B, P%tb, A,    ta,     .true., R%B)
          call mat_gemm(R%f, P%A, P%ta, R%B, .false., .true., R%A)
          call mat_nullify(R%B)
       end if
       if (atmp) call mat_remove(A)
       call mat_remove(P%A)
       if (P%hb) call mat_remove(P%B)
    ! nonzero R = f * A^ta * P%A^P%ta
    else if (.not.rev) then
       if (.not.atmp) call mat_dup(A, R%A)
       if (     atmp) call mat_move(A, R%A)
       R%ta = ta
       call mat_move(P%A, R%B)
       R%tb = P%ta
    ! nonzero R = f * P%A^P%ta * A^ta
    else if (rev) then
       call mat_move(P%A, R%A)
       R%ta = P%ta
       if (.not.atmp) call mat_dup(A, R%B)
       if (     atmp) call mat_move(A, R%B)
       R%tb = ta
    end if
    call prepare_result_proxy(R)
  end subroutine


  subroutine prx_times_prx_subr(P, Q, R)
    type(proxy), target              :: P, Q
    type(proxy), target, intent(out) :: R
    ! first reduce Q from the form f*A^ta*B^tb to 1*A
    call rectify_input_proxy(Q)
    if (P%f/=0 .and. Q%hb) then
       call mat_gemm(P%f * Q%f, Q%A, Q%ta, Q%B, Q%tb, .true., R%A)
       P%f = 1
       call mat_move(R%A, Q%A)
       Q%f = 1
       Q%ta = .false.
       call mat_remove(Q%B)
       Q%hb = .false.
       Q%tb = .false.
    end if
    ! now produce R as P*Q%A
    if (Q%f/=0) P%f = P%f*Q%f
    call mat_revtimes_prx(.true., Q%A, Q%ta, .true., P, R)
  end subroutine

  
  !--------------------------------------------
  !   dot(A,B), tr(A,B), tr(A) and norm(A)
  !--------------------------------------------

  function mat_dot_mat(A, B)
    type(matrix), target, intent(in) :: A, B
    complex(realk) mat_dot_mat
    mat_dot_mat = 0
    if (.not.iszero(A) .and. .not.iszero(B)) &
       mat_dot_mat = mat_dot(A, B, .false.)
  end function

  function mat_prodtr_mat(A, B)
    type(matrix), target :: A, B
    complex(realk) mat_prodtr_mat
    mat_prodtr_mat = 0
    if (.not.iszero(A) .and. .not.iszero(B)) &
       mat_prodtr_mat = mat_dot(A, B, .true.)
  end function

  function mat_dot_prx(A, P)
    type(matrix), target, intent(in) :: A
    type(proxy),  target, intent(in) :: P
    complex(realk) mat_dot_prx
    call mat_dot_prx_subr(A, P, .false., mat_dot_prx)
  end function

  function mat_prodtr_prx(A, P)
    type(matrix), target, intent(in) :: A
    type(proxy),  target, intent(in) :: P
    complex(realk) mat_prodtr_prx
    call mat_dot_prx_subr(A, P, .true., mat_prodtr_prx)
  end function

  function prx_dot_mat(P, A)
    type(proxy),  target, intent(in) :: P
    type(matrix), target, intent(in) :: A
    complex(realk) prx_dot_mat
    call mat_dot_prx_subr(A, P, .false., prx_dot_mat)
  end function

  function prx_prodtr_mat(P, A)
    type(proxy),  target, intent(in) :: P
    type(matrix), target, intent(in) :: A
    complex(realk) prx_prodtr_mat
    call mat_dot_prx_subr(A, P, .true., prx_prodtr_mat)
  end function

  function prx_dot_prx(P, Q)
    type(proxy),  target, intent(in) :: P, Q
    complex(realk) prx_dot_prx
    call prx_dot_prx_subr(Q, P, .false., prx_dot_prx)
  end function

  function prx_prodtr_prx(P, Q)
    type(proxy),  target, intent(in) :: P, Q
    complex(realk) prx_prodtr_prx
    call prx_dot_prx_subr(P, Q, .true., prx_prodtr_prx)
  end function

  subroutine mat_dot_prx_subr(M, P, t, r)
    type(matrix),   target      :: M
    type(proxy),    target      :: P
    logical,        intent(in)  :: t !t=F: dot(M,P), t=T: tr(M*P)
    complex(realk), intent(out) :: r
    type(matrix) AB
    call rectify_input_proxy(P)
    if (iszero(M)) P%f = 0*P%f
    if (P%hb) then
       call mat_gemm(P%f, P%A, P%ta, P%B, P%tb, .true., AB)
       if (P%f/=0) P%f = 1
       call mat_move(AB, P%A)
       P%ta = .false.
       call mat_remove(P%B)
    end if
    r = P%f * mat_dot(M, P%A, t.neqv.P%ta)
    call mat_remove(P%A)
  end subroutine

  subroutine prx_dot_prx_subr(P, Q, t, r)
    type(proxy),    target      :: P, Q
    logical,        intent(in)  :: t !t=F: dot(P,Q), t=T: tr(P*Q)
    complex(realk), intent(out) :: r
    type(matrix) PAB
    call rectify_input_proxy(P)
    if (Q%f==0) P%f = 0*P%f
    if (P%hb) then
       call mat_gemm(P%f, P%A, P%ta, P%B, P%tb, .true., PAB)
       if (P%f/=0) P%f = 1
       call mat_move(PAB, P%A)
       P%ta = .false.
       call mat_remove(P%B)
    end if
    call mat_dot_prx_subr(P%A, P, t.neqv.P%ta, r)
    r = r * P%f
  end subroutine

  function prx_trace(P)
    type(proxy), target :: P
    complex(realk) prx_trace
    call rectify_input_proxy(P)
    if (.not.P%hb) prx_trace = P%f * mat_trace(P%A)
    if (     P%hb) prx_trace = P%f * mat_dot(P%A, P%B, P%ta.neqv.P%tb)
    if (P%hb) call mat_remove(P%B)
    call mat_remove(P%A)
  end function

  function mat_norm(A)
    type(matrix), target :: A
    real(realk) mat_norm
    mat_norm = 0
    if (.not.iszero(A)) mat_norm = sqrt(real(mat_dot(A, A, .false.)))
  end function

  function prx_norm(P)
    type(proxy), target :: P
    real(realk) prx_norm
    type(matrix) AB
    if (P%hb) then
       call mat_gemm(P%f, P%A, P%ta, P%B, P%tb, .true., AB)
       if (P%f/=0) P%f = 1
       call mat_move(AB, P%A)
       P%ta = .false.
       call mat_remove(P%B)
    end if
    prx_norm = abs(P%f) * mat_norm(P%A)
    call mat_remove(P%A)
  end function


  !---------------------------------------------------------------------------
  !   operators for removing: A=0, A(:)=0, A(:,:)=0, A(:,:,:)=0, A(:,:,:,:)=0
  !---------------------------------------------------------------------------

  ! Short-hand A=0 removes matrix A if it is defined, otherwise nullifies
  subroutine mat_eq_zero(A, must_be_zero)
    type(matrix), target, intent(inout) :: A
    integer,              intent(in)    :: must_be_zero
    if (must_be_zero/=0) &
       call quit('matrix A = z error: z must be 0')
    if (isdef(A)) then
       call mat_remove(A)
    else
       call mat_nullify(A)
    end if
  end subroutine

  subroutine mat_eq_zero_1D(A, z)
    type(matrix), target, intent(inout) :: A(:)
    integer,              intent(in)    :: z
    integer i
    do i=1, size(A)
       A(i) = z
    end do
  end subroutine

  subroutine mat_eq_zero_2D(A, z)
    type(matrix), target, intent(inout) :: A(:,:)
    integer,              intent(in)    :: z
    integer j
    do j=1, size(A,2)
       A(:,j) = z
    end do
  end subroutine

  subroutine mat_eq_zero_3D(A, z)
    type(matrix), target, intent(inout) :: A(:,:,:)
    integer,              intent(in)    :: z
    integer k
    do k=1, size(A,3)
       A(:,:,k) = z
    end do
  end subroutine

  subroutine mat_eq_zero_4D(A, z)
    type(matrix), target, intent(inout) :: A(:,:,:,:)
    integer,              intent(in)    :: z
    integer l
    do l=1, size(A,4)
       A(:,:,:,l) = z
    end do
  end subroutine

  ! clear fields in proxy. Note that matrix fields A,B are not nullified
  subroutine clear_proxy(P)
    type(proxy), intent(out) :: P
    P%f  = 1
    P%ta = .false.
    P%hb = .false.
    P%tb = .false.
    P%is_delayed = .false.
  end subroutine

  ! Some compilers (eg gfortran) move type(proxy) back and forth during
  ! formula evaluation. This, however dissociates P%C%self_pointer => P%C
  ! (and the same for proxy fields P%A and P%C), making the matrices seem
  ! like aliases, and not temporaries, which results in memory leaks.
  ! To get around this, we invalidate (hide) the matrices before returning
  ! a proxy, and restore them when a proxy later is received as input
  subroutine prepare_result_proxy(P)
    type(proxy), target, intent(inout) :: P
    if (P%f/=0) call mat_hide_temp(P%A)
    if (P%hb)   call mat_hide_temp(P%B)
  end subroutine

  subroutine rectify_input_proxy(P)
    type(proxy), target, intent(inout) :: P
    type(matrix) A
    if (P%f/=0) call mat_unhide_temp(P%A)
    if (P%hb)   call mat_unhide_temp(P%B)
    ! if this proxy was delayed, the matrix 'delayed_matrix' should be added to it
    if (P%is_delayed) then
       if (.not.have_delayed_matrix) call quit( &
             'P%is_delayed .and. .not.have_delayed_matrix should be impossible')
       call mat_dup(delayed_matrix, A)
       nullify(delayed_matrix)
       have_delayed_matrix = .false.
       P%is_delayed = .false.
       if (.not.P%hb) call mat_axpy(P%f, P%A, P%ta, .false., A)
       if (     P%hb) call mat_gemm(P%f, P%A, P%ta, P%B, P%tb, .false., A)
       call mat_move(A, P%A)
       if (P%hb) call mat_remove(P%B)
       P%f = 1; P%ta = .false.; P%hb = .false.
    end if
  end subroutine

end module
