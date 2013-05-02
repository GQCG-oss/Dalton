MODULE mm_W_worker

   USE mm_global_paras_mod
 
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_get_ltsqr_W_matrix,    &
             mm_contract_Wq

CONTAINS

!-------------------------------------------------------------------------------
! routine to build the lower-triangular elements of a square W_matrix.
! only the non-zero elements of the full square matrix are written.
! LMAX is the order of translation: 
!  W_lm,jk is built for all l and j <= LMAX
! r_ab is the translation vector (final-initial).
! (note to translate a potential -r_ab should be passed, and W' used) 

   SUBROUTINE mm_get_ltsqr_W_matrix(LMAX,JMAX,r_ab,W_matrix)

      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: LMAX, JMAX
      REAL(REALK),   INTENT(IN)  :: r_ab(3)
      REAL(REALK),   INTENT(OUT) :: W_matrix(:,:)
      INTEGER :: n,m

      n = MAX(LMAX,JMAX)
      m = MIN(LMAX,JMAX)
      CALL mm_build_ltsqr_W_matrix

   CONTAINS

      SUBROUTINE mm_build_ltsqr_W_matrix
         IMPLICIT NONE
         REAL(REALK) :: R_sh(-n:n,0:n)

         CALL mm_generate_R(n,-r_ab,R_sh)
          ! note minus sign above - see TUH eqn. (9.13.58)
         CALL mm_generate_W(n,m,R_sh,W_matrix)

      END SUBROUTINE mm_build_ltsqr_W_matrix

   END SUBROUTINE mm_get_ltsqr_W_matrix

!-------------------------------------------------------------------------------

   SUBROUTINE mm_generate_R(LMAX,point,R_sh) 
 
   ! Subroutine to generate scaled regular solid harmonics
   ! See page 416, "Molecular Electronic Structure Theory"
   !-------------------------------------------------------
   ! R_sh(M,L)  == COS component of R_sh(M,L) 
   ! R_sh(-M,L) == SIN component of R_sh(-M,L) 
   ! where L and M are positive
   ! Note that R_sh(0,L) == COS component of R_sh(0,L) !
   !-------------------------------------------------------

      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: LMAX
      REAL(REALK),   INTENT(IN)  :: point(3)
      REAL(REALK),   INTENT(OUT) :: R_sh(-LMAX:LMAX,0:LMAX) 

      REAL(REALK)   :: x,y,z, r_2,r_minus2 
      REAL(REALK)   :: A(6) 
      INTEGER :: signl, l,m
    
      A(:) = zero
      x = point(1)
      y = point(2)
      z = point(3)
      r_2 = DOT_PRODUCT(point,point)
      r_minus2 = one/(r_2)

      IF (r_2 < ZERO_VECT_TOL*ZERO_VECT_TOL) THEN
         PRINT '(3E25.15)', point 
         CALL LSQUIT('ERROR: Why do we try to do zero W-vector.',-1)
      END IF

      IF (LMAX==0) THEN
         R_sh(0,0) = one
         RETURN
      END IF

      ! Start with initialization of lowest order terms

      R_sh( 0,0) = one       ! cannot assign "sin" terms for m==0
      R_sh(-1,1) = -half*y
      R_sh( 0,1) = z         ! cannot assign "sin" terms for m==0
      R_sh( 1,1) = -half*x

      ! Now iterate higher order terms

      signl = -1
      lloop: DO l = 2, LMAX
         signl = -signl
         A(1) = one/(2*l)
         A(2) = A(1)*x
         A(3) = A(1)*y*signl

         ! cos terms:
         R_sh( l,l) = A(3)*R_sh(1-l,l-1) - A(2)*R_sh(l-1,l-1)
         ! sin terms:
         R_sh(-l,l) = A(3)*R_sh(l-1,l-1) + A(2)*R_sh(1-l,l-1)

         A(5) = (2*l-1)*(z)*(r_minus2)
         DO m = 0, (l-2)
            A(6) = (r_2)/(l*l-m*m)
            ! cos terms:
            R_sh( m,l) = A(6)*( A(5)*R_sh(m,l-1) - R_sh(m,l-2) ) 
            ! sin terms:
            R_sh(-m,l) = A(6)*( A(5)*R_sh(-m,l-1) - R_sh(-m,l-2) ) 
         END DO

         ! Now do (l,l-1) terms to avoid calling elements that do not exist 
         m = l-1
         ! cos terms:
         R_sh( m,l) = z*R_sh(m,l-1)
         ! sin terms:
         R_sh(-m,l) = z*R_sh(-m,l-1)

      END DO lloop

   END SUBROUTINE mm_generate_R

!-------------------------------------------------------------------------------

   SUBROUTINE mm_generate_W(LMAX,JMAX,R_sh,W_matrix) 
 
   ! Subroutine to generate Real Translation Matrix for 
   ! given point from the Regular Solid Harmonics R_sh(m,l) 
   ! cf Helgaker et. al. pp 413
   ! Note that W-matrix is lower triangular wrt W(lm,jk) and that
   ! only NON-ZERO elements are added here.
   ! We must therefore zero W_matrix initially, or use a special contractor
   !------------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: LMAX, JMAX 
      REAL(REALK),   INTENT(INOUT) :: R_sh(-LMAX:LMAX,0:LMAX) 
      REAL(REALK),   INTENT(OUT) :: W_matrix(:,:) 

      INTEGER :: l,m,j,k, a,b, p,q, pp,qq, lm_min,lm_max
      INTEGER :: signk, phase

      ! l==j case first:
      ! diagonal = one
      ! all other off-diagonal elements identically zero
      DO j = 0, JMAX
         qq = j*(j+1) +1
         DO k = -j, j
            q = qq+k
            W_matrix(q,q) = one
         END DO
      END DO

      IF (LMAX==0) RETURN
      IF (JMAX > LMAX) CALL LSQUIT('error in mm_generate_W',-1)

      ! remaining l:
      jloop: DO j = 0, JMAX
         qq = j*(j+1) +1

         ! loop over positive k
         signk = -1
         kloop1: DO k = 0, j
            signk = -signk
            q = qq+k

            lloop1: DO l = (j+1), LMAX
               a = (l-j)
               pp = l*(l+1) +1

               ! cos-cos terms
               !--------------

               ! (m-k) terms
               phase = 1 
               DO m = (k-1), MAX((k-a),0), -1     ! reverse order 
                  ! (m-k)<0
                  phase = -phase 
                  b = -(m-k)
                  W_matrix(pp+m,q) = phase*R_sh(b,a)
               END DO
               DO m = k, (a+k)      ! (a+k)<l
                  ! (m-k)>=0
                  b = (m-k)
                  W_matrix(pp+m,q) = R_sh(b,a)
               END DO
               ! (m+k) terms
               DO m = 0, (a-k)   ! always overlaps (m-k) terms
                  ! (m+k)>=0
                  b = (m+k)
                  W_matrix(pp+m,q) = W_matrix(pp+m,q) + signk*R_sh(b,a)
               END DO

               ! sin-cos terms
               !--------------

               ! m==-k first
               m = MIN(-k,-1) 
               W_matrix(pp+m,q) = zero
               ! (m+k) terms
               DO m = -(k+a), -(k+1)
                  ! (m+k)<0
                  b = (m+k)
                  W_matrix(pp+m,q) = signk*R_sh(b,a)
               END DO
               phase = -signk
               DO m = -k+1, MIN(-1,(a-k)) 
                  ! (m+k)>0
                  phase = -phase 
                  b = -(m+k)
                  W_matrix(pp+m,q) = phase*R_sh(b,a)
               END DO
               ! (m-k) terms
               DO m = (k-a), -1   ! always overlap (m+k) terms
                  ! (m-k)<0
                  b = (m-k)
                  W_matrix(pp+m,q) = W_matrix(pp+m,q) + R_sh(b,a)
               END DO

            END DO lloop1
         END DO kloop1

         ! (half)**(delta(0,k)) correction:
         lm_min = j*(j+2) +2
         lm_max = LMAX*(LMAX+2) +1
         W_matrix(lm_min:lm_max,qq) = half*W_matrix(lm_min:lm_max,qq)   !k==0

         ! loop over negative k
         signk = 1
         kloop2: DO k = -1, -j, -1     ! reverse order
            signk = -signk
            q = qq+k

            lloop2: DO l = (j+1), LMAX
               a = (l-j)
               pp = l*(l+1) +1

               ! cos-sin terms
               !--------------
               ! (m+k) terms
               DO m = MAX(0,-(a+k)), (-k-1) 
                  ! (m+k)<0
                  b = (m+k)
                  W_matrix(pp+m,q) = signk*R_sh(b,a)
               END DO
               m = -k 
               W_matrix(pp+m,q) = zero
               phase = -signk
               DO m = -k+1, (a-k)  !<=l 
                  ! (m+k)>=0
                  phase = -phase
                  b = -(m+k)
                  W_matrix(pp+m,q) = phase*R_sh(b,a)
               END DO
               ! (m-k) terms
               phase = signk
               DO m = 0, (a+k)    ! always overlap (m+k) terms
                  ! (m-k)>0
                  phase = -phase
                  b = -(m-k)
                  W_matrix(pp+m,q) = W_matrix(pp+m,q) - phase*R_sh(b,a)
               END DO

               ! sin-sin terms
               !--------------

               ! (m-k) terms
               phase = 1
               DO m = (k-1), (k-a), -1    ! reverse order
                  ! (m-k)<0
                  phase = -phase
                  b = -(m-k)
                  W_matrix(pp+m,q) = phase*R_sh(b,a)
               END DO
               DO m = k, MIN(-1,(a+k)) 
                  ! (m-k)>=0
                  b = (m-k)
                  W_matrix(pp+m,q) = R_sh(b,a)
               END DO
               ! (m+k) terms
               phase = 1 
               DO m = -1, -(a+k), -1     ! reverse order
                  ! (m+k)<0              ! always overlap (m-k) terms
                  phase = -phase
                  b = -(m+k)
                  W_matrix(pp+m,q) = W_matrix(pp+m,q) - phase*R_sh(b,a)
               END DO

            END DO lloop2
         END DO kloop2
      END DO jloop

   END SUBROUTINE mm_generate_W

!-------------------------------------------------------------------------------

   SUBROUTINE mm_contract_Wq(N_or_T,W_mat,WLDA,vec_in,n,vec_out,m)

      IMPLICIT NONE
      CHARACTER(1),  INTENT(IN)    :: N_or_T
      INTEGER, INTENT(IN)    :: n,m,WLDA
      REAL(REALK),   INTENT(IN)    :: W_mat(WLDA,WLDA)
      REAL(REALK),   INTENT(IN)    :: vec_in(n)
      REAL(REALK),   INTENT(INOUT) :: vec_out(m)

      INTEGER :: q

      IF (N_or_T == 'N') THEN
         DO q = 1, n
            vec_out(q:m) = vec_out(q:m) + W_mat(q:m,q)*vec_in(q)
         END DO
      ELSE
         ! Perform transpose contraction
         DO q = 1, m 
            vec_out(q) = vec_out(q) + DOT_PRODUCT(W_mat(q:n,q),vec_in(q:n))
         END DO
      END IF

   END SUBROUTINE mm_contract_Wq

!-------------------------------------------------------------------------------

END MODULE mm_W_worker


!************************************************************
!* These subroutines are placed outside the module in order *
!* to be accessible from Fortran 77 code.                   *
!************************************************************


! This is a wrapper for the translation (W) matrix generating routine.
! It's only function is to enable calls to mm_get_ltsqr_W_matrix from
! the PBC code.
!
! Written by Erik Tellgren, October 2005
!
SUBROUTINE pbc_get_ltsqr_W_matrix(LMAX,r_ab,W_matrix)
  use mm_global_paras_mod
  use mm_W_worker, only: mm_get_ltsqr_W_matrix
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: LMAX
  REAL(REALK),   INTENT(IN)  :: r_ab(3)
  REAL(REALK),   INTENT(OUT) :: W_matrix((1+LMAX)**2,(1+LMAX)**2)

  integer :: jmax

  jmax = LMAX

  call mm_get_ltsqr_W_matrix(LMAX,jmax,r_ab,W_matrix)

end SUBROUTINE pbc_get_ltsqr_W_matrix



subroutine pbc_restore_W_matrix(lmax,W_NF)
  use mm_global_paras_mod
  implicit none

  integer, intent(in) :: lmax
  real(realk), intent(inout) :: W_NF((1+lmax)**2,(1+lmax)**2)

  integer :: l,j,row_min,row_max,col_min,col_max

  write(LUPRI,*)
  write(LUPRI,*) 'Transforming a T tensor from triangular to square form.'
  write(LUPRI,*)

  do l = 0,lmax
     ! determine index range
     row_min = 1 + l**2
     row_max = (1 + l)**2
     do j = l,lmax
        if (mod(j+l,2) .eq. 1) then
           ! determine index range
           col_min = 1 + j**2
           col_max = (1 + j)**2
           ! flip sign
           W_NF(row_min:row_max,col_min:col_max) = -W_NF(row_min:row_max,col_min:col_max)
        end if
     end do
  end do
end subroutine pbc_restore_W_matrix
