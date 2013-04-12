MODULE mm_T_worker

   USE mm_global_paras_mod

   IMPLICIT NONE
   PRIVATE
       REAL(REALK), SAVE :: THRESH_T = 1.d-16
   ! Public procedures
   PUBLIC :: mm_get_SPLTSQ_T_matrix,   &
             mm_get_FLTSQ_T_matrix,    &
             mm_generate_I,            &    ! added by ErikT
             mm_contract_Tq,           &
             mm_postfac_Vff,           &
             mm_scale_vec,             &
             mm_set_T_worker_threshold

CONTAINS

   SUBROUTINE mm_set_T_worker_threshold(threshold)
   implicit none
   REAL(realk) :: threshold
     THRESH_T = threshold
   END SUBROUTINE mm_set_T_worker_threshold

!------------------------------------------------------------------------------
! Get sparse (l+j <= LMAX), lower-triangular, square-based T-matrix
!FIXME: would be nice to allow possibility of building rectangular-based too

   SUBROUTINE mm_get_SPLTSQ_T_matrix(LMAX,r_ab,T_matrix)

      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: LMAX
      REAL(REALK),   INTENT(IN)  :: r_ab(3)
      REAL(REALK),   INTENT(INOUT) :: T_matrix(:,:)

      REAL(REALK) :: I_sh((1+LMAX)**2)
      INTEGER, EXTERNAL :: get_JMAX

      CALL mm_generate_I(LMAX,r_ab,I_sh)
      CALL mm_generate_T(LMAX,mm_SPLTSQ_JMAX,I_sh,T_matrix)

   END SUBROUTINE mm_get_SPLTSQ_T_matrix

   FUNCTION mm_SPLTSQ_JMAX(L,LMAX)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: L, LMAX
      INTEGER :: mm_SPLTSQ_JMAX
      mm_SPLTSQ_JMAX = LMAX-L
   END FUNCTION mm_SPLTSQ_JMAX

!------------------------------------------------------------------------------
! Get full (all elements built), lower-triangular, square T-matrix

   SUBROUTINE mm_get_FLTSQ_T_matrix(LMAX,r_ab,T_matrix)

      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: LMAX
      REAL(REALK),   INTENT(IN)  :: r_ab(3)
      REAL(REALK),   INTENT(INOUT) :: T_matrix(:,:)  ! changed "out" to "inout" / ErikT

      REAL(REALK) :: I_sh((1+2*LMAX)**2)
      INTEGER, EXTERNAL :: get_JMAX

      CALL mm_generate_I(2*LMAX,r_ab,I_sh)
      CALL mm_generate_T(LMAX,mm_FLTSQ_JMAX,I_sh,T_matrix)

   END SUBROUTINE mm_get_FLTSQ_T_matrix

   FUNCTION mm_FLTSQ_JMAX(l,LMAX)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: l, LMAX
      INTEGER :: mm_FLTSQ_JMAX
      mm_FLTSQ_JMAX = LMAX 
   END FUNCTION mm_FLTSQ_JMAX

!------------------------------------------------------------------------------

   SUBROUTINE mm_generate_I(LMAX,vector,I_sh) 
 
      ! Subroutine to generate Scaled Irregular solid harmonics
      ! See page 416, "Molecular Electronic Structure Theory".
      !--------------------------------------------------------
      ! I_sh(M,L)  == COS component of I_sh(M,L)
      ! I_sh(-M,L) == SIN component of I_sh(-M,L)
      ! where L and M are positive
      ! Note that I_sh(0,L) == COS component of I_sh(0,L) !
      ! In fact we combine (l,m) into a single index and
      ! we assume the order 0, -1,1,+1, -2,-1,0,+1,+2 ...
      !--------------------------------------------------------

      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: LMAX
      REAL(REALK),   INTENT(IN)  :: vector(3)
      REAL(REALK),   INTENT(INOUT) :: I_sh((LMAX+1)**2)  ! changed to "inout" / ErikT

      REAL(REALK)   :: tmp1, tmp2, tmp3
      REAL(REALK)   :: x,y,z, r_2, r_minus2
      INTEGER :: L,m, i,j,k,p,q,u, sign_L 

      x = vector(1)
      y = vector(2)
      z = vector(3)
      r_2 = x**2 + y**2 + z**2
      r_minus2 = one/(r_2) 

!FIXME;
      IF (r_2 < ZERO_VECT_TOL) THEN
         PRINT '(3E25.15)', vector 
         CALL LSQUIT('Why do we try to do a zero T_vector?',-1)
      END IF

      IF (LMAX==0) THEN
         I_sh(1) = SQRT(r_minus2)
         RETURN
      END IF
         
      ! Start with initialization of lowest order terms
        I_sh(1) =  SQRT(r_minus2)
        I_sh(2) = -(y)*(r_minus2)*I_sh(1)
        I_sh(3) =  (z)*(r_minus2)*I_sh(1)
        I_sh(4) = -(x)*(r_minus2)*I_sh(1)

      ! Now iterate higher order terms
      sign_L = -1
      DO L = 2, LMAX
         sign_L = -sign_L 
         i = (L+1)*(L+1) 
         j = L*L+1 
         p = j-2*L+1
         q = j-1
         tmp1 = (2*L-1)*r_minus2
         tmp2 = tmp1*y*sign_L
         tmp3 = tmp1*x
         !
         I_sh(i) = tmp2*I_sh(p) - tmp3*I_sh(q)
         I_sh(j) = tmp2*I_sh(q) + tmp3*I_sh(p)
         !
         p = p+L-1 
         q = L*(L-3)+3
         u = q+L-2
         k = i-L
         tmp2 = tmp1*z
         DO m = 0, (L-2)
            tmp3 = (u-m*m)*(r_minus2)
            !
            I_sh(k+m) = tmp2*I_sh(p+m) - tmp3*I_sh(q+m)
            I_sh(k-m) = tmp2*I_sh(p-m) - tmp3*I_sh(q-m)  
            ! I_sh(L,0)=0
            !
         END DO
         ! Now do (L, L-1) terms to avoid calling elements that do not exist 
         m = L-1
         I_sh(k+m) = tmp2*I_sh(p+m)
         I_sh(k-m) = tmp2*I_sh(p-m)
      END DO

   END SUBROUTINE mm_generate_I

!------------------------------------------------------------------------------

   SUBROUTINE mm_generate_T(LMAX,JMAX,I_sh,T_matrix)
 
   ! Subroutine to generate Real Interaction Matrix for 
   ! given point from Irregular Solid Harmonics (I_sh) 
   ! cf Helgaker et. al. pp 415
   !-----------------------------------------------------------------------
   !  T_cos_cos(lm,jk) = 2*{I_cos(l+j,m+k) + (-1)^k*I_cos(l+j,m-k)}  etc.
   ! 
   ! Note we only write the non-zero elements here, and the precise form
   ! depends on the function JMAX passed in
   !-----------------------------------------------------------------------
   !FIXME: can optimise for j=l by noting that if we build just lower
   ! triangular, we need only half the elements there

      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: LMAX
      REAL(REALK),   INTENT(IN)  :: I_sh(:)
      REAL(REALK),   INTENT(INOUT) :: T_matrix(:,:)
      INTEGER, EXTERNAL    :: JMAX

      INTEGER :: L,m,J,k, p,q, pp,qq, u
      INTEGER :: sign_k, sign_m, km

      IF (LMAX==0) THEN
         T_matrix(1,1) = two*I_sh(1)     ! NB scaling to make T symmetric
         RETURN
      END IF

      L_loop: DO L = 0, LMAX 
         qq = L*(L+1) +1
         J_loop: DO J = L, JMAX(L,LMAX)
            pp = J*(J+1) +1
            u = 1+(L+J)*(L+J+1)

            sign_m = -1
            positive_m: DO m = 0, L 
               q = qq+m
               sign_m = -sign_m

               ! cos-cos terms
               !---------------
               sign_k = -1
               DO k = 0, MIN(m,J)  
                  sign_k = -sign_k 
                  ! (m+k)>0 and (m-k)>0
                  T_matrix(pp+k,q) = I_sh(u+(m+k)) + sign_k*I_sh(u+(m-k))
               END DO

               DO k = (m+1), J
                  ! (m+k)>0 and (m-k)<0
                  T_matrix(pp+k,q) = I_sh(u+(m+k)) + sign_m*I_sh(u+(-(m-k)))
               END DO

               ! cos-sin terms
               !---------------

               ! only need to consider (m+k)=0 separately as cannot
               ! also have (m-k)=0 unless k=0

               DO k = (-J), MIN((-m-1),-1)
                  ! (m+k)<0 and (m-k)>0
                  T_matrix(pp+k,q) = I_sh(u+(m+k)) + sign_m*I_sh(u+(-(m-k)))
               END DO
 
               DO k = MAX(-J,-m), MIN(-m,-1)
                  ! (m+k)=0 and (m-k)>0
                  T_matrix(pp+k,q) = sign_m*I_sh(u+(-(m-k)))
               END DO

               sign_k = 1
               DO k = -1, MAX((-m+1),(-J)), -1 
                  sign_k = -sign_k 
                  ! (m+k)>0 and (m-k)>0
                  T_matrix(pp+k,q) = -sign_m*sign_k*I_sh(u+(-(m+k)))         &
                                     +sign_m*I_sh(u+(-(m-k)))
               END DO
 
            END DO  positive_m
 
            sign_m = 1
            negative_m: DO m = -1, -L, -1  
               q = qq+m
               sign_m = -sign_m
 
               ! sin-cos terms
               !---------------

               k = MIN(J,(-m))
               sign_k = 1
               IF (BTEST(k,0)) sign_k = -1
               ! (m+k)=0 and (m-k)<0
               T_matrix(pp+k,q) = sign_k*I_sh(u+(m-k))

               sign_k = -1 
               DO k = 0, MIN((-m-1),J)
                  sign_k = -sign_k 
                  ! (m+k)<0 and (m-k)<0
                 T_matrix(pp+k,q) = I_sh(u+(m+k)) + sign_k*I_sh(u+(m-k))
               END DO

               sign_k = -1
               IF (BTEST(J,0)) sign_k = 1
               DO k = J, (-m+1), -1 
                  sign_k = -sign_k 
                  ! (m+k)>0 and (m-k)<0
                  T_matrix(pp+k,q) = -sign_m*sign_k*I_sh(u+(-(m+k)))         &
                                     +sign_k*I_sh(u+(m-k))
               END DO

               ! sin-sin terms
               !---------------

               sign_k = -1
               IF (BTEST(J,0)) sign_k = 1
               DO k = (-J), m 
                  sign_k = -sign_k 
                  ! (m+k)<0 and (m-k)>0
                  T_matrix(pp+k,q) = -sign_m*sign_k*I_sh(u+(-(m+k)))         &
                                     +sign_k*I_sh(u+(m-k))
               END DO

               sign_k = 1 
               DO k = -1, MAX((-J),(m+1)), -1 
                  sign_k = -sign_k 
                  ! (m+k)<0 and (m-k)<0
                  T_matrix(pp+k,q) = -sign_m*sign_k*I_sh(u+(-(m+k)))         &
                                     +sign_m*I_sh(u+(-(m-k)))
               END DO

            END DO  negative_m
         END DO     J_loop
      END DO        L_loop

   END SUBROUTINE mm_generate_T

!------------------------------------------------------------------------------
! FIXME: could maybe also add other factors here to avoid storing 2 RAW_MM's
!   
!   SUBROUTINE mm_generate_U(LMAX,q_jk,T_matrix,U_matrix)
!
!      IMPLICIT NONE
!
!      INTEGER, INTENT(IN)  :: LMAX
!      REAL(REALK),   INTENT(IN)  :: q_jk((LMAX+1)**2)
!      REAL(REALK),   INTENT(IN)  :: T_matrix((LMAX+1)**2,(LMAX+1)**2)
!      REAL(REALK),   INTENT(OUT) :: U_matrix((LMAX+1)**2,LMAX+1)
!!      REAL(REALK),   INTENT(IN)  :: q_jk(:)
!!      REAL(REALK),   INTENT(IN)  :: T_matrix(:,:)
!!      REAL(REALK),   INTENT(OUT) :: U_matrix(:,:)
!
!      INTEGER :: j,k,p,q,u,v
!      INTEGER :: hi
!
!      jloop: DO j = 0, LMAX
!         hi = (LMAX-j+1)**2
!         U_matrix(1:hi,j+1) = zero
!         p = j*j +1  
!         q = (j+1)**2
!         uloop: DO u = p,q
!            IF (u>hi) EXIT uloop 
!            U_matrix(u:hi,j+1) = U_matrix(u:hi,j+1)          &
!                               + T_matrix(u:hi,u)*q_jk(u)
!         END DO uloop
!      END DO jloop
!
!      DO j = 1, LMAX
!         hi = MIN( (LMAX/2 +1)**2, (LMAX-j+1)**2 )
!         p = j*j +1
!         q = (j+1)**2
!         DO u = 1, hi 
!            IF (p<=u) p=u+1 
!            U_matrix(u,j+1) = U_matrix(u,j+1)                &
!                            + DOT_PRODUCT(T_matrix(p:q,u),q_jk(p:q))
!         END DO
!      END DO
!
!   END SUBROUTINE mm_generate_U
!
!!------------------------------------------------------------------------------
!! FIXME: could maybe also add other factors here to avoid storing 2 RAW_MM's
!   
!   SUBROUTINE mm_generate_V(LMAX,scl_fact,U_matrix,V_lm)
!
!      IMPLICIT NONE
!
!      INTEGER, INTENT(IN)  :: LMAX
!      REAL(REALK),   INTENT(IN)  :: scl_fact
!      REAL(REALK),   INTENT(IN)  :: U_matrix((LMAX+1)**2,LMAX+1) 
!      REAL(REALK),   INTENT(OUT) :: V_lm((LMAX+1)**2)
!!      REAL(REALK),   INTENT(IN)  :: U_matrix(:,:) 
!!      REAL(REALK),   INTENT(OUT) :: V_lm(:)
!
!      INTEGER :: l,m,j,u,v
!      INTEGER :: sign, hi
!      REAL(REALK)   :: prefactor, prefactor1
!
!      sign = 1
!      IF (scl_fact < zero) sign = -1
!      prefactor1 = one/(scl_fact)
!
!      ! j==0 case:
!      V_lm(:) = sign*U_matrix(:,1)
!      prefactor = sign 
!      ! remaining j:
!      DO j = 1, LMAX
!         hi = (LMAX-j+1)**2
!         prefactor = prefactor*prefactor1
!         V_lm(1:hi) = V_lm(1:hi) + prefactor*U_matrix(1:hi,j+1)
!      END DO
!
!      prefactor = one 
!      DO l = 0, LMAX
!         u = l*(l+1) +1
!         prefactor = prefactor*prefactor1
!         DO m = -l, l
!            v = u+m
!            V_lm(v) = prefactor*V_lm(v)
!         END DO
!         V_lm(u) = half*V_lm(u)         ! m==0
!      END DO
!
!   END SUBROUTINE mm_generate_V
!
!------------------------------------------------------------------------------

   FUNCTION mm_contract_Tq(LMAX,vect,T_matrix)

     IMPLICIT NONE
     INTEGER, INTENT(IN) :: LMAX
     REAL(REALK),   INTENT(IN) :: vect((LMAX+1)**2)
     REAL(REALK),   INTENT(IN) :: T_matrix(:,:)
     REAL(REALK) :: mm_contract_Tq(SIZE(vect))
     
     INTEGER :: L,m, p,q,s, u, qmin,qmax


     ! do L=0 terms first
     p = (LMAX+1)*(LMAX+1)
!Simen Room for improvement here T_matrix is very sparse
     mm_contract_Tq(1) = half*DOT_PRODUCT(vect(1:p),T_matrix(1:p,1))
     DO s = 2, p 
        mm_contract_Tq(s) = vect(1)*T_matrix(s,1)
     END DO

     contract: DO L = 1, LMAX

        u = L*(L+1) +1
        p = (LMAX-L+1)*(LMAX-L+1)

        qmin = u-L
        qmax = MIN( u+L, p )

        qloop: DO q = qmin, qmax
           mm_contract_Tq(q) = mm_contract_Tq(q)           &
                               + DOT_PRODUCT(vect(q:p),T_matrix(q:p,q))
           DO s = q+1, p 
              mm_contract_Tq(s) = mm_contract_Tq(s) + T_matrix(s,q)*vect(q)
           END DO
        END DO qloop

        ! add on extra post-factors "after" contraction
        mm_contract_Tq(u) = half*mm_contract_Tq(u)    ! m=0

     END DO contract
     
   END FUNCTION mm_contract_Tq
   
!------------------------------------------------------------------------------

   SUBROUTINE mm_scale_vec(LMAX,scl_fact,scale_vec,prefactor)

      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: LMAX
      REAL(REALK),   INTENT(IN)  :: scl_fact 
      REAL(REALK),   INTENT(INOUT) :: scale_vec(:), prefactor
      
      INTEGER :: L,L2, lo,hi, u
      REAL(REALK)   :: prefactor1, prevec
    
      prefactor1 = one/scl_fact

      prevec       = one
      scale_vec(1) = prevec
      DO L = 1, LMAX
         prevec = prefactor1*prevec
         L2 = L*L
         lo = L2+1
         hi = L2+2*L+1
         DO u = lo, hi
            scale_vec(u) = prevec
         END DO
      END DO
      
      IF (scl_fact < zero) THEN
         prefactor = -prefactor1
      ELSE 
         prefactor =  prefactor1
      END IF
      
   END SUBROUTINE mm_scale_vec
 
!------------------------------------------------------------------------------

   SUBROUTINE mm_postfac_Vff(LMAX,Vff_tmp)

      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: LMAX
      REAL(REALK),   INTENT(INOUT) :: Vff_tmp(:)

      INTEGER :: l, u
      ! add extra postfactors  (for half**(delta(0,m))
      DO l = 0, LMAX
         u = l*(l+1) +1
         Vff_tmp(u) = half*Vff_tmp(u)
      END DO

   END SUBROUTINE mm_postfac_Vff

!------------------------------------------------------------------------------

END MODULE mm_T_worker


!==============================================================================

MODULE mm_T_worker_multi

   USE mm_global_paras_mod

   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: mm_get_SPLTSQ_T_matrices,   &
             mm_get_FLTSQ_T_matrices,    &
             mm_contract_multi_Tq

CONTAINS

!------------------------------------------------------------------------------
! Get sparse (l+j <= LMAX), lower-triangular, square-based T-matrix
!FIXME: would be nice to allow possibility of building rectangular-based too

   SUBROUTINE mm_get_SPLTSQ_T_matrices(ndim,LMAX,r_ab,T_matrix)

      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: ndim, LMAX
      REAL(REALK),   INTENT(IN)  :: r_ab(:,:)
      REAL(REALK),   INTENT(INOUT) :: T_matrix(:,:,:)

      REAL(REALK) :: I_sh(ndim,(1+LMAX)**2)

      CALL mm_generate_I(ndim,LMAX,r_ab,I_sh)
      CALL mm_generate_T(ndim,LMAX,.FALSE.,I_sh,T_matrix)

   END SUBROUTINE mm_get_SPLTSQ_T_matrices

!------------------------------------------------------------------------------
! Get full (all elements built), lower-triangular, square T-matrix

   SUBROUTINE mm_get_FLTSQ_T_matrices(ndim,LMAX,r_ab,T_matrix)

      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: ndim, LMAX
      REAL(REALK),   INTENT(IN)  :: r_ab(:,:)
      REAL(REALK),   INTENT(INOUT) :: T_matrix(:,:,:)

      REAL(REALK) :: I_sh(ndim,(1+2*LMAX)**2)
     
      CALL mm_generate_I(ndim,2*LMAX,r_ab,I_sh)
      CALL mm_generate_T(ndim,LMAX,.TRUE.,I_sh,T_matrix)

   END SUBROUTINE mm_get_FLTSQ_T_matrices

!------------------------------------------------------------------------------

   SUBROUTINE mm_generate_I(ndim,LMAX,r_ab,I_sh) 
 
      ! Subroutine to generate Scaled Irregular solid harmonics
      ! See page 416, "Molecular Electronic Structure Theory".
      !--------------------------------------------------------
      ! I_sh(M,L)  == COS component of I_sh(M,L)
      ! I_sh(-M,L) == SIN component of I_sh(-M,L)
      ! where L and M are positive
      ! Note that I_sh(0,L) == COS component of I_sh(0,L) !
      ! In fact we combine (l,m) into a single index and
      ! we assume the order 0, -1,1,+1, -2,-1,0,+1,+2 ...
      !--------------------------------------------------------

      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: LMAX, NDIM
      REAL(REALK),   INTENT(IN)  :: r_ab(:,:)
      REAL(REALK),   INTENT(INOUT) :: I_sh(:,:)

      REAL(REALK)   :: r_minus2(NDIM)
      INTEGER :: L,m, i,j,k,p,q,u, sign_L , n, p1,q1, m1
      REAL(REALK)   :: var1, var2, var3, rmin2

      r_minus2(:) = one/(r_ab(:,1)**2 + r_ab(:,2)**2 + r_ab(:,3)**2)

      IF (LMAX==0) THEN
         I_sh(:,1) = SQRT(r_minus2(:))
         RETURN
      END IF
         
      I_sh(:,1) =  SQRT(r_minus2(:))
      I_sh(:,2) = -r_ab(:,2)*r_minus2(:)*I_sh(:,1)
      I_sh(:,3) =  r_ab(:,3)*r_minus2(:)*I_sh(:,1)
      I_sh(:,4) = -r_ab(:,1)*r_minus2(:)*I_sh(:,1)

      ! Now iterate higher order terms
      sign_L = -1
      DO L = 2, LMAX

         sign_L = -sign_L 
         i = (L+1)*(L+1) 
         j = L*L+1 
         p = j-2*L+1
         q = j-1
         p1 = p 
         q1 = q 
         p = p+L-1 
         q = L*(L-3)+3
         u = q+L-2
         k = i-L
         m1 = L-1

         DO n = 1, ndim

            rmin2 = r_minus2(n)
            var1 = (2*L-1)*rmin2
            var2 = var1*r_ab(n,2)*sign_L
            var3 = var1*r_ab(n,1)
            I_sh(n,i) = var2*I_sh(n,p1) - var3*I_sh(n,q1)
            I_sh(n,j) = var2*I_sh(n,q1) + var3*I_sh(n,p1)
            var2 = var1*r_ab(n,3)

            DO m = 0, L-2
               var3 = (u-m*m)*rmin2
               I_sh(n,k+m) = var2*I_sh(n,p+m) - var3*I_sh(n,q+m)
               I_sh(n,k-m) = var2*I_sh(n,p-m) - var3*I_sh(n,q-m)  
            END DO

            I_sh(n,k+m1) = var2*I_sh(n,p+m1)
            I_sh(n,k-m1) = var2*I_sh(n,p-m1)

         end do

      END DO

   END SUBROUTINE mm_generate_I

!------------------------------------------------------------------------------

   SUBROUTINE mm_generate_T(ndim,LMAX,TOLMAX,I_sh,T_matrix)
 
   ! Subroutine to generate Real Interaction Matrix for 
   ! given point from Irregular Solid Harmonics (I_sh) 
   ! cf Helgaker et. al. pp 415
   !-----------------------------------------------------------------------
   !  T_cos_cos(lm,jk) = 2*{I_cos(l+j,m+k) + (-1)^k*I_cos(l+j,m-k)}  etc.
   ! 
   ! Note we only write the non-zero elements here, and the precise form
   ! depends on the function JMAX passed in
   !-----------------------------------------------------------------------
   !FIXME: can optimise for j=l by noting that if we build just lower
   ! triangular, we need only half the elements there

      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: LMAX, NDIM
      LOGICAL,       INTENT(IN)  :: TOLMAX
      REAL(REALK),   INTENT(IN)  :: I_sh(:,:)
      REAL(REALK),   INTENT(INOUT) :: T_matrix(:,:,:)

      INTEGER :: L,m,J,k, q, pp,qq, u, Jlast
      INTEGER :: sign_m

      IF (LMAX==0) THEN
         T_matrix(1:ndim,1,1) = two*I_sh(1:ndim,1)  ! NB scaling to make T symmetric
         RETURN
      END IF

      L_loop: DO L = 0, LMAX 
         qq = L*(L+1) +1

         IF (tolmax) THEN
            Jlast = LMAX
         ELSE
            Jlast = LMAX - L
         END IF

         J_loop: DO J = L, Jlast 
            pp = J*(J+1) +1
            u = 1+(L+J)*(L+J+1)

            sign_m = -1
            positive_m: DO m = 0, L 
               q = qq+m
               sign_m = -sign_m

               ! cos-cos terms
               !---------------

               DO k = 0, MIN(m,J)  
                  IF (BTEST(k,0)) THEN
                     T_matrix(1:ndim,pp+k,q) = I_sh(1:ndim,u+m+k) - I_sh(1:ndim,u+m-k)
                  ELSE
                     T_matrix(1:ndim,pp+k,q) = I_sh(1:ndim,u+m+k) + I_sh(1:ndim,u+m-k)
                  END IF
               END DO
               IF (sign_m > 0) THEN
                  DO k = m+1, J
                     T_matrix(1:ndim,pp+k,q) = I_sh(1:ndim,u+m+k) + I_sh(1:ndim,u-m+k)
                  END DO
               ELSE
                  DO k = m+1, J
                     T_matrix(1:ndim,pp+k,q) = I_sh(1:ndim,u+m+k) - I_sh(1:ndim,u-m+k)
                  END DO
               END IF

               ! cos-sin terms
               !---------------

               IF (sign_m > 0) THEN
                  DO k = -J, MIN(-m-1,-1)
                     T_matrix(1:ndim,pp+k,q) = I_sh(1:ndim,u+(m+k)) + I_sh(1:ndim,u+(-(m-k))) 
                  END DO
                  DO k = MAX(-J,-m), MIN(-m,-1) 
                     T_matrix(1:ndim,pp+k,q) = I_sh(1:ndim,u+(-(m-k)))
                  END DO
                  DO k = -1, MAX(-m+1,-J), -1 
                     IF (BTEST(k,0)) THEN
                        T_matrix(1:ndim,pp+k,q) =   I_sh(1:ndim,u+(-(m+k)))    &
                                             + I_sh(1:ndim,u+(-(m-k)))
                     ELSE
                        T_matrix(1:ndim,pp+k,q) = - I_sh(1:ndim,u+(-(m+k)))    &
                                             + I_sh(1:ndim,u+(-(m-k)))
                     END IF
                  END DO
               ELSE
                  DO k = -J, MIN(-m-1,-1) 
                     T_matrix(1:ndim,pp+k,q) = I_sh(1:ndim,u+(m+k)) - I_sh(1:ndim,u+(-(m-k))) 
                  END DO
                  DO k = MAX(-J,-m), MIN(-m,-1)
                     T_matrix(1:ndim,pp+k,q) = - I_sh(1:ndim,u+(-(m-k)))
                  END DO
                  DO k = -1, MAX(-m+1,-J), -1 
                     IF (BTEST(k,0)) THEN
                        T_matrix(1:ndim,pp+k,q) = - I_sh(1:ndim,u+(-(m+k)))   &
                                             - I_sh(1:ndim,u+(-(m-k)))
                     ELSE
                        T_matrix(1:ndim,pp+k,q) =   I_sh(1:ndim,u+(-(m+k)))   &
                                             - I_sh(1:ndim,u+(-(m-k)))
                     END IF
                  END DO
               END IF
 
            END DO  positive_m
 
            sign_m = 1
            negative_m: DO m = -1, -L, -1  
               q = qq+m
               sign_m = -sign_m
 
               ! sin-cos terms
               !---------------

               k = MIN(J,-m)
               IF (BTEST(k,0)) THEN
                  T_matrix(1:ndim,pp+k,q) = - I_sh(1:ndim,u+(m-k))
               ELSE
                  T_matrix(1:ndim,pp+k,q) =   I_sh(1:ndim,u+(m-k))
               END IF

               DO k = 0, MIN((-m-1),J)
                  IF (BTEST(k,0)) THEN
                     T_matrix(1:ndim,pp+k,q) = I_sh(1:ndim,u+(m+k)) - I_sh(1:ndim,u+(m-k))
                  ELSE
                     T_matrix(1:ndim,pp+k,q) = I_sh(1:ndim,u+(m+k)) + I_sh(1:ndim,u+(m-k))
                  END IF
               END DO

               IF (sign_m > 0) THEN
                  DO k = J, (-m+1), -1 
                     IF (BTEST(k,0)) THEN
                        T_matrix(1:ndim,pp+k,q) =   I_sh(1:ndim,u+(-(m+k)))   &
                                             - I_sh(1:ndim,u+(m-k))
                     ELSE
                        T_matrix(1:ndim,pp+k,q) = - I_sh(1:ndim,u+(-(m+k)))   &
                                             + I_sh(1:ndim,u+(m-k))
                     END IF
                  END DO
               ELSE
                  DO k = J, (-m+1), -1 
                     IF (BTEST(k,0)) THEN
                        T_matrix(1:ndim,pp+k,q) = - I_sh(1:ndim,u+(-(m+k)))   &
                                             - I_sh(1:ndim,u+(m-k))
                     ELSE
                        T_matrix(1:ndim,pp+k,q) =   I_sh(1:ndim,u+(-(m+k)))   &
                                             + I_sh(1:ndim,u+(m-k))
                     END IF
                  END DO
               END IF

               ! sin-sin terms
               !---------------

               IF (sign_m > 0) THEN
                  DO k = -J, m 
                     IF (BTEST(k,0)) THEN
                        T_matrix(1:ndim,pp+k,q) =   I_sh(1:ndim,u+(-(m+k)))   &
                                             - I_sh(1:ndim,u+(m-k))
                     ELSE
                        T_matrix(1:ndim,pp+k,q) = - I_sh(1:ndim,u+(-(m+k)))   &
                                             + I_sh(1:ndim,u+(m-k))
                     END IF
                  END DO
                  DO k = -1, MAX(-J,m+1), -1 
                     IF (BTEST(k,0)) THEN
                        T_matrix(1:ndim,pp+k,q) =   I_sh(1:ndim,u+(-(m+k)))   &
                                             + I_sh(1:ndim,u+(-(m-k)))
                     ELSE
                        T_matrix(1:ndim,pp+k,q) = - I_sh(1:ndim,u+(-(m+k)))   &
                                             + I_sh(1:ndim,u+(-(m-k)))
                     END IF
                  END DO
               ELSE
                  DO k = -J, m 
                     IF (BTEST(k,0)) THEN
                        T_matrix(1:ndim,pp+k,q) = - I_sh(1:ndim,u+(-(m+k)))   &
                                             - I_sh(1:ndim,u+(m-k))
                     ELSE
                        T_matrix(1:ndim,pp+k,q) =   I_sh(1:ndim,u+(-(m+k)))   &
                                             + I_sh(1:ndim,u+(m-k))
                     END IF
                  END DO
                  DO k = -1, MAX(-J,m+1), -1 
                     IF (BTEST(k,0)) THEN
                        T_matrix(1:ndim,pp+k,q) = - I_sh(1:ndim,u+(-(m+k)))   &
                                             - I_sh(1:ndim,u+(-(m-k)))
                     ELSE
                        T_matrix(1:ndim,pp+k,q) =   I_sh(1:ndim,u+(-(m+k)))   &
                                             - I_sh(1:ndim,u+(-(m-k)))
                     END IF
                  END DO
               END IF

            END DO  negative_m
         END DO     J_loop
      END DO        L_loop

   END SUBROUTINE mm_generate_T

!------------------------------------------------------------------------------

   FUNCTION mm_contract_multi_Tq(LMAX,vect,T_mats,ndim)

     IMPLICIT NONE
     INTEGER, INTENT(IN) :: LMAX, ndim
     REAL(REALK),   INTENT(IN) :: vect((LMAX+1)**2)
     REAL(REALK),   INTENT(IN) :: T_mats(:,:,:)
     REAL(REALK) :: mm_contract_multi_Tq(ndim,(LMAX+1)**2)
     
     INTEGER :: L,m, p,q,r,s ,u, qmin,qmax
     REAL(REALK) :: fac 

     ! first do L=0 terms
     p = (LMAX+1)*(LMAX+1)

     fac = half*vect(1)
     DO r = 1, ndim
        mm_contract_multi_Tq(r,1) = fac*T_mats(r,1,1)
     END DO

     DO s = 2, p
     fac = half*vect(s)
        DO r = 1, ndim
           mm_contract_multi_Tq(r,1) = mm_contract_multi_Tq(r,1)    &
                                       + fac*T_mats(r,s,1)
        END DO
     END DO

     fac = vect(1)
     DO s = 2, p
        DO r = 1, ndim
           mm_contract_multi_Tq(r,s) = fac*T_mats(r,s,1)
        END DO
     END DO

     contract: DO L = 1, LMAX

        u = L*(L+1) +1
        p = (LMAX-L+1)*(LMAX-L+1)

        qmin = u-L 
        qmax = MIN( u+L, p ) 

        qloop: DO q = qmin, qmax

           DO s = q, p
           CALL DAXPY(ndim,vect(s),T_mats(1,s,q),1,mm_contract_multi_Tq(1,q),1)
           END DO

           fac = vect(q)
           DO s = q+1, p
              DO r = 1, ndim
                 mm_contract_multi_Tq(r,s) = mm_contract_multi_Tq(r,s)    &
                                             + fac*T_mats(r,s,q)
              END DO
           END DO

        END DO qloop

        mm_contract_multi_Tq(:,u) = half*mm_contract_multi_Tq(:,u)    ! m=0

     END DO contract
     
   END FUNCTION mm_contract_multi_Tq
   
!------------------------------------------------------------------------------

END MODULE mm_T_worker_multi


!************************************************************
!* These subroutines are placed outside the module in order *
!* to be accessible from Fortran 77 code.                   *
!************************************************************

! This subroutine computes the local moments loc_q resulting from
! the m2l operation represented by Tmatrix and the multipole moments
! mm_q.
!
! Written by Erik Tellgren, October 2005
!
subroutine pbc_contract_Tq(loc_q,siz_q,Tmatrix,siz_T,mm_q,lmax_q)
  use mm_global_paras_mod
  use mm_T_worker, only: mm_contract_Tq
  implicit none
  integer, intent(in) :: lmax_q, siz_q, siz_T
  real(realk), intent(in) :: Tmatrix(siz_T,siz_T), mm_q(siz_q)
  real(realk), intent(inout) :: loc_q(siz_q)

  integer :: i
  real(realk) :: dummy(siz_q)

  if (siz_q .ne. (1+lmax_q)**2) call LSquit('pbc_contract_Tq: Inconsistent sizes.',-1)

  if (siz_q .gt. siz_T) call LSquit('pbc_contract_Tq: Qlmax larger than TWlmax.',-1)

  loc_q(:) = mm_contract_Tq(lmax_q,mm_q,Tmatrix)

end subroutine pbc_contract_Tq


! This is a wrapper for the interaction (T) matrix generating routine.
! It's only function is to enable calls to mm_get_FLTSQ_T_matrix from
! the PBC code.
!
! Written by Erik Tellgren, October 2005
!
SUBROUTINE pbc_get_FLTSQ_T_matrix(LMAX,r_ab,T_matrix)
  use mm_global_paras_mod
  use mm_T_worker, only : mm_get_FLTSQ_T_matrix
  IMPLICIT NONE
  INTEGER, INTENT(in)  :: LMAX
  REAL(REALK),   INTENT(in)  :: r_ab(3)
  REAL(REALK), INTENT(inout) :: T_matrix((1+LMAX)**2,(1+LMAX)**2)

  integer :: l, j, qmin, qmax, siz

  call mm_get_FLTSQ_T_matrix(LMAX,r_ab,T_matrix)

END SUBROUTINE pbc_get_FLTSQ_T_matrix

! The FMM code generates "symmetrized" T tensors (see Eq. (67) in
! Watson et al. J Chem Phys 121:2915) with an additional factor for
! the m = 0 elements. The relation between usual T tensors and the
! ones supplied by the FMM code is
!
! Tfmm(lm,jk) = Tusual(lm,jk) (1+\delta_{m0}) (-1)^j, for jk <= lm,
! Tfmm(lm,jk) = 0,                                    for jk > lm.
!
! In this subroutine we restore the usual T tensor.
!
! Written by Erik Tellgren, November 2005
!
subroutine pbc_restore_T_matrix(lmax,T_matrix)
  use mm_global_paras_mod
  implicit none

  integer, intent(in) :: lmax
  real(realk), intent(inout) :: T_matrix((1+lmax)**2,(1+lmax)**2)

  integer :: l,m,j,k,lm,jk
  real(realk) :: fac

  write(LUPRI,*)
  write(LUPRI,*) 'Transforming a T tensor from triangular to square form.'
  write(LUPRI,*)

  rows_l: do l = lmax,0,-1
     rows_m: do m = l,-l,-1
        lm = l*(l+1)+1+m
        cols_j: do j = lmax,0,-1
           cols_k: do k = j,-j,-1
              jk = j*(j+1)+1+k

              fac = 2.0D0
              if (m .eq. 0) fac = fac / 2.0D0
              if (k .eq. 0) fac = fac / 2.0D0
              if (mod(j,2) .eq. 1) fac = -fac

              if (jk .le. lm) then
                 T_matrix(lm,jk) = T_matrix(lm,jk) * fac
              else if (mod(j-l,2) .eq. 1) then
                 T_matrix(lm,jk) = -T_matrix(jk,lm)
              else
                 T_matrix(lm,jk) = T_matrix(jk,lm)
              end if
           end do cols_k
        end do cols_j

     end do rows_m
  end do rows_l

end subroutine pbc_restore_T_matrix

subroutine pbc_generate_I(LMAX,pos,I_sh) 
  use mm_global_paras_mod
  use mm_T_worker, only : mm_generate_I
  implicit none
 
  integer, intent(in)  :: lmax
  real(realk),   intent(in)  :: pos(3)
  real(realk),intent(inout)  :: I_sh((lmax+1)**2)

  I_sh(:) = 0.0D0
  call mm_generate_I(lmax,pos,I_sh)

end subroutine pbc_generate_I
