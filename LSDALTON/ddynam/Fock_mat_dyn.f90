Module Fock_MD
! Do Fock matrix dynamics
  Use precision
  Use ls_util
  Use ls_dynamics
  Use matrix_module
  Use matrix_operations
Contains
!================================!
! Subroutine FMD_run             !
!================================!
! Basic procedure for WF dynamics with constant time-step
Subroutine FMD_run(traj,F,StepNum,NPoints,PolyOrd,Start)
!
Implicit none
Type(trajtype) :: traj
Type(Matrix), intent(inout) :: F
Integer :: StepNum, NPoints,PolyOrd
Logical :: Start
Integer :: i
! Fill in Fock_array
If (StepNum .LT. (NPoints)) then
   traj%Fock_array(StepNum+1) = F
Else
!   Do circular shift
!
!   For some reason EOSHIFT does not work here for the pointer of type(matrix)!
!   traj%Fock_array = EOSHIFT(traj%Fock_array,shift=1,boundary = F,dim=1)
    Do i = 1, (NPoints-1) 
       traj%Fock_array(i) = traj%Fock_array(i+1)
    Enddo
    traj%Fock_array(NPoints) = F 
Endif 
! Determine FMD coefficients
If (StepNum .EQ. NPoints-1) Call Calc_FMD_coef(NPoints,PolyOrd,traj%FMD_coef)
! Start extrapolation
If (StepNum .GE. (NPoints-1)) then
   Write(*,*)'Fock matrix will be extrapolated'
   Start = .TRUE.
   Call Fock_mat_dyn_const(NPoints,traj%Fock_array,traj%FMD_coef,F)
Endif
!
End subroutine FMD_run
!================================!
! Subroutine Calc_FMD_coef       !
!================================!
! Calculates expansion coefficients for Fock mat.extrapolation
! for constant time-step 
Subroutine Calc_FMD_coef(NPoints,PolyOrd,FMD_coef)
Implicit none
Integer :: NPoints, PolyOrd
Real(realk), dimension(NPoints) :: FMD_coef
Real(realk), dimension(NPoints,PolyOrd+1) :: A
Real(realk), dimension(PolyOrd+1,NPoints) :: AAA          ! Supplementary matrix
Real(realk), dimension(PolyOrd+1,PolyOrd+1) :: AA, InvAA  ! Other supplementary matrices
Integer :: i,j,p,q
Integer :: ErrorFlag
! Forming A matrix
A = 0.D0
AA = 0.D0
AAA = 0.D0
InvAA = 0.D0
Do i = 1, NPoints
   Do j = 1, PolyOrd + 1
      A(i,j) = (i - 0.5*(NPoints + 1) )**(j-1)
   Enddo
Enddo 
! Calculating FMD coefficients
If (NPoints  .EQ. (PolyOrd+1) ) then
   Call FINDInv(A,AAA,(PolyOrd+1),ErrorFlag)
Else
   AA = MATMUL(TRANSPOSE(A),A)
   Call FINDInv(AA,InvAA,(PolyOrd+1),ErrorFlag)
   AAA = MATMUL( InvAA, TRANSPOSE(A) )
Endif
FMD_coef = 0.D0
Do p = 1, Npoints
   Do q = 1, (PolyOrd+1)
      FMD_coef(p) = FMD_coef(p) + AAA(q,p)*(0.5D0*NPoints + 0.5D0)**(q-1)
   Enddo
Enddo       
!
End subroutine Calc_FMD_coef
!================================!
! Subroutine Fock_mat_dyn_const  !
!================================!
!
! Calculates the extrapolated Fock element(s+1)
!
Subroutine Fock_mat_dyn_const(NPoints,Fock_array,FMD_coef,F)
Implicit none
Integer :: NPoints, i
Type(matrix), dimension(NPoints), intent(in) :: Fock_array
Real(realk), dimension(NPoints), intent(in) :: FMD_coef   ! FMD coefficients
Type(Matrix), intent(inout) :: F
! 
Call mat_zero(F)
Write(*,*) 'coeff=',FMD_coef
Do i = 1, NPoints
   Call mat_daxpy(FMD_coef(i),Fock_array(i),F)
Enddo
!
End subroutine Fock_mat_dyn_const
!===============================================!
!  Subroutine FINdInv                           !
!===============================================!
! Finds inverse matrices
!
Subroutine FINDInv(matrix, inverse, n, ErrorFlag)
Implicit none
Integer :: n
Integer :: ErrorFlag  !Return error status. -1 for error, 0 for normal
Real(realk), dimension(n,n) :: matrix  !Input matrix
Real(realk), dimension(n,n) :: inverse !Inverted matrix
Logical :: FLAG = .TRUE.
Integer :: i, j, k, l
Real(realk) :: m
Real(realk), dimension(n,2*n) :: augmatrix !augmented matrix
!
!Augment input matrix with an identity matrix
	DO i = 1, n
		DO j = 1, 2*n
			IF (j <= n ) THEN
				augmatrix(i,j) = matrix(i,j)
			ELSE IF ((i+n) == j) THEN
				augmatrix(i,j) = 1
			Else
				augmatrix(i,j) = 0
			ENDIF
		END DO
	END DO
	
	!Reduce augmented matrix to upper traingular form
	DO k =1, n-1
		IF (augmatrix(k,k) == 0) THEN
			FLAG = .FALSE.
			DO i = k+1, n
				IF (augmatrix(i,k) /= 0) THEN
					DO j = 1,2*n
						augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
					END DO
					FLAG = .TRUE.
					EXIT
				ENDIF
				IF (FLAG .EQV. .FALSE.) THEN
					PRINT*, "Matrix is non - invertible"
					inverse = 0
					errorflag = -1
					return
				ENDIF
			END DO
		ENDIF
		DO j = k+1, n			
			m = augmatrix(j,k)/augmatrix(k,k)
			DO i = k, 2*n
				augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
			END DO
		END DO
	END DO
	
	!Test for invertibility
	DO i = 1, n
		IF (augmatrix(i,i) == 0) THEN
			PRINT*, "Matrix is non - invertible"
			inverse = 0
			errorflag = -1
			return
		ENDIF
	END DO
	
	!Make diagonal elements as 1
	DO i = 1 , n
		m = augmatrix(i,i)
		DO j = i , (2 * n)				
			   augmatrix(i,j) = (augmatrix(i,j) / m)
		END DO
	END DO
	
	!Reduced right side half of augmented matrix to identity matrix
	DO k = n-1, 1, -1
		DO i =1, k
		m = augmatrix(i,k+1)
			DO j = k, (2*n)
				augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
	   	        END DO
		END DO
	END DO				
	
	!store answer
	DO i =1, n
		DO j = 1, n
			inverse(i,j) = augmatrix(i,j+n)
		END DO
	END DO
	errorflag = 0
END Subroutine FINDinv
!
End module Fock_MD
