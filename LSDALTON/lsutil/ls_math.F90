MODULE math_fun
use precision

public :: ls_vdexp,ls_vdinv,ls_vdsqrt,FACULT,FACUL2,BINOM,NCRT
private
contains
subroutine ls_vdinv(N,A,invA)
implicit none
integer,intent(in) :: N
real(realk),intent(in) :: A(N)
real(realk),intent(inout) :: invA(N)
!
!#ifdef VAR_MKL
!call vdinv(N,A,invA)
!#else
integer :: i
DO i=1, N
   invA(i) = 1.0E0_realk/A(i)
ENDDO
!#endif
end subroutine ls_vdinv

subroutine ls_vdsqrt(N,A,sqrtA)
implicit none
integer,intent(in) :: N
real(realk),intent(in) :: A(N)
real(realk),intent(inout) :: sqrtA(N)
!
!#ifdef VAR_MKL
!   call vdsqrt(N,A,sqrtA)
!#else
integer :: i
DO i=1,N
   sqrtA(i) = SQRT(A(i))
ENDDO
!#endif
end subroutine ls_vdsqrt

subroutine ls_vdexp(N,A,expA)
implicit none
integer,intent(in) :: N
real(realk),intent(in) :: A(N)
real(realk),intent(inout) :: expA(N)
!
!#ifdef VAR_MKL
!   call vdexp(N,A,expA)
!#else
integer :: i
DO i=1,N
   expA(i) = EXP(A(i))
ENDDO
!#endif
end subroutine ls_vdexp

FUNCTION FACULT(LUPRI,N)
IMPLICIT NONE
real(realk), PARAMETER :: D1=1E0_realk
integer        :: N,I,LUPRI
real(realk)    :: FACULT
IF (N .LT. 0) THEN
   WRITE (LUPRI,'(/,A,I10,/A)')&
   &         ' Argument less than zero in FACULT:',N,&
   &         ' Program cannot continue.'
   CALL LSQUIT('Illegal argument in FACULT',lupri)
ELSE
   FACULT = D1
   DO I = 1, N
      FACULT = FACULT*I
   ENDDO
END IF
END FUNCTION FACULT

FUNCTION FACUL2(LUPRI,N)
IMPLICIT NONE
real(realk), PARAMETER :: D1=1E0_realk
real(realk)    :: FACUL2
integer :: N,I,LUPRI
IF (N .LT. 0) THEN
   FACUL2 = DFLOAT(N + 2)
   DO I = N + 4, 1, 2
      FACUL2 = FACUL2*I
   END DO
   IF (FACUL2 .EQ. 0E0_realk) THEN
      WRITE (LUPRI,'(/,A,I10,/A)')&
      &            ' Double factorial undefined for ',N,&
      &            ' Program cannot continue.'
      CALL LSQUIT('Illegal argument in FACUL2',lupri)
   ELSE
      FACUL2 = D1/FACUL2
   END IF
ELSE IF (N.EQ. 0) THEN
   FACUL2 = D1
ELSE ! N > 0
   FACUL2 = DFLOAT(N)
   DO I = N - 2, 1, -2
      FACUL2 = FACUL2*I
   END DO
END IF
END FUNCTION FACUL2

FUNCTION BINOM(LUPRI,I,J)
real(realk), PARAMETER  :: D1=1E0_realk
INTEGER :: I,J,LUPRI
real(realk)    :: BINOM

IF (I .LT. J) THEN
   WRITE (LUPRI,'(/,A,2I5,/A)')&
   &         ' Second argument larger than first argument in BINOM:',&
   &         I,J,' Program cannot continue.'
   CALL LSQUIT('Illegal arguments in BINOM',lupri)
ELSE
   BINOM = FACULT(LUPRI,I)/(FACULT(LUPRI,I-J)*FACULT(LUPRI,J))
END IF
END FUNCTION BINOM

FUNCTION NCRT(I,J,K)
IMPLICIT NONE
INTEGER  :: I,J,K,NCRT
NCRT = 1 + J + 2*K + (J + K)*(J + K - 1)/2
END FUNCTION NCRT

END MODULE math_fun


