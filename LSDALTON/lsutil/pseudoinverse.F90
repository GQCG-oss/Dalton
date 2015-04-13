!> @file
!> Contains pseudo inverse routine
MODULE pseudoinverseMod
  use precision  
  use memory_handling

CONTAINS
  !> \brief pseudo inverse routine
  !> n rows of S
  !> m columns of S
  !> S the non square matrix to inverse (Output: the inverse of input)
  SUBROUTINE PSEUDOINVERSE(n,m,S,epsilon)
  implicit none
  integer,intent(in) :: n,m
  real(realk),intent(inout) :: S(n,m)
  real(realk) :: epsilon
  !local variables
  integer                :: nmin,lwork,INFO,I,K,nn
  real(realk), pointer   :: work(:),SV(:),U(:,:),VT(:,:)
  integer,pointer        :: IPVT(:)
  real(realk)            :: RCOND, dummy(2),maxSV,svmt

  !Perform a SVD  decomposition 
  ! S = U * SIGMA * transpose(V)
  ! where SIGMA is an M-by-N matrix which is zero except for 
  ! its min(m,n) diagonal elements
  ! for efficient storage the SIGMA non zero elements are stored in 
  ! SV (non singular values)
  nmin = MIN(n,m)
  call mem_alloc(SV,nmin)   !sigular values 
  !only the first min(m,n) columns of U (the left singular
  !vectors) are returned in the array U;
  call mem_alloc(U,n,nmin)  
  !only the first min(m,n) rows of V**T (the left singular
  !vectors) are returned in the array U;
  call mem_alloc(VT,nmin,m)
  !S(n,m) = U(n,nmin) SV(nmin) VT(nmin,m)
  lwork = -1      !workspace query
  call dgesvd('S','S',n,m,S,n,SV,U,n,VT,nmin,dummy,lwork,INFO)
  lwork = dummy(1)
  call mem_alloc(work,lwork)
  call dgesvd('S','S',n,m,S,n,SV,U,n,VT,nmin,work,lwork,INFO)
  !content of S destroyed
  IF(INFO.NE.0)THEN
     print*,'dgesvd in PSEUDOINVERSE failed  INFO=',INFO
  ENDIF
  call mem_dealloc(work)
  maxSV = MAXVAL(SV)

  DO I=1,nmin
     IF(SV(I).GT.epsilon*maxSV)THEN
        svmt = 1.0E0_realk/SV(I)
        DO K=1,n
           U(K,I)=U(K,I)*svmt
        ENDDO
     ELSE
        DO K=1,n
           U(K,I)=0.0E0_realk
        ENDDO
     ENDIF
  ENDDO
  call mem_dealloc(SV)
  !S^-1(n,m) = [U(n,nmin) SV(nmin)] * VT(nmin,m)
  call dgemm('N','N',N,M,NMIN,1.0E0_realk,U,N,VT,NMIN,0.0E0_realk,S,N)

  call mem_dealloc(VT)
  call mem_dealloc(U)

end SUBROUTINE PSEUDOINVERSE

END MODULE PSEUDOINVERSEMOD
