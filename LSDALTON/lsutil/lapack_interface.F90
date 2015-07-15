!> @file 
!> \brief Contains LAPACK interface subroutines.
!>
!> These routines depend only on the module precision. 
!> UNDER NO CIRCUMSTANCES ARE YOU ALLOWED TO INTRODUCE OTHER DEPENDENCIES!
!> If your routine depends on other modules than precision, it means
!> that it belongs somewhere else than here.
!> The reason for this subroutine is mostly to avoid issues with 64 bit integers
!> and other things 
!>
!> The problem with 32/64 bit inconsistencies (assuming all integers can be treated with 32 bit integers)
!> is that when you link to 32-bit libraries, but have compiled LSDalton with 64-bit integer, then only
!> the lower 32 bits are initialized and used inside for instance DPOTRF and the rest becomes arbitrary, 
!> leading to your INFO .ne. 0. 
!> This is not a problem for the other parameters as they are intent in, and the top 32 bits are not used 
!> for them (no dimensions greater than 2G ;-) ).
!> 
!> linking 64-bit LSDALTON with 32-bit libraries will work for all intent in and intent inout variables, 
!> as long as the integer does not exceeds 32 bits. 
!> Pure intent out variables will NOT work, unless they are initialized to zero before the call.
!> 
!> This interface is the proper way of doing things but you could also just set you 64 bit INFO to zero
!> before the call to your 32 bit lapack
module lapackMod
  use precision

  contains

   !  DPOTRF computes the Cholesky factorization of a real symmetric
   !  positive definite matrix A.
   subroutine LSDPOTRF(UPLO, N, A, LDA, INFO )
     implicit none
     character,intent(in)      :: UPLO !'U':  Upper triangle of A is stored
                                       !'L':  Lower triangle of A is stored.
     Integer,intent(in)        :: N    !The order of the matrix A.  N >= 0
     real(realk),intent(inout) :: A(LDA*N)
     Integer,intent(in)        :: LDA  !first dimension of matrix A
     Integer,intent(inout)     :: INFO != 0:  successful exit
     integer(kind=INTLAPACKK)  :: INFO2
     IF(INTLAPACKK.EQ.INTK)THEN
        call DPOTRF(UPLO, N, A, LDA, INFO )
     ELSE
        call DPOTRF(UPLO, int(N,kind=INTLAPACKK),A,int(LDA,kind=INTLAPACKK),INFO2)
        INFO = int(INFO2,kind=INTK)
     ENDIF
   end subroutine LSDPOTRF
   
   !  DPOTRI computes the inverse of a real symmetric positive definite
   !  matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
   !  computed by DPOTRF.
   subroutine LSDPOTRI(UPLO, N, A, LDA, INFO )
     implicit none
     character,intent(in)      :: UPLO !'U':  Upper triangle of A is stored
                                       !'L':  Lower triangle of A is stored.
     Integer,intent(in)        :: N    !The order of the matrix A.  N >= 0
     real(realk),intent(inout) :: A(LDA*N)
     Integer,intent(in)        :: LDA  !first dimension of matrix A
     Integer,intent(inout)     :: INFO != 0:  successful exit
     integer(kind=INTLAPACKK)  :: INFO2
     IF(INTLAPACKK.EQ.INTK)THEN
        call DPOTRI(UPLO, N, A, LDA, INFO )
     ELSE
        call DPOTRI(UPLO, int(N,kind=INTLAPACKK),A,int(LDA,kind=INTLAPACKK),INFO2)
        INFO = int(INFO2,kind=INTK)
     ENDIF
   end subroutine LSDPOTRI

 end module lapackMod

