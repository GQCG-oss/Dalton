!> @file 
!> Contains the precision specifications
MODULE precision
#ifdef SYS_REAL
  INTEGER, PARAMETER :: realk = 4
#else
  INTEGER, PARAMETER :: realk = 8
#endif
#ifdef VAR_REAL_SP
  INTEGER, PARAMETER :: real_sp = 4
#endif
  INTEGER, PARAMETER :: real_pt = 8
! Long integer is equivalent to 64 bit integers (8 bytes) and ranges 
! from -9223372036854775808 to 9223372036854775807 (i.e., -2**63 to 2**63-1)
  integer, parameter :: long = 8
!integer defining the size of the MPI input (depending on wheter the mpi library
!is compiled with 64 or 32 bit integers)
#ifdef VAR_INT64
#ifdef VAR_MPI_32BIT_INT
  integer, parameter :: ls_mpik = 4
#else
  integer, parameter :: ls_mpik = 8
#endif
#else
  integer, parameter :: ls_mpik = 4
#endif
! Short integer is equivalent to 8 bit integers (1 byte) and ranges 
! from -128 to 127 (i.e., -2**7 to 2**7-1)
  INTEGER,PARAMETER :: complexk = kind((1D0,1D0))
  integer, parameter :: short = 1
  integer, parameter :: shortzero = -33
  real(realk),parameter :: shortintCRIT=1E1_realk**shortzero
!Beware: if you add 4 zeros together you get -133 and therefore integer overflow.
  ! Quick and dirty conversion of 32 bit integgers to 64 bit integers:
  ! Multiply 32 bit integer by i8
  integer(kind=8),parameter :: i8 = 1
#ifdef VAR_INT64
  INTEGER, PARAMETER :: INTK = SELECTED_INT_KIND(17)!, REALK = KIND(1D0)
#else
  INTEGER, PARAMETER :: INTK = SELECTED_INT_KIND(9)!, REALK = KIND(1D0)
#endif
  !maximum number of elements that can be held in a signed 32bit integer
  !2,147,483,648  but some compilers complain, so 2,147,483,640 is used instead
  INTEGER, PARAMETER :: MAXINT32 = 2147483640
  !the same is done for 64bit integers which can hold 9,223,372,036,854,775,808 
#ifdef VAR_INT64
  INTEGER, PARAMETER :: MAXINT = 9223372036854775800
#else
  INTEGER, PARAMETER :: MAXINT = 2147483640
#endif
  
INTERFACE Test_if_64bit_integer_required
   MODULE PROCEDURE Test_if_64bit_integer_required_2d,&
        & Test_if_64bit_integer_required_3d,&
        & Test_if_64bit_integer_required_4d 
end INTERFACE Test_if_64bit_integer_required

contains

subroutine Test_if_64bit_integer_required_2d(n1,n2)
   implicit none
   integer,intent(in) :: n1,n2
#ifndef VAR_INT64
   !local variables
   integer(kind=long) :: n   
   n = n1*n2
   IF(n.GT.MAXINT)THEN
      print*,'A 64 bit integer is required in this context'
      print*,'The two dimensions: '
      print*,'n1 = ',n1
      print*,'n2 = ',n2
      print*,'gives a combined dimension',n
      print*,'which is bigger than can be described using a 32 bit integer'
      print*,'The maximum size is ',MAXINT
      print*,'Use the flag --int64 in the setup command e.g.'
      print*,'./setup --int64 BuildUsing64int'
      call lsquit('A 64 bit integer is required, recompile using --int64',-1)
   ENDIF
#endif
 end subroutine Test_if_64bit_integer_required_2d

subroutine Test_if_64bit_integer_required_3d(n1,n2,n3)
   implicit none
   integer,intent(in) :: n1,n2,n3
#ifndef VAR_INT64
   !local variables
   integer(kind=long) :: n   
   n = n1*n2*n3
   IF(n.GT.MAXINT)THEN
      print*,'A 64 bit integer is required in this context'
      print*,'The three dimensions: '
      print*,'n1 = ',n1
      print*,'n2 = ',n2
      print*,'n3 = ',n3
      print*,'gives a combined dimension',n
      print*,'which is bigger than can be described using a 32 bit integer'
      print*,'The maximum size is ',MAXINT
      print*,'Use the flag --int64 in the setup command e.g.'
      print*,'./setup --int64 BuildUsing64int'
      call lsquit('A 64 bit integer is required, recompile using --int64',-1)
   ENDIF
#endif
 end subroutine Test_if_64bit_integer_required_3d

subroutine Test_if_64bit_integer_required_4d(n1,n2,n3,n4)
   implicit none
   integer,intent(in) :: n1,n2,n3,n4
#ifndef VAR_INT64
   !local variables
   integer(kind=long) :: n   
   n = n1*n2*n3*n4
   IF(n.GT.MAXINT)THEN
      print*,'A 64 bit integer is required in this context'
      print*,'The four dimensions: '
      print*,'n1 = ',n1
      print*,'n2 = ',n2
      print*,'n3 = ',n3
      print*,'n4 = ',n4
      print*,'gives a combined dimension',n
      print*,'which is bigger than can be described using a 32 bit integer'
      print*,'The maximum size is ',MAXINT
      print*,'Use the flag --int64 in the setup command e.g.'
      print*,'./setup --int64 BuildUsing64int'
      call lsquit('A 64 bit integer is required, recompile using --int64',-1)
   ENDIF
#endif
 end subroutine Test_if_64bit_integer_required_4d

 subroutine Obtain_CS_THRLOG(CS_THRLOG,IntegralThreshold)
   implicit none
   integer(kind=short),intent(inout) :: CS_THRLOG
   real(realk),intent(in) :: IntegralThreshold
   !local variables
   integer :: IP
   real(realk) :: TMP
   !Beware when converting from double precision to short integer 
   !If double precision is less than 10^-33 then you can run into
   !problems with short integer overflow
   IF(IntegralThreshold.GT.shortintCrit)then
      TMP=log10(IntegralThreshold)
      IP=NINT(TMP) !IP used as temp
      IF (ABS(TMP-IP).LT. 1.0E-15_realk) THEN
         !this means that the threshold is 1EX_realk X=-1,..,-19,..
         CS_THRLOG=IP
      ELSE
         CS_THRLOG=FLOOR(log10(IntegralThreshold))
      ENDIF
   ELSE
      CS_THRLOG=shortzero
   ENDIF
 end subroutine Obtain_CS_THRLOG
 
END MODULE precision
