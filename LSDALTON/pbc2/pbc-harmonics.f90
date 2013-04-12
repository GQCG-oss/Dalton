MODULE harmonics_pbc
use precision
USE matrix_module
USE pbc_matrix_operations
IMPLICIT NONE
CONTAINS

!This subroutine computes clebsch-gordan coefficients
SUBROUTINE clebsch_g(lambda,mu,ell,emm,clebsch)
IMPLICIT NONE
INTEGER, INTENT(IN):: lambda, mu,ell,emm
INTEGER, INTENT(INOUT) :: clebsch


clebsch=binomial(ell+lambda-emm+mu,lambda+mu)*&
        &binomial(ell+lambda+emm-mu,lambda-mu)*binomial(2*ell+2*lambda+1,2*lambda)
clebsch=(-1)**(lambda+mu)*(clebsch)**0.5!sqrt



END SUBROUTINE clebsch_g

!The below subroutine computes the cartesian solid harmonic in
!the points x,y,z
!This routine is according to Molecular Electronic-Structure theory
! (Helgaker, Jørgensen, Olsen) page 209 Eq. 6.4.14.
SUBROUTINE cartesian_solid_harmonic(x,y,z,l,m,clmpfunc,ifconj)
IMPLICIT NONE
real(realk), intent(IN)::x,y,z
INTEGER, intent(IN) :: l,m
COMPLEX(complexk), intent(OUT) :: clmpfunc
LOGICAL,INTENT(IN) :: ifconj
!!!Local variables
real(realk) ::func
real(realk):: Nlm,Clmt
INTEGER :: t,sgn
!real(realk) :: pi=3.14159265

   func=0.0_realk
   if(ifconj) then
     do t=0,(l-abs(m))/2
        clmt=(-0.25)**t*binomial(l-t,abs(m)+t)*binomial(l,t)
        func=func+(x**2+y**2)**t*z**(l-2.*t-abs(m))*Clmt
     ENDDO
     write(9920,*) func
     
     Nlm=(2*l+1)/(4.*pi)*factorial(l+abs(m))*factorial(l-abs(m))
     !if(m .eq. 0) Nlm=Nlm/2. ! I do not know from where I got this expression
     Nlm=(-1)**((m+abs(m))/2)/(2**abs(m)*factorial(l))*sqrt(Nlm)
     sgn=0
     if(m .gt. 0) sgn=1
     if(m .lt. 0) sgn=-1
     if(m .eq. 0) sgn=0
     clmpfunc=func*Nlm*cmplx(x,sgn*y)**abs(m)
   else
     write(*,*) 'Not implemented yet'
     stop
   endif

END SUBROUTINE cartesian_solid_harmonic


SUBROUTINE cartesian_irreg_solid_harmonic(x,y,z,l,m,clmpfunc,ifconj)
IMPLICIT NONE
real(realk), intent(IN)::x,y,z
INTEGER, intent(IN) :: l,m
COMPLEX(complexk), intent(OUT) :: clmpfunc
LOGICAL,INTENT(IN) :: ifconj
!!!Local variables
real(realk) ::func
real(realk):: Nlm,Clmt
INTEGER :: t,sgn
COMPLEX(complexk) :: phase
!real(realk) :: pi=3.14159265

   func=0.0_realk
do t=0,(l-abs(m))/2
   clmt=(-0.25)**t*binomial(l-t,abs(m)+t)*binomial(l,t)
   func=func+(x**2+y**2)**t*z**(l-2.*t-abs(m))*Clmt/(x**2+y**2+z**2)**(l+0.5)
ENDDO

Nlm=(2*l+1)/(4.*pi)*factorial(l+abs(m))*factorial(l-abs(m))
!if(m .eq. 0) Nlm=Nlm/2. ! I do not know from where I got this expression
Nlm=(-1)**((m+abs(m))/2)/(2**abs(m)*factorial(l))*sqrt(Nlm)
sgn=0
if(m .gt. 0) sgn=1
if(m .lt. 0) sgn=-1
if(m .eq. 0) sgn=0
phase=cmplx(x,-sgn*y)
if(ifconj) phase=cmplx(x,sgn*y)

clmpfunc=func*Nlm*phase**abs(m)
!if( l.eq.30 .and. abs(m) .eq. 30)then
!  write(*,*) func,(x**2+y**2+z**2)**(l+0.5),clmpfunc
!  stop
!endif

END SUBROUTINE cartesian_irreg_solid_harmonic



REAL(realk) FUNCTION binomial(n,l)
IMPLICIT NONE
INTEGER, INTENT(IN) :: n,l

binomial=factorial(n)/factorial(n-l)/factorial(l)

END FUNCTION binomial

REAL(realk) FUNCTION factorial(n)
IMPLICIT NONE
INTEGER  :: n
INTEGER :: i
REAL(realk) :: summation
!real(realk) :: pi=3.14159265

IF(n .lt. 0) THEN
  WRITE(*,*) 'INPUT ARGUMENT NEGATIVE IN factorial, not possible, n= ',n
  STOP
ENDIF

IF(n .lt. 4) THEN
summation=0
DO i=1,n
 summation =summation+log(REAL(i))
ENDDO
ELSE
summation=n*log(real(n))-n+log(n*(1.+4.*n*(1.+2.*n)))/6+log(pi)/2.
ENDIF

factorial=exp(summation)
END FUNCTION


REAL(realk) FUNCTION Nlm(l,m)
IMPLICIT NONE
INTEGER, INTENT(IN) :: l,m
!!local variables
Integer :: phase
real(realk) :: frac
!real(realk) :: pi=3.14159265

phase=(-1)**((m+abs(m))/2)
frac = 1./(2**abs(m)*factorial(l))
Nlm=phase*frac*sqrt((2*l+1)/(4.*pi)*factorial(l+abs(m))*factorial(l-abs(m)))

END FUNCTION Nlm

REAL(realk) FUNCTION Clmt(l,m,t)
IMPLICIT NONE
INTEGER, INTENT(IN) :: l,m,t
!!local variables
real(realk) :: fract

fract=(-0.25)**t
Clmt=fract*binomial(l-t,abs(m)+t)*binomial(l,t)


END FUNCTION Clmt

END MODULE harmonics_pbc
