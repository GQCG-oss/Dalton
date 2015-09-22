MODULE f12_module
use precision

TYPE GaussianGeminal
Logical :: is_set
Integer :: N
Real(realk) :: coeff(21)
Real(realk) :: exponent(21)
Real(realk) :: expProd(21)
END TYPE GaussianGeminal

INTERFACE set_GGem
  MODULE PROCEDURE set_GGem_1, set_GGem_2
END INTERFACE !set_GGem

private
public :: set_GGem, GaussianGeminal, init_GGem, free_GGem, stgfit, binominit
CONTAINS

SUBROUTINE init_GGem(GGem)
TYPE(GaussianGeminal) :: GGem
GGem%is_set = .FALSE.
END SUBROUTINE init_GGem

SUBROUTINE set_GGem_1(GGem,coeff,exponent,n)
TYPE(GaussianGeminal),intent(inout) :: GGem
integer,intent(in)                  :: n
REAL(realk),target,intent(in)       :: coeff(n)
REAL(realk),target,intent(in)       :: exponent(n)
GGem%is_set = .TRUE.
GGem%N = N
GGem%coeff(1:N) = coeff
GGem%exponent(1:N) = exponent
END SUBROUTINE set_GGem_1

SUBROUTINE set_GGem_2(GGem,coeff,exponent,expProd,n)
TYPE(GaussianGeminal),intent(inout) :: GGem
integer,intent(in)                  :: n
REAL(realk),target,intent(in)       :: coeff(n)
REAL(realk),target,intent(in)       :: exponent(n)
REAL(realk),target,intent(in)       :: expProd(n)
GGem%is_set = .TRUE.
GGem%N = N
GGem%coeff(1:N) = coeff
GGem%exponent(1:N) = exponent
GGem%expProd(1:N) = expProd
END SUBROUTINE set_GGem_2


SUBROUTINE free_GGem(GGem)
TYPE(GaussianGeminal) :: GGem
GGem%is_set = .FALSE.
END SUBROUTINE free_GGem


subroutine stgfit(slater,nlcg,lcge,lcgc)
use precision

!         slater (input): exponent of the Slater-type geminal (STG).
!         nlcg   (input): number of Gaussians used to represent the STG.
!         lcgc  (output): array with expansion coefficients.
!         lcge  (output): array with Gaussian exponents.
!    
implicit none

! input:
  integer, intent(in) :: nlcg

  real(realk), intent(in) :: slater

! output:
  real(realk), intent(out) :: lcge(nlcg), lcgc(nlcg)

! local:
  integer :: k

  real(realk) :: ggax(6,6), ggac(6,6)
    
! local data:
  data ggax / &
   & 0.6853E+00_realk,0.0000E+00_realk,0.0000E+00_realk,0.0000E+00_realk,0.0000E+00_realk,0.0000E+00_realk,& 
   & 0.4254E+00_realk,0.4520E+01_realk,0.0000E+00_realk,0.0000E+00_realk,0.0000E+00_realk,0.0000E+00_realk,& 
   & 0.3303E+00_realk,0.2321E+01_realk,0.1628E+02_realk,0.0000E+00_realk,0.0000E+00_realk,0.0000E+00_realk,& 
   & 0.2783E+00_realk,0.1591E+01_realk,0.7637E+01_realk,0.4574E+02_realk,0.0000E+00_realk,0.0000E+00_realk,& 
   & 0.2447E+00_realk,0.1225E+01_realk,0.4924E+01_realk,0.1988E+02_realk,0.1127E+03_realk,0.0000E+00_realk,& 
   & 0.2209E+00_realk,0.1004E+01_realk,0.3622E+01_realk,0.1216E+02_realk,0.4587E+02_realk,0.2544E+03_realk & 
   & /
  data ggac / & 
   & 0.7354E+00_realk,0.0000E+00_realk,0.0000E+00_realk,0.0000E+00_realk,0.0000E+00_realk,0.0000E+00_realk,&
   & 0.5640E+00_realk,0.3102E+00_realk,0.0000E+00_realk,0.0000E+00_realk,0.0000E+00_realk,0.0000E+00_realk,&
   & 0.4683E+00_realk,0.3087E+00_realk,0.1529E+00_realk,0.0000E+00_realk,0.0000E+00_realk,0.0000E+00_realk,&
   & 0.4025E+00_realk,0.3090E+00_realk,0.1570E+00_realk,0.8898E-01_realk,0.0000E+00_realk,0.0000E+00_realk,&
   & 0.3532E+00_realk,0.3072E+00_realk,0.1629E+00_realk,0.9321E-01_realk,0.5619E-01_realk,0.0000E+00_realk,&
   & 0.3144E+00_realk,0.3037E+00_realk,0.1681E+00_realk,0.9811E-01_realk,0.6024E-01_realk,0.3726E-01_realk &
   & /

! code:

  if ((nlcg.ge. 1).and.(nlcg.le. 6)) then
     do k=1,nlcg
        lcge(k)=ggax(k,nlcg)*slater**2
     end do
     do k=1,nlcg
        lcgc(k)=ggac(k,nlcg)/slater
     end do
  else
     call lsquit('This STG fit is not possible. Sorry.',-1)
  end if
end subroutine stgfit


subroutine binominit(binom,l)
use precision
! compute all binomial coefficients binom(k,n) up to binom(l,l)
! format: "n choose k" with k.le.n
! Wim Klopper@UiO/12.10.2010
integer, intent(in):: l
real(realk), intent(inout) :: binom(0:l,0:l)
integer n,n1,i
binom(0,0)=1E0_realk
if (l.eq. 0) return
do n=1,l
 n1=n+1
 binom(0,n)=1E0_realk
 do i=1,n/2
  binom(i,n)=binom(i-1,n)*(n1-i)/i
 enddo
 do i=n/2+1,n
  binom(i,n)=binom(n-i,n)
 enddo
enddo
return
end subroutine binominit


END MODULE f12_module
