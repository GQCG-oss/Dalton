!> @file
!> Contains the general McM driver 

!> \brief General McMurchie-Davidson Integral scheme
!> \author T. Kjaergaard
!> \date 2013 
MODULE IchorGaussianGeminalMod
use IchorPrecisionMod
use IchorParametersMod
use IchorMemory
integer :: nGGem
real(realk),allocatable :: GGemCoeff(:)
real(realk),allocatable :: GGemexponent(:)
real(realk),allocatable :: GGemprodexponent(:)
Real(realk),allocatable :: BinomArray(:,:)
logical :: GGemOperatorCalc
integer :: GGemOperatorSpec
private
public :: nGGem,GGemCoeff,GGemexponent,GGemprodexponent,&
     & BinomArray,GGemOperatorCalc,GGemOperatorSpec,&
     & set_GGem,free_GGem
 
!real(realk) :: ggax(6,6) !Change to 6
! local data:
!data ggax / &
!     & 0.6853E+00_realk,0.0000E+00_realk,0.0000E+00_realk,0.0000E+00_realk,0.0000E+00_realk,0.0000E+00_realk,& 
!     & 0.4254E+00_realk,0.4520E+01_realk,0.0000E+00_realk,0.0000E+00_realk,0.0000E+00_realk,0.0000E+00_realk,& 
!     & 0.3303E+00_realk,0.2321E+01_realk,0.1628E+02_realk,0.0000E+00_realk,0.0000E+00_realk,0.0000E+00_realk,& 
!     & 0.2783E+00_realk,0.1591E+01_realk,0.7637E+01_realk,0.4574E+02_realk,0.0000E+00_realk,0.0000E+00_realk,& 
!     & 0.2447E+00_realk,0.1225E+01_realk,0.4924E+01_realk,0.1988E+02_realk,0.1127E+03_realk,0.0000E+00_realk,& 
!     & 0.2209E+00_realk,0.1004E+01_realk,0.3622E+01_realk,0.1216E+02_realk,0.4587E+02_realk,0.2544E+03_realk & 
!     & /
!real(realk) :: ggac(6,6)    
!data ggac / & 
!     & 0.7354E+00_realk,0.0000E+00_realk,0.0000E+00_realk,0.0000E+00_realk,0.0000E+00_realk,0.0000E+00_realk,&
!     & 0.5640E+00_realk,0.3102E+00_realk,0.0000E+00_realk,0.0000E+00_realk,0.0000E+00_realk,0.0000E+00_realk,&
!     & 0.4683E+00_realk,0.3087E+00_realk,0.1529E+00_realk,0.0000E+00_realk,0.0000E+00_realk,0.0000E+00_realk,&
!     & 0.4025E+00_realk,0.3090E+00_realk,0.1570E+00_realk,0.8898E-01_realk,0.0000E+00_realk,0.0000E+00_realk,&
!     & 0.3532E+00_realk,0.3072E+00_realk,0.1629E+00_realk,0.9321E-01_realk,0.5619E-01_realk,0.0000E+00_realk,&
!     & 0.3144E+00_realk,0.3037E+00_realk,0.1681E+00_realk,0.9811E-01_realk,0.6024E-01_realk,0.3726E-01_realk &
!     & /
CONTAINS

  subroutine set_GGem(GGemOperatorIn,Jmax)
    implicit none
    integer,intent(in) :: GGemOperatorIn,Jmax
    !local variables
    real(realk) :: Coeff(6),expo(6)
!    real(realk),parameter :: slater = 1.0E0_realk

    real(realk) :: ggax(6),ggac(6)    
    integer :: I,J,IJ,K
    IF(GGemOperatorIn.EQ.CoulombOperator)THEN
       GGemOperatorCalc = .FALSE.
       GGemOperatorSpec = GGemOperatorIn       
       !do nothing this is a Coulomb operator no need for 
       !gaussian geminal
    ELSE
       ggax(1) = 0.2209E+00_realk
       ggax(2) = 0.1004E+01_realk
       ggax(3) = 0.3622E+01_realk
       ggax(4) = 0.1216E+02_realk
       ggax(5) = 0.4587E+02_realk
       ggax(6) = 0.2544E+03_realk 
       
       ggac(1) = 0.3144E+00_realk
       ggac(2) = 0.3037E+00_realk
       ggac(3) = 0.1681E+00_realk
       ggac(4) = 0.9811E-01_realk
       ggac(5) = 0.6024E-01_realk
       ggac(6) = 0.3726E-01_realk 
       
       GGemOperatorCalc = .TRUE.    
       GGemOperatorSpec = GGemOperatorIn
       IF(GGemOperatorSpec.EQ.GGemSqOperator.OR.GGemOperatorSpec.EQ.GGemGrdOperator)THEN
          nGGem = 21 !6*(6+1)/2
          allocate(GGemexponent(21))
          allocate(GGemprodexponent(21))
          allocate(GGemCoeff(21))
          call mem_ichor_alloc(GGemCoeff)
          call mem_ichor_alloc(GGemexponent)
          call mem_ichor_alloc(GGemprodexponent)
          do k=1,6
             expo(k)=ggax(k)!*slater**2
          end do
          do k=1,6
             Coeff(k)=ggac(k)!/slater
          end do
          IJ=0
          DO I=1,6
             DO J=1,I
                IJ = IJ + 1
                GGemCoeff(IJ) = 2E0_realk * coeff(I) * coeff(J)
                GGemprodexponent(IJ) = expo(I) * expo(J)
                GGemexponent(IJ) = expo(I) + expo(J)
             ENDDO
             GGemCoeff(IJ) = 0.5E0_realk*GGemCoeff(IJ)
          ENDDO
       ELSE 
          nGGem = 6
          allocate(GGemexponent(6))
          allocate(GGemCoeff(6))
          call mem_ichor_alloc(GGemCoeff)
          call mem_ichor_alloc(GGemexponent)
          do k=1,6
             GGemexponent(k)=ggax(k)!*slater**2
          end do
          do k=1,6
             GGemCoeff(k)=ggac(k)!/slater
          end do
       ENDIF
       IF(GGemOperatorSpec.EQ.GGemCouOperator)THEN
          allocate(BinomArray(0:JMAX,0:JMAX))
          call mem_ichor_alloc(BinomArray)
          call ichorbinominit(BinomArray,Jmax)  
       ENDIF
    ENDIF
  end subroutine set_GGem

  subroutine ichorbinominit(binom,l)
    ! compute all binomial coefficients binom(k,n) up to binom(l,l)
    integer, intent(in):: l
    real(realk), intent(inout) :: binom(0:l,0:l)
    integer n,n1,i
    binom(0,0)=1E0_realk
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
  end subroutine ichorbinominit

  subroutine free_GGem
    implicit none
    IF(GGemOperatorSpec.EQ.CoulombOperator)THEN
       GGemOperatorCalc = .FALSE.
       !do nothing this is a Coulomb operator no need for 
       !gaussian geminal
    ELSE
       GGemOperatorCalc = .FALSE.
       IF(GGemOperatorSpec.EQ.GGemSqOperator.OR.GGemOperatorSpec.EQ.GGemGrdOperator)THEN
          call mem_ichor_dealloc(GGemCoeff)
          call mem_ichor_dealloc(GGemexponent)
          call mem_ichor_dealloc(GGemprodexponent)
          deallocate(GGemexponent)
          deallocate(GGemprodexponent)
          deallocate(GGemCoeff)
       ELSE 
          call mem_ichor_dealloc(GGemCoeff)
          call mem_ichor_dealloc(GGemexponent)
          deallocate(GGemexponent)
          deallocate(GGemCoeff)
       ENDIF
       IF(GGemOperatorSpec.EQ.GGemCouOperator)THEN
          call mem_ichor_dealloc(BinomArray)
          deallocate(BinomArray)
       ENDIF
    ENDIF
  end subroutine free_GGem

end MODULE IchorGaussianGeminalMod
