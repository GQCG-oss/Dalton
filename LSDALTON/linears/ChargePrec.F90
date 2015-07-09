
!***************************************************************************
!*                                                                         *
!*  MAIN ROUTINES TO COMPUTE GRADIENT,HESSIAN ON TRIAL VECTOR AND PRECOND  *
!*                                                                         *
!* Implemented for                                                         *
!* - Pipek measure using Lowdin or Mulliken charge and power m             *
!* (Charge localization using either Mulliken or Lowdin and power m)       *
!* - Pipek-Mezey with external power m, Lowdin charge only                 *
!* - Generalized Pipek-Mezey, external power m, Lowdin charge              *
!*                                                                         *
!* Lowdin orbitals may be optimized using orbspread scheme                 *
!*                                                                         *
!* NOTE: For Pipek measure -->   .CLL , m > 0                              *
!* NOTE: For generalized Pipek measure --> .CLL,  m < 0                    *
!* NOTE: For Pipek-Mezey --> .PML m=2, Lowdin pop.analysis                 *
!* NOTE: For Pipek-Mezey --> .PMM m=2, Mulliken pop.analysis               *
!*                                                                         *
!* Author: Ida-Marie Hoeyvik                                                *
!***************************************************************************

module ChargePrecMod
use matrix_module
use typedef
use matrix_util
use matrix_operations_aux, only: mat_dmul
use lowdin_module
use loc_types
CONTAINS
subroutine charge_precond(Xout,X,mu,inp)
use loc_types
use Pipek
use matrix_operations_aux
implicit none
type(Matrix), intent(inout) :: Xout
type(Matrix), intent(in) :: X
real(realk),  intent(in) :: mu
type(PMitem), intent(in) :: inp
real(realk),pointer   :: tmp(:), tmpP(:)
integer                   :: i, ne



ne = Xout%nrow*Xout%ncol
select case(matrix_type)
case(mtype_dense)
   call mat_zero(Xout)
   do i=1,ne
     if (dabs(inp%p%elms(i)- mu) > 1d-8 ) then
         Xout%elms(i)=X%elms(i)/(inp%P%elms(i)-mu)
     else
         Xout%elms(i)=X%elms(i)
     end if
   end do
 case default
    call mem_alloc(tmp,X%nrow*X%ncol)
    call mem_alloc(tmpP,inp%P%nrow*inp%P%ncol)
    call mat_to_full(X,1E0_realk,tmp)
    call mat_to_full(inp%P,1E0_realk,tmpP)

    do i=1,ne
    if (dabs(tmpP(i) - mu)> 1d-8) tmp(i) = tmp(i)/(tmpP(i) - mu)
    enddo
  
    call mat_set_from_full(tmp,1E0_realk,Xout)
    call mem_dealloc(tmp)
    call mem_dealloc(tmpP)
 end select


end subroutine charge_precond


subroutine PMMull_LinTra(Grad,HKappa,kappa,OrbLoc)
use loc_types
use Pipek
implicit none
type(PMitem)  :: OrbLoc
type(matrix)      :: Hkappa
type(matrix)      :: kappa
type(matrix)      :: Grad
integer             :: natoms,norb
type(matrix)        :: lA,lAd,lAk
real(realk),pointer :: factor(:)
real(realk),pointer :: temp(:,:)
real(realk),pointer :: lAki(:,:),lAi(:,:)
integer :: A,i
real(realk) :: DDOT,m
call mat_zero(Hkappa)
norb   = OrbLoc%norb
natoms = OrbLoc%natoms
m      = OrbLoc%m
OrbLoc%Qkappa = 0.0d0
call mem_alloc(factor,norb)
do A=1,natoms
   do i=1,norb
      factor(i)=(OrbLoc%Q(i,A))**(int(m)-1)
   enddo
   call mem_alloc(temp,OrbLoc%nbasA(A),norb)
   call mat_init(lA,OrbLoc%nbasA(A),norb)
   call mat_init(lAd,norb,OrbLoc%nbasA(A))
   call mat_init(lAk,OrbLoc%nbasA(A),norb)
   call mat_retrieve_block(OrbLoc%CMO,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
   call mat_set_from_full(temp,1.0d0,lA)
   !-m*diag(factor)*SC(muA,:)^T C(muA,:)kappa
   call mat_mul(lA,kappa,'n','n',1d0,0d0,lAk) !CMO(muA,:)*kappa
   call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
   ! ### MAKE (QaK)_ii needed for last terms
   call mem_alloc(lAki,OrbLoc%nbasA(A),norb)
   call mat_retrieve_block(lAk,lAki,OrbLoc%nbasA(A),norb,1,1)
   do i=1,norb
       OrbLoc%Qkappa(i,A)=0.5d0*DDOT(OrbLoc%nbasA(A),temp(:,i),1,lAki(:,i),1)
   end do
   ! ### end MAKE 
   call mat_set_from_full(temp,1.0d0,lA)
   call mat_dmul(factor,la,'T',1d0,0d0,laD)  !diag(factor)*SC(muA,:)^T
   call mat_mul(laD,lAk,'n','n',-m,1d0,HKappa) !HK=HK- diag(factor)*SC(muA,:)^T C(muA,:)kappa
   !mSC(muA,:)^T C(muA,:)kappa*diag(f) 
   call mat_dmul(factor,laK,'T',1d0,0d0,lAd)  !diag(f)*(C(muA,:)*kappa)^T
   call mat_mul(lA,lAd,'T','T',m,1d0,Hkappa) ! HK = HK+ (SC)^T*C*kappa*dia(Q)
   !-m*diag(factor)*C(muA,:)^T *SC*kappa
   call mat_mul(la,kappa,'n','n',1d0,0d0,lak)! lak= SC*kappa
   call mat_retrieve_block(OrbLoc%cmo,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
   ! ### MAKE (QaK)_ii needed for last terms
   call mat_retrieve_block(lAk,lAki,OrbLoc%nbasA(A),norb,1,1)
   do i=1,norb
       OrbLoc%Qkappa(i,A)=OrbLoc%Qkappa(i,A)+0.5d0*DDOT(OrbLoc%nbasA(A),temp(:,i),1,lAki(:,i),1)
   end do
   call mem_dealloc(laki)
   ! ### end MAKE 
   call mat_set_from_full(temp,1.0d0,lA)
   call mat_dmul(factor,la,'T',1d0,0d0,laD) !diag(factor)*C(muA,:)^T
   call mat_mul(lAd,lak,'n','n',-m,1d0,Hkappa)
   !m*C(muA,:)^T *SC*kappa*diag(f)
   call mat_dmul(factor,lak,'T',1d0,0d0,lAd)
   call mat_mul(lA,lAd,'T','T',m,1d0,Hkappa)
   do i=1,norb
      factor(i)=(OrbLoc%Q(i,A))**(int(m)-2)*OrbLoc%Qkappa(i,A)
   enddo
   call mat_dmul(factor,la,'T',1d0,0d0,lAd)
   call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
   call mat_set_from_full(temp,1.0d0,lA)
   call mat_mul(lAd,lA,'n','n',-2*m*(m-1d0),1d0,HKappa)
   call mat_dmul(factor,la,'T',1d0,0d0,lAd)
   call mat_retrieve_block(OrbLoc%CMO,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
   call mat_set_from_full(temp,1.0d0,lA)
   call mat_mul(lAd,lA,'n','n',-2*m*(m-1d0),1d0,HKappa)
   call mem_dealloc(temp)
   call mat_free(lA)
   call mat_free(lAd)
   call mat_free(lAk)
end do
call mem_dealloc(factor)
 ! Add gradient term
 call mat_mul(kappa,Grad,'n','n',0.5d0,1d0,Hkappa)

 call mat_init(lA,HKappa%ncol,HKappa%nrow)
 call mat_trans(Hkappa,lA)
 call mat_daxpy(-1d0,lA,Hkappa)
 call mat_free(lA)
 call mat_scal(-1d0,Hkappa)
! call mat_scal(-1d0,Hkappa) !HACK

end subroutine PMMull_LinTra






subroutine PMLowdin_LinTra(Grad,HKappa,kappa,OrbLoc)
use loc_types
use Pipek
implicit none
type(PMitem)  :: OrbLoc
type(matrix)      :: Hkappa
type(matrix)      :: kappa
type(matrix)      :: Grad
integer             :: natoms,norb
type(matrix)        :: lA,lAd,lAk
real(realk),pointer :: factor(:)
real(realk),pointer :: temp(:,:)
real(realk),pointer :: lAki(:,:),lAi(:,:)
integer :: A,i
real(realk) :: DDOT,m
call mat_zero(Hkappa)
norb   = OrbLoc%norb
natoms = OrbLoc%natoms
m      = OrbLoc%m
OrbLoc%Qkappa = 0.0d0
call mem_alloc(factor,norb)
do A=1,natoms
   do i=1,norb
      factor(i)=OrbLoc%Q(i,A)**(int(m)-1)
   enddo
   call mem_alloc(temp,OrbLoc%nbasA(A),norb)
   call mat_init(lA,OrbLoc%nbasA(A),norb)
   call mat_init(lAd,norb,OrbLoc%nbasA(A))
   call mat_init(lAk,OrbLoc%nbasA(A),norb)
   call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
   call mat_set_from_full(temp,1.0d0,lA)
   !-m*diag(factor)*SC(muA,:)^T C(muA,:)kappa
   call mat_mul(lA,kappa,'n','n',1d0,0d0,lAk) !SC(muA,:)*kappa
   call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
   ! ### MAKE (QaK)_ii needed for last terms
   call mem_alloc(lAki,OrbLoc%nbasA(A),norb)
   call mat_retrieve_block(lAk,lAki,OrbLoc%nbasA(A),norb,1,1)
   do i=1,norb
       OrbLoc%Qkappa(i,A)=DDOT(OrbLoc%nbasA(A),temp(:,i),1,lAki(:,i),1)
   end do
   call mem_dealloc(laki)
   ! ### end MAKE 
   call mat_set_from_full(temp,1.0d0,lA)
   call mat_dmul(factor,la,'T',1d0,0d0,laD)  !diag(factor)*SC(muA,:)^T
   call mat_mul(laD,lAk,'n','n',-2*m,1d0,HKappa) !HK=HK- diag(factor)*SC(muA,:)^T C(muA,:)kappa
   !mSC(muA,:)^T C(muA,:)kappa*diag(f) 
   call mat_dmul(factor,laK,'T',1d0,0d0,lAd)  !diag(f)*(C(muA,:)*kappa)^T
   call mat_mul(lA,lAd,'T','T',2*m,1d0,Hkappa) ! HK = HK+ (SC)^T*C*kappa*dia(Q)

   do i=1,norb
      factor(i)=OrbLoc%Q(i,A)**(int(m)-2)*OrbLoc%Qkappa(i,A)
   enddo
   call mat_dmul(factor,la,'T',1d0,0d0,lAd)
   call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
   call mat_set_from_full(temp,1.0d0,lA)
   call mat_mul(lAd,lA,'n','n',-4*m*(m-1d0),1d0,HKappa)
   call mem_dealloc(temp)
   call mat_free(lA)
   call mat_free(lAd)
   call mat_free(lAk)
end do
call mem_dealloc(factor)
 ! Add gradient term
 call mat_mul(kappa,Grad,'n','n',0.25d0,1d0,Hkappa)

 call mat_init(lA,HKappa%ncol,HKappa%nrow)
 call mat_trans(Hkappa,lA)
 call mat_daxpy(-1d0,lA,Hkappa)
 call mat_free(lA)
 call mat_scal(-1d0,Hkappa)

end subroutine PMLowdin_LinTra

end module ChargePrecMod


