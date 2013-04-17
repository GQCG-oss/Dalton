
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


module loc_types
use matrix_module
use typedef
use typedeftype

TYPE PMitem
! ******************************************
! *         Essential Settings              *
! ******************************************
! Exponent m
real(realk) :: m
! Use Pipek-Mezey scheme with Lowdin charge
logical :: PipekMezeyLowdin
! Use Pipek-Mezey scheme with Lowdin charge
logical :: PipekMezeyMull
! Use Pipek measure of locality w/Mulliken charge
logical :: ChargeLocMulliken
! Use Pipek measure of locality w/Lowdin charge
logical :: ChargeLocLowdin
! use linesearch in optimization
logical :: linesearch
! if using precond
logical :: precond


! ******************************************
! * Used for construction of LT and/or grad *
! ******************************************
! SU = S_sqrt for Lowdin orbitals
! SU = S for Mulliken based Pipek measure  
type(matrix) :: SU
!SU*Coeff
type(matrix) :: SC
!CMOS
type(matrix),pointer :: CMO
! number og atoms in molecule
integer :: natoms
! number of basis functs to use
integer :: nbas
! number of orbitals in mol (CMO%ncol)
integer :: norb
! d_i^m
real(realk),pointer :: di(:)
! d_i^(m+1)
real(realk),pointer :: d_m(:)
! d_i^(m+2)
real(realk),pointer :: d_m2(:)
! Pos(A); start pos. basisfunct on atom A  
integer,pointer :: Pos(:)
! nbasA(A); # of basisfunct on atom A
integer,pointer :: nbasA(:)
! Q(i,A) --> Q_{ii}^A
real(realk), pointer :: Q(:,:)
! Qkappa(i,A)-->[Q^A Kappa]_{ii}
real(realk),pointer :: Qkappa(:,:)


! ******************************************
! *       Used in localization scheme      *
! ******************************************
! current function value
real(realk) :: funcVal
! value indicating degree of localization
real(realk) :: loc_degree
! Matrix of diagonal Hessian elements
type(matrix),pointer :: P


! ******************************************
! *       Useful output etc.               *
! ******************************************
! Array that stores sorted lowdin charges
real(realk),pointer :: SortedLow(:,:)


END TYPE PMitem

end module loc_types

module Pipek
use matrix_module
use typedef
use matrix_util
use matrix_operations_aux, only: mat_dmul
use lowdin_module
use loc_types
CONTAINS


! MAIN ROUTINES FOR GRADIENT AND LINEAR TRANSFORMATIONS
!==========================================================

subroutine Gradient_ChargeLoc(Grad,OrbLoc)
implicit none
type(PMitem) :: OrbLoc
type(matrix) :: Grad
!!##### TESTING ########
!type(matrix) :: kappa
!type(matrix) :: Hkappa,diaH
!real(realk) :: mat(1,1)
!integer :: p,q
!p=2
!q=3
!call mat_init(kappa,OrbLoc%CMO%ncol,OrbLoc%CMO%ncol)
!call mat_init(Hkappa,OrbLoc%CMO%ncol,OrbLoc%CMO%ncol)
!call mat_init(diaH,OrbLoc%CMO%ncol,OrbLoc%CMO%ncol)
!call mat_zero(kappa)
!mat(1,1)=1d0
!call mat_create_block(kappa,mat,1,1,q,p)
!mat(1,1)=-1d0
!call mat_create_block(kappa,mat,1,1,p,q)
!if (OrbLoc%PipekMezeyLowdin) then
!  call ComputeGrad_PipekMezey(OrbLoc,Grad)
!  call PMLowdin_precond(OrbLoc,diaH)
!  call PMLowdin_LinTra(Grad,HKappa,kappa,OrbLoc) 
!elseif (OrbLoc%PipekMezeyMull) then
!  call ComputeGrad_PipekMezey(OrbLoc,Grad)
!  call PMMull_LinTra(Grad,HKappa,kappa,OrbLoc)
!  call PMMull_precond(OrbLoc,diaH)
!elseif (OrbLoc%ChargeLocLowdin .or. OrbLoc%ChargeLocMulliken) then
!  call ComputeGrad_ChargeLoc(OrbLoc,Grad) 
!  call CL_precond(OrbLoc,diaH)
!  call CLLinearTrans(Grad,HKappa,kappa,OrbLoc)
!endif
!!call PMMull_precond(OrbLoc,diah)
!!call CL_precond(OrbLoc,diaH)
!print*,"**** GRADIENT ****"
!call mat_print(Grad,1,6,1,6,6)
!print*, "**** DIAGONAL HESS ****"
!call mat_print(diaH,1,6,1,6,6)
!print*, "**** LINEAR TRANS  ****"
!call mat_print(Hkappa,1,6,1,6,6)
!STOP 'TESTING'
!
!
!!!##### END TESTING ######

if (OrbLoc%ChargeLocLowdin .or. OrbLoc%ChargeLocMulliken) then
   call ComputeGrad_ChargeLoc(OrbLoc,Grad)
elseif (OrbLoc%PipekMezeyLowdin .or. OrbLoc%PipekMezeyMull) then
   call ComputeGrad_PipekMezey(OrbLoc,Grad)
else
  STOP  'Something wrong in ChargeLoc.f90. Neither option is chosen'
end if


end subroutine Gradient_ChargeLoc

!> \brief Calls correct routine for computing diagonal hessian elements
!> \param P Matrix with diagonal elements (output)
subroutine Precond_ChargeLoc(P,OrbLoc)
implicit none
type(PMitem) :: OrbLoc
type(matrix) :: P


if (OrbLoc%ChargeLocLowdin .or. OrbLoc%ChargeLocMulliken) then
   call CL_precond(OrbLoc,P)
elseif (OrbLoc%PipekMezeyLowdin) then
   call PMLowdin_precond(OrbLoc,P)
elseif (OrbLoc%PipekMezeyMull) then
   call PMMull_precond(OrbLoc,P)
endif


end subroutine Precond_ChargeLoc

!> \brief Compute gradient for charge localization (both Mull. and Lowdin)
!> \param Grad gradient (output)
subroutine ComputeGrad_ChargeLoc(OrbLoc,Grad)
implicit none
type(PMitem)      :: OrbLoc
type(matrix)            :: Grad
type(matrix)            :: lA,rA,lAd,rAd
real(realk),pointer :: temp(:,:) 
real(realk),pointer :: factor(:)
integer :: i,A,norb

call mat_zero(Grad)
norb = OrbLoc%norb

call mem_alloc(factor,norb) 
if (OrbLoc%ChargeLocMulliken) then
   ComputeGrad: do A=1,OrbLoc%natoms
      call mat_init(lA,OrbLoc%nbasA(A),norb)
      call mat_init(rA,OrbLoc%nbasA(A),norb)
      call mat_init(lAd,norb,OrbLoc%nbasA(A))
      call mat_init(rAd,norb,OrbLoc%nbasA(A))
      do i=1,norb;factor(i)=OrbLoc%Q(i,A)*OrbLoc%d_m(i);end do
      !Construct lA and rA
      call mem_alloc(temp,OrbLoc%nbasA(A),norb)
      call mat_retrieve_block(OrbLoc%CMO,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
      call mat_set_from_full(temp,1.0d0,lA)
      call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
      call mat_set_from_full(temp,1.0d0,rA)
      call mem_dealloc(temp)
      ! Make lAd=diag(f)*lA^T
      call mat_dmul(factor(:),lA,'T',1.0d0,0.0d0,lAd)
      call mat_dmul(factor(:),rA,'T',1.0d0,0.0d0,rAd)
      
      call mat_mul(rAd,lA,'n','n',-2.0d0*OrbLoc%m,1.0d0,Grad)
      call mat_mul(lAd,ra,'n','n',-2.0d0*OrbLoc%m,1.0d0,Grad)
      call mat_mul(rA,lAd,'T','T',2.0d0*OrbLoc%m,1.0d0,Grad)
      call mat_mul(lA,rAd,'T','T',2.0d0*OrbLoc%m,1.0d0,Grad)
      call mat_free(lA);call mat_free(rA)
      call mat_free(lAd);call mat_free(rAd)
   end do ComputeGrad
end if
if (Orbloc%ChargeLocLowdin) then
   ComputeGradL: do A=1,OrbLoc%natoms
      call mat_init(lA,OrbLoc%nbasA(A),norb)
      call mat_init(lAd,norb,OrbLoc%nbasA(A))
      factor=OrbLoc%Q(:,A)*OrbLoc%d_m
      !Construct lA and rA
      call mem_alloc(temp,OrbLoc%nbasA(A),norb)
      call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
      call mat_set_from_full(temp,1.0d0,lA)
      call mem_dealloc(temp)
      ! Make lAd=diag(f)*lA^T
      call mat_dmul(factor,lA,'T',1.0d0,0.0d0,lAd)
      call mat_mul(lAd,lA,'n','n',-4*OrbLoc%m,1.0d0,Grad)
      call mat_free(lA)
      call mat_free(lAd)
   end do ComputeGradL
   call mat_init(la,Grad%ncol,Grad%nrow)
   call mat_trans(Grad,la)
   call mat_daxpy(-1d0,la,Grad)
   call mat_free(la)
end if
call mem_dealloc(factor)

call mat_scal(-1.0d0,Grad)
if (OrbLoc%m < 0) then 
    call mat_scal(-1.0d0,Grad)
end if

end subroutine ComputeGrad_ChargeLoc

!> \brief COmpute gradient for PM loc (both Mull and Lowdin)
!> \param Grad  gradient (output)
subroutine ComputeGrad_PipekMezey(OrbLoc,Grad)
implicit none
type(PMitem)      :: OrbLoc
type(matrix)            :: Grad
type(matrix)            :: lA,lAd
real(realk),pointer :: temp(:,:) 
real(realk),pointer :: factor(:)
real(realk) :: m
integer     :: i,A,norb
type(matrix) :: id

norb = OrbLoc%norb
m= OrbLoc%m
call mat_zero(Grad)

call mem_alloc(factor,norb) 
if (OrbLoc%PipekMezeyLowdin) then
  ComputeGradLowdin: do A=1,OrbLoc%natoms
     call mat_init(lA,OrbLoc%nbasA(A),norb)
     call mat_init(lAd,norb,OrbLoc%nbasA(A))
     factor=OrbLoc%Q(:,A)**(int(m)-1)
     !Construct lA and rA
     call mem_alloc(temp,OrbLoc%nbasA(A),norb)
     call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
     call mat_set_from_full(temp,1.0d0,lA)
     call mem_dealloc(temp)
     ! Make lAd=diag(f)*lA^T
     call mat_dmul(factor(:),lA,'T',1.0d0,0.0d0,lAd)
     call mat_mul(lAd,lA,'n','n',-2.0d0*m,1.0d0,Grad)
     call mat_free(lA)
     call mat_free(lAd)
  end do ComputeGradLowdin
     call mat_init(lA,Grad%ncol,Grad%nrow)
     call mat_trans(Grad,lA)
     call mat_daxpy(-1d0,lA,Grad)
     call mat_free(lA)
elseif (OrbLoc%PipekMezeyMull) then
  ComputeGradMull: do A=1,OrbLoc%natoms
     call mat_init(lA,OrbLoc%nbasA(A),norb)
     call mat_init(lAd,norb,OrbLoc%nbasA(A))
     factor=(OrbLoc%Q(:,A))**(int(m)-1)
     call mem_alloc(temp,OrbLoc%nbasA(A),norb)
     call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
     call mat_set_from_full(temp,1.0d0,lA)
     call mat_dmul(factor,lA,'T',1.0d0,0.0d0,lAd)
     call mat_retrieve_block(OrbLoc%CMO,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
     call mat_set_from_full(temp,1.0d0,lA)
     call mat_mul(laD,la,'n','n',m,1d0,Grad) !factor half from mull incl.
     !call mat_mul(la,laD,'T','T',-m,1d0,Grad) !factor half from mull incl.
     call mat_dmul(factor,lA,'T',1.0d0,0.0d0,lAd)
     call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
     call mat_set_from_full(temp,1.0d0,lA)
     call mat_mul(laD,la,'n','n',m,1d0,Grad) !factor half from mull incl.
     !call mat_mul(la,lad,'T','T',-m,1d0,Grad) !factor half from mull incl.
     call mem_dealloc(temp)
     call mat_free(la)
     call mat_free(laD)
  end do ComputeGradMull
  call mat_init(lA,Grad%ncol,Grad%nrow)
  call mat_trans(Grad,lA)
  call mat_daxpy(-1d0,lA,Grad)
  call mat_free(lA)
  call mat_scal(-1d0,Grad)
end if
call mem_dealloc(factor)

end subroutine ComputeGrad_PipekMezey

!> \brief Compute preconditioner for charge localization (both Mull. and Lowdin)
!> \param diaH matrix with diagonal elements
subroutine CL_precond(OrbLoc,diaH)
use matrix_operations_aux
implicit none
type(PMitem) :: OrbLoc
type(matrix) :: diaH
real(realk),pointer :: x(:)
integer             :: natoms,norb
type(matrix)        :: lA,rA,Qtemp,QB
integer :: A,i
real(realk) :: m
real(realk),pointer :: temp(:,:)

norb   = OrbLoc%norb
natoms = OrbLoc%natoms
m      = OrbLoc%m
call mat_zero(diaH)
call mat_init(QB,diaH%nrow,diaH%ncol)
call mat_zero(QB)
do A=1,natoms
   call mat_init(lA,OrbLoc%nbasA(A),norb)
   call mat_init(rA,norb,OrbLoc%nbasA(A))
   call mem_alloc(temp,OrbLoc%nbasA(A),norb)
   if (OrbLoc%ChargeLocMulliken) then
      call mat_retrieve_block(OrbLoc%CMO,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
   elseif (OrbLoc%ChargeLocLowdin) then
      call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
   endif
   call mat_set_from_full(temp,1.0d0,lA)
   call mat_dmul(OrbLoc%Q(:,A),lA,'T',1d0,0d0,rA)
   call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
   call mat_set_from_full(temp,1.0d0,lA)
   if (OrbLoc%ChargeLocMulliken) then
      call mat_mul(ra,la,'n','n',0.5d0,1d0,QB)
      call mat_dmul(OrbLoc%Q(:,A),lA,'T',1d0,0d0,rA)
      call mat_retrieve_block(OrbLoc%CMO,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
      call mat_set_from_full(temp,1.0d0,lA)
      call mat_mul(ra,la,'n','n',0.5d0,1d0,QB)
   elseif (OrbLoc%ChargeLocLowdin)then
      call mat_mul(ra,la,'n','n',1d0,1d0,QB)
   end if
   call mat_free(la)
   call mat_free(ra)
   call mem_dealloc(temp)
enddo

call mem_alloc(x,norb)
x=1d0
do A=1,natoms
  call mat_dger(4*m,x,OrbLoc%d_m*OrbLoc%Q(:,A)*OrbLoc%Q(:,A),diaH)
  call mat_dger(-4*m,OrbLoc%Q(:,A)*OrbLoc%d_m,OrbLoc%Q(:,A),diaH)

   call mat_init(lA,OrbLoc%nbasA(A),norb)
   call mat_init(rA,OrbLoc%nbasA(A),norb)
   call mem_alloc(temp,OrbLoc%nbasA(A),norb)
   if (OrbLoc%ChargeLocMulliken) then
       call mat_retrieve_block(OrbLoc%CMO,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
   elseif (OrbLoc%ChargeLocLowdin) then
       call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
   end if
   call mat_set_from_full(temp,1.0d0,lA)
   call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
   call mat_set_from_full(temp,1.0d0,rA)
   call mem_dealloc(temp)
   call mat_init(Qtemp,norb,norb)
   if (OrbLoc%ChargeLocMulliken) then
       call mat_mul(la,ra,'T','n',0.5d0,0d0,Qtemp)
       call mat_mul(ra,la,'T','n',0.5d0,1d0,Qtemp)
   elseif (OrbLoc%ChargeLocLowdin) then
       call mat_mul(ra,la,'T','n',1.0d0,0d0,Qtemp)
   end if
   call mat_dhmul(OrbLoc%d_m,Qtemp,Qtemp,'n','n',-8*m,1d0,diaH)
   call mat_free(la)
   call mat_free(ra)
   call mat_dhmul(OrbLoc%d_m2*OrbLoc%Q(:,A),Qtemp,QB,'n','n',16*m*(m+1),1d0,diaH)
   call mat_free(Qtemp)
enddo
call mem_dealloc(x)
call mat_init(Qtemp,diaH%ncol,diaH%nrow)
call mat_trans(diaH,Qtemp)
call mat_daxpy(1d0,Qtemp,diaH)
call mat_free(Qtemp)
call mat_scal_dia(0d0,diaH)

if(OrbLoc%m< 0) call mat_scal(-1d0,diaH)


end subroutine CL_precond


!> \brief Compute precond for PM using Lowdin pop.analysis
!> \param diaH diagonal elements (output)
subroutine PMLowdin_precond(OrbLoc,diaH)
use matrix_operations_aux
implicit none
type(PMitem) :: Orbloc
type(matrix)     :: diaH
real(realk),pointer :: x(:)
integer             :: natoms,norb
type(matrix)        :: lA,rA,Qtemp
integer :: A,i
real(realk) :: m
real(realk),pointer :: temp(:,:)
norb   = OrbLoc%norb
natoms = OrbLoc%natoms
m      = OrbLoc%m
call mat_zero(diaH)

call mem_alloc(x,norb)
x=1d0
do A=1,natoms
  call mat_dger(2*m,OrbLoc%Q(:,A),OrbLoc%Q(:,A)**(int(m)-1),diaH)
  call mat_dger(-2*m,x,OrbLoc%Q(:,A)*(OrbLoc%Q(:,A)**(int(m)-1)),diaH)
  
   call mat_init(lA,OrbLoc%nbasA(A),norb)
   call mat_init(rA,OrbLoc%nbasA(A),norb)
   call mem_alloc(temp,OrbLoc%nbasA(A),norb)
   call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
   call mat_set_from_full(temp,1.0d0,lA)
   call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
   call mat_set_from_full(temp,1.0d0,rA)
   call mem_dealloc(temp)
   call mat_init(Qtemp,norb,norb)
   call mat_mul(la,ra,'T','n',1d0,0d0,Qtemp)
   call mat_dhmul(OrbLoc%Q(:,A)**(int(m)-2),Qtemp,Qtemp,'n','n',4*m*(m-1),1d0,diaH)
   call mat_free(Qtemp)
   call mat_free(la)
   call mat_free(ra)
end do
call mem_dealloc(x)

call mat_init(Qtemp,diaH%ncol,diaH%nrow)
call mat_trans(diaH,Qtemp)
call mat_daxpy(1d0,Qtemp,diaH)
call mat_free(Qtemp)
call mat_scal(-1d0,diaH)
call mat_scal_dia(0d0,diaH)
end subroutine PMLowdin_precond

!> \brief Compute preconditioner for PM using Mull. Pop.analysis
!> \param diaH matrix with diagonal elements
subroutine PMMull_precond(OrbLoc,diaH)
use matrix_operations_aux
implicit none
type(PMitem) :: Orbloc
type(matrix)     :: diaH
real(realk),pointer :: x(:)
integer             :: natoms,norb
type(matrix)        :: lA,rA,Qtemp
integer :: A,i
real(realk) :: m
real(realk),pointer :: temp(:,:)

norb   = OrbLoc%norb
natoms = OrbLoc%natoms
m      = OrbLoc%m
call mat_zero(diaH)

call mem_alloc(x,norb)
x=1d0
do A=1,natoms
  call mat_dger(2*m,OrbLoc%Q(:,A),OrbLoc%Q(:,A)**(int(m)-1),diaH)
  call mat_dger(-2*m,x,OrbLoc%Q(:,A)*(OrbLoc%Q(:,A)**(int(m)-1)),diaH)
  
   call mat_init(lA,OrbLoc%nbasA(A),norb)
   call mat_init(rA,OrbLoc%nbasA(A),norb)
   call mem_alloc(temp,OrbLoc%nbasA(A),norb)
   call mat_retrieve_block(OrbLoc%CMO,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
   call mat_set_from_full(temp,1.0d0,lA)
   call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
   call mat_set_from_full(temp,1.0d0,rA)
   call mem_dealloc(temp)
   call mat_init(Qtemp,norb,norb)
   call mat_mul(la,ra,'T','n',0.5d0,0d0,Qtemp)
   call mat_mul(ra,la,'T','n',0.5d0,1d0,Qtemp)
   call mat_dhmul(OrbLoc%Q(:,A)**(int(m)-2),Qtemp,Qtemp,'n','n',4*m*(m-1),1d0,diaH)
   call mat_free(Qtemp)
   call mat_free(la)
   call mat_free(ra)
end do
call mem_dealloc(x)
call mat_init(Qtemp,diaH%ncol,diaH%nrow)
call mat_trans(diaH,Qtemp)
call mat_daxpy(1d0,Qtemp,diaH)
call mat_free(Qtemp)
call mat_scal(-1d0,diaH)
call mat_scal_dia(0d0,diaH)


end subroutine PMMull_precond

! MISCELLANEOUS ROUTINES
!===========================



!> \brief Computes information about #bas.functions per atom and stores positions
!> \param OrbLoc%pos(A) Position from where the basis functions centered on A are placed
!> \param OrbLoc%nbasA How many basis functions are assigned to A
subroutine get_information(ls,OrbLoc)
implicit none
type(lsitem) :: ls
type(PMitem) :: OrbLoc
integer :: A

do A=1,OrbLoc%natoms
   OrbLoc%nbasA(A) = ls%setting%Molecule(1)%p%Atom(A)%nContOrbREG
   if (A==1) OrbLoc%pos(A)=1
   if (A>1) OrbLoc%pos(A)=OrbLoc%pos(A-1)+OrbLoc%nbasA(A-1)
end do

end subroutine get_information

!> \brief Computes d_m and population matrix diagonal 
!> \param Orbloc%Q(i,A) Q_ii^A (orb i, atomic center A) 
subroutine get_dm_and_Qii(CMO,OrbLoc)
implicit none
type(PMitem)      :: OrbLoc
type(matrix)            :: CMO
real(realk),pointer :: lA(:,:),rA(:,:)
integer     :: A,i,norb,natoms
real(realk) :: ddot,m,scalar

natoms=OrbLoc%natoms
norb=OrbLoc%norb
m=OrbLoc%m
Lowdin: if (OrbLoc%PipekMezeyLowdin) then
      A1:do A=1,natoms
         call mem_alloc(lA,OrbLoc%nbasA(A),norb)
         call mat_retrieve_block(OrbLoc%SC,lA,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
         do i=1,norb
                OrbLoc%Q(i,A)=DDOT(OrbLoc%nbasA(A),lA(:,i),1,lA(:,i),1)
         end do
         call mem_dealloc(lA)
      end do A1
endif Lowdin
Mulliken: if (OrbLoc%PipekMezeyMull) then
 do A=1,natoms
    call mem_alloc(lA,OrbLoc%nbasA(A),norb)
    call mem_alloc(rA,OrbLoc%nbasA(A),norb)
    call mat_retrieve_block(OrbLoc%SC,lA,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
    call mat_retrieve_block(OrbLoc%CMO,rA,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
    do i=1,norb
        scalar=DDOT(OrbLoc%nbasA(A),rA(:,i),1,lA(:,i),1)
        OrbLoc%Q(i,A)= scalar
    end do
    call mem_dealloc(lA)
    call mem_dealloc(rA)
 enddo
end if Mulliken
ChargeLocMulliken: if (OrbLoc%ChargeLocMulliken) then
OrbLoc%d_m2=0.0d0
OrbLoc%d_m=0.0d0
      A2:do A=1,natoms
         call mem_alloc(lA,OrbLoc%nbasA(A),norb)
	 call mem_alloc(rA,OrbLoc%nbasA(A),norb)
         call mat_retrieve_block(CMO,lA,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
         call mat_retrieve_block(OrbLoc%SC,rA,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
         do i=1,norb
                OrbLoc%Q(i,A)=DDOT(OrbLoc%nbasA(A),lA(:,i),1,rA(:,i),1)
                OrbLoc%d_m(i)=OrbLoc%d_m(i)+OrbLoc%Q(i,A)*OrbLoc%Q(i,A)
                if (A==natoms) then
      	          OrbLoc%di(i)= OrbLoc%d_m(i)
      	          OrbLoc%d_m2(i)= OrbLoc%d_m(i)
      	          OrbLoc%d_m(i) = OrbLoc%d_m(i)**(-(m+1))
      	          OrbLoc%d_m2(i)= OrbLoc%d_m2(i)**(-(m+2))
	          OrbLoc%di(i) = OrbLoc%di(i)**(-m)
		  if (OrbLoc%m<0) OrbLoc%di(i)=-OrbLoc%di(i)
                end if
         end do
         call mem_dealloc(lA)
	 call mem_dealloc(rA)
      end do A2
end if ChargeLocMulliken 

ChargeLocLowdin:if (OrbLoc%ChargeLocLowdin) then
! Use pipek measure with Lowdin charges
OrbLoc%d_m2=0.0d0
OrbLoc%d_m=0.0d0
      do A=1,natoms
         call mem_alloc(lA,OrbLoc%nbasA(A),norb)
         call mat_retrieve_block(OrbLoc%SC,lA,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
         do i=1,norb
                OrbLoc%Q(i,A)=DDOT(OrbLoc%nbasA(A),lA(:,i),1,lA(:,i),1)
                OrbLoc%d_m(i)=OrbLoc%d_m(i)+OrbLoc%Q(i,A)*OrbLoc%Q(i,A)
                if (A==natoms) then
                  OrbLoc%d_m2(i)= OrbLoc%d_m(i)
                  OrbLoc%di(i)= OrbLoc%d_m(i)
                  OrbLoc%d_m(i) = OrbLoc%d_m(i)**(-(m+1))
                  OrbLoc%d_m2(i)= OrbLoc%d_m2(i)**(-(m+2))
                  OrbLoc%di(i)= OrbLoc%di(i)**(-m)
		  if (OrbLoc%m<0) OrbLoc%di(i)=-OrbLoc%di(i)
                end if
         end do
         call mem_dealloc(lA)
      end do 
end if ChargeLocLowdin

end subroutine get_dm_and_Qii

subroutine compute_loc_degree(OrbLoc)
implicit none
type(PMitem) :: OrbLoc
real(realk) :: max_i(OrbLoc%norb)
real(realk) :: func(OrbLoc%natoms)
real(realk) :: temp(Orbloc%norb,Orbloc%natoms)
integer :: i,A,indx(1)

temp = OrbLoc%Q(:,:)

   do i=1,Orbloc%norb
      do A=1,Orbloc%natoms
          Orbloc%SortedLow(A,i)= maxval(temp(i,:)) 
          indx = maxloc(temp(i,1:OrbLoc%natoms))
          temp(i,indx(1))=0.0d0
      end do
   end do

if (OrbLoc%PipekMezeyLowdin .or. OrbLoc%PipekMezeyMull) then
   do i=1,OrbLoc%norb
      func=OrbLoc%Q(i,:)
      max_i(i)= maxval(func) 
   end do
   OrbLoc%loc_degree = minval(max_i)
else
   OrbLoc%loc_degree = maxval(OrbLoc%di)
end if

end subroutine

subroutine update_Orbloc(Orbloc,CMO,ls)
implicit none
type(PMitem) :: OrbLoc
type(lsitem)     :: ls
type(matrix)     :: CMO
integer :: i

call mat_mul(OrbLoc%SU,CMO,'T','n',1.0d0,0.0d0,OrbLoc%SC)

call get_information(ls,OrbLoc)
call get_dm_and_Qii(CMO,OrbLoc)
call compute_fVal(OrbLoc)
call compute_loc_degree(OrbLoc)


end subroutine update_OrbLoc

subroutine initialize_OrbLoc(CMO,OrbLoc,ls,m)
implicit none
type(PMitem) :: OrbLoc
type(lsitem)     :: ls
type(matrix)     :: CMO
real(realk)      :: m

OrbLoc%natoms = ls%setting%MOLECULE(1)%p%nAtoms
OrbLoc%norb = CMO%ncol
OrbLoc%nbas = CMO%nrow
Orbloc%m = m

call mem_alloc(OrbLoc%Pos,OrbLoc%natoms)
call mem_alloc(OrbLoc%nbasA,OrbLoc%natoms)
call mem_alloc(OrbLoc%Q,OrbLoc%norb,OrbLoc%natoms)
call mem_alloc(OrbLoc%Qkappa,OrbLoc%norb,OrbLoc%natoms)
call mem_alloc(Orbloc%SortedLow,OrbLoc%natoms,OrbLoc%norb)
if (.not. (OrbLoc%PipekMezeyLowdin.or.OrbLoc%PipekMezeyMull )) then
  call mem_alloc(OrbLoc%d_m,OrbLoc%norb)
  call mem_alloc(OrbLoc%d_m2,OrbLoc%norb)
  call mem_alloc(OrbLoc%di,OrbLoc%norb)
   OrbLoc%d_m=0.0d0;OrbLoc%d_m2=0.0d0;OrbLoc%di=0.0d0
end if

OrbLoc%Q=0.0d0
OrbLoc%Qkappa =0.0d0

call mat_init(OrbLoc%SC,OrbLoc%nbas,OrbLoc%norb)

end subroutine initialize_OrbLoc


subroutine FreeOrbLoc(OrbLoc)
implicit none
type(PMitem) :: OrbLoc

call mem_dealloc(OrbLoc%Pos)
call mem_dealloc(OrbLoc%nbasA)
call mem_dealloc(OrbLoc%Q)
call mem_dealloc(OrbLoc%Qkappa)
call mem_dealloc(OrbLoc%SortedLow)

if (.not. (OrbLoc%PipekMezeyLowdin .or. OrbLoc%PipekMezeyMull)) then
   call mem_dealloc(OrbLoc%d_m)
   call mem_dealloc(OrbLoc%d_m2)
   call mem_dealloc(OrbLoc%di)
end if
call mat_free(OrbLoc%SC)

end subroutine


subroutine compute_fVal(OrbLoc)
implicit none
type(PMitem) :: OrbLoc
real(realk) :: sum_Q
integer :: A

sum_Q = 0.0d0
if (OrbLoc%PipekMezeyLowdin .or. OrbLoc%PipekMezeyMull) then
  do A =1,Orbloc%natoms
    sum_Q = sum_Q - sum(OrbLoc%Q(:,A)**int(OrbLoc%m)) 
  end do
else
      sum_Q=sum(OrbLoc%di)
endif

OrbLoc%funcVal = sum_Q

end subroutine compute_fVal


!  SEPARATE TERMS FOR GRADIENT, LINTRANS AND PRECOND
!==========================================================


!Additional result:  Qkappa
subroutine CL_add_terms_1_4L(OrbLoc,kappa,Hkappa) 
implicit none
type(PMitem)  :: OrbLoc
integer       :: natoms,norb,m
type(matrix)  :: kappa,Hkappa !Hkappa: result
type(matrix)  :: lA,lAd,lAk
real(realk)   :: factor(OrbLoc%norb)
real(realk),pointer :: temp(:,:)
real(realk),pointer :: lAki(:,:),lAi(:,:)
integer :: A,i
real(realk) :: DDOT

norb   = OrbLoc%norb
natoms = OrbLoc%natoms
m      = OrbLoc%m
OrbLoc%Qkappa = 0.0d0

do A=1,natoms
   factor=OrbLoc%Q(:,A)*OrbLoc%d_m
   call mem_alloc(temp,OrbLoc%nbasA(A),norb)
   call mat_init(lA,OrbLoc%nbasA(A),norb)
   call mat_init(lAd,norb,OrbLoc%nbasA(A))
   call mat_init(lAk,OrbLoc%nbasA(A),norb)
   call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
   call mat_set_from_full(temp,1.0d0,lA)
   call mem_dealloc(temp)
   ! Make lAd=diag(f)*lA^T
   call mat_dmul(factor(:),lA,'T',1.0d0,0.0d0,lAd)
   ! Make lAk= lA*kappa 
   call mat_mul(lA,kappa,'n','n',1.0d0,0.0d0,lAk)
   ! Term 1
   call mat_mul(lAd,lAk,'n','n',4.0d0*m,1.0d0,Hkappa)
   ! Term 2
   call mat_mul(lAk,lAd,'T','T',-4.0d0*m,1.0d0,Hkappa)

   !Make lAkd= diag(dQ)*lAk^T, use lAd as the matrix lAkd
   call mat_dmul(factor(:),lAk,'T',1.0d0,0.0d0,lAd)

   ! Term 3
   call mat_mul(lAd,lA,'n','n',4.0d0*m,1.0d0,Hkappa)
   ! Term 4
   call mat_mul(lA,lAd,'T','T',-4.0d0*m,1.0d0,Hkappa)

   !Use lAk to make diag(Qkappa)
   call mem_alloc(lAki,OrbLoc%nbasA(A),norb)
   call mem_alloc(lAi,OrbLoc%nbasA(A),norb)
   call mat_retrieve_block(lAk,lAki,OrbLoc%nbasA(A),norb,1,1)
   call mat_retrieve_block(lA,lAi,OrbLoc%nbasA(A),norb,1,1)
   do i=1,norb
      OrbLoc%Qkappa(i,A)=DDOT(OrbLoc%nbasA(A),lAi(:,i),1,lAki(:,i),1)
   end do
 
   call mem_dealloc(lAki)
   call mem_dealloc(lAi)
   call mat_free(lA)
   call mat_free(lAd)
   call mat_free(lAk)
end do

end subroutine CL_add_terms_1_4L

!Additional result:  Qkappa
subroutine CL_add_terms_1_4M(OrbLoc,kappa,Hkappa) 
implicit none
type(PMitem)  :: OrbLoc
integer       :: natoms,norb,m
type(matrix)  :: kappa,Hkappa !Hkappa: result
type(matrix)  :: lA,rA,lAd,rAd,lAk,rAk
real(realk)   :: factor(OrbLoc%norb)
real(realk),pointer :: temp(:,:)
real(realk),pointer :: lAki(:,:),lAi(:,:),rAki(:,:),rAi(:,:)
integer :: A,i
real(realk) :: DDOT

norb   = OrbLoc%norb
natoms = OrbLoc%natoms
m      = OrbLoc%m
OrbLoc%Qkappa = 0.0d0

do A=1,natoms
   factor=OrbLoc%Q(:,A)*OrbLoc%d_m(:)
   call mem_alloc(temp,OrbLoc%nbasA(A),norb)
   call mat_init(lA,OrbLoc%nbasA(A),norb)
   call mat_init(rA,OrbLoc%nbasA(A),norb)
   call mat_init(lAd,norb,OrbLoc%nbasA(A))
   call mat_init(rAd,norb,OrbLoc%nbasA(A))
   call mat_init(lAk,OrbLoc%nbasA(A),norb)
   call mat_init(rAk,OrbLoc%nbasA(A),norb)
   call mat_retrieve_block(OrbLoc%CMO,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1) 
   call mat_set_from_full(temp,1.0d0,lA)
   call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
   call mat_set_from_full(temp,1.0d0,rA)
   call mem_dealloc(temp)
   ! Make lAd=diag(f)*lA^T
   call mat_dmul(factor(:),lA,'T',1.0d0,0.0d0,lAd)
   call mat_dmul(factor(:),rA,'T',1.0d0,0.0d0,rAd)
   ! Make lAk= lA*kappa 
   call mat_mul(lA,kappa,'n','n',1.0d0,0.0d0,lAk)
   call mat_mul(rA,kappa,'n','n',1.0d0,0.0d0,rAk)
   !Term 1
   call mat_mul(lAd,rAk,'n','n',2.0d0*m,1.0d0,Hkappa)
   call mat_mul(rAd,lAk,'n','n',2.0d0*m,1.0d0,Hkappa)
   !Term 2
   call mat_mul(lAk,rAd,'T','T',-2.0d0*m,1.0d0,Hkappa)
   call mat_mul(rAk,lAd,'T','T',-2.0d0*m,1.0d0,Hkappa)

   !Make lAkd= diag(dQ)*lAk^T, use lAd as the matrix lAkd
   call mat_dmul(factor(:),lAk,'T',1.0d0,0.0d0,lAd)
   call mat_dmul(factor(:),rAk,'T',1.0d0,0.0d0,rAd)

   ! Term 3
   call mat_mul(lAd,rA,'n','n',2.0d0*m,1.0d0,Hkappa)
   call mat_mul(rAd,lA,'n','n',2.0d0*m,1.0d0,Hkappa)
   ! Term 4
   call mat_mul(lA,rAd,'T','T',-2.0d0*m,1.0d0,Hkappa)
   call mat_mul(rA,lAd,'T','T',-2.0d0*m,1.0d0,Hkappa)

   !Use lAk and rAk to make diag(Qkappa)
   call mem_alloc(lAki,OrbLoc%nbasA(A),norb)
   call mem_alloc(lAi,OrbLoc%nbasA(A),norb)
   call mem_alloc(rAki,OrbLoc%nbasA(A),norb)
   call mem_alloc(rAi,OrbLoc%nbasA(A),norb)
   call mat_retrieve_block(lAk,lAki,OrbLoc%nbasA(A),norb,1,1)
   call mat_retrieve_block(rAk,rAki,OrbLoc%nbasA(A),norb,1,1)
   call mat_retrieve_block(lA,lAi,OrbLoc%nbasA(A),norb,1,1)
   call mat_retrieve_block(rA,rAi,OrbLoc%nbasA(A),norb,1,1)
   do i=1,norb
      OrbLoc%Qkappa(i,A)=DDOT(OrbLoc%nbasA(A),lAi(:,i),1,rAki(:,i),1)
      OrbLoc%Qkappa(i,A)=0.5d0*(OrbLoc%Qkappa(i,A)+&
      & DDOT(OrbLoc%nbasA(A),rAi(:,i),1,lAki(:,i),1))
   end do
 
   call mem_dealloc(lAki)
   call mem_dealloc(lAi)
   call mem_dealloc(rAki)
   call mem_dealloc(rAi)
   call mat_free(lA);call mat_free(rA)
   call mat_free(lAd);call mat_free(rAd)
   call mat_free(lAk);call mat_free(rAk)
end do

end subroutine CL_add_terms_1_4M


subroutine CL_add_terms_7_8L(OrbLoc,Hkappa)
implicit none
type(PMitem)      :: OrbLoc
integer                 :: natoms,norb
type(matrix)            :: Hkappa
type(matrix)            :: lA,lAd
real(realk),pointer :: temp(:,:) 
real(realk),pointer :: factor(:)
integer     :: i,A
real(realk) :: m

norb = OrbLoc%norb
natoms = OrbLoc%natoms
m = OrbLoc%m
call mem_alloc(factor,norb) 
do A=1,natoms
   call mat_init(lA,OrbLoc%nbasA(A),norb)
   call mat_init(lAd,norb,OrbLoc%nbasA(A))
   factor=OrbLoc%Qkappa(:,A)*OrbLoc%d_m
   !Construct lA and rA
   call mem_alloc(temp,OrbLoc%nbasA(A),norb)
   call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
   call mat_set_from_full(temp,1.0d0,lA)
   call mem_dealloc(temp)
   ! Make lAd=diag(f)*lA^T
   call mat_dmul(factor(:),lA,'T',1.0d0,0.0d0,lAd)
   
   call mat_mul(lAd,lA,'n','n',8.0d0*m,1.0d0,Hkappa)
   call mat_mul(lA,lAd,'T','T',-8.0d0*m,1.0d0,Hkappa)
   call mat_free(lA);call mat_free(lAd)
end do 
call mem_dealloc(factor)
end subroutine CL_add_terms_7_8L

subroutine CL_add_terms_7_8M(OrbLoc,Hkappa)
implicit none
type(PMitem)      :: OrbLoc
integer                 :: natoms,norb
type(matrix)            :: Hkappa
type(matrix)            :: lA,rA,lAd,rAd
real(realk),pointer :: temp(:,:) 
real(realk),pointer :: factor(:)
integer     :: i,A
real(realk) :: m

norb = OrbLoc%norb
natoms = OrbLoc%natoms
m = OrbLoc%m
call mem_alloc(factor,norb) 
do A=1,natoms
   call mat_init(lA,OrbLoc%nbasA(A),norb)
   call mat_init(rA,OrbLoc%nbasA(A),norb)
   call mat_init(lAd,norb,OrbLoc%nbasA(A))
   call mat_init(rAd,norb,OrbLoc%nbasA(A))
   factor=OrbLoc%Qkappa(:,A)*OrbLoc%d_m
   !Construct lA and rA
   call mem_alloc(temp,OrbLoc%nbasA(A),norb)
   call mat_retrieve_block(OrbLoc%CMO,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
   call mat_set_from_full(temp,1.0d0,lA)
   call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
   call mat_set_from_full(temp,1.0d0,rA)
   call mem_dealloc(temp)
   ! Make lAd=diag(f)*lA^T
   call mat_dmul(factor(:),lA,'T',1.0d0,0.0d0,lAd)
   call mat_dmul(factor(:),rA,'T',1.0d0,0.0d0,rAd)
   
   call mat_mul(lAd,rA,'n','n',4.0d0*m,1.0d0,Hkappa)
   call mat_mul(rAd,lA,'n','n',4.0d0*m,1.0d0,Hkappa)
   call mat_mul(lA,rAd,'T','T',-4.0d0*m,1.0d0,Hkappa)
   call mat_mul(rA,lAd,'T','T',-4.0d0*m,1.0d0,Hkappa)
   call mat_free(lA);call mat_free(rA)
   call mat_free(lAd);call mat_free(rAd)
end do 
call mem_dealloc(factor)
end subroutine CL_add_terms_7_8M

subroutine CL_add_terms_9_10L(OrbLoc,Hkappa)
implicit none
type(PMitem)      :: OrbLoc
integer                 :: natoms,norb
type(matrix)            :: Hkappa
type(matrix)            :: lA,lAd
real(realk),pointer :: temp(:,:) 
real(realk),pointer :: factor(:)
real(realk) :: m
real(realk) :: prefactor
integer     :: i,A,B

norb = OrbLoc%norb
natoms= OrbLoc%natoms
m = OrbLoc%m
prefactor= 16.0d0*m*(m+1.0d0)

call mem_alloc(factor,norb)
do A=1,natoms
   do B=1,natoms
        call mat_init(lA,OrbLoc%nbasA(A),norb)
        call mat_init(lAd,norb,OrbLoc%nbasA(A))
        factor=OrbLoc%d_m2*OrbLoc%Q(:,B)*OrbLoc%Q(:,A)*OrbLoc%Qkappa(:,B)
        !Construct lA and rA
        call mem_alloc(temp,OrbLoc%nbasA(A),norb)
        call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
        call mat_set_from_full(temp,1.0d0,lA)
        call mem_dealloc(temp)
        ! Make lAd=diag(f)*lA^T
        call mat_dmul(factor(:),lA,'T',1.0d0,0.0d0,lAd)
        
        call mat_mul(lAd,la,'n','n',-prefactor,1.0d0,Hkappa)
        call mat_mul(lA,lAd,'T','T',prefactor,1.0d0,Hkappa)
        call mat_free(lA)
        call mat_free(lAd)
    end do
end do
call mem_dealloc(factor)

end subroutine CL_add_terms_9_10L


subroutine CL_add_terms_9_10M(OrbLoc,Hkappa)
implicit none
type(PMitem)      :: OrbLoc
integer                 :: natoms,norb
type(matrix)            :: Hkappa
type(matrix)            :: lA,rA,lAd,rAd
real(realk),pointer :: temp(:,:) 
real(realk),pointer :: factor(:)
real(realk) :: m
real(realk) :: prefactor
integer     :: i,A,B

norb = OrbLoc%norb
natoms= OrbLoc%natoms
m = OrbLoc%m
prefactor= 8.0d0*m*(m+1.0d0)

call mem_alloc(factor,norb)
do A=1,natoms
   do B=1,natoms
        call mat_init(lA,OrbLoc%nbasA(A),norb)
        call mat_init(rA,OrbLoc%nbasA(A),norb)
        call mat_init(lAd,norb,OrbLoc%nbasA(A))
        call mat_init(rAd,norb,OrbLoc%nbasA(A))
        factor=OrbLoc%d_m2*OrbLoc%Q(:,B)*OrbLoc%Q(:,A)*OrbLoc%Qkappa(:,B)
        !Construct lA and rA
        call mem_alloc(temp,OrbLoc%nbasA(A),norb)
        call mat_retrieve_block(OrbLoc%CMO,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
        call mat_set_from_full(temp,1.0d0,lA)
        call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
        call mat_set_from_full(temp,1.0d0,rA)
        call mem_dealloc(temp)
        ! Make lAd=diag(f)*lA^T
        call mat_dmul(factor(:),lA,'T',1.0d0,0.0d0,lAd)
        call mat_dmul(factor(:),rA,'T',1.0d0,0.0d0,rAd)
        
        call mat_mul(lAd,ra,'n','n',-prefactor,1.0d0,Hkappa)
        call mat_mul(rAd,lA,'n','n',-prefactor,1.0d0,Hkappa)
        call mat_mul(lA,rAd,'T','T',prefactor,1.0d0,Hkappa)
        call mat_mul(rA,lAd,'T','T',prefactor,1.0d0,Hkappa)
        call mat_free(lA);call mat_free(rA)
        call mat_free(lAd);call mat_free(rAd)
    end do
end do
call mem_dealloc(factor)

end subroutine CL_add_terms_9_10M



subroutine Lowdin_add_terms_1_4(OrbLoc,kappa,Hkappa)
implicit none
type(PMitem)  :: OrbLoc
integer             :: natoms,norb,m
type(matrix)        :: kappa,Hkappa !Hkappa: result
type(matrix)        :: lA,lAd,lAk
real(realk),pointer :: factor(:)
real(realk),pointer :: temp(:,:)
real(realk),pointer :: lAki(:,:),lAi(:,:)
integer :: A,i
real(realk) :: DDOT

norb   = OrbLoc%norb
natoms = OrbLoc%natoms
m      = OrbLoc%m
OrbLoc%Qkappa = 0.0d0

call mem_alloc(factor,norb)
do A=1,natoms
   do i=1,norb;factor(i)=OrbLoc%Q(i,A)**(m-1.0d0); end do
   call mem_alloc(temp,OrbLoc%nbasA(A),norb)
   call mat_init(lA,OrbLoc%nbasA(A),norb)
   call mat_init(lAd,norb,OrbLoc%nbasA(A))
   call mat_init(lAk,OrbLoc%nbasA(A),norb)
   call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1) 
   call mat_set_from_full(temp,1.0d0,lA)
   call mem_dealloc(temp)
   ! Make lAd=diag(f)*lA^T
   call mat_dmul(factor(:),lA,'T',1.0d0,0.0d0,lAd)
   ! Make lAk= lA*kappa 
   call mat_mul(lA,kappa,'n','n',1.0d0,0.0d0,lAk)
   
   ! Term 1
   call mat_mul(lAd,lAk,'n','n',2.0d0*m,1.0d0,Hkappa)
   ! Term 2
   call mat_mul(lAk,lAd,'T','T',-2.0d0*m,1.0d0,Hkappa)
   !Make lAkd= diag(Q)*lAk^T, use lAd as the matrix lAkd
   call mat_dmul(factor(:),lAk,'T',1.0d0,0.0d0,lAd)

   ! Term 3
   call mat_mul(lA,lAd,'T','T',-2.0d0*m,1.0d0,Hkappa)
   ! Term 4
   call mat_mul(lAd,lA,'n','n',2.0d0*m,1.0d0,Hkappa)

   !Use lAk and rAk to make diag(Qkappa)
   call mem_alloc(lAki,OrbLoc%nbasA(A),norb)
   call mem_alloc(lAi,OrbLoc%nbasA(A),norb)
   call mat_retrieve_block(lAk,lAki,OrbLoc%nbasA(A),norb,1,1)
   call mat_retrieve_block(lA,lAi,OrbLoc%nbasA(A),norb,1,1)
   do i=1,norb
       OrbLoc%Qkappa(i,A)=DDOT(OrbLoc%nbasA(A),lAi(:,i),1,lAki(:,i),1)
   end do
 
   call mem_dealloc(lAki)
   call mem_dealloc(lAi)
   call mat_free(lA)
   call mat_free(lAd)
   call mat_free(lAk)
end do

call mem_dealloc(factor)
end subroutine Lowdin_add_terms_1_4

subroutine add_grad_terms(Grad,kappa,Hkappa,Lowdin)
implicit none
type(matrix) :: Grad,kappa,Hkappa
real(realk)  :: factor
logical :: lowdin

if (Lowdin) then
    factor = -0.5d0
else
    factor = -0.5d0
end if

call mat_mul(kappa,Grad,'n','n',-factor,1.0d0,Hkappa)
call mat_mul(Grad,kappa,'T','T',factor,1.0d0,Hkappa)

end subroutine add_grad_terms

subroutine Lowdin_add_terms_7_8(OrbLoc,Hkappa)
implicit none
type(PMitem)      :: OrbLoc
integer                 :: natoms,norb
type(matrix)            :: Hkappa
type(matrix)            :: lA,lAd,lAk
real(realk),pointer :: temp(:,:) 
real(realk),pointer :: factor(:)
real(realk) :: prefac
integer     :: i,A
real(realk) :: m

norb = OrbLoc%norb
natoms = OrbLoc%natoms
m = OrbLoc%m

call mem_alloc(factor,norb) 
prefac=4.0d0*m*(m-1.0d0)

do A=1,natoms
   call mat_init(lA,OrbLoc%nbasA(A),norb)
   call mat_init(lAd,norb,OrbLoc%nbasA(A))
   do i=1,norb;factor(i)=OrbLoc%Qkappa(i,A)*(OrbLoc%Q(i,A)**(m-2.0d0));end do
   !Construct lA and rA
   call mem_alloc(temp,OrbLoc%nbasA(A),norb)
   call mat_retrieve_block(OrbLoc%SC,temp,OrbLoc%nbasA(A),norb,OrbLoc%Pos(A),1)
   call mat_set_from_full(temp,1.0d0,lA)
   call mem_dealloc(temp)
   ! Make lAd=diag(f)*lA^T
   call mat_dmul(factor(:),lA,'T',1.0d0,0.0d0,lAd)
   
   call mat_mul(lAd,lA,'n','n',prefac,1.0d0,Hkappa)
   call mat_mul(lA,lAd,'T','T',-prefac,1.0d0,Hkappa)
   call mat_free(lA);call mat_free(lAd)
end do 
call mem_dealloc(factor)

end subroutine Lowdin_add_terms_7_8


end module Pipek


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



subroutine CLLinearTrans(Grad,HKappa,kappa,OrbLoc)
use loc_types
use Pipek
implicit none
type(PMitem)   :: OrbLoc
type(matrix)   :: S,SC
type(matrix)   :: Hkappa
type(matrix)   :: kappa
type(matrix)   :: Grad
type(lsitem)   :: ls
integer :: k,l

call mat_zero(Hkappa)

if (Orbloc%ChargeLocMulliken) then
  call CL_add_terms_1_4M(OrbLoc,kappa,Hkappa) 
  call CL_add_terms_7_8M(OrbLoc,Hkappa)
  call CL_add_terms_9_10M(OrbLoc,Hkappa)
elseif (OrbLoc%ChargeLocLowdin) then
  call CL_add_terms_1_4L(OrbLoc,kappa,Hkappa) 
  call CL_add_terms_7_8L(OrbLoc,Hkappa)
  call CL_add_terms_9_10L(OrbLoc,Hkappa)
end if

if (OrbLoc%m<0) call mat_scal(-1.0d0,Hkappa)

call add_grad_terms(Grad,kappa,Hkappa,OrbLoc%ChargeLocLowdin)


end subroutine CLLinearTrans

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
   factor(:)=(OrbLoc%Q(:,A))**(int(m)-1)
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

   factor(:)=(OrbLoc%Q(:,A))**(int(m)-2)*OrbLoc%Qkappa(:,A)
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
   factor(:)=OrbLoc%Q(:,A)**(int(m)-1)
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

   factor(:)=OrbLoc%Q(:,A)**(int(m)-2)*OrbLoc%Qkappa(:,A)
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

