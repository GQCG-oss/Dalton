
!***************************************************************************
!*                                                                         *
!*  MAIN ROUTINES TO COMPUTE GRADIENT,HESSIAN ON TRIAL VECTOR AND PRECOND  *
!*                                                                         *
!* Implemented for                                                         *
!* - Pipek measure using Lowdin or Mulliken charge and power m             *
!* - Pipek-Mezey with external power m, Lowdin charge only                 *
!* - Generalized Pipek-Mezey, external power m, Lowdin charge              *
!*                                                                         *
!* Lowdin orbitals may be optimized using orbspread scheme                 *
!*                                                                         *
!* NOTE: For Pipek-Mezey --> .PML m=2, Lowdin pop.analysis                 *
!* NOTE: For Pipek-Mezey --> .PMM m=2, Mulliken pop.analysis               *
!*                                                                         *
!* Author: Ida-Marie Hoeyvik                                                *
!***************************************************************************

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


!> \brief Calls correct routine for computing diagonal hessian elements
!> \param P Matrix with diagonal elements (output)
subroutine Precond_ChargeLoc(P,OrbLoc)
implicit none
type(PMitem) :: OrbLoc
type(matrix) :: P


if (OrbLoc%PipekMezeyLowdin) then
   call PMLowdin_precond(OrbLoc,P)
elseif (OrbLoc%PipekMezeyMull) then
   call PMMull_precond(OrbLoc,P)
endif


end subroutine Precond_ChargeLoc


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

subroutine add_grad_terms(Grad,kappa,Hkappa)
implicit none
type(matrix) :: Grad,kappa,Hkappa
real(realk)  :: factor

factor = -0.5d0

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


