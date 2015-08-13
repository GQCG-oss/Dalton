module charge_module
use precision
use loc_utils
use Pipek
use matrix_module
use matrix_operations
use davidson_settings
use orbspread_module
use memory_handling
use TYPEDEF
use TYPEDEFTYPE
use davidson_solv_mod
use LSTIMING
use loc_types
!use linesearch
!##########################################################
!#              CHARGE MODULE                             #
!# Routines that are specific for charge localization.    #
!# Routines called by solver (lin.trans. and precond.)    #
!# are not included in module.                            #
!#                                                        #
!##########################################################
contains


subroutine get_correct_S(CFG,ls,nbas)
implicit none
type(lsitem) :: ls
type(matrix) :: S,Sminussqrt,Stemp
integer      :: nbas,natoms
type(redspaceitem) :: CFG
real(realk),dimension(nbas,nbas) :: Smat,Smat_minussqrt,Smat_sqrt

natoms = ls%setting%MOLECULE(1)%p%nAtoms

call mat_init(S,nbas,nbas)
call mat_init(Sminussqrt,nbas,nbas)
call II_get_overlap(6,6,ls%setting,S)
call mat_to_full(S,1.0_realk,Smat)
call lowdin_diag(nbas,Smat,Smat_sqrt,Smat_minussqrt,6)


if (CFG%PM_input%PipekMezeyLowdin) then
     call mat_set_from_full(Smat_sqrt,1.0_realk,CFG%PM_input%SU)
elseif (CFG%PM_input%PipekMezeyMull) then
     call mat_copy(1d0,S,CFG%PM_input%SU)
end if    

call mat_free(S)
call mat_free(Sminussqrt)

end subroutine get_correct_S


!> \brief ROutine that drives macro iterations for charge localization
!> \author Ida-Marie Hoeyvik
subroutine charge_localize_davidson(CFG,CMOall,m,ls,norb)
implicit none
type(RedSpaceItem)           :: CFG
type(matrix) :: CMOall
type(Matrix) , target        :: CMO
TYPE(lsitem) , intent(in)    :: ls
integer      , intent(in)    :: m
type(Matrix)          :: CMOsav
type(Matrix), target  ::  X, P, G
integer     :: norb, i, imx, A,j,nbas
real(realk) :: nrmG, oVal,old_oVal, r
real(realk) :: fVal,old_fVal !function values
real(realk) :: stepsize
real(realk) :: max_loc,orig_Eval
integer, external :: idamax
integer     :: nel
real(realk) :: minel
integer     :: minel_pos(2)
real(realk), pointer :: tmp(:)

  nbas=CMOall%nrow
  ! Extract coefficients to be localized
  call mem_alloc(tmp,nbas*norb)
  call mat_init(CMO,nbas,norb)
  ! extract matrix from CMOall(1,offset)
  call mat_retrieve_block(CMOall,tmp,nbas,norb,1,CFG%offset)
  call mat_set_from_full(tmp,1.0_realk,CMO)
  call mem_dealloc(tmp)
   


  CFG%PM_input%cmo=>CMO
  r=0.d0
  norb=CMO%ncol
  call mat_init(X,norb,norb)
  call mat_zero(X)
  call mat_init(G,norb,norb)
  call mat_init(P,norb,norb)
  call mat_init(CMOsav,CMO%nrow,CMO%ncol)
  call initialize_OrbLoc(CMO,CFG%PM_input,ls,dble(m))
  call update_OrbLoc(CFG%PM_input,CMO,ls)
  call ComputeGrad_PipekMezey(CFG%PM_input,G)
  call Precond_ChargeLoc(P,CFG%PM_input)
  CFG%PM_input%P => P
  CFG%P=>CFG%PM_input%P
  
  fVal=CFG%PM_input%funcVal
  CFG%mu = 0.0_realk
  stepsize = CFG%stepsize
  if (norb < 10) CFG%macro_thresh=CFG%macro_thresh*10.0d0

  do i=1,CFG%max_macroit
    CFG%old_mu = CFG%mu
    old_fVal = fVal
    nrmG = dsqrt(mat_sqnorm2(G))/real(norb)
    max_loc = fVal/real(norb)

  write (ls%lupri,'(1X,A,I3,A,ES8.1,A,ES8.1,A,ES8.1,A,I2,A,f5.2,A,f5.2)') &
 & ' %LOC% ',i, ' max_loc = ',max_loc, &
 & ' mu = ',CFG%mu,' grd = ', nrmG, ' it = ',CFG%it, ' trust-region =', CFG%stepsize, ' step =',stepsize


   if( nrmG .le. CFG%macro_thresh .and. i.gt.1) exit
   
   call davidson_solver(CFG,G,X)
   if (CFG%PM_input%PipekMezeyMull) call mat_scal(-1d0,X)
   if (CFG%PM_input%PipekMezeyLowdin) call mat_scal(-1d0,X)

   !Dynamic convergence thresh
   if (dabs(CFG%mu)> 1.0) CFG%conv_thresh=CFG%global_conv_thresh
   if (dabs(CFG%mu)< 1.0)  CFG%conv_thresh=CFG%local_conv_thresh


   call mat_copy(1d0,CMO,CMOsav)

   call linesearch_charge(CFG%PM_input,cmo,X,stepsize,orig_eval,ls)
   nullify(CFG%PM_input%CMO)
   CFG%PM_input%CMO=>CMO
  ! Compute r
  !call updatecmo(CMO,X)
  call update_OrbLoc(CFG%PM_input,CMO,ls)
  fval = CFG%PM_input%funcVal
  r=2.0d0*(orig_eval-old_fVal)/CFG%r_denom

  if (r<0.0_realk) then
       write(ls%lupri,*) ' %LOC% Step not accepted. Go back'
       call mat_copy(1.0d0,CMOsav,CMO)
       call update_OrbLoc(CFG%PM_input,CMO,ls)
       fval = CFG%PM_input%funcVal
       CFG%stepsize = CFG%stepsize/2.0_realk
   else
       CFG%stepsize = min(CFG%stepsize*2.5_realk,CFG%max_stepsize)
   endif

   if (CFG%stepsize < 0.0001_realk) then
         write(CFG%lupri,'(a)') ' %LOC% WARNING: Too many rejections for localization. We exit..' 
         exit
    end if

    !new gradient
    call ComputeGrad_PipekMezey(CFG%PM_input,G)
    call Precond_ChargeLoc(P,CFG%PM_input)
    CFG%PM_input%P => P
    CFG%P=>CFG%PM_input%P

  !END OF LOCALIZATION LOOP
  enddo

  !Put localized block into full CMO matrix
  call mem_alloc(tmp,nbas*norb)
  call mat_to_full(CMO,1.0_realk,tmp)
  call mat_free(CMO)
  call mat_create_block(CMOall,tmp,nbas,norb,1,CFG%offset)
  call mem_dealloc(tmp)



  call FreeOrbLoc(CFG%PM_input)
  call mat_free(X)
  call mat_free(G)
  call mat_free(P)
  call mat_free(CMOsav)
end subroutine charge_localize_davidson


subroutine linesearch_charge(OrbLoc,cmo,X,stepsize,orig_eival,ls)
implicit none
type(PMitem) :: OrbLoc
type(lsitem) :: ls
type(matrix)  :: cmo,X
integer :: i,nmats,numb
type(matrix),target  :: cmotemp(6)
type(matrix) :: Xtemp(6)
real(realk) :: old_funcval,factor,step(6),stepsize,oval
real(realk) :: orig_eival

numb=6
old_funcval = OrbLoc%funcval

factor = 2.5d0
step = 1.0d0
step(1)=0.0d0
if (OrbLoc%orb_debug) write(ls%lupri,'(a,I4,a,f15.1)') &
&'Linesearch number :', 0, ' Original function value: ', old_funcval
do i=1,numb
    call mat_init(Xtemp(i),X%nrow,X%ncol)
    call mat_copy(1.0d0,X,Xtemp(i))
    call mat_init(cmotemp(i),cmo%nrow,cmo%ncol)
    call mat_copy(1.0d0,cmo,cmotemp(i))
    factor = factor + step(i)
    call mat_scal(factor,Xtemp(i))
    call updatecmo(CMOtemp(i),Xtemp(i))
    nullify(OrbLoc%cmo)
    OrbLoc%cmo=>CMOtemp(i)
    call update_OrbLoc(OrbLoc,CMOtemp(i),ls)
    oVal = OrbLoc%funcVal
    if (i==1) orig_eival = oval
    if (OrbLoc%orb_debug) write(ls%lupri,'(a,I4,a,f15.4)')&
    &'Linesearch number :', i, ' Change ', oVal-old_funcval
    if (i==1 .and. oVal > old_funcVal) then
       nmats=i
       stepsize = sqrt(mat_sqnorm2(Xtemp(i)))
       exit
    endif
    if (oVal > old_funcVal) then
           call mat_copy(1d0,cmotemp(i-1),cmo)
           stepsize = sqrt(mat_sqnorm2(Xtemp(i-1)))
	   nmats=i
           exit
    end if
    if (i==numb ) then
      stepsize = sqrt(mat_sqnorm2(Xtemp(i)))
      call mat_copy(1d0,cmotemp(i),cmo)
      nmats=i
      exit
    end if
    old_funcval=oVal
end do


do i=1,nmats
  call  mat_free(CMOtemp(i))
  call  mat_free(Xtemp(i))
end do

end subroutine linesearch_charge


end module charge_module
