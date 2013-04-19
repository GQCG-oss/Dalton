module orbspread_module
!##########################################################
!#              ORBSPREAD MODULE                          #
!# Routines that are specific for orbspread localization  #
!# Routines called by solver (lin.trans. and precond.)    #
!# are not included in module.                            #
!#                                                        #
!##########################################################
use precision
use loc_utils
use typedef
!use typedeftype
use loc_types
use matrix_module, only: matrix
use matrix_operations 
use matrix_operations_aux
use arhDensity
use kurtosis
use davidson_settings
use driver
use LSTIMING
use ARHmodule
use trustradius_mod, only: update_trustradius
use integralInterfaceMod
use decompMod !orbspread_data
use trustradius_mod


contains

!> \brief Routine that drives macro iterations for localizing orbitals using FM
!> \author Ida-Marie Hoeyvik
subroutine PFM_localize_davidson(CFG,CMO,m,ls)
implicit none
type(RedSpaceItem)           :: CFG
type(Matrix) , intent(inout ):: CMO
TYPE(lsitem) , intent(inout) :: ls
integer      , intent(in)    :: m
type(orbspread_data) :: orbspread_input
logical :: lower2
type(Matrix) :: Xsav,CMOsav
type(Matrix), target  ::  X, P, G,expX
integer :: norb, i, nbas
real(realk) :: nrmG, oVal,old_oVal, r, max_step,max_FM
real(realk) :: nrm_thresh,stepsize
real(realk) :: trial(1,1)
  
  write(CFG%lupri,'(a)') '***** PFM SETTINGS *****'
  write(CFG%lupri,'(a,i4)') ' m =        ', m
  write(CFG%lupri,'(a,l4)') ' Crossterms ', CFG%PFM_input%crossterms
  write(CFG%lupri,'(a)') '*****************************'
  
  
  r=0.d0
  norb=CMO%ncol
  nbas=CMO%nrow
  call mat_init(X,norb,norb)
  call mat_init(G,norb,norb)
  call mat_init(P,norb,norb)
  call mat_init(CMOsav,CMO%nrow,CMO%ncol)
  call mat_init(expX,CMO%ncol,CMO%ncol)
  call kurt_initMO(CFG%PFM_input,cmo)
  call kurt_value(CFG%PFM_input)
  call compute_gradient(CFG%PFM_input,G,norb)
  oVal = CFG%PFM_input%kurt_val

  call kurtosis_precond_matrix(CFG%PFM_input,P)
  CFG%PFM_input%P => P
  CFG%P => CFG%PFM_input%P
  
  CFG%mu = 0.0_realk
  if (norb < 15) CFG%macro_thresh = CFG%macro_thresh*10.0d0 
  
  lower2= .true.
  stepsize=0d0
  do i=1,200
    CFG%old_mu = CFG%mu
    old_oVal = oVal
    nrmG = dsqrt(mat_sqnorm2(G))/real(norb)
    max_FM=maxval(CFG%PFM_input%omega)
write(ls%lupri,'(I3,A,ES10.2,A,ES10.3,A,ES10.2,A,ES10.2,A,ES10.3,A,I3,A,f10.3,A,f10.3)') &
     &i, ' Pred=',CFG%r_denom,' max_FM=',max_FM,&
     &   ' r=',r,' mu=',CFG%mu,' grd=', nrmG, ' it=',CFG%it, ' step ', CFG%stepsize&
     & ,' totstep',stepsize
    
    if( nrmG.le. CFG%macro_thresh*100.0) then
       exit
    end if
   
   call solver(CFG,G,X)

   ! global and local thresholds defined in CFG settings
   if (dabs(CFG%mu)> 1.0) CFG%conv_thresh=CFG%global_conv_thresh
   if (dabs(CFG%mu)< 1.0)  CFG%conv_thresh=CFG%local_conv_thresh

   if (dabs(CFG%mu)<100.0 .and. lower2) then
      CFG%max_stepsize = 0.5d0*CFG%max_stepsize
      lower2 = .false.
   end if
   call mat_copy(1d0,cmo,cmosav)

   call linesearch_kurtosis(CFG,cmo,X,stepsize,oval) 
   
   ! COMPUTE r FOR value where factor =1d0
   r=2.0d0*(oVal-old_oVal)/CFG%r_denom
    oVal=CFG%PFM_input%kurt_val
    
    if (oVal-old_oVal < 0) then 
       write(CFG%lupri,*) "Pred: step accepted"
       CFG%stepsize=min(CFG%stepsize*1.2,CFG%max_stepsize)
       !CMOS are updated in linesearch
    else
       write(CFG%lupri,*) "Pred: step rejected"
       call mat_copy(1d0,cmosav,cmo)
       call kurt_updateAO(CFG%PFM_input,cmo)
       call kurt_value(CFG%PFM_input)
       old_oVal=CFG%PFM_input%kurt_val
       oVal=CFG%PFM_input%kurt_val
       CFG%stepsize=CFG%stepsize*0.5
       if (CFG%stepsize < 0.001) then
           write(CFG%lupri,'(a)') 'WARNING: Too many rejections for localization. We exit..' 
	   exit
       end if
       cycle
    end if
    if (CFG%stepsize < 0.001) then
         write(CFG%lupri,'(a)') 'WARNING: Stepsize very small --> cannot converge gradient norm further' 
	 exit
    end if
    call update_trustradius_david(CFG,r,ls,i) 
    !new gradient
    call compute_gradient(CFG%PFM_input,G,norb)

    !new preconditioning matrix 
    call kurtosis_precond_matrix(CFG%PFM_input,P)
    CFG%PFM_input%P => P
    CFG%P => CFG%PFM_input%P


  enddo


  call kurt_freeMO(CFG%PFM_input)
  call mat_free(X)
  call mat_free(expX)
  call mat_free(G)
  call mat_free(P)
  call mat_free(CMOsav)

end subroutine PFM_localize_davidson

!> \brief Routine that drives macro iterations for localizing using SM
!> \author Ida-Marie Hoeyvik
subroutine orbspread_localize_davidson(CFG,CMO,m,orbspread_input,ls)
implicit none
type(RedSpaceItem)           :: CFG
type(Matrix) , intent(inout ):: CMO
TYPE(lsitem) , intent(inout) :: ls
integer      , intent(in)    :: m
type(orbspread_data), target :: orbspread_input
type(Matrix) :: CMOsav
type(Matrix), target  ::  X, P, G
integer :: norb, i,imx,idamax
real(realk) :: nrmG, oVal,old_oVal, r 
real(realk) :: nrm_thresh,stepsize,orig_Eval


  r=0.d0
  norb=CMO%ncol
  CFG%orbspread_input=>orbspread_input

  call mat_init(X,norb,norb)
  call mat_init(G,norb,norb)
  call mat_init(P,norb,norb)
  call mat_init(CMOsav,CMO%nrow,CMO%ncol)


  call orbspread_init(orbspread_input,m,norb)
  call orbspread_update(orbspread_input,CMO)
  call orbspread_gradx(G,norb,orbspread_input)
  call orbspread_value(oVal,orbspread_input)


  call orbspread_precond_matrix2(orbspread_input,P,norb)
  orbspread_input%P => P
  CFG%P => orbspread_input%P

  stepsize = CFG%stepsize
  CFG%mu = 0.0_realk
  if (norb < 10) CFG%macro_thresh = CFG%macro_thresh*10.0d0 
  stepsize=0d0
  CFG%r_denom =1.0_realk
  CFG%it = 1
  do i=1,200
    CFG%old_mu = CFG%mu
    old_oVal = oVal
    imx  =  idamax(norb,orbspread_input%spread2,1)
    nrmG = sqrt(mat_sqnorm2(G))/real(norb)

    write (ls%lupri,'(I3,A,ES9.2,A,f8.3,A,f8.3,A,ES9.2,A,ES9.2,A,I3,A,f6.3,A,f6.3)') &
         &i, 'Pred=',CFG%r_denom,' max_o=',sqrt(orbspread_input%spread2(imx)),&
         &   ' r=',r,' mu=',CFG%mu,' grd=', nrmG, ' it=',CFG%it, ' step',CFG%stepsize,' tot step ', stepsize
    if(nrmG .le. CFG%macro_thresh) exit

   
   call solver(CFG,G,X)

   ! global and local thresholds defined in CFG settings
   if (dabs(CFG%mu)> 1.0) CFG%conv_thresh=CFG%global_conv_thresh
   if (dabs(CFG%mu)< 1.0)  CFG%conv_thresh=CFG%local_conv_thresh

    call mat_copy(1.0_realk,CMO,CMOsav)
 
    stepsize = CFG%stepsize
    call linesearch_orbspread2(CFG,cmo,X,stepsize,old_oval,orig_Eval)
    call orbspread_value(oVal,orbspread_input)
    
    r=2.0d0*(orig_Eval-old_oVal)/CFG%r_denom

    if (orig_Eval-old_oVal > 0) then
       write(ls%lupri,*) 'Step not accepted. Go back'
       call mat_copy(1.0d0,CMOsav,CMO)
       call orbspread_update(orbspread_input,CMO)
       call orbspread_value(oVal,orbspread_input)
       CFG%Stepsize = CFG%Stepsize/2.0d0
       if (CFG%stepsize < 0.001) then
           write(CFG%lupri,'(a)') 'WARNING: Too many rejections for localization. We exit..' 
	   exit
       end if
	cycle
    else
      CFG%Stepsize = min(CFG%Stepsize*2.5d0,CFG%max_stepsize) 
    endif

   if (CFG%stepsize < 0.001) then
           write(CFG%lupri,'(a)') 'WARNING: Stepsize too small.  We exit..' 
	   exit
   end if
    !new gradient
    call orbspread_gradx(G,norb,orbspread_input)


   call orbspread_precond_matrix2(orbspread_input,P,norb)
   orbspread_input%P => P
   CFG%P =>orbspread_input%P


  enddo

  call orbspread_free(orbspread_input)
  call mat_free(X)
  call mat_free(G)
  call mat_free(P)
  call mat_free(CMOsav)

end subroutine orbspread_localize_davidson




function orbspread_select_gt(CMO,eps,inp)
implicit none
integer, pointer :: orbspread_select_gt(:)
Type(Matrix), intent(in) :: CMO
real(realk), intent(in)  :: eps
type(orbspread_data), intent(inout) :: inp
integer :: ncol,norb, i, j

    ncol=CMO%ncol
    call orbspread_init(inp,-1,ncol)
    call orbspread_update(inp,CMO)

    norb=0
    do i=1,ncol
     if (sqrt(inp%spread2(i)).gt.eps) norb = norb+1
    enddo

    call mem_alloc(orbspread_select_gt,norb)

    j=0
    do i=1,ncol
     if (sqrt(inp%spread2(i)).gt.eps) then
      j=j+1
      orbspread_select_gt(j) = i
     endif
    enddo

    call orbspread_free(inp)

    return
end function orbspread_select_gt

function orbspread_build_selection(selected,CMO)
implicit none
Type(Matrix), pointer :: orbspread_build_selection
integer, pointer      :: selected(:)
Type(Matrix)          :: CMO
    
integer               :: i, nrow, ncol
real(realk),pointer     :: tmp(:)

    nrow = CMO%nrow
    ncol = size(selected)
    allocate(orbspread_build_selection)
    call mem_alloc(tmp,nrow)
    call mat_init(orbspread_build_selection,nrow,ncol)

    do i=1,ncol    
       call mat_retrieve_block(CMO,tmp,nrow,1,1,selected(i))
       call mat_create_block(orbspread_build_selection,tmp,nrow,1,1,i)
    enddo

    call mem_dealloc(tmp)
 
return
end function orbspread_build_selection

subroutine orbspread_set_selection(CMO,CMOs,selected)
implicit none
Type(Matrix), intent(inout) :: CMO
Type(Matrix), intent(inout) :: CMOs
integer, pointer            :: selected(:) 

integer                     :: i, nrow, ncol
real(realk),pointer     :: tmp(:)

    nrow = CMO%nrow
    ncol = size(selected)

    call mem_alloc(tmp,nrow)

    do i=1,ncol    
       call mat_retrieve_block(CMOs,tmp,nrow,1,1,i)
       call mat_create_block(CMO,tmp,nrow,1,1,selected(i))
    enddo

    call mem_dealloc(tmp)

end subroutine orbspread_set_selection


subroutine orbspread_free(orbspread_input)
implicit none
type(orbspread_data) :: orbspread_input
integer :: i

     call mem_dealloc(orbspread_input%spread2) 
 
     call mat_free(orbspread_input%Q)

     do i=1,3
      call mat_free(orbspread_input%R(i))
     enddo

     do i=1,4
      call mat_free(orbspread_input%tmpM(i))
     enddo

end subroutine orbspread_free

subroutine orbspread_propint_free(orbspread_input)
implicit none
type(orbspread_data) :: orbspread_input
integer :: i

     do i=2,5
      call mat_free(orbspread_input%propint(i))
     enddo

end subroutine orbspread_propint_free


subroutine orbspread_propint(orbspread_input,ls)
implicit none
type(orbspread_data), intent(inout) :: orbspread_input
TYPE(lsitem), intent(inout) :: ls !ls%setting is changed in the II_something routines

integer, parameter  :: nderiv=2, nMAT=10
integer i,n

! init and compute propint property integrals
! propint(2:4) -> DIPX,DIPY,DIPZ
! propint(5)   -> SECX+SECY+SECZ 

!     n   = ls%setting%BASIS(1)%p%REGULAR%nbast
     n = getNbasis(AORdefault,Contractedinttype,ls%SETTING%MOLECULE(1)%p,ls%lupri)
     DO i=1,nMAT
        call mat_init(orbspread_input%propint(I),n,n)
     ENDDO

     IF(.FALSE.)THEN !speed optimized version
        call II_get_carmom(6,6,ls%setting,orbspread_input%propint,nMAT,nderiv,0E0_realk,0E0_realk,0E0_realk)
     ELSE !memory optimized version
        do I=1,nMat
           call II_get_single_carmom(6,6,ls%setting,orbspread_input%propint(i),i,nderiv,0E0_realk,0E0_realk,0E0_realk)
        enddo
     ENDIF
     call mat_free(orbspread_input%propint(1))
     call mat_free(orbspread_input%propint(6))
     call mat_free(orbspread_input%propint(7))
     call mat_free(orbspread_input%propint(9))

! propint(5)   -> SECX+SECY+SECZ 
     call mat_daxpy(1E0_realk,orbspread_input%propint(8),orbspread_input%propint(5))
     call mat_free(orbspread_input%propint(8))
     call mat_daxpy(1E0_realk,orbspread_input%propint(10),orbspread_input%propint(5))
     call mat_free(orbspread_input%propint(10))

end subroutine orbspread_propint



subroutine orbspread_init(orbspread_input,m,norb)
implicit none
type(orbspread_data), intent(inout) :: orbspread_input
integer             , intent(in) :: norb, m
TYPE(lsitem) :: ls
integer   :: i

! init R Q and tmpM matrices
     call mat_init(orbspread_input%Q,norb,norb)
     do i=1,3
     call mat_init(orbspread_input%R(i),norb,norb)
     enddo
     do i=1,4
     call mat_init(orbspread_input%tmpM(i),norb,norb)
     enddo

! set norb
     orbspread_input%norb =  norb

! set m power
     orbspread_input%m = m

! allocate spread2

     call mem_alloc(orbspread_input%spread2,norb)

end subroutine orbspread_init


subroutine orbspread_update(orbspread_input,CMO)
implicit none
type(orbspread_data), intent(inout) :: orbspread_input
type(Matrix), intent(in) :: CMO

type(Matrix) :: tmp
integer :: nbas,norb,i
real(realk), pointer :: tmpv(:)

     nbas=CMO%nrow
     norb=CMO%ncol

     call mat_init(tmp,nbas,norb)

!   R(1:3)
    do i=1,3
     call mat_mul(orbspread_input%propint(i+1),CMO,'n','n',1E0_realk,0E0_realk,tmp)
     call mat_mul(CMO,tmp,'t','n',1E0_realk,0E0_realk,orbspread_input%R(i))
    enddo

!   Q
     call mat_mul(orbspread_input%propint(5),CMO,'n','n',1E0_realk,0E0_realk,tmp)
     call mat_mul(CMO,tmp,'t','n',1E0_realk,0E0_realk,orbspread_input%Q)

     call mat_free(tmp)

!   spread2
    call mem_alloc(tmpv,norb)

    call mat_extract_diagonal(orbspread_input%spread2,orbspread_input%Q)

    do i=1,3
    call mat_extract_diagonal(tmpv,orbspread_input%R(i))
    tmpv=tmpv**2
    call daxpy(norb,-1E0_realk,tmpv,1,orbspread_input%spread2,1)
    enddo
   
    !write(6,*) sqrt(orbspread_input%spread2)
    call mem_dealloc(tmpv)

end subroutine orbspread_update

subroutine orbspread_value(oVal, inp)
real(realk), intent(out) :: oVal
type(orbspread_data), intent(in)  :: inp

integer  :: i

 oVal = 0E0_realk
 do i=1,inp%Q%ncol
   oVal = oVal + inp%spread2(i)**inp%m
 enddo


end subroutine orbspread_value




subroutine orbspread_update_trustradius(arh,r)
type(solverItem), intent(inout)   :: arh
real(realk)     , intent(in)      :: r

real(realk)                       :: nom, denom

 !nom =   (oVal - arh%old_energy)
 !denom=  arh%denom

 !r = nom/denom 

 if (r.lt. 0E0_realk) then 
    arh%set_max_element = arh%set_max_element/2.0E0_realk
    arh%set_max_step = arh%set_max_step/2.0E0_realk
    arh%step_accepted = .false.
    write (6,*) 'Reject and contract *0.5'
    return
 endif

 if(r.gt. 0.75E0_realk) then
    arh%set_max_element = arh%set_max_element*1.2E0_realk
    arh%set_max_step = arh%set_max_step*1.2E0_realk
    write(6,*) 'Expand *1.2', r
 else if(r.gt. 0.25) then
    write(6,*) 'Keep', r
    continue
 else 
    arh%set_max_element = arh%set_max_element*0.7E0_realk
    arh%set_max_step = arh%set_max_step*0.7E0_realk
    write(6,*) 'Contract *0.7', r
 endif

 arh%step_accepted = .true.

end subroutine orbspread_update_trustradius




subroutine orbspread_gradx(G,norb,inp)
implicit none
Type(Matrix) , target, intent(inout) :: G
integer                              :: norb
type(orbspread_data), intent(inout)  :: inp

real(realk)  :: diagR(norb,3), tmp(norb)
integer      :: x, m

  m=inp%m

  do x=1, 3
  call mat_extract_diagonal(diagR(:,x),inp%R(x))
  enddo

  tmp =  (inp%spread2**(m-1))
  call mat_zero(G)
  call mat_dmul(tmp,inp%Q,'n',-2E0_realk*m,0E0_realk,G)

  do x=1,3
  tmp = diagR(:,x)*(inp%spread2**(m-1))
  call mat_dmul(tmp,inp%R(x),'n',4E0_realk*m,1E0_realk,G)
  enddo

  call mat_trans(G,inp%tmpM(1))
  call mat_daxpy(-1E0_realk,inp%tmpM(1),G)

  !call mat_scal(0.5E0_realk,G)

  inp%G => G

end subroutine  orbspread_gradx
  
!> \brief Compute diagonal hessian elements for variance localization
!> \author Ida-Marie Hoeyvik
subroutine orbspread_precond_matrix2(inp,P,norb)
implicit none
type(orbspread_data), intent(in) :: inp
type(matrix) :: P
integer :: m,norb
real(realk) :: diagQ(norb),diagR(norb,3)
real(realk) :: spm1(norb),spm2(norb)
integer :: i,x
type(matrix) :: tmp,tmp1
real(realk),pointer :: tmpvec(:)

call mat_zero(P)
m=inp%m

call mat_extract_diagonal(diagQ,inp%Q)
do x=1,3
  call mat_extract_diagonal(diagR(:,x),inp%R(x))
end do

spm1= inp%spread2**(m-1)
if (m>1) spm2=inp%spread2**(m-2)
call mat_dger(dble(2*m),diagQ,spm1,P)

do x=1,3
  call mat_dger(-dble(4*m),spm1*diagR(:,x),diagR(:,x),P)
end do

call mat_init(tmp,P%nrow,P%ncol)
call mat_zero(tmp)
if (m>1) then
  call mat_hmul(1d0,inp%Q,inp%Q,0d0,tmp)
  call mat_dmul(spm2,tmp,'n',dble(4*m*(m-1)),1d0,P)
end if
do x=1,3
  call mat_hmul(1d0,inp%R(x),inp%R(x),0d0,tmp)
  call mat_dmul(spm1,tmp,'n',-dble(8*m),1d0,P)
  if (m>1) then
    call mat_hmul(1d0,inp%Q,inp%R(x),0d0,tmp)
    call mat_dmul(spm2*diagR(:,x),tmp,'n',-dble(16*m*(m-1)),1d0,P)
  end if
end do

call mat_zero(tmp)
if (m>1) then
   do x=1,3
     call mat_dmul(diagR(:,x),inp%R(x),'n',1d0,1d0,tmp)
   end do
   call mat_init(tmp1,P%nrow,P%ncol)
   call mat_dmul(spm2,tmp,'n',dble(16*m*(m-1)),0d0,tmp1)
   call mat_hmul(1d0,tmp,tmp1,1d0,P)
   call mat_free(tmp1)
end if

call mem_alloc(tmpvec,norb)
tmpvec=1.0d0
call mat_dger(-dble(2*m),diagQ*spm1,tmpvec,P)

do x=1,3
   call mat_dger(dble(4*m),spm1*diagR(:,x)*diagR(:,x),tmpvec,P)
end do
call mem_dealloc(tmpvec)

call mat_trans(P,tmp)
call mat_daxpy(1d0,tmp,P)
call mat_free(tmp)
call mat_scal_dia(0d0,P)

end subroutine orbspread_precond_matrix2


subroutine linesearch_kurtosis(CFG,cmo,X,stepsize,oval)
implicit none
type(RedSpaceItem) :: CFG
type(matrix)  :: cmo,X
integer :: i,numb,nmats
type(matrix)  :: cmotemp(15),Xtemp(15)
real(realk) :: old_funcval,factor,step(15),stepsize,oval

numb=15
factor = 1.0d0
step = 1.0d0
step(1)=0.0d0
step(11)=1.50
step(12)=1.50
step(13)=1.50
step(14)=1.50
step(15)=1.50
    old_funcval = CFG%PFM_input%kurt_val
    
    if (CFG%orb_debug) write(CFG%lupri,'(a,I4,a,ES13.3)') &
    &'Linesearch number :', 0, ' Original function value: ', CFG%PFM_input%kurt_val
do i=1,numb
    call mat_init(Xtemp(i),X%nrow,X%ncol)
    call mat_assign(Xtemp(i),X)
    call mat_init(cmotemp(i),cmo%nrow,cmo%ncol)
    call mat_assign(cmotemp(i),cmo)
    factor = factor + step(i)
    call mat_scal(factor,Xtemp(i))
    call updatecmo(CMOtemp(i),Xtemp(i))      
    call kurt_updateAO(CFG%PFM_input,CMOtemp(i))
    call kurt_value(CFG%PFM_input)
    if (CFG%orb_debug) write(CFG%lupri,'(a,I4,a,ES13.3)') &
    &'Linesearch number :', i, ' Change ', CFG%PFM_input%kurt_val-old_funcval
    if (i==1) oVal= CFG%PFM_input%kurt_val
    if ((CFG%PFM_input%kurt_val > old_funcVal) .and. i>1) then
           call mat_assign(cmo,cmotemp(i-1))
	   stepsize = dsqrt(mat_dotproduct(xtemp(i-1),xtemp(i-1)))
	   nmats=i
           exit
    end if
    if (i==numb .or. dabs(CFG%PFM_input%kurt_val-old_funcval)<1d0) then
      call mat_assign(cmo,cmotemp(i))
      stepsize = dsqrt(mat_dotproduct(xtemp(i),xtemp(i)))
      nmats=i
      exit
    end if
    old_funcval=CFG%PFM_input%kurt_val
end do

call kurt_updateAO(CFG%PFM_input,CMO)
call kurt_value(CFG%PFM_input)
do i=1,nmats
  call  mat_free(CMOtemp(i))
  call  mat_free(Xtemp(i))
end do
end subroutine linesearch_kurtosis



subroutine linesearch_orbspread2(CFG,cmo,X,stepsize,value_last_macro,orig_eival)
implicit none
type(RedSpaceItem) :: CFG
type(matrix)  :: cmo,X
real(realk),intent(in) :: value_last_macro
integer :: i,numb=15,nmats
type(matrix)  :: cmotemp(15),Xtemp(15)
real(realk) :: old_funcval,factor,step(15),stepsize,oval
real(realk) :: orig_eival

old_funcval = value_last_macro

factor = 1.0d0
step = 1.0d0
step(1)=0.0d0
step(11)=1.50
step(12)=1.50
step(13)=1.50
step(14)=1.50
step(15)=1.50
if (CFG%orb_debug) write(CFG%lupri,'(a,I4,a,f15.1)') &
&'Linesearch number :', 0, ' Original function value: ', old_funcval
do i=1,numb
    call mat_init(Xtemp(i),X%nrow,X%ncol)
    call mat_copy(1.0d0,X,Xtemp(i))
    call mat_init(cmotemp(i),cmo%nrow,cmo%ncol)
    call mat_copy(1.0d0,cmo,cmotemp(i))
    factor = factor + step(i)
    call mat_scal(factor,Xtemp(i))
    call updatecmo(CMOtemp(i),Xtemp(i))
    call orbspread_update(CFG%orbspread_input,CMOtemp(i))
    call orbspread_value(oVal,CFG%orbspread_input)
    if (i==1) orig_eival = oval
    if (CFG%orb_debug) write(CFG%lupri,'(a,I4,a,f15.4)') &
    &'Linesearch number :', i, ' Change ', oVal-old_funcval
    if (i==1 .and. oVal > old_funcVal) then
       nmats=i
       exit
    endif
    if (oVal > old_funcVal) then
           call mat_copy(1d0,cmotemp(i-1),cmo)
           stepsize = dsqrt(mat_dotproduct(xtemp(i-1),xtemp(i-1)))
	   nmats=i
           exit
    end if
    if (i==numb .or. dabs(oVal-old_funcval)< 1.0) then
      call mat_copy(1d0,cmotemp(i),cmo)
      stepsize = dsqrt(mat_dotproduct(xtemp(i),xtemp(i)))
      nmats=i
      exit
    end if
    old_funcval=oVal
end do


call orbspread_update(CFG%orbspread_input,CMO)
call orbspread_value(oVal,CFG%orbspread_input)

do i=1,nmats
  call  mat_free(CMOtemp(i))
  call  mat_free(Xtemp(i))
end do
end subroutine linesearch_orbspread2

end module orbspread_module
