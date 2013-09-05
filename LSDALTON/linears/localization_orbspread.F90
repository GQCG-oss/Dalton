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
use loc_types
use matrix_module, only: matrix
use matrix_operations 
use matrix_operations_aux
use kurtosis
use davidson_settings
use davidson_solv_mod !davidson_solver
use LSTIMING
use integralInterfaceMod
use decompMod !orbspread_data
use orbspread_utilMod

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
integer :: norb, i, nbas,iter_number
real(realk) :: nrmG, oVal,old_oVal, max_step,max_FM
real(realk) :: nrm_thresh,stepsize
real(realk) :: trial(1,1)
real(realk),pointer :: max_orbspreads(:)  
  
  
  norb=CMO%ncol
  nbas=CMO%nrow
  call mem_alloc(max_orbspreads,CFG%max_macroit)
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
  
  lower2= .true.
  stepsize=0d0
  do i=1,CFG%max_macroit
    iter_number = i
    CFG%old_mu = CFG%mu
    old_oVal = oVal
    nrmG = dsqrt(mat_sqnorm2(G))/real(norb)
    max_FM=sqrt(sqrt(maxval(CFG%PFM_input%omega)))
    max_orbspreads(i) =  max_FM
   write(ls%lupri,'(I3,A,ES8.1,A,f6.2,A,ES8.1,A,ES8.1,A,I3,A,f5.2,A,f5.2)') &
     &i, ' Pred= ',CFG%r_denom,'  sigma_4 =',max_FM,&
     &   '  mu = ',CFG%mu,'  grd =', nrmG, '  it =',CFG%it, '  trust-region = ', CFG%stepsize&
     & ,' step= ',stepsize
    

    if (i>10) then
        if ( abs(max_orbspreads(i)-max_orbspreads(i-10)) < 0.05 ) then
           write(ls%lupri,*) '  '
           write(ls%lupri,*) '    ********* Orbital localization converged ************'
           write(ls%lupri,*) '    *                                                   *'
           write(ls%lupri,*) '    * There are insignificant changes in the  locality  *'
           write(ls%lupri,*) '    * of the least local orbital, and procedure is      *'
           write(ls%lupri,*) '    * exited irresptective of gradient norm.            *'
           write(ls%lupri,*) '    *                                                   *'
           write(ls%lupri,*) '    *****************************************************'
           write(ls%lupri,*) '  '
           exit
        endif
    endif
    if( nrmG.le. CFG%macro_thresh*10.0) then
        write(ls%lupri,*) '  '
        write(ls%lupri,*) '   ********* Orbital localization converged ************'
        write(ls%lupri,*) '   *                                                   *'
        write(ls%lupri,*) '   * The gradient norm for the orbital localization    *'
        write(ls%lupri,*) '   * function is below the threshold, and we exit      *'
        write(ls%lupri,*) '   * the localization procedure.                       *'
        write(ls%lupri,*) '   *                                                   *'
        write(ls%lupri,*) '   *****************************************************'
        write(ls%lupri,*) '  '
        exit
    end if
   
   call davidson_solver(CFG,G,X)

   ! global and local thresholds defined in CFG settings
   if (dabs(CFG%mu)> 1.0) CFG%conv_thresh=CFG%global_conv_thresh
   if (dabs(CFG%mu)< 1.0)  CFG%conv_thresh=CFG%local_conv_thresh

   if (dabs(CFG%mu)<100.0 .and. lower2) then
      CFG%max_stepsize = 0.5d0*CFG%max_stepsize
      lower2 = .false.
   end if
   call mat_copy(1d0,cmo,cmosav)

   call linesearch_kurtosis(CFG,cmo,X,stepsize,oval) 
   
    oVal=CFG%PFM_input%kurt_val
    
    if (oVal-old_oVal < 0) then 
       write(CFG%lupri,*) "Pred: step accepted"
       CFG%stepsize=min(CFG%stepsize*2.0,CFG%max_stepsize)
       !CMOS are updated in linesearch
    else
       write(CFG%lupri,*) "Pred: step rejected"
       call mat_copy(1d0,cmosav,cmo)
       call kurt_updateAO(CFG%PFM_input,cmo)
       call kurt_value(CFG%PFM_input)
       old_oVal=CFG%PFM_input%kurt_val
       oVal=CFG%PFM_input%kurt_val
       CFG%stepsize=CFG%stepsize*0.5
    endif    

   if (CFG%stepsize < 0.001) then
        write(CFG%lupri,*) ''
        write(CFG%lupri,'(a)') 'WARNING: Stepsize too small. ' 
        if (i>5) then
           if  (abs(max_orbspreads(i)-max_orbspreads(i-5))< 0.1) then 
               write(CFG%lupri,*) ' However, the locality of the least local orbital       ' 
               write(CFG%lupri,*) ' has not changed significantly the last five iterations ' 
               write(CFG%lupri,*) ' and the generated orbitals are localized, and will      ' 
               write(CFG%lupri,*) ' be written to file.   '
               write(CFG%lupri,*) ''
	       exit
           endif
        endif
        write(CFG%lupri,*) ' Cannot proceed with localization due to issues with    ' 
        write(CFG%lupri,*) ' solving the level-shifted Newton equations. You may    ' 
        write(CFG%lupri,*) ' try to restart calculation and lower the residual norm ' 
        write(CFG%lupri,*) ' threshold for the micro iterations as described in     '
        write(CFG%lupri,*) ' the user manual under section **LOCALIZE ORBITALS      '
        write(CFG%lupri,*) ' and keyword .LOOSE MICRO THRESH                              ' 
        call lsquit('Cannot converge micro iterations. ', CFG%lupri)
   elseif (oVal-old_oVal < 0) then
            cycle
   endif

    !new gradient
    call compute_gradient(CFG%PFM_input,G,norb)

    !new preconditioning matrix 
    call kurtosis_precond_matrix(CFG%PFM_input,P)
    CFG%PFM_input%P => P
    CFG%P => CFG%PFM_input%P

  enddo
    if (iter_number==CFG%max_macroit) then
        write(CFG%lupri,*) ''
        write(CFG%lupri,*) ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(CFG%lupri,*) ' %%        LOCALIZATION PROCEDURE NOT CONVERGED!!!     %% '
        write(CFG%lupri,*) ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(CFG%lupri,*) ''
        write(CFG%lupri,*) ' Localization is not converged in the maximum number of     '
        write(CFG%lupri,*) ' iterations. Restart calculation by renaming orbital file'
        write(CFG%lupri,*) ' localized_orbitals.u to orbitals_in.u, and run calculation '
        write(CFG%lupri,*) ' using keyword .ONLY LOC in **LOCALIZE ORBITALS section'
        write(CFG%lupri,*) ' or restart calculation from scratch and increase number of '
        write(CFG%lupri,*) ' macro iterations allowed using keyword .MACRO IT as described '
        write(CFG%lupri,*) ' in user manual.'
        write(CFG%lupri,*) '  '
        write(CFG%lupri,*) '  '
    endif   


  call mem_dealloc(max_orbspreads)
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
integer :: norb, i,imx,idamax,iter_number
real(realk) :: nrmG, oVal,old_oVal
real(realk) :: nrm_thresh,stepsize,orig_Eval
real(realk),pointer :: max_orbspreads(:)  


  norb=CMO%ncol
  CFG%orbspread_input=>orbspread_input
  call mem_alloc(max_orbspreads,CFG%max_macroit)
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
  stepsize=0d0
  CFG%it = 1
  do i=1,CFG%max_macroit
    iter_number= i
    CFG%old_mu = CFG%mu
    old_oVal = oVal
    imx  =  idamax(norb,orbspread_input%spread2,1)
    nrmG = sqrt(mat_sqnorm2(G))/real(norb)
    max_orbspreads(i)=sqrt(orbspread_input%spread2(imx))
    write (ls%lupri,'(I3,A,ES8.1,A,f6.2,A,ES8.1,A,ES8.1,A,I2,A,f5.2,A,f5.2)') &
         &i, ' Pred=',CFG%r_denom,' sigma_2 =',sqrt(orbspread_input%spread2(imx)),&
         &  ' mu = ',CFG%mu,' grd = ', nrmG, ' it = ',CFG%it, ' trust-region = ',CFG%stepsize,' step =', stepsize
    
  if (i>10) then
     if ( abs(max_orbspreads(i)-max_orbspreads(i-10)) < 0.05 ) then
        write(ls%lupri,*) '  '
        write(ls%lupri,*) '    ********* Orbital localization converged ************'
        write(ls%lupri,*) '    *                                                   *'
        write(ls%lupri,*) '    * There are insignificant changes in the  locality  *'
        write(ls%lupri,*) '    * of the least local orbital, and procedure is      *'
        write(ls%lupri,*) '    * exited irresptective of gradient norm.            *'
        write(ls%lupri,*) '    *                                                   *'
        write(ls%lupri,*) '    *****************************************************'
        write(ls%lupri,*) '  '
        exit
     endif
  endif
  if( nrmG.le. CFG%macro_thresh) then
        write(ls%lupri,*) '  '
        write(ls%lupri,*) '   ********* Orbital localization converged ************'
        write(ls%lupri,*) '   *                                                   *'
        write(ls%lupri,*) '   * The gradient norm for the orbital localization    *'
        write(ls%lupri,*) '   * function is below the threshold, and we exit      *'
        write(ls%lupri,*) '   * the localization procedure.                       *'
        write(ls%lupri,*) '   *                                                   *'
        write(ls%lupri,*) '   *****************************************************'
        write(ls%lupri,*) '  '
        exit
    end if

   
   call davidson_solver(CFG,G,X)

   ! global and local thresholds defined in CFG settings
   if (dabs(CFG%mu)> 1.0) CFG%conv_thresh=CFG%global_conv_thresh
   if (dabs(CFG%mu)< 1.0)  CFG%conv_thresh=CFG%local_conv_thresh
print*, "trhesh", CFG%local_conv_thresh
    call mat_copy(1.0_realk,CMO,CMOsav)
 
    stepsize = CFG%stepsize
    call linesearch_orbspread2(CFG,cmo,X,stepsize,old_oval,orig_Eval)
    call orbspread_value(oVal,orbspread_input)
    

    if (orig_Eval-old_oVal > 0) then
       write(ls%lupri,*) 'Step not accepted. Go back'
       call mat_copy(1.0d0,CMOsav,CMO)
       call orbspread_update(orbspread_input,CMO)
       call orbspread_value(oVal,orbspread_input)
       CFG%Stepsize = CFG%Stepsize/2.0d0
    else
      CFG%Stepsize = min(CFG%Stepsize*2.5d0,CFG%max_stepsize) 
    endif
     
   if (CFG%stepsize < 0.001) then
           write(CFG%lupri,*) ''
           write(CFG%lupri,'(a)') 'WARNING: Stepsize too small. ' 
           if (i>5) then
              if  (abs(max_orbspreads(i)-max_orbspreads(i-5))< 0.1) then 
                 write(CFG%lupri,*) ' However, the locality of the least local orbital       ' 
                 write(CFG%lupri,*) ' has not changed significantly the last five iterations ' 
                 write(CFG%lupri,*) ' and the generated orbitals are localized, and will      ' 
                 write(CFG%lupri,*) ' be written to file.   '
                 write(CFG%lupri,*) ''
	         exit
               endif
           endif
           write(CFG%lupri,*) ' Cannot proceed with localization due to issues with    ' 
           write(CFG%lupri,*) ' solving the level-shifted Newton equations. You may    ' 
           write(CFG%lupri,*) ' try to restart calculation and lower the residual norm ' 
           write(CFG%lupri,*) ' threshold for the micro iterations as described in     '
           write(CFG%lupri,*) ' the user manual under section **LOCALIZE ORBITALS      '
           write(CFG%lupri,*) ' and keyword .LOOSE MICRO THRESH                              ' 
           call lsquit('Cannot converge micro iterations. ', CFG%lupri)
   elseif (orig_Eval-old_oVal > 0) then
            cycle
   endif
    !new gradient
    call orbspread_gradx(G,norb,orbspread_input)


   call orbspread_precond_matrix2(orbspread_input,P,norb)
   orbspread_input%P => P
   CFG%P =>orbspread_input%P


  enddo
    if (iter_number==CFG%max_macroit) then
        write(CFG%lupri,*) ''
        write(CFG%lupri,*) ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(CFG%lupri,*) ' %%        LOCALIZATION PROCEDURE NOT CONVERGED!!!     %% '
        write(CFG%lupri,*) ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(CFG%lupri,*) ''
        write(CFG%lupri,*) ' Localization is not converged in the maximum number of     '
        write(CFG%lupri,*) ' iterations. Restart calculation by renaming orbital file'
        write(CFG%lupri,*) ' localized_orbitals.u to orbitals_in.u, and run calculation '
        write(CFG%lupri,*) ' using keyword .ONLY LOC in **LOCALIZE ORBITALS section'
        write(CFG%lupri,*) ' or restart calculation from scratch and increase number of '
        write(CFG%lupri,*) ' macro iterations allowed using keyword .MACRO IT as described '
        write(CFG%lupri,*) ' in user manual.'
        write(CFG%lupri,*) ''
    endif   

  call mem_dealloc(max_orbspreads)
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
