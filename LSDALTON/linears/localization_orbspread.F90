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
   write(ls%lupri,'(A,I3,A,f6.2,A,ES8.1,A,ES8.1,A,I3,A,f5.2,A,f5.2)') &
     &'  %LOC% ',i,' sigma_4 =',max_FM,&
     &   '  mu = ',CFG%mu,'  grd =', nrmG, '  it =',CFG%it, '  trust-region = ', CFG%stepsize&
     & ,' step= ',stepsize
    

    if (i>10) then
      if ( abs(max_orbspreads(i)-max_orbspreads(i-10)) < 0.02 .and. &
      &   abs(max_orbspreads(i)-max_orbspreads(i-9)) < 0.02  .and. &
      &   abs(max_orbspreads(i)-max_orbspreads(i-8)) < 0.02  .and. &
      &   abs(max_orbspreads(i)-max_orbspreads(i-7)) < 0.02  .and. &
      &   abs(max_orbspreads(i)-max_orbspreads(i-6)) < 0.02  .and. &
      &   abs(max_orbspreads(i)-max_orbspreads(i-5)) < 0.02) then
           write(ls%lupri,'(a)') '  %LOC%  '
           write(ls%lupri,'(a)') '  %LOC%   ********* Orbital localization converged ************'
           write(ls%lupri,'(a)') '  %LOC%   *                                                   *'
           write(ls%lupri,'(a)') '  %LOC%   * There are insignificant changes in the  locality  *'
           write(ls%lupri,'(a)') '  %LOC%   * of the least local orbital, and procedure is      *'
           write(ls%lupri,'(a)') '  %LOC%   * exited irresptective of gradient norm.            *'
           write(ls%lupri,'(a)') '  %LOC%   *                                                   *'
           write(ls%lupri,'(a)') '  %LOC%   *****************************************************'
           write(ls%lupri,'(a)') '  %LOC% '
           exit
        endif
    endif
    if( nrmG.le. CFG%macro_thresh*10.0) then
        write(ls%lupri,'(a)') '  %LOC% '
        write(ls%lupri,'(a)') '  %LOC%  ********* Orbital localization converged ************'
        write(ls%lupri,'(a)') '  %LOC%  *                                                   *'
        write(ls%lupri,'(a)') '  %LOC%  * The gradient norm for the orbital localization    *'
        write(ls%lupri,'(a)') '  %LOC%  * function is below the threshold, and we exit      *'
        write(ls%lupri,'(a)') '  %LOC%  * the localization procedure.                       *'
        write(ls%lupri,'(a)') '  %LOC%  *                                                   *'
        write(ls%lupri,'(a)') '  %LOC%  *****************************************************'
        write(ls%lupri,'(a)') '  %LOC% '
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
       write(CFG%lupri,'(a)') "  %LOC% Step accepted"
       CFG%stepsize=min(CFG%stepsize*2.0,CFG%max_stepsize)
       !CMOS are updated in linesearch
    else
       write(CFG%lupri,'(a)') "  %LOC% Step rejected"
       call mat_copy(1d0,cmosav,cmo)
       call kurt_updateAO(CFG%PFM_input,cmo)
       call kurt_value(CFG%PFM_input)
       old_oVal=CFG%PFM_input%kurt_val
       oVal=CFG%PFM_input%kurt_val
       CFG%stepsize=CFG%stepsize*0.5
    endif    

   if (CFG%stepsize < 0.001) then
        write(CFG%lupri,'(a)') '  %LOC%'
        write(CFG%lupri,'(a)') '  %LOC% WARNING: Stepsize too small. ' 
        if (i>5) then
           if  (abs(max_orbspreads(i)-max_orbspreads(i-5))< 0.1) then 
               write(CFG%lupri,'(a)') '  %LOC% However, the locality of the least local orbital       ' 
               write(CFG%lupri,'(a)') '  %LOC% has not changed significantly the last five iterations ' 
               write(CFG%lupri,'(a)') '  %LOC% and the generated orbitals are localized, and will      ' 
               write(CFG%lupri,'(a)') '  %LOC% be written to file.   '
               write(CFG%lupri,'(a)') '  %LOC%'
	       exit
           endif
        endif
        write(CFG%lupri,'(a)') '  %LOC% Cannot proceed with localization due to issues with    ' 
        write(CFG%lupri,'(a)') '  %LOC% solving the level-shifted Newton equations. You may    ' 
        write(CFG%lupri,'(a)') '  %LOC% try to restart calculation and lower the residual norm ' 
        write(CFG%lupri,'(a)') '  %LOC% threshold for the micro iterations as described in     '
        write(CFG%lupri,'(a)') '  %LOC% the user manual under section **LOCALIZE ORBITALS      '
        write(CFG%lupri,'(a)') '  %LOC% and keyword .LOOSE MICRO THRESH                              ' 
        call lsquit(' %LOC% Cannot converge micro iterations. ', CFG%lupri)
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
        write(CFG%lupri,'(a)') '  %LOC% '
        write(CFG%lupri,'(a)') '  %LOC%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(CFG%lupri,'(a)') '  %LOC%  %%        LOCALIZATION PROCEDURE NOT CONVERGED!!!     %% '
        write(CFG%lupri,'(a)') '  %LOC%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(CFG%lupri,'(a)') '  %LOC% '
        write(CFG%lupri,'(a)') '  %LOC%  Localization is not converged in the maximum number of     '
        write(CFG%lupri,'(a)') '  %LOC%  iterations. Restart calculation by renaming orbital file'
        write(CFG%lupri,'(a)') '  %LOC%  localized_orbitals.u to orbitals_in.u, and run calculation '
        write(CFG%lupri,'(a)') '  %LOC%  using keyword .ONLY LOC in **LOCALIZE ORBITALS section'
        write(CFG%lupri,'(a)') '  %LOC%  or restart calculation from scratch and increase number of '
        write(CFG%lupri,'(a)') '  %LOC%  macro iterations allowed using keyword .MACRO IT as described '
        write(CFG%lupri,'(a)') '  %LOC%  in user manual.'
        write(CFG%lupri,'(a)') '  %LOC%   '
        write(CFG%lupri,'(a)') '  '
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
subroutine orbspread_localize_davidson(CFG,CMO,m,ls)
implicit none
type(RedSpaceItem)           :: CFG
type(Matrix) , intent(inout ):: CMO
TYPE(lsitem) , intent(inout) :: ls
integer      , intent(in)    :: m
type(Matrix) :: CMOsav
type(Matrix), target  ::  X, P, G
integer :: norb, i,imx,idamax,iter_number
real(realk) :: nrmG, oVal,old_oVal
real(realk) :: nrm_thresh,stepsize,orig_Eval
real(realk),pointer :: max_orbspreads(:)  


  norb=CMO%ncol
  call mem_alloc(max_orbspreads,CFG%max_macroit)
  call mat_init(X,norb,norb)
  call mat_init(G,norb,norb)
  call mat_init(P,norb,norb)
  call mat_init(CMOsav,CMO%nrow,CMO%ncol)


  call orbspread_init(CFG%orbspread_inp,m,norb)
  call orbspread_update(CFG%orbspread_inp,CMO)
  call orbspread_gradx(G,norb,CFG%orbspread_inp)
  call orbspread_value(oVal,CFG%orbspread_inp)


  call orbspread_precond_matrix2(CFG%orbspread_inp,P,norb)
  CFG%orbspread_inp%P => P
  CFG%P=> P

  stepsize = CFG%stepsize
  CFG%mu = 0.0_realk
  stepsize=0.0_realk
  CFG%it = 1
  do i=1,CFG%max_macroit
    iter_number= i
    CFG%old_mu = CFG%mu
    old_oVal = oVal
    imx  =  idamax(norb,CFG%orbspread_inp%spread2,1)
    nrmG = sqrt(mat_sqnorm2(G))/real(norb)
    max_orbspreads(i)=sqrt(CFG%orbspread_inp%spread2(imx))
    write (ls%lupri,'(A,I3,A,f6.2,A,ES8.1,A,ES8.1,A,I2,A,f5.2,A,f5.2)') &
         &'  %LOC%',i,' sigma_2 =',max_orbspreads(i),&
         &  ' mu = ',CFG%mu,' grd = ', nrmG, ' it = ',CFG%it, ' trust-region = ',CFG%stepsize,' step =', stepsize
  

  if (i>10) then
        if ( abs(max_orbspreads(i)-max_orbspreads(i-10)) < 0.01 .and. &
         &   abs(max_orbspreads(i)-max_orbspreads(i-9)) < 0.01  .and. &
         &   abs(max_orbspreads(i)-max_orbspreads(i-8)) < 0.01  .and. &
         &   abs(max_orbspreads(i)-max_orbspreads(i-7)) < 0.01  .and. &
         &   abs(max_orbspreads(i)-max_orbspreads(i-6)) < 0.01  .and. &
         &   abs(max_orbspreads(i)-max_orbspreads(i-5)) < 0.01) then
            write(ls%lupri,'(a)') '  %LOC% '
            write(ls%lupri,'(a)') '  %LOC%   ********* Orbital localization converged ************'
            write(ls%lupri,'(a)') '  %LOC%   *                                                   *'
            write(ls%lupri,'(a)') '  %LOC%   * There are insignificant changes in the  locality  *'
            write(ls%lupri,'(a)') '  %LOC%   * of the least local orbital, and procedure is      *'
            write(ls%lupri,'(a)') '  %LOC%   * exited irresptective of gradient norm.            *'
            write(ls%lupri,'(a)') '  %LOC%   *                                                   *'
            write(ls%lupri,'(a)') '  %LOC%   *****************************************************'
            write(ls%lupri,'(a)') '  %LOC% '
            exit
     endif
  endif
  if( nrmG.le. CFG%macro_thresh) then
        write(ls%lupri,'(a)') '  %LOC% '
        write(ls%lupri,'(a)') '  %LOC%  ********* Orbital localization converged ************'
        write(ls%lupri,'(a)') '  %LOC%  *                                                   *'
        write(ls%lupri,'(a)') '  %LOC%  * The gradient norm for the orbital localization    *'
        write(ls%lupri,'(a)') '  %LOC%  * function is below the threshold, and we exit      *'
        write(ls%lupri,'(a)') '  %LOC%  * the localization procedure.                       *'
        write(ls%lupri,'(a)') '  %LOC%  *                                                   *'
        write(ls%lupri,'(a)') '  %LOC%  *****************************************************'
        write(ls%lupri,'(a)') '  %LOC% '
        exit
    end if

   
   call davidson_solver(CFG,G,X)

   ! global and local thresholds defined in CFG settings
   if (dabs(CFG%mu)> 1.0_realk) CFG%conv_thresh=CFG%global_conv_thresh
   if (dabs(CFG%mu)< 1.0_realk)  CFG%conv_thresh=CFG%local_conv_thresh
    call mat_copy(1.0_realk,CMO,CMOsav)
 
    stepsize = CFG%stepsize
    call linesearch_orbspread(CFG,cmo,X,stepsize,old_oval,orig_Eval,nrmG,i)
    call orbspread_value(oVal,CFG%orbspread_inp)

    if (orig_Eval-old_oVal > 0.0_realk) then
       write(ls%lupri,'(a)') '  %LOC% Step not accepted. Go back'
       call mat_copy(1.0d0,CMOsav,CMO)
       call orbspread_update(CFG%orbspread_inp,CMO)
       call orbspread_value(oVal,CFG%orbspread_inp)
       CFG%Stepsize = CFG%Stepsize/2.0_realk
    else
      CFG%Stepsize = min(CFG%Stepsize*2.5_realk,CFG%max_stepsize) 
    endif
     
   if (CFG%stepsize < 0.001_realk) then
           write(CFG%lupri,'(a)') '  %LOC%'
           write(CFG%lupri,'(a)') '  %LOC% WARNING: Stepsize too small. ' 
           if (i>5) then
              if  (abs(max_orbspreads(i)-max_orbspreads(i-5))< 0.1_realk) then 
                 write(CFG%lupri,'(a)') '  %LOC% However, the locality of the least local orbital       ' 
                 write(CFG%lupri,'(a)') '  %LOC% has not changed significantly the last five iterations ' 
                 write(CFG%lupri,'(a)') '  %LOC% and the generated orbitals are localized, and will      ' 
                 write(CFG%lupri,'(a)') '  %LOC% be written to file.   '
                 write(CFG%lupri,'(a)') '  %LOC%'
	         exit
               endif
           endif
           write(CFG%lupri,'(a)') '  %LOC% Cannot proceed with localization due to issues with    ' 
           write(CFG%lupri,'(a)') '  %LOC% solving the level-shifted Newton equations. You may    ' 
           write(CFG%lupri,'(a)') '  %LOC% try to restart calculation and lower the residual norm ' 
           write(CFG%lupri,'(a)') '  %LOC% threshold for the micro iterations as described in     '
           write(CFG%lupri,'(a)') '  %LOC% the user manual under section **LOCALIZE ORBITALS      '
           write(CFG%lupri,'(a)') '  %LOC% and keyword .LOOSE MICRO THRESH                              ' 
           call lsquit('%LOC% Cannot converge micro iterations. ', CFG%lupri)
   elseif (orig_Eval-old_oVal > 0.0_realk) then
            cycle
   endif
    !new gradient
   call orbspread_gradx(G,norb,CFG%orbspread_inp)
   call orbspread_precond_matrix2(CFG%orbspread_inp,P,norb)

  enddo
    if (iter_number==CFG%max_macroit) then
        write(CFG%lupri,'(a)') '  %LOC% '
        write(CFG%lupri,'(a)') '  %LOC%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(CFG%lupri,'(a)') '  %LOC%  %%        LOCALIZATION PROCEDURE NOT CONVERGED!!!     %% '
        write(CFG%lupri,'(a)') '  %LOC%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(CFG%lupri,'(a)') '  %LOC% '
        write(CFG%lupri,'(a)') '  %LOC%  Localization is not converged in the maximum number of     '
        write(CFG%lupri,'(a)') '  %LOC%  iterations. Restart calculation by renaming orbital file'
        write(CFG%lupri,'(a)') '  %LOC%  localized_orbitals.u to orbitals_in.u, and run calculation '
        write(CFG%lupri,'(a)') '  %LOC%  using keyword .ONLY LOC in **LOCALIZE ORBITALS section'
        write(CFG%lupri,'(a)') '  %LOC%  or restart calculation from scratch and increase number of '
        write(CFG%lupri,'(a)') '  %LOC%  macro iterations allowed using keyword .MACRO IT as described '
        write(CFG%lupri,'(a)') '  %LOC%  in user manual.'
        write(CFG%lupri,'(a)') '  %LOC% '
    endif   

  call mem_dealloc(max_orbspreads)
  call orbspread_free(CFG%orbspread_inp)
  call mat_free(X)
  call mat_free(G)
  call mat_free(P)
  call mat_free(CMOsav)

end subroutine orbspread_localize_davidson


subroutine linesearch_kurtosis(CFG,cmo,X,stepsize,oval)
implicit none
type(RedSpaceItem) :: CFG
type(matrix)  :: cmo,X
integer :: i,numb,nmats
type(matrix)  :: cmotemp(4),Xtemp(4)
real(realk) :: old_funcval,factor,step(4),stepsize,oval

numb=4
step(1) = 1.0d0
step(2)= 2.0d0
step(3)= 4.0d0
step(4)= 8.0d0
    old_funcval = CFG%PFM_input%kurt_val
    
    if (CFG%orb_debug) write(CFG%lupri,'(a,I4,a,ES13.3)') &
    &'Linesearch number :', 0, ' Original function value: ', CFG%PFM_input%kurt_val
do i=1,numb
    call mat_init(Xtemp(i),X%nrow,X%ncol)
    call mat_assign(Xtemp(i),X)
    call mat_init(cmotemp(i),cmo%nrow,cmo%ncol)
    call mat_assign(cmotemp(i),cmo)
    call mat_scal(step(i),Xtemp(i))
    call updatecmo(CMOtemp(i),Xtemp(i))      
    call kurt_updateAO(CFG%PFM_input,CMOtemp(i))
    call kurt_value(CFG%PFM_input)
    if (CFG%orb_debug) write(CFG%lupri,'(a,I4,a,ES13.3)') &
    &'Linesearch number :', i, ' Change ', CFG%PFM_input%kurt_val-old_funcval
    if (i==1 .and. CFG%PFM_input%kurt_val> old_funcVal) then
      oVal = CFG%PFM_input%kurt_val
      nmats=i 
      exit
    end if
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



subroutine linesearch_orbspread(CFG,cmo,X,stepsize,value_last_macro,orig_eival,nrmg,macroit)
implicit none
type(RedSpaceItem) :: CFG
type(matrix)  :: cmo,X
real(realk),intent(in) :: value_last_macro,nrmg
integer :: i,numb=5,nmats,macroit
type(matrix)  :: cmotemp(5),Xtemp(5)
real(realk) :: old_funcval,factor(6),step(5),stepsize,oval,d(6)
real(realk) :: orig_eival

   old_funcval = value_last_macro
   if (CFG%orb_debug) write(CFG%lupri,'(a,I4,a,f15.1)') &
   &'Linesearch number :', 0, ' Original function value: ', old_funcval

   factor(2)=1.0_realk
   factor(3)=2.0_realk
   factor(4)=4.0_realk
   factor(5)=8.0_realk
   do i=2,5
       call mat_init(Xtemp(i),X%nrow,X%ncol)
       call mat_copy(factor(i),X,Xtemp(i))
       call mat_init(cmotemp(i),cmo%nrow,cmo%ncol)
       call mat_copy(1.0d0,cmo,cmotemp(i))
       call updatecmo(CMOtemp(i),Xtemp(i))
       call orbspread_update(CFG%orbspread_inp,CMOtemp(i))
       call orbspread_value(oVal,CFG%orbspread_inp)
       d(i)=oVal
       if (CFG%orb_debug) write(CFG%lupri,'(a,I4,a,f15.4,a,f7.2)') &
       &'Linesearch number :', i, ' Change ', d(i)-d(i-1), '  factor  ', factor(i)
       if (i==2 .and. oVal > old_funcVal) then
          orig_eival = oVal
          nmats=i
          exit
       endif
       if (oVal > old_funcVal) then
              call mat_assign(cmo,cmotemp(i-1))
              call orbspread_update(CFG%orbspread_inp,CMO)
              call orbspread_value(oVal,CFG%orbspread_inp)
              stepsize = dsqrt(mat_dotproduct(xtemp(i-1),xtemp(i-1)))
              nmats=i
              orig_eival = old_funcVal
              exit
       end if
       if (i==5 .or. dabs(oVal-old_funcval)< 1.0) then
         call mat_assign(cmo,cmotemp(i))
         call orbspread_update(CFG%orbspread_inp,CMO)
         call orbspread_value(oVal,CFG%orbspread_inp)
         stepsize = dsqrt(mat_dotproduct(xtemp(i),xtemp(i)))
         nmats=i
         orig_eival= oVal
         exit
       end if
       old_funcval=oVal
   end do
    do i=2,nmats
       call mat_free(cmotemp(i))
       call mat_free(xtemp(i))
    enddo


end subroutine linesearch_orbspread






 end module orbspread_module
