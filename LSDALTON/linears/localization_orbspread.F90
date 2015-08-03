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
use matrix_util, only: matrix_exponential

contains

!> \brief Routine that drives macro iterations for localizing orbitals using FM
!> \author Ida-Marie Hoeyvik
  subroutine PFM_localize_davidson(CFG,CMOall,m,ls,norb)
    implicit none
    type(RedSpaceItem)           :: CFG
    type(Matrix) , intent(inout ):: CMOall
    TYPE(lsitem) , intent(inout) :: ls
    integer      , intent(in)    :: m,norb
    logical :: lower2
    type(Matrix) :: Xsav,CMOsav,CMO
    type(Matrix), target  ::  X, P, G,expX
    integer :: i, nbas,iter_number
    real(realk) :: nrmG, oVal,old_oVal, max_step,max_FM
    real(realk) :: nrm_thresh,stepsize
    real(realk) :: trial(1,1)
    real(realk),pointer :: max_orbspreads(:)  
    real(realk), pointer :: tmp(:)
    integer :: lun, counter
    logical :: OnMaster

    !initializations 
    OnMaster=.true.
    nbas=CMOall%nrow
    counter=0
    
    ! Extract coefficients to be localized
    call mem_alloc(tmp,nbas*norb)
    call mat_init(CMO,nbas,norb)
    ! extract matrix from CMOall(1,offset)
    call mat_retrieve_block(CMOall,tmp,nbas,norb,1,CFG%offset)
    call mat_set_from_full(tmp,1.0_realk,CMO)
    call mem_dealloc(tmp)


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

#ifdef COMPILER_UNDERSTANDS_FORTRAN_2003
       FLUSH(ls%lupri)
#endif

       if( nrmG.le. CFG%macro_thresh .and. i.gt.1) then
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
          write(CFG%lupri,'(a)') "   Step accepted"
          CFG%stepsize=min(CFG%stepsize*2.0,CFG%max_stepsize)
          !CMOS are updated in linesearch
       else
          write(CFG%lupri,'(a)') "   Step rejected"
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
       elseif (oVal-old_oVal > 0.0) then
          cycle
       endif

       !new gradient
       call compute_gradient(CFG%PFM_input,G,norb)

       !new preconditioning matrix 
       call kurtosis_precond_matrix(CFG%PFM_input,P)
       CFG%PFM_input%P => P
       CFG%P => CFG%PFM_input%P

       !make restart file
       counter = counter + 1
       if (counter == CFG%orbital_save_interval) then
          call mem_alloc(tmp,nbas*norb)
          call mat_to_full(CMO,1.0_realk,tmp)
          call mat_create_block(CMOall,tmp,nbas,norb,1,CFG%offset)
          call mem_dealloc(tmp)
          lun = -1
          call lsopen(lun,'localized_orbitals.restart','unknown','UNFORMATTED')
          call mat_write_to_disk(lun,CMOall,OnMaster)
          call LSclose(LUN,'KEEP')
          write(CFG%lupri,'(a)') '  %LOC% temporary orbitals written to localized_orbitals.restart'
          counter = 0
       endif


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

    
    !Put localized block into full CMO matrix
    call mem_alloc(tmp,nbas*norb)
    call mat_to_full(CMO,1.0_realk,tmp)
    call mat_free(CMO)
    call mat_create_block(CMOall,tmp,nbas,norb,1,CFG%offset)
    call mem_dealloc(tmp)


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
  subroutine orbspread_localize_davidson(CFG,CMOall,m,ls,norb)
    implicit none
    type(RedSpaceItem)           :: CFG
    type(Matrix) , intent(inout ):: CMOall
    TYPE(lsitem) , intent(inout) :: ls
    integer      , intent(in)    :: m,norb
    type(Matrix) :: CMOsav, CMO
    type(Matrix), target  ::  X, P, G
    integer ::  i,imx,idamax,iter_number
    real(realk) :: nrmG, oVal,old_oVal
    real(realk) :: nrm_thresh,stepsize
    real(realk),pointer :: max_orbspreads(:)  
    real(realk),pointer :: tmp(:)
    integer :: lun, counter,nbas
    logical :: onmaster

        !initializations 
    OnMaster=.true.
    nbas=CMOall%nrow
    counter=0

    ! Extract coefficients to be localized
    call mem_alloc(tmp,nbas*norb)
    call mat_init(CMO,nbas,norb)
    ! extract matrix from CMOall(1,offset)
    call mat_retrieve_block(CMOall,tmp,nbas,norb,1,CFG%offset)
    call mat_set_from_full(tmp,1.0_realk,CMO)
    call mem_dealloc(tmp)


    call mem_alloc(max_orbspreads,CFG%max_macroit)
    call mat_init(X,norb,norb)
    call mat_init(G,norb,norb)
    call mat_init(P,norb,norb)
    call mat_init(CMOsav,CMO%nrow,CMO%ncol)


    call orbspread_init(CFG%orbspread_inp,m,norb)
    call orbspread_update(CFG%orbspread_inp,CMO)
    call orbspread_gradx(G,norb,CFG%orbspread_inp)
    !check for symmetry minimum, make random rotation
    ! in case such a minimum is detected
    nrmG = sqrt(mat_sqnorm2(G))
    if (nrmG.lt. 1.0E-8_realk) then
       write(CFG%lupri,'(a,ES10.2)') '  %LOC% False minimum detected. Gradient norm = ', nrmG 
       call mat_assign(X,G)
       call normalize(X) 
       call mat_trans(X,G)
       call mat_daxpy(-1.0_realk,G,X)
       call normalize(X)
       call mat_scal(0.001_realk,X)
       call updatecmo(cmo,X) 
       call orbspread_update(CFG%orbspread_inp,CMO)
       call orbspread_gradx(G,norb,CFG%orbspread_inp)
       write(CFG%lupri,'(a,ES10.2)') '  %LOC% Norm after random rotation = ', sqrt(mat_sqnorm2(G)) 
    endif

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

#ifdef COMPILER_UNDERSTANDS_FORTRAN_2003
       FLUSH(ls%lupri)
#endif


       if( nrmG.le. CFG%macro_thresh .and. i.gt.1) then
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


       call linesearch_orbspread(CFG,cmo,X,stepsize,oVal)
       call orbspread_value(oVal,CFG%orbspread_inp)


       if (oVal-old_oVal > 0.0_realk) then
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
       elseif (oval-old_oVal > 0.0_realk) then
          cycle
       endif
       !new gradient
       call orbspread_gradx(G,norb,CFG%orbspread_inp)
       call orbspread_precond_matrix2(CFG%orbspread_inp,P,norb)
    
       !make restart file
       counter = counter + 1
       if (counter == CFG%orbital_save_interval) then
          call mem_alloc(tmp,nbas*norb)
          call mat_to_full(CMO,1.0_realk,tmp)
          call mat_create_block(CMOall,tmp,nbas,norb,1,CFG%offset)
          call mem_dealloc(tmp)
          lun = -1
          call lsopen(lun,'localized_orbitals.restart','unknown','UNFORMATTED')
          call mat_write_to_disk(lun,CMOall,OnMaster)
          call LSclose(LUN,'KEEP')
          counter = 0
          write(CFG%lupri,'(a)') '  %LOC% temporary orbitals written to localized_orbitals.restart'
       endif


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
   
    !Put localized block into full CMO matrix
    call mem_alloc(tmp,nbas*norb)
    call mat_to_full(CMO,1.0_realk,tmp)
    call mat_free(CMO)
    call mat_create_block(CMOall,tmp,nbas,norb,1,CFG%offset)
    call mem_dealloc(tmp)


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




  subroutine linesearch_orbspread(CFG,cmo,X,stepsize,fval)
    implicit none
    type(RedSpaceItem) :: CFG
    type(matrix),intent(inout)  :: cmo
    type(matrix) :: X
    type(matrix) :: expX,scr,cmosave
    !> fval  function value before taking step X 
    real(realk),intent(inout)  :: fval
    real(realk)                :: old_fval,new_fval
    real(realk),intent(inout)  :: stepsize
    real(realk)                :: normX
    integer :: i, max_iter

    max_iter = 5
    normX = sqrt(mat_sqnorm2(X))
    old_fval = fval
    write(CFG%lupri,'(a,I4,a,f15.1)') &
         &'Linesearch  :', 0, ' Original function value: ', old_fval

    call mat_init(expX,X%nrow,X%ncol);    call mat_zero(expX)  
    call mat_init(scr,X%nrow,X%ncol);     call mat_zero(scr)  
    call mat_init(cmosave,CMO%nrow,CMO%ncol)
    call matrix_exponential(X,expX,1E-12_realk)
    call mat_assign(scr,expX)

    !For large systems, increment of factor 2 
    if (X%ncol > 600) then
       call mat_assign(X,expX)
       call mat_mul(X,expX,'n','n',1.0_realk,0.0_realk,scr)
       call mat_assign(expX,scr)
    endif
    do i=1,max_iter
       stepsize = dble(i)*normX 
       if (X%ncol > 600) stepsize = 2.0_realk*stepsize 
       call mat_assign(cmosave,cmo)
       call mat_mul(CMOsave,scr,'n','n',1.0_realk,0.0_realk,cmo)
       call orbspread_update(CFG%orbspread_inp,cmo)  
       call orbspread_value(new_fval,CFG%orbspread_inp) 
       write(CFG%lupri,'(a,I4,a,f15.1)') &
            &'Linesearch  :', i, ' Diff. in function value: ', new_fval-old_fval
       if (new_fval > old_fval) then
          write(CFG%lupri,'(a)') 'Linesearch done, we choose previous point '
          call mat_assign(cmo,cmosave)
          call orbspread_update(CFG%orbspread_inp,cmosave) 
          fval = old_fval 
          exit
       endif
       if (i< max_iter) then
          call mat_assign(X,scr) ! use X as scratch for exp(i*X)
          call mat_mul(X,expX,'n','n',1.0_realk,0.0_realk,scr) 
          old_fval = new_fval 
       elseif (i == max_iter) then
          write(CFG%lupri,'(a)')  ' Linesearch max iter reached, we choose last point'
       endif
    enddo

    call mat_free(scr)
    call mat_free(expX)
    call mat_free(cmosave)


  end subroutine linesearch_orbspread





 end module orbspread_module
