!> @file
!> Contains level shift module.

!> \brief Dynamic determination of level shift used by Augmented Roothaan-Hall and Direct Density optimization.
!> \author S. Host
!> \date 2007
module levelshift_mod
use queue_module
use queue_ops
use precision
use memory_handling
use matrix_operations
use matrix_module
!> \brief Contains setting for level shift.
!> \author S. Host
!> \date March 2010
type lshiftItem
   !> Logical unit number for output file
   integer     :: lupri
   !> True if level shift should have a fixed value, i.e. nogen dynamic level shifting
   logical     :: fixed_shift      
   !> If fixed_shift=true, use this value for level shifting
   real(realk) :: fixed_shift_param
   !> True is CROP solver is used for linear equations
   logical     :: arh_crop
   !> True is truncated CROP scheme is used for linear equations
   logical     :: arh_truncate
   !> New damping scheme better suited for truncated scheme. NOT FULLY TESTED!
   logical     :: arh_newdamp
   !> Do not allow level shift to go below this value (0.0E0_realk => standard scheme)
   real(realk) :: min_lshift
   !> If true, step size restriction is on max element instead of Frobenius norm
   logical     :: optxelm
   !> Maximum matrix element size allowed (i.e. trust radius)
   real(realk) :: max_element
   !> Maximum Frobenius norm allowed (i.e. trust radius)
   real(realk) :: max_step
   !> Level shift from previous SCF/localization iteration
   real(realk) :: current_mu
   !> Level shift by calculating HOMO-LUMO gap in each SCF iteration
   logical     :: lshift_by_hlgap
   !INFO/DEBUG VARIABLES
   !--------------------
   !> Print detailed info from level shift
   logical     :: info_levelshift
end type lshiftItem 

contains

   !> \brief Wrapper for level shift routine.
   !> \author S. Host
   !> \date 2007
   subroutine main_levelshift(lshift,Ared,Sred,Gred,redsize,lub,red_space_dim,rowdim,coldim,mu,xF,sigmaF,vectorsubspace)
   !FIXME: Collect all these bloody arguments in a structure!!!
   implicit none
        !> Contains setting for level shift.
        type(lshiftItem),intent(in)   :: lshift
        !> Max allowed size of reduced space (allocated size)
        integer, intent(in)           :: redsize
        !> Reduced Hessian
        real(realk),intent(in)        :: Ared(redsize,redsize)
        !> Reduced overlap of trial vectors
        real(realk),intent(in)        :: Sred(redsize,redsize)
        !> Reduced gradient
        real(realk),intent(in)        :: Gred(redsize)
        !> Logical unit number for file containing previous trial vectors (if arh_truncate=false)
        integer, intent(in)           :: lub
        !> Current size of reduced space
        integer, intent(in)           :: red_space_dim
        !> Number of basis functions
        integer, intent(in)           :: rowdim
        !> Number of basis functions (or 1 if vector instead of matrix)
        integer, intent(in)           :: coldim
        !> Level shift
        real(realk), intent(inout)    :: mu
        !> Matrix used for new level shift scheme. NOT FULLY TESTED
        type(matrix), intent(in),optional :: xF
        !> Matrix used for new level shift scheme. NOT FULLY TESTED
        type(matrix), intent(in),optional :: sigmaF
        !> Contains subspace of previous trial vectors (if arh_truncate=true)
        TYPE(modFIFO), intent(inout),optional :: vectorsubspace
        type(matrix)                  :: res, x, x5(5) 
        type(matrix)                  :: scrmat
        real(realk)                   :: xnorm, maxelm, dummy
        real(realk)                   :: best_mu_yet, &
                                       & b_coeff(5), c_coeff(5), d_coeff(5)
        integer                       :: i, maxit, fullsize
        logical                       :: extradim
        real(realk), pointer          :: A(:,:), S(:,:), RHS(:)

   !write(lupri,*) 'redsize, red_space_dim', redsize, red_space_dim
   !write(lupri,*) 'vectorsubspace%offset, vectorsubspace%queuesize', vectorsubspace%offset, vectorsubspace%queuesize
   extradim = .false.
   !matdim = G%nrow
   if (lshift%fixed_shift) then
      mu = -lshift%fixed_shift_param  
   else if (lshift%lshift_by_hlgap) then
      !Do nothing
   else if (lshift%arh_crop) then 
      do i = 1, 5
         call mat_init(x5(i),rowdim,coldim)
      enddo
      if (lshift%arh_truncate) then
         if (lshift%arh_newdamp .and. vectorsubspace%offset == vectorsubspace%queuesize-1 .and. &
            & red_space_dim > vectorsubspace%queuesize-1) then
            !If subspace if full AND if it's not the first time subspace is
            !full....:
            fullsize = vectorsubspace%offset+1
            extradim = .true.
         else
            fullsize = vectorsubspace%offset
         endif
      else
         fullsize = red_space_dim
      endif
      if (lshift%arh_newdamp .and. extradim) then
         if (.not. present(xF) .or. .not. present(sigmaF)) then
            WRITE(lshift%LUPRI,'(/A)') &
            &     'xF, sigmaF must be present with new damping scheme'
            CALL LSQUIT('xF, sigmaF must be present with new damping scheme',lshift%lupri)
         endif
      endif

      call mem_alloc(A,fullsize,fullsize)
      call mem_alloc(S,fullsize,fullsize)
      call mem_alloc(RHS,fullsize)

      A(1:fullsize,1:fullsize) = Ared(1:fullsize,1:fullsize)
      S(1:fullsize,1:fullsize) = Sred(1:fullsize,1:fullsize)
      RHS(1:fullsize) = Gred(1:fullsize)

      call crop_levelshift(lshift,vectorsubspace,extradim,A,S,RHS,lub,fullsize,x5,mu,xF,sigmaF)

      call mem_dealloc(A)
      call mem_dealloc(S)
      call mem_dealloc(RHS)
      do i = 1, 5
         call mat_free(x5(i))
      enddo
   endif

   if (abs(mu) < lshift%min_lshift) then
      mu = -lshift%min_lshift
   endif

   end subroutine main_levelshift

   !> \brief Determine level shift for CROP solver.
   !> \author S. Host
   !> \date 2007
   subroutine crop_levelshift(lshift,vectorsubspace,extradim,A,S,RHS,lub,reddim,x,mu,xF,sigmaF)
   implicit none
        !> Contains setting for level shift.
        type(lshiftItem),intent(in) :: lshift
        !> Contains subspace of previous trial vectors (if arh_truncate=true)
        TYPE(modFIFO), intent(inout)   :: vectorsubspace
        !> True if new level shift scheme (extra dimension on reduced space). NOT FULLY TESTED
        logical, intent(in)         :: extradim
        !> Size of reduced space     
        integer, intent(in)         :: reddim
        !> Reduced Hessian
        real(realk),intent(in)      :: A(reddim,reddim)
        !> Reduced overlap of trial vectors
        real(realk),intent(in)      :: S(reddim,reddim)
        !> Reduced gradient/right hand side
        real(realk),intent(in)      :: RHS(reddim)
        !> Logical unit number for file containing previous trial vectors (if arh_truncate=false)
        integer, intent(in)         :: lub
        !> Work space
        type(matrix), intent(inout) :: x(5)              
        !> Level shift
        real(realk),intent(inout)   :: mu
        !> Matrix used for new level shift scheme. NOT FULLY TESTED
        type(matrix), intent(in)    :: xF
        !> Matrix used for new level shift scheme. NOT FULLY TESTED
        type(matrix), intent(in)    :: sigmaF
        real(realk)                 :: best_mu, min_eival, dev
        real(realk)                 :: intervals(5)
        integer                     :: i, acceptedsteps, n, IERR, min_eival_index(1)
        logical                     :: step_accepted, get_out
        real(realk), pointer    :: eignum(:), eignum_im(:),eigenvec(:,:), &
                                     & eigdenom(:), Atemp(:,:), Stemp(:,:) 

   !Determine initial intervals:
   !if (abs(mu) < 1.0E-2_realk) then
   !   intervals(1) = -16.0E0_realk ; intervals(2) = -8.0E0_realk ; intervals(3) = -4.0E0_realk 
   !   intervals(4) =  -2.0E0_realk ; intervals(5) =  0.0E0_realk 
   !else
   !   intervals(1) = 2.0E0_realk*mu  ; intervals(2) = 1.5E0_realk*mu ; intervals(3) = mu
   !   intervals(4) = 0.75E0_realk*mu ; intervals(5) = 0.5E0_realk*mu  
   !endif

   if (lshift%info_levelshift) then
      write (lshift%lupri,*) 'A in levelshift:'
      call LS_OUTPUT(A, 1, reddim, 1, reddim, reddim, reddim, 1, lshift%lupri)
      write (lshift%lupri,*) 'S in levelshift:'
      call LS_OUTPUT(S, 1, reddim, 1, reddim, reddim, reddim, 1, lshift%lupri)
   endif

   !Determine lowest eigenval in reduced space:
   call mem_alloc(eignum,reddim)
   call mem_alloc(eignum_im,reddim)
   call mem_alloc(eigdenom,reddim)
   call mem_alloc(eigenvec,reddim,reddim)
   call mem_alloc(Atemp,reddim,reddim)
   call mem_alloc(Stemp,reddim,reddim)
   Atemp = A !A is altered in RG
   Stemp = S
   call RGG(reddim,reddim,Atemp,Stemp,eignum,eignum_im,eigdenom,0,eigenvec,IERR) !Hess is not always symmetric
   if (IERR /= 0) then
      WRITE(lshift%LUPRI,'(/A, i4)') &
      &     'Problem in RGG, IERR =', IERR
      CALL LSQUIT(' Problem in RGG',lshift%lupri)
   endif
   do i = 1, reddim
      if (abs(eignum_im(i)/eigdenom(i)) > 1.0E-6_realk) &
         & write (lshift%lupri,*) 'WARNING: Imaginary eigenvalue in reduced Hessian, size:', abs(eignum_im(i)/eigdenom(i))
   enddo
   !Create eigenvalues from numerators eignum and denominators eigdenom
   do i = 1, reddim
      !if (lshift%info_levelshift) WRITE(lshift%LUPRI,*) 'eignum(i): ', eignum(i)
      !if (lshift%info_levelshift) WRITE(lshift%LUPRI,*) 'eigdenom(i): ', eigdenom(i)
      if (abs(eigdenom(i))>10E-12_realk .and. abs(eignum(i))>10E-12_realk) then
         eignum(i) = eignum(i)/eigdenom(i) 
      else
         !Extremely small divided by extreme small gives random results!
         !Set it to something large...
         eignum(i) = 1E3_realk
      endif
      !if (lshift%info_levelshift) WRITE(lshift%LUPRI,*) 'eigenvalue: ', eignum(i)
      !if (lshift%info_levelshift) WRITE(lshift%LUPRI,*)
   enddo
   call find_min_eival(reddim,eignum,min_eival,min_eival_index)
   if (lshift%info_levelshift) write(lshift%lupri,*) 'Smallest eigenvalue/upper bound:', min_eival

   if (min_eival >= 0.0E0_realk .or. abs(min_eival) < 1.0E-2_realk)  then
      intervals(1) = -16.0E0_realk ; intervals(2) = -8.0E0_realk ; intervals(3) = -4.0E0_realk 
      intervals(4) =  -2.0E0_realk ; intervals(5) =  0.0E0_realk 
   else if (reddim > 1) then
      intervals(1) = min_eival*10.0E0_realk  ; intervals(2) = min_eival*7.5E0_realk ; intervals(3) = min_eival*5.0E0_realk
      intervals(4) = min_eival*2.5E0_realk ; intervals(5) = min_eival
   else !quickfix: If reddim = 1, we cannot set mu=eigenval, because then H = A - mu*S = 0 
      intervals(1) = min_eival*10.0E0_realk  ; intervals(2) = min_eival*7.5E0_realk ; intervals(3) = min_eival*5.0E0_realk
      intervals(4) = min_eival*2.5E0_realk ; intervals(5) = intervals(4)
   endif
   call mem_dealloc(eignum)
   call mem_dealloc(eignum_im)
   call mem_dealloc(eigdenom)
   call mem_dealloc(eigenvec)
   call mem_dealloc(Atemp)
   call mem_dealloc(Stemp)

   acceptedsteps = 0
   n = 0
   mu = 0.0E0_realk
   do 
      n = n + 1
      if (n == 100) call LSquit('Failed to determine level shift!!!',lshift%lupri)
      call crop_levelshift_getx(lshift,vectorsubspace,extradim,intervals,A,S,RHS,lub,reddim,x,xF,sigmaF)
      step_accepted = .false.
      dev = abs((intervals(1)-intervals(2))/intervals(1))
      call crop_levelshift_check(lshift,x,step_accepted,intervals,best_mu,get_out)
      if (step_accepted) then
         !if (n > 1) then
         !   diff = abs((mu-best_mu)/mu)
         !   if (info_levelshift) write(lupri, '("Optimal mu changed by", F8.2, "%")') diff*100
         !endif
         acceptedsteps = acceptedsteps + 1
         mu = best_mu
      endif
      !dev: difference between intervals must be max 5% of largest mu
      if ((acceptedsteps >= 3 .and. dev < 0.05E0_realk) .or. get_out) exit  
   enddo
   !if (debug_arh_iter_hlgap < 0.0E0_realk .and. mu > debug_arh_iter_hlgap) then
   !   !Don't accept smaller damping than HOMO-LUMO gap
   !   mu = debug_arh_iter_hlgap
   !   if (info_levelshift) write(lshift%lupri,*) 'Damping determined from HOMO-LUMO gap: mu =', mu
   !endif
   end subroutine crop_levelshift

   !> \brief Get five trial solutions for five trial level shifts.
   !> \author S. Host
   !> \date 2007
   subroutine crop_levelshift_getx(lshift,vectorsubspace,extradim,intervals,A,S,G,lub,reddim,x,xF,sigmaF)
   implicit none
        !> Contains setting for level shift.
        type(lshiftItem),intent(in) :: lshift
        !> Contains subspace of previous trial vectors (if arh_truncate=true)
        TYPE(modFIFO), intent(inout)   :: vectorsubspace
        !> True if new level shift scheme (extra dimension on reduced space). NOT FULLY TESTED
        logical, intent(in)         :: extradim
        !> Contains the five trial level shifts
        real(realk),intent(in)      :: intervals(5)
        !> Size of reduced space
        integer, intent(in)         :: reddim 
        !> Reduced Hessian
        real(realk),intent(in)        :: A(reddim,reddim)
        !> Reduced overlap of trial vectors
        real(realk),intent(in)        :: S(reddim,reddim)
        !> Reduced gradient
        real(realk),intent(in)        :: G(reddim)
        !> Logical unit number for file containing previous trial vectors (if arh_truncate=false)
        integer, intent(in)         :: lub
        !> The five trial solutions (output)
        type(matrix), intent(inout) :: x(5)
        !> Matrix used for new level shift scheme. NOT FULLY TESTED
        type(matrix), intent(in)    :: xF
        !> Matrix used for new level shift scheme. NOT FULLY TESTED
        type(matrix), intent(in)    :: sigmaF
        real(realk)                 :: mu
        integer                     :: i, j, IERR, nsize, nvecs
        real(realk), pointer    :: H(:,:), xtrial(:,:), RHS(:)
        integer, pointer        :: IPIV(:)
        type(matrix)                :: bn
        type(matrix), pointer       :: xpointer, dummypointer
        logical                     :: firstwrong,OnMaster
        real(realk),external :: DNRM2
        OnMaster = .TRUE.
        IERR=0
   !if (extradim) then
   !   nsize = reddim-1
   !else
      nsize = reddim
   !endif
   if (extradim) then
      nvecs = reddim-1
   else
      nvecs = reddim
   endif

   firstwrong = .false.

   call mem_alloc(H,nsize,nsize)
   call mem_alloc(IPIV,nsize)
   call mem_alloc(xtrial,nsize,5)
   call mem_alloc(RHS,nsize)
   if (.not. lshift%arh_truncate) then
      call mat_init(bn,x(1)%nrow,x(1)%nrow)
   endif

   if (lshift%info_levelshift) write(lshift%lupri, '("Trial mus: ", 5F10.5)') &
   & intervals(1),intervals(2),intervals(3),intervals(4),intervals(5)

   !write (lshift%lupri,*) 'A in levelshift:'
   !call LS_OUTPUT(A, 1, nsize, 1, nsize, reddim, reddim, 1, lshift%lupri)
   !write (lshift%lupri,*) 'S in levelshift:'
   !call LS_OUTPUT(S, 1, nsize, 1, nsize, reddim, reddim, 1, lshift%lupri)

   !Solve in reduced space for different values of mu:
   do i = 1, 5
      RHS = G !RHS is altered in DGESV
      mu = intervals(i)
      H(1:nsize,1:nsize) = A(1:nsize,1:nsize) - mu*S(1:nsize,1:nsize)
!      write (lshift%lupri,*) 'H in levelshift, i =', i, ' and mu =', mu
!      call LS_OUTPUT(H, 1, nsize, 1, nsize, nsize, nsize, 1, lshift%lupri)
!      write (lshift%lupri,*) 'RHS in levelshift:'
!      call LS_OUTPUT(RHS, 1, nsize, 1, 1, nsize, 1, 1, lshift%lupri)
      IERR = 0
      call DGESV(nsize, 1, H, nsize, IPIV, RHS, nsize, IERR) !Solution vector is found in RHS.
!      write (lshift%lupri,*) 'RHS in levelshift:',IERR
!      call LS_OUTPUT(RHS, 1, nsize, 1, 1, nsize, 1, 1, lshift%lupri)
      !additional test of sigularity
      if (IERR .EQ. 0) then
         !norm of RHS should be smaller than 1.0E+10 
         if(DNRM2(nsize,RHS,1).GT.1.0E+10)THEN
            IERR = nsize+1
         ENDIF
      endif
      if (IERR /= 0) then
         if (i > 1) then
            !if IERR /= 0, it probably just means that we by accident hit a mu that
            !makes H singular - in that case, skip solution
            WRITE(lshift%LUPRI,'(/A, i4)') &
            &     'Problem in DGESV, IERR = ', IERR
            !CALL LSQUIT(' Problem in DGESV',lshift%lupri)
            write(lshift%lupri,*) 'Use previous solution'
            xtrial(1:nsize,i) = xtrial(1:nsize,i-1)
         else
            firstwrong = .true.
         endif
      else
         xtrial(1:nsize,i) = RHS
      endif
      call mat_zero(x(i))
   enddo

   if (firstwrong) then !It the first solution failed, set it equal to the second
      xtrial(1:nsize,1) = xtrial(1:nsize,2)
   endif

   !Obtain solutions in real space:
   if (.not. lshift%arh_truncate) rewind(lub)
   do i = 1, nvecs
      if (lshift%arh_truncate) then
         call get_from_modFIFO(vectorsubspace, i, xpointer, dummypointer)
         do j = 1, 5
            call mat_daxpy(xtrial(i,j), xpointer, x(j))
         enddo
      else
         call mat_read_from_disk(lub,bn,OnMaster)
         do j = 1, 5
            call mat_daxpy(xtrial(i,j), bn, x(j))
         enddo
      endif
   enddo

   if (extradim) then
      do j = 1, 5
         call mat_daxpy(xtrial(nvecs+1,j), xF, x(j))
      enddo      
   endif

   call mem_dealloc(H)
   call mem_dealloc(IPIV)
   call mem_dealloc(xtrial)
   call mem_dealloc(RHS)
   if (.not. lshift%arh_truncate) then
      call mat_free(bn)
   endif
   end subroutine crop_levelshift_getx

   !> \brief Check size of five trial solutions corresponding to five trial level shifts.
   !> \author S. Host
   !> \date 2007
   subroutine crop_levelshift_check(lshift,x,step_accepted,intervals,best_mu,get_out)
   implicit none
        !> Contains setting for level shift.
        type(lshiftItem),intent(in) :: lshift
        !> The five trial solutions
        type(matrix), intent(in)    :: x(5)
        !>  Was acceptable step found?
        logical, intent(out)        :: step_accepted
        !> New trial level shifts are determined
        real(realk),intent(inout)   :: intervals(5)
        !> The mu corresponding to the acceptable step (if found)
        real(realk),intent(inout)   :: best_mu
        !True if it doesn't make sense to continue with level shift algorithm
        logical, intent(inout)      :: get_out
        real(realk)                 :: intervalsize, xsize, xsize_old, TR !trust radius
        real(realk)                 :: upper_bound, lower_bound
        integer                     :: j
        xsize_old = 0.0E0_realk

   if (lshift%optxelm) then
      TR = lshift%max_element
   else
      TR = lshift%max_step
   endif

   get_out = .false.
   !we have 5 solutions. determine xnorm/xmax of these:
   do j = 1,5
      if (lshift%optxelm) then
         call mat_max_elm(x(j), xsize)
! intervals used before defined
!         if (lshift%info_levelshift) then
!            write(lshift%lupri, '("Point info: max elm = ", F10.5, " and mu = ", &
!                  & F10.5)') xsize, intervals(j)
!         endif
      else 
         xsize = sqrt(mat_sqnorm2(x(j)))
! intervals used before defined
!         if (lshift%info_levelshift) then
!            write(lshift%lupri, '("Point info: Xnorm = ", F10.5, " and mu = ", &
!                  & F10.5)') xsize, intervals(j)
!         endif
      endif
      if (j > 1 .and. xsize_old > xsize) then
         !previous mu was good, but now we passed the Hessian eigenvalue, and mu is too small
         step_accepted = .true.
         best_mu = intervals(j-1)
         get_out = .true.
      else if (j == 1 .and. xsize > TR) then
         !All mu's were too small, try 5 new intervals
         if (lshift%info_levelshift) write(lshift%lupri,*) "All levelshift were too small, try 5 new intervals"
         !STOP 'all levelshifts too small, this looks unhealthy!' 
         intervals(5) = intervals(1)
         intervals(4) = 2*intervals(5)
         intervals(3) = 3*intervals(5)
         intervals(2) = 4*intervals(5)
         intervals(1) = 5*intervals(5)
         exit
      else if (xsize > TR) then
         !Correct interval found, divide in smaller intervals
         step_accepted = .true.
         best_mu = intervals(j-1)
         upper_bound = intervals(j)
         lower_bound = intervals(j-1)
         intervalsize = (upper_bound-lower_bound)/4
         intervals(1) = lower_bound
         intervals(2) = intervals(1) + intervalsize
         intervals(3) = intervals(2) + intervalsize
         intervals(4) = intervals(3) + intervalsize
         intervals(5) = upper_bound
         exit
      else if (j == 5 .and. xsize < TR) then
         !Tested mu's all give smaller step than trust radius
         !so we accept mu to be the lowest Hessian eigenvalue in the reduced space
         step_accepted = .true.
         best_mu = intervals(j)
         get_out = .true.
         !if (abs(intervals(j)) < 1.0E-2_realk) then
         !   best_mu = 0.0E0_realk
         !   get_out = .true.
         !else
         !   lower_bound = intervals(5)
         !   intervalsize = -lower_bound/4
         !   intervals(1) = lower_bound
         !   intervals(2) = intervals(1) + intervalsize
         !   intervals(3) = intervals(2) + intervalsize
         !   intervals(4) = intervals(3) + intervalsize
         !   intervals(5) = 0.0E0_realk
         !endif
         exit
      endif
      xsize_old = xsize
   enddo

   end subroutine crop_levelshift_check

end module levelshift_mod

