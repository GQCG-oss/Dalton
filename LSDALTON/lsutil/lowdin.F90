module lowdin_module
! Implements Lowdin decomposition of S to S^1/2 and S^-1/2
! by diagonalization or in linearly scaling Schulz iterations.
! B. Jansik, Arhus, Jun 2006
use precision
use matrix_module
use matrix_operations
use matrix_operations_aux
use memory_handling
!use lstiming!, only: lstimer
contains

   subroutine lowdin_schulz(S, S_sqrt, S_minus_sqrt,lupri)
     use lstiming
     implicit none
     type(Matrix)                 :: S, S_sqrt, S_minus_sqrt, T1, T2
     real (realk)                 :: l2, alpha, beta, t1_converged, testx, emax, emin
     real (realk)                 :: tstart, tend, dummy
     integer                      :: iter, lupri
     logical                      :: converged = .false., diverging = .false.
#define MAX_ITER 50

   ! temporary arrays
     call mat_init(T1,S%nrow,S%ncol)
     call mat_init(T2,S%nrow,S%ncol)
 

   ! initial guess
     call mat_assign(S_sqrt,S)
     call mat_assign(T1,S)
     call mat_identity(S_minus_sqrt)

   ! get extremal eigenvalues

     if (.not.ASSOCIATED(S%raux))THEN
        call mat_eigenvalues_to_aux(.FALSE.,S)
     endif

     emax = S%raux(1)
     emin = S%raux(2)



   ! T1 should converge to this parameter
     t1_converged=sqrt(1.0E0_realk*S%ncol)

   ! do iterations
     do iter=1, MAX_ITER

       call ls_gettim(dummy,tstart)

       !update emax, emin
       if(iter.ne. 1) then
          emax=1; emin=(1.0E0_realk/4.0E0_realk)*emin*(3-emin)*(3-emin)
       end if

       l2   = 2.0E0_realk/(emax+emin)
       emax = l2*emax
       emin = l2*emin

       !intermediate scaling parameters
       alpha = -0.5E0_realk*l2*sqrt(l2)
       beta  =  1.5E0_realk*sqrt(l2)

       call mat_mul(S_sqrt,T1,'n','n',alpha,0.0E0_realk,T2)
       call mat_scal(beta,S_sqrt)
       call mat_daxpy(1.0E0_realk,T2, S_sqrt)

       call mat_mul(T1,S_minus_sqrt,'n','n',alpha,0.0E0_realk,T2)
       call mat_scal(beta,S_minus_sqrt)
       call mat_daxpy(1.00E0_realk,T2, S_minus_sqrt)

       call mat_mul(S_minus_sqrt,S_sqrt,'n','n',1.0E0_realk,0.0E0_realk, T1)

       !check convergence
       testx = abs(sqrt(mat_sqnorm2(T1))-t1_converged)
       converged = testx .le.  1.0E-7_realk 
       call ls_gettim(dummy,tend)

       write(*,'(I2,A,E12.4,A,F12.4,A,F12.4,A,E12.4,A,F8.1,A)')&
      &iter,    " Norm= ", testx, " l2= ", l2," emax= ", emax," emin= ",&
      &emin, " time= ", tend-tstart, " sec"
       if (converged) exit


     enddo

 
     ! release temporary memory
      call mat_free(T1)
      call mat_free(T2)


      if (converged) then
         write(lupri,'(1X,A,1X,I2,1X,A/,1X,A,E12.4)') &
        &"Iterative Lowdin decomposition converged in", iter, "itarations",&
        &"Residual norm         :", testx
!        write(*,'(1X,A,1X,I2,1X,A/,1X,A,E12.4)') &
!       &"Iterative Lowdin decomposition converged in", iter, "itarations",&
!       &"Residual norm         :", test
      else
         !iterative scheme not converged in MAX_ITER iterations
         write (lupri,'(/A,I4,1X,A//,A/,A/,A/,A/,A,1X,E12.4)') &
        &"Iterative LOWDIN DECOMPOSITION did NOT "//&
        &"CONVERGE in",MAX_ITER,"iterations","Please use one of",&
        &".LWDIAG", ".CHOLESK", "options in DALTON input",&
        &"Residual norm    :", testx
         call lsquit("Iterative LOWDIN DECOMPOSITION did NOT CONVERGED",lupri)
      end if



   end subroutine lowdin_schulz

!BRANO PLEASE ADD LAPACK ROUTINES TO DALTONS OWN MATH LIB SO THAT PEOPLE CAN COMPILE 
!WITHOUT LINKING TO MKL!!!!!
! Quadruple precision not supported by gfortran-4.1.1
!!$#ifndef VAR_OPEN64
!!$#ifndef GFORTRAN
!!$#ifndef VAR_PGF90
!!$   subroutine lowdin_qschulz(cfg_unres, lupri, S, S_sqrt, S_minus_sqrt)
!!$     implicit none
!!$     logical,intent(in)           :: cfg_unres
!!$     integer,intent(in)           :: lupri
!!$     integer, parameter           :: quadk = selected_real_kind(16,99)
!!$
!!$     type(Matrix)                 :: S, S_sqrt, S_minus_sqrt
!!$     real(realk),  allocatable    :: dY(:), dZ(:)
!!$     real(quadk),  allocatable    :: Y(:),Z(:,:),T(:),T1(:)
!!$
!!$     integer      :: i, info, iter, n, n2
!!$     logical      :: converged
!!$     real(quadk)  :: qdot,l2,testx,ALPHA,BETA
!!$     real(realk)  :: emin,emax,tstart,tend,dummy
!!$     external qdot
!!$
!!$      if (cfg_unres) STOP 'lowdin_qschulz not implemented for open shell'
!!$      n = S%ncol; n2 = n*n
!!$
!!$      ! get extremal  eigenvalues
!!$      emax = S%raux(1); emin = S%raux(2)
!!$      !call dsyev('N','U',n,dY,n,dZ,Z,n2,info)
!!$      !emin = dZ(1); emax = dZ(n)
!!$
!!$      ! allocate temporary arrays
!!$      allocate(Y(n2),T(n2),T1(n2),dY(n2))
!!$
!!$      !convert data to quadruple precision
!!$      call mat_to_full(S,1.0E0_realk,dY)
!!$      call d2qconv(n2,Y,dY)
!!$      call d2qconv(n2,T,dY)
!!$      deallocate(dY)
!!$
!!$     !set Z to unit
!!$
!!$      allocate(Z(n,n))
!!$      call qzero(Z,n2)
!!$      DO i=1,n
!!$         Z(i,i) = 1.0Q0
!!$      END DO
!!$
!!$      !iterations
!!$      DO iter=1, MAX_ITER
!!$
!!$        call ls_gettim(dummy,tstart)
!!$
!!$        if(iter.ne. 1) then
!!$           emax=1; emin=(1.0E0_realk/4.0E0_realk)*emin*(3-emin)*(3-emin)
!!$        end if
!!$
!!$        l2   = 2.0Q0/(emax+emin)
!!$        emax = l2*emax
!!$        emin = l2*emin
!!$
!!$
!!$        ALPHA= -0.5Q0*l2*sqrt(l2)
!!$        BETA = 1.5Q0*sqrt(l2)
!!$
!!$        call qgemm('n','n',n,n,n,ALPHA,Y,n,T,n,0.0Q0,T1,n)
!!$        call qscal(n2,BETA,Y,1)
!!$        call qaxpy(n2,1.0Q0,T1,1,Y,1)
!!$
!!$        call qgemm('n','n',n,n,n,ALPHA,T,n,Z,n,0.0Q0,T1,n)
!!$        call qscal(n2,BETA,Z,1)
!!$        call qaxpy(n2,1.0Q0,T1,1,Z,1)
!!$
!!$        call qgemm('n','n',n,n,n,1.0Q0,Z,n,Y,n,0.0Q0,T,n)
!!$
!!$        ! check convergence
!!$        testx = abs(sqrt(qdot(n2,T,1,T,1)) - sqrt(1.0Q0*n))
!!$        converged = testx.le. 1Q-10
!!$        call ls_gettim(dummy,tend)
!!$
!!$        write(*,'(I2,A,E10.4,A,F6.4,A,F6.4,A,E10.4,A,F8.1,A)')&
!!$       &iter," Norm= ", testx, " l2= ", l2," emax= ", emax," emin= ",&
!!$       &emin, " time= ", tend-tstart, " sec"
!!$        if (converged) exit
!!$
!!$
!!$      END DO
!!$
!!$      deallocate(T,T1)      
!!$
!!$!     convert back
!!$      allocate(dY(n2))
!!$      call q2dconv(n2,dY,Y)
!!$      call mat_set_from_full(dY,1.0E0_realk,S_sqrt)
!!$
!!$      call q2dconv(n2,dY,Z)
!!$      call mat_set_from_full(dY,1.0E0_realk,S_minus_sqrt)
!!$
!!$      deallocate(dY,Y,Z)
!!$
!!$      if (converged) then
!!$         write(lupri,'(1X,A,1X,I2,1X,A/,1X,A,E12.4)') &
!!$        &"Quadruple precision iterative Lowdin decomposition converged in", iter, "itarations",&
!!$        &"Residual norm         :", testx
!!$!        write(*,'(1X,A,1X,I2,1X,A/,1X,A,E12.4)') &
!!$!       &"Iterative Lowdin decomposition converged in", iter, "itarations",&
!!$!       &"Residual norm         :", test
!!$      else
!!$         !iterative scheme not converged in MAX_ITER iterations
!!$         write (lupri,'(/A,I4,1X,A//,A/,A/,A/,A/,A,1X,E12.4)') &
!!$        &"Iterative QUADRUPLE PRECISION LOWDIN DECOMPOSITION did NOT "//&
!!$        &"CONVERGE in",MAX_ITER,"iterations","Please use one of",&
!!$        &".LWDIAG", ".CHOLESK", "options in DALTON input",&
!!$        &"Residual norm    :", testx
!!$         call lsquit("Iterative LOWDIN DECOMPOSITION did NOT CONVERGED",lupri)
!!$      end if
!!$
!!$
!!$
!!$   end subroutine lowdin_qschulz
!!$#endif
!!$#endif
!!$#endif

   ! Given matrix S, computes sqrt(s) and S^{-1/2}.
   subroutine lowdin_diag(n, S,S_sqrt, S_minus_sqrt, lupri)
     implicit none
     integer, intent(in)          ::n
     real(realk), intent(in)      :: S(n,n)
     real(realk), intent(out)     :: S_sqrt(n,n), S_minus_sqrt(n,n)
     real(realk),allocatable      :: T1(:,:), T2(:,:), T3(:,:)
     real(realk),allocatable      :: eigen_sqrt(:), eigen_minus_sqrt(:)
     integer                      :: i, infdiag, lupri
     character*70                 :: msg
     integer                      :: lwork
     real(realk), dimension(:), allocatable :: work
#ifdef VAR_LSESSL
     integer :: liwork
     integer, pointer :: iwork(:)
     liwork = -1
#endif
     infdiag=0

     allocate(T1(n,n),T2(n,n))
     allocate(eigen_sqrt(n),eigen_minus_sqrt(n))

     t1 = S

!     lwork = max(n*n,5*n)
!============================================================
! we inquire the size of lwork (NOT max(n*n,5*n)
     lwork = -1
     allocate(work(5))
#ifdef VAR_LSESSL
     call mem_alloc(iwork,5)
     call dsyevd('V','U',n,T1,n,eigen_sqrt,work,lwork,iwork,liwork,infdiag)
     liwork = iwork(1)
     call mem_dealloc(iwork)
     call mem_alloc(iwork,liwork)
#else
     call dsyev('V','U',n,T1,n,eigen_sqrt,work,lwork,infdiag)
#endif
     lwork = NINT(work(1))
     deallocate(work)
!=============================================================     

     allocate(work(lwork))
     !diagonalization
#ifdef VAR_LSESSL
     call dsyevd('V','U',n,T1,n,eigen_sqrt,work,lwork,iwork,liwork,infdiag)
     call mem_dealloc(iwork)
#else
     call dsyev('V','U',n,T1,n,eigen_sqrt,work,lwork,infdiag)
#endif
     deallocate(work)

     if(infdiag.ne. 0) then
        write(lupri,*) 'lowdin_diag: dsyev failed, info=',infdiag
        call lsquit('lowdin_diag: diagonalization failed.',lupri)
     end if

     !compute squareroot S^1/2 = V*E^1/2*V'
     ! E^1/2, -1/2

     do i=1,n
        if (eigen_sqrt(i).le. 0E0_realk) then
           write(msg,'(A,1X,I4,A,1X,E14.6)') &
                &'Matrix not positive definite! Eigenvalue(',i,') =',&
                &eigen_sqrt(i)
           call lsquit(msg,lupri)
        endif
        eigen_sqrt(i)       = sqrt(eigen_sqrt(i))
        eigen_minus_sqrt(i) = 1.0E0_realk/eigen_sqrt(i)
     enddo


     ! V*E^1/2
     do i=1,n
        t2(:,i) = t1(:,i)*eigen_sqrt(i)
     enddo

     ! S^1/2
     call dgemm('n','t', n,n,n, 1.0E0_realk, t2,n, t1,n, 0.0E0_realk,S_sqrt,n)

     !V*E^-1/2
     do i=1,n
        t2(:,i) = t1(:,i)*eigen_minus_sqrt(i)
     enddo

     ! S^-1/2
     call dgemm('n','t', n,n,n, 1.0E0_realk, t2,n, t1,n, 0.0E0_realk,S_minus_sqrt,n)

     !deallocate
     deallocate(eigen_sqrt,eigen_minus_sqrt,t1,t2)

!     !verify
!     allocate(T3(n,n))
!     call dgemm('n','n',n,n,n,1.0E0_realk,S_sqrt,n,S_minus_sqrt,n,0.0E0_realk,T3,n)
!     print*,'S_sqrt,n*S_minus_sqrt'
!     call ls_output(T3,1,n,1,n,n,n,1,6)
!     deallocate(T3)

   end subroutine lowdin_diag

   ! Given matrix S, computes S^{-1/2}. OVERWRITE S
   subroutine lowdin_diag_S_minus_sqrt(n, S, S_minus_sqrt, lupri)
     implicit none
     integer, intent(in)          :: n
     real(realk), intent(inout)   :: S(n,n)
     real(realk), intent(inout)   :: S_minus_sqrt(n,n)
     real(realk),pointer          :: T(:,:)
     real(realk)                  :: eigen_minus_sqrt(n)
     integer                      :: i, infdiag, lupri
     character*70                 :: msg
     integer                      :: lwork,j
!     real(realk) :: TS,TE
     real(realk), pointer :: work(:)
#ifdef VAR_LSESSL
     integer :: liwork
     integer, pointer :: iwork(:)
     liwork = -1
#endif
     infdiag=0
!     call LSTIMER('START ',TS,TE,lupri,.TRUE.)

     ! we inquire the size of lwork (NOT max(n*n,5*n)
     lwork = -1
     call mem_alloc(work,5)
#ifdef VAR_LSESSL
     call mem_alloc(iwork,5)
     call dsyevd('V','U',n,S,n,eigen_minus_sqrt,work,lwork,iwork,liwork,infdiag)
     liwork = iwork(1)
     call mem_dealloc(iwork)
     call mem_alloc(iwork,liwork)
#else
     call dsyev('V','U',n,S,n,eigen_minus_sqrt,work,lwork,infdiag)
#endif
     lwork = NINT(work(1))
     call mem_dealloc(work)
!=============================================================
     call mem_alloc(work,lwork)
     !diagonalization
#ifdef VAR_LSESSL
     call dsyevd('V','U',n,S,n,eigen_minus_sqrt,work,lwork,iwork,liwork,infdiag)
     call mem_dealloc(iwork)
#else
     call dsyev('V','U',n,S,n,eigen_minus_sqrt,work,lwork,infdiag)
#endif
     call mem_dealloc(work)
!     call LSTIMER('LS: DSYEV ',TS,TE,lupri,.TRUE.)

     if(infdiag.ne. 0) then
        write(lupri,*) 'lowdin_diag: dsyev failed, info=',infdiag
        call lsquit('lowdin_diag: diagonalization failed.',lupri)
     end if

     !compute squareroot S^(-1/2) = V*E^(-1/2)*V'
     call mem_alloc(T,n,n)

     !$OMP PARALLEL DO PRIVATE(i,j) SHARED(eigen_minus_sqrt,S,T)
     do i=1,n
!        if (eigen_minus_sqrt(i).le. 0E0_realk) then
!           write(msg,'(A,1X,I4,A,1X,E14.6)') &
!                &'Matrix not positive definite! Eigenvalue(',i,') =',&
!                &eigen_minus_sqrt(i)
!           call lsquit(msg,lupri)
!        endif
        eigen_minus_sqrt(i) = 1.0E0_realk/sqrt(eigen_minus_sqrt(i))
        !T = V*E^-1/2
        do j=1,n
           T(j,i) = S(j,i)*eigen_minus_sqrt(i)
        enddo
     enddo
     !$OMP END PARALLEL DO 
!     call LSTIMER('LS: E^(-1/2) ',TS,TE,lupri,.TRUE.)

     ! S^-1/2 = T*V'
     call dgemm('n','t',n,n,n,1.0E0_realk,T,n,S,n,0.0E0_realk,S_minus_sqrt,n)
     !deallocate
     call mem_dealloc(T)
!     call LSTIMER('LS: DGEMM ',TS,TE,lupri,.TRUE.)

   end subroutine lowdin_diag_S_minus_sqrt

   ! Given matrix S, computes S^{-1}. OVERWRITE S
   subroutine lowdin_diag_S_minus1(n, S, S_minus1, lupri)
     implicit none
     integer, intent(in)          :: n
     real(realk), intent(inout)   :: S(n,n)
     real(realk), intent(inout)   :: S_minus1(n,n)
     real(realk),pointer          :: T(:,:)
     real(realk)                  :: eigen_minus1(n)
     integer                      :: i, infdiag, lupri
     character*70                 :: msg
     integer                      :: lwork,j
!     real(realk) :: TS,TE
     real(realk), pointer :: work(:)
#ifdef VAR_LSESSL
     integer :: liwork
     integer, pointer :: iwork(:)
     liwork = -1
#endif
     infdiag=0
!     call LSTIMER('START ',TS,TE,lupri,.TRUE.)

     ! we inquire the size of lwork (NOT max(n*n,5*n)
     lwork = -1
     call mem_alloc(work,5)
#ifdef VAR_LSESSL
     call mem_alloc(iwork,5)
     call dsyevd('V','U',n,S,n,eigen_minus1,work,lwork,iwork,liwork,infdiag)
     liwork = iwork(1)
     call mem_dealloc(iwork)
     call mem_alloc(iwork,liwork)
#else
     call dsyev('V','U',n,S,n,eigen_minus1,work,lwork,infdiag)
#endif
     lwork = NINT(work(1))
     call mem_dealloc(work)
!=============================================================
     call mem_alloc(work,lwork)
     !diagonalization
#ifdef VAR_LSESSL
     call dsyevd('V','U',n,S,n,eigen_minus1,work,lwork,iwork,liwork,infdiag)
     call mem_dealloc(iwork)
#else
     call dsyev('V','U',n,S,n,eigen_minus1,work,lwork,infdiag)
#endif
     call mem_dealloc(work)

     if(infdiag.ne. 0) then
        write(lupri,*) 'lowdin_diag: dsyev failed, info=',infdiag
        call lsquit('lowdin_diag: diagonalization failed.',lupri)
     end if

     !compute squareroot S^(-1) = V*E^(-1)*V'
     call mem_alloc(T,n,n)

     !$OMP PARALLEL DO PRIVATE(i,j) SHARED(eigen_minus1,S,T)
     do i=1,n
        eigen_minus1(i) = 1.0E0_realk/eigen_minus1(i)
        !T = V*E^-1/2
        do j=1,n
           T(j,i) = S(j,i)*eigen_minus1(i)
        enddo
     enddo
     !$OMP END PARALLEL DO 
     ! S^-1 = T*V'
     call dgemm('n','t',n,n,n,1.0E0_realk,T,n,S,n,0.0E0_realk,S_minus1,n)
     !deallocate
     call mem_dealloc(T)

   end subroutine lowdin_diag_S_minus1

end module lowdin_module
