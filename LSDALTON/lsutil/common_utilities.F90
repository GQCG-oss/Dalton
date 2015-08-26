!> @file 
!> \brief Contains common utility routines for real(realk) matrices.
!>
!> These routines depend only on the module precision. 
!> UNDER NO CIRCUMSTANCES ARE YOU ALLOWED TO INTRODUCE OTHER DEPENDENCIES!
!> If your routine depends on other modules than precision, it means
!> that it belongs somewhere else than here.
!> 

   !> \brief Interface to RGG for diagonalization of real general matrix, A*x = mu*S*X. Return eigenvalues and eigenvectors.
   !> \author S. Host
   !> \date 2005
   !> \param A Matrix to be diagonalized
   !> \param S Overlap matrix to be diagonalized
   !> \param ndim A, S are ndim x ndim matrices
   !> \param lupri Logical unit number for printing output
   !> \param eival The eigenvalues of A
   !> \param eigenvec The eigenvectors of A
   subroutine rgg_interface(A,S,ndim,lupri,eival,eigenvec)
   use precision
   implicit none
        integer, intent(in)      :: ndim, lupri
        real(realk), intent(in)  :: A(ndim,ndim), S(ndim,ndim)
        real(realk), intent(out) :: eival(ndim), eigenvec(ndim,ndim)
        real(realk), allocatable :: eival_r(:), eival_i(:), eival_denom(:)
        integer                  :: IERR, i

   allocate(eival_r(ndim),eival_i(ndim),eival_denom(ndim))

   call RGG(ndim,ndim,A,S,eival_r,eival_i,eival_denom,1,eigenvec,IERR)

   if (IERR /= 0) then
     WRITE(LUPRI,*) &
      &     'Problem in RGG (rgg_interface), IERR =', IERR
      CALL LSQUIT('Problem in RGG (rgg_interface)',lupri)
   endif

   do i = 1, ndim
      if (abs(eival_i(i)/eival_denom(i)) > 1.0E-5_realk) then
         write(lupri,*) 'WARNING: Imaginary eigenvalue in rgg_interface, size =', eival_i(i)/eival_denom(i)
      endif
      eival(i) = eival_r(i)/eival_denom(i)
      !write(lupri,*) 'eigenvalue no.', i, ' is:', eival(i)
   enddo

   deallocate(eival_r,eival_i,eival_denom)
   end subroutine rgg_interface

   !> \brief Interface to DSYEVX for diagonalization of real symmetric A*x = mu*X
   !> \author S. Host
   !> \date 2005
   !> \param A Matrix to be diagonalized
   !> \param ndim A is ndim x ndim matrix
   !> \param print_eivecs If true, eigenvectors are printed in output file lupri
   !> \param lupri Logical unit number for printing output
   !> \param string Character string for identifying matrix in output, if print_eivecs=.true.
   subroutine dsyevx_interface(A,ndim,print_eivecs,lupri,string)
   use precision
   implicit none
        character(6), intent(in), optional :: string
        integer, intent(in)      :: ndim, lupri
        real(realk), intent(in)  :: A(ndim,ndim)
        real(realk), allocatable :: eigenval(:), eigenvec(:,:)
        real(realk), allocatable :: temp(:)
        integer, allocatable     :: itemp(:), IFAIL(:)
        real(realk)              :: VL, VU
        integer                  :: IL, IU, neig, Ltemp, INFO, m
        logical, intent(in)      :: print_eivecs 
   INFO=0
   Ltemp = 8*ndim 
   allocate(eigenval(ndim))
   allocate(eigenvec(ndim,ndim),temp(Ltemp),Itemp(5*ndim),IFAIL(ndim))

   call DSYEVX('V', 'A', 'U', ndim, A, ndim, VL, VU, IL, IU, &
     &  0.0E0_realk, neig, eigenval, eigenvec, ndim, temp, Ltemp, Itemp, &
     &  IFAIL, INFO )

   if (info /= 0) STOP 'Problem in DSYEVX (diag_old)'

   if (present(string)) then
      write(lupri,'("Dim of ", a6, " is: ", i5, " and number of eigenvalues found is ", i5)') string, ndim, neig
   endif
   do m = 1, neig
      write(lupri,'(i6,D15.7)') m, eigenval(m)
   enddo

   if (print_eivecs) then
      write (lupri,*) 'Eigenvectors of A:'
      call LS_OUTPUT(eigenvec, 1, ndim, 1, ndim, ndim, ndim, 1, lupri)
   endif

   deallocate(eigenval)
   deallocate(eigenvec)
   deallocate(temp)
   deallocate(Itemp)
   deallocate(IFAIL)
   end subroutine dsyevx_interface

!!$   SUBROUTINE my_dcopy8(N,DX,INCX,DY,INCY)
!!$     use precision
!!$     implicit none
!!$     integer(kind=8),intent(in)::N
!!$     integer,intent(in)::INCX,INCY
!!$     real(realk),intent(in) :: DX(*)
!!$     real(realk),intent(out) :: DY(*)
!!$     integer(kind=8)::nel,i
!!$     
!!$     do i=1,N,MAXINT
!!$       nel=MAXINT
!!$       if(((N-i)<MAXINT).and.&
!!$         &(mod(N-i+1,MAXINT)/=0))nel=mod(N,MAXINT)
!!$       call dcopy(nel,DX(i),INCX,DY(i),INCY)
!!$     enddo
!!$   end SUBROUTINE my_dcopy8
   !> \brief Interface to DSYGV. Computes eigenvector and eigenvalues of A*x = (lambda)*B*x
   !> \author L. Thogersen
   !> \date 2003
   SUBROUTINE my_DSYGV(N,A,B,eigval,DESC)
   use precision
     implicit none
   ! eigenvectors are in A.
     !> Dimension of A
     integer,          intent(in)    :: N
     !> Input: The symmetric matrix A. Output: The eigenvectors
     real(realk)     , intent(inout) :: A(N,N)
     !> Input: The symmetric positive definite matrix B
     real(realk)     , INTENT(INOUT) :: B(N,N)
     !> The eigenvalues in ascending order
     real(realk)     , INTENT(OUT)   :: eigval(N)
     !> Should contain name of routine from which my_DSYGV is called
     character(20),     INTENT(IN)    :: desc
     !
#ifdef VAR_LSESSL
     integer :: ifail(N),iwrk(5*N),nfound
     real(realk) :: no_ref,tol
     real(realk), external :: DLAMCH
     real(realk) :: Z(N,N)
#endif
     real(realk),allocatable :: wrk(:)
     integer    :: ierr,lwrk
     ierr=0
     allocate(wrk(5))
     lwrk = -1
#ifdef VAR_LSESSL
     no_ref = 0.0E0_realk
     tol    = 0.0E0_realk
     call DSYGVX( 1,'V','A','L', N, A, N, B, N, no_ref,no_ref,no_ref,&
      &no_ref,tol,nfound,eigval, Z, N, wrk, lwrk, iwrk, ifail, ierr)
#else
     call DSYGV(1,'V','L',N,A,N,B,N,eigval,wrk,lwrk,ierr)
#endif 
     if(ierr.ne. 0) THEN
        print *, "DSYGV failed, N = ",N," ierr=", ierr," IN ", DESC
        call lsquit("programming error in my_DSYGV input. workarray inquiry",-1)
     endif
!#ifdef VAR_LSDEBUG
!     ! sometimes the optimal batch sizes do not always work, especially when
!     ! compiled with ifort --debug and --check so I introduced this, PE
!     print *,"WARNING(my_sygv): using minimal lwrk instead of optimal in debug build"
!     lwrk = 3*N-1
!#else
     lwrk = NINT(wrk(1))
!#endif
     deallocate(wrk)
     allocate(wrk(lwrk))
#ifdef VAR_LSESSL
     call DSYGVX( 1,'V','A','L', N, A, N, B, N, no_ref,no_ref,no_ref,&
      &no_ref,tol,nfound,eigval, Z, N, wrk, lwrk, iwrk, ifail, ierr)
     call dcopy(N*N,Z,1,A,1)
#else
     call DSYGV(1,'V','L',N,A,N,B,N,eigval,wrk,lwrk,ierr)
#endif 
     if(ierr.ne. 0) THEN
        print *, "DSYGV failed, N = ",N," ierr=", ierr," IN ", DESC
        IF(ierr.LT.0) print *, "illegal value of argument number ",ABS(ierr)
        IF(ierr.LE.N)then
           print *, 'SYEV failed to converge;'
           print *,  ierr, ' off-diagonal elements of an intermediate'
           print *, 'tridiagonal form did not converge to zero;'
        ENDIF
        IF(ierr.GE.N.AND.ierr.LE.2*N)then
           print *, 'the leading minor of order ',ierr-N,' of B is not positive definite.'
           print *, 'The factorization of B could not be completed and'
           print *, 'no eigenvalues or eigenvectors were computed '
        endif
        call lsquit("programming error in my_DSYGV input.",-1)
     endif
     deallocate(wrk)
   END SUBROUTINE my_DSYGV

   !> \brief Interface to DSYEV. Computes all eigenvalues and, optionally, eigenvectors of a real symmetric matrix A.
   !> \author L. Thogersen
   !> \date 2003
   subroutine my_dsyev(JOBZ, UPLO, N, A, W)
       implicit none
       !> If 'N', compute eigenvalues only; if 'V', compute eigenvalues and eigenvectors.
       CHARACTER,intent(in)          :: JOBZ
       !> If 'U', upper triangle of A is stored; if 'L', lower triangle of A is stored.
       CHARACTER,intent(in)          :: UPLO
       !> Dimension of A
       INTEGER,intent(in)            :: N
       !> Input: The symmetric matrix A(N,N). Output: If JOBZ='V', the orthonormal eigenvectors of the matrix A.
       double precision              :: A(N,N)
       !> The eigenvalues in ascending order (dimension N).
       double precision              :: W(N)
       character(len=70)             :: MSG
       INTEGER                       :: INFO, LWORK, NB, ILAENV
       double precision, allocatable     :: WORK(:)
#ifdef VAR_LSESSL
       integer :: liwork
       integer, allocatable :: iwork(:)
       liwork = -1
       if(JOBZ=='V'.or.JOBZ=='v')then
         liwork = 3+5*N
         allocate(iwork(3+5*N))
       elseif(JOBZ=='N'.or.JOBZ=='n')then
         liwork = 1
         allocate(iwork(1))
       endif
#endif
       INFO=0
   
       !find out optimal work memory size
       NB = ilaenv(1,'DSYTRD',UPLO,N,N,N,N)
       LWORK = MAX((NB+2)*N,3*N-1)
    
       !allocate work memory
       allocate(WORK(LWORK))
   
       !run
#ifdef VAR_LSESSL
       call dsyev(JOBZ,UPLO,N,A,N,W,WORK,LWORK,iwork,liwork,INFO)
       deallocate(iwork)
#else
       call dsyev(JOBZ,UPLO,N,A,N,W,WORK,LWORK,INFO)
#endif
       if (info.ne. 0) then
         write(*,'(1X,A,I3)') 'Error in dsyev in my_DSYGV, info=',info
         call LSquit('Error in dsyev in my_DSYGV',-1)
       endif
   
       deallocate(WORK)
   
       if (INFO .ne. 0) then
          write(MSG,'(A,1X,I4,1X,A)') &
          & 'DSYEV routine failed. INFO=',INFO,'See man dsyev for more information.'
             call LSquit(msg,-1)
       endif
   end subroutine my_dsyev


   !> \brief Locates smallest element in eival array
   !> \author S. Host
   !> \date 2005
   subroutine find_min_eival(ndim,eival,min_eival,min_eival_index)
   use precision
   implicit none
        !> Size of array
        integer, intent(in)      :: ndim
        !> Array for which minimum element should be located
        real(realk), intent(in)  :: eival(ndim)
        !> Minimum value of eival
        real(realk), intent(out) :: min_eival
        !> Position of minimum value i eival
        integer, intent(out)     :: min_eival_index(1)

      min_eival = minval(eival)
      min_eival_index = minloc(eival)

   end subroutine find_min_eival

   !> \brief Calculate HOMO-LUMO gap from orbital energies.
   !> \author L. Thogersen
   !> \date 2003
   !> \param unres True if unrestricted SCF 
   !> \param nocc  Number of occupied orbitals (closed shell)
   !> \param nocca Number of occupied alpha orbitals (if open shell)
   !> \param noccb Number of occupied beta orbitals  (if open shell)
   !> \param eival Array containing orbital energies
   !> \param sz Size of array with orbital energies (eival)
   function HOMO_LUMO_gap(unres,nocc,nocca,noccb,eival,sz)
   use precision
     implicit none
     logical,intent(in)      :: unres
     integer,intent(in)      :: nocc,nocca,noccb,sz
     real(realk), intent(in) :: eival(sz)
     real(realk) :: HOMO_LUMO_gap,gap_a,gap_b
     integer :: ndim

     if(unres) then
        ndim=sz/2
        IF (nocca.GE.1) gap_a=eival(nocca+1)-eival(nocca)
        IF (noccb.GE.1) gap_b=eival(noccb+1 + ndim)-eival(noccb + ndim)
        IF (nocca.EQ.0) THEN
          IF (noccb.EQ.0) THEN
            CALL LSQUIT('Error in HOMO_LUMO_gap, both nocca and noccb are zero',-1)
          ELSE
            HOMO_LUMO_gap = gap_b
          ENDIF
        ELSE IF (noccb.EQ.0) THEN
          HOMO_LUMO_gap = gap_a
        ELSE
          HOMO_LUMO_gap=MIN(gap_a,gap_b)
        ENDIF
     else
        ndim=size(eival)
        HOMO_LUMO_gap=eival(nocc+1)-eival(nocc)
     end if
   end function HOMO_LUMO_gap

   !> \brief provides the HOMO energy from orbital energies.
   !> \author T. Kjaergaard
   !> \date 2012
   !> \param unres True if unrestricted SCF 
   !> \param nocc  Number of occupied orbitals (closed shell)
   !> \param nocca Number of occupied alpha orbitals (if open shell)
   !> \param noccb Number of occupied beta orbitals  (if open shell)
   !> \param eival Array containing orbital energies
   !> \param sz Size of array with orbital energies (eival)
   function HOMO_energy(unres,nocc,nocca,noccb,eival,sz)
   use precision
     implicit none
     logical,intent(in)      :: unres
     integer,intent(in)      :: nocc,nocca,noccb,sz
     real(realk), intent(in) :: eival(sz)
     real(realk) :: HOMO_energy,gap_a,gap_b
     integer :: ndim

     if(unres) then
        ndim=sz/2
        IF (nocca.EQ.0) THEN
          IF (noccb.EQ.0) THEN
            CALL LSQUIT('Error in HOMO_energy, both nocca and noccb are zero',-1)
          ENDIF
          HOMO_energy = eival(noccb + ndim)
        ELSE IF (noccb.EQ.0) THEN
          HOMO_energy = eival(nocca)
        ELSE
          HOMO_energy=MAX(eival(nocca),eival(noccb + ndim))
        ENDIF
     else
        ndim=size(eival)
        HOMO_energy=eival(nocc)
     end if
   end function HOMO_energy

   !> \brief ** Makes the direct product matric C = AxB: [AxB]pq,rs = Apr*Bqs
   !> \author L. Thogersen
   !> \date 2003
   !> \param narow Number of rows in A
   !> \param nacol Number of columns in A
   !> \param B The matrix B
   !> \param nbrow Number of rows in B
   !> \param nbcol Number of columns in B
   !> \param C The output matrix
   subroutine dense_DIRPROD(A,narow,nacol,B,nbrow,nbcol,C)
   use precision
     implicit none
     integer, intent(in) :: narow,nacol,nbrow,nbcol
     real(realk), intent(in) :: A(narow*nacol),B(nbrow*nbcol)
     real(realk), intent(out) :: C(narow*nbrow*nacol*nbcol)
     integer :: i, j, k, l, pq, rs, ncrow
   
     do i = 1,narow
       do j = 1,nacol
         do k = 1,nbrow
           do l = 1,nbcol
             pq = (i-1)*nbrow + k
             rs = (j-1)*nbcol + l
             ncrow = narow*nbrow
             C((rs-1)*ncrow+pq) = A(narow*(j-1)+i)*B(nbrow*(l-1)+k)
            enddo
          enddo
        enddo
      enddo
   
   end subroutine dense_DIRPROD

  !> \brief Set up finite difference Hessian and gradient 
  !> \author L. Thogersen
  !> \date 2003
  subroutine debug_get_fd_hes_grad(func,x,delta,grad,hes)
  use precision
    implicit none
    !> The function of which we want Hessian and gradient
    real(realk), external :: func
    !> Center finite difference around this point
    real(realk), intent(in) :: x
    !> Use x-interval of this size
    real(realk), intent(in) :: delta
    !> Finite difference gradient
    real(realk), intent(out) :: grad
    !> Finite difference Hessian
    real(realk), intent(out) :: hes
    real(realk) :: fd,fmd,f0,f2d,fm2d,hes2,grad2
  
    f0 = func(x)
    fd = func(x+delta)
    fmd = func(x-delta)
    f2d = func(x+2E0_realk*delta)
    fm2d = func(x-2E0_realk*delta)
  
    grad2 = (fd-fmd)/(2.0E0_realk*delta)
    hes2 = (fd+fmd-2.0E0_realk*f0)/delta**2
    grad = (8.0E0_realk*fd-8.0E0_realk*fmd-f2d+fm2d)/(12.0E0_realk*delta)
    hes = (-1.0E0_realk*f2d+16.0E0_realk*fd-30.0E0_realk*f0+16.0E0_realk*fmd-1.0E0_realk*fm2d)/(12.0E0_realk*delta**2)
  
  end subroutine debug_get_fd_hes_grad

  !> \brief Naive implementation of frobenius norm 
  !> \author S. Host
  !> \date 2007
  function frob_norm(A,nrow,ncol)
  use precision
  implicit none
       real(realk)              :: frob_norm
       !> Number of rows in A
       integer, intent(in)      :: nrow
       !> Number of columns in A
       integer, intent(in)      :: ncol
       !> Find Frobenius norm of this matrix
       real(realk), intent(in)  :: A(nrow,ncol)
       integer                  :: i, j

   frob_norm = 0.0E0_realk
   do i = 1,nrow
      do j = 1, ncol
         frob_norm = frob_norm + A(i,j)*A(i,j)
      enddo
   enddo
   frob_norm = sqrt(frob_norm)

  end function frob_norm   

  !> \brief Calculate the square norm of a vector.
  !> \author S. Host
  !> \date 2005
  subroutine vecnorm(c, m, N)
  use precision
  implicit none
       !> Dimension of input vector
       integer,intent(in)       :: m
       !> Calculate norm of this vector
       real(realk), intent(in)  :: c(m)
       !> Norm of vector
       real(realk), intent(inout) :: N
       integer :: i
     N=0
     do i = 1, m
        N = N + c(i)**2
     enddo
     N=SQRT(N)
  end subroutine VECNORM

  !> \brief Sets a real array of length *LENGTH* to zero.
  !> \author H. J. Aa. Jensen. F90'fied by S. Host
  !> \date May 5, 1984
  subroutine ls_dzero(dx, length)
  use precision
  implicit none
       !Length of array
       integer, intent(in)      :: length
       !Array to be nullified
       real(realk), intent(inout) :: dx(length)
       integer                  :: i

     if (length < 0) then
        !do nothing
     else
        do i = 1, length
           dx(i) = 0.0E0_realk
        enddo
     endif
   end subroutine ls_dzero

  !> \brief Sets a real array of length *LENGTH* to zero.
  !> \author H. J. Aa. Jensen. F90'fied by S. Host
  !> \date May 5, 1984
  subroutine ls_dzero8(dx, length)
  use precision
  implicit none
       !Length of array
       integer(kind=long), intent(in) :: length
       !Array to be nullified
       real(realk), intent(inout) :: dx(length)
       integer(kind=long)             :: i

     if (length < 0) then
        !do nothing
     else
        do i = 1, length
           dx(i) = 0.0E0_realk
        enddo
     endif
   end subroutine ls_dzero8

  !> \brief Sets a short integer array of length *LENGTH* to shortzero.
  !> \author T. Kjaergaard
  !> \date May 5, 2011
  subroutine ls_sizero(dx, length)
  use precision
  implicit none
       !Length of array
       integer, intent(in)      :: length
       !Array to be nullified
       integer(kind=short), intent(inout) :: dx(length)
       integer                  :: i
       if (length < 0) then
          !do nothing
       else
          do i = 1, length
             dx(i) = shortzero
          enddo
       endif
     end subroutine ls_sizero

!> \brief Sets a short integer array of length *LENGTH* to shortzero.
!> \author T. Kjaergaard
!> \date May 5, 2011
subroutine ls_sicopy(length,dx,dy)
  use precision
  implicit none
  !Length of array
  integer(kind=long), intent(in)      :: length
  integer(kind=short), intent(in) :: dx(length)
  integer(kind=short), intent(inout):: dy(length)
  integer                  :: i
  do i = 1, length
     dy(i) = dx(i)
  enddo
end subroutine ls_sicopy

!> \brief Sets an integer array of length *LENGTH* to shortzero.
!> \author T. Kjaergaard
!> \date May 5, 2011
subroutine ls_icopy(length,dx,dy)
  use precision
  implicit none
  !Length of array
  integer(kind=long), intent(in)      :: length
  integer, intent(in) :: dx(length)
  integer, intent(inout):: dy(length)
  integer                  :: i
  do i = 1, length
     dy(i) = dx(i)
  enddo
end subroutine ls_icopy

!> \brief Sets a short integer array of length *LENGTH* to shortzero.
!> \author T. Kjaergaard
!> \date May 5, 2011
subroutine ls_dcopy(length,dx,dy)
  use precision
  implicit none
  !Length of array
  integer(kind=long), intent(in)      :: length
  real(realk), intent(in) :: dx(length)
  real(realk), intent(inout):: dy(length)
  integer                  :: i
  integer(kind=long)       :: i2
  do i2 = 1, length
     dy(i2) = dx(i2)
  enddo
end subroutine ls_dcopy

  !> \brief Sets a real array of length *LENGTH* to zero.
  !> \author T. Kjaergaard
  !> \date 2011
   subroutine ls_SetToOne(dx, length)
     use precision
     implicit none
     !Length of array
     integer, intent(in)      :: length
     !Array to be nullified
     real(realk), intent(out) :: dx(length)
     integer                  :: i
     
     if (length < 0) then
        !do nothing
     else
        do i = 1, length
           dx(i) = 1.0E0_realk
        enddo
     endif
   end subroutine ls_SetToOne

  !> \brief Sets an integer array of length *LENGTH* to zero.
  !> \author H. J. Aa. Jensen. F90'fied by S. Host
  !> \date May 5, 1984
  subroutine ls_izero(int, length)
  implicit none
       !Length of array
       integer, intent(in)  :: length
       !Array to be nullified
       integer, intent(out) :: int(length)
       integer              :: i

     if (length < 0) then
        !do nothing
     else
        do i = 1, length
           int(i) = 0
        enddo
     endif
  end subroutine ls_izero

  !> \brief End program with a text message
  !> \author T. Helgaker. F90'fied by S. Host
  !> \date 1984
  subroutine lsquit(text,lupri)
  use precision
#ifdef VAR_MPI
  use infpar_module
  use lsmpi_type
#endif
#ifdef VAR_IFORT
#ifndef VAR_INT64
  use IFCORE
#endif
#endif
  implicit none
      !> Text string to be printed
      CHARACTER(len=*), intent(in) :: TEXT
      !> Logical unit number for output
      integer, intent(in) :: lupri
      integer             :: luprin
      real(realk) :: CTOT,WTOT
      integer :: user_exit_code,qqstatus
!
!     Stamp date and time and hostname to output
!
      if (lupri > 0) then
         luprin = lupri

         CALL ls_TSTAMP('  --- SEVERE ERROR, PROGRAM WILL BE ABORTED ---', LUPRIN)
         WRITE (LUPRIN,'(/2A/)') ' Reason: ',TEXT
         CALL ls_TSTAMP('  --- SEVERE ERROR, PROGRAM WILL BE ABORTED ---', 6)
         WRITE (6,'(/2A/)') ' Reason: ',TEXT
      else
         luprin = 6
         CALL ls_TSTAMP('  --- SEVERE ERROR, PROGRAM WILL BE ABORTED ---', LUPRIN)
         WRITE (LUPRIN,'(/2A/)') ' Reason: ',TEXT
      endif


#if  defined (SYS_AIX) || defined (SYS_LINUX)
!     Write to stderr
      WRITE (0,'(/A/1X,A)') '  --- SEVERE ERROR, PROGRAM WILL BE ABORTED ---',TEXT
#endif
#ifdef VAR_IFORT
      qqstatus = -1
      user_exit_code = -1
#ifndef VAR_INT64
      CALL TRACEBACKQQ("TRACEBACKQQ INFO:",USER_EXIT_CODE,qqSTATUS)
#endif
#endif

      CALL ls_GETTIM(CTOT,WTOT)
      CALL ls_TIMTXT('>>>> Total CPU  time used in LSDALTON:',CTOT,LUPRIN)
      CALL ls_TIMTXT('>>>> Total wall time used in LSDALTON:',WTOT,LUPRIN)
      CALL ls_FLSHFO(LUPRIN)
!It may seem like a good idea to wake up the slaves so that the slaves can all call mpi_finalize and quit. However in the case of MPI the lsquit can be called in many ways.
! Option 1: The Master is awake and the slaves are sleeping. Master calls lsquit
!           Here it can make sense to wake up the slaves and have the slaves 
!           call mpi_finalize and quit
! Option 2: The Master and the Slaves are awake. Master calls lsquit
!           In this case Master should not broadcast a wake up call as the 
!           slaves are already sleeping - so it does not make sense 
!           and the calculation will hang in the MPI broadcast routine. 
!           The MPI slaves will wait in some reduction routine or something
!           While the Master i waiting in the bcast routine 
! Option 3: The Master and the Slaves are awake. Slave calls lsquit
!           Clearly it should not wake up the other slaves
!
!If master and slaves calls EXIT directly the mpiexec should kill the slaves!
!#ifdef VAR_MPI
!      IF(infpar%mynum.EQ.infpar%master)call lsmpi_finalize(lupri,.FALSE.)
!#endif
      !TRACEBACK INFO TO SEE WHERE IT CRASHED!!
#if defined (SYS_LINUX)
      CALL EXIT(100)
#else
      STOP 100
#endif
      end subroutine lsquit

      !> \brief Print Stack
      !> \author T. Kjaergaard
      !> \date Jan 2014
      subroutine LsTraceBack(text)
        use precision
#ifdef VAR_IFORT
#ifndef VAR_INT64
        use IFCORE
        implicit none
        !> Text string to be printed
        CHARACTER(len=*), intent(in) :: TEXT
        WRITE (*,'(/2A/)')"TRACEBACKQQ INFO:",TEXT
        CALL TRACEBACKQQ("TRACEBACKQQ INFO:",USER_EXIT_CODE=-1)
#else
        WRITE (*,'(/2A/)')"LsTraceBack do not work using -int64"
#endif
#else
        WRITE (*,'(/2A/)')"LsTraceBack do not work unless ifort is used"
#endif
      end subroutine LsTraceBack

      !> \brief Print a header. Based on HEADER by T. Helgaker
      !> \author S. Host
      !> \date June 2010
      subroutine LS_HEADER(HEAD,IN,lupri)
      implicit none
         CHARACTER(len=*), intent(in) :: HEAD
         integer, intent(in)          :: in
         integer, intent(in)          :: lupri
         integer                      :: length, indent, i

      LENGTH = LEN(HEAD)
      IF (IN > 0) THEN
         INDENT = IN + 1
      ELSE
         INDENT = (72 - LENGTH)/2 + 1
      END IF
      WRITE (LUPRI, '(//,80A)') (' ',I=1,INDENT), HEAD
      WRITE (LUPRI, '(80A)') (' ',I=1,INDENT), ('-',I=1,LENGTH)
      WRITE (LUPRI, '()')
      
      end subroutine LS_HEADER

  !> \brief To stamp as many as possible of text, date, time, computer, and hostname to LUPRIN
  !> \author H. J. Aa. Jensen. F90'fied by S. Host
  !> \date July 9, 1990
  subroutine ls_tstamp(TEXT,LUPRIN)
  use precision
  implicit none
      CHARACTER(*), intent(in) :: TEXT
      integer, intent(in)      :: luprin
      integer :: ltext
#if  defined (SYS_LINUX)
      CHARACTER(40) :: HSTNAM
      CHARACTER(24) :: FDATE
#endif
#if defined (SYS_AIX)
      CHARACTER(24) :: fdate
      CHARACTER(32) :: HSTNAM
#endif

      LTEXT = LEN(TEXT)
      IF (LTEXT .GT. 0) THEN
         WRITE (LUPRIN,'(/A)') TEXT
      ELSE
         WRITE (LUPRIN,'()')
      END IF

#if defined (SYS_LINUX)
      WRITE (LUPRIN,'(T6,2A)') 'Date and time (Linux)  : ',FDATE()
      CALL HOSTNM(HSTNAM)
      WRITE (LUPRIN,'(T6,2A)') 'Host name              : ',HSTNAM
#endif
#if defined (SYS_AIX)
! 930414-hjaaj: apparent error IBM's xlf library routines:
!    when 'T6' then column 1-5 not blanked but contains text
!    from a previous print statement!
!     AIX XL FORTRAN version 2.3+
      WRITE (LUPRIN,'(2A)') '     Date and time (IBM-AIX): ',fdate()
!     WRITE (LUPRIN,'(T6,2A)') 'Date and time (IBM-AIX): ',fdate_()
!     CALL hostnm_(HSTNAM)
      CALL hostnm(HSTNAM)
      WRITE (LUPRIN,'(2A)') '     Host name              : ',HSTNAM
!     WRITE (LUPRIN,'(T6,2A)') 'Host name              : ',HSTNAM
#endif
   end subroutine ls_tstamp

  !> \brief TIMTXT based on TIMER by TUH.
  !> \author H. J. Aa. Jensen. F90'fied by S. Host
  !> \date July 9, 1990
  subroutine ls_timtxt(TEXT,TIMUSD,LUPRIN)
   use precision
   implicit none 
      CHARACTER(*),intent(in) :: TEXT
      real(realk), intent(in) :: timusd
      integer, intent(in)     :: luprin
      CHARACTER :: AHOUR*6, ASEC*8, AMIN*8
      integer :: minute, isecnd, ihours

      ISECND = NINT(TIMUSD)
      IF (ISECND > 60) THEN
         MINUTE = ISECND/60
         IHOURS = MINUTE/60
         MINUTE = MINUTE - 60*IHOURS
         ISECND = ISECND - 3600*IHOURS - 60*MINUTE
         IF (IHOURS == 1) THEN
            AHOUR = ' hour '
         ELSE
            AHOUR = ' hours'
         END IF
         IF (MINUTE == 1) THEN
            AMIN = ' minute '
         ELSE
            AMIN = ' minutes'
         END IF
         IF (ISECND == 1) THEN
            ASEC = ' second '
         ELSE
            ASEC = ' seconds'
         END IF
         IF (IHOURS /= 0) THEN
            WRITE(LUPRIN,"(1X,A,I4,A,I3,A,I3,A)") TEXT, IHOURS, AHOUR, MINUTE, AMIN, ISECND, ASEC
         ELSE
            WRITE(LUPRIN,"(1X,A,     I3,A,I3,A)") TEXT, MINUTE, AMIN, ISECND, ASEC
         END IF
      ELSE
         WRITE(LUPRIN,"(1X,A,F7.2,' seconds')") TEXT,TIMUSD
      END IF

      CALL ls_FLSHFO(LUPRIN)
  end subroutine ls_timtxt

  !> \brief Return elapsed CPU time and elapsed real time.
  !> \author H. J. Aa. Jensen. F90'fied by S. Host
  !> \date December 18, 1984
  subroutine ls_gettim(cputime,walltime)
  use precision
  implicit none
      real(realk), intent(out) :: cputime, walltime

      real(realk),PARAMETER :: D0 = 0.0E0_realk
      logical    :: first = .true.
      real(realk), save :: TCPU0, twall0
      real(realk)       :: tcpu1, twall1
      integer        :: dateandtime0(8), dateandtime1(8)
      dateandtime0 = 0
      dateandtime1 = 0

      if (first) then
         first = .false.
         call cpu_time(TCPU0)
         call date_and_time(values=dateandtime0)
         call ls_get_walltime(dateandtime0,twall0)
      end if
      call cpu_time(tcpu1)
      call date_and_time(values=dateandtime1)
      call ls_get_walltime(dateandtime1,twall1)

      cputime = tcpu1 - TCPU0
      walltime = twall1 - twall0

!#if defined (SYS_AIX)
!!     SGI IRIX  etc. code
!      external time
!!     IRIX: otherwise "subroutine time", an intrinsic fu. will be used
!#endif
!      integer :: time
!      real    :: etime, tarray(2), timwl
!      logical, save :: first = .true.
!      real, save    :: wall0,TCPU0
!
!      if (first) then
!         first = .false.
!         TCPU0 = etime(tarray)
!         wall0 = time()
!      end if
!      TIMCPU = etime(tarray) - TCPU0
!      timwl = time() - wall0
!      TIMWAL = timwl
!#else
!#if (defined (SYS_LINUX) && defined (VAR_IFC))
!      integer      :: time,darr(9)
!      real         :: etime, tarray(2),err
!      logical,save :: first = .true.
!      real, save   :: TCPU0, wall0
!
!      if (first) then
!         first = .false.
!         TCPU0 = etime(tarray)
!          call gettimeofday(darr,err)
!          wall0 = darr(1)*1E0_realk + darr(2)/1E6_realk
!      end if
!      TIMCPU = etime(tarray) - TCPU0
!      call gettimeofday(darr,err)
!      TIMWAL = darr(1)*1E0_realk + darr(2)/1E6_realk - wall0
!#else
!      real(realk),PARAMETER :: D0 = 0.0E0_realk
!      logical    :: first = .true.
!      real, save :: TWALL0,TCPU0
!
!      if (first) then
!         first = .false.
!         TCPU0  = second()
!      end if
!      TIMCPU = SECOND() - TCPU0
!print *, 'first, TCPU0:', first, TCPU0
!! insert appropriate routine to get elapsed real time (wall time)
!! here
!      TIMWAL = D0
!#endif  /* VAR_IFC */
!#endif  /* SYS_AIX, VAR_ABSOFT etc. */
!
      end subroutine ls_gettim

!> \brief Get elapsed walltime in seconds since 1/1-2010 00:00:00
!> \author S. Host
!> \date October 2010
!>
!> Years that are evenly divisible by 4 are leap years. 
!> Exception: Years that are evenly divisible by 100 are not leap years, 
!> unless they are also evenly divisible by 400. Source: Wikipedia
!>
subroutine ls_get_walltime(dateandtime,walltime)
use precision
implicit none
   !> "values" output from fortran intrinsic subroutine date_and_time
   integer, intent(in) :: dateandtime(8)
   !> Elapsed wall time in seconds
   real(realk), intent(out) :: walltime
   integer                  :: month, year
   
! The output from the fortran intrinsic date_and_time
! gives the following values:
! 1. Year
! 2. Month
! 3. Day of the month
! 4. Time difference in minutes from Greenwich Mean Time (GMT)
! 5. Hour
! 6. Minute
! 7. Second
! 8. Millisecond

! Count seconds, minutes, hours, days, months and years and sum up seconds:
! We don't count milliseconds.

   walltime = 0.0E0_realk

   walltime = dble(dateandtime(8))/1.0E3_realk                !Seconds counted

   walltime = walltime + dble(dateandtime(7))                 !Seconds counted

   walltime = walltime + 60E0_realk*dateandtime(6)            !Minutes counted

   walltime = walltime + 3600E0_realk*dateandtime(5)          !Hours counted

   walltime = walltime + 24E0_realk*3600E0_realk*(dateandtime(3)-1) !Days counted (substract 1 to count only whole days)

   !Months are special, since they are not equally long:

   do month = 1, dateandtime(2)-1 !substract 1 to count only whole months
      if (month == 1 .or. month == 3 .or. month == 5 .or. month == 7 .or. &
       &  month == 8 .or. month == 10) then !Since we subtract 1, month can never be 12
         walltime = walltime + 31E0_realk*24E0_realk*3600E0_realk
      else if (month == 2) then
         if (.false.) then !insert exception for if current year is a leap year
            walltime = walltime + 29E0_realk*24E0_realk*3600E0_realk
         else
            walltime = walltime + 28E0_realk*24E0_realk*3600E0_realk
         endif
      else if (month == 4 .or. month == 6 .or. month == 9 .or. month == 11) then
         walltime = walltime + 30E0_realk*24E0_realk*3600E0_realk
      else
         stop 'Unknown month (ls_get_walltime)'
      endif
   enddo

   !Years are special, since leap years are one day longer:

   do year = 2010, dateandtime(1) 
      if (mod(year,400)==0) then
         walltime = walltime + 366*24*3600 !Leap year
      else if (mod(year,100)==0) then
         walltime = walltime + 365*24*3600 !Not leap year
      else if (mod(year,4)==0) then
         walltime = walltime + 366*24*3600 !Leap year
      else
         walltime = walltime + 365*24*3600 !Not leap year
      endif
   enddo

end subroutine ls_get_walltime

!> \brief Flush formatted output unit (empty buffers). 
!> \author H. J. Aa. Jensen. F90'fied by S. Host
!> \date February 10, 1989
!>
!> If no flush utility, flush is achieved by CLose and reOPen Formatted Output
!> 
subroutine ls_flshfo(IUNIT)
!
! Calls to this subroutine makes it possible to read the output
! up to the moment of the last call while the program continues
! executing (provided the computer allows shared access).
! This subroutine may be a dummy routine.
!
      IF (IUNIT .GE. 0) THEN
!     ... do not try to flush unused units (e.g. LUPRI on a slave) /hjaaj
#if  defined (SYS_AIX)  || defined (SYS_IRIX) || defined (SYS_NEC)  \
  || defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_LINUX) \
  || defined (SYS_SUN)  || defined (SYS_HAL) || defined (SYS_T90)   \
  || defined (SYS_HPUX) || defined (SYS_SX)
!
!        Force transfer of all buffered output to the file or device
!        associated with logical unit IUNIT.
!
         CALL FLUSH(IUNIT)
#endif
      END IF
end subroutine ls_flshfo

!> \brief OUTPUT PRINTS A REAL MATRIX IN FORMATTED FORM WITH NUMBERED ROWS AND COLUMNS
!> 
!> THE PROGRAM IS SET UP TO HANDLE 5 COLUMNS/PAGE WITH A 1P,5E24_realk.15 FORMAT
!> FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED,
!> CHANGE FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER
!> OF COLUMNS.
!>
!> \author NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
!>         FLORIDA, GAINESVILLE, FEBRUARY 26, 1971; 
!>         Revised 15-Dec-1983 by Hans Jorgen Aa. Jensen.
!>         16-Aug-2010 f90'fied by S. Host
!> \date August 2010
subroutine ls_output(AMATRX,ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,NCTL,LUPRI)
use precision
implicit none 
      !> ROW NUMBER AT WHICH OUTPUT IS TO BEGIN
      integer,intent(in) :: ROWLOW
      !> ROW NUMBER AT WHICH OUTPUT IS TO END
      integer,intent(in) :: ROWHI
      !> COLUMN NUMBER AT WHICH OUTPUT IS TO BEGIN
      integer,intent(in) :: COLLOW
      !> COLUMN NUMBER AT WHICH OUTPUT IS TO END
      integer,intent(in) :: COLHI
      !> ROW DIMENSION OF AMATRX(',')
      integer,intent(in) :: ROWDIM
      !> COLUMN DIMENSION OF AMATRX(',')
      integer,intent(in) :: COLDIM
      !> CARRIAGE CONTROL FLAG; (1 FOR SINGLE, 2 FOR DOUBLE, 3 FOR TRIPLE SPACE)
      integer,intent(in) :: NCTL
      !> Logical unit number for output file
      integer,intent(in) :: lupri
      !> MATRIX TO BE OUTPUT
      real(realk), intent(in) :: AMATRX(ROWDIM,COLDIM)
      integer :: BEGIN, KCOL, j, i, k, mctl, last
      CHARACTER(len=1) :: ASA(3), BLANK, CTL
      CHARACTER :: PFMT*20, COLUMN*8
      real(realk) :: amax
      real(realk), PARAMETER :: ZERO=0E00_realk
      integer, parameter     :: KCOLP=4, KCOLN=6
      real(realk), PARAMETER :: FFMIN=1E-3_realk, FFMAX = 1E3_realk
      DATA COLUMN/'Column  '/, BLANK/' '/, ASA/' ', '0', '-'/
 
      IF (ROWHI.LT.ROWLOW) RETURN
      IF (COLHI.LT.COLLOW) RETURN
 
      AMAX = ZERO
      DO J = COLLOW,COLHI
        DO I = ROWLOW,ROWHI
            AMAX = MAX( AMAX, ABS(AMATRX(I,J)) )
        enddo
      enddo
      IF (AMAX == ZERO) THEN
         WRITE (LUPRI,'(/T6,A)') 'Zero matrix.'
         RETURN
      END IF
      IF (FFMIN <= AMAX .AND. AMAX <= FFMAX) THEN
!        use F output format
         PFMT = '(A1,I7,2X,8F15.8)'
      ELSE
!        use 1PD output format
         PFMT = '(A1,I7,2X,1P,8D15.6)'
      END IF
!
      IF (NCTL .LT. 0) THEN
         KCOL = KCOLN
      ELSE
         KCOL = KCOLP
      END IF
      MCTL = ABS(NCTL)
      IF ((MCTL.LE. 3).AND.(MCTL.GT. 0)) THEN
         CTL = ASA(MCTL)
      ELSE
         CTL = BLANK
      END IF
!
      LAST = MIN(COLHI,COLLOW+KCOL-1)
      DO 2 BEGIN = COLLOW,COLHI,KCOL
         WRITE (LUPRI,1000) (COLUMN,I,I = BEGIN,LAST)
         DO 1 K = ROWLOW,ROWHI
            DO 4 I = BEGIN,LAST
               IF (AMATRX(K,I).NE.ZERO) GO TO 5
    4       CONTINUE
         GO TO 1
    5       WRITE (LUPRI,PFMT) CTL,K,(AMATRX(K,I), I = BEGIN,LAST)
    1    CONTINUE
    2 LAST = MIN(LAST+KCOL,COLHI)
 1000 FORMAT (/12X,6(3X,A6,I4,2X),(3X,A6,I4))
      END subroutine ls_output

      subroutine shortint_output(MAT2,dim1,dim2,lupri)
        implicit none
        integer             :: dim1,dim2,lupri,dimT2,dimT
        integer(kind=selected_int_kind(1)) :: MAT2(dim1,dim2)
        integer             :: I,J,K
        dimT=20*(dim2/20)
        DO J=1,dimT,20
           WRITE(lupri,'(A,20I4)')'Column:',&
                &J,J+1,J+2,J+3,J+4,J+5,J+6,J+7,J+8,J+9,&
                &J+10,J+11,J+12,J+13,J+14,J+15,J+16,J+17,J+18,J+19
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,20I4)')I,'   ',(MAT2(I,J+K),K=0,19)
           ENDDO
        ENDDO
        dimT2=dim2-dimT
        J=dimT
        IF(dimT2.EQ.1)THEN        
           WRITE(lupri,'(A,I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,I4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.2)THEN
           WRITE(lupri,'(A,2I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,2I4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.3)THEN
           WRITE(lupri,'(A,3I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,3I4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.4)THEN
           WRITE(lupri,'(A,4I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,4I4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.5)THEN
           WRITE(lupri,'(A,5I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,5I4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.6)THEN
           WRITE(lupri,'(A,6I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,6I4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.7)THEN
           WRITE(lupri,'(A,7I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,7I4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.8)THEN
           WRITE(lupri,'(A,8I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,8I4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.9)THEN
           WRITE(lupri,'(A,9I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,9I4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.10)THEN
           WRITE(lupri,'(A,10I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,10I4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.11)THEN
           WRITE(lupri,'(A,11I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,11I4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.12)THEN
           WRITE(lupri,'(A,12I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,12I4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.13)THEN
           WRITE(lupri,'(A,13I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,13I4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.14)THEN
           WRITE(lupri,'(A,14I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,14I4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.15)THEN
           WRITE(lupri,'(A,15I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,15I4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.16)THEN
           WRITE(lupri,'(A,16I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,16I4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.17)THEN
           WRITE(lupri,'(A,17I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,17I4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.18)THEN
           WRITE(lupri,'(A,18I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,18I4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.19)THEN
           WRITE(lupri,'(A,19I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,19I4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ENDIF
        
      end subroutine shortint_output

      subroutine loutput(MAT2,dim1,dim2,lupri)
        implicit none
        integer             :: dim1,dim2,lupri,dimT2,dimT
        logical             :: MAT2(dim1,dim2)
        integer             :: I,J,K
        dimT=20*(dim2/20)
        DO J=1,dimT,20
           WRITE(lupri,'(A,20I4)')'Column:',&
                &J,J+1,J+2,J+3,J+4,J+5,J+6,J+7,J+8,J+9,&
                &J+10,J+11,J+12,J+13,J+14,J+15,J+16,J+17,J+18,J+19
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,20L4)')I,'   ',(MAT2(I,J+K),K=0,19)
           ENDDO
        ENDDO
        dimT2=dim2-dimT
        J=dimT
        IF(dimT2.EQ.1)THEN        
           WRITE(lupri,'(A,I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,L4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.2)THEN
           WRITE(lupri,'(A,2I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,2L4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.3)THEN
           WRITE(lupri,'(A,3I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,3L4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.4)THEN
           WRITE(lupri,'(A,4I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,4L4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.5)THEN
           WRITE(lupri,'(A,5I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,5L4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.6)THEN
           WRITE(lupri,'(A,6I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,6L4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.7)THEN
           WRITE(lupri,'(A,7I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,7L4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.8)THEN
           WRITE(lupri,'(A,8I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,8L4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.9)THEN
           WRITE(lupri,'(A,9I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,9L4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.10)THEN
           WRITE(lupri,'(A,10I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,10L4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.11)THEN
           WRITE(lupri,'(A,11I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,11L4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.12)THEN
           WRITE(lupri,'(A,12I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,12L4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.13)THEN
           WRITE(lupri,'(A,13I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,13L4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.14)THEN
           WRITE(lupri,'(A,14I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,14L4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.15)THEN
           WRITE(lupri,'(A,15I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,15L4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.16)THEN
           WRITE(lupri,'(A,16I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,16L4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.17)THEN
           WRITE(lupri,'(A,17I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,17L4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.18)THEN
           WRITE(lupri,'(A,18I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,18L4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ELSEIF(dimT2.EQ.19)THEN
           WRITE(lupri,'(A,19I4)')'Column:',(J+K,K=1,dimT2)
WRITE(lupri,'(A)')'======================================================================================='
           DO I=1,dim1
              WRITE(lupri,'(I4,A3,19L4)')I,'   ',(MAT2(I,J+K),K=1,dimT2)
           ENDDO
        ENDIF
        
      end subroutine loutput

      subroutine ls_transpose(array_in,array_out,ndim)
use precision
        implicit none
        !>  the dimension
        integer, intent(in) :: ndim
        !> array to be reordered
        real(realk),intent(in) :: array_in(ndim,ndim)
        !> reordered array
        real(realk),intent(inout) :: array_out(ndim,ndim)
        integer :: a,b,ba,bb,da2,bcntr
        logical :: moda,modb
        ! assuming available cache memory is 8 MB and we need to store two blocks at a time 
        ! bs = 700 !int(SQRT((8000.0*1000.0)/(8.0*2.0))) roughly equal to 700
        integer,parameter :: bs=700
!!$        da2=(ndim/bs)*bs
!!$        moda=(mod(ndim,bs)>0)
!!$!        bcntr=bs-1    
!!$        bcntr=699
!!$        !$OMP PARALLEL DEFAULT(NONE),PRIVATE(a,b,ba,bb),SHARED(array_in,array_out,&
!!$        !$OMP& da2,bcntr,moda,ndim)
!!$        if (da2 .gt. 0) then
!!$           print*,'option1'
!!$           !$OMP DO
!!$           do ba=1,da2,bs
!!$              do bb=1,da2,bs
!!$                 do a=0,bcntr
!!$                    do b=0,bcntr
!!$                       array_out(bb+b,ba+a)=array_in(ba+a,bb+b)
!!$                    enddo
!!$                 enddo
!!$              enddo
!!$           enddo
!!$           !$OMP END DO NOWAIT
!!$        endif
!!$        if (moda .and. da2 .gt. 0) then
!!$           print*,'option2'
!!$           !$OMP DO
!!$           do bb=1,da2,bs
!!$              do a=da2+1,ndim
!!$                 do b=0,bcntr
!!$                    array_out(bb+b,a)=array_in(a,bb+b)
!!$                 enddo
!!$              enddo
!!$           enddo
!!$           !$OMP END DO NOWAIT
!!$           !$OMP DO
!!$           do ba=1,da2,bs     
!!$              do b=da2+1,ndim
!!$                 do a=0,bcntr
!!$                    array_out(b,ba+a)=array_in(ba+a,b)
!!$                 enddo
!!$              enddo
!!$           enddo
!!$           !$OMP END DO NOWAIT
!!$        endif
!!$        if (moda) then
!!$           print*,'option3'
!!$           !$OMP DO
!!$           do a=da2+1,ndim
!!$              do b=da2+1,ndim
!!$                 array_out(b,a)=array_in(a,b)
!!$!$OMP CRITICAL
!!$                 print*,'array_out(',b,',',a,')'
!!$!$OMP END CRITICAL
!!$              enddo
!!$           enddo
!!$           !$OMP END DO NOWAIT
!!$        endif
!!$        !$OMP END PARALLEL

           do a=1,ndim
              do b=1,ndim
                 array_out(b,a)=array_in(a,b)
              enddo
           enddo


      end subroutine ls_transpose


  module ls_env
    use precision

    contains
    subroutine ls_getenv(varname,leng,output_string,output_int,output_real,output_bool)
      implicit none
      integer, intent(in) :: leng
      character(leng),intent(in) :: varname
      character(*),intent(inout), optional :: output_string
      integer,intent(inout), optional      :: output_int
      real(realk),intent(inout),optional   :: output_real
      logical,intent(inout),optional       :: output_bool
      character(120) :: val
      integer :: stat,reqlen
      logical :: trim_var 
      stat = 0
      trim_var = .true.
      
      !this call is fortran 2003 standard and might not be present everywhere
      call GET_ENVIRONMENT_VARIABLE(varname,VALUE=val,LENGTH=reqlen,STATUS=stat,TRIM_NAME=trim_var)

      if(stat/=0.and.stat/=1)then
        if(stat==-1)then
          print *,"ERROR(ls_getenv): the internal val is too short to hold the&
          & value of the variable, please adjust the length and recompile"
          stop 0
        else
          print *,"ERROR(ls_getenv):error in GET_ENVIRONMENT_VARIABLE, status:",stat
          stop 0
        endif
      endif

      !if stat==1 then the variable is just not set and we skip overwriting the
      !value, if stat==0 then we can read it and overwrite
      if(stat==0)then 
        if(present(output_string))then
          output_string=TRIM(val)
        endif
        if(present(output_int))then
          read(val,*) output_int
        endif
        if(present(output_real))then
          read(val,*) output_real
        endif
        if(present(output_bool))then
          read(val,*) stat
          if(stat==1)then
            output_bool=.true.
          elseif(stat==0)then
            output_bool=.false.
          else
            print *,"ERROR(ls_getenv):the value of the variable you try to read&
            & is clearly not of boolean(logical) character, value:",stat,val
            stop 0
          endif
        endif
      endif
  
    end subroutine ls_getenv
  end module ls_env
