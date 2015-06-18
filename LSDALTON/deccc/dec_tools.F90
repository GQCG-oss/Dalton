!this file contains the lower-order tools as opposed to the dec_utils which
!conains higher routines
module dec_tools_module
  use precision
  use dec_typedef_module
  use fundamental
  use memory_handling
  use files
  use BUILDAOBATCH  

  contains

  !> \brief Print all elements of four-dimensional array to LSDALTON.OUT.
  !> Only to be used for testing purposes!
  !> \author Kasper Kristensen
  !> \date November 2011
  subroutine print_4_dimensional_array(dims,A,label)

    implicit none
    !> Dimensions of 4-dimensional array
    integer,dimension(4),intent(in) :: dims
    !> 4-dimensional array to be printed
    real(realk),intent(in) :: A(dims(1),dims(2),dims(3),dims(4))
    !> Label for array
    character(*), intent(in) :: label
    integer :: i,j,k,l

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '***********************************************************'
    write(DECinfo%output,*) '             ARRAY LABEL: ', label
    write(DECinfo%output,*) '***********************************************************'
    write(DECinfo%output,*)
    write(DECinfo%output,'(8X,a,8X,a,8X,a,8X,a,12X,a)') 'i','j','k','l', 'value'

    do i=1,dims(1)
       do j=1,dims(2)
          do k=1,dims(3)
             do l=1,dims(4)
                write(DECinfo%output,'(4i9,5X,g18.10)') i,j,k,l,A(i,j,k,l)
             end do
          end do
       end do
    end do

    write(DECinfo%output,*)
    write(DECinfo%output,*)


  end subroutine print_4_dimensional_array


  !> \brief Solve eigenvalue problem: F*C = S*C*eival
  !> \author Kasper Kristensen
  !> \date February 2011
  subroutine solve_eigenvalue_problem(n,F,S,eival,C)
    implicit none

    !> Dimension of matrices (Number of eigenvalues)
    integer,intent(in) :: n
    !> F matrix in eigenvalue problem (typically Fock matrix)
    real(realk),intent(in) :: F(n,n)
    !> Overlap matrix
    real(realk),intent(in) :: S(n,n)
    !> Eigenvectors
    real(realk),intent(inout) :: C(n,n)
    real(realk),intent(inout) :: eival(n)
    real(realk), pointer :: tmp(:,:)

    ! We must use temporary matrix as overlap matrix input because it is overwritten
    call mem_alloc(tmp,n,n)
    tmp(1:n,1:n) = S(1:n,1:n)

    ! Copy F elements to C, then use C as "F input".
    ! At the end, C is then overwritten by the eigenvectors.
    C(1:n,1:n) = F(1:n,1:n)

    ! Solve eigenvalue problem
    call my_DSYGV(N,C,tmp,eival,"DEC_SOLVE_EIGENVALUE")

    ! Free stuff
    call mem_dealloc(tmp)

  end subroutine solve_eigenvalue_problem



  !> \brief Solve eigenvalue problem: F*C = C*eival   (overlap matrix is the unit matrix)
  !> \author Kasper Kristensen
  !> \date February 2011
  subroutine solve_eigenvalue_problem_unitoverlap(n,F,eival,C)
    implicit none

    !> Dimension of matrices (Number of eigenvalues)
    integer,intent(in) :: n
    !> F matrix in eigenvalue problem (typically Fock matrix)
    real(realk),intent(in) :: F(n,n)
    !> Eigenvectors
    real(realk),intent(inout) :: C(n,n)
    real(realk),intent(inout) :: eival(n)
    real(realk), pointer :: tmp(:,:)
    integer :: i

    ! Overlap matrix is unit matrix
    call mem_alloc(tmp,n,n)
    tmp = 0.0_realk
    do i=1,n
       tmp(i,i) = 1.0_realk
    end do

    ! Copy F elements to C, then use C as "F input".
    ! At the end, C is then overwritten by the eigenvectors.
    C(1:n,1:n) = F(1:n,1:n)

    ! Solve eigenvalue problem
    call my_DSYGV(n,C,tmp,eival,"DEC_SOLVE_EIGENVALU2")

    ! Free stuff
    call mem_dealloc(tmp)

  end subroutine solve_eigenvalue_problem_unitoverlap


  !> \brief Solve nonsymmetric eigenvalue problem: 
  !> A R = eival R   
  !> L A = L eival
  !> where A is a nonsymmetrix matrix and R and L are the left and right
  !> eigenvectors, respectively.
  !> Note: Currently we quit if the eigenvalues are complex. This could be generalized.
  !> \author Kasper Kristensen
  !> \date February 2011
  subroutine solve_nonsymmetric_eigenvalue_problem_unitoverlap(n,A,eival,R,L)
    implicit none

    !> Dimension of matrices (Number of eigenvalues)
    integer,intent(in) :: n
    !> A matrix in eigenvalue problem 
    real(realk),intent(in) :: A(n,n)
    !> Right (R) and left (L) eigenvectors
    real(realk),intent(inout) :: R(n,n), L(n,n)
    !> Eigenvalues
    real(realk),intent(inout) :: eival(n)
    real(realk), pointer :: B(:,:),Atmp(:,:)
    integer :: i,j,lwork,info
    real(realk),pointer :: alphaR(:), alphaI(:), beta(:),work(:)
    real(realk),parameter :: thr=1.0e-9_realk

    ! Overlap matrix is unit matrix
    call mem_alloc(B,n,n)
    B = 0.0_realk
    do i=1,n
       B(i,i) = 1.0_realk
    end do

    ! Copy A matrix to avoid overwriting it
    call mem_alloc(Atmp,n,n)
    do j=1,n
       do i=1,n
          Atmp(i,j) = A(i,j)
       end do
    end do

    ! Allocate stuff
    call mem_alloc(alphaR,n)
    call mem_alloc(alphaI,n)
    call mem_alloc(beta,n)
    call mem_alloc(work,1)

    ! Determine optimal workspace
    lwork=-1
    info=0
    call DGGEV('V', 'V', n, Atmp, n, B, n, ALPHAR, ALPHAI,&
         & BETA, L, n, R, n, WORK, LWORK, INFO )
    lwork = int(work(1))

    ! Allocate work space
    call mem_dealloc(work)
    call mem_alloc(work,lwork)

    ! Solve eigenvalue problem
    call DGGEV('V', 'V', n, Atmp, n, B, n, ALPHAR, ALPHAI,&
         & BETA, L, n, R, n, WORK, LWORK, INFO )

    if(info/=0) then
       print *, 'INFO = ', INFO
       call lsquit('solve_nonsymmetric_eigenvalue_problem_unitoverlap: &
            & Error in DGGEV!',-1)
    end if

    ! Check that eigenvalues are real and well-defined
    do i=1,n
       if( (abs(alphaI(i)) > thr) .or. (abs(beta(i))<thr )  ) then
          print *, 'Eigenvalue number ',i 
          print *, 'Eigenvalue, real part ',alphaR(i)
          print *, 'Eigenvalue, imag part ',alphaI(i)
          print *, 'Division factor beta  ',beta(i)
          call lsquit('solve_nonsymmetric_eigenvalue_problem_unitoverlap: &
               & Complex or ill-defined eigenvalues!',-1)          
       end if
       eival(i) = alphaR(i) / beta(i)
    end do

    ! Free stuff
    call mem_dealloc(alphaR)
    call mem_dealloc(alphaI)
    call mem_dealloc(beta)
    call mem_dealloc(work)
    call mem_dealloc(B)
    call mem_dealloc(Atmp)

  end subroutine solve_nonsymmetric_eigenvalue_problem_unitoverlap


  subroutine init_batch_info(mylsitem,batch,max_allowed,nb)
     implicit none
     type(lsitem), intent(inout) :: mylsitem
     type(int_batch),intent(inout) :: batch
     integer, intent(in)  :: max_allowed,nb
     integer :: i,iorb,k
     logical :: master
     master = .true.
#ifdef VAR_MPI
     master = (infpar%lg_mynum == infpar%master)
#endif

     ! Orbital to batch information
     ! ----------------------------
     call mem_alloc(batch%orb2batch,nb)
     call build_batchesofAOS(DECinfo%output,mylsitem%setting,max_allowed,&
        & nb,batch%max_dim,batch%batchsize,batch%batchdim,batch%batchindex,&
        &batch%nbatches,batch%orb2batch,'R')

     ! Translate batchindex to orbital index
     ! -------------------------------------
     call mem_alloc(batch%batch2orb,batch%nbatches)
     do i=1,batch%nbatches
        call mem_alloc(batch%batch2orb(i)%orbindex,batch%batchdim(i))
        batch%batch2orb(i)%orbindex = 0
        batch%batch2orb(i)%norbindex = 0
     end do
     do iorb=1,nb
        i = batch%orb2batch(iorb)
        batch%batch2orb(i)%norbindex = batch%batch2orb(i)%norbindex+1
        K = batch%batch2orb(i)%norbindex
        batch%batch2orb(i)%orbindex(K) = iorb
     end do
  end subroutine init_batch_info

  subroutine free_batch_info(batch)
     type(int_batch),intent(inout) :: batch
     integer :: i
     call mem_dealloc(batch%orb2batch)
     call mem_dealloc(batch%batchdim)
     call mem_dealloc(batch%batchsize)
     call mem_dealloc(batch%batchindex)
     do i=1,batch%nbatches
        call mem_dealloc(batch%batch2orb(i)%orbindex)
        batch%batch2orb(i)%orbindex => null()
     end do
     call mem_dealloc(batch%batch2orb)
  end subroutine free_batch_info
end module dec_tools_module
