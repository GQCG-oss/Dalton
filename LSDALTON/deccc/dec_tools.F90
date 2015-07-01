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
  !> Copied (and modified slightly) from dggev.f reference file (Google it!): 
  !> If the j-th eigenvalue is real, then 
  !> the left eigenvectors u(j) are stored one
  !> after another in the columns of L, in the same order as
  !> their eigenvalues.If the j-th eigenvalue is real, then
  !> u(j) = L(:,j), the j-th column of L. If the j-th and
  !> (j+1)-th eigenvalues form a complex conjugate pair, then
  !> u(j) = L(:,j)+i*L(:,j+1) and u(j+1) = L(:,j)-i*L(:,j+1).
  !> Each eigenvector is scaled so the largest component has
  !> abs(real part)+abs(imag. part)=1.
  !> (and similarly for right eigenvectors R)
  !> \author Kasper Kristensen
  !> \date February 2011
  subroutine solve_nonsymmetric_eigenvalue_problem_unitoverlap(n,A,eivalREAL,eivalIMAG,R,L)
    implicit none

    !> Dimension of matrices (Number of eigenvalues)
    integer,intent(in) :: n
    !> A matrix in eigenvalue problem 
    real(realk),intent(in) :: A(n,n)
    !> Real and imaginary components of eigenvalues
    real(realk),intent(inout) :: eivalREAL(n),eivalIMAG(n)
    !> Right (R) and left (L) eigenvectors
    !> The right/left eigenvectors are stored one after another in columns of R/L.
    real(realk),intent(inout) :: R(n,n), L(n,n)
    real(realk), pointer :: B(:,:),Atmp(:,:),lambdaREAL(:),lambdaIMAG(:)
    integer :: i,j,lwork,info,idx
    real(realk),pointer :: work(:),Rtmp(:,:),Ltmp(:,:)
    real(realk),parameter :: thr=1.0e-10_realk
    integer,pointer :: tracklist(:)

    ! Initialization
    ! --------------

    ! Overlap matrix is assumed to be a unit matrix
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
    call mem_alloc(tracklist,n)
    call mem_alloc(work,1)
    call mem_alloc(Rtmp,n,n)
    call mem_alloc(Ltmp,n,n)
    call mem_alloc(lambdaREAL,n)
    call mem_alloc(lambdaIMAG,n)



    ! Solving eigenvalue problem
    ! --------------------------

    ! Determine optimal workspace
    lwork=-1
    info=0
    call DGEEV('V', 'V', n, Atmp, n, lambdaREAL, lambdaIMAG,&
         &  Ltmp, n, Rtmp, n, WORK, LWORK, INFO )

    lwork = int(work(1))

    if(info/=0) then
       print *, 'INFO = ', INFO
       call lsquit('solve_nonsymmetric_eigenvalue_problem_unitoverlap: &
            & Error1 in DGGEV!',-1)
    end if

    ! Allocate work space
    call mem_dealloc(work)
    call mem_alloc(work,lwork)

    ! Solve eigenvalue problem
    call DGEEV('V', 'V', n, Atmp, n, lambdaREAL, lambdaIMAG,&
         &  Ltmp, n, Rtmp, n, WORK, LWORK, INFO )


    if(info/=0) then
       print *, 'INFO = ', INFO
       call lsquit('solve_nonsymmetric_eigenvalue_problem_unitoverlap: &
            & Error2 in DGGEV!',-1)
    end if



    ! Sort eigenvalues according to size of real part (smallest first)
    ! ----------------------------------------------------------------

    ! Order eigenvalues with largest ones first
    call real_inv_sort_with_tracking(lambdaREAL,tracklist,n)

    ! Now we want the smallest first so we take the opposite order,
    ! we also set the output eigenvectors accordingly.
    do i=1,n
       eivalREAL(i) = lambdaREAL(n-i+1)

       ! Index in original eigenvalue order
       idx=tracklist(n-i+1)

       ! Set imaginary eigenvalue and eigenvectors according to new order
       eivalIMAG(i) = lambdaIMAG(idx)
       do j=1,n
          R(j,i) = Rtmp(j,idx)
          L(j,i) = Ltmp(j,idx)
       end do
    end do


    ! Free stuff
    call mem_dealloc(lambdaREAL)
    call mem_dealloc(lambdaIMAG)
    call mem_dealloc(tracklist)
    call mem_dealloc(work)
    call mem_dealloc(B)
    call mem_dealloc(Atmp)
    call mem_dealloc(Rtmp)
    call mem_dealloc(Ltmp)


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


  !> \brief Sort first 'n' elements of vector a with 'm' elements
  subroutine int_sort(a,n,m)

    implicit none
    integer, dimension(m), intent(inout) :: a
    integer, intent(in) :: m,n
    integer :: i
    integer :: tmp
    logical :: swp

    swp=.true.
    do while (swp)
       swp=.false.
       do i=1,n-1
          if(a(i)>a(i+1)) then

             tmp=a(i+1)
             a(i+1)=a(i)
             a(i)=tmp

             swp=.true.
          endif
       end do
    end do

    return
  end subroutine int_sort


  !> \brief Sort real vector, keeping track of the original indices.
  !> Note: Largest elements first!
  subroutine real_inv_sort_with_tracking(to_sort,to_track,n)

    implicit none
    integer, intent(in) :: n
    real(realk), dimension(n), intent(inout) :: to_sort
    integer, dimension(n), intent(inout) :: to_track
    real(realk) :: tmp
    integer :: tmp1,i
    logical :: swp

    ! Set original track order
    do i=1,n
       to_track(i)=i
    end do

    swp=.true.
    do while (swp)
       swp=.false.
       do i=1,n-1
          if(to_sort(i) < to_sort(i+1)) then ! reverse order

             tmp = to_sort(i+1)
             to_sort(i+1) = to_sort(i)
             to_sort(i) = tmp

             tmp1 = to_track(i+1)
             to_track(i+1) = to_track(i)
             to_track(i) = tmp1

             swp=.true.
          end if
       end do
    end do
    return
  end subroutine real_inv_sort_with_tracking


  !> \brief Sort integer vector, keeping track of the original indices.
  !> Note: Largest elements first!
  !> \author Kasper Kristensen (based on real_inv_sort_with_tracking)
  !> \date February 2011
  subroutine integer_inv_sort_with_tracking(to_sort,to_track,n)

    implicit none
    !> Dimension of vector to sort
    integer, intent(in) :: n
    !> Vector to sort
    integer, dimension(n), intent(inout) :: to_sort
    !> List of sorted original indices
    integer, dimension(n), intent(inout) :: to_track
    integer :: tmp,tmp1,i
    logical :: swp

    ! Set original track order
    do i=1,n
       to_track(i)=i
    end do

    swp=.true.
    do while (swp)
       swp=.false.
       do i=1,n-1
          if(to_sort(i) < to_sort(i+1)) then ! reverse order

             tmp = to_sort(i+1)
             to_sort(i+1) = to_sort(i)
             to_sort(i) = tmp

             tmp1 = to_track(i+1)
             to_track(i+1) = to_track(i)
             to_track(i) = tmp1

             swp=.true.
          end if
       end do
    end do

  end subroutine integer_inv_sort_with_tracking


  !> \brief Get list of atoms sorted according to distance from "MyAtom".
  !> \author Ida-Marie Hoeyvik
  subroutine GetSortedList(ListMyAtom,ListTrack,ToSort,&
       & n1,natoms,MyAtom)
    implicit none
    integer,intent(in)     :: n1,natoms,MyAtom
    real(realk),intent(in) :: ToSort(n1,natoms)
    real(realk)            :: ListMyAtom(n1), TempList(n1)
    integer                :: ListTrack(n1), TempTrack(n1)
    integer                :: counter,i

    ListMyAtom(:)=ToSort(:,MyAtom)
    ! Sort large--> small
    call real_inv_sort_with_tracking(ListMyAtom,ListTrack,n1)

    TempList = 0.0E0_realk
    TempTrack = 0
    counter = 1

    ! change to small-->large
    do i=n1,1,-1
       TempList(counter) = ListMyAtom(i)
       TempTrack(counter)= ListTrack(i)
       counter = counter + 1
    end do

    ListMyAtom=TempList
    ListTrack=TempTrack

  end subroutine GetSortedList




end module dec_tools_module
