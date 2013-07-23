!> @file
!> Contains routines needed for CROP solver in CC driver
!> \author Marcin Ziolkowski

!> Set of subroutines for DIIS or CROP solver for linear or nonlinear set of
!> equations in CC methods.
module crop

    use precision

  contains

  !> \brief Solve system of linear equations using DGESV lapack routine 
  !> \author Marcin Ziolkowski
  !> \param Amat Square matrix A in Ax=b
  !> \param Bvec RHS vector (b)
  !> \param n Size of the problem (A(n,n))
  subroutine SolveLinearEquations(Amat,Bvec,n)
  
    implicit none
    real(realk), dimension(n,n), intent(in) :: Amat
    real(realk), dimension(n), intent(inout) :: Bvec
    integer, dimension(n) :: ipiv
    integer :: n,infoLAPACK
    external dgesv

    infoLAPACK = 0
  
    call dgesv(n,1,Amat,n,ipiv,Bvec,n,infoLAPACK)

    if(infoLAPACK /= 0) then
       print *,"DGESV INFO:",infoLAPACK
       call lsquit('SolveLinearEquations: error in LAPACK DGESV routine',-1)
    end if
  
  end subroutine SolveLinearEquations

  !> \brief Solve CROP/DIIS system of equations to get mixing coefficients 
  !> \author Marcin Ziolkowski
  !> \param maxDIIS Maximum size of iterative subspace
  !> \param maxIter Maximum number of iterations
  !> \param iter Current iteration number
  !> \param Bmat B matrix in CROP/DIIS
  !> \param c Mixing coefficients
  !> \param verbose Verbose parameter
  subroutine CalculateDIIScoefficients(maxDIIS,maxIter,iter,Bmat,c,verbose)
  
    implicit none
    integer, intent(in) :: maxDIIS,maxIter,iter
    integer :: i,j,k,l
    real(realk), dimension(maxIter,maxIter), intent(in) :: Bmat
    real(realk), dimension(maxIter), intent(inout) :: c
    logical, intent(in) :: verbose
  
    ! Maximal size of the subspace equations should be 4x4 (3 from vectors and 1
    ! from Lagrange coefficient). 
    real(realk), dimension(min(maxDIIS+1,iter+1),min(maxDIIS+1,iter+1)) :: bb
    real(realk), dimension(min(maxDIIS+1,iter+1)) :: cc
  
    bb=1.0E0_realk
  
    k=1
    do i=max(iter-maxDIIS+1,1),iter
      l=1
      do j=max(iter-maxDIIS+1,1),iter
        bb(k,l)=Bmat(i,j)
        l=l+1
      end do 
      k=k+1
    end do
   
    bb(min(maxDIIS+1,iter+1),min(maxDIIS+1,iter+1))=0E0_realk
    cc=0E0_realk; cc(min(maxDIIS+1,iter+1))=1E0_realk
   
    if(verbose) call PrintMatrix(bb,min(maxDIIS+1,iter+1),'B matrix  ')
   
    ! --- solve system of the linear equations ---   
    call SolveLinearEquations(bb,cc,min(maxDIIS+1,iter+1))
    ! ---

    k=1
    c=0E0_realk
    do i=max(iter-maxDIIS+1,1),iter
      c(i)=cc(k)
      k=k+1
    end do

    return
  end subroutine CalculateDIIScoefficients

  !> \brief Print simple fortran matrix in nice form, just for debuging purpose 
  !> \author Marcin Ziolkowski
  !> \param my_matrix Matrix to print
  !> \param m Size of the matrix
  !> \param title Title for the matrix
  subroutine PrintMatrix(my_matrix,m,title)
  
    implicit none
    real(realk), dimension(m,m) :: my_matrix
    integer :: m,i,j
    character(len=10), intent(in), optional :: title
  
    print *,''
    if(present(title)) print *,' *** ',title,' ***'
    do i=1,m
      print 100,(my_matrix(i,j),j=1,m)
    enddo
    print *,''
  
100 format(15f10.5) 
  
    return
  end subroutine PrintMatrix

end module crop
