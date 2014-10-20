!> @file
!> Contains linear-equation solver routines used for overlap density-fitting
MODULE linsolvdf
use precision
use TYPEDEF  
use matrix_operations
use ls_Integral_Interface
use LSparameters
use lstiming
CONTAINS

!> \brief  Driver for the density-fitting linear-equation solver, A c = b. Includes the calculation
!>         of the metric matrix (A)   
!> \author S. Reine
!> \date   2008
!> \param  calpha    The density-fitting expansion coefficients (c)
!> \param  rhs       The right-hand side vector of the linear-equations (b)
!> \param  basisType Specifying the basis set to employ
!> \param  oper      Specifying the operator used for the linear-equations (Default is Coulomb)
!> \param  naux      The dimension of the basis-set employed
!> \param  ndmat     The number of rhs-vectors (and fitting coefficient vectors to be returned)
!> \param  SETTING   Settings for integral evaluation, used for calculation of the metric matrix (A)
!> \param  LUPRI     Print unit number
!> \param  LUERR     Error unit number
SUBROUTINE linsolv_df(calpha,rhs,basisType,oper,naux,ndmat,SETTING,LUPRI,LUERR)
implicit none
TYPE(LSSETTING)          :: SETTING
INTEGER                  :: LUPRI,LUERR,naux,ndmat
real(realk)              :: calpha(naux*ndmat)
real(realk)              :: rhs(naux*ndmat)
integer                  :: basisType,oper
TYPE(MATRIX),target      :: alphabeta2
TYPE(MATRIXP)            :: Intmat(1)

logical                  :: iterative_minushalf
character*80             :: filename
real(realk),pointer      :: alphabeta(:,:,:,:,:)
real(realk),pointer      :: work(:)
integer                  :: info,usemat
!
integer :: i1,i2,n,nsig(12)
logical :: IntegralTransformGC
!INTERFACE
!   SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
!     CHARACTER     ::     UPLO
!     INTEGER       ::   INFO, LDA, N
!     DOUBLE PRECISION  :: A( LDA, * )
!   END SUBROUTINE DPOTRF
!END INTERFACE
IntegralTransformGC = .FALSE.
iterative_minushalf=.false.
!set threshold
SETTING%SCHEME%intTHRESHOLD = SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
call Test_if_64bit_integer_required(naux,naux)
call mem_alloc(alphabeta,naux,1,naux,1,1)
! Solver based on a Lowdin decomposition, c = A^{-0.5} A^{-0.5} b
IF(ITERATIVE_MINUSHALF) THEN
    call io_get_filename(filename,'ALBEminushalf',basisType,AOEmpty,basisType,AOEmpty,0,0,&
     &                   oper,Contractedinttype,.FALSE.,LUPRI,LUERR)
    if(.not.io_file_exist(filename,SETTING%IO)) then
!      Caclulate matric matrix A
       call initIntegralOutputDims(setting%Output,naux,1,naux,1,1)
       call ls_getIntegrals(basisType,AOEmpty,basisType,AOEmpty,oper,RegularSpec,Contractedinttype,&
                             SETTING,LUPRI,LUERR)
       call retrieve_Output(lupri,setting,alphabeta,IntegralTransformGC)
!      Get decomposition A^{-0.5}
       CALL matrix_minushalf_df(alphabeta,naux,SETTING,LUPRI)
!      Write decomposed matrix to file
       call io_add_filename(SETTING%IO,filename,LUPRI)
       call io_write(alphabeta,naux,naux,1,1,1,filename,SETTING%IO,LUPRI,LUERR)
    else
!      Retrieve decomposed matrix (A^{-0.5}) from file
       call io_read(alphabeta,naux,naux,1,1,1,filename,SETTING%IO,LUPRI,LUERR)
    endif
! Solve linear-equation system using A^{-0.5} matrix twice to reduce errors
    call mem_alloc(work,naux*ndmat)
    CALL DGEMM('N','N',naux,ndmat,naux,1E0_realk,alphabeta,naux,rhs,naux,0E0_realk,work,naux)
    CALL DGEMM('N','N',naux,ndmat,naux,1E0_realk,alphabeta,naux,work,naux,0E0_realk,calpha,naux)
    call mem_dealloc(work)

! Solve linear equations using standard library routines based on Cholesky-factorization
ELSE

    call io_get_filename(filename,'ALBElu',basisType,AOEmpty,basisType,AOEmpty,0,0,&
     &                   oper,Contractedinttype,.FALSE.,LUPRI,LUERR)
    if(.not.io_file_exist(filename,SETTING%IO)) then
!       Caclulate matric matrix A
        call initIntegralOutputDims(setting%Output,naux,1,naux,1,1)
        call ls_getIntegrals(basisType,AOEmpty,basisType,AOEmpty,oper,RegularSpec,Contractedinttype,&
                             SETTING,LUPRI,LUERR)
        call retrieve_Output(lupri,setting,alphabeta,IntegralTransformGC)
        info = 0
!       Perform Cholesky factors C from matric matrix A (A = C^T C)
        CALL DPOTRF('U',naux,alphabeta,naux,info)
        IF (info.ne. 0) THEN
           WRITE(LUPRI,'(1X,A,I5)') 'DPOTRF error in linsolv_df. Info =',info
           CALL LSQUIT('DPOTRF error in linsolv_df',lupri)
        ENDIF
!       Write Cholesky factor C to file
        call io_add_filename(SETTING%IO,filename,LUPRI)
        call io_write(alphabeta,naux,naux,1,1,1,filename,SETTING%IO,LUPRI,LUERR)
    else
!       Retrieve Cholesky factor C from file
        call io_read(alphabeta,naux,naux,1,1,1,filename,SETTING%IO,LUPRI,LUERR)
    endif

!   Solve linear-equation system using the Cholesky factor and forward-backward substitution
    call DCOPY(naux*ndmat,rhs,1,calpha,1)
    CALL DPOTRS('U',naux,ndmat,alphabeta,naux,calpha,naux,info)
    IF (info.ne. 0) THEN
       WRITE(LUPRI,'(1X,A,I5)') 'DPOTRS error in linsolv_df. Info =',info
       CALL LSQUIT('DPOTRS error in linsolv_df',lupri)
    ENDIF
    
endif

call mem_dealloc(alphabeta)

END SUBROUTINE linsolv_df

!> \brief  Lowdin decomposition interface
!> \author S. Reine
!> \date   2008
!> \param  A         The square matrix to decompose on entry, the decomposed matrix on exit
!> \param  n         The matrix dimensions
!> \param  SETTING   Settings for integral evaluation, used for calculation of the metric matrix (A)
!> \param  LUPRI     Print unit number
SUBROUTINE matrix_minushalf_df(A,n,SETTING,LUPRI)
implicit none
TYPE(LSSETTING) :: SETTING
INTEGER         :: LUPRI
integer         :: n
real(realk)     :: A(n,n)
!
TYPE(matrix)    :: S, S_sqrt, S_minus_sqrt
real(realk)     :: emin,emax
integer(8)      :: operm, nperm
call mat_init(S,n,n)
call mat_init(S_sqrt,n,n)
call mat_init(S_minus_sqrt,n,n)
call mat_set_from_full(A,1E0_realk,S)
!if (.not.ASSOCIATED(S%raux)) ALLOCATE(S%raux(2))
!S%raux(1)=sqrt(mat_sqnorm2(S))     !emax
!S%raux(2)=1e-6                     !emin
call lowdin_schulz(S,S_sqrt,S_minus_sqrt,LUPRI)
call mat_to_full(S_minus_sqrt,1E0_realk,A)
call mat_free(S)
call mat_free(S_sqrt)
call mat_free(S_minus_sqrt)
END SUBROUTINE matrix_minushalf_df

!> \brief  Lowdin-Schulz decomposition
!> \author S. Reine
!> \date   2008
!> \param  S         The square matrix to decompose
!> \param  S_sqrt    The decomposed matrix
!> \param  LUPRI     Print unit number
SUBROUTINE lowdin_schulz(S, S_sqrt, S_minus_sqrt,lupri)
implicit none
type(Matrix)              :: S, S_sqrt, S_minus_sqrt, T1, T2
real(realk)               :: l2, alpha, beta, t1_converged, testx, emax, emin
real(realk)               :: tstart, tend
integer                   :: iter, lupri,nnz
logical                   :: converged = .false.
integer,parameter         :: MAX_ITER = 50

CALL LSTIMER('START ',TSTART,TEND,LUPRI)

call mat_init(T1,S%nrow,S%ncol)
call mat_init(T2,S%nrow,S%ncol)

S_sqrt=S
call mat_identity(S_minus_sqrt)
T1=S

t1_converged=sqrt(1E0_realk*S%ncol)
do iter=1, MAX_ITER

    if(iter.eq. 1) then
       emax = S%raux(1)
       emin = S%raux(2)
    else
       emax=1E0_realk
       emin=(1E0_realk/4E0_realk)*emin*(3E0_realk-emin)**2
    end if

    l2   = 2E0_realk/(emax+emin)
    emax = l2*emax
    emin = l2*emin

    alpha = -0.5E0_realk*l2*sqrt(l2)
    beta  =  1.5E0_realk*sqrt(l2)

    call mat_mul(S_sqrt,T1,'n','n',alpha,0E0_realk,T2)
    call mat_scal(beta,S_sqrt)
    call mat_daxpy(1E0_realk,T2, S_sqrt)

    call mat_mul(T1,S_minus_sqrt,'n','n',alpha,0.0E0_realk,T2)
    call mat_scal(beta,S_minus_sqrt)
    call mat_daxpy(1.0E0_realk,T2, S_minus_sqrt)

    call mat_mul(S_minus_sqrt,S_sqrt,'n','n',1.0E0_realk,0.0E0_realk, T1)

    testx = abs(sqrt(mat_sqnorm2(T1))-t1_converged)
    converged = testx .le.  1.0E-8_realk

    write(*,'(I2,A9,E14.6,A7,F14.6,A9,F14.6,A9,E14.6)')&
    &iter,    " Norm= ", testx, " l2= ", l2," emax= ", emax," emin= ", emin
    
    if (converged) exit
enddo

call mat_free(T1)
call mat_free(T2)

if (converged) then
    write(lupri,'(//1X,A,1X,I2,1X,A/,1X,A,E12.4/)') &
    &"Iterative Lowdin decomposition converged in", iter, "itarations",&
    &"Residual norm         :", testx
else
    write(lupri,'(//1X,A,1X,I2,1X,A/,1X,A,E12.4/)') &
    &"Iterative Lowdin decomposition did NOT CONVERGE in",MAX_ITER,"iterations",&
    &"Residual norm    :", testx
    call lsquit("Iterative Lowdin decomposition did NOT CONVERGED",lupri)
end if

CALL LSTIMER('MINUS HALF',TSTART,TEND,LUPRI)

END SUBROUTINE lowdin_schulz

END MODULE linsolvdf
