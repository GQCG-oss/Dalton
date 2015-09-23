!> @file
!> Contains linear-equation solver routines used for overlap density-fitting
MODULE linsolvdf
use precision
use TYPEDEF  
use matrix_operations
use ls_Integral_Interface
use LSparameters
use lstiming
use io
use io_type
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
!set threshold
SETTING%SCHEME%intTHRESHOLD = SETTING%SCHEME%THRESHOLD*SETTING%SCHEME%ONEEL_THR
call Test_if_64bit_integer_required(naux,naux)
call mem_alloc(alphabeta,naux,1,naux,1,1)
! Solve linear equations using standard library routines based on Cholesky-factorization
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

call mem_dealloc(alphabeta)

END SUBROUTINE linsolv_df

END MODULE linsolvdf
