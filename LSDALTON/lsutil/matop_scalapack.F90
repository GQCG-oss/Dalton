!> @file
!> Contains scalapack matrix module.
!> 
!> All functions are derived from their dense counterparts (mat_dense.f90) 
!> and adapted to the scalapack matrix format. The module requires scalapack.

module matrix_operations_scalapack

  use memory_handling
  use matrix_module
  use LSmatrix_type
  use lsmatrix_operations_dense
  use precision
#ifdef VAR_MPI
  use infpar_module
  use lsmpi_type
#endif

!  use module_scalapack_aux, if mpi is on 32bits, then blacs is most probably
!  too
  integer(kind=ls_mpik) :: scalapack_nodtot
  integer(kind=ls_mpik) :: scalapack_mynum
  integer(kind=ls_mpik) :: scalapack_comm
  logical :: scalapack_mpi_set
  logical :: scalapack_member
  Type Grid
     integer(kind=ls_mpik) :: ictxt !Blacs context
     integer(kind=ls_mpik) :: Masterictxt !Blacs system context
     integer(kind=ls_mpik) :: nprow !nprocess row = Number of process rows in the current process grid
     integer(kind=ls_mpik) :: npcol !nprocess col = Number of process columns in the current process grid
     integer(kind=ls_mpik) :: myrow !Row coordinate of the calling process in the process grid.
     integer(kind=ls_mpik) :: mycol !Column coordinate of the calling process in the process grid.
     integer(kind=ls_mpik) :: nprocs !number of MPI processors
     integer(kind=ls_mpik) :: mynum  !rank or process ID
     integer(kind=ls_mpik) :: comm  !grid MPI communicator
  end type Grid
  save
  integer,parameter  :: Job_gridinit = 1
  integer,parameter  :: Job_init = 2
  integer,parameter  :: Job_free = 3
  integer,parameter  :: Job_from_full = 4
  integer,parameter  :: Job_to_full = 5
  integer,parameter  :: Job_assign = 6
  integer,parameter  :: Job_trans = 7
  integer,parameter  :: Job_trace = 8
  integer,parameter  :: Job_set = 9
  integer,parameter  :: Job_dotproduct = 10
  integer,parameter  :: Job_scal = 11
  integer,parameter  :: Job_daxpy = 12
  integer,parameter  :: Job_mul = 13
  integer,parameter  :: Job_max = 14
  integer,parameter  :: Job_maxdiag = 15
  integer,parameter  :: Job_outdia_sqnrm2 = 16 
  integer,parameter  :: Job_extract_diag = 17
  integer,parameter  :: Job_ao_precond = 18
  integer,parameter  :: Job_add_identity = 19
  integer,parameter  :: Job_scal_diag = 20
  integer,parameter  :: Job_zerohalf = 21
  integer,parameter  :: Job_get_column = 22
  integer,parameter  :: Job_set_column = 23
  integer,parameter  :: Job_write_unit = 24
  integer,parameter  :: Job_read_unit = 25
  integer,parameter  :: Job_gpopen = 26
  integer,parameter  :: Job_gpclose = 27
  integer,parameter  :: Job_gprewind = 28
  integer,parameter  :: Job_eig = 29
  integer,parameter  :: Job_gridexit = 30
  integer,parameter  :: Job_absmax = 31
  integer,parameter  :: Job_diag_f = 32
  integer,parameter  :: Job_rand = 33
  integer,parameter  :: Job_retrieve_block = 34
  integer,parameter  :: Job_density_from_orbs = 35
  integer,parameter  :: Job_dsyev = 36
  integer,parameter  :: Job_dpotrf = 37
  integer,parameter  :: Job_dpotrs = 38
  integer,parameter  :: Job_test_to_full = 39
  integer,parameter  :: Job_add_block = 40
  integer,parameter  :: Job_print_global = 41
  integer,parameter  :: Job_dsyevx = 42
  integer,parameter  :: Job_setlowertriangular_zero=43
  integer,parameter  :: Job_scal_dia=44
  integer,parameter  :: Job_dmul = 45
  integer,parameter  :: Job_hmul = 46
  integer,parameter  :: Job_dger = 47
  integer,parameter  :: Job_create_block = 48
  integer,parameter  :: Job_hdiv = 49
  integer,parameter  :: Job_min = 50
  integer,parameter  :: Job_scal_dia_vec = 51
  integer,parameter  :: Job_to_full3D = 52
  Type(Grid) SLGrid
  
  Type(Matrix) :: darray(500)
  integer  :: BLOCK_SIZE
  integer, parameter :: BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1
  integer, parameter :: CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6
  integer, parameter :: RSRC_ = 7, CSRC_ = 8, LLD_ = 9
!
  integer :: gridinfo_nBlocks,gridinfo_nBlockssqrt
  integer,pointer :: gridinfo_blocks(:,:) !7,nBlocks the 7:globalrowstart,globalcolstart,localrowstart,localcolstart,owner,ownerrow,ownercol
  integer :: gridinfo_mynblocks
  integer,pointer :: gridinfo_myblocks(:) 

  contains
    SUBROUTINE DARRAY_NULLIFY
      IMPLICIT NONE
      INTEGER :: I
#ifdef VAR_SCALAPACK  
      DO I=1, UBOUND(darray,1)
         NULLIFY(darray(I)%p, darray(I)%addr_on_grid)
      ENDDO
#endif
    END SUBROUTINE DARRAY_NULLIFY

    SUBROUTINE DARRAY_FREEALL
      IMPLICIT NONE
      INTEGER :: I
#ifdef VAR_SCALAPACK  
      DO I=1, UBOUND(darray,1)
         IF (ASSOCIATED(darray(I)%p)) &
              & DEALLOCATE(darray(I)%p, darray(I)%addr_on_grid)
      ENDDO
#endif      
    END SUBROUTINE DARRAY_FREEALL

    INTEGER FUNCTION ALLOC_IN_DARRAY(A)
      IMPLICIT NONE
      Type(Matrix), target :: A
      INTEGER I, NUMROC
#ifdef VAR_SCALAPACK  
      EXTERNAL NUMROC
      ALLOC_IN_DARRAY = 0
      
      !NUMROC computes the NUMber of Rows Or Columns of a distributed
      !matrix owned by the process indicated by SLGrid%myrow/SLGrid%mycol.
      !infpar%master possesses the first row or column of the distributed matrix
      A%localnrow=NUMROC(A%nrow,BLOCK_SIZE,SLGrid%myrow,infpar%master,SLGrid%nprow)
      A%localncol=NUMROC(A%ncol,BLOCK_SIZE,SLGrid%mycol,infpar%master,SLGrid%npcol)

      IF(SLGrid%mynum.eq.infpar%master) THEN
         ALLOCATE(A%addr_on_grid(0:SLGrid%nprow-1,&
              &0:SLGrid%npcol-1))
      ELSE
         ALLOCATE(A%addr_on_grid(&
              & SLGrid%myrow:SLGrid%myrow,&
              & SLGrid%mycol:SLGrid%mycol))
      ENDIF
      A%addr_on_grid = 0       
      IF (A%localnrow*A%localncol.gt. 0)THEN
         ALLOCATE(A%p(A%localnrow,A%localncol))
      else
         ALLOCATE(A%p(1,1))
      endif
      DO I=1, UBOUND(darray,1)
         IF (.NOT.ASSOCIATED(darray(I)%p)) THEN
            EXIT
         ENDIF
      ENDDO
      A%addr_on_grid(SLGrid%myrow,SLGrid%mycol)=I
      darray(I) = A
      ALLOC_IN_DARRAY = I
#else
      ALLOC_IN_DARRAY = 0
#endif

     END FUNCTION ALLOC_IN_DARRAY
     
     SUBROUTINE FREE_IN_DARRAY(A)
       IMPLICIT NONE
       Type(Matrix), target :: A
       INTEGER I, NUMROC
#ifdef VAR_SCALAPACK  
       I=A%addr_on_grid(SLGrid%myrow,SLGrid%mycol)
       DEALLOCATE(A%p)
       DEALLOCATE(A%addr_on_grid)
       NULLIFY(A%p,darray(I)%p,&
            &A%addr_on_grid,darray(I)%addr_on_grid)
#endif       
     END SUBROUTINE FREE_IN_DARRAY     

#ifdef VAR_SCALAPACK  
      SUBROUTINE PDLADIAG( N, A, IA, JA, DESCA, FUNC, FDATA )
!     .. Scalar Arguments ..                                                                   
      INTEGER            IA, JA, N
!     ..                                                                                       
!     .. Array Arguments ..                                                                    
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A( * )
      DOUBLE PRECISION   FDATA( * )
!     ..                                                                                       
!                                                                                              
!  Purpose                                                                                     
!  =======                                                                                     
!                                                                                              
!  PDLADIAG returns diagonal of an N-by-N distributed matrix sub( A )                          
!  denoting A( IA:IA+N-1, JA:JA+N-1 ). The result is left on every                             
!  process of the grid.                                                                        
!                                                                                              
!  Arguments                                                                                   
!  =========                                                                                   
!                                                                                              
!  N       (global input) INTEGER                                                              
!          The number of rows and columns to be operated on i.e the                            
!          order of the distributed submatrix sub( A ).  N >= 0.                               
!                                                                                              
!  A       (local input) DOUBLE PRECISION pointer into the local memory                        
!          to an array of dimension ( LLD_A, LOCc(JA+N-1) ). This array                        
!          contains the local pieces of the distributed matrix the trace                       
!          is to be computed.                                                                  
!                                                                                              
!  IA      (global input) INTEGER                                                              
!          The row index in the global array A indicating the first                            
!          row of sub( A ).                                                                    
!                                                                                              
!  JA      (global input) INTEGER                                                              
!          The column index in the global array A indicating the                               
!          first column of sub( A ).                                                           
!                                                                                              
!  DESCA   (global and local input) INTEGER array of dimension DLEN_.                          
!          The array descriptor for the distributed matrix A.                                  
!                                                                                              
!  DIAG    output dimension (N)                                                                
!                                                                                              
!  ====================================================================                        
!                                                                                              
!     .. Parameters ..                                                                         
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0E+0_realk )
!     ..                                                                                       
!     .. Local Scalars ..                                                                      
      INTEGER            ICURCOL, ICURROW, II, IOFFA, J, JB, JJ, JN,&
     &                   LDA, LL, MYCOL, MYROW, NPCOL, NPROW, ID
      DOUBLE PRECISION   TRACE
!     ..                                                                                       
!     .. External Subroutines ..                                                               
      EXTERNAL           BLACS_GRIDINFO, DGSUM2D, INFOG2L
!     ..                                                                                       
!     .. External Functions ..                                                                 
      INTEGER            ICEIL
      EXTERNAL           ICEIL
!     ..                                                                                       
!     .. Intrinsic Functions ..                                                                
      INTRINSIC          MIN, MOD
!     ..                                                                                       
!     .. Executable Statements ..                                                              
!                                                                                              
!     Get grid parameters                                                                      
!                                                                                              
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
!                                                                                              
      DIAG = ZERO

      IF( N.EQ. 0 ) THEN
         RETURN
      END IF
!                                                                                              
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, II, JJ,&
     &              ICURROW, ICURCOL )
!                                                                                              
      JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
      JB = JN-JA+1
      LDA = DESCA( LLD_ )
      IOFFA = II + ( JJ - 1 ) * LDA
!                                                                                              
!     Handle first diagonal block separately                                                   
!                                                                                        
      IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
         DO 10 LL = IOFFA, IOFFA + (JB-1)*(LDA+1), LDA+1
!            ID=INDXL2G(ICEIL(LL,LDA),DESCA( NB_ ),MYROW, DESCA(RSRC_),NPROW)                  
!            DIAG(ID) = A( LL )                                                                
             CALL FUNC(DESCA,A,LL,FDATA)
   10    CONTINUE
      END IF
      IF( MYROW.EQ.ICURROW )&
     &   IOFFA = IOFFA + JB
      IF( MYCOL.EQ.ICURCOL )&
     &   IOFFA = IOFFA + JB*LDA
      ICURROW = MOD( ICURROW+1, NPROW )
      ICURCOL = MOD( ICURCOL+1, NPCOL )
!                                                                                              
!     Loop over the remaining block of columns                                                 
!                                                                                              
      DO 30 J = JN+1, JA+N-1, DESCA( NB_ )
         JB = MIN( JA+N-J, DESCA( NB_ ) )
!                                                                                              
         IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
            DO 20 LL = IOFFA, IOFFA + (JB-1)*(LDA+1), LDA+1
!            ID=INDXL2G(ICEIL(LL,LDA),DESCA( NB_ ),MYROW, DESCA(RSRC_),NPROW)                  
!            DIAG(ID) = A( LL )                                                                
            CALL FUNC(DESCA,A,LL,FDATA)
   20       CONTINUE
         END IF
         IF( MYROW.EQ.ICURROW )&
     &      IOFFA = IOFFA + JB
         IF( MYCOL.EQ.ICURCOL )&
     &      IOFFA = IOFFA + JB*LDA
         ICURROW = MOD( ICURROW+1, NPROW )
         ICURCOL = MOD( ICURCOL+1, NPCOL )
   30 CONTINUE
!                                                                                              
!      CALL DGSUM2D( DESCA( CTXT_ ), 'All', 'I', N, 1, DIAG, N, -1,&                           
!     &              MYCOL )                                                                   
!                                                                                              
      END SUBROUTINE

! Plugins for PDLADIAG                                                                      

     SUBROUTINE PLUGIN_EXTRACT(DESCA,A,LL,FDATA)
     IMPLICIT NONE
     INTEGER  :: DESCA(DLEN_), LL
     REAL(8)  :: A(*), FDATA(DESCA(M_))
     INTEGER  :: ID, INDXL2G, ICEIL
     EXTERNAL :: INDXL2G, ICEIL
     ID=INDXL2G(ICEIL(LL,DESCA(LLD_)),DESCA( NB_ ),&
          &        SLGrid%mycol,DESCA(RSRC_),SLGrid%npcol)
     FDATA(ID) = A( LL )
     END SUBROUTINE

     SUBROUTINE PLUGIN_MAX(DESCA,A,LL,FDATA)
     IMPLICIT NONE
     INTEGER  :: DESCA(DLEN_), LL
     REAL(REALK)  :: A(*), FDATA
     FDATA = MAX(FDATA,A(LL))
     END SUBROUTINE

     SUBROUTINE PLUGIN_SQ2(DESCA,A,LL,FDATA)
     IMPLICIT NONE
     INTEGER  :: DESCA(DLEN_), LL
     REAL(REALK)  :: A(*), FDATA
     FDATA = FDATA + A(LL)*A(LL)
     END SUBROUTINE

     SUBROUTINE PLUGIN_ADD(DESCA,A,LL,FDATA)
     IMPLICIT NONE
     INTEGER  :: DESCA(DLEN_), LL
     REAL(REALK)  :: A(*), FDATA
     A(LL) = A(LL) + FDATA
     END SUBROUTINE

     SUBROUTINE PLUGIN_SCAL(DESCA,A,LL,FDATA)
     IMPLICIT NONE
     INTEGER  :: DESCA(DLEN_), LL
     REAL(REALK)  :: A(*), FDATA
     A(LL) = A(LL)*FDATA
     END SUBROUTINE
#endif

     SUBROUTINE PDM_SYNC(JOB,A,B,C)
       IMPLICIT NONE
       INTEGER :: TMP(10), i,j, JOB, IERR
       TYPE(Matrix), optional :: A, B, C
#ifdef VAR_SCALAPACK  
      
       TMP(1:10) = 0      
       IF(scalapack_mynum.eq.infpar%master) then !code for master                              
          !Wake up ALL slaves 
          call ls_mpibcast(PDMSLAVE,infpar%master,MPI_COMM_LSDALTON)
          IF (PRESENT(A)) THEN
             TMP(1)=A%nrow ; TMP(2)=A%ncol
          ENDIF
          IF (PRESENT(B)) THEN
             TMP(4)=B%nrow ; TMP(5)=B%ncol
          ENDIF
          IF (PRESENT(C)) THEN
             TMP(7)=C%nrow ; TMP(8)=C%ncol
          ENDIF
          TMP(10)=JOB;
          ! Sends info about Job, and Matrices involved in Job
          ! IGESD2D Takes a general rectangular matrix [here a (10,1)] and sends it to 
          ! the destination process
          ! i is the The process row coordinate of the process to send the message to
          ! j is the The process col coordinate of the process to send the message to
         do i=0,SLGrid%nprow-1
            do j=0,SLGrid%npcol-1
               IF (PRESENT(A)) TMP(3) = A%addr_on_grid(i,j)
               IF (PRESENT(B)) TMP(6) = B%addr_on_grid(i,j)
               IF (PRESENT(C)) TMP(9) = C%addr_on_grid(i,j)
               if (.not.((i.eq. 0).and.(j.eq. 0))) then
                  CALL IGESD2D(SLGrid%ictxt,10,1,TMP,10,i,j)
               endif
            enddo
         enddo
      else !code for slaves
         !recieves TMP array
         CALL IGERV2D(SLGrid%ictxt,10,1,TMP,10,0,0)
         JOB = TMP(10) !slaves needs to know what to do
         IF (TMP(3).gt. 0) THEN
            A = darray(TMP(3))
         ELSE
            A%nrow = TMP(1); A%ncol = TMP(2)
         ENDIF
         IF (TMP(6).gt. 0) THEN
            B = darray(TMP(6))
         ELSE
            B%nrow = TMP(4); B%ncol = TMP(5)
         ENDIF
         IF (TMP(9).gt. 0) THEN
            C = darray(TMP(9))
         ELSE
            C%nrow = TMP(7); C%ncol = TMP(8)
         ENDIF
      endif
#endif
    END SUBROUTINE PDM_SYNC

    !this is for when the slaves is already awake
     SUBROUTINE PDM_MATRIXSYNC(A,B,C)
       IMPLICIT NONE
       INTEGER :: TMP(10), i,j, IERR
       TYPE(Matrix), optional :: A, B, C
#ifdef VAR_SCALAPACK  
      
       TMP(1:10) = 0      
       IF(scalapack_mynum.eq.infpar%master) then !code for master                              
          IF (PRESENT(A)) THEN
             TMP(1)=A%nrow ; TMP(2)=A%ncol
          ENDIF
          IF (PRESENT(B)) THEN
             TMP(4)=B%nrow ; TMP(5)=B%ncol
          ENDIF
          IF (PRESENT(C)) THEN
             TMP(7)=C%nrow ; TMP(8)=C%ncol
          ENDIF
          ! Sends info about Job, and Matrices involved in Job
          ! IGESD2D Takes a general rectangular matrix [here a (10,1)] and sends it to 
          ! the destination process
          ! i is the The process row coordinate of the process to send the message to
          ! j is the The process col coordinate of the process to send the message to
         do i=0,SLGrid%nprow-1
            do j=0,SLGrid%npcol-1
               IF (PRESENT(A)) TMP(3) = A%addr_on_grid(i,j)
               IF (PRESENT(B)) TMP(6) = B%addr_on_grid(i,j)
               IF (PRESENT(C)) TMP(9) = C%addr_on_grid(i,j)
               if (.not.((i.eq. 0).and.(j.eq. 0))) then
                  CALL IGESD2D(SLGrid%ictxt,10,1,TMP,10,i,j)
               endif
            enddo
         enddo
      else !code for slaves
         !recieves TMP array
         CALL IGERV2D(SLGrid%ictxt,10,1,TMP,10,0,0)
         IF(present(A))THEN
            IF (TMP(3).gt. 0) THEN
               A = darray(TMP(3))
            ELSE
               A%nrow = TMP(1); A%ncol = TMP(2)
            ENDIF
         ENDIF
         IF(present(B))THEN
            IF (TMP(6).gt. 0) THEN
               B = darray(TMP(6))
            ELSE
               B%nrow = TMP(4); B%ncol = TMP(5)
            ENDIF
         ENDIF
         IF(present(C))THEN
            IF (TMP(9).gt. 0) THEN
               C = darray(TMP(9))
            ELSE
               C%nrow = TMP(7); C%ncol = TMP(8)
            ENDIF
         ENDIF
      endif
#endif
    END SUBROUTINE PDM_MATRIXSYNC

    SUBROUTINE PDM_DSCINIT(DESC, A, BSnrow, BSncol)
      IMPLICIT NONE
      INTEGER           :: DESC(9), LLD, INFO, NUMROC
      INTEGER, optional :: BSnrow, BSncol
      INTEGER           :: bnrow, bncol 
      Type(Matrix)      :: A
#ifdef VAR_SCALAPACK  
      EXTERNAL NUMROC
      bnrow = BLOCK_SIZE
      bncol = BLOCK_SIZE
      IF (PRESENT(BSnrow)) bnrow = BSnrow
      IF (PRESENT(BSncol)) bncol = BSncol
      LLD=MAX(1,NUMROC(A%nrow, bnrow, SLGrid%myrow,0,SLGrid%nprow))
      !DESCINIT initializes the descriptor vector DESC
      !WARNING I think the two zeros should be infpar%master
      CALL DESCINIT(DESC,A%nrow,A%ncol,bnrow,bncol,0,0,SLGrid%ictxt,LLD,INFO)
      IF(INFO.NE. 0)THEN
         print*,'ERROR IN DESCINIT   INFO = ',INFO
      ENDIF
#endif
    END SUBROUTINE PDM_DSCINIT

    SUBROUTINE PDM_SCAL(ALPHA,A)
      IMPLICIT NONE
      Type(Matrix) :: A
      REAL(REALK)      :: ALPHA
#ifdef VAR_SCALAPACK  
      
      CALL PDM_SYNC(Job_scal,A)
      CALL DGEBS2D(SLGrid%ictxt,'A','I',1,1,ALPHA,1)
      if (A%localnrow*A%localncol.NE. 0)THEN
         CALL DSCAL(A%localnrow*A%localncol,ALPHA, A%p, 1)
      endif
#endif
    END SUBROUTINE PDM_SCAL

    !assigns B to be equal to A (B=A)
    SUBROUTINE PDM_ASSIGN(B,A)
      IMPLICIT NONE
      TYPE(Matrix), INTENT(INOUT) :: B
      TYPE(Matrix), INTENT(IN)    :: A
      INTEGER      :: DESC_AB(9)
#ifdef VAR_SCALAPACK  
      CALL PDM_SYNC(Job_assign,A,B)
      CALL PDM_DSCINIT(DESC_AB,A)
      ! copy A to B
      ! PDLACPY copies all or part of a distributed matrix A to another
      ! distributed matrix B.
      CALL PDLACPY('All',A%nrow,A%ncol,A%p,1,1,DESC_AB,&
           &B%p,1,1,DESC_AB)
#endif      
    END SUBROUTINE PDM_ASSIGN
    
    SUBROUTINE PDM_TRANS(A,B)
      IMPLICIT NONE
      Type(Matrix) :: A, B
      INTEGER      :: DESC_A(9), DESC_B(9)
#ifdef VAR_SCALAPACK  
      CALL PDM_SYNC(Job_trans,A,B)
      CALL PDM_DSCINIT(DESC_A,A)
      CALL PDM_DSCINIT(DESC_B,B)
      CALL PDTRAN(A%ncol,A%nrow,1.0E0_realk,A%p,1,1,DESC_A,&
           & 0.0E0_realk,B%p,1,1,DESC_B)
#endif
    END SUBROUTINE PDM_TRANS
      
    SUBROUTINE PDM_SET(A, ALPHA, BETA)
      IMPLICIT NONE
      Type(Matrix) :: A
      INTEGER      :: DESC_A(9)
      REAL(REALK)      :: ALPHA, BETA, AB(2)
#ifdef VAR_SCALAPACK  
      CALL PDM_SYNC(Job_set,A)
      CALL PDM_DSCINIT(DESC_A,A)
      AB(1)=ALPHA; AB(2)=BETA
      CALL DGEBS2D(SLGrid%ictxt,'A','I',2,1,AB,2)
      CALL PDLASET('A',A%nrow,A%ncol,ALPHA,BETA,A%p,&
           &       1,1,DESC_A)
#endif
    END SUBROUTINE PDM_SET

    SUBROUTINE PDM_IDENTITY(A)
      IMPLICIT NONE
      Type(Matrix) :: A
      CALL PDM_SET(A,0.0E0_realk,1.0E0_realk)
      
    END SUBROUTINE PDM_IDENTITY

!> \brief See mat_init in mat-operations.f90
  subroutine mat_scalapack_init(A,nrow,ncol)
     implicit none
     TYPE(Matrix) :: A
     integer, intent(in) :: nrow, ncol
     integer :: i,j,k
     integer(kind=long) :: nsizeFULL,nsizeLOCAL
#ifdef VAR_SCALAPACK        
     A%nrow = nrow
     A%ncol = ncol
     K = ALLOC_IN_DARRAY(A)
     !call ls_izero(A%addr_on_grid(0:SLGrid%nprow-1,0:SLGrid%npcol-1),SLGrid%nprow*SLGrid%npcol)
     call ls_izero(A%addr_on_grid,SLGrid%nprow*SLGrid%npcol)
     A%addr_on_grid(0,0) = K
     CALL PDM_SYNC(Job_init,A)
#if 0
     !Receives a message from the process into the general rectangular matrix.
     do i=0,SLGrid%nprow-1
        do j=0,SLGrid%npcol-1
           if (.not.((i.eq. 0).and.(j.eq. 0)))then
              CALL IGERV2D(SLGrid%ictxt,1,1,A%addr_on_grid(i,j),1,i,j)
           endif
        enddo
     enddo
#else
!    CALL lsmpi_reduction(A%addr_on_grid(0:SLGrid%nprow-1,0:SLGrid%npcol-1),SLGrid%nprow,&
!         & SLGrid%npcol,infpar%master,scalapack_comm)
     CALL lsmpi_reduction(A%addr_on_grid,SLGrid%nprow,SLGrid%npcol,infpar%master,scalapack_comm)
#endif
     nsizeFULL = A%nrow*A%ncol*mem_realsize
     nsizeLOCAL= A%localnrow*A%localncol*mem_realsize
     call mem_allocated_mem_type_matrix(nsizeLOCAL,nsizeFULL)
#endif     
   end subroutine mat_scalapack_init

  subroutine mat_scalapack_random(A)
     implicit none
     TYPE(Matrix) :: A
     integer :: i,j
#ifdef VAR_SCALAPACK        
     do i=1,A%localnrow
        do j=1,A%localncol
           A%p(i,j)=scalapack_mynum+1
        enddo
     enddo
     CALL PDM_SYNC(Job_rand,A)
#endif
   end subroutine mat_scalapack_random

!> \brief See mat_free in mat-operations.f90
  subroutine mat_scalapack_free(A)
     implicit none
     TYPE(Matrix) :: A
     integer(kind=long) :: nsizeFULL,nsizeLOCAL
#ifdef VAR_SCALAPACK        
     nsizeFULL = A%nrow*A%ncol*mem_realsize
     nsizeLOCAL= A%localnrow*A%localncol*mem_realsize
     call mem_deallocated_mem_type_matrix(nsizeLOCAL,nsizeFULL)
     CALL PDM_SYNC(Job_free,A)
     CALL FREE_IN_DARRAY(A)
#endif
   end subroutine mat_scalapack_free

!> \brief See mat_set_from_full in mat-operations.f90
   subroutine mat_scalapack_set_from_full(afull,alpha,a)
     implicit none
     real(realk), INTENT(IN)    :: afull(*)
     real(realk), intent(in) :: alpha
     TYPE(Matrix)               :: a 
     INTEGER      :: DESC_AF(9), DESC_A(9), I,nsize,n1,n2,shift
     INTEGER :: J,IP,IQ,P,Q,ip2,iq2,NBJ,NBI,iproc,preproc,nprow,npcol
#ifdef VAR_SCALAPACK        
     integer :: localnrow(0:infpar%nodtot-1)
     integer :: localncol(0:infpar%nodtot-1)
     TYPE(matrix),pointer :: localA(:)

     CALL PDM_SYNC(Job_from_full,A)
     ! Set up descriptors                                                      
     CALL PDM_DSCINIT(DESC_AF,A,A%nrow,A%ncol)
     CALL PDM_DSCINIT(DESC_A,A)
     !    DESC_AF(CTXT_) = SLGrid%Masterictxt
     ! Master Part: Distributes a full matrix AFull (on master) to SCALAPACK distributed matrix A
     IF(.NOT.infpar%ScalapackWorkAround)THEN
        CALL PDGEMR2D(A%nrow,A%ncol,AFull,1,1,DESC_AF,&
             & A%p,1,1,DESC_A,SLGrid%ictxt)
     ELSE
        shift=scalapack_nodtot/SLGrid%nprow
        J=0
        localnrow = 0
        localncol = 0
        jloop : do
           preProc = 0
           do npcol=1,SLGrid%npcol
              do NBJ=1,BLOCK_SIZE
                 J=J+1
                 IF(J.GT.A%ncol)EXIT jloop
                 Iproc = 0
                 Ploop : DO
                    IF(preProc+Iproc*shift.GT.scalapack_nodtot-1)EXIT Ploop
                    localncol(preProc+Iproc*shift) = localncol(preProc+Iproc*shift)+1 
                    Iproc = Iproc+1
                 ENDDO Ploop
                 I=0
                 iloop3 : do
                    Iproc=preProc
                    do nprow=1,SLGrid%nprow
                       localnrow(Iproc) = 0
                       do NBI=1,BLOCK_SIZE
                          I=I+1
                          IF(I.GT.A%nrow)EXIT iloop3               
                       enddo
                       Iproc=Iproc+shift
                    enddo
                 enddo iloop3
                 I=0
                 iloop : do
                    Iproc=preProc
                    do nprow=1,SLGrid%nprow
                       do NBI=1,BLOCK_SIZE
                          I=I+1
                          IF(I.GT.A%nrow)EXIT iloop               
                          localnrow(Iproc) = localnrow(Iproc)+1 
                       enddo
                       Iproc=Iproc+shift
                    enddo
                 enddo iloop
              enddo
              preProc = preProc+1
           enddo
        enddo jloop

        allocate(localA(0:scalapack_nodtot-1))
        do I=0,scalapack_nodtot-1
           localA(I)%nrow = localnrow(I)
           localA(I)%ncol = localncol(I)
           allocate(localA(I)%p(localnrow(I),localncol(I)))
        enddo

        shift=scalapack_nodtot/SLGrid%nprow
        J=0
        localnrow = 0
        localncol = 0
        jloop2 : do
           preProc = 0
           do npcol=1,SLGrid%npcol
              do NBJ=1,BLOCK_SIZE
                 J=J+1
                 IF(J.GT.A%ncol)EXIT jloop2
                 Iproc = 0
                 Ploop2 : DO
                    IF(preProc+Iproc*shift.GT.scalapack_nodtot-1)EXIT Ploop2
                    localncol(preProc+Iproc*shift) = localncol(preProc+Iproc*shift)+1 
                    Iproc = Iproc+1
                 ENDDO Ploop2
                 I=0
                 localnrow = 0
                 iloop2 : do
                    Iproc=preProc
                    do nprow=1,SLGrid%nprow
                       do NBI=1,BLOCK_SIZE
                          I=I+1
                          IF(I.GT.A%nrow)EXIT iloop2               
                          localnrow(Iproc) = localnrow(Iproc)+1 
                          localA(Iproc)%p(localnrow(Iproc),localncol(Iproc)) = AFULL(I+(J-1)*A%nrow) 
                          !                      write(6,'(A,I2,A,I2,A,I2,A,I2,A,I2,A,F16.8)')'Ap(',localnrow(Iproc),',',localncol(Iproc),',',iproc,') = AFULL(',&
                          !                           & I,',',J,')=',AFULL(I+(J-1)*A%nrow) 
                       enddo
                       Iproc=Iproc+shift
                    enddo
                 enddo iloop2
              enddo
              preProc = preProc+1
           enddo
        enddo jloop2

        do NBJ=1,A%localncol
           do NBI=1,A%localnrow
              A%p(NBI,NBJ) = localA(0)%p(NBI,NBJ)
           enddo
        enddo
        DO Iproc = 1,scalapack_nodtot-1
           IF(localA(Iproc)%nrow*localA(Iproc)%ncol.GT. 0)THEN
              call ls_mpisendrecv(localA(Iproc)%p,localA(Iproc)%nrow,localA(Iproc)%ncol,scalapack_comm,infpar%master,Iproc)
           ENDIF
        ENDDO
        do I=0,scalapack_nodtot-1
           deallocate(localA(I)%p)
        enddo
     ENDIF
     if (ALPHA.ne. 1.0E0_realk)call mat_scalapack_scal(alpha,a)
#endif
   end subroutine mat_scalapack_set_from_full

!> \brief See mat_to_full in mat-operations.f90
  subroutine mat_scalapack_to_full(a,alpha,afull)
     implicit none
    TYPE(Matrix),intent(in) :: a 
    real(realk), intent(in) :: alpha
    real(realk)  :: afull(*)
    INTEGER      :: DESC_AF(9), DESC_A(9), I,nsize,n1,n2,shift
    INTEGER :: J,IP,IQ,P,Q,ip2,iq2,NBJ,NBI,iproc,preproc,nprow,npcol
#ifdef VAR_SCALAPACK
    integer :: localnrow(0:infpar%nodtot-1)
    integer :: localncol(0:infpar%nodtot-1)
    TYPE(matrix),pointer :: localA(:)
    CALL PDM_SYNC(Job_to_full,A)
    CALL PDM_DSCINIT(DESC_AF,A,A%nrow,A%ncol)
    CALL PDM_DSCINIT(DESC_A, A)
    IF(.NOT.infpar%ScalapackWorkAround)THEN
       ! Master Part: Collects full matrix AFull (on master) from SCALAPACK distributed matrix A
       CALL PDGEMR2D(A%nrow,A%ncol,A%p,1,1,DESC_A,&
            &AFull,1,1,DESC_AF,SLGrid%ictxt)
    ELSE
       shift=scalapack_nodtot/SLGrid%nprow
       J=0
       localnrow = 0
       localncol = 0
       jloop : do
          preProc = 0
          do npcol=1,SLGrid%npcol
             do NBJ=1,BLOCK_SIZE
                J=J+1
                IF(J.GT.A%ncol)EXIT jloop
                Iproc = 0
                Ploop : DO
                   IF(preProc+Iproc*shift.GT.scalapack_nodtot-1)EXIT Ploop
                   localncol(preProc+Iproc*shift) = localncol(preProc+Iproc*shift)+1 
                   Iproc = Iproc+1
                ENDDO Ploop
                I=0
                iloop3 : do
                   Iproc=preProc
                   do nprow=1,SLGrid%nprow
                      localnrow(Iproc) = 0
                      do NBI=1,BLOCK_SIZE
                         I=I+1
                         IF(I.GT.A%nrow)EXIT iloop3               
                      enddo
                      Iproc=Iproc+shift
                   enddo
                enddo iloop3
                I=0
                iloop : do
                   Iproc=preProc
                   do nprow=1,SLGrid%nprow
                      do NBI=1,BLOCK_SIZE
                         I=I+1
                         IF(I.GT.A%nrow)EXIT iloop               
                         localnrow(Iproc) = localnrow(Iproc)+1 
                      enddo
                      Iproc=Iproc+shift
                   enddo
                enddo iloop
             enddo
             preProc = preProc+1
          enddo
       enddo jloop
       
       allocate(localA(0:scalapack_nodtot-1))
       do I=0,scalapack_nodtot-1
          localA(I)%nrow = localnrow(I)
          localA(I)%ncol = localncol(I)
          allocate(localA(I)%p(localnrow(I),localncol(I)))
       enddo
       DO Iproc = 1,scalapack_nodtot-1
          IF(localA(Iproc)%nrow*localA(Iproc)%ncol.GT. 0)THEN
             call ls_mpisendrecv(localA(Iproc)%p,localA(Iproc)%nrow,localA(Iproc)%ncol,scalapack_comm,Iproc,infpar%master)
          ENDIF
       ENDDO       
       do NBJ=1,A%localncol
          do NBI=1,A%localnrow
             localA(0)%p(NBI,NBJ) = A%p(NBI,NBJ)
          enddo
       enddo

       shift=scalapack_nodtot/SLGrid%nprow
       J=0
       localnrow = 0
       localncol = 0
       jloop2 : do
          preProc = 0
          do npcol=1,SLGrid%npcol
             do NBJ=1,BLOCK_SIZE
                J=J+1
                IF(J.GT.A%ncol)EXIT jloop2
                Iproc = 0
                Ploop2 : DO
                   IF(preProc+Iproc*shift.GT.scalapack_nodtot-1)EXIT Ploop2
                   localncol(preProc+Iproc*shift) = localncol(preProc+Iproc*shift)+1 
                   Iproc = Iproc+1
                ENDDO Ploop2
                I=0
                localnrow = 0
                iloop2 : do
                   Iproc=preProc
                   do nprow=1,SLGrid%nprow
                      do NBI=1,BLOCK_SIZE
                         I=I+1
                         IF(I.GT.A%nrow)EXIT iloop2               
                         localnrow(Iproc) = localnrow(Iproc)+1 
                         AFULL(I+(J-1)*A%nrow) = localA(Iproc)%p(localnrow(Iproc),localncol(Iproc))
!                         write(6,'(A,I2,A,I2,A,I2,A,I2,A,I2,A,F16.8)')'Ap(',localnrow(Iproc),',',localncol(Iproc),',',iproc,') = AFULL(',&
!                              & I,',',J,')=',AFULL(I+(J-1)*A%nrow) 
                      enddo
                      Iproc=Iproc+shift
                   enddo
                enddo iloop2
             enddo
             preProc = preProc+1
          enddo
       enddo jloop2

       do I=0,scalapack_nodtot-1
          deallocate(localA(I)%p)
       enddo
    ENDIF
    nsize = a%nrow*a%ncol
    if (alpha.ne. 1.0E0_realk) call dscal(nsize,alpha,afull,1)
#endif

  end subroutine mat_scalapack_to_full

!> \brief See mat_to_full in mat-operations.f90
  subroutine mat_scalapack_to_full3D(a,alpha,afull,n1,n2,n3,i3)
     implicit none
    integer, INTENT(IN)     :: n1,n2,n3,i3
    TYPE(Matrix),intent(in) :: a 
    real(realk), intent(in) :: alpha
    real(realk)  :: afull(n1,n2,n3)
    INTEGER      :: DESC_AF(9), DESC_A(9), I,nsize
    INTEGER      :: J,n,m,mp1,OFFSET
    real(realk),pointer :: Afulltmp(:)
#ifdef VAR_SCALAPACK
    nsize = a%nrow*a%ncol
    IF(n1*n2.NE.nsize)call lsquit('dim mismatch mat_scalapack_to_full3D',-1)
!    integer :: localnrow(0:infpar%nodtot-1)
!    integer :: localncol(0:infpar%nodtot-1)
!    TYPE(matrix),pointer :: localA(:)
    CALL PDM_SYNC(Job_to_full3D,A)
    CALL PDM_DSCINIT(DESC_AF,A,A%nrow,A%ncol)
    CALL PDM_DSCINIT(DESC_A, A)
    ! Master Part: Collects full matrix AFull (on master) from SCALAPACK distributed matrix A
    call mem_alloc(Afulltmp,n1*n2)
    CALL PDGEMR2D(A%nrow,A%ncol,A%p,1,1,DESC_A,&
         &AFulltmp,1,1,DESC_AF,SLGrid%ictxt)    

    N = a%nrow    !change diff
    M = MOD(N,7)  !change diff
    IF (M.NE.0) THEN
       do j = 1,a%ncol
          offset = (j-1)*N
          DO I = 1,M
             afull(i,j,i3) = alpha*afulltmp(i+offset)              
          ENDDO
       enddo
       MP1 = M + 1
       IF (N.GE.7)THEN
          do j = 1,a%ncol
             offset = (j-1)*N
             DO I = MP1,N,7
                afull(i,j,i3) = alpha*afulltmp(i+offset)
                afull(i+1,j,i3) = alpha*afulltmp(i+1+offset)
                afull(i+2,j,i3) = alpha*afulltmp(i+2+offset)
                afull(i+3,j,i3) = alpha*afulltmp(i+3+offset)
                afull(i+4,j,i3) = alpha*afulltmp(i+4+offset)
                afull(i+5,j,i3) = alpha*afulltmp(i+5+offset)
                afull(i+6,j,i3) = alpha*afulltmp(i+6+offset)
             END DO
          enddo
       ENDIF
    ELSE
       do j = 1,a%ncol
          offset = (j-1)*N
          DO I = 1,N,7
             afull(i,j,i3) = alpha*afulltmp(i+offset)
             afull(i+1,j,i3) = alpha*afulltmp(i+1+offset)
             afull(i+2,j,i3) = alpha*afulltmp(i+2+offset)
             afull(i+3,j,i3) = alpha*afulltmp(i+3+offset)
             afull(i+4,j,i3) = alpha*afulltmp(i+4+offset)
             afull(i+5,j,i3) = alpha*afulltmp(i+5+offset)
             afull(i+6,j,i3) = alpha*afulltmp(i+6+offset)
          END DO
       ENDDO
    ENDIF
    call mem_dealloc(Afulltmp)
#endif

  end subroutine mat_scalapack_to_full3D

  subroutine GRIDINFO_SETUP(nbast)
    implicit none
    integer,intent(in) :: nbast
    !
#ifdef VAR_SCALAPACK
    integer :: nRowBlocks,mynum,iBlock
    integer :: local_j,local_i,i,j,prowcoord,pcolcoord

    mynum = scalapack_mynum

    IF(MOD(nbast,BLOCK_SIZE).EQ.0)THEN
       nRowBlocks = (nbast/BLOCK_SIZE)
    ELSE
       nRowBlocks = (nbast/BLOCK_SIZE)+1
    ENDIF
    gridinfo_nBlocks = nRowBlocks*nRowBlocks
    gridinfo_nBlockssqrt = nRowBlocks
!    print*,'nBlocks:',gridinfo_nBlocks,'mynum',scalapack_mynum
    allocate(gridinfo_blocks(7,gridinfo_nBlocks))
    iBlock = 0
!    print*,'nbast',nbast
    do j=1,nbast,BLOCK_SIZE
       pcolcoord = INDXG2P(j,BLOCK_SIZE,0,0,slgrid%npcol)
       do i=1,nbast,BLOCK_SIZE
          prowcoord = INDXG2P(i,BLOCK_SIZE,0,0,slgrid%nprow)
!          print*,'global index',i,j
          local_i = INDXG2L(i,BLOCK_SIZE,0,0,slgrid%nprow)
          local_j = INDXG2L(j,BLOCK_SIZE,0,0,slgrid%nprow)
          iBlock = iBlock+1
!          print*,'iblock',iblock
          gridinfo_blocks(1,iBlock) = i
          gridinfo_blocks(2,iBlock) = j
          gridinfo_blocks(3,iBlock) = local_i
          gridinfo_blocks(4,iBlock) = local_j
          gridinfo_blocks(5,iBlock) = (pcolcoord+1+prowcoord*slgrid%npcol)-1
          gridinfo_blocks(6,iBlock) = prowcoord
          gridinfo_blocks(7,iBlock) = pcolcoord
!          print*,'local index',local_i,local_j,' on grid ',prowcoord,pcolcoord,'mynum=',(pcolcoord+1+prowcoord*slgrid%npcol)-1
       enddo
    enddo
    if(iBlock.NE.gridinfo_nBlocks)call lsquit('error in determining gridinfo in GRIDINFO_SETUP',-1)
    gridinfo_mynblocks = 0
    do iBlock = 1,gridinfo_nBlocks
       if(gridinfo_blocks(5,iBlock).EQ.mynum) gridinfo_mynblocks = gridinfo_mynblocks+1
    enddo
    allocate(gridinfo_myblocks(gridinfo_mynBlocks))
    gridinfo_mynblocks = 0
    do iBlock = 1,gridinfo_nBlocks
       if(gridinfo_blocks(5,iBlock).EQ.mynum)then
          gridinfo_mynblocks = gridinfo_mynblocks+1
          gridinfo_myblocks(gridinfo_mynblocks) = iBlock
       endif
    enddo
!    print*,'mynum',mynum,'gridinfo_mynblocks',gridinfo_mynblocks,'gridinfo_myblocks',gridinfo_myblocks
#endif
  end subroutine GRIDINFO_SETUP

  subroutine GRIDINFO_FREE()
    implicit none
    deallocate(gridinfo_blocks)
    deallocate(gridinfo_myblocks)
  end subroutine GRIDINFO_FREE

  subroutine mat_scalapack_add_block(A,Afullblock,fullrow,fullcol,insertrow,insertcol)
    implicit none
    integer, intent(in) :: fullrow,fullcol,insertrow,insertcol
    real(Realk),target :: Afullblock(fullrow,fullcol)
    type(Matrix) :: A  !, intent(inout)
    INTEGER      :: DESC_A(9), DESC_AF(9)
    INTEGER      :: LLD, INFO, TMP(4)
    integer, external :: NUMROC
#ifdef VAR_SCALAPACK

    CALL PDM_SYNC(Job_add_block,A)
    CALL PDM_DSCINIT(DESC_A, A)
    TMP(1) = fullrow
    TMP(2) = fullcol
    TMP(3) = insertrow
    TMP(4) = insertcol
    call ls_mpibcast(TMP,4,infpar%master,scalapack_comm)

    LLD =  MAX(1,NUMROC(fullrow, fullrow, SLGrid%myrow,0,SLGrid%nprow))
    CALL DESCINIT(DESC_AF,fullrow,fullcol,fullrow,fullcol,0,0,SLGrid%ictxt,LLD,INFO)

    CALL PDGEADD('N',fullrow,fullcol,1.0_realk,Afullblock,1,1,DESC_AF, &
               & 1.0_realk,A%p,insertrow,insertcol,DESC_A)
#endif
  end subroutine mat_scalapack_add_block




!> \brief See mat_retrieve_block in mat-operations.f90
  subroutine mat_scalapack_retrieve_block(A,Afullblock,fullrow,fullcol,insertrow,insertcol)
    implicit none
    integer, intent(in) :: fullrow,fullcol,insertrow,insertcol
    real(Realk),target :: Afullblock(fullrow,fullcol)
    type(Matrix) :: A  !, intent(inout)
    INTEGER      :: DESC_A(9), DESC_AF(9)
    INTEGER      :: LLD, INFO, TMP(4),i,j
    integer, external :: NUMROC
    real(realk), pointer :: ABigfull(:,:)
#ifdef VAR_SCALAPACK

    IF(.NOT.infpar%ScalapackWorkAround)THEN
       CALL PDM_SYNC(Job_retrieve_block,A)
       CALL PDM_DSCINIT(DESC_A, A)
       TMP(1) = fullrow
       TMP(2) = fullcol
       TMP(3) = insertrow
       TMP(4) = insertcol
       call ls_mpibcast(TMP,4,infpar%master,scalapack_comm)
       
       LLD =  MAX(1,NUMROC(fullrow, fullrow, SLGrid%myrow,0,SLGrid%nprow))
       CALL DESCINIT(DESC_AF,fullrow,fullcol,fullrow,fullcol,0,0,SLGrid%ictxt,LLD,INFO)
       
       CALL PDGEMR2D(fullrow,fullcol,A%p,insertrow,insertcol,DESC_A, &
            &Afullblock,1,1,DESC_AF,SLGrid%ictxt)
    ELSE
       !Work around for PDGEMR2D issue 
       call mem_alloc(ABigfull,A%nrow,A%ncol)
       call mat_scalapack_to_full(A,1.0E0_realk,aBigfull)
       !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i,j) SHARED(ABigfull,&
       !$OMP fullrow,fullcol,insertrow,insertcol,Afullblock)
       do j = insertcol, insertcol+fullcol-1
          do i = insertrow, insertrow+fullrow-1
             Afullblock(i-insertrow+1,j-insertcol+1) = ABigfull(i,j)
          enddo
       enddo
       !$OMP END PARALLEL DO
       call mat_scalapack_set_from_full(aBigfull,1.0E0_realk,A)
       call mem_dealloc(ABigfull)
    ENDIF

#endif
  end subroutine mat_scalapack_retrieve_block

  subroutine mat_scalapack_create_block(A,Afull,fullrow,fullcol,insertrow,insertcol)
    implicit none
    integer, intent(in) :: fullrow,fullcol,insertrow,insertcol
    real(Realk),target :: Afull(fullrow,fullcol)
    type(Matrix) :: A  !, intent(inout)
    INTEGER      :: DESC_A(9), DESC_AF(9)
    INTEGER      :: LLD, INFO, TMP(4),i,j
    integer, external :: NUMROC
    real(realk), pointer :: ABigfull(:,:)
#ifdef VAR_SCALAPACK

    IF(.NOT.infpar%ScalapackWorkAround)THEN
       CALL PDM_SYNC(Job_create_block,A)
       CALL PDM_DSCINIT(DESC_A, A)
       TMP(1) = fullrow
       TMP(2) = fullcol
       TMP(3) = insertrow
       TMP(4) = insertcol
       call ls_mpibcast(TMP,4,infpar%master,scalapack_comm)
       
       LLD =  MAX(1,NUMROC(fullrow, fullrow, SLGrid%myrow,0,SLGrid%nprow))
       CALL DESCINIT(DESC_AF,fullrow,fullcol,fullrow,fullcol,0,0,SLGrid%ictxt,LLD,INFO)
       
       CALL PDGEMR2D(fullrow,fullcol,Afull,1,1,DESC_AF, &
                 &A%p,insertrow,insertcol,DESC_A,SLGrid%ictxt)
    ELSE
       !Work around for PDGEMR2D issue 
       call mem_alloc(ABigfull,A%nrow,A%ncol)
       call mat_scalapack_to_full(A,1.0E0_realk,aBigfull)
       !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i,j) SHARED(ABigfull,&
       !$OMP fullrow,fullcol,insertrow,insertcol,Afull)
       do j = insertcol, insertcol+fullcol-1
          do i = insertrow, insertrow+fullrow-1
             ABigfull(i,j) = Afull(i-insertrow+1,j-insertcol+1)
          enddo
       enddo
       !$OMP END PARALLEL DO
       call mat_scalapack_set_from_full(aBigfull,1.0E0_realk,A)
       call mem_dealloc(ABigfull)
    ENDIF
#endif
  end subroutine  

  subroutine scalapack_transform_aux(A,Afullblock,fulldim1,fulldim2,fullrow,fullcol,insertrow,insertcol,&
       & CommunicationNodes,localnrow2,localncol,Abuffer,nbuffer,bufferoffset,debug,bufferproc,nodes,funccall)
    implicit none
    integer, intent(in) :: fullrow,fullcol,insertrow,insertcol,nbuffer,fulldim1,fulldim2,bufferproc,nodes
    real(Realk)  :: Afullblock(fulldim1,fulldim2)
    type(Matrix) :: A  !, intent(inout)
    type(lsMatrix) :: Abuffer(0:nbuffer)  !, intent(inout)
    logical :: CommunicationNodes(0:nodes),debug    
    integer :: localnrow2(0:nodes),bufferoffset(0:nodes)
    integer :: localncol(0:nodes)
    EXTERNAL funccall !NAME OF SUBROUTINE TO CALL    
    !
#ifdef VAR_SCALAPACK
    integer :: localnrow(0:infpar%nodtot-1)
    INTEGER      :: DESC_AF(9), DESC_A(9), I,nsize,n1,n2,shift
    INTEGER :: J,IP,IQ,P,Q,ip2,iq2,NBJ,NBI,iproc,preproc,nprow,npcol
    integer :: ndim1,ndim2,TMP(4),RowSize,ColSize
    integer :: nRowBlocks,nColBlocks,iBlock,local_i,local_j,RowSource,ColSource
    integer :: pcolcoord,prowcoord,mynum,target_ndim,startblock,nBlocksRow,nBlocksCol,jBlock,nblocks,II,tmpI,nunique,srow2,srow1,scol2,scol1
    integer,pointer :: blocks(:),blocks2(:),blocks3(:)
    TYPE(matrix),pointer :: localA(:)
    real(realk),pointer :: Bfull(:,:)
    integer :: VI,VJ,L,K,iblocksRow,iBlocksCol,Vndim1,Vndim2,iBlockOffset,jBlockOffset
    logical :: AddBlock(0:infpar%nodtot-1)    
    logical :: AddRowBlock(0:infpar%nodtot-1)    
    logical :: AddRowBlock2(0:infpar%nodtot-1)    
    integer :: Global_I,Global_J,Buffer_I,Buffer_J,startRowSize,startColSize
    integer :: localnrow3(0:infpar%nodtot-1),localncol3(0:infpar%nodtot-1)
!    real(realk),pointer :: VISUAL(:,:)
    !initialization
    CommunicationNodes=.FALSE.
    localnrow3 = 0
    localncol3 = 0
    localnrow = 0
    localncol = 0
    shift=infpar%nodtot/SLGrid%nprow
    Vndim1 = 0
    Vndim2 = 0
    TMP(1) = insertrow
    TMP(2) = insertrow+fullrow-1
    TMP(3) = insertcol
    TMP(4) = insertcol+fullcol-1
    ndim1 = TMP(2)-TMP(1)+1
    ndim2 = TMP(4)-TMP(3)+1
!    if(debug)then
!       if(fulldim1.NE.ndim1.AND.fulldim2.NE.ndim2)call lsquit('error visual1',-1)
!       call mem_alloc(VISUAL,fulldim1,fulldim2)
!    endif

    IF(MOD(TMP(1),BLOCK_SIZE).EQ.0)THEN
       startRowSize = 1
    ELSE
       startRowsize = BLOCK_SIZE-MOD(TMP(1),BLOCK_SIZE)+1
    ENDIF
    IF(MOD(TMP(3),BLOCK_SIZE).EQ.0)THEN
       startColSize = 1
    ELSE
       startColsize = BLOCK_SIZE-MOD(TMP(3),BLOCK_SIZE)+1
    ENDIF
    
    !step 1 determine start block 
    iBlock=0
    do i=1,TMP(1),BLOCK_SIZE
       iBlock = iBlock+1
    enddo
    jBlock=0
    do i=1,TMP(3),BLOCK_SIZE
       jBlock = jBlock+1
    enddo
    !this is the Matrix block that have the element
    !corresponding to the element Afullblock(1,1) 
    !which correspond to Aglobalfull(insertrow,insertcol)
    startblock = iBlock+(jBlock-1)*gridinfo_nBlockssqrt
    iBlockOffset = iBlock-1
    jBlockOffset = jBlock-1
    !step 2 determine number of row blocks and number of collum blocks 
    IF(gridinfo_blocks(1,startblock).EQ.TMP(1))THEN
       !start of block
       IF(MOD(ndim1,BLOCK_SIZE).EQ.0)THEN
          nBlocksRow = ndim1/BLOCK_SIZE
       ELSE
          nBlocksRow = ndim1/BLOCK_SIZE+1
       ENDIF
    ELSE
       nBlocksRow=1
       I = gridinfo_blocks(1,startblock)
       DO 
          IF(I+BLOCK_SIZE.LE.TMP(2))THEN
             !include one more block
             I = I+BLOCK_SIZE
             nBlocksRow = nBlocksRow+1
          ELSE
             EXIT
          ENDIF
       ENDDO
    ENDIF
    IF(gridinfo_blocks(2,startblock).EQ.TMP(3))THEN
       !start of block
       IF(MOD(ndim2,BLOCK_SIZE).EQ.0)THEN
          nBlocksCol = ndim2/BLOCK_SIZE
       ELSE
          nBlocksCol = ndim2/BLOCK_SIZE+1
       ENDIF
    ELSE
       nBlocksCol=1
       I = gridinfo_blocks(2,startblock)
       DO 
          IF(I+BLOCK_SIZE.LE.TMP(4))THEN
             !include one more block
             I = I+BLOCK_SIZE
             nBlocksCol = nBlocksCol+1
          ELSE
             EXIT
          ENDIF
       ENDDO
    ENDIF
!    if(debug)then
!       print*,'nBlocksRow',nBlocksRow
!       print*,'nBlocksCol',nBlocksCol
!    endif
    IF(nBlocksRow.EQ.1.AND.nBlocksCol.EQ.1)THEN
       VI=1
       VJ=1
       iproc = gridinfo_blocks(5,startblock)
       Global_I = insertrow!gridinfo_blocks(1,startblock)
       Global_J = insertcol!gridinfo_blocks(2,startblock)
       RowSize = ndim1
       ColSize = ndim2
       Buffer_I = 1
       Buffer_J = 1
!       print*,'A GI,GJ,BI,BJ',Global_I,Global_J,Buffer_I,Buffer_J,'iproc',iproc
!       if(debug)call fillVisual(VISUAL,Global_I,Global_J,RowSize,ColSize,iproc,fulldim1,fulldim2)
       CALL funccall(A,Afullblock,fulldim1,fulldim2,RowSize,ColSize,&
            & Abuffer,nbuffer,iproc,Global_I,Global_J,Buffer_I,Buffer_J,VI,VJ,localnrow2,bufferoffset,bufferproc,debug)!,VISUAL)
       CommunicationNodes(iproc)=.TRUE.
       localnrow(iproc) = ndim1
       localncol(iproc) = ndim2
    ELSE
!       print*,'FIRST ROW BLOCK '
!       print*,'iBlocksRow',1,'iBlocksCol',1
       AddRowBlock=.TRUE.
       VI=1
       VJ=1
       iproc = gridinfo_blocks(5,startblock)
       Global_I = insertrow!gridinfo_blocks(1,startblock)
       Global_J = insertcol!gridinfo_blocks(2,startblock)
       IF(nBlocksRow.EQ.1)THEN
          RowSize = ndim1-VI+1
       ELSE
          RowSize = startRowSize
       ENDIF
       IF(nBlocksCol.EQ.1)THEN
          ColSize = ndim2-VJ+1
       ELSE
          ColSize = startColSize
       ENDIF

       Buffer_I = 1
       Buffer_J = 1
       !First row and collom Block
!       print*,'B GI,GJ,BI,BJ',Global_I,Global_J,Buffer_I,Buffer_J,iproc
!       if(debug)call fillVisual(VISUAL,Global_I,Global_J,RowSize,ColSize,iproc,fulldim1,fulldim2)
       CALL funccall(A,Afullblock,fulldim1,fulldim2,RowSize,ColSize,&
            & Abuffer,nbuffer,iproc,Global_I,Global_J,Buffer_I,Buffer_J,VI,VJ,&
            & localnrow2,bufferoffset,bufferproc,debug)

       CommunicationNodes(iproc)=.TRUE.
       localnrow(iproc) = RowSize
       localnrow3(iproc) = RowSize
       AddRowBlock(iproc)=.FALSE.
       localncol(iproc) = ColSize

       VI=VI+RowSize
       Global_I = Global_I + RowSize
       IF(nBlocksRow.GT.1)THEN
          RowSize = BLOCK_SIZE
          DO iBlocksRow = 2,nBlocksRow-1
!             print*,'iBlocksRow',iBlocksRow
             iproc = iproc + shift
             IF(iproc.GT.infpar%nodtot-1)iproc = iproc-infpar%nodtot

             Buffer_I = localnrow3(iproc)+1
             Buffer_J = 1

             !full row block and first collom Block
!             if(debug)call fillVisual(VISUAL,Global_I,Global_J,RowSize,ColSize,iproc,fulldim1,fulldim2)
!             print*,'C GI,GJ,BI,BJ',Global_I,Global_J,Buffer_I,Buffer_J,iproc
             CALL funccall(A,Afullblock,fulldim1,fulldim2,RowSize,ColSize,&
                  & Abuffer,nbuffer,iproc,Global_I,Global_J,Buffer_I,Buffer_J,VI,VJ,&
                  & localnrow2,bufferoffset,bufferproc,debug)

             CommunicationNodes(iproc)=.TRUE.
             localnrow(iproc) = localnrow(iproc) + RowSize             
             localnrow3(iproc) = localnrow3(iproc) + RowSize             
             AddRowBlock(iproc)=.FALSE.
             localncol(iproc) = ColSize

             VI=VI+RowSize
             Global_I = Global_I + RowSize
          ENDDO
          iBlocksRow = nBlocksRow
!          print*,'Final iBlocksRow',iBlocksRow
          RowSize = ndim1-VI+1
          iproc = iproc + shift
          IF(iproc.GT.infpar%nodtot-1)iproc = iproc-infpar%nodtot

          Buffer_I = localnrow3(iproc)+1
          Buffer_J = 1

          !last row block and first collom Block
!          if(debug)call fillVisual(VISUAL,Global_I,Global_J,RowSize,ColSize,iproc,fulldim1,fulldim2)
!          print*,'D GI,GJ,BI,BJ',Global_I,Global_J,Buffer_I,Buffer_J,iproc

          CALL funccall(A,Afullblock,fulldim1,fulldim2,RowSize,ColSize,&
               & Abuffer,nbuffer,iproc,Global_I,Global_J,Buffer_I,Buffer_J,VI,VJ,localnrow2,bufferoffset,bufferproc,debug)

          CommunicationNodes(iproc)=.TRUE.
          localnrow3(iproc) = localnrow3(iproc) + RowSize             
          localnrow(iproc) = localnrow(iproc) + RowSize             
          AddRowBlock(iproc)=.FALSE.
          localncol(iproc) = ColSize
       ENDIF
       IF(nBlocksCol.GT.1)THEN          
          VJ=VJ+Colsize
          Global_J = Global_J + Colsize
          ColSize = BLOCK_SIZE
          DO iBlocksCol = 2,nBlocksCol-1
!             print*,'iBlocksCol',iBlocksCol
             localncol3 = localncol
             AddRowBlock2=AddRowBlock
             AddBlock=.TRUE.
             iBlocksRow=1
!             print*,'START iBlocksRow',iBlocksRow
             localnrow3 = 0
             Global_I = insertrow!gridinfo_blocks(1,startblock)
             VI=1
             IF(nBlocksRow.EQ.1)THEN
                RowSize = ndim1-VI+1
             ELSE
                RowSize = startRowSize
             ENDIF
             iproc = gridinfo_blocks(5,iBlocksRow+iBlockOffset+(iBlocksCol+jBlockOffset-1)*gridinfo_nBlockssqrt)
             Buffer_I = 1
             Buffer_J = localncol3(iproc)+1
             !first row block and full collom Block
!             if(debug)call fillVisual(VISUAL,Global_I,Global_J,RowSize,ColSize,iproc,fulldim1,fulldim2)
!             print*,'E GI,GJ,BI,BJ',Global_I,Global_J,Buffer_I,Buffer_J,iproc
             CALL funccall(A,Afullblock,fulldim1,fulldim2,RowSize,ColSize,&
                  & Abuffer,nbuffer,iproc,Global_I,Global_J,Buffer_I,Buffer_J,&
                  & VI,VJ,localnrow2,bufferoffset,bufferproc,debug)
             CommunicationNodes(iproc)=.TRUE.
             localnrow3(iproc) = localnrow3(iproc) + RowSize             
             IF(AddRowBlock2(iproc))THEN
                localnrow(iproc) = localnrow(iproc) + RowSize             
                AddRowBlock(iproc)=.FALSE.
             ENDIF
             IF(AddBlock(iProc))THEN
                localncol(iproc) = localncol(iproc) + ColSize
                AddBlock(iProc)=.FALSE.
             ENDIF

             VI=VI+Rowsize
             Global_I = Global_I + Rowsize
             IF(nBlocksRow.GT.1)THEN
                Rowsize = BLOCK_SIZE
                DO iBlocksRow = 2,nBlocksRow-1
!                   print*,'iBlocksRow',iBlocksRow
                   iproc = iproc + shift
                   IF(iproc.GT.infpar%nodtot-1)iproc = iproc-infpar%nodtot

                   Buffer_I = localnrow3(iproc)+1
                   Buffer_J = localncol3(iproc)+1
                   !full row block and full collom Block
!                   if(debug)call fillVisual(VISUAL,Global_I,Global_J,RowSize,ColSize,iproc,fulldim1,fulldim2)
!                   print*,'F GI,GJ,BI,BJ',Global_I,Global_J,Buffer_I,Buffer_J,iproc
                   CALL funccall(A,Afullblock,fulldim1,fulldim2,RowSize,ColSize,&
                        & Abuffer,nbuffer,iproc,Global_I,Global_J,Buffer_I,Buffer_J,VI,VJ,localnrow2,bufferoffset,bufferproc,debug)
                   CommunicationNodes(iproc)=.TRUE.
                   localnrow3(iproc) = localnrow3(iproc) + RowSize             
                   IF(AddRowBlock2(iproc))THEN
                      localnrow(iproc) = localnrow(iproc) + RowSize             
                      AddRowBlock(iproc)=.FALSE.
                   ENDIF
                   IF(AddBlock(iProc))THEN
                      localncol(iproc) = localncol(iproc) + ColSize
                      AddBlock(iProc)=.FALSE.
                   ENDIF

                   VI=VI+Rowsize
                   Global_I = Global_I + Rowsize
                ENDDO
                iBlocksRow = nBlocksRow
!                print*,'FINAL iBlocksRow',iBlocksRow
                iproc = iproc + shift
                IF(iproc.GT.infpar%nodtot-1)iproc = iproc-infpar%nodtot
                Rowsize = ndim1-VI+1

                Buffer_I = localnrow3(iproc)+1
                Buffer_J = localncol3(iproc)+1
                !last row block and full collom Block
!                if(debug)call fillVisual(VISUAL,Global_I,Global_J,RowSize,ColSize,iproc,fulldim1,fulldim2)
!                print*,'G GI,GJ,BI,BJ',Global_I,Global_J,Buffer_I,Buffer_J,iproc
                CALL funccall(A,Afullblock,fulldim1,fulldim2,RowSize,ColSize,&
                     & Abuffer,nbuffer,iproc,Global_I,Global_J,Buffer_I,Buffer_J,VI,VJ,localnrow2,bufferoffset,bufferproc,debug)
                CommunicationNodes(iproc)=.TRUE.
                localnrow3(iproc) = localnrow3(iproc) + RowSize             
                IF(AddRowBlock2(iproc))THEN
                   localnrow(iproc) = localnrow(iproc) + RowSize             
                   AddRowBlock(iproc)=.FALSE.
                ENDIF
                IF(AddBlock(iProc))THEN
                   localncol(iproc) = localncol(iproc) + ColSize
                   AddBlock(iProc)=.FALSE.
                ENDIF

             ENDIF
             VJ=VJ+ColSize
             Global_J = Global_J + Colsize
          ENDDO
          iBlocksCol = nBlocksCol
!          print*,'FINAL iBlocksCol',iBlocksCol
          localncol3 = localncol
          AddRowBlock2=AddRowBlock
          AddBlock=.TRUE.
          ColSize = ndim2-VJ+1
          iBlocksRow=1
!          print*,'FIRST iBlocksRow',iBlocksRow
          localnrow3 = 0
          VI=1
          Global_I = insertrow!gridinfo_blocks(1,startblock) 
          IF(nBlocksRow.EQ.1)THEN
             RowSize = ndim1-VI+1
          ELSE
             RowSize = startRowSize
!             Rowsize = BLOCK_SIZE-MOD(TMP(1),BLOCK_SIZE)+1
!             RowSize = (BLOCK_SIZE-(TMP(1)-Global_I))
          ENDIF
          iproc = gridinfo_blocks(5,iBlocksRow+iBlockOffset+(iBlocksCol+jBlockOffset-1)*gridinfo_nBlockssqrt)

          Buffer_I = 1
          Buffer_J = localncol3(iproc)+1
!          print*,'Buffer_J',Buffer_J,'localncol3(iproc)',localncol3(iproc)
!          print*,'Global_I,Global_J',Global_I,Global_J
!          print*,'RowSize,ColSize',RowSize,ColSize
          !fist row block and last collom Block
!          if(debug)call fillVisual(VISUAL,Global_I,Global_J,RowSize,ColSize,iproc,fulldim1,fulldim2)
!          print*,'H GI,GJ,BI,BJ',Global_I,Global_J,Buffer_I,Buffer_J,iproc
          CALL funccall(A,Afullblock,fulldim1,fulldim2,RowSize,ColSize,&
               & Abuffer,nbuffer,iproc,Global_I,Global_J,Buffer_I,Buffer_J,VI,VJ,localnrow2,bufferoffset,bufferproc,debug)
          CommunicationNodes(iproc)=.TRUE.          
          localnrow3(iproc) = localnrow3(iproc) + RowSize             
          IF(AddRowBlock2(iproc))THEN
             localnrow(iproc) = localnrow(iproc) + RowSize             
             AddRowBlock(iproc)=.FALSE.
          ENDIF
          IF(AddBlock(iProc))THEN
             localncol(iproc) = localncol(iproc) + ColSize
             AddBlock(iProc)=.FALSE.
          ENDIF

          VI=VI+RowSize
          Global_I = Global_I + RowSize
          IF(nBlocksRow.GT.1)THEN
             Rowsize = BLOCK_SIZE
             DO iBlocksRow = 2,nBlocksRow-1
!                print*,'iBlocksRow',iBlocksRow
                iproc = iproc + shift
                IF(iproc.GT.infpar%nodtot-1)iproc = iproc-infpar%nodtot

                Buffer_I = localnrow3(iproc)+1
                Buffer_J = localncol3(iproc)+1
                !full row block and last collom Block
!                if(debug)call fillVisual(VISUAL,Global_I,Global_J,RowSize,ColSize,iproc,fulldim1,fulldim2)
!                print*,'I GI,GJ,BI,BJ',Global_I,Global_J,Buffer_I,Buffer_J,iproc
                CALL funccall(A,Afullblock,fulldim1,fulldim2,RowSize,ColSize,&
                     & Abuffer,nbuffer,iproc,Global_I,Global_J,Buffer_I,Buffer_J,VI,VJ,localnrow2,bufferoffset,bufferproc,debug)
                CommunicationNodes(iproc)=.TRUE.
                localnrow3(iproc) = localnrow3(iproc) + RowSize             
                IF(AddRowBlock2(iproc))THEN
                   localnrow(iproc) = localnrow(iproc) + RowSize             
                   AddRowBlock(iproc)=.FALSE.
                ENDIF
                IF(AddBlock(iProc))THEN
                   localncol(iproc) = localncol(iproc) + ColSize
                   AddBlock(iProc)=.FALSE.
                ENDIF

                VI=VI+RowSize
                Global_I = Global_I + RowSize
             ENDDO
             iBlocksRow = nBlocksRow
!             print*,'Final iBlocksRow',iBlocksRow
             RowSize = ndim1-VI+1
             iproc = iproc + shift
             IF(iproc.GT.infpar%nodtot-1)iproc = iproc-infpar%nodtot

             Buffer_I = localnrow3(iproc)+1
             Buffer_J = localncol3(iproc)+1
             !last row block and last collom Block
!             if(debug)call fillVisual(VISUAL,Global_I,Global_J,RowSize,ColSize,iproc,fulldim1,fulldim2)
!             print*,'J GI,GJ,BI,BJ',Global_I,Global_J,Buffer_I,Buffer_J,iproc
             CALL funccall(A,Afullblock,fulldim1,fulldim2,RowSize,ColSize,&
                  & Abuffer,nbuffer,iproc,Global_I,Global_J,Buffer_I,Buffer_J,VI,VJ,localnrow2,bufferoffset,bufferproc,debug)
             CommunicationNodes(iproc)=.TRUE.
             localnrow3(iproc) = localnrow3(iproc) + RowSize             
             IF(AddRowBlock2(iproc))THEN
                localnrow(iproc) = localnrow(iproc) + RowSize             
                AddRowBlock(iproc)=.FALSE.
             ENDIF
             IF(AddBlock(iProc))THEN
                localncol(iproc) = localncol(iproc) + ColSize
                AddBlock(iProc)=.FALSE.
             ENDIF

          ENDIF
       ENDIF
    ENDIF
!    do iproc=0,infpar%nodtot-1
!       print*,'iproc',iproc,'localnrow,localncol:',localnrow(iproc),localncol(iproc),CommunicationNodes(iproc)
!    enddo
!    if(debug)then
!       print*,'VISUAL'
!       call ls_output(VISUAL,1,ndim1,1,ndim2,ndim1,ndim2,1,6)
!    endif
    localnrow2 = localnrow
#endif
  end subroutine scalapack_transform_aux

!========================================================================================
! Scalapack_transform_aux (STA) funccalls
! STA_null
!           Does nothing
! STA_BuildAlocalFromFull
!           build up info in Alocal(iproc)%p for all procs 
!========================================================================================
#ifdef VAR_SCALAPACK

  subroutine fillVisual(VISUAL,GI,GJ,RowSize,ColSize,iproc,fulldim1,fulldim2)
    implicit none    
    integer, intent(in) :: RowSize,ColSize,iproc,fulldim1,fulldim2
    integer,intent(in)  :: Gi,Gj     !global row and collum index 
    real(Realk),optional :: VISUAL(fulldim1,fulldim2)
    logical :: debug
    !
    integer :: I,J,API,APJ,AFJ2,APJ2,ABJ2,OFFSET

    DO J=0,ColSize-1
       DO I=0,RowSize-1
          VISUAL(GI+I,GJ+J) = iproc
       ENDDO
    ENDDO
  end subroutine fillVisual
  
  !different funccalls
  subroutine STA_null(A,Afullblock,fulldim1,fulldim2,RowSize,ColSize,&
       & Abuffer,nbuffer,iproc,GI,GJ,ABI,ABJ,AFI,AFJ,ABndim1,bufferoffset,bufferproc,debug)
    implicit none    
    integer, intent(in) :: fulldim1,fulldim2,RowSize,ColSize,nbuffer,iproc,bufferproc
    integer,intent(in)  :: Gi,Gj     !global row and collum index 
    integer,intent(in)  :: ABI,ABJ !row and collum index in the buffer  
    integer,intent(in)  :: ABndim1(0:infpar%nodtot-1) !row dimension of buffer 
    integer,intent(in)  :: AFI,AFJ !row and collum index in the full Afullblock matrix
    integer,intent(in)  :: bufferoffset(0:infpar%nodtot-1)  !offset in Abuffer (0)
    real(Realk) :: Afullblock(fulldim1,fulldim2)
!    real(Realk),optional :: VISUAL(fulldim1,fulldim2)
    type(Matrix)        :: A 
    logical :: debug
    type(lsMatrix)      :: Abuffer(0:nbuffer)  
    !
    integer :: I,J,API,APJ,AFJ2,APJ2,ABJ2,OFFSET
    !do nothing
  end subroutine STA_null

  subroutine STA_BuildBufferFromFull(A,Afullblock,fulldim1,fulldim2,RowSize,ColSize,&
       & Abuffer,nbuffer,iproc,GI,GJ,ABI,ABJ,AFI,AFJ,ABndim1,bufferoffset,bufferproc,debug)
    implicit none    
    integer, intent(in) :: fulldim1,fulldim2,RowSize,ColSize,nbuffer,iproc,bufferproc
    integer,intent(in)  :: Gi,Gj     !global row and collum index 
    integer,intent(in)  :: ABI,ABJ !row and collum index in the buffer  
    integer,intent(in)  :: ABndim1(0:infpar%nodtot-1) !row dimension of buffer 
    integer,intent(in)  :: AFI,AFJ !row and collum index in the full Afullblock matrix
    integer,intent(in)  :: bufferoffset(0:infpar%nodtot-1)  !offset in Abuffer (0)
    real(Realk) :: Afullblock(fulldim1,fulldim2)
!    real(Realk),optional :: VISUAL(fulldim1,fulldim2)
    type(Matrix)        :: A 
    type(lsMatrix)      :: Abuffer(0:nbuffer)  
    logical :: debug
    !
    integer :: I,J,API,APJ,AFJ2,APJ2,ABJ2,OFFSET,boffset,ABndim2
    boffset = bufferoffset(iproc)
    ABndim2 = ABndim1(bufferproc)
    DO J=0,ColSize-1
       AFJ2 = AFJ+J
       ABJ2 = ABJ+J
       OFFSET = ABI+(ABJ2-1)*ABndim2 + boffset
       DO I=0,RowSize-1
          Abuffer(iproc)%elms(I+OFFSET) = AfullBlock(AFI+I,AFJ2)
       ENDDO
    ENDDO
  end  subroutine STA_BuildBufferFromFull

  subroutine STA_BuildAllBufferFromFull(A,Afullblock,fulldim1,fulldim2,RowSize,ColSize,&
       & Abuffer,nbuffer,iproc,GI,GJ,ABI,ABJ,AFI,AFJ,ABndim1,bufferoffset,bufferproc,debug)
    implicit none    
    integer, intent(in) :: fulldim1,fulldim2,RowSize,ColSize,nbuffer,iproc,bufferproc
    integer,intent(in)  :: Gi,Gj     !global row and collum index 
    integer,intent(in)  :: ABI,ABJ !row and collum index in the buffer  
    integer,intent(in)  :: ABndim1(0:infpar%nodtot-1) !row dimension of buffer 
    integer,intent(in)  :: AFI,AFJ !row and collum index in the full Afullblock matrix
    integer,intent(inout):: bufferoffset(0:infpar%nodtot-1)  !offset in Abuffer (0)
    real(Realk) :: Afullblock(fulldim1,fulldim2)
!    real(Realk),optional :: VISUAL(fulldim1,fulldim2)
    type(Matrix)        :: A 
    type(lsMatrix)      :: Abuffer(0:nbuffer)  
    logical :: debug
    !
    integer :: I,J,API,APJ,AFJ2,APJ2,ABJ2,OFFSET,boffset,ABndim2
    boffset = bufferoffset(iproc)
    ABndim2 = ABndim1(iproc)
    DO J=0,ColSize-1
       AFJ2 = AFJ+J
       ABJ2 = ABJ+J
       OFFSET = ABI+(ABJ2-1)*ABndim2 + boffset
       DO I=0,RowSize-1
!          print*,'ALL2 Abuffer(',iproc,')%elms(',&
!               & I+OFFSET,'), OFFSET',OFFSET,'I=',I,'boffset',boffset,'iproc',iproc
          Abuffer(iproc)%elms(I+OFFSET) = AfullBlock(AFI+I,AFJ2)
!          print*,'TOALL2 AfullBlock(',AFI+I,',',&
!               & AFJ2,')=',AfullBlock(AFI+I,AFJ2)
!          print*,'TOALL2 Abuffer(',iproc,')%elms(',&
!               & I+OFFSET,')=',Abuffer(iproc)%elms(I+OFFSET)
       ENDDO
    ENDDO
!    print*,'ALL bufferOffset(iproc)',bufferOffset(iproc)
!    bufferOffset(iproc) = bufferOffset(iproc) + RowSize*ColSize
!    print*,'ALL bufferOffset(',iproc,')=',bufferOffset(iproc)

  end  subroutine STA_BuildAllBufferFromFull

  subroutine STA_BuildFullFromBuffer(A,Afullblock,fulldim1,fulldim2,RowSize,ColSize,&
       & Abuffer,nbuffer,iproc,GI,GJ,ABI,ABJ,AFI,AFJ,ABndim1,bufferoffset,bufferproc,debug)!,VISUAL)
    implicit none    
    integer, intent(in) :: fulldim1,fulldim2,RowSize,ColSize,nbuffer,iproc,bufferproc
    integer,intent(in)  :: Gi,Gj     !global row and collum index 
    integer,intent(in)  :: ABI,ABJ !row and collum index in the buffer  
    integer,intent(in)  :: ABndim1(0:infpar%nodtot-1) !row dimension of buffer 
    integer,intent(in)  :: AFI,AFJ !row and collum index in the full Afullblock matrix
    integer,intent(in)  :: bufferoffset(0:infpar%nodtot-1)  !offset in Abuffer (0)
    real(Realk) :: Afullblock(fulldim1,fulldim2)
!    real(Realk),optional :: VISUAL(fulldim1,fulldim2)
    type(Matrix)        :: A 
    type(lsMatrix)      :: Abuffer(0:nbuffer)  
    logical :: debug
    !
    integer :: I,J,API,APJ,AFJ2,APJ2,ABJ2,OFFSET,boffset,ABndim2
    boffset = bufferoffset(iproc)
!    ABndim2 = ABndim1(bufferproc) !is this not wrong????
    ABndim2 = ABndim1(iproc)
    DO J=0,ColSize-1
       AFJ2 = AFJ+J
       ABJ2 = ABJ+J
       OFFSET = ABI+(ABJ2-1)*ABndim2 + boffset
       DO I=0,RowSize-1
          AfullBlock(AFI+I,AFJ2) = Abuffer(iproc)%elms(I+OFFSET)
!          print*,'AfullBlock(',AFI+I,',',AFJ2,')=',AfullBlock(AFI+I,AFJ2)
!          print*,'Abuffer(',iproc,')%elms(',I+OFFSET,')=',Abuffer(iproc)%elms(I+OFFSET)
       ENDDO
    ENDDO
  end  subroutine STA_BuildFullFromBuffer

  subroutine STA_BuildFullFromAllBuffer(A,Afullblock,fulldim1,fulldim2,RowSize,ColSize,&
       & Abuffer,nbuffer,iproc,GI,GJ,ABI,ABJ,AFI,AFJ,ABndim1,bufferoffset,bufferproc,debug)!,VISUAL)
    implicit none    
    integer, intent(in) :: fulldim1,fulldim2,RowSize,ColSize,nbuffer,iproc,bufferproc
    integer,intent(in)  :: Gi,Gj     !global row and collum index 
    integer,intent(in)  :: ABI,ABJ !row and collum index in the buffer  
    integer,intent(in)  :: ABndim1(0:infpar%nodtot-1) !row dimension of buffer 
    integer,intent(in)  :: AFI,AFJ !row and collum index in the full Afullblock matrix
    integer,intent(inout)  :: bufferoffset(0:infpar%nodtot-1)  !offset in Abuffer (0)
    real(Realk) :: Afullblock(fulldim1,fulldim2)
!    real(Realk),optional :: VISUAL(fulldim1,fulldim2)
    type(Matrix)        :: A 
    type(lsMatrix)      :: Abuffer(0:nbuffer)  
    logical :: debug
    !
    integer :: I,J,API,APJ,AFJ2,APJ2,ABJ2,OFFSET,boffset,ABndim2
    !iproc is the processor that owns this matrix block 
    IF(debug)print*,'iproc owns the matrix block',iproc
!    print*,'not used here : bufferproc',bufferproc
    boffset = bufferoffset(iproc)
    ABndim2 = ABndim1(iproc)
    IF(debug)print*,'ABndim2',ABndim2,'bufferoffset(iproc)',bufferoffset(iproc)
    DO J=0,ColSize-1
       AFJ2 = AFJ+J
       ABJ2 = ABJ+J
       OFFSET = ABI+(ABJ2-1)*ABndim2 + boffset
       IF(debug)print*,'OFFSET=',OFFSET
       IF(debug)print*,'ABI',ABI,'ABJ',ABJ,'ABndim2',ABndim2
       DO I=0,RowSize-1
          AfullBlock(AFI+I,AFJ2) = Abuffer(iproc)%elms(I+OFFSET)
          IF(debug)print*,'FROMALL AfullBlock(',AFI+I,',',AFJ2,')=',AfullBlock(AFI+I,AFJ2)
          IF(debug)print*,'FROMALL Abuffer(',iproc,')%elms(',I+OFFSET,')=',Abuffer(iproc)%elms(I+OFFSET)
       ENDDO
    ENDDO
!    print*,'ALL bufferOffset(iproc)',bufferOffset(iproc)
!    bufferOffset(iproc) = bufferOffset(iproc) + RowSize*ColSize
!    print*,'ALL bufferOffset(',iproc,')=',bufferOffset(iproc)

  end  subroutine STA_BuildFullFromAllBuffer

!!$  subroutine STA_BuildFullFromBuffer2(A,Afullblock,fulldim1,fulldim2,RowSize,ColSize,&
!!$       & Abuffer,nbuffer,iproc,GI,GJ,ABI,ABJ,AFI,AFJ,ABndim1,bufferoffset,bufferproc,debug)!,VISUAL)
!!$    implicit none    
!!$    integer, intent(in) :: fulldim1,fulldim2,RowSize,ColSize,nbuffer,iproc,bufferproc
!!$    integer,intent(in)  :: Gi,Gj     !global row and collum index 
!!$    integer,intent(in)  :: ABI,ABJ !row and collum index in the buffer  
!!$    integer,intent(in)  :: ABndim1(0:infpar%nodtot-1) !row dimension of buffer 
!!$    integer,intent(in)  :: AFI,AFJ !row and collum index in the full Afullblock matrix
!!$    integer,intent(in)  :: bufferoffset(0:infpar%nodtot-1)  !offset in Abuffer (0)
!!$    real(Realk) :: Afullblock(fulldim1,fulldim2)
!!$!    real(Realk),optional :: VISUAL(fulldim1,fulldim2)
!!$    type(Matrix)        :: A 
!!$    type(lsMatrix)      :: Abuffer(0:nbuffer)  
!!$    !
!!$    integer :: I,J,API,APJ,AFJ2,APJ2,ABJ2,OFFSET,boffset,ABndim2
!!$    boffset = bufferoffset(bufferproc)
!!$    ABndim2 = ABndim1(bufferproc)
!!$    DO J=0,ColSize-1
!!$       AFJ2 = AFJ+J
!!$       ABJ2 = ABJ+J
!!$       OFFSET = ABI+(ABJ2-1)*ABndim2 + boffset
!!$       DO I=0,RowSize-1
!!$          AfullBlock(AFI+I,AFJ2) = Abuffer(bufferproc)%elms(I+OFFSET)
!!$          if(scalapack_mynum.EQ.1)print*,'AfullBlock(',AFI+I,',',AFJ2,')=',AfullBlock(AFI+I,AFJ2)
!!$          if(scalapack_mynum.EQ.1)print*,'Abuffer(',iproc,')%elms(',I+OFFSET,')=',Abuffer(iproc)%elms(I+OFFSET)
!!$       ENDDO
!!$    ENDDO
!!$  end  subroutine STA_BuildFullFromBuffer2

  subroutine STA_BuildFullFromSingleBuffer(A,Afullblock,fulldim1,fulldim2,RowSize,ColSize,&
       & Abuffer,nbuffer,iproc,GI,GJ,ABI,ABJ,AFI,AFJ,ABndim1,bufferoffset,bufferproc,debug)!,VISUAL)
    implicit none    
    integer, intent(in) :: fulldim1,fulldim2,RowSize,ColSize,nbuffer,iproc,bufferproc
    integer,intent(in)  :: Gi,Gj     !global row and collum index 
    integer,intent(in)  :: ABI,ABJ !row and collum index in the buffer  
    integer,intent(in)  :: ABndim1(0:infpar%nodtot-1) !row dimension of buffer 
    integer,intent(in)  :: AFI,AFJ !row and collum index in the full Afullblock matrix
    integer,intent(in)  :: bufferoffset(0:infpar%nodtot-1)  !offset in Abuffer (0)
    real(Realk) :: Afullblock(fulldim1,fulldim2)
!    real(Realk),optional :: VISUAL(fulldim1,fulldim2)
    type(Matrix)        :: A 
    logical :: debug
    type(lsMatrix)      :: Abuffer(0:nbuffer)  
    !
    integer :: I,J,API,APJ,AFJ2,APJ2,ABJ2,OFFSET,boffset
    IF(iproc.EQ.scalapack_mynum)THEN 
       !iproc is the node that has the matrix block
       !bufferproc the buffer index that has the contribution to the 
       !replace iproc with bufferproc to retrieve from Abuffer(bufferproc) instead of Abuffer(iproc)
       CALL STA_BuildFullFromBuffer(A,Afullblock,fulldim1,fulldim2,RowSize,ColSize,&
            & Abuffer,nbuffer,bufferproc,GI,GJ,ABI,ABJ,AFI,AFJ,ABndim1,bufferoffset,bufferproc,debug)
    ENDIF

  end  subroutine STA_BuildFullFromSingleBuffer

  subroutine STA_BuildAlocalFromBuffer(A,Afullblock,fulldim1,fulldim2,RowSize,ColSize,&
       & Abuffer,nbuffer,iproc,GI,GJ,ABI,ABJ,AFI,AFJ,ABndim1,bufferoffset,bufferproc,debug)
    implicit none    
    integer, intent(in) :: fulldim1,fulldim2,RowSize,ColSize,nbuffer,iproc,bufferproc
    integer,intent(in)  :: Gi,Gj   !global row and collum index 
    integer,intent(in)  :: ABI,ABJ !row and collum index in the buffer  
    integer,intent(in)  :: ABndim1(0:infpar%nodtot-1) !row dimension of buffer 
    integer,intent(in)  :: AFI,AFJ !row and collum index in the full Afullblock matrix
    integer,intent(in)  :: bufferoffset(0:infpar%nodtot-1)  !offset in Abuffer (0)
    real(Realk) :: Afullblock(fulldim1,fulldim2)
!    real(Realk),optional :: VISUAL(fulldim1,fulldim2)
    type(Matrix)        :: A 
    type(lsMatrix)      :: Abuffer(0:nbuffer)  
    logical :: debug
    !
    integer :: I,J,API,APJ,AFJ2,APJ2,ABJ2,OFFSET,boffset,ABndim2
    IF(iproc.EQ.scalapack_mynum)THEN 
       !if the owner of the block (iproc) is the same as my rank (if I own the block)
       boffset = bufferoffset(bufferproc)
       ABndim2 = ABndim1(bufferproc)
       API = INDXG2L(GI,BLOCK_SIZE,0,0,slgrid%nprow)  !local index 
       APJ = INDXG2L(GJ,BLOCK_SIZE,0,0,slgrid%npcol)  !local index 
       DO J=0,ColSize-1
          APJ2 = APJ+J
          ABJ2 = ABJ+J
          OFFSET = ABI+(ABJ2-1)*ABndim2 + boffset
          DO I=0,RowSize-1
             A%p(API+I,APJ2) = Abuffer(bufferproc)%elms(I+OFFSET)
!             if(scalapack_mynum.EQ.1)print*,'A%p(',API+I,',',APJ2,')=',A%p(API+I,APJ2)
!             if(scalapack_mynum.EQ.1)print*,'Abuffer(',bufferproc,')%elms(',I+OFFSET,')=',Abuffer(bufferproc)%elms(I+OFFSET)
          ENDDO
       ENDDO
    ENDIF
  end subroutine STA_BuildAlocalFromBuffer

  subroutine STA_BuildBufferFromAlocal(A,Afullblock,fulldim1,fulldim2,RowSize,ColSize,&
       & Abuffer,nbuffer,iproc,GI,GJ,ABI,ABJ,AFI,AFJ,ABndim1,bufferoffset,bufferproc,debug)
    implicit none    
    integer, intent(in) :: fulldim1,fulldim2,RowSize,ColSize,nbuffer,iproc,bufferproc
    integer,intent(in)  :: Gi,Gj     !global row and collum index 
    integer,intent(in)  :: ABI,ABJ !row and collum index in the buffer  
    integer,intent(in)  :: ABndim1(0:infpar%nodtot-1) !row dimension of buffer 
    integer,intent(in)  :: AFI,AFJ !row and collum index in the full Afullblock matrix
    integer,intent(in)  :: bufferoffset(0:infpar%nodtot-1) 
    real(Realk) :: Afullblock(fulldim1,fulldim2)
    logical :: debug
!    real(Realk),optional :: VISUAL(fulldim1,fulldim2)
    type(Matrix)        :: A 
    type(lsMatrix)      :: Abuffer(0:nbuffer)  
    !
    integer :: I,J,API,APJ,AFJ2,APJ2,ABJ2,OFFSET,boffset,ABndim2
!    print*,'iproc',iproc,'scalapack_mynum',scalapack_mynum
    IF(iproc.EQ.scalapack_mynum)THEN 
       boffset = bufferoffset(bufferproc)
!       print*,'bufferoffset(bufferproc)',bufferoffset(bufferproc)
!       print*,'RowSize',RowSize,'ColSize',ColSize
       ABndim2 = ABndim1(bufferproc)
       !if the owner of the block (iproc) is the same as my rank (if I own the block)
       API = INDXG2L(GI,BLOCK_SIZE,0,0,slgrid%nprow)  !local index 
       APJ = INDXG2L(GJ,BLOCK_SIZE,0,0,slgrid%npcol)  !local index 
       DO J=0,ColSize-1
          APJ2 = APJ+J
          ABJ2 = ABJ+J
          OFFSET = ABI+(ABJ2-1)*ABndim2 + boffset
          DO I=0,RowSize-1
             Abuffer(bufferproc)%elms(I+OFFSET) = A%p(API+I,APJ2)
!             print*,'TK1 Abuffer(',bufferproc,')%elms(',I+OFFSET,')=',Abuffer(bufferproc)%elms(I+OFFSET)
!             print*,'TK1 A%p(',API+I,',',APJ2,')=',A%p(API+I,APJ2)
          ENDDO
       ENDDO
!    ELSE
!       print*,'no cont iproc',iproc,'scalapack_mynum',scalapack_mynum
    ENDIF
  end  subroutine STA_BuildBufferFromAlocal

#endif
!========================================================================================
  subroutine mat_scalapack_setlowertriangular_zero(A)
    implicit none
    type(Matrix) :: A  
    INTEGER      :: DESC_A(9)
#ifdef VAR_SCALAPACK
    CALL PDM_SYNC(Job_setlowertriangular_zero,A)
    CALL PDM_DSCINIT(DESC_A, A)
    call scalapack_setlowertriangular_zero_aux(A,DESC_A)
#endif
  end subroutine mat_scalapack_setlowertriangular_zero

  subroutine scalapack_setlowertriangular_zero_aux(A,DESC_A)
    implicit none
    type(Matrix) :: A  
    INTEGER      :: DESC_A(9)
#ifdef VAR_SCALAPACK
    integer :: imyblock,iBlock,iglobalstartrow,iglobalstartcol,API,APJ,RowSize,ColSize
    integer :: j,i,apj2,GI,GJ
    do imyblock=1,gridinfo_mynblocks
       iBlock = gridinfo_myblocks(imyblock)
       iglobalstartrow = gridinfo_blocks(1,iBlock)
       iglobalstartcol = gridinfo_blocks(2,iBlock)
       API = INDXG2L(iglobalstartrow,BLOCK_SIZE,0,0,slgrid%nprow)  !local index 
       APJ = INDXG2L(iglobalstartcol,BLOCK_SIZE,0,0,slgrid%npcol)  !local index 
       RowSize = MIN(BLOCK_SIZE,(A%nrow-(iglobalstartrow-1)))
       ColSize = MIN(BLOCK_SIZE,(A%ncol-(iglobalstartcol-1)))
       DO J=0,ColSize-1
          APJ2 = APJ+J
          GJ = INDXL2G( APJ2, BLOCK_SIZE, SLgrid%mycol, 0, SLgrid%npcol ) 
          DO I=0,RowSize-1
             GI = INDXL2G( API+I, BLOCK_SIZE, SLgrid%myrow, 0, SLgrid%nprow ) 
             IF(GI.GT.GJ)THEN
                A%p(API+I,APJ2) = 0.0E0_realk
             ENDIF
          ENDDO
       ENDDO
    enddo
#endif
  end subroutine scalapack_setlowertriangular_zero_aux

  subroutine mat_scalapack_scal_dia(alpha,A)
    implicit none
    type(Matrix) :: A  
    INTEGER      :: DESC_A(9)
    real(realk) :: alpha
#ifdef VAR_SCALAPACK
    CALL PDM_SYNC(Job_scal_dia,A)
    CALL PDM_DSCINIT(DESC_A, A)
    call ls_mpibcast(alpha,infpar%master,scalapack_comm)
    call scalapack_scal_dia_aux(alpha,A,DESC_A)
#endif
  end subroutine mat_scalapack_scal_dia

  subroutine scalapack_scal_dia_aux(alpha,A,DESC_A)
    implicit none
    type(Matrix) :: A  
    INTEGER      :: DESC_A(9)
    real(realk) :: alpha
#ifdef VAR_SCALAPACK
    integer :: imyblock,iBlock,iglobalstartrow,iglobalstartcol,API,APJ,RowSize,ColSize
    integer :: j,i,apj2,GI,GJ
    do imyblock=1,gridinfo_mynblocks
       iBlock = gridinfo_myblocks(imyblock)
       iglobalstartrow = gridinfo_blocks(1,iBlock)
       iglobalstartcol = gridinfo_blocks(2,iBlock)
       API = INDXG2L(iglobalstartrow,BLOCK_SIZE,0,0,slgrid%nprow)  !local index 
       APJ = INDXG2L(iglobalstartcol,BLOCK_SIZE,0,0,slgrid%npcol)  !local index 
       RowSize = MIN(BLOCK_SIZE,(A%nrow-(iglobalstartrow-1)))
       ColSize = MIN(BLOCK_SIZE,(A%ncol-(iglobalstartcol-1)))
       DO J=0,ColSize-1
          APJ2 = APJ+J
          GJ = INDXL2G( APJ2, BLOCK_SIZE, SLgrid%mycol, 0, SLgrid%npcol ) 
          DO I=0,RowSize-1
             GI = INDXL2G( API+I, BLOCK_SIZE, SLgrid%myrow, 0, SLgrid%nprow ) 
             IF(GI.EQ.GJ)THEN
                A%p(API+I,APJ2) = alpha*A%p(API+I,APJ2)
             ENDIF
          ENDDO
       ENDDO
    enddo
#endif
  end subroutine scalapack_scal_dia_aux

  subroutine mat_scalapack_scal_dia_vec(alpha,A,ndim)
    implicit none
    integer      :: ndim
    type(Matrix) :: A  
    INTEGER      :: DESC_A(9)
    real(realk) :: alpha(ndim)
#ifdef VAR_SCALAPACK
    IF(ndim.NE.A%nrow)call lsquit('mat_scalapack_scal_dia_vec error',-1)
    CALL PDM_SYNC(Job_scal_dia_vec,A)
    CALL PDM_DSCINIT(DESC_A, A)
    call ls_mpibcast(alpha,ndim,infpar%master,scalapack_comm)
    call scalapack_scal_dia_vec_aux(alpha,A,DESC_A,ndim)
#endif
  end subroutine mat_scalapack_scal_dia_vec

  subroutine scalapack_scal_dia_vec_aux(alpha,A,DESC_A,ndim)
    implicit none
    integer      :: ndim
    type(Matrix) :: A  
    INTEGER      :: DESC_A(9)
    real(realk) :: alpha(ndim)
#ifdef VAR_SCALAPACK
    integer :: imyblock,iBlock,iglobalstartrow,iglobalstartcol,API,APJ,RowSize,ColSize
    integer :: j,i,apj2,GI,GJ
    do imyblock=1,gridinfo_mynblocks
       iBlock = gridinfo_myblocks(imyblock)
       iglobalstartrow = gridinfo_blocks(1,iBlock)
       iglobalstartcol = gridinfo_blocks(2,iBlock)
       API = INDXG2L(iglobalstartrow,BLOCK_SIZE,0,0,slgrid%nprow)  !local index 
       APJ = INDXG2L(iglobalstartcol,BLOCK_SIZE,0,0,slgrid%npcol)  !local index 
       RowSize = MIN(BLOCK_SIZE,(A%nrow-(iglobalstartrow-1)))
       ColSize = MIN(BLOCK_SIZE,(A%ncol-(iglobalstartcol-1)))
       DO J=0,ColSize-1
          APJ2 = APJ+J
          GJ = INDXL2G( APJ2, BLOCK_SIZE, SLgrid%mycol, 0, SLgrid%npcol ) 
          DO I=0,RowSize-1
             GI = INDXL2G( API+I, BLOCK_SIZE, SLgrid%myrow, 0, SLgrid%nprow ) 
             IF(GI.EQ.GJ)THEN
                A%p(API+I,APJ2) = alpha(GI)*A%p(API+I,APJ2)
             ENDIF
          ENDDO
       ENDDO
    enddo
#endif
  end subroutine scalapack_scal_dia_vec_aux

!> \brief See mat_scal in mat-operations.f90
  subroutine mat_scalapack_scal(alpha,a)
    real(realk), intent(in) :: alpha
    type(Matrix), intent(inout) :: A
    call PDM_SCAL(ALPHA,A)
  end subroutine mat_scalapack_scal

!> \brief See mat_copy in mat-operations.f90
  subroutine mat_scalapack_copy(alpha,a,b)
     implicit none
    real(realk), intent(in) :: alpha
    TYPE(Matrix), INTENT(IN)    :: a
    TYPE(Matrix), INTENT(INOUT) :: b
    INTEGER      :: DESC_AF(9), DESC_A(9), I,nsize
    CALL mat_scalapack_assign(B,A)
    if (alpha.ne. 1.0E0_realk) call mat_scalapack_scal(alpha,b)
  end subroutine mat_scalapack_copy

  !> \brief See mat_assign in mat-operations.f90
  subroutine mat_scalapack_assign(a,b)
     implicit none
     TYPE(Matrix), INTENT(INOUT) :: a
     TYPE(Matrix), INTENT(IN)    :: b
     CALL PDM_ASSIGN(A,B)
   end subroutine mat_scalapack_assign

  !> \brief See mat_trans in mat-operations.f90
  subroutine mat_scalapack_trans(a,b)
     implicit none
     TYPE(Matrix), INTENT(IN)    :: a
     TYPE(Matrix), INTENT(INOUT) :: b
     CALL PDM_TRANS(A,B)
   end subroutine mat_scalapack_trans

  !> \brief See mat_daxpy in mat-operations.f90
  subroutine mat_scalapack_daxpy(alpha,A,B)
     implicit none
     real(realk),intent(in)      :: alpha
     TYPE(Matrix), INTENT(IN)    :: A
     TYPE(Matrix), INTENT(INOUT) :: B
     !
     INTEGER      :: DESC_A(DLEN_), DESC_B(DLEN_)
     REAL(REALK)  :: AB(2)
#ifdef VAR_SCALAPACK

     CALL PDM_SYNC(Job_daxpy,A,B)

     CALL PDM_DSCINIT(DESC_A,A)
     CALL PDM_DSCINIT(DESC_B,B)

     AB(1)=ALPHA; AB(2)=1.0E0_realk
     CALL DGEBS2D(SLGrid%ictxt,'A','I',2,1,AB,2)
     CALL PDMATADD(A%nrow,A%ncol,ALPHA,A%p,1,1,DESC_A,&
          &        1.0E0_realk, B%p,1,1,DESC_B )
#endif
   end subroutine mat_scalapack_daxpy


  !> \brief print the local part of the matrix. This routine have no MPI sync
  subroutine mat_scalapack_print_local(a)
     implicit none
     TYPE(Matrix), intent(IN) :: a
     INTEGER      :: DESC_A(DLEN_)
#ifdef VAR_SCALAPACK
     CALL PDM_DSCINIT(DESC_A,A)
     write(*,'("I am number ",I3," at postion (",I3,",",I3,")")')scalapack_mynum,SLGrid%myrow,SLGrid%mycol
     write(*,'("and this is my block:A%localnrow,A%localncol:",I3,",",I3)')A%localnrow,A%localncol
     IF(A%localnrow*A%localncol.GT.0)THEN
        call ls_output(A%p,1,A%localnrow,1,A%localncol,A%localnrow,A%localncol,1,6)
     ELSE
        write(*,*)'I have no block'
     ENDIF
#endif
   end subroutine mat_scalapack_print_local

  !> \brief all nodes print their local parts of the matrix
   subroutine mat_scalapack_print_global(a)
     implicit none
     TYPE(Matrix), intent(IN) :: a
     INTEGER      :: DESC_A(DLEN_)
#ifdef VAR_SCALAPACK
     CALL PDM_SYNC(Job_print_global,A)
     call mat_scaLapack_print_local(A)
     call lsmpi_barrier(scalapack_comm)
#endif
   end subroutine mat_scalapack_print_global

  !> \brief See mat_daxpy in mat-operations.f90
  subroutine mat_scalapack_add(alpha, a, beta, b, c)
     implicit none
     TYPE(Matrix), intent(IN) :: a, b
     REAL(realk), INTENT(IN)  :: alpha, beta
     TYPE(Matrix)             :: c
     !
     IF(beta.EQ. 0E0_realk)then
        call mat_scalapack_copy(alpha,A,C)
     ELSE
        call mat_scalapack_copy(beta,B,C)
        IF(alpha.NE. 0E0_realk)call mat_scalapack_daxpy(alpha,A,C)
     ENDIF
!     call mat_scalapack_copy(1E0_realk,B,C)
!     CALL PDM_SYNC(Job_daxpy,A,C)
!     CALL PDM_DSCINIT(DESC_ABC,A)
!     AB(1)=ALPHA; AB(2)=BETA
!     CALL DGEBS2D(SLGrid%ictxt,'A','I',2,1,AB,2)
!     CALL PDMATADD(A%nrow,A%ncol,ALPHA,A%p,1,1,DESC_ABC,&
!          &        BETA, C%p,1,1,DESC_ABC )

   end subroutine mat_scalapack_add

  !> \brief See mat_mul in mat-operations.f90
   subroutine mat_scalapack_mul(a, b, ta, tb, alpha, beta, c)
     implicit none
     TYPE(Matrix), intent(IN) :: a, b
     character, intent(in)    :: ta, tb
     REAL(realk), INTENT(IN)  :: alpha, beta
     TYPE(Matrix), intent(inout):: c
     !
     INTEGER      :: DESC_A(9),DESC_B(9),DESC_C(9)
     REAL(REALK)  :: AB(2)
     CHARACTER    :: T(2)
     integer      :: m,k,n
     if (ta == 'n' .or. ta == 'N') then
       m = a%nrow
       k = a%ncol
     elseif (ta == 't' .or. ta == 'T') then
       m = a%ncol
       k = a%nrow
     endif
     
     if (tb == 'n' .or. tb == 'N') then
       n = b%ncol
     elseif (tb == 't' .or. tb == 'T') then
       n = b%nrow
     endif

#ifdef VAR_SCALAPACK

     CALL PDM_SYNC(Job_mul,A,B,C)
     
     CALL PDM_DSCINIT(DESC_A,A)
     CALL PDM_DSCINIT(DESC_B,B)
     CALL PDM_DSCINIT(DESC_C,C)
     
     AB(1)=ALPHA; AB(2)=BETA
     CALL DGEBS2D(SLGrid%ictxt,'A','I',2,1,AB,2)
     
     T(1) = TA; T(2) = TB
     call ls_mpibcast(T,2,infpar%master,scalapack_comm)
     CALL PDGEMM(TA,TB,m,n,k,ALPHA,A%p,1,1,DESC_A,&
          B%p,1,1,DESC_B,BETA,C%p,1,1,DESC_C)
#endif     
   end subroutine mat_scalapack_mul
   !> \brief calculate the outer product of two vectors and store them in a matrix
   !> \autor Patrick Ettenhuber
   !> \date May, 2012
   subroutine mat_scalapack_dger(alpha,x,y,A)
     implicit none
     TYPE(Matrix), intent(IN) :: A
     REAL(realk), INTENT(IN)  :: alpha
     real(realk) :: x(A%nrow),y(A%ncol)
     INTEGER      :: DESC_A(9),DESC_X(9),DESC_Y(9),ierr,LLD,NUMROC
     REAL(REALK)  :: AB(2)
     EXTERNAL NUMROC
#ifdef VAR_SCALAPACK
call lsquit('TKTKTKTKTKTKTKTKTKTK',-1)
!     CALL PDM_SYNC(Job_dger,A) 
!     AB(1)=alpha
!     AB(2)=0.0E0_realk
!     call PDM_DSCINIT(DESC_A,A)
!     CALL DGEBS2D(SLGrid%ictxt,'A','I',2,1,AB,2)
!     LLD=MAX(1,NUMROC(A%nrow,A%nrow,0,0,1))
!     CALL DESCINIT(DESC_X,A%nrow,1,A%nrow,1,0,0,SLGrid%ictxt,LLD,ierr)
!     LLD=MAX(1,NUMROC(A%ncol,A%ncol,0,0,1))
!     CALL DESCINIT(DESC_Y,A%ncol,1,A%ncol,1,0,0,SLGrid%ictxt,LLD,ierr)
!     CALL IGEBS2D(SLGrid%ictxt,'A','I',9,1,DESC_X,9)
!     CALL IGEBS2D(SLGrid%ictxt,'A','I',9,1,DESC_Y,9)
!     call pdger(A%nrow,A%ncol,AB(1),x,1,1,DESC_X,1,y,1,1,DESC_Y,1,A%p,1,1,DESC_A,ierr)
#endif     
   end subroutine mat_scalapack_dger
   
   !> \brief calculate the hadamard product for pdms
   !> \autor Patrick Ettenhuber
   !> \date May, 2012
   subroutine mat_scalapack_hmul(alpha,a,b,beta,c)
     implicit none
     TYPE(Matrix), intent(IN) :: a,b
     REAL(realk), INTENT(IN)  :: alpha,beta
     TYPE(Matrix), intent(inout):: c
     INTEGER      :: DESC_A(9),DESC_B(9),DESC_C(9)
     REAL(REALK)  :: AB(2)
     logical :: check
#ifdef VAR_SCALAPACK
     !check the sizes of the matrices
     check =  (((A%nrow/=B%nrow).or.(A%ncol/=B%ncol)).or.((B%nrow/=C%nrow).or.(B%ncol/=C%ncol)))
     if (check) then
       print *,"salapack hadamard prodcut (hmul), matrix sizes are not matching to each other"
       call lsquit("salapack hadamard prodcut (hmul), matrix sizes are not matching to each other")
     endif
     CALL PDM_SYNC(Job_hmul,A,B,C) 
     AB(1)=alpha
     AB(2)=beta
     CALL DGEBS2D(SLGrid%ictxt,'A','I',2,1,AB,2)
     C%p=alpha*A%p*B%p+beta*C%p
#endif     
   end subroutine mat_scalapack_hmul

   subroutine mat_scalapack_hdiv(a,b,mu)
     implicit none
     TYPE(Matrix), intent(inout):: a
     TYPE(Matrix), intent(IN) :: b
     REAL(realk), INTENT(IN)  :: mu
     REAL(REALK)  :: AB(2), denom
     integer      :: i,j
#ifdef VAR_SCALAPACK
     CALL PDM_SYNC(Job_hdiv,A,B) 
     AB(1)=mu
     CALL DGEBS2D(SLGrid%ictxt,'A','I',2,1,AB,2)
     do j=1,A%localncol
      do i=1,A%localnrow
      denom = B%p(i,j) - mu
      if (denom.ne.0.0_realk) A%p(i,j) = A%p(i,j)/denom
      enddo
     enddo
#endif     
   end subroutine mat_scalapack_hdiv


   !> \brief matrix multiply with one diagonal matrix
   subroutine mat_scalapack_dmul(a,b,tb,alpha,beta,c)
     implicit none
     real(realk) :: a(:)
     TYPE(Matrix), intent(IN) :: b
     character, intent(in)    :: tb
     REAL(realk), INTENT(IN)  :: alpha, beta
     TYPE(Matrix), intent(inout):: c
     TYPE(Matrix) :: BT
     INTEGER      :: DESC_B(9),DESC_C(9),DESC_BT(9),a_dim,i,j,loc,glb,l,max_block_idx
     REAL(REALK)  :: AB(2),val
     CHARACTER    :: T(2)
     logical :: check
#ifdef VAR_SCALAPACK
     !check the sizes of the matrices
     check = (((tb=='n').or.(tb=='N')).and.((B%nrow/=C%nrow).or.(B%ncol/=C%ncol)))
     check = check.or.(((tb=='t').or.(tb=='T')).and.((B%ncol/=C%nrow).or.(B%nrow/=C%ncol)))
     if (check) then
       print *,"salapack dmul, matrix sizes are not matching to each other with current transp"
       call lsquit("salapack dmul, matrix sizes are not matching to each other with current transp")
     endif

     ! create storing matrix for transposition and dummy if not transposed
     if ((tb=='t').or.(tb=='T'))then
      call mat_scalapack_init(BT,C%nrow,C%ncol)
     else
      call mat_scalapack_init(BT,1,1)
     endif
     !synchronize the slaves
     CALL PDM_SYNC(Job_dmul,B,C,BT) 
     AB(1)=alpha
     AB(2)=beta
     T(1)=tb
     T(2)='X'
     CALL DGEBS2D(SLGrid%ictxt,'A','I',2,1,AB,2)
     CALL SGEBS2D(SLGrid%ictxt,'A','I',2,1,T,2)
     !CALL DGEBS2D(SLGrid%ictxt,'A','I',B%nrow,1,a,B%nrow)
     call ls_mpibcast(a,B%nrow,infpar%master,SLGrid%comm)    
 
     !scale c if needed
     if ((C%localnrow*C%localncol.NE. 0) .and. (AB(2)/= 1.0E0_realk))THEN
       CALL DSCAL(C%localnrow*C%localncol,AB(2),C%p,1)
     endif

     !transpose in a halfway efficient way, thus it is recommended to rewrite the
     !code to avoid transpositions
     if ((tb=='t').or.(tb=='T')) then
       CALL PDM_DSCINIT(DESC_B,B)
       CALL PDM_DSCINIT(DESC_BT,BT)
       CALL PDTRAN(BT%nrow,BT%ncol,1.0E0_realk,B%p,1,1,DESC_B,&
           & 0.0E0_realk,BT%p,1,1,DESC_BT)
     endif
     ! calculate the multiplication locally
     l=0
     do i=1,C%localnrow,BLOCK_SIZE
       if(i*BLOCK_SIZE > C%localnrow)then
         max_block_idx = mod(C%localnrow,BLOCK_SIZE)
       else
         max_block_idx=BLOCK_SIZE
       endif
       do j=1,max_block_idx
         glb=(l*SLGrid%nprow+SLGrid%myrow)*BLOCK_SIZE+j
         loc=l*BLOCK_SIZE + j
         if((T(1)=='n') .or. (T(1)=='N'))then
           CALL daxpy(C%localncol,a(glb)*AB(1),B%p(loc,1),B%localnrow,C%p(loc,1),C%localnrow)
         else if((T(1)=='t').or.(T(1)=='T')) then
           CALL daxpy(C%localncol,a(glb)*AB(1),BT%p(loc,1),BT%localnrow,C%p(loc,1),C%localnrow)
         endif
       enddo
       l=l+1
     enddo
     call mat_scalapack_free(BT)
#endif     
   end subroutine mat_scalapack_dmul

  !> \brief See mat_density_from_orbs in mat-operations.f90
   subroutine mat_scalapack_density_from_orbs(a, b,ndim,nocc)
     implicit none
     TYPE(Matrix), intent(IN) :: a
     TYPE(Matrix), intent(inout):: b
     integer,intent(in) :: ndim,nocc
     !
     INTEGER      :: DESC_A(9),DESC_B(9)
#ifdef VAR_SCALAPACK

     CALL PDM_SYNC(Job_density_from_orbs,A,B)
     CALL PDM_DSCINIT(DESC_A,A)
     CALL PDM_DSCINIT(DESC_B,B)

     call ls_mpibcast(nocc,infpar%master,scalapack_comm)
     CALL PDGEMM('N','T',A%nrow,A%ncol,nocc,1E0_realk,A%p,1,1,DESC_A,&
          A%p,1,1,DESC_A,0E0_realk,B%p,1,1,DESC_B)
#endif     
   end subroutine mat_scalapack_density_from_orbs

!> \brief See mat_tr in mat-operations.f90
  function mat_scalapack_tr(A)
     implicit none
     TYPE(Matrix), intent(IN) :: A
     REAL(realk) :: mat_scalapack_tr 
     INTEGER      :: DESC_A(9),I
#ifdef VAR_SCALAPACK
     REAL(REALK)      :: PDLATRA, OUTP
     EXTERNAL     :: PDLATRA

     CALL PDM_SYNC(Job_trace,A)
     CALL PDM_DSCINIT(DESC_A,A)
     mat_scalapack_tr = PDLATRA(A%nrow,A%p,1,1,DESC_A)
#else
     mat_scalapack_tr = 0
#endif
   end function mat_scalapack_tr

!> \brief See mat_trAB in mat-operations.f90
  function mat_scalapack_trAB(A,B)
     implicit none
     TYPE(Matrix), intent(IN) :: A,B
     TYPE(Matrix) :: BT
     REAL(realk)  :: mat_scalapack_trAB 
     CALL mat_scalapack_init(BT,A%ncol,A%nrow)
     CALL mat_scalapack_trans(B,BT)
     mat_scalapack_trAB = mat_scalapack_DOTPRODUCT(A,BT)
     CALL mat_scalapack_free(BT)
   end function mat_scalapack_trAB

!> \brief See mat_dotproduct in mat-operations.f90
  function mat_scalapack_dotproduct(A,B)
     implicit none
     TYPE(Matrix), intent(IN) :: A,B
     REAL(REALK) :: mat_scalapack_dotproduct
     real(realk) :: DOT
     real(realk), external :: ddot
     DOT = 0.0E0_realk
#ifdef VAR_SCALAPACK
     CALL PDM_SYNC(Job_dotproduct,A,B)
     IF (A%localnrow*A%localncol.NE. 0) THEN
        DOT = DDOT(A%localnrow*A%localncol, A%p,1, B%p,1)
     ENDIF
     CALL DGSUM2D( SLGrid%ictxt,'A','I',1,1,DOT,1,0,0)
#endif
     mat_scalapack_dotproduct = DOT
   end function mat_scalapack_dotproduct

!> \brief See mat_trAB in mat-operations.f90
  function mat_scalapack_sqnorm2(A)
     implicit none
     TYPE(Matrix), intent(IN) :: A
     REAL(realk)  :: mat_scalapack_sqnorm2
     mat_scalapack_sqnorm2 = mat_scalapack_DOTPRODUCT(A,A)
   end function mat_scalapack_sqnorm2

!> \brief See mat_trAB in mat-operations.f90
  function mat_scalapack_outdia_sqnorm2(A)
     implicit none
     TYPE(Matrix), intent(IN) :: A
     REAL(realk)  :: mat_scalapack_outdia_sqnorm2
     !
     INTEGER      :: DESC_A(9), i,j,li,lj,lpi,lpj
     REAL(realk)  :: TR2(1)
#ifdef VAR_SCALAPACK
     CALL PDM_SYNC(Job_outdia_sqnrm2,A)
     CALL PDM_DSCINIT(DESC_A,A)
     TR2=0.0E0_realk
     CALL PDLADIAG(A%nrow,A%p,1,1,DESC_A,PLUGIN_SQ2,TR2)
     CALL DGSUM2D( SLGrid%ictxt,'A','I',1,1,TR2,1,0,0)
     mat_scalapack_outdia_sqnorm2 = mat_scalapack_sqnorm2(A) - TR2(1)
#else
     mat_scalapack_outdia_sqnorm2 = 0
#endif
   end function mat_scalapack_outdia_sqnorm2

  !> \brief See mat_abs_max_elms in mat-operations.f90
  subroutine mat_scalapack_abs_max_elm(A,val)
     implicit none
     TYPE(Matrix), INTENT(IN)    :: A
     real(realk),intent(inout)   :: val
     !
     INTEGER      :: DESC_A(9), i,j
     REAL(REALK)  :: MAXA
#ifdef VAR_SCALAPACK

     CALL PDM_SYNC(Job_absmax,A)
     val=0.0E0_realk
     !max localy
     DO j=1,A%localncol
        DO i=1,A%localnrow
           val = MAX(val,ABS(A%p(i,j)))
        ENDDO
     ENDDO
     CALL lsmpi_max_realk_reduction(val,infpar%master,scalapack_comm)
#endif
   end subroutine mat_scalapack_abs_max_elm

  !> \brief See mat_max_elms in mat-operations.f90
  subroutine mat_scalapack_max_elm(A,val,pos)
     implicit none
     TYPE(Matrix), INTENT(IN)    :: A
     real(realk),intent(inout)     :: val
     integer, intent(inout)        :: pos(2)
     !
     real(realk), allocatable    :: vals(:)
     integer(kind=ls_mpik), allocatable        :: poss(:)
     INTEGER      :: DESC_A(9), i,j
     REAL(REALK)  :: MAXA
     !mpi stuff:
     integer(kind=ls_mpik) :: ierr, comm, one,two,zero,pos_mpi(2)
#ifdef VAR_SCALAPACK
     one=1
     two=2
     zero=0

     CALL PDM_SYNC(Job_max,A)

     val=-huge(0.0_realk)

     !max localy
     DO j=1,A%localncol
        DO i=1,A%localnrow
           if (A%p(i,j) > val) then
              val = A%p(i,j)
              pos_mpi(1)=i; pos_mpi(2)=j
           endif
        ENDDO
     ENDDO

     !convert local positions into global positions
     pos_mpi(1) = indxl2g(pos_mpi(1),BLOCK_SIZE,SLGrid%myrow,zero,SLGrid%nprow)
     pos_mpi(2) = indxl2g(pos_mpi(2),BLOCK_SIZE,SLGrid%mycol,zero,SLGrid%npcol)

     !allocate space
     allocate(vals(SLGrid%nprocs),poss(2*SLGrid%nprocs))
 
     !gather local max values and ther position
     call MPI_GATHER(val,one,MPI_DOUBLE_PRECISION,vals,one,MPI_DOUBLE_PRECISION,zero,SLGrid%comm,ierr)
     call MPI_GATHER(pos_mpi,two,MPI_INTEGER,poss,two,MPI_INTEGER,zero,SLGrid%comm,ierr)

     !find max element
     val=-huge(0.0_realk)
     do i=1,SLGrid%nprocs
       if (vals(i) > val) then
          val = vals(i)
          pos_mpi = poss((i*2) -1: i*2);
       endif
     enddo
     
     !free arrays
     deallocate(vals,poss)
     pos=pos_mpi

#endif
   end subroutine mat_scalapack_max_elm

  subroutine mat_scalapack_min_elm(A,val,pos)
     implicit none
     TYPE(Matrix), INTENT(IN)    :: A
     real(realk),intent(inout)     :: val
     integer, intent(inout)        :: pos(2)
     !
     real(realk), allocatable    :: vals(:)
     integer(kind=ls_mpik), allocatable :: poss(:)
     INTEGER      :: DESC_A(9), i,j
     REAL(REALK)  :: MAXA
     integer(kind=ls_mpik) :: ierr, comm, one,two,zero,pos_mpi(2)
#ifdef VAR_SCALAPACK
     one=1
     two=2
     zero=0

     CALL PDM_SYNC(Job_min,A)

     val=huge(0.0_realk)

     !max localy
     DO j=1,A%localncol
        DO i=1,A%localnrow
           if (A%p(i,j) < val) then
              val = A%p(i,j)
              pos_mpi(1)=i; pos_mpi(2)=j
           endif
        ENDDO
     ENDDO

     !convert local positions into global positions
     pos_mpi(1) = indxl2g(pos_mpi(1),BLOCK_SIZE,SLGrid%myrow,zero,SLGrid%nprow)
     pos_mpi(2) = indxl2g(pos_mpi(2),BLOCK_SIZE,SLGrid%mycol,zero,SLGrid%npcol)

     !allocate space
     allocate(vals(SLGrid%nprocs),poss(2*SLGrid%nprocs))
 
     !gather local max values and ther position
     call MPI_GATHER(val,one,MPI_DOUBLE_PRECISION,vals,one,MPI_DOUBLE_PRECISION,zero,SLGrid%comm,ierr)
     call MPI_GATHER(pos_mpi,two,MPI_INTEGER,poss,two,MPI_INTEGER,zero,SLGrid%comm,ierr)
     !find max element
     val=huge(0.0_realk)
     do i=1,SLGrid%nprocs
       if (vals(i) < val) then
          val = vals(i)
          pos_mpi = poss((i*2) -1: i*2);
       endif
     enddo
     
     !free arrays
     deallocate(vals,poss)
     pos=pos_mpi

#endif
   end subroutine mat_scalapack_min_elm


  !> \brief See mat_max_diag_elm in mat-operations.f90
  subroutine mat_scalapack_max_diag_elm(A,pos,val)
     implicit none
     TYPE(Matrix), INTENT(IN)    :: A
     real(realk),intent(inout)   :: val
     integer,intent(inout)       :: pos
     !
     INTEGER      :: DESC_A(9),i,j
     REAL(REALK)  :: MAXA(1)
#ifdef VAR_SCALAPACK

     CALL PDM_SYNC(Job_maxdiag,A)

     CALL PDM_DSCINIT(DESC_A,A)
     
     CALL PDLADIAG(A%nrow,A%p,1,1,DESC_A,PLUGIN_MAX,MAXA)
     val=MAXA(1)
     CALL lsmpi_max_realk_reduction(val,infpar%master,scalapack_comm)
     pos = 0
#endif
   end subroutine mat_scalapack_max_diag_elm

  !> \brief See mat_extract_diagonal in mat-operations.f90
  subroutine mat_scalapack_extract_diagonal(diag_a,A)
     implicit none
     type(Matrix), intent(in)   :: A
     real(realk), intent(inout) :: diag_a(A%nrow)
     INTEGER :: DESC_A(9)
     !
#ifdef VAR_SCALAPACK
     CALL PDM_SYNC(Job_extract_diag,A)
     CALL PDM_DSCINIT(DESC_A, A)
     diag_a = 0.0E0_realk     
     CALL PDLADIAG(A%nrow,A%p,1,1,DESC_A, PLUGIN_EXTRACT, diag_a)
     CALL DGSUM2D(SLGrid%ictxt, 'A', 'I', A%nrow, 1, diag_a, A%nrow, 0,0)
#endif
   end subroutine mat_scalapack_extract_diagonal

  !> \brief See mat_add_identity in mat-operations.f90
  !> A = alpha*I + beta*B
  subroutine mat_scalapack_add_identity(alpha, beta, B, A)
    implicit none
    TYPE(Matrix), intent(IN) :: B
    REAL(realk), INTENT(IN)  :: alpha, beta
    TYPE(Matrix)             :: A
    !
    INTEGER      :: DESC_A(DLEN_)
    REAL(realk)      :: AL(2)
#ifdef VAR_SCALAPACK

    IF(beta.EQ. 0E0_realk)then
       call mat_scalapack_zero(A)
    ELSE
       call mat_scalapack_copy(beta,B,A)
    ENDIF
    IF(alpha.NE. 0E0_realk)THEN
       CALL PDM_SYNC(Job_add_identity,A)
       CALL PDM_DSCINIT(DESC_A,A)
       AL(1) = ALPHA
       CALL DGEBS2D(SLGrid%ictxt,'A','I',2,1,AL,2)
       CALL PDLADIAG(A%nrow,A%p,1,1,DESC_A,PLUGIN_ADD,AL)
    ENDIF
#endif
  end subroutine mat_scalapack_add_identity

  !> \brief See mat_zero in mat-operations.f90
  subroutine mat_scalapack_zero(A)
    implicit none
    TYPE(Matrix) :: A
#ifdef VAR_SCALAPACK
    CALL PDM_SET(A,0.0E0_realk,0.0E0_realk)
#endif
  end subroutine mat_scalapack_zero

  !> \brief See mat_diag_f in mat-operations.f90
  subroutine mat_scalapack_diag_f(A,B,eival,C)
    !solves AC = eBC   (A=F, B=S)
    implicit none
    TYPE(Matrix), intent(IN)  :: A,B
    type(matrix)              :: C          !output
    real(realk),intent(INOUT) :: eival(:) !output
    !
    real(realk),pointer :: WORK(:),GAP(:)
    integer,pointer :: IWORK(:),IFAIL(:),ICLUSTR(:)
    real(realk) :: ABSTOL,ORFAC,DDUMMY
    integer :: LWORK,LIWORK,INFO,IDUMMY,DESC_A(9),DESC_B(9),DESC_C(9)
    integer :: neigenvalues,neigenvectors
#ifdef VAR_SCALAPACK
    REAL(REALK),EXTERNAL :: PDLAMCH

    CALL PDM_SYNC(Job_diag_f,A,B,C)

    CALL PDM_DSCINIT(DESC_A,A)
    CALL PDM_DSCINIT(DESC_B,B)
    CALL PDM_DSCINIT(DESC_C,C)

    call mat_scalapack_diag_f_aux(A,B,eival,C,DESC_A,DESC_B,DESC_C)
#endif
  end subroutine mat_scalapack_diag_f

  subroutine mat_scalapack_diag_f_aux(A,B,eival,C,DESC_A,DESC_B,DESC_C)
    !solves AC = eBC   (A=F, B=S)
    implicit none
    TYPE(Matrix), intent(IN)  :: A,B
    type(matrix)              :: C          !output
    real(realk),intent(INOUT) :: eival(:) !output
    integer,intent(IN) :: DESC_A(9),DESC_B(9),DESC_C(9)
    !
    real(realk),pointer :: WORK(:),GAP(:),Atmp(:,:),Btmp(:,:)
    integer,pointer :: IWORK(:),IFAIL(:),ICLUSTR(:)
    real(realk) :: ABSTOL,ORFAC,DDUMMY
    integer :: LWORK,LIWORK,INFO,IDUMMY,LWORK_SAVE
    integer :: neigenvalues,neigenvectors,LIWORK_SAVE
!    INTEGER :: NUMROC
#ifdef VAR_SCALAPACK
!    INTEGER            ICEIL
!    EXTERNAL           ICEIL
    REAL(REALK),EXTERNAL :: PDLAMCH
!    EXTERNAL NUMROC

    call mem_alloc(WORK,10)
    call mem_alloc(iWORK,10)
    call mem_alloc(iFAIL,A%nrow)
    call mem_alloc(iCLUSTR,2*SLgrid%nprow*SLgrid%npcol)
    call mem_alloc(GAP,SLgrid%nprow*SLgrid%npcol)
    DDUMMY=0E0_realk
    IDUMMY=0
    ORFAC=1.0E-3_realk
    ABSTOL = PDLAMCH(SLGrid%ictxt, 'U')
    LWORK = -1  !A workspace query is assumed
    LIWORK = -1 !A workspace query is assumed
    CALL PDSYGVX(1,'V','A','L',A%nrow,A%p,1,1,DESC_A,&
               & B%p,1,1,DESC_B,DDUMMY,DDUMMY,IDUMMY,IDUMMY,&
               & ABSTOL,neigenvalues,neigenvectors,eival,&
               & ORFAC,C%p,1,1,DESC_C,WORK,LWORK,IWORK,LIWORK,&
               & ifail,iclustr,gap,INFO)
    IF(INFO.NE. 0)THEN
       CALL LSQUIT('matop_scalapack mat_diag_f: PDSYGVX Failed A ',-1)
    ENDIF
    LWORK  = INT(WORK(1))
    LIWORK = IWORK(1)
    LWORK_SAVE  = INT(WORK(1))
    LIWORK_SAVE = IWORK(1)
    call mem_dealloc(WORK)
    call mem_dealloc(IWORK)
    call mem_alloc(WORK,LWORK)
    call mem_alloc(IWORK,LIWORK)

    call mem_alloc(Atmp,size(A%p,1),size(A%p,2))
    call mem_alloc(Btmp,size(B%p,1),size(B%p,2))
    Atmp = A%p
    Btmp = B%p

    CALL PDSYGVX(1,'V','A','L',A%nrow,Atmp,1,1,DESC_A,&
               & Btmp,1,1,DESC_B,DDUMMY,DDUMMY,IDUMMY,IDUMMY,&
               & ABSTOL,neigenvalues,neigenvectors,eival,&
               & ORFAC,C%p,1,1,DESC_C,WORK,LWORK,IWORK,LIWORK,&
               & ifail,iclustr,gap,INFO)
    call mem_dealloc(WORK)
    call mem_dealloc(IWORK)
    
    IF(INFO.NE. 0)THEN
       print*,'PDSYGVX Failed: Increasing WorkArray INFO',INFO
       INFO = 0 
       !For PDSYGVX, the computed eigenvectors may not be orthogonal 
       !if the minimum workspace is supplied and ortol is too small; 
       !therefore, if you want to guarantee orthogonality 
       !(at the cost of potentially compromising performance), 
       !you should add the following to lwork: (clustersize-1)(n)
       !where clustersize is the number of eigenvalues in the largest cluster
       !, where a cluster is defined as a set of close eigenvalues: 
       !       LWORK  = LWORK + A%nrow*A%nrow
       Atmp = A%p
       Btmp = B%p

       LWORK  = LWORK_SAVE + A%nrow*A%nrow
       LIWORK = LIWORK_SAVE
       call mem_alloc(WORK,LWORK)
       call mem_alloc(IWORK,LIWORK)
       CALL PDSYGVX(1,'V','A','L',A%nrow,Atmp,1,1,DESC_A,&
            & Btmp,1,1,DESC_B,DDUMMY,DDUMMY,IDUMMY,IDUMMY,&
            & ABSTOL,neigenvalues,neigenvectors,eival,&
            & ORFAC,C%p,1,1,DESC_C,WORK,LWORK,IWORK,LIWORK,&
            & ifail,iclustr,gap,INFO)
       call mem_dealloc(WORK)
       call mem_dealloc(IWORK)
    ENDIF
    IF(INFO.NE. 0)THEN
       print*,'PDSYGVX Failed: INFO',INFO
       IF(INFO.GT. 0)THEN
          IF(MOD(INFO,2).NE. 0)THEN
             print*,'one or more eigenvectors failed to converge.'
             print*,'Their indices are stored in IFAIL.'
          ELSEIF(MOD(INFO/2,2).NE. 0)THEN
             print*,'eigenvectors corresponding to one or more clusters of eigenvalues'
             print*,'could not be reorthogonalized because of insufficient workspace.'
             print*,'The indices of the clusters are stored in the array ICLUSTR.'
          ELSEIF(MOD(INFO/4,2).NE. 0)then
             print*,'space limit prevented PDSYGVX from computing all of' 
             print*,'the eigenvectors between VL and VU. The number of eigenvectors'
             print*,'computed is returned in NZ.'
          ELSEIF(MOD(INFO/8,2).NE. 0)Then 
             print*,'PDSTEBZ failed to compute eigenvalues.'
          ELSEIF(MOD(INFO/16,2).NE. 0)then 
             print*,'B was not positive definite.'  
             print*,'IFAIL(1) indicates the order of'
             print*,'the smallest minor which is not positive definite.'
          ELSE
             print*,'Unknown reason'
          ENDIF
       ENDIF
       CALL LSQUIT('matop_scalapack mat_diag_f: PDSYGVX Failed B2 ',-1)
    ENDIF
    call mem_dealloc(Atmp)
    call mem_dealloc(Btmp)
    call mem_dealloc(IFAIL)
    call mem_dealloc(ICLUSTR)
    call mem_dealloc(GAP)
#endif
  end subroutine mat_scalapack_diag_f_aux

  !> \brief See mat_dsyev in mat-operations.f90
  subroutine mat_scalapack_dsyev(A,B,eival,ndim)
    !B will be eigenvectors and eival the eigenvalues  
    implicit none
    TYPE(Matrix), intent(IN)  :: A
    TYPE(Matrix), intent(INOUT)  :: B
    real(realk),intent(INOUT) :: eival(ndim) !output
    integer,intent(in) :: ndim
    !
    real(realk),pointer :: WORK(:),GAP(:)
    integer,pointer :: IWORK(:),IFAIL(:),ICLUSTR(:)
    real(realk) :: ABSTOL,ORFAC,DDUMMY
    integer :: LWORK,LIWORK,INFO,IDUMMY,DESC_A(9),DESC_B(9),DESC_C(9)
    integer :: neigenvalues,neigenvectors
#ifdef VAR_SCALAPACK
    REAL(REALK),EXTERNAL :: PDLAMCH

    CALL PDM_SYNC(Job_dsyev,A,B)
    CALL PDM_DSCINIT(DESC_A,A)
    CALL PDM_DSCINIT(DESC_B,B)
    call mat_scalapack_dsyev_aux(A,B,eival,ndim,DESC_A,DESC_B)
#endif
  end subroutine mat_scalapack_dsyev

  subroutine mat_scalapack_dsyev_aux(A,B,eival,ndim,DESC_A,DESC_B)
    !B will be eigenvectors and eival the eigenvalues  
    implicit none
    TYPE(Matrix), intent(IN)     :: A
    TYPE(Matrix), intent(INOUT)  :: B
    real(realk),intent(INOUT) :: eival(ndim) !output
    integer,intent(IN) :: DESC_A(9),DESC_B(9),ndim
    !
    real(realk),pointer :: WORK(:),Atmp(:,:)
    integer :: LWORK,INFO
#ifdef VAR_SCALAPACK
    REAL(REALK),EXTERNAL :: PDLAMCH

    call mem_alloc(WORK,10)
    LWORK = -1  !A workspace query is assumed
    CALL PDSYEV('V','U',A%nrow,A%p,1,1,DESC_A,&
               & eival,B%p,1,1,DESC_B,WORK,LWORK,INFO)
    IF(INFO.NE. 0)THEN
       CALL LSQUIT('matop_scalapack mat_dsyev: PDSYEV Failed A ',-1)
    ENDIF
    LWORK = INT(WORK(1))
    call mem_dealloc(WORK)
    call mem_alloc(WORK,LWORK)
    call mem_alloc(Atmp,size(A%p,1),size(A%p,2))
    Atmp = A%p
    CALL PDSYEV('V','U',A%nrow,Atmp,1,1,DESC_A,&
               & eival,B%p,1,1,DESC_B,WORK,LWORK,INFO)
    call mem_dealloc(Atmp)
    IF(INFO.NE. 0)THEN
       print*,'INFO',INFO
       CALL LSQUIT('matop_scalapack mat_dsyev: PDSYEV Failed B',-1)
    ENDIF
    call mem_dealloc(WORK)
#endif
  end subroutine mat_scalapack_dsyev_aux

  !> \brief See mat_dsyevx in mat-operations.f90
  subroutine mat_scalapack_dsyevx(A,eival,ieig)
    !B will be eigenvectors and eival the eigenvalues  
    implicit none
    TYPE(Matrix), intent(INOUT)  :: A
    real(realk),intent(INOUT)    :: eival
    integer,intent(in)           :: ieig
    !
#ifdef VAR_SCALAPACK
    integer :: ndim
    integer :: DESC_A(9)

    CALL PDM_SYNC(Job_dsyevx,A)
    CALL PDM_DSCINIT(DESC_A,A)
    CALL ls_mpibcast(ieig,infpar%master,scalapack_comm)
    call mat_scalapack_dsyevx_aux(A,eival,ieig,DESC_A)
#endif
  end subroutine mat_scalapack_dsyevx

  subroutine mat_scalapack_dsyevx_aux(A,eival,ieig,DESC_A)
    !B will be eigenvectors and eival the eigenvalues  
    implicit none
    TYPE(Matrix), intent(INOUT) :: A
    real(realk),intent(INOUT)   :: eival
    integer,intent(IN) :: DESC_A(9),ieig
    !
#ifdef VAR_SCALAPACK
    real(realk),pointer :: choltemp(:),eivec(:)
    integer,pointer     :: icholtemp(:),IFAIL(:),ICLUSTR(:)
    integer             :: ndim,neig,info,VL,VU
    real(realk),pointer :: eivalTmp(:),GAP(:),Atmp(:,:)
    integer :: lwork,liwork


    ndim = A%nrow
    neig = 1
    VL = 0
    VU = 0
    lwork = -1
    liwork = -1
    call mem_alloc(choltemp,10)
    call mem_alloc(icholtemp,10)
    call mem_alloc(ifail,ndim)
    call mem_alloc(eivalTmp,ndim)
    call mem_alloc(eivec,ndim)
    call mem_alloc(GAP,ndim)
    call mem_alloc(ICLUSTR,ndim)
!   Inquire to get the size of lwork and liwork
    call mem_alloc(Atmp,size(A%p,1),size(A%p,2))
    Atmp = A%p
    call PDSYEVX('N', 'I', 'U', ndim, Atmp, 1, 1, DESC_A, VL, VU, ieig, ieig, &
       &  0.0E0_realk, neig, 0, eivalTmp, 0.0E0_realk, eivec, 1, 1, DESC_A, &
       &  choltemp, lwork, icholtemp, liwork, &
       &  IFAIL, ICLUSTR, GAP, INFO )
    call mem_dealloc(Atmp)
    lwork =  INT(choltemp(1))
    liwork = icholtemp(1)
    call mem_dealloc(choltemp)
    call mem_dealloc(icholtemp)
    call mem_alloc(choltemp,lwork)
    call mem_alloc(icholtemp,liwork)
!   Actual calculation of the given eigenvalue
    call PDSYEVX('N', 'I', 'U', ndim, A%p, 1, 1, DESC_A, VL, VU, ieig, ieig, &
       &  0.0E0_realk, neig, 0, eivalTmp, 0.0E0_realk, eivec, 1, 1, DESC_A, &
       &  choltemp, lwork, icholtemp, liwork, &
       &  IFAIL, ICLUSTR, GAP, INFO )
    eival = eivalTmp(1)
    if(info.ne. 0) then
       print*,'mat_scalapack_dsyevx_aux: pdsyevx failed, info=',info
       call lsquit('mat_scalapack_dsyevx_aux: diagonalization failed.',-1)
    end if
    call mem_dealloc(choltemp)
    call mem_dealloc(icholtemp)
    call mem_dealloc(ifail)
    call mem_dealloc(eivalTmp)
    call mem_dealloc(eivec)
    call mem_dealloc(GAP)
    call mem_dealloc(ICLUSTR)
#endif
  end subroutine mat_scalapack_dsyevx_aux

  !> \brief See mat_dpotrf in mat-operations.f90
  subroutine mat_scalapack_dpotrf(A)
    !B will be eigenvectors and eival the eigenvalues  
    implicit none
    TYPE(Matrix), intent(INOUT) :: A
    !
    integer :: DESC_A(9)
#ifdef VAR_SCALAPACK

    CALL PDM_SYNC(Job_dpotrf,A)
    CALL PDM_DSCINIT(DESC_A,A)
    call mat_scalapack_dpotrf_aux(A,DESC_A)
#endif
  end subroutine mat_scalapack_dpotrf

  subroutine mat_scalapack_dpotrf_aux(A,DESC_A)
    !B will be eigenvectors and eival the eigenvalues  
    implicit none
    TYPE(Matrix), intent(INOUT) :: A
    integer,intent(IN)          :: DESC_A(9)
    !
    integer :: INFO
#ifdef VAR_SCALAPACK

!   Finds the Cholesky factorization of A 
    CALL PDPOTRF('U',A%nrow,A%p,1,1,DESC_A,INFO)
    IF(INFO.NE. 0)THEN
       CALL LSQUIT('matop_scalapack mat_dpotrf_aux: PDPOTRF Failed',-1)
    ENDIF
#endif
  end subroutine mat_scalapack_dpotrf_aux

  !> \brief See mat_dpotrf in mat-operations.f90
  subroutine mat_scalapack_dpotrs(A,B)
    !B will be eigenvectors and eival the eigenvalues  
    implicit none
    TYPE(Matrix), intent(IN)    :: A
    TYPE(Matrix), intent(INOUT) :: B
    !
    integer :: DESC_A(9),DESC_B(9)
#ifdef VAR_SCALAPACK

    CALL PDM_SYNC(Job_dpotrs,A,B)
    CALL PDM_DSCINIT(DESC_A,A)
    CALL PDM_DSCINIT(DESC_B,B)
    call mat_scalapack_dpotrs_aux(A,B,DESC_A,DESC_B)
#endif
  end subroutine mat_scalapack_dpotrs

  subroutine mat_scalapack_dpotrs_aux(A,B,DESC_A,DESC_B)
    !B will be eigenvectors and eival the eigenvalues  
    implicit none
    TYPE(Matrix), intent(IN)    :: A
    TYPE(Matrix), intent(INOUT) :: B
    integer,intent(IN)          :: DESC_A(9),DESC_B(9)
    !
    integer :: INFO
#ifdef VAR_SCALAPACK

!   Finds the Cholesky factorization of A 
    CALL PDPOTRS('U',A%nrow,B%ncol,A%p,1,1,DESC_A,B%p,1,1,DESC_B,INFO)
    IF(INFO.NE. 0)THEN
       CALL LSQUIT('matop_scalapack mat_dpotrs_aux: PDPOTRS Failed',-1)
    ENDIF
#endif
  end subroutine mat_scalapack_dpotrs_aux

  !> \brief See mat_ao_precond in mat-operations-aux.f90
!> \param symm Symmetry indicator: symmetric = 1, antisymmetric = 2, nonsymmetric = 0
!> \param omega Level shift
!> \param P Fock/KS matrix in OAO basis, occupied part (virtual part projected out)
!> \param Q Fock/KS matrix in OAO basis, virtual part (occupied part projected out)
!> \param U Density matrix in OAO basis
!> \param X_AO Matrix to be preconditioned
!> 
!> Preconditioning with orbital energy difference. This is done by using the occupied and virtual 
!> parts of the Fock/KS matrix: 
!> X_prec(i,j) = X(i,j) / [Q(j,j)-P(j,j) + Q(i,i)-P(i,i) - omega*(U(j,j)-U(i,i))]
!> taking care not to divide by zero and exploiting symmetry if X is symm or antisymm
  subroutine mat_scalapack_ao_precond(symm,omega,P,Q,U,A)
    implicit none
    integer, intent(in)         :: symm
    real(realk), intent(in)     :: omega
    type(Matrix), intent(in)    :: P,Q,U
    type(Matrix), intent(inout) :: A
    !
#ifdef VAR_SCALAPACK
    integer :: DESC_A(9)
    Real(REALK),pointer :: diagQP(:),diagU(:)
    integer :: i,isymm
    real(realk) :: romega
    call mem_alloc(diagU,A%nrow)
    call mem_alloc(diagQP,A%nrow)

    !Master only
    call mat_scalapack_extract_diagonal(diagU,P)
    call mat_scalapack_extract_diagonal(diagQP,Q)
    DO i=1,A%nrow
      diagQP(i) = diagQP(i) - diagU(i)
    ENDDO
    call mat_scalapack_extract_diagonal(diagU,U)

    isymm  = symm   !Due to intent(IN) of symm
    romega = omega  !Due to intent(IN) of omega
    !END Master only

    CALL PDM_SYNC(Job_ao_precond,A)
    CALL PDM_DSCINIT(DESC_A,A)
    CALL mat_scalapack_ao_precond_aux(A,A%nrow,DESC_A,isymm,romega,diagQP,diagU)

    call mem_dealloc(diagQP)
    call mem_dealloc(diagU)
#endif
  end subroutine mat_scalapack_ao_precond

 
  SUBROUTINE mat_scalapack_ao_precond_aux(A,nrow,DESC_A,symm,omega,diagQP,diagU)
  implicit none
  integer, intent(inout)      :: symm
  real(realk), intent(inout)  :: omega
  type(Matrix), intent(inout) :: A
  integer, intent(IN)         :: DESC_A(9)
  integer,intent(IN)          :: nrow
  Real(REALK),intent(INOUT)   :: diagQP(nrow),diagU(nrow)
  !
#ifdef VAR_SCALAPACK
    integer :: i
    Real(REALK),pointer :: vrow(:),vcol(:),wrow(:),wcol(:)
    integer :: irow,icol,prow,pcol,grow,gcol
    real(realk) :: denom

    call ls_mpibcast(symm,infpar%master,scalapack_comm)
    call ls_mpibcast(omega,infpar%master,scalapack_comm)
    call ls_mpibcast(diagQP,nrow,infpar%master,scalapack_comm)
    call ls_mpibcast(diagU,nrow,infpar%master,scalapack_comm)

    call mem_alloc(vrow,A%localnrow)
    call mem_alloc(wrow,A%localnrow)
    call mem_alloc(vcol,A%localncol)
    call mem_alloc(wcol,A%localncol)
    DO irow=1,A%localnrow
      grow = INDXL2G( irow, BLOCK_SIZE, SLgrid%myrow, 0, SLgrid%nprow ) !find global row index
      vrow(irow) = diagQP(grow)
      wrow(irow) = diagU(grow)
    ENDDO
    DO icol=1,A%localncol
      gcol = INDXL2G( icol, BLOCK_SIZE,SLgrid%mycol , 0, SLgrid%npcol ) !find global col index
      vcol(icol) = diagQP(gcol)
      wcol(icol) = diagU(gcol)
    ENDDO
!   Symmetric or Anti-symmetric case
    IF ((symm == 1) .OR. (symm == 2)) THEN
      DO icol=1,A%localncol
        DO irow=1,A%localnrow
          denom = vcol(icol)+vrow(irow) - omega
          if (ABS(denom) > 1.0E-10_realk) then
              A%p(irow,icol) = A%p(irow,icol)/(denom)
           endif
        ENDDO
      ENDDO
!   Non-symmetric case (I do not see why this case is different, but I here copy the mat_dense_ao_precond SR)
    ELSE
      DO icol=1,A%localncol
        DO irow=1,A%localnrow
          denom = vrow(irow)+vcol(icol) - omega*(wcol(icol)-wrow(irow))
          if (abs(denom) < 1.0E-1_realk) then !Do not divide by too small elements
             denom = denom*1.0E-1_realk/(abs(denom)) !Keep the sign on denominator
          endif
          A%p(irow,icol) = A%p(irow,icol)/(denom)
        ENDDO
      ENDDO
    ENDIF

    call mem_dealloc(vrow)
    call mem_dealloc(wrow)
    call mem_dealloc(vcol)
    call mem_dealloc(wcol)
#endif
  END SUBROUTINE mat_scalapack_ao_precond_aux

#ifdef VAR_SCALAPACK  
      INTEGER FUNCTION INDXL2G( INDXLOC, NB, IPROC, ISRCPROC, NPROCS )
!
!  -- ScaLAPACK tools routine (version 1.7) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 1, 1997
!
!     .. Scalar Arguments ..
      INTEGER            INDXLOC, IPROC, ISRCPROC, NB, NPROCS
!     ..
!
!  Purpose
!  =======
!
!  INDXL2G computes the global index of a distributed matrix entry
!  pointed to by the local index INDXLOC of the process indicated by
!  IPROC.
!
!  Arguments
!  =========
!
!  INDXLOC   (global input) INTEGER
!            The local index of the distributed matrix entry.
!
!  NB        (global input) INTEGER
!            Block size, size of the blocks the distributed matrix is
!            split into.
!
!  IPROC     (local input) INTEGER
!            The coordinate of the process whose local array row or
!            column is to be determined.
!
!  ISRCPROC  (global input) INTEGER
!            The coordinate of the process that possesses the first
!            row/column of the distributed matrix.
!
!  NPROCS    (global input) INTEGER
!            The total number processes over which the distributed
!            matrix is distributed.
!
!  =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC            MOD
!     ..
!     .. Executable Statements ..
!
      INDXL2G = NPROCS*NB*((INDXLOC-1)/NB) + MOD(INDXLOC-1,NB) + &
     &          MOD(NPROCS+IPROC-ISRCPROC, NPROCS)*NB + 1
!
      RETURN
!
!     End of INDXL2G
!
      END FUNCTION INDXL2G

      INTEGER FUNCTION INDXG2P( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS )
!
!  -- ScaLAPACK tools routine (version 1.7) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 1, 1997
!
!     .. Scalar Arguments ..
      INTEGER            INDXGLOB, IPROC, ISRCPROC, NB, NPROCS
!     ..
!
!  Purpose
!  =======
!
!  INDXG2P computes the process coordinate which posseses the entry of a
!  distributed matrix specified by a global index INDXGLOB.
!
!  Arguments
!  =========
!
!  INDXGLOB  (global input) INTEGER
!            The global index of the element.
!
!  NB        (global input) INTEGER
!            Block size, size of the blocks the distributed matrix is
!            split into.
!
!  IPROC     (local dummy) INTEGER
!            Dummy argument in this case in order to unify the calling
!            sequence of the tool-routines.
!
!  ISRCPROC  (global input) INTEGER
!            The coordinate of the process that possesses the first
!            row/column of the distributed matrix.
!
!  NPROCS    (global input) INTEGER
!            The total number processes over which the matrix is
!            distributed.
!
!  =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          MOD
!     ..
!     .. Executable Statements ..
!
      INDXG2P = MOD( ISRCPROC + (INDXGLOB - 1) / NB, NPROCS )
!
      RETURN
!
!     End of INDXG2P
!
      END FUNCTION INDXG2P

      INTEGER FUNCTION INDXG2L( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS )
!
!  -- ScaLAPACK tools routine (version 1.7) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 1, 1997
!
!     .. Scalar Arguments ..
      INTEGER            INDXGLOB, IPROC, ISRCPROC, NB, NPROCS
!     ..
!
!  Purpose
!  =======
!
!  INDXG2L computes the local index of a distributed matrix entry
!  pointed to by the global index INDXGLOB.
!
!  Arguments
!  =========
!
!  INDXGLOB  (global input) INTEGER
!            The global index of the distributed matrix entry.
!
!  NB        (global input) INTEGER
!            Block size, size of the blocks the distributed matrix is
!            split into.
!
!  IPROC     (local dummy) INTEGER
!            Dummy argument in this case in order to unify the calling
!            sequence of the tool-routines.
!
!  ISRCPROC  (local dummy) INTEGER
!            Dummy argument in this case in order to unify the calling
!            sequence of the tool-routines.
!
!  NPROCS    (global input) INTEGER
!            The total number processes over which the distributed
!            matrix is distributed.
!
!  =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          MOD
!     ..
!     .. Executable Statements ..
!
      INDXG2L = NB*((INDXGLOB-1)/(NB*NPROCS))+MOD(INDXGLOB-1,NB)+1
!
      RETURN
!
!     End of INDXG2L
!
    END FUNCTION INDXG2L
#endif
 end module matrix_operations_scalapack

! SUbroutines that must be places OUTSIDE the module and calls the module.
! These are MPI routines that activate the slaves. 

 SUBROUTINE PDM_GRIDINIT(PR,PC,NBAST)
   use matrix_operations_scalapack
   use memory_handling
   use matrix_module
#ifdef VAR_MPI
   use infpar_module
   use lsmpi_type
#endif
   IMPLICIT NONE
   INTEGER :: PR, PC, TMP(4), IERR, I, NBAST
#ifdef VAR_SCALAPACK
   integer, external :: blacs2sys_handle
   call ls_mpibcast(GRIDINIT,infpar%master,MPI_COMM_LSDALTON)   
   TMP(1)=PR; TMP(2)=PC; TMP(3)=NBAST; TMP(4)=infpar%inputblocksize
  
  
   ENTRY PDM_GRIDINIT_SLAVE
   SLGrid%mynum = scalapack_mynum
   SLGrid%nprocs = infpar%nodtot   
   call ls_mpibcast(TMP,4,infpar%master,scalapack_comm)
   !SLGrid%ICTXT specifies the BLACS context handle identifying the

   !created process grid.
   ! In MPIBLACS the BLACS system contexts are MPI communicators so 
   ! we can assign them directly.
   ! SL_INIT (more specifically blacs_get) should not be called
   ! as it can only return a global context (i.e. world communicator)
   ! The calculation will therefore crash in the call to SL_INIT 
   ! if scalapack_comm is not equal to MPI_COMM_WORLD.
   ! In the case of CHEM_SHELL, scalapack_comm is not equal to
   ! MPI_COMM_WORLD and SL_INIT should not be called directly
   SLGrid%ictxt = scalapack_comm
   ! Set up the process grid - sets ictxt to a BLACS handle
   CALL BLACS_GRIDINIT(SLGrid%ictxt, 'R', TMP(1), TMP(2))
   ! Again, do not call blacs_get. The BLACS system context
   ! and corresponding communicator are simply MPI communicators.
   SLGrid%Masterictxt = scalapack_comm
   SLGrid%comm = scalapack_comm

   IF(TMP(4).EQ.0)THEN
      IF(TMP(3).GT.200)THEN  !nbast > 200
         !set 80 to be minimum block size 
         BLOCK_SIZE = MAX(80,TMP(3)/(MAX(TMP(1),TMP(2))))
         !     print*,'BLOCK_SIZE=',BLOCK_SIZE
         !     print*,'NBAST     =',TMP(3)
      ELSE
         !for debugging purposes it is nice to not have a minumum 
         !block size
         BLOCK_SIZE = TMP(3)/(MAX(TMP(1),TMP(2)))
      ENDIF
   ELSE
      BLOCK_SIZE = TMP(4)
!     print*,'TMP4 BLOCK_SIZE=',BLOCK_SIZE
!     print*,'NBAST     =',TMP(3)
   ENDIF
   !OUTPUT from BLACS_GRIDINFO:
   !SLGrid%nprow, SLGrid%npcol, SLGrid%myrow, SLGrid%mycol
   CALL BLACS_GRIDINFO(SLGrid%ictxt ,SLGrid%nprow,&
        & SLGrid%npcol, SLGrid%myrow, SLGrid%mycol )
   CALL DARRAY_NULLIFY

   CALL GRIDINFO_SETUP(TMP(3))
#endif      
 END SUBROUTINE PDM_GRIDINIT

 SUBROUTINE PDM_SLAVE
   use matrix_operations_scalapack
   use memory_handling
   use matrix_module
#ifdef VAR_MPI
   use infpar_module
   use lsmpi_type
#endif
   IMPLICIT NONE
   TYPE(Matrix) :: A, B, C, AUX
   CHARACTER :: T(2)
   INTEGER JOB, i, J
   INTEGER DESC_A(DLEN_), DESC_B(DLEN_), DESC_C(DLEN_), DESC_AF(DLEN_)
   real(REALK) :: AF, AB(2)
   real(REALK),allocatable :: AF2(:)
   real(REALK),allocatable :: diag(:),diag2(:)
   real(REALK),pointer :: diag3(:)
   logical :: logi
   INTEGER :: n1,n2,TMP(4),nocc,symm,mynum,nbuffer,bufferproc
   real(REALK) :: PDLATRA, DDOT,omega
#ifdef VAR_SCALAPACK
   integer :: localnrow(0:infpar%nodtot-1),loc,glb,l,nelemst,max_block_idx
   integer :: localncol(0:infpar%nodtot-1)
   integer :: bufferOffset(0:infpar%nodtot-1)
   integer :: fullrow,fullcol,insertrow,insertcol,nsize,fulldim1,fulldim2,ieig
   logical :: Communicationnodes(0:infpar%nodtot-1),debug
   real(realk),pointer :: Afullblock(:,:)   
   real(realk) :: eival
   type(lsmatrix) :: Abuffer(0:0)
   integer,pointer :: address_on_grid(:,:)
   integer :: m,n,k
   integer(kind=long) :: nsize1,nsize2
   integer(kind=ls_mpik) :: tmpi(4),ierr,one=1,two=2,zero=0
   INTEGER :: LLD, INFO
   integer, external :: numroc
   EXTERNAL PDLATRA, DDOT

   CALL PDM_SYNC(JOB,A,B,C) !Job is output

!   CALL PDM_DSCINIT(DESC_A,A)
!   CALL PDM_DSCINIT(DESC_B,B)
!   CALL PDM_DSCINIT(DESC_C,C)
   
   SELECT CASE(JOB)
   CASE(Job_init)
      CALL PDM_DSCINIT(DESC_A,A)
      I = ALLOC_IN_DARRAY(A)
#if 0
      !Takes a general rectangular matrix and sends it to the destination process
      CALL IGESD2D(SLGrid%ictxt,1,1,I,1,0,0)
#else
      call mem_alloc(address_on_grid,SLGrid%nprow,SLGrid%npcol)
      call ls_izero(address_on_grid,SLGrid%nprow*SLGrid%npcol)
      address_on_grid(SLGrid%myrow+1,SLGrid%mycol+1) = I
      CALL lsmpi_reduction(address_on_grid,SLGrid%nprow,SLGrid%npcol,infpar%master,scalapack_comm)
      call mem_dealloc(address_on_grid)
#endif
      nsize1=A%localnrow*A%localncol*mem_realsize
      nsize2=A%nrow*A%ncol*mem_realsize
      call mem_allocated_mem_type_matrix(nsize1,nsize2)
   CASE(Job_rand)
      CALL PDM_DSCINIT(DESC_A,A)
      do i=1,A%localnrow
         do j=1,A%localncol
            A%p(i,j)=scalapack_mynum+1
         enddo
      enddo
   CASE(Job_from_full)
      CALL PDM_DSCINIT(DESC_AF,A,A%nrow,A%ncol)
!      DESC_AF(1:9)=0
!      DESC_AF(2)=-1
      CALL PDM_DSCINIT(DESC_A,A)
      ! Slave Part: Distributes a full matrix AFull (on master) to SCALAPACK distributed matrix A
      IF(.NOT.infpar%ScalapackWorkAround)THEN
         CALL PDGEMR2D(A%nrow,A%ncol,AF,1,1,DESC_AF,&
              &A%p,1,1,DESC_A,SLGrid%ictxt)
      ELSE
         IF(A%localnrow*A%localncol.GT. 0)THEN
            call ls_mpisendrecv(A%p,A%localnrow,A%localncol,scalapack_comm,infpar%master,scalapack_mynum)
         ENDIF
      ENDIF
   CASE(Job_to_full3D)
      CALL PDM_DSCINIT(DESC_AF,A,A%nrow,A%ncol)
      CALL PDM_DSCINIT(DESC_A,A)
      ! Slave Part: Collects full matrix AFull (on master) from SCALAPACK distributed matrix A
      CALL PDGEMR2D(A%nrow,A%ncol,A%p,1,1,DESC_A,&
           &AF,1,1,DESC_AF,SLGrid%ictxt)
   CASE(Job_to_full)
      CALL PDM_DSCINIT(DESC_AF,A,A%nrow,A%ncol)
      CALL PDM_DSCINIT(DESC_A,A)
      ! Slave Part: Collects full matrix AFull (on master) from SCALAPACK distributed matrix A
      IF(.NOT.infpar%ScalapackWorkAround)THEN
         CALL PDGEMR2D(A%nrow,A%ncol,A%p,1,1,DESC_A,&
              &AF,1,1,DESC_AF,SLGrid%ictxt)
      ELSE
         IF(A%localnrow*A%localncol.GT. 0)THEN
            call ls_mpisendrecv(A%p,A%localnrow,A%localncol,scalapack_comm,scalapack_mynum,infpar%master)
         ENDIF
      ENDIF
   CASE(Job_scal)
      CALL DGEBR2D(SLGrid%ictxt,'A','I',1,1,AF,1,&
           &SLGrid%myrow,SLGrid%mycol)
      if (A%localnrow*A%localncol.NE. 0)then
         CALL DSCAL(A%localnrow*A%localncol,AF, A%p, 1)
      endif
   CASE(Job_assign)
      CALL PDM_DSCINIT(DESC_A,A)
      CALL PDM_DSCINIT(DESC_B,B)
      CALL PDLACPY('All',A%nrow,A%ncol,A%p,1,1,DESC_A,&
                     &B%p,1,1,DESC_B)
   CASE(Job_trans)
      CALL PDM_DSCINIT(DESC_A,A)
      CALL PDM_DSCINIT(DESC_B,B)
      CALL PDTRAN(A%ncol,A%nrow,1.0E0_realk,A%p,1,1,DESC_A,&
                  &0.0E0_realk,B%p,1,1,DESC_B)
   CASE(Job_trace)
      CALL PDM_DSCINIT(DESC_A,A)
      AF = PDLATRA(A%nrow,A%p,1,1,DESC_A)
   CASE(Job_dotproduct)
      AF = 0.0E0_realk
      IF (A%localnrow*A%localncol.NE. 0)then
         AF = DDOT(A%localnrow*A%localncol, A%p,1, B%p,1)
      endif
      CALL DGSUM2D( SLGrid%ictxt,'A','I',1,1,AF,1,0,0)
   CASE(Job_daxpy)
      CALL PDM_DSCINIT(DESC_A,A)
      CALL PDM_DSCINIT(DESC_B,B)
      CALL DGEBR2D(SLGrid%ictxt,'A','I',2,1,AB,2,&
           &       SLGrid%myrow,SLGrid%mycol)
      CALL PDMATADD(A%nrow,A%ncol,AB(1),A%p,1,1,DESC_A,&
           &        AB(2), B%p,1,1,DESC_B )
   CASE(Job_mul)
      CALL PDM_DSCINIT(DESC_A,A)
      CALL PDM_DSCINIT(DESC_B,B)
      CALL PDM_DSCINIT(DESC_C,C)
      CALL DGEBR2D(SLGrid%ictxt,'A','I',2,1,AB,2,&
           &SLGrid%myrow,SLGrid%mycol)
      call ls_mpibcast(T,2,infpar%master,scalapack_comm)
     if (T(1) == 'n' .or. T(1) == 'N') then
       m = a%nrow
       k = a%ncol
     elseif (T(1) == 't' .or. T(1) == 'T') then
       m = a%ncol
       k = a%nrow
     endif

     if (T(2) == 'n' .or. T(2) == 'N') then
       n = b%ncol
     elseif (T(2) == 't' .or. T(2) == 'T') then
       n = b%nrow
     endif

      CALL PDGEMM(T(1),T(2),m,n,k,AB(1),A%p,1,1,DESC_A,&
           &B%p,1,1,DESC_B,AB(2),C%p,1,1,DESC_C)
   CASE(Job_dger)
call lsquit('TKTKTKTKTKTKTKTKTKTK',-1)
! This is wrong diag is not allocated!!!!!
!      call PDM_DSCINIT(DESC_A,A)
!      CALL DGEBR2D(SLGrid%ictxt,'A','I',2,1,AB,2,&
!           &SLGrid%myrow,SLGrid%mycol)
!      CALL IGEBR2D(SLGrid%ictxt,'A','I',9,1,DESC_B,9,&
!           &SLGrid%myrow,SLGrid%mycol)
!      CALL IGEBR2D(SLGrid%ictxt,'A','I',9,1,DESC_C,9,&
!           &SLGrid%myrow,SLGrid%mycol)
!      !allocate(diag(A%nrow))
!      !allocate(diag2(A%ncol))
!      call pdger(A%nrow,A%ncol,AB(1),diag,1,1,DESC_B&
!            ,1,diag2,1,1,DESC_C,1,A%p,1,1,DESC_A,ierr)
!      !deallocate(diag,diag2)
   CASE(Job_hmul)
      CALL DGEBR2D(SLGrid%ictxt,'A','I',2,1,AB,2,&
           &SLGrid%myrow,SLGrid%mycol)
      C%p=AB(1)*A%p*B%p+AB(2)*C%p

   CASE(Job_hdiv)
      CALL DGEBR2D(SLGrid%ictxt,'A','I',2,1,AB,2,&
           &SLGrid%myrow,SLGrid%mycol)

      do j=1,A%localncol
       do i=1,A%localnrow
       AB(2) = B%p(i,j) - AB(1)
       if (AB(2).ne.0.0_realk) A%p(i,j) = A%p(i,j)/AB(2)
       enddo
      enddo


   CASE(Job_dmul)
      CALL DGEBR2D(SLGrid%ictxt,'A','I',2,1,AB,2,&
           &SLGrid%myrow,SLGrid%mycol)
      CALL SGEBR2D(SLGrid%ictxt,'A','I',2,1,T,2,&
           &SLGrid%myrow,SLGrid%mycol)
      allocate(diag(B%nrow))
      !CALL DGEBR2D(SLGrid%ictxt,'A','I',B%nrow,1,diag,B%nrow,&
      !   &SLGrid%myrow,SLGrid%mycol)
      call ls_mpibcast(diag,B%nrow,infpar%master,SLGrid%comm)
      if (B%localnrow*B%localncol.NE. 0 .and. (AB(2) /= 1.0E0_realk))THEN
        CALL DSCAL(B%localnrow*B%localncol,AB(2),B%p,1)
      endif
      if ((T(1)=='t').or.(T(1)=='T')) then
        CALL PDM_DSCINIT(DESC_A,A)
        CALL PDM_DSCINIT(DESC_C,C)
        CALL PDTRAN(C%nrow,C%ncol,1.0E0_realk,A%p,1,1,DESC_A,&
            & 0.0E0_realk,C%p,1,1,DESC_C)
      endif
      l=0
      do i=1,B%localnrow,BLOCK_SIZE
        if(i*BLOCK_SIZE > B%localnrow)then
          max_block_idx = mod(B%localnrow,BLOCK_SIZE)
        else
          max_block_idx=BLOCK_SIZE
        endif
        do j=1,max_block_idx
          glb=(l*SLGrid%nprow+SLGrid%myrow)*BLOCK_SIZE+j
          loc=l*BLOCK_SIZE + j
          if((T(1)=='n') .or. (T(1)=='N'))then
            CALL daxpy(B%localncol,diag(glb)*AB(1),A%p(loc,1),A%localnrow,B%p(loc,1),B%localnrow)
          else if((T(1)=='t').or.(T(1)=='T')) then
            CALL daxpy(B%localncol,diag(glb)*AB(1),C%p(loc,1),C%localnrow,B%p(loc,1),B%localnrow)
          endif
        enddo
        l=l+1
      enddo
      deallocate(diag)
   CASE(Job_density_from_orbs)
      CALL PDM_DSCINIT(DESC_A,A)
      CALL PDM_DSCINIT(DESC_B,B)
      call ls_mpibcast(nocc,infpar%master,scalapack_comm)
      CALL PDGEMM('N','T',A%nrow,A%ncol,nocc,1E0_realk,A%p,1,1,DESC_A,&
           A%p,1,1,DESC_A,0E0_realk,B%p,1,1,DESC_B)

   CASE(Job_max)
   ! using AF for val and tmp(1:2) for pos to reuse variables
   ! tmp(3) for comm and tmp(4) for ierr 
     AF=-huge(0.0_realk)

     !max localy
     DO j=1,A%localncol
        DO i=1,A%localnrow
           if (A%p(i,j) > AF) then
              AF = A%p(i,j)
              tmp(1)=i; tmp(2)=j
           endif
        ENDDO
     ENDDO

     !convert local positions into global positions
     tmp(1) = indxl2g(tmp(1),BLOCK_SIZE,SLGrid%myrow,0,SLGrid%nprow)
     tmp(2) = indxl2g(tmp(2),BLOCK_SIZE,SLGrid%mycol,0,SLGrid%npcol)
     tmpi=tmp

     !gather local max values and ther position
     call MPI_GATHER(AF,one,MPI_DOUBLE_PRECISION,AB,one,MPI_DOUBLE_PRECISION,zero,SLGrid%comm,tmpi(4))
     call MPI_GATHER(tmpi,two,MPI_INTEGER,n1,two,MPI_INTEGER,zero,SLGrid%comm,tmpi(4))
   CASE(Job_min)
   ! using AF for val and tmp(1:2) for pos to reuse variables
   ! tmp(3) for comm and tmp(4) for ierr 
     AF=huge(0.0_realk)

     !min localy
     DO j=1,A%localncol
        DO i=1,A%localnrow
           if (A%p(i,j) < AF) then
              AF = A%p(i,j)
              tmp(1)=i; tmp(2)=j
           endif
        ENDDO
     ENDDO

     !convert local positions into global positions
     tmp(1) = indxl2g(tmp(1),BLOCK_SIZE,SLGrid%myrow,0,SLGrid%nprow)
     tmp(2) = indxl2g(tmp(2),BLOCK_SIZE,SLGrid%mycol,0,SLGrid%npcol)
     tmpi=tmp

     !gather local max values and ther position
     call MPI_GATHER(AF,one,MPI_DOUBLE_PRECISION,AB,one,MPI_DOUBLE_PRECISION,zero,SLGrid%comm,tmpi(4))
     call MPI_GATHER(tmpi,two,MPI_INTEGER,n1,two,MPI_INTEGER,zero,SLGrid%comm,tmpi(4))

   CASE(Job_absmax)
      AF = 0.0E0_realk
      DO j=1,A%localncol
         DO i=1,A%localnrow
            AF = MAX(AF,ABS(A%p(i,j)))
         ENDDO
      ENDDO
      CALL lsmpi_max_realk_reduction(AF,infpar%master,scalapack_comm)
   CASE(Job_maxdiag)
      CALL PDM_DSCINIT(DESC_A,A)
      AB = 0.0E0_realk
      CALL PDLADIAG(A%nrow,A%p,1,1,DESC_A,PLUGIN_MAX,AB)
      CALL lsmpi_max_realk_reduction(AF,infpar%master,scalapack_comm)
   CASE(Job_outdia_sqnrm2)
      CALL PDM_DSCINIT(DESC_A,A)
      AB = 0.0E0_realk
      CALL PDLADIAG(A%nrow,A%p,1,1,DESC_A,PLUGIN_SQ2,AB)
      CALL DGSUM2D( SLGrid%ictxt,'A','I',1,1,AB,1,0,0)
   CASE(Job_extract_diag)
      CALL PDM_DSCINIT(DESC_A,A)
      ALLOCATE(diag(A%nrow))
      diag = 0.0E0_realk
      CALL PDLADIAG(A%nrow,A%p,1,1,DESC_A,PLUGIN_EXTRACT,diag)
      CALL DGSUM2D( SLGrid%ictxt,'A','I',A%nrow,1,diag,A%ncol,0,0)
      DEALLOCATE(diag)
   CASE(Job_add_identity)
      CALL PDM_DSCINIT(DESC_A,A)
      CALL DGEBR2D(SLGrid%ictxt,'A','I',2,1,AB,2,&
           &SLGrid%myrow,SLGrid%mycol)
      CALL PDLADIAG(A%nrow,A%p,1,1,DESC_A,PLUGIN_ADD,AB)
   CASE(Job_set)
      CALL PDM_DSCINIT(DESC_A,A)
      CALL DGEBR2D(SLGrid%ictxt,'A','I',2,1,AB,2,&
           &SLGrid%myrow,SLGrid%mycol)
      CALL PDLASET('A',A%nrow,A%ncol,AB(1),AB(2),A%p,1,1,DESC_A)
   CASE(Job_diag_f)
      CALL PDM_DSCINIT(DESC_A,A)
      CALL PDM_DSCINIT(DESC_B,B)
      CALL PDM_DSCINIT(DESC_C,C)
      ALLOCATE(diag(A%nrow))
      call mat_scalapack_diag_f_aux(A,B,diag,C,DESC_A,DESC_B,DESC_C)
      DEALLOCATE(diag)
   CASE(Job_dsyev)
      CALL PDM_DSCINIT(DESC_A,A)
      CALL PDM_DSCINIT(DESC_B,B)
      ALLOCATE(diag(A%nrow))
      call mat_scalapack_dsyev_aux(A,B,diag,A%nrow,DESC_A,DESC_B)
      DEALLOCATE(diag)
   CASE(Job_dsyevx)
      CALL PDM_DSCINIT(DESC_A,A)
      CALL ls_mpibcast(ieig,infpar%master,scalapack_comm)
      call mat_scalapack_dsyevx_aux(A,eival,ieig,DESC_A)
   CASE(Job_dpotrf)
      CALL PDM_DSCINIT(DESC_A,A)
      call mat_scalapack_dpotrf_aux(A,DESC_A)
   CASE(Job_dpotrs)
      CALL PDM_DSCINIT(DESC_A,A)
      CALL PDM_DSCINIT(DESC_B,B)
      call mat_scalapack_dpotrs_aux(A,B,DESC_A,DESC_B)
   CASE(Job_ao_precond)
      ALLOCATE(diag(A%nrow))
      ALLOCATE(diag2(A%nrow))
      CALL PDM_DSCINIT(DESC_A,A)
      CALL mat_scalapack_ao_precond_aux(A,A%nrow,DESC_A,symm,omega,diag,diag2)
      DEALLOCATE(diag)
      DEALLOCATE(diag2)
   CASE(Job_add_block)
      CALL PDM_DSCINIT(DESC_A, A)
      call ls_mpibcast(TMP,4,infpar%master,scalapack_comm)
      fullrow = TMP(1)
      fullcol = TMP(2)
      insertrow = TMP(3)
      insertcol = TMP(4)

      LLD =  MAX(1,NUMROC(fullrow, fullrow, SLGrid%myrow,0,SLGrid%nprow))
      CALL DESCINIT(DESC_AF,fullrow,fullcol,fullrow,fullcol,0,0,SLGrid%ictxt,LLD,INFO)

      CALL PDGEADD('N',fullrow,fullcol,1.0_realk,AF,1,1,DESC_AF, &
                 & 1.0_realk,A%p,insertrow,insertcol,DESC_A)

   CASE(Job_retrieve_block)
      CALL PDM_DSCINIT(DESC_A,A)
      call ls_mpibcast(TMP,4,infpar%master,scalapack_comm)
      fullrow = TMP(1)
      fullcol = TMP(2)
      insertrow = TMP(3)
      insertcol = TMP(4)

      LLD =  MAX(1,NUMROC(fullrow, fullrow, SLGrid%myrow,0,SLGrid%nprow))
      CALL DESCINIT(DESC_AF,fullrow,fullcol,fullrow,fullcol,0,0,SLGrid%ictxt,LLD,INFO)

      CALL PDGEMR2D(fullrow,fullcol,A%p,insertrow,insertcol,DESC_A, &
                   &AF,1,1,DESC_AF,SLGrid%ictxt)

   CASE(Job_create_block)
      CALL PDM_DSCINIT(DESC_A,A)
      call ls_mpibcast(TMP,4,infpar%master,scalapack_comm)
      fullrow = TMP(1)
      fullcol = TMP(2)
      insertrow = TMP(3)
      insertcol = TMP(4)

      LLD =  MAX(1,NUMROC(fullrow, fullrow, SLGrid%myrow,0,SLGrid%nprow))
      CALL DESCINIT(DESC_AF,fullrow,fullcol,fullrow,fullcol,0,0,SLGrid%ictxt,LLD,INFO)

      CALL PDGEMR2D(fullrow,fullcol,AF,1,1,DESC_AF, &
                   &A%p,insertrow,insertcol,DESC_A,SLGrid%ictxt)

   CASE(Job_print_global)
      call sleep(scalapack_mynum*5) !so first slave wait 5 sec, second slave wait for 10 sec
                                  !hopefully the matrix is then printet one block at a time
      call mat_scaLapack_print_local(A)
      call lsmpi_barrier(scalapack_comm)
   CASE(Job_setlowertriangular_zero)
      CALL PDM_DSCINIT(DESC_A, A)
      call scalapack_setlowertriangular_zero_aux(A,DESC_A)
   CASE(Job_scal_dia)
      CALL PDM_DSCINIT(DESC_A, A)
      call ls_mpibcast(AF,infpar%master,scalapack_comm)
      call scalapack_scal_dia_aux(AF,A,DESC_A)
   CASE(Job_scal_dia_vec)
      CALL PDM_DSCINIT(DESC_A, A)
      call mem_alloc(diag3,A%nrow)
      call ls_mpibcast(diag3,A%nrow,infpar%master,scalapack_comm)
      call scalapack_scal_dia_vec_aux(diag3,A,DESC_A,A%nrow)
      call mem_dealloc(diag3)
   CASE(Job_free)
      CALL FREE_IN_DARRAY(A)      
      nsize1=A%localnrow*A%localncol*mem_realsize
      nsize2=A%nrow*A%ncol*mem_realsize
      call mem_deallocated_mem_type_matrix(nsize1,nsize2)
   CASE DEFAULT
   END SELECT
   
#endif
 END SUBROUTINE PDM_SLAVE

 SUBROUTINE PDM_GRIDEXIT()
   use matrix_operations_scalapack
   use memory_handling
   use matrix_module
#ifdef VAR_MPI
   use infpar_module
   use lsmpi_type
#endif
   IMPLICIT NONE
   INTEGER :: IERR   
#ifdef VAR_SCALAPACK  
   call ls_mpibcast(GRIDEXIT,infpar%master,MPI_COMM_LSDALTON)
   ENTRY PDM_GRIDEXIT_SLAVE
   CALL BLACS_GRIDEXIT( SLGrid%ictxt )
   CALL BLACS_EXIT( 1 )
   CALL DARRAY_FREEALL
   CALL DARRAY_NULLIFY
#endif      
   call GRIDINFO_FREE()
 END SUBROUTINE PDM_GRIDEXIT

