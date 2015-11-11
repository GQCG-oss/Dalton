!> @file 
!> Contains essential matrix operations module.

!> Contains wrapper routines that branch out to matrix routine for chosen matrix type.
!> \author L. Thogersen
!> \date 2003
!>
!> General rules:
!> NEVER put a type(Matrix) as intent(out), this will on some platforms
!>       make the pointer
!>       disassociated entering the routine, and memory already
!>       allocated for the matrix will be lost. \n
!> NEVER think that e.g. A%elms = matrix will copy the matrix elements from
!>       matrix to A%elms, it will only associate the pointer with the array
!>       matrix. \n
!> Use mat_assign for A = B operations
!> ALWAYS and ONLY call mat_free on a matrix you have initialized with mat_init.
!>
MODULE matrix_operations
!FIXME: order routines alphabetically
!   use lstiming
   use matrix_module
   use matrix_operations_dense
   use matrix_operations_scalapack
   use matrix_operations_pdmm
   use LSTIMING
#ifdef VAR_MPI
   use lsmpi_type, only: MATRIXTY
#endif
   use matrix_operations_csr
   use matrix_op_unres_dense

   private
   public ::  matrixfiletype2, matrixfiletype, &
        & mtype_symm_dense, mtype_dense, mtype_unres_dense, mtype_csr,&
        & mtype_scalapack, mtype_pdmm, matrix_type, &
        & SET_MATRIX_DEFAULT, mat_select_type, mat_finalize, mat_pass_info,&
        & mat_timings, mat_no_of_matmuls, mat_select_tmp_type, &
        & mat_init, mat_free, allocated_memory,&
        & stat_deallocated_memory,mat_set_from_full,mat_to_full,mat_print,&
        & mat_trans,mat_chol,mat_dpotrf,mat_dpotrs,mat_dpotri,mat_inv,&
        & mat_clone,mat_assign,mat_mpicopy,mat_copy,mat_mul,mat_add,mat_daxpy,&
        & mat_dposv,mat_abs_max_elm,mat_max_elm,mat_min_elm,mat_max_diag_elm,&
        & mat_diag_f,mat_dsyev,mat_dsyevx,mat_section,mat_insert_section,&
        & mat_identity,mat_add_identity,mat_create_block,mat_add_block,&
        & mat_retrieve_block,mat_scal,mat_scal_dia,mat_scal_dia_vec,mat_zero,&
        & mat_setlowertriangular_zero,set_lowertriangular_zero,&
        & mat_write_to_disk,mat_write_info_to_disk,mat_read_from_disk,&
        & mat_read_info_from_disk,mat_extract_diagonal,&
        & no_of_matmuls, mat_tr, mat_trab, mat_dotproduct, mat_sqnorm2, &
        & mat_outdia_sqnorm2, info_memory, max_no_of_matrices, no_of_matrices,&
        & mat_to_full3D, string11_mtype,mat_add_to_fullunres

!> Matrices are symmetric and dense (not implemented)
   integer, parameter :: mtype_symm_dense  = 1
!> Matrices are dense (default) 
   integer, parameter :: mtype_dense       = 2
!> Matrices are dense and have both alpha and beta part (default for open shell)
   integer, parameter :: mtype_unres_dense = 5
!> Matrices are compressed sparse row (CSR) 
   integer, parameter :: mtype_csr         = 7
!> Matrices are MPI memory distributed using scalapack
   integer, parameter :: mtype_scalapack   = 8
!> Matrices are MPI memory distributed using TK scheme
   integer, parameter :: mtype_pdmm        = 9
!*****************
!Possible matrix types - 
!(Exploiting symmetry when operating sparse matrices is probably a vaste of effort - 
!therefore these combinations are removed)
!******************
!mtype_dense, matop_sparse2
!mtype_cplx_dense, mtype_cplx_sparse2
!mtype_unres_dense, mtype_unres_sparse2
!mtype_cplx_unres_dense, mtype_cplx_unres_sparse2
!mtype_symm_dense
!mtype_cplx_symm_dense
!mtype_unres_symm_dense
!mtype_cplx_unres_symm_dense

!> Counts the total number of matrix multiplications used throughout calculation
   integer, save      :: no_of_matmuls! = 0
!> Counts the number of matrices allocated, if requested
   integer, save      :: no_of_matrices! = 0 
!> Tracks the maximum number of allocated matrices throughout calculation, if requested
   integer, save      :: max_no_of_matrices! = 0
!> Start time for timing of matrix routines (only if requested)
   real(realk), save  :: mat_TSTR
!> End time for timing of matrix routines (only if requested)
   real(realk), save  :: mat_TEN
!> This is set to one of the mtype... variables to indicate chosen matrix type
   integer,save :: matrix_type ! = mtype_dense !default dense
!> True if timings for matrix operations are requested
   logical,save :: INFO_TIME_MAT! = .false. !default no timings
   logical,save :: INFO_memory! = .false. !default no memory printout

   contains
     subroutine string11_mtype(matrix_type_val,String)
       implicit none
       integer,intent(in) :: matrix_type_val
       character(len=11) :: String
       select case(matrix_type_val)
       case(mtype_symm_dense)
          String = 'symm_dense '
       case(mtype_dense)
          String = 'dense      '
       case(mtype_unres_dense)
          String = 'unres_dense'
       case(mtype_csr)
          String = 'csr        '
       case(mtype_scalapack)
          String = 'scalapack  '
#ifdef VAR_ENABLE_TENSORS
       case(mtype_pdmm)
          String = 'pdmm       '
#endif
       case default
          call lsquit("Unknown type of matrix",-1)
       end select
     end subroutine string11_mtype

!*** is called from LSDALTON.f90
!> \brief Sets the global variables
!> \author T. Kjaergaard
!> \date 2011
     SUBROUTINE SET_MATRIX_DEFAULT()
       implicit none

       no_of_matmuls = 0
       no_of_matrices = 0 
       max_no_of_matrices = 0
       matrix_type = mtype_dense !default dense
       INFO_TIME_MAT = .false. !default no timings
       INFO_memory = .false. !default no memory printout  
     END SUBROUTINE SET_MATRIX_DEFAULT

!*** is called from config.f90
!> \brief Sets the global variable matrix_type that determines the matrix type
!> \author L. Thogersen
!> \date 2003
!> \param a Indicates the matrix type (see module documentation) 
     SUBROUTINE mat_select_type(a,lupri,nbast)
#ifdef VAR_MPI
       use infpar_module
       use lsmpi_type
#endif
       implicit none
       INTEGER, INTENT(IN) :: a,lupri
       INTEGER, OPTIONAL :: nbast
       integer :: nrow,ncol,tmpcol,tmprow,nproc,K,I,nblocks
       WRITE(lupri,'(A)') ' '
       if(matrix_type.EQ.mtype_unres_dense.AND.a.EQ.mtype_scalapack)then
          WRITE(6,*)'mat_select_type: FALLBACK WARNING'
          WRITE(6,*)'SCALAPACK type matrices is not implemented for unrestricted calculations'
          WRITE(6,*)'We therefore use the dense unrestricted type - which do not use memory distribution'
       else
          matrix_type = a
#ifdef VAR_MPI
          IF (infpar%mynum.EQ.infpar%master) THEN
#endif
            select case(matrix_type)             
            case(mtype_dense)
               WRITE(lupri,'(A)') 'Matrix type: mtype_dense'
            case(mtype_unres_dense)
               WRITE(lupri,'(A)') 'Matrix type: mtype_unres_dense'
            case(mtype_csr)
               WRITE(lupri,'(A)') 'Matrix type: mtype_csr'
            case(mtype_scalapack)
               WRITE(lupri,'(A)') 'Matrix type: mtype_scalapack'
#ifdef VAR_ENABLE_TENSORS
            case(mtype_pdmm)
               WRITE(lupri,'(A)') 'Matrix type: mtype_pdmm'
#endif
            case default
               call lsquit("Unknown type of matrix",-1)
            end select
#ifdef VAR_MPI
          ENDIF
          IF (infpar%mynum.EQ.infpar%master) THEN
             call ls_mpibcast(MATRIXTY,infpar%master,MPI_COMM_LSDALTON)
             call lsmpi_set_matrix_type_master(a)
             if(matrix_type.EQ.mtype_scalapack)then
                IF(.NOT.present(nbast))then
                   call lsquit('scalapack error in mat_select_type',-1)
                ENDIF
                print*,'TYPE SCALAPACK HAVE BEEN SELECTED' 
                
                print*,'infpar%ScalapackGroupSize',infpar%ScalapackGroupSize
                print*,'infpar%inputblocksize',infpar%inputblocksize
                IF(infpar%ScalapackGroupSize.EQ.-1)THEN !auto 
                   print*,'Automatic determination of how many scalapack groupsize' 
                   IF(infpar%inputblocksize.EQ.0)THEN
                      !blocksize not specified. 
                      !determine block size. 
                      !The optimal is a blocksize big enough to reduce communication and optimize lapack
                      call DetermineBlockSize(infpar%inputblocksize,infpar%nodtot,nbast,lupri)
                      WRITE(lupri,'(A,I5)')'Automatic determined Scalapack BlockSize = ',infpar%inputblocksize
                   ELSE
                      WRITE(lupri,'(A,I5)')'Scalapack BlockSize From Input = ',infpar%inputblocksize
                   ENDIF
                   nblocks = (nbast/infpar%inputblocksize)*(nbast/infpar%inputblocksize)
                   IF(MOD(nbast,infpar%inputblocksize).NE.0)THEN
                      nblocks = nblocks + 2*(nbast/infpar%inputblocksize) + 1
                   ENDIF
                   WRITE(lupri,'(A,I5)')'Number of Scalapack Blocks = ',nblocks
                   if(matrix_type.EQ.mtype_scalapack)then
                      IF(nblocks .GE. 2*infpar%nodtot)THEN
                         !2 or more blocks per node, we can use all nodes
                         infpar%ScalapackGroupSize = infpar%nodtot
                      ELSE
                         !not all nodes have 2 blocks. We cannot use all nodes efficiently
                         !We decide on a groupsize so that all nodes have 2 blocks
                         IF(nblocks/2.GT.0)THEN
                            infpar%ScalapackGroupSize = nblocks/2
                         ELSE
                            infpar%ScalapackGroupSize = nblocks
                         ENDIF
                      ENDIF
                      WRITE(lupri,'(A,I5)')'Automatic determined Scalapack GroupSize = ',infpar%ScalapackGroupSize
                      nproc = infpar%ScalapackGroupSize
                   else
                      nproc = infpar%nodtot
                   endif
                ELSE !infpar%ScalapackGroupSize.NE.-1 
                   WRITE(lupri,'(A,I5)')'Scalapack GroupSize From Input = ',infpar%ScalapackGroupSize
                   print*,'Scalapack GroupSize From Input = ',infpar%ScalapackGroupSize
                   print*,'infpar%inputblocksize',infpar%inputblocksize
                   nproc = infpar%ScalapackGroupSize
                ENDIF
                nrow = nproc
                ncol = 1
                K=1
                do 
                   K=K+1
                   IF(nproc+1.LE.K*K)EXIT
                   tmprow = nproc/K
                   tmpcol = K
                   IF(tmprow*tmpcol.EQ.nproc)THEN
                      IF(tmprow+tmpcol.LE.nrow+nrow)THEN
                         nrow = tmprow
                         ncol = tmpcol
                      ENDIF
                   ENDIF
                enddo
                if(matrix_type.EQ.mtype_scalapack)then
                   print*,'nrow=',nrow,'ncol=',ncol,'nodtot=',infpar%ScalapackGroupSize
                endif
                print*,'call PDM_GRIDINIT(',nrow,',',ncol,')'                
                !make possible subset 
                print*,'call ls_mpibcast with infpar%ScalapackGroupSize',infpar%ScalapackGroupSize
                call ls_mpibcast(infpar%ScalapackGroupSize,infpar%master,MPI_COMM_LSDALTON)                
                call ls_mpibcast(infpar%ScalapackWorkaround,infpar%master,MPI_COMM_LSDALTON)                
                print*,'master done bcast'
                IF(infpar%ScalapackGroupSize.NE.infpar%nodtot)THEN
                   IF(scalapack_mpi_set)THEN
                      !free communicator 
                      call LSMPI_COMM_FREE(scalapack_comm)
                      scalapack_mpi_set = .FALSE.
                   ENDIF
                   call init_mpi_subgroup(scalapack_nodtot,&
                        & scalapack_mynum,scalapack_comm,scalapack_member,&
                        & infpar%ScalapackGroupSize,MPI_COMM_LSDALTON,lupri)
                   scalapack_mpi_set = .TRUE.
                ELSE
                   scalapack_nodtot = infpar%nodtot
                   scalapack_mynum = infpar%mynum
                   scalapack_comm = MPI_COMM_LSDALTON
                   scalapack_member = .TRUE.
                ENDIF
                print*,'call PDM_GRIDINIT'
                CALL PDM_GRIDINIT(nrow,ncol,nbast)
                WRITE(lupri,'(A,I5)')'Scalapack Grid initiation Block Size = ',BLOCK_SIZE
                WRITE(lupri,'(A,I5)')'Scalapack Grid initiation nprow      = ',nrow
                WRITE(lupri,'(A,I5)')'Scalapack Grid initiation npcol      = ',ncol
#ifdef VAR_ENABLE_TENSORS
             elseif(matrix_type.EQ.mtype_pdmm)then
                IF(.NOT.present(nbast))then
                   call lsquit('pdmm error in mat_select_type',-1)
                ENDIF
                print*,'TYPE PDMM HAVE BEEN SELECTED' 
                !make possible subset 
                call ls_mpibcast(infpar%PDMMGroupSize,infpar%master,MPI_COMM_LSDALTON)                
                IF(infpar%PDMMGroupSize.NE.infpar%nodtot)THEN
                   IF(pdmm_mpi_set)THEN
                      !free communicator 
                      call LSMPI_COMM_FREE(pdmm_comm)
                      pdmm_mpi_set = .FALSE.
                   ENDIF
                   call init_mpi_subgroup(pdmm_nodtot,&
                        & pdmm_mynum,pdmm_comm,pdmm_member,&
                        & infpar%PDMMGroupSize,MPI_COMM_LSDALTON,lupri)
                   pdmm_mpi_set = .TRUE.
                ELSE
                   pdmm_nodtot = infpar%nodtot
                   pdmm_mynum = infpar%mynum
                   pdmm_comm = MPI_COMM_LSDALTON
                   pdmm_member = .TRUE.
                ENDIF
                CALL PDMM_GRIDINIT(nbast)
                WRITE(lupri,'(A,I5)')'PDMM initiation Block Size = ',BLOCK_SIZE_PDM
#endif
             endif
          ELSE
             !slave
             if(matrix_type.EQ.mtype_scalapack)then
                call ls_mpibcast(infpar%ScalapackGroupSize,infpar%master,&
                     & MPI_COMM_LSDALTON)
                call ls_mpibcast(infpar%ScalapackWorkaround,infpar%master,&
                     & MPI_COMM_LSDALTON)
                infpar%inputblocksize = infpar%ScalapackGroupSize
                IF(infpar%ScalapackGroupSize.NE.infpar%nodtot)THEN
                   IF(scalapack_mpi_set)THEN
                      !free communicator 
                      call LSMPI_COMM_FREE(scalapack_comm)
                      scalapack_mpi_set = .FALSE.
                   ENDIF
                   call init_mpi_subgroup(scalapack_nodtot,&
                        & scalapack_mynum,scalapack_comm,scalapack_member,&
                        & infpar%ScalapackGroupSize,MPI_COMM_LSDALTON,lupri)
                   scalapack_mpi_set = .TRUE.
                ELSE
                   scalapack_nodtot = infpar%nodtot
                   scalapack_mynum = infpar%mynum
                   scalapack_comm = MPI_COMM_LSDALTON
                   scalapack_member = .TRUE.
                ENDIF
#ifdef VAR_ENABLE_TENSORS
             elseif(matrix_type.EQ.mtype_pdmm)then
                call ls_mpibcast(infpar%PDMMGroupSize,infpar%master,&
                     & MPI_COMM_LSDALTON)
                IF(infpar%PDMMGroupSize.NE.infpar%nodtot)THEN
                   IF(PDMM_mpi_set)THEN
                      !free communicator 
                      call LSMPI_COMM_FREE(pdmm_comm)
                      PDMM_mpi_set = .FALSE.
                   ENDIF
                   call init_mpi_subgroup(pdmm_nodtot,&
                        & pdmm_mynum,pdmm_comm,pdmm_member,&
                        & infpar%PDMMGroupSize,MPI_COMM_LSDALTON,lupri)
                   PDMM_mpi_set = .TRUE.
                ELSE
                   pdmm_nodtot = infpar%nodtot
                   pdmm_mynum = infpar%mynum
                   pdmm_comm = MPI_COMM_LSDALTON
                   pdmm_member = .TRUE.
                ENDIF
#endif
             endif
          ENDIF
#endif
       endif
     END SUBROUTINE mat_select_type


!> \brief Sets the global variable matrix_type that determines the matrix type
!> It does not init SCALAPACK or PDMM module
!> \author T. Kjaergaard
!> \date 2015
!> \param a Indicates the matrix type (see module documentation) 
     SUBROUTINE mat_select_tmp_type(a,lupri,nbast)
#ifdef VAR_MPI
       use infpar_module
       use lsmpi_type
#endif
       implicit none
       INTEGER, INTENT(IN) :: a,lupri
       INTEGER, OPTIONAL :: nbast
       integer :: nrow,ncol,tmpcol,tmprow,nproc,K,I,nblocks
       WRITE(lupri,'(A)') ' '
       if(matrix_type.EQ.mtype_unres_dense.AND.a.EQ.mtype_scalapack)then
          WRITE(6,*)'mat_select_type: FALLBACK WARNING'
          WRITE(6,*)'SCALAPACK type matrices is not implemented for unrestricted calculations'
          WRITE(6,*)'We therefore use the dense unrestricted type - which do not use memory distribution'
       else
          matrix_type = a
#ifdef VAR_MPI
          IF (infpar%mynum.EQ.infpar%master) THEN
             call ls_mpibcast(MATRIXTY2,infpar%master,MPI_COMM_LSDALTON)
             call lsmpi_set_matrix_type_master(a)
          ENDIF
#endif
       endif
     END SUBROUTINE mat_select_tmp_type

     subroutine DetermineBlockSize(blocksize,nodtot,nbast,lupri)
       implicit none
       integer,intent(inout) :: blocksize
       integer(kind=ls_mpik),intent(in) :: nodtot
       integer,intent(in) :: nbast,lupri
       integer,parameter :: blocklist(4)=(/1024,512,256,128/)
       !
       integer :: I,nblocks
                        
       do I =1,size(blocklist)
          blocksize = blocklist(I)
          IF(blocksize.GT.nbast)CYCLE
          nblocks = (nbast/blocksize)*(nbast/blocksize)
          IF(MOD(nbast,blocksize).NE.0)THEN
             nblocks = nblocks + 2*nbast/blocksize + 1
          ENDIF
          WRITE(lupri,'(A,I5,A,I6,A)')'A BlockSize of ',blocksize,' gives ',nblocks,' blocks'
!          WRITE(lupri,'(A,F9.1,A)') 'Resulting in ',(nblocks*1.0E0_realk)/infpar%nodtot,' blocks per node'
          IF(nblocks .GE. 2*nodtot)THEN
             !more than 2 blocks per node, we can use all nodes efficiently
             exit
          ELSE
             !too few blocks not all nodes have 2 blocks
             !We use a smaller block size, to reduce load imbalance
          ENDIF
          IF(I.EQ.size(blocklist).AND.(nblocks .GE. 2*nodtot))THEN
             WRITE(lupri,'(A)')'The Minimum BlockSize = 128 Chosen.'
             WRITE(lupri,'(A)')'This means that not all processes can take part in'
             WRITE(lupri,'(A)')'the matrix operation parallelization. '
             WRITE(lupri,'(A)')'You can manually set the BlockSize using the .SCALAPACKBLOCKSIZE keyword'
             WRITE(lupri,'(A)')'Under the **GENERAL section. A smaller BlockSize will include more nodes'
             WRITE(lupri,'(A)')'and improve load imbalance,'
             WRITE(lupri,'(A)')'but will reduce the efficiency of the underlying Lapack Lib.'
          ENDIF
       enddo
       IF(blocksize.GT.nbast)THEN
          WRITE(lupri,'(A)')'Warning: Due to the small size of matrices Scalapack is not recommended!'
          WRITE(lupri,'(A)')'We set the blocksize equal to the number of basis functions'
          blocksize = nbast
       ENDIF
     end subroutine DetermineBlockSize

     SUBROUTINE mat_finalize()
#ifdef VAR_MPI
       use infpar_module
       use lsmpi_type
#endif
       implicit none
#ifdef VAR_MPI
       if(matrix_type.EQ.mtype_scalapack)then
          CALL PDM_GRIDEXIT
       endif
#ifdef VAR_ENABLE_TENSORS
       if(matrix_type.EQ.mtype_pdmm)then
          CALL PDMM_GRIDEXIT
       endif
#endif

#endif
     END SUBROUTINE mat_finalize

!> \brief Pass info about e.g. logical unit number for LSDALTON.OUT to matrix module 
!> \author L. Thogersen
!> \date 2003
!> \param lu_info Logical unit number for LSDALTON.OUT
!> \param info_info True if various info from matrix module should be printed
!> \param mem_monitor True if number of allocated matrices should be monitored
      SUBROUTINE mat_pass_info(lu_info,info_info,mem_monitor)
         implicit none
         integer, intent(in) :: lu_info
         logical, intent(in) :: info_info, mem_monitor
         mat_lu   = lu_info
         mat_info = info_info
         mat_mem_monitor = mem_monitor
      END SUBROUTINE mat_pass_info
!> \brief If called, timings from matrix routines will be printed 
!> \author L. Thogersen
!> \date 2003
      SUBROUTINE mat_timings
         implicit none
         INFO_TIME_MAT = .true. 
      END SUBROUTINE mat_timings
!***
!> \brief Returns the number of matrix multiplications used so far 
!> \author S. Host
!> \date 2009
!> \param n Number of matrix muliplications
      SUBROUTINE mat_no_of_matmuls(n)
         implicit none
         INTEGER, INTENT(out) :: n
         n = no_of_matmuls
      END SUBROUTINE mat_no_of_matmuls
!> \brief Initialize a type(matrix)
!> \author L. Thogersen
!> \date 2003
!> \param a type(matrix) that should be initialized
!> \param nrow Number of rows for a
!> \param ncol Number of columns for a
      SUBROUTINE mat_init(a,nrow,ncol,complex)
         implicit none
         TYPE(Matrix), TARGET :: a 
         INTEGER, INTENT(IN)  :: nrow, ncol
         LOGICAL, INTENT(IN), OPTIONAL :: complex
         ! Always start by nullifying matrix components
         call mat_nullify(a)

         !if 'a' has init tag AND self pointer, it means that it is already initialized
         !and in its original location. Re-initializing would leak memory, so err
!         if (a%init_magic_tag.EQ.mat_init_magic_value &
!             & .and. associated(a%init_self_ptr,a)) THEN
!            print*,'associated(a%init_self_ptr,a)',associated(a%init_self_ptr,a)
!            print*,'a%init_magic_tag',a%init_magic_tag
!            print*,'mat_init_magic_value',mat_init_magic_value
!            print*,'a%init_magic_tag.EQ.mat_init_magic_value',a%init_magic_tag.EQ.mat_init_magic_value
!            call lsQUIT('Error in mat_init: matrix is already initialized',-1)
!         endif
         a%init_magic_tag = mat_init_magic_value
         a%init_self_ptr => a
         !process optional complex
         if (present(complex)) then
            if (complex .and. matrix_type.NE.mtype_dense) &
               & call lsQUIT('Error in mat_init: complex only implemented for matrix type dense',-1)
            a%complex = complex
         else
            a%complex = .false.
         endif
         !record this mat_init in statistics
         if (mat_mem_monitor) then
            no_of_matrices = no_of_matrices + 1
            !write(mat_lu,*) 'Init: matrices allocated:', no_of_matrices
            if (no_of_matrices > max_no_of_matrices)then
               max_no_of_matrices = no_of_matrices!
!               WRITE(mat_lu,*)'increase max_no_of_matrices to ',no_of_matrices
!               call LsTraceBack('increase max_no_of_matrices')
            endif
         endif
         nullify(A%elms)
         nullify(A%elmsb)
         if (info_memory) write(mat_lu,*) 'Before mat_init: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_init(a,nrow,ncol)
         case(mtype_unres_dense)
             call mat_unres_dense_init(a,nrow,ncol)
         case(mtype_csr)
             call mat_csr_init(a,nrow,ncol)            
         case(mtype_scalapack)
             call mat_scalapack_init(a,nrow,ncol)            
#ifdef VAR_ENABLE_TENSORS
         case(mtype_pdmm)
             call mat_pdmm_init(a,nrow,ncol)            
#endif
         case default
              call lsquit("mat_init not implemented for this type of matrix",-1)
         end select

         NULLIFY(a%iaux, a%raux)
         if (info_memory) write(mat_lu,*) 'After mat_init: mem_allocated_global =', mem_allocated_global
      END SUBROUTINE mat_init

!> \brief Free a type(matrix) that has been initialized with mat_init
!> \author L. Thogersen
!> \date 2003
!> \param a type(matrix) that should be freed
      SUBROUTINE mat_free(a,matype)
         implicit none
         TYPE(Matrix), TARGET :: a 
         integer,intent(in),optional :: matype
         integer   :: matrix_type_case

         !to be free'ed, the matrix must be in the same location where it was init'ed.
         !If not, it is probably a duplicate (like 'a' in (/a,b,c/)), in which case
         !we may end up double-free'ing, so err
         !if (.not.ASSOCIATED(a%init_self_ptr,a)) &
         !    & call lsQUIT('Error in mat_free: matrix moved or duplicated',-1)
         nullify(a%init_self_ptr)
         !look at magic tag to verify matrix is initialized, then clear tag
         if (a%init_magic_tag.NE.mat_init_magic_value) &
             & call lsQUIT('Error in mat_free: matrix was not initialized',-1)
         a%init_magic_tag = 0
         !record this mat_free in statistics
         if (mat_mem_monitor) then
            no_of_matrices = no_of_matrices - 1
            !write(mat_lu,*) 'Free: matrices allocated:', no_of_matrices
         endif
         if (info_memory) write(mat_lu,*) 'Before mat_free: mem_allocated_global =', mem_allocated_global
         IF(present(matype))THEN
            matrix_type_case = matype
         ELSE
            matrix_type_case = matrix_type
         ENDIF
         select case(matrix_type_case)
         case(mtype_dense)
             call mat_dense_free(a)
         case(mtype_unres_dense)
             call mat_unres_dense_free(a)
         case(mtype_csr)
             call mat_csr_free(a)
         case(mtype_scalapack)
             call mat_scalapack_free(a)
#ifdef VAR_ENABLE_TENSORS
         case(mtype_pdmm)
             call mat_pdmm_free(a)
#endif
         case default
            call lsquit("mat_free not implemented for this type of matrix",-1)
         end select

      !free auxaliary data
      if (ASSOCIATED(a%iaux)) deallocate(a%iaux)
      if (ASSOCIATED(a%raux)) deallocate(a%raux)

      a%nrow = -1; a%ncol = -1
      if (info_memory) write(mat_lu,*) 'After mat_free: mem_allocated_global =', mem_allocated_global

      END SUBROUTINE mat_free

!> \brief Count allocated memory for type(matrix)
!> \author L. Thogersen
!> \date 2003
!> \param nsize Number of real(realk) elements that have been allocated
      SUBROUTINE stat_allocated_memory(nsize)
         implicit none
         integer, intent(in) :: nsize 
         integer(kind=long) :: nsize2 
         select case(matrix_type)
         case(mtype_dense)
            nsize2 = nsize
            call mem_allocated_mem_type_matrix(nsize2)
         case(mtype_unres_dense)
             call unres_dens_stat_allocated_mem(nsize)
         case default
              call lsquit("stat_allocated_memory not implemented for this type of matrix",-1)
         end select

      END SUBROUTINE stat_allocated_memory

!> \brief Count deallocated memory for type(matrix)
!> \author L. Thogersen
!> \date 2003
!> \param nsize Number of real(realk) elements that have been deallocated
      SUBROUTINE stat_deallocated_memory(nsize)
         implicit none
         integer, intent(in) :: nsize 
         integer(kind=long) :: nsize2 
         select case(matrix_type)
         case(mtype_dense)
!             call dens_stat_deallocated_memory(nsize)
            nsize2 = nsize
            call mem_deallocated_mem_type_matrix(nsize2)
         case(mtype_unres_dense)
             call unres_dens_stat_deallocated_mem(nsize)
         case default
              call lsquit("stat_deallocated_memory not implemented for this type of matrix",-1)
         end select

      END SUBROUTINE stat_deallocated_memory

!> \brief Convert a standard fortran matrix to a type(matrix) - USAGE DISCOURAGED!
!> \author L. Thogersen
!> \date 2003
!> \param afull Standard fortran matrix that should be converted (n x n)
!> \param alpha The output type(matrix) is multiplied by alpha
!> \param a The output type(matrix) ((2n x 2n) if unrestricted, (n x n) otherwise)
!> \param mat_label If the character label is present, sparsity will be printed if using block-sparse matrices
!> \param unres3 If present and true, mat_unres_dense_set_from_full3 will be called instead of mat_unres_dense_set_from_full
!>  
!> BE VERY CAREFUL WHEN USING mat_set_from_full AND mat_to_full!!!!!!
!> Usage of these routines should be avoided whenever possible, since
!> you have to hardcode an interface to make them work with unrestriced
!> matrices (see e.g. di_get_fock in dalton_interface.f90) This is because
!> usually, a and afull should have the same dimensions, but for unrestricted
!> a is (n x n) and afull is (2n x 2n). The exception is if
!> unres3 = true. In that case, both a and afull are (n x n), and both 
!> the alpha and beta parts of a will be set equal to the afull.
!> 
      SUBROUTINE mat_set_from_full(afull,alpha, a, mat_label,unres3)
         implicit none
         real(realk), INTENT(IN) :: afull(*)
         real(realk), intent(in) :: alpha
         TYPE(Matrix)            :: a  !output
         character(*), INTENT(IN), OPTIONAL :: mat_label
         logical,intent(in) , OPTIONAL :: unres3
         real(realk)             :: sparsity
         real(realk),pointer     :: full(:,:)
         call time_mat_operations1
         !write(mat_lu,*) "Usage of mat_set_from_full discouraged!!!"
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'Before mat_set_from_full: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_set_from_full(afull,alpha,a)
         case(mtype_csr)
             call mat_csr_set_from_full(afull,alpha,a)
         case(mtype_scalapack)
            call mat_scalapack_set_from_full(afull,alpha,a)
#ifdef VAR_ENABLE_TENSORS
         case(mtype_pdmm)
            call mat_pdmm_set_from_full(afull,alpha,a)
#endif
         case(mtype_unres_dense)
            if(PRESENT(unres3))then
               if(unres3)then
                  call mat_unres_dense_set_from_full3(afull,alpha,a)
               else
                  call mat_unres_dense_set_from_full(afull,alpha,a)
               endif
            else
               call mat_unres_dense_set_from_full(afull,alpha,a)
            endif
         case default
              call lsquit("mat_set_from_full not implemented for this type of matrix",-1)
         end select
         if (info_memory) write(mat_lu,*) 'After mat_set_from_full: mem_allocated_global =', mem_allocated_global
         !if (INFO_TIME_MAT) CALL LSTIMER('F_FULL',mat_TSTR,mat_TEN,mat_lu)
         call time_mat_operations2(JOB_mat_set_from_full)

      END SUBROUTINE mat_set_from_full


!> \brief Convert a type(matrix) to a standard fortran matrix - USAGE DISCOURAGED!
!> \author L. Thogersen
!> \date 2003
!> \param a The type(matrix) that should be converted (n x n)
!> \param afull The output standard fortran matrix ((2n x 2n) if unrestricted, (n x n) otherwise)
!> \param alpha The output standard fortran matrix is multiplied by alpha
!> \param mat_label If the character label is present, sparsity will be printed if using block-sparse matrices
!>  
!> BE VERY CAREFUL WHEN USING mat_set_from_full AND mat_to_full!!!!!!
!> Usage of these routines should be avoided whenever possible, since
!> you have to hardcode an interface to make them work with unrestriced
!> matrices (see e.g. di_get_fock in dalton_interface.f90)
!>
     SUBROUTINE mat_to_full(a, alpha, afull,mat_label,matype)

         implicit none
         TYPE(Matrix), intent(in):: a
         real(realk), intent(in) :: alpha
         real(realk), intent(inout):: afull(*)  !output
         character(*), INTENT(IN), OPTIONAL :: mat_label
         integer,intent(in),optional :: matype
         integer                 :: matrix_type_case
         real(realk)             :: sparsity
         call time_mat_operations1

         !write(mat_lu,*) "Usage of mat_to_full discouraged!!!"
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         !if (SIZE(afull) < a%nrow*a%ncol) then
         !  call lsquit('too small full array in mat_to_full',-1)
         !endif
         if (info_memory) write(mat_lu,*) 'Before mat_to_full: mem_allocated_global =', mem_allocated_global
         IF(present(matype))THEN
            matrix_type_case = matype
         ELSE
            matrix_type_case = matrix_type
         ENDIF
         select case(matrix_type_case)
         case(mtype_dense)
             call mat_dense_to_full(a, alpha, afull)
         case(mtype_csr)
             call mat_csr_to_full(a, alpha, afull)
         case(mtype_scalapack)
            call mat_scalapack_to_full(a, alpha, afull)
#ifdef VAR_ENABLE_TENSORS
         case(mtype_pdmm)
            call mat_pdmm_to_full(a, alpha, afull)
#endif
         case(mtype_unres_dense)
             call mat_unres_dense_to_full(a, alpha, afull)
         case default
              call lsquit("mat_to_full not implemented for this type of matrix",-1)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('TOFULL',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_to_full: mem_allocated_global =', mem_allocated_global
         call time_mat_operations2(JOB_mat_to_full)
      END SUBROUTINE mat_to_full


!> \brief Adds an alpha or beta part of a unrestricted matrix type to 
!>  a standard fortran matrix
!> \author T. Kjaergaard
!> \date 2015
!> \param a The type(matrix) that should be converted (n x n)
!> \param afull The output standard fortran matrix (n x n) 
!> \param alpha The output standard fortran matrix is multiplied by alpha
      SUBROUTINE mat_add_to_fullunres(a, alpha, afull,spin)
        implicit none
        TYPE(Matrix), intent(in):: a
        real(realk), intent(in) :: alpha
        real(realk), intent(inout):: afull(*)  !output
        integer,intent(in) :: spin
        call time_mat_operations1
        if (info_memory) write(mat_lu,*) 'Before mat_to_full: mem_allocated_global =', mem_allocated_global
        select case(matrix_type)
        case(mtype_unres_dense) 
           call mat_unres_dense_add_to_fullunres(a, alpha, afull,spin)
           !FIXME BYKOV GENERALIZE TO THE OTHER TYPES!
        case default
           call lsquit("mat_to_full2 not implemented for this type of matrix",-1)
        end select
        !if (INFO_TIME_MAT) CALL LSTIMER('TOFULL',mat_TSTR,mat_TEN,mat_lu)
        if (info_memory) write(mat_lu,*) 'After mat_to_full: mem_allocated_global =', mem_allocated_global
        call time_mat_operations2(JOB_mat_to_full)
      END SUBROUTINE mat_add_to_fullunres

!> \brief Convert a type(matrix) to a standard 3D fortran matrix - USAGE DISCOURAGED!
!> \author L. Thogersen
!> \date 2003
!> \param a The type(matrix) that should be converted (n x n)
!> \param afull The output standard fortran matrix ((2n x 2n) if unrestricted, (n x n) otherwise)
!> \param alpha The output standard fortran matrix is multiplied by alpha
!> \param mat_label If the character label is present, sparsity will be printed if using block-sparse matrices
!>  
!> BE VERY CAREFUL WHEN USING mat_set_from_full AND mat_to_full!!!!!!
!> Usage of these routines should be avoided whenever possible, since
!> you have to hardcode an interface to make them work with unrestriced
!> matrices (see e.g. di_get_fock in dalton_interface.f90)
!>
     SUBROUTINE mat_to_full3D(a, alpha, afull,n1,n2,n3,i3,i4)

         implicit none
         integer, INTENT(IN)           :: n1,n2,n3,i3
         integer, INTENT(IN), OPTIONAL :: i4 !for unres
         TYPE(Matrix), intent(in):: a
         real(realk), intent(in) :: alpha
         real(realk), intent(inout):: afull(n1,n2,n3)  !output
         integer :: local_i4
         IF(present(i4))THEN
            local_i4=i4 
            !this means place alpha in A(:,:,i3) place beta in A(:,:,i4)
         ELSE
            local_i4=i3 !this means add alpha and beta for unres 
         ENDIF
         call time_mat_operations1
         if (info_memory) write(mat_lu,*) 'Before mat_to_full: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_to_full3d(a, alpha, afull,n1,n2,n3,i3)
         case(mtype_csr)
             call mat_csr_to_full3d(a, alpha, afull,n1,n2,n3,i3)
         case(mtype_scalapack)
            call mat_scalapack_to_full3d(a, alpha, afull,n1,n2,n3,i3)
         case(mtype_unres_dense)
             call mat_unres_dense_to_full3d(a, alpha, afull,n1,n2,n3,i3,local_i4)
         case default
              call lsquit("mat_to_full3D not implemented for this type of matrix",-1)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('TOFULL',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_to_full: mem_allocated_global =', mem_allocated_global
         call time_mat_operations2(JOB_mat_to_full)
       END SUBROUTINE mat_to_full3D

!> \brief Print a type(matrix) to file in pretty format
!> \author L. Thogersen
!> \date 2003
!> \param a The type(matrix) that should be printed
!> \param i_row1 Print starting from this row
!> \param i_rown Print ending at this row
!> \param j_col1 Print starting from this column
!> \param j_coln Print ending at this column
!> \param lu Print to file with this logical unit number
      SUBROUTINE mat_print(a, i_row1, i_rown, j_col1, j_coln, lu)
         implicit none
         TYPE(Matrix),intent(in) :: a
         integer, intent(in)     :: i_row1, i_rown, j_col1, j_coln, lu 
         REAL(REALK), ALLOCATABLE :: afull(:,:)
         real(realk)              :: sparsity
         if (i_row1 < 1 .or. j_col1 < 1 .or. a%nrow < i_rown .or. a%ncol < j_coln) then
           CALL LSQUIT( 'subsection out of bounds in mat_print',-1)
         endif
         if (info_memory) write(mat_lu,*) 'Before mat_print: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_print(a, i_row1, i_rown, j_col1, j_coln, lu)
#ifdef VAR_ENABLE_TENSORS
         case(mtype_pdmm)
            call mat_pdmm_print(a, i_row1, i_rown, j_col1, j_coln, lu)
#endif
         case(mtype_scalapack)
#ifdef VAR_SCALAPACK
            print*,'FALLBACK scalapack print'
            ALLOCATE (afull(a%nrow,a%ncol))
            call mat_scalapack_to_full(a, 1E0_realk,afull)
            CALL LS_OUTPUT(afull, i_row1, i_rown, j_col1, j_coln,A%nrow,A%ncol,1, lu)
            DEALLOCATE(afull)
#endif
         case(mtype_unres_dense)
             call mat_unres_dense_print(a, i_row1, i_rown, j_col1, j_coln, lu)
         case(mtype_csr)
             call mat_csr_print(a, lu)
         case default
              call lsquit("mat_print not implemented for this type of matrix",-1)
         end select
         if (info_memory) write(mat_lu,*) 'After mat_print: mem_allocated_global =', mem_allocated_global
      END SUBROUTINE mat_print

!> \brief Transpose a type(matrix).
!> \author L. Thogersen
!> \date 2003
!> \param a The type(matrix) that should be transposed
!> \param b The transposed output type(matrix).
!>
!> Usage discouraged! If what you want is to multiply your transposed
!> matrix with something else, you should instead use mat_mul with the
!> transpose flag 'T'. This is much more efficient than transposing first 
!> and then multiplying.
!>
      SUBROUTINE mat_trans(a, b) !USAGE DISCOURAGED!!
         implicit none
         TYPE(Matrix),intent(in)     :: a
         TYPE(Matrix)                :: b !output
         REAL(REALK), ALLOCATABLE :: afull(:,:)
         call time_mat_operations1
                  
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (b%nrow /= a%ncol .or. b%ncol /= a%nrow) then
           CALL LSQUIT( 'wrong dimensions in mat_trans',-1)
         endif
         if (info_memory) write(mat_lu,*) 'Before mat_trans: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_trans(a,b)
         case(mtype_csr)
             call mat_csr_trans(a,b)
         case(mtype_scalapack)
             call mat_scalapack_trans(a,b)
#ifdef VAR_ENABLE_TENSORS
         case(mtype_pdmm)
             call mat_pdmm_trans(a,b)
#endif
         case(mtype_unres_dense)
             call mat_unres_dense_trans(a,b)
         case default
              call lsquit("mat_trans not implemented for this type of matrix",-1)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('TRANS ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_trans: mem_allocated_global =', mem_allocated_global
         call time_mat_operations2(JOB_mat_trans)

      END SUBROUTINE mat_trans

!> \brief Compute the Cholesky decomposition factors of a positive definite matrix
!> \author T. Kjærgaard
!> \date 2012
!> \param a The type(matrix) that should be decomposed
!> \param b The output type(matrix) that contains the cholesky factors.
      SUBROUTINE mat_chol(a, b) 
        implicit none
        TYPE(Matrix),intent(in)     :: a
        TYPE(Matrix),intent(inout)  :: b !output
        real(realk), pointer   :: work1(:)
        real(realk), pointer   :: U_full(:,:)
        integer,pointer    :: IPVT(:)
        real(realk)            :: RCOND, dummy(2), tmstart, tmend
        integer                :: IERR, i, j, fulldim, ndim
        call time_mat_operations1
        if (b%nrow /= a%ncol .or. b%ncol /= a%nrow) then
           CALL LSQUIT( 'wrong dimensions in mat_trans',-1)
        endif
        if (info_memory) write(mat_lu,*) 'Before mat_inv: mem_allocated_global =', mem_allocated_global
        select case(matrix_type)
        case(mtype_unres_dense)
           fulldim = 2*a%nrow
        case(mtype_dense)
           fulldim = a%nrow
        case default
           fulldim = a%nrow
        end select

        select case(matrix_type)
        case(mtype_dense)
           call mat_dense_chol(a,b)
        case(mtype_scalapack)
!           call mat_scalapack_chol(a,b)
           print*,'mat_chol FALLBACK scalapack'
           call mem_alloc(U_full,fulldim,fulldim) 
           call mem_alloc(work1,fulldim)
           call mem_alloc(IPVT,fulldim)
           call mat_scalapack_to_full(A,1.0E0_realk,U_full)
           !Set lower half of U = 0:
           do i = 1, fulldim
              do j = 1, i-1
                 U_full(i,j) = 0.0E0_realk
              enddo
           enddo           
           call dchdc(U_full,fulldim,fulldim,work1,ipvt,0,IERR)           
           call mat_scalapack_set_from_full(U_full,1.0E0_realk,B)
           call mem_dealloc(U_full) 
           call mem_dealloc(work1)
           call mem_dealloc(IPVT)
        case(mtype_csr)
!           call mat_scalapack_chol(a,b)
           print*,'mat_chol FALLBACK csr'
           call mem_alloc(U_full,fulldim,fulldim) 
           call mem_alloc(work1,fulldim)
           call mem_alloc(IPVT,fulldim)
           call mat_csr_to_full(A,1.0E0_realk,U_full)
           !Set lower half of U = 0:
           do i = 1, fulldim
              do j = 1, i-1
                 U_full(i,j) = 0.0E0_realk
              enddo
           enddo           
           call dchdc(U_full,fulldim,fulldim,work1,IPVT,0,IERR)           
           call mat_csr_set_from_full(U_full,1.0E0_realk,B)
           call mem_dealloc(U_full) 
           call mem_dealloc(work1)
           call mem_dealloc(IPVT)
         case(mtype_unres_dense)
!             call mat_unres_dense_inv(a,b)
           print*,'mat_chol FALLBACK unres dense'
           call mem_alloc(U_full,fulldim,fulldim) 
           call mem_alloc(work1,fulldim)
           call mem_alloc(IPVT,fulldim)
           call mat_unres_dense_to_full(A,1.0E0_realk,U_full)
           !Set lower half of U = 0:
           do i = 1, fulldim
              do j = 1, i-1
                 U_full(i,j) = 0.0E0_realk
              enddo
           enddo           
           call dchdc(U_full,fulldim,fulldim,work1,IPVT,0,IERR)           
           call mat_unres_dense_set_from_full(U_full,1.0E0_realk,B)
           call mem_dealloc(U_full) 
           call mem_dealloc(work1)
           call mem_dealloc(IPVT)
#ifdef VAR_ENABLE_TENSORS
        case(mtype_pdmm)
           call mat_pdmm_chol(a, b) 
#endif
        case default
           call lsquit('error mat_chol not implemented',-1)
!!$           call mem_alloc(U_full,fulldim,fulldim) 
!!$           call mem_alloc(work1,fulldim)
!!$           call mem_alloc(IPVT,fulldim)
!!$           call mat_to_full(A,1.0E0_realk,U_full)
!!$           !Set lower half of U = 0:
!!$           do i = 1, fulldim
!!$              do j = 1, i-1
!!$                 U_full(i,j) = 0.0E0_realk
!!$              enddo
!!$           enddo           
!!$           call dchdc(U_full,fulldim,fulldim,work1,0,0,IERR)           
!!$           call mat_set_from_full(U_full,1.0E0_realk,B)
!!$           call mem_dealloc(U_full) 
!!$           call mem_dealloc(work1)
!!$           call mem_dealloc(IPVT)
         end select
         if (info_memory) write(mat_lu,*) 'After mat_trans: mem_allocated_global =', mem_allocated_global
         call time_mat_operations2(JOB_mat_chol)
       END SUBROUTINE mat_chol

!> \brief Compute the Cholesky decomposition factors of a positive definite matrix
!> \author T. Kjærgaard
!> \date 2012
!> \param a The type(matrix) that should be decomposed
      SUBROUTINE mat_dpotrf(a) 
        implicit none
        TYPE(Matrix),intent(inout)     :: a !overwrites matrix A
        integer :: info
!        call time_mat_operations1
        if (info_memory) write(mat_lu,*) 'Before mat_inv: mem_allocated_global =', mem_allocated_global
        select case(matrix_type)
        case(mtype_dense)
           INFO = 0
           CALL DPOTRF('U',A%nrow,A%elms,A%nrow,INFO)
           IF(INFO.NE.0)CALL LSQUIT('DPOTRF ERROR',-1)
        case(mtype_scalapack)
           call mat_scalapack_dpotrf(a)
#ifdef VAR_ENABLE_TENSORS
        case(mtype_pdmm)
           call mat_pdmm_dpotrf(a)
#endif
        case(mtype_unres_dense)
           INFO = 0
           CALL DPOTRF('U',A%nrow,A%elms,A%nrow,INFO)
           IF(INFO.NE.0)CALL LSQUIT('DPOTRF ERROR',-1)
           INFO = 0
           CALL DPOTRF('U',A%nrow,A%elmsb,A%nrow,INFO)
           IF(INFO.NE.0)CALL LSQUIT('DPOTRF ERROR',-1)
        case default
           call lsquit("mat_dpotrf not implemented for this type of matrix",-1)
         end select
         if (info_memory) write(mat_lu,*) 'After mat_trans: mem_allocated_global =', mem_allocated_global
!         call time_mat_operations2(JOB_mat_dpotrf)
       END SUBROUTINE mat_dpotrf

!> \brief Solve the system A*X = B overwriting B with X
!> \author T. Kjærgaard
!> \date 2012
!> \param a The type(matrix) that contain cholesky factors from mat_dpotrf
!> \param b input B output X 
     SUBROUTINE mat_dpotrs(A,B) 
        implicit none
        TYPE(Matrix),intent(in)        :: a
        TYPE(Matrix),intent(inout)     :: b
        integer :: INFO
!        call time_mat_operations1
        if (info_memory) write(mat_lu,*) 'Before mat_inv: mem_allocated_global =', mem_allocated_global
        select case(matrix_type)
        case(mtype_dense)
           CALL DPOTRS('U',A%nrow,B%ncol,A%elms,A%nrow,B%elms,B%nrow,INFO)
        case(mtype_scalapack)
           call mat_scalapack_dpotrs(a,b)
#ifdef VAR_ENABLE_TENSORS
        case(mtype_pdmm)
           call mat_pdmm_dpotrs(a,b)
#endif
        case(mtype_unres_dense)
           CALL DPOTRS('U',A%nrow,B%ncol,A%elms,A%nrow,B%elms,B%nrow,INFO)
           CALL DPOTRS('U',A%nrow,B%ncol,A%elmsb,A%nrow,B%elmsb,B%nrow,INFO)
        case default
           call lsquit("mat_dpotrf not implemented for this type of matrix",-1)
         end select
         if (info_memory) write(mat_lu,*) 'After mat_trans: mem_allocated_global =', mem_allocated_global
!         call time_mat_operations2(JOB_mat_dpotrs)
       END SUBROUTINE mat_dpotrs

!> \brief compute the inverse matrix of real sym positive definite matrix A
!> \author T. Kjærgaard
!> \date 2012
!> \param a The type(matrix) that contain cholesky factors from mat_dpotrf
!> on entry the A matrix is the cholesky factors from mat_dpotrf on exit
!> is the upper triangluar matrix of the symmetric inverse of A 
     SUBROUTINE mat_dpotri(A) 
        implicit none
        TYPE(Matrix),intent(inout)        :: a
        integer :: INFO
!        call time_mat_operations1
        if (info_memory) write(mat_lu,*) 'Before mat_inv: mem_allocated_global =', mem_allocated_global
        select case(matrix_type)
        case(mtype_dense)
           CALL DPOTRI('U',A%nrow,A%elms,A%nrow,INFO)           
           IF(INFO.NE.0)CALL LSQUIT('DPOTRI ERROR',-1)
           !note that depending how you use this you will need to copy the  
!        case(mtype_scalapack)
           !           call mat_scalapack_chol(a,b)
#ifdef VAR_ENABLE_TENSORS
        case(mtype_pdmm)
           call mat_pdmm_dpotri(a)
#endif
        case(mtype_unres_dense)
           CALL DPOTRI('U',A%nrow,A%elms,A%nrow,INFO)           
           IF(INFO.NE.0)CALL LSQUIT('DPOTRI ERROR',-1)
           CALL DPOTRI('U',A%nrow,A%elmsb,A%nrow,INFO)           
           IF(INFO.NE.0)CALL LSQUIT('DPOTRI ERROR',-1)
        case default
           call lsquit("mat_dpotrf not implemented for this type of matrix",-1)
         end select
         if (info_memory) write(mat_lu,*) 'After mat_trans: mem_allocated_global =', mem_allocated_global
!         call time_mat_operations2(JOB_mat_dpotri)
       END SUBROUTINE mat_dpotri

!> \brief creates the inverse matrix of type(matrix).
!> \author T. Kjærgaard
!> \date 2012
!> \param a The type(matrix) that should be inversed
!> \param chol The type(matrix) that contains cholesky factors (from mat_chol)
!> \param c The inverse output type(matrix).
      SUBROUTINE mat_inv(A, A_inv) 
         implicit none
         TYPE(Matrix),intent(in)     :: A
         TYPE(Matrix)                :: A_inv !output
         real(realk), pointer   :: work1(:)
         real(realk), pointer   :: A_inv_full(:,:) 
         integer,pointer    :: IPVT(:)
         real(realk)            :: RCOND, dummy(2), tmstart, tmend
         integer                :: IERR, i, j, fulldim, ndim

         call time_mat_operations1
                  
        select case(matrix_type)
        case(mtype_unres_dense)
           fulldim = 2*a%nrow
        case(mtype_dense)
           fulldim = a%nrow
        case default
           fulldim = a%nrow
        end select

         if (info_memory) write(mat_lu,*) 'Before mat_inv: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_inv(a,a_inv)
#ifdef VAR_ENABLE_TENSORS
         case(mtype_pdmm)
             call mat_pdmm_inv(a,a_inv)
#endif
         case(mtype_scalapack)
!             call mat_scalapack_inv(a,b)
            call mem_alloc(A_inv_full,fulldim,fulldim) 
            call mem_alloc(work1,fulldim)
            call mem_alloc(IPVT,fulldim)
            !Invert U and Ut:
            IPVT = 0 ; RCOND = 0.0E0_realk  
            call mat_scalapack_to_full(A,1.0E0_realk,A_inv_full)
            call DGECO(A_inv_full,fulldim,fulldim,IPVT,RCOND,work1)
            call DGEDI(A_inv_full,fulldim,fulldim,IPVT,dummy,work1,01)
            !Convert framework:
            call mat_scalapack_set_from_full(A_inv_full,1.0E0_realk,A_inv)
            call mem_dealloc(A_inv_full) 
            call mem_dealloc(work1)
            call mem_dealloc(IPVT)
         case(mtype_csr)
!             call mat_csr_inv(a,b)
            call mem_alloc(A_inv_full,fulldim,fulldim) 
            call mem_alloc(work1,fulldim)
            call mem_alloc(IPVT,fulldim)
            !Invert U and Ut:
            IPVT = 0 ; RCOND = 0.0E0_realk  
            call mat_csr_to_full(A,1.0E0_realk,A_inv_full)
            call DGECO(A_inv_full,fulldim,fulldim,IPVT,RCOND,work1)
            call DGEDI(A_inv_full,fulldim,fulldim,IPVT,dummy,work1,01)
            !Convert framework:
            call mat_csr_set_from_full(A_inv_full,1.0E0_realk,A_inv)
            call mem_dealloc(A_inv_full) 
            call mem_dealloc(work1)
            call mem_dealloc(IPVT)
         case(mtype_unres_dense)
             call mat_unres_dense_inv(a,a_inv)
         case default
            call lsquit('mat_inv not implemented',-1)
!!$            call mem_alloc(A_inv_full,fulldim,fulldim) 
!!$            call mem_alloc(work1,fulldim)
!!$            call mem_alloc(IPVT,fulldim)
!!$            !Invert U and Ut:
!!$            IPVT = 0 ; RCOND = 0.0E0_realk  
!!$            call mat_to_full(A,1.0E0_realk,A_inv_full)
!!$            call DGECO(A_inv_full,fulldim,fulldim,IPVT,RCOND,work1)
!!$            call DGEDI(A_inv_full,fulldim,fulldim,IPVT,dummy,work1,01)
!!$            !Convert framework:
!!$            call mat_set_from_full(A_inv_full,1.0E0_realk,A_inv)
!!$            call mem_dealloc(A_inv_full) 
!!$            call mem_dealloc(work1)
!!$            call mem_dealloc(IPVT)
         end select
         if (info_memory) write(mat_lu,*) 'After mat_trans: mem_allocated_global =', mem_allocated_global
         call time_mat_operations2(JOB_mat_inv)

      END SUBROUTINE mat_inv

!> \brief Clone a type(matrix).
!> \author B. Jansik
!> \date 2009
!> \param dest The destination type(matrix) 
!> \param src The source type(matrix) that should be cloned.
!>
!> This makes clone of matrix src to dest, similar to mat_assign, except that it is still in same
!> memory. Matrix src can then be accessed also by the dest name. All of this could be much easier, just dest=src;
!> if not for that bloody '=' operator overload with mat_assign!!
!> 
      SUBROUTINE mat_clone(dest,src)
      implicit none
      type(Matrix) :: src, dest
            dest%ncol=src%ncol; dest%nrow=src%nrow
            dest%elms=>src%elms;! dest%idata=>src%idata
!            dest%permutation => src%permutation
            dest%elmsb=>src%elmsb; !dest%celms=>src%celms
!            dest%celmsb=>src%celmsb
            dest%iaux => src%iaux; dest%raux => src%raux
            dest%complex = src%complex
            dest%val => src%val; dest%col => src%col
            dest%row => src%row; dest%nnz = src%nnz
            dest%PDMID = src%PDMID
#ifdef VAR_SCALAPACK
            dest%localncol=src%localncol; dest%localnrow=src%localnrow
            dest%addr_on_grid => src%addr_on_grid
            dest%p => src%p
#endif
      END SUBROUTINE mat_clone
   
!> \brief Copy a type(matrix).
!> \author L. Thogersen
!> \date 2003
!> \param a The copy output type(matrix)
!> \param b The type(matrix) that should be copied.
      SUBROUTINE mat_assign(a, b)
         implicit none
         TYPE(Matrix), INTENT(INOUT) :: a
         TYPE(Matrix), INTENT(IN)    :: b
         call time_mat_operations1
         
         if (a%nrow /= b%nrow .or. a%ncol /= b%ncol) then
            call lsquit('wrong dimensions in mat_assign',-1)
         endif
         if (info_memory) write(mat_lu,*) 'Before mat_assign: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_assign(a,b)
          case(mtype_csr)
             call mat_csr_assign(a,b)
          case(mtype_scalapack)
             call mat_scalapack_assign(a,b)
#ifdef VAR_ENABLE_TENSORS
          case(mtype_pdmm)
             call mat_pdmm_assign(a,b)
#endif
         case(mtype_unres_dense)
             call mat_unres_dense_assign(a,b)
         case default
              call lsquit("mat_assign not implemented for this type of matrix",-1)
         end select
         if (info_memory) write(mat_lu,*) 'After mat_assign: mem_allocated_global =', mem_allocated_global
         call time_mat_operations2(JOB_mat_assign)

       END SUBROUTINE mat_assign

#ifndef UNITTEST
!> \brief MPI broadcast a type(matrix).
!> \author T. Kjaergaard
!> \date 2010
!> \param a The type(matrix) that should be copied
!> \param slave , true if slave process 
!> \param master integer of master process
      SUBROUTINE mat_mpicopy(a, slave, master)
         implicit none
         TYPE(Matrix), INTENT(INOUT) :: a
         integer(kind=ls_mpik),intent(in) :: master 
         logical,intent(in) :: slave

         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_mpicopy(a,slave, master)
         case(mtype_scalapack)
            call lsquit('mat_mpicopy scalapack error',-1)
#ifdef VAR_ENABLE_TENSORS
         case(mtype_pdmm)
            call lsquit('mat_mpicopy pdmm error',-1)
#endif
          case(mtype_csr)
             call mat_mpicopy_fallback(a,slave, master)
         case(mtype_unres_dense)
             call mat_unres_dense_mpicopy(a,slave, master)
         case default
              call lsquit("mpicopy_typematrix not implemented for this type of matrix",-1)
         end select
       contains
         subroutine mat_mpicopy_fallback(A,slave,master)
           use lsmpi_type
          implicit none
          type(Matrix), intent(inout) :: A
          logical                     :: slave
          integer(kind=ls_mpik)       :: master
          
          real(realk), allocatable :: Afull(:,:)
          integer                  :: i, j
          
          CALL LS_MPI_BUFFER(A%nrow,Master)
          CALL LS_MPI_BUFFER(A%ncol,Master)
          allocate(Afull(A%nrow,A%ncol))
          IF(.NOT.SLAVE)call mat_to_full(A,1.0E0_realk,Afull)
          CALL LS_MPI_BUFFER(Afull,A%nrow,A%ncol,Master)
          
          IF(SLAVE)THEN
             call mat_init(A,A%nrow,A%ncol)
             call mat_set_from_full(Afull,1.0E0_realk,A) 
          ENDIF
          deallocate(Afull)
          
        end subroutine mat_mpicopy_fallback
       END SUBROUTINE mat_mpicopy
#endif

!> \brief Copy and scale a type(matrix).
!> \author L. Thogersen
!> \date 2003
!> \param alpha The scaling parameter
!> \param a The type(matrix) that should be copied
!> \param b The scaled output type(matrix).
      SUBROUTINE mat_copy(alpha,a, b) ! USAGE DISCOURAGED!
         implicit none
         REAL(REALK),  INTENT(IN)    :: alpha
         TYPE(Matrix), INTENT(IN)    :: a
         TYPE(Matrix), INTENT(INOUT) :: b
         call time_mat_operations1
         
         if (b%nrow /= a%nrow .or. b%ncol /= a%ncol) then
           CALL LSQUIT( 'wrong dimensions in mat_copy',-1)
         endif
         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_copy(alpha,a,b)
         case(mtype_csr)
             call mat_csr_copy(alpha,a,b)
         case(mtype_scalapack)
             call mat_scalapack_copy(alpha,a,b)
#ifdef VAR_ENABLE_TENSORS
         case(mtype_pdmm)
             call mat_pdmm_copy(alpha,a,b)
#endif
         case(mtype_unres_dense)
             call mat_unres_dense_copy(alpha,a,b)
         case default
              call lsquit("mat_copy not implemented for this type of matrix",-1)
         end select
         call time_mat_operations2(JOB_mat_copy)

      END SUBROUTINE mat_copy

!> \brief Makes the trace of a square type(matrix).
!> \param a The type(matrix) we want the trace of
!> \return The trace of a
!> \author L. Thogersen
!> \date 2003
      FUNCTION mat_tr(a)
         implicit none
         TYPE(Matrix), intent(IN) :: a
         REAL(realk) :: mat_tr
         call time_mat_operations1

         if (a%nrow /= a%ncol) then
           print *, 'a%nrow, a%ncol =', a%nrow, a%ncol
           CALL LSQUIT( 'Trace is only defined for a square matrix!',-1)
         endif
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'Before mat_tr: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
         case(mtype_dense)
             mat_Tr = mat_dense_Tr(a)
         case(mtype_csr)
            mat_tr = mat_csr_Tr(a)
         case(mtype_scalapack)
            mat_tr = mat_scalapack_Tr(a)
         case(mtype_unres_dense)
             mat_Tr = mat_unres_dense_Tr(a)
         case default
              call lsquit("mat_Tr not implemented for this type of matrix",-1)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('TRACE ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_tr: mem_allocated_global =', mem_allocated_global
         call time_mat_operations2(JOB_mat_tr)

      END FUNCTION mat_tr

!> \brief Make the trace of the product of type(matrix) A and B.
!> \author L. Thogersen
!> \date 2003
!> \param a The first type(matrix) factor
!> \param b The second type(matrix) factor
!> \return Tr(a*b)
      FUNCTION mat_trAB(a,b)
         implicit none
         TYPE(Matrix), intent(IN) :: a,b
         REAL(realk) :: mat_trAB
         call time_mat_operations1

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (a%ncol /= b%nrow .or. a%nrow /= b%ncol) then
           CALL LSQUIT( 'wrong dimensions in mat_trAB',-1)
         endif
         if (info_memory) write(mat_lu,*) 'Before mat_trAB: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
         case(mtype_dense)
             mat_TrAB = mat_dense_TrAB(a,b)
         case(mtype_csr)
             mat_TrAB = mat_csr_TrAB(a,b)
         case(mtype_scalapack)
             mat_TrAB = mat_scalapack_TrAB(a,b)
#ifdef VAR_ENABLE_TENSORS
         case(mtype_pdmm)
             mat_TrAB = mat_pdmm_TrAB(a,b)
#endif
         case(mtype_unres_dense)
             mat_TrAB = mat_unres_dense_TrAB(a,b)
         case default
              call lsquit("mat_TrAB not implemented for this type of matrix",-1)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('TR_AB ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_trAB: mem_allocated_global =', mem_allocated_global
         call time_mat_operations2(JOB_mat_trAB)

      END FUNCTION mat_trAB

!=======================================================================

!> \brief Make c = alpha*ab + beta*c, where a,b,c are type(matrix) and alpha,beta are parameters
!> \author L. Thogersen
!> \date 2003
!> \param a The first type(matrix) factor
!> \param b The second type(matrix) factor
!> \param transa 'T'/'t' if a should be transposed, 'N'/'n' otherwise
!> \param transb 'T'/'t' if b should be transposed, 'N'/'n' otherwise
!> \param alpha The alpha parameter
!> \param beta The beta parameter
!> \param c The output type(matrix)
      SUBROUTINE mat_mul(a, b, transa, transb, alpha, beta, c)
         !c = alpha*ab + beta*c
         !transa = 'T'/'t' - transposed, 'N'/'n' - normal
         implicit none
         TYPE(Matrix), intent(IN) :: a, b
         character, intent(in)    :: transa, transb
         REAL(realk), INTENT(IN)  :: alpha, beta
         TYPE(Matrix), intent(inout):: c
         integer :: ak, bk, ci, cj
         call time_mat_operations1

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         no_of_matmuls = no_of_matmuls + 1
         if (transa == 'n' .or. transa == 'N') then
           ak = a%ncol
           ci = a%nrow
         elseif (transa == 't' .or. transa == 'T') then
           ak = a%nrow
           ci = a%ncol
         endif
         if (transb == 'n' .or. transb == 'N') then
           bk = b%nrow
           cj = b%ncol
         elseif (transb == 't' .or. transb == 'T') then
           bk = b%ncol
           cj = b%nrow
         endif
         if (ak /= bk .or. ci /= c%nrow .or. cj /= c%ncol) then
           print*,'ak',ak,'bk',bk,'ci',ci,'cj',cj
           Call lsquit('wrong dimensions in mat_mul or unknown trans possibility',-1)
         endif
         if (info_memory) write(mat_lu,*) 'Before mat_mul: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
         case(mtype_dense)
            call mat_dense_mul(a,b,transa, transb,alpha,beta,c)
         case(mtype_unres_dense)
            call mat_unres_dense_mul(a,b,transa, transb,alpha,beta,c)
         case(mtype_csr)
            call mat_csr_mul(a,b,transa, transb,alpha,beta,c)
         case(mtype_scalapack)
            call mat_scalapack_mul(a,b,transa, transb,alpha,beta,c)
#ifdef VAR_ENABLE_TENSORS
         case(mtype_pdmm)
            call mat_pdmm_mul(a, b, transa, transb, alpha, beta, c)
#endif
         case default
            call lsquit("mat_mul not implemented for this type of matrix",-1)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('MATMUL ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_mul: mem_allocated_global =', mem_allocated_global
         call time_mat_operations2(JOB_mat_mul)

      END SUBROUTINE mat_mul

!> \brief Make c = alpha*a + beta*b, where a,b are type(matrix) and alpha,beta are parameters
!> \author L. Thogersen
!> \date 2003
!> \param a The first type(matrix) 
!> \param alpha The alpha parameter
!> \param b The second type(matrix) 
!> \param beta The beta parameter
!> \param c The output type(matrix)
      SUBROUTINE mat_add(alpha, a, beta, b, c)
         implicit none
         TYPE(Matrix), intent(IN) :: a, b
         REAL(realk), INTENT(IN)  :: alpha, beta
         TYPE(Matrix)             :: c
         call time_mat_operations1

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (a%nrow /= b%nrow .or. a%ncol /= b%ncol .or. a%nrow /= c%nrow &
            &.or. a%ncol /= c%ncol) then
           print *, 'a%nrow, a%ncol, b%nrow, b%ncol, c%nrow, c%ncol', a%nrow, a%ncol, b%nrow, b%ncol, c%nrow, c%ncol
           CALL LSQUIT( 'wrong dimensions in mat_add',-1)
         endif
         if (info_memory) write(mat_lu,*) 'Before mat_add: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_add(alpha,a,beta,b,c)
         case(mtype_csr)
             call mat_csr_add(alpha,a,beta,b,c)
         case(mtype_scalapack)
             call mat_scalapack_add(alpha,a,beta,b,c)
#ifdef VAR_ENABLE_TENSORS
         case(mtype_pdmm)
             call mat_pdmm_add(alpha,a,beta,b,c)
#endif
         case(mtype_unres_dense)
             call mat_unres_dense_add(alpha,a,beta,b,c)
         case default
              call lsquit("mat_add not implemented for this type of matrix",-1)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('ADD   ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_add: mem_allocated_global =', mem_allocated_global
         call time_mat_operations2(JOB_mat_add)

      END SUBROUTINE mat_add

!> \brief Make Y = alpha*X + Y where X,Y are type(matrix) and a is a parameter
!> \author L. Thogersen
!> \date 2003
!> \param alpha The alpha parameter
!> \param X The input type(matrix) 
!> \param Y The input/output type(matrix) 
      SUBROUTINE mat_daxpy(alpha, X, Y)
         implicit none
         real(realk),intent(in)       :: alpha
         TYPE(Matrix), intent(IN)     :: X
         TYPE(Matrix), intent(INOUT)  :: Y
         call time_mat_operations1

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (x%nrow /= y%nrow .or. x%ncol /= y%ncol) then
            call lsquit('wrong dimensions in mat_daxpy',-1)
         endif
         if (info_memory) write(mat_lu,*) 'Before mat_daxpy: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_daxpy(alpha,x,y)
         case(mtype_csr)
             call mat_csr_daxpy(alpha,x,y)
         case(mtype_scalapack)
             call mat_scalapack_daxpy(alpha,x,y)
#ifdef VAR_ENABLE_TENSORS
         case(mtype_pdmm)
             call mat_pdmm_daxpy(alpha,x,y)
#endif
         case(mtype_unres_dense)
             call mat_unres_dense_daxpy(alpha,x,y)
         case default
              call lsquit("mat_daxpy not implemented for this type of matrix",-1)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('DAXPY ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_daxpy: mem_allocated_global =', mem_allocated_global
         call time_mat_operations2(JOB_mat_daxpy)

      END SUBROUTINE mat_daxpy

!> \brief solves Ax=b where A,x,b are type(matrix) 
!> \author J. Rekkedal
!> \date 2003
!> \param alpha The alpha parameter
!> \param A The input type(matrix) 
!> \param x The input/output type(matrix) 
!> \param b The input type(matrix) 
      SUBROUTINE mat_dposv(A,b,lupri)
         implicit none
         TYPE(Matrix), intent(INOUT)  :: A
         TYPE(Matrix), intent(INOUT)  :: b
         INTEGER,INTENT(IN)           :: lupri

         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_dposv(A,b,lupri)
         case default
              call lsquit("mat_dposv not implemented for this type of matrix",-1)
         end select

      END SUBROUTINE mat_dposv

!> \brief Make the dot product of type(matrix) a and b.
!> \author L. Thogersen
!> \date 2003
!> \param a The first type(matrix) factor
!> \param b The second type(matrix) factor
!> \return The dot product of a and b
      function mat_dotproduct(a,b)
         implicit none
         TYPE(Matrix), intent(IN) :: a,b
         REAL(realk) :: mat_dotproduct
         call time_mat_operations1

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (a%nrow*a%ncol /= b%nrow*b%ncol) then
           CALL LSQUIT( 'wrong dimensions in mat_dotproduct',-1)
         endif
         if (info_memory) write(mat_lu,*) 'Before mat_dotproduct: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
         case(mtype_dense)
            mat_dotproduct = mat_dense_dotproduct(a,b)
         case(mtype_csr)
            mat_dotproduct = mat_csr_dotproduct(a,b)
         case(mtype_scalapack)
            mat_dotproduct = mat_scalapack_dotproduct(a,b)
#ifdef VAR_ENABLE_TENSORS
         case(mtype_pdmm)
            mat_dotproduct = mat_pdmm_dotproduct(a,b)
#endif
         case(mtype_unres_dense)
             mat_dotproduct = mat_unres_dense_dotproduct(a,b)
         case default
              call lsquit("mat_dotproduct not implemented for this type of matrix",-1)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('DOTPRO',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_dotproduct: mem_allocated_global =', mem_allocated_global
         call time_mat_operations2(JOB_mat_dotproduct)

      END FUNCTION mat_dotproduct

!> \brief Make the dot product of type(matrix) a with itself.
!> \author L. Thogersen
!> \date 2003
!> \param a The type(matrix) input
!> \return The dot product of a with itself
      FUNCTION mat_sqnorm2(a)
         implicit none
         TYPE(Matrix), intent(IN) :: a
         REAL(realk) :: mat_sqnorm2
         call time_mat_operations1

         if (info_memory) write(mat_lu,*) 'Before mat_sqnorm2: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
         case(mtype_dense)
            mat_sqnorm2 = mat_dense_sqnorm2(a)
         case(mtype_csr)
            mat_sqnorm2 = mat_csr_sqnorm2(a)
         case(mtype_scalapack)
            mat_sqnorm2 = mat_scalapack_sqnorm2(a)
#ifdef VAR_ENABLE_TENSORS
         case(mtype_pdmm)
            mat_sqnorm2 = mat_pdmm_sqnorm2(a)
#endif
         case(mtype_unres_dense)
             mat_sqnorm2 = mat_unres_dense_sqnorm2(a)
         case default
              call lsquit("mat_sqnorm2 not implemented for this type of matrix",-1)
         end select
         if (info_memory) write(mat_lu,*) 'After mat_sqnorm2: mem_allocated_global =', mem_allocated_global
         call time_mat_operations2(JOB_mat_sqnorm2)

      END FUNCTION mat_sqnorm2

!> \brief Find the absolute largest element of a type(matrix).
!> \author S. Host
!> \date 2005
!> \param a The type(matrix) input
!> \param val The absolute largest element of a
      SUBROUTINE mat_abs_max_elm(a, val) 
         implicit none
         REAL(REALK),  INTENT(OUT)   :: val
         TYPE(Matrix), INTENT(IN)    :: a

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'Before mat_abs_max_elm: mem_allocated_global =', mem_allocated_global
         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_abs_max_elm(a,val)
          case(mtype_csr)
             call mat_csr_abs_max_elm(a,val)
          case(mtype_scalapack)
             call mat_scalapack_abs_max_elm(a,val)
         case(mtype_unres_dense)
             call mat_unres_dense_abs_max_elm(a,val)
         case default
              call lsquit("mat_abs_max_elm not implemented for this type of matrix",-1)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('MAXELM',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_abs_max_elm: mem_allocated_global =', mem_allocated_global
      END SUBROUTINE mat_abs_max_elm

!> \brief Find the largest element of a type(matrix).
!> \author S. Host
!> \date 2005
!> \param a The type(matrix) input
!> \param val The largest element of a
      SUBROUTINE mat_max_elm(a, val, pos) 
         implicit none
         REAL(REALK),  INTENT(OUT)   :: val
         TYPE(Matrix), INTENT(IN)    :: a
         integer, optional           :: pos(2)
         integer                     :: tmp(2)

         if (info_memory) write(mat_lu,*) 'Before mat_max_elm: mem_allocated_global =', mem_allocated_global
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_max_elm(a,val,tmp)
          case(mtype_csr)
             if (present(pos)) call lsquit('mat_max_elm(): position parameter not implemented!',-1)
             call mat_csr_max_elm(a,val)
          case(mtype_scalapack)
             call mat_scalapack_max_elm(a,val,tmp)
#ifdef VAR_ENABLE_TENSORS
          case(mtype_pdmm)
             call mat_pdmm_max_elm(a,val,tmp)
#endif
         case(mtype_unres_dense)
             if (present(pos)) call lsquit('mat_max_elm(): position parameter not implemented!',-1)
             call mat_unres_dense_max_elm(a,val)
         case default
              call lsquit('mat_max_elm not implemented for this type of matrix',-1)
         end select


         if (present(pos)) pos=tmp

         !if (INFO_TIME_MAT) CALL LSTIMER('MAXELM',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_max_elm: mem_allocated_global =', mem_allocated_global
      END SUBROUTINE mat_max_elm

      SUBROUTINE mat_min_elm(a, val, pos) 
         implicit none
         REAL(REALK),  INTENT(OUT)   :: val
         TYPE(Matrix), INTENT(IN)    :: a
         integer, optional           :: pos(2)
         integer                     :: tmp(2)

         if (info_memory) write(mat_lu,*) 'Before mat_min_elm: mem_allocated_global =', mem_allocated_global

         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_min_elm(a,val,tmp)
          case(mtype_scalapack)
             call mat_scalapack_min_elm(a,val,tmp)
#ifdef VAR_ENABLE_TENSORS
          case(mtype_pdmm)
             call mat_pdmm_min_elm(a,val,tmp)
#endif
         case default
              call lsquit("mat_min_elm not implemented for this type of matrix",-1)
         end select


         if (present(pos)) pos=tmp

         if (info_memory) write(mat_lu,*) 'After mat_min_elm: mem_allocated_global =', mem_allocated_global
      END SUBROUTINE mat_min_elm


!> \brief Find the largest element on the diagonal of a type(matrix).
!> \author S. Host
!> \date 2005
!> \param a The type(matrix) input
!> \param pos The position of the diagonal for the largest element of a
!> \param val The largest on element on the diagonal of a
      SUBROUTINE mat_max_diag_elm(a, pos, val) 
         implicit none
         TYPE(Matrix), INTENT(IN)    :: a
         INTEGER, INTENT(OUT)        :: pos
         REAL(REALK),  INTENT(OUT)   :: val

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (a%nrow /= a%ncol) then
           CALL LSQUIT( 'matrix must be symmetric in mat_max_diag_elm',-1)
         endif
         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_max_diag_elm(a,pos,val)
         case(mtype_scalapack)
!             call mat_scalapack_max_diag_elm(a,pos,val)
!             print*,'the maximum element is',val
!             print*,'but the position is more complicated and '
!             print*,'not implemented.'
             call lsquit('mat_max_diag_elm not fully implemented for scalapack type matrix',-1)
         case(mtype_unres_dense)
             call mat_unres_dense_max_diag_elm(a,pos,val)
         case default
              call lsquit("mat_max_diag_elm not implemented for this type of matrix",-1)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('MAXDIA',mat_TSTR,mat_TEN,mat_lu)
      END SUBROUTINE mat_max_diag_elm

!> \brief Squares all the off-diagonal elements in a type(matrix), and returns the sum
!> \author L. Thogersen
!> \date 2003
!> \param a The type(matrix) input
!> \return Sum of the squares of off-diagonal elements in a
      FUNCTION mat_outdia_sqnorm2(a)
         implicit none
         TYPE(Matrix), intent(IN) :: a
         REAL(realk) :: mat_outdia_sqnorm2
         select case(matrix_type)
         case(mtype_dense)
             mat_outdia_sqnorm2 = mat_dense_outdia_sqnorm2(a)
         case(mtype_csr)
             mat_outdia_sqnorm2 = mat_csr_outdia_sqnorm2(a)
         case(mtype_scalapack)
             mat_outdia_sqnorm2 = mat_scalapack_outdia_sqnorm2(a)
         case(mtype_unres_dense)
             mat_outdia_sqnorm2 = mat_unres_dense_outdia_sqnorm2(a)
         case default
              call lsquit("mat_outdia_sqnorm2 not implemented for this type of matrix",-1)
         end select
!         print *, "outdia got ", mat_outdia_sqnorm2
      END FUNCTION mat_outdia_sqnorm2

!> \brief General diagonalization F*C = S*C*e
!> \author L. Thogersen
!> \date 2003
!> \param F Fock/Kohn-Sham matrix
!> \param S Overlap matrix
!> \param eival Eigenvalues
!> \param Cmo C coefficients
      SUBROUTINE mat_diag_f(F,S,eival,Cmo)
         !solves FC = SCe 
         implicit none
         TYPE(Matrix), intent(IN) :: F,S
         type(matrix)             :: Cmo  !output
         real(realk),intent(INOUT)  :: eival(:)
!
         TYPE(MATRIX)             :: A,B 
         real(realk), allocatable :: tmp(:), eval(:), cmod(:), wrk(:)
         integer                  :: ndim

         select case(matrix_type)
         case(mtype_dense)
            call time_mat_operations1
            call mat_dense_diag_f(F,S,eival,Cmo)
            call time_mat_operations2(JOB_mat_diag_f)
#ifdef VAR_ENABLE_TENSORS
         case(mtype_pdmm)
            call time_mat_operations1
            call mat_pdmm_diag_f(F,S,eival,Cmo)
            call time_mat_operations2(JOB_mat_diag_f)
#endif
         case(mtype_scalapack)
            call mat_init(A,F%nrow,F%ncol)
            call mat_init(B,S%nrow,S%ncol)
            call mat_copy(1E0_realk,F,A)
            call mat_copy(1E0_realk,S,B)
            call time_mat_operations1
            call mat_scalapack_diag_f(A,B,eival,Cmo)
            call time_mat_operations2(JOB_mat_diag_f)
            call mat_free(A)
            call mat_free(B)
#ifndef UNITTEST
#endif
         case(mtype_unres_dense)
             call mat_unres_dense_diag_f(F,S,eival,Cmo)
         case default
            print *, "FALLBACK diag_f...", S%nrow
            ndim = s%nrow
            ALLOCATE(tmp(Ndim*Ndim),eval(Ndim),cmod(Ndim*Ndim))
            call mat_to_full(S, 1E0_realk, tmp)
            call mat_to_full(F, 1E0_realk, cmod)
            call time_mat_operations1
            call my_DSYGV(ndim,cmod,tmp,eival,"mat_diag_f          ")
            call time_mat_operations2(JOB_mat_diag_f)
            call mat_set_from_full(cmod, 1E0_realk, Cmo)
            DEALLOCATE(tmp,eval,cmod)
         end select
     END SUBROUTINE mat_diag_f

!> \brief computes all eigenvalues and, eigenvectors of a real symmetric matrix S.
!> \author T. Kjaergaard
!> \date 2012
!> \param S matrix (input and output, output is the eigenvectors)
!> \param eival Eigenvalues (output vector)
      SUBROUTINE mat_dsyev(S,eival,ndim)
         implicit none
         TYPE(Matrix), intent(INOUT) :: S
         real(realk), intent(INOUT) :: eival(ndim)
!
         type(matrix) :: B
         real(realk),pointer :: work(:),S_full(:,:)
         integer :: infdiag,ndim,lwork
         infdiag=0

         select case(matrix_type)
         case(mtype_dense)
            call time_mat_operations1
            call mat_dense_dsyev(S,eival,ndim)
            call time_mat_operations2(JOB_mat_dsyev)
         case(mtype_scalapack)
            call mat_init(B,S%nrow,S%ncol)
            call time_mat_operations1
            call mat_scalapack_dsyev(S,B,eival,ndim)
            call time_mat_operations2(JOB_mat_dsyev)
            call mat_assign(S,B)
            call mat_free(B)
#ifdef VAR_ENABLE_TENSORS
         case(mtype_pdmm)
            call time_mat_operations1
            call mat_pdmm_dsyev(S,eival,ndim)
            call time_mat_operations2(JOB_mat_dsyev)
#endif
         case(mtype_csr)
            call mem_alloc(S_full,ndim,ndim)
            call mat_csr_to_full(S,1.0E0_realk,S_full)
            !============================================================
            ! we inquire the size of lwork
            lwork = -1
            call mem_alloc(work,5)
            call time_mat_operations1
            call dsyev('V','U',ndim,S_FULL,ndim,eival,work,lwork,infdiag)
            lwork = NINT(work(1))
            call mem_dealloc(work)
            !=============================================================     
            call mem_alloc(work,lwork)
            !diagonalization
            call dsyev('V','U',ndim,S_FULL,ndim,eival,work,lwork,infdiag)
            call time_mat_operations2(JOB_mat_dsyev)
            call mem_dealloc(work)
            call mat_csr_set_from_full(S_FULL, 1E0_realk,S)
            call mem_dealloc(S_full)
            if(infdiag.ne. 0) then
               print*,'lowdin_diag: dsyev failed, info=',infdiag
               call lsquit('lowdin_diag: diagonalization failed.',-1)
            end if
         case(mtype_unres_dense)
            call mat_unres_dense_dsyev(S,eival,ndim)
         case default
            call lsquit('mat_dsyev not implemented for this type',-1)
!            call mem_alloc(S_full,ndim,ndim)
!            call mat_to_full(S,1.0E0_realk,S_full)
!            !============================================================
!            ! we inquire the size of lwork
!            lwork = -1
!            call mem_alloc(work,5)
!            call time_mat_operations1
!            call dsyev('V','U',ndim,S_FULL,ndim,eival,work,lwork,infdiag)
!            lwork = NINT(work(1))
!            call mem_dealloc(work)
!            !=============================================================     
!            call mem_alloc(work,lwork)
!            !diagonalization
!            call dsyev('V','U',ndim,S_FULL,ndim,eival,work,lwork,infdiag)
!            call time_mat_operations2(JOB_mat_dsyev)
!            call mem_dealloc(work)
!            call mat_set_from_full(S_FULL, 1E0_realk,S)
!            call mem_dealloc(S_full)
!            if(infdiag.ne. 0) then
!               print*,'lowdin_diag: dsyev failed, info=',infdiag
!               call lsquit('lowdin_diag: diagonalization failed.',-1)
!            end if
         end select
       END SUBROUTINE mat_dsyev
 
!> \brief computes one specific eigenvalue of a real symmetric matrix S.
!> \author S. Reine
!> \date May 2012
!> \param S matrix (input and output, output is the eigenvectors)
!> \param eival Eigenvalue
!> \param ieig Index of the eigenvalue to be returned (in acending order)
  SUBROUTINE mat_dsyevx(S,eival,ieig)
    implicit none
    TYPE(Matrix), intent(INOUT) :: S
    real(realk), intent(INOUT)  :: eival
    integer,intent(in)          :: ieig
!
    integer             :: ndim
    real(realk),pointer :: Sfull(:,:)
!
    select case(matrix_type)
    case(mtype_dense)
       call time_mat_operations1           
       CALL mat_dense_dsyevx(S,eival,ieig)
       call time_mat_operations2(JOB_mat_dsyevx)
    case(mtype_scalapack)
       call time_mat_operations1           
       CALL mat_scalapack_dsyevx(S,eival,ieig)
       call time_mat_operations2(JOB_mat_dsyevx)
#ifdef VAR_ENABLE_TENSORS
    case(mtype_pdmm)
       call time_mat_operations1           
       CALL mat_pdmm_dsyevx(S,eival,ieig)
       call time_mat_operations2(JOB_mat_dsyevx)
#endif
    case(mtype_unres_dense)
       call lsquit('mat_dsyevx mtype_unres_dense not implemented',-1)
       call time_mat_operations1           
       ndim = S%nrow
       call mem_alloc(Sfull,ndim,ndim)
       call mat_unres_dense_to_full(S,1.0E0_realk,Sfull)
       CALL mat_dense_dsyevx_aux(Sfull,eival,ndim,ieig)
       call mem_dealloc(Sfull)
       call time_mat_operations2(JOB_mat_dsyevx)
    case(mtype_csr)
       ndim = S%nrow
       call mem_alloc(Sfull,ndim,ndim)
       call mat_csr_to_full(S,1.0E0_realk,Sfull)
       call time_mat_operations1           
       CALL mat_dense_dsyevx_aux(Sfull,eival,ndim,ieig)
       call time_mat_operations2(JOB_mat_dsyevx)
       call mem_dealloc(Sfull)
    case default
       call lsquit('mat_dsyevx not implemented',-1)
!      write(*,*) 'FALLBACK: mat_dsyevx'
!      ndim = S%nrow
!      call mem_alloc(Sfull,ndim,ndim)
!      call mat_to_full(S,1.0E0_realk,Sfull)
!      call time_mat_operations1           
!      CALL mat_dense_dsyevx_aux(Sfull,eival,ndim,ieig)
!      call time_mat_operations2(JOB_mat_dsyevx)
!      call mem_dealloc(Sfull)
    end select
  END SUBROUTINE mat_dsyevx

!> \brief Returns a section of a matrix
!> \author L. Thogersen
!> \date 2003
!> \param A Input type(matrix)
!> \param from_row Begin at this row
!> \param to_row End at this row
!> \param from_col Begin at this column
!> \param to_col End at this column
!> \param Asec The section of the type(matrix)
    subroutine mat_section(A,from_row,to_row,from_col,to_col,Asec)
      implicit none
      type(Matrix), intent(in) :: A
      integer, intent(in) :: from_row, to_row, from_col, to_col
      type(Matrix), intent(inout) :: Asec  !output
      type(Matrix) :: B,Bsec
      real(realk),pointer :: Afull(:,:),Bsecfull(:,:)
      
      !Check if Asec is inside A
      if (to_row > A%nrow .or. from_row < 1 .or. from_col < 1 .or. A%ncol < to_col) then
         call lsquit('Asec not inside A in mat_section',-1)
      endif
      !Check if the section size is positive
      if (from_row > to_row .or. from_col > to_col) then
         call lsquit('from_row or from_col > to_row or to_col in mat_section',-1)
      endif
      !Check if allocated space for section is the right size
      if (Asec%nrow /= to_row - from_row + 1 .or.&
           & Asec%ncol /= to_col - from_col + 1) then
         CALL LSQUIT( 'Wrong dimensions in mat_section',-1)
      endif
      
      select case(matrix_type)
      case(mtype_dense)
         call mat_dense_section(A,from_row,to_row,from_col,to_col,Asec)
      case(mtype_unres_dense)
         call mat_unres_dense_section(A,from_row,to_row,from_col,to_col,Asec)
      case(mtype_scalapack)
#ifdef VAR_SCALAPACK
         write(*,'(A)') 'Fallback mat_section for mtype_scalapack'
         call mem_alloc(Afull,A%nrow,A%ncol)
         call mat_dense_init(B,A%nrow,A%ncol)
         !transform from type to dense type
         call mat_to_full(A,1E0_realk,Afull)
         call mat_dense_set_from_full(Afull,1E0_realk,B)
         call mem_dealloc(Afull)
         !build section from dense mat
         call mat_dense_init(Bsec,Asec%nrow,Asec%ncol)
         call mat_dense_section(B,from_row,to_row,from_col,to_col,Bsec)
         call mat_dense_free(B)
         !build full from dense section mat
         call mem_alloc(Bsecfull,Asec%nrow,Asec%ncol)
         call mat_dense_to_full(Bsec, 1E0_realk, Bsecfull)
         call mat_dense_free(Bsec)
         call mat_set_from_full(Bsecfull,1E0_realk,Asec)
         call mem_dealloc(Bsecfull)
         write(*,'(A)') 'debug: fallback mat_section finished'
#else
         call lsquit('matrix type scalapack requires VAR_SCALAPACK',-1)
#endif
      case default
         call lsquit("mat_section not implemented for this type of matrix",-1)
      end select
    END SUBROUTINE mat_section
   
!> \brief Inserts (overwrites) a matrix-section into at matrix (like mat_section but opposite)
!> \author C. Nygaard
!> \date June 11 2012
!> \param Asec The section to be inserted
!> \param from_row Begin at this row
!> \param to_row End at this row
!> \param from_col Begin at this column
!> \param to_col End at this column
!> \param A The matrix into which the section is inserted
subroutine mat_insert_section (Asec, from_row, to_row, from_col, to_col, A)

implicit none

type(matrix), intent(in)    :: Asec
integer, intent(in)         :: from_row, to_row, from_col, to_col
type(matrix), intent(inout) :: A

!Check if Asec is inside A
if (to_row > A%nrow .or. from_row < 1 .or. from_col < 1 .or. A%ncol < to_col) then
   CALL LSQUIT( 'Asec not inside A in mat_insert_section',-1)
endif
!Check if the section size is positive
if (from_row > to_row .or. from_col > to_col) then
   CALL LSQUIT( 'from_row or from_col > to_row or to_col in mat_insert_section',-1)
endif
!Check if allocated space for section is the right size
if (Asec%nrow /= to_row - from_row + 1 .or.&
     & Asec%ncol /= to_col - from_col + 1) then
   CALL LSQUIT( 'Wrong dimensions in mat_insert_section',-1)
endif

select case(matrix_type)
case (mtype_dense)
  call mat_dense_insert_section (Asec, from_row, to_row, from_col, to_col, A)
case(mtype_unres_dense)
  call mat_unres_dense_insert_section (Asec, from_row, to_row, from_col, to_col, A)
case default
  call lsquit("mat_insert_section not implemented for this type of matrix",-1)
end select

end subroutine mat_insert_section
 
!> \brief Set a type(matrix) to identity, i.e. I(i,j) = 1 for i = j, 0 otherwise
!> \author L. Thogersen
!> \date 2003
!> \param I Matrix to be set equal to identity
      subroutine mat_identity(I)
         implicit none
         type(Matrix), intent(inout) :: I
         real(realk), ALLOCATABLE    :: ifull(:,:)
         integer                     :: j
         !
         type(Matrix) :: TMP

         if (info_memory) write(mat_lu,*) 'Before mat_identity: mem_allocated_global =', mem_allocated_global
         !print *, "mat_identity inefficient, use mat_add_identity instead!"
         if (I%nrow /= I%ncol) then
           CALL LSQUIT( 'cannot make identity matrix with different ncol and nrow',-1)
         endif
         select case(matrix_type)
         case(mtype_dense)
             call time_mat_operations1
             call mat_dense_identity(I)
             call time_mat_operations2(JOB_mat_identity)
         case(mtype_csr)
             call time_mat_operations1
             call mat_csr_identity(I)
             call time_mat_operations2(JOB_mat_identity)
         case(mtype_scalapack)
             call time_mat_operations1
             call mat_scalapack_add_identity(1E0_realk,0E0_realk,TMP,I)
             call time_mat_operations2(JOB_mat_identity)
#ifdef VAR_ENABLE_TENSORS
         case(mtype_pdmm)
             call time_mat_operations1
             call mat_pdmm_identity(I)
             call time_mat_operations2(JOB_mat_identity)
#endif
         case(mtype_unres_dense)
             call time_mat_operations1
             call mat_unres_dense_identity(I)
             call time_mat_operations2(JOB_mat_identity)
         case default
            call time_mat_operations1
            print *, "FALLBACK: mat_identity"
            allocate(ifull(I%nrow, I%ncol))
            DO j = 1, I%ncol
               ifull(1:j-1,j) = 0E0_realk
               ifull(j,j)     = 1E0_realk
               ifull(j+1:I%nrow,j) = 0E0_realk
            END DO
            call time_mat_operations2(JOB_mat_identity)
            call mat_set_from_full(ifull,1E0_realk,I)
            deallocate(ifull)
         end select
         if (info_memory) write(mat_lu,*) 'After mat_identity: mem_allocated_global =', mem_allocated_global
      END SUBROUTINE mat_identity

!> \brief Add identity to a type(matrix), i.e. C = alpha*I + beta*B  >>> NOTE: ALLOCATES A MATRIX! <<<
!> \author L. Thogersen
!> \date 2003
!> \param alpha Alpha parameter
!> \param beta Beta parameter
!> \param B Input matrix B
!> \param C Output matrix C
      SUBROUTINE mat_add_identity(alpha, beta, B, C)
         implicit none
         TYPE(Matrix), intent(IN) :: B
         REAL(realk), INTENT(IN)  :: alpha, beta
         TYPE(Matrix)             :: C
         type(matrix)             :: I

         if (info_memory) write(mat_lu,*) 'Before mat_add_identity: mem_allocated_global =', mem_allocated_global
         if (b%nrow /= c%nrow .or. b%ncol /= c%ncol) then
           CALL LSQUIT( 'wrong dimensions in mat_add_identity',-1)
         endif
         select case(matrix_type)
         case(mtype_dense)
            call mat_init(I, b%nrow, b%ncol)
            call mat_dense_identity(I)
            call mat_dense_add(alpha,I,beta,b,c)
            call mat_free(I)
         case(mtype_csr)
            call mat_csr_add_identity(alpha, beta, B, C)
         case(mtype_scalapack)
            call mat_scalapack_add_identity(alpha, beta, B, C)
#ifdef VAR_ENABLE_TENSORS
         case(mtype_pdmm)
            call mat_pdmm_add_identity(alpha, beta, B, C)
#endif
         case(mtype_unres_dense)
            call mat_init(I, b%nrow, b%ncol)
            call mat_unres_dense_identity(I)
            call mat_unres_dense_add(alpha,I,beta,b,c)
            call mat_free(I)
         case default
              call lsquit("mat_add_identity not implemented for this type of matrix",-1)
         end select
         if (info_memory) write(mat_lu,*) 'After mat_add_identity: mem_allocated_global =', mem_allocated_global
      END SUBROUTINE mat_add_identity

!> \brief Create or overwrite block in a type(matrix)
!> \author S. Host
!> \date 2009
!> \param A Input/output matrix where we want to create a block
!> \param fullmat Standard fortran matrix containing the block to put into A
!> \param fullrow Number of rows in fullmat
!> \param fullcol Number of columns in fullmat
!> \param insertrow Insert block in A beginning at this row
!> \param insertcol Insert block in A beginning at this col
      subroutine mat_create_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         implicit none
         integer, intent(in) :: fullrow,fullcol,insertrow,insertcol
         real(Realk), intent(in) :: fullmat(fullrow,fullcol)
         type(Matrix), intent(inout) :: A
         call time_mat_operations1

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (insertrow+fullrow > A%nrow+1 .or. &
          &  insertcol+fullcol > A%ncol+1 .or. fullrow < 1 .or. fullcol < 1 .or. &
          & insertrow < 1 .or. insertcol < 1) then
           WRITE(mat_lu,*) 'Cannot create block, the indexes', &
           & fullrow,fullcol,insertrow,insertcol, &
           & 'are out of the bounds - nrow, ncol =',A%nrow,A%ncol
           CALL lsQUIT('Cannot create block (subroutine mat_create_block)',mat_lu)
         endif
         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_create_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         case(mtype_unres_dense)
             call mat_unres_dense_create_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         case(mtype_scalapack)
              call mat_scalapack_create_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
#ifdef VAR_ENABLE_TENSORS
         case(mtype_pdmm)
              call mat_pdmm_create_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
#endif
         case default
              call lsquit("mat_create_block not implemented for this type of matrix",-1)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('CREATE',mat_TSTR,mat_TEN,mat_lu)

         call time_mat_operations2(JOB_mat_create_block)
      END SUBROUTINE mat_create_block

!> \brief Add block to type(matrix) - add to existing elements, don't overwrite
!> \author T. Kjaergaard
!> \date 2009
!> \param A Input/output matrix where we want to add a block
!> \param fullmat Standard fortran matrix containing the block to add to A
!> \param fullrow Number of rows in fullmat
!> \param fullcol Number of columns in fullmat
!> \param insertrow Add block to A beginning at this row
!> \param insertcol Add block to A beginning at this col
      subroutine mat_add_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         implicit none
         integer, intent(in) :: fullrow,fullcol,insertrow,insertcol
         real(Realk), intent(in) :: fullmat(fullrow,fullcol)
         type(Matrix), intent(inout) :: A
         call time_mat_operations1

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (insertrow+fullrow > A%nrow+1 .or. &
          &  insertcol+fullcol > A%ncol+1 .or. fullrow < 1 .or. fullcol < 1 .or. &
          & insertrow < 1 .or. insertcol < 1) then
           WRITE(mat_lu,*) 'Cannot add block, the indexes', &
           & fullrow,fullcol,insertrow,insertcol, &
           & 'are out of the bounds - nrow, ncol =',A%nrow,A%ncol
           CALL lsQUIT('Cannot add block (subroutine mat_add_block)',mat_lu)
         endif
         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_add_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         case(mtype_scalapack)
            call mat_scalapack_add_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         case(mtype_unres_dense)
             call mat_unres_dense_add_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         case default
              call lsquit("mat_add_block not implemented for this type of matrix",-1)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('CREATE',mat_TSTR,mat_TEN,mat_lu)

         call time_mat_operations2(JOB_mat_add_block)
     END SUBROUTINE mat_add_block

!> \brief Retrieve block from type(matrix) 
!> \author T. Kjaergaard
!> \date 2009
!> \param A Input matrix from which we want to retrive a block
!> \param fullmat Return the desired block in this standard fortran matrix
!> \param fullrow Number of rows in fullmat
!> \param fullcol Number of columns in fullmat
!> \param insertrow Retrive block from A beginning at this row
!> \param insertcol Retrive block from A beginning at this col
      subroutine mat_retrieve_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         implicit none
         integer, intent(in) :: fullrow,fullcol,insertrow,insertcol
         real(Realk), intent(inout) :: fullmat(fullrow,fullcol)
         type(Matrix), intent(inout) :: A
         call time_mat_operations1

         if (info_memory) write(mat_lu,*) 'Before mat_retrieve_block: mem_allocated_global =', mem_allocated_global
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (insertrow+fullrow > A%nrow+1 .or. &
          &  insertcol+fullcol > A%ncol+1 .or. fullrow < 1 .or. fullcol < 1 .or. &
          & insertrow < 1 .or. insertcol < 1) then
           WRITE(mat_lu,*) 'Cannot retrieve block, the indexes', &
           & fullrow,fullcol,insertrow,insertcol, &
           & 'are out of the bounds - nrow, ncol =',A%nrow,A%ncol
           CALL lsQUIT('Cannot retrieve block (subroutine mat_retrieve_block)',mat_lu)
         endif
         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_retrieve_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         case(mtype_csr)
             call mat_csr_retrieve_block_full(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         case(mtype_scalapack)
            call mat_scalapack_retrieve_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
#ifdef VAR_ENABLE_TENSORS
         case(mtype_pdmm)
            call mat_pdmm_retrieve_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
#endif
         case(mtype_unres_dense)
             call mat_unres_dense_retrieve_block(A,fullmat,fullrow,fullcol,insertrow,insertcol)
         case default
              call lsquit("mat_retrieve_block not implemented for this type of matrix",-1)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('CREATE',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_retrieve_block: mem_allocated_global =', mem_allocated_global
         call time_mat_operations2(JOB_mat_retrieve_block)

     END SUBROUTINE mat_retrieve_block

!> \brief Scale a type(matrix) A by a scalar alpha
!> \author L. Thogersen
!> \date 2003
!> \param alpha Scaling parameter
!> \param A Input/output matrix which we want to scale
      subroutine mat_scal(alpha,A)
         implicit none
         real(realk), intent(in) :: alpha
         type(Matrix), intent(inout) :: A
         call time_mat_operations1

         if (info_memory) write(mat_lu,*) 'Before mat_scal: mem_allocated_global =', mem_allocated_global
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_scal(alpha,A)
         case(mtype_csr)
             call mat_csr_scal(alpha, A)
         case(mtype_scalapack)
             call mat_scalapack_scal(alpha, A)
#ifdef VAR_ENABLE_TENSORS
         case(mtype_pdmm)
             call mat_pdmm_scal(alpha, A)
#endif
         case(mtype_unres_dense)
             call mat_unres_dense_scal(alpha,A)
         case default
              call lsquit("mat_scal not implemented for this type of matrix",-1)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('SCAL  ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_scal: mem_allocated_global =', mem_allocated_global
         call time_mat_operations2(JOB_mat_scal)

      end subroutine mat_scal

!> \brief Scale the diagonal of a type(matrix) A by a scalar alpha
!> \author L. Thogersen
!> \date 2003
!> \param alpha Scaling parameter
!> \param A Input/output matrix which we want to scale
      subroutine mat_scal_dia(alpha,A)
         implicit none
         real(realk), intent(in) :: alpha
         type(Matrix), intent(inout) :: A
         real(realk), allocatable :: afull(:,:)
         integer i
         call time_mat_operations1

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (A%nrow /= A%ncol) then
           CALL LSQUIT( 'cannot scale diagonal since ncol /= nrow',-1)
         endif

         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_scal_dia(alpha,A)
         case(mtype_unres_dense)
             call mat_unres_dense_scal_dia(alpha,A)
         case(mtype_scalapack)
            call mat_scalapack_scal_dia(alpha,A)
         case default
            print *, "FALLBACK scale_dia"
            allocate(afull(a%nrow, a%ncol))
            call mat_to_full(a,1E0_realk,afull)
            do i = 1,A%nrow
               afull(i,i) = afull(i,i) * alpha
            enddo
            call mat_set_from_full(afull, 1E0_realk, a)
            DEALLOCATE(afull)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('SCADIA',mat_TSTR,mat_TEN,mat_lu)
         call time_mat_operations2(JOB_mat_scal_dia)
      end subroutine mat_scal_dia

!> \brief Scale the diagonal with a vector
!> \author T. Kjaergaard
!> \date 2012
!> \param alpha Scaling parameter
!> \param A Input/output matrix which we want to scale
      subroutine mat_scal_dia_vec(alpha,A,ndim)
         implicit none
         integer, intent(in)     :: ndim
         real(realk), intent(in) :: alpha(ndim)
         type(Matrix), intent(inout) :: A
         real(realk), allocatable :: afull(:,:)
         integer i
         call time_mat_operations1

         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         if (A%nrow /= A%ncol) then
           CALL LSQUIT( 'cannot scale diagonal since ncol /= nrow',-1)
         endif

         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_scal_dia_vec(alpha,A,ndim)
         case(mtype_scalapack)
            call mat_scalapack_scal_dia_vec(alpha,A,ndim)
         case default
            print *, "FALLBACK scale_dia"
            allocate(afull(a%nrow, a%ncol))
            call mat_to_full(a,1E0_realk,afull)
            do i = 1,A%nrow
               afull(i,i) = afull(i,i) * alpha(i)
            enddo
            call mat_set_from_full(afull, 1E0_realk, a)
            DEALLOCATE(afull)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('SCADIA',mat_TSTR,mat_TEN,mat_lu)
         call time_mat_operations2(JOB_mat_scal_dia_vec)
      end subroutine mat_scal_dia_vec

!> \brief Set a type(matrix) A to zero
!> \author L. Thogersen
!> \date 2003
!> \param A Input/output matrix which should be set to zero
      subroutine mat_zero(A)
         implicit none
         type(Matrix), intent(inout) :: A
         call time_mat_operations1

         if (info_memory) write(mat_lu,*) 'Before mat_zero: mem_allocated_global =', mem_allocated_global
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         select case(matrix_type)
         case(mtype_dense)
             call mat_dense_zero(A)
         case(mtype_scalapack)
             call mat_scalapack_zero(A)
#ifdef VAR_ENABLE_TENSORS
         case(mtype_pdmm)
             call mat_pdmm_zero(A)
#endif
         case(mtype_unres_dense)
             call mat_unres_dense_zero(A)
         case(mtype_csr)
             call mat_csr_zero(A)
         case default
              call lsquit("mat_zero not implemented for this type of matrix",-1)
         end select
         !if (INFO_TIME_MAT) CALL LSTIMER('ZERO  ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_zero: mem_allocated_global =', mem_allocated_global
         call time_mat_operations2(JOB_mat_zero)

      end subroutine mat_zero

!> \brief Set the lower triangular part of type(matrix) A to zero
!> \author T. Kjaergaard
!> \date 2012
!> \param A Input/output matrix which should be set to zero
subroutine mat_setlowertriangular_zero(A)
  implicit none
  type(Matrix), intent(inout) :: A
!
  real(realk),pointer :: Afull(:,:)
  call time_mat_operations1
  
  if (info_memory) write(mat_lu,*) 'Before mat_zero: mem_allocated_global =',&
       &mem_allocated_global
  select case(matrix_type)
  case(mtype_dense)
     call set_lowertriangular_zero(A%elms,A%nrow,A%ncol)
  case(mtype_unres_dense)
     call set_lowertriangular_zero(A%elms,A%nrow,A%ncol)
     call set_lowertriangular_zero(A%elmsb,A%nrow,A%ncol)
  case(mtype_scalapack)
     call mat_scalapack_setlowertriangular_zero(A)
  case(mtype_csr)
     allocate(afull(a%nrow, a%ncol))
     call mat_csr_to_full(a,1E0_realk,afull)     
     call set_lowertriangular_zero(Afull,A%nrow,A%ncol)
     call mat_csr_set_from_full(afull, 1E0_realk, a)
     deallocate(afull)
  case default
     call lsquit('Fallback mat_setlowertriangular_zero not impl',-1)
!     print*,'Fallback mat_setlowertriangular_zero'
!     allocate(afull(a%nrow, a%ncol))
!     call mat_to_full(a,1E0_realk,afull)     
!     call set_lowertriangular_zero(Afull,A%nrow,A%ncol)
!     deallocate(afull)
  end select
  if (info_memory) write(mat_lu,*) 'After mat_zero: mem_allocated_global =',&
       & mem_allocated_global
  call time_mat_operations2(JOB_mat_setlowertriangular_zero)  
end subroutine mat_setlowertriangular_zero

subroutine set_lowertriangular_zero(elms,dimenA,dimenB)
  implicit none
  integer,intent(in) :: dimenA,dimenB
  real(realk) :: elms(dimenA,dimenB)
  !
  integer :: A,B
  DO B=1,dimenB
     DO A=B+1,dimenA
        elms(A,B) = 0.0E0_realk
     ENDDO
  ENDDO
end subroutine set_lowertriangular_zero

!> \brief Write a type(matrix) to disk.
!> \author L. Thogersen
!> \date 2003
!> \param iunit Logical unit number of file which matrix should be written to
!> \param A Matrix which should be written on disk
      subroutine mat_write_to_disk(iunit,A,OnMaster)
         implicit none
         integer, intent(in) :: iunit
         type(Matrix), intent(in) :: A
         logical,optional :: OnMaster !obsolete
         !
         integer(kind=long) :: ncol,nrow 
         real(realk), allocatable :: afull(:,:)
         call time_mat_operations1

         if (info_memory) write(mat_lu,*) 'Before mat_write_to_disk: mem_allocated_global =', mem_allocated_global
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         select case(matrix_type)
         case(mtype_dense)
            call mat_dense_write_to_disk(iunit,A)
         case(mtype_csr)
            call mat_csr_write_to_disk(iunit,A)
         case(mtype_scalapack)
            !The master collects the info and write to disk
            allocate(afull(a%nrow, a%ncol))
            call mat_to_full(a,1E0_realk,afull)
            nrow = A%Nrow
            ncol = A%Ncol
            write(iunit) Nrow, Ncol
            write(iunit) afull
            deallocate(afull)
         case(mtype_unres_dense)
            call mat_unres_dense_write_to_disk(iunit,A)
         case default
            print *, "FALLBACK: mat_write_to_disk"
            allocate(afull(a%nrow, a%ncol))
            call mat_to_full(a,1E0_realk,afull)
            nrow = A%Nrow
            ncol = A%Ncol
            write(iunit) Nrow, Ncol
            write(iunit) afull
            deallocate(afull)
         end select

         !if (INFO_TIME_MAT) CALL LSTIMER('WRITE ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_write_to_disk: mem_allocated_global =', mem_allocated_global
         call time_mat_operations2(JOB_mat_write_to_disk)

      end subroutine mat_write_to_disk


!> \brief Add some logical auxiliary information to a type(matrix) on disk.
!> \author S. Host
!> \date June 2010
!> \param iunit Logical unit number of file containing the matrix
!> \param info Info to be written
!>
!> This is needed because the dens.restart that is dumped after calculatation
!> has ended can be either in the standard AO basis or the grand-canonical
!> basis. A dens.restart obtained from in the AO basis will not work in the
!> grand-canonical basis and vice versa. For the moment, it is only possible
!> to put true or false here, but that could easily be changed to a character
!> string saying e.g. 'AOBASIS', 'GCBASIS', or whatever other basis you could
!> think of. Should be independent of matrix type.
!> Must only be called after the matrix has been written and before 
!> file is rewinded!
!>
      subroutine mat_write_info_to_disk(iunit,info)
         implicit none
         integer, intent(in) :: iunit
         logical, intent(in) :: info
         logical(kind=8) :: info8
         !always read (and write) using 64 bit logicals to ensure uniform
         !read and write format independent of compilation
         info8 = info
         write(iunit) info8
      end subroutine mat_write_info_to_disk

!> \brief Read a type(matrix) from disk.
!> \author L. Thogersen
!> \date 2003
!> \param iunit Logical unit number of file from which matrix should be read
!> \param A Output matrix which is read from disk
      subroutine mat_read_from_disk(iunit,A,OnMaster)
         implicit none
         integer, intent(in) :: iunit
         type(Matrix), intent(inout) :: A  !output
         logical,optional :: OnMaster !obsolete
         !
         real(realk), allocatable :: afull(:,:)
         integer(kind=long)       :: nrow, ncol
         call time_mat_operations1
         if (info_memory) write(mat_lu,*) 'Before mat_read_from_disk: mem_allocated_global =', mem_allocated_global
         !if (INFO_TIME_MAT) CALL LSTIMER('START ',mat_TSTR,mat_TEN,mat_lu)
         select case(matrix_type)
         case(mtype_dense)
            call mat_dense_read_from_disk(iunit,A)
         case(mtype_csr)
            call mat_csr_read_from_disk(iunit,A)
         case(mtype_scalapack)
            !The master Read full matrix from disk 
            print *, "FALLBACK: mat_read_from_disk"
            allocate(afull(a%nrow, a%ncol))
            READ(iunit) Nrow, Ncol
            if(Nrow /= A%nrow) call lsquit( 'mat_read_from_disk: Nrow /= A%nrow',-1)
            if(Ncol /= A%ncol) call lsquit( 'mat_read_from_disk: Ncol /= A%ncol',-1)
            read(iunit) afull
            call mat_set_from_full(afull,1E0_realk,a)
            deallocate(afull)
         case(mtype_unres_dense)
            call mat_unres_dense_read_from_disk(iunit,A)
         case default
            print *, "FALLBACK: mat_read_from_disk"
            allocate(afull(a%nrow, a%ncol))
            READ(iunit) Nrow, Ncol
            if(Nrow /= A%nrow) call lsquit( 'mat_read_from_disk: Nrow /= A%nrow',-1)
            if(Ncol /= A%ncol) call lsquit( 'mat_read_from_disk: Ncol /= A%ncol',-1)
            read(iunit) afull
            call mat_set_from_full(afull,1E0_realk,a)
            deallocate(afull)
         end select

         !if (INFO_TIME_MAT) CALL LSTIMER('READ  ',mat_TSTR,mat_TEN,mat_lu)
         if (info_memory) write(mat_lu,*) 'After mat_read_from_disk: mem_allocated_global =', mem_allocated_global
         call time_mat_operations2(JOB_mat_read_from_disk)

      end subroutine mat_read_from_disk

!> \brief Read some logical auxiliary information from a type(matrix) on disk.
!> \author S. Host
!> \date June 2010
!> \param iunit Logical unit number of file containing the matrix
!> \param info Info to be read
!>
!> See description in mat_write_info_to_disk.
!>
      subroutine mat_read_info_from_disk(iunit,info)
        implicit none
        integer, intent(in)  :: iunit
        logical, intent(out) :: info
        logical(kind=8) :: info8
        !always read (and write) using 64 bit logicals to ensure uniform
        !read and write format independent of compilation
        read(iunit) info8
        info = info8
      end subroutine mat_read_info_from_disk

!> \brief Extract diagonal of A, store in dense vector vec.
!> \author B. Jansik
!> \date 2010
!> \param A The type(matrix) input
!> \param diag vector to hold the diagonal
    subroutine mat_extract_diagonal (diag,A)
      implicit none
      type(Matrix), intent(in) :: A
      real(realk), intent(inout) :: diag(A%nrow)

      select case(matrix_type)
      case(mtype_dense)
           call mat_dense_extract_diagonal(diag,A)
      case(mtype_scalapack)
           call mat_scalapack_extract_diagonal(diag,A)
#ifdef VAR_ENABLE_TENSORS
      case(mtype_pdmm)
           call mat_pdmm_extract_diagonal(diag,A)
#endif
      case default
            call lsquit("mat_extract_diagonal not implemented for this type of matrix",-1)
      end select

    end subroutine mat_extract_diagonal

END MODULE Matrix_Operations

#ifdef VAR_MPI
      !> \brief Pass matrix_type to other MPI nodes
      !> \author T. Kjaergaard
      !> \date 2011
      !> \param a the matrix type
      subroutine lsmpi_set_matrix_type_master(a)
  use infpar_module
  use lsmpi_type
        implicit none
        integer,intent(in) :: a
        call ls_mpibcast(a,infpar%master,MPI_COMM_LSDALTON)
      end subroutine lsmpi_set_matrix_type_master

      !> \brief obtains the matrix_type from master MPI nodes
      !> \author T. Kjaergaard
      !> \date 2011
      !> \param the matrix type
      subroutine lsmpi_set_matrix_type_slave()
  use matrix_operations
  use infpar_module
  use lsmpi_type
        implicit none
        integer :: a
        call ls_mpibcast(a,infpar%master,MPI_COMM_LSDALTON)

        call mat_select_type(a,6)

      end subroutine lsmpi_set_matrix_type_slave

      !> \brief obtains the matrix_type from master MPI nodes
      !> \author T. Kjaergaard
      !> \date 2011
      !> \param the matrix type
      subroutine lsmpi_set_matrix_tmp_type_slave()
        use matrix_operations
        use infpar_module
        use lsmpi_type
        implicit none
        integer :: a
        call ls_mpibcast(a,infpar%master,MPI_COMM_LSDALTON)
        call mat_select_tmp_type(a,6)

      end subroutine lsmpi_set_matrix_tmp_type_slave
#endif
