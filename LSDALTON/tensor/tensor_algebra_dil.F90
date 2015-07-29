!This module provides an infrastructure for distributed tensor algebra
!that avoids loading full tensors into RAM of a single node.
!AUTHOR: Dmitry I. Lyakh: quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2015/04/10 (started 2014/09/01).
!DISCLAIMER:
! This code was developed in support of the INCITE project CHP100
! at the National Center for Computational Sciences at
! the Oak Ridge National Laboratory (Oak Ridge TN, USA)
! managed by UT-Battelle for the US Department of Energy.
! The following code is open source, free of any liabilities.
! The author and UT-Battelle keep their rights on reusing
! any part of this code elsewhere.
!DETAILS:
! * Two types of tensors are considered here:
!   1) Local dense tensors (stored fully in local RAM);
!   2) Distributed tensors (stored in terms of dense tiles on multiple nodes);
! * In general, a tensor contraction will operate on subtensors, that is,
!   the operands do not have to be the full tensors. The corresponding slicing,
!   as well as the contraction pattern, are specified via <dil_tens_contr_t>.
! * A general tensor contraction operating on (sub)tensors (input) is
!   split into multiple smaller tensor contractions (tasks) operating
!   on parts of those (sub)tensors. Each part of each (sub)tensor
!   (a) is small enough to fit into the local RAM buffer;
!   (b) consists of an integer number of whole tiles the tensor is stored in terms of
!       (this applies to distributed tensors only).
!   Ultimately, a part of a distributed tensor can consist of a single tile, but no less.
!   No such restrictions are imposed on the parts of local dense tensors.
! * Dimensions of a distributed tensor are split into segments,
!   thus defining the tiles the tensor is stored in terms of.
!   Grouping multiple segments into supersegments defines tensor parts
!   participating in a particular elementary tensor contraction (task).
!   Thus, a (sub)tensor can be processed in terms of its parts,
!   each (sub)tensor part consisting of one or more whole tiles.
! * Six contraction cases (<alpha> is a scalar):
!   1) LLL: Local += Local * Local * alpha;
!   2) LLD: Local += Local * Distributed * alpha;
!   3) DLL: Distributed += Local * Local * alpha;
!   4) LDD: Local += Distributed * Distributed * alpha;
!   5) DLD: Distributed += Local * Distributed * alpha;
!   6) DDD: Distributed += Distributed * Distributed * alpha.
! * Only the destination tensor can be of rank 0 (scalar);
!   Both left and right tensor arguments must always be true tensors (rank>0)!
!   Thus, a simple multiplication of a tensor by a scalar is not considered here.
! * Tensor dimensions start from 1 (Fortran style). Thus, a tensor dimension
!   will cover the following range: [1:dim_extent]. Also, it is assumed that
!   seniority of tensor dimensions increases from left to right (Fortran-like).
!   Each tensor dimension is divided into same-length segments (the last
!   terminating segment may have a different length). Distinct dimensions of a tensor
!   may have different segment lengths. The tensor tiles are enumerated in Fortran style
!   (1st dimenstion, 2nd dimension, and so on), flat (global) tile numeration starts from 1.
!NOTES:
! * The code assumes Fortran-2003/2008 & MPI-3 (defined COMPILER_UNDERSTANDS_FORTRAN_2003, VAR_PTR_RESHAPE, VAR_MPI)!
! * The number of OMP threads spawned on CPU or MIC must not exceed the MAX_THREADS parameter!
! * In order to activate the debugging mode, define the DIL_DEBUG_ON macro below.
!PREPROCESSOR:
! * VAR_OMP: use OpenMP;
! * USE_OMP_MOD: use omp_lib module;
! * USE_BASIC_ALLOC: disable MPI_ALLOC_MEM() calls and stick to basic allocate()/malloc();
! * DIL_DEBUG_ON: enable debugging information;
! * USE_MIC: enable Intel MIC accelerators (not yet implemented).
       module tensor_algebra_dil
        use lspdm_tensor_operations_module
#ifdef COMPILER_UNDERSTANDS_FORTRAN_2003
#ifdef VAR_PTR_RESHAPE
#ifdef VAR_MPI
#define DIL_ACTIVE
#define DIL_DEBUG_ON
#endif
#endif
#endif

#ifdef DIL_ACTIVE
        use, intrinsic:: ISO_C_BINDING
#ifdef VAR_OMP
#ifdef USE_OMP_MOD
        use omp_lib
        implicit none
#else
        implicit none
        integer(4), external, private:: omp_get_max_threads,omp_get_num_threads,omp_get_thread_num
        real(8), external, private:: omp_get_wtime
#endif
#else
        implicit none
#endif
!PARAMETERS:
        integer(4), parameter, public:: INTD=4                          !default integer kind (size)
        integer(4), parameter, public:: INTL=8                          !long integer kind (size)
        integer(INTD), parameter, private:: BLAS_INT=INTD               !default integer size for BLAS/LAPACK
        integer(INTL), parameter, private:: MIN_BUF_MEM=128*1048576_INTL!min allowed local memory limit in bytes for buffer space
        integer(INTL), parameter, private:: ALIGNMENT=4096              !buffer alignment in bytes (must be multiple of 16)
        integer(INTD), parameter, private:: MAX_TILES_PER_PART=1024     !max allowed number of tiles per tensor part
        integer(INTD), parameter, private:: MIN_LOC_DIM_EXT=16          !minimal tiling length for local tensors
        integer(INTD), parameter, public:: MAX_TENSOR_RANK=16           !max allowed tensor rank
        integer(INTD), parameter, public:: IND_NUM_START=1              !number from which index range numeration starts
        integer(INTD), parameter, private:: MAX_THREADS=1024            !max number of CPU threads on multi-core CPU (or MIC)
        integer(INTD), parameter, private:: BUFS_PER_DEV=7              !number of buffers per device (3 prefetch, 3 compute, 1 upload)
        integer(INTD), parameter, private:: TASK_NULL=0                 !status "empty task"
        integer(INTD), parameter, private:: TASK_SET=1                  !status "existing task" (fully or partially set)
        integer(INTD), parameter, private:: TASK_SCHEDULED=2            !status "scheduled task"
        integer(INTD), parameter, private:: TASK_COMPLETED=3            !status "completed task" (successfully)
        integer(INTD), parameter, private:: TASK_ERR_LDS=-1             !error status "load start failed"
        integer(INTD), parameter, private:: TASK_ERR_LDC=-2             !error status "load complete failed"
        integer(INTD), parameter, private:: TASK_ERR_PRI=-3             !error status "prepare input failed"
        integer(INTD), parameter, private:: TASK_ERR_CMT=-4             !error status "compute failed"
        integer(INTD), parameter, private:: TASK_ERR_PRO=-5             !error status "prepare output failed"
        integer(INTD), parameter, private:: TASK_ERR_STS=-6             !error status "store start failed"
        integer(INTD), parameter, private:: TASK_ERR_STC=-7             !error status "store complete failed"
        integer(INTD), parameter, private:: DIL_FIRST_CALL=1            !first call of a subroutine
        integer(INTD), parameter, private:: DIL_DONE=-1                 !iterations are over
        integer(INTD), parameter, private:: DIL_NO_WORK=-1              !no work for an MPI processes
        integer(INTD), parameter, public:: DIL_SUCCESS=0                !success
        integer(INTD), parameter, public:: DIL_ALLOC_NOT=-1             !not allocated
        integer(INTD), parameter, public:: DIL_ALLOC_BASIC=0            !basic memory allocation
        integer(INTD), parameter, public:: DIL_ALLOC_PINNED=1           !pinned memory allocation
        integer(INTD), parameter, public:: DIL_ALLOC_MPI=2              !memory allocation by MPI_ALLOC_MEM
        integer(INTD), parameter, public:: DIL_ALLOC_EXT=3              !pointer associated with an external allocation
        integer(INTD), parameter, public:: DEV_HOST_CPU=0               !Regular CPU device (multi-core)
        integer(INTD), parameter, public:: DEV_NVIDIA_GPU=1             !NVidia GPU device
        integer(INTD), parameter, public:: DEV_INTEL_MIC=2              !Intel MIC device
        integer(INTD), parameter, public:: MAX_GPUS=4                   !max number of GPUs per node (MPI process)
        integer(INTD), parameter, public:: MAX_MICS=4                   !max number of Intel MICs per node (MPI process)
        integer(INTD), parameter, public:: MAX_DEVS=1+MAX_GPUS+MAX_MICS !max number of devices on a node (DEV#0 is always a (multicore) CPU)
        logical, parameter, public:: DIL_TC_EACH=.false.                !each MPI process will have its own tensor contraction
        logical, parameter, public:: DIL_TC_ALL=.true.                  !MPI processes work collectively on a tensor contraction
!VARIABLES:
#ifndef USE_BASIC_ALLOC
        integer(INTD), private:: DIL_ALLOC_TYPE=DIL_ALLOC_MPI           !default allocator type for communication buffers
#else
        integer(INTD), private:: DIL_ALLOC_TYPE=DIL_ALLOC_BASIC         !default allocator type for communication buffers
#endif
        integer(INTD), private:: CONS_OUT=6,CONS_OUT_SAVED=0            !console output device (defaults to screen)
        integer(INTD), public:: DIL_CONS_OUT=6                          !console output device for external use (defaults to screen)
        logical, private:: VERBOSE=.true.                               !verbosity (for errors)
#ifdef DIL_DEBUG_ON
        logical, public:: DIL_DEBUG=.true.                              !debugging
#else
        logical, public:: DIL_DEBUG=.false.                             !debugging
#endif
        integer(INTD), private:: DIL_DEBUG_FILE=666                     !debug file handle
        integer(INTD), private:: DIL_TMP_FILE1=1043                     !temporary file handle
        integer(INTD), private:: DIL_TMP_FILE2=1044                     !temporary file handle
#ifdef USE_MIC
!DIR$ ATTRIBUTES OFFLOAD:mic:: INTD,INTL,BLAS_INT,CONS_OUT,VERBOSE,DIL_DEBUG,MAX_TENSOR_RANK,MAX_THREADS,IND_NUM_START
!DIR$ ATTRIBUTES ALIGN:128:: INTD,INTL,BLAS_INT,CONS_OUT,VERBOSE,DIL_DEBUG,MAX_TENSOR_RANK,MAX_THREADS,IND_NUM_START
#endif
!TYPES:
 !Rank/Window descriptor:
        type, private:: rank_win_t
         integer(INTD), private:: rank
         integer(INTL), private:: window
        end type rank_win_t
 !Rank/window container:
        type, private:: rank_win_cont_t
         integer(INTD), private:: num_entries=0
         integer(INTD), private:: first_entry=-1
         type(rank_win_t), private:: rw_entry(MAX_TILES_PER_PART)
         integer(INTD), private:: next_win(MAX_TILES_PER_PART)
         integer(INTD), private:: next_rank(MAX_TILES_PER_PART)
        end type rank_win_cont_t
 !Tensor contraction specification:
        type, private:: contr_spec_t
         logical, private:: dest_zero !if .true., the LOCAL destination (sub)tensor will be set to zero before tensor contraction
         integer(INTD), private:: ndims_left               !number of uncontracted dimensions coming from the left tensor
         integer(INTD), private:: ndims_right              !number of uncontracted dimensions coming from the right tensor
         integer(INTD), private:: ndims_contr              !number of contracted dimensions
         integer(INTD), private:: dprmn(1:MAX_TENSOR_RANK) !O2N permutation required for the matricization of the destination (sub)tensor
         integer(INTD), private:: lprmn(1:MAX_TENSOR_RANK) !O2N permutation required for the matricization of the left (sub)tensor
         integer(INTD), private:: rprmn(1:MAX_TENSOR_RANK) !O2N permutation required for the matricization of the right (sub)tensor
         integer(INTD), private:: ddims(1:MAX_TENSOR_RANK) !dimension extents of the destination (sub)tensor (in order they are stored)
         integer(INTD), private:: ldims(1:MAX_TENSOR_RANK) !dimension extents of the left (sub)tensor (in order they are stored)
         integer(INTD), private:: rdims(1:MAX_TENSOR_RANK) !dimension extents of the right (sub)tensor (in order they are stored)
         integer(INTD), private:: dbase(1:MAX_TENSOR_RANK) !dimension base offsets for the destination (sub)tensor (signature)
         integer(INTD), private:: lbase(1:MAX_TENSOR_RANK) !dimension base offsets for the left (sub)tensor (signature)
         integer(INTD), private:: rbase(1:MAX_TENSOR_RANK) !dimension base offsets for the right (sub)tensor (signature)
        end type contr_spec_t
 !Locally stored tensor:
        type, private:: tens_loc_t
         integer(INTD), private:: rank                                !tensor rank (number of dimensions)
         integer(INTD), private:: dims(1:MAX_TENSOR_RANK)             !tensor dimension extents: dims(1:rank)
         integer(INTD), private:: base(1:MAX_TENSOR_RANK)             !offsets specifying a subtensor (similar to cspec%bases)
!         real(tensor_dp), pointer, contiguous, private:: elems(:)=>NULL() !tensor elements (1:*)
         real(tensor_dp), pointer, contiguous :: elems(:)=>NULL() !tensor elements (1:*)
        end type tens_loc_t
 !Tensor argument:
        type, private:: tens_arg_t
         type(tens_loc_t), private:: tens_loc                  !description of a local tensor
         type(tensor), pointer, private:: tens_distr_p=>NULL() !description of a distributed tensor (pointer)
         character(1), private:: store_type=' '                !storage type: {"L":local; "D":distributed}
        end type tens_arg_t
 !User-level complete tensor contraction specification:
        type, public:: dil_tens_contr_t
         integer(INTD), private:: contr_id         !unique tensor contraction id (`Reserved)
         type(contr_spec_t), private:: contr_spec  !formal tensor contraction specification
         type(tens_arg_t), private:: dest_arg      !destination tensor argument
         type(tens_arg_t), private:: left_arg      !left tensor argument
         type(tens_arg_t), private:: right_arg     !right tensor argument
         real(tensor_dp), private:: alpha              !multiplication prefactor (can be explicit)
         real(tensor_dp), private:: beta               !scaling factor for the destination tensor (always implicit)
         integer(INTD), private:: num_async        !number of asynchronous outstanding MPI uploads left after the contraction
         type(rank_win_t), private:: list_async(1:MAX_TILES_PER_PART)  !asynchronous outstanding MPI uploads left after the contraction
#ifdef VAR_PGF90
         real(tensor_dp), pointer, contiguous :: buffer(:)=>NULL() !work buffer
#else
         real(tensor_dp), pointer, contiguous, private:: buffer(:)=>NULL() !work buffer
#endif
         integer(INTD), private:: alloc_type=DIL_ALLOC_NOT             !allocation type for the work buffer
        end type dil_tens_contr_t
 !Subtensor part specification:
        type, public:: subtens_t
         integer(INTD), private:: rank                    !subtensor rank (number of dimensions)
         integer(INTD), private:: lbnd(1:MAX_TENSOR_RANK) !dimension lower bounds
         integer(INTD), private:: dims(1:MAX_TENSOR_RANK) !dimension extents
        end type subtens_t
 !Tensor contraction task:
        type, private:: contr_task_t
         type(subtens_t), private:: dest_arg          !destination subtensor specification
         type(subtens_t), private:: left_arg          !left subtensor specification
         type(subtens_t), private:: right_arg         !right subtensor specification
         real(8), private:: prefac=1d0                !prefactor
         integer(INTL), private:: flops               !number of Flops in the task
         integer(INTD), private:: dev_kind            !device kind
         integer(INTD), private:: dev_id              !device id
         integer(INTD), private:: task_stat=TASK_NULL !task status (negative means error)
         real(8), private:: time_start                !start time of task execution
         real(8), private:: time_finish               !finish time of task execution
        end type contr_task_t
 !Task list (list of tensor contraction tasks):
        type, private:: contr_task_list_t
         integer(INTD), private:: num_tasks=0                          !number of tasks in the list
         logical, private:: alloc=.false.                              !if .true., contr_tasks(:) was allocated (not just associated)
         logical, private:: reordered=.false.                          !if .true., the tasks have been reordered
         type(contr_task_t), pointer, private:: contr_tasks(:)=>NULL() !list of tensor contraction tasks (numeration starts from 1)
         integer(INTD), allocatable, private:: task_order(:)           !task order (if reordering has been done on the task list)
        end type contr_task_list_t
 !Argument buffer:
        type, private:: arg_buf_t
         integer(INTL), private:: buf_vol=0_INTL                        !buffer volume (number of elements)
#ifdef VAR_PGF90
         real(tensor_dp), pointer, contiguous :: buf_ptr(:)=>NULL() !buffer pointer
#else
         real(tensor_dp), pointer, contiguous, private:: buf_ptr(:)=>NULL() !buffer pointer
#endif
        end type arg_buf_t
 !Device buffers:
        type, private:: dev_buf_t
         integer(INTD), private:: num_bufs=BUFS_PER_DEV     !number of buffers belonging to device
         type(arg_buf_t), private:: arg_buf(1:BUFS_PER_DEV) !buffers belonging to device
        end type dev_buf_t
!DATA:

!INTERFACES:
 !CPU pointer memory allocation:
        interface cpu_ptr_alloc
         module procedure cpu_ptr_alloc_r !deafult real kind (tensor_dp)
        end interface cpu_ptr_alloc
 !CPU pointer memory deallocation:
        interface cpu_ptr_free
         module procedure cpu_ptr_free_r !deafult real kind (tensor_dp)
        end interface cpu_ptr_free
 !Integer to String conversion:
        interface int2str
         module procedure int2str_i4
         module procedure int2str_i8
        end interface int2str
 !Merge sort:
        interface merge_sort_key_int
         module procedure merge_sort_key_int4
         module procedure merge_sort_key_int8
        end interface merge_sort_key_int
 !Dividing n-dimensional space into blocks:
        interface dil_divide_space_int
         module procedure dil_divide_space_int4
         module procedure dil_divide_space_int8
        end interface dil_divide_space_int
!VISIBILITY:
        private dil_tens_contr_spec_setup
        private dil_tens_contr_spec_check
        private dil_tens_arg_clean
        private dil_tens_arg_set
        public dil_clean_tens_contr
        public dil_set_tens_contr_args
        private dil_get_arg_tile_vol
        public dil_get_min_buf_size
        public dil_prepare_buffer
        public dil_set_tens_contr_spec
        public dil_subtensor_set
        private dil_subtensor_vol
        private dil_subtensor_copy
        private dil_subtensor_print
        private dil_contr_task_set_arg
        private dil_contr_task_set_flops
        private dil_contr_task_copy
        private dil_contr_task_print
        private dil_contr_task_list_destroy
        private dil_contr_task_list_create
        private dil_contr_task_list_shuffle
        private dil_contr_task_list_print
        private dil_dev_num
        private dil_dev_kind
        private dil_dev_buf_destroy
        private dil_arg_buf_clean
        public dil_set_alloc_type
        private cpu_ptr_alloc
        private cpu_ptr_alloc_r
        private cpu_ptr_free
        private cpu_ptr_free_r
        private my_mpi_size
        private my_mpi_rank
        private divide_segment_i8
        public merge_sort_key_int
        private merge_sort_key_int4
        private merge_sort_key_int8
        public merge_sort_real8
        private str_parse
        public str2int
        public int2str
        private int2str_i4
        private int2str_i8
        private printf
        private dil_rank_window_clean
        private dil_rank_window_new
        private dil_get_contr_pattern
        public thread_wtime
        public process_wtime
        public permutation_invert
        public permutation_trivial
        public dil_tensor_slice
        public dil_tensor_insert
        public dil_tensor_transpose
        private dil_get_next_tile_signa
        private dil_tensor_prefetch_start
        private dil_tensor_prefetch_complete
        private dil_tensor_upload_start
        private dil_tensor_upload_complete
        private dil_tens_unpack_from_tiles
        private dil_tens_pack_into_tiles
        public dil_tens_fetch_start
        public dil_tens_fetch_finish_prep
!        public dil_tens_fetch
        public dil_tens_prep_upload_start
        public dil_tens_upload_finish
!        public dil_tens_upload
        private dil_divide_space_int
        private dil_divide_space_int4
        private dil_divide_space_int8
        private dil_tens_contr_distribute
        private dil_tens_contr_partition
        private dil_tensor_contract_pipe
        public dil_tensor_contract
        public dil_tensor_contract_finalize
        public dil_tensor_norm1
        public dil_array_norm1
        public dil_tensor_init
        public dil_array_init
        public dil_debug_to_file_start
        public dil_debug_to_file_finish
        public dil_array_print
        public dil_mm_pipe_efficient
        public dil_will_malloc_succeed
        public dil_test

       contains
!-------------------------------------------------------------------------------------------------------------------------------
        subroutine dil_tens_contr_spec_setup(cspec,drank,lrank,rrank,dprm,lprm,rprm,ddim,ldim,rdim,ierr,dbas,lbas,rbas,dst_zero) !SERIAL
!This subroutine sets up a tensor contraction specification <cspec>.
!Optional bases specify the offsets defining a subtensor:
!Each tensor dimension starts at 1; subtensor dimensions start at 1+base.
!Subtensor dimensions: [1+base:1+base+dim-1].
        implicit none
        type(contr_spec_t), intent(out):: cspec             !out: contraction specification
        integer(INTD), intent(in):: drank                   !in: destination tensor rank
        integer(INTD), intent(in):: lrank                   !in: left tensor rank
        integer(INTD), intent(in):: rrank                   !in: right tensor rank
        integer(INTD), intent(in):: dprm(1:drank)           !in: O2N that matricize the destination tensor
        integer(INTD), intent(in):: lprm(1:lrank)           !in: O2N that matricize the left tensor
        integer(INTD), intent(in):: rprm(1:rrank)           !in: O2N that matricize the right tensor
        integer(INTD), intent(in):: ddim(1:drank)           !in: destination subtensor dimension extents (as it is stored)
        integer(INTD), intent(in):: ldim(1:lrank)           !in: left subtensor dimension extents (as it is stored)
        integer(INTD), intent(in):: rdim(1:rrank)           !in: right subtensor dimension extents (as it is stored)
        integer(INTD), intent(inout):: ierr                 !out: error code (0:success)
        integer(INTD), intent(in), optional:: dbas(1:drank) !in: destination subtensor index bases (default is 0)
        integer(INTD), intent(in), optional:: lbas(1:lrank) !in: left subtensor index bases (default is 0)
        integer(INTD), intent(in), optional:: rbas(1:rrank) !in: right subtensor index bases (default is 0)
        logical, intent(in), optional:: dst_zero !in: if .true., the initial state of the destination tensor is assumed zero
        logical, parameter:: NO_CHECK=.false.    !argument check
        integer(INTD):: i

        ierr=0
        if(drank.ge.0.and.lrank.ge.0.and.rrank.ge.0.and.lrank+rrank.ge.drank.and.mod(lrank+rrank-drank,2).eq.0.and. &
           drank.le.MAX_TENSOR_RANK.and.lrank.le.MAX_TENSOR_RANK.and.rrank.le.MAX_TENSOR_RANK) then
         cspec%ndims_contr=(lrank+rrank-drank)/2
         cspec%ndims_left=lrank-cspec%ndims_contr
         cspec%ndims_right=rrank-cspec%ndims_contr
         cspec%dprmn(1:drank)=dprm(1:drank)
         cspec%lprmn(1:lrank)=lprm(1:lrank)
         cspec%rprmn(1:rrank)=rprm(1:rrank)
         cspec%ddims(1:drank)=ddim(1:drank)
         cspec%ldims(1:lrank)=ldim(1:lrank)
         cspec%rdims(1:rrank)=rdim(1:rrank)
         if(present(dbas)) then
          cspec%dbase(1:drank)=dbas(1:drank)
         else
          cspec%dbase(1:drank)=0
         endif
         if(present(lbas)) then
          cspec%lbase(1:lrank)=lbas(1:lrank)
         else
          cspec%lbase(1:lrank)=0
         endif
         if(present(rbas)) then
          cspec%rbase(1:rrank)=rbas(1:rrank)
         else
          cspec%rbase(1:rrank)=0
         endif
         if(present(dst_zero)) then
          cspec%dest_zero=dst_zero
         else
          cspec%dest_zero=.false. !by default, a tensor contraction assumes accumulation
         endif
         if(.not.NO_CHECK) then
          i=dil_tens_contr_spec_check(cspec)
          if(i.ne.0) then
           if(VERBOSE) write(CONS_OUT,'("#ERROR(dil_tens_contr_spec_setup): invalid tensor contraction specification: ",i9)') i
           ierr=1
          endif
         endif
        else
         ierr=2
        endif
        return
        end subroutine dil_tens_contr_spec_setup
!--------------------------------------------------------------
        integer(INTD) function dil_tens_contr_spec_check(cspec) !SERIAL
!This function checks a tensor contraction specification <cspec>.
        implicit none
        type(contr_spec_t), intent(in):: cspec !in: tensor contraction specification
        integer(INTD):: i,j,k,nd,nl,nr,prm0(1:MAX_TENSOR_RANK),prm1(1:MAX_TENSOR_RANK),prm2(1:MAX_TENSOR_RANK)

        dil_tens_contr_spec_check=0
        if(cspec%ndims_left.lt.0) then; dil_tens_contr_spec_check=1; return; endif
        if(cspec%ndims_right.lt.0) then; dil_tens_contr_spec_check=2; return; endif
        if(cspec%ndims_contr.lt.0) then; dil_tens_contr_spec_check=3; return; endif
        nd=cspec%ndims_left+cspec%ndims_right  !destination tensor rank
        nl=cspec%ndims_contr+cspec%ndims_left  !left tensor rank
        nr=cspec%ndims_contr+cspec%ndims_right !right tensor rank
        if(nd.gt.MAX_TENSOR_RANK) then; dil_tens_contr_spec_check=4; return; endif
        if(nl.gt.MAX_TENSOR_RANK) then; dil_tens_contr_spec_check=5; return; endif
        if(nr.gt.MAX_TENSOR_RANK) then; dil_tens_contr_spec_check=6; return; endif
        do i=1,nd; if(cspec%ddims(i).le.0) then; dil_tens_contr_spec_check=7; return; endif; enddo
        do i=1,nl; if(cspec%ldims(i).le.0) then; dil_tens_contr_spec_check=8; return; endif; enddo
        do i=1,nr; if(cspec%rdims(i).le.0) then; dil_tens_contr_spec_check=9; return; endif; enddo
        prm0(1:nd)=0
        do i=1,nd
         j=cspec%dprmn(i)
         if(j.le.0.or.j.gt.nd) then; dil_tens_contr_spec_check=10; return; endif
         if(prm0(j).gt.0) then; dil_tens_contr_spec_check=11; return; endif; prm0(j)=prm0(j)+1
        enddo
        prm0(1:nl)=0
        do i=1,nl
         j=cspec%lprmn(i)
         if(j.le.0.or.j.gt.nl) then; dil_tens_contr_spec_check=12; return; endif
         if(prm0(j).gt.0) then; dil_tens_contr_spec_check=13; return; endif; prm0(j)=prm0(j)+1
        enddo
        prm0(1:nr)=0
        do i=1,nr
         j=cspec%rprmn(i)
         if(j.le.0.or.j.gt.nr) then; dil_tens_contr_spec_check=14; return; endif
         if(prm0(j).gt.0) then; dil_tens_contr_spec_check=15; return; endif; prm0(j)=prm0(j)+1
        enddo
        if(nd.gt.0) call permutation_invert(nd,cspec%dprmn,prm0,i)
        if(nl.gt.0) call permutation_invert(nl,cspec%lprmn,prm1,i)
        if(nr.gt.0) call permutation_invert(nr,cspec%rprmn,prm2,i)
        if(DIL_DEBUG) then
         write(CONS_OUT,'("#DEBUG(dil_tens_contr_spec_check): ddims:",64(1x,i6))') cspec%ddims(1:nd)
         write(CONS_OUT,'("#DEBUG(dil_tens_contr_spec_check): ldims:",64(1x,i6))') cspec%ldims(1:nl)
         write(CONS_OUT,'("#DEBUG(dil_tens_contr_spec_check): rdims:",64(1x,i6))') cspec%rdims(1:nr)
         write(CONS_OUT,'("#DEBUG(dil_tens_contr_spec_check): dbase:",64(1x,i6))') cspec%dbase(1:nd)
         write(CONS_OUT,'("#DEBUG(dil_tens_contr_spec_check): lbase:",64(1x,i6))') cspec%lbase(1:nl)
         write(CONS_OUT,'("#DEBUG(dil_tens_contr_spec_check): rbase:",64(1x,i6))') cspec%rbase(1:nr)
        endif
        do i=1,cspec%ndims_left
         j=prm0(i); k=prm1(cspec%ndims_contr+i)
         if(cspec%dbase(j).ne.cspec%lbase(k).or.cspec%ddims(j).ne.cspec%ldims(k)) then
          dil_tens_contr_spec_check=16; return
         endif
        enddo
        do i=1,cspec%ndims_right
         j=prm0(cspec%ndims_left+i); k=prm2(cspec%ndims_contr+i)
         if(cspec%dbase(j).ne.cspec%rbase(k).or.cspec%ddims(j).ne.cspec%rdims(k)) then
          dil_tens_contr_spec_check=17; return
         endif
        enddo
        do i=1,cspec%ndims_contr
         j=prm1(i); k=prm2(i)
         if(cspec%lbase(j).ne.cspec%rbase(k).or.cspec%ldims(j).ne.cspec%rdims(k)) then
          dil_tens_contr_spec_check=18; return
         endif
        enddo
        return
        end function dil_tens_contr_spec_check
!---------------------------------------------------
        subroutine dil_tens_arg_clean(tens_arg,ierr) !SERIAL
!This subroutine cleans a tensor argument.
        implicit none
        type(tens_arg_t), intent(inout):: tens_arg    !inout: tensor argument
        integer(INTD), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INTD):: errc

        errc=0
        nullify(tens_arg%tens_loc%elems); tens_arg%tens_loc%rank=-1
        nullify(tens_arg%tens_distr_p)
        tens_arg%store_type=' '
        if(present(ierr)) ierr=errc
        return
        end subroutine dil_tens_arg_clean
!------------------------------------------------------------------------------------------------------------
        subroutine dil_tens_arg_set(tens_arg,ttype,ierr,tens_rank,tens_dims,tens_elems,tens_bases,tens_distr) !SERIAL
!This subroutine sets up a tensor argument: either "LOCAL" or "DISTRIBUTED".
        implicit none
        type(tens_arg_t), intent(inout):: tens_arg                  !inout: tensor argument
        character(1), intent(in):: ttype                            !in: tensor storage type: {"L","D"}: local or distributed
        integer(INTD), intent(inout):: ierr                         !out: error code (0:success)
        integer(INTD), intent(in), optional:: tens_rank             !in: tensor rank (for local tensors only)
        integer(INTD), intent(in), optional:: tens_dims(1:*)        !in: tensor dimension extents (for local tensors only)
        real(tensor_dp), intent(in), target, optional:: tens_elems(1:*) !in: tensor elements (for local tensors only)
        integer(INTD), intent(in), optional:: tens_bases(1:*)       !in: base offsets (if a subtensor)
        type(tensor), intent(in), target, optional:: tens_distr     !in: distributed tensor instance (for distributed tensors only)
        integer(INTD):: i
        integer(INTL):: vol

        ierr=0; tens_arg%store_type=' '
        select case(ttype)
        case('l','L') !local tensor
         if(present(tens_rank).and.present(tens_dims).and.present(tens_elems).and.(.not.present(tens_distr))) then
          if(tens_rank.ge.0.and.tens_rank.le.MAX_TENSOR_RANK) then
           vol=1_INTL; do i=1,tens_rank; if(tens_dims(i).le.0) then; ierr=1; exit; endif; vol=vol*tens_dims(i); enddo
           if(ierr.eq.0) then
            tens_arg%tens_loc%rank=tens_rank
            tens_arg%tens_loc%dims(1:tens_rank)=tens_dims(1:tens_rank)
            tens_arg%tens_loc%elems=>tens_elems(1:vol)
            if(present(tens_bases)) then
             tens_arg%tens_loc%base(1:tens_rank)=tens_bases(1:tens_rank) !offsets from IND_NUM_START
            else
             tens_arg%tens_loc%base(1:tens_rank)=0_INTD
            endif
            tens_arg%store_type=ttype
           endif
          else
           ierr=2
          endif
         else
          ierr=3
         endif
        case('d','D') !distributed tensor
         if(present(tens_distr).and.&
         &(.not.(present(tens_rank).or.present(tens_dims).or.present(tens_elems).or.present(tens_bases)))) then
          tens_arg%tens_distr_p=>tens_distr
          tens_arg%store_type=ttype
         else
          ierr=4
         endif
        case default
         ierr=5
        end select
        return
        end subroutine dil_tens_arg_set
!---------------------------------------------------
        subroutine dil_clean_tens_contr(tcontr,ierr) !SERIAL
!This user-level subroutine cleans a full tensor contraction specification.
        implicit none
        type(dil_tens_contr_t), intent(inout), target:: tcontr !inout: full tensor contraction specification
        integer(INTD), intent(inout), optional:: ierr          !out: error code (0:success)
        integer(INTD):: errc

        errc=0
        tcontr%contr_spec%ndims_left=-1; tcontr%contr_spec%ndims_right=-1; tcontr%contr_spec%ndims_contr=-1
        call dil_tens_arg_clean(tcontr%dest_arg)
        call dil_tens_arg_clean(tcontr%left_arg)
        call dil_tens_arg_clean(tcontr%right_arg)
        tcontr%alpha=1E0_tensor_dp; tcontr%beta=0E0_tensor_dp
        tcontr%num_async=-1 !-1 means undefined because 0 will have a special meaning later!
        if(associated(tcontr%buffer)) call cpu_ptr_free(tcontr%buffer,errc,attr=tcontr%alloc_type)
        if(errc.eq.0) tcontr%alloc_type=DIL_ALLOC_NOT
        if(present(ierr)) ierr=errc
        return
        end subroutine dil_clean_tens_contr
!---------------------------------------------------------------------------------------------------------------
        subroutine dil_set_tens_contr_args(tcontr,arg,ierr,tens_rank,tens_dims,tens_elems,tens_bases,tens_distr) !SERIAL
!This user-level subroutine sets up arguments of a tensor contraction.
!More precisely, it specifies the storage layout for each argument of a tensor contraction.
!For distributed tiled tensors, everything is specified in an instance of type(tensor) <tens_distr>.
!For local dense tensors, one needs:
! (a) tensor rank <tens_rank>;
! (b) tensor dimension extents <tens_dims>;
! (c) reference to a 1d array containing tensor elements <tens_elems>;
! (d) optionally, <tens_bases> may contain positive values specifying
!     offsets in numeration for each tensor dimension, thus allowing
!     local storage of tensor slices (incomplete tensors). These offsets
!     are added to IND_NUM_START.
        implicit none
        type(dil_tens_contr_t), intent(inout), target:: tcontr      !inout: full tensor contraction specification
        character(1), intent(in):: arg                              !in: argument selector: {'d','l','r'}
        integer(INTD), intent(inout):: ierr                         !out: error code (0:success)
        integer(INTD), intent(in), optional:: tens_rank             !in: tensor rank (for local tensors only)
        integer(INTD), intent(in), optional:: tens_dims(1:*)        !in: tensor dimension extents (for local tensors only)
        real(tensor_dp), intent(in), target, optional:: tens_elems(1:*) !in: tensor elements (for local tensors only)
        integer(INTD), intent(in), optional:: tens_bases(1:*)       !in: base offsets (if a subtensor)
        type(tensor), intent(in), target, optional:: tens_distr     !in: distributed tensor instance (for distributed tensors only)
        type(tens_arg_t), pointer:: targ

        ierr=0
        select case(arg)
        case('d','D')
         targ=>tcontr%dest_arg
        case('l','L')
         targ=>tcontr%left_arg
        case('r','R')
         targ=>tcontr%right_arg
        case default
         ierr=1; return
        end select
        if(present(tens_distr).and.&
        &(.not.(present(tens_rank).or.present(tens_dims).or.present(tens_elems).or.present(tens_bases)))) then !distributed tensor
         call dil_tens_arg_set(targ,'d',ierr,tens_distr=tens_distr)
         if(ierr.ne.0) then; ierr=2; return; endif
        elseif((.not.present(tens_distr)).and.present(tens_rank).and.present(tens_dims).and.present(tens_elems)) then
         if(present(tens_bases)) then
          call dil_tens_arg_set(targ,'l',ierr,tens_rank,tens_dims,tens_elems,tens_bases)
          if(ierr.ne.0) then; ierr=3; return; endif
         else
          call dil_tens_arg_set(targ,'l',ierr,tens_rank,tens_dims,tens_elems)
          if(ierr.ne.0) then; ierr=4; return; endif
         endif
        else
         ierr=5
        endif
        return
        end subroutine dil_set_tens_contr_args
!-------------------------------------------------------------
        integer(INTL) function dil_get_arg_tile_vol(targ,ierr) !SERIAL
!This function returns the volume of a tile for a give tensor argument.
        implicit none
        type(tens_arg_t), intent(in):: targ           !in: tensor argument
        integer(INTD), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INTD):: i,errc

        errc=0; dil_get_arg_tile_vol=1_INTL
        select case(targ%store_type)
        case('l','L')
         do i=1,targ%tens_loc%rank
          dil_get_arg_tile_vol=dil_get_arg_tile_vol*MIN_LOC_DIM_EXT
         enddo
        case('d','D')
         if(associated(targ%tens_distr_p)) then
          do i=1,targ%tens_distr_p%mode
           dil_get_arg_tile_vol=dil_get_arg_tile_vol*targ%tens_distr_p%tdim(i)
          enddo
         else
          errc=1
         endif
        case default
         errc=2
        end select
        if(present(ierr)) ierr=errc
        return
        end function dil_get_arg_tile_vol
!-----------------------------------------------------------------------------
        integer(INTL) function dil_get_min_buf_size(tcontr,ierr,scale,num_dev) !SERIAL
!This function returns the minimal size of the work buffer space in bytes
!required for executing a given tensor contraction. If more than one computing device
!will be used for performing the tensor contraction, <num_dev> must be specified!
!A computing device is a set of one or more computing cores sharing the same physical memory.
        implicit none
        type(dil_tens_contr_t), intent(in):: tcontr   !in: full tensor contraction specification
        integer(INTD), intent(inout):: ierr           !out: error code (0:success)
        real(8), intent(in), optional:: scale         !in: scaling factor
        integer(INTD), intent(in), optional:: num_dev !in: number of distinct computing devices (multi-core CPU is a single device)
!----------------------------------------------
        real(8), parameter:: SCALE_FOR_SURE=2d0 !overestimate the min buf size to be sure
!----------------------------------------------
        integer(INTD):: i
        integer(INTL):: sz,ndev
        real(8):: scl

        ierr=0; dil_get_min_buf_size=0_INTL
        if(present(scale)) then; scl=max(scale,SCALE_FOR_SURE); else; scl=SCALE_FOR_SURE; endif
        if(present(num_dev)) then; ndev=int(num_dev,INTL); else; ndev=1_INTL; endif
        if(ndev.ge.1.and.ndev.le.MAX_DEVS) then
         sz=dil_get_arg_tile_vol(tcontr%dest_arg,i); if(i.ne.0) then; ierr=1; return; endif
         dil_get_min_buf_size=max(dil_get_min_buf_size,sz)
         sz=dil_get_arg_tile_vol(tcontr%left_arg,i); if(i.ne.0) then; ierr=2; return; endif
         dil_get_min_buf_size=max(dil_get_min_buf_size,sz)
         sz=dil_get_arg_tile_vol(tcontr%right_arg,i); if(i.ne.0) then; ierr=3; return; endif
         dil_get_min_buf_size=max(dil_get_min_buf_size,sz)
         dil_get_min_buf_size=(int(real(dil_get_min_buf_size,8)*scl,INTL)/ALIGNMENT+1_INTL)*ALIGNMENT*BUFS_PER_DEV*ndev*tensor_dp
         dil_get_min_buf_size=max(dil_get_min_buf_size,MIN_BUF_MEM)
        else
         ierr=4
        endif
        return
        end function dil_get_min_buf_size
!----------------------------------------------------------------------------
        subroutine dil_prepare_buffer(tcontr,mem_lim,ierr,ext_buf,alloc_type) !SERIAL
!This subroutine prepares an internal work buffer for a giving tensor contraction.
        implicit none
        type(dil_tens_contr_t), intent(inout):: tcontr           !inout: full tensor contraction specification
        integer(INTL), intent(in):: mem_lim                      !in: size of the buffer in bytes
        integer(INTD), intent(inout):: ierr                      !out: error code (0:success)
        real(tensor_dp), intent(in), target, optional:: ext_buf(1:)  !in: external buffer space
        integer(INTD), intent(in), optional:: alloc_type         !in: allocation type
        integer(INTL):: nelems,sreal
        real(tensor_dp):: val

        ierr=0
        if(associated(tcontr%buffer)) call cpu_ptr_free(tcontr%buffer,ierr,attr=tcontr%alloc_type)
        if(ierr.eq.0) then
         tcontr%alloc_type=DIL_ALLOC_NOT
         if(mem_lim.ge.MIN_BUF_MEM) then
          val=0E0_tensor_dp; sreal=int(sizeof(val),INTL)
          if(present(ext_buf)) then !associate to a preallocated external buffer
           nelems=mem_lim/sreal !associate to at most <mem_lim> bytes
           tcontr%buffer(1:nelems)=>ext_buf; tcontr%alloc_type=DIL_ALLOC_EXT
          else !allocate a buffer
           nelems=(mem_lim-1_INTL)/sreal+1_INTL !allocate at least <mem_lim> bytes
           if(present(alloc_type)) then
            call cpu_ptr_alloc(tcontr%buffer,nelems,ierr,attr=alloc_type)
            if(ierr.eq.0) tcontr%alloc_type=alloc_type
           else
            call cpu_ptr_alloc(tcontr%buffer,nelems,ierr)
            if(ierr.eq.0) tcontr%alloc_type=DIL_ALLOC_TYPE
           endif
          endif
         else
          ierr=1
         endif
        else
         if(VERBOSE) write(CONS_OUT,'("#ERROR(dil_prepare_buffer): Work buffer deallocation failed: ",i11)') ierr
         ierr=2
        endif
        return
        end subroutine dil_prepare_buffer
!--------------------------------------------------------------------------------------------------------------
        subroutine dil_set_tens_contr_spec(tcontr,tcs,ierr,ddims,ldims,rdims,dbase,lbase,rbase,alpha,dest_zero) !SERIAL
!This user-level subroutine sets up a formal tensor contraction specification.
!Note that the tensor arguments must have been set already (prior to this call)!
        implicit none
        type(dil_tens_contr_t), intent(inout):: tcontr   !inout: full tensor contraction specification
        character(*), intent(inout):: tcs                !inout: tensor contraction specification string (spaces will be removed)
        integer(INTD), intent(inout):: ierr              !out: error code (0:success)
        integer(INTD), intent(in), optional:: ddims(1:*) !in: subtensor dimension extents for the destination tensor argument
        integer(INTD), intent(in), optional:: ldims(1:*) !in: subtensor dimension extents for the left tensor argument
        integer(INTD), intent(in), optional:: rdims(1:*) !in: subtensor dimension extents for the right tensor argument
        integer(INTD), intent(in), optional:: dbase(1:*) !in: subtensor base offsets for the destination tensor argument
        integer(INTD), intent(in), optional:: lbase(1:*) !in: subtensor base offsets for the left tensor argument
        integer(INTD), intent(in), optional:: rbase(1:*) !in: subtensor base offsets for the right tensor argument
        real(tensor_dp), intent(in), optional:: alpha        !in: alpha prefactor
        logical, intent(in), optional:: dest_zero        !in: destination zero flag
        integer(INTD):: i,j,k,l,m,n,nd,nl,nr
        integer(INTD):: dprm(1:MAX_TENSOR_RANK),lprm(1:MAX_TENSOR_RANK),rprm(1:MAX_TENSOR_RANK)
        integer(INTD):: ddim(1:MAX_TENSOR_RANK),ldim(1:MAX_TENSOR_RANK),rdim(1:MAX_TENSOR_RANK)
        integer(INTD):: dbas(1:MAX_TENSOR_RANK),lbas(1:MAX_TENSOR_RANK),rbas(1:MAX_TENSOR_RANK)
        logical:: dz

        ierr=0
        k=len_trim(tcs); if(k.le.0) then; ierr=1; return; endif
        l=0; do j=1,k; if(tcs(j:j).ne.' ') then; l=l+1; tcs(l:l)=tcs(j:j); endif; enddo
        call dil_get_contr_pattern(tcs(1:l),l,nd,nl,nr,dprm,lprm,rprm,ierr)
        if(ierr.eq.0) then
         if(nd.le.MAX_TENSOR_RANK.and.nl.le.MAX_TENSOR_RANK.and.nr.le.MAX_TENSOR_RANK) then
 !Destination tensor argument:
          if(present(ddims)) then
           ddim(1:nd)=ddims(1:nd)
           if(present(dbase)) then
            dbas(1:nd)=dbase(1:nd)
           else
            select case(tcontr%dest_arg%store_type)
            case('l','L')
             dbas(1:nd)=tcontr%dest_arg%tens_loc%base(1:nd)
            case('d','D')
             dbas(1:nd)=0_INTD
            case default
             ierr=2; return
            end select
           endif
          else
           if(.not.present(dbase)) then
            select case(tcontr%dest_arg%store_type)
            case('l','L')
             ddim(1:nd)=tcontr%dest_arg%tens_loc%dims(1:nd)
             dbas(1:nd)=tcontr%dest_arg%tens_loc%base(1:nd)
            case('d','D')
             ddim(1:nd)=tcontr%dest_arg%tens_distr_p%dims(1:nd)
             dbas(1:nd)=0_INTD
            case default
             ierr=3; return
            end select
           else
            ierr=4; return
           endif
          endif
 !Left tensor argument:
          if(present(ldims)) then
           ldim(1:nl)=ldims(1:nl)
           if(present(lbase)) then
            lbas(1:nl)=lbase(1:nl)
           else
            select case(tcontr%left_arg%store_type)
            case('l','L')
             lbas(1:nl)=tcontr%left_arg%tens_loc%base(1:nl)
            case('d','D')
             lbas(1:nl)=0_INTD
            case default
             ierr=5; return
            end select
           endif
          else
           if(.not.present(lbase)) then
            select case(tcontr%left_arg%store_type)
            case('l','L')
             ldim(1:nl)=tcontr%left_arg%tens_loc%dims(1:nl)
             lbas(1:nl)=tcontr%left_arg%tens_loc%base(1:nl)
            case('d','D')
             ldim(1:nl)=tcontr%left_arg%tens_distr_p%dims(1:nl)
             lbas(1:nl)=0_INTD
            case default
             ierr=6; return
            end select
           else
            ierr=7; return
           endif
          endif
 !Right tensor argument:
          if(present(rdims)) then
           rdim(1:nr)=rdims(1:nr)
           if(present(rbase)) then
            rbas(1:nr)=rbase(1:nr)
           else
            select case(tcontr%right_arg%store_type)
            case('l','L')
             rbas(1:nr)=tcontr%right_arg%tens_loc%base(1:nr)
            case('d','D')
             rbas(1:nr)=0_INTD
            case default
             ierr=8; return
            end select
           endif
          else
           if(.not.present(rbase)) then
            select case(tcontr%right_arg%store_type)
            case('l','L')
             rdim(1:nr)=tcontr%right_arg%tens_loc%dims(1:nr)
             rbas(1:nr)=tcontr%right_arg%tens_loc%base(1:nr)
            case('d','D')
             rdim(1:nr)=tcontr%right_arg%tens_distr_p%dims(1:nr)
             rbas(1:nr)=0_INTD
            case default
             ierr=9; return
            end select
           else
            ierr=10; return
           endif
          endif
 !Other stuff:
          if(present(alpha)) then; tcontr%alpha=alpha; else; tcontr%alpha=1E0_tensor_dp; endif
          if(tcontr%dest_arg%store_type.eq.'l'.or.tcontr%dest_arg%store_type.eq.'L') then
           tcontr%beta=1E0_tensor_dp !local destination arguments use accumulation in GEMM
          else
           tcontr%beta=0E0_tensor_dp !distributed destination arguments use assignment in GEMM (accumulated via MPI)
          endif
          if(present(dest_zero)) then; dz=dest_zero; else; dz=.false.; endif
 !Set up a formal tensor contraction specification:
          call dil_tens_contr_spec_setup(tcontr%contr_spec,nd,nl,nr,dprm,lprm,rprm,ddim,ldim,rdim,ierr,dbas,lbas,rbas,dz)
          if(ierr.ne.0) then
           if(VERBOSE) write(CONS_OUT,'("#ERROR(tensor_algebra_dil::dil_set_tens_contr_spec): Setup failed: ",i9)') ierr
           ierr=11; return
          endif
         else
          ierr=12
         endif
        else
         ierr=13
        endif
        return
        end subroutine dil_set_tens_contr_spec
!-------------------------------------------------------------
        subroutine dil_subtensor_set(subtens,srank,lb,ub,ierr) !SERIAL
!This subroutine sets dimensions of a contraction task argument (subtensor).
        implicit none
        type(subtens_t), intent(out):: subtens  !out: subtensor to set
        integer(INTD), intent(in):: srank       !in: subtensor rank
        integer(INTD), intent(in):: lb(1:srank) !in: dimension lower bounds
        integer(INTD), intent(in):: ub(1:srank) !in: dimension upper bounds
        integer(INTD), intent(inout):: ierr     !out: error code (0:success)
        integer(INTD):: i

        ierr=0
        if(srank.ge.0) then
         subtens%rank=srank
         do i=1,srank
          if(lb(i).le.ub(i)) then
           subtens%lbnd(i)=lb(i); subtens%dims(i)=ub(i)-lb(i)+1
          else
           ierr=2; exit
          endif
         enddo
        else
         ierr=1
        endif
        return
        end subroutine dil_subtensor_set
!--------------------------------------------------------
        function dil_subtensor_vol(subtens) result(stvol) !SERIAL
!This function computes the volume of a subtensor.
!No error checks.
        implicit none
        type(subtens_t), intent(in):: subtens !in: subtensor specification
        integer(INTL):: stvol                 !out: subtensor volume
        integer(INTD):: i

        stvol=1_INTL
        do i=1,subtens%rank
         stvol=stvol*subtens%dims(i)
        enddo
        return
        end function dil_subtensor_vol
!----------------------------------------------------
        subroutine dil_subtensor_copy(obji,objo,ierr) !SERIAL
        implicit none
        type(subtens_t), intent(in):: obji            !in: original object
        type(subtens_t), intent(out):: objo           !out: clone
        integer(INTD), intent(inout), optional:: ierr !error code (0:success)
        integer(INTD):: errc

        errc=0
        if(obji%rank.le.MAX_TENSOR_RANK) then
         objo%rank=obji%rank
         objo%lbnd(1:obji%rank)=obji%lbnd(1:obji%rank)
         objo%dims(1:obji%rank)=obji%dims(1:obji%rank)
        else
         errc=1
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine dil_subtensor_copy
!-----------------------------------------------------
        subroutine dil_subtensor_print(subtens,str,sl) !SERIAL
!This subroutine prints a subtensor into a string.
        implicit none
        type(subtens_t), intent(in):: subtens !in: subtensor specification
        character(*), intent(inout):: str     !out: output string
        integer(INTD), intent(out):: sl       !out: length of the output string
        integer(INTD):: i,j

        sl=0
        if(subtens%rank.ge.0) then
         str(sl+1:sl+1)='('; sl=sl+1
         do i=1,subtens%rank
          call int2str(subtens%lbnd(i),str(sl+1:),j); sl=sl+j
          sl=sl+1; str(sl:sl)=':'
          call int2str(subtens%lbnd(i)+subtens%dims(i)-1_INTD,str(sl+1:),j); sl=sl+j
          sl=sl+1; str(sl:sl)=','
         enddo
         if(str(sl:sl).eq.',') sl=sl-1
         sl=sl+1; str(sl:sl)=')'
        else !negative rank (error)
         str(sl+1:sl+3)='???'; sl=sl+3
        endif
        return
        end subroutine dil_subtensor_print
!----------------------------------------------------------------------------
        subroutine dil_contr_task_set_arg(contr_task,arg_ch,srank,lb,ub,ierr) !SERIAL
!This subroutine sets up an argument of a tensor contraction task.
        implicit none
        type(contr_task_t), intent(inout):: contr_task !inout: contraction task
        character(1), intent(in):: arg_ch              !in: argument selector: {'d','l','r'}
        integer(INTD), intent(in):: srank              !in: tensor rank
        integer(INTD), intent(in):: lb(1:srank)        !in: dimension lower bounds
        integer(INTD), intent(in):: ub(1:srank)        !in: dimension upper bounds
        integer(INTD), intent(inout):: ierr            !out: error code (0:success)

        ierr=0
        if(srank.ge.0) then
         select case(arg_ch)
         case('d','D')
          call dil_subtensor_set(contr_task%dest_arg,srank,lb,ub,ierr); if(ierr.ne.0) then; ierr=1; return; endif
          contr_task%task_stat=TASK_SET
         case('l','L')
          call dil_subtensor_set(contr_task%left_arg,srank,lb,ub,ierr); if(ierr.ne.0) then; ierr=2; return; endif
          contr_task%task_stat=TASK_SET
         case('r','R')
          call dil_subtensor_set(contr_task%right_arg,srank,lb,ub,ierr); if(ierr.ne.0) then; ierr=3; return; endif
          contr_task%task_stat=TASK_SET
         case default
          if(VERBOSE) write(CONS_OUT,'("#ERROR(tensor_algebra_dil::dil_contr_task_set_arg): invalid argument selector: ",'&
          &//'A1)') arg_ch
          ierr=4
         end select
        else
         if(VERBOSE) write(CONS_OUT,'("#ERROR(tensor_algebra_dil::dil_contr_task_set_arg): invalid tensor rank: ",i11)') srank
         ierr=5
        endif
        return
        end subroutine dil_contr_task_set_arg
!-----------------------------------------------------------------
        subroutine dil_contr_task_set_flops(cspec,contr_task,ierr) !SERIAL
!This subroutine counts and saves the number of Flops required
!for a given tensor contraction task <contr_task>.
        implicit none
        type(contr_spec_t), intent(in):: cspec         !in: general (sub)tensor contraction specification
        type(contr_task_t), intent(inout):: contr_task !inout: specific (sub)tensor contraction task
        integer(INTD), intent(inout):: ierr            !out: error code (0:success)
        integer(INTL):: fl
        integer(INTD):: i,j,n

        ierr=0; fl=1_INTL
        n=contr_task%left_arg%rank
        if(n.ge.0) then
         do i=1,n
          fl=fl*int(contr_task%left_arg%dims(i),INTL) !count all dimensions of the left tensor argument (contr + left)
         enddo
        else
         ierr=1; return
        endif
        n=contr_task%right_arg%rank
        if(n.ge.0) then
         do i=1,n
          j=cspec%rprmn(i)
          if(j.ge.1.and.j.le.n) then
           if(j.gt.cspec%ndims_contr) fl=fl*int(contr_task%right_arg%dims(i),INTL) !count only uncontracted dims of the right tensor
          else
           ierr=2; return
          endif
         enddo
         contr_task%flops=fl !number of Flops required to perform this contraction task
         contr_task%task_stat=TASK_SET
        else
         ierr=3; return
        endif
        return
        end subroutine dil_contr_task_set_flops
!-----------------------------------------------------
        subroutine dil_contr_task_copy(obji,objo,ierr) !SERIAL
        implicit none
        type(contr_task_t), intent(in):: obji         !in: original object
        type(contr_task_t), intent(out):: objo        !out: clone
        integer(INTD), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INTD):: errc

        errc=0
        call dil_subtensor_copy(obji%dest_arg,objo%dest_arg,errc)
        if(errc.eq.0) then
         call dil_subtensor_copy(obji%left_arg,objo%left_arg,errc)
         if(errc.eq.0) then
          call dil_subtensor_copy(obji%right_arg,objo%right_arg,errc)
          if(errc.eq.0) then
           objo%prefac=obji%prefac
           objo%flops=obji%flops
           objo%dev_kind=obji%dev_kind
           objo%dev_id=obji%dev_id
           objo%task_stat=obji%task_stat
          else
           errc=1
          endif
         else
          errc=2
         endif
        else
         errc=3
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine dil_contr_task_copy
!----------------------------------------------------------------
        subroutine dil_contr_task_print(cspec,contr_task,dev_out) !SERIAL
!This subroutine prints a tensor contraction task.
        implicit none
        type(contr_spec_t), intent(in):: cspec        !in: general tensor contraction specification
        type(contr_task_t), intent(in):: contr_task   !in: tensor contraction task
        integer(INTD), intent(in), optional:: dev_out !in: output device (default is 6: screen)
        character(1024):: task_str
        integer(INTD):: i,l,dvo

        l=0
        task_str(l+1:l+1)='D'; l=l+1
        call dil_subtensor_print(contr_task%dest_arg,task_str(l+1:),i); l=l+i
        task_str(l+1:l+2)='=L'; l=l+2
        call dil_subtensor_print(contr_task%left_arg,task_str(l+1:),i); l=l+i
        task_str(l+1:l+2)='*R'; l=l+2
        call dil_subtensor_print(contr_task%right_arg,task_str(l+1:),i); l=l+i
        if(present(dev_out)) then; dvo=dev_out; else; dvo=CONS_OUT; endif
        call printf(task_str(1:l),dvo)
        return
        end subroutine dil_contr_task_print
!-------------------------------------------------------------
        subroutine dil_contr_task_list_destroy(task_list,ierr) !SERIAL
!This subroutine destroys a tensor contraction task list.
        implicit none
        type(contr_task_list_t), intent(inout):: task_list !inout: (sub)tensor contraction task list
        integer(INTD), intent(inout):: ierr                !out: error code (0:success)

        ierr=0
        if(associated(task_list%contr_tasks)) then
         if(task_list%alloc) then
          deallocate(task_list%contr_tasks)
         else
          nullify(task_list%contr_tasks)
         endif
        endif
        if(allocated(task_list%task_order)) deallocate(task_list%task_order)
        task_list%alloc=.false.; task_list%reordered=.false.; task_list%num_tasks=0
        return
        end subroutine dil_contr_task_list_destroy
!-------------------------------------------------------------------
        subroutine dil_contr_task_list_create(task_list,length,ierr) !SERIAL
!This subroutine creates a tensor contraction tasks list.
        implicit none
        type(contr_task_list_t), intent(inout):: task_list !out: tensor contraction task list
        integer(INTD), intent(in):: length                 !in: length of the list
        integer(INTD), intent(inout):: ierr                !out: error code (0:success)
        integer(INTD):: i

        ierr=0
        if(length.gt.0) then
         if(associated(task_list%contr_tasks)) then
          if(size(task_list%contr_tasks).lt.length) then
           call dil_contr_task_list_destroy(task_list,ierr); if(ierr.ne.0) then; ierr=1; return; endif
          endif
         endif
         if(.not.associated(task_list%contr_tasks)) then
          allocate(task_list%contr_tasks(1:length),STAT=ierr); if(ierr.ne.0) then; ierr=2; return; endif
          task_list%alloc=.true.
         endif
         do i=1,length; task_list%contr_tasks(i)%task_stat=TASK_NULL; enddo
         task_list%num_tasks=length; task_list%reordered=.false.
        else
         ierr=3
        endif
        return
        end subroutine dil_contr_task_list_create
!-------------------------------------------------------------------
        subroutine dil_contr_task_list_shuffle(task_list,ierr,shift) !SERIAL
!This subroutine reshuffles a list of tensor contraction tasks
!in order to minimize simultaneous data accessess coming from different MPI ranks.
        implicit none
        type(contr_task_list_t), intent(inout):: task_list !inout: tensor contraction task list
        integer(INTD), intent(inout):: ierr                !out: error code (0:success)
        integer(INTD), intent(in), optional:: shift        !in: if present, the task list will be simply shifted
        integer(INTD):: shf,i,j

        ierr=0
        if(task_list%num_tasks.gt.0.and.associated(task_list%contr_tasks)) then
         if(present(shift)) then !simple shift
          shf=mod(shift,task_list%num_tasks)
          if(shf.ne.0) then
           if(allocated(task_list%task_order)) then
            if(size(task_list%task_order).lt.task_list%num_tasks) then
             deallocate(task_list%task_order)
             allocate(task_list%task_order(1:task_list%num_tasks),STAT=ierr); if(ierr.ne.0) ierr=1
            endif
           else
            allocate(task_list%task_order(1:task_list%num_tasks),STAT=ierr); if(ierr.ne.0) ierr=2
           endif
           if(ierr.eq.0) then
            if(shf.gt.0) then
             do i=1,task_list%num_tasks
              j=i+shf; if(j.gt.task_list%num_tasks) j=j-task_list%num_tasks
              task_list%task_order(i)=j
             enddo
            elseif(shf.lt.0) then
             do i=1,task_list%num_tasks
              j=i+shf; if(j.le.0_INTD) j=j+task_list%num_tasks
              task_list%task_order(i)=j
             enddo
            endif
            task_list%reordered=.true.
           endif
          endif
         else !random permutation
          !`Write
         endif
        else
         ierr=3
        endif
        return
        end subroutine dil_contr_task_list_shuffle
!--------------------------------------------------------------------
        subroutine dil_contr_task_list_print(cspec,task_list,ierr,fh) !SERIAL
!This subroutine prints a tensor contraction task list.
        implicit none
        type(contr_spec_t), intent(in):: cspec          !in: general tensor contraction specification
        type(contr_task_list_t), intent(in):: task_list !in: tensor contraction task list
        integer(INTD), intent(inout):: ierr             !out: error code (0:success)
        integer(INTD), intent(in), optional:: fh        !in: file handle (default is 6: screen)
        integer(INTD):: dev_out,i,j

        ierr=0
        if(present(fh)) then; dev_out=fh; else; dev_out=CONS_OUT; endif
        write(dev_out,'("#Printing a tensor contraction task list:")')
        if(task_list%num_tasks.gt.0) then
         if(task_list%reordered) then
          do i=1,task_list%num_tasks
           j=task_list%task_order(i)
           call dil_contr_task_print(cspec,task_list%contr_tasks(j),dev_out)
          enddo
         else
          do i=1,task_list%num_tasks
           call dil_contr_task_print(cspec,task_list%contr_tasks(i),dev_out)
          enddo
         endif
         write(dev_out,'("#End of printing.")')
        else
         write(dev_out,'("#Tensor contraction task list is empty!")')
         ierr=DIL_NO_WORK
        endif
        return
        end subroutine dil_contr_task_list_print
!---------------------------------------------------------------
        integer(INTD) function dil_dev_num(dev_kind,dev_id,ierr) !SERIAL
!Given a device kind and device id within its kind,
!this function returns a flat (uniform) device number.
        implicit none
        integer(INTD), intent(in):: dev_kind           !in: device kind
        integer(INTD), intent(in):: dev_id             !in: device id within its kind (0:...)
        integer(INTD), intent(inout), optional:: ierr  !out: error code (0:success)
        integer(INTD):: errc

        errc=0
        select case(dev_kind)
        case(DEV_HOST_CPU)
         if(dev_id.eq.0) then
          dil_dev_num=0
         else
          errc=1
         endif
        case(DEV_NVIDIA_GPU)
         if(dev_id.ge.0.and.dev_id.lt.MAX_GPUS) then
          dil_dev_num=1+dev_id
         else
          errc=2
         endif
        case(DEV_INTEL_MIC)
         if(dev_id.ge.0.and.dev_id.lt.MAX_MICS) then
          dil_dev_num=(1+MAX_GPUS)+dev_id
         else
          errc=3
         endif
        case default
         errc=4
        end select
        if(present(ierr)) ierr=errc
        return
        end function dil_dev_num
!------------------------------------------------------------
        subroutine dil_dev_kind(dev_num,dev_kind,dev_id,ierr) !SERIAL
!Given a flat device number, this subroutine returns device kind
!and device id within its kind (0..max).
        implicit none
        integer(INTD), intent(in):: dev_num   !in: flat device number
        integer(INTD), intent(out):: dev_kind !out: device kind
        integer(INTD), intent(out):: dev_id   !out: device id (within its kind)
        integer(INTD), intent(inout):: ierr   !out: error code (0:success)

        ierr=0
        if(dev_num.eq.0) then
         dev_kind=DEV_HOST_CPU; dev_id=0
        elseif(dev_num.ge.1.and.dev_num.le.MAX_GPUS) then
         dev_kind=DEV_NVIDIA_GPU; dev_id=dev_num-1
        elseif(dev_num.ge.(1+MAX_GPUS).and.dev_num.le.(MAX_GPUS+MAX_MICS)) then
         dev_kind=DEV_INTEL_MIC; dev_id=dev_num-(1+MAX_GPUS)
        else
         ierr=1
        endif
        return
        end subroutine dil_dev_kind
!---------------------------------------------------
        subroutine dil_dev_buf_destroy(dev_buf,ierr) !SERIAL
!This subroutine destroys a device buffer.
        implicit none
        type(dev_buf_t), intent(inout):: dev_buf !out: device buffer
        integer(INTD), intent(inout):: ierr      !out: error code (0:success)
        integer(INTD):: i

        ierr=0
        if(dev_buf%num_bufs.gt.0.and.dev_buf%num_bufs.le.BUFS_PER_DEV) then
         do i=1,dev_buf%num_bufs
          dev_buf%arg_buf(i)%buf_vol=0_INTL
          nullify(dev_buf%arg_buf(i)%buf_ptr) !buffers were associated to heap regions
         enddo
        else
         ierr=1
        endif
        end subroutine dil_dev_buf_destroy
!---------------------------------------------------------
        subroutine dil_arg_buf_clean(arg_buf,ierr,arg_vol) !PARALLEL (OMP)
!This subroutine zeroes out an argument buffer (or its part).
        implicit none
        type(arg_buf_t), intent(inout):: arg_buf      !out: argument buffer
        integer(INTD), intent(inout):: ierr           !out: error code (0:success)
        integer(INTL), intent(in), optional:: arg_vol !in: number of elements to zero out
        integer(INTL):: i0,i1,ivol

        ierr=0
        if(arg_buf%buf_vol.gt.0_INTL) then
         if(present(arg_vol)) then; ivol=arg_vol; else; ivol=arg_buf%buf_vol; endif
         if(ivol.ge.0_INTL.and.ivol.le.arg_buf%buf_vol) then
          if(associated(arg_buf%buf_ptr)) then
           if(size(arg_buf%buf_ptr).eq.arg_buf%buf_vol) then
            i1=lbound(arg_buf%buf_ptr,1)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i0) SCHEDULE(GUIDED)
            do i0=i1,i1+ivol-1_INTL
             arg_buf%buf_ptr(i0)=0E0_tensor_dp
            enddo
!$OMP END PARALLEL DO
           else
            ierr=1
           endif
          else
           ierr=2
          endif
         else
          ierr=3
         endif
        else
         ierr=4
        endif
        return
        end subroutine dil_arg_buf_clean
!------------------------------------------------
        subroutine dil_set_alloc_type(alloc_type) !SERIAL
        implicit none
        integer(INTD), intent(in):: alloc_type !in: Allocation type (see DIL_ALLOC_XXX parameters on top)
        select case(alloc_type)
        case(DIL_ALLOC_BASIC,DIL_ALLOC_PINNED,DIL_ALLOC_MPI) !only these are allowed
         DIL_ALLOC_TYPE=alloc_type
        end select
        return
        end subroutine dil_set_alloc_type
!------------------------------------------------------------
        subroutine cpu_ptr_alloc_r(arr,nelems,ierr,base,attr) !SERIAL
!This subroutine allocates memory for a 1d pointer array (default real).
        implicit none
        real(tensor_dp), pointer, contiguous, intent(inout):: arr(:) !out: 1d array
        integer(INTL), intent(in):: nelems                       !in: number of elements to allocate
        integer(INTD), intent(inout):: ierr                      !out: error code (0:success)
        integer(INTL), intent(in), optional:: base               !in: index numeration start offset (default is 1)
        integer(INTD), intent(in), optional:: attr               !in: attributes (pinned, etc.)
        integer(INTL):: bs
        integer(INTD):: flags
        integer(C_SIZE_T):: csize
        type(C_PTR):: caddr
        real(tensor_dp), pointer, contiguous:: fptr(:)
        real(tensor_dp):: val
        integer(MPI_ADDRESS_KIND):: mpi_size
        integer(tensor_mpi_kind):: mpi_err

        ierr=0
        if(nelems.gt.0_INTL) then
         if(present(attr)) then; flags=attr; else; flags=DIL_ALLOC_TYPE; endif
         if(present(base)) then; bs=base; else; bs=1_INTL; endif
         select case(flags)
 !Basic malloc:
         case(DIL_ALLOC_BASIC)
          allocate(arr(bs:bs+nelems-1_INTL),STAT=ierr)
          if(ierr.ne.0) then
           if(VERBOSE) write(CONS_OUT,'("#ERROR(tensor_algebra_dil::cpu_ptr_alloc_r): memory allocation failed: ",i11)') ierr
           ierr=1
          endif
 !Pinned malloc:
         case(DIL_ALLOC_PINNED)
          val=0E0_tensor_dp; csize=int(nelems*sizeof(val),C_SIZE_T); caddr=C_NULL_PTR
          !`Write (call C wrapper for cudaMallocHost)
          call c_f_pointer(caddr,fptr,[nelems]); arr(bs:)=>fptr; nullify(fptr)
 !MPI malloc:
         case(DIL_ALLOC_MPI)
          caddr=C_NULL_PTR
          val=0E0_tensor_dp; mpi_size=int(nelems*sizeof(val),MPI_ADDRESS_KIND)
          call MPI_ALLOC_MEM(mpi_size,MPI_INFO_NULL,caddr,mpi_err)
          if(mpi_err.eq.0) then
           call c_f_pointer(caddr,fptr,[nelems]); arr(bs:)=>fptr; nullify(fptr)
          else
           if(VERBOSE) write(CONS_OUT,'("#ERROR(tensor_algebra_dil::cpu_ptr_alloc_r): MPI memory allocation failed: ",i11)')&
           &mpi_err
           ierr=3
          endif
         case default
          if(VERBOSE) write(CONS_OUT,'("#ERROR(tensor_algebra_dil::cpu_ptr_alloc_r): invalid allocation attributes: ",i11)') flags
          ierr=4
         end select
        else
         ierr=5
        endif
        return
        end subroutine cpu_ptr_alloc_r
!-----------------------------------------------
        subroutine cpu_ptr_free_r(arr,ierr,attr) !SERIAL
!This subroutine deallocates memory assigned to a 1d pointer array (default real).
!WARNING: This subroutine assumes that MPI_FREE_MEM has an explicitly declared interface in MPI.mod!
!If not, it may NOT work correctly!!!
        implicit none
        real(tensor_dp), pointer, contiguous, intent(inout):: arr(:) !inout: 1d array
        integer(INTD), intent(inout):: ierr                      !out: error code (0:success)
        integer(INTD), intent(in), optional:: attr               !in: attributes (pinned, etc.)
        type(C_PTR):: caddr
        integer(INTD):: flags
        integer(tensor_mpi_kind):: errc

!#ifdef PROTO_MPI_FREE_MEM
!        interface
!         subroutine MPI_FREE_MEM(base_ptr,ierr)
!          import
!          type(C_PTR), value:: base_ptr
!          integer(INTD), intent(out):: ierr
!         end subroutine MPI_FREE_MEM
!        end interface
!#endif

        ierr=0
        if(associated(arr)) then
         if(present(attr)) then; flags=attr; else; flags=DIL_ALLOC_TYPE; endif
         select case(flags)
         case(DIL_ALLOC_NOT)
          ierr=1
         case(DIL_ALLOC_BASIC)
          deallocate(arr,STAT=ierr)
          if(ierr.ne.0) ierr=2
         case(DIL_ALLOC_PINNED)
          caddr=c_loc(arr(1));
          !`Write (call C wrapper for cudaMallocHost)
          if(ierr.ne.0) ierr=3
          nullify(arr)
         case(DIL_ALLOC_MPI)
          errc=0
!         caddr=c_loc(arr); call MPI_FREE_MEM(caddr,errc) !<caddr> must be passed by value!!!
          call MPI_FREE_MEM(arr,errc)
          if(errc.ne.0) ierr=4
          nullify(arr)
         case(DIL_ALLOC_EXT)
          nullify(arr)
         case default
          ierr=5
         end select
        else
         ierr=6
        endif
        return
        end subroutine cpu_ptr_free_r
!--------------------------------------------------
        integer(INTD) function my_mpi_size(my_comm) !SERIAL
!Returns the rank of an MPI process.
        implicit none
        integer(tensor_mpi_kind), intent(in), optional:: my_comm !in: MPI communicator
        integer(tensor_mpi_kind):: i,ierr
        i=0
        if(present(my_comm)) then !explicit communicator
         call MPI_COMM_SIZE(my_comm,i,ierr); if(ierr.ne.0) i=-1
        else !default communicator
         call MPI_COMM_SIZE(MPI_COMM_WORLD,i,ierr); if(ierr.ne.0) i=-1
        endif
        my_mpi_size=i
        return
        end function my_mpi_size
!--------------------------------------------------
        integer(INTD) function my_mpi_rank(my_comm) !SERIAL
!Returns the rank of an MPI process.
        implicit none
        integer(tensor_mpi_kind), intent(in), optional:: my_comm !in: MPI communicator
        integer(tensor_mpi_kind):: i,ierr
        i=-1
        if(present(my_comm)) then !explicit communicator
         call MPI_COMM_RANK(my_comm,i,ierr); if(ierr.ne.0) i=-1
        else !default communicator
         call MPI_COMM_RANK(MPI_COMM_WORLD,i,ierr); if(ierr.ne.0) i=-1
        endif
        my_mpi_rank=i
        return
        end function my_mpi_rank
!---------------------------------------------------------------------------
#ifdef USE_MIC
!DIR$ ATTRIBUTES OFFLOAD:mic:: divide_segment_i8
#endif
        subroutine divide_segment_i8(seg_range,subseg_num,subseg_sizes,ierr) !SERIAL
!A segment of range <seg_range> will be divided into <subseg_num> subsegments maximally uniformly.
!The length of each subsegment will be returned in the array <subseg_sizes(1:subseg_num)>.
!Any two subsegments will not differ in length by more than 1, longer subsegments preceding the shorter ones.
        implicit none
        integer(INTL), intent(in):: seg_range,subseg_num
        integer(INTL), intent(out):: subseg_sizes(1:subseg_num)
        integer(INTD), intent(inout):: ierr
        integer(INTL) i,j,k,l,m,n
        ierr=0_INTD
        if(seg_range.gt.0_INTL.and.subseg_num.gt.0_INTL) then
         n=seg_range/subseg_num; m=mod(seg_range,subseg_num)
         do i=1_INTL,m; subseg_sizes(i)=n+1_INTL; enddo
         do i=m+1_INTL,subseg_num; subseg_sizes(i)=n; enddo
        else
         ierr=-1_INTD
        endif
        return
        end subroutine divide_segment_i8
!-------------------------------------------------
        subroutine merge_sort_key_int4(ni,key,trn) !SERIAL
!This subroutine sorts an array of NI items in a non-descending order according to their keys.
!The algorithm is due to Johann von Neumann, implemented by Dmitry I. Lyakh.
!INPUT:
! - ni - number of items;
! - key(1:ni) - item keys (retrieved by old item numbers!): arbitrary integers;
! - trn(0:ni) - initial permutation of ni items (a sequence of old numbers), trn(0) is the initial sign;
!OUTPUT:
! - trn(0:ni) - sorted permutation (new sequence of old numbers) of ni items (according to their keys), the sign is in trn(0).
        implicit none
        integer(INTD), intent(in):: ni,key(1:ni)
        integer(INTD), intent(inout):: trn(0:ni)
        integer(INTD), parameter:: max_in_mem=1024
        integer(INTD):: i,j,k,l,m,n,k1,k2,k3,k4,kf
        integer(INTD), target:: prms(1:max_in_mem)
        integer(INTD), pointer:: prm(:)

        if(ni.gt.1) then
         if(ni.le.max_in_mem) then; prm=>prms; else; allocate(prm(1:ni)); endif
         n=1
         do while(n.lt.ni)
          m=n*2
          do i=1,ni,m
           k1=i; k2=i+n
           if(k2.gt.ni) then
            k2=ni+1; k3=0; k4=0 !no right block, only left block
           else
            k3=i+n; k4=min(ni+1,i+m) !right block present
           endif
           kf=min(ni+1,i+m)-i; l=0
           do while(l.lt.kf)
            if(k3.ge.k4) then !right block is over
             prm(i+l:i+kf-1)=trn(k1:k2-1); l=kf
            elseif(k1.ge.k2) then !left block is over
             prm(i+l:i+kf-1)=trn(k3:k4-1); l=kf
            else
             if(key(trn(k1))-key(trn(k3)).gt.0) then
              prm(i+l)=trn(k3); k3=k3+1; trn(0)=(1-2*mod(k2-k1,2))*trn(0)
             else
              prm(i+l)=trn(k1); k1=k1+1
             endif
             l=l+1
            endif
           enddo
          enddo
          trn(1:ni)=prm(1:ni)
          n=m
         enddo
         if(ni.le.max_in_mem) then; nullify(prm); else; deallocate(prm); endif
        endif
        return
        end subroutine merge_sort_key_int4
!-------------------------------------------------
        subroutine merge_sort_key_int8(ni,key,trn) !SERIAL
!This subroutine sorts an array of NI items in a non-descending order according to their keys.
!The algorithm is due to Johann von Neumann, implemented by Dmitry I. Lyakh.
!INPUT:
! - ni - number of items;
! - key(1:ni) - item keys (retrieved by old item numbers!): arbitrary integers;
! - trn(0:ni) - initial permutation of ni items (a sequence of old numbers), trn(0) is the initial sign;
!OUTPUT:
! - trn(0:ni) - sorted permutation (new sequence of old numbers) of ni items (according to their keys), the sign is in trn(0).
        implicit none
        integer(INTL), intent(in):: ni,key(1:ni)
        integer(INTL), intent(inout):: trn(0:ni)
        integer(INTL), parameter:: max_in_mem=1024
        integer(INTL):: i,j,k,l,m,n,k1,k2,k3,k4,kf
        integer(INTL), target:: prms(1:max_in_mem)
        integer(INTL), pointer:: prm(:)

        if(ni.gt.1) then
         if(ni.le.max_in_mem) then; prm=>prms; else; allocate(prm(1:ni)); endif
         n=1
         do while(n.lt.ni)
          m=n*2
          do i=1,ni,m
           k1=i; k2=i+n
           if(k2.gt.ni) then
            k2=ni+1; k3=0; k4=0 !no right block, only left block
           else
            k3=i+n; k4=min(ni+1,i+m) !right block present
           endif
           kf=min(ni+1,i+m)-i; l=0
           do while(l.lt.kf)
            if(k3.ge.k4) then !right block is over
             prm(i+l:i+kf-1)=trn(k1:k2-1); l=kf
            elseif(k1.ge.k2) then !left block is over
             prm(i+l:i+kf-1)=trn(k3:k4-1); l=kf
            else
             if(key(trn(k1))-key(trn(k3)).gt.0) then
              prm(i+l)=trn(k3); k3=k3+1; trn(0)=(1-2*mod(k2-k1,2))*trn(0)
             else
              prm(i+l)=trn(k1); k1=k1+1
             endif
             l=l+1
            endif
           enddo
          enddo
          trn(1:ni)=prm(1:ni)
          n=m
         enddo
         if(ni.le.max_in_mem) then; nullify(prm); else; deallocate(prm); endif
        endif
        return
        end subroutine merge_sort_key_int8
!----------------------------------------------
        subroutine merge_sort_real8(ni,trn,dir)
!This subroutine sorts an array of NI items.
!The algorithm is due to Johann von Neumann.
!INPUT:
! - ni - number of items;
! - trn(1:ni) - items (arbitrary real*8 numbers);
! - dir - direction of sort (negative - descending, positive - ascending order);
!OUTPUT:
! - trn(0:ni) - sorted items.
        implicit none
        integer(INTL), intent(in):: ni
        real(8), intent(inout):: trn(1:ni)
        integer(INTD), intent(in), optional:: dir
        integer(INTL):: i,j,k,l,m,n,k1,k2,k3,k4,kf
        real(8), allocatable:: prm(:)
        real(8):: ds

        if(ni.gt.1) then
         ds=+1d0
         if(present(dir)) then
          if(dir.lt.0) ds=-1d0
         endif
         allocate(prm(1:ni))
         n=1
         do while(n.lt.ni)
          m=n*2
          do i=1,ni,m
           k1=i; k2=i+n
           if(k2.gt.ni) then
            k2=ni+1; k3=0; k4=0 !no right block, only left block
           else
            k3=i+n; k4=min(ni+1,i+m) !right block present
           endif
           kf=min(ni+1,i+m)-i; l=0
           do while(l.lt.kf)
            if(k3.ge.k4) then !right block is over
             prm(i+l:i+kf-1)=trn(k1:k2-1); l=kf
            elseif(k1.ge.k2) then !left block is over
             prm(i+l:i+kf-1)=trn(k3:k4-1); l=kf
            else
             if((trn(k1)-trn(k3))*ds.gt.0d0) then
              prm(i+l)=trn(k3); k3=k3+1
             else
              prm(i+l)=trn(k1); k1=k1+1
             endif
             l=l+1
            endif
           enddo
          enddo
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
          do i=1,ni
           trn(i)=prm(i)
          enddo
!$OMP END PARALLEL DO
          n=m
         enddo
         deallocate(prm)
        endif
        return
        end subroutine merge_sort_real8
!----------------------------------------------------------------------
        subroutine str_parse(str,sep,nwords,words,ierr,str_len,sep_len) !SERIAL
!This subroutine extracts separable "words" from a string.
!It allows for multiple word separators, but all of them must have the same length!
        implicit none
        character(*), intent(in):: str                !in: string of words separated by allowed separators
        character(*), intent(in):: sep                !in: concatenated allowed word separators (default length of a separator is 1)
        integer(INTD), intent(out):: nwords           !out: number of words found in the string
        integer(INTD), intent(inout):: words(2,*)     !out: beginning (1) and end (2) position of each word in the string
        integer(INTD), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INTD), intent(in), optional:: str_len !in: string length (default will use the entire <str>)
        integer(INTD), intent(in), optional:: sep_len !in: individual separator length (all separators must have the same length)
        integer(INTD):: k,l,errc,strl,sepl,spl

        errc=0; nwords=0
        sepl=len(sep) !total length of the string containing allowed separators
        if(present(str_len)) then; strl=str_len; else; strl=len(str); endif !length of the analyzed string
        if(present(sep_len)) then; spl=sep_len; else; spl=1; endif !length of each individual separator
        if(strl.gt.0.and.spl.gt.0.and.sepl.ge.spl.and.mod(sepl,spl).eq.0) then
         k=0; l=1
         do while(l.le.strl)
          if(this_is_separator(l)) then !separator found
           if(k.gt.0) then; nwords=nwords+1; words(1:2,nwords)=(/k,l-1_INTD/); endif
           k=0; l=l+spl
          else
           if(k.eq.0) k=l !beginning of a word
           l=l+1
          endif
         enddo
         if(k.gt.0) then; nwords=nwords+1; words(1:2,nwords)=(/k,strl/); endif
        else
         errc=1
        endif
        if(present(ierr)) ierr=errc
        return
        contains

         logical function this_is_separator(spos)
         integer(INTD), intent(in):: spos !examined position in str
         integer(INTD):: j0,j1,j2
         this_is_separator=.false.
         if(strl-spos+1.ge.spl) then
          j1=0; j2=0
          do j0=1,sepl
           if(str(spos+j1:spos+j1).eq.sep(j0:j0)) j2=j2+1
           j1=j1+1
           if(j1.eq.spl) then
            if(j2.eq.spl) then; this_is_separator=.true.; return; endif
            j1=0; j2=0
           endif
          enddo
         endif
         return
         end function this_is_separator

        end subroutine str_parse
!-------------------------------------------------
        function str2int(str,sl,ierr) result(intg) !SERIAL
!This function converts an integer number given as a string into an integer.
        implicit none
        character(*), intent(in):: str          !in: string containing an integer
        integer(INTD), intent(in):: sl          !in: string length: str[1:sl]
        integer(INTD):: intg                    !out: integer
        integer(INTD), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INTD):: i,j,n,i0,i9,ie,sgn

        ie=0; intg=0
        if(sl.gt.0) then
         i=1; do while(str(i:i).eq.' '); i=i+1; if(i.gt.sl) exit; enddo
         if(i.le.sl) then
          sgn=+1
          if(str(i:i).eq.'-') then
           sgn=-1; i=i+1
          elseif(str(i:i).eq.'+') then
           sgn=+1; i=i+1
          endif
          n=0; i0=iachar('0'); i9=iachar('9')
          do while(i.le.sl)
           if(str(i:i).ne.' ') then
            j=iachar(str(i:i))
            if(j.ge.i0.and.j.le.i9) then
             intg=intg*10+(j-i0); n=n+1
            else
             ie=1; exit
            endif
           endif
           i=i+1
          enddo
          if(ie.eq.0) then
           if(n.gt.0) then
            intg=intg*sgn
           else
            ie=2
           endif
          endif
         else
          ie=3
         endif
        else
         ie=4
        endif
        if(present(ierr)) ierr=ie
        return
        end function str2int
!-----------------------------------------
        subroutine int2str_i4(intg,str,sl) !SERIAL
!This subroutine converts an integer into a string.
        implicit none
        integer(4), intent(in):: intg     !in: integer
        character(*), intent(inout):: str !out: string
        integer(INTD), intent(out):: sl   !out: string length
        character(1):: ch
        integer(4):: i,j,i0

        sl=0
        if(intg.eq.0) then
         sl=sl+1; str(sl:sl)='0'
        else
         if(intg.lt.0) then; sl=sl+1; str(sl:sl)='-'; endif
         i=abs(intg); i0=iachar('0')
         do while(i.gt.0)
          j=mod(i,10); i=i/10
          sl=sl+1; str(sl:sl)=achar(i0+j)
         enddo
         do i=1,sl/2
          j=sl+1-i; ch=str(i:i); str(i:i)=str(j:j); str(j:j)=ch
         enddo
        endif
        return
        end subroutine int2str_i4
!-----------------------------------------
        subroutine int2str_i8(intg,str,sl) !SERIAL
!This subroutine converts an integer into a string.
        implicit none
        integer(8), intent(in):: intg     !in: integer
        character(*), intent(inout):: str !out: string
        integer(INTD), intent(out):: sl   !out: string length
        character(1):: ch
        integer(8):: i
        integer(INTD):: j,i0

        sl=0
        if(intg.eq.0_8) then
         sl=sl+1; str(sl:sl)='0'
        else
         if(intg.lt.0_8) then; sl=sl+1; str(sl:sl)='-'; endif
         i=abs(intg); i0=iachar('0')
         do while(i.gt.0_8)
          j=mod(i,10_8); i=i/10_8
          sl=sl+1; str(sl:sl)=achar(i0+j)
         enddo
         do i0=1,sl/2
          j=sl+1-i0; ch=str(i0:i0); str(i0:i0)=str(j:j); str(j:j)=ch
         enddo
        endif
        return
        end subroutine int2str_i8
!------------------------------------
        subroutine printf(str,fh,adv) !SERIAL
!This subroutine prints a string str(1:*) into file#FH
        implicit none
        character(*), intent(in):: str            !in: string to print
        integer(INTD), intent(in), optional:: fh  !file number (6: screen)
        logical, intent(in), optional:: adv       !carriage return control (.false. - no carriage return)
        character(32):: frm
        integer(INTD):: l,k,od

        if(present(fh)) then; od=fh; else; od=CONS_OUT; endif
        l=len_trim(str)
        if(l.gt.0) then
         frm(1:2)='(A'
         call int2str(l,frm(3:32),k)
         k=3+k; frm(k:k)=')'
         if(present(adv)) then
          if(adv) then
           write(od,frm(1:k)) str(1:l)
          else
           write(od,frm(1:k),advance='no') str(1:l)
          endif
         else
          write(od,frm(1:k)) str(1:l)
         endif
        else
         if(present(adv)) then
          if(adv) write(od,*)
         endif
        endif
        return
        end subroutine printf
!-------------------------------------------------
        subroutine dil_rank_window_clean(rwc,ierr) !SERIAL
!This function cleans a rank/window container.
        implicit none
        type(rank_win_cont_t), intent(inout):: rwc !inout: rank_window container
        integer(INTD), intent(inout), optional:: ierr !out: error code (0:success)
        if(present(ierr)) ierr=0
        rwc%num_entries=0; rwc%first_entry=-1
        return
        end subroutine dil_rank_window_clean
!----------------------------------------------------------------
        logical function dil_rank_window_new(rwc,nrank,nwin,ierr) !SERIAL
!This function looks up a given pair {rank;window} in a rank_window container.
!If it is already there, it returns .FALSE. If not, it appends it and returns .TRUE.
        implicit none
        type(rank_win_cont_t), intent(inout):: rwc !inout: rank_window container
        integer(INTD), intent(in):: nrank          !in: MPI rank
        integer, intent(in):: nwin                 !in: MPI Window
        integer(INTD), intent(inout):: ierr        !out: error code (0:success)
        integer(INTD):: k,l,m,n

        ierr=0; dil_rank_window_new=.false.
        if(rwc%num_entries.gt.0) then
         n=rwc%first_entry; m=-1
         do while(n.gt.0)
          l=rwc%rw_entry(n)%window
          if(nwin.le.l) exit
          m=n; n=rwc%next_win(n)
         enddo
         if(nwin.eq.l) then
          m=-1
          do while(n.gt.0)
           k=rwc%rw_entry(n)%rank
           if(nrank.le.k) exit
           m=n; n=rwc%next_rank(n)
          enddo
          if(nrank.eq.k) then !entry already exists (return .false.)
           return
          else !nrank<k or end of sublist: append
           if(rwc%num_entries.ge.MAX_TILES_PER_PART) then
            if(VERBOSE) write(CONS_OUT,'("#ERROR(tensor_algebra_dil::dil_rank_window_clean): MAX_TILES_PER_PART exceeded: ",i9)')&
            &MAX_TILES_PER_PART
            ierr=1; return
           endif
           rwc%num_entries=rwc%num_entries+1
           rwc%rw_entry(rwc%num_entries)%rank=nrank; rwc%rw_entry(rwc%num_entries)%window=nwin
           if(m.gt.0) then
            rwc%next_rank(m)=rwc%num_entries; rwc%next_win(rwc%num_entries)=rwc%next_win(m)
           else
            rwc%first_entry=rwc%num_entries; rwc%next_win(rwc%num_entries)=rwc%next_win(n)
           endif
           rwc%next_rank(rwc%num_entries)=n
           dil_rank_window_new=.true.
          endif
         else !nwin<l or end of list: append
          if(rwc%num_entries.ge.MAX_TILES_PER_PART) then
           if(VERBOSE) write(CONS_OUT,'("#ERROR(tensor_algebra_dil::dil_rank_window_clean): MAX_TILES_PER_PART exceeded: ",i9)')&
           &MAX_TILES_PER_PART
           ierr=2; return
          endif
          rwc%num_entries=rwc%num_entries+1
          rwc%rw_entry(rwc%num_entries)%rank=nrank; rwc%rw_entry(rwc%num_entries)%window=nwin
          if(m.gt.0) then; rwc%next_win(m)=rwc%num_entries; else; rwc%first_entry=rwc%num_entries; endif
          rwc%next_win(rwc%num_entries)=n; rwc%next_rank(rwc%num_entries)=-1
          dil_rank_window_new=.true.
         endif
        else
         rwc%num_entries=1; rwc%first_entry=1
         rwc%rw_entry(1)%rank=nrank; rwc%rw_entry(1)%window=nwin
         rwc%next_win(1)=-1; rwc%next_rank(1)=-1
         dil_rank_window_new=.true.
        endif
        return
        end function dil_rank_window_new
!-----------------------------------------------------------------------------
        subroutine dil_get_contr_pattern(tcs,tcl,nd,nl,nr,dprm,lprm,rprm,ierr) !SERIAL
!This subroutine reads a mnemonic tensor contraction specification and
!generates the tensor argument ranks and matricization permutations.
        implicit none
        character(*), intent(in):: tcs            !in: tensor contraction specification string
        integer(INTD), intent(in):: tcl           !in: length of the tensor contraction specification string
        integer(INTD), intent(out):: nd           !out: destination tensor rank
        integer(INTD), intent(out):: nl           !out: left tensor rank
        integer(INTD), intent(out):: nr           !out: right tensor rank
        integer(INTD), intent(inout):: dprm(1:*)  !out: matricization permutation for the destination tensor argument
        integer(INTD), intent(inout):: lprm(1:*)  !out: matricization permutation for the left tensor argument
        integer(INTD), intent(inout):: rprm(1:*)  !out: matricization permutation for the right tensor argument
        integer(INTD), intent(inout):: ierr       !out: error code (0:success)
        integer(INTD):: i,j,k,l,m,n,nc
        character(MAX_TENSOR_RANK):: dts,lts,rts
        integer(INTD):: trn(0:MAX_TENSOR_RANK),key(1:MAX_TENSOR_RANK),pos(1:MAX_TENSOR_RANK)

        ierr=0
        if(DIL_DEBUG) call printf('#DEBUG(DIL): Tensor contraction pattern: '//tcs(1:tcl))
        if(tcl.gt.0) then
         l=1
 !Destination tensor argument:
         nd=0
         i=index(tcs(l:tcl),'(')+(l-1); if(i.lt.l) then; ierr=1; return; endif
         if(i.lt.tcl) then
          j=index(tcs(i+1:tcl),')')+i; if(j.lt.i+1) then; ierr=2; return; endif
          do k=i+1_INTD,j-1_INTD,2_INTD
           if(.not.this_is_letter(tcs(k:k))) then; ierr=3; return; endif
           nd=nd+1; dts(nd:nd)=tcs(k:k)
           if(tcs(k+1:k+1).ne.','.and.tcs(k+1:k+1).ne.')') then; ierr=4; return; endif
          enddo
         else
          ierr=5; return
         endif
         l=j+1; if(l.gt.tcl) then; ierr=6; return; endif
 !Left tensor argument:
         nl=0
         i=index(tcs(l:tcl),'(')+(l-1); if(i.lt.l) then; ierr=7; return; endif
         if(i.lt.tcl) then
          j=index(tcs(i+1:tcl),')')+i; if(j.lt.i+1) then; ierr=8; return; endif
          do k=i+1_INTD,j-1_INTD,2_INTD
           if(.not.this_is_letter(tcs(k:k))) then; ierr=9; return; endif
           nl=nl+1; lts(nl:nl)=tcs(k:k)
           if(tcs(k+1:k+1).ne.','.and.tcs(k+1:k+1).ne.')') then; ierr=10; return; endif
          enddo
         else
          ierr=11; return
         endif
         l=j+1; if(l.gt.tcl) then; ierr=12; return; endif
 !Right tensor argument:
         nr=0
         i=index(tcs(l:tcl),'(')+(l-1); if(i.lt.l) then; ierr=13; return; endif
         if(i.lt.tcl) then
          j=index(tcs(i+1:tcl),')')+i; if(j.lt.i+1) then; ierr=14; return; endif
          do k=i+1_INTD,j-1_INTD,2_INTD
           if(.not.this_is_letter(tcs(k:k))) then; ierr=15; return; endif
           nr=nr+1; rts(nr:nr)=tcs(k:k)
           if(tcs(k+1:k+1).ne.','.and.tcs(k+1:k+1).ne.')') then; ierr=16; return; endif
          enddo
         else
          ierr=17; return
         endif
 !Determine matricization permutations:
         if(mod(nl+nr-nd,2).eq.0) then
          nc=(nl+nr-nd)/2; lprm(1:nl)=0; rprm(1:nr)=0
  !Uncontracted dimensions:
          if(nd.gt.0) then
           do i=1,nd
            if(nl.gt.0) then; j=index(lts(1:nl),dts(i:i)); else; j=0; endif
            if(j.le.0) then
             if(nr.gt.0) then; j=index(rts(1:nr),dts(i:i)); else; j=0; endif
             if(j.gt.0) then
              key(i)=1; pos(i)=j !right index
             else
              ierr=18; return
             endif
            else
             key(i)=0; pos(i)=j !left index
            endif
           enddo
           trn(0:nd)=(/+1_INTD,(i,i=1_INTD,nd)/)
           call merge_sort_key_int(nd,key,trn)
           do i=1,nd
            j=trn(i) !old index number
            dprm(j)=i
            if(key(j).eq.0) then
             lprm(pos(j))=nc+i
            else
             rprm(pos(j))=nc*2_INTD+i-nl
            endif
           enddo
          endif
  !Contracted dimensions:
          if(nc.gt.0) then
           k=0
           do i=1,nl
            if(lprm(i).le.0) then
             k=k+1; lprm(i)=k
             j=index(rts(1:nr),lts(i:i)); if(j.le.0) then; ierr=19; return; endif
             if(rprm(j).le.0) then
              rprm(j)=k
             else
              ierr=20; return
             endif
            endif
           enddo
          endif
         else
          ierr=21; return
         endif
        else
         ierr=22; return
        endif
        if(DIL_DEBUG) then
         write(CONS_OUT,'("#DEBUG(DIL): Destination matricization permutation:",32(1x,i2))') dprm(1:nd)
         write(CONS_OUT,'("#DEBUG(DIL): Left matricization permutation:",32(1x,i2))') lprm(1:nl)
         write(CONS_OUT,'("#DEBUG(DIL): Right matricization permutation:",32(1x,i2))') rprm(1:nr)
        endif
        return

        contains

         logical function this_is_letter(ch)
         character(1), intent(in):: ch
         integer(INTD):: j0
         j0=iachar(ch)
         if((j0.ge.iachar('a').and.j0.le.iachar('z')).or.(j0.ge.iachar('A').and.j0.le.iachar('Z'))) then
          this_is_letter=.true.
         else
          this_is_letter=.false.
         endif
         return
         end function this_is_letter

        end subroutine dil_get_contr_pattern
!----------------------------------------------------
#ifdef USE_MIC
!DIR$ ATTRIBUTES OFFLOAD:mic:: thread_wtime
#endif
        function thread_wtime(time_start) result(thw) !SERIAL, threadsafe
!Returns time in seconds.
        implicit none
        real(8), intent(in), optional:: time_start !in: clock start time
        real(8):: thw
#ifdef VAR_OMP
        thw=omp_get_wtime()
#else
        call cpu_time(thw)
#endif
        if(present(time_start)) thw=thw-time_start
        return
        end function thread_wtime
!-----------------------------------------------------
        function process_wtime(time_start) result(thw) !SERIAL
!Returns time in seconds.
        implicit none
        real(8), intent(in), optional:: time_start !in: clock start time
        real(8):: thw
        thw=MPI_WTIME()
        if(present(time_start)) thw=thw-time_start
        return
        end function process_wtime
!--------------------------------------------------------
#ifdef USE_MIC
!DIR$ ATTRIBUTES OFFLOAD:mic:: permutation_invert
#endif
        subroutine permutation_invert(plen,pin,pout,ierr) !SERIAL
!This subroutine returns an inverse for a given permutation.
        implicit none
        integer(INTD), intent(in):: plen          !in: length of the permutation
        integer(INTD), intent(in):: pin(1:plen)   !in: input permutation (must be valid!)
        integer(INTD), intent(out):: pout(1:plen) !out: inverse of the input permutation
        integer(INTD), intent(inout):: ierr       !out: error code (0:success)
        integer(INTD):: i

        ierr=0
        if(plen.gt.0) then
         do i=1,plen; pout(pin(i))=i; enddo
        elseif(plen.lt.0) then
         ierr=1
        endif
        return
        end subroutine permutation_invert
!-------------------------------------------------------------------
#ifdef USE_MIC
!DIR$ ATTRIBUTES OFFLOAD:mic:: permutation_trivial
#endif
        function permutation_trivial(plen,prm,ierr) result(prm_triv) !SERIAL
!This function returns .true. if the given permutation is trivial, .false. otherwise.
        implicit none
        integer(INTD), intent(in):: plen        !in: permutation length
        integer(INTD), intent(in):: prm(1:plen) !in: permutation
        integer(INTD), intent(inout):: ierr     !out: error code (0:success)
        logical:: prm_triv                      !out: result
        integer(INTD):: i

        ierr=0; prm_triv=.true.
        if(plen.ge.0) then
         do i=1,plen; if(prm(i).ne.i) then; prm_triv=.false.; exit; endif; enddo
        else
         ierr=1
        endif
        return
        end function permutation_trivial
!--------------------------------------------------------------------------------------
#ifdef USE_MIC
!DIR$ ATTRIBUTES OFFLOAD:mic:: dil_tensor_slice
#endif
        subroutine dil_tensor_slice(dim_num,tens,tens_ext,slice,slice_ext,ext_beg,ierr) !PARALLEL (OMP)
!This subroutine extracts a slice from a tensor.
        implicit none
!-------------------------------------------------
        integer(INTD), parameter:: real_kind=tensor_dp
!-------------------------------------------------
        integer(INTD), intent(in):: dim_num              !in: number of tensor dimensions
        real(real_kind), intent(in):: tens(0:*)          !in: tensor
        integer(INTD), intent(in):: tens_ext(1:dim_num)  !in: dimension extents for <tens>
        real(real_kind), intent(inout):: slice(0:*)      !out: slice
        integer(INTD), intent(in):: slice_ext(1:dim_num) !in: dimension extents for <slice>
        integer(INTD), intent(in):: ext_beg(1:dim_num)   !in: beginning dimension offsets for <tens> (numeration starts at 0)
        integer(INTD), intent(inout):: ierr              !out: error code (0:success)
        integer(INTD):: i,j,k,l,m,n,ks,kf,im(1:dim_num)
        integer(INTL):: lts,lss,l_in,l_out,lb,le,ll,bases_in(1:dim_num),bases_out(1:dim_num),segs(0:MAX_THREADS)
        real(8) time_beg
#ifdef USE_MIC
!DIR$ ATTRIBUTES OFFLOAD:mic:: real_kind
!DIR$ ATTRIBUTES ALIGN:128:: real_kind,im,bases_in,bases_out,segs
#endif
        ierr=0
!        time_beg=thread_wtime() !debug
        if(dim_num.gt.0) then
         lts=1_INTL; do i=1,dim_num; bases_in(i)=lts; lts=lts*tens_ext(i); enddo   !tensor indexing bases
         lss=1_INTL; do i=1,dim_num; bases_out(i)=lss; lss=lss*slice_ext(i); enddo !slice indexing bases
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,m,n,im,l_in,l_out,lb,le,ll)
#ifdef VAR_OMP
         n=omp_get_thread_num(); m=omp_get_num_threads()
#else
         n=0; m=1
#endif
!$OMP MASTER
         segs(0)=0_INTL; call divide_segment_i8(lss,int(m,INTL),segs(1:),ierr); do i=2,m; segs(i)=segs(i)+segs(i-1); enddo
!$OMP END MASTER
!$OMP BARRIER
!$OMP FLUSH(segs)
         l_out=segs(n); do i=dim_num,1,-1; im(i)=l_out/bases_out(i); l_out=l_out-im(i)*bases_out(i); enddo
         l_in=ext_beg(1); do i=2,dim_num; l_in=l_in+(ext_beg(i)+im(i))*bases_in(i); enddo
         lb=int(im(1),INTL); le=int(slice_ext(1)-1,INTL); l_out=segs(n)-lb
         sloop: do while(l_out+lb.lt.segs(n+1))
          le=min(le,segs(n+1)-1_INTL-l_out)
          do ll=lb,le; slice(l_out+ll)=tens(l_in+ll); enddo
          l_out=l_out+le+1_INTL; lb=0_INTL
          do i=2,dim_num
           l_in=l_in+bases_in(i); im(i)=im(i)+1
           if(im(i).ge.slice_ext(i)) then
            l_in=l_in-im(i)*bases_in(i); im(i)=0
           else
            exit
           endif
          enddo
         enddo sloop
!$OMP END PARALLEL
        else
         ierr=1 !zero-rank tensor
        endif
!        write(CONS_OUT,'("#DEBUG(dil_tensor_algebra::dil_tensor_slice): kernel time/error code: ",F10.4,1x,i3)')&
!        &thread_wtime(time_beg),ierr !debug
        return
        end subroutine dil_tensor_slice
!---------------------------------------------------------------------------------------
#ifdef USE_MIC
!DIR$ ATTRIBUTES OFFLOAD:mic:: dil_tensor_insert
#endif
        subroutine dil_tensor_insert(dim_num,tens,tens_ext,slice,slice_ext,ext_beg,ierr) !PARALLEL (OMP)
!This subroutine inserts a slice into a tensor.
        implicit none
!-------------------------------------------------
        integer(INTD), parameter:: real_kind=tensor_dp
!-------------------------------------------------
        integer(INTD), intent(in):: dim_num              !in: number of tensor dimensions
        real(real_kind), intent(inout):: tens(0:*)       !inout: tensor
        integer(INTD), intent(in):: tens_ext(1:dim_num)  !in: dimension extents for <tens>
        real(real_kind), intent(in):: slice(0:*)         !in: slice
        integer(INTD), intent(in):: slice_ext(1:dim_num) !in: dimension extents for <slice>
        integer(INTD), intent(in):: ext_beg(1:dim_num)   !in: beginning dimension offsets for <tens> (numeration starts at 0)
        integer(INTD), intent(inout):: ierr              !out: error code (0:success)
        integer(INTD):: i,j,k,l,m,n,ks,kf,im(1:dim_num)
        integer(INTL):: lts,lss,l_in,l_out,lb,le,ll,bases_in(1:dim_num),bases_out(1:dim_num),segs(0:MAX_THREADS)
        real(8) time_beg
#ifdef USE_MIC
!DIR$ ATTRIBUTES OFFLOAD:mic:: real_kind
!DIR$ ATTRIBUTES ALIGN:128:: real_kind,im,bases_in,bases_out,segs
#endif
        ierr=0
!        time_beg=thread_wtime() !debug
        if(dim_num.gt.0) then
         lts=1_INTL; do i=1,dim_num; bases_out(i)=lts; lts=lts*tens_ext(i); enddo !tensor indexing bases
         lss=1_INTL; do i=1,dim_num; bases_in(i)=lss; lss=lss*slice_ext(i); enddo !slice indexing bases
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,m,n,im,l_in,l_out,lb,le,ll)
#ifdef VAR_OMP
         n=omp_get_thread_num(); m=omp_get_num_threads()
#else
         n=0; m=1
#endif
!$OMP MASTER
         segs(0)=0_INTL; call divide_segment_i8(lss,int(m,INTL),segs(1:),ierr); do i=2,m; segs(i)=segs(i)+segs(i-1); enddo
!$OMP END MASTER
!$OMP BARRIER
!$OMP FLUSH(segs)
         l_in=segs(n); do i=dim_num,1,-1; im(i)=l_in/bases_in(i); l_in=l_in-im(i)*bases_in(i); enddo
         l_out=ext_beg(1); do i=2,dim_num; l_out=l_out+(ext_beg(i)+im(i))*bases_out(i); enddo
         lb=int(im(1),INTL); le=int(slice_ext(1)-1,INTL); l_in=segs(n)-lb
         sloop: do while(l_in+lb.lt.segs(n+1))
          le=min(le,segs(n+1)-1_INTL-l_in)
          do ll=lb,le; tens(l_out+ll)=slice(l_in+ll); enddo
          l_in=l_in+le+1_INTL; lb=0_INTL
          do i=2,dim_num
           l_out=l_out+bases_out(i); im(i)=im(i)+1
           if(im(i).ge.slice_ext(i)) then
            l_out=l_out-im(i)*bases_out(i); im(i)=0
           else
            exit
           endif
          enddo
         enddo sloop
!$OMP END PARALLEL
        else
         ierr=1 !zero-rank tensor
        endif
!        write(CONS_OUT,'("#DEBUG(dil_tensor_algebra::dil_tensor_insert): kernel time/error code: ",F10.4,1x,i3)')&
!        &thread_wtime(time_beg),ierr !debug
        return
        end subroutine dil_tensor_insert
!--------------------------------------------------------------------------------------------
#ifdef USE_MIC
!DIR$ ATTRIBUTES OFFLOAD:mic:: dil_tensor_transpose
#endif
        subroutine dil_tensor_transpose(dim_num,dim_extents,dim_transp,tens_in,tens_out,ierr) !PARALLEL (OMP)
!This subroutine reorders dimensions of a locally stored tensor.
        implicit none
!-------------------------------------------------
        integer(INTD), parameter:: real_kind=tensor_dp
        logical, parameter:: cache_efficiency=.true.
        integer(INTL), parameter:: cache_line_len=64/real_kind     !cache line length (words)
        integer(INTL), parameter:: cache_line_min=cache_line_len*2 !lower bound for the input/output minor volume: => L1_cache_line*2
        integer(INTL), parameter:: cache_line_lim=cache_line_len*4 !upper bound for the input/output minor volume: <= SQRT(L1_size)
        integer(INTL), parameter:: small_tens_size=2**10 !up to this size it is useless to apply cache efficiency (fully fits in L1)
        integer(INTL), parameter:: vec_size=2**8 !loop reorganization parameter for direct copy
#ifdef USE_MIC
!DIR$ ATTRIBUTES OFFLOAD:mic:: real_kind,cache_efficiency,cache_line_len,cache_line_min,cache_line_lim,small_tens_size,vec_size
!DIR$ ATTRIBUTES ALIGN:128:: real_kind,cache_efficiency,cache_line_len,cache_line_min,cache_line_lim,small_tens_size,vec_size
#endif
!------------------------------------------
        integer(INTD), intent(in):: dim_num                !in: tensor rank
        integer(INTD), intent(in):: dim_extents(1:dim_num) !in: tensor dimension extents
        integer(INTD), intent(in):: dim_transp(1:dim_num)  !in: dimension permutation (O2N)
        real(real_kind), intent(in):: tens_in(0:*)         !in: input tensor
        real(real_kind), intent(out):: tens_out(0:*)       !out: output (permuted) tensor
        integer(INTD), intent(inout):: ierr                !out: error code (0:success)
        integer(INTD):: i,j,k,l,m,n,k1,k2,ks,kf,split_in,split_out
        integer(INTD):: im(1:dim_num),n2o(0:dim_num+1),ipr(1:dim_num+1),dim_beg(1:dim_num),dim_end(1:dim_num)
        integer(INTL) bases_in(1:dim_num+1),bases_out(1:dim_num+1),bases_pri(1:dim_num+1),segs(0:MAX_THREADS)
        integer(INTL) bs,l0,l1,l2,l3,ll,lb,le,ls,l_in,l_out,seg_in,seg_out,vol_min,vol_ext
        logical trivial
        real(8) time_beg,tm
#ifdef USE_MIC
!DIR$ ATTRIBUTES ALIGN:128:: im,n2o,ipr,dim_beg,dim_end,bases_in,bases_out,bases_pri,segs
#endif
        ierr=0
        if(DIL_DEBUG) write(CONS_OUT,'(3x,"#DEBUG(DIL): Transposing ... ")',ADVANCE='NO')
        time_beg=thread_wtime() !debug
        if(dim_num.lt.0) then; ierr=1; return; elseif(dim_num.eq.0) then; tens_out(0)=tens_in(0); return; endif
!Check the index permutation:
        trivial=.true.; do i=1,dim_num; if(dim_transp(i).ne.i) then; trivial=.false.; exit; endif; enddo
        if(trivial.and.cache_efficiency) then
!Trivial index permutation (no permutation):
 !Compute indexing bases:
         bs=1_INTL; do i=1,dim_num; bases_in(i)=bs; bs=bs*dim_extents(i); enddo
 !Copy input to output:
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l0,l1)
!$OMP DO SCHEDULE(GUIDED)
         do l0=0_INTL,bs-1_INTL-mod(bs,vec_size),vec_size
          do l1=0_INTL,vec_size-1_INTL; tens_out(l0+l1)=tens_in(l0+l1); enddo
         enddo
!$OMP END DO NOWAIT
!$OMP SINGLE
         do l0=bs-mod(bs,vec_size),bs-1_INTL; tens_out(l0)=tens_in(l0); enddo
!$OMP END SINGLE
!$OMP END PARALLEL
        else
!Non-trivial index permutation:
 !Compute indexing bases:
         do i=1,dim_num; n2o(dim_transp(i))=i; enddo; n2o(dim_num+1)=dim_num+1 !get the N2O
         bs=1_INTL; do i=1,dim_num; bases_in(i)=bs; bs=bs*dim_extents(i); enddo; bases_in(dim_num+1)=bs
         bs=1_INTL; do i=1,dim_num; bases_out(n2o(i))=bs; bs=bs*dim_extents(n2o(i)); enddo; bases_out(dim_num+1)=bs
 !Configure cache-efficient algorithm:
         if(bs.le.small_tens_size.or.(.not.cache_efficiency)) then !tensor block is too small to think hard about it
          ipr(1:dim_num+1)=(/(j,j=1,dim_num+1)/); kf=dim_num !trivial priorities, all indices are minor
          split_in=kf; seg_in=dim_extents(split_in); split_out=kf; seg_out=dim_extents(split_out)
         else
          do k1=1,dim_num; if(bases_in(k1+1).ge.cache_line_lim) exit; enddo; k1=k1-1
          do k2=1,dim_num; if(bases_out(n2o(k2+1)).ge.cache_line_lim) exit; enddo; k2=k2-1
          do j=k1+1,dim_num; if(dim_transp(j).le.k2) then; k1=k1+1; else; exit; endif; enddo
          do j=k2+1,dim_num; if(n2o(j).le.k1) then; k2=k2+1; else; exit; endif; enddo
          if(bases_in(k1+1).lt.cache_line_min.and.bases_out(n2o(k2+1)).ge.cache_line_min) then !split the last minor input dim
           k1=k1+1; split_in=k1; seg_in=(cache_line_lim-1_INTL)/bases_in(split_in)+1_INTL
           split_out=n2o(k2); seg_out=dim_extents(split_out)
          elseif(bases_in(k1+1).ge.cache_line_min.and.bases_out(n2o(k2+1)).lt.cache_line_min) then !split the last minor output dim
           k2=k2+1; split_in=n2o(k2); seg_in=(cache_line_lim-1_INTL)/bases_out(split_in)+1_INTL
           split_out=k1; seg_out=dim_extents(split_out)
          elseif(bases_in(k1+1).lt.cache_line_min.and.bases_out(n2o(k2+1)).lt.cache_line_min) then !split both
           k1=k1+1; k2=k2+1
           if(k1.eq.n2o(k2)) then
            split_in=k1; seg_in=(cache_line_lim-1_INTL)/min(bases_in(split_in),bases_out(split_in))+1_INTL
            split_out=k1; seg_out=dim_extents(split_out)
           else
            split_in=k1; seg_in=(cache_line_lim-1_INTL)/bases_in(split_in)+1_INTL
            split_out=n2o(k2); seg_out=(cache_line_lim-1_INTL)/bases_out(split_out)+1_INTL
           endif
          else !split none
           split_in=k1; seg_in=dim_extents(split_in)
           split_out=n2o(k2); seg_out=dim_extents(split_out)
          endif
          vol_min=1_INTL
          if(seg_in.lt.dim_extents(split_in)) vol_min=vol_min*seg_in
          if(seg_out.lt.dim_extents(split_out)) vol_min=vol_min*seg_out
          if(vol_min.gt.1_INTL) then
           do j=1,k1
            if(j.ne.split_in.and.j.ne.split_out) vol_min=vol_min*dim_extents(j)
           enddo
           do j=1,k2
            l=n2o(j)
            if(l.gt.k1.and.l.ne.split_in.and.l.ne.split_out) vol_min=vol_min*dim_extents(l)
           enddo
           l=int((cache_line_lim*cache_line_lim)/vol_min,4)
           if(l.ge.2) then
            if(split_in.eq.split_out) then
             seg_in=seg_in*l
            else
             if(l.gt.4) then
              l=int(sqrt(float(l)),4)
              seg_in=min(seg_in*l,int(dim_extents(split_in),INTL))
              seg_out=min(seg_out*l,int(dim_extents(split_out),INTL))
             else
              seg_in=min(seg_in*l,int(dim_extents(split_in),INTL))
             endif
            endif
           endif
          endif
          l=0
          do while(l.lt.k1)
           l=l+1; ipr(l)=l; if(bases_in(l+1).ge.cache_line_min) exit
          enddo
          m=l+1
          j=0
          do while(j.lt.k2)
           j=j+1; n=n2o(j)
           if(n.ge.m) then; l=l+1; ipr(l)=n; endif
           if(bases_out(n2o(j+1)).ge.cache_line_min) exit
          enddo
          n=j+1
          do j=m,k1; if(dim_transp(j).ge.n) then; l=l+1; ipr(l)=j; endif; enddo
          do j=n,k2; if(n2o(j).gt.k1) then; l=l+1; ipr(l)=n2o(j); endif; enddo
          kf=l
          do j=k2+1,dim_num; if(n2o(j).gt.k1) then; l=l+1; ipr(l)=n2o(j); endif; enddo !kf is the length of the combined minor set
          ipr(dim_num+1)=dim_num+1 !special setting
         endif
         vol_ext=1_INTL; do j=kf+1,dim_num; vol_ext=vol_ext*dim_extents(ipr(j)); enddo !external volume
!         write(CONS_OUT,'("#DEBUG(dil_tensor_algebra::dil_tensor_transpose): extents:",99(1x,i5))') dim_extents(1:dim_num) !debug
!         write(CONS_OUT,'("#DEBUG(dil_tensor_algebra::dil_tensor_transpose): permutation:",99(1x,i2))') dim_transp(1:dim_num) !debug
!         write(CONS_OUT,'("#DEBUG(dil_tensor_algebra::dil_tensor_transpose): minor ",i3,": priority:",99(1x,i2))')&
!         &kf,ipr(1:dim_num) !debug
!         write(CONS_OUT,'("#DEBUG(dil_tensor_algebra::dil_tensor_transpose): vol_ext ",i11,": segs:",4(1x,i5))')&
!         &vol_ext,split_in,split_out,seg_in,seg_out !debug
 !Transpose:
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,m,n,ks,l0,l1,l2,l3,ll,lb,le,ls,l_in,l_out,vol_min,im,dim_beg,dim_end)
#ifdef VAR_OMP
         n=omp_get_thread_num(); m=omp_get_num_threads() !multi-threaded execution
#else
         n=0; m=1 !serial execution
#endif
!         if(n.eq.0) write(CONS_OUT,'("#DEBUG(dil_tensor_algebra::dil_tensor_transpose): number of threads = ",i5)') m !debug
         if(kf.lt.dim_num) then !external indices present
!$OMP MASTER
          segs(0)=0_INTL; call divide_segment_i8(vol_ext,int(m,INTL),segs(1:),i); do j=2,m; segs(j)=segs(j)+segs(j-1); enddo
          l0=1_INTL; do i=kf+1,dim_num; bases_pri(ipr(i))=l0; l0=l0*dim_extents(ipr(i)); enddo !priority bases
!$OMP END MASTER
!$OMP BARRIER
!$OMP FLUSH(segs,bases_pri)
          dim_beg(1:dim_num)=0; dim_end(1:dim_num)=dim_extents(1:dim_num)-1
          l2=dim_end(split_in); l3=dim_end(split_out); ls=bases_out(1)
          loop0: do l1=0_INTL,l3,seg_out !output dimension
           dim_beg(split_out)=l1; dim_end(split_out)=min(l1+seg_out-1_INTL,l3)
           do l0=0_INTL,l2,seg_in !input dimension
            dim_beg(split_in)=l0; dim_end(split_in)=min(l0+seg_in-1_INTL,l2)
            ll=segs(n); do i=dim_num,kf+1,-1; j=ipr(i); im(j)=ll/bases_pri(j); ll=ll-im(j)*bases_pri(j); enddo
            vol_min=1_INTL; do i=1,kf; j=ipr(i); vol_min=vol_min*(dim_end(j)-dim_beg(j)+1); im(j)=dim_beg(j); enddo
            l_in=0_INTL; do j=1,dim_num; l_in=l_in+im(j)*bases_in(j); enddo
            l_out=0_INTL; do j=1,dim_num; l_out=l_out+im(j)*bases_out(j); enddo
            le=dim_end(1)-dim_beg(1); lb=(segs(n+1)-segs(n))*vol_min; ks=0
            loop1: do while(lb.gt.0_INTL)
             do ll=0_INTL,le
              tens_out(l_out+ll*ls)=tens_in(l_in+ll)
             enddo
             lb=lb-(le+1_INTL)
             do i=2,dim_num
              j=ipr(i) !old index number
              if(im(j).lt.dim_end(j)) then
               im(j)=im(j)+1; l_in=l_in+bases_in(j); l_out=l_out+bases_out(j)
               ks=ks+1; exit
              else
               l_in=l_in-(im(j)-dim_beg(j))*bases_in(j); l_out=l_out-(im(j)-dim_beg(j))*bases_out(j); im(j)=dim_beg(j)
              endif
             enddo !i
             ks=ks-1; if(ks.lt.0) exit loop1
            enddo loop1
            if(lb.ne.0_INTL) then
             if(VERBOSE) write(CONS_OUT,'("ERROR(dil_tensor_algebra::dil_tensor_transpose): invalid remainder: ",i11,1x,i4)') lb,n
             ierr=2
             exit loop0
            endif
           enddo !l0
          enddo loop0 !l1
         else !external indices absent
!$OMP MASTER
          l0=1_INTL; do i=kf+1,dim_num; bases_pri(ipr(i))=l0; l0=l0*dim_extents(ipr(i)); enddo !priority bases
!$OMP END MASTER
!$OMP BARRIER
!$OMP FLUSH(bases_pri)
          dim_beg(1:dim_num)=0; dim_end(1:dim_num)=dim_extents(1:dim_num)-1
          l2=dim_end(split_in); l3=dim_end(split_out); ls=bases_out(1)
!$OMP DO SCHEDULE(DYNAMIC) COLLAPSE(2)
          do l1=0_INTL,l3,seg_out !output dimension
           do l0=0_INTL,l2,seg_in !input dimension
            dim_beg(split_out)=l1; dim_end(split_out)=min(l1+seg_out-1_INTL,l3)
            dim_beg(split_in)=l0; dim_end(split_in)=min(l0+seg_in-1_INTL,l2)
            vol_min=1_INTL; do i=1,kf; j=ipr(i); vol_min=vol_min*(dim_end(j)-dim_beg(j)+1); im(j)=dim_beg(j); enddo
            l_in=0_INTL; do j=1,dim_num; l_in=l_in+im(j)*bases_in(j); enddo
            l_out=0_INTL; do j=1,dim_num; l_out=l_out+im(j)*bases_out(j); enddo
            le=dim_end(1)-dim_beg(1); lb=vol_min; ks=0
            loop2: do while(lb.gt.0_INTL)
             do ll=0_INTL,le
              tens_out(l_out+ll*ls)=tens_in(l_in+ll)
             enddo
             lb=lb-(le+1_INTL)
             do i=2,dim_num
              j=ipr(i) !old index number
              if(im(j).lt.dim_end(j)) then
               im(j)=im(j)+1; l_in=l_in+bases_in(j); l_out=l_out+bases_out(j)
               ks=ks+1; exit
              else
               l_in=l_in-(im(j)-dim_beg(j))*bases_in(j); l_out=l_out-(im(j)-dim_beg(j))*bases_out(j); im(j)=dim_beg(j)
              endif
             enddo !i
             ks=ks-1; if(ks.lt.0) exit loop2
            enddo loop2
           enddo !l0
          enddo !l1
!$OMP END DO
         endif
!$OMP END PARALLEL
        endif !trivial or not
        tm=thread_wtime(time_beg) !debug
        if(DIL_DEBUG) write(CONS_OUT,'("Done: ",F10.4," s, ",F10.4," GB/s: Status ",i9)')&
        &tm,dble(2_INTL*bs*real_kind)/(tm*1048576d0*1024d0),ierr !debug
!        write(CONS_OUT,'("#DEBUG(dil_tensor_algebra::dil_tensor_transpose): Done: ",F10.4," sec, ",F10.4," GB/s, error ",i3)')&
!        &tm,dble(2_INTL*bs*real_kind)/(tm*1024d0*1024d0*1024d0),ierr !debug
        return
        end subroutine dil_tensor_transpose
!--------------------------------------------------------------------------------------
        subroutine dil_get_next_tile_signa(tens_full,subtens,signa,sdims,tile_num,ierr) !SERIAL
!This subroutine iterates over the tiles forming a given tensor part.
!At the first call, one must set ierr > 0. On success, ierr will be cleared to zero and
!each subsequent iteration will return ierr=0 (on success). At the end, ierr=-1 will be
!returned, meaning that the iterations are over. Positive ierr on return means an error.
        implicit none
        type(tensor), intent(in):: tens_full        !in: full tensor specification
        type(subtens_t), intent(in):: subtens       !in: tensor part specification
        integer(INTD), intent(inout):: signa(1:*)   !inout: current tile signature (dimension bases)
        integer(INTD), intent(out):: sdims(1:*)     !out: current tile dimension extents
        integer, intent(out):: tile_num             !out: tile number (numeration starts from 1)
        integer(INTD), intent(inout):: ierr         !inout: begin flag (in), error code (out)
!---------------------------------------------
        logical, parameter:: NO_CHECK=.false. !argument check
!---------------------------------------------
        integer(INTD):: i,m,n

!Get tile signature:
        n=tens_full%mode
        if(n.ge.0) then
         if(ierr.eq.0) then !subsequent iterations
          i=1
          do while(i.le.n)
           signa(i)=signa(i)+tens_full%tdim(i)
           if(signa(i).lt.subtens%lbnd(i)+subtens%dims(i)) then
            exit
           else
            signa(i)=subtens%lbnd(i); i=i+1
           endif
          enddo
          if(i.gt.n) then; ierr=-1; return; endif !iterations are over
         elseif(ierr.gt.0) then !first call (positive <ierr>)
          ierr=0
          if(.not.NO_CHECK) then
           if(subtens%rank.eq.n.and.n.le.MAX_TENSOR_RANK) then
            do i=1,n
             if(subtens%lbnd(i).lt.1.or.subtens%lbnd(i).gt.tens_full%dims(i)) then; ierr=1; return; endif
             if(subtens%dims(i).lt.1) then; ierr=2; return; endif
             if(subtens%lbnd(i)+subtens%dims(i)-1_INTD.gt.tens_full%dims(i)) then; ierr=3; return; endif
             if(tens_full%ntpm(i).le.0) then; ierr=4; return; endif
            enddo
           else
            ierr=5; return
           endif
          endif
          signa(1:n)=subtens%lbnd(1:n) !signature of the 1st tile (lower bounds)
         else !iterations already over (negative <ierr>)
          return
         endif
!Get tile dimensions:
         do i=1,n; sdims(i)=min(tens_full%dims(i)-signa(i)+1,tens_full%tdim(i)); enddo
!Get the global tile number (`Keep consistent with get_cidx):
         tile_num=1; m=1 !Assumes that Global Tile Numeration starts from 1
         do i=1,n
          tile_num=tile_num+((signa(i)-1_INTD)/tens_full%tdim(i))*m
          m=m*tens_full%ntpm(i)
         enddo
        else
         ierr=6
        endif
        return
        end subroutine dil_get_next_tile_signa
!-------------------------------------------------------------------------------
        subroutine dil_tensor_prefetch_start(tens_arr,tens_part,buf,ierr,locked) !SERIAL (MPI)
!This subroutine starts collection of all tiles necessary for constructing a given tensor part.
        implicit none
        type(tensor), intent(in):: tens_arr     !in: tensor stored distributively in terms of tiles
        type(subtens_t), intent(in):: tens_part !in: tensor part specification
        type(arg_buf_t), intent(inout):: buf    !out: local buffer where the tiles will be put in
        integer(INTD), intent(inout):: ierr     !out: error code (0:success)
        logical, intent(in), optional:: locked  !in: if .TRUE., MPI windows are assumed already locked
        integer(INTD):: i,k,tile_host,signa(1:MAX_TENSOR_RANK),tile_dims(1:MAX_TENSOR_RANK)
        integer(INTL):: tile_vol,buf_end
        integer:: tile_num,tile_win,dpos,didx
        type(rank_win_cont_t):: rwc
        logical:: new_rw,win_lck

        ierr=0; call dil_rank_window_clean(rwc)
        if(present(locked)) then; win_lck=locked; else; win_lck=.false.; endif
        if(DIL_DEBUG) write(CONS_OUT,'(2x,"#DEBUG(DIL): Prefetching (",16(1x,i6,":",i6,","),")")')&
        &(/(tens_part%lbnd(i),tens_part%lbnd(i)+tens_part%dims(i)-1_INTD,i=1,tens_arr%mode)/) !debug
        buf_end=0_INTL; k=DIL_FIRST_CALL
        do while(k.ge.0) !k<0: iterations are over
         call dil_get_next_tile_signa(tens_arr,tens_part,signa,tile_dims,tile_num,k)
         if(k.eq.0) then
          call get_residence_of_tile(tens_arr,tile_num,tile_host,dpos,didx,tile_win)
          tile_vol=1_INTL; do i=1,tens_arr%mode; tile_vol=tile_vol*tile_dims(i); enddo
          new_rw=dil_rank_window_new(rwc,tile_host,tile_win,i); if(i.ne.0) ierr=ierr+1
          if(DIL_DEBUG) write(CONS_OUT,'(3x,"#DEBUG(DIL): Lock+Get on ",i9,"(",l1,"): ",i7,"/",i11)',ADVANCE='NO')&
          &tile_num,new_rw,tile_host,tens_arr%wi(tile_win)
          if((.not.win_lck).and.new_rw) call tensor_mpi_win_lock(int(tile_host,tensor_mpi_kind),tens_arr%wi(tile_win),'s')
          call tensor_get_tile(tens_arr,int(tile_num,kind=tensor_standard_int),&
             &buf%buf_ptr(buf_end+1_INTL:),tile_vol,lock_set=.true.)
          if(DIL_DEBUG) write(CONS_OUT,'(" [Ok]:",16(1x,i6))') signa(1:tens_arr%mode)
          buf_end=buf_end+tile_vol
         elseif(k.gt.0) then
          ierr=-1; return
         endif
        enddo
        return
        end subroutine dil_tensor_prefetch_start
!----------------------------------------------------------------------------------
        subroutine dil_tensor_prefetch_complete(tens_arr,tens_part,buf,ierr,locked) !SERIAL (MPI)
!This subroutine completes collection of all tiles necessary for constructing a given tensor part.
        implicit none
        type(tensor), intent(in):: tens_arr     !in: tensor stored distributively in terms of tiles
        type(subtens_t), intent(in):: tens_part !in: tensor part specification
        type(arg_buf_t), intent(inout):: buf    !out: local buffer where the tiles will be put in
        integer(INTD), intent(inout):: ierr     !out: error code (0:success)
        logical, intent(in), optional:: locked  !in: if .TRUE., MPI windows are assumed already locked
        integer(INTD):: i,k,tile_host,signa(1:MAX_TENSOR_RANK),tile_dims(1:MAX_TENSOR_RANK)
        integer:: tile_num,tile_win,dpos,didx
        type(rank_win_cont_t):: rwc
        logical:: new_rw,win_lck

        ierr=0; call dil_rank_window_clean(rwc)
        if(present(locked)) then; win_lck=locked; else; win_lck=.false.; endif
        if(DIL_DEBUG) write(CONS_OUT,'(2x,"#DEBUG(DIL): Prefetching (",16(1x,i6,":",i6,","),")")')&
        &(/(tens_part%lbnd(i),tens_part%lbnd(i)+tens_part%dims(i)-1_INTD,i=1,tens_arr%mode)/) !debug
        k=DIL_FIRST_CALL
        do while(k.ge.0) !k<0: iterations are over
         call dil_get_next_tile_signa(tens_arr,tens_part,signa,tile_dims,tile_num,k)
         if(k.eq.0) then
          call get_residence_of_tile(tens_arr,tile_num,tile_host,dpos,didx,tile_win)
          new_rw=dil_rank_window_new(rwc,tile_host,tile_win,i); if(i.ne.0) ierr=ierr+1
          if(DIL_DEBUG) write(CONS_OUT,'(3x,"#DEBUG(DIL): Unlock(Get) on ",i9,"(",l1,"): ",i7,"/",i11)',ADVANCE='NO')&
          &tile_num,new_rw,tile_host,tens_arr%wi(tile_win)
          if(new_rw) then
           if(win_lck) then
            call tensor_mpi_win_flush(tens_arr%wi(tile_win),int(tile_host,tensor_mpi_kind))
           else
            call tensor_mpi_win_unlock(int(tile_host,tensor_mpi_kind),tens_arr%wi(tile_win))
           endif
          endif
          if(DIL_DEBUG) write(CONS_OUT,'(" [Ok]",16(1x,i6))') signa(1:tens_arr%mode)
         elseif(k.gt.0) then
          ierr=-1; return
         endif
        enddo
        return
        end subroutine dil_tensor_prefetch_complete
!-----------------------------------------------------------------------------
        subroutine dil_tensor_upload_start(tens_arr,tens_part,buf,ierr,locked) !SERIAL (MPI)
!This subroutine starts upload of all tiles constituting a given tensor part.
        implicit none
        type(tensor), intent(in):: tens_arr     !in: tensor stored distributively in terms of tiles
        type(subtens_t), intent(in):: tens_part !in: tensor part specification
        type(arg_buf_t), intent(inout):: buf    !in: local buffer containing the tiles to upload
        integer(INTD), intent(inout):: ierr     !out: error code (0:success)
        logical, intent(in), optional:: locked  !in: if .TRUE., MPI windows are assumed already locked
        integer(INTD):: i,k,tile_host,signa(1:MAX_TENSOR_RANK),tile_dims(1:MAX_TENSOR_RANK)
        integer(INTL):: tile_vol,buf_end
        integer:: tile_num,tile_win,dpos,didx
        type(rank_win_cont_t):: rwc
        logical:: new_rw,win_lck

        ierr=0; call dil_rank_window_clean(rwc)
        if(present(locked)) then; win_lck=locked; else; win_lck=.false.; endif
        if(DIL_DEBUG) write(CONS_OUT,'(2x,"#DEBUG(DIL): Uploading (",16(1x,i6,":",i6,","),")")')&
        &(/(tens_part%lbnd(i),tens_part%lbnd(i)+tens_part%dims(i)-1_INTD,i=1,tens_arr%mode)/) !debug
        buf_end=0_INTL; k=DIL_FIRST_CALL
        do while(k.ge.0)
         call dil_get_next_tile_signa(tens_arr,tens_part,signa,tile_dims,tile_num,k)
         if(k.eq.0) then
          call get_residence_of_tile(tens_arr,tile_num,tile_host,dpos,didx,tile_win)
          tile_vol=1_INTL; do i=1,tens_arr%mode; tile_vol=tile_vol*tile_dims(i); enddo
          new_rw=dil_rank_window_new(rwc,tile_host,tile_win,i); if(i.ne.0) ierr=ierr+1
          if(DIL_DEBUG) write(CONS_OUT,'(3x,"#DEBUG(DIL): Lock+Accumulate on ",i9,"(",l1,"): ",i7,"/",i11)',ADVANCE='NO')&
          &tile_num,new_rw,tile_host,tens_arr%wi(tile_win)
          if((.not.win_lck).and.new_rw) call tensor_mpi_win_lock(int(tile_host,tensor_mpi_kind),tens_arr%wi(tile_win),'s')
          call tensor_accumulate_tile(tens_arr,tile_num,buf%buf_ptr(buf_end+1_INTL:),tile_vol,lock_set=.true.)
          if(DIL_DEBUG) write(CONS_OUT,'(" [Ok]",16(1x,i6))') signa(1:tens_arr%mode)
          buf_end=buf_end+tile_vol
         elseif(k.gt.0) then
          ierr=-1; return
         endif
        enddo
        return
        end subroutine dil_tensor_upload_start
!-----------------------------------------------------------------------------------------------------
        subroutine dil_tensor_upload_complete(tens_arr,tens_part,buf,ierr,locked,num_async,list_async) !SERIAL (MPI)
!This subroutine completes upload of all tiles constituting a given tensor part.
!If both <num_async> and <list_async> are present, no MPI passive synchronization will be done.
!Instead, the list of outstanding MPI uploads will be returned for later finalization.
        implicit none
        type(tensor), intent(inout):: tens_arr                     !inout: tensor stored distributively in terms of tiles
        type(subtens_t), intent(in):: tens_part                    !in: tensor part specification
        type(arg_buf_t), intent(in):: buf                          !in: local buffer containing the tiles to upload
        integer(INTD), intent(inout):: ierr                        !out: error code (0:success)
        logical, intent(in), optional:: locked                     !in: if .TRUE., MPI windows are assumed already locked
        integer(INTD), intent(inout), optional:: num_async         !inout: number of the outstanding MPI uploads left
        type(rank_win_t), intent(inout), optional:: list_async(1:) !out: list of the outstanding MPI uploads
        integer(INTD):: i,k,tile_host,signa(1:MAX_TENSOR_RANK),tile_dims(1:MAX_TENSOR_RANK),max_async
        integer:: tile_num,tile_win,dpos,didx
        type(rank_win_cont_t):: rwc
        logical:: new_rw,win_lck,async

        ierr=0; call dil_rank_window_clean(rwc)
        if(present(num_async).and.present(list_async)) then !both arguments must be present for asynchronous scenario
         async=.true.; max_async=ubound(list_async,1)
        else
         async=.false.
        endif
        if(present(locked)) then; win_lck=locked; else; win_lck=.false.; endif
        if(DIL_DEBUG) write(CONS_OUT,'(2x,"#DEBUG(DIL): Uploading (",16(1x,i6,":",i6,","),")")')&
        &(/(tens_part%lbnd(i),tens_part%lbnd(i)+tens_part%dims(i)-1_INTD,i=1,tens_arr%mode)/) !debug
        k=DIL_FIRST_CALL
        do while(k.ge.0)
         call dil_get_next_tile_signa(tens_arr,tens_part,signa,tile_dims,tile_num,k)
         if(k.eq.0) then
          call get_residence_of_tile(tens_arr,tile_num,tile_host,dpos,didx,tile_win)
          new_rw=dil_rank_window_new(rwc,tile_host,tile_win,i); if(i.ne.0) ierr=ierr+1
          if(DIL_DEBUG) write(CONS_OUT,'(3x,"#DEBUG(DIL): Unlock(Accumulate) on ",i9,"(",l1,"): ",i7,"/",i11)',ADVANCE='NO')&
          &tile_num,new_rw,tile_host,tens_arr%wi(tile_win)
          if(new_rw) then
           if(.not.async) then
            if(win_lck) then
             call tensor_mpi_win_flush(tens_arr%wi(tile_win),int(tile_host,tensor_mpi_kind))
            else
             call tensor_mpi_win_unlock(int(tile_host,tensor_mpi_kind),tens_arr%wi(tile_win))
            endif
           else
            num_async=num_async+1
            if(num_async.le.max_async) then
             list_async(num_async)%rank=tile_host
             list_async(num_async)%window=tens_arr%wi(tile_win)
            else
             ierr=1000
            endif
           endif
          endif
          if(DIL_DEBUG) write(CONS_OUT,'(" [Ok]",16(1x,i6))') signa(1:tens_arr%mode)
         elseif(k.gt.0) then
          ierr=-1; return
         endif
        enddo
        if(DIL_DEBUG.and.async) write(CONS_OUT,'(2x,"#DEBUG(DIL): Recorded for later finalization.")')
        return
        end subroutine dil_tensor_upload_complete
!-------------------------------------------------------------------------------
        subroutine dil_tens_unpack_from_tiles(tens_arr,tens_part,bufi,bufo,ierr) !PARALLEL (OMP)
!This subroutine unpacks tiles into a dense tensor slice.
        implicit none
        type(tensor), intent(in):: tens_arr     !in: tensor stored distributively in terms of tiles
        type(subtens_t), intent(in):: tens_part !in: tensor part specification
        type(arg_buf_t), intent(in):: bufi      !in: local buffer containing the tiles
        type(arg_buf_t), intent(inout):: bufo   !out: local buffer that will contain the tensor slice
        integer(INTD), intent(inout):: ierr     !out: error code (0:success)
        integer(INTD):: i,k,n,signa(1:MAX_TENSOR_RANK),tile_dims(1:MAX_TENSOR_RANK)
        integer(INTL):: tile_vol,buf_end
        integer:: tile_num
        real(8):: time_beg,tm

        ierr=0
        buf_end=0_INTL; n=tens_arr%mode; k=DIL_FIRST_CALL
        do while(k.ge.0)
         call dil_get_next_tile_signa(tens_arr,tens_part,signa,tile_dims,tile_num,k)
         if(k.eq.0) then
          if(DIL_DEBUG) time_beg=thread_wtime()
          tile_vol=1_INTL; do i=1,n; tile_vol=tile_vol*tile_dims(i); enddo
          if(DIL_DEBUG) write(CONS_OUT,'(3x,"#DEBUG(DIL): Unpacking:",16(1x,i6))',ADVANCE='NO') signa(1:n)
          call dil_tensor_insert(n,bufo%buf_ptr,tens_part%dims,bufi%buf_ptr(buf_end+1_INTL:),tile_dims,&
                                &signa(1:n)-tens_part%lbnd(1:n),i)
          buf_end=buf_end+tile_vol
          if(DIL_DEBUG) tm=thread_wtime(time_beg)
          if(DIL_DEBUG) write(CONS_OUT,'(": Unpacked: ",F10.4," s: ",F10.4," GB/s: Status ",i9)')&
           &tm,dble(2_INTL*tile_vol*tensor_dp)/(tm*1024d0*1024d0*1024d0),i
          if(i.ne.0) then; ierr=1; return; endif
         elseif(k.gt.0) then
          ierr=2; return
         endif
        enddo
        return
        end subroutine dil_tens_unpack_from_tiles
!-----------------------------------------------------------------------------
        subroutine dil_tens_pack_into_tiles(tens_arr,tens_part,bufi,bufo,ierr) !PARALLEL (OMP)
!This subroutine packs a dense tensor slice into tiles.
        implicit none
        type(tensor), intent(in):: tens_arr     !in: tensor stored distributively in terms of tiles
        type(subtens_t), intent(in):: tens_part !in: tensor part specification
        type(arg_buf_t), intent(in):: bufi      !in: local buffer containing the tensor slice
        type(arg_buf_t), intent(inout):: bufo   !out: local buffer that will contain the tiles
        integer(INTD), intent(inout):: ierr     !out: error code (0:success)
        integer(INTD):: i,k,n,signa(1:MAX_TENSOR_RANK),tile_dims(1:MAX_TENSOR_RANK)
        integer(INTL):: tile_vol,buf_end
        integer:: tile_num
        real(8):: time_beg,tm

        ierr=0
        buf_end=0_INTL; n=tens_arr%mode; k=DIL_FIRST_CALL
        do while(k.ge.0)
         call dil_get_next_tile_signa(tens_arr,tens_part,signa,tile_dims,tile_num,k)
         if(k.eq.0) then
          if(DIL_DEBUG) time_beg=thread_wtime()
          tile_vol=1_INTL; do i=1,n; tile_vol=tile_vol*tile_dims(i); enddo
          if(DIL_DEBUG) write(CONS_OUT,'(3x,"#DEBUG(DIL): Packing:",16(1x,i6))',ADVANCE='NO') signa(1:n)
          call dil_tensor_slice(n,bufi%buf_ptr,tens_part%dims,bufo%buf_ptr(buf_end+1_INTL:),tile_dims,&
                               &signa(1:n)-tens_part%lbnd(1:n),i)
          buf_end=buf_end+tile_vol
          if(DIL_DEBUG) tm=thread_wtime(time_beg)
          if(DIL_DEBUG) write(CONS_OUT,'(": Packed: ",F10.4," s: ",F10.4," GB/s: Status ",i9)')&
           &tm,dble(2_INTL*tile_vol*tensor_dp)/(tm*1024d0*1024d0*1024d0),i
          if(i.ne.0) then; ierr=1; return; endif
         elseif(k.gt.0) then
          ierr=2; return
         endif
        enddo
        return
        end subroutine dil_tens_pack_into_tiles
!---------------------------------------------------------------------------
        subroutine dil_tens_fetch_start(tens_arr,tens_part,bufi,ierr,locked) !SERIAL (MPI)
!This subroutine starts fetching all tiles necessary for constructing a given tensor part.
        implicit none
        type(tensor), intent(in):: tens_arr     !in: tensor stored distributively in terms of tiles
        type(subtens_t), intent(in):: tens_part !in: tensor part specification
        real(tensor_dp), intent(inout):: bufi(1:*)  !out: local buffer where the tiles will be put in
        integer(INTD), intent(inout):: ierr     !out: error code (0:success)
        logical, intent(in), optional:: locked  !in: if .TRUE., MPI windows are assumed already locked
        integer(INTD):: i,k,tile_host,signa(1:MAX_TENSOR_RANK),tile_dims(1:MAX_TENSOR_RANK)
        integer(INTL):: tile_vol,buf_end
        integer:: tile_num,tile_win,dpos,didx
        type(rank_win_cont_t):: rwc
        logical:: new_rw,win_lck

        ierr=0; call dil_rank_window_clean(rwc)
        if(present(locked)) then; win_lck=locked; else; win_lck=.false.; endif
!        if(DIL_DEBUG) write(CONS_OUT,'(2x,"#DEBUG(DIL): Fetching (",16(1x,i6,":",i6,","),")")')&
!         &(/(tens_part%lbnd(i),tens_part%lbnd(i)+tens_part%dims(i)-1_INTD,i=1,tens_arr%mode)/) !debug
        buf_end=0_INTL; k=DIL_FIRST_CALL
        do while(k.ge.0) !k<0: iterations are over
         call dil_get_next_tile_signa(tens_arr,tens_part,signa,tile_dims,tile_num,k)
         if(k.eq.0) then
          call get_residence_of_tile(tens_arr,tile_num,tile_host,dpos,didx,tile_win)
          tile_vol=1_INTL; do i=1,tens_arr%mode; tile_vol=tile_vol*tile_dims(i); enddo
          new_rw=dil_rank_window_new(rwc,tile_host,tile_win,i); if(i.ne.0) ierr=ierr+1
!          if(DIL_DEBUG) write(CONS_OUT,'(3x,"#DEBUG(DIL): Lock+Get on ",i9,"(",l1,"): ",i7,"/",i11)',ADVANCE='NO')&
!           &tile_num,new_rw,tile_host,tens_arr%wi(tile_win)
          if((.not.win_lck).and.new_rw) call tensor_mpi_win_lock(int(tile_host,tensor_mpi_kind),tens_arr%wi(tile_win),'s')
          call tensor_get_tile(tens_arr,int(tile_num,kind=tensor_standard_int),&
             &bufi(buf_end+1_INTL:buf_end+tile_vol),tile_vol,lock_set=.true.)
!          if(DIL_DEBUG) write(CONS_OUT,'(" [Ok]:",16(1x,i6))') signa(1:tens_arr%mode)
          buf_end=buf_end+tile_vol
         elseif(k.gt.0) then
          ierr=-1; return
         endif
        enddo
        return
        end subroutine dil_tens_fetch_start
!--------------------------------------------------------------------------------------
        subroutine dil_tens_fetch_finish_prep(tens_arr,tens_part,bufi,bufo,ierr,locked) !PARALLEL (OMP)
!This subroutine finishes tile fetching and unpacks the tiles into a dense tensor slice.
        implicit none
        type(tensor), intent(in):: tens_arr     !in: tensor stored distributively in terms of tiles
        type(subtens_t), intent(in):: tens_part !in: tensor part specification
        real(tensor_dp), intent(inout):: bufi(1:*)  !in: local buffer containing the tiles
        real(tensor_dp), intent(inout):: bufo(1:*)  !out: local buffer that will contain the dense tensor slice
        integer(INTD), intent(inout):: ierr     !out: error code (0:success)
        logical, intent(in), optional:: locked  !in: if .TRUE., MPI windows are assumed already locked
        integer(INTD):: i,k,n,tile_host,signa(1:MAX_TENSOR_RANK),tile_dims(1:MAX_TENSOR_RANK)
        integer(INTL):: tile_vol,buf_end
        integer:: tile_num,tile_win,dpos,didx
        type(rank_win_cont_t):: rwc
        logical:: new_rw,win_lck

        ierr=0; call dil_rank_window_clean(rwc)
        if(present(locked)) then; win_lck=locked; else; win_lck=.false.; endif
        buf_end=0_INTL; n=tens_arr%mode; k=DIL_FIRST_CALL
        do while(k.ge.0) !k<0: iterations are over
         call dil_get_next_tile_signa(tens_arr,tens_part,signa,tile_dims,tile_num,k)
         if(k.eq.0) then
          call get_residence_of_tile(tens_arr,tile_num,tile_host,dpos,didx,tile_win)
          tile_vol=1_INTL; do i=1,n; tile_vol=tile_vol*tile_dims(i); enddo
          new_rw=dil_rank_window_new(rwc,tile_host,tile_win,i); if(i.ne.0) ierr=ierr+1
!          if(DIL_DEBUG) write(CONS_OUT,'(3x,"#DEBUG(DIL): Unlock(Get) on ",i9,"(",l1,"): ",i7,"/",i11)',ADVANCE='NO')&
!           &tile_num,new_rw,tile_host,tens_arr%wi(tile_win)
          if(new_rw) then
           if(win_lck) then
            call tensor_mpi_win_flush(tens_arr%wi(tile_win),int(tile_host,tensor_mpi_kind))
           else
            call tensor_mpi_win_unlock(int(tile_host,tensor_mpi_kind),tens_arr%wi(tile_win))
           endif
          endif
!          if(DIL_DEBUG) write(CONS_OUT,'(" [Ok]",16(1x,i6))') signa(1:n)
!          if(DIL_DEBUG) write(CONS_OUT,'(3x,"#DEBUG(DIL): Unpacking:",16(1x,i6))',ADVANCE='NO') signa(1:n)
          call dil_tensor_insert(n,bufo,tens_part%dims,bufi(buf_end+1_INTL:buf_end+tile_vol),tile_dims,&
                                &signa(1:n)-tens_part%lbnd(1:n),i)
          buf_end=buf_end+tile_vol
!          if(DIL_DEBUG) write(CONS_OUT,'(": Unpacked: Status ",i9)') i
          if(i.ne.0) then; ierr=-2; return; endif
         elseif(k.gt.0) then
          ierr=-1; return
         endif
        enddo
        return
        end subroutine dil_tens_fetch_finish_prep
!--------------------------------------------------------------------------------------
        subroutine dil_tens_prep_upload_start(tens_arr,tens_part,bufi,bufo,ierr,locked) !SERIAL (MPI)
!This subroutine packs a dense tensor slice into tiles and starts uploading them.
        implicit none
        type(tensor), intent(inout):: tens_arr  !inout: tensor stored distributively in terms of tiles
        type(subtens_t), intent(in):: tens_part !in: tensor part specification
        real(tensor_dp), intent(in):: bufi(1:*)     !in: local buffer containing the dense tensor slice
        real(tensor_dp), intent(inout):: bufo(1:*)  !tmp: temporary buffer from where the tiles will be uploaded
        integer(INTD), intent(inout):: ierr     !out: error code (0:success)
        logical, intent(in), optional:: locked  !in: if .TRUE., MPI windows are assumed already locked
        integer(INTD):: i,k,n,tile_host,signa(1:MAX_TENSOR_RANK),tile_dims(1:MAX_TENSOR_RANK)
        integer(INTL):: tile_vol,buf_end
        integer:: tile_num,tile_win,dpos,didx
        type(rank_win_cont_t):: rwc
        logical:: new_rw,win_lck
        real(8):: time_beg,tm

        ierr=0; call dil_rank_window_clean(rwc); n=tens_arr%mode
        if(present(locked)) then; win_lck=locked; else; win_lck=.false.; endif
!        if(DIL_DEBUG) write(CONS_OUT,'(2x,"#DEBUG(DIL): Uploading (",16(1x,i6,":",i6,","),")")')&
!         &(/(tens_part%lbnd(i),tens_part%lbnd(i)+tens_part%dims(i)-1_INTD,i=1,n)/) !debug
        buf_end=0_INTL; k=DIL_FIRST_CALL
        do while(k.ge.0) !k<0: iterations are over
         call dil_get_next_tile_signa(tens_arr,tens_part,signa,tile_dims,tile_num,k)
         if(k.eq.0) then
!          if(DIL_DEBUG) time_beg=thread_wtime()
          tile_vol=1_INTL; do i=1,n; tile_vol=tile_vol*tile_dims(i); enddo
!          if(DIL_DEBUG) write(CONS_OUT,'(3x,"#DEBUG(DIL): Packing:",16(1x,i6))',ADVANCE='NO') signa(1:n)
          call dil_tensor_slice(n,bufi,tens_part%dims,bufo(buf_end+1_INTL:buf_end+tile_vol),tile_dims,&
                               &signa(1:n)-tens_part%lbnd(1:n),i)
!          if(DIL_DEBUG) tm=thread_wtime(time_beg)
!          if(DIL_DEBUG) write(CONS_OUT,'(": Packed: ",F10.4," s: ",F10.4," GB/s: Status ",i9)')&
!           &tm,dble(2_INTL*tile_vol*tensor_dp)/(tm*1024d0*1024d0*1024d0),i
          if(i.ne.0) then; ierr=1; return; endif
          call get_residence_of_tile(tens_arr,tile_num,tile_host,dpos,didx,tile_win)
          new_rw=dil_rank_window_new(rwc,tile_host,tile_win,i); if(i.ne.0) ierr=ierr+1
!          if(DIL_DEBUG) write(CONS_OUT,'(3x,"#DEBUG(DIL): Lock+Accumulate on ",i9,"(",l1,"): ",i7,"/",i11)',ADVANCE='NO')&
!           &tile_num,new_rw,tile_host,tens_arr%wi(tile_win)
          if((.not.win_lck).and.new_rw) call tensor_mpi_win_lock(int(tile_host,tensor_mpi_kind),tens_arr%wi(tile_win),'s')
          call tensor_accumulate_tile(tens_arr,tile_num,bufo(buf_end+1_INTL:buf_end+tile_vol),tile_vol,lock_set=.true.)
!          if(DIL_DEBUG) write(CONS_OUT,'(" [Ok]:",16(1x,i6))') signa(1:tens_arr%mode)
          buf_end=buf_end+tile_vol
         elseif(k.gt.0) then
          ierr=-1; return
         endif
        enddo
        return
        end subroutine dil_tens_prep_upload_start
!-----------------------------------------------------------------------------
        subroutine dil_tens_upload_finish(tens_arr,tens_part,bufo,ierr,locked) !SERIAL (MPI)
!This subroutine completes upload of all tiles constituting a given dense tensor slice.
        implicit none
        type(tensor), intent(inout):: tens_arr  !inout: tensor stored distributively in terms of tiles
        type(subtens_t), intent(in):: tens_part !in: tensor part specification
        real(tensor_dp), intent(inout):: bufo(1:*)  !in: local buffer containing the tiles to upload
        integer(INTD), intent(inout):: ierr     !out: error code (0:success)
        logical, intent(in), optional:: locked  !in: if .TRUE., MPI windows are assumed already locked
        integer(INTD):: i,k,tile_host,signa(1:MAX_TENSOR_RANK),tile_dims(1:MAX_TENSOR_RANK)
        integer:: tile_num,tile_win,dpos,didx
        type(rank_win_cont_t):: rwc
        logical:: new_rw,win_lck

        ierr=0; call dil_rank_window_clean(rwc)
        if(present(locked)) then; win_lck=locked; else; win_lck=.false.; endif
        k=DIL_FIRST_CALL
        do while(k.ge.0)
         call dil_get_next_tile_signa(tens_arr,tens_part,signa,tile_dims,tile_num,k)
         if(k.eq.0) then
          call get_residence_of_tile(tens_arr,tile_num,tile_host,dpos,didx,tile_win)
          new_rw=dil_rank_window_new(rwc,tile_host,tile_win,i); if(i.ne.0) ierr=ierr+1
!          if(DIL_DEBUG) write(CONS_OUT,'(3x,"#DEBUG(DIL): Unlock(Accumulate) on ",i9,"(",l1,"): ",i7,"/",i11)',ADVANCE='NO')&
!          &tile_num,new_rw,tile_host,tens_arr%wi(tile_win)
          if(new_rw) then
           if(win_lck) then
            call tensor_mpi_win_flush(tens_arr%wi(tile_win),int(tile_host,tensor_mpi_kind))
           else
            call tensor_mpi_win_unlock(int(tile_host,tensor_mpi_kind),tens_arr%wi(tile_win))
           endif
          endif
!          if(DIL_DEBUG) write(CONS_OUT,'(" [Ok]",16(1x,i6))') signa(1:tens_arr%mode)
         elseif(k.gt.0) then
          ierr=-1; return
         endif
        enddo
        return
        end subroutine dil_tens_upload_finish
!------------------------------------------------------------------------
        subroutine dil_divide_space_int4(ndim,dims,subvol,segs,ierr,algn) !SERIAL
!This subroutine divides an ndim-dimensional block with extents dims(1:ndim)
!into smaller ndim-dimensional blocks with extents segs(1:ndim) such that
!the volume of each such a block approximately equals to <subvol> while
!the ratio Volume/Surface is maximized.
        implicit none
        integer(INTD), intent(in):: ndim                !in: number of dimensions
        integer(INTD), intent(in):: dims(1:*)           !in: dimension extents
        integer(INTL), intent(in):: subvol              !in: target subblock volume
        integer(INTD), intent(inout):: segs(1:*)        !out: dimension segmentation
        integer(INTD), intent(inout):: ierr             !out: error code (0:success)
        integer(INTD), intent(in), optional:: algn(1:*) !in: dimension tiling (default is 1: no tiling)
!---------------------------------------------------------
        logical, parameter:: ARG_CHECK=.true. !argument check
!--------------------------------------------
        integer(INTD):: i,j,m,n,alg(1:ndim),trn(0:ndim)
        real(8):: vol
        logical:: algf

        ierr=0
        if(subvol.gt.0) then
         if(ndim.gt.0) then
          if(ARG_CHECK) then
           do i=1,ndim; if(dims(i).le.0) then; ierr=1; return; endif; enddo
          endif
          if(present(algn)) then
           alg(1:ndim)=algn(1:ndim); algf=.true.
          else
           alg(1:ndim)=1_INTD; algf=.false.
          endif
          trn(0:ndim)=(/+1_INTD,(i,i=1_INTD,ndim)/)
          if(algf) then
           call merge_sort_key_int(ndim,alg,trn)
          else
           call merge_sort_key_int(ndim,dims,trn)
          endif
          n=ndim; vol=dble(subvol)
          do i=ndim,1,-1
           j=trn(i)
           m=nint(vol**(1d0/dble(n))) !uniform segment length (ideal case: cubic)
           if(m.le.dims(j)) then
            segs(j)=max(m-mod(m,alg(j)),alg(j))
            if(mod(segs(j),alg(j)).ne.0) then; ierr=2; return; endif !trap
           else
            segs(j)=dims(j)
           endif
           if(segs(j).le.0) then; ierr=3; return; endif !trap
           vol=vol/dble(segs(j)); n=n-1
          enddo
         elseif(ndim.lt.0) then
          ierr=4
         endif
        else
         ierr=5
        endif
        return
        end subroutine dil_divide_space_int4
!------------------------------------------------------------------------
        subroutine dil_divide_space_int8(ndim,dims,subvol,segs,ierr,algn) !SERIAL
!This subroutine divides an ndim-dimensional block with extents dims(1:ndim)
!into smaller ndim-dimensional blocks with extents segs(1:ndim) such that
!the volume of each such a block approximately equals to <subvol> while
!the ratio Volume/Surface is maximized.
        implicit none
        integer(INTD), intent(in):: ndim                !in: number of dimensions
        integer(INTL), intent(in):: dims(1:*)           !in: dimension extents
        integer(INTL), intent(in):: subvol              !in: target subblock volume
        integer(INTL), intent(inout):: segs(1:*)        !out: dimension segmentation
        integer(INTD), intent(inout):: ierr             !out: error code (0:success)
        integer(INTL), intent(in), optional:: algn(1:*) !in: dimension tiling (default is 1: no tiling)
!---------------------------------------------------------
        logical, parameter:: ARG_CHECK=.true. !argument check
!--------------------------------------------
        integer(INTD):: i,j,n
        integer(INTL):: m,alg(1:ndim),trn(0:ndim)
        real(8):: vol
        logical:: algf

        ierr=0
        if(subvol.gt.0) then
         if(ndim.gt.0) then
          if(ARG_CHECK) then
           do i=1,ndim; if(dims(i).le.0) then; ierr=1; return; endif; enddo
          endif
          if(present(algn)) then
           alg(1:ndim)=algn(1:ndim); algf=.true.
          else
           alg(1:ndim)=1_INTL; algf=.false.
          endif
          trn(0:ndim)=(/+1_INTL,(m,m=1_INTL,int(ndim,INTL))/)
          if(algf) then
           call merge_sort_key_int(int(ndim,INTL),alg,trn)
          else
           call merge_sort_key_int(int(ndim,INTL),dims,trn)
          endif
          n=ndim; vol=dble(subvol)
          do i=ndim,1,-1
           j=trn(i)
           m=nint(vol**(1d0/dble(n))) !uniform segment length (ideal case: cubic)
           if(m.le.dims(j)) then
            segs(j)=max(m-mod(m,alg(j)),alg(j))
            if(mod(segs(j),alg(j)).ne.0) then; ierr=2; return; endif !trap
           else
            segs(j)=dims(j)
           endif
           if(segs(j).le.0) then; ierr=3; return; endif !trap
           vol=vol/dble(segs(j)); n=n-1
          enddo
         elseif(ndim.lt.0) then
          ierr=4
         endif
        else
         ierr=5
        endif
        return
        end subroutine dil_divide_space_int8
!--------------------------------------------------------------------
        subroutine dil_tens_contr_distribute(tcontr,impis,impir,ierr) !SERIAL
!This subroutine splits a global tensor contraction space into subspaces,
!each assigned to an MPI process for further processing.
!tcontr%contr_spec contains the global tensor contraction specification at entrance.
!It will be replaced by an individual tensor contraction specification for process #<impir>.
!If an MPI process does not get any work to do, ierr=DIL_NO_WORK will be returned.
        implicit none
        type(dil_tens_contr_t), intent(inout):: tcontr !inout: full tensor contraction specification (out: %contr_spec)
        integer(INTD), intent(in):: impis              !in: number of MPI processes
        integer(INTD), intent(in):: impir              !in: current MPI process
        integer(INTD), intent(inout):: ierr            !out: error code (0:success; DIL_NO_WORK:no work)
        integer(INTD):: i,j,k,l,m,n,nd,nl,nr,ni
        integer(INTL):: tcvol,sbvol,npieces,ll,mdim(1:3),mseg(1:3)
        integer(INTD):: dn2o(1:MAX_TENSOR_RANK),ln2o(1:MAX_TENSOR_RANK),rn2o(1:MAX_TENSOR_RANK)
        integer(INTD):: tcc(1:MAX_TENSOR_RANK*3),tct(1:MAX_TENSOR_RANK*3),tcu(1:MAX_TENSOR_RANK*3),tcs(1:MAX_TENSOR_RANK*3)
        integer(INTD):: key(1:MAX_TENSOR_RANK*3),trn(0:MAX_TENSOR_RANK*3+1)
        real(8):: tmb,tm,val
        logical:: tiled

        ierr=0
        tmb=thread_wtime()
        if(DIL_DEBUG) write(CONS_OUT,'("#DEBUG(dil_tens_contr_distribute) [",i5,"]: Entered ...")') impir
        if(impis.gt.0.and.impir.ge.0.and.impir.lt.impis) then
!Compute tensor ranks:
         nd=tcontr%contr_spec%ndims_left+tcontr%contr_spec%ndims_right
         nl=tcontr%contr_spec%ndims_contr+tcontr%contr_spec%ndims_left
         nr=tcontr%contr_spec%ndims_contr+tcontr%contr_spec%ndims_right
         ni=tcontr%contr_spec%ndims_contr+tcontr%contr_spec%ndims_left+tcontr%contr_spec%ndims_right
         if(nd.lt.0.or.nl.lt.0.or.nr.lt.0.or.ni.le.0) then; ierr=1; return; endif
!Get N2O (inverse) permutations:
         call permutation_invert(nd,tcontr%contr_spec%dprmn,dn2o,i); if(i.ne.0) then; ierr=2; return; endif
         call permutation_invert(nl,tcontr%contr_spec%lprmn,ln2o,i); if(i.ne.0) then; ierr=3; return; endif
         call permutation_invert(nr,tcontr%contr_spec%rprmn,rn2o,i); if(i.ne.0) then; ierr=4; return; endif
!Compute the global tensor contraction volume:
         tcvol=1_INTL
         do i=1,nl !contracted + left dims
          tcvol=tcvol*tcontr%contr_spec%ldims(i)
         enddo
         do i=1,nr !right dims
          if(tcontr%contr_spec%rprmn(i).gt.tcontr%contr_spec%ndims_contr) then
           tcvol=tcvol*tcontr%contr_spec%rdims(i)
          endif
         enddo
!Collect unique tensor contraction dimensions:
 !Contracted + Left:
         do i=1,nl
          j=tcontr%contr_spec%lprmn(i)
          tcc(j)=tcontr%contr_spec%ldims(i); tct(j)=1_INTD; tiled=.false.
          tcu(j)=IND_NUM_START+tcontr%contr_spec%lbase(i)+tcontr%contr_spec%ldims(i)-1_INTD
          if(j.le.tcontr%contr_spec%ndims_contr) then !contracted position
           if(tcontr%left_arg%store_type.eq.'d'.or.tcontr%left_arg%store_type.eq.'D') then
            if(associated(tcontr%left_arg%tens_distr_p)) then
             tct(j)=tcontr%left_arg%tens_distr_p%tdim(i); tiled=.true.
            else
             ierr=5; return
            endif
           endif
           if(tcontr%right_arg%store_type.eq.'d'.or.tcontr%right_arg%store_type.eq.'D') then
            if(associated(tcontr%right_arg%tens_distr_p)) then
             k=tcontr%right_arg%tens_distr_p%tdim(rn2o(j))
             if(tiled.and.k.ne.tct(j)) then; ierr=6; return; endif !trap: corresponding dimensions have different tiling
             tct(j)=k; tiled=.true.
            else
             ierr=7; return
            endif
           endif
          else !left position
           if(tcontr%left_arg%store_type.eq.'d'.or.tcontr%left_arg%store_type.eq.'D') then
            if(associated(tcontr%left_arg%tens_distr_p)) then
             tct(j)=tcontr%left_arg%tens_distr_p%tdim(i); tiled=.true.
            else
             ierr=8; return
            endif
           endif
           if(tcontr%dest_arg%store_type.eq.'d'.or.tcontr%dest_arg%store_type.eq.'D') then
            if(associated(tcontr%dest_arg%tens_distr_p)) then
             k=tcontr%dest_arg%tens_distr_p%tdim(dn2o(j-tcontr%contr_spec%ndims_contr))
             if(tiled.and.k.ne.tct(j)) then; ierr=9; return; endif !trap: corresponding dimensions have different tiling
             tct(j)=k; tiled=.true.
            else
             ierr=10; return
            endif
           endif
          endif
         enddo
 !Right:
         do i=1,nr
          j=tcontr%contr_spec%rprmn(i)
          if(j.gt.tcontr%contr_spec%ndims_contr) then !right position
           j=j-tcontr%contr_spec%ndims_contr
           tcc(nl+j)=tcontr%contr_spec%rdims(i); tct(nl+j)=1_INTD; tiled=.false.
           tcu(nl+j)=IND_NUM_START+tcontr%contr_spec%rbase(i)+tcontr%contr_spec%rdims(i)-1_INTD
           if(tcontr%right_arg%store_type.eq.'d'.or.tcontr%right_arg%store_type.eq.'D') then
            if(associated(tcontr%right_arg%tens_distr_p)) then
             tct(nl+j)=tcontr%right_arg%tens_distr_p%tdim(i); tiled=.true.
            else
             ierr=11; return
            endif
           endif
           if(tcontr%dest_arg%store_type.eq.'d'.or.tcontr%dest_arg%store_type.eq.'D') then
            if(associated(tcontr%dest_arg%tens_distr_p)) then
             k=tcontr%dest_arg%tens_distr_p%tdim(dn2o(tcontr%contr_spec%ndims_left+j))
             if(tiled.and.k.ne.tct(nl+j)) then; ierr=12; return; endif !trap: corresponding dimensions have different tiling
             tct(nl+j)=k; tiled=.true.
            else
             ierr=13; return
            endif
           endif
          endif
         enddo
         if(DIL_DEBUG) then
          write(CONS_OUT,'(1x,"#DEBUG(DIL): markers:",64(6x,A1))') (/('c',i=1,tcontr%contr_spec%ndims_contr)/),&
          &(/('l',i=1,tcontr%contr_spec%ndims_left)/),(/('r',i=1,tcontr%contr_spec%ndims_right)/) !markers
          write(CONS_OUT,'(1x,"#DEBUG(DIL): dims   :",64(1x,i6))') tcc(1:ni)
          write(CONS_OUT,'(1x,"#DEBUG(DIL): tiling :",64(1x,i6))') tct(1:ni)
          write(CONS_OUT,'(1x,"#DEBUG(dil_tens_contr_distribute) [",i5,"]: Global TC Volume ",i24," for ",i6, " procs.")')&
          &impir,tcvol,impis
         endif
!Compute matrix dimensions:
 !Contracted:
         ll=1_INTL; do i=1,tcontr%contr_spec%ndims_contr; ll=ll*tcc(i); enddo; mdim(1)=ll
 !Left
         ll=1_INTL; do i=1,tcontr%contr_spec%ndims_left; ll=ll*tcc(tcontr%contr_spec%ndims_contr+i); enddo; mdim(2)=ll
 !Right:
         ll=1_INTL; do i=1,tcontr%contr_spec%ndims_right; ll=ll*tcc(nl+i); enddo; mdim(3)=ll
!Divide the matricized 3d tensor contraction space:
         sbvol=tcvol/int(impis,INTL)
         if(sbvol.eq.0) then !not enough work: n is the number of active MPI processes
          n=int(tcvol,INTD); if(impir.lt.n) sbvol=sbvol+1_INTL
         else
          n=impis !all processes are active for now
         endif
         if(impir.lt.n) then !this MPI process will have work to do
          call dil_divide_space_int(3_INTD,mdim,sbvol,mseg,i); if(i.ne.0) then; ierr=14; return; endif
          if(DIL_DEBUG) write(CONS_OUT,'(1x,"#DEBUG(dil_tens_contr_distribute) [",i5,"]: 3d subblock dims:",3(1x,i7))')&
          &impir,mseg(1:3)
!Divide the space encapsulated in each matricized dimension:
 !Contracted:
          k=tcontr%contr_spec%ndims_contr
          if(k.gt.0) then
           call dil_divide_space_int(k,tcc(1:k),mseg(1),tcs(1:k),i,tct(1:k))
           if(i.ne.0) then; ierr=15; return; endif
          endif
 !Left:
          j=tcontr%contr_spec%ndims_left
          if(j.gt.0) then
           call dil_divide_space_int(j,tcc(k+1:k+j),mseg(2),tcs(k+1:k+j),i,tct(k+1:k+j))
           if(i.ne.0) then; ierr=16; return; endif
          endif
 !Right:
          j=tcontr%contr_spec%ndims_right
          if(j.gt.0) then
           call dil_divide_space_int(j,tcc(nl+1:nl+j),mseg(3),tcs(nl+1:nl+j),i,tct(nl+1:nl+j))
           if(i.ne.0) then; ierr=17; return; endif
          endif
          sbvol=1_INTL; do i=1,ni; sbvol=sbvol*tcs(i); enddo !subblock volume
          npieces=1_INTL; do i=1,ni; npieces=npieces*((tcc(i)-1_INTD)/tcs(i)+1_INTD); enddo !number of work pieces
          if(DIL_DEBUG) write(CONS_OUT,'(1x,"#DEBUG(dil_tens_contr_distribute) [",i5,"]: initial subblock segs:",64(1x,i6))')&
          &impir,tcs(1:ni)
          if(DIL_DEBUG) write(CONS_OUT,'(1x,"#DEBUG(dil_tens_contr_distribute) [",i5,"]: initial subblock volume = ",i12)')&
          &impir,sbvol
          if(DIL_DEBUG) write(CONS_OUT,'(1x,"#DEBUG(dil_tens_contr_distribute) [",i5,"]: initial work pieces = ",i12)')&
          &impir,npieces
!Adjust the initial work distribution (cover all space):
          if(npieces.lt.n) then !subblock size may need to be decreased
           trn(0:ni)=(/+1_INTD,(i,i=1_INTD,ni)/)
           do i=1,ni; key(i)=(tcs(i)-1_INTD)/tct(i)+1_INTD; enddo !number of tiling segments per work segment
           call merge_sort_key_int(ni,key,trn)
           do i=ni,1,-1
            j=trn(i)
            if(npieces.lt.n.and.key(j).gt.1) then 
             val=dble(npieces)/dble(n)
             k=int(dble(key(j))*val,INTD)+1_INTD
             npieces=npieces/((tcc(j)-1_INTD)/tcs(j)+1_INTD)
             if(k.gt.1) then
              tcs(j)=min(tct(j)*(k-1_INTD),tcc(j))
              if(npieces*((tcc(j)-1_INTD)/tcs(j)+1_INTD).gt.n) tcs(j)=min(tct(j)*k,tcc(j))
             else
              tcs(j)=min(tct(j),tcc(j))
             endif
             npieces=npieces*((tcc(j)-1_INTD)/tcs(j)+1_INTD)
            else
             exit
            endif
           enddo
           sbvol=1_INTL; do i=1,ni; sbvol=sbvol*tcs(i); enddo !subblock volume
           if(DIL_DEBUG) write(CONS_OUT,'(1x,"#DEBUG(dil_tens_contr_distribute) [",i5,"]: adjusted subblock segs:",64(1x,i6))')&
           &impir,tcs(1:ni)
           if(DIL_DEBUG) write(CONS_OUT,'(1x,"#DEBUG(dil_tens_contr_distribute) [",i5,"]: adjusted subblock volume = ",i12)')&
           &impir,sbvol
           if(DIL_DEBUG) write(CONS_OUT,'(1x,"#DEBUG(dil_tens_contr_distribute) [",i5,"]: adjusted work pieces = ",i12)')&
           &impir,npieces
          endif
          l=0
          do while(npieces.gt.n) !subblock size may need to be increased
           l=l+1
           if(l.gt.5) then !trap
            if(VERBOSE)&
            &write(CONS_OUT,'("#ERROR(dil_tens_contr_distribute)[",i5,"]: Unable to properly increase the subblock size!")') impir
            ierr=18; return
           endif
           trn(0:ni)=(/+1_INTD,(i,i=1_INTD,ni)/)
           do i=1,ni; key(i)=(tcs(i)-1_INTD)/tct(i)+1_INTD; enddo !number of tiling segments per work segment
           call merge_sort_key_int(ni,key,trn) 
           do i=ni,1,-1
            if(npieces.gt.n) then
             j=trn(i)
             if(tcs(j).lt.tcc(j)) then
              val=dble(npieces)/dble(n)
              k=int(dble(key(j))*val,INTD)+l
              npieces=npieces/((tcc(j)-1_INTD)/tcs(j)+1_INTD)
              tcs(j)=min(tct(j)*k,tcc(j))
              npieces=npieces*((tcc(j)-1_INTD)/tcs(j)+1_INTD)
             endif
            else
             exit
            endif
           enddo
          enddo
          sbvol=1_INTL; do i=1,ni; sbvol=sbvol*tcs(i); enddo !subblock volume
          if(DIL_DEBUG) write(CONS_OUT,'(1x,"#DEBUG(dil_tens_contr_distribute) [",i5,"]: adjusted subblock segs:",64(1x,i6))')&
          &impir,tcs(1:ni)
          if(DIL_DEBUG) write(CONS_OUT,'(1x,"#DEBUG(dil_tens_contr_distribute) [",i5,"]: adjusted subblock volume = ",i12)')&
          &impir,sbvol
          if(DIL_DEBUG) write(CONS_OUT,'(1x,"#DEBUG(dil_tens_contr_distribute) [",i5,"]: adjusted work pieces = ",i12)')&
          &impir,npieces
          if(DIL_DEBUG) write(CONS_OUT,'(1x,"#DEBUG(dil_tens_contr_distribute) [",i5,"]: work starving ratio (>=1) = ",F10.4)')&
          &impir,dble(n)/dble(npieces)
!Distribute the work subspaces among MPI processes:
 !Compute the starting point in the tensor contraction space for each MPI process:
          do i=1,ni; trn(i)=(tcc(i)-1_INTD)/tcs(i)+1_INTD; enddo !number of work segments per dimension
          trn(0)=1_INTD; do i=1,ni; trn(i)=trn(i)*trn(i-1); enddo !division bases: trn(ni)=npieces
          if(DIL_DEBUG) write(CONS_OUT,'(1x,"#DEBUG(dil_tens_contr_distribute) [",i5,"]: ",i12," work pieces for ",i6," procs.")')&
          &impir,trn(ni),n
          if(trn(ni).eq.npieces) then !number of work pieces <= number of active MPI processes
  !Coarse-grain setup:
           j=mod(impir,trn(ni))
           do i=ni,1,-1
            key(i)=(j/trn(i-1))*tcs(i); j=mod(j,trn(i-1)) !key(1:ni): starting offsets for each MPI process
           enddo
   !Substitute the global tensor contraction specification with the individual ones (for each MPI process):
    !Replace contracted dimensions:
           k=0
           do i=1,tcontr%contr_spec%ndims_contr
            l=ln2o(i); tcontr%contr_spec%lbase(l)=tcontr%contr_spec%lbase(l)+key(k+i)
            tcontr%contr_spec%ldims(l)=min(tcs(k+i),tcu(k+i)-(IND_NUM_START+tcontr%contr_spec%lbase(l))+1_INTD)
            j=rn2o(i); tcontr%contr_spec%rbase(j)=tcontr%contr_spec%rbase(j)+key(k+i)
            tcontr%contr_spec%rdims(j)=tcontr%contr_spec%ldims(l)
           enddo
    !Replace left dimensions:
           k=tcontr%contr_spec%ndims_contr
           do i=1,tcontr%contr_spec%ndims_left
            l=dn2o(i); tcontr%contr_spec%dbase(l)=tcontr%contr_spec%dbase(l)+key(k+i)
            tcontr%contr_spec%ddims(l)=min(tcs(k+i),tcu(k+i)-(IND_NUM_START+tcontr%contr_spec%dbase(l))+1_INTD)
            j=ln2o(k+i); tcontr%contr_spec%lbase(j)=tcontr%contr_spec%lbase(j)+key(k+i)
            tcontr%contr_spec%ldims(j)=tcontr%contr_spec%ddims(l)
           enddo
    !Replace right dimensions:
           k=nl
           do i=1,tcontr%contr_spec%ndims_right
            l=dn2o(tcontr%contr_spec%ndims_left+i); tcontr%contr_spec%dbase(l)=tcontr%contr_spec%dbase(l)+key(k+i)
            tcontr%contr_spec%ddims(l)=min(tcs(k+i),tcu(k+i)-(IND_NUM_START+tcontr%contr_spec%dbase(l))+1_INTD)
            j=rn2o(tcontr%contr_spec%ndims_contr+i); tcontr%contr_spec%rbase(j)=tcontr%contr_spec%rbase(j)+key(k+i)
            tcontr%contr_spec%rdims(j)=tcontr%contr_spec%ddims(l)
           enddo
  !Fine-grain setup of individual tensor contraction specifications:
           m=n/trn(ni); if(mod(impir,trn(ni)).lt.mod(n,trn(ni))) m=m+1 !m is the number of parts the subblock should be split into
           if(m.gt.1) then
            k=impir/trn(ni) !part number assigned to this MPI process
            if(k.le.1) then !2-part splitting at most
             m=2 !force 2-part splitting
             tiled=.false.
   !Try to split a right dimension into two parts:
             do i=1,tcontr%contr_spec%ndims_right
              j=rn2o(tcontr%contr_spec%ndims_contr+i)
              if(tcontr%contr_spec%rdims(j).ge.tct(nl+i)*2) then
               l=tcontr%contr_spec%rdims(j)/2_INTD; l=l-mod(l,tct(nl+i))
               if(k.eq.0) then
                tcontr%contr_spec%rdims(j)=l
               elseif(k.eq.1) then
                tcontr%contr_spec%rbase(j)=tcontr%contr_spec%rbase(j)+l
                tcontr%contr_spec%rdims(j)=tcontr%contr_spec%rdims(j)-l
               endif
               l=dn2o(tcontr%contr_spec%ndims_left+i)
               tcontr%contr_spec%dbase(l)=tcontr%contr_spec%rbase(j)
               tcontr%contr_spec%ddims(l)=tcontr%contr_spec%rdims(j)
               tiled=.true.; exit
              endif
             enddo
   !Try to split either contracted or left dimension into two parts:
             if(.not.tiled) then
              do i=1,nl
               j=ln2o(i)
               if(tcontr%contr_spec%ldims(j).ge.tct(i)*2) then
                l=tcontr%contr_spec%ldims(j)/2_INTD; l=l-mod(l,tct(i))
                if(k.eq.0) then
                 tcontr%contr_spec%ldims(j)=l
                elseif(k.eq.1) then
                 tcontr%contr_spec%lbase(j)=tcontr%contr_spec%lbase(j)+l
                 tcontr%contr_spec%ldims(j)=tcontr%contr_spec%ldims(j)-l
                endif
                if(i.le.tcontr%contr_spec%ndims_contr) then
                 l=rn2o(i)
                 tcontr%contr_spec%rbase(l)=tcontr%contr_spec%lbase(j)
                 tcontr%contr_spec%rdims(l)=tcontr%contr_spec%ldims(j)
                else
                 l=dn2o(i-tcontr%contr_spec%ndims_contr)
                 tcontr%contr_spec%dbase(l)=tcontr%contr_spec%lbase(j)
                 tcontr%contr_spec%ddims(l)=tcontr%contr_spec%ldims(j)
                endif
                tiled=.true.; exit
               endif
              enddo
             endif
            else
             tcontr%contr_spec%ddims(1:nd)=0; tcontr%contr_spec%ldims(1:nl)=0; tcontr%contr_spec%rdims(1:nr)=0
             ierr=DIL_NO_WORK !flag: NO WORK for this MPI process (#impir)
             if(DIL_DEBUG) write(CONS_OUT,'(1x,"#DEBUG(dil_tens_contr_distribute) [",i5,"]: 3d subblock dims:",3(1x,i7))')&
             &impir,0,0,0
            endif
           endif
          else
           if(VERBOSE) write(CONS_OUT,'(1x,"#ERROR(dil_tens_contr_distribute) [",i5,"]: Work pieces number mismatch:",3(1x,i12))')&
           &impir,n,trn(ni),npieces
           ierr=19; return
          endif
         else !this MPI process will have no work
          tcontr%contr_spec%ddims(1:nd)=0; tcontr%contr_spec%ldims(1:nl)=0; tcontr%contr_spec%rdims(1:nr)=0
          ierr=DIL_NO_WORK !flag: NO WORK for this MPI process (#impir)
          if(DIL_DEBUG) write(CONS_OUT,'(1x,"#DEBUG(dil_tens_contr_distribute) [",i5,"]: 3d subblock dims:",3(1x,i7))')&
          &impir,0,0,0
         endif
        else
         ierr=20
        endif
        tm=thread_wtime(tmb)
        if(DIL_DEBUG) write(CONS_OUT,'("#DEBUG(dil_tens_contr_distribute) [",i5,"]: Exited with status ",i9," ( ",F10.4," s)")')&
        &impir,ierr,tm
        return
        end subroutine dil_tens_contr_distribute
!------------------------------------------------------------------------------------------------------
        subroutine dil_tens_contr_partition(cspec,max_arg_vol,contr_case,task_list,ierr,darr,larr,rarr) !SERIAL
!This subroutine partitions a distributed tensor contraction into smaller parts that will fit into local RAM.
!If any of the tensor contraction arguments is distributed tiled, it must be explicitly passed here.
!IMPORTANT: All distributed tiled (sub)tensors in the tensor contraction specification
!must be aligned to the natural dimension segmentation boundaries (tiling boundaries)!
!NOTES: The arguments of the tensor contraction are split into parts such that
!the corresponding matrices are close to squares as possible (leading to max compute intensity).
        implicit none
        type(contr_spec_t), intent(in):: cspec             !in: tensor contraction specification
        integer(INTL), intent(in):: max_arg_vol            !in: local memory limit per argument (real words)
        character(3), intent(in):: contr_case              !in: tensor contraction case: {lll,lld,dll,ldd,dld,ddd}
        type(contr_task_list_t), intent(inout):: task_list !out: task list
        integer(INTD), intent(inout):: ierr                !out: error code (0:success)
        type(tensor), optional, intent(in):: darr          !in: destination tensor argument (if distributed tiled)
        type(tensor), optional, intent(in):: larr          !in: left tensor argument (if distributed tiled)
        type(tensor), optional, intent(in):: rarr          !in: right tensor argument (if distributed tiled)
!------------------------------------------------
        logical, parameter:: NO_CHECK=.false. !argument check
!------------------------------------------------
        integer(INTD):: i,j,k,l,m,n,nd,nl,nr,ni,ntasks,impir !,index_code,index_host,index_pos
        integer(INTD):: prmn(1:MAX_TENSOR_RANK,0:2),lb(1:MAX_TENSOR_RANK*3),ub(1:MAX_TENSOR_RANK*3),sb(1:MAX_TENSOR_RANK*3)
        integer(INTD):: ip(1:MAX_TENSOR_RANK*3),im(1:MAX_TENSOR_RANK*3),ts(1:MAX_TENSOR_RANK*3),ns(1:MAX_TENSOR_RANK*3)
        integer(INTD):: cl(1:MAX_TENSOR_RANK),cu(1:MAX_TENSOR_RANK)
        integer(INTL):: lld,lrd,lcd,lsm,ll,lr,lc,ml,mr,mc
        logical:: more_tasks,arg_loc(0:2)
        real(8):: tmb,tm

!        index_code(i,j)=i*MAX_TENSOR_RANK+j          !{host{0,1,2},position{1..MAX_TENSOR_RANK}} --> index code
!        index_host(i)=(abs(i)-1)/MAX_TENSOR_RANK     !index code --> host argument: {0,1,2}={d,l,r}
!        index_pos(i)=mod(abs(i)-1,MAX_TENSOR_RANK)+1 !index code --> index (dimension) position

        ierr=0; tmb=thread_wtime(); impir=0
        impir=my_mpi_rank(infpar%lg_comm)
        if(DIL_DEBUG) write(CONS_OUT,'("#DEBUG(tensor_algebra_dil::dil_tens_contr_partition)[",i2,"]: Entered ...")') impir !debug
!Argument check:
        do i=1,3
         if(contr_case(i:i).eq.'l'.or.contr_case(i:i).eq.'L') then
          arg_loc(i-1)=.true.
         elseif(contr_case(i:i).eq.'d'.or.contr_case(i:i).eq.'D') then
          arg_loc(i-1)=.false.
         else
          ierr=1; return
         endif
        enddo
        if(.not.NO_CHECK) then
         if((.not.arg_loc(0)).and.(.not.present(darr))) then; ierr=2; return; endif
         if((.not.arg_loc(1)).and.(.not.present(larr))) then; ierr=3; return; endif
         if((.not.arg_loc(2)).and.(.not.present(rarr))) then; ierr=4; return; endif
         i=dil_tens_contr_spec_check(cspec); if(i.ne.0) then; ierr=5; return; endif
        endif
!Get (sub)tensor ranks:
        nd=cspec%ndims_left+cspec%ndims_right  !destination (sub)tensor rank
        nl=cspec%ndims_contr+cspec%ndims_left  !left (sub)tensor rank
        nr=cspec%ndims_contr+cspec%ndims_right !right (sub)tensor rank
!Get inverse permutations:
        call permutation_invert(nd,cspec%dprmn,prmn(1:,0),i); if(i.ne.0) then; ierr=6; return; endif !N2O for the destination tensor
        call permutation_invert(nl,cspec%lprmn,prmn(1:,1),i); if(i.ne.0) then; ierr=7; return; endif !N2O for left tensor
        call permutation_invert(nr,cspec%rprmn,prmn(1:,2),i); if(i.ne.0) then; ierr=8; return; endif !N2O for right tensor
!Clean the task list:
        call dil_contr_task_list_destroy(task_list,i); if(i.ne.0) then; ierr=9; return; endif
!Init the global index sequence:
        ni=cspec%ndims_contr+cspec%ndims_left+cspec%ndims_right !total number of unique dimensions in the tensor contraction
        if(ni.le.0) then; ierr=10; return; endif !tensor contraction arguments must not be all scalars
        lld=1_INTL; lrd=1_INTL; lcd=1_INTL; ml=1_INTL; mr=1_INTL; mc=1_INTL; l=0
        do i=nl,1,-1 !old number in left tensor
         j=cspec%lprmn(i)
         if(j.le.cspec%ndims_contr) then !contracted multi-index
          l=l+1; lb(l)=IND_NUM_START+cspec%lbase(i); ub(l)=lb(l)+cspec%ldims(i)-1; ip(l)=j
          if(.not.arg_loc(1)) then
           sb(l)=larr%tdim(i) !segmented dim
           mc=mc*sb(l)
          else
           if(.not.arg_loc(2)) then
            sb(l)=rarr%tdim(prmn(j,2)) !segmented dim
            mc=mc*sb(l)
           else
            sb(l)=-cspec%ldims(i) !non-segmented dim (mark it negative)
            mc=mc*MIN_LOC_DIM_EXT
           endif
          endif
          lcd=lcd*cspec%ldims(i)
         endif
        enddo
        k=cspec%ndims_contr; l=nl
        do i=nd,1,-1 !old number in destination tensor
         j=cspec%dprmn(i)
         if(j.le.cspec%ndims_left) then !left multi-index
          k=k+1; lb(k)=IND_NUM_START+cspec%dbase(i); ub(k)=lb(k)+cspec%ddims(i)-1; ip(k)=j
          if(.not.arg_loc(0)) then
           sb(k)=darr%tdim(i) !segmented dim
           ml=ml*sb(k)
          else
           if(.not.arg_loc(1)) then
            sb(k)=larr%tdim(prmn(cspec%ndims_contr+ip(k),1)) !segmented dim
            ml=ml*sb(k)
           else
            sb(k)=-cspec%ddims(i) !non-segmented dim (mark it negative)
            ml=ml*MIN_LOC_DIM_EXT
           endif
          endif
          lld=lld*cspec%ddims(i)
         else !right multi-index
          l=l+1; lb(l)=IND_NUM_START+cspec%dbase(i); ub(l)=lb(l)+cspec%ddims(i)-1; ip(l)=j-cspec%ndims_left
          if(.not.arg_loc(0)) then
           sb(l)=darr%tdim(i) !segmented dim
           mr=mr*sb(l)
          else
           if(.not.arg_loc(2)) then
            sb(l)=rarr%tdim(prmn(cspec%ndims_contr+ip(l),2)) !segmented dim
            mr=mr*sb(l)
           else
            sb(l)=-cspec%ddims(i) !non-segmented dim (mark it negative)
            mr=mr*MIN_LOC_DIM_EXT
           endif
          endif
          lrd=lrd*cspec%ddims(i)
         endif
        enddo
!Proceed:
        lsm=int(dsqrt(real(max_arg_vol,8)),INTL) !dimension of the largest square matrix fitting in the buffer
 !Decide on splitting matrix dimensions:
        lc=lcd; ll=lld; lr=lrd
        if(lcd.le.lld) then !lcd<=lld
         if(lld.le.lrd) then !*lcd<=lld<=lrd*
          if(lld*lrd.gt.lsm*lsm) then
           if(lld.gt.lsm.and.lrd.gt.lsm) then
            ll=lsm; lr=lsm; lc=min(lcd,lsm)
           elseif(lld.le.lsm.and.lrd.gt.lsm) then
            lr=lsm*lsm/lld; ll=lld; lc=lcd
           endif
          endif
         else !lrd<lld
          if(lrd.le.lcd) then !*lrd<=lcd<=lld*
           if(lcd*lld.gt.lsm*lsm) then
            if(lcd.gt.lsm.and.lld.gt.lsm) then
             lc=lsm; ll=lsm; lr=min(lrd,lsm)
            elseif(lcd.le.lsm.and.lld.gt.lsm) then
             ll=lsm*lsm/lcd; lc=lcd; lr=lrd
            endif
           endif
          else !*lcd<lrd<lld*
           if(lrd*lld.gt.lsm*lsm) then
            if(lrd.gt.lsm.and.lld.gt.lsm) then
             lr=lsm; ll=lsm; lc=min(lcd,lsm)
            elseif(lrd.le.lsm.and.lld.gt.lsm) then
             ll=lsm*lsm/lrd; lr=lrd; lc=lcd
            endif
           endif
          endif
         endif
        else !lld<lcd
         if(lrd.lt.lld) then !*lrd<lld<lcd*
          if(lld*lcd.gt.lsm*lsm) then
           if(lld.gt.lsm.and.lcd.gt.lsm) then
            ll=lsm; lc=lsm; lr=min(lrd,lsm)
           elseif(lld.le.lsm.and.lcd.gt.lsm) then
            lc=lsm*lsm/lld; ll=lld; lr=lrd
           endif
          endif
         else !lld<=lrd
          if(lcd.lt.lrd) then !*lld<lcd<lrd*
           if(lcd*lrd.gt.lsm*lsm) then
            if(lcd.gt.lsm.and.lrd.gt.lsm) then
             lc=lsm; lr=lsm; ll=min(lld,lsm)
            elseif(lcd.le.lsm.and.lrd.gt.lsm) then
             lr=lsm*lsm/lcd; lc=lcd; ll=lld
            endif
           endif
          else !*lld<=lrd<=lcd*
           if(lrd*lcd.gt.lsm*lsm) then
            if(lrd.gt.lsm.and.lcd.gt.lsm) then
             lr=lsm; lc=lsm; ll=min(lld,lsm)
            elseif(lrd.le.lsm.and.lcd.gt.lsm) then
             lc=lsm*lsm/lrd; lr=lrd; ll=lld
            endif
           endif
          endif
         endif
        endif
        if(lr.lt.mr.and.lrd.ge.mr) then !minimal possible part volume for right dims
         lr=mr
         if(lr*ll.gt.max_arg_vol) ll=int(dble(max_arg_vol)/dble(lr))
         if(lr*lc.gt.max_arg_vol) lc=int(dble(max_arg_vol)/dble(lr))
        endif
        if(ll.lt.ml.and.lld.ge.ml) then !minimal possible part volume for left dims
         ll=ml
         if(ll*lr.gt.max_arg_vol) lr=int(dble(max_arg_vol)/dble(ll))
         if(ll*lc.gt.max_arg_vol) lc=int(dble(max_arg_vol)/dble(ll))
        endif
        if(lc.lt.mc.and.lcd.ge.mc) then !minimal possible part volume for contracted dims
         lc=mc
         if(lc*lr.gt.max_arg_vol) lr=int(dble(max_arg_vol)/dble(lc))
         if(lc*ll.gt.max_arg_vol) ll=int(dble(max_arg_vol)/dble(lc))
        endif
        if(DIL_DEBUG) then
         write(CONS_OUT,'("#DEBUG(tensor_algebra_dil::dil_tens_contr_partition)[",i2,"]: Max Vol/Dims: ",i11,3(1x,i8))')&
         &impir,max_arg_vol,ll,lr,lc !debug
        endif
        if(lc.le.0_INTL.or.ll.le.0_INTL.or.lr.le.0_INTL) then; ierr=11; return; endif
        if(lc.gt.lcd.or.ll.gt.lld.or.lr.gt.lrd) then; ierr=12; return; endif
        if(ll*lr.gt.max_arg_vol.or.lc*ll.gt.max_arg_vol.or.lc*lr.gt.max_arg_vol) then
         if(VERBOSE) write(CONS_OUT,'("#ERROR(tensor_algebra_dil::dil_tens_contr_partition)[",i2'//&
         &',"]: partitioning logic failed, contact Dmitry Lyakh!")') impir
         ierr=13; return
        endif
 !Split (sub)tensor dimensions:
        i=split_argument(1_INTD,cspec%ndims_contr,lc) !contracted dims
        if(i.ne.0) then
         if(VERBOSE) write(CONS_OUT,'("#ERROR(tensor_algebra_dil::dil_tens_contr_partition)[",i2,"]: Argument split failed: ",i9)')&
         &impir,i
         ierr=14; return
        endif
        i=split_argument(cspec%ndims_contr+1_INTD,nl,ll) !left dims
        if(i.ne.0) then
         if(VERBOSE) write(CONS_OUT,'("#ERROR(tensor_algebra_dil::dil_tens_contr_partition)[",i2,"]: Argument split failed: ",i9)')&
         &impir,i
         ierr=15; return
        endif
        i=split_argument(nl+1_INTD,ni,lr) !right dims
        if(i.ne.0) then
         if(VERBOSE) write(CONS_OUT,'("#ERROR(tensor_algebra_dil::dil_tens_contr_partition)[",i2,"]: Argument split failed: ",i9)')&
         &impir,i
         ierr=16; return
        endif
        if(DIL_DEBUG) then
         write(CONS_OUT,'(1x,"#DEBUG(DIL): XX:",64(6x,A1))') (/('c',i=1,cspec%ndims_contr)/),&
         &(/('l',i=1,cspec%ndims_left)/),(/('r',i=1,cspec%ndims_right)/) !markers
         write(CONS_OUT,'(1x,"#DEBUG(DIL): LB:",64(1x,i6))') lb(1:ni) !index lower bounds
         write(CONS_OUT,'(1x,"#DEBUG(DIL): UB:",64(1x,i6))') ub(1:ni) !index upper bounds
         write(CONS_OUT,'(1x,"#DEBUG(DIL): SB:",64(1x,i6))') sb(1:ni) !index natural segments (for storage)
         write(CONS_OUT,'(1x,"#DEBUG(DIL): TS:",64(1x,i6))') ts(1:ni) !index segments for work partitioning
        endif
 !Check segmentation:
        lsm=1_INTL; do i=1,nl; lsm=lsm*ts(i); enddo
        if(lsm.gt.max_arg_vol) then; ierr=17; return; endif
        lsm=1_INTL; do i=cspec%ndims_contr+1,ni; lsm=lsm*ts(i); enddo
        if(lsm.gt.max_arg_vol) then; ierr=18; return; endif
        lsm=1_INTL; do i=1,cspec%ndims_contr; lsm=lsm*ts(i); enddo; do i=nl+1,ni; lsm=lsm*ts(i); enddo
        if(lsm.gt.max_arg_vol) then; ierr=19; return; endif
 !Generate the task list:
        ntasks=1; do i=1,ni; ntasks=ntasks*ns(i); enddo !total number of tasks
        call dil_contr_task_list_create(task_list,ntasks,i); if(i.ne.0) then; ierr=20; return; endif        
        im(1:ni)=lb(1:ni); more_tasks=.true.; n=1
        do while(more_tasks) !generate tensor contraction tasks operating on (sub)tensor parts
  !Set destination (sub)tensor part:
         do i=cspec%ndims_contr+1,nl
          j=ip(i); k=prmn(j,0); cl(k)=im(i); cu(k)=min(im(i)+ts(i)-1,ub(i))
         enddo
         do i=nl+1,ni
          j=ip(i); k=prmn(cspec%ndims_left+j,0); cl(k)=im(i); cu(k)=min(im(i)+ts(i)-1,ub(i))
         enddo
         call dil_contr_task_set_arg(task_list%contr_tasks(n),'d',nd,cl,cu,i); if(i.ne.0) then; ierr=21; return; endif
  !Set left (sub)tensor part:
         do i=1,cspec%ndims_contr
          j=ip(i); k=prmn(j,1); cl(k)=im(i); cu(k)=min(im(i)+ts(i)-1,ub(i))
         enddo
         do i=cspec%ndims_contr+1,nl
          j=ip(i); k=prmn(cspec%ndims_contr+j,1); cl(k)=im(i); cu(k)=min(im(i)+ts(i)-1,ub(i))
         enddo
         call dil_contr_task_set_arg(task_list%contr_tasks(n),'l',nl,cl,cu,i); if(i.ne.0) then; ierr=22; return; endif
  !Set right (sub)tensor part:
         do i=1,cspec%ndims_contr
          j=ip(i); k=prmn(j,2); cl(k)=im(i); cu(k)=min(im(i)+ts(i)-1,ub(i))
         enddo
         do i=nl+1,ni
          j=ip(i); k=prmn(cspec%ndims_contr+j,2); cl(k)=im(i); cu(k)=min(im(i)+ts(i)-1,ub(i))
         enddo
         call dil_contr_task_set_arg(task_list%contr_tasks(n),'r',nr,cl,cu,i); if(i.ne.0) then; ierr=23; return; endif
  !Count task Flops:
         call dil_contr_task_set_flops(cspec,task_list%contr_tasks(n),i); if(i.ne.0) then; ierr=24; return; endif
         if(DIL_DEBUG) write(CONS_OUT,'("#DEBUG(tensor_algebra_dil::dil_tens_contr_partition): generated task # ",i9)') n !debug
         if(DIL_DEBUG) call dil_contr_task_print(cspec,task_list%contr_tasks(n)) !debug
         more_tasks=get_next_mlndx(n) !next task
        enddo
        if(n.ne.ntasks) then
         if(VERBOSE) write(CONS_OUT,'("#ERROR(tensor_algebra_dil::dil_tens_contr_partition): task number mismatch: ",'//&
         &'i9,1x,i9)') ntasks,n
         ierr=25; return
        endif
        tm=thread_wtime(tmb)
        if(DIL_DEBUG) write(CONS_OUT,'("#DEBUG(tensor_algebra_dil::dil_tens_contr_partition)[",i2,"]: Done in ",F10.4," s.")')&
        &impir,tm !debug
        return

        contains

         integer(INTD) function split_argument(ds,df,bnd) !splits (sub)tensor dimensions to fit it into the buffer
         implicit none                   !NOTE: indices {ds:df} must belong to the same index group (contracted, left, right)
         integer(INTD), intent(in):: ds  !first index of the group: {contracted,left,right}
         integer(INTD), intent(in):: df  !last index of the group: {contracted,left,right}
         integer(INTL), intent(in):: bnd !upper bound for the (sub)tensor part volume (real words)
         integer(INTL):: pvl
         integer(INTD):: j0,je
         split_argument=0; pvl=1_INTL
 !Compute the minimal volume of the tensor part:
         if(sb(ds).gt.0) then !index group consists of segmented indices --> minimal volume exists
          do j0=ds,df
           if(sb(j0)*sb(ds).gt.0) then !all indices must be segmented
            ts(j0)=min(sb(j0),ub(j0)-lb(j0)+1_INTD); pvl=pvl*ts(j0)
           else
            split_argument=1; return
           endif
          enddo
         else !index group consists of non-segmented indices --> set some minimal volume
          do j0=ds,df
           if(sb(j0)*sb(ds).gt.0) then !all indices must be non-segmented
            ts(j0)=min(MIN_LOC_DIM_EXT,ub(j0)-lb(j0)+1_INTD); pvl=pvl*ts(j0)
           else
            split_argument=2; return
           endif
          enddo
         endif
 !Increase the volume of the (sub)tensor part, if possible:
         if(pvl.le.bnd) then
          do j0=df,ds,-1
           je=ub(j0)-lb(j0)+1_INTD; pvl=pvl/ts(j0)
           if(pvl*je.gt.bnd) then
            if(sb(j0).gt.0) then !segmented dim
             ts(j0)=min((bnd/(pvl*sb(j0)))*sb(j0),je)
            else !non-segmented dim
             ts(j0)=min((bnd/(pvl*MIN_LOC_DIM_EXT))*MIN_LOC_DIM_EXT,je)
            endif
            pvl=pvl*ts(j0)
            exit
           else
            ts(j0)=je; pvl=pvl*ts(j0)
           endif
          enddo
         endif
 !Compute <ns(:)>, number of supersegments per dimension:
         do j0=ds,df; ns(j0)=(ub(j0)-lb(j0))/ts(j0)+1_INTD; enddo
         return
         end function split_argument

         logical function get_next_mlndx(jn) !multi-index iterator
         implicit none                       !imports: ni,im,ip,lb,ub,ts
         integer(INTD), intent(inout):: jn   !inout: sequential number of the multi-index
         integer(INTD):: jj
         get_next_mlndx=.false.; jj=1
         do while(jj.le.ni)
          im(jj)=im(jj)+ts(jj)
          if(im(jj).le.ub(jj)) then
           get_next_mlndx=.true.; jn=jn+1; exit
          else
           im(jj)=lb(jj); jj=jj+1
          endif
         enddo
         return
         end function get_next_mlndx

        end subroutine dil_tens_contr_partition
!----------------------------------------------------------------------------------------------
        subroutine dil_tensor_contract_pipe(cspec,darg,larg,rarg,alpha,beta,mem_lim,hbuf,ierr,&
                                           &locked,nasync,lasync,num_gpus,num_mics)             !PARALLEL (MPI+OMP+CUDA+MIC)
!This subroutine implements pipelined tensor contractions for CPU/GPU/MIC.
!Details:
! * Each Device is assigned several buffers:
!   (a) 6 argument buffers in Device RAM (maximally large): 3 prefetch buffers + 3 compute buffers;
!   (b) 1 upload MPI buffer in Host RAM, same size as (a);
!   (c) Accelerators only (No GPU Direct): 3 smaller MPI buffers in Host RAM (should fit the largest tile at least).
!Notes:
! * For GPU pipelining, <hbuf> must be pinned and its starting address must be properly aligned!
! * There is no MPI barriers in this subroutine: Each process performs its own work and returns independently.
! * If <locked> is present and is .TRUE., MPI windows for distributed tensor arguments are assumed locked,
!   such that MPI_WIN_FLUSH operation will be used for progressing the RMA. Otherwise, MPI_WIN_LOCK/UNLOCK will be used.
! * buf_conf(BUFS_PER_DEV,0:MAX_DEVS-1) is used for dynamic switching between available device buffers such that
!   (a) Argument load (prefetch) is always done into buffers referred to as 1-3 (D,L,R);
!   (b) Computation always uses buffers referred to as 4-6 (D,L,R);
!   (c) Result upload is always done from buffer referred to as 7.
! * If <nasync> and <lasync> are present, the last series of MPI uploads will not be finalized here,
!   but postponed for later: lasync(1:nasync) will contain the needed information.
        implicit none
        type(contr_spec_t), intent(in):: cspec                    !in: tensor contraction specification
        type(tens_arg_t), intent(inout):: darg                    !inout: destination tensor argument
        type(tens_arg_t), intent(in):: larg                       !in: left tensor argument
        type(tens_arg_t), intent(in):: rarg                       !in: right tensor argument
        real(tensor_dp), intent(in):: alpha                           !in: tensor contraction prefactor (GEMM alpha)
        real(tensor_dp), intent(in):: beta                            !in: GEMM beta
        integer(INTL), intent(in):: mem_lim                       !in: local memory limit (bytes): buffer space
        real(tensor_dp), intent(inout), target, contiguous:: hbuf(1:) !inout: existing external buffer (to avoid mem allocations)
        integer(INTD), intent(inout):: ierr                       !out: error code (0:success)
        logical, intent(in), optional:: locked                    !in: if .true., MPI wins for tens-args are assumed locked
        integer(INTD), intent(out), optional:: nasync             !out: number of outstanding async MPI uploads
        type(rank_win_t), intent(inout), optional:: lasync(1:)    !out: outstanding async MPI uploads
        integer(INTD), intent(in), optional:: num_gpus            !in: number of NVidia GPUs available on the node: (0..max)
        integer(INTD), intent(in), optional:: num_mics            !in: number of Intel MICs available on the node: (0..max)
!-----------------------------------------------------
        logical, parameter:: NO_CHECK=.false.         !argument check
        logical, parameter:: PREP_AND_COMM=.false.    !communication/tensor_preparation overlap
        logical, parameter:: ARGS_REUSE=.true.        !argument reuse in tensor contractions
        logical, parameter:: TASK_RESHUFFLE=.true.    !task reshuffling (to reduce the number of MPI collisions)
!-------------------------------------------------
        integer(INTD):: i,j,k,l,m,n
        type(dev_buf_t):: buf(0:MAX_DEVS-1)         !Host buffers for all devices (mapped to the Host buffer space)
        type(contr_task_list_t), target:: task_list !`Make it global threadsafe to allow reuse and avoid unnecessary allocations
        character(3):: contr_case,arg_reuse,arg_keep
        integer(INTL):: i0,i1,i2,size_of_real,tot_buf_vol,dev_buf_vol,arg_buf_vol,dvol,lvol,rvol
        integer(INTD):: impir,nd,nl,nr,num_dev,dev,buf_conf(1:BUFS_PER_DEV,0:MAX_DEVS-1),prmn(1:MAX_TENSOR_RANK,0:2)
        integer(INTD):: tasks_done(0:MAX_DEVS-1),task_curr(0:MAX_DEVS-1),task_prev(0:MAX_DEVS-1),task_next(0:MAX_DEVS-1)
        integer(INTD):: first_avail(0:MAX_DEVS-1)
        logical:: gpu_on,mic_on,first_task,next_load,prev_store,args_here
        logical:: win_lck,err_curr,err_prev,err_next,triv_d,triv_l,triv_r,async
        real(8):: tmb,tms,tm,tmm,tc_flops,mm_flops
        real(tensor_dp):: val

        ierr=0; tmb=thread_wtime()
        impir=my_mpi_rank(infpar%lg_comm) !rank in local MPI communicator
        if(DIL_DEBUG) write(CONS_OUT,'("#DEBUG(tensor_algebra_dil::dil_tensor_contract_pipe)[",i2,"]: Entered ...")') impir !debug
!Init:
        val=0E0_tensor_dp; size_of_real=sizeof(val)
        tc_flops=0d0; mm_flops=0d0; tmm=0d0
        num_dev=1; gpu_on=.false.; mic_on=.false.
        async=.false.
        if(present(nasync)) then
         if(present(lasync)) then
          nasync=0; async=.true.
         else
          call cleanup(1_INTD); return
         endif
        else
         if(present(lasync)) then; call cleanup(2_INTD); return; endif
        endif
        if(present(locked)) then; win_lck=locked; else; win_lck=.false.; endif
        if(DIL_DEBUG) write(CONS_OUT,'("#DEBUG(tensor_algebra_dil::dil_tensor_contract_pipe)[",i2,"]: OUTSIDE LOCK = ",l1'//&
        &',": ASYNC = ",l1)') impir,win_lck,async !debug
        contr_case=darg%store_type//larg%store_type//rarg%store_type !contraction case
!Check input arguments:
        if(mem_lim.lt.MIN_BUF_MEM) then; call cleanup(3_INTD); return; endif
        if(.not.NO_CHECK) then
         i=dil_tens_contr_spec_check(cspec)
         if(i.ne.0) then
          if(VERBOSE) write(CONS_OUT,'("#ERROR(tensor_algebra_dil::dil_tensor_contract_pipe)[",i2,"]: ContrSpec error ",i6)')&
          &impir,i
          call cleanup(4_INTD); return
         endif
        endif
        if(darg%store_type.eq.'d'.or.darg%store_type.eq.'D') then
         if(associated(darg%tens_distr_p)) then
          if(darg%tens_distr_p%tsize.le.0) then; call cleanup(5_INTD); return; endif
         else
          call cleanup(6_INTD); return
         endif
        endif
        if(larg%store_type.eq.'d'.or.larg%store_type.eq.'D') then
         if(associated(larg%tens_distr_p)) then
          if(larg%tens_distr_p%tsize.le.0) then; call cleanup(7_INTD); return; endif
         else
          call cleanup(8_INTD); return
         endif
        endif
        if(rarg%store_type.eq.'d'.or.rarg%store_type.eq.'D') then
         if(associated(rarg%tens_distr_p)) then
          if(rarg%tens_distr_p%tsize.le.0) then; call cleanup(9_INTD); return; endif
         else
          call cleanup(10_INTD); return
         endif
        endif
!Get tensor rank/size and other info:
        nd=cspec%ndims_left+cspec%ndims_right  !destination (sub)tensor rank
        nl=cspec%ndims_contr+cspec%ndims_left  !left (sub)tensor rank
        nr=cspec%ndims_contr+cspec%ndims_right !right (sub) tensor rank
        dvol=1_INTL; do i=1,nd; dvol=dvol*cspec%ddims(i); enddo !destination (sub)tensor volume
        lvol=1_INTL; do i=1,nl; lvol=lvol*cspec%ldims(i); enddo !left (sub)tensor volume
        rvol=1_INTL; do i=1,nr; rvol=rvol*cspec%rdims(i); enddo !right (sub)tensor volume
        triv_d=permutation_trivial(nd,cspec%dprmn,i) !triviality of the destination tensor permutation
        triv_l=permutation_trivial(nl,cspec%lprmn,i) !triviality of the left tensor permutation
        triv_r=permutation_trivial(nr,cspec%rprmn,i) !triviality of the right tensor permutation
        call permutation_invert(nd,cspec%dprmn,prmn(1:,0),i) !N2O for the destination tensor
        call permutation_invert(nl,cspec%lprmn,prmn(1:,1),i) !N2O for the left tensor
        call permutation_invert(nr,cspec%rprmn,prmn(1:,2),i) !NwO for the right tensor
        if(nd.lt.0.or.nl.le.0.or.nr.le.0.or.dvol.lt.1_INTL.or.lvol.lt.1_INTL.or.rvol.lt.1_INTL) then
         call cleanup(11_INTD); return
        endif
!Check dimension range consistency:
        if(DIL_DEBUG) then
         select case(darg%store_type)
         case('l','L')
          write(CONS_OUT,'("#DEBUG(DIL): LOC DEST STORED LBND:",16(1x,i6))') darg%tens_loc%base(1:nd)+IND_NUM_START
          write(CONS_OUT,'("#DEBUG(DIL): LOC DEST STORED UBND:",16(1x,i6))') darg%tens_loc%base(1:nd)+darg%tens_loc%dims(1:nd)&
          &+IND_NUM_START-1
         case('d','D')
          if(associated(darg%tens_distr_p)) then
           write(CONS_OUT,'("#DEBUG(DIL): DISTR DEST STORED LBND:",16(1x,i6))') (/(IND_NUM_START,i=1,nd)/)
           write(CONS_OUT,'("#DEBUG(DIL): DISTR DEST STORED UBND:",16(1x,i6))') darg%tens_distr_p%dims(1:nd)+IND_NUM_START-1
          endif
         end select
         write(CONS_OUT,'("#DEBUG(DIL): DEST PROCESSED LBND:",16(1x,i6))') cspec%dbase(1:nd)+IND_NUM_START
         write(CONS_OUT,'("#DEBUG(DIL): DEST PROCESSED UBND:",16(1x,i6))') cspec%dbase(1:nd)+cspec%ddims(1:nd)+IND_NUM_START-1
         select case(larg%store_type)
         case('l','L')
          write(CONS_OUT,'("#DEBUG(DIL): LOC LEFT STORED LBND:",16(1x,i6))') larg%tens_loc%base(1:nl)+IND_NUM_START
          write(CONS_OUT,'("#DEBUG(DIL): LOC LEFT STORED UBND:",16(1x,i6))') larg%tens_loc%base(1:nl)+larg%tens_loc%dims(1:nl)&
          &+IND_NUM_START-1
         case('d','D')
          if(associated(larg%tens_distr_p)) then
           write(CONS_OUT,'("#DEBUG(DIL): DISTR LEFT STORED LBND:",16(1x,i6))') (/(IND_NUM_START,i=1,nl)/)
           write(CONS_OUT,'("#DEBUG(DIL): DISTR LEFT STORED UBND:",16(1x,i6))') larg%tens_distr_p%dims(1:nl)+IND_NUM_START-1
          endif
         end select
         write(CONS_OUT,'("#DEBUG(DIL): LEFT PROCESSED LBND:",16(1x,i6))') cspec%lbase(1:nl)+IND_NUM_START
         write(CONS_OUT,'("#DEBUG(DIL): LEFT PROCESSED UBND:",16(1x,i6))') cspec%lbase(1:nl)+cspec%ldims(1:nl)+IND_NUM_START-1
         select case(rarg%store_type)
         case('l','L')
          write(CONS_OUT,'("#DEBUG(DIL): LOC RIGT STORED LBND:",16(1x,i6))') rarg%tens_loc%base(1:nr)+IND_NUM_START
          write(CONS_OUT,'("#DEBUG(DIL): LOC RIGT STORED UBND:",16(1x,i6))') rarg%tens_loc%base(1:nr)+rarg%tens_loc%dims(1:nr)&
          &+IND_NUM_START-1
         case('d','D')
          if(associated(rarg%tens_distr_p)) then
           write(CONS_OUT,'("#DEBUG(DIL): DISTR RIGT STORED LBND:",16(1x,i6))') (/(IND_NUM_START,i=1,nr)/)
           write(CONS_OUT,'("#DEBUG(DIL): DISTR RIGT STORED UBND:",16(1x,i6))') rarg%tens_distr_p%dims(1:nr)+IND_NUM_START-1
          endif
         end select
         write(CONS_OUT,'("#DEBUG(DIL): RIGT PROCESSED LBND:",16(1x,i6))') cspec%rbase(1:nr)+IND_NUM_START
         write(CONS_OUT,'("#DEBUG(DIL): RIGT PROCESSED UBND:",16(1x,i6))') cspec%rbase(1:nr)+cspec%rdims(1:nr)+IND_NUM_START-1
        endif
 !Destination tensor argument:
        if(nd.gt.0) then
         j=0
         select case(darg%store_type)
         case('l','L')
          do i=1,nd; if(cspec%dbase(i).lt.darg%tens_loc%base(i)) then; j=1; exit; endif; enddo
          if(j.ne.0) then; call cleanup(12_INTD); return; endif
          do i=1,nd; if(cspec%dbase(i)+cspec%ddims(i).gt.darg%tens_loc%base(i)+darg%tens_loc%dims(i)) then; j=2; exit; endif; enddo
          if(j.ne.0) then; call cleanup(13_INTD); return; endif
         case('d','D')
          if(associated(darg%tens_distr_p)) then
           do i=1,nd; if(cspec%dbase(i).lt.0) then; j=1; exit; endif; enddo
           if(j.ne.0) then; call cleanup(14_INTD); return; endif
           do i=1,nd; if(cspec%dbase(i)+cspec%ddims(i)-1.gt.darg%tens_distr_p%dims(i)) then; j=2; exit; endif; enddo
           if(j.ne.0) then; call cleanup(15_INTD); return; endif
          else
           call cleanup(16_INTD); return
          endif
         case default
          call cleanup(17_INTD); return
         end select
        endif
 !Left tensor argument:
        if(nl.gt.0) then
         j=0
         select case(larg%store_type)
         case('l','L')
          do i=1,nl; if(cspec%lbase(i).lt.larg%tens_loc%base(i)) then; j=1; exit; endif; enddo
          if(j.ne.0) then; call cleanup(18_INTD); return; endif
          do i=1,nl; if(cspec%lbase(i)+cspec%ldims(i).gt.larg%tens_loc%base(i)+larg%tens_loc%dims(i)) then; j=2; exit; endif; enddo
          if(j.ne.0) then; call cleanup(19_INTD); return; endif
         case('d','D')
          if(associated(larg%tens_distr_p)) then
           do i=1,nl; if(cspec%lbase(i).lt.0) then; j=1; exit; endif; enddo
           if(j.ne.0) then; call cleanup(20_INTD); return; endif
           do i=1,nl; if(cspec%lbase(i)+cspec%ldims(i)-1.gt.larg%tens_distr_p%dims(i)) then; j=2; exit; endif; enddo
           if(j.ne.0) then; call cleanup(21_INTD); return; endif
          else
           call cleanup(22_INTD); return
          endif
         case default
          call cleanup(23_INTD); return
         end select
        endif
 !Right tensor argument:
        if(nr.gt.0) then
         j=0
         select case(rarg%store_type)
         case('l','L')
          do i=1,nr; if(cspec%rbase(i).lt.rarg%tens_loc%base(i)) then; j=1; exit; endif; enddo
          if(j.ne.0) then; call cleanup(24_INTD); return; endif
          do i=1,nr; if(cspec%rbase(i)+cspec%rdims(i).gt.rarg%tens_loc%base(i)+rarg%tens_loc%dims(i)) then; j=2; exit; endif; enddo
          if(j.ne.0) then; call cleanup(25_INTD); return; endif
         case('d','D')
          if(associated(rarg%tens_distr_p)) then
           do i=1,nr; if(cspec%rbase(i).lt.0) then; j=1; exit; endif; enddo
           if(j.ne.0) then; call cleanup(26_INTD); return; endif
           do i=1,nr; if(cspec%rbase(i)+cspec%rdims(i)-1.gt.rarg%tens_distr_p%dims(i)) then; j=2; exit; endif; enddo
           if(j.ne.0) then; call cleanup(27_INTD); return; endif
          else
           call cleanup(28_INTD); return
          endif
         case default
          call cleanup(29_INTD); return
         end select         
        endif
!Count available computing devices: 
 !Count GPUs:
        if(present(num_gpus)) then
         if(num_gpus.gt.0) then
          if(num_gpus.le.MAX_GPUS) then
           num_dev=num_dev+num_gpus; gpu_on=.true.
          else
           call cleanup(30_INTD); return
          endif
         endif
        endif
 !Count MICs:
        if(present(num_mics)) then
         if(num_mics.gt.0) then
          if(num_mics.le.MAX_MICS) then
           num_dev=num_dev+num_mics; mic_on=.true.
          else
           call cleanup(31_INTD); return
          endif
         endif
        endif
!Allocate/associate buffers:`Accelerators should allocate only an upload buffer (full size) and 3 MPI buffers (tile size)
 !Allocate global Host buffer space:
        tot_buf_vol=(mem_lim-mod(mem_lim,ALIGNMENT*BUFS_PER_DEV*num_dev))/size_of_real !total buffer volume for all devices
        dev_buf_vol=tot_buf_vol/num_dev      !total buffer volume for each device
        arg_buf_vol=dev_buf_vol/BUFS_PER_DEV !max buffer volume for each tensor argument (tensor part)
        if(mod(arg_buf_vol*size_of_real,ALIGNMENT).ne.0_INTL) then; call cleanup(32_INTD); return; endif !trap
 !Associate CPU buffers:
        i0=0_INTL; k=dil_dev_num(DEV_HOST_CPU,0_INTD)
        do i=1,buf(k)%num_bufs
         buf(k)%arg_buf(i)%buf_ptr=>hbuf(i0+1_INTL:i0+arg_buf_vol)
         buf(k)%arg_buf(i)%buf_vol=arg_buf_vol
         i0=i0+arg_buf_vol
        enddo
 !Associate GPU buffers:`Only an upload buffer and 3 small MPI buffers on Host (+6 buffers on device)
        if(gpu_on) then
         do j=0,num_gpus-1
          k=dil_dev_num(DEV_NVIDIA_GPU,j)
          do i=1,buf(k)%num_bufs
           buf(k)%arg_buf(i)%buf_ptr=>hbuf(i0+1_INTL:i0+arg_buf_vol)
           buf(k)%arg_buf(i)%buf_vol=arg_buf_vol
           i0=i0+arg_buf_vol
          enddo
         enddo
        endif
 !Associate MIC buffers:`Only an upload buffer and 3 small MPI buffers on Host (+6 buffers on device)
        if(mic_on) then
         do j=0,num_mics-1
          k=dil_dev_num(DEV_INTEL_MIC,j)
          do i=1,buf(k)%num_bufs
           buf(k)%arg_buf(i)%buf_ptr=>hbuf(i0+1_INTL:i0+arg_buf_vol)
           buf(k)%arg_buf(i)%buf_vol=arg_buf_vol
           i0=i0+arg_buf_vol
          enddo
         enddo
        endif
        if(i0.gt.tot_buf_vol) then; call cleanup(33_INTD); return; endif
        if(DIL_DEBUG) write(CONS_OUT,'("#DEBUG(tensor_algebra_dil::dil_tensor_contract_pipe): Buffer volumes (W):",3(1x,i12))')&
        &tot_buf_vol,dev_buf_vol,arg_buf_vol !debug
!Partition tensor contraction into smaller parts (tasks) if argument(s) do not fit into local buffers:
        select case(contr_case)
        case('lll','LLL')
         call dil_tens_contr_partition(cspec,arg_buf_vol,contr_case,task_list,i)
        case('lld','LLD')
         call dil_tens_contr_partition(cspec,arg_buf_vol,contr_case,task_list,i,rarr=rarg%tens_distr_p)
        case('dll','DLL')
         call dil_tens_contr_partition(cspec,arg_buf_vol,contr_case,task_list,i,darr=darg%tens_distr_p)
        case('ldd','LDD')
         call dil_tens_contr_partition(cspec,arg_buf_vol,contr_case,task_list,i,larr=larg%tens_distr_p,rarr=rarg%tens_distr_p)
        case('dld','DLD')
         call dil_tens_contr_partition(cspec,arg_buf_vol,contr_case,task_list,i,darr=darg%tens_distr_p,rarr=rarg%tens_distr_p)
        case('ddd','DDD')
         call dil_tens_contr_partition(cspec,arg_buf_vol,contr_case,task_list,i,&
                                      &darr=darg%tens_distr_p,larr=larg%tens_distr_p,rarr=rarg%tens_distr_p)
        case default
         call cleanup(34_INTD); return
        end select
        if(i.ne.0) then
         if(VERBOSE) write(CONS_OUT,'("#ERROR(tensor_algebra_dil::dil_tensor_contract_pipe): Partitioning failed: ",i9)') i
         call cleanup(35_INTD); return
        endif
        if(DIL_DEBUG) write(CONS_OUT,'("#DEBUG(tensor_algebra_dil::dil_tensor_contract_pipe)[",i2,"]: Case: ",'//&
        &'A3,": Task amount: ",i9)') impir,contr_case,task_list%num_tasks !debug
!Reshuffle the task list:
        if(TASK_RESHUFFLE) then
         j=impir
         call dil_contr_task_list_shuffle(task_list,i,shift=j)
         if(i.ne.0) then
          if(VERBOSE) write(CONS_OUT,'("#ERROR(tensor_algebra_dil::dil_tensor_contract_pipe): Task reshuffling failed: ",i9)') i
          call cleanup(36_INTD); return
         endif
         if(DIL_DEBUG) write(CONS_OUT,'("DEBUG(tensor_algebra_dil::dil_tensor_contract_pipe)[",i2'//&
         &',"]: Task list reshuffled: Shift = ",i6)') impir,j
        endif
!Perform partitioned tensor contraction via pipelining:
        if(task_list%num_tasks.gt.0) then
         if(.not.associated(task_list%contr_tasks)) then; call cleanup(37_INTD); return; endif !trap
         if(lbound(task_list%contr_tasks,1).eq.1) then !first available task is #1
          first_avail(:)=1_INTD
         else
          call cleanup(38_INTD); return
         endif !trap
 !CPU (multicore): !`Implement REUSE for tensor arguments
         dev=dil_dev_num(DEV_HOST_CPU,0_INTD); buf_conf(1:BUFS_PER_DEV,dev)=(/(i,i=1,BUFS_PER_DEV)/)
         tasks_done(dev)=0; task_curr(dev)=0; task_prev(dev)=0; task_next(dev)=0
         task_curr(dev)=dil_get_next_task(DEV_HOST_CPU,0_INTD); first_task=.true.
         next_load=.false.; prev_store=.false.; args_here=.false.; err_curr=.false.; err_prev=.false.; err_next=.false.
         tloop: do while(task_prev(dev).ge.0_INTD) !loop over the tensor contraction tasks
          task_next(dev)=dil_get_next_task(DEV_HOST_CPU,0_INTD)
          if(DIL_DEBUG) then !debug begin
           write(CONS_OUT,'("#DEBUG(tensor_algebra_dil::dil_tensor_contract_pipe): DEVICE: ",i2,": Tasks(c,p,n):",3(1x,l1,1x,i7))')&
           &dev,(.not.err_curr),task_curr(dev),(.not.err_prev),task_prev(dev),(.not.err_next),task_next(dev) !debug
           if(task_curr(dev).ge.1.and.task_curr(dev).le.task_list%num_tasks)&
           &call dil_contr_task_print(cspec,task_list%contr_tasks(task_curr(dev))) !debug
          endif !debug end
          arg_reuse=dil_mark_arg_reuse(task_curr(dev),task_prev(dev)) !determine arguments to be reused from the previous task
          arg_keep=dil_mark_arg_reuse(task_curr(dev),task_next(dev))  !determine arguments to be reused in the next task
          if(first_task) then
           if(DIL_DEBUG) write(CONS_OUT,'(1x,"#DEBUG(DIL): First task arg load initiation started ... ")')
           tms=thread_wtime()
           call dil_args_load_start(task_list%contr_tasks(task_curr(dev)),i) !non-blocking
           if(i.ne.0) then; task_list%contr_tasks(task_curr(dev))%task_stat=TASK_ERR_LDS; err_curr=.true.; endif
           next_load=.true.; first_task=.false.
           tm=thread_wtime(tms)
           if(DIL_DEBUG) write(CONS_OUT,'(1x,"Done (",l1,"): ",F10.4," s")') (.not.err_curr),tm
          endif
          if(next_load) then
           if(DIL_DEBUG) write(CONS_OUT,'(1x,"#DEBUG(DIL): Current task arg load completion started ... ")')
           tms=thread_wtime()
           if(.not.err_curr) then
            call dil_args_load_complete(task_list%contr_tasks(task_curr(dev)),i) !blocking
            if(i.ne.0) then; task_list%contr_tasks(task_curr(dev))%task_stat=TASK_ERR_LDC; err_curr=.true.; endif
           endif
           next_load=.false.; args_here=.true.
           tm=thread_wtime(tms)
           if(DIL_DEBUG) write(CONS_OUT,'(1x,"Done (",l1,"): ",F10.4," s")') (.not.err_curr),tm
          endif
          if(args_here) then
           if(DIL_DEBUG) write(CONS_OUT,'(1x,"#DEBUG(DIL): Current task input prepare started ... ")')
           tms=thread_wtime()
           if(.not.err_curr) then
            call dil_args_prepare_input(task_list%contr_tasks(task_curr(dev)),i) !blocking
            if(i.ne.0) then; task_list%contr_tasks(task_curr(dev))%task_stat=TASK_ERR_PRI; err_curr=.true.; endif
           endif
           tm=thread_wtime(tms)
           if(DIL_DEBUG) write(CONS_OUT,'(1x,"Done (",l1,"): ",F10.4," s")') (.not.err_curr),tm
          endif
          if(task_next(dev).gt.0) then
           if(DIL_DEBUG) write(CONS_OUT,'(1x,"#DEBUG(DIL): Next task arg load initiation started ... ")')
           tms=thread_wtime()
           call dil_args_load_start(task_list%contr_tasks(task_next(dev)),i) !non-blocking
           if(i.ne.0) then; task_list%contr_tasks(task_next(dev))%task_stat=TASK_ERR_LDS; err_next=.true.; endif
           next_load=.true.
           tm=thread_wtime(tms)
           if(DIL_DEBUG) write(CONS_OUT,'(1x,"Done (",l1,"): ",F10.4," s")') (.not.err_next),tm
          endif
          if(args_here) then
           if(DIL_DEBUG) write(CONS_OUT,'(1x,"#DEBUG(DIL): Current task computation started ... ")')
           tms=thread_wtime()
           if(.not.err_curr) then
            call dil_mm_compute(task_list%contr_tasks(task_curr(dev)),i) !blocking
            if(i.ne.0) then; task_list%contr_tasks(task_curr(dev))%task_stat=TASK_ERR_CMT; err_curr=.true.; endif
           endif
           tm=thread_wtime(tms); tmm=tmm+tm
           if(DIL_DEBUG) write(CONS_OUT,'(1x,"Done (",l1,"): ",F10.4," s")') (.not.err_curr),tm
          endif
          if(prev_store) then
           if(DIL_DEBUG) write(CONS_OUT,'(1x,"#DEBUG(DIL): Previous task arg store completion started ... ")')
           tms=thread_wtime()
           if(.not.err_prev) then
            call dil_args_store_complete(task_list%contr_tasks(task_prev(dev)),i) !blocking
            if(i.eq.0) then
             task_list%contr_tasks(task_prev(dev))%task_stat=TASK_COMPLETED
             tasks_done(dev)=tasks_done(dev)+1
            else
             task_list%contr_tasks(task_prev(dev))%task_stat=TASK_ERR_STC; err_prev=.true.
            endif
           endif
           prev_store=.false.
           tm=thread_wtime(tms)
           if(DIL_DEBUG) write(CONS_OUT,'(1x,"Done (",l1,"): ",F10.4," s")') (.not.err_prev),tm
          endif
          if(args_here) then
           if(.not.err_curr) then
            if(DIL_DEBUG) write(CONS_OUT,'(1x,"#DEBUG(DIL): Current task output prepare started ... ")')
            tms=thread_wtime()
            call dil_args_prepare_output(task_list%contr_tasks(task_curr(dev)),i) !blocking
            tm=thread_wtime(tms)
            if(i.eq.0) then
             if(DIL_DEBUG) write(CONS_OUT,'(1x,"Done (",l1,"): ",F10.4," s")') (.not.err_curr),tm
             if(DIL_DEBUG) write(CONS_OUT,'(1x,"#DEBUG(DIL): Current task arg store initiation started ... ")')
             tms=thread_wtime()
             call dil_args_store_start(task_list%contr_tasks(task_curr(dev)),i) !non-blocking
             if(i.ne.0) then; task_list%contr_tasks(task_curr(dev))%task_stat=TASK_ERR_STS; err_curr=.true.; endif
             tm=thread_wtime(tms)
             if(DIL_DEBUG) write(CONS_OUT,'(1x,"Done (",l1,"): ",F10.4," s")') (.not.err_curr),tm
            else
             task_list%contr_tasks(task_curr(dev))%task_stat=TASK_ERR_PRO; err_curr=.true.
             if(DIL_DEBUG) write(CONS_OUT,'(1x,"Done (",l1,"): ",F10.4," s")') (.not.err_curr),tm
            endif
           endif
           prev_store=.true.; args_here=.false.
          endif
          if(task_curr(dev).ge.1.and.task_curr(dev).le.task_list%num_tasks) then !count Flops
           tc_flops=tc_flops+dble(task_list%contr_tasks(task_curr(dev))%flops)
          endif
          task_prev(dev)=task_curr(dev); task_curr(dev)=task_next(dev)
          err_prev=err_curr; err_curr=err_next; err_next=.false.
         enddo tloop
        else
         if(DIL_DEBUG) write(CONS_OUT,'("#DEBUG(tensor_algebra_dil::dil_tensor_contract_pipe)[",i2,"]: No tasks!")') impir !debug
        endif
        tm=thread_wtime(tmb)
        if(DIL_DEBUG) write(CONS_OUT,'("#DEBUG(tensor_algebra_dil::dil_tensor_contract_pipe)[",i2,"]: Done in ",F10.4'&
        &//'," s ( ",F15.4," GFlops/s VS MM ",F15.4," GFlops/s). Ok")') impir,tm,tc_flops/(tm*1024d0*1024d0*1024d0),&
        &mm_flops/(tmm*1024d0*1024d0*1024d0) !debug
        call cleanup(0_INTD)
        return

        contains

         subroutine dil_args_load_start(tsk,errc) !Non-blocking
 !This subroutine starts load of input tensor arguments for a task <tsk>:
 ! * Distributed tensor: starts loading necessary tiles into a local buffer;
 ! * Local tensor: does nothing.
         type(contr_task_t), intent(in):: tsk !in: tensor contraction task
         integer(INTD), intent(out):: errc    !out: error code (0:success)
         integer(INTD):: jdev,je

         errc=0
         jdev=dil_dev_num(tsk%dev_kind,tsk%dev_id)
!         if((darg%store_type.eq.'d'.or.darg%store_type.eq.'D').and.(.not.cspec%dest_zero)) then !distributed tensors only
!          if(arg_reuse(1:1).ne.'R') then
!           call dil_tensor_prefetch_start(darg%tens_distr_p,tsk%dest_arg,buf(jdev)%arg_buf(buf_conf(1,jdev)),je,win_lck)
!           if(je.ne.0) then; errc=1; return; endif
!          endif
!         endif
         if(larg%store_type.eq.'d'.or.larg%store_type.eq.'D') then !distributed tensors only
          if(arg_reuse(2:2).ne.'R') then
           call dil_tensor_prefetch_start(larg%tens_distr_p,tsk%left_arg,buf(jdev)%arg_buf(buf_conf(2,jdev)),je,win_lck)
           if(je.ne.0) then; errc=2; return; endif
          endif
         endif
         if(rarg%store_type.eq.'d'.or.rarg%store_type.eq.'D') then !distributed tensors only
          if(arg_reuse(3:3).ne.'R') then
           call dil_tensor_prefetch_start(rarg%tens_distr_p,tsk%right_arg,buf(jdev)%arg_buf(buf_conf(3,jdev)),je,win_lck)
           if(je.ne.0) then; errc=3; return; endif
          endif
         endif
         return
         end subroutine dil_args_load_start

         subroutine dil_args_load_complete(tsk,errc) !Blocking
 !This subroutine completes load of input tensor arguments for a task <tsk>:
 ! * Distributed tensor: completes loading necessary tiles into a local buffer;
 ! * Local tensor: does nothing.
         type(contr_task_t), intent(in):: tsk !in: tensor contraction task
         integer(INTD), intent(out):: errc    !out: error code (0:success)
         integer(INTD):: jdev,je

         errc=0
         jdev=dil_dev_num(tsk%dev_kind,tsk%dev_id)
!         if((darg%store_type.eq.'d'.or.darg%store_type.eq.'D').and.(.not.cspec%dest_zero)) then !distributed tensors only
!          if(arg_reuse(1:1).ne.'R') then
!           call dil_tensor_prefetch_complete(darg%tens_distr_p,tsk%dest_arg,buf(jdev)%arg_buf(buf_conf(1,jdev)),je,win_lck)
!           if(je.ne.0) then; errc=1; return; endif
!          endif
!         endif
         if(larg%store_type.eq.'d'.or.larg%store_type.eq.'D') then !distributed tensors only
          if(arg_reuse(2:2).ne.'R') then
           call dil_tensor_prefetch_complete(larg%tens_distr_p,tsk%left_arg,buf(jdev)%arg_buf(buf_conf(2,jdev)),je,win_lck)
           if(je.ne.0) then; errc=2; return; endif
          endif
         endif
         if(rarg%store_type.eq.'d'.or.rarg%store_type.eq.'D') then !distributed tensors only
          if(arg_reuse(3:3).ne.'R') then
           call dil_tensor_prefetch_complete(rarg%tens_distr_p,tsk%right_arg,buf(jdev)%arg_buf(buf_conf(3,jdev)),je,win_lck)
           if(je.ne.0) then; errc=3; return; endif
          endif
         endif
         return
         end subroutine dil_args_load_complete

         subroutine dil_args_store_start(tsk,errc) !Non-blocking
 !This subroutine starts uploading the output tensor argument for a task <tsk>:
 ! * Distributed tensor: starts uploading tiles from a local buffer;
 ! * Local tensor: does nothing.
         type(contr_task_t), intent(in):: tsk !in: tensor contraction task
         integer(INTD), intent(out):: errc    !out: error code (0:success)
         integer(INTD):: jdev,je

         errc=0
         jdev=dil_dev_num(tsk%dev_kind,tsk%dev_id)
         if(darg%store_type.eq.'d'.or.darg%store_type.eq.'D') then !distributed tensors only
          if(arg_keep(1:1).ne.'R') then
           call dil_tensor_upload_start(darg%tens_distr_p,tsk%dest_arg,buf(jdev)%arg_buf(buf_conf(7,jdev)),je,win_lck)
           if(je.ne.0) then; errc=1; return; endif
          endif
         endif
         return
         end subroutine dil_args_store_start

         subroutine dil_args_store_complete(tsk,errc) !Blocking
 !This subroutine completes uploading the output tensor argument for a task <tsk>:
 ! * Distributed tensor: completes uploading tiles from a local buffer;
 ! * Local tensor: does nothing.
         type(contr_task_t), intent(in):: tsk !in: tensor contraction task
         integer(INTD), intent(out):: errc    !out: error code (0:success)
         integer(INTD):: jdev,je

         errc=0
         jdev=dil_dev_num(tsk%dev_kind,tsk%dev_id)
         if(darg%store_type.eq.'d'.or.darg%store_type.eq.'D') then !distributed tensors only
          if(task_curr(jdev).ge.0.or.(.not.async)) then !not the last upload for this device (finalize it)
           if(arg_keep(1:1).ne.'R') then
            call dil_tensor_upload_complete(darg%tens_distr_p,tsk%dest_arg,buf(jdev)%arg_buf(buf_conf(7,jdev)),je,win_lck)
           endif
          else !last upload for this device (asynchronous)
           call dil_tensor_upload_complete(darg%tens_distr_p,tsk%dest_arg,buf(jdev)%arg_buf(buf_conf(7,jdev)),je,win_lck,&
                                          &nasync,lasync)
          endif
          if(je.ne.0) then; errc=1; return; endif
         endif
         return
         end subroutine dil_args_store_complete

         subroutine dil_args_prepare_input(tsk,errc) !Blocking
 !This subroutine prepares input tensor arguments for a tensor contraction task <tsk>:
 ! * Distributed tensor: unpack tiles (if multiple), permute dimensions (if needed);
 ! * Local tensor: permute dimensions (if needed).
         type(contr_task_t), intent(in):: tsk !in: tensor contraction task
         integer(INTD), intent(out):: errc    !out: error code (0:success)
         integer(INTD):: jdev,jb,jf,je
         integer(INTL):: jvol

         errc=0
         jdev=dil_dev_num(tsk%dev_kind,tsk%dev_id)
  !Destination tensor:
         if(arg_reuse(1:1).ne.'R') then
          jb=buf_conf(1,jdev); jf=buf_conf(4,jdev); jvol=dil_subtensor_vol(tsk%dest_arg)
          if(darg%store_type.eq.'l'.or.darg%store_type.eq.'L') then !local tensor
           if(nd.gt.0) then !slice
            call dil_tensor_slice(nd,darg%tens_loc%elems,darg%tens_loc%dims,buf(jdev)%arg_buf(jb)%buf_ptr,tsk%dest_arg%dims,&
                                 &tsk%dest_arg%lbnd(1:nd)-(darg%tens_loc%base(1:nd)+IND_NUM_START),je)
            if(je.ne.0) then; errc=1; return; endif
           elseif(nd.eq.0) then !scalar
            buf(jdev)%arg_buf(jb)%buf_ptr(1)=darg%tens_loc%elems(1)
           endif
           if(.not.triv_d) then !permute (if needed)
            call dil_tensor_transpose(nd,tsk%dest_arg%dims,cspec%dprmn,&
                                     &buf(jdev)%arg_buf(jb)%buf_ptr,buf(jdev)%arg_buf(jf)%buf_ptr,je)
            if(je.ne.0) then; errc=2; return; endif
            je=jb; jb=jf; jf=je
           endif
          elseif(darg%store_type.eq.'d'.or.darg%store_type.eq.'D') then !distributed tensor
!          call dil_arg_buf_clean(buf(jdev)%arg_buf(jb),je,jvol) !zero out buffer
!          if(je.ne.0) then; errc=3; return; endif
!          if(nd.gt.0) then
!           call dil_tens_unpack_from_tiles(darg%tens_distr_p,tsk%dest_arg,buf(jdev)%arg_buf(jb),buf(jdev)%arg_buf(jf),je) !unpack tile(s)
!           if(je.ne.0) then; errc=4; return; endif
!           je=jb; jb=jf; jf=je
!          endif
          else
           errc=5; return
          endif
          buf_conf(4,jdev)=jb; buf_conf(1,jdev)=jf
         endif
  !Left tensor:
         if(arg_reuse(2:2).ne.'R') then
          jb=buf_conf(2,jdev); jf=buf_conf(5,jdev); jvol=dil_subtensor_vol(tsk%left_arg)
          if(larg%store_type.eq.'d'.or.larg%store_type.eq.'D') then !distributed tensor
           if(nl.gt.0) then
            if(.not.one_tile_only(larg%tens_distr_p,tsk%left_arg)) then !more than one tile
             call dil_tens_unpack_from_tiles(larg%tens_distr_p,tsk%left_arg,buf(jdev)%arg_buf(jb),buf(jdev)%arg_buf(jf),je) !unpack tiles
             if(je.ne.0) then; errc=6; return; endif
             je=jb; jb=jf; jf=je
            endif
           else
            errc=7; return
           endif
          elseif(larg%store_type.eq.'l'.or.larg%store_type.eq.'L') then !local tensor
           if(nl.gt.0) then !slice
            call dil_tensor_slice(nl,larg%tens_loc%elems,larg%tens_loc%dims,buf(jdev)%arg_buf(jb)%buf_ptr,tsk%left_arg%dims,&
                                 &tsk%left_arg%lbnd(1:nl)-(larg%tens_loc%base(1:nl)+IND_NUM_START),je)
            if(je.ne.0) then; errc=8; return; endif
           elseif(nl.eq.0) then !scalar
            buf(jdev)%arg_buf(jb)%buf_ptr(1)=larg%tens_loc%elems(1)
           endif
          else
           errc=9; return
          endif
          if(.not.triv_l) then !permute (if needed)
           call dil_tensor_transpose(nl,tsk%left_arg%dims,cspec%lprmn,&
                                    &buf(jdev)%arg_buf(jb)%buf_ptr,buf(jdev)%arg_buf(jf)%buf_ptr,je)
           if(je.ne.0) then; errc=10; return; endif
           je=jb; jb=jf; jf=je
          endif
          buf_conf(5,jdev)=jb; buf_conf(2,jdev)=jf
         endif
  !Right tensor:
         if(arg_reuse(3:3).ne.'R') then
          jb=buf_conf(3,jdev); jf=buf_conf(6,jdev); jvol=dil_subtensor_vol(tsk%right_arg)
          if(rarg%store_type.eq.'d'.or.rarg%store_type.eq.'D') then !distributed tensor
           if(nr.gt.0) then
            if(.not.one_tile_only(rarg%tens_distr_p,tsk%right_arg)) then !more than one tile
             call dil_tens_unpack_from_tiles(rarg%tens_distr_p,tsk%right_arg,buf(jdev)%arg_buf(jb),buf(jdev)%arg_buf(jf),je) !unpack tiles
             if(je.ne.0) then; errc=11; return; endif
             je=jb; jb=jf; jf=je
            endif
           else
            errc=12; return
           endif
          elseif(rarg%store_type.eq.'l'.or.rarg%store_type.eq.'L') then !local tensor
           if(nr.gt.0) then !slice
            call dil_tensor_slice(nr,rarg%tens_loc%elems,rarg%tens_loc%dims,buf(jdev)%arg_buf(jb)%buf_ptr,tsk%right_arg%dims,&
                                 &tsk%right_arg%lbnd(1:nr)-(rarg%tens_loc%base(1:nr)+IND_NUM_START),je)
            if(je.ne.0) then; errc=13; return; endif
           elseif(nr.eq.0) then !scalar
            buf(jdev)%arg_buf(jb)%buf_ptr(1)=rarg%tens_loc%elems(1)
           endif
          else
           errc=14; return
          endif
          if(.not.triv_r) then !permute (if needed)
           call dil_tensor_transpose(nr,tsk%right_arg%dims,cspec%rprmn,&
                                    &buf(jdev)%arg_buf(jb)%buf_ptr,buf(jdev)%arg_buf(jf)%buf_ptr,je)
           if(je.ne.0) then; errc=15; return; endif
           je=jb; jb=jf; jf=je
          endif
          buf_conf(6,jdev)=jb; buf_conf(3,jdev)=jf
         endif
         return
         end subroutine dil_args_prepare_input

         subroutine dil_args_prepare_output(tsk,errc) !Blocking
 !This subroutine prepares the output tensor argument for a tensor contraction task <tsk>:
 ! * Distributed tensor: permute dimensions (if needed), pack tiles (if multiple);
 ! * Local tensor: permute dimensions (if needed).
         type(contr_task_t), intent(in):: tsk !in: tensor contraction task
         integer(INTD), intent(out):: errc    !out: error code (0:success)
         integer(INTD):: jdev,jb,jf,je
         integer(INTL):: jvol
         real(8):: jtb,jtm

         errc=0
         jdev=dil_dev_num(tsk%dev_kind,tsk%dev_id)
  !Destination tensor:
         if(arg_keep(1:1).ne.'R') then
          jb=buf_conf(4,jdev); jf=buf_conf(7,jdev); jvol=dil_subtensor_vol(tsk%dest_arg)
          if(.not.triv_d) then !permute (if needed)
           call dil_tensor_transpose(nd,tsk%dest_arg%dims(prmn(1:nd,0)),prmn(1:,0),&
                                    &buf(jdev)%arg_buf(jb)%buf_ptr,buf(jdev)%arg_buf(jf)%buf_ptr,je)
           if(je.ne.0) then; errc=1; return; endif
           je=jb; jb=jf; jf=je
          endif
          if(darg%store_type.eq.'d'.or.darg%store_type.eq.'D') then !distributed tensor
           if(nd.gt.0) then
            if(.not.one_tile_only(darg%tens_distr_p,tsk%dest_arg)) then !more than one tile
             call dil_tens_pack_into_tiles(darg%tens_distr_p,tsk%dest_arg,buf(jdev)%arg_buf(jb),buf(jdev)%arg_buf(jf),je) !pack tiles
             if(je.ne.0) then; errc=2; return; endif
             je=jb; jb=jf; jf=je
            endif
           else
            errc=3; return
           endif
          elseif(darg%store_type.eq.'l'.or.darg%store_type.eq.'L') then !local tensor
           if(DIL_DEBUG) then
            jtb=thread_wtime()
            write(CONS_OUT,'(2x,"#DEBUG(DIL): Updating local destination ...")',ADVANCE='NO')
           endif
           je=0
           if(nd.gt.0) then !insert tensor slice
            call dil_tensor_insert(nd,darg%tens_loc%elems,darg%tens_loc%dims,buf(jdev)%arg_buf(jb)%buf_ptr,tsk%dest_arg%dims,&
                                  &tsk%dest_arg%lbnd(1:nd)-(darg%tens_loc%base(1:nd)+IND_NUM_START),je)
           elseif(nd.eq.0) then
            darg%tens_loc%elems(1)=buf(jdev)%arg_buf(jb)%buf_ptr(1)
           endif
           if(DIL_DEBUG) then
            jtm=thread_wtime(jtb)
            write(CONS_OUT,'(" Done: ",F10.4," s: ",F10.4," GB/s: Status ",i9)')&
            &jtm,dble(2_INTL*jvol*tensor_dp)/(jtm*1024d0*1024d0*1024d0),je
           endif
           if(je.ne.0) then; errc=4; return; endif
          else
           errc=5; return
          endif
          buf_conf(7,jdev)=jb; buf_conf(4,jdev)=jf
         endif
         return
         end subroutine dil_args_prepare_output

         subroutine dil_mm_compute(tsk,errc) !Blocking
 !This subroutine performs a matrix-matrix multiplication for a tensor contraction task <tsk>.
         type(contr_task_t), intent(in):: tsk !in: tensor contraction task
         integer(INTD), intent(out):: errc    !out: error code (0:success)
         integer(INTL):: lld,lrd,lcd
         integer(INTD):: jdev,j0,jp
         real(8):: jtb,jtm,bts

         errc=0
         jdev=dil_dev_num(tsk%dev_kind,tsk%dev_id)
         lld=1_INTL; lrd=1_INTL; lcd=1_INTL
         do j0=1_INTD,cspec%ndims_left !combined extent of the left uncontracted dimensions for this task
          jp=prmn(j0,0) !original position
          lld=lld*tsk%dest_arg%dims(jp)
         enddo
         do j0=cspec%ndims_left+1_INTD,nd !combined extent of the right uncontracted dimensions for this task
          jp=prmn(j0,0) !original position
          lrd=lrd*tsk%dest_arg%dims(jp)
         enddo
         do j0=1_INTD,cspec%ndims_contr !combined extent of the contracted dimensions for this task
          jp=prmn(j0,1) !original position
          lcd=lcd*tsk%left_arg%dims(jp)
         enddo
         if(DIL_DEBUG) write(CONS_OUT,'(2x,"#DEBUG(DIL)[",i2,"]: Matrix dims:",3(1x,i11))',ADVANCE='NO') impir,lld,lrd,lcd
         mm_flops=mm_flops+dble(lld)*dble(lrd)*dble(lcd) !count Matrix Multiplication Flops
         if(lld*lrd.gt.buf(jdev)%arg_buf(buf_conf(4,jdev))%buf_vol) then; errc=1; return; endif !trap
         if(lcd*lld.gt.buf(jdev)%arg_buf(buf_conf(5,jdev))%buf_vol) then; errc=2; return; endif !trap
         if(lcd*lrd.gt.buf(jdev)%arg_buf(buf_conf(6,jdev))%buf_vol) then; errc=3; return; endif !trap
         if(arg_reuse(1:1).eq.'R') then; bts=1d0; else; bts=beta; endif !switch to accumulation if reusing
         if(DIL_DEBUG) jtb=thread_wtime()
         if(tensor_dp.eq.8) then
          call dgemm('T','N',int(lld,BLAS_INT),int(lrd,BLAS_INT),int(lcd,BLAS_INT),alpha,&
                    &buf(jdev)%arg_buf(buf_conf(5,jdev))%buf_ptr,int(lcd,BLAS_INT),&
                    &buf(jdev)%arg_buf(buf_conf(6,jdev))%buf_ptr,int(lcd,BLAS_INT),bts,&
                    &buf(jdev)%arg_buf(buf_conf(4,jdev))%buf_ptr,int(lld,BLAS_INT))
         elseif(tensor_dp.eq.4) then
          call sgemm('T','N',int(lld,BLAS_INT),int(lrd,BLAS_INT),int(lcd,BLAS_INT),alpha,&
                    &buf(jdev)%arg_buf(buf_conf(5,jdev))%buf_ptr,int(lcd,BLAS_INT),&
                    &buf(jdev)%arg_buf(buf_conf(6,jdev))%buf_ptr,int(lcd,BLAS_INT),bts,&
                    &buf(jdev)%arg_buf(buf_conf(4,jdev))%buf_ptr,int(lld,BLAS_INT))
         else
          errc=4
         endif
         if(DIL_DEBUG) then
          jtm=thread_wtime(jtb)
          write(CONS_OUT,'(": ",F10.4," GFlops/s")') (dble(lld)*dble(lrd)*dble(lcd))/(jtm*1024d0*1024d0*1024d0)
         endif
         return
         end subroutine dil_mm_compute

         function one_tile_only(tdistr,tslice) result(oto)
         type(tensor), intent(in):: tdistr    !in: distributed tensor
         type(subtens_t), intent(in):: tslice !in: tensor slice specification
         logical:: oto                        !out: result (?the slice consists only of one tile?)
         integer(INTD):: j0
         oto=.true.
         do j0=1,tdistr%mode
          if(tslice%dims(j0).gt.tdistr%tdim(j0)) then; oto=.false.; exit; endif
         enddo
         return
         end function one_tile_only

         function dil_mark_arg_reuse(tsc,tsn) result(reuse)
 !This function determines which arguments (tensor parts) will be reused among two tasks.
         integer(INTD), intent(in):: tsc !in: first task number
         integer(INTD), intent(in):: tsn !in: second task number
         character(3):: reuse            !out: triplet of letters XXX: {X="N":new arg | X="R":arg reuse}: DLR
         character(1):: jch
         type(subtens_t), pointer:: jstc,jstn
         integer(INTD):: j0
         reuse='NNN' !All three tensor arguments are assumed new by default
         if(ARGS_REUSE) then
          if(tsc.gt.0.and.tsn.gt.0) then
           jstc=>task_list%contr_tasks(tsc)%dest_arg; jstn=>task_list%contr_tasks(tsn)%dest_arg
           if(jstc%rank.eq.jstn%rank) then
            jch='R'
            do j0=1,jstc%rank
             if(jstc%lbnd(j0).ne.jstn%lbnd(j0).or.jstc%dims(j0).ne.jstn%dims(j0)) then; jch='N'; exit; endif
            enddo
            reuse(1:1)=jch
           endif
           jstc=>task_list%contr_tasks(tsc)%left_arg; jstn=>task_list%contr_tasks(tsn)%left_arg
           if(jstc%rank.eq.jstn%rank) then
            jch='R'
            do j0=1,jstc%rank
             if(jstc%lbnd(j0).ne.jstn%lbnd(j0).or.jstc%dims(j0).ne.jstn%dims(j0)) then; jch='N'; exit; endif
            enddo
            reuse(2:2)=jch
           endif
           jstc=>task_list%contr_tasks(tsc)%right_arg; jstn=>task_list%contr_tasks(tsn)%right_arg
           if(jstc%rank.eq.jstn%rank) then
            jch='R'
            do j0=1,jstc%rank
             if(jstc%lbnd(j0).ne.jstn%lbnd(j0).or.jstc%dims(j0).ne.jstn%dims(j0)) then; jch='N'; exit; endif
            enddo
            reuse(3:3)=jch
           endif
          endif
         endif
         return
         end function dil_mark_arg_reuse

         integer(INTD) function dil_get_next_task(dvk,dvn) !MASTER THREAD only!
 !This function selects the next tensor contraction task from the global task list for execution on a specific device.
 !Negative value on return means the tasks are over.
         integer(INTD), intent(in):: dvk !in: device kind
         integer(INTD), intent(in):: dvn !in: device id (within its kind)
         integer(INTD):: jn,jlo,jf,jdev
         dil_get_next_task=-1_INTD; jdev=dil_dev_num(dvk,dvn)
  !Find an appropriate task for this device:
         if(first_avail(jdev).gt.0_INTD) then
   !Forward loop:
          do jn=first_avail(jdev),task_list%num_tasks
           if(task_list%reordered) then; jlo=task_list%task_order(jn); else; jlo=jn; endif
           if(task_list%contr_tasks(jlo)%task_stat.eq.TASK_SET) then
            dil_get_next_task=jlo; jf=jn; exit
           endif
          enddo
   !Backward loop:
          if(dil_get_next_task.le.0_INTD) then
           do jn=first_avail(jdev)-1_INTD,1_INTD,-1_INTD
            if(task_list%reordered) then; jlo=task_list%task_order(jn); else; jlo=jn; endif
            if(task_list%contr_tasks(jlo)%task_stat.eq.TASK_SET) then
             dil_get_next_task=jlo; jf=jn; exit
            endif
           enddo
          endif
    !Set the task if found:
          if(dil_get_next_task.gt.0_INTD) then
           task_list%contr_tasks(dil_get_next_task)%dev_kind=dvk
           task_list%contr_tasks(dil_get_next_task)%dev_id=dvn
           task_list%contr_tasks(dil_get_next_task)%task_stat=TASK_SCHEDULED
           first_avail(jdev)=jf+1_INTD
           if(first_avail(jdev).gt.task_list%num_tasks) first_avail(jdev)=-1_INTD
          else
           first_avail(jdev)=-1_INTD
          endif
         endif
         return
         end function dil_get_next_task

         subroutine cleanup(errc) !MASTER THREAD only!
 !This function cleans up the local data used in <dil_tensor_contract_pipe>.
         integer(INTD), intent(in):: errc !in: error code
         integer(INTD):: j0,je
         ierr=errc
         if(errc.ne.0.and.VERBOSE) write(CONS_OUT,'("#ERROR(tensor_algebra_dil::dil_tensor_contract_pipe)[",i2,"]: Error ",i9)')&
         &impir,errc
 !CPU buffers:
         call dil_dev_buf_destroy(buf(dil_dev_num(DEV_HOST_CPU,0_INTD)),je)
 !GPU buffers:
         if(gpu_on) then
          do j0=0,min(num_gpus,MAX_GPUS)-1
           call dil_dev_buf_destroy(buf(dil_dev_num(DEV_NVIDIA_GPU,j0)),je)
          enddo
         endif
 !MIC buffers:
         if(mic_on) then
          do j0=0,min(num_mics,MAX_MICS)-1
           call dil_dev_buf_destroy(buf(dil_dev_num(DEV_INTEL_MIC,j0)),je)
          enddo
         endif
         return
         end subroutine cleanup

        end subroutine dil_tensor_contract_pipe
!---------------------------------------------------------------------------------------------------
        subroutine dil_tensor_contract(tcontr,globality,mem_lim,ierr,locked,async,num_gpus,num_mics) !PARALLEL (MPI)
!This is a user-level API subroutine for performing pipelined tensor contractions.
! # Argument <globality> determines the globality kind of the tensor contraction:
!    .FALSE.: Each MPI process entering here is assumed to have its own unique piece of work
!             specified via its own <tcontr%contr_spec>;
!     .TRUE.: <tcontr%contr_spec> specifies the full tensor contraction that will be
!             split into parts here, each part (piece of work) assigned to an MPI process.
! # If <async> is present and TRUE, one will need to call <dil_tensor_contract_finalize> later
!   to finalize this tensor contraction (for each MPI process). In this case, it is errorneous
!   to free or reuse <tcontr%buffer(:)> until the later call to <dil_tensor_contract_finalize> returns.
!   The full tensor contraction handle <tcontr> cannot be reused until that point!
        implicit none
        type(dil_tens_contr_t), target, intent(inout):: tcontr !in: full tensor contraction specification
        logical, intent(in):: globality                        !in: globality kind of the tensor contraction
        integer(INTL), intent(in):: mem_lim                    !in: local buffer memory limit in bytes
        integer(INTD), intent(inout):: ierr                    !out: error (0:success)
        logical, intent(in), optional:: locked                 !in: if .true., MPI windows for tensor arguments are assumed locked
        logical, intent(in), optional:: async                  !in: if .TRUE., the tensor contraction will not be finalized here
        integer(INTD), intent(in), optional:: num_gpus         !in: number of Nvidia GPUs to utilize (0..num_gpus-1)
        integer(INTD), intent(in), optional:: num_mics         !in: number of Intel MICs to utilize (0..num_mics-1)
        type(contr_spec_t):: cspec
        integer(INTD):: i,j,k,l,m,n,ngpus,nmics,nd,nl,nr,impis,impir,impir_world
        integer(INTL):: tcbv
        logical:: win_lck,asncr

        ierr=0
        impis=my_mpi_size(infpar%lg_comm); if(impis.le.0) then; ierr=1; return; endif !size of the local MPI communicator
        impir=my_mpi_rank(infpar%lg_comm); if(impir.lt.0) then; ierr=2; return; endif !rank in the local MPI communicator
        impir_world=my_mpi_rank() !rank in MPI_COMM_WORLD
        if(present(locked)) then; win_lck=locked; else; win_lck=.false.; endif
        if(present(async)) then; asncr=async; else; asncr=.false.; endif
        if(present(num_gpus)) then; ngpus=max(num_gpus,0); else; ngpus=0; endif
        if(present(num_mics)) then; nmics=max(num_mics,0); else; nmics=0; endif
        if(DIL_DEBUG)&
        &write(CONS_OUT,'("#DEBUG(dil_tensor_contract): Entered: Process ",i6," of ",i6,": Locked ",l1,": Async ",l1,": ",'//&
        &'i2," GPUs, ",i2," MICs ...")') impir,impis,win_lck,asncr,ngpus,nmics
!Allocate the work buffer, if needed:
        if(tcontr%alloc_type.eq.DIL_ALLOC_NOT) then
         call dil_prepare_buffer(tcontr,mem_lim,ierr)
         if(ierr.eq.0) then
          if(DIL_DEBUG) write(CONS_OUT,'("#DEBUG(dil_tensor_contract): Work buffer allocated: Volume = ",i12)')&
          &size(tcontr%buffer)
         else
          if(VERBOSE) write(CONS_OUT,'("#ERROR(dil_tensor_contract): Work buffer allocation failed: ",i9,1x,i12)') ierr,mem_lim
          ierr=3; return
         endif
        else
         if(associated(tcontr%buffer)) then
          tcbv=size(tcontr%buffer)
          if(tcbv*tensor_dp.ge.mem_lim) then
           if(DIL_DEBUG) write(CONS_OUT,'("#DEBUG(dil_tensor_contract): Preallocated work buffer volume = ",i12)') tcbv
          else
           if(VERBOSE) write(CONS_OUT,'("#ERROR(dil_tensor_contract): Work buffer is not large enough: ",i12,1x,i12)')&
           &tcbv*tensor_dp,mem_lim
           ierr=4; return
          endif
         else
          if(VERBOSE) write(CONS_OUT,'("#ERROR(dil_tensor_contract): Work buffer is expected to be associated, but not!")')
          ierr=5; return
         endif
        endif
!Global tensor contractions require work splitting:
        if(globality) then
 !Compute ranks of tensors:
         nd=tcontr%contr_spec%ndims_left+tcontr%contr_spec%ndims_right
         nl=tcontr%contr_spec%ndims_contr+tcontr%contr_spec%ndims_left
         nr=tcontr%contr_spec%ndims_contr+tcontr%contr_spec%ndims_right
 !Save global tensor dimensions/bases:
         cspec%ddims(1:nd)=tcontr%contr_spec%ddims(1:nd)
         cspec%ldims(1:nl)=tcontr%contr_spec%ldims(1:nl)
         cspec%rdims(1:nr)=tcontr%contr_spec%rdims(1:nr)
         cspec%dbase(1:nd)=tcontr%contr_spec%dbase(1:nd)
         cspec%lbase(1:nl)=tcontr%contr_spec%lbase(1:nl)
         cspec%rbase(1:nr)=tcontr%contr_spec%rbase(1:nr)
 !Partition the global tensor contraction space:
         call dil_tens_contr_distribute(tcontr,impis,impir,ierr)
         if(ierr.gt.0) then !ierr=-1 is OK, meaning NO WORK
          if(VERBOSE) write(CONS_OUT,'("#ERROR(dil_tensor_contract): Tensor contraction distribution failed: ",i9)') ierr
          ierr=6
         endif
        endif
!Execute tensor contraction:
        if(ierr.eq.0) then
         if(asncr) then !asynchronous execution (will need a call to <dil_tensor_contract_finalize> later)
          call dil_tensor_contract_pipe(tcontr%contr_spec,tcontr%dest_arg,tcontr%left_arg,tcontr%right_arg,&
                &tcontr%alpha,tcontr%beta,mem_lim,tcontr%buffer,ierr,locked=win_lck,&
                &nasync=tcontr%num_async,lasync=tcontr%list_async,num_gpus=ngpus,num_mics=nmics)
          if(ierr.ne.0) ierr=7
         else !blocking execution (will be finalized here)
          call dil_tensor_contract_pipe(tcontr%contr_spec,tcontr%dest_arg,tcontr%left_arg,tcontr%right_arg,&
                &tcontr%alpha,tcontr%beta,mem_lim,tcontr%buffer,ierr,locked=win_lck,num_gpus=ngpus,num_mics=nmics)
          if(ierr.ne.0) ierr=8
         endif
        elseif(ierr.lt.0) then !no work for this MPI process: Ok
         ierr=0
        endif
!Restore the global tensor contraction specification, if needed:
        if(globality) then
         tcontr%contr_spec%ddims(1:nd)=cspec%ddims(1:nd)
         tcontr%contr_spec%ldims(1:nl)=cspec%ldims(1:nl)
         tcontr%contr_spec%rdims(1:nr)=cspec%rdims(1:nr)
         tcontr%contr_spec%dbase(1:nd)=cspec%dbase(1:nd)
         tcontr%contr_spec%lbase(1:nl)=cspec%lbase(1:nl)
         tcontr%contr_spec%rbase(1:nr)=cspec%rbase(1:nr)
        endif
!Deallocate work buffer memory:
        if(.not.asncr) then
         if(tcontr%alloc_type.ne.DIL_ALLOC_EXT.and.associated(tcontr%buffer)) then
          call cpu_ptr_free(tcontr%buffer,ierr,attr=tcontr%alloc_type)
          if(ierr.eq.0) then
           tcontr%alloc_type=DIL_ALLOC_NOT
           if(DIL_DEBUG) write(CONS_OUT,'("#DEBUG(dil_tensor_contract): Work buffer freed: ",l1)') associated(tcontr%buffer)
          else
           if(VERBOSE) write(CONS_OUT,'("#ERROR(dil_tensor_contract): Work buffer deallocation failed: ",i11)') ierr
           ierr=9
          endif
         endif
        endif
        if(DIL_DEBUG)&
        &write(CONS_OUT,'("#DEBUG(dil_tensor_contract): Exited: Process ",i6," of ",i6,": Status ",i9)') impir,impis,ierr
        return
        end subroutine dil_tensor_contract
!------------------------------------------------------------------------------
        subroutine dil_tensor_contract_finalize(tcontr,ierr,locked,free_buffer) !PARALLEL (MPI)
!This subroutine finalizes any outstanding MPI communications (uploads) associated
!with a non-blocking tensor contraction represented by <tcontr>. It is safe
!to finalize blocking tensor contractions here as well.
        implicit none
        type(dil_tens_contr_t), intent(inout):: tcontr !inout: full tensor contraction specification
        integer(INTD), intent(inout):: ierr            !out: error (0:success)
        logical, intent(in), optional:: locked         !in: if .true., MPI windows for tensor arguments are assumed locked
        logical, intent(in), optional:: free_buffer    !in: if .true., the internal work buffer will be deallocated
        integer(INTD):: i
        logical:: win_lck

        ierr=0
        if(present(locked)) then; win_lck=locked; else; win_lck=.false.; endif
        do i=1,tcontr%num_async
         if(win_lck) then
          call tensor_mpi_win_flush(int(tcontr%list_async(i)%window,tensor_mpi_kind),&
             &int(tcontr%list_async(i)%rank,tensor_mpi_kind))
         else
          call tensor_mpi_win_unlock(int(tcontr%list_async(i)%rank,tensor_mpi_kind),&
             &int(tcontr%list_async(i)%window,tensor_mpi_kind))
         endif
        enddo
        if(DIL_DEBUG) write(CONS_OUT,'("#DEBUG(dil_tensor_contract_finalize): Number of finalized MPI uploads = ",i6)')&
                      &tcontr%num_async
        tcontr%num_async=0
        if(present(free_buffer)) then
         if(free_buffer) then
          if(associated(tcontr%buffer)) call cpu_ptr_free(tcontr%buffer,ierr,attr=tcontr%alloc_type)
          if(ierr.eq.0) then
           tcontr%alloc_type=DIL_ALLOC_NOT
           if(DIL_DEBUG) write(CONS_OUT,'("#DEBUG(dil_tensor_contract_finalize): Work buffer freed: ",l1)')&
           &associated(tcontr%buffer)
          else
           if(VERBOSE) write(CONS_OUT,'("#ERROR(dil_tensor_contract_finalize): Unable to free the work buffer: ",i11)') ierr
          endif
         endif
        endif
        return
        end subroutine dil_tensor_contract_finalize
!-------------------------------------------------------
        real(tensor_dp) function dil_tensor_norm1(tens,ierr) !PARALLEL (MPI+OMP)
!This function computes 1-norm of a tensor.
        implicit none
        type(tensor), intent(in):: tens               !in: tensor
        integer(INTD), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INTL):: i,l
        real(tensor_dp):: nrm1
        integer(INTD):: mpi_dtyp,errc

        errc=0; dil_tensor_norm1=0E0_tensor_dp; nrm1=0E0_tensor_dp
        call tensor_mpi_barrier(infpar%lg_comm) !complete all outstanding communications
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,l) REDUCTION(+:nrm1)
        do i=1,tens%nlti
!$OMP DO SCHEDULE(GUIDED)
         do l=1,tens%ti(i)%e
          nrm1=nrm1+abs(tens%ti(i)%t(l))
         enddo
!$OMP END DO
        enddo
!$OMP END PARALLEL
        if(tensor_dp.eq.8) then
         mpi_dtyp=MPI_REAL8
        elseif(tensor_dp.eq.4) then
         mpi_dtyp=MPI_REAL4
        else
         errc=1
        endif
        if(errc.eq.0) call MPI_ALLREDUCE(nrm1,dil_tensor_norm1,1_INTD,mpi_dtyp,MPI_SUM,infpar%lg_comm,errc)
        if(present(ierr)) ierr=errc
        return
        end function dil_tensor_norm1
!--------------------------------------------------------------
        real(tensor_dp) function dil_array_norm1(arr,arr_size,ierr) !PARALLEL (OMP)
        implicit none
        real(tensor_dp), intent(in):: arr(1:*)            !in: array
        integer(INTL), intent(in):: arr_size          !in: volume of the array (number of elements)
        integer(INTD), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INTL):: i
        integer(INTD):: errc
        real(tensor_dp):: v

        errc=0; dil_array_norm1=0E0_tensor_dp
        if(arr_size.gt.0) then
         v=0E0_tensor_dp
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED) REDUCTION(+:v)
         do i=1,arr_size
          v=v+abs(arr(i))
         enddo
!$OMP END PARALLEL DO
         dil_array_norm1=v
        else
         errc=1
        endif
        if(present(ierr)) ierr=errc
        return
        end function dil_array_norm1
!----------------------------------------
        subroutine dil_tensor_init(a,val) !SERIAL
!Each MPI process initializes its own tiles.
        implicit none
        type(tensor), intent(inout) :: a
        real(tensor_dp), intent(in), optional:: val
        integer(INTD):: lti
        real(tensor_dp):: vlu
 !get the slaves here
        if(a%access_type==AT_MASTER_ACCESS.and.infpar%lg_mynum==infpar%master) then
         call pdm_tensor_sync(infpar%lg_comm,JOB_tensor_ZERO,a)
        endif
 !loop over local tiles and zero them individually
        vlu=0E0_tensor_dp; if(present(val)) vlu=val
        do lti=1,a%nlti
!$OMP WORKSHARE
         a%ti(lti)%t=vlu
!$OMP END WORKSHARE
        enddo
        return
        end subroutine dil_tensor_init
!------------------------------------------
        subroutine dil_array_init(a,sz,val) !SERIAL
        implicit none
        real(tensor_dp), intent(inout):: a(1:*)     !out: array (contiguous)
        integer(INTL), intent(in):: sz          !in: number of elements to initialize
        real(tensor_dp), intent(in), optional:: val !in: initialization value (defaults to zero)
        integer(INTL):: i
        real(tensor_dp):: v

        if(present(val)) then; v=val; else; v=0E0_tensor_dp; endif
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
        do i=1,sz
         a(i)=v
        enddo
!$OMP END PARALLEL DO
        return
        end subroutine dil_array_init
!-----------------------------------------------
        subroutine dil_debug_to_file_start(ierr)
!This subroutine starts redirecting process's output to file.
        implicit none
        integer(INTD), intent(out), optional:: ierr
        character(128):: deb_fname
        integer(INTD):: i,impir_world,impir

        if(present(ierr)) ierr=0
        if(CONS_OUT_SAVED.le.0) then
         impir_world=my_mpi_rank(); impir=my_mpi_rank(infpar%lg_comm)
         deb_fname='dil_debug.'; call int2str(impir_world,deb_fname(11:),i); deb_fname(11+i:11+i+3)='.log'; i=11+i+3
         open(DIL_DEBUG_FILE,file=deb_fname(1:i),form='FORMATTED',status='UNKNOWN')
         CONS_OUT_SAVED=CONS_OUT; CONS_OUT=DIL_DEBUG_FILE; DIL_CONS_OUT=DIL_DEBUG_FILE
         write(CONS_OUT,'("### DEBUG BEGIN: Global Rank ",i7," (Local Rank ",i7,")")') impir_world,impir
        else
         if(present(ierr)) ierr=1 !debug file is already opened
        endif
        return
        end subroutine dil_debug_to_file_start
!---------------------------------------------
        subroutine dil_debug_to_file_finish(ierr)
!This subroutine starts redirecting process's output to file.
        implicit none
        integer(INTD), intent(out), optional:: ierr
        integer(INTD):: impir_world,impir

        if(present(ierr)) ierr=0
        if(CONS_OUT_SAVED.gt.0) then
         impir_world=my_mpi_rank(); impir=my_mpi_rank(infpar%lg_comm)
         write(CONS_OUT,'("### DEBUG END: Global Rank ",i7," (Local Rank ",i7,")")') impir_world,impir
         CONS_OUT=CONS_OUT_SAVED; DIL_CONS_OUT=CONS_OUT_SAVED; CONS_OUT_SAVED=0; close(DIL_DEBUG_FILE)
        else
         if(present(ierr)) ierr=1 !debug file has not been opened
        endif
        return
        end subroutine dil_debug_to_file_finish
!-------------------------------------------------------------------------------------
        subroutine dil_array_print(arr,ierr,fname,fhandle,print_first,sort,by_modulus) !SERIAL (I/O)
!This subroutine prints a real array into a file, with an optional ordering done to it.
        implicit none
        real(tensor_dp), intent(in), target:: arr(1:)         !inout: array to print
        integer(INTD), intent(inout), optional:: ierr     !out: error code (0:success)
        character(*), intent(in), optional:: fname        !in: EITHER a file name where to print
        integer(INTD), intent(in), optional:: fhandle     !in: OR a file handle of an already opened file
        integer(INTL), intent(in), optional:: print_first !in: print only the first <print_first> elements
        logical, intent(in), optional:: sort              !in: if .true., the array will be printed in a sorted order
        logical, intent(in), optional:: by_modulus        !in: if .true., the array will be sorted by modulus
        real(tensor_dp), pointer, contiguous:: sar(:)
        integer(INTL):: i,arr_size,pfe
        integer(INTD):: fh,errc
        logical:: alloc

        errc=0; arr_size=size(arr)
        if(arr_size.gt.0) then
         if(present(fhandle).and.(.not.present(fname))) then
          if(fhandle.gt.8) then !the first few handles are reserved
           fh=fhandle
          else
           errc=1; if(present(ierr)) ierr=errc
           return
          endif
         elseif(present(fname).and.(.not.present(fhandle))) then
          fh=DIL_TMP_FILE2
          open(fh,file=fname(1:len_trim(fname)),form='FORMATTED',status='UNKNOWN')
         elseif((.not.present(fname)).and.(.not.present(fhandle))) then
          fh=DIL_CONS_OUT
         else
          errc=2
         endif
         if(errc.eq.0) then
          alloc=.false.
          if(present(sort)) then
           if(sort.and.arr_size.gt.1) then
            allocate(sar(1:arr_size),STAT=errc)
            if(errc.eq.0) then
             alloc=.true.
             if(present(by_modulus)) then
              if(by_modulus) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
               do i=1,arr_size
                sar(i)=abs(arr(i))
               enddo
!$OMP END PARALLEL DO
              else
!$OMP WORKSHARE
               sar(1:arr_size)=arr(1:arr_size)
!$OMP END WORKSHARE
              endif
             else
!$OMP WORKSHARE
              sar(1:arr_size)=arr(1:arr_size)
!$OMP END WORKSHARE
             endif
             call merge_sort_real8(arr_size,sar,dir=-1_INTD) !sort in an descending order
            else
             errc=3
            endif
           else
            sar=>arr
           endif
          else
           sar=>arr
          endif
          if(errc.eq.0) then
           pfe=arr_size
           if(present(print_first)) pfe=print_first
           do i=1,pfe
            write(fh,'(D22.14)') sar(i)
           enddo
           if(alloc) deallocate(sar)
          endif
          if(present(fname)) close(fh)
          nullify(sar)
         endif
        else
         errc=4
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine dil_array_print
!--------------------------------------------------------------------------------------------------
        logical function dil_mm_pipe_efficient(ll,lr,lc,comp_bandwidth,data_bandwidth,data_latency) !SERIAL
!Given the tile dimensions, this function decides whether a matrix multiplication can be efficiently pipelined,
!that is, whether a matrix multiplication with the given tile dimensions can be efficiently
!pipelined on a given architecture (data transfer latency/bandwidth, PU Flops/s).
!Matrix tile multiplication: D(1:ll,1:lr)+=L(1:lc,1:ll)*R(1:lc,1:lr)
!No argument validity checks!
        implicit none
        integer(INTL), intent(in):: ll,lr,lc !in: matrix tile dimensions (left, right, contracted)
        real(8), intent(in):: comp_bandwidth !in: PU computational throughput (Flops/s)
        real(8), intent(in):: data_bandwidth !in: data transfer bandwidth (Words/s)
        real(8), intent(in):: data_latency   !in: data transfer latency (s)
        real(8), parameter:: much_more=5d0   !defines what "much more" exactly means
        real(8):: l,r,c,v

        l=real(ll,8); r=real(lr,8); c=real(lc,8)
        if(l*r*c.gt.(2d0*data_latency*comp_bandwidth)*much_more) then
         v=comp_bandwidth/data_bandwidth
         if(c.ge.v.and.l*r/(l+r).ge.v) then
          dil_mm_pipe_efficient=.true.
         else
          dil_mm_pipe_efficient=.false.
         endif
        else
         dil_mm_pipe_efficient=.false.
        endif
        return
        end function dil_mm_pipe_efficient
!---------------------------------------------------------------------------------------------------
        logical function dil_will_malloc_succeed(mem_bytes,page_size,hugepage_size,max_huge,max_mem) !SERIAL
!This function checks whether a given malloc request has a chance for success.
!If the arguments passed to this function are invalid, .FALSE. will be returned (no error status).
!NOTES:
! # Because of using the same file handle, this subroutine is not threadsafe,
!   that is, it cannot be called from multiple threads simulateneously.
! # The result returned is only a probable success (not 100% reliable) because
!   (a) the previously allocated memory might not have been touched yet (touch it, then call);
!   (b) the buddyinfo can become outdated due to a concurrent malloc() (do not call concurrent mallocs);
!   (c) the malloc() implementation may not be able to use all the pages
!       to produce an allocation of a requested size (see parameter RELIABLE_PART below).
        implicit none
        integer(INTL), intent(in):: mem_bytes               !in: number of bytes to be allocated
        integer(INTL), intent(in), optional:: page_size     !in: basic page size in bytes (default is 4K)
        integer(INTL), intent(in), optional:: hugepage_size !in: huge page size in bytes (defaults to 2M)
        integer(INTL), intent(out), optional:: max_huge     !out: maximum available memory (bytes) backed with huge pages
        integer(INTL), intent(out), optional:: max_mem      !out: total maximum available memory (bytes)
!-------------------------------------------------------
        real(8), parameter:: RELIABLE_PART=1d0 !empiric parameter to account for the non-ideality of the buddy malloc()
        integer(INTL), parameter:: DEFAULT_PAGE=4096 !default basic page size in bytes
        integer(INTL), parameter:: DEFAULT_HUGEPAGE=2097152 !default hugepage size in bytes
        integer(INTD), parameter:: MAX_BUDDY_LEVELS=128 !max anticipated number of buddy levels
        integer(INTL), parameter:: MEMINFO_UNIT=1024 !number of bytes in a /proc/meminfo memory measurement unit
!------------------------------------------------------
        integer(INTL):: psz,hsz,pls,ahpm,totm,buds(0:MAX_BUDDY_LEVELS-1)
        character(512):: str
        integer(INTD):: i,k,l,m,n,words(2,MAX_BUDDY_LEVELS+16)

        dil_will_malloc_succeed=.false.; ahpm=-1; totm=-1
        psz=DEFAULT_PAGE; hsz=DEFAULT_HUGEPAGE
        if(present(page_size)) psz=page_size
        if(present(hugepage_size)) hsz=hugepage_size
        if(psz.le.0.or.hsz.lt.psz.or.mod(hsz,psz).ne.0) return
!Read current /proc/buddyinfo (may change at any time):
        open(DIL_TMP_FILE1,file='/proc/buddyinfo',form='FORMATTED',status='OLD',ERR=999)
        buds(:)=0; m=0; i=0; str=' '
        do
         read(DIL_TMP_FILE1,'(A512)',END=100) str; l=len_trim(str)
         if(l.gt.0) then
          call str_parse(str,' ,',n,words,ierr=i,str_len=l); if(i.ne.0) exit
          call fill_buddy_info(k,i); if(i.ne.0) exit
          m=max(m,k)
          str(1:l)=' '
         endif
        enddo
100     close(DIL_TMP_FILE1)
!Compute the free memory amount:
        if(i.eq.0) then
 !Compute the amount of hugepage backed free memory (only for the mmap path):
         pls=psz; ahpm=0 !ahpm: available hugepage memory in bytes
         do l=0,m-1
          if(pls.ge.hsz) ahpm=ahpm+pls*buds(l) !count only memory chunks larger or equal to the hugepage size
          pls=pls*2
         enddo
!         if(DIL_DEBUG) write(CONS_OUT,'("#DEBUG(DIL): malloc requested ",i12," B from hugepage backed RAM of ",i12)')&
!         &mem_bytes,ahpm
         if(mem_bytes.lt.int(real(ahpm,8)*RELIABLE_PART,INTL)) dil_will_malloc_succeed=.true.
 !Compute the total amount of free memory (only for the mmap path):
         totm=0 !totm: total available memory in bytes
         open(DIL_TMP_FILE1,file='/proc/meminfo',form='FORMATTED',status='OLD',ERR=888)
         i=0; str=' '
         do
          read(DIL_TMP_FILE1,'(A512)',END=200) str; l=len_trim(str)
          if(l.gt.0) then
           call str_parse(str,' ',n,words,ierr=i,str_len=l); if(i.ne.0) exit
           call collect_total_mem(i); if(i.ne.0) exit
           str(1:l)=' '
          endif
         enddo
200      close(DIL_TMP_FILE1)
        endif
        if(present(max_huge)) max_huge=ahpm
        if(present(max_mem)) max_mem=totm
        return
888     if(VERBOSE) write(CONS_OUT,'("#ERROR(tensor_algebra_dil::dil_will_malloc_succeed): unable to open /proc/meminfo!")')
        if(present(max_huge)) max_huge=ahpm
        if(present(max_mem)) max_mem=totm
        return
999     if(VERBOSE) write(CONS_OUT,'("#ERROR(tensor_algebra_dil::dil_will_malloc_succeed): unable to open /proc/buddyinfo!")')
        if(present(max_huge)) max_huge=ahpm
        if(present(max_mem)) max_mem=totm
        return
        contains

         subroutine collect_total_mem(errc)
         integer(INTD), intent(out):: errc
         integer(INTD):: jb,je
         errc=0_INTD
         if(n.eq.3) then !three fields must be present
          jb=words(1,1); je=words(2,1)
          if(je-jb+1.eq.len('MemFree:')) then
           if(str(jb:je).eq.'MemFree:') then
            totm=totm+str2int(str(jb:je),je-jb+1_INTD,errc)*MEMINFO_UNIT; if(errc.ne.0_INTD) totm=-1
            return
           endif
          endif
          if(je-jb+1.eq.len('Buffers:')) then
           if(str(jb:je).eq.'Buffers:') then
            totm=totm+str2int(str(jb:je),je-jb+1_INTD,errc)*MEMINFO_UNIT; if(errc.ne.0_INTD) totm=-1
            return
           endif
          endif
          if(je-jb+1.eq.len('Cached:')) then
           if(str(jb:je).eq.'Cached:') then
            totm=totm+str2int(str(jb:je),je-jb+1_INTD,errc)*MEMINFO_UNIT; if(errc.ne.0_INTD) totm=-1
            return
           endif
          endif
          if(je-jb+1.eq.len('SwapFree:')) then
           if(str(jb:je).eq.'SwapFree:') then
            totm=totm+str2int(str(jb:je),je-jb+1_INTD,errc)*MEMINFO_UNIT; if(errc.ne.0_INTD) totm=-1
            return
           endif
          endif
         else
          errc=-1_INTD
         endif
         return
         end subroutine collect_total_mem

         subroutine fill_buddy_info(jl,errc)
         integer(INTD), intent(out):: jl
         integer(INTD), intent(out):: errc
         integer(INTD):: j0,jb,je,js
         errc=0; jl=-1
         do j0=1,n
          jb=words(1,j0); je=words(2,j0); js=je-jb+1_INTD
          if(jl.ge.0) then
           if(jl.ge.MAX_BUDDY_LEVELS) then
            if(VERBOSE) write(CONS_OUT,'("#ERROR(tensor_algebra_dil::dil_will_malloc_succeed): MAX_BUDDY_LEVELS exceeded!")')
            errc=-1; return
           endif
           buds(jl)=buds(jl)+str2int(str(jb:je),js,errc); if(errc.ne.0) return
           jl=jl+1
          else
           if(js.eq.len('Normal')) then
            if(str(jb:je).eq.'Normal') jl=0 !start recording
           endif
          endif
         enddo
         return
         end subroutine fill_buddy_info

        end function dil_will_malloc_succeed
!=============================================
        subroutine dil_test(dtens,ltens,rtens)
        implicit none
        type(tensor), intent(inout), target:: dtens
        type(tensor), intent(inout), target:: ltens
        type(tensor), intent(inout), target:: rtens
        integer(INTD):: i,j,k,l,m,n,impir,impir_world,drank,lrank,rrank,errc
        integer(INTL):: mem_lim,dvol,lvol,rvol
        integer(tensor_mpi_kind):: mpi_err
        integer(INTD):: dbas(1:MAX_TENSOR_RANK),lbas(1:MAX_TENSOR_RANK),rbas(1:MAX_TENSOR_RANK)
        integer(INTD):: ddim(1:MAX_TENSOR_RANK),ldim(1:MAX_TENSOR_RANK),rdim(1:MAX_TENSOR_RANK)
        integer(INTD):: dful(1:MAX_TENSOR_RANK),lful(1:MAX_TENSOR_RANK),rful(1:MAX_TENSOR_RANK)
        real(tensor_dp), pointer, contiguous:: darr(:),larr(:),rarr(:),barr(:)
        type(dil_tens_contr_t):: tcr
        character(128):: tcs
        real(tensor_dp):: val0,val1,val2
        real(8):: tmb,tm

        errc=0
        impir=my_mpi_rank(infpar%lg_comm); impir_world=my_mpi_rank()
        call dil_debug_to_file_start()
!Check tensor kernels:
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: Checking tensor algebra kernels ...")') impir
        drank=5; lrank=5; rrank=5
        ddim(1:drank)=(/25_INTD,35_INTD,15_INTD,20_INTD,10_INTD/)
        ldim(1:lrank)=(/11_INTD,23_INTD,7_INTD,4_INTD,3_INTD/)
        rdim(1:rrank)=(/3_INTD,4_INTD,7_INTD,23_INTD,11_INTD/)
        dvol=1_INTL; do i=1,drank; dvol=dvol*ddim(i); enddo
        lvol=1_INTL; do i=1,lrank; lvol=lvol*ldim(i); enddo
        rvol=1_INTL; do i=1,rrank; rvol=rvol*rdim(i); enddo
        allocate(barr(1:dvol),darr(1:dvol),larr(1:lvol),rarr(1:rvol))
        call random_number(barr(1:dvol)); darr(1:dvol)=barr(1:dvol)
        print *,'Two Random Numbers: ',barr(1),barr(dvol)
        call dil_tensor_slice(drank,darr,ddim,larr,ldim,(/2_INTD,5_INTD,1_INTD,0_INTD,3_INTD/),errc)
        print *,'SLICE ERR = ',errc
        call dil_tensor_transpose(drank,ldim,(/5_INTD,4_INTD,3_INTD,2_INTD,1_INTD/),larr,rarr,errc)
        print *,'TRN 1 ERR = ',errc
        larr(1:lvol)=5d0
        call dil_tensor_insert(drank,darr,ddim,larr,ldim,(/2_INTD,5_INTD,1_INTD,0_INTD,3_INTD/),errc)
        print *,'INSERT 1 ERR = ',errc
        call dil_tensor_transpose(drank,rdim,(/5_INTD,4_INTD,3_INTD,2_INTD,1_INTD/),rarr,larr,errc)
        print *,'TRN 2 ERR = ',errc
        call dil_tensor_insert(drank,darr,ddim,larr,ldim,(/2_INTD,5_INTD,1_INTD,0_INTD,3_INTD/),errc)
        print *,'INSERT 2 ERR = ',errc
        val0=0d0; do i=1,dvol; val0=val0+abs(darr(i)-barr(i)); enddo
        val1=0d0; do i=1,dvol; val1=val1+abs(barr(i)); enddo
        val2=0d0; do i=1,dvol; val2=val2+abs(darr(i)); enddo
        deallocate(barr); deallocate(darr); deallocate(larr); deallocate(rarr)
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: Done: Result =",3(1x,D25.15))') impir,val0,val1,val2
        flush(CONS_OUT)
!Init tensors:
        call tensor_mpi_barrier(infpar%lg_comm)
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: Input tensor init started ...")') impir
        call dil_tensor_init(ltens,1d-2)
        call dil_tensor_init(rtens,1d-3)
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: Input tensor init finished.")') impir
        call tensor_mpi_barrier(infpar%lg_comm)
!Compute input tensor norms:
        val1=tensor_tiled_pdm_get_nrm2(ltens)
        call tensor_mpi_barrier(infpar%lg_comm)
        val2=tensor_tiled_pdm_get_nrm2(rtens)
        call tensor_mpi_barrier(infpar%lg_comm)
        write(CONS_OUT,'("#DEBUG(DIL): Left/Right 2-norms:",2(1x,D25.15))') val1,val2
        tcs='D(a,b,c,d)=L(e,f,d,b)*R(f,a,c,e)'
!ALL 'DDD':
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: All-DDD tensor contraction setup started ...")') impir
        call dil_tensor_init(dtens)
        call tensor_mpi_barrier(infpar%lg_comm)
        val0=tensor_tiled_pdm_get_nrm2(dtens)
        call tensor_mpi_barrier(infpar%lg_comm)
        write(CONS_OUT,'("#DEBUG(DIL): Destination norm at the beginning: ",D25.15)') val0
        call dil_clean_tens_contr(tcr)
        call dil_set_tens_contr_args(tcr,'d',errc,tens_distr=dtens); print *,'DTENS ERR = ',errc
        call dil_set_tens_contr_args(tcr,'l',errc,tens_distr=ltens); print *,'LTENS ERR = ',errc
        call dil_set_tens_contr_args(tcr,'r',errc,tens_distr=rtens); print *,'RTENS ERR = ',errc
        call dil_set_tens_contr_spec(tcr,tcs,errc); print *,'CSPEC ERR = ',errc
        mem_lim=dil_get_min_buf_size(tcr,errc); print *,'BUF SIZE ERR = ',errc
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: Tensor contraction setup finished: Memory limit = ",i12," bytes")') impir,mem_lim
        call tensor_mpi_barrier(infpar%lg_comm)
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: Tensor contraction started:")') impir
        if(impir.eq.0) write(*,'("Global DDD started ... ")',ADVANCE='NO')
        tmb=process_wtime()
        call dil_tensor_contract(tcr,DIL_TC_ALL,mem_lim,errc)
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: Tensor contraction finished: Status ",i9)') impir,errc
        call tensor_mpi_barrier(infpar%lg_comm)
        tm=process_wtime(tmb); if(impir.eq.0) write(*,'("Done in ",F10.4," s.")') tm
        val0=tensor_tiled_pdm_get_nrm2(dtens)
        call tensor_mpi_barrier(infpar%lg_comm)
        write(CONS_OUT,'("#DEBUG(DIL): Destination norm at the end: ",D25.15)') val0
        if(impir.eq.0) write(*,'("Destination norm at the end: ",D25.15)') val0
!        goto 999
!General settings:
        drank=dtens%mode; lrank=ltens%mode; rrank=rtens%mode
        print *,'Tensor ranks: ',drank,lrank,rrank
        dful(1:drank)=(/180_INTD,180_INTD,40_INTD,40_INTD/)
        lful(1:lrank)=(/40_INTD,180_INTD,40_INTD,180_INTD/)
        rful(1:rrank)=(/180_INTD,180_INTD,40_INTD,40_INTD/)
        select case(impir)
        case(0)
         dbas(1:drank)=(/0_INTD,0_INTD,0_INTD,0_INTD/)
         lbas(1:lrank)=(/0_INTD,0_INTD,0_INTD,0_INTD/)
         rbas(1:rrank)=(/0_INTD,0_INTD,0_INTD,0_INTD/)
         ddim(1:drank)=(/60_INTD,180_INTD,28_INTD,40_INTD/)
         ldim(1:lrank)=(/40_INTD,180_INTD,40_INTD,180_INTD/)
         rdim(1:rrank)=(/180_INTD,60_INTD,28_INTD,40_INTD/)
        case(1)
         dbas(1:drank)=(/0_INTD,0_INTD,28_INTD,0_INTD/)
         lbas(1:lrank)=(/0_INTD,0_INTD,0_INTD,0_INTD/)
         rbas(1:rrank)=(/0_INTD,0_INTD,28_INTD,0_INTD/)
         ddim(1:drank)=(/60_INTD,180_INTD,12_INTD,40_INTD/)
         ldim(1:lrank)=(/40_INTD,180_INTD,40_INTD,180_INTD/)
         rdim(1:rrank)=(/180_INTD,60_INTD,12_INTD,40_INTD/)
        case(2)
         dbas(1:drank)=(/60_INTD,0_INTD,0_INTD,0_INTD/)
         lbas(1:lrank)=(/0_INTD,0_INTD,0_INTD,0_INTD/)
         rbas(1:rrank)=(/0_INTD,60_INTD,0_INTD,0_INTD/)
         ddim(1:drank)=(/120_INTD,180_INTD,28_INTD,40_INTD/)
         ldim(1:lrank)=(/40_INTD,180_INTD,40_INTD,180_INTD/)
         rdim(1:rrank)=(/180_INTD,120_INTD,28_INTD,40_INTD/)
        case(3)
         dbas(1:drank)=(/60_INTD,0_INTD,28_INTD,0_INTD/)
         lbas(1:lrank)=(/0_INTD,0_INTD,0_INTD,0_INTD/)
         rbas(1:rrank)=(/0_INTD,60_INTD,28_INTD,0_INTD/)
         ddim(1:drank)=(/120_INTD,180_INTD,12_INTD,40_INTD/)
         ldim(1:lrank)=(/40_INTD,180_INTD,40_INTD,180_INTD/)
         rdim(1:rrank)=(/180_INTD,120_INTD,12_INTD,40_INTD/)
        end select
        dvol=1_INTL; do i=1,drank; dvol=dvol*dful(i); enddo
        lvol=1_INTL; do i=1,lrank; lvol=lvol*lful(i); enddo
        rvol=1_INTL; do i=1,rrank; rvol=rvol*rful(i); enddo
        allocate(darr(1:dvol),STAT=i)
        if(i.ne.0) then
         write(*,'("MEM ALLOC 1 failed!")')
         call MPI_ABORT(infpar%lg_comm,0_tensor_mpi_kind,mpi_err)
        endif
        allocate(larr(1:lvol),STAT=i)
        if(i.ne.0) then
         write(*,'("MEM ALLOC 2 failed!")')
         call MPI_ABORT(infpar%lg_comm,0_tensor_mpi_kind,mpi_err)
        endif
        allocate(rarr(1:rvol),STAT=i)
        if(i.ne.0) then
         write(*,'("MEM ALLOC 3 failed!")')
         call MPI_ABORT(infpar%lg_comm,0_tensor_mpi_kind,mpi_err)
        endif
        do i=1,lvol; larr(i)=1d-2; enddo
        do i=1,rvol; rarr(i)=1d-3; enddo
        call tensor_mpi_barrier(infpar%lg_comm)
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: Initialization finished.")') impir
        flush(CONS_OUT)
!Tensor contraction 'DDD':
        mem_lim=384*1048576_INTL
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: DDD tensor contraction setup started (MEM_LIM = ",i12," B) ...")') impir,mem_lim
        call dil_tensor_init(dtens)
        call tensor_mpi_barrier(infpar%lg_comm)
        val0=tensor_tiled_pdm_get_nrm2(dtens)
        call tensor_mpi_barrier(infpar%lg_comm)
        write(CONS_OUT,'("#DEBUG(DIL): Destination norm: ",D25.15)') val0
        call dil_clean_tens_contr(tcr)
        call dil_set_tens_contr_args(tcr,'d',errc,tens_distr=dtens); print *,'DTENS ERR = ',errc
        call dil_set_tens_contr_args(tcr,'l',errc,tens_distr=ltens); print *,'LTENS ERR = ',errc
        call dil_set_tens_contr_args(tcr,'r',errc,tens_distr=rtens); print *,'RTENS ERR = ',errc
        call dil_set_tens_contr_spec(tcr,tcs,errc,ddim,ldim,rdim,dbas,lbas,rbas); print *,'CSPEC ERR = ',errc
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: Tensor contraction setup finished.")') impir
        call tensor_mpi_barrier(infpar%lg_comm)
        call dil_tensor_contract(tcr,DIL_TC_EACH,mem_lim,errc)
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: Tensor contraction finished: Status ",i9)') impir,errc
        call tensor_mpi_barrier(infpar%lg_comm)
        val0=tensor_tiled_pdm_get_nrm2(dtens)
        call tensor_mpi_barrier(infpar%lg_comm)
        write(CONS_OUT,'("#DEBUG(DIL): Destination norm: ",D25.15)') val0
!Tensor contraction 'LLL':
        mem_lim=1024*1048576_INTL
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: LLL tensor contraction setup started (MEM_LIM = ",i12," B) ...")') impir,mem_lim
        do i=1,dvol; darr(i)=0d0; enddo
        call tensor_mpi_barrier(infpar%lg_comm)
        call dil_clean_tens_contr(tcr)
        call dil_set_tens_contr_args(tcr,'d',errc,drank,dful,darr); print *,'DTENS ERR = ',errc
        call dil_set_tens_contr_args(tcr,'l',errc,lrank,lful,larr); print *,'DTENS ERR = ',errc
        call dil_set_tens_contr_args(tcr,'r',errc,rrank,rful,rarr); print *,'DTENS ERR = ',errc
        call dil_set_tens_contr_spec(tcr,tcs,errc,ddim,ldim,rdim,dbas,lbas,rbas); print *,'CSPEC ERR = ',errc
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: Tensor contraction setup finished.")') impir
        call tensor_mpi_barrier(infpar%lg_comm)
        call dil_tensor_contract(tcr,DIL_TC_EACH,mem_lim,errc)
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: Tensor contraction finished: Status ",i9)') impir,errc
        call tensor_mpi_barrier(infpar%lg_comm)
        val0=0d0; do i=1,dvol; val0=val0+darr(i)**2; enddo
        print *,'Local DTENS (2-norm)^2 = ',val0,impir
        call tensor_mpi_barrier(infpar%lg_comm)
!Tensor contraction 'LDD':
        mem_lim=1024*1048576_INTL
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: LDD tensor contraction setup started (MEM_LIM = ",i12," B) ...")') impir,mem_lim
        do i=1,dvol; darr(i)=0d0; enddo
        call dil_clean_tens_contr(tcr)
        call dil_set_tens_contr_args(tcr,'d',errc,drank,dful,darr); print *,'DTENS ERR = ',errc
        call dil_set_tens_contr_args(tcr,'l',errc,tens_distr=ltens); print *,'LTENS ERR = ',errc
        call dil_set_tens_contr_args(tcr,'r',errc,tens_distr=rtens); print *,'RTENS ERR = ',errc
        call dil_set_tens_contr_spec(tcr,tcs,errc,ddim,ldim,rdim,dbas,lbas,rbas); print *,'CSPEC ERR = ',errc
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: Tensor contraction setup finished.")') impir
        call tensor_mpi_barrier(infpar%lg_comm)
        call dil_tensor_contract(tcr,DIL_TC_EACH,mem_lim,errc)
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: Tensor contraction finished: Status ",i9)') impir,errc
        call tensor_mpi_barrier(infpar%lg_comm)
        val0=0d0; do i=1,dvol; val0=val0+darr(i)**2; enddo
        print *,'Local DTENS (2-norm)^2 = ',val0,impir
        call tensor_mpi_barrier(infpar%lg_comm)
!Tensor contraction 'DLL':
        mem_lim=1024*1048576_INTL
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: DLL tensor contraction setup started (MEM_LIM = ",i12," B) ...")') impir,mem_lim
        call dil_tensor_init(dtens)
        call tensor_mpi_barrier(infpar%lg_comm)
        val0=tensor_tiled_pdm_get_nrm2(dtens)
        call tensor_mpi_barrier(infpar%lg_comm)
        write(CONS_OUT,'("#DEBUG(DIL): Destination norm: ",D25.15)') val0
        call dil_clean_tens_contr(tcr)
        call dil_set_tens_contr_args(tcr,'d',errc,tens_distr=dtens); print *,'DTENS ERR = ',errc
        call dil_set_tens_contr_args(tcr,'l',errc,lrank,lful,larr); print *,'DTENS ERR = ',errc
        call dil_set_tens_contr_args(tcr,'r',errc,rrank,rful,rarr); print *,'DTENS ERR = ',errc
        call dil_set_tens_contr_spec(tcr,tcs,errc,ddim,ldim,rdim,dbas,lbas,rbas); print *,'CSPEC ERR = ',errc
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: Tensor contraction setup finished.")') impir
        call tensor_mpi_barrier(infpar%lg_comm)
        call dil_tensor_contract(tcr,DIL_TC_EACH,mem_lim,errc)
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: Tensor contraction finished: Status ",i9)') impir,errc
        call tensor_mpi_barrier(infpar%lg_comm)
        val0=tensor_tiled_pdm_get_nrm2(dtens)
        call tensor_mpi_barrier(infpar%lg_comm)
        write(CONS_OUT,'("#DEBUG(DIL): Destination norm: ",D25.15)') val0
!Tensor contraction 'LLD':
        mem_lim=1024*1048576_INTL
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: LLD tensor contraction setup started (MEM_LIM = ",i12," B) ...")') impir,mem_lim
        do i=1,dvol; darr(i)=0d0; enddo
        call dil_clean_tens_contr(tcr)
        call dil_set_tens_contr_args(tcr,'d',errc,drank,dful,darr); print *,'DTENS ERR = ',errc
        call dil_set_tens_contr_args(tcr,'l',errc,lrank,lful,larr); print *,'DTENS ERR = ',errc
        call dil_set_tens_contr_args(tcr,'r',errc,tens_distr=rtens); print *,'RTENS ERR = ',errc
        call dil_set_tens_contr_spec(tcr,tcs,errc,ddim,ldim,rdim,dbas,lbas,rbas); print *,'CSPEC ERR = ',errc
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: Tensor contraction setup finished.")') impir
        call tensor_mpi_barrier(infpar%lg_comm)
        call dil_tensor_contract(tcr,DIL_TC_EACH,mem_lim,errc)
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: Tensor contraction finished: Status ",i9)') impir,errc
        call tensor_mpi_barrier(infpar%lg_comm)
        val0=0d0; do i=1,dvol; val0=val0+darr(i)**2; enddo
        print *,'Local DTENS (2-norm)^2 = ',val0,impir
        call tensor_mpi_barrier(infpar%lg_comm)
!Tensor contraction 'DLD':
        mem_lim=1024*1048576_INTL
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: DLD tensor contraction setup started (MEM_LIM = ",i12," B) ...")') impir,mem_lim
        call dil_tensor_init(dtens)
        call tensor_mpi_barrier(infpar%lg_comm)
        val0=tensor_tiled_pdm_get_nrm2(dtens)
        call tensor_mpi_barrier(infpar%lg_comm)
        write(CONS_OUT,'("#DEBUG(DIL): Destination norm: ",D25.15)') val0
        call dil_clean_tens_contr(tcr)
        call dil_set_tens_contr_args(tcr,'d',errc,tens_distr=dtens); print *,'DTENS ERR = ',errc
        call dil_set_tens_contr_args(tcr,'l',errc,lrank,lful,larr); print *,'DTENS ERR = ',errc
        call dil_set_tens_contr_args(tcr,'r',errc,tens_distr=rtens); print *,'RTENS ERR = ',errc
        call dil_set_tens_contr_spec(tcr,tcs,errc,ddim,ldim,rdim,dbas,lbas,rbas); print *,'CSPEC ERR = ',errc
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: Tensor contraction setup finished.")') impir
        call tensor_mpi_barrier(infpar%lg_comm)
        call dil_tensor_contract(tcr,DIL_TC_EACH,mem_lim,errc)
        write(CONS_OUT,'("#DEBUG(DIL)[",i2,"]: Tensor contraction finished: Status ",i9)') impir,errc
        call tensor_mpi_barrier(infpar%lg_comm)
        val0=tensor_tiled_pdm_get_nrm2(dtens)
        call tensor_mpi_barrier(infpar%lg_comm)
        write(CONS_OUT,'("#DEBUG(DIL): Destination norm: ",D25.15)') val0
        val0=dil_tensor_norm1(dtens,errc); print *,'NORM1 ERR = ',errc
        write(CONS_OUT,'("#DEBUG(DIL): Destination norm1: ",D25.15)') val0
!Exit:
999     call dil_debug_to_file_finish()
        if(associated(darr)) deallocate(darr)
        if(associated(larr)) deallocate(larr)
        if(associated(rarr)) deallocate(rarr)
        call tensor_mpi_barrier(infpar%lg_comm)
        call MPI_ABORT(infpar%lg_comm,0_tensor_mpi_kind,mpi_err)
        return
        end subroutine dil_test
!DIL_ACTIVE (assumes Fortran-2003/2008, MPI-3):
#else
       contains
        subroutine tensor_alg_dil_dummy()
        end subroutine tensor_alg_dil_dummy
!DIL_ACTIVE (assumes Fortran-2003/2008, MPI-3):
#endif
       end module tensor_algebra_dil
