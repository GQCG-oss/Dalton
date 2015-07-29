module tensor_parameters_and_counters

   !Define atomic types used in the tensor module
   integer, parameter :: tensor_dp           = 8
   integer, parameter :: tensor_sp           = 4
#ifdef VAR_INT64
   integer, parameter :: tensor_int          = 8
   integer, parameter :: tensor_log          = 8
#else
   integer, parameter :: tensor_int          = 4
   integer, parameter :: tensor_log          = 4
#endif
   integer, parameter :: tensor_standard_int = 4
   integer, parameter :: tensor_standard_log = 4
   integer, parameter :: tensor_long_int     = 8
   integer, parameter :: tensor_long_log     = 8
   integer, parameter :: tensor_char_size    = 1
   integer, parameter :: tensor_char_4       = 4
!#ifdef VAR_MPI_32BIT_INT
   integer, parameter :: tensor_mpi_kind     = 4
!#else
!   integer, parameter :: tensor_mpi_kind     = 8
!#endif
   !
   !L2_CACHE_SIZE = 256000 ! lower estimate of cache size in bytes:
   !BS_2D = floor(sqrt(L2_CACHE_SIZE/(2*8.0E0_realk)))
   !BS_3D = floor((L2_CACHE_SIZE/(2*8.0E0_realk))**(1/3.0E0_realk))
   !BS_4D = floor((L2_CACHE_SIZE/(2*8.0E0_realk))**(1/4.0E0_realk))
   !
   !We write the explicit values to avoid internal compiler error 
   !with the gnu compiler 4.4.7 and the power function:
   integer,parameter :: BS_2D = 126
   integer,parameter :: BS_3D =  25
   integer,parameter :: BS_4D =  11

   !MPI SIGNAL, GET SLAVES, MAKE SURE THE APPLICATION HAS NO OVERLAPPING SIGNAL
   integer,parameter :: TENSOR_SLAVES_TO_SLAVE_ROUTINE_STD =  -121
   integer ::           TENSOR_SLAVES_TO_SLAVE_ROUTINE     = TENSOR_SLAVES_TO_SLAVE_ROUTINE_STD
   !MPI COMM TO USE IN TENSOR OPERATIONS, THIS IS UPDATED AT RUNTIME
   integer(kind=tensor_mpi_kind), parameter :: tensor_comm_null = -124
   integer(kind=tensor_mpi_kind) ::            tensor_work_comm = tensor_comm_null

   !parameters to define the data distribution in the tensor type
   integer(kind=tensor_standard_int), parameter :: TT_DENSE        = 1
   integer(kind=tensor_standard_int), parameter :: TT_REPLICATED   = 2
   integer(kind=tensor_standard_int), parameter :: TT_TILED        = 3
   integer(kind=tensor_standard_int), parameter :: TT_TILED_DIST   = 4
   integer(kind=tensor_standard_int), parameter :: TT_TILED_REPL   = 5
   
   !parameters for ACCESS TYPE:
   integer(kind=tensor_standard_int),parameter :: AT_NO_PDM_ACCESS = 0
   integer(kind=tensor_standard_int),parameter :: AT_MASTER_ACCESS = 1
   integer(kind=tensor_standard_int),parameter :: AT_ALL_ACCESS    = 2
   
   !other parameters
   integer,parameter :: TENSOR_MSG_LEN = 30
   integer,parameter :: DEFAULT_TDIM   = 10
   
   integer,parameter :: lspdm_stdout  = 6
   integer,parameter :: lspdm_errout  = 0
   
   !parameter to reduce the amount of windows per array
#ifdef VAR_HAVE_MPI3
   logical,parameter :: alloc_in_dummy = .true.
#else
   logical,parameter :: alloc_in_dummy = .false.
#endif
   !STANDARD MPI MESSAGE LENGTH SET TO ~10MB
   integer(kind=tensor_standard_int) :: TENSOR_MPI_MSG_LEN = 10000000
   
   !> execution time variables
   logical :: lspdm_use_comm_proc
   logical :: tensor_debug_mode           = .false.
   logical :: tensor_always_sync          = .false.
   logical :: tensor_contract_dil_backend = .false.
   logical :: tensor_segment_length_set   = .false.
   integer(kind=tensor_long_int) :: tensor_segment_length = -1
   
   !ALL COUNTERS IN BYTES!!
   !> Max allocated memory
   integer(kind=tensor_long_int) :: tensor_counter_max_hp_mem = 0_tensor_long_int
   integer(kind=tensor_long_int) :: tensor_counter_max_bg_mem = 0_tensor_long_int

   !> Allocated memory of dense array
   integer(kind=tensor_long_int) :: tensor_counter_dense_a_mem     = 0_tensor_long_int
   integer(kind=tensor_long_int) :: tensor_counter_dense_bg_a_mem  = 0_tensor_long_int
   !> Deallocated memory of dense array
   integer(kind=tensor_long_int) :: tensor_counter_dense_f_mem     = 0_tensor_long_int
   integer(kind=tensor_long_int) :: tensor_counter_dense_bg_f_mem  = 0_tensor_long_int
   !> Allocated memory of tiled array
   integer(kind=tensor_long_int) :: tensor_counter_tiled_a_mem     = 0_tensor_long_int
   integer(kind=tensor_long_int) :: tensor_counter_tiled_bg_a_mem  = 0_tensor_long_int
   !> Deallocated memory of tiled array
   integer(kind=tensor_long_int) :: tensor_counter_tiled_f_mem     = 0_tensor_long_int
   integer(kind=tensor_long_int) :: tensor_counter_tiled_bg_f_mem  = 0_tensor_long_int
   !> Allocated auxiliary memory of array
   integer(kind=tensor_long_int) :: tensor_counter_aux_a_mem       = 0_tensor_long_int
   integer(kind=tensor_long_int) :: tensor_counter_aux_bg_a_mem    = 0_tensor_long_int
   !> Deallocated auxiliary memory of array
   integer(kind=tensor_long_int) :: tensor_counter_aux_f_mem       = 0_tensor_long_int
   integer(kind=tensor_long_int) :: tensor_counter_aux_bg_f_mem    = 0_tensor_long_int
   !> Currently allocated memory on node
   integer(kind=tensor_long_int) :: tensor_counter_memory_in_use      = 0_tensor_long_int
   integer(kind=tensor_long_int) :: tensor_counter_memory_in_use_heap = 0_tensor_long_int
   integer(kind=tensor_long_int) :: tensor_counter_memory_in_use_bg   = 0_tensor_long_int

   integer(kind=tensor_long_int), pointer :: tensor_counter_ext_mem => null()

   
   
end module tensor_parameters_and_counters
