module tensor_parameters_and_counters
   use precision

   !Define atomic types used in the tensor module
   integer, parameter :: tensor_real     = 8
#ifdef VAR_INT64
   integer, parameter :: tensor_int      = 8
#else
   integer, parameter :: tensor_int      = 4
#endif
   integer, parameter :: tensor_standard_int = 4
   integer, parameter :: tensor_standard_log = 4
   integer, parameter :: tensor_long_int = 8
   integer, parameter :: tensor_char_4   = 4
!#ifdef VAR_MPI_32BIT_INT
   integer, parameter :: tensor_mpi_kind = 4
!#else
!   integer, parameter :: tensor_mpi_kind = 8
!#endif

   !parameters to define the data distribution in the tensor type
   integer, parameter :: TT_DENSE        = 1
   integer, parameter :: TT_REPLICATED   = 2
   integer, parameter :: TT_TILED        = 3
   integer, parameter :: TT_TILED_DIST   = 4
   integer, parameter :: TT_TILED_REPL   = 5
   
   !parameters for ACCESS TYPE:
   integer,parameter :: AT_NO_PDM_ACCESS = 0
   integer,parameter :: AT_MASTER_ACCESS = 1
   integer,parameter :: AT_ALL_ACCESS    = 2
   
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
   
   !> execution time variables
   logical :: lspdm_use_comm_proc
   logical :: tensor_debug_mode           = .false.
   logical :: tensor_always_sync          = .false.
   logical :: tensor_contract_dil_backend = .false.
   logical :: tensor_segment_length_set   = .false.
   integer(kind=tensor_long_int) :: tensor_segment_length = -1
   
   !ALL COUNTERS IN BYTES!!

   !size of a tile in bytes, this needs to be updated whenever the tile structure is changed
   ! right now, we assume 8 bytes for each member of the struct
   integer, parameter :: tensor_bytes_per_tile = 7*8

   !Counters for total memory usage, a = allocd, f = freed, bg = background
   real(tensor_long_int) :: tensor_counter_total_a_mem = 0_tensor_long_int
   real(tensor_long_int) :: tensor_counter_total_f_mem = 0_tensor_long_int
   real(tensor_long_int) :: tensor_counter_stack_a_mem = 0_tensor_long_int
   real(tensor_long_int) :: tensor_counter_stack_f_mem = 0_tensor_long_int
   real(tensor_long_int) :: tensor_counter_bgmem_a_mem = 0_tensor_long_int
   real(tensor_long_int) :: tensor_counter_bgmem_f_mem = 0_tensor_long_int
   !> Max allocated memory
   real(tensor_long_int) :: tensor_counter_max_memory         = 0_tensor_long_int
   !counters for the individual ATOMIC types which may be allocated
   real(tensor_long_int) :: tensor_counter_tensor_real_a_mem  = 0_tensor_long_int
   real(tensor_long_int) :: tensor_counter_tensor_real_f_mem  = 0_tensor_long_int
   real(tensor_long_int) :: tensor_counter_tensor_real_bg_a_mem = 0_tensor_long_int
   real(tensor_long_int) :: tensor_counter_tensor_real_bg_f_mem = 0_tensor_long_int

   !> Allocated memory of dense array
   real(tensor_long_int) :: tensor_counter_dense_a_mem     = 0_tensor_long_int
   real(tensor_long_int) :: tensor_counter_dense_bg_a_mem  = 0_tensor_long_int
   !> Deallocated memory of dense array
   real(tensor_long_int) :: tensor_counter_dense_f_mem     = 0_tensor_long_int
   real(tensor_long_int) :: tensor_counter_dense_bg_f_mem  = 0_tensor_long_int
   !> Allocated memory of tiled array
   real(tensor_long_int) :: tensor_counter_tiled_a_mem     = 0_tensor_long_int
   real(tensor_long_int) :: tensor_counter_tiled_bg_a_mem  = 0_tensor_long_int
   !> Deallocated memory of tiled array
   real(tensor_long_int) :: tensor_counter_tiled_f_mem     = 0_tensor_long_int
   real(tensor_long_int) :: tensor_counter_tiled_bg_f_mem  = 0_tensor_long_int
   !> Allocated auxiliary memory of array
   real(tensor_long_int) :: tensor_counter_aux_a_mem       = 0_tensor_long_int
   real(tensor_long_int) :: tensor_counter_aux_bg_a_mem    = 0_tensor_long_int
   !> Deallocated auxiliary memory of array
   real(tensor_long_int) :: tensor_counter_aux_f_mem   = 0_tensor_long_int
   real(tensor_long_int) :: tensor_counter_aux_bg_f_mem    = 0_tensor_long_int
   !> Currently allocated memory on node
   real(tensor_long_int) :: tensor_counter_memory_in_use      = 0_tensor_long_int
   real(tensor_long_int) :: tensor_counter_memory_in_use_heap = 0_tensor_long_int
   real(tensor_long_int) :: tensor_counter_memory_in_use_bg   = 0_tensor_long_int
   
   
end module tensor_parameters_and_counters
