
module tensor_type_def_module
  use,intrinsic :: iso_c_binding, only:c_ptr,c_null_ptr

  use lsmpi_module
  use tensor_parameters_and_counters
  use tensor_error_handler

  !tile structure
  ! -> always compare with the sizeof function
  integer, parameter :: tensor_bytes_per_tile = 1344
  type tile
    type(c_ptr)           :: c    =  c_null_ptr
#ifdef VAR_PTR_RESHAPE
    real(tensor_dp),pointer, contiguous :: t(:) => null()         !data in tiles
#else
    real(tensor_dp),pointer     :: t(:) => null()         !data in tiles
#endif
    integer(kind=tensor_int), pointer  :: d(:) => null()  !actual dimension of the tiles
    integer(kind=tensor_long_int)     :: e               !number of elements in current tile
    integer(kind=tensor_standard_int) :: gt              !global ti nr.
    integer(kind=tensor_mpi_kind)     :: wi              !window in pc group
  end type tile

  ! -> always compare with the sizeof function
  integer, parameter :: tensor_bytes_per_tensor = 120
  type tensor
     !mode=number of modes of the array or order of the corresponding tensor,
     !nelms=number of elements in the array
     !atype= format or distribution in which the array is stored --> dense, distributed --> see parameters in tensor_operations.f90
     integer(kind=tensor_standard_int)      :: mode
     integer(kind=tensor_long_int) :: nelms
     integer(kind=tensor_standard_int)      :: itype                                              !integer type specification
     character(tensor_char_4)      :: atype                                              !more precise character type specification

     !> Dimensions
     integer(kind=tensor_int), pointer :: dims(:)                   => null ()

     !> Data, only allocate the first for the elements and use the others just
     !to reference the data in the first pointer
     integer(kind=tensor_mpi_kind):: w1                                    ! windows for local chunk of memory
     type(c_ptr)                :: e1c                   =  c_null_ptr   ! cpointer for local chunk, maily as fail-check
     real(tensor_dp), pointer :: elm1(:)               => null()       ! local chunk of memory

     ! the following should just point to elm1
     real(tensor_dp), pointer :: elm2(:,:)             => null()
     real(tensor_dp), pointer :: elm3(:,:,:)           => null()
     real(tensor_dp), pointer :: elm4(:,:,:,:)         => null()
     real(tensor_dp), pointer :: elm5(:,:,:,:,:)       => null()
     real(tensor_dp), pointer :: elm6(:,:,:,:,:,:)     => null()
     real(tensor_dp), pointer :: elm7(:,:,:,:,:,:,:)   => null()

     !in order to have only one array type the tile information is always there
     integer(kind=tensor_mpi_kind) :: dummyw
     type(c_ptr)                   :: dummyc   =  c_null_ptr
     real(tensor_dp),pointer     :: dummy(:) => null()       !for the creation of mpi windows a dummy is required
     type(tile),pointer    :: ti(:)   => null()       !tiles, if matrix should be distributed
     integer(kind=tensor_mpi_kind),pointer  :: wi(:)     => null()       !windows for tiles, if matrix should be distributed, there are ntiles windows to be inited
     integer(kind=tensor_int)               :: nwins                = 0             !number of windows allocated
     integer(kind=tensor_standard_int),pointer       :: ntpm(:)              => null()       !dimensions in the modes, number of tiles per mode, 
     integer(kind=tensor_standard_int),pointer       :: tdim(:)              => null()       !dimension of the tiles per mode(per def symmetric, but needed)
     integer(kind=tensor_standard_int),pointer       :: addr_p_arr(:)        => null()       !address of array in persistent array "p_arr" on each compute node
     logical(kind=tensor_standard_log),pointer       :: lock_set(:)          => null()
     integer(kind=tensor_standard_int)               :: local_addr           = -1            !local address in p_arr
     !global tile information
     integer(kind=tensor_standard_int) :: ntiles,tsize                         !batching of tiles in one mode, number of tiles, tilesize (ts^mode), amount of modes of the array
     integer(kind=tensor_standard_int) :: nlti                                 !number of local tiles
     integer(kind=tensor_standard_int) :: offset                               !use offset in nodes for the distribution of arrays
     integer(kind=tensor_mpi_kind) :: nnod,comm              !number of nodes and communicator
     integer(kind=tensor_standard_int),pointer :: access_type                  !type of access to the array
     logical(kind=tensor_standard_log) :: zeros=.false.                        !use zeros in tiles --> it is at the moment not recommended to use .true. here
     !logical :: allocd_w_c_p                        ! allocated with comm_threads or not?
     logical(kind=tensor_standard_log) :: initialized = .false.                !check variable if array is initialized
     logical(kind=tensor_standard_log) :: bg_alloc

  end type tensor

  !Type counter
  type tensor_counter_type
     integer(kind=tensor_long_int) :: size_ = 0_tensor_long_int
     integer(kind=tensor_long_int) :: curr_ = 0_tensor_long_int
     integer(kind=tensor_long_int) :: high_ = 0_tensor_long_int
  end type tensor_counter_type

  !make sure the following types all have the same structure:
  integer(kind=tensor_standard_int), parameter :: tensor_mem_idx_tensor_dp           = 1
  integer(kind=tensor_standard_int), parameter :: tensor_mem_idx_tensor_dp_mpi       = 2
  integer(kind=tensor_standard_int), parameter :: tensor_mem_idx_tensor_standard_int = 3
  integer(kind=tensor_standard_int), parameter :: tensor_mem_idx_tensor_long_int     = 4
  integer(kind=tensor_standard_int), parameter :: tensor_mem_idx_tensor_standard_log = 5
  integer(kind=tensor_standard_int), parameter :: tensor_mem_idx_tensor_long_log     = 6
  integer(kind=tensor_standard_int), parameter :: tensor_mem_idx_tile                = 7
  integer(kind=tensor_standard_int), parameter :: tensor_mem_idx_tensor              = 8

  !MAKE SURE THIS INDEX CORRESPONDS TO THE HIGHEST COUNTER IN THE UPPER LIST
  integer(kind=tensor_standard_int), parameter :: tensor_nmem = 8
  type(tensor_counter_type) :: counters(tensor_nmem)
  type(tensor_counter_type) :: counters_bg(tensor_nmem)


  contains

  subroutine tensor_init_counters()
     implicit none
     call tensor_init_specific_counter(counters,tensor_nmem)
     call tensor_init_specific_counter(counters_bg,tensor_nmem)
  end subroutine tensor_init_counters

  subroutine tensor_init_specific_counter(c,l)
     implicit none
     integer(kind=tensor_standard_int), intent(in) :: l
     type(tensor_counter_type) :: c(l)
     integer(kind=tensor_standard_int) :: i
     type(tile)   :: ref_tile
     type(tensor) :: ref_tensor
#ifdef VAR_MPI
     integer(kind=MPI_ADDRESS_KIND) :: sze,lb
     integer(kind=tensor_mpi_kind)  :: ierr
#endif
     do i=1,l

        !DEFINE THE SIZE OF ONE UNIT
        select case(i)
        !ATOMIC TYPES
        case(tensor_mem_idx_tensor_dp)
           c(i)%size_ = tensor_dp
        case(tensor_mem_idx_tensor_dp_mpi)
#ifdef VAR_MPI
           call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,sze,ierr)
           c(i)%size_ = sze
#else
           c(i)%size_ = 0
#endif
        case(tensor_mem_idx_tensor_standard_int)
           c(i)%size_ = tensor_standard_int
        case(tensor_mem_idx_tensor_long_int)
           c(i)%size_ = tensor_long_int
        case(tensor_mem_idx_tensor_standard_log)
           c(i)%size_ = tensor_standard_log
        case(tensor_mem_idx_tensor_long_log)
           c(i)%size_ = tensor_long_log

        !DERIVED TYPES
        case(tensor_mem_idx_tile)
#ifdef VAR_SIZEOF_DEFINED
           c(i)%size_ = sizeof(ref_tile)
#else
           c(i)%size_ = tensor_bytes_per_tile
#endif
        case(tensor_mem_idx_tensor)
#ifdef VAR_SIZEOF_DEFINED
           c(i)%size_ = sizeof(ref_tensor)
#else
           c(i)%size_ = tensor_bytes_per_tensor
#endif
        case default 
           call tensor_status_quit("ERROR(tensor_init_basix): wrong index in&
           & setting up the memory counters. This is a coding error in&
           & tensor_type_def_module",202)
        end select

        !SET COUNTERS TO 0
        c(i)%curr_ = 0_tensor_long_int
        c(i)%high_ = 0_tensor_long_int

     enddo
  end subroutine tensor_init_specific_counter

  
  subroutine tensor_free_counters()
     implicit none
     integer(kind=tensor_standard_int) :: i

     write (*,*) "Tensor memory finalization:"
     write (*,*) "---------------------------"

     do i=1,tensor_nmem

        write (*,'(a)',advance='no') "Currently allocated bytes of type"
        call write_string_for_type(i)
        write (*,*) counters(i)%curr_
        write (*,'(a)',advance='no') "Max       allocated bytes of type"
        call write_string_for_type(i)
        write (*,*) counters(i)%high_
        write (*,*) "--------------------------------------------------------------------------"

     enddo
  end subroutine tensor_free_counters

  subroutine write_string_for_type(i)
     implicit none
     integer(kind=tensor_standard_int) :: i
     !DEFINE THE SIZE OF ONE UNIT
     select case(i)
     !ATOMIC TYPES
     case(tensor_mem_idx_tensor_dp)
        write (*,'(a)',advance='no') " tensor_dp           "
     case(tensor_mem_idx_tensor_dp_mpi)
        write (*,'(a)',advance='no') " tensor_dp_mpi       "
     case(tensor_mem_idx_tensor_standard_int)
        write (*,'(a)',advance='no') " tensor_standard_int "
     case(tensor_mem_idx_tensor_long_int)
        write (*,'(a)',advance='no') " tensor_long_int     "
     case(tensor_mem_idx_tensor_standard_log)
        write (*,'(a)',advance='no') " tensor_standard_log "
     case(tensor_mem_idx_tensor_long_log)
        write (*,'(a)',advance='no') " tensor_long_log     "

     !DERIVED TYPES
     case(tensor_mem_idx_tile)
        write (*,'(a)',advance='no') " tile                "
     case(tensor_mem_idx_tensor)
        write (*,'(a)',advance='no') " tensor              "
     case default 
        call tensor_status_quit("ERROR(tensor_init_basix): wrong index in&
        & setting up the memory counters. This is a coding error in&
        & tensor_type_def_module",202)
     end select
  end subroutine write_string_for_type

end module tensor_type_def_module
