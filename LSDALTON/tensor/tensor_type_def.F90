
module tensor_type_def_module
  use,intrinsic :: iso_c_binding, only:c_ptr,c_null_ptr
  use tensor_parameters_and_counters

  !tile structure
  type tile
    type(c_ptr)           :: c    =  c_null_ptr
#ifdef VAR_PTR_RESHAPE
    real(tensor_real),pointer, contiguous :: t(:) => null()         !data in tiles
#else
    real(tensor_real),pointer     :: t(:) => null()         !data in tiles
#endif
    integer(kind=tensor_standard_int), pointer  :: d(:) => null()  !actual dimension of the tiles
    integer(kind=tensor_long_int)     :: e               !number of elements in current tile
    integer(kind=tensor_standard_int) :: gt              !global ti nr.
    integer(kind=tensor_mpi_kind)     :: wi              !window in pc group
  end type tile

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
     real(tensor_real), pointer :: elm1(:)               => null()       ! local chunk of memory

     ! the following should just point to elm1
     real(tensor_real), pointer :: elm2(:,:)             => null()
     real(tensor_real), pointer :: elm3(:,:,:)           => null()
     real(tensor_real), pointer :: elm4(:,:,:,:)         => null()
     real(tensor_real), pointer :: elm5(:,:,:,:,:)       => null()
     real(tensor_real), pointer :: elm6(:,:,:,:,:,:)     => null()
     real(tensor_real), pointer :: elm7(:,:,:,:,:,:,:)   => null()

     !in order to have only one array type the tile information is always there
     integer(kind=tensor_mpi_kind) :: dummyw
     type(c_ptr)                   :: dummyc   =  c_null_ptr
     real(tensor_real),pointer     :: dummy(:) => null()       !for the creation of mpi windows a dummy is required
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



end module tensor_type_def_module
