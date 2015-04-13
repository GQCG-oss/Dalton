
module tensor_type_def_module
  use,intrinsic :: iso_c_binding, only:c_ptr,c_null_ptr
  use precision

  !tile structure
  type tile
    type(c_ptr)           :: c    =  c_null_ptr
#ifdef VAR_PTR_RESHAPE
    real(realk),pointer, contiguous :: t(:) => null()         !data in tiles
#else
    real(realk),pointer   :: t(:) => null()         !data in tiles
#endif
    integer,pointer       :: d(:) => null()         !actual dimension of the tiles
    integer(kind=8)       :: e                      !number of elements in current tile
    integer               :: gt                     !global ti nr.
    integer(kind=ls_mpik) :: wi                     !window in pc group
  end type tile

  type tensor
     !mode=number of modes of the array or order of the corresponding tensor,
     !nelms=number of elements in the array
     !atype= format or distribution in which the array is stored --> dense, distributed --> see parameters in tensor_operations.f90
     integer :: mode
     integer(kind=8) :: nelms
     integer      :: itype                                              !integer type specification
     character(4) :: atype                                              !more precise character type specification

     !> Dimensions
     integer, pointer :: dims(:)                   => null ()

     !> Data, only allocate the first for the elements and use the others just
     !to reference the data in the first pointer
     integer(kind=ls_mpik):: w1                                    ! windows for local chunk of memory
     type(c_ptr)          :: e1c                   =  c_null_ptr   ! cpointer for local chunk, maily as fail-check
     real(realk), pointer :: elm1(:)               => null()       ! local chunk of memory

     ! the following should just point to elm1
     real(realk), pointer :: elm2(:,:)             => null()
     real(realk), pointer :: elm3(:,:,:)           => null()
     real(realk), pointer :: elm4(:,:,:,:)         => null()
     real(realk), pointer :: elm5(:,:,:,:,:)       => null()
     real(realk), pointer :: elm6(:,:,:,:,:,:)     => null()
     real(realk), pointer :: elm7(:,:,:,:,:,:,:)   => null()

     !in order to have only one array type the tile information is always there
     integer(kind=ls_mpik) :: dummyw
     type(c_ptr)           :: dummyc               =  c_null_ptr
     real(realk),pointer   :: dummy(:)             => null()       !for the creation of mpi windows a dummy is required
     type(tile),pointer    :: ti(:)                => null()       !tiles, if matrix should be distributed
     integer(kind=ls_mpik),pointer    :: wi(:)     => null()       !windows for tiles, if matrix should be distributed, there are ntiles windows to be inited
     integer               :: nwins                = 0             !number of windows allocated
     integer,pointer       :: ntpm(:)              => null()       !dimensions in the modes, number of tiles per mode, 
     integer,pointer       :: tdim(:)              => null()       !dimension of the tiles per mode(per def symmetric, but needed)
     integer,pointer       :: addr_p_arr(:)        => null()       !address of array in persistent array "p_arr" on each compute node
     logical,pointer       :: lock_set(:)          => null()
     integer               :: local_addr           = -1            !local address in p_arr
     !global tile information
     integer :: ntiles,tsize                         !batching of tiles in one mode, number of tiles, tilesize (ts^mode), amount of modes of the array
     integer :: nlti                                 !number of local tiles
     integer :: offset                               !use offset in nodes for the distribution of arrays
     integer(kind=ls_mpik) :: nnod,comm              !number of nodes and communicator
     integer,pointer :: access_type                  !type of access to the array
     logical :: zeros=.false.                        !use zeros in tiles --> it is at the moment not recommended to use .true. here
     !logical :: allocd_w_c_p                        ! allocated with comm_threads or not?
     logical :: initialized = .false.                !check variable if array is initialized
     logical :: bg_alloc

  end type tensor

  !> Allocated memory of dense array
  real(realk) :: tensor_dense_allocd_mem   = 0.0E0_realk
  !> Deallocated memory of dense array
  real(realk) :: tensor_dense_deallocd_mem = 0.0E0_realk
  !> Currently allocated memory on node
  real(realk) :: tensor_memory_in_use      = 0.0E0_realk
  !> Max allocated memory
  real(realk) :: tensor_max_memory         = 0.0E0_realk
  !> Allocated memory of tiled array
  real(realk) :: tensor_tiled_allocd_mem   = 0.0E0_realk
  !> Deallocated memory of tiled array
  real(realk) :: tensor_tiled_deallocd_mem = 0.0E0_realk
  !> Allocated auxiliary memory of array
  real(realk) :: tensor_aux_allocd_mem     = 0.0E0_realk
  !> Deallocated auxiliary memory of array
  real(realk) :: tensor_aux_deallocd_mem   = 0.0E0_realk


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



end module tensor_type_def_module
