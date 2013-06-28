
module tensor_type_def_module
  use,intrinsic :: iso_c_binding, only:c_ptr
  use precision

  !tile structure

  type tile
    type(c_ptr) :: c
    real(realk),pointer :: t(:) => null()         !data in tiles
    integer,pointer :: d(:)     => null()         !actual dimension of the tiles
    integer :: e,gt                               !number of elements in current tile, global ti nr
  end type tile

  type array
     !mode=number of modes of the array or order of the corresponding tensor,
     !nelms=number of elements in the array
     !atype= format or distribution in which the array is stored --> dense, distributed --> see parameters in array_operations.f90
     integer :: mode
     integer(kind=8) :: nelms
     integer :: atype
     !> Dimensions
     integer, pointer :: dims(:)     => null ()
     !> Data, only allocate the first for the elements and use the others just
     !to reference the data in the first pointer
     real(realk), pointer :: elm1(:) => null()
     ! the following should just point to elm1
     real(realk), pointer :: elm2(:,:) => null()
     real(realk), pointer :: elm3(:,:,:) => null()
     real(realk), pointer :: elm4(:,:,:,:) => null()
     real(realk), pointer :: elm5(:,:,:,:,:) => null()
     real(realk), pointer :: elm6(:,:,:,:,:,:) => null()
     real(realk), pointer :: elm7(:,:,:,:,:,:,:) => null()


     !in order to have only one array type the tile information is always there
     type(c_ptr)        :: dummyc
     real(realk),pointer:: dummy(:)  => null()       !for the creation of mpi windows a dummy is required
     type(tile),pointer :: ti(:)     => null()       !tiles, if matrix should be distributed
     integer(kind=ls_mpik),pointer    :: wi(:)     => null()       !windows for tiles, if matrix should be distributed, there are ntiles windows to be inited
     integer,pointer    :: ntpm(:)   => null()       !dimensions in the modes, number of tiles per mode, 
     integer,pointer    :: tdim(:)   => null()       !dimension of the tiles per mode(per def symmetric, but needed)
     integer,pointer    :: addr_p_arr(:)   => null() !address of array in persistent array "p_arr" on each node
     logical, pointer   :: lock_set(:)     => null() !contains information about whether a lock has been set on the window with the current index
     !global tile information
     integer :: ntiles,tsize                         !batching of tiles in one mode, number of tiles, tilesize (ts^mode), amount of modes of the array
     integer :: nlti                                 !number of local tiles
     integer :: offset                               !use offset in nodes for the distribution of arrays
     integer :: init_type                            !type of initializtation
     logical :: zeros=.false.                        !use zeros in tiles --> it is at the moment not recommended to use .true. here

  end type array

  !> Allocated memory of dense array
  real(realk) :: array_dense_allocd_mem = 0.0E0_realk
  !> Deallocated memory of dense array
  real(realk) :: array_dense_deallocd_mem = 0.0E0_realk
  !> Currently allocated memory on node
  real(realk) :: array_memory_in_use = 0.0E0_realk
  !> Max allocated memory
  real(realk) :: array_max_memory = 0.0E0_realk
  !> Allocated memory of tiled array
  real(realk) :: array_tiled_allocd_mem = 0.0E0_realk
  !> Deallocated memory of tiled array
  real(realk) :: array_tiled_deallocd_mem = 0.0E0_realk
  !> Allocated auxiliary memory of array
  real(realk) :: array_aux_allocd_mem = 0.0E0_realk
  !> Deallocated auxiliary memory of array
  real(realk) :: array_aux_deallocd_mem = 0.0E0_realk


  !parameters to define the data distribution in the array type
  integer, parameter :: DENSE=1
  integer, parameter :: REPLICATED=2
  integer, parameter :: TILED=3
  integer, parameter :: TILED_DIST=4

  !parameters for PDMTYPE:
  integer,parameter :: NO_PDM=0
  integer,parameter :: MASTER_INIT=1
  integer,parameter :: ALL_INIT=2

  !other parameters
  integer,parameter :: ARR_MSG_LEN=30
  integer,parameter :: DEFAULT_TDIM=10
  
  integer,parameter :: lspdm_stdout=6
  integer,parameter :: lspdm_errout=0




end module tensor_type_def_module
