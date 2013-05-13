
module tensor_type_def_module
  use precision
  use,intrinsic :: iso_c_binding, only:c_ptr

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

  interface get_midx
    module procedure get_mode_idx8,&
                   & get_mode_idx4
  end interface get_midx

  interface get_cidx
    module procedure get_comp_idx
  end interface get_cidx

  interface get_tile_dim
    module procedure get_tileinfo_nels_frombas,&
                    &get_tileinfo_nelspmode_frombas,&
                    &get_tileinfo_nels_fromarr8,&
                    &get_tileinfo_nels_fromarr4,&
                    &get_tileinfo_nelspermode_fromarr4,&
                    &get_tileinfo_nelspermode_fromarr8
  end interface get_tile_dim
  contains
  !> \author Patrick Ettenhuber
  !> \date June 2012
  !> \brief get mode index from composite index
  subroutine get_mode_idx8(a,inds,dims,modes)
    implicit none
    integer(kind=8),intent(in):: a
    integer,intent(in)        :: dims(*),modes
    integer,intent(inout)     :: inds(*)
    integer(kind=8)           :: i,cind,ndim
    select case(modes)
    case default
      cind=a
      do i=1,modes
        ndim = dims(i)
        inds(i)=mod(cind-1,ndim)+1
        cind=(cind-inds(i))/dims(i) + 1
      enddo
    end select
  end subroutine get_mode_idx8
  !> \author Patrick Ettenhuber
  !> \date June 2012
  !> \brief get mode index from composite index
  subroutine get_mode_idx4(a,inds,dims,modes)
    implicit none
    integer(kind=4),intent(in):: a
    integer,intent(in)        :: dims(*),modes
    integer,intent(inout)     :: inds(*)
    integer(kind=4)           :: i,cind
    select case(modes)
    case default
      cind=a
      do i=1,modes
        inds(i)=mod(cind-1,dims(i))+1
        cind=(cind-inds(i))/dims(i) + 1
      enddo
    end select
  end subroutine get_mode_idx4

  !> \author Patrick Ettenhuber
  !> \date June 2012
  !> \brief get composite index from mode index
  function get_comp_idx(inds,dims,modes) result(a)
    implicit none
    integer,intent(in) :: inds(*),dims(*),modes
    integer :: i,j,cdim
    integer :: a
    select case(modes)
   ! case(4)
   !   a=inds(1)+(inds(2)-1)*dims(1)+(inds(3)-1)*dims(1)*dims(2)+(inds(4)-1)*inds(1)*inds(2)*inds(3)
    case default
      a=1
      do i=1,modes
        cdim=1
        do j=i-1,1,-1
          cdim=cdim*dims(j)
        enddo
        a=a+(inds(i)-1)*cdim
      enddo
    end select
  end function get_comp_idx

  subroutine get_tileinfo_nels_frombas(sze,tileidx,dims,tdim,mode,offset)
    implicit none
    integer, intent(out) :: sze
    integer,intent(in) :: tileidx,mode,dims(mode),tdim(mode)
    integer,intent(in),optional :: offset
    integer :: j,orig_addr(mode),offs,ntpm(mode)
    offs=1
    if(present(offset))offs=offset
    do j=1,mode
      ntpm(j)=dims(j)/tdim(j)
      if(mod(dims(j),tdim(j))>0)ntpm(j)=ntpm(j)+1
    enddo
    call get_midx(tileidx,orig_addr,ntpm,mode)
    sze=1
    do j=offs, mode
      if(((dims(j)-(orig_addr(j)-1)*tdim(j))/tdim(j))>=1)then
        sze=sze*tdim(j)
      else
        sze=sze*mod(dims(j),tdim(j))
      endif
    enddo
  end subroutine get_tileinfo_nels_frombas
  !> \brief this function returns the number of elements of a tile where the tile index is
  ! a global tile index
  !> \author Patrick Ettenhuber
  subroutine get_tileinfo_nels_fromarr8(nels,arr,tnumber)
    implicit none
    !> array for which nels shoulb be calculated
    type(array),intent(in) :: arr
    !> global tile index for which nels should be calculated
    integer(kind=long), intent(in) :: tnumber
    !> return value, number of elements in the desired tile
    integer :: nels
    integer ::orig_addr(arr%mode),j
    call get_midx(tnumber,orig_addr,arr%ntpm,arr%mode)
    nels=1
    do j=1, arr%mode
      if(((arr%dims(j)-(orig_addr(j)-1)*arr%tdim(j))/arr%tdim(j))>=1)then
        nels=nels*arr%tdim(j)
      else
        nels=nels*mod(arr%dims(j),arr%tdim(j))
      endif
    enddo
  end subroutine get_tileinfo_nels_fromarr8
  subroutine get_tileinfo_nels_fromarr4(nels,arr,tnumber)
    implicit none
    !> array for which nels shoulb be calculated
    type(array),intent(in) :: arr
    !> global tile index for which nels should be calculated
    integer(kind=4), intent(in) :: tnumber
    !> return value, number of elements in the desired tile
    integer :: nels
    integer ::orig_addr(arr%mode),j
    call get_midx(tnumber,orig_addr,arr%ntpm,arr%mode)
    nels=1
    do j=1, arr%mode
      if(((arr%dims(j)-(orig_addr(j)-1)*arr%tdim(j))/arr%tdim(j))>=1)then
        nels=nels*arr%tdim(j)
      else
        nels=nels*mod(arr%dims(j),arr%tdim(j))
      endif
    enddo
  end subroutine get_tileinfo_nels_fromarr4
  subroutine get_tileinfo_nelspermode_fromarr8(nels,arr,tnumber)
    implicit none
    !> array for which nels shoulb be calculated
    type(array),intent(in) :: arr
    !> global tile index for which nels should be calculated
    integer(kind=long), intent(in) :: tnumber
    !> return value, number of elements in the desired tile
    integer :: nels(arr%mode)
    integer ::orig_addr(arr%mode),j
    call get_midx(tnumber,orig_addr,arr%ntpm,arr%mode)
    nels=1
    do j=1, arr%mode
      if(((arr%dims(j)-(orig_addr(j)-1)*arr%tdim(j))/arr%tdim(j))>=1)then
        nels(j)=arr%tdim(j)
      else
        nels(j)=mod(arr%dims(j),arr%tdim(j))
      endif
    enddo
  end subroutine get_tileinfo_nelspermode_fromarr8
  subroutine get_tileinfo_nelspermode_fromarr4(nels,arr,tnumber)
    implicit none
    !> array for which nels shoulb be calculated
    type(array),intent(in) :: arr
    !> global tile index for which nels should be calculated
    integer(kind=4), intent(in) :: tnumber
    !> return value, number of elements in the desired tile
    integer :: nels(arr%mode)
    integer ::orig_addr(arr%mode),j
    call get_midx(tnumber,orig_addr,arr%ntpm,arr%mode)
    nels=1
    do j=1, arr%mode
      if(((arr%dims(j)-(orig_addr(j)-1)*arr%tdim(j))/arr%tdim(j))>=1)then
        nels(j)=arr%tdim(j)
      else
        nels(j)=mod(arr%dims(j),arr%tdim(j))
      endif
    enddo
  end subroutine get_tileinfo_nelspermode_fromarr4
  subroutine get_tileinfo_nelspmode_frombas(sze,tileidx,dims,tdim,mode)
    implicit none
    integer,intent(in) :: tileidx,mode,dims(mode),tdim(mode)
    integer, dimension(mode),intent(out) :: sze
    integer :: j,orig_addr(mode),ntpm(mode)
    do j=1,mode
      ntpm(j)=dims(j)/tdim(j)
      if(mod(dims(j),tdim(j))>0)ntpm(j)=ntpm(j)+1
    enddo
    call get_midx(tileidx,orig_addr,ntpm,mode)
    do j=1, mode
      if(((dims(j)-(orig_addr(j)-1)*tdim(j))/tdim(j))>=1)then
        sze(j)=tdim(j)
      else
        sze(j)=mod(dims(j),tdim(j))
      endif
    enddo
  end subroutine get_tileinfo_nelspmode_frombas

end module tensor_type_def_module
