!> @file
!> here the one sided wrappers should go, so that interfacing with
!different rma vendors is possible
!> \author Patrick Ettenhuber
!> \date April 2013
module lspdm_basic_module

  use background_buffer_module, only: mem_is_background_buf_init
  use memory_handling, only: mem_pseudo_alloc, mem_pseudo_dealloc

  use tensor_allocator
  use tensor_type_def_module
  use LSTIMING!,only:lstimer
#ifdef VAR_MPI
  use infpar_module
  use lsmpi_type
#endif
  use reorder_frontend_module

  interface get_tile_dim
    module procedure get_tileinfo_nels_frombas,&
                    &get_tileinfo_nelspmode_frombas,&
                    &get_tileinfo_nels_fromarr88,&
                    &get_tileinfo_nels_fromarr84,&
                    &get_tileinfo_nels_fromarr48,&
                    &get_tileinfo_nels_fromarr44,&
                    &get_tileinfo_nels4_fromarr8mode,&
                    &get_tileinfo_nels8_fromarr8mode,&
                    &get_tileinfo_nels4_fromarr4mode,&
                    &get_tileinfo_nels8_fromarr4mode,&
                    &get_tileinfo_nelspermode_fromarr4mode,&
                    &get_tileinfo_nelspermode_fromarr8mode,&
                    &get_tileinfo_nelspermode_fromarr4,&
                    &get_tileinfo_nelspermode_fromarr8
  end interface get_tile_dim

  interface get_residence_of_tile
     module procedure get_residence_of_tile44,&
                     &get_residence_of_tile48,&
                     &get_residence_of_tile84,&
                     &get_residence_of_tile88
  end interface get_residence_of_tile

  interface tensor_get_ntpm
     module procedure tensor_get_ntpm8888,&
                     &tensor_get_ntpm8488,&
                     &tensor_get_ntpm4444
  end interface tensor_get_ntpm

  !interface get_tile_idx
  !  module procedure get_tile_idx_from_global_idx
  !end interface get_tile_idx

  contains
  !subroutine get_tile_idx_from_global_idx()
  !  implicit none
  !end subroutine get_tile_idx_from_global_idx



  subroutine get_tileinfo_nels_frombas(sze,tileidx,dims,tdim,mode)
     implicit none
     integer, intent(out) :: sze
     integer,intent(in) :: tileidx,mode,dims(mode),tdim(mode)
     !integer,intent(in),optional :: offset
     integer :: j,orig_addr(mode),offs,ntpm(mode)
     integer(kind=tensor_int) :: dims_long(mode)
     integer(kind=tensor_standard_int) :: tdim_std(mode), ntpm_std(mode)
     dims_long = dims
     tdim_std  = tdim
     offs=1
     !if(present(offset))offs=offset
     call tensor_get_ntpm(dims_long,tdim_std,int(mode,kind=tensor_standard_int),ntpm_std)
     ntpm = ntpm_std
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
  subroutine get_tileinfo_nelspmode_frombas(sze,tileidx,dims,tdim,mode)
     implicit none
     integer, intent(in)  :: tileidx,mode,dims(mode),tdim(mode)
     integer, intent(out) :: sze(mode)
     integer :: j,orig_addr(mode),ntpm(mode)
     integer(kind=tensor_int) :: dims_long(mode)
     integer(kind=tensor_standard_int) :: tdim_std(mode),ntpm_std(mode)
     dims_long = dims
     tdim_std  = tdim
     call tensor_get_ntpm(dims_long,tdim_std,int(mode,kind=tensor_standard_int),ntpm_std)
     ntpm = ntpm_std
     call get_midx(tileidx,orig_addr,ntpm,mode)
     do j=1, mode
        if(((dims(j)-(orig_addr(j)-1)*tdim(j))/tdim(j))>=1)then
           sze(j)=tdim(j)
        else
           sze(j)=mod(dims(j),tdim(j))
        endif
     enddo
  end subroutine get_tileinfo_nelspmode_frombas

  !> \author Patrick Ettenhuber
  subroutine get_tileinfo_nels4_fromarr8mode(nels,arr,orig_addr)
     implicit none
     !> array for which nels shoulb be calculated
     type(tensor),intent(in) :: arr
     !> global mode index of the tile
     integer(kind=tensor_long_int), intent(in) :: orig_addr(arr%mode)
     !> return value, number of elements in the desired tile
     integer(kind=tensor_standard_int) :: nels
     integer ::j
     include "get_nels.inc"
  end subroutine get_tileinfo_nels4_fromarr8mode
  subroutine get_tileinfo_nels8_fromarr8mode(nels,arr,orig_addr)
     implicit none
     !> array for which nels shoulb be calculated
     type(tensor),intent(in) :: arr
     !> global mode index of the tile
     integer(kind=tensor_long_int), intent(in) :: orig_addr(arr%mode)
     !> return value, number of elements in the desired tile
     integer(kind=tensor_long_int) :: nels
     integer ::j
     include "get_nels.inc"
  end subroutine get_tileinfo_nels8_fromarr8mode
  !> \author Patrick Ettenhuber
  subroutine get_tileinfo_nels8_fromarr4mode(nels,arr,orig_addr)
     implicit none
     !> array for which nels shoulb be calculated
     type(tensor),intent(in) :: arr
     !> global mode index of the tile
     integer(kind=tensor_standard_int), intent(in) :: orig_addr(arr%mode)
     !> return value, number of elements in the desired tile
     integer(kind=tensor_long_int) :: nels
     integer ::j
     include "get_nels.inc"
  end subroutine get_tileinfo_nels8_fromarr4mode
  subroutine get_tileinfo_nels4_fromarr4mode(nels,arr,orig_addr)
     implicit none
     !> array for which nels shoulb be calculated
     type(tensor),intent(in) :: arr
     !> global mode index of the tile
     integer(kind=tensor_standard_int), intent(in) :: orig_addr(arr%mode)
     !> return value, number of elements in the desired tile
     integer(kind=tensor_standard_int) :: nels
     integer ::j
     include "get_nels.inc"
  end subroutine get_tileinfo_nels4_fromarr4mode

  !> \brief this function returns the number of elements of a tile where the tile index is
  ! a global tile index
  !> \author Patrick Ettenhuber
  subroutine get_tileinfo_nels_fromarr88(nels,arr,tnumber)
     implicit none
     !> array for which nels shoulb be calculated
     type(tensor),intent(in) :: arr
     !> global tile index for which nels should be calculated
     integer(kind=tensor_long_int), intent(in) :: tnumber
     !> return value, number of elements in the desired tile
     integer(kind=tensor_long_int) :: nels
     integer ::orig_addr(arr%mode),j
     call get_midx(tnumber,orig_addr,arr%ntpm,arr%mode)
     include "get_nels.inc"
  end subroutine get_tileinfo_nels_fromarr88
  subroutine get_tileinfo_nels_fromarr84(nels,arr,tnumber)
     implicit none
     !> array for which nels shoulb be calculated
     type(tensor),intent(in) :: arr
     !> global tile index for which nels should be calculated
     integer(kind=tensor_standard_int), intent(in) :: tnumber
     !> return value, number of elements in the desired tile
     integer(kind=tensor_long_int) :: nels
     integer ::orig_addr(arr%mode),j
     call get_midx(tnumber,orig_addr,arr%ntpm,arr%mode)
     include "get_nels.inc"
  end subroutine get_tileinfo_nels_fromarr84
  subroutine get_tileinfo_nels_fromarr48(nels,arr,tnumber)
     implicit none
     !> array for which nels shoulb be calculated
     type(tensor),intent(in) :: arr
     !> global tile index for which nels should be calculated
     integer(kind=tensor_long_int), intent(in) :: tnumber
     !> return value, number of elements in the desired tile
     integer(kind=tensor_standard_int) :: nels
     integer ::orig_addr(arr%mode),j
     call get_midx(tnumber,orig_addr,arr%ntpm,arr%mode)
     include "get_nels.inc"
  end subroutine get_tileinfo_nels_fromarr48
  subroutine get_tileinfo_nels_fromarr44(nels,arr,tnumber)
     implicit none
     !> array for which nels shoulb be calculated
     type(tensor),intent(in) :: arr
     !> global tile index for which nels should be calculated
     integer(kind=tensor_standard_int), intent(in) :: tnumber
     !> return value, number of elements in the desired tile
     integer(kind=tensor_standard_int) :: nels
     integer ::orig_addr(arr%mode),j
     call get_midx(tnumber,orig_addr,arr%ntpm,arr%mode)
     include "get_nels.inc"
  end subroutine get_tileinfo_nels_fromarr44


  subroutine get_tileinfo_nelspermode_fromarr8(nels,arr,tnumber)
     implicit none
     !> array for which nels shoulb be calculated
     type(tensor),intent(in) :: arr
     !> global tile index for which nels should be calculated
     integer(kind=tensor_long_int), intent(in) :: tnumber
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

  subroutine get_tileinfo_nelspermode_fromarr4mode(nels,arr,orig_addr)
     implicit none
     !> array for which nels shoulb be calculated
     type(tensor),intent(in) :: arr
     !> global tile index for which nels should be calculated
     integer(kind=tensor_standard_int), intent(in) :: orig_addr(arr%mode)
     !> return value, number of elements in the desired tile
     integer :: nels(arr%mode)
     integer ::j
     nels=1
     do j=1, arr%mode
        if(((arr%dims(j)-(orig_addr(j)-1)*arr%tdim(j))/arr%tdim(j))>=1)then
           nels(j)=arr%tdim(j)
        else
           nels(j)=mod(arr%dims(j),arr%tdim(j))
        endif
     enddo
  end subroutine get_tileinfo_nelspermode_fromarr4mode
  subroutine get_tileinfo_nelspermode_fromarr8mode(nels,arr,orig_addr)
     implicit none
     !> array for which nels shoulb be calculated
     type(tensor),intent(in) :: arr
     !> global tile index for which nels should be calculated
     integer(kind=long), intent(in) :: orig_addr(arr%mode)
     !> return value, number of elements in the desired tile
     integer :: nels(arr%mode)
     integer ::j
     nels=1
     do j=1, arr%mode
        if(((arr%dims(j)-(orig_addr(j)-1)*arr%tdim(j))/arr%tdim(j))>=1)then
           nels(j)=arr%tdim(j)
        else
           nels(j)=mod(arr%dims(j),arr%tdim(j))
        endif
     enddo
  end subroutine get_tileinfo_nelspermode_fromarr8mode
  subroutine get_tileinfo_nelspermode_fromarr4(nels,arr,tnumber)
     implicit none
     !> array for which nels shoulb be calculated
     type(tensor),intent(in) :: arr
     !> global tile index for which nels should be calculated
     integer(kind=tensor_standard_int), intent(in) :: tnumber
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


  !> \brief calculate the number of tiles per mode
  !> \author Patrick Ettenhuber
  !> \date march 2013
  subroutine tensor_get_ntpm8888(dims,tdim,mode,ntpm,ntiles)
     implicit none
     integer(kind=tensor_long_int), intent(in)  :: mode
     integer(kind=tensor_long_int), intent(out) :: ntpm(mode)
     integer(kind=tensor_int),     intent(in)  :: dims(mode)
     integer(kind=tensor_long_int), intent(in)  :: tdim(mode)
     integer(kind=tensor_long_int), optional, intent(out) :: ntiles
     include "get_ntpm.inc"
  end subroutine tensor_get_ntpm8888
  subroutine tensor_get_ntpm8488(dims,tdim,mode,ntpm,ntiles)
     implicit none
     integer(kind=tensor_long_int), intent(in)  :: mode
     integer(kind=tensor_long_int), intent(out) :: ntpm(mode)
     integer(kind=tensor_int),     intent(in)  :: dims(mode)
     integer(kind=tensor_standard_int), intent(in)  :: tdim(mode)
     integer(kind=tensor_long_int), optional, intent(out) :: ntiles
     include "get_ntpm.inc"
  end subroutine tensor_get_ntpm8488
  subroutine tensor_get_ntpm4444(dims,tdim,mode,ntpm,ntiles)
     implicit none
     integer(kind=tensor_standard_int), intent(in)  :: mode
     integer(kind=tensor_standard_int), intent(out) :: ntpm(mode)
     integer(kind=tensor_int),     intent(in)  :: dims(mode)
     integer(kind=tensor_standard_int), intent(in)  :: tdim(mode)
     integer(kind=tensor_standard_int), optional, intent(out) :: ntiles
     include "get_ntpm.inc"
  end subroutine tensor_get_ntpm4444

  subroutine memory_deallocate_window(arr)
     implicit none
     type(tensor),intent(inout) :: arr
     integer(kind=tensor_long_int) :: vector_size
     real(tensor_dp) :: tcpu1,twall1,tcpu2,twall2
     integer :: i,ierr

     call LSTIMER('START',tcpu1,twall1,lspdm_stdout)

     !call memory_deallocate_array(arr)
     if(associated(arr%wi)) then

#ifdef VAR_MPI
        do i=1,arr%nwins
           call lsmpi_win_free(arr%wi(i))
        enddo
#else
        print *,"THIS MESSAGE SHOULD NEVER APPEAR,WHY ARE MPI_WINs ALLOCD?"
#endif

        vector_size = int(arr%nwins*tensor_mpi_kind,kind=tensor_long_int)
        call tensor_free_mem(arr%wi)
        !$OMP CRITICAL
        tensor_counter_aux_f_mem     = tensor_counter_aux_f_mem     + vector_size
        tensor_counter_memory_in_use = tensor_counter_memory_in_use - vector_size
        !$OMP END CRITICAL

        arr%nwins = 0

     endif

     call LSTIMER('START',tcpu2,twall2,lspdm_stdout)

  end subroutine memory_deallocate_window

  subroutine memory_allocate_window(arr,nwins)
     implicit none
     type(tensor),intent(inout) :: arr
     integer, intent(in), optional :: nwins
     integer(kind=tensor_long_int) :: vector_size
     real(tensor_dp) :: tcpu1,twall1,tcpu2,twall2
     integer :: n

     call LSTIMER('START',tcpu1,twall1,lspdm_stdout)

     !call memory_deallocate_array(arr)
     if(associated(arr%wi)) then
        call lsquit("ERROR(memory_allocate_window):array already initialized, please free first",-1)
     endif

     if(present(nwins))then
        arr%nwins = nwins
     else
        arr%nwins = arr%ntiles
     endif

     vector_size = int( arr%nwins * tensor_mpi_kind, kind=tensor_long_int ) 

     call tensor_alloc_mem( arr%wi, arr%nwins )

     !$OMP CRITICAL
     tensor_counter_aux_a_mem     = tensor_counter_aux_a_mem     + vector_size
     tensor_counter_memory_in_use = tensor_counter_memory_in_use + vector_size
     !$OMP END CRITICAL

     call LSTIMER('START',tcpu2,twall2,lspdm_stdout)

  end subroutine memory_allocate_window

  !> \brief Allocate memory for general arrays with memory statistics and tiled
  !structure of the data
  !> \author Patrick Ettenhuber 
  subroutine memory_allocate_tiles(arr,bg)
     implicit none
     type(tensor) :: arr
     logical, intent(in) :: bg
     integer(kind=tensor_long_int) :: vector_size
     real(tensor_dp) :: tcpu1,twall1,tcpu2,twall2
     integer(kind=long) :: i,counter
     integer :: j,loc_idx
     integer, pointer :: idx(:)
     integer(kind=tensor_mpi_kind) :: ibuf(2)
     logical :: doit=.true.,lg_master
     integer(kind=tensor_mpi_kind) :: lg_me,lg_nnod,pc_me,pc_nnod
     integer(kind=tensor_long_int)       :: ne,tooo
     integer               :: res,from
     integer(kind=tensor_standard_int) :: tidx
     lg_master = .true.
     lg_nnod   = 1
     lg_me     = 0


#ifdef VAR_MPI
     lg_nnod   = infpar%lg_nodtot
     lg_me     = infpar%lg_mynum
     lg_master = (infpar%master==lg_me)
#else
     call lsquit("Do not use MPI-WINDOWS without MPI",lspdm_errout)
#endif

     !get zero dummy matrix in size of largest tile --> size(dummy)=tsize
     if( alloc_in_dummy )then
        call memory_allocate_dummy(arr,bg, nel = int(arr%tsize*arr%nlti,kind=tensor_long_int))
#ifdef VAR_MPI
        call memory_allocate_window(arr,nwins = 1)
        call lsmpi_win_create(arr%dummy,arr%wi(1),arr%tsize*arr%nlti,infpar%lg_comm) 
#endif
     else
        call memory_allocate_dummy( arr,bg, nel = 1_tensor_long_int)
        !prepare the integer window in the array --> ntiles windows should be created
        call memory_allocate_window(arr)
     endif

     call LSTIMER('START',tcpu1,twall1,lspdm_stdout)
     call tensor_alloc_mem(idx,arr%mode)

     !write(*,'(I2," in here and nlti ",I5)'),infpar%lg_mynum,arr%nlti
     call tensor_alloc_mem(arr%ti,arr%nlti)

     vector_size = int(arr%nlti*tensor_bytes_per_tile,kind=tensor_long_int)
     !$OMP CRITICAL
     tensor_counter_aux_a_mem     = tensor_counter_aux_a_mem     + vector_size
     tensor_counter_memory_in_use = tensor_counter_memory_in_use + vector_size
     !$OMP END CRITICAL

     counter = 1 
     !allocate tiles with zeros wherever there is mod --> this is experimental
     !and not recommended
     if(arr%zeros)then
        do i=1,arr%ntiles
           doit=.true.
#ifdef VAR_MPI
           if(.not.mod(i+arr%offset,infpar%lg_nodtot)==infpar%lg_mynum)doit=.false.
#endif
           if(doit)then
              call tensor_alloc_mem(arr%ti(counter)%t,arr%tsize)
              call tensor_alloc_mem(arr%ti(counter)%d,arr%mode)
              if( tensor_debug_mode )then
                 arr%ti(counter)%t=0.0E0_tensor_dp
              endif
              arr%ti(counter)%e=arr%tsize
              vector_size = int((arr%ti(counter)%e)*tensor_dp,kind=tensor_long_int)
              !$OMP CRITICAL
              tensor_counter_tiled_a_mem   = tensor_counter_tiled_a_mem + vector_size
              tensor_counter_memory_in_use  = tensor_counter_memory_in_use    + vector_size
              vector_size                  = int(size(arr%ti(counter)%d)*tensor_standard_int,kind=tensor_long_int)
              tensor_counter_aux_a_mem     = tensor_counter_aux_a_mem    + vector_size
              tensor_counter_memory_in_use  = tensor_counter_memory_in_use + vector_size
              !tensor_counter_max_memory    = max(tensor_counter_max_memory,tensor_counter_memory_in_use)
              !$OMP END CRITICAL
              do j=1,arr%mode
                 arr%ti(counter)%d(j)=arr%tdim(j)
              enddo
              counter = counter +1
           endif
        enddo
        !allocate tiles with the rim-tiles of mod dimensions
     else

        do i=1,arr%ntiles

           tidx = i

           !Check if the current tile resides on the current node
           doit=.true.
           call get_residence_of_tile(res,tidx,arr,idx_on_node = from)
#ifdef VAR_MPI
           if(arr%itype==TT_TILED_DIST.and.res/=lg_me) doit = .false.
#endif
           !convert global tile index i to local tile index loc_idx
           loc_idx=((i-1)/lg_nnod) + 1

           if (arr%itype==TT_TILED_REPL.or.arr%itype==TT_TILED)then
              loc_idx = i 
              from    = 1 + ( loc_idx - 1 ) * arr%tsize
           endif

           if(doit)then
              !only do these things if tile belongs to node
              !save global tile number
              arr%ti(loc_idx)%gt=i

              call tensor_alloc_mem(arr%ti(loc_idx)%d,arr%mode)
              vector_size = int(size(arr%ti(loc_idx)%d)*tensor_standard_int,kind=tensor_long_int)
              !$OMP CRITICAL
              tensor_counter_aux_a_mem = tensor_counter_aux_a_mem + vector_size
              tensor_counter_memory_in_use  = tensor_counter_memory_in_use  + vector_size
              !$OMP END CRITICAL

              !get the actual tile size of current tile, here the index is per mode
              call get_tile_dim(arr%ti(loc_idx)%d,arr,i)
              !calculate the number of elements from the tile dimensions
              arr%ti(loc_idx)%e=1
              do j=1,arr%mode
                 arr%ti(loc_idx)%e=arr%ti(loc_idx)%e*arr%ti(loc_idx)%d(j)
              enddo

#ifdef VAR_MPI
              if( alloc_in_dummy )then
                 tooo = from + arr%ti(loc_idx)%e - 1
                 arr%ti(loc_idx)%t => arr%dummy(from:tooo)
              else
                 if(bg)then
                    call mem_pseudo_alloc(arr%ti(loc_idx)%t,arr%ti(loc_idx)%e)
                 else
                    call tensor_alloc_mem(arr%ti(loc_idx)%t,arr%ti(loc_idx)%c,arr%ti(loc_idx)%e)
                 endif
                 vector_size = int(arr%ti(loc_idx)%e*tensor_dp,kind=tensor_long_int)
                 !$OMP CRITICAL
                 tensor_counter_tiled_a_mem = tensor_counter_tiled_a_mem + vector_size
                 tensor_counter_memory_in_use    = tensor_counter_memory_in_use + vector_size
                 !tensor_counter_max_memory       = max(tensor_counter_max_memory,tensor_counter_memory_in_use)
                 !$OMP END CRITICAL
                 call lsmpi_win_create(arr%ti(loc_idx)%t,arr%wi(i),arr%ti(loc_idx)%e,infpar%lg_comm) 
              endif


              if( tensor_debug_mode )then
                 arr%ti(loc_idx)%t=0.0E0_tensor_dp
              endif
#endif

              counter = counter +1
           else

              !open a window of size zero on the nodes where the tile does not
              !reside
#ifdef VAR_MPI
              if( .not. alloc_in_dummy  )call lsmpi_win_create(arr%dummy,arr%wi(i),0,infpar%lg_comm)
#endif
           endif
#ifdef VAR_MPI
           ! fence the window in preparation for comm
           if( .not. alloc_in_dummy )call lsmpi_win_fence(arr%wi(i),.true.)
#endif
        enddo
     endif
     call tensor_free_mem(idx)

     if(counter-1/=arr%nlti)then
        print*," counted wrong of node",lg_me,lg_nnod,counter-1,arr%nlti,arr%ntiles,arr%offset
        call lsquit("something went wrong with the numbering of tiles on the nodes",lspdm_errout)
     endif
     call LSTIMER('START',tcpu2,twall2,lspdm_stdout)

  end subroutine memory_allocate_tiles

  subroutine memory_allocate_dummy(arr,bg,nel)
     implicit none
     type(tensor),intent(inout) :: arr
     logical, intent(in) :: bg
     integer(kind=tensor_long_int), intent(in),optional :: nel
     integer(kind=tensor_long_int) :: nelms
     integer(kind=tensor_long_int) :: vector_size
     real(tensor_dp)     :: tcpu1,twall1,tcpu2,twall2
     integer(kind=tensor_long_int) :: ne

     call LSTIMER('START',tcpu1,twall1,lspdm_stdout)

     if(bg.and..not.mem_is_background_buf_init())then
        call lsquit("ERROR(memory_allocate_dummy): allocation in bg buffer requested, but not bg buffer is allocated",-1)
     endif

     if(associated(arr%dummy)) then
        call lsquit("ERROR(memory_allocate_dummy):array already initialized, please free first",lspdm_errout)
     endif
     if(present(nel))then
        nelms=nel
     else
        nelms=arr%tsize
     endif

     vector_size = int(nelms*tensor_dp,kind=tensor_long_int)

     call tensor_alloc_mem(arr%dummy,arr%dummyc,nelms,bg=bg)

     if( tensor_debug_mode )then
        arr%dummy = 0.0E0_tensor_dp
     endif

     !$OMP CRITICAL
     tensor_counter_aux_a_mem = tensor_counter_aux_a_mem + vector_size
     tensor_counter_memory_in_use = tensor_counter_memory_in_use + vector_size
     !$OMP END CRITICAL

     call LSTIMER('START',tcpu2,twall2,lspdm_stdout)

  end subroutine memory_allocate_dummy

  subroutine memory_deallocate_dummy(arr)
     implicit none
     type(tensor) :: arr
     integer(kind=tensor_long_int) :: vector_size,dim1
     real(tensor_dp) :: tcpu1,twall1,tcpu2,twall2
     logical :: bg

     call LSTIMER('START',tcpu1,twall1,lspdm_stdout)

     if(associated(arr%dummy)) then

        vector_size = int(size(arr%dummy)*tensor_dp,kind=tensor_long_int)

        bg = arr%bg_alloc
        call tensor_free_mem(arr%dummy,arr%dummyc,bg=bg)

        !$OMP CRITICAL
        tensor_counter_aux_f_mem     = tensor_counter_aux_f_mem     + vector_size
        tensor_counter_memory_in_use = tensor_counter_memory_in_use - vector_size
        !$OMP END CRITICAL
     end if

     call LSTIMER('START',tcpu2,twall2,lspdm_stdout)


  end subroutine memory_deallocate_dummy


  subroutine tensor_pdm_free_special_aux(arr)
     implicit none
     type(tensor), intent(inout) :: arr
     call memory_deallocate_window(arr)
     call memory_deallocate_dummy(arr)
  end subroutine tensor_pdm_free_special_aux


  subroutine get_residence_of_tile44(rank_of_node,globaltilenumber,arr,pos_on_node,idx_on_node,window_index)
     implicit none
     integer(kind=tensor_standard_int), intent(out) :: rank_of_node
     type(tensor), intent(in) :: arr
     integer(kind=tensor_standard_int),intent(in) :: globaltilenumber
     integer, intent(out), optional :: pos_on_node, idx_on_node, window_index
     integer :: nnod, pos, idx,widx
     integer(kind=tensor_long_int) :: rank8, gt

     gt = globaltilenumber

     call get_residence_of_tile88(rank8,gt,arr,&
        &pos_on_node=pos_on_node,idx_on_node=idx_on_node,window_index=window_index)

     rank_of_node = int(rank8,kind=4)

  end subroutine get_residence_of_tile44
  subroutine get_residence_of_tile48(rank_of_node,globaltilenumber,arr,pos_on_node,idx_on_node,window_index)
     implicit none
     integer(kind=tensor_standard_int), intent(out) :: rank_of_node
     type(tensor), intent(in) :: arr
     integer(kind=tensor_long_int),intent(in) :: globaltilenumber
     integer, intent(out), optional :: pos_on_node, idx_on_node, window_index
     integer :: nnod, pos, idx,widx
     integer(kind=tensor_long_int) :: rank8


     call get_residence_of_tile88(rank8,globaltilenumber,arr,&
        &pos_on_node=pos_on_node,idx_on_node=idx_on_node,window_index=window_index)

     rank_of_node = int(rank8,kind=4)

  end subroutine get_residence_of_tile48

  subroutine get_residence_of_tile84(rank_of_node,globaltilenumber,arr,pos_on_node,idx_on_node,window_index)
     implicit none
     integer(kind=tensor_long_int), intent(out) :: rank_of_node
     type(tensor), intent(in) :: arr
     integer(kind=tensor_standard_int),intent(in) :: globaltilenumber
     integer, intent(out), optional :: pos_on_node, idx_on_node, window_index
     integer :: nnod, pos, idx,widx
     integer(kind=tensor_long_int) :: gt

     gt = globaltilenumber

     call get_residence_of_tile88(rank_of_node,gt,arr,&
        &pos_on_node=pos_on_node,idx_on_node=idx_on_node,window_index=window_index)

  end subroutine get_residence_of_tile84

  subroutine get_residence_of_tile88(rank_of_node,globaltilenumber,arr,pos_on_node,idx_on_node,window_index)
     implicit none
     integer(kind=tensor_long_int), intent(out) :: rank_of_node
     type(tensor), intent(in) :: arr
     integer(kind=tensor_long_int),intent(in) :: globaltilenumber
     integer, intent(out), optional :: pos_on_node, idx_on_node,window_index
     integer :: nnod, pos, idx,widx

     nnod=arr%nnod
     rank_of_node = mod(globaltilenumber-1+arr%offset,nnod)

     if( alloc_in_dummy ) then
        widx = 1
        pos  = (globaltilenumber-1)/nnod + 1
        idx  = 1 + ( pos - 1 ) * arr%tsize
     else
        widx = globaltilenumber
        pos  = 1
        idx  = 1 
     endif

     !Return the node local index of the tile
     if(present(pos_on_node)) pos_on_node = pos
     if(present(idx_on_node)) idx_on_node = idx
     if(present(window_index)) window_index = widx
  end subroutine get_residence_of_tile88


end module lspdm_basic_module
