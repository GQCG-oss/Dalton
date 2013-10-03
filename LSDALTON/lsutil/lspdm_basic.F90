!> @file
!> here the one sided wrappers should go, so that interfacing with
!different rma vendors is possible
!> \author Patrick Ettenhuber
!> \date April 2013
module lspdm_basic_module
  use precision
  use LSTIMING!,only:lstimer
  use memory_handling
#ifdef VAR_MPI
  use infpar_module
  use lsmpi_type
#endif
  use reorder_frontend_module
  use tensor_type_def_module

  interface get_tile_dim
    module procedure get_tileinfo_nels_frombas,&
                    &get_tileinfo_nelspmode_frombas,&
                    &get_tileinfo_nels_fromarr8,&
                    &get_tileinfo_nels_fromarr4,&
                    &get_tileinfo_nelspermode_fromarr4,&
                    &get_tileinfo_nelspermode_fromarr8
  end interface get_tile_dim

  !interface get_tile_idx
  !  module procedure get_tile_idx_from_global_idx
  !end interface get_tile_idx

  contains
  !subroutine get_tile_idx_from_global_idx()
  !  implicit none
  !end subroutine get_tile_idx_from_global_idx
  

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

  subroutine memory_deallocate_window(arr)
    implicit none
    type(array),intent(inout) :: arr
    real(realk) :: vector_size
    real(realk) :: tcpu1,twall1,tcpu2,twall2
    integer :: i,ierr

    call LSTIMER('START',tcpu1,twall1,lspdm_stdout)

    !call memory_deallocate_array(arr)
    if(associated(arr%wi)) then
       do i=1,arr%ntiles
#ifdef VAR_MPI
         !call lsmpi_win_fence(arr%wi(i),.false.)
         call lsmpi_win_free(arr%wi(i))
#else
         print *,"THIS MESSAGE SHOULD NEVER APPEAR,WHY ARE MPI_WINs ALLOCD?"
#endif
       enddo
       vector_size = dble(size(arr%wi(:)))
       call mem_dealloc(arr%wi)
!$OMP CRITICAL
       array_aux_deallocd_mem = array_aux_deallocd_mem + vector_size
       array_memory_in_use = array_memory_in_use - vector_size
!$OMP END CRITICAL

    endif

    call LSTIMER('START',tcpu2,twall2,lspdm_stdout)

  end subroutine memory_deallocate_window

  subroutine memory_allocate_window(arr)
    implicit none
    type(array),intent(inout) :: arr
    real(realk) :: vector_size
    real(realk) :: tcpu1,twall1,tcpu2,twall2

    call LSTIMER('START',tcpu1,twall1,lspdm_stdout)

    !call memory_deallocate_array(arr)
    if(associated(arr%wi)) then
      call lsquit("ERROR(memory_allocate_window):array already initialized, please free first",-1)
    endif

    vector_size = dble(arr%ntiles)
   
    call mem_alloc(arr%wi,arr%ntiles)

    !print *,infpar%lg_mynum,"mem update"
!$OMP CRITICAL
    array_aux_allocd_mem = array_aux_allocd_mem + vector_size
    array_memory_in_use = array_memory_in_use + vector_size
    array_max_memory = max(array_max_memory,array_memory_in_use)
!$OMP END CRITICAL

    call LSTIMER('START',tcpu2,twall2,lspdm_stdout)

  end subroutine memory_allocate_window

  !> \brief Allocate memory for general arrays with memory statistics and tiled
  !structure of the data
  !> \author Patrick Ettenhuber 
  subroutine memory_allocate_tiles(arr)
    implicit none
    type(array) :: arr
    real(realk) :: vector_size
    real(realk) :: tcpu1,twall1,tcpu2,twall2
    integer(kind=long) :: i,counter
    integer :: j,loc_idx,nnod,me,mpi_realk,ierr
    integer, pointer :: idx(:)
    logical :: doit=.true.,master
    master = .true.
    nnod=1
    me=0
#ifdef VAR_MPI
    !print *,infpar%lg_mynum,"in routine"
    nnod=infpar%lg_nodtot
    me=infpar%lg_mynum
    if(infpar%master/=me)master=.false.
    nnod=infpar%lg_nodtot
#else
    call lsquit("Do not use MPI-WINDOWS wihout MPI",lspdm_errout)
#endif
    !get zero dummy matrix in size of largest tile --> size(dummy)=tsize
    call memory_allocate_dummy(arr)
    !prepare the integer window in the array --> ntiles windows should be created
    call memory_allocate_window(arr)

    call LSTIMER('START',tcpu1,twall1,lspdm_stdout)
    !call memory_deallocate_array(arr)
    call mem_alloc(idx,arr%mode)

    !write(*,'(I2," in here and nlti ",I5)'),infpar%lg_mynum,arr%nlti
    allocate(arr%ti(arr%nlti))

    vector_size = dble(size(arr%ti))
!$OMP CRITICAL
    array_aux_allocd_mem = array_aux_allocd_mem + vector_size
    array_memory_in_use  = array_memory_in_use  + vector_size
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
          call mem_alloc(arr%ti(counter)%t,arr%tsize)
          call mem_alloc(arr%ti(counter)%d,arr%mode)
          arr%ti(counter)%t=0.0d0
          arr%ti(counter)%e=arr%tsize
          vector_size = dble(arr%ti(counter)%e)*realk
!$OMP CRITICAL
          array_tiled_allocd_mem = array_tiled_allocd_mem + vector_size
          array_memory_in_use    = array_memory_in_use    + vector_size
          vector_size            = dble(size(arr%ti(counter)%d))
          array_aux_allocd_mem   = array_aux_allocd_mem   + vector_size
          array_memory_in_use    = array_memory_in_use    + vector_size
          array_max_memory       = max(array_max_memory,array_memory_in_use)
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

        !Check if the current tile resides on the current node
        doit=.true.
#ifdef VAR_MPI
        if(arr%itype==TILED_DIST.and.&
        &.not.mod(i-1+arr%offset,infpar%lg_nodtot)==infpar%lg_mynum)doit=.false.
#endif
        !convert global tile index i to local tile index loc_idx
        loc_idx=((i-1)/nnod) + 1
         
        if(doit)then
          !only do these things if tile belongs to node
          !save global tile number
          arr%ti(loc_idx)%gt=i
          !if(loc_idx>arr%nlti)then
          !  write( *,'(I2," has wrong index:",I3," of ",I3," in ",I3)'),infpar%lg_mynum,loc_idx,arr%nlti,i
          !endif
          call mem_alloc(arr%ti(loc_idx)%d,arr%mode)
          vector_size = dble(size(arr%ti(loc_idx)%d))
!$OMP CRITICAL
          array_aux_allocd_mem = array_aux_allocd_mem + vector_size
          array_memory_in_use  = array_memory_in_use  + vector_size
!$OMP END CRITICAL
          !get the actual tile size of current tile, here the index is per
          !mode
          call get_tile_dim(arr%ti(loc_idx)%d,arr,i)
          !calculate the number of elements from the tile dimensions
          arr%ti(loc_idx)%e=1
          do j=1,arr%mode
            arr%ti(loc_idx)%e=arr%ti(loc_idx)%e*arr%ti(loc_idx)%d(j)
          enddo

#ifdef VAR_MPI
          !if(act_ts/=arr%ti(loc_idx)%e)print*,"something wrong"
          call mem_alloc(arr%ti(loc_idx)%t,arr%ti(loc_idx)%c,arr%ti(loc_idx)%e)
          call lsmpi_win_create(arr%ti(loc_idx)%t,arr%wi(i),arr%ti(loc_idx)%e,infpar%lg_comm)
#endif

          vector_size = dble(arr%ti(loc_idx)%e)*realk
!$OMP CRITICAL
          array_tiled_allocd_mem = array_tiled_allocd_mem + vector_size
          array_memory_in_use = array_memory_in_use + vector_size
          array_max_memory = max(array_max_memory,array_memory_in_use)
!$OMP END CRITICAL
          !print *,"actual tilesizes intile on:",i,loc_idx,counter,act_ts,infpar%lg_mynum
          counter = counter +1
        else
         !open a window of size zero on the nodes where the tile does not
         !reside
#ifdef VAR_MPI
          call lsmpi_win_create(arr%dummy,arr%wi(i),0,infpar%lg_comm)
#endif
        endif
#ifdef VAR_MPI
        call lsmpi_win_fence(arr%wi(i),.true.)
#endif
      enddo
    endif
    call mem_dealloc(idx)

    if(counter-1/=arr%nlti)then
      print*,counter-1,arr%nlti,arr%ntiles
      call lsquit("something went wrong with the numbering of tiles on the nodes",lspdm_errout)
    endif
    call LSTIMER('START',tcpu2,twall2,lspdm_stdout)

  end subroutine memory_allocate_tiles

    subroutine memory_allocate_dummy(arr,nel)
      implicit none
      type(array),intent(inout) :: arr
      integer, intent(in),optional :: nel
      integer :: nelms
      real(realk) :: vector_size
      real(realk) :: tcpu1,twall1,tcpu2,twall2

      call LSTIMER('START',tcpu1,twall1,lspdm_stdout)

      !call memory_deallocate_array(arr)
      if(associated(arr%dummy)) then
        call lsquit("ERROR(memory_allocate_dummy):array already initialized, please free first",lspdm_errout)
      endif
      if(present(nel))then
        nelms=nel
      else
        nelms=arr%tsize
      endif

      vector_size = dble(nelms)*realk

#ifdef VAR_MPI
      call mem_alloc(arr%dummy,arr%dummyc,nelms)
#endif

!$OMP CRITICAL
      array_aux_allocd_mem = array_aux_allocd_mem + vector_size
      array_memory_in_use = array_memory_in_use + vector_size
      array_max_memory = max(array_max_memory,array_memory_in_use)
!$OMP END CRITICAL
      arr%dummy=0.0E0_realk

      call LSTIMER('START',tcpu2,twall2,lspdm_stdout)

    end subroutine memory_allocate_dummy

    subroutine memory_deallocate_dummy(arr)
      implicit none
      type(array) :: arr
      real(realk) :: vector_size
      real(realk) :: dim1
      real(realk) :: tcpu1,twall1,tcpu2,twall2

      call LSTIMER('START',tcpu1,twall1,lspdm_stdout)

      if(associated(arr%dummy)) then
         dim1 = dble(size(arr%dummy(:)))
         vector_size = dim1*realk
#ifdef VAR_MPI
         call mem_dealloc(arr%dummy,arr%dummyc)
#endif
         nullify(arr%dummy)
!$OMP CRITICAL
         array_aux_deallocd_mem = array_aux_deallocd_mem + vector_size
         array_memory_in_use = array_memory_in_use - vector_size
!$OMP END CRITICAL
      end if

      call LSTIMER('START',tcpu2,twall2,lspdm_stdout)


    end subroutine memory_deallocate_dummy
    subroutine array_pdm_free_special_aux(arr)
      implicit none
      type(array), intent(inout) :: arr
      call memory_deallocate_window(arr)
      call memory_deallocate_dummy(arr)
    end subroutine array_pdm_free_special_aux
    
  
end module lspdm_basic_module
