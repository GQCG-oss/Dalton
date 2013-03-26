!> @file
!> Memory manager for four dimensional arrays
module array_memory_manager

  use precision
  use ptr_assoc_module!, only: ass_D1to2,ass_D1to3, &
!         &ass_D1to4, ass_D1to5, ass_D1to6, ass_D1to7
  use LSTIMING!,only:lstimer
  use memory_handling!, only: mem_alloc,mem_dealloc
#ifdef VAR_LSMPI
  use infpar_module
  use lsmpi_type!, only: lsmpi_localwin_create_realk,&
!       & lsmpi_fence, lsmpi_win_free, lsmpi_barrier,&
!       & lsmpi_first_fence, lsmpi_last_fence
#endif
  use dec_typedef_module

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
  integer, parameter :: SCALAPACK=5

  !parameters for PDMTYPE:
  integer,parameter :: NO_PDM=0
  integer,parameter :: MASTER_INIT=1
  integer,parameter :: ALL_INIT=2

  !other parameters
  integer,parameter :: ARR_MSG_LEN=30
  integer,parameter :: DEFAULT_TDIM=10

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

    !> \brief to accurately account for mem on each node
    ! count mem allocated in each of the arrays many pointers
    !> \author Patrick Ettenhuber
    subroutine arr_set_dims(arr,dims,mode)
      implicit none
      type(array),intent(inout) :: arr
      integer,intent(in) :: dims(*)
      integer, optional :: mode
      real(realk) :: vector_size
      integer :: i
      if (present(mode))then
        if((arr%mode/=0 .and. arr%mode/=mode).or.arr%mode==0)then
          print *,"mode",mode,"arr%mode",arr%mode      
          call lsquit("wrong use of arr_set_dims",DECinfo%output)
        else
          arr%mode=mode
        endif
      endif
      if (associated(arr%dims))then
!$OMP CRITICAL
        vector_size = dble(size(arr%dims))
        array_aux_deallocd_mem = array_aux_deallocd_mem + vector_size
        array_memory_in_use    = array_memory_in_use  - vector_size
!$OMP END CRITICAL
        call mem_dealloc(arr%dims)
      endif
      if(.not.associated(arr%dims))then
        call mem_alloc(arr%dims,arr%mode)
!$OMP CRITICAL
        vector_size = dble(size(arr%dims))
        array_aux_allocd_mem = array_aux_allocd_mem + vector_size
        array_memory_in_use  = array_memory_in_use  + vector_size
!$OMP END CRITICAL
      endif
      do i=1,arr%mode
        arr%dims(i)=dims(i)
      enddo
    end subroutine arr_set_dims

    !> \brief to accurately account for mem on each node
    ! count mem allocated in each of the arrays many pointers
    !> \author Patrick Ettenhuber
    subroutine arr_set_tdims(arr,tdims,mode)
      implicit none
      type(array),intent(inout) :: arr
      integer,intent(in) :: tdims(*)
      integer, optional :: mode
      real(realk) :: vector_size
      integer :: i
      if (present(mode))then
        if((arr%mode/=0 .and. arr%mode/=mode).or.arr%mode==0)then
          print *,"mode",mode,"arr%mode",arr%mode      
          call lsquit("wrong use of arr_set_tdim",DECinfo%output)
        else
          arr%mode=mode
        endif
      endif
      if (associated(arr%tdim))then
!$OMP CRITICAL
        vector_size = dble(size(arr%tdim))
        array_aux_deallocd_mem = array_aux_deallocd_mem + vector_size
        array_memory_in_use    = array_memory_in_use  - vector_size
!$OMP END CRITICAL
        call mem_dealloc(arr%tdim)
      endif
      if(.not.associated(arr%tdim))then
        call mem_alloc(arr%tdim,arr%mode)
!$OMP CRITICAL
        vector_size = dble(size(arr%tdim))
        array_aux_allocd_mem = array_aux_allocd_mem + vector_size
        array_memory_in_use  = array_memory_in_use  + vector_size
!$OMP END CRITICAL
      endif
      do i=1,arr%mode
        arr%tdim(i)=tdims(i)
      enddo
    end subroutine arr_set_tdims

    !> \brief to accurately account for mem on each node
    ! count mem allocated in each of the arrays many pointers
    !> \author Patrick Ettenhuber
    subroutine arr_set_ntpm(arr,ntpm,mode)
      implicit none
      type(array),intent(inout) :: arr
      integer,intent(in) :: ntpm(*)
      integer, optional :: mode
      real(realk) :: vector_size
      integer :: i
      if (present(mode))then
        if((arr%mode/=0 .and. arr%mode/=mode).or.arr%mode==0)then
          print *,"mode",mode,"arr%mode",arr%mode      
          call lsquit("wrong use of arr_set_ntpm",DECinfo%output)
        else
          arr%mode=mode
        endif
      endif
      if (associated(arr%ntpm))then
!$OMP CRITICAL
        vector_size = dble(size(arr%ntpm))
        array_aux_deallocd_mem = array_aux_deallocd_mem + vector_size
        array_memory_in_use    = array_memory_in_use  - vector_size
!$OMP END CRITICAL
        call mem_dealloc(arr%ntpm)
      endif
      if(.not.associated(arr%ntpm))then
        call mem_alloc(arr%ntpm,arr%mode)
!$OMP CRITICAL
        vector_size = dble(size(arr%ntpm))
        array_aux_allocd_mem = array_aux_allocd_mem + vector_size
        array_memory_in_use  = array_memory_in_use  + vector_size
!$OMP END CRITICAL
      endif
      do i=1,arr%mode
        arr%ntpm(i)=ntpm(i)
      enddo
    end subroutine arr_set_ntpm

    subroutine arr_set_addr(arr,addr,nnodes)
      implicit none
      type(array),intent(inout) :: arr
      integer,intent(in) :: addr(*)
      integer(kind=ls_mpik) :: nnodes
      real(realk) :: vector_size
      integer :: i
      if (associated(arr%addr_p_arr))then
!$OMP CRITICAL
        vector_size = dble(size(arr%addr_p_arr))
        array_aux_deallocd_mem = array_aux_deallocd_mem + vector_size
        array_memory_in_use    = array_memory_in_use  - vector_size
!$OMP END CRITICAL
        call mem_dealloc(arr%addr_p_arr)
      endif
      if(.not.associated(arr%addr_p_arr))then
        call mem_alloc(arr%addr_p_arr,nnodes)
!$OMP CRITICAL
        vector_size = dble(size(arr%addr_p_arr))
        array_aux_allocd_mem = array_aux_allocd_mem + vector_size
        array_memory_in_use  = array_memory_in_use  + vector_size
!$OMP END CRITICAL
      endif
      do i=1,nnodes
        arr%addr_p_arr(i)=addr(i)
      enddo
    end subroutine arr_set_addr


    subroutine arr_free_aux(arr)
      implicit none
      type(array),intent(inout) :: arr
      real(realk) :: vector_size
      if (associated(arr%dims))then
!$OMP CRITICAL
        vector_size = dble(size(arr%dims))
        array_aux_deallocd_mem = array_aux_deallocd_mem + vector_size
        array_memory_in_use    = array_memory_in_use  - vector_size
!$OMP END CRITICAL
        call mem_dealloc(arr%dims)
      endif
      if (associated(arr%tdim))then
!$OMP CRITICAL
        vector_size = dble(size(arr%tdim))
        array_aux_deallocd_mem = array_aux_deallocd_mem + vector_size
        array_memory_in_use    = array_memory_in_use  - vector_size
!$OMP END CRITICAL
        call mem_dealloc(arr%tdim)
      endif
      if (associated(arr%ntpm))then
!$OMP CRITICAL
        vector_size = dble(size(arr%ntpm))
        array_aux_deallocd_mem = array_aux_deallocd_mem + vector_size
        array_memory_in_use    = array_memory_in_use  - vector_size
!$OMP END CRITICAL
        call mem_dealloc(arr%ntpm)
      endif
      if (associated(arr%addr_p_arr))then
!$OMP CRITICAL
        vector_size = dble(size(arr%addr_p_arr))
        array_aux_deallocd_mem = array_aux_deallocd_mem + vector_size
        array_memory_in_use    = array_memory_in_use  - vector_size
!$OMP END CRITICAL
        call mem_dealloc(arr%addr_p_arr)
      endif
      if(array_aux_deallocd_mem > array_aux_allocd_mem) &
      &print *,"WARNING(arr_free_aux)more memory deallocated than allocated"
      call memory_deallocate_window(arr)
      call memory_deallocate_dummy(arr)
    end subroutine arr_free_aux






  !> \brief Allocate memory for general arrays with memory statistics
  !> \author Patrick Ettenhuber adapted from Marcin Ziolkowski
    subroutine memory_allocate_array_dense(arr)
      implicit none
      type(array),intent(inout) :: arr
      real(realk) :: vector_size
      real(realk) :: tcpu1,twall1,tcpu2,twall2

      call LSTIMER('START',tcpu1,twall1,DECinfo%output)

      !call memory_deallocate_array(arr)
      if(associated(arr%elm1)) then
        call lsquit("ERROR(memory_allocate_array):array already initialized, please free first",DECinfo%output)
      endif
      vector_size = dble(arr%nelms)*realk

      call mem_alloc(arr%elm1,arr%nelms)

!$OMP CRITICAL
      array_dense_allocd_mem = array_dense_allocd_mem + vector_size
      array_memory_in_use = array_memory_in_use + vector_size
      array_max_memory = max(array_max_memory,array_memory_in_use)
!$OMP END CRITICAL
      
      call assoc_ptr_arr(arr)

      call LSTIMER('START',tcpu2,twall2,DECinfo%output)
      DECinfo%memallo_time_cpu = DECinfo%memallo_time_cpu + (tcpu2-tcpu1)
      DECinfo%memallo_time_wall = DECinfo%memallo_time_wall + (twall2-twall1)

    end subroutine memory_allocate_array_dense

  !> \brief Deallocate memory for 4d arrays with memory statistics
  !> \author Marcin Ziolkowski
  !> \param mat Four dimensional pointer to be deallocated
    subroutine memory_deallocate_array_dense(arr)
      implicit none
      type(array) :: arr
      real(realk) :: vector_size
      real(realk) :: dim1
      real(realk) :: tcpu1,twall1,tcpu2,twall2

      call LSTIMER('START',tcpu1,twall1,DECinfo%output)

      if(associated(arr%elm1)) then
         dim1 = dble(size(arr%elm1(:)))
         vector_size = dim1*realk
         call deassoc_ptr_arr(arr)
         call mem_dealloc(arr%elm1)
         arr%elm1 => null()
!$OMP CRITICAL
         array_dense_deallocd_mem = array_dense_deallocd_mem + vector_size
         array_memory_in_use = array_memory_in_use - vector_size
!$OMP END CRITICAL
      end if

      call LSTIMER('START',tcpu2,twall2,DECinfo%output)

      DECinfo%memallo_time_cpu = DECinfo%memallo_time_cpu + (tcpu2-tcpu1)
      DECinfo%memallo_time_wall = DECinfo%memallo_time_wall + (twall2-twall1)

    end subroutine memory_deallocate_array_dense


    function get_residence_of_tile(globaltilenumber,arr) result(rankofnode)
      implicit none
      type(array), intent(in) :: arr
      integer,intent(in) :: globaltilenumber
      integer :: rankofnode,nnod
      nnod=1
#ifdef VAR_LSMPI
      nnod=infpar%lg_nodtot
#endif      
      rankofnode=mod(globaltilenumber-1+arr%offset,nnod)
    end function get_residence_of_tile


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
#ifdef VAR_LSMPI
      !print *,infpar%lg_mynum,"in routine"
      nnod=infpar%lg_nodtot
      me=infpar%lg_mynum
      if(infpar%master/=me)master=.false.
      nnod=infpar%lg_nodtot
#else
      call lsquit("Do not use MPI-WINDOWS wihout MPI",DECinfo%output)
#endif
      !get zero dummy matrix in size of largest tile --> size(dummy)=tsize
      call memory_allocate_dummy(arr)
      !prepare the integer window in the array --> ntiles windows should be created
      call memory_allocate_window(arr)

      call LSTIMER('START',tcpu1,twall1,DECinfo%output)
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
#ifdef VAR_LSMPI
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
#ifdef VAR_LSMPI
          if(arr%atype==TILED_DIST.and.&
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

#ifdef VAR_LSMPI
            !if(act_ts/=arr%ti(loc_idx)%e)print*,"something wrong"
            call mem_alloc(arr%ti(loc_idx)%t,arr%ti(loc_idx)%c,arr%ti(loc_idx)%e)
            call lsmpi_win_create_realk(arr%ti(loc_idx)%t,arr%wi(i),arr%ti(loc_idx)%e,infpar%lg_comm)
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
#ifdef VAR_LSMPI
            call lsmpi_win_create_realk(arr%dummy,arr%wi(i),0,infpar%lg_comm)
#endif
          endif
#ifdef VAR_LSMPI
          call lsmpi_win_fence(arr%wi(i),.true.)
#endif
        enddo
      endif
      call mem_dealloc(idx)

      if(counter-1/=arr%nlti)then
        print*,counter-1,arr%nlti,arr%ntiles
        call lsquit("something went wrong with the numbering of tiles on the nodes",DECinfo%output)
      endif
      call LSTIMER('START',tcpu2,twall2,DECinfo%output)
      DECinfo%memallo_time_cpu = DECinfo%memallo_time_cpu + (tcpu2-tcpu1)
      DECinfo%memallo_time_wall = DECinfo%memallo_time_wall + (twall2-twall1)

    end subroutine memory_allocate_tiles

  !> \brief deallocate tiles and keep track of memory
  !> \author Patrick Ettenhuber
  !> \date September 2012
    subroutine memory_deallocate_tile(arr)
      implicit none
      type(array),intent(inout) :: arr
      real(realk) :: vector_size
      real(realk) :: dim1,dim2,dim3,dim4
      real(realk) :: tcpu1,twall1,tcpu2,twall2
      integer :: i

      call LSTIMER('START',tcpu1,twall1,DECinfo%output)

      do i=1,arr%nlti
        if(associated(arr%ti(i)%t)) then
 
           dim1 = dble(size(arr%ti(i)%t))*realk
           if(.false.)then
             call mem_dealloc(arr%ti(i)%t)
#ifdef VAR_LSMPI
           else
             call mem_dealloc(arr%ti(i)%t,arr%ti(i)%c)
#endif
           endif

           dim2 = dble(size(arr%ti(i)%d(:)))
           nullify(arr%ti(i)%t)
           call mem_dealloc(arr%ti(i)%d)
!$OMP CRITICAL
           array_tiled_deallocd_mem = array_tiled_deallocd_mem + dim1 
           array_memory_in_use      = array_memory_in_use      - dim1
           array_aux_deallocd_mem   = array_aux_deallocd_mem   + dim2 
           array_memory_in_use      = array_memory_in_use      - dim2
!$OMP END CRITICAL
        end if

      enddo
      vector_size = dble(size(arr%ti))
!$OMP CRITICAL
      array_aux_deallocd_mem = array_aux_deallocd_mem + vector_size
      array_memory_in_use    = array_memory_in_use    - vector_size
!$OMP END CRITICAL
      deallocate(arr%ti)
      nullify(arr%ti)
      call LSTIMER('START',tcpu2,twall2,DECinfo%output)

      DECinfo%memallo_time_cpu = DECinfo%memallo_time_cpu + (tcpu2-tcpu1)
      DECinfo%memallo_time_wall = DECinfo%memallo_time_wall + (twall2-twall1)

    end subroutine memory_deallocate_tile

    subroutine memory_allocate_dummy(arr,nel)
      implicit none
      type(array),intent(inout) :: arr
      integer, intent(in),optional :: nel
      integer :: nelms
      real(realk) :: vector_size
      real(realk) :: tcpu1,twall1,tcpu2,twall2

      call LSTIMER('START',tcpu1,twall1,DECinfo%output)

      !call memory_deallocate_array(arr)
      if(associated(arr%dummy)) then
        call lsquit("ERROR(memory_allocate_dummy):array already initialized, please free first",DECinfo%output)
      endif
      if(present(nel))then
        nelms=nel
      else
        nelms=arr%tsize
      endif

      vector_size = dble(nelms)*realk

#ifdef VAR_LSMPI
      call mem_alloc(arr%dummy,arr%dummyc,nelms)
#endif

!$OMP CRITICAL
      array_aux_allocd_mem = array_aux_allocd_mem + vector_size
      array_memory_in_use = array_memory_in_use + vector_size
      array_max_memory = max(array_max_memory,array_memory_in_use)
!$OMP END CRITICAL
      arr%dummy=0.0E0_realk

      call LSTIMER('START',tcpu2,twall2,DECinfo%output)
      DECinfo%memallo_time_cpu = DECinfo%memallo_time_cpu + (tcpu2-tcpu1)
      DECinfo%memallo_time_wall = DECinfo%memallo_time_wall + (twall2-twall1)

    end subroutine memory_allocate_dummy

    subroutine memory_deallocate_dummy(arr)
      implicit none
      type(array) :: arr
      real(realk) :: vector_size
      real(realk) :: dim1
      real(realk) :: tcpu1,twall1,tcpu2,twall2

      call LSTIMER('START',tcpu1,twall1,DECinfo%output)

      if(associated(arr%dummy)) then
         dim1 = dble(size(arr%dummy(:)))
         vector_size = dim1*realk
#ifdef VAR_LSMPI
         call mem_dealloc(arr%dummy,arr%dummyc)
#endif
         nullify(arr%dummy)
!$OMP CRITICAL
         array_aux_deallocd_mem = array_aux_deallocd_mem + vector_size
         array_memory_in_use = array_memory_in_use - vector_size
!$OMP END CRITICAL
      end if

      call LSTIMER('START',tcpu2,twall2,DECinfo%output)

      DECinfo%memallo_time_cpu = DECinfo%memallo_time_cpu + (tcpu2-tcpu1)
      DECinfo%memallo_time_wall = DECinfo%memallo_time_wall + (twall2-twall1)

    end subroutine memory_deallocate_dummy

    subroutine memory_allocate_window(arr)
      implicit none
      type(array),intent(inout) :: arr
      real(realk) :: vector_size
      real(realk) :: tcpu1,twall1,tcpu2,twall2

      call LSTIMER('START',tcpu1,twall1,DECinfo%output)

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

      call LSTIMER('START',tcpu2,twall2,DECinfo%output)
      DECinfo%memallo_time_cpu = DECinfo%memallo_time_cpu + (tcpu2-tcpu1)
      DECinfo%memallo_time_wall = DECinfo%memallo_time_wall + (twall2-twall1)

    end subroutine memory_allocate_window

    subroutine memory_deallocate_window(arr)
      implicit none
      type(array),intent(inout) :: arr
      real(realk) :: vector_size
      real(realk) :: tcpu1,twall1,tcpu2,twall2
      integer :: i,ierr

      call LSTIMER('START',tcpu1,twall1,DECinfo%output)

      !call memory_deallocate_array(arr)
      if(associated(arr%wi)) then
         do i=1,arr%ntiles
#ifdef VAR_LSMPI
           call lsmpi_win_fence(arr%wi(i),.false.)
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

      call LSTIMER('START',tcpu2,twall2,DECinfo%output)
      DECinfo%memallo_time_cpu = DECinfo%memallo_time_cpu + (tcpu2-tcpu1)
      DECinfo%memallo_time_wall = DECinfo%memallo_time_wall + (twall2-twall1)

    end subroutine memory_deallocate_window


    !> \brief the following collection of routines is a workaround the internal
    !fortran routines to be able to associate a one dimensional array with a
    !multidimensional pointer
    !> \author Patrick Ettenhuber
    !> \date September 2012
    subroutine assoc_ptr_arr(arr)
      implicit none
      type(array)::arr
      select case(arr%mode)
        case(2)
                call ass_D1to2(arr%elm1,arr%elm2,arr%dims)
        case(3)
                call ass_D1to3(arr%elm1,arr%elm3,arr%dims)
        case(4)
                call ass_D1to4(arr%elm1,arr%elm4,arr%dims)
        case(5)
                call ass_D1to5(arr%elm1,arr%elm5,arr%dims)
        case(6)
                call ass_D1to6(arr%elm1,arr%elm6,arr%dims)
        case(7)
                call ass_D1to7(arr%elm1,arr%elm7,arr%dims)
        case default
                return
      end select
    end subroutine assoc_ptr_arr
    
    !\brief deassociate all pointers again
    !\author Patrick Ettenhuber
    !\date September 2012
    subroutine deassoc_ptr_arr(arr)
      implicit none
      type(array)::arr
      arr%elm2 => null()
      arr%elm3 => null()
      arr%elm4 => null()
      arr%elm5 => null()
      arr%elm6 => null()
      arr%elm7 => null()
    end subroutine deassoc_ptr_arr



  !> \brief Print statistics of array4 objects
  !> \author Marcin Ziolkowski
  !> \param output File unit for output 
    subroutine print_memory_statistics(output)

      implicit none
      integer, intent(in) :: output

      write(DECinfo%output,'(/,a)')    '  Array memory statistics    '
      write(DECinfo%output,'(a)')      ' =================================================='
      write(DECinfo%output,'(a,f12.2,a)') ' Allocated memory    : ',array_dense_allocd_mem/2**20,' MB'
      write(DECinfo%output,'(a,f12.2,a)') ' Deallocated memory  : ',array_dense_deallocd_mem/2**20,' MB'
      write(DECinfo%output,'(a,f12.2,a)') ' Alloc-dealloc mem   : ', &
           (array_dense_allocd_mem-array_dense_deallocd_mem)/2**20,' MB'
      write(DECinfo%output,'(a,f12.2,a)') ' Memory in use       : ',array_memory_in_use/2**20,' MB'
      write(DECinfo%output,'(a,f12.2,a)') ' Max memory in use   : ',array_max_memory/2**20,' MB'

      write(DECinfo%output,'(a,/,a)') '  Time ', &
           ' ======='
      write(DECinfo%output,'(a,g18.2)') ' CPU Time (s) :', DECinfo%memallo_time_cpu
      write(DECinfo%output,'(a,g18.2)') ' Wall Time (s):', DECinfo%memallo_time_wall

    end subroutine print_memory_statistics
  
  !> \brief Print currenly used memory and max allocated memory so far
  !> \author Marcin Ziolkowski, modified by Patrick Ettenhuber
  !> \param output File unit for output
  subroutine array_print_memory_currents(output,retour)
    implicit none
    !> selects the output
    integer, intent(in) :: output
    !> selects to redcue it on master for checking, or outut instead
    real(realk),intent(inout),optional :: retour(8)
    
    if(.not.present(retour))then
      write(DECinfo%output,'(a,g12.4,a)') ' Allocated memory for dense array   :',&
           & array_dense_allocd_mem/(1.0E9_realk),' GB'
      write(DECinfo%output,'(a,g12.4,a)') ' Deallocated memory for dense array :',&
           & array_dense_deallocd_mem/(1.0E9_realk),' GB'
      write(DECinfo%output,'(a,g12.4,a)') ' Allocated memory for tiled array   :',&
           & array_tiled_allocd_mem/(1.0E9_realk),' GB'
      write(DECinfo%output,'(a,g12.4,a)') ' Deallocated memory for tiled array :',&
           & array_tiled_deallocd_mem/(1.0E9_realk),' GB'
      write(DECinfo%output,'(a,g12.4,a)') ' Allocated aux memory for array     :',&
           & array_aux_allocd_mem/(1.0E9_realk),' GB'
      write(DECinfo%output,'(a,g12.4,a)') ' Dellocated aux memory for array    :',&
           & array_aux_deallocd_mem/(1.0E9_realk),' GB'
      write(DECinfo%output,'(a,g12.4,a)') ' Memory in use for array      :',&
           & array_memory_in_use/(1.0E9_realk),' GB'
      write(DECinfo%output,'(a,g12.4,a)') ' Max memory in use for array  :',&
           & array_max_memory/(1.0E9_realk),' GB'
    else
      retour(1)=array_dense_allocd_mem
      retour(2)=array_dense_deallocd_mem
      retour(3)=array_tiled_allocd_mem
      retour(4)=array_tiled_deallocd_mem
      retour(5)=array_aux_allocd_mem
      retour(6)=array_aux_deallocd_mem
      retour(7)=array_memory_in_use
      retour(8)=array_max_memory
    endif

  end subroutine array_print_memory_currents

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
    !case(4)
    !  a=inds(1)+(inds(2)-1)*dims(1)+(inds(3)-1)*dims(1)*dims(2)+(inds(4)-1)*inds(1)*inds(2)*inds(3)
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


  subroutine copy_array(arr_in,arr_out)
    implicit none
    type(array), intent(in) :: arr_in
    type(array), intent(inout) :: arr_out
    integer :: i
    arr_out%mode = arr_in%mode
    arr_out%nlti = arr_in%nlti
    arr_out%tsize = arr_in%tsize
    arr_out%atype = arr_in%atype
    arr_out%nelms = arr_in%nelms
    arr_out%ntiles = arr_in%ntiles
    arr_out%init_type = arr_in%init_type
    if(associated(arr_in%dims))call arr_set_dims(arr_out,arr_in%dims,arr_out%mode)
    if(associated(arr_in%ntpm))call arr_set_ntpm(arr_out,arr_in%ntpm,arr_out%mode)
    if(associated(arr_in%tdim))call arr_set_tdims(arr_out,arr_in%tdim,arr_out%mode)
    if(associated(arr_in%addr_p_arr))call arr_set_addr(arr_out,arr_in%addr_p_arr,size(arr_in%addr_p_arr,kind=ls_mpik))
    !arr_out%dims = arr_in%dims
    !arr_out%tdim = arr_in%tdim
    !arr_out%ntpm = arr_in%ntpm
    if(associated(arr_in%elm1))then
      call memory_allocate_array_dense(arr_out)
      arr_out%elm1=arr_in%elm1
    endif
    if(associated(arr_in%ti))then
      call memory_allocate_tiles(arr_out)
      do i=1,arr_in%nlti
        arr_out%ti(i)%t=arr_in%ti(i)%t
      enddo
    endif
  end subroutine copy_array
  
  !> \author Patrick Ettenhuber adpted from Marcin Ziolkowski
  !> \date September 2012
  !> \brief free the array structure
  subroutine array_free_basic(arr)
    implicit none
    type(array), intent(inout) :: arr
    call arr_free_aux(arr)
    if(associated(arr%elm1))call memory_deallocate_array_dense(arr)
    if(associated(arr%ti))  call memory_deallocate_tile(arr)
  end subroutine array_free_basic

  !> \author Patrick Ettenhuber
  !> \date September 2012
  !> \brief nullify all pointers in the array
  subroutine array_nullify_pointers(arr)
    implicit none
    type(array), intent(inout) :: arr
    NULLIFY(arr%dims)     
    NULLIFY(arr%elm1)     
    NULLIFY(arr%elm2)     
    NULLIFY(arr%elm3)     
    NULLIFY(arr%elm4)     
    NULLIFY(arr%elm5)     
    NULLIFY(arr%elm6)     
    NULLIFY(arr%elm7)     
    NULLIFY(arr%dummy)    
    NULLIFY(arr%ti)       
    NULLIFY(arr%wi)       
    NULLIFY(arr%ntpm)     
    NULLIFY(arr%tdim)     
    NULLIFY(arr%addr_p_arr)
  end subroutine array_nullify_pointers

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
end module array_memory_manager
