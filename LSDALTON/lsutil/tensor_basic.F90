!> @file
!> Memory manager for four dimensional arrays
module tensor_basic_module

  use precision
  use ptr_assoc_module!, only: ass_D1to2,ass_D1to3, &
!         &ass_D1to4, ass_D1to5, ass_D1to6, ass_D1to7
  use LSTIMING!,only:lstimer
  use memory_handling!, only: mem_alloc,mem_dealloc
!#ifdef VAR_MPI
!  use infpar_module
!  use lsmpi_type!, only: lsmpi_localwin_create_realk,&
!       & lsmpi_fence, lsmpi_win_free, lsmpi_barrier,&
!       & lsmpi_first_fence, lsmpi_last_fence
!#endif
  use lspdm_basic_module
  use tensor_type_def_module


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
          call lsquit("wrong use of arr_set_dims",lspdm_errout)
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
          call lsquit("wrong use of arr_set_tdim",lspdm_errout)
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
  
    subroutine arr_init_lock_set(arr)
      implicit none
      type(array),intent(inout) :: arr
      real(realk) :: vector_size

      if(.not.associated(arr%lock_set))then
        call mem_alloc(arr%lock_set,arr%ntiles)
!$OMP CRITICAL
        vector_size = dble(size(arr%lock_set))
        array_aux_allocd_mem = array_aux_allocd_mem + vector_size
        array_memory_in_use  = array_memory_in_use  + vector_size
!$OMP END CRITICAL
        arr%lock_set = .false.
      endif
      
    end subroutine arr_init_lock_set

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
          call lsquit("wrong use of arr_set_ntpm",lspdm_errout)
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

    subroutine arr_set_addr(arr,addr,nnodes,loc_addr)
      implicit none
      type(array),intent(inout) :: arr
      integer,intent(in) :: addr(*)
      integer(kind=ls_mpik) :: nnodes
      logical,intent(in),optional :: loc_addr
      real(realk) :: vector_size
      integer :: i
      logical :: loc

      loc=.false.
      if(present(loc_addr)) loc = loc_addr
   
      if(loc)then
        if (associated(arr%addr_loc))then
!$OMP CRITICAL
          vector_size = dble(size(arr%addr_loc))
          array_aux_deallocd_mem = array_aux_deallocd_mem + vector_size
          array_memory_in_use    = array_memory_in_use  - vector_size
!$OMP END CRITICAL
          call mem_dealloc(arr%addr_loc)
        endif
        if(.not.associated(arr%addr_loc))then
          call mem_alloc(arr%addr_loc,nnodes)
!$OMP CRITICAL
          vector_size = dble(size(arr%addr_loc))
          array_aux_allocd_mem = array_aux_allocd_mem + vector_size
          array_memory_in_use  = array_memory_in_use  + vector_size
!$OMP END CRITICAL
        endif
        do i=1,nnodes
          arr%addr_loc(i)=addr(i)
        enddo
      else
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
      endif
    end subroutine arr_set_addr








  !> \brief Allocate memory for general arrays with memory statistics
  !> \author Patrick Ettenhuber adapted from Marcin Ziolkowski
    subroutine memory_allocate_array_dense(arr)
      implicit none
      type(array),intent(inout) :: arr
      real(realk) :: vector_size
      real(realk) :: tcpu1,twall1,tcpu2,twall2
      integer(kind=8) :: ne
      logical :: loc,parent

      call LSTIMER('START',tcpu1,twall1,lspdm_stdout)
      !call memory_deallocate_array(arr)
      if(associated(arr%elm1)) then
        call lsquit("ERROR(memory_allocate_array):array already initialized, please free first",lspdm_errout)
      endif
      vector_size = dble(arr%nelms)*realk

      if( lspdm_use_comm_proc )then
#ifdef VAR_MPI
        ne = 0_long
        if(infpar%pc_mynum==infpar%pc_nodtot-1) ne = arr%nelms
        call mem_alloc(arr%elm1,arr%e1c,ne,arr%w1,infpar%pc_comm,arr%nelms)
#else
        call lsquit("ERROR(memory_allocate_array_dense) not possible without MPI",-1)
#endif
      else
        call mem_alloc(arr%elm1,arr%nelms)
      endif

      if( array_debug_mode )then
         arr%elm1 = 0.0E0_realk
      endif

!$OMP CRITICAL
      array_dense_allocd_mem = array_dense_allocd_mem + vector_size
      array_memory_in_use = array_memory_in_use + vector_size
      array_max_memory = max(array_max_memory,array_memory_in_use)
!$OMP END CRITICAL
      
      call assoc_ptr_arr(arr)

      call LSTIMER('START',tcpu2,twall2,lspdm_stdout)

    end subroutine memory_allocate_array_dense

  !> \brief deallocate dense part
  !> \author Patrick Ettenhuber
  !> \param array for which elm1 should be deallocated
    subroutine memory_deallocate_array_dense(arr)
      implicit none
      type(array) :: arr
      real(realk) :: vector_size
      real(realk) :: dim1
      real(realk) :: tcpu1,twall1,tcpu2,twall2

      call LSTIMER('START',tcpu1,twall1,lspdm_stdout)

      if(associated(arr%elm1)) then
         dim1 = dble(size(arr%elm1(:)))
         vector_size = dim1*realk
         call deassoc_ptr_arr(arr)
         if( lspdm_use_comm_proc )then
           call mem_dealloc(arr%elm1,arr%e1c,arr%w1)
         else
           call mem_dealloc(arr%elm1)
         endif
         arr%elm1 => null()
!$OMP CRITICAL
         array_dense_deallocd_mem = array_dense_deallocd_mem + vector_size
         array_memory_in_use = array_memory_in_use - vector_size
!$OMP END CRITICAL
      end if

      call LSTIMER('START',tcpu2,twall2,lspdm_stdout)


    end subroutine memory_deallocate_array_dense





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

      call LSTIMER('START',tcpu1,twall1,lspdm_stdout)

      do i=1,arr%nlti
        if(associated(arr%ti(i)%t)) then
 
           dim1 = dble(size(arr%ti(i)%t))*realk
           if( lspdm_use_comm_proc )then
             call mem_dealloc(arr%ti(i)%t,arr%ti(i)%c,arr%ti(i)%wi)
#ifdef VAR_MPI
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
      call LSTIMER('START',tcpu2,twall2,lspdm_stdout)


    end subroutine memory_deallocate_tile





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

      write(lspdm_stdout,'(/,a)')    '  Array memory statistics    '
      write(lspdm_stdout,'(a)')      ' =================================================='
      write(lspdm_stdout,'(a,f12.2,a)') ' Allocated memory    : ',array_dense_allocd_mem/2**20,' MB'
      write(lspdm_stdout,'(a,f12.2,a)') ' Deallocated memory  : ',array_dense_deallocd_mem/2**20,' MB'
      write(lspdm_stdout,'(a,f12.2,a)') ' Alloc-dealloc mem   : ', &
           (array_dense_allocd_mem-array_dense_deallocd_mem)/2**20,' MB'
      write(lspdm_stdout,'(a,f12.2,a)') ' Memory in use       : ',array_memory_in_use/2**20,' MB'
      write(lspdm_stdout,'(a,f12.2,a)') ' Max memory in use   : ',array_max_memory/2**20,' MB'

      write(lspdm_stdout,'(a,/,a)') '  Time ', &
           ' ======='

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
      write(*,'(a,g12.4,a)') ' Allocated memory for dense array   :',&
           & array_dense_allocd_mem/(1.0E9_realk),' GB'
      write(*,'(a,g12.4,a)') ' Deallocated memory for dense array :',&
           & array_dense_deallocd_mem/(1.0E9_realk),' GB'
      write(*,'(a,g12.4,a)') ' Allocated memory for tiled array   :',&
           & array_tiled_allocd_mem/(1.0E9_realk),' GB'
      write(*,'(a,g12.4,a)') ' Deallocated memory for tiled array :',&
           & array_tiled_deallocd_mem/(1.0E9_realk),' GB'
      write(*,'(a,g12.4,a)') ' Allocated aux memory for array     :',&
           & array_aux_allocd_mem/(1.0E9_realk),' GB'
      write(*,'(a,g12.4,a)') ' Dellocated aux memory for array    :',&
           & array_aux_deallocd_mem/(1.0E9_realk),' GB'
      write(*,'(a,g12.4,a)') ' Memory in use for array      :',&
           & array_memory_in_use/(1.0E9_realk),' GB'
      write(*,'(a,g12.4,a)') ' Max memory in use for array  :',&
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
      if (associated(arr%addr_loc))then
!$OMP CRITICAL
        vector_size = dble(size(arr%addr_loc))
        array_aux_deallocd_mem = array_aux_deallocd_mem + vector_size
        array_memory_in_use    = array_memory_in_use  - vector_size
!$OMP END CRITICAL
        call mem_dealloc(arr%addr_loc)
      endif
      if (associated(arr%addr_p_arr))then
!$OMP CRITICAL
        vector_size = dble(size(arr%addr_p_arr))
        array_aux_deallocd_mem = array_aux_deallocd_mem + vector_size
        array_memory_in_use    = array_memory_in_use  - vector_size
!$OMP END CRITICAL
        call mem_dealloc(arr%addr_p_arr)
      endif
      if (associated(arr%lock_set))then
!$OMP CRITICAL
        vector_size = dble(size(arr%lock_set))
        array_aux_deallocd_mem = array_aux_deallocd_mem + vector_size
        array_memory_in_use    = array_memory_in_use  - vector_size
!$OMP END CRITICAL
        call mem_dealloc(arr%lock_set)
      endif
      if(array_aux_deallocd_mem > array_aux_allocd_mem) &
      &print *,"WARNING(arr_free_aux)more memory deallocated than allocated"
      call array_pdm_free_special_aux(arr)
  end subroutine arr_free_aux
  
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
    NULLIFY(arr%addr_loc)
    NULLIFY(arr%addr_p_arr)
  end subroutine array_nullify_pointers

end module  tensor_basic_module
