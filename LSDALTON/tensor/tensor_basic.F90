!> @file
!> basic tensor functionality
module tensor_basic_module
  use,intrinsic :: iso_c_binding, only:c_loc,c_f_pointer

  use tensor_parameters_and_counters
  use tensor_allocator
  use LSTIMING!,only:lstimer
!#ifdef VAR_MPI
!  use infpar_module
!  use lsmpi_type!, only: lsmpi_localwin_create_tensor_dp,&
!       & lsmpi_fence, lsmpi_win_free, lsmpi_barrier,&
!       & lsmpi_first_fence, lsmpi_last_fence
!#endif
  use lspdm_basic_module
  use tensor_type_def_module

  interface tensor_set_tdims
     module procedure tensor_set_tdims4,tensor_set_tdims8
  end interface tensor_set_tdims
  interface tensor_set_dims
     module procedure tensor_set_dims88,tensor_set_dims44
  end interface tensor_set_dims


  contains

    !> \brief to accurately account for mem on each node
    ! count mem allocated in each of the arrays many pointers
    !> \author Patrick Ettenhuber
    subroutine tensor_set_dims88(arr,dims,mode)
      implicit none
      type(tensor),intent(inout) :: arr
      integer(kind=tensor_long_int),intent(in) :: dims(*)
      integer(kind=tensor_long_int), optional :: mode
      include "tensor_set_dims.inc"
    end subroutine tensor_set_dims88
    !subroutine tensor_set_dims84(arr,dims,mode)
    !  implicit none
    !  type(tensor),intent(inout) :: arr
    !  integer(kind=tensor_long_int),intent(in) :: dims(*)
    !  integer(kind=tensor_standard_int), optional :: mode
    !  include "tensor_set_dims.inc"
    !end subroutine tensor_set_dims84
    subroutine tensor_set_dims44(arr,dims,mode)
      implicit none
      type(tensor),intent(inout) :: arr
      integer(kind=tensor_standard_int),intent(in) :: dims(*)
      integer(kind=tensor_standard_int), optional :: mode
      include "tensor_set_dims.inc"
    end subroutine tensor_set_dims44

    !> \brief to accurately account for mem on each node
    ! count mem allocated in each of the arrays many pointers
    !> \author Patrick Ettenhuber
    subroutine tensor_set_tdims4(arr,tdims,mode)
      implicit none
      type(tensor),intent(inout) :: arr
      integer(kind=tensor_standard_int),intent(in) :: tdims(arr%mode)
      integer, optional :: mode
      include "tensor_set_tdims.inc"
    end subroutine tensor_set_tdims4
    subroutine tensor_set_tdims8(arr,tdims,mode)
      implicit none
      type(tensor),intent(inout) :: arr
      integer(kind=tensor_long_int),intent(in) :: tdims(arr%mode)
      integer, optional :: mode
      include "tensor_set_tdims.inc"
    end subroutine tensor_set_tdims8
  
    subroutine tensor_init_lock_set(arr)
      implicit none
      type(tensor),intent(inout) :: arr
      integer(kind=tensor_long_int) :: vector_size

      if(.not.associated(arr%lock_set))then
        call tensor_alloc_mem(arr%lock_set,arr%ntiles)
        !$OMP CRITICAL
        vector_size = int(size(arr%lock_set)*tensor_standard_log,kind=tensor_long_int)
        tensor_counter_aux_a_mem = tensor_counter_aux_a_mem + vector_size
        tensor_counter_memory_in_use  = tensor_counter_memory_in_use  + vector_size
        !$OMP END CRITICAL
        arr%lock_set = .false.
      endif
      
    end subroutine tensor_init_lock_set

    !> \brief to accurately account for mem on each node
    ! count mem allocated in each of the arrays many pointers
    !> \author Patrick Ettenhuber
    subroutine tensor_set_ntpm(arr,ntpm,mode)
       implicit none
       type(tensor),intent(inout) :: arr
       integer(kind=tensor_standard_int),intent(in) :: ntpm(arr%mode)
       integer, optional :: mode
       integer(kind=tensor_long_int) :: vector_size
       integer :: i
       if (present(mode))then
          if((arr%mode/=0 .and. arr%mode/=mode).or.arr%mode==0)then
             print *,"mode",mode,"arr%mode",arr%mode      
             call lsquit("wrong use of tensor_set_ntpm",lspdm_errout)
          else
             arr%mode=mode
          endif
       endif
       vector_size = int(arr%mode*tensor_standard_int,kind=tensor_long_int)
       if (associated(arr%ntpm))then
          !$OMP CRITICAL
          tensor_counter_aux_f_mem = tensor_counter_aux_f_mem + vector_size
          tensor_counter_memory_in_use    = tensor_counter_memory_in_use  - vector_size
          !$OMP END CRITICAL
          call tensor_free_mem(arr%ntpm)
       endif
       if(.not.associated(arr%ntpm))then
          call tensor_alloc_mem(arr%ntpm,arr%mode)
          !$OMP CRITICAL
          tensor_counter_aux_a_mem = tensor_counter_aux_a_mem + vector_size
          tensor_counter_memory_in_use  = tensor_counter_memory_in_use  + vector_size
          !$OMP END CRITICAL
       endif
       do i=1,arr%mode
          arr%ntpm(i)=ntpm(i)
       enddo
    end subroutine tensor_set_ntpm

    subroutine tensor_set_addr(arr,addr,nnodes)
      implicit none
      type(tensor),intent(inout) :: arr
      integer,intent(in) :: addr(*)
      integer(kind=tensor_mpi_kind) :: nnodes
      integer(kind=tensor_standard_int) :: vector_size
      integer :: i

      vector_size = int(nnodes*tensor_standard_int,kind=tensor_long_int)
      if (associated(arr%addr_p_arr))then
         !$OMP CRITICAL
         tensor_counter_aux_f_mem     = tensor_counter_aux_f_mem      + vector_size
         tensor_counter_memory_in_use = tensor_counter_memory_in_use  - vector_size
         !$OMP END CRITICAL
         call tensor_free_mem(arr%addr_p_arr)
      endif
      if(.not.associated(arr%addr_p_arr))then
         call tensor_alloc_mem(arr%addr_p_arr,nnodes)
         !$OMP CRITICAL
         tensor_counter_aux_a_mem     = tensor_counter_aux_a_mem     + vector_size
         tensor_counter_memory_in_use = tensor_counter_memory_in_use + vector_size
         !$OMP END CRITICAL
      endif
      do i=1,nnodes
         arr%addr_p_arr(i)=addr(i)
      enddo
    end subroutine tensor_set_addr


  !> \brief Allocate memory for general arrays with memory statistics
  !> \author Patrick Ettenhuber adapted from Marcin Ziolkowski
    subroutine memory_allocate_tensor_dense(arr,bg)
      implicit none
      type(tensor),intent(inout) :: arr
      logical, intent(in) :: bg
      integer(kind=tensor_long_int) :: vector_size
      real(tensor_dp) :: tcpu1,twall1,tcpu2,twall2
      integer(kind=8) :: ne
      logical :: loc,parent

      call LSTIMER('START',tcpu1,twall1,lspdm_stdout)
      !call memory_deallocate_array(arr)
      if(associated(arr%elm1)) then
        call lsquit("ERROR(memory_allocate_array):array already initialized, please free first",lspdm_errout)
      endif
      vector_size = int((arr%nelms)*tensor_dp,kind=tensor_long_int)

      call tensor_alloc_mem(arr%elm1,arr%nelms,bg=bg)
      arr%bg_alloc = bg

      if( tensor_debug_mode )then
         arr%elm1 = 0.0E0_tensor_dp
      endif

!$OMP CRITICAL
      tensor_counter_dense_a_mem = tensor_counter_dense_a_mem + vector_size
      tensor_counter_memory_in_use = tensor_counter_memory_in_use + vector_size
      !tensor_counter_max_memory = max(tensor_counter_max_memory,tensor_counter_memory_in_use)
!$OMP END CRITICAL
      
      call assoc_ptr_arr(arr)

      call LSTIMER('START',tcpu2,twall2,lspdm_stdout)

    end subroutine memory_allocate_tensor_dense

  !> \brief deallocate dense part
  !> \author Patrick Ettenhuber
  !> \param array for which elm1 should be deallocated
    subroutine memory_deallocate_tensor_dense(arr)
      implicit none
      type(tensor) :: arr
      integer(kind=tensor_long_int) :: vector_size
      real(tensor_dp) :: dim1
      real(tensor_dp) :: tcpu1,twall1,tcpu2,twall2
      logical :: bg

      call LSTIMER('START',tcpu1,twall1,lspdm_stdout)

      if(associated(arr%elm1)) then

         vector_size = int(size(arr%elm1)*tensor_dp,kind=tensor_long_int)

         bg=arr%bg_alloc
         call tensor_free_mem(arr%elm1,bg=bg)

         !endif
         call deassoc_ptr_arr(arr)

!$OMP CRITICAL
         tensor_counter_dense_f_mem   = tensor_counter_dense_f_mem   + vector_size
         tensor_counter_memory_in_use = tensor_counter_memory_in_use - vector_size
!$OMP END CRITICAL
      end if

      call LSTIMER('START',tcpu2,twall2,lspdm_stdout)


    end subroutine memory_deallocate_tensor_dense





  !> \brief deallocate tiles and keep track of memory
  !> \author Patrick Ettenhuber
  !> \date September 2012
    subroutine memory_deallocate_tile(arr)
      implicit none
      type(tensor),intent(inout) :: arr
      integer(kind=tensor_long_int) :: vector_size
      real(tensor_dp) :: dim1,dim2,dim3,dim4
      real(tensor_dp) :: tcpu1,twall1,tcpu2,twall2
      integer :: i
      logical :: bg

      call LSTIMER('START',tcpu1,twall1,lspdm_stdout)

      do i=arr%nlti,1,-1
        if(associated(arr%ti(i)%t)) then
 
           if( .not. alloc_in_dummy )then

              vector_size = int(size(arr%ti(i)%t)*tensor_dp,kind=tensor_long_int)

              bg = arr%bg_alloc
              call tensor_free_mem(arr%ti(i)%t,arr%ti(i)%c,bg=bg)

              !$OMP CRITICAL
              tensor_counter_tiled_f_mem   = tensor_counter_tiled_f_mem   + vector_size
              tensor_counter_memory_in_use = tensor_counter_memory_in_use - vector_size
              !$OMP END CRITICAL
           endif

           arr%ti(i)%t => null()

           vector_size = int(size(arr%ti(i)%d)*tensor_int,kind=tensor_long_int)
           call tensor_free_mem(arr%ti(i)%d)
           !$OMP CRITICAL
           tensor_counter_aux_f_mem     = tensor_counter_aux_f_mem     + vector_size
           tensor_counter_memory_in_use = tensor_counter_memory_in_use - vector_size
           !$OMP END CRITICAL
        end if

      enddo

      vector_size = int(arr%nlti*tensor_bytes_per_tile,kind=tensor_long_int)
!$OMP CRITICAL
      tensor_counter_aux_f_mem     = tensor_counter_aux_f_mem     + vector_size
      tensor_counter_memory_in_use = tensor_counter_memory_in_use - vector_size
!$OMP END CRITICAL

      call tensor_free_mem(arr%ti)

      call LSTIMER('START',tcpu2,twall2,lspdm_stdout)


    end subroutine memory_deallocate_tile





    !> \brief the following collection of routines is a workaround the internal
    !fortran routines to be able to associate a one dimensional array with a
    !multidimensional pointer
    !> \author Patrick Ettenhuber
    !> \date September 2012
    subroutine assoc_ptr_arr(arr)
       implicit none
       type(tensor)::arr
       select case(arr%mode)
       case(2)
#ifdef VAR_PTR_RESHAPE
          arr%elm2(1:arr%dims(1),1:arr%dims(2)) => arr%elm1
#else
          call c_f_pointer(c_loc(arr%elm1(1)),arr%elm2,arr%dims)
#endif
       case(3)
#ifdef VAR_PTR_RESHAPE
          arr%elm3(1:arr%dims(1),1:arr%dims(2),1:arr%dims(3)) => arr%elm1
#else
          call c_f_pointer(c_loc(arr%elm1(1)),arr%elm3,arr%dims)
#endif
       case(4)
#ifdef VAR_PTR_RESHAPE
          arr%elm4(1:arr%dims(1),1:arr%dims(2),1:arr%dims(3),1:arr%dims(4)) => arr%elm1
#else
          call c_f_pointer(c_loc(arr%elm1(1)),arr%elm4,arr%dims)
#endif
       case(5)
#ifdef VAR_PTR_RESHAPE
          arr%elm5(1:arr%dims(1),1:arr%dims(2),1:arr%dims(3),1:arr%dims(4),1:arr%dims(5)) => arr%elm1
#else
          call c_f_pointer(c_loc(arr%elm1(1)),arr%elm5,arr%dims)
#endif
       case(6)
#ifdef VAR_PTR_RESHAPE
          arr%elm6(1:arr%dims(1),1:arr%dims(2),1:arr%dims(3),1:arr%dims(4),1:arr%dims(5),1:arr%dims(6)) => arr%elm1
#else
          call c_f_pointer(c_loc(arr%elm1(1)),arr%elm6,arr%dims)
#endif
       case(7)
#ifdef VAR_PTR_RESHAPE
          arr%elm7(1:arr%dims(1),1:arr%dims(2),1:arr%dims(3),1:arr%dims(4),1:arr%dims(5),1:arr%dims(6),1:arr%dims(7)) => arr%elm1
#else
          call c_f_pointer(c_loc(arr%elm1(1)),arr%elm7,arr%dims)
#endif
       case default
          return
       end select
    end subroutine assoc_ptr_arr
    
    !\brief deassociate all pointers again
    !\author Patrick Ettenhuber
    !\date September 2012
    subroutine deassoc_ptr_arr(arr)
      implicit none
      type(tensor)::arr
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
      write(lspdm_stdout,'(a,f12.2,a)') ' Allocated memory    : ',tensor_counter_dense_a_mem/2**20,' MB'
      write(lspdm_stdout,'(a,f12.2,a)') ' Deallocated memory  : ',tensor_counter_dense_f_mem/2**20,' MB'
      write(lspdm_stdout,'(a,f12.2,a)') ' Alloc-dealloc mem   : ', &
           (tensor_counter_dense_a_mem-tensor_counter_dense_f_mem)/2**20,' MB'
      write(lspdm_stdout,'(a,f12.2,a)') ' Memory in use       : ',tensor_get_total_mem()/2**20,' MB'
      write(lspdm_stdout,'(a,f12.2,a)') ' Max memory in use   : ',tensor_counter_max_hp_mem/2**20,' MB'

      write(lspdm_stdout,'(a,/,a)') '  Time ', &
           ' ======='

    end subroutine print_memory_statistics
  
  !> \brief Print currenly used memory and max allocated memory so far
  !> \author Marcin Ziolkowski, modified by Patrick Ettenhuber
  !> \param output File unit for output
  subroutine tensor_print_memory_currents(output,retour)
    implicit none
    !> selects the output
    integer, intent(in) :: output
    !> selects to redcue it on master for checking, or outut instead
    integer(kind=tensor_long_int), parameter :: nmeminfo = 9
    integer(kind=tensor_long_int),intent(inout),optional :: retour(nmeminfo)
    integer(kind=tensor_long_int) :: ret(nmeminfo)

    ret(1)=tensor_counter_dense_a_mem
    ret(2)=tensor_counter_dense_f_mem
    ret(3)=tensor_counter_tiled_a_mem
    ret(4)=tensor_counter_tiled_f_mem
    ret(5)=tensor_counter_aux_a_mem
    ret(6)=tensor_counter_aux_f_mem
    ret(7)=tensor_counter_memory_in_use_heap
    ret(8)=tensor_counter_max_hp_mem
    ret(9)=tensor_counter_max_bg_mem
    
    if(.not.present(retour))then
      write(*,'(a,i12,a)') ' Allocated memory for dense array   :',ret(1),' bytes'
      write(*,'(a,i12,a)') ' Deallocated memory for dense array :',ret(2),' bytes'
      write(*,'(a,i12,a)') ' Allocated memory for tiled array   :',ret(3),' bytes'
      write(*,'(a,i12,a)') ' Deallocated memory for tiled array :',ret(4),' bytes'
      write(*,'(a,i12,a)') ' Allocated aux memory for array     :',ret(5),' bytes'
      write(*,'(a,i12,a)') ' Dellocated aux memory for array    :',ret(6),' bytes'
      write(*,'(a,i12,a)') ' Memory in use for array            :',ret(7),' bytes'
      write(*,'(a,i12,a)') ' Max heap memory in use for array   :',ret(8),' bytes'
      write(*,'(a,i12,a)') ' Max bg memory in use for array     :',ret(9),' bytes'
    else
       retour = ret
    endif

  end subroutine tensor_print_memory_currents



  subroutine tensor_free_aux(arr)
     implicit none
     type(tensor),intent(inout) :: arr
     integer(kind=tensor_long_int) :: vector_size
     if (associated(arr%dims))then
        !$OMP CRITICAL
        vector_size = int(size(arr%dims)*tensor_int,kind=tensor_long_int)
        tensor_counter_aux_f_mem     = tensor_counter_aux_f_mem     + vector_size
        tensor_counter_memory_in_use = tensor_counter_memory_in_use - vector_size
        !$OMP END CRITICAL
        call tensor_free_mem(arr%dims)
     endif
     if (associated(arr%tdim))then
        !$OMP CRITICAL
        vector_size = int(size(arr%tdim)*tensor_standard_int,kind=tensor_long_int)
        tensor_counter_aux_f_mem     = tensor_counter_aux_f_mem     + vector_size
        tensor_counter_memory_in_use = tensor_counter_memory_in_use - vector_size
        !$OMP END CRITICAL
        call tensor_free_mem(arr%tdim)
     endif
     if (associated(arr%ntpm))then
        !$OMP CRITICAL
        vector_size = int(size(arr%ntpm)*tensor_standard_int,kind=tensor_long_int)
        tensor_counter_aux_f_mem     = tensor_counter_aux_f_mem     + vector_size
        tensor_counter_memory_in_use = tensor_counter_memory_in_use - vector_size
        !$OMP END CRITICAL
        call tensor_free_mem(arr%ntpm)
     endif
     if (associated(arr%addr_p_arr))then
        !$OMP CRITICAL
        vector_size = int(size(arr%addr_p_arr)*tensor_standard_int,kind=tensor_long_int)
        tensor_counter_aux_f_mem     = tensor_counter_aux_f_mem     + vector_size
        tensor_counter_memory_in_use = tensor_counter_memory_in_use - vector_size
        !$OMP END CRITICAL
        call tensor_free_mem(arr%addr_p_arr)
     endif
     if (associated(arr%lock_set))then
        !$OMP CRITICAL
        vector_size = int(size(arr%lock_set)*tensor_standard_log,kind=tensor_long_int)
        tensor_counter_aux_f_mem     = tensor_counter_aux_f_mem     + vector_size
        tensor_counter_memory_in_use = tensor_counter_memory_in_use - vector_size
        !$OMP END CRITICAL
        call tensor_free_mem(arr%lock_set)
     endif
     if(tensor_counter_aux_f_mem > tensor_counter_aux_a_mem) &
        &print *,"WARNING(tensor_free_aux)more memory deallocated than allocated"
     call tensor_pdm_free_special_aux(arr)
  end subroutine tensor_free_aux
  
  !> \author Patrick Ettenhuber adpted from Marcin Ziolkowski
  !> \date September 2012
  !> \brief free the array structure
  subroutine tensor_free_basic(arr)
    implicit none
    type(tensor), intent(inout) :: arr
    if(associated(arr%elm1))call memory_deallocate_tensor_dense(arr)
    if(associated(arr%ti))  call memory_deallocate_tile(arr)
    if(associated(arr%access_type)) deallocate( arr%access_type )
    call tensor_free_aux(arr)
  end subroutine tensor_free_basic

  !> \author Patrick Ettenhuber
  !> \date October 2014
  !> \brief set all values to invalid
  subroutine tensor_reset_value_defaults(arr)
     implicit none
     type(tensor) :: arr
     arr%local_addr  = -1
     arr%ntiles      = -1
     arr%nlti        = -1
     arr%tsize       = -1
     arr%offset      = -1
     arr%zeros       = .false.
     arr%initialized = .false.
     arr%nnod        = -1
     arr%comm        = -1
  end subroutine tensor_reset_value_defaults

  !> \author Patrick Ettenhuber
  !> \date September 2012
  !> \brief nullify all pointers in the array
  subroutine tensor_nullify_pointers(arr)
    implicit none
    type(tensor), intent(inout) :: arr
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
    NULLIFY(arr%access_type)
    arr%dummyc = c_null_ptr
  end subroutine tensor_nullify_pointers

end module  tensor_basic_module
