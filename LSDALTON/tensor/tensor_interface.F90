!> @file
!> Operations for general arrays
!> \author Patrick Ettenhuber

module tensor_interface_module

!`DIL backend (requires Fortran-2003/2008, MPI-3):
#ifdef COMPILER_UNDERSTANDS_FORTRAN_2003
#ifdef VAR_PTR_RESHAPE
#ifdef VAR_MPI
#define DIL_ACTIVE
#define DIL_DEBUG_ON
#endif
#endif
#endif

  ! Outside DEC directory
  use tensor_parameters_and_counters
  use tensor_mpi_interface_module
  use files!,only: lsopen,lsclose
  use LSTIMING!,only:lstimer
  use reorder_frontend_module
  use lspdm_tensor_operations_module
  use matrix_module
  use dec_workarounds_module
#ifdef DIL_ACTIVE
  use tensor_algebra_dil   !`DIL: Tensor Algebra
  public INTD,INTL         !integer sizes for DIL tensor algebra (default, long)
  public MAX_TENSOR_RANK   !max allowed tensor rank for DIL tensor algebra
  public DIL_TC_EACH       !parameter for <tensor_contract>: Each MPI process performs its own tensor contraction
  public DIL_TC_ALL        !parameter for <tensor_contract>: All MPI processes work on the same tensor contraction
  public DIL_ALLOC_NOT     !status "NOT ALLOCATED"
  public DIL_ALLOC_BASIC   !Fortran allocate() will be used for buffer allocation in <tensor_algebra_dil>
  public DIL_ALLOC_PINNED  !cudaMallocHost() will be used for buffer allocation in <tensor_algebra_dil>
  public DIL_ALLOC_MPI     !MPI_ALLOC_MEM() will be used for buffer allocation in <tensor_algebra_dil> (default for MPI)
  public DIL_ALLOC_EXT     !external buffer will be used in <tensor_algebra_dil>
  public DIL_CONS_OUT      !output for DIL messages
  public DIL_DEBUG         !DIL debugging switch
  public dil_tens_contr_t             !tensor contraction specification
  public subtens_t                    !subtensor (tensor slice) specification for Janus
  public dil_subtensor_set            !subtensor (tensor slice) setting method for Janus
  public dil_set_alloc_type           !set default memory allocation flags (BASIC,MPI_ALLOC,PINNED,etc.)
  public dil_clean_tens_contr         !clean the tensor contraction handle
  public dil_set_tens_contr_args      !set up an argument for a tensor contraction
  public dil_set_tens_contr_spec      !define the tensor contraction specification
  public dil_get_min_buf_size         !get the minimal buffer size needed to perform the tensor contraction
  public dil_prepare_buffer           !prepare a work buffer for a tensor contraction
  public dil_tensor_contract          !contract tensors (pipelined)
  public dil_tensor_contract_finalize !finalize a non-blocking tensor contraction
  public dil_debug_to_file_start      !start redirecting debugging information to a file
  public dil_debug_to_file_finish     !finish redirecting debugging information to a file
  public thread_wtime                 !OMP thread wall time
  public process_wtime                !MPI process wall time
  public dil_array_print              !print a whole (local) array or its part
  public dil_array_init               !initialize a local array
  public dil_tensor_init              !initialize a distributed tensor
  public dil_array_norm1              !compute the 1-norm of a local array
  public dil_tensor_norm1             !compute the 1-norm of a distributed tensor (blocking)
  public dil_tens_fetch_start         !tensor slice fetching for Janus (start)
  public dil_tens_fetch_finish_prep   !tensor slice fetching for Janus (finish)
  public dil_tens_prep_upload_start   !tensor slice uploading for Janus (start)
  public dil_tens_upload_finish       !tensor slice uploading for Janus  (finish)
  public dil_will_malloc_succeed      !tells whether a given malloc() request can succeed if issued
  public int2str                      !converts integers to strings
#endif

  !CALL THESE FUNCTION PRIOR TO ANY OTHER AND AS THE VERY LAST FUNCTIONS
  public tensor_initialize_interface, tensor_finalize_interface
  public tensor_set_comm, tensor_comm_null

  ! MODIFY THE BEHAVIOUR OF THE TENSOR LIB
  public tensor_set_mpi_msg_len              ! set the maximum message length of a call to MPI
  public tensor_set_always_sync_true         ! force synchronization after parallel operations
  public tensor_set_debug_mode_true          ! switch on debugging mode
  public tensor_set_dil_backend_true         ! switch on the dil backend for tensor contractions
  public tensor_set_dil_backend

  !This defines the public interface to the tensors
  !The tensor type itself
  public tensor
  !The different types of tensors in order to steer user level subroutines accordingly
  public TT_DENSE, TT_REPLICATED, TT_TILED, TT_TILED_DIST
  !The different access types in order to swith between master mediated and all at the same time accesses
  public AT_NO_PDM_ACCESS, AT_MASTER_ACCESS, AT_ALL_ACCESS
  ! USE THIS SIGNAL TO 
  public TENSOR_SLAVES_TO_SLAVE_ROUTINE_STD
  !Other parameters that may be useful for a user
  public alloc_in_dummy, TENSOR_MSG_LEN
  !Tensor timing
  public tensor_time_init
  ! User-level subroutines for the initialization
  public tensor_init, tensor_minit, tensor_ainit, tensor_free

  ! User-level subroutines for tensor operations
  public tensor_convert, print_norm, tensor_print_tile_norm
  public tensor_add, tensor_contract
  public tensor_transform_basis, tensor_ddot
  public tensor_reorder, tensor_cp_data, tensor_zero, tensor_scale, tensor_random
  public tensor_allocate_dense, tensor_deallocate_dense, tensor_hmul
  public tensor_print_norm_nrm


  ! PDM interface to the tensor structure
  public pdm_tensor_sync, new_group_reset_persistent_array
  public tensor_get_tile, tensor_put_tile, tensor_accumulate_tile
  public tensor_scatter, tensor_gather
  public tensor_lock_win, tensor_lock_wins, tensor_lock_local_wins
  public tensor_unlock_win, tensor_unlock_wins, tensor_unlock_local_wins
  public get_tensor_from_parr
  !subroutines changing the tensor distribution
  public tensor_sync_replicated
  public tensor_mv_dense2tiled, tensor_change_atype_to_d
  public tensor_cp_tiled2dense, tensor_change_atype_to_rep

  ! Special operations with tensors
  public tensor_extract_eos_indices, tensor_extract_decnp_indices
  public get_fragment_cc_energy_parallel, get_cc_energy_parallel
  public lspdm_get_combined_SingleDouble_amplitudes, get_info_for_mpi_get_and_reorder_t1 
  public get_rpa_energy_parallel, get_sosex_cont_parallel, get_starting_guess
  public precondition_doubles_parallel
  public tensor_dmul

  ! Only for testing and debugging
  public tensor_print_mem_info
  public lspdm_start_up_comm_procs, lspdm_shut_down_comm_procs

  ! Auxiliary functions on the user level
  public get_symm_tensor_segmenting_simple
  public tensor_get_ntpm, get_tile_dim
  public tensor_set_global_segment_length
  public check_if_new_instance_needed, find_free_pos_in_buf, find_tile_pos_in_buf
  public assoc_ptr_to_buf, lspdm_init_global_buffer, lspdm_free_global_buffer
  public tensor_flush_win

  private


  !> Number of created arrays
  integer(kind=long) :: ArraysCreated      = 0
  !> Number of destroyed arrays
  integer(kind=long) :: ArraysDestroyed    = 0
  !> Number of created arrays
  integer(kind=long) :: CreatedPDMArrays   = 0
  integer(kind=long) :: DestroyedPDMArrays = 0


  
  !> TIMINGS
  real(tensor_dp) :: tensor_time_init = 0


  !> convert arrays, the idea is for a general conversion only the interface
  !should be called
  interface tensor_convert
     module procedure tensor_convert_fort2tensor_wrapper1,&
        &tensor_convert_fort2tensor_wrapper2,tensor_convert_fort2tensor_wrapper3,&
        &tensor_convert_fort2tensor_wrapper4,tensor_convert_array22array,&
        &tensor_convert_tensor2fort_wrapper1,tensor_convert_tensor2fort_wrapper2,&
        &tensor_convert_tensor2fort_wrapper3,tensor_convert_tensor2fort_wrapper4
  end interface tensor_convert
  
  !> print norms of array, array2 array3, array4 and fortran arrays
  interface print_norm
     module procedure print_norm_fort_wrapper1_nrm,&
        &print_norm_fort_wrapper2_nrm,&
        &print_norm_fort_wrapper3_nrm,&
        &print_norm_fort_wrapper4_nrm,&
        &tensor_print_norm_nrm,&
        &array2_print_norm_nrm,&
        &array4_print_norm_nrm,&
        &matrix_print_norm_nrm,&
        &print_norm_fort_wrapper1_customprint,&
        &print_norm_fort_wrapper2_customprint,&
        &print_norm_fort_wrapper3_customprint,&
        &print_norm_fort_wrapper4_customprint,&
        &tensor_print_norm_customprint,&
        &array2_print_norm_customprint,&
        &array4_print_norm_customprint,&
        &print_norm_fort_nolen1_customprint,&
        &print_norm_fort_nolen2_customprint,&
        &print_norm_fort_nolen3_customprint,&
        &print_norm_fort_nolen4_customprint
  end interface print_norm


  interface tensor_add
     module procedure tensor_add_normal, tensor_add_arr2fullfort,tensor_add_fullfort2arr
  end interface tensor_add

  interface tensor_ainit
     module procedure tensor_ainit88,&
                     &tensor_ainit84,&
                     &tensor_ainit48,&
                     &tensor_ainit44
  end interface tensor_ainit

  !interface tensor_contract
  !  module procedure tensor_contract_pref
  !end interface tensor_contract


contains

  subroutine tensor_initialize_interface(comm,mem_ctr,pdm_slaves_signal)
     implicit none
     integer(kind=tensor_mpi_kind) :: comm
     !use an external counter for memory counting
     integer(kind=tensor_long_int), target, optional :: mem_ctr
     integer, intent(in), optional :: pdm_slaves_signal
     call tensor_set_comm(comm)
     if(present(mem_ctr)) call set_external_mem_ctr(mem_ctr)
     if(present(pdm_slaves_signal)) call set_signal_for_slaves(pdm_slaves_signal)
     call tensor_init_counters()
     call init_persistent_array()
  end subroutine tensor_initialize_interface
  subroutine tensor_finalize_interface()
     implicit none
     call free_persistent_array()
     call tensor_free_counters()
     if(associated(tensor_counter_ext_mem))call unset_external_mem_ctr()
  end subroutine tensor_finalize_interface
  subroutine tensor_set_mpi_msg_len(len)
     implicit none
     integer(kind=tensor_long_int) :: len
     TENSOR_MPI_MSG_LEN = len
  end subroutine tensor_set_mpi_msg_len

  subroutine tensor_set_global_segment_length(comm,seg_len)
     implicit none
     integer(kind=tensor_mpi_kind), intent(in) :: comm
     integer(kind=tensor_long_int), intent(in) :: seg_len
     integer(kind=tensor_mpi_kind) :: me, master
     integer(kind=tensor_long_int) :: seg
     me     = 0
     master = 0
     seg    = seg_len

#ifdef VAR_MPI
     call tensor_get_rank_for_comm(comm,me)

     if( me == 0 )then
        call pdm_tensor_sync(comm,JOB_SET_TENSOR_SEG_LENGTH)
     endif
     call tensor_mpi_bcast(seg,master,comm)
#endif

     if( seg<=0 )then
        call lsquit("ERROR(tensor_set_global_segment_length): invalid length",-1)
     endif

     print *,"SETTING LENGTH TO",seg
     tensor_segment_length_set  = .true.
     tensor_segment_length      = seg

  end subroutine tensor_set_global_segment_length

  subroutine tensor_set_debug_mode_true(comm,call_slaves)
     implicit none
     integer(kind=tensor_mpi_kind), intent(in) :: comm
     logical, intent(in) :: call_slaves
     integer(kind=tensor_mpi_kind) :: me
     me = 0
#ifdef VAR_MPI
     call tensor_get_rank_for_comm(comm,me)
     if( me == 0 .and. call_slaves )then
        call pdm_tensor_sync(comm,JOB_SET_TENSOR_DEBUG_TRUE)
     endif
#endif

     tensor_debug_mode  = .true.
     tensor_always_sync = .true.
  end subroutine tensor_set_debug_mode_true

  subroutine tensor_set_always_sync_true(comm,call_slaves)
     implicit none
     integer(kind=tensor_mpi_kind), intent(in) :: comm
     logical, intent(in) :: call_slaves
     integer(kind=tensor_mpi_kind) :: me
     me = 0
#ifdef VAR_MPI
     call tensor_get_rank_for_comm(comm,me)
     if( me == 0 .and. call_slaves )then
        call pdm_tensor_sync(comm,JOB_SET_TENSOR_ALWAYS_SYNC_TRUE)
     endif
#endif

     tensor_always_sync = .true.
  end subroutine tensor_set_always_sync_true

  subroutine tensor_set_dil_backend_true(comm,call_slaves)
     implicit none
     integer(kind=tensor_mpi_kind), intent(in) :: comm
     logical, intent(in) :: call_slaves
     integer(kind=tensor_mpi_kind) :: me
     me = 0
#ifdef VAR_MPI
     call tensor_get_rank_for_comm(comm,me)
     if( me == 0.and. call_slaves )then
        call pdm_tensor_sync(comm,JOB_SET_TENSOR_BACKEND_TRUE)
     endif
#endif
     tensor_contract_dil_backend = alloc_in_dummy !works only with MPI-3
  end subroutine tensor_set_dil_backend_true

  subroutine tensor_set_dil_backend(lv)
   implicit none
   logical, intent(in):: lv
   tensor_contract_dil_backend=(lv.and.alloc_in_dummy) !works only with MPI-3
  end subroutine tensor_set_dil_backend

  subroutine tensor_allocate_dense(T,bg,change)
     implicit none
     type(tensor), intent(inout) :: T
     logical, optional, intent(in) :: bg, change
     logical :: bg_int, change_int

     bg_int = .false.
     if(present(bg))bg_int = bg
     change_int = .false.
     if(present(change))change_int = change

     call memory_allocate_tensor_dense(T, bg_int)

     if(change_int)T%itype=TT_DENSE

  end subroutine tensor_allocate_dense


  subroutine copy_array(tensor_in,tensor_out,bg)
    implicit none
    type(tensor), intent(in) :: tensor_in
    type(tensor), intent(inout) :: tensor_out
    logical, intent(in), optional :: bg
    integer :: i
    logical :: bg_int

    bg_int = .false.
    if(present(bg))bg_int = bg

    tensor_out%mode = tensor_in%mode
    tensor_out%nlti = tensor_in%nlti
    tensor_out%tsize = tensor_in%tsize
    tensor_out%itype = tensor_in%itype
    tensor_out%nelms = tensor_in%nelms
    tensor_out%ntiles = tensor_in%ntiles
    tensor_out%access_type => tensor_in%access_type
    if(associated(tensor_in%dims))call tensor_set_dims(tensor_out,tensor_in%dims)
    if(associated(tensor_in%ntpm))call tensor_set_ntpm(tensor_out,tensor_in%ntpm,int(tensor_out%mode))
    if(associated(tensor_in%tdim))call tensor_set_tdims(tensor_out,tensor_in%tdim,int(tensor_out%mode))
    if(associated(tensor_in%addr_p_arr))call tensor_set_addr(tensor_out,int(tensor_in%addr_p_arr),&
       &size(tensor_in%addr_p_arr,kind=tensor_mpi_kind))
    !tensor_out%dims = tensor_in%dims
    !tensor_out%tdim = tensor_in%tdim
    !tensor_out%ntpm = tensor_in%ntpm
    if(associated(tensor_in%elm1))then
      call memory_allocate_tensor_dense(tensor_out,bg_int)
      tensor_out%elm1=tensor_in%elm1
    endif
    if(associated(tensor_in%ti))then
      call memory_allocate_tiles(tensor_out,bg_int)
      do i=1,tensor_in%nlti
        tensor_out%ti(i)%t=tensor_in%ti(i)%t
      enddo
    endif
  end subroutine copy_array


  ! x = (a *) x + b * y
  !> \brief add a scaled array to another array. The data may have different
  !distributions in the two arrays to be added
  !> \author Patrick Ettenhuber
  !> \date late 2012
  subroutine tensor_add_normal(x,b,y,a,order)
     implicit none
     !> array input, this is the result array with overwritten data
     type(tensor),intent(inout) :: x
     !> array to add
     type(tensor),intent(in) :: y
     !> scaling factor for array y
     real(tensor_dp),intent(in) :: b
     !> order the second array such that it fits the first
     integer, intent(in), optional :: order(x%mode)
     !> optional argument to scale x on the fly
     real(tensor_dp), intent(in), optional :: a
     real(tensor_dp),pointer :: buffer(:)
     real(tensor_dp) :: pre2
     integer :: ti,i,nel,o(x%mode)
     call time_start_phase( PHASE_WORK )

     pre2 = 1.0E0_tensor_dp
     if(present(a))pre2 = a

     if(x%mode/=y%mode)call lsquit("ERROR(tensor_add_normal): modes of arrays not compatible",-1)

     do i=1,x%mode
        if(present(order))then
           o(i) = order(i)
        else
           o(i) = i
        endif
        if(x%dims(i) /= y%dims(o(i)))call lsquit("ERROR(tensor_add_normal): dims of arrays not &
           &compatible (with the given order)",-1)
     enddo


     select case(x%itype)

     case(TT_DENSE,TT_REPLICATED)

        select case(y%itype)
        case(TT_DENSE,TT_REPLICATED)

           if( x%mode == 1 .or. .not.present(order))then
              if(present(a))then
                 if(a==0.0E0_tensor_dp)then
                    x%elm1 = 0.0E0_tensor_dp
                 else if(a/=1.0E0_tensor_dp)then
                    call dscal(int(x%nelms),a,x%elm1,1)
                 endif
              endif
              call daxpy(int(x%nelms),b,y%elm1,1,x%elm1,1)
           else
              select case(x%mode)
              case(2)
                 call array_reorder_2d(b,y%elm1,y%dims(1),y%dims(2),o,pre2,x%elm1)
              case(3)
                 call array_reorder_3d(b,y%elm1,y%dims(1),y%dims(2),y%dims(3),o,pre2,x%elm1)
              case(4)
                 call array_reorder_4d(b,y%elm1,y%dims(1),y%dims(2),y%dims(3),y%dims(4),o,pre2,x%elm1)
              case default
                 call lsquit("ERROR(tensor_add_normal): mode not implemented",-1)
              end select
           end if

           if(x%itype==TT_REPLICATED)call tensor_sync_replicated(x)

        case(TT_TILED_DIST)

           call tensor_alloc_mem(buffer,y%tsize)
           !TODO:IMPLEMENT MULTIPLE BUFFERING AND MOVE TO lspdm_tensor_operations!!!!!!
           do ti=1,y%ntiles
              call get_tile_dim(nel,y,ti)

              call time_start_phase( PHASE_COMM )
              call tensor_get_tile(y,int(ti,kind=tensor_standard_int),buffer,nel)
              call time_start_phase( PHASE_WORK )

              call tile_in_fort(b,buffer,ti,int(y%tdim),pre2,x%elm1,x%dims,int(x%mode),o)
           enddo
           call tensor_free_mem(buffer)

           call time_start_phase( PHASE_COMM )
           if(x%itype==TT_REPLICATED)call tensor_sync_replicated(x)
           call time_start_phase( PHASE_WORK )

        case default
           print *,x%itype,y%itype
           call lsquit("ERROR(tensor_add):not yet implemented y%itype 1",DECinfo%output)
        end select

     case(TT_TILED_DIST)

        select case(y%itype)
        case(TT_TILED_DIST)

           call tensor_add_par(pre2,x,b,y,o)

        case default
           print *,x%itype,y%itype
           call lsquit("ERROR(tensor_add):not yet implemented y%itype 2",DECinfo%output)
        end select

     case default
           print *,x%itype,y%itype
           call lsquit("ERROR(tensor_add_normal):not yet implemented x%itype",DECinfo%output)
     end select

     call time_start_phase( PHASE_WORK )
  end subroutine tensor_add_normal
  ! x = a * x + b * d * y[order]
  !> \brief add a scaled array to another array. The data may have different
  !distributions in the two arrays to be added
  !> \author Patrick Ettenhuber
  !> \date late 2012
  subroutine tensor_dmul(x,b,d,y,a,order)
     implicit none
     !> array input, this is the result array with overwritten data
     type(tensor),intent(inout) :: x
     !> array to add
     type(tensor),intent(in) :: y
     !> scaling factor for array y
     real(tensor_dp),intent(in) :: b,d(:)
     !> order the second array such that it fits the first
     integer, intent(in), optional :: order(x%mode)
     !> optional argument to scale x on the fly
     real(tensor_dp), intent(in), optional :: a
     real(tensor_dp),pointer :: buffer(:)
     real(tensor_dp) :: pre2
     integer :: ti,i,nel,o(x%mode),m,n
     call time_start_phase( PHASE_WORK )

     pre2 = 1.0E0_tensor_dp
     if(present(a))pre2 = a

     if(x%mode /= 2) call lsquit("ERROR(tensor_dmul): only implemented for mode 2 tensors",-1)

     if(x%mode/=y%mode)call lsquit("ERROR(tensor_dmul): modes of arrays not compatible",-1)

     do i=1,x%mode
        if(present(order))then
           o(i) = order(i)
        else
           o(i) = i
        endif
        if(x%dims(i) /= y%dims(o(i)))call lsquit("ERROR(tensor_dmul): dims of arrays not &
           &compatible (with the given order)",-1)
     enddo


     select case(x%itype)

     case(TT_DENSE,TT_REPLICATED)

        select case(y%itype)
        case(TT_DENSE,TT_REPLICATED)

           n = x%dims(1)
           m = x%dims(2)

           if(abs(pre2)<1.0E-15)then
              x%elm1 = 0.0E0_tensor_dp
           else if(abs(pre2-1.0E0_tensor_dp)>1.0E-15)then
              call dscal(x%nelms,pre2,x%elm1,1)
           endif

           if (o(1) == 1 .and. o(2) == 2) then

              do i=1,n
                 call daxpy(m,b*d(i),y%elm1(i),n,x%elm1(i),n)
              enddo

           else if (o(1)==2 .and. o(2)==1) then

              do i=1,m
                 call daxpy(n,b*d(i),y%elm1(n*(i-1)+1),1,x%elm1(i),m)
              enddo
           else
              call lsquit("ERROR(tensor_dmul): wrong order",-1)
           end if

           if(x%itype==TT_REPLICATED)call tensor_sync_replicated(x)

        case default
           print *,x%itype,y%itype
           call lsquit("ERROR(tensor_add):not yet implemented y%itype 1",DECinfo%output)
        end select

     case(TT_TILED_DIST)

        select case(y%itype)
        case(TT_TILED_DIST)

           call tensor_dmul_par(pre2,x,b,d,y,o)

        case default
           print *,x%itype,y%itype
           call lsquit("ERROR(tensor_add):not yet implemented y%itype 2",DECinfo%output)
        end select

     case default
           print *,x%itype,y%itype
           call lsquit("ERROR(tensor_dmul):not yet implemented x%itype",DECinfo%output)
     end select

     call time_start_phase( PHASE_WORK )
  end subroutine tensor_dmul

  subroutine tensor_transform_basis(U,nus,tens,whichU,t,maxtensmode,ntens,bg)
     implicit none
     !> specify the number of thensors that should be transformed
     integer, intent(in) :: ntens,nus,maxtensmode
     !list of which index of the trafo matrices to contract with the tensor index
     integer, intent(in) :: t(maxtensmode,ntens)
     ! list of which u to use to contract with the tensor index
     integer, intent(in) :: whichU(maxtensmode,ntens)
     !this contains the transformation matrices
     type(tensor), intent(in) :: U(nus)
     !this contains the tensors
     type(tensor), intent(in) :: tens(ntens)
     !use bg buf for temp alloc
     logical, intent(in), optional :: bg

     !internal variables
     integer :: itens, imode, it_mode, isort
     integer :: ord(maxtensmode)
     type(tensor) :: AUX, AUX1, AUX2
     real(tensor_dp), parameter :: p10 = 1.0E0_tensor_dp
     real(tensor_dp), parameter :: p00 = 0.0E0_tensor_dp

     do itens=1,ntens

        it_mode = tens(itens)%mode


        !initialize a tensor that has exactly the same type and distribution as
        !the one that should be tansformed
        call tensor_init(AUX, tens(itens)%dims, it_mode, &
           & pdm         = tens(itens)%access_type, &
           & tensor_type = tens(itens)%itype, &
           & tdims       = int(tens(itens)%tdim), &
           & fo          = int(tens(itens)%offset), &
           & bg          = bg ) 

        do imode = 1, it_mode

           do isort = 1, it_mode
              if(isort == imode)then
                 ord(isort) = 1
              else if(isort < imode )then
                 ord(isort) = isort+1
              else if(isort > imode )then
                 ord(isort) = isort
              endif
           enddo

           !USE ALIASES FOR THE ACTUAL CONTRACTION
           if( mod(imode,2) == 0)then
              AUX1 = AUX
              AUX2 = tens(itens)
           else
              AUX1 = tens(itens)
              AUX2 = AUX
           endif

           call tensor_contract(p10,U(whichU(imode,itens)),AUX1,[t(imode,itens)],&
              &[imode],1,p00,AUX2,ord(1:it_mode),force_sync=.true.)

        enddo

        call tensor_free(AUX)
     enddo

  end subroutine tensor_transform_basis


  ! x = x + b * y
  !> \brief add a scaled fortran-array to an array. The data may have arbitrary
  !distribution for the array
  !> \author Patrick Ettenhuber
  !> \date late 2012
  subroutine tensor_add_fullfort2arr(arrx,b,fortarry,order,wrk,iwrk)
    implicit none
    !> full fortan arra´y, this corresponds to y
    real(tensor_dp), intent(in) :: fortarry(*)
    !> scaling factor for fortran array
    real(tensor_dp), intent(in) :: b
    !> array which is overwritten
    type(tensor), intent(inout) :: arrx
    !> order of the fortran array with respect to the array
    integer, intent(in),optional :: order(arrx%mode)
    !> optinally workspace can be passed, the size is defined as iwrk
    integer(kind=8), intent(in),optional :: iwrk
    real(tensor_dp), intent(inout),optional :: wrk(*)
    integer :: o(arrx%mode)
    !> check if there is enough memory to send a full tile, this will die out
    integer :: i
    real(tensor_dp) :: MemFree,tilemem
    call time_start_phase( PHASE_WORK )

    do i=1,arrx%mode
      o(i) = i
    enddo
    if(present(order))o=order
    select case(arrx%itype)
      case(TT_DENSE)
        if(.not.present(order))then
          call daxpy(int(arrx%nelms),b,fortarry,1,arrx%elm1,1)
        else
          call lsquit("ERROR(tensor_add_fullfort2arr1):not implemented",-1)
        endif
      case(TT_TILED)
        call lsquit("ERROR(tensor_add_fullfort2arr):not implemented",-1)
      case(TT_TILED_DIST)
        call tensor_scatter(b,fortarry,1.0E0_tensor_dp,arrx,arrx%nelms,oo=order,wrk=wrk,iwrk=iwrk)
    end select

    call time_start_phase( PHASE_WORK )
  end subroutine tensor_add_fullfort2arr

  ! x = x + b * y
  !> \brief add a scaled array to a fortran- array. The data may have arbitrary
  !distribution for the array
  !> \author Patrick Ettenhuber
  !> \date late 2012
  subroutine tensor_add_arr2fullfort(fortarrx,b,arry,order,wrk,iwrk)
    implicit none
    !> full fortan arra´y, this corresponds to x and is overwritten
    real(tensor_dp), intent(inout) :: fortarrx(*)
    !> scaling factor for the array y
    real(tensor_dp), intent(in) :: b
    !> array to add --> y
    type(tensor), intent(in) :: arry
    !> order of the fortran array with respect to the array
    integer, intent(in),optional :: order(arry%mode)
    !> optinally workspace can be passed, the size is defined as iwrk
    integer(kind=8), intent(in),optional :: iwrk
    real(tensor_dp), intent(inout),optional :: wrk(*)
    !> check if there is enough memory to send a full tile, this will die out
    integer :: i
    real(tensor_dp) :: MemFree,tilemem
    call time_start_phase( PHASE_WORK )

    select case(arry%itype)
      case(TT_DENSE)
        call daxpy(int(arry%nelms),b,arry%elm1,1,fortarrx,1)
      case(TT_TILED)
        call lsquit("ERROR(tensor_add_fullfort2arr):not implemented",-1)
      !case(TT_TILED_DIST)
      !  !call add_tileddata2fort(arry,b,fortarrx,arry%nelms,.true.,order = order)
      !  call tensor_gather(b,arry,1.0E0_tensor_dp,fortarrx,arry%nelms,oo=order,wrk=wrk,iwrk=iwrk)
     case default
        call lsquit("ERROR(tensor_add_arr2fullfort) not implemented",-1)
    end select

    call time_start_phase( PHASE_WORK )
  end subroutine tensor_add_arr2fullfort

  !> \brief Hadamard product Cij = alpha*Aij*Bij+beta*Cij
  !> \author Thomas Kjaergaard
  !> \date 2015
  subroutine tensor_hmul(alpha,A,B,beta,C)
     implicit none
     !> array input, this is the result array with overwritten data
     type(tensor),intent(inout) :: C
     type(tensor),intent(in) :: A,B
     !> scaling factor for array C
     real(tensor_dp),intent(in) :: beta
     !> scaling factor for array A and B
     real(tensor_dp),intent(in) :: alpha
     call time_start_phase( PHASE_WORK )
     if(A%mode/=B%mode)call lsquit("ERROR(tensor_hmul_normal): modes of arrays not compatible",-1)
     if(A%mode/=C%mode)call lsquit("ERROR(tensor_hmul_normal): modes of arrays not compatible",-1)

     select case(C%itype)
        
     case(TT_TILED_DIST)
        
        select case(A%itype)

        case(TT_TILED_DIST)
           
           select case(B%itype)

           case(TT_TILED_DIST)

              call tensor_hmul_par(alpha,A,B,beta,C)
           case default
              print *,A%itype,B%itype,C%itype
              call lsquit("ERROR(tensor_hmul_normal):not yet implemented B%itype",DECinfo%output)
           end select
        case default
           print *,A%itype,B%itype,C%itype
           call lsquit("ERROR(tensor_hmul_normal):not yet implemented A%itype",DECinfo%output)
        end select
     case default
        print *,A%itype,B%itype,C%itype
        call lsquit("ERROR(tensor_hmul_normal):not yet implemented C%itype",DECinfo%output)
     end select 
    
     call time_start_phase( PHASE_WORK )
   end subroutine tensor_hmul


  !> \brief simple general tensor conraction of the type C = pre1 * A * B + pre2 * C
  !> \author Patrick Ettenhuber, Dmitry I. Lyakh (MPI-3 DIL backend)
  subroutine tensor_contract(pre1,A,B,m2cA,m2cB,nmodes2c,pre2,C,order,mem,wrk,iwrk,force_sync)
     implicit none
     real(tensor_dp), intent(in):: pre1,pre2                  !prefactors
     type(tensor), intent(in):: A,B                       !left and right tensors
     integer, intent(in):: nmodes2c                       !number of contracted dims
     integer, intent(in):: m2cA(nmodes2c), m2cB(nmodes2c) !contracted dims in left and right tensors
     type(tensor), intent(inout):: C                      !destination tensor
     integer, intent(inout):: order(C%mode)               !C dim x maps onto (A*B) dim order(x)
     real(tensor_dp), intent(in), optional:: mem              !???
     real(tensor_dp), intent(inout), optional:: wrk(:)        !external buffer
     integer(kind=long), intent(in), optional:: iwrk      !???
     logical, intent(in), optional:: force_sync           !???
     !internal variables (PETT)
     integer:: i,j,k
     logical:: contraction_mode
     integer:: rorder(C%mode)
#ifdef DIL_ACTIVE
     !internal variables (DIL)
     character(256):: tcs
     type(dil_tens_contr_t):: tch
     integer(INTL):: dil_mem
     integer(INTD):: i0,i1,i2,i3,tcl,errc,tcm(MAX_TENSOR_RANK*2)
     integer(INTD):: tens_rank,tens_dims(MAX_TENSOR_RANK),tens_bases(MAX_TENSOR_RANK)
     integer(INTD):: ddims(MAX_TENSOR_RANK),ldims(MAX_TENSOR_RANK),rdims(MAX_TENSOR_RANK)
     integer(INTD):: dbase(MAX_TENSOR_RANK),lbase(MAX_TENSOR_RANK),rbase(MAX_TENSOR_RANK)
#else
     integer(4):: i0,i1,i2,tcm(128)
#endif
     character(26), parameter:: elett='abcdefghijklmnopqrstuvwxyz'

     call time_start_phase( PHASE_WORK )

!Argument check:
     if( (A%mode-nmodes2c) + (B%mode-nmodes2c) /= C%mode) then
        call lsquit("ERROR(tensor_contract): invalid contraction pattern",-1)
     endif

     do i = 1,C%mode
        rorder(order(i)) = i
     enddo

     do i = 1,nmodes2c
        if(A%dims(m2cA(i))/=B%dims(m2cB(i))) then
           call lsquit("ERROR(tensor_contract): Contracted modes in A and B incompatible",-1)
        endif
     enddo

     i0 = 0
     i1 = 0
     k = 1
     do i = 1, A%mode
        i0 = i0 + 1
        contraction_mode=.false.
        do j=1,nmodes2c
!           contraction_mode = contraction_mode.or.(m2cA(j) == i)
           if(m2cA(j) == i) then
            contraction_mode = .true.
            i2=j
            exit
           endif
        enddo
        if(.not.contraction_mode) then
           i1 = i1 + 1
           tcm(i0) = i1
           if(A%dims(i) /= C%dims(rorder(k))) then
              call lsquit("ERROR(tensor_contract): Uncontracted modes in A and C incompatible",-1)
           endif
           k=k+1
        else
           tcm(i0) = -i2
        endif
     enddo
     do i = 1, B%mode
        i0 = i0 + 1
        contraction_mode=.false.
        do j=1,nmodes2c
!           contraction_mode = contraction_mode.or.(m2cB(j) == i)
           if(m2cB(j) == i) then
            contraction_mode = .true.
            i2=j
            exit
           endif
        enddo
        if(.not.contraction_mode) then
           i1 = i1 + 1
           tcm(i0) = i1
           if(B%dims(i) /= C%dims(rorder(k))) then
              call lsquit("ERROR(tensor_contract): Uncontracted modes in B and C incompatible",-1)
           endif
           k=k+1
        else
           tcm(i0) = -i2
        endif
     enddo

     if(k-1/=C%mode) then
        call lsquit("ERROR(tensor_contract): this should have resulted in a seg fault earlier",-1)
     endif

!Execute tensor contraction:
     if(.not.tensor_contract_dil_backend) then !PETT backend

      select case(A%itype)
      case(TT_DENSE,TT_REPLICATED)
        select case(B%itype)
        case(TT_DENSE,TT_REPLICATED)
           select case(C%itype)
           case(TT_DENSE,TT_REPLICATED)
              call tensor_contract_dense_simple(pre1,A,B,m2cA,m2cB,nmodes2c,pre2,C,order,&
                 &mem=mem,wrk=wrk,iwrk=iwrk)
              if(C%itype==TT_REPLICATED)call tensor_sync_replicated(C)
           case default
              call lsquit("ERROR(tensor_contract_simple): C%itype not implemented",-1)
           end select
        case default
           call lsquit("ERROR(tensor_contract_simple): B%itype not implemented",-1)
        end select
      case(TT_TILED_DIST)

        select case(B%itype)

        case(TT_DENSE,TT_REPLICATED)

           select case(C%itype)
           case(TT_TILED_DIST)
              call lspdm_tensor_contract_simple(pre1,A,B,m2cA,m2cB,nmodes2c,pre2,C,order,&
                 & mem=mem,wrk=wrk,iwrk=iwrk,force_sync=force_sync)
           case default
              call lsquit("error(tensor_contract_simple): c%itype not implemented",-1)
           end select

        case(TT_TILED_DIST)
           select case(C%itype)
           case(TT_TILED_DIST)
              call lspdm_tensor_contract_simple(pre1,A,B,m2cA,m2cB,nmodes2c,pre2,C,order,&
                 & mem=mem,wrk=wrk,iwrk=iwrk,force_sync=force_sync)
           case default
              call lsquit("error(tensor_contract_simple): c%itype not implemented",-1)
           end select

        case default
           call lsquit("ERROR(tensor_contract_simple): B%itype not implemented",-1)
        end select
      case default
        call lsquit("ERROR(tensor_contract_simple): A%itype not implemented",-1)
      end select

     else !`DIL backend (Fortran-2008 & MPI-3)
#ifdef DIL_ACTIVE
        !Get the symbolic tensor contraction pattern:
        tcs(1:2)='D('; tcl=2; i1=C%mode
        do i0=1,i1; tcs(tcl+1:tcl+2)=elett(i0:i0)//','; tcl=tcl+2; enddo
        if(tcs(tcl:tcl).ne.',') tcl=tcl+1; tcs(tcl:tcl)=')'
        tcs(tcl+1:tcl+4)='+=L('; tcl=tcl+4; i2=0
        do i0=1,A%mode
           i2=i2+1; i3=abs(tcm(i2))
           if(tcm(i2).gt.0) then !uncontracted index
             tcs(tcl+1:tcl+2)=elett(i3:i3)//','; tcl=tcl+2
           elseif(tcm(i2).lt.0) then !contracted index
             tcs(tcl+1:tcl+2)=elett(i1+i3:i1+i3)//','; tcl=tcl+2
           else
             call lsquit('ERROR(tensor_contract): DIL backend: symbolic part failed (A)!',-1)
           endif
        enddo
        if(tcs(tcl:tcl).ne.',') tcl=tcl+1; tcs(tcl:tcl)=')'
        tcs(tcl+1:tcl+3)='*R('; tcl=tcl+3
        do i0=1,B%mode
           i2=i2+1; i3=abs(tcm(i2))
           if(tcm(i2).gt.0) then !uncontracted index
             tcs(tcl+1:tcl+2)=elett(i3:i3)//','; tcl=tcl+2
           elseif(tcm(i2).lt.0) then !contracted index
             tcs(tcl+1:tcl+2)=elett(i1+i3:i1+i3)//','; tcl=tcl+2
           else
             call lsquit('ERROR(tensor_contract): DIL backend: symbolic part failed (B)!',-1)
           endif
        enddo
        if(tcs(tcl:tcl).ne.',') tcl=tcl+1; tcs(tcl:tcl)=')'
        if(DIL_DEBUG) write(*,*) '#DEBUG(DIL): symbolic: '//tcs(1:tcl)
        !Set tensor arguments:
        call dil_clean_tens_contr(tch)

        !Set the formal tensor contraction specification:

        !Perform the tensor contraction:

#else
        call lsquit('ERROR(tensor_contract_simple): DIL backend requires Fortran-2008 and MPI-3 at least!',-1)
#endif
     endif

     call time_start_phase( PHASE_WORK )
  end subroutine tensor_contract

  subroutine tensor_contract_dense_simple(pre1,A,B,m2cA,m2cB,nmodes2c,pre2,C,order,mem,wrk,iwrk)
     implicit none
     real(tensor_dp), intent(in)    :: pre1,pre2
     type(tensor), intent(in)    :: A,B
     integer, intent(in)        :: nmodes2c
     integer, intent(in)        :: m2cA(nmodes2c), m2cB(nmodes2c)
     type(tensor), intent(inout) :: C
     integer, intent(inout)     :: order(C%mode)
     real(tensor_dp), intent(in),    optional :: mem
     real(tensor_dp), intent(inout), target, optional :: wrk(:)
     integer(kind=long), intent(in),     optional :: iwrk
     !internal variables
     real(tensor_dp), pointer :: wA(:),  wB(:), wC(:)
     integer :: ordA(A%mode), ordB(B%mode), ro(C%mode), dims_product(C%mode)
     integer :: m_gemm, n_gemm, k_gemm, ldA,ldB
     integer :: i,j,k,l
     logical :: contraction_mode, use_wrk_space
     character :: tA, tB

     !in the case of simle matrices do not allocate buffer space
     ! TODO: this can be generalized for cases where the modes to contract are
     ! the outer modes of a tensor
     if(A%mode == 2 .and. B%mode == 2 .and. C%mode == 2)then

        m_gemm = C%dims(1)
        n_gemm = C%dims(2)
        wC => C%elm1

        if(order(1) == 1 .and. order(2) == 2)then

           select case (m2cA(1))
           case(1)
              tA = 't'
              k_gemm = A%dims(1)
           case(2)
              tA = 'n'
              k_gemm = A%dims(2)
           case default
              call lsquit("ERROR(tensor_contract_dense_simple): undefined contraction mode A",-1)
           end select
           select case (m2cB(1))
           case(1)
              tB = 'n'
           case(2)
              tB = 't'
           case default
              call lsquit("ERROR(tensor_contract_dense_simple): undefined contraction mode B",-1)
           end select


           wA => A%elm1
           ldA = A%dims(1)
           wB => B%elm1
           ldB = B%dims(1)


        else if(order(1) == 2 .and. order(2) == 1)then

           select case (m2cA(1))
           case(1)
              tA = 'n'
              k_gemm = A%dims(1)
           case(2)
              tA = 't'
              k_gemm = A%dims(2)
           case default
              call lsquit("ERROR(tensor_contract_dense_simple): undefined contraction mode A",-1)
           end select
           select case (m2cB(1))
           case(1)
              tB = 't'
           case(2)
              tB = 'n'
           case default
              call lsquit("ERROR(tensor_contract_dense_simple): undefined contraction mode B",-1)
           end select

           wA => B%elm1
           ldA = B%dims(1)
           wB => A%elm1
           ldB = A%dims(1)

        endif

        call dgemm(tA,tB,m_gemm,n_gemm,k_gemm,pre1,wA,ldA,wB,ldB,pre2,wC,m_gemm)

        wA => null()
        wB => null()
        wC => null()

     !GENERAL TENSOR CONTRACTION
     else

        do i = 1, C%mode
           ro(order(i)) = i
        enddo

        if(present(mem))then
           !use provided memory information to allcate space
           if(mem < (((A%nelms + B%nelms + C%nelms) * 8.0 )/ 1024.0**3))then
              print *,"WARNING(tensor_contract_dense_simple): too little memory &
                 &given, will try to allocate needed space anyways"
           endif
           use_wrk_space = .false.
        else if(present(wrk).and.present(iwrk))then
           !just assoctiate pointers to work space provided
           if(A%nelms + B%nelms + C%nelms > iwrk)then
              print *,"WARNING(tensor_contract_dense_simple): too small work space &
                 &given, will try to allocate needed space"
              use_wrk_space = .false.
           else
              use_wrk_space = .true.
           endif
        else
           use_wrk_space = .false.
        endif

        if(use_wrk_space)then
           wA => wrk(1:A%nelms)
           wB => wrk(A%nelms + 1 : A%nelms + B%nelms )
           wC => wrk(A%nelms + B%nelms + 1 : A%nelms + B%nelms + C%nelms )
        else
           call tensor_alloc_mem(wA,A%nelms)
           call tensor_alloc_mem(wB,B%nelms)
           call tensor_alloc_mem(wC,C%nelms)
        endif

        m_gemm = 1
        n_gemm = 1
        !get the uncontracted mode indices of the A and B arrays
        k = 1
        do i = 1, A%mode
           contraction_mode=.false.
           do j=1,nmodes2c
              contraction_mode = contraction_mode.or.(m2cA(j) == i)
           enddo
           if(.not.contraction_mode)then
              ordA(k) = i
              m_gemm  = m_gemm * C%dims(ro(k))
              k=k+1
           endif
        enddo

        if(k-1/=A%mode-nmodes2c)then
           call lsquit("ERROR(tensor_contract_dense_simple): something wrong in ordering",-1)
        endif

        k_gemm = 1
        do i = 1,nmodes2c
           ordA(k-1+i) = m2cA(i)
           ordB(i)     = m2cB(i)
           k_gemm      = k_gemm * A%dims(m2cA(i))
        end do

        l = 1
        do i = 1, B%mode
           contraction_mode=.false.
           do j=1,nmodes2c
              contraction_mode = contraction_mode.or.(m2cB(j) == i)
           enddo
           if(.not.contraction_mode)then
              ordB(nmodes2c+l) = i
              n_gemm           = n_gemm * C%dims(ro(k))
              k=k+1
              l=l+1
           endif
        enddo

        select case (A%mode)
        case(2)
           call array_reorder_2d(1.0E0_tensor_dp,A%elm1,A%dims(1),A%dims(2),ordA,0.0E0_tensor_dp,wA)
        case(3)
           call array_reorder_3d(1.0E0_tensor_dp,A%elm1,A%dims(1),A%dims(2),A%dims(3),ordA,0.0E0_tensor_dp,wA)
        case(4)
           call array_reorder_4d(1.0E0_tensor_dp,A%elm1,A%dims(1),A%dims(2),A%dims(3),A%dims(4),ordA,0.0E0_tensor_dp,wA)
        case default
           call lsquit("ERROR(tensor_contract_dense_simple): sorting A not implemented",-1)
        end select

        select case (B%mode)
        case(2)
           call array_reorder_2d(1.0E0_tensor_dp,B%elm1,B%dims(1),B%dims(2),ordB,0.0E0_tensor_dp,wB)
        case(3)
           call array_reorder_3d(1.0E0_tensor_dp,B%elm1,B%dims(1),B%dims(2),B%dims(3),ordB,0.0E0_tensor_dp,wB)
        case(4)
           call array_reorder_4d(1.0E0_tensor_dp,B%elm1,B%dims(1),B%dims(2),B%dims(3),B%dims(4),ordB,0.0E0_tensor_dp,wB)
        case default
           call lsquit("ERROR(tensor_contract_dense_simple): sorting B not implemented",-1)
        end select


        call dgemm('n','n',m_gemm,n_gemm,k_gemm,1.0E0_tensor_dp,wA,m_gemm,wB,k_gemm,0.0E0_tensor_dp,wC,m_gemm)


        do i=1,C%mode
           dims_product(i) = C%dims(ro(i))
        enddo

        !ADD THE FINALIZED TILE TO THE LOCAL TILE IN THE CORRECT ORDER
        select case (C%mode)
        case(2)
           call array_reorder_2d(pre1,wC,dims_product(1),dims_product(2),order,pre2,C%elm1)
        case(3)
           call array_reorder_3d(pre1,wC,dims_product(1),dims_product(2),dims_product(3),order,pre2,C%elm1)
        case(4)
           call array_reorder_4d(pre1,wC,dims_product(1),dims_product(2),dims_product(3),dims_product(4),order,pre2,C%elm1)
        case default
           call lsquit("ERROR(tensor_contract_dense_simple): sorting C not implemented",-1)
        end select

        if(use_wrk_space)then
           wA => null()
           wB => null()
           wC => null()
        else
           call tensor_free_mem(wA)
           call tensor_free_mem(wB)
           call tensor_free_mem(wC)
        endif
     endif
  end subroutine tensor_contract_dense_simple


  !> \brief array dotprduct, the arrays may have different distributions
  !> \author Patrick Ettenhuber
  !> \date some time at the end 2012
  function tensor_ddot(arr1,arr2,opt_par) result(res)
    implicit none
    !> the two arrays to calculate the dotproduct from
    type(tensor),intent(in) :: arr1,arr2
    !> optional integer specifying on which node the result should be stored
    integer,optional,intent(in) :: opt_par
    real(tensor_dp) :: res
    integer :: dest
    real(tensor_dp), external :: ddot
    
    if(arr1%nelms/=arr2%nelms)then
      call lsquit("ERROR(tensor_ddot):operation not defined for arrays with&
      & different nelms",DECinfo%output)
    endif

    !get the destination of the contraction
    dest = -1
    if(arr1%access_type==AT_MASTER_ACCESS)dest=0
    if(present(opt_par))dest=opt_par

    select case(arr1%itype)
    case(TT_DENSE)
      select case(arr2%itype)
      case(TT_DENSE)
        res=ddot(int(arr1%nelms),arr1%elm1,1,arr2%elm1,1)
      case default
        call lsquit("ERROR(tensor_ddot):operation not yet&
        & implemented",DECinfo%output)
      end select
            
    case(TT_TILED_DIST)
      select case(arr2%itype)
      case(TT_TILED_DIST)
        res=tensor_ddot_par(arr1,arr2,dest)
      case default
        call lsquit("ERROR(tensor_ddot):operation not yet&
        & implemented",DECinfo%output)
      end select
    case default
      call lsquit("ERROR(tensor_ddot):operation not yet&
      & implemented",DECinfo%output)
    end select

  end function tensor_ddot

  !> \brief Extract EOS indices from array4 for both occupied and virtual partitioning schemes:
  !> 1. tensor_occEOS: The occupied orbitals not assigned to the central atom are removed while
  !>                the virtual indices are unchanged.
  !> 2. tensor_virtEOS: The virtual orbitals not assigned to the central atom are removed while
  !>                 the occupied indices are unchanged.
  !> \author Patrick Ettenhuber, adapted from Kasper Kristensen
  !> \date August 2014
  subroutine tensor_extract_eos_indices(tensor_orig,MyFragment,tensor_occEOS,tensor_virtEOS)


    implicit none
    !> Array where occupied EOS indices are extracted
    type(tensor),intent(inout),optional :: tensor_occEOS
    !> Array where virtual EOS indices are extracted
    type(tensor),intent(inout),optional :: tensor_virtEOS
    !> Original array with AOS fragment indices for both occ and virt spaces
    type(tensor),intent(in) :: tensor_orig
    !> Atomic fragment
    type(decfrag),intent(inout) :: MyFragment
    integer :: nocc, nvirt

    ! Number of occ and virt orbitals on central atom in fragment
    nocc  = MyFragment%noccEOS
    nvirt = MyFragment%nvirtEOS

    ! Extract virtual EOS indices and leave occupied indices untouched
    ! ****************************************************************

    if(present(tensor_virtEOS))then
       call tensor_extract_eos_indices_virt(tensor_virtEOS,tensor_orig,&
          & nvirt,MyFragment%idxu(1:nvirt))
    endif

    ! Extract occupied EOS indices and leave virtual indices untouched
    ! ****************************************************************

    if(present(tensor_occEOS))then
       call tensor_extract_eos_indices_occ(tensor_occEOS,tensor_orig,&
          & nocc, MyFragment%idxo(1:nocc))
    endif

  end subroutine tensor_extract_eos_indices


  subroutine tensor_extract_eos_indices_virt(Arr,tensor_full,nEOS,EOS_idx)

     implicit none
     !> Array output where EOS indices are extracted
     type(tensor),intent(inout) :: Arr
     !> Original array
     type(tensor),intent(in) :: tensor_full
     !> Number of EOS indices
     integer,intent(in) :: nEOS
     !> List of EOS indices in the total (EOS+buffer) list of orbitals
     integer, dimension(nEOS),intent(in) :: EOS_idx
     integer :: nocc,nvirt,i,a,b,j,ax,bx
     integer, dimension(4) :: new_dims

     ! Initialize stuff
     ! ****************
     nocc     = tensor_full%dims(2)   ! Total number of occupied orbitals
     nvirt    = tensor_full%dims(1)   ! Total number of virtual orbitals
     new_dims = [nEOS,nocc,nEOS,nocc] ! nEOS=Number of virtual EOS orbitals


     ! Sanity checks
     ! *************

     ! 1. Positive number of orbitals
     if( (nocc<1) .or. (nvirt<1) ) then
        write(DECinfo%output,*) 'nocc = ', nocc
        write(DECinfo%output,*) 'nvirt = ', nvirt
        call lsquit('tensor_extract_eos_indices_virt: &
           & Negative or zero number of orbitals!',DECinfo%output)
     end if

     ! 2. Array structure is (virt,occ,virt,occ)
     if( (nvirt/=tensor_full%dims(3)) .or. (nocc/=tensor_full%dims(4)) ) then
        write(DECinfo%output,*) 'tensor_full%dims(1) = ', tensor_full%dims(1)
        write(DECinfo%output,*) 'tensor_full%dims(2) = ', tensor_full%dims(2)
        write(DECinfo%output,*) 'tensor_full%dims(3) = ', tensor_full%dims(3)
        write(DECinfo%output,*) 'tensor_full%dims(4) = ', tensor_full%dims(4)
        call lsquit('tensor_extract_eos_indices_virt: &
           & Arr dimensions does not match (virt,occ,virt,occ) structure!',DECinfo%output)
     end if

     ! 3. EOS dimension must be smaller than (or equal to) total number of virt orbitals
     if(nEOS > nvirt) then
        write(DECinfo%output,*) 'nvirt = ', nvirt
        write(DECinfo%output,*) 'nEOS  = ', nEOS
        call lsquit('array4_extract_eos_indices_virt_memory: &
           & Number of EOS orbitals must be smaller than (or equal to) total number of &
           & virtual orbitals!',DECinfo%output)
     end if

     ! 4. EOS indices must not exceed total number of virtual orbitals
     do i=1,nEOS
        if(EOS_idx(i) > nvirt) then
           write(DECinfo%output,'(a,i6,a)') 'EOS index number ', i, ' is larger than nvirt!'
           write(DECinfo%output,*) 'nvirt   = ', nvirt
           write(DECinfo%output,*) 'EOS_idx = ', EOS_idx(i)
           call lsquit('array4_extract_eos_indices_virt_memory: &
              & EOS index value larger than nvirt!',DECinfo%output)
        end if
     end do

     ! Extract virtual EOS indices and store in Arr
     ! ********************************************

     ! Initiate Arr with new dimensions (nvirt_EOS,nocc,nvirt_EOS,nocc) on the
     ! local node, as the tensor will be small enough to store locally
     call tensor_init(Arr,new_dims,4)

     select case ( tensor_full%itype )
     case( TT_DENSE, TT_REPLICATED )

        ! Set Arr equal to the EOS indices of the original Arr array (tensor_full)
        do j=1,nocc
           do b=1,nEOS
              bx=EOS_idx(b)
              do i=1,nocc
                 do a=1,nEOS
                    ax=EOS_idx(a)
                    Arr%elm4(a,i,b,j) = tensor_full%elm4(ax,i,bx,j)
                 end do
              end do
           end do
        end do

     case( TT_TILED_DIST )

        call tensor_zero(Arr)

        call lspdm_extract_eos_indices_virt(Arr,tensor_full,nEOS,EOS_idx)

     case default
        call lsquit("ERROR(tensor_extract_eos_indices_virt): NO PDM version implemented yet",-1)
     end select


  end subroutine tensor_extract_eos_indices_virt

  subroutine tensor_extract_eos_indices_occ(Arr,tensor_full,nEOS,EOS_idx)
     implicit none
     !> Array where EOS indices where are extracted
     type(tensor),intent(inout) :: Arr
     !> Original array
     type(tensor),intent(in) :: tensor_full
     !> Number of EOS indices
     integer,intent(in) :: nEOS
     !> List of EOS indices in the total (EOS+buffer) list of orbitals
     integer, dimension(nEOS),intent(in) :: EOS_idx
     integer :: nocc,nvirt,i,a,b,j,ix,jx
     integer, dimension(4) :: new_dims

     ! Initialize stuff
     ! ****************
     nocc     = tensor_full%dims(2)  ! Total number of occupied orbitals
     nvirt    = tensor_full%dims(1)  ! Total number of virtual orbitals
     new_dims = [nvirt,nEOS,nvirt,nEOS] ! nEOS=Number of occupied EOS orbitals

     ! Sanity checks
     ! *************
     if( tensor_full%mode /= 4)then
        call lsquit("ERROR(tensor_extract_eos_indices_occ): wrong mode of tensor_full",-1)
     endif

     ! 1. Positive number of orbitals
     if( (nocc<1) .or. (nvirt<1) ) then
        write(DECinfo%output,*) 'nocc = ', nocc
        write(DECinfo%output,*) 'nvirt = ', nvirt
        call lsquit('tensor_extract_eos_indices_occ: &
           & Negative or zero number of orbitals!',DECinfo%output)
     end if

     ! 2. Array structure is (virt,occ,virt,occ)
     if( (nvirt/=tensor_full%dims(3)) .or. (nocc/=tensor_full%dims(4)) ) then
        write(DECinfo%output,*) 'tensor_full%dims(1) = ', tensor_full%dims(1)
        write(DECinfo%output,*) 'tensor_full%dims(2) = ', tensor_full%dims(2)
        write(DECinfo%output,*) 'tensor_full%dims(3) = ', tensor_full%dims(3)
        write(DECinfo%output,*) 'tensor_full%dims(4) = ', tensor_full%dims(4)
        call lsquit('tensor_extract_eos_indices_occ: &
           & Arr dimensions does not match (virt,occ,virt,occ) structure!',DECinfo%output)
     end if

     ! 3. EOS dimension must be smaller than (or equal to) total number of occ orbitals
     if(nEOS > nocc) then
        write(DECinfo%output,*) 'nocc = ', nocc
        write(DECinfo%output,*) 'nEOS = ', nEOS
        call lsquit('tensor_extract_eos_indices_occ: &
           & Number of EOS orbitals must be smaller than (or equal to) total number of &
           & occupied orbitals!',DECinfo%output)
     end if

     ! 4. EOS indices must not exceed total number of occupied orbitals
     do i=1,nEOS
        if(EOS_idx(i) > nocc) then
           write(DECinfo%output,'(a,i6,a)') 'EOS index number ', i, ' is larger than nocc!'
           write(DECinfo%output,*) 'nocc = ', nocc
           write(DECinfo%output,*) 'EOS_idx = ', EOS_idx(i)
           call lsquit('tensor_extract_eos_indices_occ: &
              & EOS index value larger than nocc!',DECinfo%output)
        end if
     end do


     ! Extract occupied EOS indices and store in Arr
     ! *********************************************

     ! Initiate Arr with new dimensions (nvirt,nocc_EOS,nvirt,nocc_EOS)
     call tensor_init(Arr,new_dims,4)

     select case( tensor_full%itype )
     case( TT_DENSE, TT_REPLICATED )

        ! Set Arr equal to the EOS indices of the original Arr array (tensor_full)
        do j=1,nEOS
           jx=EOS_idx(j)
           do b=1,nvirt
              do i=1,nEOS
                 ix=EOS_idx(i)
                 do a=1,nvirt
                    Arr%elm4(a,i,b,j) = tensor_full%elm4(a,ix,b,jx)
                 end do
              end do
           end do
        end do

     case( TT_TILED_DIST )

        call tensor_zero(Arr)

        call lspdm_extract_eos_indices_occ(Arr,tensor_full,nEOS,EOS_idx)

     case default
        call lsquit("ERROR(tensor_extract_eos_indices_occ): NO PDM version implemented yet",-1)
     end select


  end subroutine tensor_extract_eos_indices_occ


  !> Purpose: Extract energy indices for DECNP calculation, based on Patrick routine
  !
  !> Author:  Pablo Baudin
  !> Date:    Feb. 2015
  subroutine tensor_extract_decnp_indices(tensor_full,myfragment,ArrOcc,ArrVir)

     implicit none

     !> Original array
     type(tensor),intent(in) :: tensor_full
     !> Atomic fragment
     type(decfrag), target, intent(inout) :: MyFragment
     !> Array where EOS indices where are extracted
     type(tensor),intent(inout) :: ArrOcc, ArrVir

     !> Number of EOS indices
     integer :: nEOS
     !> List of EOS indices in the total (EOS+buffer) list of orbitals
     integer, pointer :: EOS_idx(:)
     integer :: nocc,nvirt,i,a,b,j,ix,ax
     integer, dimension(4) :: new_dims

     !---------------------------------------------------------------------
     !                 EXTRACT OCCUPIED PARTITIONING ARRAY
     !---------------------------------------------------------------------

     ! Initialize stuff
     ! ****************
     nEOS     = myfragment%noccEOS
     EOS_idx  => myFragment%idxo(1:nEOS)
     nocc     = tensor_full%dims(2)     ! Total number of occupied orbitals
     nvirt    = tensor_full%dims(1)     ! Total number of virtual orbitals
     new_dims = [nvirt,nEOS,nvirt,nocc] ! nEOS=Number of occupied EOS orbitals

     ! Sanity checks
     ! *************
     if( tensor_full%mode /= 4)then
        call lsquit("ERROR(tensor_extract_decnp_indices): wrong mode of tensor_full",-1)
     endif

     ! 1. Positive number of orbitals
     if( (nocc<1) .or. (nvirt<1) ) then
        write(DECinfo%output,*) 'nocc = ', nocc
        write(DECinfo%output,*) 'nvirt = ', nvirt
        call lsquit('tensor_extract_decnp_indices: &
           & Negative or zero number of orbitals!',DECinfo%output)
     end if

     ! 2. Array structure is (virt,occ,virt,occ)
     if( (nvirt/=tensor_full%dims(3)) .or. (nocc/=tensor_full%dims(4)) ) then
        write(DECinfo%output,*) 'tensor_full%dims(1) = ', tensor_full%dims(1)
        write(DECinfo%output,*) 'tensor_full%dims(2) = ', tensor_full%dims(2)
        write(DECinfo%output,*) 'tensor_full%dims(3) = ', tensor_full%dims(3)
        write(DECinfo%output,*) 'tensor_full%dims(4) = ', tensor_full%dims(4)
        call lsquit('tensor_extract_decnp_indices: &
           & ArrOcc dimensions does not match (virt,occ,virt,occ) structure!',DECinfo%output)
     end if

     ! 3. EOS dimension must be smaller than (or equal to) total number of occ orbitals
     if(nEOS > nocc) then
        write(DECinfo%output,*) 'nocc = ', nocc
        write(DECinfo%output,*) 'nEOS = ', nEOS
        call lsquit('tensor_extract_decnp_indices: &
           & Number of EOS orbitals must be smaller than (or equal to) total number of &
           & occupied orbitals!',DECinfo%output)
     end if

     ! 4. EOS indices must not exceed total number of occupied orbitals
     do i=1,nEOS
        if(EOS_idx(i) > nocc) then
           write(DECinfo%output,'(a,i6,a)') 'EOS index number ', i, ' is larger than nocc!'
           write(DECinfo%output,*) 'nocc = ', nocc
           write(DECinfo%output,*) 'EOS_idx = ', EOS_idx(i)
           call lsquit('tensor_extract_decnp_indices: &
              & EOS index value larger than nocc!',DECinfo%output)
        end if
     end do


     ! Extract occupied EOS indices and store in ArrOcc
     ! ************************************************

     ! Initiate ArrOcc with new dimensions (nvirt,nocc,nvirt,nocc_EOS)
     call tensor_init(ArrOcc,new_dims,4)

     select case( tensor_full%itype )
     case( TT_DENSE, TT_REPLICATED )

        ! Set ArrOcc equal to the EOS indices of the original ArrOcc array (tensor_full)
        do j=1,nocc
           do b=1,nvirt
              do i=1,nEOS
                 ix=EOS_idx(i)
                 do a=1,nvirt
                    ArrOcc%elm4(a,i,b,j) = tensor_full%elm4(a,ix,b,j)
                 end do
              end do
           end do
        end do

     case( TT_TILED_DIST )

        call tensor_zero(ArrOcc)

        call lspdm_extract_decnp_indices_occ(ArrOcc,tensor_full,nEOS,EOS_idx)

     case default
        call lsquit("ERROR(tensor_extract_decnp_indices): NO PDM version implemented yet",-1)
     end select

     !---------------------------------------------------------------------
     !                 EXTRACT VIRTUAL PARTITIONING ARRAY
     !---------------------------------------------------------------------

     ! Initialize stuff
     ! ****************
     nEOS     = myfragment%nvirtEOS
     EOS_idx  => myFragment%idxu(1:nEOS)
     new_dims = [nEOS,nocc,nvirt,nocc]  ! nEOS=Number of virtual EOS orbitals


     ! Sanity checks
     ! *************

     ! 5. EOS dimension must be smaller than (or equal to) total number of virt orbitals
     if(nEOS > nvirt) then
        write(DECinfo%output,*) 'nvirt = ', nvirt
        write(DECinfo%output,*) 'nEOS  = ', nEOS
        call lsquit('tensor_extract_decnp_indices: &
           & Number of EOS orbitals must be smaller than (or equal to) total number of &
           & virtual orbitals!',DECinfo%output)
     end if

     ! 6. EOS indices must not exceed total number of virtual orbitals
     do i=1,nEOS
        if(EOS_idx(i) > nvirt) then
           write(DECinfo%output,'(a,i6,a)') 'EOS index number ', i, ' is larger than nvirt!'
           write(DECinfo%output,*) 'nvirt   = ', nvirt
           write(DECinfo%output,*) 'EOS_idx = ', EOS_idx(i)
           call lsquit('tensor_extract_decnp_indices: &
              & EOS index value larger than nvirt!',DECinfo%output)
        end if
     end do

     ! Extract virtual EOS indices and store in ArrVir
     ! ***********************************************

     ! Initiate ArrVir with new dimensions (nvirt_EOS,nocc,nvirt_EOS,nocc) on the
     ! local node, as the tensor will be small enough to store locally
     call tensor_init(ArrVir,new_dims,4)

     select case ( tensor_full%itype )
     case( TT_DENSE, TT_REPLICATED )

        ! Set ArrVir equal to the EOS indices of the original ArrVir array (tensor_full)
        do j=1,nocc
           do b=1,nvirt
              do i=1,nocc
                 do a=1,nEOS
                    ax=EOS_idx(a)
                    ArrVir%elm4(a,i,b,j) = tensor_full%elm4(ax,i,b,j)
                 end do
              end do
           end do
        end do

     case( TT_TILED_DIST )

        call tensor_zero(ArrVir)

        call lspdm_extract_decnp_indices_virt(ArrVir,tensor_full,nEOS,EOS_idx)

     case default
        call lsquit("ERROR(tensor_extract_decnp_indices): NO PDM version implemented yet",-1)
     end select

     EOS_idx  => null()

  end subroutine tensor_extract_decnp_indices


  subroutine get_starting_guess(iajb,t2, oof, vvf, local, spec, prec)
     implicit none
     type(tensor), intent(inout) :: iajb, t2, oof, vvf
     logical, intent(in) :: local, prec
     character(*), intent(in) :: spec
     real(tensor_dp), pointer :: o2v2(:)
     real(tensor_dp), pointer :: wrk(:)
     integer(kind=long) :: iwrk
     integer :: no,nv,i,j,a,b, specint
     real(tensor_dp), pointer :: elm4(:,:,:,:)

     no = t2%dims(4)
     nv = t2%dims(1)

     if( local )then


        select case(spec)
        case("MP2AMP")

           call array_reorder_4d(1.0E0_tensor_dp,iajb%elm1,iajb%dims(1),iajb%dims(2),&
              &iajb%dims(3),iajb%dims(4),[2,4,1,3],0.0E0_tensor_dp,t2%elm1)

        case ("CCSD_LAG_RHS")

           call array_reorder_4d(-4.0E0_tensor_dp,iajb%elm1,iajb%dims(1),iajb%dims(2),&
              &iajb%dims(3),iajb%dims(4),[2,4,1,3],0.0E0_tensor_dp,t2%elm1)
           call array_reorder_4d(2.0E0_tensor_dp,iajb%elm1,iajb%dims(1),iajb%dims(2),&
              &iajb%dims(3),iajb%dims(4),[2,4,3,1],1.0E0_tensor_dp,t2%elm1)

        case default
           call lsquit("ERROR(get_starting_guess): unknown spec",-1)
        end select

        if( prec )then
           !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(NONE) PRIVATE(i,a,j,b) SHARED(no,nv,t2,oof,vvf)
           do j = 1, no
              do i = 1, no
                 do b = 1, nv
                    do a = 1, nv
                       t2%elm4(a,b,i,j) = t2%elm4(a,b,i,j) / &
                          &(oof%elm2(i,i) - vvf%elm2(a,a) + oof%elm2(j,j) - vvf%elm2(b,b) )
                    enddo
                 enddo
              enddo
           enddo
           !$OMP END PARALLEL DO
        endif


     else

        select case(spec)
        case("MP2AMP")
           specint = 1
        case ("CCSD_LAG_RHS")
           specint = 2
        case default
           call lsquit("ERROR(get_starting_guess): unknown spec par",-1)
        end select
        call lspdm_get_starting_guess(iajb,t2,oof,vvf,specint,prec)

     endif
  end subroutine get_starting_guess

  !> \brief Reorder indices with additional memory allocation
  !> \author Patrick Ettenhuber
  subroutine tensor_reorder(arr,order)

     implicit none
     type(tensor), intent(inout) :: arr
     integer, dimension(arr%mode), intent(in) :: order
     integer, dimension(arr%mode) :: new_dims,order1,order2
     real(tensor_dp), pointer :: new_data(:)
     integer :: a,b,c,d
     integer :: dim1,dim2,dim3,dim4
     integer :: i,j
     integer :: aa,bb,cc,dd
     integer :: order_type,m,n
     real(tensor_dp) :: tcpu1,twall1,tcpu2,twall2
     integer(kind=long) :: nelms
     logical :: bg


     call LSTIMER('START',tcpu1,twall1,DECinfo%output)


     nelms = arr%nelms
     do i=1,arr%mode
        new_dims(i) = arr%dims(order(i))
     end do
     bg=(mem_is_background_buf_init().and.nelms<=mem_get_bg_buf_free())

     if( arr%itype == TT_DENSE )then

        ! Allocate space for reordered data

        call deassoc_ptr_arr(arr)

        if(bg)then
           call mem_pseudo_alloc( new_data,nelms )
        else
           call tensor_alloc_mem( new_data,nelms )
        endif

        select case(arr%mode)
        case(2)
           call array_reorder_2d(1.0E0_tensor_dp,arr%elm1,arr%dims(1),arr%dims(2),&
              & order,0.0E0_tensor_dp,new_data)
        case(3)
           call array_reorder_3d(1.0E0_tensor_dp,arr%elm1,arr%dims(1),arr%dims(2),&
              &arr%dims(3),order,0.0E0_tensor_dp,new_data)
        case(4)
           call array_reorder_4d(1.0E0_tensor_dp,arr%elm1,arr%dims(1),arr%dims(2),&
              &arr%dims(3),arr%dims(4),order,0.0E0_tensor_dp,new_data)
        case default
           call lsquit("ERROR(tensor_reorder) no default for arbitrary modes",-1)
        end select

        arr%dims=new_dims

#ifdef VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN
        call assign_in_subblocks(arr%elm1,'=',new_data,nelms)
#else
        !$OMP WORKSHARE
        arr%elm1 = new_data
        !$OMP END WORKSHARE
#endif

        if(bg)then
           call mem_pseudo_dealloc(new_data)
        else
           call tensor_free_mem(new_data)
        endif

        call assoc_ptr_arr(arr)

     else
        call lsquit("ERROR(tensor_reorder) only implemented for dense arrs yet",-1)
     endif

     call LSTIMER('START',tcpu2,twall2,DECinfo%output)

  end subroutine tensor_reorder


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!   ARRAY de-/init ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> \author Patrick Ettenhuber
  !> \date January 2013
  subroutine tensor_minit(arr, dims, nmodes, local, atype, tdims, fo, bg)
    !> the output array
    type(tensor),intent(inout) :: arr
    !> nmodes=order of the array, dims=dimensions in each mode
    integer, intent(in)              :: nmodes, dims(nmodes)
    integer, intent(in),optional     :: tdims(nmodes)
    logical, intent(in),optional     :: local, bg
    character(4),intent(in),optional :: atype
    integer,intent(in),optional :: fo
    character(4)  :: at
    integer(kind=tensor_standard_int) :: it
    logical :: loc, bg_int
    real(tensor_dp) :: time_minit
    call time_start_phase(PHASE_WORK, twall = time_minit )

    bg_int = .false.
    if(present(bg))bg_int = bg

    ! Sanity check
    if(arr%initialized)call lsquit("ERROR(tensor_minit):array already initialized",-1) 

    do i=1, nmodes
      if (dims(i) == 0) call lsquit("ERROR(tensor_minit): 0 dimendion not allowed",-1)
    end do

    !set defaults
    loc = .true.
    at  = 'LDAR'
    if(present(atype)) at  = atype


    !Default is to check at, but forcable with local
    if(present(local))then
       loc = local
    else
       select case(at)
       case('LDAR')
          loc = .true.
       case('REAR','REPD','TDAR','TDPD')
          loc = .false.
       end select
    endif

#ifdef VAR_MPI
    if(loc) then
      select case(at)
      case('LDAR','REAR','REPD','TDAR','TDPD','RTAR')
        call tensor_init_standard(arr,dims,nmodes,AT_NO_PDM_ACCESS,bg_int)
        arr%atype='LDAR'
      !case('TDAR','TDPD')
      !  arr=tensor_init_tiled(dims,nmodes,pdm=AT_NO_PDM_ACCESS)
      !  arr%atype='LTAR'
      case default
        call lsquit("ERROR(tensor_minit): atype not known",-1)
      end select
    else
      select case(at)
      case('LDAR')
        !INITIALIZE a Local Dense ARray
        call tensor_init_standard(arr,dims,nmodes,AT_MASTER_ACCESS,bg_int)
        arr%atype        = 'LDAR'
      case('TDAR')
        !INITIALIZE a Tiled Distributed ARray
        it               = TT_TILED_DIST
        call tensor_init_tiled(arr, dims,nmodes,at,it,AT_MASTER_ACCESS,bg_int,tdims=tdims,force_offset = fo)
        CreatedPDMArrays = CreatedPDMArrays+1
      case('RTAR')
        !INITIALIZE a Replicated Tiled ARray (all nodes have all tiles)
        it               = TT_TILED_REPL
        call tensor_init_tiled(arr,dims,nmodes,at,it,AT_MASTER_ACCESS,bg_int,tdims=tdims,force_offset = fo)
        CreatedPDMArrays = CreatedPDMArrays+1
      case('REAR')
        !INITIALIZE a REplicated ARray
        call tensor_init_replicated(arr,dims,nmodes,AT_MASTER_ACCESS,bg_int)
        CreatedPDMArrays = CreatedPDMArrays+1
        arr%itype        = TT_REPLICATED
        arr%atype        = 'REAR'
      case('TDPD')
        !INITIALIZE a Tiled Distributed Pseudo Dense array
        it               = TT_TILED_DIST ! for tensor_init_tiled routine
        call tensor_init_tiled(arr,dims,nmodes,at,it,AT_MASTER_ACCESS,bg_int,tdims=tdims,ps_d=.true.,force_offset=fo)
        arr%itype        = TT_DENSE ! back to dense after init
        CreatedPDMArrays = CreatedPDMArrays+1
      case('REPD')
        !INITIALIZE a REplicated Pseudo Dense array
        call tensor_init_replicated(arr,dims,nmodes,AT_MASTER_ACCESS,bg_int)
        CreatedPDMArrays = CreatedPDMArrays+1
        arr%itype        = TT_DENSE
        arr%atype        = 'REPD'
      case default 
        call lsquit("ERROR(tensor_minit): atype not known",-1)
      end select
    endif
#else
    call tensor_init(arr,dims,nmodes,bg=bg)
    arr%atype='LDAR'
#endif
    arr%initialized=.true.

    call time_start_phase(PHASE_WORK, ttot = time_minit )
    tensor_time_init = tensor_time_init + time_minit

  end subroutine tensor_minit

  subroutine tensor_ainit88(arr, dims, nmodes, local, atype, tdims, fo, bg )
    !> the output array
    type(tensor),intent(inout) :: arr
    !> nmodes=order of the array, dims=dimensions in each mode
    integer(kind=tensor_long_int), intent(in)              :: nmodes
    integer(kind=tensor_long_int), intent(in)              :: dims(nmodes)
    integer(kind=tensor_int), intent(in),optional     :: tdims(nmodes)
    logical, intent(in),optional     :: local
    logical, intent(in),optional     :: bg
    character(4),intent(in),optional :: atype
    integer(kind=tensor_int),intent(in),optional :: fo
    call tensor_ainit_central(arr, &
       &int(dims,kind=tensor_long_int), &
       &int(nmodes,kind=tensor_long_int), local=local, atype=atype, &
       &tdims=tdims, fo=fo, bg=bg )
  end subroutine tensor_ainit88
  subroutine tensor_ainit84(arr, dims, nmodes, local, atype, tdims, fo, bg )
    !> the output array
    type(tensor),intent(inout) :: arr
    !> nmodes=order of the array, dims=dimensions in each mode
    integer(kind=tensor_standard_int), intent(in)              :: nmodes
    integer(kind=tensor_long_int), intent(in)              :: dims(nmodes)
    integer(kind=tensor_int), intent(in),optional     :: tdims(nmodes)
    logical, intent(in),optional     :: local
    logical, intent(in),optional     :: bg
    character(4),intent(in),optional :: atype
    integer(kind=tensor_int),intent(in),optional :: fo
    call tensor_ainit_central(arr, &
       &int(dims,kind=tensor_long_int), &
       &int(nmodes,kind=tensor_long_int), local=local, atype=atype, &
       &tdims=tdims, fo=fo, bg=bg )
  end subroutine tensor_ainit84
  subroutine tensor_ainit48(arr, dims, nmodes, local, atype, tdims, fo, bg )
    !> the output array
    type(tensor),intent(inout) :: arr
    !> nmodes=order of the array, dims=dimensions in each mode
    integer(kind=tensor_long_int), intent(in)              :: nmodes
    integer(kind=tensor_standard_int), intent(in)              :: dims(nmodes)
    integer(kind=tensor_int), intent(in),optional     :: tdims(nmodes)
    logical, intent(in),optional     :: local
    logical, intent(in),optional     :: bg
    character(4),intent(in),optional :: atype
    integer(kind=tensor_int),intent(in),optional :: fo
    call tensor_ainit_central(arr, &
       &int(dims,kind=tensor_long_int), &
       &int(nmodes,kind=tensor_long_int), local=local, atype=atype, &
       &tdims=tdims, fo=fo, bg=bg )
  end subroutine tensor_ainit48
  subroutine tensor_ainit44(arr, dims, nmodes, local, atype, tdims, fo, bg )
    !> the output array
    type(tensor),intent(inout) :: arr
    !> nmodes=order of the array, dims=dimensions in each mode
    integer(kind=tensor_standard_int), intent(in)                  :: nmodes
    integer(kind=tensor_standard_int), intent(in)              :: dims(nmodes)
    integer(kind=tensor_int), intent(in),optional     :: tdims(nmodes)
    logical, intent(in),optional     :: local
    logical, intent(in),optional     :: bg
    character(4),intent(in),optional :: atype
    integer(kind=tensor_int),intent(in),optional :: fo
    call tensor_ainit_central(arr, &
       &int(dims,kind=tensor_long_int), &
       &int(nmodes,kind=tensor_long_int), local=local, atype=atype, &
       &tdims=tdims, fo=fo, bg=bg )
  end subroutine tensor_ainit44

  subroutine tensor_ainit_central(arr, dims_in, nmodes_in, local, atype, tdims, fo, bg )
    !> the output array
    type(tensor),intent(inout) :: arr
    !> nmodes=order of the array, dims=dimensions in each mode
    integer(kind=tensor_long_int), intent(in)              :: nmodes_in
    integer(kind=tensor_long_int), intent(in)              :: dims_in(nmodes_in)
    integer, intent(in),optional     :: tdims(nmodes_in)
    logical, intent(in),optional     :: local, bg
    character(4),intent(in),optional :: atype
    integer(kind=tensor_int),intent(in),optional :: fo
    character(4)  :: at
    integer(kind=tensor_standard_int):: it
    logical :: loc, bg_int
    real(tensor_dp) :: time_ainit
    integer :: dims(nmodes_in),nmodes
    call time_start_phase(PHASE_WORK, twall = time_ainit )
    dims = dims_in
    nmodes = nmodes_in
 
    bg_int = .false.
    if(present(bg))bg_int = bg

    ! Sanity check
    if(arr%initialized)call lsquit("ERROR(tensor_ainit_central):tensor already initialized",-1) 
    do i=1, nmodes
      if (dims(i) == 0) call lsquit("ERROR(tensor_ainit_central): 0 dimendion not allowed",-1)
    end do
 
    !set defaults
    loc = .true.
    at  = 'LDAR'
    if(present(atype)) at  = atype

    !Default is to check at, but forcable with local
    if(present(local))then
       loc = local
    else
       select case(at)
       case('LDAR')
          loc = .true.
       case('REAR','REPD','TDAR','TDPD')
          loc = .false.
       end select
    endif

#ifdef VAR_MPI
    if(loc) then
      select case(at)
      case('LDAR','REAR','REPD','TDAR','TDPD')
        !if local recast to a local dense array
        call tensor_init_standard(arr,int(dims),int(nmodes),AT_NO_PDM_ACCESS,bg_int)
        arr%atype='LDAR'
      !case('TDAR','TDPD')
      !  arr=tensor_init_tiled(dims,nmodes,pdm=AT_NO_PDM_ACCESS)
      !  arr%atype='LTAR'
      case default
        call lsquit("ERROR(tensor_minit): atype not known",-1)
      end select
    else
      select case(at)
      case('LDAR')
        !INITIALIZE a Local Dense ARray
        call tensor_init_standard(arr,dims,nmodes,AT_ALL_ACCESS,bg_int)
        arr%atype        = 'LDAR'
      case('TDAR')
        !INITIALIZE a Tiled Distributed ARray
        it               = TT_TILED_DIST
        call tensor_init_tiled(arr,dims,nmodes,at,it,AT_ALL_ACCESS,bg_int,&
           &tdims=int(tdims,kind=tensor_int),force_offset=fo)
        CreatedPDMArrays = CreatedPDMArrays+1
      case('REAR')
        !INITIALIZE a REplicated ARray
        call tensor_init_replicated(arr,dims,nmodes,AT_ALL_ACCESS,bg_int)
        CreatedPDMArrays = CreatedPDMArrays+1
        arr%itype        = TT_REPLICATED
        arr%atype        = 'REAR'
      case('TDPD')
        !INITIALIZE a Tiled Distributed Pseudo Dense array
        it               = TT_TILED_DIST ! for tensor_init_tiled routine
        call tensor_init_tiled(arr,dims,nmodes,at,it,AT_ALL_ACCESS,bg_int,&
           &tdims=int(tdims,kind=tensor_int),ps_d=.true.,force_offset=fo)
        arr%itype        = TT_DENSE ! back to dense after init
        CreatedPDMArrays = CreatedPDMArrays+1
      case('REPD')
        !INITIALIZE a REplicated Pseudo Dense array
        call tensor_init_replicated(arr,dims,nmodes,AT_ALL_ACCESS,bg_int)
        CreatedPDMArrays = CreatedPDMArrays+1
        arr%itype        = TT_DENSE
        arr%atype        = 'REPD'
      case default 
        call lsquit("ERROR(tensor_ainit_central): atype not known",-1)
      end select
    endif
#else
    call tensor_init_standard(arr,dims,nmodes,AT_NO_PDM_ACCESS,bg_int)
    arr%atype='LDAR'
#endif
    arr%initialized=.true.

    call time_start_phase(PHASE_WORK, ttot = time_ainit )
    tensor_time_init = tensor_time_init + time_ainit
  end subroutine tensor_ainit_central

  !> \author Patrick Ettenhuber
  !> \date September 2012
  !> \brief MAIN ARRAY INITIALIZATION ROUTINE
  subroutine  tensor_init(arr,dims,nmodes,tensor_type,pdm,tdims,fo,bg)
    implicit none
    !> output array
    type(tensor),intent(inout) :: arr
    !> nmodes=order of the array, dims=dimensions in each mode
    integer, intent(in) :: nmodes, dims(nmodes)
    !> integer specifying the type of array (a list of possible types is found
    !> at the beginning of tensor_memory.f90
    integer(kind=tensor_standard_int), optional :: tensor_type
    !> if tiled then the size of the tile in each mode can be specified explicitly 
    integer, optional :: tdims(nmodes)
    !> specifies the type of access to the array (AT_NO_PDM_ACCESS,AT_MASTER_ACCESS,AT_ALL_ACCESS)
    integer(kind=tensor_standard_int), optional :: pdm
    integer, optional :: fo
    logical, optional :: bg 
    integer :: sel_type
    integer(kind=tensor_standard_int) :: pdmtype,it
    logical :: zeros_in_tiles,wcps, bg_int
    real(tensor_dp) :: time_init

    bg_int = .false.
    if(present(bg))bg_int = bg

    !choose which kind of array
    call time_start_phase(PHASE_WORK, twall = time_init )

    !if(arr%initialized)call lsquit("ERROR(tensor_init):array already initialized",-1) 

    !DEFAULTS
    it      = TT_DENSE
    pdmtype = AT_NO_PDM_ACCESS !NO PDM

    !OPTIONAL SPECIFICATIONS
    if(present(tensor_type)) it      = tensor_type
    if(present(pdm))         pdmtype = pdm

    !CHECK INPUT
    if(pdmtype>=3)call lsquit("ERROR(tensor_init):WRONG CHOICE IN PDMTYPE",DECinfo%output)

    !EXPERIMENTAL, THIS IS NOT RECOMMENDED!!!!!!!!!!!!!!:
    !instead of modulo dimensions in the rims of an array
    !use the same dimensions as in full tiles and fill with
    !zeros
    zeros_in_tiles = .false.
    ArraysCreated = ArraysCreated+1
    
    !select corresponding routine
    select case(it)
      case(TT_DENSE)
        call tensor_init_standard(arr,dims,nmodes,pdmtype,bg_int)
        arr%atype = 'LDAR'
      case(TT_REPLICATED)
        call tensor_init_replicated(arr,dims,nmodes,pdmtype,bg_int)
        arr%atype = 'REAR'
        CreatedPDMArrays = CreatedPDMArrays+1
      case(TT_TILED)
        call tensor_init_tiled(arr,dims,nmodes,'TIAR',it,pdmtype,bg_int,tdims=tdims,force_offset=fo)
      case(TT_TILED_DIST)
        call tensor_init_tiled(arr,dims,nmodes,'TDAR',it,pdmtype,bg_int,tdims=tdims,force_offset=fo)
        CreatedPDMArrays = CreatedPDMArrays+1
    end select
    arr%access_type   = pdmtype
    arr%itype         = it
    arr%initialized   = .true.

    call time_start_phase(PHASE_WORK, ttot = time_init )
    tensor_time_init = tensor_time_init + time_init
  end subroutine tensor_init


  !> \author Patrick Ettenhuber
  !> \date September 2012
  !> \brief get mode index from composite index
  subroutine tensor_init_standard(arr,dims,nmodes,pdm,bg)
    implicit none
    integer, intent(in)   :: nmodes,dims(nmodes)
    integer(kind=tensor_standard_int), intent(in)   :: pdm
    type(tensor),intent(inout) :: arr
    logical, intent(in)   :: bg
    logical               :: master
    integer               :: i,addr,tdimdummy(nmodes)
    integer,pointer       :: buf(:)
    integer(kind=long)    :: nelms
    integer(kind=tensor_mpi_kind) :: pc_nnodes,me

    master    = .true.
    me        = 0
    
    
    !find space in the persistent array
    p_arr%curr_addr_on_node   = get_free_address(.true.)
    addr                      = p_arr%curr_addr_on_node
    p_arr%arrays_in_use       = p_arr%arrays_in_use + 1
    p_arr%a(addr)%local_addr  = addr
    p_arr%a(addr)%initialized = .true.
    !set to invalid, since not used here
    p_arr%a(addr)%nnod        = -1
    p_arr%a(addr)%comm        = -1

    !SET MODE
    p_arr%a(addr)%mode      = nmodes

    !SET DIMS
    call tensor_set_dims(p_arr%a(addr),dims,nmodes)

    !SET ARRAY TYPE
    p_arr%a(addr)%itype     = TT_DENSE

    !SET INIT TYPE
    !default
    allocate( p_arr%a(addr)%access_type )
    p_arr%a(addr)%access_type = AT_NO_PDM_ACCESS
    !if one uses comm threads the following replace the access_type
    !if( pdm == AT_MASTER_ACCESS .and. lspdm_use_comm_proc )&
    !& p_arr%a(addr)%access_type = AT_MASTER_ACCESS
    !if( pdm == AT_ALL_ACCESS .and. lspdm_use_comm_proc )&
    !& p_arr%a(addr)%access_type = AT_ALL_ACCESS

    !SET IF ALLOCATED WITH COMM PROCS
    !p_arr%a(addr)%allocd_w_c_p = lspdm_use_comm_proc

    !SET NELMS
    nelms=1
    do i=1,nmodes
      nelms=nelms*dims(i)
    enddo
    p_arr%a(addr)%nelms=nelms

    !put 0 in tdim, since for the replicated array it is not important
    tdimdummy=0
    call tensor_set_tdims(p_arr%a(addr),tdimdummy,nmodes)

    !ALLOCATE STORAGE SPACE FOR THE ARRAY
    call memory_allocate_tensor_dense(p_arr%a(addr),bg)

    !RETURN THE CURRENLY ALLOCATE ARRAY
    arr=p_arr%a(addr)

  end subroutine tensor_init_standard

  

  subroutine tensor_free_standard(arr)
    implicit none
    type(tensor), intent(inout) :: arr
    integer :: me
    logical :: parent
 
    me     = 0
    parent = .true.
    p_arr%free_addr_on_node(arr%local_addr)=.true.
    p_arr%arrays_in_use = p_arr%arrays_in_use - 1 
    call tensor_free_basic(p_arr%a(arr%local_addr)) 
    call tensor_reset_value_defaults(p_arr%a(arr%local_addr)) 
    call tensor_nullify_pointers(arr)

  end subroutine tensor_free_standard

  !> \author Patrick Ettenhuber
  !> \date January 2013
  !> \brief free a replicated matrix wich is pseudo dense for the routines,
  !> wrapper for use without mpi
  !> without mpi
  subroutine tensor_free_rpseudo_dense(arr,local)
    !> the array to free
    type(tensor) :: arr
    logical, intent(in) :: local
#ifdef VAR_MPI
    if(.not.local)then
      arr%itype=TT_REPLICATED
    endif
#endif
    call tensor_free(arr)
  end subroutine tensor_free_rpseudo_dense

  !> \author Patrick Ettenhuber
  !> \date January 2013
  !> \brief free a td matrix wich is pseudo dense for the routines,
  !> wrapper for use without mpi
  !> without mpi
  subroutine tensor_free_tdpseudo_dense(arr,local)
    !> the array to free
    type(tensor) :: arr
    logical, intent(in) :: local
#ifdef VAR_MPI
    if( .not. local ) then
      arr%itype=TT_TILED_DIST
    endif
#endif
    call tensor_free(arr)
  end subroutine tensor_free_tdpseudo_dense

  !> \brief array freeing routine, give an arbitrary array and all allocated
  !memory associated with the array will be freed
  !> \author Patrick Ettenhuber
  !> \Date probably september 2012
  subroutine tensor_free(arr)
    implicit none
    !> array to free
    type(tensor),intent(inout) :: arr

    if(.not.arr%initialized)call lsquit("ERROR(tensor_free):array not initialized",-1) 

    select case(arr%atype)
    case('LDAR')
      call tensor_free_standard(arr)
    case('TDAR','REAR','TDPD','REPD','RTAR')
      call tensor_free_pdm(arr)
      DestroyedPDMArrays = DestroyedPDMArrays + 1
    case('TIAR')
        call lsquit("ERROR(tensor_free): local tiled not maintained",-1)
    case default 
        call lsquit("ERROR(tensor_free): atype not known",-1)
    end select
    arr%initialized = .false. 
  end subroutine tensor_free


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!   ARRAY CONVERSION ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> \brief deallocate the dense part of an array and change the type to tiled
  !distributed array. in the case of non mpi builds, this routine does nothing,
  !so that implementations still run
  !> \author Patrick Ettenhuber
  !> \date January 2013
  subroutine tensor_change_itype_to_td(arr,local)
    implicit none
    !> array to change the array type
    type(tensor),intent(inout) :: arr
    logical, intent(in) :: local
#ifdef VAR_MPI
    if( .not. local )then
      arr%itype=TT_TILED_DIST
      if(associated(arr%elm1))then
        call tensor_deallocate_dense(arr)
      endif
    endif
#else
    return
#endif
  end subroutine tensor_change_itype_to_td

  !> \brief change the array type to replicated, no action in case of
  !non-mpi-build
  !> \author Patrick Ettenhuber
  !> \date January 2012
  subroutine tensor_change_atype_to_rep(arr,local)
    implicit none
    !> array to change the array type
    type(tensor),intent(inout) :: arr
    logical, intent(in) :: local
#ifdef VAR_MPI
    if(.not.local)then
      arr%itype=TT_REPLICATED
    endif
#else
    return
#endif
  end subroutine tensor_change_atype_to_rep

  !> \brief change the array type to replicated, no action in case of
  !non-mpi-build
  !> \author Patrick Ettenhuber
  !> \date January 2012
  subroutine tensor_change_atype_to_d(arr)
    implicit none
    !> array to change the array type
    type(tensor),intent(inout) :: arr
#ifdef VAR_MPI
    arr%itype=TT_DENSE
#else
    return
#endif
  end subroutine tensor_change_atype_to_d

  !> \brief copy the tiled data of an array to its dense part, if change is
  !true, also change the %itype to TT_DENSE, but NO deallocation of the tiled
  !distributed part
  !> \author Patrick Ettenhuber
  !> \date January 2012
  subroutine tensor_cp_tiled2dense(arr,change,order,bg)
    implicit none
    !> array to copy data from the tiled to its dense part
    type(tensor),intent(inout) :: arr
    !> logical to specify whether to change the atype
    logical :: change
    !> if order is given the dense part will be reordered with respect to the
    !tiled distributed part
    integer,intent(in),optional:: order(arr%mode)
    logical, intent(in),optional :: bg
    logical :: pdm, bg_int
    pdm=.false.
    bg_int = .false.
    if(present(bg)) bg_int = bg

    if(arr%itype/=TT_DENSE.and.arr%itype/=TT_REPLICATED)then
      if(.not.associated(arr%elm1))then
        call memory_allocate_tensor_dense(arr,bg_int)
      else
        call lsquit("ERROR(tensor_cp_tiled2dense):dense is already allocated,&
        & please make sure you are not doing someting stupid",DECinfo%output)
      endif
      if(arr%access_type>0)pdm=.true.
      if(.not.present(order))call cp_tileddata2fort(arr,arr%elm1,arr%nelms,pdm)
      if(present(order))call cp_tileddata2fort(arr,arr%elm1,arr%nelms,pdm,order)
      if(change)arr%itype=TT_DENSE
    endif
  end subroutine tensor_cp_tiled2dense

  !> \brief copy the dense part of an array to its tiled distributed part and
  !delallocate the dense part afterwards, if desired change the atype. no
  !action if non-mpi build.
  !> \author Patrick Ettenhuber
  !> \date January 2012
  subroutine tensor_mv_dense2tiled(arr,change,dealloc_local)
    implicit none
    !> array to copy dense to tiled and deallocate dense part
    type(tensor),intent(inout) :: arr
    logical, intent(in) :: change
    logical, intent(in),optional :: dealloc_local
    logical :: pdm,dl
#ifdef VAR_MPI
    pdm=.false.
    if(.not.associated(arr%elm1))then
      call lsquit("ERROR(tensor_cp_dense2tiled):dense is NOT allocated,&
      & please make sure you are not doing someting stupid",DECinfo%output)
    endif

    dl = .true.
    if(present(dealloc_local))dl = dealloc_local

    if(arr%access_type>0)pdm=.true.
    if(change)arr%itype=TT_TILED_DIST
    call tensor_convert_fort2arr(arr%elm1,arr,arr%nelms)
    if(dl)call tensor_deallocate_dense(arr)
#else
    return
#endif
  end subroutine tensor_mv_dense2tiled


  !\brief all the following wrappers are necessary to use the conversion routine
  !in an interface for different shapes of the fortran array  
  subroutine tensor_convert_tensor2fort_wrapper1(arr,fort,order,wrk,iwrk)
    implicit none
    type(tensor), intent(inout) :: arr
    real(tensor_dp), intent(inout) :: fort(:)
    integer, intent(in), optional :: order(arr%mode)
    real(tensor_dp),intent(inout),target,optional :: wrk(*)
    integer(kind=8),intent(in),optional,target:: iwrk
    call tensor_convert_tensor2fort(arr,fort,arr%nelms,order=order,wrk=wrk,iwrk=iwrk)
  end subroutine tensor_convert_tensor2fort_wrapper1
  subroutine tensor_convert_tensor2fort_wrapper2(arr,fort,order,wrk,iwrk)
    implicit none
    type(tensor), intent(inout) :: arr
    real(tensor_dp), intent(inout) :: fort(:,:)
    integer, intent(in), optional :: order(arr%mode)
    real(tensor_dp),intent(inout),target,optional :: wrk(*)
    integer(kind=8),intent(in),optional,target:: iwrk
    call tensor_convert_tensor2fort(arr,fort,arr%nelms,order=order,wrk=wrk,iwrk=iwrk)
  end subroutine tensor_convert_tensor2fort_wrapper2
  subroutine tensor_convert_tensor2fort_wrapper3(arr,fort,order,wrk,iwrk)
    implicit none
    type(tensor), intent(inout) :: arr
    real(tensor_dp), intent(inout) :: fort(:,:,:)
    integer, intent(in), optional :: order(arr%mode)
    real(tensor_dp),intent(inout),target,optional :: wrk(*)
    integer(kind=8),intent(in),optional,target:: iwrk
    call tensor_convert_tensor2fort(arr,fort,arr%nelms,order=order,wrk=wrk,iwrk=iwrk)
  end subroutine tensor_convert_tensor2fort_wrapper3
  subroutine tensor_convert_tensor2fort_wrapper4(arr,fort,order,wrk,iwrk)
    implicit none
    type(tensor), intent(inout) :: arr
    real(tensor_dp), intent(inout) :: fort(:,:,:,:)
    integer, intent(in), optional :: order(arr%mode)
    real(tensor_dp),intent(inout),target,optional :: wrk(*)
    integer(kind=8),intent(in),optional,target:: iwrk
    call tensor_convert_tensor2fort(arr,fort,arr%nelms,order=order,wrk=wrk,iwrk=iwrk)
  end subroutine tensor_convert_tensor2fort_wrapper4

  !\brief all the following wrappers are necessary to use the conversion routine
  !in an interface for different shapes of the fortran array  
  subroutine tensor_convert_fort2tensor_wrapper1(fortarr,arr,order,wrk,iwrk)
    implicit none
    type(tensor), intent(inout) :: arr
    real(tensor_dp), intent(in) :: fortarr(arr%nelms)
    integer, intent(in),optional :: order(arr%mode)
    real(tensor_dp),intent(inout),target,optional :: wrk(*)
    integer(kind=8),intent(in),optional,target:: iwrk
    call tensor_convert_fort2arr(fortarr,arr,arr%nelms,order=order,wrk=wrk,iwrk=iwrk)
  end subroutine tensor_convert_fort2tensor_wrapper1
  subroutine tensor_convert_fort2tensor_wrapper2(fortarr,arr,order,wrk,iwrk)
    implicit none
    real(tensor_dp), intent(in) :: fortarr(:,:)
    type(tensor), intent(inout) :: arr
    real(tensor_dp),intent(inout),target,optional :: wrk(*)
    integer(kind=8),intent(in),optional,target:: iwrk
    integer, intent(in),optional :: order(arr%mode)
    call tensor_convert_fort2arr(fortarr,arr,arr%nelms,order=order,wrk=wrk,iwrk=iwrk)
  end subroutine tensor_convert_fort2tensor_wrapper2
  subroutine tensor_convert_fort2tensor_wrapper3(fortarr,arr,order,wrk,iwrk)
    implicit none
    real(tensor_dp), intent(in) :: fortarr(:,:,:)
    type(tensor), intent(inout) :: arr
    integer, intent(in),optional :: order(arr%mode)
    real(tensor_dp),intent(inout),target,optional :: wrk(*)
    integer(kind=8),intent(in),optional,target:: iwrk
    call tensor_convert_fort2arr(fortarr,arr,arr%nelms,order=order,wrk=wrk,iwrk=iwrk)
  end subroutine tensor_convert_fort2tensor_wrapper3
  subroutine tensor_convert_fort2tensor_wrapper4(fortarr,arr,order,wrk,iwrk)
    implicit none
    real(tensor_dp), intent(in) :: fortarr(:,:,:,:)
    type(tensor), intent(inout) :: arr
    integer, intent(in),optional :: order(arr%mode)
    real(tensor_dp),intent(inout),target,optional :: wrk(*)
    integer(kind=8),intent(in),optional,target:: iwrk
    call tensor_convert_fort2arr(fortarr,arr,arr%nelms,order=order,wrk=wrk,iwrk=iwrk)
  end subroutine tensor_convert_fort2tensor_wrapper4


  !> \brief put data of a fortan array into an arbitrary array
  !> \author Patrick Ettenhuber
  !> \date late 2012
  subroutine tensor_convert_fort2arr(fortarr,arr,nelms,order,wrk,iwrk)
    implicit none
    !> the fortran array with the data
    real(tensor_dp), intent(in) :: fortarr(*)
    !> the array which should contain the data after the operation
    type(tensor), intent(inout) :: arr
    !> number of elements to copy from the fortan array to the array
    integer(kind=8), intent(in) :: nelms
    !> if the array should have a different ordering than the fortran array,
    ! this can be specified with order
    integer, intent(in),optional :: order(arr%mode)
    real(tensor_dp),intent(inout),target,optional :: wrk(*)
    integer(kind=8),intent(in),optional,target:: iwrk
    real(tensor_dp) :: tilemem,MemFree
    integer :: i,o(arr%mode),fullfortdims(arr%mode)
    real(tensor_dp) :: nrm
    logical :: simpleord

    simpleord = .true.

    do  i=1,arr%mode
      o(i)=i
    enddo
    if(present(order))o=order
    do  i=1,arr%mode
      fullfortdims(o(i))=arr%dims(i)
    enddo

    do  i=1,arr%mode
      simpleord = simpleord.and.(o(i)==i)
    enddo
    
    if(nelms/=arr%nelms)call lsquit("ERROR(tensor_convert_fort2arr):array&
       &dimensions are not the same",-1)

    select case(arr%itype)
    case(TT_DENSE,TT_REPLICATED)
       if(simpleord)then
#ifdef VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN
          call assign_in_subblocks(arr%elm1,'=',fortarr,nelms)
#else
          call dcopy(int(nelms),fortarr,1,arr%elm1,1)
#endif
       else
          select case(arr%mode)
          case(2)
             call array_reorder_2d(1.0E0_tensor_dp,fortarr,fullfortdims(1),fullfortdims(2),&
                &o,0.0E0_tensor_dp,arr%elm1)
          case(3)
             call array_reorder_3d(1.0E0_tensor_dp,fortarr,fullfortdims(1),fullfortdims(2),&
                &fullfortdims(3),o,0.0E0_tensor_dp,arr%elm1)
          case(4)
             call array_reorder_4d(1.0E0_tensor_dp,fortarr,fullfortdims(1),fullfortdims(2),&
                &fullfortdims(3),fullfortdims(4),o,0.0E0_tensor_dp,arr%elm1)
          case default
             call lsquit("ERROR(tensor_convert_fort2arr): mode not implemented",-1)
          end select
       endif

       if(arr%itype==TT_REPLICATED)call tensor_sync_replicated(arr)

    case(TT_TILED)

       call cp_data2tiled_lowmem(arr,fortarr,arr%dims,int(arr%mode))

    case(TT_TILED_DIST)

       if(arr%access_type==AT_ALL_ACCESS)then
          do i=1,arr%nlti
             call tile_from_fort(1.0E0_tensor_dp,fortarr,fullfortdims,int(arr%mode),&
                &0.0E0_tensor_dp,arr%ti(i)%t,int(arr%ti(i)%gt),int(arr%tdim),o)
          enddo
       else
          call tensor_scatter(1.0E0_tensor_dp,fortarr,0.0E0_tensor_dp,arr,nelms,oo=o,wrk=wrk,iwrk=iwrk)
       endif

    case default
       call lsquit("ERROR(tensor_convert_fort2arr) the array type is not implemented",-1)
    end select
 end subroutine tensor_convert_fort2arr

 !> \brief change the init type for a fortan array
 !> \author Patrick Ettenhuber
 !> \date late 2012
 subroutine change_access_type(arr,totype)
    implicit none
    !> array to chage the init type
    type(tensor),intent(inout) :: arr
    !> type to change it to
    integer,intent(in) :: totype
    if(arr%itype==TT_TILED_DIST.or.arr%itype==TT_REPLICATED.or. &
       & totype==TT_TILED_DIST.or.totype==TT_REPLICATED)then
       call change_access_type_td(arr,totype)
    else
       call lsquit("ERROR(change_access_type): what you want to do is not implemented",-1)
    endif
  end subroutine change_access_type

  
  !> \brief put data of an arbitrary array into a basic fortan type array
  !> \author Patrick Ettenhuber
  !> \date late 2012
  subroutine tensor_convert_tensor2fort(arr,fort,nelms,order,wrk,iwrk)
    implicit none
    !> array with the data at the beginning
    type(tensor), intent(inout) :: arr
    !> fortan array to contain the data in the end
    real(tensor_dp), intent(inout) :: fort(*)
    !> number of elements to convert, must be the same as elements in the array
    integer(kind=8), intent(in) :: nelms
    !> if the fortan array has a different order than the array this specifies
    !the reordering
    integer, intent(in), optional :: order(arr%mode)
    real(tensor_dp),intent(inout),target,optional :: wrk(*)
    integer(kind=8),intent(in),optional,target:: iwrk
    real(tensor_dp) :: tilemem,MemFree
    integer :: i,o(arr%mode),fullfortdims(arr%mode)
    real(tensor_dp) :: nrm
    logical :: simpleord

    simpleord = .true.

    do  i=1,arr%mode
      o(i)=i
    enddo
    if(present(order))o=order
    do  i=1,arr%mode
      fullfortdims(o(i))=arr%dims(i)
    enddo

    do  i=1,arr%mode
      simpleord = simpleord.and.(o(i)==i)
    enddo


    if(nelms/=arr%nelms)call lsquit("ERROR(tensor_convert_tensor2fort):array&
    &dimensions are not the same",DECinfo%output)


    select case(arr%itype)

    case(TT_DENSE,TT_REPLICATED)

       if(simpleord)then
#ifdef VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN
          call assign_in_subblocks(fort,'=',arr%elm1,nelms)
#else
          call dcopy(int(nelms),arr%elm1,1,fort,1)
#endif
       else

          select case(arr%mode)
          case(2)
             call array_reorder_2d(1.0E0_tensor_dp,arr%elm1,arr%dims(1),arr%dims(2),&
                &o,0.0E0_tensor_dp,fort)
          case(3)
             call array_reorder_3d(1.0E0_tensor_dp,arr%elm1,arr%dims(1),arr%dims(2),&
                &arr%dims(3),o,0.0E0_tensor_dp,fort)
          case(4)
             call array_reorder_4d(1.0E0_tensor_dp,arr%elm1,arr%dims(1),arr%dims(2),&
                &arr%dims(3),arr%dims(4),o,0.0E0_tensor_dp,fort)
          case default
             call lsquit("ERROR(tensor_convert_fort2arr): mode not implemented",-1)
          end select

       endif

    case(TT_TILED)
       call cp_tileddata2fort(arr,fort,nelms,.false.,order=order)
    case(TT_TILED_DIST)
       call tensor_gather(1.0E0_tensor_dp,arr,0.0E0_tensor_dp,fort,nelms,oo=order,wrk=wrk,iwrk=iwrk)
    end select
  end subroutine tensor_convert_tensor2fort


  
  !> \brief convert an old array2 structure to an array structure
  !> \author Patrick Ettenhuber
  !> \date late 2012
  subroutine tensor_convert_array22array(arraytwo,arr)
    implicit none
    !> array2 input
    type(array2),intent(inout) :: arraytwo
    !> array output
    type(tensor), intent(inout) :: arr
    integer(kind=8) :: nel
    if(arr%mode/=2)then
      call lsquit("ERROR(tensor_convert_array22array):wrong mode in arr&
      &input)",DECinfo%output)
    endif
    nel = int(arraytwo%dims(1)*arraytwo%dims(2),kind=8)
    if(arr%nelms/=nel)then
      call lsquit("ERROR(tensor_convert_array22array):number of elements in arr&
      &input)",DECinfo%output)
    endif
    call tensor_convert_fort2arr(arraytwo%val,arr,nel)
  end subroutine tensor_convert_array22array


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!   ARRAY UTILITIES !!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tensor_cp_data(from_arr,to_arr,order,wrk,iwrk)
    implicit none
    type(tensor),intent(inout) :: from_arr
    type(tensor),intent(inout) :: to_arr
    integer, intent(in),optional :: order(to_arr%mode)
    real(tensor_dp),intent(inout),target,optional :: wrk(*)
    integer(kind=8),intent(in),optional,target:: iwrk
    integer :: i
    real(tensor_dp) :: tilemem,MemFree
    integer :: o(from_arr%mode)
    if(from_arr%nelms/=to_arr%nelms)then
      call lsquit("ERROR(tensor_cp_data):arrays need the same number of& 
      & elements",DECinfo%output)
    endif

    do i=1,from_arr%mode
       o(i) = i
    enddo
    if(present(order)) o = order
    
    select case(from_arr%itype)


    case(TT_DENSE,TT_REPLICATED)
       select case(to_arr%itype)

       case(TT_DENSE,TT_REPLICATED)

          call tensor_convert(from_arr, to_arr%elm1, order = o)
          if(to_arr%itype == TT_REPLICATED)call tensor_sync_replicated(to_arr)

       case(TT_TILED_DIST)

          if(to_arr%access_type==AT_ALL_ACCESS)then
             do i=1,to_arr%nlti
                call tile_from_fort(1.0E0_tensor_dp,from_arr%elm1,from_arr%dims,int(to_arr%mode),&
                   &0.0E0_tensor_dp,to_arr%ti(i)%t,int(to_arr%ti(i)%gt),int(to_arr%tdim),o)
             enddo
          else
             call tensor_scatter(1.0E0_tensor_dp,from_arr%elm1,0.0E0_tensor_dp,to_arr,from_arr%nelms,oo=o,wrk=wrk,iwrk=iwrk)
          endif

       case default
          call lsquit("ERROR(tensor_cp_data):operation not yet implemented",DECinfo%output)
       end select



    case(TT_TILED_DIST)

       select case(to_arr%itype)
       case(TT_DENSE)

          call cp_tileddata2fort(from_arr,to_arr%elm1,from_arr%nelms,.true.,order = o)

       case(TT_TILED_DIST)

          call tensor_cp_tiled(from_arr,to_arr,o)

       case default

          call lsquit("ERROR(tensor_cp_data):operation not yet  implemented",DECinfo%output)

       end select


    case default
        call lsquit("ERROR(tensor_cp_data):operation not yet  implemented",DECinfo%output)
    end select

  end subroutine tensor_cp_data

  subroutine tensor_zero(zeroed)
     implicit none
     type(tensor) :: zeroed
     integer :: i

     select case(zeroed%itype)
     case(TT_DENSE)
        zeroed%elm1=0.0E0_tensor_dp
     case(TT_REPLICATED)
        zeroed%elm1=0.0E0_tensor_dp
        call tensor_sync_replicated(zeroed)
     case(TT_TILED)
        if (zeroed%atype=='RTAR') then
           call tensor_zero_tiled_dist(zeroed)
        else
           do i=1,zeroed%ntiles
              zeroed%ti(i)%t=0.0E0_tensor_dp
           enddo
        end if
     case(TT_TILED_DIST,TT_TILED_REPL)
        call tensor_zero_tiled_dist(zeroed)
     case default
        call lsquit("ERROR(tensor_zero):not yet implemented",-1)
     end select

  end subroutine tensor_zero

  subroutine tensor_random(zeroed)
     implicit none
     type(tensor) :: zeroed
     integer :: i

     select case(zeroed%itype)
     case(TT_DENSE)
        call random_seed()
        call random_number(zeroed%elm1)
     case(TT_REPLICATED)
        call random_seed()
        call random_number(zeroed%elm1)
        call tensor_sync_replicated(zeroed)
     case(TT_TILED)
        if (zeroed%atype=='RTAR') then
           call tensor_rand_tiled_dist(zeroed)
        else
           call random_seed()
           do i=1,zeroed%ntiles
              call random_number(zeroed%ti(i)%t)
              
           enddo
        end if
     case(TT_TILED_DIST,TT_TILED_REPL)
        call tensor_rand_tiled_dist(zeroed)
     case default
        call lsquit("ERROR(tensor_rand):not yet implemented",-1)
     end select

  end subroutine tensor_random

  subroutine tensor_print_tile_norm(arr,globtinr,nrm,returnsquared)
    implicit none
    type(tensor),intent(in) :: arr
    integer,intent(in) :: globtinr
    real(tensor_dp),intent(inout),optional::nrm
    logical,intent(in),optional :: returnsquared
    real(tensor_dp)::norm
    integer :: loctinr,i,j,on
    logical :: squareback
    squareback=.false.
    if(present(returnsquared))squareback=returnsquared
    norm=0.0d0
    if(globtinr>arr%ntiles)then
      call lsquit("ERROR(tensor_print_tile_norm):tile does not exist",DECinfo%output)
    endif 
    select case(arr%itype)
      case(TT_DENSE)
        print *,"WARNING INVALID OPTION(tensor_print_tile_norm):no tiles in dense array"     
        do i=1,arr%nelms
          norm=norm+arr%elm1(i)*arr%elm1(i)
        enddo
        on = 1
      case(TT_TILED)
        do j=1,arr%ti(globtinr)%e
          norm=norm + arr%ti(globtinr)%t(j) * arr%ti(globtinr)%t(j)
        enddo
        on = 1
      case(TT_TILED_DIST)
        call tensor_tiled_pdm_print_ti_nrm(arr,globtinr,on,norm)
    end select
    if(.not.squareback)norm=sqrt(norm)
    if(present(nrm))nrm=norm
    if(.not.present(nrm))then
    if(.not.squareback)write(DECinfo%output,'("LOCAL TILE NORM ON",I3,f20.15)')on,norm
    if(squareback)write(DECinfo%output,'("LOCAL TILE NORM^2 ON",I3,f20.15)')on,norm
    endif
  end subroutine tensor_print_tile_norm

  subroutine tensor_print_norm_nrm(arr,nrm,returnsquared)
    implicit none
    real(tensor_dp),intent(inout),optional :: nrm
    type(tensor),intent(in) :: arr
    logical,intent(in),optional :: returnsquared
    real(tensor_dp)::norm
    integer(kind=8) :: i,j
    logical :: squareback
    squareback=.false.
    if(present(returnsquared))squareback=returnsquared
    norm=0.0d0
    select case(arr%itype)
      case(TT_DENSE)
        do i=1,arr%nelms
          norm=norm+arr%elm1(i)*arr%elm1(i)
        enddo
      case(TT_REPLICATED)
        norm = tensor_print_norm_repl(arr)
      case(TT_TILED)
        do i=1,arr%nlti
          do j=1,arr%ti(i)%e
            !ISNAN is gfort intrinsic, so this options should be commented out in
            !general
            !if(ISNAN(arr%ti(i)%t(j)))then
            !  write(*,'("NaN detected in norm_t, tile:",I5," element: ",I5)')i,j
            !  stop 1
            !endif
            norm=norm + arr%ti(i)%t(j) * arr%ti(i)%t(j)
          enddo
        enddo
      case(TT_TILED_DIST)
        norm=tensor_tiled_pdm_get_nrm2(arr)
    end select
    if(.not.squareback)norm = sqrt(norm)
    if(.not.squareback.and..not.present(nrm))print *,"NORM:",norm
    if(squareback.and..not.present(nrm))print *,"NORM^2:",norm
    if(present(nrm))nrm=norm
  end subroutine tensor_print_norm_nrm

  subroutine tensor_print_norm_customprint(arr,msg,returnsquared,print_)
    implicit none
    character*(*),intent(in) :: msg
    type(tensor),intent(in) :: arr
    logical,intent(in),optional :: returnsquared
    logical,intent(in),optional :: print_
    real(tensor_dp)::norm
    integer(kind=8) :: i,j
    logical :: squareback
    integer(kind=tensor_mpi_kind) :: me

    squareback=.false.

    if(present(returnsquared))squareback=returnsquared

    norm=0.0d0
    select case(arr%itype)
    case(TT_DENSE)
       do i=1,arr%nelms
          norm=norm+arr%elm1(i)*arr%elm1(i)
       enddo
    case(TT_REPLICATED)
       norm = tensor_print_norm_repl(arr)
    case(TT_TILED)
       do i=1,arr%nlti
          do j=1,arr%ti(i)%e
             norm=norm + arr%ti(i)%t(j) * arr%ti(i)%t(j)
          enddo
       enddo
    case(TT_TILED_DIST,TT_TILED_REPL)
       norm=tensor_tiled_pdm_get_nrm2(arr)
    end select

    if(.not.squareback)norm = sqrt(norm)

    if(print_)print *,msg,norm

  end subroutine tensor_print_norm_customprint


  !> \brief Master routine for getting memory information in different shapes
  !> \autonr Patrick Ettenhuber
  subroutine tensor_print_mem_info(output,print_all_nodes,allaccess,reducetocheck)
    implicit none
    !> integer controling the output, if there should be any
    integer, intent(in) :: output
    !> optional logigal whether all nodes should print their information
    logical, intent(in), optional :: print_all_nodes
    !> optional logical if all nodes access this routine at the same time
    logical, intent(in), optional :: allaccess
    !> optional logigal stating that the informaion should be gathered on master
    integer, intent(inout), optional :: reducetocheck
    logical :: alln,red,master
    integer :: nnod, nod
    integer(kind=tensor_long_int),pointer :: red_info(:),narr(:)
    alln   = .false.
    red    = .false.
    master =.true.
    nnod   = 1
#ifdef VAR_MPI
    if(infpar%lg_mynum/=0)master=.false.
    nnod=infpar%lg_nodtot
#endif
    if(present(print_all_nodes))alln=print_all_nodes
    if(present(reducetocheck))then
      red=.true.
      call tensor_alloc_mem(red_info,nnod*9)
      call tensor_alloc_mem(narr,nnod)
      red_info = 0
      narr     = 0
    endif
    if(alln)then
      if(present(allaccess).and.red)call print_mem_per_node(output,allaccess,red_info,narr)
      if(.not.present(allaccess).and.red)call print_mem_per_node(output,.false.,red_info,narr)
      if(present(allaccess).and..not.red)call print_mem_per_node(output,allaccess)
      if(.not.present(allaccess).and..not.red)call print_mem_per_node(output,.false.)
    else
      call tensor_print_memory_currents(output)
    endif

    if(red)then
       !set test status
      reducetocheck=0
      do nod=1,nnod
         !check position 7 for each node, this ist the total memory currently
         !allocated in the tensor structure
         if(red_info(7+(nod-1)*9) /= 0) reducetocheck = 1
         if(narr(nod) /= 0)             reducetocheck = 1
      enddo
      call tensor_free_mem(red_info)
      call tensor_free_mem(narr)
    endif

  end subroutine tensor_print_mem_info

  subroutine print_norm_fort_wrapper1_nrm(fort,nelms,nrm,square)
    implicit none
    real(tensor_dp),intent(in) :: fort(:)
    integer(kind=8),intent(in) ::  nelms
    real(tensor_dp),intent(out),optional :: nrm
    logical,intent(in),optional :: square
    if(present(nrm))call print_norm_fort_nrm(fort,nelms,nrm)
    if(present(nrm).and.present(square))call print_norm_fort_nrm(fort,nelms,nrm,square)
    if(.not.present(nrm).and..not.present(square))call print_norm_fort_nrm(fort,nelms)
  end subroutine print_norm_fort_wrapper1_nrm
  subroutine print_norm_fort_wrapper2_nrm(fort,nelms,nrm,square)
    implicit none
    real(tensor_dp),intent(in) :: fort(:,:)
    integer(kind=8),intent(in) ::  nelms
    real(tensor_dp),intent(out),optional :: nrm
    logical,intent(in),optional :: square
    if(present(nrm))call print_norm_fort_nrm(fort,nelms,nrm)
    if(present(nrm).and.present(square))call print_norm_fort_nrm(fort,nelms,nrm,square)
    if(.not.present(nrm).and..not.present(square))call print_norm_fort_nrm(fort,nelms)
  end subroutine print_norm_fort_wrapper2_nrm
  subroutine print_norm_fort_wrapper3_nrm(fort,nelms,nrm,square)
    implicit none
    real(tensor_dp),intent(in) :: fort(:,:,:)
    integer(kind=8),intent(in) ::  nelms
    real(tensor_dp),intent(out),optional :: nrm
    logical,intent(in),optional :: square
    if(present(nrm))call print_norm_fort_nrm(fort,nelms,nrm)
    if(present(nrm).and.present(square))call print_norm_fort_nrm(fort,nelms,nrm,square)
    if(.not.present(nrm).and..not.present(square))call print_norm_fort_nrm(fort,nelms)
  end subroutine print_norm_fort_wrapper3_nrm
  subroutine print_norm_fort_wrapper4_nrm(fort,nelms,nrm,square)
    implicit none
    real(tensor_dp),intent(in) :: fort(:,:,:,:)
    integer(kind=8),intent(in) ::  nelms
    real(tensor_dp),intent(out),optional :: nrm
    logical,intent(in),optional :: square
    if(present(nrm))call print_norm_fort_nrm(fort,nelms,nrm)
    if(present(nrm).and.present(square))call print_norm_fort_nrm(fort,nelms,nrm,square)
    if(.not.present(nrm).and..not.present(square))call print_norm_fort_nrm(fort,nelms)
  end subroutine print_norm_fort_wrapper4_nrm
  subroutine print_norm_fort_nolen1_customprint(fort,msg)
    implicit none
    real(tensor_dp),intent(in) :: fort(:)
    character*(*),intent(in) :: msg
    integer(kind=8) ::  nelms
    nelms = size(fort)
    call print_norm_fort_customprint(fort,nelms,msg)
  end subroutine print_norm_fort_nolen1_customprint
  subroutine print_norm_fort_nolen2_customprint(fort,msg)
    implicit none
    real(tensor_dp),intent(in) :: fort(:,:)
    character*(*),intent(in) :: msg
    integer(kind=8) ::  nelms
    nelms = size(fort)
    call print_norm_fort_customprint(fort,nelms,msg)
  end subroutine print_norm_fort_nolen2_customprint
  subroutine print_norm_fort_nolen3_customprint(fort,msg)
    implicit none
    real(tensor_dp),intent(in) :: fort(:,:,:)
    character*(*),intent(in) :: msg
    integer(kind=8) ::  nelms
    nelms = size(fort)
    call print_norm_fort_customprint(fort,nelms,msg)
  end subroutine print_norm_fort_nolen3_customprint
  subroutine print_norm_fort_nolen4_customprint(fort,msg)
    implicit none
    real(tensor_dp),intent(in) :: fort(:,:,:,:)
    character*(*),intent(in) :: msg
    integer(kind=8) ::  nelms
    nelms = size(fort)
    call print_norm_fort_customprint(fort,nelms,msg)
  end subroutine print_norm_fort_nolen4_customprint
  subroutine print_norm_fort_wrapper1_customprint(fort,nelms,msg,square)
    implicit none
    real(tensor_dp),intent(in) :: fort(:)
    integer(kind=8),intent(in) ::  nelms
    character*(*),intent(in) :: msg
    logical,intent(in),optional :: square
    if(.not.present(square))call print_norm_fort_customprint(fort,nelms,msg)
    if(present(square))call print_norm_fort_customprint(fort,nelms,msg,square)
  end subroutine print_norm_fort_wrapper1_customprint
  subroutine print_norm_fort_wrapper2_customprint(fort,nelms,msg,square)
    implicit none
    real(tensor_dp),intent(in) :: fort(:,:)
    integer(kind=8),intent(in) ::  nelms
    character*(*),intent(in) :: msg
    logical,intent(in),optional :: square
    if(.not.present(square))call print_norm_fort_customprint(fort,nelms,msg)
    if(present(square))call print_norm_fort_customprint(fort,nelms,msg,square)
  end subroutine print_norm_fort_wrapper2_customprint
  subroutine print_norm_fort_wrapper3_customprint(fort,nelms,msg,square)
    implicit none
    real(tensor_dp),intent(in) :: fort(:,:,:)
    integer(kind=8),intent(in) ::  nelms
    character*(*),intent(in) :: msg
    logical,intent(in),optional :: square
    if(.not.present(square))call print_norm_fort_customprint(fort,nelms,msg)
    if(present(square))call print_norm_fort_customprint(fort,nelms,msg,square)
  end subroutine print_norm_fort_wrapper3_customprint
  subroutine print_norm_fort_wrapper4_customprint(fort,nelms,msg,square)
    implicit none
    real(tensor_dp),intent(in) :: fort(:,:,:,:)
    integer(kind=8),intent(in) ::  nelms
    character*(*),intent(in) :: msg
    logical,intent(in),optional :: square
    if(.not.present(square))call print_norm_fort_customprint(fort,nelms,msg)
    if(present(square))call print_norm_fort_customprint(fort,nelms,msg,square)
  end subroutine print_norm_fort_wrapper4_customprint

  subroutine array2_print_norm_nrm(arrtwo,nrm,square)
    implicit none
    type(array2),intent(in) :: arrtwo
    real(tensor_dp),intent(inout), optional :: nrm
    integer(kind=8) :: nelms
    logical,intent(in),optional :: square
    nelms = int(arrtwo%dims(1)*arrtwo%dims(2),kind=8)
    if(present(nrm).and..not.present(square))call print_norm_fort_nrm(arrtwo%val,nelms,nrm)
    if(present(nrm).and.present(square))call print_norm_fort_nrm(arrtwo%val,nelms,nrm,square)
    if(.not.present(nrm).and..not.present(square))call print_norm_fort_nrm(arrtwo%val,nelms)
  end subroutine array2_print_norm_nrm
  subroutine array2_print_norm_customprint(arrtwo,msg,square)
    implicit none
    type(array2),intent(in) :: arrtwo
    character*(*), intent(in) :: msg
    integer(kind=8) :: nelms
    logical,intent(in),optional :: square
    nelms = int(arrtwo%dims(1)*arrtwo%dims(2),kind=8)
    if(.not.present(square))call print_norm_fort_customprint(arrtwo%val,nelms,msg)
    if(present(square))call print_norm_fort_customprint(arrtwo%val,nelms,msg,square)
  end subroutine array2_print_norm_customprint
  subroutine array4_print_norm_nrm(arrf,nrm,square)
    implicit none
    type(array4),intent(in) :: arrf
    real(tensor_dp),intent(inout), optional :: nrm
    logical,intent(in),optional :: square
    integer(kind=8) :: nelms
    nelms = int(arrf%dims(1)*arrf%dims(2)*arrf%dims(3)*arrf%dims(4),kind=8)
    if(present(nrm).and..not.present(square))call print_norm_fort_nrm(arrf%val,nelms,nrm)
    if(present(nrm).and.present(square))call print_norm_fort_nrm(arrf%val,nelms,nrm,square)
    if(.not.present(nrm).and..not.present(square))call print_norm_fort_nrm(arrf%val,nelms)
  end subroutine array4_print_norm_nrm
  subroutine array4_print_norm_customprint(arrf,msg,square)
    implicit none
    type(array4),intent(in) :: arrf
    character*(*),intent(in):: msg
    logical,intent(in),optional :: square
    integer(kind=8) :: nelms
    nelms = int(arrf%dims(1)*arrf%dims(2)*arrf%dims(3)*arrf%dims(4),kind=8)
    if(.not.present(square))call print_norm_fort_customprint(arrf%val,nelms,msg)
    if(present(square))call print_norm_fort_customprint(arrf%val,nelms,msg,square)
  end subroutine array4_print_norm_customprint
  subroutine matrix_print_norm_nrm(mat,nrm,square)
    implicit none
    type(matrix),intent(in) :: mat
    real(tensor_dp),intent(inout), optional :: nrm
    logical,intent(in),optional :: square
    integer(kind=8) :: nelms
    nelms = int(mat%nrow*mat%ncol,kind=8)
    if(present(nrm).and..not.present(square))call print_norm_fort_nrm(mat%elms,nelms,nrm)
    if(present(nrm).and.present(square))call print_norm_fort_nrm(mat%elms,nelms,nrm,square)
    if(.not.present(nrm).and..not.present(square))call print_norm_fort_nrm(mat%elms,nelms)
  end subroutine matrix_print_norm_nrm

  subroutine print_norm_fort_nrm(fort,nelms,nrm,returnsquared)
    implicit none
    real(tensor_dp),intent(in) :: fort(*)
    integer(kind=8),intent(in) ::  nelms
    real(tensor_dp),intent(out),optional :: nrm
    logical,intent(in),optional :: returnsquared
    integer(kind=8) :: i
    real(tensor_dp) :: norm
    logical :: squareback
    squareback=.false.
    if(present(returnsquared))squareback=returnsquared
    norm=0.0E0_tensor_dp
    do i=1,nelms
      norm = norm + fort(i) * fort(i)
    enddo
    if(.not.squareback)norm = sqrt(norm)
    if(.not.present(nrm).and..not.squareback)print *,"NORM:",norm
    if(.not.present(nrm).and.squareback)print *,"NORM^2:",norm
    if(present(nrm))nrm=norm
  end subroutine print_norm_fort_nrm
  subroutine print_norm_fort_customprint(fort,nelms,string,returnsquared)
    implicit none
    real(tensor_dp),intent(in) :: fort(*)
    integer(kind=8),intent(in) ::  nelms
    character*(*),intent(in) :: string
    logical,intent(in),optional :: returnsquared
    integer(kind=8) :: i
    real(tensor_dp) :: norm
    logical :: squareback
    squareback=.false.
    if(present(returnsquared))squareback=returnsquared
    norm=0.0E0_tensor_dp
    do i=1,nelms
      norm = norm + fort(i) * fort(i)
    enddo
    if(.not.squareback)norm = sqrt(norm)
    print *,string,norm
  end subroutine print_norm_fort_customprint



  subroutine tensor_scale(arr,sc)
    implicit none
    type(tensor) :: arr
    real(tensor_dp) :: sc
    
    select case(arr%itype)
    case(TT_DENSE)
      call dscal(int(arr%nelms),sc,arr%elm1,1)
    case(TT_REPLICATED)
      call dscal(int(arr%nelms),sc,arr%elm1,1)
      call tensor_sync_replicated(arr)
    case(TT_TILED_DIST)
      call tensor_scale_td(arr,sc)
    case default
      call lsquit("ERROR(tensor_scale):not yet implemented",DECinfo%output)
    end select
  end subroutine tensor_scale
  subroutine get_symm_tensor_segmenting_simple(nnodes,a,b,a_seg,b_seg)
     implicit none
     integer, intent(in)  :: a,b,nnodes
     integer, intent(out) :: a_seg, b_seg
     integer :: counter
     integer :: modtilea, modtileb
     real(tensor_dp) :: max_mem_p_tile_in_GB
     !get segmenting for tensors, divide dimensions until tiles are less than
     !100MB and/or until enough tiles are available such that each node gets at
     !least one and as long as a_seg>=2 and b_seg>=2
     b_seg    = b
     a_seg    = a
     modtilea = 0
     modtileb = 0

     if( tensor_segment_length_set )then
        b_seg    = tensor_segment_length
        a_seg    = tensor_segment_length
     else
        max_mem_p_tile_in_GB = DECinfo%cc_solver_tile_mem

        select case(DECinfo%tensor_segmenting_scheme)
        case (1)

           counter  = 1

           !FIRST a then b

           do while(   ( ( b_seg**2*a_seg**2)*8.0E0_tensor_dp/(1024.0E0_tensor_dp**3) > max_mem_p_tile_in_GB &
                 & .or.((b/b_seg+modtileb)**2*(a/a_seg+modtilea)**2<nnodes)                  )&
                 & .and. (b_seg>=1.or.a_seg>=1) .and. (a/a_seg+modtilea) <= 4 )

              a_seg = a / counter + mod(a,counter)

              counter = counter + 1

              modtilea = 0
              if(mod(a,a_seg)/=0)modtilea = 1

           enddo

           counter  = 1

           do while(   ( ( b_seg**2*a_seg**2)*8.0E0_tensor_dp/(1024.0E0_tensor_dp**3) > max_mem_p_tile_in_GB &
                 &  .or. ((b/b_seg+modtileb)**2*(a/a_seg+modtilea)**2<nnodes)      )&
                 & .and. (b_seg>=1.or.a_seg>=1)  .and. (b/b_seg+modtileb) <= 4   )

              b_seg = b / counter + mod(b,counter)

              counter = counter + 1

              modtileb = 0
              if(mod(b,b_seg)/=0)modtileb = 1

           enddo

        case(2)

           !FIRST b then a
           counter  = 1

           do while(   ( ( b_seg**2*a_seg**2)*8.0E0_tensor_dp/(1024.0E0_tensor_dp**3) > max_mem_p_tile_in_GB &
                 &  .or. ((b/b_seg+modtileb)**2*(a/a_seg+modtilea)**2<nnodes)      )&
                 & .and. (b_seg>=1.or.a_seg>=1)   .and. (b/b_seg+modtileb) <= 4   )

              b_seg = b / counter + mod(b,counter)

              counter = counter + 1

              modtileb = 0
              if(mod(b,b_seg)/=0)modtileb = 1

           enddo

           counter  = 1

           do while(   ( ( b_seg**2*a_seg**2)*8.0E0_tensor_dp/(1024.0E0_tensor_dp**3) > max_mem_p_tile_in_GB &
                 &  .or. ((b/b_seg+modtileb)**2*(a/a_seg+modtilea)**2<nnodes)      )&
                 & .and. (b_seg>=1.or.a_seg>=1)  .and. (a/a_seg+modtilea) <= 4  )

              a_seg = a / counter + mod(a,counter)

              counter = counter + 1

              modtilea = 0
              if(mod(a,a_seg)/=0)modtilea = 1

           enddo
        case (3)

           counter  = 1

           !FIRST a then b

           do while(   ( ( b_seg**2*a_seg**2)*8.0E0_tensor_dp/(1024.0E0_tensor_dp**3) > max_mem_p_tile_in_GB &
                 & .or.((b/b_seg+modtileb)**2*(a/a_seg+modtilea)**2<nnodes)                  )&
                 & .and. (b_seg>=1.or.a_seg>=1) .and. (a/a_seg+modtilea) <= 4 )

              a_seg = min(a,b) / counter + mod(min(a,b),counter)
              b_seg = a_seg

              counter = counter + 1

              modtilea = 0
              if(mod(a,a_seg)/=0)modtilea = 1

           enddo

        case (4)

           counter  = 1

           !BOTH a and b
           do while(   ( ( b_seg**2*a_seg**2)*8.0E0_tensor_dp/(1024.0E0_tensor_dp**3) > max_mem_p_tile_in_GB &
                 &  .or. ((b/b_seg+modtileb)**2*(a/a_seg+modtilea)**2<nnodes)      )&
                 & .and. (b_seg>=10.or.a_seg>=10)  .and. (a/a_seg+modtilea) <= 4 .and. (b/b_seg+modtileb) <= 4  )

              b_seg = b / counter + mod(b,counter)
              a_seg = a / counter + mod(a,counter)

              counter = counter + 1

              modtilea = 0
              if(mod(a,a_seg)/=0)modtilea = 1

              modtileb = 0
              if(mod(b,b_seg)/=0)modtileb = 1

           enddo

           !then make sure that pure virtual batches have a size < thr

           do while(   ( ( b_seg**4)*8.0E0_tensor_dp/(1024.0E0_tensor_dp**3) > max_mem_p_tile_in_GB &
                 &  .or. ((b/b_seg+modtileb)**4 < nnodes) )&
                 & .and. ( b_seg>=10 )  )

              b_seg = b / counter + mod(b,counter)

              counter = counter + 1

              modtileb = 0
              if(mod(b,b_seg)/=0)modtileb = 1

           enddo

        case default

           counter  = 1

           !BOTH a and b
           do while(   ( ( b_seg**2*a_seg**2)*8.0E0_tensor_dp/(1024.0E0_tensor_dp**3) > max_mem_p_tile_in_GB &
                 &  .or. ((b/b_seg+modtileb)**2*(a/a_seg+modtilea)**2<nnodes)      )&
                 & .and. (b_seg>=2.or.a_seg>=2)  .and. (a/a_seg+modtilea) <= 4 .and. (b/b_seg+modtileb) <= 4  )

              b_seg = b / counter + mod(b,counter)
              a_seg = a / counter + mod(a,counter)

              counter = counter + 1

              modtilea = 0
              if(mod(a,a_seg)/=0)modtilea = 1

              modtileb = 0
              if(mod(b,b_seg)/=0)modtileb = 1

           enddo

        end select
        endif

     if(DECinfo%PL>2)then
        print *,"SPLITTING OF DIMS IN A^2B^2"
        print *,"#a",a,"seg:",a_seg
        print *,"#b",b,"seg:",b_seg
     endif

  end subroutine get_symm_tensor_segmenting_simple

end module tensor_interface_module

