!> @file
!> This is the file where all higher order pdm operations should go, especially,
!everything that has to do with arrays, allocating, feeing and so on, the one
!sided wrappers should go to lspdm_base_module
!> \author Patrick Ettenhuber
!> \date April 2013
module lspdm_tensor_operations_module
  use,intrinsic :: iso_c_binding,only:c_f_pointer,c_loc


!#ifdef VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN
!  use dec_workarounds_module
!#endif
#ifdef TENSORS_IN_LSDALTON
  use LSTIMING
#endif

  use get_idx_mod
  use tensor_parameters_and_counters
  use tensor_mpi_binding_module
  use tensor_bg_buf_module
  use tensor_mpi_interface_module
  use tensor_mpi_operations_module

  use tensor_basic_module
  use lspdm_basic_module


  !INTERFACES
  !**********
  interface tensor_get_tile
     module procedure tensor_gett44,&
        &tensor_gett48,&
        &tensor_gett84,&
        &tensor_gett88,&
        &tensor_gettile_modeidx
  end interface tensor_get_tile

  interface tensor_put_tile
     module procedure tensor_putt44,&
        &tensor_putt48,&
        &tensor_putt84,&
        &tensor_putt88,&
        &tensor_puttile_modeidx
  end interface tensor_put_tile

  interface tensor_acc_tile
     module procedure tensor_acct44,&
        &tensor_acct48,&
        &tensor_acct84,&
        &tensor_acct88,&
        &tensor_acctile_modeidx
  end interface tensor_acc_tile

  interface tensor_lock_win
     module procedure tensor_lock_win8,tensor_lock_win4
  end interface tensor_lock_win

#ifdef COMPILER_UNDERSTANDS_FORTRAN_2003
  abstract interface
  subroutine put_acc_tile(arr,globtilenr,fort,nelms,lock_set,flush_it,req)
     use tensor_parameters_and_counters
     import
     implicit none
     type(tensor),intent(in) :: arr
     integer(kind=tensor_long_int),intent(in) :: globtilenr
     integer(kind=tensor_long_int),intent(in) :: nelms
     real(tensor_dp),intent(inout) :: fort(nelms)
     logical, optional, intent(in) :: lock_set,flush_it
     integer(kind=tensor_mpi_kind),intent(inout), optional :: req
  end subroutine put_acc_tile
  end interface
#endif

  !interface tensor_accumulate_tile_nobuff
  !  module procedure tensor_accumulate_tile_combidx4_nobuff,&
  !                  &tensor_accumulate_tile_combidx8_nobuff,&
  !                  &tensor_accumulate_tile_modeidx_nobuff
  !end interface tensor_accumulate_tile_nobuff


  !Persistent array type definition
  !--------------------------------
  !> \brief the persistent array is a collection of n=500 arrays on each node
  !with some additional information. Here the storage fort tiled distributed and
  !replicated arrays is allocated, if 500  is not enough, please change here
  !> amount of arrays which are storable in the persistent array
  integer, parameter :: n_arrays = 500
  !> persistent array type-def
  type persistent_array
     !> collection of arrays
     type(tensor),pointer :: a(:)
     !> current address on node
     integer :: curr_addr_on_node  = 1
     !> counter for how many arrays were allocated in the persisten array
     integer :: arrays_allocated   = 0
     !> counter for how many arrays were deallocated
     integer :: arrays_deallocated = 0
     !> conter for the arrays currently in use
     integer :: arrays_in_use      = 0
     !> offset for the first tile allocation to get a better load distribution
     integer :: new_offset         = 0
     !> list of n logicals as indicator wheter an adress is free to allocate a
     !new array
     logical,pointer :: free_addr_on_node(:) => null()
  endtype persistent_array

  ! Global communication buffer in this module
  type global_module_buffer
     real(tensor_dp), pointer :: buf(:)
     integer(kind=tensor_long_int)      :: n    = 0
     logical              :: init = .false.
  end type global_module_buffer

  save

  ! job parameters for pdm jobs
  integer,parameter :: JOB_TEST_FRAMEWORK           =  1
  integer,parameter :: JOB_BARRIER                  =  2
  integer,parameter :: JOB_FREE_TENSOR_STD          =  3
  integer,parameter :: JOB_INIT_TENSOR_TILED        =  4
  integer,parameter :: JOB_FREE_TENSOR_PDM          =  5
  integer,parameter :: JOB_INIT_TENSOR_REPLICATED   =  6
  integer,parameter :: JOB_PRINT_MEM_INFO1          =  7
  integer,parameter :: JOB_PRINT_MEM_INFO2          =  8
  integer,parameter :: JOB_GET_NRM2_TILED           =  9
  integer,parameter :: JOB_DATA2TILED_DIST          = 10
  integer,parameter :: JOB_GET_TILE_SEND            = 11
  integer,parameter :: JOB_PRINT_TI_NRM             = 12
  integer,parameter :: JOB_SYNC_REPLICATED          = 13
  integer,parameter :: JOB_GET_NORM_REPLICATED      = 14
  integer,parameter :: JOB_PREC_DOUBLES_PAR         = 15
  integer,parameter :: JOB_DDOT_PAR                 = 16
  integer,parameter :: JOB_ADD_PAR                  = 17
  integer,parameter :: JOB_CP_ARR                   = 18
  integer,parameter :: JOB_TENSOR_ZERO              = 19
  integer,parameter :: JOB_GET_CC_ENERGY            = 20
  integer,parameter :: JOB_GET_FRAG_CC_ENERGY       = 21
  integer,parameter :: JOB_CHANGE_ACCESS_TYPE       = 22
  integer,parameter :: JOB_TENSOR_SCALE             = 23

  integer,parameter :: JOB_GET_MP2_ENERGY           = 25
  integer,parameter :: JOB_GET_RPA_ENERGY           = 26
  integer,parameter :: JOB_GET_SOS_ENERGY           = 27
  integer,parameter :: JOB_TENSOR_CONTRACT_SIMPLE   = 28
  integer,parameter :: JOB_TENSOR_CONTRACT_BDENSE   = 29
  integer,parameter :: JOB_TENSOR_EXTRACT_VEOS      = 30
  integer,parameter :: JOB_TENSOR_EXTRACT_OEOS      = 31
  integer,parameter :: JOB_TENSOR_EXTRACT_ODECNP    = 32
  integer,parameter :: JOB_TENSOR_EXTRACT_VDECNP    = 33
  integer,parameter :: JOB_GET_COMBINEDT1T2_1       = 34
  integer,parameter :: JOB_GET_COMBINEDT1T2_2       = 35
  integer,parameter :: JOB_GET_MP2_ST_GUESS         = 36
  integer,parameter :: JOB_tensor_rand              = 37
  integer,parameter :: JOB_HMUL_PAR                 = 38
  integer,parameter :: JOB_DMUL_PAR                 = 39
  integer,parameter :: JOB_LSPDM_INIT_GLOBAL_BUFFER = 40
  integer,parameter :: JOB_LSPDM_FREE_GLOBAL_BUFFER = 41
  integer,parameter :: JOB_SET_TENSOR_SEG_LENGTH    = 42
  integer,parameter :: JOB_SET_TENSOR_DEBUG_TRUE    = 43
  integer,parameter :: JOB_SET_TENSOR_ALWAYS_SYNC_TRUE = 44
  integer,parameter :: JOB_SET_TENSOR_BACKEND_TRUE  = 45

  !> definition of the persistent array 
  type(persistent_array) :: p_arr

  !define the buffer
  type(global_module_buffer) :: gm_buf

  !> timing and measuring variables
  real(tensor_dp) :: time_pdm_acc          = 0.0E0_tensor_dp
  integer(kind=tensor_long_int) :: bytes_transferred_acc = 0
  integer(kind=tensor_long_int) :: nmsg_acc = 0
  real(tensor_dp) :: time_pdm_put          = 0.0E0_tensor_dp
  integer(kind=tensor_long_int) :: bytes_transferred_put = 0
  integer(kind=tensor_long_int) :: nmsg_put = 0
  real(tensor_dp) :: time_pdm_get          = 0.0E0_tensor_dp
  integer(kind=tensor_long_int) :: bytes_transferred_get = 0
  integer(kind=tensor_long_int) :: nmsg_get = 0

#ifdef VAR_MPI
#ifdef COMPILER_UNDERSTANDS_FORTRAN_2003
  procedure(put_acc_tile),pointer :: acc_tlong1
  procedure(put_acc_tile),pointer :: put_tlong1 
#endif
#endif

  contains

  !> \brief intitialize storage room for the tiled distributed arrays
  !> \author Patrick Ettenhuber
  !> \date May 2013
  subroutine init_persistent_array()
     implicit none
     integer :: i

     call tensor_alloc_mem(p_arr%a,n_arrays)
     call tensor_alloc_mem(p_arr%free_addr_on_node,n_arrays)


     do i = 1, n_arrays
        p_arr%free_addr_on_node(i) = .true.
        call tensor_reset_value_defaults(p_arr%a(i)) 
        call tensor_nullify_pointers(p_arr%a(i)) 
     end do

     !if( lspdm_use_comm_proc ) call tensor_status_quit("ERROR(init_persistent_array)&
     !& lspdm_use_comm_proc cannot be true at startup",-1)
     lspdm_use_comm_proc = .false.

     !set defaults for comm buffer
     gm_buf%n    = 0
     gm_buf%init = .false.
  end subroutine init_persistent_array
  !>  \brief free storage room for the tiled distributed arrays
  !> \author Patrick Ettenhuber
  !> \date May 2013
  subroutine free_persistent_array()
     implicit none
     integer :: i
     if(associated(p_arr%a))then

        do i = 1, n_arrays
           call tensor_reset_value_defaults(p_arr%a(i)) 
           call tensor_nullify_pointers(p_arr%a(i)) 
        end do

        call tensor_free_mem(p_arr%a)

     endif

     if(associated(p_arr%free_addr_on_node))then
        call tensor_free_mem(p_arr%free_addr_on_node)
     endif

     if( lspdm_use_comm_proc ) call tensor_status_quit("ERROR(free_persistent_array) &
        & lspdm_use_comm_proc has to be disabled at shutdown, otherwise there &
        & still might be processes running",-1)
  end subroutine free_persistent_array

  subroutine new_group_reset_persistent_array
     implicit none
     p_arr%new_offset = 0
  end subroutine new_group_reset_persistent_array


  ! initalizing the global buffer array
  subroutine lspdm_init_global_buffer(call_slaves_from_slaveroutine)
     implicit none
     logical, intent(in) :: call_slaves_from_slaveroutine
     integer(kind=tensor_mpi_kind) :: me
     integer :: i, checked
     integer(kind=tensor_long_int) :: nelms

     me = 0
#ifdef VAR_MPI
     call tensor_get_rank(me)

     if( me == 0 .and. call_slaves_from_slaveroutine )then
        call pdm_tensor_sync(JOB_LSPDM_INIT_GLOBAL_BUFFER)
     endif

#endif

     if( gm_buf%init )then
        call tensor_status_quit("ERROR(lspdm_init_global_buffer): global buffer already initialized",-1)
     endif

     !Find the buffer size based on the allocated arrays
     checked  = 0
     nelms    = 0
     gm_buf%n = 0
     LoopAllAllocd: do i=1,n_arrays

        if(.not.p_arr%free_addr_on_node(i))then
           if(p_arr%a(i)%itype==TT_TILED_DIST.or.p_arr%a(i)%itype==TT_TILED_REPL)then

              nelms = max(nelms,long1*2*p_arr%a(i)%tsize)

              !counter for fast exit
              checked = checked + 1
           endif
        endif

        if(checked>=p_arr%arrays_in_use) exit LooPAllAllocd

     enddo LoopAllAllocd

     call lspdm_reinit_global_buffer(nelms)

     gm_buf%init = .true.

  end subroutine lspdm_init_global_buffer

  subroutine lspdm_reinit_global_buffer(nelms)
     implicit none
     integer(kind=tensor_long_int), intent(in) :: nelms

     if(nelms > gm_buf%n)then

        if((8.0E0_tensor_dp*nelms)/(1024.0**2) > 500.0E0_tensor_dp)then
           print *,"WARNING(lspdm_reinit_global_buffer): background buffer more than&
              & 500MB. Smaller tiles or switching off the background buffer prevent&
              & this warning"
        endif

        if( associated(gm_buf%buf) ) call tensor_free_mem(gm_buf%buf)

        gm_buf%n = nelms

        call tensor_alloc_mem(gm_buf%buf,gm_buf%n)

     endif

  end subroutine lspdm_reinit_global_buffer

  ! freeing the global buffer array
  subroutine lspdm_free_global_buffer(call_slaves_from_slaveroutine)
     implicit none
     logical, intent(in) :: call_slaves_from_slaveroutine
     integer(kind=tensor_mpi_kind) :: me

     me = 0
#ifdef VAR_MPI
     call tensor_get_rank(me)

     if( me == 0 .and. call_slaves_from_slaveroutine )then
        call pdm_tensor_sync(JOB_LSPDM_FREE_GLOBAL_BUFFER)
     endif

#endif

     if( .not. gm_buf%init )then
        call tensor_status_quit("ERROR(lspdm_free_global_buffer): global buffer not initialized",-1)
     endif

     if(associated(gm_buf%buf))then
        call tensor_free_mem(gm_buf%buf)
     endif

     gm_buf%n    = 0
     gm_buf%init = .false.
  end subroutine lspdm_free_global_buffer


  subroutine set_signal_for_slaves(signal_idx)
     implicit none
     integer, intent(in) :: signal_idx
     TENSOR_SLAVES_TO_SLAVE_ROUTINE = signal_idx
  end subroutine set_signal_for_slaves

  !> \brief main subroutine for the communication of nodes on grid handling 
  !arr structures. this routine assumes, that always the process with rank 0 in
  !the given communicator is the "manager"
  !> \author Patrick Ettenhuber
  !> \date May 2012
  subroutine pdm_tensor_sync(job,a,b,c,d)
     implicit none
     !> job is input for master and output for slaves, the arguments have to be
     !in the job paramenters list in top of this file
     integer                          :: job
     !> the communicator on which this routine should work
     !integer(kind=tensor_mpi_kind),intent(in) :: comm
     !the array(s) to be passed to the slaves for which the operation is
     !performed
     type(tensor),optional             :: a,b,c,d

     !> comm arrays
     integer,pointer                  :: TMPI(:), dims(:)
     !character :: TMPC(12)
     integer :: i, j, context,modes(3),counter, stat,ierr,basic
     integer(kind=tensor_mpi_kind)            :: sendctr,me,nn, comm
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     logical                          :: loc
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif

     modes=0
#ifdef VAR_MPI
     comm = tensor_work_comm

     call tensor_get_rank(me)
     call tensor_get_size(nn)

     basic = 12

     loc = .false.
     if(loc) call tensor_status_quit("ERROR(pdm_tensor_sync): this feature has been deactivated",-1)

     if( me == root) then
        !**************************************************************************************
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!code for MASTER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !**************************************************************************************
        !Wake up slaves
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_COMM )
#endif
        call tensor_mpi_bcast( TENSOR_SLAVES_TO_SLAVE_ROUTINE, me, tensor_work_comm)
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_WORK )
#endif
        !1     = JOB
        !2-5   = address in slot a-c
        !5-8   = modes a-c
        !9-13  = zero -> bool passed as integer for a-c
        !rest specifies dimensions
        counter = basic
        if (present(A)) then
           counter     = counter+2*A%mode
        endif
        if (present(B)) then
           counter     = counter+2*B%mode
        endif
        if (present(C)) then
           counter     = counter+2*C%mode
        endif
        if (present(D)) then
           counter     = counter+2*D%mode
        endif

#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_COMM )
#endif
        call tensor_mpi_bcast(counter,root,comm)
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_WORK )
#endif

        call tensor_alloc_mem(TMPI,counter)

        !change counter and basic for checking in the end
        TMPI(1) = counter
        counter = basic
        basic   = TMPI(1) 

        !get comm vector done
        TMPI    = 0
        TMPI(1) = job
        if (present(A)) then
           TMPI(6)                        = A%mode
           if(A%zeros) TMPI(10)           = 1
           TMPI(counter+1:counter+A%mode) = A%dims
           counter = counter + A%mode
           TMPI(counter+1:counter+A%mode) = A%tdim
           counter = counter + A%mode
        endif
        if (present(B)) then
           TMPI(7)                        = B%mode
           if(B%zeros) TMPI(11)            = 1
           TMPI(counter+1:counter+B%mode) = B%dims
           counter = counter + B%mode
           TMPI(counter+1:counter+B%mode) = B%tdim
           counter = counter + B%mode
        endif
        if (present(C)) then
           TMPI(8)                        = C%mode
           if(C%zeros) TMPI(12)           = 1
           TMPI(counter+1:counter+C%mode) = C%dims
           counter = counter + C%mode
           TMPI(counter+1:counter+C%mode) = C%tdim
           counter = counter + C%mode
        endif
        if (present(D)) then
           TMPI(10)                       = D%mode
           if(D%zeros) TMPI(13)           = 1
           TMPI(counter+1:counter+D%mode) = D%dims
           counter = counter + D%mode
           TMPI(counter+1:counter+D%mode) = D%tdim
           counter = counter + D%mode
        endif

        if(counter/=basic)call tensor_status_quit("ERROR(pdm_tensor_sync):different number of&
           & elements for MASTER",lspdm_stdout)
        !if(loc)then
        !if(nn>1.and.present(A).and..not.associated(A%addr_loc))&
        !&call tensor_status_quit("ERROR(pdm_tensor_sync):addr_loc for array A not associated",lspdm_stdout)
        !if(nn>1.and.present(B).and..not.associated(B%addr_loc))&
        !&call tensor_status_quit("ERROR(pdm_tensor_sync):addr_loc for array B not associated",lspdm_stdout)
        !if(nn>1.and.present(C).and..not.associated(C%addr_loc))&
        !&call tensor_status_quit("ERROR(pdm_tensor_sync):addr_loc for array C not associated",lspdm_stdout)
        !if(nn>1.and.present(D).and..not.associated(D%addr_loc))&
        !&call tensor_status_quit("ERROR(pdm_tensor_sync):addr_loc for array D not associated",lspdm_stdout)
        !else
        if(nn>1.and.present(A))then
           if(.not.associated(A%addr_p_arr))&
              &call tensor_status_quit("ERROR(pdm_tensor_sync):addr_p_arr for array A not associated",lspdm_stdout)
        endif
        if(nn>1.and.present(B))then
           if(.not.associated(B%addr_p_arr))&
              &call tensor_status_quit("ERROR(pdm_tensor_sync):addr_p_arr for array B not associated",lspdm_stdout)
        endif
        if(nn>1.and.present(C))then
           if(.not.associated(C%addr_p_arr))&
              &call tensor_status_quit("ERROR(pdm_tensor_sync):addr_p_arr for array C not associated",lspdm_stdout)
        endif
        if(nn>1.and.present(D))then
           if(.not.associated(D%addr_p_arr))&
              &call tensor_status_quit("ERROR(pdm_tensor_sync):addr_p_arr for array D not associated",lspdm_stdout)
        endif
        !endif


        do sendctr=1,nn-1
           if (present(A)) TMPI(2)  = A%addr_p_arr(sendctr+1)
           if (present(B)) TMPI(3)  = B%addr_p_arr(sendctr+1)
           if (present(C)) TMPI(4)  = C%addr_p_arr(sendctr+1)
           if (present(D)) TMPI(5)  = D%addr_p_arr(sendctr+1)

#ifdef TENSORS_IN_LSDALTON
           call time_start_phase( PHASE_COMM )
#endif
           call tensor_mpi_sendrecv( TMPI, counter, comm, root, sendctr)
#ifdef TENSORS_IN_LSDALTON
           call time_start_phase( PHASE_WORK )
#endif
        enddo
        call tensor_free_mem(TMPI)


     else  

        !**************************************************************************************
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!code for SLAVES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !**************************************************************************************
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_COMM )
#endif
        call tensor_mpi_bcast( counter, root, comm )
        call tensor_alloc_mem( TMPI, counter )
        call tensor_mpi_sendrecv( TMPI, counter, comm, root, me)
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_WORK )
#endif

        !get data from info vector; THIS COUNTER CONSTRUCTION HAS TO BE REWRITTEN
        !if NEEDED FOR NOW IT IS CONVENIENT, BECAUSE IT IS SIMPLE
        counter = basic
        job = TMPI(1) !slaves needs to know what to do
        !1     = JOB
        !2-5   = address in slot a-c
        !6-9   = modes a-c
        !10-13 = zero -> bool passed as integer for a-c
        !rest specifies dimensions
        if (TMPI(2).gt.0) then
           A = p_arr%a(TMPI(2))
        else
           if(TMPI(6).gt.0)then
              A%mode                = TMPI(6)
              if(TMPI(10)==1) A%zeros= .true.
              call tensor_set_dims(A,TMPI(counter+1:counter+A%mode),int(A%mode))
              counter = counter + A%mode
              call tensor_set_tdims(A,TMPI(counter+1:counter+A%mode),int(A%mode))
              counter = counter + A%mode
           endif
        endif
        if (TMPI(3).gt.0) then
           B = p_arr%a(TMPI(3))
        else
           if(TMPI(7).gt.0)then
              B%mode                = TMPI(7)
              if(TMPI(11)==1)B%zeros= .true.
              call tensor_set_dims(B,TMPI(counter+1:counter+B%mode),int(B%mode))
              counter = counter + B%mode
              call tensor_set_tdims(B,TMPI(counter+1:counter+B%mode),int(B%mode))
              counter = counter + B%mode
           endif
        endif
        if (TMPI(4).gt. 0) then
           C = p_arr%a(TMPI(4))
        else
           if(TMPI(8).gt.0)then
              C%mode                = TMPI(8)
              if(TMPI(12)==1)C%zeros= .true.
              call tensor_set_dims(C,TMPI(counter+1:counter+C%mode),int(C%mode))
              counter = counter + C%mode
              call tensor_set_tdims(C,TMPI(counter+1:counter+C%mode),int(C%mode))
              counter = counter + C%mode
           endif
        endif
        if (TMPI(5).gt. 0) then
           !C = associate_to_p_arr(TMPI(4))
           D = p_arr%a(TMPI(5))
        else
           if(TMPI(9).gt.0)then
              D%mode                = TMPI(9)
              if(TMPI(13)==1)D%zeros= .true.
              call tensor_set_dims(D,TMPI(counter+1:counter+D%mode),int(D%mode))
              counter = counter + D%mode
              call tensor_set_tdims(D,TMPI(counter+1:counter+D%mode),int(D%mode))
              counter = counter + D%mode
           endif
        endif
        call tensor_free_mem(TMPI)
     endif

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
#endif
  end subroutine pdm_tensor_sync

  !> \author Patrick Ettenhuber
  !> \date December 2012
  !> \brief get an array from the persistent array by specifying its address
  function get_tensor_from_parr(addr) result(arr)
     implicit none
     !> the address of the array to extract
     integer,intent(in) :: addr
     !> array extracted from persisten array 
     type(tensor) :: arr
     arr=p_arr%a(addr)
  end function get_tensor_from_parr



  !> \author Patrick Ettenhuber
  !> \date December 2012
  !> \brief calculate fragment eos cc energy in parallel (PDM)
  function get_fragment_cc_energy_parallel(t1,t2,gmo,occ_num,virt_num,occ_idx,virt_idx) result(fEc)
     implicit none
     !> singles amplitudes
     type(tensor), intent(inout) :: t1
     !> two electron integrals in the mo-basis
     type(tensor), intent(inout) :: gmo
     !> doubles amplitudes
     type(tensor), intent(in) :: t2
     !> number of occupied indices
     integer, intent(in) :: occ_num
     !> number of virtual indices
     integer, intent(in) :: virt_num
     !> referencing the occupied indices of the fragment to the full basis
     integer, intent(in) :: occ_idx(occ_num)
     !> referencing the virtueal indices of the fragment to the full basis
     integer, intent(in) :: virt_idx(virt_num)
     !> return-calue fEc contains the fragment correlation energy
     real(tensor_dp) :: Evirt,Eocc,fEc
     real(tensor_dp),pointer :: t(:,:,:,:)
     integer :: lt,i,j,a,b,o(t2%mode),fr_i,fr_j,fr_a,fr_b
     integer :: i_high,j_high,a_high,b_high
     logical :: use_bg
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: rank, nnod

     Eocc  = 0.0E0_tensor_dp
     Evirt = 0.0E0_tensor_dp
     fEc   = 0.0E0_tensor_dp

     call tensor_get_rank(rank)
     call tensor_get_size(nnod)

#ifdef VAR_MPI
     !Get the slaves to this routine
     if(rank == root)then
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase(PHASE_COMM)
#endif

        call pdm_tensor_sync(JOB_GET_FRAG_CC_ENERGY,t1,t2,gmo)

        call tensor_buffer(occ_num,root=root,comm=tensor_work_comm)
        call tensor_buffer(occ_idx,occ_num)
        call tensor_buffer(virt_num)
        call tensor_buffer(virt_idx,virt_num,finalize=.true.)

#ifdef TENSORS_IN_LSDALTON
        call time_start_phase(PHASE_WORK)
#endif
     endif
     use_bg = (tensor_bg_init()) .and. (tensor_bg_free() > gmo%nelms)
     call memory_allocate_tensor_dense(gmo,use_bg)

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_COMM)
#endif
     call cp_tileddata2fort(gmo,gmo%elm1,gmo%nelms,.true.)
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif


     do lt=1,t2%nlti
#ifdef VAR_PTR_RESHAPE
        t(1:t2%ti(lt)%d(1),1:t2%ti(lt)%d(2),1:t2%ti(lt)%d(3),1:t2%ti(lt)%d(4)) => t2%ti(lt)%t
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
        call c_f_pointer(c_loc(t2%ti(lt)%t(1)),t,t2%ti(lt)%d)
#else
        call tensor_status_quit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
        !get offset for global indices
        call get_midx(t2%ti(lt)%gt,o,t2%ntpm,t2%mode)

        do j=1,t2%mode
           o(j)=(o(j)-1)*t2%tdim(j)
        enddo

        !find the limits of the current tile
        a_high= t2%ti(lt)%d(1)
        b_high= t2%ti(lt)%d(2)
        i_high= t2%ti(lt)%d(3)
        j_high= t2%ti(lt)%d(4)


        ! Energy using occupied scheme
        ! ****************************
        do j=1,occ_num
           fr_j = occ_idx(j)-o(4)   ! occupied EOS index in occupied AOS list
           if(j_high>=fr_j.and.fr_j>0)then 
              do i=1,occ_num
                 fr_i = occ_idx(i)-o(3)
                 if(i_high>=fr_i.and.fr_i>0)then
                    do b=1,t2%ti(lt)%d(2)
                       do a=1,t2%ti(lt)%d(1)

                          Eocc = Eocc + &
                             & ( t(a,b,fr_i,fr_j) + t1%elm2(a+o(1),fr_i+o(3))*t1%elm2(b+o(2),fr_j+o(4)) ) * &
                             & ( 2.0E0_tensor_dp*gmo%elm4(fr_i+o(3),a+o(1),fr_j+o(4),b+o(2)) &
                             & - gmo%elm4(fr_i+o(3),b+o(2),fr_j+o(4),a+o(1)) )

                       end do
                    end do
                 endif
              end do
           endif
        end do

        ! Energy using virtual scheme
        do a=1,virt_num
           fr_a = virt_idx(a)-o(1)
           if(a_high>=fr_a.and.fr_a>0)then
              do j=1,t2%ti(lt)%d(4)
                 do b=1,virt_num
                    fr_b = virt_idx(b)-o(2)  ! virtual EOS index in occupied AOS list
                    if(b_high>=fr_b.and.fr_b>0)then
                       do i=1,t2%ti(lt)%d(3)

                          Evirt = Evirt + &
                             & ( t(fr_a,fr_b,i,j) + t1%elm2(fr_a+o(1),i+o(3))*t1%elm2(fr_b+o(2),j+o(4)) ) * &
                             & ( 2.0E0_tensor_dp*gmo%elm4(i+o(3),fr_a+o(1),j+o(4),fr_b+o(2)) &
                             &- gmo%elm4(i+o(3),fr_b+o(2),j+o(4),fr_a+o(1)) )

                       end do
                    endif
                 end do
              end do
           endif
        end do

        nullify(t)

     enddo

     ! Hybrid scheme: Mixture of occupied and virtual partitioning schemes
     ! Ehybrid = 1/2* (Eocc + Evirt)

     call tensor_deallocate_dense(gmo)

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_COMM)
#endif
     call tensor_mpi_reduce(Eocc,  root, tensor_work_comm)
     call tensor_mpi_reduce(Evirt, root, tensor_work_comm)
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     fEc = 0.50E0_tensor_dp*(Eocc + Evirt)

     if( tensor_always_sync ) call tensor_mpi_barrier(tensor_work_comm)

#endif
  end function get_fragment_cc_energy_parallel

  !> \author Patrick Ettenhuber
  !> \date December 2012, modified several times afterwards
  !> \brief calculate aos cc energy in parallel (PDM)
  function get_cc_energy_parallel(t2,gmo,t1) result(Ec)
     implicit none
     !> two electron integrals in the mo-basis
     type(tensor), intent(inout) :: gmo
     !> doubles amplitudes
     type(tensor), intent(in) :: t2
     !> singles amplitudes, optional so that an MP2 contribution can be calculated
     type(tensor), intent(inout),optional :: t1
     !> on return Ec contains the correlation energy
     real(tensor_dp) :: E1,E2,Ec
     real(tensor_dp),pointer :: t2tile(:,:,:,:),gmotile(:,:,:,:),gmotile1d(:)
     real(tensor_dp),pointer :: gmo_tile(:)
     integer :: lt,i,j,a,b,o(t2%mode),da,db,di,dj,gmo_ts
     integer :: order_c(gmo%mode), gmo_ctidx(gmo%mode), gmo_ctdim(gmo%mode), cbuf
     real(tensor_dp), pointer :: gmo_ctile_buf(:,:),gmo_ctile(:,:,:,:)
     integer :: order_e(gmo%mode), gmo_etidx(gmo%mode), gmo_etdim(gmo%mode), ebuf
     integer(kind=tensor_standard_int) :: gmo_ccidx, gmo_ecidx
     real(tensor_dp), pointer :: gmo_etile(:,:,:,:)
#ifdef VAR_PTR_RESHAPE
     real(tensor_dp), contiguous, pointer :: gmo_tile_buf(:,:)
#else
     real(tensor_dp), pointer :: gmo_tile_buf(:,:)
#endif
     integer :: nt,nbuffs,nbuffs_c, nbuffs_e
     integer(kind=tensor_mpi_kind) :: mode
     integer(kind=tensor_long_int) :: tiledim
     real(tensor_dp), external :: ddot
     integer, pointer :: table_iajb(:,:), table_ibja(:,:)
     integer(kind=tensor_mpi_kind), pointer :: reqC(:),reqE(:)
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: me, nnod
#ifdef VAR_MPI
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif

     call tensor_get_rank(me)

     mode = MPI_MODE_NOCHECK

     !reorder gmo to have the same order as t2 - both coulomb and exchange parts
     order_c   = [2,4,1,3]
     order_e   = [4,2,1,3]

     do i=1,gmo%mode
        if(gmo%dims(order_c(i)) /= t2%dims(i))then
           call tensor_status_quit("ERROR(get_cc_energy_parallel):the assumed sorting of the&
              & gmos and t2 amplitudes is incorrect",443)
        endif
     enddo
     !Get the slaves to this routine
     if(me == root)then
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_COMM )
#endif
        if(present(t1))then
           call pdm_tensor_sync(JOB_GET_CC_ENERGY,t2,gmo,t1)
        else
           call pdm_tensor_sync(JOB_GET_MP2_ENERGY,t2,gmo)
        endif
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_WORK )
#endif
     endif

     nbuffs = 6
     if(mod(nbuffs,2)/=0)call tensor_status_quit("ERROR(get_cc_energy_parallel): nbuffs must be an even number",3)
     nbuffs_c = nbuffs/2
     nbuffs_e = nbuffs/2

     call tensor_alloc_mem(gmo_tile_buf,int(gmo%tsize),int(nbuffs))

     E1=0.0E0_tensor_dp
     E2=0.0E0_tensor_dp
     Ec=0.0E0_tensor_dp

     if( alloc_in_dummy )then
        call tensor_lock_wins(gmo,'s',all_nodes = .true.)
        call tensor_alloc_mem(reqC,nbuffs_c)
        call tensor_alloc_mem(reqE,nbuffs_e)
     endif

     !Preload nbuffs_c tiles
     do lt=1,min(nbuffs_c-1,t2%nlti)

        !get offset for global indices and tile indices for gmo
        call get_midx(t2%ti(lt)%gt,o,t2%ntpm,t2%mode)
        do j=1,t2%mode
           gmo_ctidx(order_c(j)) = o(j)
           gmo_etidx(order_e(j)) = o(j)
        enddo

        call get_tile_dim(gmo_ctdim,gmo,gmo_ctidx)
        call get_tile_dim(gmo_etdim,gmo,gmo_etidx)


        gmo_ts = 1
        do j=1,gmo%mode
#ifdef VAR_LSDEBUG
           if(gmo_ctdim(order_c(j)) /= gmo_etdim(order_e(j)))then
              print *,me,j,"wrong",gmo_ctdim(order_c(j)),gmo_etdim(order_e(j)),"tdim",gmo_ctdim,",",gmo_etdim
              call tensor_status_quit("ERROR(get_cc_energy_parallel): something wrong with the gmo tiles",-1)
           endif
#endif
           gmo_ts = gmo_ts * gmo_ctdim(j)
        enddo

        cbuf = mod(lt,nbuffs_c) + 1

        gmo_ccidx = get_cidx(gmo_ctidx,gmo%ntpm,gmo%mode)
        gmo_ecidx = get_cidx(gmo_etidx,gmo%ntpm,gmo%mode)

        !GET COULOMB TILE
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_COMM )
#endif
        if( alloc_in_dummy )then

           call tensor_get_tile(gmo,gmo_ccidx,gmo_tile_buf(:,cbuf),gmo_ts,lock_set=.true.,req=reqC(cbuf))

        else

           call tensor_lock_win(gmo,gmo_ccidx,'s',assert = mode)
           call tensor_get_tile(gmo,gmo_ccidx,gmo_tile_buf(:,cbuf),gmo_ts,lock_set=.true.)

        endif
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_WORK )
#endif

        !GET EXCHANGE TILE
        if(gmo_ccidx/=gmo_ecidx)then
           ebuf = nbuffs_c + mod(lt,nbuffs_c) + 1
#ifdef TENSORS_IN_LSDALTON
           call time_start_phase( PHASE_COMM )
#endif
           if( alloc_in_dummy )then

              call tensor_get_tile(gmo,gmo_ecidx,gmo_tile_buf(:,ebuf),gmo_ts,lock_set=.true.,req=reqE(cbuf))

           else
              call tensor_lock_win(gmo,gmo_ecidx,'s',assert = mode)
              call tensor_get_tile(gmo,gmo_ecidx,gmo_tile_buf(:,ebuf),gmo_ts,lock_set=.true.)

           endif
#ifdef TENSORS_IN_LSDALTON
           call time_start_phase( PHASE_WORK )
#endif
        else
           ebuf = cbuf
        endif


     enddo


     !DO LOOP OVER LOCAL TILES OF T2
     do lt=1,t2%nlti


        !PREFETCH TILES
        nt = lt + nbuffs_c - 1
        if( nt <= t2%nlti )then
           !get offset for global indices and tile indices for gmo
           call get_midx(t2%ti(nt)%gt,o,t2%ntpm,t2%mode)
           do j=1,t2%mode
              gmo_ctidx(order_c(j)) = o(j)
              gmo_etidx(order_e(j)) = o(j)
           enddo

           call get_tile_dim(gmo_ctdim,gmo,gmo_ctidx)
           call get_tile_dim(gmo_etdim,gmo,gmo_etidx)


           gmo_ts = 1
           do j=1,gmo%mode
#ifdef VAR_LSDEBUG
              if(gmo_ctdim(order_c(j)) /= gmo_etdim(order_e(j)))then
                 print *,me,j,"wrong",gmo_ctdim(order_c(j)),gmo_etdim(order_e(j)),"tdim",gmo_ctdim,",",gmo_etdim
                 call tensor_status_quit("ERROR(get_cc_energy_parallel): something wrong with the gmo tiles",-1)
              endif
#endif
              gmo_ts = gmo_ts * gmo_ctdim(j)
           enddo

           cbuf = mod(nt,nbuffs_c) + 1

           gmo_ccidx = get_cidx(gmo_ctidx,gmo%ntpm,gmo%mode)
           gmo_ecidx = get_cidx(gmo_etidx,gmo%ntpm,gmo%mode)

           !GET COULOMB TILE
#ifdef TENSORS_IN_LSDALTON
           call time_start_phase( PHASE_COMM )
#endif
           if( alloc_in_dummy )then

              call tensor_get_tile(gmo,gmo_ccidx,gmo_tile_buf(:,cbuf),gmo_ts,lock_set=.true.,req=reqC(cbuf))

           else
              call tensor_lock_win(gmo,gmo_ccidx,'s',assert = mode)
              call tensor_get_tile(gmo,gmo_ccidx,gmo_tile_buf(:,cbuf),gmo_ts,lock_set=.true.)
           endif
#ifdef TENSORS_IN_LSDALTON
           call time_start_phase( PHASE_WORK )
#endif

           !GET EXCHANGE TILE
           if(gmo_ccidx/=gmo_ecidx)then
              ebuf = nbuffs_c + mod(nt,nbuffs_c) + 1
#ifdef TENSORS_IN_LSDALTON
              call time_start_phase( PHASE_COMM )
#endif
              if( alloc_in_dummy )then

                 call tensor_get_tile(gmo,gmo_ecidx,gmo_tile_buf(:,ebuf),gmo_ts,lock_set=.true.,req=reqE(cbuf))

              else
                 call tensor_lock_win(gmo,gmo_ecidx,'s',assert = mode)
                 call tensor_get_tile(gmo,gmo_ecidx,gmo_tile_buf(:,ebuf),gmo_ts,lock_set=.true.)
              endif
#ifdef TENSORS_IN_LSDALTON
              call time_start_phase( PHASE_WORK )
#endif
           else
              ebuf = cbuf
           endif

        endif

        !get offset for global indices and tile indices for gmo
        call get_midx(t2%ti(lt)%gt,o,t2%ntpm,t2%mode)
        do j=1,t2%mode
           gmo_ctidx(order_c(j)) = o(j)
           gmo_etidx(order_e(j)) = o(j)
           o(j)=(o(j)-1)*t2%tdim(j)
        enddo

        call get_tile_dim(gmo_ctdim,gmo,gmo_ctidx)
        call get_tile_dim(gmo_etdim,gmo,gmo_etidx)


        gmo_ts = 1
        do j=1,gmo%mode
#ifdef VAR_LSDEBUG
           if(gmo_ctdim(order_c(j)) /= gmo_etdim(order_e(j)))then
              print *,me,j,"wrong",gmo_ctdim(order_c(j)),gmo_etdim(order_e(j)),"tdim",gmo_ctdim,",",gmo_etdim
              call tensor_status_quit("ERROR(get_cc_energy_parallel): something wrong with the gmo tiles",-1)
           endif
#endif
           gmo_ts = gmo_ts * gmo_ctdim(j)
        enddo

        cbuf = mod(lt,nbuffs_c) + 1

        gmo_ccidx = get_cidx(gmo_ctidx,gmo%ntpm,gmo%mode)
        gmo_ecidx = get_cidx(gmo_etidx,gmo%ntpm,gmo%mode)

#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_COMM )
#endif
        if( alloc_in_dummy )then
           call tensor_mpi_wait(reqC(cbuf))
        else
           call tensor_unlock_win(gmo,int(gmo_ccidx))
        endif
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_WORK )
#endif

        if(gmo_ccidx/=gmo_ecidx)then
           ebuf = nbuffs_c + mod(lt,nbuffs_c) + 1
#ifdef TENSORS_IN_LSDALTON
           call time_start_phase( PHASE_COMM )
#endif
           if( alloc_in_dummy )then
              call tensor_mpi_wait(reqE(cbuf))
           else
              call tensor_unlock_win(gmo,int(gmo_ecidx))
           endif
#ifdef TENSORS_IN_LSDALTON
           call time_start_phase( PHASE_WORK )
#endif
        else
           ebuf = cbuf
        endif

#if defined(VAR_PTR_RESHAPE) && !defined(VAR_PGF90)
     gmo_ctile(1:gmo_ctdim(1),1:gmo_ctdim(2),1:gmo_ctdim(3),1:gmo_ctdim(4)) => gmo_tile_buf(1:,cbuf)
     gmo_etile(1:gmo_etdim(1),1:gmo_etdim(2),1:gmo_etdim(3),1:gmo_etdim(4)) => gmo_tile_buf(1:,ebuf)
     t2tile(1:t2%ti(lt)%d(1),1:t2%ti(lt)%d(2),1:t2%ti(lt)%d(3),1:t2%ti(lt)%d(4)) => t2%ti(lt)%t
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
     call c_f_pointer(c_loc(gmo_tile_buf(1,cbuf)),gmo_ctile,gmo_ctdim)
     call c_f_pointer(c_loc(gmo_tile_buf(1,ebuf)),gmo_etile,gmo_etdim)
     call c_f_pointer(c_loc(t2%ti(lt)%t(1)),t2tile,t2%ti(lt)%d)
#else
     call tensor_status_quit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

     da = t2%ti(lt)%d(1)
     db = t2%ti(lt)%d(2)
     di = t2%ti(lt)%d(3)
     dj = t2%ti(lt)%d(4)
     !count over local indices
     !$OMP  PARALLEL DO DEFAULT(NONE) SHARED(o,t2tile,gmo_ctile,gmo_etile,&
     !$OMP  da,db,di,dj) PRIVATE(i,j,a,b) REDUCTION(+:E1,E2) COLLAPSE(3)
     do j=1,dj
        do i=1,di
           do b=1,db
              do a=1,da

                 E2 = E2 + t2tile(a,b,i,j)*&
                    & (2.0E0_tensor_dp*  gmo_ctile(i,a,j,b) - gmo_etile(i,b,j,a))
              enddo 
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
     if(present(t1))then
        !$OMP  PARALLEL DO DEFAULT(NONE) SHARED(o,t1,gmo_ctile,gmo_etile,&
        !$OMP  da,db,di,dj) PRIVATE(i,j,a,b) REDUCTION(+:E1,E2) COLLAPSE(3)
        do j=1,dj
           do i=1,di
              do b=1,db
                 do a=1,da

                    E1 = E1 + ( t1%elm2(a+o(1),i+o(3))*t1%elm2(b+o(2),j+o(4)) ) * &
                       (2.0E0_tensor_dp*gmo_ctile(i,a,j,b)-gmo_etile(i,b,j,a))

                 enddo 
              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO
     endif
     t2tile    => null()
     gmo_ctile => null()
     gmo_etile => null()
  enddo

  if( alloc_in_dummy )then
     call tensor_unlock_wins(gmo,all_nodes = .true.)
     call tensor_free_mem(reqC)
     call tensor_free_mem(reqE)
  endif

#ifdef TENSORS_IN_LSDALTON
  call time_start_phase( PHASE_COMM )
#endif
  if(present(t1))call tensor_mpi_reduce(E1,root,tensor_work_comm)
  call tensor_mpi_reduce(E2,root,tensor_work_comm)
#ifdef TENSORS_IN_LSDALTON
  call time_start_phase( PHASE_WORK )
#endif

  if(present(t1))then
     Ec=E1+E2
  else
     Ec=E2
  endif

  call tensor_free_mem(gmo_tile_buf)

  if( tensor_always_sync ) call tensor_mpi_barrier(tensor_work_comm)
#else
  Ec = 0.0E0_tensor_dp
#endif
  end function get_cc_energy_parallel

  subroutine lspdm_extract_eos_indices_virt(Arr,tensor_full,nEOS,EOS_idx)
     implicit none
     !> Array where EOS indices where are extracted
     type(tensor),intent(inout) :: Arr
     !> Original array in the order nv,no,nv,no
     type(tensor),intent(in) :: tensor_full
     !> Number of EOS indices
     integer,intent(in) :: nEOS
     !> List of EOS indices in the total (EOS+buffer) list of orbitals
     integer, dimension(nEOS),intent(in) :: EOS_idx
     integer :: nocc,nvirt,i,a,b,j,ix,jx
     integer, dimension(4) :: new_dims, o
     integer, pointer :: idxatil(:), idxbtil(:)
     integer :: lt, di, da, dj, db, nidxa, nidxb, a_eos,b_eos
     real(tensor_dp), pointer :: tile(:,:,:,:)
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: me, nnod
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif

     call tensor_get_rank(me)
     call tensor_get_size(nnod)

     ! Initialize stuff
     ! ****************
     nocc     = tensor_full%dims(2)  ! Total number of occupied orbitals
     nvirt    = tensor_full%dims(1)  ! Total number of virtual orbitals
     new_dims = [nEOS,nocc,nEOS,nocc] ! nEOS=Number of occupied EOS orbitals
#ifdef VAR_MPI
     if(me == root.and.tensor_full%access_type==AT_MASTER_ACCESS)then

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_COMM )
#endif
     call pdm_tensor_sync(JOB_TENSOR_EXTRACT_VEOS,tensor_full)

     call tensor_buffer(nEOS,root=root,comm=tensor_work_comm)
     call tensor_buffer(EOS_idx,nEOS)
     call tensor_buffer(4)
     call tensor_buffer(new_dims,4,finalize=.true.)
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
  endif

  call tensor_alloc_mem(idxatil,nEOS)
  call tensor_alloc_mem(idxbtil,nEOS)

  do lt=1,tensor_full%nlti

     call get_midx(tensor_full%ti(lt)%gt,o,tensor_full%ntpm,tensor_full%mode)

     !#ifdef VAR_PTR_RESHAPE
     !        tile(1:tensor_full%ti(lt)%d(1),1:tensor_full%ti(lt)%d(2),&
     !        &1:tensor_full%ti(lt)%d(3),1:tensor_full%ti(lt)%d(4)) => tensor_full%ti(lt)%t
     !#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
     call c_f_pointer(c_loc(tensor_full%ti(lt)%t(1)),tile,tensor_full%ti(lt)%d)
     !#else
     !        call tensor_status_quit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
     !#endif

     !get offset for tile counting
     do j=1,tensor_full%mode
        o(j)=(o(j)-1)*tensor_full%tdim(j)
     enddo

     da = tensor_full%ti(lt)%d(1)
     di = tensor_full%ti(lt)%d(2)
     db = tensor_full%ti(lt)%d(3)
     dj = tensor_full%ti(lt)%d(4)

     ! GET EOS mapping to the tile
     nidxa = 0
     do a_eos = 1, nEOS
        do a = 1, da
           if( o(1) + a == EOS_idx(a_eos))then
              idxatil(nidxa+1) = a
              nidxa = nidxa + 1
           endif
        enddo
     enddo

     nidxb = 0
     do b_eos = 1, nEOS
        do b = 1, db
           if( o(3) + b == EOS_idx(b_eos))then
              idxbtil(nidxb+1) = b
              nidxb = nidxb + 1
           endif
        enddo
     enddo

     if(nidxa > 0 .and. nidxb>0)then

        do j=1,dj
           do b=1,nidxb
              do i=1,di
                 do a=1,nidxa
                    Arr%elm4(o(1)+idxatil(a),o(2)+i,o(3)+idxbtil(b),o(4)+j) = tile(idxatil(a),i,idxbtil(b),j)
                 end do
              end do
           end do
        end do

     endif

     tile => null()
  enddo

  call tensor_free_mem(idxatil)
  call tensor_free_mem(idxbtil)

#ifdef TENSORS_IN_LSDALTON
  call time_start_phase( PHASE_COMM )
#endif
  if(tensor_full%access_type==AT_MASTER_ACCESS)then
     call tensor_mpi_reduce(Arr%elm1,Arr%nelms,root,tensor_work_comm)
  else if(tensor_full%access_type==AT_ALL_ACCESS)then
     call tensor_mpi_allreduce(Arr%elm1,Arr%nelms,tensor_work_comm)
  endif
#ifdef TENSORS_IN_LSDALTON
  call time_start_phase( PHASE_WORK )
#endif

  if( tensor_always_sync ) call tensor_mpi_barrier(tensor_work_comm)
#endif

  end subroutine lspdm_extract_eos_indices_virt

  subroutine lspdm_extract_eos_indices_occ(Arr,tensor_full,nEOS,EOS_idx)
     implicit none
     !> Array where EOS indices where are extracted
     type(tensor),intent(inout) :: Arr
     !> Original array in the order nv,no,nv,no
     type(tensor),intent(in) :: tensor_full
     !> Number of EOS indices
     integer,intent(in) :: nEOS
     !> List of EOS indices in the total (EOS+buffer) list of orbitals
     integer, dimension(nEOS),intent(in) :: EOS_idx
     integer :: nocc,nvirt,i,a,b,j,ix,jx
     integer, dimension(4) :: new_dims, o
     integer, pointer :: idxitil(:), idxjtil(:)
     integer :: lt, di, da, dj, db, nidxi, nidxj, i_eos, j_eos
     real(tensor_dp), pointer :: tile(:,:,:,:)
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: me, nnod
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif

     call tensor_get_rank(me)

     ! Initialize stuff
     ! ****************
     nocc     = tensor_full%dims(2)  ! Total number of occupied orbitals
     nvirt    = tensor_full%dims(1)  ! Total number of virtual orbitals
     new_dims = [nvirt,nEOS,nvirt,nEOS] ! nEOS=Number of occupied EOS orbitals

#ifdef VAR_MPI
     if(me == root.and.tensor_full%access_type==AT_MASTER_ACCESS)then
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_COMM )
#endif
     call pdm_tensor_sync(JOB_TENSOR_EXTRACT_OEOS,tensor_full)

     call tensor_buffer(nEOS,root=root,comm=tensor_work_comm)
     call tensor_buffer(EOS_idx,nEOS)
     call tensor_buffer(4)
     call tensor_buffer(new_dims,4,finalize=.true.)
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
  endif

  call tensor_alloc_mem(idxitil,nEOS)
  call tensor_alloc_mem(idxjtil,nEOS)

  do lt=1,tensor_full%nlti

     call get_midx(tensor_full%ti(lt)%gt,o,tensor_full%ntpm,tensor_full%mode)

     !#ifdef VAR_PTR_RESHAPE
     !        tile(1:tensor_full%ti(lt)%d(1),1:tensor_full%ti(lt)%d(2),&
     !        &1:tensor_full%ti(lt)%d(3),1:tensor_full%ti(lt)%d(4)) => tensor_full%ti(lt)%t
     !#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
     call c_f_pointer(c_loc(tensor_full%ti(lt)%t(1)),tile,tensor_full%ti(lt)%d)
     !#else
     !        call tensor_status_quit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
     !#endif

     !get offset for tile counting
     do j=1,tensor_full%mode
        o(j)=(o(j)-1)*tensor_full%tdim(j)
     enddo

     da = tensor_full%ti(lt)%d(1)
     di = tensor_full%ti(lt)%d(2)
     db = tensor_full%ti(lt)%d(3)
     dj = tensor_full%ti(lt)%d(4)

     ! GET EOS mapping to the tile
     nidxi = 0
     do i_eos = 1, nEOS
        do i = 1, di
           if( o(2) + i == EOS_idx(i_eos))then
              idxitil(nidxi+1) = i
              nidxi = nidxi + 1
           endif
        enddo
     enddo

     nidxj = 0
     do j_eos = 1, nEOS
        do j = 1, dj
           if( o(4) + j == EOS_idx(j_eos))then
              idxjtil(nidxj+1) = j
              nidxj = nidxj + 1
           endif
        enddo
     enddo

     if(nidxi > 0 .and. nidxj>0)then

        do j=1,nidxj
           do b=1,db
              do i=1,nidxi
                 do a=1,da
                    Arr%elm4(o(1)+a,o(2)+idxitil(i),o(3)+b,o(4)+idxjtil(j)) = tile(a,idxitil(i),b,idxjtil(j))
                 end do
              end do
           end do
        end do

     endif

     tile => null()
  enddo

  call tensor_free_mem(idxitil)
  call tensor_free_mem(idxjtil)

#ifdef TENSORS_IN_LSDALTON
  call time_start_phase( PHASE_COMM )
#endif
  if(tensor_full%access_type==AT_MASTER_ACCESS)then
     call tensor_mpi_reduce(Arr%elm1,Arr%nelms,root,tensor_work_comm)
  else if(tensor_full%access_type==AT_ALL_ACCESS)then
     call tensor_mpi_allreduce(Arr%elm1,Arr%nelms,tensor_work_comm)
  endif
#ifdef TENSORS_IN_LSDALTON
  call time_start_phase( PHASE_WORK )
#endif

  if( tensor_always_sync ) call tensor_mpi_barrier(tensor_work_comm)
#endif

  end subroutine lspdm_extract_eos_indices_occ

  subroutine lspdm_extract_decnp_indices_virt(Arr,tensor_full,nEOS,EOS_idx)
     implicit none
     !> Array where EOS indices where are extracted
     type(tensor),intent(inout) :: Arr
     !> Original array in the order nv,no,nv,no
     type(tensor),intent(in) :: tensor_full
     !> Number of EOS indices
     integer,intent(in) :: nEOS
     !> List of EOS indices in the total (EOS+buffer) list of orbitals
     integer, dimension(nEOS),intent(in) :: EOS_idx
     integer :: nocc,nvirt,i,a,b,j
     integer, dimension(4) :: new_dims, o
     integer, pointer :: idxatil(:)
     integer :: lt, di, da, dj, db, nidxa, a_eos
     real(tensor_dp), pointer :: tile(:,:,:,:)
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: me, nnod
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif

     call tensor_get_rank(me)

     ! Initialize stuff
     ! ****************
     nocc     = tensor_full%dims(2)    ! Total number of occupied orbitals
     nvirt    = tensor_full%dims(1)    ! Total number of virtual orbitals
     new_dims = [nEOS,nocc,nvirt,nocc] ! nEOS=Number of occupied EOS orbitals

#ifdef VAR_MPI
     if(me == root.and. tensor_full%access_type==AT_MASTER_ACCESS)then

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_COMM )
#endif
     call pdm_tensor_sync(JOB_TENSOR_EXTRACT_VDECNP,tensor_full)

     call tensor_buffer(nEOS,root=root,comm=tensor_work_comm)
     call tensor_buffer(EOS_idx,nEOS)
     call tensor_buffer(4)
     call tensor_buffer(new_dims,4,finalize=.true.)
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
  endif

  call tensor_alloc_mem(idxatil,nEOS)

  do lt=1,tensor_full%nlti

     call get_midx(tensor_full%ti(lt)%gt,o,tensor_full%ntpm,tensor_full%mode)

     !#ifdef VAR_PTR_RESHAPE
     !        tile(1:tensor_full%ti(lt)%d(1),1:tensor_full%ti(lt)%d(2),&
     !        &1:tensor_full%ti(lt)%d(3),1:tensor_full%ti(lt)%d(4)) => tensor_full%ti(lt)%t
     !#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
     call c_f_pointer(c_loc(tensor_full%ti(lt)%t(1)),tile,tensor_full%ti(lt)%d)
     !#else
     !        call tensor_status_quit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
     !#endif

     !get offset for tile counting
     do j=1,tensor_full%mode
        o(j)=(o(j)-1)*tensor_full%tdim(j)
     enddo

     da = tensor_full%ti(lt)%d(1)
     di = tensor_full%ti(lt)%d(2)
     db = tensor_full%ti(lt)%d(3)
     dj = tensor_full%ti(lt)%d(4)

     ! GET EOS mapping to the tile
     nidxa = 0
     do a_eos = 1, nEOS
        do a = 1, da
           if( o(1) + a == EOS_idx(a_eos))then
              idxatil(nidxa+1) = a
              nidxa = nidxa + 1
           endif
        enddo
     enddo

     if(nidxa > 0)then

        do j=1,dj
           do b=1,db
              do i=1,di
                 do a=1,nidxa
                    Arr%elm4(o(1)+idxatil(a),o(2)+i,o(3)+b,o(4)+j) = tile(idxatil(a),i,b,j)
                 end do
              end do
           end do
        end do

     endif

     tile => null()
  enddo

  call tensor_free_mem(idxatil)

#ifdef TENSORS_IN_LSDALTON
  call time_start_phase( PHASE_COMM )
#endif
  if(tensor_full%access_type==AT_MASTER_ACCESS)then
     call tensor_mpi_reduce(Arr%elm1,Arr%nelms,root,tensor_work_comm)
  else if(tensor_full%access_type==AT_ALL_ACCESS)then
     call tensor_mpi_allreduce(Arr%elm1,Arr%nelms,tensor_work_comm)
  endif
#ifdef TENSORS_IN_LSDALTON
  call time_start_phase( PHASE_WORK )
#endif

  if( tensor_always_sync ) call tensor_mpi_barrier(tensor_work_comm)
#endif

  end subroutine lspdm_extract_decnp_indices_virt


  subroutine lspdm_extract_decnp_indices_occ(Arr,tensor_full,nEOS,EOS_idx)
     implicit none
     !> Array where EOS indices where are extracted
     type(tensor),intent(inout) :: Arr
     !> Original array in the order nv,no,nv,no
     type(tensor),intent(in) :: tensor_full
     !> Number of EOS indices
     integer,intent(in) :: nEOS
     !> List of EOS indices in the total (EOS+buffer) list of orbitals
     integer, dimension(nEOS),intent(in) :: EOS_idx
     integer :: nocc,nvirt,i,a,b,j
     integer, dimension(4) :: new_dims, o
     integer, pointer :: idxitil(:)
     integer :: lt, di, da, dj, db, nidxi, i_eos
     real(tensor_dp), pointer :: tile(:,:,:,:)
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: me, nnod
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif

     call tensor_get_rank(me)

     ! Initialize stuff
     ! ****************
     nocc     = tensor_full%dims(2)     ! Total number of occupied orbitals
     nvirt    = tensor_full%dims(1)     ! Total number of virtual orbitals
     new_dims = [nvirt,nEOS,nvirt,nocc] ! nEOS=Number of occupied EOS orbitals

#ifdef VAR_MPI
     if(me == root.and.tensor_full%access_type==AT_MASTER_ACCESS)then

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_COMM )
#endif
     call pdm_tensor_sync(JOB_TENSOR_EXTRACT_ODECNP,tensor_full)

     call tensor_buffer(nEOS,root=root,comm=tensor_work_comm)
     call tensor_buffer(EOS_idx,nEOS)
     call tensor_buffer(4)
     call tensor_buffer(new_dims,4,finalize=.true.)
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
  endif

  call tensor_alloc_mem(idxitil,nEOS)

  do lt=1,tensor_full%nlti

     call get_midx(tensor_full%ti(lt)%gt,o,tensor_full%ntpm,tensor_full%mode)

     !#ifdef VAR_PTR_RESHAPE
     !        tile(1:tensor_full%ti(lt)%d(1),1:tensor_full%ti(lt)%d(2),&
     !        &1:tensor_full%ti(lt)%d(3),1:tensor_full%ti(lt)%d(4)) => tensor_full%ti(lt)%t
     !#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
     call c_f_pointer(c_loc(tensor_full%ti(lt)%t(1)),tile,tensor_full%ti(lt)%d)
     !#else
     !        call tensor_status_quit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
     !#endif

     !get offset for tile counting
     do j=1,tensor_full%mode
        o(j)=(o(j)-1)*tensor_full%tdim(j)
     enddo

     da = tensor_full%ti(lt)%d(1)
     di = tensor_full%ti(lt)%d(2)
     db = tensor_full%ti(lt)%d(3)
     dj = tensor_full%ti(lt)%d(4)

     ! GET EOS mapping to the tile
     nidxi = 0
     do i_eos = 1, nEOS
        do i = 1, di
           if( o(2) + i == EOS_idx(i_eos))then
              idxitil(nidxi+1) = i
              nidxi = nidxi + 1
           endif
        enddo
     enddo

     if(nidxi > 0)then

        do j=1,dj
           do b=1,db
              do i=1,nidxi
                 do a=1,da
                    Arr%elm4(o(1)+a,o(2)+idxitil(i),o(3)+b,o(4)+j) = tile(a,idxitil(i),b,j)
                 end do
              end do
           end do
        end do

     endif

     tile => null()
  enddo

  call tensor_free_mem(idxitil)

#ifdef TENSORS_IN_LSDALTON
  call time_start_phase( PHASE_COMM )
#endif
  if(tensor_full%access_type==AT_MASTER_ACCESS)then
     call tensor_mpi_reduce(Arr%elm1,Arr%nelms,root,tensor_work_comm)
  else if(tensor_full%access_type==AT_ALL_ACCESS)then
     call tensor_mpi_allreduce(Arr%elm1,Arr%nelms,tensor_work_comm)
  endif
#ifdef TENSORS_IN_LSDALTON
  call time_start_phase( PHASE_WORK )
#endif

  if( tensor_always_sync ) call tensor_mpi_barrier(tensor_work_comm)
#endif

  end subroutine lspdm_extract_decnp_indices_occ

  subroutine lspdm_get_combined_SingleDouble_amplitudes(t1,t2,u)
     implicit none
     !> Singles amplitudes t1(a,i)
     type(tensor),intent(in) :: t1
     !> Doubles amplitudes t2(a,i,b,j)
     type(tensor),intent(in) :: t2
     !> Combined single+double amplitudes
     type(tensor),intent(inout) :: u
     integer :: i,j,a,b,nocc,nvirt,da,db,di,dj,gtnr,lt,nelt
     integer :: o(4)
     real(tensor_dp), pointer :: tt(:,:,:,:), ut(:,:,:,:)
     real(tensor_dp), pointer :: ttile(:)
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: me, nnod

     call tensor_get_rank(me)

#ifdef VAR_MPI
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_COMM )
#endif
     if( t2%access_type == AT_MASTER_ACCESS .and. me == root)then
        if(t1%itype == TT_REPLICATED.or.t1%itype==TT_TILED_DIST)then
           call pdm_tensor_sync(JOB_GET_COMBINEDT1T2_1,t1,t2,u)
        else if(t1%itype == TT_DENSE)then
           call pdm_tensor_sync(JOB_GET_COMBINEDT1T2_2,t2,u)
           call tensor_buffer(t1%dims,2,root=root,comm=tensor_work_comm)
           call tensor_buffer(t1%elm1,t1%nelms,finalize=.true.)
        else
           call tensor_status_quit("ERROR(lspdm_get_combined_SingleDouble_amplitudes):no valid t1%itype",-1)
        endif
     endif
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif

     select case(t1%itype)
     case(TT_DENSE,TT_REPLICATED)

        call tensor_alloc_mem(ttile,u%tsize)

        do lt = 1, u%nlti

           gtnr = u%ti(lt)%gt

           nelt = u%ti(lt)%e

           call get_midx(gtnr,o,u%ntpm,u%mode)

           !This is just a copy since we enforced u to have the same
           !distribution as t, but for the sake of generality we use the MPI_GET
           !to perform the copy. For now, prefetching is not necessary
#ifdef TENSORS_IN_LSDALTON
           call time_start_phase( PHASE_COMM )
#endif
           call tensor_get_tile(t2,int(gtnr,kind=tensor_standard_int),ttile,nelt)
#ifdef TENSORS_IN_LSDALTON
           call time_start_phase( PHASE_WORK )
#endif

           !Facilitate access
#ifdef VAR_PTR_RESHAPE
           ut(1:u%ti(lt)%d(1),1:u%ti(lt)%d(2),1:u%ti(lt)%d(3),1:u%ti(lt)%d(4)) => u%ti(lt)%t
           tt(1:u%ti(lt)%d(1),1:u%ti(lt)%d(2),1:u%ti(lt)%d(3),1:u%ti(lt)%d(4)) => ttile
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
           call c_f_pointer( c_loc(u%ti(lt)%t(1)), ut, u%ti(lt)%d )
           call c_f_pointer( c_loc(ttile(1)),      tt, u%ti(lt)%d )
#else
           call tensor_status_quit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

           !get offset for tile counting
           do j=1,u%mode
              o(j)=(o(j)-1)*u%tdim(j)
           enddo

           da = u%ti(lt)%d(1)
           di = u%ti(lt)%d(2)
           db = u%ti(lt)%d(3)
           dj = u%ti(lt)%d(4)

           do j=1,dj
              do b=1,db
                 do i=1,di
                    do a=1,da
                       ut(a,i,b,j) = tt(a,i,b,j) + t1%elm2(o(1)+a,o(2)+i) * t1%elm2(o(3)+b,o(4)+j)
                    end do
                 end do
              end do
           end do


           ut => null()
           tt => null()
        enddo

        call tensor_free_mem(ttile)

     case default
        call tensor_status_quit("ERROR(lspdm_get_combined_SingleDouble_amplitudes): not yet&
           & implemented for tiled t1, but it should be simple :)",-1)
     end select

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_IDLE )
#endif
     call tensor_mpi_barrier(tensor_work_comm)
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
#endif

  end subroutine lspdm_get_combined_SingleDouble_amplitudes


  subroutine get_info_for_mpi_get_and_reorder_t1(arr,table_iajb,table_ibja, &
        & dims,ord,t1,t1tile)
     implicit none

     type(tensor), intent(in) :: arr, t1
     integer, intent(inout) :: table_iajb(:,:), table_ibja(:,:)   
     integer, intent(in) :: dims(4), ord(4)
     real(tensor_dp), intent(inout) :: t1tile(:)

     !> mode and combined idices of the tile:
     integer :: timode(4), ticomb
     !> mode index of the elmt in the tile:
     integer :: modeinti(4)
     !> mode idex of the elmt in dense array
     integer ::  modeinde(4)

     integer :: source, pos, k, idx, da, db, di, dj, a, b, i, j, dpos, didx, dwidx
#ifdef VAR_MPI
     da = dims(1)
     db = dims(2)
     di = dims(3)
     dj = dims(4)

     ! Make table of indices:
     !$OMP  PARALLEL DO DEFAULT(NONE) SHARED(arr,ord,t1,t1tile, &
     !$OMP  table_iajb,table_ibja,da,db,di,dj) PRIVATE(i,j,k,a,b,idx,pos, &
     !$OMP  source,modeinde,modeinti,timode,ticomb,dpos,didx,dwidx) COLLAPSE(3)
     do j=1,dj
        do i=1,di
           do b=1,db
              do a=1,da
                 ! get combined index for tiles:
                 idx = a + (b-1)*da + (i-1)*da*db + (j-1)*da*db*di

                 ! GET G_IAJB ELMT FROM PDM ARRAY:    
                 modeinde = [i+ord(3),a+ord(1),j+ord(4),b+ord(2)]
                 ! get tile mode index:
                 do k=1,arr%mode
                    timode(k)   = (modeinde(k)-1)/arr%tdim(k) + 1
                    modeinti(k) = mod((modeinde(k)-1) , arr%tdim(k)) + 1
                 end do
                 ticomb = get_cidx(timode,arr%ntpm,arr%mode)
                 pos    = get_cidx(modeinti,arr%tdim,arr%mode)
                 call get_residence_of_tile(arr, ticomb, source, dpos, didx, dwidx)

                 ! get tile elmts from source:
                 table_iajb(idx,:) = [pos,source,ticomb]

                 ! GET G_IBJA ELMT FROM PDM ARRAY:    
                 modeinde = [i+ord(3),b+ord(2),j+ord(4),a+ord(1)]
                 ! get tile mode index:
                 do k=1,arr%mode
                    timode(k)   = (modeinde(k)-1)/arr%tdim(k) + 1
                    modeinti(k) = mod((modeinde(k)-1) , arr%tdim(k)) + 1
                 end do
                 ticomb = get_cidx(timode,arr%ntpm,arr%mode)
                 pos    = get_cidx(modeinti,arr%tdim,arr%mode)
                 call get_residence_of_tile(arr, ticomb, source, dpos, didx, dwidx)

                 ! get tile elmts from source:
                 table_ibja(idx,:) = [pos,source,ticomb]

                 ! reorder t1 contributions into one big tile:
                 t1tile(idx) = t1%elm2(a+ord(1),i+ord(3))*t1%elm2(b+ord(2),j+ord(4))
              enddo 
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
#endif
  end subroutine get_info_for_mpi_get_and_reorder_t1


  function get_rpa_energy_parallel(t2,gmo) result(Ec)
     implicit none
     !> two electron integrals in the mo-basis
     type(tensor), intent(inout) :: gmo
     !> doubles amplitudes
     type(tensor), intent(in) :: t2
     !> on return Ec contains the correlation energy
     real(tensor_dp) :: E1,E2,Ec
     real(tensor_dp),pointer :: t(:,:,:,:)
     integer :: lt,i,j,a,b,o(t2%mode),da,db,di,dj
     logical :: use_bg
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: me, nnod

     call tensor_get_rank(me)

#ifdef VAR_MPI
     !Get the slaves to this routine
     if(me==root)then
        call pdm_tensor_sync(JOB_GET_RPA_ENERGY,t2,gmo)
     endif

     use_bg = (tensor_bg_init()) .and. (tensor_bg_free() > gmo%nelms)

     call memory_allocate_tensor_dense(gmo, use_bg)
     call cp_tileddata2fort(gmo,gmo%elm1,gmo%nelms,.true.)

     E2=0.0E0_tensor_dp
     Ec=0.0E0_tensor_dp
     do lt=1,t2%nlti

#ifdef VAR_PTR_RESHAPE
        t(1:t2%ti(lt)%d(1),1:t2%ti(lt)%d(2),1:t2%ti(lt)%d(3),1:t2%ti(lt)%d(4)) => t2%ti(lt)%t
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
        call c_f_pointer(c_loc(t2%ti(lt)%t(1)),t,t2%ti(lt)%d)
#else
        call tensor_status_quit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

        !get offset for global indices
        call get_midx(t2%ti(lt)%gt,o,t2%ntpm,t2%mode)
        do j=1,t2%mode
           o(j)=(o(j)-1)*t2%tdim(j)
        enddo

        da = t2%ti(lt)%d(1)
        db = t2%ti(lt)%d(2)
        di = t2%ti(lt)%d(3)
        dj = t2%ti(lt)%d(4)
        !count over local indices
        !count over local indices
        !$OMP  PARALLEL DO DEFAULT(NONE) SHARED(gmo,o,t,&
        !$OMP  da,db,di,dj) PRIVATE(i,j,a,b) REDUCTION(+:E2) COLLAPSE(3)
        do j=1,dj
           do i=1,di
              do b=1,db
                 do a=1,da

                    E2 = E2 + t(a,b,i,j)*&
                       & (1.0E0_tensor_dp*  gmo%elm4(i+o(3),a+o(1),j+o(4),b+o(2)))

                 enddo 
              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO
        nullify(t)
     enddo

     call tensor_deallocate_dense(gmo)

     call tensor_mpi_reduce(E2,root,tensor_work_comm)

     Ec = E2
     if( tensor_always_sync ) call tensor_mpi_barrier(tensor_work_comm)
#else
     Ec = 0.0E0_tensor_dp
#endif
  end function get_rpa_energy_parallel


  function get_sosex_cont_parallel(t2,gmo) result(Ec)
     implicit none
     !> two electron integrals in the mo-basis
     type(tensor), intent(inout) :: gmo
     !> doubles amplitudes
     type(tensor), intent(in) :: t2
     !> on return Ec contains the correlation energy
     real(tensor_dp) :: E2,Ec
     real(tensor_dp),pointer :: t(:,:,:,:)
     integer :: lt,i,j,a,b,o(t2%mode),da,db,di,dj
     logical :: use_bg
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: me, nnod

     call tensor_get_rank(me)

#ifdef VAR_MPI
     !Get the slaves to this routine
     if(me==root)then
        call pdm_tensor_sync(JOB_GET_SOS_ENERGY,t2,gmo)
     endif

     use_bg = (tensor_bg_init()) .and. (tensor_bg_free() > gmo%nelms)

     call memory_allocate_tensor_dense(gmo,use_bg)
     call cp_tileddata2fort(gmo,gmo%elm1,gmo%nelms,.true.)

     E2=0.0E0_tensor_dp
     Ec=0.0E0_tensor_dp
     do lt=1,t2%nlti

#ifdef VAR_PTR_RESHAPE
        t(1:t2%ti(lt)%d(1),1:t2%ti(lt)%d(2),1:t2%ti(lt)%d(3),1:t2%ti(lt)%d(4)) => t2%ti(lt)%t
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
        call c_f_pointer(c_loc(t2%ti(lt)%t(1)),t,t2%ti(lt)%d)
#else
        call tensor_status_quit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

        !get offset for global indices
        call get_midx(t2%ti(lt)%gt,o,t2%ntpm,t2%mode)
        do j=1,t2%mode
           o(j)=(o(j)-1)*t2%tdim(j)
        enddo

        da = t2%ti(lt)%d(1)
        db = t2%ti(lt)%d(2)
        di = t2%ti(lt)%d(3)
        dj = t2%ti(lt)%d(4)
        !count over local indices
        !$OMP  PARALLEL DO DEFAULT(NONE) SHARED(gmo,o,t,&
        !$OMP  da,db,di,dj) PRIVATE(i,j,a,b) REDUCTION(+:E2) COLLAPSE(3)
        do j=1,dj
           do i=1,di
              do b=1,db
                 do a=1,da

                    E2 = E2 + t(a,b,i,j)*&
                       & (-0.5E0_tensor_dp* gmo%elm4(i+o(3),b+o(2),j+o(4),a+o(1)) )

                 enddo 
              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO
        nullify(t)
     enddo

     call tensor_deallocate_dense(gmo)

     call tensor_mpi_reduce(E2,root,tensor_work_comm)

     Ec = E2
     if( tensor_always_sync ) call tensor_mpi_barrier(tensor_work_comm)
#else
     Ec = 0.0E0_tensor_dp
#endif
  end function get_sosex_cont_parallel

  subroutine lspdm_get_starting_guess(iajb,t2,oof,vvf,spec,prec) 
     implicit none
     type(tensor), intent(inout) :: iajb,t2,oof,vvf
     real(tensor_dp), pointer :: buf_c(:),buf_e(:),t(:,:,:,:),c(:,:,:,:), e(:,:,:,:)
     integer, intent(in) :: spec
     logical, intent(in) :: prec
     integer :: gtidx,lt,o(t2%mode),da,db,di,dj,a,b,i,j,nelms,bcast_spec
     logical :: bcast_prec
     integer :: order_c(4), order_e(4), gmo_ctidx(4), gmo_etidx(4)
     integer(kind=tensor_standard_int) :: gmo_ccidx, gmo_ecidx
     integer, parameter :: MP2AMP       = 1
     integer, parameter :: CCSD_LAG_RHS = 2
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: me, nnod
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)
#ifdef VAR_MPI

     order_c = [3,1,4,2]
     order_e = [3,2,4,1]

     !sanity checks
     if(t2%access_type /= iajb%access_type)then
        call tensor_status_quit("ERROR(lspdm_get_mp2_starting_guess): pdm tensors should have the same access types",-1)
     endif
     if(oof%itype /= vvf%itype)then
        call tensor_status_quit("ERROR(lspdm_get_mp2_starting_guess): Fock matrices should have the same types",-1)
     endif

     do i=1,t2%mode
        if( iajb%dims(i) /= t2%dims(order_c(i)))then
           call tensor_status_quit("ERROR(lspdm_get_mp2_starting_guess): dimensions of iajb and t2 not in the assumed order",-1)
        endif
        if( iajb%tdim(i) /= t2%tdim(order_c(i)))then
           call tensor_status_quit("ERROR(lspdm_get_mp2_starting_guess): tiling of iajb and t2 not as expected",-1)
        endif
     enddo

     if(t2%access_type == AT_MASTER_ACCESS .and. me == root)then
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase(PHASE_COMM)
#endif
        call pdm_tensor_sync(JOB_GET_MP2_ST_GUESS,iajb,t2,oof,vvf)
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase(PHASE_WORK)
#endif
        bcast_spec = spec
        bcast_prec = prec
        call tensor_mpi_bcast(bcast_spec,root,tensor_work_comm)
        call tensor_mpi_bcast(bcast_prec,root,tensor_work_comm)
     endif

     if( oof%itype == TT_DENSE .or. oof%itype == TT_REPLICATED )then

        !TODO: introduce prefetching of tiles and adapt to alloc_in_dummy
        call tensor_alloc_mem(buf_c,iajb%tsize)
        if( spec == CCSD_LAG_RHS )then
           call tensor_alloc_mem(buf_e,iajb%tsize)
        endif

        do lt=1,t2%nlti

           gtidx = t2%ti(lt)%gt
           nelms = t2%ti(lt)%e

           da = t2%ti(lt)%d(1)
           db = t2%ti(lt)%d(2)
           di = t2%ti(lt)%d(3)
           dj = t2%ti(lt)%d(4)

           !get offset for global indices
           call get_midx(gtidx,o,t2%ntpm,t2%mode)

           !get tile numbers for e and c tiles
           do j=1,t2%mode
              gmo_ctidx(j) = o(order_c(j))
              gmo_etidx(j) = o(order_e(j))
           enddo
           gmo_ccidx = get_cidx(gmo_ctidx,iajb%ntpm,iajb%mode)
           gmo_ecidx = get_cidx(gmo_etidx,iajb%ntpm,iajb%mode)

#ifdef TENSORS_IN_LSDALTON
           call time_start_phase(PHASE_COMM)
#endif
           call tensor_get_tile(iajb,gmo_ccidx,buf_c,nelms)
           if( spec == CCSD_LAG_RHS .and. gmo_ccidx/=gmo_ecidx)then
              call tensor_get_tile(iajb,gmo_ecidx,buf_e,nelms)
           endif
#ifdef TENSORS_IN_LSDALTON
           call time_start_phase(PHASE_WORK)
#endif

#ifdef VAR_PTR_RESHAPE
           t(1:t2%ti(lt)%d(1),1:t2%ti(lt)%d(2),1:t2%ti(lt)%d(3),1:t2%ti(lt)%d(4)) => t2%ti(lt)%t
           c(1:di,1:da,1:dj,1:db) => buf_c
           if( spec == CCSD_LAG_RHS )then
              if(gmo_ccidx/=gmo_ecidx)then
                 e(1:di,1:db,1:dj,1:da) => buf_e
              else
                 e(1:di,1:db,1:dj,1:da) => buf_c
              endif
           endif
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
           call c_f_pointer( c_loc(t2%ti(lt)%t(1)), t, t2%ti(lt)%d )
           call c_f_pointer( c_loc(buf_c(1)), c, [di,da,dj,db] )
           if( spec == CCSD_LAG_RHS )then
              if(gmo_ccidx/=gmo_ecidx)then
                 call c_f_pointer( c_loc(buf_e(1)), e, [di,db,dj,da] )
              else
                 call c_f_pointer( c_loc(buf_c(1)), e, [di,db,dj,da] )
              endif
           endif
#else
           call tensor_status_quit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

           do j=1,t2%mode
              o(j)=(o(j)-1)*t2%tdim(j)
           enddo

           !count over local indices
           select case(spec)
           case(MP2AMP)
              !$OMP  PARALLEL DO DEFAULT(NONE) SHARED(c,o,t,oof,vvf,&
              !$OMP  da,db,di,dj) PRIVATE(i,j,a,b) COLLAPSE(3)
              do j=1,dj
                 do i=1,di
                    do b=1,db
                       do a=1,da

                          t(a,b,i,j) = c(i,a,j,b) 

                       enddo 
                    enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
           case(CCSD_LAG_RHS)
              !$OMP  PARALLEL DO DEFAULT(NONE) SHARED(c,e,o,t,oof,vvf,&
              !$OMP  da,db,di,dj) PRIVATE(i,j,a,b) COLLAPSE(3)
              do j=1,dj
                 do i=1,di
                    do b=1,db
                       do a=1,da

                          t(a,b,i,j) = -4.0E0_tensor_dp * c(i,a,j,b) + 2.0E0_tensor_dp * e(i,b,j,a) 

                       enddo 
                    enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
           case default
              call tensor_status_quit("ERROR(lspdm_get_starting_guess): wrong coice of spec",-1)
           end select

           if(prec)then
              !count over local indices
              !$OMP  PARALLEL DO DEFAULT(NONE) SHARED(o,t,oof,vvf,&
              !$OMP  da,db,di,dj) PRIVATE(i,j,a,b) COLLAPSE(3)
              do j=1,dj
                 do i=1,di
                    do b=1,db
                       do a=1,da

                          t(a,b,i,j) = t(a,b,i,j) / &
                             & (oof%elm2(o(3)+i,o(3)+i) - vvf%elm2(o(1)+a,o(1)+a) &
                             &+ oof%elm2(o(4)+j,o(4)+j) - vvf%elm2(o(2)+b,o(2)+b))

                       enddo 
                    enddo
                 enddo
              enddo
              !$OMP END PARALLEL DO
           endif

           t => null()
           c => null()
           e => null()

        enddo

        call tensor_free_mem(buf_c)
        if( spec == CCSD_LAG_RHS )then
           call tensor_free_mem(buf_e)
        endif

     else
        call tensor_status_quit("ERROR(lspdm_get_mp2_starting_guess): the routine does not accept this type of fock matrix",-1)
     end if

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_IDLE)
#endif
     call tensor_mpi_barrier(tensor_work_comm)
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif
#endif
  end subroutine lspdm_get_starting_guess



  !> \brief doubles preconditionning routine for pdm distributed doubles
  !amplitudes
  !> \author Patrick Ettenhuber
  !> \date december 2012
  subroutine precondition_doubles_parallel(omega2,ppfock,qqfock,prec)
     implicit none
     !> doubles residual, occupied and virtual blocks of the fock matrix
     type(tensor), intent(in) :: omega2,ppfock,qqfock
     !> output is the preconditioned doubles residual
     type(tensor), intent(inout) :: prec
     integer :: lt,a, b, i, j, dims(4)
     real(tensor_dp),pointer :: om(:,:,:,:),pp(:,:),qq(:,:),p(:,:,:,:)
     real(tensor_dp) :: nrm
     integer :: t(4),da,db,di,dj
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: me, nnod
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)

#ifdef VAR_MPI

     !CHECK if the distributions are the same, if it becomes necessary, that they
     !are not, then this routine has to be rewritten
     if(omega2%tdim(1)/=prec%tdim(1).or.omega2%tdim(2)/=prec%tdim(2).or.&
        &omega2%tdim(3)/=prec%tdim(3).or.omega2%tdim(4)/=prec%tdim(4))then
     call tensor_status_quit("ERROR(precondition_doubles_parallel):omega2 and prec have&
        &different distributions", lspdm_stdout)
  endif
  !Get the slaves to this routine
  if(me==root)then
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_COMM)
#endif
     call pdm_tensor_sync(JOB_PREC_DOUBLES_PAR,omega2,ppfock,qqfock,prec)
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif
  endif

  dims=prec%dims

  !TODO: introduce prefetching of tiles, but not so important if the offset
  !for prec and omega2 is chosen to be the same

  !do a loop over the local tiles of the preconditioned matrix and get the
  !corresponding tiles of the residual to form the preconditioned residual
  do lt=1,prec%nlti

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_COMM)
#endif

     call tensor_get_tile(omega2,prec%ti(lt)%gt,prec%ti(lt)%t,prec%ti(lt)%e)

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif


#ifdef VAR_PTR_RESHAPE
     om(1:prec%ti(lt)%d(1),1:prec%ti(lt)%d(2),1:prec%ti(lt)%d(3),1:prec%ti(lt)%d(4)) => prec%ti(lt)%t
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
     call c_f_pointer(c_loc(prec%ti(lt)%t(1)),om,prec%ti(lt)%d)
#else
     call tensor_status_quit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif

     !get offset for global indices
     call get_midx(prec%ti(lt)%gt,dims,prec%ntpm,prec%mode)
     do j=1,prec%mode
        dims(j)=(dims(j)-1)*prec%tdim(j)
     enddo

     !workaround for intel OMP error
     da = prec%ti(lt)%d(1)
     db = prec%ti(lt)%d(2)
     di = prec%ti(lt)%d(3)
     dj = prec%ti(lt)%d(4)

     !count over local indices
     !$OMP  PARALLEL DO DEFAULT(NONE)  SHARED(om,dims,&
     !$OMP  ppfock,qqfock,da,db,di,dj) PRIVATE(i,j,a,b) COLLAPSE(3)
     do j=1,dj
        do i=1,di
           do b=1,db
              do a=1,da

                 om(a,b,i,j) = om(a,b,i,j) / &
                    ( ppfock%elm2(i+dims(3),i+dims(3)) - qqfock%elm2(a+dims(1),a+dims(1)) + &
                    ppfock%elm2(j+dims(4),j+dims(4)) - qqfock%elm2(b+dims(2),b+dims(2)) )

              enddo 
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
     nullify(om)
  enddo

  !crucial barrier, wait for all slaves to finish their jobs
#ifdef TENSORS_IN_LSDALTON
  call time_start_phase(PHASE_IDLE)
#endif

  call tensor_mpi_barrier(tensor_work_comm)

#ifdef TENSORS_IN_LSDALTON
  call time_start_phase(PHASE_WORK)
#endif
#endif
  end subroutine precondition_doubles_parallel

  !> \brief calculate the dot product of two parallel distributed arrays. the
  !arrays must have the same tiling parameters, otherwise it is not implemented
  !> \author Patrick Ettenhuber
  !> \date december 2012
  function tensor_ddot_par(arr1,arr2,dest) result(res)
     implicit none
     !> the two arrays to calculate the dot-product from
     type(tensor),intent(in) :: arr1, arr2
     !> rank of the node to collect the result, -1 means all
     integer, intent(in) :: dest
     !> result
     real(tensor_dp) :: res
     real(tensor_dp),pointer :: buffer(:)
     integer :: lt,rem_els,mode
     real(tensor_dp), external :: ddot
     integer(kind=tensor_mpi_kind) :: dest_mpi
     logical :: distribution_ok
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: me, nnod
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)

#ifdef VAR_MPI
     !check if the init-types are the same
     if(arr1%access_type/=arr2%access_type)then
        call tensor_status_quit("ERROR(tensor_ddot_par):different init types of the&
           & arrays is not possible",lspdm_stdout)
     endif

     !check if the destination to collet the resut makes sense in connection with
     !the access_type
     if(arr1%access_type==AT_MASTER_ACCESS.and.dest/=0)then
        call tensor_status_quit("ERROR(tensor_ddot_par): the choice of destnation is&
           & useless",lspdm_stdout)
     endif

     !get the slaves to this routine
     if(arr1%access_type==AT_MASTER_ACCESS.and.me==root)then
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_COMM )
#endif
        call pdm_tensor_sync(JOB_DDOT_PAR,arr1,arr2)
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_WORK )
#endif
     endif

     distribution_ok = (arr1%mode==arr2%mode)
     do mode=1,arr1%mode
        if( arr1%tdim(mode) /= arr2%tdim(mode) )then
           distribution_ok = .false.
        endif
     enddo


     !check for the same distribution of the arrays
     if( distribution_ok )then

        !zeroing the result
        res    = 0.0E0_tensor_dp

        !allocate buffer for the tiles
        !TODO: introduce prefetching make preftching dependent on wrk and iwrk on input
        call tensor_alloc_mem(buffer,arr1%tsize)
        buffer=0.0E0_tensor_dp

        !loop over local tiles of array2  and get the corresponding tiles of
        !array1
        do lt=1,arr2%nlti
#ifdef TENSORS_IN_LSDALTON
           call time_start_phase( PHASE_COMM )
#endif
           call tensor_get_tile(arr1,arr2%ti(lt)%gt,buffer,arr2%ti(lt)%e)
#ifdef TENSORS_IN_LSDALTON
           call time_start_phase( PHASE_WORK )
#endif
           res = res + ddot(arr2%ti(lt)%e,arr2%ti(lt)%t,1,buffer,1)
        enddo

        call tensor_free_mem(buffer)

     else

        call tensor_status_quit("ERROR(tensor_ddot_par):NOT YET IMPLEMENTED, if the arrays have&
           & different distributions",lspdm_stdout)

     endif

     !get result on the specified node/s
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_COMM )
#endif
     if(dest==-1)then
        call tensor_mpi_allreduce(res,tensor_work_comm)
     else
        dest_mpi=dest
        call tensor_mpi_reduce(res,dest_mpi,tensor_work_comm)
     endif
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif

     if( tensor_always_sync ) call tensor_mpi_barrier(tensor_work_comm)
#else
     res = 0.0E0_tensor_dp
#endif
  end function tensor_ddot_par

  !> x = a * x + b * y
  !> \brief array addition routine for TT_TILED_DIST arrays
  !> \author Patrick Ettenhuber
  !> \date January 2013
  subroutine tensor_add_par(a,x,b,y,order)
     implicit none
     !> array to collect the result in
     type(tensor), intent(inout) :: x
     !> array to add to x
     type(tensor), intent(in) :: y
     !> scale factor without intent, because it might be overwiritten for the slaves
     real(tensor_dp),intent(in) :: a,b
     !> order y to adapt to dims of b
     integer, intent(in) :: order(x%mode)
     real(tensor_dp),pointer :: buffer(:)
     real(tensor_dp) :: prex, prey
     integer :: i,lt,nbuffs,ibuf,ibuf_idx,fbuf_idx,cmidy,buffer_lt
     integer :: xmidx(x%mode), ymidx(y%mode), ytdim(y%mode), ynels
     integer(kind=tensor_mpi_kind),pointer :: req(:)
#ifdef VAR_MPI
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: me, nnod
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)

     prex = a
     prey = b

     !check if the access_types are the same
     if(x%access_type/=y%access_type)then
        call tensor_status_quit("ERROR(tensor_add_par):different init types&
           & impossible",lspdm_stdout)
     endif

     !IF NOT AT_MASTER_ACCESS all processes should know b on call-time, else b is
     !broadcasted here
     if(x%access_type==AT_MASTER_ACCESS.and.me==root)then
        call pdm_tensor_sync(JOB_ADD_PAR,x,y)
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase(PHASE_COMM)
#endif

        call tensor_buffer(order,x%mode,root=root,comm=tensor_work_comm)
        call tensor_buffer(prex)
        call tensor_buffer(prey,finalize=.true.)
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase(PHASE_WORK)
#endif
     endif

     do i=1,x%mode
        if(x%tdim(i) /= y%tdim(order(i))) then
           call tensor_status_quit("ERROR(tensor_add_par): tdims of arrays not &
              &compatible (with the given order)",-1)
        endif
     enddo

     ! now set to two and that ought to be enough, but should work with any
     ! number >0
     if( gm_buf%init )then
        nbuffs = gm_buf%n / x%tsize
        !allocate buffer for the tiles
        buffer => gm_buf%buf
     else
        nbuffs = 2
        !allocate buffer for the tiles
        call tensor_alloc_mem(buffer,x%tsize*nbuffs)
     endif

     if( alloc_in_dummy )then
        call tensor_alloc_mem(req,nbuffs)
        call tensor_lock_wins(y,'s',all_nodes=.true.)
     endif


     !fill buffer
     do lt=1,min(nbuffs-1,x%nlti)

        call get_midx(x%ti(lt)%gt,xmidx,x%ntpm,x%mode)

        do i=1,x%mode
           ymidx(order(i)) = xmidx(i)
        enddo

        call get_tile_dim(ytdim,y,ymidx)

        ynels = 1
        do i=1,y%mode
           ynels = ynels * ytdim(i)
        enddo

        if(ynels /= x%ti(lt)%e)call tensor_status_quit("ERROR(tensor_add_par): #elements in tiles mismatch",-1)

        ibuf = mod(lt-1,nbuffs)
        ibuf_idx = ibuf*x%tsize + 1
        fbuf_idx = ibuf_idx + ynels - 1

        cmidy = get_cidx(ymidx,y%ntpm,y%mode)

#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_COMM )
#endif
        if(alloc_in_dummy )then

           call tensor_get_tile(y,ymidx,buffer(ibuf_idx:fbuf_idx),ynels,lock_set=.true.,req=req(ibuf+1))

        else
           call tensor_lock_win(y,cmidy,'s')
           call tensor_get_tile(y,ymidx,buffer(ibuf_idx:fbuf_idx),ynels,lock_set=.true.)
        endif
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_WORK )
#endif

     enddo


     !lsoop over local tiles of array x
     do lt=1,x%nlti

        !buffer last element
        buffer_lt = lt + nbuffs - 1
        if(buffer_lt <= x%nlti)then
           call get_midx(x%ti(buffer_lt)%gt,xmidx,x%ntpm,x%mode)

           do i=1,x%mode
              ymidx(order(i)) = xmidx(i)
           enddo

           call get_tile_dim(ytdim,y,ymidx)

           ynels = 1
           do i=1,y%mode
              ynels = ynels * ytdim(i)
           enddo

           if(ynels /= x%ti(buffer_lt)%e)call tensor_status_quit("ERROR(tensor_add_par): #elements in tiles mismatch",-1)

           ibuf = mod(buffer_lt-1,nbuffs)
           ibuf_idx = ibuf*x%tsize + 1
           fbuf_idx = ibuf_idx + ynels - 1 


           cmidy = get_cidx(ymidx,y%ntpm,y%mode)

#ifdef TENSORS_IN_LSDALTON
           call time_start_phase( PHASE_COMM )
#endif

           if(alloc_in_dummy )then

              call tensor_get_tile(y,ymidx,buffer(ibuf_idx:fbuf_idx),ynels,lock_set=.true.,req=req(ibuf+1))

           else

              call tensor_lock_win(y,cmidy,'s')
              call tensor_get_tile(y,ymidx,buffer(ibuf_idx:fbuf_idx),ynels,lock_set=.true.)

           endif

#ifdef TENSORS_IN_LSDALTON
           call time_start_phase( PHASE_WORK )
#endif
        endif

        call get_midx(x%ti(lt)%gt,xmidx,x%ntpm,x%mode)

        do i=1,x%mode
           ymidx(order(i)) = xmidx(i)
        enddo

        call get_tile_dim(ytdim,y,ymidx)

        ynels = 1
        do i=1,y%mode
           ynels = ynels * ytdim(i)
        enddo

        if(ynels /= x%ti(lt)%e)call tensor_status_quit("ERROR(tensor_add_par): #elements in tiles mismatch",-1)

        ibuf = mod(lt-1,nbuffs)
        ibuf_idx = ibuf*x%tsize + 1
        fbuf_idx = ibuf_idx + ynels - 1

        cmidy = get_cidx(ymidx,y%ntpm,y%mode)

#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_COMM )
#endif
        if(alloc_in_dummy )then
           call tensor_mpi_wait(req(ibuf+1))
        else
           call tensor_unlock_win(y,cmidy)
        endif
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_WORK )
#endif

        select case(x%mode)
        case(1)
           if(prex==0.0E0_tensor_dp)then

              x%ti(lt)%t = 0.0E0_tensor_dp

           else if(prex /= 1.0E0_tensor_dp) then

              call dscal(x%ti(lt)%e,prex,x%ti(lt)%t,1)

           endif

           call daxpy(x%ti(lt)%e,prey,buffer(ibuf_idx:fbuf_idx),1,x%ti(lt)%t,1)

        case(2)
           call array_reorder_2d(prey,buffer(ibuf_idx:fbuf_idx),ytdim(1),ytdim(2),order,prex,x%ti(lt)%t)
        case(3)
           call array_reorder_3d(prey,buffer(ibuf_idx:fbuf_idx),ytdim(1),ytdim(2),ytdim(3),order,prex,x%ti(lt)%t)
        case(4)
           call array_reorder_4d(prey,buffer(ibuf_idx:fbuf_idx),ytdim(1),ytdim(2),ytdim(3),ytdim(4),order,prex,x%ti(lt)%t)
        case default
           call tensor_status_quit("ERROR(tensor_add_par): mode>4 not yet implemented",-1)
        end select
     enddo

     if( alloc_in_dummy )then
        call tensor_free_mem(req)
        call tensor_unlock_wins(y,all_nodes=.true.)
     endif

     if(.not.gm_buf%init)then
        call tensor_free_mem(buffer)
     endif

     !crucial barrier, because direct memory access is used
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_IDLE )
#endif
     call tensor_mpi_barrier(tensor_work_comm)
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif

     if( tensor_always_sync ) call tensor_mpi_barrier(tensor_work_comm)
#endif
  end subroutine tensor_add_par

  !> x = a * x + b diag(d)* y[order]
  !> \brief tensor dmul for pdm tensors
  !> \author Patrick Ettenhuber
  !> \date January 2013
  subroutine tensor_dmul_par(a,x,b,d,y,order)
     implicit none
     !> array to collect the result in
     type(tensor), intent(inout) :: x
     !> array to add to x
     type(tensor), intent(in) :: y
     !> scale factor without intent, because it might be overwiritten for the slaves
     real(tensor_dp),intent(in) :: a,b,d(:)
     !> order y to adapt to dims of b
     integer, intent(in) :: order(x%mode)
     real(tensor_dp),pointer :: buffer(:)
     real(tensor_dp) :: prex, prey
     integer :: i,lt,nbuffs,ibuf,ibuf_idx,cmidy,buffer_lt, d_len,m,n,mt,nt,of
     integer :: xmidx(x%mode), ymidx(y%mode), ytdim(y%mode), ynels
     integer(kind=tensor_mpi_kind),pointer :: req(:)
#ifdef VAR_MPI
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: me, nnod
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)

     !check if the access_types are the same
     if(x%access_type/=y%access_type)then
        call tensor_status_quit("ERROR(tensor_dmul_par):different init types&
           & impossible",lspdm_stdout)
     endif
     if(x%mode/=2)then
        call tensor_status_quit("ERROR(tensor_dmul_par):not implemented for mode>2",lspdm_stdout)
     endif
     if(x%mode/=y%mode)then
        call tensor_status_quit("ERROR(tensor_dmul_par):not implemented for x%mode /= y%mode",lspdm_stdout)
     endif

     prex = a
     prey = b

     n = x%dims(1)
     m = x%dims(2)

     if (order(1) == 1 .and. order(2) == 2) then
        d_len = m
     else if (order(1) == 2 .and. order(2) == 1) then
        d_len = n
     endif


     !IF NOT AT_MASTER_ACCESS all processes should know b on call-time, else b is
     !broadcasted here
     if(x%access_type==AT_MASTER_ACCESS.and.me==root)then
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase(PHASE_COMM)
#endif
        call pdm_tensor_sync(JOB_DMUL_PAR,x,y)

        call tensor_buffer(order,x%mode,root=root,comm=tensor_work_comm)
        call tensor_buffer(prex)
        call tensor_buffer(prey)
        call tensor_buffer(d_len)
        call tensor_buffer(d(1:d_len),d_len,finalize=.true.)
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase(PHASE_WORK)
#endif
     endif

     do i=1,x%mode
        if(x%tdim(i) /= y%tdim(order(i))) then
           call tensor_status_quit("ERROR(tensor_dmul_par): tdims of arrays not &
              &compatible (with the given order)",-1)
        endif
     enddo

     ! now set to two and that ought to be enough, but should work with any
     ! number >0
     if( gm_buf%init )then
        nbuffs = gm_buf%n / x%tsize
        !allocate buffer for the tiles
        buffer => gm_buf%buf
     else
        nbuffs = 2
        !allocate buffer for the tiles
        call tensor_alloc_mem(buffer,x%tsize*nbuffs)
     endif

     if( alloc_in_dummy )then
        call tensor_alloc_mem(req,nbuffs)
        call tensor_lock_wins(y,'s',all_nodes=.true.)
     endif


     !fill buffer
     do lt=1,min(nbuffs-1,x%nlti)

        call get_midx(x%ti(lt)%gt,xmidx,x%ntpm,x%mode)

        do i=1,x%mode
           ymidx(order(i)) = xmidx(i)
        enddo

        call get_tile_dim(ytdim,y,ymidx)

        ynels = 1
        do i=1,y%mode
           ynels = ynels * ytdim(i)
        enddo

        if(ynels /= x%ti(lt)%e)call tensor_status_quit("ERROR(tensor_dmul_par): #elements in tiles mismatch",-1)

        ibuf = mod(lt-1,nbuffs)
        ibuf_idx = ibuf*x%tsize + 1

        cmidy = get_cidx(ymidx,y%ntpm,y%mode)

#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_COMM )
#endif
        if(alloc_in_dummy )then

           call tensor_get_tile(y,ymidx,buffer(ibuf_idx:),ynels,lock_set=.true.,req=req(ibuf+1))

        else
           call tensor_lock_win(y,cmidy,'s')
           call tensor_get_tile(y,ymidx,buffer(ibuf_idx:),ynels,lock_set=.true.)
        endif
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_WORK )
#endif

     enddo


     !lsoop over local tiles of array x
     do lt=1,x%nlti

        !buffer last element
        buffer_lt = lt + nbuffs - 1
        if(buffer_lt <= x%nlti)then
           call get_midx(x%ti(buffer_lt)%gt,xmidx,x%ntpm,x%mode)

           do i=1,x%mode
              ymidx(order(i)) = xmidx(i)
           enddo

           call get_tile_dim(ytdim,y,ymidx)

           ynels = 1
           do i=1,y%mode
              ynels = ynels * ytdim(i)
           enddo

           if(ynels /= x%ti(buffer_lt)%e)call tensor_status_quit("ERROR(tensor_dmul_par): #elements in tiles mismatch",-1)

           ibuf = mod(buffer_lt-1,nbuffs)
           ibuf_idx = ibuf*x%tsize + 1

           cmidy = get_cidx(ymidx,y%ntpm,y%mode)

#ifdef TENSORS_IN_LSDALTON
           call time_start_phase( PHASE_COMM )
#endif

           if(alloc_in_dummy )then

              call tensor_get_tile(y,ymidx,buffer(ibuf_idx:),ynels,lock_set=.true.,req=req(ibuf+1))

           else

              call tensor_lock_win(y,cmidy,'s')
              call tensor_get_tile(y,ymidx,buffer(ibuf_idx:),ynels,lock_set=.true.)

           endif

#ifdef TENSORS_IN_LSDALTON
           call time_start_phase( PHASE_WORK )
#endif
        endif

        call get_midx(x%ti(lt)%gt,xmidx,x%ntpm,x%mode)

        do i=1,x%mode
           ymidx(order(i)) = xmidx(i)
        enddo

        call get_tile_dim(ytdim,y,ymidx)

        ynels = 1
        do i=1,y%mode
           ynels = ynels * ytdim(i)
        enddo

        if(ynels /= x%ti(lt)%e)call tensor_status_quit("ERROR(tensor_dmul_par): #elements in tiles mismatch",-1)

        ibuf = mod(lt-1,nbuffs)
        ibuf_idx = ibuf*x%tsize + 1

        cmidy = get_cidx(ymidx,y%ntpm,y%mode)

#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_COMM )
#endif
        if(alloc_in_dummy )then
           call tensor_mpi_wait(req(ibuf+1))
        else
           call tensor_unlock_win(y,cmidy)
        endif
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_WORK )
#endif

        if(abs(prex)<1.0E-15)then
           x%ti(lt)%t = 0.0E0_tensor_dp
        else if(abs(prex-1.0E0_tensor_dp)>1.0E-15)then
           call dscal(x%ti(lt)%e,prex,x%ti(lt)%t,1)
        endif

        nt = x%ti(lt)%d(1)
        mt = x%ti(lt)%d(2)

        of = (xmidx(order(1))-1)*x%tdim(order(1)) 

        if (order(1) == 1 .and. order(2) == 2) then

           do i=1,nt
              call daxpy(mt,prey*d(of+i),buffer(ibuf_idx+i-1),nt,x%ti(lt)%t(i),nt)
           enddo

        else if (order(1)==2 .and. order(2)==1) then

           do i=1,mt
              call daxpy(nt,prey*d(of+i),buffer(ibuf_idx+nt*(i-1)),1,x%ti(lt)%t(i),mt)
           enddo 

        else                                                                                                                        
           call tensor_status_quit("ERROR(tensor_dmul_par): wrong order",-1)                                                                        
        end if

     enddo

     if( alloc_in_dummy )then
        call tensor_free_mem(req)
        call tensor_unlock_wins(y,all_nodes=.true.)
     endif

     if(.not.gm_buf%init)then
        call tensor_free_mem(buffer)
     endif

     !crucial barrier, because direct memory access is used
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_IDLE )
#endif
     call tensor_mpi_barrier(tensor_work_comm)
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
#endif
  end subroutine tensor_dmul_par

  !> \brief Hadamard product Cij = alpha*Aij*Bij+beta*Cij for TT_TILED_DIST arrays
  !> \author Thomas Kjærgaard and Patrick Ettenhuber
  !> \date 2015
  subroutine tensor_hmul_par(alpha,A,B,beta,C)
     implicit none
     !> array input, this is the result array with overwritten data
     type(tensor),intent(inout) :: C
     type(tensor),intent(in) :: A,B
     !> scaling factor for array C
     real(tensor_dp),intent(in) :: beta
     !> scaling factor for array A and B
     real(tensor_dp),intent(in) :: alpha
     integer :: i,lt
#ifdef VAR_MPI
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: me, nnod
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)

     !check if the access_types are the same
     if(A%access_type/=C%access_type)then
        call tensor_status_quit("ERROR(tensor_hmul_par):different init types&
           & impossible",lspdm_stdout)
     endif
     if(B%access_type/=C%access_type)then
        call tensor_status_quit("ERROR(tensor_hmul_par):different init types&
           & impossible",lspdm_stdout)
     endif

     !IF AT_MASTER_ACCESS all processes should know b on call-time, else b is
     !broadcasted here
     if(C%access_type==AT_MASTER_ACCESS.and.me==root)then
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase(PHASE_COMM)
#endif
        call pdm_tensor_sync(JOB_HMUL_PAR,A,B,C)

        call tensor_buffer(alpha,root=root,comm=tensor_work_comm)
        call tensor_buffer(beta, finalize=.true.)
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase(PHASE_WORK)
#endif
     endif

     IF(A%offset.NE.C%offset)THEN
        call tensor_status_quit('Offset not same in tensor_hmul_par',-1)
     ENDIF
     IF(B%offset.NE.C%offset)THEN
        call tensor_status_quit('Offset not same in tensor_hmul_par',-1)
     ENDIF

     !lsoop over local tiles of array C
     do lt=1,C%nlti
        call HMULSUB(C%ti(lt)%e,alpha,A%ti(lt)%t,B%ti(lt)%t,beta,C%ti(lt)%t)
     enddo

     !   if( tensor_always_sync )then
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_IDLE )
#endif
     call tensor_mpi_barrier(tensor_work_comm)
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
     !   endif
#endif
  end subroutine tensor_hmul_par

  subroutine HMULSUB(N,alpha,A,B,beta,C)
     implicit none
     integer(kind=tensor_long_int),intent(in) :: N 
     real(tensor_dp),intent(in) :: alpha,beta
     real(tensor_dp),intent(in) :: A(N),B(N)
     real(tensor_dp),intent(inout) :: C(N)
     integer(kind=tensor_long_int) :: i        
     !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i) SHARED(N,C,alpha,A,B,beta)
     do i = 1,N
        C(i) = alpha*A(i)*B(i) + beta*C(i)
     enddo
     !$OMP END PARALLEL DO
  end subroutine HMULSUB

  !> \brief array copying routine for TT_TILED_DIST arrays
  !> \author Patrick Ettenhuber
  !> \date January 2013
  subroutine tensor_cp_tiled(from,to_ar, order )
     implicit none
     !> source, array to copy
     type(tensor), intent(in) :: from
     !> drain, the copied array
     type(tensor), intent(inout) :: to_ar
     integer, intent(in) :: order(to_ar%mode)
#ifdef VAR_PTR_RESHAPE
     real(tensor_dp), pointer, contiguous :: buffer(:,:)
#else
     real(tensor_dp), pointer :: buffer(:,:)
#endif
     real(tensor_dp), parameter :: prex = 0.0E0_tensor_dp
     real(tensor_dp), parameter :: prey = 1.0E0_tensor_dp
     integer :: i,lt,nbuffs,ibuf,cmidy,buffer_lt
     integer :: xmidx(from%mode), ymidx(to_ar%mode), ytdim(to_ar%mode), ynels
#ifdef VAR_MPI
     integer(kind=tensor_mpi_kind), pointer :: req(:)
     integer :: order_comm(to_ar%mode)
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: me, nnod
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)

     !if(present(order)) call tensor_status_quit("ERROR(tensor_cp_tiled): order not yet implemented",-1)

     !check for the same access_types
     if(from%access_type/=to_ar%access_type)then
        call tensor_status_quit("ERROR(tensor_cp_tiled):different init types&
           & impossible",lspdm_stdout)
     endif

     order_comm = order

     !get the slaves
     if(from%access_type==AT_MASTER_ACCESS.and.me==root)then
        call pdm_tensor_sync(JOB_CP_ARR,from,to_ar)
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase(PHASE_COMM)
#endif
        call tensor_mpi_bcast(order_comm,from%mode,root,tensor_work_comm)
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase(PHASE_WORK)
#endif
     endif

     ! now set to two and that ought to be enough, but should work with any
     ! number >0
     nbuffs = 2

     !allocate buffer for the tiles
     call tensor_alloc_mem(buffer,to_ar%tsize,nbuffs)
     if( alloc_in_dummy )then
        call tensor_alloc_mem(req,nbuffs)
        call tensor_lock_wins(from,'s',all_nodes=.true.)
     endif

     !fill buffer
     do lt=1,min(nbuffs-1,to_ar%nlti)

        call get_midx(to_ar%ti(lt)%gt,xmidx,to_ar%ntpm,to_ar%mode)

        do i=1,to_ar%mode
           ymidx(order(i)) = xmidx(i)
        enddo

        call get_tile_dim(ytdim,from,ymidx)

        ynels = 1
        do i=1,from%mode
           ynels = ynels * ytdim(i)
        enddo

        if(ynels /= to_ar%ti(lt)%e)call tensor_status_quit("ERROR(tensor_cp_tiled): #elements in tiles mismatch",-1)

        ibuf = mod(lt-1,nbuffs)+1

        cmidy = get_cidx(ymidx,from%ntpm,from%mode)

        if( alloc_in_dummy )then

           call tensor_get_tile(from,ymidx,buffer(:,ibuf),ynels,lock_set=.true.,req=req(ibuf))

        else
           call tensor_lock_win(from,cmidy,'s')
           call tensor_get_tile(from,ymidx,buffer(:,ibuf),ynels,lock_set=.true.)
        endif
     enddo

     !lsoop over local tiles of array to_ar
     do lt=1,to_ar%nlti

        !buffer last element
        buffer_lt = lt + nbuffs - 1
        if(buffer_lt <= to_ar%nlti)then

           call get_midx(to_ar%ti(buffer_lt)%gt,xmidx,to_ar%ntpm,to_ar%mode)

           do i=1,to_ar%mode
              ymidx(order(i)) = xmidx(i)
           enddo

           call get_tile_dim(ytdim,from,ymidx)

           ynels = 1
           do i=1,from%mode
              ynels = ynels * ytdim(i)
           enddo

           if(ynels /= to_ar%ti(buffer_lt)%e)call tensor_status_quit("ERROR(tensor_cp_tiled): #elements in tiles mismatch",-1)

           ibuf = mod(buffer_lt-1,nbuffs)+1


           cmidy = get_cidx(ymidx,from%ntpm,from%mode)

           if( alloc_in_dummy )then

              call tensor_get_tile(from,ymidx,buffer(:,ibuf),ynels,lock_set=.true.,req=req(ibuf))

           else
              call tensor_lock_win(from,cmidy,'s')
              call tensor_get_tile(from,ymidx,buffer(:,ibuf),ynels,lock_set=.true.)
           endif
        endif

        call get_midx(to_ar%ti(lt)%gt,xmidx,to_ar%ntpm,to_ar%mode)

        do i=1,to_ar%mode
           ymidx(order(i)) = xmidx(i)
        enddo

        call get_tile_dim(ytdim,from,ymidx)

        ynels = 1
        do i=1,from%mode
           ynels = ynels * ytdim(i)
        enddo

        if(ynels /= to_ar%ti(lt)%e)call tensor_status_quit("ERROR(tensor_cp_tiled): #elements in tiles mismatch",-1)

        ibuf = mod(lt-1,nbuffs)+1

        cmidy = get_cidx(ymidx,from%ntpm,from%mode)

        if( alloc_in_dummy )then
           call tensor_mpi_wait(req(ibuf))
        else
           call tensor_unlock_win(from,cmidy)
        endif

        select case(to_ar%mode)
        case(1)
           call dcopy(to_ar%ti(lt)%e,buffer(:,ibuf),1,to_ar%ti(lt)%t,1)
        case(2)
           call array_reorder_2d(prey,buffer(:,ibuf),ytdim(1),ytdim(2),order,prex,to_ar%ti(lt)%t)
        case(3)
           call array_reorder_3d(prey,buffer(:,ibuf),ytdim(1),ytdim(2),ytdim(3),order,prex,to_ar%ti(lt)%t)
        case(4)
           call array_reorder_4d(prey,buffer(:,ibuf),ytdim(1),ytdim(2),ytdim(3),ytdim(4),order,prex,to_ar%ti(lt)%t)
        case default
           call tensor_status_quit("ERROR(tensor_cp_tiled): mode>4 not yet implemented",-1)
        end select
     enddo

     call tensor_free_mem(buffer)

     if( alloc_in_dummy )then
        call tensor_free_mem(req)
        call tensor_unlock_wins(from,all_nodes=.true.)
     endif

     !crucial barrier, because direct memory access is used
     call tensor_mpi_barrier(tensor_work_comm)
     if( tensor_always_sync ) call tensor_mpi_barrier(tensor_work_comm)
#endif
  end subroutine tensor_cp_tiled


  !> \brief zeroing routine for tiled distributed arrays
  !> \author Patrick Ettenhuber
  !> \date late 2012
  subroutine tensor_zero_tiled_dist(a)
     implicit none
     !> array to zero
     type(tensor),intent(inout) :: a
     integer :: lt
#ifdef VAR_MPI
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: me, nnod
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)
     !get the slaves here
     if(a%access_type==AT_MASTER_ACCESS.and.me==root)then
        call pdm_tensor_sync(JOB_tensor_ZERO,a)
     endif

     !loop over local tiles and zero them individually
     do lt=1,a%nlti
        a%ti(lt)%t=0.0E0_tensor_dp
     enddo

     if( tensor_always_sync ) call tensor_mpi_barrier(tensor_work_comm)
#endif
  end subroutine tensor_zero_tiled_dist
  !> \brief randomizing routine for tiled distributed arrays
  !> \author Patrick Ettenhuber
  !> \date beginning 2015
  subroutine tensor_rand_tiled_dist(a)
     implicit none
     !> array to zero
     type(tensor),intent(inout) :: a
     integer :: lt
#ifdef VAR_MPI
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: me, nnod
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)
     !get the slaves here
     if(a%access_type==AT_MASTER_ACCESS.and.me==root)then
        call pdm_tensor_sync(JOB_tensor_rand,a)
     endif

     call random_seed()
     !loop over local tiles and zero them individually
     do lt=1,a%nlti
        call random_number(a%ti(lt)%t)
     enddo

     if( tensor_always_sync ) call tensor_mpi_barrier(tensor_work_comm)
#endif
  end subroutine tensor_rand_tiled_dist


  !> \author Patrick Ettenhuber
  !> \date January 2013
  !> \brief initialized a replicated matrix on each node
  subroutine tensor_init_replicated(arr,dims,nmodes,pdm,bg)
     implicit none
     !> array to be initialilzed
     type(tensor),intent(inout) :: arr
     !> number of modes and the dimensions of the array
     integer,intent(in) :: nmodes,dims(nmodes)
     !> integer specifying the access_type of the array
     integer(kind=tensor_standard_int),intent(in) :: pdm
     logical, intent(in) :: bg
     integer(kind=tensor_long_int) :: i,j
     integer :: addr
     integer(kind=tensor_standard_int) :: tdimdummy(nmodes)
     integer :: nelms
     integer, pointer :: lg_buf(:),pc_buf(:)
     logical :: master, bg_int
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: nnod
     integer(kind=tensor_mpi_kind) :: me
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)
     call tensor_get_size(nnod)

     !set the initial values and overwrite them later
     master = (me == root)
     bg_int = bg

     !allocate all pdm in p_arr therefore get free address and associate it with
     !the array, and increment the array counter
     p_arr%curr_addr_on_node   = get_free_address(.true.)
     addr                      = p_arr%curr_addr_on_node
     p_arr%arrays_in_use       = p_arr%arrays_in_use + 1
     p_arr%a(addr)%local_addr  = addr
     p_arr%a(addr)%initialized = .true.
#ifdef VAR_MPI
     p_arr%a(addr)%nnod        = nnod
     p_arr%a(addr)%comm        = tensor_work_comm
#else
     !set to invalid, since not used here
     p_arr%a(addr)%nnod        = -1
     p_arr%a(addr)%comm        = -1
#endif

     allocate( p_arr%a(addr)%access_type )
     p_arr%a(addr)%access_type = pdm

     !SET MODE
     p_arr%a(addr)%mode = nmodes

     !SET DIMS
     call tensor_set_dims(p_arr%a(addr),dims,nmodes)

     !SET ARRAY TYPE
     p_arr%a(addr)%itype=TT_REPLICATED

     !SET NELMS
     nelms=1
     do i=1,nmodes
        nelms=nelms*dims(i)
     enddo
     p_arr%a(addr)%nelms=nelms

     !put 0 in tdim, since for the replicated array it is not important
     tdimdummy=0
     call tensor_set_tdims(p_arr%a(addr),tdimdummy,nmodes)
     !SET NELMS

     !In the initialization the addess has to be set, since pdm_tensor_sync
     !depends on the  adresses, but setting them correctly is done later
     call tensor_alloc_mem(lg_buf,nnod)
     lg_buf = 0

     !if master init only master has to init the addresses addresses before
     !pdm syncronization
     if(master .and. p_arr%a(addr)%access_type==AT_MASTER_ACCESS)then
        call tensor_set_addr(p_arr%a(addr),lg_buf,nnod)
#ifdef VAR_MPI
        call pdm_tensor_sync(JOB_INIT_TENSOR_REPLICATED,p_arr%a(addr))
#endif
     endif

     !if AT_ALL_ACCESS all have to have the addresses allocated
     if(p_arr%a(addr)%access_type==AT_ALL_ACCESS)then
        call tensor_set_addr(p_arr%a(addr),lg_buf,nnod)
        !if(lspdm_use_comm_proc)call tensor_set_addr(p_arr%a(addr),pc_buf,pc_nnodes,.true.)
     endif

#ifdef VAR_MPI
     !SET THE ADDRESSES ON ALL NODES     
     lg_buf(me+1)=addr 
     call tensor_mpi_allreduce(lg_buf,nnod,tensor_work_comm)

     if( p_arr%a(addr)%access_type==AT_MASTER_ACCESS)then
        call tensor_mpi_bcast(bg_int,root,tensor_work_comm)
     endif

     call tensor_set_addr(p_arr%a(addr),lg_buf,nnod)
#endif

     !ALLOCATE STORAGE SPACE FOR THE ARRAY
     call memory_allocate_tensor_dense(p_arr%a(addr),bg_int)

     !RETURN THE CURRENLY ALLOCATE ARRAY
     arr=p_arr%a(addr)

     call tensor_free_mem(lg_buf)
     !if(lspdm_use_comm_proc)call tensor_free_mem(pc_buf)
  end subroutine tensor_init_replicated


  !> \brief print the norm of a replicated array from each node, just a
  !debugging routine
  !> \author Patrick Ettenhuber
  !> \date January 2012
  function tensor_print_norm_repl(arr) result(nrm)
     implicit none
     !> replicated array to print the norm from
     type(tensor), intent(in) :: arr
     !return-value is the norm
     real(tensor_dp) :: nrm
     integer :: i
#ifdef VAR_MPI
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: nnod
     integer(kind=tensor_mpi_kind) :: me
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)
     call tensor_get_size(nnod)

     !get the slaves
     if(me==root.and.arr%access_type==AT_MASTER_ACCESS)then
        call pdm_tensor_sync(JOB_GET_NORM_REPLICATED,arr)
     endif

     !zero the norm an calculate it
     nrm =0.0E0_tensor_dp
     do i=1,arr%nelms
        nrm=nrm+arr%elm1(i)*arr%elm1(i)
     enddo
     print *,"on nodes",me,sqrt(nrm)
#else
     nrm = 0.0E0_tensor_dp
#endif
  end function tensor_print_norm_repl


  !> \brief synchronize a replicated array from a source
  !> \author Patrick Ettenhuber
  !> \date cannot remember, 2012
  subroutine tensor_sync_replicated(arr,fromnode)
     implicit none
     !> array to synchronize
     type(tensor), intent(inout) :: arr
     !> specify the node which holds the original data that should be
     !synchronized to all nodes
     integer,optional, intent(in) :: fromnode
     integer(kind=tensor_mpi_kind) :: source
#ifdef VAR_MPI
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: nnod
     integer(kind=tensor_mpi_kind) :: me
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)
     call tensor_get_size(nnod)

     !give meaningful quit statement for useless input
     if(present(fromnode).and.arr%access_type==AT_MASTER_ACCESS)then
        call tensor_status_quit("ERROR(tensor_sync_replicated): This combintion of input&
           &elements does not give sense",lspdm_stdout)
        ! why would you want to collect the data on a node you cannot direcly
        ! access, or if you can access the data in the calling subroutine on the
        ! specified node, why is the init_tyep AT_MASTER_ACCESS?
     endif

     ! get slaves
     if(me==root.and.arr%access_type==AT_MASTER_ACCESS)then
        call pdm_tensor_sync(JOB_SYNC_REPLICATED,arr)
     endif


     !specify the source of the data, by default master
     source = root
     if(present(fromnode))source=fromnode

     !do the synchronization
     call tensor_mpi_bcast(arr%elm1,arr%nelms,source,tensor_work_comm)
#endif    
  end subroutine tensor_sync_replicated


  !> \brief calculate the default tile-dimensions for the tiled dirtributed
  !array. default tile dimensions are currently a compromize between data
  !distribution and communication cost, it is cheaper to communicate large
  !chunks than many of them
  !> \author Patrick Ettenhuber
  !> \date march 2013
  subroutine tensor_default_batches(dims,nmodes,tdim,div)
     implicit none
     !> mode of the array
     integer :: nmodes
     !> dimensions in the modes
     integer :: dims(nmodes)
     !> divisor the last dimension whic is slict
     integer,intent(inout) :: div
     !> tdim output 
     integer(kind=tensor_standard_int) :: tdim(nmodes)
     integer :: i,j
     integer(kind=tensor_mpi_kind) :: nlocalnodes
     integer :: cdims
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: me
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)
     call tensor_get_size(nlocalnodes)

     tdim = DEFAULT_TDIM

  end subroutine tensor_default_batches


  !> \author Patrick Ettenhuber
  !> \date September 2012
  !> \brief initialized a distributed tiled array
  subroutine tensor_init_tiled(arr,dims,nmodes,at,it,pdm,bg,tdims,ps_d,force_offset)
     implicit none
     type(tensor),intent(inout) :: arr
     integer,intent(in) :: nmodes,dims(nmodes)
     character(4) :: at
     integer(kind=tensor_standard_int) :: it
     integer(kind=tensor_standard_int) :: pdm
     logical, intent(in) :: bg
     integer,optional :: tdims(nmodes)
     logical, optional :: ps_d
     integer,intent(in), optional :: force_offset
     integer(kind=tensor_long_int) :: i,j
     integer ::addr,pdmt,k,div
     integer(kind=tensor_standard_int) :: dflt(nmodes)
     integer :: cdims
     integer, pointer :: lg_buf(:)
     logical :: master,defdims, pseudo_dense
     integer :: infobuf(2),fo
     logical,parameter :: zeros_in_tiles=.false.
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: nnod
     integer(kind=tensor_mpi_kind) :: me
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)
     call tensor_get_size(nnod)

     !set the initial values and overwrite them later
     master = (me==root)

     pseudo_dense = .false.
     if(present(ps_d))pseudo_dense=ps_d

     fo = -1
     if(present(force_offset))then
        fo = force_offset
     endif

     !allocate all tiled arrays in p_arr, get free
     p_arr%curr_addr_on_node   = get_free_address(.true.)
     addr                      = p_arr%curr_addr_on_node
     p_arr%arrays_in_use       = p_arr%arrays_in_use + 1
     p_arr%a(addr)%local_addr  = addr
     p_arr%a(addr)%initialized = .true.
     p_arr%a(addr)%bg_alloc    = bg
#ifdef VAR_MPI
     p_arr%a(addr)%nnod        = nnod
     p_arr%a(addr)%comm        = tensor_work_comm
#else
     !set to invalid, since not used here
     p_arr%a(addr)%nnod        = -1
     p_arr%a(addr)%comm        = -1
#endif

     allocate( p_arr%a(addr)%access_type )
     p_arr%a(addr)%access_type = pdm

     !INITIALIZE TILE STRUCTURE, if master from basics, if slave most is already
     !there
     defdims = .false.
     p_arr%a(addr)%zeros = zeros_in_tiles

     !SET MODE
     p_arr%a(addr)%mode = nmodes

     !SET DIMS
     call tensor_set_dims(p_arr%a(addr),dims,nmodes)

     if(present(tdims))then
        dflt=tdims
     else
        defdims=.true.
        !insert a routine for estimation of ts according to mem here
     endif

     !check if invalid numbers occur and fall back to default if so
     if(.not.defdims.and.present(tdims))then
        do i=1,nmodes
           if(dflt(i)<=0)then
              call tensor_status_quit("ERROR(tensor_init_tiled): invalid tile request, one of the dimensions is <= 0",33)
           endif
        enddo
     endif

     !if needed, get default batch sizes, which are chosen such, that the best
     !distribution in terms of transfer speed and even distribution occur 
     !-> lots of consecutive !elements, big tiles, enough tiles
     if(defdims)then
        call tensor_default_batches(dims,nmodes,dflt,div)
     endif
     call tensor_set_tdims(p_arr%a(addr),dflt,int(p_arr%a(addr)%mode))
     if (master) then 
        p_arr%a(addr)%itype=it
        p_arr%a(addr)%atype=at
     end if    

     !divide A into tiles, according to dimensions
     !begin with counting the number of tiles needed in each mode
     dflt=0
     p_arr%a(addr)%nelms=1
     call tensor_get_ntpm(p_arr%a(addr)%dims,p_arr%a(addr)%tdim,p_arr%a(addr)%mode,dflt)
     do i=1,p_arr%a(addr)%mode
        p_arr%a(addr)%nelms = p_arr%a(addr)%nelms * p_arr%a(addr)%dims(i)
     enddo
     call tensor_set_ntpm(p_arr%a(addr),dflt,int(p_arr%a(addr)%mode))
     !count the total number of tiles for the array and allocate in structure
     !calculate tilesize

     p_arr%a(addr)%ntiles = 1
     p_arr%a(addr)%tsize  = 1
     do i=1,p_arr%a(addr)%mode
        p_arr%a(addr)%ntiles = p_arr%a(addr)%ntiles * p_arr%a(addr)%ntpm(i)
        p_arr%a(addr)%tsize  = p_arr%a(addr)%tsize  * min(p_arr%a(addr)%tdim(i),p_arr%a(addr)%dims(i))
     enddo

     !Adapt background buffer
     if( gm_buf%init) call lspdm_reinit_global_buffer(2*long1*p_arr%a(addr)%tsize)

     !In the initialization the addess has to be set, since pdm_tensor_sync
     !depends on the  adresses, but setting them correctly is done later
     call tensor_alloc_mem(lg_buf,2*nnod)
     lg_buf = 0

     !if master init only master has to get addresses
     if(master .and. p_arr%a(addr)%access_type==AT_MASTER_ACCESS)then
        call tensor_set_addr(p_arr%a(addr),lg_buf,nnod)
#ifdef VAR_MPI
        call pdm_tensor_sync(JOB_INIT_TENSOR_TILED,p_arr%a(addr))
#endif
     endif

#ifdef VAR_MPI
     call tensor_buffer(p_arr%a(addr)%itype,&
        & init_size = tensor_int + tensor_standard_int + len(p_arr%a(addr)%atype)*tensor_char_size + tensor_log,&
        & root = root, comm = tensor_work_comm)
     call tensor_buffer( fo )
     call tensor_buffer( p_arr%a(addr)%atype )
     call tensor_buffer( p_arr%a(addr)%bg_alloc, finalize=.true. )
#endif

     !if AT_ALL_ACCESS only all have to know the addresses
     if(p_arr%a(addr)%access_type==AT_ALL_ACCESS)call tensor_set_addr(p_arr%a(addr),lg_buf,nnod)

     if(fo==-1)then
        call get_distribution_info(p_arr%a(addr))
     else
        call get_distribution_info(p_arr%a(addr),force_offset = fo)
     endif

#ifdef VAR_MPI
     lg_buf(me+1)=addr 
     lg_buf(nnod+me+1)=p_arr%a(addr)%offset
     call tensor_mpi_allreduce(lg_buf,2*nnod,tensor_work_comm)

     call tensor_set_addr(p_arr%a(addr),lg_buf,nnod)

     do i=1,nnod
        if(lg_buf(nnod+i)/=p_arr%a(addr)%offset)then
           print * ,me,"found",lg_buf(nnod+i),p_arr%a(addr)%offset,i,fo
           call tensor_status_quit("ERROR(tensor_init_tiled):offset &
              &is not the same on all nodes",lspdm_stdout)
        endif
     enddo
#endif
     call tensor_free_mem(lg_buf)

     call tensor_init_lock_set(p_arr%a(addr))
     call memory_allocate_tiles(p_arr%a(addr),p_arr%a(addr)%bg_alloc)

     if(pseudo_dense .and. (master.or.p_arr%a(addr)%access_type==AT_ALL_ACCESS))then
        call memory_allocate_tensor_dense(p_arr%a(addr),p_arr%a(addr)%bg_alloc)
     endif

     arr = p_arr%a(addr)

     !if(lspdm_use_comm_proc)call tensor_free_mem(pc_buf)
  end subroutine tensor_init_tiled

  subroutine lspdm_tensor_contract_simple(pre1,A,B,m2cA,m2cB,nmodes2c,pre2,C,order,mem,wrk,iwrk,force_sync)
     implicit none
     real(tensor_dp), intent(in)    :: pre1,pre2
     type(tensor), intent(in)    :: A,B
     integer, intent(in)        :: nmodes2c
     integer, intent(in)        :: m2cA(nmodes2c),m2cB(nmodes2c)
     type(tensor), intent(inout) :: C
     integer, intent(inout)     :: order(C%mode)
     real(tensor_dp), intent(in),    optional :: mem !in GB
     real(tensor_dp), intent(inout), target, optional :: wrk(:)
     integer(kind=tensor_long_int), intent(in),     optional :: iwrk
     logical, intent(in),        optional :: force_sync
     !internal variables
     integer :: use_wrk_space
     logical :: test_all_master_access,test_all_all_access,master,contraction_mode,sync
     real(tensor_dp), pointer :: buffA(:,:),buffB(:,:),wA(:),wB(:),wC(:),tA(:),tB(:),tC(:)
     integer :: ibufA, ibufB, nbuffsA,nbuffsB, nbuffs, buffer_cm, ntens_to_get_from, tsizeB
     integer :: gc, gm(C%mode), ro(C%mode), locC
     integer :: mA(A%mode), mB(B%mode), tdimA(A%mode), tdimB(B%mode), ordA(A%mode), ordB(B%mode)
     integer :: cmidA, cmidB
     integer(kind=tensor_standard_int) :: fBtdim(B%mode), ntpmB(B%mode)
     integer :: tdim_ord(B%mode)
     integer :: tdimC(C%mode),tdim_product(C%mode)
     integer :: nelmsTA, nelmsTB
     integer :: i,j,k,l, sq, cci, max_mode_ci(nmodes2c),cm,current_mode(nmodes2c)
     integer :: m_gemm, n_gemm, k_gemm
     logical :: B_dense,locA,locB, found
     integer(kind=tensor_mpi_kind),pointer :: reqA(:),reqB(:)
     integer, pointer :: buf_posA(:), buf_posB(:)
     logical, pointer :: buf_logA(:), buf_logB(:), nfA(:), nfB(:)
     real(tensor_dp), pointer :: w(:)
     integer(kind=tensor_long_int) :: itest
     integer, parameter :: USE_INPUT_WORK      = 334
     integer, parameter :: USE_GLOBAL_BUFFER   = 335
     integer, parameter :: USE_INTERNAL_ALLOC  = 336
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: nnod
     integer(kind=tensor_mpi_kind) :: me
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)
     call tensor_get_size(nnod)

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif

     sync = .false.
     if(present(force_sync))sync = force_sync

     B_dense = (B%itype == TT_DENSE) .or. (B%itype == TT_REPLICATED)

     if(B_dense)then
        ntens_to_get_from = 1
     else
        ntens_to_get_from = 2
     endif

#ifdef VAR_MPI
     master = (me == root)

     test_all_master_access = (A%access_type == AT_MASTER_ACCESS).and.&
        &((B%access_type == AT_MASTER_ACCESS) .or. (B%itype == TT_DENSE)).and.&
        &(C%access_type == AT_MASTER_ACCESS)
     test_all_all_access = (A%access_type == AT_ALL_ACCESS).and.&
        &((B%access_type == AT_ALL_ACCESS) .or. (B%itype == TT_DENSE)).and.&
        &(C%access_type == AT_ALL_ACCESS)

     if(  (.not.test_all_master_access.and..not.test_all_all_access) .or. &
        & (     test_all_master_access.and.     test_all_all_access)  )then
     call tensor_status_quit("ERROR(lspdm_tensor_contract_simple):: Invalid access types",-1)
  endif

  if(test_all_master_access.and.(present(wrk).or.present(iwrk)))then
     print *,"WARNING(lspdm_tensor_contract_simple): in master access ignoring, wrk and iwrk"
  endif

  !calculate reverse order
  do i = 1, C%mode
     ro(order(i)) = i
  enddo


  if(B_dense)then

     !Set contraction modes
     do i=1,nmodes2c
        fBtdim(m2cB(i)) = A%tdim(m2cA(i))
     enddo

     !Set uncontracted modes
     k = A%mode - nmodes2c + 1
     do i = 1, B%mode
        contraction_mode=.false.
        do j=1,nmodes2c
           contraction_mode = contraction_mode.or.(m2cB(j) == i)
        enddo
        if(.not.contraction_mode)then
           fBtdim(i) = C%tdim(ro(k))
           k = k + 1
        endif
     enddo

     tsizeB = 1
     call tensor_get_ntpm(B%dims,fBtdim,B%mode,ntpmB)
     do i=1,B%mode
        tsizeB = tsizeB * fBtdim(i)
     enddo

  else

     ntpmB  = B%ntpm
     fBtdim = B%tdim
     tsizeB = B%tsize

  endif

  !calculate the combined contraction index by looping over the contraction
  !modes and the respective number of tiles in these
  cci = 1
  max_mode_ci = 0
  do cm=1,nmodes2c
     if(A%ntpm(m2cA(cm))/=ntpmB(m2cB(cm)))then
        call tensor_status_quit("ERROR(lspdm_tensor_contract_simple):: A and B do not have the same &
           &ntpm in the given contraction modes",-1)
     else
        cci = cci * A%ntpm(m2cA(cm))
        max_mode_ci(cm) = A%ntpm(m2cA(cm))
     endif
  enddo


#ifdef TENSORS_IN_LSDALTON
  call time_start_phase( PHASE_COMM )
#endif
  if(master.and.test_all_master_access)then
     if(B%itype == TT_TILED_DIST .or. B%itype == TT_REPLICATED)then
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase(PHASE_COMM)
#endif
        call pdm_tensor_sync(JOB_TENSOR_CONTRACT_SIMPLE,A,B,C)


        call tensor_buffer(nmodes2c,root=root,comm=tensor_work_comm)
        call tensor_buffer(m2cA,nmodes2c)
        call tensor_buffer(m2cB,nmodes2c)
        call tensor_buffer(order,C%mode)
        call tensor_buffer(pre1)
        call tensor_buffer(pre2)
        call tensor_buffer(sync,finalize = .true.)

#ifdef TENSORS_IN_LSDALTON
        call time_start_phase(PHASE_WORK)
#endif
     else
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase(PHASE_COMM)
#endif
        call pdm_tensor_sync(JOB_TENSOR_CONTRACT_BDENSE,A,C)

        call tensor_buffer(nmodes2c,root=root,comm=tensor_work_comm)
        call tensor_buffer(m2cA,nmodes2c)
        call tensor_buffer(m2cB,nmodes2c)
        call tensor_buffer(order,C%mode)
        call tensor_buffer(pre1)
        call tensor_buffer(pre2)
        call tensor_buffer(sync)
        call tensor_buffer(B%mode)
        call tensor_buffer(B%dims,B%mode)
        call tensor_buffer(B%elm1,B%nelms, finalize=.true.)
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase(PHASE_WORK)
#endif
     endif
  endif
#ifdef TENSORS_IN_LSDALTON
  call time_start_phase( PHASE_WORK )
#endif


  !TODO: Optimize number of buffers for A and B
  nbuffsA = 2
  if( B_dense )then
     nbuffsB = 0
  else
     nbuffsB = 2
  endif
  !Desired buffer size
  itest = long1*A%tsize*nbuffsA+long1*tsizeB*nbuffsB+long1*A%tsize+long1*tsizeB+long1*C%tsize

  !Find possible buffer size
  if(test_all_all_access.and.present(mem))then
     !use provided memory information to allcate space
     itest = max(int(mem*1024.0E0**3)/8,itest)
     if( gm_buf%init )then
        use_wrk_space = USE_GLOBAL_BUFFER
     else
        use_wrk_space = USE_INTERNAL_ALLOC
     endif
  else if(test_all_all_access.and.present(wrk).and.present(iwrk))then
     if( gm_buf%init )then
        if( gm_buf%n > iwrk .or. itest > iwrk )then
           use_wrk_space = USE_GLOBAL_BUFFER
           itest = max(gm_buf%n,itest)
        else
           use_wrk_space = USE_INPUT_WORK
           itest = iwrk
        endif
     else
        use_wrk_space = USE_INPUT_WORK
        itest = iwrk
     endif
  else
     if( gm_buf%init )then
        use_wrk_space = USE_GLOBAL_BUFFER
        itest = max(gm_buf%n,itest)
     else
        use_wrk_space = USE_INTERNAL_ALLOC
     endif
  endif

  !calculate number of buffers from possible buffer size
  nbuffsA = ((itest-(A%tsize + tsizeB + C%tsize))/ntens_to_get_from)/A%tsize
  if(B_dense)then
     nbuffsB = 0
  else
     nbuffsB = ((itest-(A%tsize + tsizeB + C%tsize))/ntens_to_get_from)/tsizeB
  endif

  if((nbuffsA==0.or.(nbuffsB==0.and..not.B_dense)).and.use_wrk_space == USE_INPUT_WORK)then
     print *,"WARNING(tensor_contract_par): the specified work space is too small, switching to allocations"
     use_wrk_space = USE_INTERNAL_ALLOC
     nbuffsA = 1
     if(.not.B_dense) nbuffsB = 1
  endif

  if(B_dense)then
     nbuffs  = max(2,nbuffsA)
     nbuffsB = 0
  else
     nbuffs  = max(2,min(nbuffsA, nbuffsB))
     nbuffsB = nbuffs
  endif
  nbuffsA = nbuffs

  select case(use_wrk_space)
  case( USE_INPUT_WORK )
     w => wrk
  case( USE_GLOBAL_BUFFER )
     call lspdm_reinit_global_buffer(itest)
     w => gm_buf%buf
  case( USE_INTERNAL_ALLOC )
     call tensor_alloc_mem(w,itest)
  case default 
     call tensor_status_quit("ERROR(tensor_contract_par): unknown buffering scheme",-1)
  end select


#ifdef VAR_PTR_RESHAPE
  buffA(1:A%tsize,1:nbuffsA) => w
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
  call c_f_pointer(c_loc(w(1)),buffA,[A%tsize,int(nbuffsA,kind=tensor_standard_int)])
#else
  call tensor_status_quit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
  if(B_dense)then
     buffB => null()
  else
#ifdef VAR_PTR_RESHAPE
     buffB(1:tsizeB,1:nbuffsB) => w(nbuffsA*A%tsize+1:)
#elif defined(COMPILER_UNDERSTANDS_FORTRAN_2003)
     call c_f_pointer(c_loc(w(nbuffsA*A%tsize+1)),buffB,[tsizeB,nbuffsB])
#else
     call tensor_status_quit("ERROR, YOUR COMPILER IS NOT F2003 COMPATIBLE",-1)
#endif
  endif

  wA => w(nbuffsA*A%tsize+nbuffsB*tsizeB+1:nbuffsA*A%tsize+nbuffsB*tsizeB+A%tsize)
  wB => w(nbuffsA*A%tsize+nbuffsB*tsizeB+A%tsize+1:nbuffsA*A%tsize+nbuffsB*tsizeB+A%tsize+tsizeB)
  wC => w(nbuffsA*A%tsize+nbuffsB*tsizeB+A%tsize+tsizeB+1:nbuffsA*A%tsize+nbuffsB*tsizeB+A%tsize+tsizeB+C%tsize)

  if( alloc_in_dummy )then

     call tensor_alloc_mem(reqA,nbuffsA)
     locA = A%lock_set(1)
     if(.not.locA)call tensor_lock_wins(A,'s',all_nodes=.true.)

     if(.not.B_dense)then
        call tensor_alloc_mem(reqB,nbuffsB)
        locB = B%lock_set(1)
        if(.not.locB)call tensor_lock_wins(B,'s',all_nodes=.true.)
     endif

  endif

  !initialize buffer book keeping
  call tensor_alloc_mem(buf_posA,nbuffsA)
  call tensor_alloc_mem(buf_logA,nbuffsA)
  call tensor_alloc_mem(nfA,     nbuffsA)

  if(.not.B_dense)then
     call tensor_alloc_mem(buf_posB,nbuffsB)
     call tensor_alloc_mem(buf_logB,nbuffsB)
     call tensor_alloc_mem(nfB,     nbuffsA)
  endif

  buf_posA = -1
  buf_logA = .false.
  call fill_buffer_for_tensor_contract_simple(0,nbuffsA,buf_logA,buf_posA,nfA,A,buffA,&
     &A,B,m2cA,m2cB,nmodes2c,C,order,.true.)
  if(.not.B_dense)then
     buf_posB = -1
     buf_logB = .false.
     call fill_buffer_for_tensor_contract_simple(0,nbuffsB,buf_logB,buf_posB,nfB,B,buffB,&
        &A,B,m2cA,m2cB,nmodes2c,C,order,.true.)
  endif

  !loop over local tiles of C and contract corresponding 
  LocalTiles: do locC = 1, C%nlti
     !get the global combined and global mode indices of the current C tile
     gc = C%ti(locC)%gt
     call get_midx(gc,gm,C%ntpm,C%mode)

     !determine gemm parameters m and n
     m_gemm = 1
     n_gemm = 1

     mA = -1
     mB = -1

     !get the uncontracted mode indices of the A and B arrays
     k = 1
     do i = 1, A%mode
        contraction_mode=.false.
        do j=1,nmodes2c
           contraction_mode = contraction_mode.or.(m2cA(j) == i)
        enddo
        if(.not.contraction_mode)then
           mA(i)   = gm(ro(k))
           ordA(k) = i
           m_gemm  = m_gemm * C%ti(locC)%d(ro(k))
           k=k+1
        endif
     enddo

     if(k-1/=A%mode-nmodes2c)then
        call tensor_status_quit("ERROR(lspdm_tensor_contract_simple): something wrong in ordering",-1)
     endif

     do i = 1,nmodes2c
        ordA(k-1+i) = m2cA(i)
        ordB(i)     = m2cB(i)
     end do

     l = 1
     do i = 1, B%mode
        contraction_mode=.false.
        do j=1,nmodes2c
           contraction_mode = contraction_mode.or.(m2cB(j) == i)
        enddo
        if(.not.contraction_mode)then
           mB(i)            = gm(ro(k))
           ordB(nmodes2c+l) = i
           n_gemm           = n_gemm * C%ti(locC)%d(ro(k))
           k=k+1
           l=l+1
        endif
     enddo

     if(B_dense)then
        do i = 1, B%mode
           tdim_ord(i)   = fBtdim(ordB(i))
        enddo
     endif

     !zero local wC and accumulate all contributions therein
#ifdef VAR_LSDEBUG
     wA = 0.0E0_tensor_dp
     wB = 0.0E0_tensor_dp
#endif
     wC = 0.0E0_tensor_dp

     !loop over all tiles in the contraction modes via a combined contraction index
     do cm = 1, cci

        sq = cm + (locC - 1) * cci

        !build full mode index for A and B
        call get_midx(cm,current_mode,max_mode_ci,nmodes2c)
        do i=1,nmodes2c
           mA(m2cA(i)) = current_mode(i)
           mB(m2cB(i)) = current_mode(i)
        enddo

        !get number of elements in tiles for A and B
        call get_tile_dim(tdimA,A,mA)
        nelmsTA = 1
        do i = 1, A%mode
           nelmsTA = nelmsTA * tdimA(i)
        end do

        if(.not.B_dense)then
           call get_tile_dim(tdimB,B,mB)
           nelmsTB = 1
           do i = 1, B%mode
              nelmsTB = nelmsTB * tdimB(i)
           end do
        endif

        !get the tiles into the local buffer, insert multiple buffering here
        cmidA = get_cidx(mA,A%ntpm,A%mode)
        cmidB = get_cidx(mB,ntpmB, B%mode)


        !get buffer positions for A and B bufs
        call find_tile_pos_in_buf(int(cmidA,kind=tensor_standard_int),buf_posA,nbuffsA,ibufA,found)
        if( .not. found) call tensor_status_quit("ERROR(lspdm_tensor_contract_simple): A tile must be present",-1)  

#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_COMM )
#endif
        call make_sure_tile_is_here(A,cmidA,nfA,buf_posA,ibufA,nbuffs)
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_WORK )
#endif

        ! sort for the contraction such that in gemm the arguments are always 'n' and 'n', 
        ! always sort such, that the contraction modes are in the order of m2CA, something smarter could be done here!!
        ! > determine the dgemm parameter k_gemm

        k_gemm = 1
        do i = 1,nmodes2c
           k_gemm = k_gemm * tdimA(m2cA(i))
        end do

        select case(A%mode)
        case(2)
           call array_reorder_2d(1.0E0_tensor_dp,buffA(:,ibufA),tdimA(1),tdimA(2),ordA,0.0E0_tensor_dp,wA)
        case(3)
           call array_reorder_3d(1.0E0_tensor_dp,buffA(:,ibufA),tdimA(1),tdimA(2),tdimA(3),ordA,0.0E0_tensor_dp,wA)
        case(4)
           call array_reorder_4d(1.0E0_tensor_dp,buffA(:,ibufA),tdimA(1),tdimA(2),tdimA(3),tdimA(4),ordA,0.0E0_tensor_dp,wA)
        case default
           call tensor_status_quit("ERROR(lspdm_tensor_contract_simple): sorting A not implemented",-1)
        end select

        buf_logA(ibufA) = .false.
        call fill_buffer_for_tensor_contract_simple(sq,nbuffsA,buf_logA,buf_posA,nfA,A,buffA,&
           &A,B,m2cA,m2cB,nmodes2c,C,order,.true.)

        if(B_dense)then
           call tile_from_fort(1.0E0_tensor_dp,B%elm1,B%dims,int(B%mode),0.0E0_tensor_dp,wB,cmidB,tdim_ord,ordB)
        else

           call find_tile_pos_in_buf(int(cmidB,kind=tensor_standard_int),buf_posB,nbuffsB,ibufB,found)

           if( .not. found) call tensor_status_quit("ERROR(lspdm_tensor_contract_simple): B tile must be present",-1)  

#ifdef TENSORS_IN_LSDALTON
           call time_start_phase( PHASE_COMM )
#endif
           call make_sure_tile_is_here(B,cmidB,nfB,buf_posB,ibufB,nbuffs)
#ifdef TENSORS_IN_LSDALTON
           call time_start_phase( PHASE_WORK )
#endif

           select case(B%mode)
           case(2)
              call array_reorder_2d(1.0E0_tensor_dp,buffB(:,ibufB),tdimB(1),tdimB(2),ordB,0.0E0_tensor_dp,wB)
           case(3)
              call array_reorder_3d(1.0E0_tensor_dp,buffB(:,ibufB),tdimB(1),tdimB(2),tdimB(3),ordB,0.0E0_tensor_dp,wB)
           case(4)
              call array_reorder_4d(1.0E0_tensor_dp,buffB(:,ibufB),tdimB(1),tdimB(2),tdimB(3),tdimB(4),ordB,&
                 &0.0E0_tensor_dp,wB)
           case default
              call tensor_status_quit("ERROR(lspdm_tensor_contract_simple): sorting B not implemented",-1)
           end select

           buf_logB(ibufB) = .false.
           call fill_buffer_for_tensor_contract_simple(sq,nbuffsB,buf_logB,buf_posB,nfB,B,buffB,&
              &A,B,m2cA,m2cB,nmodes2c,C,order,.true.)

        endif

        !carry out the contraction
        call dgemm('n','n',m_gemm,n_gemm,k_gemm,1.0E0_tensor_dp,wA,m_gemm,wB,k_gemm,1.0E0_tensor_dp,wC,m_gemm)


     end do

     call get_tile_dim(tdimC,C,gm)

     do i=1,C%mode
        tdim_product(i) = tdimC(ro(i))
     enddo

     !ADD THE FINALIZED TILE TO THE LOCAL TILE IN THE CORRECT ORDER
     select case (C%mode)
     case(2)
        call array_reorder_2d(pre1,wC,tdim_product(1),tdim_product(2),order,pre2,C%ti(locC)%t)
     case(3)
        call array_reorder_3d(pre1,wC,tdim_product(1),tdim_product(2),tdim_product(3),order,pre2,C%ti(locC)%t)
     case(4)
        call array_reorder_4d(pre1,wC,tdim_product(1),tdim_product(2),tdim_product(3),tdim_product(4),order,pre2,C%ti(locC)%t)
     case default
        call tensor_status_quit("ERROR(lspdm_tensor_contract_simple): sorting C not implemented",-1)
     end select

  enddo LocalTiles

  if(count(buf_logA) > 0 ) call tensor_status_quit("ERROR(lspdm_tensor_contract_simple): there should be no tiles left A",-1)
  call tensor_free_mem(buf_posA)
  call tensor_free_mem(buf_logA)
  call tensor_free_mem(nfA     )

  if(.not.B_dense)then
     if(count(buf_logB) > 0 ) call tensor_status_quit("ERROR(lspdm_tensor_contract_simple): there should be no tiles left B",-1)
     call tensor_free_mem(buf_posB)
     call tensor_free_mem(buf_logB)
     call tensor_free_mem(nfB)
  endif

  if(use_wrk_space==USE_INTERNAL_ALLOC)then
     call tensor_free_mem(w)
  endif
  buffA => null()
  buffB => null()
  wA    => null()
  wB    => null()
  wC    => null()


  if( alloc_in_dummy )then

     call tensor_free_mem(reqA)
     if(.not.locA)call tensor_unlock_wins(A,all_nodes=.true.)

     if(.not.B_dense)then
        call tensor_free_mem(reqB)
        if(.not.locB)call tensor_unlock_wins(B,all_nodes=.true.)
     endif

  endif

  !critical barrier if synchronization is not achieved by other measures
#ifdef TENSORS_IN_LSDALTON
  call time_start_phase( PHASE_IDLE )
#endif
  if( sync .or. tensor_always_sync ) call tensor_mpi_barrier(tensor_work_comm)
#ifdef TENSORS_IN_LSDALTON
  call time_start_phase( PHASE_WORK )
#endif
  !stop 0
#else
  call tensor_status_quit("ERROR(lspdm_tensor_contract_simple): cannot be called without MPI",-1)
#endif
  end subroutine lspdm_tensor_contract_simple

  subroutine make_sure_tile_is_here(T,cmidT,nfT,intT,ibufT,nbuffs)
     implicit none
     type(tensor), intent(in) :: T
     integer, intent(in)      :: cmidT,nbuffs,ibufT,intT(nbuffs)
     logical, intent(inout)   :: nfT(nbuffs)
     integer :: i, flushed_node, another_node, ne2, dpos, didx, dwidx

#ifdef VAR_MPI
     if( alloc_in_dummy )then

        !call tensor_mpi_wait(reqA(ibufA))
        if( nfT(ibufT) )then

           call tensor_flush_win(T, gtidx=cmidT, only_owner=.true., local = .true.)

           call get_residence_of_tile(T, intT(ibufT), flushed_node, dpos, didx, dwidx)

           do i = 1, nbuffs

              if( nfT(i) )then

                 call get_residence_of_tile(T, intT(i), another_node, dpos, didx, dwidx)

                 if(another_node == flushed_node)then

                    call get_tile_dim(ne2,T,intT(i))

                    nfT(i) = .false.

                 endif

              endif
           enddo

        endif

     else

        if(T%lock_set(cmidT))call tensor_unlock_win(T,cmidT)

     endif
#endif

  end subroutine make_sure_tile_is_here

  subroutine fill_buffer_for_tensor_contract_simple(current,nbuffs,loglist,intlist,&
        &nf,T,buf,A,B,m2cA,m2cB,nm2c,C,order,load_it)
     implicit none
     type(tensor), intent(in) :: T
     type(tensor), intent(in) :: A,B,C
     integer, intent(in)         :: current,nbuffs,order(C%mode),nm2c
     real(tensor_dp), intent(inout)  :: buf(:,:)
     integer, intent(in)         :: m2cA(nm2c),m2cB(nm2c)
     logical, intent(inout)      :: loglist(nbuffs), nf(nbuffs)
     integer, intent(inout)      :: intlist(nbuffs)
     logical, intent(in)         :: load_it
     integer, parameter :: one = 1
     integer, parameter :: two = 2
     integer :: next, nn, cc, gc, cci, tt(two), mm(two), ro(C%mode), gm(C%mode)
     integer :: maxCI(nm2c), cm(nm2c)
     integer :: mA(A%mode)
     integer :: mB(B%mode)
     integer :: mT(max(A%mode,B%mode))
     integer :: i,j,k, cmidA, cmidB, cmidT, ibufT, nelmsT
     logical :: contM, TisA, TisB, get_new, found, break_condition, BDense
     integer :: maxcntr
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: nnod
     integer(kind=tensor_mpi_kind) :: me
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)
     call tensor_get_size(nnod)

#ifdef VAR_MPI
     ! find the identity of T
     TisA   = ( A%addr_p_arr(me+1) == T%addr_p_arr(me+1) )
     TisB   = ( B%addr_p_arr(me+1) == T%addr_p_arr(me+1) )
     BDense = (B%itype == TT_DENSE) .or. (B%itype == TT_REPLICATED)

     if( TisA .eqv. TisB )then
        call tensor_status_quit("ERROR(fill_buffer_for_tensor_contract_simple): both cannot be T, one has to be T",-1)
     endif

     if( TisB .and. BDense )then
        call tensor_status_quit("ERROR(fill_buffer_for_tensor_contract_simple): T must not be B if B is dense",-1)
     endif

     do i = 1, C%mode
        ro(order(i)) = i
     enddo

     cci = 1
     do i=1,nm2c
        cci      = cci * A%ntpm(m2cA(i))
        maxCI(i) = A%ntpm(m2cA(i))
     enddo

     maxcntr         = cci + (C%nlti-1)*cci
     break_condition = ( count(loglist) < nbuffs ) .and. (current < maxcntr)
     next            = current + 1
     mm              = [cci,int(C%nlti)]

     do while (break_condition)

        call get_midx(next,tt,mm,two)
        cc = tt(one)
        nn = tt(two)

        gc = C%ti(nn)%gt

        call get_midx(gc,gm,C%ntpm,C%mode)

        mA = -1
        mB = -1

        !get the uncontracted mode indices of the A and B arrays
        k = 1
        do i = 1, A%mode
           contM=.false.
           do j=1,nm2c
              contM = contM.or.(m2cA(j) == i)
           enddo
           if(.not.contM)then
              mA(i) = gm(ro(k))
              k=k+1
           endif
        enddo

        if(k-1/=A%mode-nm2c)then
           call tensor_status_quit("ERROR(fill_buffer_for_tensor_contract_simple): something wrong in ordering",-1)
        endif

        do i = 1, B%mode
           contM=.false.
           do j=1,nm2c
              contM = contM.or.(m2cB(j) == i)
           enddo
           if(.not.contM)then
              mB(i) = gm(ro(k))
              k=k+1
           endif
        enddo

        call get_midx(cc,cm,maxCI,nm2c)

        do i=1,nm2c
           mA(m2cA(i)) = cm(i)
           mB(m2cB(i)) = cm(i)
        enddo


        cmidA = get_cidx(mA,A%ntpm,A%mode)

        if(.not. BDense) then
           cmidB = get_cidx(mB,B%ntpm,B%mode)
        endif


        if( TisA ) cmidT = cmidA
        if( TisB ) cmidT = cmidB


        call check_if_new_instance_needed(cmidT,intlist,nbuffs,get_new,pos=ibufT,set_needed=loglist)

        if( get_new .and. load_it )then

           call find_free_pos_in_buf(loglist,nbuffs,ibufT,found)

           if( found )then

              call get_tile_dim(nelmsT,T,cmidT)

#ifdef TENSORS_IN_LSDALTON
              call time_start_phase( PHASE_COMM )
#endif

              if( alloc_in_dummy )then

                 !call tensor_get_tile(T,cmidT,buf(:,ibufT),nelmsTA,lock_set=.true.,req=reqA(ibufA))
                 call tensor_get_tile(T,cmidT,buf(:,ibufT),nelmsT,lock_set=.true.)

              else
                 call tensor_lock_win(T,cmidT,'s')
                 call tensor_get_tile(T,cmidT,buf(:,ibufT),nelmsT,lock_set=.true.)
              endif

#ifdef TENSORS_IN_LSDALTON
              call time_start_phase( PHASE_WORK )
#endif

              intlist(ibufT) = cmidT
              loglist(ibufT) = .true.
              nf(ibufT)      = .true.

           endif

        endif

        next = next + 1
        break_condition = ( (count(loglist) < nbuffs) .and. ( next <= maxcntr ) )

     enddo
#endif

  end subroutine fill_buffer_for_tensor_contract_simple

  !> \brief add tiled distributed data to a basic fortran type array
  !> \author Patrick Ettenhuber
  !> date march 2013
  subroutine add_tileddata2fort(arr,b,fort,nelms,pdm,order)
     implicit none
     !> array to add to the input
     type(tensor),intent(in) :: arr
     !> nuber of elements in the array
     integer(kind=tensor_long_int), intent(in) :: nelms
     !> basic fotran type array to which arr is added
     real(tensor_dp),intent(inout) :: fort(nelms)
     !> scaling factor for arr
     real(tensor_dp),intent(in) :: b
     !> logical specifying whether the tiles are in pdm
     logical, intent(in) :: pdm
     !> reorder if the array is reorder with respect to the fortran array
     integer, intent(in), optional :: order(arr%mode)
     integer :: i,j,k,tmdidx(arr%mode), o(arr%mode), mode
     integer :: l,nelintile,tdim(arr%mode),fullfortdim(arr%mode)
     real(tensor_dp), pointer :: tmp(:)
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif

     tdim = arr%tdim
     mode = arr%mode

     !check nelms
     if(nelms/=arr%nelms)call tensor_status_quit("ERROR(cp_tileddate2fort):array&
        &dimensions are not the same",lspdm_stdout)

     !allocate space 
     if(pdm)then
        call tensor_alloc_mem(tmp,arr%tsize)
     endif

     do i=1,arr%mode
        o(i)=i
     enddo

     if(present(order))o=order

     do i = 1, arr%mode
        fullfortdim(i) = arr%dims(o(i))
     enddo

     !TODO: use async buffering

     do i=1,arr%ntiles
        call get_midx(i,tmdidx,arr%ntpm,arr%mode)
        call get_tile_dim(nelintile,arr,i)
        if(pdm)then
#ifdef VAR_MPI
#ifdef TENSORS_IN_LSDALTON
           call time_start_phase( PHASE_COMM )
#endif
           call tensor_get_tile(arr,int(i,kind=tensor_standard_int),tmp,nelintile)
#ifdef TENSORS_IN_LSDALTON
           call time_start_phase( PHASE_WORK )
#endif
#endif
        else
           tmp => arr%ti(i)%t
        endif
        call tile_in_fort(b,tmp,i,tdim,1.0E0_tensor_dp,fort,fullfortdim,mode,o)
     enddo

     if(pdm)then
        call tensor_free_mem(tmp)
     else
        nullify(tmp)
     endif

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
  end subroutine add_tileddata2fort

  subroutine cp_tileddata2fort(arr,fort,nelms,pdm,order)
     implicit none
     type(tensor),intent(in) :: arr
     integer(kind=tensor_long_int), intent(in) :: nelms
     real(tensor_dp),intent(inout) :: fort(nelms)
     logical, intent(in) :: pdm
     integer, intent(in), optional :: order(arr%mode)
     integer :: i,j,k,minimode(arr%mode),o(arr%mode), mode
     integer :: glbmodeidx(arr%mode),glbidx,nelintile,tdim(arr%mode),fullfortdim(arr%mode)
     real(tensor_dp), pointer :: tmp(:)

     tdim = arr%tdim
     mode = arr%mode


     if(nelms/=arr%nelms)call tensor_status_quit("ERROR(cp_tileddate2fort):array&
        &dimensions are not the same",lspdm_stdout)
     if(pdm)then
        call tensor_alloc_mem(tmp,arr%tsize)
     endif
     do i=1,arr%mode
        o(i)=i
     enddo
     if(present(order))o=order

     do i = 1, arr%mode
        fullfortdim(i) = arr%dims(o(i))
     enddo

     do i=1,arr%ntiles
        !call get_midx(i,tmdidx,arr%ntpm,arr%mode)
        call get_tile_dim(nelintile,arr,i)
        if(pdm)then
#ifdef VAR_MPI
           call tensor_get_tile(arr,int(i,kind=tensor_standard_int),tmp,nelintile)
#endif
        else
           tmp => arr%ti(i)%t
        endif
        call tile_in_fort(1.0E0_tensor_dp,tmp,i,tdim,0.0E0_tensor_dp,fort,fullfortdim,mode,o)
     enddo

     if(pdm)then
        call tensor_free_mem(tmp)
     else
        nullify(tmp)
     endif
  end subroutine cp_tileddata2fort




  ! arr = pre1 * fort + pre2 * arr
  subroutine tensor_scatter(pre1,fort,pre2,arr,nelms,oo,wrk,iwrk)
     implicit none
     real(tensor_dp),intent(in)          :: pre1,pre2
     type(tensor),intent(in)         :: arr
     integer(kind=tensor_long_int), intent(in)  :: nelms
     real(tensor_dp),intent(in)          :: fort(nelms)
     integer(kind=tensor_mpi_kind)           :: nod
     integer, intent(in), optional             :: oo(arr%mode)
     real(tensor_dp),intent(inout),target,optional :: wrk(*)
     integer(kind=tensor_long_int),intent(in),optional,target:: iwrk
     integer(kind=tensor_mpi_kind) :: src
     integer :: mode, tdim(arr%mode)
     integer               :: i,ltidx,o(arr%mode)
     integer               :: fullfortdim(arr%mode)
     real(tensor_dp), pointer  :: tmp(:)
     integer               :: tmps, elms_sent,last_flush_i,j
     logical               :: internal_alloc,lock_outside,ls
     integer               :: maxintmp,b,e,minstart,bidx
     integer(kind=tensor_mpi_kind),pointer   :: req(:)
     integer(kind=tensor_long_int) :: itest,iti, nelintile
     logical :: lock_was_not_set
#ifdef VAR_MPI
#ifdef COMPILER_UNDERSTANDS_FORTRAN_2003
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: nnod
     integer(kind=tensor_mpi_kind) :: me
     procedure(put_acc_tile), pointer :: put_acc => null()
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)
     call tensor_get_size(nnod)

     acc_tlong1 => tensor_acct88
     put_tlong1 => tensor_putt88

     mode = arr%mode
     tdim = arr%tdim

     do i=1,arr%mode
        o(i)=i
     enddo
     if(present(oo))o=oo

     if(arr%mode==4)then
        if(o(1) == 2.and.o(2)==4.and.o(3)==1.and.o(4)==3)&
           &print *,"WARNING(tensor_scatter)this sort may be implemented,&
           & wrongly, plese check your results"
     endif

     if(pre2==0.0E0_tensor_dp) put_acc => tensor_putt88
     if(pre2/=0.0E0_tensor_dp) put_acc => tensor_acct88

     if(pre2/=0.0E0_tensor_dp.and.pre2/=1.0E0_tensor_dp)then
        call tensor_scale_td(arr,pre2)
     endif

#ifdef VAR_LSDEBUG
     if((present(wrk).and..not.present(iwrk)).or.(.not.present(wrk).and.present(iwrk)))then
        call tensor_status_quit("ERROR(tensor_scatter):both or neither wrk and iwrk have to &
           &be given",-1)
     endif
#endif

     !CHECK WHICH BUFFERING METHOD TO USE
     internal_alloc = .not. gm_buf%init
     itest          = 0
     if(present(wrk).and.present(iwrk))then
        if( iwrk > arr%tsize )then
           internal_alloc = .false.
        endif
        itest = iwrk
     endif

     if(internal_alloc)then
#ifdef VAR_LSDEBUG
        print *,'WARNING(tensor_scatter):Allocating internally'
#endif
        tmps = arr%tsize
        call tensor_alloc_mem(tmp,tmps)
     else
        if( itest > gm_buf%n )then
           tmps =  itest
           tmp  => wrk(1:tmps)
        else
           tmps =  gm_buf%n
           tmp  => gm_buf%buf(1:tmps)
        endif
     endif

#ifdef VAR_LSDEBUG
     if(nelms/=arr%nelms)call tensor_status_quit("ERROR(tensor_scatter):array&
        &dimensions are not the same",lspdm_stdout)
#endif

     do i = 1, arr%mode
        fullfortdim(i) = arr%dims(o(i))
     enddo

     maxintmp = tmps / arr%tsize

     call tensor_alloc_mem(req,maxintmp)

     lock_was_not_set = .not.arr%lock_set(1)
     if( alloc_in_dummy .and. lock_was_not_set)call tensor_lock_wins(arr,'s', all_nodes = .true. )
     elms_sent    = 0
     last_flush_i = 0

     do i=1,arr%ntiles

        !set the buffer index
        bidx = mod(i-1,maxintmp)+1

        if( i>maxintmp )then
           if(alloc_in_dummy)then
              call get_tile_dim(nelintile,arr,i-maxintmp)

              call tensor_mpi_wait(req(bidx))

           else
              if(arr%lock_set(i-maxintmp)) call tensor_unlock_win(arr,i-maxintmp)
           endif
        endif

        call get_tile_dim(nelintile,arr,i)

        !ADDRESSING IN TMP BUFFER ALWAYS WITH FULL TILE SIZES
        b = 1 +( bidx - 1 ) * arr%tsize
        e = b + nelintile - 1

        if( alloc_in_dummy )then
           ls = arr%lock_set(1)
        else
           ls = arr%lock_set(i)
        endif

        call tile_from_fort(pre1,fort,fullfortdim,mode,0.0E0_tensor_dp,tmp(b:e),i,tdim,o)

        iti = i

        if( alloc_in_dummy )then
           call put_acc(arr,iti,tmp(b:e),nelintile,lock_set=ls,req=req(bidx))
        else
           call put_acc(arr,iti,tmp(b:e),nelintile,lock_set=ls)
        endif

        elms_sent = elms_sent + nelintile

     enddo

     if(arr%ntiles - maxintmp >= 0)then
        minstart = arr%ntiles - maxintmp + 1
     else
        minstart = 1
     endif

     if( alloc_in_dummy )then
        do i=minstart, arr%ntiles
           call get_tile_dim(nelintile,arr,i)
           bidx = mod(i-1,maxintmp)+1
           call tensor_mpi_wait(req(bidx))
        enddo
        if(arr%lock_set(1).and.lock_was_not_set)call tensor_unlock_wins(arr, all_nodes = .true.)
     else
        do i=minstart, arr%ntiles
           if(arr%lock_set(i))call tensor_unlock_win(arr,i)
        enddo
     endif

     if(internal_alloc)then
        if(.not.lock_was_not_set.and.alloc_in_dummy)call tensor_mpi_win_flush(arr%wi(1), local = .true.)
        call tensor_free_mem(tmp)
     else
        tmp  => null()
     endif

     call tensor_free_mem(req)
#else
     call tensor_status_quit("ERROR(tensor_scatter):this routine is FORTRAN 2003 only",-1)
#endif
#else
     call tensor_status_quit("ERROR(tensor_scatter):this routine is MPI only",-1)
#endif
  end subroutine tensor_scatter


  ! fort = pre1 * arr + pre2 *fort
  subroutine tensor_gather(pre1,arr,pre2,fort,nelms,oo,wrk,iwrk)
     implicit none
     real(tensor_dp),intent(in)             :: pre1,pre2
     type(tensor),intent(in)            :: arr
     integer(kind=tensor_long_int), intent(in)     :: nelms
     real(tensor_dp),intent(inout)          :: fort(nelms)
     integer(kind=tensor_mpi_kind)              :: nod
     integer, intent(in), optional             :: oo(arr%mode)
     real(tensor_dp),intent(inout),target,optional :: wrk(*)
     integer(kind=tensor_long_int),intent(in),optional,target:: iwrk
     integer :: tdim(arr%mode), mode
     integer               :: i,ltidx,o(arr%mode)
     integer               :: nelintile,fullfortdim(arr%mode)
     real(tensor_dp), pointer  :: tmp(:)
     integer               :: tmps, elms_sent,last_flush_i,j,bidx, itest
     logical               :: internal_alloc,lock_outside,so,consecutive,ff,ls,lock_was_not_set
     integer               :: maxintmp,b,e,minstart
#ifdef VAR_MPI
     integer(kind=tensor_mpi_kind),pointer :: req(:)
     real(tensor_dp)   :: nrm

     tdim = arr%tdim
     mode = arr%mode

     do i=1,arr%mode
        o(i)=i
     enddo
     if(present(oo))o=oo

     so = .true.
     do i=1,arr%mode
        if(o(i)/=i)so = .false.
     enddo

#ifdef VAR_LSDEBUG
     if((present(wrk).and..not.present(iwrk)).or.(.not.present(wrk).and.present(iwrk)))then
        call tensor_status_quit('ERROR(tensor_gather):both or neither wrk and iwrk have to &
           &be given',-1)
     endif
#endif


     !CHECK WHICH BUFFERING METHOD TO USE
     internal_alloc = .not. gm_buf%init
     itest          = 0
     if(present(wrk).and.present(iwrk))then
        if( iwrk > arr%tsize )then
           internal_alloc = .false.
        endif
        itest = iwrk
     endif

     if(.not. internal_alloc)then
        if( itest > gm_buf%n )then
           tmps =  itest
           tmp  => wrk(1:tmps)
        else
           tmps =  gm_buf%n
           tmp  => gm_buf%buf(1:tmps)
        endif
     endif

#ifdef VAR_LSDEBUG
     if(nelms/=arr%nelms)call tensor_status_quit('ERROR(tensor_gather):array&
        &dimensions are not the same',lspdm_stdout)
#endif

     do i = 1, arr%mode
        fullfortdim(i) = arr%dims(o(i))
     enddo

     consecutive = .true.
     ff = .false.
     do i = 1, arr%mode
        fullfortdim(i) = arr%dims(o(i))
        if( arr%dims(i) /= arr%tdim(i) .and. .not. ff) ff = .true.
        if( arr%dims(i) /= arr%tdim(i) .and. arr%tdim(i) /= 1 .and. ff) consecutive = .false.
     enddo

     elms_sent    = 0
     last_flush_i = 0

     lock_was_not_set = .not.arr%lock_set(1)
     if( lock_was_not_set .or. .not. alloc_in_dummy )call tensor_lock_wins(arr,'s', all_nodes = alloc_in_dummy,check =.true. )


     if(so.and.pre1==1.0E0_tensor_dp.and.pre2==0.0E0_tensor_dp.and.consecutive)then

        b=1
        do i=1,arr%ntiles

           call get_tile_dim(nelintile,arr,i)
           e = b + nelintile - 1

           if( alloc_in_dummy ) then
              ls = arr%lock_set(1)
           else
              ls = arr%lock_set(i)
           endif

           call tensor_get_tile(arr,i,fort(b:e),nelintile,lock_set=ls)
           b = e + 1
           elms_sent = elms_sent + nelintile

        enddo

     else

        if(internal_alloc)then
#ifdef VAR_LSDEBUG
           print *,'WARNING(tensor_gather):Allocating internally'
#endif
           tmps = arr%tsize
           call tensor_alloc_mem(tmp,tmps)
        endif

        maxintmp = tmps / arr%tsize

        call tensor_alloc_mem(req,maxintmp)

        do i=1,arr%ntiles

           !set the buffer index
           bidx = mod(i-1,maxintmp)+1

           call get_tile_dim(nelintile,arr,i)

           !ADDRESSING IN TMP BUFFER ALWAYS WITH FULL TILE SIZES
           b = 1 + (bidx - 1) * arr%tsize
           !b = 1 + mod(i - 1, maxintmp) * arr%tsize
           e = b + arr%tsize - 1


           if(i>maxintmp)then

              call get_tile_dim(nelintile,arr,i-maxintmp)

              if( alloc_in_dummy ) then
                 call tensor_mpi_wait(req(bidx))
              else
                 if(.not.arr%lock_set(i-maxintmp)) call tensor_status_quit("ERROR(tensor_gather): lock has to be set here 1!",-1)
                 call tensor_unlock_win(arr,i-maxintmp)
              endif
              !call pn(tmp(b:e),nelintile,norm=nrm)
              !write(*,'("have tile",I9," :",g9.3)')i-maxintmp,nrm
              call tile_in_fort(pre1,tmp(b:e),i-maxintmp,tdim,pre2,fort,fullfortdim,mode,o)
           endif

           call get_tile_dim(nelintile,arr,i)

           if( alloc_in_dummy) then
              call tensor_get_tile(arr,int(i,kind=tensor_standard_int),tmp(b:e),nelintile,lock_set=.true.,req = req(bidx))
           else
              call tensor_get_tile(arr,int(i,kind=tensor_standard_int),tmp(b:e),nelintile,lock_set=.true.)
           endif

           elms_sent = elms_sent + nelintile

        enddo

        minstart = max(arr%ntiles - maxintmp + 1,1)

        do i=minstart, arr%ntiles

           bidx = mod(i-1,maxintmp)+1

           b = 1 + (bidx - 1) * arr%tsize
           e = b + arr%tsize -1

           call get_tile_dim(nelintile,arr,i)

           if( alloc_in_dummy )then
              call tensor_mpi_wait(req(bidx))
           else
              if(.not.arr%lock_set(i)) call tensor_status_quit("ERROR(tensor_gather): lock has to be set here 2!",-1)
              call tensor_unlock_win(arr,i)
           endif

           !call pn(tmp(b:e),nelintile,norm=nrm)
           call tile_in_fort(pre1,tmp(b:e),i,tdim,pre2,fort,fullfortdim,mode,o)

        enddo

        if(internal_alloc)then
           call tensor_free_mem(tmp)
        else
           tmp  => null()
        endif
        call tensor_free_mem(req)
     endif
     if( alloc_in_dummy .and. lock_was_not_set )call tensor_unlock_wins(arr, all_nodes = .true. )
     if( .not. alloc_in_dummy )then
        do i = 1, arr%nwins
           if(arr%lock_set(i))call tensor_status_quit("ERROR(tensor_gather) a window is not freed",-1)
        enddo
     endif
#else
     call tensor_status_quit('ERROR(tensor_gather):this routine is MPI only',-1)
#endif
  end subroutine tensor_gather

  !gather as 2 coulomb minus exchange
  subroutine tensor_gather_2cme(arr,fort,nelms,pos,oo,wrk,iwrk)
     implicit none
     type(tensor),intent(in)             :: arr
     integer(kind=tensor_long_int), intent(in)     :: nelms
     real(tensor_dp),intent(inout)          :: fort(nelms)
     integer,intent(in)                 :: pos(2)
     integer(kind=tensor_mpi_kind)              :: nod
     integer, intent(in), optional             :: oo(arr%mode)
     real(tensor_dp),intent(inout),target,optional :: wrk(*)
     integer(kind=tensor_long_int),intent(in),optional,target:: iwrk
     integer(kind=tensor_mpi_kind) :: src
     integer :: tdim(arr%mode), mode
     integer               :: i,ltidx,o(arr%mode),excho(arr%mode)
     integer               :: nelintile,fullfortdim(arr%mode)
     real(tensor_dp), pointer  :: tmp(:)
     integer               :: tmps, elms_sent,last_flush_i,j
     logical               :: internal_alloc,lock_outside,ls
     integer               :: maxintmp,b,e,minstart
     real(tensor_dp)           :: pre1,pre2
#ifdef VAR_MPI
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: nnod
     integer(kind=tensor_mpi_kind) :: me
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)
     call tensor_get_size(nnod)

     tdim = arr%tdim
     mode = arr%mode

     do i=1,arr%mode
        o(i)=i
     enddo
     if(present(oo))o=oo
     excho = o
     excho(pos(2)) = o(pos(1))
     excho(pos(1)) = o(pos(2))

!#ifdef VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN
!     call assign_in_subblocks(fort(1:nelms),'=',fort(1:nelms),nelms,scal2=0.0E0_tensor_dp)
!#else
     fort(1:nelms) = 0.0E0_tensor_dp
!#endif

#ifdef VAR_LSDEBUG
     if((present(wrk).and..not.present(iwrk)).or.(.not.present(wrk).and.present(iwrk)))then
        call tensor_status_quit('ERROR(tensor_gather):both or neither wrk and iwrk have to &
           &be given',-1)
     endif
#endif


     !CHECK IF INTERNAL MEMORY ALLOCATION IS NEEDED
     internal_alloc = .true.
     if(present(wrk).and.present(iwrk))then
        if(iwrk>arr%tsize)then
           internal_alloc=.false.
#ifdef VAR_LSDEBUG
        else
           print *,'WARNING(tensor_gather):allocating internally, given buffer not large enough'
#endif
        endif
     endif

#ifdef VAR_LSDEBUG
     if(nelms/=arr%nelms)call tensor_status_quit('ERROR(tensor_gather):array&
        &dimensions are not the same',lspdm_stdout)
#endif

     do i = 1, arr%mode
        fullfortdim(i) = arr%dims(o(i))
     enddo


     elms_sent    = 0
     last_flush_i = 0


     if(internal_alloc)then
#ifdef VAR_LSDEBUG
        print *,'WARNING(tensor_gather):Allocating internally'
#endif
        tmps = arr%tsize
        call tensor_alloc_mem(tmp,tmps)
     else
        tmps =  iwrk
        tmp  => wrk(1:tmps)
     endif

     maxintmp = tmps / arr%tsize

     do i=1,arr%ntiles

        if(i>maxintmp)then
           b = 1 + mod(i - maxintmp - 1, maxintmp) * arr%tsize
           e = b + arr%tsize -1
           if( alloc_in_dummy ) then
              call tensor_flush_win(arr, gtidx = i-maxintmp, local = .false., only_owner=.true.)
           else
              if(arr%lock_set(i-maxintmp))call tensor_unlock_win(arr,i-maxintmp)
           endif
           call tile_in_fort(2.0E0_tensor_dp,tmp(b:e),i-maxintmp,tdim,1.0E0_tensor_dp,fort,fullfortdim,mode,o)
           call tile_in_fort(-1.0E0_tensor_dp,tmp(b:e),i-maxintmp,tdim,1.0E0_tensor_dp,fort,fullfortdim,mode,excho)
        endif

        call get_tile_dim(nelintile,arr,i)

        !ADDRESSING IN TMP BUFFER ALWAYS WITH FULL TILE SIZES
        b = 1 + mod(i - 1, maxintmp) * arr%tsize
        e = b + arr%tsize - 1

        if( alloc_in_dummy ) then
           ls = arr%lock_set(1)
        else
           ls = arr%lock_set(i)
        endif

        call tensor_get_tile(arr,int(i,kind=tensor_standard_int),tmp(b:e),nelintile,lock_set=ls)

     enddo

     if(arr%ntiles - maxintmp >= 0)then
        minstart = arr%ntiles - maxintmp + 1
     else
        minstart = 1
     endif

     do i=minstart, arr%ntiles
        b = 1 + mod(i - 1, maxintmp) * arr%tsize
        e = b + arr%tsize -1
        if( alloc_in_dummy.and.arr%lock_set(1))call tensor_flush_win(arr, gtidx = i, local = .false., only_owner=.false.)
        if(.not.alloc_in_dummy.and.arr%lock_set(i))call tensor_unlock_win(arr,i)
        call tile_in_fort(2.0E0_tensor_dp,tmp(b:e),i,tdim,1.0E0_tensor_dp,fort,fullfortdim,mode,o)
        call tile_in_fort(-1.0E0_tensor_dp,tmp(b:e),i,tdim,1.0E0_tensor_dp,fort,fullfortdim,mode,excho)
     enddo

     if(internal_alloc)then
        call tensor_free_mem(tmp)
     else
        tmp  => null()
     endif

#else
     call tensor_status_quit('ERROR(tensor_gather):this routine is MPI only',-1)
#endif
  end subroutine tensor_gather_2cme



  subroutine print_mem_per_node(output,allaccs,infoonmaster,narr_allocd)
     implicit none
     integer, intent(in) :: output
     logical,intent(in)  :: allaccs
     integer(kind=tensor_long_int), parameter :: nmeminfo = 9
     integer(kind=tensor_long_int),optional  :: infoonmaster(:)
     integer(kind=tensor_long_int),optional  :: narr_allocd(:)
     integer(kind=tensor_long_int) :: get_mem(nmeminfo)
     integer(kind=tensor_mpi_kind) :: i
     integer :: allallocd
     logical :: master
     real(tensor_dp) :: mb_acc,mb_put,mb_get,total(9),speed_acc,speed_get,speed_put
#ifdef VAR_MPI
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: nnod
     integer(kind=tensor_mpi_kind) :: me
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)
     call tensor_get_size(nnod)

     master = (me==root)

     if(present(infoonmaster))then
        if(size(infoonmaster) < nmeminfo*nnod)then
           call tensor_status_quit("ERROR(print_mem_per_node): infoonmaster too small",83)
        endif
     endif
     if(present(narr_allocd))then
        if(size(narr_allocd) < nnod)then
           call tensor_status_quit("ERROR(print_mem_per_node): narr_allocd too small",83)
        endif
     endif
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !NODE SPECIFIC ONE-SIDED INFO!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     mb_acc = (bytes_transferred_acc*1.0E0_tensor_dp)/(1024.0**2)
     mb_put = (bytes_transferred_put*1.0E0_tensor_dp)/(1024.0**2) 
     mb_get = (bytes_transferred_get*1.0E0_tensor_dp)/(1024.0**2)
     speed_acc=0.0E0_tensor_dp
     speed_put=0.0E0_tensor_dp
     speed_get=0.0E0_tensor_dp
     if(time_pdm_acc/=0.0E0_tensor_dp)speed_acc=mb_acc/time_pdm_acc
     if(time_pdm_put/=0.0E0_tensor_dp)speed_put=mb_put/time_pdm_put
     if(time_pdm_get/=0.0E0_tensor_dp)speed_get=mb_get/time_pdm_get

     if(.not.present(infoonmaster))then
        if(master.and..not.allaccs)then
           call pdm_tensor_sync(JOB_PRINT_MEM_INFO1)
        endif
        do i=1,nnod
           if(me+1==i)then
              write(*,'("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")')
              write(*,'("Printing memory information for rank",I3)') me
              write(*,'("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")')
              call tensor_print_memory_currents(output)
              write(*,'("")')
              write(*,'(" Printing one-sided transfer information for rank",I3)') me
              write(*,'(" ***************************************************")')
              write(*,'(I6," acc: ",f15.4," MB in ",f15.4," s, bandwidth ",f15.4," MB/s")') &
                 &nmsg_acc,mb_acc,time_pdm_acc,speed_acc
              write(*,'(I6," put: ",f15.4," MB in ",f15.4," s, bandwidth ",f15.4," MB/s")') &
                 &nmsg_put,mb_put,time_pdm_put,speed_put
              write(*,'(I6," get: ",f15.4," MB in ",f15.4," s, bandwidth ",f15.4," MB/s")') &
                 &nmsg_get,mb_get,time_pdm_get,speed_get
              write(*,'(" currently",I4," arrays allocated")')p_arr%arrays_in_use
              write(*,'("")')
           endif
           call tensor_mpi_barrier(tensor_work_comm)
        enddo
     else
        if(master.and..not.allaccs)then
           call pdm_tensor_sync(JOB_PRINT_MEM_INFO2)
        endif

        call tensor_print_memory_currents(output,get_mem)

        infoonmaster(me*nmeminfo+1:me*nmeminfo+nmeminfo)=get_mem(1:nmeminfo)

        if(present(narr_allocd))then
           narr_allocd(me+1) = p_arr%arrays_in_use
           call tensor_mpi_allreduce(narr_allocd,nnod,tensor_work_comm)
        endif
        call tensor_mpi_allreduce(infoonmaster,nmeminfo*nnod,tensor_work_comm)

     endif


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !TOTAL ONE-SIDED INFO!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     total(1) = mb_acc
     total(2) = mb_put
     total(3) = mb_get
     total(4) = time_pdm_acc
     total(5) = time_pdm_put
     total(6) = time_pdm_get
     total(7) = nmsg_acc*1.0E0_tensor_dp
     total(8) = nmsg_put*1.0E0_tensor_dp
     total(9) = nmsg_get*1.0E0_tensor_dp
     call tensor_mpi_reduce(total,9,root,tensor_work_comm)
     if(master)then
        speed_acc=0.0E0_tensor_dp
        speed_put=0.0E0_tensor_dp
        speed_get=0.0E0_tensor_dp
        if(total(4)/=0.0E0_tensor_dp)speed_acc=total(1)/total(4)
        if(total(5)/=0.0E0_tensor_dp)speed_put=total(2)/total(5)
        if(total(6)/=0.0E0_tensor_dp)speed_get=total(3)/total(6)
        if(.not.present(infoonmaster))then
           write(*,'("")')
           write(*,'("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")')
           write(*,'("Printing one-sided transfer information summary")')
           write(*,'("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")')
           write(*,'(I9," Total acc: ",f15.4," MB in ",f15.4," s, bandwidth ",f15.8," MB/s")') &
              &int(total(7)),total(1),total(4),speed_acc
           write(*,'(I9," Total put: ",f15.4," MB in ",f15.4," s, bandwidth ",f15.8," MB/s")') &
              &int(total(8)),total(2),total(5),speed_put
           write(*,'(I9," Total get: ",f15.4," MB in ",f15.4," s, bandwidth ",f15.8," MB/s")') &
              &int(total(9)),total(3),total(6),speed_get
        endif
     endif
#endif
  end subroutine print_mem_per_node


  subroutine cp_data2tiled_lowmem(arr,A,dims,mode)
     implicit none
     type(tensor),intent(inout) :: arr
     real(tensor_dp),intent(in) :: A(*)
     integer,intent(in) :: mode, dims(mode)
     integer :: fib,lt,ce,j,step,mod_step,iter,nccblocks,st
     integer(kind=tensor_mpi_kind) :: dest, assert,ierr, act_step
     integer :: loc_ti,i,nelms,fe_in_block
     integer, pointer :: elm_in_tile(:),in_tile_mode(:),orig_addr(:),remote_td(:)
     integer(kind=tensor_long_int) :: widx, comp_ti, comp_el, dpos,didx,dwidx
     logical :: pdm
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: nnod
     integer(kind=tensor_mpi_kind) :: me
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)
     call tensor_get_size(nnod)

     pdm=(arr%itype==TT_TILED_DIST)

     !TRY TO INCLUDE MPI_PUT FOR THAT OPERATION,SO MAYBE MASTER_SLAVE DEPENDENCE
     !IN THIS ROUTINE IS GONE

     !begin with sanity checks
     if(arr%mode/=mode)then
        print *,"ERROR(cp_data2tiled_lowmem):mode of array does not match mode of tiled_array"
        stop 1
     endif
     do i=1,mode
        if(arr%dims(i)/=dims(i))then
           print *,"ERROR(cp_data2tiled_lowmem):dims in input do not match dims of tiled_array"
           stop 1
        endif
     enddo
     ! corresponding elements
     call tensor_alloc_mem(elm_in_tile,arr%mode)
     call tensor_alloc_mem(in_tile_mode,arr%mode)
     call tensor_alloc_mem(orig_addr,arr%mode)
     call tensor_alloc_mem(remote_td,arr%mode)

     !find consecutive elements in both full and tiled matrix according to tile
     !dimensions --> determines largest memory blocks that can be transferred in
     !one go --> step
     step=arr%tdim(1)
     do i=1,arr%mode
        if(arr%tdim(i)==arr%dims(i).and.i<arr%mode)then
           step=step*arr%tdim(i+1)
        else
           exit
        endif
     enddo
     !determine the truncated block size --> mod_step is first the number of
     !elements in dimensions 1..i, this is used to determine the number of
     !iterations needed with the determined step size, and the mod gives the
     !remainder 
     mod_step=1
     do j=1,i
        mod_step = mod_step * arr%dims(j)
     enddo
     ! determine how many blocks, including truncated blocks fit into the
     ! (eventually comnbined) dimensions of the full matrix for if none
     ! of the following modes are considered
     iter = mod_step/step 
     mod_step = mod(mod_step,step)
     if(mod_step>0)iter=iter+1
     if(mod_step==0)mod_step=step
     ! determine how many consecutive blocks are in the full matrix, i.e.
     ! multiplying iter=number of blocks in the (combined) dimension by the
     ! number of elements in the non-consecutive modes of the array
     nccblocks = iter
     do j=i+1,arr%mode
        nccblocks = nccblocks * arr%dims(j)
     enddo
     !print *,step,mod_step,iter,nccblocks,arr%nelms

     fe_in_block=1
     do i=1,nccblocks
        if(mod(i,iter)>0) act_step=step
        if(mod(i,iter)==0)act_step=mod_step
        ! get the position for the tile of the first element in the block, the previous
        ! treatment ensures that all the following act_step elements are
        ! consecutive in the same tile
        call get_midx(fe_in_block,orig_addr,arr%dims,arr%mode)

        do j=1,arr%mode
           in_tile_mode(j)=(orig_addr(j)-1)/arr%tdim(j) + 1
        enddo
        comp_ti=get_cidx(in_tile_mode,arr%ntpm,arr%mode)

        !check where the current tile resides and jump the following steps if not
        !master where the full matrix resides or the destination slave
        !dest = mod(comp_ti-1+arr%offset,nnod) 
        if(pdm)call get_residence_of_tile(arr,comp_ti,dest,dpos,didx,dwidx) 

        !get the dimensions of the remote tile
        call get_tile_dim(remote_td,arr,comp_ti)

        !now get position of the first element of the batch in the current tile
        do j=1,arr%mode
           elm_in_tile(j) = mod(orig_addr(j)-1,arr%tdim(j)) + 1
        enddo

        !get the one index element number for the remote tile
        comp_el=get_cidx(elm_in_tile,remote_td,arr%mode)

        !copy data to the identified places
        if(pdm)then
           if( alloc_in_dummy ) then
              widx = 1
           else
              widx = comp_ti
           endif
#ifdef VAR_MPI
#ifdef TENSORS_IN_LSDALTON
           call time_start_phase( PHASE_COMM )
#endif
           call tensor_mpi_win_lock(dest,arr%wi(comp_ti),'e')
           call tensor_mpi_put(A(fe_in_block:fe_in_block+act_step-1),act_step,comp_el,dest,arr%wi(widx))
           call tensor_mpi_win_unlock(dest,arr%wi(comp_ti))
#ifdef TENSORS_IN_LSDALTON
           call time_start_phase( PHASE_WORK )
#endif
#endif
        else
           call dcopy(act_step,A(fe_in_block),1,arr%ti(comp_ti)%t(comp_el),1)
        endif
        fe_in_block=fe_in_block + act_step
     enddo
     call tensor_free_mem(remote_td)
     call tensor_free_mem(elm_in_tile)
     call tensor_free_mem(in_tile_mode)
     call tensor_free_mem(orig_addr)

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
  end subroutine cp_data2tiled_lowmem



  !\> \brief lock all windows of a tensor from the current node with the
  !specified lock and assertion
  !\> \author Patrick Ettenhuber
  !\> \date July 2013
  subroutine tensor_lock_win8(arr,ti_idx,locktype,assert)
     implicit none
     type(tensor) :: arr
     integer(kind=tensor_long_int),intent(in) :: ti_idx
     character, intent(in) :: locktype
     integer(kind=tensor_mpi_kind), optional,intent(in) :: assert
     integer(kind=tensor_mpi_kind) ::node
     integer(kind=tensor_long_int) :: wi_idx,dpos,didx,dwidx
#ifdef VAR_MPI
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_COMM )
#endif

     call get_residence_of_tile(arr,ti_idx,node,dpos,didx,wi_idx)
     call tensor_mpi_win_lock(node,arr%wi(wi_idx),locktype,ass=assert)
     arr%lock_set(wi_idx)=.true.

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
#endif
  end subroutine tensor_lock_win8
  subroutine tensor_lock_win4(arr,ti_idx,locktype,assert)
     implicit none
     type(tensor) :: arr
     integer(kind=tensor_standard_int),intent(in) :: ti_idx
     character, intent(in) :: locktype
     integer(kind=tensor_mpi_kind), optional,intent(in) :: assert
     integer(kind=tensor_mpi_kind) ::node
     integer(kind=tensor_long_int) :: wi_idx,dpos,didx,dwidx
#ifdef VAR_MPI
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_COMM )
#endif

     call get_residence_of_tile(arr,ti_idx,node,dpos,didx,wi_idx)
     call tensor_mpi_win_lock(node,arr%wi(wi_idx),locktype,ass=assert)
     arr%lock_set(wi_idx)=.true.

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
#endif
  end subroutine tensor_lock_win4

  subroutine tensor_lock_win_on_all_nodes(arr,ti_idx,locktype,assert)
     implicit none
     type(tensor) :: arr
     integer,intent(in) :: ti_idx
     character, intent(in) :: locktype
     integer(kind=tensor_mpi_kind), optional,intent(in) :: assert
     integer(kind=tensor_mpi_kind) ::node
     integer(kind=tensor_long_int) :: wi_idx,dpos,didx,dwidx
#ifdef VAR_MPI
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_COMM )
#endif
     call get_residence_of_tile(arr,ti_idx,node,dpos,didx,wi_idx)
     call tensor_mpi_win_lock(node,arr%wi(wi_idx),locktype,ass=assert)
     arr%lock_set(wi_idx)=.true.

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
#endif
  end subroutine tensor_lock_win_on_all_nodes

  subroutine tensor_unlock_win(arr,ti_idx)
     implicit none
     type(tensor) :: arr
     integer,intent(in) :: ti_idx
     integer(kind=tensor_mpi_kind) :: node
     integer(kind=tensor_long_int) :: wi_idx,dpos,didx,dwidx
#ifdef VAR_MPI
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_COMM )
#endif


     call get_residence_of_tile(arr,ti_idx,node,dpos,didx,wi_idx)
     call tensor_mpi_win_unlock(node,arr%wi(wi_idx))
     arr%lock_set(wi_idx) = .false.

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
#endif
  end subroutine tensor_unlock_win

  subroutine tensor_lock_wins(arr,locktype,assert,all_nodes,check)
     implicit none
     type(tensor) :: arr
     character, intent(in) :: locktype
     integer(kind=tensor_mpi_kind), optional,intent(in) :: assert
     logical, optional,intent(in) :: all_nodes, check
     integer(kind=tensor_mpi_kind) :: node
     integer :: i
     logical :: an,ch
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_COMM )
#endif
#ifdef VAR_MPI
     !Per default this routine only locks the window on the specific nodes where
     !the tiles reside
     an = .false.
     if(present(all_nodes))an = all_nodes
     ch = .false.
     if(present(check)) ch = check

     if(an)then

        if(locktype == 'e') print *,"WARNING(tensor_lock_wins): you called a &
           & lock_win on all windows on all nodees with an exclusive lock, &
           &that should not be done"

        do i=1,arr%nwins
           call tensor_mpi_win_lock_all(arr%wi(i),ass=assert)
        enddo

     else
        if(ch)then
           do i=1,arr%nwins
              if( .not.arr%lock_set(i) )call tensor_lock_win(arr,i,locktype,assert=assert)
           enddo
        else
           do i=1,arr%nwins
              if( arr%lock_set(i) ) call tensor_status_quit("ERROR(tensor_lock_wins): no lock should be set here",-1)
              call tensor_lock_win(arr,i,locktype,assert=assert)
           enddo
        endif

     endif

     arr%lock_set = .true.

#endif
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
  end subroutine tensor_lock_wins

  subroutine tensor_lock_local_wins(arr,locktype,assert)
     implicit none
     type(tensor) :: arr
     character, intent(in) :: locktype
     integer(kind=tensor_mpi_kind), optional,intent(in) :: assert
     integer(kind=tensor_mpi_kind) :: node
     integer :: i,gt
#ifdef VAR_MPI
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(node)

     if( alloc_in_dummy )then

        call tensor_mpi_win_lock(node,arr%wi(1),locktype,ass=assert)
        arr%lock_set(1) = .true.

     else

        do i=1,arr%nlti
           gt = arr%ti(i)%gt
           call tensor_lock_win(arr,gt,locktype,assert=assert)
        enddo

     endif

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
#endif
  end subroutine tensor_lock_local_wins
  subroutine tensor_unlock_local_wins(arr)
     implicit none
     type(tensor) :: arr
     integer(kind=tensor_mpi_kind) :: node
     integer :: i,gt
#ifdef VAR_MPI
     call tensor_get_rank(node)
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_COMM )
#endif

     if( alloc_in_dummy )then

        call tensor_mpi_win_unlock(node,arr%wi(1))
        arr%lock_set(1) = .false.

     else

        do i=1,arr%nlti
           gt = arr%ti(i)%gt
           call tensor_unlock_win(arr,gt)
        enddo

     endif

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
#endif
  end subroutine tensor_unlock_local_wins

  !\> \brief unlock all windows of a tensor 
  !\> \author Patrick Ettenhuber
  !\> \date July 2013
  subroutine tensor_unlock_wins(arr,check,all_nodes)
     implicit none
     type(tensor) :: arr
     logical, intent(in),optional:: check
     logical, optional,intent(in) :: all_nodes
     integer(kind=tensor_mpi_kind) :: node
     integer :: i,dpos,didx,dwidx
     logical :: ch,an

#ifdef VAR_MPI
     ch = .false.
     if(present(check))     ch = check
     an = .false.
     if(present(all_nodes)) an = all_nodes

     if(ch.and.an)call tensor_status_quit("ERROR(tensor_unlock_wins): input error, an all &
        &node unlock can only be requested without check, since an MPI_UNLOCK_ALL &
        &needs a preceding MPI_LOCK_ALL",-1)

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_COMM )
#endif

     if(arr%itype/=TT_DENSE.or.arr%atype=="TDAR".or.arr%atype=="TDPD")then

        if(an)then
           !UNLOCK ALL WINDOWS ON ALL NODES
           do i=1,arr%nwins
              call tensor_mpi_win_unlock_all(arr%wi(i))
           enddo

           arr%lock_set =.false.

        else

           if(ch)then

              !UNLOCK ALL WINDOWS THAT ARE MARKED AS LOCKED
              do i=1,arr%nwins
                 if(arr%lock_set(i))then
                    call get_residence_of_tile(arr,i,node,dpos,didx,dwidx)
                    call tensor_mpi_win_unlock(node,arr%wi(i))
                    arr%lock_set(i)=.false.
                 endif
              enddo

           else

              !UNLOCK ALL WINDOWS 
              do i=1,arr%nwins
                 call get_residence_of_tile(arr,i,node,dpos,didx,dwidx)
                 call tensor_mpi_win_unlock(node,arr%wi(i))
                 arr%lock_set(i)=.false.
              enddo

           endif
        endif
     endif

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif

#endif
  end subroutine tensor_unlock_wins


  subroutine pn(a,n,norm)
     implicit none
     integer,intent(in) :: n
     real(tensor_dp), intent(in) :: a(n)
     real(tensor_dp), intent(out), optional :: norm
     integer :: i
     real(tensor_dp) :: nrm
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
     nrm = 0.0E0_tensor_dp
     do i=1,n
        nrm=nrm+a(i)*a(i)
     enddo
     nrm = sqrt(nrm)
     if(present(norm)) then
        norm = nrm
     else
        print *,"NORM:",nrm
     endif
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
  end subroutine

  subroutine tensor_deallocate_dense(arr, change)
     implicit none
     type(tensor), intent(inout) :: arr
     logical, intent(in), optional :: change
     integer :: a
     logical :: change_int
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: nnod
     integer(kind=tensor_mpi_kind) :: me
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)
     call tensor_get_size(nnod)

     call memory_deallocate_tensor_dense(arr)
     a = me + 1
     change_int = .true.
     if(present(change)) change_int = change

     if(associated(p_arr%a(arr%addr_p_arr(a))%elm1))then
        p_arr%a(arr%addr_p_arr(a))%elm1 => null()
     endif

     !this is a simple solution, but will not work if itype is changed to be a
     !pointer, which it should
     if(change_int)then
        arr%itype = p_arr%a(arr%addr_p_arr(a))%itype
     endif

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
  end subroutine tensor_deallocate_dense


  subroutine tensor_free_pdm(arr)
     implicit none
     type(tensor) :: arr
     logical     :: parent
#ifdef VAR_MPI
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: nnod
     integer(kind=tensor_mpi_kind) :: me
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)
     call tensor_get_size(nnod)

     if( arr%access_type==AT_MASTER_ACCESS .and. me==root )then

     call pdm_tensor_sync(JOB_FREE_tensor_PDM,arr)

  endif

  p_arr%free_addr_on_node(arr%local_addr)=.true.
  p_arr%arrays_in_use = p_arr%arrays_in_use - 1 
  call tensor_free_basic(p_arr%a(arr%local_addr)) 
  call tensor_reset_value_defaults(p_arr%a(arr%local_addr)) 
  call tensor_nullify_pointers(arr)
#endif
#ifdef TENSORS_IN_LSDALTON
  call time_start_phase( PHASE_WORK )
#endif
  end subroutine tensor_free_pdm

  subroutine get_distribution_info(arr,force_offset)
     implicit none
     type(tensor),intent(inout) :: arr
     integer, intent(in), optional :: force_offset
     integer :: i,ntiles2dis
     logical :: parent
     integer(kind=tensor_mpi_kind) :: me,nnod,pc_me,pc_nnod,buf(2)
#ifdef VAR_MPI
     integer(kind=tensor_mpi_kind), parameter :: root = 0
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)
     call tensor_get_size(nnod)

     if(arr%access_type==AT_NO_PDM_ACCESS)then

        arr%offset       = 0
        p_arr%new_offset = 0
        arr%nlti         = arr%ntiles

     else

        select case(arr%itype)
        case(TT_TILED, TT_TILED_REPL)

           arr%offset       = 0
           p_arr%new_offset = 0
           arr%nlti         = arr%ntiles

        case(TT_TILED_DIST)

           if(present(force_offset))then
              arr%offset       = force_offset
           else
              arr%offset       = p_arr%new_offset
              p_arr%new_offset = mod(p_arr%new_offset+arr%ntiles,nnod)
           endif

           arr%nlti         = arr%ntiles/nnod

           if(mod(arr%ntiles,nnod)>mod(me+nnod-arr%offset,nnod))then
              arr%nlti=arr%nlti+1
           endif

        case default
           call tensor_status_quit("ERROR(get_distribution_info): unknown array%itype",-1)
        end select

     end if
#endif
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
  end subroutine get_distribution_info

  !> \brief routine to get a free address in the persisten array
  !> \autor Patrick Ettenhuber
  function get_free_address(occ_addr) result(addr)
     implicit none
     !> retrurn value with the address
     integer :: addr
     !> logical which tells the routine to set the value of the found address to occupied
     logical, intent(in) :: occ_addr
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif

     if(p_arr%arrays_in_use==n_arrays)then
        call tensor_status_quit("ERROR(get_free_address):max number of arrays in p_arr allocated, change&
           & the parameter n_arrays in lsutil/lspdm_tensor_operations.F90 and recompile",-1)
     endif

     do addr=1,p_arr%arrays_in_use+1
        if(p_arr%free_addr_on_node(addr))then
           if(occ_addr)p_arr%free_addr_on_node(addr) = .false.
           return
        endif
     enddo

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
  end function get_free_address

  !> \brief debugging routine to check the norms of individual tiles
  !> \author Patrick Ettenhuber
  subroutine tensor_tiled_pdm_print_ti_nrm(arr,globtinr,whichnode,nrm) 
     implicit none
     !> input array for which to check the tile
     type(tensor), intent(in) :: arr
     !> global index number of the tile
     integer, intent(in) :: globtinr
     !> optional input, return value for the destination of the tile
     integer, intent(inout), optional :: whichnode
     !> optional input, return value for the norm
     real(tensor_dp), intent(inout), optional :: nrm
     real(tensor_dp) :: norm
     integer :: i,j,loctinr
     integer(kind=tensor_standard_int) :: gtnr
     integer(kind=tensor_mpi_kind) :: dest,dpos,didx,dwidx
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: nnod
     integer(kind=tensor_mpi_kind) :: me
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)
     call tensor_get_size(nnod)

#ifdef VAR_MPI
     gtnr=globtinr
     if(arr%access_type==AT_MASTER_ACCESS.and.me==root)then
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_COMM )
#endif
        call pdm_tensor_sync(JOB_PRINT_TI_NRM,arr)
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_WORK )
#endif
     endif
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_COMM )
#endif
     call tensor_mpi_bcast(gtnr,root,tensor_work_comm)
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif

     call get_residence_of_tile(arr,gtnr,dest,dpos,didx,dwidx)
     if(present(whichnode))whichnode=dest

     if(dest==me)then
        loctinr=(gtnr-1)/nnod + 1
        norm=0.0E0_tensor_dp
        do j=1,arr%ti(loctinr)%e
           norm = norm + arr%ti(loctinr)%t(j) * arr%ti(loctinr)%t(j)
        enddo

#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_COMM )
#endif
        call tensor_mpi_sendrecv(norm,tensor_work_comm,me,root)
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_WORK )
#endif
     endif

     if(me==root.and.me/=dest)then
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_COMM )
#endif
        call tensor_mpi_sendrecv(norm,tensor_work_comm,dest,root)
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_WORK )
#endif
     endif

     !if nrm is present return the squared norm, else print the norm
     if(me==root.and.present(nrm))then
        nrm = norm
     else if(me==root)then
        write(lspdm_stdout,'("LOCAL TILE NORM ON",I3,f20.15)') dest,sqrt(norm)
     endif
#endif
  end subroutine tensor_tiled_pdm_print_ti_nrm

  function tensor_tiled_pdm_get_nrm2(arr) result(nrm)
     implicit none
     type(tensor), intent(in) :: arr
     real(tensor_dp) :: nrm
     integer :: i,j,should
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: nnod
     integer(kind=tensor_mpi_kind) :: me
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)
     call tensor_get_size(nnod)
#ifdef VAR_MPI
     if(me==root.and.arr%access_type==AT_MASTER_ACCESS) then
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_COMM )
#endif
        call pdm_tensor_sync(JOB_GET_NRM2_TILED,arr)
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_WORK )
#endif
     endif
     nrm=0.0E0_tensor_dp

     do i=1,arr%nlti
        do j=1,arr%ti(i)%e
           nrm = nrm +(arr%ti(i)%t(j) * arr%ti(i)%t(j))
           !nrm = arr%ti(i)%t(j) 
        enddo
     enddo

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_COMM )
#endif
     if(arr%access_type==AT_MASTER_ACCESS) call tensor_mpi_reduce(nrm,root,tensor_work_comm)
     if(arr%access_type==AT_ALL_ACCESS)    call tensor_mpi_allreduce(nrm,tensor_work_comm)
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif

#else
     nrm = 0.0E0_tensor_dp
#endif
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
  end function tensor_tiled_pdm_get_nrm2

  subroutine change_access_type_td(arr,totype)
     implicit none
     type(tensor),intent(inout) :: arr
     integer,intent(in) :: totype
#ifdef VAR_MPI
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: nnod
     integer(kind=tensor_mpi_kind) :: me
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)
     call tensor_get_size(nnod)

     if(totype/=TT_REPLICATED.and.totype/=TT_DENSE.and.totype/=TT_TILED_DIST.and.totype/=TT_TILED)then
        call tensor_status_quit("ERROR(change_access_type_td): wrong type given",-1)
     endif
     if(me==root.and.arr%access_type==AT_MASTER_ACCESS) then
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_COMM )
#endif
        call pdm_tensor_sync(JOB_CHANGE_access_type,arr)
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_WORK )
#endif
     endif
     arr%access_type=totype
#endif
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
  end subroutine change_access_type_td



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!                 ACCUMULATE TILES
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> \brief direct communication routine for the accumulation of arrays,
  !> interface to the combined index routine
  !> \author Patrick Ettenhuber
  subroutine tensor_acctile_modeidx(arr,modidx,fort,nelms,lock_set,flush_it,req)
     implicit none
     !> input array for which a tile should be accumulated
     type(tensor),intent(in) ::arr
     !> input, the index of the tile in modular form and the number of elements
     integer,intent(in) :: modidx(arr%mode),nelms
     !> input the fortan array which should be transferred to the tile
     real(tensor_dp),intent(inout) :: fort(nelms)
     logical, optional, intent(in) :: lock_set,flush_it
     integer(kind=tensor_mpi_kind),intent(inout), optional :: req
     integer(kind=tensor_long_int) :: gt
     integer(kind=tensor_long_int) :: n
     gt = get_cidx(modidx,arr%ntpm,arr%mode)
     n  = nelms
     call tensor_acc_tile(arr,gt,fort,n,lock_set=lock_set,flush_it=flush_it,req=req)
  end subroutine tensor_acctile_modeidx
  subroutine tensor_acct44(arr,globtilenr,fort,nelms,lock_set,flush_it,req)
     implicit none
     type(tensor),intent(in) :: arr
     integer(kind=tensor_standard_int),intent(in) :: globtilenr
     integer(kind=tensor_standard_int),intent(in) :: nelms
     real(tensor_dp),intent(inout) :: fort(nelms)
     logical, optional, intent(in) :: lock_set,flush_it
     integer(kind=tensor_mpi_kind),intent(inout), optional :: req
     integer(kind=tensor_long_int) :: gt
     integer(kind=tensor_long_int) :: n
     gt = globtilenr
     n  = nelms
     call tensor_accumulate_tile_combidx(arr,gt,fort,n,lock_set=lock_set,flush_it=flush_it,req=req)
  end subroutine tensor_acct44
  subroutine tensor_acct48(arr,globtilenr,fort,nelms,lock_set,flush_it,req)
     implicit none
     type(tensor),intent(in) :: arr
     integer(kind=tensor_standard_int),intent(in) :: globtilenr
     integer(kind=tensor_long_int),intent(in) :: nelms
     real(tensor_dp),intent(inout) :: fort(nelms)
     logical, optional, intent(in) :: lock_set,flush_it
     integer(kind=tensor_mpi_kind),intent(inout), optional :: req
     integer(kind=tensor_long_int) :: gt
     integer(kind=tensor_long_int) :: n
     gt = globtilenr
     n  = nelms
     call tensor_accumulate_tile_combidx(arr,gt,fort,n,lock_set=lock_set,flush_it=flush_it,req=req)
  end subroutine tensor_acct48
  subroutine tensor_acct84(arr,globtilenr,fort,nelms,lock_set,flush_it,req)
     implicit none
     type(tensor),intent(in) :: arr
     integer(kind=tensor_long_int),intent(in) :: globtilenr
     integer(kind=tensor_standard_int),intent(in) :: nelms
     real(tensor_dp),intent(inout) :: fort(nelms)
     logical, optional, intent(in) :: lock_set,flush_it
     integer(kind=tensor_mpi_kind),intent(inout), optional :: req
     integer(kind=tensor_long_int) :: gt
     integer(kind=tensor_long_int) :: n
     gt = globtilenr
     n  = nelms
     call tensor_accumulate_tile_combidx(arr,gt,fort,n,lock_set=lock_set,flush_it=flush_it,req=req)
  end subroutine tensor_acct84
  subroutine tensor_acct88(arr,globtilenr,fort,nelms,lock_set,flush_it,req)
     implicit none
     type(tensor),intent(in) :: arr
     integer(kind=tensor_long_int),intent(in) :: globtilenr
     integer(kind=tensor_long_int),intent(in) :: nelms
     real(tensor_dp),intent(inout) :: fort(nelms)
     logical, optional, intent(in) :: lock_set,flush_it
     integer(kind=tensor_mpi_kind),intent(inout), optional :: req
     integer(kind=tensor_long_int) :: gt
     integer(kind=tensor_long_int) :: n
     gt = globtilenr
     n  = nelms
     call tensor_accumulate_tile_combidx(arr,gt,fort,n,lock_set=lock_set,flush_it=flush_it,req=req)
  end subroutine tensor_acct88
  subroutine tensor_accumulate_tile_combidx(arr,gt,fort,nelms,lock_set,flush_it,req)
     implicit none
     type(tensor),intent(in) :: arr
     integer(kind=tensor_long_int),intent(in) :: gt
     integer(kind=tensor_long_int),intent(in) :: nelms
     !> input the fortan array which should be transferred to the tile
     real(tensor_dp),intent(inout) :: fort(nelms)
     logical, optional, intent(in) :: lock_set, flush_it
     integer(kind=tensor_mpi_kind),intent(inout), optional :: req
     integer(kind=tensor_mpi_kind) :: dest
     logical :: ls
     integer(kind=tensor_standard_int) :: dpos
#ifdef VAR_MPI
     integer(kind=tensor_long_int) :: p,pos,widx
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_COMM )
#endif

     ls = .false.
     if(present(lock_set))ls=lock_set

     call get_residence_of_tile( arr, gt, dest, dpos, p, widx)

     time_pdm_acc = time_pdm_acc - mpi_wtime()

     if(.not.ls)call tensor_mpi_win_lock(dest,arr%wi(widx),'s')

     if(present(req))then
        call tensor_mpi_acc(fort,nelms,p,dest,arr%wi(widx),req, fl = flush_it)
     else
        call tensor_mpi_acc(fort,nelms,p,dest,arr%wi(widx), fl = flush_it)
     endif

     if(.not.ls)CALL tensor_mpi_win_unlock(dest, arr%wi(widx))

     time_pdm_acc = time_pdm_acc + mpi_wtime()
     bytes_transferred_acc = bytes_transferred_acc + nelms * tensor_dp
     nmsg_acc     = nmsg_acc + 1

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
#endif
  end subroutine tensor_accumulate_tile_combidx



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!                   PUT TILES
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tensor_puttile_modeidx(arr,modidx,fort,nelms,lock_set,flush_it,req)
     implicit none
     type(tensor),intent(in) ::arr
     integer,intent(in) :: modidx(arr%mode),nelms
     real(tensor_dp),intent(inout) :: fort(nelms)
     logical, optional, intent(in) :: lock_set,flush_it
     integer(kind=tensor_mpi_kind), intent(inout), optional :: req
     integer(kind=tensor_long_int) :: gt
     integer(kind=tensor_long_int) :: n

     gt= get_cidx(modidx,arr%ntpm,arr%mode)
     n = nelms

     call tensor_put_tile(arr,gt,fort,n,lock_set=lock_set,flush_it=flush_it,req=req)
  end subroutine tensor_puttile_modeidx

  subroutine tensor_putt88(arr,globtilenr,fort,nelms,lock_set,flush_it,req)
     implicit none
     type(tensor),intent(in) :: arr
     integer(kind=tensor_long_int), intent(in) :: globtilenr
     integer(kind=tensor_long_int), intent(in) :: nelms
     real(tensor_dp),intent(inout) :: fort(nelms)
     logical, optional, intent(in) :: lock_set,flush_it
     integer(kind=tensor_mpi_kind), intent(inout), optional :: req
     integer(kind=tensor_long_int) :: gt
     integer(kind=tensor_long_int) :: n
     gt = globtilenr
     n  = nelms
     call tensor_puttile_combidx(arr,gt,fort,n,lock_set=lock_set,flush_it=flush_it,req=req)
  end subroutine tensor_putt88
  subroutine tensor_putt84(arr,globtilenr,fort,nelms,lock_set,flush_it,req)
     implicit none
     type(tensor),intent(in) :: arr
     integer(kind=tensor_long_int), intent(in) :: globtilenr
     integer(kind=tensor_standard_int), intent(in) :: nelms
     real(tensor_dp),intent(inout) :: fort(nelms)
     logical, optional, intent(in) :: lock_set,flush_it
     integer(kind=tensor_mpi_kind), intent(inout), optional :: req
     integer(kind=tensor_long_int) :: gt
     integer(kind=tensor_long_int) :: n
     gt = globtilenr
     n  = nelms
     call tensor_puttile_combidx(arr,gt,fort,n,lock_set=lock_set,flush_it=flush_it,req=req)
  end subroutine tensor_putt84
  subroutine tensor_putt48(arr,globtilenr,fort,nelms,lock_set,flush_it,req)
     implicit none
     type(tensor),intent(in) :: arr
     integer(kind=tensor_standard_int), intent(in) :: globtilenr
     integer(kind=tensor_long_int), intent(in) :: nelms
     real(tensor_dp),intent(inout) :: fort(nelms)
     logical, optional, intent(in) :: lock_set,flush_it
     integer(kind=tensor_mpi_kind), intent(inout), optional :: req
     integer(kind=tensor_long_int) :: gt
     integer(kind=tensor_long_int) :: n
     gt = globtilenr
     n  = nelms
     call tensor_puttile_combidx(arr,gt,fort,n,lock_set=lock_set,flush_it=flush_it,req=req)
  end subroutine tensor_putt48
  subroutine tensor_putt44(arr,globtilenr,fort,nelms,lock_set,flush_it,req)
     implicit none
     type(tensor),intent(in) :: arr
     integer(kind=tensor_standard_int), intent(in) :: globtilenr
     integer(kind=tensor_standard_int), intent(in) :: nelms
     real(tensor_dp),intent(inout) :: fort(nelms)
     logical, optional, intent(in) :: lock_set,flush_it
     integer(kind=tensor_mpi_kind), intent(inout), optional :: req
     integer(kind=tensor_long_int) :: gt
     integer(kind=tensor_long_int) :: n
     gt = globtilenr
     n  = nelms
     call tensor_puttile_combidx(arr,gt,fort,n,lock_set=lock_set,flush_it=flush_it,req=req)
  end subroutine tensor_putt44

  subroutine tensor_puttile_combidx(arr,globtilenr,fort,nelms,lock_set,flush_it,req)
     implicit none
     type(tensor),intent(in) :: arr
     integer(kind=tensor_long_int),intent(in) :: globtilenr
     integer(kind=tensor_long_int),intent(in) :: nelms
     real(tensor_dp),intent(inout) :: fort(nelms)
     logical, optional, intent(in) :: lock_set,flush_it
     integer(kind=tensor_mpi_kind), intent(inout), optional :: req
     logical :: ls
     integer(kind=tensor_mpi_kind) :: dest
     integer :: pos,dpos
     integer(kind=tensor_long_int) :: p,widx
     integer(kind=tensor_standard_int) :: gt
#ifdef VAR_MPI
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_COMM )
#endif

     gt = globtilenr

     ls = .false.
     if(present(lock_set))ls=lock_set

     call get_residence_of_tile(arr,gt,dest,dpos,p,widx)

     time_pdm_put          = time_pdm_put - mpi_wtime()

     if(.not.ls)call tensor_mpi_win_lock(dest,arr%wi(widx),'s')

     if(present(req))then
        call tensor_mpi_put(fort,nelms,p,dest,arr%wi(widx),req, fl = flush_it)
     else
        call tensor_mpi_put(fort,nelms,p,dest,arr%wi(widx), fl = flush_it)
     endif

     if(.not.ls)call tensor_mpi_win_unlock(dest,arr%wi(widx))

     time_pdm_put          = time_pdm_put + mpi_wtime()
     bytes_transferred_put = bytes_transferred_put + nelms * tensor_dp
     nmsg_put              = nmsg_put + 1

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
#endif
  end subroutine tensor_puttile_combidx

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!                   GET TILES
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !interface to the tensor_gettile_combidx
  subroutine tensor_gettile_modeidx(arr,modidx,fort,nelms,lock_set,flush_it,req)
     implicit none
     type(tensor),intent(in) ::arr
     integer,intent(in) :: modidx(arr%mode),nelms
     real(tensor_dp),intent(inout) :: fort(nelms)
     logical, optional, intent(in) :: lock_set,flush_it
     integer(kind=tensor_mpi_kind),intent(inout),optional :: req
     logical :: ls
     integer(kind=tensor_standard_int) :: cidx
     ls = .false.
     if(present(lock_set))ls=lock_set
     cidx=get_cidx(modidx,arr%ntpm,arr%mode)
     call tensor_get_tile(arr,cidx,fort,nelms,lock_set=ls,flush_it=flush_it,req=req)
  end subroutine tensor_gettile_modeidx
  subroutine tensor_gett88(arr,globtilenr,fort,nelms,lock_set,flush_it,req)
     implicit none
     type(tensor),intent(in) :: arr
     integer(kind=tensor_long_int),intent(in) :: globtilenr
     integer(kind=tensor_long_int),intent(in) :: nelms
     real(tensor_dp),intent(inout) :: fort(nelms)
     logical, optional, intent(in) :: lock_set,flush_it
     integer(kind=tensor_mpi_kind),intent(inout),optional :: req
     call tensor_gettile_combidx(arr,globtilenr,&
        &fort,nelms,lock_set=lock_set,flush_it=flush_it,req=req)
  end subroutine tensor_gett88
  subroutine tensor_gett48(arr,globtilenr,fort,nelms,lock_set,flush_it,req)
     implicit none
     type(tensor),intent(in) :: arr
     integer(kind=tensor_standard_int),intent(in) :: globtilenr
     integer(kind=tensor_long_int),intent(in) :: nelms
     real(tensor_dp),intent(inout) :: fort(nelms)
     logical, optional, intent(in) :: lock_set,flush_it
     integer(kind=tensor_mpi_kind),intent(inout),optional :: req
     call tensor_gettile_combidx(arr,int(globtilenr,kind=tensor_long_int),&
        &fort,nelms,lock_set=lock_set,flush_it=flush_it,req=req)
  end subroutine tensor_gett48
  subroutine tensor_gett84(arr,globtilenr,fort,nelms,lock_set,flush_it,req)
     implicit none
     type(tensor),intent(in) :: arr
     integer(kind=tensor_long_int),intent(in) :: globtilenr
     integer(kind=tensor_standard_int),intent(in) :: nelms
     real(tensor_dp),intent(inout) :: fort(nelms)
     logical, optional, intent(in) :: lock_set,flush_it
     integer(kind=tensor_mpi_kind),intent(inout),optional :: req
     call tensor_gettile_combidx(arr,globtilenr,&
        &fort,int(nelms,kind=tensor_long_int),lock_set=lock_set,flush_it=flush_it,req=req)
  end subroutine tensor_gett84
  subroutine tensor_gett44(arr,globtilenr,fort,nelms,lock_set,flush_it,req)
     implicit none
     type(tensor),intent(in) :: arr
     integer(kind=tensor_standard_int),intent(in) :: globtilenr
     integer(kind=tensor_standard_int),intent(in) :: nelms
     real(tensor_dp),intent(inout) :: fort(nelms)
     logical, optional, intent(in) :: lock_set,flush_it
     integer(kind=tensor_mpi_kind),intent(inout),optional :: req
     call tensor_gettile_combidx(arr,int(globtilenr,kind=tensor_long_int),&
        &fort,int(nelms,kind=tensor_long_int),lock_set=lock_set,flush_it=flush_it,req=req)
  end subroutine tensor_gett44

  subroutine tensor_gettile_combidx(arr,globtilenr,fort,nelms,lock_set,flush_it,req)
     implicit none
     type(tensor),intent(in) :: arr
     integer(kind=tensor_long_int),intent(in) :: globtilenr
     integer(kind=tensor_long_int),intent(in) :: nelms
     real(tensor_dp),intent(inout) :: fort(nelms)
     logical, optional, intent(in) :: lock_set,flush_it
     integer(kind=tensor_mpi_kind),intent(inout),optional :: req
     integer(kind=tensor_mpi_kind) :: source,r
     integer(kind=tensor_standard_int) :: gt
     logical :: ls
#ifdef VAR_MPI
     integer(kind=tensor_long_int) :: p, dpos, widx
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_COMM )
#endif

     gt = globtilenr

     ls = .false.
     if(present(lock_set))ls=lock_set

     call get_residence_of_tile(arr,gt,source,dpos,p,widx)

     time_pdm_get = time_pdm_get - mpi_wtime()

     if(.not.ls)call tensor_mpi_win_lock(source,arr%wi(widx),'s')

     if(present(req))then
        call tensor_mpi_get(fort,nelms,p,source,arr%wi(widx),req, fl = flush_it)
     else
!#ifdef VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN
!        call tensor_mpi_get(fort,nelms,int(p),source,arr%wi(widx),r,fl=.true.)
!        call tensor_mpi_wait(r)
!        call tensor_mpi_win_flush(arr%wi(widx), rank=source, local=.false.)
!#else
        call tensor_mpi_get(fort,nelms,p,source,arr%wi(widx), fl = flush_it)
     endif

     if(.not.ls)call tensor_mpi_win_unlock(source,arr%wi(widx))

     time_pdm_get = time_pdm_get + mpi_wtime()
     bytes_transferred_get = bytes_transferred_get + nelms * tensor_dp
     nmsg_get = nmsg_get + 1

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
#endif
  end subroutine tensor_gettile_combidx

  subroutine get_int_dist_info(o2v2,firstintel,nintel,remoterank)
     implicit none
     integer(kind=tensor_long_int), intent(in) :: o2v2
     integer, intent(inout) :: firstintel,nintel
     integer(kind=tensor_mpi_kind), intent(in), optional :: remoterank
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: nnod
     integer(kind=tensor_mpi_kind) :: me
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase(PHASE_WORK)
#endif

     call tensor_get_rank(me)
     call tensor_get_size(nnod)

     !overwirte to calculate the info for another node
     if(present(remoterank))then
        me = remoterank
     endif

     nintel = o2v2/nnod
     firstintel = me*nintel + 1
     if(me<int(mod(o2v2,int(nnod,kind=tensor_long_int)),kind=tensor_mpi_kind))then
        nintel = nintel + 1
        firstintel = firstintel + int(me) 
     else if(me>=int(mod(o2v2,int(nnod,kind=tensor_long_int)),kind=tensor_mpi_kind))then
        firstintel = firstintel + int(mod(o2v2,int(nnod,kind=tensor_long_int))) 
     endif

#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif
  end subroutine get_int_dist_info


  subroutine tensor_scale_td(arr,sc)
     implicit none
     type(tensor) :: arr
     real(tensor_dp) :: sc
#ifdef VAR_MPI
     integer     :: i
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: me

     call tensor_get_rank(me)

     if(arr%access_type==AT_MASTER_ACCESS.and.me==root)then
#ifdef TENSORS_IN_LSDALTON
        call time_start_phase( PHASE_COMM )
#endif
        call pdm_tensor_sync(JOB_tensor_SCALE,arr)
        call tensor_mpi_bcast(sc,root,tensor_work_comm)
     endif
#ifdef TENSORS_IN_LSDALTON
     call time_start_phase( PHASE_WORK )
#endif

     do i=1,arr%nlti
        call dscal(int(arr%ti(i)%e),sc,arr%ti(i)%t,1)
     enddo

     if( tensor_always_sync ) call tensor_mpi_barrier(tensor_work_comm)

#endif
  end subroutine tensor_scale_td


  subroutine tensor_flush_win(T,node,gtidx,local,only_owner)
     implicit none
     type(tensor) :: T
     integer(kind=tensor_mpi_kind), intent(in), optional :: node
     integer, intent(in), optional :: gtidx
     logical, intent(in), optional :: local,only_owner
     integer :: tidx,n
     logical :: all_tiles,oo
     integer :: widx, pos, idx
     integer(kind=tensor_mpi_kind)     :: node2
     integer(kind=tensor_standard_int) :: gt


     all_tiles = .not.present(gtidx)
     if(present(gtidx)) gt = gtidx

     if(present(node)) n = node

     oo = .false.
     if(present(only_owner)) oo = only_owner

#ifdef VAR_MPI
#ifdef VAR_HAVE_MPI3
     if( oo )then

        if( all_tiles )then

           do widx = 1, T%nwins

              call tensor_mpi_win_flush(T%wi(widx),rank=node,local=local)

           enddo

        else

           call get_residence_of_tile(T,gt, node2, pos, idx, widx)

           call tensor_mpi_win_flush(T%wi(widx),rank=node2,local=local)

        endif
     else
        if( all_tiles )then

           do widx = 1, T%nwins
              call tensor_mpi_win_flush(T%wi(widx),rank=node,local=local)
           enddo

        else

           if( alloc_in_dummy ) then
              widx = 1
           else
              widx = gt
           endif

           call tensor_mpi_win_flush(T%wi(widx),rank=node,local=local)

        endif

     endif
#else

     print *,"WARNING(tensor_flush_win): it is not recommended to use this &
        &without mpi 3, if you have mpi 3 make sure you comile with VAR_HAVE_MPI3"

     if( all_tiles )then

        do widx = 1,T%nwins

           if(T%lock_set(widx))call tensor_unlock_win(T,widx)
           call tensor_lock_win(T,widx,'s')

        enddo

     else

        if( alloc_in_dummy ) then
           widx = 1
        else
           widx = gt
        endif

        if(T%lock_set(widx))call tensor_unlock_win(T,widx)
        call tensor_lock_win(T,widx,'s')

     endif

     !ENDIF VAR_HAVE_MPI3
#endif 

     !ENDIF VAR_MPI
#endif
  end subroutine tensor_flush_win
  subroutine assoc_ptr_to_buf(tilenr,arr,nbuffs,buf_pos,buf_log,ptr,bg_buf,pos,req)
     implicit none
     integer, intent(in):: tilenr,nbuffs
     type(tensor), intent(inout) :: arr
     integer, intent(inout):: buf_pos(nbuffs)
     logical, intent(inout):: buf_log(nbuffs)
     real(tensor_dp), intent(out),   pointer :: ptr(:)
     real(tensor_dp), intent(inout), pointer :: bg_buf(:,:)
     integer, intent(out) :: pos
     integer(kind=tensor_mpi_kind), intent(inout) :: req(nbuffs)
     integer :: i_search_buf
     integer(kind=tensor_long_int) :: ts
     logical :: found
     integer(kind=tensor_mpi_kind) :: mode
     integer(kind=tensor_standard_int) :: tnr
     pos = 0
     tnr = tilenr
#ifdef VAR_MPI
     mode = MPI_MODE_NOCHECK

     call find_tile_pos_in_buf(tnr,buf_pos,nbuffs,pos,found)

     if(found)then

        ptr => bg_buf(:,pos)

     else

        call find_free_pos_in_buf(buf_log,nbuffs,pos,found)

        if( .not. found)then

           call tensor_status_quit("ERROR(assoc_ptr_to_buf):tile not found in buf and no&
              & free position available to load",-1)

        endif

        if( .not. alloc_in_dummy ) call tensor_lock_win(arr,tnr,'s',assert=mode)

        call get_tile_dim(ts,arr,tnr)

        if( alloc_in_dummy )then
           call tensor_get_tile(arr,tnr,bg_buf(:,pos),ts,lock_set=.true.,req=req(pos))
        else
           call tensor_get_tile(arr,tnr,bg_buf(:,pos),ts,lock_set=.true.,flush_it=.true.)
        endif

        buf_pos(pos) = tnr
        buf_log(pos) = .true.
        ptr          => bg_buf(:,pos)

     endif
#endif
  end subroutine assoc_ptr_to_buf

  subroutine find_tile_pos_in_buf(tilenr,buf,nbuffs,pos,found)
     implicit none
     integer(kind=tensor_standard_int), intent(in)  :: tilenr
     integer, intent(in)  :: nbuffs
     integer, intent(in)  :: buf(nbuffs)
     logical, intent(out) :: found
     integer, intent(out) :: pos
     integer :: i

     found = .false.
     do i=1,nbuffs
        if(buf(i)==tilenr)then
           pos   = i
           found = .true.
           exit
        endif
     enddo

  end subroutine find_tile_pos_in_buf

  subroutine find_free_pos_in_buf(buf,nbuffs,pos,found)
     implicit none
     integer, intent(in)  :: nbuffs
     logical, intent(in)  :: buf(nbuffs)
     logical, intent(out) :: found
     integer, intent(out) :: pos
     integer :: i

     found = .false.
     do i=1,nbuffs
        if(.not.buf(i))then
           pos   = i
           found = .true.
           exit
        endif
     enddo

  end subroutine find_free_pos_in_buf

  subroutine check_if_new_instance_needed(tilenr,buf,nbuffs,NewN,pos,set_needed)
     implicit none
     integer, intent(in)  :: tilenr
     integer, intent(in)  :: nbuffs
     integer, intent(in)  :: buf(nbuffs)
     logical, intent(out) :: NewN
     integer, intent(out), optional :: pos
     logical, intent(inout), optional :: set_needed(nbuffs)
     logical :: found
     integer :: pos_int
     integer(kind=tensor_standard_int) :: tnr
     tnr = tilenr

     call find_tile_pos_in_buf(tnr,buf,nbuffs,pos_int,found) 

     NewN = ( .not. found )

     if(present(pos)) pos = pos_int

     if(present(set_needed))then
        if( found ) set_needed(pos_int) = .true.
     endif

  end subroutine check_if_new_instance_needed

  subroutine tensor_barrier(call_slaves)
     implicit none
     logical, intent(in) :: call_slaves
     integer(kind=tensor_mpi_kind), parameter :: root = 0
     integer(kind=tensor_mpi_kind) :: nnod
     integer(kind=tensor_mpi_kind) :: me
#ifdef VAR_MPI
     call tensor_get_rank(me)
     call tensor_get_size(nnod)
     if(me==root.and.call_slaves) then
        call pdm_tensor_sync(JOB_BARRIER)
     endif

     call tensor_mpi_barrier(tensor_work_comm)
#endif
  end subroutine tensor_barrier
end module lspdm_tensor_operations_module




