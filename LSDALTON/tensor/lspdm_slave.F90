
subroutine pdm_tensor_slave(comm)
  use tensor_parameters_and_counters
  use tensor_allocator
  use lstiming
  !use matrix_operations_scalapack, only: BLOCK_SIZE, SLGrid, DLEN_
  use matrix_operations, only: mtype_scalapack, matrix_type
  use dec_typedef_module
#ifdef VAR_MPI
  use infpar_module
  use tensor_mpi_interface_module
#endif
  use tensor_mpi_operations_module

  !this is a special case since this slave routine is associated with the tensor
  !structure it is allwed to use lspdm_tensor_operaions, otherwise this is not
  !allowed
  use lspdm_tensor_operations_module
  use tensor_interface_module

   IMPLICIT NONE
   INTEGER(kind=tensor_mpi_kind),intent(in) :: comm
   type(tensor)  :: A, B, C, D, AUX
   CHARACTER    :: T(2)
   INTEGER      :: JOB
   real(tensor_dp),pointer :: realar1(:)
   integer(kind=tensor_long_int), pointer  :: lintar1(:), lintar2(:)
   integer, pointer  :: intarr1(:), intarr2(:), intarr3(:), intarr4(:)
   integer(kind=tensor_long_int) :: LIN1
   integer           :: INT1,       INT2,       INT3,       INT4
   real(tensor_dp)   :: REAL1,      REAL2
   logical           :: LOG1
   logical :: loc
   character (4) :: at 
#ifdef VAR_MPI

   call time_start_phase(PHASE_COMM)
   CALL PDM_tensor_SYNC(comm,JOB,A,B,C,D) !Job is output
   call time_start_phase(PHASE_WORK)

   SELECT CASE(JOB)
   !CASE(JOB_PC_DEALLOC_DENSE)
   !   !call memory_deallocate_tensor_dense(A)
   !   print *,"this is deprecated 1"
   !   stop 0
   !CASE(JOB_PC_ALLOC_DENSE)
   !   !call memory_allocate_tensor_dense(A)
   !   print *,"this is deprecated 2"
   !   stop 0

   CASE(JOB_INIT_TENSOR_TILED)

      call tensor_alloc_mem(intarr2,A%mode)
      call tensor_alloc_mem(intarr1,A%mode)
      intarr2=A%tdim
      intarr1 =A%dims
      call tensor_free_aux(A)

      call tensor_mpi_bcast(INT1,infpar%master,comm)

      if(INT1==-1)then
         call tensor_init_tiled(A,intarr1,int(A%mode),at,int(INT1,kind=tensor_standard_int),&
            &AT_MASTER_ACCESS,.false.,tdims=intarr2) 
      else
         call tensor_init_tiled(A,intarr1,int(A%mode),at,int(INT1,kind=tensor_standard_int),&
            &AT_MASTER_ACCESS,.false.,tdims=intarr2,force_offset=INT1) 
      endif

      call tensor_free_mem(intarr2)
      call tensor_free_mem(intarr1)

   CASE(JOB_FREE_TENSOR_STD)
      call tensor_free_pdm(A) 
   CASE(JOB_FREE_tensor_PDM)
      call tensor_free_pdm(A) 
   CASE(JOB_INIT_tensor_REPLICATED)
      call tensor_alloc_mem(intarr1,A%mode)
      intarr1 =A%dims
      call tensor_free_aux(A)
      call tensor_init_replicated(A,intarr1,int(A%mode),AT_MASTER_ACCESS,.false.) 
      call tensor_free_mem(intarr1)
   CASE(JOB_PRINT_MEM_INFO1)
      call print_mem_per_node(DECinfo%output,.false.)
   CASE(JOB_PRINT_MEM_INFO2)
      call tensor_alloc_mem(lintar1,9*infpar%lg_nodtot)
      call tensor_alloc_mem(lintar2,infpar%lg_nodtot)
      lintar1 = 0
      lintar2 = 0
      call print_mem_per_node(DECinfo%output,.false.,lintar1,lintar2)
      call tensor_free_mem(lintar1)
      call tensor_free_mem(lintar2)
   CASE(JOB_GET_NRM2_TILED)
      REAL1 = tensor_tiled_pdm_get_nrm2(A)
   CASE(JOB_DATA2TILED_DIST)
      !realar1 has to be allocated, otherwise seg faults might occur with 
      !some compilers
      call tensor_alloc_mem(realar1,1)
      !call cp_data2tiled(A,realar1,A%dims,A%mode,.true.)
      print *,"not necessary"
      call tensor_free_mem(realar1)
   CASE(JOB_GET_TILE_SEND)
      !INT1,realar1 and INT2 are just dummy arguments
      call tensor_get_tile(A,int(INT1,kind=tensor_standard_int),realar1,INT2)
   CASE(JOB_PRINT_TI_NRM)
      call tensor_tiled_pdm_print_ti_nrm(A,0)
   CASE(JOB_SYNC_REPLICATED)
      call tensor_sync_replicated(A)
   CASE(JOB_GET_NORM_REPLICATED)
      REAL1 = tensor_print_norm_repl(A)
   CASE(JOB_PREC_DOUBLES_PAR)
      call precondition_doubles_parallel(A,B,C,D)
   CASE(JOB_DDOT_PAR)
      REAL1 = tensor_ddot_par(A,B,0)
   CASE(JOB_ADD_PAR)
      INT1 = A%mode
      call tensor_alloc_mem(intarr1,INT1)
      call time_start_phase(PHASE_COMM)
      call tensor_buffer(intarr1,INT1,root=infpar%master,comm=comm)
      call tensor_buffer(REAL1)
      call tensor_buffer(REAL2,finalize=.true.)
      call time_start_phase(PHASE_WORK)

      call tensor_add_par(REAL1,A,REAL2,B,intarr1)

      call tensor_free_mem(intarr1)

   CASE(JOB_DMUL_PAR)
      INT1 = A%mode
      call tensor_alloc_mem(intarr1,INT1)
      call time_start_phase(PHASE_COMM)
      call tensor_buffer(intarr1,INT1,root=infpar%master,comm=comm)
      call tensor_buffer(REAL1)
      call tensor_buffer(REAL2)
      call tensor_buffer(INT2)
      call tensor_alloc_mem(realar1,INT2)
      call tensor_buffer(realar1,INT2,finalize=.true.)
      call time_start_phase(PHASE_WORK)

      call tensor_dmul_par(REAL1,A,REAL2,realar1,B,intarr1)

      call tensor_free_mem(intarr1)
      call tensor_free_mem(realar1)

   CASE(JOB_HMUL_PAR)
      call time_start_phase(PHASE_COMM)
      call tensor_buffer(REAL1,root=infpar%master,comm=comm)
      call tensor_buffer(REAL2,finalize=.true.)
      call time_start_phase(PHASE_WORK)
      call tensor_hmul_par(REAL1,A,B,REAL2,C)
   CASE(JOB_CP_ARR)
      INT1 = A%mode
      call tensor_alloc_mem(intarr1,INT1)
      call tensor_mpi_bcast(intarr1,INT1,infpar%master,infpar%lg_comm)
      call tensor_cp_tiled(A,B,intarr1)
      call tensor_free_mem(intarr1)
   CASE(JOB_tensor_ZERO)
      call tensor_zero_tiled_dist(A)
   CASE(JOB_GET_CC_ENERGY)
      REAL1 = get_cc_energy_parallel(A,B,C)
   CASE(JOB_GET_MP2_ENERGY)
      REAL1 = get_cc_energy_parallel(A,B)
   CASE(JOB_GET_FRAG_CC_ENERGY)
      !the counterpart to this buffer is in get_fragment_cc_energy
      call time_start_phase(PHASE_COMM)
      call tensor_buffer(INT1,root=infpar%master,comm=comm)
      call tensor_alloc_mem(intarr1,INT1)
      call tensor_buffer(intarr1,INT1)
      call tensor_buffer(INT2)
      call tensor_alloc_mem(intarr2,INT2)
      call tensor_buffer(intarr2,INT2,finalize=.true.)
      call time_start_phase(PHASE_WORK)

      REAL1 = get_fragment_cc_energy_parallel(A,B,C,INT1,INT2,intarr1,intarr2)

      call tensor_free_mem(intarr1)
      call tensor_free_mem(intarr2)
   CASE(JOB_CHANGE_ACCESS_TYPE)
      call change_access_type_td(A,INT1)
   CASE(JOB_tensor_SCALE)
      call tensor_mpi_bcast(REAL1,infpar%master,infpar%lg_comm)
      call tensor_scale_td(A,REAL1)
   CASE(JOB_GET_RPA_ENERGY)
      write(*,*) 'Johannes RPA ene'
      REAL1 = get_rpa_energy_parallel(A,B)
   CASE(JOB_GET_SOS_ENERGY)
      write(*,*) 'Johannes SOS ene'
      REAL1 = get_sosex_cont_parallel(A,B)
   CASE(JOB_TENSOR_CONTRACT_SIMPLE)
      call time_start_phase(PHASE_COMM)
      call tensor_buffer(INT1,root=infpar%master,comm=comm)

      INT2 = INT1
      INT3 = C%mode
      call tensor_alloc_mem(intarr1,INT1)
      call tensor_alloc_mem(intarr2,INT2)
      call tensor_alloc_mem(intarr3,INT3)

      call tensor_buffer(intarr1,INT1)
      call tensor_buffer(intarr2,INT2)
      call tensor_buffer(intarr3,INT3)
      call tensor_buffer(REAL1)
      call tensor_buffer(REAL2)
      call tensor_buffer(LOG1, finalize = .true.)
      call time_start_phase(PHASE_WORK)

      if(LOG1)then
         call lspdm_tensor_contract_simple(REAL1,A,B,intarr1,intarr2,INT1,REAL2,C,intarr3,force_sync=LOG1)
      else
         call lspdm_tensor_contract_simple(REAL1,A,B,intarr1,intarr2,INT1,REAL2,C,intarr3)
      endif

      call tensor_free_mem(intarr1)
      call tensor_free_mem(intarr2)
      call tensor_free_mem(intarr3)
   CASE(JOB_TENSOR_CONTRACT_BDENSE)
      call time_start_phase(PHASE_COMM)
      call tensor_buffer(INT1,root=infpar%master,comm=comm)
      INT2 = INT1
      INT3 = B%mode
      call tensor_alloc_mem(intarr1,INT1)
      call tensor_alloc_mem(intarr2,INT2)
      call tensor_alloc_mem(intarr3,INT3)
      call tensor_buffer(intarr1,INT1)
      call tensor_buffer(intarr2,INT2)
      call tensor_buffer(intarr3,INT3)
      call tensor_buffer(REAL1)
      call tensor_buffer(REAL2)
      call tensor_buffer(LOG1)
      call tensor_buffer(INT4)
      call tensor_alloc_mem(intarr4,INT4)
      call tensor_buffer(intarr4,INT4)

      call tensor_init(AUX,intarr4,int4)

      call tensor_buffer(AUX%elm1,AUX%nelms, finalize=.true.)

      call time_start_phase(PHASE_WORK)

      if(LOG1)then
         call lspdm_tensor_contract_simple(REAL1,A,AUX,intarr1,intarr2,INT1,REAL2,B,intarr3,force_sync=LOG1)
      else
         call lspdm_tensor_contract_simple(REAL1,A,AUX,intarr1,intarr2,INT1,REAL2,B,intarr3)
      endif

      call tensor_free(AUX)

      call tensor_free_mem(intarr1)
      call tensor_free_mem(intarr2)
      call tensor_free_mem(intarr3)
      call tensor_free_mem(intarr4)
   CASE(JOB_TENSOR_EXTRACT_VEOS)
      call tensor_buffer(INT1,root=infpar%master,comm=comm)
      call tensor_alloc_mem(intarr1,INT1)
      call tensor_buffer(intarr1,INT1)
      call tensor_buffer(INT2)
      call tensor_alloc_mem(intarr2,INT2)
      call tensor_buffer(intarr2,INT2,finalize=.true.)

      call tensor_init(AUX,intarr2,INT2)

      call tensor_free_mem(intarr2)

      call tensor_zero(AUX)

      call lspdm_extract_eos_indices_virt(AUX,A,INT1,intarr1)

      call tensor_free_mem(intarr1)
      call tensor_free(AUX)

   CASE(JOB_TENSOR_EXTRACT_OEOS)
      call tensor_buffer(INT1,root=infpar%master,comm=comm)
      call tensor_alloc_mem(intarr1,INT1)
      call tensor_buffer(intarr1,INT1)
      call tensor_buffer(INT2)
      call tensor_alloc_mem(intarr2,INT2)
      call tensor_buffer(intarr2,INT2,finalize=.true.)
      
      call tensor_init(AUX,intarr2,INT2)

      call tensor_free_mem(intarr2)

      call tensor_zero(AUX)

      call lspdm_extract_eos_indices_occ(AUX,A,INT1,intarr1)

      call tensor_free_mem(intarr1)
      call tensor_free(AUX)

   CASE(JOB_TENSOR_EXTRACT_VDECNP)
      call tensor_buffer(INT1,root=infpar%master,comm=comm)
      call tensor_alloc_mem(intarr1,INT1)
      call tensor_buffer(intarr1,INT1)
      call tensor_buffer(INT2)
      call tensor_alloc_mem(intarr2,INT2)
      call tensor_buffer(intarr2,INT2,finalize=.true.)
      
      call tensor_init(AUX,intarr2,INT2)

      call tensor_free_mem(intarr2)

      call tensor_zero(AUX)

      call lspdm_extract_decnp_indices_virt(AUX,A,INT1,intarr1)

      call tensor_free_mem(intarr1)
      call tensor_free(AUX)

   CASE(JOB_TENSOR_EXTRACT_ODECNP)
      call tensor_buffer(INT1,root=infpar%master,comm=comm)
      call tensor_alloc_mem(intarr1,INT1)
      call tensor_buffer(intarr1,INT1)
      call tensor_buffer(INT2)
      call tensor_alloc_mem(intarr2,INT2)
      call tensor_buffer(intarr2,INT2,finalize=.true.)
      
      call tensor_init(AUX,intarr2,INT2)

      call tensor_free_mem(intarr2)

      call tensor_zero(AUX)

      call lspdm_extract_decnp_indices_occ(AUX,A,INT1,intarr1)

      call tensor_free_mem(intarr1)
      call tensor_free(AUX)

   CASE(JOB_GET_COMBINEDT1T2_1)

      call lspdm_get_combined_SingleDouble_amplitudes(A,B,C)

   CASE(JOB_GET_COMBINEDT1T2_2)

      call tensor_alloc_mem(intarr1,2)

      call tensor_buffer(intarr1,2,root=infpar%master,comm=comm)
      call tensor_init(AUX,intarr1,2)
      call tensor_buffer(AUX%elm1,AUX%nelms,finalize=.true.)

      call tensor_free_mem(intarr1)

      call lspdm_get_combined_SingleDouble_amplitudes(AUX,A,B)

      call tensor_free(AUX)
   CASE(JOB_GET_MP2_ST_GUESS)
      call tensor_mpi_bcast(INT1,infpar%master,infpar%lg_comm)
      call tensor_mpi_bcast(LOG1,infpar%master,infpar%lg_comm)
      call lspdm_get_starting_guess(A,B,C,D,INT1,LOG1)
   CASE(JOB_tensor_rand)
      call tensor_rand_tiled_dist(A)

   case(JOB_LSPDM_INIT_GLOBAL_BUFFER)
      call lspdm_init_global_buffer(comm,.false.)
   case(JOB_LSPDM_FREE_GLOBAL_BUFFER)
      call lspdm_free_global_buffer(comm,.false.)
   case(JOB_SET_TENSOR_SEG_LENGTH)
      call tensor_set_global_segment_length(comm,LIN1)
   case(JOB_SET_TENSOR_DEBUG_TRUE)
      call tensor_set_debug_mode_true(comm,.false.)
   case(JOB_SET_TENSOR_ALWAYS_SYNC_TRUE)
      call tensor_set_always_sync_true(comm,.false.)
   case(JOB_SET_TENSOR_BACKEND_TRUE)
      call tensor_set_dil_backend_true(comm,.false.)

   CASE DEFAULT
        call lsquit("ERROR(pdm_tensor_slave): Unknown job",-1)
   END SELECT
#endif
end subroutine pdm_tensor_slave
