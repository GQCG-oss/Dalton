
subroutine pdm_tensor_slave(comm)
  use precision
  use lstiming
  !use matrix_operations_scalapack, only: BLOCK_SIZE, SLGrid, DLEN_
  use memory_handling, only: mem_alloc,mem_dealloc
  use matrix_operations, only: mtype_scalapack, matrix_type
  use dec_typedef_module
#ifdef VAR_MPI
  use infpar_module
  use lsmpi_type
#endif

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
   integer, pointer    :: intarr1(:), intarr2(:), intarr3(:), intarr4(:)
   integer             :: INT1,       INT2,       INT3,       INT4
   real(tensor_dp)   :: REAL1,      REAL2
   logical             :: LOG1
   logical :: loc
   character (4) :: at 
#ifdef VAR_MPI
   loc = (infpar%parent_comm /= MPI_COMM_NULL)


   call time_start_phase(PHASE_COMM)
   CALL PDM_tensor_SYNC(comm,JOB,A,B,C,D,loc_addr=loc) !Job is output
   call time_start_phase(PHASE_WORK)

   SELECT CASE(JOB)
   CASE(JOB_PC_DEALLOC_DENSE)
      call memory_deallocate_tensor_dense_pc(A)
   CASE(JOB_PC_ALLOC_DENSE)
      call memory_allocate_tensor_dense_pc(A)
   CASE(JOB_INIT_tensor_PC)
      call mem_alloc(intarr1,A%mode)
      intarr1 =A%dims
      call tensor_free_aux(A)
      print *,"this is deprecated"
      stop 0
      !call tensor_init_standard(A,intarr1,A%mode,AT_MASTER_ACCESS) 
      call mem_dealloc(intarr1)

   CASE(JOB_INIT_TENSOR_TILED)

      call mem_alloc(intarr2,A%mode)
      call mem_alloc(intarr1,A%mode)
      intarr2=A%tdim
      intarr1 =A%dims
      call tensor_free_aux(A)

      call ls_mpibcast(INT1,infpar%master,infpar%lg_comm)

      if(INT1==-1)then
         call tensor_init_tiled(A,intarr1,int(A%mode),at,int(INT1,kind=tensor_standard_int),&
            &AT_MASTER_ACCESS,.false.,tdims=intarr2) 
      else
         call tensor_init_tiled(A,intarr1,int(A%mode),at,int(INT1,kind=tensor_standard_int),&
            &AT_MASTER_ACCESS,.false.,tdims=intarr2,force_offset=INT1) 
      endif

      call mem_dealloc(intarr2)
      call mem_dealloc(intarr1)

   CASE(JOB_FREE_TENSOR_STD)
      call tensor_free_pdm(A) 
   CASE(JOB_FREE_tensor_PDM)
      call tensor_free_pdm(A) 
   CASE(JOB_INIT_tensor_REPLICATED)
      call mem_alloc(intarr1,A%mode)
      intarr1 =A%dims
      call tensor_free_aux(A)
      call tensor_init_replicated(A,intarr1,int(A%mode),AT_MASTER_ACCESS,.false.) 
      call mem_dealloc(intarr1)
   CASE(JOB_PRINT_MEM_INFO1)
      call print_mem_per_node(DECinfo%output,.false.)
   CASE(JOB_PRINT_MEM_INFO2)
      call mem_alloc(realar1,1)
      call print_mem_per_node(DECinfo%output,.false.,realar1)
      call mem_dealloc(realar1)
   CASE(JOB_GET_NRM2_TILED)
      REAL1 = tensor_tiled_pdm_get_nrm2(A)
   CASE(JOB_DATA2TILED_DIST)
      !realar1 has to be allocated, otherwise seg faults might occur with 
      !some compilers
      call mem_alloc(realar1,1)
      !call cp_data2tiled(A,realar1,A%dims,A%mode,.true.)
      print *,"not necessary"
      call mem_dealloc(realar1)
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
      call mem_alloc(intarr1,INT1)
      call time_start_phase(PHASE_COMM)
      call ls_mpiinitbuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call ls_mpi_buffer(intarr1,INT1,infpar%master)
      call ls_mpi_buffer(REAL1,infpar%master)
      call ls_mpi_buffer(REAL2,infpar%master)
      call ls_mpifinalizebuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call time_start_phase(PHASE_WORK)

      call tensor_add_par(REAL1,A,REAL2,B,intarr1)

      call mem_dealloc(intarr1)

   CASE(JOB_DMUL_PAR)
      INT1 = A%mode
      call mem_alloc(intarr1,INT1)
      call time_start_phase(PHASE_COMM)
      call ls_mpiinitbuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call ls_mpi_buffer(intarr1,INT1,infpar%master)
      call ls_mpi_buffer(REAL1,infpar%master)
      call ls_mpi_buffer(REAL2,infpar%master)
      call ls_mpi_buffer(INT2,infpar%master)
      call mem_alloc(realar1,INT2)
      call ls_mpi_buffer(realar1,INT2,infpar%master)
      call ls_mpifinalizebuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call time_start_phase(PHASE_WORK)

      call tensor_dmul_par(REAL1,A,REAL2,realar1,B,intarr1)

      call mem_dealloc(intarr1)
      call mem_dealloc(realar1)

   CASE(JOB_HMUL_PAR)
      call time_start_phase(PHASE_COMM)
      call ls_mpiinitbuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call ls_mpi_buffer(REAL1,infpar%master)
      call ls_mpi_buffer(REAL2,infpar%master)
      call ls_mpifinalizebuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call time_start_phase(PHASE_WORK)
      call tensor_hmul_par(REAL1,A,B,REAL2,C)
   CASE(JOB_CP_ARR)
      INT1 = A%mode
      call mem_alloc(intarr1,INT1)
      call ls_mpibcast(intarr1,INT1,infpar%master,infpar%lg_comm)
      call tensor_cp_tiled(A,B,intarr1)
      call mem_dealloc(intarr1)
   CASE(JOB_tensor_ZERO)
      call tensor_zero_tiled_dist(A)
   CASE(JOB_GET_CC_ENERGY)
      REAL1 = get_cc_energy_parallel(A,B,C)
   CASE(JOB_GET_MP2_ENERGY)
      REAL1 = get_cc_energy_parallel(A,B)
   CASE(JOB_GET_FRAG_CC_ENERGY)
      !the counterpart to this buffer is in get_fragment_cc_energy
      call time_start_phase(PHASE_COMM)
      call ls_mpiinitbuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call ls_mpi_buffer(INT1,infpar%master)
      call mem_alloc(intarr1,INT1)
      call ls_mpi_buffer(intarr1,INT1,infpar%master)
      call ls_mpi_buffer(INT2,infpar%master)
      call mem_alloc(intarr2,INT2)
      call ls_mpi_buffer(intarr2,INT2,infpar%master)
      call ls_mpifinalizebuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call time_start_phase(PHASE_WORK)

      REAL1 = get_fragment_cc_energy_parallel(A,B,C,INT1,INT2,intarr1,intarr2)

      call mem_dealloc(intarr1)
      call mem_dealloc(intarr2)
   CASE(JOB_CHANGE_ACCESS_TYPE)
      call change_access_type_td(A,INT1)
   CASE(JOB_tensor_SCALE)
      call ls_mpibcast(REAL1,infpar%master,infpar%lg_comm)
      call tensor_scale_td(A,REAL1)
   CASE(JOB_GET_RPA_ENERGY)
      write(*,*) 'Johannes RPA ene'
      REAL1 = get_rpa_energy_parallel(A,B)
   CASE(JOB_GET_SOS_ENERGY)
      write(*,*) 'Johannes SOS ene'
      REAL1 = get_sosex_cont_parallel(A,B)
   CASE(JOB_tensor_CONTRACT_SIMPLE)
      call time_start_phase(PHASE_COMM)
      call ls_mpiinitbuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call ls_mpi_buffer(INT1,infpar%master)
      INT2 = INT1
      INT3 = C%mode
      call mem_alloc(intarr1,INT1)
      call mem_alloc(intarr2,INT2)
      call mem_alloc(intarr3,INT3)
      call ls_mpi_buffer(intarr1,INT1,infpar%master)
      call ls_mpi_buffer(intarr2,INT2,infpar%master)
      call ls_mpi_buffer(intarr3,INT3,infpar%master)
      call ls_mpi_buffer(REAL1,infpar%master)
      call ls_mpi_buffer(REAL2,infpar%master)
      call ls_mpi_buffer(LOG1, infpar%master)
      call ls_mpifinalizebuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call time_start_phase(PHASE_WORK)

      if(LOG1)then
         call lspdm_tensor_contract_simple(REAL1,A,B,intarr1,intarr2,INT1,REAL2,C,intarr3,force_sync=LOG1)
      else
         call lspdm_tensor_contract_simple(REAL1,A,B,intarr1,intarr2,INT1,REAL2,C,intarr3)
      endif

      call mem_dealloc(intarr1)
      call mem_dealloc(intarr2)
      call mem_dealloc(intarr3)
   CASE(JOB_tensor_CONTRACT_BDENSE)
      call time_start_phase(PHASE_COMM)
      call ls_mpiinitbuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call ls_mpi_buffer(INT1,infpar%master)
      INT2 = INT1
      INT3 = B%mode
      call mem_alloc(intarr1,INT1)
      call mem_alloc(intarr2,INT2)
      call mem_alloc(intarr3,INT3)
      call ls_mpi_buffer(intarr1,INT1,infpar%master)
      call ls_mpi_buffer(intarr2,INT2,infpar%master)
      call ls_mpi_buffer(intarr3,INT3,infpar%master)
      call ls_mpi_buffer(REAL1,infpar%master)
      call ls_mpi_buffer(REAL2,infpar%master)
      call ls_mpi_buffer(LOG1, infpar%master)
      call ls_mpi_buffer(INT4, infpar%master)
      call mem_alloc(intarr4,INT4)
      call ls_mpi_buffer(intarr4,INT4, infpar%master)

      call tensor_init(AUX,intarr4,int4)

      call ls_mpi_buffer(AUX%elm1,AUX%nelms,infpar%master)

      call ls_mpifinalizebuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call time_start_phase(PHASE_WORK)

      if(LOG1)then
         call lspdm_tensor_contract_simple(REAL1,A,AUX,intarr1,intarr2,INT1,REAL2,B,intarr3,force_sync=LOG1)
      else
         call lspdm_tensor_contract_simple(REAL1,A,AUX,intarr1,intarr2,INT1,REAL2,B,intarr3)
      endif

      call tensor_free(AUX)

      call mem_dealloc(intarr1)
      call mem_dealloc(intarr2)
      call mem_dealloc(intarr3)
      call mem_dealloc(intarr4)
   CASE(JOB_TENSOR_EXTRACT_VEOS)
      call ls_mpiinitbuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call ls_mpi_buffer(INT1,infpar%master)
      call mem_alloc(intarr1,INT1)
      call ls_mpi_buffer(intarr1,INT1,infpar%master)
      call ls_mpi_buffer(INT2,infpar%master)
      call mem_alloc(intarr2,INT2)
      call ls_mpi_buffer(intarr2,INT2,infpar%master)
      call ls_mpifinalizebuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)

      call tensor_init(AUX,intarr2,INT2)

      call mem_dealloc(intarr2)

      call tensor_zero(AUX)

      call lspdm_extract_eos_indices_virt(AUX,A,INT1,intarr1)

      call mem_dealloc(intarr1)
      call tensor_free(AUX)

   CASE(JOB_TENSOR_EXTRACT_OEOS)
      call ls_mpiinitbuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call ls_mpi_buffer(INT1,infpar%master)
      call mem_alloc(intarr1,INT1)
      call ls_mpi_buffer(intarr1,INT1,infpar%master)
      call ls_mpi_buffer(INT2,infpar%master)
      call mem_alloc(intarr2,INT2)
      call ls_mpi_buffer(intarr2,INT2,infpar%master)
      call ls_mpifinalizebuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      
      call tensor_init(AUX,intarr2,INT2)

      call mem_dealloc(intarr2)

      call tensor_zero(AUX)

      call lspdm_extract_eos_indices_occ(AUX,A,INT1,intarr1)

      call mem_dealloc(intarr1)
      call tensor_free(AUX)

   CASE(JOB_TENSOR_EXTRACT_VDECNP)
      call ls_mpiinitbuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call ls_mpi_buffer(INT1,infpar%master)
      call mem_alloc(intarr1,INT1)
      call ls_mpi_buffer(intarr1,INT1,infpar%master)
      call ls_mpi_buffer(INT2,infpar%master)
      call mem_alloc(intarr2,INT2)
      call ls_mpi_buffer(intarr2,INT2,infpar%master)
      call ls_mpifinalizebuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      
      call tensor_init(AUX,intarr2,INT2)

      call mem_dealloc(intarr2)

      call tensor_zero(AUX)

      call lspdm_extract_decnp_indices_virt(AUX,A,INT1,intarr1)

      call mem_dealloc(intarr1)
      call tensor_free(AUX)

   CASE(JOB_TENSOR_EXTRACT_ODECNP)
      call ls_mpiinitbuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call ls_mpi_buffer(INT1,infpar%master)
      call mem_alloc(intarr1,INT1)
      call ls_mpi_buffer(intarr1,INT1,infpar%master)
      call ls_mpi_buffer(INT2,infpar%master)
      call mem_alloc(intarr2,INT2)
      call ls_mpi_buffer(intarr2,INT2,infpar%master)
      call ls_mpifinalizebuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      
      call tensor_init(AUX,intarr2,INT2)

      call mem_dealloc(intarr2)

      call tensor_zero(AUX)

      call lspdm_extract_decnp_indices_occ(AUX,A,INT1,intarr1)

      call mem_dealloc(intarr1)
      call tensor_free(AUX)

   CASE(JOB_GET_COMBINEDT1T2_1)

      call lspdm_get_combined_SingleDouble_amplitudes(A,B,C)

   CASE(JOB_GET_COMBINEDT1T2_2)

      call mem_alloc(intarr1,2)

      call ls_mpiinitbuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call ls_mpi_buffer(intarr1,2,infpar%master)
      call tensor_init(AUX,intarr1,2)
      call ls_mpi_buffer(AUX%elm1,AUX%nelms,infpar%master)
      call ls_mpifinalizebuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)

      call mem_dealloc(intarr1)

      call lspdm_get_combined_SingleDouble_amplitudes(AUX,A,B)

      call tensor_free(AUX)
   CASE(JOB_GET_MP2_ST_GUESS)
      call ls_mpibcast(INT1,infpar%master,infpar%lg_comm)
      call ls_mpibcast(LOG1,infpar%master,infpar%lg_comm)
      call lspdm_get_starting_guess(A,B,C,D,INT1,LOG1)
   CASE(JOB_tensor_rand)
      call tensor_rand_tiled_dist(A)

   CASE DEFAULT
        call lsquit("ERROR(pdm_tensor_slave): Unknown job",-1)
   END SELECT
#endif
end subroutine pdm_tensor_slave
