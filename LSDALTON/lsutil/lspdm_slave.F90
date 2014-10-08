
subroutine pdm_array_slave(comm)
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

  use tensor_interface_module

   IMPLICIT NONE
   INTEGER(kind=ls_mpik),intent(in) :: comm
   TYPE(array)  :: A, B, C, D, AUX
   CHARACTER    :: T(2)
   INTEGER      :: JOB, i, j
   real(REALK)  :: AF, AB(2), REAL1, REAL2
   real(REALK),pointer :: dummy(:)
   logical :: logi
   integer, pointer  :: intarr1(:),   intarr2(:),   intarr3(:)
   integer           :: intarr1_size, intarr2_size, intarr3_size
   logical :: loc
   character (4) :: at 
   integer       :: it
#ifdef VAR_MPI
   loc = (infpar%parent_comm /= MPI_COMM_NULL)

   call time_start_phase(PHASE_COMM)
   CALL PDM_ARRAY_SYNC(comm,JOB,A,B,C,D,loc_addr=loc) !Job is output
   call time_start_phase(PHASE_WORK)

   SELECT CASE(JOB)
   CASE(JOB_PC_DEALLOC_DENSE)
      call memory_deallocate_array_dense_pc(A)
   CASE(JOB_PC_ALLOC_DENSE)
      call memory_allocate_array_dense_pc(A)
   CASE(JOB_INIT_ARR_PC)
      call mem_alloc(intarr1,A%mode)
      intarr1 =A%dims
      call arr_free_aux(A)
      call array_init_standard(A,intarr1,A%mode,MASTER_ACCESS) 
      call mem_dealloc(intarr1)
   CASE(JOB_INIT_ARR_TILED)
      call mem_alloc(intarr2,A%mode)
      call mem_alloc(intarr1,A%mode)
      intarr2=A%tdim
      intarr1 =A%dims
      call arr_free_aux(A)
      call array_init_tiled(A,intarr1,A%mode,at,it,MASTER_ACCESS,intarr2,A%zeros) 
      call mem_dealloc(intarr2)
      call mem_dealloc(intarr1)
   CASE(JOB_FREE_ARR_STD)
      call array_free_pdm(A) 
   CASE(JOB_FREE_ARR_PDM)
      call array_free_pdm(A) 
   CASE(JOB_INIT_ARR_REPLICATED)
      call mem_alloc(intarr1,A%mode)
      intarr1 =A%dims
      call arr_free_aux(A)
      call array_init_replicated(A,intarr1,A%mode,MASTER_ACCESS) 
      call mem_dealloc(intarr1)
   CASE(JOB_PRINT_MEM_INFO1)
      call print_mem_per_node(DECinfo%output,.false.)
   CASE(JOB_PRINT_MEM_INFO2)
      call mem_alloc(dummy,1)
      call print_mem_per_node(DECinfo%output,.false.,dummy)
      call mem_dealloc(dummy)
   CASE(JOB_GET_NRM2_TILED)
      AF=array_tiled_pdm_get_nrm2(A)
   CASE(JOB_DATA2TILED_DIST)
      !dummy has to be allocated, otherwise seg faults might occur with 
      !some compilers
      call mem_alloc(dummy,1)
      !call cp_data2tiled(A,dummy,A%dims,A%mode,.true.)
      print *,"not necessary"
      call mem_dealloc(dummy)
   CASE(JOB_GET_TILE_SEND)
      !i,dummy and j are just dummy arguments
      call array_get_tile(A,i,dummy,j)
   CASE(JOB_PRINT_TI_NRM)
      call array_tiled_pdm_print_ti_nrm(A,0)
   CASE(JOB_SYNC_REPLICATED)
      call array_sync_replicated(A)
   CASE(JOB_GET_NORM_REPLICATED)
      AF=array_print_norm_repl(A)
   CASE(JOB_PREC_DOUBLES_PAR)
      call precondition_doubles_parallel(A,B,C,D)
   CASE(JOB_DDOT_PAR)
      AF=array_ddot_par(A,B,0)
   CASE(JOB_ADD_PAR)
      intarr1_size = A%mode
      call mem_alloc(intarr1,intarr1_size)
      call time_start_phase(PHASE_COMM)
      call ls_mpiinitbuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call ls_mpi_buffer(intarr1,intarr1_size,infpar%master)
      call ls_mpi_buffer(REAL1,infpar%master)
      call ls_mpi_buffer(REAL2,infpar%master)
      call ls_mpifinalizebuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call time_start_phase(PHASE_WORK)

      call array_add_par(REAL1,A,REAL2,B,intarr1)

      call mem_dealloc(intarr1)
   CASE(JOB_CP_ARR)
      call array_cp_tiled(A,B)
   CASE(JOB_ARRAY_ZERO)
      call array_zero_tiled_dist(A)
   CASE(JOB_GET_CC_ENERGY)
      AF = get_cc_energy_parallel(A,B,C)
   CASE(JOB_GET_MP2_ENERGY)
      AF = get_mp2_energy_parallel(A,B)
   CASE(JOB_GET_FRAG_CC_ENERGY)
      !the counterpart to this buffer is in get_fragment_cc_energy
      call time_start_phase(PHASE_COMM)
      call ls_mpiinitbuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call ls_mpi_buffer(intarr1_size,infpar%master)
      call mem_alloc(intarr1,intarr1_size)
      call ls_mpi_buffer(intarr1,intarr1_size,infpar%master)
      call ls_mpi_buffer(intarr2_size,infpar%master)
      call mem_alloc(intarr2,intarr2_size)
      call ls_mpi_buffer(intarr2,intarr2_size,infpar%master)
      call ls_mpifinalizebuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call time_start_phase(PHASE_WORK)

      AF = get_fragment_cc_energy_parallel(A,B,C,i,j,intarr1,intarr2)

      call mem_dealloc(intarr1)
      call mem_dealloc(intarr2)
   CASE(JOB_CHANGE_ACCESS_TYPE)
      call change_access_type_td(A,i)
   CASE(JOB_ARRAY_SCALE)
      call ls_mpibcast(AF,infpar%master,infpar%lg_comm)
      call array_scale_td(A,AF)
   CASE(JOB_GET_RPA_ENERGY)
      write(*,*) 'Johannes RPA ene'
      AF = get_rpa_energy_parallel(A,B)
   CASE(JOB_GET_SOS_ENERGY)
      write(*,*) 'Johannes SOS ene'
      AF = get_sosex_cont_parallel(A,B)
   CASE(JOB_ARR_CONTRACT_SIMPLE)
      call time_start_phase(PHASE_COMM)
      call ls_mpiinitbuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call ls_mpi_buffer(intarr1_size,infpar%master)
      intarr2_size = intarr1_size
      intarr3_size = C%mode
      call mem_alloc(intarr1,intarr1_size)
      call mem_alloc(intarr2,intarr2_size)
      call mem_alloc(intarr3,intarr3_size)
      call ls_mpi_buffer(intarr1,intarr1_size,infpar%master)
      call ls_mpi_buffer(intarr2,intarr2_size,infpar%master)
      call ls_mpi_buffer(intarr3,intarr3_size,infpar%master)
      call ls_mpi_buffer(REAL1,infpar%master)
      call ls_mpi_buffer(REAL2,infpar%master)
      call ls_mpifinalizebuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call time_start_phase(PHASE_WORK)

      call lspdm_array_contract_simple(REAL1,A,B,intarr1,intarr2,intarr1_size,REAL2,C,intarr3)

      call mem_dealloc(intarr1)
      call mem_dealloc(intarr2)

   CASE DEFAULT
        call lsquit("ERROR(pdm_array_slave): Unknown job",-1)
   END SELECT
#endif
end subroutine pdm_array_slave
